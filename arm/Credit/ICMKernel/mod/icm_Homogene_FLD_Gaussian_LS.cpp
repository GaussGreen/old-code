

#include "ARMKernel\glob\firsttoinc.h" 
#include "ICMKernel\mod\icm_Homogene_FastLossDistrib_Gaussian.h"
#include "ICMKernel\glob\icm_maths.h"


//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

//----------------------------------------------------------------------------------------------------------------------
// Long Short
//----------------------------------------------------------------------------------------------------------------------
ICM_Gauss1FLossDistrib_LJ::ICM_Gauss1FLossDistrib_LJ(const int& nbnames,
							const ARM_Vector& pdefault,
							const ARM_Vector& beta,
							// const std::vector<double>& LossRates,
							const ARM_Vector& LossRates,
							const std::vector<int>& SortedIndice,
							const int& discretizationstep,
							const int& CopulaType,
							const qIntegratorChoice&	IntegrationMethod,
							const int& IntegrationStep)
{
	Init();
	Set(nbnames,pdefault,beta,LossRates, SortedIndice, discretizationstep,CopulaType,IntegrationMethod,IntegrationStep);
}


// -----------------------------------------------------------------------
//LongShort
// -----------------------------------------------------------------------
void ICM_Gauss1FLossDistrib_LJ::Set(const int& nbnames,
										   const ARM_Vector& pdef,
										   const ARM_Vector&  beta,
										   const ARM_Vector & LossRates,
										   const std::vector<int>& SortedIndice,
										   const int& discretizationstep ,
										   const int& CopulaType ,
										   const qIntegratorChoice&	IntegrationMethod ,
										   const int& IntegrationStep)
{
	Set(nbnames,pdef,beta,LossRates,discretizationstep ,CopulaType ,IntegrationMethod ,IntegrationStep);
	ComputeNbNameShort(LossRates);
	SetSortedIndices(SortedIndice);
}


// -------------------------------------------------------------------
//	Long Short CDO
// -------------------------------------------------------------------

double
ICM_Gauss1FLossDistrib_LJ::compute_expectedlosstranche_LongShort(const double& tranche_up, 
												  const double& tranche_down, 
												  const double& lossunit
												  // const double& minloss
//												  ICM_QMatrix<double>* ShiftMatrix,
												  // const int& Tenor,
												  // const int& Issuer
												  )
{
	double exp_loss_tranche_down = 0.;
	double exp_loss_tranche_up = 0.;
	 
	int lup=floor(tranche_up/lossunit);
	int ldown=floor(tranche_down/lossunit);
	int lmin= floor(fabs(its_LSminloss));

	int i=0;	
	
	its_ProbCond->ResizeWithInitialize(lmin+lup+1,0.);
	its_ProbCond_Down->ResizeWithInitialize(lmin+ldown+1,0.);

	its_lossdistrib.resize(lmin+lup+1); 
	its_lossdistrib_Down.resize(lmin+ldown+1);
	
	// Compute the Distribution (by default, ldown = 0)
	if (tranche_down)
		compute_distrib_LongShort(lup, ldown, lmin);
	else
		compute_distrib_LongShort(lup, lmin);
	
	//View the loss ditribution
	/* #ifdef _DEBUG
	FILE *stream = fopen("c:\\temp\\LossDistribLongShort.txt", "w+");
	this->View("",stream);
	fclose(stream);
	#endif */
 

	int l	=	0;
	double	loss_level;

	// first loop, deal with beta_down and up
	for (l=1+lmin; l<=ldown+lmin; l++) 
	{
		loss_level	=	(l-lmin)*lossunit;
		exp_loss_tranche_down	+=	its_lossdistrib_Down[l] * loss_level;		
		exp_loss_tranche_up		+=	its_lossdistrib[l] * loss_level;
	}
		
	// second loop, deal only with beat_up
	for (l=lmin+ldown+1;l<=lup+lmin;l++)
	{
		loss_level	=	(l-lmin)*lossunit;
		exp_loss_tranche_up		+=	its_lossdistrib[l] * loss_level;
	}

	// add the tail
	if (tranche_down)
		exp_loss_tranche_down	+=	tranche_down * its_taildistrib_Down;
	exp_loss_tranche_up		+=	tranche_up * its_taildistrib;

	return (exp_loss_tranche_up - exp_loss_tranche_down);
}

//Distribution

double ICM_Gauss1FLossDistrib_LJ::Compute_pk(const int& ind_x, const int& ind_sortedname)
{
	double tmp_barrier=0.;
	double pk=0.;
	double x=TheIntegrator.GetAbscissa(ind_x)*SQRT2;
	
	//Barrier
	tmp_barrier = its_coeff_a[ind_sortedname] - its_coeff_b[ind_sortedname] * x;
	
	// conditional default probability
	pk=NAG_cumul_normal(tmp_barrier);

	return pk;

}

double ICM_Gauss1FLossDistrib_LJ::Compute_pk_down(const int& ind_x, const int& ind_sortedname)
{
	double tmp_barrier=0.;
	double pk=0.;
	double x=TheIntegrator.GetAbscissa(ind_x)*SQRT2;
	
	//Barrier
	tmp_barrier = its_coeff_a_down[ind_sortedname] - its_coeff_b_down[ind_sortedname] * x;
	
	// conditional default probability
	pk=NAG_cumul_normal(tmp_barrier);

	return pk;

}

void
ICM_Gauss1FLossDistrib_LJ::compute_distrib_LongShort(const int& lup, const int& lmin)
{
	double tmpPrbLoss=0.,cumul_distrib=0.;
	int min_lup =0,max_lup =0;
	int l = 0;

	int ind_x = 0;
	int ind_name = 0;
	int ind_loss = 0;
	int ind_sortedname = 0;
	max_lup		=	lup+lmin;
	its_ind_Loss = 0;

	switch (itsCopulaType)
	{
		case qNO_COPULA:	
		case qGAUSSIAN:

			//Prob Cond computation

			//First Loop : underlying factor (integration)
			for (ind_x=0; ind_x<TheIntegrator.GetIntegrationStep(); ind_x++)
			{
				//Second Loop : portfolio size
				InitProbCondPortSize0(lmin, ind_x);
				double pk=0.;

				//Short Issuer Part
				for (ind_name=1; ind_name<=this->its_nbnames_short; ind_name++)
				{
					//Sorted Index
					ind_sortedname = this->its_sorted_indices[ind_name-1];
					pk=0.;
					pk = Compute_pk(ind_x,ind_sortedname);

					//Conditionnal loss distribution on short part
					for (ind_loss=0; ind_loss<=lmin; ind_loss++)
					{
						ComputeProbCondShort(ind_x, ind_loss, ind_name, lmin, ind_sortedname, pk);
					}
				}
			
				//Long Issuer Part
				for (ind_name=its_nbnames_short+1; ind_name<=this->its_nbnames; ind_name++)
				{
					//Sorted Index
					ind_sortedname = this->its_sorted_indices[ind_name-1];
					pk=0.;
					pk = Compute_pk(ind_x,ind_sortedname);

					for (ind_loss=0; ind_loss<=max_lup; ind_loss++)
					{
						ComputeProbCond(ind_x, ind_loss, ind_name, ind_sortedname, pk);
					}
				}					
			}

			//Loos distribution computation
			for (ind_loss=0; ind_loss<=max_lup; ind_loss++)
			{
				its_lossdistrib[ind_loss] = 0.;
				for (ind_x=0; ind_x<TheIntegrator.GetIntegrationStep(); ind_x++)
				{
					its_lossdistrib[ind_loss] += its_ProbCond->Elt(its_nbnames,ind_x,ind_loss)*TheIntegrator.GetWeight(ind_x);
				}	
				its_lossdistrib[ind_loss]*=ONEOVERSQRTPI;
			}
	}

	for (l=0;l<=max_lup;l++)
	{
		cumul_distrib += its_lossdistrib[l];
	}

	its_taildistrib			=	1.0 - cumul_distrib;
}

//Util Long Short
void ICM_Gauss1FLossDistrib_LJ::InitProbCondPortSize0(const int& lmin, const int& ind_x)
{
	its_ProbCond->Elt(0,ind_x,lmin) = 1;
	its_ProbCond_Down->Elt(0,ind_x,lmin) = 1;
}

void ICM_Gauss1FLossDistrib_LJ::ComputeProbCond(const int& ind_x, const int& ind_loss, const int& ind_name, const int& ind_sortedname, const double& pk)
	{
	double tmp1=0.;
	double tmp2=0.;
	//double tmp_barrier=0.;
	//double pk=0.;
	
	//double x=TheIntegrator.GetAbscissa(ind_x)*SQRT2;
	// int ind_lastloss = ind_loss - GetIndLossRates()[ind_sortedname];
	int ind_lastloss = ind_loss - GetIntLossRates()[ind_sortedname];

	//Barrier
	//tmp_barrier = its_coeff_a[ind_sortedname] - its_coeff_b[ind_sortedname] * x;
	
	// conditional default probability
	tmp1=its_ProbCond->Elt(ind_name-1,ind_x,ind_loss);
	//pk=NAG_cumul_normal(tmp_barrier);
	tmp2		=	tmp1 * (1.-pk);

	if (ind_lastloss >= 0)
		tmp2	+=	pk * its_ProbCond->Elt(ind_name-1,ind_x,ind_lastloss);

	its_ProbCond->Elt(ind_name,ind_x,ind_loss)=tmp2;
}

void ICM_Gauss1FLossDistrib_LJ::ComputeProbCondDown(const int& ind_x, const int& ind_loss, const int& ind_name, const int& ind_sortedname, const double& pk_down)
{
	double tmp1_down=0.;
	double tmp2_down=0;
	//double tmp_barrier_down=0.;
	//double pk_down=0.;
	//double x=TheIntegrator.GetAbscissa(ind_x)*SQRT2;

	int ind_lastloss = ind_loss - GetIntLossRates()[ind_sortedname];
	//Barrier
	//tmp_barrier_down = its_coeff_a_down[ind_sortedname] - its_coeff_b_down[ind_sortedname] * x;		
	
	// conditional default probability
	tmp1_down=its_ProbCond_Down->Elt(ind_name-1,ind_x,ind_loss);
	//pk_down=NAG_cumul_normal(tmp_barrier_down);
	tmp2_down	=	tmp1_down * (1.-pk_down);

	if (ind_lastloss >= 0)
		tmp2_down	+=	pk_down * its_ProbCond_Down->Elt(ind_name-1,ind_x,ind_lastloss);

	its_ProbCond_Down->Elt(ind_name,ind_x,ind_loss)=tmp2_down;
}

void ICM_Gauss1FLossDistrib_LJ::ComputeProbCondShort(const int& ind_x, const int& ind_loss, const int& ind_name, const int& lmin, const int& ind_sortedname, const double& pk)
{
	double tmp1=0.;
	double tmp2=0.;
//	double tmp_barrier=0.;
//	double pk=0.;
	
//	double x=TheIntegrator.GetAbscissa(ind_x)*SQRT2;
	int ind_lastloss = ind_loss - GetIntLossRates()[ind_sortedname];

	//Barrier
//	tmp_barrier = its_coeff_a[ind_sortedname] - its_coeff_b[ind_sortedname] * x;
	
	// conditional default probability
	tmp1=its_ProbCond->Elt(ind_name-1,ind_x,ind_loss);
//	pk=NAG_cumul_normal(tmp_barrier);
	tmp2		=	tmp1 * (1.-pk);

	if (ind_lastloss <= lmin)
		tmp2	+=	pk * its_ProbCond->Elt(ind_name-1,ind_x,ind_lastloss);

	its_ProbCond->Elt(ind_name,ind_x,ind_loss)=tmp2;
}
	
void ICM_Gauss1FLossDistrib_LJ::ComputeProbCondDownShort(const int& ind_x, const int& ind_loss, const int& ind_name, const int& lmin, const int& ind_sortedname,const double& pk_down)
{
	double tmp1_down=0.;
	double tmp2_down=0;
//	double tmp_barrier_down=0.;
//	double pk_down=0.;
//	double x=TheIntegrator.GetAbscissa(ind_x)*SQRT2;

	int ind_lastloss = ind_loss - GetIntLossRates()[ind_sortedname];
	//Barrier
//	tmp_barrier_down = its_coeff_a_down[ind_sortedname] - its_coeff_b_down[ind_sortedname] * x;		
	
	// conditional default probability
	tmp1_down=its_ProbCond_Down->Elt(ind_name-1,ind_x,ind_loss);
//	pk_down=NAG_cumul_normal(tmp_barrier_down);
	tmp2_down	=	tmp1_down * (1.-pk_down);

	if (ind_lastloss <= lmin)
		tmp2_down	+=	pk_down * its_ProbCond_Down->Elt(ind_name-1,ind_x,ind_lastloss);

	its_ProbCond_Down->Elt(ind_name,ind_x,ind_loss)=tmp2_down;
}

void
ICM_Gauss1FLossDistrib_LJ::compute_distrib_LongShort(const int& lup,const int& ldown, const int& lmin)
{
	double tmpPrbLoss=0.,cumul_distrib=0.;
	double tmpPrbLoss_down=0.,cumul_distrib_down=0.;
	int min_lup =0,max_lup =0, max_ldown=0;
	int l = 0;

	int ind_x = 0;
	int ind_name = 0;
	int ind_loss = 0;
	int ind_sortedname = 0;
	max_ldown	=	ldown+lmin;
	max_lup		=	lup+lmin;
	its_ind_Loss = 0;

	switch (itsCopulaType)
	{
		case qNO_COPULA:	
		case qGAUSSIAN:

			//Prob Cond computation

			//First Loop : underlying factor (integration)
			for (ind_x=0; ind_x<TheIntegrator.GetIntegrationStep(); ind_x++)
			{
				//Second Loop : portfolio size
				InitProbCondPortSize0(lmin, ind_x);
				double pk=0.;
				double pk_down=0.;

				//Short Issuer Part
				for (ind_name=1; ind_name<=this->its_nbnames_short; ind_name++)
				{
					//Sorted Index
					ind_sortedname = this->its_sorted_indices[ind_name-1];
					pk=0.;
					pk_down=0.;
					pk = Compute_pk(ind_x,ind_sortedname);
					pk_down = Compute_pk_down(ind_x,ind_sortedname);

					//Conditionnal loss distribution on short part
					for (ind_loss=0; ind_loss<=lmin; ind_loss++)
					{
						ComputeProbCondShort(ind_x, ind_loss, ind_name, lmin, ind_sortedname, pk);
						ComputeProbCondDownShort(ind_x, ind_loss, ind_name, lmin, ind_sortedname, pk_down);
					}
				}
			
				//Long Issuer Part
				for (ind_name=its_nbnames_short+1; ind_name<=this->its_nbnames; ind_name++)
				{
					//Sorted Index
					ind_sortedname = this->its_sorted_indices[ind_name-1];
					pk=0.;
					pk_down=0.;
					pk = Compute_pk(ind_x,ind_sortedname);
					pk_down = Compute_pk_down(ind_x,ind_sortedname);

					for (ind_loss=0; ind_loss<=max_ldown; ind_loss++)
					{
						ComputeProbCond(ind_x, ind_loss, ind_name, ind_sortedname, pk);
						ComputeProbCondDown(ind_x, ind_loss, ind_name, ind_sortedname, pk_down);
					}
					for (ind_loss=max_ldown+1; ind_loss<=max_lup; ind_loss++)
					{
						ComputeProbCond(ind_x, ind_loss, ind_name, ind_sortedname, pk);
					}
				}					
			}

			//Loos distribution computation
			for (ind_loss=0; ind_loss<=max_ldown; ind_loss++)
			{
				its_lossdistrib[ind_loss] = 0.;
				its_lossdistrib_Down[ind_loss] = 0.;
				for (ind_x=0; ind_x<TheIntegrator.GetIntegrationStep(); ind_x++)
				{
					its_lossdistrib[ind_loss] += its_ProbCond->Elt(its_nbnames,ind_x,ind_loss)*TheIntegrator.GetWeight(ind_x);
					its_lossdistrib_Down[ind_loss] += its_ProbCond_Down->Elt(its_nbnames,ind_x,ind_loss)*TheIntegrator.GetWeight(ind_x);
				}	
				its_lossdistrib[ind_loss]*=ONEOVERSQRTPI;
				its_lossdistrib_Down[ind_loss]*=ONEOVERSQRTPI;
			}
			for (ind_loss=max_ldown+1; ind_loss<=max_lup; ind_loss++)
			{
				its_lossdistrib[ind_loss] = 0.;
				for (ind_x=0; ind_x<TheIntegrator.GetIntegrationStep(); ind_x++)
					its_lossdistrib[ind_loss] += its_ProbCond->Elt(its_nbnames,ind_x,ind_loss)*TheIntegrator.GetWeight(ind_x);
				its_lossdistrib[ind_loss]*=ONEOVERSQRTPI;
			}
	}

	for (l=0;l<=max_ldown;l++)
	{
		cumul_distrib += its_lossdistrib[l];
		cumul_distrib_down += its_lossdistrib_Down[l];
	}
	for (l=max_ldown+1;l<=max_lup;l++)
	{
		cumul_distrib += its_lossdistrib[l];
	}

	its_taildistrib			=	1.0 - cumul_distrib;
	its_taildistrib_Down	=	1.0 - cumul_distrib_down;
	
}