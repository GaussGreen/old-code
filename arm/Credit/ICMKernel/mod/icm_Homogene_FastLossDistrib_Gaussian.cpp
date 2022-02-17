#include "ARMKernel\glob\firsttoinc.h" 
#include "ICMKernel\mod\icm_Homogene_FastLossDistrib_Gaussian.h"
#include "ICMKernel\glob\icm_maths.h"
#include "ARMKernel\util\interpol.h"
#include "ICMKernel\util\icm_HermiteIntegration.h"

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////


void
ICM_Gauss1FLossDistrib_LJ::Init()
{
	itsIntegrationMethod	=	qGAUSS_HERMITE;		// Gauss-Legendre, Hermite, etc.
	itsIntegrationStep		=	20;

	itsNbStrikes		=	0;

	itsCopulaType		=	qNO_COPULA;
	itsFreedomDegree	=	4;

	itsStrikes.clear();
	itsBetasMatrix		=	NULL;

	its_ProbCond_Down = NULL;
	its_ProbCond_Perturb_Down = NULL;
	its_lossdistrib_perturb_Down = NULL;
	
	its_taildistrib_Down	=	0.0;

	its_beta_Down.clear();		 
	its_lossdistrib_Down.clear(); 
	its_taildistrib_perturb_Down.clear();
	its_coeff_a.clear();
	its_coeff_b.clear();
	its_coeff_a_down.clear();
	its_coeff_b_down.clear();
	its_coeff_a_perturb.clear();
	its_coeff_a_down_perturb.clear();
	itsStrikes.clear();

	//Forward Collateral --------------------------------------------------------------
	its_collat_fwd_taildistrib_Down=0.;

	its_collat_fwd_ProbCond_Down=NULL;
	its_collat_fwd_ProbCond_Perturb_Down=NULL;
	its_collat_fwd_lossdistrib_perturb_Down=NULL;

	its_collat_fwd_beta_Down.clear();	
	its_collat_fwd_lossdistrib_Down.clear();
	its_collat_fwd_taildistrib_perturb_Down.clear();

	its_collat_fwd_coeff_a.clear();
	its_collat_fwd_coeff_b.clear();
	its_collat_fwd_coeff_a_down.clear();
	its_collat_fwd_coeff_b_down.clear();
	its_collat_fwd_coeff_a_perturb.clear();
	its_collat_fwd_coeff_a_down_perturb.clear();
	its_collat_fwdStrikes.clear();

	its_beta_Down_t2.clear();		 
	its_collat_fwd_beta_Down_t2.clear();
	its_coeff_b_perturb.clear();
	its_coeff_b_down_perturb.clear();
	its_coeff_a_t2.clear();
	its_coeff_b_t2.clear();
	its_coeff_a_down_t2.clear();
	its_coeff_b_down_t2.clear();
	its_coeff_a_perturb_t2.clear();
	its_coeff_a_down_perturb_t2.clear();
	its_coeff_b_perturb_t2.clear();
	its_coeff_b_down_perturb_t2.clear();	

	its_ProbCond_stepup = NULL;
}

//----------------------------------------------------------------------------------------------------------------------
// Constructor
//----------------------------------------------------------------------------------------------------------------------
ICM_Gauss1FLossDistrib_LJ::ICM_Gauss1FLossDistrib_LJ(const int& nbnames,
										   const ARM_Vector& pdef,
										   const ARM_Vector& beta,
										   const ARM_Vector& LossRates,
										   const int& discretizationstep ,
										   const int& CopulaType ,
										   const qIntegratorChoice&	IntegrationMethod ,
										   const int& IntegrationStep)
{
	Init();

	Set(nbnames,pdef,beta,LossRates,discretizationstep,CopulaType,IntegrationMethod,IntegrationStep);

}


void ICM_Gauss1FLossDistrib_LJ::Set(const int& nbnames,

										   const ARM_Vector&  pdef,
										   const ARM_Vector&  beta,
										   const ARM_Vector & LossRates,
										   const int& discretizationstep ,
										   const int& CopulaType ,
										   const qIntegratorChoice&	IntegrationMethod ,
										   const int& IntegrationStep)
{
	int size;
	if (discretizationstep	==	0) size	=	20;
	else size	=	discretizationstep;

	if (nbnames != GetNbNames())
	{
		ICM_Distribution::Set(nbnames);
	
		if (its_lossdistrib_perturb_Down) delete its_lossdistrib_perturb_Down;
		its_lossdistrib_perturb_Down= new ICM_QMatrix<double>(1,nbnames);

		its_taildistrib_perturb_Down.resize(nbnames);
		its_beta_Down.resize(nbnames);
		its_coeff_a.resize(nbnames);
		its_coeff_b.resize(nbnames);
		its_coeff_a_down.resize(nbnames);
		its_coeff_b_down.resize(nbnames);
		its_coeff_a_perturb.resize(nbnames);
		its_coeff_a_down_perturb.resize(nbnames);
	}

	if ((nbnames != its_nbnames) || (itsIntegrationStep != size) || 
		// 17783 (!its_ProbCond)||(!its_ProbCond_Perturb)||(!its_ProbCond_Perturb_Down)||(!its_ProbCond_Down))
		(!its_ProbCond)||(!its_ProbCond_Perturb_Down)||(!its_ProbCond_Down))
	{
		if (its_ProbCond) delete its_ProbCond;
		its_ProbCond= new ICM_QCubix<double>(nbnames+1,size,1,0.);

		// 17783 if (its_ProbCond_Perturb) delete its_ProbCond_Perturb;
		// 17783 its_ProbCond_Perturb= new ICM_QCubix<double>(nbnames+1,size,1,0.);

		if (its_ProbCond_Perturb_Down) delete its_ProbCond_Perturb_Down;
		its_ProbCond_Perturb_Down= new ICM_QCubix<double>(nbnames+1,size,1,0.);

		if (its_ProbCond_Down) delete its_ProbCond_Down;
		its_ProbCond_Down= new ICM_QCubix<double>(nbnames+1,size,1,0.);
	}

	SetNbNames(nbnames);

	itsCopulaType = CopulaType;
	
	itsIntegrationStep = IntegrationStep;
	itsIntegrationMethod = IntegrationMethod;
	// Si le nombre de pas et la méthode utilisés ne sont pas compatibles, la parité du nombre de pas est prioritaire
	if ( (itsIntegrationStep % 2 == 0) && (itsIntegrationMethod == qGAUSS_LEGENDRE))
		itsIntegrationMethod = qGAUSS_HERMITE;
	else if ( (itsIntegrationStep % 2 == 1) && (itsIntegrationMethod == qGAUSS_HERMITE))
		itsIntegrationMethod = qGAUSS_LEGENDRE;
	

	SetUniqueBeta(beta); 
	SetPdefAtMaturity(pdef); 
	compute_min_pdef(GetPdefAtMaturity());

	SetIntLossRates(LossRates);
	check_homogeneous(LossRates);
	compute_barrier();		// Copula dependant

	its_lossdistrib_Down.resize(1);
}


void ICM_Gauss1FLossDistrib_LJ::BitwiseCopy(const ARM_Object* src)
{

    ICM_Gauss1FLossDistrib_LJ* srcdistrib = (ICM_Gauss1FLossDistrib_LJ *) src;

	itsCopulaType = srcdistrib->itsCopulaType;
	itsIntegrationStep = srcdistrib->itsIntegrationStep;
	itsFreedomDegree = srcdistrib->itsFreedomDegree;
	itsNbStrikes = srcdistrib->itsNbStrikes;
	its_taildistrib_Down = srcdistrib->its_taildistrib_Down;
	itsIntegrationMethod = srcdistrib->itsIntegrationMethod;
	
	if (srcdistrib->itsBetasMatrix)
		itsBetasMatrix = (ICM_QMatrix<double>*) srcdistrib->itsBetasMatrix->Clone();

	if (srcdistrib->its_ProbCond_Down)
		its_ProbCond_Down = (ICM_QCubix<double>*) srcdistrib->its_ProbCond_Down->Clone();

	if (srcdistrib->its_ProbCond_Perturb_Down)
		its_ProbCond_Perturb_Down = (ICM_QCubix<double>*) srcdistrib->its_ProbCond_Perturb_Down->Clone();

	if (srcdistrib->its_lossdistrib_perturb_Down)
		its_lossdistrib_perturb_Down = (ICM_QMatrix<double>*) srcdistrib->its_lossdistrib_perturb_Down->Clone();

	its_beta_Down=srcdistrib->its_beta_Down;
	its_lossdistrib_Down=srcdistrib->its_lossdistrib_Down;
	its_taildistrib_perturb_Down=srcdistrib->its_taildistrib_perturb_Down;
	
	its_coeff_a=srcdistrib->its_coeff_a;
	its_coeff_b=srcdistrib->its_coeff_b;
	its_coeff_a_down=srcdistrib->its_coeff_a_down;
	its_coeff_b_down=srcdistrib->its_coeff_b_down;
	its_coeff_a_perturb=srcdistrib->its_coeff_a_perturb;
	its_coeff_a_down_perturb=srcdistrib->its_coeff_a_down_perturb;
	itsStrikes=srcdistrib->itsStrikes;

	//fwd collateral
	its_collat_fwd_taildistrib_Down = srcdistrib->its_collat_fwd_taildistrib_Down;
	
	if (srcdistrib->its_collat_fwd_ProbCond_Down)
		its_collat_fwd_ProbCond_Down = (ICM_QCubix<double>*) srcdistrib->its_collat_fwd_ProbCond_Down->Clone();

	if (srcdistrib->its_collat_fwd_ProbCond_Perturb_Down)
		its_collat_fwd_ProbCond_Perturb_Down = (ICM_QCubix<double>*) srcdistrib->its_collat_fwd_ProbCond_Perturb_Down->Clone();

	if (srcdistrib->its_lossdistrib_perturb_Down)
		its_collat_fwd_lossdistrib_perturb_Down = (ICM_QMatrix<double>*) srcdistrib->its_collat_fwd_lossdistrib_perturb_Down->Clone();

	its_collat_fwd_beta_Down=srcdistrib->its_collat_fwd_beta_Down;
	its_collat_fwd_lossdistrib_Down=srcdistrib->its_collat_fwd_lossdistrib_Down;
	its_collat_fwd_taildistrib_perturb_Down=srcdistrib->its_collat_fwd_taildistrib_perturb_Down;
	
	its_collat_fwd_coeff_a=srcdistrib->its_collat_fwd_coeff_a;
	its_collat_fwd_coeff_b=srcdistrib->its_collat_fwd_coeff_b;
	its_collat_fwd_coeff_a_down=srcdistrib->its_collat_fwd_coeff_a_down;
	its_collat_fwd_coeff_b_down=srcdistrib->its_collat_fwd_coeff_b_down;
	its_collat_fwd_coeff_a_perturb=srcdistrib->its_collat_fwd_coeff_a_perturb;
	its_collat_fwd_coeff_a_down_perturb=srcdistrib->its_collat_fwd_coeff_a_down_perturb;
	its_collat_fwdStrikes=srcdistrib->its_collat_fwdStrikes;

}

// -------------
//	Copy Method 
// -------------
void ICM_Gauss1FLossDistrib_LJ::Copy(const ARM_Object* src)
{
	ICM_Distribution::Copy(src);
    BitwiseCopy(src);
}


ARM_Object* ICM_Gauss1FLossDistrib_LJ::Clone(void)
{
     ICM_Gauss1FLossDistrib_LJ* theClone = new ICM_Gauss1FLossDistrib_LJ();

     theClone->Copy(this);
 
     return(theClone);
}


//----------------------------------------------------------------------------------------------------------------------

void ICM_Gauss1FLossDistrib_LJ::SetSmileParameters(const int& NbStrikes,
												const vector<double>&	Strikes,
												ICM_QMatrix<double>* BetaMatrix)
{
	// Nb Strikes
	itsNbStrikes = NbStrikes;

	// Strikes
	itsStrikes = Strikes;
	
	// Smile Betas Matrix	
	if (itsBetasMatrix)
		delete itsBetasMatrix;
	itsBetasMatrix = (ICM_QMatrix<double>*) BetaMatrix->Clone();
}


void
ICM_Gauss1FLossDistrib_LJ::compute_distrib(const int& lup)
{
	switch (itsIntegrationMethod)
	{
	case qGAUSS_LEGENDRE:
	case qGAUSS_HERMITE:
	case qTRAPEZE:
	default:
		compute_distrib_Gauss_Legendre(lup);
	break;
	//compute_distrib_Hermite(lup);
	}
}

// -------------------------------------------------------------------
// GAUSSIAN COPULA
// -------------------------------------------------------------------

//	Loss Level is 0.0

extern "C" void 
compute_cond_distribZeroLoss_Gaussian_Integrator(void* Param, double x, double& res)
{
	// get parameters from stack
	ICM_Gauss1FLossDistrib_LJ*	TheModel;
	TheModel	=	(ICM_Gauss1FLossDistrib_LJ*)((*(AddressVector*)Param)[0]);

	res	=	TheModel->compute_cond_distribZeroLoss_Gaussian(x,Param);
}


double
ICM_Gauss1FLossDistrib_LJ::compute_cond_distribZeroLoss_Gaussian(double x, void* params)
{
	double	pk=0.;
	double	tmp_barrier=0.;

	its_ind_x+=1;

	unsigned long sizemax=-1;
	(*(AddressVector*)params).GetCount(sizemax);
	ICM_QMatrix<double>* matrix = (ICM_QMatrix<double>*)((*(AddressVector*)params)[sizemax-2]);
	int* nointegrationcoef = (int*)((*(AddressVector*)params)[sizemax-1]);

	double CondProbaZeroLoss=1.;

	its_ProbCond->SetElt(0,its_ind_x,0,CondProbaZeroLoss);

	if (itsIntegrationMethod ==	qGAUSS_HERMITE)
		x	*=	SQRT2;

	for (int k=1;k<=its_nbnames;k++) 
	{

		if (*nointegrationcoef>=0)
		{
			if ((*matrix)(k-1,*nointegrationcoef) !=-999.)
			{   pk = (*matrix)(k-1,*nointegrationcoef);}
			else
			{
				tmp_barrier=(its_barrier[k-1]-its_unique_beta[k-1]*x)/sqrt(1.-its_unique_beta[k-1]*its_unique_beta[k-1]);
				(*matrix)(k-1,*nointegrationcoef) = pk=NAG_cumul_normal(tmp_barrier);

				//add ------------------
				if (IsTSR() && GetT2())
				{
					double tmp_barrier_t2=0.;
					double tmp_barrier_end=0.;

					tmp_barrier_t2	= (its_barrier_t2[k-1]-its_unique_beta_t2[k-1]*x)/sqrt(1.-its_unique_beta_t2[k-1]*its_unique_beta_t2[k-1]);
					pk += NAG_cumul_normal(tmp_barrier_t2);
	
					pk -=	NAG_bivariate_normal_dist(tmp_barrier,tmp_barrier_t2,RHO_DIST(GetT1(),GetT2()));
					(*matrix)(k-1,*nointegrationcoef) = pk;
				}
				//add -------------------

			} 
		}
		else
		{
		tmp_barrier=(its_barrier[k-1]-its_unique_beta[k-1]*x)/sqrt(1.-its_unique_beta[k-1]*its_unique_beta[k-1]);
		pk=NAG_cumul_normal(tmp_barrier);

		//add ------------------
		if (IsTSR() && GetT2())
		{
			double tmp_barrier_t2=0.;
			double tmp_barrier_end=0.;

			tmp_barrier_t2	= (its_barrier_t2[k-1]-its_unique_beta_t2[k-1]*x)/sqrt(1.-its_unique_beta_t2[k-1]*its_unique_beta_t2[k-1]);
			pk += NAG_cumul_normal(tmp_barrier_t2);
	
			pk -=	NAG_bivariate_normal_dist(tmp_barrier,tmp_barrier_t2,RHO_DIST(GetT1(),GetT2()));
		}
		//add -------------------

		}

		CondProbaZeroLoss *= (1.-pk);
		its_ProbCond->SetElt(k,its_ind_x,0,CondProbaZeroLoss);
	}

	return its_ProbCond->Elt(its_nbnames,its_ind_x,0);	
}


//	Loss Level is 0.0

extern "C" void 
compute_cond_distribZeroLoss_Gaussian_Smile(void* Param, double x, double& ProbCond_Up, double& ProbCond_Down)
{
	int k;
	double	pk=0.;
	double	tmp_barrier=0.;

	double CondProbaZeroLoss	=	1.;
	double CondProbaZeroLoss_Down=	1.;

	// get parameters from stack
	ICM_Gauss1FLossDistrib_LJ*	TheModel;
	TheModel	=	(ICM_Gauss1FLossDistrib_LJ*)((*(AddressVector*)Param)[0]);

	int*	its_ind_x;
	its_ind_x	=	(int*)((*(AddressVector*)Param)[1]);
	(*its_ind_x) += 1;

	int*	its_ind_Loss;
	its_ind_Loss	=	(int*)((*(AddressVector*)Param)[2]);

	int itsIntegrationMethod = TheModel->GetIntegrationMethod();

	if (itsIntegrationMethod	==	qGAUSS_HERMITE)
		x	*=	SQRT2;

	TheModel->SetProbCond_Elt(0, *its_ind_x,0,CondProbaZeroLoss);
	TheModel->SetProbCond_Down_Elt(0, *its_ind_x,0,CondProbaZeroLoss_Down);

	for (k=0;k<TheModel->GetNbNames();k++) 
	{
		tmp_barrier	=	TheModel->GetCoeff_a(k) - TheModel->GetCoeff_b(k) * x;
		pk=NAG_cumul_normal(tmp_barrier);

		if (TheModel->IsTSR() && TheModel->GetT2())
		{
			double tmp_barrier_t2=0.;
			double tmp_barrier_end=0.;

			tmp_barrier_t2	=	TheModel->GetCoeff_a_t2(k) - TheModel->GetCoeff_b_t2(k) * x;
			pk += NAG_cumul_normal(tmp_barrier_t2);

			pk -=	NAG_bivariate_normal_dist(tmp_barrier,tmp_barrier_t2,RHO_DIST(TheModel->GetT1(),TheModel->GetT2()));
		}

		CondProbaZeroLoss *= (1.-pk);
		TheModel->SetProbCond_Elt(k+1, *its_ind_x,0,CondProbaZeroLoss);

		tmp_barrier	=	TheModel->GetCoeff_a_down(k) - TheModel->GetCoeff_b_down(k) * x;
		pk=NAG_cumul_normal(tmp_barrier);

		if (TheModel->IsTSR() && TheModel->GetT2())
		{
			double tmp_barrier_t2=0.;
			double tmp_barrier_end=0.;

			tmp_barrier_t2	=	TheModel->GetCoeff_a_down_t2(k) - TheModel->GetCoeff_b_down_t2(k) * x;
			pk += NAG_cumul_normal(tmp_barrier_t2);

			pk -=	NAG_bivariate_normal_dist(tmp_barrier,tmp_barrier_t2,RHO_DIST(TheModel->GetT1(),TheModel->GetT2()));
		}

		CondProbaZeroLoss_Down *= (1.-pk);
		TheModel->SetProbCond_Down_Elt(k+1, *its_ind_x,0,CondProbaZeroLoss_Down);
	}

	ProbCond_Up		=	CondProbaZeroLoss;	
	ProbCond_Down	=	CondProbaZeroLoss_Down;	
}


// -------------------------------------------------------------------
//	Strictly Positive Loss Level
// -------------------------------------------------------------------

extern "C" void 
compute_cond_distrib_Gaussian_Integrator(void* Param, double x, double& res)
{
	// get parameters from stack
	ICM_Gauss1FLossDistrib_LJ*	TheModel;
	TheModel	=	(ICM_Gauss1FLossDistrib_LJ*)((*(AddressVector*)Param)[0]);

	res	=	TheModel->compute_cond_distrib_Gaussian(x,Param);
}


double
ICM_Gauss1FLossDistrib_LJ::compute_cond_distrib_Gaussian(double x,void* params)
{

	its_ind_x+=1;

	int	k;
	double tmp1=0.; // proba de loss=L pour un panier de taille k
	double tmp2=0.,tmp3=0.; // proba de loss=L-1 pour un panier de taille k
	double pk=0.;     // proba cond de défaut du nom k
	int wk=0,ind_lastloss=0;

	double tmp_barrier=0.,tmp4=0.;

	unsigned long sizemax=-1;
	(*(AddressVector*)params).GetCount(sizemax);
	ICM_QMatrix<double>* matrix = (ICM_QMatrix<double>*)((*(AddressVector*)params)[sizemax-2]);
	int* nointegrationcoef = (int*)((*(AddressVector*)params)[sizemax-1]);

	its_ProbCond->SetElt(0,its_ind_x,its_ind_Loss,0.);

	// its_ind_Loss: Loss Distribution for a basket of size (0) knowing factor value its_ind_x
	tmp1	=	its_ProbCond->Elt(0, its_ind_x, its_ind_Loss);


	// HERMITE
	if (itsIntegrationMethod	==	qGAUSS_HERMITE)
		x	*=	its_sqrt2;
	{
		for (k=1;k<=its_nbnames;k++) 
		{
			// current level of loss - loss rate of names k
			double LossRates_inf = its_int_lossrates[k-1];
			double LossRates_delta = its_dbl_lossrates[k-1] - LossRates_inf;

			ind_lastloss = its_ind_Loss - its_int_lossrates[k-1];

			if (*nointegrationcoef>=0)
			{
				// 17783 
				pk = (*matrix)(k-1,*nointegrationcoef);
				if (pk==-999.)
				{
					tmp_barrier=(its_barrier[k-1]-its_unique_beta[k-1]*x)/sqrt(1.-its_unique_beta[k-1]*its_unique_beta[k-1]);
					pk = NAG_cumul_normal(tmp_barrier);

					if (IsTSR() && GetT2())
					{
						double tmp_barrier_t2=0.;
						double tmp_barrier_end=0.;
						tmp_barrier_t2	=(its_barrier_t2[k-1]-its_unique_beta_t2[k-1]*x)/sqrt(1.-its_unique_beta_t2[k-1]*its_unique_beta_t2[k-1]);
						pk += NAG_cumul_normal(tmp_barrier_t2);
						pk -=	NAG_bivariate_normal_dist(tmp_barrier,tmp_barrier_t2,RHO_DIST(GetT1(),GetT2()));
						
					}
					(*matrix)(k-1,*nointegrationcoef) = pk;
				}
				/** 
				if ((*matrix)(k-1,*nointegrationcoef) !=-999.)
				{   pk = (*matrix)(k-1,*nointegrationcoef);}
				else
				{
					tmp_barrier=(its_barrier[k-1]-its_unique_beta[k-1]*x)/sqrt(1.-its_unique_beta[k-1]*its_unique_beta[k-1]);
					(*matrix)(k-1,*nointegrationcoef) = pk = NAG_cumul_normal(tmp_barrier);

					if (IsTSR() && GetT2())
					{
						double tmp_barrier_t2=0.;
						double tmp_barrier_end=0.;

						tmp_barrier_t2	=(its_barrier_t2[k-1]-its_unique_beta_t2[k-1]*x)/sqrt(1.-its_unique_beta_t2[k-1]*its_unique_beta_t2[k-1]);
						pk += NAG_cumul_normal(tmp_barrier_t2);

						pk -=	NAG_bivariate_normal_dist(tmp_barrier,tmp_barrier_t2,RHO_DIST(GetT1(),GetT2()));
						(*matrix)(k-1,*nointegrationcoef) = pk;
					}
				} **/ 
			}
			else
			{
			// barrier
			tmp_barrier=(its_barrier[k-1]-its_unique_beta[k-1]*x)/sqrt(1.-its_unique_beta[k-1]*its_unique_beta[k-1]);
			// conditional default probability
			pk=NAG_cumul_normal(tmp_barrier);

			if (IsTSR() && GetT2())
			{
				double tmp_barrier_t2=0.;
				double tmp_barrier_end=0.;

				tmp_barrier_t2	=(its_barrier_t2[k-1]-its_unique_beta_t2[k-1]*x)/sqrt(1.-its_unique_beta_t2[k-1]*its_unique_beta_t2[k-1]);
				pk += NAG_cumul_normal(tmp_barrier_t2);

				pk -=	NAG_bivariate_normal_dist(tmp_barrier,tmp_barrier_t2,RHO_DIST(GetT1(),GetT2()));
			}
			}

			// its_ind_Loss: Loss Distribution for a basket of size (k-1) knowing factor value its_ind_x
			// tmp1

			// Loss Distribution for a basket of size k knowing factor value its_ind_x
			tmp2	=	tmp1 * (1-pk);
			
			if (ind_lastloss >= 0)
			{
				tmp2	+=	pk * its_ProbCond->Elt(k-1,its_ind_x,ind_lastloss)* (1.-LossRates_delta);

				if (ind_lastloss > 0)
				{
					tmp2	+=	pk * its_ProbCond->Elt(k-1,its_ind_x,ind_lastloss-1)* LossRates_delta;
				}
			}

			its_ProbCond->SetElt(k,its_ind_x,its_ind_Loss, tmp2);

			tmp1 = tmp2;
		}
	}

	return tmp1;
}


extern "C" void 
compute_cond_distrib_Gaussian_Smile(void* Param, double x, double& ProbCond_Up, double& ProbCond_Down)
{
	// get parameters from stack
	ICM_Gauss1FLossDistrib_LJ*	TheModel;
	TheModel	=	(ICM_Gauss1FLossDistrib_LJ*)((*(AddressVector*)Param)[0]);

	int*	its_ind_x;
	its_ind_x	=	(int*)((*(AddressVector*)Param)[1]);
	(*its_ind_x) += 1;

	int*	its_ind_Loss;
	its_ind_Loss	=	(int*)((*(AddressVector*)Param)[2]);

	ICM_QMatrix<double>* matrix_up = (ICM_QMatrix<double>*)((*(AddressVector*)Param)[3]);
	ICM_QMatrix<double>& matrix_down = *(ICM_QMatrix<double>*)((*(AddressVector*)Param)[4]);
	int* nointegrationcoef = (int*)((*(AddressVector*)Param)[5]);

	int	k;
	double tmp1, tmp1_down; // proba de loss=L pour un panier de taille k
	double tmp2, tmp2_down; // proba de loss=L-1 pour un panier de taille k
	double pk, pk_down;     // proba cond de défaut du nom k
	int ind_lastloss;

	double tmp_barrier, tmp_barrier_down;

	int itsIntegrationMethod = TheModel->GetIntegrationMethod();
	if (itsIntegrationMethod	==	qGAUSS_HERMITE)
		x	*=	SQRT2;

	TheModel->SetProbCond_Elt(0,*its_ind_x,*its_ind_Loss,0.);
	TheModel->SetProbCond_Down_Elt(0,*its_ind_x,*its_ind_Loss,0.);

	// its_ind_Loss: Loss Distribution for a basket of size (0) knowing factor value its_ind_x
	tmp1		=	TheModel->GetProbCond_Elt(0, *its_ind_x, *its_ind_Loss);
	tmp1_down	=	TheModel->GetProbCond_Down_Elt(0, *its_ind_x, *its_ind_Loss);

	{
		for (k=0;k<TheModel->GetNbNames();k++) 
		{
			// current level of loss - loss rate of names k
			double LossRates_inf = TheModel->GetIntLossRates()[k];
			double LossRates_delta = TheModel->GetDblLossRates()[k] - LossRates_inf;

			ind_lastloss = *its_ind_Loss - TheModel->GetIntLossRates()[k];

			if (*nointegrationcoef>=0) 
			{
				// 17783: same with less matrix access. 
				pk =(*matrix_up)(k,*nointegrationcoef) ; 
				if (pk==-999.)
				{
					tmp_barrier = TheModel->GetCoeff_a(k) - TheModel->GetCoeff_b(k) * x;
					pk = NAG_cumul_normal(tmp_barrier);

					if (TheModel->IsTSR() && TheModel->GetT2())
					{
						double tmp_barrier_t2=0.;
						double tmp_barrier_end=0.;

						tmp_barrier_t2	=	TheModel->GetCoeff_a_t2(k) - TheModel->GetCoeff_b_t2(k) * x;
						pk += NAG_cumul_normal(tmp_barrier_t2);
						pk -=	NAG_bivariate_normal_dist(tmp_barrier,tmp_barrier_t2,RHO_DIST(TheModel->GetT1(),TheModel->GetT2()));

					}
					(*matrix_up)(k,*nointegrationcoef) = pk;
				}
				/** 
				if ((*matrix_up)(k,*nointegrationcoef) !=-999.)
				{  pk = (*matrix_up)(k,*nointegrationcoef);}
				else
				{
					tmp_barrier = TheModel->GetCoeff_a(k) - TheModel->GetCoeff_b(k) * x;
					(*matrix_up)(k,*nointegrationcoef) = pk = NAG_cumul_normal(tmp_barrier);

					if (TheModel->IsTSR() && TheModel->GetT2())
					{
						double tmp_barrier_t2=0.;
						double tmp_barrier_end=0.;

						tmp_barrier_t2	=	TheModel->GetCoeff_a_t2(k) - TheModel->GetCoeff_b_t2(k) * x;
						pk += NAG_cumul_normal(tmp_barrier_t2);

						pk -=	NAG_bivariate_normal_dist(tmp_barrier,tmp_barrier_t2,RHO_DIST(TheModel->GetT1(),TheModel->GetT2()));

						(*matrix_up)(k,*nointegrationcoef) = pk;
					}
				}
				**/ 
				//	17783
				pk_down = matrix_down(k,*nointegrationcoef); 
				if (pk_down ==-999.) 
				{
					tmp_barrier_down = TheModel->GetCoeff_a_down(k) - TheModel->GetCoeff_b_down(k) * x;
					pk_down = NAG_cumul_normal(tmp_barrier_down);

					if (TheModel->IsTSR() && TheModel->GetT2())
					{		
						double tmp_barrier_t2=0.;
						double tmp_barrier_end=0.;

						tmp_barrier_t2	=	TheModel->GetCoeff_a_down_t2(k) - TheModel->GetCoeff_b_down_t2(k) * x;
						pk_down += NAG_cumul_normal(tmp_barrier_t2);

						pk_down -=	NAG_bivariate_normal_dist(tmp_barrier_down,tmp_barrier_t2,RHO_DIST(TheModel->GetT1(),TheModel->GetT2()));
					}
					matrix_down(k,*nointegrationcoef) = pk_down;
				}

				/** 
				if ((*matrix_down)(k,*nointegrationcoef) !=-999.)
				{ pk_down = (*matrix_down)(k,*nointegrationcoef);}
				else
				{
					tmp_barrier_down = TheModel->GetCoeff_a_down(k) - TheModel->GetCoeff_b_down(k) * x;
					(*matrix_down)(k,*nointegrationcoef) = pk_down = NAG_cumul_normal(tmp_barrier_down);

					if (TheModel->IsTSR() && TheModel->GetT2())
					{		
						double tmp_barrier_t2=0.;
						double tmp_barrier_end=0.;

						tmp_barrier_t2	=	TheModel->GetCoeff_a_down_t2(k) - TheModel->GetCoeff_b_down_t2(k) * x;
						pk_down += NAG_cumul_normal(tmp_barrier_t2);

						pk_down -=	NAG_bivariate_normal_dist(tmp_barrier_down,tmp_barrier_t2,RHO_DIST(TheModel->GetT1(),TheModel->GetT2()));
						(*matrix_down)(k,*nointegrationcoef) = pk_down;
					}
				}**/ 

			}
			else
			{
				// barrier
				tmp_barrier = TheModel->GetCoeff_a(k) - TheModel->GetCoeff_b(k) * x;
				tmp_barrier_down = TheModel->GetCoeff_a_down(k) - TheModel->GetCoeff_b_down(k) * x;

				// conditional default probability
				pk=NAG_cumul_normal(tmp_barrier);
				pk_down=NAG_cumul_normal(tmp_barrier_down);

				if (TheModel->IsTSR() && TheModel->GetT2())
				{
					double tmp_barrier_t2=0.;
					double tmp_barrier_end=0.;

					tmp_barrier_t2	=	TheModel->GetCoeff_a_t2(k) - TheModel->GetCoeff_b_t2(k) * x;
					pk += NAG_cumul_normal(tmp_barrier_t2);

					pk -=	NAG_bivariate_normal_dist(tmp_barrier,tmp_barrier_t2,RHO_DIST(TheModel->GetT1(),TheModel->GetT2()));
				}

				if (TheModel->IsTSR() && TheModel->GetT2())
				{
					double tmp_barrier_t2=0.;
					double tmp_barrier_end=0.;

					tmp_barrier_t2	=	TheModel->GetCoeff_a_down_t2(k) - TheModel->GetCoeff_b_down_t2(k) * x;
					pk_down += NAG_cumul_normal(tmp_barrier_t2);

					pk_down -=	NAG_bivariate_normal_dist(tmp_barrier_down,tmp_barrier_t2,RHO_DIST(TheModel->GetT1(),TheModel->GetT2()));
				}
			}

			// its_ind_Loss: Loss Distribution for a basket of size (k-1) knowing factor value its_ind_x
			// tmp1

			// Loss Distribution for a basket of size k knowing factor value its_ind_x
			tmp2		=	tmp1 * (1.-pk);
			tmp2_down	=	tmp1_down * (1.-pk_down);
			
			if (ind_lastloss >= 0)
			{
				tmp2	+=	pk * TheModel->GetProbCond_Elt(k, *its_ind_x, ind_lastloss) * (1.-LossRates_delta);
				tmp2_down	+=	pk_down * TheModel->GetProbCond_Down_Elt(k, *its_ind_x, ind_lastloss) * (1.-LossRates_delta);

				if (ind_lastloss >0)
				{
					tmp2	+=	pk * TheModel->GetProbCond_Elt(k, *its_ind_x, ind_lastloss-1) * (LossRates_delta);
					tmp2_down	+=	pk_down * TheModel->GetProbCond_Down_Elt(k, *its_ind_x, ind_lastloss-1) * (LossRates_delta);
				}

			}
			// else: the loss of name k is too big

			TheModel->SetProbCond_Elt(k+1, *its_ind_x, *its_ind_Loss, tmp2);
			TheModel->SetProbCond_Down_Elt(k+1, *its_ind_x, *its_ind_Loss, tmp2_down);

			tmp1 = tmp2;
			tmp1_down = tmp2_down;
		}
	}

	ProbCond_Up		=	 tmp1;
	ProbCond_Down	=	 tmp1_down;
}

// -------------------------------------------------------------------
// STUDENT COPULA
// -------------------------------------------------------------------

//	Loss Level is 0.0

double
ICM_Gauss1FLossDistrib_LJ::compute_cond_distribZeroLoss_Student(const double& x)
{
	double	pk=0.;
	double	tmp_barrier=0.;

	its_ind_x+=1;

	double CondProbaZeroLoss=1.;

	its_ProbCond->SetElt(0,its_ind_x,0,CondProbaZeroLoss);

	for (int k=1;k<=its_nbnames;k++) 
	{
		tmp_barrier=(its_barrier[k-1]-its_unique_beta[k-1]*x)/sqrt(1.-its_unique_beta[k-1]*its_unique_beta[k-1]);

		// Cumulative Student distribution: nag_t_prob
		pk = NAG_prob_students_t( tmp_barrier, itsFreedomDegree );

		CondProbaZeroLoss *= (1.-pk);
		its_ProbCond->SetElt(k,its_ind_x,0,CondProbaZeroLoss);
	}

	return its_ProbCond->Elt(its_nbnames,its_ind_x,0);
}


// -------------------------------------------------------------------
//	Strictly Positive Loss Level
// -------------------------------------------------------------------

double
ICM_Gauss1FLossDistrib_LJ::compute_cond_distrib_Student(const double& x)
{

	its_ind_x+=1;

	int	k;
	double tmp1=0.; // proba de loss=L pour un panier de taille k
	double tmp2=0.,tmp3=0.; // proba de loss=L-1 pour un panier de taille k
	double pk=0.;     // proba cond de défaut du nom k
	int wk=0,ind_lastloss=0;

	double tmp_barrier=0.,tmp4=0.;

	its_ProbCond->SetElt(0,its_ind_x,its_ind_Loss,0.);

	// its_ind_Loss: Loss Distribution for a basket of size (0) knowing factor value its_ind_x
	tmp1	=	its_ProbCond->Elt(0, its_ind_x, its_ind_Loss);

/*	if (its_ishomogeneous)
	{
		for (k=1;k<=its_nbnames;k++) 
		{
			// barrier
			tmp_barrier=(its_barrier[k-1]-its_unique_beta[k-1]*x)/sqrt(1.-its_unique_beta[k-1]*its_unique_beta[k-1]);

			// conditional default probability
			// Cumulative Student distribution
			pk = NAG_prob_students_t( tmp_barrier, itsFreedomDegree);

			// its_ind_Loss: Loss Distribution for a basket of size (k-1) knowing factor value its_ind_x
			// tmp1

			// Loss Distribution for a basket of size k knowing factor value its_ind_x
			tmp2	=	tmp1 * (1.-pk) + pk * its_ProbCond->Elt(k-1,its_ind_x,its_ind_Loss-1);

			its_ProbCond->SetElt(k, its_ind_x, its_ind_Loss, tmp2);

			tmp1 = tmp2;
		}
	}
	else
*/	{
		for (k=1;k<=its_nbnames;k++) 
		{
			// current level of loss - loss rate of names k
			double LossRates_inf = its_int_lossrates[k-1] ;
			double LossRates_delta = its_dbl_lossrates[k-1] - LossRates_inf;

			ind_lastloss = its_ind_Loss - its_int_lossrates[k-1] ;

			// barrier
			tmp_barrier=(its_barrier[k-1]-its_unique_beta[k-1]*x)/sqrt(1.-its_unique_beta[k-1]*its_unique_beta[k-1]);

			// conditional default probability
			// Cumulative Student distribution
			pk = NAG_prob_students_t( tmp_barrier, itsFreedomDegree);

			// its_ind_Loss: Loss Distribution for a basket of size (k-1) knowing factor value its_ind_x
			// tmp1

			// Loss Distribution for a basket of size k knowing factor value its_ind_x
			tmp2	=	tmp1 * (1.-pk);
			
			if (ind_lastloss >= 0)
			{
				tmp2	+=	pk * its_ProbCond->Elt(k-1,its_ind_x,ind_lastloss)*(1.-LossRates_delta);

				if (ind_lastloss >= 0)
				{
					tmp2	+=	pk * its_ProbCond->Elt(k-1,its_ind_x,ind_lastloss-1)*(LossRates_delta);
				}
			}
			// else: the loss of name k is too big

			its_ProbCond->SetElt(k,its_ind_x,its_ind_Loss, tmp2);

			tmp1 = tmp2;
		}
	}

	return tmp1;
}


void
ICM_Gauss1FLossDistrib_LJ::compute_distrib(const int& lup,const int& ldown)
{
	double tmpPrbLoss=0.,cumul_distrib=0.;
	double tmpPrbLoss_down=0.,cumul_distrib_down=0.;
	int min_lup =0,max_lup =0;
	int l = 0;

	// ---------------------
	its_lup	=	0;
//	max_lup = MAX(lup,its_lup);
	
	// ----------------------
	// Parameters
	AddressVector	TheParameterVector;

	TheParameterVector.Append(this);
	TheParameterVector.Append(&its_ind_x);
	TheParameterVector.Append(&its_ind_Loss);

	ICM_QMatrix<double> COEFS_up(its_nbnames,itsIntegrationStep,-999.);
	TheParameterVector.Append(&COEFS_up);

	ICM_QMatrix<double> COEFS_down(its_nbnames,itsIntegrationStep,-999.);
	TheParameterVector.Append(&COEFS_down);

	int nointegrationcoef=-1;
	TheParameterVector.Append(&nointegrationcoef);

	switch (itsCopulaType)
	{
		case qNO_COPULA:	
		case qGAUSSIAN:

			for (l=its_lup;l<=ldown;l++) 
			{
				if (l ==0)
				{
					its_ind_x=-1;

					TheIntegrator.Integrate(-6.0, 6.0, compute_cond_distribZeroLoss_Gaussian_Smile, &TheParameterVector, tmpPrbLoss, tmpPrbLoss_down); 

					its_lossdistrib[0] = tmpPrbLoss;
					its_lossdistrib_Down[0] = tmpPrbLoss_down;
					its_ind_Loss=1;
				}
				else
				{
					its_ind_x=-1;

					TheIntegrator.Integrate(-6.0, 6.0, compute_cond_distrib_Gaussian_Smile, &TheParameterVector, tmpPrbLoss, tmpPrbLoss_down); 

					its_lossdistrib[l] = tmpPrbLoss;
					its_lossdistrib_Down[l] = tmpPrbLoss_down;

					its_ind_Loss+=1;
				}
			}

			COEFS_down.RAZ(-999.);

			for (l=ldown+1;l<=lup;l++) 
			{
				its_ind_x=-1;

				TheIntegrator.Integrate(-6.0, 6.0, compute_cond_distrib_Gaussian_Integrator, &TheParameterVector, tmpPrbLoss); 

				its_lossdistrib[l] = tmpPrbLoss;

				its_ind_Loss+=1;
			}
			break;

		case qSTUDENT:

			break;
	}

	its_ind_Loss--;

	for (l=0;l<=ldown;l++)
	{
		cumul_distrib += its_lossdistrib[l];
		cumul_distrib_down += its_lossdistrib_Down[l];
	}
	for (l=ldown+1;l<=lup;l++)
	{
		cumul_distrib += its_lossdistrib[l];
	}

	its_taildistrib			=	1.0 - cumul_distrib;
	its_taildistrib_Down	=	1.0 - cumul_distrib_down;

}


void
ICM_Gauss1FLossDistrib_LJ::compute_distrib_Gauss_Legendre(const int& lup)
{ 
	double tmpPrbLoss=0.,cumul_distrib=0.;
	int min_lup =0,max_lup =0;
	int l = 0;

	// ------------------------
	its_lup		=	0;
	max_lup		=	lup;
//	max_lup = MAX(lup,its_lup);
	// ------------------------

	// ------------------------
	// Parameters
	AddressVector	TheParameterVector;

	TheParameterVector.Append(this);
	TheParameterVector.Append(&its_ind_x);
	TheParameterVector.Append(&its_ind_Loss);
	// ------------------------

	ICM_QMatrix<double> COEFS(its_nbnames,itsIntegrationStep,-999.);
	TheParameterVector.Append(&COEFS);

	int nointegrationcoef=-1;
	TheParameterVector.Append(&nointegrationcoef);

	switch (itsCopulaType)
	{
		case qNO_COPULA:	
		case qGAUSSIAN:

			for (l=its_lup;l<=max_lup;l++) 
			{
				if (l ==0)
				{
					its_ind_x=-1;
					TheIntegrator.Integrate(-6.0, 6.0, compute_cond_distribZeroLoss_Gaussian_Integrator, &TheParameterVector, tmpPrbLoss); 

					its_lossdistrib[0] = tmpPrbLoss;
					its_ind_Loss=1;

				}
				else
				{
					its_ind_x=-1;
					TheIntegrator.Integrate(-6.0, 6.0, compute_cond_distrib_Gaussian_Integrator, &TheParameterVector, tmpPrbLoss); 
					its_lossdistrib[l] = tmpPrbLoss;

					its_ind_Loss+=1;
				}
			}

			break;

		case qSTUDENT:
			break;

	}

	its_ind_Loss--;

	for (l=0;l<=lup;l++)
		cumul_distrib += its_lossdistrib[l];

	its_taildistrib			=	1.0 - cumul_distrib;

}


void
ICM_Gauss1FLossDistrib_LJ::compute_barrier()
{
	double Limit_case_Minus = -10.;
	double Limit_case_Plus	= 10.;

	int k;

	its_barrier.resize(its_nbnames);

	switch (itsCopulaType)
	{
	case qNO_COPULA:	
	case qGAUSSIAN:
		for (k=0;k<its_nbnames;k++)
		{
			if (fabs(its_pdef_at_maturity[k]) < DB_TOL)
				its_barrier[k] = Limit_case_Minus;
			else if (fabs(its_pdef_at_maturity[k]-1.0) < DB_TOL)
				its_barrier[k] = Limit_case_Plus;
			else
				its_barrier[k] =	NAG_deviates_normal_dist(its_pdef_at_maturity[k]);
		}
		break;
	case qSTUDENT:
		for (k=0;k<its_nbnames;k++)
			// Inverse Cumulative Student: nag_t_deviate
			its_barrier[k]	=	NAG_deviates_students_t(its_pdef_at_maturity[k], itsFreedomDegree );
		break;
	}

	//---------------------------------------------------------------
	// Term structure review
	//---------------------------------------------------------------
	if (IsTSR())
	{
		its_barrier_t2.resize(its_nbnames);

		switch (itsCopulaType)
		{
		case qNO_COPULA:	
		case qGAUSSIAN:
		for (k=0;k<its_nbnames;k++)
		{
			if (fabs(its_pdef_at_maturity_t2[k]) < DB_TOL)
				its_barrier_t2[k] = Limit_case_Minus;
			else if (fabs(its_pdef_at_maturity_t2[k]-1.0) < DB_TOL)
				its_barrier_t2[k] = Limit_case_Plus;
			else
				its_barrier_t2[k] =	NAG_deviates_normal_dist(its_pdef_at_maturity_t2[k]);
		}
		break;
		case qSTUDENT:
		for (k=0;k<its_nbnames;k++)
			// Inverse Cumulative Student: nag_t_deviate
			its_barrier_t2[k]	=	NAG_deviates_students_t(its_pdef_at_maturity_t2[k], itsFreedomDegree );
		break;
		}

	}
}


void ICM_Gauss1FLossDistrib_LJ::UpdateSmileCorrelation(const double& tranche_down, 
													   const double& tranche_up, 
													   const bool&	HedgeFlag)
{
	int i = 0;
	double	base_correl_down = 0.;
	double	base_correl_up = 0.;

	// maybe to be generalized with vectors of Betas.
	// for each Credit, linear interpolate both Base Correlation Down and Up

	for (i=0; i< its_nbnames; i++)
	{			
		double* V = itsBetasMatrix->DupRowAsTab(i);
		base_correl_down	=	linInterpol(&(*itsStrikes.begin()), itsNbStrikes, tranche_down, V);
		if (base_correl_down == K_HUGE_DOUBLE)
		{
			throw Exception(__LINE__, __FILE__, ERR_INVALID_DATA,
				"ERROR: unable to compute Base Correlation for lower bound!");
		}

		base_correl_up	=	linInterpol(&(*itsStrikes.begin()), itsNbStrikes, tranche_up, V);
		if (base_correl_up == K_HUGE_DOUBLE)
		{
			throw Exception(__LINE__, __FILE__, ERR_INVALID_DATA,
				"ERROR: unable to compute Base Correlation for upper bound!");
		}

		if (V) delete V;

		its_unique_beta[i] = base_correl_up;
		its_beta_Down[i] = base_correl_down;
	}

	precompute_coeffs();

	if (HedgeFlag)
		precompute_coeffs_perturb();
}



double
ICM_Gauss1FLossDistrib_LJ::compute_expectedlosstranche(const double& tranche_up, 
												  const double& tranche_down, 
												  const double& lossunit,
												  // ICM_QMatrix<double>* ShiftMatrix,
												  // const int& Tenor,
												  // const int& Issuer,
												  vector<double>& losses,
												  vector<double>& ProbasLossesUp,
												  vector<double>& ProbasLossesDown)
{
	double exp_loss_tranche_down = 0.;
	double exp_loss_tranche_up = 0.;
	
	int lup=floor(tranche_up/lossunit);
	int ldown=floor(tranche_down/lossunit);

	int i=0;	
	
	its_ProbCond->ResizeWithCopy(lup+1,0.);
	its_ProbCond_Down->ResizeWithCopy(ldown+1,0.);

	its_lossdistrib.resize(lup+1);
	its_lossdistrib_Down.resize(ldown+1);
	// Compute the Distribution (by default, ldown = 0)
	if (tranche_down)
		compute_distrib(lup, ldown);
	else
		compute_distrib(lup);

	ProbasLossesUp.resize(lup+1,0.);
	ProbasLossesDown.resize(lup+1,0.);
	
	int l	=	0;
	double	loss_level;
	losses.clear();

	// first loop, deal with beta_down and up
	ProbasLossesDown[0]=its_lossdistrib_Down[0];
	ProbasLossesUp[0]=its_lossdistrib[0];

	for (l=1; l<=ldown; l++) 
	{
		loss_level	=	l*lossunit;
		exp_loss_tranche_down	+=	its_lossdistrib_Down[l] * loss_level;		
		exp_loss_tranche_up		+=	its_lossdistrib[l] * loss_level;
		ProbasLossesDown[l] = its_lossdistrib_Down[l];
		ProbasLossesUp[l]  = its_lossdistrib[l];
	}
		
	// second loop, deal only with beat_up
	for (l=ldown+1;l<=lup;l++)
	{
		loss_level	=	l*lossunit;
		exp_loss_tranche_up		+=	its_lossdistrib[l] * loss_level;
		ProbasLossesUp[l] = its_lossdistrib[l];
	}

	for (l=0;l<its_lossdistrib.size();l++)
	{	losses.push_back(its_lossdistrib[l]); }

	// add the tail
	if (tranche_down)
		exp_loss_tranche_down	+=	tranche_down * its_taildistrib_Down;
	exp_loss_tranche_up		+=	tranche_up * its_taildistrib;

	return (exp_loss_tranche_up - exp_loss_tranche_down);
}


void ICM_Gauss1FLossDistrib_LJ :: precompute_coeffs()
{
	int k;
	double	den, beta_value;

	//term structure review
	if (IsTSR())
		{precompute_coeffs_TSR();}

	for (k=0; k<its_nbnames; k++)
	{
		if (fabs(beta_value = its_unique_beta[k]) == 1.0)
		{
			throw Exception(__LINE__, __FILE__, ERR_INVALID_DATA,
				"Beta values must be strictly between -1.0 and 1.0");
		}
		
		den = 1.0 / sqrt(1.0 - beta_value * beta_value);
		
		its_coeff_a[k]	=	den;
		its_coeff_b[k]	=	den;
		
		its_coeff_b[k]	*=	beta_value;
		its_coeff_a[k]	*=	its_barrier[k];


		if (fabs(beta_value = its_beta_Down[k]) == 1.0)
		{
			throw Exception(__LINE__, __FILE__, ERR_INVALID_DATA,
				"Beta values must be strictly between -1.0 and 1.0");
		}
		
		den = 1.0 / sqrt(1.0 - beta_value * beta_value);
		
		its_coeff_a_down[k]	=	den;
		its_coeff_b_down[k]	=	den;
		
		its_coeff_b_down[k]	*=	beta_value;
		its_coeff_a_down[k]	*=	its_barrier[k];
	}

	//case forward collateral in TS mode
	if (GetFwdCollatTScase())
	{
		for (k=0; k<its_nbnames; k++)
		{
		if (fabs(beta_value = its_collat_fwd_unique_beta[k]) == 1.0)
		{
			throw Exception(__LINE__, __FILE__, ERR_INVALID_DATA,
				"Beta values must be strictly between -1.0 and 1.0");
		}
		
		den = 1.0 / sqrt(1.0 - beta_value * beta_value);
		
		its_collat_fwd_coeff_a[k]	=	den;
		its_collat_fwd_coeff_b[k]	=	den;
		
		its_collat_fwd_coeff_b[k]	*=	beta_value;
		its_collat_fwd_coeff_a[k]	*=	its_barrier[k];


		if (fabs(beta_value = its_collat_fwd_beta_Down[k]) == 1.0)
		{
			throw Exception(__LINE__, __FILE__, ERR_INVALID_DATA,
				"Beta values must be strictly between -1.0 and 1.0");
		}
		
		den = 1.0 / sqrt(1.0 - beta_value * beta_value);
		
		its_collat_fwd_coeff_a_down[k]	=	den;
		its_collat_fwd_coeff_b_down[k]	=	den;
		
		its_collat_fwd_coeff_b_down[k]	*=	beta_value;
		its_collat_fwd_coeff_a_down[k]	*=	its_barrier[k];
		}
	}
}


void ICM_Gauss1FLossDistrib_LJ :: precompute_coeffs_perturb()
{
	int k;
	double	den, beta_value;

	//term structure review
	if (IsTSR())
		{precompute_coeffs_perturb_TSR();}

	if (GetCtxt()){
	if ((GetCtxt()->m_IsUsed) && (GetCtxt()->m_DistribType==qDISTRIB_STD_TSR))
	
	{
	for (k=0; k<its_nbnames; k++)
	{
		if (fabs(beta_value = GetCtxt()->m_ts_beta_up_maturity_t1[k]) == 1.0)
		{
			throw Exception(__LINE__, __FILE__, ERR_INVALID_DATA,
				"Beta values must be strictly between -1.0 and 1.0");
		}
		
		den = 1.0 / sqrt(1.0 - beta_value * beta_value);
		
		its_coeff_a_perturb[k]	=	den;
		its_coeff_a_perturb[k]	*=	its_barrier_perturb[k];

		if (its_coeff_b_perturb.size()>0){
		its_coeff_b_perturb[k]	=	den;
		its_coeff_b_perturb[k]	*=	beta_value;
		}

		if (fabs(beta_value = GetCtxt()->m_ts_beta_dw_maturity_t1[k]) == 1.0)
		{
			throw Exception(__LINE__, __FILE__, ERR_INVALID_DATA,
				"Beta values must be strictly between -1.0 and 1.0");
		}
		
		den = 1.0 / sqrt(1.0 - beta_value * beta_value);
		
		its_coeff_a_down_perturb[k]	=	den;
		its_coeff_a_down_perturb[k]	*=	its_barrier_perturb[k];

		if (its_coeff_b_down_perturb.size()>0){
		its_coeff_b_down_perturb[k]	=	den;
		its_coeff_b_down_perturb[k]	*=	beta_value;
		}
	} 
	
	}}
	else
	{
	for (k=0; k<its_nbnames; k++)
	{
		if (fabs(beta_value = its_unique_beta[k]) == 1.0)
		{
			throw Exception(__LINE__, __FILE__, ERR_INVALID_DATA,
				"Beta values must be strictly between -1.0 and 1.0");
		}
		
		den = 1.0 / sqrt(1.0 - beta_value * beta_value);
		
		its_coeff_a_perturb[k]	=	den;
		its_coeff_a_perturb[k]	*=	its_barrier_perturb[k];

		if (its_coeff_b_perturb.size()>0){
		its_coeff_b_perturb[k]	=	den;
		its_coeff_b_perturb[k]	*=	beta_value;
		}

		if (fabs(beta_value = its_beta_Down[k]) == 1.0)
		{
			throw Exception(__LINE__, __FILE__, ERR_INVALID_DATA,
				"Beta values must be strictly between -1.0 and 1.0");
		}
		
		den = 1.0 / sqrt(1.0 - beta_value * beta_value);
		
		its_coeff_a_down_perturb[k]	=	den;
		its_coeff_a_down_perturb[k]	*=	its_barrier_perturb[k];

		if (its_coeff_b_down_perturb.size()>0){
		its_coeff_b_down_perturb[k]	=	den;
		its_coeff_b_down_perturb[k]	*=	beta_value;
		}
	} 
	}

	//case forward collateral in TS mode
	if (GetFwdCollatTScase())
	{
		for (k=0; k<its_nbnames; k++)
		{	
		if (fabs(beta_value = its_collat_fwd_unique_beta[k]) == 1.0)
		{
			throw Exception(__LINE__, __FILE__, ERR_INVALID_DATA,
				"Beta values must be strictly between -1.0 and 1.0");
		}
		
		den = 1.0 / sqrt(1.0 - beta_value * beta_value);
		
		its_collat_fwd_coeff_a_perturb[k]	=	den;
		its_collat_fwd_coeff_a_perturb[k]	*=	its_barrier_perturb[k];


		if (fabs(beta_value = its_collat_fwd_beta_Down[k]) == 1.0)
		{
			throw Exception(__LINE__, __FILE__, ERR_INVALID_DATA,
				"Beta values must be strictly between -1.0 and 1.0");
		}
		
		den = 1.0 / sqrt(1.0 - beta_value * beta_value);
		
		its_collat_fwd_coeff_a_down_perturb[k]	=	den;
		its_collat_fwd_coeff_a_down_perturb[k]	*=	its_barrier_perturb[k];

		}
	} 

	if (IsTSR() && GetFwdCollatTScase()) 
	{
		for (k=0; k<its_nbnames; k++)
		{	
		if (fabs(beta_value = its_collat_fwd_unique_beta[k]) == 1.0)
		{
			throw Exception(__LINE__, __FILE__, ERR_INVALID_DATA,
				"Beta values must be strictly between -1.0 and 1.0");
		}
		
		den = 1.0 / sqrt(1.0 - beta_value * beta_value);
		
		its_coeff_a_perturb[k]	=	den;
		its_coeff_a_perturb[k]	*=	its_barrier_perturb[k];

		its_coeff_b_perturb[k]	=	den;
		its_coeff_b_perturb[k]	*=	beta_value;

		if (fabs(beta_value = its_collat_fwd_beta_Down[k]) == 1.0)
		{
			throw Exception(__LINE__, __FILE__, ERR_INVALID_DATA,
				"Beta values must be strictly between -1.0 and 1.0");
		}
		
		den = 1.0 / sqrt(1.0 - beta_value * beta_value);
		
		its_coeff_a_down_perturb[k]	=	den;
		its_coeff_a_down_perturb[k]	*=	its_barrier_perturb[k];

		its_coeff_b_down_perturb[k]	=	den;
		its_coeff_b_down_perturb[k]	*=	beta_value;

		}
	}

}


void ICM_Gauss1FLossDistrib_LJ :: SetIntegratorType(const qIntegratorChoice&	TheIntegratorType, 
													const int&	TheStep)
{
	itsIntegrationMethod	=	TheIntegratorType;
	itsIntegrationStep		=	TheStep;

	TheIntegrator.SetIntegrationType(TheIntegratorType);
	
	// Rajout du || sur la méthode qTrapeze afin que le step de l'integrator soit setter a la bonne valeur
	if ((TheIntegratorType == qGAUSS_LEGENDRE) || (TheIntegratorType == qGAUSS_HERMITE) || (TheIntegratorType == qTRAPEZE))
		TheIntegrator.SetIntegrationStep(TheStep);
}




void ICM_Gauss1FLossDistrib_LJ::UpdateCorrelation(const double& beta_down, 
												  const double& beta_up,
												  const bool& HedgeFlag)
{
	int i;

	// maybe to be generalized with vectors of Betas.
	// for each Credit, linear interpolate both Base Correlation Down and Up

	for (i=0; i< its_nbnames; i++)
	{			
		its_unique_beta[i] = beta_up;
		its_beta_Down[i] = beta_down;
	}

	precompute_coeffs();

	if (HedgeFlag)
		precompute_coeffs_perturb();
}


void ICM_Gauss1FLossDistrib_LJ::UpdateCorrelation(const ARM_Vector& V_beta_down,
												  const ARM_Vector&V_beta_up, 
												  bool HedgeFlag)
{
	int i = 0;

	// maybe to be generalized with vectors of Betas.
	// for each Credit, linear interpolate both Base Correlation Down and Up

	for (i=0; i< its_nbnames; i++)
	{			
		its_unique_beta[i] = V_beta_up[i];
		its_beta_Down[i] = V_beta_down[i];
	}

	precompute_coeffs();

	if (HedgeFlag)
		precompute_coeffs_perturb();
}

//----------------------------------------------//
//					View	
//----------------------------------------------//
void ICM_Gauss1FLossDistrib_LJ::View(char* id, FILE* ficOut)
{	
	FILE* fOut;
	char  fOutName[200];

	if ( ficOut == NULL )
	{
	ARM_GetViewFile(ARM_VIEW_FILE, id, fOutName);
	
	   (void) unlink(fOutName);

       fOut = fopen(fOutName, "w"); 
    }
	else
	{
	fOut = ficOut;
	} 

	int size =0; 

    fprintf(fOut, "\t\t\t ----------------- Homogeneous Fast Loss Dsitrib Gaussian  ----------------- \n");
	int ind_loss=0;
	for (ind_loss=0; ind_loss<its_lossdistrib_Down.size(); ind_loss++)
		fprintf(fOut,"ind_x : %i\t%f\t%f\n",ind_loss, its_lossdistrib[ind_loss], its_lossdistrib_Down[ind_loss]);
	for (ind_loss=its_lossdistrib_Down.size(); ind_loss<its_lossdistrib.size(); ind_loss++)
	{
		if (ind_loss==its_lossdistrib_Down.size())
			fprintf(fOut,"ind_x : %i\t%f\t%f\n",ind_loss, its_lossdistrib[ind_loss], its_taildistrib_Down);
		else
			fprintf(fOut,"ind_x : %i\t%f\n",ind_loss, its_lossdistrib[ind_loss]);
	}
	fprintf(fOut,"tail : \t%f\n",its_taildistrib);

	/*fprintf(fOut, "\t\t\t ----------------- Prob Cond  ----------------- \n");
	for (int ind_x=0; ind_x<TheIntegrator.GetIntegrationStep(); ind_x++)
	{
		for (ind_loss=0; ind_loss<its_lossdistrib.size(); ind_loss++)
		{
			for (int ind_name=0; ind_name<=this->its_nbnames; ind_name++)
				fprintf(fOut,"%f\t",its_ProbCond->Elt(ind_name,ind_x,ind_loss));
			fprintf(fOut,"\n");
		}
		
		fprintf(fOut,"\n\n");
	}*/
	
	
	
	
	ICM_Distribution::View(id, fOut);

	if ( ficOut == NULL )
	{
		fclose(fOut);
	}
}


//--------------------------------------------------------------------------------------
// Methode Set cadre Term Structure Review
//--------------------------------------------------------------------------------------
ICM_Gauss1FLossDistrib_LJ::ICM_Gauss1FLossDistrib_LJ(const double& t1_corr,
							const double& t2_corr,
							const int& nbnames,
							const ARM_Vector& pdefault_maturity_up_t1,
							const ARM_Vector& pdefault_maturity_up_t2,
							const ARM_Vector& pdefault_maturity_dw_t1,
							const ARM_Vector& pdefault_maturity_dw_t2,
							const ARM_Vector& beta_t1,
							const ARM_Vector& beta_t2,
							const ARM_Vector & LossRates,
							const int& discretizationstep ,
							const int& CopulaType ,
							const qIntegratorChoice&	IntegrationMethod ,
							const int& IntegrationStep)
{
	Init();

	Set(t1_corr,t2_corr,nbnames,pdefault_maturity_up_t1,pdefault_maturity_up_t2,pdefault_maturity_dw_t1,pdefault_maturity_dw_t2,beta_t1,beta_t2,LossRates,discretizationstep,CopulaType,IntegrationMethod,IntegrationStep);

}

//--------------------------------------------------------------------------------------
// Methode Set cadre Term Structure Review
//--------------------------------------------------------------------------------------
void ICM_Gauss1FLossDistrib_LJ::Set(const double& t1_corr,
							const double& t2_corr,
							const int& nbnames,
							const ARM_Vector&  pdefault_maturity_up_t1,
							const ARM_Vector&  pdefault_maturity_up_t2,
							const ARM_Vector&  pdefault_maturity_dw_t1,
							const ARM_Vector&  pdefault_maturity_dw_t2,
							const ARM_Vector&  beta_t1,
							const ARM_Vector&  beta_t2,
							const ARM_Vector & LossRates,
							const int& discretizationstep ,
							const int& CopulaType ,
							const qIntegratorChoice&	IntegrationMethod ,
							const int& IntegrationStep)
{
	SetTSR(true);
	SetT1(t1_corr);
	SetT2(t2_corr);

	int size;
	if (discretizationstep	==	0) size	=	20;
	else size	=	discretizationstep;

	if (nbnames != GetNbNames())
	{
		ICM_Distribution::Set(nbnames);
	
		if (its_lossdistrib_perturb_Down) delete its_lossdistrib_perturb_Down;
		its_lossdistrib_perturb_Down= new ICM_QMatrix<double>(1,nbnames);

		its_taildistrib_perturb_Down.resize(nbnames);

		its_beta_Down.resize(nbnames);
		its_coeff_a.resize(nbnames);
		its_coeff_b.resize(nbnames);
		its_coeff_a_down.resize(nbnames);
		its_coeff_b_down.resize(nbnames);

		//term structure review
		its_beta_Down_t2.resize(nbnames);
		its_coeff_a_t2.resize(nbnames);
		its_coeff_b_t2.resize(nbnames);
		its_coeff_a_down_t2.resize(nbnames);
		its_coeff_b_down_t2.resize(nbnames);
	}

	if ((nbnames != its_nbnames) || (itsIntegrationStep != size) || 
		// 17783 (!its_ProbCond)||(!its_ProbCond_Perturb)||(!its_ProbCond_Perturb_Down)||(!its_ProbCond_Down))
		(!its_ProbCond)||(!its_ProbCond_Perturb_Down)||(!its_ProbCond_Down))
	{
		if (its_ProbCond) delete its_ProbCond;
		its_ProbCond= new ICM_QCubix<double>(nbnames+1,size,1,0.);

		// 17783 if (its_ProbCond_Perturb) delete its_ProbCond_Perturb;
		// 17783 its_ProbCond_Perturb= new ICM_QCubix<double>(nbnames+1,size,1,0.);

		if (its_ProbCond_Perturb_Down) delete its_ProbCond_Perturb_Down;
		its_ProbCond_Perturb_Down= new ICM_QCubix<double>(nbnames+1,size,1,0.);

		if (its_ProbCond_Down) delete its_ProbCond_Down;
		its_ProbCond_Down= new ICM_QCubix<double>(nbnames+1,size,1,0.);
	}

	SetNbNames(nbnames);

	itsCopulaType = CopulaType;
	
	itsIntegrationStep = IntegrationStep;
	itsIntegrationMethod = IntegrationMethod;
	// Si le nombre de pas et la méthode utilisés ne sont pas compatibles, la parité du nombre de pas est prioritaire
	if ( (itsIntegrationStep % 2 == 0) && (itsIntegrationMethod == qGAUSS_LEGENDRE))
		itsIntegrationMethod = qGAUSS_HERMITE;
	else if ( (itsIntegrationStep % 2 == 1) && (itsIntegrationMethod == qGAUSS_HERMITE))
		itsIntegrationMethod = qGAUSS_LEGENDRE;
	
	SetUniqueBeta(beta_t1); 
	SetPdefAtMaturity(pdefault_maturity_up_t1); 
	SetPdefPerturb(pdefault_maturity_dw_t1); 
	compute_min_pdef(GetPdefAtMaturity());

	//term structure review
	SetUniqueBeta_t2(beta_t2); 
	SetPdefAtMaturity_t2(pdefault_maturity_up_t1); 
	SetPdefPerturb_t2(pdefault_maturity_dw_t2); 

	SetIntLossRates(LossRates);
	check_homogeneous(LossRates);
	compute_barrier();		// Copula dependant

	its_lossdistrib_Down.resize(1);
}

//--------------------------------------------------------------------------------------
// Methode UpdateCorrelation_TSR cadre Term Structure Review
//--------------------------------------------------------------------------------------
void ICM_Gauss1FLossDistrib_LJ::UpdateCorrelation_TSR(const ARM_Vector & V_beta_down_t1, 
													const ARM_Vector & V_beta_down_t2, 
													const ARM_Vector & V_beta_up_t1, 
													const ARM_Vector & V_beta_up_t2)
{
	int i = 0;

	// maybe to be generalized with vectors of Betas.
	// for each Credit, linear interpolate both Base Correlation Down and Up

	for (i=0; i< its_nbnames; i++)
	{			
		its_unique_beta[i] = V_beta_up_t1[i];
		its_beta_Down[i] = V_beta_down_t1[i];
		its_unique_beta_t2[i] = V_beta_up_t2[i];
		its_beta_Down_t2[i] = V_beta_down_t2[i];
	}

	precompute_coeffs();

}

//--------------------------------------------------------------------------------------
// Methode precompute_coeffs_TSR cadre Term Structure Review
//--------------------------------------------------------------------------------------

void ICM_Gauss1FLossDistrib_LJ :: precompute_coeffs_TSR()
{
	int k;
	double	den, beta_value;

	for (k=0; k<its_nbnames; k++)
	{
		if (fabs(beta_value = its_unique_beta_t2[k]) == 1.0)
		{
			throw Exception(__LINE__, __FILE__, ERR_INVALID_DATA,
				"Beta values must be strictly between -1.0 and 1.0");
		}
		
		den = 1.0 / sqrt(1.0 - beta_value * beta_value);
		
		its_coeff_a_t2[k]	=	den;
		its_coeff_b_t2[k]	=	den;
		
		its_coeff_b_t2[k]	*=	beta_value;
		its_coeff_a_t2[k]	*=	its_barrier_t2[k];


		if (fabs(beta_value = its_beta_Down_t2[k]) == 1.0)
		{
			throw Exception(__LINE__, __FILE__, ERR_INVALID_DATA,
				"Beta values must be strictly between -1.0 and 1.0");
		}
		
		den = 1.0 / sqrt(1.0 - beta_value * beta_value);
		
		its_coeff_a_down_t2[k]	=	den;
		its_coeff_b_down_t2[k]	=	den;
		
		its_coeff_b_down_t2[k]	*=	beta_value;
		its_coeff_a_down_t2[k]	*=	its_barrier_t2[k];
	}
}



void ICM_Gauss1FLossDistrib_LJ :: precompute_coeffs_perturb_TSR()
{
	int k;
	double	den, beta_value;

	if (GetFwdCollatTScase())
	{
	//case forward collateral in TS mode
	for (k=0; k<its_nbnames; k++)
	{	
		if (fabs(beta_value = its_collat_fwd_unique_beta_t2[k]) == 1.0)
		{
			throw Exception(__LINE__, __FILE__, ERR_INVALID_DATA,
				"Beta values must be strictly between -1.0 and 1.0");
		}
		
		den = 1.0 / sqrt(1.0 - beta_value * beta_value);
		
		its_coeff_a_perturb_t2[k]	=	den;
		its_coeff_a_perturb_t2[k]	*=	its_barrier_perturb_t2[k];

		its_coeff_b_perturb_t2[k]	=	den;
		its_coeff_b_perturb_t2[k]	*=	beta_value;

		if (fabs(beta_value = its_collat_fwd_beta_Down_t2[k]) == 1.0)
		{
			throw Exception(__LINE__, __FILE__, ERR_INVALID_DATA,
				"Beta values must be strictly between -1.0 and 1.0");
		}
		
		den = 1.0 / sqrt(1.0 - beta_value * beta_value);
		
		its_coeff_a_down_perturb_t2[k]	=	den;
		its_coeff_a_down_perturb_t2[k]	*=	its_barrier_perturb_t2[k];

		its_coeff_b_down_perturb_t2[k]	=	den;
		its_coeff_b_down_perturb_t2[k]	*=	beta_value;
	}
	}
	// else if ((GetCtxt()->m_IsUsed) && (GetCtxt()->m_DistribType==qDISTRIB_STD_TSR))
	else if (GetCtxt()!=NULL && (GetCtxt()->m_IsUsed) && (GetCtxt()->m_DistribType==qDISTRIB_STD_TSR))
	{
	//case forward collateral in TS mode
	for (k=0; k<its_nbnames; k++)
	{	
		if (fabs(beta_value = GetCtxt()->m_ts_beta_up_maturity_t2[k]) == 1.0)
		{
			throw Exception(__LINE__, __FILE__, ERR_INVALID_DATA,
				"Beta values must be strictly between -1.0 and 1.0");
		}
		
		den = 1.0 / sqrt(1.0 - beta_value * beta_value);
		
		its_coeff_a_perturb_t2[k]	=	den;
		its_coeff_a_perturb_t2[k]	*=	its_barrier_perturb_t2[k];

		its_coeff_b_perturb_t2[k]	=	den;
		its_coeff_b_perturb_t2[k]	*=	beta_value;

		if (fabs(beta_value = GetCtxt()->m_ts_beta_dw_maturity_t2[k]) == 1.0)
		{
			throw Exception(__LINE__, __FILE__, ERR_INVALID_DATA,
				"Beta values must be strictly between -1.0 and 1.0");
		}
		
		den = 1.0 / sqrt(1.0 - beta_value * beta_value);
		
		its_coeff_a_down_perturb_t2[k]	=	den;
		its_coeff_a_down_perturb_t2[k]	*=	its_barrier_perturb_t2[k];

		its_coeff_b_down_perturb_t2[k]	=	den;
		its_coeff_b_down_perturb_t2[k]	*=	beta_value;
	}
	}

}

// ---------------------------------------------------------------------
// Generic constructor
// ---------------------------------------------------------------------
ICM_Gauss1FLossDistrib_LJ::ICM_Gauss1FLossDistrib_LJ(CtxtDistrib* c)
{
	Init();

	Set(c);
}

void ICM_Gauss1FLossDistrib_LJ::Set(CtxtDistrib* c)
{

	double output=0.;

	SetCtxt(c);

	switch (c->m_DistribType)
	{
		case qDISTRIB_STD_NEW :
			{
				Set_STD(c);
				break;
			}
		case qDISTRIB_APPROX :
			{
				Set_APPROX(c);
				break;
			}
		case qDISTRIB_VN :
			{
				Set_VN(c);
				break;
			}
		case qDISTRIB_COLLFWD_TSR:
		case qDISTRIB_COLLFWD:
		case qDISTRIB_STD_TSR:
			{
				Set(c->m_ts_t1_corr,
					c->m_ts_t2_corr,
					c->m_nbnames,
					c->m_ts_pdef_up_maturity_t1,
					c->m_ts_pdef_up_maturity_t2,
					c->m_ts_pdef_dw_maturity_t1,
					c->m_ts_pdef_dw_maturity_t2,
					c->m_ts_beta_up_maturity_t1,
					c->m_ts_beta_up_maturity_t2,
					c->m_LossRates, 
					c->m_IntegrationStep,
					c->m_CopulaType,
					c->m_IntegrationMethod,
					c->m_IntegrationStep);
				break;
			}
		default :
		case qDISTRIB_STD:
			{
				Set(c->m_nbnames,
					c->m_pdef_maturity,
					c->m_beta_up_maturity,
					c->m_LossRates, 
					c->m_IntegrationStep,
					c->m_CopulaType,
					c->m_IntegrationMethod,
					c->m_IntegrationStep);
				break;
			}
	}
}


//-------------------------------------------------------------------
// Generic EL computation
//-------------------------------------------------------------------
double ICM_Gauss1FLossDistrib_LJ::ComputeEL(CtxtDistrib* c)
{
	double ELT = 0.;

	SetCtxt(c);
	GetCtxt()->ComputeAll();
	SetIntegrationMethod(c->m_IntegrationMethod);

	switch (c->m_DistribType)
	{
		case qDISTRIB_STD_NEW :
			{
				SetCopulaType(qGAUSSIAN);
				SetIntegratorType_STD(c->m_IntegrationMethod,c->m_IntegrationStep);
				SetIntegrationStep(c->m_IntegrationStep);

				ELT =compute_expectedlosstranche_STD(c); 

				break;
			}
		case qDISTRIB_APPROX :
			{
				SetCopulaType(qGAUSSIAN);
				SetIntegratorType_APPROX(c->m_IntegrationMethod,c->m_IntegrationStep);
				SetIntegrationStep(c->m_IntegrationStep);

				ELT =compute_expectedlosstranche_APPROX(c); 

				break;
			}
		case qDISTRIB_VN :
			{
				SetCopulaType(qGAUSSIAN);
				SetIntegratorType_VN(c->m_IntegrationMethod,c->m_IntegrationStep);
				SetIntegrationStep(c->m_IntegrationStep);

				ELT =compute_expectedlosstranche_VN(c); 

				break;
			}
		case qDISTRIB_COLLFWD_TSR:
		case qDISTRIB_COLLFWD:
		case qDISTRIB_STD_TSR:
			{
				UpdateCorrelation_TSR(c->m_ts_beta_dw_maturity_t1,c->m_ts_beta_dw_maturity_t2,
									  c->m_ts_beta_up_maturity_t1,c->m_ts_beta_up_maturity_t2);

				SetCopulaType(qGAUSSIAN);
				SetIntegratorType(c->m_IntegrationMethod,c->m_IntegrationStep);
				SetIntegrationStep(c->m_IntegrationStep);

				ELT =compute_expectedlosstranche(c->m_TrancheUp,c->m_TrancheDown,c->m_LossUnit,
															/**-1.,-1.,*/c->m_OutLosses,c->m_ProbasLossesUp,c->m_ProbasLossesDown); 

				break;
			}
		case qDISTRIB_STD:
		default :
			{
				UpdateCorrelation(c->m_beta_dw_maturity,c->m_beta_up_maturity);

				SetCopulaType(qGAUSSIAN);
				SetIntegratorType(c->m_IntegrationMethod,c->m_IntegrationStep);
				SetIntegrationStep(c->m_IntegrationStep);

				ELT =compute_expectedlosstranche(c->m_TrancheUp,c->m_TrancheDown,c->m_LossUnit,
															/**-1.,-1.,*/c->m_OutLosses,c->m_ProbasLossesUp,c->m_ProbasLossesDown); 

				break;
			}
	}

	return (ELT);
}

extern "C" void 
ZeroLoss_TSR_stepup_Integrator(void* Param, double x, double& res)
{
	// get parameters from stack
	ICM_Gauss1FLossDistrib_LJ*	TheModel;
	TheModel	=	(ICM_Gauss1FLossDistrib_LJ*)((*(AddressVector*)Param)[0]);

	res	=	TheModel->ZeroLoss_TSR_stepup(x,Param);
}


//-------------------------------------------------------------------
// Step Up Subordination
//-------------------------------------------------------------------
double ICM_Gauss1FLossDistrib_LJ::ComputeStepUp(CtxtDistrib* ctxt)
{
	SetCtxt(ctxt);

	vector<double> tmpPrbLoss_up;
	vector<double> tmpPrbLoss_down;
	double CumtmpPrbLoss=0.;
	double CumtmpPrbLoss_=0.;
	double Loss = 0.,Loss1 = 0.,Loss2 = 0.;
	int i=0;

	//-------------------------------------------------------------------
	// partie non step up
	//-------------------------------------------------------------------
	int lup=MIN(floor(GetCtxt()->m_TrancheUp/GetCtxt()->m_LossUnit),GetCtxt()->m_nbnames);
	int ldown=MIN(floor(GetCtxt()->m_TrancheDown/GetCtxt()->m_LossUnit),GetCtxt()->m_nbnames);

	double stepup = GetCtxt()->m_TrancheUp_stepup;
	double stepdw = GetCtxt()->m_TrancheDown_stepup;

	int lup_stepup=MIN(floor(GetCtxt()->m_TrancheUp_stepup/GetCtxt()->m_LossUnit),GetCtxt()->m_nbnames);
	int ldown_stepup=MIN(floor(GetCtxt()->m_TrancheDown_stepup/GetCtxt()->m_LossUnit),GetCtxt()->m_nbnames);

	double TrancheUp = GetCtxt()->m_TrancheUp;
	double TrancheDw = GetCtxt()->m_TrancheDown;
	const ARM_Vector&  Pdefault_t1 = GetCtxt()->m_ts_pdef_up_maturity_t1;
	const ARM_Vector&  Pdefault_t2 = GetCtxt()->m_ts_pdef_up_maturity_t2;
	const ARM_Vector&  Pdefault_down_t1 = GetCtxt()->m_ts_pdef_dw_maturity_t1;
	const ARM_Vector&  Pdefault_down_t2 = GetCtxt()->m_ts_pdef_dw_maturity_t2;

	GetCtxt()->m_TrancheUp = GetCtxt()->m_TrancheUp_stepup;
	GetCtxt()->m_TrancheDown = GetCtxt()->m_TrancheDown_stepup;

	const ARM_Vector& ts_beta_dw_maturity_t1 = GetCtxt()->m_ts_beta_dw_maturity_t1;
	const ARM_Vector& ts_beta_dw_maturity_t2 = GetCtxt()->m_ts_beta_dw_maturity_t2;
	const ARM_Vector& ts_beta_up_maturity_t1 = GetCtxt()->m_ts_beta_up_maturity_t1;
	const ARM_Vector& ts_beta_up_maturity_t2 = GetCtxt()->m_ts_beta_up_maturity_t2;

	GetCtxt()->m_ts_pdef_up_maturity_t1=GetCtxt()->m_ts_stepup_pdef_up_maturity_t1;
	GetCtxt()->m_ts_pdef_up_maturity_t2=GetCtxt()->m_ts_stepup_pdef_up_maturity_t2;
	GetCtxt()->m_ts_pdef_dw_maturity_t1=GetCtxt()->m_ts_stepup_pdef_dw_maturity_t1;
	GetCtxt()->m_ts_pdef_dw_maturity_t2=GetCtxt()->m_ts_stepup_pdef_dw_maturity_t2;

	GetCtxt()->m_ts_beta_dw_maturity_t1=GetCtxt()->m_ts_stepup_beta_dw_maturity_t1;
	GetCtxt()->m_ts_beta_dw_maturity_t2=GetCtxt()->m_ts_stepup_beta_dw_maturity_t2;
	GetCtxt()->m_ts_beta_up_maturity_t1=GetCtxt()->m_ts_stepup_beta_up_maturity_t1;
	GetCtxt()->m_ts_beta_up_maturity_t2=GetCtxt()->m_ts_stepup_beta_up_maturity_t2;

	GetCtxt()->ResetAll();
	GetCtxt()->ComputeALLBarriers();

	SetPdefAtMaturity(GetCtxt()->m_ts_pdef_up_maturity_t1); 
	SetPdefPerturb(GetCtxt()->m_ts_pdef_dw_maturity_t1); 
	SetPdefAtMaturity_t2(GetCtxt()->m_ts_pdef_up_maturity_t2); 
	SetPdefPerturb_t2(GetCtxt()->m_ts_pdef_dw_maturity_t2); 
	compute_barrier();
	compute_barrier_perturb();

	//calcul de l'expected loss
	ComputeEL(GetCtxt());
	tmpPrbLoss_up=GetCtxt()->m_ProbasLossesUp;
	tmpPrbLoss_down=GetCtxt()->m_ProbasLossesDown;

	GetCtxt()->m_TrancheUp = TrancheUp;
	GetCtxt()->m_TrancheDown = TrancheDw;

	GetCtxt()->m_ts_pdef_up_maturity_t1=Pdefault_t1;
	GetCtxt()->m_ts_pdef_up_maturity_t2=Pdefault_t2;
	GetCtxt()->m_ts_pdef_dw_maturity_t1=Pdefault_down_t1;
	GetCtxt()->m_ts_pdef_dw_maturity_t2=Pdefault_down_t2;

	GetCtxt()->m_ts_beta_dw_maturity_t1=ts_beta_dw_maturity_t1;
	GetCtxt()->m_ts_beta_dw_maturity_t2=ts_beta_dw_maturity_t2;
	GetCtxt()->m_ts_beta_up_maturity_t1=ts_beta_up_maturity_t1;
	GetCtxt()->m_ts_beta_up_maturity_t2=ts_beta_up_maturity_t2;

	GetCtxt()->ResetAll();
	GetCtxt()->ComputeAll();

	//-------------------------------------------------------------------
	// partie step up
	//-------------------------------------------------------------------
	double LossUp1 = 0.,LossDown1 = 0.,value=0.;
	vector<double> tmpPrbLoss_up_stepup;
	vector<double> tmpPrbLoss_down_stepup;
	double CumtmpPrbLossUp_stepup=0.;
	double CumtmpPrbLossDown_stepup=0.;

	if (its_ProbCond_stepup)
		delete its_ProbCond_stepup;
	its_ProbCond_stepup= new ICM_QCubix<double>(GetCtxt()->m_nbnames+1,GetCtxt()->m_nbnames+1,GetCtxt()->m_nbnames+1,0.);
	its_ProbCond_stepup->SetElt(0,0,0,1.);

	vector<double> CheckLoss1;
	vector<double> CheckLoss2;

	for (i=0;i<=lup_stepup;i++)
	{
		CumtmpPrbLoss+=tmpPrbLoss_up[i];
		if (i<=lup) {CumtmpPrbLoss_=CumtmpPrbLoss;}
	}

	for (i=0;i<=lup;i++)
	{
		double proba_dw = tmpPrbLoss_down[i];
		if (i>ldown) 
			{proba_dw=0.;}
		Loss2 +=(tmpPrbLoss_up[i]-proba_dw)*(i*GetCtxt()->m_LossUnit-GetCtxt()->m_TrancheDown);
		CheckLoss2.push_back(Loss2);
	}
	Loss2 += (1.-CumtmpPrbLoss_)*(GetCtxt()->m_TrancheUp-GetCtxt()->m_TrancheDown);
	Loss2 /= (GetCtxt()->m_TrancheUp-GetCtxt()->m_TrancheDown);

	//parametres d'integration
	AddressVector	TheParameterVector;

	TheParameterVector.Append(this);
	TheParameterVector.Append(&its_ind_x);
	TheParameterVector.Append(&its_ind_Loss);
	// ------------------------

	ICM_QMatrix<double> COEFS_UP_tu(GetCtxt()->m_nbnames,GetCtxt()->m_IntegrationStep,-999.);
	TheParameterVector.Append(&COEFS_UP_tu);
	ICM_QMatrix<double> COEFS_UP_ts(GetCtxt()->m_nbnames,GetCtxt()->m_IntegrationStep,-999.);
	TheParameterVector.Append(&COEFS_UP_ts);

	ICM_QMatrix<double> COEFS_DW_tu(GetCtxt()->m_nbnames,GetCtxt()->m_IntegrationStep,-999.);
	TheParameterVector.Append(&COEFS_DW_tu);
	ICM_QMatrix<double> COEFS_DW_ts(GetCtxt()->m_nbnames,GetCtxt()->m_IntegrationStep,-999.);
	TheParameterVector.Append(&COEFS_DW_ts);

	int nointegrationcoef=-1;
	TheParameterVector.Append(&nointegrationcoef);

	for (GetCtxt()->m_losses_reset_beg=0;GetCtxt()->m_losses_reset_beg<=lup_stepup;GetCtxt()->m_losses_reset_beg++)
	{
		CumtmpPrbLossUp_stepup=0.;
		CumtmpPrbLossDown_stepup=0.;
		LossDown1 = LossUp1 =0.;

		double CDOloss = MAX((GetCtxt()->m_TrancheUp - GetCtxt()->m_TrancheDown)-
						MAX(GetCtxt()->m_losses_reset_beg*GetCtxt()->m_LossUnit-GetCtxt()->m_TrancheDown,0.),0.)/
						(GetCtxt()->m_TrancheUp - GetCtxt()->m_TrancheDown);

		double stepdw_ = MAX(stepdw-GetCtxt()->m_losses_reset_beg*GetCtxt()->m_LossUnit,0.);

		double stepup_ = stepdw_ 
				+ MAX(
						(GetCtxt()->m_TrancheUp_stepup - GetCtxt()->m_TrancheDown_stepup)-
						MAX(GetCtxt()->m_losses_reset_beg*GetCtxt()->m_LossUnit-GetCtxt()->m_TrancheDown_stepup,0.)
				,0.);

		if ((!stepup_)||(!CDOloss))  
			{continue;}

		double LossProbaUp = tmpPrbLoss_up[GetCtxt()->m_losses_reset_beg];
		double LossProbaDown = tmpPrbLoss_down[GetCtxt()->m_losses_reset_beg];

		tmpPrbLoss_up_stepup.resize(lup_stepup+1,0.);
		tmpPrbLoss_down_stepup.resize(lup_stepup+1,0.);
	
		its_IsUp=true;
		for (i=GetCtxt()->m_losses_reset_beg;i<=lup_stepup;i++)
		{
		GetCtxt()->m_losses_reset_end=i;
		//tmpPrbLoss_up_stepup[i]=HermiteIntegration(ff1::mem_call(&ICM_Gauss1FLossDistrib_LJ::ZeroLoss_TSR_stepup, (*this))).Integrate(GetCtxt()->m_IntegrationStep);
		//tmpPrbLoss_up_stepup[i]*=its_normPi;
		TheIntegrator.Integrate(-6.0, 6.0, ZeroLoss_TSR_stepup_Integrator, &TheParameterVector, value); 
		tmpPrbLoss_up_stepup[i]=value;
		CumtmpPrbLossUp_stepup+=tmpPrbLoss_up_stepup[i];
		}

		its_IsUp=false;
		for (i=GetCtxt()->m_losses_reset_beg;i<=ldown_stepup;i++)
		{
		GetCtxt()->m_losses_reset_end=i;
		//tmpPrbLoss_down_stepup[i]=HermiteIntegration(ff1::mem_call(&ICM_Gauss1FLossDistrib_LJ::ZeroLoss_TSR_stepup, (*this))).Integrate(GetCtxt()->m_IntegrationStep);
		//tmpPrbLoss_down_stepup[i]*=its_normPi;
		TheIntegrator.Integrate(-6.0, 6.0, ZeroLoss_TSR_stepup_Integrator, &TheParameterVector, value); 
		tmpPrbLoss_down_stepup[i]=value;
		CumtmpPrbLossDown_stepup+=tmpPrbLoss_down_stepup[i];
		}

		for (i=GetCtxt()->m_losses_reset_beg;i<=ldown_stepup;i++) 
		{
			LossDown1+=tmpPrbLoss_down_stepup[i]*(i-GetCtxt()->m_losses_reset_beg)*GetCtxt()->m_LossUnit; 
			LossUp1+=tmpPrbLoss_up_stepup[i]*(i-GetCtxt()->m_losses_reset_beg)*GetCtxt()->m_LossUnit;
		}

		for (i=MAX(GetCtxt()->m_losses_reset_beg,ldown_stepup+1);i<=lup_stepup;i++) 
		{LossUp1+=tmpPrbLoss_up_stepup[i]*(i-GetCtxt()->m_losses_reset_beg)*GetCtxt()->m_LossUnit;}

		LossDown1+=(LossProbaDown-CumtmpPrbLossDown_stepup)*stepdw_;
		LossUp1+=(LossProbaUp-CumtmpPrbLossUp_stepup)*stepup_;

		double StepUpLoss = LossUp1-LossDown1;

		Loss1 += CDOloss*StepUpLoss/(stepup_-stepdw_);

		//forward case
		//Loss1 += (LossUp1-LossDown1)/(GetCtxt()->m_TrancheUp_stepup - GetCtxt()->m_TrancheDown_stepup);
		CheckLoss1.push_back(CDOloss*StepUpLoss/(stepup_-stepdw_));

	}

	if (Loss1<0.)
		{Loss1 = 0.;}

	Loss = Loss2 + Loss1 ;

	return (Loss);
}



double ICM_Gauss1FLossDistrib_LJ::ZeroLoss_TSR_stepup(double x, void* params)
{
	its_ind_x+=1;

	unsigned long sizemax=-1;
	(*(AddressVector*)params).GetCount(sizemax);
	ICM_QMatrix<double>* matrix_up_tu = (ICM_QMatrix<double>*)((*(AddressVector*)params)[sizemax-5]);
	ICM_QMatrix<double>* matrix_up_ts = (ICM_QMatrix<double>*)((*(AddressVector*)params)[sizemax-4]);
	ICM_QMatrix<double>* matrix_dw_tu = (ICM_QMatrix<double>*)((*(AddressVector*)params)[sizemax-3]);
	ICM_QMatrix<double>* matrix_dw_ts = (ICM_QMatrix<double>*)((*(AddressVector*)params)[sizemax-2]);
	int* nointegrationcoef = (int*)((*(AddressVector*)params)[sizemax-1]);

/*	if (GetCtxt()->m_IntegrationMethod ==	qGAUSS_HERMITE)
		x	*=	SQRT2;
*/
	double Z_v1_v2_j=0.,Pi_j_u=0.,Pi_j_s=0.,Pu=0.,Ps=0.,result = 0.,value = 0.,cum=0.;
	int j=0;
	int nbnames = GetCtxt()->m_nbnames;
	double tmp_barrier = 0.;

	for (j=1;j<nbnames+1;j++)
	{
		if (its_IsUp)
		{
			if (((*matrix_up_tu)(j-1,*nointegrationcoef) !=-999.) && (*nointegrationcoef>=0))
			{	Pi_j_u=(*matrix_up_tu)(j-1,*nointegrationcoef);	}
			else
			{
				tmp_barrier = GetCtxt()->m_ts_a_up_maturity_t1[j-1] - GetCtxt()->m_ts_b_up_maturity_t1[j-1] * x*its_sqrt2;
				Pi_j_u=NAG_cumul_normal(tmp_barrier);
				if (IsTSR() && GetT2())
					{
					double tmp_barrier_t2=0.;
					double tmp_barrier_end=0.;

					tmp_barrier_t2	=	GetCtxt()->m_ts_a_up_maturity_t2[j-1] - GetCtxt()->m_ts_b_up_maturity_t2[j-1] * x*its_sqrt2;
					Pi_j_u += NAG_cumul_normal(tmp_barrier_t2);
					Pi_j_u -=	NAG_bivariate_normal_dist(tmp_barrier,tmp_barrier_t2,RHO_DIST(GetT1(),GetT2()));
					}
				(*matrix_up_tu)(j-1,*nointegrationcoef)=Pi_j_u;
			}

			if (((*matrix_up_ts)(j-1,*nointegrationcoef) !=-999.) && (*nointegrationcoef>=0))
			{	Pi_j_s=(*matrix_up_ts)(j-1,*nointegrationcoef);	}
			else
			{
				tmp_barrier = GetCtxt()->m_ts_stepup_a_up_maturity_t1[j-1] - GetCtxt()->m_ts_stepup_b_up_maturity_t1[j-1] * x*its_sqrt2;
				Pi_j_s=NAG_cumul_normal(tmp_barrier);
				if (IsTSR() && GetT2())
					{
					double tmp_barrier_t2=0.;
					double tmp_barrier_end=0.;

					tmp_barrier_t2	=	GetCtxt()->m_ts_stepup_a_up_maturity_t2[j-1] - GetCtxt()->m_ts_stepup_b_up_maturity_t2[j-1] * x*its_sqrt2;
					Pi_j_s += NAG_cumul_normal(tmp_barrier_t2);
					Pi_j_s -=	NAG_bivariate_normal_dist(tmp_barrier,tmp_barrier_t2,RHO_DIST(GetT1(),GetT2()));
					}
				(*matrix_up_ts)(j-1,*nointegrationcoef)=Pi_j_s;
			}
		}
		else
		{
			if (((*matrix_dw_tu)(j-1,*nointegrationcoef) !=-999.) && (*nointegrationcoef>=0))
			{	Pi_j_u=(*matrix_dw_tu)(j-1,*nointegrationcoef);	}
			else
			{
				tmp_barrier = GetCtxt()->m_ts_a_dw_maturity_t1[j-1] - GetCtxt()->m_ts_b_dw_maturity_t1[j-1] * x*its_sqrt2;
				Pi_j_u=NAG_cumul_normal(tmp_barrier);
				if (IsTSR() && GetT2())
					{
					double tmp_barrier_t2=0.;
					double tmp_barrier_end=0.;

					tmp_barrier_t2	=	GetCtxt()->m_ts_a_dw_maturity_t2[j-1] - GetCtxt()->m_ts_b_dw_maturity_t2[j-1] * x*its_sqrt2;
					Pi_j_u += NAG_cumul_normal(tmp_barrier_t2);
					Pi_j_u -=	NAG_bivariate_normal_dist(tmp_barrier,tmp_barrier_t2,RHO_DIST(GetT1(),GetT2()));
					}	
				(*matrix_dw_tu)(j-1,*nointegrationcoef)=Pi_j_u;
			}

			if (((*matrix_dw_ts)(j-1,*nointegrationcoef) !=-999.) && (*nointegrationcoef>=0))
			{	Pi_j_s=(*matrix_dw_ts)(j-1,*nointegrationcoef);	}
			else
			{
				tmp_barrier = GetCtxt()->m_ts_stepup_a_dw_maturity_t1[j-1] - GetCtxt()->m_ts_stepup_b_dw_maturity_t1[j-1] * x*its_sqrt2;
				Pi_j_s=NAG_cumul_normal(tmp_barrier);
				if (IsTSR() && GetT2())
					{
					double tmp_barrier_t2=0.;
					double tmp_barrier_end=0.;

					tmp_barrier_t2	=	GetCtxt()->m_ts_stepup_a_dw_maturity_t2[j-1] - GetCtxt()->m_ts_stepup_b_dw_maturity_t2[j-1] * x*its_sqrt2;
					Pi_j_s += NAG_cumul_normal(tmp_barrier_t2);
					Pi_j_s -=	NAG_bivariate_normal_dist(tmp_barrier,tmp_barrier_t2,RHO_DIST(GetT1(),GetT2()));
					}
				(*matrix_dw_ts)(j-1,*nointegrationcoef)=Pi_j_s;
			}
		}

		Pu=Pi_j_u;
		Ps=Pi_j_s;

		for (int v1=0;v1<= MIN(j,GetCtxt()->m_losses_reset_beg);v1++)
			for (int v2=v1;v2<= MIN(j,GetCtxt()->m_losses_reset_end) ;v2++)
		{
			Z_v1_v2_j = 0.;

			int llr = (int) GetCtxt()->m_LossRates[j-1];

			if (v1>=llr) 
			{	Z_v1_v2_j = (1.-Pu)*its_ProbCond_stepup->Elt(v1,v2,j-1) + 
						Ps*its_ProbCond_stepup->Elt(v1-llr,v2-llr,j-1) +
						(Pu - Ps)*(value=its_ProbCond_stepup->Elt(v1,v2-llr,j-1));
			}
			else if (v2<llr)
			{	Z_v1_v2_j = (1.-Pu)*its_ProbCond_stepup->Elt(v1,v2,j-1);
			}	
			else if ((v1<llr) && (llr<=v2))
			{	Z_v1_v2_j = (1.-Pu)*its_ProbCond_stepup->Elt(v1,v2,j-1) + 
						(Pu - Ps)*its_ProbCond_stepup->Elt(v1,v2-llr,j-1);
			}

		its_ProbCond_stepup->SetElt(v1,v2,j,Z_v1_v2_j);
		}
	} 

	result = its_ProbCond_stepup->Elt(GetCtxt()->m_losses_reset_beg,GetCtxt()->m_losses_reset_end,nbnames);

	return (result);

}



