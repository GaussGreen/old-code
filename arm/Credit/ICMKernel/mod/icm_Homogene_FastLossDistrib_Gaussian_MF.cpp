
#include "ARMKernel\glob\firsttoinc.h" 
#include "ICMKernel\mod\icm_Homogene_FastLossDistrib_Gaussian_MF.h"
#include "ICMKernel\glob\icm_maths.h"

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

void
ICM_Gaussian_LossDistrib_MF::Init()
{
	itsIntegrationMethod	=	qGAUSS_HERMITE;		// Gauss-Legendre, Hermite, etc.
	itsIntegrationStep		=	20;

	itsNbStrikes		=	0;

	itsCopulaType		=	qNO_COPULA;
	itsFreedomDegree	=	4;
	itsNbSectors=0;
	its_CorrelationType = TFCT_FULL;

	itsBetasMatrix		=	NULL;

	its_ProbCond = NULL;
	its_ProbCond_Down = NULL;
	its_SectorLoss = NULL;
	its_SectorLoss_Down = NULL;
	
	its_taildistrib_Down	=	0.0;

	itsIndicSectors.clear();
	its_beta_Down.clear();
	its_lambda.clear();
	its_lambda_Down.clear();		 
	its_lossdistrib.clear(); 
	its_lossdistrib_Down.clear(); 
	
	its_coeff_a.clear();
	its_coeff_b.clear();
	its_coeff_c.clear();
	its_coeff_a_down.clear();
	its_coeff_b_down.clear();
	its_coeff_c_down.clear();
	itsStrikes.clear();
	itsIsShortName.clear();
}

//----------------------------------------------------------------------------------------------------------------------
// Constructor
//----------------------------------------------------------------------------------------------------------------------
ICM_Gaussian_LossDistrib_MF::ICM_Gaussian_LossDistrib_MF(const int& nbnames,
										   const ARM_Vector&  pdef,
										   const ARM_Vector&  beta,
										   const ARM_Vector&  LossRates,
										   const int& NbSec,
										   const std::vector<int>& IndSec,
										   const int& CopulaType ,
										   const qIntegratorChoice&	IntegrationMethod ,
										   const int& IntegrationStep)
{
	Init();

	Set(nbnames,pdef,beta,LossRates,NbSec,IndSec,CopulaType,IntegrationMethod,IntegrationStep);

}

//LongShort
/*********************************************************************/
ICM_Gaussian_LossDistrib_MF::ICM_Gaussian_LossDistrib_MF(const int& nbnames,
										   const ARM_Vector&  pdef,
										   const ARM_Vector&  beta,
										   const ARM_Vector&  LossRates,
										   const int& NbSec,
										   const std::vector<int>& IndSec,
										   const std::vector<int>& SortedIndice,
										   const int& CopulaType ,
										   const qIntegratorChoice&	IntegrationMethod ,
										   const int& IntegrationStep)
{
	Init();
	Set(nbnames,pdef,beta,LossRates,NbSec,IndSec,SortedIndice,CopulaType,IntegrationMethod,IntegrationStep);
}

void ICM_Gaussian_LossDistrib_MF::Set(const int& nbnames,
										   const ARM_Vector&  pdef,
										   const ARM_Vector&  beta,
										   const ARM_Vector&  LossRates,
										   const int& NbSec,
										   const std::vector<int>& IndSec,
										   const int& CopulaType ,
										   const qIntegratorChoice&	IntegrationMethod ,
										   const int& IntegrationStep)
{
	int size;
	if (IntegrationStep	==	0) size	=	20;
	else size	=	IntegrationStep;
	
	itsNbSectors = NbSec;

	if (nbnames != GetNbNames())
	{
		ICM_Distribution::Set(nbnames);
		
		its_beta_Down.resize(nbnames);
		its_lambda.resize(nbnames);
		its_lambda_Down.resize(nbnames);

		its_coeff_a.resize(nbnames);
		its_coeff_b.resize(nbnames);
		its_coeff_c.resize(nbnames);
		its_coeff_a_down.resize(nbnames);
		its_coeff_b_down.resize(nbnames);
		its_coeff_c_down.resize(nbnames);
	}

	if ((nbnames != its_nbnames) || (itsIntegrationStep != size) || 
		(!its_ProbCond)||(!its_ProbCond_Down)||(!its_SectorLoss)||(!its_SectorLoss_Down))
	{
		if (its_ProbCond) delete its_ProbCond;
		its_ProbCond= new ICM_QCubix<double>(nbnames+1,size*size,1,0.);

		if (its_ProbCond_Down) delete its_ProbCond_Down;
		its_ProbCond_Down= new ICM_QCubix<double>(nbnames+1,size*size,1,0.);

		if (its_SectorLoss) delete its_SectorLoss;
		its_SectorLoss= new ICM_QCubix<double>(NbSec+1,size,1,0.);

		if (its_SectorLoss_Down) delete its_SectorLoss_Down;
		its_SectorLoss_Down= new ICM_QCubix<double>(NbSec+1,size,1,0.);		
	}

	SetNbNames(nbnames);

	//MaJ du vecteur d'indic du secteur
	if (IndSec.size() == nbnames)
		itsIndicSectors=IndSec;
	else
		itsIndicSectors.resize(0);

	SetUniqueBeta(beta); 

	itsCopulaType = CopulaType;
	itsIntegrationStep = IntegrationStep;
	itsIntegrationMethod = IntegrationMethod;
	SetIntegratorType(itsIntegrationMethod, itsIntegrationStep);
	// Si le nombre de pas et la méthode utilisés ne sont pas compatibles, la parité du nombre de pas est prioritaire
	if ( (itsIntegrationStep % 2 == 0) && (itsIntegrationMethod == qGAUSS_LEGENDRE))
		itsIntegrationMethod = qGAUSS_HERMITE;
	else if ( (itsIntegrationStep % 2 == 1) && (itsIntegrationMethod == qGAUSS_HERMITE))
		itsIntegrationMethod = qGAUSS_LEGENDRE;
	
	SetPdefAtMaturity(pdef); 
	compute_min_pdef(GetPdefAtMaturity());

	SetIntLossRates(LossRates);
	check_homogeneous(LossRates);
	compute_barrier();

	its_lossdistrib_Down.resize(1);
}



//LongShort
/*********************************************************************/
void ICM_Gaussian_LossDistrib_MF::Set(const int& nbnames,
										   const ARM_Vector&  pdef,
										   const ARM_Vector&  beta,
										   const ARM_Vector&  LossRates,
										   const int& NbSec,
										   const std::vector<int>& IndSec,
										   const std::vector<int>& SortedIndice,
										   const int& CopulaType ,
										   const qIntegratorChoice&	IntegrationMethod ,
										   const int& IntegrationStep)
{
	Set(nbnames,pdef,beta,LossRates,NbSec,IndSec,CopulaType,IntegrationMethod,IntegrationStep);
	ComputeNbNameShort(LossRates);
	SetSortedIndices(SortedIndice);
	SetIsShortName(LossRates); 
}

void ICM_Gaussian_LossDistrib_MF::BitwiseCopy(const ARM_Object* src)
{

    ICM_Gaussian_LossDistrib_MF* srcdistrib = (ICM_Gaussian_LossDistrib_MF *) src;

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

	if (srcdistrib->its_SectorLoss)
		its_SectorLoss = (ICM_QCubix<double>*) srcdistrib->its_SectorLoss->Clone();

	if (srcdistrib->its_SectorLoss_Down)
		its_SectorLoss_Down = (ICM_QCubix<double>*) srcdistrib->its_SectorLoss_Down->Clone();

	itsNbSectors=srcdistrib->itsNbSectors;
	itsIndicSectors=srcdistrib->itsIndicSectors;
	its_beta_Down=srcdistrib->its_beta_Down;
	its_lambda=srcdistrib->its_lambda;
	its_lambda_Down=srcdistrib->its_lambda_Down;
	its_lossdistrib=srcdistrib->its_lossdistrib;
	its_lossdistrib_Down=srcdistrib->its_lossdistrib_Down;
	
	
	its_coeff_a=srcdistrib->its_coeff_a;
	its_coeff_b=srcdistrib->its_coeff_b;
	its_coeff_c=srcdistrib->its_coeff_c;
	its_coeff_a_down=srcdistrib->its_coeff_a_down;
	its_coeff_b_down=srcdistrib->its_coeff_b_down;
	its_coeff_c_down=srcdistrib->its_coeff_c_down;
	
	itsStrikes=srcdistrib->itsStrikes;
	itsIsShortName=srcdistrib->itsIsShortName;
}

// -------------
//	Copy Method 
// -------------
void ICM_Gaussian_LossDistrib_MF::Copy(const ARM_Object* src)
{
	ICM_Distribution::Copy(src);
    BitwiseCopy(src);
}


ARM_Object* ICM_Gaussian_LossDistrib_MF::Clone(void)
{
     ICM_Gaussian_LossDistrib_MF* theClone = new ICM_Gaussian_LossDistrib_MF();

     theClone->Copy(this);
 
     return(theClone);
}


//----------------------------------------------------------------------------------------------------------------------

void ICM_Gaussian_LossDistrib_MF::SetSmileParameters(const int& NbStrikes,
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

double ICM_Gaussian_LossDistrib_MF::compute_expectedlosstranche(const double& tranche_up, 
									   const double& tranche_down, 
									   const double& lossunit,
									   //ICM_QMatrix<double>* ShiftMatrix,
									   const int& Tenor,
									   const int& Issuer)
{
	double exp_loss_tranche_down = 0.;
	double exp_loss_tranche_up = 0.;
	
	int lup=floor(tranche_up/lossunit);
	int ldown=floor(tranche_down/lossunit);

	int i=0;	
	
	its_ProbCond->ResizeWithCopy(lup+1,0.);
	its_ProbCond_Down->ResizeWithCopy(ldown+1,0.);

	//Loss sectorielle
	its_SectorLoss->ResizeWithCopy(lup+1,0.);
	its_SectorLoss_Down->ResizeWithCopy(ldown+1,0.);

	its_lossdistrib.resize(lup+1);
	its_lossdistrib_Down.resize(ldown+1);
	
	// Compute the Distribution (by default, ldown = 0)
	if (tranche_down)
		compute_distrib(lup, ldown);
	else
		compute_distrib(lup);
	
	//View the loss ditribution
	/*#ifdef _DEBUG

	FILE *stream = fopen("c:\\temp\\LossDistrib.txt", "w+");
	this->View("",stream);
	fclose(stream);
	#endif */

	int l	=	0;
	double	loss_level;

	// first loop, deal with beta_down and up
	for (l=1; l<=ldown; l++) 
	{
		loss_level	=	l*lossunit;
		exp_loss_tranche_down	+=	its_lossdistrib_Down[l] * loss_level;		
		exp_loss_tranche_up		+=	its_lossdistrib[l] * loss_level;
	}
		
	// second loop, deal only with beat_up
	for (l=ldown+1;l<=lup;l++)
	{
		loss_level	=	l*lossunit;
		exp_loss_tranche_up		+=	its_lossdistrib[l] * loss_level;
	}

	// add the tail
	if (tranche_down)
		exp_loss_tranche_down	+=	tranche_down * its_taildistrib_Down;
	exp_loss_tranche_up		+=	tranche_up * its_taildistrib;

	return (exp_loss_tranche_up - exp_loss_tranche_down);
}

double ICM_Gaussian_LossDistrib_MF::compute_expectedlosstranche_LongShort(const double& tranche_up, 
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

	//Loss sectorielle
	its_SectorLoss->ResizeWithInitialize(lmin+lup+1,0.);
	its_SectorLoss_Down->ResizeWithInitialize(lmin+ldown+1,0.);

 	its_lossdistrib.resize(lmin+lup+1,0.); 
	its_lossdistrib_Down.resize(lmin+ldown+1,0.);
	
	// Compute the Distribution (by default, ldown = 0)
	if (tranche_down)
		compute_distrib_LongShort(lup, ldown, lmin);
	else
		compute_distrib_LongShort(lup, lmin);
	
	//View the loss ditribution
	 #ifdef _DEBUG
	FILE *stream = fopen("c:\\temp\\LossDistribLongShort.txt", "w+");
	this->ViewS("",stream,39);
	fclose(stream);
	#endif 
 

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


void ICM_Gaussian_LossDistrib_MF::InitProbCondPortSize0(const int& lmin, const int& ind_x)
{
	its_ProbCond->Elt(0,ind_x,lmin) = 1;
	its_ProbCond_Down->Elt(0,ind_x,lmin) = 1;
}

double ICM_Gaussian_LossDistrib_MF::Compute_pk(const int& ind_x, const int& ind_sortedname)
{
	double tmp_barrier=0.;
	double pk=0.;
	double x1,x2;
	int nb_factor_cond = TheIntegrator.GetIntegrationStep();
	
	if (itsIntegrationMethod	==	qGAUSS_LEGENDRE)
	{
		//Les 2 facteurs sont stockés linéairement x1 -> V et x2 -> Vi
		x1=TheIntegrator.GetAbscissa(ind_x/nb_factor_cond);
		x2=TheIntegrator.GetAbscissa(ind_x%nb_factor_cond);
	}
	else
	{
		x1=TheIntegrator.GetAbscissa(ind_x/nb_factor_cond)*SQRT2;
		x2=TheIntegrator.GetAbscissa(ind_x%nb_factor_cond)*SQRT2;
	}

	//Barrier
	tmp_barrier = its_coeff_a[ind_sortedname] - its_coeff_b[ind_sortedname] * x1 - its_coeff_c[ind_sortedname] * x2;
	
	// conditional default probability
	pk=NAG_cumul_normal(tmp_barrier);

	return pk;

}

double ICM_Gaussian_LossDistrib_MF::Compute_pk_down(const int& ind_x, const int& ind_sortedname)
{
	double tmp_barrier_down=0.;
	double pk_down=0.;
	double x1,x2;
	int nb_factor_cond = TheIntegrator.GetIntegrationStep();
	
	if (itsIntegrationMethod	==	qGAUSS_LEGENDRE)
	{
		//Les 2 facteurs sont stockés linéairement x1 -> V et x2 -> Vi
		x1=TheIntegrator.GetAbscissa(ind_x/nb_factor_cond);
		x2=TheIntegrator.GetAbscissa(ind_x%nb_factor_cond);
	}
	else
	{
		x1=TheIntegrator.GetAbscissa(ind_x/nb_factor_cond)*SQRT2;
		x2=TheIntegrator.GetAbscissa(ind_x%nb_factor_cond)*SQRT2;
	}	
	//Barrier
	tmp_barrier_down = its_coeff_a_down[ind_sortedname] - its_coeff_b_down[ind_sortedname] * x1 - its_coeff_c_down[ind_sortedname] * x2;		;		
	
	// conditional default probability
	pk_down=NAG_cumul_normal(tmp_barrier_down);

	return pk_down;

}

void ICM_Gaussian_LossDistrib_MF::ComputeProbCond(const int& ind_x, const int& ind_loss, const int& ind_name, const int& ind_sortedname, const double& pk)
{
	double tmp1=0.;
	double tmp2=0.;
	//double tmp_barrier=0.;
	//double pk=0.;
	//double x1,x2;
	//int nb_factor_cond = TheIntegrator.GetIntegrationStep();
	
	/* (itsIntegrationMethod	==	qGAUSS_LEGENDRE)
	{
		//Les 2 facteurs sont stockés linéairement x1 -> V et x2 -> Vi
		x1=TheIntegrator.GetAbscissa(ind_x/nb_factor_cond);
		x2=TheIntegrator.GetAbscissa(ind_x%nb_factor_cond);
	}
	else
	{
		x1=TheIntegrator.GetAbscissa(ind_x/nb_factor_cond)*SQRT2;
		x2=TheIntegrator.GetAbscissa(ind_x%nb_factor_cond)*SQRT2;
	}*/
							
	int ind_lastloss = ind_loss - GetIntLossRates()[ind_sortedname];

	//Barrier Multi factor
	//tmp_barrier = its_coeff_a[ind_sortedname] - its_coeff_b[ind_sortedname] * x1 - its_coeff_c[ind_sortedname] * x2;
	
	// conditional default probability
	tmp1=its_ProbCond->Elt(ind_name-1,ind_x,ind_loss);
	//pk=NAG_cumul_normal(tmp_barrier);
	tmp2		=	tmp1 * (1.-pk);

	if (ind_lastloss >= 0)
		tmp2	+=	pk * its_ProbCond->Elt(ind_name-1,ind_x,ind_lastloss);

	its_ProbCond->Elt(ind_name,ind_x,ind_loss)=tmp2;
}

void ICM_Gaussian_LossDistrib_MF::ComputeProbCondDown(const int& ind_x, const int& ind_loss, const int& ind_name, const int& ind_sortedname, const double& pk_down)
{
	double tmp1_down=0.;
	double tmp2_down=0;
	//double tmp_barrier_down=0.;
	//double pk_down=0.;
	//double x1,x2;
	//int nb_factor_cond = TheIntegrator.GetIntegrationStep();
	
	int ind_lastloss = ind_loss - GetIntLossRates()[ind_sortedname];
	
	/*if (itsIntegrationMethod	==	qGAUSS_LEGENDRE)
	{
		//Les 2 facteurs sont stockés linéairement x1 -> V et x2 -> Vi
		x1=TheIntegrator.GetAbscissa(ind_x/nb_factor_cond);
		x2=TheIntegrator.GetAbscissa(ind_x%nb_factor_cond);
	}
	else
	{
		x1=TheIntegrator.GetAbscissa(ind_x/nb_factor_cond)*SQRT2;
		x2=TheIntegrator.GetAbscissa(ind_x%nb_factor_cond)*SQRT2;
	}*/
	
	//Barrier
	//tmp_barrier_down = its_coeff_a_down[ind_sortedname] - its_coeff_b_down[ind_sortedname] * x1 - its_coeff_c_down[ind_sortedname] * x2;		
	
	// conditional default probability
	tmp1_down=its_ProbCond_Down->Elt(ind_name-1,ind_x,ind_loss);
	//pk_down=NAG_cumul_normal(tmp_barrier_down);
	tmp2_down	=	tmp1_down * (1.-pk_down);

	if (ind_lastloss >= 0)
		tmp2_down	+=	pk_down * its_ProbCond_Down->Elt(ind_name-1,ind_x,ind_lastloss);

	its_ProbCond_Down->Elt(ind_name,ind_x,ind_loss)=tmp2_down;
}

void ICM_Gaussian_LossDistrib_MF::ComputeProbCondShort(const int& ind_x, const int& ind_loss, const int& ind_name, const int& lmin, const int& ind_sortedname, const double& pk)
{
	double tmp1=0.;
	double tmp2=0.;
	//double tmp_barrier=0.;
	//double pk=0.;
	//double x1,x2;
	//int nb_factor_cond = TheIntegrator.GetIntegrationStep();

	/*if (itsIntegrationMethod	==	qGAUSS_LEGENDRE)
	{
		//Les 2 facteurs sont stockés linéairement x1 -> V et x2 -> Vi
		x1=TheIntegrator.GetAbscissa(ind_x/nb_factor_cond);
		x2=TheIntegrator.GetAbscissa(ind_x%nb_factor_cond);
	}
	else
	{
		x1=TheIntegrator.GetAbscissa(ind_x/nb_factor_cond)*SQRT2;
		x2=TheIntegrator.GetAbscissa(ind_x%nb_factor_cond)*SQRT2;
	}*/
	
	int ind_lastloss = ind_loss - GetIntLossRates()[ind_sortedname];

	//Barrier
	//tmp_barrier = its_coeff_a[ind_sortedname] - its_coeff_b[ind_sortedname] * x1 - its_coeff_c[ind_sortedname] * x2;
	
	// conditional default probability
	tmp1=its_ProbCond->Elt(ind_name-1,ind_x,ind_loss);
	//pk=NAG_cumul_normal(tmp_barrier);
	tmp2		=	tmp1 * (1.-pk);

	if (ind_lastloss <= lmin)
		tmp2	+=	pk * its_ProbCond->Elt(ind_name-1,ind_x,ind_lastloss);

	its_ProbCond->Elt(ind_name,ind_x,ind_loss)=tmp2;
}

void ICM_Gaussian_LossDistrib_MF::ComputeProbCondDownShort(const int& ind_x, const int& ind_loss, const int& ind_name, const int& lmin, const int& ind_sortedname, const double& pk_down)
{
	double tmp1_down=0.;
	double tmp2_down=0;
	//double tmp_barrier_down=0.;
	//double pk_down=0.;
	//double x1,x2;
	//int nb_factor_cond = TheIntegrator.GetIntegrationStep();
	
	/*if (itsIntegrationMethod	==	qGAUSS_LEGENDRE)
	{
		//Les 2 facteurs sont stockés linéairement x1 -> V et x2 -> Vi
		x1=TheIntegrator.GetAbscissa(ind_x/nb_factor_cond);
		x2=TheIntegrator.GetAbscissa(ind_x%nb_factor_cond);
	}
	else
	{
		x1=TheIntegrator.GetAbscissa(ind_x/nb_factor_cond)*SQRT2;
		x2=TheIntegrator.GetAbscissa(ind_x%nb_factor_cond)*SQRT2;
	}*/

	int ind_lastloss = ind_loss - GetIntLossRates()[ind_sortedname];
	//Barrier
	//tmp_barrier_down = its_coeff_a_down[ind_sortedname] - its_coeff_b_down[ind_sortedname] * x1 - its_coeff_c_down[ind_sortedname] * x2;		;		
	
	// conditional default probability
	tmp1_down=its_ProbCond_Down->Elt(ind_name-1,ind_x,ind_loss);
	//pk_down=NAG_cumul_normal(tmp_barrier_down);
	tmp2_down	=	tmp1_down * (1.-pk_down);

	if (ind_lastloss <= lmin)
		tmp2_down	+=	pk_down * its_ProbCond_Down->Elt(ind_name-1,ind_x,ind_lastloss);

	its_ProbCond_Down->Elt(ind_name,ind_x,ind_loss)=tmp2_down;
}
	
// Estimation de la distribution de perte conditionnelle
void ICM_Gaussian_LossDistrib_MF::compute_distrib(const int& lup)
{
	double tmpPrbLoss=0.;
	
	int l = 0;
	
	int ind_x = 0;
	int ind_x2 = 0;
	int ind_name = 0;
	int ind_sortedname = 0;
	int ind_GlobalName =0;
	int its_indSect = 0;
	int nb_name_secteur = 0;
	int ind_loss = 0;
	int ind_loss2 = 0;
	double res = 0.;
	double cumul_distrib = 0.;

	int nb_fact_cond = TheIntegrator.GetIntegrationStep();


	switch (itsCopulaType)
	{
		case qNO_COPULA:	
		case qGAUSSIAN:

			//Premiere etape : estimation des probas conditionnelles (V, Vi) pour chaque secteur
			
			
			//Boucle sur les différents secteurs
			for (its_indSect=0; its_indSect<itsNbSectors; its_indSect++)
			{
				nb_name_secteur = GetNbNameSector(its_indSect);
				double pk=0.;
				//Boucle sur les facteurs sous-jacent
				//On a N*N points d'intégration
				//On peut faire tous les secteurs en même temps
				for (ind_x=0; ind_x<nb_fact_cond*nb_fact_cond ; ind_x++)
				{
					//Init : Port size = 0
					InitProbCondPortSize0(0, ind_x);
					
					//Boucle sur les émetteurs : on suppose que le veceur SortedIndice est trié selon 
					//le long short et le secteur correspondant
					for (ind_name=1; ind_name<=nb_name_secteur; ind_name++)
					{
						ind_sortedname = this->its_sorted_indices[ind_GlobalName+ind_name-1];
						pk=0.;
						pk = Compute_pk(ind_x,ind_sortedname);

						//Boucle sur les loss
						for (ind_loss=0; ind_loss<=lup; ind_loss++)
						{
							ComputeProbCond(ind_x, ind_loss, ind_name, ind_sortedname, pk);
						}
					}
				}
				//ind_GlobalName = ind_name-2;
				ind_GlobalName = ind_GlobalName + ind_name - 1;
				
				//Sector Loss distribution computation : integration par rapport à chaque Vi
				for (ind_x=0; ind_x<nb_fact_cond ; ind_x++)
				{
					//Boucle sur les loss
					integrate_SectorLoss(lup,nb_fact_cond,nb_name_secteur,ind_x,its_indSect);

				}

			}

			//Initialisation de its_SectorLoss(itsNbSectors)
			initialize_SectorLoss(lup, nb_fact_cond);

			//Boucle sur les secteurs de manière à estimer la distribution conditionnelle par rapport a V
			for (its_indSect=1; its_indSect<itsNbSectors; its_indSect++)
			{	
				for (ind_x=0; ind_x<nb_fact_cond ; ind_x++)
				{	vector<double> SumElt(lup+1,0.);
					for (ind_loss=0; ind_loss<=lup; ind_loss++)
					{				
						for(ind_loss2=0; ind_loss2<=ind_loss; ind_loss2++)
							SumElt[ind_loss] += its_SectorLoss->Elt(itsNbSectors,ind_x,ind_loss-ind_loss2)*its_SectorLoss->Elt(its_indSect,ind_x,ind_loss2);
					
					}
					for (ind_loss=0; ind_loss<=lup; ind_loss++)
						its_SectorLoss->SetElt(itsNbSectors,ind_x,ind_loss,SumElt[ind_loss]);
					
				}
			}
					
			//Integration de la distribution totale 
			//Boucle sur chaque loss 
			integrate_LossDistrib(lup,nb_fact_cond);

	}

	for (l=0;l<=lup;l++)
	{
		cumul_distrib += its_lossdistrib[l];
	}
	its_taildistrib			=	1.0 - cumul_distrib;

}
	
void ICM_Gaussian_LossDistrib_MF::compute_distrib(const int& lup, const int& ldown)
{

	double tmpPrbLoss=0.;
	
	int l = 0;
	
	int ind_x = 0;
	int ind_x2 = 0;
	int ind_name = 0;
	int ind_sortedname = 0;
	int ind_GlobalName =0;
	int its_indSect = 0;
	int nb_name_secteur = 0;
	int ind_loss = 0;
	int ind_loss2 = 0;
	double res = 0.;
	double resdown = 0.;
	double cumul_distrib = 0.;
	double cumul_distrib_down = 0.;

	int nb_fact_cond = TheIntegrator.GetIntegrationStep();


	switch (itsCopulaType)
	{
		case qNO_COPULA:	
		case qGAUSSIAN:

			//Premiere etape : estimation des probas conditionnelles (V, Vi) pour chaque secteur
			
			
			//Boucle sur les différents secteurs
			for (its_indSect=0; its_indSect<itsNbSectors; its_indSect++)
			{
				nb_name_secteur = GetNbNameSector(its_indSect);
				double pk=0.;
				double pk_down=0.;
				//Boucle sur les facteurs sous-jacent
				//On a N*N points d'intégration
				//On peut faire tous les secteurs en même temps
				for (ind_x=0; ind_x<nb_fact_cond*nb_fact_cond ; ind_x++)
				{
					//Init : Port size = 0
					InitProbCondPortSize0(0, ind_x);
					
					//Boucle sur les émetteurs : on suppose que le veceur SortedIndice est trié selon 
					//le long short et le secteur correspondant
					for (ind_name=1; ind_name<=nb_name_secteur; ind_name++)
					{
						ind_sortedname = this->its_sorted_indices[ind_GlobalName+ind_name-1];
						pk =0.;
						pk_down=0.;
						pk = Compute_pk(ind_x,ind_sortedname);
						pk_down = Compute_pk_down(ind_x,ind_sortedname);

						//Boucle sur les loss
						for (ind_loss=0; ind_loss<=ldown; ind_loss++)
						{
							ComputeProbCond(ind_x, ind_loss, ind_name, ind_sortedname, pk);
							ComputeProbCondDown(ind_x, ind_loss, ind_name, ind_sortedname, pk_down);
						}
						for (ind_loss=ldown+1; ind_loss<=lup; ind_loss++)
						{
							ComputeProbCond(ind_x, ind_loss, ind_name, ind_sortedname, pk);
						}
					}
				}
				ind_GlobalName = ind_GlobalName + ind_name - 1;
				
				//Sector Loss distribution computation : integration par rapport à chaque Vi
				for (ind_x=0; ind_x<nb_fact_cond ; ind_x++)
				{
					//Boucle sur les loss
					integrate_SectorLoss(lup,nb_fact_cond,nb_name_secteur,ind_x,its_indSect);
					integrate_SectorLossDown(ldown,nb_fact_cond,nb_name_secteur,ind_x,its_indSect);
				}
			}

			//Initialisation de its_SectorLoss(itsNbSectors)
			initialize_SectorLossDown(ldown, nb_fact_cond);
			initialize_SectorLoss(lup, nb_fact_cond);


			//Boucle sur les secteurs de manière à estimer la distribution conditionnelle par rapport a V
			for (its_indSect=1; its_indSect<itsNbSectors; its_indSect++)
			{	
				for (ind_x=0; ind_x<nb_fact_cond ; ind_x++)
				{	vector<double> SumElt(lup+1,0.);
					for (ind_loss=0; ind_loss<=lup; ind_loss++)
					{				
						for(ind_loss2=0; ind_loss2<=ind_loss; ind_loss2++)
							SumElt[ind_loss] += its_SectorLoss->Elt(itsNbSectors,ind_x,ind_loss-ind_loss2)*its_SectorLoss->Elt(its_indSect,ind_x,ind_loss2);
					
					}
					for (ind_loss=0; ind_loss<=lup; ind_loss++)
						its_SectorLoss->SetElt(itsNbSectors,ind_x,ind_loss,SumElt[ind_loss]);
					
				}
				for (ind_x=0; ind_x<nb_fact_cond ; ind_x++)
				{	vector<double> SumElt(ldown+1,0.);
					for (ind_loss=0; ind_loss<=ldown; ind_loss++)
					{				
						for(ind_loss2=0; ind_loss2<=ind_loss; ind_loss2++)
							SumElt[ind_loss] += its_SectorLoss_Down->Elt(itsNbSectors,ind_x,ind_loss-ind_loss2)*its_SectorLoss_Down->Elt(its_indSect,ind_x,ind_loss2);
					
					}
					for (ind_loss=0; ind_loss<=ldown; ind_loss++)
						its_SectorLoss_Down->SetElt(itsNbSectors,ind_x,ind_loss,SumElt[ind_loss]);
					
				}			
			}
					
			//Integration de la distribution totale 
			//Boucle sur chaque loss 
			integrate_LossDistrib(lup,nb_fact_cond);
			integrate_LossDistribDown(ldown,nb_fact_cond);

	}

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

void ICM_Gaussian_LossDistrib_MF::compute_distrib_LongShort(const int& lup, const int& lmin)
{
	double tmpPrbLoss=0.;
	int min_lup =0,max_lup =0;
	int maxloss; 
	max_lup		=	lup+lmin;

	int l = 0;
	
	int ind_x = 0;
	int ind_x2 = 0;
	int ind_name = 0;
	int ind_sortedname = 0;
	int ind_GlobalName =0;
	int its_indSect = 0;
	int nb_name_secteur = 0;
	int ind_loss = 0;
	int ind_loss2 = 0;
	double res = 0.;
	double cumul_distrib = 0.;
	bool IsShortName=false; 
	std::deque<bool> /* FIXMEFRED: mig.vc8 (23/05/2007 14:27:19): dont use std::vector<bool>. see item 18 in "Effective STL" from Scott Meyers*/ IsShortSector(itsNbSectors);

	int nb_fact_cond = TheIntegrator.GetIntegrationStep();


	switch (itsCopulaType)
	{
		case qNO_COPULA:	
		case qGAUSSIAN:

			//Premiere etape : estimation des probas conditionnelles (V, Vi) pour chaque secteur
			
			
			//Boucle sur les différents secteurs
			for (its_indSect=0; its_indSect<itsNbSectors; its_indSect++)
			{
				its_ProbCond->ResizeWithInitialize(max_lup+1,0.);
				nb_name_secteur = GetNbNameSector(its_indSect);
				//Boucle sur les facteurs sous-jacent
				//On a N*N points d'intégration
				//On peut faire tous les secteurs en même temps
				for (ind_x=0; ind_x<nb_fact_cond*nb_fact_cond ; ind_x++)
				{
					//Init : Port size = 0
					double pk=0.;
					IsShortName = itsIsShortName[this->its_sorted_indices[ind_GlobalName]];
					if (IsShortName) 
						InitProbCondPortSize0(lmin, ind_x);
					else 
						InitProbCondPortSize0(0, ind_x);
					
					//Boucle sur les émetteurs : on suppose que le veceur SortedIndice est trié selon 
					//le long short et le secteur correspondant
					for (ind_name=1; ind_name<=nb_name_secteur; ind_name++)
					{
						ind_sortedname = this->its_sorted_indices[ind_GlobalName+ind_name-1];
						pk=0.;
						pk = Compute_pk(ind_x,ind_sortedname);

						//Boucle sur les loss
						if (IsShortName) 
						{	for (ind_loss=0; ind_loss<=lmin; ind_loss++)
							{
								ComputeProbCondShort(ind_x, ind_loss, ind_name, lmin, ind_sortedname, pk);
							}
						} 
						else 
						{	for (ind_loss=0; ind_loss<=max_lup; ind_loss++)
							{
								ComputeProbCond(ind_x, ind_loss, ind_name, ind_sortedname, pk);
							}
						} 

					}
					/*FILE *stream = fopen("c:\\temp\\LossDistrib.txt", "w+");
					this->ViewF("",stream,ind_x);
					fclose(stream);*/

					IsShortSector[its_indSect]=IsShortName;
				}
				ind_GlobalName = ind_GlobalName + ind_name - 1;
				
				//Sector Loss distribution computation : integration par rapport à chaque Vi
				for (ind_x=0; ind_x<nb_fact_cond ; ind_x++)
				{
				
					if (IsShortName) 
						maxloss = lmin;
					else
						maxloss = max_lup;

					//Boucle s0.ur les loss
					integrate_SectorLoss(maxloss,nb_fact_cond,nb_name_secteur,ind_x,its_indSect);
					/*for (ind_loss=0; ind_loss<=maxloss; ind_loss++)
					{	
						res=0.;
						//Intgération sur le facteur sectoriel Vi
						for (ind_x2=0; ind_x2< ; ind_x2++)
							res += its_ProbCond->Elt(nb_name_secteur,ind_x*nb_fact_cond+ind_x2,ind_loss)*TheIntegrator.GetWeight(ind_x2);
						
						if (itsIntegrationMethod	==	qGAUSS_HERMITE)
							its_SectorLoss->SetElt(its_indSect,ind_x,ind_loss,res*ONEOVERSQRTPI);
						else if (itsIntegrationMethod	==	qGAUSS_LEGENDRE)
							its_SectorLoss->SetElt(its_indSect,ind_x,ind_loss,res/SQRT2PI);
					}	*/
					/*FILE *stream = fopen("c:\\temp\\SecLossDistrib.txt", "w+");
					this->ViewS("",stream,ind_x);
					fclose(stream);*/
				}
			}

			//Initialisation de its_SectorLoss(itsNbSectors)
			initialize_SectorLoss(lmin,nb_fact_cond);
			initialize_SectorLoss0(lmin,max_lup, nb_fact_cond);

			//Boucle sur les secteurs de manière à estimer la distribution conditionnelle par rapport a V
			for (its_indSect=1; its_indSect<itsNbSectors; its_indSect++)
			{	if (IsShortSector[its_indSect]) 
				{	for (ind_x=0; ind_x<nb_fact_cond ; ind_x++)
						{	
							vector<double> SumElt(lmin+1,0.);
							for (ind_loss=lmin; ind_loss>=0; ind_loss--)
							{				
								for(ind_loss2=lmin; ind_loss2>=ind_loss; ind_loss2--)
									SumElt[ind_loss] += its_SectorLoss->Elt(itsNbSectors,ind_x,ind_loss-ind_loss2+lmin)*its_SectorLoss->Elt(its_indSect,ind_x,ind_loss2);
							
							}
						for (ind_loss=0; ind_loss<=lmin; ind_loss++)
							its_SectorLoss->SetElt(itsNbSectors,ind_x,ind_loss,SumElt[ind_loss]);

						for (ind_loss=lmin+1; ind_loss<=max_lup; ind_loss++)
							its_SectorLoss->SetElt(itsNbSectors,ind_x,ind_loss,0.);

						}
				}
				else
				{	for (ind_x=0; ind_x<nb_fact_cond ; ind_x++)
						{	
							vector<double> SumElt(max_lup+1,0.);
							for (ind_loss=0; ind_loss<=max_lup; ind_loss++)
							{				
								for(ind_loss2=0; ind_loss2<=ind_loss; ind_loss2++)
									SumElt[ind_loss] += its_SectorLoss->Elt(itsNbSectors,ind_x,ind_loss-ind_loss2)*its_SectorLoss->Elt(its_indSect,ind_x,ind_loss2);
							
							}
						for (ind_loss=0; ind_loss<=max_lup; ind_loss++)
							its_SectorLoss->SetElt(itsNbSectors,ind_x,ind_loss,SumElt[ind_loss]);
					
						/*FILE *stream = fopen("c:\\temp\\SecLossDistrib.txt", "w+");
						this->ViewS("",stream,ind_x);
						fclose(stream);*/
						}
				}
			}
					
			//Integration de la distribution totale 
			//Boucle sur chaque loss 
			integrate_LossDistrib(max_lup,nb_fact_cond);

	}

	for (l=0;l<=max_lup;l++)
	{
		cumul_distrib += its_lossdistrib[l];
	}
	its_taildistrib			=	1.0 - cumul_distrib;


}

void ICM_Gaussian_LossDistrib_MF::compute_distrib_LongShort(const int& lup, const int& ldown, const int& lmin)
{
	double tmpPrbLoss=0.;

	int min_lup =0,max_lup =0,max_ldown=0;
	int maxloss; 
	max_lup		=	lup+lmin;
	max_ldown	=	ldown+lmin;

	int l = 0;
	
	int ind_x = 0;
	int ind_x2 = 0;
	int ind_name = 0;
	int ind_sortedname = 0;
	int ind_GlobalName =0;
	int its_indSect = 0;
	int nb_name_secteur = 0;
	int ind_loss = 0;
	int ind_loss2 = 0;
	double res = 0.;
	double resdown = 0.;
	double cumul_distrib = 0.;
	double cumul_distrib_down = 0.;
	bool IsShortName=false; 
	std::deque<bool> /* FIXMEFRED: mig.vc8 (23/05/2007 14:27:19): dont use std::vector<bool>. see item 18 in "Effective STL" from Scott Meyers*/ IsShortSector(itsNbSectors);
	

	int nb_fact_cond = TheIntegrator.GetIntegrationStep();


	switch (itsCopulaType)
	{
		case qNO_COPULA:	
		case qGAUSSIAN:

			//Premiere etape : estimation des probas conditionnelles (V, Vi) pour chaque secteur
			
			
			//Boucle sur les différents secteurs
			for (its_indSect=0; its_indSect<itsNbSectors; its_indSect++)
			{
				its_ProbCond->ResizeWithInitialize(max_lup+1,0.);
				its_ProbCond_Down->ResizeWithInitialize(max_ldown+1,0.);
				nb_name_secteur = GetNbNameSector(its_indSect);
				double pk=0.;
				double pk_down=0.;
				//Boucle sur les facteurs sous-jacent
				//On a N*N points d'intégration
				//On peut faire tous les secteurs en même temps
				for (ind_x=0; ind_x<nb_fact_cond*nb_fact_cond ; ind_x++)
				{
					//Init : Port size = 0
					IsShortName = itsIsShortName[this->its_sorted_indices[ind_GlobalName]];
					if (IsShortName) 
						InitProbCondPortSize0(lmin, ind_x);
					else 
						InitProbCondPortSize0(0, ind_x);
					
					//Boucle sur les émetteurs : on suppose que le veceur SortedIndice est trié selon 
					//le long short et le secteur correspondant
					for (ind_name=1; ind_name<=nb_name_secteur; ind_name++)
					{
						ind_sortedname = this->its_sorted_indices[ind_GlobalName+ind_name-1];
						pk =0.;
						pk_down=0.;
						pk = Compute_pk(ind_x,ind_sortedname);
						pk_down = Compute_pk_down(ind_x,ind_sortedname);

						//Boucle sur les loss
						if (IsShortName) 
						{	for (ind_loss=0; ind_loss<=lmin; ind_loss++)
							{	
								ComputeProbCondShort(ind_x, ind_loss, ind_name, lmin, ind_sortedname, pk);
								ComputeProbCondDownShort(ind_x, ind_loss, ind_name, lmin, ind_sortedname, pk_down);
							}
						} 
						else 
						{	for (ind_loss=0; ind_loss<=max_ldown; ind_loss++)
							{
								ComputeProbCond(ind_x, ind_loss, ind_name, ind_sortedname, pk);
								ComputeProbCondDown(ind_x, ind_loss, ind_name, ind_sortedname, pk_down);
							}
							for (ind_loss=max_ldown+1; ind_loss<=max_lup; ind_loss++)
								ComputeProbCond(ind_x, ind_loss, ind_name, ind_sortedname, pk);
						}
						
					}
					
					IsShortSector[its_indSect]=IsShortName;
				}
				ind_GlobalName = ind_GlobalName + ind_name - 1;
				
				//Sector Loss distribution computation : integration par rapport à chaque Vi
				for (ind_x=0; ind_x<nb_fact_cond ; ind_x++)
				{
					if (IsShortName) 
						maxloss = lmin;
					else
						maxloss = max_lup;

					integrate_SectorLoss(maxloss,nb_fact_cond,nb_name_secteur,ind_x,its_indSect);

					if (IsShortName) 
						maxloss = lmin;
					else
						maxloss = max_ldown;

					integrate_SectorLossDown(maxloss,nb_fact_cond,nb_name_secteur,ind_x,its_indSect);

				}
			}

			//Initialisation de its_SectorLoss(itsNbSectors) and its_SectorLoss_Down(itsNbSectors)
			initialize_SectorLoss(lmin,nb_fact_cond);
			initialize_SectorLoss0(lmin,max_lup, nb_fact_cond);
			
			initialize_SectorLossDown(lmin,nb_fact_cond);
			initialize_SectorLoss0Down(lmin,max_ldown, nb_fact_cond);
			
			//Boucle sur les secteurs de manière à estimer la distribution conditionnelle par rapport a V
			for (its_indSect=1; its_indSect<itsNbSectors; its_indSect++)
			{	if (IsShortSector[its_indSect])
				{
					for (ind_x=0; ind_x<nb_fact_cond ; ind_x++)
					{	vector<double> SumElt(lmin+1,0.);
						for (ind_loss=lmin; ind_loss>=0; ind_loss--)
						{				
							for(ind_loss2=lmin; ind_loss2>=ind_loss; ind_loss2--)
								SumElt[ind_loss] += its_SectorLoss->Elt(itsNbSectors,ind_x,ind_loss-ind_loss2+lmin)*its_SectorLoss->Elt(its_indSect,ind_x,ind_loss2);
						
						}
						for (ind_loss=0; ind_loss<=lmin; ind_loss++)
							its_SectorLoss->SetElt(itsNbSectors,ind_x,ind_loss,SumElt[ind_loss]);
						
						for (ind_loss=lmin+1; ind_loss<=max_lup; ind_loss++)
							its_SectorLoss->SetElt(itsNbSectors,ind_x,ind_loss,0.);
						
					}
					for (ind_x=0; ind_x<nb_fact_cond ; ind_x++)
					{	vector<double> SumElt(lmin+1,0.);
						for (ind_loss=lmin; ind_loss>=0; ind_loss--)
						{				
							for(ind_loss2=lmin; ind_loss2>=ind_loss; ind_loss2--)
								SumElt[ind_loss] += its_SectorLoss_Down->Elt(itsNbSectors,ind_x,ind_loss-ind_loss2+lmin)*its_SectorLoss_Down->Elt(its_indSect,ind_x,ind_loss2);
						
						}
						for (ind_loss=0; ind_loss<=lmin; ind_loss++)
							its_SectorLoss_Down->SetElt(itsNbSectors,ind_x,ind_loss,SumElt[ind_loss]);

						for (ind_loss=lmin+1; ind_loss<=max_ldown; ind_loss++)
							its_SectorLoss_Down->SetElt(itsNbSectors,ind_x,ind_loss,0.);
						
					}
				}
				else
				{
					for (ind_x=0; ind_x<nb_fact_cond ; ind_x++)
					{	vector<double> SumElt(max_lup+1,0.);
						for (ind_loss=0; ind_loss<=max_lup; ind_loss++)
						{				
							for(ind_loss2=0; ind_loss2<=ind_loss; ind_loss2++)
								SumElt[ind_loss] += its_SectorLoss->Elt(itsNbSectors,ind_x,ind_loss-ind_loss2)*its_SectorLoss->Elt(its_indSect,ind_x,ind_loss2);
						
						}
						for (ind_loss=0; ind_loss<=max_lup; ind_loss++)
							its_SectorLoss->SetElt(itsNbSectors,ind_x,ind_loss,SumElt[ind_loss]);
					/*FILE *stream = fopen("c:\\temp\\SecLossDistrib.txt", "w+");
					this->ViewS("",stream,ind_x);
					fclose(stream);	*/
					}

					for (ind_x=0; ind_x<nb_fact_cond ; ind_x++)
					{	vector<double> SumElt(max_ldown+1,0.);
						for (ind_loss=0; ind_loss<=max_ldown; ind_loss++)
						{				
							for(ind_loss2=0; ind_loss2<=ind_loss; ind_loss2++)
								SumElt[ind_loss] += its_SectorLoss_Down->Elt(itsNbSectors,ind_x,ind_loss-ind_loss2)*its_SectorLoss_Down->Elt(its_indSect,ind_x,ind_loss2);
						
						}
						for (ind_loss=0; ind_loss<=max_ldown; ind_loss++)
							its_SectorLoss_Down->SetElt(itsNbSectors,ind_x,ind_loss,SumElt[ind_loss]);
					/*FILE *stream = fopen("c:\\temp\\SecLossDistrib.txt", "w+");
					this->ViewS("",stream,ind_x);
					fclose(stream);	*/
					}	
					

				}
			}
					
			//Integration de la distribution totale 
			//Boucle sur chaque loss 
			
			integrate_LossDistrib(max_lup,nb_fact_cond);
			integrate_LossDistribDown(max_ldown,nb_fact_cond);


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


void ICM_Gaussian_LossDistrib_MF::integrate_SectorLoss(const int& lmax, const int& nb_fact_cond, const int& nb_name_secteur, const int& ind_x, const int& indSect)
{
	int ind_loss; 
	int ind_x2;
	double res=0.;

	//Boucle sur les loss
	for (ind_loss=0; ind_loss<=lmax; ind_loss++)
	{	res=0.;
		//Intgération sur le facteur sectoriel Vi
		for (ind_x2=0; ind_x2<nb_fact_cond ; ind_x2++)
			res += its_ProbCond->Elt(nb_name_secteur,ind_x*nb_fact_cond+ind_x2,ind_loss)*TheIntegrator.GetWeight(ind_x2);
		
		if (itsIntegrationMethod	==	qGAUSS_HERMITE)
			its_SectorLoss->SetElt(indSect,ind_x,ind_loss,res*ONEOVERSQRTPI);
		else if (itsIntegrationMethod	==	qGAUSS_LEGENDRE)
			its_SectorLoss->SetElt(indSect,ind_x,ind_loss,res/SQRT2PI);
	}	

}

void ICM_Gaussian_LossDistrib_MF::integrate_SectorLossDown(const int& lmax, const int& nb_fact_cond, const int& nb_name_secteur, const int& ind_x, const int& indSect)
{
	int ind_loss; 
	int ind_x2;
	double resDown=0.;
	
	//Boucle sur les loss
	for (ind_loss=0; ind_loss<=lmax; ind_loss++)
	{	resDown=0.;
		//Intgération sur le facteur sectoriel Vi
		for (ind_x2=0; ind_x2<nb_fact_cond ; ind_x2++)
			resDown += its_ProbCond_Down->Elt(nb_name_secteur,ind_x*nb_fact_cond+ind_x2,ind_loss)*TheIntegrator.GetWeight(ind_x2);
		
		if (itsIntegrationMethod	==	qGAUSS_HERMITE)
			its_SectorLoss_Down->SetElt(indSect,ind_x,ind_loss,resDown*ONEOVERSQRTPI);
		else if (itsIntegrationMethod	==	qGAUSS_LEGENDRE)
			its_SectorLoss_Down->SetElt(indSect,ind_x,ind_loss,resDown/SQRT2PI);
	}

}

void ICM_Gaussian_LossDistrib_MF::integrate_LossDistrib(const int& lmax, const int& nb_fact_cond)
{
	int ind_loss; 
	int ind_x;
	
	for (ind_loss=0; ind_loss<=lmax; ind_loss++) 
	{	
		its_lossdistrib[ind_loss] = 0.;
		//Intégration sur le facteur V 
		for (ind_x=0; ind_x<nb_fact_cond ; ind_x++)
			its_lossdistrib[ind_loss]+=its_SectorLoss->Elt(itsNbSectors,ind_x,ind_loss)*TheIntegrator.GetWeight(ind_x);

		if (itsIntegrationMethod	==	qGAUSS_HERMITE)
			its_lossdistrib[ind_loss]*=ONEOVERSQRTPI;
		else if (itsIntegrationMethod	==	qGAUSS_LEGENDRE)
			its_lossdistrib[ind_loss]/=SQRT2PI;
	}


}

void ICM_Gaussian_LossDistrib_MF::integrate_LossDistribDown(const int& lmax, const int& nb_fact_cond)
{
	int ind_loss; 
	int ind_x;
	
	for (ind_loss=0; ind_loss<=lmax; ind_loss++) 
	{	
		its_lossdistrib_Down[ind_loss] = 0.;
		//Intégration sur le facteur V 
		for (ind_x=0; ind_x<nb_fact_cond ; ind_x++)
			its_lossdistrib_Down[ind_loss]+=its_SectorLoss_Down->Elt(itsNbSectors,ind_x,ind_loss)*TheIntegrator.GetWeight(ind_x);

		if (itsIntegrationMethod	==	qGAUSS_HERMITE)
			its_lossdistrib_Down[ind_loss]*=ONEOVERSQRTPI;
		else if (itsIntegrationMethod	==	qGAUSS_LEGENDRE)
			its_lossdistrib_Down[ind_loss]/=SQRT2PI;
	}

}

void ICM_Gaussian_LossDistrib_MF::initialize_SectorLoss(const int& lmax, const int& nb_fact_cond)
{
	int ind_loss;
	int ind_x;

	for (ind_loss=0; ind_loss<=lmax; ind_loss++)
	{	for (ind_x=0; ind_x<nb_fact_cond ; ind_x++)
			its_SectorLoss->SetElt(itsNbSectors,ind_x,ind_loss,its_SectorLoss->Elt(0,ind_x,ind_loss));
	}
}

void ICM_Gaussian_LossDistrib_MF::initialize_SectorLoss0(const int& lmin, const int& lmax, const int& nb_fact_cond)
{
	int ind_loss;
	int ind_x;

	for (ind_loss=lmin+1; ind_loss<=lmax; ind_loss++)
	{	for (ind_x=0; ind_x<nb_fact_cond ; ind_x++)
			its_SectorLoss->SetElt(itsNbSectors,ind_x,ind_loss,0.);
	}
}

void ICM_Gaussian_LossDistrib_MF::initialize_SectorLossDown(const int& lmax, const int& nb_fact_cond)
{
	int ind_loss;
	int ind_x;

	for (ind_loss=0; ind_loss<=lmax; ind_loss++)
	{	for (ind_x=0; ind_x<nb_fact_cond ; ind_x++)
			its_SectorLoss_Down->SetElt(itsNbSectors,ind_x,ind_loss,its_SectorLoss_Down->Elt(0,ind_x,ind_loss));
	}
}

void ICM_Gaussian_LossDistrib_MF::initialize_SectorLoss0Down(const int& lmin, const int& lmax, const int& nb_fact_cond)
{
	int ind_loss;
	int ind_x;

	for (ind_loss=lmin+1; ind_loss<=lmax; ind_loss++)
	{	for (ind_x=0; ind_x<nb_fact_cond ; ind_x++)
			its_SectorLoss_Down->SetElt(itsNbSectors,ind_x,ind_loss,0.);
	}
}



//Fonctions intermédiaires 
void ICM_Gaussian_LossDistrib_MF::compute_barrier()
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
	}
}

void ICM_Gaussian_LossDistrib_MF::precompute_coeffs()
{
	int k;
	double	den, beta_value, lambda_value;

	for (k=0; k<its_nbnames; k++)
	{
		//Strike up
		if (fabs(beta_value = its_unique_beta[k]) == 1.0)
		{
			throw Exception(__LINE__, __FILE__, ERR_INVALID_DATA,
				"Beta values must be strictly between -1.0 and 1.0");
		}
		
		if (fabs(lambda_value = its_lambda[k]) == 1.0)
		{
			throw Exception(__LINE__, __FILE__, ERR_INVALID_DATA,
				"Lambda values must be strictly between -1.0 and 1.0");
		}

		den = 1.0 / sqrt(1.0 - beta_value * beta_value);
		
		its_coeff_a[k]	=	den;
		its_coeff_b[k]	=	den;
		its_coeff_c[k]	=	den;

		its_coeff_a[k]	*=	its_barrier[k];
		its_coeff_b[k]	*=	beta_value * lambda_value;
		its_coeff_c[k]	*=	beta_value * sqrt(1.0 - lambda_value * lambda_value);
		
		//Strike Down-
		if (fabs(beta_value = its_beta_Down[k]) == 1.0)
		{
			throw Exception(__LINE__, __FILE__, ERR_INVALID_DATA,
				"Beta values must be strictly between -1.0 and 1.0");
		}
		if (fabs(lambda_value = its_lambda_Down[k]) == 1.0)
		{
			throw Exception(__LINE__, __FILE__, ERR_INVALID_DATA,
				"Lambda values must be strictly between -1.0 and 1.0");
		}
		
		den = 1.0 / sqrt(1.0 - beta_value * beta_value);
		
		its_coeff_a_down[k]	=	den;
		its_coeff_b_down[k]	=	den;
		its_coeff_c_down[k]	=	den;
		
		its_coeff_a_down[k]	*=	its_barrier[k];
		its_coeff_b_down[k]	*=	beta_value * lambda_value;
		its_coeff_c_down[k]	*=  beta_value * sqrt(1.0 - lambda_value * lambda_value);
	}
} 
	

void ICM_Gaussian_LossDistrib_MF::UpdateCorrelation(const double& beta_down, 
												  const double& beta_up,
												  const double& lambda_down, 
												  const double& lambda_up,
												  const bool& HedgeFlag											  
													)
{
	int i;

	// maybe to be generalized with vectors of Betas.
	// for each Credit, linear interpolate both Base Correlation Down and Up

	for (i=0; i< its_nbnames; i++)
	{			
		its_unique_beta[i] = beta_up;
		its_beta_Down[i] = beta_down;
		its_lambda[i] = lambda_up;
		its_lambda_Down[i] = lambda_down;

	}

	precompute_coeffs();

	//if (HedgeFlag)
		//precompute_coeffs_perturb();
}


void ICM_Gaussian_LossDistrib_MF::UpdateCorrelation(std::vector<double>& V_beta_down,
												  std::vector<double>& V_beta_up,
												  std::vector<double>& V_lambda_down,
												  std::vector<double>& V_lambda_up, 
												  bool HedgeFlag 												  
												  )
{
	int i = 0;

	// maybe to be generalized with vectors of Betas.
	// for each Credit, linear interpolate both Base Correlation Down and Up

	for (i=0; i< its_nbnames; i++)
	{			
		its_unique_beta[i] = V_beta_up[i];
		its_beta_Down[i] = V_beta_down[i];
		its_lambda[i] = V_lambda_up[i];
		its_lambda_Down[i] = V_lambda_down[i];
	}

	precompute_coeffs();

	//if (HedgeFlag)
		//precompute_coeffs_perturb();
}

void ICM_Gaussian_LossDistrib_MF::SetIntegratorType(const qIntegratorChoice&	TheIntegratorType, 
													const int&	TheStep)
{
	itsIntegrationMethod	=	TheIntegratorType;
	itsIntegrationStep		=	TheStep;

	TheIntegrator.SetIntegrationType(TheIntegratorType);
	
	// Rajout du || sur la méthode qTrapeze afin que le step de l'integrator soit setter a la bonne valeur
	if ((TheIntegratorType == qGAUSS_LEGENDRE) || (TheIntegratorType == qGAUSS_HERMITE) || (TheIntegratorType == qTRAPEZE))
		TheIntegrator.SetIntegrationStep(TheStep);
}

void ICM_Gaussian_LossDistrib_MF::ViewF(char* id, FILE* ficOut,const int& ind_x)
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

    fprintf(fOut, "\t\t\t ----------------- Homogeneous Fast Loss Dsitrib Gaussian MF ----------------- \n");
/*	int ind_loss=0;
	for (ind_loss=0; ind_loss<its_lossdistrib_Down.size(); ind_loss++)
		fprintf(fOut,"ind_x : %i\t%f\t%f\n",ind_loss, its_lossdistrib[ind_loss], its_lossdistrib_Down[ind_loss]);
	for (ind_loss=its_lossdistrib_Down.size(); ind_loss<its_lossdistrib.size(); ind_loss++)
	{
		if (ind_loss==its_lossdistrib_Down.size())
			fprintf(fOut,"ind_x : %i\t%f\t%f\n",ind_loss, its_lossdistrib[ind_loss], its_taildistrib_Down);
		else
			fprintf(fOut,"ind_x : %i\t%f\n",ind_loss, its_lossdistrib[ind_loss]);
	}
	fprintf(fOut,"tail : \t%f\n",its_taildistrib);*/

	fprintf(fOut, "\t\t\t ----------------- Prob Cond  ----------------- \n");

		for (int ind_loss=0; ind_loss<its_lossdistrib.size(); ind_loss++)
		{
			for (int ind_name=0; ind_name<=this->its_nbnames; ind_name++)
				fprintf(fOut,"%f\t",its_ProbCond->Elt(ind_name,ind_x,ind_loss));
			fprintf(fOut,"\n");
		}

		for (ind_loss=0; ind_loss<its_lossdistrib_Down.size(); ind_loss++)
		{
			for (int ind_name=0; ind_name<=this->its_nbnames; ind_name++)
				fprintf(fOut,"%f\t",its_ProbCond_Down->Elt(ind_name,ind_x,ind_loss));
			fprintf(fOut,"\n");
		}
		
		fprintf(fOut,"\n\n");
	

	if ( ficOut == NULL )
	{
		fclose(fOut);
	}
}


void ICM_Gaussian_LossDistrib_MF::ViewS(char* id, FILE* ficOut,const int& ind_x)
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

    fprintf(fOut, "\t\t\t ----------------- Homogeneous Fast Loss Dsitrib Gaussian MF ----------------- \n");
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

	fprintf(fOut, "\t\t\t ----------------- Prob Cond  ----------------- \n");

/*		for (int ind_loss=0; ind_loss<its_lossdistrib.size(); ind_loss++)
		{
			for (int ind_sec=0; ind_sec<=this->itsNbSectors; ind_sec++)
				fprintf(fOut,"%f\t",its_SectorLoss->Elt(ind_sec,ind_x,ind_loss));
			fprintf(fOut,"\n");
		}

		for (ind_loss=0; ind_loss<its_lossdistrib_Down.size(); ind_loss++)
		{
			for (int ind_sec=0; ind_sec<=this->itsNbSectors; ind_sec++)
				fprintf(fOut,"%f\t",its_SectorLoss_Down->Elt(ind_sec,ind_x,ind_loss));
			fprintf(fOut,"\n");
		}		*/

		fprintf(fOut,"\n\n");
	

	if ( ficOut == NULL )
	{
		fclose(fOut);
	}
}
