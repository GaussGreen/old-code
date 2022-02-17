#include "ARMKernel\glob\firsttoinc.h" 
#include "ICMKernel\mod\icm_Homogene_FastLossDistrib_Gaussian.h"
#include "ICMKernel\glob\icm_maths.h"


//----------------------------------------------------------------------------------------
//Collat Forward
//----------------------------------------------------------------------------------------
ICM_Gauss1FLossDistrib_LJ::ICM_Gauss1FLossDistrib_LJ(const int& nbnames,
										   const ARM_Vector& pdef,
										   const ARM_Vector& pdef_start,
										   const ARM_Vector& beta,
										   const ARM_Vector& LossRates,
										   const int& discretizationstep ,
										   const int& CopulaType ,
										   const qIntegratorChoice&	IntegrationMethod ,
										   const int& IntegrationStep)
{
	Init();
	Set(nbnames,pdef,pdef_start,beta,LossRates,discretizationstep,CopulaType,IntegrationMethod,IntegrationStep);

}

//-----------------------------------------------------------------------
//Set Forward collateral
//-----------------------------------------------------------------------
void ICM_Gauss1FLossDistrib_LJ::SetFwdCollat(const int& nbnames,
										   const ARM_Vector&  pdef,
										   const ARM_Vector&  beta,
										   const ARM_Vector& LossRates,
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
		its_coeff_b_perturb.resize(nbnames);
		its_coeff_a_down_perturb.resize(nbnames);
		its_coeff_b_down_perturb.resize(nbnames);

		//fwd collateral
		if (its_collat_fwd_lossdistrib_perturb_Down) delete its_collat_fwd_lossdistrib_perturb_Down;
		its_collat_fwd_lossdistrib_perturb_Down= new ICM_QMatrix<double>(1,nbnames);

		its_collat_fwd_taildistrib_perturb_Down.resize(nbnames);
		its_collat_fwd_beta_Down.resize(nbnames);
		its_collat_fwd_coeff_a.resize(nbnames);
		its_collat_fwd_coeff_b.resize(nbnames);
		its_collat_fwd_coeff_a_down.resize(nbnames);
		its_collat_fwd_coeff_b_down.resize(nbnames);
		its_collat_fwd_coeff_a_perturb.resize(nbnames);
		its_collat_fwd_coeff_a_down_perturb.resize(nbnames);
	}

	if ((nbnames != its_nbnames) || (itsIntegrationStep != size) || 
		// 17783  (!its_ProbCond)||(!its_ProbCond_Perturb)||(!its_ProbCond_Perturb_Down)||(!its_ProbCond_Down))
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


//---------------------------------------------------------------------------
//Collat Fwd
//---------------------------------------------------------------------------
void ICM_Gauss1FLossDistrib_LJ::Set(const int& nbnames,
										   const ARM_Vector&  pdef,
										   const ARM_Vector&  pdef_start,
										   const ARM_Vector&  beta,
										   const ARM_Vector& LossRates,
										   const int& discretizationstep ,
										   const int& CopulaType ,
										   const qIntegratorChoice&	IntegrationMethod ,
										   const int& IntegrationStep)
{
	SetFwdCollat(nbnames,pdef,beta,LossRates,discretizationstep ,CopulaType ,IntegrationMethod ,IntegrationStep);
	//MaJ des Pdef Start dans le vecteur PdefPerturb
	SetPdefPerturb(pdef_start);
	//Estimation des barrières correspondantes
	compute_barrier_perturb();
}

void ICM_Gauss1FLossDistrib_LJ::Collat_fwd_UpdateCorrelation(std::vector<double>& V_beta_down,
															std::vector<double>& V_beta_up, 
															bool HedgeFlag)
{
	int i = 0;
	SetFwdCollatTScase(true);

	// maybe to be generalized with vectors of Betas.
	// for each Credit, linear interpolate both Base Correlation Down and Up

	its_collat_fwd_unique_beta.resize(its_nbnames);
	its_collat_fwd_beta_Down.resize(its_nbnames);

	for (i=0; i< its_nbnames; i++)
	{			
		its_collat_fwd_unique_beta[i] = V_beta_up[i];
		its_collat_fwd_beta_Down[i] = V_beta_down[i];
	}

	precompute_coeffs();

	if (HedgeFlag)
		precompute_coeffs_perturb();
}

//----------------------------------------------//
//				Forward Collateral
//----------------------------------------------//

//Utilities
extern "C" void 
InitProbCondPortSize0CollatFwd(void* Param, double x, double& res)
{
	int k;
	double	pk=0.;
	double	tmp_barrier=0., tmp_barrier_start=0.;

	double CondProbaZeroLoss	=	1.;

	// get parameters from stack
	ICM_Gauss1FLossDistrib_LJ*	TheModel;
	TheModel	=	(ICM_Gauss1FLossDistrib_LJ*)((*(AddressVector*)Param)[0]);

	int*	its_ind_x;
	its_ind_x	=	(int*)((*(AddressVector*)Param)[1]);
	(*its_ind_x) += 1;

	int itsIntegrationMethod = TheModel->GetIntegrationMethod();

	if (itsIntegrationMethod	==	qGAUSS_HERMITE)
		x	*=	SQRT2;

	TheModel->SetProbCond_Elt(0, *its_ind_x,0,CondProbaZeroLoss);

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

		//case forward collateral in TS mode
		if (ICM_Distribution::itsFwdCollatTScase)
			tmp_barrier_start	=	TheModel->Get_collat_fwdCoeff_Perturb_a(k) - TheModel->Get_collat_fwdCoeff_b(k) * x;
		else
			tmp_barrier_start	=	TheModel->GetCoeff_Perturb_a(k) - TheModel->GetCoeff_Perturb_b(k) * x;

		pk -= NAG_cumul_normal(tmp_barrier_start);

		if (TheModel->IsTSR() && TheModel->GetCollatT2())
		{
			double tmp_barrier_t2=0.;
			double tmp_barrier_end=0.;

			tmp_barrier_t2	=	TheModel->GetCoeff_Perturb_a_t2(k) - TheModel->GetCoeff_Perturb_b_t2(k) * x;
			pk -= NAG_cumul_normal(tmp_barrier_t2);
			pk +=	NAG_bivariate_normal_dist(tmp_barrier_start,tmp_barrier_t2,RHO_DIST(TheModel->GetCollatT1(),TheModel->GetCollatT2()));
		}

		if (pk<0.) pk=0.;

		CondProbaZeroLoss *= (1.-pk);
		TheModel->SetProbCond_Elt(k+1, *its_ind_x,0,CondProbaZeroLoss); 
	}

	res		=	CondProbaZeroLoss;		
}

extern "C" void 
InitProbCondPortSize0CollatFwdSmile(void* Param, double x, double& ProbCond_Up, double& ProbCond_Down)
{
	int k;
	double	pk=0.;
	double	tmp_barrier=0., tmp_barrier_start=0.;

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

		//case forward collateral in TS mode
		if (ICM_Distribution::itsFwdCollatTScase)
			tmp_barrier_start	=	TheModel->Get_collat_fwdCoeff_Perturb_a(k) - TheModel->Get_collat_fwdCoeff_b(k) * x;
		else
			tmp_barrier_start	=	TheModel->GetCoeff_Perturb_a(k) - TheModel->GetCoeff_Perturb_b(k) * x;

		pk -= NAG_cumul_normal(tmp_barrier_start);

		if (TheModel->IsTSR() && TheModel->GetCollatT2())
		{
			double tmp_barrier_t2=0.;
			double tmp_barrier_end=0.;

			tmp_barrier_t2	=	TheModel->GetCoeff_Perturb_a_t2(k) - TheModel->GetCoeff_Perturb_b_t2(k) * x;
			pk -= NAG_cumul_normal(tmp_barrier_t2);

			pk +=	NAG_bivariate_normal_dist(tmp_barrier_start,tmp_barrier_t2,RHO_DIST(TheModel->GetCollatT1(),TheModel->GetCollatT2()));
		}

		if (pk<0.) pk=0.;

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

		//case forward collateral in TS mode
		if (ICM_Distribution::itsFwdCollatTScase)
			tmp_barrier_start =	TheModel->Get_collat_fwdCoeff_Perturb_a_down(k) - TheModel->Get_collat_fwdCoeff_b_down(k) * x;
		else
			tmp_barrier_start =	TheModel->GetCoeff_Perturb_a_down(k) - TheModel->GetCoeff_Perturb_b_down(k) * x;

		pk -= NAG_cumul_normal(tmp_barrier_start);

		if (TheModel->IsTSR() && TheModel->GetCollatT2())
		{
			double tmp_barrier_t2=0.;
			double tmp_barrier_end=0.;

			tmp_barrier_t2	=	TheModel->GetCoeff_Perturb_a_down_t2(k) - TheModel->GetCoeff_Perturb_b_down_t2(k) * x;
			pk -= NAG_cumul_normal(tmp_barrier_t2);

			pk +=	NAG_bivariate_normal_dist(tmp_barrier_start,tmp_barrier_t2,RHO_DIST(TheModel->GetCollatT1(),TheModel->GetCollatT2()));
		}

		if (pk<0.) pk=0.;

		CondProbaZeroLoss_Down *= (1.-pk);
		TheModel->SetProbCond_Down_Elt(k+1, *its_ind_x,0,CondProbaZeroLoss_Down);
	}

	ProbCond_Up		=	CondProbaZeroLoss;	
	ProbCond_Down	=	CondProbaZeroLoss_Down;	
}

extern "C" void 
ComputeProbCondCollatFwd(void* Param, double x, double& res)
{
	// get parameters from stack
	ICM_Gauss1FLossDistrib_LJ*	TheModel;
	TheModel	=	(ICM_Gauss1FLossDistrib_LJ*)((*(AddressVector*)Param)[0]);

	int*	its_ind_x;
	its_ind_x	=	(int*)((*(AddressVector*)Param)[1]);
	(*its_ind_x) += 1;

	int*	its_ind_Loss;
	its_ind_Loss	=	(int*)((*(AddressVector*)Param)[2]);

	int	k;
	double tmp1; // proba de loss=L pour un panier de taille k
	double tmp2; // proba de loss=L-1 pour un panier de taille k
	double pk;     // proba cond de défaut du nom k
	int ind_lastloss;

	double tmp_barrier, tmp_barrier_start;

	int itsIntegrationMethod = TheModel->GetIntegrationMethod();
	if (itsIntegrationMethod	==	qGAUSS_HERMITE)
		x	*=	SQRT2;

	TheModel->SetProbCond_Elt(0,*its_ind_x,*its_ind_Loss,0.);

	// its_ind_Loss: Loss Distribution for a basket of size (0) knowing factor value its_ind_x
	tmp1		=	TheModel->GetProbCond_Elt(0, *its_ind_x, *its_ind_Loss);

	if (TheModel->IsHomogeneous())
	{
		for (k=0;k<TheModel->GetNbNames();k++) 
		{
			// barrier : Les pdef associées à la partie Start sont stockées dans les objets peerturb
			tmp_barrier = TheModel->GetCoeff_a(k) - TheModel->GetCoeff_b(k) * x;
			pk=NAG_cumul_normal(tmp_barrier);

			if (TheModel->IsTSR() && TheModel->GetT2())
			{
				double tmp_barrier_t2=0.;
				double tmp_barrier_end=0.;

				tmp_barrier_t2	=	TheModel->GetCoeff_a_t2(k) - TheModel->GetCoeff_b_t2(k) * x;
				pk += NAG_cumul_normal(tmp_barrier_t2);

				pk -=	NAG_bivariate_normal_dist(tmp_barrier,tmp_barrier_t2,RHO_DIST(TheModel->GetT1(),TheModel->GetT2()));
			}

			//case forward collateral in TS mode
			if (ICM_Distribution::itsFwdCollatTScase)
				tmp_barrier_start = TheModel->Get_collat_fwdCoeff_Perturb_a(k) - TheModel->Get_collat_fwdCoeff_b(k) * x;
			else
				tmp_barrier_start = TheModel->GetCoeff_Perturb_a(k) - TheModel->GetCoeff_Perturb_b(k) * x;

			// conditional default probability
			pk -= NAG_cumul_normal(tmp_barrier_start);

			if (TheModel->IsTSR() && TheModel->GetCollatT2())
			{
				double tmp_barrier_t2=0.;
				double tmp_barrier_end=0.;

				tmp_barrier_t2	=	TheModel->GetCoeff_Perturb_a_t2(k) - TheModel->GetCoeff_Perturb_b_t2(k) * x;
				pk -= NAG_cumul_normal(tmp_barrier_t2);
				pk +=	NAG_bivariate_normal_dist(tmp_barrier_start,tmp_barrier_t2,RHO_DIST(TheModel->GetCollatT1(),TheModel->GetCollatT2()));
			}

			if (pk<0.) pk = 0.;

			// its_ind_Loss: Loss Distribution for a basket of size (k-1) knowing factor value its_ind_x
			// tmp1

			// Loss Distribution for a basket of size k knowing factor value its_ind_x
			tmp2		=	tmp1 * (1.-pk) + pk * TheModel->GetProbCond_Elt(k,*its_ind_x,*its_ind_Loss-1);
			TheModel->SetProbCond_Elt(k+1, *its_ind_x, *its_ind_Loss, tmp2);
			tmp1 = tmp2;
		}
	}
	else
	{
		for (k=0;k<TheModel->GetNbNames();k++) 
		{
			// current level of loss - loss rate of names k
			// double LossRates_inf = floor(TheModel->GetIndLossRates()[k]);
			double LossRates_inf = TheModel->GetIntLossRates()[k];
			double LossRates_delta = TheModel->GetDblLossRates()[k] - LossRates_inf;
			
			// ind_lastloss = *its_ind_Loss - (int)LossRates_inf;
			ind_lastloss = *its_ind_Loss - TheModel->GetIntLossRates()[k];

			// barrier : Les pdef associées à la partie Start sont stockées dans les objets peerturb
			tmp_barrier = TheModel->GetCoeff_a(k) - TheModel->GetCoeff_b(k) * x;
			pk=NAG_cumul_normal(tmp_barrier);

			if (TheModel->IsTSR() && TheModel->GetT2())
			{
				double tmp_barrier_t2=0.;
				double tmp_barrier_end=0.;

				tmp_barrier_t2	=	TheModel->GetCoeff_a_t2(k) - TheModel->GetCoeff_b_t2(k) * x;
				pk += NAG_cumul_normal(tmp_barrier_t2);

				pk -=	NAG_bivariate_normal_dist(tmp_barrier,tmp_barrier_t2,RHO_DIST(TheModel->GetT1(),TheModel->GetT2()));
			}

			//case forward collateral in TS mode
			if (ICM_Distribution::itsFwdCollatTScase)
				tmp_barrier_start = TheModel->Get_collat_fwdCoeff_Perturb_a(k) - TheModel->Get_collat_fwdCoeff_b(k) * x;
			else
				tmp_barrier_start = TheModel->GetCoeff_Perturb_a(k) - TheModel->GetCoeff_Perturb_b(k) * x;

			// conditional default probability
			pk -= NAG_cumul_normal(tmp_barrier_start);

			if (TheModel->IsTSR() && TheModel->GetCollatT2())
			{
				double tmp_barrier_t2=0.;
				double tmp_barrier_end=0.;

				tmp_barrier_t2	=	TheModel->GetCoeff_Perturb_a_t2(k) - TheModel->GetCoeff_Perturb_b_t2(k) * x;
				pk -= NAG_cumul_normal(tmp_barrier_t2);
				pk +=	NAG_bivariate_normal_dist(tmp_barrier_start,tmp_barrier_t2,RHO_DIST(TheModel->GetCollatT1(),TheModel->GetCollatT2()));
			}

			if (pk<0.) pk = 0.;

			// its_ind_Loss: Loss Distribution for a basket of size (k-1) knowing factor value its_ind_x
			// tmp1

			// Loss Distribution for a basket of size k knowing factor value its_ind_x
			tmp2		=	tmp1 * (1.-pk);
			
			if (ind_lastloss >= 0)
			{
				tmp2	+=	pk * TheModel->GetProbCond_Elt(k, *its_ind_x, ind_lastloss)* (1.-LossRates_delta);

				if (ind_lastloss > 0)
				{	
						tmp2	+=	pk * TheModel->GetProbCond_Elt(k, *its_ind_x, ind_lastloss-1)* LossRates_delta;
				}
			}
			// else: the loss of name k is too big
			TheModel->SetProbCond_Elt(k+1, *its_ind_x, *its_ind_Loss, tmp2);
			tmp1 = tmp2;
		}
	}
	res		=	 tmp1;
}

extern "C" void 
ComputeProbCondCollatFwdSmile(void* Param, double x, double& ProbCond_Up, double& ProbCond_Down)
{
	// get parameters from stack
	ICM_Gauss1FLossDistrib_LJ*	TheModel;
	TheModel	=	(ICM_Gauss1FLossDistrib_LJ*)((*(AddressVector*)Param)[0]);

	int*	its_ind_x;
	its_ind_x	=	(int*)((*(AddressVector*)Param)[1]);
	(*its_ind_x) += 1;

	int*	its_ind_Loss;
	its_ind_Loss	=	(int*)((*(AddressVector*)Param)[2]);

	int	k;
	double tmp1, tmp1_down; // proba de loss=L pour un panier de taille k
	double tmp2, tmp2_down; // proba de loss=L-1 pour un panier de taille k
	double pk, pk_down;     // proba cond de défaut du nom k
	int ind_lastloss;

	double tmp_barrier, tmp_barrier_down,tmp_barrier_start, tmp_barrier_down_start;

	int itsIntegrationMethod = TheModel->GetIntegrationMethod();
	if (itsIntegrationMethod	==	qGAUSS_HERMITE)
		x	*=	SQRT2;

	TheModel->SetProbCond_Elt(0,*its_ind_x,*its_ind_Loss,0.);
	TheModel->SetProbCond_Down_Elt(0,*its_ind_x,*its_ind_Loss,0.);

	// its_ind_Loss: Loss Distribution for a basket of size (0) knowing factor value its_ind_x
	tmp1		=	TheModel->GetProbCond_Elt(0, *its_ind_x, *its_ind_Loss);
	tmp1_down	=	TheModel->GetProbCond_Down_Elt(0, *its_ind_x, *its_ind_Loss);

	if (TheModel->IsHomogeneous())
	{
		for (k=0;k<TheModel->GetNbNames();k++) 
		{
			// barrier : Les pdef associées à la partie Start sont stockées dans les objets peerturb
			tmp_barrier = TheModel->GetCoeff_a(k) - TheModel->GetCoeff_b(k) * x;
			pk=NAG_cumul_normal(tmp_barrier);

			if (TheModel->IsTSR() && TheModel->GetT2())
			{
				double tmp_barrier_t2=0.;
				double tmp_barrier_end=0.;

				tmp_barrier_t2	=	TheModel->GetCoeff_a_t2(k) - TheModel->GetCoeff_b_t2(k) * x;
				pk += NAG_cumul_normal(tmp_barrier_t2);

				pk -=	NAG_bivariate_normal_dist(tmp_barrier,tmp_barrier_t2,RHO_DIST(TheModel->GetT1(),TheModel->GetT2()));
			}

			//case forward collateral in TS mode
			if (ICM_Distribution::itsFwdCollatTScase)
				tmp_barrier_start = TheModel->Get_collat_fwdCoeff_Perturb_a(k) - TheModel->Get_collat_fwdCoeff_b(k) * x;
			else
				tmp_barrier_start = TheModel->GetCoeff_Perturb_a(k) - TheModel->GetCoeff_Perturb_b(k) * x;

			pk -= NAG_cumul_normal(tmp_barrier_start);

			if (TheModel->IsTSR() && TheModel->GetCollatT2())
			{
				double tmp_barrier_t2=0.;
				double tmp_barrier_end=0.;

				tmp_barrier_t2	=	TheModel->GetCoeff_Perturb_a_t2(k) - TheModel->GetCoeff_Perturb_b_t2(k) * x;
				pk -= NAG_cumul_normal(tmp_barrier_t2);

				pk +=	NAG_bivariate_normal_dist(tmp_barrier_start,tmp_barrier_t2,RHO_DIST(TheModel->GetCollatT1(),TheModel->GetCollatT2()));
			}

			tmp_barrier_down = TheModel->GetCoeff_a_down(k) - TheModel->GetCoeff_b_down(k) * x;
			pk_down=NAG_cumul_normal(tmp_barrier_down);

			if (TheModel->IsTSR() && TheModel->GetT2())
			{
				double tmp_barrier_t2=0.;
				double tmp_barrier_end=0.;

				tmp_barrier_t2	=	TheModel->GetCoeff_a_down_t2(k) - TheModel->GetCoeff_b_down_t2(k) * x;
				pk_down += NAG_cumul_normal(tmp_barrier_t2);

				pk_down -=	NAG_bivariate_normal_dist(tmp_barrier_down,tmp_barrier_t2,RHO_DIST(TheModel->GetT1(),TheModel->GetT2()));
			}

			//case forward collateral in TS mode
			if (ICM_Distribution::itsFwdCollatTScase)
				tmp_barrier_down_start = TheModel->Get_collat_fwdCoeff_Perturb_a_down(k) - TheModel->Get_collat_fwdCoeff_b_down(k) * x;
			else
				tmp_barrier_down_start = TheModel->GetCoeff_Perturb_a_down(k) - TheModel->GetCoeff_Perturb_b_down(k) * x;

			// conditional default probability
			pk_down -= NAG_cumul_normal(tmp_barrier_down_start);

			if (TheModel->IsTSR() && TheModel->GetCollatT2())
			{
				double tmp_barrier_t2=0.;
				double tmp_barrier_end=0.;

				tmp_barrier_t2	=	TheModel->GetCoeff_Perturb_a_down_t2(k) - TheModel->GetCoeff_Perturb_b_down_t2(k) * x;
				pk_down -= NAG_cumul_normal(tmp_barrier_t2);

				pk_down +=	NAG_bivariate_normal_dist(tmp_barrier_down_start,tmp_barrier_t2,RHO_DIST(TheModel->GetCollatT1(),TheModel->GetCollatT2()));
			}

			if (pk<0.) pk = 0.;
			if (pk_down<0.) pk_down = 0.;

			// its_ind_Loss: Loss Distribution for a basket of size (k-1) knowing factor value its_ind_x
			// tmp1

			// Loss Distribution for a basket of size k knowing factor value its_ind_x
			tmp2		=	tmp1 * (1.-pk) + pk * TheModel->GetProbCond_Elt(k,*its_ind_x,*its_ind_Loss-1);
			tmp2_down	=	tmp1_down * (1.-pk_down) + pk_down * TheModel->GetProbCond_Down_Elt(k,*its_ind_x,*its_ind_Loss-1);

			TheModel->SetProbCond_Elt(k+1, *its_ind_x, *its_ind_Loss, tmp2);
			TheModel->SetProbCond_Down_Elt(k+1, *its_ind_x, *its_ind_Loss, tmp2_down);

			tmp1 = tmp2;
			tmp1_down = tmp2_down;
		}
	}
	else
	{
		for (k=0;k<TheModel->GetNbNames();k++) 
		{
			// current level of loss - loss rate of names k
			double LossRates_inf = TheModel->GetIntLossRates()[k];
			double LossRates_delta = TheModel->GetDblLossRates()[k] - LossRates_inf;
			
			ind_lastloss = *its_ind_Loss - TheModel->GetIntLossRates()[k];

			// barrier : Les pdef associées à la partie Start sont stockées dans les objets peerturb
			tmp_barrier = TheModel->GetCoeff_a(k) - TheModel->GetCoeff_b(k) * x;
			pk=NAG_cumul_normal(tmp_barrier);

			if (TheModel->IsTSR() && TheModel->GetT2())
			{
				double tmp_barrier_t2=0.;
				double tmp_barrier_end=0.;

				tmp_barrier_t2	=	TheModel->GetCoeff_a_t2(k) - TheModel->GetCoeff_b_t2(k) * x;
				pk += NAG_cumul_normal(tmp_barrier_t2);

				pk -=	NAG_bivariate_normal_dist(tmp_barrier,tmp_barrier_t2,RHO_DIST(TheModel->GetT1(),TheModel->GetT2()));
			}

			//case forward collateral in TS mode
			if (ICM_Distribution::itsFwdCollatTScase)
				tmp_barrier_start = TheModel->Get_collat_fwdCoeff_Perturb_a(k) - TheModel->Get_collat_fwdCoeff_b(k) * x;
			else
				tmp_barrier_start = TheModel->GetCoeff_Perturb_a(k) - TheModel->GetCoeff_Perturb_b(k) * x;
			pk -= NAG_cumul_normal(tmp_barrier_start);

			if (TheModel->IsTSR() && TheModel->GetCollatT2())
			{
				double tmp_barrier_t2=0.;
				double tmp_barrier_end=0.;

				tmp_barrier_t2	=	TheModel->GetCoeff_Perturb_a_t2(k) - TheModel->GetCoeff_Perturb_b_t2(k) * x;
				pk -= NAG_cumul_normal(tmp_barrier_t2);

				pk +=	NAG_bivariate_normal_dist(tmp_barrier_start,tmp_barrier_t2,RHO_DIST(TheModel->GetCollatT1(),TheModel->GetCollatT2()));
			}

			tmp_barrier_down = TheModel->GetCoeff_a_down(k) - TheModel->GetCoeff_b_down(k) * x;
			pk_down=NAG_cumul_normal(tmp_barrier_down);

			if (TheModel->IsTSR() && TheModel->GetT2())
			{
				double tmp_barrier_t2=0.;
				double tmp_barrier_end=0.;

				tmp_barrier_t2	=	TheModel->GetCoeff_a_down_t2(k) - TheModel->GetCoeff_b_down_t2(k) * x;
				pk_down += NAG_cumul_normal(tmp_barrier_t2);

				pk_down -=	NAG_bivariate_normal_dist(tmp_barrier_down,tmp_barrier_t2,RHO_DIST(TheModel->GetT1(),TheModel->GetT2()));
			}

			//case forward collateral in TS mode
			if (ICM_Distribution::itsFwdCollatTScase)
				tmp_barrier_down_start = TheModel->Get_collat_fwdCoeff_Perturb_a_down(k) - TheModel->Get_collat_fwdCoeff_b_down(k) * x;
			else
				tmp_barrier_down_start = TheModel->GetCoeff_Perturb_a_down(k) - TheModel->GetCoeff_Perturb_b_down(k) * x;

			// conditional default probability
			pk_down -= NAG_cumul_normal(tmp_barrier_down_start);

			if (TheModel->IsTSR() && TheModel->GetCollatT2())
			{
				double tmp_barrier_t2=0.;
				double tmp_barrier_end=0.;

				tmp_barrier_t2	=	TheModel->GetCoeff_Perturb_a_down_t2(k) - TheModel->GetCoeff_Perturb_b_down_t2(k) * x;
				pk_down -= NAG_cumul_normal(tmp_barrier_t2);

				pk_down +=	NAG_bivariate_normal_dist(tmp_barrier_down_start,tmp_barrier_t2,RHO_DIST(TheModel->GetCollatT1(),TheModel->GetCollatT2()));
			}

			if (pk<0.) pk = 0.;
			if (pk_down<0.) pk_down = 0.;

			// its_ind_Loss: Loss Distribution for a basket of size (k-1) knowing factor value its_ind_x
			// tmp1

			// Loss Distribution for a basket of size k knowing factor value its_ind_x
			tmp2		=	tmp1 * (1.-pk);
			tmp2_down	=	tmp1_down * (1.-pk_down);
			
			if (ind_lastloss >= 0)
			{
				tmp2	+=	pk * TheModel->GetProbCond_Elt(k, *its_ind_x, ind_lastloss)* (1.-LossRates_delta);
				tmp2_down	+=	pk_down * TheModel->GetProbCond_Down_Elt(k, *its_ind_x, ind_lastloss)* (1.-LossRates_delta);

				if (ind_lastloss > 0)
				{	
					tmp2	+=	pk * TheModel->GetProbCond_Elt(k, *its_ind_x, ind_lastloss-1)* (LossRates_delta);
					tmp2_down	+=	pk_down * TheModel->GetProbCond_Down_Elt(k, *its_ind_x, ind_lastloss-1)* (LossRates_delta);
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

//Expected Loss

double ICM_Gauss1FLossDistrib_LJ::compute_expectedlosstranche_CollatFwd(const double& tranche_up, 
												  const double& tranche_down, 
												  const double& lossunit
												  //ICM_QMatrix<double>* ShiftMatrix,
												  // const int& Tenor,
												  // const int& Issuer
												  )
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
		compute_distrib_CollatFwd(lup, ldown);
	else
		compute_distrib_CollatFwd(lup);
	
	//View the loss ditribution
	/*#ifdef _DEBUG
	//FILE *stream = fopen("c:\\temp\\LossDistrib.txt", "w+");
	//this->View("",stream);
	//fclose(stream);
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


//Equity Distribution
void ICM_Gauss1FLossDistrib_LJ::compute_distrib_CollatFwd(const int& lup)
{
	double tmpPrbLoss=0.,cumul_distrib=0.;
	int min_lup =0,max_lup =0;
	int l = 0;

	its_lup		=	0;
	max_lup		=	lup;
	
	// Parameters
	AddressVector	TheParameterVector;

	TheParameterVector.Append(this);
	TheParameterVector.Append(&its_ind_x);
	TheParameterVector.Append(&its_ind_Loss);

	int nointegrationcoef=-1;
	TheParameterVector.Append(&nointegrationcoef);
	// ------------------------

	switch (itsCopulaType)
	{
		case qNO_COPULA:	
		case qGAUSSIAN:

			for (l=its_lup;l<=max_lup;l++) 
			{
				if (l ==0)
				{
					its_ind_x=-1;
					TheIntegrator.Integrate(-6.0, 6.0, InitProbCondPortSize0CollatFwd, &TheParameterVector, tmpPrbLoss); 

					its_lossdistrib[0] = tmpPrbLoss;
					its_ind_Loss=1;

				}
				else
				{
					its_ind_x=-1;
					TheIntegrator.Integrate(-6.0, 6.0, ComputeProbCondCollatFwd, &TheParameterVector, tmpPrbLoss); 
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

//Mezzanine Distribution
void ICM_Gauss1FLossDistrib_LJ::compute_distrib_CollatFwd(const int& lup, const int& ldown)
{
	double tmpPrbLoss=0.,cumul_distrib=0.;
	double tmpPrbLoss_down=0.,cumul_distrib_down=0.;
	int min_lup =0,max_lup =0;
	int l = 0;

	its_lup	=	0;
	
	// Parameters
	AddressVector	TheParameterVector;

	TheParameterVector.Append(this);
	TheParameterVector.Append(&its_ind_x);
	TheParameterVector.Append(&its_ind_Loss);

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

					TheIntegrator.Integrate(-6.0, 6.0, InitProbCondPortSize0CollatFwdSmile, &TheParameterVector, tmpPrbLoss, tmpPrbLoss_down); 

					its_lossdistrib[0] = tmpPrbLoss;
					its_lossdistrib_Down[0] = tmpPrbLoss_down;
					its_ind_Loss=1;
				}
				else
				{
					its_ind_x=-1;

					TheIntegrator.Integrate(-6.0, 6.0, ComputeProbCondCollatFwdSmile, &TheParameterVector, tmpPrbLoss, tmpPrbLoss_down); 

					its_lossdistrib[l] = tmpPrbLoss;
					its_lossdistrib_Down[l] = tmpPrbLoss_down;

					its_ind_Loss+=1;
				}
			}
			for (l=ldown+1;l<=lup;l++) 
			{
				its_ind_x=-1;

				TheIntegrator.Integrate(-6.0, 6.0, ComputeProbCondCollatFwd, &TheParameterVector, tmpPrbLoss); 

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

//-----------------------------------------------------------------------
// Forward collateral : Term Structure Review
//-----------------------------------------------------------------------
void ICM_Gauss1FLossDistrib_LJ::SetFwdCollat_TSR(const int& nbnames,
													 const ARM_Vector&  pdef,
													 const ARM_Vector&  pdef_t2,
													 const ARM_Vector&  beta_t1,
													 const ARM_Vector&  beta_t2,
													 const ARM_Vector&  start_beta_t1,
													 const ARM_Vector&  start_beta_t2,
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

		its_beta_Down_t2.resize(nbnames);
		its_coeff_a_t2.resize(nbnames);
		its_coeff_b_t2.resize(nbnames);
		its_coeff_a_down_t2.resize(nbnames);
		its_coeff_b_down_t2.resize(nbnames);
		its_coeff_a_perturb_t2.resize(nbnames);
		its_coeff_a_down_perturb_t2.resize(nbnames);

		//fwd collateral
		if (its_collat_fwd_lossdistrib_perturb_Down) delete its_collat_fwd_lossdistrib_perturb_Down;
		its_collat_fwd_lossdistrib_perturb_Down= new ICM_QMatrix<double>(1,nbnames);

		its_collat_fwd_taildistrib_perturb_Down.resize(nbnames);
		its_collat_fwd_beta_Down.resize(nbnames);
		its_collat_fwd_coeff_a.resize(nbnames);
		its_collat_fwd_coeff_b.resize(nbnames);
		its_collat_fwd_coeff_a_down.resize(nbnames);
		its_collat_fwd_coeff_b_down.resize(nbnames);
		its_collat_fwd_coeff_a_perturb.resize(nbnames);
		its_collat_fwd_coeff_a_down_perturb.resize(nbnames);

		its_beta_Down_t2.resize(nbnames);
		its_collat_fwd_beta_Down_t2.resize(nbnames);
		its_coeff_b_perturb.resize(nbnames);
		its_coeff_b_down_perturb.resize(nbnames);
		its_coeff_a_t2.resize(nbnames);
		its_coeff_b_t2.resize(nbnames);
		its_coeff_a_down_t2.resize(nbnames);
		its_coeff_b_down_t2.resize(nbnames);
		its_coeff_a_perturb_t2.resize(nbnames);
		its_coeff_a_down_perturb_t2.resize(nbnames);
		its_coeff_b_perturb_t2.resize(nbnames);
		its_coeff_b_down_perturb_t2.resize(nbnames);	
	}

	if ((nbnames != its_nbnames) || (itsIntegrationStep != size) || 
		// 17783 (!its_ProbCond)||(!its_ProbCond_Perturb)||(!its_ProbCond_Perturb_Down)||(!its_ProbCond_Down))
		(!its_ProbCond)||(!its_ProbCond_Perturb_Down)||(!its_ProbCond_Down))
	{
		if (its_ProbCond) delete its_ProbCond;
		its_ProbCond= new ICM_QCubix<double>(nbnames+1,size,1,0.);

		// 17783	if (its_ProbCond_Perturb) delete its_ProbCond_Perturb;
		// 17783	its_ProbCond_Perturb= new ICM_QCubix<double>(nbnames+1,size,1,0.);

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
	SetUniqueBeta_t2(beta_t2); 

	SetCollatUniqueBeta(start_beta_t1); 
	SetCollatUniqueBeta_t2(start_beta_t2); 

	SetPdefAtMaturity(pdef); 
	SetPdefAtMaturity_t2(pdef_t2); 

	compute_min_pdef(GetPdefAtMaturity());

	SetIntLossRates(LossRates);
	check_homogeneous(LossRates);
	compute_barrier();		// Copula dependant

	its_lossdistrib_Down.resize(1);
}

//-----------------------------------------------------------------------
// Forward collateral : Term Structure Review
//-----------------------------------------------------------------------

void ICM_Gauss1FLossDistrib_LJ::Collat_fwd_UpdateCorrelation_TSR(const ARM_Vector & V_beta_down,
																 const ARM_Vector & V_beta_down_t2,
																 const ARM_Vector & V_beta_up,
																 const ARM_Vector & V_beta_up_t2,
																 const ARM_Vector & V_collat_beta_down,
																 const ARM_Vector & V_collat_beta_down_t2,
																 const ARM_Vector & V_collat_beta_up,
																 const ARM_Vector & V_collat_beta_up_t2)
{
	int i = 0;
	//SetFwdCollatTScase(true);

	UpdateCorrelation_TSR(V_beta_down,V_beta_down_t2,V_beta_up,V_beta_up_t2);

	// maybe to be generalized with vectors of Betas.
	// for each Credit, linear interpolate both Base Correlation Down and Up

	its_collat_fwd_unique_beta.resize(its_nbnames);
	its_collat_fwd_beta_Down.resize(its_nbnames);

	its_collat_fwd_unique_beta_t2.resize(its_nbnames);
	its_collat_fwd_beta_Down_t2.resize(its_nbnames);


	for (i=0; i< its_nbnames; i++)
	{			
		its_collat_fwd_unique_beta[i] = V_collat_beta_up[i];
		its_collat_fwd_beta_Down[i] = V_collat_beta_down[i];
		its_collat_fwd_unique_beta_t2[i] = V_collat_beta_up_t2[i];
		its_collat_fwd_beta_Down_t2[i] = V_collat_beta_down_t2[i];
	}

	precompute_coeffs_perturb();

}


//-----------------------------------------------------------------------
// Forward collateral : Term Structure Review
//-----------------------------------------------------------------------
ICM_Gauss1FLossDistrib_LJ::ICM_Gauss1FLossDistrib_LJ(const double& t1_corr,
													 const double& t2_corr,
													 const double& start_t1_corr,
													 const double& start_t2_corr,
													 const int& nbnames,
													 const ARM_Vector& pdef,
													 const ARM_Vector& pdef_t2,
													 const ARM_Vector& pdef_start,
													 const ARM_Vector& pdef_start_t2,
													 const ARM_Vector& beta_t1,
													 const ARM_Vector& beta_t2,
													 const ARM_Vector& start_beta_t1,
													 const ARM_Vector& start_beta_t2,
													 // const std::vector<double>& LossRates,
													 const ARM_Vector& LossRates,
													 const int& discretizationstep ,
													 const int& CopulaType ,
													 const qIntegratorChoice&	IntegrationMethod ,
													 const int& IntegrationStep)
{
	Init();
	Set(t1_corr,t2_corr,start_t1_corr,start_t2_corr,nbnames,pdef,pdef_t2,pdef_start,
		pdef_start_t2,beta_t1,beta_t2,start_beta_t1,start_beta_t2,LossRates,discretizationstep ,CopulaType ,
		IntegrationMethod,IntegrationStep);

}

//-----------------------------------------------------------------------
// Forward collateral : Term Structure Review
//-----------------------------------------------------------------------
void ICM_Gauss1FLossDistrib_LJ::Set(const double& t1_corr,
													 const double& t2_corr,
													 const double& start_t1_corr,
													 const double& start_t2_corr,
													 const int& nbnames,
													 const ARM_Vector&  pdef,
													 const ARM_Vector&  pdef_t2,
													 const ARM_Vector&  pdef_start,
													 const ARM_Vector&  pdef_start_t2,
													 const ARM_Vector&  beta_t1,
													 const ARM_Vector&  beta_t2,
													 const ARM_Vector&  start_beta_t1,
													 const ARM_Vector&  start_beta_t2,
													 const ARM_Vector & LossRates,
													 const int& discretizationstep ,
													 const int& CopulaType ,
													 const qIntegratorChoice&	IntegrationMethod ,
													 const int& IntegrationStep)
{
	SetTSR(true);
	SetT1(t1_corr);
	SetT2(t2_corr);
	SetCollatT1(start_t1_corr);
	SetCollatT2(start_t2_corr);

	SetFwdCollat_TSR(nbnames,pdef,pdef_t2,beta_t1,beta_t2,start_beta_t1,start_beta_t2,LossRates,
				discretizationstep ,CopulaType ,IntegrationMethod ,IntegrationStep);
	//MaJ des Pdef Start dans le vecteur PdefPerturb
	SetPdefPerturb(pdef_start);
	SetPdefPerturb_t2(pdef_start_t2);
	//Estimation des barrières correspondantes
	compute_barrier_perturb();
}
