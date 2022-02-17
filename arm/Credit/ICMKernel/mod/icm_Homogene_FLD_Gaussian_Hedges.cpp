
#include "ARMKernel\glob\firsttoinc.h" 
#include "ICMKernel\mod\icm_Homogene_FastLossDistrib_Gaussian.h"
#include "ICMKernel\glob\icm_maths.h"
#include "ICMKernel\util\icm_HermiteIntegration.h"

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

 

/**
17783 
extern "C" void 
compute_cond_distribZeroLoss_shift_k_Gaussian_Integrator(void* Param, double x, double& res)
{
	// get parameters from stack
	ICM_Gauss1FLossDistrib_LJ*	TheModel;
	TheModel	=	(ICM_Gauss1FLossDistrib_LJ*)((*(AddressVector*)Param)[0]);

	res	=	TheModel->compute_cond_distribZeroLoss_shift_k_Gaussian(x);
}

extern "C" void 
compute_cond_distrib_Gaussian_shift_k_Integrator(void* Param, double x, double& res)
{
	// get parameters from stack
	ICM_Gauss1FLossDistrib_LJ*	TheModel;
	TheModel	=	(ICM_Gauss1FLossDistrib_LJ*)((*(AddressVector*)Param)[0]);

	res	=	TheModel->compute_cond_distrib_shift_k_Gaussian(x);
}

void
ICM_Gauss1FLossDistrib_LJ::compute_distrib_perturb(const int& lup)
{
	switch (itsIntegrationMethod)
	{
	case qGAUSS_LEGENDRE:
	case qGAUSS_HERMITE:
	case qTRAPEZE:
		compute_distrib_perturb_Gauss_Legendre(lup);
		break;
	default:
		compute_distrib_perturb_Hermite(lup);
		break;
	}
}


void
ICM_Gauss1FLossDistrib_LJ::compute_distrib_perturb_Hermite(const int& lup)
{
	double tmpPrbLoss=0.,cumul_distrib=0.;
	int min_lup =0,max_lup =0, max_ldown =0;
	int k, l = 0;

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

	switch (itsCopulaType)
	{
		case qNO_COPULA:	
		case qGAUSSIAN:

		for (k=0;k<its_nbnames;k++) 
		{
			its_ind_x=-1;
			its_ind_name=k;
			
			tmpPrbLoss = HermiteIntegration(ff1::mem_call(&ICM_Gauss1FLossDistrib_LJ::compute_cond_distribZeroLoss_shift_k_Gaussian, (*this))).Integrate_20();	
			tmpPrbLoss*=its_normPi;

			its_lossdistrib_perturb->SetValue(0,k,tmpPrbLoss);
			cumul_distrib=tmpPrbLoss;

			its_ind_Loss=1;

			for (l=its_lup+1;l<=max_lup;l++) 
			{
				 its_ind_x=-1;

				 tmpPrbLoss = HermiteIntegration(ff1::mem_call(&ICM_Gauss1FLossDistrib_LJ::compute_cond_distrib_shift_k_Gaussian, (*this))).Integrate_20();
				 tmpPrbLoss*=its_normPi;
				 cumul_distrib+=tmpPrbLoss;
	 			 its_lossdistrib_perturb->SetValue(l,k,tmpPrbLoss);
				 its_ind_Loss+=1;
			}

			its_taildistrib_perturb[k]=1.-cumul_distrib;
		}
	
		break;
	}
}


// -------------------------------------------------------------------
// GAUSSIAN COPULA
// -------------------------------------------------------------------

//	Loss Level is 0.0


void
ICM_Gauss1FLossDistrib_LJ::compute_distrib_perturb_Gauss_Legendre(const int& lup)
{
	double tmpPrbLoss=0.,cumul_distrib=0.;
	int min_lup =0,max_lup =0;
	int k;

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

	int nointegrationcoef=-1;
	TheParameterVector.Append(&nointegrationcoef);
	// ------------------------

	switch (itsCopulaType)
	{
		case qNO_COPULA:	
		case qGAUSSIAN:

			for (k=0;k<its_nbnames;k++) 
			{
				its_ind_x=-1;
				its_ind_name=k;
				TheIntegrator.Integrate(-6.0, 6.0, compute_cond_distribZeroLoss_shift_k_Gaussian_Integrator, &TheParameterVector, tmpPrbLoss); 

				its_lossdistrib_perturb->SetValue(0,k,tmpPrbLoss);
				cumul_distrib=tmpPrbLoss;

				its_ind_Loss=1;

				for (int l=its_lup+1;l<=max_lup;l++) 
				{
					its_ind_x=-1;
					TheIntegrator.Integrate(-6.0, 6.0, compute_cond_distrib_Gaussian_shift_k_Integrator, &TheParameterVector, tmpPrbLoss); 
					cumul_distrib+=tmpPrbLoss;

	 				its_lossdistrib_perturb->SetValue(l,k,tmpPrbLoss);
					its_ind_Loss+=1;
				}

				its_taildistrib_perturb[k]=1.-cumul_distrib;
			}
		
			break;

		case qSTUDENT:
			break;
	}

}
**/ 
/**
// 17783 
void
ICM_Gauss1FLossDistrib_LJ::compute_expectedlosstranche_fast_spread_hedge(const double& tranche_up, 
																		const double& tranche_down, 
																		const double& lossunit,
																		ICM_QMatrix<double>* ShiftMatrix)
{
	int lup=floor(tranche_up/lossunit);
	int ldown=floor(tranche_down/lossunit);
	
	its_lup	=	0;
	
	its_ProbCond->ResizeWithCopy(lup+1,0.);
	its_ProbCond_Down->ResizeWithCopy(ldown+1,0.);
	its_lossdistrib.resize(lup+1);
	its_lossdistrib_Down.resize(ldown+1);

	// Compute the Distribution (by default, ldown = 0)
	if (tranche_down)
		compute_distrib(lup, ldown);
	else
		compute_distrib(lup);

	compute_expectedlosstranche_perturb(tranche_up,tranche_down,lossunit,ShiftMatrix);
}


void
ICM_Gauss1FLossDistrib_LJ::compute_expectedlosstranche_perturb(const double& tranche_up, 
														  const double& tranche_down, 
														  const double& lossunit,
														  ICM_QMatrix<double>* ShiftMatrix)
{

	int nbscenario=ShiftMatrix->Getnbrows();		// inversion
	ICM_QMatrix<double>* result= new ICM_QMatrix<double>(its_nbnames,nbscenario);

	int lup=floor(tranche_up/lossunit);
	int ldown=floor(tranche_down/lossunit);
	int s,k,l;	

	double tmp=0.;
	double	exp_loss_tranche_down, exp_loss_tranche_up;
	its_lup	=	0;

	for (s=0;s<nbscenario;s++)  //boucle sur les scenarios de hedges par buckets
	{
		SetPdefPerturb(*ShiftMatrix->RowAsVector(s));
		compute_barrier_perturb();
		precompute_coeffs_perturb();

		its_ProbCond_Perturb->ResizeWithCopy(lup+1,0.);
		its_lossdistrib_perturb->Resize(lup+1,its_nbnames);
		its_ProbCond_Perturb_Down->ResizeWithCopy(ldown+1,0.);
		its_lossdistrib_perturb_Down->Resize(ldown+1,its_nbnames);

		if (tranche_down)
			compute_distrib_perturb(lup, ldown);
		else
			compute_distrib_perturb(lup);

		for (k=0;k<its_nbnames;k++)  //boucle sur les scenarios de hedges par labels
		{
			result->SetValue(k,s,0.);
			
			exp_loss_tranche_down	= 0.;
			exp_loss_tranche_up		= 0.;

			double	loss_level;

			// first loop, deal with beta_down and up
			for (l=1; l<=ldown; l++) 
			{
				loss_level	=	l*lossunit;
				exp_loss_tranche_down	+=	(its_lossdistrib_perturb_Down->Getvalue(l,k)) * loss_level;		
				exp_loss_tranche_up		+=	(its_lossdistrib_perturb->Getvalue(l,k)) * loss_level;
			}
				
			// second loop, deal only with beat_up
			for (l=ldown+1;l<=lup;l++)
			{
				loss_level	=	l*lossunit;
				exp_loss_tranche_up		+=	(its_lossdistrib_perturb->Getvalue(l,k)) * loss_level;
			}

			// add the tail
			if (tranche_down)
				exp_loss_tranche_down	+=	tranche_down * its_taildistrib_perturb_Down[k];
			exp_loss_tranche_up		+=	tranche_up * its_taildistrib_perturb[k];

			result->SetValue(k,s,exp_loss_tranche_up - exp_loss_tranche_down);
		}
	}

	// 17783  SetShifts(result);
}
**/ 

/** 
17783 
double
ICM_Gauss1FLossDistrib_LJ::compute_cond_distribZeroLoss_shift_k_Gaussian(double x)
{
	double	pk, pshiftk;
	double	tmp_barrier_perturb, tmp_barrier;
	double	CondProbaZeroLossshiftk = 1.;

	its_ind_x+=1;

	double	beta_value;
	double	SQRT_one_minus_beta_square;

	if (itsIntegrationMethod	==	qGAUSS_HERMITE)
		x	*=	its_sqrt2;

	beta_value	=	its_unique_beta[its_ind_name];
	SQRT_one_minus_beta_square	=	1.0	- beta_value * beta_value;
	SQRT_one_minus_beta_square	=	sqrt(SQRT_one_minus_beta_square);

	tmp_barrier_perturb=(its_barrier_perturb[its_ind_name]-beta_value*x)/SQRT_one_minus_beta_square;
	pshiftk=NAG_cumul_normal(tmp_barrier_perturb);

	tmp_barrier=(its_barrier[its_ind_name]-beta_value*x)/SQRT_one_minus_beta_square;
	pk=NAG_cumul_normal(tmp_barrier);

	CondProbaZeroLossshiftk=its_ProbCond->Elt(its_nbnames,its_ind_x,0) / (1.-pk);

	its_ProbCond_Perturb->SetElt(its_ind_name,its_ind_x,0,CondProbaZeroLossshiftk);

	return (CondProbaZeroLossshiftk*(1.-pshiftk));
}


//	Loss Level is 0.0

extern "C" void 
compute_cond_distribZeroLoss_Gaussian_shift_k_Smile(void* Param, double x, double& ProbCond_Up, double& ProbCond_Down)
{
	double	pk=0.;
	double	tmp_barrier=0.;
	double	tmp_barrier_perturb;
	double	pshiftk, pshiftk_down;

	double CondProbaZeroLossshiftk		=	1.;
	double CondProbaZeroLoss_Down_shiftk;

	// get parameters from stack
	ICM_Gauss1FLossDistrib_LJ*	TheModel;
	TheModel	=	(ICM_Gauss1FLossDistrib_LJ*)((*(AddressVector*)Param)[0]);

	int*	its_ind_x;
	its_ind_x	=	(int*)((*(AddressVector*)Param)[1]);
	(*its_ind_x) += 1;

	int*	its_ind_Loss;
	its_ind_Loss	=	(int*)((*(AddressVector*)Param)[2]);

	int*	its_ind_name;
	its_ind_name	=	(int*)((*(AddressVector*)Param)[3]);

	int itsIntegrationMethod = TheModel->GetIntegrationMethod();

	if (itsIntegrationMethod	==	qGAUSS_HERMITE)
		x	*=	SQRT2;

	int	its_nbnames	= TheModel->GetNbNames();

	tmp_barrier_perturb	=	TheModel->GetCoeff_Perturb_a(*its_ind_name) - TheModel->GetCoeff_b(*its_ind_name) * x;
	pshiftk	=	NAG_cumul_normal(tmp_barrier_perturb);

	tmp_barrier	=	TheModel->GetCoeff_a(*its_ind_name) - TheModel->GetCoeff_b(*its_ind_name) * x;
	pk		=	NAG_cumul_normal(tmp_barrier);

	CondProbaZeroLossshiftk = TheModel->GetProbCond_Elt(its_nbnames, *its_ind_x, 0) / (1.-pk);
	TheModel->SetProbCond_Perturb_Elt(*its_ind_name, *its_ind_x, 0, CondProbaZeroLossshiftk);

	// Down
	tmp_barrier_perturb	=	TheModel->GetCoeff_Perturb_a_down(*its_ind_name) - TheModel->GetCoeff_b_down(*its_ind_name) * x;
	pshiftk_down	=	NAG_cumul_normal(tmp_barrier_perturb);

	tmp_barrier	=	TheModel->GetCoeff_a_down(*its_ind_name) - TheModel->GetCoeff_b_down(*its_ind_name) * x;
	pk	=	NAG_cumul_normal(tmp_barrier);

	CondProbaZeroLoss_Down_shiftk = TheModel->GetProbCond_Down_Elt(its_nbnames, *its_ind_x, 0) / (1.-pk);
	TheModel->SetProbCond_Perturb_Down_Elt(*its_ind_name, *its_ind_x, 0, CondProbaZeroLoss_Down_shiftk);

	ProbCond_Up		=	CondProbaZeroLossshiftk * (1.0 - pshiftk);
	ProbCond_Down	=	CondProbaZeroLoss_Down_shiftk * (1.0 - pshiftk_down);
}


// -------------------------------------------------------------------
//	Strictly Positive Loss Level
// -------------------------------------------------------------------


double
ICM_Gauss1FLossDistrib_LJ::compute_cond_distrib_shift_k_Gaussian(double x)
{

	its_ind_x+=1;

	double tmp1=0.; // proba de loss=L pour un panier de taille k
	double tmp2=0.,tmp3=0.; // proba de loss=L-1 pour un panier de taille k
	double pk=0.;     // proba cond de défaut du nom k
	int wk=0,ind_lastloss=0;

	double	tmp_barrier_perturb, pk_perturb;
	double	tmp;
	double	Prb_Cond_l_shiftk;
	double tmp_barrier=0.;
	double	tmp_prev;

	if (itsIntegrationMethod	==	qGAUSS_HERMITE)
		x	*=	its_sqrt2;

	if (its_ishomogeneous)
	{
		// barrier perturb
		tmp_barrier_perturb=(its_barrier_perturb[its_ind_name]-its_unique_beta[its_ind_name]*x)/sqrt(1.-its_unique_beta[its_ind_name]*its_unique_beta[its_ind_name]);
		// conditional perturb default probability
		pk_perturb=NAG_cumul_normal(tmp_barrier_perturb);

		// barrier
		tmp_barrier=(its_barrier[its_ind_name]-its_unique_beta[its_ind_name]*x)/sqrt(1.-its_unique_beta[its_ind_name]*its_unique_beta[its_ind_name]);
		// conditional default probability
		pk=NAG_cumul_normal(tmp_barrier);

		tmp_prev	=	its_ProbCond_Perturb->Elt(its_ind_name,its_ind_x,its_ind_Loss-1);

		tmp	= (its_ProbCond->Elt(its_nbnames,its_ind_x,its_ind_Loss)-tmp_prev*pk)/(1-pk);

		its_ProbCond_Perturb->SetElt(its_ind_name,its_ind_x,its_ind_Loss,tmp);

		Prb_Cond_l_shiftk = tmp*(1.-pk_perturb)+tmp_prev*pk_perturb;

	}
	else
	{
		// current level of loss - loss rate of names k (global variable k = its_ind_name)
		ind_lastloss = its_ind_Loss - its_int_lossrates[its_ind_name];

		// barrier perturb
		tmp_barrier_perturb=(its_barrier_perturb[its_ind_name]-its_unique_beta[its_ind_name]*x)/sqrt(1.-its_unique_beta[its_ind_name]*its_unique_beta[its_ind_name]);
		// conditional perturb default probability
		pk_perturb=NAG_cumul_normal(tmp_barrier_perturb);

		// barrier
		tmp_barrier=(its_barrier[its_ind_name]-its_unique_beta[its_ind_name]*x)/sqrt(1.-its_unique_beta[its_ind_name]*its_unique_beta[its_ind_name]);
		// conditional default probability
		pk=NAG_cumul_normal(tmp_barrier);

		tmp	= its_ProbCond->Elt(its_nbnames,its_ind_x,its_ind_Loss);
		if (ind_lastloss >= 0)
			tmp	-= its_ProbCond_Perturb->Elt(its_ind_name,its_ind_x,ind_lastloss)*pk;
		tmp	/= (1.-pk);

		its_ProbCond_Perturb->SetElt(its_ind_name,its_ind_x,its_ind_Loss,tmp);

		tmp2	=	tmp*(1.-pk_perturb);
		if (ind_lastloss >= 0)
			tmp2	+=	its_ProbCond_Perturb->Elt(its_ind_name,its_ind_x,ind_lastloss)*pk_perturb;

		Prb_Cond_l_shiftk	=	tmp2;
		
	}

	return Prb_Cond_l_shiftk;
}


extern "C" void 
compute_cond_distrib_Gaussian_shift_k_Smile(void* Param, double x, double& ProbCond_Up, double& ProbCond_Down)
{

	// get parameters from stack
	ICM_Gauss1FLossDistrib_LJ*	TheModel;
	TheModel	=	(ICM_Gauss1FLossDistrib_LJ*)((*(AddressVector*)Param)[0]);

	int*	its_ind_x;
	its_ind_x	=	(int*)((*(AddressVector*)Param)[1]);
	(*its_ind_x) += 1;

	int*	its_ind_Loss;
	its_ind_Loss	=	(int*)((*(AddressVector*)Param)[2]);

	int*	its_ind_name;
	its_ind_name	=	(int*)((*(AddressVector*)Param)[3]);

	double pk, pk_down;     // proba cond de défaut du nom k
	int ind_lastloss;

	double tmp, tmp_barrier, tmp_barrier_down, tmp_barrier_perturb_down;
	double	pk_perturb_down, tmp_down, tmp_barrier_perturb, pk_perturb;

	int itsIntegrationMethod = TheModel->GetIntegrationMethod();

	if (itsIntegrationMethod	==	qGAUSS_HERMITE)
		x	*=	SQRT2;

	int	its_nbnames	=	TheModel->GetNbNames();

	double	prev_ProbCond_Perturb;
	double	prev_ProbCond_Perturb_down;

	if (TheModel->IsHomogeneous())
	{
		// barriers perturb
		tmp_barrier_perturb	=	TheModel->GetCoeff_Perturb_a(*its_ind_name) - TheModel->GetCoeff_b(*its_ind_name) * x;
		tmp_barrier_perturb_down	=	TheModel->GetCoeff_Perturb_a_down(*its_ind_name) - TheModel->GetCoeff_b_down(*its_ind_name) * x;

		// conditional perturb default probability
		pk_perturb	=	NAG_cumul_normal(tmp_barrier_perturb);
		pk_perturb_down	=	NAG_cumul_normal(tmp_barrier_perturb_down);

		// barrier
		tmp_barrier = TheModel->GetCoeff_a(*its_ind_name) - TheModel->GetCoeff_b(*its_ind_name) * x;
		tmp_barrier_down = TheModel->GetCoeff_a_down(*its_ind_name) - TheModel->GetCoeff_b_down(*its_ind_name) * x;

		// conditional default probability
		pk	=	NAG_cumul_normal(tmp_barrier);
		pk_down	=	NAG_cumul_normal(tmp_barrier_down);

		tmp	=	(TheModel->GetProbCond_Elt(its_nbnames, *its_ind_x, *its_ind_Loss) - TheModel->GetProbCond_Perturb_Elt(*its_ind_name, *its_ind_x, *its_ind_Loss -1) * pk) / (1-pk);
		tmp_down	=	(TheModel->GetProbCond_Down_Elt(its_nbnames, *its_ind_x, *its_ind_Loss) - TheModel->GetProbCond_Perturb_Down_Elt(*its_ind_name, *its_ind_x, *its_ind_Loss -1) * pk_down) / (1-pk_down);

		TheModel->SetProbCond_Perturb_Elt(*its_ind_name, *its_ind_x, *its_ind_Loss, tmp);
		TheModel->SetProbCond_Perturb_Down_Elt(*its_ind_name, *its_ind_x, *its_ind_Loss, tmp_down);

		ProbCond_Up		=	tmp*(1.-pk_perturb) + TheModel->GetProbCond_Perturb_Elt(*its_ind_name, *its_ind_x, *its_ind_Loss-1)*pk_perturb;
		ProbCond_Down	=	tmp_down*(1.-pk_perturb_down) + TheModel->GetProbCond_Perturb_Down_Elt(*its_ind_name, *its_ind_x, *its_ind_Loss-1)*pk_perturb_down;

	}
	else
	{
		// current level of loss - loss rate of names k
		ind_lastloss = *its_ind_Loss - TheModel->GetIntLossRates()[*its_ind_name];

		// barriers perturb
		tmp_barrier_perturb	=	TheModel->GetCoeff_Perturb_a(*its_ind_name) - TheModel->GetCoeff_b(*its_ind_name) * x;
		tmp_barrier_perturb_down	=	TheModel->GetCoeff_Perturb_a_down(*its_ind_name) - TheModel->GetCoeff_b_down(*its_ind_name) * x;

		// conditional perturb default probability
		pk_perturb	=	NAG_cumul_normal(tmp_barrier_perturb);
		pk_perturb_down	=	NAG_cumul_normal(tmp_barrier_perturb_down);

		// barrier
		tmp_barrier = TheModel->GetCoeff_a(*its_ind_name) - TheModel->GetCoeff_b(*its_ind_name) * x;
		tmp_barrier_down = TheModel->GetCoeff_a_down(*its_ind_name) - TheModel->GetCoeff_b_down(*its_ind_name) * x;

		// conditional default probability
		pk	=	NAG_cumul_normal(tmp_barrier);
		pk_down	=	NAG_cumul_normal(tmp_barrier_down);

		// Loss Distribution for a basket of size k knowing factor value its_ind_x
		tmp	=	TheModel->GetProbCond_Elt(its_nbnames, *its_ind_x, *its_ind_Loss);
		tmp_down	=	TheModel->GetProbCond_Down_Elt(its_nbnames, *its_ind_x, *its_ind_Loss);
		if (ind_lastloss >= 0)
		{
			prev_ProbCond_Perturb	=	TheModel->GetProbCond_Perturb_Elt(*its_ind_name, *its_ind_x, ind_lastloss);
			prev_ProbCond_Perturb_down	=	TheModel->GetProbCond_Perturb_Down_Elt(*its_ind_name, *its_ind_x, ind_lastloss);

			tmp	-=	prev_ProbCond_Perturb * pk;
			tmp_down	-=	prev_ProbCond_Perturb_down * pk_down;
		}
		tmp	/=	(1. - pk);
		tmp_down	/=	(1. - pk_down);

		TheModel->SetProbCond_Perturb_Elt(*its_ind_name, *its_ind_x,*its_ind_Loss, tmp);
		TheModel->SetProbCond_Perturb_Down_Elt(*its_ind_name, *its_ind_x, *its_ind_Loss, tmp_down);

		ProbCond_Up		=	tmp*(1.-pk_perturb);
		ProbCond_Down	=	tmp_down*(1.-pk_perturb_down);
		if (ind_lastloss >= 0)
		{
			ProbCond_Up		+=	prev_ProbCond_Perturb * pk_perturb;
			ProbCond_Down	+=	prev_ProbCond_Perturb_down * pk_perturb_down;
		}

	}
}


void
ICM_Gauss1FLossDistrib_LJ::compute_distrib_perturb(const int& lup, const int& ldown)
{
	double tmpPrbLoss=0.,cumul_distrib=0.;
	double tmpPrbLoss_down=0.,cumul_distrib_down=0.;
	int min_lup =0,max_lup =0;
	int k, l;

	// ---------------------
	its_lup	=	0;
//	max_lup = MAX(lup,its_lup);
	
	// ----------------------
	// Parameters
	AddressVector	TheParameterVector;

	TheParameterVector.Append(this);
	TheParameterVector.Append(&its_ind_x);
	TheParameterVector.Append(&its_ind_Loss);
	TheParameterVector.Append(&its_ind_name);

	int nointegrationcoef=-1;
	TheParameterVector.Append(&nointegrationcoef);

	switch (itsCopulaType)
	{
		case qNO_COPULA:	
		case qGAUSSIAN:

			for (k=0;k<its_nbnames;k++) 
			{
				its_ind_x=-1;
				its_ind_name=k;

				max_lup	=	lup;
				
				TheIntegrator.Integrate(-6.0, 6.0, compute_cond_distribZeroLoss_Gaussian_shift_k_Smile, &TheParameterVector, tmpPrbLoss, tmpPrbLoss_down); 

				its_lossdistrib_perturb->SetValue(0,k,tmpPrbLoss);
				its_lossdistrib_perturb_Down->SetValue(0,k,tmpPrbLoss_down);

				cumul_distrib	=	tmpPrbLoss;
				cumul_distrib_down	=	tmpPrbLoss_down;

				its_ind_Loss=1;

				for (l=its_lup+1;l<=ldown;l++) 
				{
					its_ind_x=-1;

					TheIntegrator.Integrate(-6.0, 6.0, compute_cond_distrib_Gaussian_shift_k_Smile, &TheParameterVector, tmpPrbLoss, tmpPrbLoss_down); 

					its_lossdistrib_perturb->SetValue(l,k,tmpPrbLoss);
					its_lossdistrib_perturb_Down->SetValue(l,k,tmpPrbLoss_down);

					cumul_distrib +=	tmpPrbLoss;
					cumul_distrib_down +=	tmpPrbLoss_down;

					its_ind_Loss+=1;
				}
				for (l=ldown+1;l<=lup;l++) 
				{
					its_ind_x=-1;

					TheIntegrator.Integrate(-6.0, 6.0, compute_cond_distrib_Gaussian_shift_k_Integrator, &TheParameterVector, tmpPrbLoss); 

					its_lossdistrib_perturb->SetValue(l,k,tmpPrbLoss);
					cumul_distrib +=	tmpPrbLoss;

					its_ind_Loss+=1;
				}

				its_taildistrib_perturb[k]	=	1.0 - cumul_distrib;
				its_taildistrib_perturb_Down[k]	=	1.0 - cumul_distrib_down;
			}

			break;
	}
}

 17783 
 **/ 