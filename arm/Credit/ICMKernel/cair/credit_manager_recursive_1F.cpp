#include "ARMKernel\glob\firsttoinc.h"
/* -----------------------------------------------------------------
   -----------------------------------------------------------------

	FILE:		CREDIT_MANAGER_RECURSIVE_1F.CPP
	PROJECT:	CAIR
	
	DESCRIPTION:	implementation based on 1 Factor recursive algorithms

   -----------------------------------------------------------------
   
	CAIR Library

		version 1.0
		developped by Laurent Jacquel

  -----------------------------------------------------------------
  ----------------------------------------------------------------- */

#include "ICMKernel\cair\credit_manager.h"

#include "ICMKernel\glob\icm_maths.h"

#include <nags.h>		//	s15abc
#include <nagg01.h>		//	g01fac


// ------------------------------------------------------------------------
//
//	Computes Probability that at time T we have 
//		LossMIN < Loss < LossMAX
//
//	Algorithm is  FFT
// 
// ------------------------------------------------------------------------

void	CreditManager::ComputeProbabilityOrExpectationLossesBeforeT_RECURSIVE_1F(double T, CreditLossComputation ProbOrExpectationFlag, double LossMin, double LossMax, double& Result)
{
	// ---------------------------------
	// GO TO HEDGES
	if (IsActivateShift())
	{
		ComputeProbabilityOrExpectationLossesBeforeT_RECURSIVE_1F_Hedges_Spread(T, ProbOrExpectationFlag, LossMin, LossMax, Result);
		return;
	}
	// ---------------------------------

	int		i,j; 
	double	z, pp, zmin, zmax; 
	double	TheValue;
	double	TheCorrelation;
	double	Loss1F;	

	DoubleVector	X;
	DoubleVector	ProbNtD;
	DoubleVector	ProbNtD_z;

	double	LossVMin, LossVMax, DLoss;
	// resize arrays
	its_Current_DefProb.clear();
	its_Current_DefProb.resize(Get_NbCredits());

	X.resize(Get_NbCredits());

	// 0 to Get_NbCredits() included

	if (T <= 0.)
	{
		// should be improved in order to take into account instantaneous defaults?
		Result = 0.0;
		return;
	}

	Loss1F = 0.;

	// ----------------------------------
	// According to Copula Type
	// ----------------------------------
	
	// ----------------------------------
	// COMPUTE BARRIERS, Copula dependant
	ComputeBarriers(DAYSTIME(T));
	X	=	itsBarriers_Standard;
	// ----------------------------------

	// RAZ to 0.0
	fill(ProbNtD.begin(), ProbNtD.end(), 0.0);
	
	switch (its_CopulaType)
	{
	case CCT_GAUSSIAN:
		
		// integration limits (standard error)
		zmin	=	-6.0;
		zmax	=	6.0;
		double	integration_size;
		integration_size	=	(zmax - zmin);

		for (i=1 ;i<=its_NIntegration_1F; i++)
		{
			// global variable
			its_index_z	=	i-1;
			z = integration_size * its_Xi[i] + zmin;
			its_factor_value	=	z;
			
			// Knowing Factor value z, Evaluate Default Prob for each name at time T
			for (j=0; j<Get_NbCredits(); j++)
			{
				switch (its_CorrelationType)
				{
					case CT_FLAT:
						TheCorrelation	=	its_CorrelationValue;
						break;
					
					case CT_BETA:
						TheCorrelation	=	its_Beta[j] * its_Beta[j];
						break;
				}

				pp = X[j]  - sqrt(TheCorrelation) * z;

				if (! CHECK_EQUAL(fabs(TheCorrelation), 1.0))
				{
					//	Normal Cumulative Function
					double	tmpVale	=	its_CopulaChoice->Cumulative_Density_Function(pp / sqrt(1.0 - TheCorrelation));
					its_Current_DefProb[j] = NAG_cumul_normal(pp / sqrt(1.0 - TheCorrelation));
				}
				else
					if (pp > 0.0)	// + infinite => N(+infinite)
						its_Current_DefProb[j] = 1.0;
					else			// - infinite => N(-infinite)
						its_Current_DefProb[j] = 0.0;
			}

			if (ProbOrExpectationFlag == CLC_EXPECTATION)
			{
				LossVMin	=	Compute_Expected_LossTranche(LossMin, LossMax);
				LossVMax	=	0.0;

			}
			else if (ProbOrExpectationFlag == CLC_PROBABILITY)
			{
				LossVMin	=	Compute_Probability_LossTranche(LossMin);
				LossVMax	=	0.0;
			}
			else if (ProbOrExpectationFlag == CLC_EXPECTATION_INTERPOLATION)
			{
				LossVMin = LossCallInterp(LossMin);

				if (LossMax < its_TotalLoss) 
					LossVMax = LossCallInterp(LossMax);
				else
					LossVMax = 0.0;
			}
			else if (ProbOrExpectationFlag == CLC_PROBABILITY_INTERPOLATION)
			{
				LossVMin = LossDigitInterp(LossMin);

//				if (LossMax < its_TotalLoss) 
//					LossVMax = LossDigitInterp(LossMax);
//				else
					LossVMax = 0.0;
			}

			// ----------------------------------------------------------------
			DLoss = LossVMin - LossVMax;
			// ----------------------------------------------------------------

			// ----------------------------------------------------------------
			// GaussLegendre Integration for each Nth to Default
			TheValue = its_CopulaChoice->Density_Function(z) * its_Wi[i] * integration_size;
			// ----------------------------------------------------------------
			Loss1F += TheValue * DLoss; 
			// ----------------------------------------------------------------

		}
	
	break;

	case CCT_STUDENT:

		ICMTHROW(ERR_INVALID_DATA,"Student Copula not yet implemented for 1F Model!");
		break;

	default:
		ICMTHROW(ERR_INVALID_DATA,"Unknown Copula Type!");
		break;
	}

/*
	if (T == 30.0)
	{
		// -------------------------------------------
		// FOR A GIVEN T
		FILE *stream = fopen("c:\\test\\ProbCond.txt", "w+");
		// -------------------------------------------
		// VIEW
		fprintf(stream, "\tTime T:\t%lf\n", T);

		its_ProbCond->View("",stream);

		fclose(stream);
	}
	
*/
	Result	=	Loss1F;

}




// ------------------------------------------------------------------------
//
//	LOSS CALL its_With INTERPOLATION
//
//	Algorithm is RECURSIVE
// 
// ------------------------------------------------------------------------



double CreditManager::LossProbabilityDistribution(double Strike)
{
	return	LossCallInterp(Strike);
/*
	switch (its_Recursive_1F_Loss_Unit_Choice)
	{
	case COFRLUC_PGCD:
//return	LossProb();
		break;
		
	case COFRLUC_LOSS_UNIT_MIN:
		return	LossCallInterp_LossUnitMin(Strike);
		break;

	case COFRLUC_NB_LOSS_STEP:
		break;
	}

	return -1.0;
	*/
}


// --------------------------------------------------------
// this may be improved
// when the discretization matches exactly the loss amount
// --------------------------------------------------------

double CreditManager::LossCallInterp(double Strike)
{
	int		i, j, Index, Dim ;
	double	A, E, M, DX, lambda, Result;
	double	def_prob;

	double	CurrentLossAmount;

	if (Strike < 0.0)
		return	1.0;

	if (Strike >= its_TotalLoss)
		return	0.0;
	
	Dim = its_Recursive_NbLossSteps;

	M = 0.0;
	E = 0.0;
	DX = its_TotalLoss / ((double)Dim);
	Dim = ((int) (Strike / DX)) + 1;

	DoubleVector TMP1(Dim+1);
	DoubleVector TMP2(Dim+1);

	if (its_Algo_Type == 1)
	{
		// RAZ to 0.0
		fill(TMP1.begin(), TMP1.end(), 0.0);
		TMP1[Dim] = DX * ((double)Dim) - Strike;

		for (i=0; i < Get_NbCredits(); i++)
		{
			CurrentLossAmount	=	its_LossesAmount[i];
			A = CurrentLossAmount;

//			Index = (int) (CurrentLossAmount / DX);
			
//			lambda = CurrentLossAmount / DX - ((double)Index);
			
			// Current Def Prob
			def_prob	=	its_Current_DefProb[i];
			
			for (j=0; j<=Dim; j++)
			{
				if (A>Strike)
					TMP2[j] = (A + E - Strike) * def_prob + (1.0 - def_prob) * TMP1[j];
				else if (A<= Strike - M)
					TMP2[j] = 0.0 ;
				else
					TMP2[j] = TMP1[j + 1] * def_prob + (1.0 - def_prob) * TMP1[j];
//					TMP2[j] = (TMP1[j + Index] * (1.0 - lambda) + TMP1[j + Index + 1] * lambda) * def_prob + (1.0 - def_prob) * TMP1[j];
				
				A += DX;
			}

			M += CurrentLossAmount ;
			E += def_prob * CurrentLossAmount ;
			
			for (j=0; j<=Dim; j++)
				TMP1[j] = TMP2[j];
		}
	}
	else if (its_Algo_Type == 2)
	{

		// SECOND ALGORITHM
		int	jj;

		// RAZ to 0.0
		M = 0.0;
		E = 0.0;
		fill(TMP1.begin(), TMP1.end(), 0.0);
		fill(TMP2.begin(), TMP2.end(), 0.0);

		// last point, not null if interpolation is required
		TMP1[Dim] = DX * ((double)Dim) - Strike;

		for (i=0; i<Get_NbCredits(); i++)
		{
			CurrentLossAmount	=	its_LossesAmount[i];
			A = CurrentLossAmount;

//			Index = (int) (CurrentLossAmount / DX);
			
//			lambda = CurrentLossAmount / DX - ((double)Index);
			
			// Current Def Prob
			def_prob	=	its_Current_DefProb[i];
			
			jj	=	(int) ((Strike - M - CurrentLossAmount) / DX);
			jj	=	FMAX(jj, -1);
			A	+=	(jj+1) * DX;

			for (j=jj+1; j<=Dim; j++)
			{
				if (A>Strike)
					TMP2[j] = (A + E - Strike) * def_prob + (1.0 - def_prob) * TMP1[j];
				else if (A<= Strike - M)
					TMP2[j] = 0.0 ;
				else
					TMP2[j] = TMP1[j + 1] * def_prob + (1.0 - def_prob) * TMP1[j];
//					TMP2[j] = (TMP1[j + Index] * (1.0 - lambda) + TMP1[j + Index + 1] * lambda) * def_prob + (1.0 - def_prob) * TMP1[j];
				
				A += DX;
			}

			M += CurrentLossAmount ;
			E += def_prob * CurrentLossAmount ;
			
			for (j=0; j<=Dim; j++)
				TMP1[j] = TMP2[j];
		}
	}
	else if (its_Algo_Type == 3)
	{
		// RAZ to 0.0
		fill(TMP1.begin(), TMP1.end(), 0.0);
		TMP1[Dim] = DX * ((double)Dim) - Strike;

		for (i=0; i < Get_NbCredits(); i++)
		{
			CurrentLossAmount	=	its_LossesAmount[i];
			A = CurrentLossAmount;

			Index = (int) (CurrentLossAmount / DX);
			
			lambda = CurrentLossAmount / DX - ((double)Index);
			
			// Current Def Prob
			def_prob	=	its_Current_DefProb[i];
			
			for (j=0; j<=Dim; j++)
			{
				if (A>Strike)
					TMP2[j] = (A + E - Strike) * def_prob + (1.0 - def_prob) * TMP1[j];
				else if (A<= Strike - M)
					TMP2[j] = 0.0 ;
				else
					TMP2[j] = (TMP1[j + Index] * (1.0 - lambda) + TMP1[j + Index + 1] * lambda) * def_prob + (1.0 - def_prob) * TMP1[j];
				
				A += DX;
			}

			M += CurrentLossAmount ;
			E += def_prob * CurrentLossAmount ;
			
			for (j=0; j<=Dim; j++)
				TMP1[j] = TMP2[j];
		}
	}
	else if (its_Algo_Type == 4)
	{

		// SECOND ALGORITHM
		int	jj;

		// RAZ to 0.0
		M = 0.0;
		E = 0.0;
		fill(TMP1.begin(), TMP1.end(), 0.0);
		fill(TMP2.begin(), TMP2.end(), 0.0);

		// last point, not null if interpolation is required
		TMP1[Dim] = DX * ((double)Dim) - Strike;

		for (i=0; i<Get_NbCredits(); i++)
		{
			CurrentLossAmount	=	its_LossesAmount[i];
			A = CurrentLossAmount;

			Index = (int) (CurrentLossAmount / DX);
			
			lambda = CurrentLossAmount / DX - ((double)Index);
			
			// Current Def Prob
			def_prob	=	its_Current_DefProb[i];
			
			jj	=	(int) ((Strike - M - CurrentLossAmount) / DX);
			jj	=	FMAX(jj, -1);
			A	+=	(jj+1) * DX;

			for (j=jj+1; j<=Dim; j++)
			{
				if (A>Strike)
					TMP2[j] = (A + E - Strike) * def_prob + (1.0 - def_prob) * TMP1[j];
				else if (A<= Strike - M)
					TMP2[j] = 0.0 ;
				else
					TMP2[j] = (TMP1[j + Index] * (1.0 - lambda) + TMP1[j + Index + 1] * lambda) * def_prob + (1.0 - def_prob) * TMP1[j];
				
				A += DX;
			}

			M += CurrentLossAmount ;
			E += def_prob * CurrentLossAmount ;
			
			for (j=0; j<=Dim; j++)
				TMP1[j] = TMP2[j];
		}
	}

	Result = TMP1[0];

	return Result;
}


// --------------------------------------------------------
// this may be improved
// when the discretization matches exactly the loss amount
// --------------------------------------------------------

double CreditManager::LossDigitInterp(double Strike)
{
	int		i, j, Index, Dim;
	double	A, M, DX, lambda, Result;
	double	def_prob;

	double	CurrentLossAmount;

	if (Strike < 0.0)
		return	1.0;

	if (Strike >= its_TotalLoss)
		return	0.0;
	
	Dim = its_Recursive_NbLossSteps;

	M = 0.0;
	DX = its_TotalLoss / ((double)Dim);
	Dim = ((int) (Strike / DX)) + 1;

	DoubleVector TMP1(Dim+1);
	DoubleVector TMP2(Dim+1);

	// RAZ to 0.0
	fill(TMP1.begin(), TMP1.end(), 0.0);
	TMP1[Dim] = DX * ((double)Dim) - Strike;

	for (i=0; i < Get_NbCredits(); i++)
	{
		CurrentLossAmount	=	its_LossesAmount[i];
		A = CurrentLossAmount;

		Index = (int) (CurrentLossAmount / DX);
		
		lambda = CurrentLossAmount / DX - ((double)Index);
		
		// Current Def Prob
		def_prob	=	its_Current_DefProb[i];
		
		for (j=0; j<=Dim; j++)
		{
			if (A >= Strike)
				TMP2[j] = def_prob + (1.0 - def_prob) * TMP1[j];
			else if (A < Strike - M)
				TMP2[j] = 0.0 ;
			else
				TMP2[j] = (TMP1[j + Index] * (1.0 - lambda) + TMP1[j + Index + 1] * lambda) * def_prob + (1.0 - def_prob) * TMP1[j];
			
			A += DX;
		}

		M += CurrentLossAmount ;
		
		for (j=0; j<=Dim; j++)
			TMP1[j] = TMP2[j];
	}

	Result = TMP1[0];

	return Result;
}


double	CreditManager::Compute_Expected_LossTranche(double	Loss_Min, double Loss_Max)
{
	int	lup, ldown;
	double exp_loss_tranche_up;

	// -------------------------------------------------------
	// in order to keep some consistency its_With notations
	double	tranche_down;
	double	tranche_up;
	double	lossunit;

	tranche_down	=	Loss_Min;
	tranche_up		=	Loss_Max;
	
	lossunit		=	its_LossUnit;
	// -------------------------------------------------------

	ldown	=	floor(tranche_down/lossunit);
	lup		=	floor(tranche_up/lossunit);
	
	// Resize
	its_ProbCond->ResizeWithCopy(lup+1, 0.);
	its_lossdistrib.resize(lup+1);

	// Compute of the Loss Distribution for this given LossUp
	LossProbabilityDistribution(lup);

	int		l;
	double	loss_level;

	exp_loss_tranche_up	=	0.0;
	// first loop, deal its_With beta_down and up
/*	for (l=1; l<=ldown; l++) 
	{
		loss_level	=	l*lossunit;
		exp_loss_tranche_down	+=	its_lossdistrib_Down[l] * loss_level;		
		exp_loss_tranche_up		+=	its_lossdistrib[l] * loss_level;
	}
*/		
	// second loop, deal only its_With beat_up
	for (l=ldown+1; l<=lup; l++)
	{
		loss_level	=	l*lossunit;
		exp_loss_tranche_up		+=	its_lossdistrib[l] * (loss_level - tranche_down);
	}

	// add the tail
//	if (tranche_down)
//		exp_loss_tranche_down	+=	tranche_down * its_taildistrib_Down;
//	exp_loss_tranche_up		+=	tranche_up * its_taildistrib;

//	return (exp_loss_tranche_up - exp_loss_tranche_down);

	exp_loss_tranche_up		+=	(tranche_up - tranche_down) * its_taildistrib;
	
	return	exp_loss_tranche_up;
}


double	CreditManager::Compute_Probability_LossTranche(double Loss)
{
	int	lup;

	// -------------------------------------------------------
	// in order to keep some consistency its_With notations
	double	tranche_up;
	double	lossunit;

	// -------------------------------------------------------
	tranche_up	=	Loss;
	
	lossunit		=	its_LossUnit;
	
	// -------------------------------------------------------

	if (CHECK_EQUAL(Loss, TheBasketNotional))
		return 0.0;
	else if (CHECK_EQUAL(Loss, 0.0))
		return 1.0;
	else
	{
		// -------------------------------------------------------
		lup		=	floor(tranche_up/lossunit);
		// -------------------------------------------------------
		
		// -------------------------------------------------------
		// Resize
		its_ProbCond->ResizeWithCopy(lup+1, 0.);
		its_lossdistrib.resize(lup+1);
		// -------------------------------------------------------

		// -------------------------------------------------------
		// Compute of the Loss Distribution for this given LossUp
		LossProbabilityDistribution(lup);
		// -------------------------------------------------------
	}

	return	its_taildistrib;
}


void CreditManager::LossProbabilityDistribution(const int& lup)
{
	int		l;

	double	cumul_distrib;
	
	// global variable
	its_index_loss	=	0;
	its_lossdistrib[its_index_loss]	=	compute_cond_distribZeroLoss();
	
	for (its_index_loss=1; its_index_loss<=lup; its_index_loss++)
		its_lossdistrib[its_index_loss]	=	compute_cond_distrib();

	cumul_distrib	=	0.0;
	for (l=0; l<=lup; l++)
		cumul_distrib	+=	its_lossdistrib[l];

	its_taildistrib	=	1.0 - cumul_distrib;

}


double CreditManager::compute_cond_distribZeroLoss()
{
	double pk;

	int		k;
	double	CondProbaZeroLoss;
	
	CondProbaZeroLoss	=	1.0;

	// its_With no elements in the basket
	// its_index_loss = 0
	its_ProbCond->SetElt(0, its_index_z, 0, CondProbaZeroLoss);

	for (k=1; k<=Get_NbCredits(); k++) 
	{
		// Conditional default probability for Credit k
		pk	=	its_Current_DefProb[k-1];

		// loss level l = 0.0
		// no credit defaults!

		CondProbaZeroLoss	*=	(1-pk);

		// its_index_loss = 0
		its_ProbCond->SetElt(k, its_index_z, 0, CondProbaZeroLoss);
	}

	return CondProbaZeroLoss;
}


double CreditManager::compute_cond_distrib()
{
	int	k;
	
	// for the current level of loss
	double	loss_Basket_size_k; // proba de loss=L pour un panier de taille k
	
	double	tmp2; // proba de loss=L-1 pour un panier de taille k
	double	pk;     // proba cond de défaut du nom k
	int		ind_lastloss;

	// start its_With Loss = 0.0
	ind_lastloss	=	0;
	its_ProbCond->SetElt(0, its_index_z, its_index_loss, 0.);

	// its_index_Loss: Loss Distribution for a basket of size (0) Knowing factor value its_index_z
	loss_Basket_size_k	=	its_ProbCond->Elt(0, its_index_z, its_index_loss);

	if (IsHomogeneousBasketLossUnit())
	{
		for (k=1; k<=Get_NbCredits(); k++) 
		{
			// Conditional default probability for Credit k
			pk	=	its_Current_DefProb[k-1];

			// its_index_Loss: Loss Distribution for a basket of size (k-1) Knowing factor value its_index_z
			// tmp1

			// Loss Distribution for a basket of size k Knowing factor value its_index_z
			tmp2	=	loss_Basket_size_k * (1-pk) + pk * its_ProbCond->Elt(k-1, its_index_z, its_index_loss-1);

			its_ProbCond->SetElt(k, its_index_z, its_index_loss, tmp2);

			loss_Basket_size_k = tmp2;
		}
	}
	else
	{
		for (k=1; k<=Get_NbCredits(); k++) 
		{
			// current level of loss - loss rate of names k
			ind_lastloss = its_index_loss - its_LossRate[k-1];

			// Conditional default probability for Credit k
			pk	=	its_Current_DefProb[k-1];

			// its_index_Loss: Loss Distribution for a basket of size (k-1) Knowing factor value its_index_z
			// tmp1

			// Loss Distribution for a basket of size k Knowing factor value its_index_z
			tmp2	=	loss_Basket_size_k * (1-pk);
			
			if (ind_lastloss >= 0)
				tmp2	+=	pk * its_ProbCond->Elt(k-1, its_index_z, ind_lastloss);
			// else: the loss of name k is too big

			its_ProbCond->SetElt(k,its_index_z, its_index_loss, tmp2);

			loss_Basket_size_k = tmp2;
		}
	}

	return tmp2;
}


// --------------------------------------------------------------------------------------
//	HEDGES
// --------------------------------------------------------------------------------------

double	CreditManager::Compute_Expected_LossTranche_Fast_Spread_Hedge(const double& Loss_Min, const double& Loss_Max)
{
	int	lup, ldown;
	double	lossunit;
	double	tranche_down;
	double	tranche_up;

	tranche_down	=	Loss_Min;
	tranche_up		=	Loss_Max;
	
	lossunit		=	its_LossUnit;
	// -------------------------------------------------------

	ldown	=	floor(tranche_down/lossunit);
	lup		=	floor(tranche_up/lossunit);
	
//	its_ProbCond->ResizeWithCopy(lup+1,0.);
//	its_lossdistrib.resize(lup+1);

	// Compute the Distribution (by default, ldown = 0)
//	if (tranche_down)
//		LossProbabilityDistribution(lup, ldown);
//	else

	// DONE WITH PRICING
//	LossProbabilityDistribution(lup);

	return	Compute_Expected_LossTranche_Perturb(Loss_Min, Loss_Max);
}


double	CreditManager::Compute_Expected_LossTranche_Perturb(double	Loss_Min, double Loss_Max)
{
	int	lup, ldown;

	// -------------------------------------------------------
	// in order to keep some consistency its_With notations
	double	tranche_down;
	double	tranche_up;
	double	lossunit;
	double	tmp;

	tranche_down	=	Loss_Min;
	tranche_up		=	Loss_Max;
	
	lossunit		=	its_LossUnit;
	// -------------------------------------------------------

	ldown	=	floor(tranche_down/lossunit);
	lup		=	floor(tranche_up/lossunit);

	// -------------------------------------------------------
	int	k_id, s_id;

	s_id	=	GetTenorShift();
	k_id	=	GetIssuerShift();
	// -------------------------------------------------------
	
	int	l;
//	int	l, k, s;

//	for (s=0; s<its_nb_scenarii; s++)  //boucle sur les scenarios de hedges par buckets
//	{
		// ---------------------------------------------		
//		its_ProbCond_Perturb->ResizeWithCopy(lup+1, 0.0);
//		its_lossdistrib_perturb->ResizeWithCopy(lup+1, 0.0);

	// -------------------------------------------
		// Barrier Perturb
		Compute_Barrier_Perturb();
		
		// because  its_Current_DefProb has been computed (depends on s, k and factor integration)
		Compute_Distrib_Perturb(lup);

//		for (k=0; k<Get_NbCredits(); k++)  //boucle sur les scenarios de hedges par labels
//		{
//			its_AllMatrixNPVS->SetValue(k_id, s_id, 0.0);
			
			tmp	=	0.0;

			for (l=ldown+1; l<=lup; l++) 
				tmp	+=	(its_lossdistrib_perturb->Elt(k_id, s_id, l)) * (l * lossunit-tranche_down);

			tmp	+=	(tranche_up-tranche_down) * its_taildistrib_perturb->Getvalue(k_id, s_id);
			its_Shift_Results->SetValue(k_id, s_id, tmp);
//		}
//
//	}


			return	tmp;
//	SetShifts(result);

}


void CreditManager::Compute_Distrib_Perturb(const int& lup)
{
	int		l;
	double	tmpPrbLoss;
	double	cumul_distrib;
	
	int	min_lup;
	int	max_lup;

	double	its_lup;

	its_lup	=	0;
	min_lup	=	0;
	max_lup	=	lup;

	tmpPrbLoss		=	0.0;
	cumul_distrib	=	0.0;

	// Loop for each name
//	for (k=0; k<Get_NbCredits(); k++) 
//	{
	// -------------------------------------------------------
	int	k_id, s_id;

	s_id	=	GetTenorShift();
	k_id	=	GetIssuerShift();
	// -------------------------------------------------------

		its_ind_name		=	k_id;
		its_ind_scenario	=	s_id;
		
		max_lup	=	lup;
		
		its_index_loss	=	0;
		tmpPrbLoss		=	Compute_Cond_DistribZeroLoss_Shift_One_Name();
		
		its_lossdistrib_perturb->SetElt(k_id, s_id, 0, tmpPrbLoss);
		
		cumul_distrib	=	tmpPrbLoss;
		
		its_index_loss	=	1;

		for (l=its_lup+1; l<=max_lup; l++) 
		{
			 tmpPrbLoss = Compute_Cond_Distrib_Shift_One_Name();
			 
			 cumul_distrib	+=	tmpPrbLoss;
	 		 
			 its_lossdistrib_perturb->SetElt(k_id, s_id, l, tmpPrbLoss);
			 
			 its_index_loss	+=	1;
		}

		(*its_taildistrib_perturb)(k_id, s_id)	=	1.0 - cumul_distrib;
//	}
}



// for hedge computation
double CreditManager::Compute_Cond_DistribZeroLoss_Shift_One_Name()
{
	double	cond_def_prob;
	double	tmp_barrier_perturb;
	double	cond_def_prob_perturb;
	double	tmp_barrier;
	double	CondProbaZeroLossshiftk;

	// --------------------------------------------------
	// Correlation dealing

	double	beta_value;
	double	SQRT_one_minus_beta_square;

	beta_value	=	its_Used_Beta[its_ind_name];
	SQRT_one_minus_beta_square	=	its_Used_SQRT_OneMinusBetaSquare[its_ind_name];
	// --------------------------------------------------

	tmp_barrier_perturb	=	(its_barrier_perturb->Getvalue(its_ind_name, its_ind_scenario) - beta_value * its_factor_value);
	tmp_barrier_perturb	/=	SQRT_one_minus_beta_square;

	cond_def_prob_perturb	=	nag_cumul_normal(tmp_barrier_perturb);

	tmp_barrier			=	(itsBarriers_Standard[its_ind_name] - beta_value * its_factor_value);
	tmp_barrier			/=	SQRT_one_minus_beta_square;

	cond_def_prob		=	nag_cumul_normal(tmp_barrier);

	// here when I divide I should be more cautious...
	CondProbaZeroLossshiftk	=	its_ProbCond->Elt(Get_NbCredits(), its_index_z, 0) / (1.0 - cond_def_prob);

	its_ProbCond_Perturb->SetElt(its_ind_name, its_index_z, 0, CondProbaZeroLossshiftk);

	return (CondProbaZeroLossshiftk * (1.0 - cond_def_prob_perturb));
}



double CreditManager::Compute_Cond_Distrib_Shift_One_Name()
{
	double	cond_def_prob;
	double	cond_def_prob_perturb;

	double tmp;     // proba cond de défaut du nom k
	
	double	tmp_barrier;
	double	tmp_barrier_perturb;

	double Prb_Cond_l_k=0.;  // proba cond de loss=l pour un panier sans le name k
	double Prb_Cond_l_shiftk=0.;  // proba cond de loss=l pour un panier où le name k a été shifté
	double Prb_Cond_last_l_k=0.;  // proba cond de loss=l-1 pour un panier sans le name k
	int ind_lastloss=0;

	// --------------------------------------------------
	// Correlation dealing

	double	beta_value;
	double	SQRT_one_minus_beta_square;

	beta_value	=	its_Used_Beta[its_ind_name];
	SQRT_one_minus_beta_square	=	its_Used_SQRT_OneMinusBetaSquare[its_ind_name];
	// --------------------------------------------------

	if (IsHomogeneousBasketLossUnit())
	{
		tmp_barrier_perturb	=	(its_barrier_perturb->Getvalue(its_ind_name, its_ind_scenario) - beta_value * its_factor_value);
		tmp_barrier_perturb	/=	SQRT_one_minus_beta_square;

		cond_def_prob_perturb	=	nag_cumul_normal(tmp_barrier_perturb);

		tmp_barrier	=	(itsBarriers_Standard[its_ind_name] - beta_value * its_factor_value);
		tmp_barrier	/=	SQRT_one_minus_beta_square;

		cond_def_prob	=	nag_cumul_normal(tmp_barrier);

		tmp	=	its_ProbCond->Elt(Get_NbCredits(), its_index_z, its_index_loss);
		tmp	-=	its_ProbCond_Perturb->Elt(its_ind_name, its_index_z, its_index_loss-1) * cond_def_prob;
		tmp	/=	(1.0 - cond_def_prob);

		its_ProbCond_Perturb->SetElt(its_ind_name, its_NIntegration_1F, its_index_loss, tmp);

		Prb_Cond_l_shiftk	=	tmp * (1.0 - cond_def_prob_perturb);
		Prb_Cond_l_shiftk	+=	its_ProbCond_Perturb->Elt(its_ind_name, its_index_z, its_index_loss-1) * cond_def_prob_perturb;
	}
	else
	{
		ind_lastloss	=	its_index_loss - its_LossRate[its_ind_name];

		tmp_barrier_perturb		=	(its_barrier_perturb->Getvalue(its_ind_name, its_ind_scenario) - beta_value * its_factor_value);
		tmp_barrier_perturb		/=	SQRT_one_minus_beta_square;
		
		cond_def_prob_perturb	=	nag_cumul_normal(tmp_barrier_perturb);

		// pk = Conditional default probability for Credit k
		cond_def_prob	=	its_Current_DefProb[its_ind_name-1];

		if (ind_lastloss < 0)
		{
			tmp	=	its_ProbCond->Elt(Get_NbCredits(), its_index_z, its_index_loss);
			tmp	/=	(1.0 - cond_def_prob);

			its_ProbCond_Perturb->SetElt(its_ind_name, its_index_z, its_index_loss, tmp);
			
			Prb_Cond_l_shiftk	=	tmp * (1.0 - cond_def_prob_perturb);

		}
		else
		{
			tmp	=	its_ProbCond->Elt(Get_NbCredits(), its_index_z, its_index_loss);
			tmp	-=	its_ProbCond_Perturb->Elt(its_ind_name, its_index_z, ind_lastloss) * cond_def_prob;
			tmp	/=	(1.0 - cond_def_prob);

			its_ProbCond_Perturb->SetElt(its_ind_name, its_index_z, its_index_loss, tmp);

			Prb_Cond_l_shiftk	=	tmp * (1.0 - cond_def_prob_perturb);
			Prb_Cond_l_shiftk	+=	its_ProbCond_Perturb->Elt(its_ind_name, its_index_z, ind_lastloss) * cond_def_prob_perturb;
		}

	}

	return (Prb_Cond_l_shiftk);
}


// ---------------------------------------------
// Computes Barrier after a shift of Spreads
// ---------------------------------------------

void CreditManager::Compute_Barrier_Perturb()
{
	int k;
	double	TheValue;
	double Limit_case_Minus =	-10.;
	double Limit_case_Plus	=	10.;

//	its_barrier_perturb.resize(Get_NbCredits());	

	for (k=0; k<Get_NbCredits(); k++)
	{
		if (fabs(its_Current_DefProb[k]) < DB_TOL)
			TheValue	=	Limit_case_Minus;
		else if (fabs(its_Current_DefProb[k]-1.0) < DB_TOL)
			TheValue	=	Limit_case_Plus;
		else
			TheValue	=	NAG_deviates_normal_dist(its_Current_DefProb[k]);
		
		its_barrier_perturb->SetValue(k, its_ind_scenario, Limit_case_Plus);		
	}
}


// ---------------------------------------------
// FAST HEDGES
// ---------------------------------------------

// Spread Sensitivity

void CreditManager::HedgesFastSpreadParallel_Recursive(double& NPV)
{

	int	i, s;
	DoubleVector	Outputs;

	if (its_CurrentHedgesIndex != -1)
	{
		NPV	=	(*its_AllMatrixNPVS)(its_CurrentHedgesIndex, its_CurrentHedgesScenario);
		return;
	}

	// --------------------------------------------------------------
	// Activate Hedge Flag Computation
	ActivateShift();
	// --------------------------------------------------------------
		
	SetFirstComputationForHedgeFlag(true);

	// Clear vector
//	ClearShifts();

/*
FILE	*fOut;
if ((fOut = fopen("c:\\test\\hedge_times.txt", "a+")) == NULL) return;
DWORD	 t1, t2;
DWORD	 tb, te;

tb = GetTickCount();
*/

	// --------------------------------------------------------------
	// ALLOCATIONS
	// --------------------------------------------------------------
	if (its_AllMatrixNPVS)
		delete	its_AllMatrixNPVS;
	its_AllMatrixNPVS	=	new ICM_QMatrix<double>(Get_NbCredits(), its_nb_scenarii);
	
	if (its_Shift_Results)
		delete	its_Shift_Results;
	its_Shift_Results	=	new ICM_QMatrix<double>(Get_NbCredits(), its_nb_scenarii);
	
	if (its_ProbCond_Perturb)
		delete	its_ProbCond_Perturb;
	its_ProbCond_Perturb	=	new ICM_QCubix<double>(Get_NbCredits() + 1, its_NIntegration_1F, 1, 0.);	// some resize will be done

	if (its_lossdistrib_perturb)
		delete its_lossdistrib_perturb;
	its_lossdistrib_perturb = new ICM_QCubix<double>(Get_NbCredits(), its_nb_scenarii, 1, 0.0);	// some resize will be done l_up

	if (its_taildistrib_perturb)
		delete	its_taildistrib_perturb;
	its_taildistrib_perturb	=	new ICM_QMatrix<double>(Get_NbCredits(), its_nb_scenarii);

	if (its_barrier_perturb)
		delete	its_barrier_perturb;
	its_barrier_perturb	=	new ICM_QMatrix<double>(Get_NbCredits(), its_nb_scenarii);

	// --------------------------------------------------------------
	// COMPUTATION
	// --------------------------------------------------------------
	for (s=0; s<its_nb_scenarii; s++)
	{
		// Keep track of the current scenario
		SetTenorShift(s);	// Tenor
		its_ind_scenario	=	s; // GetTenorShift();

		for (i=0; i<Get_NbCredits(); i++)	//	The given Name
		{
			if (i)
				SetFirstComputationForHedgeFlag(false);

			SetIssuerShift(i);	// Label
			its_ind_name	=	i;

//t1 = GetTickCount();

			// -------------------------------------------------------
			// Reset Outputs
			// -------------------------------------------------------
			ResetOutputs();	// NPV to 0.0
				
			// Pricing
			PriceBasket_1F(Outputs);

			// Outputs
			(*its_AllMatrixNPVS)(i, s)		=	Outputs[0];
//			its_AllDefLegs[i]	=	Outputs[1];
//			its_AllPremLegs[i]	=	Outputs[2];
						

//t2 = GetTickCount();

//fprintf(fOut, "\nModified Price for label %u and tenor %u:\t%6.3f\n", i,j,(float) (t2-t1) / (float) 1000.0);

		}

		SetFirstComputationForHedgeFlag(false);
	}

	DeActivateShift();
			
/*
te = GetTickCount();
fprintf(fOut, "\nTOTAL TIME :\t%6.3f\n", (float) (te-tb) / (float) 1000.0);
fclose(fOut);
*/
}



// -------------------------------------------------------------------------------------------------
// COMPUTE HEDGES for SPREADS
// -------------------------------------------------------------------------------------------------

void	CreditManager::ComputeProbabilityOrExpectationLossesBeforeT_RECURSIVE_1F_Hedges_Spread(double T, CreditLossComputation ProbOrExpectationFlag, double LossMin, double LossMax, double& Result)
{
	// ---------------------------------
	// GO TO HEDGES
	if (! IsActivateShift())
		ICMTHROW(ERR_INVALID_DATA,"Compute Probability Or Expectation Losses Before T RECURSIVE 1F Hedges Spread!");
	// ---------------------------------

	// if First Passing, then I have to compute it all
	// otherwise, just pick up the result
	if (!IsFirstComputationForHedgeFlag())
	{
		// Already computed
		int	i_id, j_id;

		// and if 
		if (CHECK_EQUAL(T, 0.0))
			Result	=	0.0;
		else
		{
			j_id	=	GetTenorShift();
			i_id	=	GetIssuerShift();

			Result	=	GetShifts(T)->Getvalue(i_id, j_id);
		}

		return;
	}
	
	int		i,j; 
	double	z, pp, zmin, zmax; 
	double	TheValue, DefProb;
	double	Loss1F;	

	DoubleVector	X;
	DoubleVector	ProbNtD;
	DoubleVector	ProbNtD_z;

	double	LossVMin, LossVMax, DLoss;

	ICM_DefaultCurve*	TheShiftedDefaultCurve;
	ICM_DefaultCurve*	TheOriginalDefaultCurve;

	X.resize(Get_NbCredits());

	// 0 to Get_NbCredits() included

	if (T == 0.)
	{
		Result	=	0.0;
		// RAZ
//		its_Shift_Results->RAZ(0.0);
		return;
	}

	// -------------------------------------------------------
	int	lup, ldown;
	double	lossunit;
	double	tranche_down;
	double	tranche_up;

	tranche_down	=	LossMin;
	tranche_up		=	LossMax;
	
	lossunit		=	its_LossUnit;

	ldown	=	floor(tranche_down/lossunit);
	lup		=	floor(tranche_up/lossunit);

	// ---------------------------------------------		
	its_ProbCond_Perturb->ResizeWithCopy(lup+1, 0.0);
	its_lossdistrib_perturb->ResizeWithCopy(lup+1, 0.0);
	// -------------------------------------------------------

	ldown	=	floor(tranche_down/lossunit);
	lup		=	floor(tranche_up/lossunit);
	int	s, k;

	double	TheOriginalBarrier;
	double	TheBeta;
	double	TheSQRTOneMinusBetaSquare;

	for (s=0; s<its_nb_scenarii; s++)  //boucle sur les scenarios de hedges par buckets
	{
		for (k=0; k<Get_NbCredits(); k++)  // boucle sur les scenarios de hedges par labels
		{
			// -----------------------------------
			// CHANGE Prob DEF, just for this Credit and this current Scenario.
			TheShiftedDefaultCurve	=	(*its_MatrixShiftedDefaultCrv)(k, s);

			TheOriginalDefaultCurve	=	its_ArrayDefaultCrv[k];

			its_ArrayDefaultCrv[k]	=	TheShiftedDefaultCurve;

			TheOriginalBarrier		=	itsBarriers_Standard[k];
			// -----------------------------------

			// ----------------------------------
			// According to Copula Type
			// ----------------------------------
			switch (its_CopulaType)
			{
			case CCT_GAUSSIAN:
	
				// Computes The Shifted Barrier for the Current Credit along with the scenario
				
				// Default Proba expects only yearterms		// SHIFTED CURVES 
				DefProb = TheShiftedDefaultCurve->DefaultProba(T);	// T is already in year fraction

				// should test if DefProb = 0.0 or = 1.0
				if (DefProb == 0.0)
					TheValue	=	_MINUS_INFINITY_;	//	minus infinity
				else if (DefProb == 1.0)
					TheValue	=	_PLUS_INFINITY_;	//	plus infinity
				else
				{
					TheValue	=	its_CopulaChoice->Inverse_Cumulative_Density_Function(DefProb);
//					TheValue	=	g01fac(Nag_LowerTail, DefProb, NAGERR_DEFAULT);
				}
						
				itsBarriers_Standard[k]	=	TheValue;
	
				// RAZ to 0.0
				fill(ProbNtD.begin(), ProbNtD.end(), 0.0);
				
				// -----------------------------------
				// INTEGRATION
				// -----------------------------------
				
				// integration limits (standard error)
				zmin	=	-6.0;
				zmax	=	6.0;
				double	integration_size;
				integration_size	=	(zmax - zmin);

				X	=	itsBarriers_Standard;

				Loss1F	=	0.0;

				for (i=1 ;i<=its_NIntegration_1F; i++)
				{
					// global variable
					its_index_z	=	i-1;
					z = integration_size * its_Xi[i] + zmin;
					its_factor_value	=	z;
					
					// ----------------------------------
					// CONDITIONAL DEFAULT PROBABILITIES
					// ----------------------------------
					
					// Knowing Factor value z, Evaluate Default Prob for each name at time T
					for (j=0; j<Get_NbCredits(); j++)
					{
						TheBeta	=	its_Used_Beta[j];
						TheSQRTOneMinusBetaSquare	=	its_Used_SQRT_OneMinusBetaSquare[j];

						pp = X[j]  - TheBeta * z;

						if (! CHECK_EQUAL(fabs(TheBeta), 1.0))
						{
							//	Normal Cumulative Function
							double	tmpVale	=	its_CopulaChoice->Cumulative_Density_Function(pp / TheSQRTOneMinusBetaSquare);
							its_Current_DefProb[j] = NAG_cumul_normal(pp / TheSQRTOneMinusBetaSquare);
						}
						else
							if (pp > 0.0)	// + infinite => N(+infinite)
								its_Current_DefProb[j] = 1.0;
							else			// - infinite => N(-infinite)
								its_Current_DefProb[j] = 0.0;
					}
					// once, its_Current_DefProb has been computed, the Barrier can also

					if (ProbOrExpectationFlag == CLC_EXPECTATION)
					{
						Compute_Expected_LossTranche_Fast_Spread_Hedge(LossMin, LossMax);
						LossVMin	=	(*its_Shift_Results)(k, s);
						LossVMax	=	0.0;

					}				
		/*			else if (ProbOrExpectationFlag == CLC_PROBABILITY)
					{
						LossVMin	=	Compute_Probability_LossTranche(LossMin);
						LossVMax	=	0.0;
					}
					else if (ProbOrExpectationFlag == CLC_EXPECTATION_INTERPOLATION)
					{
						LossVMin = LossCallInterp(LossMin);

						if (LossMax < its_TotalLoss) 
							LossVMax = LossCallInterp(LossMax);
						else
							LossVMax = 0.0;
					}
					else if (ProbOrExpectationFlag == CLC_PROBABILITY_INTERPOLATION)
					{
						LossVMin = LossDigitInterp(LossMin);

		//				if (LossMax < its_TotalLoss) 
		//					LossVMax = LossDigitInterp(LossMax);
		//				else
							LossVMax = 0.0;
					}
*/
					
					// ----------------------------------------------------------------
					DLoss = LossVMin - LossVMax;
					// ----------------------------------------------------------------

					// ----------------------------------------------------------------
					// GaussLegendre Integration for each Nth to Default
					TheValue = its_CopulaChoice->Density_Function(z) * its_Wi[i] * integration_size;
					// ----------------------------------------------------------------
					Loss1F += TheValue * DLoss; 
					// ----------------------------------------------------------------
				}

				// STORE
				its_Shift_Results->SetValue(k, s, Loss1F);

				// -----------------------------------
				// RESTORE DEFAULT CURVE
				its_ArrayDefaultCrv[k]	=	TheOriginalDefaultCurve;

				itsBarriers_Standard[k]	=	TheOriginalBarrier;
				// -----------------------------------
	
				break;

			case CCT_STUDENT:

				ICMTHROW(ERR_INVALID_DATA,"Student Copula not yet implemented for 1F Model!");
				break;

			default:
				ICMTHROW(ERR_INVALID_DATA,"Unknown Copula Type!");
				break;
			}
		}
	}
/*
	if (T == 30.0)
	{
		FILE *stream = fopen("c:\\test\\ProbCond_Perturb.txt", "a+");
		// -------------------------------------------
		// VIEW
		fprintf(stream, "\tTime T:\t%lf\n", T);
	//	fprintf(stream,"Scenario:\t%u\t\tCredit:\t%u\n", s_id, k_id);

		its_ProbCond_Perturb->View("",stream);

		fclose(stream);
	}
*/	
	// -----------------------------------
	// push the ShiftMatrix in the vector (Clone???)
	AppendShifts(T, (ICM_QMatrix<double> *) its_Shift_Results->Clone());
	// -----------------------------------
	
	Result	=	(*its_Shift_Results)(0, 0);
}



// ----------------------------------------------------------------
// HEDGES: ZERO SHIFT for ALL NAMES
// ----------------------------------------------------------------
/*
double CreditManager::Compute_Cond_DistribZeroLoss_Shift_All_Names()
{
	double	cond_def_prob;
	double	tmp_barrier_perturb;
	double	cond_def_prob_perturb;
	double	tmp_barrier;
	double	CondProbaZeroLossshiftk;

	// --------------------------------------------------
	// Correlation dealing

	for (k=0; k<Get_NbCredits(); k++)
	{
		tmp_barrier			=	(itsBarriers_Standard[k] - beta_value * its_factor_value);
		tmp_barrier			/=	SQRT_one_minus_beta_square;

		cond_def_prob		=	nag_cumul_normal(tmp_barrier);

		for (s=0; s<its_nb_scenarii; s++)
		{
			tmp_barrier_perturb	=	(its_barrier_perturb->Getvalue(its_ind_name, its_ind_scenario) - beta_value * its_factor_value);
			tmp_barrier_perturb	/=	SQRT_one_minus_beta_square;

	double	beta_value;
	double	SQRT_one_minus_beta_square;

	beta_value	=	its_Used_Beta[its_ind_name];
	SQRT_one_minus_beta_square	=	its_Used_SQRT_OneMinusBetaSquare[its_ind_name];
	// --------------------------------------------------


	cond_def_prob_perturb	=	nag_cumul_normal(tmp_barrier_perturb);


	// here when I divide I should be more cautious...
	CondProbaZeroLossshiftk	=	its_ProbCond->Elt(Get_NbCredits(), its_index_z, 0) / (1.0 - cond_def_prob);

	its_ProbCond_Perturb->SetElt(its_ind_name, its_index_z, 0, CondProbaZeroLossshiftk);

	return (CondProbaZeroLossshiftk * (1.0 - cond_def_prob_perturb));
}

  */




void CreditManager::Hedges_Recursive(double& dummyNPV)
{
	int	i;

	ICM_DefaultCurve*	TheDefaultCurve	=	NULL;
	ICM_DefaultCurve*	TheShiftedDefaultCurve	=	NULL;

	int iscenario;

	// PARALLEL
	iscenario	=	0;

	// First passage
	if (its_CurrentHedgesIndex	==	-1)
	{
		for (i=0; i<Get_NbCredits(); i++)
		{
			// Just modify one curve, and reprice it all
			// I will have to improve it to take into account the Homogeneous case

			//	Retrieve Default Curve from Id
//			TheDefaultCurve	=	its_ModelMultiCurves->GetDefaultCurve(i);
			TheDefaultCurve	=	its_ArrayDefaultCrv[i];
			if (TheDefaultCurve == NULL)
				ICMTHROW(ERR_INVALID_DATA,"Parameters:  The Default Curve is NULL!");

			// Now deal with shifted Curve
//			TheShiftedDefaultCurve	=	(*its_MatrixShiftedDefaultCrv)(i, iscenario);
			TheShiftedDefaultCurve	=	(*its_MatrixShiftedDefaultCrv)(i, iscenario);
			if (TheShiftedDefaultCurve == NULL)
				ICMTHROW(ERR_INVALID_DATA,"Parameters:  The Shifted Default Curve is NULL!");

//			its_ModelMultiCurves->SetDefaultCurve(i, TheShiftedDefaultCurve);
			its_ArrayDefaultCrv[i]	=	TheShiftedDefaultCurve;

			// PRICE
			Price();

			// GET RESULTS
			its_AllNPVS[i]		=	NPV;
			its_AllDefLegs[i]	=	DefaultLegPV;
			its_AllPremLegs[i]	=	PremiumLegPV;

			// RESTORE
//			its_ModelMultiCurves->SetDefaultCurve(i, TheDefaultCurve);
			its_ArrayDefaultCrv[i]	=	TheDefaultCurve;

			if (IsHomogeneousBasketLossUnit())
			{
				for (int j=1; j<Get_NbCredits(); j++)
				{
					its_AllNPVS[j]		=	NPV;
					its_AllDefLegs[j]	=	DefaultLegPV;
					its_AllPremLegs[j]	=	PremiumLegPV;
				}

				return;
			}
		}
	}
	else
	{
		// simply get the right value
		dummyNPV		=	its_AllNPVS[its_CurrentHedgesIndex];

//		PremLeg	=	its_AllPremLegs[its_CurrentHedgesIndex];

//		DefLeg	=	its_AllDefLegs[its_CurrentHedgesIndex];
	}

}