#include "ARMKernel\glob\firsttoinc.h"
/* -----------------------------------------------------------------
   -----------------------------------------------------------------

	FILE:		CREDIT_MANAGER_LHP.CPP
	PROJECT:	CAIR
	
#include "ARMKernel\glob\firsttoinc.h"
	DESCRIPTION:	implementation based on Large Homogeneous Portfolio

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

#include "ICMKernel/util/icm_gaussian.h"
#include "ICMKernel/util/icm_nig.h"


// ------------------------------------------------------------------------
//
//	Computes Probability that at time T we have 
//		LossMIN < Loss < LossMAX
//
//	Algorithm is  LARGE PORTFOLIO
// 
// ------------------------------------------------------------------------

void	CreditManager::ComputeLossesBeforeT_LHP(double T, CreditLossComputation ProbOrExpectationFlag, double LossMin, double LossMax, double& Result)
{
	double	TheValue;
	double	Recovery, DefProb;
	double	K_1, K_2;
	double	ND2_LossMin, ND2_LossMax;
	double	InvNorm_LossMin, InvNorm_LossMax;
	double	ExpectedLoss;

	// ----------------------------------
	if ((T <= 0.) || (CHECK_EQUAL(LossMin, LossMax)))
	{
		Result = 0.0;
		return;
	}
	// ----------------------------------

	double	MinusSQRTOneMinusCorrel;

	// FLAT CORRELATION
	MinusSQRTOneMinusCorrel	=	- sqrt(1.0 - its_CorrelationValue);

	// EXCEL CONSTRUCTION --> the last one of the Portfolio
//	i = Get_NbCredits() - 1;

	// RECOVERY
	Recovery	=	its_LHP_Recovery;
	K_1	=	LossMin	/ TheBasketNotional / (1.0 - Recovery);
	K_2	=	LossMax	/ TheBasketNotional / (1.0 - Recovery);

	// ----------------------------------
	// According to Copula Type
	// ----------------------------------
	switch (its_CopulaType)
	{
	case CCT_GAUSSIAN:
	
		// Computes N-1 (SurvivalProb(T))	(Inverse Cumulative Gaussian Distribution)
		
		// Curve i
		// Default Proba expects only yearterms
		DefProb = its_DefCurve->DefaultProba(T);	// T is already in year fraction

		// should test if DefProb = 0.0 or = 1.0
		if (CHECK_EQUAL(DefProb, 0.0))
			TheValue	=	_MINUS_INFINITY_;	//	minus infinity
		else if (CHECK_EQUAL(DefProb, 1.0))
			TheValue	=	_PLUS_INFINITY_;	//	plus infinity
		else
		{
			TheValue	=	its_CopulaChoice->Inverse_Cumulative_Density_Function(DefProb);
			TheValue	=	NAG_deviates_normal( DefProb);
		}

		// for Loss Min
		if (CHECK_EQUAL(K_1, 0.0))
			InvNorm_LossMin	=	_MINUS_INFINITY_;	//	minus infinity
		else if (CHECK_EQUAL(K_1, 1.0))
			InvNorm_LossMin	=	_PLUS_INFINITY_;	//	plus infinity
		else
			InvNorm_LossMin	=	NAG_deviates_normal( K_1);

		// for Loss Max
		if (CHECK_EQUAL(K_2, 0.0))
			InvNorm_LossMax	=	_MINUS_INFINITY_;	//	minus infinity
		else if (CHECK_EQUAL(K_2, 1.0))
			InvNorm_LossMax	=	_PLUS_INFINITY_;	//	plus infinity
		else
			InvNorm_LossMax	=	NAG_deviates_normal( K_2);

		// The Barrier = TheValue
		ND2_LossMin	=	ND2(-InvNorm_LossMin, TheValue, MinusSQRTOneMinusCorrel);
		ND2_LossMax	=	ND2(-InvNorm_LossMax, TheValue, MinusSQRTOneMinusCorrel);

		// EXPECTED LOSS
		ExpectedLoss	=	(ND2_LossMin - ND2_LossMax) / (K_2 - K_1);

		if (ProbOrExpectationFlag == CLC_EXPECTATION)
		{
			ExpectedLoss	*=	(LossMax - LossMin);
		}
		else if (ProbOrExpectationFlag == CLC_PROBABILITY)
		{
			;
		}
		else
			ICMTHROW(ERR_INVALID_DATA,"Computes either PROB or EXPECTATION with LHP!");

		break;
	
	case CCT_STUDENT:

		ICMTHROW(ERR_INVALID_DATA,"Student Copula not yet implemented for 1F Model!");
		break;
	
	case CCT_NIG:
		ICMTHROW(ERR_INVALID_DATA,"NIG Copula not yet implemented for 1F Model!");
		break;
		

	default:
		ICMTHROW(ERR_INVALID_DATA,"Unknown Copula Type!");
		break;
	}

	Result	=	ExpectedLoss;
}


// ------------------------------------------------------------------------
//
//	Computes Probability that at time T we have 
//		LossMIN < Loss < LossMAX
//
//	Algorithm is  LARGE PORTFOLIO
// 
// ------------------------------------------------------------------------

void	CreditManager::ComputeLossesBeforeT_LHP_JPM(double T, CreditLossComputation ProbOrExpectationFlag, double LossMin, double LossMax, double& Result)
{
	int		i; 
	double	z, pp, zmin, zmax; 
	double	TheValue, DefProb;
	double	TheCorrelation;
	double	Loss1F;	
	double	tmpVal;

	DoubleVector	X;
	DoubleVector	ProbNtD;
	DoubleVector	ProbNtD_z;

	double	LossVMin, LossVMax, DLoss;

	// resize arrays
	its_Current_DefProb.clear();
	its_Current_DefProb.resize(Get_NbCredits());

	X.resize(Get_NbCredits());

	// 0 to Get_NbCredits() included

	// ----------------------------------
	if ((T <= 0.) || (CHECK_EQUAL(LossMin, LossMax)))
	{
		Result = 0.0;
		return;
	}
	// ----------------------------------

	Loss1F = 0.;
	
	double	DiffLoss	=	(LossMax-LossMin);

	//------------------------------------------------
	// If NIG Copula
	//------------------------------------------------
	
	double	alpha, betanig, mu, delta, rho ;
	alpha = Get_NIG_Alpha();
	betanig = Get_NIG_Beta();
	mu = - alpha * betanig/sqrt(alpha*alpha-betanig*betanig);
	delta = alpha ;
	rho = Get_NIG_Rho() ;

	ICM_Nig* NigCopula = new ICM_Nig(alpha, betanig, rho) ;

	// ----------------------------------
	// According to Copula Type
	// ----------------------------------
	switch (its_CopulaType)
	{
	case CCT_GAUSSIAN:
	
		// Computes N-1 (SurvivalProb(T))	(Inverse Cumulative Gaussian Distribution)
		// Default Proba expects only yearterms
		DefProb = its_DefCurve->DefaultProba(T);	// T is already in year fraction

		// should test if DefProb = 0.0 or = 1.0
		if (DefProb == 0.0)
			TheValue	=	_MINUS_INFINITY_;	//	minus infinity
		else if (DefProb == 1.0)
			TheValue	=	_PLUS_INFINITY_;	//	plus infinity
		else
			TheValue	=	NAG_deviates_normal( DefProb);
	

		// RAZ to 0.0
		fill(ProbNtD.begin(), ProbNtD.end(), 0.0);
		
		// integration limits (standard error)
		zmin	=	-6.0;
		zmax	=	6.0;
		double	integration_size;
		integration_size	=	(zmax - zmin);

		for (i=1 ;i<=its_NIntegration_1F; i++)
		{
			// global variable
//			its_index_z	=	i-1;
			z = integration_size * its_Xi[i] + zmin;	
			
			// Knowing Factor value z, Evaluate Default Prob for each name at time T

//			case CT_FLAT:
			TheCorrelation	=	its_CorrelationValue;

			pp = TheValue  - sqrt(TheCorrelation) * z;

			if (! CHECK_EQUAL(fabs(TheCorrelation), 1.0))
			{
				//	Normal Cumulative Function
				double	tmpValue	=	its_CopulaChoice->Cumulative_Density_Function(pp / sqrt(1.0 - TheCorrelation));
				DefProb = NAG_cumul_normal(pp / sqrt(1.0 - TheCorrelation));
			}
			else
				if (pp > 0.0)	// + infinite => N(+infinite)
					DefProb = 1.0;
				else			// - infinite => N(-infinite)
					DefProb = 0.0;

			if (ProbOrExpectationFlag == CLC_EXPECTATION)
			{
				tmpVal	=	DefProb * (1.0 - its_LHP_Recovery) * TheBasketNotional;
				LossVMin	=	FMIN(FMAX(tmpVal - LossMin, 0.0), DiffLoss);
				LossVMax	=	0.0;
			}
			else if (ProbOrExpectationFlag == CLC_PROBABILITY)
			{
				LossVMin	=	DefProb * (1.0 - its_LHP_Recovery);
				LossVMax	=	0.0;

			}

			// ----------------------------------------------------------------
			DLoss	=	LossVMin - LossVMax;

			// rescale by tranche size
//			DLoss	*=	its_BasketNotional;

			// ----------------------------------------------------------------
			// GaussLegendre Integration for each Nth to Default
			tmpVal = its_CopulaChoice->Density_Function(z) * its_Wi[i] * integration_size;
			// ----------------------------------------------------------------
			Loss1F += tmpVal * DLoss; 
			// ----------------------------------------------------------------

		}
	
	break;

	case CCT_NIG:
		
			
		// Computes N-1 (SurvivalProb(T))	(Inverse Cumulative NIG Distribution)
		// Default Proba expects only yearterms
		DefProb = its_DefCurve->DefaultProba(T);	// T is already in year fraction

		// should test if DefProb = 0.0 or = 1.0
		if (DefProb == 0.0)
			TheValue	=	_MINUS_INFINITY_;	//	minus infinity
		else if (DefProb == 1.0)
			TheValue	=	_PLUS_INFINITY_;	//	plus infinity
		else
			//TheValue	=	g01fac(Nag_LowerTail, DefProb, NAGERR_DEFAULT);
			TheValue	=	NigCopula->GetDefaultBarrier(DefProb);
		
		// RAZ to 0.0
		fill(ProbNtD.begin(), ProbNtD.end(), 0.0);
		
		// integration limits (standard error)
		zmin	=	-6.0;
		zmax	=	6.0;
		//double	integration_size;
		integration_size	=	(zmax - zmin);

		for (i=1 ;i<=its_NIntegration_1F; i++)
		{
			// global variable
			//			its_index_z	=	i-1;
			z = integration_size * its_Xi[i] + zmin;	
			
			// Knowing Factor value z, Evaluate Default Prob for each name at time T

			//			case CT_FLAT:
			
			pp = TheValue  - rho * z;

			if (! CHECK_EQUAL(fabs(rho), 1.0))
			{
				//	NIG Cumulative Function
				//DefProb = s15abc(pp / sqrt(1.0 - TheCorrelation));
				double fact = sqrt(1-rho*rho)/rho ;
				ICM_Nig* NigCopula2 = new ICM_Nig(fact*alpha, fact*betanig, fact* mu, fact*delta, rho);
				DefProb = NigCopula2->GetDistribution(pp / sqrt(1.0 - rho*rho));
				delete NigCopula2;
			}
			else
				if (pp > 0.0)	// + infinite => N(+infinite)
					DefProb = 1.0;
				else			// - infinite => N(-infinite)
					DefProb = 0.0;
				
			if (ProbOrExpectationFlag == CLC_EXPECTATION)
			{
				tmpVal	=	DefProb * (1.0 - its_LHP_Recovery) * TheBasketNotional;
				LossVMin	=	FMIN(FMAX(tmpVal - LossMin, 0.0), DiffLoss);
				LossVMax	=	0.0;
			}
			else 
				if (ProbOrExpectationFlag == CLC_PROBABILITY)
				{
					LossVMin	=	DefProb * (1.0 - its_LHP_Recovery);
					LossVMax	=	0.0;
				}
				
				// ----------------------------------------------------------------
				DLoss	=	LossVMin - LossVMax;

				// rescale by tranche size
				//			DLoss	*=	its_BasketNotional;

				// ----------------------------------------------------------------
				// GaussLegendre Integration for each Nth to Default
				//tmpVal = StandardGaussianDensity(z) * its_Wi[i] * integration_size;
				tmpVal = NigCopula->GetDensity(z) * its_Wi[i] * integration_size;
				// ----------------------------------------------------------------
				Loss1F += tmpVal * DLoss; 
			// ----------------------------------------------------------------
		}
		
		delete NigCopula ;
		
		break;

	case CCT_STUDENT:

		ICMTHROW(ERR_INVALID_DATA,"Student Copula not yet implemented for 1F Model!");
		break;

	default:
		ICMTHROW(ERR_INVALID_DATA,"Unknown Copula Type!");
		break;
	}

	Result	=	Loss1F;
}



// ------------------------------------------------------------------------
//
//	Computes Probability that at time T we have 
//		LossMIN < Loss < LossMAX
//
//	Algorithm is  LARGE PORTFOLIO
// 
// ------------------------------------------------------------------------

void	CreditManager::ComputeLossesBeforeT_LHP_PLUS(double T, CreditLossComputation ProbOrExpectationFlag, double LossMin, double LossMax, double& Result)
{
	int	i;

	double	TheValue;
	double	Recovery, DefProb;
	double	K_1, K_2;
//	double	InvNorm_LossMax;
	double	ExpectedLoss;

	// ----------------------------------
	if ((T <= 0.) || (CHECK_EQUAL(LossMin, LossMax)))
	{
		Result = 0.0;
		return;
	}
	// ----------------------------------

	// EXCEL CONSTRUCTION --> the last one of the Portfolio
//	i = Get_NbCredits() - 1;

	// RECOVERY
	Recovery	=	its_LHP_Recovery;
	K_1	=	LossMin	/ TheBasketNotional / (1.0 - Recovery);
	K_2	=	LossMax	/ TheBasketNotional / (1.0 - Recovery);

	// ----------------------------------
	// LEHMAN notations
	// ----------------------------------
	
	double	/*A, */C, C0;
	double	A1, B1, A2, B2;
	double	Beta;
	double	SQRTOneMinusBetaSquare;
	double	Beta0;
	double	SQRTOneMinusBeta0Square;
	double	BetaSquare;
	double	Beta0Square;
	
	double	TheValueA, TheValueB, TheValueC;
	double	term1, term2, term3, term4;
	double	rho;

	rho	=	0.5;

	Beta	=	sqrt(its_CorrelationValue);
	BetaSquare	=	its_CorrelationValue;
	SQRTOneMinusBetaSquare	=	sqrt(1.0 - BetaSquare);

	double	New_LossMin;
	double	New_LossMax;

	DoubleVector	Bs;
//	double	B;
	double	TheLossAmount;

	Bs.resize(Get_NbCredits());

	double *Sigma = new double[4];
	double *LimitA = new double[4];
	double *LimitB = new double[4];

	LimitA[0]	=	0.0;
	LimitB[0]	=	0.0;

	Sigma[0]	=	0.0;
	Sigma[3]	=	Beta;

	// ----------------------------------
	// According to Copula Type
	// ----------------------------------
	switch (its_CopulaType)
	{
	case CCT_GAUSSIAN:

		if (ProbOrExpectationFlag == CLC_EXPECTATION_DERIV_SPREAD)
		{
			// Curve for Homogeneous Portfolio
			// Default Proba expects only yearterms
			DefProb = its_DefCurve->DefaultProba(T);	// T is already in year fraction

			// should test if DefProb = 0.0 or = 1.0
			if (CHECK_EQUAL(DefProb, 0.0))
				TheValue	=	_MINUS_INFINITY_;	//	minus infinity
			else if (CHECK_EQUAL(DefProb, 1.0))
				TheValue	=	_PLUS_INFINITY_;	//	plus infinity
			else
				TheValue	=	its_CopulaChoice->Inverse_Cumulative_Density_Function(DefProb);

			C	=	TheValue;

			//	maybe it is not optimal, deal with one credit spread
			//	its_CurrentHedgesIndex

			i	=	its_CurrentHedgesIndex;

			// --------------------------------------------------------------------------------
			// A Value
			// --------------------------------------------------------------------------------
			if (CHECK_EQUAL(LossMin, 0.0))
				TheValue	=	_MINUS_INFINITY_;	//	minus infinity
			else
				if (CHECK_EQUAL(Recovery, 1.0))
					TheValue	=	_PLUS_INFINITY_;
				else
				{
					K_1	=	LossMin	/ TheBasketNotional / (1.0 - Recovery);
					TheValue	=	its_CopulaChoice->Inverse_Cumulative_Density_Function(K_1);
				}
				
			A1	=	(C - SQRTOneMinusBetaSquare * TheValue);
			A1	/=	Beta;

			// --------------------------------------------------------------------------------
			// B Value
			// --------------------------------------------------------------------------------
			// Default Proba expects only yearterms
			DefProb = (its_ArrayDefaultCrv[i])->DefaultProba(T);	// T is already in year fraction

			// should test if DefProb = 0.0 or = 1.0
			if (DefProb == 0.0)
				C0	=	_MINUS_INFINITY_;	//	minus infinity
			else if (DefProb == 1.0)
				C0	=	_PLUS_INFINITY_;	//	plus infinity
			else
				C0	=	NAG_deviates_normal( DefProb);
			
			switch (its_CorrelationType)
			{
				case CT_FLAT:
					Beta0	=	sqrt(its_CorrelationValue);
					Beta0Square	=	its_CorrelationValue;
					break;
				
				case CT_BETA:
					Beta0	=	its_Beta[i];
					Beta0Square	=	Beta0 * Beta0;
					break;
			}

			TheLossAmount	=	its_LossesAmount[i];
			
			SQRTOneMinusBeta0Square	=	sqrt(1.0 - Beta0Square);

			New_LossMin		=	LossMin - TheLossAmount;

			if (CHECK_EQUAL(New_LossMin, 0.0))
				TheValue	=	_MINUS_INFINITY_;	//	minus infinity
			else
				if (CHECK_EQUAL(Recovery, 1.0))
					if (New_LossMin > 0.0)
						TheValue	=	_PLUS_INFINITY_;
					else
						// should never happen ! inconsistent
						TheValue	=	_MINUS_INFINITY_;	//	minus infinity
				else
				{
					K_1	=	New_LossMin	/ TheBasketNotional / (1.0 - Recovery);
					TheValue	=	its_CopulaChoice->Inverse_Cumulative_Density_Function(K_1);
				}
				
			B1	=	(C - SQRTOneMinusBetaSquare * TheValue);
			B1	/=	Beta0;

			// --------------------------------------------------------------------------------
			// SOME VALUES
			TheValueA	=	(A1 - Beta0 * C0)/SQRTOneMinusBeta0Square;
			TheValueB	=	(B1 - Beta0 * C0)/SQRTOneMinusBeta0Square;
			TheValueC	=	(C - Beta * Beta0 * C0) / sqrt(1.0 - BetaSquare * Beta0Square);

			rho	=	Beta	*	(1.0 - Beta0Square) / SQRTOneMinusBetaSquare / SQRTOneMinusBeta0Square;

			// --------------------------------------------------------------------------------
			// --------------------------------------------------------------------------------
			// TERM 1
			// --------------------------------------------------------------------------------
			if (CHECK_EQUAL(LossMin, 0.0))
				term1	=	0.0;
			else
			{
				term1	=	its_CopulaChoice->Cumulative_Density_Function(TheValueA);
				term1	*=	LossMin;	
			}

			// --------------------------------------------------------------------------------
			// TERM 2
			// --------------------------------------------------------------------------------
			if (CHECK_EQUAL(New_LossMin, 0.0))
				term2	=	0.0;
			else
			{
				term2	=	its_CopulaChoice->Cumulative_Density_Function(TheValueB);
				term2	*=	New_LossMin;
			}

			// --------------------------------------------------------------------------------
			// TERM 3
			// --------------------------------------------------------------------------------
			if (CHECK_EQUAL(Recovery, 1.0))
				term3	=	0.0;
			else
			{
				term3	=	ND2(TheValueC, TheValueB, rho);
				term3	*=	TheBasketNotional * (1.0 - Recovery);
			}

			// --------------------------------------------------------------------------------
			// TERM 4
			// --------------------------------------------------------------------------------
			term4	=	ND2(TheValueC, TheValueA, rho);

			// --------------------------------------------------------------------------------
			// FINALLY
			// --------------------------------------------------------------------------------
			ExpectedLoss	=	term1 + term2 + term3 - term4;


			// --------------------------------------------------------------------------------
			// A Value: LOSS MAX
			// --------------------------------------------------------------------------------
			if (CHECK_EQUAL(LossMax, 0.0))
				TheValue	=	_MINUS_INFINITY_;	//	minus infinity
			else
				if (CHECK_EQUAL(Recovery, 1.0))
					TheValue	=	_PLUS_INFINITY_;
				else
				{
					K_1	=	LossMax	/ TheBasketNotional / (1.0 - Recovery);
					TheValue	=	its_CopulaChoice->Inverse_Cumulative_Density_Function(K_1);
				}
				
			A2	=	(C - SQRTOneMinusBetaSquare * TheValue);
			A2	/=	Beta;

			// --------------------------------------------------------------------------------
			// B Value
			// --------------------------------------------------------------------------------
			New_LossMax		=	LossMax - TheLossAmount;

			if (CHECK_EQUAL(New_LossMax, 0.0))
				TheValue	=	_MINUS_INFINITY_;	//	minus infinity
			else
				if (CHECK_EQUAL(Recovery, 1.0))
					if (New_LossMax > 0.0)
						TheValue	=	_PLUS_INFINITY_;
					else
						// should never happen ! inconsistent
						TheValue	=	_MINUS_INFINITY_;	//	minus infinity
				else
				{
					K_1	=	New_LossMax	/ TheBasketNotional / (1.0 - Recovery);
					TheValue	=	its_CopulaChoice->Inverse_Cumulative_Density_Function(K_1);
				}
				
			B2	=	(C - SQRTOneMinusBetaSquare * TheValue);
			B2	/=	Beta0;

			// --------------------------------------------------------------------------------
			// SOME VALUES
			TheValueA	=	(A2 - Beta0 * C0)/SQRTOneMinusBeta0Square;
			TheValueB	=	(B2 - Beta0 * C0)/SQRTOneMinusBeta0Square;

			// --------------------------------------------------------------------------------
			// --------------------------------------------------------------------------------
			// TERM 1
			// --------------------------------------------------------------------------------
			if (CHECK_EQUAL(LossMax, 0.0))
				term1	=	0.0;
			else
			{
				term1	=	its_CopulaChoice->Cumulative_Density_Function(TheValueA);
				term1	*=	LossMax;	
			}

			// --------------------------------------------------------------------------------
			// TERM 2
			// --------------------------------------------------------------------------------
			if (CHECK_EQUAL(New_LossMax, 0.0))
				term2	=	0.0;
			else
			{
				term2	=	its_CopulaChoice->Cumulative_Density_Function(TheValueB);
				term2	*=	New_LossMax;
			}

			// --------------------------------------------------------------------------------
			// TERM 3
			// --------------------------------------------------------------------------------
			if (CHECK_EQUAL(Recovery, 1.0))
				term3	=	0.0;
			else
			{
				term3	=	ND2(TheValueC, TheValueB, rho);
				term3	*=	TheBasketNotional * (1.0 - Recovery);
			}

			// --------------------------------------------------------------------------------
			// TERM 4
			// --------------------------------------------------------------------------------
			term4	=	ND2(TheValueC, TheValueA, rho);

			// --------------------------------------------------------------------------------
			// FINALLY: LOSS MAX
			// --------------------------------------------------------------------------------
			ExpectedLoss	+=	term1 + term2 + term3 - term4;

			// --------------------------------------------------------------------------------
			// FINALLY:
			// --------------------------------------------------------------------------------
			ExpectedLoss	*=	T * (1.0 - DefProb) / (1.0 - Recovery);

		}
		else if (ProbOrExpectationFlag == CLC_EXPECTATION)
		{
			// NOT USED FOR PRICING, ONLY FOR HEDGES ANALYSIS

/*			
			// Computes N-1 (SurvivalProb(T))	(Inverse Cumulative Gaussian Distribution)
			
			// Curve for Homogeneous Portfolio
			// Default Proba expects only yearterms
			DefProb = its_DefCurve->DefaultProba(T);	// T is already in year fraction

			// should test if DefProb = 0.0 or = 1.0
			if (CHECK_EQUAL(DefProb, 0.0))
				TheValue	=	_MINUS_INFINITY_;	//	minus infinity
			else if (CHECK_EQUAL(DefProb, 1.0))
				TheValue	=	_PLUS_INFINITY_;	//	plus infinity
			else
			{
				TheValue	=	its_CopulaChoice->Inverse_Cumulative_Density_Function(DefProb);
	//			TheValue	=	g01fac(Nag_LowerTail, DefProb, NAGERR_DEFAULT);
			}

			C	=	TheValue;
			LimitA[2]	=	C;
			LimitB[2]	=	C;

			// --------------------------------------------------------------------------------
			// A Value
			// --------------------------------------------------------------------------------
			if (CHECK_EQUAL(LossMin, 0.0))
				TheValue	=	_MINUS_INFINITY_;	//	minus infinity
			else
				if (CHECK_EQUAL(Recovery, 1.0))
					TheValue	=	_PLUS_INFINITY_;
				else
				{
					K_1	=	LossMin	/ TheBasketNotional / (1.0 - Recovery);
					TheValue	=	its_CopulaChoice->Inverse_Cumulative_Density_Function(K_1);
				}
				
			A	=	(C - SQRTOneMinusBetaSquare * TheValue);
			A	/=	Beta;

			// --------------------------------------------------------------------------------
			// Bs' Values
			// --------------------------------------------------------------------------------
			for (i=0; i<Get_NbCredits(); i++)
			{
				// Default Proba expects only yearterms
				DefProb = (its_ArrayDefaultCrv[i])->DefaultProba(T);	// T is already in year fraction

				// should test if DefProb = 0.0 or = 1.0
				if (DefProb == 0.0)
					C0	=	_MINUS_INFINITY_;	//	minus infinity
				else if (DefProb == 1.0)
					C0	=	_PLUS_INFINITY_;	//	plus infinity
				else
					C0	=	g01fac(Nag_LowerTail, DefProb, NAGERR_DEFAULT);
				
				switch (its_CorrelationType)
				{
					case CT_FLAT:
						Beta0	=	sqrt(its_CorrelationValue);
						break;
					
					case CT_BETA:
						Beta0	=	its_Beta[i];
						break;
				}

				TheLossAmount	=	its_LossesAmount[i];
				SQRTOneMinusBeta0Square	=	sqrt(1.0 - Beta0 * Beta0);

				New_LossMin		=	LossMin - TheLossAmount;

				if (CHECK_EQUAL(New_LossMin, 0.0))
					TheValue	=	_MINUS_INFINITY_;	//	minus infinity
				else
					if (CHECK_EQUAL(Recovery, 1.0))
						if (New_LossMin > 0.0)
							TheValue	=	_PLUS_INFINITY_;
						else
							// should never happen ! inconsistent
							TheValue	=	_MINUS_INFINITY_;	//	minus infinity
					else
					{
						K_1	=	New_LossMin	/ TheBasketNotional / (1.0 - Recovery);
						TheValue	=	its_CopulaChoice->Inverse_Cumulative_Density_Function(K_1);
					}
					
				B	=	(C - SQRTOneMinusBetaSquare * TheValue);
				B	/=	Beta0;

				Bs[i]	=	B;

				// -----------------------------------------------------------------
				// for Loss Min

				// Correlation Matrix
				Sigma[1]	=	Beta * Beta0;
				Sigma[2]	=	Beta0;

				LimitB[1]	=	C0;
				LimitB[3]	=	B;
				LimitA[1]	=	C0;
				LimitA[3]	=	A;

				// first term
				ExpectedLoss	=	LossMin * (ND2(C0, A, Beta0) - its_CopulaChoice->Cumulative_Density_Function(A));
				ExpectedLoss	-=	New_LossMin * ND2(C0, B, Beta0);
	//			ExpectedLoss	+=	(1.0 - Recovery) * TheBasketNotional * (ND2(C, A, Beta) + ND3(LimitB, Sigma) - ND3(LimitA, Sigma));

			}	// loop Credits
			// -----------------------------------------------------------------

			// for Loss Max
			if (CHECK_EQUAL(K_2, 0.0))
				InvNorm_LossMax	=	_MINUS_INFINITY_;	//	minus infinity
			else if (CHECK_EQUAL(K_2, 1.0))
				InvNorm_LossMax	=	_PLUS_INFINITY_;	//	plus infinity
			else
				InvNorm_LossMax	=	g01fac(Nag_LowerTail, K_2, NAGERR_DEFAULT);

			// The Barrier = TheValue
	//		ND2_LossMin	=	ND2(-InvNorm_LossMin, TheValue, MinusSQRTOneMinusCorrel);
	//		ND2_LossMax	=	ND2(-InvNorm_LossMax, TheValue, MinusSQRTOneMinusCorrel);

			// EXPECTED LOSS
	//		ExpectedLoss	=	(ND2_LossMin - ND2_LossMax) / (K_2 - K_1);

			ExpectedLoss	=	0.0;

			ExpectedLoss	*=	(LossMax - LossMin);
*/
			ICMTHROW(ERR_INVALID_DATA,"LHP Plus only available for Hedges Analysis!");
		}
		else if (ProbOrExpectationFlag == CLC_PROBABILITY)
		{
//			ExpectedLoss	=	0.0;
			ICMTHROW(ERR_INVALID_DATA,"LHP Plus only available for Hedges Analysis!");
		}
		else
			ICMTHROW(ERR_INVALID_DATA,"LHP Plus only available for Hedges Analysis!");

		break;
	
	case CCT_STUDENT:

		ICMTHROW(ERR_INVALID_DATA,"Student Copula not yet implemented for LHP Plus Model!");
		break;
	
	case CCT_NIG:
		ICMTHROW(ERR_INVALID_DATA,"NIG Copula not yet implemented for LHP Plus Model!");
		break;
		

	default:
		ICMTHROW(ERR_INVALID_DATA,"Unknown Copula Type!");
		break;
	}

	delete[] Sigma;
	delete[] LimitA;
	delete[] LimitB;

	Result	=	ExpectedLoss;
}
