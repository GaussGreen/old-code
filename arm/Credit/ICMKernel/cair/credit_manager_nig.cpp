/* -----------------------------------------------------------------
   -----------------------------------------------------------------

	FILE:		CREDIT_MANAGER_LHP.CPP
	PROJECT:	CAIR
	
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
	if ((T == 0.) || (CHECK_EQUAL(LossMin, LossMax)))
	{
		Result = 0.0;
		return;
	}
	// ----------------------------------

	double	MinusSQRTOneMinusCorrel;

	// FLAT CORRELATION
	MinusSQRTOneMinusCorrel	=	- sqrt(1.0 - its_CorrelationValue);

	// EXCEL CONSTRUCTION --> the last one of the Portfolio
//	i = its_NbCredits - 1;

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
			TheValue	=	NAG_deviates_normal(DefProb);

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
	its_Current_DefProb.resize(its_NbCredits);

	X.resize(its_NbCredits);

	// 0 to its_NbCredits included

	// ----------------------------------
	if ((T == 0.) || (CHECK_EQUAL(LossMin, LossMax)))
	{
		Result = 0.0;
		return;
	}
	// ----------------------------------

	Loss1F = 0.;
	
	double	DiffLoss	=	(LossMax-LossMin);

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
				//	Normal Cumulative Function
				DefProb = NAG_cumul_normal(pp / sqrt(1.0 - TheCorrelation));
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

	case CCT_STUDENT:

		ICMTHROW(ERR_INVALID_DATA,"Student Copula not yet implemented for 1F Model!");
		break;

	default:
		ICMTHROW(ERR_INVALID_DATA,"Unknown Copula Type!");
		break;
	}

	Result	=	Loss1F;
}
