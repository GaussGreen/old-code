#include "ARMKernel\glob\firsttoinc.h"
/* -----------------------------------------------------------------
   -----------------------------------------------------------------

	FILE:		CREDIT_MANAGER_RECURSIVE_STOCHASTIC_CORRELATION_BERNOULLI_1F.CPP
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
//	Algorithm is  RECURSION
// 
// ------------------------------------------------------------------------

void	CreditManager::ComputeProbabilityOrExpectationLossesBeforeT_RECURSIVE_STOCHASTIC_CORRELATION_BERNOULLI_1F(double T, CreditLossComputation ProbOrExpectationFlag, double LossMin, double LossMax, double& Result)
{

	// ---------------------------------
	// GO TO HEDGES
	if (IsActivateShift())
	{
		ICMTHROW(ERR_INVALID_DATA,"Fast Spread Hedges not done for Stochastic Correlation 1F Model!");

//		ComputeProbabilityOrExpectationLossesBeforeT_RECURSIVE_STOCHASTIC_CORRELATION_BERNOULLI_1F(T, ProbOrExpectationFlag, LossMin, LossMax, Result);
		return;
	}
	// ---------------------------------

	int		i,j; 
	double	z, pp, zmin, zmax; 
	double	TheValue;
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

	int	k;
	double	TheBeta;
	double	TheSQRTOneMinusBeta;
	double	TheDefProb;
	double	CurrDefProb;
	
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
				TheDefProb	=	0.0;

				for (k=0; k<its_SC_NbFactors; k++)
				{
					TheBeta				=	its_Used_Beta_SC_Vector[k];
					TheSQRTOneMinusBeta	=	its_Used_SQRT_OneMinusBetaSquare_SC_Vector[k];
				
					pp = X[j]  - TheBeta * z;

					if (! CHECK_EQUAL(fabs(TheBeta), 1.0))
						//	Normal Cumulative Function
						CurrDefProb	=	its_CopulaChoice->Cumulative_Density_Function(pp / TheSQRTOneMinusBeta);
					else
						if (pp > 0.0)	// + infinite => N(+infinite)
							CurrDefProb = 1.0;
						else			// - infinite => N(-infinite)
							CurrDefProb = 0.0;
					TheDefProb	+=	CurrDefProb * its_SC_Coefficients[j];
				}
				its_Current_DefProb[j]	=	TheDefProb;
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

	Result	=	Loss1F;
}




// ------------------------------------------------------------------------
//
//	STOCHASTIC CORRELATION ALLOCATIONS and INITIALISATIONS
// 
// ------------------------------------------------------------------------

void	CreditManager::Allocation_Stochastic_Correlation_Data()
{
	int		NbFactors;
	double	tmpValue;
	DoubleVector	tmpVector;
	DoubleVector	tmpVectorCoeff;

	switch (its_CreditModelType)
	{
	case CMT_ANALYTIC_RECURSIVE_STOCHASTIC_CORRELATION_BERNOULLI_1F:
		
		NbFactors	=	2;
		
		tmpVector.resize(NbFactors);
		
		tmpVector[0]	=	its_SC_Rho1;
		tmpVector[1]	=	its_SC_Rho2;

		tmpVectorCoeff.resize(NbFactors);

		tmpVectorCoeff[0]	=	its_SC_Prob;
		tmpVectorCoeff[1]	=	1.0 - its_SC_Prob;

		break;

	default:
		ICMTHROW(ERR_INVALID_DATA,"Allocation_Stochastic_Correlation_Data: Unknown Model Type!");
	}

	its_SC_NbFactors	=	NbFactors;

	// ALLOCATIONS
	its_Used_Beta_SC_Vector.resize(its_SC_NbFactors);
	its_Used_SQRT_OneMinusBetaSquare_SC_Vector.resize(its_SC_NbFactors);

	// FILL them
	int j;
	double	tmpCorrel;

	for (j=0; j<its_SC_NbFactors; j++)
	{
		// correl
		tmpCorrel	=	tmpVector[0];

		// beta
		tmpValue	=	sqrt(tmpCorrel);
		its_Used_Beta_SC_Vector[j]	=	tmpValue;

		// SQRT(1.0 - correl)
		tmpValue	=	sqrt(1.0 - tmpCorrel);
		its_Used_SQRT_OneMinusBetaSquare_SC_Vector[j]	=	tmpValue;

		its_SC_Coefficients[j]	=	tmpVectorCoeff[j];
	}

}