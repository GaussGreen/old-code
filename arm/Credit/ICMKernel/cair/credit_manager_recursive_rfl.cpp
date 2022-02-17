/* -----------------------------------------------------------------
   -----------------------------------------------------------------

	FILE:		CREDIT_MANAGER_RECURSIVE_RFL.CPP
	PROJECT:	CAIR
	
	DESCRIPTION:	implementation based on Random Factor Loading Model

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

void	CreditManager::ComputeProbabilityOrExpectationLossesBeforeT_RECURSIVE_RFL(double T, CreditLossComputation ProbOrExpectationFlag, double LossMin, double LossMax, double& Result)
{

	// ---------------------------------
	// GO TO HEDGES
	if (IsActivateShift())
	{
		ICMTHROW(ERR_INVALID_DATA,"Fast Spread Hedges not done for Random Factor Loading Model!");
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
	its_Current_DefProb.resize(its_NbCredits);

	X.resize(its_NbCredits);

	// 0 to its_NbCredits included

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

	double	TheBeta;
	double	TheSQRTOneMinusBeta;
	double	TheDefProb;
	
	double	integration_size;

	// RAZ to 0.0
	fill(ProbNtD.begin(), ProbNtD.end(), 0.0);

	if (its_CopulaType != CCT_GAUSSIAN)
		ICMTHROW(ERR_INVALID_DATA,"RFL dos just take into account GAUSSIAN Copula!");

	if (its_CorrelationType != CT_FACTOR_LOADING_2)
		ICMTHROW(ERR_INVALID_DATA,"RFL dos just take into account GAUSSIAN Copula!");

	switch (its_Correlation_RFL_Type)
	{
	case CRFL_2F_CONSTANT:	// name by name... maybe too much
		
		// integration limits (standard error)
		zmin	=	-6.0;
		zmax	=	6.0;
		integration_size	=	(zmax - zmin);

		for (i=1 ;i<=its_NIntegration_1F; i++)
		{
			// global variable
			its_index_z	=	i-1;
			z = integration_size * its_Xi[i] + zmin;
			its_factor_value	=	z;

			// Knowing Factor value z, Evaluate Default Prob for each name at time T
			if (z < its_FL_Alpha[0])
			{				
				TheBeta	=	its_FL_Beta1[0];
				TheSQRTOneMinusBeta	=	its_FL_Beta1_Complement[0];
			}
			else
			{
				TheBeta	=	its_FL_Beta2[0];
				TheSQRTOneMinusBeta	=	its_FL_Beta2_Complement[0];
			}
			
			for (j=0; j<its_NbCredits; j++)
			{
				pp = X[j]  - TheBeta * z;
				
				if (! CHECK_EQUAL(fabs(TheBeta), 1.0))
					//	Normal Cumulative Function
					TheDefProb	=	its_CopulaChoice->Cumulative_Density_Function(pp / TheSQRTOneMinusBeta);
				else
					if (pp > 0.0)	// + infinite => N(+infinite)
						TheDefProb = 1.0;
					else			// - infinite => N(-infinite)
						TheDefProb = 0.0;			

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

	case CRFL_2F_INDIVIDUAL:	// name by name... maybe too much
		
		// integration limits (standard error)
		zmin	=	-6.0;
		zmax	=	6.0;
		integration_size	=	(zmax - zmin);

		for (i=1 ;i<=its_NIntegration_1F; i++)
		{
			// global variable
			its_index_z	=	i-1;
			z = integration_size * its_Xi[i] + zmin;
			its_factor_value	=	z;

			// Knowing Factor value z, Evaluate Default Prob for each name at time T
			for (j=0; j<its_NbCredits; j++)
			{
				if (z < its_FL_Alpha[j])
				{				
					TheBeta	=	its_FL_Beta1[j];
					TheSQRTOneMinusBeta	=	its_FL_Beta1_Complement[j];
				}
				else
				{
					TheBeta	=	its_FL_Beta2[j];
					TheSQRTOneMinusBeta	=	its_FL_Beta2_Complement[j];
				}
			
				pp = X[j]  - TheBeta * z;

				if (! CHECK_EQUAL(fabs(TheBeta), 1.0))
					//	Normal Cumulative Function
					TheDefProb	=	its_CopulaChoice->Cumulative_Density_Function(pp / TheSQRTOneMinusBeta);
				else
					if (pp > 0.0)	// + infinite => N(+infinite)
						TheDefProb = 1.0;
					else			// - infinite => N(-infinite)
						TheDefProb = 0.0;				
				
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

	case CRFL_TANH:	// name by name... maybe too much
		
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
			for (j=0; j<its_NbCredits; j++)
			{
				TheBeta	=	tanh(its_FL_TanH_alpha + its_FL_TanH_delta * tanh(its_FL_TanH_theta * z + its_FL_TanH_mu));
				TheSQRTOneMinusBeta	=	sqrt(1.0 - TheBeta * TheBeta);
			
				pp = X[j]  - TheBeta * z;

				if (! CHECK_EQUAL(fabs(TheBeta), 1.0))
					//	Normal Cumulative Function
					TheDefProb	=	its_CopulaChoice->Cumulative_Density_Function(pp / TheSQRTOneMinusBeta);
				else
					if (pp > 0.0)	// + infinite => N(+infinite)
						TheDefProb = 1.0;
					else			// - infinite => N(-infinite)
						TheDefProb = 0.0;				
				
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
	
	default:
		ICMTHROW(ERR_INVALID_DATA,"Unknown Copula Type for RFL Model!");
		break;
	}	

	Result	=	Loss1F;
}




// ------------------------------------------------------------------------
//
//	RFL ALLOCATIONS and INITIALISATIONS
// 
// ------------------------------------------------------------------------


void	CreditManager::Allocation_RFL_Correlation_Data()
{
	double	tmpValue;

	// ALLOCATIONS
	its_FL_Beta1_Complement.clear();
	its_FL_Beta2_Complement.clear();
	its_FL_Beta1_Complement.resize(its_NbCredits);
	its_FL_Beta2_Complement.resize(its_NbCredits);

	// FILL them
	int j;
	double	tmpCorrel;

	for (j=0; j<its_NbCredits; j++)
	{
		// BETA 1
		tmpValue	=	its_FL_Beta1[j];
		tmpCorrel	=	tmpValue * tmpValue;
		tmpValue	=	sqrt(1.0 - tmpCorrel);
		its_FL_Beta1_Complement[j]	=	tmpValue;

		// BETA 2
		tmpValue	=	its_FL_Beta2[j];
		tmpCorrel	=	tmpValue * tmpValue;
		tmpValue	=	sqrt(1.0 - tmpCorrel);
		its_FL_Beta2_Complement[j]	=	tmpValue;

	}

}