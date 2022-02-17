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

#include "ICMKernel\util\icm_integrator.h"
#include "ICMKernel\glob\icm_maths.h"

#include <nags.h>		//	s15abc
#include <nagg01.h>		//	g01fac
#include <numeric>


// ------------------------------------------------------------------------
//
//	Computes Probability that at time T we have 
//		LossMIN < Loss < LossMAX
//
//	Algorithm is  LARGE PORTFOLIO
// 
// ------------------------------------------------------------------------

void	CreditManager::ComputeLossesBeforeT_LARGE_PORTFOLIO_1F(double T, CreditLossComputation ProbOrExpectationFlag, double LossMin, double LossMax, double& Result)
{
	int		i,j; 
	double	z, pp, zmin, zmax; 
	double	TheValue, DefProb;
	double	TheCorrelation;
	double	Loss1F;	

	DoubleVector	X;
	DoubleVector	ProbNtD;
	DoubleVector	ProbNtD_z;

	double	LossVMin, LossVMax, DLoss;

	double	TheMean, TheVariance;

	double	TheLossAmount;
	double	TheLossFraction;
	double	TheCurrentMean;

	// ----------------------------------
	// Hermite Polynomials data
/*
	double	z_01;
	double	z_01_Square;
	double	z_01_Cube;
*/
	DoubleVector	Hi_z_01;

	// Up to _MAX_LP_Degree_ = 5
	its_CurrentHermiteCoeffs.resize(_MAX_LP_Degree_ + 1);
	Hi_z_01.resize(_MAX_LP_Degree_ + 1);

	double	HermiteSum;
	HermiteSum	=	0.0;
	// ----------------------------------
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


	TheMean		=	0.0;
	TheVariance	=	0.0;

	double	LossMinInPct;
	double	LossMaxInPct;

	LossMinInPct	=	LossMin / TheBasketNotional;
	LossMaxInPct	=	LossMax / TheBasketNotional;

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
			
			// Knowing Factor value z, Evaluate Default Prob for each name at time T
			TheMean		=	0.0;
			TheVariance	=	0.0;

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
					double	tmpValue	=	its_CopulaChoice->Cumulative_Density_Function(pp / sqrt(1.0 - TheCorrelation));
					its_Current_DefProb[j] = NAG_cumul_normal(pp / sqrt(1.0 - TheCorrelation));
				}
				else
					if (pp > 0.0)	// + infinite => N(+infinite)
						its_Current_DefProb[j] = 1.0;
					else			// - infinite => N(-infinite)
						its_Current_DefProb[j] = 0.0;

				// Mean and Variance are conditionally independent
				DefProb			=	its_Current_DefProb[j];
				TheLossAmount	=	its_LossesAmount[j];
				TheLossFraction	=	TheLossAmount / TheBasketNotional;

				TheCurrentMean	=	TheLossFraction * DefProb;
				
				TheMean			+=	TheCurrentMean;
				TheVariance		+=	TheCurrentMean * TheLossFraction * (1.0 - DefProb);
			}

			if (its_LP_Degree	==	0)
			{
				// STANDARD APPROXIMATION: NORMAL DENSITY
				if (ProbOrExpectationFlag == CLC_EXPECTATION)
				{
					LossVMin = LossLargePortfolioTranche(TheMean, TheVariance, LossMinInPct, LossMaxInPct);
					LossVMax = 0.0;
					
					// loop
				}
				else if (ProbOrExpectationFlag == CLC_PROBABILITY)
				{
					ICMTHROW(ERR_INVALID_DATA,"1F LARGE PORTFOLIO - Probability Computation not implemented yet!");
				}
			}
			else
			{
				// Check
				if (its_LP_Degree > _MAX_LP_Degree_)
					ICMTHROW(ERR_INVALID_DATA,"1F LARGE PORTFOLIO - Implemented only for degree N <= " << _MAX_LP_Degree_);

				// Computes its_CurrentHermiteCoeffs
				Compute_LargePortfolioHermitePolynomials();

				// STANDARD APPROXIMATION: NORMAL DENSITY
				if (ProbOrExpectationFlag == CLC_EXPECTATION)
				{
					LossVMin = LossLargePortfolioTranche_Hermite(TheMean, TheVariance, LossMinInPct, LossMaxInPct);
					LossVMax = 0.0;
					
					// loop
				}
				else if (ProbOrExpectationFlag == CLC_PROBABILITY)
				{
					ICMTHROW(ERR_INVALID_DATA,"1F LARGE PORTFOLIO - Probability Computation not implemented yet!");
				}


				// Pseudo-Loop
/*
				// Deal with centered and normalized zi
				z_01		=	(z - TheMean) / sqrt(TheVariance);
				z_01_Square	=	z_01 * z_01;
				z_01_Cube	=	z_01_Square * z_01;

				// k= 0
				// H0(z_01) = 1
				Hi_z_01[0]	=	1.0;

				// k= 1
				// H1(z_01) = z_01
				Hi_z_01[1]	=	z_01;

				// k= 2
				// H2(z_01) = z_01^2 - 1
				Hi_z_01[2]	=	z_01_Square - 1.0;

				// k= 3
				// H3(z_01) = z_01^3 - 3 * z_01
				Hi_z_01[3]	=	z_01_Cube - 3.0 * z_01;

				// k= 4
				// H4(z_01) = z_01^4 - 6 * z_01^2 + 3
				Hi_z_01[4]	=	z_01_Square * z_01_Square - 6.0 * z_01_Square + 3.0;

				// k= 5
				// H5(z_01) = z_01^5 - 10 * z_01^3 + 15 * z_01
				Hi_z_01[5]	=	z_01_Cube * z_01_Square - 10.0 * z_01_Cube + 15.0 * z_01;

				HermiteSum	=	inner_product(Hi_z_01.begin(), Hi_z_01.end(), its_CurrentHermiteCoeffs.begin(), 0.0);

				HermiteSum	*=	StandardGaussianDensity(z_01);
				HermiteSum	/=	sqrt(TheVariance);
*/
			}

			// ----------------------------------------------------------------
			DLoss	=	LossVMin - LossVMax;

			DLoss	*=	TheBasketNotional;

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


// --------------------------------------------------------
// Computation of TLoss(a,b)
// --------------------------------------------------------

// I avoid the 1.0 / (b-a) ratio
double CreditManager::LossLargePortfolioTranche(double mean, double variance, double LossMin, double LossMax)
{
	double	a, b;
	double	m, v;
	double	da, db, d1;
	double	Sqrt_v;

	double Integral;

	a	=	LossMin;
	b	=	LossMax;
	m	=	mean;
	v	=	variance;

	// v must be > 0
	// b-a must be > 0

	Sqrt_v	=	sqrt(v);
	da	=	(a - m) / Sqrt_v;
	db	=	(b - m) / Sqrt_v;
	d1	=	(1.0 - m) / Sqrt_v;

	double	Ncum_da;
	double	Ncum_db;
	double	Ncum_d1;

	Ncum_da	=	its_CopulaChoice->Cumulative_Density_Function(da);
	Ncum_db	=	its_CopulaChoice->Cumulative_Density_Function(db);
	Ncum_d1	=	its_CopulaChoice->Cumulative_Density_Function(d1);

	double	n_da;
	double	n_db;

	n_da	=	its_CopulaChoice->Density_Function(da);
	n_db	=	its_CopulaChoice->Density_Function(db);

	// first part: [0 ; a]
	Integral	=	0.0;
	// second part: [a ; b]
	Integral	+=	m * (Ncum_db - Ncum_da);

	//	StandardGaussianDensity(x) = 1 / SQRT(2 * PI) * exp(- x * x / 2.0)
	Integral	+=	Sqrt_v * (n_da - n_db);

	// third part: [b ; 1]	
	Integral	+=	Ncum_d1 * (b - a);
	Integral	+=	(a * Ncum_da - b * Ncum_db);

	return	Integral;
}


// --------------------------------------------------------
// Hermite Polynomials
// --------------------------------------------------------

// default computation is Degree (N) = 5

void CreditManager::Compute_LargePortfolioHermitePolynomials()
{
	int i;
	
	double TheDefProb;
	double TheDefProbSquare;
//	double TheDefProbSquareSquare;
	
	double	TheVariance;

	DoubleVector	DefProb_2;
	DoubleVector	DefProb_3;
//	DoubleVector	DefProb_4;
//	DoubleVector	DefProb_5;
	DoubleVector	TheVariances;

	DefProb_2.resize(Get_NbCredits());
	DefProb_3.resize(Get_NbCredits());
//	DefProb_4.resize(Get_NbCredits());
//	DefProb_5.resize(Get_NbCredits());
	TheVariances.resize(Get_NbCredits());

	double	current_LossAmount_Fraction;
	double	current_Variance;

	TheVariance	=	0.0;

	for (i=0; i<Get_NbCredits(); i++)
	{
		TheDefProb			=	its_Current_DefProb[i];
		TheDefProbSquare	=	TheDefProb * TheDefProb;
//		TheDefProbSquareSquare	= TheDefProbSquare * TheDefProbSquare;

		DefProb_2[i]	=	 TheDefProbSquare;
		DefProb_3[i]	=	 TheDefProbSquare * TheDefProb;
//		DefProb_4[i]	=	 TheDefProbSquareSquare;
//		DefProb_5[i]	=	 TheDefProbSquareSquare * TheDefProb;
		
		current_LossAmount_Fraction	=	its_LossesAmount[i] / TheBasketNotional;

		current_Variance	=	TheDefProb * (1.0 - TheDefProb);
		current_Variance	*=	current_LossAmount_Fraction * current_LossAmount_Fraction;

		TheVariances[i]	=	current_Variance;
		TheVariance	+=	current_Variance;
	}

	
	double	STL_TheVariance	=	accumulate(TheVariances.begin(), TheVariances.end(), 0.0);

	double	SqrtTheVariance;
	double	TheVarianceSquare;

	SqrtTheVariance	=	sqrt(TheVariance);
	TheVarianceSquare	=	TheVariance * TheVariance;

	// -----------------------------------
	// STEP 1
	// -----------------------------------
	// Credits Cumulants
	// -----------------------------------
	// kappa(i) (1,...,5)

	// -----------------------------------
	// STEP 2
	// -----------------------------------
	// Normalized Portfolio Cumulants
	// -----------------------------------
	// kappa (1,...,5)

	// -----------------------------------
	// STEP 3
	// -----------------------------------
	// Normalized Portfolio Raw Moments
	// -----------------------------------
	// mu_prime (1,...,5)

	DoubleVector	mu_primes;

	// _MAX_LP_Degree_ = 5
	mu_primes.resize(_MAX_LP_Degree_ + 1);

	// not used
	mu_primes[0]	=	0.0;

	mu_primes[1]	=	0.0;
	mu_primes[2]	=	1.0;

	double	current_Var;

	double	temp3	=	0.0;
	double	temp4	=	0.0;
	double	temp5	=	0.0;

	double	current_DefProb;
	double	current_DefProb_2;
	double	current_DefProb_3;

//	mu_prime[3]	=	inner_product(TheVariances.begin(), TheVariances.end(), 
	for (i=0; i<Get_NbCredits(); i++)
	{
		current_Var			=	TheVariances[i];
		current_LossAmount_Fraction	=	its_LossesAmount[i] / TheBasketNotional;

		current_DefProb		=	its_Current_DefProb[i];
		current_DefProb_2	=	DefProb_2[i];
		current_DefProb_3	=	DefProb_3[i];

		temp3	+=	current_Var * current_LossAmount_Fraction * (1.0 - 2 * current_DefProb);

		temp4	+=	current_Var * current_LossAmount_Fraction * current_LossAmount_Fraction *(1.0 - 6 * (current_DefProb - current_DefProb_2));

		temp5	+=	current_Var * current_LossAmount_Fraction * current_LossAmount_Fraction * current_LossAmount_Fraction * (1.0 - 14 * current_DefProb + 36 * current_DefProb_2 - 24 * current_DefProb_3);

	}

	temp3	/=	TheVariance;
	temp3	/=	SqrtTheVariance;
	mu_primes[3]	=	temp3;

	temp4	/=	TheVarianceSquare;
	temp4	+=	3;
	mu_primes[4]	=	temp4;
	
	temp5	/=	TheVarianceSquare;
	temp5	/=	SqrtTheVariance;
	temp5	+=	10 * temp3;
	mu_primes[5]	=	temp5;

	// -----------------------------------
	// STEP 4
	// -----------------------------------
	// Hermite coefficients: linear solving
	// -----------------------------------
	// c (0,...,5)

	its_CurrentHermiteCoeffs[0]	=	1.0;
	its_CurrentHermiteCoeffs[1]	=	0.0;
	its_CurrentHermiteCoeffs[2]	=	0.0;
	its_CurrentHermiteCoeffs[3]	=	mu_primes[3] / 6.0;
	its_CurrentHermiteCoeffs[4]	=	(mu_primes[4] - 3.0) / 24.0;
	its_CurrentHermiteCoeffs[5]	=	(-10.0 * mu_primes[3] + mu_primes[5]) / 120.0;

}


// --------------------------------------------------------
// Computation of TLoss(a,b)
// --------------------------------------------------------

// I avoid the 1.0 / (b-a) ratio
double CreditManager::LossLargePortfolioTranche_Hermite(double mean, double variance, double LossMin, double LossMax)
{
	double	a, b;
	double	m, v;
	double	da, db, d1;
	double	Sqrt_v;

	double Integral;

	a	=	LossMin;
	b	=	LossMax;
	m	=	mean;
	v	=	variance;

	// v must be > 0
	// b-a must be > 0

	Sqrt_v	=	sqrt(v);
	da	=	(a - m) / Sqrt_v;
	db	=	(b - m) / Sqrt_v;
	d1	=	(1.0 - m) / Sqrt_v;

	// --------------------------------
	// PRE-COMPUTATION
	// --------------------------------
	// BETA(i,da,db) and BETA(i,db,d1)
	// --------------------------------
	
	// for n(x) = StandardGaussianDensity(x)
	double	n_da;
	double	n_db;
	double	n_d1;

	double	NCum_da;
	double	NCum_db;
	double	NCum_d1;

	double	da_2;
	double	db_2;
	double	d1_2;
	double	da_4;
	double	db_4;
	double	d1_4;

	da_2	=	da * da;
	db_2	=	db * db;
	d1_2	=	d1 * d1;
	da_4	=	da_2 * da_2;
	db_4	=	db_2 * db_2;
	d1_4	=	d1_2 * d1_2;

	n_da	=	its_CopulaChoice->Density_Function(da);
	n_db	=	its_CopulaChoice->Density_Function(db);
	n_d1	=	its_CopulaChoice->Density_Function(d1);
	
	NCum_da	=	its_CopulaChoice->Cumulative_Density_Function(da);
	NCum_db	=	its_CopulaChoice->Cumulative_Density_Function(db);
	NCum_d1	=	its_CopulaChoice->Cumulative_Density_Function(d1);

	DoubleVector	Beta_i_a_b;
	DoubleVector	Beta_i_b_1;

	double	Beta_i_a_b_0;
	double	Beta_i_a_b_1;
	double	Beta_i_a_b_2;
	double	Beta_i_a_b_3;
	double	Beta_i_a_b_4;
	double	Beta_i_a_b_5;
	double	Beta_i_a_b_6;

	double	Beta_i_b_1_0;
	double	Beta_i_b_1_1;
	double	Beta_i_b_1_2;
	double	Beta_i_b_1_3;
	double	Beta_i_b_1_4;
	double	Beta_i_b_1_5;
	double	Beta_i_b_1_6;

	Beta_i_a_b.resize(_MAX_LP_Degree_ + 1);
	Beta_i_b_1.resize(_MAX_LP_Degree_ + 1);

	Beta_i_a_b_0	=	NCum_db - NCum_da;
	Beta_i_a_b_1	=	n_da - n_db;
	Beta_i_a_b_2	=	(da * n_da - db * n_db) + Beta_i_a_b_0;
	Beta_i_a_b_3	=	(da_2 * n_da - db_2 * n_db) + 2 * Beta_i_a_b_1;
	Beta_i_a_b_4	=	(da_2 * da * n_da - db_2 * db * n_db) + 3 * Beta_i_a_b_2;
	Beta_i_a_b_5	=	(da_4 * n_da - db_4 * n_db) + 4 * Beta_i_a_b_3;
	Beta_i_a_b_6	=	(da_4 * da * n_da - db_4 * db * n_db) + 5 * Beta_i_a_b_4;

	Beta_i_a_b[0]	=	Beta_i_a_b_0;
	Beta_i_a_b[1]	=	Beta_i_a_b_1;
	Beta_i_a_b[2]	=	Beta_i_a_b_2;
	Beta_i_a_b[3]	=	Beta_i_a_b_3;
	Beta_i_a_b[4]	=	Beta_i_a_b_4;
	Beta_i_a_b[5]	=	Beta_i_a_b_5;

	Beta_i_b_1_0	=	NCum_d1 - NCum_db;
	Beta_i_b_1_1	=	n_db - n_d1;
	Beta_i_b_1_2	=	(db * n_db - d1 * n_d1) + Beta_i_b_1_0;
	Beta_i_b_1_3	=	(db_2 * n_db - d1_2 * n_d1) + 2 * Beta_i_b_1_1;
	Beta_i_b_1_4	=	(db_2 * db * n_db - d1_2 * d1 * n_d1) + 3 * Beta_i_b_1_2;
	Beta_i_b_1_5	=	(db_4 * n_db - d1_4 * n_d1) + 4 * Beta_i_b_1_3;
	Beta_i_b_1_6	=	(db_4 * db * n_db - d1_4 * d1 * n_d1) + 5 * Beta_i_b_1_4;

	Beta_i_b_1[0]	=	Beta_i_b_1_0;
	Beta_i_b_1[1]	=	Beta_i_b_1_1;
	Beta_i_b_1[2]	=	Beta_i_b_1_2;
	Beta_i_b_1[3]	=	Beta_i_b_1_3;
	Beta_i_b_1[4]	=	Beta_i_b_1_4;
	Beta_i_b_1[5]	=	Beta_i_b_1_5;

	DoubleVector	Gamma_i_a_b;
	DoubleVector	Gamma_i_b_1;
	DoubleVector	Delta_i_a_b;
//	DoubleVector	Delta_i_b_1;

	Gamma_i_a_b.resize(_MAX_LP_Degree_ + 1);
	Gamma_i_b_1.resize(_MAX_LP_Degree_ + 1);
	Delta_i_a_b.resize(_MAX_LP_Degree_ + 1);
//	Delta_i_b_1.resize(_MAX_LP_Degree_ + 1);

	Gamma_i_a_b[0]	=	Beta_i_a_b_0;
	Gamma_i_a_b[1]	=	Beta_i_a_b_1;
	Gamma_i_a_b[2]	=	Beta_i_a_b_2 - Beta_i_a_b_0;
	Gamma_i_a_b[3]	=	Beta_i_a_b_3 - 3 * Beta_i_a_b_1;
	Gamma_i_a_b[4]	=	Beta_i_a_b_4 - 6 * Beta_i_a_b_2 + 3 * Beta_i_a_b_0;
	Gamma_i_a_b[5]	=	Beta_i_a_b_5 - 10 * Beta_i_a_b_3 + 15 * Beta_i_a_b_1;

	Gamma_i_b_1[0]	=	Beta_i_b_1_0;
	Gamma_i_b_1[1]	=	Beta_i_b_1_0;
	Gamma_i_b_1[2]	=	Beta_i_b_1_2 - Beta_i_b_1_0;
	Gamma_i_b_1[3]	=	Beta_i_b_1_3 - 3 * Beta_i_b_1_1;
	Gamma_i_b_1[4]	=	Beta_i_b_1_4 - 6 * Beta_i_b_1_2 + 3 * Beta_i_b_1_0;
	Gamma_i_b_1[5]	=	Beta_i_b_1_5 - 10 * Beta_i_b_1_3 + 15 * Beta_i_b_1_1;

	Delta_i_a_b[0]	=	Beta_i_a_b_1;
	Delta_i_a_b[1]	=	Beta_i_a_b_2;
	Delta_i_a_b[2]	=	Beta_i_a_b_3 - Beta_i_a_b_1;
	Delta_i_a_b[3]	=	Beta_i_a_b_4 - 3 * Beta_i_a_b_2;
	Delta_i_a_b[4]	=	Beta_i_a_b_5 - 6 * Beta_i_a_b_3 + 3 * Beta_i_a_b_1;
	Delta_i_a_b[5]	=	Beta_i_a_b_6 - 10 * Beta_i_a_b_4 + 15 * Beta_i_a_b_2;
/*
	Delta_i_b_1[0]	=	Beta_i_b_1_1;
	Delta_i_b_1[1]	=	Beta_i_b_1_1;
	Delta_i_b_1[2]	=	Beta_i_b_1_3 - Beta_i_b_1_1;
	Delta_i_b_1[3]	=	Beta_i_b_1_4 - 3 - Beta_i_b_1_2;
	Delta_i_b_1[4]	=	Beta_i_b_1_5 - 6 * Beta_i_b_1_3 + 3 * Beta_i_b_1_1;
	Delta_i_b_1[5]	=	Beta_i_b_1_6 - 10 * Beta_i_b_1_4 + 15 * Beta_i_b_1_2;
*/
	// --------------------------------
	// INTEGRALS
	// --------------------------------
	// first part: [0 ; a]
	Integral	=	0.0;

	// second part: [a ; b]
	double	Integral2;
	double	TmpIntegral;

	TmpIntegral	=	inner_product(its_CurrentHermiteCoeffs.begin(), its_CurrentHermiteCoeffs.begin() + its_LP_Degree + 1, 
						Delta_i_a_b.begin(), 0.0);

	Integral2	=	TmpIntegral;
	Integral2	*=	Sqrt_v;


	TmpIntegral	=	inner_product(its_CurrentHermiteCoeffs.begin(), its_CurrentHermiteCoeffs.begin() + its_LP_Degree + 1, 
						Gamma_i_a_b.begin(), 0.0);

	TmpIntegral	*=	(m - a);

	Integral2	+=	TmpIntegral;

	// third part: [b ; 1]	
	double	Integral3;
	
	Integral3	=	inner_product(its_CurrentHermiteCoeffs.begin(), its_CurrentHermiteCoeffs.begin() + its_LP_Degree + 1, 
						Gamma_i_b_1.begin(), 0.0);

	Integral3	*=	(b-a);

	Integral	+=	Integral2;
	Integral	+=	Integral3;

	return	Integral;
}