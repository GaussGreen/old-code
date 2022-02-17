#include "ARMKernel\glob\firsttoinc.h"

/*********************************************************************************/
/*! \class  CreditManager CreditManager_FFT_1F.cpp "CreditManager_FFT_1F.cpp"
 *  \author:	L. JACQUEL 
 *	\version:	1.0
 *	\date:		June 2005
 *	\brief:		Implementation of the 1 Factor Copula Model - JPL & JG paper:
				Paper reference Basket Default Swaps, CDO’s and Factor Copulas
/***********************************************************************************/


#include "ICMKernel\cair\credit_manager.h"

#include "ICMKernel\glob\icm_maths.h"

#include "ICMKernel\util\icm_fft.h"

#include <nags.h>		//	s15abc
#include <nagg01.h>		//	g01fac

// ------------------------------------------------------------------------
//
//	Computes Probability that at time T we have more than Nth to Default (included)
//		i.e. : (1 - Sum (i=0,NtD-1) Pi)
//		where	Pi is the 
//
// ------------------------------------------------------------------------

void	CreditManager::ProbabilityOfAtLeastNDefaultBeforeT_1F(double	T, int NtD, double& Result)
{
	int		i,j; 
	double	z, pp, zmin, zmax; 
	double	TheValue, SurvProb;
	double	TheCorrelation;
	
	DoubleVector	ProbSurvOneName;
	DoubleVector	X;
	DoubleVector	ProbNtD;
	DoubleVector	ProbNtD_z;

	// resize arrays
	ProbSurvOneName.resize(Get_NbCredits());
	X.resize(Get_NbCredits());

	// 0 to Get_NbCredits() included
	ProbNtD.resize(Get_NbCredits()+1);
	ProbNtD_z.resize(Get_NbCredits()+1);

	if (T <= 0.)
	{
		Result = 0.0;
		return;
	}

	// ----------------------------------
	// According to Copula Type
	// ----------------------------------
	switch (its_CopulaType)
	{
	case CCT_GAUSSIAN:
	
		// Computes N-1 (SurvivalProb(T))	(Inverse Cumulative Gaussian Distribution)
		for (i=0; i<Get_NbCredits(); i++)
		{
			// Curve i
			// Default Proba expects only yearterms
			SurvProb = (its_ArrayDefaultCrv[i])->SurvivalProba(T);	// T is already in year fraction

			// should test if DefProb = 0.0 or = 1.0, so SurvProb = 1.0 - DefProb
			if (SurvProb == 0.0)
				TheValue	=	_MINUS_INFINITY_;	//	minus infinity
			else if (SurvProb == 1.0)
				TheValue	=	_PLUS_INFINITY_;	//	plus infinity
			else
				TheValue	=	NAG_deviates_normal(SurvProb);

			X[i] = TheValue;
		}
	
		// RAZ to 0.0
		fill(ProbNtD.begin(), ProbNtD.end(), 0.0);
		
		// integration limits (standard error)
		zmin	=	-6.0;
		zmax	=	6.0;
		double	integration_size;
		integration_size	=	(zmax - zmin);

		for (i=1 ;i<=its_NIntegration_1F; i++)
		{
			z = integration_size * its_Xi[i] + zmin;	
			
			// Knoits_Wing Factor value z, Evaluate Survival Prob for each name at time T
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
					//	Normal Cumulative Function
					ProbSurvOneName[j] = NAG_cumul_normal(pp / sqrt(1.0 - TheCorrelation));
				else
					if (pp > 0.0)	// + infinite => N(+infinite)
						ProbSurvOneName[j] = 1.0;
					else			// - infinite => N(-infinite)
						ProbSurvOneName[j] = 0.0;
			}

			// Probabilities from 0 to Get_NbCredits() to default for given Factor value z
			CumulativeBinomialDistibutionNumbers(ProbSurvOneName, ProbNtD_z);

			// GaussLegendre Integration for each Nth to Default
			TheValue = StandardGaussianDensity(z) * its_Wi[i] * integration_size;

			for (j=0; j<NtD; j++)
				ProbNtD[j] += ProbNtD_z[j] * TheValue;
		}
	
	break;

	case CCT_STUDENT:

		ICMTHROW(ERR_INVALID_DATA,"Student Copula not yet implemented for 1F Model!");
		break;

	default:
		ICMTHROW(ERR_INVALID_DATA,"Unknown Copula Type!");
		break;
	}

	// Sum all the probabilities from 0 to N-1 default states and then 1-P 
//	if (NtD)
//		Result	=	1.0	- accumulate(ProbNtD.begin(), ProbNtD..at(NtD-1), 0);

	Result = 0.;
	for (j=0; j<NtD;j++)
		Result += ProbNtD[j];
	Result = 1. - Result;

}



// ------------------------------------------------------------------------
//
//	Computes Probability that at time T we have up to Nth Defaults (included)
//		where	Pi is the individual Default Probabilities
//
// ------------------------------------------------------------------------

void	CreditManager::ProbabilityOfNDefaultsBeforeT_1F(double	T, int NtD, DoubleVector& Result)
{
	int		i,j; 
	double	z, pp, zmin, zmax; 
	double	TheValue, SurvProb;
	double	TheCorrelation;
	
	DoubleVector	ProbSurvOneName;
	DoubleVector	X;
	DoubleVector	ProbNtD;
	DoubleVector	ProbNtD_z;

	// resize arrays
	ProbSurvOneName.resize(Get_NbCredits());
	X.resize(Get_NbCredits());

	// 0 to Get_NbCredits() included
	ProbNtD.resize(Get_NbCredits()+1);
	ProbNtD_z.resize(Get_NbCredits()+1);

/*
	Result.clear();
	Result.resize(NtD+1);
*/

	if (T <= 0.)
	{
		fill(Result.begin(), Result.end(), 0.0);
		return;
	}

	// ----------------------------------
	// According to Copula Type
	// ----------------------------------
	switch (its_CopulaType)
	{
	case CCT_GAUSSIAN:
	
		// Computes N-1 (SurvivalProb(T))	(Inverse Cumulative Gaussian Distribution)
		for (i=0; i<Get_NbCredits(); i++)
		{
			// Curve i
			// Default Proba expects only yearterms
			SurvProb = (its_ArrayDefaultCrv[i])->SurvivalProba(T);	// T is already in year fraction

			// should test if DefProb = 0.0 or = 1.0, so SurvProb = 1.0 - DefProb
			if (SurvProb == 0.0)
				TheValue	=	_MINUS_INFINITY_;	//	minus infinity
			else if (SurvProb == 1.0)
				TheValue	=	_PLUS_INFINITY_;	//	plus infinity
			else
				TheValue	=	NAG_deviates_normal( SurvProb);

			X[i] = TheValue;
		}

	
		// RAZ to 0.0
		fill(ProbNtD.begin(), ProbNtD.end(), 0.0);
		
		// integration limits (standard error)
		zmin	=	-6.0;
		zmax	=	6.0;
		double	integration_size;
		integration_size	=	(zmax - zmin);

		for (i=1 ;i<=its_NIntegration_1F; i++)
		{
			z = integration_size * its_Xi[i] + zmin;	
			
			// Knoits_Wing Factor value z, Evaluate Survival Prob for each name at time T
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
					//	Normal Cumulative Function
					ProbSurvOneName[j] = NAG_cumul_normal(pp / sqrt(1.0 - TheCorrelation));
				else
					if (pp > 0.0)	// + infinite => N(+infinite)
						ProbSurvOneName[j] = 1.0;
					else			// - infinite => N(-infinite)
						ProbSurvOneName[j] = 0.0;
			}

			// Probabilities from 0 to Get_NbCredits() to default for given Factor value z
			CumulativeBinomialDistibutionNumbers(ProbSurvOneName, ProbNtD_z);

			// GaussLegendre Integration for each Nth to Default
			TheValue = StandardGaussianDensity(z) * its_Wi[i] * integration_size;

			for (j=0; j<=NtD; j++)
				ProbNtD[j] += ProbNtD_z[j] * TheValue;
		}
	
	break;

	case CCT_STUDENT:

		ICMTHROW(ERR_INVALID_DATA,"Student Copula not yet implemented for 1F Model!");
		break;

	default:
		ICMTHROW(ERR_INVALID_DATA,"Unknown Copula Type!");
		break;
	}

	// Take all the probabilities from 0 to N default states

	double	Sum	=	0.0;

	for (j=0; j<=NtD;j++)
	{
		TheValue	=	ProbNtD[j];
		Result[j]	=	TheValue;
		Sum			+=	TheValue;
	}

	// Finally the complementary part
	Result[NtD+1]	=	1.0 - Sum;
}


// -----------------------------------------------------------------------------------------
// Hull - APPENDIX A Algorithm
//
// This function returns in the Qk vector the probability of having from k a sum of m = Dim 
// Bernouilli variables - whose probabilies are in the p vector
// -----------------------------------------------------------------------------------------

void CreditManager::CumulativeBinomialDistibutionNumbers(DoubleVector& p, DoubleVector& Qk)
{
	int	i, j, N;
	
	N = p.size();
	DoubleVector	Table;
	Table.resize(N+1);

	if (N > Get_NbCredits())
		ICMTHROW(ERR_INVALID_DATA," CumulativeBinomialDistibutionNumbers can not exceed Number of Credits! " << N << " instead of " << Get_NbCredits() << " expected!");

	fill(Qk.begin(), Qk.end(), 0.0);
	Qk[0] = 1.0;

	for (i=0; i<N; i++)
	{
		Table[0] = Qk[0] * p[i];
		for (j=1; j<=i; j++) 
		{
			Table[j] = p[i] * Qk[j] + (1.0 - p[i]) * Qk[j - 1];
			Qk[j - 1] = Table[j - 1];
		}
		Table[i + 1] = Qk[i] * (1.0 - p[i]);
		Qk[i] = Table[i];
		Qk[i + 1] = Table[i + 1];
	}

}



// -----------------------------------------------------------------------------------------
// -----------------------------------------------------------------------------------------
//
//		FFT
//
// -----------------------------------------------------------------------------------------
// -----------------------------------------------------------------------------------------

// The PortfolioLoss is defined as the sum of independent Bernouilli random variables weighted by Losses i.e.:
// Loss = L[0] * B(0) + ... + L[n-1] * B(n-1) its_With Proba[B(i) = 0] = DefProb[i]
// Return Output is in Psi real and imaginary parts of the Characteristic Function of Loss at point u.

void	CreditManager::PortfolioLoss_CharacteristicFunction(double u, DoubleVector& DefProb, double& psi_r, double& psi_i)
{
	int		i;
	double	a, b, c, d, l;


	l = its_LossesAmount[0] * u;

	a = DefProb[0] * (cos(l) - 1.0) + 1.0;
	b = DefProb[0] * sin(l);

	psi_r = a;
	psi_i = b;
	
	for (i=1 ;i<Get_NbCredits() ;i++)
	{
		l = its_LossesAmount[i] * u;
		a = DefProb[i] * (cos(l) - 1.0) + 1.0;
		b = DefProb[i] * sin(l);
		
		c = psi_r;
		d = psi_i;
		
		psi_r = a * c - b * d; 
		psi_i = a * d + b * c;
	}

}


//
// PortfolioLoss Computation, computed via FFT
//
void	CreditManager::PortfolioLoss_FFT(DoubleVector& DefProb, DoubleVector& ProbaDensityFunctionLossDistribution)
{
	long	i;		
	double	psi_r, psi_i, wr;
	double	TwoPiOverTotalLoss, OneOverNFFT;

	double	*m_datas;
	m_datas	=	new double [2 * its_N_FFT + 3];
	
	// total Loss represents the size
	TwoPiOverTotalLoss = TWOPI / its_TotalLoss;
	OneOverNFFT = 1. / its_N_FFT;
	
	wr = 0.0; 
	PortfolioLoss_CharacteristicFunction(wr, DefProb, psi_r, psi_i);

	// first point
	m_datas[1] = psi_r; 
	m_datas[2] = psi_i; 
	
	for (i=3; i<=its_N_FFT-1; i+=2)
	{
		// next point: 2 PI convention (cf. NR)
		// should be pre-computed!
		wr = wr + TwoPiOverTotalLoss;
		PortfolioLoss_CharacteristicFunction(wr, DefProb, psi_r, psi_i);

		m_datas[i]					= psi_r;   
		m_datas[i + 1]				= psi_i;
		
		// we use symmetry (the density is real), so that its FFT has the property: FFT(-u) = FFT(u)*
		// where * means the complement
		m_datas[its_N_FFT * 2 + 2 - i]	= psi_r;
		m_datas[its_N_FFT * 2 + 3 - i]	= -psi_i;
	}

	wr = wr + TwoPiOverTotalLoss;
	PortfolioLoss_CharacteristicFunction(wr, DefProb, psi_r, psi_i);
	m_datas[its_N_FFT+1] = psi_r; 
	m_datas[its_N_FFT+2] = psi_i;

	// Compute the fft inverse:
	fourier::fft(m_datas, its_N_FFT, -1);

	// Fill density vector
	for (i=1; i<=its_N_FFT; i++) 
		ProbaDensityFunctionLossDistribution[i-1] = m_datas[2*i-1] * OneOverNFFT;

	delete [] m_datas;
}



// ------------------------------------------------------------------------
//
//	Computes Probability that at time T we have 
//		LossMIN < Loss < LossMAX
//
//	Algorithm is  FFT
// 
// ------------------------------------------------------------------------

void	CreditManager::ComputeProbabilityOrExpectationLossesBeforeT_FFT_1F(double	T, CreditLossComputation ProbOrExpectationFlag, double LossMin, double LossMax, DoubleVector& Result)
{
	int		i,j; 
	double	z, pp, zmin, zmax; 
	double	TheValue, SurvProb;
	double	TheCorrelation;
	
	// Gauss Legendre Integration
	DoubleVector	ProbSurvOneName;
	DoubleVector	X;
	
	// FFT Integration
	DoubleVector	PDF_Loss_Distribution;
	DoubleVector	PDF_Loss_Distribution_z;

	// resize arrays
	ProbSurvOneName.resize(Get_NbCredits());
	X.resize(Get_NbCredits());

	PDF_Loss_Distribution.resize(its_N_FFT);
	PDF_Loss_Distribution_z.resize(its_N_FFT);

	if (T <= 0.)
	{
		fill(Result.begin(), Result.end(), 0.0);
		return;
	}

	string	TheLabel;
	char	buffer[20];

	// ----------------------------------
	// According to Copula Type
	// ----------------------------------
	switch (its_CopulaType)
	{
	case CCT_GAUSSIAN:
	
		// Computes N-1 (SurvivalProb(T))	(Inverse Cumulative Gaussian Distribution)
		for (i=0; i<Get_NbCredits(); i++)
		{
			// Curve i
			// Default Proba expects only yearterms
			SurvProb = (its_ArrayDefaultCrv[i])->SurvivalProba(T);	// T is already in year fraction

			// should test if DefProb = 0.0 or = 1.0, so SurvProb = 1.0 - DefProb
			if (SurvProb == 0.0)
				TheValue	=	_MINUS_INFINITY_;	//	minus infinity
			else if (SurvProb == 1.0)
				TheValue	=	_PLUS_INFINITY_;	//	plus infinity
			else
				TheValue	=	NAG_deviates_normal( SurvProb);

			X[i] = TheValue;
		}

	
		// RAZ to 0.0
		fill(PDF_Loss_Distribution.begin(), PDF_Loss_Distribution.end(), 0.0);

		// integration limits (standard error)
		zmin	=	-6.0;
		zmax	=	6.0;
		double	integration_size;
		integration_size	=	(zmax - zmin);

		for (i=1 ;i<=its_NIntegration_1F; i++)
		{
			z = integration_size * its_Xi[i] + zmin;	
			
			// Knoits_Wing Factor value z, Evaluate Survival Prob for each name at time T
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
					//	Normal Cumulative Function
					ProbSurvOneName[j] = NAG_cumul_normal(pp / sqrt(1.0 - TheCorrelation));
				else
					if (pp > 0.0)	// + infinite => N(+infinite)
						ProbSurvOneName[j] = 1.0;
					else			// - infinite => N(-infinite)
						ProbSurvOneName[j] = 0.0;
			}

			// Compute Loss Probability
			PortfolioLoss_FFT(ProbSurvOneName, PDF_Loss_Distribution_z);
			
			// GaussLegendre Integration for each Nth to Default
			TheValue = StandardGaussianDensity(z) * its_Wi[i] * integration_size;

			// integration at point z for each discretized FFT
			for (j=0; j<its_N_FFT; j++)
				PDF_Loss_Distribution[j] += PDF_Loss_Distribution_z[j] * TheValue;

			TheLabel	=	"Integration step: ";
			itoa(i, buffer, 10);
			TheLabel.append(buffer);

			Display_Vector(TheLabel, PDF_Loss_Distribution);
		}
	
	break;

	case CCT_STUDENT:

		ICMTHROW(ERR_INVALID_DATA,"Student Copula not yet implemented for 1F Model!");
		break;

	default:
		ICMTHROW(ERR_INVALID_DATA,"Unknown Copula Type!");
		break;
	}

	// According to the Choice, either Maturity (one date)
	// or density
	
	// Extract from PDF (Probability Density Function)
	// the probability of Loss lying between LossMIN and LossMax
	CreditCalibratorSelection	TheCreditCalibratorChoice;
	GetCreditCalibratorChoice(TheCreditCalibratorChoice);

	switch (TheCreditCalibratorChoice)
	{
	case CCS_NLOSS:
		GetLossFromPDF_Distribution(ProbOrExpectationFlag, LossMin, LossMax, PDF_Loss_Distribution, Result[0]);

		break;

	case CCS_DENSITY_LOSS:

		double	LowBound	=	LossMin;
		double	UpBound		=	LowBound;
		double	TrancheSize;
		double	TheStep;

		TrancheSize	=	(LossMax - LossMin);
		TheStep		=	its_CreditCalibrator_Density_Loss_Step * TrancheSize;

		j	=	0;
		if (CHECK_EQUAL(LossMin, 0.0))
			// LossMin is equal to 0.0
			Result[j]	=	0.0;
		else
			GetLossFromPDF_Distribution(ProbOrExpectationFlag, 0.0, LowBound, PDF_Loss_Distribution, Result[j]);
		
		j	=	size_DENSITY_LOSS + 1;
		if (CHECK_EQUAL(LossMax, its_TotalLoss))
			// LossMax is equal to 1000.0
			Result[j]	=	0.0;
		else
			GetLossFromPDF_Distribution(ProbOrExpectationFlag, LossMax, its_TotalLoss, PDF_Loss_Distribution, Result[j]);

		// all remaining ones
		j	=	0;
		GetLossFromPDF_Distribution(ProbOrExpectationFlag, LowBound, LowBound, PDF_Loss_Distribution, Result[j]);

		while (UpBound <= LossMax)
		{
			UpBound	+=	TheStep;
			GetLossFromPDF_Distribution(ProbOrExpectationFlag, LowBound, UpBound, PDF_Loss_Distribution, Result[j]);
			LowBound	=	UpBound;
			j++;
		}

		break;
	}
}


void	CreditManager::GetLossFromPDF_Distribution(CreditLossComputation ProbOrExpectationFlag, double LossMin, double LossMax, DoubleVector& PDF_Loss_Distribution, double& Prob_Loss)
{
	int	i, iStartNumInteg, iEndNumInteg;
	double	TheIntegral;
	double	The_EPS	=	1.0e-6;	// just to be sure to include or exclude the ranges

	// computes according to LossMin and LossMax, the numerical indices (should be done outside the integration loop)
	
	// no more than the Total Loss
	double	CutOff	=	FMIN(LossMax, its_TotalLoss);

	// the step (divide the [0; TotalLoss] interval its_With its_N_FFT points
	double	dL	=	its_TotalLoss / (double) (its_N_FFT-1);

	iStartNumInteg	=	(int) ((LossMin - The_EPS) / dL);
	iEndNumInteg	=	(int) ((CutOff - The_EPS) / dL);

	TheIntegral	=	0.0;

	// Numerical Integration: Sum (Loss-LossMin) * PDFLossDistribution * dLi

	switch (ProbOrExpectationFlag)
	{
	case CLC_EXPECTATION:

		for (i=iStartNumInteg; i<iEndNumInteg; i++)
			TheIntegral	+=	PDF_Loss_Distribution[i] * (dL * i - LossMin);
		// and the final one
		TheIntegral	+= PDF_Loss_Distribution[iEndNumInteg] * (CutOff - LossMin);

		break;

	case CLC_PROBABILITY:

		if (iStartNumInteg	==	iEndNumInteg)
			TheIntegral	+=	PDF_Loss_Distribution[iStartNumInteg];
		else
			for (i=iStartNumInteg+1; i<=iEndNumInteg; i++)
				TheIntegral	+=	PDF_Loss_Distribution[i];

		break;

	default:
		
		ICMTHROW(ERR_INVALID_DATA,"Unknown Loss Computation: Only Probability or Expectation are yet implemented for 1F Model!");
	}

	Prob_Loss	=	TheIntegral;
}