#include "ARMKernel\glob\firsttoinc.h"
/* -----------------------------------------------------------------
   -----------------------------------------------------------------

	FILE:		ICM_Probability_Density_Gaussian.H
	PROJECT:	CAIR
	
	DESCRIPTION:	this class provides the ancestor for Gaussian Probability Density


   -----------------------------------------------------------------
   
	ICM KERNEL Library

		version 1.0
		developped by Laurent Jacquel

  -----------------------------------------------------------------
  ----------------------------------------------------------------- */

#include "ICMKernel\cair\ICM_Probability_Density_Gaussian.h"

# include "math.h"
#include "ICMKernel\glob\icm_maths.h"


void ICM_Probability_Density_Gaussian :: Init()
{
}

ICM_Probability_Density* ICM_Probability_Density_Gaussian :: Clone() const
{
	ICM_Probability_Density_Gaussian *New;
	New = new ICM_Probability_Density_Gaussian(*this);
	return New;
}

void ICM_Probability_Density_Gaussian :: Copy(const ICM_Probability_Density& Data)
{
	mean		=	((ICM_Probability_Density_Gaussian&)Data).mean;
	variance	=	((ICM_Probability_Density_Gaussian&)Data).variance;
}


void ICM_Probability_Density_Gaussian :: SetParameters(double data_mean, double data_variance)
{
	mean		=	data_mean;
	variance	=	data_variance;
}


		
// ----------------------------------------------------------- 
// Density Function
// ----------------------------------------------------------- 

double	ICM_Probability_Density_Gaussian :: Density_Function(double x)
{
	double	Sqrt_variance;
	double	y;


	Sqrt_variance = sqrt(variance);
	y	=	(x - mean) / Sqrt_variance;
	
	return (StandardGaussianDensity(y) / Sqrt_variance);
}


// -----------------------------------------------------------
// Cumulative Density Function
// -----------------------------------------------------------

double	ICM_Probability_Density_Gaussian :: Cumulative_Density_Function(double x)
{
	double	Sqrt_variance;
	double	y;

	Sqrt_variance = sqrt(variance);
	y	=	(x - mean) / Sqrt_variance;

	return NAG_cumul_normal(y);
}


// -----------------------------------------------------------
// Inverse Cumulative Density Function
// -----------------------------------------------------------

double	ICM_Probability_Density_Gaussian :: Inverse_Cumulative_Density_Function(double x)
{
	double	Sqrt_variance;
	double	y;

	Sqrt_variance = sqrt(variance);

	y	=	NAG_deviates_normal( x);

	return 	(y * Sqrt_variance + mean);
}
	
