#include "ARMKernel\glob\firsttoinc.h"
/* -----------------------------------------------------------------
   -----------------------------------------------------------------

	FILE:		ICM_Probability_Density_Standard_Gaussian.H
	PROJECT:	CAIR
	
	DESCRIPTION:	this class provides the ancestor for Gaussian Probability Density


   -----------------------------------------------------------------
   
	ICM KERNEL Library

		version 1.0
		developped by Laurent Jacquel

  -----------------------------------------------------------------
  ----------------------------------------------------------------- */

#include "ICMKernel\cair\ICM_Probability_Density_Standard_Gaussian.h"

#include "ICMKernel\glob\icm_maths.h"



void ICM_Probability_Density_Standard_Gaussian :: Init()
{
}

ICM_Probability_Density* ICM_Probability_Density_Standard_Gaussian :: Clone() const
{
	ICM_Probability_Density_Standard_Gaussian *New;
	New = new ICM_Probability_Density_Standard_Gaussian(*this);
	return New;
}

void ICM_Probability_Density_Standard_Gaussian :: Copy(const ICM_Probability_Density& Data)
{

}


void ICM_Probability_Density_Standard_Gaussian :: SetParameters(double data_mean, double data_variance)
{
}


		
// ----------------------------------------------------------- 
// Density Function
// ----------------------------------------------------------- 

double	ICM_Probability_Density_Standard_Gaussian :: Density_Function(double x)
{
	return StandardGaussianDensity(x);
}


// -----------------------------------------------------------
// Cumulative Density Function
// -----------------------------------------------------------

double	ICM_Probability_Density_Standard_Gaussian :: Cumulative_Density_Function(double x)
{
	return NAG_cumul_normal(x);
}


// -----------------------------------------------------------
// Inverse Cumulative Density Function
// -----------------------------------------------------------

double	ICM_Probability_Density_Standard_Gaussian :: Inverse_Cumulative_Density_Function(double x)
{
	return NAG_deviates_normal( x );
}
	
