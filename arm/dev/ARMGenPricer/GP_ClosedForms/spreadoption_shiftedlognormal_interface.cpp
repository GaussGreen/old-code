/*!
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 * Version initiale 01/15/2004
 *
 *  basic functions for the closed form framework 
 *
 *	\file spreadoption_shiftedlognormal.cpp
 *
 *  \brief
 *
 *	\author  O. Croissant
 *	\version 1.0
 *	\date January 2004
 */
#include <glob/firsttoinc.h>

#include "gpbase/port.h"

#include "gpclosedforms/spreadoption_shiftedlognormal_formula.h"


/// uses the template PowerSpreadOption_Pricing_With_Limits which declared in powerspreadoption.h



CC_BEGIN_NAMESPACE(ARM)


///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
///  
///			Begining Exportable  Pricing Functions 
///
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////



/// callput =  1 (K_CALL) for call
/// callput = -1 (K_PUT) for put

double Export_Gaussian_ShiftedLN_Power_SpreadOption(double S1,double S2,
								   double sigma1, double alpha1,
								   double sigma2, double alpha2,
								   double copula_corr,double t,
								   double a10,double b10,double k10,double a20,double b20,double k20,int n)
{
	ArgumentList a( S1, S2,
					sigma1,alpha1,
					sigma2,alpha2,
					copula_corr, t,
					a10, b10, k10, a20, b20, k20, n);

	Power_Expression<ARM_CF_ShiftedLN_Gaussian_PowerSpreadOption_Formula> y;
	return y(a);
}

///////////////////////////////////////////////////////////////////////
///  
///			1st Derivatives  
///
///////////////////////////////////////////////////////////////////////
double Export_Gaussian_ShiftedLN_Power_SpreadOption(int i,
								   double S1,double S2,
								   double sigma1, double alpha1,
								   double sigma2, double alpha2,
								   double copula_corr,double t,
								   double a10,double b10,double k10,double a20,double b20,double k20,int n)
{
	ArgumentList a( S1, S2,
					sigma1,alpha1,
					sigma2,alpha2,
					copula_corr, t,
					a10, b10, k10, a20, b20, k20, n);
	
	Power_Expression<ARM_CF_ShiftedLN_Gaussian_PowerSpreadOption_Formula> y;
	return y(i,a);
}

///////////////////////////////////////////////////////////////////////
///  
///			2nd Derivatives  
///
///////////////////////////////////////////////////////////////////////
double Export_Gaussian_ShiftedLN_Power_SpreadOption(int i,int j,
								   double S1,double S2,
								   double sigma1, double alpha1,
								   double sigma2, double alpha2,
								   double copula_corr,double t,
								   double a10,double b10,double k10,double a20,double b20,double k20,int n)
{
	ArgumentList a( S1, S2,
					sigma1,alpha1,
					sigma2,alpha2,
					copula_corr, t,
					a10, b10, k10, a20, b20, k20, n);
	
	Power_Expression<ARM_CF_ShiftedLN_Gaussian_PowerSpreadOption_Formula> y;
	return y(i,j,a);
}

///////////////////////////////////////////////////////////////////////
///  
///			Certitude  
///
///////////////////////////////////////////////////////////////////////
double Export_Gaussian_ShiftedLN_Power_SpreadOption_Certitude(
										 double copula_corr,double t,int n)
{
ArgumentList a(copula_corr);
return ARM_CF_PowerSpreadOption_Formula<ShiftedLogNormal_Smile,GaussianCopula>::Certitude(a,t,n);

}


CC_END_NAMESPACE()
 


/*---------------------------------------------------------------------------*/
/*---- End of file ----*/
