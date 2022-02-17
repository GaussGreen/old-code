/*!
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 * Version initiale 01/15/2004
 *
 *  basic functions for the closed form framework 
 *
 *	\file change_numeraire.cpp
 *
 *  \brief
 *
 *	\author  O. Croissant
 *	\version 1.0
 *	\date January 2004
 */
#include <glob/firsttoinc.h>
#include "gpbase/port.h"

#include "gpclosedforms/spreadoption_sabr.h"
#include "gpclosedforms/spreadoption_sabr_formula.h"
#include "gpclosedforms/spreadoption_sabr_student_formula.h"



CC_BEGIN_NAMESPACE(ARM)
//////////////////////////////////////////////////////////////////////////////////////////
///
///   Pricing Functions
///
//////////////////////////////////////////////////////////////////////////////////////////


///////////////////////////////////////////////////////////////////////
///  
///		Gaussian Copula, Pricing Functions
///
///////////////////////////////////////////////////////////////////////


/// callput =  1  for call
/// callput = -1  for put

double Export_Gaussian_SABR_Power_SpreadOption(double S1,double S2,
								   double alpha1, double beta1, double rho1, double nu1,
								   double alpha2, double beta2, double rho2, double nu2,
								   double copula_corr,double t,int flag,
								   double a10,double b10,double k10,double a20,double b20,double k20,int n,double alpha_exp,double alpha_tanh,double kb_tanh)
{
	ArgumentList a( S1, S2,
					alpha1,  beta1,  rho1,  nu1,
					alpha2,  beta2,  rho2,  nu2,
					copula_corr, t,
					a10, b10, k10, a20, b20, k20,alpha_exp,alpha_tanh ,kb_tanh,flag,n,GENERIC_SPREADOPTION);

	Power_Expression<ARM_CF_SABR_Gaussian_PowerSpreadOption_Formula> y;
	return y(a);
}

double Export_Gaussian_SABR_Power_Digital_SpreadOption(double S1,double S2,
								   double alpha1, double beta1, double rho1, double nu1,
								   double alpha2, double beta2, double rho2, double nu2,
								   double copula_corr,double t,int flag,
								  double k10,double a20,double b20,double k20,int n,double alpha_exp,double alpha_tanh,double kb_tanh)
{
	ArgumentList a( S1, S2,
					alpha1,  beta1,  rho1,  nu1,
					alpha2,  beta2,  rho2,  nu2,
					copula_corr, t,
					0, 0, k10, a20, b20, k20, alpha_exp,alpha_tanh,kb_tanh,flag,n,DIGITAL_SPREADOPTION);

	Power_Expression<ARM_CF_SABR_Gaussian_PowerSpreadOption_Formula> y;
	return y(a);
}



///////////////////////////////////////////////////////////////////////
///  
///			Gaussian Copula, 	1st Derivatives  
///
///////////////////////////////////////////////////////////////////////
double Export_Gaussian_SABR_Power_SpreadOption(int i,
								   double S1,double S2,
								   double alpha1, double beta1, double rho1, double nu1,
								   double alpha2, double beta2, double rho2, double nu2,
								   double copula_corr,double t,int flag,
								   double a10,double b10,double k10,double a20,double b20,double k20,int n,double alpha_exp,double alpha_tanh,double kb_tanh)
{
	ArgumentList a( S1, S2,
		alpha1,  beta1,  rho1,  nu1,
		alpha2,  beta2,  rho2,  nu2,
		copula_corr, t,
		a10, b10, k10, a20, b20, k20,alpha_exp,alpha_tanh,kb_tanh,flag, n,GENERIC_SPREADOPTION);
	
	Power_Expression<ARM_CF_SABR_Gaussian_PowerSpreadOption_Formula> y;
	return y(i,a);
}

double Export_Gaussian_SABR_Power_Digital_SpreadOption(int i,
								   double S1,double S2,
								   double alpha1, double beta1, double rho1, double nu1,
								   double alpha2, double beta2, double rho2, double nu2,
								   double copula_corr,double t,int flag,
								   double k10,double a20,double b20,double k20,int n,double alpha_exp,double alpha_tanh,double kb_tanh)
{
	ArgumentList a( S1, S2,
		alpha1,  beta1,  rho1,  nu1,
		alpha2,  beta2,  rho2,  nu2,
		copula_corr, t,
		0, 0, k10, a20, b20, k20,alpha_exp,alpha_tanh,kb_tanh,flag,n,DIGITAL_SPREADOPTION);
	
	Power_Expression<ARM_CF_SABR_Gaussian_PowerSpreadOption_Formula> y;
	return y(i,a);
}

///////////////////////////////////////////////////////////////////////
///  
///				Gaussian Copula,  2nd Derivatives  
///
///////////////////////////////////////////////////////////////////////
double Export_Gaussian_SABR_Power_SpreadOption(int i,int j,
								   double S1,double S2,
								   double alpha1, double beta1, double rho1, double nu1,
								   double alpha2, double beta2, double rho2, double nu2,
								   double copula_corr,double t,int flag,
								   double a10,double b10,double k10,double a20,double b20,double k20,int n,double alpha_exp,double alpha_tanh,double kb_tanh)
{
	ArgumentList a( S1, S2,
		alpha1,  beta1,  rho1,  nu1,
		alpha2,  beta2,  rho2,  nu2,
		copula_corr, t,
		a10, b10, k10, a20, b20, k20, alpha_exp,alpha_tanh,kb_tanh,flag,n,GENERIC_SPREADOPTION);
	
	Power_Expression<ARM_CF_SABR_Gaussian_PowerSpreadOption_Formula> y;
	return y(i,j,a);
}

double Export_Gaussian_SABR_Power_Digital_SpreadOption(int i,int j,
								   double S1,double S2,
								   double alpha1, double beta1, double rho1, double nu1,
								   double alpha2, double beta2, double rho2, double nu2,
								   double copula_corr,double t,int flag,
								   double k10,double a20,double b20,double k20,int n,double alpha_exp,double alpha_tanh,double kb_tanh)
{
	ArgumentList a( S1, S2,
		alpha1,  beta1,  rho1,  nu1,
		alpha2,  beta2,  rho2,  nu2,
		copula_corr, t,
		0, 0, k10, a20, b20, k20, alpha_exp,alpha_tanh,kb_tanh,flag,n,DIGITAL_SPREADOPTION);
	
	Power_Expression<ARM_CF_SABR_Gaussian_PowerSpreadOption_Formula> y;
	return y(i,j,a);
}

///////////////////////////////////////////////////////////////////////
///  
///			Gaussian Copula,	Certitude  
///
///////////////////////////////////////////////////////////////////////
double Export_Gaussian_SABR_Power_SpreadOption_Certitude(
										 double copula_corr,double t,int n)
{
ArgumentList a(copula_corr);
return ARM_CF_PowerSpreadOption_Formula<SABR_smile,GaussianCopula>::Certitude(a,t,n);

}



///////////////////////////////////////////////////////////////////////
///  
///		Student Copula, Pricing Functions
///
///////////////////////////////////////////////////////////////////////


/// callput =  1  for call
/// callput = -1  for put

double Export_Student_SABR_Power_SpreadOption(double S1,double S2,
								   double alpha1, double beta1, double rho1, double nu1,
								   double alpha2, double beta2, double rho2, double nu2,
								   double copula_corr,double copula_degre,double t,int flag,
								   double a10,double b10,double k10,double a20,double b20,double k20,int n,double alpha_exp,double alpha_tanh,double kb_tanh)
{
	ArgumentList a( S1, S2,
					alpha1,  beta1,  rho1,  nu1,
					alpha2,  beta2,  rho2,  nu2,
					copula_corr,copula_degre, t,
					a10, b10, k10, a20, b20, k20, alpha_exp,alpha_tanh,kb_tanh,flag,n,GENERIC_SPREADOPTION);

	Power_Expression<ARM_CF_SABR_Student_PowerSpreadOption_Formula> y;
	return y(a);
}

double Export_Student_SABR_Power_Digital_SpreadOption(double S1,double S2,
								   double alpha1, double beta1, double rho1, double nu1,
								   double alpha2, double beta2, double rho2, double nu2,
								   double copula_corr,double copula_degre,double t,int flag,
								  double k10,double a20,double b20,double k20,int n,double alpha_exp,double alpha_tanh,double kb_tanh)
{
	ArgumentList a( S1, S2,
					alpha1,  beta1,  rho1,  nu1,
					alpha2,  beta2,  rho2,  nu2,
					copula_corr,copula_degre, t,
					0, 0, k10, a20, b20, k20,alpha_exp,alpha_tanh,kb_tanh, flag,n,DIGITAL_SPREADOPTION);

	Power_Expression<ARM_CF_SABR_Student_PowerSpreadOption_Formula> y;
	return y(a);
}


///////////////////////////////////////////////////////////////////////
///  
///		Gaussian Copulas, uses the direct copula approach
///
///////////////////////////////////////////////////////////////////////

double Export_GaussianSABRDigitalCall(double f1,
							   double alpha1,
							   double beta1,
							   double rho1,
							   double nu1,
							   int SABRFlag1,
							   double f2,
							   double alpha2,
							   double beta2,
							   double rho2,
							   double nu2,
							   int SABRFlag2,
							   double rho,
							   double K,
							   double T,
							   int legendreNb ,double alpha_exp,double alpha_tanh,double kb_tanh
							   )
{
	return GaussianSABRDigitalCall( f1,
							    alpha1,
							    beta1,
							    rho1,
							    nu1,
							    SABRFlag1,
							    f2,
							    alpha2,
							    beta2,
							    rho2,
							    nu2,
							    SABRFlag2,
							    rho,
							    K,
							    T,
							    legendreNb, alpha_exp, alpha_tanh, kb_tanh
							   );
}


double Export_GaussianSABRDigitalCallPayingS1(double f1,
							   double alpha1,
							   double beta1,
							   double rho1,
							   double nu1,
							   int SABRFlag1,
							   double f2,
							   double alpha2,
							   double beta2,
							   double rho2,
							   double nu2,
							   int SABRFlag2,
							   double rho,
							   double K,
							   double T,
							   int legendreNb,double alpha_exp,double alpha_tanh,double kb_tanh
							   )
{		
	return GaussianSABRDigitalCallPayingS1( f1,
							    alpha1,
							    beta1,
							    rho1,
							    nu1,
							    SABRFlag1,
							    f2,
							    alpha2,
							    beta2,
							    rho2,
							    nu2,
							    SABRFlag2,
							    rho,
							    K,
							    T,
							    legendreNb, alpha_exp, alpha_tanh, kb_tanh
							   );
}





double Export_GaussianSABRDigitalCallPayingS2(double f1,
							   double alpha1,
							   double beta1,
							   double rho1,
							   double nu1,
							   int SABRFlag1,
							   double f2,
							   double alpha2,
							   double beta2,
							   double rho2,
							   double nu2,
							   int SABRFlag2,
							   double rho,
							   double K,
							   double T,
							   int legendreNb,double alpha_exp,double alpha_tanh,double kb_tanh
							   )
{
	return GaussianSABRDigitalCallPayingS2( f1,
							    alpha1,
							    beta1,
							    rho1,
							    nu1,
							    SABRFlag1,
							    f2,
							    alpha2,
							    beta2,
							    rho2,
							    nu2,
							    SABRFlag2,
							    rho,
							    K,
							    T,
							    legendreNb,alpha_exp,alpha_tanh,kb_tanh
							   );
}

double Export_GaussianSABRDigitalCallPayingS3(double f1,
							   double alpha1,
							   double beta1,
							   double rho1,
							   double nu1,
							   int SABRFlag1,
							   double f2,
							   double alpha2,
							   double beta2,
							   double rho2,
							   double nu2,
							   int SABRFlag2,
							   double f3,
							   double sigma3,
							   double rho12,
							   double rho13,
							   double rho23,
							   double K,
							   double T,
							   int legendreNb,double alpha_exp,double alpha_tanh,double kb_tanh
							   )
{
	return GaussianSABRDigitalCallPayingS3( f1,
							    alpha1,
							    beta1,
							    rho1,
							    nu1,
							    SABRFlag1,
							    f2,
							    alpha2,
							    beta2,
							    rho2,
							    nu2,
							    SABRFlag2,
							    f3,
							    sigma3,
							    rho12,
							    rho13,
							    rho23,
							    K,
							    T,
							    legendreNb,alpha_exp,alpha_tanh,kb_tanh
							   );
}



CC_END_NAMESPACE()
 


/*---------------------------------------------------------------------------*/
/*---- End of file ----*/