/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 * Version initiale 01/15/2004
 *
 *  basic functions for the closed form framework 
 *
 *	\file spreadoption_lognormal.h
 *
 *  \brief
 *
 *	\author  O. Croissant
 *	\version 1.0
 *	\date January 2004
 */
 
#ifndef _GP_CF_SPREADOPTION_SABR_INTERFACE_H
#define _GP_CF_SPREADOPTION_SABR_INTERFACE_H

#include <glob/firsttoinc.h>
#include "gpbase/port.h"



CC_BEGIN_NAMESPACE(ARM)



/////////////////////////////////////////////////////////////////////////////////////////////
///  
///		Gaussian Copula,  Pricing Functions 
///
//////////////////////////////////////////////////////////////////////////////////////////////


/// callput =  1 (K_CALL) for call
/// callput = -1  (K_PUT) for put

double Export_Gaussian_SABR_Power_SpreadOption(double S1,double S2,
								   double alpha1, double beta1, double rho1, double nu1,
								   double alpha2, double beta2, double rho2, double nu2,
								   double copula_corr,double t,int flag,
								   double a10,double b10,double k10,double a20,double b20,double k20,int n,double alpha_exp,double alpha_tanh,double kb_tanh);

double Export_Gaussian_SABR_Power_Digital_SpreadOption(double S1,double S2,
								   double alpha1, double beta1, double rho1, double nu1,
								   double alpha2, double beta2, double rho2, double nu2,
								   double copula_corr,double t,int flag,
								   double k10,double a20,double b20,double k20,int n,double alpha_exp,double alpha_tanh,double kb_tanh);

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
							   int legendreNb,double alpha_exp,double alpha_tanh,double kb_tanh
							   );

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
							   );

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
							   );

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
							   );




///////////////////////////////////////////////////////////////////////
///  
///				Gaussian Copula,  1st Derivatives  
///
///////////////////////////////////////////////////////////////////////
double Export_Gaussian_SABR_Power_SpreadOption(int i,
								   double S1,double S2,
								   double alpha1, double beta1, double rho1, double nu1,
								   double alpha2, double beta2, double rho2, double nu2,
								   double copula_corr,double t,int flag,
								   double a10,double b10,double k10,double a20,double b20,double k20,int n,double alpha_exp,double alpha_tanh,double kb_tanh);

double Export_Gaussian_SABR_Power_Digital_SpreadOption(int i,
								   double S1,double S2,
								   double alpha1, double beta1, double rho1, double nu1,
								   double alpha2, double beta2, double rho2, double nu2,
								   double copula_corr,double t,int flag,
								   double k10,double a20,double b20,double k20,int n,double alpha_exp,double alpha_tanh,double kb_tanh);

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
								   double a10,double b10,double k10,double a20,double b20,double k20,int n,double alpha_exp,double alpha_tanh,double kb_tanh);

double Export_Gaussian_SABR_Power_Digital_SpreadOption(int i,int j,
								   double S1,double S2,
								   double alpha1, double beta1, double rho1, double nu1,
								   double alpha2, double beta2, double rho2, double nu2,
								   double copula_corr,double t,int flag,
								   double k10,double a20,double b20,double k20,int n,double alpha_exp,double alpha_tanh,double kb_tanh);

///////////////////////////////////////////////////////////////////////
///  
///				Gaussian Copula,  Certitude  
///
///////////////////////////////////////////////////////////////////////

double Export_Gaussian_SABR_Power_SpreadOption_Certitude(
										 double copula_corr,double t,int n);

/////////////////////////////////////////////////////////////////////////////////////////////
///  
///		Student Copula,  Pricing Functions 
///
//////////////////////////////////////////////////////////////////////////////////////////////


/// callput =  1 (K_CALL) for call
/// callput = -1  (K_PUT) for put

double Export_Student_SABR_Power_SpreadOption(double S1,double S2,
								   double alpha1, double beta1, double rho1, double nu1,
								   double alpha2, double beta2, double rho2, double nu2,
								   double copula_corr,double copula_degre,double t,int flag,
								   double a10,double b10,double k10,double a20,double b20,double k20,int n,double alpha_exp,double alpha_tanh,double kb_tanh);

double Export_Student_SABR_Power_Digital_SpreadOption(double S1,double S2,
								   double alpha1, double beta1, double rho1, double nu1,
								   double alpha2, double beta2, double rho2, double nu2,
								   double copula_corr,double copula_degre,double t,int flag,
								   double k10,double a20,double b20,double k20,int n,double alpha_exp,double alpha_tanh,double kb_tanh);




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
							   int legendreNb,double alpha_exp,double alpha_tanh,double kb_tanh
							   );


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
							   );

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
							   );

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
							   );



CC_END_NAMESPACE()



#endif

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/

