/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 * Version initiale 01/15/2004
 *
 *  basic functions for the closed form framework 
 *
 *	\file extended_sabr.h
 *
 *  \brief
 *
 *	\author  O. Croissant
 *	\version 1.0
 *	\date January 2004
 */
 
#ifndef _GP_CF_EXTENDED_SABR_INTERFACE_H
#define _GP_CF_EXTENDED_SABR_INTERFACE_H

#include "firsttoinc.h"
#include "gpbase/port.h"

#include "gpclosedforms/sabr_calibration.h"


CC_BEGIN_NAMESPACE(ARM)


///////////////////////////////////////////////////////////////////////
///  
///			Begining Exportable  Pricing Functions 
///
///////////////////////////////////////////////////////////////////////


/// callput =  1  for call
/// callput = -1  for put

double Export_SABR_ImplicitVol(double f, 
							   double K, 
							   double tex,
							   double alpha, 
							   double beta, 
							   double rho, 
							   double nu,
							   int flag,
							   int nbsteps,
							   double C_Alpha_Exp,
							   double C_Alpha_Tanh,
							   double C_Kb_Tanh);

double Export_SABR_ImplicitVol(double f, 
							   double K, 
							   double tex,
							   double alpha, 
							   double beta, 
							   double rho, 
							   double nu,
							   int flag,
							   int nbsteps);


double Export_SABR2B_ImplicitVol(double f, 
							   double K, 
							   double tex,
							   double alpha, 
							   double beta1, 
							   double beta2, 
							   double rho, 
							   double nu,
							   double zero,
							   double lambda);


///////////////////////////////////////////////////////////////////////
///  
///			1st Derivatives  
///
///////////////////////////////////////////////////////////////////////
double Export_SABR_ImplicitVol(int i,
							   double f, 
							   double K, 
							   double tex,
							   double alpha, 
							   double beta, 
							   double rho, 
							   double nu,
							   int flag,
							   int nbsteps,
							   double C_Alpha_Exp,
							   double C_Alpha_Tanh,
							   double C_Kb_Tanh);

double Export_SABR_ImplicitVol(int i,
							   double f, 
							   double K, 
							   double tex,
							   double alpha, 
							   double beta, 
							   double rho, 
							   double nu,
							   int flag,
							   int nbsteps);


///////////////////////////////////////////////////////////////////////
///  
///			2nd Derivatives  
///
///////////////////////////////////////////////////////////////////////
double Export_SABR_ImplicitVol(int i,int j,
							   double f, 
							   double K, 
							   double tex,
							   double alpha, 
							   double beta, 
							   double rho, 
							   double nu,
							   int flag,
							   int nbsteps,
							   double C_Alpha_Exp,
							   double C_Alpha_Tanh,
							   double C_Kb_Tanh);

double Export_SABR_ImplicitVol(int i,int j,
							   double f, 
							   double K, 
							   double tex,
							   double alpha, 
							   double beta, 
							   double rho, 
							   double nu,
							   int flag,
							   int nbsteps);


///////////////////////////////////////////////////////////////////////
///  
///			Begining Exportable  Pricing Functions 
///
///////////////////////////////////////////////////////////////////////


/// callput =  1  for call
/// callput = -1  for put

double Export_SABR_VanillaOption(double f, 
							   double K, 
							   double tex,
							   double alpha, 
							   double beta, 
							   double rho, 
							   double nu,
							   int callorput,
							   int flag,
							   int nbsteps,
							   double C_Alpha_Exp,
							   double C_Alpha_Tanh,
							   double C_Kb_Tanh);

double Export_SABR_VanillaOption(double f, 
							   double K, 
							   double tex,
							   double alpha, 
							   double beta, 
							   double rho, 
							   double nu,
							   int callorput,
							   int flag,
							   int nbsteps);

double Export_SABR_FromGaussianToDistribution(double f, 
							   double K, 
							   double tex,
							   double alpha, 
							   double beta, 
							   double rho, 
							   double nu,
							   int flag,
							   int nbsteps,double st1,double st2,double st3);

double Export_SABR_FromDistributionToGaussian(double f, 
							   double K, 
							   double tex,
							   double alpha, 
							   double beta, 
							   double rho, 
							   double nu,
							   int flag,
							   int nbsteps,double st1,double st2,double st3);

///////////////////////////////////////////////////////////////////////
///  
///			1st Derivatives  
///
///////////////////////////////////////////////////////////////////////
double Export_SABR_VanillaOption(int i,
							   double f, 
							   double K, 
							   double tex,
							   double alpha, 
							   double beta, 
							   double rho, 
							   double nu,
							   int callorput,
							   int flag,
							   int nbsteps,
							   double C_Alpha_Exp,
							   double C_Alpha_Tanh,
							   double C_Kb_Tanh);

double Export_SABR_VanillaOption(int i,
							   double f, 
							   double K, 
							   double tex,
							   double alpha, 
							   double beta, 
							   double rho, 
							   double nu,
							   int callorput,
							   int flag,
							   int nbsteps);


///////////////////////////////////////////////////////////////////////
///  
///			2nd Derivatives  
///
///////////////////////////////////////////////////////////////////////
double Export_SABR_VanillaOption(int i,int j,
							   double f, 
							   double K, 
							   double tex,
							   double alpha, 
							   double beta, 
							   double rho, 
							   double nu,
							   int callorput,
							   int flag,
							   int nbsteps,
							   double C_Alpha_Exp,
							   double C_Alpha_Tanh,
							   double C_Kb_Tanh);

double Export_SABR_VanillaOption(int i,int j,
							   double f, 
							   double K, 
							   double tex,
							   double alpha, 
							   double beta, 
							   double rho, 
							   double nu,
							   int callorput,
							   int flag,
							   int nbsteps);



////////////////////////////////////////////////////////////////////////////////////////////
///
///				computation of sum{0,t} of exp(-a/x+b*x)/sqrt(x)
///
////////////////////////////////////////////////////////////////////////////////////////////


double SABR_expintegral(double a,double b, double t, int nbsteps);

////////////////////////////////////////////////////////////////////////////////////////////
///
///				Calibration of SABR
///
////////////////////////////////////////////////////////////////////////////////////////////
SABR_ParameterSet*  Export_SABR_Model_Calibrate(std::vector<double>* K_Vec,std::vector<double>* ImpVol_Vec,std::vector<double>* Weigth_Vec,double f,double tex,int flag,int nbsteps,int algorithm,
												double alpha0,double beta0,double rho0,double nu0);/// 12 params

SABR_ParameterSet*  Export_SABR_Model_Calibrate(std::vector<double>* K_Vec,std::vector<double>* ImpVol_Vec,std::vector<double>* Weigth_Vec,double tex,int flag,int nbsteps,int algorithm,
												double alpha0,double beta0,double rho0,double nu0,double f0); /// 12 params 


///  f and beta fixed
SABR_ParameterSet*  Export_SABR_Model_Calibrate(std::vector<double>* K_Vec,std::vector<double>* ImpVol_Vec,std::vector<double>* Weigth_Vec,double f,double beta,double tex,
												int flag,double alpha_exp,double alpha_tanh,double kb_tanh,int nbsteps,int algorithm,
												double alpha0,double rho0,double nu0);


SABR_ParameterSet*  Export_SABR_Model_Calibrate(std::vector<double>* K_Vec,std::vector<double>* ImpVol_Vec,std::vector<double>* Weigth_Vec,double f,double beta,double tex,
												int flag,int nbsteps,int algorithm,
												double alpha0,double rho0,double nu0,double alphap,double rhop, double nup, double rweight_alpha,double rweight_rho,double rweight_nu);  /// 18 params





///////////////////////////////////////////////////////////////////////
///  
///			conversion from sigma(K) to alpha  
///
///////////////////////////////////////////////////////////////////////
double Export_SABR_From_Sigma_To_Alpha (
							   double f, 
							   double K, 
							   double tex,
							   double sigma, 
							   double beta, 
							   double rho, 
							   double nu,
							   int flag,
							   int nbsteps);

CC_END_NAMESPACE()

#endif

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/

