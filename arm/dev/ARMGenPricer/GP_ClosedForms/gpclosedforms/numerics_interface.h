/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 * Version initiale 01/15/2004
 *
 *  basic functions for the closed form framework 
 *
 *	\file numerics_interface.h
 *
 *  \brief
 *
 *	\author  O. Croissant
 *	\version 1.0
 *	\date January 2004
 */
///////////////////////////////////////////////////////////////////////////////
///
///					Process :
///			dS= S(rdt+V^(1/2) dW1+ Jdq
///			dV=(omega-theta*V)dt +ksi*V^(1/2) dW2 
///			dW1.dW2=rho*dt
///
///			Jump J : probability lambda, volatility sigmaJ, log-size muJ
///
/////////////////////////////////////////////////////////////////////////////:

 
#ifndef _GP_CF_NUMERICS_INTERFACE_H
#define _GP_CF_NUMERICS_INTERFACE_H


#include <glob/firsttoinc.h>
#include "gpbase/port.h"
#include "gpclosedforms/gaussian_integrals.h"
//#include "gpclosedforms/sabr_calibration.h"
//#include "heston_calibration.h"
#include "gpclosedforms/eigenvalues.h"
#include "gpbase/gpvector.h"
#include "gpbase/typedef.h"



CC_BEGIN_NAMESPACE(ARM)

class GeneralizedHeston_ParameterSet;
class SABR_ParameterSet;

/////////////////////////////////////////////////////////////////////////////////////////////////////
///  
///			Begining Exportable  Pricing Functions 
///
/////////////////////////////////////////////////////////////////////////////////////////////////////
double Export_Bivariate(double x, double y, double rho, int p = 0, int q = 0);

double Export_Gamma(double a);

double Export_Lambert_Function(double a);

double Export_IncompleteBeta(const double& a, const double& b, const double& x);

double Export_InverseIncompleteBeta(const double& a, const double& b, const double& x0, const double& x);

double  Export_Hypergeometric2F1( double  a, double  b, double  c, double  z);

GaussLegendre_Coefficients* Export_GaussLegendre_Coefficients(double a,double b,int n);

GaussHermite_Coefficients* Export_GaussHermite_Coefficients(int n);

GaussLaguerre_Coefficients* Export_GaussLaguerre_Coefficients(double a,int n);


GeneralizedHeston_ParameterSet*   Export_GeneralizedHeston_Model_Calibrate(ARM_GP_Vector* K_Vec,ARM_GP_Vector* ImpVol_Vec,double f,double t,int nbsteps,int algorithm,
												double V0,double omega,double theta,double ksi,double rho,double muJ,double sigmaJ,double lambda);

GeneralizedHeston_ParameterSet*   Export_GeneralizedHeston_Model_Calibrate_NoJump(ARM_GP_Vector* K_Vec,ARM_GP_Vector* ImpVol_Vec,double f,double t,double muJ0,double sigmaJ0,double lambda0,int nbsteps,int algorithm,
												double V0,double omega,double theta,double ksi,double rho);

SABR_ParameterSet*  Export_SABR_Model_Calibrate(ARM_GP_Vector* K_Vec,ARM_GP_Vector* ImpVol_Vec,double f,double beta,double tex,int flag,int nbsteps,int algorithm,
												double alpha0,double rho0,double nu0,double alphap,double rhop, double nup,double rweight_alpha,double rweight_rho,double rweight_nu);


double Export_ImcompleteBeta_Inverse(double a, double b, double x);

double Export_Student_QIntegral(double a, double b, double x);

double  Export_Hypergeometric_Whittaker_W (double  a,double  b,double  z);

double  Export_Hypergeometric_Whittaker_M (double  a,double  b,double  z);

double Export_Bessel_I(double nu, double x);

double Export_Bessel_J(double nu, double x);

double Export_Bessel_K(double nu, double x);

double Export_Bessel_Y(double nu, double x);

double Export_Hypergeometric_Appell(double a,double b1,double b2,double c ,double x,double xim, double y,double yim,int nb);

void Export_Eigenvalues4(double rho12,double rho13, double rho14,double rho23,double rho24,double rho34,double* e1,double* e2,double* e3,double* e4);

void Export_Eigenvalues3(double rho12,double rho13, double rho23,double* e1,double* e2,double* e3);



CC_END_NAMESPACE()



#endif

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/
