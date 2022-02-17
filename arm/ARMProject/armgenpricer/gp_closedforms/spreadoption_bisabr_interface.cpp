/*!
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 * Version initiale 01/15/2004
 *
 *  basic functions for the closed form framework 
 *
 *	\file spreadoption_bisabr_interface.cpp
 *
 *  \brief
 *
 *	\author  O. Croissant
 *	\version 1.0
 *	\date January 2006
 */
#include "firsttoinc.h"
#include "gpbase/port.h"

#include <cmath>
#include <complex>

#include "gpclosedforms/bisabr_spreadoption.h"
#include "gpclosedforms/bisabr_digital_spreadoption.h"
#include "gpclosedforms/basic_distributions.h"
#include "gpclosedforms/bisabr_spreadoption_formula.h"
#include "gpclosedforms/bisabr_s3_spreadoption_formula.h"
#include "gpclosedforms/spreadoption_bisabr_interface.h"

#include "expt.h"


using std::complex;

CC_BEGIN_NAMESPACE(ARM)


/*
double Export_BiSABR_SpreadOption(double F1,double alpha1,double beta1,double rho1,double nu1,
								  double F2,double alpha2,double beta2,double rho2,double nu2,
								  double K,double T,int CallPut,double rhos,double rhov,double rhoc12,double rhoc21,
								  int flag,int nbsteps)
{
	return Packaged_BiSABR_SpreadOption( F1, alpha1, beta1, rho1, nu1, F2, alpha2, beta2, rho2, nu2,
								   K, T, CallPut, rhos, rhov, rhoc12, rhoc21, flag);

}
*/

double Export_BiSABR_SpreadOption(
									 double F1,double alpha1,double beta1,double rho1,double nu1,
									 double F2,double alpha2,double beta2,double rho2,double nu2,
									 double rhos,double rhov,double rhoc12,double rhoc21,
									 double K,double T,int CallPut,int flag)
{
	ArgumentList a(
		F1, alpha1, beta1, rho1, nu1,
		F2, alpha2, beta2, rho2, nu2,
		rhos, rhov, rhoc12, rhoc21,
		K, T, CallPut, flag
		);

	Power_Expression<ARM_CF_BiSABR_SpreadOption_Formula> y;
	return y(a);
}


double Export_BiSABR_S3_SpreadOption(
									 double F1,double alpha1,double beta1,double rho1,double nu1,
									 double F2,double alpha2,double beta2,double rho2,double nu2,
									 double F3,double alpha3,double beta3,double rho3,double nu3,
									 double rhos,double rhov,double rhoc12,double rhoc21,double copularho,
									 double T,double a1,double b1,double k1,double a2,double b2,double k2,int flag,int nbsteps)
{
	ArgumentList a(
		F1, alpha1, beta1, rho1, nu1,
		F2, alpha2, beta2, rho2, nu2,
		F3, alpha3, beta3, rho3, nu3,
		rhos, rhov, rhoc12, rhoc21, copularho,
		T,a1,b1,k1,a2,b2,k2, flag, nbsteps
		);

	Power_Expression<ARM_CF_BiSABR_S3_SpreadOption_Formula> y;
	return y(a);
}







/// the following extension is mainly used for using negative F1 and/or  F2 inside the computation of digitale paying float 



double Export_Complexified_BiSABR_SpreadOption(double F1,double alpha1,double beta1,double rho1,double nu1,
								  double F2,double alpha2,double beta2,double rho2,double nu2,
					double K,double T,int callput,double rhos,double rhov,double rhoc12,double rhoc21,int flag)
{
	return Packaged_Complexified_BiSABR_SpreadOption( F1, alpha1, beta1, rho1, nu1,
								   F2, alpha2, beta2, rho2, nu2,
					 K, T, callput, rhos, rhov, rhoc12, rhoc21, flag);
}


double Export_SABR_BetaEqualZero_Option (double f,double K,double T,double mu,double alpha,double rho,double nu,int callput)
{
	return  BetaEqualZeroSABR( f, K, T, mu, alpha, rho, nu, callput);
}


CC_END_NAMESPACE()
 


/*---------------------------------------------------------------------------*/
/*---- End of file ----*/