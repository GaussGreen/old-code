/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 * Version initiale 01/15/2004
 *
 *  basic functions for the closed form framework 
 *
 *	\file heston.h
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

 
#ifndef _GP_CF_NORMAL_HESTON_H
#define _GP_CF_NORMAL_HESTON_H


#include "firsttoinc.h"
#include "gpbase/port.h"

#include <complex>
using std::complex;

CC_BEGIN_NAMESPACE(ARM)


complex<long double> GaussianHestonFondamentalTransform(double rho,double lambdaV,double thetaV,
						double kappaV,double V0, complex<long double> omega, double T);

complex<long double> FourierPayoff(complex<long double> omega, complex<long double> k, double lambdaB);


double NormalHeston(double rho,double lambdaV,double thetaV,
					double kappaV,double V0,  double S0,double k,double T,int callput,double lambdaB = 0.1, 
					int nbfirst = 40,int nb = 40,
					int NbStage = -1, int NbOscill = 0,double prec = 1e-8);

double Complexified_NormalHeston(complex<double> rho,complex<double> lambdaV,complex<double> thetaV,
						complex<double> kappaV,complex<double> V0,  complex<double> S0,complex<double> k,complex<double> T,
						int callput,double lambdaB,int nbfirst,int nb,
						int NbStage, int NbOscill,double prec);

double NormalHestonTry(double rho,double lambdaV,double thetaV,
						double kappaV,double V0,  double S0,double k,double T,int callput,int nb = 20, double lamdaB = 0.1);

complex<long double> SuperNormalHestonFondamentalTransform(double rho1,double lambda1,double theta1,double kappa1,double V01,
														   double rho2,double lambda2,double theta2,double kappa2,double V02,
														   complex<long double> omega, double T);

double SuperNormalHeston(double rho1, double lambdaV1, double thetaV1, double kappV1, double V01,
						 double rho2, double lambdaV2, double thetaV2, double kappV2, double V02,
						 double S0, double k, double T, int callput, int nb = 20);
						


CC_END_NAMESPACE()



#endif

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/
