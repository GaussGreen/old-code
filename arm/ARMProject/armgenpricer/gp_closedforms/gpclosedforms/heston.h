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

 
#ifndef _GP_CF_HESTON_H
#define _GP_CF_HESTON_H


#include "firsttoinc.h"
#include "gpbase/port.h"


CC_BEGIN_NAMESPACE(ARM)

double GeneralizedHeston_ImplicitVol(double F,double K,double V0, double t,double omega,double theta,double ksi,
						 double rho, double muJ, double sigmaJ, double lambda,int nb);


double GeneralizedHeston(double F,double K,double sig, double t,double omega,double theta,
				double ksi,double rho, double muJ, double sigmaJ, double lambda,int nb);

double GeneralizedHeston2(double F,double K,double sig, double t,double omega,double theta,double ksi,
						 double rho, double muJ, double sigmaJ, double lambda,int nb);

double GeneralizedHeston2(double F,double K,double sig, double t,double lgtvol,double theta,double ksi,
						 double rho, double muJ, double sigmaJ, double lambda,
						 int nb1, int nb,int NbStage, int NbOscill,double prec);

long double GeneralizedHeston_Integrand(double F,double K,double k, double V0, double t,double omega,
				double theta,double ksi,double rho, double muJ, double sigmaJ, double lambda);

long double GeneralizedHeston_Integrand_OscillatoryComponent(double F,double K,double k, double V0, double t,
				double omega,double theta,double ksi,double rho, double muJ, double sigmaJ, double lambda);


CC_END_NAMESPACE()



#endif

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/
