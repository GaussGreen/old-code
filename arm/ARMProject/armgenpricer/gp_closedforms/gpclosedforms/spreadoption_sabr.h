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
 
#ifndef _GP_CF_SPREADOPTION_SABR_H
#define _GP_CF_SPREADOPTION_SABR_H

////////////////////////////////////////////////////////////////////////////////////

/// templated approach:

////////////////////////////////////////////////////////////////////////////////////

/// we use generic pricing mechanism based on structure::PowerSpreadOption_Pricing and 
/// PowerSpreadOption_Pricing_With_Limits for the case where we can compute the starting point of integration
/// called DistributionZaInverseLimit<smile,copula>


////////////////////////////////////////////////////////////////////////////////////

/// Direct approach from a gaussian copula perspective: 

////////////////////////////////////////////////////////////////////////////////////
#include "firsttoinc.h"
#include "gpbase/port.h"



CC_BEGIN_NAMESPACE(ARM)


double GaussianSABRDigitalCall(double f1,
							   double alpha1,
							   double beta1,
							   double rho1,
							   double nu1,
							   int		SABRFlag1,
							   double f2,
							   double alpha2,
							   double beta2,
							   double rho2,
							   double nu2,
							   int		SABRFlag2,
							   double rho,
							   double K,
							   double T,
							   int legendreNb,double alpha_exp,double alpha_tanh,double kb_tanh
							   );

double GaussianSABRDigitalCallPayingS1(double f1,
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

double GaussianSABRDigitalCallPayingS2(double f1,
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

double GaussianSABRDigitalCallPayingS3(double f1,
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

