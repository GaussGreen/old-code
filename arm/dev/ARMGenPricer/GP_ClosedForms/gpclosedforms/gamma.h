/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 * Version initiale 01/15/2004
 *
 *  basic functions for the closed form framework 
 *
 *	\file gamma.h
 *
 *  \brief
 *
 *	\author  O. Croissant
 *	\version 1.0
 *	\date January 2004
 */
 
#ifndef _GP_CF_GAMMA_H
#define _GP_CF_GAMMA_H

#include <glob/firsttoinc.h>
#include "gpbase/port.h"

#include "long_double.h"
#include <complex>

using std::complex;
//using std::complex<double>;

CC_BEGIN_NAMESPACE(ARM)


/////////////////////////////////////////////////////////////////////////
///
///    Gamma(x)  ( Gamma(x) = Factorial(x-1) for x integer )
///		return a long_double to handle large numbers
///  the factorial and gamma functions are here for small numbers, gammaLD handling large numbers
///
/////////////////////////////////////////////////////////////////////////


long_double  gammaLD(const double& xx);

double gamma(const double& xx);

double gammalog(const double& xx);

long double factorial(int n);

complex<double> GammaLog(complex<double> z);

complex<double> Gamma(complex<double> z);

void gser(double *gamser, double a, double x, double *gln);

void gcf(double *gammcf, double a, double x, double *gln);

double gammp(double a, double x);

CC_END_NAMESPACE()



#endif

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/

