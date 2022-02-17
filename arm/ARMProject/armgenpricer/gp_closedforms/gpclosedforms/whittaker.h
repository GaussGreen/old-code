/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 * Version initiale 01/15/2004
 *
 *  basic functions for the closed form framework 
 *
 *	\file whittaker.h
 *
 *  \brief
 *
 *	\author  O. Croissant
 *	\version 1.0
 *	\date January 2004
 */
 
#ifndef _GP_CF_WHITTAKER_H
#define _GP_CF_WHITTAKER_H

#include "gpbase/port.h"
#include "firsttoinc.h"

#include <complex>
#include "long_double.h"

///
/// Implementation of the whittaker functions (generalization of the bessel functions, an avatar of the confluent hypergeometric functions) 
/// 
CC_USING_NS(std,complex)

CC_BEGIN_NAMESPACE(ARM)

complex<double> GammaLog(complex<double> z);
complex<double> Gamma(complex<double> z);

complex<double>  Hypergeometric_Whittaker_M ( complex<double>  a, complex<double>  b, complex<double>  z);
complex<double>  Hypergeometric_Whittaker_W ( complex<double>  a, complex<double>  b, complex<double>  z);

double  Hypergeometric_Whittaker_M (double  a,double  b,double  z);
double  Hypergeometric_Whittaker_W (double  a,double  b,double  z);

/// the following function returns structure oftype long_double which are meant to transport numbers with huge exponents...
/// this is because the C++ standard long double do not work
long_double  Hypergeometric_Whittaker_M_L (double  a,double  b,double  z);
long_double  Hypergeometric_Whittaker_W_L (double  a,double  b,double  z);

long double exp_long(double x);

CC_END_NAMESPACE()

#endif

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/

