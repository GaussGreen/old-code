/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 * Version initiale 01/15/2004
 *
 *  basic functions for the closed form framework 
 *
 *	\file hypergeometric.h
 *
 *  \brief
 *
 *	\author  O. Croissant
 *	\version 1.0
 *	\date January 2004
 */
 
#ifndef _GP_CF_HYPERGEOMETRIC_H
#define _GP_CF_HYPERGEOMETRIC_H

#include <glob/firsttoinc.h>
#include "gpbase/port.h"


#include "gpclosedforms/long_double.h"
#include <complex>


using std::complex;

CC_BEGIN_NAMESPACE(ARM)


complex<double>  Hypergeometric1F1 ( complex<double>  a, complex<double>  b, complex<double>  z);
complex<double>  Hypergeometric2F1 ( complex<double>  a, complex<double>  b, complex<double>  c, complex<double>  z);
complex<double>  Hypergeometric2F0 ( complex<double>  a, complex<double>  b, complex<double>  z);
complex<double>  HypergeometricU ( complex<double>  a, complex<double>  b, complex<double>  z);
complex<double>  HypergeometricAppellF1 ( complex<double>  a, complex<double>  b1,complex<double>  b2, complex<double>  c, complex<double>  x,complex<double> y, int Nb);



double  Hypergeometric2F1 ( double  a, double  b, double  c, double  z);

long double   Hypergeometric2F0_serie ( long double   a, long double   b, long double   z);
long double  Hypergeometric2F0 ( long double  a, long double  b, long double  z);
long double   HypergeometricU ( long double   a, long double   b, long double   z);



CC_END_NAMESPACE()



#endif

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/

