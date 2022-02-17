/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 * Version initiale 01/15/2004
 *
 *  basic functions for the closed form framework 
 *
 *	\file erf.h
 *
 *  \brief
 *
 *	\author  O. Croissant
 *	\version 1.0
 *	\date January 2004
 */
 
#ifndef _GP_CF_ERF_H
#define _GP_CF_ERF_H

#include <glob/firsttoinc.h>
#include "gpbase/port.h"

#include <complex>
CC_USING_NS(std,complex)

CC_BEGIN_NAMESPACE(ARM)

////////////////////////////////////////////////////////////////////////////////////////////
///
///				computation of real error function erf(x) 
///					erf(t) =   2/Sqrt[Pi]  Integrate[Exp[-x^2], {x, 0, t}]
///
////////////////////////////////////////////////////////////////////////////////////////////

double erf(double x);
 
////////////////////////////////////////////////////////////////////////////////////////////
///
///				computation of imaginary error function erfi(x) 
///					erfi(x) = erf(i x) / i  =   2/Sqrt[Pi]  Integrate[Exp[x^2], {x, 0, t}]
///
////////////////////////////////////////////////////////////////////////////////////////////

double erfi(const double x);

////////////////////////////////////////////////////////////////////////////////////////////
///
///				computation of complex error function erf(z) 
///					
////////////////////////////////////////////////////////////////////////////////////////////

complex<double> cerf(complex<double> z,int nbsteps);

complex<double> cnormal(complex<double> z,int nbsteps);


////////////////////////////////////////////////////////////////////////////////////////////
///
///			Exportation of	computation of complex error function erf(z) 
///					
////////////////////////////////////////////////////////////////////////////////////////////

double Export_RealPart_ComplexErf(double realpart, double imaginarypart, int nbterm);
double Export_ImaginaryPart_ComplexErf(double realpart, double imaginarypart, int nbterm);


CC_END_NAMESPACE()



#endif

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/

