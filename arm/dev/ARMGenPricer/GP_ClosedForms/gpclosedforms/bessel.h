/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 * Version initiale 01/15/2004
 *
 *  basic functions for the closed form framework 
 *
 *	\file bessel.h
 *
 *  \brief
 *
 *	\author  O. Croissant
 *	\version 1.0
 *	\date January 2004
 */
 
#ifndef _GP_CF_BESSEL_H
#define _GP_CF_BESSEL_H

#include <glob/firsttoinc.h>
#include "gpbase/port.h"

CC_BEGIN_NAMESPACE(ARM)

/// Implementation of the subroutine given in Numerical Recipes  in C++ for the computation ofthe Bessel Function J and Y and I and K 
/// for integer dimension and non integer (fractional) dimensions
/// the structure of the subroutine have been made a little bit more transparent and natural

double bessel_integer_J(const int n, const double x);
double bessel_integer_Y(const int n, const double x);
double bessel_integer_I(const int n, const double x);
double bessel_integer_K(const int n, const double x);
double bessel_fractional_I(double nu, double x);
double bessel_fractional_K(double nu, double x);
double bessel_fractional_I_Derx(double nu, double x);	// derivative with respect ot x
double bessel_fractional_K_Derx(double nu, double x);	// derivative with respect ot x
double bessel_fractional_J(double nu, double x);
double bessel_fractional_Y(double nu, double x);
double bessel_fractional_J_Derx(double nu, double x);	// derivative with respect ot x
double bessel_fractional_Y_Derx(double nu, double x);	// derivative with respect ot x
double bessel_fractional_I_LargeZ_Log(double nu,double z, int n);
double Hankel_Symbol(double alpha,double n);

CC_END_NAMESPACE()



#endif

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/

