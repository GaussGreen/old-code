/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 * Version initiale 01/15/2004
 *
 *  basic functions for the closed form framework 
 *
 *	\file normal.h
 *
 *  \brief
 *
 *	\author  O. Croissant
 *	\version 1.0
 *	\date January 2004
 */
 
#ifndef _GP_CF_NORMAL_H
#define _GP_CF_NORMAL_H

#include "firsttoinc.h"
#include "gpbase/port.h"


CC_BEGIN_NAMESPACE(ARM)

double NormalCDF(double x);

double NormalCDF(double x, double y, double rho,
				 double nGaussLegendre_1=10.0, double nGaussLegendre_2=4.0, double lowerBound=-8.0, double upperBound=8.0);

double Normal_X_Expectation(double X, double y, double rho,
				 double nGaussLegendre_1=10.0, double nGaussLegendre_2=4.0, double lowerBound=-8.0, double upperBound=8.0);

double Normal_Y_Expectation(double x, double Y, double rho,
				 double nGaussLegendre_1=10.0, double nGaussLegendre_2=4.0, double lowerBound=-8.0, double upperBound=8.0);

double Normal_XY_Expectation(double x, double Y, double rho,
				 double nGaussLegendre_1=10.0, double nGaussLegendre_2=4.0, double lowerBound=-8.0, double upperBound=8.0);

double Normal_XX_Expectation(double x, double Y, double rho,
				 double nGaussLegendre_1=10.0, double nGaussLegendre_2=4.0, double lowerBound=-8.0, double upperBound=8.0);

double Normal_YY_Expectation(double x, double Y, double rho,
				 double nGaussLegendre_1=10.0, double nGaussLegendre_2=4.0, double lowerBound=-8.0, double upperBound=8.0);


double NormalCDF2(double x);

double NormalCDF3(double x);

double NormalPDF(double x);

double NormalCDFInverse(double x);

double RegularNormalCDF2(double d);

/// Integral of the bivariate with respect to y (derivative of the biavriate cummulative w.r. to x)
double DxBivariateCummulative2D(double x , double y, double rho);

/// Integral of the bivariate with respect to x (derivative of the biavriate cummulative w.r. to y)
double DyBivariateCummulative2D(double x , double y, double rho);


		/// Derivative of the gaussian 2-copula with respect to x

double DxBivariateCopula(double x,double y, double rho);

		/// Derivative of the gaussian 2-copula with respect to y

double DyBivariateCopula(double x,double y, double rho);




CC_END_NAMESPACE()



#endif

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/