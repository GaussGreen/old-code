/*!
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 * Version initiale 01/15/2004
 *
 *  basic functions for the closed form framework 
 *
 *	\file student_copula.cpp
 *
 *  \brief
 *
 *	\author  O. Croissant
 *	\version 1.0
 *	\date January 2004
 */
#include "firsttoinc.h"
#include "gpbase/port.h"

#include <cmath>

#include "gpclosedforms/student_copula.h"
#include "gpbase/numericconstant.h"
#include "gpclosedforms/incompletebeta.h"
#include "gpclosedforms/gamma.h"
#include "gpclosedforms/student_copula_digital_spreadoption.h"
#include "gpclosedforms/hypergeometric.h"
#include "gpclosedforms/gaussian_integrals.h"

#include "expt.h"

using namespace std;
CC_BEGIN_NAMESPACE(ARM)



///////////////////////////////////////////////////////////////////////////////////////
///
/// 
/// here the copula is the the student copula

/// take gaussian points and add the copula dependency to return the density
///			Contens of the argumentlist for the gaussian copula:
///
///		Underlying[0],	//Correlation
///		Underlying[1],	//Degre
///
///////////////////////////////////////////////////////////////////////////////////////

double StudentIntegralQ(double nu,double H,double x)
{
	return H*pow(nu/(nu+x*x),1.+nu/2.)*Hypergeometric2F1(0.5,1.+nu/2.,1.5,-H*H/(x*x+nu))/(2.*ARM_NumericConstants::ARM_PI);
}


double StudentIntegralQT(double nu,double x)
{
	return pow(nu,nu/2.)*pow(x*x+nu,-(1.+nu)/2.)*gamma((1.+nu)/2.)/(2.*gamma(nu/2.))/ARM_NumericConstants::ARM_SQRT_PI;
}


double Student_Copula_Digital_Spreadoption(double rho,double nu,
										   double l1a,double l2a,double l3a,double l4a,double l5a,double l6a,
										   double l1b,double l2b,double l3b,double l4b,double l5b,double l6b,
										   double K,
										   GaussLegendre_Coefficients* coefs)
{
	return 0;
}


CC_END_NAMESPACE()
 


/*---------------------------------------------------------------------------*/
/*---- End of file ----*/
