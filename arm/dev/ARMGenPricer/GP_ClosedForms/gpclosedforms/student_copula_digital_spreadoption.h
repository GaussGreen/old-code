/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 * Version initiale 01/15/2004
 *
 *  basic functions for the closed form framework 
 *
 *	\file student_copula_digital_spreadoption.h
 *
 *  \brief
 *
 *	\author  O. Croissant
 *	\version 1.0
 *	\date January 2004
 */
 
#ifndef _GP_CF_STUDENTCOPULA_DIGITAL_SPREADOPTION_H
#define _GP_CF_STUDENTCOPULA_DIGITAL_SPREADOPTION_H

#include <glob/firsttoinc.h>
#include "gpbase/port.h"

#include <math.h>

#include "gpclosedforms/basic_distributions.h"
#include "gpclosedforms/generic_copula.h"


CC_BEGIN_NAMESPACE(ARM)
class GaussLegendre_Coefficients;

///////////////////////////////////////////////////////////////////////////////////////
///
///	Class  : Student Copula : two args  : correlation and degré
/// 
///////////////////////////////////////////////////////////////////////////////////////

double StudentIntegralQ(double nu,double H,double x);

double StudentIntegralQT(double nu,double x);

double Student_Copula_Digital_Spreadoption(double rho,double nu,
										   double l1a,double l2a,double l3a,double l4a,double l5a,double l6a,
										   double l1b,double l2b,double l3b,double l4b,double l5b,double l6b,
										   double K,
										   GaussLegendre_Coefficients* coefs);



CC_END_NAMESPACE()



#endif

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/

