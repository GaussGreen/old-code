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
 
#ifndef _GP_CF_BI_SPREADOPTION_LOGNORMAL_H
#define _GP_CF_BI_SPREADOPTION_LOGNORMAL_H

 
#include "firsttoinc.h"
#include "gpbase/port.h"


CC_BEGIN_NAMESPACE(ARM)




//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////
///
/// Integral{phi(x)*phi(y)*1{a1*exp(d*y)+a2*exp(b2*y)+a3*exp(b3*y)>0}*1{k+l*exp(m x)>0}
///
//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////


double BilognormalIntegralType1_compute(double d,double a1,double a2,double b2,
								   double a3,double b3,
								   double k,double l,double m,
								  GaussLegendre_Coefficients* glcoeffs_ptr);

/////////////////////////////////////////////////////////////////////////////////
///
///              implementation of the generalized digital: 
///                  payoff =1	only if 	a0 + a1*S1^g1+a2*S2^g2 >0
///											and b0+b1*S1^gb>0
///
/////////////////////////////////////////////////////////////////////////////////

double Generalized_BiSpreadDigitalCall_aux(double l1,double l2,
							  double m1,double m2,
							  double sig1,double sig2,
							  double r12,
							  double alpha0,double alpha1,double alpha2,
							  double beta0,double beta1,double T,
							  double g1,double g2,double gb,
							  GaussLegendre_Coefficients* glcoeffs_ptr);

///////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////
///
///
///  double digital option where the condition are : 
///			a0 + a1*S1+a2*S2 >0
///		and b0+b1*S1>0
///
///
//////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////



double BiSpreadDigitalOption(double l1,double l2,
								   double v1,double v2,
								   double m1,double m2,
								   double r12,
								   double a0,double a1,double a2,
								   double b0,double b1,
								   double T,int callput,int n);


CC_END_NAMESPACE()



#endif

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/

