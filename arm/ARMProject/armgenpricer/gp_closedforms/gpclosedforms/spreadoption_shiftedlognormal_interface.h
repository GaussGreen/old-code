/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 * Version initiale 01/15/2004
 *
 *  basic functions for the closed form framework 
 *
 *	\file spreadoption_shiftedlognormal.h
 *
 *  \brief
 *
 *	\author  O. Croissant
 *	\version 1.0
 *	\date January 2004
 */
 
#ifndef _GP_CF_SPREADOPTION_SHIFTEDLOGNORMAL_INTERFACE_H
#define _GP_CF_SPREADOPTION_SHIFTEDLOGNORMAL_INTERFACE_H

#include "firsttoinc.h"
#include "gpbase/port.h"


CC_BEGIN_NAMESPACE(ARM)



//////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////

///  
///			Begining Exportable  Pricing Functions 
///
//////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////
///  
///			Value  
///
///////////////////////////////////////////////////////////////////////


/// callput =  1 (K_CALL) for call
/// callput = -1  (K_PUT) for put

double Export_Gaussian_ShiftedLN_Power_SpreadOption(double S1,double S2,
								   double sigma1, double alpha1, 
								   double sigma2, double alpha2, 
								   double copula_corr,double t,
								   double a10,double b10,double k10,double a20,double b20,double k20,int n=64);

///////////////////////////////////////////////////////////////////////
///  
///			1st Derivatives  
///
///////////////////////////////////////////////////////////////////////
double Export_Gaussian_ShiftedLN_Power_SpreadOption(int i,
								   double S1,double S2,
								   double sigma1, double alpha1, 
								   double sigma2, double alpha2, 
								   double copula_corr,double t,
								   double a10,double b10,double k10,double a20,double b20,double k20,int n=64);

///////////////////////////////////////////////////////////////////////
///  
///			2nd Derivatives  
///
///////////////////////////////////////////////////////////////////////
double Export_Gaussian_ShiftedLN_Power_SpreadOption(int i,int j,
								   double S1,double S2,
								   double sigma1, double alpha1, 
								   double sigma2, double alpha2, 
								   double copula_corr,double t,
								   double a10,double b10,double k10,double a20,double b20,double k20,int n=64);

///////////////////////////////////////////////////////////////////////
///  
///			Certitude  
///
///////////////////////////////////////////////////////////////////////

double Export_Gaussian_ShiftedLN_Power_SpreadOption_Certitude(
										 double copula_corr,double t,int n);

///////////////////////////////////////////////////////////////////////
///  
///			End Exportable Pricing Functions
///
///////////////////////////////////////////////////////////////////////

#undef ARM_CF_SQRT_2_PI 
#undef ARM_CF_INVSQRTPI 

CC_END_NAMESPACE()



#endif

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/

