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
 
#ifndef _GP_CF_TRISPREADOPTION_LOGNORMAL_INTERFACE_H
#define _GP_CF_TRISPREADOPTION_LOGNORMAL_INTERFACE_H

#include <glob/firsttoinc.h>
#include "gpbase/port.h"



CC_BEGIN_NAMESPACE(ARM)

/////////////////////////////////////////////////////////////////////////////////////////////////////
///
///   Forward value Tri Spread Option
///
///
/////////////////////////////////////////////////////////////////////////////////////////////////////


double Export_LogNormal_TriSpreadOption(double S1,double S2,double S3,double sig1,double sig2,double sig3,
										double rho12,double rho13,double rho23,double mu1,double mu2,double mu3,
										double a0,double a1,double a2,double a3,double t,int n);




///  
///			1st Derivatives  Tri Spread Option
///

double Export_LogNormal_TriSpreadOption(int i,double S1,double S2,double S3,double sig1,double sig2,double sig3,
										double rho12,double rho13,double rho23,double mu1,double mu2,double mu3,
										double a0,double a1,double a2,double a3,double t,int n);



///  
///			2nd Derivatives  Tri Spread Option
///

double Export_LogNormal_TriSpreadOption(int i,int j,double S1,double S2,double S3,double sig1,double sig2,double sig3,
										double rho12,double rho13,double rho23,double mu1,double mu2,double mu3,double a0,
										double a1,double a2,double a3,double t,int n);

/////////////////////////////////////////////////////////////////////////////////////////////////////
///
///   Forward value  Digital Tri Spread Digital option
///
///
/////////////////////////////////////////////////////////////////////////////////////////////////////


double Export_LogNormal_TriSpreadDigitalOption(double S1,double S2,double S3,double sig1,double sig2,double sig3,
											   double rho12,double rho13,double rho23,double mu1,double mu2,double mu3,
											   double a0,double a1,double a2,double a3,double t,int calput,int n);





///  
///			1st Derivatives  Digital Tri Spread Digital option
///

double Export_LogNormal_TriSpreadDigitalOption(int i,double S1,double S2,double S3,double sig1,double sig2,double sig3,
											   double rho12,double rho13,double rho23,double mu1,double mu2,double mu3,
											   double a0,double a1,double a2,double a3,double t,int calput,int n);



///  
///			2nd Derivatives  Digital Tri Spread Digital option
///

double Export_LogNormal_TriSpreadDigitalOption(int i,int j,double S1,double S2,double S3,double sig1,double sig2,double sig3,
											   double rho12,double rho13,double rho23,double mu1,double mu2,double mu3,
											   double a0,double a1,double a2,double a3,double t,int calput,int n);


/////////////////////////////////////////////////////////////////////////////////////////////////////
///
///   Forward value  Digital Tri Spread Digital option2 
///		which pays  (B0+B1*S1)^+ if A0+A2*S2+A3*S3>0
///
///
/////////////////////////////////////////////////////////////////////////////////////////////////////


double Export_LogNormal_TriSpreadDigitalOption2(double S1,double S2,double S3,double sig1,double sig2,double sig3,
											   double rho12,double rho13,double rho23,double mu1,double mu2,double mu3,
											   double a0,double a2,double a3,double b0,double b1,double t,int n);

///  
///			1st Derivatives  Digital Tri Spread Digital option
///

double Export_LogNormal_TriSpreadDigitalOption2(int i,double S1,double S2,double S3,double sig1,double sig2,double sig3,
											   double rho12,double rho13,double rho23,double mu1,double mu2,double mu3,
											   double a0,double a2,double a3,double b0,double b1,double t,int n);


///  
///			2nd Derivatives  Digital Tri Spread Digital option
///

double Export_LogNormal_TriSpreadDigitalOption2(int i,int j,double S1,double S2,double S3,double sig1,double sig2,double sig3,
											   double rho12,double rho13,double rho23,double mu1,double mu2,double mu3,
											   double a0,double a2,double a3,double b0,double b1,double t,int n);



CC_END_NAMESPACE()



#endif

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/

