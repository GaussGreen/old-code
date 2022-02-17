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
 
#ifndef _GP_CF_SPREADOPTION_LOGNORMAL_H
#define _GP_CF_SPREADOPTION_LOGNORMAL_H

#include <glob/firsttoinc.h>
#include "gpbase/port.h"



CC_BEGIN_NAMESPACE(ARM)


double SpreadDigitalCall(double S1,double S2,double sig1,double sig2,double rho,double k,double t,int n);

double SpreadDigitalOption(double S1,double S2,double sig1,double sig2,double rho,double k,double t,int callput,int n);

double Vega1SpreadOption(double S1,double S2,double sig1,double sig2,double rho,double k,double t,int callput,int n);

double Vega2SpreadOption(double S1,double S2,double sig1,double sig2,double rho,double k,double t,int callput,int n);

double  LogNormal_SpreadOption_Calibrate_Correlation(double S1,double S2,double sig1,double sig2,double optionprice,
											 double k,double t,int callput,int optiontype,int n);
double  Smiled_LogNormal_SpreadOption_Calibrate_Correlation(double S1,double S2,double sig1,double sig2,double optionprice,
											 double k,double t,double slope1,double slope2,int callput,int optiontype,int n);


double CorrelationRobust_SpreadDigitalCall(double S10,double S20,double sig10,double sig20,double rho0,double k0,double t0,int n);



CC_END_NAMESPACE()



#endif

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/

