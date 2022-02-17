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
 
#ifndef _GP_CF_SPREADOPTION_LOGNORMAL_INTERFACE_H
#define _GP_CF_SPREADOPTION_LOGNORMAL_INTERFACE_H

#include "firsttoinc.h"
#include "gpbase/port.h"



CC_BEGIN_NAMESPACE(ARM)


///////////////////////////////////////////////////////////////////////
///  
///			Begining Exportable  Pricing Functions 
///
///////////////////////////////////////////////////////////////////////


/// callput =  1  (K_CALL) for call
/// callput = -1  (K_PUT) for put

///  optiontype=0	Digital			(condition S1-S2-k>=0)
///  optiontype=1	pays S1-S2-k	(condition S1-S2-k>=0)
///  optiontype=2	pays S1			(condition S1-S2-k>=0)
///  optiontype=3	pays S2	(condition S1-S2-k>=0)

/////////////////////////////////////////////////////////////////////////////////////////////////////
///
///   Forward value
///
///
/////////////////////////////////////////////////////////////////////////////////////////////////////


double Export_LogNormal_SpreadOption(double S1,double S2,double sig1,double sig2,double rho,double k,double t,int callput,int optiontype,int n=64);


/////////////////////////////////////////////////////////////////////////////////////////////////////
///
///   Forward value : the smile is brought by the its slope (2 additional parameters slope1 and slope2)
///
///
/////////////////////////////////////////////////////////////////////////////////////////////////////
double Export_Smiled_LogNormal_SpreadOption(double S1,double S2,double sig1,double sig2,double rho,double k,
										 double t,double slope1, double slope2,int callput,int optiontype,int n=64);

///////////////////////////////////////////////////////////////////////
///  
///			1st Derivatives  
///
///////////////////////////////////////////////////////////////////////
double Export_LogNormal_SpreadOption(int i,double S1,double S2,double sig1,double sig2,double rho,double k,double t,int callput,int optiontype,int n=64);

/////////////////////////////////////////////////////////////////////////////////////////////////////
///
///   1st Derivatives : the smile is brought by the its slope (2 additional parameters slope1 and slope2)
///
///
/////////////////////////////////////////////////////////////////////////////////////////////////////
double Export_Smiled_LogNormal_SpreadOption(int i,double S1,double S2,double sig1,double sig2,double rho,double k,
										 double t,double slope1, double slope2 ,int callput,int optiontype,int n=64);

///////////////////////////////////////////////////////////////////////
///  
///			2nd Derivatives  
///
///////////////////////////////////////////////////////////////////////
double Export_LogNormal_SpreadOption(int i,int j,double S1,double S2,double sig1,double sig2,double rho,double k,double t,int callput,int optiontype,int n=64);


/////////////////////////////////////////////////////////////////////////////////////////////////////
///
///   2nd Derivatives : the smile is brought by the its slope (2 additional parameters slope1 and slope2)
///
///
/////////////////////////////////////////////////////////////////////////////////////////////////////

double Export_Smiled_LogNormal_SpreadOption(int i,int j,double S1,double S2,double sig1,double sig2,double rho,double k,
								  double t,double slope1, double slope2 ,int callput,int optiontype,int n=64);

////////////////////////////////////////////////////////////////////////////////////////////////////////
///  
///			Calibration of the implicit correlation for option 
///
////////////////////////////////////////////////////////////////////////////////////////////////////////



double Export_LogNormal_SpreadOption_Calibrate_Correlation(double S1,double S2,double sig1,double sig2,double optionprice,
											 double k,double t,int callput,int optiontype,int n=64);


double Export_Smiled_LogNormal_SpreadOption_Calibrate_Correlation(double S1,double S2,double sig1,double sig2,double optionprice,
											 double k,double t,double slope1, double slope2,int callput,int optiontype,int n);
///////////////////////////////////////////////////////////////////////
///  
///			End Exportable Pricing Functions
///
///////////////////////////////////////////////////////////////////////


CC_END_NAMESPACE()



#endif

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/

