/*!
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 * Version initiale 01/15/2004
 *
 *  basic functions for the closed form framework 
 *
 *	\file change_numeraire.cpp
 *
 *  \brief
 *
 *	\author  O. Croissant
 *	\version 1.0
 *	\date January 2004
 */

#include "firsttoinc.h"
#include "gpbase/port.h"

#include "gpclosedforms/spreadoption_normal_formula.h"


CC_BEGIN_NAMESPACE(ARM)


///////////////////////////////////////////////////////////////////////
///  
///			Export Pricing Functions 

/// callput =  1  for call
/// callput = -1  for put

///  optiontype=0	Digital			(condition S1-S2-k>=0)
///  optiontype=1	pays S1			(condition S1-S2-k>=0)
///  optiontype=2	pays S2			(condition S1-S2-k>=0)
///  optiontype=3	pays S1-S2-k	(condition S1-S2-k>=0)


///
///////////////////////////////////////////////////////////////////////
double Export_Normal_SpreadOption(double S1,double S2,double sig1,double sig2,double rho,double k,double t,int callput,int optiontype)
{
	ArgumentList a(S1,S2,sig1,sig2,rho,k,t,callput,optiontype);
	Power_Expression<ARM_CF_SpreadDigitalOption_N_Formula> y;
	return y(a);
}

double Export_Smiled_Normal_SpreadOption(double S1,double S2,double sig1,double sig2,double rho,double k,double t,
										 double slope1, double slope2,int callput,int optiontype)
{
	ArgumentList a(S1,S2,sig1,sig2,rho,k,t,callput,optiontype,slope1,slope2);
			Power_Expression<ARM_CF_Smiled_SpreadDigitalOption_N_Formula> y;
			return y(a);
}

///////////////////////////////////////////////////////////////////////
///  
///			1st Derivatives  
///
///////////////////////////////////////////////////////////////////////
double Export_Normal_SpreadOption(int i, double S1,double S2,double sig1,double sig2,double rho,double k,double t,int callput,int optiontype)
{
	ArgumentList a(S1,S2,sig1,sig2,rho,k,t,callput,optiontype);
	Power_Expression<ARM_CF_SpreadDigitalOption_N_Formula> y;
	return y(i,a);
}


double Export_Smiled_Normal_SpreadOption(int i, double S1,double S2,double sig1,double sig2,double rho,double k,double t,
										 double slope1, double slope2,int callput,int optiontype)
{
	ArgumentList a(S1,S2,sig1,sig2,rho,k,t,callput,optiontype,slope1,slope2);
			Power_Expression<ARM_CF_Smiled_SpreadDigitalOption_N_Formula> y;
			return y(i,a);
}



///////////////////////////////////////////////////////////////////////
///  
///			2nd Derivatives  
///
///////////////////////////////////////////////////////////////////////
double Export_Normal_SpreadOption(int i, int j, double S1,double S2,double sig1,double sig2,double rho,double k,double t,int callput,int optiontype)
{
	ArgumentList a(S1,S2,sig1,sig2,rho,k,t,callput,optiontype);
	Power_Expression<ARM_CF_SpreadDigitalOption_N_Formula> y;
	return y(i,j,a);
}


double Export_Smiled_Normal_SpreadOption(int i, int j,double S1,double S2,double sig1,double sig2,double rho,double k,double t,
										 double slope1, double slope2,int callput,int optiontype)
{
	ArgumentList a(S1,S2,sig1,sig2,rho,k,t,callput,optiontype,slope1,slope2);
	Power_Expression<ARM_CF_Smiled_SpreadDigitalOption_N_Formula> y;
	return y(i,j,a);
}



CC_END_NAMESPACE()
 

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/

