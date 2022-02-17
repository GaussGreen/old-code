/*!
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 * Version initiale 01/15/2004
 *
 *  basic functions for the closed form framework 
 *
 *	\file spreadoption_lognormal.cpp
 *
 *  \brief
 *
 *	\author  O. Croissant
 *	\version 1.0
 *	\date January 2004
 */
#include "firsttoinc.h"
#include "gpbase/port.h"

#include "gpclosedforms/spreadoption_lognormal_formula.h"
#include "gpclosedforms/spreadoption_lognormal.h"
#include "gpclosedforms/vanille_bs_formula.h"
#include "gpclosedforms/vanilla_bs.h"


inline double cmin(double x, double y) {return ((x<=y) ? x : y);}
inline double cmax(double x, double y) {return ((x<=y) ? y : x);}



CC_BEGIN_NAMESPACE(ARM)

///////////////////////////////////////////////////////////////////////
///  
///			Export Pricing Functions 
///
///////////////////////////////////////////////////////////////////////
double Export_LogNormal_SpreadOption(double S1,double S2,double sig1,double sig2,double rho,double k,double t,int callput,int optiontype,int n)
{
	ArgumentList a(t,S1,S2,sig1,sig2,rho,k,callput,n);
	switch (optiontype) 
	{
	case ARM_CF_SpreadDigitalOption_Formula::DIGITALOPTION :
		{
			Power_Expression<ARM_CF_SpreadDigitalOption_Formula> y;
			return cmax(y(a),0);
			break;
		}
	case ARM_CF_SpreadDigitalOption_Formula::SPREADOPTION :
		{
			if(S2<1e-10*S1)
			{
				if(k<=0)
				{
					return cmax(S1-k,0);
				}
				return BlackSholes_Formula(	S1,sig1*sqrt(t),1.0,k,callput);
			}
			if(S1<1e-10*S2)
			{
				if(k>=0)
				{
					return 0;
				}
				return BlackSholes_Formula(	S2,sig2*sqrt(t),1.0,-k,-callput);
			}
			if(fabs(sig2)<0.00000001)
			{
				ArgumentList a(S1,sig1*sqrt(t),1.0,k+S2,callput);
				Power_Expression<ARM_CF_BS_Formula> y;
				return cmax(y(a),0);
			}
			else
			if(fabs(sig1)<0.00000001)
			{
				ArgumentList a(S2,sig2*sqrt(t),1.0,S1-k,-callput);
				Power_Expression<ARM_CF_BS_Formula> y;
				return cmax(y(a),0);
			}
			else
			{
				Power_Expression<ARM_CF_SpreadOption_Formula> y;
				return cmax(y(a),0);
			}
			break;
		}
	case ARM_CF_SpreadDigitalOption_Formula::PAYFIRST :
		{
			Power_Expression<ARM_CF_Index1Paying_SpreadDigitalOption_Formula> y;
			return cmax(y(a),0);
			break;
		}
	case ARM_CF_SpreadDigitalOption_Formula::PAYSECOND :
		{
			Power_Expression<ARM_CF_Index2Paying_SpreadDigitalOption_Formula> y;
			return cmax(y(a),0);
			break;
		}
	default :
		{
			throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"Export_LogNormal_SpreadOption : incorrect optiontype toggle ");
		}
	}
}

double Export_Smiled_LogNormal_SpreadOption(double S1,double S2,double sig1,double sig2,double rho,double k,double t,
										 double slope1, double slope2,int callput,int optiontype,int n)
{
	ArgumentList a(t,S1,S2,sig1,sig2,rho,k,callput,n,slope1,slope2);
	ArgumentList a0(t,S1,S2,sig1,sig2,rho,k,callput,n);
	switch (optiontype) 
	{
	case ARM_CF_SpreadDigitalOption_Formula::DIGITALOPTION :
		{
			Power_Expression<ARM_CF_Smiled_SpreadDigitalOption_Formula> y;
			return y(a);
			break;
		}
	case ARM_CF_SpreadDigitalOption_Formula::SPREADOPTION :
		{
			Power_Expression<ARM_CF_SpreadOption_Formula> y;
			return y(a0);
			break;
		}
	case ARM_CF_SpreadDigitalOption_Formula::PAYFIRST :
		{
			Power_Expression<ARM_CF_Index1Paying_Smiled_SpreadDigitalOption_Formula> y;
			return y(a);
			break;
		}
	case ARM_CF_SpreadDigitalOption_Formula::PAYSECOND :
		{
			Power_Expression<ARM_CF_Index2Paying_Smiled_SpreadDigitalOption_Formula> y;
			return y(a);
			break;
		}
	default :
		{
			throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"Export_Smiled_LogNormal_SpreadOption : incorrect optiontype toggle ");
		}
	}
}

///////////////////////////////////////////////////////////////////////
///  
///			1st Derivatives  
///
///////////////////////////////////////////////////////////////////////
double Export_LogNormal_SpreadOption(int i,double S1,double S2,double sig1,double sig2,double rho,double k,double t,int callput,int optiontype,int n)
{
	ArgumentList a(t,S1,S2,sig1,sig2,rho,k,callput,n);
	switch (optiontype) 
	{
	case ARM_CF_SpreadDigitalOption_Formula::DIGITALOPTION :
		{
			Power_Expression<ARM_CF_SpreadDigitalOption_Formula> y;
			return y(i,a);
			break;
		}
	case ARM_CF_SpreadDigitalOption_Formula::SPREADOPTION :
		{
			Power_Expression<ARM_CF_SpreadOption_Formula> y;
			return y(i,a);
			break;
		}
	case ARM_CF_SpreadDigitalOption_Formula::PAYFIRST :
		{
			Power_Expression<ARM_CF_Index1Paying_SpreadDigitalOption_Formula> y;
			return y(i,a);
			break;
		}
	case ARM_CF_SpreadDigitalOption_Formula::PAYSECOND :
		{
			Power_Expression<ARM_CF_Index2Paying_SpreadDigitalOption_Formula> y;
			return y(i,a);
			break;
		}
	default :
		{
			throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"Export_LogNormal_SpreadOption : incorrect optiontype toggle ");
		}
	}
}

double Export_Smiled_LogNormal_SpreadOption(int i, double S1,double S2,double sig1,double sig2,double rho,double k,double t,
										 double slope1, double slope2,int callput,int optiontype,int n)
{
	ArgumentList a(t,S1,S2,sig1,sig2,rho,k,callput,n,slope1,slope2);
	ArgumentList a0(t,S1,S2,sig1,sig2,rho,k,callput,n);
	switch (optiontype) 
	{
	case ARM_CF_SpreadDigitalOption_Formula::DIGITALOPTION :
		{
			Power_Expression<ARM_CF_Smiled_SpreadDigitalOption_Formula> y;
			return y(i,a);
			break;
		}
	case ARM_CF_SpreadDigitalOption_Formula::SPREADOPTION :
		{
			Power_Expression<ARM_CF_SpreadOption_Formula> y;
			return y(i,a0);
			break;
		}
	case ARM_CF_SpreadDigitalOption_Formula::PAYFIRST :
		{
			Power_Expression<ARM_CF_Index1Paying_Smiled_SpreadDigitalOption_Formula> y;
			return y(i,a);
			break;
		}
	case ARM_CF_SpreadDigitalOption_Formula::PAYSECOND :
		{
			Power_Expression<ARM_CF_Index2Paying_Smiled_SpreadDigitalOption_Formula> y;
			return y(i,a);
			break;
		}
	default :
		{
			throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"Export_Smiled_LogNormal_SpreadOption : incorrect optiontype toggle ");
		}
	}
}

///////////////////////////////////////////////////////////////////////
///  
///			2nd Derivatives  
///
///////////////////////////////////////////////////////////////////////
double Export_LogNormal_SpreadOption(int i,int j,double S1,double S2,double sig1,double sig2,double rho,double k,double t,int callput,int optiontype,int n)
{
	ArgumentList a(t,S1,S2,sig1,sig2,rho,k,callput,n);
	switch (optiontype) 
	{
	case ARM_CF_SpreadDigitalOption_Formula::DIGITALOPTION :
		{
			Power_Expression<ARM_CF_SpreadDigitalOption_Formula> y;
			return y(i,j,a);
			break;
		}
	case ARM_CF_SpreadDigitalOption_Formula::SPREADOPTION :
		{
			Power_Expression<ARM_CF_SpreadOption_Formula> y;
			return y(i,j,a);
			break;
		}
	case ARM_CF_SpreadDigitalOption_Formula::PAYFIRST :
		{
			Power_Expression<ARM_CF_Index1Paying_SpreadDigitalOption_Formula> y;
			return y(i,j,a);
			break;
		}
	case ARM_CF_SpreadDigitalOption_Formula::PAYSECOND :
		{
			Power_Expression<ARM_CF_Index2Paying_SpreadDigitalOption_Formula> y;
			return y(i,j,a);
			break;
		}
	default :
		{
			throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"Export_LogNormal_SpreadOption : incorrect optiontype toggle ");
		}
	}
}

double Export_Smiled_LogNormal_SpreadOption(int i, int j,double S1,double S2,double sig1,double sig2,double rho,double k,double t,
										 double slope1, double slope2,int callput,int optiontype,int n)
{
	ArgumentList a(t,S1,S2,sig1,sig2,rho,k,callput,n,slope1,slope2);
	ArgumentList a0(t,S1,S2,sig1,sig2,rho,k,callput,n);
	switch (optiontype) 
	{
	case ARM_CF_SpreadDigitalOption_Formula::DIGITALOPTION :
		{
			Power_Expression<ARM_CF_Smiled_SpreadDigitalOption_Formula> y;
			return y(i,j,a);
			break;
		}
	case ARM_CF_SpreadDigitalOption_Formula::SPREADOPTION :
		{
			Power_Expression<ARM_CF_SpreadOption_Formula> y;
			return y(i,j,a0);
			break;
		}
	case ARM_CF_SpreadDigitalOption_Formula::PAYFIRST :
		{
			Power_Expression<ARM_CF_Index1Paying_Smiled_SpreadDigitalOption_Formula> y;
			return y(i,j,a);
			break;
		}
	case ARM_CF_SpreadDigitalOption_Formula::PAYSECOND :
		{
			Power_Expression<ARM_CF_Index2Paying_Smiled_SpreadDigitalOption_Formula> y;
			return y(i,j,a);
			break;
		}
	default :
		{
			throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"Export_Smiled_LogNormal_SpreadOption : incorrect optiontype toggle ");
		}
	}
}

///////////////////////////////////////////////////////////////////////
///  
///			Calibration of the implicit correlation for option 
///
///////////////////////////////////////////////////////////////////////



double Export_LogNormal_SpreadOption_Calibrate_Correlation(double S1,double S2,double sig1,double sig2,double optionprice,
											 double k,double t,int callput,int optiontype,int n)
{

	double borneinf=Export_LogNormal_SpreadOption( S1, S2, sig1, sig2, -0.995, k, t, callput, optiontype, n);

	double bornesup=Export_LogNormal_SpreadOption( S1, S2, sig1, sig2, 0.995, k, t, callput, optiontype, n);


	if (((optionprice<borneinf) && (optionprice<bornesup)) || ((optionprice>borneinf) && (optionprice>bornesup)))
	{
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
			"Export_LogNormal_SpreadOption_Calibrate_Correlation : optionprice  outside raisonable bound (corr=+ or -1)");
	}
	
	return LogNormal_SpreadOption_Calibrate_Correlation( S1, S2, sig1, sig2, optionprice, k, t, callput, optiontype, n);
}


double Export_Smiled_LogNormal_SpreadOption_Calibrate_Correlation(double S1,double S2,double sig1,double sig2,double optionprice,
											 double k,double t,double slope1, double slope2,int callput,int optiontype,int n)
{

	/// la valeur de l'option en fonction de rho n'est pas monotone donc 
	/// toute tentative de borner la valeur de l'option avec des correllation extermes
	/// n'est pas valable

	return Smiled_LogNormal_SpreadOption_Calibrate_Correlation( S1, S2, sig1, sig2, optionprice, k, t, slope1,  slope2, callput, optiontype, n);
}
CC_END_NAMESPACE()
 
#undef ARM_CF_SQRT_2_PI 
#undef ARM_CF_SQRT2 
#undef ARM_CF_INVSQRTPI

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/