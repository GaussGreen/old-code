
/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 * Version initiale 01/15/2004
 *
 *  Functions for   BlackSholes 
 *
 *	\file vanille_bs_interface.cpp
 *
 *  \brief
 *
 *	\author  O. Croissant
 *	\version 1.0
 *	\date January 2004
 */

#include <glob/firsttoinc.h>
#include "gpbase/port.h"
#include <math.h>

#include "gpclosedforms/vanille_bs_formula.h"
#include "gpclosedforms/vanilla_bs.h"

CC_BEGIN_NAMESPACE( ARM )



///////////////////////////////////////////////////////////////////////
///  
///			Export Pricing Functions 
///
///////////////////////////////////////////////////////////////////////
double Export_BlackSholes(double forward,
							double totalvolatility,
							double bondprice,
							double strike,
							double CallPut)
{
	ArgumentList a(forward,totalvolatility,bondprice,strike,CallPut);
	
	Power_Expression<ARM_CF_BS_Formula> y;
	return y(a);
}

double Export_BlackSholes(int i, double forward,
							double totalvolatility,
							double bondprice,
							double strike,
							double CallPut)
{
	ArgumentList a(forward,totalvolatility,bondprice,strike,CallPut);
	
	Power_Expression<ARM_CF_BS_Formula> y;
	return y(i,a);
}


double Export_BlackSholes(int i,int j, double forward,
							double totalvolatility,
							double bondprice,
							double strike,
							double CallPut)
{
	ArgumentList a(forward,totalvolatility,bondprice,strike,CallPut);
	
	Power_Expression<ARM_CF_BS_Formula> y;
	return y(i,j,a);
}


double Export_BlackSholes_ImplicitVol(
							double forward,
							double discount,
							double strike,
							double CallPut,
							double optprice,
							int algorithm = 1)
{
	double intrinsecvalue;
	if (CallPut==K_CALL) 
	{
		intrinsecvalue=(forward-strike)>0.0?forward-strike:0.0;
	}
	else
	{
		intrinsecvalue=(strike-forward)>0.0?strike-forward:0.0;
	}
	if((optprice-intrinsecvalue)<forward/1000000.)
	{
		return Asymptotic_BlackSholesTimeValue_ImplicitVol(forward,strike,(optprice-intrinsecvalue)/discount);
	}
	
	else
	{
		
		if(algorithm==1)
		{
			return ARM_CF_BS_Formula::callimplicit_totalvolatility(forward, discount, strike, CallPut,optprice,1e-10);
		}
		else if(algorithm==2)
		{
			return ARM_CF_BS_Formula::callimplicit_totalvolatility1(forward, discount, strike, CallPut,optprice,1e-10);
		}
		else if(algorithm==3)
		{
			return ARM_CF_BS_Formula::callimplicit_totalvolatility2(forward, discount, strike, CallPut,optprice,1e-10);
		}
		else 
		{
			if(fabs(strike/forward-1.0)<0.1)
			{
				return ARM_CF_BS_Formula::callimplicit_totalvolatility2(forward, discount, strike, CallPut,optprice,1e-10);
			}
			else
			{
				return ARM_CF_BS_Formula::callimplicit_totalvolatility(forward, discount, strike, CallPut,optprice,1e-10);
				
			}
			
		}
	}
}



double Export_BlackSholesTimeValue(double forward,
							double totalvolatility,
							double bondprice,
							double strike,
							double CallPut)
{
	if(fabs(log(forward/strike)/totalvolatility)<5.)
	{
		double intrinsec_value=forward>strike ?(forward-strike)*bondprice : 0;
		return Export_BlackSholes(forward,totalvolatility,bondprice,strike,K_CALL)-intrinsec_value;
	}
	else
	{
		return bondprice*Asymptotic_BlackSholesTimeValue(forward,strike,totalvolatility);
	}

}


double Export_BlackSholesTimeValue_ImplicitVol(
							double forward,
							double discount,
							double strike,
							double CallPut,
							double optprice)
{
	if((optprice/forward)>0.0001)
	{
		double intrinsec_value=forward>strike ?(forward-strike)*discount : 0;
		return Export_BlackSholes_ImplicitVol(forward,discount,strike,K_CALL,optprice+intrinsec_value);
	}
		else
	{
		return Asymptotic_BlackSholesTimeValue_ImplicitVol(forward,strike,optprice/discount);
	}
	
}


double Export_BlackScholesDigitalOption(double forward, double strike, double maturity, int callput,double volatility)
{
	return 	DigitalBlackSholes_Formula(	 forward,
							 volatility*sqrt(maturity),
							 1.0,
							 strike,
							 callput);
}

double Export_LN_RatioOption(double S1, double Mu1, double Sigma1,double S2, double Mu2, double Sigma2, double Rho, double K,double T, int CallPut)
{
	return 	BS_RatioOption(	 S1,Mu1,Sigma1,S2,Mu2,Sigma2,Rho,K,T,CallPut
							 );
}

double Export_LN_ProductOption(double S1, double Mu1, double Sigma1,double S2, double Mu2, double Sigma2, double Rho, double K,double T, int CallPut)
{
	return 	BS_ProductOption(	 S1,Mu1,Sigma1,S2,Mu2,Sigma2,Rho,K,T,CallPut
							 );
}
///////////////////////////////////////////////////////////////////////
///  
///			End of Export Pricing Functions 
///
///////////////////////////////////////////////////////////////////////



CC_END_NAMESPACE()


/*---------------------------------------------------------------------------*/
/*---- End of file ----*/