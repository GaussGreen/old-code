/*!
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 * Version initiale 01/15/2004
 *
 *  basic functions for the closed form framework 
 *
 *	\file vanille_normal_interface.cpp
 *
 *  \brief
 *
 *	\author  O. Croissant
 *	\version 1.0
 *	\date January 2004
 */
#include "firsttoinc.h"
#include "gpbase/port.h"


#include "gpclosedforms/vanille_normal_formula.h"
#include "gpclosedforms/vanilla_normal.h"

CC_BEGIN_NAMESPACE(ARM)



 ///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
///  
///			Begining Exportable  Pricing Functions 
///
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////



/// callput =  1 (K_CALL) for call
/// callput = -1 (K_PUT) for put

double Export_normal_VanillaOption(
									double F,
									double K,
									double V0,
									double t,
									int callorput
									)
{
	ArgumentList a(F,V0,K,t,callorput);
	
	Power_Expression<ARM_CF_Vanille_Normal_Formula> y;
	return y(a);
}

///////////////////////////////////////////////////////////////////////
///  
///			1st Derivatives  
///
///////////////////////////////////////////////////////////////////////
double Export_normal_VanillaOption(int i,
									double F,
									double K,
									double V0,
									double t,
									int callorput
									)
{
	ArgumentList a(F,V0,K,t,callorput);
	
	Power_Expression<ARM_CF_Vanille_Normal_Formula> y;
	return y(i,a);
}

///////////////////////////////////////////////////////////////////////
///  
///			2nd Derivatives  
///
///////////////////////////////////////////////////////////////////////
double Export_normal_VanillaOption(int i,int j,
									double F,
									double K,
									double V0,
									double t,
									int callorput
									)
{
	ArgumentList a(F,V0,K,t,callorput);
	
	Power_Expression<ARM_CF_Vanille_Normal_Formula> y;
	return y(i,j,a);
}

///////////////////////////////////////////////////////////////////////
///  
///			implicit volatility 
///
///////////////////////////////////////////////////////////////////////

double Export_Normal_ImpliedVol(double F,double K,double opt, int CallPut)
{
	double val_intrinseque=(F-K)*CallPut;
	val_intrinseque=(val_intrinseque>0) ? val_intrinseque : 0;
	if(opt<val_intrinseque)
	{
			throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"Export_Normal_ImpliedVol: Val intrinseque> option !");
	
	}
	if(CallPut==K_CALL)
	{
		return VanillaCall_N_BSphi_Inverse(F-K,opt,1e-12);
	}
	else
	{
		double call=opt+(F-K);
		if(call<=0)
		{
			throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"Export_Normal_ImpliedVol : call put parity problem! ");
		}
		return VanillaCall_N_BSphi_Inverse(F-K,call,1e-12);
	}	
}

double Export_Normal_Digital_ImpliedVol(double F,double K,double opt, int CallPut)
{
	return DigitalCall_N_ImpliedVol( F, K,  opt, CallPut);
}


double Export_Normal_DoubleDigital(double fwd1, 
								   double fwd2,
								   double maturity,
								   double K1, double spread1,
								   double K2, double spread2,
								   double vol1plus, double vol1minus,
								   double vol2plus, double vol2minus,
								   double correl,
								   int callorput1,
								   int callorput2)
{
	return DoubleDigital_N( maturity, K1, spread1, K2, spread2, fwd1, vol1plus, vol1minus, fwd2, vol2plus, vol2minus, correl, callorput1, callorput2);
}

CC_END_NAMESPACE()
 


/*---------------------------------------------------------------------------*/
/*---- End of file ----*/
