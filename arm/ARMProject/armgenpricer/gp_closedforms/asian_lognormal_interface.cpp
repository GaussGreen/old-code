/*!
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 * Version initiale 01/15/2004
 *
 *  basic functions for the closed form framework 
 *
 *	\file heston.cpp
 *
 *  \brief
 *
 *	\author  O. Croissant
 *	\version 1.0
 *	\date January 2004
 */
#include "firsttoinc.h"
#include "gpbase/port.h"

#include "gpclosedforms/asian_lognormal_formula.h"



///////////////////////////////////////////////////////////////////////////////
///
///					Process :
///			dS= S(rdt+V^(1/2) dW1+ Jdq
///			dV=(omega-theta*V)dt +ksi*V^(1/2) dW2 
///			dW1.dW2=rho*dt
///
///			Jump J : probability lambda, volatility sigmaJ, log-size muJ
///
/////////////////////////////////////////////////////////////////////////////:


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

double Export_Lognormal_Asian_VanillaOption(
											double S,
											double k,
											double T,
											double r,
											double v,
											double alpha,
											int callput,
											int n)
{
	ArgumentList a(S,k,T,r,v,callput,alpha,n);

	Power_Expression<ARM_CF_Asian_Vanilla_Lognormal_Formula> y;
	return y(a);
}

///////////////////////////////////////////////////////////////////////
///  
///			1st Derivatives  
///
///////////////////////////////////////////////////////////////////////
double Export_Lognormal_Asian_VanillaOption(
											int i,
											double S,
											double k,
											double T,
											double r,
											double v,
											double alpha,
											int callput,
											int n)
{
	ArgumentList a(S,k,T,r,v,callput,alpha,n);
	Power_Expression<ARM_CF_Asian_Vanilla_Lognormal_Formula> y;
	return y(i,a);
}

///////////////////////////////////////////////////////////////////////
///  
///			2nd Derivatives  
///
///////////////////////////////////////////////////////////////////////
double Export_Lognormal_Asian_VanillaOption(int i,
											int j,
											double S,
											double k,
											double T,
											double r,
											double v,
											double alpha,
											int callput,
											int n)
{
	ArgumentList a(S,k,T,r,v,callput,alpha,n);
	Power_Expression<ARM_CF_Asian_Vanilla_Lognormal_Formula> y;
	return y(i,j,a);
}


CC_END_NAMESPACE()
 


/*---------------------------------------------------------------------------*/
/*---- End of file ----*/
