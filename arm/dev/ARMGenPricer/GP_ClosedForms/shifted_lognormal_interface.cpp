/*!
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 * Version initiale 01/15/2004
 *
 *  basic functions for the closed form framework 
 *
 *	\file shifted_lognormal_interface.cpp
 *
 *  \brief
 *
 *	\author  O. Croissant
 *	\version 1.0
 *	\date January 2004
 */
#include <glob/firsttoinc.h>
#include "gpbase/port.h"


#include "gpclosedforms/shifted_lognormal_formula.h"


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

double Export_shifted_lognormal_VanillaOption(
									double F,
									double K,
									double V0,
									double t,
									double alpha,
									int callorput
									)
{
	ArgumentList a(F,K,t,V0,alpha,callorput);
	
	Power_Expression<ARM_CF_Shifted_Lognormal_Formula> y;
	return y(a);
}

///////////////////////////////////////////////////////////////////////
///  
///			1st Derivatives  
///
///////////////////////////////////////////////////////////////////////
double Export_shifted_lognormal_VanillaOption(int i,
									double F,
									double K,
									double V0,
									double t,
									double alpha,
									int callorput
									)
{
ArgumentList a(F,K,t,V0,alpha,callorput);
	
	Power_Expression<ARM_CF_Shifted_Lognormal_Formula> y;
	return y(i,a);
}

///////////////////////////////////////////////////////////////////////
///  
///			2nd Derivatives  
///
///////////////////////////////////////////////////////////////////////
double Export_shifted_lognormal_VanillaOption(int i,int j,
									double F,
									double K,
									double V0,
									double t,
									double alpha,
									int callorput
									)
{
	ArgumentList a(F,K,t,V0,alpha,callorput);
	
	Power_Expression<ARM_CF_Shifted_Lognormal_Formula> y;
	return y(i,j,a);
}




CC_END_NAMESPACE()
 


/*---------------------------------------------------------------------------*/
/*---- End of file ----*/
