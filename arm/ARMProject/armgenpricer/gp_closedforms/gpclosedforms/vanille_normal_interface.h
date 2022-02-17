/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 * Version initiale 01/15/2004
 *
 *  basic functions for the closed form framework 
 *
 *	\file vanille_normal_interface.h
 *
 *  \brief
 *
 *	\author  O. Croissant
 *	\version 1.0
 *	\date January 2004
 */
///////////////////////////////////////////////////////////////////////////////
///
///					Process :
///			dS= rdt+ sigma*dW
///			
///
/////////////////////////////////////////////////////////////////////////////:

 
#ifndef _GP_CF_VANILLE_NORMAL_INTERFACE_H
#define _GP_CF_VANILLE_NORMAL_INTERFACE_H


#include "firsttoinc.h"
#include "gpbase/port.h"

CC_BEGIN_NAMESPACE(ARM)

/////////////////////////////////////////////////////////////////////////////////////////////////////
///  
///			Begining Exportable  Pricing Functions 
///
/////////////////////////////////////////////////////////////////////////////////////////////////////
double Export_normal_VanillaOption(
									double F,
									double K,
									double V0,
									double t,
									int callorput
									);

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
									);

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
									);


///////////////////////////////////////////////////////////////////////
///  
///			implicit volatility 
///
///////////////////////////////////////////////////////////////////////
double Export_Normal_ImpliedVol(double F,double K,double opt, int CallPut);
double Export_Normal_Digital_ImpliedVol(double F,double K,double opt, int CallPut);


///////////////////////////////////////////////////////////////////////
///  
///			Double digitale
///
///////////////////////////////////////////////////////////////////////
double Export_Normal_DoubleDigital(double fwd1, 
								   double fwd2,
								   double maturity,
								   double K1, double spread1,
								   double K2, double spread2,
								   double vol1plus, double vol1minus,
								   double vol2plus, double vol2minus,
								   double correl,
								   int callorput1,
								   int callorput2);


CC_END_NAMESPACE()



#endif

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/
