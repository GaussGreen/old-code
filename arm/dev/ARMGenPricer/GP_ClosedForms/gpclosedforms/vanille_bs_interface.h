/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 * Version initiale 01/15/2004
 *
 *  Header for the   BlackSholes templates and related
 *
 *	\file vanille_bs_interface.h
 *
 *  \brief
 *
 *	\author  O. Croissant
 *	\version 1.0
 *	\date January 2004
 */

#ifndef _GP_CF_VANILLE_BS_INTERFACE_H
#define _GP_CF_VANILLE_BS_INTERFACE_H

#include <glob/firsttoinc.h>
#include "gpbase/port.h"

CC_BEGIN_NAMESPACE( ARM )

///////////////////////////////////////////////////////////////////////
///  
///			Utility Functions 
///
///////////////////////////////////////////////////////////////////////

//double BS(double f, double k, double t, double v);

double VanillaOptionFromCall(double forward,double k,double call,int CallOrPutFlag);

///////////////////////////////////////////////////////////////////////
///  
///			Export Pricing Functions 
///
///////////////////////////////////////////////////////////////////////
double Export_BlackSholes(double forward,
							double totalvolatility,
							double bondprice,
							double strike,
							double CallPut);


double Export_BlackSholes(int i, double forward,
							double totalvolatility,
							double bondprice,
							double strike,
							double CallPut);



double Export_BlackSholes(int i,int j, double forward,
							double totalvolatility,
							double bondprice,
							double strike,
							double CallPut);



double Export_BlackSholes_ImplicitVol(double forward,
							double discount,
							double strike,
							double CallPut,
							double optprice,
							int algorithm=1);

double Export_BlackSholesTimeValue(double forward,
							double totalvolatility,
							double bondprice,
							double strike,
							double CallPut);


double Export_BlackSholesTimeValue_ImplicitVol(
							double forward,
							double discount,
							double strike,
							double CallPut,
							double optprice);

double Export_BlackScholesDigitalOption(
							double forward,
							double strike,
							double maturity,
							int callput,
							double volatility);

double Export_LN_RatioOption(double S1, double Mu1, double Sigma1,double S2, double Mu2, double Sigma2, double Rho, double K,double T, int CallPut);


double Export_LN_ProductOption(double S1, double Mu1, double Sigma1,double S2, double Mu2, double Sigma2, double Rho, double K,double T, int CallPut);




CC_END_NAMESPACE()

#endif

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/
