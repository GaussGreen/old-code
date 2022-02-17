/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 * Version initiale 01/15/2004
 *
 *  Header for the   BlackSholes templates and related
 *
 *	\file stochasticvol_ln_interface.h
 *
 *  \brief
 *
 *	\author  O. Croissant
 *	\version 1.0
 *	\date January 2004
 */

#ifndef _GP_CF_STOCHASTICVOL_LN_INTERFACE_H
#define _GP_CF_STOCHASTICVOL_LN_INTERFACE_H

#include "gpbase/port.h"
#include <gpbase/argconv.h>
#include "basic_distributions.h"

CC_BEGIN_NAMESPACE( ARM )

double Export_StochasticVol_LN_Geometric_VanillaOption_with_Reset(	double forward,
							double strike,
							double maturity,
							double drift,
							double volatility,
							double voldrift,
							double volvol,
							double averaging,
							double reset,
							double CallPut,
							int nbhermite
							);
double Export_StochasticVol_LN_Arithmetic_VanillaOption_with_Reset(	double forward,
							double strike,
							double maturity,
							double drift,
							double volatility,
							double voldrift,
							double volvol,
							double averaging,
							double reset,
							double CallPut,
							int nbhermite
							);

double Export_StochasticVol_LN_Geometric_VanillaOption_with_Reset(int i,double forward,
							double strike,
							double maturity,
							double drift,
							double volatility,
							double voldrift,
							double volvol,
							double averaging,
							double reset,
							double CallPut,
							int nbhermite
							);

double Export_StochasticVol_LN_Arithmetic_VanillaOption_with_Reset(int i,double forward,
							double strike,
							double maturity,
							double drift,
							double volatility,
							double voldrift,
							double volvol,
							double averaging,
							double reset,
							double CallPut,
							int nbhermite
							);

double Export_StochasticVol_LN_Arithmetic_VanillaOption_with_Reset(int i,int j,double forward,
							double strike,
							double maturity,
							double drift,
							double volatility,
							double voldrift,
							double volvol,
							double averaging,
							double reset,
							double CallPut,
							int nbhermite
							);

double Export_StochasticVol_LN_Geometric_VanillaOption_with_Reset(int i,int j,double forward,
							double strike,
							double maturity,
							double drift,
							double volatility,
							double voldrift,
							double volvol,
							double averaging,
							double reset,
							double CallPut,
							int nbhermite
							);


CC_END_NAMESPACE()

#endif

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/
