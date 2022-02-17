/*!
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 * Version initiale 01/15/2004
 *
 *  basic functions for the closed form framework 
 *
 *	\Barriere_bs.cpp
 *
 *  \brief
 *
 *	\author  O. Croissant
 *	\version 1.0
 *	\date January 2004
 */

#include <glob/firsttoinc.h>
#include "gpbase/port.h"


#include <cmath>

#include"gpclosedforms/vanilla_bs.h"


CC_BEGIN_NAMESPACE(ARM)

double shifted_lognormal_vanilla_call(
									  const double index,
									  const double strike,
									  const double timetomaturity,
									  const double volatility,
									  const double alpha)
{
	return BS(index+alpha,strike+alpha,timetomaturity,volatility);
}

double RelativeShiftBSVanilla(double fwd, double strike, double ttm, int callPut, double shift, double vol)
{
	double F = fwd/fabs(shift);
	double K = shift < 0. ? F * (1. + shift) - strike : strike - (shift - 1.) * F;
	
	return BS(F, K, ttm, vol, callPut);
}



CC_END_NAMESPACE()
 

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/
