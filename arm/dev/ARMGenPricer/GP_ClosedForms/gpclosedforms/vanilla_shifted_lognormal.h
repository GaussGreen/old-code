/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 * Version initiale 01/15/2004
 *
 *  basic functions for the closed form framework 
 *
 *	\file vanilla_shiftedlognormal.h
 *
 *  \brief
 *
 *	\author  O. Croissant
 *	\version 1.0
 *	\date January 2004
 */
 
#ifndef _GP_CF_VANILLA_SHIFTEDLOGNORMAL_H
#define _GP_CF_VANILLA_SHIFTEDLOGNORMAL_H

#include <glob/firsttoinc.h>
#include "gpbase/port.h"


CC_BEGIN_NAMESPACE(ARM)


double shifted_lognormal_vanilla_call(
									  const double index,
									  const double strike,
									  const double timetomaturity,
									  const double volatility,
									  const double alpha);



double RelativeShiftBSVanilla(double fwd, double strike, double ttm, int callPut, double shift, double vol);

CC_END_NAMESPACE()



#endif

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/

