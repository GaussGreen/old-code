/*!
 *
 * Copyright (c) CDC IXIS CM November 2003 Paris
 *
 *	\file normalinvcum.h
 *  \brief General file for the inverse of the cumulative normal distribution
 *
 *	\author  E. Benhamou
 *	\version 1.0
 *	\date March 2004
 */

#ifndef _INGPNUMLIB_NORMALINVCUM_H
#define _INGPNUMLIB_NORMALINVCUM_H

#include "gpbase/port.h"
#include "gpbase/env.h"

CC_BEGIN_NAMESPACE( ARM )

//////////////////////////////
/// \class ARM_
//////////////////////////////
struct ARM_NormalInvCum
{
	static double HUGE_NB;
	static double Inverse_erf_Moro( double x );
	static double Inverse_erf( double x );
	static double Inverse_erf_conditional( double x, double a );
};

CC_END_NAMESPACE()

#endif
//-----------------------------------------------------------------------------
/*---- End of file ----*/
