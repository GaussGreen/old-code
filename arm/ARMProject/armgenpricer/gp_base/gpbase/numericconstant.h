/*!
 *
 * Copyright (c) CDC IXIS CM November 2004 Paris
 *
 *	\file numericconstant.h
 *
 *  \brief General file for the various numerical constants
 *	\author  E. Benhamou
 *	\version 1.0
 *	\date September 2004
 */

#ifndef _INGPNUMLIB_NUMERICCONSTANT_H
#define _INGPNUMLIB_NUMERICCONSTANT_H

#include "gpbase/port.h"

CC_BEGIN_NAMESPACE( ARM )

struct ARM_NumericConstants
{
	static const double ARM_ONE_THIRD;
	static const double ARM_PI;
	static const double ARM_2_PI;
	static const double ARM_PI_BY_2;
	static const double ARM_PI_BY_4;
	static const double ARM_SQRT_PI;
	static const double ARM_SQRT_2;
	static const double ARM_SQRT_3;
	static const double ARM_TWO_DIVIDED_BY_SQRT_3;
	static const double ARM_SQRT_2_DIVIDED_BY_SQRT_3;
	static const double ARM_SQRT_2_DIVIDED_BY_SQRT_PI;
	static const double ARM_SQRT_2_PI;
	static const double ARM_SQRT_8_PI;
	static const double ARM_INVSQRTPI;
	static const double ARM_INVSQRT2PI;
	static const double ARM_INVSQRTSQRTPI;

	static const double ARM_HALF_LOG_2PI;
	static const double ARM_LOG_2PI;
	static const double ARM_HALF_LOG_PI;

	static const double ARM_INVERSE_E;

	static const double ARM_BIGGEST_POSITIVE_NUMBER;
    static const double ARM_LOWEST_POSITIVE_NUMBER;
    static const double ARM_SQRT_LOWEST_POSITIVE_NUMBER;
    static const double ARM_INFINITY;
    static const double ARM_TOLERENCE;
    static const double ARM_STD_DOUBLE_TOLERANCE;
    static const double ARM_DOUBLE_TOLERENCE;
    static const size_t ARM_GP_MAX_ITER;
    static const double ARM_GP_GRAD_TOL;

};

CC_END_NAMESPACE()

#endif
//-----------------------------------------------------------------------------
/*---- End of file ----*/
