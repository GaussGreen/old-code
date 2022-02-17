/*!
 *
 * Copyright (c) CDC IXIS CM November 2003 Paris
 *
 *	\file numericconstant.cpp
 *  \brief 
 *	\author  E. Benhamou
 *	\version 1.0
 *	\date March 2003
 */


#include "gpbase/numericconstant.h"
#include <cmath>


CC_BEGIN_NAMESPACE( ARM )

const double ARM_NumericConstants::ARM_ONE_THIRD				= 0.3333333333333333;
const double ARM_NumericConstants::ARM_PI						= 3.1415926535897932;
const double ARM_NumericConstants::ARM_2_PI						= 2*ARM_NumericConstants::ARM_PI;
const double ARM_NumericConstants::ARM_PI_BY_4					= 0.7853981633974483;
const double ARM_NumericConstants::ARM_PI_BY_2					= 0.5*ARM_NumericConstants::ARM_PI;
const double ARM_NumericConstants::ARM_SQRT_PI					= 1.772453850905516;
const double ARM_NumericConstants::ARM_SQRT_2					= 1.4142135623730950;
const double ARM_NumericConstants::ARM_SQRT_3					= sqrt(3.0);
const double ARM_NumericConstants::ARM_TWO_DIVIDED_BY_SQRT_3	= 2.0/sqrt(3.0);
const double ARM_NumericConstants::ARM_SQRT_2_DIVIDED_BY_SQRT_3	= sqrt(2.0/3.0);
const double ARM_NumericConstants::ARM_SQRT_2_DIVIDED_BY_SQRT_PI= sqrt(2.0/3.1415926535897932);
const double ARM_NumericConstants::ARM_SQRT_2_PI				= 2.506628274631001;
const double ARM_NumericConstants::ARM_SQRT_8_PI				= 5.013256549262001;
const double ARM_NumericConstants::ARM_INVSQRTPI				= 0.5641895835477563;
const double ARM_NumericConstants::ARM_INVSQRT2PI				= 0.398942280401433;
const double ARM_NumericConstants::ARM_INVSQRTSQRTPI			= 0.7511255444649424;

const double ARM_NumericConstants::ARM_HALF_LOG_2PI				= 0.91893853320467274;
const double ARM_NumericConstants::ARM_LOG_2PI					= 1.83787706640934548;
const double ARM_NumericConstants::ARM_HALF_LOG_PI				= 0.57236494292470008;

const double ARM_NumericConstants::ARM_INVERSE_E				= exp(-1.);

const double ARM_NumericConstants::ARM_BIGGEST_POSITIVE_NUMBER	    = 1e+308;
const double ARM_NumericConstants::ARM_LOWEST_POSITIVE_NUMBER	    = 1e-308;
const double ARM_NumericConstants::ARM_SQRT_LOWEST_POSITIVE_NUMBER	= 1e-154;
const double ARM_NumericConstants::ARM_INFINITY				        = 1.0e+66;
const double ARM_NumericConstants::ARM_TOLERENCE				    = 1.0e-12;
const double ARM_NumericConstants::ARM_STD_DOUBLE_TOLERANCE		    = 1.0e-14;
const double ARM_NumericConstants::ARM_DOUBLE_TOLERENCE			    = 1.0e-17;
const size_t ARM_NumericConstants::ARM_GP_MAX_ITER                  = 100;
const double ARM_NumericConstants::ARM_GP_GRAD_TOL                  = 1.0e-12;
CC_END_NAMESPACE()

///---------------------------------------------------------------------------
///---- End of file ----