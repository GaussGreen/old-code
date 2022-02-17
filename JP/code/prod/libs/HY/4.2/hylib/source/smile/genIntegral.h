// -----------------------------------------------------------------------
// This proprietary software has been developed strictly for J.P. Morgan's
// internal use.  Any use or misuse, intentional or otherwise, which
// contradicts or places this policy in jeopardy is strictly forbidden.
//
// Copyright 1999 J.P. Morgan & Co. Incorporated.  All rights reserved.
// -----------------------------------------------------------------------
//
// 1.0 10/19/99 Neil Yang
// 
//

#if ! defined _GEN_INTEGRAL_
#define _GEN_INTEGRAL_

#include "BaseFunction.h"
#include "integral.h"

#define MAX_INTEGRAL_ERROR 1.0E-6


int dummyFunction(double x, void  *data, double *result);




template<class Functor>
double Integral(
	Functor&  function,          /* (I) Function te be integrated */
	double    boundLo,           /* (I) Lower integral limit */
	double    boundHi,           /* (I) Higher integral limit */
	double    maxError)          /* (I) Maximum error */
{
	double result;
	ConcreteFunctorWrapper<Functor> funcWrapper(function);
	
	GtoIntegral(dummyFunction,
				boundLo,
				boundHi,
				&funcWrapper,
				maxError,
				&result);
	return result;

}

template<class Functor>
double Integral(
	Functor&  function,          /* (I) Function te be integrated */
	double    boundLo,           /* (I) Lower integral limit */
	double    boundHi)           /* (I) Higher integral limit */
{
	return Integral(function,
					boundLo,
					boundHi,
					MAX_INTEGRAL_ERROR);
}

#endif


