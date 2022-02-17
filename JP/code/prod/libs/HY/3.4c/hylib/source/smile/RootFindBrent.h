// -----------------------------------------------------------------------
// This proprietary software has been developed strictly for J.P. Morgan's
// internal use.  Any use or misuse, intentional or otherwise, which
// contradicts or places this policy in jeopardy is strictly forbidden.
//
// Copyright 1999 J.P. Morgan & Co. Incorporated.  All rights reserved.
// -----------------------------------------------------------------------
//
// 1.0 10/7/99 Afshin Bayrooti
// $Header$
//

#if ! defined(_DLIB_BU_BLOCKS_ROOT_FIND_BRENT_)
#define _DLIB_BU_BLOCKS_ROOT_FIND_BRENT_

//#include "General/RuntimeError.h"
//using namespace Dlib;
#include "kexception.h"

#include "BaseFunction.h"
//#include "systemh.h"
#include "rtbrent.h"


// Set up default values.
const double SOLVER_INITIAL_X_STEP  = 1E-8;
const int    SOLVER_INITIAL_F_DERIV = 0;
const double SOLVER_X_TOLERANCE     = 1E300;
const double SOLVER_F_TOLERANCE     = 1E-12;
const int    SOLVER_MAX_ITERATIONS  = 50;
const double SOLVER_LOWER_BOUND     = -0.999;
const double SOLVER_UPPER_BOUND     = 10000.00;

const double ROOTLO = 0;
const double ROOTHI = 100000;



int dummyRootFindFunction(double x, void  *data, double *result);

template<class Functor>
double RootFindBrent(
	Functor&  function,			 /*(I) Function to be solved */
	double          guess,               /*(I) Initial guess */
	double          boundLo,             /*(I) Lower bound of legal X */
	double          boundHi,             /*(I) Upper bound of legal X */
	int             numIterations,       /*(I) Maximum number of iterations */
	double          initialXStep,        /*(I) Size of step in x */
	double          initialFDeriv,       /*(I) Derivative, defaults to 0 */
	double          xacc,                /*(I) X accuracy tolerance */
	double          facc)                /*(I) Function accuracy tolerance */
{
	double result;
	ConcreteFunctorWrapper<Functor> funcWrapper(function);
	
	if(GtoRootFindBrent(dummyRootFindFunction,
					 &funcWrapper,
					 boundLo,
					 boundHi,
					 numIterations,
					 guess,
					 initialXStep,
					 initialFDeriv,
					 xacc,
					 facc,
					 &result)==FAILURE)
	{

	     KException( "root not found" );
	}
	return result;

}


template<class Functor>
double RootFindBrent(
	Functor&        function,			 /*(I) Function to be solved */
    double          guess,               /*(I) Initial guess */
	double          boundLo,             /*(I) Lower bound of legal X */
	double          boundHi)             /*(I) Upper bound of legal X */
{
    return RootFindBrent(function,
                         guess,
                         boundLo,
                         boundHi,
                         SOLVER_MAX_ITERATIONS,
                         SOLVER_INITIAL_X_STEP,
                         SOLVER_INITIAL_F_DERIV,
                         SOLVER_X_TOLERANCE,
                         SOLVER_F_TOLERANCE);              

}


#endif

