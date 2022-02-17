//----------------------------------------------------------------------------
//
//   Group       : Core Analytics Team
//
//   Filename    : ZeroSysEqn.cpp
//
//   Description : Wrapper around IMSL_ZERO_SYS_EQN
//
//   Author      : Richard Appleton
//
//   Date        : 19th June 2006
//
//----------------------------------------------------------------------------
#include "edginc/config.hpp"
#include "edginc/ZeroSysEqn.hpp"
#include "edginc/IMSLException.hpp"


DRLIB_BEGIN_NAMESPACE


// callback mechanism is not thread safe
static ZeroSysEqnFunc* func;
static ModelException* error;


static void imsl_callback(int n, double* x, double* errors)
{
    try
    {
        (*func)(n, x, errors);
    }
    catch(exception& e)
    {
        // get iteration to stop
        for (int i = 0 ; i < n ; i++)
            errors[i] = 0.0;

        // record error
        error = new ModelException(e);
    }
}


ZeroSysEqnSolver::ZeroSysEqnSolver(double pTolerance, int pMaxIterations)
: tolerance(pTolerance), maxIterations(pMaxIterations)
{
}


ZeroSysEqnSolver::~ZeroSysEqnSolver()
{
}


DoubleArraySP ZeroSysEqnSolver::solve(ZeroSysEqnFunc& f, DoubleArray& guess) // throws IMSLException
{
    // setup callback
    func = &f;

    // create result array
    int n = guess.size();
	DoubleArraySP outx(new DoubleArray(n));

    // invoke IMSL function
    error = NULL;
    imsl_d_zeros_sys_eqn(imsl_callback,
                         n,
                         IMSL_XGUESS , &guess[0],
                         IMSL_ERR_REL, tolerance,
                         IMSL_MAX_ITN, maxIterations,
                         IMSL_RETURN_USER, &(*outx)[0],
                         0 /* terminates imsl vararg list */);

    if (error)
    {
        ModelException x(*error);
        delete error;
        error = NULL;
        throw x;
    }

	// throw IMSLException if IMSL function fails
	IMSLException::throwIMSLExceptionIfError();

    return outx;
}


DRLIB_END_NAMESPACE
