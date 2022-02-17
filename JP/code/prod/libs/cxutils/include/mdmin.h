/*
***************************************************************************
** FILENAME: mdmin.h
**
** Multi-dimensional minimization routine - based on Powell from NumRec
***************************************************************************
*/

#ifndef _CX_MDMIN_H
#define _CX_MDMIN_H

#include "cxutils.h"

#include <alib/gtomat.h>

#ifdef __cplusplus
extern "C"
{
#endif

typedef int (*CxTMultiDimObjFunc) (int n, double *x, void* para, double *f);

/**
***************************************************************************
** Performs multi-dimensional minimization without constraints.
**
** If you just provide a function evaluator, then it will use Powell's
** method.
**
** With a derivative function as well then it will use the Fletcher-Reeves-
** Polak-Ribiere method, but this has not been implemented yet. The idea
** is that the function prototype will not need to change if we decide
** to implement that as well.
**
** The idea is that the minimization has a state, and you call the function
** with the current state and then at the end of the call it has a new
** state which you can use as the input for other function calls.
**
** The state contains the values of the variables you are optimizing, the
** total number of iterations, the current target value and the latest
** degree of convergence. In addition it contains a matrix of search
** directions which is useful for ensuring that you start at the same
** point when you re-start the optimization.
**
** When you initialize the state, you only need to set the initial guess
** of the variables you are optimizing.
**
** This can be quite a slow function if you have a large number of
** iterations or small relative tolerance, which is why we provide both
** the ability to call the function in a loop (via the state variables)
** and the ability to watch the progress via a logfile.
***************************************************************************
*/
CxTMultiDimMinState* CxMultiDimMinimization(
    /** Prior state of the minimization. Not modified by this call - the
        returned state is a separate instance. */
    CxTMultiDimMinState *initState,
    /** Function which computes the single valued f(x) */
    CxTMultiDimObjFunc   func,
    /** Function which computes the vector valued df/dx. Optional. */
    CxTMultiDimObjFunc   dfunc,
    /** Data required by func and dfunc for evaluating f and df/dx. */
    void                *data,
    /** Relative tolerance - when the minimum value does not change by more
        than the relative tolerance, then the iteration will stop. */
    double               ftol,
    /** Maximum number of iterations at this step. */
    int                  maxiter,
    /** Maximum time allowed (in seconds). Time=0 => no limit on the time. */
    double               maxtime,
    /** Optional logfilename for recording progress. NULL or "" if you
        don't need want to write the progress. The file is started again
        for each function call. */
    char                *logfilename);

int CxMultiDimMinStateValidate (CxTMultiDimMinState *state);

#ifdef __cplusplus
}
#endif

#endif

