/* --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---  --
Function: irxRootFindBrent

created by: 4/1/96 Julia Tolpin

description: Defines the new root finder, irxRootFindBrent

inputs: 
--- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---  -- 
Proprietary software, whether developed for Morgan by in-house
staff or external contractors, must not be copied or removed from Morgan
premises without the approval of a Senior Vice President and Audit.

This proprietary software has been developed strictly for Morgan's own
internal use.  Any use or misuse, intentional or otherwise which contradicts
or places this policy in jeopardy is strictly forbidden.  Also, any actions or 
inactions resulting in a breach of this goal will be dealt with harshly.

Do Not Distribute this to anyone outside the Analytics Library
Group without authorization from its management. 

Copyright 1995 J.P. Morgan & Co. Incorporated.   All rights reserved.
-------------------------------------------------------------------------  */

#ifndef IRX_RTBRENT_H
#define IRX_RTBRENT_H

#include "cgeneral.h"

#ifdef __cplusplus
extern "C"
{
#endif

/** Objective function. Has one variable input (x) which is changed by
    the root finder. Should return in f the value of f(x). Parameters 
    para are passed through from the call to the root finder unchanged
    into this function. */
typedef int (*IrxTObjectFunc) (double x, void * para, double *f);

/**
 * Finds the root of f(x) =  0 using a combination of 
 * secant,bisection and an inverse quadratic interpolation method.
 *
 * Failure can occur if the objective function fails at any point, or
 * if a solution cannot be found to the desired accuracy within the maximum
 * number of iterations.
 */
extern int irxRootFindBrent(
   /** Function to be solved. */
   IrxTObjectFunc funcd,           
   /** Passed unchanged to the objective function (as the second parameter) */ 
   void           *data,         
   /** Lower bound on legal values of x */
   double         boundLo,   
   /** Upper bound on legal values of x */
   double         boundHi,  
   /** Maximum number of iterations */
   int            numIterations, 
   /** Initial guess */
   double         guess,        
   /** Initial step size in x. If set to zero, then will use one percent of
       (boundHi-boundLo). */
   double         initialXStep, 
   /** Initial derivative if known. If not zero, then the second guess point
       is (guess - f(guess) / initialFDeriv   */
   double         initialFDeriv,
   /** X-accuracy tolerance. The solution will be less than this value from
       the previously tested value. If you are not interested in x-accuracy
       then set this to a high value. */
   double         xacc,        
   /** Function accuracy tolerance. The function value for this solution will
       be less than (in absolute terms) this accuracy. If you are not
       interested in function accuracy, then set this to a high value.
 
       It would not make sense to have both xacc and facc as high values. */
   double         facc,       
   /** Output. The solution is returned. */
   double         *solution);

#ifdef __cplusplus
}
#endif

#endif    /* IRX_RTBRENT_H */
