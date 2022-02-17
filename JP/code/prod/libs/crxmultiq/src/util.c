/******************************************************************************
 * Module:      Q3
 * Submodule:
 * File: util.c 
 * Function:    
 * Author:      Interest Rates DR
 * Revision:    $Header$
 *****************************************************************************/
#include <math.h>
#include <stdio.h>

#include "crmultiq.h"


/*f----------------------------------------------------------------------------
 * DrlSimpsIntegral    
 * 
 * Numerical integration : 1D Simpson (from data points).
 * <br><br>
 * Integrate function given by its values
 * using extended Simpson's rule
 * For less than 8 points use central Riemann sums
 * For a single point use value * step
 */

double  Q3SimpsIntegral
(
    double            step,             /* (I) */
    long              numPoints,        /* (I) */
    double           *fValues)          /* (I) */ 
{
    double c1 = 17./48.;
    double c2 = 59./48.;
    double c3 = 43./48.;
    double c4 = 49./48.;

    long i;
    long numEndPts; /*  for each end of the interval */
    double  intValue = 0.;
    double *fValue = NULL;

    if(numPoints >= 8) 
    {
        numEndPts = 4;
        intValue =
            (fValues[0] + fValues[numPoints-1]) * c1 +
            (fValues[1] + fValues[numPoints-2]) * c2 +
            (fValues[2] + fValues[numPoints-3]) * c3 +
            (fValues[3] + fValues[numPoints-4]) * c4;
    } 
    else if (numPoints >= 2) 
    {
        numEndPts = 1;
        intValue = 
            (fValues[0] + fValues[numPoints-1]) * 0.5;
    } 
    else 
    {
        numEndPts = 0;
    }
        
    for (i=2*numEndPts,fValue=fValues+numEndPts; i<numPoints;
         i++,fValue++) 
    {
        intValue += *fValue;
    }

    intValue *= step;

    return (intValue);
} /* Q3SimpsIntegral */
        
/*f----------------------------------------------------------------------------
 * Math : smooth Heaviside function.
 * 
 * <br><br>
 * Returns smoothed step function value of x.
 *                  Heaviside function is defined as:
 *                  F(X) = 1 if X >= 0.
 *                  F(X) = 0 if X <  0.
 */
double  Q3SmoothStepFcn(
    double  x,                         /* (I) Argument */
    double  step)                      /* (I) 1/2 step */
{
    double smoothValue, xSquare;

    if (x <= -step) 
    {                 /* If both x & step = 0, return 0 */ 
        smoothValue = 0.;
    }
    else if (x >= step || fabs(step) < TINY) 
    {
        smoothValue = 1.;        
    } 
    else 
    {
        x /=step;
        xSquare = x * x;
        smoothValue = ((3.*xSquare-10.) * xSquare + 15.) * x *.0625 +.5;
    } 

    return (smoothValue);
} /* Q3SmoothStepFcn */

/*f----------------------------------------------------------------------------
 * Math : smooth max function.
 * 
 * <br><br>
 * Q3SmoothMAX(): Returns smoothed MAX(x,0).
 *                 Is given by the integral of the smooth step fcn.
 */
double  Q3SmoothMAX(
    double  x,                         /* (I) Argument */
    double  step)                      /* (I) 1/2 step */
{
    double smoothValue, xSquare;

    if (x <= -step) 
    {                 /* If both x & step = 0, return 0 */
        smoothValue = 0.;
    } 
    else if (x >= step || fabs(step) < TINY) 
    {
        smoothValue = x;        
    } 
    else 
    {
        x /=step;
        xSquare = x * x;
        smoothValue = ((xSquare-5.) * xSquare + 15.) * xSquare * 0.03125 
            + 0.5 * x + 0.15625;
        smoothValue *= step;
    } 
    return (smoothValue);
} /* Q3SmoothMAX */


/*-----------------------------------------------------------------------------
 * Q3ProbInBarriers
 *
 * Returns probability that loBar < x < hiBar, where the barriers are smoothed
 * by the epsilon parameters loEps and hiEps respectively.
 *
 */
double Q3ProbInBarriers(
    double x,     /* (I) input value                      */
    double loBar, /* (I) low barrier                      */
    double hiBar, /* (I) high barrier                     */
    double loEps, /* (I) low barrier smoothing parameter  */
    double hiEps, /* (I) high barrier smoothing parameter */
    double smooth /* (I) smoothing step size              */
    )
{

    double prob = 0;
    double probLo, probHi;

    /* avoid warnings about unused inputs */
    smooth=smooth;

    /* outside with certainty */
    if ((x < (loBar - loEps)) || (x > (hiBar + hiEps))) goto RETURN;

    /* prob above low barrier */
    if (fabs(loEps) < TINY)
    {
	if (x > loBar)
	{
	    probLo = 1.;
	}
	else
	{
	    probLo = 0.;
	}
    }
    else
    {
	probLo = (MAX(x - (loBar - loEps), 0.) - 
		  MAX(x - loBar, 0.))/loEps;
    }
    
    /* prob above high barrier */
    if (fabs(hiEps) < TINY)
    {
	if (x > hiBar)
	{
	    probHi = 1.;
	}
	else
	{
	    probHi = 0.;
	}
    }
    else
    {
	probHi = (MAX(x - hiBar, 0.) - 
		  MAX(x - (hiBar + hiEps), 0.))/hiEps;
    }


    /* prob within barriers */
    prob = probLo - probHi;


 RETURN:

    return prob;

}


/*-----------------------------------------------------------------------------
 * Q3ProbOutBarriers
 *
 * Returns probability that x < loBar && x > hiBar, where the barriers are smoothed
 * by the epsilon parameters loEps and hiEps respectively.
 *
 */
double Q3ProbOutBarriers(
    double x,     /* (I) input value                      */
    double loBar, /* (I) low barrier                      */
    double hiBar, /* (I) high barrier                     */
    double loEps, /* (I) low barrier smoothing parameter  */
    double hiEps, /* (I) high barrier smoothing parameter */
    double smooth /* (I) smoothing step size              */
    )
{

    double prob;

    /* prob inside barriers */
    prob = Q3ProbInBarriers( 
	x,
	loBar,
	hiBar,
	loEps,
	hiEps,
	smooth);
    
    prob = 1. - prob;

    return prob;

}

/*---------------------------------------------------------------------------
 * Q3ProbInBarriersSmooth
 *
 * Returns probability that loBar < x < hiBar, where the barriers are 
 * smoothed by both the smooth window and epsilons.
 *
 */
double Q3ProbInBarriersSmooth(
    double x,     /* (I) input value                      */
    double loBar, /* (I) low barrier                      */
    double hiBar, /* (I) high barrier                     */
    double loEps, /* (I) low barrier smoothing parameter  */
    double hiEps, /* (I) high barrier smoothing parameter */
    double smooth /* (I) smoothing step size              */
    )
{

    double prob = 0;
    double probLo, probHi;

    /* prob above low barrier */
    if (fabs(loEps) < TINY)
    {
        probLo = Q3SmoothStepFcn(x-loBar,smooth);
    }
    else
    {
        probLo = (Q3SmoothMAX(x - loBar + loEps, smooth) -
                  Q3SmoothMAX(x - loBar, smooth)) / loEps;
    }
    
    /* prob above high barrier */
    if (fabs(hiEps) < TINY)
    {
        probHi = Q3SmoothStepFcn(x-hiBar,smooth);
    }
    else
    {
        probHi = (Q3SmoothMAX(x - hiBar, smooth) -
                  Q3SmoothMAX(x - hiBar - hiEps, smooth)) / hiEps;
    }


    /* prob within barriers */
    prob = probLo - probHi;

    return prob;

}

/*------------------------------------------------------------------------
 * Q3ProbOutBarriersSmooth
 *
 * Returns probability that x < loBar && x > hiBar, where the barriers are 
 * smoothed by both the smooth window and epsilons.
 *
 */
double Q3ProbOutBarriersSmooth(
    double x,     /* (I) input value                      */
    double loBar, /* (I) low barrier                      */
    double hiBar, /* (I) high barrier                     */
    double loEps, /* (I) low barrier smoothing parameter  */
    double hiEps, /* (I) high barrier smoothing parameter */
    double smooth /* (I) smoothing step size              */
    )
{

    double prob;

    /* prob inside barriers */
    prob = Q3ProbInBarriersSmooth( 
    x,
    loBar,
    hiBar,
    loEps,
    hiEps,
    smooth);
    
    prob = 1. - prob;

    return prob;
}


/* --------------------------------------------------------------------------
 * Swop two elements 
 */

int Q3SwapDouble(double *x, double *y)
{
    double Tmp;
    Tmp = *x;
    *x  = *y;
    *y  = Tmp;
    return SUCCESS;
}

int Q3SwapInt(int *x, int *y)
{
    int Tmp;
    Tmp = *x;
    *x  = *y;
    *y  = Tmp;
    return SUCCESS;
}

