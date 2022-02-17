/************************************************************************
 * Module:	DRL - MATH
 * Function:	Mathematical Functions
 * Function:	
 * Author:	J. Chislenko
 * Revision:	$Header$
 ************************************************************************/
#include "drlstd.h"		/* platform compatibility */
#include "drlerr.h"		/* DrlErrMsg */
#include <math.h>

#include "drlmath.h"         

/*f-----------------------------------------------------
 * Numerical integration : 1D Simpson (from data points).
 * 
 * <br><br>
 * Integrate function given by its values
 * using extended Simpson's rule
 * For less than 8 points use central Riemann sums
 * For a single point use value * step
 */

double  DrlSimpsIntegral(
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
    else if(numPoints >= 2)
    {
	numEndPts = 1;
	intValue = 
	    (fValues[0] + fValues[numPoints-1]) * 0.5;
    }
    else 
    {
	numEndPts = 0;
    }
	
    for (i=2*numEndPts, fValue=fValues+numEndPts; i<numPoints; i++, fValue++)
    {
	intValue += *fValue;
    }

    intValue *= step;

    return (intValue);
}
	

/*f-----------------------------------------------------
 * Math : smooth Heaviside function.
 * 
 * <br><br>
 * Returns smoothed step function value of x.
 *                  Heaviside function is defined as:
 *                  F(X) = 1 if X >= 0.
 *                  F(X) = 0 if X <  0.
 */
double  DrlSmoothStepFcn(
     double  x,                         /* (I) Argument */
     double  step)                     /* (I) 1/2 step */
{
    double smoothValue, xSquare;

    if (x <= -step)                  /* If both x & step = 0, return 0 */
    {
        smoothValue = 0.;
    }
    else if (x >= step ||
             IS_ALMOST_ZERO(step))
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
}

/*f-----------------------------------------------------
 * Math : smooth max function.
 * 
 * <br><br>
 * DrlSmoothMAX(): Returns smoothed MAX(x,0).
 *                 Is given by the integral of the smooth step fcn.
 */
double  DrlSmoothMAX(
     double  x,                         /* (I) Argument */
     double  step)                      /* (I) 1/2 step */
{
    double smoothValue, xSquare;

    if (x <= -step)                  /* If both x & step = 0, return 0 */
    {
        smoothValue = 0.;
    }
    else if (x >= step ||
             IS_ALMOST_ZERO(step))
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
}

