//----------------------------------------------------------------------------
//
//   Group       : Quantitative Research & Development
//
//   Filename    : SecantBrentRootFinder.cpp
//
//   Description : Port of ALIB root finder.
//
//   Author      : Richard Appleton
//
//   Date        : 10th June 2005
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/SecantBrentRootFinder.hpp"
#include "edginc/Atomic.hpp"


DRLIB_BEGIN_NAMESPACE


SecantBrentRootFinder::SecantBrentRootFinder(double pTolerance, int pMaxIterations)
: tolerance(pTolerance), maxIterations(pMaxIterations)
{
}


// ALIB: rtbrent.c#83
double SecantBrentRootFinder::solve(
    const Func1D::NoDeriv& fn, 
    const double           guess, 
    const double           lowerBound, 
    const double           upperBound,
    double                 initialXStep)
{
    // NB. not static as we want to know which iteration failed
    string method = string(typeid(fn).name()) + "::solve";
    static const double ONE_PERCENT = 0.01;

    try
    {
        // make sure lower bound is below higher bound
        if (lowerBound >= upperBound)
        {
            string msg = Format::toString(
                "Lower bound (%2.6e) >= higher bound (%2.6e)",
                lowerBound, upperBound);
            throw ModelException(method, msg);
        }

        // check to make sure the guess is within the bounds
        if (guess < lowerBound || guess > upperBound)
        {
            string msg = Format::toString(
                "Guess (%2.6e) is out of range [%2.6e, %2.6e]", 
                guess, lowerBound, upperBound);
            throw ModelException(method, msg);
        }

		double xPoints[3];
		double yPoints[3];

		// check if guess is the root
		xPoints[0] = guess;
		yPoints[0] = fn(xPoints[0]);
		if (Maths::areEqualWithinTol(0.0, yPoints[0], tolerance))
		{
			return xPoints[0];
		}

		// set initial step size
    	double boundSpread = (upperBound - lowerBound);
        if (Maths::isZero(initialXStep))
        {
            initialXStep = ONE_PERCENT * boundSpread;
        }

		// as derivative is not passed in take step of initial step size
		xPoints[2] = xPoints[0] + initialXStep;

		// now check xPoints[2] is within bounds -if not adjust it so that it is
		if (xPoints[2] < lowerBound || xPoints[2] > upperBound)
		{
			// switch step direction
			xPoints[2] = xPoints[0] - initialXStep;
			if (xPoints[2] < lowerBound)
			{
				// xPoints[2] is still too small, so make it lower bound
				xPoints[2] = lowerBound;
			}

			if (xPoints[2] > upperBound)
			{
				// xPoints[2] is too large, so make it upper bound
				xPoints[2] = upperBound;
			}

			if (xPoints[2] == xPoints[0])
			{
				// cannot have xPoints[0] and xPoints[2] be the same
				if (xPoints[2] == lowerBound)
				{
					xPoints[2] = lowerBound + ONE_PERCENT * boundSpread;
				}
				else
				{
					xPoints[2] = upperBound - ONE_PERCENT * boundSpread;
				}
			}
		}

		// finally, check function can return value at xPoints[2]
		yPoints[2] = fn(xPoints[2]);

		// ... check if new point meets tolerance requirements
		if (Maths::areEqualWithinTol(0.0, yPoints[2], tolerance))
		{
			return xPoints[2];
		}

		// call secant method to find root, or to get a 3rd point, so that 2 of
		// the 3 points bracket the root
		bool bracketed = false;
		bool foundIt = false;
		double result = secantMethod(fn, xPoints, yPoints, lowerBound, upperBound, foundIt, bracketed);
		if (foundIt)
		{
			return result;
		}
        else if (!bracketed)
		{
			// root was not bracketed, so try at the bounds
			double fLo = fn(lowerBound);
			if (Maths::areEqualWithinTol(0.0, fLo, tolerance))
			{
				return lowerBound;
			}

			// if these points bracket the root assign points so that 
			// xPoints[0] < xPoints[2]
			if (yPoints[0] * fLo < 0.0)
			{
				xPoints[2] = xPoints[0];
				xPoints[0] = lowerBound;
				yPoints[2] = yPoints[0];
				yPoints[0] = fLo;
			}
			else
			{
				// root is still not bracketed, so try at upper bound
				double fHi = fn(upperBound);
				if (Maths::areEqualWithinTol(0.0, fHi, tolerance))
				{
					return upperBound;
				}

				// if points bracket root assign xPoints[2] to upper bound
				if (yPoints[0] * fHi < 0.0)
				{
					xPoints[2] = upperBound;
					yPoints[2] = fHi;
				}
				else
				{
					// root could not be bracketed at the bounds
                    string msg = Format::toString(
						"Function values (%2.6e,%2.6e) at bounds (%2.6e, %2.6e) imply no root exists",
						fLo, fHi, lowerBound, upperBound);
					throw ModelException(method, msg);
				}
			}

			// xPoints[0] and xPoints[2] bracket the root, but we need a 3rd point
			// to do Brent method.  Take the midpoint.
			xPoints[1] = 0.5 * (xPoints[0] + xPoints[2]);
			yPoints[1] = fn(xPoints[1]);
			if (Maths::areEqualWithinTol(0.0, yPoints[1], tolerance))
			{
				return xPoints[1];
			}
		}

		return brentMethod(fn, xPoints, yPoints);
    }
    catch(exception& e)
    {
        throw ModelException(e, method, "failed to find root");
    }
}


void SWITCH(double& first, double& second)
{
	double tmp = first;
	first = second;
	second = tmp;
}

double SecantBrentRootFinder::brentMethod(
     const Func1D::NoDeriv& fn, 
     double                 *xPoints,      /* (I) Array of x values */
     double                 *yPoints)      /* (I) Array of y values */
{
	static const string method = "SecantBrentRootFinder::brentMethod";

    int         j;                      /* Index */
    double      ratio;                  /* (x3-x1)/(x2-x1) */
    double      x31;                    /* x3-x1*/
    double      x21;                    /* x2-x1*/
    double      f21;                    /* f2-f1 */
    double      f31;                    /* f3-f1 */
    double      f32;                    /* f3-f2 */
    double      xm;                     /* New point found using Brent method*/
    double      fm;                     /* f(xm) */

    double      x1 = xPoints[0];        /* xN short hand for xPoints[n] */
    double      x2 = xPoints[1];
    double      x3 = xPoints[2];

    double      f1 = yPoints[0];
    double      f2 = yPoints[1];
    double      f3 = yPoints[2];

    for (j=1; j<=maxIterations; j++) 
    {
        /* Always want to be sure that f1 and f2 have opposite signs,
         * and f2 and f3 have the same sign.
         */
        if (f2*f1>0.0)
        {   
            SWITCH(x1,x3);
            SWITCH(f1,f3);
        }
        f21 = f2-f1;
        f32 = f3-f2;
        f31 = f3-f1;
        x21 = x2-x1;
        x31 = x3-x1;
        /* Check whether this is suitable for interpolation. When checking
         * for f21,f31,f32 = 0, we don't use IS_ALMOST_ZERO for efficiency
         * reasons. If the objective function has been truncated to a 
         * certain number of digits of accuracy, f21,f31,or f32 could be
         * (were in one case) zero. In this case we need to protect against
         * division by zero. So we use bisection instead of brent.
         */
        ratio = (x3-x1)/(x2-x1);
        if (f3*f31<ratio*f2*f21 || f21 == 0. || f31 == 0. || f32 == 0.) 
        {
            /* This is not suitable, do bisection 
             */
            x3 = x2;
            f3 = f2; 

        }
        else 
        {
            xm = x1 - (f1/f21)*x21 + ((f1*f2)/(f31*f32))*x31 - 
                ((f1*f2)/(f21*f32))*x21;
			fm = fn(xm);
			if (Maths::areEqualWithinTol(0.0, fm, tolerance))
            {
                return xm;
            }
            /* If the new point and point1 bracket the root,
               replace point3 with the new point */
            if (fm*f1<0.0)
            {
                x3=xm;
                f3=fm;
            }
            /* If the root is not bracketed, replace point1 with new point,
               and point3 with point2 */
            else
            {
                x1=xm;
                f1=fm;
                x3=x2;
                f3=f2;
            }
        }             
        x2 = 0.5*(x1+x3);
		f2 = fn(x2);
		if (Maths::areEqualWithinTol(0.0, f2, tolerance))
        {
            return x2;
        }
    }

    throw ModelException(method, "Maximum number of iterations exceeded.");
}    
       

#ifndef ABS
#define ABS(x) ((x) > 0 ? (x) : -(x))
#endif


double SecantBrentRootFinder::secantMethod(
     const Func1D::NoDeriv& fn, 
     double*                xPoints,       // (I/O) Array of x points
     double*                yPoints,       // (I/O) Array of y points
     double                 boundLo,       // (I) Lower bound
     double                 boundHi,       // (I) Upper bound
     bool&                  foundIt,       // (O) If solution was found
     bool&                  bracketed)     // (O) if root was bracketed
{
	static const string method = "SecantBrentRootFinder::secantMethod";

    int           j=maxIterations;      /* Index */
    double        dx;                   /* Delta x used for secant */       

    foundIt = false;           // Until solution is found
    bracketed = false;         // Until bracketed

    while (j--)
    {
        /* Swap points so that yPoints[0] is smaller than yPoints[2] */
        if (ABS(yPoints[0]) > ABS(yPoints[2]))
        {   
            SWITCH(xPoints[0],xPoints[2]);
            SWITCH(yPoints[0],yPoints[2]);
        }

        /* Make sure that you do not divide by a very small value */
        if (ABS(yPoints[0]-yPoints[2]) <= tolerance)
        {
            if (yPoints[0] - yPoints[2] > 0)
            {
                dx = -yPoints[0]*(xPoints[0] - xPoints[2])/tolerance;
            }
            else
            {
                dx = yPoints[0]*(xPoints[0] - xPoints[2])/tolerance;
            }
        }
        else
        {
            dx= (xPoints[2]-xPoints[0])* yPoints[0]/(yPoints[0]-yPoints[2]);
        }
        xPoints[1] = xPoints[0] + dx;

        /* Make sure that the point is within bounds 
         */
        if (xPoints[1] < boundLo || xPoints[1] > boundHi)
        {
            return 0.0;     // Not bracketed, not found
        }

		yPoints[1] = fn(xPoints[1]);
		if (Maths::areEqualWithinTol(0.0, yPoints[1], tolerance))
        {
			foundIt = true;
			bracketed = true;
            return xPoints[1];
        }

        if ((yPoints[0] < 0 && yPoints[1] < 0 && yPoints[2] < 0) ||
            (yPoints[0] > 0 && yPoints[1] > 0 && yPoints[2] > 0))
        {
            /* Swap points so that yPoints[0] is always smallest 
             */
            if (ABS(yPoints[0]) > ABS(yPoints[1]))
            {
                xPoints[2] = xPoints[0];
                yPoints[2] = yPoints[0];
                xPoints[0] = xPoints[1];
                yPoints[0] = yPoints[1];
            }
            else 
            {
                xPoints[2] = xPoints[1];
                yPoints[2] = yPoints[1];
            }
            continue;
        }
        else
        { 
            /* Root was bracketed. 
             * Swap points so that yPoints[0]*yPoints[2] < 0 
             */
            if (yPoints[0]*yPoints[2] > 0)
            {
                /* Make sure that you swap so that 
                 * xPoints[0]<xPoints[1]<xPoints[2] 
                 */
                if (xPoints[1] < xPoints[0])
                {
                    SWITCH(xPoints[0], xPoints[1]);
                    SWITCH(yPoints[0], yPoints[1]);
                }
                else
                {
                    SWITCH (xPoints[1], xPoints[2]);
                    SWITCH (yPoints[1], yPoints[2]);
                }
            }
            /* Root was bracketed, but not found.
             */
            bracketed = true;
            return 0.0;
        }
    } /* while */


    /* Root not bracketed or found.
     */
    return 0.0;
}


#undef ABS


DRLIB_END_NAMESPACE
