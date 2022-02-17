//----------------------------------------------------------------------------
//
//   Group       : Convertibles DR
//
//   Filename    : FDUtils.cpp
//
//   Description : Some utility functions useful for FD. Mainly interpolation.
//
//   Author      : Tycho von Rosenvinge
//
//   Date        : June 24, 2002
//
//----------------------------------------------------------------------------
#include "edginc/config.hpp" // needed for precompiled headers
#include "edginc/FDUtils.hpp"
#include "edginc/Maths.hpp"
#include "edginc/UtilFuncs.hpp"
//#include "C:/fldll194z/headers/nagmk19.hpp"
#include "edginc/Addin.hpp"
#ifdef _MSC_VER
#include <malloc.h>
#else
#include <alloca.h>
#endif

DRLIB_BEGIN_NAMESPACE

#define FDMAXSIZE 10000



static int FDGetSplineCoefficients(
    int        size,            /* (I) size of array */
    double *xvec,            /* (I) array of x values */
    double *yvec,            /* (I) array of y=f(x) values */
    double *fpp)            /* (O) array of y's second derivatives */
/*-------------------------------------------------------------------
** FUNCTION:    FDGetSplineCoefficients
**
** AUTHOR:        Milan Kovacevic, February 24, 1998
**
** DESCRIPTION:    Cubic spline interpolator (derived from
**                SplineCoefficients11)
**
** RETURNS:        Array of second derivatives
**
** AMENDMENTS:
**-----------------------------------------------------------------*/
{
    int ii6;
    double *uu,sig,pp;
    double ubuf[FDMAXSIZE];

    /* Key to program variables: */
    /* errorBuffer: place to store error messages */
    /* fpp: second derivative array */
    /* ii6: index variable for xi */
    /* size: number of grid cells for xi */
    /* yvec: yvec */
    /* xvec: solution variable */
    /* uu: work space */
    /* pp, sig: temporary */

    // MK: Use memory on stack for performance reasons if the array is not too big
    if (size <= FDMAXSIZE)
        uu = ubuf;
    else
    {
            uu = (double *)calloc(size, sizeof(double));
            if (uu == NULL)
                return FAILURE;
    }
    /* Calculate the 2nd derivatives. */
    /* The natural end condition at the left. */
    fpp[0] = 0.;
    uu[0] = 0.;
    /* The forward sweep of tridiagonal solver. */
    for (ii6=1; ii6<size - 1; ii6++) {
        sig = (xvec[ii6] - xvec[ii6 - 1]) / (xvec[ii6 + 1] - xvec[ii6 - 1]);
        pp = sig*fpp[ii6 - 1] + 2;
        fpp[ii6] = (sig - 1) / pp;
        uu[ii6] = ((6.*((yvec[ii6 - 1] - yvec[ii6]) / (xvec[ii6] - xvec[ii6 - 1
           ]) + (yvec[ii6 + 1] - yvec[ii6]) / (xvec[ii6 + 1] - xvec[ii6]))) / (xvec
           [ii6 + 1] - xvec[ii6 - 1]) - sig*uu[ii6 - 1]) / pp;
    }
    /* The natural end condition at the right. */
    fpp[size-1] = 0.;
    /* The backward sweep of tridiagonal solver. */
    for (ii6=size - 2; ii6>=0; ii6--) {
        fpp[ii6] = fpp[ii6]*fpp[ii6 + 1] + uu[ii6];
    }
    if (size > FDMAXSIZE)
        ::free(uu);
    return SUCCESS;
}

double FDInterpolateLinear(
    int        tabsize,    /* (I) size of the table to interpolate from */
    double*    xtab,        /* (I) x values of the table to interpolate from */
    double*    ytab,        /* (I) y values of the table to interpolate from */
    double    x            /* (I) x value at which to interpolate */
    )
/*-------------------------------------------------------------------
** FUNCTION:    FDInterpolateLinear
**
** AUTHOR:        Milan Kovacevic, February 24, 1998
**
** DESCRIPTION:    Linear interpolator (derived from
**                InterpolateLinear2)
**
** RETURNS:        Interpolated value at x
**
** AMENDMENTS:
**-----------------------------------------------------------------*/
{
    double    tt1;    /* tt1: fraction of interval */
    int        loc1;    /* loc1: raw interpolation index */
    int        inew1;    /* inew1: adjusted interpolation index */
    int        imid;    /* imid: average of ilow and ihigh */
    int        ilow;    /* ilow: index below interpolation index */
    int        ihigh;    /* ihigh: index above the interpolation index */

    if (tabsize == 1)
        return ytab[0];
    else if (tabsize == 0)
        return 0;
    /* Find the index of the data point. */
    ilow = -1;
    ihigh = tabsize;
    while (ihigh - ilow>1)
    {
        imid = (int)(floor((ihigh + ilow) / 2.0));
        if (x>xtab[imid]&&xtab[tabsize-1]<=xtab[0]||
            xtab[tabsize-1]>xtab[0]&&x<=xtab[imid])
        {
            ihigh = imid;
        }
        else
        {
            ilow = imid;
        }
    }
    loc1 = ilow;
    /* Adjust the indices of the end points. */
    if (loc1<=-1)
    {
        inew1 = 0;
    }
    else
    {
        if (loc1>=tabsize-1)
        {
            inew1 = tabsize - 2;
        }
        else
        {
            inew1 = loc1;
        }
    }
    /* Compute the fraction of the interval. */
    tt1 = (x - xtab[inew1]) / (xtab[inew1 + 1] - xtab[inew1]);
    return ytab[inew1] + tt1*(ytab[inew1 + 1] - ytab[inew1]);
}

double FDInterpolateLinearD(
    int        tabsize,    /* (I) size of the table to interpolate from */
    double*    xtab,        /* (I) x values of the table to interpolate from */
    double*    ytab,        /* (I) y values of the table to interpolate from */
    double    x,            /* (I) x value at which to interpolate */
    double*    y1            /* (O) first derivative */
    )
/*-------------------------------------------------------------------
** FUNCTION:    FDInterpolateLinear
**
** AUTHOR:        Milan Kovacevic, February 24, 1998
**
** DESCRIPTION:    Linear interpolator (derived from
**                InterpolateLinear2)
**
** RETURNS:        Interpolated value at x
**
** AMENDMENTS:
**-----------------------------------------------------------------*/
{
    double    tt1;    /* tt1: fraction of interval */
    int        loc1;    /* loc1: raw interpolation index */
    int        inew1;    /* inew1: adjusted interpolation index */
    int        imid;    /* imid: average of ilow and ihigh */
    int        ilow;    /* ilow: index below interpolation index */
    int        ihigh;    /* ihigh: index above the interpolation index */

    if (tabsize == 1)
        return ytab[0];
    else if (tabsize == 0)
        return 0;
    /* Find the index of the data point. */
    ilow = -1;
    ihigh = tabsize;
    while (ihigh - ilow>1)
    {
        imid = (int)(floor((ihigh + ilow) / 2.0));
        if (x>xtab[imid]&&xtab[tabsize-1]<=xtab[0]||
            xtab[tabsize-1]>xtab[0]&&x<=xtab[imid])
        {
            ihigh = imid;
        }
        else
        {
            ilow = imid;
        }
    }
    loc1 = ilow;
    /* Adjust the indices of the end points. */
    if (loc1<=-1)
    {
        inew1 = 0;
    }
    else
    {
        if (loc1>=tabsize-1)
        {
            inew1 = tabsize - 2;
        }
        else
        {
            inew1 = loc1;
        }
    }
    /* Compute the fraction of the interval. */
    tt1 = (x - xtab[inew1]) / (xtab[inew1 + 1] - xtab[inew1]);
    *y1 = (ytab[inew1 + 1] - ytab[inew1])/(xtab[inew1 + 1] - xtab[inew1]);
    return ytab[inew1] + tt1*(ytab[inew1 + 1] - ytab[inew1]);
}


int FDInterpolation(
    int        nOld,        /* (I) size of the table to interpolate from*/
    double*    xOld,        /* (I) x values of the table to interpolate from */
    double*    yOld,        /* (I) y values of the table to interpolate from */
    int        n,            /* (I) size of the output array */
    double*    x,        /* (I) array of x values at which to interpolate */
    double*    y)        /* (O) array of interpolated values for xvec*/
/*-------------------------------------------------------------------
** FUNCTION:    FDInterpolation
**
** AUTHOR:        Milan Kovacevic, December 10, 2001
**
** DESCRIPTION:    Interpolates linearly linear segments and with cubic spline
**                a nonlinear segement
**
** RETURNS:        Array of interpolated values
**
** AMENDMENTS:
**-----------------------------------------------------------------*/
{
    int    i;
    int    iMin,iMax;
    double    delta,deltaPrev,deltaNext;
//    double tol=1e-6;

    if (nOld < 3)
    {
        // Do linear interpolation here
        for (i=0; i<n;i++)
            FDInterpolateLinear(n,xOld,yOld,x[i]);
        return SUCCESS;
    }
    iMin = nOld-1;
    iMax = 0;
    // Determine the beginning of the nonlinear segment
    deltaPrev = (yOld[1]-yOld[0])/(xOld[1]-xOld[0]);
    delta = (yOld[2]-yOld[1])/(xOld[2]-xOld[1]);
    if (!DBL_EQUAL(deltaPrev,delta))
//    if (fabs(deltaPrev/delta - 1) > tol)
        iMin = 0;
    else
    {
        for (i=2; i<nOld-1; i++)
        {
            deltaPrev = delta;
            delta = (yOld[i+1]-yOld[i])/(xOld[i+1]-xOld[i]);
            if (!DBL_EQUAL(delta,deltaPrev) && i < nOld-2)
//            if (fabs(deltaPrev/delta - 1) > tol && i < nOld-2)
            {
                // Check whether a non-linear segment starts
                deltaNext = (yOld[i+2]-yOld[i+1])/(xOld[i+2]-xOld[i+1]);
                if (!DBL_EQUAL(deltaNext,delta))
//                if (fabs(deltaNext/delta - 1) > tol)
                {
                    // The non linear segment starts at i
                    iMin = i;
                    break;
                }                
            }
        }
    }
    // Find now the end of the nonlinear segment
    deltaNext = (yOld[nOld-1]-yOld[nOld-2])/(xOld[nOld-1]-xOld[nOld-2]);
    delta = (yOld[nOld-2]-yOld[nOld-3])/(xOld[nOld-2]-xOld[nOld-3]);
    if (!DBL_EQUAL(deltaNext,delta))
//    if (fabs(deltaNext/delta - 1) > tol)
        iMax = nOld-1;
    else
    {
        for (i=nOld-4; i>iMin; i--)
        {
            deltaNext = delta;
            delta = (yOld[i+1]-yOld[i])/(xOld[i+1]-xOld[i]);

            if (!DBL_EQUAL(delta,deltaNext) && i > 1)
            {
                // Check whether a non-linear segment ends
                deltaPrev = (yOld[i]-yOld[i-1])/(xOld[i]-xOld[i-1]);
                if (!DBL_EQUAL(deltaPrev,delta))
                {
                    // The non linear segment ends at i
                    iMax = i+1;
                    break;
                }                
            }
        }
    }

    if (iMax > iMin + 1)
    {
        // Do cubic spline interpolation for values between iMin and iMax
        // We interpolate for all x's since they may not be sorted. The unneeded ones will be overwritten later
        FDCubicSplineInterp(iMax-iMin+1,&(xOld[iMin]),&(yOld[iMin]),n,x,y);
    }

    // Do linear interpolation for the remaining points
    for (i=0; i<n;i++)
    {
        if (x[i] < xOld[iMin] && iMin > 0)
            y[i] = FDInterpolateLinear(iMin,xOld,yOld,x[i]);
        else if (x[i] > xOld[iMax] && iMax < nOld-1)
            y[i] = FDInterpolateLinear(nOld-iMax-1,&(xOld[iMax]),&(yOld[iMax]),x[i]);
    }
    return SUCCESS;
}

// Interpolates linearly linear segments and with cubic spline a nonlinear segement
int FDInterpolation1F(
    int        nOld,        // (I) size of the table to interpolate from
    double*    xOld,        // (I) x values of the table to interpolate from
    double*    yOld,        // (I) y values of the table to interpolate from
    int        n,            // (I) size of the output array
    double*    x,            // (I) array of x values at which to interpolate
    double*    y)            // (O) array of interpolated values for xvec
{
    int    i;
    int    iMin,iMax;
    double    delta,deltaPrev,deltaNext;
#ifdef FDUTILS_NEW
    double tol=1e-6;
#endif

    if (nOld < 3) {
        // Do linear interpolation here
        for (i=0; i<n;i++)
            y[i] = FDInterpolateLinear(n,xOld,yOld,x[i]);
        return SUCCESS;
    }
    iMin = nOld-1;
    iMax = 0;
    // Determine the beginning of the nonlinear segment
    deltaPrev = (yOld[1]-yOld[0])/(xOld[1]-xOld[0]);
    delta = (yOld[2]-yOld[1])/(xOld[2]-xOld[1]);
    #ifndef FDUTILS_NEW
       if (!Maths::equals(deltaPrev,delta))
    #else
       if (fabs(deltaPrev/delta - 1) > tol)
   #endif
        iMin = 0;
    else
    {
        for (i=2; i<nOld-1; i++)
        {
            deltaPrev = delta;
            delta = (yOld[i+1]-yOld[i])/(xOld[i+1]-xOld[i]);
#ifndef FDUTILS_NEW
            if (!Maths::equals(delta,deltaPrev) && i < nOld-2)
#else
            if (fabs(deltaPrev/delta - 1) > tol && i < nOld-2)
#endif
            {
                // Check whether a non-linear segment starts
                deltaNext = (yOld[i+2]-yOld[i+1])/(xOld[i+2]-xOld[i+1]);
#ifndef FDUTILS_NEW
                if (!Maths::equals(deltaNext,delta))
#else
                if (fabs(deltaNext/delta - 1) > tol)
#endif
                {
                    // The non linear segment starts at i
                    iMin = i;
                    break;
                }                
            }
        }
    }
    // Find now the end of the nonlinear segment
    deltaNext = (yOld[nOld-1]-yOld[nOld-2])/(xOld[nOld-1]-xOld[nOld-2]);
    delta = (yOld[nOld-2]-yOld[nOld-3])/(xOld[nOld-2]-xOld[nOld-3]);
#ifndef FDUTILS_NEW
    if (!Maths::equals(deltaNext,delta))
#else
    if (fabs(deltaNext/delta - 1) > tol)
#endif
        iMax = nOld-1;
    else
    {
        for (i=nOld-4; i>iMin; i--)
        {
            deltaNext = delta;
            delta = (yOld[i+1]-yOld[i])/(xOld[i+1]-xOld[i]);

#ifndef FDUTILS_NEW
            if (!Maths::equals(delta,deltaNext) && i > 1)
#else
            if (fabs(deltaNext/delta - 1) > tol && i > 1)
#endif
            {
                // Check whether a non-linear segment ends
                deltaPrev = (yOld[i]-yOld[i-1])/(xOld[i]-xOld[i-1]);
#ifndef FDUTILS_NEW
                if (!Maths::equals(deltaPrev,delta))
#else
                if (fabs(deltaPrev/delta - 1) > tol)
#endif
                {
                    // The non linear segment ends at i
                    iMax = i+1;
                    break;
                }                
            }
        }
    }

    if (iMax > iMin + 1)
    {
        // Do cubic spline interpolation for values between iMin and iMax
        // We interpolate for all x's since they may not be sorted. The unneeded ones will be overwritten later
        FDCubicSplineInterp(iMax-iMin+1,&(xOld[iMin]),&(yOld[iMin]),n,x,y);
    }

    // Do linear interpolation for the remaining points
    for (i=0; i<n;i++)
    {
        if (x[i] < xOld[iMin] && iMin > 0)
            y[i] = FDInterpolateLinear(iMin,xOld,yOld,x[i]);
        else if (x[i] > xOld[iMax] && iMax < nOld-1)
            y[i] = FDInterpolateLinear(nOld-iMax-1,&(xOld[iMax]),&(yOld[iMax]),x[i]);
    }
    return SUCCESS;
}

int FDInterpolationD(
    int        nOld,        /* (I) size of the table to interpolate from*/
    double*    xOld,        /* (I) x values of the table to interpolate from */
    double*    yOld,        /* (I) y values of the table to interpolate from */
    int        n,            /* (I) size of the output array */
    double*    x,        /* (I) array of x values at which to interpolate */
    double*    y,        /* (O) array of interpolated values for xvec*/
    double*    y1,        /* (O) array of first derivatives */
    double*    y2)        /* (O) array of second derivatives */
/*-------------------------------------------------------------------
** FUNCTION:    FDInterpolation
**
** AUTHOR:        Milan Kovacevic, December 10, 2001
**
** DESCRIPTION:    Interpolates linearly linear segments and with cubic spline
**                a nonlinear segement
**
** RETURNS:        Array of interpolated values
**
** AMENDMENTS:
**-----------------------------------------------------------------*/
{
    int    i;
    int    iMin,iMax;
    double    delta,deltaPrev,deltaNext;
//    double tol=1e-6;

    if (nOld < 3)
    {
        // Do linear interpolation here
        for (i=0; i<n;i++)
            FDInterpolateLinear(n,xOld,yOld,x[i]);
        return SUCCESS;
    }
    iMin = nOld-1;
    iMax = 0;
    // Determine the beginning of the nonlinear segment
    deltaPrev = (yOld[1]-yOld[0])/(xOld[1]-xOld[0]);
    delta = (yOld[2]-yOld[1])/(xOld[2]-xOld[1]);
    if (!DBL_EQUAL(deltaPrev,delta))
//    if (fabs(deltaPrev/delta - 1) > tol)
        iMin = 0;
    else
    {
        for (i=2; i<nOld-1; i++)
        {
            deltaPrev = delta;
            delta = (yOld[i+1]-yOld[i])/(xOld[i+1]-xOld[i]);
            if (!DBL_EQUAL(delta,deltaPrev) && i < nOld-2)
            {
                // Check whether a non-linear segment starts
                deltaNext = (yOld[i+2]-yOld[i+1])/(xOld[i+2]-xOld[i+1]);
                if (!DBL_EQUAL(deltaNext,delta))
//                if (fabs(deltaNext/delta - 1) > tol)
                {
                    // The non linear segment starts at i
                    iMin = i;
                    break;
                }                
            }
        }
    }
    // Find now the end of the nonlinear segment
    deltaNext = (yOld[nOld-1]-yOld[nOld-2])/(xOld[nOld-1]-xOld[nOld-2]);
    delta = (yOld[nOld-2]-yOld[nOld-3])/(xOld[nOld-2]-xOld[nOld-3]);
    if (!DBL_EQUAL(deltaNext,delta))
//    if (fabs(deltaNext/delta - 1) > tol)
        iMax = nOld-1;
    else
    {
        for (i=nOld-4; i>iMin; i--)
        {
            deltaNext = delta;
            delta = (yOld[i+1]-yOld[i])/(xOld[i+1]-xOld[i]);

            if (!DBL_EQUAL(delta,deltaNext) && i > 1)
//            if (fabs(deltaNext/delta - 1) > tol && i > 1)
            {
                // Check whether a non-linear segment ends
                deltaPrev = (yOld[i]-yOld[i-1])/(xOld[i]-xOld[i-1]);
                if (!DBL_EQUAL(deltaPrev,delta))
                {
                    // The non linear segment ends at i
                    iMax = i+1;
                    break;
                }                
            }
        }
    }

    if (iMax > iMin + 1)
    {
        // Do cubic spline interpolation for values between iMin and iMax
        // We interpolate for all x's since they may not be sorted. The unneeded ones will be overwritten later
        FDCubicSplineInterpD(iMax-iMin+1,&(xOld[iMin]),&(yOld[iMin]),n,x,y,y1,y2);
    }

    // Do linear interpolation for the remaining points
    for (i=0; i<n;i++)
    {
        if (x[i] < xOld[iMin] && iMin > 0)
        {
            y[i] = FDInterpolateLinearD(iMin,xOld,yOld,x[i],&(y1[i]));
            y2[i] = 0;
        }
        else if (x[i] > xOld[iMax] && iMax < nOld-1)
        {
            y[i] = FDInterpolateLinearD(nOld-iMax-1,&(xOld[iMax]),&(yOld[iMax]),x[i],&(y1[i]));
            y2[i] = 0;
        }
    }
    return SUCCESS;
}

// Interpolates linearly linear segments and with cubic spline
// a nonlinear segement
int FDInterpolationD1F(
    int        nOld,        // (I) size of the table to interpolate from
    double*    xOld,        // (I) x values of the table to interpolate from
    double*    yOld,        // (I) y values of the table to interpolate from
    int        n,            // (I) size of the output array
    double*    x,            // (I) array of x values at which to interpolate
    double*    y,            // (O) array of interpolated values for xvec
    double*    y1,            // (O) array of first derivatives
    double*    y2)            // (O) array of second derivatives
{
    int    i;
    int    iMin,iMax;
    double    delta,deltaPrev,deltaNext;
#ifdef FDUTILS_NEW
    double tol=1e-6;
#endif

    if (nOld < 3)
    {
        // Do linear interpolation here
        for (i=0; i<n;i++)
            y[i] = FDInterpolateLinear(n,xOld,yOld,x[i]);
        return SUCCESS;
    }
    iMin = nOld-1;
    iMax = 0;
    // Determine the beginning of the nonlinear segment
    deltaPrev = (yOld[1]-yOld[0])/(xOld[1]-xOld[0]);
    delta = (yOld[2]-yOld[1])/(xOld[2]-xOld[1]);
#ifndef FDUTILS_NEW
    if (!Maths::equals(deltaPrev,delta))
#else
    if (fabs(deltaPrev/delta - 1) > tol)
#endif
        iMin = 0;
    else
    {
        for (i=2; i<nOld-1; i++)
        {
            deltaPrev = delta;
            delta = (yOld[i+1]-yOld[i])/(xOld[i+1]-xOld[i]);
            if (!Maths::equals(delta,deltaPrev) && i < nOld-2)
            {
                // Check whether a non-linear segment starts
                deltaNext = (yOld[i+2]-yOld[i+1])/(xOld[i+2]-xOld[i+1]);
#ifndef FDUTILS_NEW
                if (!Maths::equals(deltaNext,delta))
#else
                if (fabs(deltaNext/delta - 1) > tol)
#endif
                {
                    // The non linear segment starts at i
                    iMin = i;
                    break;
                }                
            }
        }
    }
    // Find now the end of the nonlinear segment
    deltaNext = (yOld[nOld-1]-yOld[nOld-2])/(xOld[nOld-1]-xOld[nOld-2]);
    delta = (yOld[nOld-2]-yOld[nOld-3])/(xOld[nOld-2]-xOld[nOld-3]);
#ifndef FDUTILS_NEW
    if (!Maths::equals(deltaNext,delta))
#else
    if (fabs(deltaNext/delta - 1) > tol)
#endif
        iMax = nOld-1;
    else
    {
        for (i=nOld-4; i>iMin; i--)
        {
            deltaNext = delta;
            delta = (yOld[i+1]-yOld[i])/(xOld[i+1]-xOld[i]);

#ifndef FDUTILS_NEW
            if (!Maths::equals(delta,deltaNext) && i > 1)
#else
            if (fabs(deltaNext/delta - 1) > tol && i > 1)
#endif
            {
                // Check whether a non-linear segment ends
                deltaPrev = (yOld[i]-yOld[i-1])/(xOld[i]-xOld[i-1]);
#ifndef FDUTILS_NEW
                if (!Maths::equals(deltaPrev,delta))
#else
                if (fabs(deltaPrev/delta - 1) > tol && i > 1)
#endif
                {
                    // The non linear segment ends at i
                    iMax = i+1;
                    break;
                }                
            }
        }
    }

    if (iMax > iMin + 1)
    {
        // Do cubic spline interpolation for values between iMin and iMax
        // We interpolate for all x's since they may not be sorted. The unneeded ones will be overwritten later
        FDCubicSplineInterpD(iMax-iMin+1,&(xOld[iMin]),&(yOld[iMin]),n,x,y,y1,y2);
    }

    // Do linear interpolation for the remaining points
    for (i=0; i<n;i++)
    {
        if (x[i] < xOld[iMin] && iMin > 0)
        {
            y[i] = FDInterpolateLinearD(iMin,xOld,yOld,x[i],&(y1[i]));
            y2[i] = 0;
        }
        else if (x[i] > xOld[iMax] && iMax < nOld-1)
        {
            y[i] = FDInterpolateLinearD(nOld-iMax-1,&(xOld[iMax]),&(yOld[iMax]),x[i],&(y1[i]));
            y2[i] = 0;
        }
    }
    return SUCCESS;
}

int FDCubicSplineInterp(
    int        tabsize,    /* (I) size of the table to interpolate from*/
    double*    xtab,        /* (I) x values of the table to interpolate from */
    double*    ytab,        /* (I) y values of the table to interpolate from */
    int        vecsize,    /* (I) size of the output array */
    double*    xvec,        /* (I) array of x values at which to interpolate */
    double*    yvec)        /* (O) array of interpolated values for xvec*/
/*-------------------------------------------------------------------
** FUNCTION:    FDCubicSplineInterp
**
** AUTHOR:        Milan Kovacevic, February 24, 1998
**
** DESCRIPTION:    Cubic spline interpolator (derived from
**                CubicSplineLoop11)
**
** RETURNS:        Array of interpolated values
**
** AMENDMENTS:
**-----------------------------------------------------------------*/
{
    int        i;            /* index variable for x1 */
    int        loc12;        /* raw interpolation index */
    int        inew12;        /* adjusted interpolation index */
    int        imid2;        /* average of ilow and ihigh */
    int        ilow2;        /* index below interpolation index */
    int        ihigh2;        /* index above the interpolation index */
    double    x1;            /* value to interpolate at */
    double    tt12;        /* fraction of interval */
    double*    fpp;        /* second derivative array */
    double    ans1;        /* temporary */
    double fbuf[FDMAXSIZE];

    if (tabsize == 0)
    {
        for (i=0; i<vecsize; i++)
            yvec[i] = 0;
        return FAILURE;
    }
    else if (tabsize == 1)
    {
        for (i=0; i<vecsize; i++)
            yvec[i] = ytab[0];
        return SUCCESS;
    }

    if (tabsize <= FDMAXSIZE)
        fpp = fbuf;
    else
    {
        fpp = (double *)calloc(tabsize, sizeof(double));
        if (fpp == NULL)
        {
            return FAILURE;
        }
    }

    /* Find the cubic spline interpolants for a list of values. */
    if (FDGetSplineCoefficients(tabsize,xtab,ytab,fpp)==FAILURE)
        return FAILURE;
    /* Loop over the data values. */
    for (i=0; i<vecsize; i++)
    {
        x1 = xvec[i];
        /* Find the index of the data point. */
        ilow2 = -1;
        ihigh2 = tabsize;
        while (ihigh2 - ilow2>1)
        {
            imid2 = (int)(floor((ihigh2 + ilow2) / 2.0));
            if (x1>xtab[imid2]&&xtab[tabsize-1]<=xtab[0]||
                xtab[tabsize-1]>xtab[0]&&x1<=xtab[imid2])
            {
                ihigh2 = imid2;
            }
            else
            {
                ilow2 = imid2;
            }
        }
        loc12 = ilow2;
        /* Adjust the indices of the end points. */
        if (loc12<=-1)
        {
            inew12 = 0;
        }
        else
        {
            if (loc12>=tabsize-1)
            {
                inew12 = tabsize - 2;
            }
            else
            {
                inew12 = loc12;
            }
        }
        /* Compute the fraction of the interval. */
        tt12 = (x1 - xtab[inew12]) / (xtab[inew12 + 1] - xtab[inew12]);
        /* The linear interpolation. */
        ans1 = (1 - tt12)*ytab[inew12] + tt12*ytab[inew12 + 1];
        /* Cubic correction on the interior of the data. */
        if (tt12>=0.&&tt12<=1.)
            {
            ans1 = ans1 + (((::pow(1 - tt12,2) - 1)*(1 - tt12)*fpp[inew12]
                + tt12*(::pow(tt12,2) - 1)*fpp[inew12 + 1])*::pow(xtab[inew12 + 1]
                - xtab[inew12],2)) / 6;
            }
        yvec[i] = ans1;
    }
    if (tabsize > FDMAXSIZE)
        ::free(fpp);
    return SUCCESS;
}


int FDCubicSplineInterpD(
    int tabsize,            /* (I) size of the table to interpolate from */
    double *xtab,            /* (I) x values of the table to interpolate from */
    double *ytab,            /* (I) y values of the table to interpolate from */
    int vecsize,            /* (I) size of the output array */
    double *xvec,            /* (I) array of x values at which to interpolate */
    double *yvec,            /* (O) array of interpolated values */
    double *deriv1,            /* (O) first derivative */
    double *deriv2)            /* (0) second derivative */
/*-------------------------------------------------------------------
** FUNCTION:    FDCubicSplineInterpD
**
** AUTHOR:        Milan Kovacevic, February 24, 1998
**
** DESCRIPTION:    Cubic spline interpolator (derived from
**                CubicSplineLoop11). It is the same function
**                as FDCubicSplineInterp except that it also
**                returns first and second derivative
**
** RETURNS:        Array of interpolated values
**
** AMENDMENTS:
**-----------------------------------------------------------------*/
{
    int i,loc12,inew12,imid2,ilow2,ihigh2;
    double x1,tt12,*fpp,ans1;
    double fbuf[FDMAXSIZE];
    
    double xx1, xx2, a, b, y1, y2;
    double y2ndderiv1, y2ndderiv2;
    /* Key to program variables: */
    /* ans1: temporary */
    /* errorBuffer: place to store error messages */
    /* fpp: second derivative array */
    /* i: index variable for xi */
    /* ihigh2: index above the interpolation index */
    /* ilow2: index below interpolation index */
    /* tabsize: number of grid cells for xi */
    /* imid2: average of ilow and ihigh */
    /* inew12: adjusted interpolation index */
    /* loc12: raw interpolation index */
    /* ytab: ytab */
    /* xtab: solution variable */
    /* tt12: fraction of interval */
    /* x1: value to interpolate at */
    /* xvec: updated var */
    /* yvec: adjusted for discrete events */
    
    if (tabsize == 0)
    {
        for (i=0; i<vecsize; i++)
            yvec[i] = 0;
        return FAILURE;
    }
    else if (tabsize == 1)
    {
        for (i=0; i<vecsize; i++)
            yvec[i] = ytab[0];
        return SUCCESS;
    }
    
    if (tabsize <= FDMAXSIZE)
        fpp = fbuf;
    else
    {
        fpp = (double *)calloc(tabsize, sizeof(double));
        if (fpp == NULL)
        {
            return FAILURE;
        }
    }
    /* Find the cubic spline interpolants for a list of values. */
    if (FDGetSplineCoefficients(tabsize,xtab,ytab,fpp)==FAILURE)
        return FAILURE;
    /* Loop over the data values. */
    for (i=0; i<vecsize; i++) {
        x1 = xvec[i];
        /* Find the index of the data point. */
        ilow2 = -1;
        ihigh2 = tabsize;
        while (ihigh2 - ilow2>1) {
            imid2 = (int)(floor((ihigh2 + ilow2) / 2.0));
            if (x1>xtab[imid2]&&xtab[tabsize-1]<=xtab[0]||xtab[tabsize-1]>xtab[0]&&x1<=xtab[imid2])
                {
                ihigh2 = imid2;
                }
            else
                {
                ilow2 = imid2;
                }
        }
        loc12 = ilow2;
        /* Adjust the indices of the end points. */
        if (loc12<=-1)
            {
            inew12 = 0;
            }
        else
            {
            if (loc12>=tabsize-1)
                {
                inew12 = tabsize - 2;
                }
            else
                {
                inew12 = loc12;
                }
            }
        /* Compute the fraction of the interval. */
        tt12 = (x1 - xtab[inew12]) / (xtab[inew12 + 1] - xtab[inew12]);

        a = 1 - tt12;
        b = tt12;
        xx1 = xtab[inew12];
        xx2 = xtab[inew12 + 1];
        y1 = ytab[inew12];
        y2 = ytab[inew12 + 1];
        y2ndderiv1 = fpp[inew12];
        y2ndderiv2 = fpp[inew12 + 1];

        /* The linear interpolation. */
        ans1 = (1 - tt12)*ytab[inew12] + tt12*ytab[inew12 + 1];
        /* Cubic correction on the interior of the data. */
        if (tt12>=0.&&tt12<=1.)
        {
            ans1 = ans1 + (((::pow(1 - tt12,2) - 1)*(1 - tt12)*fpp[inew12
                ] + tt12*(::pow(tt12,2) - 1)*fpp[inew12 + 1])*::pow(xtab[inew12
                + 1] - xtab[inew12],2)) / 6;
        }
        yvec[i] = ans1;
        deriv1[i] = (y2 - y1)/(xx2 - xx1) -
            (3*a*a-1)*(xx2-xx1)*y2ndderiv1/6. +
            (3*b*b-1)*(xx2-xx1)*y2ndderiv2/6.;
        deriv2[i] = a*y2ndderiv1 + b*y2ndderiv2;
    }
    if (tabsize > FDMAXSIZE)
        ::free(fpp);
    return SUCCESS;
}

// Interpolates option values for the new FD grid points 
int    GridInterpolate1F(
    int        mOld,
    double*    Sold,
    double*    Vold,
    int        m,
    double*    S,
    double*    V,                // (O) vector of option values
    double    upBarrier,
    double    valueAtUpBar,
    double    valueAtUpBarDelta,
    double    downBarrier,
    double    valueAtDownBar,
    double    valueAtDownBarDelta)
//    int        iMaxOld,
//    int*    iMax,    // (O) maximum spot for the grid
{
    int        numSaved=0;
    double    Vsaved=0;
    double    upB=0;
    int        iMaxOld=mOld;
    double    alpha;
    int        i;

    // If we have the same grid, we don't need to do anything. We assume here that
    // type of coodinate transformation does not change..
    if (mOld != m || Sold[m] != S[m] || Sold[0] != S[0]) 
    {
        if (upBarrier > 0)
        {
            upB = upBarrier*0.9999999999;
            while (Sold[iMaxOld] >= upB)
                iMaxOld--;

            if (iMaxOld < mOld)    
            {
                iMaxOld++;

                if (Sold[iMaxOld] > upBarrier*1.0000001)
                {
                    // Adjust Vold[iMaxOld] to reflect the fact that the barrier is between grid points
                    double upPay;
                    alpha = (Sold[iMaxOld]-Sold[iMaxOld-1])/(upBarrier-Sold[iMaxOld-1]);
                    upPay = valueAtUpBar + valueAtUpBarDelta*Sold[iMaxOld];
                    Vsaved= Vold[iMaxOld];
                    numSaved++;
                    Vold[iMaxOld] = alpha*upPay+(1-alpha)*Vold[iMaxOld-1];
//                    V[iMaxOld] = alpha*upPay+(1-alpha)*Vold[iMaxOld-1]; // Bug in the production version
                }
            }
            else    // S[mOld] < upB so the barier will be ignored
                upBarrier = -1;
        }

        // Interpolate option values for the new spots
        FDInterpolation1F(iMaxOld+1,Sold,Vold,m+1,S,V);

        // Now apply barrier constraints to the interpolated value

        // Correct values that are above the barrier
        if (upBarrier>0)
        {
            if (numSaved > 0)
                Vold[iMaxOld] = Vsaved;    // Restore the old value
            
    //        upB = upBarrier*0.9999999999;
            iMaxOld=m;
            while (S[iMaxOld] >= upB)
                iMaxOld--;
            if (iMaxOld<m)
            {
                iMaxOld++;
                for (i=iMaxOld; i<=m; i++)
                    V[i] = valueAtUpBar + valueAtUpBarDelta*S[i];
            }
        }
    }
    else
        memcpy(V, Vold, (size_t)((m+1) * sizeof(double)));

    return SUCCESS;
}






////////////////////////////////////////////////////////////////////////////////////
//Create a new class for common fts used in FD1F and FD2F solver
//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////
//The following codes are used for solve FD1F in the framework of FD2F ADI
//ie: in the form of matrix for the purpose to speed the program.
//avoid the operations from extract lines or cols from matrix,
//call solverFD1F
//then, copy res back to the matrix
/////////////////////////////////////////////////////////////////////////////////

//-----------------------------------------------------------------------
                
void FDUtils::euroOneStep(
    int bot, int top,
    double * alpha, 
    double * beta, 
    double * gamma,
    double * A, 
    double * B, 
    double * C,
    double* sol_n,     
    double* sol_nplus1,
    bool solveByLine ,
    int which_i,
    int dim2)
{
    //dim will be the original size of FD, not necessary the current FD size
    //if solveByLine = True, dim = dim of col of sol_n
    double * rhs = (double *)::alloca( ( top - bot + 1 ) * sizeof( double ) ) - bot;

    //check dim
//    if (top - bot + 1 > dim){
//          throw ModelException("FDUtils::euroOneStep", 
//                        "Mismatch dim!");
//    }

    if (top - bot  <= 2){
          throw ModelException("FDUtils::euroOneStep", 
                        "Only 3 points!");
    }

    computeRhs(bot, top, 
        alpha, beta, gamma, A, B, C, rhs, sol_n, sol_nplus1, solveByLine, which_i, dim2);

    if (solveByLine == true) {
        //alphaY, betaY, gammaY
        triDiag2D (bot, top, alpha,  beta,  gamma, 
             rhs, sol_n, solveByLine ,which_i, dim2);
    
    }else{
        //alphaX, betaX, gammaX
        triDiag2D (bot, top, alpha,  beta,  gamma, 
             rhs, sol_n, solveByLine ,which_i, dim2);
    }
}

//-----------------------------------------------------------------------

void FDUtils::computeRhs(
    int bot, int top,
    double * alpha,
    double * beta,
    double * gamma,
    double * A,
    double * B,
    double * C,
    double * rhs,
    double* sol_n,
    double* sol_nplus1,
    bool solveByLine,
    int which_i,
    int dim2 )
{
    //same issue, maybe need to check the consistence of dim
    int ix;
    int l = bot + 1;
    int t = top -1 ;
    int k;

    // maybe need to add check dim for dim of sol_nplus1

    if (solveByLine == true) {
        for (ix = l ; ix <= t; ix++){
            //aY, By, Cy
//            rhs[ix] = A[ix] * sol_nplus1[which_i][ix-1] 
//                    + B[ix] * sol_nplus1[which_i][ix] 
//                    + C[ix] * sol_nplus1[which_i][ix+1];

//            k = which_i + (ix) * (dim1 +1);  //index on dim 1 goes from 0 to dim1
//
//            rhs[ix] = A[ix] * sol_nplus1[k - (dim1+1)] 
//                    + B[ix] * sol_nplus1[k] 
//                    + C[ix] * sol_nplus1[k + (dim1+1)];


            k = ix + which_i * dim2;  //index on dim 2 goes from 0 to dim2

            rhs[ix] = A[ix] * sol_nplus1[k - 1]
                    + B[ix] * sol_nplus1[k]
                    + C[ix] * sol_nplus1[k + 1];

        
        }

        k = (l-1) + which_i * dim2;
        rhs[l] -= alpha[l] * sol_n[k]; 

        k = (t+1) + which_i * dim2;
        rhs[t] -= gamma[t] * sol_n[k] ;
    }else{
        //aX, bX, cX
        for (ix = l; ix <= t; ix++){
//            rhs[ix] = A[ix] * sol_nplus1[ix-1][which_i] 
//                    + B[ix] * sol_nplus1[ix][which_i] 
//                    + C[ix] * sol_nplus1[ix+1][which_i];

//            k = ix + which_i * (dim1 +1);

            k = which_i + ix * dim2;

            rhs[ix] = A[ix] * sol_nplus1[k - dim2]
                    + B[ix] * sol_nplus1[k]
                    + C[ix] * sol_nplus1[k + dim2];


        }
        k = which_i + (l-1) * dim2;

        rhs[l] -= alpha[l] * sol_n[k] ;

        k = which_i + (t+1) * dim2;
        rhs[t] -= gamma[t] * sol_n[k];

//        rhs[l] -= alpha[l] * sol_n[l-1][which_i] ;
//        rhs[t] -= gamma[t] * sol_n[t+1][which_i];

/*        
        for (ix = l; ix <= t; ix++){
            rhs[idx[ix]] = A[idx[ix]] * sol_nplus1[idx[ix]-1][which_i] 
                    + B[idx[ix]] * sol_nplus1[idx[ix]][which_i] 
                    + C[idx[ix]] * sol_nplus1[idx[ix]+1][which_i];

        }
        rhs[l] -= alpha[l] * sol_n[l-1][which_i] ;
        rhs[t] -= gamma[t] * sol_n[t+1][which_i];
*/

    }
}

//-----------------------------------------------------------------------------    

void FDUtils::euroOneStepWithSource(
    int bot, int top, 
    double * alpha, 
    double * beta, 
    double * gamma,
    double * A, 
    double * B, 
    double * C,                         
    double* sol_n, 
    double* sol_nplus1, 
    double* source,
    bool solveByLine,
    int which_i, int dim2)
{
    double * rhs = (double *)::alloca( ( top - bot + 1 ) * sizeof( double ) ) - bot;

    computeRhsWithSource(bot, top, 
                        alpha, beta, gamma, A, B, C,
                        rhs,  sol_n, sol_nplus1, 
                        source, solveByLine ,which_i, dim2);  

    if (solveByLine == true){
        //alphaY, betaY, gammaY
        triDiag2D (bot, top, alpha,  beta,  gamma, 
             rhs, sol_n, solveByLine ,which_i, dim2);

    }else{
        triDiag2D (bot, top, alpha,  beta,  gamma, 
             rhs, sol_n, solveByLine ,which_i, dim2);
    }
}

//-----------------------------------------------------------------------

void FDUtils::computeRhsWithSource(
    int bot, int top,
    double * alpha,
    double * beta,
    double * gamma,
    double * A,
    double * B,
    double * C,
    double * rhs,
    double* sol_n,
    double* sol_nplus1,
    double* source,
    bool solveByLine,
    int which_i,
    int dim2)
{
    int ix;
    int l = bot +1;
    int t = top -1;
    int k;

    computeRhs(bot, top, 
        alpha, beta, gamma, A, B, C,
        rhs, sol_n, sol_nplus1,solveByLine, which_i, dim2);
    if (solveByLine == true){
        for (ix = l; ix <= t ; ix++){
            k = ix + which_i * dim2;
//            rhs[ix] += source[which_i][ix];
            rhs[ix] += source[k];
        }
    }else{
        for (ix = l; ix <= t ; ix++){
            k = which_i + ix * dim2;
//            rhs[ix] += source[ix][which_i];            
            rhs[ix] += source[k];            
        }
    }
}

// triDiag2D is solving just for one dimension... right?
// copied from cmgb_invert_tridia.cpp comments and replace the varianble
/*
            |beta(1)    gamma(1)                                            | | sol(1)    |    | rhs(1)   |
            |alpha(2)    beta(2)        gamma(2)                                | | sol(2)    |    | rhs(1)   |
            |                                                                | |            |=    |          |
            |                                                                | |            |    |          |
            |                                    beta(n-1)    gamma(n-1)        | | sol(n-1)|    | rhs(n-1) |
            |                                    alpha(n)    beta(n)            | | sol(n)    |    | rhs(n)   |

                = M                                                                =X                =Y

                B(n) = beta(n)
                Y(n) = rhs(n)
                
                for i = n-1 --> 1
                    B(i) = beta(i) - Sup(i).Inf(i+1)/B(i+1)
                    Y(i) = rhs(i) - Sup(i).Y(i+1)/B(i+1)

                sol(1) = Y(1)/B(1)

                for i=2 --> n
                    sol(i) = (Y(i) - Inf(i).sol(i-1))/B(i)

*/

//-----------------------------------------------------------------------
//bot: the bot index of FD 
//top: the top index of FD
//Dim = top - bot + 1 is the nb of points
//Dim -2 will be the unknown v.a, (bot+1, top -1)
 
void FDUtils::triDiag2D(
    int bot, int top,                                  
    double * alpha, 
    double * beta,
    double * gamma,
    double * rhs, 
    double* sol,
    bool solveByLine,
    int which_i,
    int dim2)
{
    //to be safe, sould check the dim of vector and matrix, to do...
    //check dim error
    if (bot > top) {
        // hard coded for now !!!
        throw ModelException("FDUtils::triDiag2D ", 
        "solver does not support 9 coefficients.");
    }

    /******************************************************************/
    /*   alpha:  lowdiag;                                          
    //   beta:  diag;                                                   
    //   gamma:  updiag;                                                
    //   rhs:  member at right of eq;                                            
    //   Ndim: dimension of the matrix;                             
    //////////////////////////////////////////////////////////////////*/

    int ix;
    int t = top - 1;
    int l = bot +1;
    int k;
    //t = Ndim -1, l=0
    //to set the same size for easy dealing with index

    double * betatemp = (double *)::alloca( ( top - bot + 1 ) * sizeof( double ) ) - bot;
    double * rhstemp = (double *)::alloca( ( top - bot + 1 ) * sizeof( double ) ) - bot;

    betatemp[t] = beta[t];
    for (ix= t-1; ix >= l;  ix--)
        betatemp[ix] = beta[ix] - gamma[ix]*alpha[ix+1]/betatemp[ix+1] ;

    /*************/
    /* backward  */
    /*************/
    rhstemp[t] = rhs[t]; 
    for (ix = t-1; ix >= l;  ix--){
        rhstemp[ix] = rhs[ix] 
            - gamma[ix]*rhstemp[ix+1]/betatemp[ix+1] ;
    }
    
    /*************/
    /* forward  */
    /*************/

    if (solveByLine == true ) {
        k = l + which_i * dim2;
//        sol[which_i][l] = rhstemp[l]/betatemp[l]    ;
        sol[k] = rhstemp[l]/betatemp[l]    ;
        for (ix = l + 1; ix <= t ; ix++){
            k = ix + which_i * dim2;

            sol[k] = rhstemp[ix]/betatemp[ix] 
                - alpha[ix]*sol[k - 1]/betatemp[ix]    ;      

//            sol[which_i][ix] = rhstemp[ix]/betatemp[ix] 
//                - alpha[ix]*sol[which_i][ix-1]/betatemp[ix]    ;      
        }
    }
    else{//solveByCol
        k = which_i + l * dim2;
        sol[k] = rhstemp[l]/betatemp[l]    ;
        for (ix = l + 1; ix <= t ; ix++){    
            k = which_i + ix * dim2;

            sol[k] = rhstemp[ix]/betatemp[ix] 
                - alpha[ix]*sol[k - dim2]/betatemp[ix]    ;              
        }

//        sol[l][which_i] = rhstemp[l]/betatemp[l]    ;
//        for (ix = l + 1; ix <= t ; ix++){    
//            sol[ix][which_i] = rhstemp[ix]/betatemp[ix] 
//                - alpha[ix]*sol[ix-1][which_i]/betatemp[ix]    ;              
//        }
    }
}

//-------------------------------------------------------------
//some utility ft to set segments for barrier option
//add +- 1d around barrier date
//--------------------------------------------------------------

void FDUtils::calcSegBarDates(const DateTime& valueDate,
                                const DateTime& matDate,
                                TimeMetricConstSP metric,
                                const DateTimeArray& barDates,
                                DateTimeArray& segDates, 
                                IntArray* isAddedSeg){
    //double diffD;
    double oneD = 1.0 / 250.0;
    DateTime lastD; 
    DateTime nextD; 

    //segDate[0]=valuedate
    segDates.push_back(valueDate);

    for (int i=0; i < barDates.size(); i++)
    {
        if (barDates[i] > valueDate && barDates[i] <= matDate)
        {
            //bar dates as seg date, additionally, added +-1d as additional seg.
            double dummy;

            //left
            DateTime lastDin; 
            if(i==0){
                lastDin = valueDate;
            }else{
                lastDin = barDates[i-1];
            }

            DateTime date = metric->impliedTime(barDates[i], - oneD, dummy); 

            if ( date <= lastDin){
                if (i ==0){
                    (*isAddedSeg)[0]=0; //overwrite the default value (-1) at index 0 to 0
                }else{
                    isAddedSeg->push_back(0);                    
                }
            }else{//
                // add at -1 trading day
                lastD = segDates[segDates.size() -1];

                if(date <= lastD){
                    if (i==0){
                        (*isAddedSeg)[0] = 0; //overwrite the default value (-1) at index 0 to 0
                    }else{
                        isAddedSeg->push_back(0);
                    }
                }else{

                    if (i ==0){
                        (*isAddedSeg)[0]=1; //overwrite the default value (-1) at index 0 to 0
                    }else{
                        isAddedSeg->push_back(1);                    
                    }

                    //isAddedSeg->push_back(1);
                    segDates.push_back(date); 
                    isAddedSeg->push_back(0);
                }
            }
            
            //add bar date to critDate
            segDates.push_back(barDates[i]);

            //right
            // add at 1 trading day
            date = metric->impliedTime(barDates[i], oneD, dummy); 

            if(i == barDates.size() - 1){
                nextD =matDate;
            }else{
                nextD = barDates[i+1];
            }

            if (date < nextD){
                segDates.push_back(date); 
                isAddedSeg->push_back(0);                    
            }
        }
    }

    //for last seg
    if (segDates[segDates.size()-1] < matDate){
        isAddedSeg->push_back(1);    
        segDates.push_back(matDate);
    }
}

void FDUtils::SetSegBarDates1Bar(const DateTime& valueDate,
                                    const DateTime& matDate,
                                    TimeMetricConstSP metric,
                                    const ScheduleSP bar, 
                                    DateTimeArray& segDates, IntArray* isAddedSeg){


    if ( bar->getInterp() != Schedule::INTERP_NONE){
        segDates.push_back(valueDate);

        if (bar->firstDate() > valueDate){
            segDates.push_back(bar->firstDate());
        }
        if (bar->lastDate() < matDate){
            segDates.push_back(bar->lastDate());
        }

        segDates.push_back(matDate);

    }else{
        calcSegBarDates(valueDate, matDate, metric, bar->getDates(), segDates, isAddedSeg);
    }
}

void FDUtils::SetSegBarDates2Bar(const DateTime& valueDate,
                                       const DateTime& matDate,
                                       TimeMetricConstSP metric,
                                       const ScheduleSP ubar, const ScheduleSP lbar, 
                                       DateTimeArray& segDates, IntArray* isAddedSeg){
    //segDates //O
    //isAddedSeg //O

    if (( ubar->getInterp() != Schedule::INTERP_NONE) && ( lbar->getInterp() != Schedule::INTERP_NONE) ){

        // add monitor start and monitor end to segment dates list if needed

        //segDates.push_back(valueDate);
        segDates.push_back(ubar->firstDate());
        segDates.push_back(lbar->firstDate());
        segDates.push_back(ubar->lastDate());
        segDates.push_back(lbar->lastDate());

        //seg[0] = valueDate, and last will be matDates
        metric->SortDate(valueDate , matDate, true, segDates);
        
        //wants to insert valueDate to seg[0]
        segDates.insert(segDates.begin(),  valueDate);

        segDates.push_back(matDate);
    }
    else{

        if(( ubar->getInterp() == Schedule::INTERP_NONE) && ( lbar->getInterp() != Schedule::INTERP_NONE) ){
            //upBarrier discrete
            const DateTimeArray& barDates = ubar->getDates();

            calcSegBarDates(valueDate, matDate, metric, barDates, segDates, isAddedSeg);

        }else if(( ubar->getInterp() != Schedule::INTERP_NONE) && ( lbar->getInterp() == Schedule::INTERP_NONE) ){
            const DateTimeArray& barDates = lbar->getDates();

            calcSegBarDates(valueDate, matDate, metric, barDates, segDates, isAddedSeg);

        }else{
            const DateTimeArray& barDates = DateTime::merge(ubar->getDates(), 
                                        lbar->getDates());  
            calcSegBarDates(valueDate, matDate, metric, barDates, segDates, isAddedSeg);
        }
    }
}

void FDUtils::SetSegDates(const DateTime& valueDate,
                                        const DateTime& matDate,
                                        TimeMetricConstSP metric,
                                        const ScheduleSP ubar, const ScheduleSP lbar, 
                                        const string& uBarType, const string& lBarType,
                                        DateTimeArray& segDates, IntArray* isAddedSeg){
    //segDates //O
    //isAddedSeg //O

    //isAddedSeg[0]=-1 indicates that we didn't add any additional seg.
    //timeline will use this info
    //needed otherwise TimeLine willn't work when up and down barrier both = NA
    (*isAddedSeg).push_back(-1);

    if ( (uBarType != "NA" ) && (lBarType != "NA" )  ){
        //has up and low barriers

        SetSegBarDates2Bar(valueDate, matDate,
                            metric,
                            ubar, 
                            lbar, 
                            segDates, isAddedSeg);

    }else if( (uBarType != "NA" ) && (lBarType == "NA" )  ){
        //has upper barrier
        SetSegBarDates1Bar(valueDate, matDate,
                            metric,
                           ubar , 
                           segDates, isAddedSeg);

    }else if( (uBarType == "NA" ) && (lBarType != "NA" )  ){
        //has low barrier
        SetSegBarDates1Bar(valueDate, matDate,
                            metric,
                           lbar, 
                           segDates, isAddedSeg);

    }else{//no barrier
        segDates.push_back(valueDate);
        segDates.push_back(matDate);
    }
}

/** more than 2 barriers */
void FDUtils::SetSegDatesGen(const DateTime& valueDate,
                                const DateTime& matDate,
                                TimeMetricConstSP metric,
                                const vector<ScheduleSP> bar, 
                                const vector<string>& barType, 
                                DateTimeArray& segDates, IntArray* isAddedSeg){
    //segDates //O
    //isAddedSeg //O

    int size = barType.size();
    int i;
    bool needSpecialSegSet = false;
    bool hasBarrier = false;
    bool firstTime = false;
    DateTimeArray barDates ;

    //isAddedSeg[0]=-1 indicates that we didn't add any additional seg.
    //timeline will use this info
    //needed otherwise TimeLine willn't work when up and down barrier both = NA
    (*isAddedSeg).push_back(-1);

    for(i = 0; i < size; i++){
        if (barType[i] != "NA"){
            if(bar[i]->getInterp() == Schedule::INTERP_NONE){//continous barrier case
                if (firstTime == false){
                    barDates = bar[i]->getDates();
                    firstTime = true;
                }else{
                    barDates = DateTime::merge(barDates, bar[i]->getDates());   
                }
                needSpecialSegSet = true;
            }
            hasBarrier = true; //has at least 1 barrier
        }
    }
    
    if(needSpecialSegSet == true){//at least 1 discrete barrier
        calcSegBarDates(valueDate, matDate, metric, barDates, segDates, isAddedSeg);
    }else{
        if(hasBarrier){//all continous barriers
            for(i = 0; i < size; i++){
                if (barType[i] != "NA"){
                    if(bar[i]->getInterp() != Schedule::INTERP_NONE){//NA or continous barrier case
                        segDates.push_back(bar[i]->firstDate());
                        segDates.push_back(bar[i]->lastDate());
                    }
                }
            }
            //seg[0] = valueDate, and last will be matDates
            metric->SortDate(valueDate , matDate, true, segDates);            
            //wants to insert valueDate to seg[0]
            segDates.insert(segDates.begin(),  valueDate);
            segDates.push_back(matDate);

        }else{//no barriers
            segDates.push_back(valueDate);
            segDates.push_back(matDate);
        }
    }
}


//collecting barrier info and set up average critical points
//return critical points value
double FDUtils::setCritSpacePts(const DateTime& valueDate,
                                const DateTime& matDate,
                                const ScheduleSP in){
    const DateTimeArray& inDates = in->getDates();
    const DoubleArray& inValues = in->getValues();

    int j;
    double FP_MIN = 1.0e-10; //copy from tree1f, temporary

    double maxV = - 2.0*FP_MIN;
    double minV = 1000000.0 / FP_MIN;
    double avgV = 0.0;

    int count = 0;    
    for (j=0; j < in->length(); j++){
        if ( (inDates[j] >= valueDate) && 
            (inDates[j]<= matDate) ){

            minV = Maths::min(minV, inValues[j]);
            maxV = Maths::max(maxV, inValues[j]);
            avgV += inValues[j];
            count++;
        }
    }

    if (count > 2){
        //take out thebiggest and smallest
        avgV = (avgV - minV - maxV)/(count - 2);
    }else if (count != 0){
        avgV = avgV / count;
    }

    return avgV;
}

/** more than 2 barriers */
void FDUtils::setCritSpacePtsAllGen(const DateTime& valueDate,
                                    const DateTime& matDate,
                                    const  vector<ScheduleSP> bar,
                                    const  vector<string>& barType,
                                    DoubleArray& outCritPts){

    int size = barType.size();
    int i;
    int count = 0;

    outCritPts.resize(size);
    for (i = 0; i < size; i++){
        if(barType[i] != "NA"){
            outCritPts[count] = setCritSpacePts( valueDate, matDate, bar[i]); 
            count = count +1;
        }
    }
    outCritPts.resize(count);
}

void FDUtils::setCritSpacePtsAll(const DateTime& valueDate,
                                const DateTime& matDate,
                            const ScheduleSP ubar, const ScheduleSP lbar, 
                            const string& uBarType, const string& lBarType,
                            DoubleArray& outCritPts){

    if ( (uBarType != "NA" ) && (lBarType != "NA" )  ){
        //has up and low barriers
        outCritPts.resize(2);
        outCritPts[0] = setCritSpacePts( valueDate, matDate, lbar);
        outCritPts[1] = setCritSpacePts( valueDate, matDate, ubar);

    }else if( (uBarType != "NA" ) && (lBarType == "NA" )  ){
        //has upper barrier
        outCritPts.resize(1);
        outCritPts[0] = setCritSpacePts( valueDate, matDate, ubar);


    }else if( (uBarType == "NA" ) && (lBarType != "NA" )  ){
        //has low barrier
        outCritPts.resize(1);
        outCritPts[0] = setCritSpacePts( valueDate, matDate, lbar);
    }else{//no barrier
        outCritPts.resize(0);
    }
}

//----------------------------------------------------------------------------

void FDUtils::addMoreSpacePts(vector<double>& v_dxM, 
                              DoubleArray& critSpacePts, int dim1, int minNe,
                                double w5, double w3, double w2,
                                double& lBound5, double& uBound5, 
                                double& lBound3, double& uBound3, 
                                double& lBound2, double& uBound2, 
                                double oneDVolT){

    int nbCritSpts = critSpacePts.size(); 

    double lB = lBound2;
    double uB = uBound2;
    double dxTemp = (uB - lB)/(dim1);
    int impNe = int ((2 * oneDVolT) / dxTemp );


    if (nbCritSpts ==2){
        //total, allocate min minNe pts, here has 2 pts
        impNe = Maths::max(impNe, minNe / 2);
    }else if (nbCritSpts ==1){
        impNe = Maths::max(impNe, minNe);
    }else{
        impNe = 0;
    }

    if (nbCritSpts ==2){
        double barLow = critSpacePts[0];
        double barUp = critSpacePts[1];

        double barLowL = log(barLow) - oneDVolT;
        double barLowR = log(barLow) + oneDVolT;

        double barUpL = log(barUp) - oneDVolT;
        double barUpR = log(barUp) + oneDVolT;

        addMoreSpacePtsAt2Pts(v_dxM, dim1, impNe, w5, w3, w2,
                        lBound5, uBound5, 
                        lBound3, uBound3, 
                        lBound2, uBound2, 
                        barLowL, barLowR,
                        barUpL, barUpR);
    }else if(nbCritSpts ==1){
        double bar = critSpacePts[0];
        double barL = log(bar) - oneDVolT;
        double barR = log(bar) + oneDVolT;

        addMoreSpacePtsAt1Pt(v_dxM, dim1, impNe, w5, w3, w2,
                        lBound5, uBound5, 
                        lBound3, uBound3, 
                        lBound2, uBound2, 
                        barL, barR);
    }
}

void FDUtils::addMoreSpacePtsAt2Pts(vector<double>& v_dxM, int dim1, 
                                int impNe, double w5, double w3, double w2,
                                double& lBound5, double& uBound5, 
                                double& lBound3, double& uBound3, 
                                double& lBound2, double& uBound2, 
                                double barLowL, double barLowR,
                                double barUpL, double barUpR){
    int l, h;

    if (barLowR > barUpL){
        throw ModelException("FDUtils::addMoreSpacePts", "not implemented yet! --up and low barriers are too close!");
    }else{
        
        int dim1R = Maths::max(dim1 - 2 * impNe, 0);
        int n5 = int (w5 * dim1R);
        int n3 = int (w3 * dim1R);
        int n2 = dim1R - 2 * n5 - 2 * n3 ;

        if ((barLowL >= lBound2) && (barUpL <= uBound2)){//barLow, barUp both within [lbound2, uBound2]
            uBound2 = Maths::max(uBound2, barUpR);
            
            //lBound5, lBound3, n5
            l = 0;
            h =  l + n5;
            calcLocaldx(v_dxM, lBound5, lBound3, l , h);

            //lbound3, lBound2, n3
            l = n5;
            h =  l + n3;
            calcLocaldx(v_dxM, lBound3, lBound2, l , h);

            //lBound2, uBound2, n2 + 2* impNe
            l = n3 + n5;
            h = n3 + n5 + n2 + 2 * impNe;
            addMoreSpacePtsConst(v_dxM, impNe, lBound2, uBound2, l, h, 
                                barLowL, barLowR, barUpL, barUpR);

            //uBound2, uBound3, n3,
            l = n5 + n3 + n2 + 2 * impNe;
            h =  l + n3;
            calcLocaldx(v_dxM, uBound2, uBound3, l , h);

            //uBound3, uBound5, n5
            l = n5 + n3 + n2 + 2 * impNe + n3;
            h =  l + n5;
            calcLocaldx(v_dxM, uBound3, uBound5, l , h);

        }else if((barLowL >= lBound3) && (barUpR <= uBound3)){    //barLow, barUp both within [lbound3, uBound3]
            uBound3 = Maths::max(uBound3, barUpR);

            //lBound5, lBound3, n5
            l = 0;
            h =  l + n5;
            calcLocaldx(v_dxM, lBound5, lBound3, l , h);

            //lBound3, uBound3, 2*(n3+impNe) + n2
            l = n5;
            h = n5 + 2 * (n3 + impNe) + n2;
            addMoreSpacePtsConst(v_dxM, impNe, lBound3, uBound3, l, h, 
                                barLowL, barLowR, barUpL, barUpR);

            //uBound3, uBound5, n5
            l = n5 + 2 * (n3 + impNe) + n2;
            h =  l + n5;
            calcLocaldx(v_dxM, uBound3, uBound5, l , h);

        }else if((barLowL >= lBound5) && (barLowL < lBound3) ){
            //barLow within lBound5 lBound3, set grid based on barUp, then ajust to barLow

            lBound3 = Maths::max(lBound3, barLowR);

            //set grid based on barUp
            //since we don't add impNe pts around barLow
            //we need to add impNe pts back to barUp

            addMoreSpacePtsAt1Pt(v_dxM, dim1, 2*impNe, w5, w3, w2,
                        lBound5, uBound5, 
                        lBound3, uBound3, 
                        lBound2, uBound2, barUpL, barUpR);

            //ajust based on barLow, pb: n5 may < impNe???,
            //to do
            
        }else if ((barUpL > uBound3) && (barUpR <= uBound5) ){
            //barUp within uBound3, uBound5, set grid based on barLow, then adjust to barUp

            //since we don't add impNe pts around barUp
            //we need to add impNe pts back to barLow

            addMoreSpacePtsAt1Pt(v_dxM, dim1, 2*impNe, w5, w3, w2,
                        lBound5, uBound5, 
                        lBound3, uBound3, 
                        lBound2, uBound2, barLowL, barLowR);

        }else{

            //barLow in (lBound5, lBound3), and barUp in (uBound3, uBound5)
            //add points at thses 2 crit pts, and const grid otherwise
            addMoreSpacePtsConst(v_dxM, impNe, lBound5, uBound5, 0, dim1, 
                                barLowL, barLowR, barUpL, barUpR);
        }
    }
}

void FDUtils::addMoreSpacePtsAt1Pt(vector<double>& v_dxM, int dim1, int impNe, double w5, double w3, double w2,
                                double& lBound5, double& uBound5, 
                                double& lBound3, double& uBound3, 
                                double& lBound2, double& uBound2, 
                                double barL, double barR ){

    int n5, n3, n2;
    double wL;
    int l, h;

    double oneDVolT = (barR - barL) / 2.0;

    int dim1R = Maths::max(dim1 - impNe, 0);
    n5 = int(w5 * dim1R);
    n3 = int(w3 * dim1R);
    n2 = int(w2 * dim1R);

    if ((barL >= lBound2) && (barR <= uBound2)){

        wL = (barL - lBound2)/(uBound2 - lBound2 - 2*oneDVolT);
        int n2L = int(wL * n2);

        //lbound5, lbound3, n5
        l = 0;
        h =  l + n5;
        calcLocaldx(v_dxM, lBound5, lBound3, l , h);

        //lbound3, lbound2, n3
        l = n5;
        h =  l + n3;
        calcLocaldx(v_dxM, lBound3, lBound2, l , h);

        //lbound2, barL, n2*weight
        l = n5 + n3;
        h =  l + n2L;
        calcLocaldx(v_dxM, lBound2, barL, l , h);

        //barL, barR, impNe
        l = n5 + n3 + n2L;
        h =  l + impNe;
        calcLocaldx(v_dxM, barL, barR, l , h);

        //barR, uBound2, n2*(1-weight)
        l = n5 + n3 + n2L + impNe;
        h =  l + (n2 - n2L);
        calcLocaldx(v_dxM, barR, uBound2, l , h);

        //uBound2, uBound3, n3
        l = n5 + n3 + n2L + impNe + (n2 - n2L);
        h =  l + n3;
        calcLocaldx(v_dxM, uBound2, uBound3, l , h);

        //uBound3, uBound5, n5
        l = n5 + n3 + n2L + impNe + (n2 - n2L) + n3;
        h =  l + n5;
        calcLocaldx(v_dxM, uBound3, uBound5, l , h);

     }else if ((barL >= lBound3) && (barL < lBound2) ){
         //ajust lbound2 = min(lBound2, barR)

        lBound2 = Maths::max(lBound2, barR);

        wL = (barL - lBound3)/(lBound2 - lBound3 - 2*oneDVolT);
        int n3L = int(wL * n3);

     //lbound5, lbound3, n5
        l = 0;
        h =  l + n5;
        calcLocaldx(v_dxM, lBound5, lBound3, l , h);


     //lbound3, barL, n3 * weight
        l = n5;
        h =  l + n3L;
        calcLocaldx(v_dxM, lBound3, barL, l , h);

     //barL, barR, impNe
        l = n5 + n3L;
        h =  l + impNe;
        calcLocaldx(v_dxM, barL, barR, l , h);

     //barR, lBound2, n3*(1-weight)
        l = n5 + n3L + impNe;
        h =  l + (n3 - n3L);
        calcLocaldx(v_dxM, barR, lBound2, l , h);

     //lBound2, uBound2, n2
        l = n5 + n3L + impNe + (n3 - n3L);
        h =  l + n2;
        calcLocaldx(v_dxM, lBound2, uBound2, l , h);

     //uBound2, uBound3, n3
        l = n5 + n3L + impNe + (n3 - n3L) + n2;
        h =  l + n3;
        calcLocaldx(v_dxM, uBound2, uBound3, l , h);

     //uBound3, uBound5, n5
        l = n5 + n3L + impNe + (n3 - n3L) + n2 + n3;
        h =  l + n5;
        calcLocaldx(v_dxM, uBound3, uBound5, l , h);

     }else if((barR > uBound2) && (barR <= uBound3)){
         //adjust uBound2 = max(uBound2, barL)
        uBound2 = Maths::min(uBound2, barL);

        wL = (barL - uBound2)/(uBound3 - uBound2 - 2*oneDVolT);
        int n3L = int(wL * n3);

         //lbound5, lbound3, n5
        l = 0;
        h =  l + n5;
        calcLocaldx(v_dxM, lBound5, lBound3, l , h);

         //lbound3, lbound2, n3
        l = n5;
        h =  l + n3;
        calcLocaldx(v_dxM, lBound3, lBound2, l , h);

         //lbound2, uBound2, n2
        l = n5 + n3;
        h =  l + n2;
        calcLocaldx(v_dxM, lBound2, uBound2, l , h);

         //uBound2, barL, n3 * w
        l = n5 + n3 + n2;
        h =  l + n3L;
        calcLocaldx(v_dxM, uBound2, barL, l , h);

         //barL, barR, impNe
        l = n5 + n3 + n2 + n3L;
        h =  l + impNe;
        calcLocaldx(v_dxM, barL, barR, l , h);

         //barR, uBound3, n2*(1-weight)
        l = n5 + n3 + n2 + n3L +impNe;
        h =  l + (n3 - n3L);
        calcLocaldx(v_dxM, barR, uBound3, l , h);

         //uBound3, uBound5, n5
        l = n5 + n3 + n2 + n3L +impNe + (n3-n3L);
        h =  l + n5;
        calcLocaldx(v_dxM, uBound3, uBound5, l , h);

     }else{ // (lbound5, lBound3) or (uBound3, uBound5) 
         // add pts to (barL, barR), const for remaining

        wL = (barL - lBound5) / (uBound5 - lBound5 -2 * oneDVolT);
        n5 = int (wL * dim1R);
        n2 = dim1R - n5 ;

         //lBound5, barL, n5
        l = 0;
        h =  l + n5;
        calcLocaldx(v_dxM, lBound5, barL, l , h);

         //barL, barR, impNe
        l = n5;
        h =  l + impNe;
        calcLocaldx(v_dxM, barL, barR, l , h);

         //barR, uBound5, n2
        l = n5 + impNe;
        h =  l + n2;
        calcLocaldx(v_dxM, barR, uBound5, l , h);
     }
}


void FDUtils::addMoreSpacePtsConst(vector<double>& v_dxM, 
                                   int impNe, const double lB, const double uB, const int lN, const int uN, 
                            double barLowL, double barLowR, double barUpL, double barUpR){
    
    if (barLowR > barUpL){
        throw ModelException("FDUtils::addMoreSpacePtsAt2Crit", "two barriers points are too close!, not implemented yet!");
    }

    double oneDVolT = (barLowR - barLowL) /2.0;

    int n = uN - lN;

    double dxTemp = (uB - lB)/(n);
    //double dxTemp = (uB - lB)/(n + 2 * impNe);  //to think

    int impNeNew = int ((barLowR - barLowL) / dxTemp );

    //to check
    impNeNew= Maths::max(impNe, impNeNew);

    //have 2 crit pts
    int n2Bis = n - 2 *impNeNew;

    double wL = (barLowL - lB)/ (uB - lB - 4.0 * oneDVolT);
    double wM = (barUpL - barLowR)/ (uB - lB - 4.0 * oneDVolT);

    int n2L = int(wL * n2Bis);
    int n2M = int(wM * n2Bis);
    int n2R = n2Bis - n2L - n2M;

    int l ;
    int h ;

    //lB, barLOwL, n2L
    l = lN ;
    h =  l + n2L ;
    calcLocaldx(v_dxM, lB, barLowL, l, h);

    //barLowL, barLowR, impNeNew
    l = lN + n2L ;
    h =  l + impNeNew ;
    calcLocaldx(v_dxM, barLowL, barLowR, l, h);

    //barLOwR, barUpL, n2M
    l = lN + n2L + impNeNew;
    h =  l + n2M ;
    calcLocaldx(v_dxM, barLowR, barUpL, l, h);


    //barUpL, barUpR, impNeNew
    l = lN + n2L + impNeNew + n2M;
    h =  l + impNeNew ;
    calcLocaldx(v_dxM, barUpL, barUpR, l, h);

    //barUpR, uB, n2R
    l = lN + n2L + impNeNew + n2M + impNeNew;
    h =  l + n2R;
    calcLocaldx(v_dxM, barUpR, uB, l, h);
}


void FDUtils::calcLocaldx(vector<double>& v_dxM, double lowB, double upB, int lDim, int uDim){

    int i;
    int size = uDim - lDim;

    double dxTemp = (upB - lowB) / (size);

    for( i = lDim+1; i <= uDim; i++){
        v_dxM[i] = dxTemp;            
    }
}

//----------------------------------------------------------------------------

FDUtilsBarrier::FDUtilsBarrier() : CObject(TYPE){

    needSpecialFDUpBar = false;
    needSpecialFDDownBar = false;

    hasUpBarAtStep = false;
    hasDownBarAtStep = false;

    upBarrier = -1;
    valueAtUpBar = -1;

    downBarrier = -1;
    valueAtDownBar = -1;
}

FDUtilsBarrier::~FDUtilsBarrier(){

}

//Neighbour 
//mode = 1,  position of the smallest item larger than or equal to target 
//mode = -1,  position of the largest item smaller than or equal to target
void FDUtilsBarrier::calcClosestIndex(int& iDown, int& iUp, 
                    double* s, int sSize){

    int last = sSize -1;
    int mode = -1;

    hasDownBarAtStep = false;
    hasUpBarAtStep = false;

    if(downBarrier != -1){ // prod has overwritten the default value -1
       int iD = Neighbour(downBarrier, s, 0, last, mode);
        if (iD <= 0){ 
            //no barriers
            iD = 0;
        }else{
            iDown = iD;
            
            //can be solved in a smaller grids
            if(needSpecialFDDownBar == true){
                hasDownBarAtStep = true;
            }
        }
    }

    if(upBarrier !=-1){
        int iU = Neighbour(upBarrier, s, 0, last, mode);
        if ((iU <= 0) || (iU == iUp)){ //don't need special barrier treatment
            iUp = Maths::max(0, iU);
        }else{
            iUp = iU;

            //can be solved in a smaller grids
            if(needSpecialFDUpBar ==true){
                hasUpBarAtStep = true;
            }
        }
    }
    
    botDimBar = iDown;
    topDimBar = iUp;

    //if we only have below 3 points,
    //don't need to call solver
    if (iUp - iDown < 2){
        //to do
        //think where to control this and how to assign the value for thses points
        throw ModelException("FDUtilsBarrier::calcCloestIndex", "Barriers are too close !");
    }
}


void FDUtilsBarrier::adjustCBDs(const vector<double>& gridLevel1, 
                    int iDown, int iUp,  //these two index are the cloest index to the critical points
                    double* price, double* lastP) {

    double alpha;
    int size = gridLevel1.size() -1;
    double rebate ;

    if (hasDownBarAtStep == true) {
        //down side, should check s1[iDown + 1] - downBarrier <> 0??
        if (iDown +1 <= size) {
            rebate = valueAtDownBar;

//            alpha = (gridLevel1[iDown +1] - gridLevel1[iDown])/(gridLevel1[iDown + 1] - downBarrierGrid);
//            price[iDown] = alpha * rebate + (1 - alpha) * lastP[iDown + 1];

            alpha = (lastP[iDown + 1] - rebate) / (gridLevel1[iDown + 1] - downBarrierGrid );
            price[iDown] = rebate + alpha * (gridLevel1[iDown] - downBarrierGrid );
        }else{

        }
    }else{

    }

    if (hasUpBarAtStep == true){
//        if (iUp + 1 <= size) {
//            rebate =     valueAtUpBar;
//            //linear
//            alpha = (lastP[iUp + 1] - rebate) / (gridLevel1[iUp+1] - upBarrierGrid );
//            price[iUp] = rebate + alpha * (gridLevel1[iUp] - upBarrierGrid );
//        }

        //up side
        if (iUp -1 >= 0) {
            rebate =     valueAtUpBar;
            //linear
            alpha = (lastP[iUp - 1] - rebate) / (gridLevel1[iUp-1] - upBarrierGrid );
            price[iUp] = rebate + alpha * (gridLevel1[iUp] - upBarrierGrid );
        }
    }else{

    }   
}

// solve (u^n+1 - u^n)/da + (theta A^n u^n + (1-theta) A^n+1 u^n+1) + g^n+1 = 0
void FDUtils::euroOneStepWithSource(
    int bot, 
	int top, 
    SparseCollection aCurr, 
    SparseCollection aPrev,                       
    double * solCurr, 
    double * solPrev, 
    double * source,
	double interval,
	double theta)
{
	int iRow;

	double tmp1 = 1.;
	double tmp2 = (1. - theta) * interval;

	aPrev.linearTransform(tmp1, tmp2);
	aPrev.leftMultiplyTo(solPrev);

	for (iRow = bot; iRow <= top; ++iRow){
		source[iRow] = source[iRow] * interval + solPrev[iRow];
	}

	tmp2 = -theta * interval;
	aCurr.linearTransform(tmp1, tmp2);

	sparseLinearSystem(bot, top, aCurr, source, solCurr);
}

void FDUtils::sparseLinearSystem(
	int bot, 
	int top, 
	SparseCollection & coeff, 
	double * source, 
	double * solCurr)
{}

/*
void FDUtils::sparseLinearSystem(
	int bot, 
	int top, 
	SparseCollection & coeff, 
	double * source, 
	double * solCurr)
{
	// parameters ...
	int nin = 5;
	int nout = 6;
	int nmax = 5000;
	int liwork = 7 * nmax + 2;
	int lwork = 6 * nmax;

	// local scalars ...
	double			anorm = 0.0;
	double			sigmax = 0.0;
	double			stplhs, stprhs;
	double			dtol = 0.;
	double			tol = 1.0e-10;

	int				iterm = 1;
	int				lfill = 0;
	int				m = 1;
	int				maxitn = 20;
	int				monit = 1;
	int				n = top - bot + 1;
	int				nnz = coeff.size();
	int				la = 10 * nnz;
	int				i, ifail, ifail1, ifailx, irevcm, itn;
	int				lwreq, nnzc, npivm;

	bool			loop;

	char			weight[] = "N";
	int				length_weight = 1;

	char			method[] = "BICGSTAB";
	int				length_method = 8;

	char			norm[] = "1";
	int				length_norm = 1;

	char			precon[] = "P";
	int				length_precon = 1;

	char			milu[] = "N";
	int				length_milu = 1;

	char			pstrat[] = "C";
	int				length_pstrat = 1;

	// local arrays ...
	DoubleArray		a(la), b(nmax), wgt(nmax), work(lwork), x(nmax);
	IntArray		icol(la), idiag(nmax), ipivp(nmax), ipivq(nmax);
	IntArray		irow(la), istr(nmax+1), iwork(liwork);

	if (n < nmax) {
		for (i = 0; i < nnz; ++i){
			coeff.convert(i, irow[i], icol[i], a[i]);
			irow[i]++;
			icol[i]++;
		}

		for (i = 0; i < n; ++i){
			b[i] = source[bot+i];
			x[0] = source[bot+i];
		}

		// Calc incomplete LU factorization
		ifail = -1;
		F11DAF(n, nnz, &a[0], la, &irow[0], 
				&icol[0], lfill, dtol, pstrat, length_pstrat, 
				milu, length_milu, &ipivp[0], &ipivq[0], &istr[0], 
				&idiag[0], nnzc, npivm, &iwork[0], liwork, 
				ifail);

		// call F11BDF to initialize the solver
		ifail = -1;
		F11BDF(method, length_method, 
				precon, length_precon,
				norm, length_norm,
				weight, length_weight,
				iterm, n, m, tol, maxitn, 
				anorm, sigmax, monit, lwreq, &work[0], lwork, ifail);

		// call repeatedly F11BEF to solve the equations. 
		// Note that the arrays B and X are overwritten
		// On final exit, X will contain the solution and B the residual vector

		ifail = -1;
		irevcm = 0;
		loop = true;
		lwreq = lwork;
		
		while (loop) {
			F11BEF(irevcm, &x[0], &b[0], &wgt[0], &work[0], lwreq, ifail);
			if (irevcm == -1) {
				ifailx = -1;
				F11XAF("T", 1, n, nnz, &a[0], &irow[0], &icol[0], "N", 1, 
					&x[0], &b[0], ifailx);
			} else if (irevcm == 1) {
				ifailx = -1;
				F11XAF("N", 1, n, nnz, &a[0], &irow[0], &icol[0], "N", 1, 
					&x[0], &b[0], ifailx);
			} else if (irevcm == 2) {
				ifail1 = -1;
				F11DBF("N", 1, n, &a[0], la, &irow[0], &icol[0], &ipivp[0], &ipivq[0], 
					&istr[0], &idiag[0], "N", 1, &x[0], &b[0], ifail1);
				if (ifail1 != 0) irevcm = 6;
			} else if(irevcm == 3) {
				ifail1 = -1;
				F11BFF(itn, stplhs, stprhs, anorm, sigmax, &work[0], lwreq, ifail1);
			} else if (irevcm == 4) {
				loop = false;
			}
		}
		

		// Obtain information about the computation
		ifail1 = -1;
		F11BFF(itn, stplhs, stprhs, anorm, sigmax, &work[0], lwreq, ifail1);

		// Print the output data

	}            
	if (0 != ifail){
        throw ModelException("NAG routine F11BEF : Failed with IFAIL = " + Format::toString(ifail));
    } 

	for (i = 0; i < n; ++i){
		solCurr[bot+i] = x[i];
	}
}

class SparseLinearEquations: public CObject{
public:
	static CClassConstSP const TYPE;

	IObjectSP solve(){
		// parameters ...
		int nin = 5;
		int nout = 6;
		int nmax = 1000;
		int la = 10000;
		int liwork = 7 * nmax + 2;
		int lwork = 6 * nmax;

		// local scalars ...
		double			anorm = 0.;
		double			sigmax = 0.;
		double			stplhs, stprhs;

		char			weight[] = "N";
		int				length_weight = 1;

		char			method[] = "BICGSTAB";
		int				length_method = 8;

		char			norm[] = "1";
		int				length_norm = 1;

		char			precon[] = "P";
		int				length_precon = 1;

		char			milu[] = "N";
		int				length_milu = 1;

		char			pstrat[] = "C";
		int				length_pstrat = 1;

		int				i, ifail, ifail1, ifailx, irevcm, itn;
		int				lwreq, nnzc, npivm;

		bool			loop;

		// local arrays ...
		DoubleArray		a(la), b(nmax), wgt(nmax), work(lwork), x(nmax);
		IntArray		icol(la), idiag(nmax), ipivp(nmax), ipivq(nmax);
		IntArray		irow(la), istr(nmax+1), iwork(liwork);

		if (n < nmax) {

			for (i = 0; i < nnz; ++i){
				a[i] = aInput[i];
				irow[i] = irowInput[i];
				icol[i] = icolInput[i];
			}

			for (i = 0; i < n; ++i){
				b[i] = bInput[i];
				x[i] = 0.;
			}

			// Calc incomplete LU factorization
			ifail = 0;
			F11DAF(n, nnz, &a[0], la, &irow[0], 
					&icol[0], lfill, dtol, pstrat, length_pstrat, 
					milu, length_milu, &ipivp[0], &ipivq[0], &istr[0], 
					&idiag[0], nnzc, npivm, &iwork[0], liwork, 
					ifail);

			// call F11BDF to initialize the solver
			ifail = -1;
			F11BDF(method, length_method, 
					precon, length_precon,
					norm, length_norm,
					weight, length_weight,
					iterm, n, m, tol, maxitn, 
					anorm, sigmax, monit, lwreq, &work[0], lwork, ifail);

			// call repeatedly F11BEF to solve the equations. 
			// Note that the arrays B and X are overwritten
			// On final exit, X will contain the solution and B the residual vector

			ifail = 0;
			irevcm = 0;
			loop = true;
			lwreq = lwork;
			
			while (loop) {
				F11BEF(irevcm, &x[0], &b[0], &wgt[0], &work[0], lwreq, ifail);
				if (irevcm == -1) {
					ifailx = -1;
					F11XAF("Transpose", 9, n, nnz, &a[0], &irow[0], &icol[0], "No checking", 11, 
						&x[0], &b[0], ifailx);
				} else if (irevcm == 1) {
					ifailx = -1;
					F11XAF("No transpose", 12, n, nnz, &a[0], &irow[0], &icol[0], "No checking", 11, 
						&x[0], &b[0], ifailx);
				} else if (irevcm == 2) {
					ifail1 = -1;
					F11DBF("No transpose", 12, n, &a[0], la, &irow[0], &icol[0], &ipivp[0], &ipivq[0], 
						&istr[0], &idiag[0], "No checking", 11, &x[0], &b[0], ifail1);
					if (ifail1 != 0) irevcm = 6;
				} else if(irevcm == 3) {
					ifail1 = 0;
					F11BFF(itn, stplhs, stprhs, anorm, sigmax, &work[0], lwreq, ifail1);
				} else if (irevcm == 4) {
					loop = false;
				}
			}
			
			// Obtain information about the computation
			ifail1 = 0;
			F11BFF(itn, stplhs, stprhs, anorm, sigmax, &work[0], lwreq, ifail1);

			// Print the output data

		}            
		if (0 != ifail){
			throw ModelException("NAG routine F11BEF : Failed with IFAIL = " + Format::toString(ifail));
		} 

		DoubleArrayArray	obj(2);

		obj[0].resize(n);
		obj[1].resize(n);

		for (i = 0; i < n; ++i){
			obj[0][i] = x[i];
			obj[1][i] = b[i];
		}
		DoubleArrayArraySP	objSP(new DoubleArrayArray(obj));

		return objSP;
	}

private:
	int		n;
	int		iterm;
	int		m;
	int		maxitn;
	int		monit;
	int		lfill;
	int		nnz;

	double	tol;
	double	dtol;

	DoubleArray		aInput;
	IntArray		irowInput;
	IntArray		icolInput;

	DoubleArray		bInput;

    static void load(CClassSP& clazz){
        REGISTER(SparseLinearEquations, clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultCtor);
       
        FIELD(aInput, "value of non-zero element of sparse coeff matrix");        
		FIELD(irowInput, "row of non-zero element of sparse coeff matrix");   
		FIELD(icolInput, "col of non-zero element of sparse coeff matrix");   
        FIELD(bInput, "source");
              
        Addin::registerClassObjectMethod("SPARSE_LINEAR_EQUATIONS",
                                          Addin::UTILITIES,
                                          "Solve linear equations, with sparse coefficient matrix",
                                          TYPE,
                                          true,
                                          Addin::returnHandle,
                                          (Addin::ObjMethod*)output);
    }

	static IObjectSP output(SparseLinearEquations* params) {
        return params->solve();
    }

	static IObject* defaultCtor(){
        return new SparseLinearEquations();
    }
	
    SparseLinearEquations():
    CObject(TYPE){}

	virtual void validatePop2Object(){
        static const string method = "SparseLinearEqations::validatePop2Object";
        try{
			n = bInput.size();
			if(aInput.size()> n * n){
                throw ModelException(method,
                                     "size of aInput can not be larger than square of size of bInput.");
            }

			nnz = aInput.size();
			if(nnz != irowInput.size()){
                throw ModelException(method,
                                     "size of aInput equals to size of irowInput.");
            }

			if(nnz != icolInput.size()){
                throw ModelException(method,
                                     "size of aInput equals to size of icolInput.");
            }

			int i;
			for (i = 1; i < irowInput.size(); ++i){
				if (irowInput[i-1] > irowInput[i]) {
					throw ModelException(method,
										"value of irowInput vector should be increasing.");
				}
			}

			for (i = 1; i < icolInput.size(); ++i){
				if ((icolInput[i-1] > icolInput[i])&&(irowInput[i-1] == irowInput[i] )){
					throw ModelException(method,
										"value of icolInput vector should be increasing.");
				}
			}

			iterm	= 1;
			m		= 1;
			maxitn	= 20;
			monit	= 1;
			lfill	= 0;

			dtol	= 0.;
			tol		= 1.e-8;
	
        } catch (exception& e){
            throw ModelException(e, method);
        }
    }
};

CClassConstSP const SparseLinearEquations::TYPE =
			CClass::registerClassLoadMethod("SparseLinearEquations", typeid(SparseLinearEquations), load);
*/
DRLIB_END_NAMESPACE

