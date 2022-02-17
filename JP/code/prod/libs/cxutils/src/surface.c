/*
***************************************************************************
** SOURCE FILE: surface.c
**
** Generic surface functions.
**
** $Header$
***************************************************************************
*/

#include "surface.h"

#include "bsearch.h"
#include "cxmacros.h"

#include <alib/lintrp.h>
#include <alib/sintrp.h>

#include <alib/cversion.h>

#include <math.h>

static void spliced(int wantY, int wantSlope, int wantCurve, int n,
                    const double *xData, const double *yData,
                    double gamma1, double gamman, int npoints, 
                    const double *x, int *err,
                    double *y, double *slope, double *curve);

#undef ELEMENT
#define ELEMENT(array,skip,i) (*((double*)(((char*)array) + (i)*(skip))))

/*f
***************************************************************************
** Interpolation routine for base correlation.
**
** This is a spline internally but with strange behaviour around the
** edges. Extrapolation requests are handled by linear extrapolation.
***************************************************************************
*/
static int BaseCorrInterpDoublePoint1
(double           *x,                /* (I) Ordered Array of X values */
 int               xskip,            /* (I) # bytes between x values */
 int               N,                /* (I) Length of X & F arrays */
 double           *f,                /* (I) Ordered Array of F values */
 int               fskip,            /* (I) # bytes between f values */
 double            xDesired,         /* (I) X for which F is desired */
 TMetricDoubleFunc mfunc,            /* (I) Metric Function */
 double           *fInterp)          /* (O) Interpolated value */
{
    static char routine[] = "BaseCorrInterpDoublePoint1";
    int         status    = FAILURE;

    double xLinInterp[2];
    double fLinInterp[2];

    TBoolean linInterp = FALSE;
    int      splineStart = 0;

    REQUIRE (N > 0);
    REQUIRE (x != NULL);
    REQUIRE ((size_t)xskip >= sizeof(double));
    REQUIRE (f != NULL);
    REQUIRE ((size_t)fskip >= sizeof(double));
    REQUIRE (fInterp != NULL);

    if (ELEMENT(x,xskip,0) < 0.0)
    {
        /* linearly interpolate if xDesired <= x[2] or if xDesired >= x[N-1] */
        if (N >= 2 &&
            (xDesired <= ELEMENT(x,xskip,2) ||
             xDesired >= ELEMENT(x,xskip,N-1)))
        {
            linInterp = TRUE;
            if (xDesired <= ELEMENT(x,xskip,2))
            {
                xLinInterp[0] = 0.0;
                xLinInterp[1] = ELEMENT(x,xskip,2);
                fLinInterp[0] = ELEMENT(f,fskip,0);
                fLinInterp[1] = ELEMENT(f,fskip,2);
            }
            else
            {
                ASSERT (xDesired >= ELEMENT(x,xskip,N-1));
                xLinInterp[0] = ELEMENT(x,xskip,N-2);
                xLinInterp[1] = ELEMENT(x,xskip,N-1);
                fLinInterp[0] = ELEMENT(f,fskip,N-2);
                fLinInterp[1] = ELEMENT(f,fskip,N-1);
            }
        }
        else
        {
            splineStart = 1;
            linInterp   = FALSE;
        }
    }
    else
    {
        /* linear extrapolation */
        if (xDesired <= ELEMENT(x,xskip,0))
        {
            linInterp = TRUE;
            xLinInterp[0] = ELEMENT(x,xskip,0);
            xLinInterp[1] = ELEMENT(x,xskip,1);
            fLinInterp[0] = ELEMENT(f,fskip,0);
            fLinInterp[1] = ELEMENT(f,fskip,1);
        }
        else if (xDesired >= ELEMENT(x,xskip,N-1))
        {
            linInterp = TRUE;
            xLinInterp[0] = ELEMENT(x,xskip,N-2);
            xLinInterp[1] = ELEMENT(x,xskip,N-1);
            fLinInterp[0] = ELEMENT(f,fskip,N-2);
            fLinInterp[1] = ELEMENT(f,fskip,N-1);
        }
        else
        {
            linInterp = FALSE;
            splineStart = 0;
        }
    }

    if (linInterp)
    {
        if (GtoLinInterpDoublePoint1 (xLinInterp,
                                      sizeof(double),
                                      2,
                                      fLinInterp,
                                      sizeof(double),
                                      xDesired,
                                      mfunc,
                                      fInterp) != SUCCESS)
            goto done; /* failure */
    }
    else
    {
        if (GtoSplineInterpDoublePoint1 (&(ELEMENT(x,xskip,splineStart)),
                                         xskip,
                                         N-splineStart,
                                         &(ELEMENT(f,fskip,splineStart)),
                                         fskip,
                                         xDesired,
                                         mfunc,
                                         fInterp) != SUCCESS)
            goto done; /* failure */
    }

    status = SUCCESS;

done:
    
    if (status != SUCCESS)
        GtoErrMsgFailure (routine);

    return status;
}

/*f
***************************************************************************
** Interpolation routine using spliced parabolas.
**
** For a particular point this involves averaging two parabolas. Consider
** four points - x1,x2,x3,x4 - and we are interested in a point between
** x2 and x3. Then we find the parabola that fits x1,x2,x3 and we find the
** parabola that fits x2,x3,x4. Then we linearly interpolate the result
** based on the distance between x2 and x3. Thus near point x2 we are close
** to using the x1,x2,x3 parabola and near point x3 we are close to using
** the x2,x3,x4 parabola.
**
** Gives a degree of smoothness but also guarantees relatively local
** interpolation.
**
** Extrapolation requests are an error.
***************************************************************************
*/
static int GtoSpliceInterpDoublePoint1
(double           *x,                /* (I) Ordered Array of X values */
 int               xskip,            /* (I) # bytes between x values */
 int               N,                /* (I) Length of X & F arrays */
 double           *f,                /* (I) Ordered Array of F values */
 int               fskip,            /* (I) # bytes between f values */
 double            xDesired,         /* (I) X for which F is desired */
 TMetricDoubleFunc mfunc,            /* (I) Metric Function */
 double           *fInterp)          /* (O) Interpolated value */
{
    static char routine[] = "GtoSpliceInterpDoublePoint1";
    int         status    = FAILURE;

    double     *xData = NULL;
    double     *yData = NULL;
    
    int         i;
    int         err;
    double      loSlope;
    double      hiSlope;

    REQUIRE (N > 0);
    REQUIRE (x != NULL);
    REQUIRE ((size_t)xskip >= sizeof(double));
    REQUIRE (f != NULL);
    REQUIRE ((size_t)fskip >= sizeof(double));
    REQUIRE (fInterp != NULL);

    /* relax this requirement */
    REQUIRE (mfunc == NULL);

    /* call the spliced routine */
    xData = NEW_ARRAY(double, N);
    yData = NEW_ARRAY(double, N);

    for (i = 0; i < N; ++i)
    {
        xData[i] = *((double*)(((char*)x) + i*xskip));
        yData[i] = *((double*)(((char*)f) + i*fskip));
    }

    if (N > 1)
    {
        loSlope = (yData[1] - yData[0]) / (xData[1] - xData[0]);
        hiSlope = (yData[N-1] - yData[N-2]) / (xData[N-1] - xData[N-2]);
    }
    else
    {
        loSlope = 0.0;
        hiSlope = 0.0;
    }

    /* for the moment do not allow extrapolation */
    if (xDesired < xData[0] || xDesired > xData[N-1])
    {
        GtoErrMsg ("%s: Requested value %f is not in range [%f,%f]\n",
                 routine, xDesired, xData[0], xData[N-1]);
        goto done; /* failure */
    }

    spliced (1, 0, 0, N, xData, yData, loSlope, hiSlope, 1, &xDesired, &err,
             fInterp, NULL, NULL);
    if (err)
    {
        GtoErrMsg ("%s: Error in splice interpolation\n", routine);
        goto done; /* failure */
    }

    status = SUCCESS;

done:
    
    FREE (xData);
    FREE (yData);

    if (status != SUCCESS)
        GtoErrMsgFailure (routine);

    return status;
}

/*
    NAME:       spliced.c 
    SUMMARY:    interpolation method using spliced parabolas. 
    AUTHOR:     G. A. Gatoff.
    HISTORY:    created on 23rd August 1995.
                Added slope and curve outputs on March 4th, 1996.
                May 2003: merged with fourWeights.c and added getSpliceWeights

    RETURNS:    void.

    CAUTION:    None.

    INPUTS:

            wantY        -     1 if you want y(x) interpolated, else 0.
            wantSlope    -     1 if you want y'(x) interpolated, else 0.
            wantCurve    -     1 if you want y''(x) interpolated, else 0.
            n            -    number of points in data
            xData          -    array of size n (x_1,... x_n)
            yData        -    array of size n (y_1,... y_n)
            gamma1        -    slope for x < x1
            gamman        -    slope for x >= xn
            npoints     -    number of interpolted points.
            x            -    x points for which y needs interpolating.
    OUTPUTS:
            err         -    error number 
            y            -    interpolated y points
            slope         -    interploated slopes (deltas)
            curve         -   interpolated curves (gammas)

*/

#define REPETITION_ERR 1
#define SORTING_ERR    2

static void spliced(int wantY, int wantSlope, int wantCurve, int n, 
                    const double *xData, const double *yData,
                    double gamma1, double gamman, 
                    int npoints, const double *x, 
                    int *err, double *y, double *slope, double *curve)
{
    int i,j;
    double pj=0.0,pjPlus1=0.0,pjdash=0.0;
    double pjPlus1dash=0.0,pjdashdash=0.0,pjPlus1dashdash=0.0;
    double x1 = xData[0];
    double x2 = xData[1];
    double x3 = xData[2];
    double xnMin2 = xData[n-3];
    double xnMin1 = xData[n-2];
    double xn    = xData[n-1];
    double y1 = yData[0];
    double y2 = yData[1];
    double y3 = yData[2];
    double ynMin2 = yData[n-3];
    double ynMin1 = yData[n-2];
    double yn    = yData[n-1];
    double A1,B1,An,Bn,p2dash,pnMin1dash;        

    int errCode = 0;
    
    /* checking input data */
    
    for (j=1 ; j < n ; j++)
    {
        if (IS_EQUAL(xData[j], xData[j-1]))
        {
            errCode = REPETITION_ERR;
            goto done; /* failure */
        }
        if (xData[j] < xData[j-1])
        {
            errCode = SORTING_ERR;
            goto done; /* failure */
        }
    }

    for (i=0; i<npoints ; i++)
    {
        
        if      (x[i]  < x1)
        {
            if (wantY)             y[i] =    y1-gamma1*(x1-x[i]);
            if (wantSlope)    slope[i] =     gamma1;
            if (wantCurve)     curve[i] =     0.0;
        }
        else if (xn <= x[i])
        {
            if (wantY)          y[i] =    yn+gamman*(x[i]- xn);
            if (wantSlope)     slope[i] =    gamman;
            if (wantCurve)     curve[i] =    0.0;
        }
        else if ((x1 <= x[i]) && (x[i] < x2))
        {
            p2dash = (x2-x3)*y1/((x1-x2)*(x1-x3)) + y2/(x2-x1)
                + y2/(x2-x3) + (x2-x1)*y3/ ((x3-x1)*(x3-x2));
            A1    =    ((x2-x1)*(p2dash+gamma1)-2*(y2-y1))/pow((x2-x1),3.0);
            B1    =    (3*(y2-y1)-(x2-x1)*
                        (p2dash+2*gamma1))/((x2-x1)*(x2-x1));
            if (wantY)
                y[i]    =    (x[i]-x1)*(x[i]-x1)*(x[i]- x1)*A1+
                    (x[i]-x1)*(x[i]-x1)*B1+(x[i]-x1)*gamma1+y1;
            if (wantSlope)
                slope[i] = 3*(x[i]-x1)*(x[i]-x1)*A1+2*(x[i]-x1)*B1+gamma1;
            if (wantCurve)
                curve[i] = 6*(x[i]-x1)*A1+2*B1;
            
        }
        else if ((xnMin1 <= x[i]) && (x[i] < xn))
        {
            pnMin1dash = 
                (xnMin1-xn)*ynMin2/((xnMin2-xnMin1)*(xnMin2-xn))
                + ynMin1/(xnMin1-xnMin2) + ynMin1/(xnMin1-xn) 
                + (xnMin1-xnMin2)*yn /((xn-xnMin2)*(xn-xnMin1));
            An    =    ((xn-xnMin1)*(pnMin1dash+gamman)-2*(yn-ynMin1))/
                pow((xn-xnMin1),3.0);
            Bn    =    ((xn-xnMin1)*(pnMin1dash+2*gamman)-3*(yn-ynMin1))/
                ((xn-xnMin1)*(xn-xnMin1));
            if (wantY)
                y[i]    =    (x[i]-xn)*(x[i]-xn)*(x[i]-xn)*An+
                    (x[i]-xn)*(x[i]-xn)*Bn+(x[i]-xn)*gamman+yn;
            if (wantSlope)
                slope[i] = 3*(x[i]-xn)*(x[i]-xn)*An+2*(x[i]-xn)*Bn+gamman;
            if (wantCurve)
                curve[i] = 6*(x[i]-xn)*An+2*Bn;
            
        }
        else 
        {
            for (j=1; j<(n-2); j++)
            {
                if ((xData[j] <= x[i])    && (x[i] < xData[j+1])) 
                {
                    break;
                }
            }
            
            if ((wantY) || (wantSlope))
            {
                pj     =      ((x[i]- xData[j])*(x[i]- xData[j+1]))*yData[j-1]/
                    ((xData[j-1]-xData[j])*(xData[j-1]-xData[j+1]))    
                    +((x[i]- xData[j-1])*(x[i]- xData[j+1]))*yData[j]/
                    ((xData[j]-xData[j-1])*(xData[j]-xData[j+1]))
                    +((x[i]- xData[j-1])*(x[i]- xData[j]))*yData[j+1]/
                    ((xData[j+1]-xData[j-1])*(xData[j+1]-xData[j]));
                pjPlus1 = ((x[i]- xData[j+1])*(x[i]- xData[j+2]))*yData[j]/
                    ((xData[j]-xData[j+1])*(xData[j]-xData[j+2]))    
                    +((x[i]- xData[j])*(x[i]- xData[j+2]))*yData[j+1]/
                    ((xData[j+1]-xData[j])*(xData[j+1]-xData[j+2]))
                    +((x[i]- xData[j])*(x[i]- xData[j+1]))*yData[j+2]/
                    ((xData[j+2]-xData[j])*(xData[j+2]-xData[j+1]));
            }
            if ((wantSlope) || (wantCurve))
            {
                pjdash         =    
                    (2*x[i]-xData[j]-xData[j+1])*yData[j-1]/
                    ((xData[j-1]-xData[j])*(xData[j-1]-xData[j+1]))    
                    +(2*x[i]-xData[j-1]-xData[j+1])*yData[j]/
                    ((xData[j]-xData[j-1])*(xData[j]-xData[j+1]))
                    +(2*x[i]-xData[j-1]-xData[j])*yData[j+1]/
                    ((xData[j+1]-xData[j-1])*(xData[j+1]-xData[j]));
                pjPlus1dash  =
                    (2*x[i]- xData[j+1]- xData[j+2])*yData[j]/
                    ((xData[j]-xData[j+1])*(xData[j]-xData[j+2]))    
                    +(2*x[i]- xData[j]-xData[j+2])*yData[j+1]/
                    ((xData[j+1]-xData[j])*(xData[j+1]-xData[j+2]))
                    +(2*x[i]- xData[j]-xData[j+1])*yData[j+2]/
                    ((xData[j+2]-xData[j])*(xData[j+2]-xData[j+1]));
            }
            
            if (wantCurve)
            {
                pjdashdash         =
                    2*yData[j-1]/
                    ((xData[j-1]-xData[j])*(xData[j-1]-xData[j+1]))    
                    +2*yData[j]/
                    ((xData[j]-xData[j-1])*(xData[j]-xData[j+1]))
                    +2*yData[j+1]/
                    ((xData[j+1]-xData[j-1])*(xData[j+1]-xData[j]));
                pjPlus1dashdash    =
                    2*yData[j]/
                    ((xData[j]-xData[j+1])*(xData[j]-xData[j+2]))    
                    +2*yData[j+1]/
                    ((xData[j+1]-xData[j])*(xData[j+1]-xData[j+2]))
                    +2*yData[j+2]/
                    ((xData[j+2]-xData[j])*(xData[j+2]-xData[j+1]));
                
            }
            if (wantY)     
                y[i] =      (xData[j+1]-x[i])*pj/(xData[j+1]-xData[j])
                    + (x[i]- xData[j])*pjPlus1/(xData[j+1]-xData[j]);
            if (wantSlope)
                slope[i] = -pj/(xData[j+1]-xData[j])
                    + pjPlus1/(xData[j+1]-xData[j])
                    + (xData[j+1]-x[i])*pjdash/(xData[j+1]-xData[j])
                    +(x[i]- xData[j])*pjPlus1dash/(xData[j+1]-xData[j]);
            if (wantCurve)
                curve[i] = -2*pjdash/(xData[j+1]-xData[j])
                    + 2*pjPlus1dash/(xData[j+1]-xData[j])
                    + (xData[j+1]-x[i])*pjdashdash/(xData[j+1]-xData[j])
                    + (x[i]- xData[j])*pjPlus1dashdash/
                    (xData[j+1]-xData[j]);
        }
    }
                    
 done:
    
    *err = errCode;
    return;
}



/*f
***************************************************************************
** Constructor for CXSurface
***************************************************************************
*/
CxTSurface* CxSurfaceMakeFromInputs(
int             numItems,            /* (I) */
double*         x,                   /* (I) [numItems] */
double*         y,                   /* (I) [numItems] */
double*         z,                   /* (I) [numItems] */
TBoolean*       include,             /* (I) [numItems] optional */
TBoolean        xRegular,
TBoolean        yRegular)
{
    static char routine[] = "CxSurfaceMakeFromInputs";
    int         status    = FAILURE;
    
    CxTSurface  *surface = NULL;
    double     *xIncluded = NULL;
    double     *yIncluded = NULL;
    double     *zIncluded = NULL;
    CxTSurfaceRegularity regularity;

    int numIncludes = 0;
    int i;
    int j;

    REQUIRE (numItems > 0);
    REQUIRE (x != NULL);
    REQUIRE (y != NULL);
    REQUIRE (z != NULL);

    /* TBD: relax these requirements */
    REQUIRE (xRegular);
    REQUIRE (!yRegular);

    /* so the first step is to only include what we need */
    if (include == NULL)
    {
        /* include everything */
        xIncluded = NEW_ARRAY (double, numItems);
        yIncluded = NEW_ARRAY (double, numItems);
        zIncluded = NEW_ARRAY (double, numItems);

        COPY_ARRAY (xIncluded, x, double, numItems);
        COPY_ARRAY (yIncluded, y, double, numItems);
        COPY_ARRAY (zIncluded, z, double, numItems);

        numIncludes = numItems;
    }
    else
    {
        /* only include what was requested */
        for (i = 0; i < numItems; ++i)
        {
            if (include[i])
                numIncludes++;
        }
        
        if (numIncludes == 0)
        {
            GtoErrMsg ("%s: No items requested for inclusion\n", routine);
            goto done; /* failure */
        }

        xIncluded = NEW_ARRAY (double, numIncludes);
        yIncluded = NEW_ARRAY (double, numIncludes);
        zIncluded = NEW_ARRAY (double, numIncludes);

        for (i = 0, j = 0; i < numItems; ++i)
        {
            if (include[i])
            {
                xIncluded[j] = x[i];
                yIncluded[j] = y[i];
                zIncluded[j] = z[i];
                ++j;
            }
        }
        ASSERT (j == numIncludes);
    }

    if (xRegular && yRegular)
        regularity = CX_SURFACE_REGULAR_GRID;
    else if (xRegular)
        regularity = CX_SURFACE_REGULAR_X;
    else if (yRegular)
        regularity = CX_SURFACE_REGULAR_Y;
    else
        regularity = CX_SURFACE_IRREGULAR;

    /* TBD: add some tedious sorting routines later */
    /* for grid - sort by x and y - validate that all the y's are the same */
    /* for x-regular - sort by x then y */
    /* for y-regular - sort by y then x */
    /* for irregular - sort by x then y - but perhaps we don't need to sort */

    /* now call the real constructor */
    surface = CxSurfaceMake (numIncludes,
                             xIncluded,
                             yIncluded,
                             zIncluded,
                             regularity);
    if (surface == NULL)
        goto done; /* failure */

    status = SUCCESS;

 done:

    FREE (xIncluded);
    FREE (yIncluded);
    FREE (zIncluded);

    if (status != SUCCESS)
    {
        GtoErrMsgFailure (routine);
        CxSurfaceFree (surface);
        surface = NULL;
    }
    
    return surface;
}



/*
***************************************************************************
** Interpolate a single value from a surface
***************************************************************************
*/
int CxSurfaceInterp(
    CxTSurface*          surface,
    CxTSurfaceInterp*    interp,
    double               x,
    double               y,
    double              *result)
{
    static char routine[] = "CxSurfaceInterp";
    int         status    = FAILURE;

    double      z = 0.0;

    long        exact;
    long        loBound;
    long        hiBound;

    typedef int (*INTERP_FUNC) (
        double           *x,                /* (I) Ordered Array of X values */
        int               xskip,            /* (I) # bytes between x values */
        int               N,                /* (I) Length of X & F arrays */
        double           *f,                /* (I) Ordered Array of F values */
        int               fskip,            /* (I) # bytes between f values */
        double            xDesired,         /* (I) X for which F is desired */
        TMetricDoubleFunc mfunc,            /* (I) Metric Function */
        double           *fInterp);         /* (O) Interpolated value */

    INTERP_FUNC interpFunc;

    REQUIRE (surface != NULL);
    REQUIRE (interp != NULL);
    REQUIRE (result != NULL);

    /* TBD - relax requirements */
    REQUIRE (surface->regularity == CX_SURFACE_REGULAR_X);
    REQUIRE (interp->interpSurface == FALSE);
    REQUIRE (interp->xInterp == CX_SURFACE_INTERP_LINEAR);
    REQUIRE (interp->xExtrap == CX_SURFACE_EXTRAP_FLAT);
    REQUIRE (interp->yInterp == CX_SURFACE_INTERP_SPLINE ||
             interp->yInterp == CX_SURFACE_INTERP_SPLICE ||
             interp->yInterp == CX_SURFACE_INTERP_LINEAR ||
             interp->yInterp == CX_SURFACE_INTERP_BASE_CORR);
    REQUIRE (interp->yExtrap == CX_SURFACE_EXTRAP_NONE);

    /* given our rather restricted set of inputs, this becomes very simple */

    /* binary search on x */
    /* spline interp on y for lo and hi values of x */
    /* linear interp on x */

    if (CxBinarySearchDouble (x,
                              surface->x,
                              sizeof(double),
                              surface->numItems,
                              &exact,
                              &loBound,
                              &hiBound) != SUCCESS)
        goto done; /* failure */

    switch (interp->yInterp)
    {
    case CX_SURFACE_INTERP_SPLINE:
        interpFunc = GtoSplineInterpDoublePoint1;
        break;
    case CX_SURFACE_INTERP_SPLICE:
        interpFunc = GtoSpliceInterpDoublePoint1;
        break;
    case CX_SURFACE_INTERP_LINEAR:
        interpFunc = GtoLinInterpDoublePoint1;
        break;
    case CX_SURFACE_INTERP_BASE_CORR:
        interpFunc = BaseCorrInterpDoublePoint1;
        break;
    default: 
        ASSERT(0);
    }
        
    if (exact >= 0 || hiBound >= surface->numItems || loBound < 0)
    {
        if (exact < 0)
        {
            double tmp;

            if (hiBound >= surface->numItems)
                tmp = surface->x[surface->numItems-1];
            else
                tmp = surface->x[0];

            if (CxBinarySearchDouble (tmp,
                                      surface->x,
                                      sizeof(double),
                                      surface->numItems,
                                      &exact,
                                      &loBound,
                                      &hiBound) != SUCCESS)
                goto done; /* failure */
            ASSERT (exact >= 0);
        }

        /* no linear interpolation on x required */
        if (interpFunc (surface->y + loBound + 1,
                        sizeof(double),
                        hiBound - loBound - 1,
                        surface->z + loBound + 1,
                        sizeof(double),
                        y,
                        NULL,
                        &z) != SUCCESS)
            goto done; /* failure */
    }
    else
    {
        double xPrev;
        double xNext;
        double zPrev;
        double zNext;

        ASSERT (hiBound < surface->numItems);
        ASSERT (loBound >= 0);

        xPrev = surface->x[loBound];
        xNext = surface->x[hiBound];

        if (CxBinarySearchDouble (xPrev,
                                  surface->x,
                                  sizeof(double),
                                  surface->numItems,
                                  &exact,
                                  &loBound,
                                  &hiBound) != SUCCESS)
            goto done; /* failure */
        ASSERT (exact >= 0);

        if (interpFunc (surface->y + loBound + 1,
                        sizeof(double),
                        hiBound - loBound - 1,
                        surface->z + loBound + 1,
                        sizeof(double),
                        y,
                        NULL,
                        &zPrev) != SUCCESS)
            goto done; /* failure */
        
        if (CxBinarySearchDouble (xNext,
                                  surface->x,
                                  sizeof(double),
                                  surface->numItems,
                                  &exact,
                                  &loBound,
                                  &hiBound) != SUCCESS)
            goto done; /* failure */
        ASSERT (exact >= 0);

        if (interpFunc (surface->y + loBound + 1,
                        sizeof(double),
                        hiBound - loBound - 1,
                        surface->z + loBound + 1,
                        sizeof(double),
                        y,
                        NULL,
                        &zNext) != SUCCESS)
            goto done; /* failure */

        z = zPrev + (x - xPrev) * (zNext - zPrev) / (xNext - xPrev);
    }

    status = SUCCESS;
    *result = z;

 done:

    if (status != SUCCESS)
        GtoErrMsgFailure (routine);

    return status;
}

/*
***************************************************************************
** Interpolate a matrix from a surface
***************************************************************************
*/
TMatrix2D* CxSurfaceInterpMatrix(
    CxTSurface*           surface,
    CxTSurfaceInterp*     interp,
    int                  xSize,
    double              *x,
    int                  ySize,
    double              *y)
{
    static char routine[] = "CxSurfaceInterpMatrix";
    int         status    = FAILURE;

    TMatrix2D  *matrix = NULL;
    int i;
    int j;

    REQUIRE (xSize > 0);
    REQUIRE (x != NULL);
    REQUIRE (ySize > 0);
    REQUIRE (y != NULL);

    matrix = GtoMatrixNewEmpty (xSize, ySize);
    if (matrix == NULL)
        goto done; /* failure */

    for (i = 0; i < xSize; ++i)
    {
        for (j = 0; j < ySize; ++j)
        {
            double value;
            if (CxSurfaceInterp (surface, interp, x[i], y[j], &value))
                goto done; /* failure */

            matrix->data[i][j] = value;
        }
    }

    status = SUCCESS;

 done:

    if (status != SUCCESS)
    {
        GtoErrMsgFailure (routine);
        GtoMatrixFree (matrix);
        matrix = NULL;
    }

    return matrix;
}

/*
***************************************************************************
** Interpolate a data table from a surface
***************************************************************************
*/
TDataTable* CxSurfaceInterpTable(
    CxTSurface*           surface,
    CxTSurfaceInterp*     interp,
    int                  xSize,
    double              *x,
    char               **xHeading,
    int                  ySize,
    double              *y,
    char               **yHeading)
{
    static char routine[] = "CxSurfaceInterpTable";
    int         status    = FAILURE;

    TMatrix2D  *matrix = NULL;
    TDataTable *table = NULL;

    REQUIRE (xHeading != NULL);
    REQUIRE (yHeading != NULL);

    matrix = CxSurfaceInterpMatrix (surface, interp, xSize, x, ySize, y);
    if (matrix == NULL)
        goto done; /* failure */

    table = GtoDataTableMake (matrix, xHeading, yHeading);
    if (table == NULL)
        goto done; /* failure */

    status = SUCCESS;

 done:

    GtoMatrixFree (matrix);
    if (status != SUCCESS)
    {
        GtoErrMsgFailure (routine);
        GtoDataTableDelete (table);
        table = NULL;
    }

    return table;
}
