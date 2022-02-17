//----------------------------------------------------------------------------
//
//   Group       : Derivatives Research
//
//   Filename    : CCMLossUnit.cpp
//
//   Description : Calculates loss unit
//
//   Date        : Aug 2004
//
//
//----------------------------------------------------------------------------
#include "edginc/config.hpp"
#include "edginc/CCMLossUnit.hpp"

DRLIB_BEGIN_NAMESPACE

#define MAX_DISCRETISATION_UNITS 100000

/* Euclide's algorithm for Greatest Common Denominator */
static int GCD(int x,int y )
{
    while( x != 0 ) {
        int oldY = y;
        y = x;
        x = oldY % x;
    }
    return y;
}

/* Least Common Multiple */
static int LCM(int x, int y)
{
    return x/GCD(x,y) * y;
}

/* express x as a simple continued fraction a0+ 1/(a1+1/(a2+1/...))) */
static void simpleContinuedFraction(
    double    x,    /* (I) number to be decomposed */
    int       nmax, /* (I) max size of the continued fraction expansion */
    IntArray& a){   /* (O) array[n] with the resulting decomposition    */
    a.reserve(nmax);
    double r=x;
    ASSERT(r>=0.);
    for(int i=0; i<nmax ; ++i) {
        int b = (int)floor(r);
        a.push_back(b);
        if (fabs(r-b)<1.e-9 ) break;
        r = 1./(r-b);
    }
}

#define MAX_CONTINUED_FRACTION 100

/* express x as a fraction A/B where A<Nmax */
static void rationalConvergent(
    double x, 
    int&   A, 
    int&   B, 
    int    maxDenominator){
    IntArray a;
    simpleContinuedFraction(x, 100, a);

    int Am1, Am2, Bm1, Bm2, i;
    A = a[0], B = 1; /* handle x=int */
    for (i = 1,Am2=1,Am1=a[0],Bm2=0,Bm1=1; i < a.size();
         i++, Bm2=Bm1, Am2=Am1, Am1=A, Bm1=B) {
        A = a[i]*Am1 + Am2;
        B = a[i]*Bm1 + Bm2;
        if (B > maxDenominator) {
            A = Am1;
            B = Bm1;
            break;
        }
    }
}

/** calculate loss unit for a given array of loss amt */
static void lossUnit(
    const DoubleArray&  l,         /* (I) loss amount        */
    int                 maxDenom,  /* (I) max denominator        */
    double&             lossUnit,  /* (O) loss unit              */
    double&             sum){      /* (O) sum of loss amt        */
    const char routine[] = "lossunit";
    try{
        int nl = l.size(); // for ease
        sum = 0.;
        int j;
        for (j=0; j<nl; ++j) {
            ASSERT(l[j] >= 0.);
            sum += l[j];
        }
        
        /* bail out with failure if all names have 0 notional (all default) */
        if (sum <= 0.){
            throw ModelException(routine, "All names have 0 notional");
        }
        
        int m = 1;
        for (j=0; j<nl; ++j) {
            int p,q;
            rationalConvergent(l[j]/sum,p,q,maxDenom);
            m = LCM(m,q);
            if (m > maxDenom) {
                lossUnit = sum / maxDenom;
                return;
            }
        }
        
        lossUnit = sum / m;
    } catch (exception& e){
        throw ModelException(e, routine);
    }
}
#if 0
/** we think this might be handy for tranche of tranche */
int ccmLossunitMultiArray(
    long           nv, /* (I) nb vector              */
    const long    *nl, /* (I) nb loss amt[nv]        */
    const double **l,  /* (I) loss amount[nv][nl[nv]]*/
    long    maxDenom,  /* (I) max denominator        */
    double *lossUnit   /* (O) loss unit              */
    )
{
    const char routine[] = "ccmLossunitMultiArray";
    int status = -1;
    int i;
    double sumMax, tmp;
    double *lu  = NULL;
    double *sum = NULL;
    MALLOC(double, lu,  nv);
    MALLOC(double, sum, nv);

    sumMax = -1;
    for (i=0; i<nv; ++i)
    {
        if (ccmLossunit(nl[i],l[i],maxDenom, &lu[i], &sum[i]) != SUCCESS) 
            goto done;
        sumMax = MAX(sumMax,sum[i]);
    }

    if (ccmLossunit(nv,lu,maxDenom, lossUnit, &tmp) != SUCCESS) 
        goto done;

    if (sumMax / *lossUnit > maxDenom)
        *lossUnit = sumMax / maxDenom;
    status = SUCCESS;
done:
    if (status != SUCCESS) ccmError("%s: failed", routine);
    FREE(lu,free);
    return status;
}
#endif

/**
 * Calculate loss unit (bin size for loss distribution discretization)
 * lossAmt and lossAmtCata must already be allocated with size n
 * if a name is defaulted, it is still in the lossAmt list with amt 0.
 */
double CCMLossUnit::calcLossUnit(
    const DoubleArray& notionals,    /* (I) name notionals */
    int                maxSlice,     /* (I) max nb of slice              */
    int                subdivision){ /* (I) subdivion of the GCD */
    const char routine[] = "CCMLossUnit::calcLossUnit";
    try{
        int n = notionals.size();
        DoubleArray l(n);
        double total = 0.;
        for (int i=0; i<n; ++i) {
            l[i] = fabs(notionals[i]);
            total += l[i];
        }
        double tmp;
        double result;
        lossUnit(l, maxSlice, result, tmp);
        result /= subdivision; 
        
        if (total /result > maxSlice){
            result = total / maxSlice;
        }
        return result;
    } catch (exception& e){
        throw ModelException(e, routine);
    }
}

    


/* Calcualtes bounds for the density arrays */ 
void CCMLossUnit::calcLossDensityBounds(
    const DoubleArray& amt,           /* (I) amounts as a double, 
                                         usually loss Cata */
    int&               maxLongLoss,   /* (O) sum(long loss amt) */
    int&               maxShortLoss){ /* (O) sum(short loss amount) */
    static char routine[] = "ccmCalcLossDensityBounds";
    int maxL = 0;
    int maxS = 0;
    
    // Use a "double" version of maxL and maxS in case computation
    // of maxL and maxS leads to numbers greater than the greatest integer (in which
    // case we have random behaviours potentially causing crashes or bad
    // allocation errors depending on the value of maxS+maxL+1)
    double maxLdb = 0.0;
    double maxSdb = 0.0;
    
    for (int i=0; i<amt.size(); ++i) {
        if (amt[i] > 0) {
            maxL += (long)(ceil(amt[i]));
            maxLdb += ceil(amt[i]);
        }
        if (amt[i] < 0) {
            maxS -= (long)(floor(amt[i]));
            maxSdb -= floor(amt[i]);
        }
        /* and if amt[i], do nothing */
    }    

    // check maxL and maxS won't lead to too many discretisation units for
    // the loss distribution
    double nbDiscretisationUnits = maxLdb + maxSdb + 1.0;
    if (nbDiscretisationUnits >= MAX_DISCRETISATION_UNITS)
    {
        throw ModelException(routine,
            "Unable to discretise loss distribution: "
            "too many discretisation units (" +
            Format::toString(nbDiscretisationUnits) +
            "). Check if the loss unit override is not too small.");
    }
    
    maxLongLoss  = maxL;
    maxShortLoss = maxS;
    if (maxL+maxS==0){
        throw ModelException(routine, "Sum of long and short loss is zero");
    }
}

#undef MAX_DISCRETISATION_UNITS

DRLIB_END_NAMESPACE

