/*
***************************************************************************
** SOURCE FILE: zerocurve.c
**
** An alternative interface to the ALIB TCurve interface - with some
** extra features.
**
** Based on the requirement that we need to combine CX and CRX.
** CRX uses TCurve, so rather than change CRX we decided to change CX from
** using CXZeroCurve so that it uses TCurve instead.
**
** $Header$
***************************************************************************
*/

#include "zerocurve.h"

#include <assert.h>
#include <ctype.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "bsearch.h"
#include "date.h"
#include "datelist.h"
#include "dateutils.h"
#include "cxmacros.h"

#include <alib/datelist.h>
#include <alib/ldate.h>
#include <alib/tcurve.h>
#include <alib/zr2fwd.h>

#define NaN sqrt(-1.0)

static int zcInterpRate (TCurve*, TDate, long, long, double*);
static int zcRateCC (TCurve*, int, double*);

/*f
***************************************************************************
** Constructor for TCurve
***************************************************************************
*/
TCurve* CxZeroCurveMake(
TDate           baseDate,            /* (I) */
int             numItems,            /* (I) */
TDate*          zeroDate,            /* (I) [numItems] */
double*         zeroRate,            /* (I) [numItems] */
TRateType*      rateType,            /* (I) */
CxTDayCountConv  dcc                  /* (I) */
)
{
    static char routine[] = "CxZeroCurveMake";
    
    TCurve *out = NULL;
    TRateType annualBasis = GTO_ANNUAL_BASIS;

    if (rateType == NULL)
        rateType = &annualBasis;

    out = GtoMakeTCurve (baseDate,
                         zeroDate,
                         zeroRate,
                         numItems,
                         *rateType,
                         (long)dcc);

    if (out == NULL)
        GtoErrMsgFailure (routine);

    return out;
}

/*f
***************************************************************************
** Adds a long rate to a TCurve so that ALIB interpolation will work
** beyond the last date of the curve.
**
** The resulting curve will be ANNUAL ACT/365F with a point 100 years
** beyond the last point in the curve such that the forward rate for the
** last 100 years matches the last forward rate in the input curve.
***************************************************************************
*/
TCurve* CxZeroCurveAddLongRate (TCurve* in)
{
    static char routine[] = "CxZeroCurveAddLongRate";
    int         status    = SUCCESS;

    TCurve     *tc = NULL;
    int         i;

    REQUIRE (in != NULL);

    tc = GtoNewTCurve (in->fBaseDate,
                       in->fNumItems+1,
                       GTO_ANNUAL_BASIS,
                       GTO_ACT_365F);

    if (tc == NULL)
        goto done;

    COPY_ARRAY (tc->fArray, in->fArray, TRatePt, in->fNumItems);
    if (IS_NOT_EQUAL(in->fBasis, GTO_ANNUAL_BASIS) || 
        in->fDayCountConv != GTO_ACT_365F)
    {
        /* rate conversions are necessary */
        /* do not use GtoRateToRate since that fails for rates at the
           baseDate and also does too much maths */
        for (i = 0; i < in->fNumItems; ++i)
        {
            if (CxConvertCompoundRate (in->fArray[i].fRate,
                                       in->fBasis,
                                       in->fDayCountConv,
                                       GTO_ANNUAL_BASIS,
                                       GTO_ACT_365F,
                                       &tc->fArray[i].fRate) != SUCCESS)
                goto done; /* failure */
        }
    }

    /* add the extra point 100 years past the last date for extrapolation */
    i = in->fNumItems;
    tc->fArray[i].fDate = tc->fArray[i-1].fDate + 36524;
    if (CxConvertCompoundRate (CxZeroRate (in, tc->fArray[i].fDate),
                               GTO_CONTINUOUS_BASIS,
                               GTO_ACT_365F,
                               GTO_ANNUAL_BASIS,
                               GTO_ACT_365F,
                               &tc->fArray[i].fRate) != SUCCESS)
        goto done; /* failure */

    status = SUCCESS;

  done:

    if (status != SUCCESS)
    {
        GtoFreeTCurve (tc);
        tc = NULL;
        GtoErrMsgFailure (routine);
    }

    return tc;
}
    


/*f
***************************************************************************
** Constructs a zero curve from a list of zero prices.
**
** The zero dates must be in ascending order.
** You can include the value date in the list of zero dates, but the
** corresponding zero price must be 1.0 in that case.
***************************************************************************
*/
TCurve* CxZeroCurveMakeFromZeroPrices
(TDate   valueDate,   /* (I) */
 int     numItems,    /* (I) */
 TDate  *zeroDate,    /* (I) [numItems] */
 double *zeroPrice)   /* (I) [numItems] */
{
    static char routine[] = "CxZeroCurveMakeFromZeroPrices";

    TCurve *zc = NULL; /* output curve */

    char    dateBuf[16];
    int     i;
    int     j = -1;
    double *zeroRate = NEW_ARRAY(double, numItems);

    for (i = 0; i < numItems; ++i)
    {
        if (zeroDate[i] != valueDate)
        {
            if (GtoDiscountToRate (zeroPrice[i],
                                   valueDate,
                                   zeroDate[i],
                                   GTO_ACT_365F,
                                   GTO_ANNUAL_BASIS,
                                   &zeroRate[i]) != SUCCESS)
                goto done; /* failure */
        }
        else
        {
            if (IS_NOT_EQUAL(zeroPrice[i], 1.0))
            {
                GtoErrMsg ("%s: Zero price (%f) on value date (%s) is not 1.0\n",
                         routine, zeroPrice[i],
                         CxDateFormat(valueDate, "DD-MMM-YYYY", dateBuf));
                goto done; /* failure */
            }
            j = i;
        }
    }

    if (j >= 0)
    {
        /* value date included in the curve - use the next date for this
           zero rate - the zero rate is somewhat meaningless for this date,
           but the existence of the date in the curve may make a difference
           if we have dates before the value date - also if we extrapolate
           from before the start of the curve, then the rule is that we use
           the first rate in the curve, so this is why we use the next
           rate */
        i = j+1;
        if ((j+1) < numItems)
            zeroRate[j] = zeroRate[j+1];
        else if (j > 0)
            zeroRate[j] = zeroRate[j-1];
        else
        {
            char d1[16];
            assert (numItems == 1);
            GtoErrMsg ("%s: Only one date provided - and it is the value date "
                     "%s\n", routine,
                     CxDateFormat(valueDate, "DD-MMM-YYYY", d1));
            goto done; /* failure */
        }
    }

    /* CxZeroCurveMake ensures that the dates are sorted */
    zc = GtoMakeTCurve (valueDate, zeroDate, zeroRate, numItems,
                        GTO_ANNUAL_BASIS, GTO_ACT_365F);
    if (zc == NULL) goto done; /* failure */

 done:

    if (zc == NULL) GtoErrMsgFailure(routine);

    FREE(zeroRate);
    return zc;
}


/*f
***************************************************************************
** The risky zero curve is the product of two zero curves. The first curve
** is typically the risk free curve. The second curve is the pure default
** curve.
**
** The resulting curve will have the minimum of the two base dates.
***************************************************************************
*/
TCurve* CxZeroCurveMakeRisky
(TCurve* zc1,
 TCurve* zc2)
{
    static char routine[] = "CxZeroCurveMakeRisky";

    TCurve      *zcRisky = NULL;
    int          i;
    TDateList   *dlzc1 = NULL;
    TDateList   *dlzc2 = NULL;
    TDateList   *dl = NULL;
    double      *zeroPrices = NULL;
    TDate        baseDate;

    REQUIRE (zc1 != NULL);
    REQUIRE (zc2 != NULL);

    baseDate = MIN(zc1->fBaseDate, zc2->fBaseDate);

    dlzc1 = GtoNewDateListFromTCurve (zc1);
    dlzc2 = GtoNewDateListFromTCurve (zc2);
    dl    = GtoMergeDateLists (dlzc1, dlzc2);
    if (zc1->fBaseDate > baseDate)
        dl = CxDateListAddDatesFreeOld (dl, 1, (TDate*)&zc1->fBaseDate);
    if (zc2->fBaseDate > baseDate)
        dl = CxDateListAddDatesFreeOld (dl, 1, (TDate*)&zc2->fBaseDate);

    zeroPrices = NEW_ARRAY(double, dl->fNumItems);

    for (i = 0; i < dl->fNumItems; ++i)
    {
        TDate date = dl->fArray[i];
		zeroPrices[i] = CxForwardZeroPrice(zc1, baseDate, date) * 
            CxForwardZeroPrice(zc2, baseDate, date);
	}

    zcRisky = CxZeroCurveMakeFromZeroPrices (baseDate,
                                             dl->fNumItems,
                                             dl->fArray,
                                             zeroPrices);

 done:
	
    if (zcRisky == NULL) GtoErrMsgFailure(routine);
    FREE (zeroPrices);
    GtoFreeDateList (dl);
    GtoFreeDateList (dlzc1);
    GtoFreeDateList (dlzc2);
    
    return zcRisky;
}

/*f
***************************************************************************
** Takes a copy of the zero curve with a new value of today.
**
** This is needed when we have a zero curve from another system which
** might not use today as the base date of the curve.
**
** The rules are as follows:
**
** 1. If today < old value of today, then the old value of today is
**    retained in the resulting curve.
** 2. If today = old value of today, then the resulting curve is a copy
**    of the input.
** 3. If today > old value of today, then we consider this to be an error.
***************************************************************************
*/
TCurve* CxZeroCurveMoveToToday
(TCurve *zc,
 TDate   today)
{
    static char routine[] = "CxZeroCurveMoveToToday";

    TCurve* out = NULL;

    char         buf1[20];
    char         buf2[20];
    TDateList*   dl = NULL;
    double*      zeroPrices = NULL;

    REQUIRE (zc != NULL);

    if (today > zc->fBaseDate)
    {
        GtoErrMsg ("%s: Cannot move to today (%s) greater than original value "
                 "of baseDate (%s)\n", routine,
                 CxDateFormat (today, "DD-MMM-YYYY", buf1),
                 CxDateFormat (zc->fBaseDate, "DD-MMM-YYYY", buf2));
        goto done; /* failure */
    }
    
    if (today == zc->fBaseDate)
    {
        out = GtoCopyCurve (zc);
    }
    else
    {
        double      p0;
        int         i;

        dl = GtoNewDateListFromTCurve (zc);
        if (dl == NULL) goto done; /* failure */
        
        dl  = CxDateListAddDatesFreeOld (dl, 1, &zc->fBaseDate);
        if (dl == NULL) goto done; /* failure */

        p0 = CxZeroPrice (zc, today); /* expect p0 > 1 */
        zeroPrices = NEW_ARRAY(double, dl->fNumItems);
        if (zeroPrices == NULL) goto done; /* failure */

        for (i = 0; i < dl->fNumItems; ++i)
            zeroPrices[i] = CxZeroPrice (zc, dl->fArray[i]) / p0;

        out = CxZeroCurveMakeFromZeroPrices (today,
                                             dl->fNumItems,
                                             dl->fArray,
                                             zeroPrices);
    }

 done:

    if (out == NULL) GtoErrMsgFailure(routine);

    FREE (zeroPrices);
    GtoFreeDateList (dl);

    return out;
}
        

/*f
***************************************************************************
** Calculates par forward yield for a specific maturity from a zero coupon
** curve (without currency basis).
**
** Makes a number of simplifying assumptions (no holidays etc).
**
** Returns SUCCESS or FAILURE.
**
** You can use NULL for parYield or annuity if you do not care about that
** particular value.
***************************************************************************
*/
int CxParYield
(TCurve*        zeroCurve,    /* (I) Zero curve */
 TDate          startDate,    /* (I) Start date */
 CxTDayCountConv dcc,          /* (I) Day count convention */
 long           numPeriods,   /* (I) Number of coupon periods */
 long           frequency,    /* (I) Frequency of coupon payments */
 double*        parYield,     /* (O) Par yield */
 double*        annuity)      /* (O) Annuity */
{
    static char routine[] = "CxParYield";
    int         status    = FAILURE;

    double      myParYield;
    double      myAnnuity;
    long        numMonths;
    long        i;
    TDate       prevDate;
    double      zeroPrice = 1.0;
    double      startZeroPrice;
    double      endZeroPrice;

    REQUIRE (numPeriods > 0);
    REQUIRE (frequency > 0);
    REQUIRE (12 % frequency == 0);
    
    numMonths = 12 / frequency;
    
    myAnnuity = 0.0;

    prevDate = startDate;
    for (i = 1; i <= numPeriods; ++i)
    {
        /* get date i periods ahead of start date (end of month adjusted) */
        TDate date = CxDateAddMonths (startDate, i * numMonths, TRUE);
        double time = CxDayCountFraction (prevDate, date, dcc);
        
        zeroPrice = CxZeroPrice (zeroCurve, date);
        myAnnuity += zeroPrice * time;

        prevDate = date;
    }

    startZeroPrice = CxZeroPrice (zeroCurve, startDate);
    endZeroPrice   = zeroPrice; /* last value calculated in the coupon loop */
    myParYield     = (startZeroPrice - endZeroPrice) / myAnnuity;

    if (parYield != NULL) *parYield = myParYield;
    if (annuity != NULL)  *annuity  = myAnnuity;
    status  = SUCCESS;

 done:
    
    if (status != SUCCESS) GtoErrMsgFailure(routine);
    return status;
}
 
/*f
***************************************************************************
** Calculates par forward yield for a specific maturity from a zero coupon
** curve (with currency basis).
**
** Makes a number of simplifying assumptions (no holidays etc).
**
** Note that we do not need to know the day count convention of the
** floating leg, but we do not need to know its frequency (which must be
** valid vis-a-vis the numPeriod and frequency parameters).
**
** Returns SUCCESS or FAILURE.
**
** You can use NULL for parYield or annuity if you do not care about that
** particular value.
***************************************************************************
*/
int CxParYieldCB
(TCurve*        indexCurve,    /* (I) Index curve */
 TCurve*        discountCurve, /* (I) Discount curve */
 TDate          startDate,     /* (I) Start date */
 CxTDayCountConv dcc,           /* (I) Day count convention */
 long           numPeriods,    /* (I) Number of coupon periods */
 long           frequency,     /* (I) Frequency of coupon payments */
 long           floatFreq,     /* (I) Frequency of floating payments */
 double*        parYield,      /* (O) Par yield */
 double*        annuity)       /* (O) Annuity */
{
    indexCurve=indexCurve;
    discountCurve=discountCurve;
    startDate=startDate;
    dcc=dcc;
    numPeriods=numPeriods;
    frequency=frequency;
    floatFreq=floatFreq;
    parYield=parYield;
    annuity=annuity;
    return FAILURE;
}

/*f
***************************************************************************
** Calculates a simple forward rate between two dates from a zero coupon
** curve. Returns NaN for errors.
***************************************************************************
*/
double CxForwardRate
(TCurve*        zeroCurve,    /* (I) Zero curve */
 TDate          startDate,    /* (I) Start date */
 TDate          maturityDate, /* (I) Maturity date */
 CxTDayCountConv dcc,          /* (I) Day count convention */
 long           frequency)    /* (I) Frequency of rate - 0 = simple rate */
{
    static char routine[] = "CxForwardRate";
    int         status    = FAILURE;

    double      rate;

    if (GtoForwardFromZCurve (zeroCurve,
                              GTO_FLAT_FORWARDS,
                              startDate,
                              maturityDate,
                              (long)dcc,
                              frequency,
                              &rate) != SUCCESS)
    {
        goto done; /* failure */
    }

    status = SUCCESS;
    
 done:

    if (status != SUCCESS)
    {
        GtoErrMsgFailure(routine);
        return NaN;
    }
    return rate;
}

/*f
***************************************************************************
** Extracts the times and zero prices from the zero curve.
***************************************************************************
*/
int CxZeroCurveTimesAndPrices
(TCurve*    zeroCurve,
 int*       numItems,
 double**   times,
 double**   zeroPrices)
{
    static char routine[] = "CxZeroCurveTimesAndPrices";
    int         status    = FAILURE;

    int     i;
    double *myTimes = NULL;
    double *myZeroPrices = NULL;

    REQUIRE(zeroCurve != NULL);
    REQUIRE(numItems != NULL);
    REQUIRE(times != NULL);
    REQUIRE(zeroPrices != NULL);

    myTimes = NEW_ARRAY(double, zeroCurve->fNumItems);
    myZeroPrices = NEW_ARRAY(double, zeroCurve->fNumItems);
    
    for (i = 0; i < zeroCurve->fNumItems; ++i)
    {
        myTimes[i] = (zeroCurve->fArray[i].fDate - zeroCurve->fBaseDate)/365.0;
        if (GtoRateToDiscount (zeroCurve->fArray[i].fRate,
                               zeroCurve->fBaseDate,
                               zeroCurve->fArray[i].fDate,
                               (long)zeroCurve->fBasis,
                               zeroCurve->fDayCountConv,
                               &myZeroPrices[i]) != SUCCESS)
            goto done; /* failure */
    }

    *numItems    = zeroCurve->fNumItems;
    *times       = myTimes;
    *zeroPrices  = myZeroPrices;
    myTimes      = NULL;
    myZeroPrices = NULL;
    status       = SUCCESS;

 done:

    if (status != SUCCESS) GtoErrMsgFailure(routine);
    
    FREE(myTimes);
    FREE(myZeroPrices);

    return status;
}


/*f
***************************************************************************
** Calculates the zero price for a given start date and maturity date.
** Returns NaN for errors.
***************************************************************************
*/
double CxForwardZeroPrice
(TCurve* zeroCurve,
 TDate   startDate,
 TDate   maturityDate)
{
    double startPrice    = CxZeroPrice(zeroCurve, startDate);
    double maturityPrice = CxZeroPrice(zeroCurve, maturityDate);
    return maturityPrice / startPrice;
}

/*f
***************************************************************************
** Calculates the zero price for a given date. Returns NaN for errors.
***************************************************************************
*/
double CxZeroPrice
(TCurve* zeroCurve,
 TDate   date)
{
    double      zeroPrice = 0.0;
    double      rate;
    double      time;

    rate = CxZeroRate (zeroCurve, date);

    /* 
    ** rate is continuously compounded calculated between valueDate of
    ** the zeroCurve and the required date
    */
    time = (date - zeroCurve->fBaseDate) / 365.0;
    zeroPrice = exp(-rate * time);
    return zeroPrice;
}

/*f
***************************************************************************
** Calculates the zero rate for a given date using ACT/365F and continously
** compounded rates.
**
** We do not trust the ALIB to do this because it does not extrapolate
** beyond the last date of the curve.
***************************************************************************
*/
double CxZeroRate
(TCurve* zeroCurve,
 TDate   date)
{
    static char routine[] = "CxZeroRate";
    int         status    = FAILURE;

    long        exact;
    long        lo;
    long        hi;
    double      rate = 0.0;

    REQUIRE (zeroCurve != NULL);
    REQUIRE (zeroCurve->fNumItems > 0);
    REQUIRE (zeroCurve->fArray != NULL);

    if (CxBinarySearchLong (date,
                            &zeroCurve->fArray[0].fDate,
                            sizeof(TRatePt),
                            zeroCurve->fNumItems,
                            &exact,
                            &lo,
                            &hi) != SUCCESS) 
        goto done; /* failure */

    if (exact >= 0)
    {
        /* date found in zeroDates */
        if (zcRateCC (zeroCurve, exact, &rate) != SUCCESS)
            goto done; /* failure */
    }
    else if (lo < 0)
    {
        /* date before start of zeroDates */
        if (zcRateCC (zeroCurve, 0, &rate) != SUCCESS)
            goto done; /* failure */
    }
    else if (hi >= zeroCurve->fNumItems)
    {
        /* date after end of zeroDates */
        if (zeroCurve->fNumItems == 1)
        {
            if (zcRateCC (zeroCurve, 0, &rate) != SUCCESS)
                goto done; /* failure */
        }
        else
        {
            /* extrapolate using last flat segment of the curve */
            lo = zeroCurve->fNumItems-2;
            hi = zeroCurve->fNumItems-1;
            if (zcInterpRate (zeroCurve, date, lo, hi, &rate) != SUCCESS)
                goto done; /* failure */
        }
    }
    else
    {
        /* date between start and end of zeroDates */
        if (zcInterpRate (zeroCurve, date, lo, hi, &rate) != SUCCESS)
            goto done; /* failure */
    }

    status    = SUCCESS;

 done:

    if (status != SUCCESS)
    {
        GtoErrMsgFailure (routine);
        return NaN;
    }
    return rate;
}

/*
***************************************************************************
** Interpolates a rate segment of a zero curve expressed with continuously
** compounded rates using flat forwards.
**
** Always returns a continously compounded ACT/365F rate.
***************************************************************************
*/
static int zcInterpRate 
(TCurve* zc, TDate date, long lo, long hi, double *rate)
{
    static char routine[] = "zcInterpRate";
    int         status    = FAILURE;

    long   t1;
    long   t2;
    long   t;
    double zt;
    double z1t1;
    double z2t2;
    double z1;
    double z2;

    t1   = zc->fArray[lo].fDate - zc->fBaseDate;
    t2   = zc->fArray[hi].fDate - zc->fBaseDate;
    t    = date - zc->fBaseDate;

    assert (t > t1);
    assert (t2 > t1);

    if (zcRateCC (zc, lo, &z1) != SUCCESS)
        goto done; /* failure */

    if (zcRateCC (zc, hi, &z2) != SUCCESS)
        goto done; /* failure */

    /* rates are continuously compounded, i.e. exp(-rt) */
    /* flat forwards => (zt) should be linear in t */
    z1t1 = z1 * t1;
    z2t2 = z2 * t2;
    if (t == 0)
    {
        /* If the date equals the base date, then the zero rate
           is essentially undefined and irrelevant - so let us
           get the rate for the following day which is in thr
           right ballpark at any rate. */
        /* An exception to this approach is when t2 = 0 as well.
           In this case, rate = z2 */
        if (t2 == 0)
        {
            *rate = z2;
            goto success;
        }
        t = 1;
    }
    zt   = z1t1 + (z2t2 - z1t1) * (double)(t - t1) / (double)(t2 - t1);
    *rate = zt / t;

success:
    status = SUCCESS;

 done:

    if (status != SUCCESS)
        GtoErrMsgFailure (routine);

    return status;
}



/*f
***************************************************************************
** Create a smoother version of a zero curve.
**
** All the dates of the existing curve are preserved.
**
** We need constraints for the beginning and the end. The end constraint
** will be dr/dt = 0 at the end date. The start constraint is either than
** r is linear on the first segment, but if this leads to a negative rate
** at the start date, then we switch to r = 0 at the start date.
**
** The user provides a timeline of extra dates that is required.
**
** Continuous forwards are used with maturity weighting (as opposed to
** duration weighting which is not available)
***************************************************************************
*/
TCurve* CxZeroCurveMakeSmoother
(TCurve*    zeroCurve,      /* (I) */
 TDateList* newDates)       /* (I) */
{
    static char routine[] = "CxZeroCurveMakeSmoother";
    int        status    = FAILURE;

    TCurve      *zcSmooth = NULL;
    double      *a = NULL; /* parabola is att + bt + c */
    double      *b = NULL;
    double      *c = NULL;
    double      *r = NULL; /* r is flat forward rate on each segment */
    double      *t = NULL; /* t is time for each segment */
    double      *z = NULL; /* z is discount at each zero date */
    int          i;

    double       rs; /* rate at start of segment */
    double       re; /* rate at end of segment */

    TDate       *zcDates = NULL; /* for original curve */
    TDateList   *allDates = NULL;
    double      *allZeros = NULL; /* discounts at all dates */

    REQUIRE (zeroCurve != NULL);
    REQUIRE (zeroCurve->fNumItems > 0);
    REQUIRE (newDates != NULL);
    REQUIRE (newDates->fNumItems > 0);

    a = NEW_ARRAY(double, zeroCurve->fNumItems);
    b = NEW_ARRAY(double, zeroCurve->fNumItems);
    c = NEW_ARRAY(double, zeroCurve->fNumItems);
    r = NEW_ARRAY(double, zeroCurve->fNumItems);
    t = NEW_ARRAY(double, zeroCurve->fNumItems);
    z = NEW_ARRAY(double, zeroCurve->fNumItems);

    zcDates = GtoDatesFromCurve (zeroCurve);
    if (zcDates == NULL)
        goto done; /* failure */

    for (i = 0; i < zeroCurve->fNumItems; ++i)
    {
        double tz;
        double rate;

        if (zcRateCC (zeroCurve, i, &rate) != SUCCESS)
            goto done; /* failure */
        tz = (double)(zcDates[i] - zeroCurve->fBaseDate)/365.0;
        z[i] = exp(-rate * tz);
        if (i > 0)
        {
            t[i] = (double)(zcDates[i] - zcDates[i-1]) / 365.0;
            r[i] = log(z[i-1]/z[i]) / t[i];
        }
        else
        {
            t[i] = tz;
            r[i] = rate;
        }
    }

    if (zeroCurve->fNumItems == 1)
    {
        /* we really have no choice but to be flat throughout */
        a[0] = 0.0;
        b[0] = 0.0;
        c[0] = r[0];
    }
    else
    {
        /* first segment */
        re = (r[1] * t[0] + r[0] * t[1]) / (t[0] + t[1]);
        /* constraint is that 3att + 2bt + c = re at t = t[0] */
        /* constraint is that att + bt + c = r at t = t[0] */
        /* these give us
           att = c + re - 2r
           bt  = 3r - re - 2c
        */
        /* linear in first segment */
        /* constraint is that a = 0 */
        /* we also require c >= 0 */
        a[0] = 0.0;
        b[0] = (re - r[0]) / t[0];
        c[0] = 2*r[0] - re;
        if (c[0] < 0.0)
        {
            c[0] = 0.0;
            b[0] = (3.0 * r[0] - re) / t[0];
            a[0] = (re - 2*r[0]) / (t[0] * t[0]);
        }
        /* middle segments */
        for (i = 1; i < zeroCurve->fNumItems-1; ++i)
        {
            rs = re;
            re = (r[i+1] * t[i] + r[i] * t[i+1]) / (t[i] + t[i+1]);
            /* constraints are:
               3att + 2bt + c = rs at t = 0
               3att + 2bt + c = re at t = t[i]
               att + bt + c = r at t = t[i]
            */
            c[i] = rs;
            b[i] = (3.0 * r[i] - re - 2.0 * rs) / t[i];
            a[i] = (re + rs - 2*r[i]) / (t[i] * t[i]);
        }
        /* last segment */
        /* constraints are:
           3att + 2bt + c = rs at t = 0
           att + bt + c = r at t = t[i]
           3at + b = 0 at t = t[i]
        */
        /* the last constraint is that dr/dt should be zero at the end */
        assert (i == zeroCurve->fNumItems-1);
        rs = re;
        c[i] = rs;
        b[i] = 1.5 * (r[i] - rs) / t[i];
        a[i] = 0.5 * (rs - r[i]) / (t[i] * t[i]);
    }

    allDates = CxDateListCopyUnique ((TDateList*)newDates); /*TBD - const */
    allDates = CxDateListAddDatesFreeOld (allDates, 
                                          zeroCurve->fNumItems,
                                          zcDates);
    allDates = CxDateListTruncate (allDates,
                                   zeroCurve->fBaseDate,
                                   FALSE, /* exclude base date */
                                   TRUE,  /* exclude before base date */
                                   TRUE); /* modify in place */
    allDates = CxDateListTruncate (allDates,
                                   zcDates[zeroCurve->fNumItems-1],
                                   TRUE,  /* include last date */
                                   FALSE, /* exclude after base date */
                                   TRUE); /* modify in place */
    if (allDates == NULL)
        goto done; /* failure */
    
    allZeros = NEW_ARRAY(double, allDates->fNumItems);
    if (allZeros == NULL)
        goto done; /* failure */

    for (i = 0; i < allDates->fNumItems; ++i)
    {
        long exact;
        long lo;
        long hi;
        TDate date = allDates->fArray[i];

        if (CxBinarySearchLong (date,
                                zcDates,
                                sizeof(TDate),
                                zeroCurve->fNumItems,
                                &exact,
                                &lo,
                                &hi) != SUCCESS) 
            goto done; /* failure */

        if (exact >= 0)
        {
            /* date found in zeroDates */
            allZeros[i] = z[exact]; /* we worked it out earlier */
        }
        else if (lo < 0)
        {
            /* we are in first segment */
            double t;

            t = (double)(date - zeroCurve->fBaseDate) / 365.0;
            allZeros[i] = exp(-(((a[0]*t + b[0]) * t + c[0]) * t));
        }
        else
        {
            /* we are not in the first segment */
            double t;
            
            /* can assert this because we truncated allDates */
            assert (hi < zeroCurve->fNumItems);
            /* requirement that zero curve has distinct dates */
            assert (hi == lo+1); 

            t = (double)(date - zcDates[lo]) / 365.0;
            allZeros[i] = z[lo] * exp(-(((a[hi]*t + b[hi]) * t + c[hi]) * t));
        }
    }

    zcSmooth = CxZeroCurveMakeFromZeroPrices (zeroCurve->fBaseDate,
                                              allDates->fNumItems,
                                              allDates->fArray,
                                              allZeros);
    if (zcSmooth == NULL)
        goto done; /* failure */

    status = SUCCESS;

 done:

    if (status != SUCCESS)
    {
        GtoErrMsgFailure (routine);
        GtoFreeTCurve (zcSmooth);
        zcSmooth = NULL;
    }

    FREE (allZeros);
    GtoFreeDateList (allDates);
    FREE (zcDates);
    FREE (z);
    FREE (t);
    FREE (r);
    FREE (c);
    FREE (b);
    FREE (a);

    return zcSmooth;
}

/*
***************************************************************************
** Gets a continuously compounded rate from a TCurve for element idx.
***************************************************************************
*/
static int zcRateCC 
(TCurve *tc,
 int     idx,
 double *ccRate)
{
    int status;

    status =  CxConvertCompoundRate (tc->fArray[idx].fRate,
                                     tc->fBasis,
                                     tc->fDayCountConv,
                                     GTO_CONTINUOUS_BASIS,
                                     GTO_ACT_365F,
                                     ccRate);
    return status;
}

/*
***************************************************************************
** Converts a compound rate from one frequency to another.
** Can also convert between ACT-style day count conventions.
***************************************************************************
*/
int CxConvertCompoundRate 
(double  inRate,
 double  inBasis,
 long    inDayCountConv,
 double  outBasis,
 long    outDayCountConv,
 double *outRate)
{
    static char routine[] = "CxConvertCompoundRate";
    int         status    = FAILURE;

    double      ccRate;

    /* this routine is a hotspot and was taking too long for the case where we
       do nothing */

    if (IS_EQUAL(inBasis,outBasis))
    {
        if (inDayCountConv == outDayCountConv)
        {
            *outRate = inRate;
        }
        else if (inDayCountConv == GTO_ACT_365F && outDayCountConv == GTO_ACT_360)
        {
            *outRate = inRate * 360.0/365.0;
        }
        else if (inDayCountConv == GTO_ACT_360 && outDayCountConv == GTO_ACT_365F)
        {
            *outRate = inRate * 365.0/360.0;
        }
        else
        {
            GtoErrMsg ("%s: Can only convert between ACT/360 and ACT/365F day count "
                       "conventions\n", routine);
            goto done; /* failure */
        }
    }
    else
    {
        double dayFactor;

        if (inDayCountConv == outDayCountConv)
        {
            dayFactor = 1.0;
        }
        else if (inDayCountConv == GTO_ACT_365F && outDayCountConv == GTO_ACT_360)
        {
            dayFactor = 360.0/365.0;
        }
        else if (inDayCountConv == GTO_ACT_360 && outDayCountConv == GTO_ACT_365F)
        {
            dayFactor = 365.0/360.0;
        }
        else
        {
            GtoErrMsg ("%s: Can only convert between ACT/360 and ACT/365F day count "
                       "conventions\n", routine);
            goto done; /* failure */
        }

        /* convert inRate to ccRate, then convert to outRate */
        if (IS_EQUAL(inBasis, GTO_CONTINUOUS_BASIS))
        {
            ccRate = inRate * dayFactor;
        }
        else if (inBasis >= 1.0 && inBasis <= 365.0)
        {
            ccRate = dayFactor * inBasis * log (1.0 + inRate / inBasis);
        }
        else
        {
            GtoErrMsg ("%s: Input basis %f is not a compounding frequency\n",
                     routine, inBasis);
            goto done; /* failure */
        }

        if (IS_EQUAL(outBasis, GTO_CONTINUOUS_BASIS))
        {
            *outRate = ccRate;
        }
        else if (outBasis >= 1.0 && outBasis <= 365.0)
        {
            *outRate = outBasis * (exp (ccRate/outBasis) - 1.0);
        }
        else
        {
            GtoErrMsg ("%s: Output basis %f is not a compounding frequency\n",
                     routine, outBasis);
            goto done; /* failure */
        }
    }
    status = SUCCESS;

 done:

    if (status != SUCCESS)
        GtoErrMsgFailure (routine);

    return status;
}

        
