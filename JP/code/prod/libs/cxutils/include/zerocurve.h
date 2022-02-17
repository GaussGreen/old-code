/*
***************************************************************************
** HEADER FILE: zerocurve.h
**
** Defines an extremely simple zero curve structure.
**
** $Header$
***************************************************************************
*/

#ifndef CX_ZEROCURVE_H
#define CX_ZEROCURVE_H

#include "cxutils.h"  /* basic data types */

#ifdef __cplusplus
extern "C"
{
#endif

/*f
***************************************************************************
** Constructor for TCurve for add-in users
***************************************************************************
*/
TCurve* CxZeroCurveMake(
TDate           baseDate,            /* (I) */
int             numItems,            /* (I) */
TDate*          zeroDate,            /* (I) [numItems] */
double*         zeroRate,            /* (I) [numItems] */
TRateType*      rateType,            /* (I) */
CxTDayCountConv  dcc                  /* (I) */
);

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
TCurve* CxZeroCurveAddLongRate (TCurve* in);

/*f
***************************************************************************
** Constructs a zero curve from a list of zero prices.
**
** The zero dates must be in ascending order.
** You can include today in the list of zero dates, but the
** corresponding zero price must be 1.0 in that case.
***************************************************************************
*/
TCurve* CxZeroCurveMakeFromZeroPrices
(TDate   today,       /* (I) */
 int     numItems,    /* (I) */
 TDate  *zeroDate,    /* (I) [numItems] */
 double *zeroPrice);  /* (I) [numItems] */

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
 TCurve* zc2);

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
 TDate   today);

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
 double*        annuity);     /* (O) Annuity */

/*f
***************************************************************************
** Calculates par forward yield for a specific maturity from a zero coupon
** curve (with currency basis).
**
** Makes a number of simplifying assumptions (no holidays etc).
**
** Note that we do not need to know the day count convention of the
** floating leg, but we do need to know its frequency (which must be
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
 double*        annuity);      /* (O) Annuity */

/*f
***************************************************************************
** Calculates a simple forward rate between two dates from a zero coupon
** curve.
***************************************************************************
*/
double CxForwardRate
(TCurve*        zeroCurve,    /* (I) Zero curve */
 TDate          startDate,    /* (I) Start date */
 TDate          maturityDate, /* (I) Maturity date */
 CxTDayCountConv dcc,          /* (I) Day count convention */
 long           frequency);   /* (I) Frequency of rate - 0 = simple rate */

/*f
***************************************************************************
** Calculates the zero price for a given date. Returns NaN for errors.
***************************************************************************
*/
double CxZeroPrice
(TCurve* zeroCurve,
 TDate   date);

/*f
***************************************************************************
** Calculates the zero price for a given start date and maturity date.
** Returns NaN for errors.
***************************************************************************
*/
double CxForwardZeroPrice
(TCurve* zeroCurve,
 TDate   startDate,
 TDate   maturityDate);

/*f
***************************************************************************
** Calculates the zero rate for a given date using ACT/365F and continously
** compounded rates. If you want another convention, then use CxForwardRate
** instead.
***************************************************************************
*/
double CxZeroRate
(TCurve* zeroCurve,
 TDate   date);

/*f
***************************************************************************
** Extracts all the times and zero prices from the zero curve.
**
** This is useful for analytics which just use time and discount factor
** and know nothing about dates.
***************************************************************************
*/
int CxZeroCurveTimesAndPrices
(TCurve*  zeroCurve,   /* (I) */
 int*     numItems,    /* (O) */
 double** times,       /* (O) */
 double** zeroPrices); /* (O) */

/*f
***************************************************************************
** Create a smoothed version of a zero curve.
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
 TDateList* newDates);      /* (I) */

/*f
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
 double *outRate);


#ifdef __cplusplus
}
#endif

#endif /* CX_ZERO_CURVE_D_H */

