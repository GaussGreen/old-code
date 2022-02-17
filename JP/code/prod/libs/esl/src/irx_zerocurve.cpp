#include "irx/zerocurve.h"

#include <math.h>

extern "C" IrxTBootstrapMethod ZeroInterpTypeFlag;

#include "irx/rate.h"

#include <irx/bsearch.h>
#include <irx/error.h>
#include <irx/macros.h>
#include <irx/irxutilio.h>
#include <irx/convert.h>

static int validateZeroDates (
    IrxTDate        baseDate,
    int             numDates,
    const IrxTDate *zeroDates,
    int            *baseDateIdx);

static int computeFwdRates (
    IrxTZeroCurve *zeroCurve);


IrxTZeroCurve* irxZeroCurveMakeFromRates(
    IrxTDate          baseDate,
    int               numDates,
    const IrxTDate   *zeroDates,
    const double     *zeroRates,
    IrxTRateType      rateType,
    IrxTDayCountConv  dcc)
{
    static char routine[] = "irxZeroCurveMakeFromRates";
    IrxTZeroCurve* zc = irxZeroCurveMakeEmpty(numDates > 0 ? numDates : 1);
    if (zc == NULL || irxZeroCurveConstructFromRates(zc,
                                                     baseDate,
                                                     numDates,
                                                     zeroDates,
                                                     zeroRates,
                                                     rateType,
                                                     dcc) != SUCCESS)
    {
        irxZeroCurveFree(zc);
        zc = NULL;
        irxErrorFailure (routine);
    }

    return zc;

}


int irxZeroCurveConstructFromRates(
    IrxTZeroCurve*    zc,
    IrxTDate          baseDate,
    int               numDates,
    const IrxTDate   *zeroDates,
    const double     *zeroRates,
    IrxTRateType      rateType,
    IrxTDayCountConv  dcc)
{
    static char routine[] = "irxZeroCurveConstructFromRates";
    int         status    = FAILURE;

    status = irxZeroCurveConstructFromRates2(zc,
                                             baseDate,
                                             numDates,
                                             zeroDates,
                                             zeroRates,
                                             rateType,
                                             ZeroInterpTypeFlag,
                                             dcc);

    if (status != SUCCESS)
        irxErrorFailure (routine);

    return status;
}


IrxTZeroCurve* irxZeroCurveMakeFromRates2(
    IrxTDate             baseDate,
    int                  numDates,
    const IrxTDate*      zeroDates,
    const double*        zeroRates,
    IrxTRateType         rateType,
    IrxTBootstrapMethod  interpMethod,
    IrxTDayCountConv     dcc)
{
    static char routine[] = "irxZeroCurveMakeFromRates2";
    IrxTZeroCurve* zc = irxZeroCurveMakeEmpty(numDates > 0 ? numDates : 1);
    if (zc == NULL || irxZeroCurveConstructFromRates2(zc,
                                                      baseDate,
                                                      numDates,
                                                      zeroDates,
                                                      zeroRates,
                                                      rateType,
                                                      interpMethod,
                                                      dcc) != SUCCESS)
    {
        irxZeroCurveFree(zc);
        zc = NULL;
        irxErrorFailure (routine);
    }

    return zc;
}

int irxZeroCurveConstructFromRates2(
    IrxTZeroCurve*       zc,
    IrxTDate             baseDate,
    int                  numDates,
    const IrxTDate*      zeroDates,
    const double*        zeroRates,
    IrxTRateType         rateType,
    IrxTBootstrapMethod  interpMethod,
    IrxTDayCountConv     dcc)
{
    static char routine[] = "irxZeroCurveConstructFromRates2";
    int         status    = FAILURE;

    int baseDateIdx = 0;
    int i;
    int offset;
    int numItems;
    double baseRate;
    IrxTBool addBaseDate;
    IrxTBool addLastDate = FALSE;

    IrxTDate lastDate = 0;

    REQUIRE(baseDate > 0);
    REQUIRE(numDates >= 0);
    if (numDates > 0)
    {
        REQUIRE(zeroDates != NULL);
        REQUIRE(zeroRates != NULL);

        if (validateZeroDates(baseDate, numDates, zeroDates,
                              &baseDateIdx) != SUCCESS)
            goto RETURN; /* failure */

        numItems = numDates;
        addBaseDate = zeroDates[baseDateIdx] != baseDate;
        if (addBaseDate)
            ++numItems;
        else
            baseRate = zeroRates[baseDateIdx];

        /* We extend the curve to at least 100 years for
         * swaption vol bootstrapping
         */
        lastDate = irxDateAddMonths (baseDate, 1200L, TRUE);

        if (zeroDates[numDates-1] < lastDate)
        {
            ++numItems;
            addLastDate = TRUE;
        }
    }
    else
    {
        addBaseDate = TRUE;
        numItems = 1;
    }

    /* check to see whether we need to reconstruct with bigger numItems */
    if (zc->numItems != numItems)
    {
        irxZeroCurveDestroy(zc);
        irxZeroCurveConstructEmpty(zc, numItems);
    }

    zc->baseDate = baseDate;
    zc->baseRate = baseRate;
    zc->maxDate  = 0;

    offset = 0;
    for (i = 0; i < baseDateIdx; ++i)
    {
        zc->startDates[i]      = zeroDates[i];
        if (irxRateToDiscount (zeroRates[i],
                               baseDate,
                               zeroDates[i],
                               dcc,
                               rateType,
                               &zc->prices[i]) != SUCCESS)
            goto RETURN; /* failure */
    }

    if (addBaseDate)
    {
        offset = 1;
        zc->startDates[baseDateIdx]      = baseDate;
        zc->prices[baseDateIdx] = 1.0;
    }

    for (i = baseDateIdx; i < numDates; ++i)
    {
        zc->startDates[i+offset]      = zeroDates[i];
        if (irxRateToDiscount (zeroRates[i],
                               baseDate,
                               zeroDates[i],
                               dcc,
                               rateType,
                               &zc->prices[i+offset]) != SUCCESS)
            goto RETURN; /* failure */
    }

    if (addBaseDate)
    {
        if (baseDateIdx > 0)
            zc->baseRate = zeroRates[baseDateIdx - 1];
        else if (baseDateIdx < numDates)
            zc->baseRate = zeroRates[baseDateIdx];
    }

    if (addLastDate)
    {
        /* We extend the curve to at least 100 years using flat zero curve
         * for swaption vol bootstrapping (FIX3)
         */
        zc->startDates[numDates+offset] = lastDate;
        if (irxRateToDiscount (zeroRates[numDates-1],
                               baseDate,
                               lastDate,
                               dcc,
                               rateType,
                               &zc->prices[numDates+offset]) != SUCCESS)
            goto RETURN; /* failure */

    }

    zc->interpMethod = interpMethod;
    switch (zc->interpMethod)
    {
        case IRX_BOOTSTRAP_LINEAR_ZEROES:
            status = SUCCESS;
            break;

        case IRX_BOOTSTRAP_FLAT:
            status = computeFwdRates(zc);
            break;

        case IRX_BOOTSTRAP_SMOOTH:
            irxError ("%s: Interp method IRX_BOOTSTRAP_SMOOTH not yet supported.\n", routine);
            status = FAILURE;
            break;
   
    }

 RETURN:
    if (status != SUCCESS)
        irxErrorFailure (routine);

    return status;
}


IrxTZeroCurve* irxZeroCurveMakeFromPrices(
IrxTDate          baseDate,
int               numDates,
const IrxTDate         *zeroDates,
const double           *zeroPrices)
{
    static char routine[] = "irxZeroCurveMakeFromPrices";
    int         status    = FAILURE;

    IrxTZeroCurve *zc = NULL; /* to be returned */

    int baseDateIdx = 0;
    int i;
    int offset;
    int numItems;
    IrxTBool addBaseDate;
    IrxTBool addLastDate = FALSE;

    IrxTDate lastDate = 0;
    double   lastZeroRate;

    IrxTRateType IRX_ANNUAL_RATE = { IRX_RT_FREQUENCY, 1 };

    REQUIRE(baseDate > 0);
    REQUIRE(numDates >= 0);
    if (numDates > 0)
    {
        REQUIRE(zeroDates != NULL);
        REQUIRE(zeroPrices != NULL);

        if (validateZeroDates(baseDate, numDates, zeroDates,
                              &baseDateIdx) != SUCCESS)
            goto RETURN; /* failure */

        numItems = numDates;
        if (zeroDates[baseDateIdx] == baseDate)
        {
            addBaseDate = FALSE;
            if (IS_NOT_EQUAL(zeroPrices[baseDateIdx],1.0))
            {
                irxError ("%s: Discount factor %f for base date is not 1.0\n",
                          routine, zeroPrices[baseDateIdx]);
                goto RETURN; /* failure */
            }
        }
        else
        {
            addBaseDate = TRUE;
            ++numItems;
        }

        /* We extend the curve to at least 100 years for
         * swaption vol bootstrapping
         */
        lastDate = irxDateAddMonths (baseDate, 1200L, TRUE);

        if (zeroDates[numDates-1] < lastDate)
        {
            ++numItems;
            addLastDate = TRUE;
        }
    }
    else
    {
        addBaseDate = TRUE;
        numItems = 1;
    }

    zc = irxZeroCurveMakeEmpty(numItems);
    if (zc == NULL)
        goto RETURN; /* failure */

    zc->baseDate = baseDate;
    zc->maxDate  = 0;

    offset = 0;
    for (i = 0; i < baseDateIdx; ++i)
    {
        zc->startDates[i]      = zeroDates[i];
        zc->prices[i] = zeroPrices[i];
    }

    if (addBaseDate)
    {
        offset = 1;
        zc->startDates[baseDateIdx]      = baseDate;
        zc->prices[baseDateIdx] = 1.0;
    }

    for (i = baseDateIdx; i < numDates; ++i)
    {
        zc->startDates[i+offset]      = zeroDates[i];
        zc->prices[i+offset] = zeroPrices[i];
    }

    if (addLastDate)
    {
        /* We extend the curve to at least 100 years using flat zero curve
         * for swaption vol bootstrapping (FIX3)
         */
        zc->startDates[numDates+offset] = lastDate;
        if (irxDiscountToRate (zc->prices[numDates+offset-1],
                               baseDate,
                               zc->startDates[numDates+offset-1],
                               IRX_ACT_365F,
                               IRX_ANNUAL_RATE,
                               &lastZeroRate) != SUCCESS)
            goto RETURN; /* failure */

        if (irxRateToDiscount (lastZeroRate,
                               baseDate,
                               lastDate,
                               IRX_ACT_365F,
                               IRX_ANNUAL_RATE,
                               &zc->prices[numDates+offset]) != SUCCESS)
            goto RETURN; /* failure */
    }

    zc->interpMethod = ZeroInterpTypeFlag;
    switch (zc->interpMethod)
    {
        case IRX_BOOTSTRAP_LINEAR_ZEROES:
            status = SUCCESS;
            break;

        case IRX_BOOTSTRAP_FLAT:
          status = computeFwdRates(zc);
            break;

        case IRX_BOOTSTRAP_SMOOTH:
            irxError ("%s: Interp method IRX_BOOTSTRAP_SMOOTH not yet supported.\n", routine);
            status = FAILURE;
            break;
    }

 RETURN:
    if (status != SUCCESS)
    {
        irxZeroCurveFree(zc);
        zc = NULL;
        irxErrorFailure (routine);
    }

    return zc;
}


/**---------------------------------------------------------
 * I/O: read IrxTZeroCurve from wrapper file (skip comments).
     * (1) London format:
     * # Start date
     * 19970411
     * # Money Market basis (360 or 365)
     * 360
     * # Annual or semi-annual curve ("A" or "S")
     * S
     * # Year basis for benchmark swaps ("ACT", "365" or "360")
     * ACT
     * # No of entries
     * 84
     * #zero maturity yyyymmdd rates (ACT/365F annual)
     * 19970412 5.734349908886
     * 19970511 5.913551526763
     * etc
 */
int  irxZeroCurveConstructFromDRWFile(IrxTZeroCurve *zc, FILE *fp)
{
    static  char    routine[] = "irxZeroCurveConstructFromDRWFile";
    int status      = FAILURE;

    char            c;
    IrxTDate        baseDate;
    int             MMB;
    char            SwapFreq, SwapDCC[4], *SwapDCC1 = SwapDCC;

    int             nPts;
    IrxTDate        *zeroDate = NULL;    /* Zero maturity dates               */
    double          *zeroRate = NULL;    /* Zero rates (Annual ACT/365 basis) */

    /*
     * Identify format
     */
    if ((c = fgetc(fp)) == EOF) goto RETURN;
    ungetc(c, fp);


    READ_DATA(IRX_TDATE_T,   &baseDate,   "baseDate");
    READ_DATA(IRX_INT_T,     &MMB,        "MMB");
    READ_DATA(IRX_CHAR_T,    &SwapFreq,   "SwapFreq");
    READ_DATA(IRX_STRING_T,  &SwapDCC1,   "SwapDCC");

    READ_DATA(IRX_INT_T,     &nPts,       "numPoints");

    if (nPts>0)
    {
        if(irxVectArrayFpReadV(
                  fp,
                  nPts,
                  IRX_TDATE_T,   (void*) &zeroDate,
                  IRX_PERCENT_T, (void*) &zeroRate,
                  IRX_NULL_T) == FAILURE)
        {
            irxError("%s: Cannot read zero curve schedule.\n", routine);
            goto RETURN;
        }

    }

    if (irxZeroCurveConstructFromRates(
                        zc,
                        baseDate,
                        nPts,
                        zeroDate,
                        zeroRate,
                        IRX_ANNUAL_RATE,   /* Annual/365F */
                        IRX_ACT_365F) != SUCCESS)
        goto RETURN;


    /* Legacy Fix3 stuffs */
    zc->ValueDate = baseDate;
    zc->Today     = baseDate;
    zc->SpotDays  = 0;

    zc->MMB       = MMB;
    zc->SwapFreq  = SwapFreq;
    strcpy(zc->SwapDCC, SwapDCC);

    /* made it through */
    status = SUCCESS;


RETURN:
    if (status != SUCCESS) {
        irxError("%s: failed.\n", routine);
    }

    FREE(zeroDate);
    FREE(zeroRate);

    return(status);

}



/*
 * Validates the zero dates to ensure that they are in ascending order.
 * Also returns the insertion point for the baseDate.
 *
 * zeroDates[*baseDateIdx]   >= baseDate;
 * zeroDates[*baseDateIdx-1] < baseDate;
 *
 * Obviously if there is equality, then there is no need to insert the
 * baseDate.
 */
static int validateZeroDates (
    IrxTDate        baseDate,
    int             numDates,
    const IrxTDate *zeroDates,
    int            *baseDateIdx)
{
    static char routine[] = "validateZeroDates";

    int i;

    *baseDateIdx = 0;

    for (i = 0; i < numDates; ++i)
    {
        if (i > 1 && (zeroDates[i] <= zeroDates[i-1]))
        {
            char buf1[16];
            char buf2[16];
            irxError ("%s: zero dates %s and %s not in ascending order",
                      routine,
                      irxDateFormat(zeroDates[i-1],"DD-MMM-YYYY", buf1),
                      irxDateFormat(zeroDates[i],"DD-MMM-YYYY", buf2));
            return FAILURE;
        }
        if (zeroDates[i] < baseDate)
            *baseDateIdx = i+1;
    }
    return SUCCESS;
}

/*
 * Computes the forward rates for a zero curve given its discount factors.
 * Assumes we are using flat forwards and that a0,a1,a2 are not yet
 * populated.
 */
static int computeFwdRates (
    IrxTZeroCurve *zc)
{
    static char routine[] = "computeFwdRates";
    int         status    = FAILURE;

    int i;
    double rate;

    REQUIRE(zc != NULL);

    rate = 0.0;
    for (i = 0; i < zc->numItems-1; ++i)
    {
        double discount;
        discount = zc->prices[i+1] / zc->prices[i];
        if (irxDiscountToRate (discount,
                               zc->startDates[i],
                               zc->startDates[i+1],
                               IRX_ACT_365F,
                               IRX_CONTINUOUS_RATE,
                               &rate) != SUCCESS)
            goto RETURN; /* failure */
        zc->a0[i] = rate;
    }
    zc->a0[zc->numItems-1] = rate;

    status = SUCCESS;

 RETURN:

    if (status != SUCCESS)
        irxErrorFailure (routine);

    return status;
}


IrxTDate irxZeroCurveBaseDate(const IrxTZeroCurve *zeroCurve)
{
    return zeroCurve != NULL ? zeroCurve->baseDate : -1;
}

IrxTDate irxZeroCurveLastDate(const IrxTZeroCurve *zeroCurve)
{
    return zeroCurve != NULL ? zeroCurve->startDates[zeroCurve->numItems-1] : -1;
}


/**
 * Calculates the zero price for a given date. This is the price as seen
 * from the base date of the zero curve.
 */
int irxZeroPrice(
    const IrxTZeroCurve *zc,
    IrxTDate             date,
    double              *price)
{
    static char routine[] = "irxZeroPrice";
    int         status    = FAILURE;

    long   lo;
    long   hi;
    long   exact;
    long   i;

    IrxTDate startDate;
    double   a0;
    double   a1;
    double   a2;
    double   startPrice;

    REQUIRE(zc != NULL);

    if (irxBinarySearchLong (date,
                             zc->startDates,
                             sizeof(IrxTDate),
                             zc->numItems,
                             &exact,
                             &lo,
                             &hi) != SUCCESS)
        goto RETURN; /* failure */

    if (exact >= 0)
    {
        i = exact;
    }
    else if (lo < 0)
    {
        i = 0;
    }
    else if (hi >= zc->numItems)
    {
        /* date after end of zeroDates - depends on maxDate */
        if (date > zc->maxDate)
        {
            char buf1[16];
            char buf2[16];
            IrxTDate lastDate = zc->startDates[zc->numItems-1];
            if (zc->maxDate > lastDate)
                lastDate = zc->maxDate;

            irxError("%s: Cannot compute price for date %s beyond last "
                     "available date in curve %s.\n", routine,
                     irxDateFormat(date, "DD-MMM-YYYY", buf1),
                     irxDateFormat(lastDate, "DD-MMM-YYYY", buf2));
            goto RETURN; /* failure */
        }
        i = zc->numItems-1;
    }
    else
    {
        i = lo;
    }

    startDate  = zc->startDates[i];
    startPrice = zc->prices[i];

    if (date == startDate)
    {
        *price = startPrice;
    }
    else
    {
        switch (zc->interpMethod)
    	{
    		case IRX_BOOTSTRAP_FLAT:
    		case IRX_BOOTSTRAP_SMOOTH:
    		{
				double t;  /* time since startDate */
				double rt; /* integral of r.dt from 0 to t */
				a0 = zc->a0[i];
				t = (double)(date-startDate)/365.0;
				if (date < startDate)
				{
					ASSERT (lo < 0);
					/* before start of curve use the first spot rate available */
					/* this is a0 - and we ignore a1 and a2 */
					rt = a0 * t;
				}
				else
				{
					a1 = zc->a1[i];
					a2 = zc->a2[i];
					rt = ((a2 * t / 3.0 + a1 / 2.0) * t + a0) * t;
				}
				*price = startPrice * exp(-rt);
				break;
			}
			case IRX_BOOTSTRAP_LINEAR_ZEROES:
			{
				double rzlo;
				double rzhi;
				double rzmid = 0;

				if (lo < 0)
				{
					lo = 0;
				}

				if (hi >= zc->numItems)
				{
					hi = zc->numItems - 1;
					
					if (hi < 0)
					{
						irxError("%s : No curve components to perform interpolation", routine);
					}
				}

                /* handle case where all we have is DF = 1 at baseDate therefore, use
                 * the baseRate which we cleverly (hackily) stored during construction */
                if (zc->startDates[lo] == zc->baseDate)
                {
                    rzlo = zc->baseRate;
                }
                else
                {
				    if (irxDiscountToRate(zc->prices[lo],
					    				  zc->baseDate,
						    			  zc->startDates[lo],
							    		  IRX_ACT_365F,
								    	  IRX_ANNUAL_RATE,
									      &rzlo) == FAILURE) goto RETURN;
                }

				if (irxDiscountToRate(zc->prices[hi],
									  zc->baseDate,
									  zc->startDates[hi],
									  IRX_ACT_365F,
									  IRX_ANNUAL_RATE,
									  &rzhi) == FAILURE) goto RETURN;

				if (zc->startDates[lo] > date)
				{
					rzmid = rzlo;
				}
				else if (date > zc->startDates[hi])
				{
					rzmid = rzhi;
				}
				else
				{
					rzmid = rzlo + (rzhi - rzlo) * (date - zc->startDates[lo]) / (zc->startDates[hi] - zc->startDates[lo]);
				}

				if (irxRateToDiscount(rzmid,
									  zc->baseDate,
									  date,
									  IRX_ACT_365F,
									  IRX_ANNUAL_RATE,
									  price) == FAILURE) goto RETURN;
				break;
			}
			default:
			{
				irxError("%s : Unknown interpolation type", routine);
            	goto RETURN; /* failure */
            }
        }
    }

    status = SUCCESS;

 RETURN:

    if (status != SUCCESS)
        irxErrorFailure (routine);

    return status;
}


/**
 * Calculates the zero price for a given start date and maturity date.
 */
int irxFwdZeroPrice(
    const IrxTZeroCurve *zc,
    IrxTDate             startDate,
    IrxTDate             endDate,
    double              *fwdPrice)
{
    static char routine[] = "irxForwardZeroPrice";
    int         status    = FAILURE;

    double startPrice;
    double endPrice;

    if (irxZeroPrice (zc, startDate, &startPrice) != SUCCESS)
        goto RETURN; /* failure */

    if (irxZeroPrice (zc, endDate, &endPrice) != SUCCESS)
        goto RETURN; /* failure */

    *fwdPrice = endPrice/startPrice;
    status    = SUCCESS;

 RETURN:

    if (status != SUCCESS)
        irxErrorFailure (routine);

    return status;
}

/**
 * Calculates the spot zero rate from the base date to the given date
 */
int irxZeroRate(
    const IrxTZeroCurve *zc,
    IrxTDate             date,
    IrxTDayCountConv     dcc,
    IrxTRateType         rateType,
    double              *rate)
{
    static char routine[] = "irxZeroRate";
    int         status    = FAILURE;

    double discount;

    if (irxZeroPrice (zc, date, &discount) != SUCCESS)
        goto RETURN; /* failure */

    if (irxDiscountToRate (discount, zc->baseDate, date, dcc, rateType,
                           rate) != SUCCESS)
        goto RETURN; /* failure */

    status = SUCCESS;

 RETURN:

    if (status != SUCCESS)
        irxErrorFailure (routine);

    return status;
}



/**
 * Calculates a forward zero rate for a given start date and end date.
 */
int irxFwdZeroRate(
    const IrxTZeroCurve *zc,
    IrxTDate             startDate,
    IrxTDate             endDate,
    IrxTDayCountConv     dcc,
    IrxTRateType         rateType,
    double              *fwdRate)
{
    static char routine[] = "irxFwdZeroRate";
    int         status    = FAILURE;

    double discount;

    if (irxFwdZeroPrice (zc, startDate, endDate, &discount) != SUCCESS)
        goto RETURN; /* failure */

    if (irxDiscountToRate (discount, startDate, endDate, dcc, rateType,
                           fwdRate) != SUCCESS)
        goto RETURN; /* failure */

    status = SUCCESS;

 RETURN:

    if (status != SUCCESS)
        irxErrorFailure (routine);

    return status;
}

/**
 * Sets the zero curve's notion of Today and adjusts the zero curve to start
 * from Today (where it previously started from baseDate).
 * Returns the discount factor between today and baseDate
 */
double irxZeroCurveSetAndExtendToToday(IrxTZeroCurve* zc, IrxTDate newToday) {
    static char routine[] = "irxZeroCurveExtendToToday";
    int status = FAILURE;

    int             i;
    double          discount = 0.0;
    int             dgap0, dgap1;
    IrxTDate*       startDates = NULL;
    double*         prices = NULL;
    double*         a0 = NULL;
    double*         a1 = NULL;
    double*         a2 = NULL;

    REQUIRE(zc != NULL);
    REQUIRE(zc->numItems >= 2);
    REQUIRE(zc->startDates[0] == zc->baseDate);
    REQUIRE(newToday <= zc->Today);

    zc->Today = newToday;

    dgap0 = zc->baseDate - zc->Today;

    REQUIRE(dgap0 >= 0);

    if (dgap0 == 0) return 1.0;

    dgap1 = zc->startDates[1] -  zc->startDates[0];

    discount = pow(zc->prices[1], ((double)dgap0)/dgap1);

    startDates = zc->startDates;
    prices = zc->prices;
    a0 = zc->a0;
    a1 = zc->a1;
    a2 = zc->a2;

    zc->startDates = NULL;
    zc->prices = NULL;
    zc->a0 = NULL;
    zc->a1 = NULL;
    zc->a2 = NULL;

    ++(zc->numItems);

    zc->startDates = NEW_ARRAY(IrxTDate, zc->numItems);
    if (zc->startDates == NULL) goto RETURN;

    zc->prices = NEW_ARRAY(double, zc->numItems);
    if (zc->prices == NULL) goto RETURN;

    zc->a0 = NEW_ARRAY(double, zc->numItems);
    if (zc->a0 == NULL) goto RETURN;

    zc->a1 = NEW_ARRAY(double, zc->numItems);
    if (zc->a1 == NULL) goto RETURN;

    zc->a2 = NEW_ARRAY(double, zc->numItems);
    if (zc->a2 == NULL) goto RETURN;

    zc->startDates[0] = zc->Today;
    zc->prices[0] = 1.0;
    zc->a0[0] = -log(prices[1])/(dgap1/365.0);
    zc->a1[0] = 0.0;
    zc->a2[0] = 0.0;

    COPY_ARRAY(zc->startDates+1, startDates, IrxTDate, zc->numItems-1);
    COPY_ARRAY(zc->prices+1, prices, double, zc->numItems-1);
    COPY_ARRAY(zc->a0+1, a0, double, zc->numItems-1);
    COPY_ARRAY(zc->a1+1, a1, double, zc->numItems-1);
    COPY_ARRAY(zc->a2+1, a2, double, zc->numItems-1);

    for (i=1; i<zc->numItems; ++i)
        zc->prices[i] *= discount;

    zc->baseDate = zc->Today;
    status = SUCCESS;

RETURN:
    IRX_FREE(startDates);
    IRX_FREE(prices);
    IRX_FREE(a0);
    IRX_FREE(a1);
    IRX_FREE(a2);

    if (status != SUCCESS)
    {
        irxError("%s: Failed!\n", routine);
        IRX_FREE(zc->startDates);
        IRX_FREE(zc->prices);
        IRX_FREE(zc->a0);
        IRX_FREE(zc->a1);
        IRX_FREE(zc->a2);
        zc->numItems = 0;
    }


    return discount;
}


double irxZeroCurveExtendToToday(IrxTZeroCurve* zc) {
    return irxZeroCurveSetAndExtendToToday(zc, zc->Today);
}
