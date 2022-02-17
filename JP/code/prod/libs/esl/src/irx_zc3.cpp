/*
 * zero curve bootstrapping
 *
 * this file is far too big - consider some refactoring
 */

#include "irx/zc3.h"

#include <math.h>

#include "irx/cfl.h"
#include "irx/rate.h"
#include "irx/swap.h"
#include <irx/calendar.h>
#include <irx/dateutils.h>
#include <irx/macros.h>
#include <irx/rtbrent.h>

#define LOWER_BOUND      (-10000.00)
#define UPPER_BOUND      (10000.00)
#define GUESS_RATE       (0.06)
#define INITIAL_X_STEP   (0.0001)
#define INITIAL_F_DERIV  (0)
#define X_TOLERANCE      (DBL_MAX)
#define F_TOLERANCE      (1E-10)
#define MAX_ITERATIONS   (50)

/* Validates that the dates are in ascending order and that the instrument
   types make sense */
static int validateDates
(char*               routine,
 IrxTDate            today,
 IrxTDate            baseDate,
 int                 numInstruments,
 char**              instTypes,
 IrxTDate*           maturityDates,
 long*               includeFlags,
 int                 numShortRates,
 IrxTDate*           shortMaturityDates,
 long*               shortIncludeFlags);

/* Builds the short part of the curve - with or without currency basis */
static IrxTZeroCurve* buildShortEnd
(char*               routine,
 IrxTMarketConv     *marketConv,
 IrxTDate            today,
 IrxTDate            baseDate,
 int                 numInstruments,
 int                 numShortRates,
 IrxTDate*           shortMaturityDates,
 double*             shortRates,
 double*             shortCbsSpreads,
 long*               shortIncludeFlags);

/* Count a set of rates given by a number of include flags */
static int countRates
(int   numRates,
 long *includeFlags,
 long  maxInclude);



/* Adds a cash rate (with or without currency basis) */
static int addCashRate
(char*               routine,
 IrxTZeroCurve      *zc,
 int                 j,
 IrxTBool            smoothing,
 IrxTBool            first,
 IrxTBool            last,
 IrxTMarketConv     *marketConv,
 IrxTDate            baseDate,
 IrxTDate            cashDate,
 IrxTCalendar       *calendar,
 double              rate,
 double              cbsSpread);



/* Adds a swap rate (with or without currency basis) 
   For currency basis case this will be for calculating the indexCurve */
static int addSwapRate
(char*               routine,
 IrxTZeroCurve      *zc,
 IrxTZeroCurve const*zcDiscount,
 int                 j,
 IrxTBool            smoothing,
 IrxTBool            first,
 IrxTBool            last,
 IrxTMarketConv     *marketConv,
 IrxTDate            baseDate,
 IrxTDate            swapDate,
 IrxTCalendar       *calendar,
 IrxTBool            valueFloating,
 IrxTBool            firstFloatFixed,
 double              firstFloatFixedRate,
 IrxTStubLocation    stubLocation,
 double              rate,
 double              price);

/* Adds a swap rate (with currency basis) to the discount curve */
static int addSwapRateCB
(char*               routine,
 IrxTZeroCurve      *zcDiscount,
 int                 j,
 IrxTBool            smoothing,
 IrxTBool            first,
 IrxTBool            last,
 IrxTMarketConv     *marketConv,
 IrxTDate            baseDate,
 IrxTDate            swapDate,
 IrxTCalendar       *calendar,
 IrxTCalendar       *cbsCalendar,
 IrxTStubLocation    stubLocation,
 double              rate,
 double              cbsSpread,
 double              price);


/* Bootstraps a zero curve using flat forwards without currency basis */
static IrxTZeroCurve* buildFlatForwards
(char*               routine,
 IrxTMarketConv     *marketConv,
 IrxTDate            today,
 IrxTDate            baseDate,
 int                 numInstruments,
 char**              instTypes,
 IrxTDate*           maturityDates,
 double*             rates,
 double*             prices,
 long*               includeFlags,
 IrxTCalendar*       calendar,
 IrxTBool            valueFloating,
 int                 numShortRates,
 IrxTDate*           shortMaturityDates,
 double*             shortRates,
 long*               shortIncludeFlags,
 IrxTBool            firstFloatFixed,
 double              firstFloatFixedRate,
 IrxTStubLocation    stubLocation,
 IrxTDate            extrapDate);

/* Bootstraps a zero curve using smooth forwards without currency basis */
static IrxTZeroCurve* buildSmoothForwards
(char*               routine,
 IrxTZeroCurve      *zcFlat,
 IrxTMarketConv     *marketConv,
 IrxTDate            today,
 IrxTDate            baseDate,
 int                 numInstruments,
 char**              instTypes,
 IrxTDate*           maturityDates,
 double*             rates,
 double*             prices,
 long*               includeFlags,
 double*             adjustments,
 IrxTCalendar*       calendar,
 IrxTBool            valueFloating,
 int                 numShortRates,
 IrxTDate*           shortMaturityDates,
 double*             shortRates,
 long*               shortIncludeFlags,
 IrxTBool            firstFloatFixed,
 double              firstFloatFixedRate,
 IrxTStubLocation    stubLocation,
 IrxTDate            extrapDate);

/* Bootstraps a zero curve using flat forwards with currency basis */
static int buildFlatForwardsCB
(char*               routine,
 IrxTBool            reuseDiscountCurve,
 IrxTMarketConv     *marketConv,
 IrxTDate            today,
 IrxTDate            baseDate,
 int                 numInstruments,
 char**              instTypes,
 IrxTDate*           maturityDates,
 double*             rates,
 double*             cbsSpreads,
 double*             prices,
 long*               includeFlags,
 IrxTCalendar*       calendar,
 IrxTCalendar*       cbsCalendar,
 int                 numShortRates,
 IrxTDate*           shortMaturityDates,
 double*             shortRates,
 double*             shortCbsSpreads,
 long*               shortIncludeFlags,
 IrxTBool            firstFloatFixed,
 double              firstFloatFixedRate,
 IrxTStubLocation    stubLocation,
 IrxTDate            extrapDate,
 IrxTZeroCurve**     zcDiscount,
 IrxTZeroCurve**     zcIndex);

/* Bootstraps a zero curve using smooth forwards with currency basis */
static int buildSmoothForwardsCB
(char*               routine,
 IrxTZeroCurve      *zcDiscountFlat,
 IrxTZeroCurve      *zcIndexFlat,
 IrxTBool            reuseDiscountCurve,
 IrxTMarketConv     *marketConv,
 IrxTDate            today,
 IrxTDate            baseDate,
 int                 numInstruments,
 char**              instTypes,
 IrxTDate*           maturityDates,
 double*             rates,
 double*             cbsSpreads,
 double*             prices,
 long*               includeFlags,
 double*             adjustments,
 IrxTCalendar*       calendar,
 IrxTCalendar*       cbsCalendar,
 int                 numShortRates,
 IrxTDate*           shortMaturityDates,
 double*             shortRates,
 double*             shortCbsSpreads,
 long*               shortIncludeFlags,
 IrxTBool            firstFloatFixed,
 double              firstFloatFixedRate,
 IrxTStubLocation    stubLocation,
 IrxTDate            extrapDate,
 IrxTZeroCurve**     zcDiscount,
 IrxTZeroCurve**     zcIndex);

static int initSmoothCurve
(IrxTZeroCurve*       zcSmooth,
 IrxTZeroCurve const* zcFlat,
 int                  numInstruments,
 IrxTDate*            maturityDate,
 long*                includeFlags,
 int                  baseDateIndex);

/* Calculate durations corresponding to the maturity dates.
   The method chosen is somewhat artificial and agrees with the existing
   implementation. */
static int calcDuration (
    IrxTZeroCurve const *zc,
    IrxTDate             maturityDate,
    double              *duration);

/* For the first point of the curve, computes a1[j] and prices[j+1] given
   a0[j] */
static void smoothCalcFirst(double a0, double a0next, double time,
                            double *a2, double *discount);

/* For middle point of the curve, computes a2[j] and prices[j+1] given
   a1[j] */
static void smoothCalcMiddle(double a0, double a0next, double a1, double time,
                             double *a2, double *discount);

/* For the last point of the curve, computes a2[j] and prices[j+1] given
   a2[j] */
static void smoothCalcLast(double a0, double a1, double time,
                           double *a2, double *discount);

typedef struct
{
    int               j;
    IrxTBool          valueFloating;
    IrxTZeroCurve    *zc;
    /* it is possible that zcDisc = zc */
    IrxTZeroCurve const *zcDisc;
    IrxTCashFlowList *fixedFlows;
    IrxTSwapPayments *floatPayments;
    double            pv;
    double            time;
} SOLVER_CONTEXT;

static int flatObjFunc (double x, SOLVER_CONTEXT *context,
                        double *pvDiff);
static int smoothObjFuncFirstPoint (double x, SOLVER_CONTEXT *context,
                                    double *pvDiff);
static int smoothObjFuncMiddlePoint (double x, SOLVER_CONTEXT *context,
                                     double *pvDiff);
static int smoothObjFuncLastPoint (double x, SOLVER_CONTEXT *context,
                                   double *pvDiff);
static int calcPvDiff (SOLVER_CONTEXT *context, double *pvDiff);

static int flatSolver (double rate, SOLVER_CONTEXT *context);
static int smoothSolverFirstPoint (double rate, SOLVER_CONTEXT *context);
static int smoothSolverMiddlePoint (SOLVER_CONTEXT *context);
static int smoothSolverLastPoint (SOLVER_CONTEXT *context);


/**
 * Zero curve bootstrapping without currency basis.
 */
IrxTSwapZeroCurve* irxZeroCurve3
(IrxTBootstrapMethod method,
 IrxTMarketConv     *marketConv,
 IrxTDate            today,
 IrxTDate            baseDate,
 int                 numInstruments,
 char**              instTypes,
 IrxTDate*           maturityDates,
 double*             rates,
 double*             prices,
 long*               includeFlags,
 double*             adjustments,
 IrxTCalendar*       calendar,
 IrxTBool            valueFloating,
 int                 numShortRates,
 IrxTDate*           shortMaturityDates,
 double*             shortRates,
 long*               shortIncludeFlags,
 IrxTBool            firstFloatFixed,
 double              firstFloatFixedRate,
 IrxTStubLocation    stubLocation,
 IrxTDate            extrapDate)
{
    static char routine[] = "irxZeroCurve3";
    int         status    = FAILURE;

    IrxTZeroCurve *zc = NULL;
    IrxTSwapZeroCurve *swapCurve = NULL;

    if (validateDates (routine,
                       today,
                       baseDate,
                       numInstruments,
                       instTypes,
                       maturityDates,
                       includeFlags,
                       numShortRates,
                       shortMaturityDates,
                       shortIncludeFlags) != SUCCESS)
        goto RETURN; /* failure */

    REQUIRE(marketConv != NULL);
    REQUIRE(rates != NULL);
    if (numShortRates > 0)
    {
        REQUIRE (shortRates != NULL);
    }

    zc = buildFlatForwards(routine,
                           marketConv,
                           today,
                           baseDate,
                           numInstruments,
                           instTypes,
                           maturityDates,
                           rates,
                           prices,
                           includeFlags,
                           calendar,
                           valueFloating,
                           numShortRates,
                           shortMaturityDates,
                           shortRates,
                           shortIncludeFlags,
                           firstFloatFixed,
                           firstFloatFixedRate,
                           stubLocation,
                           extrapDate);
    if (zc == NULL)
        goto RETURN; /* failure */

    if (method == IRX_BOOTSTRAP_SMOOTH)
    {
        IrxTZeroCurve *tmp = buildSmoothForwards(routine,
                                                 zc,
                                                 marketConv,
                                                 today,
                                                 baseDate,
                                                 numInstruments,
                                                 instTypes,
                                                 maturityDates,
                                                 rates,
                                                 prices,
                                                 includeFlags,
                                                 adjustments,
                                                 calendar,
                                                 valueFloating,
                                                 numShortRates,
                                                 shortMaturityDates,
                                                 shortRates,
                                                 shortIncludeFlags,
                                                 firstFloatFixed,
                                                 firstFloatFixedRate,
                                                 stubLocation,
                                                 extrapDate);

        if (tmp == NULL)
            goto RETURN; /* failure */

        irxZeroCurveFree(zc);
        zc = tmp;
    }

    /********************************************************************
     * For transitional backward compatibility with wrapper and Fix3 
     * In principle, we should start with IrxTSwapZeroCurve rather than
     * IrxTZeroCurve.
     ********************************************************************/
    if (IrxTZeroCurveSetFromMarket(zc,
                                   today,
                                   marketConv) != SUCCESS)
        goto RETURN; /* failure */


    swapCurve = irxSwapZeroCurveMake (today,
                                      marketConv,
                                      calendar,
                                      zc,
                                      NULL);
    if (swapCurve == NULL)
        goto RETURN; /* failure */

    status = SUCCESS;

 RETURN:

    irxZeroCurveFree(zc);
    if (status != SUCCESS)
    {
        irxSwapZeroCurveFree(swapCurve);
        swapCurve = NULL;
        irxErrorFailure(routine);
    }

    return swapCurve;
}
 

/**
 * Zero curve bootstrapping with currency basis.
 *
 * Same as ZeroCurve3 but with two extra parameters (cbsSpreads, cbsCalendar)
 * and one less parameter (valueFloating).
 *
 * Also uses more information from the marketConv.
 */
IrxTSwapZeroCurve* irxZeroCurve3CB
(IrxTBootstrapMethod method,
 IrxTBool            reuseDiscountCurve,
 IrxTMarketConv     *marketConv,
 IrxTDate            today,
 IrxTDate            baseDate,
 int                 numInstruments,
 char**              instTypes,
 IrxTDate*           maturityDates,
 double*             rates,
 double*             cbsSpreads,
 double*             prices,
 long*               includeFlags,
 double*             adjustments,
 IrxTCalendar*       calendar,
 IrxTCalendar*       cbsCalendar,
 int                 numShortRates,
 IrxTDate*           shortMaturityDates,
 double*             shortRates,
 double*             shortCbsSpreads,
 long*               shortIncludeFlags,
 IrxTBool            firstFloatFixed,
 double              firstFloatFixedRate,
 IrxTStubLocation    stubLocation,
 IrxTDate            extrapDate)
{
    static char routine[] = "irxZeroCurve3CB";
    int         status    = FAILURE;

    IrxTZeroCurve     *zcDiscount = NULL;
    IrxTZeroCurve     *zcIndex = NULL;
    IrxTSwapZeroCurve *swapCurve = NULL;

    if (validateDates (routine,
                       today,
                       baseDate,
                       numInstruments,
                       instTypes,
                       maturityDates,
                       includeFlags,
                       numShortRates,
                       shortMaturityDates,
                       shortIncludeFlags) != SUCCESS)
        goto RETURN; /* failure */

    REQUIRE(marketConv != NULL);
    REQUIRE(rates != NULL);
    if (numShortRates > 0)
    {
        REQUIRE (shortRates != NULL);
    }

    if (buildFlatForwardsCB(routine,
                            reuseDiscountCurve,
                            marketConv,
                            today,
                            baseDate,
                            numInstruments,
                            instTypes,
                            maturityDates,
                            rates,
                            cbsSpreads,
                            prices,
                            includeFlags,
                            calendar,
                            cbsCalendar,
                            numShortRates,
                            shortMaturityDates,
                            shortRates,
                            shortCbsSpreads,
                            shortIncludeFlags,
                            firstFloatFixed,
                            firstFloatFixedRate,
                            stubLocation,
                            extrapDate,
                            &zcDiscount,
                            &zcIndex) != SUCCESS)
        goto RETURN; /* failure */


    if (method == IRX_BOOTSTRAP_SMOOTH)
    {
        IrxTZeroCurve *zcDiscount2;
        IrxTZeroCurve *zcIndex2;
        if (buildSmoothForwardsCB (routine,
                                   zcDiscount,
                                   zcIndex,
                                   reuseDiscountCurve,
                                   marketConv,
                                   today,
                                   baseDate,
                                   numInstruments,
                                   instTypes,
                                   maturityDates,
                                   rates,
                                   cbsSpreads,
                                   prices,
                                   includeFlags,
                                   adjustments,
                                   calendar,
                                   cbsCalendar,
                                   numShortRates,
                                   shortMaturityDates,
                                   shortRates,
                                   shortCbsSpreads,
                                   shortIncludeFlags,
                                   firstFloatFixed,
                                   firstFloatFixedRate,
                                   stubLocation,
                                   extrapDate,
                                   &zcDiscount2,
                                   &zcIndex2) != SUCCESS)
            goto RETURN; /* failure */

        irxZeroCurveFree(zcDiscount);
        irxZeroCurveFree(zcIndex);

        zcDiscount = zcDiscount2;
        zcIndex    = zcIndex2;
    }


    /********************************************************************
     * For transitional backward compatibility with wrapper and Fix3 
     * In principle, we should start with IrxTSwapZeroCurve rather than
     * IrxTZeroCurve.
     ********************************************************************/
    if (IrxTZeroCurveSetFromMarket(zcIndex,
                                   today,
                                   marketConv) != SUCCESS ||
        IrxTZeroCurveSetFromMarket(zcDiscount,
                                   today,
                                   marketConv) != SUCCESS)
        goto RETURN; /* failure */


    swapCurve = irxSwapZeroCurveMake (today,
                                      marketConv,
                                      calendar,
                                      zcDiscount,
                                      zcIndex);
    if (swapCurve == NULL)
        goto RETURN; /* failure */

    status = SUCCESS;

 RETURN:

    irxZeroCurveFree(zcIndex);
    irxZeroCurveFree(zcDiscount);
    if (status != SUCCESS)
    {
        irxSwapZeroCurveFree(swapCurve);
        swapCurve = NULL;
        irxErrorFailure(routine);
    }

    return swapCurve;
}
 


/* Validates that the dates are in ascending order and that the instrument
   types make sense */
static int validateDates
(char*               routine,
 IrxTDate            today,
 IrxTDate            baseDate,
 int                 numInstruments,
 char**              instTypes,
 IrxTDate*           maturityDates,
 long*               includeFlags,
 int                 numShortRates,
 IrxTDate*           shortMaturityDates,
 long*               shortIncludeFlags)
{
    int status = FAILURE;
    
    int i = 0;

    IrxTDate prevMaturityDate = 0;
    char     prevInstType     = 'M';
    char     dateBuf[2][16];

    REQUIRE (today <= baseDate);
    REQUIRE (numInstruments > 0);
    REQUIRE (numShortRates >= 0);
    REQUIRE (instTypes != NULL);
    REQUIRE (maturityDates != NULL);

    if (numShortRates > 0)
    {
        REQUIRE (shortMaturityDates != NULL);
    }

    /* Instrument dates must be in increasing order with the M's first
       and the S's second. */
    for (i = 0; i < numInstruments; ++i)
    {
        IrxTDate maturityDate;
        char*    instType;

        if (includeFlags && !includeFlags[i])
            continue; /* ignore when includeFlags[i] == 0 */

        maturityDate = maturityDates[i];
        instType     = instTypes[i];

        if (maturityDate <= baseDate)
        {
            irxError ("%s: Maturity date (%s) for instrument #%d must be "
                      "after base date (%s)\n", routine,
                      irxDateFormat(maturityDate, NULL, dateBuf[0]),
                      i+1,
                      irxDateFormat(baseDate, NULL, dateBuf[1]));
            goto RETURN; /* failure */
        }

        if (maturityDate <= prevMaturityDate)
        {
            irxError ("%s: Maturity date (%s) for instrument #%d must be "
                      "after previous date (%s)\n", routine,
                      irxDateFormat(maturityDate, NULL, dateBuf[0]),
                      i+1,
                      irxDateFormat(prevMaturityDate, NULL, dateBuf[1]));
            goto RETURN; /* failure */
        }
        prevMaturityDate = maturityDate;

        if (instType == NULL)
        {
            irxError ("%s: Undefined instrument type for instrument #%d\n",
                      routine, i+1);
            goto RETURN; /* failure */
        }

        switch (*instType)
        {
        case 'M':
            if (prevInstType != 'M')
            {
                irxError ("%s: Cannot define money market instrument #%d "
                          "after starting to define swaps\n", routine, i+1);
                goto RETURN; /* failure */
            }
            break;
        case 'S':
            prevInstType = 'S';
            break;
        default:
            irxError ("%s: Bad instrument type (%s) for instrument #%d\n",
                      routine, instType, i+1);
            goto RETURN; /* failure */
        }
    }

    /* Short rates must be increasing order between today and baseDate.
       Last date must equal baseDate. */

    prevMaturityDate = today;

    for (i = 0; i < numShortRates; ++i)
    {
        IrxTDate maturityDate;

        if (shortIncludeFlags && !shortIncludeFlags[i])
            continue; /* ignore when shortIncludeFlags[i] == 0 */

        maturityDate = shortMaturityDates[i];

        if (maturityDate <= prevMaturityDate)
        {
            irxError ("%s: Maturity date (%s) for short rate #%d must be "
                      "after previous date (%s)\n", routine,
                      irxDateFormat(maturityDate, NULL, dateBuf[0]),
                      i+1,
                      irxDateFormat(prevMaturityDate, NULL, dateBuf[1]));
            goto RETURN; /* failure */
        }
        prevMaturityDate = maturityDate;
    }

    if (prevMaturityDate != today)
    {
        if (prevMaturityDate != baseDate)
        {
            irxError ("%s: Maturity date (%s) for last short rate must be "
                      "equal to base date (%s)\n", routine,
                      irxDateFormat(prevMaturityDate, NULL, dateBuf[0]),
                      irxDateFormat(baseDate, NULL, dateBuf[1]));
            goto RETURN; /* failure */
        }
    }

    status = SUCCESS;

 RETURN:

    return status;
}

/* Count a set of rates given by a number of include flags */
static int countRates
(int   numRates,
 long *includeFlags,
 long  maxInclude)
{
    int i = 0;
    int count = 0;

    if (includeFlags == NULL)
        return numRates;

    for (i = 0; i < numRates; ++i)
    {
        if (includeFlags[i] > 0 && includeFlags[i] <= maxInclude)
            ++count;
    }

    return count;
}


/* Builds the short end of the curve (with or without currency basis) */
static IrxTZeroCurve* buildShortEnd
(char*               routine,
 IrxTMarketConv     *marketConv,
 IrxTDate            today,
 IrxTDate            baseDate,
 int                 numInstruments,
 int                 numShortRates,
 IrxTDate*           shortMaturityDates,
 double*             shortRates,
 double*             shortCbsSpreads,
 long*               shortIncludeFlags)
{
    int status = FAILURE;

    int numShortRatesIncluded = countRates(numShortRates,
                                           shortIncludeFlags,
                                           1);
                                           
    int numItems = numInstruments + 1 + numShortRatesIncluded;
    int i;
    int j;
    int k;
    double   discount;
    double   cbsDiscount;
    IrxTDate startDate;

    IrxTZeroCurve* zc = irxZeroCurveMakeEmpty(numItems);

    if (zc == NULL)
        goto RETURN; /* failure */

    zc->baseDate = baseDate;
    j = numShortRatesIncluded;
    i = numShortRates;

    zc->startDates[j] = baseDate;
    zc->prices[j]     = 1.0;

    while (i > 0)
    {
        --i;
        if (shortIncludeFlags != NULL && shortIncludeFlags[i] != 1)
            continue;

        --j;
        k = i;
        startDate = today;
        while (k > 0)
        {
            --k;
            if (shortIncludeFlags != NULL && shortIncludeFlags[k] != 1)
                continue;
            startDate = shortMaturityDates[k];
            break;
        }

        if (irxRateToDiscount (shortRates[i],
                               startDate,
                               shortMaturityDates[i],
                               marketConv->mmDcc,
                               IRX_SIMPLE_RATE,
                               &discount) != SUCCESS)
            goto RETURN; /* failure */

        if (shortCbsSpreads != NULL)
        {
            if (irxRateToDiscount (shortCbsSpreads[i],
                                   startDate,
                                   shortMaturityDates[i],
                                   marketConv->cbsDcc,
                                   IRX_SIMPLE_RATE,
                                   &cbsDiscount) != SUCCESS)
                goto RETURN; /* failure */
        }
        else
        {
            cbsDiscount = 1.0;
        }

        zc->startDates[j] = startDate;
        zc->prices[j]     = zc->prices[j+1] / discount / cbsDiscount;
        
        if (irxDiscountToRate (discount * cbsDiscount,
                               startDate,
                               shortMaturityDates[i],
                               IRX_ACT_365F,
                               IRX_CONTINUOUS_RATE,
                               &zc->a0[j]) != SUCCESS)
            goto RETURN; /* failure */
    }

    ASSERT(j == 0);

    status = SUCCESS;

 RETURN:

    if (status != SUCCESS)
    {
        irxZeroCurveFree(zc);
        zc = NULL;
    }
    return zc;
}







/* Bootstraps a zero curve using flat forwards without currency basis */
static IrxTZeroCurve* buildFlatForwards
(char*               routine,
 IrxTMarketConv     *marketConv,
 IrxTDate            today,
 IrxTDate            baseDate,
 int                 numInstruments,
 char**              instTypes,
 IrxTDate*           maturityDates,
 double*             rates,
 double*             prices,
 long*               includeFlags,
 IrxTCalendar*       calendar,
 IrxTBool            valueFloating,
 int                 numShortRates,
 IrxTDate*           shortMaturityDates,
 double*             shortRates,
 long*               shortIncludeFlags,
 IrxTBool            firstFloatFixed,
 double              firstFloatFixedRate,
 IrxTStubLocation    stubLocation,
 IrxTDate            extrapDate)
{
    int status = FAILURE;

    IrxTZeroCurve *zc = NULL;
    int i;
    int j;

    int numInstrumentsToInclude = countRates(numInstruments,
                                             includeFlags,
                                             1);
    IrxTDate maxDate;

    zc = buildShortEnd (routine,
                        marketConv,
                        today,
                        baseDate,
                        numInstrumentsToInclude,
                        numShortRates,
                        shortMaturityDates,
                        shortRates,
                        NULL, /* shortCbsSpreads */
                        shortIncludeFlags);
    if (zc == NULL)
        goto RETURN; /* failure */

    j = zc->numItems - numInstrumentsToInclude - 1;
    ASSERT (j >= 0);
    ASSERT (zc->startDates[j] == baseDate);
    zc->numItems = j+1;

    for (i = 0; i < numInstruments; ++i)
    {
        IrxTDate maturityDate;
        char*    instType;
        double   price;

        if (includeFlags != NULL && includeFlags[i] != 1)
            continue;

        maturityDate = maturityDates[i];
        instType     = instTypes[i];
        switch (*instType)
        {
        case 'M':
            if (addCashRate(routine, zc, j, FALSE, FALSE, FALSE,
                            marketConv, baseDate,
                            maturityDate, calendar,
                            rates[i], 0.0) != SUCCESS)
                goto RETURN; /* failure */
            break;
        case 'S':
            if (prices != NULL)
                price = prices[i];
            else
                price = 1.0;
            if (addSwapRate(routine, zc, zc, j, FALSE, FALSE, FALSE,
                            marketConv, baseDate,
                            maturityDate, calendar, valueFloating,
                            firstFloatFixed, firstFloatFixedRate, stubLocation,
                            rates[i], price) != SUCCESS)
                goto RETURN; /* failure */
            break;
        default:
            PROGRAM_BUG();
            goto RETURN; /* failure */
        }
        ++j;
        zc->numItems = j+1;
    }

    if (zc->numItems > 1)
        zc->a0[zc->numItems-1] = zc->a0[zc->numItems-2];

    maxDate = MAX(zc->startDates[zc->numItems-1], extrapDate);
    zc->maxDate = maxDate;

    status = SUCCESS;

 RETURN:

    if (status != SUCCESS)
    {
        irxZeroCurveFree (zc);
        zc = NULL;
    }
    return zc;
}

/* This does the following:

   discountCurve = add spread to rate, compute curve based on these swaps
   indexCurve    = compute curve such that fixed=floating using discountCurve
                   to discount them.
*/

static int buildFlatForwardsCB
(char*               routine,
 IrxTBool            reuseDiscountCurve,
 IrxTMarketConv     *marketConv,
 IrxTDate            today,
 IrxTDate            baseDate,
 int                 numInstruments,
 char**              instTypes,
 IrxTDate*           maturityDates,
 double*             rates,
 double*             cbsSpreads,
 double*             prices,
 long*               includeFlags,
 IrxTCalendar*       calendar,
 IrxTCalendar*       cbsCalendar,
 int                 numShortRates,
 IrxTDate*           shortMaturityDates,
 double*             shortRates,
 double*             shortCbsSpreads,
 long*               shortIncludeFlags,
 IrxTBool            firstFloatFixed,
 double              firstFloatFixedRate,
 IrxTStubLocation    stubLocation,
 IrxTDate            extrapDate,
 IrxTZeroCurve     **out_zcDisc,
 IrxTZeroCurve     **out_zcIndex
)
{
    int status = FAILURE;

    IrxTZeroCurve *zcDisc = NULL;
    IrxTZeroCurve *zcIndex = NULL;
    IrxTZeroCurve *zcDiscUsed;
    int i;
    int j;
    int k;

    int numInstrumentsToInclude = countRates(numInstruments,
                                             includeFlags,
                                             1);
    IrxTDate maxDate;

    zcDisc = buildShortEnd (routine,
                            marketConv,
                            today,
                            baseDate,
                            numInstrumentsToInclude,
                            numShortRates,
                            shortMaturityDates,
                            shortRates,
                            shortCbsSpreads,
                            shortIncludeFlags);
    if (zcDisc == NULL)
        goto RETURN; /* failure */

    zcIndex = buildShortEnd (routine,
                             marketConv,
                             today,
                             baseDate,
                             numInstrumentsToInclude,
                             numShortRates,
                             shortMaturityDates,
                             shortRates,
                             NULL,
                             shortIncludeFlags);
    if (zcIndex == NULL)
        goto RETURN; /* failure */

    j = zcDisc->numItems - numInstrumentsToInclude - 1;
    k = zcIndex->numItems - numInstrumentsToInclude - 1;
    ASSERT (j >= 0);
    ASSERT (j == k);
    ASSERT (zcDisc->startDates[j] == baseDate);
    ASSERT (zcIndex->startDates[k] == baseDate);

    /* Do the discount curve first */
    for (i = 0; i < numInstruments; ++i)
    {
        IrxTDate maturityDate;
        char*    instType;
        double   price;

        if (includeFlags != NULL && includeFlags[i] != 1)
            continue;

        maturityDate = maturityDates[i];
        instType     = instTypes[i];
        switch (*instType)
        {
        case 'M':
            if (addCashRate(routine, zcDisc, j, FALSE, FALSE, FALSE,
                            marketConv, baseDate,
                            maturityDate, calendar,
                            rates[i], cbsSpreads[i]) != SUCCESS)
                goto RETURN; /* failure */
            break;
        case 'S':
            if (prices != NULL)
                price = prices[i];
            else
                price = 1.0;
            if (addSwapRateCB(routine, zcDisc, j, FALSE, FALSE, FALSE,
                              marketConv, baseDate,
                              maturityDate, calendar, cbsCalendar,
                              stubLocation, rates[i],
                              cbsSpreads[i], price) != SUCCESS)
                goto RETURN; /* failure */
            break;
        default:
            PROGRAM_BUG();
            goto RETURN; /* failure */
        }
        ++j;
        zcDisc->numItems = j+1;
    }

    /* Just in case we have some little glitch at the end regarding
       date differences between zcDisc and zcIndex ensure that we can
       extrapolate the zcDisc */
    zcDisc->maxDate = zcDisc->startDates[zcDisc->numItems-1] + 100;

    if (reuseDiscountCurve)
        zcDiscUsed = zcDisc;
    else
        zcDiscUsed = zcIndex;

    /* Do the index curve second */
    for (i = 0; i < numInstruments; ++i)
    {
        IrxTDate maturityDate;
        char*    instType;

        if (includeFlags != NULL && includeFlags[i] != 1)
            continue;

        maturityDate = maturityDates[i];
        instType     = instTypes[i];
        switch (*instType)
        {
        case 'M':
            if (addCashRate(routine, zcIndex, k, FALSE, FALSE, FALSE,
                            marketConv, baseDate,
                            maturityDate, calendar,
                            rates[i], 0.0) != SUCCESS)
                goto RETURN; /* failure */
            break;
        case 'S':
            if (addSwapRate(routine, zcIndex, zcDiscUsed, k, 
                            FALSE, FALSE, FALSE, marketConv,
                            baseDate, maturityDate, calendar, TRUE,
                            firstFloatFixed, firstFloatFixedRate, 
                            stubLocation, rates[i], 0.0) != SUCCESS)
                goto RETURN; /* failure */
            break;
        default:
            PROGRAM_BUG();
            goto RETURN; /* failure */
        }
        ++k;
        zcIndex->numItems = k+1;
    }

    ASSERT(zcDisc->numItems == zcIndex->numItems);

    maxDate = MAX(zcDisc->startDates[zcDisc->numItems-1],
                  zcIndex->startDates[zcIndex->numItems-1]);
    maxDate = MAX(maxDate, extrapDate);

    zcDisc->maxDate  = maxDate;
    zcIndex->maxDate = maxDate;

    *out_zcDisc  = zcDisc;
    *out_zcIndex = zcIndex;

    zcDisc  = NULL;
    zcIndex = NULL;
        
    status = SUCCESS;

 RETURN:

    irxZeroCurveFree (zcIndex);
    irxZeroCurveFree (zcDisc);
    return status;
}



/* Bootstraps a zero curve using smooth forwards without currency basis */
static IrxTZeroCurve* buildSmoothForwards
(char*               routine,
 IrxTZeroCurve      *zcFlat,
 IrxTMarketConv     *marketConv,
 IrxTDate            today,
 IrxTDate            baseDate,
 int                 numInstruments,
 char**              instTypes,
 IrxTDate*           maturityDates,
 double*             rates,
 double*             prices,
 long*               includeFlags,
 double*             adjustments,
 IrxTCalendar*       calendar,
 IrxTBool            valueFloating,
 int                 numShortRates,
 IrxTDate*           shortMaturityDates,
 double*             shortRates,
 long*               shortIncludeFlags,
 IrxTBool            firstFloatFixed,
 double              firstFloatFixedRate,
 IrxTStubLocation    stubLocation,
 IrxTDate            extrapDate)
{
    int status = FAILURE;

    IrxTZeroCurve *zc = NULL;
    IrxTZeroCurve *zcFlatPhaseTwo = NULL;
    long   *includeFlagsPhaseTwo  = NULL;
    double *ratesPhaseTwo         = NULL;

    double  price;

    int i;
    int j;
    int lastPoint;
    IrxTBool first;
    IrxTBool last;
    IrxTDate maxDate;

    int numPhaseOneInstruments = countRates(numInstruments,
                                            includeFlags,
                                            1);
    int numPhaseTwoInstruments = countRates(numInstruments,
                                            includeFlags,
                                            2);
    int numInstrumentsToInclude;

    if (numPhaseTwoInstruments > numPhaseOneInstruments)
    {
        ASSERT (includeFlags != NULL);
        REQUIRE (adjustments != NULL);
        includeFlagsPhaseTwo = NEW_ARRAY(long, numInstruments);
        ratesPhaseTwo        = NEW_ARRAY(double, numInstruments);
        COPY_ARRAY (includeFlagsPhaseTwo, includeFlags, long, numInstruments);
        COPY_ARRAY (ratesPhaseTwo, rates, double, numInstruments);
        for (i = 0; i < numInstruments; ++i)
        {
            if (includeFlags[i] == 2)
            {
                includeFlagsPhaseTwo[i] = 1;
                ratesPhaseTwo[i]        = rates[i] + adjustments[i];
            }
        }
        zcFlatPhaseTwo = buildFlatForwards(routine,
                                           marketConv,
                                           today,
                                           baseDate,
                                           numInstruments,
                                           instTypes,
                                           maturityDates,
                                           ratesPhaseTwo,
                                           prices,
                                           includeFlagsPhaseTwo,
                                           calendar,
                                           valueFloating,
                                           numShortRates,
                                           shortMaturityDates,
                                           shortRates,
                                           shortIncludeFlags,
                                           firstFloatFixed,
                                           firstFloatFixedRate,
                                           stubLocation,
                                           extrapDate);
        if (zcFlatPhaseTwo == NULL)
            goto RETURN; /* failure */

        zcFlat       = zcFlatPhaseTwo;
        includeFlags = includeFlagsPhaseTwo;
        rates        = ratesPhaseTwo;
    }

    numInstrumentsToInclude = countRates(numInstruments,
                                         includeFlags,
                                         1);

    ASSERT(numInstrumentsToInclude == numPhaseTwoInstruments);

    /* now we have a curve - we calculate the annuity of each maturity date */

    /* the short end of the curve does not have any smoothing */
    /* after all we don't expect any valid business dates except those in
       the short end of the curve */
    zc = buildShortEnd (routine,
                        marketConv,
                        today,
                        baseDate,
                        numInstrumentsToInclude,
                        numShortRates,
                        shortMaturityDates,
                        shortRates,
                        NULL, /* shortCbsSpreads */
                        shortIncludeFlags);
    if (zc == NULL)
        goto RETURN; /* failure */

    j = zc->numItems - numInstrumentsToInclude - 1;
    if (initSmoothCurve (zc, 
                         zcFlat,
                         numInstruments,
                         maturityDates,
                         includeFlags,
                         j) != SUCCESS)
        goto RETURN; /* failure */

    ASSERT (j >= 0);
    ASSERT (zc->startDates[j] == baseDate);
    
    /* When we process point j we fill in a0[j-1], a1[j-1], a2[j-1]
       and prices[j].

       We use knowledge of a0[j-1] and a0[j] to apply constraints. */

    first = TRUE;
    last  = FALSE;
    lastPoint = zc->numItems-2; /* tends to change as we go along */
    for (i = 0; i < numInstruments; ++i)
    {
        IrxTDate maturityDate;
        char*    instType;

        if (includeFlags != NULL && includeFlags[i] != 1)
            continue;

        maturityDate = maturityDates[i];
        instType     = instTypes[i];

        last = (j == lastPoint);
        switch (*instType)
        {
        case 'M':
            if (addCashRate(routine, zc, j, TRUE, first, last,
                            marketConv, baseDate,
                            maturityDate, calendar,
                            rates[i], 0.0) != SUCCESS)
                goto RETURN; /* failure */
            break;
        case 'S':
            if (prices != NULL)
                price = prices[i];
            else
                price = 1.0;
            if (addSwapRate(routine, zc, zc, j, TRUE, first, last,
                            marketConv, baseDate,
                            maturityDate, calendar, valueFloating,
                            firstFloatFixed, firstFloatFixedRate,
                            stubLocation, rates[i], price) != SUCCESS)
                goto RETURN; /* failure */
            break;
        default:
            PROGRAM_BUG();
            goto RETURN; /* failure */
        }
        ++j;
        first = FALSE;
        zc->numItems = j+1;
    }

    maxDate = MAX(zc->startDates[zc->numItems-1], extrapDate);
    zc->maxDate = maxDate;

    status = SUCCESS;

 RETURN:

    FREE (includeFlagsPhaseTwo);
    FREE (ratesPhaseTwo);
    irxZeroCurveFree (zcFlatPhaseTwo);

    if (status != SUCCESS)
    {
        irxZeroCurveFree (zc);
        zc = NULL;
    }
    return zc;
}

/* Bootstraps a zero curve using smooth forwards with currency basis */
static int buildSmoothForwardsCB
(char*               routine,
 IrxTZeroCurve      *zcDiscFlat,
 IrxTZeroCurve      *zcIndexFlat,
 IrxTBool            reuseDiscountCurve,
 IrxTMarketConv     *marketConv,
 IrxTDate            today,
 IrxTDate            baseDate,
 int                 numInstruments,
 char**              instTypes,
 IrxTDate*           maturityDates,
 double*             rates,
 double*             cbsSpreads,
 double*             prices,
 long*               includeFlags,
 double*             adjustments,
 IrxTCalendar*       calendar,
 IrxTCalendar*       cbsCalendar,
 int                 numShortRates,
 IrxTDate*           shortMaturityDates,
 double*             shortRates,
 double*             shortCbsSpreads,
 long*               shortIncludeFlags,
 IrxTBool            firstFloatFixed,
 double              firstFloatFixedRate,
 IrxTStubLocation    stubLocation,
 IrxTDate            extrapDate,
 IrxTZeroCurve**     out_zcDisc,
 IrxTZeroCurve**     out_zcIndex)
{
    int status = FAILURE;

    IrxTZeroCurve *zcDisc = NULL;
    IrxTZeroCurve *zcIndex = NULL;
    IrxTZeroCurve *zcDiscFlatPhaseTwo = NULL;
    IrxTZeroCurve *zcIndexFlatPhaseTwo = NULL;
    IrxTZeroCurve *zcDiscUsed;
    long   *includeFlagsPhaseTwo  = NULL;
    double *ratesPhaseTwo         = NULL;

    IrxTBool first;
    IrxTBool last;
    int      lastPoint;

    int i;
    int j;
    int k;
    IrxTDate maxDate;

    int numPhaseOneInstruments = countRates(numInstruments,
                                            includeFlags,
                                            1);
    int numPhaseTwoInstruments = countRates(numInstruments,
                                            includeFlags,
                                            2);
    int numInstrumentsToInclude;

    if (numPhaseTwoInstruments > numPhaseOneInstruments)
    {
        ASSERT (includeFlags != NULL);
        REQUIRE (adjustments != NULL);
        includeFlagsPhaseTwo = NEW_ARRAY(long, numInstruments);
        ratesPhaseTwo        = NEW_ARRAY(double, numInstruments);
        COPY_ARRAY (includeFlagsPhaseTwo, includeFlags, long, numInstruments);
        COPY_ARRAY (ratesPhaseTwo, rates, double, numInstruments);
        for (i = 0; i < numInstruments; ++i)
        {
            if (includeFlags[i] == 2)
            {
                includeFlagsPhaseTwo[i] = 1;
                ratesPhaseTwo[i]        = rates[i] + adjustments[i];
            }
        }
        if (buildFlatForwardsCB(routine,
                                reuseDiscountCurve,
                                marketConv,
                                today,
                                baseDate,
                                numInstruments,
                                instTypes,
                                maturityDates,
                                ratesPhaseTwo,
                                cbsSpreads,
                                prices,
                                includeFlagsPhaseTwo,
                                calendar,
                                cbsCalendar,
                                numShortRates,
                                shortMaturityDates,
                                shortRates,
                                shortCbsSpreads,
                                shortIncludeFlags,
                                firstFloatFixed,
                                firstFloatFixedRate,
                                stubLocation,
                                extrapDate,
                                &zcDiscFlatPhaseTwo,
                                &zcIndexFlatPhaseTwo) != SUCCESS)
            goto RETURN; /* failure */

        zcDiscFlat   = zcDiscFlatPhaseTwo;
        zcIndexFlat  = zcIndexFlatPhaseTwo;
        includeFlags = includeFlagsPhaseTwo;
        rates        = ratesPhaseTwo;
    }

    numInstrumentsToInclude = countRates(numInstruments,
                                         includeFlags,
                                         1);

    ASSERT(numInstrumentsToInclude == numPhaseTwoInstruments);

    /* now we have a curve - we calculate the annuity of each maturity date */

    /* the short end of the curve does not have any smoothing */
    /* after all we don't expect any valid business dates except those in
       the short end of the curve */
    zcDisc = buildShortEnd (routine,
                            marketConv,
                            today,
                            baseDate,
                            numInstrumentsToInclude,
                            numShortRates,
                            shortMaturityDates,
                            shortRates,
                            shortCbsSpreads,
                            shortIncludeFlags);
    if (zcDisc == NULL)
        goto RETURN; /* failure */

    zcIndex = buildShortEnd (routine,
                             marketConv,
                             today,
                             baseDate,
                             numInstrumentsToInclude,
                             numShortRates,
                             shortMaturityDates,
                             shortRates,
                             NULL,
                             shortIncludeFlags);
    if (zcIndex == NULL)
        goto RETURN; /* failure */

    j = zcDisc->numItems - numInstrumentsToInclude - 1;
    k = zcIndex->numItems - numInstrumentsToInclude - 1;

    if (initSmoothCurve (zcDisc,
                         zcDiscFlat,
                         numInstruments,
                         maturityDates,
                         includeFlags,
                         j) != SUCCESS)
        goto RETURN; /* failure */

    if (initSmoothCurve (zcIndex,
                         zcIndexFlat,
                         numInstruments,
                         maturityDates,
                         includeFlags,
                         k) != SUCCESS)
        goto RETURN; /* failure */

    ASSERT (j >= 0);
    ASSERT (j == k);
    ASSERT (zcDisc->startDates[j] == baseDate);
    ASSERT (zcIndex->startDates[k] == baseDate);

    /* now we have the following state - we have a0 set for all 
       intermediate points - i.e. for the baseDate and last point the
       value of a0 is not set */

    /* Do the discount curve first */
    first = TRUE;
    last  = FALSE;
    lastPoint = zcDisc->numItems-2; /* tends to change as we go along */

    for (i = 0; i < numInstruments; ++i)
    {
        IrxTDate maturityDate;
        char*    instType;
        double   price;

        if (includeFlags != NULL && includeFlags[i] != 1)
            continue;

        maturityDate = maturityDates[i];
        instType     = instTypes[i];
        last         = j == lastPoint;
        switch (*instType)
        {
        case 'M':
            if (addCashRate(routine, zcDisc, j, TRUE, first, last,
                            marketConv, baseDate,
                            maturityDate, calendar,
                            rates[i], cbsSpreads[i]) != SUCCESS)
                goto RETURN; /* failure */
            break;
        case 'S':
            if (prices != NULL)
                price = prices[i];
            else
                price = 1.0;
            if (addSwapRateCB(routine, zcDisc, j, TRUE, first, last,
                              marketConv, baseDate,
                              maturityDate, calendar, cbsCalendar,
                              stubLocation, rates[i],
                              cbsSpreads[i], price) != SUCCESS)
                goto RETURN; /* failure */
            break;
        default:
            PROGRAM_BUG();
            goto RETURN; /* failure */
        }
        first = FALSE;
        ++j;
        zcDisc->numItems = j+1;
    }

    /* Just in case we have some little glitch at the end regarding
       date differences between zcDisc and zcIndex ensure that we can
       extrapolate the zcDisc */
    zcDisc->maxDate = zcDisc->startDates[zcDisc->numItems-1] + 100;

    if (reuseDiscountCurve)
        zcDiscUsed = zcDisc;
    else
        zcDiscUsed = zcIndex;

    /* Do the index curve second */
    first = TRUE;
    last  = FALSE;
    lastPoint = zcIndex->numItems-2; /* tends to change as we go along */

    for (i = 0; i < numInstruments; ++i)
    {
        IrxTDate maturityDate;
        char*    instType;

        if (includeFlags != NULL && includeFlags[i] != 1)
            continue;

        maturityDate = maturityDates[i];
        instType     = instTypes[i];
        last         = k == lastPoint;
        switch (*instType)
        {
        case 'M':
            if (addCashRate(routine, zcIndex, k, TRUE, first, last,
                            marketConv, baseDate,
                            maturityDate, calendar,
                            rates[i], 0.0) != SUCCESS)
                goto RETURN; /* failure */
            break;
        case 'S':
            if (addSwapRate(routine, zcIndex, zcDiscUsed, k, 
                            TRUE, first, last, marketConv,
                            baseDate, maturityDate, calendar, TRUE,
                            firstFloatFixed, firstFloatFixedRate, 
                            stubLocation, rates[i], 0.0) != SUCCESS)
                goto RETURN; /* failure */
            break;
        default:
            PROGRAM_BUG();
            goto RETURN; /* failure */
        }
        first = FALSE;
        ++k;
        zcIndex->numItems = k+1;
    }

    maxDate = MAX(zcDisc->startDates[zcDisc->numItems-1],
                  zcIndex->startDates[zcIndex->numItems-1]);
    maxDate = MAX(maxDate, extrapDate);

    zcDisc->maxDate  = maxDate;
    zcIndex->maxDate = maxDate;

    *out_zcDisc  = zcDisc;
    *out_zcIndex = zcIndex;
            
    zcDisc  = NULL;
    zcIndex = NULL;

    status = SUCCESS;

 RETURN:

    FREE (includeFlagsPhaseTwo);
    FREE (ratesPhaseTwo);
    irxZeroCurveFree (zcDiscFlatPhaseTwo);
    irxZeroCurveFree (zcIndexFlatPhaseTwo);
    irxZeroCurveFree (zcDisc);
    irxZeroCurveFree (zcIndex);

    return status;
}


/* Calculates the duration for a given maturity date */
/* Based on a mythical market which has particular market conventions 
   and an existing zero curve. */
static int calcDuration (
    IrxTZeroCurve const *zc,
    IrxTDate             maturityDate,
    double              *duration)
{
    static char routine[] = "calcDuration";
    int status = FAILURE;

    double        years;
    double        freq = 2.0;
    long          days30360;

    IrxTSwap      swap;
    double        swapRate;

    static IrxTMarketConv usdMarketConvNoBadDays =
    {
        {6, 'M', FALSE},  /* fixedIvl */
        {3, 'M', FALSE},  /* floatIvl */
        {6, 'M', FALSE},  /* cbsIvl */
        IRX_B30_360,      /* fixedDcc */
        IRX_ACT_360,      /* floatDcc */
        IRX_B30_360,      /* cbsDcc */
        IRX_ACT_360,      /* mmDcc */
        IRX_BAD_DAY_NONE, /* paymentBdc */
        IRX_BAD_DAY_NONE, /* accrualBdc */
        IRX_BAD_DAY_NONE, /* resetBdc */
        IRX_BAD_DAY_NONE, /* mmBdc */
        2                 /* daysToSpot */
    };

    IrxTMarketConv *marketConv = &usdMarketConvNoBadDays;

    if (irxDayCountDays(zc->baseDate, maturityDate, IRX_B30_360,
                        &days30360) != SUCCESS)
        goto RETURN; /* failure */

    if (days30360 == 0)
    {
        *duration = 0.0;
    }
    else
    {
        memset(&swap, 0, sizeof(IrxTSwap));
        swap.marketConv       = marketConv;
        swap.couponRate       = 0.0;
        swap.spread           = 0.0;
        swap.startDate        = zc->baseDate;
        swap.maturityDate     = maturityDate;
        swap.rollDate         = 0;
        swap.stubLocation     = IRX_SHORT_FRONT_STUB;
        swap.fixedStubPayment = IRX_STUB_SIMPLE;
        swap.floatStubPayment = IRX_STUB_SIMPLE;

        if (irxSwapRate (&swap,
                         zc,
                         FALSE, /* valueFloating */
                         NULL,  /* indexCurve */
                         1.0,   /* pv */
                         NULL,  /* calendar */
                         0,     /* firstFloatFixedDate */
                         0.0,   /* firstFloatFixedRate */
                         &swapRate) != SUCCESS)
            goto RETURN; /* failure */

        if (irxDayCountFraction (zc->baseDate,
                                 maturityDate,
                                 IRX_ACT_365F,
                                 &years) != SUCCESS)
            goto RETURN; /* failure */
        
        if (IS_ALMOST_ZERO(swapRate))
        {
            *duration = years;
        }
        else
        {
            *duration = (1.0 - pow(1.0+swapRate/freq, -years*freq)) / swapRate;
        }
    }
    status = SUCCESS;

 RETURN:

    if (status != SUCCESS)
        irxErrorFailure (routine);

    return status;
}

static int initSmoothCurve
(IrxTZeroCurve*       zcSmooth,
 IrxTZeroCurve const* zcFlat,
 int                  numInstruments,
 IrxTDate*            maturityDates,
 long*                includeFlags,
 int                  baseDateIndex)
{
    static char routine[] = "initSmoothCurve";
    int         status    = FAILURE;
    
    int     i,j,k;
    double  durationM2;
    double  durationM1;
    double  duration;

    j = baseDateIndex;
    k = -1; /* the previous point */
    durationM1 = 0.0;
    durationM2 = 0.0;

    /* fill in the values of a0 for all intermediate points */
    for (i = 0; i < numInstruments; ++i)
    {
        IrxTDate maturityDate;

        if (includeFlags != NULL && includeFlags[i] != 1)
            continue;

        ++j;
        maturityDate = maturityDates[i];

        if (calcDuration (zcFlat, maturityDate, &duration) != SUCCESS)
            goto RETURN; /* failure */

        zcSmooth->startDates[j] = zcFlat->startDates[j];
        if (k > 0)
        {
            double d2 = duration-durationM1;
            double d1 = durationM1-durationM2;
            double f2 = zcFlat->a0[k];
            double f1 = zcFlat->a0[k-1];
            double a0 = (f2*d1 + f1*d2) / (d1+d2);
            zcSmooth->a0[k] = a0;
        }
        k = j;
        
        /* shuffle up the variables */
        durationM2   = durationM1;
        durationM1   = duration;
    }

    /* now we have the following state - we have a0 set for all 
       intermediate points - i.e. for the baseDate and last point the
       value of a0 is not set */

    ASSERT (zcSmooth->startDates[baseDateIndex] == zcSmooth->baseDate);
    ASSERT (IS_ZERO(zcSmooth->a0[baseDateIndex]));
    ASSERT (IS_ZERO(zcSmooth->a0[zcSmooth->numItems-1]));

    status = SUCCESS;

 RETURN:

    return status;
}




/* Adds a cash rate (with or without currency basis) */

static int flatSolver (double rate, SOLVER_CONTEXT *context)
{
    double a0;

    int j             = context->j;
    IrxTZeroCurve *zc = context->zc;
    double time       = context->time;

    if (irxRootFindBrent ((IrxTObjectFunc)flatObjFunc,
                          context,
                          LOWER_BOUND,
                          UPPER_BOUND,
                          MAX_ITERATIONS,
                          rate, /* a not very good guess */
                          INITIAL_X_STEP,
                          INITIAL_F_DERIV,
                          X_TOLERANCE,
                          F_TOLERANCE,
                          &a0) != SUCCESS)
        return FAILURE;

    zc->a0[j]       = a0;
    zc->a1[j]       = 0.0;
    zc->a2[j]       = 0.0;
    zc->a0[j+1]     = a0;
    zc->prices[j+1] = zc->prices[j] * exp(-a0*time);

    return SUCCESS;
}
    
static int smoothSolverFirstPoint (double rate, SOLVER_CONTEXT *context)
{
    double a0;
    double a2;
    double discount;
    double guess;

    int j             = context->j;
    IrxTZeroCurve *zc = context->zc;
    double time       = context->time;

    if (!context->valueFloating && context->fixedFlows->numItems == 1)
    {
        /* we can solve directly */
        /* in case there is a date difference we will use this as a guess */
        discount = context->pv / context->fixedFlows->amounts[0];
        guess = -0.5 * (zc->a0[j+1] + 3.0 * log(discount)/context->time);
    }
    else
    {
        guess = rate;
    }

    if (irxRootFindBrent ((IrxTObjectFunc)smoothObjFuncFirstPoint,
                          context,
                          LOWER_BOUND,
                          UPPER_BOUND,
                          MAX_ITERATIONS,
                          guess,
                          INITIAL_X_STEP,
                          INITIAL_F_DERIV,
                          X_TOLERANCE,
                          F_TOLERANCE,
                          &a0) != SUCCESS)
        return FAILURE;
    
    smoothCalcFirst (a0, zc->a0[j+1], time, &a2, &discount);
    zc->a0[j]       = a0;
    zc->a1[j]       = 0;
    zc->a2[j]       = a2;
    zc->prices[j+1] = zc->prices[j] * discount;
    
    return SUCCESS;
}

static int smoothSolverMiddlePoint (SOLVER_CONTEXT *context)
{
    double a1;
    double a2;
    double discount;
    double guess;

    int j             = context->j;
    IrxTZeroCurve *zc = context->zc;
    double time       = context->time;

    if (!context->valueFloating && context->fixedFlows->numItems == 1)
    {
        /* we can solve directly */
        /* in case there is a date difference we will use this as a guess */
        discount =  context->pv / context->fixedFlows->amounts[0];
        discount /= zc->prices[j];
        guess    = -2.0 * (3.0*log(discount)/context->time +
                           zc->a0[j+1] + 2*zc->a0[j]) / context->time;
    }
    else
    {
        /* a linear guess */
        guess = (zc->a0[j+1] - zc->a0[j]) / context->time;
    }

    if (irxRootFindBrent ((IrxTObjectFunc)smoothObjFuncMiddlePoint,
                          context,
                          LOWER_BOUND,
                          UPPER_BOUND,
                          MAX_ITERATIONS,
                          guess,
                          INITIAL_X_STEP,
                          INITIAL_F_DERIV,
                          X_TOLERANCE,
                          F_TOLERANCE,
                          &a1) != SUCCESS)
        return FAILURE;
    
    smoothCalcMiddle (zc->a0[j], zc->a0[j+1], a1, time, &a2, &discount);
    zc->a1[j]       = a1;
    zc->a2[j]       = a2;
    zc->prices[j+1] = zc->prices[j] * discount;
    
    return SUCCESS;
}


static int smoothSolverLastPoint (SOLVER_CONTEXT *context)
{
    double a1;
    double a2;
    double discount;
    double guess;

    int j             = context->j;
    IrxTZeroCurve *zc = context->zc;
    double time       = context->time;

    if (!context->valueFloating && context->fixedFlows->numItems == 1)
    {
        /* we can solve directly */
        /* in case there is a date difference we will use this as a guess */
        discount =  context->pv / context->fixedFlows->amounts[0];
        discount /= zc->prices[j];
        guess    = -3.0 * (log(discount)/context->time + zc->a0[j]) / 
            context->time;
    }
    else
    {
        /* a flat guess */
        guess = 0.0;
    }

    if (irxRootFindBrent ((IrxTObjectFunc)smoothObjFuncLastPoint,
                          context,
                          LOWER_BOUND,
                          UPPER_BOUND,
                          MAX_ITERATIONS,
                          guess,
                          INITIAL_X_STEP,
                          INITIAL_F_DERIV,
                          X_TOLERANCE,
                          F_TOLERANCE,
                          &a1) != SUCCESS)
        return FAILURE;
    
    smoothCalcLast (zc->a0[j], a1, time, &a2, &discount);
    zc->a1[j]       = a1;
    zc->a2[j]       = a2;
    zc->prices[j+1] = zc->prices[j] * discount;
    zc->a0[j+1]     = zc->a0[j] + time * (a1 + time * a2);

    return SUCCESS;
}

/* Solves for forward rate from j to j+1 for a flat curve */
static int flatObjFunc (double x, SOLVER_CONTEXT *context,
                        double *pvDiff)
{
    int            j    = context->j;
    IrxTZeroCurve *zc   = context->zc;
    double         time = context->time;

    /* it is possible that zc and zcDisc are the same - in which
       case changing zc also changes zcDisc - despite the const status
       of zcDisc */
    zc->numItems    = j+2;
    zc->a0[j]       = x;
    zc->prices[j+1] = zc->prices[j] * exp(-x*time);

    return calcPvDiff (context, pvDiff);
}

/* Solves for first point of the smooth curve.
   Solves a0[j] using a0[j+1].
   Implies a1[j] and prices[j+1].
*/
static int smoothObjFuncFirstPoint (double x, SOLVER_CONTEXT *context,
                                    double *pvDiff)
{
    int            j    = context->j;
    IrxTZeroCurve *zc   = context->zc;
    double         time = context->time;

    double         a2;
    double         discount;

    smoothCalcFirst (x, zc->a0[j+1], time, &a2, &discount);
    zc->numItems    = j+2;
    zc->a0[j]       = x;
    zc->a1[j]       = 0;
    zc->a2[j]       = a2;
    zc->prices[j+1] = zc->prices[j] * discount;

    return calcPvDiff (context, pvDiff);
}

/* Solves for middle point of the curve.
   Solves a1[j] using a0[j], a0[j+1].
   a2[j] and prices[j+1] are implied.
 */
static int smoothObjFuncMiddlePoint (double x, SOLVER_CONTEXT *context,
                                    double *pvDiff)
{
    int            j    = context->j;
    IrxTZeroCurve *zc   = context->zc;
    double         time = context->time;

    double         a2;
    double         discount;

    /* it is possible that zc and zcDisc are the same - in which
       case changing zc also changes zcDisc - despite the const status
       of zcDisc */
    smoothCalcMiddle (zc->a0[j], zc->a0[j+1], x, time, &a2, &discount);
    zc->numItems    = j+2;
    zc->a1[j]       = x;
    zc->a2[j]       = a2;
    zc->prices[j+1] = zc->prices[j] * discount;

    return calcPvDiff (context, pvDiff);
}

/* Solves for last point of the curve.
   Solves a1[j] using a0[j].
   a2[j] and prices[j+1] are implied.
 */
static int smoothObjFuncLastPoint (double x, SOLVER_CONTEXT *context,
                                   double *pvDiff)
{
    int            j    = context->j;
    IrxTZeroCurve *zc   = context->zc;
    double         time = context->time;

    double         a2;
    double         discount;

    /* it is possible that zc and zcDisc are the same - in which
       case changing zc also changes zcDisc - despite the const status
       of zcDisc */
    smoothCalcLast (zc->a0[j], x, time, &a2, &discount);
    zc->numItems    = j+2;
    zc->a1[j]       = x;
    zc->a2[j]       = a2;
    zc->prices[j+1] = zc->prices[j] * discount;

    return calcPvDiff (context, pvDiff);
}


/* After adjusting the zc according to the rules for the particular thing
   we are trying to solve, this routine actually computes the pvDifference */
static int calcPvDiff (SOLVER_CONTEXT *context, double *pvDiff)
{
    int status = FAILURE;

    double flowsPv = 0.0;
    double fixedPv = 0.0;

    if (context->valueFloating)
    {
        IrxTCashFlowList *tmp = irxSwapPaymentsFlows (
            context->floatPayments,
            1,
            (IrxTZeroCurve const**)(&context->zc),
            0);
        if (tmp == NULL)
            goto RETURN; /* failure */

        if (irxCashFlowListPV (tmp, context->zcDisc, &flowsPv) != SUCCESS)
        {
            irxCashFlowListFree (tmp);
            goto RETURN; /* failure */
        }

        irxCashFlowListFree (tmp);
    }

    if (irxCashFlowListPV (context->fixedFlows, context->zcDisc, 
                           &fixedPv) != SUCCESS)
        goto RETURN; /* failure */

    *pvDiff = fixedPv - flowsPv - context->pv;
    status  = SUCCESS;

 RETURN:

    return status;
}


/* Note: we ignore the possible complication that the date for the currency
   basis is different from the date for the regular swap (due to different
   holiday calendars). */
static int addCashRate
(char*               routine,
 IrxTZeroCurve      *zc,
 int                 j,
 IrxTBool            smoothing,
 IrxTBool            first,
 IrxTBool            last,
 IrxTMarketConv     *marketConv,
 IrxTDate            baseDate,
 IrxTDate            cashDate,
 IrxTCalendar       *calendar,
 double              rate,
 double              cbsSpread)
{
    int status = FAILURE;

    double time;
    double payment;
    double discount;
    /* we seem to be ignoring the input price here */
    double pv = 1.0;

    IrxTDate adjDate;

    if (irxBusinessDay (cashDate,
                        marketConv->mmBdc,
                        calendar,
                        &adjDate) != SUCCESS)
        goto RETURN; /* failure */
    
    /* special case for 1D rate near the end of the month */
    if (adjDate <= baseDate)
    {
        if (irxBusinessDay (cashDate,
                            IRX_BAD_DAY_FOLLOW,
                            calendar,
                            &adjDate) != SUCCESS)
            goto RETURN; /* failure */
    }

    if (irxDayCountFraction (baseDate, adjDate, marketConv->mmDcc,
                             &time) != SUCCESS)
        goto RETURN; /* failure */

/* ALIB ZC3 assumed the same day count convention for cash currency basis */
    payment = 1.0 + (rate+cbsSpread)*time;

    if (!smoothing)
    {
        discount = pv / payment;
        zc->startDates[j+1] = adjDate;
        zc->prices[j+1]     = discount;
        if (irxDiscountToRate (zc->prices[j+1]/zc->prices[j],
                               zc->startDates[j],
                               zc->startDates[j+1],
                               IRX_ACT_365F,
                               IRX_CONTINUOUS_RATE,
                               &zc->a0[j]) != SUCCESS)
            goto RETURN; /* failure */
    }
    else
    {
        IrxTCashFlowList cfl;
        SOLVER_CONTEXT context;

        cfl.numItems = 1;
        cfl.dates    = &adjDate;
        cfl.amounts  = &payment;
        
        time = (zc->startDates[j+1]-zc->startDates[j]) / 365.0;

        context.j             = j;
        context.valueFloating = FALSE;
        context.zc            = zc;
        context.zcDisc        = zc;
        context.fixedFlows    = &cfl;
        context.floatPayments = NULL;
        context.pv            = pv;
        context.time          = time;

        if (first)
        {
            if (smoothSolverFirstPoint (rate, &context) != SUCCESS)
                goto RETURN; /* failure */
        }
        else if (last)
        {
            if (smoothSolverLastPoint (&context) != SUCCESS)
                goto RETURN; /* failure */
        }
        else
        {
            if (smoothSolverMiddlePoint (&context) != SUCCESS)
                goto RETURN; /* failure */
        }        
    }
    status = SUCCESS;

 RETURN:

    if (status != SUCCESS)
    {
        char buf[16];
        irxError ("%s: Could not add cash rate %.3f%% with date %s to curve\n",
                  routine, 1e2*rate,
                  irxDateFormat(cashDate, "DD-MMM-YYYY", buf));
    }
    return status;
}

/* Adds a swap rate (without currency basis) */
static int addSwapRate
(char*               routine,
 IrxTZeroCurve      *zc,
 IrxTZeroCurve const*zcDiscount,
 int                 j,
 IrxTBool            smoothing,
 IrxTBool            first,
 IrxTBool            last,
 IrxTMarketConv     *marketConv,
 IrxTDate            baseDate,
 IrxTDate            swapDate,
 IrxTCalendar       *calendar,
 IrxTBool            valueFloating,
 IrxTBool            firstFloatFixed,
 double              firstFloatFixedRate,
 IrxTStubLocation    stubLocation,
 double              rate,
 double              price)
{
    int status = FAILURE;

    IrxTCashFlowList *fixedFlows = NULL;
    IrxTCashFlowList *fixedBefore = NULL;
    IrxTCashFlowList *fixedAfter = NULL;
    IrxTSwapPayments *floatPayments = NULL;
    IrxTSwapPayments *floatBefore = NULL;
    IrxTSwapPayments *floatAfter = NULL;
    IrxTSwap         *swap = NULL;
    IrxTDate          maxDate;
    double            time;
    double            pv;
    double            pvBefore;
    IrxTDate          splitDate = zc->startDates[j];

    SOLVER_CONTEXT context;

    swap = irxSwapMake (marketConv,
                        rate,
                        0.0,
                        baseDate,
                        swapDate,
                        0,
                        stubLocation,
                        IRX_STUB_SIMPLE,
                        IRX_STUB_SIMPLE);

    fixedFlows = irxSwapFixedFlows (swap,
                                    1.0,
                                    baseDate,
                                    calendar,
                                    !valueFloating,
                                    FALSE);
    if (fixedFlows == NULL)
        goto RETURN; /* failure */

    maxDate = fixedFlows->dates[fixedFlows->numItems-1];

    if (irxCashFlowListSplit (fixedFlows,
                              splitDate,
                              &fixedBefore,
                              &fixedAfter) != SUCCESS)
        goto RETURN; /* failure */

    if (irxCashFlowListPV (fixedBefore, zcDiscount, &pvBefore) != SUCCESS)
        goto RETURN; /* failure */

    if (valueFloating)
    {
        int i;
        IrxTDate maxFloatDate;

        floatPayments = irxSwapPaymentsGenerate (
            swap,
            0.0,
            1.0,
            calendar,
            FALSE,
            FALSE,
            firstFloatFixed ? baseDate : 0,
            firstFloatFixed ? firstFloatFixedRate : 0.0);

        if (floatPayments == NULL)
            goto RETURN; /* failure */

        if (irxSwapPaymentsSplit (floatPayments,
                                  splitDate,
                                  &floatBefore,
                                  &floatAfter) != SUCCESS)
            goto RETURN; /* failure */

        pv = -pvBefore;
        if (floatBefore->numFlows > 0)
        {
            irxCashFlowListFree(fixedBefore);
            fixedBefore = irxSwapPaymentsFlows (
                floatBefore, 
                1,
                (IrxTZeroCurve const**)&zc,
                0);
            if (fixedBefore == NULL)
                goto RETURN; /* failure */
            
            if (irxCashFlowListPV (fixedBefore,
                                   zcDiscount,
                                   &pvBefore) != SUCCESS)
                goto RETURN; /* failure */
            pv += pvBefore;
        }
        
        for (i = 0; i <= 1; ++i)
        {
            /* i=0 for pay dates; i=1 for rate estimation dates */
            if (irxSwapPaymentsLastDate(floatAfter, i,
                                        &maxFloatDate) != SUCCESS)
                goto RETURN; /* failure */
            
            if (maxFloatDate > maxDate)
                maxDate = maxFloatDate;
        }
    }
    else
    {
        pv = price - pvBefore;
    }

    zc->startDates[j+1] = maxDate;

    time = (double)(zc->startDates[j+1] - zc->startDates[j])/365.0;

    context.zc            = zc;
    context.zcDisc        = zcDiscount;
    context.j             = j;
    context.fixedFlows    = fixedAfter;
    context.pv            = pv;
    context.valueFloating = valueFloating;
    context.floatPayments = floatAfter;
    context.time          = time;

    if (!smoothing)
    {
        if (flatSolver (rate, &context) != SUCCESS)
            goto RETURN; /* failure */
    }
    else if (first)
    {
        if (smoothSolverFirstPoint (rate, &context) != SUCCESS)
            goto RETURN; /* failure */
    } 
    else if (last)
    {
        if (smoothSolverLastPoint (&context) != SUCCESS)
            goto RETURN; /* failure */
    }
    else
    {
        if (smoothSolverMiddlePoint (&context) != SUCCESS)
            goto RETURN; /* failure */
    }        
    
    status = SUCCESS;
    
 RETURN:
    
    irxCashFlowListFree (fixedFlows);
    irxCashFlowListFree (fixedBefore);
    irxCashFlowListFree (fixedAfter);
    irxSwapPaymentsFree (floatPayments);
    irxSwapPaymentsFree (floatBefore);
    irxSwapPaymentsFree (floatAfter);
    irxSwapFree (swap);
    
    return status;
}

/* Adds a swap rate (with currency basis) to the discount curve */
static int addSwapRateCB
(char*               routine,
 IrxTZeroCurve      *zcDisc,
 int                 j,
 IrxTBool            smoothing,
 IrxTBool            first,
 IrxTBool            last,
 IrxTMarketConv     *marketConv,
 IrxTDate            baseDate,
 IrxTDate            swapDate,
 IrxTCalendar       *calendar,
 IrxTCalendar       *cbsCalendar,
 IrxTStubLocation    stubLocation,
 double              rate,
 double              cbsSpread,
 double              price)
{
    int status = FAILURE;

    IrxTCashFlowList *fixedFlows = NULL;
    IrxTCashFlowList *basisFlows = NULL;
    IrxTCashFlowList *combinedFlows = NULL;
    IrxTCashFlowList *fixedBefore = NULL;
    IrxTCashFlowList *fixedAfter = NULL;
    IrxTSwap         *swap = NULL;
    IrxTSwap         *basisSwap = NULL;
    IrxTMarketConv    basisMktConv;
    IrxTDate          maxDate;
    double            time;
    double            pvBefore;
    IrxTDate          splitDate = zcDisc->startDates[j];

    SOLVER_CONTEXT context;

    /* Adjust the basisMktConv so that the basis conventions are correct */
    basisMktConv = *marketConv;
    basisMktConv.fixedIvl = marketConv->cbsIvl;
    basisMktConv.fixedDcc = marketConv->cbsDcc;

    swap = irxSwapMake (marketConv,
                        rate,
                        0.0,
                        baseDate,
                        swapDate,
                        0,
                        stubLocation,
                        IRX_STUB_SIMPLE,
                        IRX_STUB_SIMPLE);

    basisSwap = irxSwapMake (&basisMktConv,
                             cbsSpread,
                             0.0,
                             baseDate,
                             swapDate,
                             0,
                             stubLocation,
                             IRX_STUB_SIMPLE,
                             IRX_STUB_SIMPLE);

    fixedFlows = irxSwapFixedFlows (swap,
                                    1.0,
                                    baseDate,
                                    calendar,
                                    TRUE,
                                    FALSE);
    if (fixedFlows == NULL)
        goto RETURN; /* failure */

    basisFlows = irxSwapFixedFlows (basisSwap,
                                    1.0,
                                    baseDate,
                                    cbsCalendar,
                                    FALSE,
                                    FALSE);
    if (basisFlows == NULL)
        goto RETURN; /* failure */

    combinedFlows = irxCashFlowListMerge (fixedFlows, basisFlows);
    if (combinedFlows == NULL)
        goto RETURN; /* failure */

    maxDate = combinedFlows->dates[combinedFlows->numItems-1];

    if (irxCashFlowListSplit (combinedFlows,
                              splitDate,
                              &fixedBefore,
                              &fixedAfter) != SUCCESS)
        goto RETURN; /* failure */

    if (irxCashFlowListPV (fixedBefore, zcDisc, &pvBefore) != SUCCESS)
        goto RETURN; /* failure */

    zcDisc->startDates[j+1] = maxDate;

    time = (double)(zcDisc->startDates[j+1] - zcDisc->startDates[j])/365.0;

    context.zc            = zcDisc;
    context.zcDisc        = zcDisc;
    context.j             = j;
    context.fixedFlows    = fixedAfter;
    context.pv            = price - pvBefore;
    context.valueFloating = FALSE;
    context.floatPayments = NULL;
    context.time          = time;

    if (!smoothing)
    {
        if (flatSolver (rate, &context) != SUCCESS)
            goto RETURN; /* failure */
    }
    else if (first)
    {
        if (smoothSolverFirstPoint (rate, &context) != SUCCESS)
            goto RETURN; /* failure */
    } 
    else if (last)
    {
        if (smoothSolverLastPoint (&context) != SUCCESS)
            goto RETURN; /* failure */
    }
    else
    {
        if (smoothSolverMiddlePoint (&context) != SUCCESS)
            goto RETURN; /* failure */
    }        
    status = SUCCESS;

 RETURN:

    irxCashFlowListFree (fixedAfter);
    irxCashFlowListFree (fixedBefore);
    irxCashFlowListFree (combinedFlows);
    irxCashFlowListFree (basisFlows);
    irxCashFlowListFree (fixedFlows);
    irxSwapFree (basisSwap);
    irxSwapFree (swap);
    
    return status;
}

#if 0
/* Adds a swap rate (with currency basis) to the index curve */
static int addSwapRateCBIndex
(char*               routine,
 IrxTZeroCurve      *zcIndex,
 IrxTZeroCurve const*zcDisc,
 int                 j,
 IrxTBool            smoothing,
 IrxTBool            first,
 IrxTBool            last,
 IrxTMarketConv     *marketConv,
 IrxTDate            baseDate,
 IrxTDate            swapDate,
 IrxTCalendar       *calendar,
 IrxTBool            firstFloatFixed,
 double              firstFloatFixedRate,
 IrxTStubLocation    stubLocation,
 double              rate)
{
    int status = FAILURE;

    IrxTCashFlowList *fixedFlows = NULL;
    IrxTSwapPayments *floatPayments = NULL;
    IrxTSwap         *swap = NULL;
    IrxTDate          maxDate;
    double            fwdRate;
    double            time;
    int               i;
    IrxTDate          maxFloatDate;

    SOLVER_CONTEXT context;

    swap = irxSwapMake (marketConv,
                        rate,
                        0.0,
                        baseDate,
                        swapDate,
                        0,
                        stubLocation,
                        IRX_STUB_SIMPLE,
                        IRX_STUB_SIMPLE);

    fixedFlows = irxSwapFixedFlows (swap,
                                    1.0,
                                    baseDate,
                                    calendar,
                                    FALSE,
                                    FALSE);
    if (fixedFlows == NULL)
        goto RETURN; /* failure */

    maxDate = fixedFlows->dates[fixedFlows->numItems-1];

    floatPayments = irxSwapPaymentsGenerate (
        swap,
        0.0,
        1.0,
        calendar,
        FALSE,
        FALSE,
        firstFloatFixed ? baseDate : 0,
        firstFloatFixed ? firstFloatFixedRate : 0.0);
    
    if (floatPayments == NULL)
        goto RETURN; /* failure */
    
    for (i = 0; i <= 1; ++i)
    {
        /* i=0 for pay dates; i=1 for rate estimation dates */
        if (irxSwapPaymentsLastDate(floatPayments, i,
                                    &maxFloatDate) != SUCCESS)
            goto RETURN; /* failure */
        
        if (maxFloatDate > maxDate)
            maxDate = maxFloatDate;
    }

    zcIndex->startDates[j+1] = maxDate;

    time = (double)(zcIndex->startDates[j+1] - zcIndex->startDates[j])/365.0;

    context.zc            = zcIndex;
    context.zcDisc        = zcDisc;
    context.j             = j;
    context.fixedFlows    = fixedFlows;
    context.pv            = 0.0;
    context.valueFloating = TRUE;
    context.floatPayments = floatPayments;
    context.time          = time;

    if (irxRootFindBrent ((IrxTObjectFunc)flatObjFunc,
                          &context,
                          LOWER_BOUND,
                          UPPER_BOUND,
                          MAX_ITERATIONS,
                          rate, /* a not very good guess */
                          INITIAL_X_STEP,
                          INITIAL_F_DERIV,
                          X_TOLERANCE,
                          F_TOLERANCE,
                          &fwdRate) != SUCCESS)
        goto RETURN; /* failure */

    zcIndex->numItems    = j+1;
    zcIndex->a0[j]       = fwdRate;
    zcIndex->prices[j+1] = zcIndex->prices[j] * exp(-fwdRate*time);

    status = SUCCESS;

 RETURN:

    irxCashFlowListFree (fixedFlows);
    irxSwapPaymentsFree (floatPayments);
    irxSwapFree (swap);

    return status;
}
#endif

/* For the first point of the curve, computes a1[j] and prices[j+1] given
   a0[j], a0[j+1] */
static void smoothCalcFirst(double a0, double a0next, double time,
                            double *a2, double *discount)
{
    *a2 = (a0next - a0) / (time*time);
    *discount = exp(-time*(a0 + (*a2) * time * time / 3.0));
}

/* For middle point of the curve, computes a2[j] and prices[j+1] given
   a0[j], a0[j+1], a1[j] */
static void smoothCalcMiddle(double a0, double a0next, double a1, double time,
                             double *a2, double *discount)
{
    *a2 = (a0next - a0 - time * a1) / (time*time);
    *discount = exp(-time*(a0 + time * (a1/2.0 + (*a2) * time/3.0)));
}

/* For the last point of the curve, computes a2[j] and prices[j+1] given
   a0[j], a1[j] */
static void smoothCalcLast(double a0, double a1, double time,
                           double *a2, double *discount)
{
    *a2 = -a1 / (2.0*time); /* flat at t=time */
    *discount = exp(-time*(a0 + time * (a1/2.0 + (*a2) * time/3.0)));
}



