/*
***************************************************************************
** SOURCE FILE: cdscurvedateadj.c
**
** Implementation file for function which converts market spreads of a CDS
** curve into spreads on a different set of dates such that when you
** bootstrap from the new set of dates you back-out the results on the
** market dates.
***************************************************************************
*/

#include "cdscurvedateadj.h"

#include "cdsbootstrap.h"
#include "cds.h"
#include "recovery.h"

#include <alib/mgpimsl.h>
#include <alib/cmemory.h>
#include <alib/ldate.h>
#include <alib/tcurve.h>

#include <cxutils/include/date.h>
#include <cxutils/include/cxmacros.h>
#include <cxutils/include/zerocurve.h>

typedef struct _SOLVER_DATA
{
    TDate           today;
    TCurve         *discCurve;
    TDate           startDate;
    TDate           valueDate;
    int             nbDate;
    TDate          *endDates;
    double         *couponRates;
    double         *prices;
    TDate          *adjEndDates;
    double         *recalcSpreads;
    CxTRecoveryCurve *recoveryCurve;
    TBoolean        payAccOnDefault;
    TDateInterval  *couponInterval;
    CxTDayCountConv  paymentDcc;
    CxTStubType      stubType;
    CxTCreditCurveType curveType;
    TDateInterval  *timestep;
    TBoolean        protectStart;
    TBoolean        isPriceClean;
    long            delay;
    CxTBadDayConv   badDayConv;
    CxTCalendar    *calendar;
    TBoolean        optimize;
    int             status;
} SOLVER_DATA;

static SOLVER_DATA global;

static TMatrix2D* computeTweakMatrix
(int     nbDate,       /* (I) */
 double *spotSpreads); /* (I) */

static void matchMarketSpreads
(int     n,             /* (I) */
 double *adjSpreads,    /* (I) */ 
 double *marketErrors); /* (O) */

/**
 * This function takes a curve of CDS rates at the market maturity dates
 * and creates a curve of CDS rates at slightly different dates which
 * might typically be used by risk management systems.
 *
 * When you bootstrap using the resulting curve of CDS rates then you
 * will obtain the CDS rates for the market maturity dates.
 *
 * Typically this is required due to the market convention of using ISDA
 * dates (20th of the months of March, June, September, December) whereas
 * some systems assume the CDSs have maturity exactly N years after the
 * start date of the deal.
 *
 * In addition we can return on request a Jacobian matrix (dAdj/dMkt) to
 * enable us to convert positions returned by the risk management system.
 *
 * Risk starts at the end of today. Where relevant, the PV is computed
 * for value date. The CDS's all start at a single startDate, and
 * end at the given endDates. The last date is always protected - the start
 * date is only protected if protectStart=True.
 *
 * Interest accrues for the same number of days as there is protection.
 * Thus if protectStart=True you get one extra day of accrued interest in
 * comparison with an interest rate swap. This extra day is assumed to be
 * the last day of the CDS and means that the last period is one day longer
 * than for an interest rate swap.
 *
 * Maximum compatibility with CMLib is obtained by setting protectStart=False
 * and curveType=CX_CURVE_TYPE_FLOW.
 */
CxTCdsCurveDateAdj* CxCdsCurveDateAdjustment(
    TDate           today,
    TCurve         *discCurve,
    TDate           startDate,
    TDate           valueDate,
    long            nbDate,
    TDate          *endDates,
    double         *couponRates,
    double         *prices, 
    TBoolean       *includes,
    TDate          *adjEndDates,
    CxTRecoveryCurve *recoveryCurve,
    TBoolean        payAccOnDefault,
    TDateInterval  *couponInterval,
    CxTDayCountConv paymentDcc,
    CxTStubType     stubType,
    CxTCreditCurveType curveType,
    TDateInterval  *timestep,   
    TBoolean        protectStart,
    TBoolean        isPriceClean,
    long            delay,
    CxTBadDayConv   badDayConv,
    CxTCalendar    *calendar,
    TBoolean        computeTweakMatrixFlag
)
{
    static char routine[] = "CxCdsCurveDateAdjustment";
    int         status    = FAILURE;

    CxTCdsCurveDateAdj *curveDateAdj = NULL;
    double             *solveSpreads = NULL;
    double             *guessSpreads = NULL;
    double             *recalcSpreads = NULL;

    TDate              *includeEndDates = NULL;
    TDate              *includeAdjEndDates = NULL;
    double             *includeCouponRates = NULL;
    double             *includePrices = NULL;

    TBoolean            optimize = FALSE;

    int                i;
    CxTCreditCurve     *marketCreditCurve = NULL;
    CxTCreditCurve     *adjCreditCurve = NULL;
    //CxTRecoveryCurve   *recoveryCurve = NULL;

    char dateBuf[16]; /* used for date formatting */

    REQUIRE (discCurve != NULL);
    REQUIRE (startDate >= today);
    REQUIRE (valueDate >= today);
    REQUIRE (nbDate > 0);
    REQUIRE (endDates != NULL);
    REQUIRE (couponRates != NULL);
    REQUIRE (adjEndDates != NULL);
    //REQUIRE (recoveryRate >= 0.0 && recoveryRate <= 1.0);

    //recoveryCurve = CxRecoveryCurveMakeFromRecoveryRate(recoveryRate);
    REQUIRE(recoveryCurve !=NULL);
    // if (recoveryCurve == NULL)
    //  goto done; /* failure */

    /*
     * 1. Bootstrap the Cds curve using the market rates.
     * 2. Imply the adjSpreads from this curve and use this as the initial
     *    guess.
     * 3. Solve for adjSpreads such that the curve bootstrapped from this
     *    will return the market rates.
     */
    
    marketCreditCurve = CxCdsBootstrap (today,
                                        discCurve,
                                        startDate,
                                        valueDate,
                                        nbDate,
                                        endDates,
                                        couponRates,
                                        prices,
                                        includes,
                                        recoveryCurve,
                                        payAccOnDefault,
                                        couponInterval,
                                        paymentDcc,
                                        stubType,
                                        curveType,
                                        timestep,
                                        NULL,
                                        protectStart,
                                        isPriceClean,
                                        delay,
                                        badDayConv,
                                        calendar);
    if (marketCreditCurve == NULL)
    {
        GtoErrMsg ("%s: Could not bootstrap market spread curve\n", routine);
        goto done; /* failure */
    }

    optimize = (curveType == CX_CURVE_TYPE_FLOW);

    /* when we are in the solver we only want to tweak the spreads which
       have an effect */
    if (includes != NULL)
    {
        int nbIncludes = 0;
        int j = 0;
        for (i = 0; i < nbDate; ++i)
        {
            if (includes[i]) ++nbIncludes;
        }
        if (nbIncludes == 0)
        {
            GtoErrMsg ("%s: No dates included\n", routine);
            goto done; /* failure */
        }
        includeEndDates    = NEW_ARRAY(TDate, nbIncludes);
        includeCouponRates = NEW_ARRAY(double, nbIncludes);
        if (prices != NULL)
            includePrices  = NEW_ARRAY(double, nbIncludes);
        includeAdjEndDates = NEW_ARRAY(TDate, nbIncludes);
        j = 0;
        for (i = 0; i < nbDate; ++i)
        {
            if (includes[i])
            {
                includeEndDates[j]    = endDates[i];
                includeCouponRates[j] = couponRates[i];
                if (prices != NULL)
                    includePrices[j]  = prices[i];
                includeAdjEndDates[j] = adjEndDates[i];
                ++j;
            }
        }
        nbDate      = nbIncludes;
        endDates    = includeEndDates;
        couponRates = includeCouponRates;
        prices      = includePrices;
        adjEndDates = includeAdjEndDates;
    }

    guessSpreads  = NEW_ARRAY(double, nbDate);
    recalcSpreads = NEW_ARRAY(double, nbDate);

    if (optimize)
    {
        /* Instead of going through the par spreads we will use the clean
           spreads internally. This should enable us to skip a whole bunch
           of bootstrap routines. The guess spreads will be expressed as
           continuously compounded rates.

           This relies on knowing that the dates for the clean spread curve
           are the same dates as the adjusted spread curve. */

        for (i = 0; i < nbDate; ++i)
        {
            guessSpreads[i] = CxZeroRate (marketCreditCurve->tc,
                                          adjEndDates[i]);
        }
    }
    else
    {
        for (i = 0; i < nbDate; ++i)
        {
            if (CxCdsParSpread (today,
                                valueDate,
                                startDate,
                                adjEndDates[i],
                                delay, /* delay */
                                0.0, /* price */
                                payAccOnDefault,
                                couponInterval,
                                stubType,
                                paymentDcc,
                                badDayConv,
                                calendar,
                                discCurve,
                                marketCreditCurve,
                                recoveryCurve,
                                protectStart,
                                isPriceClean,
                                &guessSpreads[i]) != SUCCESS)
            {
                GtoErrMsg ("%s: Could not compute initial adjusted spread for "
                           "%s\n", routine,
                           CxDateFormat (adjEndDates[i], 
                                         "YYYY-MM-DD", 
                                         dateBuf));
                goto done; /* failure */
            }
        }
    }

    global.today            = today;
    global.discCurve        = discCurve;
    global.startDate        = startDate;
    global.valueDate        = valueDate;
    global.nbDate           = nbDate;
    global.endDates         = endDates;
    global.couponRates      = couponRates;
    global.prices           = prices;
    global.adjEndDates      = adjEndDates;
    global.recalcSpreads    = recalcSpreads;
    global.recoveryCurve    = recoveryCurve;
    global.payAccOnDefault  = payAccOnDefault;
    global.couponInterval   = couponInterval;
    global.paymentDcc       = paymentDcc;
    global.stubType         = stubType;
    global.curveType        = curveType;
    global.timestep         = timestep;
    global.protectStart     = protectStart;
    global.delay            = delay;
    global.isPriceClean     = isPriceClean;
    global.badDayConv       = badDayConv;
    global.calendar         = calendar;
    global.optimize         = optimize;
    global.status           = SUCCESS;

/*     GtoImslSetErrorHandling(); */
    if (GtoImslZerosSysEqn (matchMarketSpreads,
                            nbDate,
                            guessSpreads,
                            1.0e-8, /* errBound */
                            100,    /* maxIter */
                            &solveSpreads) != SUCCESS ||
        global.status != SUCCESS)
    {
        GtoErrMsg ("%s: Could not solve for adjusted dates\n", routine);
        goto done; /* failure */
    }

    if (optimize)
    {
        /* need to convert solveSpreads from clean spreads to par spreads */
        TCurve *tmp = NULL;

        tmp = GtoMakeTCurve (today,
                             adjEndDates,
                             solveSpreads,
                             nbDate,
                             GTO_CONTINUOUS_BASIS,
                             GTO_ACT_365F);
        if (tmp == NULL)
            goto done; /* failure */

        adjCreditCurve = CxCreditCurveMake(curveType,
                                           tmp,
                                           timestep);
        GtoFreeTCurve (tmp);
        if (adjCreditCurve == NULL)
            goto done; /* failure */

        for (i = 0; i < nbDate; ++i)
        {
            if (CxCdsParSpread (today,
                                valueDate,
                                startDate,
                                adjEndDates[i],
                                delay, /* delay */
                                0.0, /* price */
                                payAccOnDefault,
                                couponInterval,
                                stubType,
                                paymentDcc,
                                badDayConv,
                                calendar,
                                discCurve,
                                adjCreditCurve,
                                recoveryCurve,
                                protectStart,
                                isPriceClean,
                                &solveSpreads[i]) != SUCCESS)
            {
                GtoErrMsg ("%s: Could not compute adjusted spread for %s\n",
                           routine,
                           CxDateFormat (adjEndDates[i], 
                                         "YYYY-MM-DD", 
                                         dateBuf));
                goto done; /* failure */
            }
        }
    }

    /* sanity checks on solve rates */
    for (i = 0; i < nbDate; ++i)
    {
        if (solveSpreads[i] < 0.0)
        {
            GtoErrMsg ("%s: Negative spread %f%% for date %s\n",
                       routine, solveSpreads[i]*100.0,
                       CxDateFormat (adjEndDates[i], "YYYY-MM-DD", dateBuf));
            goto done; /* failure */
        }
        /* might re-check market spreads */
    }

    curveDateAdj = CxCdsCurveDateAdjMake (nbDate, 
                                          adjEndDates,
                                          solveSpreads,
                                          NULL);
    if (curveDateAdj == NULL) 
        goto done; /* failure */

    /* Compute tweak matrix */
    if (computeTweakMatrixFlag && prices == NULL)
    {
        curveDateAdj->tweakMatrix = computeTweakMatrix (nbDate, solveSpreads);
        if (curveDateAdj->tweakMatrix == NULL)
            goto done; /* failure */
    }

    status = SUCCESS;

 done:

    if (status != SUCCESS)
    {
        CxCdsCurveDateAdjFree (curveDateAdj);
        curveDateAdj = NULL;
        GtoErrMsgFailure (routine);
    }

    CxCreditCurveFree (marketCreditCurve);
    CxCreditCurveFree (adjCreditCurve);
    FREE (includeAdjEndDates);
    FREE (includeEndDates);
    FREE (includePrices);
    FREE (includeCouponRates);
    FREE (guessSpreads);
    FREE (recalcSpreads);
    GtoFreeSafe (solveSpreads); /* allocated within ALIB */

    return curveDateAdj;
}

/**
 * First tweak spot rates to get dMarket/dAdj, then invert the matrix using 
 * IMSL.
 */
static TMatrix2D* computeTweakMatrix
(int     nbDate,      /* (I) */
 double *adjSpreads)  /* (I) */
{
#define TWEAK_SIZE 1e-5

    static char routine[] = "computeTweakMatrix";
    int         status    = FAILURE;

    double *dmat = NULL;
    double *mktChanges = NULL;
    double *adjTweaks  = NULL;
    double *inv = NULL;
    TMatrix2D *tweakMatrix = NULL;
    int size = nbDate*nbDate;

    int i,j;
    
    dmat       = NEW_ARRAY(double, size);
    adjTweaks  = NEW_ARRAY(double, nbDate);
    mktChanges = NEW_ARRAY(double, nbDate);
    if(dmat==NULL || adjTweaks==NULL || mktChanges==NULL)
        goto done;

    COPY_ARRAY (adjTweaks, adjSpreads, double, nbDate);

    /* since we are using matchMarketSpreads we need to set global.optimize
       to FALSE otherwise we won't understand the type of spreads passed in */
    global.optimize = FALSE;

    /* j-th column are the tweaks of all fwds w/respect to j-th spot
     */
    for (j=0; j<nbDate; j++)
    {
        /* tweak up */
        adjTweaks[j] = adjSpreads[j] + TWEAK_SIZE;
        matchMarketSpreads(nbDate,
                           adjTweaks,
                           mktChanges);

        if(global.status == FAILURE)
            goto done;

        for (i=0; i<nbDate; i++) /* fill in j-th column */
        {
            /* at this step we just store dM for tweak up */
            dmat[j+i*nbDate] = mktChanges[i];
        }

        /* tweak down */
        adjTweaks[j] = adjSpreads[j] - TWEAK_SIZE;
        matchMarketSpreads(nbDate,
                           adjTweaks,
                           mktChanges);
        if(global.status == FAILURE)
            goto done;
        
        for (i=0; i<nbDate; i++) /* fill in j-th column */
        {
            int k = j + i*nbDate;

            /* we have dM for tweak up stored in dmat and dM for tweak
               down calculated here - so we now calculate dmat */
            
            dmat[k] = (dmat[k] - mktChanges[i]) / (2.0 * TWEAK_SIZE);
        }
        
        /* no tweak - put adjTweaks back to original value */
        adjTweaks[j] = adjSpreads[j];
    }
    
    /* Invert matrix */
    if (GtoImslInvertMatrix(nbDate, dmat, &inv) ISNT SUCCESS )
        goto done; /* failure */

    /* convert the inverted matrix into a TMatrix2D structure */
    tweakMatrix = GtoMatrix1DTo2DNew (nbDate, nbDate, inv);
    if (tweakMatrix == NULL)
        goto done; /* failure */

    status = SUCCESS;

 done: 

    FREE (adjTweaks);
    FREE (mktChanges);
    FREE (dmat);
    GtoFreeSafe (inv); /* allocated within ALIB */
                          
    if (status != SUCCESS)
    {
        GtoErrMsgFailure (routine);
        GtoMatrixFree (tweakMatrix);
        tweakMatrix = NULL;
    }

    return tweakMatrix;
}



static void matchMarketSpreads
(int     n,             /* (I) */
 double *adjSpreads,    /* (I) */ 
 double *marketErrors)  /* (O) */
{
    static char routine[] = "matchMarketSpreads";
    int         status    = FAILURE;

    int                i;
    CxTCreditCurve     *adjCreditCurve = NULL;
    char dateBuf[16]; /* used for date formatting */

    for (i = 0; i < n; ++i)
    {
        if (adjSpreads[i] < 0.0)
        {
            GtoErrMsg ("%s: Root solver implies negative spread (%.2f%%) for "
                       "date %s\n",
                       routine, 1e2 * adjSpreads[i],
                       CxDateFormat (global.adjEndDates[i], "YYYY-MM-DD",
                                     dateBuf));
            goto done; /* failure */
        }
    }

    if (global.optimize)
    {
        TCurve *tmp = NULL;

        tmp = GtoMakeTCurve (global.today,
                             global.adjEndDates,
                             adjSpreads,
                             global.nbDate,
                             GTO_CONTINUOUS_BASIS,
                             GTO_ACT_365F);
        if (tmp == NULL)
            goto done; /* failure */

        adjCreditCurve = CxCreditCurveMake(global.curveType,
                                           tmp,
                                           global.timestep);
        GtoFreeTCurve (tmp);
        if (adjCreditCurve == NULL)
            goto done; /* failure */
    }
    else
    {
        adjCreditCurve = CxCdsBootstrap (global.today,
                                         global.discCurve,
                                         global.startDate,
                                         global.valueDate,
                                         global.nbDate,
                                         global.adjEndDates,
                                         adjSpreads,
                                         NULL,
                                         NULL,
                                         global.recoveryCurve,
                                         global.payAccOnDefault,
                                         global.couponInterval,
                                         global.paymentDcc,
                                         global.stubType,
                                         global.curveType,
                                         global.timestep,
                                         NULL,
                                         global.protectStart,
                                         global.isPriceClean,
                                         global.delay,
                                         global.badDayConv,
                                         global.calendar);
        if (adjCreditCurve == NULL)
        {
            GtoErrMsg ("%s: Could not bootstrap adjusted spread curve\n",
                       routine);
            goto done; /* failure */
        }
    }
        
    for (i = 0; i < global.nbDate; ++i)
    {
        double price;
        
        if (global.prices != NULL)
            price = global.prices[i];
        else
            price = 0.0;

        if (CxCdsParSpread (global.today,
                            global.valueDate,
                            global.startDate,
                            global.endDates[i],
                            global.delay,
                            price,
                            global.payAccOnDefault,
                            global.couponInterval,
                            global.stubType,
                            global.paymentDcc,
                            global.badDayConv,
                            global.calendar,
                            global.discCurve,
                            adjCreditCurve,
                            global.recoveryCurve,
                            global.protectStart,
                            global.isPriceClean,
                            &global.recalcSpreads[i]) != SUCCESS)
        {
            GtoErrMsg ("%s: Could not recalculate market spread for %s\n",
                       routine,
                       CxDateFormat (global.endDates[i], "YYYY-MM-DD",
                                     dateBuf));
            goto done; /* failure */
        }
        marketErrors[i] = global.recalcSpreads[i] - global.couponRates[i];
    }

    status = SUCCESS;

 done:

    if (status != SUCCESS)
    {
        /* setting marketErrors to zero will end the iteration */
        /* setting global.status to FAILURE will tell the outside routine
           that there was a failure in here */
        global.status = FAILURE;
        for (i = 0; i < global.nbDate; ++i)
            marketErrors[i] = 0.0;
        GtoErrMsgFailure (routine);
    }

    CxCreditCurveFree (adjCreditCurve);
}
    

