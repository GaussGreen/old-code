/*
***************************************************************************
** SOURCE FILE: cdsbootstrap.c
**
** CDS bootstrap routines.
**
** $Header$
***************************************************************************
*/

#include "cdsbootstrap.h"

#include "cds.h"
#include "contingentleg.h"
#include "creditcurve.h"
#include "feeleg.h"
#include "recovery.h"

#include <cxutils/include/strutils.h>
#include <cxutils/include/zerocurve.h>
#include <cxutils/include/cxmacros.h>
#include <cxutils/include/datelist.h>
#include <cxutils/include/date.h>
#include <crxflow/include/crcrv.h>
#include <alib/mgpimsl.h>
#include <alib/rtbrent.h>

typedef struct
{
    int              i;
    TCurve          *discountCurve;
    CxTCreditCurve   *cdsCurve;
    CxTRecoveryCurve *recoveryCurve;
    double           pvPrice;
    TBoolean         isPriceClean;
    CxTContingentLeg *cl;
    CxTFeeLeg        *fl;
} CDS_BOOTSTRAP_CONTEXT;


/* static function declarations */
static int cdsBootstrapPointFunction
(double   cleanSpread,
 void    *data,
 double  *pv);

static CxTCreditCurve* CdsBootstrapFlowType
(TDate           today,           /* (I) Mostly ignored                     */
 TCurve         *discountCurve,   /* (I) Risk-free discount curve           */
 TDate           startDate,       /* (I) Start of CDS for accrual and risk  */
 TDate           valueDate,       /* (I) Settlement date for price          */
 long            nbDate,          /* (I) Number of benchmark dates          */
 TDate          *endDates,        /* (I) Maturity dates of CDS to bootstrap */
 double         *couponRates,     /* (I) Coupons (e.g. 0.05 = 5% = 500bp)   */ 
 double         *prices,          /* (I) Prices (a.k.a. upfront charge)
                                     Can be NULL if no upfront charges.     */
 CxTRecoveryCurve *recoveryCurve, /* (I) Recovery curve                     */
 TBoolean        payAccOnDefault, /* (I) Pay accrued on default             */
 TDateInterval  *couponInterval,  /* (I) Interval between fee payments      */
 CxTDayCountConv paymentDCC,      /* (I) DCC for fee payments and accrual   */
 CxTStubType     stubType,        /* (I) Stub type for fee leg              */
 TDateInterval  *timestep,        /* (I) Integration timestep .             */
 TDateInterval  *smoothInterval,  /* (I) Smoothing interval. NULL if no
                                     smoothing.                             */
 TBoolean        protectStart,    /* (I) Protected on both start and end
                                     date of the CDS */
 TBoolean        isPriceClean,  /* (I) Price is quoted as clean           */
 long            delay,
 CxTBadDayConv   badDayConv,
 CxTCalendar    *calendar);

static CxTCreditCurve* CdsBootstrapExoticType
(TDate           today,           /* (I) Mostly ignored                     */
 TCurve         *discountCurve,   /* (I) Risk-free discount curve           */
 TDate           startDate,       /* (I) Start of CDS for accrual and risk  */
 TDate           valueDate,       /* (I) Settlement date for price          */
 long            nbDate,          /* (I) Number of benchmark dates          */
 TDate          *endDates,        /* (I) Maturity dates of CDS to bootstrap */
 double         *couponRates,         /* (I) CouponRates (e.g. 0.05 = 5% = 500bp)   */ 
 double         *prices,          /* (I) Prices (a.k.a. upfront charge)
                                     Can be NULL if no upfront charges.     */
 CxTRecoveryCurve *recoveryCurve, /* (I) Recovery curve                     */
 TBoolean        payAccOnDefault, /* (I) Pay accrued on default        */
 TDateInterval  *couponInterval,  /* (I) Interval between fee payments      */
 CxTDayCountConv paymentDCC,      /* (I) DCC for fee payments and accrual   */
 CxTStubType     stubType,        /* (I) Stub type for fee leg              */
 TDateInterval  *timestep,        /* (I) Integration timestep .             */
 TDateInterval  *smoothInterval,  /* (I) Smoothing interval. NULL if no
                                     smoothing.                             */
 TBoolean        protectStart,    /* (I) Protected on both start and end
                                         date of the CDS */
 TBoolean        isPriceClean,  /* (I) Price is quoted as clean           */
 long            delay,
 CxTBadDayConv   badDayConv,
 CxTCalendar    *calendar);

/*
***************************************************************************
** We use an IMSL routine for multi-dimensional solving. This does not
** support the concept of user defined data, so we have to pass the data
** via a file global variable.
**
** This is the SOLVER_DATA global defined below.
***************************************************************************
*/
typedef struct _SOLVER_DATA
{
    TDate           today;
    TCurve         *discountCurve;
    TDate           startDate;
    TDate           valueDate;
    int             nbDate;
    TDate          *endDates;
    double         *couponRates;
    double         *prices;
    CxTRecoveryCurve *recoveryCurve;
    TBoolean        payAccOnDefault;
    TDateInterval  *couponInterval;
    CxTDayCountConv paymentDCC;
    CxTStubType     stubType;
    int             status;
    TDateList      *allDates;
    TBoolean        protectStart;
    TBoolean        isPriceClean;
    long            delay;
    CxTBadDayConv   badDayConv;
    CxTCalendar    *calendar;
} SOLVER_DATA;

static SOLVER_DATA global;

static void matchBenchmarkPrices
(int     n,             /* (I) */
 double *cleanSpreads,  /* (I) */ 
 double *priceDiffs);   /* (O) */

static CxTCreditCurve* CxCdsBootstrapSmoothing
(CxTCreditCurve *cdsCurve,      /* (I) unsmoothed curve */
 TCurve         *discountCurve,
 TDate           startDate,
 TDate           valueDate,
 long            nbDate,
 TDate          *endDates,
 double         *couponRates,
 double         *prices,
 CxTRecoveryCurve *recoveryCurve,
 TBoolean        payAccOnDefault,
 TDateInterval  *couponInterval,
 CxTDayCountConv paymentDCC,
 CxTStubType     stubType,
 TDateInterval  *smoothInterval,
 TBoolean        protectStart,
 TBoolean        isPriceClean,
 long            delay,
 CxTBadDayConv   badDayConv,
 CxTCalendar    *calendar);


/*f
***************************************************************************
** Performs a smoothing algorithm to provide relatively continuous forwards
** for a credit curve.
**
** Probably not a complete solution. The mechanism has some problems:
** 1. The end condition is not quite right. It suggests the curve goes 
**    flat at the last date of the curve which does not seem correct.
** 2. The multi-dimensional solver introduces potential hedge leaks between
**    different parts of the curve.
**
** There is an alternative but this method was so easy to implement that
** we did this one first and the alternative remains on the back burner.
***************************************************************************
*/
static CxTCreditCurve* CxCdsBootstrapSmoothing
(CxTCreditCurve *cdsCurve,      /* (I) unsmoothed curve */
 TCurve         *discountCurve,
 TDate           startDate,
 TDate           valueDate,
 long            nbDate,
 TDate          *endDates,
 double         *couponRates,
 double         *prices,
 CxTRecoveryCurve *recoveryCurve,
 TBoolean        payAccOnDefault,
 TDateInterval  *couponInterval,
 CxTDayCountConv paymentDCC,
 CxTStubType     stubType,
 TDateInterval  *smoothInterval,
 TBoolean        protectStart,
 TBoolean        isPriceClean,
 long            delay,
 CxTBadDayConv   badDayConv,
 CxTCalendar    *calendar)
{
    static char routine[] = "CxCdsBootstrapSmoothing";
    int         status    = FAILURE;

    TCurve        *sparseCurve = NULL;
    CxTCreditCurve *smoothCurve = NULL;
    double        *cleanSpreads = NULL;
    TDateList     *allDates = NULL;
    TDate         *zcDates = NULL;
    double        *zcRates = NULL;

    REQUIRE (cdsCurve->tc->fNumItems == nbDate);
    REQUIRE (cdsCurve->type == CX_CURVE_TYPE_FLOW);

    zcDates  = GtoDatesFromCurve (cdsCurve->tc);
    zcRates  = GtoRatesFromCurve (cdsCurve->tc);
    allDates = CxDateListMakeRegular (cdsCurve->tc->fBaseDate,
                                      cdsCurve->tc->fArray[nbDate-1].fDate,
                                      smoothInterval,
                                      CX_LONG_FRONT_STUB);
    if (allDates == NULL)
        goto done; /* failure */

    global.today         = cdsCurve->tc->fBaseDate;
    global.discountCurve = discountCurve;
    global.startDate     = startDate;
    global.valueDate     = valueDate;
    global.nbDate        = nbDate;
    global.endDates      = endDates;
    global.couponRates   = couponRates;
    global.prices        = prices;
    global.recoveryCurve = recoveryCurve;
    global.payAccOnDefault = payAccOnDefault;
    global.couponInterval = couponInterval;
    global.paymentDCC    = paymentDCC;
    global.stubType      = stubType;
    global.status        = SUCCESS;
    global.allDates      = allDates;
    global.protectStart  = protectStart;
    global.isPriceClean  = isPriceClean;
    global.delay         = delay;
    global.badDayConv    = badDayConv;
    global.calendar      = calendar;

    if (GtoImslZerosSysEqn (matchBenchmarkPrices,
                            global.nbDate,
                            zcRates, /* initial guess */
                            1.0e-8, /* errBound */
                            100,    /* maxIter */
                            &cleanSpreads) != SUCCESS ||
        global.status != SUCCESS)
    {
        GtoErrMsg ("%s: Could not solve for smoothed curve\n", routine);
        goto done; /* failure */
    }

    sparseCurve = GtoMakeTCurve (global.today,
                                 global.endDates,
                                 cleanSpreads,
                                 global.nbDate,
                                 GTO_CONTINUOUS_BASIS,
                                 GTO_ACT_365F);
    if (sparseCurve == NULL)
        goto done; /* failure */

    smoothCurve           = CxCreditCurveMakeEmpty();
    smoothCurve->type     = cdsCurve->type;
    smoothCurve->tc       = CxZeroCurveMakeSmoother (sparseCurve,
                                                     global.allDates);
    smoothCurve->timestep = GtoDateIntervalMakeCopy(cdsCurve->timestep);
    if (CxCreditCurveValidate (smoothCurve) != SUCCESS)
        goto done; /* failure */

    status = SUCCESS;

 done:

    if (status != SUCCESS)
    {
        GtoErrMsgFailure (routine);
        CxCreditCurveFree (smoothCurve);
        smoothCurve = NULL;
    }

    GtoFreeTCurve (sparseCurve);
    GtoFreeSafe (cleanSpreads); /* allocated by ALIB */
    GtoFreeDateList (allDates);
    FREE (zcDates);
    FREE (zcRates);

    return smoothCurve;
}


/*
***************************************************************************
** The main bootstrap routine.
**
** In fact this routine immediately re-directs itself to either the flow
** curve method (a.k.a. CX) or the exotic curve method (a.k.a. CRX).
***************************************************************************
*/
CxTCreditCurve* CxCdsBootstrap
(TDate           today,           /* (I) Mostly ignored                     */
 TCurve         *discountCurve,   /* (I) Risk-free discount curve           */
 TDate           startDate,       /* (I) Start of CDS for accrual and risk  */
 TDate           valueDate,       /* (I) Settlement date for price          */
 long            nbDate,          /* (I) Number of benchmark dates          */
 TDate          *endDates,        /* (I) Maturity dates of CDS to bootstrap */
 double         *couponRates,     /* (I) CouponRates (e.g. 0.05 = 5% = 500bp)   */ 
 double         *prices,          /* (I) Prices (a.k.a. upfront charge)
                                     Can be NULL if no upfront charges.     */
 TBoolean       *includes,        /* (I) Include this date. Can be NULL if
                                     all are included.                      */
 CxTRecoveryCurve *recoveryCurve, /* (I) Recovery curve                     */
 TBoolean        payAccOnDefault, /* (I) Pay accrued on default             */
 TDateInterval  *couponInterval,  /* (I) Interval between fee payments      */
 CxTDayCountConv paymentDCC,     /* (I) DCC for fee payments and accrual    */
 CxTStubType     stubType,       /* (I) Stub type for fee leg               */
 CxTCreditCurveType curveType,    /* (I) Type of curve to create            */
 TDateInterval  *timestep,        /* (I) Integration timestep. NULL =
                                     couponInterval */
 TDateInterval  *smoothInterval,  /* (I) Smoothing interval. NULL if no
                                     smoothing.                             */
 TBoolean        protectStart,    /* (I) Protected on both start and end
                                     date of the CDS */
 TBoolean        isPriceClean,
 long            delay,
 CxTBadDayConv   badDayConv,
 CxTCalendar    *calendar
)
{
    static char routine[] = "CxCdsBootstrap";
    CxTCreditCurve *out = NULL;

    TDate           *includeEndDates = NULL;
    double          *includeCouponRates = NULL;
    double          *includePrices = NULL;

    TDateInterval ivl3M;

    SET_TDATE_INTERVAL(ivl3M,3,'M');
    if (couponInterval == NULL)
        couponInterval = &ivl3M;

    /* put common requirements between two routines here */
    REQUIRE (discountCurve != NULL);
    REQUIRE (nbDate > 0);
    REQUIRE (endDates != NULL);
    REQUIRE (couponRates != NULL);
    REQUIRE (recoveryCurve != NULL);
    REQUIRE (valueDate != 0 || prices == NULL);

    if (timestep == NULL)
        timestep = couponInterval;

    if (includes != NULL)
    {
        /* need to pick and choose which names appear */
        int  i;
        long nbInclude = 0;
        long j;
        for (i = 0; i < nbDate; ++i) 
        {
            if (includes[i]) ++nbInclude;
        }
        REQUIRE (nbInclude > 0);

        includeEndDates    = NEW_ARRAY(TDate, nbInclude);
        includeCouponRates = NEW_ARRAY(double, nbInclude);
        if (prices != NULL) 
            includePrices = NEW_ARRAY(double, nbInclude);
        
        j = 0;
        for (i = 0; i < nbDate; ++i)
        {
            if (includes[i])
            {
                includeEndDates[j]    = endDates[i];
                includeCouponRates[j] = couponRates[i];
                if (prices != NULL)
                    includePrices[j] = prices[i];
                ++j;
            }
        }
        ASSERT (j == nbInclude);
        
        nbDate      = nbInclude;
        endDates    = includeEndDates;
        couponRates = includeCouponRates;
        prices      = includePrices;
    }

    switch (curveType)
    {
    case CX_CURVE_TYPE_FLOW:
        out = CdsBootstrapFlowType (today,
                                    discountCurve,
                                    startDate,
                                    valueDate,
                                    nbDate,
                                    endDates,
                                    couponRates,
                                    prices,
                                    recoveryCurve,
                                    payAccOnDefault,
                                    couponInterval,
                                    paymentDCC,
                                    stubType,
                                    timestep,
                                    smoothInterval,
                                    protectStart,
                                    isPriceClean,
                                    delay,
                                    badDayConv,
                                    calendar);
        break;
    case CX_CURVE_TYPE_EXOTIC:
        out = CdsBootstrapExoticType (today,
                                      discountCurve,
                                      startDate,
                                      valueDate,
                                      nbDate,
                                      endDates,
                                      couponRates,
                                      prices,
                                      recoveryCurve,
                                      payAccOnDefault,
                                      couponInterval,
                                      paymentDCC,
                                      stubType,
                                      timestep,
                                      smoothInterval,
                                      protectStart,
                                      isPriceClean,
                                      delay,
                                      badDayConv,
                                      calendar);
        break;
    default:
        GtoErrMsg ("%s: Invalid curve type %d\n", routine, (int)curveType);
    }

 done:
    FREE(includeEndDates);
    FREE(includeCouponRates);
    FREE(includePrices);
    if (out == NULL)
        GtoErrMsgFailure (routine);

    return out;
}

/*
***************************************************************************
** This is the CDS bootstrap routine using CX analytics.
**
** Very little attempt has been made at extreme optimisation - this is quite
** a basic bootstrap routine which simply calls the underlying CDS pricer
** for each benchmark instrument while it changes the CDS zero rate at the
** maturity date of the benchmark instrument.
**
** After the first pass, it will optionally use a smoothing algorithm to
** create a curve with relatively continuous forward hazard rates.
***************************************************************************
*/
static CxTCreditCurve* CdsBootstrapFlowType
(TDate           today,           /* (I) Mostly ignored                     */
 TCurve         *discountCurve,   /* (I) Risk-free discount curve           */
 TDate           startDate,       /* (I) Start of CDS for accrual and risk  */
 TDate           valueDate,       /* (I) Settlement date for price          */
 long            nbDate,          /* (I) Number of benchmark dates          */
 TDate          *endDates,        /* (I) Maturity dates of CDS to bootstrap */
 double         *couponRates,     /* (I) CouponRates (e.g. 0.05 = 5% = 500bp)   */ 
 double         *prices,          /* (I) Prices (a.k.a. upfront charge)
                                     Can be NULL if no upfront charges.     */
 CxTRecoveryCurve *recoveryCurve, /* (I) Recovery curve                     */
 TBoolean        payAccOnDefault, /* (I) Pay accrued on default         */
 TDateInterval  *couponInterval,  /* (I) Interval between fee payments      */
 CxTDayCountConv  paymentDCC,     /* (I) DCC for fee payments and accrual   */
 CxTStubType      stubType,       /* (I) Stub type for fee leg              */
 TDateInterval  *timestep,        /* (I) Integration timestep .             */
 TDateInterval  *smoothInterval,  /* (I) Smoothing interval. NULL if no
                                     smoothing.                             */
 TBoolean        protectStart,    /* (I) Protected on both start and end
                                         date of the CDS */
 TBoolean        isPriceClean,    /* (I) Pay accrued at start date      */
 long            delay,
 CxTBadDayConv   badDayConv,
 CxTCalendar    *calendar)
{
    static char routine[] = "CdsBootstrapFlowType";
    int         status    = FAILURE;

    CxTCreditCurve *cdsCurve = NULL;
    int            i;
    
    CDS_BOOTSTRAP_CONTEXT context;
    CxTContingentLeg *cl = NULL;
    CxTFeeLeg        *fl = NULL;
    double           settleDiscount = 0.0;

    /* the spreads are of course complete rubbish at this point */
    /* however we will have validated that endDates are in increasing
       order */
    cdsCurve       = CxCreditCurveMakeEmpty ();
    cdsCurve->type = CX_CURVE_TYPE_FLOW;
    /* we work with a continuously compounded curve since that is faster -
       but we will convert to annual compounded since that is traditional */
    cdsCurve->tc   = GtoMakeTCurve (today,
                                    endDates,
                                    couponRates,
                                    nbDate,
                                    GTO_CONTINUOUS_BASIS,
                                    GTO_ACT_365F);
    if (cdsCurve->tc == NULL)
        goto done; /* failure */
    if (timestep != NULL)
        cdsCurve->timestep = GtoDateIntervalMakeCopy(timestep);
    if (CxCreditCurveValidate (cdsCurve) != SUCCESS)
        goto done; /* failure */

    context.discountCurve = discountCurve;
    context.cdsCurve      = cdsCurve;
    context.recoveryCurve = recoveryCurve;
    context.isPriceClean  = isPriceClean;
    
    if (valueDate != 0)
    {
        settleDiscount = CxForwardZeroPrice (discountCurve,
                                             cdsCurve->tc->fBaseDate,
                                             valueDate);
    }

    for (i = 0; i < nbDate; ++i)
    {
        double recovery;
        double guess;
        double spread;
        double price;

        if (CxRecoveryCurveInterp (recoveryCurve, endDates[i],
                                   &recovery) != SUCCESS)
            goto done; /* failure */

        /* note that the guess is not very good if the price is non-zero 
           or if the recovery is not constant */
        /* c'est la vie */
        guess = couponRates[i] / (1.0 - recovery);

        if (prices == NULL) price = 0.0;
        else price = prices[i];

        cl = CxCdsContingentLegMake (startDate,
                                     endDates[i],
                                     1.0,
                                     0,
                                     protectStart);
        if (cl == NULL)
            goto done; /* failure */

        fl = CxCdsFeeLegMake (startDate,
                              endDates[i],
                              payAccOnDefault,
                              couponInterval,
                              stubType,
                              1.0,
                              couponRates[i],
                              paymentDCC,
                              badDayConv,
                              calendar,
                              protectStart);
        if (fl == NULL)
            goto done; /* failure */

        context.i  = i;
        context.cl = cl;
        context.fl = fl;
        context.pvPrice = price * settleDiscount;

        if (GtoRootFindBrent ((TObjectFunc)cdsBootstrapPointFunction,
                              (void*) &context,
                              0.0,    /* boundLo */
                              1e10,   /* boundHi */
                              100,    /* numIterations */
                              guess,
                              0.0005, /* initialXstep */
                              0,      /* initialFDeriv */
                              1e-10,  /* xacc */
                              1e-10,  /* facc */
                              &spread) != SUCCESS)
        {
            char dateBuf[16];
            GtoErrMsg ("%s: Could not add CDS maturity %s spread %.2fbp\n",
                     routine,
                     CxDateFormat (endDates[i], "DD-MMM-YYYY", dateBuf),
                     1e4 * couponRates[i]);
            goto done; /* failure */
        }
        cdsCurve->tc->fArray[i].fRate = spread;

        CxContingentLegFree (cl);
        CxFeeLegFree (fl);
        cl = NULL;
        fl = NULL;
    }

    if (smoothInterval != NULL)
    {
        CxTCreditCurve* cdsSmoothCurve = CxCdsBootstrapSmoothing (
            cdsCurve,
            discountCurve,
            startDate,
            valueDate,
            nbDate,
            endDates,
            couponRates,
            prices,
            recoveryCurve,
            payAccOnDefault,
            couponInterval,
            paymentDCC,
            stubType,
            smoothInterval,
            protectStart,
            isPriceClean,
            delay,
            badDayConv,
            calendar);
        if (cdsSmoothCurve == NULL)
            goto done; /* failure */

        CxCreditCurveFree (cdsCurve);
        cdsCurve = cdsSmoothCurve;
    }

    if (CxCreditCurveConvertRateType (cdsCurve, GTO_ANNUAL_BASIS) != SUCCESS)
        goto done; /* failure */

    status = SUCCESS;

 done:

    if (status != SUCCESS)
    {
        CxCreditCurveFree (cdsCurve);
        cdsCurve = NULL;
        GtoErrMsgFailure (routine);
    }

    CxContingentLegFree (cl);
    CxFeeLegFree (fl);
        
    return cdsCurve;
}


static CxTCreditCurve* CdsBootstrapExoticType
(TDate           today,           /* (I) Mostly ignored                     */
 TCurve         *discountCurve,   /* (I) Risk-free discount curve           */
 TDate           startDate,       /* (I) Start of CDS for accrual and risk  */
 TDate           valueDate,       /* (I) Settlement date for price          */
 long            nbDate,          /* (I) Number of benchmark dates          */
 TDate          *endDates,        /* (I) Maturity dates of CDS to bootstrap */
 double         *couponRates,     /* (I) CouponRates (e.g. 0.05 = 5% = 500bp)   */ 
 double         *prices,          /* (I) Prices (a.k.a. upfront charge)
                                     Can be NULL if no upfront charges.     */
 CxTRecoveryCurve *recoveryCurve, /* (I) Recovery curve                     */
 TBoolean        payAccOnDefault,  /* (I) Pay accrued on default        */
 TDateInterval  *couponInterval,     /* (I) Interval between fee payments      */
 CxTDayCountConv  paymentDCC,      /* (I) DCC for fee payments and accrual   */
 CxTStubType      stubType,        /* (I) Stub type for fee leg              */
 TDateInterval  *timestep,        /* (I) Integration timestep .             */
 TDateInterval  *smoothInterval,  /* (I) Smoothing interval. NULL if no
                                     smoothing.                             */
 TBoolean        protectStart, /* (I) Protected on both start and end
                                         date of the CDS */
 TBoolean        isPriceClean,   /* (I) Pay accrued at start date      */
 long            delay,
 CxTBadDayConv   badDayConv,
 CxTCalendar    *calendar)
{
    static char routine[] = "CdsBootstrapExoticType";
    int         status    = FAILURE;

    TCurve        *ccStub = NULL;
    TCurve        *cc     = NULL;
    TCurve        *tcDisc = NULL;
    TDateInterval  ivlDelay;
    CxTCreditCurve *out = NULL;
    KFeeLeg_D     *fl = NULL;
    int            i;
    double         recovery;
    
    KAccrualConv     payAccr;

    SET_TDATE_INTERVAL(ivlDelay, delay, 'D');

    /* avoid unused parameters warning */
    valueDate = valueDate;
    smoothInterval = smoothInterval;
    protectStart = protectStart;

    REQUIRE (isPriceClean == FALSE);
    REQUIRE (recoveryCurve != NULL);

    if (recoveryCurve->numItems != 1)
    {
        GtoErrMsg ("%s: Recovery curve has more than one recovery rate. "
                   "This is not handled by exotic bootstrapper.\n", routine);
        goto done; /* failure */
    }
    recovery = recoveryCurve->recoveryRates[0];

    if (payAccOnDefault)
        payAccr = ACCRUAL_PAY_ALL;
    else
        payAccr = ACCRUAL_PAY_NONE;

    if (prices != NULL)
    {
        for (i = 0; i < nbDate; ++i)
            REQUIRE (IS_ALMOST_ZERO(prices[i]));
    }

    /*
    ** Convert the input discount curve into a TCurve with baseDate equal
    ** to today. Otherwise the CRX bootstrapping adjusts the credit curve
    ** with unfortunate results.
    */
    if (discountCurve->fBaseDate != today)
    {
        tcDisc = GtoZCShift2 (discountCurve, today, GTO_FLAT_FORWARDS);
        if (tcDisc == NULL)
            goto done; /* failure */
        discountCurve = tcDisc;
    }
    
    /*
    ** We add one at a time and hence we don't use the capacity of
    ** CreditZCBootstrap to add a number of on-cycle CDS all at the
    ** same time.
    */
    for (i = 0; i < nbDate; ++i)
    {
        int maturityIdx[1];

        fl = CrxFeeLegCreateFromFreq (startDate,
                                      endDates[i],
                                      *couponInterval,
                                      (KStubLocation)stubType,
                                      1.0, /* notional */
                                      couponRates[i],
                                      paymentDCC,
                                      payAccr,
                                      *timestep);
		if (fl == NULL)
			goto done; /* failure */
	

        maturityIdx[0] = fl->mNbCF - 1;

        cc = CreditZCBootstrap (today,
                                today,
                                paymentDCC,
                                payAccr,
                                PAY_DEF, /* pay protection at default */
                                ivlDelay,
                                fl->mNbCF,
                                fl->mAccStDates,
                                fl->mAccEndDates,
                                fl->mPayDates,
                                1,
                                maturityIdx,
                                &couponRates[i],
                                recovery,
                                *timestep,
                                discountCurve,
                                ccStub);
        GtoFreeTCurve (ccStub);
        ccStub = cc;
        if (cc == NULL)
            goto done; /* failure */

        CrxFeeLegFree (fl);
        fl = NULL;
    }

    out = CxCreditCurveMake (CX_CURVE_TYPE_EXOTIC, cc, timestep);
    status = SUCCESS;

 done:

    if (status != SUCCESS)
    {
        CxCreditCurveFree (out);
        out = NULL;
        GtoErrMsgFailure (routine);
    }

    GtoFreeTCurve (tcDisc);
    GtoFreeTCurve (cc);
    CrxFeeLegFree (fl);
  
    return out;
}



/*
 * objective function for root-solver
 *
 * returns the PV as at the startDate of the CDS curve
 */
static int cdsBootstrapPointFunction
(double   cleanSpread,
 void    *data,
 double  *pv)
{
    static char routine[] = "cdsBootstrapPointFunction";
    int         status    = FAILURE;

    CDS_BOOTSTRAP_CONTEXT *context = (CDS_BOOTSTRAP_CONTEXT*)data;

    int               i             = context->i;
    TCurve           *discountCurve = context->discountCurve;
    CxTCreditCurve   *cdsCurve      = context->cdsCurve;
    CxTRecoveryCurve *recoveryCurve = context->recoveryCurve;
    CxTContingentLeg *cl            = context->cl;
    CxTFeeLeg        *fl            = context->fl;
    double            pvPrice       = context->pvPrice;
    TBoolean          isPriceClean  = context->isPriceClean;
    TDate             cdsBaseDate   = cdsCurve->tc->fBaseDate;

    double           pvC; /* PV of contingent leg */
    double           pvF; /* PV of fee leg */

    cdsCurve->tc->fArray[i].fRate = cleanSpread;

    if (CxContingentLegPV (cl,
                           cdsBaseDate,
                           cdsBaseDate,
                           discountCurve,
                           cdsCurve,
                           recoveryCurve,
                           &pvC) != SUCCESS)
        goto done; /* failure */
                              
    if (CxFeeLegPV (fl,
                    cdsBaseDate,
                    cdsBaseDate,
                    discountCurve,
                    cdsCurve,
                    isPriceClean,
                    &pvF) != SUCCESS)
        goto done; /* failure */

    /* Note: price is discounted to cdsBaseDate */
    *pv = pvC - pvF - pvPrice;
    status = SUCCESS;

 done:

    if (status != SUCCESS)
        GtoErrMsgFailure (routine);

    return status;
}

/**
 * Objective function for multi-dimensional solver.
 */
static void matchBenchmarkPrices
(int     n,             /* (I) */
 double *cleanSpreads,  /* (I) */ 
 double *priceDiffs)   /* (O) */
{
    static char routine[] = "matchBenkmarkPrices";
    int         status    = FAILURE;

/**
 * 1. Create sparse zero curve with the given cleanSpreads.
 * 2. Convert this into a smoother and denser zero curve.
 * 3. Compute the couponRates for each benchmark instrument.
 * 4. Return the spread differences.
 */
    CxTCreditCurve *sparseCurve = NULL;
    CxTCreditCurve *smoothCurve = NULL;

    int  i;
    char dateBuf[16]; /* used for date formatting */
    TDate valueDate;

    REQUIRE (n == global.nbDate);

    if (global.valueDate == 0)
        valueDate = global.startDate;
    else
        valueDate = global.valueDate;

    sparseCurve = CxCreditCurveNew (CX_CURVE_TYPE_FLOW,
                                    global.today,
                                    global.nbDate,
                                    global.endDates,
                                    cleanSpreads,
                                    GTO_CONTINUOUS_BASIS,
                                    CX_ACT_365F,
                                    NULL);
    if (sparseCurve == NULL)
        goto done; /* failure */

    smoothCurve = CxCreditCurveMakeSmoother (sparseCurve,
                                             global.allDates);
    if (smoothCurve == NULL)
        goto done; /* failure */

    for (i = 0; i < global.nbDate; ++i)
    {
        double price;
        if (CxCdsPrice (global.today,
                        valueDate,
                        global.startDate,
                        global.endDates[i],
                        global.delay,
                        global.couponRates[i],
                        global.payAccOnDefault,
                        global.couponInterval,
                        global.stubType,
                        global.paymentDCC,
                        global.badDayConv,
                        global.calendar,
                        global.discountCurve,
                        smoothCurve,
                        global.recoveryCurve,
                        global.protectStart,
                        global.isPriceClean,
                        &price) != SUCCESS)
        {
            GtoErrMsg ("%s: Could not calculate price for %s in solver\n",
                     routine,
                     CxDateFormat (global.endDates[i], "YYYY-MM-DD",
                                   dateBuf));
            goto done; /* failure */
        }
        if (global.prices == NULL)
            priceDiffs[i] = price;
        else
            priceDiffs[i] = price - global.prices[i];
    }

    status = SUCCESS;

 done:

    if (status != SUCCESS)
    {
        /* setting priceDiffs to zero will end the iteration */
        /* setting global.status to FAILURE will tell the outside routine
           that there was a failure in here */
        global.status = FAILURE;
        for (i = 0; i < global.nbDate; ++i)
            priceDiffs[i] = 0.0;
        GtoErrMsgFailure (routine);
    }
    
    CxCreditCurveFree (smoothCurve);
    CxCreditCurveFree (sparseCurve);
}

