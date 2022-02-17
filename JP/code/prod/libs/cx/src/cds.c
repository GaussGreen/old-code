/*
***************************************************************************
** FILE NAME: cds.c
**
** CDS functions.
***************************************************************************
*/

#include "cds.h"

#include "contingentleg.h"
#include "feeleg.h"

#include <cxutils/include/calendar.h>
#include <cxutils/include/datelist.h>
#include <cxutils/include/cxmacros.h>
#include <cxutils/include/zerocurve.h>
#include <alib/dtivlo.h>

/*
 * Makes a contingent leg for a vanilla CDS
 *
 * Protection starts either at the beginning of startDate (protectStart=True)
 * or at the end of startDate.
 *
 * Protection ends at the end of endDate.
 *
 * Notional is the amount of notional protected.
 * Delay is the delay in payment after default measured in days.
 */
CxTContingentLeg* CxCdsContingentLegMake
(TDate     startDate,
 TDate     endDate,  
 double    notional, 
 long      delay,    
 TBoolean  protectStart)
{
    static char routine[] = "CxCdsContingentLegMake";
    int         status    = FAILURE;
    
    CxTContingentLeg *cl = NULL;
    
    REQUIRE (endDate > startDate);
    REQUIRE (delay >= 0);

    cl = CxContingentLegMakeEmpty(1);
    if (cl == NULL)
        goto done; /* failure */

    /* cl->startDate is defined as giving protection from end of startDate.
       So if we want to protect on the start date, we need to move this
       date forward by one. */
    if (protectStart)
    {
        cl->startDate = startDate-1;
    }
    else
    {
        cl->startDate = startDate;
    }
    cl->dates[0]     = endDate;
    cl->notionals[0] = notional;
    cl->payType      = CX_PROT_PAY_DEF;
    cl->payDelay     = delay;

    status = SUCCESS;

 done:

    if (status != SUCCESS)
    {
        CxContingentLegFree(cl);
        cl = NULL;
        GtoErrMsgFailure(routine);
    }

    return cl;
}

/*
 * computes the PV for a contingent leg for a vanilla CDS
 */
int CxCdsContingentLegPV
(TDate             today,
 TDate             valueDate,
 TDate             startDate,
 TDate             endDate,
 double            notional,
 long              delay,
 TCurve           *discCurve,
 CxTCreditCurve   *spreadCurve,
 CxTRecoveryCurve *recoveryCurve,
 TBoolean          protectStart,
 double           *pv)
{
    static char routine[] = "CxCdsContingentLegPV";
    int         status    = FAILURE;

    CxTContingentLeg *cl = NULL;

    cl = CxCdsContingentLegMake (startDate, endDate, notional, delay,
                                 protectStart);
    if (cl == NULL)
        goto done; /* failure */

    if (CxContingentLegPV (cl, today, valueDate, discCurve, spreadCurve,
                           recoveryCurve, pv) != SUCCESS)
        goto done; /* failure */

    status = SUCCESS;

 done:

    if (status != SUCCESS)
        GtoErrMsgFailure (routine);

    CxContingentLegFree (cl);

    return status;
}

/*
 * makes a fixed fee leg for a vanilla CDS
 *
 * note that you are protected on both startDate and endDate
 */
CxTFeeLeg* CxCdsFeeLegMake
(TDate           startDate,
 TDate           endDate,
 TBoolean        payAccOnDefault,
 TDateInterval  *dateInterval,
 CxTStubType     stubType,
 double          notional,
 double          couponRate,
 CxTDayCountConv paymentDcc,
 CxTBadDayConv   badDayConv,
 CxTCalendar    *calendar,
 TBoolean        protectStart)
{
    static char routine[] = "CxCdsFeeLegMake";
    int         status    = FAILURE;

    TDateList  *dl = NULL;
    CxTFeeLeg  *fl = NULL;
    TDateInterval ivl3M;
    TDate         prevDate;
    TDate         prevDateAdj;

    int i;

    SET_TDATE_INTERVAL(ivl3M,3,'M');
    if (dateInterval == NULL)
        dateInterval = &ivl3M;

    REQUIRE (endDate > startDate);

    dl = CxDateListMakeRegular (startDate, endDate, dateInterval, stubType);
    if (dl == NULL)
        goto done; /* failure */

    /* the datelist includes both start date and end date */
    /* therefore it has one more element than the fee leg requires */

    fl = CxFeeLegMakeEmpty (dl->fNumItems-1);
    if (fl == NULL)
        goto done; /* failure */

    if (payAccOnDefault)
    {
        fl->accrualPayConv = CX_ACCRUAL_PAY_ALL;
    }
    else
    {
        fl->accrualPayConv = CX_ACCRUAL_PAY_NONE;
        /* and we will assume that it observes at end of the period */
    }
    fl->dcc = paymentDcc;

    prevDate = dl->fArray[0];
    prevDateAdj = prevDate; /* first date is not bad day adjusted */

    for (i = 0; i < fl->nbDates; ++i)
    {
        TDate nextDate = dl->fArray[i+1];
        TDate nextDateAdj;

        if (CxBusinessDay (nextDate, badDayConv, calendar, 
                           &nextDateAdj) != SUCCESS)
            goto done; /* failure */

        fl->accStartDates[i] = prevDateAdj;
        fl->accEndDates[i]   = nextDateAdj;
        fl->payDates[i]      = nextDateAdj;
        fl->notionals[i]     = notional;
        fl->couponRates[i]   = couponRate;

        prevDate    = nextDate;
        prevDateAdj = nextDateAdj;
    }

    /* the last accrual date is not adjusted */
    /* also we may have one extra day of accrued interest */
    if (protectStart)
    {
        fl->accEndDates[fl->nbDates-1] = prevDate+1;
        fl->obsStartOfDay              = TRUE;
    }
    else
    {
        fl->accEndDates[fl->nbDates-1] = prevDate;
        fl->obsStartOfDay              = FALSE;
    }
    status = SUCCESS;

 done:

    if (status != SUCCESS)
    {
        GtoErrMsgFailure(routine);
        CxFeeLegFree (fl);
        fl = NULL;
    }

    GtoFreeDateList (dl);

    return fl;
}

/*
 * computes the PV for a fee leg for a vanilla CDS
 *
 * note that you are protected on both start date and end date
 */
int CxCdsFeeLegPV
(TDate           today,
 TDate           valueDate,
 TDate           startDate,
 TDate           endDate,
 TBoolean        payAccOnDefault,
 TDateInterval  *dateInterval,
 CxTStubType      stubType,
 double          notional,
 double          couponRate,
 CxTDayCountConv  paymentDcc,
 CxTBadDayConv    badDayConv,
 CxTCalendar     *calendar,
 TCurve         *discCurve,
 CxTCreditCurve  *spreadCurve,
 TBoolean        protectStart,
 TBoolean        cleanPrice,
 double         *pv)
{
    static char routine[] = "CxCdsFeeLegPV";
    int         status    = FAILURE;

    CxTFeeLeg *fl = NULL;

    fl = CxCdsFeeLegMake (startDate, endDate, payAccOnDefault,
                          dateInterval, stubType, notional,
                          couponRate, paymentDcc, badDayConv,
                          calendar, protectStart);
    if (fl == NULL)
        goto done; /* failure */

    if (CxFeeLegPV (fl, today, valueDate, discCurve, spreadCurve, cleanPrice,
                    pv) != SUCCESS)
        goto done; /* failure */

    status = SUCCESS;

 done:

    if (status != SUCCESS)
        GtoErrMsgFailure (routine);

    CxFeeLegFree (fl);
    return status;
}

/*
 * computes the par spread for a vanilla CDS
 *
 * note that you are protected on both start date and end date
 */
int CxCdsParSpread
(TDate           today,
 TDate           settleDate,
 TDate           startDate,
 TDate           endDate,
 long            delay,
 double          price,
 TBoolean        payAccOnDefault,
 TDateInterval  *dateInterval,
 CxTStubType      stubType,
 CxTDayCountConv  paymentDcc,
 CxTBadDayConv    badDayConv,
 CxTCalendar     *calendar,
 TCurve         *discCurve,
 CxTCreditCurve  *spreadCurve,
 CxTRecoveryCurve *recoveryCurve,
 TBoolean        protectStart,
 TBoolean        isPriceClean,
 double         *parSpread)
{
    static char routine[] = "CxCdsParSpread";
    int         status    = FAILURE;

    double      annuityPV;
    double      contingentLegPV;
    double      pricePV;
    TDate      valueDate;

    REQUIRE(delay == 0); /* because delay not yet applied to fee legs */
    REQUIRE(parSpread != NULL);
    REQUIRE(IS_ZERO(price) || settleDate > 0);
    /* all other requirements can be handled by the routines we call */

    valueDate = startDate;

    if (CxCdsFeeLegPV (today,
                       valueDate,
                       startDate,
                       endDate,
                       payAccOnDefault,
                       dateInterval,
                       stubType,
                       1.0, /* notional */
                       1.0, /* couponRate */
                       paymentDcc,
                       badDayConv,
                       calendar,
                       discCurve,
                       spreadCurve,
                       protectStart,
                       isPriceClean,
                       &annuityPV) != SUCCESS)
        goto done; /* failure */

    if (CxCdsContingentLegPV (today,
                              valueDate,
                              startDate,
                              endDate,
                              1.0, /* notional */
                              delay,
                              discCurve,
                              spreadCurve,
                              recoveryCurve,
                              protectStart,
                              &contingentLegPV) != SUCCESS)
        goto done; /* failure */

    if (IS_NOT_ZERO(price))
    {
        /* the upfront charge payment is not conditional on default */
        /* therefore we use the riskless zero curve to discount it  */
        /* typically the discounting is for one or two days         */
        double discount;

        discount = CxForwardZeroPrice(discCurve, valueDate, settleDate);
        pricePV = price * discount;
    }
    else
    {
        pricePV = 0.0;
    }

    *parSpread = (contingentLegPV - pricePV) / annuityPV;
    status = SUCCESS;

 done:

    if (status != SUCCESS)
        GtoErrMsgFailure (routine);

    return status;
}


/*
 * computes the price for a vanilla CDS
 *
 * note that you are protected on both start date and end date
 */
int CxCdsPrice
(TDate           today,
 TDate           settleDate,
 TDate           startDate,
 TDate           endDate,
 long            delay,
 double          couponRate,
 TBoolean        payAccOnDefault,
 TDateInterval  *dateInterval,
 CxTStubType      stubType,
 CxTDayCountConv  paymentDcc,
 CxTBadDayConv    badDayConv,
 CxTCalendar     *calendar,
 TCurve         *discCurve,
 CxTCreditCurve  *spreadCurve,
 CxTRecoveryCurve *recoveryCurve,
 TBoolean        protectStart,
 TBoolean        isPriceClean,
 double         *price)
{
    static char routine[] = "CxCdsPrice";
    int         status    = FAILURE;

    double      feeLegPV;
    double      contingentLegPV;
    TDate      valueDate;

    REQUIRE(delay == 0); /* because delay not yet applied to fee legs */
    REQUIRE(price != NULL);
    /* all other requirements can be handled by the routines we call */

    valueDate = settleDate;

    if (CxCdsFeeLegPV (today,
                       valueDate,
                       startDate,
                       endDate,
                       payAccOnDefault,
                       dateInterval,
                       stubType,
                       1.0, /* notional */
                       couponRate,
                       paymentDcc,
                       badDayConv,
                       calendar,
                       discCurve,
                       spreadCurve,
                       protectStart,
                       isPriceClean,
                       &feeLegPV) != SUCCESS)
        goto done; /* failure */

    if (CxCdsContingentLegPV (today,
                              valueDate,
                              startDate,
                              endDate,
                              1.0, /* notional */
                              delay,
                              discCurve,
                              spreadCurve,
                              recoveryCurve,
                              protectStart,
                              &contingentLegPV) != SUCCESS)
        goto done; /* failure */

    *price = contingentLegPV - feeLegPV;
    status = SUCCESS;

 done:

    if (status != SUCCESS)
        GtoErrMsgFailure (routine);

    return status;
}



/**
 * Computes the non-contingent cash flows for a fee leg. These are the
 * cash flows you will receive if there is no default.
 */
TCashFlowList* CxCdsFeeLegFlows
(TDate           startDate,
 TDate           endDate,
 TBoolean        payAccOnDefault,
 TDateInterval  *dateInterval,
 CxTStubType     stubType,
 double          notional,
 double          couponRate,
 CxTDayCountConv paymentDcc,
 CxTBadDayConv   badDayConv,
 CxTCalendar    *calendar,
 TBoolean        protectStart)
{
    static char routine[] = "CxCdsFeeLegFlows";

    TCashFlowList *cfl = NULL;
    CxTFeeLeg     *fl  = NULL;

    fl = CxCdsFeeLegMake (startDate,
                          endDate,
                          payAccOnDefault,
                          dateInterval,
                          stubType,
                          notional,
                          couponRate,
                          paymentDcc,
                          badDayConv,
                          calendar,
                          protectStart);

    if (fl == NULL)
        goto done; /* failure */

    cfl = CxFeeLegFlows (fl);

 done:

    CxFeeLegFree (fl);
    if (cfl == NULL)
        GtoErrMsgFailure (routine);

    return cfl;
}



/**
 * Computes the expected cash flows for a fee leg. This returns the cash
 * flows on the fee leg payment dates taking into account the survival
 * probability implied by the credit curve.
 *
 * The calculation involves computing the PV of each fee period, and then
 * forward valuing that value to the cash flow payment date. Although this
 * is not very sensitive to the discount curve, we therefore need a
 * discount curve.
 *
 * The guarantee is that if you subsequently PV these flows using the
 * discount curve, then you should get the same result as CxFeeLegPV.
 */
TCashFlowList* CxCdsFeeLegExpectedFlows
(TDate           today,
 TDate           startDate,
 TDate           endDate,
 TBoolean        payAccOnDefault,
 TDateInterval  *dateInterval,
 CxTStubType     stubType,
 double          notional,
 double          couponRate,
 CxTDayCountConv paymentDcc,
 CxTBadDayConv   badDayConv,
 CxTCalendar    *calendar,
 TCurve         *discCurve,
 CxTCreditCurve *creditCurve,
 TBoolean        protectStart)
{
    static char routine[] = "CxCdsFeeLegExpectedFlows";

    TCashFlowList *cfl = NULL;
    CxTFeeLeg     *fl  = NULL;

    fl = CxCdsFeeLegMake (startDate,
                          endDate,
                          payAccOnDefault,
                          dateInterval,
                          stubType,
                          notional,
                          couponRate,
                          paymentDcc,
                          badDayConv,
                          calendar,
                          protectStart);

    if (fl == NULL)
        goto done; /* failure */

    cfl = CxFeeLegExpectedFlows (fl,
                                 today,
                                 discCurve,
                                 creditCurve);

 done:

    CxFeeLegFree (fl);
    if (cfl == NULL)
        GtoErrMsgFailure (routine);

    return cfl;
}


/**
 * Computes the expected cash flows for a contingent leg. This returns the
 * cash flows on the corresponding payment dates taking into account
 * the survival probability implied by the credit curve.
 *
 * The calculation involves computing the PV of each fee period, and then
 * forward valuing that value to the fee cash flow payment date. Although
 * this is not very sensitive to the discount curve, we therefore need
 * a discount curve. We also need a fee leg to provide the dates.
 *
 * The guarantee is that if you subsequently PV these flows using the
 * discount curve, then you should get the same results as
 * CxContingentLegPV.
 */
TCashFlowList* CxCdsContingentLegExpectedFlows
(TDate             today,
 TDate             startDate,
 TDate             endDate,  
 TDateInterval    *dateInterval,
 CxTStubType       stubType,
 double            notional, 
 long              delay,    
 CxTBadDayConv     badDayConv,
 CxTCalendar      *calendar,
 TCurve           *discCurve,
 CxTCreditCurve   *creditCurve,
 CxTRecoveryCurve *recoveryCurve,
 TBoolean          protectStart)
{
    static char routine[] = "CxCdsContingentLegExpectedFlows";

    TCashFlowList    *cfl = NULL;
    CxTFeeLeg        *fl  = NULL;
    CxTContingentLeg *cl  = NULL;
    TDateList        *dl  = NULL;

    fl = CxCdsFeeLegMake (startDate,
                          endDate,
                          FALSE,
                          dateInterval,
                          stubType,
                          notional,
                          1.0,
                          CX_ACT_360,
                          badDayConv,
                          calendar,
                          protectStart);

    if (fl == NULL)
        goto done; /* failure */

    dl = GtoNewDateListFromDates (fl->payDates,
                                  fl->nbDates);
    if (dl == NULL)
        goto done; /* failure */

    cl = CxCdsContingentLegMake (startDate,
                                 endDate,
                                 notional,
                                 delay,
                                 protectStart);
    if (cl == NULL)
        goto done; /* failure */

    cfl = CxContingentLegExpectedFlows (cl,
                                        dl,
                                        today,
                                        discCurve,
                                        creditCurve,
                                        recoveryCurve);

 done:

    CxContingentLegFree (cl);
    CxFeeLegFree (fl);
    GtoFreeDateList (dl);

    if (cfl == NULL)
        GtoErrMsgFailure (routine);

    return cfl;
}


TCashFlowList* CxCdsExpectedFlows
(TDate             today,
 TDate             startDate,
 TDate             endDate,
 long              delay,
 double            notional,
 double            couponRate,
 TBoolean          payAccOnDefault,
 TDateInterval    *dateInterval,
 CxTStubType       stubType,
 CxTDayCountConv   paymentDcc,
 CxTBadDayConv     badDayConv,
 CxTCalendar      *calendar,
 TCurve           *discCurve,
 CxTCreditCurve   *creditCurve,
 CxTRecoveryCurve *recoveryCurve,
 TBoolean          protectStart)
{
    static char routine[] = "CxCdsExpectedFlows";

    TCashFlowList    *cfl = NULL;

    TCashFlowList    *cflF = NULL;
    TCashFlowList    *cflC = NULL;
    CxTFeeLeg        *fl   = NULL;
    CxTContingentLeg *cl   = NULL;
    TDateList        *dl   = NULL;

    fl = CxCdsFeeLegMake (startDate,
                          endDate,
                          payAccOnDefault,
                          dateInterval,
                          stubType,
                          notional,
                          couponRate,
                          paymentDcc,
                          badDayConv,
                          calendar,
                          protectStart);

    if (fl == NULL)
        goto done; /* failure */

    dl = GtoNewDateListFromDates (fl->payDates,
                                  fl->nbDates);
    if (dl == NULL)
        goto done; /* failure */

    cl = CxCdsContingentLegMake (startDate,
                                 endDate,
                                 -notional,
                                 delay,
                                 protectStart);
    if (cl == NULL)
        goto done; /* failure */

    cflC = CxContingentLegExpectedFlows (cl,
                                         dl,
                                         today,
                                         discCurve,
                                         creditCurve,
                                         recoveryCurve);
    if (cflC == NULL)
        goto done; /* failure */

    cflF = CxFeeLegExpectedFlows (fl,
                                  today,
                                  discCurve,
                                  creditCurve);
    if (cflF == NULL)
        goto done; /* failure */

    cfl = GtoMergeCFL (cflF, cflC);

 done:

    CxContingentLegFree (cl);
    CxFeeLegFree (fl);
    GtoFreeDateList (dl);
    GtoFreeCFL (cflF);
    GtoFreeCFL (cflC);

    if (cfl == NULL)
        GtoErrMsgFailure (routine);

    return cfl;
}
