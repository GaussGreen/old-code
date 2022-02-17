/*
***************************************************************************
** FILE NAME: feeLeg.c
**
** Analytics for a fee leg.
***************************************************************************
*/

#include "feeleg.h"

#include <math.h>

#include <cxutils/include/alibconv.h>
#include "creditcurve.h"
#include "crxconv.h"
#include "timeline.h"

#include <cxutils/include/cxmacros.h>
#include <cxutils/include/datelist.h>
#include <cxutils/include/dateutils.h>
#include <cxutils/include/zerocurve.h>

static int FeeLegPV
(CxTFeeLeg      *fl,
 TDate           today,
 TDate           valueDate,
 TCurve         *discCurve,
 CxTCreditCurve *spreadCurve,
 TBoolean        payAccruedAtStart,
 double         *pv);

static int FeePaymentPVWithTimeLine
(CxTAccrualPayConv accrualPayConv,
 TDate             today,
 TDate             accStartDate,
 TDate             accEndDate,
 TDate             payDate,
 CxTDayCountConv   accrueDCC,
 double            notional,
 double            couponRate,
 TCurve           *discCurve,
 CxTCreditCurve   *spreadCurve,
 TDateList        *tl,
 TBoolean          obsStartOfDay,
 double           *pv);



/*f
** Calculates the PV of a fee leg with fixed fee payments.
*/
int CxFeeLegPV
(CxTFeeLeg      *fl,
 TDate           today,
 TDate           valueDate,
 TCurve         *discCurve,
 CxTCreditCurve *spreadCurve,
 TBoolean        payAccruedAtStart,
 double         *pv)
{
    static char routine[] = "CxFeeLegPV";
    int         status    = FAILURE;

    REQUIRE (spreadCurve != NULL);

    switch (spreadCurve->type)
    {
    case CX_CURVE_TYPE_FLOW:
    {
        status = FeeLegPV (fl,
                           today,
                           valueDate,
                           discCurve,
                           spreadCurve,
                           payAccruedAtStart,
                           pv);
        break;
    }
    case CX_CURVE_TYPE_EXOTIC:
    {
        KFeeLeg_D *crxfl     = CxFeeLegConvert (fl, spreadCurve->timestep);
        KStubType  priceConv = payAccruedAtStart ? BOND : NONE;
        if (crxfl != NULL)
        {
            status = RiskyFeePV_O (pv,
                                   today,
                                   valueDate,
                                   crxfl,
                                   priceConv,
                                   discCurve,
                                   spreadCurve->tc);
        }
        CrxFeeLegFree (crxfl);
        break;
    }
    default:
        GtoErrMsg ("%s: Invalid curve type %d\n", routine,
                   (int)spreadCurve->type);
    }

 done:

    if (status != SUCCESS)
        GtoErrMsgFailure (routine);

    return status;
}

/*f
** Calculates the PV of a fee leg with fixed fee payments.
*/
static int FeeLegPV
(CxTFeeLeg      *fl,
 TDate           today,
 TDate           valueDate,
 TCurve         *discCurve,
 CxTCreditCurve *spreadCurve,
 TBoolean        payAccruedAtStart,
 double         *pv)
{
    static char routine[] = "CxFeeLegPV";
    int         status    = FAILURE;

    int    i;
    double myPv;
    double valueDatePv;
    TDateList  *tl = NULL;

    REQUIRE (fl != NULL);
    REQUIRE (discCurve != NULL);
    REQUIRE (spreadCurve != NULL);
    REQUIRE (pv != NULL);
    REQUIRE (valueDate >= today);

    REQUIRE (payAccruedAtStart == FALSE); /* TBD */

    myPv = 0.0;

    if (fl->nbDates > 1)
    {
        /* it is more efficient to compute the timeLine just the once
           and truncate it for each payment */
        TDate startDate = fl->accStartDates[0];
        TDate endDate   = fl->accEndDates[fl->nbDates-1];
        
        tl = CxRiskyTimeLine (startDate,
                              endDate,
                              discCurve,
                              spreadCurve,
                              NULL, /* lossCurve */
                              0); /* delay */
        if (tl == NULL)
            goto done; /* failure */
    }


    for (i = 0; i < fl->nbDates; ++i)
    {
        double thisPv;

        if (FeePaymentPVWithTimeLine (fl->accrualPayConv,
                                      today,
                                      fl->accStartDates[i],
                                      fl->accEndDates[i],
                                      fl->payDates[i],
                                      fl->dcc,
                                      fl->notionals[i],
                                      fl->couponRates[i],
                                      discCurve,
                                      spreadCurve,
                                      tl,
                                      fl->obsStartOfDay,
                                      &thisPv) != SUCCESS)
            goto done; /* failure */

        myPv += thisPv;
    }

    valueDatePv = CxForwardZeroPrice (discCurve, today, valueDate);
    
    status = SUCCESS;
    *pv = myPv / valueDatePv;

 done:

    if (status != SUCCESS)
        GtoErrMsgFailure (routine);

    GtoFreeDateList (tl);

    return status;
}

/**
 * Calculates the PV of a single fee payment.
 *
 * Calling this function repeatedly with sensible inputs and summing the
 * result would give the same answer as using CxFeeLegPV directly
 * with one difference - the pv returned by this function is for today,
 * whereas the pv returned by CxFeeLegPV is for valueDate.
 *
 * The conversion is a matter of dividing by the discount factor between
 * today and valueDate using the risk-free curve.
 */
int CxFeePaymentPV
(CxTAccrualPayConv accrualPayConv,
 TDate             today,
 TDate             accStartDate,
 TDate             accEndDate,
 TDate             payDate,
 CxTDayCountConv   accrueDCC,
 double            notional,
 double            couponRate,
 TCurve           *discCurve,
 CxTCreditCurve   *spreadCurve,
 TBoolean          obsStartOfDay,
 double           *pv)
{
    static char routine[] = "CxFeePaymentPV";
    int         status    = FAILURE;

    REQUIRE (spreadCurve != NULL);
    REQUIRE (spreadCurve->type == CX_CURVE_TYPE_FLOW); /* TBD */

    if (FeePaymentPVWithTimeLine (accrualPayConv,
                                  today,
                                  accStartDate,
                                  accEndDate,
                                  payDate,
                                  accrueDCC,
                                  notional,
                                  couponRate,
                                  discCurve,
                                  spreadCurve,
                                  NULL, /* tl */
                                  obsStartOfDay,
                                  pv) != SUCCESS)
        goto done; /* failure */
    
    status = SUCCESS;
    
 done:

    if (status != SUCCESS)
        GtoErrMsgFailure (routine);

    return status;
}

/**
 * Calculates the PV of a single fee payment.
 * Uses a pre-calculated timeline for efficiency.
 *
 * Calling this function repeatedly with sensible inputs and summing the
 * result would give the same answer as using CxFeeLegPV directly
 * with one difference - the pv returned by this function is for today,
 * whereas the pv returned by CxFeeLegPV is for valueDate.
 *
 * The conversion is a matter of dividing by the discount factor between
 * today and valueDate using the risk-free curve.
 */
static int FeePaymentPVWithTimeLine
(CxTAccrualPayConv accrualPayConv,
 TDate             today,
 TDate             accStartDate,
 TDate             accEndDate,
 TDate             payDate,
 CxTDayCountConv   accrueDCC,
 double            notional,
 double            couponRate,
 TCurve           *discCurve,
 CxTCreditCurve   *spreadCurve,
 TDateList        *tl,
 TBoolean          obsStartOfDay,
 double           *pv)
{
    static char routine[] = "CxFeePaymentPV";
    int         status    = FAILURE;

    double myPv = 0.0;

    /*
     * Because survival is calculated at the end of the day, then if
     * we observe survival at the start of the day, we need to subtract
     * one from the date.
     */
    int    obsOffset = obsStartOfDay ? -1 : 0;

    REQUIRE (discCurve != NULL);
    REQUIRE (spreadCurve != NULL);
    REQUIRE (pv != NULL);

    switch (accrualPayConv)
    {
    case CX_ACCRUAL_PAY_NONE:
    {
        /* fee leg pays at pay date if it has survived to accrual end date */
        double accTime;
        double amount;
        double survival;
        double discount;
        
        accTime  = CxDayCountFraction (accStartDate,
                                       accEndDate,
                                       accrueDCC);
        amount   = notional * couponRate * accTime;
        survival = CxForwardZeroPrice (spreadCurve->tc,
                                       today,
                                       accEndDate + obsOffset);
        discount = CxForwardZeroPrice (discCurve,
                                       today,
                                       payDate);
        myPv = amount * survival * discount;
        break;
    }
    case CX_ACCRUAL_PAY_ALL:
    {
        /* fee leg pays accrual on default - otherwise it pays at pay date
           if it has survived to accrual end date */
        double accTime;
        double amount;
        double survival;
        double discount;
        double accrual;
        
        accTime  = CxDayCountFraction (accStartDate,
                                       accEndDate,
                                       accrueDCC);
        amount   = notional * couponRate * accTime;
        survival = CxForwardZeroPrice (spreadCurve->tc,
                                       today,
                                       accEndDate + obsOffset);
        discount = CxForwardZeroPrice (discCurve,
                                       today,
                                       payDate);
        myPv = amount * survival * discount;
        
        /* also need to calculate accrual PV
           currently there is no delay inside the fixed fee leg,
           so it has to be paid immediately */
        
        if (AccrualOnDefaultPVWithTimeLine (today,
                                            accStartDate + obsOffset,
                                            accEndDate + obsOffset,
                                            0, /* delay */
                                            amount,
                                            discCurve,
                                            spreadCurve,
                                            tl,
                                            &accrual) != SUCCESS)
            goto done; /* failure */
        
        myPv += accrual;
        break;
    }
    default:
        GtoErrMsg ("%s: Invalid accrual payment type %d\n",  
                 routine, (int)accrualPayConv);
        goto done; /* failure */
    }

    status = SUCCESS;
    *pv = myPv;

 done:

    if (status != SUCCESS)
        GtoErrMsgFailure (routine);

    return status;
}

/*f
** Calculates the PV of the accruals which occur on default with delay.
*/
int CxAccrualOnDefaultPV
(TDate           today,
 TDate           startDate,
 TDate           endDate,
 long            delay,
 double          amount,
 TCurve         *discCurve,
 CxTCreditCurve *spreadCurve,
 TBoolean        obsStartOfDay,
 double         *pv)
{
    static char routine[] = "CxAccrualOnDefaultPV";
    int         status    = FAILURE;

    /*
     * Because survival is calculated at the end of the day, then if
     * we observe survival at the start of the day, we need to subtract
     * one from the date.
     */
    int    obsOffset = obsStartOfDay ? -1 : 0;

    REQUIRE (spreadCurve != NULL);
    REQUIRE (spreadCurve->type == CX_CURVE_TYPE_FLOW); /* TBD */

    if (AccrualOnDefaultPVWithTimeLine (today,
                                        startDate + obsOffset,
                                        endDate + obsOffset,
                                        delay,
                                        amount,
                                        discCurve,
                                        spreadCurve,
                                        NULL, /* tl */
                                        pv) != SUCCESS)
        goto done; /* failure */

    status = SUCCESS;

 done:

    if (status != SUCCESS)
        GtoErrMsgFailure (routine);

    return status;
}

/*f
** Calculates the PV of the accruals which occur on default with delay.
** Uses a pre-calculated timeline for efficiency.
*/
int AccrualOnDefaultPVWithTimeLine
(TDate           today,
 TDate           startDate,
 TDate           endDate,
 long            delay,
 double          amount,
 TCurve         *discCurve,
 CxTCreditCurve *spreadCurve,
 TDateList      *criticalDates,
 double         *pv)
{
    static char routine[] = "CxAccrualOnDefaultPVWithTimeLine";
    int         status    = FAILURE;

    double  myPv = 0.0;
    int     i;

    double t;
    double s0;
    double s1;
    double df0;
    double df1;
    double accRate;

    TDateList  *tl = NULL;

    REQUIRE (endDate > startDate);
    REQUIRE (discCurve != NULL);
    REQUIRE (spreadCurve != NULL);
    REQUIRE (pv != NULL);
    REQUIRE (delay >= 0);
    
    /*
    ** Timeline is points on the spreadCurve between startDate and endDate,
    ** combined with delay adjusted points from the discCurve, plus
    ** the startDate and endDate.
    */
    if (criticalDates != NULL)
    {
        tl = CxTruncateTimeLine (criticalDates, startDate, endDate);
    }
    else
    {
        tl = CxRiskyTimeLine (startDate, endDate, discCurve, spreadCurve,
                                 NULL, delay);
    }
    if (tl == NULL)
        goto done; /* failure */

    /* the integration - we can assume flat forwards between points on
       the timeline - this is true for both curves 

       we are integrating -Zt dS/dt where Z is the discount factor and
       S is the survival probability and t is the accrual time

       assuming flat forwards on each part of the integration, this is an
       exact integral

       we then make some adjustments to account for the fact that in reality
       we are sampling only once per business day, which on average over the
       year means there is about one day of delay in the integral
    */
    
    t       = (double)(endDate-startDate)/365.0;
    accRate = amount/t;
    s1      = CxForwardZeroPrice (spreadCurve->tc, today, startDate);
    df1     = CxForwardZeroPrice (discCurve, today, startDate+delay);
    for (i = 1; i < tl->fNumItems; ++i)
    {
        double lambda;
        double fwdRate;
        double thisPv;
        double t0;
        double t1;

        /* adj, adj2 and adjDays are all to do with sampling frequency error */
        double adj;
        double adj2;
        double adjDays=0.0; /* 0.5 for all days + weekends + holidays */
        TBoolean adjSamplingError = FALSE;
        
        s0  = s1;
        df0 = df1;
        s1  = CxForwardZeroPrice (spreadCurve->tc, today, tl->fArray[i]);
        df1 = CxForwardZeroPrice (discCurve, today, tl->fArray[i]+delay);
        t0  = (double)(tl->fArray[i-1] - startDate)/365.0;
        t1  = (double)(tl->fArray[i] - startDate)/365.0;
        t   = t1-t0;

        lambda  = log(s0/s1)/t;
        fwdRate = log(df0/df1)/t;

        thisPv  = lambda * accRate * s0 * df0 * (
            (t0 + 1.0/(lambda+fwdRate))/(lambda+fwdRate) -
            (t1 + 1.0/(lambda+fwdRate))/(lambda+fwdRate) * 
            s1/s0 * df1/df0);

        if (adjSamplingError)
        {
            /* it is not clear whether we ever want to do this - and this
               little bit of code is a noticeable performance drain
               given the number of times this routine is called */
            adj  = exp(-fwdRate*adjDays/365.0); 
            adj2 = lambda*adjDays/365.0 * 
                (1.0 - exp(-(lambda+fwdRate)*t)) / (lambda+fwdRate);
            myPv += (thisPv+adj2) * adj;
        }
        else
        {
            myPv += thisPv;
        }
    }

    status = SUCCESS;
    *pv = myPv;
        
 done:

    if (status != SUCCESS)
        GtoErrMsgFailure (routine);

    GtoFreeDateList (tl);
    return status;
}





/*
***************************************************************************
** Computes the non-contingent cash flows for a fee leg.
***************************************************************************
*/
TCashFlowList* CxFeeLegFlows
(CxTFeeLeg      *fl)
{
    static char routine[] = "CxFeeLegFlows";
    int         status    = FAILURE;

    TCashFlowList *cfl = NULL;
    int i;

    cfl = GtoNewEmptyCFL(fl->nbDates);
    if (cfl == NULL)
        goto done; /* failure */

    for (i = 0; i < fl->nbDates; ++i)
    {
        double amount;
        double time;

        time = CxDayCountFraction(fl->accStartDates[i],
                                  fl->accEndDates[i],
                                  fl->dcc);

        amount = time * fl->couponRates[i] * fl->notionals[i];

        cfl->fArray[i].fDate   = fl->payDates[i];
        cfl->fArray[i].fAmount = amount;
    }

    status = SUCCESS;

 done:

    if (status != SUCCESS)
    {
        GtoFreeCFL(cfl);
        cfl = NULL;
        GtoErrMsgFailure(routine);
    }
    return cfl;
}




/*
***************************************************************************
** Computes the expected cash flows for a fee leg. This returns the cash
** flows on the fee leg payment dates which are expected (taking account
** of default etc).
**
** The calculation involves computing the PV of each fee period, and then
** forward valuing that value to the cash flow payment date. Although this
** is not very sensitive to the discount curve, we therefore need a
** discount curve.
***************************************************************************
*/
TCashFlowList* CxFeeLegExpectedFlows
(CxTFeeLeg      *fl,
 TDate           today,
 TCurve         *discCurve,
 CxTCreditCurve *spreadCurve)
{
    static char routine[] = "CxFeeLegContingentFlows";
    int         status    = FAILURE;

    int    i;
    TCashFlowList *cfl = NULL;

    TDateList  *tl = NULL;

    REQUIRE (fl != NULL);
    REQUIRE (discCurve != NULL);
    REQUIRE (spreadCurve != NULL);

    cfl = GtoNewEmptyCFL(fl->nbDates);
    if (cfl == NULL)
        goto done; /* failure */

    if (fl->nbDates > 1)
    {
        /* it is more efficient to compute the timeLine just the once
           and truncate it for each payment */
        TDate startDate = fl->accStartDates[0];
        TDate endDate   = fl->accEndDates[fl->nbDates-1];
        
        tl = CxRiskyTimeLine (startDate,
                              endDate,
                              discCurve,
                              spreadCurve,
                              NULL, /* lossCurve */
                              0); /* delay */
        if (tl == NULL)
            goto done; /* failure */
    }

    for (i = 0; i < fl->nbDates; ++i)
    {
        double pv;
        double fv;
        double zeroPrice;

        if (FeePaymentPVWithTimeLine (fl->accrualPayConv,
                                      today,
                                      fl->accStartDates[i],
                                      fl->accEndDates[i],
                                      fl->payDates[i],
                                      fl->dcc,
                                      fl->notionals[i],
                                      fl->couponRates[i],
                                      discCurve,
                                      spreadCurve,
                                      tl,
                                      fl->obsStartOfDay,
                                      &pv) != SUCCESS)
            goto done; /* failure */

        zeroPrice = CxForwardZeroPrice (discCurve, today, fl->payDates[i]);
        fv = pv / zeroPrice;

        cfl->fArray[i].fDate   = fl->payDates[i];
        cfl->fArray[i].fAmount = fv;
    }

    status = SUCCESS;

 done:

    GtoFreeDateList (tl);
    if (status != SUCCESS)
    {
        GtoFreeCFL(cfl);
        cfl = NULL;
        GtoErrMsgFailure(routine);
    }
    return cfl;
}





