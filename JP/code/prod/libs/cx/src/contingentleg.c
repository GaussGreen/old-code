/*
***************************************************************************
** FILE NAME: contingentLeg.c
**
** Analytics for a contingent leg.
**
** $Header$
***************************************************************************
*/

#include "contingentleg.h"

#include <math.h>

#include <cxutils/include/alibconv.h>
#include "creditcurve.h"
#include "crxconv.h"
#include "recovery.h"
#include "timeline.h"

#include <cxutils/include/cxmacros.h>
#include <cxutils/include/datelist.h>
#include <cxutils/include/zerocurve.h>

#include <crxflow/include/crcrv.h>

static int ContingentLegPV
(CxTContingentLeg *cl,               /* (I) Contingent leg                  */
 TDate             today,            /* (I) No observations before today    */
 TDate             valueDate,        /* (I) Value date for discounting      */
 TCurve           *discountCurve,    /* (I) Risk-free curve                 */
 CxTCreditCurve   *spreadCurve,      /* (I) Spread curve                    */
 CxTRecoveryCurve *recoveryCurve,    /* (I) Recovery curve                  */
 double           *pv);              /* (O) Present value of contingent leg */

static int onePeriodIntegralWithDelay
(TDate             today,
 TDate             startDate,
 TDate             endDate,
 long              delay,
 TCurve           *discountCurve,
 CxTCreditCurve   *spreadCurve,
 CxTRecoveryCurve *recoveryCurve,
 double           *pv);

static int onePeriodIntegralAtPayDate
(TDate             today,
 TDate             startDate,
 TDate             endDate,
 TDate             payDate, 
 TCurve           *discountCurve,
 CxTCreditCurve   *spreadCurve,
 CxTRecoveryCurve *recoveryCurve,
 double           *pv);

/**
 * Computes the PV of a contingent leg as a whole
 */
int CxContingentLegPV
(CxTContingentLeg *cl,               /* (I) Contingent leg                  */
 TDate             today,            /* (I) No observations before today    */
 TDate             valueDate,        /* (I) Value date for discounting      */
 TCurve           *discountCurve,    /* (I) Risk-free curve                 */
 CxTCreditCurve   *spreadCurve,      /* (I) Spread curve                    */
 CxTRecoveryCurve *recoveryCurve,    /* (I) Recovery curve                  */
 double           *pv)               /* (O) Present value of contingent leg */
{
    static char routine[] = "CxContingentLegPV";
    int         status    = FAILURE;

    REQUIRE (spreadCurve != NULL);
    
    switch (spreadCurve->type)
    {
    case CX_CURVE_TYPE_FLOW:
        status = ContingentLegPV (cl,
                                  today,
                                  valueDate,
                                  discountCurve,
                                  spreadCurve,
                                  recoveryCurve,
                                  pv);
        break;
    case CX_CURVE_TYPE_EXOTIC:
    {
        KProtLeg_D *pl = CxProtectionLegConvert (cl, recoveryCurve,
                                                 spreadCurve->timestep);
        if (pl != NULL)
        {
            status = ProtectionPV_O (pv,
                                     today,
                                     valueDate,
                                     pl,
                                     discountCurve,
                                     spreadCurve->tc);
        }
        CrxProtectionFree (pl);
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

/**
 * Computes the PV of a contingent leg as a whole
 */
int ContingentLegPV
(CxTContingentLeg *cl,               /* (I) Contingent leg                  */
 TDate             today,            /* (I) No observations before today    */
 TDate             valueDate,        /* (I) Value date for discounting      */
 TCurve           *discountCurve,    /* (I) Risk-free curve                 */
 CxTCreditCurve   *spreadCurve,      /* (I) Spread curve                    */
 CxTRecoveryCurve *recoveryCurve,    /* (I) Recovery curve                  */
 double          *pv)                /* (O) Present value of contingent leg */
{
    static char routine[] = "ContingentLegPV";
    int         status    = FAILURE;

    int    i;
    double myPv = 0.0;
    double valueDatePv;

    TDate startDate;

    REQUIRE (cl != NULL);
    REQUIRE (discountCurve != NULL);
    REQUIRE (spreadCurve != NULL);
    REQUIRE (pv != NULL);
    REQUIRE (spreadCurve->type = CX_CURVE_TYPE_FLOW);

    startDate = cl->startDate;

    switch (cl->payType)
    {
    case CX_PROT_PAY_MAT:
        for (i = 0; i < cl->nbDates; ++i)
        {
            double tmp;
            
            if (onePeriodIntegralAtPayDate (today,
                                            startDate,
                                            cl->dates[i],
                                            cl->dates[i] + cl->payDelay,
                                            discountCurve,
                                            spreadCurve,
                                            recoveryCurve,
                                            &tmp) != SUCCESS)
                goto done; /* failure */
            
            startDate = cl->dates[i];
            myPv += tmp * cl->notionals[i];
        }
        break;
    case CX_PROT_PAY_DEF:
        for (i = 0; i < cl->nbDates; ++i)
        {
            double tmp;
            
            if (onePeriodIntegralWithDelay (today,
                                            startDate,
                                            cl->dates[i],
                                            cl->payDelay,
                                            discountCurve,
                                            spreadCurve,
                                            recoveryCurve,
                                            &tmp) != SUCCESS)
                goto done; /* failure */
            
            startDate = cl->dates[i];
            myPv += tmp * cl->notionals[i];
        }
        break;
    default:
        GtoErrMsg ("%s: Invalid payment type %d\n", routine, (int)cl->payType);
        goto done; /* failure */
    }

    /* myPv has been calculated as at today - need it at valueDate */
    valueDatePv = CxForwardZeroPrice (discountCurve, today, valueDate);

    status = SUCCESS;
    *pv    = myPv / valueDatePv;

 done:

    if (status != SUCCESS)
        GtoErrMsgFailure (routine);

    return status;
}

/**
 * Computes the PV of a single period of a contingent leg.
 *
 * Calling this function repeatedly with the correct sensible inputs would
 * give the same result as CxContingentLegPV on the whole contingent leg.
 */
int CxContingentPaymentPV
(CxTProtPayConv    payType,          /* (I) Contingent leg delay type       */
 TDate             today,            /* (I) No observations before today    */
 double            notional,         /* (I) Notional value                  */
 TDate             startDate,        /* (I) Observation start date          */
 TDate             endDate,          /* (I) Observation end date            */
 long              payDelay,         /* (I) Delay in payment                */
 TCurve           *discountCurve,    /* (I) Risk-free curve                 */
 CxTCreditCurve   *spreadCurve,      /* (I) Spread curve                    */
 CxTRecoveryCurve *recoveryCurve,    /* (I) Recovery curve                  */
 double           *pv)               /* (O) Present value of contingent leg */
{
    static char routine[] = "CxContingentPaymentPV";
    int         status    = FAILURE;

    REQUIRE (spreadCurve != NULL);
    REQUIRE (spreadCurve->type == CX_CURVE_TYPE_FLOW); /* TBD */

    switch (payType)
    {
    case CX_PROT_PAY_MAT:
        if (onePeriodIntegralAtPayDate (today,
                                        startDate,
                                        endDate,
                                        endDate + payDelay,
                                        discountCurve,
                                        spreadCurve,
                                        recoveryCurve,
                                        pv) != SUCCESS)
            goto done; /* failure */
        break;
    case CX_PROT_PAY_DEF:
        if (onePeriodIntegralWithDelay (today,
                                        startDate,
                                        endDate,
                                        payDelay,
                                        discountCurve,
                                        spreadCurve,
                                        recoveryCurve,
                                        pv) != SUCCESS)
            goto done; /* failure */
        break;
    default:
        GtoErrMsg ("%s: Invalid payment type %d\n", routine, (int)payType);
        goto done; /* failure */
    }

    *pv   *= notional;
    status = SUCCESS;

 done:

    if (status != SUCCESS)
        GtoErrMsgFailure (routine);

    return status;
}

/*
 * Computes a one period integral with a delay period.
 */
static int onePeriodIntegralWithDelay
(TDate             today,
 TDate             startDate,
 TDate             endDate,
 long              delay,
 TCurve           *discCurve,
 CxTCreditCurve   *spreadCurve,
 CxTRecoveryCurve *recoveryCurve,
 double           *pv)
{
    static char routine[] = "onePeriodIntegralAtPayDate";
    int         status    = FAILURE;

    double  myPv = 0.0;
    int     i;

    double t;
    double s0;
    double s1;
    double df0;
    double df1;
    double loss0;
    double loss1;
    double recovery;

    TDateList *tl = NULL;

    REQUIRE (endDate > startDate);
    REQUIRE (discCurve != NULL);
    REQUIRE (spreadCurve != NULL);
    REQUIRE (pv != NULL);
    REQUIRE (delay >= 0);

    if (today >= endDate)
    {
        myPv = 0.0;
        goto success;
    }

    if (today > startDate)
        startDate = today;

    tl = CxRiskyTimeLine (startDate, endDate, discCurve, spreadCurve, 
                          recoveryCurve, delay);
    if (tl == NULL) goto done; /* failure */

    /* the integration - we can assume flat forwards between points on
       the timeline - this is true for both curves 

       we are integrating -Z dS/dt where Z is the discount factor and
       S is the survival probability

       assuming flat forwards on each part of the integration, this is an
       exact integral

       we then make some adjustments to account for the fact that in reality
       we are sampling only once per business day, which on average over the
       year means there is about one day of delay in the integral
    */

    s1  = CxForwardZeroPrice (spreadCurve->tc, today, startDate);
    df1 = CxForwardZeroPrice (discCurve, today, startDate+delay);
    if (recoveryCurve != NULL)
    {
        if (CxRecoveryCurveInterp (recoveryCurve, startDate, 
                                   &recovery) != SUCCESS)
            goto done; /* failure */
        loss1 = 1.0 - recovery;
    }
    else
    {
        loss1 = 1.0;
    }

    for (i = 1; i < tl->fNumItems; ++i)
    {
        double lambda;
        double fwdRate;
        double thisPv;

        /* adj and adjDays are to do with sampling frequency error */
        double adj;
        double adjDays=0.0; /* 0.5 for all days + weekends + holidays */
        
        s0  = s1;
        df0 = df1;
        s1  = CxForwardZeroPrice (spreadCurve->tc, today, tl->fArray[i]);
        df1 = CxForwardZeroPrice (discCurve, today, tl->fArray[i]+delay);
        loss0 = loss1;
        if (recoveryCurve != NULL)
        {
            if (CxRecoveryCurveInterp(recoveryCurve, tl->fArray[i], 
                                      &recovery) != SUCCESS)
                goto done; /* failure */
            loss1 = 1.0 - recovery;
        }
        else
        {
            loss1 = 1.0;
        }
        t   = (double)(tl->fArray[i] - tl->fArray[i-1])/365.0;
        
        lambda  = log(s0/s1)/t;
        fwdRate = log(df0/df1)/t;
        
        adj     = exp(-fwdRate*adjDays/365.0);
        thisPv  = loss0 * lambda / (lambda + fwdRate) * 
            (1.0 - exp(-(lambda + fwdRate) * t)) * s0 * df0;
        
        if (IS_NOT_EQUAL(loss1,loss0))
        {
            switch (recoveryCurve->interpType)
            {
            case CX_RECOVERY_LINEAR_INTERP:
                /* linear interpolation is the tricky case */
                thisPv += (loss1-loss0)/t * lambda * s0 * df0 *
                    (1.0/(lambda+fwdRate) - exp(-(lambda+fwdRate)*t) *
                     (t + 1.0/(lambda+fwdRate))) / (lambda+fwdRate);
                break;
            case CX_RECOVERY_RIGHT_INTERP:
                /* we should have used loss1 instead of loss 0 */
                thisPv *= (loss1 / loss0); 
                break;
            case CX_RECOVERY_LEFT_INTERP:
                break; /* do nothing - the formula was correct */
            }
        }

        myPv   += thisPv * adj;
    }

 success:

    status = SUCCESS;
    *pv = myPv;
        
 done:

    if (status != SUCCESS)
        GtoErrMsgFailure (routine);

    GtoFreeDateList (tl);
    return status;
}


/*
 * Computes a one period integral with payment at a specific payment date.
 *
 * This is actually trivial.
 */
static int onePeriodIntegralAtPayDate
(TDate             today,
 TDate             startDate,
 TDate             endDate,
 TDate             payDate, 
 TCurve           *discCurve,
 CxTCreditCurve   *spreadCurve,
 CxTRecoveryCurve *recoveryCurve,
 double           *pv)
{
    static char routine[] = "onePeriodIntegralAtPayDate";
    int         status    = FAILURE;

    double df;
    double s0;
    double s1;
    double loss;

    REQUIRE (endDate > startDate);
    REQUIRE (discCurve != NULL);
    REQUIRE (spreadCurve != NULL);
    REQUIRE (pv != NULL);

    if (today >= endDate)
    {
        *pv = 0.0;
    }
    else
    {
        if (today > startDate)
            startDate = today;

        s0  = CxForwardZeroPrice (spreadCurve->tc, today, startDate);
        s1  = CxForwardZeroPrice (spreadCurve->tc, today, endDate);
        df  = CxForwardZeroPrice (discCurve, today, payDate);
        if (recoveryCurve != NULL)
        {
            double recovery;
            if (CxRecoveryCurveAverage(recoveryCurve, startDate, endDate, 
                                       &recovery) != SUCCESS)
                goto done; /* failure */
            loss = 1.0 - recovery;
        }
        else
        {
            loss = 1.0;
        }
        *pv = (s0 - s1) * df * loss;
    }

    status = SUCCESS;
        
 done:

    if (status != SUCCESS)
        GtoErrMsgFailure (routine);
    
    return status;
}




/*
***************************************************************************
** Computes the expected cash flows for a contingent leg. This returns the
** cash flows on the corresponding payment dates.
**
** The calculation involves computing the PV of each fee period, and then
** forward valuing that value to the fee cash flow payment date. Although
** this is not very sensitive to the discount curve, we therefore need
** a discount curve. We also need a fee leg to provide the dates.
***************************************************************************
*/
TCashFlowList* CxContingentLegExpectedFlows
(CxTContingentLeg *cl,
 TDateList        *flowDates,
 TDate             today,
 TCurve           *discCurve,
 CxTCreditCurve   *spreadCurve,
 CxTRecoveryCurve *recoveryCurve)
{
    static char routine[] = "CxContingentLegExpectedFlows";
    int         status    = SUCCESS;
    
    int i;
    int j;
    
    TDate startDate = cl->startDate;

    TCashFlowList *cfl = NULL;

    REQUIRE (cl != NULL);
    REQUIRE (flowDates != NULL);
    REQUIRE (discCurve != NULL);
    REQUIRE (spreadCurve != NULL);
    REQUIRE (recoveryCurve != NULL);


    cfl = GtoNewEmptyCFL(flowDates->fNumItems);
    if (cfl == NULL)
        goto done; /* failure */

    for (j = 0; j < flowDates->fNumItems; ++j)
    {
        cfl->fArray[j].fDate = flowDates->fArray[j];
    }

    for (i = 0; i < cl->nbDates; ++i)
    {
        TDate endDate = cl->dates[i];
        for (j = 0; j < flowDates->fNumItems; ++j)
        {
            TDate payDate = flowDates->fArray[j];
            TDate thisStartDate;
            TDate thisEndDate;
            double pv;
            double zeroPrice;
            double fv;

            if (payDate <= startDate)
                continue;

            if (j > 0)
            {
                thisStartDate = flowDates->fArray[j-1];
                if (thisStartDate < startDate)
                    thisStartDate = startDate;
            }
            else
            {
                thisStartDate = startDate;
            }

            if (payDate > endDate)
            {
                thisEndDate = endDate;
            }
            else
            {
                thisEndDate = payDate;
            }

            if (thisEndDate <= thisStartDate)
                continue;

            ASSERT(thisStartDate >= startDate);
            ASSERT(thisEndDate <= endDate);
            ASSERT(thisEndDate > thisStartDate);

            if (CxContingentPaymentPV (cl->payType,
                                       today,
                                       cl->notionals[i],
                                       thisStartDate,
                                       thisEndDate,
                                       cl->payDelay,
                                       discCurve,
                                       spreadCurve,
                                       recoveryCurve,
                                       &pv) != SUCCESS)
                goto done; /* failure */

            zeroPrice = CxForwardZeroPrice (discCurve, today, payDate);
            fv = pv / zeroPrice;

            cfl->fArray[j].fAmount += fv;
        }

        startDate = endDate;
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


