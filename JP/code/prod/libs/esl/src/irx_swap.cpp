#include "irx/swap.h"

#include <assert.h>

#include "irx/cfl.h"
#include "irx/rate.h"
#include "irx/zerocurve.h"
#include <irx/calendar.h>
#include <irx/datelist.h>
#include <irx/dateutils.h>
#include <irx/macros.h>

static int swapDates
(IrxTDate         startDate,      /* (I) Start date */
 IrxTDate         maturityDate,   /* (I) End date */
 IrxTDate         rollDate,       /* (I) Roll date */
 IrxTStubLocation stubLocation,   /* (I) Stub location */
 IrxTBadDayConv   badDayConv,     /* (I) Use the payment bad day convention */
 const IrxTCalendar *calendar,       
 IrxTDateInterval interval,       /* (I) Date interval */
 IrxTDateList   **unadjDates,     /* (O) Unadjusted coupon dates */
 IrxTDateList   **adjDates,       /* (O) Adjusted coupon dates */
 IrxTBool        *hasFrontStub,
 IrxTBool        *hasBackStub
);

static int guessRollDateFrontStub
(IrxTDate         startDate,
 IrxTDate         maturityDate,
 IrxTBadDayConv   badDayConv,
 const IrxTCalendar    *calendar,
 IrxTDateInterval interval,
 IrxTDate        *guessRollDate);

/**
***************************************************************************
** Makes a date list using swap-style conventions. Includes bad day
** adjustments, stub location and roll date.
**
** Two cases:
**   rollDate = 0  - all dates are on cycle with either startDate or
**                   maturityDate (depending on the stub location)
**   rollDate != 0 - all dates are on cycle with the rollDate. Then the
**                   dates are business day adjusted, and the resulting
**                   date list should include either the startDate or
**                   the maturityDate (depending on the stub location)
**
** The dates are bad day adjusted and represent the coupon payment dates.
***************************************************************************
*/
IrxTDateList* irxSwapPayDates
(IrxTDate         startDate,     /* (I) Start date */
 IrxTDate         maturityDate,  /* (I) End date */
 IrxTDate         rollDate,      /* (I) Roll date */
 IrxTStubLocation stubLocation,  /* (I) Stub location */
 IrxTBadDayConv   payBadDayConv, /* (I) Payment bad day convention */
 IrxTCalendar    *calendar,      /* (I) Holiday calendar */
 IrxTDateInterval interval)      /* (I) Date interval */
{
    static char routine[] = "irxSwapPayDates";
    int         status    = FAILURE;

    IrxTDateList *unadjDates   = NULL;
    IrxTDateList *adjDates     = NULL;
    IrxTDateList *swapPayDates = NULL;
    IrxTDate      adjMaturityDate;
    IrxTBool      hasFrontStub;
    IrxTBool      hasBackStub;

    if (irxBusinessDay (maturityDate,
                        payBadDayConv,
                        calendar,
                        &adjMaturityDate) != SUCCESS)
        goto RETURN; /* failure */

    if (swapDates (startDate,
                   maturityDate,
                   rollDate,
                   stubLocation,
                   payBadDayConv,
                   calendar,
                   interval,
                   &unadjDates,
                   &adjDates,
                   &hasFrontStub,
                   &hasBackStub) != SUCCESS)
        goto RETURN; /* failure */
    
    swapPayDates = irxDateListMake (adjDates->fNumItems-1,
                                    adjDates->fArray+1);

    if (swapPayDates == NULL)
        goto RETURN; /* failure */

    if ((swapPayDates->fNumItems < 1 ) ||
        (adjMaturityDate > swapPayDates->fArray[swapPayDates->fNumItems-1]))
    {
        PROGRAM_BUG();
        goto RETURN; /* failure */
    }
        
    swapPayDates->fArray[swapPayDates->fNumItems-1] = adjMaturityDate;
    status = SUCCESS;

 RETURN:

    irxDateListFree (adjDates);
    irxDateListFree (unadjDates);
    
    if (status != SUCCESS)
    {
        irxDateListFree (swapPayDates);
        swapPayDates = NULL;
        irxErrorFailure (routine);
    }
    return swapPayDates;
}


IrxTSwapPayments* irxSwapPaymentsGenerate (
    const IrxTSwap*     swap,
    double              fixedNotional,
    double              floatNotional,
    const IrxTCalendar* calendar,
    IrxTBool            includeFinalNotional,
    IrxTBool            includeInitialNotional,
    IrxTDate            firstFloatFixedDate,
    double              firstFloatFixedRate)
{
    static char routine[] = "irxSwapPaymentsGenerate";
    int         status    = FAILURE;

    IrxTSwapPayments  *swapPayments = NULL;
    IrxTDateList      *unadjFixedDates = NULL;
    IrxTDateList      *adjFixedDates = NULL;
    IrxTDateList      *unadjFloatDates = NULL;
    IrxTDateList      *adjFloatDates = NULL;
    IrxTDate           adjMaturityDate;
    IrxTBool           fixedHasFrontStub;
    IrxTBool           fixedHasBackStub;
    IrxTBool           floatHasFrontStub;
    IrxTBool           floatHasBackStub;
    IrxTBool           hasAccruedInterest;
    IrxTCouponPayment *payment;
    int                totalDates;
    int                i;
    int                j;

    REQUIRE (swap != NULL);
    REQUIRE (IS_NOT_ZERO(fixedNotional) || IS_NOT_ZERO(floatNotional));
   
    if (irxSwapMarketConvValidate (swap->marketConv) != SUCCESS)
        goto RETURN; /* failure */

    if (irxBusinessDay (swap->maturityDate,
                        swap->marketConv->paymentBdc,
                        calendar,
                        &adjMaturityDate) != SUCCESS)
        goto RETURN; /* failure */

    totalDates = 0;
    if (IS_NOT_ZERO(fixedNotional))
    {
        if (swapDates (swap->startDate,
                       swap->maturityDate,
                       swap->rollDate,
                       swap->stubLocation,
                       swap->marketConv->paymentBdc,
                       calendar,
                       swap->marketConv->fixedIvl,
                       &unadjFixedDates,
                       &adjFixedDates,
                       &fixedHasFrontStub,
                       &fixedHasBackStub) != SUCCESS)
            goto RETURN; /* failure */

        totalDates += adjFixedDates->fNumItems;
    }

    if (IS_NOT_ZERO(floatNotional))
    {
        if (swapDates (swap->startDate,
                       swap->maturityDate,
                       swap->rollDate,
                       swap->stubLocation,
                       swap->marketConv->paymentBdc,
                       calendar,
                       swap->marketConv->floatIvl,
                       &unadjFloatDates,
                       &adjFloatDates,
                       &floatHasFrontStub,
                       &floatHasBackStub) != SUCCESS)
            goto RETURN; /* failure */

        totalDates += adjFloatDates->fNumItems;
    }

    swapPayments = irxSwapPaymentsMakeEmpty (totalDates);

    /*
    ** Process fixed dates.
    **
    ** If there is a front stub, then we care whether we have
    ** a BOND stub (adds accrued interest), a NONE stub (adds no
    ** accrued interest) or a SIMPLE stub.
    **
    ** Ditto with a back stub, but in this case a BOND stub is
    ** not allowed.
    **
    ** If we want to add an initial notional, then we do that
    ** at the startDate. If we want to add a final notional,
    ** then we do that at the maturityDate.
    **
    ** If we have no accrued interest and no initial notional,
    ** then we can skip the startDate in the list of swap
    ** payments (in which case swapPayments will have been
    ** over allocated).
    **
    */

    j = 0;

    if (IS_NOT_ZERO(fixedNotional))
    {
        IrxTDate accrualDate;
        int      iStart; /* position in dateList of first payment date       */
        int      iEnd;   /* position in dateList of penultimate payment date */

        switch (swap->fixedStubPayment)
        {
        case IRX_STUB_NONE: /* full coupon */
        case IRX_STUB_BOND: /* full coupon - accrued interest */
            if (irxBusinessDay(unadjFixedDates->fArray[0],
                               swap->marketConv->accrualBdc,
                               calendar,
                               &accrualDate) != SUCCESS)
                goto RETURN; /* failure */
            break;
        case IRX_STUB_SIMPLE: /* part coupon */
            accrualDate = swap->startDate;
            break;
        default:
            PROGRAM_BUG();
            goto RETURN; /* failure */
        }
        
        hasAccruedInterest = fixedHasFrontStub && 
            (swap->fixedStubPayment == IRX_STUB_BOND);

        
        iStart = 1;
        if (fixedHasFrontStub &&
            swap->stubLocation == IRX_LONG_FRONT_STUB)
        {
            /* The dates generated include all on-cycle dates.    */
            /* This would mean we have a short stub in this case. */
            /* By starting at date 2 we enforce a long stub.      */
            iStart = 2; 
        }

        iEnd = adjFixedDates->fNumItems-2;
        if (fixedHasBackStub)
        {
            if (swap->fixedStubPayment != IRX_STUB_SIMPLE)
            {
                irxError ("%s: Only simple stubs are supported for back stubs",
                          routine);
                goto RETURN; /* failure */
            }
            if (swap->stubLocation == IRX_LONG_BACK_STUB)
            {
                /* The dates generated include all on-cycle dates. */
                /* With a back stub this would imply a short stub. */
                /* By ending at -3 we enforce a long stub.         */
                iEnd = adjFixedDates->fNumItems-3;
            }
        }

        if (hasAccruedInterest || includeInitialNotional)
        {
            /* Under these circumstances there are negative cash flows at
               the startDate. */
            if (hasAccruedInterest)
            {
                payment = swapPayments->couponPayments+j; /* notational ease */
                /* Use -notional rather than backwards interval */
                /* This is because dcc(a,b) not always equal -dcc(b,a) */
                payment->notional        = -fixedNotional;
                payment->couponRate      = swap->couponRate;
                payment->spread          = 0.0;
                payment->dcc             = swap->marketConv->fixedDcc;
                payment->accrueStartDate = accrualDate;
                payment->accrueEndDate   = swap->startDate;
            }
            if (includeInitialNotional)
                swapPayments->amounts[j] = -fixedNotional;

            swapPayments->dates[j] = swap->startDate;
            ++j;
        }

        for (i = iStart; i <= iEnd; ++i)
        {
            /* This is a coupon payment - excluding the final coupon */
            payment = swapPayments->couponPayments+j; /* notational ease */

            payment->notional         = fixedNotional;
            payment->couponRate       = swap->couponRate;
            payment->spread           = 0.0;
            payment->indexCurveNumber = 0;
            payment->dcc              = swap->marketConv->fixedDcc;
            payment->accrueStartDate  = accrualDate;
            if (irxBusinessDay (unadjFixedDates->fArray[i],
                                swap->marketConv->accrualBdc,
                                calendar,
                                &accrualDate) != SUCCESS)
                goto RETURN; /* failure */
            payment->accrueEndDate    = accrualDate;

            swapPayments->dates[j] = adjFixedDates->fArray[i];
            ++j;
        }

        /* Handle the final coupon separately */
        payment = swapPayments->couponPayments+j; /* notational ease */
        
        payment->notional         = fixedNotional;
        payment->couponRate       = swap->couponRate;
        payment->spread           = 0.0;
        payment->indexCurveNumber = 0;
        payment->dcc              = swap->marketConv->fixedDcc;
        payment->accrueStartDate  = accrualDate;
        if (irxBusinessDay (swap->maturityDate,
                            swap->marketConv->accrualBdc,
                            calendar,
                            &accrualDate) != SUCCESS)
            goto RETURN; /* failure */
        payment->accrueEndDate    = accrualDate;
        swapPayments->dates[j]         = adjMaturityDate;

        if (includeFinalNotional)
            swapPayments->amounts[j] += fixedNotional;

        ++j;
    }

    /*
    ** Process float dates.
    **
    ** If there is a front stub, then we care whether we have
    ** a BOND stub or a NONE stub as before. However we cannot
    ** calculate the accrued interest for the BOND stub unless
    ** it is a known fixed rate. Thus we would require for a
    ** non-trivial BOND stub that the firstFloatFixedDate was
    ** equal to the BOND stub pre-coupon date. We will therefore
    ** make the same requirement for a non-trivial NONE stub
    ** as well.
    **
    ** Even if there is no front stub, we may have the case that
    ** the firstFloatFixedDate kicks in to convert a floating
    ** rate into a fixed rate.
    **
    ** Otherwise - the initial notional and final notional are
    ** treated as for the fixed leg. Once again we may be able
    ** to skip the startDate in the list of swap payments.
    */
    if (IS_NOT_ZERO(floatNotional))
    {
        IrxTDate accrualDate;
        IrxTDate rateStartDate;
        IrxTDate rateEndDate;
        int      iStart; /* position in dateList of first payment date       */
        int      iEnd;   /* position in dateList of penultimate payment date */

        switch (swap->floatStubPayment)
        {
        case IRX_STUB_NONE: /* full coupon */
        case IRX_STUB_BOND: /* full coupon - accrued interest */
            rateStartDate = unadjFloatDates->fArray[0];
            if (irxBusinessDay(unadjFloatDates->fArray[0],
                               swap->marketConv->accrualBdc,
                               calendar,
                               &accrualDate) != SUCCESS)
                goto RETURN; /* failure */
            break;
        case IRX_STUB_SIMPLE: /* part coupon */
            rateStartDate = swap->startDate;
            accrualDate   = swap->startDate;
            break;
        default:
            PROGRAM_BUG();
            goto RETURN; /* failure */
        }
        
        hasAccruedInterest = floatHasFrontStub && 
            (swap->floatStubPayment == IRX_STUB_BOND);
        
        iStart = 1;
        if (floatHasFrontStub &&
            swap->stubLocation == IRX_LONG_FRONT_STUB)
        {
            /* The dates generated include all on-cycle dates.    */
            /* This would mean we have a short stub in this case. */
            /* By starting at date 2 we enforce a long stub.      */
            iStart = 2; 
        }

        iEnd = adjFloatDates->fNumItems-2;
        if (floatHasBackStub)
        {
            if (swap->floatStubPayment != IRX_STUB_SIMPLE)
            {
                irxError ("%s: Only simple stubs are supported for back stubs",
                          routine);
                goto RETURN; /* failure */
            }
            if (swap->stubLocation == IRX_LONG_BACK_STUB)
            {
                /* The dates generated include all on-cycle dates. */
                /* With a back stub this would imply a short stub. */
                /* By ending at -3 we enforce a long stub.         */
                iEnd = adjFloatDates->fNumItems-3;
            }
        }

        if (hasAccruedInterest || includeInitialNotional)
        {
            /* Under these circumstances there are negative cash flows at
               the startDate. */
            if (hasAccruedInterest)
            {
                payment = swapPayments->couponPayments+j; /* notational ease */
                payment->notional         = -floatNotional;
                payment->couponRate       = 0.0;
                payment->spread           = swap->spread;
                payment->indexCurveNumber = 1;
                payment->dcc              = swap->marketConv->floatDcc;
                payment->accrueStartDate  = accrualDate;
                payment->accrueEndDate    = swap->startDate;
                /* rateStartDate pre-calculated - there are special cases
                   for initial coupon which are handled by calculating
                   rateStartDate beforehand */
                if (irxBusinessDay (rateStartDate,
                                    swap->marketConv->paymentBdc,
                                    calendar,
                                    &payment->rateStartDate) != SUCCESS)
                    goto RETURN; /* failure */

                if (rateStartDate == unadjFloatDates->fArray[0])
                {
                    /* rateStartDate is on cycle - so add the normal period
                       (after adjusting for bad days) */
                    if (irxDateAddInterval (payment->rateStartDate,
                                            swap->marketConv->floatIvl,
                                            &rateEndDate) != SUCCESS)
                        goto RETURN; /* failure */
                }
                else
                {
                    /* rateStart is not on cycle - use unusual period */
                    rateEndDate = unadjFloatDates->fArray[iStart];
                }
                if (irxBusinessDay (rateEndDate,
                                    swap->marketConv->paymentBdc,
                                    calendar,
                                    &payment->rateEndDate) != SUCCESS)
                    goto RETURN; /* failure */
            }
            if (includeInitialNotional)
                swapPayments->amounts[j] = -floatNotional;

            swapPayments->dates[j] = swap->startDate;
            ++j;
        }

        for (i = iStart; i <= iEnd; ++i)
        {
            /* This is a coupon payment - excluding the final coupon */
            payment = swapPayments->couponPayments+j; /* notational ease */

            payment->notional         = floatNotional;
            payment->couponRate       = 0.0;
            payment->spread           = 0.0;
            payment->indexCurveNumber = 1;
            payment->dcc              = swap->marketConv->floatDcc;
            payment->accrueStartDate  = accrualDate;
            if (irxBusinessDay (unadjFloatDates->fArray[i],
                                swap->marketConv->accrualBdc,
                                calendar,
                                &accrualDate) != SUCCESS)
                goto RETURN; /* failure */
            payment->accrueEndDate    = accrualDate;
            /* rateStartDate pre-calculated - there are special cases
               for initial coupon which are handled by calculating
               rateStartDate beforehand */
            if (irxBusinessDay (rateStartDate,
                                swap->marketConv->paymentBdc,
                                calendar,
                                &payment->rateStartDate) != SUCCESS)
                goto RETURN; /* failure */
            
            if (rateStartDate == unadjFloatDates->fArray[i-1])
            {
                /* rateStartDate is on cycle - so add the normal period
                   (after adjusting for bad days) */
                if (irxDateAddInterval (payment->rateStartDate,
                                        swap->marketConv->floatIvl,
                                        &rateEndDate) != SUCCESS)
                    goto RETURN; /* failure */
            }
            else
            {
                /* rateStart is not on cycle - use unusual period */
                rateEndDate = unadjFloatDates->fArray[iStart];
            }
            if (irxBusinessDay (rateEndDate,
                                swap->marketConv->paymentBdc,
                                calendar,
                                &payment->rateEndDate) != SUCCESS)
                goto RETURN; /* failure */

            rateStartDate = unadjFloatDates->fArray[i]; /* for next coupon */
            swapPayments->dates[j] = adjFloatDates->fArray[i];
            ++j;
        }

        /* Handle the final coupon separately */
        payment = swapPayments->couponPayments+j; /* notational ease */
        
        payment->notional         = floatNotional;
        payment->couponRate       = 0.0;
        payment->spread           = 0.0;
        payment->indexCurveNumber = 1;
        payment->dcc              = swap->marketConv->floatDcc;
        payment->accrueStartDate  = accrualDate;
        if (irxBusinessDay (swap->maturityDate,
                            swap->marketConv->accrualBdc,
                            calendar,
                            &accrualDate) != SUCCESS)
            goto RETURN; /* failure */

        payment->accrueEndDate    = accrualDate;

        /* rateStartDate pre-calculated - there are special cases
           for initial coupon which are handled by calculating
           rateStartDate beforehand */
        if (irxBusinessDay (rateStartDate,
                            swap->marketConv->paymentBdc,
                            calendar,
                            &payment->rateStartDate) != SUCCESS)
            goto RETURN; /* failure */
            
        if (floatHasBackStub)
        {
            /* dates are not on cycle - so use the maturity date */
            rateEndDate = swap->maturityDate;
        }
        else
        {
            /* dates are on cycle - so add the normal period
               (after adjusting for bad days) */
            if (irxDateAddInterval (payment->rateStartDate,
                                    swap->marketConv->floatIvl,
                                    &rateEndDate) != SUCCESS)
                goto RETURN; /* failure */
        }
        if (irxBusinessDay (rateEndDate,
                            swap->marketConv->paymentBdc,
                            calendar,
                            &payment->rateEndDate) != SUCCESS)
            goto RETURN; /* failure */

        swapPayments->dates[j] = adjMaturityDate;
        if (includeFinalNotional)
            swapPayments->amounts[j] += floatNotional;

        ++j;
    }

    assert(j <= swapPayments->numFlows);
    swapPayments->numFlows = j;

    status = SUCCESS;

 RETURN:

    irxDateListFree(unadjFixedDates);
    irxDateListFree(adjFixedDates);
    irxDateListFree(unadjFloatDates);
    irxDateListFree(adjFloatDates);
    if (status != SUCCESS)
    {
        irxSwapPaymentsFree (swapPayments);
        swapPayments = NULL;
        irxErrorFailure (routine);
    }

    return swapPayments;
}


/*
***************************************************************************
** Returns swap dates - both unadjusted and adjusted for payment bad day
** convention. The first date in the list is the start of the first
** coupon period. There should therefore be one more date in the list
** than actual payments in the swap.
**
** The act of adjusting for bad days may remove stubs at the end of the
** set of coupons, and thus give us one less date in the datelist.
** However if for accrual we do not adjust for bad days, then we would
** be stuffed if we didn't also returns the unadjusted dates.
**
** The first date will be on or before the start date.
** The last date in the adjusted list will be on or after the adjusted
** maturity date.
***************************************************************************
*/
static int swapDates
(IrxTDate         startDate,     /* (I) Start date */
 IrxTDate         maturityDate,  /* (I) End date */
 IrxTDate         rollDate,      /* (I) Roll date */
 IrxTStubLocation stubLocation,  /* (I) Stub location */
 IrxTBadDayConv   payBadDayConv,
 const IrxTCalendar *calendar,
 IrxTDateInterval interval,
 IrxTDateList   **pUnadjDates,
 IrxTDateList   **pAdjDates,
 IrxTBool        *hasFrontStub,
 IrxTBool        *hasBackStub)
{
    static char routine[] = "swapDates";
    int         status    = FAILURE;

    IrxTDateList *unadjDates = NULL;
    IrxTDateList *adjDates   = NULL;
    IrxTDate      adjMaturityDate;
    IrxTDate      adjStartDate;
    int           i;

    char dateBuf[2][16];

    REQUIRE (maturityDate > startDate);

    if (irxBusinessDay (maturityDate, payBadDayConv, calendar,
                        &adjMaturityDate) != SUCCESS)
        goto RETURN; /* failure */

    if (irxBusinessDay (startDate, payBadDayConv, calendar,
                        &adjStartDate) != SUCCESS)
        goto RETURN; /* failure */

    if (startDate != adjStartDate)
    {
        irxError ("%s: Start date %s is not a business day\n", routine,
                  irxDateFormat(startDate, "DD-MMM-YYYY", dateBuf[0]));
        goto RETURN; /* failure */
    }

    if (rollDate == 0)
    {
        if (irxStubLocationAtFront(stubLocation))
        {
            if (guessRollDateFrontStub (startDate,
                                        maturityDate,
                                        payBadDayConv,
                                        calendar,
                                        interval,
                                        &rollDate) != SUCCESS)
                goto RETURN; /* failure */
        }
        else
        {
            rollDate = startDate;
        }
    }

    unadjDates = irxDateListMakeWithRoll (startDate,
                                          maturityDate,
                                          rollDate,
                                          interval);
    if (unadjDates == NULL)
        goto RETURN; /* failure */

    adjDates = irxDateListMakeEmpty (unadjDates->fNumItems);
    for (i = 0; i < unadjDates->fNumItems; ++i)
    {
        if (irxBusinessDay (unadjDates->fArray[i],
                            payBadDayConv,
                            calendar,
                            &adjDates->fArray[i]) != SUCCESS)
            goto RETURN; /* failure */
    }

    /* imagine we are near the end of the month
       then it is highly possible that unadjDates[N-2] < maturityDate
       whereas adjDates[N-2] = maturityDate
       
       in this case we have one more date than we really need */

    if (adjDates->fNumItems < 2)
    {
        /* startDate < maturityDate means this cannot happen */
        PROGRAM_BUG();
        goto RETURN; /* failure */
    }

    if (adjDates->fArray[adjDates->fNumItems-2] == adjMaturityDate)
    {
        --(unadjDates->fNumItems);
        --(adjDates->fNumItems);
    }

    if (adjDates->fArray[0] > startDate ||
        adjDates->fArray[adjDates->fNumItems-1] < adjMaturityDate)
    {
        PROGRAM_BUG();
        goto RETURN; /* failure */
    }

    *hasFrontStub = adjDates->fArray[0] < startDate;
    *hasBackStub  = adjDates->fArray[adjDates->fNumItems-1] > adjMaturityDate;
    *pUnadjDates  = unadjDates;
    *pAdjDates    = adjDates;

    unadjDates = NULL;
    adjDates   = NULL;

    status = SUCCESS;

 RETURN:

    irxDateListFree (adjDates);
    irxDateListFree (unadjDates);

    if (status != SUCCESS)
        irxErrorFailure (routine);
    return status;
}

/*
 * Guesses the roll date when we have a front stub.
 *
 * Assumption is that startDate is not bad day adjusted, but that
 * maturityDate might be bad day adjusted.
 *
 * For complete accuracy in pricing swaps you should provide the roll
 * date as a parameter. When building a zero curve, the extra roll date
 * parameter is a bit of a pain and is liable not to be provided.
 *
 * Does the following:
 *
 *    If maturityDate to startDate is exact number of date intervals,
 *    then it assumes maturityDate is rollDate.
 *
 *    If startDate to maturityDate is exact number of date intervals,
 *    then it assumes startDate is rollDate.
 *
 *    If startDate to maturityDate is not an exact number of date
 *    intervals, but when you adjust the last period for bad days you
 *    hit the maturity date, then it assumes that the start date is
 *    the roll date.
 *
 *    If startDate,maturityDate is nowhere near an exact number of date
 *    intervals and the date interval is monthly, then it will use
 *    the maturityDate.
 */
static int guessRollDateFrontStub
(IrxTDate         startDate,
 IrxTDate         maturityDate,
 IrxTBadDayConv   badDayConv,
 const IrxTCalendar    *calendar,
 IrxTDateInterval interval,
 IrxTDate        *guessRollDate)
{
    static char routine[] = "guessRollDateFrontStub";
    int         status    = FAILURE;

    int numIntervals;
    int remainder;

    IrxTDate         rollDate;
    IrxTDate         aDate;
    int              i;

    REQUIRE (interval.prd > 0);
    REQUIRE (maturityDate > startDate);

    /* First try - if there is no stub from the viewpoint of the maturityDate
       then clearly we can use the maturityDate as the rollDate */
    if (irxCountDateIntervals (maturityDate,
                               startDate,
                               interval,
                               &numIntervals,
                               &remainder) != SUCCESS)
        goto RETURN; /* failure */

    if (remainder == 0)
    {
        rollDate = maturityDate;
        goto success;
    }

    /* Second try - if there is no stub from the viewpoint of the startDate
       then we can use the startDate as the rollDate - slightly less clear */
    if (irxCountDateIntervals (startDate,
                               maturityDate,
                               interval,
                               &numIntervals,
                               &remainder) != SUCCESS)
        goto RETURN; /* failure */

    if (remainder == 0)
    {
        rollDate = startDate;
        goto success;
    }

    /* We try to use the startDate and find that the maturity date arises
       by bad day adjustment of the startDate implied maturity date, then
       we can use the startDate as the rollDate - again slightly less clear */
    for (i = numIntervals; i <= numIntervals+1; ++i)
    {
        if (irxDateAddMultiInterval (startDate,
                                     i,
                                     interval,
                                     &aDate) != SUCCESS)
            goto RETURN; /* failure */
        
        if (irxBusinessDay (aDate, badDayConv, calendar, &aDate) != SUCCESS)
            goto RETURN; /* failure */
        
        if (aDate == maturityDate)
        {
            rollDate = startDate;
            goto success;
        }
    }

    rollDate = maturityDate;

 success:

    *guessRollDate = rollDate;
    status         = SUCCESS;

 RETURN:

    /* Bad status arises from bad inputs - not a failed algorithm */
    if (status != SUCCESS)
        irxErrorFailure (routine);
        
    return status;
}



/*
 * Computes the last relevant date for a given curve for a set of swap
 * payments.
 */
int irxSwapPaymentsLastDate (
    const IrxTSwapPayments  *swapPayments,
    int                      curveNumber,
    IrxTDate                *lastDate)
{
    static char routine[] = "irxSwapPaymentsLastDate";
    int         status    = FAILURE;

    IrxTDate myLastDate = 0;
    int      i;

    REQUIRE (swapPayments != NULL);
    REQUIRE (curveNumber >= 0);
    REQUIRE (lastDate != NULL);

    if (curveNumber == 0)
    {
        for (i = 0; i < swapPayments->numFlows; ++i)
        {
            if (swapPayments->dates[i] > myLastDate)
                myLastDate = swapPayments->dates[i];
        }
    }
    else
    {
        for (i = 0; i < swapPayments->numFlows; ++i)
        {
            IrxTCouponPayment *payments = swapPayments->couponPayments+i;
            if (payments->indexCurveNumber == curveNumber &&
                payments->rateEndDate > myLastDate)
            {
                myLastDate = payments->rateEndDate;
            }
        }
    }

    *lastDate = myLastDate;
    status = SUCCESS;
    
 RETURN:

    if (status != SUCCESS)
        irxErrorFailure (routine);

    return status;
}





/*
 * Converts swap payments into a cash flow list using one (or more) index
 * curves. Normally one index curve is enough.
 *
 * Cash flows on the same date are netted. A net value of zero remains in
 * the cash flow list.
 */
IrxTCashFlowList* irxSwapPaymentsFlows (
    const IrxTSwapPayments *swapPayments,
    int                     numIndexCurves,
    const IrxTZeroCurve   **indexCurves,
    IrxTDate                minDate)
{
    static char routine[] = "irxSwapPaymentsFlows";
    int         status    = FAILURE;

    IrxTCashFlowList *cfl = NULL;
    IrxTCashFlowList *tmp;
    int               i;
    int               j;

    REQUIRE (swapPayments != NULL);
    REQUIRE (numIndexCurves >= 0);
    if (numIndexCurves > 0)
        REQUIRE (indexCurves != NULL);

    cfl = irxCashFlowListMakeEmpty(swapPayments->numFlows);
    if (cfl == NULL)
        goto RETURN; /* failure */

    j = 0;
    for (i = 0; i < swapPayments->numFlows; ++i)
    {
        double             amount  = swapPayments->amounts[i];
        IrxTCouponPayment *payment = swapPayments->couponPayments+i;
        IrxTDate           date    = swapPayments->dates[i];

        if (date < minDate) 
            continue;
        
        if (IS_NOT_ZERO(payment->notional))
        {
            double couponRate;
            double accrueTime;
            double couponAmount;

            if (payment->indexCurveNumber > 0)
            {
                const IrxTZeroCurve *indexCurve;

                REQUIRE(payment->indexCurveNumber <= numIndexCurves);
                indexCurve = indexCurves[payment->indexCurveNumber-1];
                REQUIRE(indexCurve != NULL);
                REQUIRE (payment->rateStartDate >= indexCurve->baseDate);

                if (irxFwdZeroRate (indexCurve,
                                    payment->rateStartDate,
                                    payment->rateEndDate,
                                    payment->dcc,
                                    IRX_SIMPLE_RATE,
                                    &couponRate) != SUCCESS)
                    goto RETURN; /* failure */
            }
            else
            {
                couponRate = payment->couponRate;
            }
            
            if (irxDayCountFraction (payment->accrueStartDate,
                                     payment->accrueEndDate,
                                     payment->dcc,
                                     &accrueTime) != SUCCESS)
                goto RETURN; /* failure */
            
            couponAmount = (couponRate + payment->spread) * accrueTime;
            
            amount += couponAmount * payment->notional;
        }
        cfl->dates[j]   = date;
        cfl->amounts[j] = amount;
        ++j;
    }
    cfl->numItems = j;

    tmp = irxCashFlowListCopySort (cfl);
    if (tmp == NULL)
        goto RETURN; /* failure */

    irxCashFlowListFree (cfl);
    cfl = tmp;

    status = SUCCESS;

 RETURN:

    if (status != SUCCESS)
    {
        irxCashFlowListFree (cfl);
        cfl = NULL;
        irxErrorFailure (routine);
    }

    return cfl;
}



/**
 * Computes the fixed payments for a swap.
 */
IrxTCashFlowList* irxSwapFixedFlows (
    const IrxTSwap     *swap,
    double        notional,
    IrxTDate      valueDate,
    const IrxTCalendar *calendar,
    IrxTBool      includeFinalNotional,
    IrxTBool      includeInitialNotional)
{
    static char routine[] = "irxSwapFixedFlows";
    int         status    = FAILURE;

    IrxTCashFlowList *cfl = NULL;
    IrxTSwapPayments *swapPayments = NULL;
    IrxTSwap         *copySwap = NULL;

    REQUIRE (swap != NULL);

    /* We handle a stub payment of BOND by changing the start date of the
       swap. Obviously we only do this on a copy of the swap. */
    if (swap->fixedStubPayment == IRX_STUB_BOND && valueDate > swap->startDate)
    {
        copySwap = irxSwapCopy (swap);
        if (copySwap == NULL)
            goto RETURN; /* failure */
        if (copySwap->rollDate == 0)
            copySwap->rollDate = copySwap->startDate;
        copySwap->startDate = valueDate;
        swap = copySwap;
    }
    
    swapPayments = irxSwapPaymentsGenerate (swap,
                                            notional,
                                            0.0,
                                            calendar,
                                            includeFinalNotional,
                                            includeInitialNotional,
                                            FALSE,
                                            0.0);
    if (swapPayments == NULL)
        goto RETURN; /* failure */

    cfl = irxSwapPaymentsFlows (swapPayments,
                                0,
                                NULL,
                                valueDate);
    if (cfl == NULL)
        goto RETURN; /* failure */

    status = SUCCESS;

 RETURN:

    irxSwapFree (copySwap);
    irxSwapPaymentsFree (swapPayments);

    if (status != SUCCESS)
    {
        irxCashFlowListFree (cfl);
        cfl = NULL;
        irxErrorFailure (routine);
    }

    return cfl;
}


/**
 * Computes the floating payments for a swap.
 */
IrxTCashFlowList* irxSwapFloatFlows (
    const IrxTSwap      *swap,
    const IrxTZeroCurve *indexCurve,
    double         notional,
    IrxTDate       valueDate,
    const IrxTCalendar  *calendar,
    IrxTBool       includeFinalNotional,
    IrxTBool       includeInitialNotional,
    IrxTDate       firstFloatFixedDate,
    double         firstFloatFixedRate)
{
    static char routine[] = "irxSwapFixedFlows";
    int         status    = FAILURE;

    IrxTCashFlowList *cfl = NULL;
    IrxTSwapPayments *swapPayments = NULL;
    IrxTSwap         *copySwap = NULL;

    REQUIRE (swap != NULL);

    /* We handle a stub payment of BOND by changing the start date of the
       swap. Obviously we only do this on a copy of the swap. */
    if (swap->floatStubPayment == IRX_STUB_BOND && valueDate > swap->startDate)
    {
        copySwap = irxSwapCopy (swap);
        if (copySwap == NULL)
            goto RETURN; /* failure */
        if (copySwap->rollDate == 0)
            copySwap->rollDate = copySwap->startDate;
        copySwap->startDate = valueDate;
        swap = copySwap;
    }
    
    swapPayments = irxSwapPaymentsGenerate (swap,
                                            0.0,
                                            notional,
                                            calendar,
                                            includeFinalNotional,
                                            includeInitialNotional,
                                            firstFloatFixedDate,
                                            firstFloatFixedRate);
    if (swapPayments == NULL)
        goto RETURN; /* failure */

    cfl = irxSwapPaymentsFlows (swapPayments,
                                1,
                                &indexCurve,
                                valueDate);
    if (cfl == NULL)
        goto RETURN; /* failure */

    status = SUCCESS;

 RETURN:

    irxSwapFree (copySwap);
    irxSwapPaymentsFree (swapPayments);

    if (status != SUCCESS)
    {
        irxCashFlowListFree (cfl);
        cfl = NULL;
        irxErrorFailure (routine);
    }

    return cfl;
}





/**
 * Computes the par rate for a swap. The coupon within the swap definition
 * is ignored.
 */
int irxSwapRate (
    const IrxTSwap      *swap,
    const IrxTZeroCurve *discountCurve,
    IrxTBool       valueFloating,
    const IrxTZeroCurve *indexCurve,
    double         pv,
    const IrxTCalendar  *calendar,
    IrxTDate       firstFloatFixedDate,
    double         firstFloatFixedRate,
    double        *swapRate)
{
    static char routine[] = "irxSwapRate";
    int         status    = FAILURE;

    IrxTCashFlowList *cflFloat = NULL;
    IrxTCashFlowList *cflFixed = NULL;
    double            pvAnnuity;
    double            pvTarget;
    IrxTSwap         *mySwap = NULL;
    IrxTDate          valueDate;

    REQUIRE (swap != NULL);

    mySwap = irxSwapCopy (swap);
    if (mySwap == NULL)
        goto RETURN; /* failure */
    mySwap->couponRate = 1.0;
    valueDate = mySwap->startDate;

    cflFixed = irxSwapFixedFlows (mySwap,
                                  1.0, /* notional */
                                  valueDate,
                                  calendar,
                                  FALSE, /* includeFinalNotional */
                                  FALSE); /* includeInitialNotional */
    if (cflFixed == NULL)
        goto RETURN; /* failure */

    if (cflFixed->numItems <= 0)
    {
        irxError ("%s: No cash flows for swap\n", routine);
        goto RETURN; /* failure */
    }

    if (irxCashFlowListFV (cflFixed,
                           discountCurve,
                           valueDate,
                           &pvAnnuity) != SUCCESS)
        goto RETURN; /* failure */

    /* It is hard to imagine how pvAnnuity == 0.0 given that the coupon
       was 1.0 and we have checked that the number of flows is not zero */
    if (IS_ZERO(pvAnnuity))
    {
        PROGRAM_BUG();
        goto RETURN; /* failure */
    }
    
    if (valueFloating)
    {
        cflFloat = irxSwapFloatFlows (mySwap,
                                      indexCurve,
                                      1.0, /* notional */
                                      valueDate,
                                      calendar,
                                      FALSE,
                                      FALSE,
                                      firstFloatFixedDate,
                                      firstFloatFixedRate);
        if (irxCashFlowListFV (cflFloat,
                               discountCurve,
                               valueDate,
                               &pvTarget) != SUCCESS)
            goto RETURN; /* failure */
    }
    else
    {
        /* pvTarget = pv - pvLastDate */
        /* since flows are sorted - the last date is easy to obtain */
        /* it might be different from mySwap->maturityDate because of
           bad day adjustments */
        IrxTDate lastDate;
        double   pvLastDate;
        lastDate = cflFixed->dates[cflFixed->numItems-1];
        
        if (irxFwdZeroPrice (discountCurve,
                             valueDate,
                             lastDate,
                             &pvLastDate) != SUCCESS)
            goto RETURN; /* failure */
        
        pvTarget = pv-pvLastDate;
    }

    *swapRate = pvTarget / pvAnnuity;
    status = SUCCESS;

 RETURN:

    irxCashFlowListFree (cflFloat);
    irxCashFlowListFree (cflFixed);
    irxSwapFree (mySwap);

    if (status != SUCCESS)
        irxErrorFailure (routine);

    return status;
}



/**
 * Validates market conventions as regards the swap market.
 *
 * We have the following rules:
 *
 *   The coupon intervals must be positive.
 *
 *   The bad day conventions must be consistent. Bad day conventions
 *   can either be trivial or consistently non-trivial. If payment
 *   bad day convention is trivial, then accrual bad day convention
 *   must also be trivial.
 */
int irxSwapMarketConvValidate (
    const IrxTMarketConv* marketConv)
{
    static char routine[] = "irxSwapMarketConvValidate";
    int         status    = FAILURE;

    IrxTBadDayConv nonTrivial = IRX_BAD_DAY_NONE;
    IrxTBadDayConv bdcList[3];
    int i;

    REQUIRE (marketConv != NULL);
    REQUIRE (marketConv->fixedIvl.prd > 0);
    REQUIRE (marketConv->floatIvl.prd > 0);
    REQUIRE (marketConv->daysToSpot >= 0);

    bdcList[0] = marketConv->paymentBdc;
    bdcList[1] = marketConv->accrualBdc;
    bdcList[2] = marketConv->resetBdc;

    for (i = 0; i < 3; ++i)
    {
        if (bdcList[i] != IRX_BAD_DAY_NONE)
        {
            if (nonTrivial == IRX_BAD_DAY_NONE)
            {
                nonTrivial = bdcList[i];
            }
            else if (bdcList[i] != nonTrivial)
            {
                irxError ("%s: Swap bad day conventions must be consistently "
                          "non-trivial. There is a mixture of %s and %s\n",
                          routine,
                          irxBadDayConvToString(nonTrivial),
                          irxBadDayConvToString(bdcList[i]));
                goto RETURN; /* failure */
            }
        }
    }

    if (marketConv->paymentBdc == IRX_BAD_DAY_NONE &&
        marketConv->accrualBdc != IRX_BAD_DAY_NONE)
    {
        irxError ("%s: Cannot have no payment bad day convention and "
                  "non-trivial accrued bad day convention %s\n",
                  routine,
                  irxBadDayConvToString(marketConv->accrualBdc));
        goto RETURN; /* failure */
    }

    status = SUCCESS;

 RETURN:

    if (status != SUCCESS)
        irxErrorFailure (routine);

    return status;
}


/**
 * Splits a swap payments stream into flows which are completely
 * on or before a given date, and flows which are after that date.
 *
 * Definition of the last critical date for an individual flow is the
 * maximum of the maturity of the effective date and the payment date.
 */
int irxSwapPaymentsSplit (
    const IrxTSwapPayments *payments,
    IrxTDate                splitDate,
    IrxTSwapPayments      **before,
    IrxTSwapPayments      **after)
{
    static char routine[] = "irxSwapPaymentsSplit";
    int         status    = FAILURE;

    int i,j,k;
    
    IrxTSwapPayments *myBefore = NULL;
    IrxTSwapPayments *myAfter  = NULL;

    REQUIRE (payments != NULL);
    REQUIRE (before != NULL);
    REQUIRE (after != NULL);

    myBefore = irxSwapPaymentsMakeEmpty (payments->numFlows);
    myAfter  = irxSwapPaymentsMakeEmpty (payments->numFlows);
    if (myBefore == NULL || myAfter == NULL)
        goto RETURN; /* failure */

    j = 0;
    k = 0;
    for (i = 0; i < payments->numFlows; ++i)
    {
        IrxTDate maxDate = payments->dates[i];
        
        if (payments->couponPayments[i].indexCurveNumber > 0)
        {
            if (payments->couponPayments[i].rateEndDate > maxDate)
                maxDate = payments->couponPayments[i].rateEndDate;
        }

        if (payments->dates[i] <= splitDate)
        {
            myBefore->dates[j]          = payments->dates[i];
            myBefore->amounts[j]        = payments->amounts[i];
            myBefore->couponPayments[j] = payments->couponPayments[i];
            ++j;
        }
        else
        {
            myAfter->dates[k]          = payments->dates[i];
            myAfter->amounts[k]        = payments->amounts[i];
            myAfter->couponPayments[k] = payments->couponPayments[i];
            ++k;
        }
    }
    myBefore->numFlows = j;
    myAfter->numFlows  = k;

    *before = myBefore;
    *after  = myAfter;

    myBefore = NULL;
    myAfter  = NULL;
    
    status = SUCCESS;

 RETURN:

    irxSwapPaymentsFree (myAfter);
    irxSwapPaymentsFree (myBefore);

    if (status != SUCCESS)
        irxErrorFailure (routine);

    return status;
}
