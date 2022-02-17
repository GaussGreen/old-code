#ifndef IRX_SWAP_H
#define IRX_SWAP_H

#include "irxflow.h"

#ifdef __cplusplus
extern "C"
{
#endif

/**
***************************************************************************
** Makes a date list using swap-style conventions. Includes bad day
** adjustments, stub location and roll date.
**
** Two cases:
**   rollDate = 0  - all dates are on cycle with either startDate or
**                   maturityDate (depending on the stub location)
**   rollDate != 0 - all dates are on cycle with the rollDate
**
** The dates are bad day adjusted and represent the coupon payment dates.
** The last date in the list is always the bad day adjusted maturityDate.
** The first date in the list is after the startDate.
** If payBadDayConv is non-trivial, then startDate must be a business day.
***************************************************************************
*/
    IrxTDateList* irxSwapPayDates (
        IrxTDate         startDate,     /* (I) Start date */
        IrxTDate         maturityDate,  /* (I) End date */
        IrxTDate         rollDate,      /* (I) Roll date */
        IrxTStubLocation stubLocation,  /* (I) Stub location */
        IrxTBadDayConv   payBadDayConv, /* (I) Payment bad day convention */
        IrxTCalendar    *calendar,      /* (I) Holiday calendar */
        IrxTDateInterval interval);     /* (I) Date interval */


/**
 * Generates swap payments for a swap.
 *
 * Note you can generate flows for both fixed and floating leg, and there will
 * be no netting in this case.
 */
    IrxTSwapPayments* irxSwapPaymentsGenerate (
        /** Swap object - defines market conventions, coupon, spread etc. */
        const IrxTSwap*     swap,
        /** Notional of fixed leg. If zero, then no fixed payments. */
        double        fixedNotional,
        /** Notional of floating leg. If zero, then no floating payments. */
        double        floatNotional,
        /** Calendar for adjusting for bad days. */
        const IrxTCalendar* calendar,
        /** Include final notional payments */
        IrxTBool      includeFinalNotional,
        /** Include initial notional payments */
        IrxTBool      includeInitialNotional,
        /** Is the first floating rate fixed?
            In this case provide a non-zero date which represents the
            date for which the rate was fixed. */
        IrxTDate      firstFloatFixedDate,
        /** If the first floating rate is fixed then what is it? */
        double        firstFloatFixedRate);

/**
 * Computes the last relevant date for a given curve for a set of swap
 * payments.
 */
    int irxSwapPaymentsLastDate (
        const IrxTSwapPayments  *swapPayments,
        /** Use 0 for the discount curve. This will give you the last payment
            date in the swap payments.
            
            Use >0 for an index curve. This will give you the last rate
            maturity date for the relevant index curve.
        */
        int                curveNumber,
        IrxTDate          *lastDate);


/**
 * Converts swap payments into a cash flow list using one (or more) index
 * curves. Normally one index curve is enough.
 *
 * Cash flows on the same date are netted. A net value of zero remains in
 * the cash flow list.
 */
    IrxTCashFlowList* irxSwapPaymentsFlows (
        const IrxTSwapPayments *swapPayments,
        int               numIndexCurves,
        /* array of size numIndexCurves */
        IrxTZeroCurve const  **indexCurves,
        IrxTDate          minDate);



/**
 * Computes the fixed payments for a swap.
 */
    IrxTCashFlowList* irxSwapFixedFlows (
        /** Definition of the swap */
        const IrxTSwap  *swap,
        /** Notional value. */
        double        notional,
        /** Value date - any flows before the value date are removed.
            Also if we have a bond stub, then the accrued interest is 
            calculated for the value date and will appear as a negative
            coupon in the cash flow list. */
        IrxTDate      valueDate,
        /** Calendar for adjusting the cash flows */
        const IrxTCalendar *calendar,
        /** Include the final notional payment (+1) */
        IrxTBool      includeFinalNotional,
        /** Include the initial notional payment (-1) */
        IrxTBool      includeInitialNotional);


/**
 * Computes the floating payments for a swap.
 */
    IrxTCashFlowList* irxSwapFloatFlows (
        /** Definition of the swap */
        const IrxTSwap      *swap,
        /** Index curve for estimating the swap payments */
        const IrxTZeroCurve *indexCurve,
        /** Notional value. */
        double         notional,
        /** Value date - any flows before the value date are removed.
            Also if we have a bond stub, then the accrued interest is 
            calculated for the value date and will appear as a negative
            coupon in the cash flow list. */
        IrxTDate       valueDate,
        /** Calendar for adjusting the cash flows */
        const IrxTCalendar  *calendar,
        /** Include the final notional payment (+1) */
        IrxTBool       includeFinalNotional,
        /** Include the initial notional payment (-1) */
        IrxTBool       includeInitialNotional,
        /** Is the first floating rate fixed?
            In this case provide a non-zero date which represents the
            date for which the rate was fixed. */
        IrxTDate      firstFloatFixedDate,
        /** If the first floating rate is fixed then what is it? */
        double        firstFloatFixedRate);



/**
 * Computes the par rate for a swap. The coupon within the swap definition
 * is ignored.
 */
    int irxSwapRate (
        /** Definition of the swap - the coupon is ignored */
        const IrxTSwap      *swap,
        /** Discount curve for discounting swap payments. */
        const IrxTZeroCurve *discountCurve,
        /** Should we value the floating leg? */
        IrxTBool       valueFloating,
        /** Index curve for estimating the swap payments.
            Ignored if we don't value the floating leg. */
        const IrxTZeroCurve *indexCurve,
        /** PV of the fixed leg we are trying to match. This is computed 
            as of the startDate of the swap. You need a good reason not
            to provide 1.0.
            Ignored if we value the floating leg. */
        double         pv,
        /** Calendar for adjusting the cash flows */
        const IrxTCalendar  *calendar,
        /** Is the first floating rate fixed?
            In this case provide a non-zero date which represents the
            date for which the rate was fixed. */
        IrxTDate      firstFloatFixedDate,
        /** If the first floating rate is fixed then what is it? */
        double        firstFloatFixedRate,
        /** Swap rate returned */
        double       *swapRate);

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
        const IrxTMarketConv* marketConv);


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
        IrxTSwapPayments      **after);

#ifdef __cplusplus
}
#endif

#endif
