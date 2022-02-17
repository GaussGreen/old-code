/*
***************************************************************************
** FILE NAME: cds.h
**
** Vanilla CDS functions.
**
** $Header$
***************************************************************************
*/

#ifndef CX_CDS_H
#define CX_CDS_H

#include "cx.h"

/**
 * Makes a contingent leg for a vanilla CDS.
 *
 * The CDS starts at startDate and ends at endDate. The last date is always
 * protected - the start date is only protected if protectStart=True.
 *
 * Maximum compatibility with CMLib is obtained by setting protectStart=False.
 */
CxTContingentLeg* CxCdsContingentLegMake(
    /** Date when protection begins. Either at start or end of day (depends
        on protectStart) */
    TDate     startDate,
    /** Date when protection ends (end of day) */
    TDate     endDate,
    /** Notional value protected */
    double    notional, 
    /** Delay in contingent payment after default (in days) */
    long      delay, 
    /** Should protection include the start date */
    TBoolean  protectStart);

/**
 * Computes the PV for a contingent leg for a vanilla CDS.
 *
 * Risk starts at the end of today. The PV is computed for a given value date.
 * The CDS starts at startDate and ends at endDate. The last date is always
 * protected - the start date is only protected if protectStart=True.
 *
 * Maximum compatibility with CMLib is obtained by setting protectStart=False.
 */
int CxCdsContingentLegPV(
    /** Risk starts at the end of today */
    TDate           today,
    /** Date for which the PV is calculated */
    TDate           valueDate,
    /** Date when protection begins. Either at start or end of day (depends
        on protectStart) */
    TDate     startDate,
    /** Date when protection ends (end of day) */
    TDate     endDate,
    /** Notional value protected */
    double    notional, 
    /** Delay in contingent payment after default (in days) */
    long      delay, 
    /** Interest rate discount curve - assumes flat forward interpolation */
    TCurve         *discCurve,
    /** Credit clean spread curve */
    CxTCreditCurve *spreadCurve,
    /** Assumed recovery curve in case of default */
    CxTRecoveryCurve *recoveryCurve,
    /** True => protection includes start date */
    TBoolean        protectStart,
    /** Output - the present value is returned */
    double         *pv);

/**
 * Makes a fixed fee leg for a vanilla CDS.
 *
 * The CDS starts at startDate and ends at endDate. The last date is always
 * protected - the start date is only protected if protectStart=True.
 *
 * Interest accrues for the same number of days as there is protection.
 * Thus if protectStart=True you get one extra day of accrued interest in
 * comparison with an interest rate swap. This extra day is assumed to be
 * the last day of the CDS and means that the last period is one day longer
 * than for an interest rate swap.
 *
 * Maximum compatibility with CMLib is obtained by setting protectStart=False.
 */
CxTFeeLeg* CxCdsFeeLegMake(
    /** Date when protection begins. Either at start or end of day (depends
        on protectStart) */
    TDate           startDate,
    /** Date when protection ends (end of day) */
    TDate           endDate,
    /** Should accrued interest be paid on default. Usually set to TRUE */
    TBoolean        payAccOnDefault,
    /** Interval between coupon payments. Can be NULL when 3M is assumed */
    TDateInterval  *couponInterval,
    /** If the startDate and endDate are not on cycle, then this parameter
        determines location of coupon dates. Normal is CX_SHORT_FRONT_STUB. */
    CxTStubType     stubType,
    /** Notional value protected */
    double          notional, 
    /** Fixed coupon rate (a.k.a. spread) for the fee leg */
    double          couponRate,
    /** Day count convention for coupon payment. Normal is CX_ACT_360 */
    CxTDayCountConv paymentDcc,
    /** Bad day convention for adjusting coupon payment dates. */
    CxTBadDayConv   badDayConv,
    /** Calendar used when adjusting coupon dates. Can be NULL which equals
        a calendar with no holidays and including weekends. */
    CxTCalendar    *calendar,
    /** Should protection include the start date */
    TBoolean  protectStart);

/**
 * Computes the PV for a fixed fee leg for a vanilla CDS.
 *
 * Risk starts at the end of today. The PV is computed for a given value date.
 * The CDS starts at startDate and ends at endDate. The last date is always
 * protected - the start date is only protected if protectStart=True.
 *
 * Interest accrues for the same number of days as there is protection.
 * Thus if protectStart=True you get one extra day of accrued interest in
 * comparison with an interest rate swap. This extra day is assumed to be
 * the last day of the CDS and means that the last period is one day longer
 * than for an interest rate swap.
 *
 * Maximum compatibility with CMLib is obtained by setting protectStart=False.
 */
int CxCdsFeeLegPV(
    /** Risk starts at the end of today */
    TDate           today,
    /** Date for which the PV is calculated */
    TDate           valueDate,
    /** Date when protection begins. Either at start or end of day (depends
        on protectStart) */
    TDate           startDate,
    /** Date when protection ends (end of day) */
    TDate           endDate,
    /** Should accrued interest be paid on default. Usually set to TRUE */
    TBoolean        payAccOnDefault,
    /** Interval between coupon payments. Can be NULL when 3M is assumed */
    TDateInterval  *couponInterval,
    /** If the startDate and endDate are not on cycle, then this parameter
        determines location of coupon dates. Normal is CX_SHORT_FRONT_STUB. */
    CxTStubType     stubType,
    /** Notional value protected */
    double          notional, 
    /** Fixed coupon rate (a.k.a. spread) for the fee leg */
    double          couponRate,
    /** Day count convention for coupon payment. Normal is CX_ACT_360 */
    CxTDayCountConv paymentDcc,
    /** Bad day convention for adjusting coupon payment dates. */
    CxTBadDayConv   badDayConv,
    /** Calendar used when adjusting coupon dates. Can be NULL which equals
        a calendar with no holidays and including weekends. */
    CxTCalendar    *calendar,
    /** Interest rate discount curve - assumes flat forward interpolation */
    TCurve         *discCurve,
    /** Credit clean spread curve */
    CxTCreditCurve *spreadCurve,
    /** Should protection include the start date */
    TBoolean        protectStart,
    /** Should the present value be computed as a clean price (removing
        accrued interest) */
    TBoolean        isPriceClean,
    /** Output - the present value is returned */
    double         *pv);
  
/**
 * Computes the par spread for a vanilla CDS
 *
 * Risk starts at the end of today. The PV is computed for a given value date.
 * The CDS starts at startDate and ends at endDate. The last date is always
 * protected - the start date is only protected if protectStart=True.
 *
 * Interest accrues for the same number of days as there is protection.
 * Thus if protectStart=True you get one extra day of accrued interest in
 * comparison with an interest rate swap. This extra day is assumed to be
 * the last day of the CDS and means that the last period is one day longer
 * than for an interest rate swap.
 *
 * Maximum compatibility with CMLib is obtained by setting protectStart=False.
 */
int CxCdsParSpread(
    /** Risk starts at the end of today */
    TDate           today,
    /** Date for which the PV is calculated */
    TDate           valueDate,
    /** Date when protection begins. Either at start or end of day (depends
        on protectStart) */
    TDate           startDate,
    /** Date when protection ends (end of day) */
    TDate           endDate,
    /** Delay in contingent payment after default (in days) */
    long            delay, 
    /** Price (a.k.a. upfront charge) for the CDS (see also isPriceClean) */
    double          price,
    /** Should accrued interest be paid on default. Usually set to TRUE */
    TBoolean        payAccOnDefault,
    /** Interval between coupon payments. Can be NULL when 3M is assumed */
    TDateInterval  *couponInterval,
    /** If the startDate and endDate are not on cycle, then this parameter
        determines location of coupon dates. Normal is CX_SHORT_FRONT_STUB. */
    CxTStubType     stubType,
    /** Day count convention for coupon payment. Normal is CX_ACT_360 */
    CxTDayCountConv paymentDcc,
    /** Bad day convention for adjusting coupon payment dates. */
    CxTBadDayConv   badDayConv,
    /** Calendar used when adjusting coupon dates. Can be NULL which equals
        a calendar with no holidays and including weekends. */
    CxTCalendar    *calendar,
    /** Interest rate discount curve - assumes flat forward interpolation */
    TCurve         *discCurve,
    /** Credit clean spread curve */
    CxTCreditCurve *spreadCurve,
    /** Assumed recovery curve in case of default */
    CxTRecoveryCurve *recoveryCurve,
    /** True => protection includes start date */
    TBoolean        protectStart,
    /** Is the price expressed as a clean price (removing accrued interest) */
    TBoolean        isPriceClean,
    /** Output - the par spread is returned */
    double         *parSpread);

/**
 * Computes the price (a.k.a. upfront charge) for a vanilla CDS
 *
 * Risk starts at the end of today. The PV is computed for a given value date.
 * The CDS starts at startDate and ends at endDate. The last date is always
 * protected - the start date is only protected if protectStart=True.
 *
 * Interest accrues for the same number of days as there is protection.
 * Thus if protectStart=True you get one extra day of accrued interest in
 * comparison with an interest rate swap. This extra day is assumed to be
 * the last day of the CDS and means that the last period is one day longer
 * than for an interest rate swap.
 *
 * Maximum compatibility with CMLib is obtained by setting protectStart=False.
 */
int CxCdsPrice(
    /** Risk starts at the end of today */
    TDate           today,
    /** Date for which the PV is calculated */
    TDate           valueDate,
    /** Date when protection begins. Either at start or end of day (depends
        on protectStart) */
    TDate           startDate,
    /** Date when protection ends (end of day) */
    TDate           endDate,
    /** Delay in contingent payment after default (in days) */
    long            delay, 
    /** Fixed coupon rate (a.k.a. spread) for the fee leg */
    double          couponRate,
    /** Should accrued interest be paid on default. Usually set to TRUE */
    TBoolean        payAccOnDefault,
    /** Interval between coupon payments. Can be NULL when 3M is assumed */
    TDateInterval  *couponInterval,
    /** If the startDate and endDate are not on cycle, then this parameter
        determines location of coupon dates. Normal is CX_SHORT_FRONT_STUB. */
    CxTStubType     stubType,
    /** Day count convention for coupon payment. Normal is CX_ACT_360 */
    CxTDayCountConv paymentDcc,
    /** Bad day convention for adjusting coupon payment dates. */
    CxTBadDayConv   badDayConv,
    /** Calendar used when adjusting coupon dates. Can be NULL which equals
        a calendar with no holidays and including weekends. */
    CxTCalendar    *calendar,
    /** Interest rate discount curve - assumes flat forward interpolation */
    TCurve         *discCurve,
    /** Credit clean spread curve */
    CxTCreditCurve *spreadCurve,
    /** Assumed recovery curve in case of default */
    CxTRecoveryCurve *recoveryCurve,
    /** True => protection includes start date */
    TBoolean        protectStart,
    /** Is the price expressed as a clean price (removing accrued interest) */
    TBoolean        isPriceClean,
    /** Output - price (a.k.a. upfront charge) for the CDS is returned 
        (see also isPriceClean) */
    double         *price);











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
 TBoolean        protectStart);

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
 TBoolean        protectStart);

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
 TBoolean          protectStart);

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
 TBoolean          protectStart);

#endif

