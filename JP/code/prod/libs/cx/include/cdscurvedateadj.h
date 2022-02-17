/*
***************************************************************************
** HEADER FILE: cdscurvedateadj.h
**
** Header file for function which converts the market spreads of a CDS
** curve into spreads on a different set of dates such that when you
** bootstrap from the new set of dates you back-out the results on the
** market dates.
**
** $Header$
***************************************************************************
*/

#ifndef CX_CDSCURVEDATEADJ_H
#define CX_CDSCURVEDATEADJ_H

#include "cx.h"

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
    /** Risk starts at the end of today */
    TDate           today,
    /** Interest rate discount curve - assumes flat forward interpolation */
    TCurve         *discCurve,
    /** Date when protection begins. Either at start or end of day (depends
        on protectStart) */
    TDate           startDate,
    /** Date for which the PV is calculated (where relevant) */
    TDate           valueDate,
    /** Number of benchmark dates */
    long            nbDate,
    /** Dates when protection ends for each benchmark (end of day).
        Array of size nbDate */
    TDate          *endDates,
    /** Coupon rates for each benchmark instrument. Array of size nbDate */
    double         *couponRates,
    /** Price (a.k.a. upfront charge) for each benchmark instrument. See also
        isPriceClean. Can be NULL if there are no upfront charges. Otherwise
        an array of size nbDate. */
    double         *prices, 
    /** Flags to denote that we include particular benchmarks. This makes it
        easy for the user to include or exclude benchmarks on a one-by-one
        basis. Can be NULL if all are included. Otherwise an array of size
        nbDate. */
    TBoolean       *includes,
    /** Dates for which we want the adjusted spreads and/or prices.
        Array of size nbDate. */
    TDate          *adjEndDates,
    /** Assumed recovery rate in case of default. */
    CxTRecoveryCurve *recoveryCurve,
    /** Should accrued interest be paid on default. Usually set to TRUE */
    TBoolean        payAccOnDefault,
    /** Interval between coupon payments. Can be NULL when 3M is assumed */
    TDateInterval  *couponInterval,
    /** Day count convention for coupon payment. Normal is CX_ACT_360 */
    CxTDayCountConv paymentDcc,
    /** If the startDate and endDate are not on cycle, then this parameter
        determines location of coupon dates. Normal is CX_SHORT_FRONT_STUB. */
    CxTStubType     stubType,
    /** Type of curve used when building the intermediate bootstrap curves.
        Use CX_CURVE_TYPE_FLOW for standard methodologies. 
        Use CX_CURVE_TYPE_EXOTIC for compatibility with exotic credit product
        pricing algorithms. */
    CxTCreditCurveType curveType,
    /** Integration timestep. NULL = use couponInterval */
    TDateInterval  *timestep,   
    /** True => protection includes start date */
    TBoolean        protectStart,
    /** Is the price expressed as a clean price (removing accrued interest) */
    TBoolean        isPriceClean,
    /** Delay in contingent payment after default (in days) */
    long            delay,
    /** Bad day convention for adjusting coupon payment dates. */
    CxTBadDayConv   badDayConv,
    /** Calendar used when adjusting coupon dates. Can be NULL which equals
        a calendar with no holidays and including weekends. */
    CxTCalendar    *calendar,
    /** Compute the tweak matrix. Note that the tweak matrix cannot be
      computed if the benchmark instruments are not at par (i.e. prices
      is defined). */
    TBoolean        computeTweakMatrixFlag
);

#endif
