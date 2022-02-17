//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : SwapTool.hpp
//
//   Description : All the nasty date generation routines you need for swaps
//
//   Author      : Andrew J Swain
//
//   Date        : 7 February 2001
//
//
//----------------------------------------------------------------------------

#ifndef SWAPTOOL_HPP
#define SWAPTOOL_HPP
#include <string>
#include "edginc/Maths.hpp"
#include "edginc/Format.hpp"
#include "edginc/DateTime.hpp"
#include "edginc/MaturityPeriod.hpp"
#include "edginc/CashFlow.hpp"
#include "edginc/DayCountConvention.hpp"
#include "edginc/Stub.hpp"
#include "edginc/BadDayConvention.hpp"
#include "edginc/Holiday.hpp"
#include "edginc/StubSimple.hpp"

using namespace std;    // string

DRLIB_BEGIN_NAMESPACE

/** All the nasty date generation routines you need for swaps */

class IRVolBase;
class StubPlacement;
class ZeroCurve;
class IYieldCurve;


class MARKET_DLL SwapTool {
public:

/** Calculates a par swap rate using passed zero curve for the swap starting at startDate,
 * maturing at maturityDate, with fixed day count convention 
 * fixedDayCountConv and fixed payments occuring at time 
 * intervals defined by interval. In other words, the routine calculates 
 * the fixed rate such that the present value of the fixed and 
 * floating sides are equal. The floating side is assumed to be at par.
 */
/* static */
static double couponRate(
    const DateTime&           startDate,       // (I) Date instrument begins at
    const DateTime&           maturityDate,    // (I) Date instrument matures at
    const MaturityPeriod&     interval,        // (I) Time between payments 
    bool                      stubAtEnd,
    const DayCountConvention* dcc,
    const ZeroCurve&          zc);

static double calculateCoupon(
    const ZeroCurve&     zc, 
    const CashFlowArray& cfl, 
    double               presentValue);

static double swapFixedPV(
    const ZeroCurve&          zeroCurve,
    const double              couponRate,
    const DateTime&           startDate,
    const MaturityPeriod&     couponInterval,
    const DateTime&           maturityDate,
    const DayCountConvention& dcc,
    const Stub&               stubType,         // SIMPLE, BOND
    const StubPlacement&      stubPlacement,    // stub is always short + simple
    const bool                subtractInitial,  // notional payment?
    const bool                addFinal,         // notional payment?
    const BadDayConvention&   accBadDayConv,    // accrual bad day convention
    const BadDayConvention&   payBadDayConv,    // payment bad day convention
    const Holiday&            holidays,
    const DateTime&           valueDate);

static double swapFloatPV(
    const ZeroCurve&          discounting,
    const ZeroCurve&          estimating,
    const double              spread,
    const bool                isAdditive,
    const DateTime&           startDate,
    const MaturityPeriod&     interval,
    const DateTime&           maturityDate,
    const DayCountConvention& dcc,
    const Stub&               stubType,
    const StubPlacement&      stubPlacement,    // stub is always short + simple
    const bool                subtractInitial,  // notional payment?
    const bool                addFinal,         // notional payment?
    const BadDayConvention&   accBadDayConv,    // accrual bad day convention
    const BadDayConvention&   payBadDayConv,    // payment bad day convention
    const BadDayConvention&   resetBadDayConv,  // reset bad day convention
    const Holiday&            holidays,
    const bool                firstFixed,
    const double              firstFixRate,
    const DateTime&           valueDate);

/**
 * Calculates a par swap rate for the swap starting at startDate, maturing
 * at endDate, with fixed day count convention fixedDcc and fixed
 * payments occurring at time intervals defined by fixedIvl. In other words, 
 * the routine calculates the fixed rate such that the present value of the
 * fixed and floating sides are equal.
 *
 * Takes account of holidays, and the possibility that the floating side is
 * not at par and cost of basis.
 *
 * The convexity/delay adjustment parameter is currently not supported,
 * and must be set to FALSE.
 *
 * The function is copied from the ALIB GtoSwapRate2 function.
 */
static double swapRate(
    const ZeroCurve&          discountCurve, 
    const DateTime&           startDate,
    const DateTime&           maturityDate,
    const MaturityPeriod&     fixedIvl,
    const DayCountConvention& fixedDayCountConv,
    const bool                valueFloating,
    const double              floatPV,
    const ZeroCurve*          estimatingCurve,
    const MaturityPeriod*     floatIvl,
    const DayCountConvention* floatDayCountConv,
    const bool                firstFixed,
    const double              firstFixRate,
    const bool                convexityDelayAdj,
    const IRVolBase*          volModelIR,
    const StubPlacement&      stubPos,
    const BadDayConvention&   accBadDayConv,    // accrual bad day convention
    const BadDayConvention&   payBadDayConv,    // payment bad day convention
    const BadDayConvention&   resetBadDayConv,  // reset bad day convention
    const Holiday&            holidays);

/**
 * Calculates a par swap rate,  but also includes a roll date input.
 * This date, with is in most cases identical to the start date is used in 
 * generating the cashFlow dates. In cases where there are front stubs, and the
 * maturity date is at the end of a non-31 day month it may be necessary to
 * enter a rollDate different from the start date.
 *
 * Amortizations (if provided) are assumed to coincide with payments.
 *
 * The function is copied from the ALIB GtoSwapRateWithRollDate function.
 */
static double swapRate(
    const ZeroCurve&          discountCurve, 
    const DateTime&           startDate,
    const DateTime&           rollDate,
    const DateTime&           maturityDate,
    const MaturityPeriod&     fixedIvl,
    const DayCountConvention& fixedDayCountConv,
    const bool                valueFloating,
    double                    floatPV,
    const ZeroCurve*          estimatingCurve,
    const double              floatSpread,
    const MaturityPeriod*     floatIvl,
    const DayCountConvention* floatDayCountConv,
    const bool                firstFixed,
    const double              firstFixRate,
    const bool                convexityDelayAdj,
    const IRVolBase*          volModelIR,
    const CashFlowArray*      amortSched,
    const StubPlacement&      stubPos,
    const BadDayConvention&   accBadDayConv,    // accrual bad day convention
    const BadDayConvention&   payBadDayConv,    // payment bad day convention
    const BadDayConvention&   resetBadDayConv,  // reset bad day convention
    const Holiday&            holidays);

/** Calculates a par swap rate for the swap starting at startDate,
 * maturing at maturityDate, with fixed day count convention 
 * fixedDayCountConv and fixed payments occurring at time 
 * intervals defined by interval. In other words, the routine calculates 
 * the fixed rate such that the present value of the fixed and 
 * floating sides are equal. The floating side is assumed to be at par.
 */
/*
    make it as template function so that it can be useful for YieldCurve and Zero curve, Daniel Ng
*/
template<class IRCURVE>
static double parSwapRate(
    const IRCURVE&          curve,
    const DateTime&           startDate,   // Date instrument begins at
    const DateTime&           maturityDate,// Date instrument matures at
    const MaturityPeriod&     period,      // Time between payments 
    const DayCountConvention& dcc,
    const Stub&               stubType,
    bool                      stubAtEnd,
    const BadDayConvention&   accBadDayConv,
    const BadDayConvention&   payBadDayConv,
    const Holiday&            holidays);

/** Makes a cash flow list for a swap instrument. */
static CashFlowArray cashflows(
    const DateTime&           startDate,
    const DateTime&           maturityDate,
    const Stub&               stubType,
    bool                      stubAtEnd,
    bool                      longStub,
    const BadDayConvention&   accrualBadDayConv,
    const BadDayConvention&   payBadDayConv,
    const Holiday&            holidays,
    bool                      keepStartDate,
    bool                      subtractInitial,
    bool                      addFinal,
    double                    rate,
    const MaturityPeriod&     period,           // as tenor, eg. 3M
    const DayCountConvention& dcc);

/** Makes a cash flow list for a swap instrument. */
static CashFlowArray cashflows(
    const DateTime&           startDate,
    const DateTime&           rollDate,
    const DateTime&           maturityDate,
    const Stub&               stubType,
    bool                      stubAtEnd,
    bool                      longStub,
    const BadDayConvention&   accrualBadDayConv,
    const BadDayConvention&   payBadDayConv,
    const Holiday&            holidays,
    bool                      keepStartDate,
    bool                      subtractInitial,
    bool                      addFinal,
    double                    rate,
    const MaturityPeriod&     period,           // as tenor, eg. 3M
    const DayCountConvention& dcc);

/** Makes a cash flow list for a swap instrument. */
static CashFlowArray cashflows(
    const DateTime&           startDate,
    const DateTime&           rollDate,
    const DateTime&           maturityDate,
    const Stub*               stub,
    bool                      stubAtEnd,
    bool                      longStub,
    const BadDayConvention*   accrualBadDayConv,
    const BadDayConvention*   payBadDayConv,
    const Holiday*            holidays,
    bool                      keepStartDate,
    bool                      subtractInitial,
    bool                      addFinal,
    double                    rate,
    int                       count,           // interval = count periods
    const string&             period,          // e.g. Y, M, W, D
    const DayCountConvention *dcc);

/** 
 * Makes a cash flow list.
 *
 * Based on the ALIB GtoMakeCFL2 function.
 */
static CashFlowArray cashflows(
    const double              rate,
    const DateTime&           startDate,
    const DateTime&           valueDate,
    const DateTime&           maturityDate,
    const Stub&               stubType,
    bool                      stubAtEnd,
    bool                      longStub,
    const BadDayConvention&   accrualBadDayConv,
    const BadDayConvention&   payBadDayConv,
    const Holiday&            holidays,
    bool                      keepStartDate,
    bool                      subtractInitial,
    bool                      addFinal,
    const MaturityPeriod&     period,           // as tenor, eg. 3M
    const DayCountConvention& dcc);

/** 
 * Makes a cashflow list for a swap floating leg.
 *
 * The function is copied from the ALIB GtoMakeCFLFloatingRoll function.
 */
static CashFlowArray cashflows(
    const double              spread,
    const ZeroCurve&          curve,              // estimating curve
    const DateTime&           startDate,
    const DateTime&           valueDate,
    const DateTime&           rollDate,
    const DateTime&           maturityDate,
    const MaturityPeriod&     interval,           // as tenor, eg. 3M
    const DayCountConvention& dayCountConv,
    const Stub&               stubType,
    const bool                stubAtEnd,
    const bool                longStub,
    const bool                subtractInitial,
    const bool                addFinal,
    const BadDayConvention&   accrualBadDayConv,
    const BadDayConvention&   payBadDayConv,
    const BadDayConvention&   resetBadDayConv,
    const Holiday&            holidays,
    const bool                firstFixed,
    const double              firstFixRate,
    const bool                isAdditive = true);   // is spread additive or percentage?

// create an array of dates separated by interval
static DateTimeArray* dateArray(
    const DateTime&       startDate,
    const DateTime&       rollDate,
    const DateTime&       maturityDate,
    const MaturityPeriod& interval,        // as tenor, eg. 3M
    bool                  stubAtEnd,
    bool                  constrainEndDates = true);

// create an array of dates separated by interval
static DateTimeArray* dateArray(
    const DateTime& startDate,
    const DateTime& rollDate,
    const DateTime& maturityDate,
    int             count,            // interval = count periods
    const string&   period,           // e.g. Y, M, W, D
    bool            stubAtEnd);

/** is a date on cycle ? */
static bool onCycle(
    const DateTime& valueDate,
    const DateTime& date,            // is this date on cycle ?
    int             count,           // interval = count periods
    const string&   period);         // e.g. Y, M, W, D

/** Checks if swap maturity dates are on the regular cycle
    from valueDate. Works only if floating and fixed freqs are equal.
*/
static bool swapDatesOnCycle(
    const DateTime&      valueDate,     
    const DateTimeArray& swapDates,      
    int                  freq);

/** Counts # intervals in a range of dates */
static void countDates(
    const DateTime& fromDate,
    const DateTime& toDate,
    int             count,            // interval = count periods
    const string&   period,           // e.g. Y, M, W, D
    int*            numIntervals,     // (O) Answer (Quotient) 
    int*            extraDays);       // (O) Days left over(remainder)

/** Returns a date list of standard euro-money swap dates.  This
    basically means dates of all coupons for a swap an integral number of
    years from now.  This differs from the coupons for the last swap, in that
    it may not be an integral number of years from the value date.
    The routine outputs a date list including all coupons for a swap 
    that's not after the given  maturity date and is after the second to
    last coupon date of the input swap.  And is an integral number of
    coupon frequencies from the value date.
*/
static DateTimeArray canonicalSwapDates(
    const DateTime& startDate,
    const DateTime& minDate,   // none before here
    const DateTime& endDate,
    int             frequency);

/** builds date list from start to end using given interval (stub @ end) */
static void simpleSwapDates(
    const DateTime& startDate,
    const DateTime& endDate,
    int             frequency,
    DateTimeArray&  dates); // (O) 

/** Makes a cash flow list for a swap instrument. */
static CashFlowArray cashflows(
    const DateTime& valueDate,
    const DateTime& maturity,
    bool            stubAtEnd,
    double          rate,
    int             count,           // interval = count periods
    const string&   period,          // e.g. Y, M, W, D
    const DayCountConvention *dcc);

/** Makes a cash flow list for a swap instrument. */
static CashFlowArray cashflows(
    const DateTime&           startDate,
    const DateTime&           maturity,
    const Stub*               stub,
    bool                      stubAtEnd,
    const BadDayConvention*   accrualBadDayConv,
    const BadDayConvention*   payBadDayConv,
    const Holiday*            holidays,
    bool                      keepStartDate,
    bool                      addPrincipal,
    double                    rate,
    int                       count,           // interval = count periods
    const string&             period,          // e.g. Y, M, W, D
    const DayCountConvention *dcc);

// build up accrue and pay dates for a swap
static void swapDates(
    const DateTime&           startDate,
    const DateTime&           maturity,
    bool                      stubAtEnd,
    const BadDayConvention*   accrualBadDayConv,
    const BadDayConvention*   payBadDayConv,
    const Holiday*            holidays,
    int                       count,           // interval = count periods
    const string&             period,          // e.g. Y, M, W, D
    DateTimeArray&            accrualDates,    // (O) 
    DateTimeArray&            payDates);       // (O)

/** create an array of swap payment dates */
static DateTimeArray* paymentDates(
    const DateTime& baseDate,        // start here
    const DateTime& maturity,        // end here
    int             count,           // interval = count periods
    const string&   period,          // e.g. Y, M, W, D
    bool            stubAtEnd);

/** create an array of dates separated by interval (given by count periods) */
static DateTimeArray* dateArray(
    const DateTime& baseDate,        // start here
    const DateTime& maturity,        // end here
    int             count,           // interval = count periods
    const string&   period,          // e.g. Y, M, W, D
    bool            stubAtEnd,
    bool            constrainEndDates = true);

/** create an array of dates separated by interval (given by count periods) */
static DateTimeArray* dateArray(
    const DateTime& startDate,       // start here
    const DateTime& endDate,         // end here
    int             count,           // interval = count periods
    const string&   period,          // e.g. Y, M, W, D
    int             startIdx,        // 0=start @ basedate, 1=start @ baseDate + interval
    bool            stubAtEnd,       // true=add dates forwards, false=add backwards
    int             numDates);       // how many dates

/** create an array of dates separated by interval (given by count periods) */
static DateTimeArray* dateArray(
    const DateTime& baseDate,        // start here
	int             count,           // interval = count periods
	const string&   period,          // e.g. Y, M, W, D
	int             startIdx,        // 0=start @ basedate, 1=start @ baseDate + interval
	int             arrayIncrement,  // Usually +1 or -1
	int             numDates);       // how many dates

private:
    SwapTool();

    static double calculateFloatLegWithFixing(
        const ZeroCurve&          discountCurve,
        const DateTime&           startDate,
        const DateTime&           rollDate,
        const DateTime&           maturityDate,
        const MaturityPeriod&     floatIvl,
        const DayCountConvention& floatDcc,
        const double              firstFixRate,
        const bool                stubAtEnd,
        const BadDayConvention&   badDayConv,
        const Holiday&            holidays);
};



template<class IRCURVE>
double SwapTool::parSwapRate(
                          const IRCURVE&          curve,
                          const DateTime&           startDate,   // Date instrument begins at
                          const DateTime&           maturityDate,// Date instrument matures at
                          const MaturityPeriod&     period,      // Time between payments 
                          const DayCountConvention& dcc,
                          const Stub&               stubType,
                          bool                      stubAtEnd,
                          const BadDayConvention&   accBadDayConv,
                          const BadDayConvention&   payBadDayConv,
                          const Holiday&            holidays)
{  // ALIB: szcbuild.c#3986 (ParSwapRate)
      static const string method = "SwapTool::parSwapRate";

        try 
        {
            DateTime startDateAdj = payBadDayConv.adjust(startDate, &holidays);
            DateTime matDateAdj = payBadDayConv.adjust(maturityDate, &holidays);

            CashFlowArray cfl = cashflows(
                startDate, 
                maturityDate, 
                StubSimple(), 
                stubAtEnd, 
                false,          // long stub
                accBadDayConv, 
                payBadDayConv, 
                holidays, 
                true, 
                false, 
                false, 
                1.0,    // coupon = 1 for annuity
                period, 
                dcc);

            if (cfl.empty())
            {
                string msg = Format::toString(
                    "no cashflows between %s and %s",
                    startDateAdj.toString().c_str(),
                    matDateAdj.toString().c_str());
                throw ModelException(method, msg);
              
            }

            // Get present value of 1 at startDate
            double startDatePV = curve.discountFactor(startDateAdj);


            // Compute coupon for the given zero-coupon rates
            double couponsPV = curve.pv(cfl);  // I think this is not correct.
                                                // The start and end dates in cfl must be adjusted . Daniel Ng

            
            if (!Maths::isPositive(couponsPV))
            {
                string msg = Format::toString("coupons with rate = 1.0 value to <= 0.0");
                throw ModelException(method, msg);
            }
            

            double lastPV = curve.discountFactor(matDateAdj);
            double couponRate = (startDatePV - lastPV) / couponsPV;
            return couponRate;
        }
        catch (exception&e )
        {
            throw ModelException(e, method);
        }

}


DRLIB_END_NAMESPACE
#endif
