//----------------------------------------------------------------------------
//
//   Group       : QR&D Core Analytics Team
//
//   Filename    : RecurseInfo.hpp
//
//   Description : Helper class for bootstrapping ALIB style zero curve
//
//   Author      : Richard Appleton
//
//   Date        : 2nd May 2005
//
//----------------------------------------------------------------------------

#ifndef ZC3_RECURSE_INFO_HPP
#define ZC3_RECURSE_INFO_HPP

#include "edginc/config.hpp"
#include "edginc/Atomic.hpp"
#include "edginc/DateTime.hpp"
#include "edginc/CashFlow.hpp"
#include "edginc/CashStream.hpp"
#include "edginc/FourPlusIZeroCurve.hpp"    // TBD!! replace with SimpleZeroCurve
#include "edginc/ZC3CurveInstrument.hpp"


DRLIB_BEGIN_NAMESPACE
using namespace std;

class BadDayConvention;
class Holiday;

class ZC3ZeroCurve;
class ZC3FuturesGapRule;
class ZC3ZeroInterpolation;

class MARKET_DLL ZC3RecurseInfo{
public:

    ZC3RecurseInfo(
        const ZC3SwapDataArray& swapsData,
        const BadDayConvention& badDayConv,
        const Holiday&          holidays);

    bool isFirstPass() const;

    /**
     * Evaluates curve to see whether recursive bootstrapping is complete. This
     * structure is modified to provide information for the next pass (if it is
     * necessary).
     */
    bool evaluateBootstrap(
        ZC3ZeroCurve&               szc,            // curve being built
        const DateTimeArray&        cashDateList,   // date list from non-swap instruments
        const ZC3SwapDataArray&     swapsData,      // swap instrument data
        const BadDayConvention&     badDayConv,     // bad day convention
        const Holiday&              holidays,       // holidays for date adjustment
        const Holiday&              basisHolidays,  // holidays for basis payments
        bool                        valueFloating,  // true if floating leg valued
        bool                        withBasis,      // true if has basis payments
        bool                        monthEndRoll,   // true if month end roll
        const MaturityPeriod*       swapFloatIvl,   // float interval
        const ZC3FuturesGapRule&    futuresGapRule, // cash/futures-swaps gap rule
        const ZC3ZeroInterpolation& zeroInterpType  // zero interpolation type
        );

    /**
     * This is used in recursive bootstrapping. Discount factors for payments
     * not corresponding to benchmark maturities are calculated from the
     * fitted curve produced in the previous iteration.
     */
    void addNonBenchmarkDiscounts(
        ZC3ZeroCurve&        discountCurve,
        const CashFlowArray& cfl) const;


    /**
     * Adds a single discount factor to a smooth zero curve. The discount
     * factor must already have been adjusted to remove the effects of
     * year-end adjustments. The interp type used is flat forwards.
     *
     * This function is used in the recursive bootstrapping process. In fact it
     * should only be used in such an iterative manner. Otherwise the absolute
     * adjustments may not be calculated correctly.
     * 
     * The date must be after the base date and the last date of the curve.
     * The call will fail if the previous segment was added with smoothing. 
     */
    void addAdjustedDiscount(
        ZC3ZeroCurve&   curve, 
        double          discount,
        const DateTime& discDate) const;

    void getCurrent(CashStream*& cs, CashFlowArray*& cfl);
    void next(const CashStream& cs, const CashFlowArray& cfl);
    void reset();

private:
    /**
     * This function is used with recursive bootstrapping methods. It is used
     * to select the set of benchmark swap end dates.
     */
    void findSwapBenchmarkDates(
        const DateTimeArray&     cashDateList,
        const DateTime&          lastCashDate,
        const ZC3SwapDataArray&  swapsData,
        bool                     valueFloating,
        bool                     withBasis,
        const MaturityPeriod&    floatIvl,
        const BadDayConvention&  badDayConv,
        const Holiday&           holidays,
        const Holiday&           basisHolidays,
        const ZC3FuturesGapRule& futuresGapRule);

    /**
     * This function is used with recursive bootstrapping methods. It is used
     * to select the set of cash flow dates which do not coincide with 
     * benchmark dates. It is at these days that the fitting process attempts
     * to minimize the difference between the discount factors resulting from
     * the bootstrapping and the discount factors predicted by interpolating
     * from the benchmarks using our curve model.
     */
    void findSwapCriticalDates(
        const ZC3ZeroCurve&     szc,
        const DateTime&         lastCashDate,
        const BadDayConvention& badDayConv,
        const Holiday&          holidays,
        bool                    valueFloating,
        bool                    withBasis,
        bool                    monthEndRoll,
        const MaturityPeriod&   floatIvl);

    /**
     * This function is used with recursive bootstrapping methods. It does
     * the following :-
     * (1) Creates array of discount factors at benchmark dates
     * (2) Creates array of discount factors at fitting dates by interpolating
     *       with a particular model (eg spline interpolation of DFs)
     * (3) Creates discount curve from above array in (2)
     * (4) Returns a single number showing net error between DF's at the fitting
     *       dates in the original curve and the new curve.
     */
    bool evaluateRecursiveBootstrap(
        ZC3ZeroCurve&               szc,
        const ZC3ZeroInterpolation& zeroInterpType);

    DateTimeArray        benchmarkDates;
    DateTimeArray        ignoreDates;
    DateTimeArray        fittingDates;
    FourPlusIZeroCurveSP fittedCurve;    // TBD!! use SimpleZeroCurve
    CashFlowArrayArray   cfls;
    CashStreamArray      cashStreams;
    DateTime             startDate;
    DateTime             lastSwapDate;
    DateTime             lastFittingDate;
    int                  passIndex;
    int                  swapIndex;
};


typedef refCountPtr<ZC3RecurseInfo> ZC3RecurseInfoSP;


DRLIB_END_NAMESPACE
#endif // ZC3_RECURSE_INFO_HPP
