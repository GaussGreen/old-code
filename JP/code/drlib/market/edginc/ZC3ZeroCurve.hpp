//----------------------------------------------------------------------------
//
//   Group       : QR Core Analytics Team
//
//   Filename    : ZC3ZeroCurve.hpp
//
//   Description : ALIB style zero curve
//
//   Author      : Richard Appleton
//
//   Date        : 22nd April 2005
//
//----------------------------------------------------------------------------

#ifndef ZC3_ZERO_CURVE_HPP
#define ZC3_ZERO_CURVE_HPP

#include "edginc/config.hpp"
#include "edginc/ZC3Interval.hpp"
#include "edginc/ZC3Interpolation.hpp"
#include "edginc/ZC3CurveSegment.hpp"
#include "edginc/ZC3FuturesStubRule.hpp"
#include "edginc/CashFlow.hpp"
#include "edginc/CashStream.hpp"
#include "edginc/ZeroCurve.hpp"
#include "edginc/YieldCurve.hpp"
#include <string>
#include <set>


using namespace std;    // string
DRLIB_BEGIN_NAMESPACE

class ZC3ZeroCurveHelper;

class DayCountConvention;
class IRVolBase;
class ZeroCurve3;
class StubPlacement;

class ZC3FuturesGapRule;
class RateData;
class ZC3RecurseInfo;


/**
 * Implementation of ALIB 'zero curve 3' style zero curve.
 *
 * The mutator functions are intended for use only during the curve 
 * bootstrapping process - after a curve has been built it should not be
 * subsequently amended.
 *
 * This class is not meant to be used directly by clients - it is intended to
 * be accessed via a yield curve which acts as a proxy onto the methods 
 * detailed below.
 */
class MARKET_DLL ZC3ZeroCurve : public ZeroCurve
{
public:
    static CClassConstSP const TYPE;

    /** Construct empty curve */
    ZC3ZeroCurve(const DateTime& baseDate, const string& stubInterpType);

    ZC3ZeroCurve(
        const ZC3ZeroCurve&           szcStub,
        const ZC3ZeroCurve*           discountCurve,
        const ZC3RateDataArray&       ratesData,
        const ZC3SwapDataArray&       swapsData,
        const ZC3TurnArray&           turnsData,
        const bool                    hasAdjustments,
        const bool                    mtmAdjustmentFlag,
        ZC3FuturesStubRuleSP&         futuresStubRule,
        const ZC3FuturesGapRule&      futuresGapRule,
        const DateTime*               futuresMMDate,
        const IRVolBase*              volModelIR,
        const StubPlacement&          stubPos,
        const BadDayConvention&       swapBadDayConvention,
        const bool                    valueFloating,
        const MaturityPeriod*         swapFloatIvl,
        const ZeroCurve*              fixingCurve,
        const ZC3CouponInterpolation* couponInterpolation,
        const ZC3ZeroInterpolation&   zeroInterpolation,
        const Holiday&                holidays,
        const Holiday&                basisHolidays,
        const DateTime*               extrapDate,
        const bool                    withBasis);

    virtual ~ZC3ZeroCurve();

    /**
     * A shifted curve should preserve forward rates between the same calendar 
     * dates.  This is achieved by not altering the original curve, but by
     * performing all discount functions using the original curve then rebasing
     * the results to the new base date.
     */
    void shift(const DateTime& newBaseDate);

    // ZeroCurve interface

    /**
     * Computes discount factor for a smooth zero curve at a given date.
     */
    double discountFactor(const DateTime& date) const;

    /**
     * Calculates an interpolated rate from a ZCurve at some date, with 
     * extrapolation for dates past end of curve.  Adjustments are taken
     * into account.
     */
    virtual double zeroCouponRate(const DateTime& date) const;

    /** how many fitted points are on curve ? */
    int length() const;
    
    const DateTime& getBaseDate() const;

    /** when is first date for which we have genuine information? */
    const DateTime& firstDate() const;

    /** when does it end ? [base date if no points on curve]*/
    const DateTime& endDate() const;
    
    /** strip out the dates */
    DateTimeArray getDates() const;

    /** strip out the rates and dates */
    CashFlowArraySP getRatesAndDates() const;


    // methods useful for bootstrapping

    void validate(bool strong) const;
    
    /** Get year fraction between two dates using curve day count convention */
    double years(const DateTime& startDate, const DateTime& endDate) const;

    /** Get year fraction between two dates assuming Act/365F day count convention */
    double years(int startDate, int endDate) const;

    /** Converts zc-style rate into a discount factor */
    double rateToDiscount(double rate, const DateTime& date) const;

    /** Converts zc-style rate into a discount factor */
    double rateToDiscount(double rate, double t) const;

    /** Converts continuous or zc-style rate into a discount factor */
    double rateToDiscount(double rate, const DateTime& startDate, const DateTime& endDate, bool continuous = false) const;

    /** Converts discount factor to zc-style rate */
    double discountToRate(double discount, const DateTime& date) const;

    /** Converts discount factor to zc-style rate */
    double discountToRate(double discount, double t) const;

    /** Converts discount factor into a zc-style or continuous rate */
    double discountToRate(double discount, const DateTime& startDate, const DateTime& endDate, bool continuous = false) const;

    /**
     * Calculates an interpolated rate from a ZCurve at some date, with 
     * extrapolation for dates past end of curve.  Adjustments are not
     * taken into account.
     */
    double unadjustedZeroCouponRate(const DateTime& date) const;

    /**
     * Computes adjustment factors for a smooth zero curve between two dates.
     *
     * If there are no adjustments, then the adjustment factor would be 1.
     *
     * If there are adjustments that cause us to expect higher than expected
     * forward rates (that is what adjustments are) then adjustment < 1.
     *
     * The adjustment is the same for both smoothed and unsmoothed discount
     * factors, unless there are absolute adjustments in the structure.
     */
    double adjustmentFactor(const DateTime& startDate, const DateTime& endDate) const;

    void extendLastRate(const DateTime& lastZeroDate, const DateTime& adjusted);

    /**
     * Set flag to indicate that unsmoothed discount factors are wanted.
     * Unsmoothed discount factors are required during curve bootstrapping.
     *
     * After getting an unsmoothed discount factor the curve state should 
     * be reset to use smoothing.
     */
    void setUseSmoothing(bool state) const;

    /**
     * Returns a key used to optimize repeated calculations of discount
     * factors/forward rate. The calc method for this key returns the 
     * natural logarithm of the discount factor.
     */
    virtual YieldCurve::IKey* logOfDiscFactorKey() const;

protected:
    ZC3ZeroCurve(CClassConstSP clazz);
    ZC3ZeroCurve(CClassConstSP clazz, const DateTime& valueDate);

private:
    ZC3ZeroCurve();
    ZC3ZeroCurve& operator=(const ZC3ZeroCurve& other); // undefined to prevent usage


    // methods used during bootstrapping

    /**
     * Adds a number of rates to create a new smooth zero curve.
     *
     * Four types of rate can be added. Each defines a rate between the start
     * date and the maturity date for that rate.
     *
     * These rate types are as follows:
     *    a) Standard forward or spot rate.
     *    b) Futures rate which will have MTM adjustments applied.
     *    c) Period-end relative adjustments.
     *    d) Period-end absolute adjustments.
     */
    void addRates(
        const ZC3ZeroCurve&         stubCurve, 
        const ZC3RateDataArray&     ratesData,
        const ZC3TurnArray&         turnsData,
        bool                        hasAdjustments,
        bool                        mtmAdjustmentFlag, // futures flag #1
        ZC3FuturesStubRuleSP&       futuresStubRule,   // futures flag #2
        const ZC3FuturesGapRule&    futuresGapRule,    // futures flags #3 & 4
        const DateTime*             futuresMMDate,
        const IRVolBase*            volModelIR,
        const ZC3ZeroInterpolation& interpolation, 
        const DateTime*             extrapDate,
        const ZC3ZeroCurve*         discountCurve,
        const Holiday&              holidays);

    /**
     * Main worker function for addRates.
     * 
     * This function does not perform any mark to market adjustments.
     *
     * The plan is to build the zero curve using the unadjusted futures rates 
     * using this function, then adjust the futures rates and call this 
     * function again. There is a potential for iteration!
     */
    void addRates(
        const ZC3ZeroCurve&         stubCurve, 
        const ZC3RateDataArray&     ratesData,
        const ZC3TurnArray&         turnsData,
        const ZC3ZeroInterpolation& interpolation, 
        const DateTime*             extrapDate,
        const ZC3ZeroCurve*         discountCurve,
        const Holiday&              holidays,
        bool                        withAdjustments);

    bool futuresMarkToMarket(
        const ZC3RateDataArray& ratesData, 
        const IRVolBase&        volModelIR) const;

    /**
     * Converts input rates to continuous rates and adjustments.
     */
    void getRatesAndAdjustments(
        const ZC3RateDataArray& ratesData, 
        const ZC3TurnArray&     turnsData,
        ZC3IntervalArray&       intervals, 
        ZC3StubArray&           stubs,
		bool                    withBasis,
        bool                    withAdjustment);

    void bootstrap(
        CashStream&                 cashStream,
        double                      presentValue,
        const ZC3ZeroInterpolation& interpolation,
        const ZC3ZeroCurve*         discountCurve,
        const IRVolBase*            volModelIR);

    /**
     * Add multiple segments using the shape of the discount curve. Effectively
     * we add multiple segments which are the rates in the discount curve plus 
     * a spread.
     */
    void addShapedSegments(
        const ZC3ZeroInterpolation& interpolation,
        const ZC3ZeroCurve&         discountCurve, 
        CashStream&                 cashStream, 
        const double                pvRemainder, 
        const IRVolBase*            volModelIR);

    
    /**
     * populates curve with rates derived from shape array and spread.
     */
    void populateShapedCurve(
        const int          oldSize, 
        const int          numSegs, 
        const DoubleArray& shapeRates, 
        const double       spread, 
        const BoolArray&   useAdjAbsolutes, 
        const bool         calcZero);


    /**
     * Performs the smoothing algorithm for the last and next-to-last segment
     * of a smooth zero curve.
     *
     * Uses just in time methods to calculate and store the duration of
     * segments of the zero curve.
     */
    double smoothing(
        const CashStream&       cashStream,     // cash stream to add (all inclusive)
        double                  presentValue,   // of cash stream at curve base date
        const ZeroCurve*        discountCurve,
        const CashStream&       cs,             // remaining cash flows
        const IRVolBase*        volModelIR,
        bool                    foundAdjAbs);

    /**
     * Under certain circumstances it is necessary to revert a zero curve to
     * the state before the previous segment was added. This is needed for
     * example when GtoZeroCurve3CB is being to build an estimating curve with
     * multiple interpolation types.
     */
    ZC3ZeroInterpolationConstSP revert();

    /**
     * Adds a number of swaps to a smooth zero curve. The zero curve is amended
     * in place.
     */
    void swaps2(
        const ZC3ZeroCurve*           discountCurve,
        const ZC3SwapDataArray&       swapsData,
        const BadDayConvention&       bdc,
        const StubPlacement&          stubPos,
        const Holiday&                holidays,
        const Holiday&                basisHolidays,
        const ZC3CouponInterpolation* couponInterp,
        const ZC3ZeroInterpolation&   zeroInterp,
        const ZC3FuturesGapRule&      gapRule,
        const DateTime*               extrapDate,
        const bool                    convDelayAdj,
        const IRVolBase*              volModelIR,
        ZC3RecurseInfo*               recurseInfo,
        const bool                    withBasis);

    /**
     * Adds a single swap to a smooth zero curve.
     *
     * The swap could be a genuine or a synthetic swap.
     */
    void addOneSwap(
        const ZC3ZeroCurve*         discountCurve,
        ZC3SwapData&                swapData,
        double                      couponRate,
        const BadDayConvention&     badDayConv,
        const StubPlacement&        stubPos,
        const Holiday&              holidays,
        const Holiday&              basisHolidays,
        const ZC3ZeroInterpolation& zeroInterp,
        const bool                  convDelayAdj,
        const IRVolBase*            volModelIR,
        bool                        recursive,
        ZC3RecurseInfo*             recurseInfo,
        const bool                  withBasis);

    /**
     * "Safe" version of bootstrap which checks that the cash stream does
     * extend the zero curve before calling it.  This will be based on the 
     * close enough principle of 5 days to protect against really unsorted
     * swaps.
     *
     * Note that since there has already been some culling of swaps before
     * we get here (based on unadjusted dates), the only chance of having
     * a swap cash stream which does not extend the zero curve is when we
     * have a holiday which forces the swap cash stream to end earlier than
     * expected (i.e. when the cash stream matures naturally at the end of the
     * month, but that it is a holiday and we are using modified-following).
     *
     * Despite the rarity of the errors, this solution is not 100% ideal.
     *
     * What we really need to do is this:
     *
     * 1. Sort the swaps.
     * 2. Change bootstrap to be configurable regarding whether an
     *    instrument extends the curve or not.
     *
     * Since the ideal solution has too many side effects, we are going to 
     * concentrate for the moment on a solution which is good enough.
     */
    void safeBootstrap(
        CashStream&                 cashStream,
        double                      presentValue,
        const ZC3ZeroInterpolation& interpolation,
        const ZC3ZeroCurve*         discountCurve,
        const IRVolBase*            volModelIR);


    /**
     * Add uncalculated segment to curve.
     */
    ZC3CurveSegment& addSegment(
        const DateTime&             date, 
        const ZC3ZeroInterpolation& interpolation, 
        ZC3CurveSegment*&           prevSegment);   // I/O

    /**
     * Adds a single rate to a smooth zero curve. The zero curve is amended in
     * place via a call to bootstrap.
     *
     * The rate cannot start before the base date of the curve, and cannot
     * mature before the current last date of the curve.
     */
    void addRate(
        double                      rate, 
        const DateTime&             startDate, 
        const DateTime&             endDate, 
        const DayCountConvention&   dayCountConv,
        int                         basis,
        const ZC3ZeroInterpolation& interpolation);

    /** 
     * Adds a number of stub rates before the base date of the zero curve in
     * such a way as to guarantee that the same results will be given for dates
     * beyond the base date.
     *
     * The zero curve is amended in place.
     */
    void addStub(
        const ZC3StubArray&       stubs, 
        const DayCountConvention& rateDayCountConv, 
        int                       rateBasis);

    /**
     * Check that the stubs are contiguous and join up to the base date of the
     * zero curve.
     */
    void checkStubs(const ZC3StubArray& stubs) const;

    /**
     * Validate the inputs for addStub. 
     *
     * We require that the start dates are in strictly ascending order and that
     * all dates are strictly before the base date of the zero curve.
     *
     * We also expect the base date of the zero curve to be EQUAL to the first
     * date in the zero curve.
     */
    // TBD!! overlaps with checkStubs ???
    void addStubValidateInputs(const ZC3StubArray& stubs, const string& method) const;

    /**
     * Add an adjustment in date order to list of adjustments.
     */
    void addAdjustmentToList(const ZC3Adjustment& adjustment);

    /**
     * Finds out whether there are any absolute adjustments in a particular
     * segment of a zero curve.
     */
    bool findAbsoluteAdjustments(
        const ZC3CurveSegment& prevSegment, 
        const ZC3CurveSegment& segment) const;

    /**
     * Modifies an adjustment in a smooth zero curve - must exist already.
     */
    void modifyAdjustment(const ZC3Adjustment& adjustment, double contRate);

    /**
     * Computes adjustment factors between two dates.
     *
     * If there are no adjustments, then the adjustment factor would be 1.
     *
     * If there are adjustments that cause us to expect higher than expected
     * forward rates (that is what adjustments are) then adjustment < 1.
     */
    double adjustmentFromList(
        const ZC3AdjustmentArray& adjustments, 
        const int                 startDate, 
        const int                 endDate) const;

	void flatFwdAdjAbsolute(
		const ZC3CurveSegment& prevSegment,
		const ZC3CurveSegment& segment);

    /**
     * Sets the dates and rates within a curve using the critical dates which 
     * are sufficient to reprice any benchmarks used when building the curve 
     * originally.
     */
    void encapsulatedCurveSetRates();

    DateTimeArraySP getSegmentDates() const;

    void addDiscountFactor(
        const DateTime&             date, 
        double                      disc, 
        const ZC3ZeroInterpolation& interpolation);

    void setExtrapolationDate(const DateTime* date);
    void addCriticalDate(const DateTime& date);

    /**
     * Validate benchmark flag.
     *
     * Returns true if instrument is to be used in bootstrapping.
     */
    bool checkBenchmarkFlag(int benchmark, const string& instrument) const;

    // delete all fitted data so far [used in recursive algorithms]
    void reset(const ZC3ZeroCurve* copy);

    // delete all temporary data
    void cleanup();
    
    // helper classes used during bootstrapping
    friend class ZeroCurve3;
    friend class ZC3FlatForwardsZeroInterpolation;
    friend class ZC3SmoothForwardsZeroInterpolation;
    friend class ZC3LinearZeroInterpolation;
    friend class FlatFwdObjFunc;
    friend class LinearObjFunc;
    friend class SmoothFunc;
    friend class ShapedFlatFwdObjFunc;
    friend class ZC3Turn;
    friend class ZC3AbsoluteTurn;
    friend class ZC3RecurseInfo;
    friend class ZC3CurveSegment;
    friend class ZC3FuturesStubRuleFlatFwds;
    friend class ZC3ZeroCurveHelper;
    
    // helper class for logOfDiscKey
    class ZC3LogOfDiscFactorKey;
    friend class ZC3LogOfDiscFactorKey;

    // for local market curves
    friend class LMParabCurve;
    friend class ZCBrzFI;

    // fields
    DateTime                  valueDate;     // value date if curve has been shifted
    double                    valueDateDiscount;
    DateTime                  baseDate;      // base date
    DateTime                  extrapDate;    // extrapolation date
    int                       basis;         // basis of rates

    DayCountConventionConstSP dcc;           // time measure for zero curve
    ZC3CurveSegmentArray      data;          // zero curve points
    DateTimeArray             criticalDates; // list of critical dates when built
    set<DateTime>             tmpCriticalDates; // set used during boostrapping for speed $unregistered
    CashFlowArray             datesAndRates; // transient field for get rates and dates (not analytics!)
    ZC3AdjustmentArray        adjustments;   // array of adjustments
    ZC3AdjustmentArray        adjAbsolute;   // absolute adjustment information
    CashStreamSP              cs;            // cash flow list for last segment $unregistered
    double                    PV;            // PV of last segment at its start date $unregistered
    bool                      noBootstrap;   // can we continue to bootstrap?

    // transient field used when building smoothed curve
    mutable bool              useSmoothing; // $unregistered

    // transient fields used for internal calculations when multiple 
    // interpolation types are needed.
    CashStreamSP              prevFull;      // previous full cashstream $unregistered
    CashStreamSP              prevCs;        // previous value of cs $unregistered
    double                    prevPv;        // previous last segment PV $unregistered
    ZC3CurveSegmentSP         prevPoint;     // previous smooth ZC point $unregistered

    // used for improving speed of interpolation
    mutable int               loBound; // $unregistered
    mutable int               hiBound; // $unregistered
};


typedef smartConstPtr<ZC3ZeroCurve>         ZC3ZeroCurveConstSP;
typedef smartPtr<ZC3ZeroCurve>              ZC3ZeroCurveSP;
typedef array<ZC3ZeroCurveSP, ZC3ZeroCurve> ZC3ZeroCurveArray;



DRLIB_END_NAMESPACE
#endif // ZC3_ZERO_CURVE_HPP
