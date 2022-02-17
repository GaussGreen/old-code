//----------------------------------------------------------------------------
//
//   Group       : QR&D Core Analytics Team
//
//   Filename    : ZC3CurveInstrument.hpp
//
//   Description : Benchmark for curve bootstrapping
//
//   Author      : Richard Appleton
//
//   Date        : 26th April 2005
//
//----------------------------------------------------------------------------

#ifndef ZC3_CURVE_INSTRUMENT_HPP
#define ZC3_CURVE_INSTRUMENT_HPP

#include "edginc/config.hpp"
#include "edginc/CashStream.hpp"
#include "edginc/DateTime.hpp"
#include "edginc/CashFlow.hpp"
#include "edginc/MaturityPeriod.hpp"
#include "edginc/CurrencyBasis.hpp"
#include "edginc/Stub.hpp"
#include <vector>



DRLIB_BEGIN_NAMESPACE
class BadDayConvention;
class DayCountConvention;
class IRVolBase;
class StubPlacement;
class ZeroCurve;
class ZC3ZeroCurve;
class ZC3Stub;


/**
 * Base class for instrument helper classes used for zero curve 3 bootstrapping
 * process.
 *
 * These classes should not be used outside of zero curve 3 bootstrapping.
 */
class MARKET_DLL ZC3CurveInstrument
{
public:

    virtual ~ZC3CurveInstrument();

    virtual double getContinuousRate(bool withBasis, bool withAdjustment = true) const = 0;

    const DateTime& getStartDate() const;
    const DateTime& getEndDate() const;

    // for sorting into maturity date ascending order using STL functions
    bool operator<(const ZC3CurveInstrument& other) const;

    bool   isStub(const DateTime& baseDate) const;
    double getRate() const;
    double getBasisRate() const;

protected:
    ZC3CurveInstrument(
        const DateTime& startDate, 
        const DateTime& endDate, 
        double          rate,
        double          basisRate);

    DateTime            startDate;
    DateTime            endDate;
    double              rate;
    double              basisRate;
};


class MARKET_DLL ZC3RateData : public ZC3CurveInstrument
{
public:
    static ZC3RateData* interpolate(
        const ZC3RateData* lastBefore, 
        const ZC3RateData& firstAfter, 
        const DateTime&    startDate,
        const DateTime&    endDate);

    ZC3RateData(
        const DateTime&           startDate, 
        const DateTime&           endDate, 
        double                    rate,
        const DayCountConvention& dcc,
        double                    basisRate);

    virtual bool futuresMarkToMarket(const ZC3ZeroCurve& szc, const IRVolBase& volModelIR) = 0;

    /**
     * See if there are any futures or MM rates that we wish to exclude.
     */
    virtual void ratesInitUsage(
        const DateTime& baseDate,
        const DateTime* futuresMMDate,
        DateTime&       futuresStartDate, 
        DateTime&       futuresEndDate, 
        DateTime&       badMMDate, 
        DateTime&       lastMMDate) = 0;

    /**
     * Prepare a stub. When we loop through the list of swaps some of the 
     * instruments might effectively be money market stubs before the start of
     * the curve.
     *
     * We will always convert these into ACT/360 simple rates, taking account
     * of currency basis.
     *
     * We need to do this because we add all the stubs in one function call,
     * and need to check for contiguity at that time.
     */
    void prepareStub(vector<ZC3Stub>& stubs, bool withBasis) const;

    virtual bool isReallyBadDate(DateTime& badMMDate) const;
    virtual bool isCashStartingOn(const DateTime& date) const;

    const DayCountConvention& getDcc() const;

protected:
    const DayCountConvention* dcc;
};


class MARKET_DLL ZC3CashData : public ZC3RateData
{
public:
    ZC3CashData(
        const DateTime&             startDate, 
        const DateTime&             endDate, 
        double                      rate, 
        const DayCountConvention&   dcc,
        double                      basisRate);

    double getContinuousRate(bool withBasis, bool withAdjustment = true) const;

    bool   futuresMarkToMarket(const ZC3ZeroCurve& szc, const IRVolBase& volModelIR);

    void ratesInitUsage(
        const DateTime& baseDate,
        const DateTime* futuresMMDate,
        DateTime&       futuresStartDate, 
        DateTime&       futuresEndDate, 
        DateTime&       badMMDate, 
        DateTime&       lastMMDate);

    bool isCashStartingOn(const DateTime& date) const;
};


class MARKET_DLL ZC3FuturesData : public ZC3RateData
{
public:
    ZC3FuturesData(
        const DateTime&             startDate, 
        const DateTime&             endDate, 
        double                      rate, 
        const DayCountConvention&   dcc,
        double                      adjustment,
        double                      basisRate);

    double getContinuousRate(bool withBasis, bool withAdjustment = true) const;

    bool   futuresMarkToMarket(const ZC3ZeroCurve& szc, const IRVolBase& volModelIR);

    void ratesInitUsage(
        const DateTime& baseDate,
        const DateTime* futuresMMDate,
        DateTime&       futuresStartDate, 
        DateTime&       futuresEndDate, 
        DateTime&       badMMDate, 
        DateTime&       lastMMDate);

private:
    bool isReallyBadDate(DateTime& badMMDate) const;

    int    basis;     // eg. CompoundBasis::SIMPLE
    string type; // futures type "1M", "3M", "AUD", "USD_TBILL"
    double adjustment;
};


class MARKET_DLL ZC3Turn : public ZC3RateData
{
public:
    ZC3Turn(
        const DateTime&           startDate,
        const DateTime&           endDate,
        double                    rate,
        const DayCountConvention& dcc);

    virtual void addAdjustment(ZC3ZeroCurve& curve) = 0;

    virtual double turnEffectInPeriod(double* turnYF, int* turnDays) const = 0;

    bool futuresMarkToMarket(const ZC3ZeroCurve& szc, const IRVolBase& volModelIR);

    void ratesInitUsage(
        const DateTime& baseDate,
        const DateTime* futuresMMDate,
        DateTime&       futuresStartDate, 
        DateTime&       futuresEndDate, 
        DateTime&       badMMDate, 
        DateTime&       lastMMDate);
};


class MARKET_DLL ZC3AbsoluteTurn : public ZC3Turn
{
public:
    ZC3AbsoluteTurn(
        const DateTime& startDate,
        const DateTime& endDate,
        double          rate,
        const DayCountConvention& dcc);

    void   addAdjustment(ZC3ZeroCurve& curve);
    double getContinuousRate(bool withBasis, bool withAdjustment = true) const;
    double turnEffectInPeriod(double* turnYF, int* turnDays) const;
};


class MARKET_DLL ZC3RelativeTurn : public ZC3Turn
{
public:
    ZC3RelativeTurn(
        const DateTime& startDate,
        const DateTime& endDate,
        double          rate,
        const DayCountConvention& dcc);

    void   addAdjustment(ZC3ZeroCurve& curve);
    double getContinuousRate(bool withBasis, bool withAdjustment = true) const;
    double turnEffectInPeriod(double* turnYF, int* turnDays) const;
};


class MARKET_DLL ZC3ImpliedTurn : public ZC3Turn
{
public:
    ZC3ImpliedTurn(
        const DateTime& startDate,
        const DateTime& endDate,
        double          rate,
        const DayCountConvention& dcc);

    void   addAdjustment(ZC3ZeroCurve& curve);
    double getContinuousRate(bool withBasis, bool withAdjustment = true) const;
    double turnEffectInPeriod(double* turnYF, int* turnDays) const;
};


class ZC3SwapData;
typedef refCountPtr<ZC3SwapData>      ZC3SwapDataSP;
typedef vector<ZC3SwapDataSP>         ZC3SwapDataArray;
typedef refCountPtr<ZC3SwapDataArray> ZC3SwapDataArraySP;



class MARKET_DLL ZC3SwapData : public ZC3CurveInstrument
{
public:
    static ZC3SwapDataArraySP getInterpRates(
            const ZC3SwapDataArray& data,
            const ZeroCurve&        discountCurve,
            const ZeroCurve*        estimatingCurve,
            const BadDayConvention& badDayConv,
            const StubPlacement&    stubPos,
            const Holiday&          holidays,
            const Holiday&          basisHolidays,
            const bool              convDelayAdj,
            const IRVolBase*        volModelIR,
            const bool              annualize,
            const bool              withBasis);

    /** 
     * Returns TRUE if a stub should be at the end. The decision is
     * based on the following.
     * 1. First the default stub position is checked. If this is a back
     *    stub then return TRUE - if a front stub return FALSE.
     * 2. If the default stub position is auto then the function will
     *    return TRUE, unless there is a stub in which case FALSE is returned.
     */
    ZC3SwapData(
        const DateTime&             baseDate,
        const DateTime&             startDate, 
        const DateTime&             endDate, 
        double                      rate, 
        const double*               adjustment,
        double                      price,
        const MaturityPeriod&       fixedIvl,
        const DayCountConvention&   fixedDcc,
        const MaturityPeriod*       floatIvl,
        const DayCountConvention*   floatDcc,
        int                         includeFlag,
        bool                        valueFloating,
        bool                        floatRateFixed,
        const Holiday&              holidays,
        const BadDayConvention&     badDayConv,
        const ZeroCurve*            fixCurve,
        const MaturityPeriod*       basisIvl,
        const DayCountConvention*   basisDcc,
        double                      basisRate);

    bool isForPass(int pass) const;
    bool isExcluded(const DateTime& lastZeroDate, const DateTime& baseDate) const;

    double getContinuousRate(bool withBasis, bool withAdjustment = true) const;

    void prepareStub(
        vector<ZC3Stub>&          stubs,
        double                    couponRate, 
        const BadDayConvention&   badDayConv, 
        const Holiday&            holidays,
        const bool                withBasis) const;

    CashFlowArraySP getKnownCashflows(
        const ZC3ZeroCurve&    szc,
        double                  couponRate,
        const StubPlacement&    stubPos,
        const BadDayConvention& badDayConv, 
        const Holiday&          holidays,
        const Holiday&          basisHolidays,
        const bool              withBasis,
        const bool              estimating) const;

    void getFloatingCashflows(
        CashStream&             cs,
        double                  couponRate,
        const StubPlacement&    stubPos,
        const BadDayConvention& badDayConv, 
        const Holiday&          holidays,
        const bool              mustValueFloating) const;

    /**
     * Get swap rate from zero curve adding any adjustment..
     */
    double getCouponRate(
        const ZeroCurve&        discountCurve, 
        const ZeroCurve*        estimatingCurve,
        const bool              convDelayAdj,
        const IRVolBase*        volModelIR,
        const StubPlacement&    stubPos, 
        const BadDayConvention& badDayConv,
        const Holiday&          holidays,
        const bool              withBasis) const;

    // used for bootstrapping curves with multiple interpolation types
    bool isRecursiveIncludedInFlatSection() const;
    void setRecursiveIncludedInFlatSection();

    void setBenchmark(int benchmark);

private:
    static DateTime removeCouponsNotNeeded(
        const DateTimeArray&    dl, 
        const DateTime&         lastCurveDate, 
        const BadDayConvention& badDayConv, 
        const Holiday&          holidays);

    static StubSP stubType;

    /*
     * Validates that coupon interpolation makes sense for a set of coupon 
     * dates.
     *
     * As a by-product returns the following:
     * 1. last curve date
     * 2. common fixed swap frequency
     * 3. common basis swap frequency (if relevant)
     */
    static void validateSwapDatesToInterp(
        const ZC3SwapDataArray& data,
        const bool withBasis, 
        const bool indexSwap);

    /** 
     * Returns TRUE if a stub should be at the end. The decision is
     * based on the following.
     * 1. First the default stub position is checked. If this is a back
     *    stub then return TRUE - if a front stub return FALSE.
     * 2. If the default stub position is auto then the function will
     *    return TRUE, unless there is a stub in which case FALSE is returned.
     */
    bool isEndStub(const StubPlacement& stubPos, const MaturityPeriod& interval) const;

    CashFlowArraySP makeCFL(
        double                    couponRate, 
        bool                      endStub, 
        const BadDayConvention&   badDayConv,
        const Holiday&            holidays,
        const MaturityPeriod&     interval,
        const DayCountConvention& dcc,
        bool                      includeNotionalFlows) const;

    DateTimeArraySP getFixedDates(bool stubAtEnd) const;

    int                       benchmark;
    double                    adjustment;
    double                    price;
    const MaturityPeriod&     fixedIvl;
    const DayCountConvention& fixedDcc;
    bool                      valueFloating;
    const MaturityPeriod&     floatIvl;
    const DayCountConvention& floatDcc;
    bool                      floatFixed;
    double                    floatFixRate;
    const MaturityPeriod&     basisIvl;
    const DayCountConvention& basisDcc;
    bool                      recursiveInFlatSection;
};


typedef refCountPtr<ZC3RateData>      ZC3RateDataSP;
typedef vector<ZC3RateDataSP>         ZC3RateDataArray;
typedef refCountPtr<ZC3RateDataArray> ZC3RateDataArraySP;

typedef refCountPtr<ZC3Turn>          ZC3TurnSP;
typedef vector<ZC3TurnSP>             ZC3TurnArray;
typedef refCountPtr<ZC3TurnArray>     ZC3TurnArraySP;

DRLIB_END_NAMESPACE

#endif // ZC3_CURVE_INSTRUMENT_HPP
