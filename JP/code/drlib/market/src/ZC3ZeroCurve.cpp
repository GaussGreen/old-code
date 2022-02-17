//----------------------------------------------------------------------------
//
//   Group       : QR Core Analytics Team
//
//   Filename    : ZC3ZeroCurve.cpp
//
//   Description : ALIB style zero curve
//
//   Author      : Richard Appleton
//
//   Date        : 22nd April 2005
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/ZC3ZeroCurve.hpp"
#include "edginc/CashStream.hpp"
#include "edginc/StubPlacement.hpp"
#include "edginc/ZC3Interval.hpp"
#include "edginc/ZC3Stub.hpp"
#include "edginc/ZC3CurveInstrument.hpp"
#include "edginc/ZC3CurveSegment.hpp"
#include "edginc/ZC3FuturesGapRule.hpp"
#include "edginc/ZC3FuturesStubRule.hpp"
#include "edginc/ZC3Interpolation.hpp"
#include "edginc/ZC3Iteration.hpp"
#include "edginc/ZC3RecurseInfo.hpp"
#include "edginc/Atomic.hpp"
#include "edginc/Maths.hpp"
#include "edginc/ModelException.hpp"
#include "edginc/CompoundBasis.hpp"
#include "edginc/RateConversion.hpp"
#include "edginc/Actual365F.hpp"
#include <math.h>
#include <algorithm>


DRLIB_BEGIN_NAMESPACE

using namespace std;

static bool tmp1 = lessThan(ZC3RateDataSP(   ), ZC3RateDataSP(   ));
static bool tmp2 = lessThan(ZC3SwapDataSP(   ), ZC3SwapDataSP(   ));
static bool tmp3 = lessThan(ZC3TurnSP(   ), ZC3TurnSP(   ));


// ALIB: zcbuild.c#997 (GtoZeroCurve3)
ZC3ZeroCurve::ZC3ZeroCurve(
    const ZC3ZeroCurve&           szcStub,
    const ZC3ZeroCurve*           discountCurve,   // non-NULL for estimating curve
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
    const bool                    withBasis)
  : ZeroCurve(TYPE), 
    /*
     * Convert the stub curve into a smooth zero curve
     */
    valueDate(szcStub.baseDate),
    valueDateDiscount(1.0),
    baseDate(szcStub.baseDate), 
    basis(szcStub.basis), 
    dcc(szcStub.dcc), 
    data(szcStub.data),
    criticalDates(0),
    tmpCriticalDates(szcStub.criticalDates.begin(), szcStub.criticalDates.end()),
    datesAndRates(szcStub.datesAndRates),
    adjustments(szcStub.adjustments),
    adjAbsolute(szcStub.adjAbsolute),
    cs(szcStub.cs),
    PV(szcStub.PV),
    noBootstrap(false),
    useSmoothing(true),
    prevFull(szcStub.prevFull),
    prevCs(szcStub.prevCs),
    prevPv(szcStub.PV),
    prevPoint(NULL),
    loBound(0), 
    hiBound(0)
{
    static const string method = "ZC3ZeroCurve::ZC3ZeroCurve";

    bool hasDiscCurve = (discountCurve != NULL);

    ZC3RecurseInfoSP recurseInfo;

    /* 
     * To handle the single step spline the bootstrapping process is recursive.
     * We recurse until the errors in the curve are below a defined level. For
     * all other bootstrapping methods we set the error to zero after the first
     * pass to make sure that the code is only executed once.
     */
    ZC3ZeroInterpolationConstSP firstInterpType(&zeroInterpolation);
    bool recursive = firstInterpType->isRecursive();
    if (recursive)
    {
        firstInterpType = ZC3ZeroInterpolation::make("Flat");

        if (couponInterpolation != NULL)
        {
            string msg = "Cannot use coupon interpolation with single step spline interpolation";
            throw ModelException(method, msg);
        }

        recurseInfo = ZC3RecurseInfoSP(
			new ZC3RecurseInfo(swapsData, swapBadDayConvention, holidays));
    }

    /*
     * Add the single rate instruments
     */
    if (!ratesData.empty() || !turnsData.empty())
    {
        addRates(szcStub, ratesData, turnsData, hasAdjustments,
            mtmAdjustmentFlag, futuresStubRule, futuresGapRule,
            futuresMMDate, volModelIR, *firstInterpType, extrapDate, 
            discountCurve, holidays);
    }

    ZC3ZeroCurveSP szcCash;
    DateTimeArraySP cashDateList;

    DateTime lastCashDate = endDate(); // OK as > 0 points on curve
    if (recursive && !hasDiscCurve)
    {
        cashDateList = getSegmentDates();
        if (cashDateList->empty())
        {
            throw ModelException(method,"There are no cash dates present.");
        }

        // Keep a copy of cash part for later iterations
        szcCash.reset(new ZC3ZeroCurve(*this));
    }

    /*
     * Add the swaps
     */
    bool monthEndRoll = false;
    for (bool morePasses = true ; morePasses ; ) //zcbuild.c#1258
    {
        if (recursive && !recurseInfo->isFirstPass() && !hasDiscCurve)
        {
            reset(szcCash.get());
        }

        if (!swapsData.empty())
        {
            // monthEndRoll = ((*swapsData)[0].prd_type == 'F');  // TBD!! flexible end of month

            // Chase 10 method [zcbuild.c#1272] not supported
            swaps2(
                discountCurve,
                swapsData,
                swapBadDayConvention,
                stubPos,
                holidays,
                basisHolidays,
                couponInterpolation,
                zeroInterpolation,
                futuresGapRule,
                extrapDate,
                false,
                volModelIR,
                recurseInfo.get(),
                withBasis);
        }

        if (recursive && !hasDiscCurve)
        {
            // zcbuild3.c#1338 & 2234

            // GCC will not pass ref to non-const temporary
            ZC3ZeroCurve tmp(discountCurve ? *discountCurve : ZC3ZeroCurve());
            
            bool complete = recurseInfo->evaluateBootstrap(
                discountCurve ? tmp : *this,
                *cashDateList,
                swapsData,
                swapBadDayConvention,
                holidays,
                holidays,
                withBasis ? false : valueFloating,
                withBasis,
                monthEndRoll,
                swapFloatIvl,
                futuresGapRule,
                zeroInterpolation);

            if (complete)
            {
                morePasses = false; // finish off curve here
            }
            else
            {
                reset(NULL);
            }
        }
        else
        {
            morePasses = false; // make sure we only do a single pass
        }
    }

    // chase spline method at zcbuild.c#1385 not supported

    encapsulatedCurveSetRates();

    /*
     * NB. cleanup of temporary data has to be initiated by ZeroCurve3, as copy
     * constructor may be used during bootstrapping (eg. for multiple segments).
     * The copy constructor copies values from the previous curve instance, so
     * they must have been preserved after the end of this method.
     */
}


// create new empty curve
ZC3ZeroCurve::ZC3ZeroCurve(const DateTime& pBaseDate, const string& pStubInterpType)
  : ZeroCurve(TYPE),
    valueDate(pBaseDate),
    valueDateDiscount(1.0),
    baseDate(pBaseDate),
    basis(CompoundBasis::ANNUAL),
    dcc(new Actual365F()),  // should restrict to Act/365F or Act/360
    criticalDates(0),
    datesAndRates(0),
    cs(new CashStream()),
    noBootstrap(false),
    useSmoothing(true),
    prevFull(NULL),
    prevCs(new CashStream()),
    prevPv(0.0),
    prevPoint(NULL),
    loBound(0), 
    hiBound(0)

{
    // ALIB: zcsmooth.c#1394
    // always have a point at base date
    addDiscountFactor(baseDate, 1.0, *ZC3ZeroInterpolation::make(pStubInterpType));
}


ZC3ZeroCurve::~ZC3ZeroCurve()
{
    adjustments.resize(0);
    adjAbsolute.resize(0);
}


void ZC3ZeroCurve::shift(const DateTime& pNewBaseDate)
{
    if (pNewBaseDate.empty())
    {
        throw ModelException("ZC3ZeroCurve::shift", "Invalid new base date");
    }

    valueDate = pNewBaseDate;
    valueDateDiscount = 1.0;

    if (baseDate != valueDate)
    {
        // NB. valueDateDiscount must be reset to 1.0 before this call
        valueDateDiscount = discountFactor(valueDate);
    }
}


// ALIB: zcsmooth.c#1973 (GtoSZCDiscount)
double ZC3ZeroCurve::discountFactor(const DateTime& date) const
{
    static const string method = "ZC3ZeroCurve::discountFactor";

    try
    {
        double adjustment = adjustmentFactor(baseDate, date);
        double discount = ZC3CurveSegment::discountFactor(data, *this, date, useSmoothing);
        discount *= adjustment;
        discount /= valueDateDiscount;
        return discount;
    }
    catch (exception& e)
    {
        throw ModelException(e, method);
    }
}


double ZC3ZeroCurve::zeroCouponRate(const DateTime& date) const
{
    static const string method = "ZC3ZeroCurve::zeroCouponRate";

    try
    {
        // going via discount factor ensures adjustments taken into account
        double discount = discountFactor(date);
        return discountToRate(discount, baseDate, date);
    }
    catch (exception& e)
    {
        throw ModelException(e, method);
    }
}


// how long is the curve ?
int ZC3ZeroCurve::length() const
{
    return data.size();
}


const DateTime& ZC3ZeroCurve::getBaseDate() const
{
    return baseDate;
}


// what is first date with real information?
const DateTime& ZC3ZeroCurve::firstDate() const
{
    return data.empty() ? baseDate : data[0].date;
}


// when does it end ?
const DateTime& ZC3ZeroCurve::endDate() const
{
    return data.empty() ? baseDate : data.back().date;
}


/** strip out the dates */
DateTimeArray ZC3ZeroCurve::getDates() const
{
    return criticalDates;
}


/** strip out the rates and dates */
CashFlowArraySP ZC3ZeroCurve::getRatesAndDates() const
{
    return CashFlowArraySP(dynamic_cast<CashFlowArray*>(datesAndRates.clone()));
}


// Converts zc-style rate into a discount factor
double ZC3ZeroCurve::rateToDiscount(double rate, const DateTime& date) const 
{
    static const string method = "ZC3ZeroCurve::rateToDiscount";

    try
    {
        return RateConversion::rateToDiscount(rate, baseDate, date, dcc.get(), basis);
    }
    catch (exception& e)
    {
        throw ModelException(e, method);
    }
}


// Converts zc-style rate into a discount factor
double ZC3ZeroCurve::rateToDiscount(double rate, double t) const 
{
    static const string method = "ZC3ZeroCurve::rateToDiscount";

    try
    {
        return RateConversion::rateToDiscountYearFrac(rate, t, basis);
    }
    catch (exception& e)
    {
        throw ModelException(e, method);
    }
}


// Converts continuous or zc-style rate into a discount factor
double ZC3ZeroCurve::rateToDiscount(
    double          rate, 
    const DateTime& startDate, 
    const DateTime& endDate, 
    bool            continuous) const 
{
    static const string method = "ZC3ZeroCurve::rateToDiscount";

    try
    {
        int basisToUse = continuous ? CompoundBasis::CONTINUOUS : basis;
        return RateConversion::rateToDiscount(rate, startDate, endDate, dcc.get(), basisToUse);
    }
    catch (exception& e)
    {
        throw ModelException(e, method);
    }
}


// Converts discount factor into a zc-style rate
double ZC3ZeroCurve::discountToRate(double discount, const DateTime& date) const 
{
    static const string method = "ZC3ZeroCurve::discountToRate";

    try
    {
        return RateConversion::discountToRate(discount, baseDate, date, dcc.get(), basis);
    }
    catch (exception& e)
    {
        throw ModelException(e, method);
    }
}


// Converts discount factor into a zc-style rate
double ZC3ZeroCurve::discountToRate(double discount, double t) const 
{
    static const string method = "ZC3ZeroCurve::discountToRate";

    try
    {
        return RateConversion::discountToRateYearFrac(discount, t, basis);
    }
    catch (exception& e)
    {
        throw ModelException(e, method);
    }
}


// Converts discount factor into a zc-style or continuous rate
double ZC3ZeroCurve::discountToRate(
    double          discount, 
    const DateTime& startDate, 
    const DateTime& endDate, 
    bool            continuous) const 
{
    static const string method = "ZC3ZeroCurve::discountToRate";

    try
    {
        int basisToUse = continuous ? CompoundBasis::CONTINUOUS : basis;
        return RateConversion::discountToRate(discount, startDate, endDate, dcc.get(), basisToUse);
    }
    catch (exception& e)
    {
        throw ModelException(e, method);
    }
}


double ZC3ZeroCurve::unadjustedZeroCouponRate(const DateTime& date) const
{
    static const string method = "ZC3ZeroCurve::unadjustedZeroCouponRate";

    int size = data.size();

    if (size < 1)
    {
        throw ModelException(method, "no points in zero curve");
    }

    return ZC3CurveSegment::unadjustedZeroCouponRate(data, *this, date, loBound, hiBound);
}


// convenience function to avoid exporting dcc
double ZC3ZeroCurve::years(const DateTime& startDate, const DateTime& endDate) const
{
    return dcc->years(startDate, endDate);
}


/** 
 * Simple Act/365F day count routine using dates as ints. 
 * 
 * Assumes ZC3ZeroCurve will always use Act/365F for dcc 
 * (as currently enforced by ZC3ZeroCurve constructors). 
 */ 
double ZC3ZeroCurve::years(int startDate, int endDate) const
{
    return (endDate - startDate) / 365.0; 
}


// ALIB: zcsmooth.c#3161
// ... used by zero curve 3 bootstrapping
void ZC3ZeroCurve::validate(bool strong) const
{
    static const string method = "ZC3ZeroCurve::validate";

    if (data.empty())
        throw ModelException(method, "No zero curve points");

    if (data[0].date.getDate() > baseDate.getDate())
    {
        string msg = Format::toString(
            "First date (%s) is after base date (%s)",
            data[0].date.toString().c_str(), 
            baseDate.toString().c_str());
        throw ModelException(method, msg);
    }

    if (data.size() > 2)
    {
        // ensure last date > first date
        if (data.back().date.getDate() <= data[0].date.getDate())
        {
            string msg = Format::toString(
                "Last date (%s) must be after first date (%s)",
                data[data.size()-1].date.toString().c_str(), 
                data[0].date.toString().c_str());
            throw ModelException(method, msg);
        }
    }

    if (strong)
    {
        string msg = Format::toString("TBD!! strong validation not yet implemented [%s line %d]", __FILE__, __LINE__);
        throw ModelException(method, msg);
    }
}


// Adds a discount factor (at a specified date) to a ZCurve
void ZC3ZeroCurve::addDiscountFactor(
    const DateTime&             date, 
    double                      disc, 
    const ZC3ZeroInterpolation& interpolation)
{
    static const string method = "ZC3ZeroCurve::addDiscountFactor";

    try
    {
        double rate = RateConversion::discountToRate(disc, baseDate, date, dcc.get(), basis);

        ZC3CurveSegment cp(date, interpolation, rate, disc, disc, data.size());
        data.push_back(cp);
    }
    catch (exception& e)
    {
        throw ModelException(e, method);
    }
}


// ALIB: szcprvt.c#590
void ZC3ZeroCurve::setExtrapolationDate(const DateTime* date)
{
    static const string method = "ZC3ZeroCurve::setExtrapolationDate";

    DateTime lastDate = data.back().date;
    extrapDate = (date && *date > lastDate) ? *date : DateTime();

    if (!extrapDate.empty())
    {
        addCriticalDate(extrapDate);
    }
}


void ZC3ZeroCurve::addCriticalDate(const DateTime& date)
{
    tmpCriticalDates.insert(date);
}


DateTimeArraySP ZC3ZeroCurve::getSegmentDates() const
{
    DateTimeArraySP dates(new DateTimeArray());
    dates->reserve(data.size());

    for (int i = 0 ; i < data.size() ; i++ )
    {
        dates->push_back(data[i].date);
    }

    return dates;
}


// ALIB: tcurve.c#545 (GtoEncapsulatedCurveSetRates)
// no additional dates added, as GtoCurveDatesAndRates does not return them
// (although  OBJECT_GET does)
void ZC3ZeroCurve::encapsulatedCurveSetRates()
{
    static const string method = "ZC3ZeroCurve::encapsulatedCurveSetRates";

    try
    {
        set<DateTime> dates(tmpCriticalDates.begin(), tmpCriticalDates.end());

        datesAndRates.resize(dates.size());
        criticalDates.resize(dates.size());
        int j = 0;
        for (set<DateTime>::iterator k = dates.begin() ; k != dates.end() ; k++)
        {
            DateTime date = *k;

            if (date == valueDate)
            {
                // not as arbitrary as seems!
                date = date.rollDate(1);
            }

            criticalDates[j] = *k;
            datesAndRates[j].date = *k;
            datesAndRates[j].amount = zeroCouponRate(date);
            j++;
        }
    }
    catch (exception& e) 
    {
        throw ModelException(e, method);
    }
}


// ALIB: szcrates.c#1059 (GtoSZCRates)
void ZC3ZeroCurve::addRates(
    const ZC3ZeroCurve&         stubCurve, 
    const ZC3RateDataArray&     ratesDataIn,
    const ZC3TurnArray&         turnsData,
    bool                        hasAdjustments,
    bool                        mtmAdjustmentFlag,    // futures flag #1
    ZC3FuturesStubRuleSP&       futuresStubRule,      // futures flag #2
    const ZC3FuturesGapRule&    futuresGapRule,       // futures flags #3&4
    const DateTime*             futuresMMDate,
    const IRVolBase*            volModelIR,
    const ZC3ZeroInterpolation& interpolation, 
    const DateTime*             extrapDate,
    const ZC3ZeroCurve*         discountCurve,
    const Holiday&              holidays)
{
    // weed out MM /futures conflicts
    ZC3RateDataArraySP ratesData = futuresStubRule->ratesUsage
        (ratesDataIn, turnsData, futuresMMDate, stubCurve, mtmAdjustmentFlag, volModelIR);

    if (mtmAdjustmentFlag && !volModelIR && hasAdjustments)
    {
        /*
         * Logic here is that if we request futures MTM adjustment but do not
         * provide a volModelIR, then we can provide explicit futures rate 
         * adjustments.
         */
        addRates(stubCurve, *ratesData, turnsData, interpolation, extrapDate, discountCurve, holidays, true);
    }
    else
    {
        /*
         * Otherwise we ignore the explicit futures rates adjustments.
         */
        addRates(stubCurve, *ratesData, turnsData, interpolation, extrapDate, discountCurve, holidays, false);

        if (mtmAdjustmentFlag && volModelIR)
        {
            if (futuresMarkToMarket(*ratesData, *volModelIR)) // changes rates in futures data
            {
                reset(NULL);
                addRates(stubCurve, *ratesData, turnsData, interpolation, extrapDate, discountCurve, holidays, false);
            }
        }
    }
}


// ALIB: szcrates.c#2713
bool ZC3ZeroCurve::futuresMarkToMarket(
    const ZC3RateDataArray& ratesData, 
    const IRVolBase&        volModelIR) const
{
    bool adjusted = false;

    for (size_t i = 0 ; i < ratesData.size() ; i++)
    {
        adjusted |= ratesData[i]->futuresMarkToMarket(*this, volModelIR);
    }

    return adjusted;
}


// ALIB: szcrates.c#2445 (SZCRates)
void ZC3ZeroCurve::addRates(
    const ZC3ZeroCurve&         stubCurve, 
    const ZC3RateDataArray&     ratesData,
    const ZC3TurnArray&         turnsData,
    const ZC3ZeroInterpolation& interpolation, 
    const DateTime*             extrapDate,
    const ZC3ZeroCurve*         discountCurve,
    const Holiday&              holidays,
    bool                        withAdjustments)
{
    static const string method = "ZC3ZeroCurve::addRates";

   /*
    * First of all we extract the intervals and adjustments from the input date.
    *
    * Then we add the adjustments into the smooth zero curve.
    *
    * Then we perform the interval merging.
    *
    * This puts us in excellent shape to generate simple cash streams and
    * bootstrap the zero curve.
    */
    ZC3IntervalArray intervals;
    ZC3StubArray stubs;

    getRatesAndAdjustments(ratesData, turnsData, intervals, stubs, discountCurve == NULL, withAdjustments);
    ZC3Interval::mergeAndSort(intervals);

    /*
     * Loop through array of inputs, with all rates expressed in a continuously 
     * compounded manner.
     */
    bool recalcLast = true;
    for (int i = 0 ; i < intervals.size() ; i++)
    {
        CashStream cs;

        if (discountCurve == NULL)
        {
            // first point already added to curve
            cs.addFixedFlow(intervals[i].getStartDate(), -1.0);
            cs.addFixedFlow(intervals[i].getEndDate(), exp(intervals[i].getContRate()));
        }
        else
        {
            /* 
             * It is possible that the discount curve has changed since the
             * previous segment was added. This may happen if we are using
             * multiple interpolation types to build a curve and we are just
             * starting a new interpolated section following a smooth segment. 
             * All we have to do is to re-add the last segment again.
             */
            if (recalcLast && prevFull.get() != NULL)
            {
                // PV = 0.0 means price is in the cash stream
                bootstrap(*prevFull, 0.0, *revert(), discountCurve, NULL);
            }

            recalcLast = false;

            /*
             * (DCC is unused as accrual period matches observation period.
             *  It is also unknown after rates have been converted into intervals)
             */
            cs.addFloatFlow(intervals[i].getEndDate(), -1.0, 0.0, NULL,
                intervals[i].getStartDate(), intervals[i].getEndDate(), 
                intervals[i].getStartDate(), intervals[i].getEndDate());
            cs.addFixedFlow
                (intervals[i].getEndDate(), exp(intervals[i].getContRate()) -1.0);
        }

        // PV = 0.0 means price is in the cash stream
        bootstrap(cs, 0.0, interpolation, discountCurve, NULL);
    }

    addStub(stubs, *ZC3Stub::STUB_DCC, ZC3Stub::STUB_BASIS);
    setExtrapolationDate(extrapDate);
}


// VC6 Release configuration fix for unresolved lessThan<ZC3Turn>
static inline bool lessThan_ZC3Turn(const refCountPtr< ZC3Turn >& first, const refCountPtr< ZC3Turn >& second)
{
    return lessThan(first, second); 
}
// ALIB: szcrates.c#2818
void ZC3ZeroCurve::getRatesAndAdjustments(
    const ZC3RateDataArray& ratesData, 
    const ZC3TurnArray&     turnsData,
    ZC3IntervalArray&       intervals, 
    ZC3StubArray&           stubs,
	bool                    withBasis,
    bool                    withAdjustment)
{
    static const string method = "ZC3ZeroCurve::getRatesAndAdjustments";

    /*
     * All this loop does is convert the single big array of inputs into
     * sub-arrays of rates and adjustments, with all rates expressed in a
     * continuously compounded manner.
     */
    for (size_t i = 0 ; i < ratesData.size() ; i++)
    {
        ZC3RateData& instrument = *ratesData[i];

        if (instrument.getStartDate() == instrument.getEndDate())
        {
            continue;   // ignore this no-information case
        }
        else if (instrument.getStartDate() > instrument.getEndDate())
        {
            string msg = Format::toString(
                "Start date (%s) after end date (%s) for element %d",
                instrument.getStartDate().toString().c_str(),
                instrument.getEndDate().toString().c_str(),
                i);
            throw ModelException(method, msg);
        }
        else
        {
            if (instrument.isStub(baseDate))
            {
                // this is a stub rate
                instrument.prepareStub(stubs, withBasis);
            }
            else
            {
                ZC3Interval interval(instrument.getStartDate(), 
					                 instrument.getEndDate(), 
					                 instrument.getContinuousRate(withBasis, withAdjustment));
                intervals.push_back(interval);
            }
        }
    }

    // handle adjustments (NB. must cast away const-ness for sort)
    ZC3TurnArray& turnsData2 = const_cast<ZC3TurnArray&>(turnsData);
    sort(turnsData2.begin(), turnsData2.end(), lessThan_ZC3Turn);

    for (size_t j = 0 ; j < turnsData.size() ; j++)
    {
        turnsData[j]->addAdjustment(*this);
    }

    if (intervals.empty() && adjustments.empty() && adjAbsolute.empty())
    {
        throw ModelException(method, "No rates specified");
    }
}


// ALIB: szcprvt.c#288 (SZCAddAdjustmentToList)
void ZC3ZeroCurve::addAdjustmentToList(const ZC3Adjustment& adjustment)
{
    static const string method = "ZC3ZeroCurve::addAdjustmentToList";

    ZC3AdjustmentArray& list = adjustment.isAbsolute() ? adjAbsolute : adjustments;

    /*
     * Adjustments will be put into increasing date order by this function.
     *
     * Adjustments must not overlap (slightly arbitrary constraint).
     */
    list.push_back(adjustment);
    sort(list.begin(), list.end(), ZC3Interval::lessThan);

    for (int i = 1 ; i < list.size() ; i++)
    {
        if (list[i-1].getEndDate() > list[i].getStartDate())
        {
            string msg = Format::toString(
                "New adjustment for %s overlaps with existing adjustment for %s", 
                list[i].toString().c_str(), list[i-1].toString().c_str());
            throw ModelException(method, msg);
        }
    }
}


// ALIB: szcprvt.c#805 (GtoSZCCheckStubs)
void ZC3ZeroCurve::checkStubs(const ZC3StubArray& stubs) const
{
    static const string method = "ZC3ZeroCurve::checkStubs";

    if (!stubs.empty())
    {
        int size = stubs.size();

        for (int i = 0 ; i < size - 1 ; i++)
        {
            // end date for this stub must be start date for next stub.
            if (!stubs[i].isContiguousWith(stubs[i+1]))
            {
                string msg = Format::toString(
                    "Rate from %s is not contiguous with rate from %s",
                    stubs[i].toString().c_str(),
                    stubs[i+1].toString().c_str());
                throw ModelException(method, msg);
            }
        }

        // end date for last stub must be base date of zero curve
        if (stubs[size-1].getEndDate() != baseDate)
        {
            string msg = Format::toString(
                    "Rate from %s is not contiguous with base date of zero curve [%s]",
                    stubs[size-1].toString().c_str(),
                    baseDate.toString().c_str());
            throw ModelException(method, msg);
        }
    }
}


// ALIB: szcbuild.c#418 (SZCAddStubValidateInputs)
void ZC3ZeroCurve::addStubValidateInputs(
    const ZC3StubArray& stubs, 
    const string&       method) const
{
    if (!stubs.empty())
    {
        int size = stubs.size();

        for (int i = 1 ; i < size ; i++)
        {
            if (stubs[i-1].getStartDate() >= stubs[i].getStartDate())
            {
                string msg = Format::toString(
                    "Start dates must be in strictly ascending order (%s) >= (%s)",
                    stubs[i-1].getStartDate().toString().c_str(),
                    stubs[i].getStartDate().toString().c_str());
                throw ModelException(method, msg);
            }
        }

        if (stubs[size-1].getStartDate() >= baseDate)
        {
            string msg = Format::toString(
                "Start dates must be before zero curve base date (%s) >= (%s)",
                stubs[size-1].getStartDate().toString().c_str(),
                baseDate.toString().c_str());
            throw ModelException(method, msg);
        }
    }

    if (data.empty())
    {
        throw ModelException(method,"No rates in zero curve");
    }

    if (baseDate != data[0].date)
    {
        string msg = Format::toString(
            "Base date must be first date in the zero curve: (%s) does not equal (%s)",
            baseDate.toString().c_str(),
            data[0].date.toString().c_str());
        throw ModelException(method, msg);
    }
}


// ALIB: szcbuild.c#227 (GtoSZCAddStub)
void ZC3ZeroCurve::addStub(
    const ZC3StubArray&       stubs, 
    const DayCountConvention& rateDayCountConv, 
    int                       rateBasis)
{
    static const string method = "ZC3ZeroCurve::addStub";

    // no rates added should be acceptable
    if (stubs.empty())
        return;

    // validate stub data
    checkStubs(stubs);
    addStubValidateInputs(stubs, method);

    /*
     * This is the only function which inserts segments at the beginning of the
     * curve.
     *
     * The modified curve has stubs.size() more dates than the original curve. 
     * These dates are at the beginning of the curve.
     */
    ZC3ZeroInterpolationConstSP interpolation = ZC3ZeroInterpolation::make("Flat");
    DateTime endDate = baseDate;
    double cumulativeDiscount = 1.0;

    for (int i = stubs.size() - 1 ; i >= 0 ; --i)
    {
        const DateTime& startDate = stubs[i].getStartDate();

        double discount = RateConversion::rateToDiscount(stubs[i].getStubRate(), 
            startDate, endDate, &rateDayCountConv, rateBasis);
        cumulativeDiscount *= discount;
        double zeroRate = discountToRate(cumulativeDiscount, startDate, baseDate);

        /*
         * The choice of flat forwards rather than linear interpolation is to
         * agree most closely with ChaseTools, which as a rule does not have 
         * rates between today and the spot date, and so as a rule these are 
         * interpolated using flat forwards - it doesn't really matter.
         *
         * (The point indexes are reset in ZC3CurveSegment::calculateAverageRates).
         */
        ZC3CurveSegment point(startDate, *interpolation, zeroRate, 1.0/cumulativeDiscount, 0.0, 0);
        data.insert(data.begin(), point);

        addCriticalDate(startDate);
        endDate = startDate;
    }

    /*
     * Go back and calculate the average rates. Needs the discount factor at
     * the start of the segment. We could have done this in the previous loop,
     * but we already have the code to do it in the segment member function
     * averageRate.
     */
    ZC3CurveSegment::calculateAverageRates(data, *this, stubs.size());
}


// VC6 Release configuration fix for unresolved lessThan<ZC3SwapData>
static inline bool lessThan_ZC3SwapData(const refCountPtr< ZC3SwapData >& first, const refCountPtr< ZC3SwapData >& second)
{
    return lessThan(first, second); 
}
// ALIB: szcswaps.c#377 (GtoSZCSwaps2 combined with GtoSZCSwapsWithBasis2 at szcswaps.c#1307)
void ZC3ZeroCurve::swaps2(
    const ZC3ZeroCurve*           discountCurve,
    const ZC3SwapDataArray&             swapsDataIn,
    const BadDayConvention&       badDayConv,
    const StubPlacement&          stubPos,
    const Holiday&                holidays,
    const Holiday&                basisHolidays,
    const ZC3CouponInterpolation* couponInterp,
    const ZC3ZeroInterpolation&   zeroInterpolation,
    const ZC3FuturesGapRule&      gapRule,
    const DateTime*               extrapDate,
    const bool                    convDelayAdj,
    const IRVolBase*              volModelIR,
    ZC3RecurseInfo*               recurseInfo,
    const bool                    withBasis)
{
    static const string method = "ZC3ZeroCurve::swaps2";

    if ( convDelayAdj)
    {
        throw ModelException(method, "Convexity & delay adjustments not yet supported.");
    }

    ZC3ZeroInterpolationConstSP zeroInterp(&zeroInterpolation);

    /*
     * For backward compatibility we currently allow "money-market" rates to be
     * added as swaps via these routines. Some of these money-market rates can
     * even be "stub" rates before the base date of the zero curve.
     */
    bool recursive = false;
    if (recurseInfo)
    {
        if (couponInterp)
        {
            throw ModelException(method, "Cannot use coupon interpolation with single step spline interpolation");
        }

        if (discountCurve)
        {
            zeroInterp = ZC3ZeroInterpolation::make("Flat", DateTime(), true);
        }
        else
        {
            zeroInterp = ZC3ZeroInterpolation::make("Flat");
            recursive = true;
            recurseInfo->reset();
        }
    }

    zeroInterpolation.validateInterpType(couponInterp);

    // Validation in ALIB szcswaps.c#485 and #582 is not required due to ZC3SwapData
    // constructor defaults

    /*
     * If we are using coupon interpolation then get a larger set of swap rates 
     * by coupon interpolation.
     *
     * Note that coupon interpolation is incompatible with using adjustments. 
     * This is validated by getInterpRates().
     */
    ZC3SwapDataArraySP swapsDataInterpolated;
    if (couponInterp != NULL)
    {
        swapsDataInterpolated = couponInterp->getInterpRates(
            swapsDataIn,
            discountCurve ? *discountCurve : *this,
            this, 
            badDayConv,
            stubPos,
            holidays,
            basisHolidays,
            convDelayAdj,
            volModelIR,
            withBasis);
    }

    ZC3SwapDataArray& swapsData = swapsDataInterpolated.get() 
        ? *swapsDataInterpolated : const_cast<ZC3SwapDataArray&>(swapsDataIn);

    const DateTime lastZeroDate = endDate();    

    sort(swapsData.begin(), swapsData.end(), lessThan_ZC3SwapData);

    /*
     * See if we need to make two loops.
     */
    bool twoLoops = false;
    for (size_t j = 0 ; j < swapsData.size() ; j++)
    {
        twoLoops |= swapsData[j]->isForPass(2);
    }

    ZC3ZeroCurveSP szcDiscFirst(const_cast<ZC3ZeroCurve*>(discountCurve));
    ZC3ZeroCurveSP szcEstFirst;
    ZC3ZeroCurve* szcCurrent = this;

    if (twoLoops)
    {
        if (withBasis && discountCurve)
        {
            szcDiscFirst = ZC3ZeroCurveSP(const_cast<ZC3ZeroCurve*>(discountCurve));
            szcEstFirst.reset(new ZC3ZeroCurve(*this));
            szcCurrent = szcEstFirst.get();
        }
        else
        {
            szcDiscFirst.reset(new ZC3ZeroCurve(*this));
            szcEstFirst = szcDiscFirst;
            szcCurrent = szcDiscFirst.get();
        }
    }
    
    bool finalLoop = !twoLoops;
    bool anotherLoop = true;
    bool secondLoop = false;
    ZC3StubArray stubs;

    /*
     * This code logic looks a little odd, so please let me explain what is going 
     * on here. 
     *
     * Essentially we can go through the list of swaps once or twice. We do it
     * twice if there are any swaps to be added in the second loop. Adding a swap
     * in the second loop involves computing its swap rate from the initial zero
     * curve and then adding the adjustment.
     *
     * We control this by using a loop which can in fact only be processed at
     * most twice. The alternative would be to have a function which added
     * swaps, and two bits of code which decided which swaps should be set-up
     * for calling that particular function. I think that the approach chosen
     * is just as good, and involves setting up less intermediate arrays.
     *
     * The double looping is controlled by the three booleans - anotherLoop,
     * finalLoop, secondLoop. Naturally anotherLoop starts off as TRUE and
     * secondLoop as FALSE. finalLoop can be TRUE or FALSE. There are some
     * things we only do on the finalLoop - namely add any rates before the
     * start of the zero curve (these are not really swaps, but supported by
     * this routine for backward compatibility with ZERO_CURVE_SWAPS2).
     *
     * There are also some things we only do on the secondLoop - namely add
     * some of the artificial swap rates.
     *
     * At the end of the loop, we re-set these three booleans, and also
     * perform some other processing. You can see this code at the end of
     * the while(anotherLoop) {...} block of code.
     *
     * Each swap (real or artificial) is added using addOneSwap (which is a
     * private function).
     *
     * In one sense, coupon interpolation is not optimized - we go through the
     * same routine each time (addOneSwap which calls bootstrap), even
     * though it might be possible to "remember" intermediate results from the
     * previous function call. However bootstrap would not have to iterate
     * to find zero rates, which is a small saving in calculation time.
     *
     * The upshot is that the new zero curve is slower than the old zero curve
     * for adding swaps using coupon interpolation.
     */
    while(anotherLoop)
    {
        for (size_t i = 0 ; i < swapsData.size() ; i++)
        {
            /*
             * for efficiency we do not calculate the estimating curve within
             * the recursive loop.
             */
            if (discountCurve && recursive)
            {
                break;
            }

            double couponRate = 0.0;

            // Exclude swaps that do not extend curve
            if (swapsData[i]->isExcluded(lastZeroDate, baseDate))
            {
                continue;
            }
            else if (swapsData[i]->isForPass(1) && !secondLoop)
            {
                // process this swap
                if (gapRule.ignore(*swapsData[i], *szcCurrent, *this, badDayConv, holidays))
                {
                    continue;
                }

                couponRate = swapsData[i]->getRate();
            }
            else if (swapsData[i]->isForPass(2)) // process this swap on the second pass only
            {
                if (!secondLoop)
                {
                    continue;
                }

                /*
                 * We get the swap rate from the first zero curve adding the adjustment.
                 */
                couponRate = swapsData[i]->getCouponRate
                    (*szcDiscFirst, szcEstFirst.get(), 
                    convDelayAdj, volModelIR, stubPos, badDayConv, holidays, withBasis);
            }

            // is it a stub?
            if (swapsData[i]->isStub(baseDate))
            {
                if (finalLoop)
                {
                    swapsData[i]->prepareStub(stubs, couponRate, badDayConv, holidays, withBasis);
                }

                continue;
            }

            // NB. Chase style interpolation not implemented
            szcCurrent->addOneSwap(discountCurve, *(swapsData[i]),
                couponRate, badDayConv, stubPos, holidays, basisHolidays, *zeroInterp, 
                convDelayAdj, volModelIR, recursive, recurseInfo, withBasis);
        }

        if (!stubs.empty() && finalLoop)
        {
             // Only add stubs on the final pass.
             szcCurrent->addStub(stubs, *ZC3Stub::STUB_DCC, ZC3Stub::STUB_BASIS);
        }

        if (finalLoop)
        {
            anotherLoop = false;
        }
        else
        {
            anotherLoop = true;
            finalLoop = true;
            secondLoop = true;

            // Set curves for the second loop.
            szcCurrent = this;
        }
    }

    if (this != szcCurrent)
    {
        throw ModelException(method, "Program bug");
    }

    setExtrapolationDate(extrapDate);
}


void ZC3ZeroCurve::extendLastRate(
    const DateTime& lastZeroDate, 
    const DateTime& adjusted)
{
    double zcLastRate = data[data.size() - 1].avgRate;

    addRate(zcLastRate, lastZeroDate, adjusted, 
        *dcc, CompoundBasis::CONTINUOUS, *ZC3ZeroInterpolation::make("Flat"));
}


// ALIB: szcswaps.c#956 (GtoSZCAddOneSwap) & 2168 (GtoSZCAddOneSwapWithBasis) & 2377 (GtoSZCAddOneSwapEstimating)
// NB. in ALIB the same BDC & stub position is used for basis leg
void ZC3ZeroCurve::addOneSwap(
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
    const bool                  withBasis)
{
    static const string method = "ZC3ZeroCurve::addOneSwap";

    try
    {
        CashStreamSP    spCashStream;
        CashFlowArraySP spCfl;
		CashStream*     cashStream;
		CashFlowArray*  cfl;

        bool estimating = withBasis && discountCurve != NULL;

        // For efficiency we cache cfl's and cash streams.
        if (!estimating && recursive && !recurseInfo->isFirstPass())
        {
            recurseInfo->getCurrent(cashStream, cfl);
        }
        else
        {
            if (convDelayAdj)
            {
                throw ModelException(method,"Convexity and delay adjustments are not yet supported");
            }

            useSmoothing = zeroInterp.isSmoothed();
            spCfl = swapData.getKnownCashflows(*this, couponRate, stubPos, 
                badDayConv, holidays, basisHolidays, withBasis, estimating);
			cfl = spCfl.get();
            useSmoothing = true;

            spCashStream = CashStreamSP(new CashStream());
    	    cashStream = spCashStream.get();
            cashStream->add(*cfl);

            if (withBasis && !estimating)
            {
                if (recursive)
                {
                    /*
                     * Note that as at this point only fixed flows are in the 
                     * cash stream the interpolation type does not matter.
                     */
                    spCfl = cashStream->flows(*this, volModelIR);
			        cfl = spCfl.get();
                }
            }
            else
            {
                swapData.getFloatingCashflows(
                    *cashStream, couponRate, stubPos, badDayConv, holidays, estimating);
            }
        }

        // szcswaps.c#1155
        if (!estimating && recursive && !recurseInfo->isFirstPass())
        {
            recurseInfo->addNonBenchmarkDiscounts(*this, *cfl);
        }

        /*
         * Add the swap to the zero curve.
         */
        safeBootstrap(*cashStream, 0.0, zeroInterp, discountCurve, NULL);

        if (recursive)
        {
            recurseInfo->next(*cashStream, *cfl);
        }
    }
    catch(exception& e)
    {
        string msg = Format::toString("Could not add swap ending on %s", 
            swapData.getEndDate().toString().c_str());
        throw ModelException(e, method, msg);
    }
}


// ALIB: szcswaps.c#5234 (SZCBootstrap)
void ZC3ZeroCurve::safeBootstrap(
    CashStream&                 cashStream,
    double                      presentValue,
    const ZC3ZeroInterpolation& interpolation,
    const ZC3ZeroCurve*         discountCurve,
    const IRVolBase*            volModelIR)
{
    int closeEnough = 5;

    DateTime lastZeroDate = data[data.size() - 1].date;
    DateTime maturityDate = cashStream.maturityDate();

    if ((maturityDate <= lastZeroDate)
     && (maturityDate >= lastZeroDate.rollDate(-1 * closeEnough)))
    {
        return;
    }
    else
    {
        /*
         * if maturityDate < lastZeroCurve - closeEnough then
         * bootstrap fails anyway - but we will get an error
         * message out of bootstrap by calling the function.
         */
        bootstrap(cashStream, presentValue, interpolation, discountCurve, volModelIR);
    }
}


// ALIB: szcbuild.c#748 (GtoSZCBootstrap)
void ZC3ZeroCurve::bootstrap(
    CashStream&                 cashStream,
    double                      presentValue,
    const ZC3ZeroInterpolation& interpolation,
    const ZC3ZeroCurve*         discountCurve,
    const IRVolBase*            volModelIR)
{
    static const string method = "ZC3ZeroCurve::bootstrap";

    validate(false);

    if (noBootstrap)
    {
        throw ModelException(method, "Cannot bootstrap further with this curve");
    }

    if (volModelIR)
    {
        string msg = "Convexity and delay adjustments not supported. "
                     "You must set the volModelIR to NULL.";
        throw ModelException(method, msg);
    }

    // inputs are validated so we can safely calculate the last date of the zero curve
    ZC3CurveSegment* prevSegment = &(data[data.size() - 1]);
    bool smoothed = !prevSegment->isSmoothed() && !interpolation.isSmoothed();

    /*
     * Handle the part of the cash flow list already covered by the zero curve.
     *
     * We do not allow cash flows *before* the zero curve base date. These can
     * be handled by using the AddStubs routine.
     *
     * We do not care about the temporal order of the cash stream that we are
     * adding.
     *
     * We do care that the maturity date of this cash stream extends the zero
     * curve.
     */

    double pvKnown = 0.0;
    DateTime maturityDate(0,0);

    /*
     * Add to critical date list so that we can re-price benchmarks by
     * simply interpolating rates from the smooth zero curve and using
     * these to construct a set of dates and rates.
     *
     * The payment date and the observed rate start and end dates are
     * both critical for this purpose.
     */
    DateTimeArraySP dates = cashStream.getCriticalDates();
    for (int j = 0 ; j < dates->size() ; j++)
    {
        addCriticalDate((*dates)[j]);
    }

    CashStreamSP csp(new CashStream());
    for (int i = 0 ; i < cashStream.size() ; i++ )
    {
        /*
         * Now we see whether the cash flow is already covered by the zero
         * curve, in which case we price it now, or else it is not covered,
         * in which case we keep track of it for the subsequent bootstrap.
         */
        const DateTime& csMaturityDate = cashStream[i].getMaturityDate();

        if (cashStream[i].getPayDate() < baseDate)
        {
            string msg = Format::toString(
                "Cash flow at date %s before zero curve base date %s",
                cashStream[i].getPayDate().toString().c_str(),
                baseDate.toString().c_str());
            throw ModelException(method, msg);
        }
        else if (csMaturityDate <= prevSegment->date 
             || (cashStream[i].isFixed() && discountCurve != NULL))
        {
            // Calculate cash stream amount without smoothing
            setUseSmoothing(smoothed);
            double amount = cashStream[i].getAmount(*this, volModelIR);
            setUseSmoothing(true);

            if (!Maths::isZero(amount))
            {
                double discount;

                if (discountCurve != NULL)
                {
                    discount = discountCurve->discountFactor(cashStream[i].getPayDate());
                }
                else
                {
                    setUseSmoothing(smoothed);
                    discount = discountFactor(cashStream[i].getPayDate());
                    setUseSmoothing(true);
                }

                pvKnown += discount * amount;
            }
        }
        else
        {
            if (!cashStream[i].isZero())
            {
                if (csMaturityDate > maturityDate)
                {
                    maturityDate = csMaturityDate;
                }
                
                csp->add(const_cast<RiskFreeCashFlow*>(&(cashStream[i])));
            }
        }
    }

    if (csp->empty())
    {
        string msg = Format::toString(
            "Cash stream does not extend beyond last zero date %s",
            prevSegment->date.toString().c_str());
        throw ModelException(method, msg);
    }

    /*
     * Calculate PV of cash stream at prevSegment date.  We can calculate PV at
     * base date and divide this by discount at prevSegment date.
     */
    setUseSmoothing(smoothed);
    double prevDiscount = discountFactor(prevSegment->date);
    setUseSmoothing(true);

    double pvRemainder;

    if (discountCurve != NULL)
    {
        /*
         * When we iterate we will be computing discount factors to the base
         * date of the curve using the discount curve.
         */
        pvRemainder = presentValue - pvKnown;
    }
    else
    {
        /*
         * When we iterate we will be computing discount factors to the start
         * of the new segment.  Hence we divide by the previous discount factor
         * to get the pvRemainder as the PV of the remainder at the start of
         * the new segment.
         */
        pvRemainder = (presentValue - pvKnown) / prevDiscount;
    }

    /*
     * Now we calculate the rates for the new segment. If there is more than
     * one cash flow in csp this requires iteration. Otherwise we can perform
     * an exact calculation. Note that only two interpolation types are supported.
     *
     * This block of code should undoubtedly be moved into its own function!
     *
     * As a first step we add an empty segment to the zero curve. Subsequently
     * when we amend segment we simultaneously amend the zero curve.
     */
    if (maturityDate <= baseDate)
    {
        string msg = Format::toString("program bug [%s line #%d]", __FILE__, __LINE__);
        throw ModelException(method, msg);
    }

    /*
     * When building a Chase style index curve (with ZeroCurve3) we use the
     * shape of the discount curve. Effectively we add multiple segments which
     * are the rates in the discount curve plus a spread.
     */
    if (interpolation.isShaped())
    {
        // NB. discount curve is not NULL if got here
        addShapedSegments(interpolation, *discountCurve, *csp, pvRemainder, volModelIR);
        return;
    }

    ZC3CurveSegment& segment = addSegment(maturityDate, interpolation, prevSegment);
    prevPoint = ZC3CurveSegmentSP(new ZC3CurveSegment(*prevSegment));

    /*
     * Check to see if there any absolute adjustments in the interval
     *    [prevSegment date, segment date]
     * If they exist, then we will perform more complex calculations in
     * the objective functions.
     *
     * Also if they exist, it prevents closed form solutions.
     *
     * Finally if they exist, then they must be included wholly within
     * the segment.
     */
    bool foundAdjAbsolute = findAbsoluteAdjustments(*prevSegment, segment);

    if (csp->size() == 1 && (*csp)[0].isFixed() && !foundAdjAbsolute)
    {
        if (discountCurve != NULL)
        {
            string msg = Format::toString("program bug [%s line %d]", __FILE__, __LINE__);
            throw ModelException(method, msg);
        }

        if ((*csp)[0].getPayDate() != maturityDate)
        {
            string msg = Format::toString("program bug [%s line %d]", __FILE__, __LINE__);
            throw ModelException(method, msg);
        }

        /*
         * adjustment is the discount adjustment between prevSegment->date and
         * maturityDate.
         *
         * discount is the discount factor between prevSegment->date and
         * maturityDate, taking adjustments into account.
         *
         * pvRemainder is the pv of the remaining cash stream at
         * prevSegment->date.
         *
         * We combine these to get segment->discount, and depending upon
         * the interpolation type we calculate the type specific data.
         */
        double adjustment = adjustmentFactor(prevSegment->date, maturityDate);
        setUseSmoothing(smoothed);
        double amount = (*csp)[0].getAmount(*this, volModelIR);
        setUseSmoothing(true);
        double discount = pvRemainder / (amount * adjustment);

        segment.discount = prevSegment->discount * discount;

        /*
         * Regardless of interpolation type, rate is the zero rate
         * corresponding to discount within this segment.
         */
        segment.rate = discountToRate(segment.discount, maturityDate);

        /*
         * The average rate is easy since we know the discount factor
         * between the last zero date and the maturity date. The only
         * trick to remember is that the average rate is continuously
         * compounded.
         */
        segment.avgRate = interpolation.discountToRate
            (*this, discount, prevSegment->date, maturityDate);
    }
    else
    {
        /*
         * Iteration is required for more than one cash flow or floating flows
         * or for absolute adjustments.
         *
         * For linear interpolation we need to get the rate.
         * For flat forwards interpolation we need to get the avgRate.
         */
        interpolation.solve(*csp, *this, discountCurve, volModelIR, *prevSegment, segment,
            presentValue, pvRemainder, pvKnown, foundAdjAbsolute, smoothed, prevSegment->discount);
    }

    if (prevSegment->date == baseDate)
    {
        /*
         * This ensures that the rate at the base date is defined,
         * even though on its own account it is meaningless.
         *
         * However we need to do this to ensure linear interpolation
         * works as expected in the interval between the base date and
         * the first date in the curve.
         */
        prevSegment->rate = segment.rate;
        prevSegment->averageRate(*this, true);
    }

    // if required, smooth zero curve
    if (prevSegment->isSmoothed() || segment.isSmoothed())
    {
        pvRemainder = smoothing(cashStream, presentValue, 
            discountCurve, *csp, volModelIR, foundAdjAbsolute);
    }

    /*
     * If we are building an estimating curve we may need to save a copy of
     * the input cash stream as well.
     */
    if (discountCurve)
    {
        prevFull.reset(new CashStream(cashStream));
        prevCs.reset(cs.get());
        prevPv = PV;
    }
    else
    {
        prevFull.reset();
        prevCs.reset(new CashStream());
        prevPv = 0.0;
    }

    cs = csp;
    PV = pvRemainder;

    // done!!
}


// ALIB: szcbuild.c#2461
void ZC3ZeroCurve::addShapedSegments(
    const ZC3ZeroInterpolation& interpolation,
    const ZC3ZeroCurve&         discountCurve, 
    CashStream&                 cashStream, 
    const double                pvRemainder, 
    const IRVolBase*            volModelIR)
{
    static const string method = "ZC3ZeroCurve::addShapedSegments";

    // first get critical points from discount curve
    int numDates = discountCurve.criticalDates.size();

    int oldSize = data.size();
    ZC3CurveSegment* prevSegment = &data[oldSize - 1];
    if (prevSegment->isSmoothed())
    {
        throw ModelException(method, "Cannot mix shaped flat forwards with smoothing");
    }

    DateTime csMatDate = cashStream[cashStream.size() - 1].getMaturityDate();
    DateTimeArray shapeStartDts(numDates);
    DateTimeArray shapeEndDts(numDates);
    DoubleArray shapeRates(numDates);
    BoolArray   useAdjAbsolutes(numDates);

    // get the shape from the discount curve
    shapeStartDts[0] = prevSegment->date;
    int numSegs = 0;
    for (int i = 0 ; i < numDates ; i++)
    {
        DateTime date = discountCurve.criticalDates[i];
        if (date > prevSegment->date && date < csMatDate)
        {
            /*
             * of course we cannot have segments ending within an absolute 
             * adjustment...
             */
            bool ignoreThisDate = false;
            if (adjAbsolute.size() > 0)
            {
                for (int idx = 0 ; idx < adjAbsolute.size() ; idx++)
                {
                    ZC3Adjustment& adjustment = adjAbsolute[idx];

                    if (adjustment.getStartDate() >= shapeStartDts[numSegs]
                     && adjustment.getStartDate() < date)
                    {
                        if (adjustment.getEndDate() > date)
                        {
                            ignoreThisDate = true;
                        }
                    }
                }
            }

            if (ignoreThisDate)
            {
                continue;
            }

            shapeEndDts[numSegs] = date;
            addCriticalDate(date);
            shapeStartDts[numSegs+1] = date;
            shapeRates[numSegs] = discountCurve.fwd(shapeStartDts[numSegs], 
                shapeEndDts[numSegs], dcc.get(), CompoundBasis::CONTINUOUS);
            numSegs++;
        }

        if (date >= csMatDate)
        {
            break;  // all done
        }
    }

    shapeEndDts[numSegs] = csMatDate;
    shapeRates[numSegs] = discountCurve.fwd(shapeStartDts[numSegs], 
        shapeEndDts[numSegs], dcc.get(), CompoundBasis::CONTINUOUS);
    numSegs++;
    
    // add the new segments to the curve
    for (int j = 0 ; j < numSegs ; j++)
    {
        ZC3CurveSegment& segment = addSegment(shapeEndDts[j], interpolation, prevSegment);
        useAdjAbsolutes[j] = findAbsoluteAdjustments(*prevSegment, segment);
    }

    prevSegment = &data[oldSize - 1];
    ZC3CurveSegment* firstSeg = &data[oldSize];

    /*
     * For flat forwards interpolation we will calculate times from the 
     * previous segment date.  This is because we are simply trying to 
     * calculate the average rate since the last time point in the zero
     * curve.  Thus within the iteration we can use the pv of the remaining
     * cash flows at the previous time point in the zero curve.
     */

    ShapedFlatFwdObjFunc func(*this, cashStream, discountCurve, pvRemainder, 
        volModelIR, useAdjAbsolutes, numSegs, shapeRates, *prevSegment);
    double spread = func.solve(0.0);    // first guess with zero spread

    // we have the spread - now fully populate the new segments
    populateShapedCurve(oldSize, numSegs, shapeRates, spread, useAdjAbsolutes, true);

    if (prevSegment->date == baseDate)
    {
        /*
         * This ensures that the rate at the base date is defined, even though
         * on its own account it is meaningless.
         *
         * However, we need to do this to ensure that linear interpolation 
         * works as expected in the interval between the base date and the 
         * first date in the curve.
         */
        prevSegment->rate = firstSeg->rate;
        prevSegment->averageRate(*this, true);
    }

    prevFull.reset(new CashStream());
    PV = pvRemainder;
    cs.reset(&cashStream);
    prevCs.reset(new CashStream());
    prevPv = 0.0;
}


// ALIB: szcbuild.c#2721 (PopulateShapedCurve)
void ZC3ZeroCurve::populateShapedCurve(
    const int          oldSize, 
    const int          numSegs, 
    const DoubleArray& shapeRates, 
    const double       spread, 
    const BoolArray&   useAdjAbsolutes, 
    const bool         calcZero)
{
    for (int i = 0; i < numSegs ; i++)
    {
        ZC3CurveSegment& segment = data[oldSize + i];
        ZC3CurveSegment& prevSegment = data[oldSize + i - 1];

        segment.avgRate = shapeRates[i] + spread;
        if (useAdjAbsolutes[i])
        {
            flatFwdAdjAbsolute(prevSegment, segment);
        }
        else
        {
            /*
             * remove relative adjustments to stop them being included twice -
             * since they are already included in the shape from the discount 
             * curve.
             */
            double adjustment = adjustmentFactor(prevSegment.date, segment.date);
            double yearFrac = dcc->years(prevSegment.date, segment.date);
            segment.avgRate += log(adjustment) / yearFrac;
        }

        double discount = rateToDiscount(segment.avgRate, prevSegment.date, segment.date, true);
        segment.discount = discount * prevSegment.discount;

        if (calcZero)
        {
            segment.rate = discountToRate(segment.discount, segment.date);
        }
    }
}


// ALIB: szcbuild.c#1731
double ZC3ZeroCurve::smoothing(
    const CashStream& cashStream,       // cash stream to add (all inclusive)
    double            presentValue,     // of cash stream at curve base date
    const ZeroCurve*  discountCurve,
    const CashStream& cs,               // remaining cash flows
    const IRVolBase*  volModelIR,
    bool              foundAdjAbs)
{
    static const string method = "ZC3ZeroCurve::smoothing";

    ZC3CurveSegment& lastSegment = data.back();
    ZC3CurveSegment& prevSegment = data[data.size() - 2];

    // We might not need to do any smoothing at all
    if (!lastSegment.isSmoothed() && !prevSegment.isSmoothed())
    {
        throw ModelException(method, "No smoothing required");
    }

    /*
     * We always need the average rate for the last segment.
     */
    double lastAvgRate = lastSegment.avgRate; // Average forward rate for last segment

    /*
     * At least one of the segments needs smoothing.
     * We always need the forward rate at the beginning of the last segment.
     */
    double  startFwdRate = 0.0;                // Forward rate at start of last segment
    if (prevSegment.date == baseDate)
    {
        /*
         * If the previous segment is the base date segment, then 
         * (a) It cannot usually compute an average rate.
         * (b) We need it to be flat in any case.
         *
         * Hence our start forward rate is simply the average rate for the last
         * segment.
         */
        startFwdRate = lastAvgRate;
    }
    else
    {
        if (data.size() < 3)
        {
            /*
             * This is impossible, because the array size must be at least 2,
             * and when it is precisely 2, the previous segment must be the
             * base date segment.
             */
            throw ModelException(method, "Program bug. Array size must be >= 3");
        }

        ZC3CurveSegment& prevSegment2 = data[data.size() - 3];

        /*
         * This is the duration weighting. The calculation of the duration is
         * somewhat arbitrary (uses USD swap conventions).
         */
        prevSegment2.segmentDuration(*this);
        prevSegment.segmentDuration(*this);
        lastSegment.segmentDuration(*this);

        double t0 = prevSegment2.duration; // Duration at beginning of previous segment
        double t1 = prevSegment.duration;  // Duration at beginning of last segment
        double t2 = lastSegment.duration;  // Duration at end of last segment

        double prevAvgRate = prevSegment.avgRate;

        startFwdRate = (lastAvgRate * (t1-t0) + prevAvgRate * (t2-t1)) / (t2-t0);
    }

    /*
     * Now complete the smoothing for the previous segment.
     * This is where we calculate the parabola.
     */
    prevSegment.smooth(*this, startFwdRate, discountCurve, volModelIR);

    /*
     * Now that we have completed the smoothing for the previous segment we
     * are in position to calculate the pv of the cash flows using the smooth
     * forwards method for all flows <= last date of the zero curve.
     *
     * Note that this loop is similar to the loop in bootstrap where
     * we are dealing with the unsmoothed bootstrap.
     */
    double pvKnown = 0.0;
    for (int i = 0 ; i < cashStream.size() ; i++)
    {
        DateTime csMaturityDate = cashStream[i].getMaturityDate();
        if (cashStream[i].getPayDate() < baseDate)
        {
            /*
             * This cannot happen since the loop in bootstrap would have
             * caught it.
             */
            throw ModelException(method, "Program bug. cashStream[i].getPayDate() < base date");
        }
        else if (csMaturityDate <= prevSegment.date || (cashStream[i].isFixed() && discountCurve))
        {
            double amount = cashStream[i].getAmount(*this, volModelIR);

            if (!Maths::isZero(amount))
            {
                double discount = 0.0;

                if (discountCurve)
                {
                    discount = discountCurve->discountFactor(cashStream[i].getPayDate());
                }
                else
                {
                    discount = discountFactor(cashStream[i].getPayDate());
                }

                pvKnown += discount * amount;
            }
        }
    }

    /*
     * Now perform partial smoothing for the last segment (if necessary).  We
     * have the start rate (this will be a0) plus a constraint on the cash flow
     * list.
     *
     * We can make an arbitrary second constraint. Originally (Sept 1998) we
     * chose to set a2=0 and get linear forwards in the final segment.
     *
     * Now (Feb 2000) we have decided to set a2=-a1/2T, where T is the length
     * of the segment. This has the advantage that the forward rate will have
     * zero derivative, and will optimize the extreme value of the forward rate
     * in the last segment.
     *
     * i.e. The maximum forward rate using this scheme is better than the
     * maximum forward rate using linear forwards (for an increasing curve).
     *
     * There is actually a question of efficiency here - on many occasions we
     * will be immediately adding another cash flow, in which case this step is
     * a waste of time. Perhaps we need an extra parameter to denote this.
     */
    double pvRemainder = 0.0;

    if (discountCurve)
    {
        /*
         * With a discount curve, pvRemainder is the remainder to the base date.
         */
        pvRemainder = presentValue - pvKnown;
    }
    else
    {
        /*
         * Without a discount curve, the pvRemainder is the remainder at the
         * previous last date in the zero curve.
         *
         * However we need to take adjustments into account when calculating PV.
         */
        double adjustment = adjustmentFactor(baseDate, prevSegment.date);
        pvRemainder = (presentValue - pvKnown) / (prevSegment.getDiscount() * adjustment);
    }

    /*
     * The rate, discount and avgRate for the segment must be recalculated
     * to take into account that the previous segment has now been smoothed.
     */
    lastSegment.recalculate(*this, prevSegment, cs, startFwdRate, presentValue,
         pvKnown, pvRemainder, foundAdjAbs, discountCurve, volModelIR);

    return pvRemainder;
}


// ALIB: zcsmooth.c#3974 (GtoSZCRevert)
ZC3ZeroInterpolationConstSP ZC3ZeroCurve::revert()
{
    static const string method = "ZC3ZeroCurve::revert";

    if (data.size() < 2)
    {
        throw ModelException(method, "Data array size is less than 2");
    }

    if (!prevFull.get())
    {
        throw ModelException(method, "Previous cashstreams missing");
    }

    ZC3ZeroInterpolationConstSP prevInterpolation = data[data.size() - 1].interpolation;
    data.resize(data.size() - 1);   // removes last element
    cs = prevCs;
    prevCs = CashStreamSP(new CashStream());
    PV = prevPv;
    prevPv = 0.0;
    data.pop_back();
    data.push_back(*prevPoint);
    return prevInterpolation;
}


// ALIB: szcbuild.c#4111 (FindAbsoluteAdjustments)
bool ZC3ZeroCurve::findAbsoluteAdjustments(
    const ZC3CurveSegment& prevSegment, 
    const ZC3CurveSegment& segment) const
{
    static const string method = "ZC3ZeroCurve::findAbsoluteAdjustments";

    bool foundAdjAbsolute = false;

    for (int i = 0 ; i < adjAbsolute.size() ; i++ )
    {
        if (adjAbsolute[i].startsWithin(prevSegment.date, segment.date))
        {
            if (adjAbsolute[i].getEndDate() > segment.date)
            {
                string msg = Format::toString(
                    "Absolute adjustment %s  not contained within segment [%s,%s]",
                    adjAbsolute[i].toString().c_str(),
                    prevSegment.date.toString().c_str(),
                    segment.date.toString().c_str());
                throw ModelException(method, msg);
            }

            foundAdjAbsolute = true;
        }
    }

    return foundAdjAbsolute;
}


// ALIB: szcprvt.c#152 (GtoSZCModifyAdjustment)
void ZC3ZeroCurve::modifyAdjustment(const ZC3Adjustment& adjustment, double contRate)
{
    static const string method = "ZC3ZeroCurve::modifyAdjustment";

    for (int i = 0 ; i < adjustments.size() ; i++)
    {
        if (adjustments[i].isForSamePeriodAs(adjustment))
        {
            adjustments[i].modify(contRate);
            return;
        }
    }

    string msg = Format::toString("Could not find interval %s", adjustment.toString().c_str());
    throw ModelException(method, msg);
}


double ZC3ZeroCurve::adjustmentFactor(
    const DateTime& startDate, 
    const DateTime& endDate) const
{
    static const string method = "ZC3ZeroCurve::adjustmentFactor";

    // ALIB: zcsmooth.c#2168 (GtoSZCAdjustment)
    validate(false);

    /*
     * Compute the adjustments.
     */
    double adjustment = adjustmentFromList(adjustments, startDate.getDate(), endDate.getDate()); 

    if (!useSmoothing)
    {
        // ALIB: zcsmooth.c#2221 (GtoSZCAdjustmentUnsmoothed)
        /*
         * This routine is used for unsmoothed discount factors, and returns a
         * different answer if there are absolute adjustments in the structure.
         */

        // compute the re-adjustment
        if (!adjAbsolute.empty())
        {
            /*
             * We need to put into adjCopy the corresponding values from
             * adjustments. In general the dates in adjCopy should be a sub-set
             * of adjustments. Also in general there is only one date in adjCopy,
             * and it will be the first date in the adjustments array.
             *
             * Q: Why not store this permanently?
             * A: This function is usually only called during bootstrap
             *    and the values in these arrays will be continuously changing.
             */
            ZC3AdjustmentArray adjCopy(adjAbsolute);

            for (int i = 0 ; i < adjAbsolute.size() ; i++)
            {
                bool found = false;

                for (int j = 0 ; j < adjustments.size() ; j++)
                {
                    if (adjCopy[i].isForSamePeriodAs(adjustments[j]))
                    {
                        found = true;
                        adjCopy[i].modify(adjustments[j].getContRate());
                    }
                }

                if (!found)
                {
                    string msg = Format::toString(
                        "Absolute adjustment %s not found.", 
                        adjCopy[i].toString().c_str());
                    throw ModelException(method, msg);
                }
            }

            /*
             * Now we have two adjustment lists for the same dates. The ones in
             * adjustments are required, and the ones in adjCopy need to be
             * excluded (having already been included in the call to 
             * adjustment above).
             */
            double reAdjustment = adjustmentFromList(adjAbsolute, startDate.getDate(), endDate.getDate());
            adjustment *= reAdjustment;
            reAdjustment = adjustmentFromList(adjCopy, startDate.getDate(), endDate.getDate());
            adjustment /= reAdjustment;
        }
    }

    return adjustment;
}


// ALIB: zcsmooth.c#2350 (GtoSZCAdjustmentFromList)
double ZC3ZeroCurve::adjustmentFromList(
    const ZC3AdjustmentArray& adjustments, 
    const int                 startDate, 
    const int                 endDate) const
{
    static const string method = "ZC3ZeroCurve::adjustmentFromList";

    double rt = 0.0;  // adjustment

    /**
     * Compute the adjustments.
     */
    for (int i = 0 ; i < adjustments.size() ; i++)
    {
        int adjStartDate = adjustments[i].getStartDate().getDate();
        int adjEndDate = adjustments[i].getEndDate().getDate();

        /*
         * Here we are assuming that the adjustments are in ascending order, so
         * that if we have a date before the adjustment start date, that none
         * of the subsequent adjustments will apply.
         */
        if (endDate <= adjStartDate)
        {
            break;
        }

        int appStartDate = max(startDate, adjStartDate);  // apportion start date
        int appEndDate = min(endDate, adjEndDate);        // apportion start date

        if (appEndDate == adjEndDate && appStartDate == adjStartDate)
        {
            /*
             * This adjustment is fully within the requested interval, and
             * hence is applied completely.
             */
            rt += adjustments[i].getContRate();     // Full adjustment
        }
        else if (appStartDate >= appEndDate)
        {
            /*
             * This adjustment is fully before the requested interval, and
             * hence is not applied at all.
             */
        }
        else
        {
            /*
             * Need to apportion the adjustment. This requires us to compute
             * the time from startDate to date and the time from startDate to
             * endDate. The apportionment is the ratio in these times.
             *
             * We multiply the rt in the adjustment by this apportionment and
             * add it to the cumulative rt for this date.
             */
            double t1 = years(appStartDate, appEndDate);
            double t2 = years(adjStartDate, adjEndDate);

            /*
             * Because of the restrictions on the dates, although it is
             * possible for t2 to be zero (some day count conventions are
             * funny) then in this case t1 will also have to be zero.
             *
             * Hence the zero divide is essentially trapped by the first case
             * (t1 >= t2)
             */
            if (t1 >= t2)
            {
                rt += adjustments[i].getContRate();     // Full adjustment
            }
            else
            {
                rt += (t1 / t2) * adjustments[i].getContRate();
            }
        }
    }

    return Maths::isZero(rt) ? 1.0 : exp(-rt);
}


// ALIB: szcbuild.c#3104 (SZCFlatFwdAdjAbsolute)
void ZC3ZeroCurve::flatFwdAdjAbsolute(
	const ZC3CurveSegment& prevSegment,
	const ZC3CurveSegment& segment)
{
    static const string method = "ZC3ZeroCurve::flatFwdAdjAbsolute";

    /*
     * We need to match one or more absolute adjustments within the current
     * segment.
     *
     * To do this we will compute the relative adjustment implied by the
     * current choice of rate.
     *
     * This calculation is independent of whether we have a discount curve or
     * curve or not.
     */
    for (int i = 0 ; i < adjAbsolute.size() ; i++)
    {
        ZC3Adjustment& adjustment = adjAbsolute[i];

        if (adjustment.startsWithin(prevSegment.date, segment.date))
        {
            // lack of overlap was validated externally
            // use of this step validated externally

            // rate in this segment
            double avgRate = segment.avgRate;

            // time at adjustment start and end dates
            double ts = years(prevSegment.date, adjustment.getStartDate());
            double te = years(prevSegment.date, adjustment.getEndDate());

            adjustment.setContRate( (ts - te) * avgRate - log(adjustment.getDiscount()) );

            // this line changes the input curve
            modifyAdjustment(adjustment, adjustment.getContRate());
        }
    }
}


// ALIB: szcbuild.c#635 (GtoSZCAddRate)
void ZC3ZeroCurve::addRate(
    double                      rate, 
    const DateTime&             startDate, 
    const DateTime&             endDate, 
    const DayCountConvention&   dayCountConv,
    int                         basis,
    const ZC3ZeroInterpolation& interpolation)
{
    static const string method = "ZC3ZeroCurve::addRate";

    if (startDate == endDate)
    {
        // ignore this no-information case
        return;
    }

    if (startDate > endDate)
    {
        string msg = Format::toString(
            "Start date (%s) after end date (%s)",
            startDate.toString().c_str(),
            endDate.toString().c_str());
        throw ModelException(method, msg);
    }

    CashStream cs;

    /**
     * We need to generate a cash stream corresponding to the forward rate.
     * This has -1 * discount factor at the start date and 1 at the end date.
     */
    double discount = RateConversion::rateToDiscount(rate, startDate, endDate, &dayCountConv, basis);
    cs.addFixedFlow(startDate, -discount);
    cs.addFixedFlow(endDate, 1.0);

    // PV = 0.0 means price is in the cash stream
    bootstrap(cs, 0.0, interpolation, NULL, NULL);
}


// ALIB: szcbuild.c#1628 (SmoothZeroCurveAddSegment)
ZC3CurveSegment& ZC3ZeroCurve::addSegment(
    const DateTime&             date, 
    const ZC3ZeroInterpolation& interpolation, 
    ZC3CurveSegment*&           prevSegment)
{
    ZC3CurveSegment point(date, interpolation, *prevSegment);
    data.push_back(point);
    prevSegment = &(data[data.size() - 2]);  // vector may have resized
    return data.back();
}


void ZC3ZeroCurve::setUseSmoothing(bool state) const
{
    useSmoothing = state;
}


// delete [reset] fitted curve data, optionally copying data from other curve
void ZC3ZeroCurve::reset(const ZC3ZeroCurve* copy)
{
    if (copy)
    {
        // ALIB: zcsmooth.c#1712 (GtoSZCCopy)
        // base date, basis and DCC do not change
        data = copy->data;
        criticalDates = copy->criticalDates;
        adjustments = copy->adjustments;
        adjAbsolute = copy->adjAbsolute;

        cs = copy->cs;
        prevFull = copy->prevFull;
        prevCs = copy->prevCs;
    }
    else
    {
        data.resize(0);
        criticalDates.resize(0);
        adjustments.resize(0);
        adjAbsolute.resize(0);

        prevFull.reset();
        prevCs.reset(new CashStream());
        prevPoint.reset();
        prevPv = 0.0;
    }

    datesAndRates.resize(0);

    loBound = 0;
    hiBound = 0;
    noBootstrap = false;
}


void ZC3ZeroCurve::cleanup()
{
    cs.reset(new CashStream());
    prevCs.reset(new CashStream());
    PV = 0.0;
    prevPv = 0.0;
    prevFull.reset();
    prevPoint.reset();
}


class ZC3ZeroCurve::ZC3LogOfDiscFactorKey: public YieldCurve::IKey
{
public:
    ZC3LogOfDiscFactorKey(const ZC3ZeroCurve* pZc): 
        zc(*pZc), loRateTT(0), hiRateTT(0), 
            valueDate(pZc->valueDate.getDate()),
            baseDate(pZc->getBaseDate().getDate()),
            endDate(pZc->endDate().getDate())
        {
			// adjustments invalidate the fast interp algorithm
			hasFastInterp = zc.adjustments.empty() 
				         && zc.adjAbsolute.empty() 
						 && ZC3CurveSegment::hasFastInterp(zc.data);
        }

    virtual double calc(const DateTime& newLoDate,
                        const DateTime& newHiDate)
    {
        // see if we can quit early
        if (newLoDate.getDate() == newHiDate.getDate())
        {
            return 0.0;
        } 

        // cannot use fast interpolation
        if (!hasFastInterp)
        {
            // pv is an attribute on ZC3 and a method on ZeroCurve
            return log(static_cast<const ZeroCurve&>(zc).pv(newLoDate, newHiDate));
        }

		// extrapolation may vary according to interpolation type and extrapDate
        if (newHiDate.getDate() < baseDate || newHiDate.getDate() > endDate)
        {
            // pv is an attribute on ZC3 and a method on ZeroCurve
            return log(static_cast<const ZeroCurve&>(zc).pv(newLoDate, newHiDate));
        }


		// fast interpolation

        if (newHiDate == loDate)
        {
            hiDate = newHiDate;
            hiRateTT = loRateTT;
            loDate = newLoDate;
            double rate = zc.unadjustedZeroCouponRate(newLoDate);
            loRateTT = zc.years(valueDate, newLoDate.getDate()) * log(1.0 + rate);
        }
        else if (newLoDate == hiDate)
        {
            loDate = newLoDate;
            loRateTT = hiRateTT;
            hiDate = newHiDate;
            double rate = zc.unadjustedZeroCouponRate(newHiDate);
            hiRateTT = zc.years(valueDate, newHiDate.getDate()) * log(1.0 + rate);
        } 
        else 
        {
            if (newHiDate != hiDate)
            {
                hiDate = newHiDate;
                double rate = zc.unadjustedZeroCouponRate(newHiDate);
                hiRateTT = zc.years(valueDate, newHiDate.getDate()) * log(1.0 + rate);
            }

            if (newLoDate != loDate)
            {
                loDate = newLoDate;
                double rate = zc.unadjustedZeroCouponRate(newLoDate);
                loRateTT = zc.years(valueDate, newLoDate.getDate()) * log(1.0 + rate);
            }
        }

        return (loRateTT-hiRateTT); // = -1* (hiRateTT - loRateTT)
    }

private:

    const ZC3ZeroCurve& zc;
	bool                hasFastInterp;
    DateTime            loDate;     // last loDate
    DateTime            hiDate;     // last hiDate
    double              loRateTT;   // rate from base date to lo date * year frac
    double              hiRateTT;   // rate from base date to hi date * year frac
    const int           valueDate;  // cached value
    const int           baseDate;   // cached value
    const int           endDate;    // cached value
};


/**
 * Returns a key used to optimise repeated calculations of forward rates.
 */
YieldCurve::IKey* ZC3ZeroCurve::logOfDiscFactorKey() const
{
    return new ZC3LogOfDiscFactorKey(this);
}


/*
 * Reflection support.
 */
ZC3ZeroCurve::ZC3ZeroCurve()
  : ZeroCurve(TYPE),
    valueDateDiscount(1.0),
    baseDate(DateTime()),
    basis(CompoundBasis::ANNUAL),
    dcc(new Actual365F()),  // should restrict to Act/365F or Act/360
    criticalDates(0),
    datesAndRates(0),
    cs(new CashStream()),
    noBootstrap(false),
    useSmoothing(true),
    prevFull(NULL),
    prevCs(new CashStream()),
    prevPv(0.0),
    prevPoint(NULL),
    loBound(0), 
    hiBound(0)
{
    // ALIB: zcsmooth.c#1394
    // always have a point at base date
    addDiscountFactor(baseDate, 1.0, *ZC3ZeroInterpolation::make("Flat"));
}


class ZC3ZeroCurveHelper
{
public:
    /**
     * Invoked when Class is 'loaded' - we only really need the reflection
     * info to serialize the class.
     */
    static void load(CClassSP& clazz)
    {
        REGISTER(ZC3ZeroCurve, clazz);
        SUPERCLASS(ZeroCurve);
        EMPTY_SHELL_METHOD(defaultZC3ZeroCurve);
        FIELD_NO_DESC(valueDate);
        FIELD_NO_DESC(valueDateDiscount);
        FIELD_NO_DESC(baseDate);
        FIELD_NO_DESC(extrapDate);
        FIELD_NO_DESC(basis);
        FIELD_NO_DESC(dcc);
        FIELD_NO_DESC(data);
        FIELD_NO_DESC(criticalDates);
        FIELD_NO_DESC(datesAndRates);
        FIELD_NO_DESC(adjustments);
        FIELD_NO_DESC(adjAbsolute);
		FIELD_NO_DESC(noBootstrap);
        FIELD_MAKE_TRANSIENT(datesAndRates);
		FIELD_MAKE_TRANSIENT(noBootstrap);
    }

    static IObject* defaultZC3ZeroCurve()
    {
        return new ZC3ZeroCurve();
    }
};

CClassConstSP const ZC3ZeroCurve::TYPE = CClass::registerClassLoadMethod(
    "ZC3ZeroCurve", typeid(ZC3ZeroCurve), ZC3ZeroCurveHelper::load);

DEFINE_TEMPLATE_TYPE(ZC3ZeroCurveArray);


DRLIB_END_NAMESPACE
