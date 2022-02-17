//----------------------------------------------------------------------------
//
//   Group       : Quantitative Research & Development
//
//   Filename    : ZC3CurveInstrument.cpp
//
//   Description : Benchmark for bootstrapping curve.
//
//   Author      : Richard Appleton
//
//   Date        : 26th April 2005
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/ZC3CurveInstrument.hpp"
#include "edginc/ZC3ZeroCurve.hpp"
#include "edginc/ZC3Stub.hpp"
#include "edginc/StubPlacement.hpp"
#include "edginc/CashStream.hpp"
#include "edginc/RateConversion.hpp"
#include "edginc/CompoundBasis.hpp"
#include "edginc/DayCountConvention.hpp"
#include "edginc/DayCountConventionFactory.hpp"
#include "edginc/BadDayConvention.hpp"
#include "edginc/BadDayConventionFactory.hpp"
#include "edginc/IRVolBase.hpp"
#include "edginc/SwapTool.hpp"
#include "edginc/B30360.hpp"
#include "edginc/BenchmarkDate.hpp"
#include "edginc/StubSimple.hpp"
#include <algorithm>



DRLIB_BEGIN_NAMESPACE
using namespace std;

// force template instantiation for linker
static bool tmp1 = lessThan(ZC3RateDataSP(), ZC3RateDataSP());
static bool tmp2 = lessThan(ZC3SwapDataSP(), ZC3SwapDataSP());
static bool tmp3 = lessThan(ZC3TurnSP(), ZC3TurnSP());

StubSP ZC3SwapData::stubType(new StubSimple());


ZC3CurveInstrument::ZC3CurveInstrument(
    const DateTime& pStartDate,
    const DateTime& pEndDate,
    double          pRate,
    double          pBasisRate)
    : startDate(pStartDate), endDate(pEndDate), rate(pRate), basisRate(pBasisRate)
{
    if (startDate > endDate)
    {
        string msg = Format::toString(
            "Instrument start date (%s) cannot be after end date (%s)",
            startDate.toString().c_str(),
            endDate.toString().c_str());
        throw ModelException("ZC3CurveInstrument::ZC3CurveInstrument", msg);
    }
}


ZC3CurveInstrument::~ZC3CurveInstrument()
{
}


bool ZC3CurveInstrument::operator<(const ZC3CurveInstrument& other) const
{
    if (endDate == other.endDate)
    {
        return startDate < other.startDate;
    }
    else
    {
        return endDate < other.endDate;
    }
}


bool ZC3CurveInstrument::isStub(const DateTime& baseDate) const
{
    return endDate <= baseDate;
}


const DateTime& ZC3CurveInstrument::getStartDate() const
{
    return startDate;
}


const DateTime& ZC3CurveInstrument::getEndDate() const
{
    return endDate;
}


double ZC3CurveInstrument::getRate() const
{
    return rate;
}


double ZC3CurveInstrument::getBasisRate() const
{
    return basisRate;
}


/* static */
// ALIB: szcrates.c#1798
/*
 * ALIB uses 2 passes, once with rates and once with rates + basis.  For QLib
 * this has been changed to a single pass, hence the basis rate interpolation.
 */
ZC3RateData* ZC3RateData::interpolate(
    const ZC3RateData* lastBefore,
    const ZC3RateData& firstAfter,
    const DateTime&    startDate,
    const DateTime&    endDate)
{
    static const string method = "ZC3RateData::interpolate";

    double rate = 0.0;
    double rateWithBasis = 0.0;
    const DayCountConvention& dcc = *firstAfter.dcc;

    if (!lastBefore)
    {
        // extrapolate flat back to baseDate
        rate = firstAfter.rate;
        rateWithBasis = rate + firstAfter.basisRate;
    }
    else
    {
        // interpolate
        if (!dcc.equalTo(lastBefore->dcc))
        {
            string msg = Format::toString(
                "Not interpolating like with like [day count convention %s is not the same as before %s]",
                dcc.toString().c_str(),
                lastBefore->dcc->toString().c_str());
            throw ModelException(method, msg);
        }

        double t1 = dcc.years(lastBefore->endDate, endDate);
        double t2 = dcc.years(lastBefore->endDate, firstAfter.endDate);
        double r0 = lastBefore->rate;
        double r2 = firstAfter.rate;
        double b0 = lastBefore->basisRate;
        double b2 = firstAfter.basisRate;

        // new rate is r1 at t=t1 where we know r0 at t=0 and r2 at t=t2
        if (Maths::isZero(t2))
        {
            // possible with strange day count conventions
            // (split the difference in this case)
            rate = (r0 + r2) * 0.5;
            rateWithBasis = ((r0 + b0) + (r2 + b2)) * 0.5;
        }
        else
        {
            rate = r0 + (t1/t2) * (r2 - r0);
            rateWithBasis = (r0 + b0) + (t1/t2) * ((r2 + b2) - (r0 + b0));
        }
    }

    ZC3CashData* result = new ZC3CashData(startDate, endDate, rate, dcc, 0.0);
    result->basisRate = rateWithBasis - rate;
    return result;
}


ZC3RateData::ZC3RateData(
    const DateTime&           pStartDate,
    const DateTime&           pEndDate,
    double                    pRate,
    const DayCountConvention& pDcc,
    double                    pBasisRate)
    : ZC3CurveInstrument(pStartDate, pEndDate, pRate, pBasisRate), dcc(&pDcc)
{
}


const DayCountConvention& ZC3RateData::getDcc() const
{
    return *dcc;
}


bool ZC3RateData::isReallyBadDate(DateTime& badMMDate) const
{
    return false;
}


bool ZC3RateData::isCashStartingOn(const DateTime& date) const
{
    return false;
}


// ALIB: szcprvt.c#639 (GtoSZCPrepareStub)
void ZC3RateData::prepareStub(vector<ZC3Stub>& stubs, bool withBasis) const
{
    static const DayCountConvention* basisDcc = DayCountConventionFactory::make("Act/360");
    static const Holiday*            holidays = Holiday::noHolidays();
    static const BadDayConvention*   badDayConv = BadDayConventionFactory::make("N");

    double stubRate = getRate();
    if (withBasis)
    {
        stubRate += basisRate;
    }

    ZC3Stub stub(startDate, endDate, stubRate, *dcc, 0.0, *basisDcc, *badDayConv, *holidays);
    stubs.push_back(stub);
    sort(stubs.begin(), stubs.end());
}


ZC3CashData::ZC3CashData(
    const DateTime&             pStartDate,
    const DateTime&             pEndDate,
    double                      pRate,
    const DayCountConvention&   pDcc,
    double                      pBasisRate)
    : ZC3RateData(pStartDate, pEndDate, pRate, pDcc, pBasisRate)
{
}


double ZC3CashData::getContinuousRate(bool withBasis, bool withAdjustment) const
{
    // NB. check for stub is in calling function
	double r = rate;
	if (withBasis)
	{
		r += basisRate;
	}

    double df = RateConversion::rateToDiscount(r, startDate, endDate, dcc, CompoundBasis::SIMPLE);
    return -log(df);
}


bool ZC3CashData::futuresMarkToMarket(const ZC3ZeroCurve& szc, const IRVolBase& volModelIR)
{
    return false;
}


// ALIB:: szcrates.c#1608 (SZCRatesInitUsage)
void ZC3CashData::ratesInitUsage(
    const DateTime& baseDate,
    const DateTime* futuresMMDate,
    DateTime&       futuresStartDate,
    DateTime&       futuresEndDate,
    DateTime&       badMMDate,
    DateTime&       lastMMDate)
{
    /*
     * If this rate starts at the base date then it might end up being excluded.
     */
    if (startDate == baseDate)
    {
        if (lastMMDate.empty() || endDate > lastMMDate)
        {
            lastMMDate = endDate;
        }
    }
}


bool ZC3CashData::isCashStartingOn(const DateTime& date) const
{
    return startDate == date;
}


ZC3FuturesData::ZC3FuturesData(
    const DateTime&             pStartDate,
    const DateTime&             pEndDate,
    double                      pRate,
    const DayCountConvention&   pDcc,
    double                      pAdjustment,
    double                      pBasisRate)
    : ZC3RateData(pStartDate, pEndDate, pRate, pDcc, pBasisRate),
    adjustment(pAdjustment)
{
    // TBD!! zcinstyp.c#1045 re futures type rather than end date
}


double ZC3FuturesData::getContinuousRate(bool withBasis, bool withAdjustment) const
{
    double zcr = rate;
	if (withBasis)
	{
		zcr += basisRate;
	}

    // add futures adjustment
    if (withAdjustment)
    {
        zcr += adjustment;
    }

    // convert rate to discount
    double df = RateConversion::rateToDiscount(zcr, startDate, endDate, dcc, CompoundBasis::SIMPLE);
    return -log(df);
}


// ALIB: szcrates.c#2747
bool ZC3FuturesData::futuresMarkToMarket(const ZC3ZeroCurve& szc, const IRVolBase& volModelIR)
{
    string msg = Format::toString("TBD!! [%s line %d]", __FILE__, __LINE__);
    throw ModelException("ZC3FuturesData::futuresMarkToMarket", msg);

    /*
    B30360 dcc2;
    double futMat = dcc2.years(startDate, endDate);

    // szcrates.c#2782
    double rate = getRate();
    // double fwdRate = volModelIR.FutureToFwd(getRate(), szc.getBaseDate(), startDate, endDate, dcc, CompoundBasis::SIMPLE, szc);
    */

    return true;
}


// ALIB:: szcrates.c#1568 (SZCRatesInitUsage)
void ZC3FuturesData::ratesInitUsage(
    const DateTime& baseDate,
    const DateTime* futuresMMDate,
    DateTime&       futuresStartDate,
    DateTime&       futuresEndDate,
    DateTime&       badMMDate,
    DateTime&       lastMMDate)
{
    if (futuresStartDate.empty() || startDate < futuresStartDate)
    {
        futuresStartDate = startDate;
    }

    if (futuresEndDate.empty() || endDate > futuresEndDate)
    {
        futuresEndDate = endDate;
    }

    if (futuresMMDate && *futuresMMDate == endDate)
    {
        badMMDate = startDate;
    }
}


// ALIB:: szcrates.c#1668
bool ZC3FuturesData::isReallyBadDate(DateTime& badMMDate) const
{
    if (badMMDate == endDate)
    {
        badMMDate = startDate;
        return true;
    }
    else
    {
        return false;
    }
}


ZC3Turn::ZC3Turn(
    const DateTime&           pStartDate,
    const DateTime&           pEndDate,
    double                    pRate,
    const DayCountConvention& pDcc)
  : ZC3RateData(pStartDate, pEndDate, pRate, pDcc, 0.0)
{
    if (endDate <= startDate)
    {
        string msg = Format::toString(
            "Adjustment dates [%s,%s] are invalid",
            startDate.toString().c_str(),
            endDate.toString().c_str());
        throw ModelException("ZC3Turn::ZC3Turn", msg);
    }
}


bool ZC3Turn::futuresMarkToMarket(const ZC3ZeroCurve& szc, const IRVolBase& volModelIR)
{
    return false;
}


void ZC3Turn::ratesInitUsage(
    const DateTime& baseDate,
    const DateTime* futuresMMDate,
    DateTime&       futuresStartDate,
    DateTime&       futuresEndDate,
    DateTime&       badMMDate,
    DateTime&       lastMMDate)
{
}


// ALIB: szcprvt.c#101 (GtoAddAdjustment)
void ZC3Turn::addAdjustment(ZC3ZeroCurve& curve)
{
    ZC3Adjustment adjustment(startDate, endDate, getContinuousRate(false), 0.0);
    curve.addAdjustmentToList(adjustment);
    curve.addCriticalDate(startDate);
    curve.addCriticalDate(endDate);
}


ZC3AbsoluteTurn::ZC3AbsoluteTurn(
        const DateTime&           pStartDate,
        const DateTime&           pEndDate,
        double                    pRate,
        const DayCountConvention& pDcc)
  : ZC3Turn(pStartDate, pEndDate, pRate, pDcc)
{
}


// ALIB: szcprvt.c#226
double ZC3AbsoluteTurn::getContinuousRate(bool withBasis, bool withAdjustment) const
{
    // absolute adjustments ignore the rate
    return 0.0;
}


// ALIB: szcprvt.c#190 & szcrates.c#2996 (GtoAddAdjustment)
void ZC3AbsoluteTurn::addAdjustment(ZC3ZeroCurve& curve)
{
    /*
     * Need to add a zero adjustment (contRate=0) to both the absolute and
     * relative adjustments.
     *
     * Need to put the discount factor into the array of absolute adjustments.
     *
     * Absolute adjustments are usually provided as simple rates.
     */
    double discount = RateConversion::rateToDiscount(rate, startDate, endDate, dcc, CompoundBasis::SIMPLE);
    ZC3Adjustment adjustment(startDate, endDate, 0.0, discount);

    ZC3Turn::addAdjustment(curve);
    curve.addAdjustmentToList(adjustment);
}


// ALIB: [13.6.1] futstub.c#474
double ZC3AbsoluteTurn::turnEffectInPeriod(double* turnYF, int* turnDays) const
{
    if (turnYF)
    {
        double yearFrac = dcc->years(startDate, endDate);
        *turnYF -= yearFrac;
    }

    if (turnDays)
    {
        *turnDays += endDate.daysDiff(startDate);
    }

    return RateConversion::rateToDiscount(rate, startDate, endDate, dcc, CompoundBasis::SIMPLE);
}


ZC3RelativeTurn::ZC3RelativeTurn(
        const DateTime&           pStartDate,
        const DateTime&           pEndDate,
        double                    pRate,
        const DayCountConvention& pDcc)
  : ZC3Turn(pStartDate, pEndDate, pRate, pDcc)
{
}


// ALIB: szcrates.c#2974
double ZC3RelativeTurn::getContinuousRate(bool withBasis, bool withAdjustment) const
{
    // spread adjustments are assumed to be provided as continuously compounded
    double yearFrac = dcc->years(startDate, endDate);
    return rate * yearFrac;
}


// ALIB: szcprvt.c#101 (GtoAddAdjustment)
void ZC3RelativeTurn::addAdjustment(ZC3ZeroCurve& curve)
{
    ZC3Turn::addAdjustment(curve);
}


// ALIB: [13.6.1] futstub.c#469
double ZC3RelativeTurn::turnEffectInPeriod(double* turnYF, int* turnDays) const
{
    double yearFrac = dcc->years(startDate, endDate);
    return rate * yearFrac;
}


ZC3ImpliedTurn::ZC3ImpliedTurn(
        const DateTime&           pStartDate,
        const DateTime&           pEndDate,
        double                    pRate,
        const DayCountConvention& pDcc)
  : ZC3Turn(pStartDate, pEndDate, pRate, pDcc)
{
}


// ALIB: szcrates.c#3016
double ZC3ImpliedTurn::getContinuousRate(bool withBasis, bool withAdjustment) const
{
    // implied adjustments ignore the rate
    return 0.0;
}


// ALIB: szcrates.c#3116 (GtoAddAdjustment)
void ZC3ImpliedTurn::addAdjustment(ZC3ZeroCurve& curve)
{
    ZC3Turn::addAdjustment(curve);
}


// ALIB: [13.6.1] futstub.c#474
double ZC3ImpliedTurn::turnEffectInPeriod(double* turnYF, int* turnDays) const
{
    if (turnYF)
    {
        double yearFrac = dcc->years(startDate, endDate);
        *turnYF -= yearFrac;
    }

    if (turnDays)
    {
        *turnDays += endDate.daysDiff(startDate);
    }

    return RateConversion::rateToDiscount(rate, startDate, endDate, dcc, CompoundBasis::SIMPLE);
}


// ALIB: zcbuild3.c#3250
ZC3SwapData::ZC3SwapData(
        const DateTime&             pBaseDate,
        const DateTime&             pStartDate,
        const DateTime&             pEndDate,
        double                      pRate,
        const double*               pAdjustment,
        double                      pPrice,
        const MaturityPeriod&       pFixedIvl,
        const DayCountConvention&   pFixedDcc,
        const MaturityPeriod*       pFloatIvl,
        const DayCountConvention*   pFloatDcc,
        int                         pIncludeFlag,
        bool                        pValueFloating,
        bool                        pFloatRateFixed,
        const Holiday&              pHolidays,
        const BadDayConvention&     pBadDayConv,
        const ZeroCurve*            pFixCurve,
        const MaturityPeriod*       pBasisIvl,
        const DayCountConvention*   pBasisDcc,
        double                      pBasisRate)
    : ZC3CurveInstrument(pStartDate, pEndDate, pRate, pBasisRate),
    benchmark(pIncludeFlag),
    adjustment(pAdjustment ? *pAdjustment : 0.0),
    price(Maths::isZero(pPrice) ? 1.0 : pPrice),
    fixedIvl(pFixedIvl),
    fixedDcc(pFixedDcc),
    valueFloating(pValueFloating),
    floatIvl(pFloatIvl ? *pFloatIvl : pFixedIvl),
    floatDcc(pFloatDcc ? *pFloatDcc : pFixedDcc),
    floatFixed(pFloatRateFixed && (pStartDate == pBaseDate)),
    floatFixRate(0.0),
    basisIvl(pBasisIvl ? *pBasisIvl : (pFloatIvl ? *pFloatIvl : pFixedIvl)),
    basisDcc(pBasisDcc ? *pBasisDcc : (pFloatDcc ? *pFloatDcc : pFixedDcc)),
	recursiveInFlatSection(false)
{
    static const string method = "ZC3SwapData::ZC3SwapData";

    if (benchmark < 1 || benchmark > 2)
    {
        string msg = Format::toString(
            "Invalid benchmark for swap: %s",
            endDate.toString().c_str());
        throw ModelException(method, msg);
    }

    if (benchmark == 2 && pAdjustment == NULL)
    {
        string msg = Format::toString(
            "Must provide adjustments for benchmark flag 2: %s",
            endDate.toString().c_str());
        throw ModelException(method, msg);
    }

    if (floatFixed)
    {
        /*
         * What rate do we need? We need the rate for the
         * floating interval after the start date. We then
         * use this date to interpolate from the fixing curve.
         *
         * In practice there will be often only be one rate in
         * the fixing curve, and this will be the correct date
         * as well!
         */
        DateTime rateEndDate = DateFwdThenAdjust(pBaseDate, floatIvl, 1, pBadDayConv, pHolidays);
        floatFixRate = pFixCurve->fwd(pBaseDate, rateEndDate, &floatDcc, CompoundBasis::SIMPLE);
    }
}


bool ZC3SwapData::isExcluded(const DateTime& lastZeroDate, const DateTime& baseDate) const
{
    return endDate <= lastZeroDate && endDate > baseDate;
}


bool ZC3SwapData::isForPass(int pass) const
{
    return benchmark == pass;
}


bool ZC3SwapData::isEndStub(const StubPlacement& stubPos, const MaturityPeriod& tenor) const
{
    return stubPos.isEndStub(tenor, startDate, endDate);
}


double ZC3SwapData::getContinuousRate(bool withBasis, bool withAdjustment) const
{
    return 0.0;
}


bool ZC3SwapData::isRecursiveIncludedInFlatSection() const
{
    return recursiveInFlatSection;
}


void ZC3SwapData::setRecursiveIncludedInFlatSection()
{
    recursiveInFlatSection = true;
}


void ZC3SwapData::setBenchmark(int pBenchmark)
{
    benchmark = pBenchmark;
}


// ALIB: szcprvt.c#639 (GtoSZCPrepareStub)
void ZC3SwapData::prepareStub(
    vector<ZC3Stub>&        stubs,
    double                  couponRate,
    const BadDayConvention& badDayConv,
    const Holiday&          holidays,
    const bool              withBasis) const
{
    static const DayCountConvention* basisDcc = DayCountConventionFactory::make("Act/360");

    ZC3Stub stub(
        startDate,
        endDate,
        couponRate,
        fixedDcc,
        withBasis ? basisRate : 0.0,
        withBasis ? *basisDcc : *ZC3Stub::STUB_DCC,
        badDayConv,
        holidays);

    stubs.push_back(stub);
    sort(stubs.begin(), stubs.end());
}


// ALIB: szcswaps.c#1025 & 2229 [GtoSZCAddOneSwapWithBasis]
CashFlowArraySP ZC3SwapData::getKnownCashflows(
    const ZC3ZeroCurve&     szc,
    double                  couponRate,
    const StubPlacement&    stubPos,
    const BadDayConvention& badDayConv,
    const Holiday&          holidays,
    const Holiday&          basisHolidays,
    const bool              withBasis,
    const bool              estimating) const
{
    static const string method = "ZC3SwapData::getKnownCashflows";

    /*
     * Fixed leg of the swap.
     */
    bool endStub = isEndStub(stubPos, fixedIvl);
    bool includeNotionalFlows = estimating ? false : (withBasis ? true : !valueFloating);

    CashFlowArraySP cfl = makeCFL(couponRate, endStub, badDayConv,
        holidays, fixedIvl, fixedDcc, includeNotionalFlows);

    if (!estimating && withBasis)
    {
        /*
         * Since we subtracted initial notional we amend the first cash
         * flow to be equal to the price.
         */
        (*cfl)[0].amount += (1.0 - price);

        // basis fixed payments to be added to the regular fixed payments
        cfl = CashFlow::merge(cfl, makeCFL(basisRate, endStub, badDayConv,
            basisHolidays, basisIvl, basisDcc, false));
    }
    else if (!valueFloating)
    {
        /*
         * We only use the price when NOT valueFloating. We used
         * GTO_SUBTRACT_INITIAL to get a cash flow of -1 at the start, so
         * we need to amend this to be equal to the price.
         */
        if (floatFixed)
        {
            if (!Maths::equals(1.0,price))
            {
                string msg = "Prices must equal 1.0 if a floating rate has "
                             "been fixed and valueFloating is false";
                throw ModelException(method, msg);
            }

            DateTime rateEndDate = DateFwdThenAdjust
                (szc.getBaseDate(), floatIvl, 1, badDayConv, holidays);
            DateTime lastZeroDate = szc.endDate();

            if (rateEndDate >= lastZeroDate)
            {
                string msg = Format::toString(
                    "Cannot calculate floating leg price for swap with first floating"
                    " rate fixed.  Coupon date (%s) is after last point in curve (%s)",
                    rateEndDate.toString().c_str(),
                    lastZeroDate.toString().c_str());
                throw ModelException(method, msg);
            }

            double yearFrac = floatDcc.years(szc.getBaseDate(), rateEndDate);
            double discount = szc.discountFactor(rateEndDate);
            (*cfl)[0].amount += 1.0 - (1.0 + floatFixRate * yearFrac) * discount;
        }
        else
        {
            (*cfl)[0].amount += (1.0 - price);
        }
    }

    return cfl;
}


// ALIB: szcswaps.c#1113
void ZC3SwapData::getFloatingCashflows(
    CashStream&             cs,
    double                  couponRate,
    const StubPlacement&    stubPos,
    const BadDayConvention& badDayConv,
    const Holiday&          holidays,
    const bool              mustValueFloating) const
{
    static const string method = "ZC3SwapData::getFloatingCashflows";

    /*
     * If necessary add the floating leg of the swap.
     */
    if (mustValueFloating || valueFloating)
    {
        bool isFloatingEndStub = isEndStub(stubPos, floatIvl);

        cs.addFloatLegVanilla(
            -1.0,     // notional
            0.0,      // spread
            true,     // spread is additive
            startDate,
            startDate,
            startDate,
            floatIvl,
            endDate,
            *stubType,
            isFloatingEndStub,
            false,    // no principal payments
            false,    // no principal payments
            badDayConv,
            badDayConv,
            badDayConv,
            floatFixed,
            floatFixRate,
            floatDcc,
            holidays);
    }
}


// ALIB: cashflow.c#949 (GtoMakeCFL)
CashFlowArraySP ZC3SwapData::makeCFL(
    double                    couponRate,
    bool                      endStub,
    const BadDayConvention&   badDayConv,
    const Holiday&            holidays,
    const MaturityPeriod&     interval,
    const DayCountConvention& dcc,
    bool                      includeNotionalFlows) const
{
    static const string method = "ZC3SwapData::makeCFL";

    CashFlowArray cashflows = SwapTool::cashflows(startDate, endDate, *stubType,
        endStub, true, badDayConv, badDayConv, holidays, true,
        includeNotionalFlows, includeNotionalFlows, couponRate, interval, dcc);

    return CashFlowArraySP(new CashFlowArray(cashflows));
}


// ALIB: szcswaps.c#767
double ZC3SwapData::getCouponRate(
    const ZeroCurve&        discountCurve,
    const ZeroCurve*        estimatingCurve,
    const bool              convDelayAdj,
    const IRVolBase*        volModelIR,
    const StubPlacement&    stubPos,
    const BadDayConvention& badDayConv,
    const Holiday&          holidays,
    const bool              withBasis) const
{
    double couponRate = SwapTool::swapRate(discountCurve,
        startDate,
        endDate,
        fixedIvl,
        fixedDcc,
        withBasis ? true : valueFloating,
        withBasis ? 0.0 : price,
        estimatingCurve,
        &floatIvl,
        &floatDcc,
        floatFixed,
        floatFixRate,
        convDelayAdj,
        volModelIR,
        stubPos,
        badDayConv,
        badDayConv,
        badDayConv,
        holidays);

    couponRate += adjustment;
    return couponRate;
}


DateTimeArraySP ZC3SwapData::getFixedDates(bool stubAtEnd) const
{
    DateTimeArraySP dates(SwapTool::dateArray(startDate, startDate, endDate, fixedIvl, stubAtEnd));
    return dates;
}


/* static */
// ALIB: szcswaps.c#3143 (GetInterpRates)
ZC3SwapDataArraySP ZC3SwapData::getInterpRates(
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
        const bool              withBasis)
{
    static const string method = "ZC3SwapData::getInterpRates";

    /*
     * No swaps to add is not a failure.
     */
    if (data.size() == 0)
    {
        return ZC3SwapDataArraySP(new ZC3SwapDataArray(data));
    }

    bool indexSwap = false;

    /*
     * Generate the coupon dates for the last swap.
     */
    DateTime startDate = data.back()->getStartDate();
    DateTime lastSwapDate = data.back()->getEndDate();
    DateTime lastCurveDate = estimatingCurve ? estimatingCurve->endDate() : discountCurve.endDate();

    if (lastSwapDate <= lastCurveDate)
    {
        throw ModelException(method, "last swap date does not extend the zero curve");
    }

    /*
     * Generate a date list with a stub at the beginning.
     *
     * This means that if startDate is not on-cycle with lastSwapDate,
     * then startDate will not be in the resulting date list, but
     * we will get an extra date earlier than lastCurveDate instead.
     */
    DateTimeArraySP dl = data.back()->getFixedDates(false);

    DateTime firstSwapDate = removeCouponsNotNeeded(*dl, lastCurveDate, badDayConv, holidays);
    validateSwapDatesToInterp(data, withBasis, indexSwap);
    double fixedFreq = data.front()->fixedIvl.approxAnnualFrequency();
    double basisFreq = withBasis ? data.front()->basisIvl.approxAnnualFrequency() : 0.0;

    /*
     * The useIndex array is a list of the swaps that are actually going to be
     * used in the interpolation.
     */
    ZC3SwapDataArray useIndex;
    for (size_t j = 0 ; j < data.size() ; j++)
    {
        DateTime maturityDate = data[j]->getEndDate();

        /*
         * Maturity dates which do not extend the curve are irrelevant and
         * hence ignored.
         */
        if (maturityDate <= lastCurveDate)
        {
            continue;
        }

        useIndex.push_back(data[j]);

        /*
         * Now add the actual swap date to the coupon date list (which is the
         * on-cycle list of maturity dates).
         */
        int found = count(dl->begin(), dl->end(), maturityDate);
        if (found == 0)
        {
            dl->push_back(maturityDate);
        }
    }

    sort(dl->begin(), dl->end());

    /*
     * We exclude all dates in the combined date list which are less than or
     * equal to the firstSwapDate calculated earlier.
     *
     * Note that the date list has no duplicates and is sorted in date order.
     */
    DateTimeArray::iterator iterator = dl->begin();
    while (iterator != dl->end())
    {
        if (*iterator > firstSwapDate)
        {
            break;
        }

        iterator++;
    }

    dl->erase(dl->begin(), iterator);

    if (dl->size() == 0)
    {
        throw ModelException(method, "None of the swaps really extends the zero curve");
    }

    /*
     * OK - now we get to work.
     * We get the swap rate for the last date covered by the curve.
     */
    bool valueFloating = withBasis ? true : data.front()->valueFloating;
    const ZeroCurve* tc1 = &discountCurve;
    const ZeroCurve* tc2 = valueFloating ? (estimatingCurve ? estimatingCurve : tc1) : NULL;

    /*
     * We need some validation to show that we have at least added some
     * money market rates into the curve. In the old n+i models, this was
     * typically explicit - something like needing the 1Y money market rate.
     *
     * We are not going to be this explicit here.
     */
    if (firstSwapDate <= startDate)
    {
        throw ModelException(method, "No rates in stub curve after swap start date "
                                     " - cannot use coupon interpolation method");
    }

    DateTime prevMaturityDate = firstSwapDate;
    double prevCouponRate = SwapTool::swapRate(*tc1,
        startDate,
        prevMaturityDate,
        data.front()->fixedIvl,
        data.front()->fixedDcc,
        valueFloating,
        data.front()->price,
        tc2,
        &data.front()->floatIvl,
        &data.front()->floatDcc,
        data.front()->floatFixed,
        data.front()->floatFixRate,
        convDelayAdj,
        volModelIR,
        stubPos,
        badDayConv,
        badDayConv,
        badDayConv,
        holidays);

    double prevBasisRate = 0.0;
    if (withBasis)
    {
        /*
         * Calculate basis swap spread corresponding to the prevCouponRate.
         *
         * This uses the formula:
         *    price = fixed leg + basis leg
         */
        double fixedLeg = SwapTool::swapFixedPV(*tc1,
                            prevCouponRate,
                            startDate,
                            data.front()->fixedIvl,
                            prevMaturityDate,
                            data.front()->fixedDcc,
                            *stubType,
                            stubPos,
                            false,      // subtract initial notional payment
                            true,       // add final notional payment
                            badDayConv,
                            badDayConv,
                            holidays,
                            startDate);

        double basisAnnuity = SwapTool::swapFixedPV(*tc1,
                            1.0,
                            startDate,
                            data.front()->basisIvl,
                            prevMaturityDate,
                            data.front()->basisDcc,
                            *stubType,
                            stubPos,
                            false,      // subtract initial notional payment
                            false,      // add final notional payment
                            badDayConv,
                            badDayConv,
                            basisHolidays,
                            startDate);

        prevBasisRate = (data.front()->price - fixedLeg) / basisAnnuity;
    }

    if (annualize)
    {
        /*
         * It is not totally clear what annualization means when currency basis
         * / cost of funds is involved. Probably this is just a stupid thing to
         * do. However the usual method is 4+i which does not involve annualization.
         */
        prevCouponRate = pow(1.0 + prevCouponRate/fixedFreq, fixedFreq) - 1.0;
        prevBasisRate  = pow(1.0 + prevBasisRate/basisFreq, basisFreq) - 1.0;
    }

    /*
     * Now perform the interpolation. At each point we will use linear
     * interpolation with the 30/360 day count convention to compute time.
     *
     * This interpolation is consistent with 3+i and 4+i methodology.
     *
     * We will interpolate on adjusted maturity dates for backward compatibility.
     */
    ZC3SwapDataArraySP swapDataOut(new ZC3SwapDataArray());
    B30360 dcc;
    int jdx = -1;
    for (size_t k = 0 ; k < useIndex.size() ; k++)
    {
        DateTime maturityDate = useIndex[k]->getEndDate();
        double nextCouponRate = useIndex[k]->getRate();

        if (maturityDate <= lastCurveDate)
        {
            string msg = Format::toString(
                "Maturity date (%s) <= last curve date (%s)",
                maturityDate.toString().c_str(),
                lastCurveDate.toString().c_str());
            throw ModelException(method, msg);
        }

        // Calculate time for this segment (using adjusted dates)
        DateTime prevMaturityDateAdj = badDayConv.adjust(prevMaturityDate, &holidays);
        DateTime maturityDateAdj = badDayConv.adjust(maturityDate, &holidays);
        double segTime = dcc.years(prevMaturityDateAdj, maturityDateAdj);

        if (annualize)
        {
            nextCouponRate = pow(1.0 + nextCouponRate / fixedFreq, fixedFreq) - 1.0;
        }

        double nextBasisRate = 0.0;
        if (withBasis)
        {
            nextBasisRate = useIndex[k]->basisRate;
            if(annualize)
            {
                nextBasisRate = pow(1.0 + nextBasisRate / basisFreq, basisFreq) - 1.0;
            }
        }

        do
        {
            double rate;
            double basisRate;

            ++jdx;
            if (jdx >= dl->size())
            {
                string msg = Format::toString("Program bug [file %s, line %d]", __FILE__, __LINE__);
                throw ModelException(method, msg);
            }

            DateTime adjDate = badDayConv.adjust((*dl)[jdx], &holidays);
            double time = dcc.years(prevMaturityDateAdj, adjDate);

            /*
             * By testing for time > 0 we avoid the possible case where
             * segTime = 0 (very unlikely) since segTime >= time in all cases.
             */
            if (Maths::isPositive(time))
            {
                rate = prevCouponRate + (time / segTime) * (nextCouponRate - prevCouponRate);
                basisRate = prevBasisRate + (time / segTime) * (nextBasisRate - prevBasisRate);
            }
            else
            {
                rate = prevCouponRate;
                basisRate = prevBasisRate;
            }

            if (annualize)
            {
                rate = fixedFreq * (pow(1.0 +rate, 1.0/fixedFreq) - 1.0);
                basisRate = basisFreq * (pow(1.0 + basisRate, 1.0/basisFreq) - 1.0);
            }

            bool valueFloating = true;
            if (!indexSwap)
            {
                valueFloating = withBasis ? false : useIndex[k]->valueFloating;
            }

            // create output SwapData and add to output array
            ZC3SwapData* swap = new ZC3SwapData(
                tc1->getBaseDate(),
                startDate,
                (*dl)[jdx],
                rate,
                NULL,
                useIndex[k]->price,
                useIndex[k]->fixedIvl,
                useIndex[k]->fixedDcc,
                &useIndex[k]->floatIvl,
                &useIndex[k]->floatDcc,
                useIndex[k]->benchmark,
                valueFloating,
                data[0]->floatFixed,
                holidays,
                badDayConv,
                NULL,
                &useIndex[k]->basisIvl,
                &useIndex[k]->basisDcc,
                basisRate);
            swap->floatFixed = useIndex[k]->floatFixed;
            swap->floatFixRate = useIndex[k]->floatFixRate;
            swapDataOut->push_back(ZC3SwapDataSP(swap));
        }
        while(jdx < dl->size() && useIndex[k]->getEndDate() > (*dl)[jdx]);

        /*
         * Did we find the date? We must have done since all the dates used in
         * maturityDates where put into dl earlier in the function.
         */
        if (useIndex[k]->getEndDate() != (*dl)[jdx])
        {
            string msg = Format::toString("Program bug [file %s at line %d]", __FILE__, __LINE__);
            throw ModelException(method, msg);
        }

        prevMaturityDate = maturityDate;
        prevCouponRate = nextCouponRate;
        prevBasisRate = nextBasisRate;
    }

    return swapDataOut;
}


/* static */
DateTime ZC3SwapData::removeCouponsNotNeeded(
    const DateTimeArray&    dl,
    const DateTime&         lastCurveDate,
    const BadDayConvention& badDayConv,
    const Holiday&          holidays)
{
    DateTime firstSwapDate = dl.front();
    DateTime nextSwapDate = firstSwapDate;
    DateTime nextSwapDateUnadj = firstSwapDate;

    // Remove coupons that are not needed
    for (int i = 0 ; nextSwapDate <= lastCurveDate ; i++)
    {
        /*
         *
         * Consider 3 use cases - the last is atypical but must be considered.
         *
         * 1. MM curve built with 4 futures. Last curve date is 1.15 years
         *    after the base date. The first swap date is 1 year after the
         *    base date. Need to use the 1Y swap rate to interpolate the 1.5Y
         *    swap rate.
         *
         * 2. MM curve has 1 year point as a holiday at the end of the month,
         *    and hence by modified-following rules is just short of 1Y from
         *    base date. Hence the first swap date could be 0 years or 0.5
         *    years (annual / semi-annual) after base date. In this case
         *    we want to calculate the 1Y swap rate rather than the 0 or 0.5
         *    year swap rates. How can we force this?
         *
         * 3. Stub curve is out to 4 years, and the 4-year point is a holiday.
         *    User may be trying to demonstrate that the 4Y curve + 5Y rate
         *    is the same as the 5Y curve built directly.
         */

        /*
         * What we do is consider the next swap date and see whether this
         * is greater than the lastCurveDate purely by virtue of bad day
         * conventions. It is possible for the badDayConv to move this date
         * backwards so that the nextSwapDate becomes <= lastCurveDate.
         */
        firstSwapDate = nextSwapDateUnadj;
        nextSwapDateUnadj = dl[i];
        nextSwapDate = badDayConv.adjust(nextSwapDateUnadj, &holidays);
    }

    return firstSwapDate;
}


/* static */
// ALIB: szcswaps.c#4863
void ZC3SwapData::validateSwapDatesToInterp(
    const ZC3SwapDataArray& data,
    const bool              withBasis,
    const bool              indexSwap)
{
    static const string method = "ZC3SwapData::validateSwapDatesToInterp";

    /*
     * Validates that coupon interpolation is valid.
     *
     * The following are errors:
     * 1. Inconsistent swap definition.
     * 2. Use of adjustments (benchmark flag > 1).
     * 3. Start dates must all be the same.
     * 4. Maturity dates must be increasing.
     * 5. Inconsistent use of valueFloating.
     * 6. Price not par.
     *
     * It does not matter whether the original swaps are on-cycle or not.
     * We also populate the genSwapDates, genSwapRates, genBasisRates
     * at this point.
     */
    bool valueFloating = true;
    if (!indexSwap)
    {
        // in this context
        valueFloating = withBasis ? false : data.front()->valueFloating;
    }

    for (size_t i = 0 ; i < data.size() ; i++)
    {
        // 1. value floating test
        bool thisValueFloating = true;
        if (!indexSwap)
        {
            // in this context
            thisValueFloating = withBasis ? false : data[i]->valueFloating;
        }

        if (thisValueFloating != valueFloating)
        {
            throw ModelException(method, "Value floating must be all TRUE or FALSE");
        }

        // 2. inconsistent swap definition test
        static const char* msg1 = "benchmark [%d] : coupon interpolation requires consistent coupon intervals";
        static const char* msg2 = "benchmark [%d] : coupon interpolation requires consistent day count conventions";
        static const char* msg3 = "benchmark [%d] : coupon interpolation requires consistent behaviour of first rate fixing";

        if (data.front()->fixedIvl.approxAnnualFrequency()
         != data[i]->fixedIvl.approxAnnualFrequency())
        {
            string msg = Format::toString(msg1, i);
            throw ModelException(method, msg);
        }

        if (!data.front()->fixedDcc.equalTo(&data[i]->fixedDcc))
        {
            string msg = Format::toString(msg2, i);
            throw ModelException(method, msg);
        }

        if (valueFloating || withBasis)
        {
            if (data.front()->floatIvl.approxAnnualFrequency()
             != data[i]->floatIvl.approxAnnualFrequency())
            {
                string msg = Format::toString(msg1, i);
                throw ModelException(method, msg);
            }

            if (!data.front()->floatDcc.equalTo(&data[i]->floatDcc))
            {
                string msg = Format::toString(msg2, i);
                throw ModelException(method, msg);
            }

            if ((data.front()->floatFixed != data[i]->floatFixed)
             || (data.front()->floatFixed && data.front()->floatFixRate != data[i]->floatFixRate))
            {
                string msg = Format::toString(msg3, i);
                throw ModelException(method, msg);
            }

            if (withBasis)
            {
                if (data.front()->basisIvl.approxAnnualFrequency()
                 != data[i]->basisIvl.approxAnnualFrequency())
                {
                    string msg = Format::toString(msg1, i);
                    throw ModelException(method, msg);
                }

                if (!data.front()->basisDcc.equalTo(&data[i]->basisDcc))
                {
                    string msg = Format::toString(msg2, i);
                    throw ModelException(method, msg);
                }
            }
        }

        // 3. use of adjustments test
        if (data[i]->isForPass(2))
        {
            throw ModelException(method, "Coupon interpolation and swap rate adjustments cannot be mixed");
        }

        // 4. start date test
        if (data.front()->startDate != data[i]->startDate)
        {
            string msg = Format::toString("benchmark [%d] : swap start dates must be the same for coupon interpolation");
            throw ModelException(method, msg);
        }

        // 5. maturity date test
        // (NB. start date > end date for each swap validated on parsing data)
        if (i > 0 && (data[i]->endDate < data[i-1]->endDate))
        {
            string msg = Format::toString(
                "benchmark [%d] : maturity date %s is not in increasing order",
                i, data[i]->endDate.toString().c_str());
            throw ModelException(method, msg);
        }

        // 6. prices test
        if (!valueFloating && !Maths::equals(data[i]->price, 1.0))
        {
            string msg = Format::toString("benchmark [%d] : cannot interpolate coupon rates unless prices is 1");
            throw ModelException(method, msg);
        }
    }
}


DRLIB_END_NAMESPACE
