//----------------------------------------------------------------------------
//
//   Group       : QR Core Analytics Team
//
//   Filename    : ZCBrzFI.cpp
//
//   Description : Curve factory implementation ported from ALIB BRZ_ZC_FI.
//
//   Author      : Richard Appleton
//
//   Date        : 24th March 2006
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/ZCBrzFI.hpp"
#include "edginc/Atomic.hpp"

#include "edginc/Addin.hpp"
#include "edginc/Actual360.hpp"
#include "edginc/Business252.hpp"
#include "edginc/DayCountConventionFactory.hpp"
#include "edginc/BadDayConventionFactory.hpp"
#include "edginc/BadDayFollowing.hpp"
#include "edginc/BadDayNone.hpp"
#include "edginc/RateConversion.hpp"
#include "edginc/CompoundBasis.hpp"
#include "edginc/Hashtable.hpp" // for hash_string
#include "edginc/LMParabCurve.hpp"
#include <algorithm>


DRLIB_BEGIN_NAMESPACE

double ZCBrzFI::BRZ_FUTURES_PRICE_UNIT = 100000.0;


ZCBrzFI::ZCBrzFI()
  : IZeroCurveFactory(TYPE)
{
}


ZCBrzFI::~ZCBrzFI()
{
    //empty
}


int ZCBrzFI::hashCode() const
{
    int hCode = 0;
    return hCode;
}


bool ZCBrzFI::equals(const IZeroCurveFactory* zcm) const
{
    //check type
    if (!ZCBrzFI::TYPE->isInstance(zcm)) return false;

    const ZCBrzFI* zcb = dynamic_cast<const ZCBrzFI*>(zcm);

    //check instance variables
    if (copomDates.get())
    {
        if (!zcb->copomDates.get() || !copomDates->equalTo(zcb->copomDates.get()))
        {
            return false;
        }
    }
    else if (zcb->copomDates.get())
    {
        return false;
    }

    return true;
}


bool ZCBrzFI::isFloatLegValued() const
{
    return false;
}


void ZCBrzFI::validatePop2Object()
{
}


void ZCBrzFI::validate(
    const DateTime&                today,
    const DateTime&                spotDate,
    const ZeroCurveBenchmarkArray& benchmarks,
    const DayCountConvention*      moneyMarketDcc,
    const DayCountConvention*      swapFixedDcc,
    const DayCountConvention*      swapFloatDcc,
    const DayCountConvention*      swapBasisDcc,
    const MaturityPeriod*          swapFixedIvl,
    const MaturityPeriod*          swapFloatIvl,
    const MaturityPeriod*          swapBasisIvl,
    const BadDayConvention*        swapBadDayConvention,
    const bool                     floatRateFixed,
    const ExpiryArray*             fixDates,
    const DoubleArray*             fixRates,
    IRVolBaseConstSP               volModelIR,
    HolidayConstSP                 holidays,
    const CurrencyBasisWrapper&    ccyBasis,
    const Expiry*                  futuresMaturity
    ) const
{
    static const string method = "ZCBrzFI::validate";

    // verify mandatory parameters

    //if (moneyMarketDcc)
    //{
    //    throw ModelException(method, "moneyMarketDcc specified");
    //}

    //if (swapFixedDcc)
    //{
    //    throw ModelException(method, "FixedDcc specified");
    //}

    if (swapFloatDcc)
    {
        throw ModelException(method, "FloatDcc specified");
    }

    if (swapBasisDcc)
    {
        throw ModelException(method, "BasisDcc specified");
    }

    //if (swapFixedIvl)
    //{
    //    throw ModelException(method, "FixedIvl specified");
    //}

    if (swapFloatIvl)
    {
        throw ModelException(method, "FloatIvl specified");
    }

    if (swapBasisIvl)
    {
        throw ModelException(method, "BasisIvl specified");
    }

    if (swapBadDayConvention)
    {
        throw ModelException(method, "BadDayConvention specified");
    }

    if (fixDates)
    {
        throw ModelException(method, "FixDates specified");
    }

    if (fixRates)
    {
        throw ModelException(method, "FixRates specified");
    }

    if (volModelIR.get())
    {
        throw ModelException(method, "IrVol specified");
    }

    //if (futuresMaturity)
    //{
    //    throw ModelException(method, "FuturesMaturity specified");
    //}

    if (!(benchmarks[0]->isCash() && benchmarks[0]->getEnd()->toString() == "ON"))
    {
        throw ModelException(method, "First instrument must be money market ON rate");
    }

    if (!(benchmarks[1]->isCash() && benchmarks[1]->getEnd()->toString() == "1D"))
    {
        throw ModelException(method, "Second instrument must be money market 1D rate");
    }
}


ZeroPairSP ZCBrzFI::bootstrap(
        const DateTime&                tradeDate,
        const DateTime&                spotDate,
        const ZeroCurveBenchmarkArray& benchmarks,
        const DayCountConvention*      moneyMarketDcc,
        const DayCountConvention*      swapFixedDcc,
        const DayCountConvention*      swapFloatDcc,
        const DayCountConvention*      swapBasisDcc,
        const MaturityPeriod*          swapFixedIvl,
        const MaturityPeriod*          swapFloatIvl,
        const MaturityPeriod*          swapBasisIvl,
        const BadDayConvention*        swapBadDayConvention,
        const bool                     floatRateFixed,
        const ExpiryArray*             fixDates,
        const DoubleArray*             fixRates,
        IRVolBaseConstSP               volModelIR,
        HolidayConstSP                 holidays,
        const CurrencyBasisWrapper&    ccyBasis,
        const Expiry*                  futuresMaturity
    ) const
{
    validate(
        tradeDate,
        spotDate,
        benchmarks,
        moneyMarketDcc,
        swapFixedDcc,
        swapFloatDcc,
        swapBasisDcc,
        swapFixedIvl,
        swapFloatIvl,
        swapBasisIvl,
        swapBadDayConvention,
        floatRateFixed,
        fixDates,
        fixRates,
        volModelIR,
        holidays,
        ccyBasis,
        futuresMaturity);

    Business252     bus252DCC(holidays);
    BadDayFollowing bdf;
    BadDayNone      bdn;
    DateTime        baseDate = spotDate;
    DateTime        date;
    double          rate;

    vector<DateTime> relDates;
    vector<double>   relSpreads;
    vector<DateTime> benchDates;
    vector<double>   benchRates;

    LMParabCurveSP basePC;
    ZC3ZeroCurveSP baseZC;

    // For backwards compatibility with an old version of ZeroCurve3
    // we need to be able to specify that cash maturity dates are bad
    // day adjusted.
    const BadDayConvention& bdc = adjustDates ?
        static_cast<BadDayConvention&>(bdf) :
        static_cast<BadDayConvention&>(bdn);
 
    // compute dates and discounts associated with ON rates
    double ONRate = benchmarks[0]->getRate();
    double SNRate = benchmarks[1]->getRate();
    double ONBus252 = pow(1. + ONRate/30., 252.) - 1.;

    date = bdc.adjust(benchmarks[1]->getBenchmarkDate(baseDate), holidays.get());
    rate = pow(1. + SNRate/30., 252.) - 1.;
    benchDates.push_back(date);
    benchRates.push_back(rate);

    // compute dates and BUS252 rates associated with swaps and futures
    int i;
    Actual360 act360;
    for (i = 2; i < benchmarks.size(); ++i)
    {
        if (adjustDates || !benchmarks[i]->getEnd())
        {
            // For backwards compatibility with an old version of ZeroCurve3,
            // or if an end expiry has not been provided, futuresMaturity is used
            // to calculate the expiry from the future's start date.
            date = benchmarks[i]->getStart()->toDate(baseDate);
            date = futuresMaturity->toDate(bdc.adjust(date, holidays.get()));
        }
        else
        {            
            date = benchmarks[i]->getEnd()->toDate(baseDate);
        }

        date = holidays->addBusinessDays(date, 1);

        if (benchmarks[i]->getIncludeFlag() == 2)
        {
            // Include flag of 2 denotes a relative benchmark
            relDates.push_back(date);
            relSpreads.push_back(benchmarks[i]->getAdjustment());
        }
        else if (benchmarks[i]->isFuture())
        {
            double price = adjustDates ?
                benchmarks[i]->getRate() :
                (1. - benchmarks[i]->getRate())*FuturesBenchmark::ZERO_PERCENT;

            rate = RateConversion::discountToRate(
                price / BRZ_FUTURES_PRICE_UNIT,
                baseDate,
                date,
                &bus252DCC,
                CompoundBasis::ANNUAL);

            benchDates.push_back(date);
            benchRates.push_back(rate);
        }
        else if (benchmarks[i]->isSwap())
        {
            rate = RateConversion::rateConvert(
                tradeDate,
                benchmarks[i]->getBenchmarkDate(baseDate),
                benchmarks[i]->getRate(),
                &act360,
                CompoundBasis::ANNUAL,
                &bus252DCC,
                CompoundBasis::ANNUAL);

            benchDates.push_back(date);
            benchRates.push_back(rate);
        }
    }

    // compute the flat dates
    for (i = 0; i < copomDates->size() && (*copomDates)[i] <= tradeDate; ++i);

    vector<DateTime> flatRateDates;
    flatRateDates.reserve(copomDates->size() + 1 - i);

    if( i < copomDates->size() && benchmarks.size() > 1)
    {
        flatRateDates.push_back(holidays->addBusinessDays((*copomDates)[i],1));

        if (flatRateDates.front() >= benchDates[1])
        {
            flatRateDates[0] = holidays->addBusinessDays(baseDate,1);
        }
        else
        {
            ++i;
        }
    }
    else
    {
        flatRateDates.push_back(holidays->addBusinessDays(baseDate,1));
    }

    for (int j = i; j < copomDates->size(); ++j)
    {
        flatRateDates.push_back(holidays->addBusinessDays((*copomDates)[j], 1));
    }

    DateTime lastFlatDate = flatRateDates.back();

    //
    // Create the first zero curve which only uses benchmark prices and dates.
    // The first curve has num futures + num swaps + 1 points.
    //

    basePC = LMParabCurve::genFromSwaps(baseDate,
                                        &benchDates,
                                        &benchRates,
                                        NULL,
                                        NULL,
                                        NULL,
                                        NULL,
                                        NULL,
                                        NULL,
                                        NULL,
                                        &relDates,
                                        &relSpreads,
                                        bus252DCC,
                                        CompoundBasis::ANNUAL,
                                        &flatRateDates,
                                        NULL,
                                        NULL,
                                        0,
                                        bus252DCC,
                                        BadDayNone(),
                                        BadDayNone(),
                                        lastFlatDate,
                                        holidays,
                                        ZC3ZeroInterpolation::SMOOTH_FORWARDS);

    baseZC = basePC->toZeroCurve();

    // add ONRate as a stub rate

    ZC3StubArray stub(1, ZC3Stub(tradeDate,
                                 baseDate,
                                 ONBus252,
                                 *ZC3Stub::STUB_DCC,
                                 0.,
                                 *ZC3Stub::STUB_DCC,
                                 BadDayNone(),
                                 *holidays.get()));

    baseZC->addStub(stub,
                    bus252DCC,
                    CompoundBasis::ANNUAL);

    // For compatibility with ALIB, insert trade date as first critical date
    baseZC->criticalDates.insert(baseZC->criticalDates.begin(), tradeDate);

    return ZeroPairSP(new ZeroPair(baseZC, baseZC));
}


/*
 * Reflection support.
 */

static IObject* defaultZCBrzFI()
{
    return new ZCBrzFI();
}


/** Invoked when Class is 'loaded' */
void ZCBrzFI::load(CClassSP& clazz)
{
    clazz->setPublic(); // make visible to EAS/spreadsheet
    REGISTER(ZCBrzFI, clazz);
    SUPERCLASS(IZeroCurveFactory);
    EMPTY_SHELL_METHOD(defaultZCBrzFI);

    FIELD(copomDates, "Copom Dates");

    Addin::registerConstructor("ZC_BRZ_FI",
                               Addin::MARKET,
                               "Creates factory for bootstrapping brazilian curves",
                               ZCBrzFI::TYPE);
}


CClassConstSP const ZCBrzFI::TYPE = 
    CClass::registerClassLoadMethod("ZCBrzFI", typeid(ZCBrzFI), load);



DRLIB_END_NAMESPACE
