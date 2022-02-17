//----------------------------------------------------------------------------
//
//   Group       : Quantitative Research
//
//   Filename    : UntweakableBasisIndexCurve.hpp
//
//   Description : Basis index curve that is already bootstrapped and ready
//                 for use by the diffusion model
//
//   Author      : 
//
//   Date        : 23rd Aug 2006
//
//
//---------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/UntweakableBasisIndexCurve.hpp"
#include "edginc/CompoundBasis.hpp"
#include "edginc/DayCountConventionFactory.hpp"
#include "edginc/SwapTool.hpp"
#include "edginc/Maths.hpp"
#include "edginc/Format.hpp"
#include "edginc/Addin.hpp"
#include "edginc/ValueDateCollector.hpp"
#include "edginc/WrapperNameCollector.hpp"
#include "edginc/BadDayConventionFactory.hpp"
#include "edginc/Equals.hpp"
/*
#include "edginc/ZeroCurve.hpp"
#include "edginc/BootstrappedYieldCurve.hpp"
#include "edginc/RateConversion.hpp"
#include "edginc/SecantBrentRootFinder.hpp"
#include "edginc/StubPlacement.hpp"
#include "edginc/StubSimple.hpp"
*/


DRLIB_BEGIN_NAMESPACE


/** overrides CObject version to allow for easy default */
bool UntweakableBasisIndexCurve::accept(ICollector* collector) const
{
    if (!CClass::invokeAcceptMethod(this, collector))
    {
        // if no method registered try  vol
        return irVol->accept(collector);
    }

    return false;
}

void UntweakableBasisIndexCurve::acceptValueDateCollector(
    const UntweakableBasisIndexCurve* ic,
    CValueDateCollector*               collector)
{
    collector->valueDateValidate(ic->today, ic->getName());
}

void UntweakableBasisIndexCurve::acceptWrapperNameCollector(
    const UntweakableBasisIndexCurve* ic,
    WrapperNameCollector*              collector)
{
    collector->addName(ic->getName());
}

/** Returns an processed vol - which combines the vol market data with the
instrument data in the volRequest */
IVolProcessed* UntweakableBasisIndexCurve::getProcessedVol(const CVolRequest* volRequest) const
{
    // (1) see whether request is for irVol or for spVol:
    if (SPCalib::recognizedVolRequest(volRequest))
    {
        if (spVol.isEmpty())
            throw ModelException("UntweakableBasisIndexCurve::getProcessedVol",
            "Pricing model requires SP Vol to be supplied for "
            "UntweakableBasisIndexCurve with name " + getName());

        return spVol->getProcessedVol(volRequest);
    }
    else // default to irVol
    {
        if (irVol.isEmpty())
            throw ModelException("UntweakableBasisIndexCurve::getProcessedVol",
            "Pricing model requires IR Vol be supplied for "
            "UntweakableBasisIndexCurve with name " + getName());

        return irVol->getProcessedVol(volRequest, this);
    }
}

/** populate from market cache */
void UntweakableBasisIndexCurve::getMarket(
    const IModel* model,
    const MarketData* market )
{
    if (today.empty())
    {
        today = market->GetReferenceDate();
    }

    holidays.getData(model, market);

    if ( !reference.isEmpty() )
    {
        reference.getData(model, market);
    }

    /*
    if ( !discount.isEmpty() )
    {
        discount.getData(model, market);
    }
    */

    // if specified, ask for the sp vol
    if ( !spVol.isEmpty())
    {
        spVol.getData(model, market);
    }

    // if specified, ask for the vol
    if ( !irVol.isEmpty())
    {
        irVol.getData(model, market); // NB This might result in a null vol
        // if the model reckons we don't need it
    }
}


/* for reflection */
UntweakableBasisIndexCurve::UntweakableBasisIndexCurve(
    const DateTime&     today, 
    const DateTime&     valueDate, 
    IYieldCurveConstSP  basisCurve,
    DoubleArraySP       parSwapRates)
    : MarketObject(TYPE), 
    today(today), 
    valueDate(valueDate),
    spotOffset(valueDate.daysDiff(today)), 
    basisCurve(basisCurve), 
    parSwapRates(parSwapRates)
{}

UntweakableBasisIndexCurve::UntweakableBasisIndexCurve(): 
MarketObject(TYPE), spotOffset(0), isAdditive(true)
{}

bool UntweakableBasisIndexCurve::equalTo(const IObject* curve) const
{
    try
    {
        if (this == curve)    // obvious first test
            return true;

        if (!curve || getClass() != curve->getClass())
            return false;

        const UntweakableBasisIndexCurve* ic2 = STATIC_CAST(UntweakableBasisIndexCurve, curve);

        if (!today.equals(ic2->today))
            return false;

        if (isAdditive != ic2->isAdditive)
            return false;

        if (spotOffset != ic2->spotOffset)
            return false;

        /*
        if (!GenericArrayEquals<ExpiryArray>(expiries, ic2->expiries))
            return false;

        if (!NumericArrayEquals<DoubleArray>(spreads, ic2->spreads))
            return false;

        if (!StringArrayEquals<StringArray>(instruments, ic2->instruments))
            return false;

        // ... NB. no comparison function on MarketDataWrapper
        if (discount.get())
        {
            if (!discount.get()->zeroCurveEquals(ic2->discount.get()))
                return false;
        }
        else
        {
            if (ic2->discount.get())
                return false;
        }
        */

        if (basisCurve.get())
        {
            if (!basisCurve.get()->zeroCurveEquals(ic2->basisCurve.get()))
                return false;
        }
        else
        {
            if (ic2->basisCurve.get())
                return false;
        }

        if (reference.get())
        {
            if (!reference.get()->zeroCurveEquals(ic2->reference.get()))
                return false;
        }
        else
        {
            if (ic2->reference.get())
                return false;
        }

        /*
        if (!equalsValue<IZeroCurveFactory>(zeroCurveFactory, ic2->zeroCurveFactory))
            return false;

        // assume day count conventions are unique by class - should really
        // put equals method on them
        if (!equalsClass<DayCountConvention>(moneyMarketDcc, ic2->moneyMarketDcc))
            return false;
            */

        if (!equalsClass<BadDayConvention>(badDayConvention, ic2->badDayConvention))
            return false;

        if (!equalsClass<DayCountConvention>(basisSwapLiborDcc, ic2->basisSwapLiborDcc))
            return false;

        if (!equalsClass<DayCountConvention>(basisSwapBasisDcc, ic2->basisSwapBasisDcc))
            return false;

        if (!equalsValue<MaturityPeriod>(basisSwapLiborIvl, ic2->basisSwapLiborIvl))
            return false;

        if (!equalsValue<MaturityPeriod>(basisSwapBasisIvl, ic2->basisSwapBasisIvl))
            return false;

        /*
        if (!equalsClass<DayCountConvention>(refSwapFixedDcc, ic2->refSwapFixedDcc))
            return false;

        if (!equalsClass<DayCountConvention>(refSwapFixedDcc, ic2->refSwapFixedDcc))
            return false;

        if (!equalsValue<MaturityPeriod>(refSwapFixedIvl, ic2->refSwapFixedIvl))
            return false;

        if (!equalsValue<MaturityPeriod>(refSwapFloatIvl, ic2->refSwapFloatIvl))
            return false;
            */

        // this code also gets triggered on construction, so may not have data
        // from market cache yet
        if (badDayConvention.get())
        {
            if (holidays.get() && !holidays->equals(ic2->holidays.get()))
            {
                return false;
            }
        }

        return true;
    }
    catch (exception& e)
    {
        throw ModelException(e, "UntweakableBasisIndexCurve::equalTo");
    }
}

void UntweakableBasisIndexCurve::validate(const string& method) const
{
    // validate no data is missing

    if (today.empty())
    {
        throw ModelException(method, "Today is not defined");
    }

    if (name.empty())
    {
        throw ModelException(method, "Name is not defined");
    }

    if (reference.isEmpty())
    {
        throw ModelException(method, "Reference curve is not defined");
    }

    if (!basisCurve.get())
    {
        throw ModelException(method, "Basis zero curve is not defined");
    }

    /*
    if (discount.isEmpty())
    {
        throw ModelException(method, "Discount curve is not defined");
    }

    if (!expiries.get() || expiries->empty())
    {
        throw ModelException(method, "Expiries are not defined");
    }

    if (!spreads.get() || spreads->empty())
    {
        throw ModelException(method, "Spreads are not defined");
    }

    if (!instruments.get() || instruments->empty())
    {
        throw ModelException(method, "Instruments are not defined");
    }

    if (!zeroCurveFactory.get())
    {
        throw ModelException(method, "Basis curve zero curve factory is not defined");
    }

    if (!moneyMarketDcc.get())
    {
        throw ModelException(method, "Money market day count convention is not defined");
    }

    if (!refSwapFixedDcc.get())
    {
        throw ModelException(method, "Reference swap fixed leg day count convention is not defined");
    }

    if (!refSwapFloatDcc.get())
    {
        throw ModelException(method, "Reference swap floating leg day count convention is not defined");
    }

    if (!refSwapFixedIvl.get())
    {
        throw ModelException(method, "Reference swap fixed leg reset period is not defined");
    }

    if (!refSwapFloatIvl.get())
    {
        throw ModelException(method, "Reference swap floating leg reset period is not defined");
    }
    */

    if (!basisSwapLiborDcc.get())
    {
        throw ModelException(method, "Basis swap libor day count convention is not defined");
    }

    if (!basisSwapBasisDcc.get())
    {
        throw ModelException(method, "Basis swap basis day count convention is not defined");
    }

    if (!basisSwapLiborIvl.get())
    {
        throw ModelException(method, "Basis swap libor reset period is not defined");
    }

    if (!basisSwapBasisIvl.get())
    {
        throw ModelException(method, "Basis swap basis reset period is not defined");
    }

    if (!badDayConvention.get())
    {
        throw ModelException(method, "Bad day convention is not defined");
    }

    /*
    // validate array sizes
    if (expiries->size() != spreads->size())
    {
        string msg = Format::toString(
            "Number of expiries (%d) differs from number of spreads (%d)",
            expiries->size(), spreads->size());
        throw ModelException(method, msg);
    }

    if (instruments->size() != spreads->size())
    {
        string msg = Format::toString(
            "Number of instruments (%d) differs from number of spreads (%d)",
            instruments->size(), spreads->size());
        throw ModelException(method, msg);
    }

    if (!reference->getCcy().empty() && !discount->getCcy().empty()
        && reference->getCcy() != discount->getCcy())
    {
        string msg = "reference curve currency [" +
            reference->getCcy() +
            "] differs from discount curve currency [" +
            discount->getCcy() + "]";
        throw ModelException(method, msg);
    }

    if (reference->getToday() != discount->getToday())
    {
        string msg = "reference curve value date [" +
            reference->getToday().toString() +
            "] differs from discount curve value date [" +
            discount->getToday().toString() + "]";
        throw ModelException(method, msg);
    }
    */

    if (!reference->getCcy().empty() && !basisCurve->getCcy().empty()
        && reference->getCcy() != basisCurve->getCcy())
    {
        string msg = "reference curve currency [" +
            reference->getCcy() +
            "] differs from basisCurve curve currency [" +
            basisCurve->getCcy() + "]";
        throw ModelException(method, msg);
    }

    if (reference->getToday() != basisCurve->getToday())
    {
        string msg = "reference curve value date [" +
            reference->getToday().toString() +
            "] differs from basisCurve curve value date [" +
            basisCurve->getToday().toString() + "]";
        throw ModelException(method, msg);
    }
}

CClassConstSP const UntweakableBasisIndexCurve::TYPE = CClass::registerClassLoadMethod(
    "UntweakableBasisIndexCurve", typeid(UntweakableBasisIndexCurve), load);

/* external symbol to allow class to be forced to be linked in */
bool UntweakableBasisIndexCurveLinkIn(){
    return true;
}


/* static */
IObject* defaultUntweakableBasisIndexCurve()
{
    return new UntweakableBasisIndexCurve();
}

void UntweakableBasisIndexCurve::load(CClassSP& clazz)
{
    clazz->setPublic(); // make visible to EAS/spreadsheet
    REGISTER(UntweakableBasisIndexCurve, clazz);
    SUPERCLASS(MarketObject);
    IMPLEMENTS(IBasisIndexCurve);
    IMPLEMENTS(IGetMarket);
    IMPLEMENTS(IMarketFactor);
    EMPTY_SHELL_METHOD(defaultUntweakableBasisIndexCurve);

    // basis index curve parameters
    FIELD(name,             "basis curve name");
    FIELD(today,            "today");
    FIELD(valueDate,        "value date");
    FIELD(isAdditive,       "true if spread basis swap type - else percentage");
    FIELD(reference,        "reference zero curve");
    FIELD       (basisCurve,        "bootstrapped basis zero curve");
    FIELD       (parSwapRates,      "bootstrapped par swap rates [optional]");
    FIELD(spotOffset,        "spot offset");
    FIELD(holidays,          "holidays");
    FIELD(irVol,             "interest rate volatility");
    FIELD(spVol,             "spread volatility and mean reversion");
    FIELD       (badDayConvention,  "bad day convention");
    FIELD       (basisSwapLiborDcc, "basis swap libor leg day count convention");
    FIELD       (basisSwapBasisDcc, "basis swap basis leg day count convention");
    FIELD       (basisSwapLiborIvl, "basis swap libor leg period");
    FIELD       (basisSwapBasisIvl, "basis swap basis leg period");

    FIELD_MAKE_OPTIONAL(parSwapRates);
    FIELD_MAKE_OPTIONAL(irVol);
    FIELD_MAKE_OPTIONAL(spVol);

    /*
    FIELD       (expiries,         "basis swap dates");
    FIELD       (spreads,          "basis swap spreads");
    FIELD       (instruments,      "instrument types - (S)wap or (M)oney market");
    FIELD(discount,         "discount zero curve");
    FIELD       (moneyMarketDcc,    "money market day count convention");
    FIELD       (zeroCurveFactory,  "zero curve factory for building basis curve");
    FIELD       (refSwapFixedDcc,   "reference swap fixed leg day count convention");
    FIELD       (refSwapFloatDcc,   "reference swap floating leg day count convention");
    FIELD       (refSwapFixedIvl,   "reference swap fixed leg period");
    FIELD       (refSwapFloatIvl,   "reference swap floating leg period");

    FIELD_MAKE_OPTIONAL(today);
    */

    ClassSetAcceptMethod(UntweakableBasisIndexCurve::acceptValueDateCollector);
    ClassSetAcceptMethod(UntweakableBasisIndexCurve::acceptWrapperNameCollector);
}


DateTime UntweakableBasisIndexCurve::getSpotDate() const
{
    return basisCurve->getSpotDate();
}

DateTime UntweakableBasisIndexCurve::getFirstCurveDate() const
{
    return getSpotDate().rollDate(spotOffset+1);
}

DateTime UntweakableBasisIndexCurve::getRefPaymentDate(
    const DateTime& resetDate ) const
{
    DateTime payDate = basisSwapLiborIvl->toDate( resetDate );
    return badDayConvention->adjust(payDate, holidays.get());
}


double UntweakableBasisIndexCurve::fwd(const DateTime& startDate) const
{
    DateTime endDate = basisSwapBasisIvl->toDate(startDate);
    return fwd(startDate, endDate, basisSwapBasisDcc.get(), CompoundBasis::SIMPLE);
}

double UntweakableBasisIndexCurve::fwd(
    const DateTime&           lodate,
    const DateTime&           hidate,
    const DayCountConvention* dcc,
    int                       basis) const
{
    return basisCurve->fwd(lodate, hidate, dcc, basis);
}

DoubleArraySP UntweakableBasisIndexCurve::getParSwapRates() const
{
    if (!parSwapRates.get())
        throw ModelException("UntweakableBasisIndexCurve::getParSwapRates()", 
            "parSwapRates not supplied by user");

    return parSwapRates;
}

double UntweakableBasisIndexCurve::parSwapRate(const DateTime& maturity) const
{
    return basisCurve->parSwapRate(maturity);
}

double UntweakableBasisIndexCurve::parSpread(const DateTime& maturity) const
{
    double stdRate = reference->parSwapRate(maturity);
    double basisRate = parSwapRate(maturity);
    return isAdditive ? (stdRate - basisRate) : (stdRate / basisRate);
}

double UntweakableBasisIndexCurve::fwdSpread(const DateTime& reset, const DateTime& maturity) const
{
    double liborRate = reference->fwd(reset, maturity, basisSwapLiborDcc.get(), CompoundBasis::SIMPLE); // TODO: is the DCC here correct?
    double basisRate = fwd(reset, maturity, basisSwapBasisDcc.get(), CompoundBasis::SIMPLE);
    return isAdditive ? (liborRate - basisRate) : (liborRate / basisRate);
}


string UntweakableBasisIndexCurve::getName() const
{
    return name;
}

int UntweakableBasisIndexCurve::hashCode() const
{
    int hCode = (size_t) getClass();
    hCode ^= today.hashCode();
    hCode ^= CBool::hashCode(isAdditive);
    hCode ^= spotOffset;

    /*
    if (expiries.get())          hCode ^= expiries->hashCode();
    if (spreads.get())           hCode ^= spreads->hashCode();
    if (instruments.get())       hCode ^= instruments->hashCode();
    if (discount.get())          hCode ^= discount->hashCode();
    if (zeroCurveFactory.get())  hCode ^= zeroCurveFactory->hashCode();
    if (moneyMarketDcc.get())    hCode ^= moneyMarketDcc->hashCode();
    if (refSwapFixedDcc.get())   hCode ^= refSwapFixedDcc->hashCode();
    if (refSwapFloatDcc.get())   hCode ^= refSwapFloatDcc->hashCode();
    if (refSwapFixedIvl.get())   hCode ^= refSwapFixedIvl->hashCode();
    if (refSwapFloatIvl.get())   hCode ^= refSwapFloatIvl->hashCode();
    */
    if (reference.get())         hCode ^= reference->hashCode();
    if (basisSwapLiborDcc.get()) hCode ^= basisSwapLiborDcc->hashCode();
    if (basisSwapBasisDcc.get()) hCode ^= basisSwapBasisDcc->hashCode();
    if (basisSwapLiborIvl.get()) hCode ^= basisSwapLiborIvl->hashCode();
    if (basisSwapBasisIvl.get()) hCode ^= basisSwapBasisIvl->hashCode();
    if (badDayConvention.get())  hCode ^= badDayConvention->hashCode();

    return hCode;
}

YieldCurveWrapper UntweakableBasisIndexCurve::getRefCurve() const
{
    return reference;
}

bool UntweakableBasisIndexCurve::getIsAdditive() const
{
    return isAdditive;
}

DayCountConventionConstSP UntweakableBasisIndexCurve::getBasisSwapLiborDcc() const
{
    return basisSwapLiborDcc;
}

DayCountConventionConstSP UntweakableBasisIndexCurve::getBasisSwapBasisDcc() const
{
    return basisSwapBasisDcc;
}

MaturityPeriodConstSP UntweakableBasisIndexCurve::getBasisSwapLiborIvl() const
{
    return basisSwapLiborIvl;
}

MaturityPeriodConstSP UntweakableBasisIndexCurve::getBasisSwapBasisIvl() const
{
    return basisSwapBasisIvl;
}

BadDayConventionConstSP UntweakableBasisIndexCurve::getBadDayConvention() const
{
    return badDayConvention;
}

const HolidayWrapper& UntweakableBasisIndexCurve::getHolidays() const
{
    return holidays;
}

DRLIB_END_NAMESPACE