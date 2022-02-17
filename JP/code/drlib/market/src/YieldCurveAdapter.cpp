#include "edginc/config.hpp"
#include "edginc/YieldCurveAdapter.hpp"

DRLIB_BEGIN_NAMESPACE


YieldCurveAdapter::YieldCurveAdapter(const IYieldCurve& pYc, bool pProjection)
: ZeroCurve(TYPE), yc(const_cast<IYieldCurve*>(&pYc)), projection(pProjection)
{
}

YieldCurveAdapter::~YieldCurveAdapter()
{
}


const IYieldCurve& YieldCurveAdapter::get() const
{
    if (!yc.get())
    {
        throw ModelException("YieldCurveAdapter::get()","Yield curve not set");
    }

    yc->setProjectionCurve(projection);
    return *yc;
}


double YieldCurveAdapter::discountFactor(const DateTime& date) const
{
    return get().pv(date);
}


double YieldCurveAdapter::zeroCouponRate(const DateTime& date) const
{
    return get().zero(date);
}


int YieldCurveAdapter::length() const
{
    throw ModelException("YieldCurveAdapter::length()", "Not implemented");
}


const DateTime& YieldCurveAdapter::firstDate() const
{
    // generally true
    firstDt = get().getSpotDate();
    return firstDt;
}


const DateTime& YieldCurveAdapter::endDate() const
{
    throw ModelException("YieldCurveAdapter::endDate()", "Not implemented");
}


DateTimeArray YieldCurveAdapter::getDates() const
{
    return get().zeroDates();
}


CashFlowArraySP YieldCurveAdapter::getRatesAndDates() const
{
    return get().getRatesAndDates();
}


const DateTime& YieldCurveAdapter::getBaseDate() const
{
    baseDate = get().getSpotDate();
    return baseDate;
}


YieldCurve::IKey* YieldCurveAdapter::logOfDiscFactorKey() const
{
    return get().logOfDiscFactorKey();
}


// for reflection

YieldCurveAdapter::YieldCurveAdapter() : ZeroCurve(TYPE), projection(false)
{
}

/* static */
void YieldCurveAdapter::load(CClassSP& clazz)
{
    REGISTER(YieldCurveAdapter, clazz);
    SUPERCLASS(ZeroCurve);
    FIELD       (yc,         "Yield curve to treat as a zero curve");
    FIELD(projection, "True if using estimating curve");
}

CClassConstSP const YieldCurveAdapter::TYPE = CClass::registerClassLoadMethod(
    "YieldCurveAdapter", typeid(YieldCurveAdapter), YieldCurveAdapter::load);


DRLIB_END_NAMESPACE
