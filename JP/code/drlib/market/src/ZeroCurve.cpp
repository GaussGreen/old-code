//----------------------------------------------------------------------------
//
//   Group       : QR&D Core Analytics Team
//
//   Filename    : ZeroCurve.cpp
//
//   Description : Abstract base class for zero curves
//
//   Author      : Richard Appleton
//
//   Date        : 4th July 2005
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/ZeroCurve.hpp"
#include "edginc/Maths.hpp"
#include "edginc/ModelException.hpp"
#include "edginc/CompoundBasis.hpp"
#include "edginc/RateConversion.hpp"
#include "edginc/Actual365F.hpp"
#include "edginc/CashFlow.hpp"
#include "edginc/Addin.hpp"
#include "edginc/DayCountConventionFactory.hpp"
#include "edginc/StubFactory.hpp"
#include "edginc/BadDayConventionFactory.hpp"
#include "edginc/SwapTool.hpp"
#include "edginc/NonPricingModel.hpp"

#include <math.h>

DRLIB_BEGIN_NAMESPACE


ZeroCurve::ZeroCurve(CClassConstSP clazz)
: CObject(clazz)
{
}


ZeroCurve::~ZeroCurve()
{
}


double ZeroCurve::pv(const DateTime& lodate, const DateTime& hidate) const
{
    double lo = discountFactor(lodate);
    double hi = discountFactor(hidate);
    return hi / lo;
}


double ZeroCurve::pv(const CashFlowArray& cfl) const
{
    double result = 0.0;
    for (int i = 0; i < cfl.size(); i++)
    {
        result += cfl[i].amount * discountFactor(cfl[i].date);
    }

    return result;
}


double ZeroCurve::fv(const CashFlowArray& cfl, const DateTime& date) const
{
    double forwardValue = pv(cfl);
    if (!Maths::isZero(forwardValue))
    {
        forwardValue /= discountFactor(date);
    }

    return forwardValue;
}


double ZeroCurve::fwd(
    const DateTime&           lodate, 
    const DateTime&           hidate,
    const DayCountConvention* dcc,
    int                       basis) const
{
    double disc = pv(lodate, hidate);
    return RateConversion::discountToRate(disc, lodate, hidate, dcc, basis);
}

double ZeroCurve::parSwapRate(
                   const DateTime&           startDate, 
                   const DateTime&           endDate,
                   const MaturityPeriod&     period, 
                   const DayCountConvention& dcc,
                   const Stub&               stubType,
                   bool                      stubAtEnd,
                   const BadDayConvention&   accBadDayConv,
                   const BadDayConvention&   payBadDayConv,
                   const Holiday&            holidays) const 
{
    return SwapTool::parSwapRate(*this, startDate, endDate, period, dcc, stubType, stubAtEnd, accBadDayConv, payBadDayConv, holidays);
}

class ZeroCurve::LogOfDiscFactorKey: public YieldCurve::IKey
{
public:
    LogOfDiscFactorKey(const ZeroCurve& zc) : zc(zc)
    {
    }

    /**
     * Default implementation for this is for no optimization - just take the 
     * log of the pv factor between the two dates.
     */
    virtual double calc(const DateTime&  newLoDate, const DateTime&  newHiDate)
    {
        return log(zc.pv(newLoDate, newHiDate));
    }

private:
    const ZeroCurve& zc;
};


/**
 * Returns a key used to optimize repeated calculations of forward rates.
 */
YieldCurve::IKey* ZeroCurve::logOfDiscFactorKey() const
{
    return new LogOfDiscFactorKey(*this);
}


class ZeroCurveHelper
{
public:
    /**
     * Invoked when Class is 'loaded' - we only really need the reflection
     * info to serialize the class.
     */
    static void load(CClassSP& clazz)
    {
        REGISTER(ZeroCurve, clazz);
        SUPERCLASS(CObject);
    }
};


CClassConstSP const ZeroCurve::TYPE = CClass::registerClassLoadMethod(
    "ZeroCurve", typeid(ZeroCurve), ZeroCurveHelper::load);



// Util: curve discount factor
class ZCDiscountFactorAddin: public CObject
{
public:
    static CClassConstSP const TYPE;
    ZeroCurveSP         zeroCrv;
    DateTime            date;
    ZCDiscountFactorAddin(): CObject(TYPE) {}
    static IObject* defaultZCDiscountFactorAddin()
    {
        return new ZCDiscountFactorAddin();
    }
    double discountFactor()
    {
        return zeroCrv->discountFactor(date);
    }
    static void load(CClassSP& clazz)
    {
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(ZCDiscountFactorAddin, clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultZCDiscountFactorAddin);
        FIELD(zeroCrv, "zero curve");
        FIELD(date, "date");
        Addin::registerDoubleMethod("ZC_DISCOUNT_FACTOR",
            Addin::UTILITIES,
            "Returns discount factor from zero curve",
            &ZCDiscountFactorAddin::discountFactor);                                    
    }
};
CClassConstSP const ZCDiscountFactorAddin::TYPE=CClass::registerClassLoadMethod(
    "ZCDiscountFactorAddin", typeid(ZCDiscountFactorAddin), load);


// Util: curve forward rate
class ZCForwardRateAddin: public CObject
{
public:
    static CClassConstSP const TYPE;

    ZeroCurveSP          zeroCrv;
    DateTime             dateStart;
    DateTime             dateEnd;
    DayCountConventionSP dcc;
    string               basis;
    MarketDataConstSP    marketData;



    ZCForwardRateAddin(): CObject(TYPE) {}

    static IObject* defaultZCForwardRateAddin()
    {
        return new ZCForwardRateAddin();
    }

    double forwardRate()
    {
        if (marketData.get())
        {
                // May need to fetch holidays for Business/252 DCC
                marketData->fetchForNonMarketObjects(dcc, 
						     IModelConstSP(new NonPricingModel()),
						     "Business252");
        }
	return zeroCrv->fwd(dateStart, dateEnd, dcc.get(), CompoundBasis::toInt(basis));
    }

    static void load(CClassSP& clazz)
    {
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(ZCForwardRateAddin, clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultZCForwardRateAddin);
        FIELD(zeroCrv, "zero curve");
        FIELD(dateStart, "start date");
        FIELD(dateEnd, "end date");
        FIELD(dcc, "day count convention");
        FIELD(basis, "compound basis");
        FIELD(marketData, "market data");
        FIELD_MAKE_OPTIONAL(marketData);

        Addin::registerDoubleMethod("ZC_FORWARD_RATE",
            Addin::UTILITIES,
            "Returns the forward rate from curve",
            &ZCForwardRateAddin::forwardRate);                                    
    }

};

CClassConstSP const ZCForwardRateAddin::TYPE=CClass::registerClassLoadMethod(
    "ZCForwardRateAddin", typeid(ZCForwardRateAddin), load);

// Util: curve par swap rate
class ZCParSwapRateAddin: public CObject
{
public:
    static CClassConstSP const TYPE;

    ZeroCurveSP          zeroCrv;
    DateTime             startDate;
    DateTime             endDate;
    string               period;
    DayCountConventionSP dcc;
    string               stubRule;
    bool                 stubAtEnd;
    string               accBadDayConv;
    string               payBadDayConv;
    HolidaySP            holidays;
    MarketDataConstSP    marketData;


    ZCParSwapRateAddin(): CObject(TYPE) {}

    static IObject* defaultZCParSwapRateAddin()
    {
        return new ZCParSwapRateAddin();
    }

    double parSwapRate()
    {
      
        if (marketData.get())
        {
                // May need to fetch holidays for Business/252 DCC
	        marketData->fetchForNonMarketObjects(dcc, 
						     IModelConstSP(new NonPricingModel()),
						     "Business252");
        }
      

        MaturityPeriod _period(period);
        StubSP _stub(StubFactory::make(stubRule));
        BadDayConventionSP _accBadDayConv_ptr (BadDayConventionFactory::make(accBadDayConv));
        BadDayConventionSP _payBadDayConv_ptr (BadDayConventionFactory::make(payBadDayConv));

	double rtn = zeroCrv->parSwapRate(startDate, endDate, _period, *dcc.get(), *_stub.get(), 
            stubAtEnd, *_accBadDayConv_ptr.get(), *_payBadDayConv_ptr.get(), *holidays.get());
        return rtn;
    }

    static void load(CClassSP& clazz)
    {
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(ZCParSwapRateAddin, clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultZCParSwapRateAddin);
        FIELD(zeroCrv, "zero curve");
        FIELD(startDate, "start date");
        FIELD(endDate,   "end date");
        FIELD(period, "payment period");
        FIELD(dcc, "day count convention");
        FIELD(stubRule, "stub rule");
        FIELD(stubAtEnd, "stub at end");
        FIELD(accBadDayConv, "accrual bad day convention");
        FIELD(payBadDayConv, "pay bad day convention");
        FIELD(holidays, "holidays");
        FIELD(marketData, "market data");
	FIELD_MAKE_OPTIONAL(marketData);

        Addin::registerDoubleMethod("ZC_PAR_SWAP_RATE",
            Addin::UTILITIES,
            "Returns the forward rate from curve",
            &ZCParSwapRateAddin::parSwapRate);                                    
    }

};

CClassConstSP const ZCParSwapRateAddin::TYPE=CClass::registerClassLoadMethod(
    "ZCParSwapRateAddin", typeid(ZCParSwapRateAddin), load);



DRLIB_END_NAMESPACE

