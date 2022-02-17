//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : GenericSimpleIR.cpp
//
//   Description : Base class for simple interest rate generic instruments
//
//   Author      : Andrew J Swain
//
//   Date        : 18 January 2002
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/GenericSimpleIR.hpp"

DRLIB_BEGIN_NAMESPACE

/* for reflection */
GenericSimpleIR::GenericSimpleIR(CClassConstSP clazz): CInstrument(clazz){
    // empty
}

GenericSimpleIR::~GenericSimpleIR(){}

/** Get the asset and discount market data */
void GenericSimpleIR::GetMarket(const IModel*          model, 
                                const CMarketDataSP    market)
{
    market->GetReferenceDate(valueDate);
    discount.getData(model, market);
    // vol is optional for now - v.weak
    if (!vol.isEmpty()) {
        vol.getData(model, market);
    }
    coupon = YieldCurveSP(discount.getSP().clone());
    coupon->setProjectionCurve(); // use 'growth' zc
}

DateTime GenericSimpleIR::getValueDate() const {
    return valueDate;
}

// roll through time 
bool GenericSimpleIR::sensShift(Theta* theta){
    // roll today 
    valueDate = theta->rollDate(valueDate);
    return true; // continue to tweak components which implement Theta
}


/** Returns the name of the instrument's discount currency. */
string GenericSimpleIR::discountYieldCurveName() const {
    return discount.getName();
}


class GenericSimpleIRHelper {
public:
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(GenericSimpleIR, clazz);
        SUPERCLASS(CInstrument);
        IMPLEMENTS(Theta::Shift);
        FIELD(valueDate,"valuation Date");
        FIELD_MAKE_OPTIONAL(valueDate);
        FIELD(discount,"discount curve");
        FIELD_NO_DESC(coupon);
        FIELD_MAKE_TRANSIENT(coupon);
        FIELD(vol,"volatility");
        FIELD_MAKE_OPTIONAL(vol); // weak - but allows futures as generics
    }
};

CClassConstSP const GenericSimpleIR::TYPE = CClass::registerClassLoadMethod(
    "GenericSimpleIR", typeid(GenericSimpleIR), GenericSimpleIRHelper::load);

DRLIB_END_NAMESPACE

