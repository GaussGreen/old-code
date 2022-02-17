//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : GenericNFactor.cpp
//
//   Description : Base class for N factor generic instruments
//                 If you want to book an equity based N factor instrument
//                 in Pyramid 'Generics' this is the mandatory starting point
//
//   Author      : Andrew J Swain
//
//   Date        : 24 October 2001
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/GenericNFactor.hpp"
#include "edginc/Control.hpp"

#include "edginc/InstrumentSettlement.hpp"
#include "edginc/MultiMarketFactors.hpp"

DRLIB_BEGIN_NAMESPACE

GenericNFactor::~GenericNFactor(){}

/* for reflection */
GenericNFactor::GenericNFactor(CClassConstSP clazz): CInstrument(clazz) {
    // empty
}

GenericNFactor::GenericNFactor(CClassConstSP clazz,
                   const DateTime&          valueDate,
                   double                   notional,         
                   InstrumentSettlementSP   instSettle, 
                   IMultiMarketFactorsSP    assets, 
                   const YieldCurveWrapper& discount) 
    : CInstrument(clazz), valueDate(valueDate),
      notional(notional), instSettle(instSettle),  
      assets(assets), discount(discount)
{
    // empty
}


/** Get the asset and discount market data */
void GenericNFactor::GetMarket(const IModel*          model, 
                               const CMarketDataSP    market)
{
    market->GetReferenceDate(valueDate);
    assets->getMarket(model, market.get());
    discount.getData(model, market);
    instSettle->getMarket(model, market.get());
}

DateTime GenericNFactor::getValueDate() const {
    return valueDate;
}

/** Validate instrument having aquired market data */
void GenericNFactor::Validate(){ 
    // in case we need to do some
}

//// roll through time (setting historic values)
bool GenericNFactor::sensShift(Theta* theta){
    valueDate = theta->rollDate(valueDate);
    return true; // continue to tweak components which implement Theta
}

/** Adds the FWD_AT_MAT if requested */
void GenericNFactor::addRequests(Control*        control,
                                 Results*        results,
                                 const DateTime& maturity) const {

    OutputRequest* request = NULL;
    if (maturity.isGreaterOrEqual(valueDate)) {
        if (control->requestsOutput(OutputRequest::FWD_AT_MAT, request)) {
            assets->recordFwdAtMat(request, results, maturity);
        }
    }
}

/** Returns the name of the instrument's discount currency. */
string GenericNFactor::discountYieldCurveName() const {
    return discount.getName();
}


class GenericNFactorHelper {
public:
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(GenericNFactor, clazz);
        SUPERCLASS(CInstrument);
        IMPLEMENTS(Theta::Shift);
        FIELD(valueDate,"valuation date");
        FIELD_MAKE_OPTIONAL(valueDate);
        FIELD(notional,"notional");
        FIELD(instSettle, "instrument settlement at maturity");
        FIELD(assets,"underlying assets");
        FIELD(discount,"discount curve");
    }

};

CClassConstSP const GenericNFactor::TYPE = CClass::registerClassLoadMethod(
    "GenericNFactor", typeid(GenericNFactor), GenericNFactorHelper::load);

DRLIB_END_NAMESPACE

