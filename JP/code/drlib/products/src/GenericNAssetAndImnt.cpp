//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : GenericNAssetAndImnt.cpp
//
//   Description : Base class for N underlying generic instruments that also
//                 have another instrument as an underlying
//                 If you want to book something like a trail fee
//                 in Pyramid 'Generics' this is the mandatory starting point
//
//   Author      : Andrew J Swain
//
//   Date        : 10 April 2002
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/GenericNAssetAndImnt.hpp"

DRLIB_BEGIN_NAMESPACE

/* for reflection */
GenericNAssetAndImnt::GenericNAssetAndImnt(CClassConstSP clazz): CInstrument(clazz) {
    // empty
}

/** Get the asset and discount market data */
void GenericNAssetAndImnt::GetMarket(const IModel*       model, 
                                     const CMarketDataSP market) {
    try {
        market->GetReferenceDate(valueDate);
        discount.getData(model, market);
        instSettle->getMarket(model, market.get());

        for (int i = 0; i < assets.size(); i++) {
            CAsset::getAssetMarketData(model,
                                       market.get(), 
                                       CAsset::CCY_TREATMENT_NONE,
                                       discount.getName(), 
                                       assets[i]);
        }

        imntModel->getInstrumentAndModelMarket(market.get(), imnt.get());
    }
    catch (exception& e) {
        throw ModelException(e, "GenericNAssetAndImnt::GetMarket");
    }
}

DateTime GenericNAssetAndImnt::getValueDate() const {
    return valueDate;
}

/** Validate instrument having aquired market data */
void GenericNAssetAndImnt::Validate(){ 
    // in case we need to do some
}

//// roll through time (setting historic values)
bool GenericNAssetAndImnt::sensShift(Theta* theta){
    valueDate = theta->rollDate(valueDate);
    return true; // continue to tweak components which implement Theta
}


/** Returns the name of the instrument's discount currency. */
string GenericNAssetAndImnt::discountYieldCurveName() const {
    return discount.getName();
}


class GenericNAssetAndImntHelper {
public:
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(GenericNAssetAndImnt, clazz);
        SUPERCLASS(CInstrument);
        IMPLEMENTS(Theta::Shift);
        FIELD(valueDate,"valuation date");
        FIELD_MAKE_OPTIONAL(valueDate);
        FIELD(instSettle, "instrument settlement at maturity");
        FIELD(assets,"underlying assets");
        FIELD(discount,"discount curve");
        FIELD(imnt,"underlying instrument");
        FIELD(imntModel,"and how to price it");
    }
};

CClassConstSP const GenericNAssetAndImnt::TYPE = CClass::registerClassLoadMethod(
    "GenericNAssetAndImnt", typeid(GenericNAssetAndImnt), GenericNAssetAndImntHelper::load);


DRLIB_END_NAMESPACE

