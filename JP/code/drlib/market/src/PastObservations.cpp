//----------------------------------------------------------------------------
//
//   Group       : QR
//
//   Filename    : PastObservations.cpp
//
//   Description : Records historic (overriding) values for an asset
//                 including source/observation type for centralised sampling
//
//   Author      : Ian Stares
//
//   Date        : 15 May 2006
//

//
//----------------------------------------------------------------------------
#include "edginc/config.hpp"
#include "edginc/PastObservations.hpp"
#include "edginc/MarketObservable.hpp"

DRLIB_BEGIN_NAMESPACE

void PastObservations::load(CClassSP& clazz){
    REGISTER(PastObservations, clazz);
    SUPERCLASS(CObject);
    EMPTY_SHELL_METHOD(defaultPastObs);
    FIELD(samples, "The overriding past samples for this underlying");
    FIELD(source, "The source for centralised sampling for this underlying");
    FIELD(obsType, "The source for centralised sampling for this underlying");
    FIELD(fxSource, "The source for FX centralised sampling for this underlying");
    FIELD_MAKE_OPTIONAL(fxSource);
    FIELD(fxObsType, "The source for FX centralised sampling for this underlying");
    FIELD_MAKE_OPTIONAL(fxObsType);
    // Register conversion magic to allow a raw CashFlowArray to be
    // transformed into an PastObservations
    registerObjectFromArrayMethod(CashFlowArray::TYPE,
                                  TYPE,
                                  &fromCashFlowArray);
    clazz->setPublic(); // make visible to EAS/spreadsheet
}

PastObservations::PastObservations(): CObject(TYPE) {}

PastObservations::PastObservations(CashFlowArray& smpls) :
                        CObject(TYPE), samples(smpls) {
    source = IMarketObservable::DEFAULT_SOURCE;
    obsType = ObservationType::make("NotUsed");
}

IObjectSP PastObservations::fromCashFlowArray(const IObjectSP& object, 
                                              CClassConstSP    requiredType) {
    CashFlowArray& samples = dynamic_cast<CashFlowArray&>(*object);
    return IObjectSP(new PastObservations(samples));
}

const CashFlowArray& PastObservations::getSamples() {
    return samples;
}

IObject* PastObservations::defaultPastObs() {
    return new PastObservations();
}

CClassConstSP const PastObservations::TYPE = CClass::registerClassLoadMethod(
    "PastObservations", typeid(PastObservations), PastObservations::load);

DEFINE_TEMPLATE_TYPE(PastObservationsArray);

DRLIB_END_NAMESPACE
