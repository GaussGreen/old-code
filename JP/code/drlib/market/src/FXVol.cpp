//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : FXVol.cpp
//
//   Description : Class to take ordinary equity vol and make it appear as
//                 an fx vol
//
//   Author      : Mark A Robson
//
//   Date        : 5 Dec 2001
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/Void.hpp"
#include "edginc/FXVol.hpp"
#include "edginc/VolParallel.hpp"
#include "edginc/VolPointwise.hpp"
#include "edginc/RiskProperty.hpp"

DRLIB_BEGIN_NAMESPACE

/** Returns name of vol */
string FXVol::getName() const{
    return name;
}


/** populate from market cache - default implementation provided */
void FXVol::getMarket(const IModel* model, const MarketData* market) {
    vol.getData(model, market);
}

/** Combines market and instrument data together to give a
    Processed Vol */
IVolProcessed* FXVol::getProcessedVol(const CVolRequest* volRequest,
                                      const CAsset*      asset) const{
    return vol->getProcessedVol(volRequest, asset);
}

/** Combines market and instrument data together to give a
    Processed Vol. Here the processed volatility is a processed
    struck volatility ie it reflects the combination of this
    CVolBase together with the supplied FX asset and the
    correlation between this CVolBase and the vol of the
    FX. */
IVolProcessed* FXVol::getProcessedVol(
    const CVolRequest* volRequest,
    const CAsset*      eqAsset,
    const FXAsset*     fxAsset,
    const Correlation* eqFXCorr) const{
    return vol->getProcessedVol(volRequest, eqAsset, fxAsset, eqFXCorr);
}

PDFCalculator* FXVol::getPDFCalculator(
    const PDFRequest* request,
    const CAsset*     asset) const {
    // see if we can pass through to the real vol
    const IPDFCalculator* pdf = dynamic_cast<const IPDFCalculator*>(vol.get());
    if (!pdf) {
        throw ModelException("FXVol::getPDFCalculator",
                             "FX vol " + vol->getName() + " does not support "
                             "pdf calculator");
    }
    return pdf->getPDFCalculator(request, asset);
}


/** Returns name identifying vol for fx vega  */
string FXVol::sensName(FXVega* shift) const{
    OutputNameArrayConstSP names(RiskProperty<VolParallel>().
                                     subjectNames(vol.getSP()));

    if (names->size() > 1) {
        // hmmm
        throw ModelException("FXVol::sensName", vol->getName() + 
                             " vol has > 1 vega name");
    }
    if (names->empty()) { 
        return "";
    }
    return (*names)[0]->toString();
}

/** Shifts the object using given shift */
bool FXVol::sensShift(FXVega* shift){
    PropertyTweakHypothesis<VolParallel>(shift->getShiftSize()).
        applyTo(vol.getSP());
    return false;  // nothing else to shift here
}


/** Returns name identifying vol for fx vega pointwise */
string FXVol::sensName(FXVegaPointwise* shift) const{
    OutputNameArrayConstSP names(
        RiskProperty<VolPointwise>().subjectNames(vol.getSP()));

    if (names->size() > 1) {
        // hmmm
        throw ModelException("FXVol::sensName", vol->getName() + 
                             " vol has > 1 vega name");
    }
    if (names->empty()) { 
        return "";
    }
    return (*names)[0]->toString();
}


/** Returns the array of expiries (ie maturities/benchmark dates) that
    need to be tweaked for this vol */
ExpiryArrayConstSP FXVol::sensExpiries(FXVegaPointwise* shift) const{
    return ExpiryWindow::expiries(
        RiskProperty<VolPointwise>().subjectQualifiers(
            vol.getSP(), shift->getMarketDataName()));
}

/** Shifts the object using given shift */
bool FXVol::sensShift(FXVegaPointwise* shift){
    PropertyTweakHypothesis<VolPointwise>(shift->getShiftSize(),
                                          shift->getMarketDataName(),
                                          shift->getExpiryWindow()).
        applyTo(vol.getSP());
    return false;  // nothing else to shift here    
}

FXVol::~FXVol(){}

FXVol::FXVol(const CVolBase* vol): FXVolBase(TYPE), name(vol->getName()),
                                   vol(copy(vol)){}

IObject* FXVol::defaultFXVol(){
    return new FXVol();
}

FXVol::FXVol(): FXVolBase(TYPE){};

/** Invoked when Class is 'loaded' */
void FXVol::load(CClassSP& clazz){
    clazz->setPublic(); // make visible to EAS/spreadsheet
    REGISTER(FXVol, clazz);
    SUPERCLASS(FXVolBase);
    IMPLEMENTS(FXVega::IShift);
    IMPLEMENTS(FXVegaPointwise::IShift);
    IMPLEMENTS(IPDFCalculator);
    EMPTY_SHELL_METHOD(defaultFXVol);
    FIELD(name, "Vol identifier");
    FIELD(vol, "Vol");
}

CClassConstSP const FXVol::TYPE = CClass::registerClassLoadMethod(
    "FXVol", typeid(FXVol), FXVol::load);


DRLIB_END_NAMESPACE
