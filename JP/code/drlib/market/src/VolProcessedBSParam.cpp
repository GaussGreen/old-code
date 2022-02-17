//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : VolProcessedBSParam.cpp
//
//   Description : Processed BS parameterised vols. 
//                 Previously, was part of VolParam.hpp
//
//   Author      : Regis Guichard
//
//   Date        : 02 Mai 2001
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/VolProcessedBSParam.hpp"

DRLIB_BEGIN_NAMESPACE




// CVolProcessedBSParam stuff ------------------------------------------------
void CVolProcessedBSParam::load(CClassSP& clazz){
    REGISTER(CVolProcessedBSParam, clazz);
    SUPERCLASS(CVolProcessedBS);
    FIELD(myParamVol, "Parameterised vol");
    FIELD(myVolRequestLN, "Vol Request");
    FIELD(myAsset, "Asset");
    FIELD(myFXAsset, "FX Asset");
    FIELD(myEqFXCorr, "Asset-FX Correlation");
    FIELD(myProcessedSurface, "Cached process vol");
    FIELD_MAKE_TRANSIENT(myProcessedSurface);
    FIELD(myVolSurf, "Cached surface");
    FIELD_MAKE_TRANSIENT(myVolSurf);
}

CClassConstSP const CVolProcessedBSParam::TYPE = 
CClass::registerClassLoadMethod(
    "VolProcessedBSParam", typeid(CVolProcessedBSParam), load);

CVolProcessedBSParam::CVolProcessedBSParam(
    const CClassConstSP& clazz, 
    const CVolBase*         vol,
    const CVolParamConstSP& volParam,
    const CVolRequestLN* volRequest, const CAsset* asset,
    const FXAsset*   fxAsset, const  Correlation* eqFXCorr):
    CVolProcessedBS(clazz),
    myVol(CVolBaseConstSP::attachToRef(vol)), myParamVol(volParam), 
    myVolRequestLN(CVolRequestLNConstSP(copyIfRef(volRequest))),
    myAsset(AssetConstSP(copyIfRef(asset))),
    myFXAsset(FXAssetConstSP(copyIfRef(fxAsset))),
    myEqFXCorr(CorrelationConstSP(copyIfRef(eqFXCorr))) {}

CVolProcessedBSParam::CVolProcessedBSParam(
    const CVolBase*         vol,
    const CVolParamConstSP& volParam,
    const CVolRequestLN*    volRequest,
    const CAsset*           asset):
    CVolProcessedBS(TYPE),
    myVol(CVolBaseConstSP::attachToRef(vol)), myParamVol(volParam),
    myVolRequestLN(CVolRequestLNConstSP(copyIfRef(volRequest))), 
    myAsset(AssetConstSP(copyIfRef(asset))) {}

CVolProcessedBSParam::CVolProcessedBSParam(
    const CVolBase*         vol,
    const CVolParamConstSP& volParam,
    const CVolRequestLN*    volRequest,
    const CAsset*           asset,
    const FXAsset*          fxAsset, 
    const Correlation*      eqFXCorr):
    CVolProcessedBS(TYPE),
    myVol(CVolBaseConstSP::attachToRef(vol)), myParamVol(volParam),
    myVolRequestLN(CVolRequestLNConstSP(copyIfRef(volRequest))), 
    myAsset(AssetConstSP(copyIfRef(asset))),
    myFXAsset(FXAssetConstSP(copyIfRef(fxAsset))), 
    myEqFXCorr(CorrelationConstSP(copyIfRef(eqFXCorr))) {}

/** updates this VolProcessedBS with the given vol request and asset.
    Returns the updated VolProcessedBS - which may or may not be the
    same as 'this'. Use smart pointers to manage memory */
CVolProcessedBSParam* CVolProcessedBSParam::update(
    const CVolRequestLN*    volRequest, 
    const CAsset*           asset){
    if (myVolRequestLN.get() == volRequest &&  myAsset.get() == asset){
        // do nothing
    }
    else if (getRefCount() == 1){
        // can recycle old one as it is not in use
        myVolRequestLN = CVolRequestLNConstSP(copyIfRef(volRequest));
        myAsset = AssetConstSP(copyIfRef(asset));
    } else {
        // construct new one - old one is in use
        return new CVolProcessedBSParam(myVol.get(), myParamVol,
                                        volRequest, asset);
    }
    // clear cache
    myProcessedSurface.reset();
    myVolSurf.reset();
    return this;
}
                      

double CVolProcessedBSParam::CalcVar(const DateTime& date1,
                                     const DateTime& date2) const {
    processToSurface();
    return myProcessedSurface->CalcVar(date1, date2);
}

void CVolProcessedBSParam::CalcVar(const DateTimeArray& dateList,
                                   TCalcType      calcType, 
                                   CDoubleArray&  vols) const {
    processToSurface();
    myProcessedSurface->CalcVar(dateList,
                                calcType,
                                vols);
}

void CVolProcessedBSParam::CalcVar(const DateTime&      dateFrom,
                                   const DateTimeArray& datesTo,
                                   TCalcType      calcType, 
                                   CDoubleArray&  vols) const {
    processToSurface();
    myProcessedSurface->CalcVar(dateFrom,
                                datesTo,
                                calcType,
                                vols);
}

double CVolProcessedBSParam::CalcVol(const DateTime& date1,
                                     const DateTime& date2) const {
    processToSurface();
    return myProcessedSurface->CalcVol(date1, date2);
}

void CVolProcessedBSParam::CalcVol(const DateTimeArray& dateList,
                                   TCalcType      calcType, 
                                   CDoubleArray&  vols) const {
    processToSurface();
    myProcessedSurface->CalcVol(dateList,
                                calcType,
                                vols);
}

void CVolProcessedBSParam::CalcVol(const DateTime&      dateFrom,
                                   const DateTimeArray& datesTo,
                                   TCalcType      calcType, 
                                   CDoubleArray&  vols) const {
    processToSurface();
    myProcessedSurface->CalcVol(dateFrom,
                                datesTo,
                                calcType,
                                vols);
}

void CVolProcessedBSParam::populateCompositeVol(
    CompositeVol* compositeVol) const{
    // just delegate (just need to supply benchmark dates and vols)
    processToSurface();
    myProcessedSurface->populateCompositeVol(compositeVol);
}

/** identifies the market data name of the volatility */
string CVolProcessedBSParam::getName() const{
    processToSurface();
    return myVolSurf->getName();
}

/** calculates the trading time between two dates */
double CVolProcessedBSParam::calcTradingTime(const DateTime &date1, 
                                             const DateTime &date2) const{
    processToSurface();
    return myProcessedSurface->calcTradingTime(date1, date2);
}

/** retieve time measure for the vol */
TimeMetricConstSP CVolProcessedBSParam::GetTimeMetric() const{
    processToSurface();
    return myProcessedSurface->GetTimeMetric();
}

void CVolProcessedBSParam::processToSurface() const {
    if (!myProcessedSurface){
        if (!strikes){
              strikes = DoubleArraySP(new DoubleArray());
              strikes->reserve(1);
        } else {
            strikes->clear();
        }
        bool isStruck = myFXAsset.get() && myEqFXCorr.get();
        CVolRequestLNConstSP volRequestForStrikes = myVolRequestLN;
        double fxSpot = 1.0;
        if (isStruck){
            // need to scale strikes by 1/spot fx, but only non fwd start ones
            // some requests have mixture of fwd and non fwd start strikes
            fxSpot = myFXAsset->getSpot();
            CVolRequestLNSP volRequestCopy(volRequestForStrikes.clone());
            volRequestCopy->scale(1.0/fxSpot);
            volRequestForStrikes = volRequestCopy;
        }
        /* to do - make asset optional when doing vol (unless it's really 
           needed). Better idea: create dummy asset class which clients can
           pass down when asset is not really needed. Currently the vol
           surface generator addin passed null in for the asset */
        volRequestForStrikes->getSensitiveStrike(myAsset.get()? 
                                                 myAsset->getSpot()/fxSpot: 
                                                 -100.0,
                                                 strikes);

        myVolSurf = 
            VolSurfaceConstSP(myParamVol->
                              spotVolSurfaceFromStrikes(myVol.get(),*strikes));

        CVolProcessed* interpVol;
        if (myFXAsset.get() && myEqFXCorr.get()){
            interpVol = myVolSurf->getProcessedVol(myVolRequestLN.get(),
                                                   myAsset.get(),
                                                   myFXAsset.get(),
                                                   myEqFXCorr.get());
        } else {
            interpVol = myVolSurf->getProcessedVol(myVolRequestLN.get(),
                                                   myAsset.get());
        }
        myProcessedSurface = CVolProcessedBSConstSP(
            dynamic_cast<CVolProcessedBS*>(interpVol));
    }
}


DRLIB_END_NAMESPACE
