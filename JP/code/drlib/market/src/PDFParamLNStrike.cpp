//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : PDFDefaultLNStrike.cpp
//
//   Description : Implementation of PDFCalculator for parameterised vols
//
//   Author      : Mark A Robson
//
//   Date        : 24 March 2002
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/PDFParamLNStrike.hpp"
#include "edginc/Black.hpp"

DRLIB_BEGIN_NAMESPACE

/** Caclulates spreads using specified lo and hi strikes. Overrides
    default implementation  in PDFDefaultLNStrike to use the 
    ComputeImpVol on the VolParam */
void PDFParamLNStrike::spreads(const CLatticeDouble&    loStrikes,
                               const CLatticeDouble&    hiStrikes,
                               const DateTimeArray&     maturities,
                               const DoubleArray&       fwds, // at maturities
                               CLatticeDouble&          spread) const {
    static const string method = "PDFParamLNStrike::probabilities";
    try{
        // calculate vols at lo and hi strikes
        CLatticeDouble loVol(loStrikes.sizes());
        CLatticeDouble hiVol(hiStrikes.sizes());
        DateTime startDate = volFwdStarting? 
            lnRequest->getVolRequest()->getStartDate(): valueDate;
        /* this class is horribly confused about which fwd start dates
           are which - so avoid 'spotAtStart' field */
        double mySpotAtStart = asset->fwdValue(startDate);
        if (volFwdStarting){
            // capture fwd starting data in required object
            CVolParam::FwdStart fwdStart(valueDate, startDate, metric,
                                         asset->getSpot(), 
                                         mySpotAtStart);
            paramVol->computeFwdStartImpVol(vol.get(), fwdStart, loStrikes,
                                            true, // strikes are %
                                            maturities, loVol);
            paramVol->computeFwdStartImpVol(vol.get(), fwdStart, hiStrikes,
                                            true, // strikes are %
                                            maturities, hiVol);
        } else {
            paramVol->ComputeImpVol(vol.get(), loStrikes, maturities, loVol);
            paramVol->ComputeImpVol(vol.get(), hiStrikes, maturities, hiVol);
        }
        for(int iStep = 0; iStep < maturities.size(); iStep++) {
            double tradTime = metric->yearFrac(startDate, maturities[iStep]);
            for (int i = 0; i < loStrikes[iStep].size(); i++) {
                // turn vol into variances
                double loVar = loVol[iStep][i] * loVol[iStep][i] * tradTime;
                double hiVar = hiVol[iStep][i] * hiVol[iStep][i] * tradTime;
                // calculate prices
                double absStrike = volFwdStarting? 
                    loStrikes[iStep][i] * mySpotAtStart: loStrikes[iStep][i];
                double lo = Black::price(true, fwds[iStep], 
                                         absStrike, 1.0, loVar);
                absStrike = volFwdStarting? 
                    hiStrikes[iStep][i] * mySpotAtStart: hiStrikes[iStep][i]; 
                double hi = Black::price(true, fwds[iStep], 
                                         absStrike, 1.0, hiVar);
                // calculate 'spread'
                spread[iStep][i] = lo-hi;
            }
        }
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }
}

PDFParamLNStrike::PDFParamLNStrike(
    const DateTime&                  valueDate,
    const CAssetConstSP&             asset,
    const TimeMetricConstSP&         metric,
    const CVolBaseConstSP&           vol,
    const CVolParamConstSP&          paramVol,
    const PDFRequestLNStrikeConstSP& lnRequest):
    PDFDefaultLNStrike(TYPE, valueDate, asset, lnRequest),
    metric(metric), vol(vol), paramVol(paramVol){
}

PDFParamLNStrike::~PDFParamLNStrike(){}

PDFParamLNStrike::PDFParamLNStrike():PDFDefaultLNStrike(TYPE) {}

class PDFParamLNStrike::Helper{
public:
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        REGISTER(PDFParamLNStrike, clazz);
        SUPERCLASS(PDFDefaultLNStrike);
        EMPTY_SHELL_METHOD(defaultPDFParamLNStrike);
        FIELD(metric, "metric");
        FIELD(vol, "vol");
        FIELD(paramVol, "parameterised vol");
    }
    
    static IObject* defaultPDFParamLNStrike(){
        return new PDFParamLNStrike();
    }    
};

CClassConstSP const PDFParamLNStrike::TYPE = CClass::registerClassLoadMethod(
    "PDFParamLNStrike", typeid(PDFParamLNStrike), 
    PDFParamLNStrike::Helper::load);

DRLIB_END_NAMESPACE
