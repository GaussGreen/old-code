//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : FourierProcessGammaOU.cpp
//
//   Date        : 09 March 04
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/VolRequestRaw.hpp"
#include "edginc/FourierProcessGammaOU.hpp"
#include "edginc/FourierEngine.hpp"
#include "edginc/MDFAssetVol.hpp"

DRLIB_BEGIN_NAMESPACE

/* Started Log Return */
Complex FourierProcessGammaOU::scalelessCumulant(
    const StFourierProductLogRtn& product, 
    const Complex& z, 
    const DateTime& matDate) const{
    static const string method = "FourierProcessGammaOU::scalelessCumulant";
    try{
        return theVol->scalelessCumulant(static_cast<const StFourierProcessLogRtn&>(*this),
                                         product, 
                                         z, 
                                         matDate);
    }
    catch(exception& e){
        throw ModelException(e, method);
    }

}

double FourierProcessGammaOU::lowerRealBound(
    const StFourierProductLogRtn& product, 
    const DateTime& matDate) const{
    if (!stFrequencyRange){
        throw ModelException("FourierProcessGammaOU::lowerRealBound", 
                             "stFrequencyRange is missing");
    }
    return stFrequencyRange->lowerRealBound;
}

double FourierProcessGammaOU::upperRealBound(
    const StFourierProductLogRtn& product, 
    const DateTime& matDate) const{
    if (!stFrequencyRange){
        throw ModelException("FourierProcessGammaOU::lowerRealBound", 
                             "stFrequencyRange is missing");
    }
    return stFrequencyRange->upperRealBound;
}

/* Fwd starting Log Return  */
Complex FourierProcessGammaOU::scalelessCumulant(
    const FwdStFourierProductLogRtn& product, 
    const Complex& z, 
    const DateTime& matDate) const{
    static const string method = "FourierProcessGammaOU::scalelessCumulant";
    try{
        return theVol->scalelessCumulant(
            static_cast<const FwdStFourierProcessLogRtn&>(*this),
            product, 
            z, 
            matDate);
    }
    catch(exception& e){
        throw ModelException(e, method);
    }

}

double FourierProcessGammaOU::lowerRealBound(
    const FwdStFourierProductLogRtn& product, 
    const DateTime& matDate) const{
    if (!fwdStFrequencyRange){
        throw ModelException("FourierProcessGammaOU::lowerRealBound", 
                             "fwdStFrequencyRange is missing");
    }
    return fwdStFrequencyRange->lowerRealBound;
}

double FourierProcessGammaOU::upperRealBound(
    const FwdStFourierProductLogRtn& product, 
    const DateTime& matDate) const{
    if (!fwdStFrequencyRange){
        throw ModelException("FourierProcessGammaOU::lowerRealBound", 
                             "fwdStFrequencyRange is missing");
    }
    return fwdStFrequencyRange->upperRealBound;
}

/* Started integrated variance */
Complex FourierProcessGammaOU::cumulant(const StFourierProductIntVar& product, 
                                        const Complex z, 
                                        const DateTime& matDate) const {
    
    static const string method = "FourierProcessGammaOU::cumulant";
    try{
        return theVol->cumulant(
            static_cast<const StFourierProcessIntVar&>(*this),
            product, 
            z, 
            matDate);
    }
    catch(exception& e){
        throw ModelException(e, method);
    }
}

double FourierProcessGammaOU::lowerRealBound(
    const StFourierProductIntVar& product, 
    const DateTime& matDate) const{
    if (!stFrequencyRangeIntVar){
        throw ModelException("FourierProcessGammaOU::lowerRealBound", 
                             "StFrequencyRangeIntVar is missing");
    }
    return stFrequencyRangeIntVar->lowerRealBound;
}

double FourierProcessGammaOU::upperRealBound(
    const StFourierProductIntVar& product, 
    const DateTime& matDate) const{
    if (!stFrequencyRangeIntVar){
        throw ModelException("FourierProcessGammaOU::lowerRealBound", 
                             "StFrequencyRangeIntVar is missing");
    }
    return stFrequencyRangeIntVar->upperRealBound;
}
    
/* Forward Starting integrated variance */
Complex FourierProcessGammaOU::cumulant(
    const FwdStFourierProductIntVar& product, 
    const Complex z, 
    const DateTime& matDate) const {        
    static const string method = "FourierProcessGammaOU::cumulant";
    try{
        return theVol->cumulant(
            static_cast<const FwdStFourierProcessIntVar&>(*this),
            product,
            z, 
            matDate);
    }
    catch(exception& e){
        throw ModelException(e, method);
    }
}    

double FourierProcessGammaOU::lowerRealBound(
    const FwdStFourierProductIntVar& product, 
    const DateTime& matDate) const{
    if (!fwdStFrequencyRangeIntVar){
        throw ModelException("FourierProcessGammaOU::lowerRealBound", 
                             "fwdStFrequencyRangeIntVar is missing");
    }
    return fwdStFrequencyRangeIntVar->lowerRealBound;
}

double FourierProcessGammaOU::upperRealBound(
    const FwdStFourierProductIntVar& product, 
    const DateTime& matDate) const{
    if (!fwdStFrequencyRangeIntVar){
        throw ModelException("FourierProcessGammaOU::lowerRealBound", 
                             "fwdStFrequencyRangeIntVar is missing");
    }
    return fwdStFrequencyRangeIntVar->upperRealBound;
}

/** Create a MarketDataFetcher which will be used by the [Fourier] model
 * for retrieving market data etc */
MarketDataFetcherSP FourierProcessGammaOU::marketDataFetcher() const {
    return MarketDataFetcherSP(new MDFAssetVol(VolGammaOU::TYPE->getName(), 
                                               VolSurface::TYPE->getName()));
}

const TimeMetric& FourierProcessGammaOU::getTimeMetric() const{
    return *timeMetric;
}

string FourierProcessGammaOU::getParameters() const {
    return FourierProcess::extractParameters(theVol.get());
}

void FourierProcessGammaOU::validate(const FourierProduct* product){
    static const string method = "FourierProcessGammaOU::validate";
    try{
        const IMultiFactors& mAsset = product->getMultiAsset();
        // validate size
        if (mAsset.NbAssets() != 1){
            throw ModelException(method, "only 1 asset suppported");
        }
        // copy the vol that's in the multiasset
        theVol = VolGammaOUSP::dynamicCast(
            VolRequestRaw::copyVolBase(mAsset, 0));
        // copy time metric
        timeMetric = TimeMetricSP(copy(&theVol->getTimeMetric()));
    }
    catch (exception& e){
        throw ModelException(e, method);
    }
}

FourierProcessGammaOU::FourierProcessGammaOU(): 
FourierProcessParametric(TYPE){}

IObject* FourierProcessGammaOU::defaultCtor(){
    return new FourierProcessGammaOU();
}

/** Invoked when Class is 'loaded' */
void FourierProcessGammaOU::load(CClassSP& clazz){
    clazz->setPublic(); // make visible to EAS/spreadsheet
    REGISTER(FourierProcessGammaOU, clazz);
    SUPERCLASS(FourierProcessParametric);
    IMPLEMENTS(StFourierProcessLogRtn);
    IMPLEMENTS(FwdStFourierProcessLogRtn);
    IMPLEMENTS(FwdStFourierProcessIntVar);
    IMPLEMENTS(StFourierProcessIntVar);
    EMPTY_SHELL_METHOD(defaultCtor);
    
    // Optional
    FIELD(stFrequencyRange, "stFrequencyRange");
    FIELD_MAKE_OPTIONAL(stFrequencyRange);
    FIELD(fwdStFrequencyRange, "fwdStFrequencyRange");
    FIELD_MAKE_OPTIONAL(fwdStFrequencyRange);
    FIELD(fwdStFrequencyRangeIntVar, "fwdStFrequencyRangeIntVar");
    FIELD_MAKE_OPTIONAL(fwdStFrequencyRangeIntVar);
    FIELD(stFrequencyRangeIntVar, "stFrequencyRangeIntVar");
    FIELD_MAKE_OPTIONAL(stFrequencyRangeIntVar);

    // transient
    FIELD(theVol, ""); 
    FIELD_MAKE_TRANSIENT_BUT_TWEAKABLE(theVol);
    FIELD(timeMetric, "");
    FIELD_MAKE_TRANSIENT(timeMetric);
}
 
CClassConstSP const FourierProcessGammaOU::TYPE =
CClass::registerClassLoadMethod("FourierProcessGammaOU", typeid(FourierProcessGammaOU), load);

/* external symbol to allow class to be forced to be linked in */
bool FourierProcessGammaOULoad(){
    return (FourierProcessGammaOU::TYPE != 0);
}

DRLIB_END_NAMESPACE
