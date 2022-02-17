//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : FourierProcessMerton.cpp
//
//   Date        : 10 March 04
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/VolRequestRaw.hpp"
#include "edginc/FourierProcessMerton.hpp"
#include "edginc/FourierEngine.hpp"
#include "edginc/MDFAssetVol.hpp"

DRLIB_BEGIN_NAMESPACE

/* Started Log Return */
Complex FourierProcessMerton::scalelessCumulant(
    const StFourierProductLogRtn& product, 
    const Complex& z, 
    const DateTime& matDate) const{
    return theVol->scalelessCumulant(
        static_cast<const StFourierProcessLogRtn&>(*this), 
        product, 
        z, matDate);
}

double FourierProcessMerton::lowerRealBound(
    const StFourierProductLogRtn& product, 
    const DateTime& matDate) const{
    if (!stFrequencyRange){
        throw ModelException("FourierProcessMerton::lowerRealBound", 
                             "stFrequencyRange is missing");
    }
    return stFrequencyRange->lowerRealBound;
}

double FourierProcessMerton::upperRealBound(
    const StFourierProductLogRtn& product, 
    const DateTime& matDate) const{
    if (!stFrequencyRange){
        throw ModelException("FourierProcessMerton::lowerRealBound", 
                             "stFrequencyRange is missing");
    }
    return stFrequencyRange->upperRealBound;
}

/* Fwd starting Log Return  */
Complex FourierProcessMerton::scalelessCumulant(
    const FwdStFourierProductLogRtn& product, 
    const Complex& z, 
    const DateTime& matDate) const{
    return theVol->scalelessCumulant(
        static_cast<const FwdStFourierProcessLogRtn&>(*this),
        product, 
        z, 
        matDate);
}

double FourierProcessMerton::lowerRealBound(
    const FwdStFourierProductLogRtn& product, 
    const DateTime& matDate) const{
    if (!fwdStFrequencyRange){
        throw ModelException("FourierProcessMerton::lowerRealBound", 
                             "fwdStFrequencyRange is missing");
    }
    return fwdStFrequencyRange->lowerRealBound;
}

double FourierProcessMerton::upperRealBound(
    const FwdStFourierProductLogRtn& product, 
    const DateTime& matDate) const{
    if (!fwdStFrequencyRange){
        throw ModelException("FourierProcessMerton::lowerRealBound", 
                             "fwdStFrequencyRange is missing");
    }
    return fwdStFrequencyRange->upperRealBound;
}

/* Started integrated variance */
Complex FourierProcessMerton::cumulant(const StFourierProductIntVar& product, 
                                    const Complex z, 
                                    const DateTime& matDate) const {
    return theVol->cumulant(static_cast<const StFourierProcessIntVar&>(*this),
                            product, 
                            z, 
                            matDate);
}

double FourierProcessMerton::lowerRealBound(
    const StFourierProductIntVar& product, 
    const DateTime& matDate) const{
    if (!stFrequencyRangeIntVar){
        throw ModelException("FourierProcessMerton::lowerRealBound", 
                             "StFrequencyRangeIntVar is missing");
    }
    return stFrequencyRangeIntVar->lowerRealBound;
}

double FourierProcessMerton::upperRealBound(
    const StFourierProductIntVar& product, 
    const DateTime& matDate) const{
    if (!stFrequencyRangeIntVar){
        throw ModelException("FourierProcessMerton::lowerRealBound", 
                             "StFrequencyRangeIntVar is missing");
    }
    return stFrequencyRangeIntVar->upperRealBound;
}

Complex FourierProcessMerton::cumulant(const StFourierProductQuadVar& product, 
                                    const Complex z, 
                                    const DateTime& matDate) const {
    return theVol->cumulant(static_cast<const StFourierProcessQuadVar&>(*this),
                            product, 
                            z, 
                            matDate);
}

double FourierProcessMerton::lowerRealBound(
    const StFourierProductQuadVar& product, 
    const DateTime& matDate) const{
    if (!stFrequencyRangeIntVar){
        throw ModelException("FourierProcessMerton::lowerRealBound", 
                             "StFrequencyRangeIntVar is missing");
    }
    return stFrequencyRangeIntVar->lowerRealBound;
}

double FourierProcessMerton::upperRealBound(
    const StFourierProductQuadVar& product, 
    const DateTime& matDate) const{
    if (!stFrequencyRangeIntVar){
        throw ModelException("FourierProcessMerton::lowerRealBound", 
                             "StFrequencyRangeIntVar is missing");
    }
    return stFrequencyRangeIntVar->upperRealBound;
}

Complex FourierProcessMerton::cumulant(const FwdStFourierProductIntVar& product, 
                                    const Complex z, 
                                    const DateTime& matDate) const {        
    return theVol->cumulant(
        static_cast<const FwdStFourierProcessIntVar&>(*this),
        product,
        z, 
        matDate);
}    

double FourierProcessMerton::lowerRealBound(
    const FwdStFourierProductIntVar& product, 
    const DateTime& matDate) const{
    if (!fwdStFrequencyRangeIntVar){
        throw ModelException("FourierProcessMerton::lowerRealBound", 
                             "fwdStFrequencyRangeIntVar is missing");
    }
    return fwdStFrequencyRangeIntVar->lowerRealBound;
}

double FourierProcessMerton::upperRealBound(
    const FwdStFourierProductIntVar& product, 
    const DateTime& matDate) const{
    if (!fwdStFrequencyRangeIntVar){
        throw ModelException("FourierProcessMerton::lowerRealBound", 
                             "fwdStFrequencyRangeIntVar is missing");
    }
    return fwdStFrequencyRangeIntVar->upperRealBound;
} 

double FourierProcessMerton::expectation(const StFourierProductQuadVar& product, 
                                      const DateTime& matDate) const {
    return theVol->expectation(static_cast<const StFourierProcessQuadVar&>(*this),
                               product, 
                               matDate);
}

/** Create a MarketDataFetcher which will be used by the [Fourier] model
 * for retrieving market data etc */
MarketDataFetcherSP FourierProcessMerton::marketDataFetcher() const {
    return MarketDataFetcherSP(new MDFAssetVol(VolMerton::TYPE->getName(), 
                                               VolSurface::TYPE->getName()));
}

const TimeMetric& FourierProcessMerton::getTimeMetric() const{
    return *timeMetric;
}

string FourierProcessMerton::getParameters() const {
    return FourierProcess::extractParameters(theVol.get());
}

void FourierProcessMerton::validate(const FourierProduct* product){
    static const string method = "FourierProcessMerton::validate";
    try{
        const IMultiFactors& mAsset = product->getMultiAsset();
        // validate size
        if (mAsset.NbAssets() != 1){
            throw ModelException(method, "only 1 asset suppported");
        }
        // copy the vol that's in the multiasset
        theVol = VolMertonSP::dynamicCast(
            VolRequestRaw::copyVolBase(mAsset, 0));
        // copy time metric
        timeMetric = TimeMetricSP(copy(&theVol->getTimeMetric()));
    }
    catch (exception& e){
        throw ModelException(e, method);
    }
}

FourierProcessMerton::FourierProcessMerton(): 
FourierProcessParametric(TYPE){}

IObject* FourierProcessMerton::defaultCtor(){
    return new FourierProcessMerton();
}

/** Invoked when Class is 'loaded' */
void FourierProcessMerton::load(CClassSP& clazz){
    clazz->setPublic(); // make visible to EAS/spreadsheet
    REGISTER(FourierProcessMerton, clazz);
    SUPERCLASS(FourierProcessParametric);
    IMPLEMENTS(StFourierProcessLogRtn);
    IMPLEMENTS(FwdStFourierProcessLogRtn);
    IMPLEMENTS(FwdStFourierProcessIntVar);
    IMPLEMENTS(StFourierProcessIntVar);
    IMPLEMENTS(StFourierProcessQuadVar);
    EMPTY_SHELL_METHOD(defaultCtor);
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
 
CClassConstSP const FourierProcessMerton::TYPE =
CClass::registerClassLoadMethod("FourierProcessMerton", 
                                typeid(FourierProcessMerton), load);

/* external symbol to allow class to be forced to be linked in */
bool FourierProcessMertonLoad(){
    return (FourierProcessMerton::TYPE != 0);
}

DRLIB_END_NAMESPACE
