//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : FourierProcessAJDSuper.cpp
//
//   Date        : 10 March 2004
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/VolRequestRaw.hpp"
#include "edginc/FourierProcessAJDSuper.hpp"
#include "edginc/FourierEngine.hpp"
#include "edginc/MDFAssetVol.hpp"

DRLIB_BEGIN_NAMESPACE

/* Started Log Return */
Complex FourierProcessAJDSuper::scalelessCumulant(
    const StFourierProductLogRtn& product, 
    const Complex& z, 
    const DateTime& matDate) const{
    return theVol->scalelessCumulant(
        static_cast<const StFourierProcessLogRtn&>(*this), 
        product, 
        z, matDate);
}

double FourierProcessAJDSuper::lowerRealBound(
    const StFourierProductLogRtn& product, 
    const DateTime& matDate) const{
    if (!stFrequencyRange){
        throw ModelException("FourierProcessAJDSuper::lowerRealBound", 
                             "stFrequencyRange is missing");
    }
    return stFrequencyRange->lowerRealBound;
}

double FourierProcessAJDSuper::upperRealBound(
    const StFourierProductLogRtn& product, 
    const DateTime& matDate) const{
    if (!stFrequencyRange){
        throw ModelException("FourierProcessAJDSuper::lowerRealBound", 
                             "stFrequencyRange is missing");
    }
    return stFrequencyRange->upperRealBound;
}

/* Fwd starting Log Return  */
Complex FourierProcessAJDSuper::scalelessCumulant(
    const FwdStFourierProductLogRtn& product, 
    const Complex& z, 
    const DateTime& matDate) const{
    return theVol->scalelessCumulant(
        static_cast<const FwdStFourierProcessLogRtn&>(*this),
        product, 
        z, 
        matDate);
}

double FourierProcessAJDSuper::lowerRealBound(
    const FwdStFourierProductLogRtn& product, 
    const DateTime& matDate) const{
    if (!fwdStFrequencyRange){
        throw ModelException("FourierProcessAJDSuper::lowerRealBound", 
                             "fwdStFrequencyRange is missing");
    }
    return fwdStFrequencyRange->lowerRealBound;
}

double FourierProcessAJDSuper::upperRealBound(
    const FwdStFourierProductLogRtn& product, 
    const DateTime& matDate) const{
    if (!fwdStFrequencyRange){
        throw ModelException("FourierProcessAJDSuper::lowerRealBound", 
                             "fwdStFrequencyRange is missing");
    }
    return fwdStFrequencyRange->upperRealBound;
}

#if 0
Complex FourierProcessAJDSuper::cumulant(const StFourierProductIntVar& product, 
                                    const Complex z, 
                                    const DateTime& matDate) const {
    return theVol->cumulant(static_cast<const StFourierProcessIntVar&>(*this),
                            product, 
                            z, 
                            matDate);
}

double FourierProcessAJDSuper::lowerRealBound(
    const StFourierProductIntVar& product, 
    const DateTime& matDate) const{
    if (!stFrequencyRangeIntVar){
        throw ModelException("FourierProcessAJDSuper::lowerRealBound", 
                             "StFrequencyRangeIntVar is missing");
    }
    return stFrequencyRangeIntVar->lowerRealBound;
}

double FourierProcessAJDSuper::upperRealBound(
    const StFourierProductIntVar& product, 
    const DateTime& matDate) const{
    if (!stFrequencyRangeIntVar){
        throw ModelException("FourierProcessAJDSuper::lowerRealBound", 
                             "StFrequencyRangeIntVar is missing");
    }
    return stFrequencyRangeIntVar->upperRealBound;
}
    
Complex FourierProcessAJDSuper::cumulant(const FwdStFourierProductIntVar& product, 
                                    const Complex z, 
                                    const DateTime& matDate) const {        
    return theVol->cumulant(
        static_cast<const FwdStFourierProcessIntVar&>(*this),
        product,
        z, 
        matDate);
}    

double FourierProcessAJDSuper::lowerRealBound(
    const FwdStFourierProductIntVar& product, 
    const DateTime& matDate) const{
    if (!fwdStFrequencyRangeIntVar){
        throw ModelException("FourierProcessAJDSuper::lowerRealBound", 
                             "fwdStFrequencyRangeIntVar is missing");
    }
    return fwdStFrequencyRangeIntVar->lowerRealBound;
}

double FourierProcessAJDSuper::upperRealBound(
    const FwdStFourierProductIntVar& product, 
    const DateTime& matDate) const{
    if (!fwdStFrequencyRangeIntVar){
        throw ModelException("FourierProcessAJDSuper::lowerRealBound", 
                             "fwdStFrequencyRangeIntVar is missing");
    }
    return fwdStFrequencyRangeIntVar->upperRealBound;
}
   
/* Started Quadratic Variation */
Complex FourierProcessAJDSuper::cumulant(const StFourierProductQuadVar& product, 
                                    const Complex z, 
                                    const DateTime& matDate) const {
    return theVol->cumulant(static_cast<const StFourierProcessQuadVar&>(*this),
                            product, 
                            z, 
                            matDate);
}

double FourierProcessAJDSuper::lowerRealBound(
    const StFourierProductQuadVar& product, 
    const DateTime& matDate) const{
    if (!stFrequencyRangeIntVar){
        throw ModelException("FourierProcessAJDSuper::lowerRealBound", 
                             "StFrequencyRangeIntVar is missing");
    }
    return stFrequencyRangeIntVar->lowerRealBound;
}

double FourierProcessAJDSuper::upperRealBound(
    const StFourierProductQuadVar& product, 
    const DateTime& matDate) const{
    if (!stFrequencyRangeIntVar){
        throw ModelException("FourierProcessAJDSuper::lowerRealBound", 
                             "StFrequencyRangeIntVar is missing");
    }
    return stFrequencyRangeIntVar->upperRealBound;
}

double FourierProcessAJDSuper::expectation(const StFourierProductQuadVar& product, 
                                      const DateTime& matDate) const {
    return theVol->expectation(static_cast<const StFourierProcessQuadVar&>(*this),
                               product, 
                               matDate);
}
#endif

IModel::WantsRiskMapping FourierProcessAJDSuper::wantsRiskMapping() const {
    return IModel::riskMappingIrrelevant;
}

/** Create a MarketDataFetcher which will be used by the [Fourier] model
 * for retrieving market data etc */
MarketDataFetcherSP FourierProcessAJDSuper::marketDataFetcher() const {
    return MarketDataFetcherSP(new MDFAssetVol(VolAJDSuper::TYPE->getName(), 
                                               VolSurface::TYPE->getName()));
}


const TimeMetric& FourierProcessAJDSuper::getTimeMetric() const{
    return *timeMetric;
}

string FourierProcessAJDSuper::getParameters() const {
    return FourierProcess::extractParameters(theVol.get());
}

void FourierProcessAJDSuper::validate(const FourierProduct* product){
    static const string method = "FourierProcessAJDSuper::validate";
    try{
        const IMultiFactors& mAsset = product->getMultiAsset();
        // validate size
        if (mAsset.NbAssets() != 1){
            throw ModelException(method, "only 1 asset suppported");
        }
        // copy the vol that's in the multiasset
        theVol = VolAJDSuperSP::dynamicCast(
            VolRequestRaw::copyVolBase(mAsset, 0));
        // copy time metric
        timeMetric = TimeMetricSP(copy(&theVol->getTimeMetric()));
    }
    catch (exception& e){
        throw ModelException(e, method);
    }
}

FourierProcessAJDSuper::FourierProcessAJDSuper(): 
FourierProcess(TYPE){}

IObject* FourierProcessAJDSuper::defaultCtor(){
    return new FourierProcessAJDSuper();
}

/** Invoked when Class is 'loaded' */
void FourierProcessAJDSuper::load(CClassSP& clazz){
    clazz->setPublic(); // make visible to EAS/spreadsheet
    REGISTER(FourierProcessAJDSuper, clazz);
    SUPERCLASS(FourierProcess);
    IMPLEMENTS(StFourierProcessLogRtn);
    IMPLEMENTS(FwdStFourierProcessLogRtn);
//    IMPLEMENTS(FwdStFourierProcessIntVar);
//    IMPLEMENTS(StFourierProcessIntVar);
//    IMPLEMENTS(StFourierProcessQuadVar);
    EMPTY_SHELL_METHOD(defaultCtor);
    FIELD(stFrequencyRange, "stFrequencyRange");
    FIELD_MAKE_OPTIONAL(stFrequencyRange);
    FIELD(fwdStFrequencyRange, "fwdStFrequencyRange");
    FIELD_MAKE_OPTIONAL(fwdStFrequencyRange);
//    FIELD(fwdStFrequencyRangeIntVar, "fwdStFrequencyRangeIntVar");
//    FIELD_MAKE_OPTIONAL(fwdStFrequencyRangeIntVar);
//    FIELD(stFrequencyRangeIntVar, "stFrequencyRangeIntVar");
//    FIELD_MAKE_OPTIONAL(stFrequencyRangeIntVar);
    // transient
    FIELD(theVol, ""); 
    FIELD_MAKE_TRANSIENT_BUT_TWEAKABLE(theVol);
    FIELD(timeMetric, "");
    FIELD_MAKE_TRANSIENT(timeMetric);
}

CClassConstSP const FourierProcessAJDSuper::TYPE =
CClass::registerClassLoadMethod("FourierProcessAJDSuper", 
                                typeid(FourierProcessAJDSuper), load);

/* external symbol to allow class to be forced to be linked in */
bool FourierProcessAJDSuperLoad(){
    return (FourierProcessAJDSuper::TYPE != 0);
}

DRLIB_END_NAMESPACE
