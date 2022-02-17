//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : FourierProcessSVJJ.hpp
//
//   Date        : 09 March 2004
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/VolRequestRaw.hpp"
#include "edginc/FourierProcessSVJJ.hpp"
#include "edginc/FourierEngine.hpp"
#include "edginc/MDFAssetVol.hpp"

DRLIB_BEGIN_NAMESPACE

/* Started Log Return */
Complex FourierProcessSVJJ::scalelessCumulant(const StFourierProductLogRtn& product, 
                                              const Complex& z, 
                                              const DateTime& matDate) const{
    return theVol->scalelessCumulant(static_cast<const StFourierProcessLogRtn&>(*this),
                                     product, 
                                     z, 
                                     matDate);
}

double FourierProcessSVJJ::lowerRealBound(const StFourierProductLogRtn& product, 
                                          const DateTime& matDate) const{
    if (!stFrequencyRange){
        throw ModelException("FourierProcessSVJJ::lowerRealBound", 
                             "stFrequencyRange is missing");
    }
    return stFrequencyRange->lowerRealBound;
}

double FourierProcessSVJJ::upperRealBound(const StFourierProductLogRtn& product, 
                                          const DateTime& matDate) const{
    if (!stFrequencyRange){
        throw ModelException("FourierProcessSVJJ::lowerRealBound", 
                             "stFrequencyRange is missing");
    }
    return stFrequencyRange->upperRealBound;
}

/* Fwd starting Log Return  */
Complex FourierProcessSVJJ::scalelessCumulant(const FwdStFourierProductLogRtn& product, 
                                              const Complex& z, 
                                              const DateTime& matDate) const{
    return theVol->scalelessCumulant(static_cast<const FwdStFourierProcessLogRtn&>(*this),
                                     product, 
                                     z, 
                                     matDate);
}

double FourierProcessSVJJ::lowerRealBound(const FwdStFourierProductLogRtn& product, 
                                          const DateTime& matDate) const{
    if (!fwdStFrequencyRange){
        throw ModelException("FourierProcessSVJJ::lowerRealBound", 
                             "fwdStFrequencyRange is missing");
    }
    return fwdStFrequencyRange->lowerRealBound;
}

double FourierProcessSVJJ::upperRealBound(const FwdStFourierProductLogRtn& product, 
                                          const DateTime& matDate) const{
    if (!fwdStFrequencyRange){
        throw ModelException("FourierProcessSVJJ::lowerRealBound", 
                             "fwdStFrequencyRange is missing");
    }
    return fwdStFrequencyRange->upperRealBound;
}

/* Started integrated variance */
Complex FourierProcessSVJJ::cumulant(const StFourierProductIntVar& product, 
                                     const Complex z, 
                                     const DateTime& matDate) const {
    return theVol->cumulant(static_cast<const StFourierProcessIntVar&>(*this),
                            product, 
                            z, 
                            matDate);
}

double FourierProcessSVJJ::lowerRealBound(const StFourierProductIntVar& product, 
                                          const DateTime& matDate) const{
    if (!stFrequencyRangeIntVar){
        throw ModelException("FourierProcessSVJJ::lowerRealBound", 
                             "StFrequencyRangeIntVar is missing");
    }
    return stFrequencyRangeIntVar->lowerRealBound;
}

double FourierProcessSVJJ::upperRealBound(const StFourierProductIntVar& product, 
                                          const DateTime& matDate) const{
    if (!stFrequencyRangeIntVar){
        throw ModelException("FourierProcessSVJJ::lowerRealBound", 
                             "StFrequencyRangeIntVar is missing");
    }
    return stFrequencyRangeIntVar->upperRealBound;
}
    
/* Forward Starting integrated variance */
Complex FourierProcessSVJJ::cumulant(const FwdStFourierProductIntVar& product, 
                                     const Complex z, 
                                     const DateTime& matDate) const {        
    return theVol->cumulant(static_cast<const FwdStFourierProcessIntVar&>(*this),
                            product,
                            z, 
                            matDate);
}    

double FourierProcessSVJJ::lowerRealBound(const FwdStFourierProductIntVar& product, 
                                          const DateTime& matDate) const{
    if (!fwdStFrequencyRangeIntVar){
        throw ModelException("FourierProcessSVJJ::lowerRealBound", 
                             "fwdStFrequencyRangeIntVar is missing");
    }
    return fwdStFrequencyRangeIntVar->lowerRealBound;
}

double FourierProcessSVJJ::upperRealBound(const FwdStFourierProductIntVar& product, 
                                          const DateTime& matDate) const{
    if (!fwdStFrequencyRangeIntVar){
        throw ModelException("FourierProcessSVJJ::lowerRealBound", 
                             "fwdStFrequencyRangeIntVar is missing");
    }
    return fwdStFrequencyRangeIntVar->upperRealBound;
}
   
/** Create a MarketDataFetcher which will be used by the [Fourier] model
 * for retrieving market data etc */
MarketDataFetcherSP FourierProcessSVJJ::marketDataFetcher() const {
    return MarketDataFetcherSP(new MDFAssetVol(VolSVJJ::TYPE->getName(), 
                                               VolSurface::TYPE->getName()));
}

const TimeMetric& FourierProcessSVJJ::getTimeMetric() const{
    return *timeMetric;
}

string FourierProcessSVJJ::getParameters() const {
    return FourierProcess::extractParameters(theVol.get());
}

void FourierProcessSVJJ::validate(const FourierProduct* product){
    static const string method = "FourierProcessSV::validate";
    try{
        const IMultiFactors& mAsset = product->getMultiAsset();
        // validate size
        if (mAsset.NbAssets() != 1){
            throw ModelException(method, "only 1 asset suppported");
        }
        // copy the vol that's in the multiasset
        theVol = VolSVJJSP::dynamicCast(VolRequestRaw::copyVolBase(mAsset, 0));
        // copy time metric
        timeMetric = TimeMetricSP(copy(&theVol->getTimeMetric()));
    }
    catch (exception& e){
        throw ModelException(e, method);
    }
}

FourierProcessSVJJ::FourierProcessSVJJ(): 
FourierProcessParametric(TYPE){}

IObject* FourierProcessSVJJ::defaultCtor(){
    return new FourierProcessSVJJ();
}

/** Invoked when Class is 'loaded' */
void FourierProcessSVJJ::load(CClassSP& clazz){
    clazz->setPublic(); // make visible to EAS/spreadsheet
    REGISTER(FourierProcessSVJJ, clazz);
    SUPERCLASS(FourierProcessParametric);
    IMPLEMENTS(StFourierProcessLogRtn);
    IMPLEMENTS(FwdStFourierProcessLogRtn);
    IMPLEMENTS(FwdStFourierProcessIntVar);
    IMPLEMENTS(StFourierProcessIntVar);
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

CClassConstSP const FourierProcessSVJJ::TYPE =
CClass::registerClassLoadMethod("FourierProcessSVJJ", typeid(FourierProcessSVJJ), load);

/* external symbol to allow class to be forced to be linked in */
bool FourierProcessSVJJLoad(){
    return (FourierProcessSVJJ::TYPE != 0);
}

DRLIB_END_NAMESPACE
