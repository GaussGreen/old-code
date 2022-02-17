//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : FourierProcessSVCJ.cpp
//
//   Date        : 27 April 2004
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/VolRequestRaw.hpp"
#include "edginc/FourierProcessSVCJ.hpp"
#include "edginc/FourierEngine.hpp"
#include "edginc/MDFAssetVol.hpp"

DRLIB_BEGIN_NAMESPACE

/* Started Log Return */
Complex FourierProcessSVCJ::scalelessCumulant(const StFourierProductLogRtn& product, 
                                              const Complex& z, 
                                              const DateTime& matDate) const{
    return theVol->scalelessCumulant(static_cast<const StFourierProcessLogRtn&>(*this),
                                     product, 
                                     z, 
                                     matDate);
}

double FourierProcessSVCJ::lowerRealBound(const StFourierProductLogRtn& product, 
                                          const DateTime& matDate) const{
    if (!stFrequencyRange){
        throw ModelException("FourierProcessSVCJ::lowerRealBound", 
                             "stFrequencyRange is missing");
    }
    return stFrequencyRange->lowerRealBound;
}

double FourierProcessSVCJ::upperRealBound(const StFourierProductLogRtn& product, 
                                          const DateTime& matDate) const{
    if (!stFrequencyRange){
        throw ModelException("FourierProcessSVCJ::lowerRealBound", 
                             "stFrequencyRange is missing");
    }
    return stFrequencyRange->upperRealBound;
}

/* Fwd starting Log Return  */
Complex FourierProcessSVCJ::scalelessCumulant(const FwdStFourierProductLogRtn& product, 
                                              const Complex& z, 
                                              const DateTime& matDate) const{
    return theVol->scalelessCumulant(static_cast<const FwdStFourierProcessLogRtn&>(*this),
                                     product, 
                                     z, 
                                     matDate);
}

double FourierProcessSVCJ::lowerRealBound(const FwdStFourierProductLogRtn& product, 
                                          const DateTime& matDate) const{
    if (!fwdStFrequencyRange){
        throw ModelException("FourierProcessSVCJ::lowerRealBound", 
                             "fwdStFrequencyRange is missing");
    }
    return fwdStFrequencyRange->lowerRealBound;
}

double FourierProcessSVCJ::upperRealBound(const FwdStFourierProductLogRtn& product, 
                                          const DateTime& matDate) const{
    if (!fwdStFrequencyRange){
        throw ModelException("FourierProcessSVCJ::lowerRealBound", 
                             "fwdStFrequencyRange is missing");
    }
    return fwdStFrequencyRange->upperRealBound;
}

/* Started integrated variance */
Complex FourierProcessSVCJ::cumulant(const StFourierProductIntVar&  product, 
                                     const Complex                  z, 
                                     const DateTime&                matDate) const {
    return theVol->cumulant(static_cast<const StFourierProcessIntVar&>(*this),
                            product, 
                            z, 
                            matDate);
}

double FourierProcessSVCJ::lowerRealBound(const StFourierProductIntVar& product, 
                                          const DateTime& matDate) const{
    if (!stFrequencyRangeIntVar){
        throw ModelException("FourierProcessSVCJ::lowerRealBound", 
                             "StFrequencyRangeIntVar is missing");
    }
    return stFrequencyRangeIntVar->lowerRealBound;
}

double FourierProcessSVCJ::upperRealBound(const StFourierProductIntVar& product, 
                                          const DateTime&               matDate) const{
    if (!stFrequencyRangeIntVar){
        throw ModelException("FourierProcessSVCJ::lowerRealBound", 
                             "StFrequencyRangeIntVar is missing");
    }
    return stFrequencyRangeIntVar->upperRealBound;
}
    
/* Forward Starting integrated variance */
Complex FourierProcessSVCJ::cumulant(const FwdStFourierProductIntVar& product, 
                                     const Complex                    z, 
                                     const DateTime&                  matDate) const {        
    return theVol->cumulant(static_cast<const FwdStFourierProcessIntVar&>(*this),
                            product,
                            z, 
                            matDate);
}    

double FourierProcessSVCJ::lowerRealBound(const FwdStFourierProductIntVar&  product, 
                                          const DateTime&                   matDate) const{
    if (!fwdStFrequencyRangeIntVar){
        throw ModelException("FourierProcessSVCJ::lowerRealBound", 
                             "fwdStFrequencyRangeIntVar is missing");
    }
    return fwdStFrequencyRangeIntVar->lowerRealBound;
}

double FourierProcessSVCJ::upperRealBound(const FwdStFourierProductIntVar&  product, 
                                          const DateTime&                   matDate) const{
    if (!fwdStFrequencyRangeIntVar){
        throw ModelException("FourierProcessSVCJ::lowerRealBound", 
                             "fwdStFrequencyRangeIntVar is missing");
    }
    return fwdStFrequencyRangeIntVar->upperRealBound;
}

/* Started quadratic variation */
Complex FourierProcessSVCJ::cumulant(const StFourierProductQuadVar& product, 
                                     const Complex                  z,
                                     const DateTime&                matDate) const {
    return theVol->cumulant(static_cast<const StFourierProcessQuadVar&>(*this),
                            product,
                            z, 
                            matDate);
}

double FourierProcessSVCJ::lowerRealBound(const StFourierProductQuadVar&    product, 
                                          const DateTime&                   matDate) const{
    if (!stFrequencyRangeIntVar){
        throw ModelException("FourierProcessSVCJ::lowerRealBound", 
                             "StFrequencyRangeIntVar is missing");
    }
    return stFrequencyRangeIntVar->lowerRealBound;
}

double FourierProcessSVCJ::upperRealBound(const StFourierProductQuadVar&    product, 
                                          const DateTime&                   matDate) const{
    if (!stFrequencyRangeIntVar){
        throw ModelException("FourierProcessSVCJ::lowerRealBound", 
                             "StFrequencyRangeIntVar is missing");
    }
    return stFrequencyRangeIntVar->upperRealBound;
}

double FourierProcessSVCJ::expectation(const StFourierProductQuadVar& product, 
                                       const DateTime& matDate) const {
    return theVol->expectation(static_cast<const StFourierProcessQuadVar&>(*this),
                               product, 
                               matDate);
}
    
/* Forward Starting quadratic variation */
Complex FourierProcessSVCJ::cumulant(const FwdStFourierProductQuadVar&  product, 
                                     const Complex                      z, 
                                     const DateTime&                    matDate) const {        
    return theVol->cumulant(static_cast<const FwdStFourierProcessQuadVar&>(*this),
                            product,
                            z, 
                            matDate);
}    

double FourierProcessSVCJ::lowerRealBound(const FwdStFourierProductQuadVar& product, 
                                          const DateTime&                   matDate) const{
    if (!fwdStFrequencyRangeIntVar){
        throw ModelException("FourierProcessSVCJ::lowerRealBound", 
                             "fwdStFrequencyRangeIntVar is missing");
    }
    return fwdStFrequencyRangeIntVar->lowerRealBound;
}

double FourierProcessSVCJ::upperRealBound(const FwdStFourierProductQuadVar& product, 
                                          const DateTime&                   matDate) const{
    if (!fwdStFrequencyRangeIntVar){
        throw ModelException("FourierProcessSVCJ::lowerRealBound", 
                             "fwdStFrequencyRangeIntVar is missing");
    }
    return fwdStFrequencyRangeIntVar->upperRealBound;
}

double FourierProcessSVCJ::expectation(const FwdStFourierProductQuadVar& product, 
                                       const DateTime& matDate) const{
    return theVol->expectation(static_cast<const FwdStFourierProcessQuadVar&>(*this),
                               product, 
                               matDate);
}

/** Create a MarketDataFetcher which will be used by the [Fourier] model
 * for retrieving market data etc */
MarketDataFetcherSP FourierProcessSVCJ::marketDataFetcher() const {
    return MarketDataFetcherSP(new MDFAssetVol(VolSVCJ::TYPE->getName(), 
                                               VolSurface::TYPE->getName()));
}


const TimeMetric& FourierProcessSVCJ::getTimeMetric() const{
    return *timeMetric;
}

string FourierProcessSVCJ::getParameters() const {
    return FourierProcess::extractParameters(theVol.get());
}

void FourierProcessSVCJ::validate(const FourierProduct* product){
    static const string method = "FourierProcessSV::validate";
    try{
        const IMultiFactors& mAsset = product->getMultiAsset();
        // validate size
        if (mAsset.NbAssets() != 1){
            throw ModelException(method, "only 1 asset suppported");
        }
        // copy the vol that's in the multiasset
        theVol = VolSVCJSP::dynamicCast(VolRequestRaw::copyVolBase(mAsset, 0));
        // copy time metric
        timeMetric = TimeMetricSP(copy(&theVol->getTimeMetric()));
    }
    catch (exception& e){
        throw ModelException(e, method);
    }
}

FourierProcessSVCJ::FourierProcessSVCJ(): 
FourierProcessParametric(TYPE){}

IObject* FourierProcessSVCJ::defaultCtor(){
    return new FourierProcessSVCJ();
}

/** Invoked when Class is 'loaded' */
void FourierProcessSVCJ::load(CClassSP& clazz){
    clazz->setPublic(); // make visible to EAS/spreadsheet
    REGISTER(FourierProcessSVCJ, clazz);
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

CClassConstSP const FourierProcessSVCJ::TYPE =
CClass::registerClassLoadMethod("FourierProcessSVCJ", typeid(FourierProcessSVCJ), load);

/* external symbol to allow class to be forced to be linked in */
bool FourierProcessSVCJLoad(){
    return (FourierProcessSVCJ::TYPE != 0);
}

DRLIB_END_NAMESPACE
