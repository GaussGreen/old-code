//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : FourierProcessSVJ.cpp
//
//   Date        : 09 March 2004
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/VolRequestRaw.hpp"
#include "edginc/FourierProcessSVJ.hpp"
#include "edginc/FourierEngine.hpp"
#include "edginc/MDFAssetVol.hpp"

DRLIB_BEGIN_NAMESPACE

/* Started Log Return */
Complex FourierProcessSVJ::scalelessCumulant(
    const StFourierProductLogRtn& product, 
    const Complex& z, 
    const DateTime& matDate) const{
    return theVol->scalelessCumulant(
        static_cast<const StFourierProcessLogRtn&>(*this), 
        product, 
        z, matDate);
}

double FourierProcessSVJ::lowerRealBound(
    const StFourierProductLogRtn& product, 
    const DateTime& matDate) const{
    if (!stFrequencyRange){
        throw ModelException("FourierProcessSVJ::lowerRealBound", 
                             "stFrequencyRange is missing");
    }
    return stFrequencyRange->lowerRealBound;
}

double FourierProcessSVJ::upperRealBound(
    const StFourierProductLogRtn& product, 
    const DateTime& matDate) const{
    if (!stFrequencyRange){
        throw ModelException("FourierProcessSVJ::lowerRealBound", 
                             "stFrequencyRange is missing");
    }
    return stFrequencyRange->upperRealBound;
}

/* Fwd starting Log Return  */
Complex FourierProcessSVJ::scalelessCumulant(
    const FwdStFourierProductLogRtn& product, 
    const Complex& z, 
    const DateTime& matDate) const{
    return theVol->scalelessCumulant(
        static_cast<const FwdStFourierProcessLogRtn&>(*this),
        product, 
        z, 
        matDate);
}

double FourierProcessSVJ::lowerRealBound(
    const FwdStFourierProductLogRtn& product, 
    const DateTime& matDate) const{
    if (!fwdStFrequencyRange){
        throw ModelException("FourierProcessSVJ::lowerRealBound", 
                             "fwdStFrequencyRange is missing");
    }
    return fwdStFrequencyRange->lowerRealBound;
}

double FourierProcessSVJ::upperRealBound(
    const FwdStFourierProductLogRtn& product, 
    const DateTime& matDate) const{
    if (!fwdStFrequencyRange){
        throw ModelException("FourierProcessSVJ::lowerRealBound", 
                             "fwdStFrequencyRange is missing");
    }
    return fwdStFrequencyRange->upperRealBound;
}


/**ARNAUD**////////////////////////////////////////////////////////////
/* Fwd Starting Expected Quad Var */
Complex FourierProcessSVJ::cumulant(
    const FwdStFourierProductExpQuadVar& product, 
    const Complex& z, 
    const DateTime& matDate) const{
    return theVol->cumulant(
        static_cast<const FwdStFourierProcessExpQuadVar&>(*this),
        product, 
        z, 
        matDate);
}

double FourierProcessSVJ::lowerRealBound(
    const FwdStFourierProductExpQuadVar& product, 
    const DateTime& matDate) const{
    if (!fwdStFrequencyRangeIntVar){
        throw ModelException("FourierProcessSVJ::lowerRealBound", 
                             "fwdStFrequencyRangeIntVar is missing");
    }
    return fwdStFrequencyRangeIntVar->lowerRealBound;
}

double FourierProcessSVJ::upperRealBound(
    const FwdStFourierProductExpQuadVar& product, 
    const DateTime& matDate) const{
    if (!fwdStFrequencyRangeIntVar){
        throw ModelException("FourierProcessSVJ::lowerRealBound", 
                             "fwdStFrequencyRangeIntVar is missing");
    }
    return fwdStFrequencyRangeIntVar->upperRealBound;
}


//////////////////////////////////////////////////////////////////////////



/* Started integrated variance */
Complex FourierProcessSVJ::cumulant(const StFourierProductIntVar& product, 
                                    const Complex z, 
                                    const DateTime& matDate) const {
    return theVol->cumulant(static_cast<const StFourierProcessIntVar&>(*this),
                            product, 
                            z, 
                            matDate);
}

double FourierProcessSVJ::lowerRealBound(
    const StFourierProductIntVar& product, 
    const DateTime& matDate) const{
    if (!stFrequencyRangeIntVar){
        throw ModelException("FourierProcessSVJ::lowerRealBound", 
                             "StFrequencyRangeIntVar is missing");
    }
    return stFrequencyRangeIntVar->lowerRealBound;
}

double FourierProcessSVJ::upperRealBound(
    const StFourierProductIntVar& product, 
    const DateTime& matDate) const{
    if (!stFrequencyRangeIntVar){
        throw ModelException("FourierProcessSVJ::lowerRealBound", 
                             "StFrequencyRangeIntVar is missing");
    }
    return stFrequencyRangeIntVar->upperRealBound;
}
    
/* Forward Starting integrated variance */
Complex FourierProcessSVJ::cumulant(const FwdStFourierProductIntVar& product, 
                                    const Complex z, 
                                    const DateTime& matDate) const {        
    return theVol->cumulant(
        static_cast<const FwdStFourierProcessIntVar&>(*this),
        product,
        z, 
        matDate);
}    

double FourierProcessSVJ::lowerRealBound(
    const FwdStFourierProductIntVar& product, 
    const DateTime& matDate) const{
    if (!fwdStFrequencyRangeIntVar){
        throw ModelException("FourierProcessSVJ::lowerRealBound", 
                             "fwdStFrequencyRangeIntVar is missing");
    }
    return fwdStFrequencyRangeIntVar->lowerRealBound;
}

double FourierProcessSVJ::upperRealBound(
    const FwdStFourierProductIntVar& product, 
    const DateTime& matDate) const{
    if (!fwdStFrequencyRangeIntVar){
        throw ModelException("FourierProcessSVJ::lowerRealBound", 
                             "fwdStFrequencyRangeIntVar is missing");
    }
    return fwdStFrequencyRangeIntVar->upperRealBound;
}
   
/* Started Quadratic Variation */
Complex FourierProcessSVJ::cumulant(const StFourierProductQuadVar& product, 
                                    const Complex z, 
                                    const DateTime& matDate) const {
    return theVol->cumulant(static_cast<const StFourierProcessQuadVar&>(*this),
                            product, 
                            z, 
                            matDate);
}

double FourierProcessSVJ::lowerRealBound(
    const StFourierProductQuadVar& product, 
    const DateTime& matDate) const{
    if (!stFrequencyRangeIntVar){
        throw ModelException("FourierProcessSVJ::lowerRealBound", 
                             "StFrequencyRangeIntVar is missing");
    }
    return stFrequencyRangeIntVar->lowerRealBound;
}

double FourierProcessSVJ::upperRealBound(
    const StFourierProductQuadVar& product, 
    const DateTime& matDate) const{
    if (!stFrequencyRangeIntVar){
        throw ModelException("FourierProcessSVJ::lowerRealBound", 
                             "StFrequencyRangeIntVar is missing");
    }
    return stFrequencyRangeIntVar->upperRealBound;
}

double FourierProcessSVJ::expectation(const StFourierProductQuadVar& product, 
                                      const DateTime& matDate) const {
    return theVol->expectation(static_cast<const StFourierProcessQuadVar&>(*this),
                               product, 
                               matDate);
}

/* Forward Starting Quadratic Variation */
Complex FourierProcessSVJ::cumulant(const FwdStFourierProductQuadVar& product, 
                                    const Complex z, 
                                    const DateTime& matDate) const {
    return theVol->cumulant(static_cast<const FwdStFourierProcessQuadVar&>(*this),
                            product, 
                            z, 
                            matDate);
}

double FourierProcessSVJ::lowerRealBound(
    const FwdStFourierProductQuadVar& product, 
    const DateTime& matDate) const{
    if (!fwdStFrequencyRangeIntVar){
        throw ModelException("FourierProcessSVJ::lowerRealBound", 
                             "FwdStFrequencyRangeIntVar is missing");
    }
    return fwdStFrequencyRangeIntVar->lowerRealBound;
}

double FourierProcessSVJ::upperRealBound(
    const FwdStFourierProductQuadVar& product, 
    const DateTime& matDate) const{
    if (!fwdStFrequencyRangeIntVar){
        throw ModelException("FourierProcessSVJ::lowerRealBound", 
                             "FwdStFrequencyRangeIntVar is missing");
    }
    return fwdStFrequencyRangeIntVar->upperRealBound;
}

double FourierProcessSVJ::expectation(const FwdStFourierProductQuadVar& product, 
                                      const DateTime& matDate) const {
    return theVol->expectation(static_cast<const FwdStFourierProcessQuadVar&>(*this),
                               product, 
                               matDate);
}

/** Create a MarketDataFetcher which will be used by the [Fourier] model
 * for retrieving market data etc */
MarketDataFetcherSP FourierProcessSVJ::marketDataFetcher() const {
    return MarketDataFetcherSP(new MDFAssetVol(VolSVJ::TYPE->getName(), 
                                               VolSurface::TYPE->getName()));
}

const TimeMetric& FourierProcessSVJ::getTimeMetric() const{
    return *timeMetric;
}

string FourierProcessSVJ::getParameters() const {
    return FourierProcess::extractParameters(theVol.get());
}

void FourierProcessSVJ::validate(const FourierProduct* product){
    static const string method = "FourierProcessSV::validate";
    try{
        const IMultiFactors& mAsset = product->getMultiAsset();
        // validate size
        if (mAsset.NbAssets() != 1){
            throw ModelException(method, "only 1 asset suppported");
        }
        // copy the vol that's in the multiasset
        theVol = VolSVJSP::dynamicCast(VolRequestRaw::copyVolBase(mAsset, 0));
        // copy time metric
        timeMetric = TimeMetricSP(copy(&theVol->getTimeMetric()));
    }
    catch (exception& e){
        throw ModelException(e, method);
    }
}

FourierProcessSVJ::FourierProcessSVJ(): 
FourierProcessParametric(TYPE){}

IObject* FourierProcessSVJ::defaultCtor(){
    return new FourierProcessSVJ();
}

/** Invoked when Class is 'loaded' */
void FourierProcessSVJ::load(CClassSP& clazz){
    clazz->setPublic(); // make visible to EAS/spreadsheet
    REGISTER(FourierProcessSVJ, clazz);
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
 
CClassConstSP const FourierProcessSVJ::TYPE =
CClass::registerClassLoadMethod("FourierProcessSVJ", 
                                typeid(FourierProcessSVJ), load);

/* external symbol to allow class to be forced to be linked in */
bool FourierProcessSVJLoad(){
    return (FourierProcessSVJ::TYPE != 0);
}

DRLIB_END_NAMESPACE
