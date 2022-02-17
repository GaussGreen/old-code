#include "edginc/config.hpp"
#include "edginc/VolRequestRaw.hpp"
#include "edginc/FourierProcessVSCurve.hpp"
#include "edginc/FourierEngine.hpp"
#include "edginc/MDFAssetVol.hpp"

DRLIB_BEGIN_NAMESPACE

/* Started Log Return */
Complex FourierProcessVSCurve::scalelessCumulant(const StFourierProductLogRtn& product, 
                                                 const Complex& z, 
                                                 const DateTime& matDate) const{
    return theVol->scalelessCumulant(static_cast<const StFourierProcessLogRtn&>(*this),
                                     product, 
                                     z, 
                                     matDate);
}

double FourierProcessVSCurve::lowerRealBound(const StFourierProductLogRtn& product, 
                                             const DateTime& matDate) const{
    if (!stFrequencyRange){
        throw ModelException("FourierProcessVSCurve::lowerRealBound", 
                             "stFrequencyRange is missing");
    }
    return stFrequencyRange->lowerRealBound;
}

double FourierProcessVSCurve::upperRealBound(const StFourierProductLogRtn& product, 
                                             const DateTime& matDate) const{
    if (!stFrequencyRange){
        throw ModelException("FourierProcessVSCurve::lowerRealBound", 
                             "stFrequencyRange is missing");
    }
    return stFrequencyRange->upperRealBound;
}

/* Fwd starting Log Return  */
Complex FourierProcessVSCurve::scalelessCumulant(const FwdStFourierProductLogRtn& product, 
                                                 const Complex& z,
                                                 const DateTime& matDate) const{
    return theVol->scalelessCumulant(static_cast<const FwdStFourierProcessLogRtn&>(*this),
                                     product, 
                                     z, 
                                     matDate);
}

double FourierProcessVSCurve::lowerRealBound(const FwdStFourierProductLogRtn& product, 
                                             const DateTime& matDate) const{
    if (!fwdStFrequencyRange){
        throw ModelException("FourierProcessVSCurve::lowerRealBound", 
                             "fwdStFrequencyRange is missing");
    }
    return fwdStFrequencyRange->lowerRealBound;
}

double FourierProcessVSCurve::upperRealBound(const FwdStFourierProductLogRtn& product, 
                                             const DateTime& matDate) const{
    if (!fwdStFrequencyRange){
        throw ModelException("FourierProcessVSCurve::lowerRealBound", 
                             "fwdStFrequencyRange is missing");
    }
    return fwdStFrequencyRange->upperRealBound;
}


/* Fwd Starting Expected Quad Var */
Complex FourierProcessVSCurve::cumulant(
    const FwdStFourierProductExpQuadVar& product, 
    const Complex& z, 
    const DateTime& matDate) const{
    return theVol->cumulant(
        static_cast<const FwdStFourierProcessExpQuadVar&>(*this),
        product, 
        z, 
        matDate);
}

double FourierProcessVSCurve::lowerRealBound(
    const FwdStFourierProductExpQuadVar& product, 
    const DateTime& matDate) const{
    if (!fwdStFrequencyRangeIntVar){
        throw ModelException("FourierProcessVSCurve::lowerRealBound", 
                             "fwdStFrequencyRangeIntVar is missing");
    }
    return fwdStFrequencyRangeIntVar->lowerRealBound;
}

double FourierProcessVSCurve::upperRealBound(
    const FwdStFourierProductExpQuadVar& product, 
    const DateTime& matDate) const{
    if (!fwdStFrequencyRangeIntVar){
        throw ModelException("FourierProcessVSCurve::lowerRealBound", 
                             "fwdStFrequencyRangeIntVar is missing");
    }
    return fwdStFrequencyRangeIntVar->upperRealBound;
}



/* Started integrated variance */
Complex FourierProcessVSCurve::cumulant(const StFourierProductIntVar& product, 
                                        const Complex z, 
                                        const DateTime& matDate) const {
    return theVol->cumulant(static_cast<const StFourierProcessIntVar&>(*this),
                            product, 
                            z, 
                            matDate);
}

double FourierProcessVSCurve::lowerRealBound(const StFourierProductIntVar& product, 
                                             const DateTime& matDate) const{
    if (!stFrequencyRangeIntVar){
        throw ModelException("FourierProcessVSCurve::lowerRealBound", 
                             "StFrequencyRangeIntVar is missing");
    }
    return stFrequencyRangeIntVar->lowerRealBound;
}

double FourierProcessVSCurve::upperRealBound(const StFourierProductIntVar& product, 
                                             const DateTime& matDate) const{
    if (!stFrequencyRangeIntVar){
        throw ModelException("FourierProcessVSCurve::lowerRealBound", 
                             "StFrequencyRangeIntVar is missing");
    }
    return stFrequencyRangeIntVar->upperRealBound;
}
    
/* Forward Starting integrated variance */
Complex FourierProcessVSCurve::cumulant(const FwdStFourierProductIntVar& product, 
                                        const Complex z, 
                                        const DateTime& matDate) const {        
    return theVol->cumulant(static_cast<const FwdStFourierProcessIntVar&>(*this),
                            product,
                            z, 
                            matDate);
}    

double FourierProcessVSCurve::lowerRealBound(const FwdStFourierProductIntVar& product, 
                                             const DateTime& matDate) const{
    if (!fwdStFrequencyRangeIntVar){
        throw ModelException("FourierProcessVSCurve::lowerRealBound", 
                             "fwdStFrequencyRangeIntVar is missing");
    }
    return fwdStFrequencyRangeIntVar->lowerRealBound;
}

double FourierProcessVSCurve::upperRealBound(const FwdStFourierProductIntVar& product, 
                                             const DateTime& matDate) const{
    if (!fwdStFrequencyRangeIntVar){
        throw ModelException("FourierProcessVSCurve::lowerRealBound", 
                             "fwdStFrequencyRangeIntVar is missing");
    }
    return fwdStFrequencyRangeIntVar->upperRealBound;
}

/* Started quad variation */
Complex FourierProcessVSCurve::cumulant(const StFourierProductQuadVar& product, 
                                        const Complex z, 
                                        const DateTime& matDate) const{
    return( theVol->cumulant(static_cast<const StFourierProcessQuadVar&>(*this),
                            product, 
                            z, 
                            matDate) );
}

double FourierProcessVSCurve::lowerRealBound(const StFourierProductQuadVar& product, 
                                             const DateTime& matDate) const{
    if (!stFrequencyRangeIntVar){
        throw ModelException("FourierProcessVSCurve::lowerRealBound", 
                             "StFrequencyRangeIntVar is missing");
    }
    return stFrequencyRangeIntVar->lowerRealBound;
}

double FourierProcessVSCurve::upperRealBound(const StFourierProductQuadVar& product, 
                                             const DateTime& matDate) const{
    if (!stFrequencyRangeIntVar){
        throw ModelException("FourierProcessVSCurve::lowerRealBound", 
                             "StFrequencyRangeIntVar is missing");
    }
    return stFrequencyRangeIntVar->upperRealBound;
}

double FourierProcessVSCurve::expectation(const StFourierProductQuadVar& product, 
                                          const DateTime& matDate) const {
    return theVol->expectation(static_cast<const StFourierProcessQuadVar&>(*this),
                               product, 
                               matDate);
}
    
/* Forward Starting quadratic variation */
Complex FourierProcessVSCurve::cumulant(const FwdStFourierProductQuadVar& product, 
                                        const Complex z, 
                                        const DateTime& matDate) const{        
    return( theVol->cumulant(static_cast<const FwdStFourierProcessQuadVar&>(*this),
                             product,
                             z, 
                             matDate) );
}

double FourierProcessVSCurve::lowerRealBound(const FwdStFourierProductQuadVar& product, 
                                             const DateTime& matDate) const{
    if (!fwdStFrequencyRangeIntVar){
        throw ModelException("FourierProcessVSCurve::lowerRealBound", 
                             "fwdStFrequencyRangeIntVar is missing");
    }
    return fwdStFrequencyRangeIntVar->lowerRealBound;
}

double FourierProcessVSCurve::upperRealBound(const FwdStFourierProductQuadVar& product, 
                                         const DateTime& matDate) const{
    if (!fwdStFrequencyRangeIntVar){
        throw ModelException("FourierProcessVSCurve::lowerRealBound", 
                             "fwdStFrequencyRangeIntVar is missing");
    }
    return fwdStFrequencyRangeIntVar->upperRealBound;
}

double FourierProcessVSCurve::expectation(const FwdStFourierProductQuadVar& product, 
                                          const DateTime& matDate) const{
    return theVol->expectation(static_cast<const FwdStFourierProcessQuadVar&>(*this),
                               product, 
                               matDate);
}
  
/** Create a MarketDataFetcher which will be used by the [Fourier] model
 * for retrieving market data etc */
MarketDataFetcherSP FourierProcessVSCurve::marketDataFetcher() const {
    return MarketDataFetcherSP(new MDFAssetVol(VolVSCurve::TYPE->getName(), 
                                               VolSurface::TYPE->getName()));
}

const TimeMetric& FourierProcessVSCurve::getTimeMetric() const{
    return *timeMetric;
}

string FourierProcessVSCurve::getParameters() const {
    return FourierProcess::extractParameters(theVol.get());
}

void FourierProcessVSCurve::validate(const FourierProduct* product){
    static const string method = "FourierProcessVSCurve::validate";
    try{
        const IMultiFactors& mAsset = product->getMultiAsset();
        // validate size
        if (mAsset.NbAssets() != 1){
            throw ModelException(method, "only 1 asset suppported");
        }
        // copy the vol that's in the multiasset
        theVol = VolVSCurveSP::dynamicCast(VolRequestRaw::copyVolBase(mAsset, 0));
        // copy time metric
        timeMetric = TimeMetricSP(copy(&theVol->getTimeMetric()));
    }
    catch (exception& e){
        throw ModelException(e, method);
    }
}

FourierProcessVSCurve::FourierProcessVSCurve(): 
FourierProcessParametric(TYPE){}

IObject* FourierProcessVSCurve::defaultCtor(){
    return new FourierProcessVSCurve();
}

/** Invoked when Class is 'loaded' */
void FourierProcessVSCurve::load(CClassSP& clazz){
    clazz->setPublic(); // make visible to EAS/spreadsheet
    REGISTER(FourierProcessVSCurve, clazz);
    SUPERCLASS(FourierProcessParametric);
    IMPLEMENTS(StFourierProcessLogRtn);
    IMPLEMENTS(FwdStFourierProcessLogRtn);
    
    IMPLEMENTS(FwdStFourierProcessIntVar);
    IMPLEMENTS(StFourierProcessIntVar);
    
    IMPLEMENTS(StFourierProcessQuadVar);
    IMPLEMENTS(FwdStFourierProcessQuadVar);
    
    EMPTY_SHELL_METHOD(defaultCtor);
    
    FIELD(stFrequencyRange, "stFrequencyRange");
    FIELD_MAKE_OPTIONAL(stFrequencyRange);
    FIELD(fwdStFrequencyRange, "fwdStFrequencyRange");
    FIELD_MAKE_OPTIONAL(fwdStFrequencyRange);
    
    FIELD(fwdStFrequencyRangeIntVar, "fwdStFrequencyRangeIntVar");
    FIELD_MAKE_OPTIONAL(fwdStFrequencyRangeIntVar);
    FIELD(stFrequencyRangeIntVar, "stFrequencyRangeIntVar");
    FIELD_MAKE_OPTIONAL(stFrequencyRangeIntVar);
    
    FIELD(fwdStFrequencyRangeQuadVar, "fwdStFrequencyRangeQuadVar");
    FIELD_MAKE_OPTIONAL(fwdStFrequencyRangeQuadVar);
    FIELD(stFrequencyRangeQuadVar, "stFrequencyRangeQuadVar");
    FIELD_MAKE_OPTIONAL(stFrequencyRangeQuadVar);
    
    // transient
    FIELD(theVol, ""); 
    FIELD_MAKE_TRANSIENT_BUT_TWEAKABLE(theVol);
    FIELD(timeMetric, "");
    FIELD_MAKE_TRANSIENT(timeMetric);
}

CClassConstSP const FourierProcessVSCurve::TYPE =
CClass::registerClassLoadMethod("FourierProcessVSCurve", typeid(FourierProcessVSCurve), load);

/* external symbol to allow class to be forced to be linked in */
bool FourierProcessVSCurveLoad(){
    return (FourierProcessVSCurve::TYPE != 0);
}

DRLIB_END_NAMESPACE
