//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : FourierProcessSV.cpp
//
//   Date        : 09 March 04
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/VolRequestRaw.hpp"
#include "edginc/FourierProcessSV.hpp"
#include "edginc/FourierEngine.hpp"
#include "edginc/MDFAssetVol.hpp"

DRLIB_BEGIN_NAMESPACE

/* Started Log Return */
Complex FourierProcessSV::scalelessCumulant(
    const StFourierProductLogRtn& product, 
    const Complex& z, 
    const DateTime& matDate) const{
    return theVol->scalelessCumulant(
        static_cast<const StFourierProcessLogRtn&>(*this), 
        product, 
        z, matDate);
}

double FourierProcessSV::lowerRealBound(
    const StFourierProductLogRtn& product, 
    const DateTime& matDate) const{
    if (!stFrequencyRange){
        throw ModelException("FourierProcessSV::lowerRealBound", 
                             "stFrequencyRange is missing");
    }
    return stFrequencyRange->lowerRealBound;
}

double FourierProcessSV::upperRealBound(
    const StFourierProductLogRtn& product, 
    const DateTime& matDate) const{
    if (!stFrequencyRange){
        throw ModelException("FourierProcessSV::lowerRealBound", 
                             "stFrequencyRange is missing");
    }
    return stFrequencyRange->upperRealBound;
}

/* Fwd starting Log Return  */
Complex FourierProcessSV::scalelessCumulant(
    const FwdStFourierProductLogRtn& product, 
    const Complex& z, 
    const DateTime& matDate) const{
    return theVol->scalelessCumulant(
        static_cast<const FwdStFourierProcessLogRtn&>(*this),
        product, 
        z, 
        matDate);
}

double FourierProcessSV::lowerRealBound(
    const FwdStFourierProductLogRtn& product, 
    const DateTime& matDate) const{
    if (!fwdStFrequencyRange){
        throw ModelException("FourierProcessSV::lowerRealBound", 
                             "fwdStFrequencyRange is missing");
    }
    return fwdStFrequencyRange->lowerRealBound;
}

double FourierProcessSV::upperRealBound(
    const FwdStFourierProductLogRtn& product, 
    const DateTime& matDate) const{
    if (!fwdStFrequencyRange){
        throw ModelException("FourierProcessSV::lowerRealBound", 
                             "fwdStFrequencyRange is missing");
    }
    return fwdStFrequencyRange->upperRealBound;
}



/**ARNAUD**////////////////////////////////////////////////////////////
/* Fwd Starting Expected Quad Var */
Complex FourierProcessSV::cumulant(
    const FwdStFourierProductExpQuadVar& product, 
    const Complex& z, 
    const DateTime& matDate) const{
    return theVol->cumulant(
        static_cast<const FwdStFourierProcessExpQuadVar&>(*this),
        product, 
        z, 
        matDate);
}

double FourierProcessSV::lowerRealBound(
    const FwdStFourierProductExpQuadVar& product, 
    const DateTime& matDate) const{
    if (!fwdStFrequencyRangeIntVar){
        throw ModelException("FourierProcessSV::lowerRealBound", 
                             "fwdStFrequencyRangeIntVar is missing");
    }
    return fwdStFrequencyRangeIntVar->lowerRealBound;
}

double FourierProcessSV::upperRealBound(
    const FwdStFourierProductExpQuadVar& product, 
    const DateTime& matDate) const{
    if (!fwdStFrequencyRangeIntVar){
        throw ModelException("FourierProcessSV::lowerRealBound", 
                             "fwdStFrequencyRangeIntVar is missing");
    }
    return fwdStFrequencyRangeIntVar->upperRealBound;
}


//////////////////////////////////////////////////////////////////////////

/* Started integrated variance */
Complex FourierProcessSV::cumulant(const StFourierProductIntVar& product, 
                                    const Complex z, 
                                    const DateTime& matDate) const {
    return theVol->cumulant(static_cast<const StFourierProcessIntVar&>(*this),
                            product, 
                            z, 
                            matDate);
}

double FourierProcessSV::lowerRealBound(
    const StFourierProductIntVar& product, 
    const DateTime& matDate) const{
    if (!stFrequencyRangeIntVar){
        throw ModelException("FourierProcessSV::lowerRealBound", 
                             "StFrequencyRangeIntVar is missing");
    }
    return stFrequencyRangeIntVar->lowerRealBound;
}

double FourierProcessSV::upperRealBound(
    const StFourierProductIntVar& product, 
    const DateTime& matDate) const{
    if (!stFrequencyRangeIntVar){
        throw ModelException("FourierProcessSV::lowerRealBound", 
                             "StFrequencyRangeIntVar is missing");
    }
    return stFrequencyRangeIntVar->upperRealBound;
}
    
/* Forward Starting integrated variance */
Complex FourierProcessSV::cumulant(const FwdStFourierProductIntVar& product, 
                                    const Complex z, 
                                    const DateTime& matDate) const {        
    return theVol->cumulant(
        static_cast<const FwdStFourierProcessIntVar&>(*this),
        product,
        z, 
        matDate);
}    

double FourierProcessSV::lowerRealBound(
    const FwdStFourierProductIntVar& product, 
    const DateTime& matDate) const{
    if (!fwdStFrequencyRangeIntVar){
        throw ModelException("FourierProcessSV::lowerRealBound", 
                             "fwdStFrequencyRangeIntVar is missing");
    }
    return fwdStFrequencyRangeIntVar->lowerRealBound;
}

double FourierProcessSV::upperRealBound(
    const FwdStFourierProductIntVar& product, 
    const DateTime& matDate) const{
    if (!fwdStFrequencyRangeIntVar){
        throw ModelException("FourierProcessSV::lowerRealBound", 
                             "fwdStFrequencyRangeIntVar is missing");
    }
    return fwdStFrequencyRangeIntVar->upperRealBound;
}
   
/* Started Quadratic Variation */
Complex FourierProcessSV::cumulant(const StFourierProductQuadVar& product, 
                                    const Complex z, 
                                    const DateTime& matDate) const {
    return theVol->cumulant(static_cast<const StFourierProcessQuadVar&>(*this),
                            product, 
                            z, 
                            matDate);
}

double FourierProcessSV::lowerRealBound(
    const StFourierProductQuadVar& product, 
    const DateTime& matDate) const{
    if (!stFrequencyRangeIntVar){
        throw ModelException("FourierProcessSV::lowerRealBound", 
                             "StFrequencyRangeIntVar is missing");
    }
    return stFrequencyRangeIntVar->lowerRealBound;
}

double FourierProcessSV::upperRealBound(
    const StFourierProductQuadVar& product, 
    const DateTime& matDate) const{
    if (!stFrequencyRangeIntVar){
        throw ModelException("FourierProcessSV::lowerRealBound", 
                             "StFrequencyRangeIntVar is missing");
    }
    return stFrequencyRangeIntVar->upperRealBound;
}

double FourierProcessSV::expectation(const StFourierProductQuadVar& product, 
                                     const DateTime& matDate) const {
    return theVol->expectation(static_cast<const StFourierProcessQuadVar&>(*this),
                               product, 
                               matDate);
}


/* Forward starting Quadratic Variation */
Complex FourierProcessSV::cumulant(const FwdStFourierProductQuadVar& product, 
                                    const Complex z, 
                                    const DateTime& matDate) const {
    return theVol->cumulant(static_cast<const FwdStFourierProcessQuadVar&>(*this),
                            product, 
                            z, 
                            matDate);
}

double FourierProcessSV::lowerRealBound(
    const FwdStFourierProductQuadVar& product, 
    const DateTime& matDate) const{
    if (!fwdStFrequencyRangeIntVar){
        throw ModelException("FourierProcessSV::lowerRealBound", 
                             "FwdStFrequencyRangeIntVar is missing");
    }
    return fwdStFrequencyRangeIntVar->lowerRealBound;
}

double FourierProcessSV::upperRealBound(
    const FwdStFourierProductQuadVar& product, 
    const DateTime& matDate) const{
    if (!fwdStFrequencyRangeIntVar){
        throw ModelException("FourierProcessSV::lowerRealBound", 
                             "FwdStFrequencyRangeIntVar is missing");
    }
    return fwdStFrequencyRangeIntVar->upperRealBound;
}

double FourierProcessSV::expectation(const FwdStFourierProductQuadVar& product, 
                                     const DateTime& matDate) const {
    return theVol->expectation(static_cast<const FwdStFourierProcessQuadVar&>(*this),
                               product, 
                               matDate);
}

/** Create a MarketDataFetcher which will be used by the [Fourier] model
 * for retrieving market data etc */
MarketDataFetcherSP FourierProcessSV::marketDataFetcher() const {
    return MarketDataFetcherSP(new MDFAssetVol(VolSV::TYPE->getName(), 
                                               VolSurface::TYPE->getName()));
}

const TimeMetric& FourierProcessSV::getTimeMetric() const{
    return *timeMetric;
}

string FourierProcessSV::getParameters() const {
    return FourierProcess::extractParameters(theVol.get());
}

void FourierProcessSV::validate(const FourierProduct* product){
    static const string method = "FourierProcessSV::validate";
    try{
        const IMultiFactors& mAsset = product->getMultiAsset();
        // validate size
        if (mAsset.NbAssets() != 1){
            throw ModelException(method, "only 1 asset suppported");
        }
        // copy the vol that's in the multiasset
        theVol = VolSVSP::dynamicCast(VolRequestRaw::copyVolBase(mAsset, 0));
        // copy time metric
        timeMetric = TimeMetricSP(copy(&theVol->getTimeMetric()));
    }
    catch (exception& e){
        throw ModelException(e, method);
    }
}

FourierProcessSV::FourierProcessSV(): 
FourierProcessParametric(TYPE){}

IObject* FourierProcessSV::defaultCtor(){
    return new FourierProcessSV();
}

/** Invoked when Class is 'loaded' */
void FourierProcessSV::load(CClassSP& clazz){
    clazz->setPublic(); // make visible to EAS/spreadsheet
    REGISTER(FourierProcessSV, clazz);
    SUPERCLASS(FourierProcessParametric);
    IMPLEMENTS(StFourierProcessLogRtn);
    IMPLEMENTS(FwdStFourierProcessLogRtn);
    IMPLEMENTS(FwdStFourierProcessIntVar);
    IMPLEMENTS(StFourierProcessIntVar);
    IMPLEMENTS(StFourierProcessQuadVar);
	IMPLEMENTS(FwdStFourierProcessExpQuadVar);
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

CClassConstSP const FourierProcessSV::TYPE =
CClass::registerClassLoadMethod("FourierProcessSV", 
                                typeid(FourierProcessSV), load);

/* external symbol to allow class to be forced to be linked in */
bool FourierProcessSVLoad(){
    return (FourierProcessSV::TYPE != 0);
}

DRLIB_END_NAMESPACE
