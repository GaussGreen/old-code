//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : FourierProcessIGOU.cpp
//
//   Date        : 09 March 04
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/VolRequestRaw.hpp"
#include "edginc/FourierProcessIGOU.hpp"
#include "edginc/FourierEngine.hpp"
#include "edginc/MDFAssetVol.hpp"

DRLIB_BEGIN_NAMESPACE

/* Started Log Return */
Complex FourierProcessIGOU::scalelessCumulant(const StFourierProductLogRtn& product, 
                                              const Complex& z, 
                                              const DateTime& matDate) const{
    return theVol->scalelessCumulant(static_cast<const StFourierProcessLogRtn&>(*this),
                                     product, 
                                     z, 
                                     matDate);
}

double FourierProcessIGOU::lowerRealBound(const StFourierProductLogRtn& product, 
                                          const DateTime& matDate) const{
    if (!stFrequencyRange){
        throw ModelException("FourierProcessIGOU::lowerRealBound", 
                             "stFrequencyRange is missing");
    }
    return stFrequencyRange->lowerRealBound;
}

double FourierProcessIGOU::upperRealBound(const StFourierProductLogRtn& product, 
                                          const DateTime& matDate) const{
    if (!stFrequencyRange){
        throw ModelException("FourierProcessIGOU::lowerRealBound", 
                             "stFrequencyRange is missing");
    }
    return stFrequencyRange->upperRealBound;
}

/* Fwd starting Log Return  */
Complex FourierProcessIGOU::scalelessCumulant(const FwdStFourierProductLogRtn& product, 
                                              const Complex& z, 
                                              const DateTime& matDate) const{
    return theVol->scalelessCumulant(static_cast<const FwdStFourierProcessLogRtn&>(*this),
                                     product, 
                                     z, 
                                     matDate);
}

double FourierProcessIGOU::lowerRealBound(const FwdStFourierProductLogRtn& product, 
                                          const DateTime& matDate) const{
    if (!fwdStFrequencyRange){
        throw ModelException("FourierProcessIGOU::lowerRealBound", 
                             "fwdStFrequencyRange is missing");
    }
    return fwdStFrequencyRange->lowerRealBound;
}

double FourierProcessIGOU::upperRealBound(const FwdStFourierProductLogRtn& product, 
                                          const DateTime& matDate) const{
    if (!fwdStFrequencyRange){
        throw ModelException("FourierProcessIGOU::lowerRealBound", 
                             "fwdStFrequencyRange is missing");
    }
    return fwdStFrequencyRange->upperRealBound;
}

/* Started integrated variance */
Complex FourierProcessIGOU::cumulant(const StFourierProductIntVar& product, 
                                     const Complex z, 
                                     const DateTime& matDate) const {
    return theVol->cumulant(static_cast<const StFourierProcessIntVar&>(*this),
                            product, 
                            z, 
                            matDate);
}

double FourierProcessIGOU::lowerRealBound(const StFourierProductIntVar& product, 
                                          const DateTime& matDate) const{
    if (!stFrequencyRangeIntVar){
        throw ModelException("FourierProcessIGOU::lowerRealBound", 
                             "StFrequencyRangeIntVar is missing");
    }
    return stFrequencyRangeIntVar->lowerRealBound;
}

double FourierProcessIGOU::upperRealBound(const StFourierProductIntVar& product, 
                                          const DateTime& matDate) const{
    if (!stFrequencyRangeIntVar){
        throw ModelException("FourierProcessIGOU::lowerRealBound", 
                             "StFrequencyRangeIntVar is missing");
    }
    return stFrequencyRangeIntVar->upperRealBound;
}
    
/* Forward Starting integrated variance */
Complex FourierProcessIGOU::cumulant(const FwdStFourierProductIntVar& product, 
                                     const Complex z, 
                                     const DateTime& matDate) const {        
    return theVol->cumulant(static_cast<const FwdStFourierProcessIntVar&>(*this),
                            product,
                            z, 
                            matDate);
}    

double FourierProcessIGOU::lowerRealBound(const FwdStFourierProductIntVar& product, 
                                          const DateTime& matDate) const{
    if (!fwdStFrequencyRangeIntVar){
        throw ModelException("FourierProcessIGOU::lowerRealBound", 
                             "fwdStFrequencyRangeIntVar is missing");
    }
    return fwdStFrequencyRangeIntVar->lowerRealBound;
}

double FourierProcessIGOU::upperRealBound(const FwdStFourierProductIntVar& product, 
                                          const DateTime& matDate) const{
    if (!fwdStFrequencyRangeIntVar){
        throw ModelException("FourierProcessIGOU::lowerRealBound", 
                             "fwdStFrequencyRangeIntVar is missing");
    }
    return fwdStFrequencyRangeIntVar->upperRealBound;
}
   

/** Create a MarketDataFetcher which will be used by the [Fourier] model
 * for retrieving market data etc */
MarketDataFetcherSP FourierProcessIGOU::marketDataFetcher() const {
    return MarketDataFetcherSP(new MDFAssetVol(VolIGOU::TYPE->getName(), 
                                               VolSurface::TYPE->getName()));
}

const TimeMetric& FourierProcessIGOU::getTimeMetric() const{
    return *timeMetric;
}

string FourierProcessIGOU::getParameters() const {
    return FourierProcess::extractParameters(theVol.get());
}

void FourierProcessIGOU::validate(const FourierProduct* product){
    static const string method = "FourierProcessIGOU::validate";
    try{
        const IMultiFactors& mAsset = product->getMultiAsset();
        // validate size
        if (mAsset.NbAssets() != 1){
            throw ModelException(method, "only 1 asset suppported");
        }
        // copy the vol that's in the multiasset
        theVol = VolIGOUSP::dynamicCast(
            VolRequestRaw::copyVolBase(mAsset, 0));
        // copy time metric
        timeMetric = TimeMetricSP(copy(&theVol->getTimeMetric()));
    }
    catch (exception& e){
        throw ModelException(e, method);
    }
}

FourierProcessIGOU::FourierProcessIGOU(): 
FourierProcessParametric(TYPE){}

IObject* FourierProcessIGOU::defaultCtor(){
    return new FourierProcessIGOU();
}

/** Invoked when Class is 'loaded' */
void FourierProcessIGOU::load(CClassSP& clazz){
    clazz->setPublic(); // make visible to EAS/spreadsheet
    REGISTER(FourierProcessIGOU, clazz);
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
 
CClassConstSP const FourierProcessIGOU::TYPE =
CClass::registerClassLoadMethod("FourierProcessIGOU", typeid(FourierProcessIGOU), load);

/* external symbol to allow class to be forced to be linked in */
bool FourierProcessIGOULoad(){
    return (FourierProcessIGOU::TYPE != 0);
}

DRLIB_END_NAMESPACE
