//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : FourierProcessCGMYHeston.cpp
//
//   Date        : 10 March 04
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/VolRequestRaw.hpp"
#include "edginc/FourierProcessCGMYHeston.hpp"
#include "edginc/FourierEngine.hpp"
#include "edginc/MDFAssetVol.hpp"

DRLIB_BEGIN_NAMESPACE

/* Started Log Return */
Complex FourierProcessCGMYHeston::scalelessCumulant(const StFourierProductLogRtn& product, 
                                                    const Complex& z, 
                                                    const DateTime& matDate) const{
    return theVol->scalelessCumulant(static_cast<const StFourierProcessLogRtn&>(*this),
                                     product, 
                                     z, 
                                     matDate);
}

double FourierProcessCGMYHeston::lowerRealBound(const StFourierProductLogRtn& product, 
                                                const DateTime& matDate) const{
    if (!stFrequencyRange){
        throw ModelException("FourierProcessCGMYHeston::lowerRealBound", 
                             "stFrequencyRange is missing");
    }
    return stFrequencyRange->lowerRealBound;
}

double FourierProcessCGMYHeston::upperRealBound(const StFourierProductLogRtn& product, 
                                                const DateTime& matDate) const{
    if (!stFrequencyRange){
        throw ModelException("FourierProcessCGMYHeston::lowerRealBound", 
                             "stFrequencyRange is missing");
    }
    return stFrequencyRange->upperRealBound;
}

/* Fwd starting Log Return  */
Complex FourierProcessCGMYHeston::scalelessCumulant(const FwdStFourierProductLogRtn& product, 
                                                    const Complex& z, 
                                                    const DateTime& matDate) const{
    return theVol->scalelessCumulant(static_cast<const FwdStFourierProcessLogRtn&>(*this),
                                     product, 
                                     z, 
                                     matDate);
}

double FourierProcessCGMYHeston::lowerRealBound(const FwdStFourierProductLogRtn& product, 
                                                const DateTime& matDate) const{
    if (!fwdStFrequencyRange){
        throw ModelException("FourierProcessCGMYHeston::lowerRealBound", 
                             "fwdStFrequencyRange is missing");
    }
    return fwdStFrequencyRange->lowerRealBound;
}

double FourierProcessCGMYHeston::upperRealBound(const FwdStFourierProductLogRtn& product, 
                                                const DateTime& matDate) const{
    if (!fwdStFrequencyRange){
        throw ModelException("FourierProcessCGMYHeston::lowerRealBound", 
                             "fwdStFrequencyRange is missing");
    }
    return fwdStFrequencyRange->upperRealBound;
}

/** Create a MarketDataFetcher which will be used by the [Fourier] model
 * for retrieving market data etc */
MarketDataFetcherSP FourierProcessCGMYHeston::marketDataFetcher() const {
    return MarketDataFetcherSP(new MDFAssetVol(VolCGMYHeston::TYPE->getName(), 
                                               VolSurface::TYPE->getName()));
}

const TimeMetric& FourierProcessCGMYHeston::getTimeMetric() const{
    return *timeMetric;
}

string FourierProcessCGMYHeston::getParameters() const {
    return FourierProcess::extractParameters(theVol.get());
}

void FourierProcessCGMYHeston::validate(const FourierProduct* product){
    static const string method = "FourierProcessCGMYHeston::validate";
    try{
        const IMultiFactors& mAsset = product->getMultiAsset();
        // validate size
        if (mAsset.NbAssets() != 1){
            throw ModelException(method, "only 1 asset suppported");
        }
        // copy the vol that's in the multiasset
        theVol = VolCGMYHestonSP::dynamicCast(
            VolRequestRaw::copyVolBase(mAsset, 0));
        // copy time metric
        timeMetric = TimeMetricSP(copy(&theVol->getTimeMetric()));
    }
    catch (exception& e){
        throw ModelException(e, method);
    }
}

FourierProcessCGMYHeston::FourierProcessCGMYHeston(): 
FourierProcessParametric(TYPE){}

IObject* FourierProcessCGMYHeston::defaultCtor(){
    return new FourierProcessCGMYHeston();
}

/** Invoked when Class is 'loaded' */
void FourierProcessCGMYHeston::load(CClassSP& clazz){
    clazz->setPublic(); // make visible to EAS/spreadsheet
    REGISTER(FourierProcessCGMYHeston, clazz);
    SUPERCLASS(FourierProcessParametric);
    IMPLEMENTS(StFourierProcessLogRtn);
    IMPLEMENTS(FwdStFourierProcessLogRtn);
    EMPTY_SHELL_METHOD(defaultCtor);
    FIELD(stFrequencyRange, "stFrequencyRange");
    FIELD_MAKE_OPTIONAL(stFrequencyRange);
    FIELD(fwdStFrequencyRange, "fwdStFrequencyRange");
    FIELD_MAKE_OPTIONAL(fwdStFrequencyRange);
    // transient
    FIELD(theVol, ""); 
    FIELD_MAKE_TRANSIENT_BUT_TWEAKABLE(theVol);
    FIELD(timeMetric, "");
    FIELD_MAKE_TRANSIENT(timeMetric);
}
 
CClassConstSP const FourierProcessCGMYHeston::TYPE =
CClass::registerClassLoadMethod("FourierProcessCGMYHeston", typeid(FourierProcessCGMYHeston), load);

/* external symbol to allow class to be forced to be linked in */
bool FourierProcessCGMYHestonLoad(){
    return (FourierProcessCGMYHeston::TYPE != 0);
}

DRLIB_END_NAMESPACE
