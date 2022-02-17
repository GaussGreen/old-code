//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : VolBaseParam.cpp
//
//   Description : A parameterised CVolBase which supports VolRequestRaw
//
//   Date        : 17 May 2005
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/VolBaseParam.hpp"
#include "edginc/VolRequestRaw.hpp"

DRLIB_BEGIN_NAMESPACE

VolBaseParam::~VolBaseParam(){}

VolBaseParam::VolBaseParam(CClassConstSP clazz): CVolBaseParamSurface(clazz){}

VolBaseParam::VolBaseParam(CClassConstSP clazz, const string& name): 
CVolBaseParamSurface(clazz, name) {}

const TimeMetric& VolBaseParam::getTimeMetric() const{
    return *timeMetric;
}

 // unfortunate - to tidy up
TimeMetricConstSP VolBaseParam::GetTimeMetric() const{
    return timeMetric;
}

double VolBaseParam::calcTradingTime(const DateTime &date1, 
                                     const DateTime &date2) const{
    return timeMetric->yearFrac(date1, date2);
}
 
/** Retrieves time metric from parent's vol surface */
void VolBaseParam::getMarket(const IModel*     model, 
                             const MarketData* market,
                             const string&     name){
    CVolBaseParamSurface::getMarket(model, market, name);
    // Dummy try/catch to avoid vc6.opt crash
	try {
        const VolSurface* backbone = getBackboneSurface();
        timeMetric = backbone->getTimeMetric();
    } catch (...) { throw; }
}

/** Retrieves time metric from parent's vol surface */
void VolBaseParam::getMarket(const IModel* model, 
                             const MarketData* market){
    getMarket(model, market, CVolBaseParamSurface::getName());
}

IVolProcessed* VolBaseParam::getProcessedVol(
    const CVolRequest* volRequest,
    const CAsset*      asset) const{
    // intercept VolRequestRaw here
    if (VolRequestRaw::TYPE->isInstance(volRequest)){
        return const_cast<VolBaseParam*>(this);
    }
    // otherwise delegate to base class
    return CVolBaseParamSurface::getProcessedVol(volRequest, asset);
}
    

IVolProcessed* VolBaseParam::getProcessedVol(
    const CVolRequest* volRequest,
    const CAsset*      eqAsset,
    const FXAsset*     fxAsset,
    const Correlation* eqFXCorr) const{
    // intercept VolRequestRaw here
    if (VolRequestRaw::TYPE->isInstance(volRequest)){
        return const_cast<VolBaseParam*>(this); // meaningless? no?
    }
    // otherwise delegate to base class
    return CVolBaseParamSurface::getProcessedVol(volRequest, 
                                                 eqAsset, 
                                                 fxAsset, 
                                                 eqFXCorr);
}

void VolBaseParam::load(CClassSP& clazz){
    REGISTER(VolBaseParam, clazz);
    SUPERCLASS(CVolBaseParamSurface);
    IMPLEMENTS(IVolProcessed);
    FIELD_NO_DESC(timeMetric);
    FIELD_MAKE_TRANSIENT(timeMetric);
}

CClassConstSP const VolBaseParam::TYPE =
CClass::registerClassLoadMethod("VolBaseParam", typeid(VolBaseParam), load);

DRLIB_END_NAMESPACE
