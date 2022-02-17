//----------------------------------------------------------------------------
//
//   Group       : Derivatives Research
//
//   Filename    : IRGridPointCache.cpp
//
//   Description : Class to help models implement 
//                 IRVegaPointwise::ISensitivePoints - this is primarily 
//                 targeted at SRM3 style models that use the swap dates as well
//                 as the swaptions dates
//
//   Author      : Mark A Robson
//
//   Date        : June 20, 2005
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#if defined(_MSC_VER)
// disable warning truncated decorated names information
#pragma warning(disable : 4503)
#endif
#include "edginc/IRGridPointCache.hpp"
#include "edginc/Hashtable.hpp"
#include ext_hash_map

DRLIB_BEGIN_NAMESPACE
/** Want to keep IRGridPoints per IR vol as well as per yield curve.
    This is so we can report what to tweak for IRVegaPointwise and we can
    correctly return an endDate for RhoPointwise (we assume it is impacted by 
    the pricing algorithm because it looks at coupons of swaps which take place
    beyond the end of instrument's end date) */
class IRGridPointCache::Imp{
public:
    typedef hash_map<string, IRGridPointAbsArraySP,
                     Hashtable::StringHash> HashMap; // for ease
    HashMap ptsByVol;
    HashMap ptsByYC;
    bool    cacheOn;
    Imp(bool cacheOn): cacheOn(cacheOn){}
};

IRGridPointCache::~IRGridPointCache(){
    // normally would use auto_ptr but reset method raises issues about 
    // different implementations
    delete my;
}

//// Simple constructor. cacheOn controls whether anything is cached
IRGridPointCache::IRGridPointCache(bool cacheOn): my(new Imp(cacheOn)){}

IRGridPointCache::IRGridPointCache(){}

/** Resets object, returns state before reset. cacheOn controls
    whether anything is cached */
void IRGridPointCache::setCachingMode(bool cacheOn){
    my->cacheOn= cacheOn;
}

void IRGridPointCache::cacheGridPoints(IYieldCurveConstSP yc,
                                        VolProcessedBSIRSP processedVol){
    if (my->cacheOn){
        const string& isoCode = yc->getCcy();
        const string& volName = processedVol->getName();
        // we assume that we treat the same currencies and vols consistently
        if (my->ptsByYC.find(isoCode) == my->ptsByYC.end() ||
            my->ptsByVol.find(volName) == my->ptsByVol.end()){
            IRGridPointAbsArraySP pts(processedVol->sensitiveIRVolPoints());
            my->ptsByYC[isoCode] = pts;
            my->ptsByVol[volName] = pts;
        }
    }
}
//// can happen if called in wrong way
static void checkForNullName(OutputNameConstSP name){
    if (!name){
        throw ModelException("IRGridPointCache::checkForNullName",
                             "Internal error - null name supplied");
    }
}

//// return the cached grid points for the specified ir vol
IRGridPointAbsArraySP IRGridPointCache::sensitiveIRVolPoints(
    OutputNameConstSP  vegaName) const{
    static const string method("QuantoCDSAlgorithm::IRGridPointCache::"
                               "sensitiveIRVolPoints");
    checkForNullName(vegaName);
    const string& name = vegaName->idGet(0);
    Imp::HashMap::const_iterator iter = my->ptsByVol.find(name);
    if (iter == my->ptsByVol.end() || vegaName->idCount() != 1){
        throw ModelException(method, "Internal error - no information for "
                             "IRVol with name "+vegaName->toString());
    }
    return iter->second;
}

/** Return when to stop tweaking for rho pointwise (as far as the swaps
    are concerned). */
DateTime IRGridPointCache::rhoEndDate(OutputNameConstSP  rhoName,
                                      const DateTime&    irVegaEndDate) const{
    static const string method("QuantoCDSAlgorithm::IRGridPointCache::"
                               "rhoEndDate");
    checkForNullName(rhoName);
    if (rhoName->idCount() != 1){
        throw ModelException(method, "Internal error - unexpected name "
                             "of yield curve: "+rhoName->toString());
    }
    const string& name = rhoName->idGet(0);
    Imp::HashMap::const_iterator iter = my->ptsByYC.find(name);
    if (iter != my->ptsByYC.end()){
        // get relevant points
        IRGridPointAbsArraySP pts = iter->second;
        // remove ones not uses
        IRGridPointAbs::trim(*pts, irVegaEndDate);
        return IRGridPointAbs::maxSwapMaturity(*pts);
    } else {
        return irVegaEndDate;
    }
}
    

DRLIB_END_NAMESPACE
