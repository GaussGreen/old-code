//----------------------------------------------------------------------------
//
//   Group       : Derivatives Research
//
//   Filename    : IRGridPointCache.hpp
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

#ifndef QLIB_IRGRIDPOINTCACHE_HPP
#define QLIB_IRGRIDPOINTCACHE_HPP
#include "edginc/IRGridPointAbs.hpp"
#include "edginc/VirtualDestructorBase.hpp"
#include "edginc/YieldCurve.hpp"
#include "edginc/VolProcessedBSIR.hpp"
#include "edginc/OutputName.hpp"

DRLIB_BEGIN_NAMESPACE
class IRGridPointCache;
typedef smartPtr<IRGridPointCache> IRGridPointCacheSP;

/** Class to help models implement IRVegaPointwise::ISensitivePoints - this is
    primarily targeted at SRM3 style models that use the swap dates as well
    as the swaptions dates.

    Want to keep IRGridPoints per IR vol as well as per yield curve.
    This is so we can report what to tweak for IRVegaPointwise and we can
    correctly return an endDate for RhoPointwise (we assume it is impacted by 
    the pricing algorithm because it looks at coupons of swaps which take place
    beyond the end of instrument's end date) */
class MCARLO_DLL IRGridPointCache: public virtual VirtualDestructorBase{
public:
    virtual ~IRGridPointCache();

    //// Simple constructor. cacheOn controls whether anything is cached
    IRGridPointCache(bool cacheOn);

    /** Resets object, returns state before reset. cacheOn controls
     * whether anything is cached */
    void setCachingMode(bool cacheOn);

    void cacheGridPoints(IYieldCurveConstSP yc,
                         VolProcessedBSIRSP processedVol);

    //// return the cached grid points for the specified ir vol
    IRGridPointAbsArraySP sensitiveIRVolPoints(
        OutputNameConstSP  vegaName) const;

    /** Return when to stop tweaking for rho pointwise (as far as the swaps
        are concerned). */
    DateTime rhoEndDate(OutputNameConstSP  rhoName,
                        const DateTime&    irVegaEndDate) const;
private:
    IRGridPointCache();
    // hide implementation in separate class
    class Imp;
    Imp*   my;
};

DRLIB_END_NAMESPACE
#endif
