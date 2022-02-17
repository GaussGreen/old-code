//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : SpotPriceProxy.hpp
//
//   Description : Spot price for fund proxies
//
//   Author      : Andrew J Swain
//
//   Date        : 24 January 2002
//
//
//----------------------------------------------------------------------------

#ifndef _SPOTPRICEPROXY_HPP
#define _SPOTPRICEPROXY_HPP
#include "edginc/Sensitivity.hpp"

DRLIB_BEGIN_NAMESPACE
class Results;

/** SpotPriceProxy Sensitivity. Effectively does a DELTA shift then reports
    the initial spot price as the sensitivity. This functionality is 
    needed by analytics */
class RISKMGR_DLL SpotPriceProxy: public Sensitivity {
public:
    static CClassConstSP const TYPE;
    const static string NAME;
   
    /** identifies the name used storing associated results in the output */
    const string& getSensOutputName() const;

    /** Is this sensitivity made using a discrete shift (ie a jump) or a
        an approximately continuous one. Returns false */
    virtual bool discreteShift() const;

    /** Combines spot prices between results packets (ie does a merge) */
    virtual void addResult(Results*           results,     // (M)
                           const Results*     resultsToAdd,
                           double             scaleFactor) const;
    
    
    /** Does a DELTA shift and records the spot price */
    virtual void calculate(TweakGroup*     tweakGroup,
                           Results*        results);

    SpotPriceProxy();

private:
    friend class SpotPriceProxyHelper;

    SpotPriceProxy(const SpotPriceProxy &rhs);
    SpotPriceProxy& operator=(const SpotPriceProxy& rhs);
};

typedef smartConstPtr<SpotPriceProxy> SpotPriceProxyConstSP;
typedef smartPtr<SpotPriceProxy> SpotPriceProxySP;


DRLIB_END_NAMESPACE

#endif

