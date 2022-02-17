//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Description : Helper class for models to get vols for 'CAssets'
//                 (ie 'spot' type assets) out of market cache
//
//   Author      : Mark A Robson
//
//----------------------------------------------------------------------------

#ifndef EDR_MDFASSETVOL_HPP
#define EDR_MDFASSETVOL_HPP

#include "edginc/MarketDataFetcher.hpp"

DRLIB_BEGIN_NAMESPACE


/** Helper class for LN models to get data out of market cache */        
class MARKET_DLL MDFAssetVol:public MarketDataFetcher {
public:
    /** key method - uses specifed type to refine search when looking for
        'equity' and fx vols. Does not default to VolSurface when VolPreferred 
        specified (cf MarketDataFetcherLN) */
    virtual MarketObjectSP fetch(const MarketData*    market,
                                 const string&        name,
                                 const CClassConstSP& type,
                                 const IModel*        model) const;

    /** fxVolType allows to specify how to build fx vol 
        eg VolSurface => look for VolSurface and then wrap it
        so that it behaves like an FX Vol */
    MDFAssetVol(const string& volType, const string& fxVolType);

    /** Same as above but defaults to looking for any sort of fx vol */
    MDFAssetVol(const string& volType);

    void setAllowVolTypeConversions(bool allowVolTypeConversions);

protected:
    /** Methods to change the way this market data fetcher works - can be used
     * to simulate (on demand) a MarketDataFetcherLN - see MDFAssetVol */
    bool setCollapseToVolSurface (bool collapse);
    string setVolType(const string& newVolType);

    /** Fetch with conversion */
    MarketObjectSP convertAndfetch(const MarketData*    market,
                                   const string&        name,
                                   const CClassConstSP& type,
                                   const IModel*        model) const;

private:
    MDFAssetVol(const MDFAssetVol& rhs);
    MDFAssetVol& operator=(const MDFAssetVol& rhs);

    void initialise();
    // fields
    string volType;
    string fxVolType; /* eg VolSurface => look for VolSurface and then wrap it
                       * so that it behaves like an FX Vol */
    bool   allowVolTypeConversions;
    mutable CClassConstSP volClass; 
    mutable CClassConstSP fxVolClass; 

    bool collapseToVolSurface; /* If trying to fetch volPreferred and not in 
                                * the cache, collapse to VolSurface (or not) */
};

typedef smartConstPtr<MDFAssetVol> MDFAssetVolConstSP;
typedef smartPtr<MDFAssetVol> MDFAssetVolSP;

DRLIB_END_NAMESPACE
#endif
