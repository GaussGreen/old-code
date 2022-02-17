//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Description : Helper class for LN models to get data out of market cache
//
//   Author      : Jose Hilera
//
//   Date        : 13 May 2005
//
//----------------------------------------------------------------------------

#ifndef MARKETDATAFETCHERLNSPLINE_HPP
#define MARKETDATAFETCHERLNSPLINE_HPP

#include "edginc/MarketDataFetcherLN.hpp"
#include <string>

using namespace std;   // string

DRLIB_BEGIN_NAMESPACE

/** Helper class for LN models to get data out of market cache */        
class MARKET_DLL MarketDataFetcherLNSpline : public MarketDataFetcherLN {
public:
    MarketDataFetcherLNSpline(const string& volType);

    /** Override the default implementation of modifyMarketData so that it
     * invokes AssetUtil::surfaceSplined */
    virtual MarketObjectSP modifyMarketData(const IModel*         model,
                                            const MarketData*     market,
                                            const CClassConstSP&  clazz,    
                                            const MarketObjectSP& mo) const;

protected:
    MarketDataFetcherLNSpline();
    MarketDataFetcherLNSpline(const MarketDataFetcherLNSpline& rhs);
    MarketDataFetcherLNSpline& operator=(const MarketDataFetcherLNSpline& rhs);
};

typedef smartConstPtr<MarketDataFetcherLNSpline> MarketDataFetcherLNSplineConstSP;
typedef smartPtr<MarketDataFetcherLNSpline> MarketDataFetcherLNSplineSP;

DRLIB_END_NAMESPACE
#endif
