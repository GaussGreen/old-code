//----------------------------------------------------------------------------
//
//   Description : Demo market data fetcher
//
//   Author      : Jay Z Wang
//
//----------------------------------------------------------------------------

#ifndef EDR_MARKETDATAFETCHERDEMO_HPP
#define EDR_MARKETDATAFETCHERDEMO_HPP

#include "edginc/MarketDataFetcher.hpp"

DRLIB_BEGIN_NAMESPACE

/** This class doesn't need to exist any more - could be replaced with
    a static constructor to MarketDataFetcher. */        
class MCARLO_DLL MarketDataFetcherDemo: public MarketDataFetcher {
public:
    ~MarketDataFetcherDemo();

    MarketDataFetcherDemo();
 
private:
};

typedef smartConstPtr<MarketDataFetcherDemo> MarketDataFetcherDemoConstSP;
typedef smartPtr<MarketDataFetcherDemo> MarketDataFetcherDemoSP;

DRLIB_END_NAMESPACE
#endif
