//----------------------------------------------------------------------------
//
//   Group       : GCCG Derivatives Research
//
//   Filename    : MarketDataFetcherE3F.hpp
//
//   Description : Helper class for E3F model to get data out of market cache
//
//   Author      : Spyridon Schismenos
//
//   Date        : 5 July 2005
//
//----------------------------------------------------------------------------

#ifndef MARKETDATAFETCHERE3F_HPP
#define MARKETDATAFETCHERE3F_HPP

#include "edginc/MarketDataFetcher.hpp"

DRLIB_BEGIN_NAMESPACE

/** Helper class for E3F Monte Carlo model to get data out of market cache.
    Unclear why this class exists */ 
class MCARLO_DLL MarketDataFetcherE3F : public MarketDataFetcher {
public:
    MarketDataFetcherE3F();
    MarketDataFetcherE3F(const string& volType);
};

typedef smartConstPtr<MarketDataFetcherE3F> MarketDataFetcherE3FConstSP;
typedef smartPtr<MarketDataFetcherE3F> MarketDataFetcherE3FSP;

DRLIB_END_NAMESPACE
#endif
