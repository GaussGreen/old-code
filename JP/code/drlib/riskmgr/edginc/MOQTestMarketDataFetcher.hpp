//------------------------------------------------------------------------------
//
//   Group       : QR - Core Analytics 
//
//   Description : Market data fetcher for testing retrieval of market data using
//                 market data qualifiers
//
//   Author      : Andrew Greene 
//
//   Date        : 13 September 2006
//
//------------------------------------------------------------------------------

#ifndef MOQ_TEST_MARKET_DATA_FETCHER_HPP
#define MOQ_TEST_MARKET_DATA_FETCHER_HPP

#include "edginc/MarketDataFetcher.hpp"
#include "edginc/MarketObjectQualifierString.hpp"

DRLIB_BEGIN_NAMESPACE

class RISKMGR_DLL MOQTestMarketDataFetcher : public MarketDataFetcher
{
public:
    MOQTestMarketDataFetcher(string qualifierChoice);

    /** Resolve an array of multiply selected market data
        into a single chosen market datum */
    MarketObjectConstSP resolveMultiData(const MarketObjectArray& marketObjs,
                                         CClassConstSP            type) const;

private:
    MarketObjectQualifierString qualifierChoice;
};

DECLARE(MOQTestMarketDataFetcher);

DRLIB_END_NAMESPACE

#endif // MOQ_TEST_MARKET_DATA_FETCHER_HPP
