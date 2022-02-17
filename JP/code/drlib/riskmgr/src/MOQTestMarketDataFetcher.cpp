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

#include "edginc/config.hpp"

#include "edginc/MOQTestMarketDataFetcher.hpp"

DRLIB_BEGIN_NAMESPACE

MOQTestMarketDataFetcher::MOQTestMarketDataFetcher(string qualifierChoice):
    qualifierChoice(qualifierChoice)
{}

/** Resolve an array of multiply selected market data
    into a single chosen market datum */
MarketObjectConstSP MOQTestMarketDataFetcher::resolveMultiData(
    const MarketObjectArray& marketObjs,
    CClassConstSP            type) const
{
    // Select market data with specific market object qualifier string
    for (MarketObjectArray::const_iterator i = marketObjs.begin();
        i != marketObjs.end(); ++i)
    {
        if ((*i)->getMarketObjectQualifier()->equals(&qualifierChoice))
            return *i;
    }

    return MarketObjectConstSP();
}

DRLIB_END_NAMESPACE
