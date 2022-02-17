//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : MarketDataFetcherLN.hpp
//
//   Description : Helper class for LN models to get data out of market cache
//
//   Author      : Andrew J Swain
//
//   Date        : 1 February 2002
//
//
//----------------------------------------------------------------------------

#ifndef MARKETDATAFETCHERLN_HPP
#define MARKETDATAFETCHERLN_HPP

#include "edginc/MarketDataFetcher.hpp"
#include <string>

using namespace std;   // string

DRLIB_BEGIN_NAMESPACE

/** Helper class for LN models to get data out of market cache */        
class MARKET_DLL MarketDataFetcherLN : public MarketDataFetcher {
public:
    MarketDataFetcherLN(const string& volType, const bool useCcyBasis = false);

    /** LN fetch method */
    virtual MarketObjectSP fetch(const MarketData*    market,
                                 const string&        name,
                                 const CClassConstSP& type,
                                 const IModel*        model) const;

protected:
    MarketDataFetcherLN();
    MarketDataFetcherLN(const MarketDataFetcherLN& rhs);
    MarketDataFetcherLN& operator=(const MarketDataFetcherLN& rhs);

    string volType;
    mutable CClassConstSP volClass; // not in registration
private:
    void initialise();
};

typedef smartConstPtr<MarketDataFetcherLN> MarketDataFetcherLNConstSP;
typedef smartPtr<MarketDataFetcherLN> MarketDataFetcherLNSP;
#ifndef QLIB_MARKETDATAFETCHERLN_CPP
EXTERN_TEMPLATE(class MARKET_DLL_SP smartPtr<MarketDataFetcherLN>);
#else
INSTANTIATE_TEMPLATE(class MARKET_DLL smartPtr<MarketDataFetcherLN>);
#endif

DRLIB_END_NAMESPACE
#endif
