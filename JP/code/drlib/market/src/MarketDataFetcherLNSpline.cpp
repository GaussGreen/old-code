//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Description : Helper class for LN models to get data out of market cache.
//                 Implements a surface spline when modifying the market data
//
//   Author      : Jose Hilera
//
//   Date        : 13 May 2005
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/AssetUtil.hpp"
#include "edginc/MarketDataFetcherLNSpline.hpp"


DRLIB_BEGIN_NAMESPACE

MarketDataFetcherLNSpline::MarketDataFetcherLNSpline(const string& volType): 
    MarketDataFetcherLN(volType) {}

MarketDataFetcherLNSpline::MarketDataFetcherLNSpline(): 
    MarketDataFetcherLN() {}

/** Gives MCImplied the chance to Spline the VolSurface */
MarketObjectSP MarketDataFetcherLNSpline::modifyMarketData(
    const IModel*         model,
    const MarketData*     market,
    const CClassConstSP&  clazz,      // what type was originally requested
    const MarketObjectSP& mo) const   /* what GetMarket returned or what was
                                       * "inline" already */
{
    try {
        return AssetUtil::surfaceSplined(market, model, mo);
    } 
    catch(exception& e) {
        throw ModelException(e, "MarketDataFetcherLNSpline::modifyMarketData");
    }
}


DRLIB_END_NAMESPACE
