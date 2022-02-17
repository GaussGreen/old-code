//----------------------------------------------------------------------------
//
//   Description : demo market data fetcher
//
//   Author      : Jay Z Wang
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/MarketDataFetcherDemo.hpp"
#include "edginc/IRVol.hpp"

DRLIB_BEGIN_NAMESPACE

MarketDataFetcherDemo::~MarketDataFetcherDemo(){}

MarketDataFetcherDemo::MarketDataFetcherDemo(){
    // if we are asked for IRVolBase (or derivatives) then choose
    // IRVol if that is compatible
    setRetrievalMode(IRVolBase::TYPE, true, IRVolCommon::TYPE);
}

DRLIB_END_NAMESPACE

