/**
 * @file MarketData_forward.hpp
 */

#include "edginc/smartPtr.hpp"

#ifndef DRLIB_MarketData_forward_H
#define DRLIB_MarketData_forward_H

DRLIB_BEGIN_NAMESPACE

class MarketData;
typedef smartConstPtr<MarketData> CMarketDataConstSP;
typedef smartPtr<MarketData> CMarketDataSP;
#ifndef QLIB_MARKETDATA_CPP
EXTERN_TEMPLATE(class RISKMGR_DLL_SP smartConstPtr<MarketData>);
EXTERN_TEMPLATE(class RISKMGR_DLL_SP smartPtr<MarketData>);
EXTERN_TEMPLATE(IObjectSP RISKMGR_DLL FieldGetSmartPtr<CMarketDataSP>(CMarketDataSP* t));
EXTERN_TEMPLATE(void RISKMGR_DLL FieldSetSmartPtr<CMarketDataSP>(CMarketDataSP* t,
                                                     IObjectSP o));
#else
INSTANTIATE_TEMPLATE(class RISKMGR_DLL smartConstPtr<MarketData>);
INSTANTIATE_TEMPLATE(class RISKMGR_DLL smartPtr<MarketData>);
INSTANTIATE_TEMPLATE(IObjectSP RISKMGR_DLL FieldGetSmartPtr<CMarketDataSP>(
                         CMarketDataSP* t));
INSTANTIATE_TEMPLATE(void RISKMGR_DLL FieldSetSmartPtr<CMarketDataSP>(CMarketDataSP* t,
                                                          IObjectSP o));
#endif
typedef CMarketDataConstSP MarketDataConstSP;
typedef CMarketDataSP MarketDataSP;

DRLIB_END_NAMESPACE

// Local Variables: ***
// c-basic-offset:4 ***
// End: ***

#endif
