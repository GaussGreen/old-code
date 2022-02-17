//	handle_marketdatacollection.h : Declaration of the CMarketDataCollection
//
/////////////////////////////////////////////////////////////////////////////

#ifndef __MARKETDATACOLLECTION_H_
#define __MARKETDATACOLLECTION_H_

#include "resource.h"
#include "siriuscomobjectcollection.h"

class ATL_NO_VTABLE CMarketDataCollection : public CSiriusComObjectCollection<IMarketData, CMarketDataCollection, IMarketDataCollection, &CLSID_MarketDataCollection, &IID_IMarketDataCollection, &LIBID_Sirius, &CLSID_MarketData, IDR_MARKETDATACOLLECTION>
{
};

#endif
