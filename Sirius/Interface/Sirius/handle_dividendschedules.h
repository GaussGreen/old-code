//	handle_dividendschedules.h : Declaration of the CDividendSchedules
//
/////////////////////////////////////////////////////////////////////////////

#ifndef __DIVIDENDSCHEDULES_H_
#define __DIVIDENDSCHEDULES_H_

#include "resource.h"
#include "comobjectcollectionserialisable.h"

class ATL_NO_VTABLE CDividendSchedules : public CComObjectCollectionSerialisable<IDividendSchedule, CDividendSchedules, IDividendSchedules, &CLSID_DividendSchedules, &IID_IDividendSchedules, &LIBID_Sirius, &CLSID_DividendSchedule, IDR_DIVIDENDSCHEDULES> 
{
public:
	static HRESULT						Load(const std::string& szDummy, DataSourceEnum ds, long nDate, CComPtr<IDividendSchedules>& spDividendSchedules);
	STDMETHOD(Reset)(void);
	STDMETHOD(Shift)(/*[in]*/ double Yield);
	STDMETHOD(Stick)(void);
};

#endif
