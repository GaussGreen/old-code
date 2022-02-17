//	Aggregators.h : Declaration of the CAggregators
//
/////////////////////////////////////////////////////////////////////////////

#ifndef __AGGREGATORS_H_
#define __AGGREGATORS_H_

#include "resource.h"
#include <comobjectcollection.h>

class ATL_NO_VTABLE CAggregators : public IDispatchImpl<IAggregatable, &IID_IAggregatable, &LIBID_SiriusRisk>, public CComObjectCollection<IAggregatable, CAggregators, IAggregators, &CLSID_Aggregators, &IID_IAggregators, &LIBID_SiriusRisk, IDR_AGGREGATORS>
{
	STDMETHOD(Aggregate)(/*[in]*/ BSTR Name, /*[in]*/ IResults* pResults);

	BEGIN_COM_MAP(CAggregators)	// This overrides the default in CComObjectCollection<>
		COM_INTERFACE_ENTRY(IAggregators)
		COM_INTERFACE_ENTRY(ISupportErrorInfo)
		COM_INTERFACE_ENTRY2(IDispatch, IAggregators)
		COM_INTERFACE_ENTRY(IAggregatable)		
	END_COM_MAP()
};

#endif
