//	risk_price.h : Declaration of the CPrice
//
/////////////////////////////////////////////////////////////////////////////

#ifndef __PRICE_H_
#define __PRICE_H_

#include "resource.h"
#include "risk.h"

class ATL_NO_VTABLE CPrice : 
	public CComObjectRootEx<CComSingleThreadModel>,
	public CComCoClass<CPrice, &CLSID_Price>,
	public ISupportErrorInfo,
	public IDispatchImpl<IPrice, &IID_IPrice, &LIBID_SiriusRisk>,
	public IDispatchImpl<IGreek, &IID_IGreek, &LIBID_SiriusRisk>,
	public CRisk
{
public:
	CPrice(){}

	DECLARE_REGISTRY_RESOURCEID(IDR_PRICE)
	DECLARE_NOT_AGGREGATABLE(CPrice)
	DECLARE_PROTECT_FINAL_CONSTRUCT()
	BEGIN_COM_MAP(CPrice)
		COM_INTERFACE_ENTRY(IPrice)
		COM_INTERFACE_ENTRY(ISupportErrorInfo)
		COM_INTERFACE_ENTRY2(IDispatch, IPrice)
		COM_INTERFACE_ENTRY(IGreek)	
	END_COM_MAP()

	STDMETHOD(Calculate)(/*[in]*/ IEvaluatable* Portfolio, /*[in*/ IResult* pVal);
	STDMETHOD(InterfaceSupportsErrorInfo)(REFIID riid);	
};

#endif
