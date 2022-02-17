// risk_VegaSkew.h : Declaration of the CVegaSkew

#ifndef __VEGASKEW_H_
#define __VEGASKEW_H_

#include "resource.h"       // main symbols
#include "risk.h"
/////////////////////////////////////////////////////////////////////////////
// CVegaSkew
class ATL_NO_VTABLE CVegaSkew : 
	public CComObjectRootEx<CComSingleThreadModel>,
	public CComCoClass<CVegaSkew, &CLSID_VegaSkew>,
	public ISupportErrorInfo,
	public IDispatchImpl<IVegaSkew, &IID_IVegaSkew, &LIBID_SiriusRisk>,
	public IDispatchImpl<IGreek, &IID_IGreek, &LIBID_SiriusRisk>,
	public CRisk
{
public:
	CVegaSkew(){}
	HRESULT								FinalConstruct();
	DECLARE_REGISTRY_RESOURCEID(IDR_VEGASKEW)
	DECLARE_NOT_AGGREGATABLE(CVegaSkew)
	DECLARE_PROTECT_FINAL_CONSTRUCT()
	BEGIN_COM_MAP(CVegaSkew)
		COM_INTERFACE_ENTRY(IVegaSkew)
		COM_INTERFACE_ENTRY(ISupportErrorInfo)
		COM_INTERFACE_ENTRY2(IDispatch, IVegaSkew)
		COM_INTERFACE_ENTRY(IGreek)
	END_COM_MAP()

	STDMETHOD(Calculate)(/*[in]*/ IEvaluatable* Portfolio, /*[in]*/ IResult* pVal);
	STDMETHOD(InterfaceSupportsErrorInfo)(REFIID riid);
	STDMETHOD(get_ExpiryDate)(/*[out, retval]*/ DATE* pVal);
	STDMETHOD(put_ExpiryDate)(/*[in]*/ DATE newVal);
		
protected:
	DATE								m_dtExpiryDate;
};

#endif //__VEGASKEW_H_
