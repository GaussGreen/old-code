// Gamma5.h : Declaration of the CGamma5

#ifndef __GAMMA5_H_
#define __GAMMA5_H_

#include "resource.h"       // main symbols
#include "risk.h"
/////////////////////////////////////////////////////////////////////////////
// CGamma5
class ATL_NO_VTABLE CGamma5 : 
	public CComObjectRootEx<CComSingleThreadModel>,
	public CComCoClass<CGamma5, &CLSID_Gamma5>,
	public ISupportErrorInfo,
	public IDispatchImpl<IGamma5, &IID_IGamma5, &LIBID_SiriusRisk>,
	public IDispatchImpl<IGreek, &IID_IGreek, &LIBID_SiriusRisk>,
	public CRisk

{
public:
	CGamma5(){}
	HRESULT								FinalConstruct();
	DECLARE_REGISTRY_RESOURCEID(IDR_GAMMA5)
	DECLARE_NOT_AGGREGATABLE(CGamma5)
	DECLARE_PROTECT_FINAL_CONSTRUCT()
	BEGIN_COM_MAP(CGamma5)
		COM_INTERFACE_ENTRY(IGamma5)
		COM_INTERFACE_ENTRY(ISupportErrorInfo)
		COM_INTERFACE_ENTRY2(IDispatch, IGamma5)
		COM_INTERFACE_ENTRY(IGreek)
	END_COM_MAP()

	STDMETHOD(Calculate)(/*[in]*/ IEvaluatable* Portfolio, /*[in]*/ IResult* pVal);
	STDMETHOD(InterfaceSupportsErrorInfo)(REFIID riid);
	STDMETHOD(get_Shift)(/*[out, retval]*/ double* pVal);
	STDMETHOD(put_Shift)(/*[in]*/ double newVal);
		
protected:
	
};

#endif //__GAMMA5_H_
