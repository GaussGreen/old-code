//	risk_gamma.h : Declaration of the CGamma
//
/////////////////////////////////////////////////////////////////////////////

#ifndef __GAMMA_H_
#define __GAMMA_H_

#include "resource.h"
#include "risk.h"

class ATL_NO_VTABLE CGamma : 
	public CComObjectRootEx<CComSingleThreadModel>,
	public CComCoClass<CGamma, &CLSID_Gamma>,
	public ISupportErrorInfo,
	public IDispatchImpl<IGamma, &IID_IGamma, &LIBID_SiriusRisk>,
	public IDispatchImpl<IGreek, &IID_IGreek, &LIBID_SiriusRisk>,
	public CRisk
{
public:	
	CGamma(){}
	HRESULT								FinalConstruct();
	DECLARE_REGISTRY_RESOURCEID(IDR_GAMMA)
	DECLARE_NOT_AGGREGATABLE(CGamma)
	DECLARE_PROTECT_FINAL_CONSTRUCT()
	BEGIN_COM_MAP(CGamma)
		COM_INTERFACE_ENTRY(IGamma)
		COM_INTERFACE_ENTRY(ISupportErrorInfo)
		COM_INTERFACE_ENTRY2(IDispatch, IGamma)
		COM_INTERFACE_ENTRY(IGreek)		
	END_COM_MAP()
	
	STDMETHOD(Calculate)(/*[in]*/ IEvaluatable* Portfolio, /*[in]*/ IResult* pVal);
	STDMETHOD(InterfaceSupportsErrorInfo)(REFIID riid);
	STDMETHOD(get_Shift)(/*[out, retval]*/ double* pVal);
	STDMETHOD(put_Shift)(/*[in]*/ double newVal);

protected:
	double								m_fShift;
};

#endif
