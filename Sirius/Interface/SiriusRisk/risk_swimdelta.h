// risk_SwimDelta.h : Declaration of the CSwimDelta

#ifndef __SWIMDELTA_H_
#define __SWIMDELTA_H_

#include "resource.h"       // main symbols
#include "risk.h"


/////////////////////////////////////////////////////////////////////////////
// CSwimDelta
class ATL_NO_VTABLE CSwimDelta : 
	public CComObjectRootEx<CComSingleThreadModel>,
	public CComCoClass<CSwimDelta, &CLSID_SwimDelta>,
	public ISupportErrorInfo,
	public IDispatchImpl<ISwimDelta, &IID_ISwimDelta, &LIBID_SiriusRisk>,
	public IDispatchImpl<IGreek, &IID_IGreek, &LIBID_SiriusRisk>,
	public CRisk

{
public:
	CSwimDelta()
	{
	}

	

	HRESULT								FinalConstruct();
	DECLARE_REGISTRY_RESOURCEID(IDR_SWIMDELTA)
	DECLARE_NOT_AGGREGATABLE(CSwimDelta)
	DECLARE_PROTECT_FINAL_CONSTRUCT()
	BEGIN_COM_MAP(CSwimDelta)
		COM_INTERFACE_ENTRY(ISwimDelta)
		COM_INTERFACE_ENTRY(ISupportErrorInfo)
		COM_INTERFACE_ENTRY2(IDispatch, ISwimDelta)
		COM_INTERFACE_ENTRY(IGreek)
	END_COM_MAP()

	STDMETHOD(Calculate)(/*[in]*/ IEvaluatable* Portfolio, /*[in]*/ IResult* pVal);
	STDMETHOD(InterfaceSupportsErrorInfo)(REFIID riid);
	STDMETHOD(get_Shift)(/*[out, retval]*/ double* pVal);
	STDMETHOD(put_Shift)(/*[in]*/ double newVal);
		
protected:
	double								m_fShift;
};

#endif //__SWIMDELTA_H_
