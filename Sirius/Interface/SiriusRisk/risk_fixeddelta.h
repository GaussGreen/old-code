// risk_FixedDelta.h : Declaration of the CFixedDelta

#ifndef __FIXEDDELTA_H_
#define __FIXEDDELTA_H_

#include "resource.h"       // main symbols
#include "risk.h"
/////////////////////////////////////////////////////////////////////////////
// CFixedDelta
class ATL_NO_VTABLE CFixedDelta : 
	public CComObjectRootEx<CComSingleThreadModel>,
	public CComCoClass<CFixedDelta, &CLSID_FixedDelta>,
	public ISupportErrorInfo,
	public IDispatchImpl<IFixedDelta, &IID_IFixedDelta, &LIBID_SiriusRisk>,
	public IDispatchImpl<IGreek, &IID_IGreek, &LIBID_SiriusRisk>,
	public CRisk

{
public:
	CFixedDelta(){}

	HRESULT								FinalConstruct();

	DECLARE_REGISTRY_RESOURCEID(IDR_FIXEDDELTA)
	DECLARE_NOT_AGGREGATABLE(CFixedDelta)
	DECLARE_PROTECT_FINAL_CONSTRUCT()
	
	BEGIN_COM_MAP(CFixedDelta)
		COM_INTERFACE_ENTRY(IFixedDelta)
		COM_INTERFACE_ENTRY(ISupportErrorInfo)
		COM_INTERFACE_ENTRY2(IDispatch, IFixedDelta)
		COM_INTERFACE_ENTRY(IGreek)
	END_COM_MAP()


	STDMETHOD(Calculate)(/*[in]*/ IEvaluatable* Portfolio, /*[in]*/ IResult* pVal);

// ISupportsErrorInfo
	STDMETHOD(InterfaceSupportsErrorInfo)(REFIID riid);
	STDMETHOD(get_Shift)(/*[out, retval]*/ double* pVal);
	STDMETHOD(put_Shift)(/*[in]*/ double newVal);
		
protected:
	double								m_fShift;
// IFixedDelta
public:
};

#endif //__FIXEDDELTA_H_
