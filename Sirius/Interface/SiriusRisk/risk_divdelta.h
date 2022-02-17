// DivDelta.h : Declaration of the CDivDelta

#ifndef __DIVDELTA_H_
#define __DIVDELTA_H_

#include "resource.h"       // main symbols
#include "risk.h"
/////////////////////////////////////////////////////////////////////////////
// CDivDelta
class ATL_NO_VTABLE CDivDelta : 
	public CComObjectRootEx<CComSingleThreadModel>,
	public CComCoClass<CDivDelta, &CLSID_DivDelta>,
	public ISupportErrorInfo,
	public IDispatchImpl<IDivDelta, &IID_IDivDelta, &LIBID_SiriusRisk>,
	public IDispatchImpl<IGreek, &IID_IGreek, &LIBID_SiriusRisk>,
	public CRisk
{
public:	
	CDivDelta(){}
	HRESULT								FinalConstruct();
	DECLARE_REGISTRY_RESOURCEID(IDR_DIVDELTA)
	DECLARE_NOT_AGGREGATABLE(CDivDelta)
	DECLARE_PROTECT_FINAL_CONSTRUCT()
	BEGIN_COM_MAP(CDivDelta)
		COM_INTERFACE_ENTRY(IDivDelta)
		COM_INTERFACE_ENTRY(ISupportErrorInfo)
		COM_INTERFACE_ENTRY2(IDispatch, IDivDelta)
		COM_INTERFACE_ENTRY(IGreek)
	END_COM_MAP()

	STDMETHOD(Calculate)(/*[in]*/ IEvaluatable* Portfolio, /*[in]*/ IResult* pVal);
	STDMETHOD(InterfaceSupportsErrorInfo)(REFIID riid);
	STDMETHOD(get_Shift)(/*[out, retval]*/ double* pVal);
	STDMETHOD(put_Shift)(/*[in]*/ double newVal);
		
protected:
	double								m_fShift;
};

#endif //__DIVDELTA_H_
