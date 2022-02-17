// risk_Theta.h : Declaration of the CTheta

#ifndef __THETA_H_
#define __THETA_H_

#include "resource.h"       // main symbols
#include "risk.h"
/////////////////////////////////////////////////////////////////////////////
// CTheta
class ATL_NO_VTABLE CTheta : 
	public CComObjectRootEx<CComSingleThreadModel>,
	public CComCoClass<CTheta, &CLSID_Theta>,
	public ISupportErrorInfo,
	public IDispatchImpl<ITheta, &IID_ITheta, &LIBID_SiriusRisk>,
	public IDispatchImpl<IGreek, &IID_IGreek, &LIBID_SiriusRisk>,
	public CRisk
{
public:
	CTheta(){}

DECLARE_REGISTRY_RESOURCEID(IDR_THETA)
DECLARE_NOT_AGGREGATABLE(CTheta)
	HRESULT								FinalConstruct();
DECLARE_PROTECT_FINAL_CONSTRUCT()

BEGIN_COM_MAP(CTheta)
		COM_INTERFACE_ENTRY(ITheta)
		COM_INTERFACE_ENTRY(ISupportErrorInfo)
		COM_INTERFACE_ENTRY2(IDispatch, ITheta)
		COM_INTERFACE_ENTRY(IGreek)		
	END_COM_MAP()
// ITheta
public:
	STDMETHOD(Calculate)(/*[in]*/ IEvaluatable* Portfolio, /*[in]*/ IResult* pVal);
	STDMETHOD(InterfaceSupportsErrorInfo)(REFIID riid);
	STDMETHOD(get_Shift)(/*[out, retval]*/ double* pVal);
	STDMETHOD(put_Shift)(/*[in]*/ double newVal);

protected:
	double								m_fShift;
};

#endif //__THETA_H_
