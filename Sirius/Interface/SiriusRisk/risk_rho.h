// risk_rho.h : Declaration of the CRHo

#ifndef __RHO_H_
#define __RHO_H_

#include "resource.h"       // main symbols
#include "risk.h"

/////////////////////////////////////////////////////////////////////////////
// CRHo
class ATL_NO_VTABLE CRho : 
	public CComObjectRootEx<CComSingleThreadModel>,
	public CComCoClass<CRho, &CLSID_Rho>,
	public ISupportErrorInfo,
	public IDispatchImpl<IRho, &IID_IRho, &LIBID_SiriusRisk>,
	public IDispatchImpl<IGreek, &IID_IGreek, &LIBID_SiriusRisk>,
	public CRisk

{
public:
	CRho()
	{
	}

DECLARE_REGISTRY_RESOURCEID(IDR_RHO)

DECLARE_PROTECT_FINAL_CONSTRUCT()

BEGIN_COM_MAP(CRho)
	COM_INTERFACE_ENTRY(IRho)
	COM_INTERFACE_ENTRY(ISupportErrorInfo)
	COM_INTERFACE_ENTRY2(IDispatch, IRho)
	COM_INTERFACE_ENTRY(IGreek)
END_COM_MAP()

	HRESULT FinalConstruct();
	STDMETHOD(Calculate)(/*[in]*/ IEvaluatable* Portfolio, /*[in]*/ IResult* pVal);
	STDMETHOD(InterfaceSupportsErrorInfo)(REFIID riid);
	STDMETHOD(get_Shift)(/*[out, retval]*/ double* pVal);
	STDMETHOD(put_Shift)(/*[in]*/ double newVal);
		
protected:
	double								m_fShift;
};

#endif //__RHO_H_
