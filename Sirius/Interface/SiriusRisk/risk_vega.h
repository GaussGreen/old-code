//	risk_vega.h : Declaration of the CVega
//
/////////////////////////////////////////////////////////////////////////////

#ifndef __VEGA_H_
#define __VEGA_H_

#include "resource.h"
#include "risk.h"

class ATL_NO_VTABLE CVega : 
	public CComObjectRootEx<CComSingleThreadModel>,
	public CComCoClass<CVega, &CLSID_Vega>,
	public ISupportErrorInfo,
	public IDispatchImpl<IVega, &IID_IVega, &LIBID_SiriusRisk>,
	public IDispatchImpl<IGreek, &IID_IGreek, &LIBID_SiriusRisk>,
	public IDispatchImpl<IPartialGreek, &IID_IPartialGreek, &LIBID_SiriusRisk>,
	public CRisk
{
public:
	CVega(){}
	HRESULT								FinalConstruct();
	DECLARE_REGISTRY_RESOURCEID(IDR_VEGA)
	DECLARE_NOT_AGGREGATABLE(CVega)
	DECLARE_PROTECT_FINAL_CONSTRUCT()
	BEGIN_COM_MAP(CVega)
		COM_INTERFACE_ENTRY(IVega)
		COM_INTERFACE_ENTRY(ISupportErrorInfo)
		COM_INTERFACE_ENTRY2(IDispatch, IPrice)
		COM_INTERFACE_ENTRY(IGreek)	
		COM_INTERFACE_ENTRY(IPartialGreek)
	END_COM_MAP()

	STDMETHOD(Calculate)(/*[in]*/ IEvaluatable* Portfolio, /*[in]*/ IResult* pVal);
	STDMETHOD(InterfaceSupportsErrorInfo)(REFIID riid);
	STDMETHOD(get_Shift)(/*[out, retval]*/ double* pVal);
	STDMETHOD(put_Shift)(/*[in]*/ double newVal);
	
	STDMETHOD(get_CalculatePartialGreek)(/*[out, retval] */ VARIANT_BOOL* pVal);
	STDMETHOD(put_CalculatePartialGreek)(/*[in]*/ VARIANT_BOOL newVal);
	STDMETHOD(GetCalculatedPartialGreekResults)(/*[in]*/ IResult* Result);

protected:
	double								m_fShift;
	bool								m_bPartialGreek; 
	std::map<std::string,double>		m_mapPartialGreeks;
};

#endif
