//	risk_delta.h : Declaration of the CDelta
//
/////////////////////////////////////////////////////////////////////////////

#ifndef __DELTA_H_
#define __DELTA_H_

#include "resource.h"
#include "risk.h"


class ATL_NO_VTABLE CDelta : 
	public CComObjectRootEx<CComSingleThreadModel>,
	public CComCoClass<CDelta, &CLSID_Delta>,
	public ISupportErrorInfo,
	public IDispatchImpl<IDelta, &IID_IDelta, &LIBID_SiriusRisk>,
	public IDispatchImpl<IGreek, &IID_IGreek, &LIBID_SiriusRisk>,
	public IDispatchImpl<IPartialGreek, &IID_IPartialGreek, &LIBID_SiriusRisk>,
	public CRisk
	
{
public:	
	CDelta(){}
	HRESULT								FinalConstruct();
	DECLARE_REGISTRY_RESOURCEID(IDR_DELTA)
	DECLARE_NOT_AGGREGATABLE(CDelta)
	DECLARE_PROTECT_FINAL_CONSTRUCT()
	BEGIN_COM_MAP(CDelta)
		COM_INTERFACE_ENTRY(IDelta)
		COM_INTERFACE_ENTRY(ISupportErrorInfo)
		COM_INTERFACE_ENTRY2(IDispatch, IDelta)
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
