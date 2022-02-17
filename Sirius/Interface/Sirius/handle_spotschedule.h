//	handle_spotschedule.h : Declaration of the CSpotSchedule
//
/////////////////////////////////////////////////////////////////////////////

#ifndef __SPOTSCHEDULE_H_
#define __SPOTSCHEDULE_H_

#include "resource.h"
#include "mleqspotschedule.h"			// for m_h

class ATL_NO_VTABLE CSpotSchedule : 
	public CComObjectRootEx<CComSingleThreadModel>,
	public CComCoClass<CSpotSchedule, &CLSID_SpotSchedule>,
	public ISupportErrorInfo,
	public IDispatchImpl<ISpotSchedule, &IID_ISpotSchedule, &LIBID_Sirius>
{
public:		
	CSpotSchedule()	: m_h(new MlEqSpotSchedule(_Module.GetUseYesterdaysClose)) {}
	DECLARE_REGISTRY_RESOURCEID(IDR_SPOTSCHEDULE)
	DECLARE_NOT_AGGREGATABLE(CSpotSchedule)
	DECLARE_PROTECT_FINAL_CONSTRUCT()
	BEGIN_COM_MAP(CSpotSchedule)
		COM_INTERFACE_ENTRY(ISpotSchedule)
		COM_INTERFACE_ENTRY(IDispatch)
		COM_INTERFACE_ENTRY(ISupportErrorInfo)
	END_COM_MAP()
	associate_analytic_object(CSpotSchedule, MlEqSpotSchedule);
	declare_serialisable;

	HRESULT								FinalConstruct(void);
	static HRESULT						Load(const std::string& szIdentifier, DataSourceEnum ds, long nDate, CComPtr<ISpotSchedule>& spSpotSchedule);

	STDMETHOD(InterfaceSupportsErrorInfo)(REFIID riid);
	STDMETHOD(get_Value)(/*[out, retval]*/ VARIANT *pVal);
	STDMETHOD(put_Value)(/*[in]*/ VARIANT newVal);
	STDMETHOD(get_ValueAt)(/*[in]*/ DATE Date, /*[out, retval]*/ double *pVal);
	STDMETHOD(put_ValueAt)(/*[in]*/ DATE Date, /*[in]*/ double newVal);
	STDMETHOD(get_LastDate)(/*[out, retval]*/ DATE* pVal);
	STDMETHOD(get_FirstDate)(/*[out, retval]*/ DATE* pVal);	
	STDMETHOD(Clear)();
	STDMETHOD(Save)(/*[out, retval]*/ BSTR* pVal);
	STDMETHOD(get_Values)(/*[out, retval]*/ VARIANT *pVal);
	STDMETHOD(get_Dates)(/*[out, retval]*/ VARIANT *pVal);
	STDMETHOD(get_LastValue)(/*[out, retval]*/ double *pVal);
	STDMETHOD(get_FirstValue)(/*[out, retval]*/ double *pVal);
	STDMETHOD(put_LastValue)(/*[in]*/ double newVal);
	STDMETHOD(put_FirstValue)(/*[in]*/ double newVal);
	STDMETHOD(get_Name)(/*[out, retval]*/ BSTR* pVal);	
	STDMETHOD(Multiply)(/*[in]*/ double Amount);
	declare_member_variable_rekey(DATE, Date);
	declare_member_variable_rekey(DataSourceEnum, DataSource);
};

#endif
