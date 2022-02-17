//	handle_dateschedule.h : Declaration of the CDateSchedule
//
/////////////////////////////////////////////////////////////////////////////

#ifndef __DATESCHEDULE_H_
#define __DATESCHEDULE_H_

#include "resource.h"
#include "mleqdateschedule.h"

class ATL_NO_VTABLE CDateSchedule : 
	public CComObjectRootEx<CComSingleThreadModel>,
	public CComCoClass<CDateSchedule, &CLSID_DateSchedule>,
	public ISupportErrorInfo,
	public IDispatchImpl<IDateSchedule, &IID_IDateSchedule, &LIBID_Sirius>
{
public:	
	CDateSchedule() : m_h(new MlEqDateSchedule){}
	DECLARE_REGISTRY_RESOURCEID(IDR_DATESCHEDULE)
	DECLARE_NOT_AGGREGATABLE(CDateSchedule)
	DECLARE_PROTECT_FINAL_CONSTRUCT()
	BEGIN_COM_MAP(CDateSchedule)
		COM_INTERFACE_ENTRY(IDateSchedule)
		COM_INTERFACE_ENTRY(IDispatch)
		COM_INTERFACE_ENTRY(ISupportErrorInfo)
	END_COM_MAP()
	associate_analytic_object(CDateSchedule, MlEqDateSchedule);
	HRESULT								FinalConstruct(void);

	STDMETHOD(InterfaceSupportsErrorInfo)(REFIID riid);
	STDMETHOD(get_Value)(/*[out, retval]*/ VARIANT *pVal);
	STDMETHOD(put_Value)(/*[in]*/ VARIANT newVal);
	STDMETHOD(get_NumberAt)(/*[in]*/ DATE Date, /*[out, retval]*/ double *pVal);
	STDMETHOD(put_NumberAt)(/*[in]*/ DATE Date, /*[in]*/ double newVal);
	STDMETHOD(get_ValueAt)(DATE Date, /*[out, retval]*/ VARIANT *pVal);
	STDMETHOD(put_ValueAt)(DATE Date, /*[in]*/ VARIANT newVal);
};

#endif
