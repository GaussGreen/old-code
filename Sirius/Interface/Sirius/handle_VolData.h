//	handle_voldata.h : Declaration of the CVolData
//
/////////////////////////////////////////////////////////////////////////////

#ifndef __VOLDATA_H_
#define __VOLDATA_H_

#include "resource.h"
#include "MlEqObjects.h"

class ATL_NO_VTABLE CVolData : 
	public CComObjectRootEx<CComSingleThreadModel>,
	public CComCoClass<CVolData, &CLSID_VolData>,
	public ISupportErrorInfo,
	public IDispatchImpl<IVolData, &IID_IVolData, &LIBID_Sirius>
{
public:
	MlEqVolDataHandle					m_h;
	
	CVolData() : m_h(new MlEqVolData) {}
	DECLARE_REGISTRY_RESOURCEID(IDR_VOLDATA)
	DECLARE_NOT_AGGREGATABLE(CVolData)
	DECLARE_PROTECT_FINAL_CONSTRUCT()
	BEGIN_COM_MAP(CVolData)
		COM_INTERFACE_ENTRY(IVolData)
		COM_INTERFACE_ENTRY(IDispatch)
		COM_INTERFACE_ENTRY(ISupportErrorInfo)
	END_COM_MAP()
	HRESULT								FinalConstruct(void);

	STDMETHOD(InterfaceSupportsErrorInfo)(REFIID riid);
	STDMETHOD(get_Value)(/*[out, retval]*/ VARIANT* pVal);
	STDMETHOD(put_Value)(/*[in]*/ VARIANT newVal);	
	
protected:
	std::vector<CAdapt<CComQIPtr<IStrikes> > >		m_aspStrike;
	std::vector<CAdapt<CComQIPtr<IInterpolator> > >	m_aspInterpolator;
};

#endif
