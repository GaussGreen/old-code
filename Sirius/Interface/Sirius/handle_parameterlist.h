//	handle_parameterlist.h : Declaration of the CParameterList
//
/////////////////////////////////////////////////////////////////////////////

#ifndef __PARAMETERLIST_H_
#define __PARAMETERLIST_H_

#include "resource.h"       // main symbols
#include "mleqparameterlist.h"

class ATL_NO_VTABLE CParameterList : 
	public CComObjectRootEx<CComSingleThreadModel>,
	public CComCoClass<CParameterList, &CLSID_ParameterList>,
	public ISupportErrorInfo,
	public IDispatchImpl<IParameterList, &IID_IParameterList, &LIBID_Sirius>	
{
public:
	CParameterList() : m_h(new MlEqParameterList){}
	HRESULT								FinalConstruct(void);
	DECLARE_REGISTRY_RESOURCEID(IDR_PARAMETERLIST)	
	DECLARE_NOT_AGGREGATABLE(CParameterList)
	DECLARE_PROTECT_FINAL_CONSTRUCT()
	BEGIN_COM_MAP(CParameterList)
		COM_INTERFACE_ENTRY(IParameterList)
		COM_INTERFACE_ENTRY(IDispatch)
		COM_INTERFACE_ENTRY(ISupportErrorInfo)
	END_COM_MAP()
	associate_analytic_object(CParameterList, MlEqParameterList);

	STDMETHOD(InterfaceSupportsErrorInfo)(REFIID riid);
	STDMETHOD(get_Value)(/*[out, retval]*/ VARIANT *pVal);
	STDMETHOD(put_Value)(/*[in]*/ VARIANT newVal);
	STDMETHOD(AddValue)(/*[in]*/ BSTR Name, /*[in]*/ VARIANT Value);
	STDMETHOD(GetValue)(/*[in]*/ BSTR Name, /*[out, retval]*/ VARIANT* pVal);
	STDMETHOD(get_Count)(/*[out, retval]*/ long *pVal);	
};

#endif
