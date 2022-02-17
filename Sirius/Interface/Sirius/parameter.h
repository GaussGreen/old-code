//	parameter.h : Declaration of the CParameter
//
//////////////////////////////////////////////////////////////////////

#ifndef __PARAMETER_H_
#define __PARAMETER_H_

#include "resource.h"

class ATL_NO_VTABLE CParameter : 
	public CComObjectRootEx<CComSingleThreadModel>,
	public CComCoClass<CParameter, &CLSID_Parameter>,
	public ISupportErrorInfo,
	public IDispatchImpl<IParameter, &IID_IParameter, &LIBID_Sirius>,
	public CParameterMap
{
public:		
	CParameter(){}
	DECLARE_REGISTRY_RESOURCEID(IDR_PARAMETER)
	DECLARE_NOT_AGGREGATABLE(CParameter)
	DECLARE_PROTECT_FINAL_CONSTRUCT()
	BEGIN_COM_MAP(CParameter)
		COM_INTERFACE_ENTRY(IParameter)
		COM_INTERFACE_ENTRY(IDispatch)
		COM_INTERFACE_ENTRY(ISupportErrorInfo)
	END_COM_MAP()
	
	STDMETHOD(InterfaceSupportsErrorInfo)(REFIID riid);
	STDMETHOD(get_Value)(/*[out, retval]*/ VARIANT *pVal);
	STDMETHOD(put_Value)(/*[in]*/ VARIANT newVal);
	STDMETHOD(get_Name)(/*[out, retval]*/ BSTR *pVal);
	STDMETHOD(put_Name)(/*[in]*/ BSTR newVal);
	STDMETHOD(GetColumn)(/*[in]*/ long Column, /*[out, retval]*/ IParameter** pVal);
	STDMETHOD(GetRow)(/*[in]*/ long Row, /*[out, retval]*/ IParameter** pVal);
	STDMETHOD(get_Element)(/*[in]*/ long Row, /*[in]*/ long Column, /*[out, retval]*/ VARIANT *pVal);
	STDMETHOD(put_Element)(/*[in]*/ long Row, /*[in]*/ long Column, /*[in]*/ VARIANT newVal);
	STDMETHOD(CopyToClipboard)();
	STDMETHOD(RemoveRow)(/*[in]*/ long Row);
	STDMETHOD(InsertRow)(/*[in]*/ long At);
	STDMETHOD(InsertColumn)(/*[in]*/ long At);
	STDMETHOD(RemoveColumn)(/*[in]*/ long Column);
	STDMETHOD(get_Columns)(/*[out, retval]*/ long *pVal);
	STDMETHOD(put_Columns)(/*[in]*/ long newVal);
	STDMETHOD(get_Rows)(/*[out, retval]*/ long *pVal);
	STDMETHOD(put_Rows)(/*[in]*/ long newVal);
	STDMETHOD(Transpose)();
	STDMETHOD(get_AutoGrow)(/*[out, retval]*/ VARIANT_BOOL *pVal);
	STDMETHOD(put_AutoGrow)(/*[in]*/ VARIANT_BOOL newVal);
	STDMETHOD(AddToEnd)(/*[in]*/ IParameter* Val);
	STDMETHOD(AddToRHS)(/*[in]*/ IParameter* Val);
	STDMETHOD(IsNumeric)(/*[out, retval]*/ VARIANT_BOOL* pVal);
	STDMETHOD(ReplaceBlanksWithSpaces)(void);
	STDMETHOD(IsBlank)(/*[out, retval]*/ VARIANT_BOOL* pVal);

protected:
	estring								m_szName;
};

#endif
