//	handle_array.h : Declaration of the CArray
//
/////////////////////////////////////////////////////////////////////////////

#ifndef __ARRAY_H_
#define __ARRAY_H_

#include "resource.h"
#include "mleqarray.h"					// for m_h

class ATL_NO_VTABLE CArray : 
	public CComObjectRootEx<CComSingleThreadModel>,
	public CComCoClass<CArray, &CLSID_Array>,
	public ISupportErrorInfo,
	public IDispatchImpl<IArray, &IID_IArray, &LIBID_Sirius>
{
public:
	CArray() : m_h(new MlEqArray){}
	DECLARE_REGISTRY_RESOURCEID(IDR_ARRAY)
	DECLARE_NOT_AGGREGATABLE(CArray)
	DECLARE_PROTECT_FINAL_CONSTRUCT()
	BEGIN_COM_MAP(CArray)
		COM_INTERFACE_ENTRY(IArray)
		COM_INTERFACE_ENTRY(IDispatch)
		COM_INTERFACE_ENTRY(ISupportErrorInfo)
	END_COM_MAP()
	associate_analytic_object(CArray, MlEqArray);

	STDMETHOD(InterfaceSupportsErrorInfo)(REFIID riid);
	STDMETHOD(get_Value)(/*[out, retval]*/ VARIANT* pVal);
	STDMETHOD(put_Value)(/*[in]*/ VARIANT newVal);
	STDMETHOD(get_Element)(long Index, /*[out, retval]*/ double *pVal);
	STDMETHOD(put_Element)(long Index, /*[in]*/ double newVal);
	STDMETHOD(get_Size)(/*[out, retval]*/ long *pVal);
	STDMETHOD(put_Size)(/*[in]*/ long newVal);
};

#endif
