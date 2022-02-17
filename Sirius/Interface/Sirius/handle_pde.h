//	handle_pde.h : Declaration of the CPde
//
/////////////////////////////////////////////////////////////////////////////

#ifndef __PDE_H_
#define __PDE_H_

#include "resource.h"
#include "MlEqHandles.h"
#include "MlEqPde.h"




class ATL_NO_VTABLE CPde : 
	public CComObjectRootEx<CComSingleThreadModel>,
	public CComCoClass<CPde, &CLSID_Pde>,
	public ISupportErrorInfo,
	public IDispatchImpl<IPde, &IID_IPde, &LIBID_Sirius>
{
public:	
	CPde() : m_h(new MlEqPdeDriver){}
	DECLARE_REGISTRY_RESOURCEID(IDR_PDE)
	DECLARE_NOT_AGGREGATABLE(CPde)
	DECLARE_PROTECT_FINAL_CONSTRUCT()
	BEGIN_COM_MAP(CPde)
		COM_INTERFACE_ENTRY(IPde)
		COM_INTERFACE_ENTRY(IDispatch)
		COM_INTERFACE_ENTRY(ISupportErrorInfo)
	END_COM_MAP()
	associate_analytic_object(CPde, MlEqPdeDriver);	

	STDMETHOD(InterfaceSupportsErrorInfo)(REFIID riid);
	STDMETHOD(get_Value)(/*[out, retval]*/ VARIANT *pVal);
	STDMETHOD(put_Value)(/*[in]*/ VARIANT newVal);
};

#endif
