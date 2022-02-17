//	HullAndWhite.h : Declaration of the CHullAndWhite
//
/////////////////////////////////////////////////////////////////////////////

#ifndef __HULLANDWHITE_H_
#define __HULLANDWHITE_H_

#include "resource.h"
#include "montecarlo.h"

class ATL_NO_VTABLE CHullAndWhite : 
	public CComObjectRootEx<CComSingleThreadModel>,
	public CComCoClass<CHullAndWhite, &CLSID_HullAndWhite>,
	public ISupportErrorInfo,
	public IDispatchImpl<IHullAndWhite, &IID_IHullAndWhite, &LIBID_Sirius>	
{
public:
	CHullAndWhite() : m_h(new MlEqHullAndWhite(NULL)){}
	DECLARE_REGISTRY_RESOURCEID(IDR_HULLANDWHITE)
	DECLARE_NOT_AGGREGATABLE(CHullAndWhite)
	DECLARE_PROTECT_FINAL_CONSTRUCT()
	BEGIN_COM_MAP(CHullAndWhite)
		COM_INTERFACE_ENTRY(IHullAndWhite)
		COM_INTERFACE_ENTRY(IDispatch)
		COM_INTERFACE_ENTRY(ISupportErrorInfo)
	END_COM_MAP()
	associate_analytic_object(CHullAndWhite, MlEqHullAndWhite);

	STDMETHOD(InterfaceSupportsErrorInfo)(REFIID riid);
	STDMETHOD(get_Value)(/*[out, retval]*/ VARIANT *pVal);
	STDMETHOD(put_Value)(/*[in]*/ VARIANT newVal);
};

#endif
