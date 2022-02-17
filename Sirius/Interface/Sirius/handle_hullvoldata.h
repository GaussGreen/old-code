//	handle_hullvoldata.h : Declaration of the CHullVolData
//
/////////////////////////////////////////////////////////////////////////////

#ifndef __HULLVOLDATA_H_
#define __HULLVOLDATA_H_

#include "resource.h"
#include "MlEqObjects.h"

class ATL_NO_VTABLE CHullVolData : 
	public CComObjectRootEx<CComSingleThreadModel>,
	public CComCoClass<CHullVolData, &CLSID_HullVolData>,
	public ISupportErrorInfo,
	public IDispatchImpl<IHullVolData, &IID_IHullVolData, &LIBID_Sirius>	
{
public:
	MLHullVolDataHandle					m_h;	
	CHullVolData() : m_h(new MLHullVolData){}
	DECLARE_REGISTRY_RESOURCEID(IDR_HULLVOLDATA)
	DECLARE_NOT_AGGREGATABLE(CHullVolData)
	DECLARE_PROTECT_FINAL_CONSTRUCT()
	BEGIN_COM_MAP(CHullVolData)
		COM_INTERFACE_ENTRY(IHullVolData)
		COM_INTERFACE_ENTRY(IDispatch)
		COM_INTERFACE_ENTRY(ISupportErrorInfo)
	END_COM_MAP()
	HRESULT								FinalConstruct(void);

	STDMETHOD(InterfaceSupportsErrorInfo)(REFIID riid);
	STDMETHOD(get_Value)(/*[out, retval]*/ VARIANT *pVal);
	STDMETHOD(put_Value)(/*[in]*/ VARIANT newVal);

protected:
	std::vector<CAdapt<CComQIPtr<IStrikes> > >	m_aspStrike;
	CParameterMap								m_pmHullCoefficients;
};

#endif
