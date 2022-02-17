
#ifndef __ARPROPFITVOLDATA_H_
#define __ARPROPFITVOLDATA_H_

#include "resource.h"       // main symbols
#include "MlEqObjects.h"

class ATL_NO_VTABLE CARPropFitVolData : 
	public CComObjectRootEx<CComSingleThreadModel>,
	public CComCoClass<CARPropFitVolData, &CLSID_ARPropFitVolData>,
	public ISupportErrorInfo,
	public IDispatchImpl<IARPropFitVolData, &IID_IARPropFitVolData, &LIBID_Sirius>
{
public:
	MlARPropFitHandle					m_h;	
	CARPropFitVolData() : m_h(new MlARPropFit){}
	DECLARE_REGISTRY_RESOURCEID(IDR_ARPROPFITVOLDATA)
	DECLARE_NOT_AGGREGATABLE(CARPropFitVolData)
	DECLARE_PROTECT_FINAL_CONSTRUCT()

	BEGIN_COM_MAP(CARPropFitVolData)
		COM_INTERFACE_ENTRY(IARPropFitVolData)
		COM_INTERFACE_ENTRY(IDispatch)
		COM_INTERFACE_ENTRY(ISupportErrorInfo)
	END_COM_MAP()

	STDMETHOD(InterfaceSupportsErrorInfo)(REFIID riid);
	STDMETHOD(get_Value)(/*[out, retval]*/ VARIANT *pVal);
	STDMETHOD(put_Value)(/*[in]*/ VARIANT newVal);
};

#endif
