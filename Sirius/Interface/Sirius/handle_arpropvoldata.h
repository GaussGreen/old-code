// ARPropVolData.h : Declaration of the CARPropVolData

#ifndef __ARPROPVOLDATA_H_
#define __ARPROPVOLDATA_H_

#include "resource.h"       // main symbols
#include "MlEqObjects.h"

class ATL_NO_VTABLE CARPropVolData : 
	public CComObjectRootEx<CComSingleThreadModel>,
	public CComCoClass<CARPropVolData, &CLSID_ARPropVolData>,
	public ISupportErrorInfo,
	public IDispatchImpl<IARPropVolData, &IID_IARPropVolData, &LIBID_Sirius>	
{
public:
	MlARPropVolDataHandle				m_h;	
	CARPropVolData() : m_h(new MlARPropVolData){}
	DECLARE_REGISTRY_RESOURCEID(IDR_ARPROPVOLDATA)
	DECLARE_NOT_AGGREGATABLE(CARPropVolData)
	DECLARE_PROTECT_FINAL_CONSTRUCT()
	BEGIN_COM_MAP(CARPropVolData)
		COM_INTERFACE_ENTRY(IARPropVolData)
		COM_INTERFACE_ENTRY(IDispatch)
		COM_INTERFACE_ENTRY(ISupportErrorInfo)
	END_COM_MAP()
	HRESULT								FinalConstruct(void);

	STDMETHOD(InterfaceSupportsErrorInfo)(REFIID riid);
	STDMETHOD(get_Value)(/*[out, retval]*/ VARIANT *pVal);
	STDMETHOD(put_Value)(/*[in]*/ VARIANT newVal);
};

#endif //__ARPROPVOLDATA_H_
