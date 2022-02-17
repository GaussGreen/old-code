//	handle_multipliervoldata.h : Declaration of the CMultiplierVolData
//
/////////////////////////////////////////////////////////////////////////////

#ifndef __MULTIPLIERVOLDATA_H_
#define __MULTIPLIERVOLDATA_H_

#include "resource.h"
#include "MlEqObjects.h"

class ATL_NO_VTABLE CMultiplierVolData : 
	public CComObjectRootEx<CComSingleThreadModel>,
	public CComCoClass<CMultiplierVolData, &CLSID_MultiplierVolData>,
	public ISupportErrorInfo,
	public IDispatchImpl<IMultiplierVolData, &IID_IMultiplierVolData, &LIBID_Sirius>	
{
public:
	MlEqVolMultiplierDataHandle			m_h;	
	CMultiplierVolData() : m_h(new MlEqVolMultiplierData){}
	DECLARE_REGISTRY_RESOURCEID(IDR_MULTIPLIERVOLDATA)
	DECLARE_NOT_AGGREGATABLE(CMultiplierVolData)
	DECLARE_PROTECT_FINAL_CONSTRUCT()
	BEGIN_COM_MAP(CMultiplierVolData)
		COM_INTERFACE_ENTRY(IMultiplierVolData)
		COM_INTERFACE_ENTRY(IDispatch)
		COM_INTERFACE_ENTRY(ISupportErrorInfo)
	END_COM_MAP()
	STDMETHOD(InterfaceSupportsErrorInfo)(REFIID riid);
	STDMETHOD(get_Value)(/*[out, retval]*/ VARIANT *pVal);
	STDMETHOD(put_Value)(/*[in]*/ VARIANT newVal);

protected:
	std::vector<CAdapt<CComQIPtr<IStrikes> > >			m_aspStrike;
	std::vector<CAdapt<CComQIPtr<IInterpolator> > >		m_aspInterpolator;
};

#endif
