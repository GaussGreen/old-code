//	handle_hermitevoldata.h : Declaration of the CHermiteVolData
//
/////////////////////////////////////////////////////////////////////////////

#ifndef __HERMITEVOLDATA_H_
#define __HERMITEVOLDATA_H_

#include "resource.h"
#include "MlEqObjects.h"

class ATL_NO_VTABLE CHermiteVolData : 
	public CComObjectRootEx<CComSingleThreadModel>,
	public CComCoClass<CHermiteVolData, &CLSID_HermiteVolData>,
	public ISupportErrorInfo,
	public IDispatchImpl<IHermiteVolData, &IID_IHermiteVolData, &LIBID_Sirius>
{
public:
	CHermiteVolData() : m_h(new MLHermiteVolData) {}
	DECLARE_REGISTRY_RESOURCEID(IDR_HERMITEVOLDATA)
	DECLARE_NOT_AGGREGATABLE(CHermiteVolData)
	DECLARE_PROTECT_FINAL_CONSTRUCT()
	BEGIN_COM_MAP(CHermiteVolData)
		COM_INTERFACE_ENTRY(IHermiteVolData)
		COM_INTERFACE_ENTRY(IDispatch)
		COM_INTERFACE_ENTRY(ISupportErrorInfo)
	END_COM_MAP()
	associate_analytic_object(CHermiteVolData, MLHermiteVolData);
	
	STDMETHOD(InterfaceSupportsErrorInfo)(REFIID riid);
	STDMETHOD(get_Value)(/*[out, retval]*/ VARIANT* pVal);
	STDMETHOD(put_Value)(/*[in]*/ VARIANT newVal);

protected:
	CComPtr<IDate>						m_spStartDate;
	std::vector<CParameterMap>			m_vpm;
};

#endif
