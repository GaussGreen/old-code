//	jumpwingvoldata.h : Declaration of the CJumpWingVolData
//
/////////////////////////////////////////////////////////////////////////////

#ifndef __JUMPWINGVOLDATA_H_
#define __JUMPWINGVOLDATA_H_

#include "resource.h"
#include "MlEqObjects.h"

class ATL_NO_VTABLE CJumpWingVolData : 
	public CComObjectRootEx<CComSingleThreadModel>,
	public CComCoClass<CJumpWingVolData, &CLSID_JumpWingVolData>,
	public ISupportErrorInfo,
	public IDispatchImpl<IJumpWingVolData, &IID_IJumpWingVolData, &LIBID_Sirius>
{
public:	
	CJumpWingVolData() : m_h(new MlEqJumpWingVolData) {}
	DECLARE_REGISTRY_RESOURCEID(IDR_JUMPWINGVOLDATA)
	DECLARE_NOT_AGGREGATABLE(CJumpWingVolData)
	DECLARE_PROTECT_FINAL_CONSTRUCT()
	BEGIN_COM_MAP(CJumpWingVolData)
		COM_INTERFACE_ENTRY(IJumpWingVolData)
		COM_INTERFACE_ENTRY(IDispatch)
		COM_INTERFACE_ENTRY(ISupportErrorInfo)
	END_COM_MAP()
	associate_analytic_object(CJumpWingVolData, MlEqJumpWingVolData);

	STDMETHOD(InterfaceSupportsErrorInfo)(REFIID riid);
	STDMETHOD(get_Value)(/*[out, retval]*/ VARIANT *pVal);
	STDMETHOD(put_Value)(/*[in]*/ VARIANT newVal);
	STDMETHOD(get_JumpWingParameters)(/*[out, retval]*/ VARIANT* pVal);

protected:
	std::vector<CAdapt<CComQIPtr<IStrikes> > >	m_aspStrike;	
	CParameterMap								m_pmJumpWingCoefficients;
};

#endif
