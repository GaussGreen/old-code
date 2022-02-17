//	handle_jumpwingfitvoldata.h : Declaration of the CJumpWingFitVolData
//
/////////////////////////////////////////////////////////////////////////////

#ifndef __JUMPWINGFITVOLDATA_H_
#define __JUMPWINGFITVOLDATA_H_

#include "resource.h"
#include "MlEqObjects.h"

class ATL_NO_VTABLE CJumpWingFitVolData : 
	public CComObjectRootEx<CComSingleThreadModel>,
	public CComCoClass<CJumpWingFitVolData, &CLSID_JumpWingFitVolData>,
	public ISupportErrorInfo,
	public IDispatchImpl<IJumpWingFitVolData, &IID_IJumpWingFitVolData, &LIBID_Sirius>
{
public:	
	CJumpWingFitVolData() : m_h(new MlEqJumpWingFitVolData) {}
	DECLARE_REGISTRY_RESOURCEID(IDR_JUMPWINGFITVOLDATA)
	DECLARE_NOT_AGGREGATABLE(CJumpWingFitVolData)
	DECLARE_PROTECT_FINAL_CONSTRUCT()
	BEGIN_COM_MAP(CJumpWingFitVolData)
		COM_INTERFACE_ENTRY(IJumpWingFitVolData)
		COM_INTERFACE_ENTRY(IDispatch)
		COM_INTERFACE_ENTRY(ISupportErrorInfo)
	END_COM_MAP()
	associate_analytic_object(CJumpWingFitVolData, MlEqJumpWingFitVolData);
	HRESULT								FinalConstruct(void);

	STDMETHOD(InterfaceSupportsErrorInfo)(REFIID riid);
	STDMETHOD(get_Value)(/*[out, retval]*/ VARIANT *pVal);
	STDMETHOD(put_Value)(/*[in]*/ VARIANT newVal);
	STDMETHOD(get_JumpWingParameters)(/*[out, retval]*/ VARIANT* pVal);

protected:
	// member variables that Jimmy will get rid of go here
	std::vector<CParameterMap>			m_vpm;
};

#endif
