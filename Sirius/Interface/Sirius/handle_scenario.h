//	handle_scenario.h : Declaration of the CScenario
//
/////////////////////////////////////////////////////////////////////////////

#ifndef __SCENARIO_H_
#define __SCENARIO_H_

#include "resource.h"
#include "mleqscenario.h"				// for m_h

class ATL_NO_VTABLE CScenario : 
	public CComObjectRootEx<CComSingleThreadModel>,
	public CComCoClass<CScenario, &CLSID_Scenario>,
	public ISupportErrorInfo,
	public IDispatchImpl<IScenario, &IID_IScenario, &LIBID_Sirius>
{
public:
	CScenario() : m_h(new MlEqScenario) {}
	DECLARE_REGISTRY_RESOURCEID(IDR_SCENARIO)
	DECLARE_NOT_AGGREGATABLE(CScenario)
	DECLARE_PROTECT_FINAL_CONSTRUCT()
	BEGIN_COM_MAP(CScenario)
		COM_INTERFACE_ENTRY(IScenario)
		COM_INTERFACE_ENTRY(IDispatch)
		COM_INTERFACE_ENTRY(ISupportErrorInfo)
	END_COM_MAP()
	associate_analytic_object(CScenario, MlEqScenario)

	STDMETHOD(InterfaceSupportsErrorInfo)(REFIID riid);
	STDMETHOD(get_Value)(/*[out, retval]*/ VARIANT *pVal);
	STDMETHOD(put_Value)(/*[in]*/ VARIANT newVal);
};

#endif
