//	handle_interpolator.h : Declaration of the CInterpolator
//
/////////////////////////////////////////////////////////////////////////////

#ifndef __INTERPOLATOR_H_
#define __INTERPOLATOR_H_

#include "resource.h"
#include "MlEqInterpolator.h"

class ATL_NO_VTABLE CInterpolator : 
	public CComObjectRootEx<CComSingleThreadModel>,
	public CComCoClass<CInterpolator, &CLSID_Interpolator>,
	public ISupportErrorInfo,
	public IDispatchImpl<IInterpolator, &IID_IInterpolator, &LIBID_Sirius>
{
public:		
	CInterpolator(){}
	HRESULT								FinalConstruct(void);
	DECLARE_REGISTRY_RESOURCEID(IDR_INTERPOLATOR)
	DECLARE_NOT_AGGREGATABLE(CInterpolator)
	DECLARE_PROTECT_FINAL_CONSTRUCT()
	BEGIN_COM_MAP(CInterpolator)
		COM_INTERFACE_ENTRY(IInterpolator)
		COM_INTERFACE_ENTRY(IDispatch)
		COM_INTERFACE_ENTRY(ISupportErrorInfo)
	END_COM_MAP()
	associate_analytic_object(CInterpolator, MlEqInterpolator);

	STDMETHOD(InterfaceSupportsErrorInfo)(REFIID riid);
	STDMETHOD(get_Value)(/*[out, retval]*/ VARIANT *pVal);
	STDMETHOD(put_Value)(/*[in]*/ VARIANT newVal);
	STDMETHOD(get_InterpolatorType)(/*[out, retval]*/ InterpolatorTypeEnum *pVal);
	STDMETHOD(GetValue)(/*[in]*/ double xValue, long whichDataOpt,/*[out, retval]*/ double* pVal);
	
public:	
	std::vector<CParameterMap>			m_vpm;
};

#endif
