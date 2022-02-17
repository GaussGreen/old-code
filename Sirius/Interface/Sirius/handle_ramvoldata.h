//	handle_ramvoldata.h : Declaration of the CRamVolData
//
/////////////////////////////////////////////////////////////////////////////

#ifndef __RAMVOLDATA_H_
#define __RAMVOLDATA_H_

#include "resource.h"
#include "MlEqObjects.h"

class ATL_NO_VTABLE CRamVolData : 
	public CComObjectRootEx<CComSingleThreadModel>,
	public CComCoClass<CRamVolData, &CLSID_RamVolData>,
	public ISupportErrorInfo,
	public IDispatchImpl<IRamVolData, &IID_IRamVolData, &LIBID_Sirius>	
{
public:
	MLRamVolDataHandle					m_h;
	CRamVolData() : m_h(new MLRamVolData){}
	DECLARE_REGISTRY_RESOURCEID(IDR_RAMVOLDATA)	
	DECLARE_NOT_AGGREGATABLE(CRamVolData)
	DECLARE_PROTECT_FINAL_CONSTRUCT()
	BEGIN_COM_MAP(CRamVolData)
		COM_INTERFACE_ENTRY(IRamVolData)
		COM_INTERFACE_ENTRY(IDispatch)
		COM_INTERFACE_ENTRY(ISupportErrorInfo)
	END_COM_MAP()
		
	STDMETHOD(InterfaceSupportsErrorInfo)(REFIID riid);
	STDMETHOD(get_Value)(/*[out, retval]*/ VARIANT *pVal);
	STDMETHOD(put_Value)(/*[in]*/ VARIANT newVal);
	STDMETHOD(get_BasicParameters)(/*[out, retval]*/ VARIANT* pVal);
	STDMETHOD(get_Parameters)(/*[out, retval]*/ VARIANT* pVal);

protected:
	std::vector<CAdapt<CComQIPtr<IStrikes> > >			m_aspStrike;
	CVector												m_BaseVolArr;
	CVector												m_SlopeArr;
	CVector												m_CurvatureArr;
	std::vector<CAdapt<CComQIPtr<IStrikes> > >			m_aspDownCutoffArr;
	std::vector<CAdapt<CComQIPtr<IStrikes> > >			m_aspUpCutoffArr;
};

#endif
