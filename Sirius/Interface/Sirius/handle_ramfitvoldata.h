//	handle_ramfitvoldata.h : Declaration of the CRamFitVolData
//
/////////////////////////////////////////////////////////////////////////////

#ifndef __RAMFITVOLDATA_H_
#define __RAMFITVOLDATA_H_

#include "resource.h"
#include "MlEqObjects.h"

class ATL_NO_VTABLE CRamFitVolData : 
	public CComObjectRootEx<CComSingleThreadModel>,
	public CComCoClass<CRamFitVolData, &CLSID_RamFitVolData>,
	public ISupportErrorInfo,
	public IDispatchImpl<IRamFitVolData, &IID_IRamFitVolData, &LIBID_Sirius>	
{
public:
	RamFitHandle						m_h;
	CRamFitVolData() : m_h(new RamFit){}
	DECLARE_REGISTRY_RESOURCEID(IDR_RAMFITVOLDATA)
	DECLARE_NOT_AGGREGATABLE(CRamFitVolData)
	DECLARE_PROTECT_FINAL_CONSTRUCT()
	BEGIN_COM_MAP(CRamFitVolData)
		COM_INTERFACE_ENTRY(IRamFitVolData)
		COM_INTERFACE_ENTRY(IDispatch)
		COM_INTERFACE_ENTRY(ISupportErrorInfo)
	END_COM_MAP()

	STDMETHOD(InterfaceSupportsErrorInfo)(REFIID riid);
	STDMETHOD(get_Value)(/*[out, retval]*/ VARIANT *pVal);
	STDMETHOD(put_Value)(/*[in]*/ VARIANT newVal);

protected:
	CComPtr<IVolatilityStructure>				m_spVolatilityStructure;
	std::vector<CAdapt<CComQIPtr<IStrikes> > >	m_aspFittingStrikes;
	std::vector<CAdapt<CComQIPtr<IStrikes> > >	m_aspUpStrikes;
	std::vector<CAdapt<CComQIPtr<IStrikes> > >	m_aspDownStrikes;
	std::vector<long>							m_nFittingDates;
};

#endif
