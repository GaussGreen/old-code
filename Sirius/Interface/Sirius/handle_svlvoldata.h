//	handle_svlvoladata.h : Declaration of the CSVLVolData
//
/////////////////////////////////////////////////////////////////////////////

#ifndef __SVLVOLDATA_H_
#define __SVLVOLDATA_H_

#include "resource.h"
#include "MlEqObjects.h"

class ATL_NO_VTABLE CSVLVolData : 
	public CComObjectRootEx<CComSingleThreadModel>,
	public CComCoClass<CSVLVolData, &CLSID_SVLVolData>,
	public ISupportErrorInfo,
	public IDispatchImpl<ISVLVolData, &IID_ISVLVolData, &LIBID_Sirius>	
{
public:
	MlSvlVolDataHandle					m_h;
	CSVLVolData() : m_h(new MlSvlVolData){}
	DECLARE_REGISTRY_RESOURCEID(IDR_SVLVOLDATA)
	DECLARE_NOT_AGGREGATABLE(CSVLVolData)
	DECLARE_PROTECT_FINAL_CONSTRUCT()
	BEGIN_COM_MAP(CSVLVolData)
		COM_INTERFACE_ENTRY(ISVLVolData)
		COM_INTERFACE_ENTRY(IDispatch)
		COM_INTERFACE_ENTRY(ISupportErrorInfo)
	END_COM_MAP()

	STDMETHOD(InterfaceSupportsErrorInfo)(REFIID riid);
	STDMETHOD(get_Value)(/*[out, retval]*/ VARIANT *pVal);
	STDMETHOD(put_Value)(/*[in]*/ VARIANT newVal);
	STDMETHOD(get_BasicSvlParameters)(/*[out, retval]*/ VARIANT* pVal);
	STDMETHOD(put_BasicSvlParameters)(/*[in]*/ VARIANT newVal);
	STDMETHOD(get_BasicValue)(/*[out, retval]*/ VARIANT *pVal);
	STDMETHOD(put_BasicValue)(/*[in]*/ VARIANT newVal);
	STDMETHOD(get_SvlParameters)(/*[out, retval]*/ VARIANT* pVal);
	STDMETHOD(put_SvlParameters)(/*[in]*/ VARIANT newVal);

protected:
	static bool							GetFromNew(const std::vector<SvlCoeffs>& svl, CParameterMap* ppm);	// Return type is defined as true for Svl function, false for Svi. Odd but useful.
	static bool							GetFromOld(const std::vector<SvlCoeffs>& svl, CParameterMap* ppm);	// Return type is defined as true for Svl function, false for Svi. Odd but useful.
	HRESULT								GetValue(VARIANT* pVal, bool (*Get)(const std::vector<SvlCoeffs>& , CParameterMap*));
	static bool							PutFromNew(const CParameterMap& pm, std::vector<SvlCoeffs>* psvl);  // Return type is defined as true for Svl function, false for Svi. Odd but useful.
	static bool							PutFromOld(const CParameterMap& pm, std::vector<SvlCoeffs>* psvl);	// Return type is defined as true for Svl function, false for Svi. Odd but useful.
	void								PutValue(const std::vector<CComVariant>& vv, bool (*Put)(const CParameterMap&, std::vector<SvlCoeffs>*));
	
	std::vector<CAdapt<CComQIPtr<IStrikes> > >	m_aspStrike;
};

#endif
