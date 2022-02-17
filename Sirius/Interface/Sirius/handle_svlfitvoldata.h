//	handle_svlfitvoldata.h : Declaration of the CSVLFitVolData
//
/////////////////////////////////////////////////////////////////////////////

#ifndef __SVLFITVOLDATA_H_
#define __SVLFITVOLDATA_H_

#include "resource.h"
#include "MlEqObjects.h"


class ATL_NO_VTABLE CSVLFitVolData : 
	public CComObjectRootEx<CComSingleThreadModel>,
	public CComCoClass<CSVLFitVolData, &CLSID_SVLFitVolData>,
	public ISupportErrorInfo,
	public IDispatchImpl<ISVLFitVolData, &IID_ISVLFitVolData, &LIBID_Sirius>
{
public:
	SvlFitHandle						m_h;
	CSVLFitVolData() : m_h(new SvlFit){}
	DECLARE_REGISTRY_RESOURCEID(IDR_SVLFITVOLDATA)
	DECLARE_NOT_AGGREGATABLE(CSVLFitVolData)
	DECLARE_PROTECT_FINAL_CONSTRUCT()
	BEGIN_COM_MAP(CSVLFitVolData)
		COM_INTERFACE_ENTRY(ISVLFitVolData)
		COM_INTERFACE_ENTRY(IDispatch)
		COM_INTERFACE_ENTRY(ISupportErrorInfo)
	END_COM_MAP()
	HRESULT								FinalConstruct(void);

	STDMETHOD(InterfaceSupportsErrorInfo)(REFIID riid);
	STDMETHOD(get_Value)(/*[out, retval]*/ VARIANT *pVal);
	STDMETHOD(put_Value)(/*[in]*/ VARIANT newVal);

protected:
	std::vector<CAdapt<CComQIPtr<IStrikes> > >	m_aspStrike;
	CMatrix										m_InitialGuess;
};

#endif
