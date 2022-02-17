//	handle_result.h : Declaration of the CResult
//
/////////////////////////////////////////////////////////////////////////////

#ifndef __RESULT_H_
#define __RESULT_H_

#include "resource.h"       // main symbols

class ATL_NO_VTABLE CResult : 
	public CComObjectRootEx<CComSingleThreadModel>,
	public CComCoClass<CResult, &CLSID_Result>,
	public ISupportErrorInfo,
	public IDispatchImpl<IResult, &IID_IResult, &LIBID_Sirius>
{
public:				
	CResult(){}
	DECLARE_REGISTRY_RESOURCEID(IDR_RESULT)
	DECLARE_NOT_AGGREGATABLE(CResult)
	DECLARE_PROTECT_FINAL_CONSTRUCT()
	BEGIN_COM_MAP(CResult)
		COM_INTERFACE_ENTRY(IResult)
		COM_INTERFACE_ENTRY(IDispatch)
		COM_INTERFACE_ENTRY(ISupportErrorInfo)
	END_COM_MAP()
	HRESULT								FinalConstruct();
	const CResult&						operator*=(double f);
	const CResult&						operator+=(const CResult& r);
	BOOL								operator==(const CResult&)	{ATLASSERT(false); return FALSE;/*do not use*/}

	STDMETHOD(InterfaceSupportsErrorInfo)(REFIID riid);
	STDMETHOD(get_Value)(/*[out, retval]*/ VARIANT *pVal);
	STDMETHOD(put_Value)(/*[in]*/ VARIANT newVal);
	STDMETHOD(GetTheta)(/*[out, retval]*/ double* pVal);
	STDMETHOD(GetPrice)(/*[out, retval]*/ double* pVal);
	STDMETHOD(AddValue)(/*[in]*/ BSTR Name, /*[in]*/ VARIANT Value);
	STDMETHOD(Add)(/*[in]*/ IResult* ResultObject);
	STDMETHOD(Multiply)(/*[in]*/ double Amount);
	STDMETHOD(GetValue)(/*[in]*/ BSTR Name, /*[out, retval]*/ VARIANT* pVal);
	STDMETHOD(SetCalculate)(/*[in]*/ BSTR Calculate);
	STDMETHOD(Calculate)(/*[in]*/ BSTR Name, /*[out, retval]*/ VARIANT_BOOL* pVal);
	STDMETHOD(CopyToClipboard)(void);
	STDMETHOD(HasParameter)(/*[in]*/ BSTR Name, /*[out, retval]*/ VARIANT_BOOL* pVal);

	HRESULT										GetParameters(VARIANT* pVal);
	const std::map<std::string, double>&		GetStringToDoubleMap(void) const {return m_mapStringToDouble;}
	const std::map<std::string, CComVariant>&	GetStringToVariantMap(void) const {return m_mapStringToVariant;}
	const std::map<std::string, std::string>&	GetLowerCaseToProperCaseMap(void) const {return m_mapLowerCaseToProperCase;}

protected:	
	std::map<std::string, double>		m_mapStringToDouble;			// handles all the numeric maps
	std::map<std::string, std::string>	m_mapLowerCaseToProperCase;		// transforms the internal lower case representation of the input values to the originally specified string
	std::map<std::string, CComVariant>	m_mapStringToVariant;			// handles all the other result handle value cases
	std::map<std::string, bool>			m_mapCalculate;					// contains a map of all values that should be calculated by a given product type
	boolean								m_bCalculateAll;				// TRUE if we calculate all values

	HRESULT								GetDouble(const std::string& szName, const std::string& szName_lc, double* pVal) const;
};

#endif
