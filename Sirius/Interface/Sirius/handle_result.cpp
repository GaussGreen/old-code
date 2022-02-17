//	handle_result.cpp : Implementation of CResult
//
/////////////////////////////////////////////////////////////////////////////

#include "stdafx.h"
#include "handle_result.h"
#include "siriusapplication.h"


const CResult& CResult::operator+=(const CResult& r)
{
	ATLASSERT(!m_mapStringToVariant.size());	// ToDo - support non-numeric cases somehow
	std::map<std::string, double>::const_iterator it_other = r.m_mapStringToDouble.begin();
	std::map<std::string, double>::iterator it_this = m_mapStringToDouble.begin();

	while (it_other != r.m_mapStringToDouble.end() || it_this != m_mapStringToDouble.end()){		
		if (it_other == r.m_mapStringToDouble.end() || (it_this != m_mapStringToDouble.end() && it_this->first < it_other->first)){
			// we don't need to add anything to it_this
			it_this++;
		} else if (it_this == m_mapStringToDouble.end() || it_this->first > it_other->first){
			// insert it_other into m_mapStringToDouble
			{
				std::map<std::string, std::string>::const_iterator it = r.m_mapLowerCaseToProperCase.find(it_other->first);
				m_mapLowerCaseToProperCase[it_other->first] = it->second;
			}
			{
				std::map<std::string, double>::const_iterator it = r.m_mapStringToDouble.find(it_other->first);
				m_mapStringToDouble[it_other->first] = it->second;
			}
			it_other++;
		} else /*it_this->first == it_other->first*/{		
			// add it_other to it_this
			it_this->second += it_other->second;
			it_this++;
			it_other++;
		}		
	}		
	return *this;
}

const CResult& CResult::operator*=(double f)
{
	ATLASSERT(!m_mapStringToVariant.size());	// ToDo - support non-numeric cases somehow
	for (std::map<std::string, double>::iterator it = m_mapStringToDouble.begin(); it != m_mapStringToDouble.end(); ++it){
		it->second *= f;
	}
	return *this;
}

STDMETHODIMP CResult::Add(IResult* ResultObject)
{
	begin_function
	CResult* pResult = dynamic_cast<CResult*>(ResultObject);
	operator+=(*pResult);
	end_function
}

STDMETHODIMP CResult::AddValue(BSTR Name, VARIANT Value)
{		
	CComVariant							vDouble;						// double equivalent of Value
	estring								szName(Name);
	estring								szName_lc(Name);				// lower case form of szName

	szName_lc.StripWhiteSpace();
	szName_lc.lc();
	if (!szName.size()) return CParameterMap::ReturnErrorR(IDS_INVALID_NAME, IID_IResult);
	m_mapLowerCaseToProperCase[szName_lc] = szName;
		
	if (!vDouble.ChangeType(VT_R8, &Value)){
		// use the string to double map
		m_mapStringToDouble[szName_lc] = vDouble.dblVal;		
	} else {
		// use the string to variant map
		m_mapStringToVariant[szName_lc] = Value;
	}
	return S_OK;
}

STDMETHODIMP CResult::CopyToClipboard(void)
{
	begin_function
	HRESULT								hr;
	CComVariant							v;
	CParameterMap						pm;
	
	if (hr = get_Value(&v)) return hr;
	if (hr = pm.SetValue(v)) return hr;
	return pm.CopyToClipboard();	
	end_function
}

HRESULT CResult::FinalConstruct()
{	
	m_bCalculateAll = false;
	return S_OK;
}

HRESULT CResult::GetDouble(const std::string& szName, const std::string& szName_lc, double* pVal) const
{	
	std::map<std::string, double>::const_iterator it = m_mapStringToDouble.find(szName_lc);
	if (it == m_mapStringToDouble.end()){
		return CParameterMap::ReturnErrorRS(IDS_DATA_ELEMENT_NOT_FOUND, szName, IID_IResult);
	}
	*pVal = it->second;
	return S_OK;
}

// returns a list of parameters in the result object
HRESULT CResult::GetParameters(VARIANT* pVal)
{
	CParameterMap						pm(m_mapStringToDouble.size() + m_mapStringToVariant.size(), 1);
	long								nRow(0);
	HRESULT								hr = S_OK;
		
	// string to double map	
	for (std::map<std::string, double>::const_iterator it = m_mapStringToDouble.begin(); it != m_mapStringToDouble.end(); ++it){								
		pm.SetValue(nRow, 0, m_mapLowerCaseToProperCase[it->first]);	
		nRow++;
	}

	// string to variant map	
	for (std::map<std::string, CComVariant>::const_iterator it = m_mapStringToVariant.begin(); it != m_mapStringToVariant.end(); ++it){
		pm.SetValue(nRow, 0, m_mapLowerCaseToProperCase[it->first]);			
		nRow++;			
	}	
	return pm.GetValue(pVal);
}

STDMETHODIMP CResult::GetPrice(double *pVal)
{	
	// it's OK to hardcode the string here since the GetPrice method name is hardcoded in the same sense	
	return GetDouble("price", "price", pVal);
}

STDMETHODIMP CResult::GetTheta(double *pVal)
{
	// it's OK to hardcode the string here since the GetTheta method name is hardcoded in the same sense
	return GetDouble("theta", "theta", pVal);
}

//	returns the value stored in a given property name
STDMETHODIMP CResult::GetValue(BSTR Name, VARIANT *pVal)
{
	estring								szName(Name);
	estring								szName_lc(Name);
	double								fValue;

	szName_lc.StripWhiteSpace();
	szName_lc.lc();

	// try the double map first
	if (!GetDouble(szName, szName_lc, &fValue)){
		// found it
		return CComVariant(fValue).Detach(pVal);
	}

	// try the variant map
	std::map<std::string, CComVariant>::const_iterator it = m_mapStringToVariant.find(szName_lc);
	if (it == m_mapStringToVariant.end()){
		return CParameterMap::ReturnErrorRS(IDS_DATA_ELEMENT_NOT_FOUND, szName, IID_IResult);
	}
	return CComVariant(it->second).Detach(pVal);	
}

STDMETHODIMP CResult::get_Value(VARIANT* pVal)
{
	begin_function
	
	CParameterMap						pm(m_mapStringToDouble.size() + m_mapStringToVariant.size(), 2);
	long								nRow(0);	
			
	// special case where the handle only contains "Price"
	if (m_mapStringToDouble.size() == 1 && m_mapStringToVariant.size() == 0){
		std::map<std::string, double>::const_iterator it = m_mapStringToDouble.find("price");
		if (it != m_mapStringToDouble.end()){
			return CComVariant(it->second).Detach(pVal);		
		}
	}
	
	// string to double map	
	for (std::map<std::string, double>::const_iterator it = m_mapStringToDouble.begin(); it != m_mapStringToDouble.end(); ++it){								
		pm.SetValue(nRow, 0, m_mapLowerCaseToProperCase[it->first]);
		pm.SetValue(nRow, 1, it->second);
		nRow++;
	}
	
	// string to variant map	
	for (std::map<std::string, CComVariant>::const_iterator it = m_mapStringToVariant.begin(); it != m_mapStringToVariant.end(); ++it){
		pm.SetValue(nRow, 0, m_mapLowerCaseToProperCase[it->first]);
		pm.SetValue(nRow, 1, it->second);			
		nRow++;			
	}
	
	if (pm.IsBlank()) throw CStringResource(IDS_NO_DATA_FOUND);
	return pm.GetValue(pVal);
	end_function
}

STDMETHODIMP CResult::HasParameter(BSTR Name, VARIANT_BOOL* pVal)
{
	begin_function
	estring szName(Name);
	szName.lc();		
	*pVal = m_mapLowerCaseToProperCase.find(szName) != m_mapLowerCaseToProperCase.end() ? VARIANT_TRUE : VARIANT_FALSE;
	end_function
}

STDMETHODIMP CResult::InterfaceSupportsErrorInfo(REFIID riid)
{
	static const IID* arr[] = {&IID_IResult};
	for (int i=0; i < sizeof(arr) / sizeof(arr[0]); i++){
		if (InlineIsEqualGUID(*arr[i],riid)) return S_OK;
	}
	return S_FALSE;
}

STDMETHODIMP CResult::Multiply(double Amount)
{
	operator*=(Amount);
	return S_OK;
}

STDMETHODIMP CResult::Calculate(BSTR Name, VARIANT_BOOL* pVal)
{
	estring								szName_lc(Name);
	szName_lc.StripWhiteSpace();
	szName_lc.lc();

	if (!szName_lc.size()) return CParameterMap::ReturnErrorR(IDS_INVALID_NAME, IID_IResult);
	if (m_bCalculateAll) {
		*pVal = true;
		return S_OK;
	}
	if (szName_lc == "price" && !m_mapCalculate.size()){
		*pVal = true;
		return S_OK;
	}
	std::map<std::string, bool>::const_iterator it = m_mapCalculate.find(szName_lc);
	*pVal = (it != m_mapCalculate.end());
	return S_OK;
}

STDMETHODIMP CResult::put_Value(VARIANT newVal)
{
	CParameterMap						pm;
	estring								szName;
	estring								szName_lc;
	HRESULT								hr;
	CParameterMap						pmElement;	
	
	if (hr = pm.SetValue(newVal)) return hr;
	if (pm.IsDouble()){
		// special case where the handle will only contain "Price" - see get_Value
		pm.InsertColumn(0);
		pm.SetValue(0, 0, std::string("Price"));
		return put_Value(pm.GetValue());
	} else {										
		if (!pm.GetRows()) return CParameterMap::ReturnErrorR(IDS_ROWS_INVALID, IID_IResult);
		if (pm.GetCols() != 2) return CParameterMap::ReturnErrorR(IDS_COLUMNS_INVALID, IID_IResult);						
		m_bCalculateAll = false;
		m_mapLowerCaseToProperCase.clear(); // ToDo - don't reassign this until the end?
		m_mapStringToDouble.clear();
		m_mapStringToVariant.clear();
		m_mapCalculate.clear();
		for (long nRow = 0; nRow < pm.GetRows(); nRow++){
			if (pm.GetValue(nRow, 0, &szName)) return CParameterMap::ReturnErrorR(IDS_INVALID_NAME, IID_IResult);
			szName_lc.assign(szName);
			szName_lc.StripWhiteSpace();
			szName_lc.lc();
			if (!szName_lc.size()) return CParameterMap::ReturnErrorR(IDS_INVALID_NAME, IID_IResult);
			if (m_mapLowerCaseToProperCase.find(szName_lc) != m_mapLowerCaseToProperCase.end()) return CParameterMap::ReturnErrorRS(IDS_NAME_ALREADY_IN_USE, szName_lc, IID_IResult);
			if (pm.GetValue(nRow, 1, &pmElement)) return CParameterMap::ReturnErrorR(IDS_INVALID_PARAMETER_VALUE, IID_IResult);
			if (pmElement.IsScalar() && pmElement.IsDouble()){
				// add the value to the string to double map
				double f;
				if (pmElement.GetValue(&f)) ATLASSERT(false);
				m_mapStringToDouble[szName_lc] = f;
			} else {
				// add the value to the string to variant map
				CComVariant v;
				if (pmElement.GetValue(&v)) ATLASSERT(false);
				m_mapStringToVariant[szName_lc] = v;
			}
			m_mapLowerCaseToProperCase[szName_lc] = szName;
		}
	}
	return S_OK;
}

STDMETHODIMP CResult::SetCalculate(BSTR Calculate)
{
	estring										szCalculate(Calculate);
	std::vector<std::string>					vectorCalculate;
	std::vector<std::string>::const_iterator	it;

	m_mapCalculate.clear();
	m_bCalculateAll = false;
	szCalculate.StripWhiteSpace();
	szCalculate.lc();
	szCalculate.Split(",", &vectorCalculate);
	
	// check to see if default has been inputted (blank input implies just a price calculation)
	if (!vectorCalculate.size()){
		// this is in an instruction to calculate just the price (i.e. the default behaviour)
		m_bCalculateAll = false;
		return S_OK;
	}
	
	// Set up m_mapCalculate. We set m_bCalculateAll if the "all" element is found
	for (it = vectorCalculate.begin(); it < vectorCalculate.end(); it++){
		if (*it == "all"){
			m_bCalculateAll = true;
			m_mapCalculate.clear();
			return S_OK;
		}
		m_mapCalculate[*it] = true;
	}
	return S_OK;
}