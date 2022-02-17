//	parameter.cpp : Implementation of CParameter
//
/////////////////////////////////////////////////////////////////////////////

#include "stdafx.h"
#include "parameter.h"

STDMETHODIMP CParameter::AddToEnd(IParameter* Val)
{
	CParameterMap*	pm = dynamic_cast<CParameterMap*>(Val);
	if (!pm) return E_POINTER;
	CParameterMap::AddToEnd(*pm);
	return S_OK;
}
STDMETHODIMP CParameter::AddToRHS(IParameter* Val)
{
	CParameterMap*	pm = dynamic_cast<CParameterMap*>(Val);
	if (!pm) return E_POINTER;	
	CParameterMap::AddToRHS(*pm);
	return S_OK;
}

STDMETHODIMP CParameter::get_AutoGrow(VARIANT_BOOL* pVal)
{
	*pVal = GetAutoGrow() ? VARIANT_TRUE : VARIANT_FALSE;	
	return S_OK;
}

STDMETHODIMP CParameter::CopyToClipboard()
{
	return CParameterMap::CopyToClipboard();
}

STDMETHODIMP CParameter::GetColumn(long Column, IParameter** pVal)
{
	CComPtr<IParameter>			spParameter;
	CParameter*					ppm;
	HRESULT						hr;

	if (hr = spParameter.CoCreateInstance(CLSID_Parameter)) return hr;
	ppm = dynamic_cast<CParameter*>(spParameter.p);
	if (hr = CParameterMap::GetColumn(Column - 1, ppm)) return hr;
	return spParameter.CopyTo(pVal);
}

STDMETHODIMP CParameter::get_Columns(long *pVal)
{
	*pVal = CParameterMap::GetCols();
	return S_OK;
}

STDMETHODIMP CParameter::get_Element(long Row, long Column, VARIANT *pVal)
{
	begin_function
	if (Row < 1 || Row > m_nRows || Column < 1 || Column > m_nCols) throw CStringResource(IDS_ARRAY_INDEX);
	return CParameterMap::GetValue(Row - 1, Column - 1, pVal);	
	end_function
}

STDMETHODIMP CParameter::get_Name(BSTR *pVal)
{
	return m_szName.GetBSTR(pVal);
}

STDMETHODIMP CParameter::GetRow(long Row, IParameter** pVal)
{
	CComPtr<IParameter>			spParameter;
	CParameter*					ppm;
	HRESULT						hr;

	if (hr = spParameter.CoCreateInstance(CLSID_Parameter)) return hr;
	ppm = dynamic_cast<CParameter*>(spParameter.p);
	if (hr = CParameterMap::GetRow(Row - 1, ppm)) return hr;
	return spParameter.CopyTo(pVal);
}

STDMETHODIMP CParameter::get_Rows(long *pVal)
{	
	*pVal = CParameterMap::GetRows();
	return S_OK;
}

STDMETHODIMP CParameter::InsertColumn(long At)
{
	return CParameterMap::InsertColumn(At - 1);
}

STDMETHODIMP CParameter::InsertRow(long At)
{
	return CParameterMap::InsertRow(At - 1);
}

STDMETHODIMP CParameter::InterfaceSupportsErrorInfo(REFIID riid)
{
	static const IID* arr[] = { &IID_IParameter };
	for (int i=0; i < sizeof(arr) / sizeof(arr[0]); i++){
		if (InlineIsEqualGUID(*arr[i],riid)) return S_OK;
	}
	return S_FALSE;
}

STDMETHODIMP CParameter::IsBlank(VARIANT_BOOL* pVal)
{
	begin_function
	*pVal = CParameterMap::IsBlank() ? VARIANT_TRUE : VARIANT_FALSE;
	end_function
}

STDMETHODIMP CParameter::IsNumeric(VARIANT_BOOL* pVal)
{
	begin_function	
	*pVal = IsDouble() ? VARIANT_TRUE : VARIANT_FALSE;
	end_function
}

STDMETHODIMP CParameter::put_Columns(long newVal)
{
	return CParameterMap::SetColumns(newVal);
}

STDMETHODIMP CParameter::put_Name(BSTR newVal)
{
	if (!::SysStringLen(newVal)) return CParameterMap::ReturnErrorR(IDS_INVALID_NAME);	
	return m_szName.Set(newVal);	
}

STDMETHODIMP CParameter::get_Value(VARIANT *pVal)
{
	if (CParameterMap::IsBlank()){
		return CComVariant().Detach(pVal);
	} else if (m_nRows == 1 && m_nCols == 1){	
		// scalar (literally - i.e. don't use the IsScalar() method in case it does something flashy!)
		return GetValue(0, 0, pVal);
	} else {
		return GetValue(pVal);
	}
}

STDMETHODIMP CParameter::put_AutoGrow(VARIANT_BOOL newVal)
{
	SetAutoGrow(newVal ? true : false);
	return S_OK;
}

STDMETHODIMP CParameter::put_Element(long Row, long Column, VARIANT newVal)
{
	return CParameterMap::SetValue(Row - 1, Column - 1, newVal);	
}

STDMETHODIMP CParameter::put_Rows(long newVal)
{
	return CParameterMap::SetRows(newVal);
}

STDMETHODIMP CParameter::put_Value(VARIANT newVal)
{	
	return SetValue(newVal);	
}

STDMETHODIMP CParameter::RemoveColumn(long Column)
{
	return CParameterMap::RemoveColumn(Column - 1);
}

STDMETHODIMP CParameter::RemoveRow(long Row)
{
	return CParameterMap::RemoveRow(Row - 1);
}

STDMETHODIMP CParameter::ReplaceBlanksWithSpaces(void)
{
	CParameterMap::ReplaceBlanksWithSpaces();
	return S_OK;
}

STDMETHODIMP CParameter::Transpose()
{
	return CParameterMap::Transpose();	
}
