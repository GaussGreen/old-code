//	handle_parameterlist.cpp : Implementation of CParameterList
//
/////////////////////////////////////////////////////////////////////////////

#include "stdafx.h"
#include "handle_parameterlist.h"

STDMETHODIMP CParameterList::AddValue(BSTR Name, VARIANT Value)
{	
	begin_function
	m_h->AddValue(Name, Value);
	end_function
}

HRESULT CParameterList::FinalConstruct(void)
{		
	return S_OK;
}

STDMETHODIMP CParameterList::get_Count(long *pVal)
{
	begin_function
	*pVal = m_h->GetCount();
	end_function
}

STDMETHODIMP CParameterList::GetValue(BSTR Name, VARIANT *pVal)
{	
	begin_function
	m_h->MlEqDictionary<CComVariant>::GetValue(Name, (CComVariant*)pVal, true);
	end_function
}

STDMETHODIMP CParameterList::get_Value(VARIANT *pVal)
{	
	begin_function
	m_h->MlEqDictionary<CComVariant>::GetValue((CComVariant*)pVal);	
	end_function
}

STDMETHODIMP CParameterList::InterfaceSupportsErrorInfo(REFIID riid)
{
	static const IID* arr[] = {&IID_IParameterList};
	for (int i=0; i < sizeof(arr) / sizeof(arr[0]); i++){
		if (InlineIsEqualGUID(*arr[i],riid)) return S_OK;
	}
	return S_FALSE;
}

STDMETHODIMP CParameterList::put_Value(VARIANT newVal)
{
	begin_function
	m_h->PutValue(newVal);
	end_function
}
