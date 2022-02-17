//	handle_scenario.cpp : Implementation of CScenario
//
/////////////////////////////////////////////////////////////////////////////

#include "stdafx.h"
#include "handle_scenario.h"

STDMETHODIMP CScenario::get_Value(VARIANT* pVal)
{	
	begin_function
	m_h->GetValue((CComVariant*)pVal);
	end_function
}

STDMETHODIMP CScenario::InterfaceSupportsErrorInfo(REFIID riid)
{
	static const IID* arr[] = { &IID_IScenario };
	for (int i=0; i < sizeof(arr) / sizeof(arr[0]); i++){
		if (InlineIsEqualGUID(*arr[i],riid)) return S_OK;
	}
	return S_FALSE;
}

STDMETHODIMP CScenario::put_Value(VARIANT newVal)
{	
	begin_function
	m_h->PutValue(newVal);
	end_function
}