//	product_generic.cpp : Implementation of CGeneric
//
/////////////////////////////////////////////////////////////////////////////

#include "stdafx.h"
#include "Products.h"
#include "product_generic.h"

STDMETHODIMP CGeneric::Evaluate(IResult* pVal)
{
	HRESULT								hr;
	
	begin_function	
	begin_calculate(pVal, L"Price")
		if (hr = pVal->AddValue(L"Price", CComVariant(m_fPrice))) return hr;
	end_calculate
	begin_calculate(pVal, L"Delta")
		if (hr = pVal->AddValue(L"Delta", CComVariant(m_fDelta))) return hr;
	end_calculate
	begin_calculate(pVal, L"Gamma")
		if (hr = pVal->AddValue(L"Gamma", CComVariant(m_fGamma))) return hr;
	end_calculate
	begin_calculate(pVal, L"Vega")
		if (hr = pVal->AddValue(L"Vega", CComVariant(m_fVega))) return hr;
	end_calculate
	begin_calculate(pVal, L"FXDelta")
		if (hr = pVal->AddValue(L"FXDelta", CComVariant(m_fFXDelta))) return hr;
	end_calculate
	end_function
}

HRESULT CGeneric::FinalConstruct(void)
{
	m_fPrice = 0.0;
	m_fDelta = 0.0;
	m_fGamma = 0.0;
	m_fVega = 0.0;
	m_fFXDelta = 0.0;
	m_szRamModel = "???";
	m_cePay = NoCurrency;
	m_szUnderlying = "???";
	return S_OK;
}

STDMETHODIMP CGeneric::InterfaceSupportsErrorInfo(REFIID riid)
{
	static const IID* arr[] = { &IID_IGeneric };
	for (int i=0; i < sizeof(arr) / sizeof(arr[0]); i++){
		if (InlineIsEqualGUID(*arr[i],riid)) return S_OK;
	}
	return S_FALSE;
}