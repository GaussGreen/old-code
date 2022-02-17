// risk_Theta.cpp : Implementation of CTheta
#include "stdafx.h"
#include "SiriusRisk.h"
#include "risk_Theta.h"

/////////////////////////////////////////////////////////////////////////////
// CTheta


STDMETHODIMP CTheta::Calculate(IEvaluatable* Portfolio, IResult* pVal)
{
	begin_function
	CComPtr<IResult> spResult (pVal);// = GetResultObject(&pVal);
	double f = Price(Portfolio, spResult);

	double fTheta = PriceTheta (m_fShift, Portfolio, spResult);

	// ToDo - Probably need to move marketData date back

	return spResult->AddValue(L"Theta", CComVariant((fTheta - f) / m_fShift ));

	end_function		
}

HRESULT CTheta::FinalConstruct()
{
	m_fShift = 1;
	return S_OK;
}

STDMETHODIMP CTheta::get_Shift(double* pVal)
{
	*pVal = m_fShift;
	return S_OK;
}

STDMETHODIMP CTheta::InterfaceSupportsErrorInfo(REFIID riid)
{
	static const IID* arr[] = { &IID_IVega };
	for (int i=0; i < sizeof(arr) / sizeof(arr[0]); i++){
		if (InlineIsEqualGUID(*arr[i],riid)) return S_OK;
	}
	return S_FALSE;
}

STDMETHODIMP CTheta::put_Shift(double newVal)
{
	m_fShift = newVal;
	return S_OK;
}
