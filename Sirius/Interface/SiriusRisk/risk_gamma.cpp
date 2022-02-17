//	risk_gamma.cpp : Implementation of CGamma
//
/////////////////////////////////////////////////////////////////////////////

#include "stdafx.h"
#include "SiriusRisk.h"
#include "risk_gamma.h"

STDMETHODIMP CGamma::Calculate(IEvaluatable* Portfolio, IResult* pVal)
{
	begin_function

	CComPtr<IResult> spResult(pVal) ;//GetResultObject(&pVal);

	double fUnderPrice =  UnderlyingPrice(Portfolio, spResult);

	
	double f = Price(Portfolio, spResult);
	double fHigh = PriceSpot(m_fShift, Portfolio, spResult, false,"");
	double fLow = PriceSpot(-m_fShift, Portfolio, spResult, false,"");

	double fGamma = (fLow + fHigh - 2 * f ) /  (fUnderPrice * m_fShift); 
	return spResult->AddValue(L"Gamma", CComVariant(fGamma));
	end_function		
}

HRESULT CGamma::FinalConstruct()
{
	m_fShift = 0.01;
	return S_OK;
}

STDMETHODIMP CGamma::get_Shift(double* pVal)
{
	*pVal = m_fShift;
	return S_OK;
}

STDMETHODIMP CGamma::InterfaceSupportsErrorInfo(REFIID riid)
{
	static const IID* arr[] = { &IID_IGamma };
	for (int i=0; i < sizeof(arr) / sizeof(arr[0]); i++){
		if (InlineIsEqualGUID(*arr[i],riid)) return S_OK;
	}
	return S_FALSE;
}

STDMETHODIMP CGamma::put_Shift(double newVal)
{
	m_fShift = newVal;
	return S_OK;
}
