//	risk_delta.cpp : Implementation of CDelta
//
/////////////////////////////////////////////////////////////////////////////

#include "stdafx.h"
#include "SiriusRisk.h"
#include "risk_gamma5.h"

STDMETHODIMP CGamma5::Calculate(IEvaluatable* Portfolio, IResult* pVal)
{
	begin_function
	CComPtr<IResult> spResult(pVal) ;// = GetResultObject(&pVal);
	
	double fUnderPrice	= UnderlyingPrice(Portfolio, spResult);
	double f			= Price(Portfolio, spResult);
	double fHigh1		= PriceSpot(0.01, Portfolio, spResult, false,"");
	double fHigh5		= PriceSpot(0.05, Portfolio, spResult, false,"");
	double fLow1		= PriceSpot(-0.01, Portfolio, spResult, false,"");
	double fHigh4		= PriceSpot(0.0395, Portfolio, spResult, false,"");
	double fHigh6		= PriceSpot(0.0605, Portfolio, spResult, false,"");

	double fHighDelta	= 	(fHigh6 - fHigh4)  /  (fUnderPrice * 0.0210 ); 
	double fLowDelta	= 	(fHigh1 - fLow1)  /  (fUnderPrice * 0.02); 

	double fGamma5 = (fHighDelta -fLowDelta	 	)  ; 
	return spResult->AddValue(L"Gamma5", CComVariant(fGamma5 ));

	end_function
}

HRESULT CGamma5::FinalConstruct()
{
	return S_OK;
}

STDMETHODIMP CGamma5::get_Shift(double* pVal)
{
	//*pVal = m_fShift;
	return S_OK;
}

STDMETHODIMP CGamma5::InterfaceSupportsErrorInfo(REFIID riid)
{
	static const IID* arr[] = { &IID_IGamma5};
	for (int i=0; i < sizeof(arr) / sizeof(arr[0]); i++){
		if (InlineIsEqualGUID(*arr[i],riid)) return S_OK;
	}
	return S_FALSE;
}

STDMETHODIMP CGamma5::put_Shift(double newVal)
{
	//m_fShift = newVal;
	return S_OK;
}
