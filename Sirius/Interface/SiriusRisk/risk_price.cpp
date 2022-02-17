//	risk_price.cpp : Implementation of CPrice
//
/////////////////////////////////////////////////////////////////////////////

#include "stdafx.h"
#include "SiriusRisk.h"
#include "risk_price.h"

STDMETHODIMP CPrice::Calculate(IEvaluatable* Portfolio, IResult* pVal)
{
	begin_function
	CComPtr<IResult> spResult(pVal) ;// = GetResultObject(&pVal);
	double f = Price(Portfolio, spResult);
	return spResult->AddValue(L"Price", CComVariant(f));
	end_function
}

STDMETHODIMP CPrice::InterfaceSupportsErrorInfo(REFIID riid)
{
	static const IID* arr[] = { &IID_IPrice };
	for (int i=0; i < sizeof(arr) / sizeof(arr[0]); i++){
		if (InlineIsEqualGUID(*arr[i],riid)) return S_OK;
	}
	return S_FALSE;
}
