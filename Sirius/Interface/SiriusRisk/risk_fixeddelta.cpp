// risk_FixedDelta.cpp : Implementation of CFixedDelta
#include "stdafx.h"
#include "SiriusRisk.h"
#include "risk_FixedDelta.h"

/////////////////////////////////////////////////////////////////////////////
// CFixedDelta

STDMETHODIMP CFixedDelta::InterfaceSupportsErrorInfo(REFIID riid)
{
	static const IID* arr[] = 
	{
		&IID_IFixedDelta
	};
	for (int i=0; i < sizeof(arr) / sizeof(arr[0]); i++)
	{
		if (InlineIsEqualGUID(*arr[i],riid))
			return S_OK;
	}
	return S_FALSE;
}

HRESULT CFixedDelta::FinalConstruct()
{
	m_fShift = 0.01;
	return S_OK;
}

STDMETHODIMP CFixedDelta::Calculate(IEvaluatable* Portfolio, IResult* pVal)
{
	begin_function
	HRESULT hr = S_OK;
	CComPtr<IResult> spResult(pVal) ;// = GetResultObject(&pVal);
	double f = Delta(m_fShift ,Portfolio, spResult, true);
	hr =  spResult->AddValue(L"FixedDelta", CComVariant(f));


	return hr ;
	end_function
}




STDMETHODIMP CFixedDelta::get_Shift(double* pVal)
{
	*pVal = m_fShift;
	return S_OK;
}


STDMETHODIMP CFixedDelta::put_Shift(double newVal)
{
	m_fShift = newVal;
	return S_OK;
}
