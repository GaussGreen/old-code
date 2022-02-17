// DivDelta.cpp : Implementation of CDivDelta
#include "stdafx.h"
#include "SiriusRisk.h"
#include "risk_DivDelta.h"

/////////////////////////////////////////////////////////////////////////////
// CDivDelta


STDMETHODIMP CDivDelta::Calculate(IEvaluatable* Portfolio, IResult* pVal)
{
	begin_function
	CComPtr<IResult> spResult(pVal) ;   
	double f = DivDelta(m_fShift, Portfolio, spResult);
	return spResult->AddValue(L"DivDelta", CComVariant(f));
	end_function
}

HRESULT CDivDelta::FinalConstruct()
{
	m_fShift = -0.01;
	return S_OK;
}

STDMETHODIMP CDivDelta::get_Shift(double* pVal)
{
	*pVal = m_fShift;
	return S_OK;
}

STDMETHODIMP CDivDelta::InterfaceSupportsErrorInfo(REFIID riid)
{
	static const IID* arr[] = { &IID_IDivDelta};
	for (int i=0; i < sizeof(arr) / sizeof(arr[0]); i++){
		if (InlineIsEqualGUID(*arr[i],riid)) return S_OK;
	}
	return S_FALSE;
}

STDMETHODIMP CDivDelta::put_Shift(double newVal)
{
	m_fShift = newVal;
	return S_OK;
}
