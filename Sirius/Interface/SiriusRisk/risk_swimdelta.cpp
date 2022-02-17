// risk_SwimDelta.cpp : Implementation of CSwimDelta
#include "stdafx.h"
#include "SiriusRisk.h"
#include "risk_SwimDelta.h"

/////////////////////////////////////////////////////////////////////////////
// CSwimDelta


STDMETHODIMP CSwimDelta::Calculate(IEvaluatable* Portfolio, IResult* pVal)
{
	begin_function
	CComPtr<IResult> spResult(pVal) ;// = GetResultObject(&pVal);
	double f = 0.0;//= Delta(m_fShift , Portfolio, spResult);
	return spResult->AddValue(L"SwimDelta", CComVariant(f));
	end_function
}

HRESULT CSwimDelta::FinalConstruct()
{
	m_fShift = 0.01;
	return S_OK;
}

STDMETHODIMP CSwimDelta::get_Shift(double* pVal)
{
	*pVal = m_fShift;
	return S_OK;
}

STDMETHODIMP CSwimDelta::InterfaceSupportsErrorInfo(REFIID riid)
{
	static const IID* arr[] = { &IID_ISwimDelta};
	for (int i=0; i < sizeof(arr) / sizeof(arr[0]); i++){
		if (InlineIsEqualGUID(*arr[i],riid)) return S_OK;
	}
	return S_FALSE;
}

STDMETHODIMP CSwimDelta::put_Shift(double newVal)
{
	m_fShift = newVal;
	return S_OK;
}
