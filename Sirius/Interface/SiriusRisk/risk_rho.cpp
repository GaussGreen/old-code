// risk_rho.cpp : Implementation of CRHo
#include "stdafx.h"
#include "SiriusRisk.h"
#include "risk_rho.h"

/////////////////////////////////////////////////////////////////////////////
// CRHo

STDMETHODIMP CRho::InterfaceSupportsErrorInfo(REFIID riid)
{
	static const IID* arr[] = 
	{
		&IID_IRho
	};
	for (int i=0; i < sizeof(arr) / sizeof(arr[0]); i++)
	{
		if (InlineIsEqualGUID(*arr[i],riid))
			return S_OK;
	}
	return S_FALSE;
}




STDMETHODIMP CRho::Calculate(IEvaluatable* Portfolio, IResult* pVal)
{
	begin_function
	CComPtr<IResult> spResult(pVal) ;// = GetResultObject(&pVal);
	double f = Rho(m_fShift ,Portfolio, spResult);
	return spResult->AddValue(L"Rho", CComVariant(f));
	end_function
}

HRESULT CRho::FinalConstruct()
{
	m_fShift = 0.01;
	return S_OK;
}

STDMETHODIMP CRho::get_Shift(double* pVal)
{
	*pVal = m_fShift;
	return S_OK;
}


STDMETHODIMP CRho::put_Shift(double newVal)
{
	m_fShift = newVal;
	return S_OK;
}

