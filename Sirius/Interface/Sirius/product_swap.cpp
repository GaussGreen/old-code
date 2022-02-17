//	product_swap.cpp : Implementation of CSwap
//
/////////////////////////////////////////////////////////////////////////////

#include "stdafx.h"
#include "product_swap.h"
#include "MlEqArray.h"
#include "MlEqZeroCurve.h"


//////////////////////////////////////////////////////////////////////////////
//	Evaluate
//
//	Required pricing function
//
STDMETHODIMP CSwap::Evaluate(IResult* pVal)
{
	begin_function
	HRESULT								hr;
	double								f;

	begin_calculate(pVal, L"Price");
		// ToDo - only allow swap price if the zero curves are the same.
		f = m_h->GetPV(A_Leg) + m_h->GetPV(B_Leg);
		if (hr = pVal->AddValue(L"Price", CComVariant(f))) return hr;
	end_calculate

	begin_calculate(pVal, L"A_PV");
		f = m_h->GetPV(A_Leg);
		if (hr = pVal->AddValue(L"A_PV", CComVariant(f))) return hr;
	end_calculate

	begin_calculate(pVal, L"B_PV");
		f = m_h->GetPV(B_Leg);
		if (hr = pVal->AddValue(L"B_PV", CComVariant(f))) return hr;
	end_calculate

	begin_calculate(pVal, L"A_FixedPV");
		f = m_h->GetFixedPV(A_Leg);
		if (hr = pVal->AddValue(L"A_FixedPV", CComVariant(f))) return hr;
	end_calculate

	begin_calculate(pVal, L"B_FixedPV");
		f = m_h->GetFixedPV(B_Leg);
		if (hr = pVal->AddValue(L"B_FixedPV", CComVariant(f))) return hr;
	end_calculate

	begin_calculate(pVal, L"A_FloatingPV");
		f = m_h->GetFloatingPV(A_Leg);
		if (hr = pVal->AddValue(L"A_FloatingPV", CComVariant(f))) return hr;
	end_calculate

	begin_calculate(pVal, L"B_FloatingPV");
		f = m_h->GetFloatingPV(B_Leg);
		if (hr = pVal->AddValue(L"B_FloatingPV", CComVariant(f))) return hr;
	end_calculate

	end_function
}


//////////////////////////////////////////////////////////////////////////////
//	Get[...]PV
//
//	Returns present value of cash flows associated with the Swap.
//
STDMETHODIMP CSwap::GetFixedPV(SwapLegEnum SwapLeg, double* pVal)
{
	begin_function
	*pVal = m_h->GetFixedPV(SwapLeg);
	end_function
}

STDMETHODIMP CSwap::GetFloatingPV(SwapLegEnum SwapLeg, double* pVal)
{
	begin_function
	*pVal = m_h->GetFloatingPV(SwapLeg);
	end_function
}

STDMETHODIMP CSwap::GetPV(SwapLegEnum SwapLeg, double* pVal)
{
	begin_function
	*pVal = m_h->GetPV(SwapLeg);
	end_function
}


//////////////////////////////////////////////////////////////////////////////
//	InterfaceSupportsErrorInfo
//
//	Allows association of an ErrorInfo object with this class
//
STDMETHODIMP CSwap::InterfaceSupportsErrorInfo(REFIID riid)
{
	static const IID* arr[] = { &IID_ISwap };
	for (int i=0; i < sizeof(arr) / sizeof(arr[0]); i++){
		if (InlineIsEqualGUID(*arr[i],riid)) return S_OK;
	}
	return S_FALSE;
}


