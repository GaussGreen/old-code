//	Puff.cpp : Implementation of CPuff
//
/////////////////////////////////////////////////////////////////////////////

#include "stdafx.h"
#include "APACProducts.h"
#include "Puff.h"

STDMETHODIMP CPuff::InterfaceSupportsErrorInfo(REFIID riid)
{
	static const IID* arr[] = { &IID_IPuff };
	for (int i=0; i < sizeof(arr) / sizeof(arr[0]); i++){
		if (InlineIsEqualGUID(*arr[i],riid)) return S_OK;
	}
	return S_FALSE;
}


STDMETHODIMP CPuff::Evaluate(IResult *pVal)
{
	begin_function

	begin_calculate(pVal, L"Price")
		double f = m_h->GetPrice();
		pVal->AddValue(L"Price", CComVariant(f));
	end_calculate

	begin_calculate(pVal, L"Horse")
		pVal->AddValue(L"Horse", CComVariant("Jimmy"));
	end_calculate

	end_function
}
