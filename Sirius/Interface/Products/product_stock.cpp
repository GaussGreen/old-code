// product_stock.cpp : Implementation of CStock
#include "stdafx.h"
#include "Products.h"
#include "product_stock.h"

/////////////////////////////////////////////////////////////////////////////
// CStock

STDMETHODIMP CStock::InterfaceSupportsErrorInfo(REFIID riid)
{
	static const IID* arr[] = 
	{
		&IID_IStock
	};
	for (int i=0; i < sizeof(arr) / sizeof(arr[0]); i++)
	{
		if (InlineIsEqualGUID(*arr[i],riid))
			return S_OK;
	}
	return S_FALSE;
}


STDMETHODIMP CStock::Evaluate(IResult* pVal)
{
	begin_function	
	if (!m_hUnderlying) throw "No underlying defined";

	HRESULT		hr;	

	long	nDate = m_hUnderlying->GetDateHandle()->GetDate();
	double	fSpot = m_hUnderlying->GetSpot(nDate);

	begin_calculate(pVal, L"Price")
		if (hr = pVal->AddValue(L"Price", CComVariant(fSpot))) return hr;			
	end_calculate	
	return S_OK;

	end_function
}
