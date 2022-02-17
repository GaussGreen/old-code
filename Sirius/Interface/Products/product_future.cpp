//	product_future.cpp : Implementation of CForward
//
/////////////////////////////////////////////////////////////////////////////

#include "stdafx.h"
#include "Products.h"
#include "product_future.h"

HRESULT CFuture::FinalConstruct(void)
{
	m_fStrike = 0.0;
	m_dateMaturity = 0.0;
	m_bTotalReturn = false;
	return S_OK;
}

STDMETHODIMP CFuture::Evaluate(IResult *pVal)
{
	begin_function
	HRESULT								hr;	

	if (!m_hUnderlying) throw "No underlying defined";	
	
	long								nToday = m_hUnderlying->GetDateHandle()->GetDate();

	if (nToday >= m_dateMaturity){
		// We discount nToday == m_dateMaturity since we only really care about precise valuation in the afternoon after which
		// back office have closed out the posisiton.
		begin_calculate(pVal, L"Price")
			if (hr = pVal->AddValue(L"Price", CComVariant(0.0))) return hr;
		end_calculate
	} else {
		double fForward = m_hUnderlying->GetQuantoForward(nToday, m_dateMaturity, m_bTotalReturn);	
		double fTime = m_hUnderlying->GetDateHandle()->GetYearFraction(m_dateMaturity);
		//double fDiscountFactor = m_hUnderlying->GetPayZeroCurve(true)->GetDiscountFactor(nToday, m_dateMaturity);
		
		begin_calculate(pVal, L"Price")
			double fPrice = (fForward - m_fStrike);
			if (hr = pVal->AddValue(L"Price", CComVariant(fPrice))) return hr;
		end_calculate	
	}
		
	end_function
}

STDMETHODIMP CFuture::InterfaceSupportsErrorInfo(REFIID riid)
{
	static const IID* arr[] = { &IID_IARPropFitVolData };
	for (int i=0; i < sizeof(arr) / sizeof(arr[0]); i++){
		if (InlineIsEqualGUID(*arr[i],riid)) return S_OK;
	}
	return S_FALSE;
}