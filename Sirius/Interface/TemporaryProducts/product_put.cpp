//	product_put.cpp : Implementation of CProductPut
//
/////////////////////////////////////////////////////////////////////////////

#include "stdafx.h"
#include "TemporaryProducts.h"
#include "product_put.h"

STDMETHODIMP CProductPut::Evaluate(IResult* pVal)
{
	begin_function	
	HRESULT								hr;	
								
	if (!m_hUnderlying) throw "No underlying defined";
	
	long nToday = m_hUnderlying->GetDateHandle()->GetDate();
	double fForward = m_hUnderlying->GetQuantoForward(nToday, m_dateMaturity, false);
	double fVolatility = m_hUnderlying->GetCompositeVolatility(MlEqStrike(m_fStrike, m_hUnderlying->GetDateHandle()), nToday, m_dateMaturity);	
	double fTime = m_hUnderlying->GetVolatilityStructure()->getDateToDouble()->GetYearFraction(m_dateMaturity);				
	double d1 = (log(fForward / m_fStrike) + 0.5 * fVolatility * fVolatility * fTime) / (fVolatility * sqrt(fTime));	
	double d2 = (log(fForward / m_fStrike) - 0.5 * fVolatility * fVolatility * fTime) / (fVolatility * sqrt(fTime));
	double fDiscountFactor = m_hUnderlying->GetPayZeroCurve(true)->GetDiscountFactor(nToday, m_dateMaturity);
					
	begin_calculate(pVal, L"Price")
		double f = fDiscountFactor * (m_fStrike * ::normal(-d2) - fForward * ::normal(-d1));
		if (hr = pVal->AddValue(L"Price", CComVariant(f))) return hr;
	end_calculate

	begin_calculate(pVal, L"Delta")
	if (hr = pVal->AddValue(L"Delta", CComVariant(::normal(d1) - 1.0))) return hr;
	end_calculate
		
	end_function
}

HRESULT CProductPut::FinalConstruct(void)
{
	m_fStrike = 1.0;
	m_dateMaturity = MlEqDate(MlEqDate::GetCurrentDate()).AddTenor("1y")->GetDate();
	return S_OK;
}

STDMETHODIMP CProductPut::InterfaceSupportsErrorInfo(REFIID riid)
{
	static const IID* arr[] = { &IID_IMattOption };
	for (int i=0; i < sizeof(arr) / sizeof(arr[0]); i++){
		if (InlineIsEqualGUID(*arr[i],riid)) return S_OK;
	}
	return S_FALSE;
}