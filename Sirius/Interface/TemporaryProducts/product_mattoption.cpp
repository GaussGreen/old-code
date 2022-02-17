//	product_mattoption.cpp : Implementation of CMattOption
//
/////////////////////////////////////////////////////////////////////////////

#include "stdafx.h"
#include "TemporaryProducts.h"
#include "product_mattoption.h"

STDMETHODIMP CMattOption::Evaluate(IResult *pVal)
{
	begin_function
	HRESULT								hr;	
						
	if (!m_hUnderlying1 || !m_hUnderlying2) throw "No underlying defined";
	long nToday = m_hUnderlying1->GetDateHandle()->GetDate();	
	double fForward1 = m_hUnderlying1->GetQuantoForward(nToday , m_dateMaturity, false);		
	double fVolatility1 = m_hUnderlying1->GetCompositeVolatility(MlEqStrike(m_fStrike1), nToday , m_dateMaturity);
	double fForward2 = m_hUnderlying2->GetQuantoForward(nToday , m_dateMaturity, false);
	double fVolatility2 = m_hUnderlying2->GetCompositeVolatility(MlEqStrike(m_fStrike2), nToday , m_dateMaturity);			
	double fTime1 = m_hUnderlying1->GetVolatilityStructure()->getDateToDouble()->GetYearFraction(m_dateMaturity);
	double fTime2 = m_hUnderlying2->GetVolatilityStructure()->getDateToDouble()->GetYearFraction(m_dateMaturity);
				
	begin_calculate(pVal, L"Price")
		double fDiscountFactor1 = m_hUnderlying1->GetPayZeroCurve(true)->GetDiscountFactor(nToday, m_datePay);
		double fDiscountFactor2 = m_hUnderlying2->GetPayZeroCurve(true)->GetDiscountFactor(nToday, m_datePay);
		double f1 = ::Bs(fForward1, fVolatility1, fTime1, m_fStrike1, fDiscountFactor1, 1);
		double f2 = ::Bs(fForward2, fVolatility2, fTime2, m_fStrike2, fDiscountFactor2, 1);		
		double fPrice = m_fNotional1 * f1 + m_fNotional2 * f2;
		if (hr = pVal->AddValue(L"Price", CComVariant(fPrice))) return hr;
	end_calculate
		
	end_function
}

HRESULT CMattOption::FinalConstruct(void)
{	
	m_fStrike1 = 0.0;	
	m_fStrike2 = 0.0;		
	m_datePay = m_dateMaturity = MlEqDate(MlEqDate::GetCurrentDate()).AddTenor("1y")->GetDate();
	m_fNotional1 = 1.0;	
	m_fNotional2 = 1.0;	
	return S_OK;
}

STDMETHODIMP CMattOption::InterfaceSupportsErrorInfo(REFIID riid)
{
	static const IID* arr[] = { &IID_IMattOption };
	for (int i=0; i < sizeof(arr) / sizeof(arr[0]); i++){
		if (InlineIsEqualGUID(*arr[i],riid)) return S_OK;
	}
	return S_FALSE;
}