
#include "stdafx.h"
#include "Products.h"
#include "product_varianceswap.h"

/////////////////////////////////////////////////////////////////////////////
// CVarianceSwapPricer

#include "CVarianceSwap.h"


STDMETHODIMP CVarianceSwapPricer::InterfaceSupportsErrorInfo(REFIID riid)
{
	static const IID* arr[] = { &IID_IARPropFitVolData };
	for (int i=0; i < sizeof(arr) / sizeof(arr[0]); i++){
		if (InlineIsEqualGUID(*arr[i],riid)) return S_OK;
	}
	return S_FALSE;
}


HRESULT CVarianceSwapPricer::FinalConstruct(void)
{	
	m_upCutoff = m_downCutoff = 0.0;
	m_perfType = LogarithmicType ;
	m_boundType = NoStrikesType ;
	m_useCurrentSpot = No;
	m_szCalendar = "EUR" ;
	m_nDaysPerYear = 252 ;

	return S_OK;
}


STDMETHODIMP CVarianceSwapPricer::Evaluate(IResult* pVal)
{
	begin_function

	if (!m_hFixingSchedule) throw "No fixing schedule defined";

	HRESULT		hr;
	CMatrix		mResults;
	std::vector<long>		vDates;

	m_hFixingSchedule->GetDates(vDates);
	long dateStart = vDates[0];
	long dateMaturity = vDates[vDates.size()-1];

	bool bUseCurrentSpot = (m_useCurrentSpot == Yes);

	if( !( m_boundType == ForwardBased || m_boundType == Normalised || m_boundType == Fixed  || m_boundType == NoStrikesType) )
	{
		throw"Wrong bounds type";
	}
	
	RCPtr< CVarianceSwap > pricer = new CVarianceSwap(	m_hUnderlying,
														dateStart,
														dateMaturity,
														m_hFixingSchedule,
														m_szCalendar,
														m_perfType,
														bUseCurrentSpot,
														m_boundType,
														m_upCutoff,
														m_downCutoff,
														m_nDaysPerYear);

	pricer->fComputePrice( mResults );

	long	nToday = m_hUnderlying->GetDateHandle()->GetDate();
	double  df = m_hUnderlying->GetPayZeroCurve(true)->GetDiscountFactor(nToday, m_datePay);


	begin_calculate(pVal, L"Price")
		if (hr = pVal->AddValue(L"Price", CComVariant(df*mResults[0][0]))) return hr;
	end_calculate

	begin_calculate(pVal, L"Var")
		if (hr = pVal->AddValue(L" Var", CComVariant(mResults[0][0]))) return hr;
	end_calculate

	begin_calculate(pVal, L"Var(today, start)")
		if (hr = pVal->AddValue(L" Var(today, start)", CComVariant(mResults[1][0]))) return hr;
	end_calculate

	begin_calculate(pVal, L"Var(today, mat)")
		if (hr = pVal->AddValue(L" Var(today, mat)", CComVariant(mResults[2][0]))) return hr;
	end_calculate

	begin_calculate(pVal, L"Cutoff Up")
		if (hr = pVal->AddValue(L"Cutoff Up", CComVariant(mResults[3][0]))) return hr;
	end_calculate

	begin_calculate(pVal, L"Cutoff Down")
		if (hr = pVal->AddValue(L"Cutoff Down", CComVariant(mResults[3][1]))) return hr;
	end_calculate


	begin_calculate(pVal, L"N perf")
		if (hr = pVal->AddValue(L"N perf", CComVariant(mResults[0][1]))) return hr;
	end_calculate
	
	begin_calculate(pVal, L"N perf(today, mat)")
		if (hr = pVal->AddValue(L"N perf(today, mat)", CComVariant(mResults[2][1]))) return hr;
	end_calculate

	begin_calculate(pVal, L"N perf(today, start)")
		if (hr = pVal->AddValue(L"N perf(today, start)", CComVariant(mResults[1][1]))) return hr;
	end_calculate

	end_function
}
















/*
STDMETHODIMP CVarianceSwapPricer::EvaluateGreeks(IResult* pVal)
{
	begin_function

	if (!m_hFixingSchedule) throw "No fixing schedule defined";

	HRESULT		hr;
	CMatrix		mResults;
	std::vector<long>		vDates;

	m_hFixingSchedule->GetDates(vDates);
	long dateStart = vDates[0];
	long dateMaturity = vDates[vDates.size()-1];

	bool bUseCurrentSpot = (m_useCurrentSpot == Yes);

	if( !( m_boundType == ForwardBased || m_boundType == Normalised || m_boundType == Fixed  || m_boundType == NoStrikesType) )
	{
		throw"Wrong bounds type";
	}
	
	RCPtr< CVarianceSwap > pricer = new CVarianceSwap(	m_hUnderlying,
														dateStart,
														dateMaturity,
														m_hFixingSchedule,
														m_szCalendar,
														m_perfType,
														bUseCurrentSpot,
														m_boundType,
														m_upCutoff,
														m_downCutoff,
														m_nDaysPerYear);

	pricer->fComputePriceAndGreeks( mResults );

	long	nToday = m_hUnderlying->GetDateHandle()->GetDate();
	double  df = m_hUnderlying->GetPayZeroCurve(true)->GetDiscountFactor(nToday, m_datePay);


	begin_calculate(pVal, L"Price")
		if (hr = pVal->AddValue(L"Price", CComVariant(df*mResults[0][0]))) return hr;
	end_calculate

	begin_calculate(pVal, L"Var")
		if (hr = pVal->AddValue(L"Var", CComVariant(mResults[0][0]))) return hr;
	end_calculate

	begin_calculate(pVal, L"Delta")
		if (hr = pVal->AddValue(L"Delta", CComVariant(df*mResults[1][0]))) return hr;
	end_calculate

	begin_calculate(pVal, L"Gamma")
		if (hr = pVal->AddValue(L"Gamma", CComVariant(df*mResults[2][0]))) return hr;
	end_calculate	
		
	begin_calculate(pVal, L"Var(today, start)")
		if (hr = pVal->AddValue(L" Var(today, start)", CComVariant(mResults[0][1]))) return hr;
	end_calculate

	begin_calculate(pVal, L"Delta(today, start)")
		if (hr = pVal->AddValue(L"Delta(today, start)", CComVariant(mResults[1][1]))) return hr;
	end_calculate

	begin_calculate(pVal, L"Gamma (today, start)")
		if (hr = pVal->AddValue(L"Gamma (today, start)", CComVariant(mResults[2][1]))) return hr;
	end_calculate

	begin_calculate(pVal, L"Var(today, mat)")
		if (hr = pVal->AddValue(L" Var(today, mat)", CComVariant(mResults[0][2]))) return hr;
	end_calculate

	begin_calculate(pVal, L"Delta(today, mat)")
		if (hr = pVal->AddValue(L"Delta(today, mat)", CComVariant(mResults[1][2]))) return hr;
	end_calculate

	begin_calculate(pVal, L"Gamma(today, mat)")
		if (hr = pVal->AddValue(L"Gamma(today, mat)", CComVariant(mResults[2][2]))) return hr;
	end_calculate


	end_function
}
*/
