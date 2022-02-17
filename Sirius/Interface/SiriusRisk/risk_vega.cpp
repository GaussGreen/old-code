//	risk_vega.cpp : Implementation of CVega
//
/////////////////////////////////////////////////////////////////////////////

#include "stdafx.h"
#include "SiriusRisk.h"
#include "risk_vega.h"

STDMETHODIMP CVega::Calculate(IEvaluatable* Portfolio, IResult* pVal)
{
	begin_function
	HRESULT hr = S_OK;
	CComPtr<IResult> spResult (pVal);// = GetResultObject(&pVal);
	double f = Price(Portfolio, spResult);
	double fHigh = PriceVolatility(m_fShift, Portfolio, spResult,"");
	hr = spResult->AddValue(L"Vega", CComVariant((fHigh - f) * (0.01/m_fShift) ));



	if(m_bPartialGreek) {
		CComPtr<IVolatilityStructures> spVolatilityStructures;
		
		if( hr = Portfolio->GetVolatilityStructures(&spVolatilityStructures)) propagate_error_ex(hr);
	
		long nCount= 0;
		
		if( hr =  spVolatilityStructures->get_Count(&nCount) ) propagate_error_ex(hr);
		for ( long idx = 1 ; idx <= nCount; idx++) {
			CComPtr<IVolatilityStructure> spVol;
			CComBSTR sItem;
			if ( hr = spVolatilityStructures->get_Item(_variant_t(idx), &spVol)) propagate_error_ex(hr);
			if ( hr = spVol->get_Name(&sItem) ) propagate_error_ex(hr);
			double pf = PriceVolatility(m_fShift ,Portfolio, spResult, estring(sItem));
			pf= ( (pf- f) * (0.01/m_fShift) );
			m_mapPartialGreeks.insert(std::pair<std::string,double>(estring(sItem), pf));
		}
	}


	return hr;


	end_function		
}

HRESULT CVega::FinalConstruct()
{
	m_fShift = 0.01;
	m_bPartialGreek = false;
	return S_OK;
}

STDMETHODIMP CVega::get_Shift(double* pVal)
{
	*pVal = m_fShift;
	return S_OK;
}

STDMETHODIMP CVega::InterfaceSupportsErrorInfo(REFIID riid)
{
	static const IID* arr[] = { &IID_IVega };
	for (int i=0; i < sizeof(arr) / sizeof(arr[0]); i++){
		if (InlineIsEqualGUID(*arr[i],riid)) return S_OK;
	}
	return S_FALSE;
}

STDMETHODIMP CVega::put_Shift(double newVal)
{
	m_fShift = newVal;
	return S_OK;
}




STDMETHODIMP CVega::get_CalculatePartialGreek(/*[out, retval] */ VARIANT_BOOL* pVal)
{
	*pVal = (m_bPartialGreek) ? VARIANT_TRUE : VARIANT_FALSE; 

	return S_OK;
}

STDMETHODIMP CVega::put_CalculatePartialGreek(/*[in]*/ VARIANT_BOOL newVal)
{
	m_bPartialGreek = (newVal == VARIANT_TRUE) ;
	return S_OK;
}

STDMETHODIMP CVega::GetCalculatedPartialGreekResults(/*[in]*/ IResult* Result)
{
	begin_function
	HRESULT					hr = S_OK;
	CComPtr<IResult>		spResult(Result);
	for ( std::map<std::string,double>::const_iterator it =m_mapPartialGreeks.begin();  it != m_mapPartialGreeks.end(); it++) {
		_bstr_t sName(it->first.c_str());
		if ( hr = Result->AddValue(sName , _variant_t(it->second))) propagate_error_ex(hr);
	}
	return S_OK;
	end_function
}