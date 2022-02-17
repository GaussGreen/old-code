//	risk_delta.cpp : Implementation of CDelta
//
/////////////////////////////////////////////////////////////////////////////

#include "stdafx.h"
#include "SiriusRisk.h"
#include "risk_delta.h"

STDMETHODIMP CDelta::Calculate(IEvaluatable* Portfolio, IResult* pVal)
{
	begin_function
	HRESULT hr = S_OK;
	m_mapPartialGreeks.clear();
	CComPtr<IResult> spResult(pVal) ;// = GetResultObject(&pVal);
	double f = Delta(m_fShift ,Portfolio, spResult, false);
	
	if(m_bPartialGreek) {
		CComPtr<IAssets> spAssets;
		CComPtr<IAssets> spBaseAssets;
		if( hr = Portfolio->GetUnderlyings(&spAssets)) propagate_error_ex(hr);
		if( hr = spAssets->GetBaseUnderlyings(&spBaseAssets)) propagate_error_ex(hr);
		long nAssets = 0;
		
		if( hr =  spBaseAssets->get_Count(&nAssets) ) propagate_error_ex(hr);
		for ( long idx = 1 ; idx <= nAssets; idx++) {
			CComPtr<IAsset> spAsset;
			CComBSTR sAsset;
			if ( hr = spBaseAssets->get_Item(_variant_t(idx), &spAsset)) propagate_error_ex(hr);
			if ( hr = spAsset->get_Name(&sAsset) ) propagate_error_ex(hr);
			double f = Delta(m_fShift ,Portfolio, spResult, false, estring(sAsset));
			m_mapPartialGreeks.insert(std::pair<std::string,double>(estring(sAsset), f));
		}
	}

	return hr;
	end_function
}

HRESULT CDelta::FinalConstruct()
{
	m_fShift = 0.01;
	m_bPartialGreek = false;

	return S_OK;
}

STDMETHODIMP CDelta::get_Shift(double* pVal)
{
	
	*pVal = m_fShift;
	return S_OK;
	

}

STDMETHODIMP CDelta::InterfaceSupportsErrorInfo(REFIID riid)
{
	static const IID* arr[] = { &IID_IDelta };
	for (int i=0; i < sizeof(arr) / sizeof(arr[0]); i++){
		if (InlineIsEqualGUID(*arr[i],riid)) return S_OK;
	}
	return S_FALSE;
}

STDMETHODIMP CDelta::put_Shift(double newVal)
{
	m_fShift = newVal;
	return S_OK;
}


STDMETHODIMP CDelta::get_CalculatePartialGreek(/*[out, retval] */ VARIANT_BOOL* pVal)
{
	*pVal = (m_bPartialGreek) ? VARIANT_TRUE : VARIANT_FALSE; 

	return S_OK;
}

STDMETHODIMP CDelta::put_CalculatePartialGreek(/*[in]*/ VARIANT_BOOL newVal)
{
	m_bPartialGreek = (newVal == VARIANT_TRUE) ;
	return S_OK;
}

STDMETHODIMP CDelta::GetCalculatedPartialGreekResults(/*[in]*/ IResult* Result)
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

