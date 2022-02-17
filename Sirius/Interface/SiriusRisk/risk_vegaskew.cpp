// risk_VegaSkew.cpp : Implementation of CVegaSkew
#include "stdafx.h"
#include "SiriusRisk.h"
#include "risk_VegaSkew.h"

/////////////////////////////////////////////////////////////////////////////
// CVegaSkew

STDMETHODIMP CVegaSkew::Calculate(IEvaluatable* Portfolio, IResult* pVal)
{
	begin_function

	CComPtr<IResult> spResult(pVal) ;   

	double f = VegaSkew(Portfolio, pVal, m_dtExpiryDate);

	return spResult->AddValue(L"VegaSkew", CComVariant(f));
	end_function
}

HRESULT CVegaSkew::FinalConstruct()
{
	CComVariant vExpiryDate;
	CComPtr<IDate> sptodaysDate ; 
	long dt;

	//Get a year form today as by  Default for ExpiryDate
	HRESULT hr = sptodaysDate.CoCreateInstance(L"Sirius.Date");
	
	hr = sptodaysDate->get_SerialNumber(&dt);
	MLEDate(CComVariant(dt ),CComVariant(L"1Y"),&vExpiryDate);
	hr = vExpiryDate.ChangeType(VT_R8);
	m_dtExpiryDate = (vExpiryDate.dblVal);
	return S_OK;
}

STDMETHODIMP CVegaSkew::get_ExpiryDate(DATE* pVal)
{
	*pVal = m_dtExpiryDate;
	return S_OK;
}

STDMETHODIMP CVegaSkew::InterfaceSupportsErrorInfo(REFIID riid)
{
	static const IID* arr[] = { &IID_IVegaSkew };
	for (int i=0; i < sizeof(arr) / sizeof(arr[0]); i++){
		if (InlineIsEqualGUID(*arr[i],riid)) return S_OK;
	}
	return S_FALSE;
}

STDMETHODIMP CVegaSkew::put_ExpiryDate(DATE newVal)
{
	m_dtExpiryDate = newVal;
	return S_OK;
}
