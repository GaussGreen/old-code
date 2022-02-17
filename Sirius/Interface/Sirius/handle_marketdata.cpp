//	handle_marketdata.cpp : Implementation of CMarketData
//
/////////////////////////////////////////////////////////////////////////////

#include "stdafx.h"
#include "handle_marketdata.h"
#include "siriusapplication.h"
#include "MlEqObjects.h"

HRESULT CMarketData::FinalConstruct(void)
{
	HRESULT								hr;
	
	// Create instances of the CComPtr<> member variables.
	if (hr = m_spCorrelationMatrices.CoCreateInstance(CLSID_CorrelationMatrices)) return hr;
	if (hr = m_spDividendSchedules.CoCreateInstance(CLSID_DividendSchedules)) return hr;
	if (hr = m_spSpotSchedules.CoCreateInstance(CLSID_SpotSchedules)) return hr;
	if (hr = m_spVolatilityStructures.CoCreateInstance(CLSID_VolatilityStructures)) return hr;		
	if (hr = m_spZeroCurves.CoCreateInstance(CLSID_ZeroCurves)) return hr;
	if (hr = m_spAssets.CoCreateInstance(CLSID_Assets)) return hr;
	return S_OK;
}

STDMETHODIMP CMarketData::InterfaceSupportsErrorInfo(REFIID riid)
{
	static const IID* arr[] = {&IID_IMarketData};
	for (int i=0; i < sizeof(arr) / sizeof(arr[0]); i++){
		if (InlineIsEqualGUID(*arr[i],riid)) return S_OK;
	}
	return S_FALSE;
}