//	handle_voldata.cpp : Implementation of CVolData
//
/////////////////////////////////////////////////////////////////////////////

#include "stdafx.h"
#include "handle_voldata.h"
#include "handle_strikes.h"
#include "handle_interpolator.h"
#include "siriusapplication.h"

HRESULT CVolData::FinalConstruct(void)
{	
	return S_OK;
}

STDMETHODIMP CVolData::InterfaceSupportsErrorInfo(REFIID riid)
{
	static const IID* arr[] = { &IID_IVolData };
	for (int i=0; i < sizeof(arr) / sizeof(arr[0]); i++){
		if (InlineIsEqualGUID(*arr[i],riid)) return S_OK;
	}
	return S_FALSE;
}

STDMETHODIMP CVolData::get_Value(VARIANT *pVal)
{	
	HRESULT								hr;
	CParameterMap						pmStrikes;
	CParameterMap						pmInterpolators;
		
	if (hr = pmStrikes.SetValue(m_aspStrike)) return hr;
	if (hr = pmInterpolators.SetValue(m_aspInterpolator)) return hr;
	return CParameterMap::VariableArgumentListToArray(pVal, 2, pmStrikes, pmInterpolators);
}
																
STDMETHODIMP CVolData::put_Value(VARIANT newVal)
{
	/*parameter list is 0 - Strikes
						1 - Interpolators*/

	HRESULT											hr;
	std::vector<CComVariant>						vv;
	
	begin_function
	if (hr = CParameterMap::ArrayToVector(newVal, &vv, NULL)) return hr;
	if (vv.size() != 2) return CParameterMap::ReturnErrorR(IDS_NUMBER_PARAMETERS, IID_IVolData);
	map_strikes_vector_parameter(vv[0], hStrikes);
	map_object_vector_parameter(vv[1], Interpolator, hInterpolators);
	m_h->initialize(hInterpolators, hStrikes);
	m_aspStrike = _shStrikes;
	m_aspInterpolator = _shInterpolators;
	end_function;
}