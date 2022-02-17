//	handle_multipliervoldata.cpp : Implementation of CMultiplierVolData
//
/////////////////////////////////////////////////////////////////////////////

#include "stdafx.h"
#include "handle_strikes.h"
#include "handle_interpolator.h"
#include "handle_multipliervoldata.h"

STDMETHODIMP CMultiplierVolData::InterfaceSupportsErrorInfo(REFIID riid)
{
	static const IID* arr[] = { &IID_IMultiplierVolData };
	for (int i=0; i < sizeof(arr) / sizeof(arr[0]); i++){
		if (InlineIsEqualGUID(*arr[i],riid)) return S_OK;
	}
	return S_FALSE;
}

STDMETHODIMP CMultiplierVolData::get_Value(VARIANT *pVal)
{	
	begin_function
	HRESULT								hr;
	CParameterMap						pmStrikes;
	CParameterMap						pmInterpolators;
	unmap_parameter(m_h->GetReferenceVol(), pmReferenceVol);
		
	if (hr = pmStrikes.SetValue(m_aspStrike)) return hr;
	if (hr = pmInterpolators.SetValue(m_aspInterpolator)) return hr;
	return CParameterMap::VariableArgumentListToArray(pVal, 3, pmStrikes, pmInterpolators, pmReferenceVol);
	end_function
}

STDMETHODIMP CMultiplierVolData::put_Value(VARIANT newVal)
{
	/*parameter list is 0 - Strikes
						1 - Interpolators
						2 - Reference volatilities */
	
	HRESULT											hr;		
	std::vector<CComVariant>						vv;
	
	begin_function
	if (hr = CParameterMap::ArrayToVector(newVal, &vv, NULL)) return hr;
	if (vv.size() != 3) return CParameterMap::ReturnErrorR(IDS_NUMBER_PARAMETERS, IID_IVolData);	
	map_strikes_vector_parameter(vv[0], hStrikes);
	map_object_vector_parameter(vv[1], Interpolator, hInterpolators);
	map_parameter(vv[2],CVector, refVols);
		
	m_h->initialize(hInterpolators, hStrikes);
	m_h->PutReferenceVol(refVols);
	
	m_aspStrike = _shStrikes;
	m_aspInterpolator = _shInterpolators;
	end_function;
}
