//	handle_arpropfitvoldata.cpp : Implementation of CARPropFitVolData
//
/////////////////////////////////////////////////////////////////////////////

#include "stdafx.h"
#include "handle_arpropfitvoldata.h"
#include "handle_volatilitystructure.h"
#include "handle_strikes.h"
#include "handle_date.h"


STDMETHODIMP CARPropFitVolData::InterfaceSupportsErrorInfo(REFIID riid)
{
	static const IID* arr[] = { &IID_IARPropFitVolData };
	for (int i=0; i < sizeof(arr) / sizeof(arr[0]); i++){
		if (InlineIsEqualGUID(*arr[i],riid)) return S_OK;
	}
	return S_FALSE;
}

STDMETHODIMP CARPropFitVolData::get_Value(VARIANT *pVal)
{
	// TODO: Add your implementation code here
	return S_OK;
}

STDMETHODIMP CARPropFitVolData::put_Value(VARIANT newVal)
{
	/*parameter list is 
			vv[0], VolatilityStructure
			vv[1], FittingStrikes
			vv[2], DateToDouble
			vv[3], timeCutoff
			vv[4], fittingDates
			vv[5], fwds   
			vv[6], curvature   */
 
	HRESULT											hr;
	std::vector<CComVariant>						vv;	

	begin_function	
	if (hr = CParameterMap::ArrayToVector(newVal, &vv, NULL)) return hr;
	if (vv.size() != 7) return CParameterMap::ReturnErrorR(IDS_NUMBER_PARAMETERS, IID_IARPropFitVolData);

	map_object_parameter(vv[0], VolatilityStructure, hVolatilityStructure);
	map_strikes_vector_parameter(vv[1], hFittingStrikes);
	map_object_parameter(vv[2], Date, hDateToDouble);
	map_parameter(vv[3], long, timeCutoff);
	map_parameter(vv[4], vector<long>, fittingDates);
	map_parameter(vv[5], CVector, fwds);
	map_parameter(vv[6], double, curvature);

	CMatrix BoundsOnConstraints;
	m_h->initialize(hVolatilityStructure, hFittingStrikes, fittingDates, fwds, hDateToDouble, BoundsOnConstraints, timeCutoff, curvature);

	end_function
	return S_OK;
}
