
#include "stdafx.h"
#include "handle_arpropvoldata.h"
#include "handle_date.h"

/////////////////////////////////////////////////////////////////////////////
//
// CARPropVolData
//
/////////////////////////////////////////////////////////////////////////////

STDMETHODIMP CARPropVolData::InterfaceSupportsErrorInfo(REFIID riid)
{
	static const IID* arr[] = { &IID_IARPropVolData };
	for (int i=0; i < sizeof(arr) / sizeof(arr[0]); i++){
		if (InlineIsEqualGUID(*arr[i],riid)) return S_OK;
	}
	return S_FALSE;
}

STDMETHODIMP CARPropVolData::get_Value(VARIANT *pVal)
{	
	return S_OK;
}

HRESULT CARPropVolData::FinalConstruct(void)
{
	return S_OK;
}

STDMETHODIMP CARPropVolData::put_Value(VARIANT newVal)
{
	/*parameter list is 0 - AtmVolLambda
						1 - AtmVolS0,
						2 - AtmVolSInfty,
						3 - SkewVolLambda,
						4 - SkewVolS0,
						5 - SkewVolInfty,	
						6 - curvature,
						7 - timeCutoff,
						8 - fwds,
						9 - dateToDouble, 
						10- dates	*/

	HRESULT											hr;
	std::vector<CComVariant>						vv;

	begin_function	
	if (hr = CParameterMap::ArrayToVector(newVal, &vv, NULL)) return hr;
	if (vv.size() != 11) return CParameterMap::ReturnErrorR(IDS_NUMBER_PARAMETERS, IID_IHullVolData);

	map_parameter(vv[0], double, AtmVolLambda);
	map_parameter(vv[1], double, AtmVolS0);
	map_parameter(vv[2], double, AtmVolSInfty);
	map_parameter(vv[3], double, SkewVolLambda);
	map_parameter(vv[4], double, SkewVolS0);
	map_parameter(vv[5], double, SkewVolInfty);
	map_parameter(vv[6], double, curvature);
	map_parameter(vv[7], long, timeCutoff);
	map_parameter(vv[8], CVector, fwds);
	map_object_parameter(vv[9],  Date, hDateToDouble);
	map_parameter(vv[10], vector<long>, dates);
		
	m_h->initialize(AtmVolLambda, AtmVolS0, AtmVolSInfty, SkewVolLambda, SkewVolS0, SkewVolInfty, curvature, timeCutoff, fwds, hDateToDouble, dates);

	end_function
	return S_OK;
}
