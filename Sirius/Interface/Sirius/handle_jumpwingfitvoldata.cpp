//	handle_jumpwingfitvoldata.cpp : Implementation of CJumpWingFitVolData
//
/////////////////////////////////////////////////////////////////////////////

#include "stdafx.h"
#include "handle_jumpwingfitvoldata.h"


HRESULT CJumpWingFitVolData::FinalConstruct(void)
{	
	return S_OK;
}

STDMETHODIMP CJumpWingFitVolData::get_JumpWingParameters(VARIANT* pVal)
{
	// ToDo - the code here is exactly the same as in CJumpWingVolData::get_JumpWingParameters(VARIANT* pVal);
	// this will disappear once this handle is a special case of CJumpWingVolData
	
	begin_function	
		
	int									nSlices = m_h->getNumberOfSlices();
	CMatrix								jwv(nSlices,7);	
	const vector<JumpWingCoeffs>		jw = m_h->getJumpWingParameters();

	for (int n = 0 ; n < nSlices; n++){
		jwv[n][0] = jw[n].ATMVolF;
		jwv[n][1] = jw[n].skew;
		jwv[n][2] = jw[n].CallWing;
		jwv[n][3] = jw[n].PutWing;
		jwv[n][4] = jw[n].volMin;
		jwv[n][5] = jw[n].a;
		jwv[n][6] = jw[n].ASym;
	}
	
	CParameterMap pm;
	pm.SetValue(jwv);
	return pm.GetValue(pVal);	
	end_function
}

STDMETHODIMP CJumpWingFitVolData::get_Value(VARIANT *pVal)
{
	return CParameterMap::VectorToArray(m_vpm, pVal);	
}

STDMETHODIMP CJumpWingFitVolData::InterfaceSupportsErrorInfo(REFIID riid)
{
	static const IID* arr[] = { &IID_IJumpWingFitVolData };
	for (int i=0; i < sizeof(arr) / sizeof(arr[0]); i++){
		if (InlineIsEqualGUID(*arr[i],riid)) return S_OK;
	}
	return S_FALSE;
}

STDMETHODIMP CJumpWingFitVolData::put_Value(VARIANT newVal)
{

	/*parameter list is vv[0] - MktVols
						vv[1] - JumpWingGuess
						vv[2] - FittingDates
						vv[3] - FitMidPoint
						vv[4] - FitSpread
						vv[5] - UseVega
						*/


	HRESULT											hr;
	std::vector<CComVariant>						vv;	

	begin_function	
	if (hr = CParameterMap::ArrayToVector(newVal, &vv, &m_vpm)) return hr;
	if (vv.size() != 6) throw CStringResource(IDS_NUMBER_PARAMETERS);
	
	map_object_parameter(vv[0], VolatilityStructure, hMktVols);
	map_object_parameter(vv[1], JumpWingVolData, hJumpWingVolDataInitialGuess);	
	map_parameter(vv[2], std::vector<long>, FittingDates);
	map_parameter(vv[3], double, FitMidPoints);
	map_parameter(vv[4], double, FitSpread);
	map_parameter(vv[5], double, UseVega);

	CMatrix BoundsOnConstraints;
	m_h->initialize(hJumpWingVolDataInitialGuess, hMktVols, FittingDates, BoundsOnConstraints, UseVega, FitSpread, FitMidPoints);
	
	end_function	
}
