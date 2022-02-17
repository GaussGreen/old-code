//	handle_ramfitvoldata.cpp : Implementation of CRamFitVolData
//
/////////////////////////////////////////////////////////////////////////////

#include "stdafx.h"
#include "handle_ramfitvoldata.h"
#include "handle_strikes.h"
#include "siriusapplication.h"
#include "handle_volatilitystructure.h"


STDMETHODIMP CRamFitVolData::InterfaceSupportsErrorInfo(REFIID riid)
{
	static const IID* arr[] = { &IID_IRamFitVolData };
	for (int i=0; i < sizeof(arr) / sizeof(arr[0]); i++){
		if (InlineIsEqualGUID(*arr[i],riid)) return S_OK;
	}
	return S_FALSE;
}

STDMETHODIMP CRamFitVolData::get_Value(VARIANT *pVal)
{	
	CParameterMap		pmVolatility;				pmVolatility.SetToObjectPointer(m_spVolatilityStructure);
	CParameterMap		pmFittingStrikes;			pmFittingStrikes.SetValue(m_aspFittingStrikes);
	CParameterMap		pmUpStrikes;				pmUpStrikes.SetValue(m_aspUpStrikes);
	CParameterMap		pmDownStrikes;				pmDownStrikes.SetValue(m_aspDownStrikes);
	CParameterMap		pmFittingDates;				pmFittingDates.SetValue(m_nFittingDates);	
	CParameterMap		pmFitBoundaries;			pmFitBoundaries.SetValue(m_h->m_fitBoundaries);
	CParameterMap		pmUseVega;					pmUseVega.SetValue(m_h->m_useVega);
	CParameterMap		pmFitSpread;				pmFitSpread.SetValue(m_h->m_fitSpread);
	CParameterMap		pmFitMidVols;				pmFitMidVols.SetValue(m_h->m_fitMidVols);
		
	return CParameterMap::VariableArgumentListToArray(pVal, 9, pmVolatility, pmFittingStrikes, pmUpStrikes, pmDownStrikes, pmFittingDates, pmFitBoundaries, pmUseVega, pmFitSpread, pmFitMidVols);
}

STDMETHODIMP CRamFitVolData::put_Value(VARIANT newVal)
{	
	/*parameter list is vv[0] - VolatilityStructure
						vv[1] - FittingStrikes
						vv[2] - UpStrikes
						vv[3] - DownStrikes
						vv[4] - FittingDates		
						vv[5] - FitBoundaries		
						vv[6] - UseVega	
						vv[7] - FitSpread		
						vv[8] - FitMidVols		*/

	HRESULT											hr;
	std::vector<CComVariant>						vv;	
		
	begin_function
	if (hr = CParameterMap::ArrayToVector(newVal, &vv, NULL)) return hr;
	if (vv.size() != 9) return CParameterMap::ReturnErrorR(IDS_NUMBER_PARAMETERS, IID_IRamVolData);

	map_object_parameter(vv[0], VolatilityStructure, hVolatilityStructure);
	map_strikes_vector_parameter(vv[1], hFittingStrikes);	
	map_strikes_vector_parameter(vv[2], hUpStrikes);	
	map_strikes_vector_parameter(vv[3], hDownStrikes);	
	map_parameter(vv[4], std::vector<long>, FittingDates);
	map_parameter(vv[5], long, FitBoundaries);
	map_parameter(vv[6], long, UseVega);
	map_parameter(vv[7], double, FitSpread);
	map_parameter(vv[8], double, FitMidVols);
		
	m_h->initialize(hVolatilityStructure, hFittingStrikes, FittingDates, CMatrix(), FitBoundaries, UseVega, FitSpread, FitMidVols, hUpStrikes, hDownStrikes);

	m_spVolatilityStructure = _shVolatilityStructure;
	m_aspFittingStrikes = _shFittingStrikes;
	m_aspUpStrikes = _shUpStrikes;
	m_aspDownStrikes = _shDownStrikes;
	m_nFittingDates = FittingDates;		
	end_function
}