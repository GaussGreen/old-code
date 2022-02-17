//	handle_hullvoldata.cpp : Implementation of CHullVolData
//
/////////////////////////////////////////////////////////////////////////////

#include "stdafx.h"
#include "handle_hullvoldata.h"
#include "handle_strikes.h"
#include "siriusapplication.h"

HRESULT CHullVolData::FinalConstruct(void)
{
	return S_OK;
}

STDMETHODIMP CHullVolData::get_Value(VARIANT *pVal)
{	
	HRESULT								hr;
	CParameterMap						pmStrikes;
	
	if (hr = pmStrikes.SetValue(m_aspStrike)) return hr;
	return CParameterMap::VariableArgumentListToArray(pVal, 2, pmStrikes, m_pmHullCoefficients);
}	

STDMETHODIMP CHullVolData::InterfaceSupportsErrorInfo(REFIID riid)
{
	static const IID* arr[] = { &IID_IHullVolData };
	for (int i=0; i < sizeof(arr) / sizeof(arr[0]); i++){
		if (InlineIsEqualGUID(*arr[i],riid)) return S_OK;
	}
	return S_FALSE;
}
	
STDMETHODIMP CHullVolData::put_Value(VARIANT newVal)
{	
	/*parameter list is vpm[0] - Strikes
						vpm[1] - HullParameters*/
	
	HRESULT											hr;
	std::vector<CComVariant>						vv;	
	vector<HullCoeffs>								hl;

	begin_function	
	if (hr = CParameterMap::ArrayToVector(newVal, &vv, NULL)) return hr;
	if (vv.size() != 2) return CParameterMap::ReturnErrorR(IDS_NUMBER_PARAMETERS, IID_IHullVolData);
	map_strikes_vector_parameter(vv[0], hStrikes);
	map_parameter(vv[1], CParameterMap, pmHullCoefficients);
		
		
	if (pmHullCoefficients.GetCols() != 6) return CParameterMap::ReturnErrorR(IDS_COLUMNS_INVALID, IID_IHullVolData);
	hl.resize(pmHullCoefficients.GetRows());
	for (int n = 0; n < hl.size(); n++){
		pmHullCoefficients.GetValue(n, 0, &hl[n].VoltMin);	 
		pmHullCoefficients.GetValue(n, 1, &hl[n].NsAtmVoltMin);			
		pmHullCoefficients.GetValue(n, 2, &hl[n].RightInflectionPoint);
		pmHullCoefficients.GetValue(n, 3, &hl[n].CallWing);
		pmHullCoefficients.GetValue(n, 4, &hl[n].LeftInflectionPoint);
		pmHullCoefficients.GetValue(n, 5, &hl[n].PutWing);
	}	
	m_h->initialize(hStrikes, hl);

	m_aspStrike = _shStrikes;
	m_pmHullCoefficients.Attach(&pmHullCoefficients);
	end_function
}