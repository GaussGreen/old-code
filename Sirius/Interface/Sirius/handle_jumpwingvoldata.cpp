//	jumpwingvoldata.cpp : Implementation of CJumpWingVolData
//
/////////////////////////////////////////////////////////////////////////////

#include "stdafx.h"
#include "handle_jumpwingvoldata.h"
#include "handle_strikes.h"
#include "siriusapplication.h"

STDMETHODIMP CJumpWingVolData::InterfaceSupportsErrorInfo(REFIID riid)
{
	static const IID* arr[] = { &IID_IJumpWingVolData };
	for (int i=0; i < sizeof(arr) / sizeof(arr[0]); i++){
		if (InlineIsEqualGUID(*arr[i],riid)) return S_OK;
	}
	return S_FALSE;
}

STDMETHODIMP CJumpWingVolData::get_JumpWingParameters(VARIANT* pVal)
{
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

STDMETHODIMP CJumpWingVolData::get_Value(VARIANT *pVal)
{
	HRESULT								hr;
	CParameterMap						pmStrikes;
	
	if (hr = pmStrikes.SetValue(m_aspStrike)) return hr;
	return CParameterMap::VariableArgumentListToArray(pVal, 2, pmStrikes, m_pmJumpWingCoefficients);
}

STDMETHODIMP CJumpWingVolData::put_Value(VARIANT newVal)
{
	/*parameter list is 0 - Strikes
						1 - JumpWingCoefficients*/
	
	HRESULT											hr;
	std::vector<CComVariant>						vv;	
	vector<JumpWingCoeffs>							jw;	

	begin_function
	if (hr = CParameterMap::ArrayToVector(newVal, &vv, NULL)) return hr;
	if (vv.size() != 2) return CParameterMap::ReturnErrorR(IDS_NUMBER_PARAMETERS, IID_IJumpWingVolData);
	map_strikes_vector_parameter(vv[0], hStrikes);
	map_parameter(vv[1], CParameterMap, pmJumpWingCoefficients);
			
	if (pmJumpWingCoefficients.GetCols() != 5) return CParameterMap::ReturnErrorR(IDS_COLUMNS_INVALID, IID_IJumpWingVolData);		
	jw.resize(pmJumpWingCoefficients.GetRows());
	for (int n = 0; n < pmJumpWingCoefficients.GetRows(); n++){
		pmJumpWingCoefficients.GetValue(n, 0, &jw[n].ATMVolF);
//		pmJumpWingCoefficients.GetValue(n, 1, &jw[n].a);
		pmJumpWingCoefficients.GetValue(n, 1, &jw[n].skew);
		pmJumpWingCoefficients.GetValue(n, 2, &jw[n].CallWing);
		pmJumpWingCoefficients.GetValue(n, 3, &jw[n].PutWing);
//		pmJumpWingCoefficients.GetValue(n, 4, &jw[n].ASym);
		pmJumpWingCoefficients.GetValue(n, 4, &jw[n].volMin);

	}
	m_h->initialize(hStrikes, jw);
	m_aspStrike = _shStrikes;
	m_pmJumpWingCoefficients.Attach(&pmJumpWingCoefficients);
	end_function
}
