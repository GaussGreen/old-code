//	handle_ramvoldata.cpp : Implementation of CRamVolData
//
/////////////////////////////////////////////////////////////////////////////

#include "stdafx.h"
#include "handle_ramvoldata.h"
#include "handle_strikes.h"
#include "siriusapplication.h"


STDMETHODIMP CRamVolData::get_BasicParameters(VARIANT* pVal)
{
	begin_function
	CParameterMap				pm;
	const vector<RamCoeffs>&	vrc = m_h->getRamParameters();
	
	// return Base, DownCutoff, DownLevel, UpLevel, UpCutoff
	pm.SetSize(vrc.size(), 5);
	for (long nRow = 0; nRow < vrc.size(); nRow++){
		pm.SetValue(nRow, 0, vrc[nRow].BaseVol);
		pm.SetValue(nRow, 1, vrc[nRow].DownCutoff->m_strike);
		pm.SetValue(nRow, 2, vrc[nRow].Vol90Strike - vrc[nRow].BaseVol);
		pm.SetValue(nRow, 3, vrc[nRow].Vol110Strike - vrc[nRow].BaseVol);
		pm.SetValue(nRow, 4, vrc[nRow].UpCutoff->m_strike);
	}
	return pm.GetValue(pVal);
	end_function
}


STDMETHODIMP CRamVolData::get_Parameters(VARIANT* pVal)
{
	begin_function
	CParameterMap				pm;
	const vector<RamCoeffs>&	vrc = m_h->getRamParameters();
	
	// return Base, Slope, Curvature, Down, Up
	pm.SetSize(vrc.size(), 5);
	for (long nRow = 0; nRow < vrc.size(); nRow++){
		pm.SetValue(nRow, 0, vrc[nRow].BaseVol);
		pm.SetValue(nRow, 1, vrc[nRow].Slope);
		pm.SetValue(nRow, 2, vrc[nRow].Curvature);
		pm.SetValue(nRow, 3, vrc[nRow].DownCutoff->m_strike);
		pm.SetValue(nRow, 4, vrc[nRow].UpCutoff->m_strike);
	}
	return pm.GetValue(pVal);
	end_function
}

STDMETHODIMP CRamVolData::get_Value(VARIANT *pVal)
{	
	CParameterMap		pmStrikes;			pmStrikes.SetValue(m_aspStrike);
	CParameterMap		pmBaseVolArr;		pmBaseVolArr.SetValue(m_BaseVolArr);
	CParameterMap		pmSlopeArr;			pmSlopeArr.SetValue(m_SlopeArr);
	CParameterMap		pmCurvatureArr;		pmCurvatureArr.SetValue(m_CurvatureArr);
	CParameterMap		pmDownCutoffArr;	pmDownCutoffArr.SetValue(m_aspDownCutoffArr);
	CParameterMap		pmUpCutoffArr;		pmUpCutoffArr.SetValue(m_aspUpCutoffArr);
			
	return CParameterMap::VariableArgumentListToArray(pVal, 6, pmStrikes, pmBaseVolArr, pmSlopeArr, pmCurvatureArr, pmDownCutoffArr, pmUpCutoffArr);
}

STDMETHODIMP CRamVolData::InterfaceSupportsErrorInfo(REFIID riid)
{
	static const IID* arr[] = { &IID_IRamVolData };
	for (int i=0; i < sizeof(arr) / sizeof(arr[0]); i++){
		if (InlineIsEqualGUID(*arr[i],riid)) return S_OK;
	}
	return S_FALSE;
}

STDMETHODIMP CRamVolData::put_Value(VARIANT newVal)
{		
	/*parameter list is vpm[0] - Strikes
						vpm[1] - BaseVolArr
						vpm[2] - SlopeArr
						vpm[3] - CurvatureArr
						vpm[4] - DownCutoffArr
						vpm[5] - UpCutoffArr*/		

	HRESULT											hr; 
	std::vector<CComVariant>						vv;
	std::vector<CParameterMap>						vpm;
	std::vector<RamCoeffs>							rm;
	
	begin_function
	if (hr = CParameterMap::ArrayToVector(newVal, &vv, &vpm)) return hr;
	if (vv.size() != 6) return CParameterMap::ReturnErrorR(IDS_NUMBER_PARAMETERS, IID_IRamVolData);
	map_strikes_vector_parameter(vv[0], hStrikes);
	map_parameter(vv[1], CVector, BaseVolArr);
	map_parameter(vv[2], CVector, SlopeArr);
	map_parameter(vv[3], CVector, CurvatureArr);			
	map_strikes_vector_parameter(vv[4], hDownStrikes);	
	map_strikes_vector_parameter(vv[5], hUpStrikes);
	int numberSlices = BaseVolArr.getsize();
	if (BaseVolArr.getsize() != SlopeArr.getsize() || BaseVolArr.getsize() != CurvatureArr.getsize() ) return CParameterMap::ReturnErrorR(IDS_VECTOR_SIZES_DIFFER, IID_IRamVolData);

	rm.resize(BaseVolArr.getsize());
	for (int i = 0; i < rm.size(); i++ )
	{
		rm[i].BaseVol		= getObjectFromCVector(BaseVolArr,i);
		rm[i].Slope			= getObjectFromCVector(SlopeArr,i);
		rm[i].Curvature		= getObjectFromCVector(CurvatureArr,i);
		rm[i].DownCutoff	= getObjectFromVectorVector(hDownStrikes, i, 0, numberSlices);
		rm[i].UpCutoff		= getObjectFromVectorVector(hUpStrikes, i, 0, numberSlices);		
		rm[i].m_strike		= getObjectFromVectorVector(hStrikes, i, 0, numberSlices);	
	}

	m_h->initialize(hStrikes, rm);
	m_aspStrike			=	_shStrikes;
	m_BaseVolArr		=	BaseVolArr;
	m_SlopeArr			=	SlopeArr;
	m_CurvatureArr		=	CurvatureArr;
	m_aspDownCutoffArr	=	_shDownStrikes;
	m_aspUpCutoffArr	=	_shUpStrikes;
	m_h->PutNumberOfSlices(numberSlices);
	end_function
}