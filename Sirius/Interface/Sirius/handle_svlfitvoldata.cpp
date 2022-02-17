//	handle_svlfitvoldata.cpp : Implementation of CSVLFitVolData
//
/////////////////////////////////////////////////////////////////////////////

#include "stdafx.h"
#include "handle_svlfitvoldata.h"
#include "handle_strikes.h"
#include "siriusapplication.h"

HRESULT CSVLFitVolData::FinalConstruct(void)
{	
	return S_OK;
}

STDMETHODIMP CSVLFitVolData::get_Value(VARIANT *pVal)
{
	CParameterMap		pmStrikes;		pmStrikes.SetValue(m_aspStrike);
	CParameterMap		pmBid;			pmBid.SetValue(m_h->GetBidVols());
	CParameterMap		pmMid;			pmMid.SetValue(m_h->GetMidVols());
	CParameterMap		pmAsk;			pmAsk.SetValue(m_h->GetAskVols());
	CParameterMap		pmTimes;		pmTimes.SetValue(m_h->GetDoubleTime());
	CParameterMap		pmInitialGuess;	pmInitialGuess.SetValue(m_InitialGuess);
	CParameterMap		pmFitIdeals;	pmFitIdeals.SetValue(m_h->m_fitIdeals);
	CParameterMap		pmFitMidPoints;	pmFitMidPoints.SetValue(m_h->m_fitMidVols);
	CParameterMap		pmFitSpread;	pmFitSpread.SetValue(m_h->m_fitSpread);
		
	return CParameterMap::VariableArgumentListToArray(pVal, 9, pmStrikes, pmBid, pmMid, pmAsk, pmTimes, pmInitialGuess, pmFitIdeals, pmFitMidPoints, pmFitSpread);
}

STDMETHODIMP CSVLFitVolData::InterfaceSupportsErrorInfo(REFIID riid)
{
	static const IID* arr[] = { &IID_ISVLFitVolData };
	for (int i=0; i < sizeof(arr) / sizeof(arr[0]); i++){
		if (InlineIsEqualGUID(*arr[i],riid)) return S_OK;
	}
	return S_FALSE;
}

STDMETHODIMP CSVLFitVolData::put_Value(VARIANT newVal)
{	
	/*parameter list is vv[0] - Strikes
						vv[1] - Bid
						vv[2] - Mid
						vv[3] - Ask
						vv[4] - Times
						vv[5] - InitialGuess
						vv[6] - FitIdeals
						vv[7] - FitMidPoints
						vv[8] - FitSpread*/

	HRESULT											hr;
	std::vector<CComVariant>						vv;	

	begin_function	
	if (hr = CParameterMap::ArrayToVector(newVal, &vv, NULL)) return hr;
	if (vv.size() != 9) return CParameterMap::ReturnErrorR(IDS_NUMBER_PARAMETERS, IID_ISVLFitVolData);
	map_strikes_vector_parameter(vv[0], hStrikes);		
	map_parameter(vv[1], vector<CVector>, Bid);
	map_parameter(vv[2], vector<CVector>, Mid);
	map_parameter(vv[3], vector<CVector>, Ask);
	map_parameter(vv[4], CVector, Times);
	map_parameter(vv[5], CMatrix, InitialGuess);
	map_parameter(vv[6], double, FitIdeals);
	map_parameter(vv[7], double, FitMidPoints);
	map_parameter(vv[8], double, FitSpread);
	
	m_h->initialize(Bid, Mid, Ask, hStrikes, Times, FitIdeals, FitMidPoints, FitSpread, InitialGuess, CMatrix());

	m_aspStrike = _shStrikes;
	m_InitialGuess = InitialGuess;
	end_function
}
