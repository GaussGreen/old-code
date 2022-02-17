// handle_hullandwhite.cpp : Implementation of CHullAndWhite
#include "stdafx.h"
#include "Sirius.h"
#include "handle_hullandwhite.h"
#include "handle_zerocurve.h"
#include "handle_date.h"
#include "MonteCarlo.h"
#include "MlEqPde.h"

#include "mleqobjects.h"

// Franz

/////////////////////////////////////////////////////////////////////////////
// CHullAndWhite

STDMETHODIMP CHullAndWhite::InterfaceSupportsErrorInfo(REFIID riid)
{
	static const IID* arr[] = 
	{
		&IID_IHullAndWhite
	};
	for (int i=0; i < sizeof(arr) / sizeof(arr[0]); i++)
	{
		if (InlineIsEqualGUID(*arr[i],riid))
			return S_OK;
	}
	return S_FALSE;
}

STDMETHODIMP CHullAndWhite::get_Value(VARIANT *pVal)
{

//  TODO: Add your implementation code here

	return S_OK;
}


STDMETHODIMP CHullAndWhite::put_Value(VARIANT newVal)
{
	HRESULT								hr;			
	std::vector<CParameterMap>			vpm;
	
	begin_function
	if (hr = CParameterMap::ArrayToVector(newVal, NULL, &vpm)) return hr;

	if (vpm.size() != 10) return CParameterMap::ReturnErrorR(IDS_NUMBER_PARAMETERS, IID_IInterpolator);

	map_parameter(vpm[0], GVector<long>, lDates);
	map_parameter(vpm[1], double, lambda);

	map_object_parameter(vpm[2], ZeroCurve, hZeroCurve);
	map_object_parameter(vpm[3], Date, hDateToDouble);
	map_parameter(vpm[4],CVector,capletVols);
	map_parameter(vpm[5],CVector,capletStrikes);
	map_parameter(vpm[6],CVector,swaptionVols);
	map_parameter(vpm[7],CVector,swaptionStrikes);
	map_parameter(vpm[8],CMatrix,swaptionDates);
	map_enum_parameter(vpm[9], DayCountConventionEnum, dc);


	/*GVector<long> lDates(Dates.getsize());
	for ( int i = 0 ; i < lDates.getsize(); i++ ){
		lDates[i] = Dates[i];
	}*/

	CVector beta(lDates.getsize());
	CVector gamma(lDates.getsize());

	for ( int i = 0 ; i < lDates.getsize(); i++ ){
		beta[i]		= 1.0;
		gamma[i]	= 0.04;
	}

	
	m_h = new MlEqHullAndWhite(hDateToDouble);
	
	m_h->initialize(lDates, beta, gamma, lambda, hZeroCurve, capletVols, capletStrikes, swaptionVols, swaptionStrikes, swaptionDates, dc);
	end_function
}
