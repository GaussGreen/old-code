//	handle_hermitevoldata.cpp : Implementation of CHermiteVolData
//
/////////////////////////////////////////////////////////////////////////////

#include "stdafx.h"
#include "handle_hermitevoldata.h"

STDMETHODIMP CHermiteVolData::get_Value(VARIANT *pVal)
{
	return CParameterMap::VectorToArray(m_vpm, pVal);
}

STDMETHODIMP CHermiteVolData::InterfaceSupportsErrorInfo(REFIID riid)
{
	static const IID* arr[] = { &IID_IHermiteVolData };
	for (int i=0; i < sizeof(arr) / sizeof(arr[0]); i++){
		if (InlineIsEqualGUID(*arr[i],riid)) return S_OK;
	}
	return S_FALSE;
}

STDMETHODIMP CHermiteVolData::put_Value(VARIANT newVal)
{
	/*parameter list is 0 - StartDateHandle
						1 - DatesArr
						2 - HermiteParameters
						3 - ForwardsArrOpt
	*/
	
	HRESULT											hr;
	std::vector<CParameterMap>						vpm;			
	vector<HermiteCoeff>							hc;

	begin_function
	if (hr = CParameterMap::ArrayToVector(newVal, NULL, &vpm)) return hr;
	if (vpm.size() < 3 || vpm.size() > 4) return CParameterMap::ReturnErrorR(IDS_NUMBER_PARAMETERS, IID_IHermiteVolData);
	vpm.resize(4);

	map_object_parameter(vpm[0], Date, StartDateHandle);
	map_parameter(vpm[1], std::vector<long>, Dates);
	map_parameter(vpm[2], CMatrix, HermiteParameters);
	map_optional_parameter(vpm[3], CVector, ForwardsArr,CVector());
	
	if ( ForwardsArr.getsize() == 0 ){
		ForwardsArr.resize(HermiteParameters.rows());
		for ( int i = 0 ; i < ForwardsArr.getsize(); i++ ){
			ForwardsArr[i] = 1e99;
		}
	}

	// check the sizes
	if (Dates.size() != ForwardsArr.getsize()) throw "The number of elements in the forwards must match the number of dates";
	if (Dates.size() != HermiteParameters.rows()) throw "The number of rows in the HermiteParameters matrix must match the number of dates";
	
	hc.resize(Dates.size());
	for (int nRow = 0; nRow < Dates.size(); nRow++){	
		hc[nRow].fwd = ForwardsArr[nRow];
		hc[nRow].nToday = StartDateHandle;
		hc[nRow].maturity = Dates[nRow];
		hc[nRow].hCoeff = HermiteParameters[nRow];	
		hc[nRow].mat = StartDateHandle->GetYearFraction(Dates[nRow]);
		
	}

	m_h->initialize(hc);
	
	// assign members of this class for get_Value implementation
	m_spStartDate = _sStartDateHandle;
	m_vpm = vpm;			
	end_function
}
