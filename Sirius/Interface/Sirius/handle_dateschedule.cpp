//	handle_dateschedule.cpp : Implementation of CDateSchedule
//
/////////////////////////////////////////////////////////////////////////////

#include "stdafx.h"
#include "handle_dateschedule.h"

HRESULT CDateSchedule::FinalConstruct(void)
{
	return S_OK;
}

STDMETHODIMP CDateSchedule::get_NumberAt(DATE Date, double* pVal)
{		
	std::vector<double>					vector;
	std::string							szDate;
				
	begin_function	
	vector = m_h->GetValueAt(Date);		
	if (vector.size() != 1){
		CParameterMap::DateToString(Date, &szDate);
		return CParameterMap::ReturnErrorRS(IDS_NOT_NUMBER_AT_DATE, szDate);
	}
	*pVal = vector[0];
	end_function
}

STDMETHODIMP CDateSchedule::put_NumberAt(DATE Date, double newVal)
{
	std::vector<double>					vector;
	
	vector.push_back(newVal);
	(*m_h)[Date] = vector;
	return S_OK;
}

STDMETHODIMP CDateSchedule::get_Value(VARIANT *pVal)
{
	CParameterMap					pmValues;		
	CParameterMap					pmRet;
	HRESULT							hr;
	long							nRow(0);	
	
	if (!m_h->size()) return S_OK;				
	pmRet.SetSize(m_h->size(), 1);
	for (std::map<long, std::vector<double> >::const_iterator itCount = m_h->begin(); itCount != m_h->end(); ++itCount){
		CParameterMap		pm;		
		// add the date
		if (hr = pmRet.SetValue(nRow, 0, itCount->first)) return hr;
		// add the value
		if (hr = pm.SetValue(itCount->second)) return hr;
		pm.Transpose();
		pmValues.AddToEnd(pm);
		nRow++;
	}		
	pmRet.AddToRHS(pmValues);
	return pmRet.GetValue(pVal);	
}

STDMETHODIMP CDateSchedule::get_ValueAt(DATE Date, VARIANT *pVal)
{
	HRESULT								hr;	
	CParameterMap						pm;	
	
	std::map<long, std::vector<double> >::const_iterator it = m_h->find(Date);
	if (it == m_h->end()){
		std::string szDate;
		CParameterMap::DateToString(Date, &szDate);
		return CParameterMap::ReturnErrorRS(IDS_NO_VALUE_AT_DATE, szDate);
	}
	if (hr = pm.SetValue(it->second)) return hr;
	return pm.GetValue(pVal);
}

STDMETHODIMP CDateSchedule::InterfaceSupportsErrorInfo(REFIID riid)
{
	static const IID* arr[] = { &IID_IDateSchedule };
	for (int i=0; i < sizeof(arr) / sizeof(arr[0]); i++){
		if (InlineIsEqualGUID(*arr[i],riid)) return S_OK;
	}
	return S_FALSE;
}

STDMETHODIMP CDateSchedule::put_Value(VARIANT newVal)
{
	HRESULT								hr;
	std::vector<CParameterMap>			vpm;	

	m_h->clear();
	if (hr = CParameterMap::ArrayToVector(newVal, NULL, &vpm, false/*the decomposition is performed below*/)) return hr;
	for (long n = 1; n < vpm.size(); n++){
		vpm[0].AddToRHS(vpm[n]);
	}

	for (long nRow = 0; nRow < vpm[0].GetRows(); nRow++){
		CParameterMap				pm;
		long						nDate;
		std::vector<double>			vector;
		if (vpm[0].GetRow(nRow, &pm)) continue;
		
		if (pm.GetValue(0, 0, &nDate)){
			std::string sz;
			pm.GetValue(0, 0, &sz);
			return CParameterMap::ReturnErrorRS(IDS_INVALID_DATE, sz, IID_IDateSchedule);
		}
		pm.RemoveColumn(0);
		if (!pm.IsBlank() && !pm.IsDouble(-1, -1, IDS_INVALID_PARAMETER_VALUE)) return E_FAIL;		
		pm.GetValue(&vector);		
		(*m_h)[nDate] = vector;
	}
	return S_OK;
}

STDMETHODIMP CDateSchedule::put_ValueAt(DATE Date, VARIANT newVal)
{
	HRESULT								hr;
	CParameterMap						pm;	
	std::vector<double>					vector;
		
	if (hr = pm.SetValue(newVal)) return hr;
	if (hr = pm.GetValue(&vector)) return hr;
	(*m_h)[Date] = vector;
	return S_OK;	
}
