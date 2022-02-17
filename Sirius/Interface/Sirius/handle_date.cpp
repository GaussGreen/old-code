//	handle_date.cpp : Implementation of CDate
//
/////////////////////////////////////////////////////////////////////////////

#include "stdafx.h"
#include "handle_date.h"
#include "siriusapplication.h"

STDMETHODIMP CDate::Add(long Days, long Months, long Years, BusinessDayConventionEnum BusinessDayConvention, VARIANT CalendarOpt, VARIANT_BOOL PreserveEoM, IDate** pVal)
{
	begin_function
	map_optional_parameter(CalendarOpt, estring, szCalendar, "");
	map_parameter(PreserveEoM, bool, bPreserveEoM);	
	map_analytic_to_com(m_h->Add(Days, Months, Years, BusinessDayConvention, szCalendar, bPreserveEoM), Date, spDate);
	return spDate.CopyTo(pVal);
	end_function
}

STDMETHODIMP CDate::AddTenor(BSTR Tenor, BusinessDayConventionEnum BusinessDayConvention, VARIANT CalendarOpt, VARIANT_BOOL PreserveEoM, IDate** pVal)
{
	begin_function	
	map_parameter(Tenor, estring, szTenor);
	map_optional_parameter(CalendarOpt, estring, szCalendar, "");
	map_parameter(PreserveEoM, bool, bPreserveEoM);	
	map_analytic_to_com(m_h->AddTenor(szTenor, BusinessDayConvention, szCalendar, bPreserveEoM), Date, spDate);
	return spDate.CopyTo(pVal);
	end_function
}

STDMETHODIMP CDate::CopyToClipboard()
{			
	HRESULT								hr;
	CComVariant							v;

	if (hr = get_Value(&v)) return hr;	
	unmap_parameter(v, pm);	
	return pm.CopyToClipboard();
}

STDMETHODIMP CDate::get_DayCountConvention(DayCountConventionEnum *pVal)
{
	begin_function
	*pVal = m_h->GetDayCountConvention();
	end_function
}

STDMETHODIMP CDate::get_DayOfMonth(long* pVal)
{
	begin_function
	*pVal = m_h->GetDayOfMonth();
	end_function
}

STDMETHODIMP CDate::get_DayOfWeek(WeekdayEnum* pVal)
{
	begin_function
	*pVal = m_h->GetDayOfWeek();
	end_function
}

STDMETHODIMP CDate::get_DaysInYear(long* pVal)
{
	begin_function
	*pVal = m_h->GetDaysInYear();
	end_function
}

STDMETHODIMP CDate::get_Month(MonthEnum* pVal)
{
	begin_function
	*pVal = m_h->GetMonth();
	end_function
}

STDMETHODIMP CDate::get_SerialNumber(long* pVal)
{
	begin_function
	*pVal = m_h->GetDate();
	end_function
}

STDMETHODIMP CDate::get_Value(VARIANT *pVal)
{
	return g_pApplication->GetObjectManager().ImplementGetValue(this, pVal);
}

STDMETHODIMP CDate::get_Year(long* pVal)
{	
	begin_function
	*pVal = m_h->GetYear();
	end_function
}

STDMETHODIMP CDate::GetYearFraction(long Date, double* pVal)
{
	begin_function
	*pVal = m_h->GetYearFraction(Date);
	end_function
}

STDMETHODIMP CDate::InterfaceSupportsErrorInfo(REFIID riid)
{
	static const IID* arr[] = { &IID_IDate };
	for (int i=0; i < sizeof(arr) / sizeof(arr[0]); i++){
		if (InlineIsEqualGUID(*arr[i],riid)) return S_OK;
	}
	return S_FALSE;
}

STDMETHODIMP CDate::put_DayCountConvention(DayCountConventionEnum newVal)
{
	begin_function
	m_h->PutDayCountConvention(newVal);	
	end_function
}

STDMETHODIMP CDate::put_SerialNumber(long newVal)
{
	begin_function
	m_h->PutDate(newVal);
	end_function
}

STDMETHODIMP CDate::put_Value(VARIANT newVal)
{
	begin_function
	HRESULT								hr;		
	std::vector<CParameterMap>			vpm;	
	
	if (hr = CParameterMap::ArrayToVector(newVal, NULL, &vpm)) return hr;
	if (vpm.size() == 1 && vpm[0].GetCols() == 2){
		// Use the general function		
		hr = g_pApplication->GetObjectManager().ImplementPutValue(newVal, this);
		return hr;
	} else {		
		if (vpm.size() != 2) return CParameterMap::ReturnErrorR(IDS_NUMBER_PARAMETERS);		
		map_parameter(vpm[0], long, Date);
		map_enum_parameter(vpm[1], DayCountConventionEnum, DayCountConvention);				
		m_h = new MlEqDate(Date, DayCountConvention);
	}
	end_function
}