//	handle_date.h : Declaration of the CDate
//
/////////////////////////////////////////////////////////////////////////////

#ifndef __DATE_H_
#define __DATE_H_

#include "resource.h"
#include "MlEqDate.h"

class ATL_NO_VTABLE CDate : 
	public CComObjectRootEx<CComSingleThreadModel>,
	public CComCoClass<CDate, &CLSID_Date>,
	public ISupportErrorInfo,
	public IDispatchImpl<IDate, &IID_IDate, &LIBID_Sirius>
{
public:	
	CDate() : m_h(new MlEqDate(MlEqDate::GetCurrentDate())) {}
	DECLARE_REGISTRY_RESOURCEID(IDR_DATE)
	DECLARE_NOT_AGGREGATABLE(CDate)
	DECLARE_PROTECT_FINAL_CONSTRUCT()
	BEGIN_COM_MAP(CDate)
		COM_INTERFACE_ENTRY(IDate)
		COM_INTERFACE_ENTRY(IDispatch)
		COM_INTERFACE_ENTRY(ISupportErrorInfo)
	END_COM_MAP()
	associate_analytic_object(CDate, MlEqDate);
	
	STDMETHOD(InterfaceSupportsErrorInfo)(REFIID riid);	
	STDMETHOD(get_Value)(/*[out, retval]*/ VARIANT* pVal);
	STDMETHOD(put_Value)(/*[in]*/ VARIANT newVal);		
	STDMETHOD(get_SerialNumber)(/*[out, retval]*/ long* pVal);
	STDMETHOD(put_SerialNumber)(/*[in]*/ long newVal);
	STDMETHOD(get_DayCountConvention)(/*[out, retval]*/ DayCountConventionEnum* pVal);
	STDMETHOD(put_DayCountConvention)(/*[in]*/ DayCountConventionEnum newVal);
	STDMETHOD(CopyToClipboard)();
	STDMETHOD(get_DayOfWeek)(/*[out, retval]*/ WeekdayEnum *pVal);
	STDMETHOD(get_DayOfMonth)(/*[out, retval]*/ long* pVal);
	STDMETHOD(get_Year)(/*[out, retval]*/ long* pVal);
	STDMETHOD(get_Month)(/*[out, retval]*/ MonthEnum *pVal);
	STDMETHOD(Add)(/*[in, defaultvalue(0)]*/ long Days, /*[in, defaultvalue(0)]*/ long Months, /*[in, defaultvalue(0)]*/ long Years, /*[in, defaultvalue(0)]*/ BusinessDayConventionEnum BusinessDayConvention, /*[in, optional]*/ VARIANT CalendarOpt, /*[in, optional]*/ VARIANT_BOOL PreserveEoM, /*[out, retval]*/ IDate** pVal);
	STDMETHOD(get_DaysInYear)(/*[out, retval]*/ long* pVal);
	STDMETHOD(AddTenor)(/*[in]*/ BSTR Tenor, /*[in, defaultvalue(0)]*/ BusinessDayConventionEnum BusinessDayConvention, /*[in, optional]*/ VARIANT CalendarOpt, /*[in, optional]*/ VARIANT_BOOL PreserveEoM, /*[out, retval]*/ IDate** pVal);
	STDMETHOD(GetYearFraction)(/*[in]*/ long Date, /*[out, retval]*/ double* pVal);
};

#endif
