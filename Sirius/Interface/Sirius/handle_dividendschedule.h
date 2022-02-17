//	handle_dividendschedule.h : Declaration of the CDividendSchedule
//
/////////////////////////////////////////////////////////////////////////////

#ifndef __DIVIDENDSCHEDULE_H_
#define __DIVIDENDSCHEDULE_H_

#include "resource.h"
#include "mleqdividendschedule.h"		// for m_h

class ATL_NO_VTABLE CDividendSchedule : 
	public CComObjectRootEx<CComSingleThreadModel>,
	public CComCoClass<CDividendSchedule, &CLSID_DividendSchedule>,
	public ISupportErrorInfo,
	public IDispatchImpl<IDividendSchedule, &IID_IDividendSchedule, &LIBID_Sirius>
{
public:				
	CDividendSchedule() : m_h(new MlEqDividendSchedule) {}
	DECLARE_REGISTRY_RESOURCEID(IDR_DIVIDENDSCHEDULE)
	DECLARE_NOT_AGGREGATABLE(CDividendSchedule)
	DECLARE_PROTECT_FINAL_CONSTRUCT()
	BEGIN_COM_MAP(CDividendSchedule)
		COM_INTERFACE_ENTRY(IDividendSchedule)
		COM_INTERFACE_ENTRY(IDispatch)
		COM_INTERFACE_ENTRY(ISupportErrorInfo)
	END_COM_MAP()
	associate_analytic_object(CDividendSchedule, MlEqDividendSchedule);
	declare_serialisable;

	STDMETHOD(InterfaceSupportsErrorInfo)(REFIID riid);
	STDMETHOD(GetYield)(/*[in]*/ DATE From, /*[in]*/ DATE To, /*[out, retval]*/ double* Yield);
	STDMETHOD(Reset)(void);
	STDMETHOD(Shift)(/*[in]*/ double Yield);
	STDMETHOD(get_Value)(/*[out, retval]*/ VARIANT* pVal);
	STDMETHOD(get_ValueAt)(/*[in]*/ DividendTypeEnum DividendType, /*[in]*/ DATE Date, /*[out, retval]*/ double* pVal);
	STDMETHOD(put_ValueAt)(/*[in]*/ DividendTypeEnum DividendType, /*[in]*/ DATE Date, /*[in]*/ double newVal);
	STDMETHOD(put_Value)(/*[in]*/ VARIANT newVal);
	STDMETHOD(GetPresentValue)(/*[in]*/ IZeroCurve* ZeroCurve, /*[in]*/ DATE From, /*[in]*/ DATE To, /*[out, retval]*/ double* pVal);
	STDMETHOD(Save)(/*[out, retval]*/ BSTR* pVal);
	STDMETHOD(get_Schedule)(/*[out, retval]*/ VARIANT* pVal);
	STDMETHOD(put_Schedule)(/*[in]*/ VARIANT newVal);
	STDMETHOD(get_FirstDate)(/*[in]*/ DividendTypeEnum DividendType, /*[out, retval]*/ DATE* pVal);
	STDMETHOD(get_LastDate)(/*[in]*/ DividendTypeEnum DividendType, /*[out, retval]*/ DATE* pVal);
	STDMETHOD(get_DayCountConvention)(/*[out, retval]*/ DayCountConventionEnum* pVal);
	STDMETHOD(put_DayCountConvention)(/*[in]*/ DayCountConventionEnum newVal);
	STDMETHOD(get_YieldCurve)(/*[out, retval]*/ IZeroCurve** pVal);
	STDMETHOD(put_YieldCurve)(/*[in]*/ IZeroCurve* newVal);
	STDMETHOD(get_Name)(/*[out, retval]*/ BSTR* pVal);
	STDMETHOD(Refresh)(/*[in, defaultvalue(0L)]*/ VARIANT_BOOL Recursive);
	STDMETHOD(DiscreteShift)(/*[in]*/ double Growth, /*[in]*/ BSTR Tenor);
	STDMETHOD(CalibrateYield)(/*[in]*/ double Spot, /*[in]*/ IZeroCurve* ZeroCurve, /*[in]*/ BSTR Tenor, /*[in]*/ BSTR AddAt);
	STDMETHOD(Stick)(void);

	static HRESULT						Load(const std::string& szIdentifier, DataSourceEnum ds, long nDate, CComPtr<IDividendSchedule>& spDividendSchedule);

protected:
	CComPtr<IZeroCurve>					m_spZeroCurve;

	declare_member_variable_rekey(DATE, Date);
	declare_member_variable_rekey(DataSourceEnum, DataSource);
};

#endif

