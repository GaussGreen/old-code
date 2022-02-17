//	handle_zerocurve.h : Declaration of CZeroCurve
//
/////////////////////////////////////////////////////////////////////////////

#ifndef __ZEROCURVE_H_
#define __ZEROCURVE_H_

#include "resource.h"
#include "mleqzerocurve.h"

class ATL_NO_VTABLE CZeroCurve :
	public CComObjectRootEx<CComSingleThreadModel>,
	public CComCoClass<CZeroCurve, &CLSID_ZeroCurve>,
	public ISupportErrorInfo,
	public IDispatchImpl<IZeroCurve, &IID_IZeroCurve, &LIBID_Sirius>
{
public:		
	CZeroCurve() : m_h(new MlEqZeroCurve) {}
	DECLARE_REGISTRY_RESOURCEID(IDR_ZEROCURVE)
	DECLARE_NOT_AGGREGATABLE(CZeroCurve)
	DECLARE_PROTECT_FINAL_CONSTRUCT()
	BEGIN_COM_MAP(CZeroCurve)
		COM_INTERFACE_ENTRY(IZeroCurve)
		COM_INTERFACE_ENTRY(IDispatch)
		COM_INTERFACE_ENTRY(ISupportErrorInfo)
	END_COM_MAP()
	associate_analytic_object(CZeroCurve, MlEqZeroCurve);
	declare_serialisable;
	HRESULT								FinalConstruct(void);
	static HRESULT						Load(const std::string& szIdentifier, DataSourceEnum ds, long nDate, CComPtr<IZeroCurve>& spZeroCurve);
	std::string							GetGDAHandle(void) const;
			
	declare_member_variable_rekey(DATE, Date);
	declare_member_variable_rekey(DataSourceEnum, DataSource);
	STDMETHOD(InterfaceSupportsErrorInfo)(REFIID riid);
	STDMETHOD(get_Value)(/*[out, retval]*/ VARIANT *pVal);
	STDMETHOD(put_Value)(/*[in]*/ VARIANT newVal);
	STDMETHOD(get_Name)(/*[out, retval]*/ BSTR* pVal);	
	STDMETHOD(GetDiscountFactor)(/*[in]*/ DATE Maturity, /*[out, retval]*/ double* pVal);
	STDMETHOD(GetForwardDiscountFactor)(/*[in]*/ DATE From, /*[in]*/ DATE To, /*[out, retval]*/ double* pVal);
	STDMETHOD(Reset)(void);
	STDMETHOD(Shift)(/*[in]*/ double Amount);	
	STDMETHOD(get_FxSpot)(/*[out, retval]*/ double *pVal);
	STDMETHOD(put_FxSpot)(/*[in]*/ double newVal);
	STDMETHOD(Save)(/*[out, retval]*/ BSTR* pVal);
	STDMETHOD(get_InterpolatorType)(/*[out, retval]*/ BSTR *pVal);
	STDMETHOD(put_InterpolatorType)(/*[in]*/ BSTR newVal);
	STDMETHOD(GetCashRate)(/*[in]*/ BSTR TermOrDate, /*[out, retval]*/ double* pVal);
	STDMETHOD(Default)(/*[in]*/ BSTR Currency, /*[in]*/ BSTR Request, /*[out, retval]*/ BSTR* pVal);	
	STDMETHOD(ShiftByInterpolator)(/*[in]*/ IInterpolator* Interpolator);
	STDMETHOD(Stick)(void);
};

#endif
