//	handle_deal.h : Declaration of the CDeal
//
/////////////////////////////////////////////////////////////////////////////

#ifndef __DEAL_H_
#define __DEAL_H_

#include "resource.h"
#include "MlEqDeal.h"

class ATL_NO_VTABLE CDeal : 
	public CComObjectRootEx<CComSingleThreadModel>,
	public CComCoClass<CDeal, &CLSID_Deal>,
	public ISupportErrorInfo,
	public IDispatchImpl<IDeal, &IID_IDeal, &LIBID_Sirius>,
	public IDispatchImpl<IEvaluatable, &IID_IEvaluatable, &LIBID_Sirius>
{
public:				
	CDeal() : m_h(new MlEqDeal){}	
	DECLARE_REGISTRY_RESOURCEID(IDR_DEAL)
	DECLARE_NOT_AGGREGATABLE(CDeal)
	DECLARE_PROTECT_FINAL_CONSTRUCT()
	BEGIN_COM_MAP(CDeal)
		COM_INTERFACE_ENTRY(IDeal)
		COM_INTERFACE_ENTRY(ISupportErrorInfo)
		COM_INTERFACE_ENTRY2(IDispatch, IDeal)
		COM_INTERFACE_ENTRY(IEvaluatable)		
	END_COM_MAP()
	associate_analytic_object(CDeal, MlEqDeal);
	declare_serialisable;
	static HRESULT						Load(const std::string& szName, DataSourceEnum ds, long nDate, CComPtr<IDeal>& spDeal);

	STDMETHOD(Evaluate)(/*[in, defaultvalue("Price")]*/ BSTR Calculate, /*[out, retval]*/ IResult** pVal);	
	STDMETHOD(GetFxUnderlyings)(/*[out, retval]*/ IAssets** pVal);
	STDMETHOD(get_PositionNotional)(/*[in]*/ IPosition* Position, /*[out, retval]*/ double *pVal);
	STDMETHOD(get_Positions)(/*[out, retval]*/ IPositions** pVal);
	STDMETHOD(GetUnderlyings)(/*[out, retval]*/ IAssets** pVal);
	STDMETHOD(get_Value)(/*[out, retval]*/ VARIANT *pVal);	
	STDMETHOD(InterfaceSupportsErrorInfo)(REFIID riid);
	STDMETHOD(put_PositionNotional)(/*[in]*/ IPosition* Position, /*[in]*/ double newVal);
	STDMETHOD(put_Positions)(/*[in]*/ IPositions* newVal);
	STDMETHOD(put_Value)(/*[in]*/ VARIANT newVal);	
	STDMETHOD(Save)(/*[out, retval]*/ BSTR* pVal);
	STDMETHOD(GetVolatilityStructures)(/*[out, retval]*/ IVolatilityStructures** pVal);
	STDMETHOD(GetZeroCurves)(/*[out, retval]*/ IZeroCurves** pVal);
	STDMETHOD(PutDataSource)(/*[in]*/ DataSourceEnum newVal);
	STDMETHOD(PutDate)(/*[in]*/ DATE newVal);
	
protected:	
	static const long					s_nMaxPositionsInDeal;
	declare_member_variable_rekey(DATE, Date);
	declare_member_variable_rekey(DataSourceEnum, DataSource);	
	declare_member_string_rekey(Name);
	declare_member_variable(double, Notional);
	declare_member_string(Book);
};

#endif
