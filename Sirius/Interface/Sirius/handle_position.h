//	handle_position.h : Declaration of the CPosition
//
/////////////////////////////////////////////////////////////////////////////

#ifndef __POSITION_H_
#define __POSITION_H_

#include "resource.h"
#include "MlEqPosition.h"

class ATL_NO_VTABLE CPosition : 
	public CComObjectRootEx<CComSingleThreadModel>,
	public CComCoClass<CPosition, &CLSID_Position>,
	public ISupportErrorInfo,
	public IDispatchImpl<IPosition, &IID_IPosition, &LIBID_Sirius>,
	public IDispatchImpl<IEvaluatable, &IID_IEvaluatable, &LIBID_Sirius>
{
public:	
	CPosition() : m_h(new MlEqPosition){}	
	DECLARE_REGISTRY_RESOURCEID(IDR_POSITION)
	DECLARE_NOT_AGGREGATABLE(CPosition)
	DECLARE_PROTECT_FINAL_CONSTRUCT()
	BEGIN_COM_MAP(CPosition)
		COM_INTERFACE_ENTRY(IPosition)
		COM_INTERFACE_ENTRY(ISupportErrorInfo)
		COM_INTERFACE_ENTRY2(IDispatch, IPosition)
		COM_INTERFACE_ENTRY(IEvaluatable)		
	END_COM_MAP()			
	associate_analytic_object(CPosition, MlEqPosition);
	declare_serialisable;
	static HRESULT						Load(const std::string& szName, DataSourceEnum ds, long nDate, CComPtr<IPosition>& spPosition);
	
	STDMETHOD(InterfaceSupportsErrorInfo)(REFIID riid);
	STDMETHOD(Evaluate)(/*[in, defaultvalue("Price")]*/ BSTR Calculate, /*[out, retval]*/ IResult** pVal);
	STDMETHOD(get_ProductNotional)(/*[in]*/ IProduct* Product, /*[out, retval]*/ double *pVal);
	STDMETHOD(GetUnderlyings)(/*[out, retval]*/ IAssets** pVal);
	STDMETHOD(get_Value)(/*[out, retval]*/ VARIANT* pVal);
	STDMETHOD(put_ProductNotional)(/*[in]*/ IProduct* Product, /*[in]*/ double newVal);
	STDMETHOD(put_Value)(/*[in]*/ VARIANT newVal);
	STDMETHOD(Save)(/*[out, retval]*/ BSTR* pVal);
	STDMETHOD(GetFxUnderlyings)(/*[out, retval]*/ IAssets** pVal);
	STDMETHOD(get_Products)(/*[out, retval]*/ IProducts** pVal);
	STDMETHOD(put_Products)(/*[in]*/ IProducts* newVal);
	STDMETHOD(get_Result)(/*[out, retval]*/ IResult** pVal);
	STDMETHOD(put_Result)(/*[in]*/ IResult* newVal);
	STDMETHOD(GetVolatilityStructures)(/*[out, retval]*/ IVolatilityStructures** pVal);
	STDMETHOD(GetZeroCurves)(/*[out, retval]*/ IZeroCurves** pVal);
	STDMETHOD(PutDate)(/*[in]*/ DATE newVal);
	STDMETHOD(PutDataSource)(/*[in]*/ DataSourceEnum newVal);

protected:	
	declare_member_variable_rekey(DATE, Date);
	declare_member_variable_rekey(DataSourceEnum, DataSource);
	declare_member_string_rekey(Name);
	declare_member_string(Description);
	declare_member_string(Counterparty);
	declare_member_string(Book);
	declare_member_variable(double, Notional);
	declare_member_string(Comment);
	declare_member_string(Strategy);
};

#endif
