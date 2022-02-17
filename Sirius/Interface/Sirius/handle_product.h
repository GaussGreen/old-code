//	handle_product.h : Declaration of the CProduct
//
/////////////////////////////////////////////////////////////////////////////

#ifndef __PRODUCT_H_
#define __PRODUCT_H_

#include "resource.h"
#include "MlEqProduct.h"

class ATL_NO_VTABLE CProduct : 
	public CComObjectRootEx<CComSingleThreadModel>,
	public CComCoClass<CProduct, &CLSID_Product>,
	public ISupportErrorInfo,
	public IDispatchImpl<IProduct, &IID_IProduct, &LIBID_Sirius>,
	public IDispatchImpl<IEvaluatable, &IID_IEvaluatable, &LIBID_Sirius>
{
public:	
	CProduct() : m_h(new MlEqProduct) {}	
	DECLARE_REGISTRY_RESOURCEID(IDR_PRODUCT)
	DECLARE_NOT_AGGREGATABLE(CProduct)
	DECLARE_PROTECT_FINAL_CONSTRUCT()
	BEGIN_COM_MAP(CProduct)
		COM_INTERFACE_ENTRY(IProduct)
		COM_INTERFACE_ENTRY(ISupportErrorInfo)
		COM_INTERFACE_ENTRY2(IDispatch, IProduct)
		COM_INTERFACE_ENTRY(IEvaluatable)
	END_COM_MAP()
	associate_analytic_object(CProduct, MlEqProduct);
		
	STDMETHOD(CopyToClipboard)();
	STDMETHOD(Evaluate)(/*[in, defaultvalue("Price")]*/ BSTR Calculate, /*[out, retval]*/ IResult** pVal);
	STDMETHOD(GetFxUnderlyings)(/*[out, retval]*/ IAssets** pVal);
	STDMETHOD(GetObject)(/*[out, retval]*/ IDispatch** pVal);
	STDMETHOD(get_Parameters)(/*[out, retval]*/ IParameters** pVal);
	STDMETHOD(GetRequiredParameters)(/*[out, retval]*/ VARIANT* pVal);
	STDMETHOD(GetUnderlyings)(/*[out, retval]*/ IAssets** pVal);
	STDMETHOD(get_Value)(/*[out, retval]*/ VARIANT *pVal);	
	STDMETHOD(InterfaceSupportsErrorInfo)(REFIID riid);	
	STDMETHOD(put_Parameters)(/*[in]*/ IParameters* newVal);
	STDMETHOD(put_Value)(/*[in]*/ VARIANT newVal);
	STDMETHOD(PutObject)(/*[in]*/ IDispatch* Product);
	STDMETHOD(HasParameter)(/*[in]*/ BSTR Name, /*[out, retval]*/ VARIANT_BOOL* pVal);
	STDMETHOD(GetPayCurrency)(/*[out, retval]*/ CurrencyEnum* pVal);
	STDMETHOD(GetVolatilityStructures)(/*[out, retval]*/ IVolatilityStructures** pVal);
	STDMETHOD(GetZeroCurves)(/*[out, retval]*/ IZeroCurves** pVal);
	STDMETHOD(PutDate)(/*[in]*/ DATE newVal);
	STDMETHOD(PutDataSource)(/*[in]*/ DataSourceEnum newVal);

	HRESULT								GetValue(const std::string& szName, VARIANT *pVal);
	HRESULT								PutValue(const std::string& szName, const CComVariant& Value);
	
protected:
	void								GetParameterMap(CParameterMap* ppm) const;

	declare_member_variable(DATE, Date);
	declare_member_variable(DataSourceEnum, DataSource);
	declare_member_string(ProductType);
};

#endif
