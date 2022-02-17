//	handle_products.h : Declaration of the CProducts
//
/////////////////////////////////////////////////////////////////////////////

#ifndef __PRODUCTS_H_
#define __PRODUCTS_H_

#include "resource.h"
#include "siriuscomobjectcollection.h"

class ATL_NO_VTABLE CProducts : 
	public CSiriusComObjectCollection<IProduct, CProducts, IProducts, &CLSID_Products, &IID_IProducts, &LIBID_Sirius, &CLSID_Product, IDR_PRODUCTS>,
	public IDispatchImpl<IEvaluatable, &IID_IEvaluatable, &LIBID_Sirius>
{
public:
	BEGIN_COM_MAP(CProducts)
		COM_INTERFACE_ENTRY(IProducts)
		COM_INTERFACE_ENTRY(ISupportErrorInfo)
		COM_INTERFACE_ENTRY2(IDispatch, IProducts)
		COM_INTERFACE_ENTRY(IEvaluatable)
	END_COM_MAP()

	STDMETHOD(Evaluate)(/*[in, defaultvalue("Price")]*/ BSTR Calculate, /*[out, retval]*/ IResult** pVal);
	STDMETHOD(GetFxUnderlyings)(/*[out, retval]*/ IAssets** pVal);
	STDMETHOD(GetUnderlyings)(/*[out, retval]*/ IAssets** pVal);
	STDMETHOD(get_Date)(/*[out, retval]*/ DATE* pVal);
	STDMETHOD(GetVolatilityStructures)(/*[out, retval]*/ IVolatilityStructures** pVal);
	STDMETHOD(GetZeroCurves)(/*[out, retval]*/ IZeroCurves** pVal);
	STDMETHOD(PutDataSource)(/*[in]*/ DataSourceEnum newVal);
	STDMETHOD(PutDate)(/*[in]*/ DATE newVal);	
};

#endif
