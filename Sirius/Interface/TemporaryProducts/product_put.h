//	product_put.h : Declaration of the CProductPut
//
/////////////////////////////////////////////////////////////////////////////

#ifndef __PRODUCTPUT_H_
#define __PRODUCTPUT_H_

#include "resource.h"
#include "mleqobjects.h"

class ATL_NO_VTABLE CProductPut : 
	public CComObjectRootEx<CComSingleThreadModel>,
	public CComCoClass<CProductPut, &CLSID_ProductPut>,
	public ISupportErrorInfo,
	public IDispatchImpl<IProductPut, &IID_IProductPut, &LIBID_TemporaryProductsLib>
{
public:
	HRESULT FinalConstruct(void);
	DECLARE_REGISTRY_RESOURCEID(IDR_PRODUCTPUT)
	DECLARE_PROTECT_FINAL_CONSTRUCT()
	BEGIN_COM_MAP(CProductPut)
		COM_INTERFACE_ENTRY(IProductPut)
		COM_INTERFACE_ENTRY(IDispatch)
		COM_INTERFACE_ENTRY(ISupportErrorInfo)
	END_COM_MAP()
	CProductPut(){}

	// member variables with COM wrappers	
	DECLARE_MEMBER_VARIABLE(double, m_fStrike, Strike)
	DECLARE_MEMBER_VARIABLE(DATE, m_dateMaturity, Maturity)
	DECLARE_MEMBER_OBJECT(Asset, m_hUnderlying, Underlying)

public:
	STDMETHOD(Evaluate)(IResult* pVal);
	STDMETHOD(InterfaceSupportsErrorInfo)(REFIID riid);
};

#endif
