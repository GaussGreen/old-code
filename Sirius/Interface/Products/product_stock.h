// product_stock.h : Declaration of the CStock

#ifndef __STOCK_H_
#define __STOCK_H_

#include "resource.h"
#include "mleqobjects.h"

/////////////////////////////////////////////////////////////////////////////
// CStock
class ATL_NO_VTABLE CStock : 
	public CComObjectRootEx<CComSingleThreadModel>,
	public CComCoClass<CStock, &CLSID_Stock>,
	public ISupportErrorInfo,
	public IDispatchImpl<IStock, &IID_IStock, &LIBID_Products>
{
public:
	CStock()
	{
	}

DECLARE_REGISTRY_RESOURCEID(IDR_STOCK)
DECLARE_NOT_AGGREGATABLE(CStock)

DECLARE_PROTECT_FINAL_CONSTRUCT()

BEGIN_COM_MAP(CStock)
	COM_INTERFACE_ENTRY(IStock)
	COM_INTERFACE_ENTRY(IDispatch)
	COM_INTERFACE_ENTRY(ISupportErrorInfo)
END_COM_MAP()

// ISupportsErrorInfo
	STDMETHOD(InterfaceSupportsErrorInfo)(REFIID riid);


public:
	DECLARE_MEMBER_OBJECT(Asset, m_hUnderlying, Underlying)
	STDMETHOD(Evaluate)(IResult* pVal);
};

#endif
