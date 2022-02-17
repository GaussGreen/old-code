//	product_future.h : Declaration of the CFuture
//
/////////////////////////////////////////////////////////////////////////////

#ifndef __FUTURE_H_
#define __FUTURE_H_

#include "resource.h"
#include "mleqobjects.h"

class ATL_NO_VTABLE CFuture : 
	public CComObjectRootEx<CComSingleThreadModel>,
	public CComCoClass<CFuture, &CLSID_Future>,
	public ISupportErrorInfo,
	public IDispatchImpl<IFuture, &IID_IFuture, &LIBID_Products>
{
public:
	HRESULT								FinalConstruct(void);
	CFuture(){}
	DECLARE_REGISTRY_RESOURCEID(IDR_FUTURE)
	DECLARE_NOT_AGGREGATABLE(CFuture)
	DECLARE_PROTECT_FINAL_CONSTRUCT()
	BEGIN_COM_MAP(CFuture)
		COM_INTERFACE_ENTRY(IFuture)
		COM_INTERFACE_ENTRY(IDispatch)
		COM_INTERFACE_ENTRY(ISupportErrorInfo)
	END_COM_MAP()	

	// member variables with COM wrappers
	DECLARE_MEMBER_VARIABLE(double, m_fStrike, Strike);
	DECLARE_MEMBER_BOOL(m_bTotalReturn, TotalReturn);
	DECLARE_MEMBER_OBJECT(Asset, m_hUnderlying, Underlying);
	DECLARE_MEMBER_VARIABLE(DATE, m_dateMaturity, Maturity);
	
public:
	STDMETHOD(Evaluate)(IResult* pVal);
	STDMETHOD(InterfaceSupportsErrorInfo)(REFIID riid);
};

#endif
