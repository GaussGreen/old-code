// product_funding.h : Declaration of the CFunding

#ifndef __FUNDING_H_
#define __FUNDING_H_

#include "resource.h"       // main symbols
#include "MlEqObjects.h"



class ATL_NO_VTABLE CFunding : 
	public CComObjectRootEx<CComSingleThreadModel>,
	public CComCoClass<CFunding, &CLSID_Funding>,
	public ISupportErrorInfo,
	public IDispatchImpl<IFunding, &IID_IFunding, &LIBID_TemporaryProductsLib>
{
public:
	CFunding(){}
	HRESULT FinalConstruct();
	DECLARE_REGISTRY_RESOURCEID(IDR_FUNDING)
	DECLARE_NOT_AGGREGATABLE(CFunding)
	DECLARE_PROTECT_FINAL_CONSTRUCT()
	BEGIN_COM_MAP(CFunding)
		COM_INTERFACE_ENTRY(IFunding)
		COM_INTERFACE_ENTRY(IDispatch)
		COM_INTERFACE_ENTRY(ISupportErrorInfo)
	END_COM_MAP()

	STDMETHOD(InterfaceSupportsErrorInfo)(REFIID riid);
	DECLARE_MEMBER_OBJECT(ZeroCurve, m_hZeroCurve, ZeroCurve);
	DECLARE_MEMBER_OBJECT(Array, m_hCoupons, Coupons);			
	DECLARE_MEMBER_OBJECT(DateSchedule, m_hFixingSchedule, DateSchedule);

public:
	STDMETHOD(Evaluate)(IResult* pVal);
};

#endif
