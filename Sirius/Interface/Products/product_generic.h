//	product_generic.h : Declaration of the CGeneric
//						This is a product that simply stores its value
//						and various hedge parameters. It's analogue is the
//                      RAM all-input.
//
//	author :            David Cuin
//
/////////////////////////////////////////////////////////////////////////////

#ifndef __GENERIC_H_
#define __GENERIC_H_

#include "resource.h"

class ATL_NO_VTABLE CGeneric : 
	public CComObjectRootEx<CComSingleThreadModel>,
	public CComCoClass<CGeneric, &CLSID_Generic>,
	public ISupportErrorInfo,
	public IDispatchImpl<IGeneric, &IID_IGeneric, &LIBID_Products>
{
public:
	HRESULT FinalConstruct(void);		
	CGeneric(){}
	DECLARE_REGISTRY_RESOURCEID(IDR_GENERIC)
	DECLARE_NOT_AGGREGATABLE(CGeneric)
	DECLARE_PROTECT_FINAL_CONSTRUCT()
	BEGIN_COM_MAP(CGeneric)
		COM_INTERFACE_ENTRY(IGeneric)
		COM_INTERFACE_ENTRY(IDispatch)
		COM_INTERFACE_ENTRY(ISupportErrorInfo)
	END_COM_MAP()
	
	// member variables with COM wrappers
	DECLARE_MEMBER_VARIABLE(double, m_fPrice, Price)
	DECLARE_MEMBER_VARIABLE(double, m_fDelta, Delta)
	DECLARE_MEMBER_VARIABLE(double, m_fGamma, Gamma)
	DECLARE_MEMBER_VARIABLE(double, m_fVega, Vega)
	DECLARE_MEMBER_VARIABLE(double, m_fFXDelta, FXDelta)
	DECLARE_MEMBER_STRING(m_szRamModel, RamModel)
	DECLARE_MEMBER_VARIABLE(CurrencyEnum, m_cePay, PayCurrency)	
	DECLARE_MEMBER_STRING(m_szUnderlying, Underlying)

public:	
	STDMETHOD(Evaluate)(IResult* pVal);
	STDMETHOD(InterfaceSupportsErrorInfo)(REFIID riid);
};

#endif
