//	product_mattoption.h : Declaration of the CMattOption.
//						   This is a sum of two calls designed to test
//						   Multi-asset all input.
//
/////////////////////////////////////////////////////////////////////////////

#ifndef __MATTOPTION_H_
#define __MATTOPTION_H_

#include "resource.h"
#include "mleqobjects.h"

class ATL_NO_VTABLE CMattOption : 
	public CComObjectRootEx<CComSingleThreadModel>,
	public CComCoClass<CMattOption, &CLSID_MattOption>,
	public ISupportErrorInfo,
	public IDispatchImpl<IMattOption, &IID_IMattOption, &LIBID_TemporaryProductsLib>
{
public:
	CMattOption(){}
	DECLARE_REGISTRY_RESOURCEID(IDR_MATTOPTION)
	DECLARE_PROTECT_FINAL_CONSTRUCT()
	BEGIN_COM_MAP(CMattOption)
		COM_INTERFACE_ENTRY(IMattOption)
		COM_INTERFACE_ENTRY(IDispatch)
		COM_INTERFACE_ENTRY(ISupportErrorInfo)
	END_COM_MAP()
	HRESULT								FinalConstruct(void);

	DECLARE_MEMBER_OBJECT(Asset, m_hUnderlying1, Underlying1)
	DECLARE_MEMBER_OBJECT(Asset, m_hUnderlying2, Underlying2)
	DECLARE_MEMBER_VARIABLE(DATE, m_datePay, PayDate)
	DECLARE_MEMBER_VARIABLE(double, m_fStrike1, Strike1)
	DECLARE_MEMBER_VARIABLE(double, m_fStrike2, Strike2)
	DECLARE_MEMBER_VARIABLE(double, m_fNotional1, Notional1)
	DECLARE_MEMBER_VARIABLE(double, m_fNotional2, Notional2)
	DECLARE_MEMBER_VARIABLE(DATE, m_dateMaturity, Maturity)

public:
	STDMETHOD(InterfaceSupportsErrorInfo)(REFIID riid);
	STDMETHOD(Evaluate)(IResult* pVal);
};

#endif