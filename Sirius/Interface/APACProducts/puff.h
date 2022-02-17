//	Puff.h : Declaration of the CPuff
//
/////////////////////////////////////////////////////////////////////////////

#ifndef __PUFF_H_
#define __PUFF_H_

#include "resource.h"
#include "MlEqPuff.h"

class MlEqPuff;
typedef RCPtr<MlEqPuff>				MlEqPuffHandle;


class ATL_NO_VTABLE CPuff : 
	public CComObjectRootEx<CComSingleThreadModel>,
	public CComCoClass<CPuff, &CLSID_Puff>,
	public ISupportErrorInfo,
	public IDispatchImpl<IPuff, &IID_IPuff, &LIBID_APACProducts>
{
public:
	CPuff() : m_h(new MlEqPuff) {}

	DECLARE_REGISTRY_RESOURCEID(IDR_PUFF)
	DECLARE_NOT_AGGREGATABLE(CPuff)
	DECLARE_PROTECT_FINAL_CONSTRUCT()
	BEGIN_COM_MAP(CPuff)
		COM_INTERFACE_ENTRY(IPuff)
		COM_INTERFACE_ENTRY(IDispatch)
		COM_INTERFACE_ENTRY(ISupportErrorInfo)
	END_COM_MAP()
	STDMETHOD(InterfaceSupportsErrorInfo)(REFIID riid);
	associate_analytic_object(CPuff, MlEqPuff);

	STDMETHOD(Evaluate)(IResult* pVal);
	
	declare_member_variable(double, Strike);
	declare_member_variable(double, Maturity);
	declare_member_variable(double, Multiplier);
	declare_member_object(Asset, Underlying);
};

#endif
