//	handle_strikes.h : Declaration of the CStrikes
//
/////////////////////////////////////////////////////////////////////////////

#ifndef __STRIKES_H_
#define __STRIKES_H_

#include "resource.h"
#include "MlEqObjects.h"

class ATL_NO_VTABLE CStrikes : 
	public CComObjectRootEx<CComSingleThreadModel>,
	public CComCoClass<CStrikes, &CLSID_Strikes>,
	public ISupportErrorInfo,
	public IDispatchImpl<IStrikes, &IID_IStrikes, &LIBID_Sirius>
{
public:
	enum StrikesTypeFromAsset // enumerators detailing methods on how to create strikes
							  // that are never passed into analytic classes.
							  // Make sure that their values do not clash with
							  // StrikesTypeEnum
	{
		NormalisedFromAsset = -100,
		ForwardBasedFromAsset,
		SpotBasedFromAsset,
		LogBasedFromAsset
	};
	
public:	
	CStrikes() : m_h(new MlEqStrikes) {}
	DECLARE_REGISTRY_RESOURCEID(IDR_STRIKES)
	DECLARE_NOT_AGGREGATABLE(CStrikes)
	DECLARE_PROTECT_FINAL_CONSTRUCT()
	BEGIN_COM_MAP(CStrikes)
		COM_INTERFACE_ENTRY(IStrikes)
		COM_INTERFACE_ENTRY(IDispatch)
		COM_INTERFACE_ENTRY(ISupportErrorInfo)
	END_COM_MAP()
	associate_analytic_object(CStrikes, MlEqStrikes);
	
	STDMETHOD(InterfaceSupportsErrorInfo)(REFIID riid);
	STDMETHOD(get_Value)(/*[out, retval]*/ VARIANT *pVal);
	STDMETHOD(put_Value)(/*[in]*/ VARIANT newVal);
	STDMETHOD(get_StrikesType)(/*[out, retval]*/ StrikesTypeEnum *pVal);	

protected:
	CComPtr<IDate>						m_spStartDate;
};

#endif
