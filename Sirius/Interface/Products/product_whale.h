// product_whale.h : Declaration of the CWhale

//#pragma once

#ifndef __WHALE_H_
#define __WHALE_H_

#include "resource.h"      
#include "MlEqObjects.h"
#include "WhalePricer.h"

/////////////////////////////////////////////////////////////////////////////
// CWhale
class ATL_NO_VTABLE CWhale : 
	public CComObjectRootEx<CComSingleThreadModel>,
	public CComCoClass<CWhale, &CLSID_Whale>,
	public ISupportErrorInfo,
	public IDispatchImpl<IWhale, &IID_IWhale, &LIBID_Products>,
	public WhalePricer
{
public:
	HRESULT		FinalConstruct(void);
	CWhale(){};
	DECLARE_REGISTRY_RESOURCEID(IDR_WHALE)
	DECLARE_PROTECT_FINAL_CONSTRUCT()
	BEGIN_COM_MAP(CWhale)
		COM_INTERFACE_ENTRY(IWhale)
		COM_INTERFACE_ENTRY(IDispatch)		
		COM_INTERFACE_ENTRY(ISupportErrorInfo)
	END_COM_MAP()

	DECLARE_MEMBER_VARIABLE(double, m_fStrike, Strike);
	DECLARE_MEMBER_VARIABLE(DATE, m_datePay, PayDate);
	DECLARE_MEMBER_VARIABLE(PayoffTypeEnum, m_cpe, CallOrPut);
	DECLARE_MEMBER_OBJECT(Array, m_hWeights, Weights);
	DECLARE_MEMBER_OBJECT(ParameterList, m_hParameterList, ModelInfo);
	DECLARE_MEMBER_OBJECT(MonteCarlo, m_hMonteCarlo, ForwardSkewHandle);
	DECLARE_MEMBER_OBJECT(Asset, m_hUnderlying, Underlying);			
	DECLARE_MEMBER_OBJECT(DateSchedule, m_hFixingSchedule, DateAndFixingSchedule);

public:	
	STDMETHOD(Evaluate)(IResult* pVal);
	STDMETHOD(InterfaceSupportsErrorInfo)(REFIID riid);
};

#endif //__WHALE_H_
