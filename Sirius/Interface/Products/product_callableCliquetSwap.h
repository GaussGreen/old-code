// product_callableCliquetSwap.h : Declaration of the CCallableCliquetSwap

#pragma once

#include "resource.h"       // main symbols
#include "CallableCliquetSwapPricer.h"
#include "MlEqObjects.h"


/////////////////////////////////////////////////////////////////////////////
// CCallableCliquetSwap
class ATL_NO_VTABLE CCallableCliquetSwap : 
	public CComObjectRootEx<CComSingleThreadModel>,
	public CComCoClass<CCallableCliquetSwap, &CLSID_CallableCliquetSwap>,
	public ISupportErrorInfo,
	public IDispatchImpl<ICallableCliquetSwap, &IID_ICallableCliquetSwap, &LIBID_Products>,
	public CallableCliquetSwapPricer
{
public:
	CCallableCliquetSwap(){};
	HRESULT		FinalConstruct(void);


	DECLARE_REGISTRY_RESOURCEID(IDR_CALLABLECLIQUETSWAP)
	DECLARE_NOT_AGGREGATABLE(CCallableCliquetSwap)

	DECLARE_PROTECT_FINAL_CONSTRUCT()

	BEGIN_COM_MAP(CCallableCliquetSwap)
		COM_INTERFACE_ENTRY(ICallableCliquetSwap)
		COM_INTERFACE_ENTRY(IDispatch)
		COM_INTERFACE_ENTRY(ISupportErrorInfo)
	END_COM_MAP()

public:
	DECLARE_MEMBER_OBJECT(DateSchedule, m_hCouponSchedule, CouponSchedule);
	DECLARE_MEMBER_OBJECT(Array, m_hGlobalCap, GlobalCap);
	DECLARE_MEMBER_OBJECT(Array, m_hGlobalFloor, GlobalFloor);
	DECLARE_MEMBER_OBJECT(Array, m_hFixedCoupon, FixedCoupon);

	DECLARE_MEMBER_VARIABLE(double, m_fCallableLevel, CallableLevel);

	DECLARE_MEMBER_OBJECT(Array, m_hCallPuts, CallPutArray);			
	DECLARE_MEMBER_OBJECT(Array, m_hLocalCaps, LocalCaps);
	DECLARE_MEMBER_OBJECT(Array, m_hLocalFloors, LocalFloors);
	DECLARE_MEMBER_OBJECT(Array, m_hStrikes, Strikes);	
	DECLARE_MEMBER_OBJECT(Array, m_hWeights, Weights);
	DECLARE_MEMBER_OBJECT(ParameterList, m_hParameterList, ModelInfo);
	DECLARE_MEMBER_OBJECT(MonteCarlo, m_hMonteCarlo, ForwardSkewHandle);
	DECLARE_MEMBER_OBJECT(Asset, m_hUnderlying, Underlying);			
	DECLARE_MEMBER_OBJECT(DateSchedule, m_hFixingSchedule, DateAndFixingSchedule);
	DECLARE_MEMBER_OBJECT(Array, m_hPeriodIdentifier, PeriodIdentifier);

//	DECLARE_MEMBER_OBJECT(ZeroCurve, m_hZeroCurve, ZeroCurve);
//	DECLARE_MEMBER_OBJECT(DateSchedule, m_hFundingSchedule, FundingSchedule);
//	DECLARE_MEMBER_VARIABLE(double, m_fFunding, Funding);

//	DECLARE_MEMBER_OBJECT(Swap, m_hSwap, SwapComponent);
	
public:		
	STDMETHOD(InterfaceSupportsErrorInfo)(REFIID riid);
	STDMETHOD(Evaluate)(IResult* pVal);
};

