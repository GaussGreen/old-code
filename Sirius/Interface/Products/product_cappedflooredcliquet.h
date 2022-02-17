//	product_cappedflooredcliquet.h : Declaration of the CCappedFlooredCliquet
//
/////////////////////////////////////////////////////////////////////////////

#ifndef __CAPPEDFLOOREDCLIQUET_H_
#define __CAPPEDFLOOREDCLIQUET_H_

#include "resource.h"
#include "cappedflooredcliquetpricer.h"
#include "MlEqObjects.h"

class ATL_NO_VTABLE CCappedFlooredCliquet : 
	public CComObjectRootEx<CComSingleThreadModel>,
	public CComCoClass<CCappedFlooredCliquet, &CLSID_CappedFlooredCliquet>,	
	public ISupportErrorInfo,
	public IDispatchImpl<ICappedFlooredCliquet, &IID_ICappedFlooredCliquet, &LIBID_Products>,
	public CappedFlooredCliquetPricer
{
public:
	HRESULT								FinalConstruct(void);
	CCappedFlooredCliquet(){}
	DECLARE_REGISTRY_RESOURCEID(IDR_CAPPEDFLOOREDCLIQUET)
	DECLARE_PROTECT_FINAL_CONSTRUCT()
	BEGIN_COM_MAP(CCappedFlooredCliquet)
		COM_INTERFACE_ENTRY(ICappedFlooredCliquet)
		COM_INTERFACE_ENTRY(IDispatch)
		COM_INTERFACE_ENTRY(ISupportErrorInfo)
	END_COM_MAP()
	
	DECLARE_MEMBER_VARIABLE(double, m_fGlobalCap, GlobalCap)
	DECLARE_MEMBER_VARIABLE(double, m_fGlobalFloor, GlobalFloor)
	DECLARE_MEMBER_VARIABLE(double, m_fGlobalRedemption, GlobalRedemption)
	DECLARE_MEMBER_VARIABLE(double, m_fGlobalGearing, GlobalGearing)
	DECLARE_MEMBER_VARIABLE(double, m_fNotional, Notional)
	DECLARE_MEMBER_OBJECT(Array, m_hCallPuts, CallPutArray);			
	DECLARE_MEMBER_OBJECT(Array, m_hLocalCaps, LocalCaps);
	DECLARE_MEMBER_OBJECT(Array, m_hLocalFloors, LocalFloors);
	DECLARE_MEMBER_OBJECT(Array, m_hStrikes, Strikes);	
	DECLARE_MEMBER_OBJECT(Array, m_hWeights, Weights);
	DECLARE_MEMBER_OBJECT(ParameterList, m_hParameterList, ModelInfo);
	DECLARE_MEMBER_OBJECT(MonteCarlo, m_hMonteCarlo, ForwardSkewHandle);
	DECLARE_MEMBER_OBJECT(Asset, m_hUnderlying, Underlying);			
	DECLARE_MEMBER_OBJECT(DateSchedule, m_hFixingSchedule, DateAndFixingSchedule);
	DECLARE_MEMBER_OBJECT(Matrix, m_hLockinSchedule, LockinScheduleOpt);
	DECLARE_MEMBER_BOOL(m_bRebalancingBasket, RebalancingBasketOpt);

public:		
	STDMETHOD(InterfaceSupportsErrorInfo)(REFIID riid);
	STDMETHOD(Evaluate)(IResult* pVal);
};

#endif
