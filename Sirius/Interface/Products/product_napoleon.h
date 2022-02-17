//	product_napoleon.h : Declaration of the CNapoleon
//
/////////////////////////////////////////////////////////////////////////////

#ifndef __NAPOLEON_H_
#define __NAPOLEON_H_

#include "resource.h"
#include "NapoleonPricer.h"
#include "MlEqObjects.h"

class ATL_NO_VTABLE CNapoleon : 
	public CComObjectRootEx<CComSingleThreadModel>,
	public CComCoClass<CNapoleon, &CLSID_Napoleon>,
	public ISupportErrorInfo,
	public IDispatchImpl<INapoleon, &IID_INapoleon, &LIBID_Products>,
	public MlEqNapoleonPricer
{
public:
	HRESULT								FinalConstruct(void);
	CNapoleon(){}
	DECLARE_REGISTRY_RESOURCEID(IDR_NAPOLEON)
	DECLARE_PROTECT_FINAL_CONSTRUCT()
	BEGIN_COM_MAP(CNapoleon)
		COM_INTERFACE_ENTRY(INapoleon)
		COM_INTERFACE_ENTRY(IDispatch)	
		COM_INTERFACE_ENTRY(ISupportErrorInfo)
	END_COM_MAP()

	DECLARE_MEMBER_OBJECT(DateSchedule, m_hPayDateScheduleOpt, PayDateScheduleOpt)
	DECLARE_MEMBER_OBJECT(DateSchedule, m_hFixingSchedule, DateAndFixingSchedule)
	DECLARE_MEMBER_OBJECT(Array, m_hStartsPeriodArr, StartsPeriodArr)
	DECLARE_MEMBER_OBJECT(Array, m_hCallPutArr, CallPutArr)
	DECLARE_MEMBER_OBJECT(Array, m_hStrikesArr, StrikesArr)
	DECLARE_MEMBER_OBJECT(Array, m_hLocalFloorsArr, LocalFloorsArr)
	DECLARE_MEMBER_OBJECT(Array, m_hLocalCapsArr, LocalCapsArr)
	DECLARE_MEMBER_OBJECT(Array, m_hWeightsArr, WeightsArr)
	DECLARE_MEMBER_OBJECT(Array, m_hCouponsArr, CouponsArr)
	DECLARE_MEMBER_OBJECT(Array, m_hGearingsArr, GearingsArr)
	DECLARE_MEMBER_OBJECT(Array, m_hPeriodFloorsArr, PeriodFloorsArr)
	DECLARE_MEMBER_OBJECT(Array, m_hPeriodCapsArr, PeriodCapsArr)
	DECLARE_MEMBER_BOOL(m_bDelayPaymentsToEnd, DelayPaymentsToEnd)
	DECLARE_MEMBER_OBJECT(ParameterList, m_hParameterList, ModelInfo)
	DECLARE_MEMBER_OBJECT(Asset, m_hUnderlying, Underlying)	

public:
	STDMETHOD(Evaluate)(IResult* pVal);
	STDMETHOD(InterfaceSupportsErrorInfo)(REFIID riid);
};

#endif
