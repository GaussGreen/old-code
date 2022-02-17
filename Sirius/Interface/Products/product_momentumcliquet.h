// product_momentumcliquet.h : Declaration of the CMomentumCliquet

#ifndef __MOMENTUMCLIQUET_H_
#define __MOMENTUMCLIQUET_H_

#include "resource.h"    
#include "MomentumCliquetPricer.h"
#include "MlEqObjects.h"

/////////////////////////////////////////////////////////////////////////////
// CMomentumCliquet
class ATL_NO_VTABLE CMomentumCliquet : 
	public CComObjectRootEx<CComSingleThreadModel>,
	public CComCoClass<CMomentumCliquet, &CLSID_MomentumCliquet>,
	public ISupportErrorInfo,
	public IDispatchImpl<IMomentumCliquet, &IID_IMomentumCliquet, &LIBID_Products>,
	public GeneralCliquetPricer
{
public:
	HRESULT								FinalConstruct(void);
	CMomentumCliquet(){};

DECLARE_REGISTRY_RESOURCEID(IDR_MOMENTUMCLIQUET)
DECLARE_NOT_AGGREGATABLE(CMomentumCliquet)

DECLARE_PROTECT_FINAL_CONSTRUCT()

BEGIN_COM_MAP(CMomentumCliquet)
	COM_INTERFACE_ENTRY(IMomentumCliquet)
	COM_INTERFACE_ENTRY(IDispatch)
	COM_INTERFACE_ENTRY(ISupportErrorInfo)
END_COM_MAP()

	DECLARE_MEMBER_OBJECT(DateSchedule, m_hPayDateSchedule, PayDateSchedule)
	DECLARE_MEMBER_OBJECT(DateSchedule, m_hFixingSchedule, DateAndFixingSchedule)
	DECLARE_MEMBER_OBJECT(Array, m_hStartsPeriodArr, StartsPeriodArr)
	DECLARE_MEMBER_OBJECT(Array, m_hStrikesArr, StrikesArr)
	DECLARE_MEMBER_OBJECT(Array, m_hLocalFloorsArr, LocalFloorsArr)
	DECLARE_MEMBER_OBJECT(Array, m_hLocalCapsArr, LocalCapsArr)
	DECLARE_MEMBER_OBJECT(Array, m_hWeightsArr, WeightsArr)
	DECLARE_MEMBER_OBJECT(Array, m_hCouponsArr, CouponsArr)
	DECLARE_MEMBER_OBJECT(Array, m_hGearingsArr, GearingsArr)
	DECLARE_MEMBER_OBJECT(Array, m_hPeriodFloorsArr, PeriodFloorsArr)
	DECLARE_MEMBER_OBJECT(Array, m_hPeriodCapsArr, PeriodCapsArr)
	DECLARE_MEMBER_VARIABLE(PayoffTypeEnum, m_PayoffType, PayoffType)
	DECLARE_MEMBER_VARIABLE(MinMaxEnum, m_MinMaxType, ExtremumType)
	DECLARE_MEMBER_OBJECT(ParameterList, m_hParameterList, ModelInfo)
	DECLARE_MEMBER_OBJECT(Asset, m_hUnderlying, Underlying)	
	DECLARE_MEMBER_OBJECT(MonteCarlo, m_hMonteCarlo, ForwardSkewHandle);

public:
	STDMETHOD(Evaluate)(IResult* pVal);
	STDMETHOD(InterfaceSupportsErrorInfo)(REFIID riid);
};

#endif //__MOMENTUMCLIQUET_H_
