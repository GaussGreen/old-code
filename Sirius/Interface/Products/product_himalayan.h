//	product_himalayan.h : Declaration of the CHimalayan
//
/////////////////////////////////////////////////////////////////////////////

#ifndef __HIMALAYAN_H_
#define __HIMALAYAN_H_

#include "resource.h"
#include "cappedflooredcliquetpricer.h"
#include "MlEqObjects.h"

class ATL_NO_VTABLE CHimalayan : 
	public CComObjectRootEx<CComSingleThreadModel>,
	public CComCoClass<CHimalayan, &CLSID_Himalayan>,
	public ISupportErrorInfo,
	public IDispatchImpl<IHimalayan, &IID_IHimalayan, &LIBID_Products>,
	public Himalaya
{
public:
	HRESULT								FinalConstruct(void);
	CHimalayan(){}
	DECLARE_REGISTRY_RESOURCEID(IDR_HIMALAYAN)
	DECLARE_NOT_AGGREGATABLE(CHimalayan)
	DECLARE_PROTECT_FINAL_CONSTRUCT()
	BEGIN_COM_MAP(CHimalayan)
		COM_INTERFACE_ENTRY(IHimalayan)
		COM_INTERFACE_ENTRY(IDispatch)
		COM_INTERFACE_ENTRY(ISupportErrorInfo)
	END_COM_MAP()
		
	DECLARE_MEMBER_OBJECT(DateSchedule, m_hFixingSchedule, DateAndFixingSchedule)
	DECLARE_MEMBER_OBJECT(Array, m_hRainbowWeights, RainbowWeights)
	DECLARE_MEMBER_VARIABLE(double, m_fNotional, Notional)
	DECLARE_MEMBER_VARIABLE(double, m_fStrike, Strike)
	DECLARE_MEMBER_VARIABLE(PayoffTypeEnum, m_CallOrPut, CallOrPut)
	DECLARE_MEMBER_OBJECT(Array, m_hCallLevelsArr, CallLevels)	
	DECLARE_MEMBER_OBJECT(Array, m_hRebates, Rebates)
	DECLARE_MEMBER_BOOL(m_bDelayRebateToEnd, DelayRebateToEnd)
	DECLARE_MEMBER_BOOL(m_bAsianFromStart, AsianFromStart)
	DECLARE_MEMBER_OBJECT(Array, m_hIsHimalayaArr, IsHimalayaArr)	
	DECLARE_MEMBER_OBJECT(ParameterList, m_hParameterList, ModelInfo)		
	declare_member_collection(Assets, Asset, m_hAssets, Underlyings)
	DECLARE_MEMBER_VARIABLE(DATE, m_dateFixing, FixingDate)
	DECLARE_MEMBER_OBJECT(Array, m_hFixingValues, FixingValues);
			
public:	
	STDMETHOD(Evaluate)(IResult* pVal);
	STDMETHOD(InterfaceSupportsErrorInfo)(REFIID riid);
};

#endif