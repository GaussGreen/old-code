// product_autocallableNote.h : Declaration of the CAutocallableNote

#ifndef __AUTOCALLABLENOTE_H_
#define __AUTOCALLABLENOTE_H_

#include "resource.h"       // main symbols
#include "mleqobjects.h"
#include "callablecliquetswappricer.h"

/////////////////////////////////////////////////////////////////////////////
// CAutocallableNote
class ATL_NO_VTABLE CAutocallableNote : 
	public CComObjectRootEx<CComSingleThreadModel>,
	public CComCoClass<CAutocallableNote, &CLSID_AutocallableNote>,
	public IDispatchImpl<IAutocallableNote, &IID_IAutocallableNote, &LIBID_Products>,
	public AutoCallableSwapPricer
{
public:
	HRESULT FinalConstruct(void);

	CAutocallableNote(){}

	DECLARE_REGISTRY_RESOURCEID(IDR_AUTOCALLABLENOTE)
	DECLARE_NOT_AGGREGATABLE(CAutocallableNote)

	DECLARE_PROTECT_FINAL_CONSTRUCT()

	BEGIN_COM_MAP(CAutocallableNote)
		COM_INTERFACE_ENTRY(IAutocallableNote)
		COM_INTERFACE_ENTRY(IDispatch)
	END_COM_MAP()

	DECLARE_MEMBER_VARIABLE(double, m_fStrike, Strike)
	DECLARE_MEMBER_VARIABLE(double, m_fGearing, Gearing)
	DECLARE_MEMBER_VARIABLE(PayoffTypeEnum, m_cpe, CallOrPut)
	DECLARE_MEMBER_VARIABLE(double, m_fFixedCoupon, FixedCoupon)

	DECLARE_MEMBER_OBJECT(Asset, m_hUnderlying, Underlying)
	DECLARE_MEMBER_OBJECT(DateSchedule, m_hResetSchedule, ResetSchedule);
	DECLARE_MEMBER_OBJECT(DateSchedule, m_hPaySchedule, PaySchedule);
	DECLARE_MEMBER_OBJECT(Array, m_hLiborFixings, CouponFixings)		
	DECLARE_MEMBER_OBJECT(ZeroCurve, m_hZeroCurve, ZeroCurve)

	DECLARE_MEMBER_VARIABLE(DayCountConventionEnum, m_dcc, DayCountConvention)
	DECLARE_MEMBER_OBJECT(Array, m_hSpread, SpreadArray)			
	DECLARE_MEMBER_OBJECT(Array, m_hWeights, WeightsArray)		

	DECLARE_MEMBER_OBJECT(DateSchedule, m_hBarrierSchedule, BarrierSchedule);
	DECLARE_MEMBER_OBJECT(Array, m_hCallLevels, CallLevels)		
	DECLARE_MEMBER_OBJECT(Array, m_hRebates, Rebates)	

	DECLARE_MEMBER_OBJECT(Array, m_hThresholds, Thresholds)		
	DECLARE_MEMBER_OBJECT(Array, m_hCouponLow, CouponLow)	
	DECLARE_MEMBER_OBJECT(Array, m_hCouponUp, CouponUp)	

	DECLARE_MEMBER_OBJECT(Array, m_hRainbowWeights, RainbowWeights)	
	
	DECLARE_MEMBER_OBJECT(ParameterList, m_hParameterList, ModelInfo);

public:	
	STDMETHOD(Evaluate)(IResult* pVal);
	STDMETHOD(InterfaceSupportsErrorInfo)(REFIID riid);
protected:
};

#endif //__AUTOCALLABLENOTE_H_
