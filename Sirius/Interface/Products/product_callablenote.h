// product_callableNote.h : Declaration of the CCallableNote

#ifndef __CALLABLENOTE_H_
#define __CALLABLENOTE_H_

#include "resource.h"       // main symbols
#include "mleqobjects.h"

/////////////////////////////////////////////////////////////////////////////
// CCallableNote
class ATL_NO_VTABLE CCallableNote : 
	public CComObjectRootEx<CComSingleThreadModel>,
	public CComCoClass<CCallableNote, &CLSID_CallableNote>,
	public IDispatchImpl<ICallableNote, &IID_ICallableNote, &LIBID_Products>
{
public:
	HRESULT FinalConstruct(void);

	DECLARE_REGISTRY_RESOURCEID(IDR_CALLABLENOTE)
	DECLARE_NOT_AGGREGATABLE(CCallableNote)
	DECLARE_PROTECT_FINAL_CONSTRUCT()

	BEGIN_COM_MAP(CCallableNote)
		COM_INTERFACE_ENTRY(ICallableNote)
		COM_INTERFACE_ENTRY(IDispatch)
	END_COM_MAP()

	CCallableNote(){}

	DECLARE_MEMBER_VARIABLE(double, m_fRefSpot, ReferenceSpot)
		
	DECLARE_MEMBER_VARIABLE(double, m_fStrike, Strike)
	DECLARE_MEMBER_VARIABLE(double, m_fGearing, Gearing)
	DECLARE_MEMBER_VARIABLE(PayoffTypeEnum, m_cpe, CallOrPut)
	DECLARE_MEMBER_VARIABLE(double, m_fFixedCoupon, FixedCoupon)
	DECLARE_MEMBER_VARIABLE(DATE, m_dateMaturity, Maturity)

	DECLARE_MEMBER_OBJECT(Asset, m_hUnderlying, Underlying)
	DECLARE_MEMBER_OBJECT(DateSchedule, m_hResetSchedule, ResetSchedule);
	DECLARE_MEMBER_OBJECT(DateSchedule, m_hPaySchedule, PaySchedule);
	DECLARE_MEMBER_OBJECT(Array, m_hLiborFixings, CouponFixings)		
	DECLARE_MEMBER_OBJECT(ZeroCurve, m_hZeroCurve, ZeroCurve)

	DECLARE_MEMBER_VARIABLE(DayCountConventionEnum, m_dcc, DayCountConvention)
	DECLARE_MEMBER_OBJECT(Array, m_hSpread, SpreadArray)			
	DECLARE_MEMBER_OBJECT(Array, m_hWeights, WeightsArray)		

public:	
	STDMETHOD(Evaluate)(IResult* pVal);
	STDMETHOD(InterfaceSupportsErrorInfo)(REFIID riid);
protected:
};

#endif //__CALLABLENOTE_H_
