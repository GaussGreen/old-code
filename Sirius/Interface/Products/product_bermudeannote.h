// product_bermudeanNote.h : Declaration of the CBermudeanNote

#ifndef __BERMUDEANNOTE_H_
#define __BERMUDEANNOTE_H_

#include "resource.h"       // main symbols
#include "mleqobjects.h"

 



class ATL_NO_VTABLE CBermudeanNote : 
	public CComObjectRootEx<CComSingleThreadModel>,
	public CComCoClass<CBermudeanNote, &CLSID_BermudeanNote>,
	public IDispatchImpl<IBermudeanNote, &IID_IBermudeanNote, &LIBID_Products>
{
public:
	HRESULT FinalConstruct(void);

	DECLARE_REGISTRY_RESOURCEID(IDR_BERMUDEANNOTE)
	DECLARE_NOT_AGGREGATABLE(CBermudeanNote)

	DECLARE_PROTECT_FINAL_CONSTRUCT()

	BEGIN_COM_MAP(CBermudeanNote)
		COM_INTERFACE_ENTRY(IBermudeanNote)
		COM_INTERFACE_ENTRY(IDispatch)
	END_COM_MAP()

	CBermudeanNote(){}

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
	
	DECLARE_MEMBER_OBJECT(DateSchedule, m_hCallDates, CallDates);
	DECLARE_MEMBER_OBJECT(Array, m_hRebates, Rebates)		


public:	
	STDMETHOD(Evaluate)(IResult* pVal);
	STDMETHOD(InterfaceSupportsErrorInfo)(REFIID riid);
};

#endif //__BERMUDEANNOTE_H_
