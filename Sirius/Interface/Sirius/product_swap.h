//	product_swap.h : Declaration of the CSwap
//
//////////////////////////////////////////////////////////////////////

#ifndef __SWAP_H_
#define __SWAP_H_

#include "resource.h"
#include "MlEqSwap.h"

class ATL_NO_VTABLE CSwap : 
	public CComObjectRootEx<CComSingleThreadModel>,
	public CComCoClass<CSwap, &CLSID_Swap>,
	public ISupportErrorInfo,
	public IDispatchImpl<ISwap, &IID_ISwap, &LIBID_Sirius>
{
public:	
	CSwap()	: m_h(new MlEqSwap) {}	
	DECLARE_REGISTRY_RESOURCEID(IDR_SWAP)
	DECLARE_NOT_AGGREGATABLE(CSwap)
	DECLARE_PROTECT_FINAL_CONSTRUCT()
	BEGIN_COM_MAP(CSwap)
		COM_INTERFACE_ENTRY(ISwap)
		COM_INTERFACE_ENTRY(IDispatch)
		COM_INTERFACE_ENTRY(ISupportErrorInfo)
	END_COM_MAP()
	associate_analytic_object(CSwap, MlEqSwap);
	STDMETHOD(InterfaceSupportsErrorInfo)(REFIID riid);
			
	STDMETHOD(Evaluate)(/*[in]*/ IResult* pVal);
	STDMETHOD(GetPV)(/*[in]*/ SwapLegEnum SwapLeg, /*[out, retval]*/ double* pVal);
	STDMETHOD(GetFixedPV)(/*[in]*/ SwapLegEnum SwapLeg, /*[out, retval]*/ double* pVal);
	STDMETHOD(GetFloatingPV)(/*[in]*/ SwapLegEnum SwapLeg, /*[out, retval]*/ double* pVal);

public:
	// A-leg
	declare_member_bool(A_FullFirst);
	declare_member_bool(A_RollForward);	

	declare_member_object(ZeroCurve, A_ZeroCurve);
	declare_member_object(Array, A_ResetRates);

	declare_member_string(A_AccrualCalendar);
	declare_member_string(A_FixedOrAccrualType);	
	declare_member_string(A_Frequency);
	declare_member_string(A_PayCalendar);
	declare_member_string(A_RateType);		
	declare_member_string(A_ResetLeadCalendar);
	declare_member_string(A_TermOrEndDate);

	declare_member_variable(BusinessDayConventionEnum, A_AccrualBusinessDayConvention);
	declare_member_variable(double, A_CouponOrMargin);
	declare_member_variable(DayCountConventionEnum, A_DayCountConvention);	
	declare_member_variable(SwapLegTypeEnum, A_LegType);	
	declare_member_variable(double, A_Notional);
	declare_member_variable(BusinessDayConventionEnum, A_PayBusinessDayConvention);
	declare_member_variable(PrincipalExchangedTypeEnum, A_PrincipalExchanged);
	declare_member_variable(long, A_ResetLead);
	declare_member_variable(DATE, A_RollDate);
	declare_member_variable(RollDateConventionEnum, A_RollDateConvention);
	declare_member_variable(DATE, A_StartDate);
	declare_member_variable(StubTypeEnum, A_LongStub);
		
	// B-leg
	declare_member_bool(B_FullFirst);
	declare_member_bool(B_RollForward);

	declare_member_object(ZeroCurve, B_ZeroCurve);
	declare_member_object(Array, B_ResetRates);

	declare_member_string(B_AccrualCalendar);
	declare_member_string(B_FixedOrAccrualType);
	declare_member_string(B_Frequency);
	declare_member_string(B_PayCalendar);
	declare_member_string(B_RateType);
	declare_member_string(B_ResetLeadCalendar);
	declare_member_string(B_TermOrEndDate);

	declare_member_variable(BusinessDayConventionEnum, B_AccrualBusinessDayConvention);
	declare_member_variable(double, B_CouponOrMargin);
	declare_member_variable(DayCountConventionEnum, B_DayCountConvention);
	declare_member_variable(SwapLegTypeEnum, B_LegType);
	declare_member_variable(double, B_Notional);
	declare_member_variable(BusinessDayConventionEnum, B_PayBusinessDayConvention);
	declare_member_variable(PrincipalExchangedTypeEnum, B_PrincipalExchanged);	
	declare_member_variable(long, B_ResetLead);
	declare_member_variable(DATE, B_RollDate);
	declare_member_variable(RollDateConventionEnum, B_RollDateConvention);
	declare_member_variable(DATE, B_StartDate);
	declare_member_variable(StubTypeEnum, B_LongStub);	
};

#endif