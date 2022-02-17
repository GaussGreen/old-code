//	mleqswap.h :	  Swap rate class. (Essentially a GDA wrapper.)
//
//	Authors :         David Cuin, Romain Barc
/////////////////////////////////////////////////////////////////////////////

#ifndef __MLEQSWAP_H_
#define __MLEQSWAP_H_

#include "smart.h"

#define _swap_parameter(Leg, DataType, VariableName, FunctionName)				\
	DataType			Get##Leg##_##FunctionName(void) const					\
	{																			\
		return m_sp##Leg##.##VariableName;										\
	}																			\
	void				Put##Leg##_##FunctionName(const DataType& d)			\
	{																			\
		m_sp##Leg##.##VariableName = d;											\
	}

#define swap_parameter(DataType, VariableName, FunctionName)					\
	_swap_parameter(A, DataType, VariableName, FunctionName)					\
	_swap_parameter(B, DataType, VariableName, FunctionName)


class MlEqSwap: public RCObject
{
public:
	MlEqSwap(void);	

	swap_parameter(BusinessDayConventionEnum, m_bdcAccrual, AccrualBusinessDayConvention);
	swap_parameter(BusinessDayConventionEnum, m_bdcPay, PayBusinessDayConvention);
	swap_parameter(bool, m_bFullFirst, FullFirst);
	swap_parameter(bool, m_bRollForward, RollForward);
	swap_parameter(long, m_dateRoll, RollDate);
	swap_parameter(long, m_nStartDate, StartDate);	
	swap_parameter(DayCountConventionEnum, m_dcc, DayCountConvention);
	swap_parameter(RollDateConventionEnum, m_dccRoll, RollDateConvention);
	swap_parameter(MlEqZeroCurveHandle, m_hZeroCurve, ZeroCurve);
	swap_parameter(double, m_fCouponOrMargin, CouponOrMargin);
	swap_parameter(double, m_fNotional, Notional);
	swap_parameter(SwapLegTypeEnum, m_LegType, LegType);
	swap_parameter(long, m_nResetLead, ResetLead);
	swap_parameter(MlEqArrayHandle, m_hResetRates, ResetRates);
	swap_parameter(PrincipalExchangedTypeEnum, m_pe, PrincipalExchanged);
	swap_parameter(StubTypeEnum, m_stLongStub, LongStub);
	swap_parameter(std::string, m_szAccrualCalendar, AccrualCalendar);
	swap_parameter(std::string, m_szFixedOrAccrualType, FixedOrAccrualType);
	swap_parameter(std::string, m_szFrequency, Frequency);
	swap_parameter(std::string, m_szPayCalendar, PayCalendar);
	swap_parameter(std::string, m_szRateType, RateType);
	swap_parameter(std::string, m_szResetCalendar, ResetLeadCalendar);
	swap_parameter(std::string, m_szTerm, TermOrEndDate);

	double								GetPV(SwapLegEnum SwapLeg) const;
	double								GetFixedPV(SwapLegEnum SwapLeg) const;
	double								GetFloatingPV(SwapLegEnum SwapLeg) const;
	
protected:
	struct SwapParameters
	{
		SwapParameters(void);
		
		BusinessDayConventionEnum		m_bdcAccrual;
		BusinessDayConventionEnum		m_bdcPay;		
		bool							m_bFullFirst;						
		bool							m_bRollForward;
		long							m_dateRoll;		
		long							m_nStartDate;
		DayCountConventionEnum			m_dcc;				
		RollDateConventionEnum			m_dccRoll;
		MlEqZeroCurveHandle				m_hZeroCurve;		
		double							m_fCouponOrMargin;
		double							m_fNotional;				
		SwapLegTypeEnum					m_LegType;
		long							m_nResetLead;		
		MlEqArrayHandle					m_hResetRates;
		PrincipalExchangedTypeEnum		m_pe;
		StubTypeEnum					m_stLongStub;		
		std::string						m_szAccrualCalendar;
		std::string						m_szFixedOrAccrualType;
		std::string						m_szFrequency;
		std::string						m_szPayCalendar;
		std::string						m_szRateType;
		std::string						m_szResetCalendar;
		std::string						m_szTerm;		
	};
	
	SwapParameters						m_spA;
	SwapParameters						m_spB;		
		
protected:
	static GDA::HDElement				TermOrDateToHDElement(const std::string& sz);		
	double								GetFixedPV(SwapLegTypeEnum SwapLegType, MlEqZeroCurveHandle hZeroCurve, double fNotional, long nStartDate, const estring& szTerm, double fCoupon, const estring& szFrequency, DayCountConventionEnum dcc, long dateRoll, RollDateConventionEnum dccRoll, const std::string& szAccrualCalendar, BusinessDayConventionEnum bdcAccrual, bool bFullFirst, bool bRollForward, StubTypeEnum stLongStub, const std::string& szPayCalendar, BusinessDayConventionEnum bdcPay, PrincipalExchangedTypeEnum pe, const estring& szFixedType) const;
	double								GetFloatingPV(SwapLegTypeEnum LegType, MlEqZeroCurveHandle hZeroCurve, double fNotional, long nStartDate, const estring& szTerm, double fMargin, const estring& szFrequency, DayCountConventionEnum dcc, long dateRoll, RollDateConventionEnum dccRoll, const std::string& szAccrualCalendar, BusinessDayConventionEnum bdcAccrual, bool bRollForward, StubTypeEnum stLongStub, const std::string& szPayCalendar, BusinessDayConventionEnum bdcPay, const estring& szRateType, long nResetLead, const std::string& szResetCalendar, MlEqArrayHandle hResetRates, PrincipalExchangedTypeEnum pe, const estring& szAccrualType) const;
	
	/*ToDo - Unused code
	double								_GetForwardValue(long nStartDate);
	double								GetForwardValue(long nStartDate, long nEndDate);
	MlEqZeroCurveHandle					m_payhZeroCurve;	// set it to leg A by default...
	double								m_fxA;
	double								m_fxB;*/

};

#undef _swap_parameter
#undef swap_parameter
	
#endif

