//	MlEqSwap.cpp :		Implementation of MlEqSwap
//
//	Authors :           David Cuin, Romain Barc
//
//  Style   :           K & R
//
/////////////////////////////////////////////////////////////////////////////

#include "stdafx.h"
#include "MlEqSwap.h"
#include "MlEqDate.h"
#include "MlEqZeroCurve.h"
#include "MlEqArray.h"


MlEqSwap::SwapParameters::SwapParameters(void)
{		
	m_bdcAccrual = NoBusinessDayConvention;
	m_bdcPay = NoBusinessDayConvention;			
	m_bFullFirst = false;
	m_bRollForward = false;
	m_dateRoll = 0.0;
	m_nStartDate = MlEqDate::GetCurrentDate();
	m_dcc = NoDayCountConvention;
	m_dccRoll = NoRollDateConvention;	
	m_fCouponOrMargin = 0.0;
	m_fNotional = 1.0;		
	m_LegType = NoSwapLegType;
	m_nResetLead = 0;	
	m_pe = NoPrincipalExchange;
	m_stLongStub = NoStubType;
	
	m_szAccrualCalendar;
	m_szFixedOrAccrualType;
	m_szFrequency = "SemiAnnual";	
}

MlEqSwap::MlEqSwap(void)
{
	// Do nothing
}

double MlEqSwap::GetFixedPV(SwapLegEnum SwapLeg) const
{
	double								f = 0.0;
	switch (SwapLeg){
	case A_Leg:
		f = GetFixedPV(m_spA.m_LegType, m_spA.m_hZeroCurve, m_spA.m_fNotional, m_spA.m_nStartDate, m_spA.m_szTerm, m_spA.m_fCouponOrMargin, m_spA.m_szFrequency, m_spA.m_dcc, m_spA.m_dateRoll, m_spA.m_dccRoll, m_spA.m_szAccrualCalendar, m_spA.m_bdcAccrual, m_spA.m_bFullFirst, m_spA.m_bRollForward, m_spA.m_stLongStub, m_spA.m_szPayCalendar, m_spA.m_bdcPay, m_spA.m_pe, m_spA.m_szFixedOrAccrualType);
		return f;
	case B_Leg:	
		f = GetFixedPV(m_spB.m_LegType, m_spB.m_hZeroCurve, m_spB.m_fNotional, m_spB.m_nStartDate, m_spB.m_szTerm, m_spB.m_fCouponOrMargin, m_spB.m_szFrequency, m_spB.m_dcc, m_spB.m_dateRoll, m_spB.m_dccRoll, m_spB.m_szAccrualCalendar, m_spB.m_bdcAccrual, m_spB.m_bFullFirst, m_spB.m_bRollForward, m_spB.m_stLongStub, m_spB.m_szPayCalendar, m_spB.m_bdcPay, m_spB.m_pe, m_spB.m_szFixedOrAccrualType);
		return f;
	default:
		throw "Invalid SwapLeg value '" + CEnumMap::GetString("SwapLegEnum", LIBID_Sirius, SwapLeg);
		return 0.0;
	}
}

double MlEqSwap::GetFixedPV(SwapLegTypeEnum SwapLegType, MlEqZeroCurveHandle hZeroCurve, double fNotional, long nStartDate, const estring& szTerm, double fCoupon, const estring& szFrequency, DayCountConventionEnum dcc, long dateRoll, RollDateConventionEnum dccRoll, const std::string& szAccrualCalendar, BusinessDayConventionEnum bdcAccrual, bool bFullFirst, bool bRollForward, StubTypeEnum stLongStub, const std::string& szPayCalendar, BusinessDayConventionEnum bdcPay, PrincipalExchangedTypeEnum pe, const estring& szFixedType) const
{
	double f = 0.0;
	if (SwapLegType == FixedLeg && fNotional){
		GDA::HDElement hdeParams(GDA::HDElement::Dictionary);
		if (!hZeroCurve) throw "No zero curve defined";		
		hdeParams.insert("YieldCurve", hZeroCurve->GetYieldCurve());
		if (fNotional) hdeParams.insert("Notional", GDA::HDElement(fNotional));
		if (nStartDate) hdeParams.insert("StartDate", GDA::HDElement(MlEqDate::ExcelDateToJulianDate(nStartDate)));
		if (szTerm.size()) hdeParams.insert("EndDate", TermOrDateToHDElement(szTerm));		
		if (fCoupon) hdeParams.insert("Coupon", GDA::HDElement(fCoupon));
		if (szFrequency.size()) hdeParams.insert("Freq", GDA::HDElement(szFrequency.data()));
		if (dcc != NoDayCountConvention) hdeParams.insert("DCM", GDA::HDElement(MlEqDate::SiriusToGDADayCountConvention(CEnumMap::GetString("DayCountConventionEnum", LIBID_Sirius, dcc)).data()));
		if (dateRoll) hdeParams.insert("RollDate", GDA::HDElement(MlEqDate::ExcelDateToJulianDate(dateRoll)));
		if (dccRoll != NoRollDateConvention) hdeParams.insert("RDC", GDA::HDElement(CEnumMap::GetString("RollDateConvention", LIBID_Sirius, dccRoll).data()));		
		if (szAccrualCalendar.size()) hdeParams.insert("AccCal", GDA::HDElement(szAccrualCalendar.data()));
		if (bdcAccrual != NoBusinessDayConvention) hdeParams.insert("AccBDC", GDA::HDElement(CEnumMap::GetString("BusinessDayConventionEnum", LIBID_Sirius, bdcAccrual).data()));		
		if (bFullFirst) hdeParams.insert("FullFirst", GDA::HDElement(TRUE));
		if (bRollForward) hdeParams.insert("RollForward", GDA::HDElement(TRUE));
		if (stLongStub != NoStubType) hdeParams.insert("LongStub", GDA::HDElement(CEnumMap::GetString("StubTypeEnum", LIBID_Sirius, stLongStub).data()));
		if (szPayCalendar.size()) hdeParams.insert("PayCal", GDA::HDElement(szPayCalendar.data()));
		if (bdcPay != NoBusinessDayConvention) hdeParams.insert("PayBDC", GDA::HDElement(CEnumMap::GetString("BusinessDayConventionEnum", LIBID_Sirius, bdcPay).data()));
		if (pe != NoPrincipalExchange) hdeParams.insert("PrinExch", GDA::HDElement(CEnumMap::GetString("PrincipalExchangedTypeEnum", LIBID_Sirius, pe).data()));	
		if (szFixedType.size()) hdeParams.insert("FixedType", GDA::HDElement(szFixedType.data()));
				
		GDA::HDElement hdeResult = GDA::GetLibraryInstance()->getFactory()->createAndApply("VanillaIR.FixedSide", hdeParams, GDA::HDElement("PV")).lookup("PV");
		if (hdeResult.isException()) hdeResult.raise();		
		f += hdeResult.asDouble();
	}
	return f;
}

double MlEqSwap::GetFloatingPV(SwapLegEnum SwapLeg) const
{	
	double								f = 0.0;
	switch (SwapLeg){
	case A_Leg:
		f = GetFloatingPV(m_spA.m_LegType, m_spA.m_hZeroCurve, m_spA.m_fNotional, m_spA.m_nStartDate, m_spA.m_szTerm, m_spA.m_fCouponOrMargin, m_spA.m_szFrequency, m_spA.m_dcc, m_spA.m_dateRoll, m_spA.m_dccRoll, m_spA.m_szAccrualCalendar, m_spA.m_bdcAccrual, m_spA.m_bRollForward, m_spA.m_stLongStub, m_spA.m_szPayCalendar, m_spA.m_bdcPay, m_spA.m_szRateType, m_spA.m_nResetLead, m_spA.m_szResetCalendar, m_spA.m_hResetRates, m_spA.m_pe, m_spA.m_szFixedOrAccrualType);
		return f;
	case B_Leg:
		f = GetFloatingPV(m_spB.m_LegType, m_spB.m_hZeroCurve, m_spB.m_fNotional, m_spB.m_nStartDate, m_spB.m_szTerm, m_spB.m_fCouponOrMargin, m_spB.m_szFrequency, m_spB.m_dcc, m_spB.m_dateRoll, m_spB.m_dccRoll, m_spB.m_szAccrualCalendar, m_spB.m_bdcAccrual, m_spB.m_bRollForward, m_spB.m_stLongStub, m_spB.m_szPayCalendar, m_spB.m_bdcPay, m_spB.m_szRateType, m_spB.m_nResetLead, m_spB.m_szResetCalendar, m_spB.m_hResetRates, m_spB.m_pe, m_spB.m_szFixedOrAccrualType);
		return f;
	default:
		throw "Invalid SwapLeg value '" + CEnumMap::GetString("SwapLegEnum", LIBID_Sirius, SwapLeg);
		return 0.0;
	}	
}
double MlEqSwap::GetFloatingPV(SwapLegTypeEnum SwapLegType, MlEqZeroCurveHandle hZeroCurve, double fNotional, long nStartDate, const estring& szTerm, double fMargin, const estring& szFrequency, DayCountConventionEnum dcc, long dateRoll, RollDateConventionEnum dccRoll, const std::string& szAccrualCalendar, BusinessDayConventionEnum bdcAccrual, bool bRollForward, StubTypeEnum stLongStub, const std::string& szPayCalendar, BusinessDayConventionEnum bdcPay, const estring& szRateType, long nResetLead, const std::string& szResetCalendar, MlEqArrayHandle hResetRates, PrincipalExchangedTypeEnum pe, const estring& szAccrualType) const
{
	double f = 0.0;
	if (SwapLegType == FloatingLeg && fNotional){
		GDA::HDElement hdeParams(GDA::HDElement::Dictionary);
		if (!hZeroCurve) throw "No zero curve defined";		
		hdeParams.insert("YieldCurve", hZeroCurve->GetYieldCurve());
		if (fNotional) hdeParams.insert("Notional", GDA::HDElement(fNotional));
		if (nStartDate) hdeParams.insert("StartDate", GDA::HDElement(MlEqDate::ExcelDateToJulianDate(nStartDate)));
		if (szTerm.size()) hdeParams.insert("EndDate", TermOrDateToHDElement(szTerm));
		if (fMargin) hdeParams.insert("Margin", GDA::HDElement(fMargin));
		if (szFrequency.size()) hdeParams.insert("Freq", GDA::HDElement(szFrequency.data()));
		if (dcc != NoDayCountConvention) hdeParams.insert("DCM", GDA::HDElement(MlEqDate::SiriusToGDADayCountConvention(CEnumMap::GetString("DayCountConventionEnum", LIBID_Sirius, dcc)).data()));
		if (dateRoll) hdeParams.insert("RollDate", GDA::HDElement(MlEqDate::ExcelDateToJulianDate(dateRoll)));
		if (dccRoll != NoRollDateConvention) hdeParams.insert("RDC", GDA::HDElement(CEnumMap::GetString("RollDateConvention", LIBID_Sirius, dccRoll).data()));		
		if (szAccrualCalendar.size()) hdeParams.insert("AccCal", GDA::HDElement(szAccrualCalendar.data()));
		if (bdcAccrual != NoBusinessDayConvention) hdeParams.insert("AccBDC", GDA::HDElement(CEnumMap::GetString("BusinessDayConventionEnum", LIBID_Sirius, bdcAccrual).data()));
		if (bRollForward) hdeParams.insert("RollForward", GDA::HDElement(TRUE));
		if (stLongStub != NoStubType) hdeParams.insert("LongStub", GDA::HDElement(CEnumMap::GetString("StubTypeEnum", LIBID_Sirius, stLongStub).data()));
		if (szPayCalendar.size()) hdeParams.insert("PayCal", GDA::HDElement(szPayCalendar.data()));
		if (bdcPay != NoBusinessDayConvention) hdeParams.insert("PayBDC", GDA::HDElement(CEnumMap::GetString("BusinessDayConventionEnum", LIBID_Sirius, bdcPay).data()));
		if (szRateType.size()) hdeParams.insert("RateType", GDA::HDElement(szRateType.data()));
		if (nResetLead) hdeParams.insert("ResetLead", GDA::HDElement(nResetLead));
		if (szResetCalendar.size()) hdeParams.insert("ResetLeadCal", GDA::HDElement(szResetCalendar.data()));
		if (!!hResetRates && hResetRates->getsize()){
			hdeParams.insert("ResetRates", *hResetRates);
		}
		if (pe != NoPrincipalExchange) hdeParams.insert("PrinExch", GDA::HDElement(CEnumMap::GetString("PrincipalExchangedTypeEnum", LIBID_Sirius, pe).data()));	
		if (szAccrualType.size()) hdeParams.insert("AccType", GDA::HDElement(szAccrualType.data()));
		GDA::HDElement hdeResult = GDA::GetLibraryInstance()->getFactory()->createAndApply("VanillaIR.FloatSide", hdeParams, GDA::HDElement("PV")).lookup("PV");
		if (hdeResult.isException()) hdeResult.raise();		
		f += hdeResult.asDouble();
	}
	return f;
}
double MlEqSwap::GetPV(SwapLegEnum SwapLeg) const
{
	return GetFixedPV(SwapLeg) + GetFloatingPV(SwapLeg);
}


//////////////////////////////////////////////////////////////////////////////
//	TermOrDateToHDElement
//
//	Maps a string input s.t. if the input is a date then we convert it
//	to a Julian Date format.
//
/*static*/ GDA::HDElement MlEqSwap::TermOrDateToHDElement(const std::string& sz)
{
	CParameterMap						pm(sz);
	DATE								date;
	GDA::HDElement						hde;

	if (!pm.GetValue(&date)){
		hde = CParameterMap::ExcelDateToJulianDate(date);		
	} else {
		hde = sz.data();
	}
	return hde;		
}


/*ToDo - unused code:
double MlEqSwap::_GetForwardValue(long nStartDate)
{
	double	f = 0.0;
	f += GetFloatingPV(m_spA.m_LegType, m_spA.m_hZeroCurve, m_spA.m_fNotional, nStartDate, m_spA.m_szTerm, m_spA.m_fCouponOrMargin, m_spA.m_szFrequency, m_spA.m_dcc, m_spA.m_dateRoll, m_spA.m_dccRoll, m_spA.m_szAccrualCalendar, m_spA.m_bdcAccrual, m_spA.m_bRollForward, m_spA.m_stLongStub, m_spA.m_szPayCalendar, m_spA.m_bdcPay, m_spA.m_szRateType, m_spA.m_nResetLead, m_spA.m_szResetCalendar, m_spA.m_hResetRates, m_spA.m_pe, m_spA.m_szFixedOrAccrualType);
	f += GetFloatingPV(m_spB.m_LegType, m_spB.m_hZeroCurve, m_spB.m_fNotional, nStartDate, m_spB.m_szTerm, m_spB.m_fCouponOrMargin, m_spB.m_szFrequency, m_spB.m_dcc, m_spB.m_dateRoll, m_spB.m_dccRoll, m_spB.m_szAccrualCalendar, m_spB.m_bdcAccrual, m_spB.m_bRollForward, m_spB.m_stLongStub, m_spB.m_szPayCalendar, m_spB.m_bdcPay, m_spB.m_szRateType, m_spB.m_nResetLead, m_spB.m_szResetCalendar, m_spB.m_hResetRates, m_spB.m_pe, m_spB.m_szFixedOrAccrualType);
	f += GetFixedPV(m_spA.m_LegType, m_spA.m_hZeroCurve, m_spA.m_fNotional, nStartDate, m_spA.m_szTerm, m_spA.m_fCouponOrMargin, m_spA.m_szFrequency, m_spA.m_dcc, m_spA.m_dateRoll, m_spA.m_dccRoll, m_spA.m_szAccrualCalendar, m_spA.m_bdcAccrual, m_spA.m_bFullFirst, m_spA.m_bRollForward, m_spA.m_stLongStub, m_spA.m_szPayCalendar, m_spA.m_bdcPay, m_spA.m_pe, m_spA.m_szFixedOrAccrualType);
	f += GetFixedPV(m_spB.m_LegType, m_spB.m_hZeroCurve, m_spB.m_fNotional, nStartDate, m_spB.m_szTerm, m_spB.m_fCouponOrMargin, m_spB.m_szFrequency, m_spB.m_dcc, m_spB.m_dateRoll, m_spB.m_dccRoll, m_spB.m_szAccrualCalendar, m_spB.m_bdcAccrual, m_spB.m_bFullFirst, m_spB.m_bRollForward, m_spB.m_stLongStub, m_spB.m_szPayCalendar, m_spB.m_bdcPay, m_spB.m_pe, m_spB.m_szFixedOrAccrualType);	
	return f;
}

double MlEqSwap::GetForwardValue(long nStartDate, long nEndDate)
{
	return _GetForwardValue(nStartDate) - _GetForwardValue(nEndDate);
}*/