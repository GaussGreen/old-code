
#pragma warning(disable :4005 4786)

#include <libCCxll\CCxll.h>

#include <ARM\libarm_local\ARM_local_glob.h>
#include <ARM\libarm_local\ARM_local_class.h>

#include <ARM\libicm_local\ICM_local_frn.h>

#include "ARM_local_interface.h"
#include "ARM_local_interglob.h"

#include "ARM_xl_trycatch_local.h"
#include <util\fromto.h>
#include "ExcelTools.h"

__declspec(dllexport) LPXLOPER WINAPI ARM_Credit_FRN (LPXLOPER XL_Spread, 
													  LPXLOPER XL_StartDate,
													  LPXLOPER XL_Maturity,
													  LPXLOPER XL_Index,
													  LPXLOPER XL_First_period_refdate,
													  LPXLOPER XL_InitialRate,
													  LPXLOPER XL_LastIndexFixing,
													  LPXLOPER XL_AccOnDef,	
													  LPXLOPER XL_NotAmount,	
													  LPXLOPER XL_DayCountFrq,
													  LPXLOPER XL_AccruedDayCount,
													  LPXLOPER XL_SettlementGap,
													  LPXLOPER XL_Currency,
													  LPXLOPER XL_ResetCal,
													  LPXLOPER XL_PayCal)

{
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{

	ARM_NOCALCIFWIZ();

	// C variable
	

	double C_Spread = 0.;
	double C_StartDate = 0.;
	double C_Maturity = 0.;

	CCString C_Index;
	double C_First_period_refdate = 0.;
	double C_First_period_refdate_default = -1.;

	double C_InitialRate = 0.;
	double C_InitialRate_default = -999.;

	double C_LastIndexFixing = 0.;
	double C_LastIndexFixing_default = 0.;


	// CCString AccOnDef;
	// int C_AccOnDef = 0;

	CCString C_AmortOnDef;
	double AmortOnDef = 0.;

	CCString C_IntOnDef;
	double IntOnDef = 0.;

	double C_NotAmount = 0.;
	double C_NotAmount_default = 1000000.;

	CCString DayCountFrq;
	int C_DayCountFrq = 0;

	CCString AccruedDayCount;
	int C_AccruedDayCount = 0;

	double C_SettlementGap;
	double C_SettlementGap_default = 2;

	CCString C_Currency;
	CCString C_ResetCal;
	CCString C_PayCal;

	// error
	static int error;
	static char* reason = "";

	XL_readNumCell(XL_Spread,C_Spread," ARM_ERR: Spread: numeric expected",C_result);
	XL_readNumCell(XL_StartDate,C_StartDate," ARM_ERR: Start date: date expected",C_result);
	XL_readNumCell(XL_Maturity,C_Maturity," ARM_ERR: Maturity date: date expected",C_result);
	XL_readStrCell(XL_Index,C_Index," ARM_ERR: IRIndex: object IRIndex expected",C_result);
	XL_readNumCellWD(XL_First_period_refdate,C_First_period_refdate,C_First_period_refdate_default," ARM_ERR: First_period_refdate Rate: date expected",C_result);
	XL_readNumCellWD(XL_InitialRate,C_InitialRate,C_InitialRate_default," ARM_ERR: InitialRate: numeric expected",C_result);
	XL_readNumCellWD(XL_LastIndexFixing,C_LastIndexFixing,C_LastIndexFixing_default," ARM_ERR: LastFixing: numeric expected",C_result);
	// XL_readStrCellWD(XL_AccOnDef,AccOnDef,"ACC"," ARM_ERR: currency: string expected",C_result);
	qPAYMENT_PREMIUM_LEG C_AccOnDef; 
	ExcelTools::econvert(XL_AccOnDef,"ACC",C_AccOnDef); 

	XL_readNumCellWD(XL_NotAmount,C_NotAmount,C_NotAmount_default," ARM_ERR: NotAmountIn: numeric expected",C_result);
	XL_readStrCellWD(XL_DayCountFrq,DayCountFrq,"ACTUAL"," ARM_ERR: DayCountFrq: numeric expected",C_result);
	XL_readStrCellWD(XL_AccruedDayCount,AccruedDayCount,"ACTUAL"," ARM_ERR: AccruedDayCount: numeric expected",C_result);
	XL_readNumCellWD(XL_SettlementGap,C_SettlementGap,C_SettlementGap_default," ARM_ERR: SettlementGap: numeric expected",C_result);
	XL_readStrCellWD(XL_Currency,C_Currency,"DEFAULT"," ARM_ERR: currency: string expected",C_result);
	XL_readStrCellWD(XL_ResetCal,C_ResetCal,"NULL"," ARM_ERR: Reset Calendar: string expected",C_result);
	XL_readStrCellWD(XL_PayCal,C_PayCal,"NULL"," ARM_ERR: Payment Calendar: string expected",C_result);


	ARM_result currencyres;
	ARMLOCAL_ARM_GetDefaultCurrency (currencyres);

	if(C_Currency == "DEFAULT")
	{
		if(currencyres.getRetCode () != ARM_OK)
		{
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}
		else
		{
			C_Currency = currencyres.getString ();
		}
	}


	long retCode;

	if((C_DayCountFrq = ARM_ConvDayCount (DayCountFrq)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	if((C_AccruedDayCount = ARM_ConvDayCount (AccruedDayCount)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}
	
	// if((C_AccOnDef = ARM_ConvAccOnDef (AccOnDef, C_result)) == ARM_DEFAULT_ERR)
	// {
	// 	ARM_ARG_ERR();
	// 	return (LPXLOPER)&XL_result;
	// }

	if (C_SettlementGap>30)
	{
		C_result.setMsg("Error: Settlement must be <= 30");
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	long objId;
	CCString prevClass;
	
	CCString curClass = LOCAL_FRN_CLASS;
	CCString stringId = GetLastCurCellEnvValue ();
	
	if(!stringId)
	{	
		retCode = ICMLOCAL_FRN(C_Spread,
							   C_StartDate,
							   C_Maturity,
							   LocalGetNumObjectId(C_Index),
							   C_InitialRate,
							   C_LastIndexFixing,
							   C_AccOnDef,
							   C_NotAmount ,
							   C_DayCountFrq,
							   C_AccruedDayCount,
							   C_SettlementGap,	
							   C_Currency,
							   C_ResetCal,
							   C_PayCal,
							   C_First_period_refdate,
							   C_result);
		

		if(retCode == ARM_OK)
		{
			objId = C_result.getLong ();

			LocalSetCurCellEnvValue (curClass, objId); 

			stringId = LocalMakeObjectId (objId, curClass);
		}
	}
	else
	{
		prevClass = LocalGetStringObjectClass (stringId);
		
		objId = LocalGetNumObjectId (stringId);
			
		if(curClass == prevClass)
		{
		retCode = ICMLOCAL_FRN(C_Spread,
							   C_StartDate,
							   C_Maturity,
							   LocalGetNumObjectId(C_Index),
							   C_InitialRate,
							   C_LastIndexFixing,
							   C_AccOnDef,
							   C_NotAmount ,
							   C_DayCountFrq,
							   C_AccruedDayCount,
							   C_SettlementGap,	
							   C_Currency,
							   C_ResetCal,
							   C_PayCal,
							   C_First_period_refdate,
							   C_result,
							   objId);
		
			if(retCode == ARM_OK)
			{
				LocalSetCurCellEnvValue (curClass, objId); 

				stringId = LocalMakeObjectId (objId, curClass);
			}
		}
		else
		{
			FreeCurCellContent ();

		retCode = ICMLOCAL_FRN(C_Spread,
							   C_StartDate,
							   C_Maturity,
							   LocalGetNumObjectId(C_Index),
							   C_InitialRate,
							   C_LastIndexFixing,
							   C_AccOnDef,
							   C_NotAmount ,
							   C_DayCountFrq,
							   C_AccruedDayCount,
							   C_SettlementGap,	
							   C_Currency,
							   C_ResetCal,
							   C_PayCal,
							   C_First_period_refdate,
							   C_result);

			
			if(retCode == ARM_OK)
			{
				objId = C_result.getLong ();
			
				LocalSetCurCellEnvValue (curClass, objId); 

				stringId = LocalMakeObjectId (objId, curClass);
			}
		}
	}

	if(retCode == ARM_OK)
	{
		FreeCurCellErr ();
		XL_result.xltype = xltypeStr;
		XL_result.val.str = XL_StrC2StrPascal (stringId);
		XL_result.xltype |= xlbitDLLFree;
	}
	else
	{
		ARM_ERR();
	}

//	ARM_END();
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in ARM_Credit_FRN" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


__declspec(dllexport) LPXLOPER WINAPI ARM_Credit_CLN (LPXLOPER XL_StartDate,
													  LPXLOPER XL_Maturity,
													  LPXLOPER XL_Spread, 
													  LPXLOPER XL_IndexObjectId,
													  LPXLOPER XL_refdate,
													  LPXLOPER XL_fstcpnrefdate,
													  LPXLOPER XL_Notional,	
													  LPXLOPER XL_AccOnDef,	
													  LPXLOPER XL_DayCountFrq,	
													  LPXLOPER XL_DecompFreq,
													  LPXLOPER XL_StubRule,
													  LPXLOPER XL_ResetGap,
													  LPXLOPER XL_Currency,
													  LPXLOPER XL_ResetCalendar,
													  LPXLOPER XL_PayCalendar,
													  LPXLOPER XL_NxChange,
													  LPXLOPER XL_IncludeMaturity,
													  LPXLOPER XL_adjstartdate,
													  LPXLOPER XL_Binary)
{
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();

	// C variable
	
	CCString C_IndexObjectId;

	double C_Spread = 0.;
	double C_StartDate = 0.;
	double C_Maturity = 0.;

	double C_refdate = 0.;
	double C_refdate_default = -1.;

	double C_fstcpnrefdate = 0.;
	double C_fstcpnrefdate_default = -1.;

	double C_resetgap = 0.;
	double C_resetgap_default = 1.;

	// CCString AccOnDef;
	// int C_AccOnDef = 0;

	double C_Notional = 0.;
	double C_Notional_default = 1000000.;

	CCString DayCountFrq;
	int C_DayCountFrq = 0;

	CCString DecompFreq;
	int C_DecompFreq = 0;

	CCString StubRule;
	int C_StubRule = 0;

	CCString NxChange;
	int C_NxChange = 0;

	CCString IncludeMaturity;
	bool C_IncludeMaturity = false;

	CCString adjstartdate;
	int C_adjstartdate=0;

	CCString C_Currency;
	CCString C_ResetCalendar;
	CCString C_PayCalendar;

	double C_Binary=0.;
	double C_Binary_default=-999.;

	// error
	static int error;
	static char* reason = "";

	 XL_readNumCell(XL_StartDate,C_StartDate," ARM_ERR: Start date: date expected",C_result);
 	 XL_readNumCell(XL_Maturity,C_Maturity," ARM_ERR: Maturity date: date expected",C_result);
	 XL_readNumCell(XL_Spread,C_Spread," ARM_ERR: Spread: numeric expected",C_result);
	 XL_readStrCellWD(XL_IndexObjectId,C_IndexObjectId,"NONE"," ARM_ERR: IRIndex: object IRIndex expected",C_result);
	 XL_readNumCellWD(XL_refdate,C_refdate,C_refdate_default," ARM_ERR: reference date : date expected",C_result);
	 XL_readNumCellWD(XL_fstcpnrefdate,C_fstcpnrefdate,C_fstcpnrefdate_default," ARM_ERR: First period refdate : date expected",C_result);
	 XL_readNumCellWD(XL_Notional,C_Notional,C_Notional_default," ARM_ERR: Notional: numeric expected",C_result);
	 // XL_readStrCellWD(XL_AccOnDef,AccOnDef,"ACC"," ARM_ERR: currency: string expected",C_result);
	qPAYMENT_PREMIUM_LEG C_AccOnDef; 
	ExcelTools::econvert(XL_AccOnDef,"ACC",C_AccOnDef); 

	 XL_readStrCellWD(XL_DayCountFrq,DayCountFrq,"ACTUAL"," ARM_ERR: DayCountFrq: numeric expected",C_result);
	 XL_readStrCellWD(XL_DecompFreq,DecompFreq,"P"," ARM_ERR: DecompFreq: string expected",C_result);
	 XL_readStrCellWD(XL_StubRule,StubRule,"SS"," ARM_ERR: Stub Rule expected",C_result);
	 XL_readNumCellWD(XL_ResetGap,C_resetgap,C_resetgap_default," ARM_ERR: reset gap: numeric expected",C_result);
	 XL_readStrCellWD(XL_Currency,C_Currency,"DEFAULT"," ARM_ERR: currency: string expected",C_result);
	 XL_readStrCellWD(XL_ResetCalendar,C_ResetCalendar,"EUR"," ARM_ERR: Reset Calendar: string expected",C_result);
	 XL_readStrCellWD(XL_PayCalendar,C_PayCalendar,"EUR"," ARM_ERR: Payment Calendar: string expected",C_result);
	 XL_readStrCellWD(XL_NxChange,NxChange,"NXNONE"," ARM_ERR: notional exchange: string expected",C_result);
	 XL_readStrCellWD(XL_IncludeMaturity,IncludeMaturity,"NO"," ARM_ERR: include maturity: string expected",C_result);
	 XL_readStrCellWD(XL_adjstartdate,adjstartdate,"NO"," ARM_ERR: adjust start date: string expected",C_result);
	 XL_readNumCellWD(XL_Binary,C_Binary,C_Binary_default," ARM_ERR: adjust start date: string expected",C_result);

	ARM_result currencyres;
	ARMLOCAL_ARM_GetDefaultCurrency (currencyres);

	if(C_Currency == "DEFAULT")
	{
		if(currencyres.getRetCode () != ARM_OK)
		{	ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;}
		else
		{	C_Currency = currencyres.getString ();	}
	}

	long retCode;

	if((C_DayCountFrq = ARM_ConvDayCount (DayCountFrq)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	// if((C_AccOnDef = ARM_ConvAccOnDef (AccOnDef, C_result)) == ARM_DEFAULT_ERR)
	// {
	// 	ARM_ARG_ERR();
	// 	return (LPXLOPER)&XL_result;
	// }

	if((C_StubRule = ARM_ConvStubRule (StubRule)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	C_DecompFreq = ARM_ConvDecompFrequency (DecompFreq);

	long nxChange = ARM_NotionalExchange(NxChange);

	if (IncludeMaturity=="Y") C_IncludeMaturity=true;
	if (adjstartdate=="Y") C_adjstartdate=1;

	//long LegType = ARM_ConvLegType ("RUNNING",C_result);
	long LegType = ARM_ConvLegType ("SWAPLEG",C_result);

	long objId;
	CCString prevClass;
	
	CCString curClass = LOCAL_CLN_CLASS;
	CCString stringId = GetLastCurCellEnvValue ();
	
	if(!stringId)
	{	
		retCode = ICMLOCAL_CLN(C_StartDate,
							   C_Maturity,
							   C_refdate,
							   C_fstcpnrefdate,
							   LocalGetNumObjectId(C_IndexObjectId),
							   C_Spread,
							   C_Notional,
							   C_AccOnDef,
							   C_DayCountFrq,
							   C_DecompFreq,
							   C_StubRule,
							   C_resetgap,
							   C_Currency,
							   C_ResetCalendar,
							   C_PayCalendar,
							   nxChange,
							   C_IncludeMaturity,	
							   C_adjstartdate,
							   LegType,
							   C_Binary,
							   C_result);
		

		if(retCode == ARM_OK)
		{
			objId = C_result.getLong ();

			LocalSetCurCellEnvValue (curClass, objId); 

			stringId = LocalMakeObjectId (objId, curClass);
		}
	}
	else
	{
		prevClass = LocalGetStringObjectClass (stringId);
		
		objId = LocalGetNumObjectId (stringId);
			
		if(curClass == prevClass)
		{
		retCode = ICMLOCAL_CLN(C_StartDate,
							   C_Maturity,
							   C_refdate,
							   C_fstcpnrefdate,
							   LocalGetNumObjectId(C_IndexObjectId),
							   C_Spread,
							   C_Notional,
							   C_AccOnDef,
							   C_DayCountFrq,
							   C_DecompFreq,
							   C_StubRule,
							   C_resetgap,
							   C_Currency,
							   C_ResetCalendar,
							   C_PayCalendar,
							   nxChange,
							   C_IncludeMaturity,	
							   C_adjstartdate,
							   LegType,
							   C_Binary,
							   C_result,
							   objId);
		
			if(retCode == ARM_OK)
			{
				LocalSetCurCellEnvValue (curClass, objId); 

				stringId = LocalMakeObjectId (objId, curClass);
			}
		}
		else
		{
			FreeCurCellContent ();

		retCode = ICMLOCAL_CLN(C_StartDate,
							   C_Maturity,
							   C_refdate,
							   C_fstcpnrefdate,
							   LocalGetNumObjectId(C_IndexObjectId),
							   C_Spread,
							   C_Notional,
							   C_AccOnDef,
							   C_DayCountFrq,
							   C_DecompFreq,
							   C_StubRule,
							   C_resetgap,
							   C_Currency,
							   C_ResetCalendar,
							   C_PayCalendar,
							   nxChange,
							   C_IncludeMaturity,	
							   C_adjstartdate,
							   LegType,
							   C_Binary,
							   C_result);

			
			if(retCode == ARM_OK)
			{
				objId = C_result.getLong ();
			
				LocalSetCurCellEnvValue (curClass, objId); 

				stringId = LocalMakeObjectId (objId, curClass);
			}
		}
	}

	if(retCode == ARM_OK)
	{
		FreeCurCellErr ();
		XL_result.xltype = xltypeStr;
		XL_result.val.str = XL_StrC2StrPascal (stringId);
		XL_result.xltype |= xlbitDLLFree;
	}
	else
	{
		ARM_ERR();
	}

//	ARM_END();
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in ARM_Credit_CLN" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}




/*---- End Of File ----*/

// EOF %M% 
