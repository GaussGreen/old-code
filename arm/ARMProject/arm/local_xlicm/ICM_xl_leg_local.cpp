
#pragma warning(disable :4005 4786)

#include <libCCxll\CCxll.h>

#include <ARM\libarm_local\ARM_local_glob.h>
#include <ARM\libarm_local\ARM_local_class.h>

#include <ARM\libicm_local\ICM_local_leg.h>

#include "ARM_local_interface.h"
#include "ARM_local_interglob.h"
#include "ARM_xl_trycatch_local.h"
#include <util\fromto.h>
#include "ExcelTools.h"


__declspec(dllexport) LPXLOPER WINAPI ARM_Credit_FIXEDLEG (LPXLOPER XL_startDate,
														   LPXLOPER XL_endDate,
														   LPXLOPER XL_fixedRate,
														   LPXLOPER XL_AccOnDefault,
														   LPXLOPER XL_AccDayCount,
														   LPXLOPER XL_LastIndexFixing,
														   LPXLOPER	XL_rcvOrPay,
														   LPXLOPER	XL_freq,
														   LPXLOPER	XL_dayCount,
														   LPXLOPER	XL_decompFreq,
														   LPXLOPER	XL_payTiming,
														   LPXLOPER	XL_intRule,
														   LPXLOPER	XL_stubRule,
														   LPXLOPER	XL_discountCcy,
														   LPXLOPER	XL_payCalName,
														   LPXLOPER	XL_nxChange,
														   LPXLOPER	XL_refDate)
{
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();

	// C variable
	double C_startDate;
	double C_endDate;
	
	CCString C_receiveOrPay;
	double receiveOrPayId;

	double   C_fixedRate_double;
	CCString C_fixedRate_str;
	long     fixedRateType;
	
	CCString C_dayCount;
	long dayCountId;

	CCString C_freq;
	long freqId;

	CCString C_decompFreq;
	long decompFreqId;

	CCString C_payTiming;
	long payTimingId;
	
	CCString C_intRule;
	long intRuleId;
	
	CCString C_stubRule;
	long stubRuleId;

	CCString C_ccy;
	long ccyId;
			
	CCString C_payCal;

    CCString C_nxChange;
    long     nxChange;

	double C_refDate;
	double C_refDate_default = -1.0;

	// Ajouts Credit ----------------------------------------------

	// CCString AccruedOnDefault;
	// int		lAccruedOnDefault;	

	CCString AccruedDayCount;
	int		AccruedDayCountId;

	double	 LastIndexFixing;
	double	 LastIndexFixing_default = 0.;

	// error
	static int error;
	static char* reason = "";

	XL_readNumCell(XL_startDate,C_startDate," ARM_ERR: start date: date expected",C_result);
	XL_readNumCell(XL_endDate,C_endDate," ARM_ERR: end date: date expected",C_result);
	XL_readStrCell(XL_rcvOrPay,C_receiveOrPay," ARM_ERR: receive or pay: string expected",C_result);
	XL_readStrOrNumCell(XL_fixedRate, C_fixedRate_str, C_fixedRate_double, fixedRateType,
		   " ARM_ERR: fixed rate: numeric or object ID string expected",C_result);
	// Ajouts Credit ----------------------------------------------

	// XL_readStrCellWD(XL_AccOnDefault, AccruedOnDefault, "ACC"," ARM_ERR: AccruedOnDefault: possibles values are ACC, CUR, NOS",C_result);

	
	qPAYMENT_PREMIUM_LEG lAccruedOnDefault; 
	ExcelTools::econvert(XL_AccOnDefault,"ACC",lAccruedOnDefault); 

	XL_readStrCellWD(XL_AccDayCount,AccruedDayCount,"30/360"," ARM_ERR: Accrued day count: string expected",C_result);
	XL_readNumCellWD(XL_LastIndexFixing,LastIndexFixing,LastIndexFixing_default," ARM_ERR: LastIndexFixing: numeric expected",C_result);

	// Fin Ajouts Credit ------------------------------------------
	XL_readStrCellWD(XL_dayCount,C_dayCount,"30/360"," ARM_ERR: day count: string expected",C_result);
	XL_readStrCellWD(XL_freq,C_freq,"1"," ARM_ERR: frequency: string expected",C_result);
	XL_readStrCellWD(XL_decompFreq,C_decompFreq,"P"," ARM_ERR: decomp frequency: string expected",C_result);
	XL_readStrCellWD(XL_payTiming,C_payTiming,"ARR"," ARM_ERR: pay timing: string expected",C_result);
	XL_readStrCellWD(XL_intRule,C_intRule,"ADJ"," ARM_ERR: interpolation rule: string expected",C_result);
	XL_readStrCellWD(XL_stubRule,C_stubRule,"SS"," ARM_ERR: stub rule: string expected",C_result);
	XL_readStrCellWD(XL_discountCcy,C_ccy,"DEFAULT"," ARM_ERR: currency id: object expected",C_result);
	XL_readStrCellWD(XL_payCalName,C_payCal,"NULL"," ARM_ERR: payment calendar name: string expected",C_result);
	XL_readStrCellWD(XL_nxChange,C_nxChange,"NXNONE"," ARM_ERR: Notional xchge: string expected",C_result);
	XL_readNumCellWD(XL_refDate,C_refDate,C_refDate_default," ARM_ERR: reference date: date expected",C_result);

	if(C_ccy == "DEFAULT")
	{
		ccyId = ARM_NULL_OBJECT;
	}
	else
	{
		ccyId = LocalGetNumObjectId (C_ccy);
	}

	if((receiveOrPayId = ARM_ConvRecOrPay (C_receiveOrPay, C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	if((dayCountId = ARM_ConvDayCount (C_dayCount)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	if((freqId = ARM_ConvFrequency (C_freq, C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	if((decompFreqId = ARM_ConvDecompFrequency (C_decompFreq)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	if((payTimingId = ARM_ConvPayResetRule (C_payTiming)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	if((intRuleId = ARM_ConvIntRule (C_intRule)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	if((stubRuleId = ARM_ConvStubRule (C_stubRule)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	if ( fixedRateType == XL_TYPE_STRING )
	{
		C_fixedRate_double = (double) LocalGetNumObjectId(C_fixedRate_str);

		fixedRateType = 1L;
	}
	else
	{
		fixedRateType = 0L;
	}

	if ((nxChange = ARM_NotionalExchange(C_nxChange)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	// Ajouts Credit -------------------------------------------------

	/**if((lAccruedOnDefault = ARM_ConvAccOnDef (AccruedOnDefault, C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}**/ 

	if((AccruedDayCountId = ARM_ConvDayCount (AccruedDayCount)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}
	
	// Fin Ajouts Credit ----------------------------------------------


	long retCode;
	long objId;
	CCString prevClass;
	
	CCString curClass = LOCAL_ICM_LEG_CLASS;
	CCString stringId = GetLastCurCellEnvValue ();
	
	if(!stringId)
	{
		retCode = ICMLOCAL_FIXEDLEG (C_startDate,
									 C_endDate,
									 C_fixedRate_double,
									 lAccruedOnDefault,
									 AccruedDayCountId,
									 LastIndexFixing,
									 receiveOrPayId,
									 freqId,
									 dayCountId,
									 decompFreqId,
									 payTimingId,
									 intRuleId,
									 stubRuleId,
									 ccyId,
									 C_payCal,
									 nxChange,
									 C_refDate,
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
		retCode = ICMLOCAL_FIXEDLEG (C_startDate,
									 C_endDate,
									 C_fixedRate_double,
									 lAccruedOnDefault,
									 AccruedDayCountId,
									 LastIndexFixing,
									 receiveOrPayId,
									 freqId,
									 dayCountId,
									 decompFreqId,
									 payTimingId,
									 intRuleId,
									 stubRuleId,
									 ccyId,
									 C_payCal,
									 //fixedRateType,
									 nxChange,
									 C_refDate,
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

		retCode = ICMLOCAL_FIXEDLEG (C_startDate,
									 C_endDate,
									 C_fixedRate_double,
									 lAccruedOnDefault,
									 AccruedDayCountId,
									 LastIndexFixing,
									 receiveOrPayId,
									 freqId,
									 dayCountId,
									 decompFreqId,
									 payTimingId,
									 intRuleId,
									 stubRuleId,
									 ccyId,
									 C_payCal,
									 //fixedRateType,
									 nxChange,
									 C_refDate,
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in ARM_Credit_FIXEDLEG" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
	
	
}


__declspec(dllexport) LPXLOPER WINAPI ARM_Credit_GenLeg (LPXLOPER XL_StartDate,
														 LPXLOPER XL_EndDate,
														 LPXLOPER XL_FixedRate,
														 LPXLOPER XL_FixedNotional,
														 LPXLOPER XL_VariableNotional,
														 LPXLOPER XL_VariableRate,
														 LPXLOPER XL_ExchangeNotional,
														 LPXLOPER XL_frequency,
														 LPXLOPER XL_dayCount,
														 LPXLOPER XL_payTiming,
														 LPXLOPER XL_intRule,
														 LPXLOPER XL_stubRule,
														 LPXLOPER XL_discountCcy,
														 LPXLOPER XL_payCalName,
														 LPXLOPER XL_refDate,
														 LPXLOPER XL_IncludeMaturity,
														 LPXLOPER XL_AdjStartDate,
														 LPXLOPER XL_LegType,
														 LPXLOPER XL_IndexId,
														 LPXLOPER XL_CreditLag ,
														 /*LPXLOPER XL_Binary,
														 LPXLOPER XL_Name,
														 LPXLOPER XL_Nxchange,*/
														 LPXLOPER XL_AccOnDefault)
{
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	
/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{

	ARM_NOCALCIFWIZ();
	// C variable
	double C_StartDate;
	double C_EndDate;
	
	double   C_FixedRate=0.;
	double   C_FixedRate_default = 0.;

	double   C_FixedNotional=0.;
	double   C_FixedNotional_default = 10000000.;
	
	CCString C_VariableNotional;
	CCString C_VariableRate;
	CCString C_ExchangeNotional;

	CCString C_frequency;
	long l_frequency;

	CCString C_dayCount;
	long l_dayCountId;

	CCString C_payTiming;
	long l_payTiming;
	
	CCString C_intRule;
	long l_intRuleId;
	
	CCString C_stubRule;
	long l_stubRuleId;

	CCString C_ccy;
	long l_ccyId;

	CCString C_IncludeMaturity;
	bool l_IncludeMaturity = false;

	CCString C_AdjStartDate;
	int l_AdjStartDate = 0;

	CCString C_LegType;
	int l_LegType = 0;

	CCString C_payCalName;

	double C_refDate;
	double C_refDate_default = -1.0;

	CCString C_IndexId;

	double   C_Binary=0.;
	double   C_Binary_default = -999.;

	CCString C_Name;
	CCString C_Name_default = "UNDEF";
	CCString C_Nxchange;
	int l_Nxchange = 0;
	int l_Nxchange_default =0;

	// error
	static int error;
	static char* reason = "";

	double   C_CreditLag=0.;
	double   C_CreditLag_default = 0.;

	// CCString AccruedOnDefault;
	// int		lAccruedOnDefault;

	XL_readNumCell(XL_StartDate,C_StartDate," ARM_ERR: start date: date expected",C_result);
	XL_readNumCell(XL_EndDate,C_EndDate," ARM_ERR: end date: date expected",C_result);
	XL_readNumCellWD(XL_FixedRate,C_FixedRate, C_FixedRate_default, " ARM_ERR: Fixed rate: numeric expected",C_result);
	XL_readNumCellWD(XL_FixedNotional,C_FixedNotional,C_FixedNotional_default, " ARM_ERR: Fixed notional: numeric expected",C_result);
	XL_readStrCellWD(XL_VariableNotional, C_VariableNotional,"NONE"," ARM_ERR: Var Notional: object ARM_ReferenceValue expected",C_result);
	XL_readStrCellWD(XL_VariableRate,C_VariableRate,"NONE"," ARM_ERR: Var Rate: object ARM_ReferenceValue expected",C_result);
	XL_readStrCellWD(XL_ExchangeNotional,C_ExchangeNotional,"NONE"," ARM_ERR: Exchange Notional: object ARM_ReferenceValue expected",C_result);
	XL_readStrCellWD(XL_frequency,C_frequency,"Q"," ARM_ERR: frequency: string expected",C_result);
	XL_readStrCellWD(XL_dayCount,C_dayCount,"A360"," ARM_ERR: day count: string expected",C_result);
	XL_readStrCellWD(XL_payTiming,C_payTiming,"ARR"," ARM_ERR: pay timing: string expected",C_result);
	XL_readStrCellWD(XL_intRule,C_intRule,"ADJ"," ARM_ERR: interpolation rule: string expected",C_result);
	XL_readStrCellWD(XL_stubRule,C_stubRule,"SS"," ARM_ERR: stub rule: string expected",C_result);
	XL_readStrCellWD(XL_discountCcy,C_ccy,"DEFAULT"," ARM_ERR: currency id: object expected",C_result);
	XL_readStrCellWD(XL_payCalName,C_payCalName,"EUR"," ARM_ERR: payment calendar name: string expected",C_result);
	XL_readNumCellWD(XL_refDate,C_refDate,C_refDate_default," ARM_ERR: reference date: date expected",C_result);
	XL_readStrCellWD(XL_IncludeMaturity,C_IncludeMaturity,"Y"," ARM_ERR: payment calendar name: string expected",C_result);
	XL_readStrCellWD(XL_AdjStartDate,C_AdjStartDate,"Y"," ARM_ERR: payment calendar name: string expected",C_result);
	XL_readStrCellWD(XL_LegType,C_LegType,"RUNNING"," ARM_ERR: payment calendar name: string expected",C_result);
	XL_readStrCellWD(XL_IndexId, C_IndexId," ","ARM_ERR: Index : Object Index expected",C_result);
	XL_readNumCellWD(XL_CreditLag,C_CreditLag,C_CreditLag_default," ARM_ERR: creditlag: numeric expected",C_result);
	/*XL_readNumCellWD(XL_Binary,C_Binary, C_Binary_default, " ARM_ERR: Binary rate: numeric expected",C_result);
	XL_readStrCellWD(XL_Name,C_Name,"UNDEF"," ARM_ERR: name: string expected",C_result);
	XL_readStrCellWD(XL_Nxchange,C_Nxchange,"NXNONE"," ARM_ERR: name: string expected",C_result);*/

	// XL_readStrCellWD(XL_AccOnDefault, AccruedOnDefault, "ACC"," ARM_ERR: AccruedOnDefault: possibles values are ACC, CUR, NOS",C_result);

	qPAYMENT_PREMIUM_LEG lAccruedOnDefault ;
	ExcelTools::econvert(XL_AccOnDefault,"ACC",lAccruedOnDefault) ;

	
	if (C_IncludeMaturity == "Y") l_IncludeMaturity = true;
	if (C_AdjStartDate == "Y") l_AdjStartDate = 1;

	if(C_ccy == "DEFAULT")
	{
		l_ccyId = ARM_NULL_OBJECT;
	}
	else
	{
		l_ccyId = LocalGetNumObjectId (C_ccy);
	}

	if((l_dayCountId = ARM_ConvDayCount (C_dayCount)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	if((l_frequency = ARM_ConvFrequency (C_frequency, C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	if((l_payTiming = ARM_ConvPayResetRule (C_payTiming)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	if((l_intRuleId = ARM_ConvIntRule (C_intRule)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	if((l_stubRuleId = ARM_ConvStubRule (C_stubRule)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	if((l_LegType = ARM_ConvLegType(C_LegType, C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	/*if((l_Nxchange = ARM_NotionalType(C_Nxchange, C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}*/
	/** if((lAccruedOnDefault = ARM_ConvAccOnDef (AccruedOnDefault, C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}**/ 

	long retCode;
	long objId;
	CCString prevClass;
	
	CCString curClass = LOCAL_ICM_LEG_CLASS;
	CCString stringId = GetLastCurCellEnvValue ();
	
	if(!stringId)
	{
		retCode = ICMLOCAL_GenLeg (C_StartDate,
										C_EndDate,
										C_FixedRate,
										C_FixedNotional,
										LocalGetNumObjectId(C_VariableNotional),
										LocalGetNumObjectId(C_VariableRate),
										LocalGetNumObjectId(C_ExchangeNotional),
										l_frequency,
										l_dayCountId,
										l_payTiming,
										l_intRuleId,
										l_stubRuleId,
										l_ccyId,
										C_payCalName,
										C_refDate,
										l_IncludeMaturity,
										l_AdjStartDate,
										l_LegType,
										LocalGetNumObjectId(C_IndexId),
										(int)C_CreditLag,
										C_Binary_default,
										C_Name_default,
										l_Nxchange_default,
										lAccruedOnDefault,
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
		retCode = ICMLOCAL_GenLeg(C_StartDate,
										C_EndDate,
										C_FixedRate,
										C_FixedNotional,
										LocalGetNumObjectId(C_VariableNotional),
										LocalGetNumObjectId(C_VariableRate),
										LocalGetNumObjectId(C_ExchangeNotional),
										l_frequency,
										l_dayCountId,
										l_payTiming,
										l_intRuleId,
										l_stubRuleId,
										l_ccyId,
										C_payCalName,
										C_refDate,
										l_IncludeMaturity,
										l_AdjStartDate,
										l_LegType,	
										LocalGetNumObjectId(C_IndexId),
										(int)C_CreditLag,
										C_Binary_default,
										C_Name_default,
										l_Nxchange_default,
										lAccruedOnDefault,
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

		retCode = ICMLOCAL_GenLeg(C_StartDate,
										C_EndDate,
										C_FixedRate,
										C_FixedNotional,
										LocalGetNumObjectId(C_VariableNotional),
										LocalGetNumObjectId(C_VariableRate),
										LocalGetNumObjectId(C_ExchangeNotional),
										l_frequency,
										l_dayCountId,
										l_payTiming,
										l_intRuleId,
										l_stubRuleId,
										l_ccyId,
										C_payCalName,
										C_refDate,
										l_IncludeMaturity,
										l_AdjStartDate,
										l_LegType,
										LocalGetNumObjectId(C_IndexId),
										(int)C_CreditLag,
										C_Binary_default,
										C_Name_default,
										l_Nxchange_default,
										lAccruedOnDefault,
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in ARM_Credit_FIXEDLEGGEN" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
	
}


__declspec(dllexport) LPXLOPER WINAPI ARM_Credit_SetVariableSpread (LPXLOPER XL_secId,
																	LPXLOPER XL_refValId)
{
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	try {

	ARM_NOCALCIFWIZ();

	// C variable
	CCString C_secId;
	CCString C_refValId;

	// error
	static int error;
	static char* reason = "";

	XL_readStrCell(XL_secId,C_secId," ARM_ERR: security id: object expected",C_result);
	XL_readStrCell(XL_refValId,C_refValId," ARM_ERR: reference value id: object expected",C_result);

	long retCode = ICMLOCAL_SetVariableSpread (LocalGetNumObjectId (C_secId),
											LocalGetNumObjectId (C_refValId),
											C_result);

	if ( retCode == ARM_OK )
	{
		FreeCurCellErr();
	
        XL_result.xltype = xltypeStr;
		XL_result.val.str = XL_StrC2StrPascal(C_secId);
		XL_result.xltype |= xlbitDLLFree;
	}
	else
	{
		ARM_ERR();
	}

//	ARM_END();
	}
	catch (Exception&e)
	{
		ExcelCaller::get().setError(e.GetErrorString()); 
		ExcelTools::convert("ARM_ERR",&XL_result) ;
	}
	catch (std::exception&e)
	{
		ExcelCaller::get().setError(e.what()); 
		ExcelTools::convert("ARM_ERR",&XL_result) ;
	}
	catch(...)
	{
		ExcelCaller::get().setError("Unknown Exception");
		ExcelTools::convert("ARM_ERR",&XL_result) ;
	}
	
	return (LPXLOPER)&XL_result;
}


__declspec(dllexport) LPXLOPER WINAPI ARM_Credit_SetLeg (LPXLOPER XL_SecId,
														 LPXLOPER XL_LegId,
														 LPXLOPER XL_LegType)
{
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	try {
	ARM_NOCALCIFWIZ();

	// C variable
	CCString C_SecId;
	CCString C_LegId;
	CCString C_LegType;

	int l_option = 0;

	// error
	static int error;
	static char* reason = "";

	XL_readStrCell(XL_SecId,C_SecId," ARM_ERR: security id: object expected",C_result);
	XL_readStrCell(XL_LegId,C_LegId," ARM_ERR: leg id: object expected",C_result);
	XL_readStrCellWD(XL_LegType,C_LegType,"FEELEG"," ARM_ERR: leg Type: string expected ",C_result);

	if((l_option = ICM_ConvCptType (C_LegType, C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	long retCode = ICMLOCAL_SetLeg (LocalGetNumObjectId (C_SecId),
									LocalGetNumObjectId (C_LegId),
									l_option,	
									C_result);

	if ( retCode == ARM_OK )
	{
		FreeCurCellErr();
	
        XL_result.xltype = xltypeStr;
		XL_result.val.str = XL_StrC2StrPascal("Done");
		XL_result.xltype |= xlbitDLLFree;
	}
	else
	{
		ARM_ERR();
	}
	}
	catch (Exception&e)
	{
		ExcelCaller::get().setError(e.GetErrorString()); 
		ExcelTools::convert("ARM_ERR",&XL_result) ;
	}
	catch (std::exception&e)
	{
		ExcelCaller::get().setError(e.what()); 
		ExcelTools::convert("ARM_ERR",&XL_result) ;
	}
	catch(...)
	{
		ExcelCaller::get().setError("Unknown Exception");
		ExcelTools::convert("ARM_ERR",&XL_result) ;
	}

//	ARM_END();
	
	return (LPXLOPER)&XL_result;
}

// Old Credit Index
LPXLOPER WINAPI ARM_Credit_IndexInt(LPXLOPER XL_IndexName,	
														LPXLOPER XL_labels,
														LPXLOPER XL_YearTerm,
														LPXLOPER XL_Spread,
														LPXLOPER XL_Method,
														LPXLOPER XL_Basis,
														LPXLOPER XL_ResetFrequency,
														LPXLOPER XL_PayFrequency,
														LPXLOPER XL_Currency,
														LPXLOPER XL_DefaultCurveId,
														LPXLOPER XL_FwdRule,
														LPXLOPER XL_ResetTiming,
														LPXLOPER XL_ResetGap,
														LPXLOPER XL_PayTiming,
														LPXLOPER XL_PayGap,
														LPXLOPER XL_IntRule,
														LPXLOPER XL_AdjCal,
														LPXLOPER XL_CM_resetWeekDay,
														LPXLOPER XL_CM_resetOccur)

{
	static XLOPER XL_result;
	ARM_result C_result;
	long prevId = ExcelCaller::get().getObjectId();
	long retCode = 0;
	
	// parameters
	// vector<string>* pvLabels = NULL;
	
	try {
		ARM_NOCALCIFWIZ();

		string IndexName;					ExcelTools::convert(XL_IndexName,IndexName);
		vector<string> Labels;			
		vector<string> vLabels_default ;	ExcelTools::convert(XL_labels,vLabels_default, Labels);
		// if(vLabels.size() > 1)
		// vector<string> pvLabels ; 
		ARM_Vector C_YearTerm;
		ARM_Vector C_YearTerm_default(1); C_YearTerm_default[0]= 5.;	
		
											ExcelTools::convert(XL_YearTerm,C_YearTerm_default, C_YearTerm);
		ARM_Date aDate;

		try {
			for (int j = 0; j < C_YearTerm.size(); j++)
			{
				ARM_Date aDate;
				Local_XLDATE2ARMDATE(C_YearTerm.Elt(j), aDate);
				C_YearTerm.Elt(j) = aDate.GetJulian();
			}
		}catch(...) {
			// nothing
		}
		
		ARM_Vector C_Spread;
		ARM_Vector C_Spread_default(1); C_Spread_default[0] = -999.;	
											ExcelTools::convert(XL_Spread,C_Spread_default, C_Spread);
		string C_Method ;
		long l_Method = 0;					ExcelTools::convert(XL_Method,"AVG", C_Method);
		int Basis =0;
		string C_Basis;						ExcelTools::convert(XL_Basis,"A360", C_Basis);
		string C_ResetFrequency;
		int ResetFrequency = 0;				ExcelTools::convert(XL_ResetFrequency,"Q",C_ResetFrequency);
		string C_PayFrequency;
		int PayFrequency = 0;				ExcelTools::convert(XL_PayFrequency,"Q",C_PayFrequency);
		string Currency;					ExcelTools::convert(XL_Currency,"DEFAULT",Currency);
		string C_DefaultCurveId;			ExcelTools::convert(XL_DefaultCurveId,"NONE",C_DefaultCurveId);
		string C_fwdRule;
		long fwdRuleId = 0;					ExcelTools::convert(XL_FwdRule,"MF",C_fwdRule);
		string C_resetTiming;
		long resetTimingId;					ExcelTools::convert(XL_ResetTiming,"ADV",C_resetTiming);
		int C_resetGap = 0.;
		int C_resetGap_default = 0.;		ExcelTools::convert(XL_ResetGap,C_resetGap_default,C_resetGap);
		string C_payTiming;
		long payTimingId = 0;				ExcelTools::convert(XL_PayTiming,"ARR",C_payTiming);
		int C_payGap = 0;
		int C_payGap_default = 0;		ExcelTools::convert(XL_PayGap,C_payGap_default,C_payGap);
		string C_intRule;
		long intRuleId = 0;					ExcelTools::convert(XL_IntRule,"ADJ",C_intRule);
		qCDS_ADJ l_AdjCal = qCredit_Adjust20;
		string C_AdjCal;					ExcelTools::convert(XL_AdjCal,"STDCDS",C_AdjCal);
		double cm_resetWeekDay = 0.,cm_resetWeek_default=5 ;	// vendredi
		ExcelTools::convert(XL_CM_resetWeekDay,cm_resetWeek_default,cm_resetWeekDay);
		double cm_resetOccur = 0.,cm_resetOccur_default=2;	// 2eme vendredi
		ExcelTools::convert(XL_CM_resetOccur,cm_resetOccur_default,cm_resetOccur );
		// conversion
		if((fwdRuleId = ARM_ConvFwdRule (C_fwdRule.c_str(), C_result)) == ARM_DEFAULT_ERR){ ARM_ARG_ERR();return (LPXLOPER)&XL_result; }
		resetTimingId = ARM_ConvPayResetRule (C_resetTiming.c_str());
		payTimingId = ARM_ConvPayResetRule (C_payTiming.c_str());
		intRuleId = ARM_ConvIntRule (C_intRule.c_str());
		//Ccy
		ARM_result currencyres;
		ARMLOCAL_ARM_GetDefaultCurrency (currencyres);
		if(Currency == "DEFAULT")
		{
			if(currencyres.getRetCode () != ARM_OK){ARM_ARG_ERR();return (LPXLOPER)&XL_result;}
			else
				Currency = currencyres.getString ();
		}

		if((l_AdjCal = ARM_ConvAdjCalCDS (C_AdjCal.c_str(), C_result)) == ARM_DEFAULT_ERR){	ARM_ARG_ERR();	return (LPXLOPER)&XL_result;}

		if((ResetFrequency = ARM_ConvFrequency (C_ResetFrequency.c_str(), C_result)) == ARM_DEFAULT_ERR){	ARM_ARG_ERR();	return (LPXLOPER)&XL_result;}

		if((PayFrequency = ARM_ConvFrequency (C_PayFrequency.c_str(), C_result)) == ARM_DEFAULT_ERR){ARM_ARG_ERR();	return (LPXLOPER)&XL_result;}

		if((Basis = ARM_ConvDayCount (C_Basis.c_str())) == ARM_DEFAULT_ERR)	{	ARM_ARG_ERR();	return (LPXLOPER)&XL_result;}

		if((l_Method = ARM_ConvIndexMethod(C_Method.c_str(), C_result)) == ARM_DEFAULT_ERR){	ARM_ARG_ERR();	return (LPXLOPER)&XL_result;}

		
		long retCode =  ICMLOCAL_Index( IndexName,
									Labels,
									Basis,
									ResetFrequency,
									PayFrequency,
									C_YearTerm,
									C_Spread,								
									Currency,
									l_Method,
									(int)LocalGetNumObjectId((CCString)C_DefaultCurveId.c_str()),
									(int)fwdRuleId,
									(int)resetTimingId,
									C_resetGap,
									(int)payTimingId,
									C_payGap,
									(int)intRuleId,
									l_AdjCal,
									(int)cm_resetWeekDay,
									(int)cm_resetOccur);
		
		
		 string objectLabel = ExcelCaller::get().setObject(retCode, LOCAL_IRINDEX_CLASS);	
		 ExcelTools::convert(objectLabel, &XL_result);

	} 
	catch (Exception& e)
	{
		ExcelCaller::get().setError(e.GetErrorString());
		ExcelTools::convert("ARM_ERR", &XL_result);
	}
	catch (std::exception& e)
	{
		ExcelCaller::get().setError(e.what() );
		ExcelTools::convert("ARM_ERR", &XL_result);
	}
	catch (...)
	{
		ExcelCaller::get().setError("unknown error");
		ExcelTools::convert("ARM_ERR", &XL_result);
	}
    return (LPXLOPER)&XL_result;
}
// oteh parameters in default case of index
__declspec(dllexport) LPXLOPER WINAPI ARM_Credit_Index(LPXLOPER XL_labels,
														LPXLOPER XL_YearTerm,
														LPXLOPER XL_Spread,
														LPXLOPER XL_Method,
														LPXLOPER XL_Basis,
														LPXLOPER XL_ResetFrequency,
														LPXLOPER XL_PayFrequency,
														LPXLOPER XL_Currency,
														LPXLOPER XL_DefaultCurveId,
														LPXLOPER XL_FwdRule,
														LPXLOPER XL_ResetTiming,
														LPXLOPER XL_ResetGap,
														LPXLOPER XL_PayTiming,
														LPXLOPER XL_PayGap,
														LPXLOPER XL_IntRule,
														LPXLOPER XL_AdjCal,
														LPXLOPER XL_CM_resetWeekDay,
														LPXLOPER XL_CM_resetOccur)												
{
														  
	return	ARM_Credit_IndexInt(XL_labels,	
						XL_labels,
						XL_YearTerm,
						XL_Spread,
						XL_Method,
						XL_Basis,
						XL_ResetFrequency,
						XL_PayFrequency,
						XL_Currency,
						XL_DefaultCurveId,
						XL_FwdRule,
						XL_ResetTiming,
						XL_ResetGap,
						XL_PayTiming,
						XL_PayGap,
						XL_IntRule,
						XL_AdjCal,
						XL_CM_resetWeekDay,
						XL_CM_resetOccur);
	
}


__declspec(dllexport) LPXLOPER WINAPI ARM_Credit_IndexCompo(LPXLOPER XL_IndexName,	
														LPXLOPER XL_labels,
														LPXLOPER XL_YearTerm,
														LPXLOPER XL_Spread,
														LPXLOPER XL_Method,
														LPXLOPER XL_Basis,
														LPXLOPER XL_ResetFrequency,
														LPXLOPER XL_PayFrequency,
														LPXLOPER XL_Currency,
														LPXLOPER XL_FwdRule,
														LPXLOPER XL_ResetTiming,
														LPXLOPER XL_ResetGap,
														LPXLOPER XL_PayTiming,
														LPXLOPER XL_PayGap,
														LPXLOPER XL_IntRule,
														LPXLOPER XL_AdjCal,
														LPXLOPER XL_CM_resetWeekDay,
														LPXLOPER XL_CM_resetOccur)

{
	return	ARM_Credit_IndexInt(XL_IndexName,	
						XL_labels,
						XL_YearTerm,
						XL_Spread,
						XL_Method,
						XL_Basis,
						XL_ResetFrequency,
						XL_PayFrequency,
						XL_Currency,
						XL_IndexName, // LPXLOPER vide !!?
						XL_FwdRule,
						XL_ResetTiming,
						XL_ResetGap,
						XL_PayTiming,
						XL_PayGap,
						XL_IntRule,
						XL_AdjCal,
						XL_CM_resetWeekDay,
						XL_CM_resetOccur);
}



__declspec(dllexport) LPXLOPER WINAPI	ARM_Credit_CorridorLeg(
														LPXLOPER XL_Name,
														LPXLOPER XL_startDate,
														LPXLOPER XL_endDate,
														LPXLOPER XL_refdate,
														LPXLOPER XL_fstcpneffdate,
														LPXLOPER XL_RefValueSpreads,
									                    LPXLOPER XL_floatingIdx,
														LPXLOPER XL_leverageFloatIdx,
														LPXLOPER XL_creditIdx,
														LPXLOPER XL_refvalueKUP,
														LPXLOPER XL_refvalueKDW,
														LPXLOPER XL_accondef,
														LPXLOPER XL_accdaycount,
														LPXLOPER XL_payfreq,
														//LPXLOPER XL_resetfreq,
														LPXLOPER XL_daycount,
														LPXLOPER XL_paytiming,
														LPXLOPER XL_intrule,
														LPXLOPER XL_stubrule,
														LPXLOPER XL_disc_ccy,
														LPXLOPER XL_paycalname)
{
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
		ARM_NOCALCIFWIZ();

		// C variable
		double C_startDate;
		double C_endDate;

		double C_refdate;
		double C_refdate_default=-1.;

		double C_fstcpneffdate;
		double C_fstcpneffdate_default=-1.;
		
		CCString C_RefValueSpreads;
		long l_RefValueSpreads;

		CCString C_floatingIdx;
		long l_floatingIdx;

		double C_leverageFloatIdx;
		double C_leverageFloatIdx_default=1.;

		CCString C_creditIdx;
		long l_creditIdx;

		CCString C_refvalueKUP;
		long l_refvalueKUP;

		CCString C_refvalueKDW;
		long l_refvalueKDW;

		// CCString C_accondef;
		// CCString C_accondef_default="ACC";
		// long l_accondef;

		double d_notional = 0.;

		CCString C_accdaycount;
		CCString C_accdaycount_default="A360";
		long l_accdaycount;

		CCString C_payfreq;
		CCString C_payfreq_default="Q";
		long l_payfreq;

		CCString C_resetfreq;
		CCString C_resetfreq_default="M";
		C_resetfreq = C_resetfreq_default;
		long l_resetfreq;

		CCString C_daycount;
		CCString C_daycount_default="A360";
		long l_daycount;

		CCString C_paytiming;
		CCString C_paytiming_default="ARR";
		long l_paytiming;

		CCString C_intrule;
		CCString C_intrule_default="ADJ";
		long l_intrule;

		CCString C_stubrule;
		CCString C_stubrule_default="LS";
		long l_stubrule;

		CCString C_disc_ccy;
		CCString C_disc_ccy_default="DEFAULT";

		CCString C_paycalname;
		CCString C_paycalname_default="DEFAULT";

		CCString C_Name = "";

		// error
		static int error;
		static char* reason = "";
		XL_readStrCell(XL_Name, C_Name, " ARM_ERR: name expected", C_result);
		XL_readNumCell(XL_startDate,C_startDate," ARM_ERR: start date: date expected",C_result);
		XL_readNumCell(XL_endDate,C_endDate," ARM_ERR: end date: date expected",C_result);
		XL_readNumCellWD(XL_refdate,C_refdate,C_refdate_default," ARM_ERR: reference date: date expected",C_result);
		XL_readNumCellWD(XL_fstcpneffdate,C_fstcpneffdate,C_fstcpneffdate_default," ARM_ERR: effective date: date expected",C_result);
		XL_readStrCell(XL_RefValueSpreads,C_RefValueSpreads," ARM_ERR: spreads: object expected",C_result);    
		XL_readStrCell(XL_floatingIdx,C_floatingIdx," ARM_ERR: float index : object expected",C_result);    
		XL_readNumCellWD(XL_leverageFloatIdx,C_leverageFloatIdx,C_leverageFloatIdx_default," ARM_ERR: leverage: numeric expected",C_result);
		XL_readStrCell(XL_creditIdx,C_creditIdx," ARM_ERR: credit index : object expected",C_result);    
		XL_readStrCell(XL_refvalueKUP,C_refvalueKUP," ARM_ERR: strike up: object expected",C_result);    
		XL_readStrCell(XL_refvalueKDW,C_refvalueKDW," ARM_ERR: strike down: object expected",C_result);   
		// XL_readStrCellWD(XL_accondef,C_accondef,C_accondef_default," ARM_ERR: accrued on def: object expected",C_result);    
		
		qPAYMENT_PREMIUM_LEG l_accondef; 
		ExcelTools::econvert(XL_accondef,"ACC",l_accondef); 
		XL_readStrCellWD(XL_accdaycount,C_accdaycount,C_accdaycount_default," ARM_ERR: accrued daycoubt: string expected",C_result);    
		XL_readStrCellWD(XL_payfreq,C_payfreq,C_payfreq_default," ARM_ERR: payfreq: string expected",C_result);    
		//XL_readStrCellWD(XL_resetfreq,C_resetfreq,C_resetfreq_default," ARM_ERR: resetfreq: string expected",C_result);    
		XL_readStrCellWD(XL_daycount,C_daycount,C_daycount_default," ARM_ERR: daycount: string expected",C_result);    
		XL_readStrCellWD(XL_paytiming,C_paytiming,C_paytiming_default," ARM_ERR: paytiming: string expected",C_result);    
		XL_readStrCellWD(XL_intrule,C_intrule,C_intrule_default," ARM_ERR: intrule: string expected",C_result);    
		XL_readStrCellWD(XL_stubrule,C_stubrule,C_stubrule_default," ARM_ERR: stub: string expected",C_result);    
		XL_readStrCellWD(XL_disc_ccy,C_disc_ccy,C_disc_ccy_default," ARM_ERR: disc ccy: string expected",C_result);    
		XL_readStrCellWD(XL_paycalname,C_paycalname,C_paycalname_default," ARM_ERR: paycal: string expected",C_result);    


		l_RefValueSpreads = LocalGetNumObjectId(C_RefValueSpreads);
		l_floatingIdx = LocalGetNumObjectId(C_floatingIdx);
		l_creditIdx = LocalGetNumObjectId(C_creditIdx);
		l_refvalueKUP = LocalGetNumObjectId(C_refvalueKUP);
		l_refvalueKDW = LocalGetNumObjectId(C_refvalueKDW);

		if ( C_payfreq == "-1" )
		{
		   l_payfreq = K_DEF_FREQ;
		}
		else
		{
			if ((l_payfreq = ARM_ConvFrequency(C_payfreq, C_result)) == ARM_DEFAULT_ERR)
			{
	  			ARM_ARG_ERR();
				return (LPXLOPER)&XL_result;
			}
		}
			
		if ( C_resetfreq == "-1" )
		{
			l_resetfreq = K_DAILY;
		}
		else
		{
			if ((l_resetfreq = ARM_ConvFrequency(C_resetfreq, C_result)) == ARM_DEFAULT_ERR)
			{
	  			ARM_ARG_ERR();
				return (LPXLOPER)&XL_result;
			}
		}

		// if ((l_accondef = ARM_ConvAccOnDef(C_accondef, C_result)) == ARM_DEFAULT_ERR)
		// {
	  	// 	ARM_ARG_ERR();
		// 	return (LPXLOPER)&XL_result;
		// }

		
		if((l_accdaycount = ARM_ConvDayCount (C_accdaycount)) == ARM_DEFAULT_ERR)
		{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
		}

		if((l_daycount = ARM_ConvDayCount (C_daycount)) == ARM_DEFAULT_ERR)
		{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
		}


		l_paytiming = ARM_ConvPayResetRule (C_paytiming);
		l_stubrule = ARM_ConvStubRule (C_stubrule);
		l_intrule = ARM_ConvIntRule (C_intrule);

		if((l_stubrule = ARM_ConvStubRule (C_stubrule)) == ARM_DEFAULT_ERR)
		{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
		}


		long retCode;
		long objId;
		CCString prevClass;
		
		CCString curClass = LOCAL_ICM_CORRIDOR_LEG_CLASS;
		CCString stringId = GetLastCurCellEnvValue ();

		if (!stringId)
		{
			retCode = ICMLOCAL_CORRIDORLEG(C_Name,
								 C_startDate,
								 C_endDate,
								 C_refdate,
								 C_fstcpneffdate,
								 l_RefValueSpreads,
								 l_floatingIdx,
								 C_leverageFloatIdx,
								 l_creditIdx,
								 l_refvalueKUP,
								 l_refvalueKDW,
								 l_accondef,
								 l_accdaycount,
								 l_payfreq,
								 l_resetfreq,
								 l_daycount,
								 l_paytiming,
								 l_intrule,
								 l_stubrule,
								 C_disc_ccy,
								 C_paycalname,
								 C_result);

			if ( retCode == ARM_OK )
			{
				objId = C_result.getLong ();

				LocalSetCurCellEnvValue(curClass, objId); 

				stringId = LocalMakeObjectId (objId, curClass);
			}
		}
		else
		{
			prevClass = LocalGetStringObjectClass (stringId);
			
			objId = LocalGetNumObjectId (stringId);
				
			if(curClass == prevClass)
			{
				retCode = ICMLOCAL_CORRIDORLEG(
								C_Name,
								C_startDate,
								 C_endDate,
								 C_refdate,
								 C_fstcpneffdate,
								 l_RefValueSpreads,
								 l_floatingIdx,
								 C_leverageFloatIdx,
								 l_creditIdx,
								 l_refvalueKUP,
								 l_refvalueKDW,
								 l_accondef,
								 l_accdaycount,
								 l_payfreq,
								 l_resetfreq,
								 l_daycount,
								 l_paytiming,
								 l_intrule,
								 l_stubrule,
								 C_disc_ccy,
								 C_paycalname,
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
				FreeCurCellContent();

				retCode = ICMLOCAL_CORRIDORLEG(C_Name, C_startDate,
								 C_endDate,
								 C_refdate,
								 C_fstcpneffdate,
								 l_RefValueSpreads,
								 l_floatingIdx,
								 C_leverageFloatIdx,
								 l_creditIdx,
								 l_refvalueKUP,
								 l_refvalueKDW,
								 l_accondef,
								 l_accdaycount,
								 l_payfreq,
								 l_resetfreq,
								 l_daycount,
								 l_paytiming,
								 l_intrule,
								 l_stubrule,
								 C_disc_ccy,
								 C_paycalname,
								 C_result);

				if(retCode == ARM_OK)
				{
					objId = C_result.getLong ();
				
					LocalSetCurCellEnvValue (curClass, objId); 

					stringId = LocalMakeObjectId (objId, curClass);
				}
			}
		}

		if (retCode == ARM_OK)
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in ARM_Credit_CorridorLeg" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}

__declspec(dllexport) LPXLOPER WINAPI ARM_Credit_CorridorLeg_Sche(LPXLOPER XL_Notional,
														LPXLOPER XL_RecieveOrPay,
														LPXLOPER XL_RefValueSpreads,
									                    LPXLOPER XL_floatingIdx,
														LPXLOPER XL_leverageFloatIdx,
														LPXLOPER XL_creditIdx,
														LPXLOPER XL_refvalueKUP,
														LPXLOPER XL_refvalueKDW,
														LPXLOPER XL_ScheduleInfoId,
														LPXLOPER XL_accondef,
														LPXLOPER XL_disc_ccy,
														LPXLOPER XL_Name)
{

	static XLOPER XL_result;
	ARM_result C_result;
	try {
		long prevId = ExcelCaller::get().getObjectId();
		string C_RefValueSpreads = "";
		long l_RefValueSpreads = 0;

		string recieveOrPay="";
		long ConvertRorP = 1;

		string C_floatingIdx="";
		long l_floatingIdx =0;

		string s_leverageFloatIdx="";
		double d_leverageFloatIdx =0;
		double C_leverageFloatIdx_default=1.;

		double notional = 0;

		string C_creditIdx= "";
		long l_creditIdx =0;

		string C_refvalueKUP ="";
		long l_refvalueKUP=0;

		string C_refvalueKDW="";
		long l_refvalueKDW=0;

		string strSchedule_info = "";
		long l_schedule_Info = 0;

		string C_disc_ccy = "";

		string C_Name = "";

		qPAYMENT_PREMIUM_LEG l_accondef; 
		ExcelTools::econvert(XL_accondef,"ACC",l_accondef);

		// error
		static int error;
		static char* reason = "";

		ExcelTools::convert(XL_Notional,notional); 
		if (notional == 0 || notional < 0) {
			C_result.setMsg("notional expected and must be positif");
			ARM_ERR();
		}
		ExcelTools::convert(XL_RecieveOrPay,recieveOrPay);
		if (recieveOrPay.empty()) {
			C_result.setMsg("recieveOrPay expected");
			ARM_ERR();
		}
		if(recieveOrPay == "R" )
			ConvertRorP = 1;
		else if (recieveOrPay == "P")
			ConvertRorP = -1;
		else {
			C_result.setMsg("recieveOrPay  : P or R ");
			ARM_ERR();
		}
		
		ExcelTools::convert(XL_RefValueSpreads,C_RefValueSpreads); 
		if (C_RefValueSpreads.empty()) {
			C_result.setMsg("RefValueSpreads expected");
			ARM_ERR();
		}
		ExcelTools::convert(XL_floatingIdx, C_floatingIdx);
		if (C_floatingIdx.empty()) {
			C_result.setMsg("floatingIdx expected");
			ARM_ERR();
		}
		ExcelTools::convert(XL_leverageFloatIdx,s_leverageFloatIdx);
		if (s_leverageFloatIdx.empty()) {
			d_leverageFloatIdx = C_leverageFloatIdx_default;
		}else {
			ExcelTools::convert(XL_leverageFloatIdx,C_leverageFloatIdx_default,d_leverageFloatIdx);
		}
		ExcelTools::convert(XL_creditIdx,C_creditIdx);
		if (C_creditIdx.empty()) {
			C_result.setMsg("creditIdx expected");
			ARM_ERR();
		}
		ExcelTools::convert(XL_refvalueKUP, C_refvalueKUP);
		if (C_refvalueKUP.empty()) {
			C_result.setMsg("refvalueKUP expected");
			ARM_ERR();
		}
		ExcelTools::convert(XL_refvalueKDW, C_refvalueKDW);
		if (C_refvalueKDW.empty()) {
			C_result.setMsg("refvalueKDW expected");
			ARM_ERR();
		}
		ExcelTools::convert(XL_ScheduleInfoId, strSchedule_info);
		if (strSchedule_info.empty()) {
			C_result.setMsg("Schedule_info object Id is expected");
			ARM_ERR();
		}
		ExcelTools::convert(XL_Name,"", C_Name);
		ExcelTools::convert(XL_disc_ccy, C_disc_ccy);
		if (C_disc_ccy.empty()) C_disc_ccy = "EUR";


		l_RefValueSpreads = LocalGetNumObjectId(CCString(C_RefValueSpreads.c_str()));
		l_floatingIdx = LocalGetNumObjectId(CCString(C_floatingIdx.c_str()));
		l_creditIdx = LocalGetNumObjectId(CCString(C_creditIdx.c_str()));
		l_refvalueKUP = LocalGetNumObjectId(CCString(C_refvalueKUP.c_str()));
		l_refvalueKDW = LocalGetNumObjectId(CCString(C_refvalueKDW.c_str()));
		l_schedule_Info = LocalGetNumObjectId(CCString(strSchedule_info.c_str()));

		long newId = ICMLOCAL_CORRIDORLEG_SCHE(C_Name, 
								 notional,
								 ConvertRorP,
								 l_RefValueSpreads,
								 l_floatingIdx,
								 d_leverageFloatIdx,
								 l_creditIdx,
								 l_refvalueKUP,
								 l_refvalueKDW,
								 l_schedule_Info,
								 l_accondef,
								 C_disc_ccy,
								 C_result);
		string objectLabel = ExcelCaller::get().setObject(newId,LOCAL_ICM_CORRIDOR_LEG_CLASS);
		//ExcelCaller::get().setObject(objectLabel);
		ExcelTools::convert(objectLabel, &XL_result);	
	} 
	catch (Exception& e)
	{
		ExcelCaller::get().setError(e.GetErrorString());
		ExcelTools::convert("ARM_ERR", &XL_result);
	}
	catch (std::exception& e)
	{
		ExcelCaller::get().setError(e.what() );
		ExcelTools::convert("ARM_ERR", &XL_result);
	}
	catch (...)
	{
		ExcelCaller::get().setError("unknown error");
		ExcelTools::convert("ARM_ERR", &XL_result);
	}

    return (LPXLOPER)&XL_result;
}

__declspec(dllexport) LPXLOPER WINAPI ARM_Credit_IRLEGTOCREDITLEG(LPXLOPER XL_LegId,
																  LPXLOPER XL_LegType,
																  LPXLOPER XL_creditindexId,
																  LPXLOPER XL_PricerId)
{
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	try {

	ARM_NOCALCIFWIZ();

	// C variable
	
	CCString C_LegId;
	CCString C_creditindexId;
	CCString C_LegType;
	CCString C_Pricer;

	CCString C_creditindexId_default = "NONE";
	CCString C_LegType_default = "SWAPLEG";
	CCString C_Pricer_default = "NONE";

	long l_LegType; 

	long retCode;

	// error
	static int error;
	static char* reason = "";

	XL_readStrCell(XL_LegId,C_LegId," ARM_ERR: Leg: string expected",C_result);
	XL_readStrCellWD(XL_LegType,C_LegType,C_LegType_default," ARM_ERR: LegType: string expected",C_result);
	XL_readStrCellWD(XL_creditindexId,C_creditindexId,C_creditindexId_default," ARM_ERR: RcvFee: string expected",C_result);
	XL_readStrCellWD(XL_PricerId,C_Pricer,C_Pricer_default," ARM_ERR: Pricer: string expected",C_result);
	
	if((l_LegType = ARM_ConvLegType(C_LegType, C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}
	long objId;
	CCString prevClass;
	
	CCString curClass =  LOCAL_ICM_LEG_CLASS;
	CCString stringId = GetLastCurCellEnvValue ();
	
	if(!stringId)
	{	
		retCode = ICMLOCAL_IRLEGTOCREDITLEG(LocalGetNumObjectId (C_LegId),
											l_LegType,
											LocalGetNumObjectId (C_creditindexId),
											LocalGetNumObjectId (C_Pricer),
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
		retCode = ICMLOCAL_IRLEGTOCREDITLEG(LocalGetNumObjectId (C_LegId),
											l_LegType,
											LocalGetNumObjectId (C_creditindexId),
											LocalGetNumObjectId (C_Pricer),
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

		retCode = ICMLOCAL_IRLEGTOCREDITLEG(LocalGetNumObjectId (C_LegId),
											l_LegType,
											LocalGetNumObjectId (C_creditindexId),
											LocalGetNumObjectId (C_Pricer),
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
	catch (Exception&e)
	{
		ExcelCaller::get().setError(e.GetErrorString()); 
		ExcelTools::convert("ARM_ERR",&XL_result) ;
	}
	catch (std::exception&e)
	{
		ExcelCaller::get().setError(e.what()); 
		ExcelTools::convert("ARM_ERR",&XL_result) ;
	}
	catch(...)
	{
		ExcelCaller::get().setError("Unknown Exception");
		ExcelTools::convert("ARM_ERR",&XL_result) ;
	}
	
	return (LPXLOPER)&XL_result;
}



/*---- End Of File ----*/

// EOF %M% 
