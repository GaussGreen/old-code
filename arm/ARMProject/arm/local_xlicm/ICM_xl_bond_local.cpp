
#pragma warning(disable :4005 4786)

#include <libCCxll\CCxll.h>

#include <ARM\libarm_local\ARM_local_glob.h>
#include <ARM\libarm_local\ARM_local_class.h>

#include <ARM\libicm_local\ICM_local_bond.h>

#include "ARM_local_interface.h"
#include "ARM_local_interglob.h"

#include "ARM_xl_trycatch_local.h"
#include "ExcelTools.h"
#include <util\fromto.h>


__declspec(dllexport) LPXLOPER WINAPI ARM_Credit_BOND ( LPXLOPER XL_StartDate,
														LPXLOPER XL_Maturity,
														LPXLOPER XL_Rate,
														LPXLOPER XL_FixingFreq,
														LPXLOPER XL_DayCountFrq,
														LPXLOPER XL_First_period_refdate,
														LPXLOPER XL_NotAmount,	
														LPXLOPER XL_AccOnDef,	
														LPXLOPER XL_Currency,
														LPXLOPER XL_PayCal,
														LPXLOPER XL_AccruedDayCount,
														LPXLOPER XL_SettlementGap,
														LPXLOPER XL_RedemptionValue)

{


	// return
	static XLOPER XL_result;
	ARM_result C_result;

	// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();

	// C variable
	

	double StartDate;
	double Maturity;
	double Rate;

	CCString C_FixingFreq;

	int FixingFreq;

	int DayCountFrq;

	double First_period_refdate;
	double First_period_refdate_default = -1.;

	double NotAmount;	
	double NotAmount_default = 1000000.;	

	CCString Currency;
	CCString C_AccruedDayCount;
	double AccruedDayCount = 0.;
	
	double SettlementGap;
	double SettlementGap_default = 2.;

	double RedemptionValue;
	double RedemptionValue_default = 100.;

	CCString C_Currency;
	CCString C_DayCountFrq;

	// CCString C_AccOnDef;
	// double AccruedOnDefault = 0.;

	CCString C_AmortOnDef;
	double AmortOnDef = 0.;

	CCString C_IntOnDef;
	double IntOnDef = 0.;

	CCString C_PayCal;

	// error
	static int error;
	static char* reason = "";


	XL_readNumCell(XL_StartDate,StartDate," ARM_ERR: Start date: date expected",C_result);
	XL_readNumCell(XL_Maturity,Maturity," ARM_ERR: Maturity date: date expected",C_result);
	XL_readNumCell(XL_Rate,Rate," ARM_ERR: Rate: numeric expected",C_result);

	XL_readStrCellWD(XL_FixingFreq,C_FixingFreq,"A"," ARM_ERR: FixingFreq: numeric expected",C_result);
	XL_readStrCellWD(XL_DayCountFrq,C_DayCountFrq,"ACTUAL"," ARM_ERR: DayCountFrq: numeric expected",C_result);
	
	XL_readNumCellWD(XL_First_period_refdate,First_period_refdate,First_period_refdate_default," ARM_ERR: First_period_refdate Rate: date expected",C_result);
	XL_readNumCellWD(XL_NotAmount,NotAmount,NotAmount_default," ARM_ERR: NotAmountIn: numeric expected",C_result);

	// XL_readStrCellWD(XL_AccOnDef,C_AccOnDef,"ACC"," ARM_ERR: AccOnDef: string expected",C_result);
	qPAYMENT_PREMIUM_LEG AccruedOnDefault ;
	ExcelTools::econvert(XL_AccOnDef,"ACC",AccruedOnDefault ); 
	

	XL_readStrCellWD(XL_Currency,C_Currency,"DEFAULT"," ARM_ERR: currency: string expected",C_result);
	XL_readStrCellWD(XL_PayCal,C_PayCal,"NULL"," ARM_ERR: Payment Calendar: string expected",C_result);

	XL_readStrCellWD(XL_AccruedDayCount,C_AccruedDayCount,"ACTUAL"," ARM_ERR: AccruedDayCount: numeric expected",C_result);
	XL_readNumCellWD(XL_SettlementGap,SettlementGap,SettlementGap_default," ARM_ERR: SettlementGap: numeric expected",C_result);
	XL_readNumCellWD(XL_RedemptionValue,RedemptionValue,RedemptionValue_default," ARM_ERR: RedemptionValue (%): numeric expected",C_result);


	ARM_result currencyres;
	ARMLOCAL_ARM_GetDefaultCurrency (currencyres);

	if((ARM_CheckDate(StartDate, C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	if((ARM_CheckDate(Maturity, C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

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

	if((FixingFreq = ARM_ConvFrequency (C_FixingFreq, C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	if((DayCountFrq = ARM_ConvDayCount (C_DayCountFrq)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	if((AccruedDayCount = ARM_ConvDayCount (C_AccruedDayCount)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}
	
	// if((AccruedOnDefault = ARM_ConvAccOnDef (C_AccOnDef, C_result)) == ARM_DEFAULT_ERR)
	// {
	// 	ARM_ARG_ERR();
	// 	return (LPXLOPER)&XL_result;
	// }

	if (SettlementGap>30)
	{
		C_result.setMsg("Error: Settlement must be <= 30");
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	if ((RedemptionValue<0.) || (RedemptionValue>1000.))
	{
		C_result.setMsg("Error: Redemption Value must be in [0;1000]");
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	long objId;
	CCString prevClass;
	
	CCString curClass = LOCAL_BONDDEF_CLASS;
	CCString stringId = GetLastCurCellEnvValue ();
	
	if(!stringId)
	{	
		retCode = ICMLOCAL_BOND   (Rate,
								   StartDate,
								   Maturity,
								   FixingFreq,
								   First_period_refdate,
								   NotAmount ,
								    AccruedOnDefault,
								   C_Currency,
								   C_PayCal,
								   (int)DayCountFrq,
								   (int)AccruedDayCount,
								   (int)SettlementGap,
								   RedemptionValue,
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
		retCode = ICMLOCAL_BOND   (Rate,
								   StartDate,
								   Maturity,
								   FixingFreq,
								   First_period_refdate,
								   NotAmount ,
								   AccruedOnDefault,
								   C_Currency,
								   C_PayCal,
								   (int)DayCountFrq,
								   (int)AccruedDayCount,
								   (int)SettlementGap,
								   RedemptionValue,
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

		retCode = ICMLOCAL_BOND   (Rate,
								   StartDate,
								   Maturity,
								   FixingFreq,
								   First_period_refdate,
								   NotAmount ,
								   AccruedOnDefault,
								   C_Currency,
								   C_PayCal,
								   (int)DayCountFrq,
								   (int)AccruedDayCount,
								   (int)SettlementGap,
								   RedemptionValue,
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in ARM_Credit_BOND" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


__declspec(dllexport) LPXLOPER WINAPI ARM_Credit_YTOPRICE (LPXLOPER XL_bondId,
													  LPXLOPER XL_settlement,
													  LPXLOPER XL_yield)
{
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	try {
	ARM_NOCALCIFWIZ();

	// C variable
	CCString C_bondId;
	double C_settlement;
	double C_yield;
	
	// error
	static int error;
	static char* reason = "";

	XL_readStrCell(XL_bondId,C_bondId," ARM_ERR: bond id: object expected",C_result);
	XL_readNumCell(XL_settlement,C_settlement," ARM_ERR: settlement: date expected",C_result);
	XL_readNumCell(XL_yield,C_yield," ARM_ERR: yield: numeric expected",C_result);
	
	if((ARM_CheckDate(C_settlement, C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	long retCode = ICMLOCAL_YTOPRICE (LocalGetNumObjectId (C_bondId), C_settlement, C_yield, C_result);

	if(retCode == ARM_OK)
	{
		FreeCurCellErr ();
		XL_result.xltype = xltypeNum;
		XL_result.val.num = C_result.getDouble ();
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


__declspec(dllexport) LPXLOPER WINAPI ARM_Credit_PTOYIELD (LPXLOPER XL_bondId,
													  LPXLOPER XL_settlement,
													  LPXLOPER XL_price)
{
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	try {

	ARM_NOCALCIFWIZ();

	// C variable
	CCString C_bondId;
	double C_settlement;
	double C_price;
	
	// error
	static int error;
	static char* reason = "";

	XL_readStrCell(XL_bondId,C_bondId," ARM_ERR: bond id: object expected",C_result);
	XL_readNumCell(XL_settlement,C_settlement," ARM_ERR: settlement: date expected",C_result);
	XL_readNumCell(XL_price,C_price," ARM_ERR: price: numeric expected",C_result);
	
	if((ARM_CheckDate(C_settlement, C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	long retCode = ICMLOCAL_PTOYIELD (LocalGetNumObjectId (C_bondId), C_settlement, C_price, C_result);

	if(retCode == ARM_OK)
	{
		FreeCurCellErr ();
		XL_result.xltype = xltypeNum;
		XL_result.val.num = C_result.getDouble ();
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



/*---- End Of File ----*/

// EOF %M% 
