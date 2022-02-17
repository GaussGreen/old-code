#include <ARM\libarm_local\firstToBeIncluded.h>

#include <libCCxll\CCxll.h>

#include <ARM\libarm_local\ARM_local_swap.h>
#include <ARM\libarm_local\ARM_local_glob.h>
#include <ARM\libarm_local\ARM_local_class.h>

#include "ARM_local_interface.h"
#include "ARM_local_interglob.h"
#include "ARM_xl_trycatch_local.h"
#include <util\fromto.h>
#include <ARM\local_xlarm\ARM_xl_wrapper_local.h>

#include <ARM\local_xlarm\ARM_xl_gp_fctorhelper.h>

#include "util\tech_macro.h"

__declspec(dllexport) LPXLOPER WINAPI Local_LiborAssetSwapMargin (LPXLOPER XL_model,
																  LPXLOPER XL_startDate,
																  LPXLOPER XL_endDate,
																  LPXLOPER XL_fixedRate,
																  LPXLOPER XL_fixDayCount,
																  LPXLOPER XL_fixFrequency,
																  LPXLOPER XL_fixDecompFrequency,
																  LPXLOPER XL_fixPayTiming,
																  LPXLOPER XL_fixIntRule,
																  LPXLOPER XL_liborType,
																  LPXLOPER XL_floatResetFreq,
																  LPXLOPER XL_floatPayFreq,
																  LPXLOPER XL_price,
																  LPXLOPER XL_discountCcy,
																  LPXLOPER XL_redemptionPrice,
																  LPXLOPER XL_supplFee)
{
	ADD_LOG("Local_LiborAssetSwapMargin ");
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
		ARM_NOCALCIFWIZ();

		// C variable
		CCString C_model;
		
		double C_startDate;
		double C_endDate;

		double C_fixedRate;
		
		CCString C_fixDayCount;
		long fixDayCountId;

		CCString C_fixFrequency;
		long fixFrequencyId;

		CCString C_fixDecompFrequency;
		long fixDecompFrequencyId;

		CCString C_fixPayTiming;
		long fixPayTimingId;

		CCString C_fixIntRule;
		long fixIntRuleId;

		CCString C_liborType;
		long liborTypeId;

		CCString C_floatResetFreq;
		long floatResetFreqId;

		CCString C_floatPayFreq;
		long floatPayFreqId;

		double C_price;

		CCString C_discountCcy;
			
		double C_redemptionPrice;
		double C_redemptionPrice_default = -1;

		double C_supplFee;
		double C_supplFee_default = 0;
		
		double C_spread = 0;
		double C_assetGap = 3;

		CCString C_receiveOrPay ("R");
		long receiveOrPayId;
			
		// error
		static int error;
		static char* reason = "";
		
		XL_readStrCell(XL_model,C_model," ARM_ERR: model id: object expected",C_result);
		XL_readNumCell(XL_startDate,C_startDate," ARM_ERR: start date: date expected",C_result);
		XL_readNumCell(XL_endDate,C_endDate," ARM_ERR: end date: date expected",C_result);
		XL_readNumCell(XL_fixedRate,C_fixedRate," ARM_ERR: fixed rate: numeric expected",C_result);
		XL_readStrCell(XL_fixDayCount,C_fixDayCount," ARM_ERR: fix day count: string expected",C_result);
		XL_readStrCell(XL_fixFrequency,C_fixFrequency," ARM_ERR: fix frequency: string expected",C_result);
		XL_readStrCell(XL_fixDecompFrequency,C_fixDecompFrequency," ARM_ERR: fix decomp frequency: string expected",C_result);
		XL_readStrCell(XL_fixPayTiming,C_fixPayTiming," ARM_ERR: fix pay timing: string expected",C_result);
		XL_readStrCell(XL_fixIntRule,C_fixIntRule," ARM_ERR: interpolation rule: string expected",C_result);
		XL_readStrCell(XL_liborType,C_liborType," ARM_ERR: libor type: string expected",C_result);
		XL_readStrCell(XL_floatResetFreq,C_floatResetFreq," ARM_ERR: float reset frequency: string expected",C_result);
		XL_readStrCell(XL_floatPayFreq,C_floatPayFreq," ARM_ERR: float pay frequency: string expected",C_result);
		XL_readNumCell(XL_price,C_price," ARM_ERR: price: numeric expected",C_result);
		XL_readStrCellWD(XL_discountCcy,C_discountCcy,"DEFAULT"," ARM_ERR: discount currency: string expected",C_result);
		XL_readNumCellWD(XL_redemptionPrice,C_redemptionPrice,C_redemptionPrice_default," ARM_ERR: redemption price: numeric expected",C_result);
		XL_readNumCellWD(XL_supplFee,C_supplFee,C_supplFee_default," ARM_ERR: supplement fee: numeric expected",C_result);
		
		if(C_discountCcy == "DEFAULT")
		{
			ARM_result currencyres;
			ARMLOCAL_ARM_GetDefaultCurrency (currencyres);
			if(currencyres.getRetCode () != ARM_OK)
			{
				ARM_ARG_ERR();
				return (LPXLOPER)&XL_result;
			}
			else
			{
				C_discountCcy = currencyres.getString ();
			}
		}

		fixDayCountId = ARM_ConvDayCount (C_fixDayCount);

		if((fixFrequencyId = ARM_ConvFrequency (C_fixFrequency, C_result)) == ARM_DEFAULT_ERR)
		{
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}

		fixDecompFrequencyId = ARM_ConvDecompFrequency (C_fixDecompFrequency);

		fixPayTimingId = ARM_ConvPayResetRule (C_fixPayTiming);

		fixIntRuleId = ARM_ConvIntRule (C_fixIntRule);

		if((liborTypeId = ARM_ConvIrIndName (C_liborType, C_result)) == ARM_DEFAULT_ERR)
		{
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}

		if((floatResetFreqId = ARM_ConvFrequency (C_floatResetFreq, C_result)) == ARM_DEFAULT_ERR)
		{
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}
		
		if((floatPayFreqId = ARM_ConvFrequency (C_floatPayFreq, C_result)) == ARM_DEFAULT_ERR)
		{
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}

		if((receiveOrPayId = ARM_ConvRecOrPay (C_receiveOrPay, C_result)) == ARM_DEFAULT_ERR)
		{
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}

		long retCode = ARM_KO;

		char* id = "";

		retCode = ARMLOCAL_LiborAssetSwapMargin (LocalGetNumObjectId (C_model),C_startDate,
											C_endDate, C_fixedRate, receiveOrPayId,
											fixDayCountId, fixFrequencyId,
											fixDecompFrequencyId, fixPayTimingId,
											fixIntRuleId, liborTypeId,
											C_spread, floatResetFreqId,
											floatPayFreqId,(long)C_assetGap,0,
											C_price, C_discountCcy,C_redemptionPrice,
											C_supplFee,0,id,C_result);

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
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_LiborAssetSwapMargin" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}



__declspec(dllexport) LPXLOPER WINAPI Local_LIBORSWAP (LPXLOPER XL_startDate,
													   LPXLOPER XL_endDate,
													   LPXLOPER XL_liborType,
													   LPXLOPER XL_receiveOrPay,
													   LPXLOPER XL_fixedRate,
													   LPXLOPER XL_spread,
													   LPXLOPER XL_ccy,
													   LPXLOPER XL_fixedDayCount,
													   LPXLOPER XL_floatingDayCount)
{
	ADD_LOG("Local_LIBORSWAP ");
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
		
		CCString C_liborType;
		long liborTypeId;

		CCString C_receiveOrPay;
		long receiveOrPayId;

		double C_fixedRate;

		CCString C_spread_str;
		double C_spread;
		double C_spread_default = 0.0;
		long spreadType;
			
		CCString C_ccy;
		bool ccyIsObject = false;
		
		CCString C_fixedDayCount;
		CCString C_floatingDayCount;
		long fixedDayCountId;
		long floatingDayCountId;

		// error
		static int error;
		static char* reason = "";

		XL_readNumCell(XL_startDate, C_startDate, " ARM_ERR: start date: date expected",C_result);
		XL_readNumCell(XL_endDate, C_endDate, " ARM_ERR: end date: date expected",C_result);
		XL_readStrCell(XL_liborType, C_liborType, " ARM_ERR: libor type: string expected",C_result);
		XL_readStrCell(XL_receiveOrPay, C_receiveOrPay, " ARM_ERR: receive or pay: string expected",C_result);
		XL_readNumCell(XL_fixedRate, C_fixedRate, " ARM_ERR: fixed rate: numeric expected",C_result);
		XL_readStrOrNumCellWD(XL_spread, C_spread_str, C_spread, C_spread_default, spreadType,
			   " ARM_ERR: spread: numeric or object ID string expected",C_result);
	//	XL_readNumCellWD(XL_spread, C_spread, C_spread_default," ARM_ERR: spread: numeric expected",C_result);
		XL_readStrCellWD(XL_ccy,C_ccy, "DEFAULT", " ARM_ERR: currency id: string expected",C_result);
		XL_readStrCellWD(XL_fixedDayCount,C_fixedDayCount,"-1"," ARM_ERR: fixedDayCount: string expected",C_result);
		XL_readStrCellWD(XL_floatingDayCount,C_floatingDayCount,"-1"," ARM_ERR: floatingDayCount: string expected",C_result);
		

		if (( C_ccy.GetLen() > 3 )
			&&
			(!( C_ccy == "DEFAULT" ))
		   )
		   ccyIsObject = true;

		if ((liborTypeId = ARM_ConvIrIndName(C_liborType, C_result)) == ARM_DEFAULT_ERR )
		{
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}

		if ((receiveOrPayId = ARM_ConvRecOrPay(C_receiveOrPay, C_result)) == ARM_DEFAULT_ERR)
		{
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}

		if ( spreadType == XL_TYPE_STRING )
		{
		   C_spread = (double) LocalGetNumObjectId(C_spread_str);

		   spreadType = 1L;
		}
		else
		{
		   spreadType = 0L;
		}
		
		fixedDayCountId = ARM_ConvDayCount(C_fixedDayCount);
		floatingDayCountId = ARM_ConvDayCount(C_floatingDayCount);

		long retCode;
		long objId;
		CCString prevClass;
		
		CCString curClass = LOCAL_SWAP_CLASS;
		CCString stringId = GetLastCurCellEnvValue ();
		
		if (!stringId)
		{
			retCode = ARMLOCAL_LIBORSWAP(C_startDate, C_endDate, liborTypeId, 
										 receiveOrPayId,C_fixedRate,spreadType,
										 C_spread, ccyIsObject, C_ccy, (int) fixedDayCountId, (int) floatingDayCountId, C_result);

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
				retCode = ARMLOCAL_LIBORSWAP(C_startDate, C_endDate, liborTypeId, 
											 receiveOrPayId, C_fixedRate,spreadType,
											 C_spread, ccyIsObject, C_ccy, (int) fixedDayCountId, (int) floatingDayCountId, C_result, objId);

				if(retCode == ARM_OK)
				{			
					LocalSetCurCellEnvValue (curClass, objId); 

					stringId = LocalMakeObjectId (objId, curClass);
				}
			}
			else
			{
				FreeCurCellContent();

				retCode = ARMLOCAL_LIBORSWAP(C_startDate, C_endDate, liborTypeId, 
											 receiveOrPayId, C_fixedRate,spreadType,
											 C_spread, ccyIsObject, C_ccy, (int) fixedDayCountId, (int) floatingDayCountId, C_result);
			
				if ( retCode == ARM_OK )
				{
					objId = C_result.getLong ();
				
					LocalSetCurCellEnvValue (curClass, objId); 

					stringId = LocalMakeObjectId (objId, curClass);
				}
			}
		}

		if ( retCode == ARM_OK )
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_LIBORSWAP" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}



__declspec(dllexport) LPXLOPER WINAPI Local_PXL_LIBORSWAP (LPXLOPER XL_startDate,
														   LPXLOPER XL_endDate,
														   LPXLOPER XL_liborType,
														   LPXLOPER XL_receiveOrPay,
														   LPXLOPER XL_fixedRate,
														   LPXLOPER XL_spread,
														   LPXLOPER XL_ccy,
														   LPXLOPER XL_fixedDayCount,
														   LPXLOPER XL_floatingDayCount)
{
	ADD_LOG("Local_PXL_LIBORSWAP ");
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
		
		CCString C_liborType;
		long liborTypeId;

		CCString C_receiveOrPay;
		long receiveOrPayId;

		double C_fixedRate;

		CCString C_spread_str;
		double C_spread;
		double C_spread_default = 0.0;
		long spreadType;
			
		CCString C_ccy;
		bool ccyIsObject = false;
		
		CCString C_fixedDayCount;
		CCString C_floatingDayCount;
		long fixedDayCountId;
		long floatingDayCountId;

		// error
		static int error;
		static char* reason = "";

		XL_readNumCell(XL_startDate,C_startDate," ARM_ERR: start date: date expected",C_result);
		XL_readNumCell(XL_endDate,C_endDate," ARM_ERR: end date: date expected",C_result);
		XL_readStrCell(XL_liborType,C_liborType," ARM_ERR: libor type: string expected",C_result);
		XL_readStrCell(XL_receiveOrPay,C_receiveOrPay," ARM_ERR: receive or pay: string expected",C_result);
		XL_readNumCell(XL_fixedRate,C_fixedRate," ARM_ERR: fixed rate: numeric expected",C_result);
		XL_readStrOrNumCellWD(XL_spread, C_spread_str, C_spread, C_spread_default, spreadType,
			   " ARM_ERR: spread: numeric or object ID string expected",C_result);
	//	XL_readNumCellWD(XL_spread,C_spread,C_spread_default," ARM_ERR: spread: numeric expected",C_result);
		XL_readStrCellWD(XL_ccy,C_ccy,"DEFAULT"," ARM_ERR: currency id: string expected",C_result);
		XL_readStrCellWD(XL_fixedDayCount,C_fixedDayCount,"-1"," ARM_ERR: fixedDayCount: string expected",C_result);
		XL_readStrCellWD(XL_floatingDayCount,C_floatingDayCount,"-1"," ARM_ERR: floatingDayCount: string expected",C_result);
		
		if (( C_ccy.GetLen() > 3 )
			&&
			(!( C_ccy == "DEFAULT" ))
		   )
		   ccyIsObject = true;

		if((liborTypeId = ARM_ConvIrIndName (C_liborType, C_result)) == ARM_DEFAULT_ERR)
		{
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}

		if((receiveOrPayId = ARM_ConvRecOrPay (C_receiveOrPay, C_result)) == ARM_DEFAULT_ERR)
		{
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}

		if ( spreadType == XL_TYPE_STRING )
		{
		   C_spread = (double) LocalGetNumObjectId(C_spread_str);

		   spreadType = 1L;
		}
		else
		{
		   spreadType = 0L;
		}
		
		fixedDayCountId = ARM_ConvDayCount(C_fixedDayCount);
		floatingDayCountId = ARM_ConvDayCount(C_floatingDayCount);

		long retCode;
		long objId;
		
		CCString curClass = LOCAL_SWAP_CLASS;
		CCString stringId;
		
		retCode = ARMLOCAL_LIBORSWAP (C_startDate, C_endDate, (long)liborTypeId, 
									 (long)receiveOrPayId, C_fixedRate,spreadType,
									 C_spread, ccyIsObject, C_ccy, (int) fixedDayCountId, (int) floatingDayCountId, C_result);

		if ( retCode == ARM_OK )
		{
			objId = C_result.getLong ();

			stringId = LocalMakeObjectId (objId, curClass);
		}

		if ( retCode == ARM_OK )
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_PXL_LIBORSWAP" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


__declspec(dllexport) LPXLOPER WINAPI Local_LIBORLEG (LPXLOPER XL_startDate,
													  LPXLOPER XL_endDate,
													  LPXLOPER XL_liborType,
													  LPXLOPER XL_receiveOrPay,
													  LPXLOPER XL_spread,
													  LPXLOPER XL_resetFreq,
													  LPXLOPER XL_payFreq,
													  LPXLOPER XL_resetTiming,
													  LPXLOPER XL_payTiming,
													  LPXLOPER XL_ccy,
													  LPXLOPER XL_intRule,
													  LPXLOPER XL_resetGap,
													  LPXLOPER XL_resetCal,
													  LPXLOPER XL_payCal,
													  LPXLOPER XL_decompPricingFlag,
													  LPXLOPER XL_nxChange,
													  LPXLOPER XL_stubRule,
													  LPXLOPER XL_refDate,
													  LPXLOPER XL_adjStartDate,
													  LPXLOPER XL_couponDayCount)
{
	ADD_LOG("Local_LIBORLEG ");
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

		CCString C_liborType;
		long liborTypeId;
		
		CCString C_receiveOrPay;
		long receiveOrPayId;

		CCString C_spread_str;
		double C_spread_double;
		double C_default_spread = 0;
		long   spreadType;

		CCString C_resetFreq;
		long resetFreqId;

		CCString C_payFreq;
		long payFreqId;

		CCString C_resetTiming;
		long resetTimingId;
		
		CCString C_payTiming;
		long payTimingId;
		
		CCString C_ccy;
		bool ccyIsObject = false;

		CCString C_intRule;
		long intRuleId;

		double C_resetGap;
		double C_resetGap_default = 10000.;

		CCString C_resetCal;
		CCString C_payCal;
		
		double C_decompPricingFlag;
		double C_decompPricingFlag_default(1.0);

		CCString C_nxChange;
		long     nxChange;

		CCString C_stubRule;
		long stubRuleId;

		double C_refDate;
		double C_refDate_default = -1.0;

		CCString C_adjStartDate;
		long     adjStartDateId;
		
		CCString C_couponDayCount;
		long dayCountId;

		// error
		static int error;
		static char* reason = "";

		XL_readNumCell(XL_startDate,C_startDate," ARM_ERR: start date: date expected",C_result);
		XL_readNumCell(XL_endDate,C_endDate," ARM_ERR: end date: date expected",C_result);
		XL_readStrCell(XL_liborType,C_liborType," ARM_ERR: libor type: string expected",C_result);
		XL_readStrCell(XL_receiveOrPay,C_receiveOrPay," ARM_ERR: receive or pay: string expected",C_result);
		XL_readStrOrNumCellWD(XL_spread, C_spread_str, C_spread_double, C_default_spread, spreadType,
			   " ARM_ERR: spread: numeric or object ID string expected",C_result);
		XL_readStrCellWD(XL_resetFreq,C_resetFreq,"-1"," ARM_ERR: reset frequency: string expected",C_result);
		XL_readStrCellWD(XL_payFreq,C_payFreq,"-1"," ARM_ERR: pay frequency: string expected",C_result);
		XL_readStrCellWD(XL_resetTiming,C_resetTiming,"ADV"," ARM_ERR: reset timing: string expected",C_result);
		XL_readStrCellWD(XL_payTiming,C_payTiming,"ARR"," ARM_ERR: pay timing: string expected",C_result);
		XL_readStrCellWD(XL_ccy,C_ccy,"DEFAULT"," ARM_ERR: currency id: string expected",C_result);
		XL_readStrCellWD(XL_intRule,C_intRule,"ADJ"," ARM_ERR: intRule: string expected",C_result);
		XL_readNumCellWD(XL_resetGap,C_resetGap,C_resetGap_default," ARM_ERR: reset Gap: numeric expected",C_result);
		XL_readStrCellWD(XL_resetCal,C_resetCal,"NULL"," ARM_ERR: reset Calendar: string expected",C_result);
		XL_readStrCellWD(XL_payCal,C_payCal,"NULL"," ARM_ERR: pay Calendar: string expected",C_result);
		XL_readNumCellWD(XL_decompPricingFlag,C_decompPricingFlag,C_decompPricingFlag_default," ARM_ERR: decompPricingFlag: numeric expected",C_result);
		XL_readStrCellWD(XL_nxChange,C_nxChange,"NXNONE"," ARM_ERR: Notional xchge: string expected",C_result);
		XL_readStrCellWD(XL_stubRule,C_stubRule,"SS"," ARM_ERR: stub Rule: string expected",C_result);
		XL_readNumCellWD(XL_refDate,C_refDate,C_refDate_default," ARM_ERR: reference date: date expected",C_result);
		XL_readStrCellWD(XL_adjStartDate,C_adjStartDate,"YES"," ARM_ERR: adjStartDate: string expected",C_result);
		XL_readStrCellWD(XL_couponDayCount,C_couponDayCount,"-1"," ARM_ERR: couponDayCount: string expected",C_result);
		

		if ((liborTypeId = ARM_ConvIrIndName (C_liborType, C_result)) == ARM_DEFAULT_ERR)
		{
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}

		if ((receiveOrPayId = ARM_ConvRecOrPay (C_receiveOrPay, C_result)) == ARM_DEFAULT_ERR)
		{
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}

		if ((payFreqId = ARM_ConvFrequency (C_payFreq, C_result)) == ARM_DEFAULT_ERR)
		{
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}

		if ((resetFreqId = ARM_ConvFrequency (C_resetFreq, C_result)) == ARM_DEFAULT_ERR)
		{
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}

		payTimingId = ARM_ConvPayResetRule (C_payTiming);

		resetTimingId = ARM_ConvPayResetRule (C_resetTiming);

		if (( C_ccy.GetLen() > 3 )
			&& 
			( !(C_ccy == "DEFAULT"))
		   )
		   ccyIsObject = true;

		intRuleId = ARM_ConvIntRule (C_intRule);

		nxChange = ARM_NotionalExchange(C_nxChange);

		stubRuleId = ARM_ConvStubRule (C_stubRule);

		if ( spreadType == XL_TYPE_STRING )
		{
		   C_spread_double = (double) LocalGetNumObjectId(C_spread_str);

		   spreadType = 1L;
		}
		else
		{
		   spreadType = 0L;
		}

		if ((adjStartDateId = ARM_ConvYesOrNo(C_adjStartDate, C_result)) == ARM_DEFAULT_ERR)
		{
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}
		
		dayCountId = ARM_ConvDayCount(C_couponDayCount);

		long retCode;
		long objId;
		CCString prevClass;
		
		CCString curClass = LOCAL_SWAPLEG_CLASS;
		CCString stringId = GetLastCurCellEnvValue ();
		
		if (!stringId)
		{
			retCode = ARMLOCAL_LIBORLEG (C_startDate,
										 C_endDate,
										 liborTypeId,
										 receiveOrPayId,
										 spreadType,
										 C_spread_double,
										 resetFreqId,
										 payFreqId,
										 resetTimingId,
										 payTimingId,
										 ccyIsObject,
										 C_ccy,
										 intRuleId,
										 (long) C_resetGap,
										 C_resetCal,
										 C_payCal,
										 (long) C_decompPricingFlag,
										 nxChange,
										 stubRuleId,
										 C_refDate,
										 adjStartDateId,
										 (int) dayCountId,
										 C_result);

			if (retCode == ARM_OK)
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
				
			if (curClass == prevClass)
			{
				retCode = ARMLOCAL_LIBORLEG (C_startDate,
											 C_endDate,
											 liborTypeId,
											 receiveOrPayId,
											 spreadType,
											 C_spread_double,
											 resetFreqId,
											 payFreqId,
											 resetTimingId,
											 payTimingId,
											 ccyIsObject,
											 C_ccy,
											 intRuleId,
											 (long) C_resetGap,
											 C_resetCal,
											 C_payCal,
											 (long) C_decompPricingFlag,
											 nxChange,
											 stubRuleId,
											 C_refDate,
											 adjStartDateId,
											 (int) dayCountId,
											 C_result,
											 objId);

				if (retCode == ARM_OK)
				{
					LocalSetCurCellEnvValue (curClass, objId); 

					stringId = LocalMakeObjectId (objId, curClass);
				}
			}
			else
			{
				FreeCurCellContent ();

				retCode = ARMLOCAL_LIBORLEG (C_startDate,
											 C_endDate,
											 liborTypeId,
											 receiveOrPayId,
											 spreadType,
											 C_spread_double,
											 resetFreqId,
											 payFreqId,
											 resetTimingId,
											 payTimingId,
											 ccyIsObject,
											 C_ccy,
											 intRuleId,
											 (long) C_resetGap,
											 C_resetCal,
											 C_payCal,
											 (long) C_decompPricingFlag,
											 nxChange,
											 stubRuleId,
											 C_refDate,
											 adjStartDateId,
											 (int) dayCountId,
											 C_result);

				if (retCode == ARM_OK)
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_LIBORLEG" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}




__declspec(dllexport) LPXLOPER WINAPI Local_PXL_LIBORLEG (LPXLOPER XL_startDate,
														  LPXLOPER XL_endDate,
														  LPXLOPER XL_liborType,
														  LPXLOPER XL_receiveOrPay,
														  LPXLOPER XL_spread,
														  LPXLOPER XL_resetFreq,
														  LPXLOPER XL_payFreq,
														  LPXLOPER XL_resetTiming,
														  LPXLOPER XL_payTiming,
														  LPXLOPER XL_ccy,
														  LPXLOPER XL_intRule,
														  LPXLOPER XL_resetGap,
														  LPXLOPER XL_resetCal,
														  LPXLOPER XL_payCal,
														  LPXLOPER XL_decompPricingFlag,
														  LPXLOPER XL_nxChange,
														  LPXLOPER XL_stubRule,
														  LPXLOPER XL_refDate,
														  LPXLOPER XL_adjStartDate,
														  LPXLOPER XL_couponDayCount)
{
	ADD_LOG("Local_PXL_LIBORLEG ");
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

		CCString C_liborType;
		long liborTypeId;
		
		CCString C_receiveOrPay;
		long receiveOrPayId;

		CCString C_spread_str;
		double C_spread_double;
		double C_default_spread = 0;
		long   spreadType;
		
		CCString C_resetFreq;
		long resetFreqId;

		CCString C_payFreq;
		long payFreqId;

		CCString C_resetTiming;
		long resetTimingId;
		
		CCString C_payTiming;
		long payTimingId;
		
		CCString C_ccy;
		bool ccyIsObject = false;
				
		CCString C_intRule;
		long intRuleId;

		double C_resetGap;
		double C_resetGap_default = 10000.;

		CCString C_resetCal;
		CCString C_payCal;
		
		double C_decompPricingFlag;
		double C_decompPricingFlag_default(1.0);

		CCString C_nxChange;
		long     nxChange;

		CCString C_stubRule;
		long stubRuleId;

		double C_refDate;
		double C_refDate_default = -1.0;

		CCString C_adjStartDate;
		long     adjStartDateId;

		CCString C_couponDayCount;
		long dayCountId;	

		// error
		static int error;
		static char* reason = "";

		XL_readNumCell(XL_startDate,C_startDate," ARM_ERR: start date: date expected",C_result);
		XL_readNumCell(XL_endDate,C_endDate," ARM_ERR: end date: date expected",C_result);
		XL_readStrCell(XL_liborType,C_liborType," ARM_ERR: libor type: string expected",C_result);
		XL_readStrCell(XL_receiveOrPay,C_receiveOrPay," ARM_ERR: receive or pay: string expected",C_result);
		XL_readStrOrNumCellWD(XL_spread, C_spread_str, C_spread_double, C_default_spread, spreadType,
			   " ARM_ERR: spread: numeric or object ID string expected",C_result);
		XL_readStrCellWD(XL_resetFreq,C_resetFreq,"-1"," ARM_ERR: reset frequency: string expected",C_result);
		XL_readStrCellWD(XL_payFreq,C_payFreq,"-1"," ARM_ERR: pay frequency: string expected",C_result);
		XL_readStrCellWD(XL_resetTiming,C_resetTiming,"ADV"," ARM_ERR: reset timing: string expected",C_result);
		XL_readStrCellWD(XL_payTiming,C_payTiming,"ARR"," ARM_ERR: pay timing: string expected",C_result);
		XL_readStrCellWD(XL_ccy,C_ccy,"DEFAULT"," ARM_ERR: currency id: string expected",C_result);
		XL_readStrCellWD(XL_intRule,C_intRule,"ADJ"," ARM_ERR: intRule: string expected",C_result);
		XL_readNumCellWD(XL_resetGap,C_resetGap,C_resetGap_default," ARM_ERR: reset Gap: numeric expected",C_result);
		XL_readStrCellWD(XL_resetCal,C_resetCal,"NULL"," ARM_ERR: reset Calendar: string expected",C_result);
		XL_readStrCellWD(XL_payCal,C_payCal,"NULL"," ARM_ERR: pay Calendar: string expected",C_result);
		XL_readNumCellWD(XL_decompPricingFlag,C_decompPricingFlag,C_decompPricingFlag_default," ARM_ERR: decompPricingFlag: numeric expected",C_result);
		XL_readStrCellWD(XL_nxChange,C_nxChange,"NXNONE"," ARM_ERR: Notional xchge: string expected",C_result);
		XL_readStrCellWD(XL_stubRule,C_stubRule,"SS"," ARM_ERR: stub Rule: string expected",C_result);
		XL_readNumCellWD(XL_refDate,C_refDate,C_refDate_default," ARM_ERR: reference date: date expected",C_result);
		XL_readStrCellWD(XL_adjStartDate,C_adjStartDate,"YES"," ARM_ERR: adjStartDate: string expected",C_result);
		XL_readStrCellWD(XL_couponDayCount,C_couponDayCount,"-1"," ARM_ERR: couponDayCount: string expected",C_result);
		

		if((liborTypeId = ARM_ConvIrIndName (C_liborType, C_result)) == ARM_DEFAULT_ERR)
		{
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}

		if((receiveOrPayId = ARM_ConvRecOrPay (C_receiveOrPay, C_result)) == ARM_DEFAULT_ERR)
		{
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}

		if((payFreqId = ARM_ConvFrequency (C_payFreq, C_result)) == ARM_DEFAULT_ERR)
		{
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}

		if((resetFreqId = ARM_ConvFrequency (C_resetFreq, C_result)) == ARM_DEFAULT_ERR)
		{
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}

		payTimingId = ARM_ConvPayResetRule (C_payTiming);

		resetTimingId = ARM_ConvPayResetRule (C_resetTiming);

		if (( C_ccy.GetLen() > 3 )
			&& 
			( !( C_ccy == "DEFAULT" ))
		   )
		   ccyIsObject = true;

		intRuleId = ARM_ConvIntRule (C_intRule);

		nxChange = ARM_NotionalExchange(C_nxChange);

		stubRuleId = ARM_ConvStubRule (C_stubRule);

		if ( spreadType == XL_TYPE_STRING )
		{
		   C_spread_double = (double) LocalGetNumObjectId(C_spread_str);

		   spreadType = 1L;
		}
		else
		{
		   spreadType = 0L;
		}

		if ((adjStartDateId = ARM_ConvYesOrNo(C_adjStartDate, C_result)) == ARM_DEFAULT_ERR)
		{
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}
		
		dayCountId = ARM_ConvDayCount(C_couponDayCount);

		long retCode;
		long objId;
		
		CCString curClass = LOCAL_SWAPLEG_CLASS;
		CCString stringId;
		
		retCode = ARMLOCAL_LIBORLEG (C_startDate,
									 C_endDate,
									 liborTypeId,
									 receiveOrPayId,
									 spreadType,
									 C_spread_double,
									 resetFreqId,
									 payFreqId,
									 resetTimingId,
									 payTimingId,
									 ccyIsObject,
									 C_ccy,
									 intRuleId,
									 (long) C_resetGap,
									 C_resetCal,
									 C_payCal,
									 (long) C_decompPricingFlag,
									 nxChange,
									 stubRuleId,
									 C_refDate,
									 adjStartDateId,
									 (int) dayCountId,
									 C_result);

		if ( retCode == ARM_OK )
		{
			objId = C_result.getLong ();

			stringId = LocalMakeObjectId (objId, curClass);
		}
		
		if ( retCode == ARM_OK )
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_PXL_LIBORLEG" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
	
}



__declspec(dllexport) LPXLOPER WINAPI Local_SWAPLEG(LPXLOPER XL_irIndex,
													LPXLOPER XL_startDate,
													LPXLOPER XL_endDate,
													LPXLOPER XL_receiveOrPay,
													LPXLOPER XL_spread,
													LPXLOPER XL_ccy,
													LPXLOPER XL_dayCount,
													LPXLOPER XL_resetGap,
													LPXLOPER XL_resetCal,
													LPXLOPER XL_payCal,
													LPXLOPER XL_decompPricingFlag,
													LPXLOPER XL_nxChange,
													LPXLOPER XL_stubRule,
													LPXLOPER XL_refDate,
													LPXLOPER XL_adjStartDate,
													LPXLOPER XL_rollDay)
{
	ADD_LOG("Local_SWAPLEG");
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
		ARM_NOCALCIFWIZ();

		// C variable
		CCString C_irIndex;

		double C_startDate;
		double C_endDate;
		
		CCString C_receiveOrPay;
		long receiveOrPayId;

		double C_spread_double;
		CCString C_spread_str;
		long   spreadType;
		
		CCString C_ccy;
		bool ccyIsObject = false;

		CCString C_dayCount;
		long dayCountId;

		double C_resetGap;
		double C_resetGap_default = 10000.;

		CCString C_resetCal;
		CCString C_payCal;
		
		double C_decompPricingFlag;
		double C_decompPricingFlag_default(1.0);

		CCString C_nxChange;
		long     nxChange;

		CCString C_stubRule;
		long stubRuleId;

		double C_refDate;
		double C_refDate_default = -1.0;

		CCString C_adjStartDate;
		long     adjStartDateId;

		double	 C_rollDay;
		double   C_rollDay_default = -1.0;
		long     rollDay;

		// error
		static int error;
		static char* reason = "";

		XL_readStrCell(XL_irIndex,C_irIndex," ARM_ERR: interest rate index id: object expected",C_result);
		XL_readNumCell(XL_startDate,C_startDate," ARM_ERR: start date: date expected",C_result);
		XL_readNumCell(XL_endDate,C_endDate," ARM_ERR: end date: date expected",C_result);
		XL_readStrCell(XL_receiveOrPay,C_receiveOrPay," ARM_ERR: receive or pay: string expected",C_result);
		XL_readStrOrNumCell(XL_spread, C_spread_str, C_spread_double, spreadType,
			   " ARM_ERR: spread: numeric or object ID string expected",C_result);
		XL_readStrCellWD(XL_ccy,C_ccy,"DEFAULT"," ARM_ERR: currency id: string expected",C_result);
		XL_readStrCellWD(XL_dayCount,C_dayCount,"-1"," ARM_ERR: day count: string expected",C_result);
		XL_readNumCellWD(XL_resetGap,C_resetGap,C_resetGap_default," ARM_ERR: reset Gap: numeric expected",C_result);
		XL_readStrCellWD(XL_resetCal,C_resetCal,"NULL"," ARM_ERR: reset Calendar: string expected",C_result);
		XL_readStrCellWD(XL_payCal,C_payCal,"NULL"," ARM_ERR: pay Calendar: string expected",C_result);
		XL_readNumCellWD(XL_decompPricingFlag,C_decompPricingFlag,C_decompPricingFlag_default," ARM_ERR: decompPricingFlag: numeric expected",C_result);
		XL_readStrCellWD(XL_nxChange,C_nxChange,"NXNONE"," ARM_ERR: Notional xchge: string expected",C_result);
		XL_readStrCellWD(XL_stubRule,C_stubRule,"SS"," ARM_ERR: stub Rule: string expected",C_result);
		XL_readNumCellWD(XL_refDate,C_refDate,C_refDate_default," ARM_ERR: reference date: date expected",C_result);
		XL_readStrCellWD(XL_adjStartDate,C_adjStartDate,"YES"," ARM_ERR: adjStartDate: string expected",C_result);
		XL_readNumCellWD(XL_rollDay,C_rollDay, C_rollDay_default," ARM_ERR: rollDay: value between 1 and 31 expected",C_result);

		if (( C_ccy.GetLen() > 3 )
			&& 
			( !(C_ccy == "DEFAULT"))
		   )
		   ccyIsObject = true;

		if ((receiveOrPayId = ARM_ConvRecOrPay (C_receiveOrPay, C_result)) == ARM_DEFAULT_ERR)
		{
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}

		if ( spreadType == XL_TYPE_STRING )
		{
		   C_spread_double = (double) LocalGetNumObjectId(C_spread_str);

		   spreadType = 1L;
		}
		else
		{
		   spreadType = 0L;
		}

		dayCountId = ARM_ConvDayCount (C_dayCount);

		nxChange = ARM_NotionalExchange(C_nxChange);

		stubRuleId = ARM_ConvStubRule (C_stubRule);

		if ((adjStartDateId = ARM_ConvYesOrNo(C_adjStartDate, C_result)) == ARM_DEFAULT_ERR)
		{
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}

		rollDay	=	(long) C_rollDay;

		long retCode;
		long objId;
		CCString prevClass;
		
		CCString curClass = LOCAL_SWAPLEG_CLASS;
		CCString stringId = GetLastCurCellEnvValue ();
		
		if (!stringId)
		{
			retCode = ARMLOCAL_SWAPLEG (LocalGetNumObjectId (C_irIndex),
										C_startDate,
										C_endDate,
										(long)receiveOrPayId,
										spreadType,
										C_spread_double,
										ccyIsObject,
										C_ccy,
										(long)dayCountId,
										(long) C_resetGap,
										C_resetCal,
										C_payCal,
										(long) C_decompPricingFlag,
										nxChange,
										stubRuleId,
										C_refDate,
										adjStartDateId,
										rollDay,
										C_result);

			if ( retCode == ARM_OK )
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
				retCode = ARMLOCAL_SWAPLEG (LocalGetNumObjectId (C_irIndex),
											C_startDate,
											C_endDate,
											(long)receiveOrPayId,
											spreadType,
											C_spread_double,
											ccyIsObject,
											C_ccy,
											(long)dayCountId,
											(long) C_resetGap,
											C_resetCal,
											C_payCal,
											(long) C_decompPricingFlag,
											nxChange,
											stubRuleId,
											C_refDate,
											adjStartDateId,
											rollDay,
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

				retCode = ARMLOCAL_SWAPLEG (LocalGetNumObjectId (C_irIndex),
											C_startDate,
											C_endDate,
											(long)receiveOrPayId,
											spreadType,
											C_spread_double,
											ccyIsObject,
											C_ccy,
											(long)dayCountId,
											(long) C_resetGap,
											C_resetCal,
											C_payCal,
											(long) C_decompPricingFlag,
											nxChange,
											stubRuleId,
											C_refDate,
											adjStartDateId,
											rollDay,
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_SWAPLEG" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


__declspec(dllexport) LPXLOPER WINAPI Local_PXL_SWAPLEG (LPXLOPER XL_irIndex,
														 LPXLOPER XL_startDate,
														 LPXLOPER XL_endDate,
														 LPXLOPER XL_receiveOrPay,
														 LPXLOPER XL_spread,
														 LPXLOPER XL_ccy,
														 LPXLOPER XL_dayCount,
														 LPXLOPER XL_resetGap,
														 LPXLOPER XL_resetCal,
														 LPXLOPER XL_payCal,
														 LPXLOPER XL_decompPricingFlag,
														 LPXLOPER XL_nxChange,
														 LPXLOPER XL_stubRule,
														 LPXLOPER XL_refDate,
														 LPXLOPER XL_adjStartDate,
														 LPXLOPER XL_rollDay)
{
	ADD_LOG("Local_PXL_SWAPLEG ");
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
		ARM_NOCALCIFWIZ();

		// C variable
		CCString C_irIndex;

		double C_startDate;
		double C_endDate;
		
		CCString C_receiveOrPay;
		long receiveOrPayId;

		double C_spread_double;
		CCString C_spread_str;
		long   spreadType;
		
		CCString C_ccy;
		bool ccyIsObject = false;

		CCString C_dayCount;
		long dayCountId;

		double C_resetGap;
		double C_resetGap_default = 10000.;

		CCString C_resetCal;
		CCString C_payCal;
		
		double C_decompPricingFlag;
		double C_decompPricingFlag_default(1.0);

		CCString C_nxChange;
		long     nxChange;

		CCString C_stubRule;
		long stubRuleId;

		double C_refDate;
		double C_refDate_default = -1.0;

		CCString C_adjStartDate;
		long     adjStartDateId;

		double	 C_rollDay;
		double   C_rollDay_default = -1.0;
		long     rollDay;

		// error
		static int error;
		static char* reason = "";

		XL_readStrCell(XL_irIndex,C_irIndex," ARM_ERR: interest rate index id: object expected",C_result);
		XL_readNumCell(XL_startDate,C_startDate," ARM_ERR: start date: date expected",C_result);
		XL_readNumCell(XL_endDate,C_endDate," ARM_ERR: end date: date expected",C_result);
		XL_readStrCell(XL_receiveOrPay,C_receiveOrPay," ARM_ERR: receive or pay: string expected",C_result);
		XL_readStrOrNumCell(XL_spread, C_spread_str, C_spread_double, spreadType,
			   " ARM_ERR: spread: numeric or object ID string expected",C_result);
		XL_readStrCellWD(XL_ccy,C_ccy,"DEFAULT"," ARM_ERR: currency id: string expected",C_result);
		XL_readStrCellWD(XL_dayCount,C_dayCount,"-1"," ARM_ERR: day count: string expected",C_result);
		XL_readNumCellWD(XL_resetGap,C_resetGap,C_resetGap_default," ARM_ERR: reset Gap: numeric expected",C_result);
		XL_readStrCellWD(XL_resetCal,C_resetCal,"NULL"," ARM_ERR: reset Calendar: string expected",C_result);
		XL_readStrCellWD(XL_payCal,C_payCal,"NULL"," ARM_ERR: pay Calendar: string expected",C_result);
		XL_readNumCellWD(XL_decompPricingFlag,C_decompPricingFlag,C_decompPricingFlag_default," ARM_ERR: decompPricingFlag: numeric expected",C_result);
		XL_readStrCellWD(XL_nxChange,C_nxChange,"NXNONE"," ARM_ERR: Notional xchge: string expected",C_result);
		XL_readStrCellWD(XL_stubRule,C_stubRule,"SS"," ARM_ERR: stub Rule: string expected",C_result);
		XL_readNumCellWD(XL_refDate,C_refDate,C_refDate_default," ARM_ERR: reference date: date expected",C_result);
		XL_readStrCellWD(XL_adjStartDate,C_adjStartDate,"YES"," ARM_ERR: adjStartDate: string expected",C_result);
		XL_readNumCellWD(XL_rollDay,C_rollDay, C_rollDay_default," ARM_ERR: rollDay: value between 1 and 31 expected",C_result);

		if ((C_ccy.GetLen() > 3)
			&& !(C_ccy == "DEFAULT"))
			ccyIsObject = true;

		if((receiveOrPayId = ARM_ConvRecOrPay (C_receiveOrPay, C_result)) == ARM_DEFAULT_ERR)
		{
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}

		if ( spreadType == XL_TYPE_STRING )
		{
		   C_spread_double = (double) LocalGetNumObjectId(C_spread_str);

		   spreadType = 1L;
		}
		else
		{
		   spreadType = 0L;
		}

		dayCountId = ARM_ConvDayCount (C_dayCount);

		nxChange = ARM_NotionalExchange(C_nxChange);

		stubRuleId = ARM_ConvStubRule (C_stubRule);

		if ((adjStartDateId = ARM_ConvYesOrNo(C_adjStartDate, C_result)) == ARM_DEFAULT_ERR)
		{
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}

		rollDay	=	(long) C_rollDay;

		long retCode;
		long objId;
		
		CCString curClass = LOCAL_SWAPLEG_CLASS;
		CCString stringId;
		
		retCode = ARMLOCAL_SWAPLEG (LocalGetNumObjectId (C_irIndex),
									C_startDate,
									C_endDate,
									(long)receiveOrPayId,
									spreadType,
									C_spread_double,
									ccyIsObject,
									C_ccy,
									(long)dayCountId,
									(long) C_resetGap,
									C_resetCal,
									C_payCal,
									(long) C_decompPricingFlag,
									nxChange,
									stubRuleId,
									C_refDate,
									adjStartDateId,
									rollDay,
									C_result);

		if(retCode == ARM_OK)
		{
			objId = C_result.getLong ();

			stringId = LocalMakeObjectId (objId, curClass);
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

///	ARM_END();
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_PXL_SWAPLEG" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


// tatata 
__declspec(dllexport) LPXLOPER WINAPI Local_FIXEDLEG (LPXLOPER XL_startDate,
													  LPXLOPER XL_endDate,
													  LPXLOPER XL_receiveOrPay,
													  LPXLOPER XL_fixRate,
													  LPXLOPER XL_dayCount,
													  LPXLOPER XL_freq,
													  LPXLOPER XL_decompFreq,
													  LPXLOPER XL_payTiming,
													  LPXLOPER XL_intRule,
													  LPXLOPER XL_stubRule,
													  LPXLOPER XL_ccy,
													  LPXLOPER XL_payCal,
													  LPXLOPER XL_nxChange,
													  LPXLOPER XL_refDate,
													  LPXLOPER XL_adjStartDate,
													  LPXLOPER XL_rollDay)
{
	ADD_LOG("Local_FIXEDLEG ");
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
		bool ccyIsObject = false;
				
		CCString C_payCal;

		CCString C_nxChange;
		long     nxChange;

		double C_refDate;
		double C_refDate_default = -1.0;

		CCString C_adjStartDate;
		long     adjStartDateId;

		double	 C_rollDay;
		double   C_rollDay_default = -1.0;
		long     rollDay;

		// error
		static int error;
		static char* reason = "";

		XL_readNumCell(XL_startDate,C_startDate," ARM_ERR: start date: date expected",C_result);
		XL_readNumCell(XL_endDate,C_endDate," ARM_ERR: end date: date expected",C_result);
		XL_readStrCell(XL_receiveOrPay,C_receiveOrPay," ARM_ERR: receive or pay: string expected",C_result);
		XL_readStrOrNumCell(XL_fixRate, C_fixedRate_str, C_fixedRate_double, fixedRateType,
			   " ARM_ERR: fixed rate: numeric or object ID string expected",C_result);
		XL_readStrCellWD(XL_dayCount,C_dayCount,"30/360"," ARM_ERR: day count: string expected",C_result);
		XL_readStrCellWD(XL_freq,C_freq,"1"," ARM_ERR: frequency: string expected",C_result);
		XL_readStrCellWD(XL_decompFreq,C_decompFreq,"P"," ARM_ERR: decomp frequency: string expected",C_result);
		XL_readStrCellWD(XL_payTiming,C_payTiming,"ARR"," ARM_ERR: pay timing: string expected",C_result);
		XL_readStrCellWD(XL_intRule,C_intRule,"ADJ"," ARM_ERR: interpolation rule: string expected",C_result);
		XL_readStrCellWD(XL_stubRule,C_stubRule,"SS"," ARM_ERR: stub rule: string expected",C_result);
		XL_readStrCellWD(XL_ccy,C_ccy,"DEFAULT"," ARM_ERR: currency id: string expected",C_result);
		XL_readStrCellWD(XL_payCal,C_payCal,"NULL"," ARM_ERR: payment calendar name: string expected",C_result);
		XL_readStrCellWD(XL_nxChange,C_nxChange,"NXNONE"," ARM_ERR: Notional xchge: string expected",C_result);
		XL_readNumCellWD(XL_refDate,C_refDate,C_refDate_default," ARM_ERR: reference date: date expected",C_result);
		XL_readStrCellWD(XL_adjStartDate,C_adjStartDate,"YES"," ARM_ERR: adjStartDate: string expected",C_result);
		XL_readNumCellWD(XL_rollDay,C_rollDay, C_rollDay_default," ARM_ERR: rollDay: value between 1 and 31 expected",C_result);

		if (( C_ccy.GetLen() > 3 )
			&& 
			( !(C_ccy == "DEFAULT" ))
		    )
			ccyIsObject = true;

		if((receiveOrPayId = ARM_ConvRecOrPay (C_receiveOrPay, C_result)) == ARM_DEFAULT_ERR)
		{
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}

		dayCountId = ARM_ConvDayCount (C_dayCount);

		if((freqId = ARM_ConvFrequency (C_freq, C_result)) == ARM_DEFAULT_ERR)
		{
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}

		decompFreqId = ARM_ConvDecompFrequency (C_decompFreq);

		payTimingId = ARM_ConvPayResetRule (C_payTiming);

		intRuleId = ARM_ConvIntRule (C_intRule);

		stubRuleId = ARM_ConvStubRule (C_stubRule);
		
		if ( fixedRateType == XL_TYPE_STRING )
		{
			C_fixedRate_double = (double) LocalGetNumObjectId(C_fixedRate_str);

			fixedRateType = 1L;
		}
		else
		{
			fixedRateType = 0L;
		}

		nxChange = ARM_NotionalExchange(C_nxChange);

		if ((adjStartDateId = ARM_ConvYesOrNo(C_adjStartDate, C_result)) == ARM_DEFAULT_ERR)
		{
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}

		rollDay	=	(long) C_rollDay;

		long retCode;
		long objId;
		CCString prevClass;
		
		CCString curClass = LOCAL_SWAPLEG_CLASS;
		CCString stringId = GetLastCurCellEnvValue ();
		
		if(!stringId)
		{
			retCode = ARMLOCAL_FIXEDLEG (C_startDate,
										 C_endDate,
										 receiveOrPayId,
										 fixedRateType,
										 C_fixedRate_double,
										 dayCountId,
										 freqId,
										 decompFreqId,
										 payTimingId,
										 intRuleId,
										 stubRuleId,
										 ccyIsObject,
										 C_ccy,
										 C_payCal,
										 nxChange,
										 C_refDate,
										 adjStartDateId,
										 rollDay,
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
				retCode = ARMLOCAL_FIXEDLEG (C_startDate,
											 C_endDate,
											 receiveOrPayId,
											 fixedRateType,
											 C_fixedRate_double,
											 dayCountId,
											 freqId,
											 decompFreqId,
											 payTimingId,
											 intRuleId,
											 stubRuleId,
											 ccyIsObject,
											 C_ccy,
											 C_payCal,
											 nxChange,
											 C_refDate,
											 adjStartDateId,
											 rollDay,
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

				retCode = ARMLOCAL_FIXEDLEG (C_startDate,
											 C_endDate,
											 receiveOrPayId,
											 fixedRateType,
											 C_fixedRate_double,
											 dayCountId,
											 freqId,
											 decompFreqId,
											 payTimingId,
											 intRuleId,
											 stubRuleId,
											 ccyIsObject,
											 C_ccy,
											 C_payCal,
											 nxChange,
											 C_refDate,
											 adjStartDateId,
											 rollDay,
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_FIXEDLEG" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}



__declspec(dllexport) LPXLOPER WINAPI Local_PXL_FIXEDLEG (LPXLOPER XL_startDate,
														  LPXLOPER XL_endDate,
														  LPXLOPER XL_receiveOrPay,
														  LPXLOPER XL_fixRate,
														  LPXLOPER XL_dayCount,
														  LPXLOPER XL_freq,
														  LPXLOPER XL_decompFreq,
														  LPXLOPER XL_payTiming,
														  LPXLOPER XL_intRule,
														  LPXLOPER XL_stubRule,
														  LPXLOPER XL_ccy,
														  LPXLOPER XL_payCal,
														  LPXLOPER XL_nxChange,
														  LPXLOPER XL_refDate,
														  LPXLOPER XL_adjStartDate,
														  LPXLOPER XL_rollDay)
{
	ADD_LOG("Local_PXL_FIXEDLEG ");
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
		bool ccyIsObject = false;
				
		CCString C_payCal;

		CCString C_nxChange;
		long     nxChange;

		double C_refDate;
		double C_refDate_default = -1.0;

		CCString C_adjStartDate;
		long     adjStartDateId;

		double	 C_rollDay;
		double   C_rollDay_default = -1.0;
		long     rollDay;

		// error
		static int error;
		static char* reason = "";

		XL_readNumCell(XL_startDate,C_startDate," ARM_ERR: start date: date expected",C_result);
		XL_readNumCell(XL_endDate,C_endDate," ARM_ERR: end date: date expected",C_result);
		XL_readStrCell(XL_receiveOrPay,C_receiveOrPay," ARM_ERR: receive or pay: string expected",C_result);
		XL_readStrOrNumCell(XL_fixRate, C_fixedRate_str, C_fixedRate_double, fixedRateType,
			   " ARM_ERR: fixed rate: numeric or object ID string expected",C_result);
		XL_readStrCellWD(XL_dayCount,C_dayCount,"30/360"," ARM_ERR: day count: string expected",C_result);
		XL_readStrCellWD(XL_freq,C_freq,"1"," ARM_ERR: frequency: string expected",C_result);
		XL_readStrCellWD(XL_decompFreq,C_decompFreq,"P"," ARM_ERR: decomp frequency: string expected",C_result);
		XL_readStrCellWD(XL_payTiming,C_payTiming,"ARR"," ARM_ERR: pay timing: string expected",C_result);
		XL_readStrCellWD(XL_intRule,C_intRule,"ADJ"," ARM_ERR: interpolation rule: string expected",C_result);
		XL_readStrCellWD(XL_stubRule,C_stubRule,"SS"," ARM_ERR: stub rule: string expected",C_result);
		XL_readStrCellWD(XL_ccy,C_ccy,"DEFAULT"," ARM_ERR: currency id: string expected",C_result);
		XL_readStrCellWD(XL_payCal,C_payCal,"NULL"," ARM_ERR: payment calendar name: object expected",C_result);
		XL_readStrCellWD(XL_nxChange,C_nxChange,"NXNONE"," ARM_ERR: Notional xchge: string expected",C_result);
		XL_readNumCellWD(XL_refDate,C_refDate,C_refDate_default," ARM_ERR: reference date: date expected",C_result);
		XL_readStrCellWD(XL_adjStartDate,C_adjStartDate,"YES"," ARM_ERR: adjStartDate: string expected",C_result);
		XL_readNumCellWD(XL_rollDay,C_rollDay, C_rollDay_default," ARM_ERR: rollDay: value between 1 and 31 expected",C_result);

		if (( C_ccy.GetLen() > 3 )
			&& 
			(!( C_ccy == "DEFAULT" ))
		   )
		   ccyIsObject = true;

		if ( (receiveOrPayId = ARM_ConvRecOrPay (C_receiveOrPay, C_result)) == ARM_DEFAULT_ERR )
		{
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}

		dayCountId = ARM_ConvDayCount (C_dayCount);

		if ( (freqId = ARM_ConvFrequency (C_freq, C_result)) == ARM_DEFAULT_ERR )
		{
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}

		decompFreqId = ARM_ConvDecompFrequency (C_decompFreq);

		payTimingId = ARM_ConvPayResetRule (C_payTiming);

		intRuleId = ARM_ConvIntRule (C_intRule);

		stubRuleId = ARM_ConvStubRule (C_stubRule);

		if ( fixedRateType == XL_TYPE_STRING )
		{
			C_fixedRate_double = (double) LocalGetNumObjectId(C_fixedRate_str);

			fixedRateType = 1L;
		}
		else
		{
			fixedRateType = 0L;
		}

		nxChange = ARM_NotionalExchange(C_nxChange);

		if ((adjStartDateId = ARM_ConvYesOrNo(C_adjStartDate, C_result)) == ARM_DEFAULT_ERR)
		{
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}

		rollDay	=	(long) C_rollDay;

		long retCode;
		long objId;
		
		CCString curClass = LOCAL_SWAPLEG_CLASS;
		CCString stringId;
		
		retCode = ARMLOCAL_FIXEDLEG (C_startDate,
									 C_endDate,
									 receiveOrPayId,
									 fixedRateType,
									 C_fixedRate_double,
									 dayCountId,
									 freqId,
									 decompFreqId,
									 payTimingId,
									 intRuleId,
									 stubRuleId,
									 ccyIsObject,
									 C_ccy,
									 C_payCal,
									 nxChange,
									 C_refDate,
									 adjStartDateId,
									 rollDay,
									 C_result);

		if (retCode == ARM_OK)
		{
			objId = C_result.getLong ();

			stringId = LocalMakeObjectId (objId, curClass);
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_PXL_FIXEDLEG" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}



__declspec(dllexport) LPXLOPER WINAPI Local_SWAP_PRICE_TO_RATE (LPXLOPER XL_swap,
																LPXLOPER XL_date,
																LPXLOPER XL_price,
																LPXLOPER XL_model)
{
	ADD_LOG("Local_SWAP_PRICE_TO_RATE ");
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
		ARM_NOCALCIFWIZ();

		// C variable
		CCString C_swap;
		double C_date;
		double C_price;
		CCString C_model;
		
		// error
		static int error;
		static char* reason = "";

		XL_readStrCell(XL_swap,C_swap," ARM_ERR: swap id: object expected",C_result);
		XL_readNumCell(XL_date,C_date," ARM_ERR: date: date expected",C_result);
		XL_readNumCell(XL_price,C_price," ARM_ERR: price: numeric expected",C_result);
		XL_readStrCell(XL_model,C_model," ARM_ERR: model id: object expected",C_result);
		
		long retCode = ARMLOCAL_SWAP_PRICE_TO_RATE (LocalGetNumObjectId (C_swap), 
													C_date, 
													C_price, 
													LocalGetNumObjectId (C_model),
													C_result);

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
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_SWAP_PRICE_TO_RATE" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


__declspec(dllexport) LPXLOPER WINAPI Local_SWAP_RATE_TO_PRICE (LPXLOPER XL_swap,
																LPXLOPER XL_date,
																LPXLOPER XL_rate,
																LPXLOPER XL_model)
{
	ADD_LOG("Local_SWAP_RATE_TO_PRICE ");
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
		ARM_NOCALCIFWIZ();

		// C variable
		CCString C_swap;
		double C_date;
		double C_rate;
		CCString C_model;
		
		// error
		static int error;
		static char* reason = "";

		XL_readStrCell(XL_swap,C_swap," ARM_ERR: swap id: object expected",C_result);
		XL_readNumCell(XL_date,C_date," ARM_ERR: date: numeric expected",C_result);
		XL_readNumCell(XL_rate,C_rate," ARM_ERR: rate: numeric expected",C_result);
		XL_readStrCell(XL_model,C_model," ARM_ERR: model id: object expected",C_result);
		
		long retCode = ARMLOCAL_SWAP_RATE_TO_PRICE (LocalGetNumObjectId (C_swap), 
													C_date, 
													C_rate, 
													LocalGetNumObjectId (C_model), 
													C_result);

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
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_SWAP_RATE_TO_PRICE" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}



__declspec(dllexport) LPXLOPER WINAPI Local_SWAP (LPXLOPER XL_swLeg1,
												  LPXLOPER XL_swLeg2,
												  LPXLOPER XL_minPay,
												  LPXLOPER XL_fixedRates)
{
	ADD_LOG("Local_SWAP ");
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
		ARM_NOCALCIFWIZ();

		// C variable
		CCString C_swLeg1;
		CCString C_swLeg2;

		double C_minPay;


		double C_minPay_default = -1;

		VECTOR<double> C_fixedRates;
		VECTOR<double> C_fixedRates_default (0);
		
		// error
		static int error;
		static char* reason = "";

		XL_readStrCell(XL_swLeg1,C_swLeg1," ARM_ERR: swap leg 1 id: object expected",C_result);
		XL_readStrCell(XL_swLeg2,C_swLeg2," ARM_ERR: swap leg 2 id: object expected",C_result);
		XL_readNumCellWD(XL_minPay,C_minPay,C_minPay_default," ARM_ERR: minimum pay: numeric expected",C_result);
		XL_readNumVectorWD(XL_fixedRates,C_fixedRates,C_fixedRates_default," ARM_ERR: past fixed rates: array of numeric expected",C_result);
		
		long retCode;
		long objId;
		CCString prevClass;
		
		CCString curClass = LOCAL_SWAP_CLASS;
		CCString stringId = GetLastCurCellEnvValue ();
		
		if (!stringId)
		{
			retCode = ARMLOCAL_SWAP (LocalGetNumObjectId (C_swLeg1),
								LocalGetNumObjectId (C_swLeg2),
								(long)C_minPay, C_fixedRates, C_result);
									 
			if (retCode == ARM_OK)
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
				
			if (curClass == prevClass)
			{
				retCode = ARMLOCAL_SWAP (LocalGetNumObjectId (C_swLeg1),
									LocalGetNumObjectId (C_swLeg2),
									(long)C_minPay, C_fixedRates, C_result, objId);

				if (retCode == ARM_OK)
				{
					LocalSetCurCellEnvValue (curClass, objId); 

					stringId = LocalMakeObjectId (objId, curClass);
				}
			}
			else
			{
				FreeCurCellContent ();
				retCode = ARMLOCAL_SWAP (LocalGetNumObjectId (C_swLeg1),
									LocalGetNumObjectId (C_swLeg2),
									(long)C_minPay, C_fixedRates, C_result);
			
				if (retCode == ARM_OK)
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_SWAP" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}



__declspec(dllexport) LPXLOPER WINAPI Local_PXL_SWAP (LPXLOPER XL_swLeg1,
													  LPXLOPER XL_swLeg2,
													  LPXLOPER XL_minPay,
													  LPXLOPER XL_fixedRates)
{
	ADD_LOG("Local_PXL_SWAP ");
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
		ARM_NOCALCIFWIZ();

		// C variable
		CCString C_swLeg1;
		CCString C_swLeg2;

		double C_minPay;
		double C_minPay_default = -1;
		
		VECTOR<double> C_fixedRates;
		VECTOR<double> C_fixedRates_default (0);
		
		// error
		static int error;
		static char* reason = "";

		XL_readStrCell(XL_swLeg1,C_swLeg1," ARM_ERR: swap leg 1 id: object expected",C_result);
		XL_readStrCell(XL_swLeg2,C_swLeg2," ARM_ERR: swap leg 2 id: object expected",C_result);
		XL_readNumCellWD(XL_minPay,C_minPay,C_minPay_default," ARM_ERR: minimum pay: numeric expected",C_result);
		XL_readNumVectorWD(XL_fixedRates,C_fixedRates,C_fixedRates_default," ARM_ERR: past fixed rates: array of numeric expected",C_result);
		
		long retCode;
		long objId;
		
		CCString curClass = LOCAL_SWAP_CLASS;
		CCString stringId;
		
		retCode = ARMLOCAL_SWAP (LocalGetNumObjectId (C_swLeg1),
								LocalGetNumObjectId (C_swLeg2),
								(long)C_minPay, C_fixedRates, C_result);
									 
		if(retCode == ARM_OK)
		{
			objId = C_result.getLong ();

			stringId = LocalMakeObjectId (objId, curClass);
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_PXL_SWAP" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}

class getswaplegfromswapFunc : public ARMResultLong2LongFunc
{
public:
	getswaplegfromswapFunc( long swapId,
							int legNumber)
    :
    C_swapId(swapId),
    C_legNumber(legNumber)
    {};
	
	long operator()( ARM_result& result, long objId ){
		return ARMLOCAL_GETSWAPLEGFROMSWAP(
            C_swapId,
            C_legNumber,
            result,
            objId);
    }

private:
	long    C_swapId;
	int    C_legNumber;
};


LPXLOPER WINAPI Local_SwapFromExpiry (LPXLOPER XL_expiry,
									  LPXLOPER XL_tenor,
									  LPXLOPER XL_liborType,
									  LPXLOPER XL_fixedRate,
									  LPXLOPER XL_ccy)
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
		CCString C_expiry;
		CCString C_tenor;
		CCString C_liborType;
		long liborTypeId;

		CCString C_ccy;
		bool ccyIsObject = false;

		double C_fixedRate;
		double C_fixedRate_default = -1;

		// error
		static int error;
		static char* reason = "";

		XL_readStrCell(XL_expiry,C_expiry," ARM_ERR: expiry : string expected",C_result);
		XL_readStrCell(XL_tenor,C_tenor," ARM_ERR: tenor : string expected",C_result);
		XL_readStrCell(XL_liborType,C_liborType," ARM_ERR: libor type : string expected",C_result);
		XL_readNumCellWD(XL_fixedRate,C_fixedRate,C_fixedRate_default," ARM_ERR: fixed rate : double expected",C_result);
		XL_readStrCellWD(XL_ccy,C_ccy,"DEFAULT"," ARM_ERR: currency : string expected",C_result);

		long retCode;
		long objId;
		
		CCString prevClass;
		CCString curClass = LOCAL_SWAP_CLASS;
		CCString stringId = GetLastCurCellEnvValue ();
		
		if((liborTypeId = ARM_ConvIrIndName (C_liborType, C_result)) == ARM_DEFAULT_ERR)
		{
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}
		
		if ((C_ccy.GetLen() > 3)
		&&
		!(C_ccy == "DEFAULT"))
		ccyIsObject = true;


		if (!stringId)
		{
			retCode = ARMLOCAL_SwapFromExpiry (C_expiry,
											   C_tenor,
											   liborTypeId,
											   C_fixedRate,
											   ccyIsObject,
											   C_ccy,
											   C_result);
									 
			if (retCode == ARM_OK)
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
				
			if (curClass == prevClass)
			{
				retCode = ARMLOCAL_SwapFromExpiry (C_expiry,
												   C_tenor,
												   liborTypeId,
												   C_fixedRate,
												   ccyIsObject,
												   C_ccy,
												   C_result,
												   objId);

				if (retCode == ARM_OK)
				{
					LocalSetCurCellEnvValue (curClass, objId); 

					stringId = LocalMakeObjectId (objId, curClass);
				}
			}
			else
			{
				FreeCurCellContent ();
				retCode = ARMLOCAL_SwapFromExpiry (C_expiry,
												   C_tenor,
												   liborTypeId,
												   C_fixedRate,
												   ccyIsObject,
												   C_ccy,
												   C_result);
			
				if (retCode == ARM_OK)
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_SWAP" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
	
}

LPXLOPER WINAPI Local_PXL_SwapFromExpiry (LPXLOPER XL_expiry,
										  LPXLOPER XL_tenor,
										  LPXLOPER XL_liborType,
										  LPXLOPER XL_fixedRate,
										  LPXLOPER XL_ccy)
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
		CCString C_expiry;
		CCString C_tenor;
		CCString C_liborType;
		long liborTypeId;

		CCString C_ccy;
		bool ccyIsObject = false;

		double C_fixedRate;
		double C_fixedRate_default = -1;

		// error
		static int error;
		static char* reason = "";

		XL_readStrCell(XL_expiry,C_expiry," ARM_ERR: expiry : string expected",C_result);
		XL_readStrCell(XL_tenor,C_tenor," ARM_ERR: tenor : string expected",C_result);
		XL_readStrCell(XL_liborType,C_liborType," ARM_ERR: libor type : string expected",C_result);
		XL_readNumCellWD(XL_fixedRate,C_fixedRate,C_fixedRate_default," ARM_ERR: fixed rate : double expected",C_result);
		XL_readStrCellWD(XL_ccy,C_ccy,"DEFAULT"," ARM_ERR: currency : string expected",C_result);
		
		long retCode;
		long objId;
		CCString prevClass;
		
		CCString curClass = LOCAL_SWAP_CLASS;
		CCString stringId;
		
		if((liborTypeId = ARM_ConvIrIndName (C_liborType, C_result)) == ARM_DEFAULT_ERR)
		{
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}
		
		if ((C_ccy.GetLen() > 3)
		&&
		!(C_ccy == "DEFAULT"))
		ccyIsObject = true;

		retCode = ARMLOCAL_SwapFromExpiry (C_expiry,
										   C_tenor,
										   liborTypeId,
										   C_fixedRate,
										   ccyIsObject,
										   C_ccy,
										   C_result);
		
		if (retCode == ARM_OK)
		{
			objId = C_result.getLong ();

			LocalSetCurCellEnvValue (curClass, objId); 

			stringId = LocalMakeObjectId (objId, curClass);
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_SWAP" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


LPXLOPER Local_GetSwapLegFromSwap_Common(
	LPXLOPER XL_swapId,
	LPXLOPER XL_legNumber,
	bool PersistentInXL )
{
	
	/// this is defined first because it is used in XL macros
	static XLOPER XL_result;
	
	/// Get the variables from the XLOper variables
	ARM_result C_result;
	
	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
		/// to avoid computation if called by the wizard
		ARM_NOCALCIFWIZ();
		
		/// this is used by macros 
		/// and therefore this has to be defined
		static int error;
		static char* reason = "";

		CCString swapStrId;
		XL_readStrCell(XL_swapId, swapStrId,"ARM_ERR: swap leg 1 id: object expected",C_result);

		double legNumber;
		XL_readNumCell(XL_legNumber, legNumber, "ARM_ERR: Leg Number: long 1 or 2 expected",C_result);
	
		getswaplegfromswapFunc ourFunc(LocalGetNumObjectId(swapStrId), (int)legNumber);

		fillXL_Result( LOCAL_SWAPLEG_CLASS, ourFunc, C_result, XL_result, PersistentInXL );
	}

	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_BermudaSwaptionGet_Common" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}



__declspec(dllexport) LPXLOPER WINAPI Local_GetSwapLegFromSwap(
	LPXLOPER XL_swapId,
	LPXLOPER XL_legNumber)
{
	ADD_LOG("Local_GetSwapLegFromSwap");
	bool PersistentInXL = true;
	return Local_GetSwapLegFromSwap_Common(
	    XL_swapId,
	    XL_legNumber,
        PersistentInXL );
}

_declspec(dllexport) LPXLOPER WINAPI Local_PXL_GetSwapLegFromSwap(
	LPXLOPER XL_swapId,
	LPXLOPER XL_legNumber)
{
	bool PersistentInXL = false;
	return Local_GetSwapLegFromSwap_Common(
	    XL_swapId,
	    XL_legNumber,
        PersistentInXL );
}

__declspec(dllexport) LPXLOPER WINAPI Local_CMSLEG (LPXLOPER XL_startDate,
													LPXLOPER XL_endDate,
													LPXLOPER XL_cmsType,
													LPXLOPER XL_receiveOrPay,
													LPXLOPER XL_spread,
													LPXLOPER XL_yieldDecompFreq,
													LPXLOPER XL_swapLegDayCount,
													LPXLOPER XL_resetFreq,
													LPXLOPER XL_intRule,
													LPXLOPER XL_ccy,
													LPXLOPER XL_resetTiming,
                                                    LPXLOPER XL_payFreq,
													LPXLOPER XL_resetGap,
													LPXLOPER XL_stubRule,
													LPXLOPER XL_adjStartDate)
{
	ADD_LOG("Local_CMSLEG ");
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

		CCString C_cmsType;
		long cmsTypeId;

		CCString C_receiveOrPay;
		long receiveOrPayId;

		double C_spread_double;
		CCString C_spread_str;
		long spreadType;
		double C_default_spread = 0.0;
		
		CCString C_yieldDecompFreq;
		long yieldDecompFreqId;

		CCString C_swapLegDayCount;
		long swapLegDayCountId;
		
		CCString C_resetFreq;
		long resetFreqId;

		CCString C_payFreq;
		long payFreqId;
		
		double C_resetGap;
		double C_resetGap_default = 10000.;

		CCString C_intRule;
		long intRuleId;

		CCString C_ccy;
		bool ccyIsObject = false;
				
		CCString C_resetTiming;
		long resetTiming;

		CCString C_stubRule;
		long stubRule;

		CCString C_adjStartDate;
		long adjStartDate = 1;

		// error
		static int error;
		static char* reason = "";

		XL_readNumCell(XL_startDate,C_startDate," ARM_ERR: start date: date expected",C_result);
		XL_readNumCell(XL_endDate,C_endDate," ARM_ERR: end date: date expected",C_result);
		XL_readStrCell(XL_cmsType,C_cmsType," ARM_ERR: cms type: string expected",C_result);
		XL_readStrCell(XL_receiveOrPay,C_receiveOrPay," ARM_ERR: receive or pay: string expected",C_result);
		XL_readStrOrNumCellWD(XL_spread,C_spread_str,C_spread_double,C_default_spread,spreadType," ARM_ERR: spread: numeric or object expected",C_result);
		XL_readStrCellWD(XL_yieldDecompFreq,C_yieldDecompFreq,"-1"," ARM_ERR: yield decomp frequency: string expected",C_result);
		XL_readStrCellWD(XL_swapLegDayCount,C_swapLegDayCount,"A360"," ARM_ERR: swap leg day count: string expected",C_result);
		XL_readStrCellWD(XL_resetFreq,C_resetFreq,"Q"," ARM_ERR: reset frequency: string expected",C_result);
		XL_readStrCellWD(XL_intRule,C_intRule,"ADJ"," ARM_ERR: interpolation rule: string expected",C_result);
		XL_readStrCellWD(XL_ccy,C_ccy,"DEFAULT"," ARM_ERR: currency id: string expected",C_result);
		XL_readStrCellWD(XL_resetTiming,C_resetTiming,"ADV"," ARM_ERR: pay timing: string expected",C_result);
		XL_readStrCellWD(XL_payFreq,C_payFreq,"DEFAULT"," ARM_ERR: pay frequency: string expected",C_result);
		XL_readNumCellWD(XL_resetGap,C_resetGap,C_resetGap_default," ARM_ERR: reset gap: numeric expected",C_result);
		XL_readStrCellWD(XL_stubRule,C_stubRule,"SS"," ARM_ERR: stub rule: string expected",C_result);
		XL_readStrCellWD(XL_adjStartDate, C_adjStartDate, "YES", " ARM_ERR: adjust Start Date: YES or NO expected", C_result);

		if((XL_spread->xltype == xltypeMissing) || (XL_spread->xltype == xltypeNil))
		{
			spreadType = 0;
		}
		else
		if ( spreadType == XL_TYPE_STRING )
		{
		   C_spread_double = (double) LocalGetNumObjectId (C_spread_str);

		   spreadType = 1;
		}
		else
		{
		   spreadType = 0;
		}

		cmsTypeId = ARM_ConvIrType (C_cmsType);

		if((receiveOrPayId = ARM_ConvRecOrPay (C_receiveOrPay, C_result)) == ARM_DEFAULT_ERR)
		{
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}

		yieldDecompFreqId = ARM_ConvDecompFrequency (C_yieldDecompFreq);

		swapLegDayCountId = ARM_ConvDayCount (C_swapLegDayCount);
		
		intRuleId = ARM_ConvIntRule (C_intRule);

		if((resetFreqId = ARM_ConvFrequency (C_resetFreq, C_result)) == ARM_DEFAULT_ERR)
		{
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}
		
		if(C_payFreq == "DEFAULT")
		{
			C_payFreq = C_resetFreq;
		}

		if((payFreqId = ARM_ConvFrequency (C_payFreq, C_result)) == ARM_DEFAULT_ERR)
		{
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}

		if (( C_ccy.GetLen() > 3 )
			&&
			( !(C_ccy == "DEFAULT") )
		   )
		   ccyIsObject = true;

		resetTiming = ARM_ConvPayResetRule (C_resetTiming);

		stubRule = ARM_ConvStubRule (C_stubRule);

		if( C_adjStartDate != CCString("YES") )
		{
			adjStartDate = 0;
		}

		long retCode;
		long objId;
		CCString prevClass;
		
		CCString curClass = LOCAL_CMSLEG_CLASS;
		CCString stringId = GetLastCurCellEnvValue ();
		
		if(!stringId)
		{
			retCode = ARMLOCAL_CMSLEG (C_startDate,
									   C_endDate,
									   cmsTypeId,
									   receiveOrPayId,
									   spreadType,
									   C_spread_double,
									   yieldDecompFreqId,
									   swapLegDayCountId,
									   resetFreqId,
									   payFreqId,
									   (long) C_resetGap,
									   intRuleId,
									   ccyIsObject,
									   C_ccy,
									   resetTiming,
									   stubRule,
									   adjStartDate,
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
				retCode = ARMLOCAL_CMSLEG (C_startDate,
										   C_endDate,
										   cmsTypeId,
										   receiveOrPayId,
										   spreadType,
										   C_spread_double,
										   yieldDecompFreqId,
										   swapLegDayCountId,
										   resetFreqId,
										   payFreqId,
										   (long) C_resetGap,
										   intRuleId,
										   ccyIsObject,
										   C_ccy,
										   resetTiming,
										   stubRule,
										   adjStartDate,
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
				retCode = ARMLOCAL_CMSLEG (C_startDate,
										   C_endDate,
										   cmsTypeId,
										   receiveOrPayId,
										   spreadType,
										   C_spread_double,
										   yieldDecompFreqId,
										   swapLegDayCountId,
										   resetFreqId,
										   payFreqId,
										   (long) C_resetGap,
										   intRuleId,
										   ccyIsObject,
										   C_ccy,
										   resetTiming,
										   stubRule,
										   adjStartDate,
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_CMSLEG" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


__declspec(dllexport) LPXLOPER WINAPI Local_PXL_CMSLEG (LPXLOPER XL_startDate,
														LPXLOPER XL_endDate,
														LPXLOPER XL_cmsType,
														LPXLOPER XL_receiveOrPay,
														LPXLOPER XL_spread,
														LPXLOPER XL_yieldDecompFreq,
														LPXLOPER XL_swapLegDayCount,
														LPXLOPER XL_resetFreq,
														LPXLOPER XL_intRule,
														LPXLOPER XL_ccy,
														LPXLOPER XL_resetTiming,
														LPXLOPER XL_payFreq,
														LPXLOPER XL_resetGap,
														LPXLOPER XL_stubRule,
														LPXLOPER XL_adjStartDate)
{
	ADD_LOG("Local_PXL_CMSLEG ");
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

		CCString C_cmsType;
		long cmsTypeId;

		CCString C_receiveOrPay;
		long receiveOrPayId;

		double C_spread_double;
		CCString C_spread_str;
		long spreadType;
		double C_default_spread = 0.0;
		
		CCString C_yieldDecompFreq;
		long yieldDecompFreqId;

		CCString C_swapLegDayCount;
		long swapLegDayCountId;
		
		CCString C_resetFreq;
		long resetFreqId;

		CCString C_payFreq;
		long payFreqId;
		
		double C_resetGap;
		double C_resetGap_default = 10000.;

		CCString C_intRule;
		long intRuleId;

		CCString C_ccy;
		bool ccyIsObject = false;
				
		CCString C_resetTiming;
		long resetTiming;

		CCString C_stubRule;
		long stubRule;

		CCString C_adjStartDate;
		long adjStartDate = 1;

		// error
		static int error;
		static char* reason = "";

		XL_readNumCell(XL_startDate,C_startDate," ARM_ERR: start date: date expected",C_result);
		XL_readNumCell(XL_endDate,C_endDate," ARM_ERR: end date: date expected",C_result);
		XL_readStrCell(XL_cmsType,C_cmsType," ARM_ERR: cms type: string expected",C_result);
		XL_readStrCell(XL_receiveOrPay,C_receiveOrPay," ARM_ERR: receive or pay: string expected",C_result);
		XL_readStrOrNumCellWD(XL_spread,C_spread_str,C_spread_double,C_default_spread,spreadType," ARM_ERR: spread: numeric or object expected",C_result);
		XL_readStrCellWD(XL_yieldDecompFreq,C_yieldDecompFreq,"-1"," ARM_ERR: yield decomp frequency: string expected",C_result);
		XL_readStrCellWD(XL_swapLegDayCount,C_swapLegDayCount,"A360"," ARM_ERR: swap leg day count: string expected",C_result);
		XL_readStrCellWD(XL_resetFreq,C_resetFreq,"Q"," ARM_ERR: reset frequency: string expected",C_result);
		XL_readStrCellWD(XL_intRule,C_intRule,"ADJ"," ARM_ERR: interpolation rule: string expected",C_result);
		XL_readStrCellWD(XL_ccy,C_ccy,"DEFAULT"," ARM_ERR: currency id: string expected",C_result);
		XL_readStrCellWD(XL_resetTiming,C_resetTiming,"ADV"," ARM_ERR: pay timing: string expected",C_result);
		XL_readStrCellWD(XL_payFreq,C_payFreq,"DEFAULT"," ARM_ERR: pay frequency: string expected",C_result);
		XL_readNumCellWD(XL_resetGap,C_resetGap,C_resetGap_default," ARM_ERR: reset gap: numeric expected",C_result);
		XL_readStrCellWD(XL_stubRule,C_stubRule,"SS"," ARM_ERR: stub rule: string expected",C_result);
		XL_readStrCellWD(XL_adjStartDate, C_adjStartDate, "YES", " ARM_ERR: adjust Start Date: YES or NO expected", C_result);

		if((XL_spread->xltype == xltypeMissing) || (XL_spread->xltype == xltypeNil))
		{
			spreadType = 0;
		}
		else if ( spreadType == XL_TYPE_STRING )
		{
		   C_spread_double = (double) LocalGetNumObjectId (C_spread_str);

		   spreadType = 1;
		}
		else
		{
		   spreadType = 0;
		}

		cmsTypeId = ARM_ConvIrType (C_cmsType);

		if((receiveOrPayId = ARM_ConvRecOrPay (C_receiveOrPay, C_result)) == ARM_DEFAULT_ERR)
		{
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}

		yieldDecompFreqId = ARM_ConvDecompFrequency (C_yieldDecompFreq);

		swapLegDayCountId = ARM_ConvDayCount (C_swapLegDayCount);

		intRuleId = ARM_ConvIntRule (C_intRule);

		if((resetFreqId = ARM_ConvFrequency (C_resetFreq, C_result)) == ARM_DEFAULT_ERR)
		{
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}

		if (C_payFreq == "DEFAULT")
		{
			C_payFreq = C_resetFreq;
		}
		if((payFreqId = ARM_ConvFrequency (C_payFreq, C_result)) == ARM_DEFAULT_ERR)
		{
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}

		if ((C_ccy.GetLen() > 3)
			&&
			!(C_ccy == "DEFAULT")
		   )
		   ccyIsObject = true;

		resetTiming = ARM_ConvPayResetRule (C_resetTiming);
		stubRule = ARM_ConvStubRule (C_stubRule);
		
		long retCode;
		long objId;
		
		CCString curClass = LOCAL_CMSLEG_CLASS;
		CCString stringId;
		
		retCode = ARMLOCAL_CMSLEG (C_startDate,
								   C_endDate,
								   cmsTypeId,
								   receiveOrPayId,
								   spreadType,
								   C_spread_double,
								   yieldDecompFreqId,
								   swapLegDayCountId,
								   resetFreqId,
								   payFreqId,
								   (long) C_resetGap,
								   intRuleId,
								   ccyIsObject,
								   C_ccy,
								   resetTiming,
								   stubRule,
								   adjStartDate,
								   C_result);

		if(retCode == ARM_OK)
		{
			objId = C_result.getLong ();

			stringId = LocalMakeObjectId (objId, curClass);
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_PXL_CMSLEG" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}




__declspec(dllexport) LPXLOPER WINAPI Local_CORRIDORLEG(LPXLOPER XL_startDate,
														LPXLOPER XL_endDate,
														LPXLOPER XL_receiveOrPay,
														LPXLOPER XL_payIndexId,
														LPXLOPER XL_payFreq,
														LPXLOPER XL_spread,
														LPXLOPER XL_refIndexId,
														LPXLOPER XL_resetFreq,
														LPXLOPER XL_paidRateResetTiming,
														LPXLOPER XL_refRateResetTiming,
														LPXLOPER XL_stubRule,
														LPXLOPER XL_levelDownId,
														LPXLOPER XL_downSpec,
														LPXLOPER XL_levelUpId,
														LPXLOPER XL_upSpec,
														LPXLOPER XL_ccy,
														LPXLOPER XL_MCFreq,
														LPXLOPER XL_MCInterp,
														LPXLOPER XL_refIndexResetGap,
                                                        LPXLOPER XL_decompPricingFlag)
{
	ADD_LOG("Local_CORRIDORLEG");
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
		long receiveOrPay;

		CCString C_payIndexId;
		long payIndexId;

		CCString C_payFreq;
		long payFreq;

		CCString C_spread_str;
		double C_spread_double;
		double C_spread_default = 0.0;
		long spreadType;

		CCString C_refIndexId;
		long refIndexId;

		CCString C_resetFreq;
		long resetFreq;

		CCString C_paidRateResetTiming;
		long paidRateResetTiming;

		CCString C_refRateResetTiming;
		long refRateResetTiming;
		
		CCString C_stubRule;
		long stubRule;

		CCString C_levelDownId;
		long levelDownId;
    						  
		CCString C_downSpec;
		long     downSpec;

		CCString C_levelUpId;
		long levelUpId;

		CCString C_upSpec;
		long     upSpec;

		CCString C_ccy;
		bool ccyIsObject = false;

		CCString C_MCFreq;
		long MCFreq;

		CCString C_MCInterp;
		long MCInterp;

		CCString C_LDPricingMethod("DIG");
		long LDPricingMethodId;

		double decompPricingFlag;
		double decompPricingFlag_default = 1.0;

        double C_RefIndexResetGap = 0.0;;
        double C_RefIndexresetGap_default = 0.0;

		// error
		static int error;
		static char* reason = "";

		XL_readNumCell(XL_startDate,C_startDate," ARM_ERR: start date: date expected",C_result);
		XL_readNumCell(XL_endDate,C_endDate," ARM_ERR: end date: date expected",C_result);
		XL_readStrCell(XL_receiveOrPay,C_receiveOrPay," ARM_ERR: receive or pay: string expected",C_result);    
		XL_readStrCell(XL_payIndexId,C_payIndexId," ARM_ERR: Object paymentIndex Id expected",C_result);
		XL_readStrCellWD(XL_payFreq,C_payFreq,"-1"," ARM_ERR: frequency: string expected",C_result);
		XL_readStrOrNumCellWD(XL_spread,C_spread_str,C_spread_double,C_spread_default,spreadType," ARM_ERR: spread: numeric or object expected",C_result);
		XL_readStrCell(XL_refIndexId,C_refIndexId," ARM_ERR: Object refindex Id expected",C_result);
		XL_readStrCellWD(XL_resetFreq,C_resetFreq,"-1"," ARM_ERR: frequency: string expected",C_result);
		XL_readStrCellWD(XL_paidRateResetTiming,C_paidRateResetTiming,"ADV"," ARM_ERR: paidRateResetTiming: string expected",C_result);
		XL_readStrCellWD(XL_refRateResetTiming,C_refRateResetTiming,"ADV"," ARM_ERR: refRateResetTiming: string expected",C_result);
		XL_readStrCellWD(XL_stubRule,C_stubRule,"SS"," ARM_ERR: stub rule: string expected",C_result);
		XL_readStrCellWD(XL_levelDownId,C_levelDownId,"DEFAULT"," ARM_ERR: Object refValue expected",C_result);
		XL_readStrCellWD(XL_downSpec,C_downSpec,"STD"," ARM_ERR: DownSpec: String expected",C_result);
		XL_readStrCellWD(XL_levelUpId,C_levelUpId,"DEFAULT"," ARM_ERR: Object RefValue expected",C_result);
		XL_readStrCellWD(XL_upSpec,C_upSpec,"STD"," ARM_ERR: DownSpec: String expected",C_result);
		XL_readStrCellWD(XL_ccy,C_ccy,"DEFAULT"," ARM_ERR: currency id: string expected",C_result);
		XL_readStrCellWD(XL_MCFreq,C_MCFreq,"-1"," ARM_ERR: frequency: string expected",C_result);
		
        XL_readStrCellWD(XL_MCInterp,C_MCInterp,"LIN"," ARM_ERR: Interpolation method: string expected",C_result);
		
        
        XL_readNumCellWD(XL_refIndexResetGap, C_RefIndexResetGap, 
                         C_RefIndexresetGap_default," ARM_ERR: Ref Index reset gap: numeric expected", C_result);
		
        XL_readNumCellWD(XL_decompPricingFlag,decompPricingFlag,decompPricingFlag_default," ARM_ERR: decompPricingFlag: string expected",C_result);
	
        

		if ((receiveOrPay = ARM_ConvRecOrPay(C_receiveOrPay, C_result)) == ARM_DEFAULT_ERR)
		{
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}

		payIndexId = LocalGetNumObjectId(C_payIndexId);
		refIndexId = LocalGetNumObjectId(C_refIndexId);

		if ( C_payFreq == "-1" )
		{
		   payFreq = K_DEF_FREQ;
		}
		else
		{
			if ((payFreq = ARM_ConvFrequency(C_payFreq, C_result)) == ARM_DEFAULT_ERR)
			{
	  			ARM_ARG_ERR();
				return (LPXLOPER)&XL_result;
			}
		}

		if ( spreadType == XL_TYPE_STRING )
		{
			C_spread_double = (double) LocalGetNumObjectId(C_spread_str);

			spreadType = 1L;
		}
		else
		{
			spreadType = 0L;
		}
			
		if ( C_resetFreq == "-1" )
		{
			resetFreq = K_DAILY;
		}
		else
		{
			if ((resetFreq = ARM_ConvFrequency(C_resetFreq, C_result)) == ARM_DEFAULT_ERR)
			{
	  			ARM_ARG_ERR();
				return (LPXLOPER)&XL_result;
			}
		}

		paidRateResetTiming = ARM_ConvPayResetRule (C_paidRateResetTiming);

		refRateResetTiming = ARM_ConvPayResetRule (C_refRateResetTiming);

		stubRule = ARM_ConvStubRule (C_stubRule);

		if ( C_levelDownId == "DEFAULT" )
		{
			levelDownId = ARM_NULL_OBJECT;
		}
		else
		{
			levelDownId = LocalGetNumObjectId(C_levelDownId);
		}

		if ( ( downSpec = ARM_ConvStdBoost(C_downSpec, C_result) ) == ARM_DEFAULT_ERR )
		{
			ARM_ARG_ERR();
			return((LPXLOPER) &XL_result);
		}

		if ( C_levelUpId == "DEFAULT" )
		{
			levelUpId = ARM_NULL_OBJECT;
		}
		else
		{
			levelUpId = LocalGetNumObjectId(C_levelUpId);
		}

		if (( upSpec = ARM_ConvStdBoost(C_upSpec, C_result) ) == ARM_DEFAULT_ERR )
		{
			ARM_ARG_ERR();
			return((LPXLOPER) &XL_result);
		}

		if (( C_ccy.GetLen() > 3 ) 
			&& 
			( !(C_ccy == "DEFAULT") )
		   )
		   ccyIsObject = true;

		if ( C_MCFreq == "-1" )
		{
			MCFreq = K_DAILY;
		}
		else
		{
			if ((MCFreq = ARM_ConvFrequency(C_MCFreq, C_result)) == ARM_DEFAULT_ERR)
			{
	  			ARM_ARG_ERR();
				return (LPXLOPER)&XL_result;
		   }
		}

		MCInterp = K_LINEAR;
		if ( C_MCInterp == "LOG" )
		{
			MCInterp = K_CONTINUOUS;
		}
		else if ( C_MCInterp == "RAC" )
		{
			MCInterp = K_SQRT;
		}
		else if ( C_MCInterp == "SL_LIN" )
		{
			MCInterp = K_SLOPELIN;
		}
		else if ( C_MCInterp == "SL_RAC" )
		{
			MCInterp = K_SLOPESQRT;
		}

		LDPricingMethodId = ARM_ConvPricCorridorLD(C_LDPricingMethod) ;

		long retCode;
		long objId;
		CCString prevClass;
		
		CCString curClass = LOCAL_CORRIDORLEG_CLASS;
		CCString stringId = GetLastCurCellEnvValue ();

		if (!stringId)
		{
			retCode = ARMLOCAL_CORRIDORLEG(C_startDate,
										   C_endDate,
										   receiveOrPay,
										   payIndexId,
										   payFreq,
										   spreadType,
										   C_spread_double,
										   refIndexId,
										   resetFreq,
										   paidRateResetTiming,
										   refRateResetTiming,
										   stubRule,
										   levelDownId,
										   downSpec,
										   levelUpId,
										   upSpec,
										   ccyIsObject,
										   C_ccy,
										   MCFreq,
										   MCInterp,
										   LDPricingMethodId,
										   long(decompPricingFlag),
										   "NULL", // reset calendar
										   "NULL", // pay calendar
                                           long(C_RefIndexResetGap),
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
				
			if ( curClass == prevClass )
			{
				retCode = ARMLOCAL_CORRIDORLEG(C_startDate,
											   C_endDate,
											   receiveOrPay,
											   payIndexId,
											   payFreq,
											   spreadType,
											   C_spread_double,
											   refIndexId,
											   resetFreq,
											   paidRateResetTiming,
											   refRateResetTiming,
											   stubRule,
											   levelDownId,
											   downSpec,
											   levelUpId,
											   upSpec,
											   ccyIsObject,
											   C_ccy,
											   MCFreq,
											   MCInterp,
											   LDPricingMethodId,
											   long(decompPricingFlag),
											   "NULL", // reset calendar
											   "NULL", // pay calendar
                                               long(C_RefIndexResetGap),
											   C_result,
											   objId);

				if ( retCode == ARM_OK )
				{
					LocalSetCurCellEnvValue (curClass, objId); 

					stringId = LocalMakeObjectId (objId, curClass);
				}
			}
			else
			{
				FreeCurCellContent();

				retCode = ARMLOCAL_CORRIDORLEG(C_startDate,
											   C_endDate,
											   receiveOrPay,
											   payIndexId,
											   payFreq,
											   spreadType,
											   C_spread_double,
											   refIndexId,
											   resetFreq,
											   paidRateResetTiming,
											   refRateResetTiming,
											   stubRule,
											   levelDownId,
											   downSpec,
											   levelUpId,
											   upSpec,
											   ccyIsObject,
											   C_ccy,
											   MCFreq,
											   MCInterp,
											   LDPricingMethodId,
											   long(decompPricingFlag),
											   "NULL", // reset calendar
											   "NULL", // pay calendar
                                               long(C_RefIndexResetGap),
											   C_result);

				if ( retCode == ARM_OK )
				{
					objId = C_result.getLong ();
				
					LocalSetCurCellEnvValue (curClass, objId); 

					stringId = LocalMakeObjectId (objId, curClass);
				}
			}
		}

		if ( retCode == ARM_OK )
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_CORRIDORLEG" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}




__declspec(dllexport) LPXLOPER WINAPI Local_PXL_CORRIDORLEG(LPXLOPER XL_startDate,
															LPXLOPER XL_endDate,
															LPXLOPER XL_receiveOrPay,
															LPXLOPER XL_payIndexId,
															LPXLOPER XL_payFreq,
															LPXLOPER XL_spread,
															LPXLOPER XL_refIndexId,
															LPXLOPER XL_resetFreq,
															LPXLOPER XL_paidRateResetTiming,
															LPXLOPER XL_refRateResetTiming,
															LPXLOPER XL_stubRule,
															LPXLOPER XL_levelDownId,
															LPXLOPER XL_downSpec,
															LPXLOPER XL_levelUpId,
															LPXLOPER XL_upSpec,
															LPXLOPER XL_ccy,
															LPXLOPER XL_MCFreq,
															LPXLOPER XL_MCInterp,
														    LPXLOPER XL_refIndexResetGap,
                                                            LPXLOPER XL_decompPricingFlag)
{
	ADD_LOG("Local_PXL_CORRIDORLEG");
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
		long receiveOrPay;

		CCString C_payIndexId;
		long payIndexId;

		CCString C_payFreq;
		long payFreq;

		CCString C_spread_str;
		double C_spread_double;
		double C_spread_default = 0.0;
		long spreadType;

		CCString C_refIndexId;
		long refIndexId;

		CCString C_resetFreq;
		long resetFreq;

		CCString C_paidRateResetTiming;
		long paidRateResetTiming;

		CCString C_refRateResetTiming;
		long refRateResetTiming;
		
		CCString C_stubRule;
		long stubRule;

		CCString C_levelDownId;
		long levelDownId;
    						  
		CCString C_downSpec;
		long     downSpec;

		CCString C_levelUpId;
		long levelUpId;

		CCString C_upSpec;
		long     upSpec;

		CCString C_ccy;
		bool ccyIsObject = false;

		CCString C_MCFreq;
		long MCFreq;

		CCString C_MCInterp;
		long MCInterp;

		CCString C_LDPricingMethod("DIG");
		long LDPricingMethodId;

		double decompPricingFlag;
		double decompPricingFlag_default = 1.0;

        double C_RefIndexResetGap = 0.0;
        double C_RefIndexresetGap_default = 0.0;

		// error
		static int error;
		static char* reason = "";

		XL_readNumCell(XL_startDate,C_startDate," ARM_ERR: start date: date expected",C_result);
		XL_readNumCell(XL_endDate,C_endDate," ARM_ERR: end date: date expected",C_result);
		XL_readStrCell(XL_receiveOrPay,C_receiveOrPay," ARM_ERR: receive or pay: string expected",C_result);    
		XL_readStrCell(XL_payIndexId,C_payIndexId, "ARM_ERR: Object PaymentIndex Id expected",C_result);
		XL_readStrCellWD(XL_payFreq,C_payFreq,"-1"," ARM_ERR: frequency: string expected",C_result);
		XL_readStrOrNumCellWD(XL_spread,C_spread_str,C_spread_double,C_spread_default,spreadType," ARM_ERR: spread: numeric or object expected",C_result);
		XL_readStrCell(XL_refIndexId,C_refIndexId," ARM_ERR: Object refIdex Id expected",C_result);
		XL_readStrCellWD(XL_resetFreq,C_resetFreq,"-1"," ARM_ERR: frequency: string expected",C_result);
		XL_readStrCellWD(XL_paidRateResetTiming,C_paidRateResetTiming,"ADV"," ARM_ERR: paidRateResetTiming: string expected",C_result);
		XL_readStrCellWD(XL_refRateResetTiming,C_refRateResetTiming,"ADV"," ARM_ERR: refRateResetTiming: string expected",C_result);
		XL_readStrCellWD(XL_stubRule,C_stubRule,"SS"," ARM_ERR: stub rule: string expected",C_result);
		XL_readStrCellWD(XL_levelDownId,C_levelDownId,"DEFAULT"," ARM_ERR: Object RefValue expected",C_result);
		XL_readStrCellWD(XL_downSpec,C_downSpec,"STD"," ARM_ERR: DownSpec: String expected",C_result);
		XL_readStrCellWD(XL_levelUpId,C_levelUpId,"DEFAULT"," ARM_ERR: Object RefValue expected",C_result);
		XL_readStrCellWD(XL_upSpec,C_upSpec,"STD"," ARM_ERR: DownSpec: String expected",C_result);
		XL_readStrCellWD(XL_ccy,C_ccy,"DEFAULT"," ARM_ERR: currency id: string expected",C_result);
		XL_readStrCellWD(XL_MCFreq,C_MCFreq,"-1"," ARM_ERR: frequency: string expected",C_result);
		XL_readStrCellWD(XL_MCInterp,C_MCInterp,"LIN"," ARM_ERR: Interpolation method: string expected",C_result);
		
        XL_readNumCellWD(XL_refIndexResetGap, C_RefIndexResetGap, 
                         C_RefIndexresetGap_default," ARM_ERR: Ref Index reset gap: numeric expected", C_result);
		
        XL_readNumCellWD(XL_decompPricingFlag,decompPricingFlag,decompPricingFlag_default," ARM_ERR: decompPricingFlag: string expected",C_result);
		
        

		if ((receiveOrPay = ARM_ConvRecOrPay(C_receiveOrPay, C_result)) == ARM_DEFAULT_ERR)
		{
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}

		payIndexId = LocalGetNumObjectId(C_payIndexId);
		refIndexId = LocalGetNumObjectId(C_refIndexId);

		if ( C_payFreq == "-1" )
		{
		   payFreq = K_DEF_FREQ;
		}
		else
		{
			if ((payFreq = ARM_ConvFrequency(C_payFreq, C_result)) == ARM_DEFAULT_ERR)
			{
	  			ARM_ARG_ERR();
				return (LPXLOPER)&XL_result;
			}
		}

		if ( spreadType == XL_TYPE_STRING )
		{
			C_spread_double = (double) LocalGetNumObjectId(C_spread_str);

			spreadType = 1L;
		}
		else
		{
			spreadType = 0L;
		}
    

		if ( C_resetFreq == "-1" )
		{
			resetFreq = K_DAILY;
		}
		else
		{
			if ((resetFreq = ARM_ConvFrequency(C_resetFreq, C_result)) == ARM_DEFAULT_ERR)
			{
	  			ARM_ARG_ERR();
				return (LPXLOPER)&XL_result;
			}
		}

		paidRateResetTiming = ARM_ConvPayResetRule (C_paidRateResetTiming);

		refRateResetTiming = ARM_ConvPayResetRule (C_refRateResetTiming);

		stubRule = ARM_ConvStubRule (C_stubRule);

		if ( C_levelDownId == "DEFAULT" )
		{
			levelDownId = ARM_NULL_OBJECT;
		}
		else
		{
			levelDownId = LocalGetNumObjectId(C_levelDownId);
		}

		if ( ( downSpec = ARM_ConvStdBoost(C_downSpec, C_result) ) == ARM_DEFAULT_ERR )
		{
			ARM_ARG_ERR();
			return((LPXLOPER) &XL_result);
		}

		if ( C_levelUpId == "DEFAULT" )
		{
			levelUpId = ARM_NULL_OBJECT;
		}
		else
		{
			levelUpId = LocalGetNumObjectId(C_levelUpId);
		}

		if (( upSpec = ARM_ConvStdBoost(C_upSpec, C_result) ) == ARM_DEFAULT_ERR )
		{
			ARM_ARG_ERR();
			return((LPXLOPER) &XL_result);
		}

		if (( C_ccy.GetLen() > 3 )
			&&
			( !(C_ccy == "DEFAULT") )
		   )
		   ccyIsObject = true;

		if ( C_MCFreq == "-1" )
		{
			MCFreq = K_DAILY;
		}
		else
		{
			if ((MCFreq = ARM_ConvFrequency(C_MCFreq, C_result)) == ARM_DEFAULT_ERR)
			{
	  			ARM_ARG_ERR();
				return (LPXLOPER)&XL_result;
		   }
		}

		MCInterp = K_LINEAR;
		if ( C_MCInterp == "LOG" )
		{
			MCInterp = K_CONTINUOUS;
		}
		else if ( C_MCInterp == "RAC" )
		{
			MCInterp = K_SQRT;
		}
		else if ( C_MCInterp == "SL_LIN" )
		{
			MCInterp = K_SLOPELIN;
		}
		else if ( C_MCInterp == "SL_RAC" )
		{
			MCInterp = K_SLOPESQRT;
		}

		LDPricingMethodId = ARM_ConvPricCorridorLD(C_LDPricingMethod) ;

		long retCode;
		long objId;
		
		CCString curClass = LOCAL_CORRIDORLEG_CLASS;
		CCString stringId;

		retCode = ARMLOCAL_CORRIDORLEG(C_startDate,
									   C_endDate,
									   receiveOrPay,
									   payIndexId,
									   payFreq,
									   spreadType,
									   C_spread_double,
									   refIndexId,
									   resetFreq,
									   paidRateResetTiming,
									   refRateResetTiming,
									   stubRule,
									   levelDownId,
									   downSpec,
									   levelUpId,
									   upSpec,
									   ccyIsObject,
									   C_ccy,
									   MCFreq,
									   MCInterp,
									   LDPricingMethodId,
                                       long(decompPricingFlag),
									   "NULL", // reset calendar
									   "NULL", // pay calendar
                                       long(C_RefIndexResetGap),
									   C_result);

		if ( retCode == ARM_OK )
		{
			objId = C_result.getLong ();

			stringId = LocalMakeObjectId (objId, curClass);
		}

		if ( retCode == ARM_OK )
		{
			FreeCurCellErr();

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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_PXL_CORRIDORLEG" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}




LPXLOPER Local_CorridorLegWithFixings_Common(  LPXLOPER		XL_startDate,
											   LPXLOPER		XL_endDate,
											   LPXLOPER		XL_payReceive,
											   LPXLOPER		XL_payIndex,
											   LPXLOPER		XL_payFreq,
											   LPXLOPER		XL_spread,
											   LPXLOPER		XL_refIndex,
											   LPXLOPER		XL_resetFreq,
											   LPXLOPER		XL_paidRateResetTiming,
											   LPXLOPER		XL_refRateResetTiming,
											   LPXLOPER		XL_stubRule,
											   LPXLOPER		XL_levelDown,
											   LPXLOPER		XL_downSpec,
											   LPXLOPER		XL_levelUp,
											   LPXLOPER		XL_upSpec,
											   LPXLOPER		XL_ccy,
											   LPXLOPER		XL_decompPricingFlag,
											   LPXLOPER		XL_refIndexResetGap,
											   LPXLOPER		XL_refPastFixings,
											   LPXLOPER		XL_payPastFixings,
											   bool			PersistentInXL)
{
	//This is defined first because it is used in XL macros
	static XLOPER XL_result;
	
	//Get the variables from the XLOper variables
	ARM_result C_result;
	
	//To make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
		//To avoid computation if called by the wizard
		ARM_NOCALCIFWIZ();
		
		//This is used by macros 
		//and therefore this has to be defined
		static int error;
		static char* reason = "";

		//StartDate
		double startDate;
		XL_readNumCell(XL_startDate, startDate, "ARM_ERR: Start Date: date expected",C_result);
		
		//EndDate
		double endDate;
		XL_readNumCell(XL_endDate, endDate, "ARM_ERR: End Date: date expected",C_result);

		//ReceiveOrPay
		CCString payReceiveStr;
		long lPayReceive;

		XL_readStrCell(XL_payReceive, payReceiveStr, "ARM_ERR: Payer / Receiver: string expected",C_result);
	
        if((lPayReceive = ARM_ConvRecOrPay (payReceiveStr, C_result)) == ARM_DEFAULT_ERR)
		{
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}
		int payReceive = (int)lPayReceive;

		//PaymentIndex
		CCString payIndexStr;
		long payIndexId;
		XL_readStrCell(XL_payIndex, payIndexStr, "ARM_ERR: Payment Index: String Id expected", C_result);
		payIndexId = LocalGetNumObjectId(payIndexStr);

		//PaymentFrequency
		CCString C_payFreq;
		long payFreq;
		XL_readStrCellWD(XL_payFreq, C_payFreq, "-1", "ARM_ERR: Payment Frequency: string expected",	C_result);
		if ( C_payFreq == "-1" )
		{
		   payFreq = K_DEF_FREQ;
		}
		else
		{
			if ((payFreq = ARM_ConvFrequency(C_payFreq, C_result)) == ARM_DEFAULT_ERR)
			{
				ARM_ARG_ERR();
				return (LPXLOPER)&XL_result;
			}
		}

		//Spread
		double spread;
		CCString spreadStr;
		long     spreadId;
		XL_readStrOrNumCell(XL_spread, spreadStr, spread, spreadId, "ARM_ERR: Spread: numeric or refValue Id expected",C_result);	
		if(spreadId == XL_TYPE_STRING)
		{
			spreadId = LocalGetNumObjectId(spreadStr);
		}
		else
		{
			spreadId = ARM_NULL_OBJECT;
		}

		//ReferenceIndex
		CCString refIndexStr;
		long refIndexId;
		XL_readStrCell(XL_refIndex, refIndexStr,"ARM_ERR: Reference Index: String Id expected", C_result);
		refIndexId = LocalGetNumObjectId(refIndexStr);

		//ResetFrequency
		CCString C_resetFreq;
		long resetFreq;
		XL_readStrCellWD(XL_resetFreq, C_resetFreq, "-1","ARM_ERR: Reset Frequency: string expected",	C_result);
		if ( C_resetFreq == "-1" )
		{
		   resetFreq = K_DEF_FREQ;
		}
		else
		{
			if ((resetFreq = ARM_ConvFrequency(C_resetFreq, C_result)) == ARM_DEFAULT_ERR)
			{
				ARM_ARG_ERR();
				return (LPXLOPER)&XL_result;
			}
		}

		//PaidRateResetTiming
		CCString C_paidRateResetTiming;
		long paidRateResetTiming;
		XL_readStrCellWD(XL_paidRateResetTiming, C_paidRateResetTiming, "ADV","ARM_ERR: Paid Rate Reset Timing: string expected",C_result);
		paidRateResetTiming = ARM_ConvPayResetRule(C_paidRateResetTiming);

		//RefRateResetTiming
		CCString C_refRateResetTiming;
		long refRateResetTiming;
		XL_readStrCellWD(XL_refRateResetTiming, C_refRateResetTiming, "ADV","ARM_ERR: Ref Rate Reset Timing: string expected",C_result);
		refRateResetTiming = ARM_ConvPayResetRule(C_refRateResetTiming);

		//StubRule
		CCString C_stubRule;
		long stubRule;
		XL_readStrCellWD(XL_stubRule, C_stubRule, "SS"," ARM_ERR: Stub Rule: string expected",C_result);
		stubRule = ARM_ConvStubRule (C_stubRule);

		//LevelDown
		double levelDown;
		CCString levelDownStr;
		long     levelDownId;
		XL_readStrOrNumCell(XL_levelDown, levelDownStr, levelDown, levelDownId, "ARM_ERR: Level Down: numeric or refValue Id expected",C_result);	
		if(levelDownId == XL_TYPE_STRING)
		{
			levelDownId = LocalGetNumObjectId(levelDownStr);
		}
		else
		{
			levelDownId = ARM_NULL_OBJECT;
		}

		//DownSpec
		CCString C_downSpec;
		long     downSpec;
		XL_readStrCellWD(XL_downSpec, C_downSpec, "STD","ARM_ERR: Down Spec: String expected",C_result);
		if ( ( downSpec = ARM_ConvStdBoost(C_downSpec, C_result) ) == ARM_DEFAULT_ERR )
		{
			ARM_ARG_ERR();
			return((LPXLOPER) &XL_result);
		}

		//LevelUp
		double levelUp;
		CCString levelUpStr;
		long     levelUpId;
		XL_readStrOrNumCell(XL_levelUp, levelUpStr, levelUp, levelUpId, "ARM_ERR: Level Up: numeric or refValue Id expected",C_result);	
		if(levelUpId == XL_TYPE_STRING)
		{
			levelUpId = LocalGetNumObjectId(levelUpStr);
		}
		else
		{
			levelUpId = ARM_NULL_OBJECT;
		}

		//UpSpec
		CCString C_upSpec;
		long     upSpec;
		XL_readStrCellWD(XL_upSpec, C_upSpec, "STD","ARM_ERR: Up Spec: String expected",C_result);
		if ( ( upSpec = ARM_ConvStdBoost(C_upSpec, C_result) ) == ARM_DEFAULT_ERR )
		{
			ARM_ARG_ERR();
			return((LPXLOPER) &XL_result);
		}

		//Currency
		CCString C_ccy;
		long ccyId;
		XL_readStrCellWD(XL_ccy, C_ccy , "DEFAULT"," ARM_ERR: Currency id: string expected",C_result);
		if (C_ccy == "DEFAULT")
		{
			ccyId =  ARM_NULL_OBJECT;
		}
		else
		{
		    ccyId = LocalGetNumObjectId(C_ccy);
		}
		
		//DecompPricingFlag
		double decompPricingFlag;
		double decompPricingFlag_default = 1.0;
		XL_readNumCellWD(XL_decompPricingFlag, decompPricingFlag, decompPricingFlag_default,"ARM_ERR: DecompPricingFlag: string expected",C_result);
		
		//RefIndexResetGap
		double refIndexResetGap;
		double refIndexResetGap_default = 0.0;
		XL_readNumCellWD(XL_refIndexResetGap, refIndexResetGap, refIndexResetGap_default,"ARM_ERR: Ref Index Reset Gap: int expected",C_result);

		//RefPastFixings
		CCString C_refPastFixings;
		long refPastFixingsId;
		XL_readStrCellWD(XL_refPastFixings, C_refPastFixings,"DEFAULT","ARM_ERR: Ref Index Past Fixings: object expected",C_result);
	    if (C_refPastFixings == "DEFAULT")
		{
			refPastFixingsId =  ARM_NULL_OBJECT;
		}
		else
		{
		    refPastFixingsId = LocalGetNumObjectId(C_refPastFixings);
		}

		//PayPastFixings
		CCString C_payPastFixings;
		long payPastFixingsId;
		XL_readStrCellWD(XL_payPastFixings, C_payPastFixings,"DEFAULT","ARM_ERR: Pay Index Past Fixings: object expected",C_result);
	    if (C_payPastFixings == "DEFAULT")
		{
			payPastFixingsId = ARM_NULL_OBJECT;
		}
		else
		{
		    payPastFixingsId = LocalGetNumObjectId(C_payPastFixings);
		}

		//Functor call
		exportFunc23Args< double, 
						  double, 
						  int, 
						  long, 
						  long, 
						  double, 
						  long, 
						  long, 
						  long, 
						  long, 
						  long, 
						  long, 
						  double, 
						  long, 
						  long, 
						  double, 
						  long, 
						  long, 
						  long, 
						  double, 
						  double, 
						  long, 
						  long > ourFunc( 
						  startDate, 
						  endDate, 
						  payReceive, 
						  payIndexId, 
						  payFreq, 
						  spread, 
						  spreadId, 
						  refIndexId, 
						  resetFreq, 
						  paidRateResetTiming, 
						  refRateResetTiming,
						  stubRule, 
						  levelDown, 
						  levelDownId, 
						  downSpec, 
						  levelUp, 
						  levelUpId, 
						  upSpec, 
						  ccyId, 
						  decompPricingFlag, 
						  refIndexResetGap, 
						  refPastFixingsId, 
						  payPastFixingsId,
						  ARMLOCAL_CORRIDORLEGWITHFIXINGS);

		//Call the general function
		fillXL_Result( LOCAL_CORRIDORLEG_CLASS, ourFunc, C_result, XL_result, PersistentInXL );
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_CorridorLegWithFixings_Common" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


__declspec(dllexport) LPXLOPER Local_PXL_CorridorLegWithFixings( LPXLOPER XL_startDate,
																 LPXLOPER XL_endDate,
																 LPXLOPER XL_receiveOrPay,
																 LPXLOPER XL_payIndex,
																 LPXLOPER XL_payFreq,
																 LPXLOPER XL_spread,
																 LPXLOPER XL_refIndex,
																 LPXLOPER XL_resetFreq,
																 LPXLOPER XL_paidRateResetTiming,
																 LPXLOPER XL_refRateResetTiming,
																 LPXLOPER XL_stubRule,
																 LPXLOPER XL_levelDown,
																 LPXLOPER XL_downSpec,
																 LPXLOPER XL_levelUp,
																 LPXLOPER XL_upSpec,
																 LPXLOPER XL_ccy,
																 LPXLOPER XL_decompPricingFlag,
																 LPXLOPER XL_refIndexResetGap,
																 LPXLOPER XL_refPastFixings,
																 LPXLOPER XL_payPastFixings)
{
	bool PersistentInXL = false;

	return Local_CorridorLegWithFixings_Common( XL_startDate,
												XL_endDate,
												XL_receiveOrPay,
												XL_payIndex,
												XL_payFreq,
												XL_spread,
												XL_refIndex,
												XL_resetFreq,
												XL_paidRateResetTiming,
												XL_refRateResetTiming,
												XL_stubRule,
												XL_levelDown,
												XL_downSpec,
												XL_levelUp,
												XL_upSpec,
												XL_ccy,
												XL_decompPricingFlag,
												XL_refIndexResetGap,
												XL_refPastFixings,
												XL_payPastFixings,
												PersistentInXL);
}



__declspec(dllexport) LPXLOPER WINAPI Local_CorridorLegWithFixings( LPXLOPER XL_startDate,
																	LPXLOPER XL_endDate,
																	LPXLOPER XL_receiveOrPay,
																	LPXLOPER XL_payIndex,
																	LPXLOPER XL_payFreq,
																	LPXLOPER XL_spread,
																	LPXLOPER XL_refIndex,
																	LPXLOPER XL_resetFreq,
																	LPXLOPER XL_paidRateResetTiming,
																	LPXLOPER XL_refRateResetTiming,
																	LPXLOPER XL_stubRule,
																	LPXLOPER XL_levelDown,
																	LPXLOPER XL_downSpec,
																	LPXLOPER XL_levelUp,
																	LPXLOPER XL_upSpec,
																	LPXLOPER XL_ccy,
																	LPXLOPER XL_decompPricingFlag,
																	LPXLOPER XL_refIndexResetGap,
																	LPXLOPER XL_refPastFixings,
																	LPXLOPER XL_payPastFixings)
{
	ADD_LOG("Local_CorridorLegWithFixings");
	bool PersistentInXL = true;

	return Local_CorridorLegWithFixings_Common( XL_startDate,
												XL_endDate,
												XL_receiveOrPay,
												XL_payIndex,
												XL_payFreq,
												XL_spread,
												XL_refIndex,
												XL_resetFreq,
												XL_paidRateResetTiming,
												XL_refRateResetTiming,
												XL_stubRule,
												XL_levelDown,
												XL_downSpec,
												XL_levelUp,
												XL_upSpec,
												XL_ccy,
												XL_decompPricingFlag,
												XL_refIndexResetGap,
												XL_refPastFixings,
												XL_payPastFixings,
												PersistentInXL);
}



__declspec(dllexport) LPXLOPER WINAPI Local_CORRIDORLEG_WITHCALENDARS(LPXLOPER XL_startDate,
																	  LPXLOPER XL_endDate,
																	  LPXLOPER XL_receiveOrPay,
																	  LPXLOPER XL_payIndexId,
																	  LPXLOPER XL_payFreq,
																	  LPXLOPER XL_spread,
																	  LPXLOPER XL_refIndexId,
																	  LPXLOPER XL_resetFreq,
																	  LPXLOPER XL_paidRateResetTiming,
																	  LPXLOPER XL_refRateResetTiming,
																	  LPXLOPER XL_stubRule,
																	  LPXLOPER XL_levelDownId,
																	  LPXLOPER XL_downSpec,
																	  LPXLOPER XL_levelUpId,
																	  LPXLOPER XL_upSpec,
																	  LPXLOPER XL_ccy,
																	  LPXLOPER XL_resetCal,
																	  LPXLOPER XL_payCal,
																	  LPXLOPER XL_decompPricingFlag,
                                                                      LPXLOPER XL_refIndexResetGap)
{
	ADD_LOG("Local_CORRIDORLEG_WITHCALENDARS");
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
		long receiveOrPay;

		CCString C_payIndexId;
		long payIndexId;

		CCString C_payFreq;
		long payFreq;

		CCString C_spread_str;
		double C_spread_double;
		double C_spread_default = 0.0;
		long spreadType;

		CCString C_refIndexId;
		long refIndexId;

		CCString C_resetFreq;
		long resetFreq;

		CCString C_paidRateResetTiming;
		long paidRateResetTiming;

		CCString C_refRateResetTiming;
		long refRateResetTiming;
		
		CCString C_stubRule;
		long stubRule;

		CCString C_levelDownId;
		long levelDownId;
    						  
		CCString C_downSpec;
		long     downSpec;

		CCString C_levelUpId;
		long levelUpId;

		CCString C_upSpec;
		long     upSpec;

		CCString C_ccy;
		bool ccyIsObject = false;

		CCString C_resetCal;
		CCString C_payCal;

		double decompPricingFlag;
		double decompPricingFlag_default = 1.0;

        double C_RefIndexResetGap;
        double C_RefIndexresetGap_default = 0.0;
        
        // error
		static int error;
		static char* reason = "";

		XL_readNumCell(XL_startDate,C_startDate," ARM_ERR: start date: date expected",C_result);
		XL_readNumCell(XL_endDate,C_endDate," ARM_ERR: end date: date expected",C_result);
		XL_readStrCell(XL_receiveOrPay,C_receiveOrPay," ARM_ERR: receive or pay: string expected",C_result);    
		XL_readStrCell(XL_payIndexId,C_payIndexId," ARM_ERR: Object paymentIndex Id expected",C_result);
		XL_readStrCellWD(XL_payFreq,C_payFreq,"-1"," ARM_ERR: frequency: string expected",C_result);
		XL_readStrOrNumCellWD(XL_spread,C_spread_str,C_spread_double,C_spread_default,spreadType," ARM_ERR: spread: numeric or object expected",C_result);
		XL_readStrCell(XL_refIndexId,C_refIndexId," ARM_ERR: Object refindex Id expected",C_result);
		XL_readStrCellWD(XL_resetFreq,C_resetFreq,"-1"," ARM_ERR: frequency: string expected",C_result);
		XL_readStrCellWD(XL_paidRateResetTiming,C_paidRateResetTiming,"ADV"," ARM_ERR: paidRateResetTiming: string expected",C_result);
		XL_readStrCellWD(XL_refRateResetTiming,C_refRateResetTiming,"ADV"," ARM_ERR: refRateResetTiming: string expected",C_result);
		XL_readStrCellWD(XL_stubRule,C_stubRule,"SS"," ARM_ERR: stub rule: string expected",C_result);
		XL_readStrCellWD(XL_levelDownId,C_levelDownId,"DEFAULT"," ARM_ERR: Object refValue expected",C_result);
		XL_readStrCellWD(XL_downSpec,C_downSpec,"STD"," ARM_ERR: DownSpec: String expected",C_result);
		XL_readStrCellWD(XL_levelUpId,C_levelUpId,"DEFAULT"," ARM_ERR: Object RefValue expected",C_result);
		XL_readStrCellWD(XL_upSpec,C_upSpec,"STD"," ARM_ERR: DownSpec: String expected",C_result);
		XL_readStrCellWD(XL_ccy,C_ccy,"DEFAULT"," ARM_ERR: currency id: string expected",C_result);
		XL_readStrCellWD(XL_resetCal,C_resetCal,"NULL"," ARM_ERR: reset Calendar: string expected",C_result);
		XL_readStrCellWD(XL_payCal,C_payCal,"NULL"," ARM_ERR: pay Calendar: string expected",C_result);
		XL_readNumCellWD(XL_decompPricingFlag,decompPricingFlag,decompPricingFlag_default," ARM_ERR: decompPricingFlag: string expected",C_result);
		
        XL_readNumCellWD(XL_refIndexResetGap, C_RefIndexResetGap, 
                         C_RefIndexresetGap_default," ARM_ERR: Ref Index reset gap: numeric expected", C_result);

		if ((receiveOrPay = ARM_ConvRecOrPay(C_receiveOrPay, C_result)) == ARM_DEFAULT_ERR)
		{
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}

		payIndexId = LocalGetNumObjectId(C_payIndexId);
		refIndexId = LocalGetNumObjectId(C_refIndexId);

		if ( C_payFreq == "-1" )
		{
		   payFreq = K_DEF_FREQ;
		}
		else
		{
			if ((payFreq = ARM_ConvFrequency(C_payFreq, C_result)) == ARM_DEFAULT_ERR)
			{
	  			ARM_ARG_ERR();
				return (LPXLOPER)&XL_result;
			}
		}

		if ( spreadType == XL_TYPE_STRING )
		{
			C_spread_double = (double) LocalGetNumObjectId(C_spread_str);

			spreadType = 1L;
		}
		else
		{
			spreadType = 0L;
		}
			
		if ( C_resetFreq == "-1" )
		{
			resetFreq = K_DAILY;
		}
		else
		{
			if ((resetFreq = ARM_ConvFrequency(C_resetFreq, C_result)) == ARM_DEFAULT_ERR)
			{
	  			ARM_ARG_ERR();
				return (LPXLOPER)&XL_result;
			}
		}

		paidRateResetTiming = ARM_ConvPayResetRule (C_paidRateResetTiming);

		refRateResetTiming = ARM_ConvPayResetRule (C_refRateResetTiming);

		stubRule = ARM_ConvStubRule (C_stubRule);

		if ( C_levelDownId == "DEFAULT" )
		{
			levelDownId = ARM_NULL_OBJECT;
		}
		else
		{
			levelDownId = LocalGetNumObjectId(C_levelDownId);
		}

		if ( ( downSpec = ARM_ConvStdBoost(C_downSpec, C_result) ) == ARM_DEFAULT_ERR )
		{
			ARM_ARG_ERR();
			return((LPXLOPER) &XL_result);
		}

		if ( C_levelUpId == "DEFAULT" )
		{
			levelUpId = ARM_NULL_OBJECT;
		}
		else
		{
			levelUpId = LocalGetNumObjectId(C_levelUpId);
		}

		if (( upSpec = ARM_ConvStdBoost(C_upSpec, C_result) ) == ARM_DEFAULT_ERR )
		{
			ARM_ARG_ERR();
			return((LPXLOPER) &XL_result);
		}

		if (( C_ccy.GetLen() > 3 ) 
			&& 
			( !(C_ccy == "DEFAULT") )
		   )
		   ccyIsObject = true;

		long retCode;
		long objId;
		CCString prevClass;
		
		CCString curClass = LOCAL_CORRIDORLEG_CLASS;
		CCString stringId = GetLastCurCellEnvValue ();

		if (!stringId)
		{
			retCode = ARMLOCAL_CORRIDORLEG(C_startDate,
										   C_endDate,
										   receiveOrPay,
										   payIndexId,
										   payFreq,
										   spreadType,
										   C_spread_double,
										   refIndexId,
										   resetFreq,
										   paidRateResetTiming,
										   refRateResetTiming,
										   stubRule,
										   levelDownId,
										   downSpec,
										   levelUpId,
										   upSpec,
										   ccyIsObject,
										   C_ccy,
										   K_DEF_FREQ,	// MC Freq
										   K_LINEAR,	// MC Interpol
										   K_DIGITALE,	// LD Pricing Method
										   long(decompPricingFlag),
   										   C_resetCal,
										   C_payCal,
                                           long(C_RefIndexResetGap),
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
				
			if ( curClass == prevClass )
			{
				retCode = ARMLOCAL_CORRIDORLEG(C_startDate,
											   C_endDate,
											   receiveOrPay,
											   payIndexId,
											   payFreq,
											   spreadType,
											   C_spread_double,
											   refIndexId,
											   resetFreq,
											   paidRateResetTiming,
											   refRateResetTiming,
											   stubRule,
											   levelDownId,
											   downSpec,
											   levelUpId,
											   upSpec,
											   ccyIsObject,
											   C_ccy,
											   K_DEF_FREQ,	// MC Freq
											   K_LINEAR,	// MC Interpol
											   K_DIGITALE,	// LD Pricing Method
											   long(decompPricingFlag),
											   C_resetCal,
											   C_payCal,
                                               long(C_RefIndexResetGap),
											   C_result,
											   objId);

				if ( retCode == ARM_OK )
				{
					LocalSetCurCellEnvValue (curClass, objId); 

					stringId = LocalMakeObjectId (objId, curClass);
				}
			}
			else
			{
				FreeCurCellContent();

				retCode = ARMLOCAL_CORRIDORLEG(C_startDate,
											   C_endDate,
											   receiveOrPay,
											   payIndexId,
											   payFreq,
											   spreadType,
											   C_spread_double,
											   refIndexId,
											   resetFreq,
											   paidRateResetTiming,
											   refRateResetTiming,
											   stubRule,
											   levelDownId,
											   downSpec,
											   levelUpId,
											   upSpec,
											   ccyIsObject,
											   C_ccy,
											   K_DEF_FREQ,	// MC Freq
											   K_LINEAR,	// MC Interpol
											   K_DIGITALE,	// LD Pricing Method
											   long(decompPricingFlag),
											   C_resetCal,
											   C_payCal,
                                               long(C_RefIndexResetGap),
											   C_result);

				if ( retCode == ARM_OK )
				{
					objId = C_result.getLong ();
				
					LocalSetCurCellEnvValue (curClass, objId); 

					stringId = LocalMakeObjectId (objId, curClass);
				}
			}
		}

		if ( retCode == ARM_OK )
		{
			FreeCurCellErr();

			XL_result.xltype  = xltypeStr;
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_CORRIDORLEG_WITHCALENDARS" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}



__declspec(dllexport) LPXLOPER WINAPI Local_PXL_CORRIDORLEG_WITHCALENDARS(LPXLOPER XL_startDate,
															LPXLOPER XL_endDate,
															LPXLOPER XL_receiveOrPay,
															LPXLOPER XL_payIndexId,
															LPXLOPER XL_payFreq,
															LPXLOPER XL_spread,
															LPXLOPER XL_refIndexId,
															LPXLOPER XL_resetFreq,
															LPXLOPER XL_paidRateResetTiming,
															LPXLOPER XL_refRateResetTiming,
															LPXLOPER XL_stubRule,
															LPXLOPER XL_levelDownId,
															LPXLOPER XL_downSpec,
															LPXLOPER XL_levelUpId,
															LPXLOPER XL_upSpec,
															LPXLOPER XL_ccy,
															LPXLOPER XL_resetCal,
															LPXLOPER XL_payCal,
                                                            LPXLOPER XL_decompPricingFlag,
                                                            LPXLOPER XL_refIndexResetGap)
{
	ADD_LOG("Local_PXL_CORRIDORLEG_WITHCALENDARS");
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
		long receiveOrPay;

		CCString C_payIndexId;
		long payIndexId;

		CCString C_payFreq;
		long payFreq;

		CCString C_spread_str;
		double C_spread_double;
		double C_spread_default = 0.0;
		long spreadType;

		CCString C_refIndexId;
		long refIndexId;

		CCString C_resetFreq;
		long resetFreq;

		CCString C_paidRateResetTiming;
		long paidRateResetTiming;

		CCString C_refRateResetTiming;
		long refRateResetTiming;
		
		CCString C_stubRule;
		long stubRule;

		CCString C_levelDownId;
		long levelDownId;
    						  
		CCString C_downSpec;
		long     downSpec;

		CCString C_levelUpId;
		long levelUpId;

		CCString C_upSpec;
		long     upSpec;

		CCString C_ccy;
		bool ccyIsObject = false;

		CCString C_resetCal;
		CCString C_payCal;

		double decompPricingFlag;
		double decompPricingFlag_default = 1.0;

        double C_RefIndexResetGap;
        double C_RefIndexresetGap_default = 0.0;

		// error
		static int error;
		static char* reason = "";

		XL_readNumCell(XL_startDate,C_startDate," ARM_ERR: start date: date expected",C_result);
		XL_readNumCell(XL_endDate,C_endDate," ARM_ERR: end date: date expected",C_result);
		XL_readStrCell(XL_receiveOrPay,C_receiveOrPay," ARM_ERR: receive or pay: string expected",C_result);    
		XL_readStrCell(XL_payIndexId,C_payIndexId, "ARM_ERR: Object PaymentIndex Id expected",C_result);
		XL_readStrCellWD(XL_payFreq,C_payFreq,"-1"," ARM_ERR: frequency: string expected",C_result);
		XL_readStrOrNumCellWD(XL_spread,C_spread_str,C_spread_double,C_spread_default,spreadType," ARM_ERR: spread: numeric or object expected",C_result);
		XL_readStrCell(XL_refIndexId,C_refIndexId," ARM_ERR: Object refIdex Id expected",C_result);
		XL_readStrCellWD(XL_resetFreq,C_resetFreq,"-1"," ARM_ERR: frequency: string expected",C_result);
		XL_readStrCellWD(XL_paidRateResetTiming,C_paidRateResetTiming,"ADV"," ARM_ERR: paidRateResetTiming: string expected",C_result);
		XL_readStrCellWD(XL_refRateResetTiming,C_refRateResetTiming,"ADV"," ARM_ERR: refRateResetTiming: string expected",C_result);
		XL_readStrCellWD(XL_stubRule,C_stubRule,"SS"," ARM_ERR: stub rule: string expected",C_result);
		XL_readStrCellWD(XL_levelDownId,C_levelDownId,"DEFAULT"," ARM_ERR: Object RefValue expected",C_result);
		XL_readStrCellWD(XL_downSpec,C_downSpec,"STD"," ARM_ERR: DownSpec: String expected",C_result);
		XL_readStrCellWD(XL_levelUpId,C_levelUpId,"DEFAULT"," ARM_ERR: Object RefValue expected",C_result);
		XL_readStrCellWD(XL_upSpec,C_upSpec,"STD"," ARM_ERR: DownSpec: String expected",C_result);
		XL_readStrCellWD(XL_ccy,C_ccy,"DEFAULT"," ARM_ERR: currency id: string expected",C_result);
		XL_readStrCellWD(XL_resetCal,C_resetCal,"NULL"," ARM_ERR: reset Calendar: string expected",C_result);
		XL_readStrCellWD(XL_payCal,C_payCal,"NULL"," ARM_ERR: pay Calendar: string expected",C_result);
		XL_readNumCellWD(XL_decompPricingFlag,decompPricingFlag,decompPricingFlag_default," ARM_ERR: decompPricingFlag: string expected",C_result);

        XL_readNumCellWD(XL_refIndexResetGap, C_RefIndexResetGap, 
                         C_RefIndexresetGap_default," ARM_ERR: Ref Index reset gap: numeric expected", C_result);

		
		if ((receiveOrPay = ARM_ConvRecOrPay(C_receiveOrPay, C_result)) == ARM_DEFAULT_ERR)
		{
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}

		payIndexId = LocalGetNumObjectId(C_payIndexId);
		refIndexId = LocalGetNumObjectId(C_refIndexId);

		if ( C_payFreq == "-1" )
		{
		   payFreq = K_DEF_FREQ;
		}
		else
		{
			if ((payFreq = ARM_ConvFrequency(C_payFreq, C_result)) == ARM_DEFAULT_ERR)
			{
	  			ARM_ARG_ERR();
				return (LPXLOPER)&XL_result;
			}
		}

		if ( spreadType == XL_TYPE_STRING )
		{
			C_spread_double = (double) LocalGetNumObjectId(C_spread_str);

			spreadType = 1L;
		}
		else
		{
			spreadType = 0L;
		}
    

		if ( C_resetFreq == "-1" )
		{
			resetFreq = K_DAILY;
		}
		else
		{
			if ((resetFreq = ARM_ConvFrequency(C_resetFreq, C_result)) == ARM_DEFAULT_ERR)
			{
	  			ARM_ARG_ERR();
				return (LPXLOPER)&XL_result;
			}
		}

		paidRateResetTiming = ARM_ConvPayResetRule (C_paidRateResetTiming);

		refRateResetTiming = ARM_ConvPayResetRule (C_refRateResetTiming);

		stubRule = ARM_ConvStubRule (C_stubRule);

		if ( C_levelDownId == "DEFAULT" )
		{
			levelDownId = ARM_NULL_OBJECT;
		}
		else
		{
			levelDownId = LocalGetNumObjectId(C_levelDownId);
		}

		if ( ( downSpec = ARM_ConvStdBoost(C_downSpec, C_result) ) == ARM_DEFAULT_ERR )
		{
			ARM_ARG_ERR();
			return((LPXLOPER) &XL_result);
		}

		if ( C_levelUpId == "DEFAULT" )
		{
			levelUpId = ARM_NULL_OBJECT;
		}
		else
		{
			levelUpId = LocalGetNumObjectId(C_levelUpId);
		}

		if (( upSpec = ARM_ConvStdBoost(C_upSpec, C_result) ) == ARM_DEFAULT_ERR )
		{
			ARM_ARG_ERR();
			return((LPXLOPER) &XL_result);
		}

		if (( C_ccy.GetLen() > 3 )
			&&
			( !(C_ccy == "DEFAULT") )
		   )
		   ccyIsObject = true;

		long retCode;
		long objId;
		
		CCString curClass = LOCAL_CORRIDORLEG_CLASS;
		CCString stringId;

		retCode = ARMLOCAL_CORRIDORLEG(C_startDate,
									   C_endDate,
									   receiveOrPay,
									   payIndexId,
									   payFreq,
									   spreadType,
									   C_spread_double,
									   refIndexId,
									   resetFreq,
									   paidRateResetTiming,
									   refRateResetTiming,
									   stubRule,
									   levelDownId,
									   downSpec,
									   levelUpId,
									   upSpec,
									   ccyIsObject,
									   C_ccy,
									   K_DEF_FREQ,	// MC Freq
									   K_LINEAR,	// MC Interpol
									   K_DIGITALE,	// LD Pricing Method
									   long(decompPricingFlag),
									   C_resetCal,
									   C_payCal,
                                       long(C_RefIndexResetGap),
									   C_result);

		if ( retCode == ARM_OK )
		{
			objId = C_result.getLong ();

			stringId = LocalMakeObjectId (objId, curClass);
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_PXL_CORRIDORLEG_WITHCALENDARS" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}



__declspec(dllexport) LPXLOPER WINAPI Local_RESTRIKABLELEG(LPXLOPER XL_startDate,
														   LPXLOPER XL_endDate,
														   LPXLOPER XL_receiveOrPay,
														   LPXLOPER XL_payIndexId,
														   LPXLOPER XL_spread,
														   LPXLOPER XL_refIndexId,
														   LPXLOPER XL_stubRule,
														   LPXLOPER XL_range,
														   LPXLOPER XL_ccy,
														   LPXLOPER XL_Alpha,
														   LPXLOPER XL_Beta,
                                                           LPXLOPER XL_decompPricingFlag,
														   LPXLOPER XL_AdjStartDateFlag)
{
	ADD_LOG("Local_RESTRIKABLELEG");
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
		long receiveOrPay;

		CCString C_payIndexId;
		long payIndexId;

		CCString C_spreadStr;
		double C_spreadDouble;
		double C_spread_default = 0.0;
		long spreadType;

		CCString C_refIndexId;
		long refIndexId;
	
		CCString C_stubRule;
		long stubRule;

		CCString C_rangeStr;
		double C_rangeDouble;
		double C_range_default = 0.0;
		long rangeType;
    						  
		CCString C_ccy;
		bool ccyIsObject = false;
				
		double C_Alpha;
		double C_Alpha_default = 0.5;

		double C_Beta;
		double C_Beta_default = 0.5;

		double C_decompPricingFlag;
		double C_decompPricingFlag_default(1.0);

		double C_AdjStartDateFlag;
		double C_AdjStartDateFlag_default(1.0);

		// error
		static int error;
		static char* reason = "";

		XL_readNumCell(XL_startDate,C_startDate," ARM_ERR: start date: date expected",C_result);
		XL_readNumCell(XL_endDate,C_endDate," ARM_ERR: end date: date expected",C_result);
		XL_readStrCell(XL_receiveOrPay,C_receiveOrPay," ARM_ERR: receive or pay: string expected",C_result);  
		XL_readStrCellWD(XL_payIndexId,C_payIndexId,"DEFAULT"," ARM_ERR: Object index ID expected",C_result);
		XL_readStrOrNumCellWD(XL_spread,C_spreadStr,C_spreadDouble,C_spread_default,spreadType," ARM_ERR: spread: numeric expected",C_result);
		XL_readStrCellWD(XL_refIndexId,C_refIndexId,"DEFAULT"," ARM_ERR: Object index ID expected",C_result);
		XL_readStrCellWD(XL_stubRule,C_stubRule,"SS"," ARM_ERR: stub rule: string expected",C_result);
		XL_readStrOrNumCellWD(XL_range,C_rangeStr,C_rangeDouble,C_range_default,rangeType," ARM_ERR: range: numeric expected",C_result);
		XL_readStrCellWD(XL_ccy,C_ccy,"DEFAULT"," ARM_ERR: currency id: string expected",C_result);
		XL_readNumCellWD(XL_Alpha,C_Alpha, C_Alpha_default," ARM_ERR: Alpha: double expected",C_result);
		XL_readNumCellWD(XL_Beta,C_Beta ,C_Beta_default," ARM_ERR: Beta: double expected",C_result);
		XL_readNumCellWD(XL_decompPricingFlag,C_decompPricingFlag,C_decompPricingFlag_default," ARM_ERR: decompPricingFlag: numeric expected",C_result);
		XL_readNumCellWD(XL_AdjStartDateFlag,C_AdjStartDateFlag,C_AdjStartDateFlag_default," ARM_ERR: AdjStartDateFlag: numeric expected",C_result);


		if ((receiveOrPay = ARM_ConvRecOrPay(C_receiveOrPay, C_result)) == ARM_DEFAULT_ERR)
		{
			ARM_ARG_ERR();		
			return (LPXLOPER)&XL_result;
		}

		if ( C_payIndexId == "DEFAULT" )
		{
			payIndexId = ARM_NULL_OBJECT;
		}
		else
		{
			payIndexId = LocalGetNumObjectId(C_payIndexId);
		}

		if ( C_refIndexId == "DEFAULT" )
		{
			refIndexId = ARM_NULL_OBJECT;
		}
		else
		{
			refIndexId = LocalGetNumObjectId(C_refIndexId);
		}

		stubRule = ARM_ConvStubRule (C_stubRule);


		if (( C_ccy.GetLen() > 3 )
			&&
			( !(C_ccy == "DEFAULT") )
		   )
		   ccyIsObject = true;

		if (spreadType == XL_TYPE_STRING)
			spreadType = LocalGetNumObjectId(C_spreadStr);
		else
			spreadType = ARM_NULL_OBJECT;

		if (rangeType == XL_TYPE_STRING)
			rangeType = LocalGetNumObjectId(C_rangeStr);
		else
			rangeType = ARM_NULL_OBJECT;

		long retCode;
		long objId;
		CCString prevClass;
		
		CCString curClass = LOCAL_CORRIDORLEG_CLASS;
		CCString stringId = GetLastCurCellEnvValue ();

		if (!stringId)
		{
			retCode = ARMLOCAL_RESTRIKABLELEG(C_startDate, C_endDate,
										 receiveOrPay, payIndexId,
										 spreadType, C_spreadDouble, 
										 refIndexId,stubRule,
										 rangeType, C_rangeDouble,
										 ccyIsObject, C_ccy,
										 C_Alpha, C_Beta,
										 (long) C_decompPricingFlag,
										 (long) C_AdjStartDateFlag,
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
				retCode = ARMLOCAL_RESTRIKABLELEG(C_startDate, C_endDate,
													receiveOrPay, payIndexId,
													spreadType, C_spreadDouble, 
													refIndexId, stubRule,
													rangeType, C_rangeDouble,
													ccyIsObject, C_ccy,
													C_Alpha, C_Beta,
													(long) C_decompPricingFlag,
													(long) C_AdjStartDateFlag,
													C_result, objId);

				if(retCode == ARM_OK)
				{
					LocalSetCurCellEnvValue (curClass, objId); 

					stringId = LocalMakeObjectId (objId, curClass);
				}
			}
			else
			{
				FreeCurCellContent();

				retCode = ARMLOCAL_RESTRIKABLELEG(C_startDate, C_endDate,
											 receiveOrPay, payIndexId,
											 spreadType, C_spreadDouble, 
											 refIndexId, stubRule,
											 rangeType, C_rangeDouble,
											 ccyIsObject, C_ccy,
											 C_Alpha, C_Beta,
											 (long) C_decompPricingFlag,
											 (long) C_AdjStartDateFlag,
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_RESTRIKABLELEG" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}



__declspec(dllexport) LPXLOPER WINAPI Local_PXL_RESTRIKABLELEG(LPXLOPER XL_startDate,
															   LPXLOPER XL_endDate,
															   LPXLOPER XL_receiveOrPay,
															   LPXLOPER XL_payIndexId,
															   LPXLOPER XL_spread,
															   LPXLOPER XL_refIndexId,
															   LPXLOPER XL_stubRule,
															   LPXLOPER XL_range,
															   LPXLOPER XL_ccy,
															   LPXLOPER XL_Alpha,
															   LPXLOPER XL_Beta,
                                                               LPXLOPER XL_decompPricingFlag,
															   LPXLOPER XL_AdjStartDateFlag)
{
	ADD_LOG("Local_PXL_RESTRIKABLELEG");
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
		long receiveOrPay;

		CCString C_payIndexId;
		long payIndexId;

		CCString C_spreadStr;
		double C_spreadDouble;
		double C_spread_default = 0.0;
		long spreadType;

		CCString C_refIndexId;
		long refIndexId;
	
		CCString C_stubRule;
		long stubRule;

		CCString C_rangeStr;
		double C_rangeDouble;
		double C_range_default = 0.0;
		long rangeType;
    						  
		CCString C_ccy;
		bool ccyIsObject = false;
				
		double C_Alpha;
		double C_Alpha_default = 0.5;

		double C_Beta;
		double C_Beta_default = 0.5;

		double C_decompPricingFlag;
		double C_decompPricingFlag_default(1.0);

		double C_AdjStartDateFlag;
		double C_AdjStartDateFlag_default(1.0);

		// error
		static int error;
		static char* reason = "";

		XL_readNumCell(XL_startDate,C_startDate," ARM_ERR: start date: date expected",C_result);
		XL_readNumCell(XL_endDate,C_endDate," ARM_ERR: end date: date expected",C_result);
		XL_readStrCell(XL_receiveOrPay,C_receiveOrPay," ARM_ERR: receive or pay: string expected",C_result);  
		XL_readStrCellWD(XL_payIndexId,C_payIndexId,"DEFAULT"," ARM_ERR: Object index ID expected",C_result);
		XL_readStrOrNumCellWD(XL_spread,C_spreadStr,C_spreadDouble,C_spread_default,spreadType," ARM_ERR: spread: numeric expected",C_result);
		XL_readStrCellWD(XL_refIndexId,C_refIndexId,"DEFAULT"," ARM_ERR: Object index ID expected",C_result);
		XL_readStrCellWD(XL_stubRule,C_stubRule,"SS"," ARM_ERR: stub rule: string expected",C_result);
		XL_readStrOrNumCellWD(XL_range,C_rangeStr,C_rangeDouble,C_range_default,rangeType," ARM_ERR: range: numeric expected",C_result);
		XL_readStrCellWD(XL_ccy,C_ccy,"DEFAULT"," ARM_ERR: currency id: string expected",C_result);
		XL_readNumCellWD(XL_Alpha,C_Alpha, C_Alpha_default," ARM_ERR: Alpha: double expected",C_result);
		XL_readNumCellWD(XL_Beta,C_Beta ,C_Beta_default," ARM_ERR: Beta: double expected",C_result);
		XL_readNumCellWD(XL_decompPricingFlag,C_decompPricingFlag,C_decompPricingFlag_default," ARM_ERR: decompPricingFlag: numeric expected",C_result);
		XL_readNumCellWD(XL_AdjStartDateFlag,C_AdjStartDateFlag,C_AdjStartDateFlag_default," ARM_ERR: AdjStartDateFlag: numeric expected",C_result);


		if ((receiveOrPay = ARM_ConvRecOrPay(C_receiveOrPay, C_result)) == ARM_DEFAULT_ERR)
		{
			ARM_ARG_ERR();		
			return (LPXLOPER)&XL_result;
		}

		if ( C_payIndexId == "DEFAULT" )
		{
			payIndexId = ARM_NULL_OBJECT;
		}
		else
		{
			payIndexId = LocalGetNumObjectId(C_payIndexId);
		}

		if ( C_refIndexId == "DEFAULT" )
		{
			refIndexId = ARM_NULL_OBJECT;
		}
		else
		{
			refIndexId = LocalGetNumObjectId(C_refIndexId);
		}

		stubRule = ARM_ConvStubRule (C_stubRule);


		if (( C_ccy.GetLen() > 3 )
			&&
			( !(C_ccy == "DEFAULT") )
		   )
		   ccyIsObject = true;

		if (spreadType == XL_TYPE_STRING)
			spreadType = LocalGetNumObjectId(C_spreadStr);
		else
			spreadType = ARM_NULL_OBJECT;

		if (rangeType == XL_TYPE_STRING)
			rangeType = LocalGetNumObjectId(C_rangeStr);
		else
			rangeType = ARM_NULL_OBJECT;

		long retCode;
		long objId;
		CCString prevClass;
		
		CCString curClass = LOCAL_CORRIDORLEG_CLASS;
		CCString stringId;

		retCode = ARMLOCAL_RESTRIKABLELEG(C_startDate, C_endDate,
									 receiveOrPay, payIndexId,
									 spreadType, C_spreadDouble, 
									 refIndexId,stubRule,
									 rangeType, C_rangeDouble,
									 ccyIsObject, C_ccy,
									 C_Alpha, C_Beta,
									 (long) C_decompPricingFlag,
									 (long) C_AdjStartDateFlag,
									 C_result);

		if(retCode == ARM_OK)
		{
			objId = C_result.getLong ();

			stringId = LocalMakeObjectId (objId, curClass);
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_PXL_RESTRIKABLELEG" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


__declspec(dllexport) LPXLOPER WINAPI Local_CMTLEG (LPXLOPER XL_startDate,
													LPXLOPER XL_endDate,
													LPXLOPER XL_cmtType,
													LPXLOPER XL_bondCouponFreq,
													LPXLOPER XL_bondDayCount,
													LPXLOPER XL_receiveOrPay,
													LPXLOPER XL_spread,
													LPXLOPER XL_yieldDecompFreq,
													LPXLOPER XL_swapLegDayCount,
													LPXLOPER XL_intRule,
													LPXLOPER XL_resetGap,
													LPXLOPER XL_resetFreq,
													LPXLOPER XL_ntlAmount,
													LPXLOPER XL_ccy,
													LPXLOPER XL_resetTiming,
													LPXLOPER XL_adjStartDate)
{
	ADD_LOG("Local_CMTLEG ");
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

		CCString C_cmtType;
		long cmtTypeId;

		CCString C_bondCouponFreq;
		long bondCouponFreqId;

		CCString C_bondDayCount;
		long bondDayCountId;

		CCString C_receiveOrPay;
		long receiveOrPayId;

		double C_spread;
		double C_default_spread = 0;
		
		CCString C_yieldDecompFreq;
		long yieldDecompFreqId;

		CCString C_swapLegDayCount;
		long swapLegDayCountId;
		
		CCString C_intRule;
		long intRuleId;

		double C_resetGap;
		double C_resetGap_default = -1;
		
		CCString C_resetFreq;
		long resetFreqId;

		double C_ntlAmount;
		double C_ntlAmount_default = 0;
		
		CCString C_ccy;
		bool ccyIsObject = false;

		CCString C_resetTiming;
		long resetTimingId;

		CCString C_adjStartDate;
		long adjStartDate = 1;

		// error
		static int error;
		static char* reason = "";

		XL_readNumCell(XL_startDate,C_startDate," ARM_ERR: start date: date expected",C_result);
		XL_readNumCell(XL_endDate,C_endDate," ARM_ERR: end date: date expected",C_result);
		XL_readStrCell(XL_cmtType,C_cmtType," ARM_ERR: cmt type: string expected",C_result);
		XL_readStrCellWD(XL_bondCouponFreq,C_bondCouponFreq,"A"," ARM_ERR: bond coupon frequnecy: string expected",C_result);
		XL_readStrCellWD(XL_bondDayCount,C_bondDayCount,"30/360"," ARM_ERR: bond day count: string expected",C_result);
		XL_readStrCellWD(XL_receiveOrPay,C_receiveOrPay,"R"," ARM_ERR: receive or pay: string expected",C_result);
		XL_readNumCellWD(XL_spread,C_spread,C_default_spread," ARM_ERR: spread: numeric expected",C_result);
		XL_readStrCellWD(XL_yieldDecompFreq,C_yieldDecompFreq,"-1"," ARM_ERR: yield decomp frequency: string expected",C_result);
		XL_readStrCellWD(XL_swapLegDayCount,C_swapLegDayCount,"30/360"," ARM_ERR: swap leg day count: string expected",C_result);
		XL_readStrCellWD(XL_intRule,C_intRule,"ADJ"," ARM_ERR: interpolation rule: string expected",C_result);
		XL_readNumCellWD(XL_resetGap,C_resetGap,C_resetGap_default," ARM_ERR: reset gap: numeric expected",C_result);
		XL_readStrCellWD(XL_resetFreq,C_resetFreq,"Q"," ARM_ERR: reset frequency: string expected",C_result);
		XL_readNumCellWD(XL_ntlAmount,C_ntlAmount,C_ntlAmount_default," ARM_ERR: ntl amount: numeric expected",C_result);
		XL_readStrCellWD(XL_ccy,C_ccy,"DEFAULT"," ARM_ERR: currency id: string expected",C_result);
		XL_readStrCellWD(XL_resetTiming,C_resetTiming,"ADV"," ARM_ERR: reset timing: string expected",C_result);
		XL_readStrCellWD(XL_adjStartDate, C_adjStartDate, "YES"," ARM_ERR: adjust Start Date: YES or NO expected", C_result);
		
		cmtTypeId = ARM_ConvIrType (C_cmtType);

		if((bondCouponFreqId = ARM_ConvFrequency (C_bondCouponFreq, C_result)) == ARM_DEFAULT_ERR)
		{
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}

		bondDayCountId = ARM_ConvDayCount (C_bondDayCount);

		if((receiveOrPayId = ARM_ConvRecOrPay (C_receiveOrPay, C_result)) == ARM_DEFAULT_ERR)
		{
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}

		yieldDecompFreqId = ARM_ConvDecompFrequency (C_yieldDecompFreq);

		swapLegDayCountId = ARM_ConvDayCount (C_swapLegDayCount);

		intRuleId = ARM_ConvIntRule (C_intRule);

		if((resetFreqId = ARM_ConvFrequency (C_resetFreq, C_result)) == ARM_DEFAULT_ERR)
		{
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}

		if (( C_ccy.GetLen() > 3 )
			&&
			( !(C_ccy == "DEFAULT") )
		   )
		   ccyIsObject = true;

		resetTimingId = ARM_ConvPayResetRule (C_resetTiming);

		if(C_adjStartDate != CCString("YES"))
			adjStartDate = 0;

		long retCode;
		long objId;
		CCString prevClass;
		
		CCString curClass = LOCAL_CMTLEG_CLASS;
		CCString stringId = GetLastCurCellEnvValue ();
		
		if(!stringId)
		{
			retCode = ARMLOCAL_CMTLEG (C_startDate, C_endDate,
								  cmtTypeId, bondCouponFreqId,
								  bondDayCountId, receiveOrPayId,
								  C_spread, yieldDecompFreqId,
								  swapLegDayCountId, intRuleId,
								  (long)C_resetGap, resetFreqId,
								  C_ntlAmount, ccyIsObject, C_ccy,
								  resetTimingId, adjStartDate, C_result);

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
				retCode = ARMLOCAL_CMTLEG (C_startDate, C_endDate,
									  cmtTypeId, bondCouponFreqId,
									  bondDayCountId, receiveOrPayId,
									  C_spread, yieldDecompFreqId,
									  swapLegDayCountId, intRuleId,
									  (long)C_resetGap, resetFreqId,
									  C_ntlAmount, ccyIsObject, C_ccy,
									  resetTimingId, adjStartDate, C_result, objId);

				if(retCode == ARM_OK)
				{
					LocalSetCurCellEnvValue (curClass, objId); 

					stringId = LocalMakeObjectId (objId, curClass);
				}
			}
			else
			{
				FreeCurCellContent ();
				retCode = ARMLOCAL_CMTLEG (C_startDate, C_endDate,
									  cmtTypeId, bondCouponFreqId,
									  bondDayCountId, receiveOrPayId,
									  C_spread, yieldDecompFreqId,
									  swapLegDayCountId, intRuleId,
									  (long)C_resetGap, resetFreqId,
									  C_ntlAmount, ccyIsObject, C_ccy,
									  resetTimingId, adjStartDate, C_result);

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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_CMTLEG" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


__declspec(dllexport) LPXLOPER WINAPI Local_PXL_CMTLEG (LPXLOPER XL_startDate,
														LPXLOPER XL_endDate,
														LPXLOPER XL_cmtType,
														LPXLOPER XL_bondCouponFreq,
														LPXLOPER XL_bondDayCount,
														LPXLOPER XL_receiveOrPay,
														LPXLOPER XL_spread,
														LPXLOPER XL_yieldDecompFreq,
														LPXLOPER XL_swapLegDayCount,
														LPXLOPER XL_intRule,
														LPXLOPER XL_resetGap,
														LPXLOPER XL_resetFreq,
														LPXLOPER XL_ntlAmount,
														LPXLOPER XL_ccy,
														LPXLOPER XL_resetTiming,
														LPXLOPER XL_adjStartDate)
{
	ADD_LOG("Local_PXL_CMTLEG ");
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

		CCString C_cmtType;
		long cmtTypeId;

		CCString C_bondCouponFreq;
		long bondCouponFreqId;

		CCString C_bondDayCount;
		long bondDayCountId;
			
		CCString C_receiveOrPay;
		long receiveOrPayId;

		double C_spread;
		double C_default_spread = 0;
		
		CCString C_yieldDecompFreq;
		long yieldDecompFreqId;

		CCString C_swapLegDayCount;
		long swapLegDayCountId;
		
		CCString C_intRule;
		long intRuleId;

		double C_resetGap;
		double C_resetGap_default = -1;
		
		CCString C_resetFreq;
		long resetFreqId;

		double C_ntlAmount;
		double C_ntlAmount_default = 0;
		
		CCString C_ccy;
		bool ccyIsObject = false;
				
		CCString C_resetTiming;
		long resetTimingId;

		CCString C_adjStartDate;
		long adjStartDate = 1;

		// error
		static int error;
		static char* reason = "";

		XL_readNumCell(XL_startDate,C_startDate," ARM_ERR: start date: date expected",C_result);
		XL_readNumCell(XL_endDate,C_endDate," ARM_ERR: end date: date expected",C_result);
		XL_readStrCell(XL_cmtType,C_cmtType," ARM_ERR: cmt type: string expected",C_result);
		XL_readStrCellWD(XL_bondCouponFreq,C_bondCouponFreq,"A"," ARM_ERR: bond coupon frequnecy: string expected",C_result);
		XL_readStrCellWD(XL_bondDayCount,C_bondDayCount,"30/360"," ARM_ERR: bond day count: string expected",C_result);
		XL_readStrCellWD(XL_receiveOrPay,C_receiveOrPay,"R"," ARM_ERR: receive or pay: string expected",C_result);
		XL_readNumCellWD(XL_spread,C_spread,C_default_spread," ARM_ERR: spread: numeric expected",C_result);
		XL_readStrCellWD(XL_yieldDecompFreq,C_yieldDecompFreq,"-1"," ARM_ERR: yield decomp frequency: string expected",C_result);
		XL_readStrCellWD(XL_swapLegDayCount,C_swapLegDayCount,"30/360"," ARM_ERR: swap leg day count: string expected",C_result);
		XL_readStrCellWD(XL_intRule,C_intRule,"ADJ"," ARM_ERR: interpolation rule: string expected",C_result);
		XL_readNumCellWD(XL_resetGap,C_resetGap,C_resetGap_default," ARM_ERR: reset gap: numeric expected",C_result);
		XL_readStrCellWD(XL_resetFreq,C_resetFreq,"Q"," ARM_ERR: reset frequency: string expected",C_result);
		XL_readNumCellWD(XL_ntlAmount,C_ntlAmount,C_ntlAmount_default," ARM_ERR: ntl amount: numeric expected",C_result);
		XL_readStrCellWD(XL_ccy,C_ccy,"DEFAULT"," ARM_ERR: currency id: string expected",C_result);
		XL_readStrCellWD(XL_resetTiming,C_resetTiming,"ADV"," ARM_ERR: reset timing: string expected",C_result);
		XL_readStrCellWD(XL_adjStartDate, C_adjStartDate, "YES", " ARM_ERR: adjust Start Date : YES or NO expected", C_result);
		
		cmtTypeId = ARM_ConvIrType (C_cmtType);

		if((bondCouponFreqId = ARM_ConvFrequency (C_bondCouponFreq, C_result)) == ARM_DEFAULT_ERR)
		{
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}

		bondDayCountId = ARM_ConvDayCount (C_bondDayCount);

		if((receiveOrPayId = ARM_ConvRecOrPay (C_receiveOrPay, C_result)) == ARM_DEFAULT_ERR)
		{
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}

		yieldDecompFreqId = ARM_ConvDecompFrequency (C_yieldDecompFreq);

		swapLegDayCountId = ARM_ConvDayCount (C_swapLegDayCount);

		intRuleId = ARM_ConvIntRule (C_intRule);

		if((resetFreqId = ARM_ConvFrequency (C_resetFreq, C_result)) == ARM_DEFAULT_ERR)
		{
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}

		if (( C_ccy.GetLen() > 3 )
			&&
			( !(C_ccy == "DEFAULT") )
		   )
		   ccyIsObject = true;

		resetTimingId = ARM_ConvPayResetRule (C_resetTiming);

		if(C_adjStartDate != CCString("YES"))
			adjStartDate = 0;

		long retCode;
		long objId;
		
		CCString curClass = LOCAL_CMTLEG_CLASS;
		CCString stringId;
		
		retCode = ARMLOCAL_CMTLEG (C_startDate, C_endDate,
								  cmtTypeId, bondCouponFreqId,
								  bondDayCountId, receiveOrPayId,
								  C_spread, yieldDecompFreqId,
								  swapLegDayCountId, intRuleId,
								  (long)C_resetGap, resetFreqId,
								  C_ntlAmount, ccyIsObject, C_ccy,
								  resetTimingId, adjStartDate, C_result);

		if(retCode == ARM_OK)
		{
			objId = C_result.getLong ();

			stringId = LocalMakeObjectId (objId, curClass);
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_PXL_CMTLEG" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


__declspec(dllexport) LPXLOPER WINAPI Local_TMLEG(LPXLOPER XL_tmIxType,
												  LPXLOPER XL_startDate,
												  LPXLOPER XL_endDate,
												  LPXLOPER XL_receiveOrPay,
												  LPXLOPER XL_spread,
                                                  LPXLOPER XL_payFreq,
                                                  LPXLOPER XL_resetFreq,
                                                  LPXLOPER XL_intRule,
                                                  LPXLOPER XL_fwdRule,
												  LPXLOPER XL_stubRule,
                                                  LPXLOPER XL_ccy)
{
	ADD_LOG("Local_TMLEG");
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
		ARM_NOCALCIFWIZ();

		// C variable
		CCString C_tmIxType;
		long tmIxTypeId;

		double C_startDate;
		double C_endDate;
		
		CCString C_receiveOrPay;
		double receiveOrPayId;

		double C_spread;
		double C_default_spread = 0;

        CCString C_payFreq;
		long payFreqId;

        CCString C_resetFreq;
		long resetFreqId;

        CCString C_ccy;
		bool ccyIsObject = false;

		CCString C_intRule;
		long intRuleId;

        CCString C_fwdRule;
		long fwdRuleId;
		
		CCString C_stubRule;
		long stubRuleId;

		// error
		static int error;
		static char* reason = "";

		XL_readStrCell(XL_tmIxType,C_tmIxType," ARM_ERR: tmIxtype: string expected",C_result);
		XL_readNumCell(XL_startDate,C_startDate," ARM_ERR: start date: date expected",C_result);
		XL_readNumCell(XL_endDate,C_endDate," ARM_ERR: end date: date expected",C_result);
		XL_readStrCell(XL_receiveOrPay,C_receiveOrPay," ARM_ERR: receive or pay: string expected",C_result);
		XL_readNumCellWD(XL_spread,C_spread,C_default_spread," ARM_ERR: spread: numeric expected",C_result);
		XL_readStrCellWD(XL_payFreq, C_payFreq, "A"," ARM_ERR: pay frequency: string expected",C_result);
        XL_readStrCellWD(XL_resetFreq, C_resetFreq, "A"," ARM_ERR: reset frequency: string expected",C_result);
        
		XL_readStrCellWD(XL_intRule,C_intRule,"ADJ"," ARM_ERR: intRule: string expected",C_result);
        XL_readStrCellWD(XL_fwdRule,C_fwdRule,"MF"," ARM_ERR: fwdRule: string expected",C_result);
        XL_readStrCellWD(XL_stubRule,C_stubRule,"SS"," ARM_ERR: stubRule: string expected",C_result);
        XL_readStrCellWD(XL_ccy, C_ccy, "DEFAULT"," ARM_ERR: currency id: string expected",C_result);

		if ((receiveOrPayId = ARM_ConvRecOrPay (C_receiveOrPay, C_result)) == ARM_DEFAULT_ERR )
		{
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}

		tmIxTypeId = ARM_ConvIrType(C_tmIxType);

        if ((payFreqId = ARM_ConvFrequency (C_payFreq, C_result)) == ARM_DEFAULT_ERR)
		{
			ARM_ARG_ERR();
		
            return (LPXLOPER)&XL_result;
		}

		if ((resetFreqId = ARM_ConvFrequency (C_resetFreq, C_result)) == ARM_DEFAULT_ERR)
		{
			ARM_ARG_ERR();
			
            return (LPXLOPER)&XL_result;
		}

        // interest Rule
        intRuleId = ARM_ConvIntRule(C_intRule);

        // Fwd Rule

        if ( (fwdRuleId = ARM_ConvFwdRule (C_fwdRule, C_result)) == ARM_DEFAULT_ERR )
		{
			ARM_ARG_ERR();
		
            return (LPXLOPER)&XL_result;
		}
		
		// stubRule
        stubRuleId = ARM_ConvStubRule(C_stubRule);

        if (( C_ccy.GetLen() > 3 )
			&& 
			( !(C_ccy == "DEFAULT"))
		   )
		   ccyIsObject = true;

		long retCode;
		long objId;
		CCString prevClass;
		
		CCString curClass = LOCAL_T4MLEG_CLASS;
		CCString stringId = GetLastCurCellEnvValue ();
		
		if(!stringId)
		{
			retCode = ARMLOCAL_TMLEG(tmIxTypeId, C_startDate, C_endDate, receiveOrPayId, 
                                     C_spread,
                                     payFreqId,
                                     resetFreqId,
                                     intRuleId,
                                     fwdRuleId,
									 stubRuleId,
                                     ccyIsObject,
								     C_ccy,
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
				retCode = ARMLOCAL_TMLEG(tmIxTypeId, C_startDate, C_endDate, receiveOrPayId, 
                                         C_spread,
                                         payFreqId,
                                         resetFreqId,
                                         intRuleId,
                                         fwdRuleId,
										 stubRuleId,
                                         ccyIsObject,
								         C_ccy,
                                         C_result, objId);
			}
			else
			{
				FreeCurCellContent ();
				
                retCode = ARMLOCAL_TMLEG(tmIxTypeId, C_startDate, C_endDate, receiveOrPayId, 
                                         C_spread,
                                         payFreqId,
                                         resetFreqId,
                                         intRuleId,
                                         fwdRuleId,
										 stubRuleId,
                                         ccyIsObject,
								         C_ccy,
                                         C_result);
			
				if ( retCode == ARM_OK )
				{
					objId = C_result.getLong ();

					LocalSetCurCellEnvValue (curClass, objId); 

					stringId = LocalMakeObjectId (objId, curClass);
				}
			}
		}

		if ( retCode == ARM_OK )
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_TMLEG" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}



__declspec(dllexport) LPXLOPER WINAPI Local_PXL_TMLEG (LPXLOPER XL_tmIxType,
													   LPXLOPER XL_startDate,
													   LPXLOPER XL_endDate,
													   LPXLOPER XL_receiveOrPay,
													   LPXLOPER XL_spread,
                                                       LPXLOPER XL_payFreq,
                                                       LPXLOPER XL_resetFreq,
                                                       LPXLOPER XL_intRule,
                                                       LPXLOPER XL_fwdRule,
													   LPXLOPER XL_stubRule,
                                                       LPXLOPER XL_ccy)
{
	ADD_LOG("Local_PXL_TMLEG ");
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
		ARM_NOCALCIFWIZ();

		// C variable
		CCString C_tmIxType;
		long tmIxTypeId;

		double C_startDate;
		double C_endDate;
		
		CCString C_receiveOrPay;
		double receiveOrPayId;

		double C_spread;
		double C_default_spread = 0;

        CCString C_payFreq;
		long payFreqId;

        CCString C_resetFreq;
		long resetFreqId;

        CCString C_ccy;
		bool ccyIsObject = false;

		CCString C_intRule;
		long intRuleId;

        CCString C_fwdRule;
		long fwdRuleId;
	
		CCString C_stubRule;
		long stubRuleId;

		// error
		static int error;
		static char* reason = "";

		XL_readStrCell(XL_tmIxType,C_tmIxType," ARM_ERR: tmIxtype: string expected",C_result);
		XL_readNumCell(XL_startDate,C_startDate," ARM_ERR: start date: date expected",C_result);
		XL_readNumCell(XL_endDate,C_endDate," ARM_ERR: end date: date expected",C_result);
		XL_readStrCell(XL_receiveOrPay,C_receiveOrPay," ARM_ERR: receive or pay: string expected",C_result);
		XL_readNumCellWD(XL_spread,C_spread,C_default_spread," ARM_ERR: spread: numeric expected",C_result);
		
        XL_readStrCellWD(XL_payFreq, C_payFreq, "A"," ARM_ERR: pay frequency: string expected",C_result);
        XL_readStrCellWD(XL_resetFreq, C_resetFreq, "A"," ARM_ERR: reset frequency: string expected",C_result);

        XL_readStrCellWD(XL_intRule,C_intRule,"ADJ"," ARM_ERR: intRule: string expected",C_result);
        XL_readStrCellWD(XL_fwdRule,C_fwdRule,"MF"," ARM_ERR: fwdRule: string expected",C_result);
        XL_readStrCellWD(XL_stubRule,C_stubRule,"SS"," ARM_ERR: stubRule: string expected",C_result);
        XL_readStrCellWD(XL_ccy, C_ccy, "DEFAULT"," ARM_ERR: currency id: string expected",C_result);

		if((receiveOrPayId = ARM_ConvRecOrPay (C_receiveOrPay, C_result)) == ARM_DEFAULT_ERR)
		{
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}

		tmIxTypeId = ARM_ConvIrType (C_tmIxType);

        if ((payFreqId = ARM_ConvFrequency (C_payFreq, C_result)) == ARM_DEFAULT_ERR)
		{
			ARM_ARG_ERR();
		
            return (LPXLOPER)&XL_result;
		}

		if ((resetFreqId = ARM_ConvFrequency (C_resetFreq, C_result)) == ARM_DEFAULT_ERR)
		{
			ARM_ARG_ERR();
			
            return (LPXLOPER)&XL_result;
		}

        // interest Rule
        intRuleId = ARM_ConvIntRule(C_intRule);

        // Fwd Rule

        if ( (fwdRuleId = ARM_ConvFwdRule (C_fwdRule, C_result)) == ARM_DEFAULT_ERR )
		{
			ARM_ARG_ERR();
		
            return (LPXLOPER)&XL_result;
		}

		// stubRule
        stubRuleId = ARM_ConvStubRule(C_stubRule);

        if (( C_ccy.GetLen() > 3 )
			&& 
			( !(C_ccy == "DEFAULT"))
		   )
		   ccyIsObject = true;

		long retCode;
		long objId;
		
		CCString curClass = LOCAL_T4MLEG_CLASS;
		CCString stringId;
		
		retCode = ARMLOCAL_TMLEG(tmIxTypeId, C_startDate, C_endDate, receiveOrPayId, 
                                 C_spread,
                                 payFreqId,
                                 resetFreqId,
                                 intRuleId,
                                 fwdRuleId,
								 stubRuleId,
                                 ccyIsObject,
								 C_ccy,
                                 C_result);

		if ( retCode == ARM_OK )
		{
			objId = C_result.getLong ();

			stringId = LocalMakeObjectId (objId, curClass);
		}

		if ( retCode == ARM_OK )
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_PXL_TMLEG" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


__declspec(dllexport) LPXLOPER WINAPI Local_REVERSEFLOAT(LPXLOPER XL_startDate,
														 LPXLOPER XL_endDate,
														 LPXLOPER XL_NotionalRatio,
														 LPXLOPER XL_receiveOrPay,
														 LPXLOPER XL_liborType,
														 LPXLOPER XL_fltCcy,
														 LPXLOPER XL_fltSpreadsId,
														 LPXLOPER XL_fltDayCount,
														 LPXLOPER XL_fxCouponsId,
														 LPXLOPER XL_extraCouponsId,
														 LPXLOPER XL_reverseCcy,
														 LPXLOPER XL_payFreq,
														 LPXLOPER XL_reverseIndexId,
														 LPXLOPER XL_fxDayCount,
														 LPXLOPER XL_multiplier,
														 LPXLOPER XL_couponFloor,
														 LPXLOPER XL_StubDate)
{
	ADD_LOG("Local_REVERSEFLOAT");
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
		ARM_NOCALCIFWIZ();

		// C variable
		double C_startDate;
		double C_endDate;
		double C_stubDate;
		double C_NotionnalRatio;

		CCString C_receiveOrPay;
		long receiveOrPayId;
		
		CCString C_liborType;
		long liborTypeId;

		CCString C_fltCcy;

		CCString C_fltSpreads;    

		CCString C_fltDc;
		long fltDcId;

		CCString C_fxCoupons;
		CCString C_extraCoupons;

		CCString C_reverseCcy;

  		CCString C_payFreq;
		long payFreqId;

  		CCString C_reverseIndexId;
		long reverseIndexId;

  		CCString C_fxDayCount;
		long fxDayCount;

		double C_multiplier;
		double C_couponFloor;

		double one  = 1.0;
		double zero = 0.0;
		double minusOne = -1.0;
		
		// error
		static int error;
		static char* reason = "";

		XL_readNumCell(XL_startDate,C_startDate," ARM_ERR: start date: date expected",C_result);
		XL_readNumCell(XL_endDate,C_endDate," ARM_ERR: end date: date expected",C_result);
		XL_readNumCell(XL_NotionalRatio,C_NotionnalRatio," ARM_ERR: notional ratio : numeric expected",C_result);
		XL_readStrCell(XL_receiveOrPay,C_receiveOrPay," ARM_ERR: receive or pay: string expected",C_result);
		XL_readStrCell(XL_liborType,C_liborType," ARM_ERR: libor type: string expected",C_result);
		XL_readStrCell(XL_fltCcy, C_fltCcy," ARM_ERR: floating currency: string expected",C_result);
		XL_readStrCell(XL_fltSpreadsId,C_fltSpreads," ARM_ERR: floating spreads : object expected",C_result);
		XL_readStrCell(XL_fltDayCount,C_fltDc," ARM_ERR: floating day count : string expected",C_result);
		XL_readStrCell(XL_fxCouponsId,C_fxCoupons," ARM_ERR: fixed coupons : object expected",C_result);
		XL_readStrCell(XL_extraCouponsId,C_extraCoupons," ARM_ERR: extra coupons : object expected",C_result);
		XL_readStrCell(XL_reverseCcy, C_reverseCcy," ARM_ERR: reverse currency: string expected",C_result);    
		XL_readStrCell(XL_payFreq,C_payFreq," ARM_ERR: pay frequency: string expected",C_result);
		XL_readStrCell(XL_reverseIndexId,C_reverseIndexId," ARM_ERR: reverse index type: string expected",C_result);
		XL_readStrCell(XL_fxDayCount,C_fxDayCount," ARM_ERR: fix day count : string expected",C_result);
		XL_readNumCellWD(XL_multiplier,C_multiplier, one, " ARM_ERR: multiplier : numeric expected",C_result);
		XL_readNumCellWD(XL_couponFloor,C_couponFloor, zero, " ARM_ERR: coupon floor: numeric expected",C_result);
		XL_readNumCellWD(XL_StubDate,C_stubDate, minusOne, " ARM_ERR: stubDate : date expected",C_result);

		if(C_fltCcy == "DEFAULT")
		{
			ARM_result currencyres;
			ARMLOCAL_ARM_GetDefaultCurrency (currencyres);
			if(currencyres.getRetCode () != ARM_OK)
			{
				ARM_ARG_ERR();
				return (LPXLOPER)&XL_result;
			}
			else
			{
				C_fltCcy = currencyres.getString ();
			}
		}

		if(C_reverseCcy == "DEFAULT")
		{
			ARM_result currencyres;
			ARMLOCAL_ARM_GetDefaultCurrency (currencyres);
			if(currencyres.getRetCode () != ARM_OK)
			{
				ARM_ARG_ERR();
				return (LPXLOPER)&XL_result;
			}
			else
			{
				C_reverseCcy = currencyres.getString ();
			}
		}

		if((receiveOrPayId = ARM_ConvRecOrPay (C_receiveOrPay, C_result)) == ARM_DEFAULT_ERR)
		{
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}

		if((liborTypeId = ARM_ConvIrIndName (C_liborType, C_result)) == ARM_DEFAULT_ERR)
		{
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}

		if((reverseIndexId = ARM_ConvIrIndName (C_reverseIndexId, C_result)) == ARM_DEFAULT_ERR)
		{
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}

		fltDcId  = ARM_ConvDayCount (C_fltDc);

		if((payFreqId = ARM_ConvFrequency (C_payFreq, C_result)) == ARM_DEFAULT_ERR)
		{
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}

		fxDayCount  = ARM_ConvDayCount (C_fxDayCount);

		long retCode;
		long objId;
		CCString prevClass;
		
		CCString curClass = LOCAL_REVERSEFLOATER_CLASS;
		CCString stringId = GetLastCurCellEnvValue ();
		
		if(!stringId)
		{
			retCode = ARMLOCAL_REVERSEFLOAT ( C_startDate, C_endDate, C_NotionnalRatio, receiveOrPayId,
								liborTypeId, C_fltCcy, LocalGetNumObjectId(C_fltSpreads),
								fltDcId, LocalGetNumObjectId(C_fxCoupons),
								LocalGetNumObjectId(C_extraCoupons),  C_reverseCcy,
								payFreqId, reverseIndexId, fxDayCount, C_multiplier,
								C_couponFloor, C_stubDate, C_result);

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
				retCode = ARMLOCAL_REVERSEFLOAT ( C_startDate, C_endDate, C_NotionnalRatio, receiveOrPayId,
								liborTypeId, C_fltCcy, LocalGetNumObjectId(C_fltSpreads),
								fltDcId, LocalGetNumObjectId(C_fxCoupons), 
								LocalGetNumObjectId(C_extraCoupons), C_reverseCcy,
								payFreqId, reverseIndexId, fxDayCount, C_multiplier,
								C_couponFloor, C_stubDate, C_result,objId);

				if(retCode == ARM_OK)
				{
					LocalSetCurCellEnvValue (curClass, objId); 

					stringId = LocalMakeObjectId (objId, curClass);
				}
			}
			else
			{
				FreeCurCellContent ();
				retCode = ARMLOCAL_REVERSEFLOAT ( C_startDate, C_endDate, C_NotionnalRatio, receiveOrPayId,
								liborTypeId, C_fltCcy, LocalGetNumObjectId(C_fltSpreads),
								fltDcId, LocalGetNumObjectId(C_fxCoupons),
								LocalGetNumObjectId(C_extraCoupons), C_reverseCcy,
								payFreqId, reverseIndexId, fxDayCount, C_multiplier,
								C_couponFloor, C_stubDate, C_result);
			
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_REVERSEFLOAT" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


__declspec(dllexport) LPXLOPER WINAPI Local_PXL_REVERSEFLOAT(LPXLOPER XL_startDate,
															 LPXLOPER XL_endDate,
															 LPXLOPER XL_NotionalRatio,
															 LPXLOPER XL_receiveOrPay,
															 LPXLOPER XL_liborType,
															 LPXLOPER XL_fltCcy,
															 LPXLOPER XL_fltSpreadsId,
															 LPXLOPER XL_fltDayCount,
															 LPXLOPER XL_fxCouponsId,
															 LPXLOPER XL_extraCouponsId,
															 LPXLOPER XL_reverseCcy,
															 LPXLOPER XL_payFreq,
															 LPXLOPER XL_reverseIndexId,
															 LPXLOPER XL_fxDayCount,
															 LPXLOPER XL_multiplier,
															 LPXLOPER XL_couponFloor,
															 LPXLOPER XL_StubDate)
{
	ADD_LOG("Local_PXL_REVERSEFLOAT");
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
		ARM_NOCALCIFWIZ();

		// C variable
		double C_startDate;
		double C_endDate;
		double C_stubDate;
		double C_NotionnalRatio;

		CCString C_receiveOrPay;
		long receiveOrPayId;
		
		CCString C_liborType;
		long liborTypeId;

		CCString C_fltCcy;

		CCString C_fltSpreads;    

		CCString C_fltDc;
		long fltDcId;

		CCString C_fxCoupons;
		CCString C_extraCoupons;    

		CCString C_reverseCcy;

  		CCString C_payFreq;
		long payFreqId;

  		CCString C_reverseIndexId;
		long reverseIndexId;

  		CCString C_fxDayCount;
		long fxDayCount;

		double C_multiplier;
		double C_couponFloor;

		double one  = 1.0;
		double zero = 0.0;
		double minusOne = -1.0;
		
		// error
		static int error;
		static char* reason = "";

		XL_readNumCell(XL_startDate,C_startDate," ARM_ERR: start date: date expected",C_result);
		XL_readNumCell(XL_endDate,C_endDate," ARM_ERR: end date: date expected",C_result);
		XL_readNumCell(XL_NotionalRatio,C_NotionnalRatio," ARM_ERR: notional ratio : numeric expected",C_result);
		XL_readStrCell(XL_receiveOrPay,C_receiveOrPay," ARM_ERR: receive or pay: string expected",C_result);
		XL_readStrCell(XL_liborType,C_liborType," ARM_ERR: libor type: string expected",C_result);
		XL_readStrCell(XL_fltCcy, C_fltCcy," ARM_ERR: floating currency: string expected",C_result);
		XL_readStrCell(XL_fltSpreadsId,C_fltSpreads," ARM_ERR: floating spreads : object expected",C_result);
		XL_readStrCell(XL_fltDayCount,C_fltDc," ARM_ERR: floating day count : string expected",C_result);
		XL_readStrCell(XL_fxCouponsId,C_fxCoupons," ARM_ERR: fixed coupons : object expected",C_result);
		XL_readStrCell(XL_extraCouponsId,C_extraCoupons," ARM_ERR: extra coupons : object expected",C_result);    
		XL_readStrCell(XL_reverseCcy, C_reverseCcy," ARM_ERR: reverse currency: string expected",C_result);    
		XL_readStrCell(XL_payFreq,C_payFreq," ARM_ERR: pay frequency: string expected",C_result);
		XL_readStrCell(XL_reverseIndexId,C_reverseIndexId," ARM_ERR: reverse index type: string expected",C_result);
		XL_readStrCell(XL_fxDayCount,C_fxDayCount," ARM_ERR: fix day count : string expected",C_result);
		XL_readNumCellWD(XL_multiplier,C_multiplier, one, " ARM_ERR: multiplier : numeric expected",C_result);
		XL_readNumCellWD(XL_couponFloor,C_couponFloor, zero, " ARM_ERR: coupon floor: numeric expected",C_result);
		XL_readNumCellWD(XL_StubDate,C_stubDate, minusOne, " ARM_ERR: stubDate : date expected",C_result);

		if(C_fltCcy == "DEFAULT")
		{
			ARM_result currencyres;
			ARMLOCAL_ARM_GetDefaultCurrency (currencyres);
			if(currencyres.getRetCode () != ARM_OK)
			{
				ARM_ARG_ERR();
				return (LPXLOPER)&XL_result;
			}
			else
			{
				C_fltCcy = currencyres.getString ();
			}
		}

		if(C_reverseCcy == "DEFAULT")
		{
			ARM_result currencyres;
			ARMLOCAL_ARM_GetDefaultCurrency (currencyres);
			if(currencyres.getRetCode () != ARM_OK)
			{
				ARM_ARG_ERR();
				return (LPXLOPER)&XL_result;
			}
			else
			{
				C_reverseCcy = currencyres.getString ();
			}
		}

		if((receiveOrPayId = ARM_ConvRecOrPay (C_receiveOrPay, C_result)) == ARM_DEFAULT_ERR)
		{
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}

		if((liborTypeId = ARM_ConvIrIndName (C_liborType, C_result)) == ARM_DEFAULT_ERR)
		{
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}

		if((reverseIndexId = ARM_ConvIrIndName (C_reverseIndexId, C_result)) == ARM_DEFAULT_ERR)
		{
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}

		fltDcId  = ARM_ConvDayCount (C_fltDc);

		if((payFreqId = ARM_ConvFrequency (C_payFreq, C_result)) == ARM_DEFAULT_ERR)
		{
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}

		fxDayCount  = ARM_ConvDayCount (C_fxDayCount);

		long retCode;
		long objId;
		
		CCString curClass = LOCAL_REVERSEFLOATER_CLASS;
		CCString stringId;
		
		retCode = ARMLOCAL_REVERSEFLOAT ( C_startDate, C_endDate, C_NotionnalRatio, receiveOrPayId,
								liborTypeId, C_fltCcy, LocalGetNumObjectId(C_fltSpreads),
								fltDcId, LocalGetNumObjectId(C_fxCoupons),
								LocalGetNumObjectId(C_extraCoupons),
								C_reverseCcy,
								payFreqId, reverseIndexId, fxDayCount, C_multiplier,
								C_couponFloor, C_stubDate, C_result);

		if(retCode == ARM_OK)
		{
			objId = C_result.getLong ();

			stringId = LocalMakeObjectId (objId, curClass);
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_PXL_REVERSEFLOAT" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;}




__declspec(dllexport) LPXLOPER WINAPI Local_IMPLIEDSPREAD (LPXLOPER XL_swap,
														   LPXLOPER XL_model,
														   LPXLOPER XL_price,
														   LPXLOPER XL_leg1Or2)
{
	ADD_LOG("Local_IMPLIEDSPREAD ");
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
		ARM_NOCALCIFWIZ();

		// C variable
		CCString C_swap;
		CCString C_model;
		double C_price;
		double C_leg1Or2;
		
		// error
		static int error;
		static char* reason = "";

		XL_readStrCell(XL_swap,C_swap," ARM_ERR: swap id: object expected",C_result);
		XL_readStrCell(XL_model,C_model," ARM_ERR: model id: object expected",C_result);
		XL_readNumCell(XL_price,C_price," ARM_ERR: price: numeric expected",C_result);
		XL_readNumCell(XL_leg1Or2,C_leg1Or2," ARM_ERR: leg 1 or 2: numeric expected",C_result);

		long retCode = ARMLOCAL_IMPLIEDSPREAD (LocalGetNumObjectId (C_swap),
											   LocalGetNumObjectId (C_model),
											   C_price,
											   (long)C_leg1Or2,
											   C_result);

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
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_IMPLIEDSPREAD" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}

__declspec(dllexport) LPXLOPER WINAPI Local_ARM_GenAmortization(LPXLOPER XL_swapleg,
																LPXLOPER XL_amortMethod,
																LPXLOPER XL_amortFrequency,
																LPXLOPER XL_amortAmount,
																LPXLOPER XL_basisDaycount,
																LPXLOPER XL_legNotional,
																LPXLOPER XL_amortRate,
																LPXLOPER XL_reducedMaturity,
																LPXLOPER XL_modelId,
																LPXLOPER XL_cleanup)
{
	ADD_LOG("Local_ARM_GenAmortization");
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
		ARM_NOCALCIFWIZ();

		// C variable
		CCString C_swapleg;
		
		CCString C_amortMethod;
		long amortMethodId;

		CCString C_amortFrequency;
		long freqId;

		double C_amortAmount;
		double C_amortAmount_default (0.);

		CCString C_basisDaycount;
		long dayCountId;

		double C_legNotional;
		double C_legNotional_default(-9999);

		double C_amortRate;
		double C_amortRate_default(0.);

		double C_reducedMaturity;
		double C_reducedMaturity_default(0.);

		CCString C_model;
		long C_modelId;

		double C_cleanup;
		double C_cleanup_default(0.);

		// error
		static int error;
		static char* reason = "";

		XL_readStrCell(XL_swapleg,C_swapleg," ARM_ERR: swapleg id: object expected",C_result);
		XL_readStrCell(XL_amortMethod,C_amortMethod," ARM_ERR: amortization Method: string expected",C_result);
		XL_readStrCellWD(XL_amortFrequency,C_amortFrequency,"-1"," ARM_ERR: amortization Frequency: string expected",C_result);
		XL_readNumCellWD(XL_amortAmount,C_amortAmount,C_amortAmount_default," ARM_ERR: amortization Amount: numeric expected",C_result);
		XL_readStrCellWD(XL_basisDaycount,C_basisDaycount,"-1"," ARM_ERR: basis daycount: string expected",C_result);
		XL_readNumCellWD(XL_legNotional,C_legNotional,C_legNotional_default," ARM_ERR: leg Notional: numeric expected",C_result);
		XL_readNumCellWD(XL_amortRate,C_amortRate,C_amortRate_default," ARM_ERR: amortization rate: numeric expected",C_result);
		XL_readNumCellWD(XL_reducedMaturity,C_reducedMaturity,C_reducedMaturity_default," ARM_ERR: reduced Maturity: numeric expected",C_result);
		XL_readStrCellWD(XL_modelId,C_model,"DEFAULT"," ARM_ERR: model Id: string expected",C_result);
		XL_readNumCellWD(XL_cleanup,C_cleanup,C_cleanup_default," ARM_ERR: clean up: numeric expected",C_result);

		dayCountId = ARM_ConvDayCount (C_basisDaycount);

		if((amortMethodId = ARM_ConvAmortMethod (C_amortMethod, C_result)) == ARM_DEFAULT_ERR)
		{
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}

		if((freqId = ARM_ConvFrequency (C_amortFrequency, C_result)) == ARM_DEFAULT_ERR)
		{
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}

		if(C_model == "DEFAULT")
		{
			C_modelId = ARM_NULL_OBJECT;
		}
		else
		{
			C_modelId = LocalGetNumObjectId (C_model);
		}

		long retCode;
		long objId;
		CCString prevClass;
		
		CCString curClass = LOCAL_REFVAL_CLASS;
		CCString stringId = GetLastCurCellEnvValue ();
		
		if(!stringId)
		{
			retCode = ARMLOCAL_ARM_GenAmortization (LocalGetNumObjectId (C_swapleg),
													amortMethodId,
													freqId,
													C_amortAmount,
													dayCountId,
													C_legNotional,
													C_amortRate,
													C_reducedMaturity,
													C_modelId,
													C_cleanup,
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
				retCode = ARMLOCAL_ARM_GenAmortization (LocalGetNumObjectId (C_swapleg),
														amortMethodId,
														freqId,
														C_amortAmount,
														dayCountId,
														C_legNotional,
														C_amortRate,
														C_reducedMaturity,
														C_modelId,
														C_cleanup,
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
				retCode = ARMLOCAL_ARM_GenAmortization (LocalGetNumObjectId (C_swapleg),
														amortMethodId,
														freqId,
														C_amortAmount,
														dayCountId,
														C_legNotional,
														C_amortRate,
														C_reducedMaturity,
														C_modelId,
														C_cleanup,
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_ARM_GenAmortization" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


__declspec(dllexport) LPXLOPER WINAPI Local_PXL_ARM_GenAmortization(LPXLOPER XL_swapleg,
																	LPXLOPER XL_amortMethod,
																	LPXLOPER XL_amortFrequency,
																	LPXLOPER XL_amortAmount,
																	LPXLOPER XL_basisDaycount,
																	LPXLOPER XL_legNotional,
																	LPXLOPER XL_amortRate,
																	LPXLOPER XL_reducedMaturity,
																	LPXLOPER XL_modelId,
																	LPXLOPER XL_cleanup)
{
	ADD_LOG("Local_PXL_ARM_GenAmortization");
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
		ARM_NOCALCIFWIZ();

		// C variable
		CCString C_swapleg;
		
		CCString C_amortMethod;
		long amortMethodId;

		CCString C_amortFrequency;
		long freqId;

		double C_amortAmount;
		double C_amortAmount_default (0.);

		CCString C_basisDaycount;
		long dayCountId;

		double C_legNotional;
		double C_legNotional_default(-9999);

		double C_amortRate;
		double C_amortRate_default(0.);

		double C_reducedMaturity;
		double C_reducedMaturity_default(0.);

		CCString C_model;
		long C_modelId;

		double C_cleanup;
		double C_cleanup_default(0.);

		// error
		static int error;
		static char* reason = "";

		XL_readStrCell(XL_swapleg,C_swapleg," ARM_ERR: swapleg id: object expected",C_result);
		XL_readStrCell(XL_amortMethod,C_amortMethod," ARM_ERR: amortization Method: string expected",C_result);
		XL_readStrCellWD(XL_amortFrequency,C_amortFrequency,"-1"," ARM_ERR: amortization Frequency: string expected",C_result);
		XL_readNumCellWD(XL_amortAmount,C_amortAmount,C_amortAmount_default," ARM_ERR: amortization Amount: numeric expected",C_result);
		XL_readStrCellWD(XL_basisDaycount,C_basisDaycount,"-1"," ARM_ERR: basis daycount: string expected",C_result);
		XL_readNumCellWD(XL_legNotional,C_legNotional,C_legNotional_default," ARM_ERR: leg Notional: numeric expected",C_result);
		XL_readNumCellWD(XL_amortRate,C_amortRate,C_amortRate_default," ARM_ERR: amortization rate: numeric expected",C_result);
		XL_readNumCellWD(XL_reducedMaturity,C_reducedMaturity,C_reducedMaturity_default," ARM_ERR: reduced Maturity: numeric expected",C_result);
		XL_readStrCellWD(XL_modelId,C_model,"DEFAULT"," ARM_ERR: reduced Maturity: numeric expected",C_result);
		XL_readNumCellWD(XL_cleanup,C_cleanup,C_cleanup_default," ARM_ERR: clean up: numeric expected",C_result);

		dayCountId = ARM_ConvDayCount (C_basisDaycount);

		if((amortMethodId = ARM_ConvAmortMethod (C_amortMethod, C_result)) == ARM_DEFAULT_ERR)
		{
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}

		if((freqId = ARM_ConvFrequency (C_amortFrequency, C_result)) == ARM_DEFAULT_ERR)
		{
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}

		if(C_model == "DEFAULT")
		{
			C_modelId = ARM_NULL_OBJECT;
		}
		else
		{
			C_modelId = LocalGetNumObjectId (C_model);
		}

		long retCode;
		long objId;
		
		CCString curClass = LOCAL_REFVAL_CLASS;
		CCString stringId;


		retCode = ARMLOCAL_ARM_GenAmortization (LocalGetNumObjectId (C_swapleg),
												amortMethodId,
												freqId,
												C_amortAmount,
												dayCountId,
												C_legNotional,
												C_amortRate,
												C_reducedMaturity,
												C_modelId,
												C_cleanup,
												C_result);

		if(retCode == ARM_OK)
		{
			objId = C_result.getLong ();

			stringId = LocalMakeObjectId (objId, curClass);
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_PXL_ARM_GenAmortization" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


__declspec(dllexport) LPXLOPER WINAPI Local_SWAP_WITH_NOTIONNAL (LPXLOPER XL_swLeg1,
																 LPXLOPER XL_swLeg2,
																 LPXLOPER XL_not,
																 LPXLOPER XL_minPay)
{
	ADD_LOG("Local_SWAP_WITH_NOTIONNAL ");
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
		ARM_NOCALCIFWIZ();

		// C variable
		CCString C_swLeg1;
		CCString C_swLeg2;
		CCString C_not;

		double C_minPay;
		double C_minPay_default = -1;
		
		// error
		static int error;
		static char* reason = "";

		XL_readStrCell(XL_swLeg1,C_swLeg1," ARM_ERR: swap leg 1 id: object expected",C_result);
		XL_readStrCell(XL_swLeg2,C_swLeg2," ARM_ERR: swap leg 2 id: object expected",C_result);
		XL_readStrCell(XL_not,C_not," ARM_ERR: notionnal id: object expected",C_result);
		XL_readNumCellWD(XL_minPay,C_minPay,C_minPay_default," ARM_ERR: minimum pay: numeric expected",C_result);
		
		long retCode;
		long objId;
		CCString prevClass;
		
		CCString curClass = LOCAL_SWAP_CLASS;
		CCString stringId = GetLastCurCellEnvValue ();
		
		if(!stringId)
		{
			retCode = ARMLOCAL_SWAP_WITH_NOTIONNAL (LocalGetNumObjectId (C_swLeg1),
													LocalGetNumObjectId (C_swLeg2),
													LocalGetNumObjectId (C_not),
													(long)C_minPay,
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
				retCode = ARMLOCAL_SWAP_WITH_NOTIONNAL (LocalGetNumObjectId (C_swLeg1),
														LocalGetNumObjectId (C_swLeg2),
														LocalGetNumObjectId (C_not),
														(long)C_minPay,
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
				retCode = ARMLOCAL_SWAP_WITH_NOTIONNAL (LocalGetNumObjectId (C_swLeg1),
														LocalGetNumObjectId (C_swLeg2),
														LocalGetNumObjectId (C_not),
														(long)C_minPay,
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_SWAP_WITH_NOTIONNAL" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}



__declspec(dllexport) LPXLOPER WINAPI Local_PXL_SWAP_WITH_NOTIONNAL (LPXLOPER XL_swLeg1,
																	 LPXLOPER XL_swLeg2,
																	 LPXLOPER XL_not,
																	 LPXLOPER XL_minPay)
{
	ADD_LOG("Local_PXL_SWAP_WITH_NOTIONNAL ");
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
		ARM_NOCALCIFWIZ();

		// C variable
		CCString C_swLeg1;
		CCString C_swLeg2;
		CCString C_not;

		double C_minPay;
		double C_minPay_default = -1;
		
		// error
		static int error;
		static char* reason = "";

		XL_readStrCell(XL_swLeg1,C_swLeg1," ARM_ERR: swap leg 1 id: object expected",C_result);
		XL_readStrCell(XL_swLeg2,C_swLeg2," ARM_ERR: swap leg 2 id: object expected",C_result);
		XL_readStrCell(XL_not,C_not," ARM_ERR: notionnal id: object expected",C_result);
		XL_readNumCellWD(XL_minPay,C_minPay,C_minPay_default," ARM_ERR: minimum pay: numeric expected",C_result);
		
		long retCode;
		long objId;
		
		CCString curClass = LOCAL_SWAP_CLASS;
		CCString stringId;

		retCode = ARMLOCAL_SWAP_WITH_NOTIONNAL (LocalGetNumObjectId (C_swLeg1),
												LocalGetNumObjectId (C_swLeg2),
												LocalGetNumObjectId (C_not),
												(long)C_minPay,
												C_result);
									 
		if(retCode == ARM_OK)
		{
			objId = C_result.getLong ();

			stringId = LocalMakeObjectId (objId, curClass);
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_PXL_SWAP_WITH_NOTIONNAL" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


__declspec(dllexport) LPXLOPER WINAPI Local_ARM_CUSTOMFIRSTDATE (LPXLOPER XL_swlegId,
																 LPXLOPER XL_date)
{
	ADD_LOG("Local_ARM_CUSTOMFIRSTDATE ");
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
		ARM_NOCALCIFWIZ();

		// C variable
		CCString C_swlegId;

		double C_date;

		// error
		static int error;
		static char* reason = "";

		XL_readStrCell(XL_swlegId,C_swlegId," ARM_ERR: swapleg id: object expected",C_result);
		XL_readNumCell(XL_date,C_date," ARM_ERR: date: numeric expected",C_result);

		long retCode = ARMLOCAL_ARM_CUSTOMFSTCPN (LocalGetNumObjectId (C_swlegId),
												  C_date,
												  C_result);

		if ( retCode == ARM_OK )
		{
			FreeCurCellErr();
		
			XL_result.xltype = xltypeStr;
			XL_result.val.str = XL_StrC2StrPascal(C_swlegId);
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_ARM_CUSTOMFIRSTDATE" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}




__declspec(dllexport) LPXLOPER WINAPI Local_REVERSESTICKYLEG (LPXLOPER XL_startDate,
															  LPXLOPER XL_endDate,
															  LPXLOPER XL_liborType,
															  LPXLOPER XL_receiveOrPay,
															  LPXLOPER XL_spreads,
															  LPXLOPER XL_strikes,
															  LPXLOPER XL_first_coupon,
															  LPXLOPER XL_resetFreq,
															  LPXLOPER XL_payFreq,
															  LPXLOPER XL_resetTiming,
															  LPXLOPER XL_payTiming,
															  LPXLOPER XL_ccy)
{
	ADD_LOG("Local_REVERSESTICKYLEG ");
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

		CCString C_liborType;
		long liborTypeId;
		
		CCString C_receiveOrPay;
		long receiveOrPayId;

		double C_default_coupon = 0;
		
		CCString C_resetFreq;
		long resetFreqId;

		CCString C_payFreq;
		long payFreqId;

		CCString C_resetTiming;
		long resetTimingId;
		
		CCString C_payTiming;
		long payTimingId;
		
		CCString C_ccy;
		bool ccyIsObject = false;

		CCString C_kRefValId;
		CCString C_SpreadRefValId;

		double firstCoupon;
		double C_firstCoupon_default=0.0;

		// error
		static int error;
		static char* reason = "";

		XL_readNumCell(XL_startDate,C_startDate," ARM_ERR: start date: date expected",C_result);
		XL_readNumCell(XL_endDate,C_endDate," ARM_ERR: end date: date expected",C_result);
		XL_readStrCell(XL_liborType,C_liborType," ARM_ERR: libor type: string expected",C_result);
		XL_readStrCell(XL_receiveOrPay,C_receiveOrPay," ARM_ERR: receive or pay: string expected",C_result);
		XL_readStrCellWD(XL_resetFreq,C_resetFreq,"-1"," ARM_ERR: reset frequency: string expected",C_result);
		XL_readStrCellWD(XL_payFreq,C_payFreq,"-1"," ARM_ERR: pay frequency: string expected",C_result);
		XL_readStrCellWD(XL_resetTiming,C_resetTiming,"ADV"," ARM_ERR: reset timing: string expected",C_result);
		XL_readStrCellWD(XL_payTiming,C_payTiming,"ARR"," ARM_ERR: pay timing: string expected",C_result);
		XL_readStrCellWD(XL_ccy,C_ccy,"DEFAULT"," ARM_ERR: currency id: string expected",C_result);

		XL_readStrCell(XL_strikes,C_kRefValId," ARM_ERR: strike reference value id: object expected",C_result);
		XL_readStrCell(XL_spreads,C_SpreadRefValId," ARM_ERR: spread reference value id: object expected",C_result);
   		XL_readNumCellWD(XL_first_coupon,firstCoupon,C_firstCoupon_default," ARM_ERR: first coupon: numeric expected",C_result);


		if((liborTypeId = ARM_ConvIrIndName (C_liborType, C_result)) == ARM_DEFAULT_ERR)
		{
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}

		if((receiveOrPayId = ARM_ConvRecOrPay (C_receiveOrPay, C_result)) == ARM_DEFAULT_ERR)
		{
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}

		if((payFreqId = ARM_ConvFrequency (C_payFreq, C_result)) == ARM_DEFAULT_ERR)
		{
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}

		if((resetFreqId = ARM_ConvFrequency (C_resetFreq, C_result)) == ARM_DEFAULT_ERR)
		{
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}

		payTimingId = ARM_ConvPayResetRule (C_payTiming);

		resetTimingId = ARM_ConvPayResetRule (C_resetTiming);

		if (( C_ccy.GetLen() > 3 )
			&&
			( !(C_ccy == "DEFAULT") )
		   )
		   ccyIsObject = true;

		long retCode;
		long objId;
		CCString prevClass;

		CCString curClass = LOCAL_REVERSESTICKYLEG_CLASS;
		CCString stringId = GetLastCurCellEnvValue ();
		
		if(!stringId)
		{
			retCode = ARMLOCAL_REVERSESTICKYLEG (C_startDate,
												 C_endDate,
												 liborTypeId,
												 receiveOrPayId,									 
												 resetFreqId,
												 payFreqId,
												 resetTimingId,
												 payTimingId,
												 ccyIsObject,
												 C_ccy,
												 LocalGetNumObjectId (C_kRefValId),
												 LocalGetNumObjectId (C_SpreadRefValId),
												 firstCoupon,
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
				retCode = ARMLOCAL_REVERSESTICKYLEG (C_startDate,
													 C_endDate,
													 liborTypeId,
													 receiveOrPayId,									 
													 resetFreqId,
													 payFreqId,
													 resetTimingId,
													 payTimingId,
													 ccyIsObject,
													 C_ccy,
													 LocalGetNumObjectId (C_kRefValId),
													 LocalGetNumObjectId (C_SpreadRefValId),
													 firstCoupon,
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

				retCode = ARMLOCAL_REVERSESTICKYLEG (C_startDate,
													 C_endDate,
													 liborTypeId,
													 receiveOrPayId,									 
													 resetFreqId,
													 payFreqId,
													 resetTimingId,
													 payTimingId,
													 ccyIsObject,
													 C_ccy,
													 LocalGetNumObjectId (C_kRefValId),
													 LocalGetNumObjectId (C_SpreadRefValId),
													 firstCoupon,
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_REVERSESTICKYLEG" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}



__declspec(dllexport) LPXLOPER WINAPI Local_PXL_REVERSESTICKYLEG (LPXLOPER XL_startDate,
																  LPXLOPER XL_endDate,
																  LPXLOPER XL_liborType,
																  LPXLOPER XL_receiveOrPay,
																  LPXLOPER XL_spreads,
																  LPXLOPER XL_strikes,
																  LPXLOPER XL_first_coupon,
																  LPXLOPER XL_resetFreq,
																  LPXLOPER XL_payFreq,
																  LPXLOPER XL_resetTiming,
																  LPXLOPER XL_payTiming,
																  LPXLOPER XL_ccy)
{
	ADD_LOG("Local_PXL_REVERSESTICKYLEG ");
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

		CCString C_liborType;
		long liborTypeId;
		
		double firstCoupon;
		double C_firstCoupon_default=0.0;

		CCString C_receiveOrPay;
		long receiveOrPayId;

		double C_default_coupon = 0;
		
		CCString C_resetFreq;
		long resetFreqId;

		CCString C_payFreq;
		long payFreqId;

		CCString C_resetTiming;
		long resetTimingId;
		
		CCString C_payTiming;
		long payTimingId;
		
		CCString C_ccy;
		bool ccyIsObject = false;
				
		CCString C_kRefValId;
		CCString C_SpreadRefValId;

		// error
		static int error;
		static char* reason = "";

		XL_readNumCell(XL_startDate,C_startDate," ARM_ERR: start date: date expected",C_result);
		XL_readNumCell(XL_endDate,C_endDate," ARM_ERR: end date: date expected",C_result);
		XL_readStrCell(XL_liborType,C_liborType," ARM_ERR: libor type: string expected",C_result);
		XL_readStrCell(XL_receiveOrPay,C_receiveOrPay," ARM_ERR: receive or pay: string expected",C_result);
   		XL_readNumCellWD(XL_first_coupon,firstCoupon,C_firstCoupon_default," ARM_ERR: first coupon: numeric expected",C_result);
		XL_readStrCellWD(XL_resetFreq,C_resetFreq,"-1"," ARM_ERR: reset frequency: string expected",C_result);
		XL_readStrCellWD(XL_payFreq,C_payFreq,"-1"," ARM_ERR: pay frequency: string expected",C_result);
		XL_readStrCellWD(XL_resetTiming,C_resetTiming,"ADV"," ARM_ERR: reset timing: string expected",C_result);
		XL_readStrCellWD(XL_payTiming,C_payTiming,"ARR"," ARM_ERR: pay timing: string expected",C_result);
		XL_readStrCellWD(XL_ccy,C_ccy,"DEFAULT"," ARM_ERR: currency id: string expected",C_result);

		XL_readStrCell(XL_strikes,C_kRefValId," ARM_ERR: strike reference value id: object expected",C_result);
		XL_readStrCell(XL_spreads,C_SpreadRefValId," ARM_ERR: spread reference value id: object expected",C_result);

		if((liborTypeId = ARM_ConvIrIndName (C_liborType, C_result)) == ARM_DEFAULT_ERR)
		{
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}

		if((receiveOrPayId = ARM_ConvRecOrPay (C_receiveOrPay, C_result)) == ARM_DEFAULT_ERR)
		{
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}

		if((payFreqId = ARM_ConvFrequency (C_payFreq, C_result)) == ARM_DEFAULT_ERR)
		{
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}

		if((resetFreqId = ARM_ConvFrequency (C_resetFreq, C_result)) == ARM_DEFAULT_ERR)
		{
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}

		payTimingId = ARM_ConvPayResetRule(C_payTiming);

		resetTimingId = ARM_ConvPayResetRule(C_resetTiming);

		if (( C_ccy.GetLen() > 3 )
			&&
			( !(C_ccy == "DEFAULT") )
		   )
		   ccyIsObject = true;

		long retCode;
		long objId;
		
		CCString curClass = LOCAL_REVERSESTICKYLEG_CLASS;
		CCString stringId;
		
		retCode = ARMLOCAL_REVERSESTICKYLEG (C_startDate,
											 C_endDate,
											 liborTypeId,
											 receiveOrPayId,									 
											 resetFreqId,
											 payFreqId,
											 resetTimingId,
											 payTimingId,
											 ccyIsObject,
											 C_ccy,
											 LocalGetNumObjectId (C_kRefValId),
											 LocalGetNumObjectId (C_SpreadRefValId),
											 firstCoupon,
											 C_result);

		if(retCode == ARM_OK)
		{
			objId = C_result.getLong ();

			stringId = LocalMakeObjectId (objId, curClass);
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_PXL_REVERSESTICKYLEG" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
	
}



__declspec(dllexport) LPXLOPER WINAPI Local_IMPLIEDSPREADWITHMODELS (LPXLOPER XL_swap,
																	 LPXLOPER XL_model1,
																	 LPXLOPER XL_model2,
																	 LPXLOPER XL_price,
																	 LPXLOPER XL_leg1Or2)
{
	ADD_LOG("Local_IMPLIEDSPREADWITHMODELS ");
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
		ARM_NOCALCIFWIZ();

		// C variable
		CCString C_swap;
		CCString C_model1;
		CCString C_model2;
		double C_price;
		double C_leg1Or2;
		
		// error
		static int error;
		static char* reason = "";

		XL_readStrCell(XL_swap,C_swap," ARM_ERR: swap id: object expected",C_result);
		XL_readStrCell(XL_model1,C_model1," ARM_ERR: model 1 id: object expected",C_result);
		XL_readStrCell(XL_model2,C_model2," ARM_ERR: model 2 id: object expected",C_result);
		XL_readNumCell(XL_price,C_price," ARM_ERR: price: numeric expected",C_result);
		XL_readNumCell(XL_leg1Or2,C_leg1Or2," ARM_ERR: leg 1 or 2: numeric expected",C_result);
			
		long retCode = ARMLOCAL_IMPLIEDSPREADWITHMODELS (LocalGetNumObjectId (C_swap),
														 LocalGetNumObjectId (C_model1),
														 LocalGetNumObjectId (C_model2),
														 C_price,
														 (long)C_leg1Or2,
														 C_result);

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
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_IMPLIEDSPREADWITHMODELS" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}




__declspec(dllexport) LPXLOPER WINAPI Local_ARM_SetLastFixing(LPXLOPER XL_securityId,
															  LPXLOPER XL_rate,
															  LPXLOPER XL_FixingBeforeLast,
															  LPXLOPER XL_asof,
															  LPXLOPER XL_resetDate)
{
	ADD_LOG("Local_ARM_SetLastFixing");
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
		ARM_NOCALCIFWIZ();

		// C variable
		CCString C_securityId;
		
		double C_rate;
		double C_FixingBeforeLast;
		double C_FixingBeforeLastDef = 0.0;

		double C_asof;
		double C_resetDate = -1.0;
		double C_resetDate_default = -1.0;
		
		// error
		static int error;
		static char* reason = "";

		XL_readStrCell(XL_securityId,C_securityId," ARM_ERR: security id: object expected",C_result);
		XL_readNumCell(XL_rate,C_rate," ARM_ERR: fixing rate: numeric expected",C_result);
		XL_readNumCellWD(XL_FixingBeforeLast,C_FixingBeforeLast,C_FixingBeforeLastDef," ARM_ERR: fixing: rate expected",C_result);
		XL_readNumCell(XL_asof,C_asof," ARM_ERR: as of date: date expected",C_result);
		XL_readNumCellWD(XL_resetDate,C_resetDate,C_resetDate_default," ARM_ERR: reset date: date expected",C_result);

		long retCode;
		long objId;
		CCString prevClass;

		CCString curClass = LocalGetStringObjectClass(C_securityId);
		CCString stringId = GetLastCurCellEnvValue ();
		
		if(!stringId)
		{
			retCode = ARMLOCAL_ARM_SetLastFixing(LocalGetNumObjectId (C_securityId),
															C_rate,
															C_FixingBeforeLast,
															C_asof,
															C_resetDate,
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
				retCode = ARMLOCAL_ARM_SetLastFixing(LocalGetNumObjectId (C_securityId),
																C_rate,
																C_FixingBeforeLast,
																C_asof,
																C_resetDate,
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

				retCode = ARMLOCAL_ARM_SetLastFixing(LocalGetNumObjectId (C_securityId),
																C_rate,
																C_FixingBeforeLast,
																C_asof,
																C_resetDate,
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_ARM_SetLastFixing" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


__declspec(dllexport) LPXLOPER WINAPI Local_PXL_ARM_SetLastFixing(LPXLOPER XL_securityId,
																  LPXLOPER XL_rate,
																  LPXLOPER XL_FixingBeforeLast,
																  LPXLOPER XL_asof,
																  LPXLOPER XL_resetDate)
{
	ADD_LOG("Local_PXL_ARM_SetLastFixing");
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
		ARM_NOCALCIFWIZ();

		// C variable
		CCString C_securityId;
		
		double C_rate;
		double C_FixingBeforeLast;
		double C_FixingBeforeLastDef = 0.0;

		double C_asof;
		double C_resetDate;
		double C_resetDate_default = -1.0;
		
		// error
		static int error;
		static char* reason = "";

		XL_readStrCell(XL_securityId,C_securityId," ARM_ERR: security id: object expected",C_result);
		XL_readNumCell(XL_rate,C_rate," ARM_ERR: fixing rate: numeric expected",C_result);
		XL_readNumCellWD(XL_FixingBeforeLast,C_FixingBeforeLast,C_FixingBeforeLastDef," ARM_ERR: fixing: rate expected",C_result);
		XL_readNumCell(XL_asof,C_asof," ARM_ERR: as of date: date expected",C_result);
		XL_readNumCellWD(XL_resetDate,C_resetDate,C_resetDate_default," ARM_ERR: reset date: date expected",C_result);

		long retCode;
		long objId;

		CCString curClass = LocalGetStringObjectClass(C_securityId);
		CCString stringId;
		
		retCode = ARMLOCAL_ARM_SetLastFixing(LocalGetNumObjectId (C_securityId),
											 C_rate,
											 C_FixingBeforeLast,
											 C_asof,
											 C_resetDate,
											 C_result);

		if(retCode == ARM_OK)
		{
			objId = C_result.getLong ();

			stringId = LocalMakeObjectId (objId, curClass);
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_ARM_SetLastFixing" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


// Add a vector of past fixings to a swap leg
// Input : 
// - securityId : the security to consider
// - fixRates : vector of fixed rates
__declspec(dllexport) LPXLOPER WINAPI Local_ARM_SetFixRates(LPXLOPER XL_securityId,
															LPXLOPER XL_fixedRates)
{
	ADD_LOG("Local_ARM_SetFixRates");
	//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
		ARM_NOCALCIFWIZ();

		// C variable
		CCString C_securityId;
		
		VECTOR<double> C_fixedRates;
		VECTOR<double> C_fixedRates_default (0);
		
		// error
		static int error;
		static char* reason = "";

		XL_readStrCell(XL_securityId,C_securityId," ARM_ERR: security id: object expected",C_result);
		XL_readNumVectorWD(XL_fixedRates,C_fixedRates,C_fixedRates_default," ARM_ERR: past fixed rates: array of numeric expected",C_result);

		long retCode;
		long objId;
		CCString prevClass;

		CCString curClass = LocalGetStringObjectClass(C_securityId);
		CCString stringId = GetLastCurCellEnvValue ();
		
		if(!stringId)
		{
			retCode = ARMLOCAL_ARM_SetFixRates(LocalGetNumObjectId (C_securityId),
															C_fixedRates,
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
				retCode = ARMLOCAL_ARM_SetFixRates(LocalGetNumObjectId (C_securityId),
																C_fixedRates,
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

				retCode = ARMLOCAL_ARM_SetFixRates(LocalGetNumObjectId (C_securityId),
																C_fixedRates,
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_ARM_SetFixRates" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}



// Add a vector of past fixings to a swap leg
// Input : 
// - securityId : the security to consider
// - fixRates : vector of fixed rates
__declspec(dllexport) LPXLOPER WINAPI Local_PXL_ARM_SetFixRates(LPXLOPER XL_securityId,
																LPXLOPER XL_fixedRates)
{
	ADD_LOG("Local_PXL_ARM_SetFixRates");
	//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
		ARM_NOCALCIFWIZ();

		// C variable
		CCString C_securityId;

		VECTOR<double> C_fixedRates;
		VECTOR<double> C_fixedRates_default (0);
		
		// error
		static int error;
		static char* reason = "";

		XL_readStrCell(XL_securityId,C_securityId," ARM_ERR: security id: object expected",C_result);
		XL_readNumVectorWD(XL_fixedRates,C_fixedRates,C_fixedRates_default," ARM_ERR: past fixed rates: array of numeric expected",C_result);

		long retCode;
		long objId;
		CCString prevClass;

		CCString curClass = LocalGetStringObjectClass(C_securityId);
		CCString stringId = GetLastCurCellEnvValue ();
		
		retCode = ARMLOCAL_ARM_SetFixRates(LocalGetNumObjectId (C_securityId),
														C_fixedRates,
														C_result);

		if(retCode == ARM_OK)
		{
			objId = C_result.getLong ();

			stringId = LocalMakeObjectId (objId, curClass);
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_ARM_SetFixRates" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}



__declspec(dllexport) LPXLOPER WINAPI Local_ARM_GENLEG (LPXLOPER XL_startDates,
														LPXLOPER XL_endDates,
														LPXLOPER XL_paymentDates,
														LPXLOPER XL_resetDates,
														LPXLOPER XL_intDays,
														LPXLOPER XL_notional,
														LPXLOPER XL_irIndex,
														LPXLOPER XL_rcvOrPay,
														LPXLOPER XL_fixRateOrSpread,
														LPXLOPER XL_fwdOrFixing,
														LPXLOPER XL_ccy,
														LPXLOPER XL_nxChange,
														LPXLOPER XL_couponDayCount)
{
	ADD_LOG("Local_ARM_GENLEG ");
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
		ARM_NOCALCIFWIZ();

		// C variable
		VECTOR<double> C_startDates;
		VECTOR<double> C_endDates;
		VECTOR<double> C_paymentDates;
		VECTOR<double> C_resetDates;
		VECTOR<double> C_intDays;

		CCString C_notional;
		long notionalId;

		CCString C_irIndex;
		long irIndexId;

		CCString C_receiveOrPay;
		long receiveOrPayId;

		double C_fixrateorspread_double;
		double C_fixrateorspread_double_default = 0.0;
		CCString C_fixrateorspread_str;
		long   fixrateorspreadType;

		CCString C_ccy;
		bool ccyIsObject = false;

		VECTOR<double> C_fwdorfixing;
		VECTOR<double> C_fwdorfixing_default;

		CCString C_nxChange;
		long NxId;

		CCString C_couponDayCount;
		long dayCountId;

		// error
		static int error;
		static char* reason = "";

		XL_readNumVector(XL_startDates,C_startDates," ARM_ERR: start Dates: array of numeric expected",C_result);
		XL_readNumVector(XL_endDates,C_endDates," ARM_ERR: end Dates: array of numeric expected",C_result);
		XL_readNumVector(XL_paymentDates,C_paymentDates," ARM_ERR: payment Dates: array of numeric expected",C_result);
		XL_readNumVector(XL_resetDates,C_resetDates," ARM_ERR: reset Dates: array of numeric expected",C_result);
		XL_readNumVector(XL_intDays,C_intDays," ARM_ERR: interest days: array of numeric expected",C_result);
		XL_readStrCellWD(XL_notional,C_notional,"NULL"," ARM_ERR: notional: object expected (refvalue)",C_result);

		XL_readStrCellWD(XL_irIndex,C_irIndex,"NULL"," ARM_ERR: interest rate index id: object expected",C_result);
		XL_readStrCellWD(XL_rcvOrPay,C_receiveOrPay,"R"," ARM_ERR: receive or pay: string expected",C_result);

		XL_readStrOrNumCellWD(XL_fixRateOrSpread, C_fixrateorspread_str, C_fixrateorspread_double, C_fixrateorspread_double_default, fixrateorspreadType,
			   " ARM_ERR: spread (for floating leg) or fixed rate (for fixed leg): numeric or object ID string expected",C_result);
		XL_readNumVectorWD(XL_fwdOrFixing,C_fwdorfixing,C_fwdorfixing_default," ARM_ERR: fwd or fixing: array of numeric expected",C_result);

		XL_readStrCellWD(XL_ccy,C_ccy,"DEFAULT"," ARM_ERR: currency id: string expected",C_result);

		XL_readStrCellWD(XL_nxChange,C_nxChange,"NXNONE"," ARM_ERR: Notional xchge: string expected",C_result);
		XL_readStrCellWD(XL_couponDayCount,C_couponDayCount,"-1"," ARM_ERR: Coupon Day Count: string expected",C_result);


		if (( C_ccy.GetLen() > 3 )
			&&
			( !(C_ccy == "DEFAULT") )
		   )
		   ccyIsObject = true;

		if((receiveOrPayId = ARM_ConvRecOrPay (C_receiveOrPay, C_result)) == ARM_DEFAULT_ERR)
		{
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}

		if ( fixrateorspreadType == XL_TYPE_STRING )
		{
		   C_fixrateorspread_double = (double) LocalGetNumObjectId(C_fixrateorspread_str);

		   fixrateorspreadType = 1L;
		}
		else
		{
		   fixrateorspreadType = 0L;
		}

		NxId = ARM_NotionalExchange(C_nxChange);

		if (C_notional == "NULL")
			notionalId = ARM_NULL_OBJECT;
		else
			notionalId = LocalGetNumObjectId(C_notional);

		if (C_irIndex == "NULL")
			irIndexId = ARM_NULL_OBJECT;
		else
			irIndexId = LocalGetNumObjectId(C_irIndex);

		dayCountId = ARM_ConvDayCount(C_couponDayCount);

		long retCode;
		long objId;
		CCString prevClass;
		
		CCString curClass = LOCAL_SWAPLEG_CLASS;
		CCString stringId = GetLastCurCellEnvValue ();
		
		if(!stringId)
		{
			retCode = ARMLOCAL_GENLEG (C_startDates,
									   C_endDates,
									   C_paymentDates,
									   C_resetDates,
									   C_intDays,
									   C_fwdorfixing,
									   notionalId,
									   irIndexId,
									   (long)receiveOrPayId,
									   fixrateorspreadType,
									   C_fixrateorspread_double,
									   ccyIsObject,
									   C_ccy,
									   (long)NxId,
									   (int) dayCountId,
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
				retCode = ARMLOCAL_GENLEG (C_startDates,
										   C_endDates,
										   C_paymentDates,
										   C_resetDates,
										   C_intDays,
										   C_fwdorfixing,
										   notionalId,
										   irIndexId,
										   (long)receiveOrPayId,
										   fixrateorspreadType,
										   C_fixrateorspread_double,
										   ccyIsObject,
										   C_ccy,
										   (long)NxId,
										   (int) dayCountId,
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

				retCode = ARMLOCAL_GENLEG (C_startDates,
										   C_endDates,
										   C_paymentDates,
										   C_resetDates,
										   C_intDays,
										   C_fwdorfixing,
										   notionalId,
										   irIndexId,
										   (long)receiveOrPayId,
										   fixrateorspreadType,
										   C_fixrateorspread_double,
										   ccyIsObject,
										   C_ccy,
										   (long)NxId,
										   (int) dayCountId,
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_ARM_GENLEG" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


__declspec(dllexport) LPXLOPER WINAPI Local_PXL_ARM_GENLEG (LPXLOPER XL_startDates,
															LPXLOPER XL_endDates,
															LPXLOPER XL_paymentDates,
															LPXLOPER XL_resetDates,
															LPXLOPER XL_intDays,
															LPXLOPER XL_notional,
															LPXLOPER XL_irIndex,
															LPXLOPER XL_rcvOrPay,
															LPXLOPER XL_fixRateOrSpread,
															LPXLOPER XL_fwdOrFixing,
															LPXLOPER XL_ccy,
															LPXLOPER XL_nxChange,
															LPXLOPER XL_couponDayCount)
{
	ADD_LOG("Local_PXL_ARM_GENLEG ");
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
		ARM_NOCALCIFWIZ();

		// C variable
		VECTOR<double> C_startDates;
		VECTOR<double> C_endDates;
		VECTOR<double> C_paymentDates;
		VECTOR<double> C_resetDates;
		VECTOR<double> C_intDays;

		CCString C_notional;
		long notionalId;

		CCString C_irIndex;
		long irIndexId;

		CCString C_receiveOrPay;
		long receiveOrPayId;

		double C_fixrateorspread_double;
		double C_fixrateorspread_double_default = 0.0;
		CCString C_fixrateorspread_str;
		long   fixrateorspreadType;

		CCString C_ccy;
		bool ccyIsObject = false;

		VECTOR<double> C_fwdorfixing;
		VECTOR<double> C_fwdorfixing_default;

		CCString C_nxChange;
		long NxId;
		
		CCString C_couponDayCount;
		long dayCountId;

		// error
		static int error;
		static char* reason = "";

		XL_readNumVector(XL_startDates,C_startDates," ARM_ERR: start Dates: array of numeric expected",C_result);
		XL_readNumVector(XL_endDates,C_endDates," ARM_ERR: end Dates: array of numeric expected",C_result);
		XL_readNumVector(XL_paymentDates,C_paymentDates," ARM_ERR: payment Dates: array of numeric expected",C_result);
		XL_readNumVector(XL_resetDates,C_resetDates," ARM_ERR: reset Dates: array of numeric expected",C_result);
		XL_readNumVector(XL_intDays,C_intDays," ARM_ERR: interest days: array of numeric expected",C_result);
		XL_readStrCellWD(XL_notional,C_notional,"NULL"," ARM_ERR: notional: object expected (refvalue)",C_result);

		XL_readStrCellWD(XL_irIndex,C_irIndex,"NULL"," ARM_ERR: interest rate index id: object expected",C_result);
		XL_readStrCellWD(XL_rcvOrPay,C_receiveOrPay,"R"," ARM_ERR: receive or pay: string expected",C_result);

		XL_readStrOrNumCellWD(XL_fixRateOrSpread, C_fixrateorspread_str, C_fixrateorspread_double, C_fixrateorspread_double_default, fixrateorspreadType,
			   " ARM_ERR: spread (for floating leg) or fixed rate (for fixed leg): numeric or object ID string expected",C_result);
		XL_readNumVectorWD(XL_fwdOrFixing,C_fwdorfixing,C_fwdorfixing_default," ARM_ERR: fwd or fixing: array of numeric expected",C_result);

		XL_readStrCellWD(XL_ccy,C_ccy,"DEFAULT"," ARM_ERR: currency id: string expected",C_result);

		XL_readStrCellWD(XL_nxChange,C_nxChange,"NXNONE"," ARM_ERR: Notional xchge: string expected",C_result);
		XL_readStrCellWD(XL_couponDayCount,C_couponDayCount,"-1"," ARM_ERR: Coupon Day Count: string expected",C_result);


		if (( C_ccy.GetLen() > 3 )
			&&
			( !(C_ccy == "DEFAULT") )
		   )
		   ccyIsObject = true;

		if((receiveOrPayId = ARM_ConvRecOrPay (C_receiveOrPay, C_result)) == ARM_DEFAULT_ERR)
		{
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}

		if ( fixrateorspreadType == XL_TYPE_STRING )
		{
		   C_fixrateorspread_double = (double) LocalGetNumObjectId(C_fixrateorspread_str);

		   fixrateorspreadType = 1L;
		}
		else
		{
		   fixrateorspreadType = 0L;
		}

		NxId = ARM_NotionalExchange(C_nxChange);

		if (C_notional == "NULL")
			notionalId = ARM_NULL_OBJECT;
		else
			notionalId = LocalGetNumObjectId(C_notional);

		if (C_irIndex == "NULL")
			irIndexId = ARM_NULL_OBJECT;
		else
			irIndexId = LocalGetNumObjectId(C_irIndex);
		
		dayCountId = ARM_ConvDayCount(C_couponDayCount);

		long retCode;
		long objId;
		CCString prevClass;
		
		CCString curClass = LOCAL_SWAPLEG_CLASS;
		CCString stringId;
		
		retCode = ARMLOCAL_GENLEG (C_startDates,
								   C_endDates,
								   C_paymentDates,
								   C_resetDates,
								   C_intDays,
								   C_fwdorfixing,
								   notionalId,
								   irIndexId,
								   (long)receiveOrPayId,
								   fixrateorspreadType,
								   C_fixrateorspread_double,
								   ccyIsObject,
								   C_ccy,
								   (long)NxId,
								   (int) dayCountId,
								   C_result
								   );

		if(retCode == ARM_OK)
		{
			objId = C_result.getLong ();

			stringId = LocalMakeObjectId (objId, curClass);
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_PXL_ARM_GENLEG" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}



__declspec(dllexport) LPXLOPER WINAPI Local_SWAPLEG2 (LPXLOPER XL_irIndex,
													 LPXLOPER XL_startDate,
													 LPXLOPER XL_endDate,
													 LPXLOPER XL_receiveOrPay,
													 LPXLOPER XL_spread,
													 LPXLOPER XL_ccy,
													 LPXLOPER XL_dayCount,
													 LPXLOPER XL_resetTiming,
                                                     LPXLOPER XL_resetGap,
													 LPXLOPER XL_payTiming,
                                                     LPXLOPER XL_payGap,
													 LPXLOPER XL_decompPricingFlag)
{
	ADD_LOG("Local_SWAPLEG2 ");
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
		ARM_NOCALCIFWIZ();
		// C variable
		CCString C_irIndex;

		double C_startDate;
		double C_endDate;
		
		CCString C_receiveOrPay;
		long receiveOrPayId;

		double C_spread_double;
		CCString C_spread_str;
		long   spreadType;
		
		CCString C_ccy;
		bool ccyIsObject = false;

		CCString C_dayCount;
		long dayCountId;
    
		CCString C_resetTiming;
		long resetTimingId;
		
		double C_resetGap;
		double C_resetGap_default = 10000.;

		CCString C_payTiming;
		long payTimingId;	

		double C_payGap;
		double C_payGap_default = 10000.;

		double C_decompPricingFlag;
		double C_decompPricingFlag_default(1.0);

    
		// default swapleg values
		long     nxChange_default = ARM_NotionalExchange("NXNONE");
		long stubRuleId_default = ARM_ConvStubRule ("SS");
		double C_refDate_default = -1.0;
		long     adjStartDateId_default = ARM_ConvYesOrNo("YES", C_result);

		// error
		static int error;
		static char* reason = "";

		XL_readStrCell(XL_irIndex,C_irIndex," ARM_ERR: interest rate index id: object expected",C_result);
		XL_readNumCell(XL_startDate,C_startDate," ARM_ERR: start date: date expected",C_result);
		XL_readNumCell(XL_endDate,C_endDate," ARM_ERR: end date: date expected",C_result);
		XL_readStrCell(XL_receiveOrPay,C_receiveOrPay," ARM_ERR: receive or pay: string expected",C_result);
		XL_readStrOrNumCell(XL_spread, C_spread_str, C_spread_double, spreadType,
			   " ARM_ERR: spread: numeric or object ID string expected",C_result);
		XL_readStrCellWD(XL_ccy,C_ccy,"DEFAULT"," ARM_ERR: currency id: string expected",C_result);
		XL_readStrCellWD(XL_dayCount,C_dayCount,"-1"," ARM_ERR: day count: string expected",C_result);
		XL_readStrCellWD(XL_resetTiming,C_resetTiming,"ADV"," ARM_ERR: reset timing: string expected",C_result);
		XL_readNumCellWD(XL_resetGap,C_resetGap,C_resetGap_default," ARM_ERR: reset Gap: numeric expected",C_result);
		XL_readStrCellWD(XL_payTiming,C_payTiming,"ARR"," ARM_ERR: pay timing: string expected",C_result);
		XL_readNumCellWD(XL_payGap,C_payGap,C_payGap_default," ARM_ERR: pay Gap: numeric expected",C_result);
		XL_readNumCellWD(XL_decompPricingFlag,C_decompPricingFlag,C_decompPricingFlag_default," ARM_ERR: decompPricingFlag: numeric expected",C_result);
		
		if (( C_ccy.GetLen() > 3 )
			&&
			( !(C_ccy == "DEFAULT") )
		   )
		   ccyIsObject = true;

		if((receiveOrPayId = ARM_ConvRecOrPay (C_receiveOrPay, C_result)) == ARM_DEFAULT_ERR)
		{
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}

		if ( spreadType == XL_TYPE_STRING )
		{
		   C_spread_double = (double) LocalGetNumObjectId(C_spread_str);

		   spreadType = 1L;
		}
		else
		{
		   spreadType = 0L;
		}

		dayCountId = ARM_ConvDayCount (C_dayCount);

		payTimingId = ARM_ConvPayResetRule (C_payTiming);

		resetTimingId = ARM_ConvPayResetRule (C_resetTiming);

		long retCode;
		long objId;
		CCString prevClass;
		
		CCString curClass = LOCAL_SWAPLEG_CLASS;
		CCString stringId = GetLastCurCellEnvValue ();
		
		if(!stringId)
		{
			retCode = ARMLOCAL_SWAPLEG (LocalGetNumObjectId (C_irIndex),
										C_startDate,
										C_endDate,
										(long)receiveOrPayId,
										spreadType,
										C_spread_double,
										ccyIsObject,
										C_ccy,
										(long)dayCountId,
										resetTimingId,
										(long) C_resetGap,
										payTimingId,
										(long) C_payGap,
										"NULL",
										"NULL",
										(long) C_decompPricingFlag,
										nxChange_default,
										stubRuleId_default,
										C_refDate_default,
										adjStartDateId_default,
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
				retCode = ARMLOCAL_SWAPLEG (LocalGetNumObjectId (C_irIndex),
											C_startDate,
											C_endDate,
											(long)receiveOrPayId,
											spreadType,
											C_spread_double,
											ccyIsObject,
											C_ccy,
											(long)dayCountId,
											resetTimingId,
											(long) C_resetGap,
											payTimingId,
											(long) C_payGap,
											"NULL",
											"NULL",
											(long) C_decompPricingFlag,
											nxChange_default,
											stubRuleId_default,
											C_refDate_default,
											adjStartDateId_default,
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

				retCode = ARMLOCAL_SWAPLEG (LocalGetNumObjectId (C_irIndex),
											C_startDate,
											C_endDate,
											(long)receiveOrPayId,
											spreadType,
											C_spread_double,
											ccyIsObject,
											C_ccy,
											(long)dayCountId,
											resetTimingId,
											(long) C_resetGap,
											payTimingId,
											(long) C_payGap,
											"NULL",
											"NULL",
											(long) C_decompPricingFlag,
											nxChange_default,
											stubRuleId_default,
											C_refDate_default,
											adjStartDateId_default,
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_SWAPLEG2" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}

__declspec(dllexport) LPXLOPER WINAPI Local_PXL_SWAPLEG2 (LPXLOPER XL_irIndex,
													     LPXLOPER XL_startDate,
													     LPXLOPER XL_endDate,
													     LPXLOPER XL_receiveOrPay,
													     LPXLOPER XL_spread,
													     LPXLOPER XL_ccy,
													     LPXLOPER XL_dayCount,
													     LPXLOPER XL_resetTiming,
                                                         LPXLOPER XL_resetGap,
													     LPXLOPER XL_payTiming,
                                                         LPXLOPER XL_payGap,
													     LPXLOPER XL_decompPricingFlag)
{
	ADD_LOG("Local_PXL_SWAPLEG2 ");
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
		ARM_NOCALCIFWIZ();
		// C variable
		CCString C_irIndex;

		double C_startDate;
		double C_endDate;
		
		CCString C_receiveOrPay;
		long receiveOrPayId;

		double C_spread_double;
		CCString C_spread_str;
		long   spreadType;
		
		CCString C_ccy;
		bool ccyIsObject = false;

		CCString C_dayCount;
		long dayCountId;
    
		CCString C_resetTiming;
		long resetTimingId;
		
		double C_resetGap;
		double C_resetGap_default = 10000.;

		CCString C_payTiming;
		long payTimingId;	

		double C_payGap;
		double C_payGap_default = 10000.;

		double C_decompPricingFlag;
		double C_decompPricingFlag_default(1.0);

    
		// default swapleg values
		long     nxChange_default = ARM_NotionalExchange("NXNONE");
		long stubRuleId_default = ARM_ConvStubRule ("SS");
		double C_refDate_default = -1.0;
		long     adjStartDateId_default = ARM_ConvYesOrNo("YES", C_result);

		// error
		static int error;
		static char* reason = "";

		XL_readStrCell(XL_irIndex,C_irIndex," ARM_ERR: interest rate index id: object expected",C_result);
		XL_readNumCell(XL_startDate,C_startDate," ARM_ERR: start date: date expected",C_result);
		XL_readNumCell(XL_endDate,C_endDate," ARM_ERR: end date: date expected",C_result);
		XL_readStrCell(XL_receiveOrPay,C_receiveOrPay," ARM_ERR: receive or pay: string expected",C_result);
		XL_readStrOrNumCell(XL_spread, C_spread_str, C_spread_double, spreadType,
			   " ARM_ERR: spread: numeric or object ID string expected",C_result);
		XL_readStrCellWD(XL_ccy,C_ccy,"DEFAULT"," ARM_ERR: currency id: string expected",C_result);
		XL_readStrCellWD(XL_dayCount,C_dayCount,"-1"," ARM_ERR: day count: string expected",C_result);
		XL_readStrCellWD(XL_resetTiming,C_resetTiming,"ADV"," ARM_ERR: reset timing: string expected",C_result);
		XL_readNumCellWD(XL_resetGap,C_resetGap,C_resetGap_default," ARM_ERR: reset Gap: numeric expected",C_result);
		XL_readStrCellWD(XL_payTiming,C_payTiming,"ARR"," ARM_ERR: pay timing: string expected",C_result);
		XL_readNumCellWD(XL_payGap,C_payGap,C_payGap_default," ARM_ERR: pay Gap: numeric expected",C_result);
		XL_readNumCellWD(XL_decompPricingFlag,C_decompPricingFlag,C_decompPricingFlag_default," ARM_ERR: decompPricingFlag: numeric expected",C_result);
		
		if (( C_ccy.GetLen() > 3 )
			&&
			( !(C_ccy == "DEFAULT") )
		   )
		   ccyIsObject = true;

		if((receiveOrPayId = ARM_ConvRecOrPay (C_receiveOrPay, C_result)) == ARM_DEFAULT_ERR)
		{
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}

		if ( spreadType == XL_TYPE_STRING )
		{
		   C_spread_double = (double) LocalGetNumObjectId(C_spread_str);

		   spreadType = 1L;
		}
		else
		{
		   spreadType = 0L;
		}

		dayCountId = ARM_ConvDayCount (C_dayCount);

		payTimingId = ARM_ConvPayResetRule (C_payTiming);

		resetTimingId = ARM_ConvPayResetRule (C_resetTiming);

		long retCode;
		long objId;
		CCString prevClass;
		
		CCString curClass = LOCAL_SWAPLEG_CLASS;
		CCString stringId;
		
		retCode = ARMLOCAL_SWAPLEG (LocalGetNumObjectId (C_irIndex),
									C_startDate,
									C_endDate,
									(long)receiveOrPayId,
									spreadType,
									C_spread_double,
									ccyIsObject,
									C_ccy,
									(long)dayCountId,
									resetTimingId,
									(long) C_resetGap,
									payTimingId,
									(long) C_payGap,
									"NULL",
									"NULL",
									(long) C_decompPricingFlag,
									nxChange_default,
									stubRuleId_default,
									C_refDate_default,
									adjStartDateId_default,
									C_result);

		if(retCode == ARM_OK)
		{
			objId = C_result.getLong ();

			stringId = LocalMakeObjectId (objId, curClass);
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_PXL_SWAPLEG2" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}

__declspec(dllexport) LPXLOPER WINAPI Local_LIVRETALEG(LPXLOPER XL_startDate,
													  LPXLOPER XL_endDate,
													  LPXLOPER XL_receiveOrPay,
													  LPXLOPER XL_spread,
													  LPXLOPER XL_resetFreq,
													  LPXLOPER XL_payFreq,
													  LPXLOPER XL_resetTiming,
													  LPXLOPER XL_payTiming,
													  LPXLOPER XL_ccy,
													  LPXLOPER XL_intRule,
													  LPXLOPER XL_resetGap,
													  LPXLOPER XL_resetCal,
													  LPXLOPER XL_payCal,
													  LPXLOPER XL_decompPricingFlag,
													  LPXLOPER XL_nxChange,
													  LPXLOPER XL_stubRule,
													  LPXLOPER XL_refDate,
													  LPXLOPER XL_adjStartDate,
													  LPXLOPER XL_daycount)
{
	ADD_LOG("Local_LIVRETALEG");
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

		CCString C_liborType="LIVRETA";
		long liborTypeId;
		
		CCString C_receiveOrPay;
		long receiveOrPayId;

		CCString C_spread_str;
		double C_spread_double;
		double C_default_spread = 0;
		long   spreadType;

		CCString C_resetFreq;
		long resetFreqId;

		CCString C_payFreq;
		long payFreqId;

		CCString C_resetTiming;
		long resetTimingId;
		
		CCString C_payTiming;
		long payTimingId;
		
		CCString C_ccy;
		bool ccyIsObject = false;

		CCString C_intRule;
		long intRuleId;

		double C_resetGap;
		double C_resetGap_default = 10000.;

		CCString C_resetCal;
		CCString C_payCal;
		
		double C_decompPricingFlag;
		double C_decompPricingFlag_default(1.0);

		CCString C_nxChange;
		long     nxChange;

		CCString C_stubRule;
		long stubRuleId;

		double C_refDate;
		double C_refDate_default = -1.0;

		CCString C_adjStartDate;
		long     adjStartDateId;

		CCString C_daycount;
		CCString C_daycountDefault = "30/360";
		long daycountId;

		// error
		static int error;
		static char* reason = "";

		XL_readNumCell(XL_startDate,C_startDate," ARM_ERR: start date: date expected",C_result);
		XL_readNumCell(XL_endDate,C_endDate," ARM_ERR: end date: date expected",C_result);
		XL_readStrCell(XL_receiveOrPay,C_receiveOrPay," ARM_ERR: receive or pay: string expected",C_result);
		XL_readStrOrNumCellWD(XL_spread, C_spread_str, C_spread_double, C_default_spread, spreadType,
			   " ARM_ERR: spread: numeric or object ID string expected",C_result);
		XL_readStrCellWD(XL_resetFreq,C_resetFreq,"-1"," ARM_ERR: reset frequency: string expected",C_result);
		XL_readStrCellWD(XL_payFreq,C_payFreq,"-1"," ARM_ERR: pay frequency: string expected",C_result);
		XL_readStrCellWD(XL_resetTiming,C_resetTiming,"ADV"," ARM_ERR: reset timing: string expected",C_result);
		XL_readStrCellWD(XL_payTiming,C_payTiming,"ARR"," ARM_ERR: pay timing: string expected",C_result);
		XL_readStrCellWD(XL_ccy,C_ccy,"DEFAULT"," ARM_ERR: currency id: string expected",C_result);
		XL_readStrCellWD(XL_intRule,C_intRule,"ADJ"," ARM_ERR: intRule: string expected",C_result);
		XL_readNumCellWD(XL_resetGap,C_resetGap,C_resetGap_default," ARM_ERR: reset Gap: numeric expected",C_result);
		XL_readStrCellWD(XL_resetCal,C_resetCal,"NULL"," ARM_ERR: reset Calendar: string expected",C_result);
		XL_readStrCellWD(XL_payCal,C_payCal,"NULL"," ARM_ERR: pay Calendar: string expected",C_result);
		XL_readNumCellWD(XL_decompPricingFlag,C_decompPricingFlag,C_decompPricingFlag_default," ARM_ERR: decompPricingFlag: numeric expected",C_result);
		XL_readStrCellWD(XL_nxChange,C_nxChange,"NXNONE"," ARM_ERR: Notional xchge: string expected",C_result);
		XL_readStrCellWD(XL_stubRule,C_stubRule,"SS"," ARM_ERR: stub Rule: string expected",C_result);
		XL_readNumCellWD(XL_refDate,C_refDate,C_refDate_default," ARM_ERR: reference date: date expected",C_result);
		XL_readStrCellWD(XL_adjStartDate,C_adjStartDate,"YES"," ARM_ERR: adjStartDate: string expected",C_result);
		XL_readStrCellWD(XL_daycount, C_daycount, C_daycountDefault, "ARM_ERR: Daycount Basis: string expected",C_result);

		if((liborTypeId = ARM_ConvIrIndName (C_liborType, C_result)) == ARM_DEFAULT_ERR)
		{
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}

		if((receiveOrPayId = ARM_ConvRecOrPay (C_receiveOrPay, C_result)) == ARM_DEFAULT_ERR)
		{
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}

		if((payFreqId = ARM_ConvFrequency (C_payFreq, C_result)) == ARM_DEFAULT_ERR)
		{
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}

		if((resetFreqId = ARM_ConvFrequency (C_resetFreq, C_result)) == ARM_DEFAULT_ERR)
		{
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}

		payTimingId = ARM_ConvPayResetRule (C_payTiming);

		resetTimingId = ARM_ConvPayResetRule (C_resetTiming);

		if (( C_ccy.GetLen() > 3 )
			&&
			( !(C_ccy == "DEFAULT") )
		   )
		   ccyIsObject = true;


		intRuleId = ARM_ConvIntRule (C_intRule);

		nxChange = ARM_NotionalExchange(C_nxChange);

		stubRuleId = ARM_ConvStubRule (C_stubRule);
		if ( spreadType == XL_TYPE_STRING )
		{
		   C_spread_double = (double) LocalGetNumObjectId(C_spread_str);

		   spreadType = 1L;
		}
		else
		{
		   spreadType = 0L;
		}

		if ((adjStartDateId = ARM_ConvYesOrNo(C_adjStartDate, C_result)) == ARM_DEFAULT_ERR)
		{
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}

		daycountId = ARM_ConvDayCount(C_daycount);

		long retCode;
		long objId;
		CCString prevClass;
		
		CCString curClass = LOCAL_SWAPLEG_CLASS;
		CCString stringId = GetLastCurCellEnvValue ();
		
		if(!stringId)
		{
			retCode = ARMLOCAL_LIVRETALEG (C_startDate,
										 C_endDate,
										 liborTypeId,
										 receiveOrPayId,
										 spreadType,
										 C_spread_double,
										 resetFreqId,
										 payFreqId,
										 resetTimingId,
										 payTimingId,
										 ccyIsObject,
										 C_ccy,
										 intRuleId,
										 (long) C_resetGap,
										 C_resetCal,
										 C_payCal,
										 (long) C_decompPricingFlag,
										 nxChange,
										 stubRuleId,
										 C_refDate,
										 adjStartDateId,
										 daycountId,
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
				retCode = ARMLOCAL_LIVRETALEG (C_startDate,
											 C_endDate,
											 liborTypeId,
											 receiveOrPayId,
											 spreadType,
											 C_spread_double,
											 resetFreqId,
											 payFreqId,
											 resetTimingId,
											 payTimingId,
											 ccyIsObject,
											 C_ccy,
											 intRuleId,
											 (long) C_resetGap,
											 C_resetCal,
											 C_payCal,
											 (long) C_decompPricingFlag,
											 nxChange,
											 stubRuleId,
											 C_refDate,
											 adjStartDateId,
											 daycountId,
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

				retCode = ARMLOCAL_LIVRETALEG (C_startDate,
											 C_endDate,
											 liborTypeId,
											 receiveOrPayId,
											 spreadType,
											 C_spread_double,
											 resetFreqId,
											 payFreqId,
											 resetTimingId,
											 payTimingId,
											 ccyIsObject,
											 C_ccy,
											 intRuleId,
											 (long) C_resetGap,
											 C_resetCal,
											 C_payCal,
											 (long) C_decompPricingFlag,
											 nxChange,
											 stubRuleId,
											 C_refDate,
											 adjStartDateId,
											 daycountId,
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_LIVRETALEG" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


__declspec(dllexport) LPXLOPER WINAPI Local_PXL_LIVRETALEG(LPXLOPER XL_startDate,
														  LPXLOPER XL_endDate,
														  LPXLOPER XL_receiveOrPay,
														  LPXLOPER XL_spread,
														  LPXLOPER XL_resetFreq,
														  LPXLOPER XL_payFreq,
														  LPXLOPER XL_resetTiming,
														  LPXLOPER XL_payTiming,
														  LPXLOPER XL_ccy,
														  LPXLOPER XL_intRule,
														  LPXLOPER XL_resetGap,
														  LPXLOPER XL_resetCal,
														  LPXLOPER XL_payCal,
														  LPXLOPER XL_decompPricingFlag,
														  LPXLOPER XL_nxChange,
														  LPXLOPER XL_stubRule,
														  LPXLOPER XL_refDate,
														  LPXLOPER XL_adjStartDate,
														  LPXLOPER XL_daycount)
{
	ADD_LOG("Local_PXL_LIVRETALEG");
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

		CCString C_liborType="LIVRETA";
		long liborTypeId;
		
		CCString C_receiveOrPay;
		long receiveOrPayId;

		CCString C_spread_str;
		double C_spread_double;
		double C_default_spread = 0;
		long   spreadType;

		CCString C_resetFreq;
		long resetFreqId;

		CCString C_payFreq;
		long payFreqId;

		CCString C_resetTiming;
		long resetTimingId;
		
		CCString C_payTiming;
		long payTimingId;
		
		CCString C_ccy;
		bool ccyIsObject = false;

		CCString C_intRule;
		long intRuleId;

		double C_resetGap;
		double C_resetGap_default = 10000.;

		CCString C_resetCal;
		CCString C_payCal;
		
		double C_decompPricingFlag;
		double C_decompPricingFlag_default(1.0);

		CCString C_nxChange;
		long     nxChange;

		CCString C_stubRule;
		long stubRuleId;

		double C_refDate;
		double C_refDate_default = -1.0;

		CCString C_adjStartDate;
		long     adjStartDateId;

		CCString C_daycount;
		CCString C_daycountDefault = "30/360";
		long daycountId;

		// error
		static int error;
		static char* reason = "";

		XL_readNumCell(XL_startDate,C_startDate," ARM_ERR: start date: date expected",C_result);
		XL_readNumCell(XL_endDate,C_endDate," ARM_ERR: end date: date expected",C_result);
		XL_readStrCell(XL_receiveOrPay,C_receiveOrPay," ARM_ERR: receive or pay: string expected",C_result);
		XL_readStrOrNumCellWD(XL_spread, C_spread_str, C_spread_double, C_default_spread, spreadType,
			   " ARM_ERR: spread: numeric or object ID string expected",C_result);
		XL_readStrCellWD(XL_resetFreq,C_resetFreq,"-1"," ARM_ERR: reset frequency: string expected",C_result);
		XL_readStrCellWD(XL_payFreq,C_payFreq,"-1"," ARM_ERR: pay frequency: string expected",C_result);
		XL_readStrCellWD(XL_resetTiming,C_resetTiming,"ADV"," ARM_ERR: reset timing: string expected",C_result);
		XL_readStrCellWD(XL_payTiming,C_payTiming,"ARR"," ARM_ERR: pay timing: string expected",C_result);
		XL_readStrCellWD(XL_ccy,C_ccy,"DEFAULT"," ARM_ERR: currency id: string expected",C_result);
		XL_readStrCellWD(XL_intRule,C_intRule,"ADJ"," ARM_ERR: intRule: string expected",C_result);
		XL_readNumCellWD(XL_resetGap,C_resetGap,C_resetGap_default," ARM_ERR: reset Gap: numeric expected",C_result);
		XL_readStrCellWD(XL_resetCal,C_resetCal,"NULL"," ARM_ERR: reset Calendar: string expected",C_result);
		XL_readStrCellWD(XL_payCal,C_payCal,"NULL"," ARM_ERR: pay Calendar: string expected",C_result);
		XL_readNumCellWD(XL_decompPricingFlag,C_decompPricingFlag,C_decompPricingFlag_default," ARM_ERR: decompPricingFlag: numeric expected",C_result);
		XL_readStrCellWD(XL_nxChange,C_nxChange,"NXNONE"," ARM_ERR: Notional xchge: string expected",C_result);
		XL_readStrCellWD(XL_stubRule,C_stubRule,"SS"," ARM_ERR: stub Rule: string expected",C_result);
		XL_readNumCellWD(XL_refDate,C_refDate,C_refDate_default," ARM_ERR: reference date: date expected",C_result);
		XL_readStrCellWD(XL_adjStartDate,C_adjStartDate,"YES"," ARM_ERR: adjStartDate: string expected",C_result);
		XL_readStrCellWD(XL_daycount, C_daycount, C_daycountDefault, "ARM_ERR: Daycount Basis: string expected",C_result);

		if((liborTypeId = ARM_ConvIrIndName (C_liborType, C_result)) == ARM_DEFAULT_ERR)
		{
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}

		if((receiveOrPayId = ARM_ConvRecOrPay (C_receiveOrPay, C_result)) == ARM_DEFAULT_ERR)
		{
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}

		if((payFreqId = ARM_ConvFrequency (C_payFreq, C_result)) == ARM_DEFAULT_ERR)
		{
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}

		if((resetFreqId = ARM_ConvFrequency (C_resetFreq, C_result)) == ARM_DEFAULT_ERR)
		{
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}

		payTimingId = ARM_ConvPayResetRule (C_payTiming);

		resetTimingId = ARM_ConvPayResetRule (C_resetTiming);

		if (( C_ccy.GetLen() > 3 )
			&&
			( !(C_ccy == "DEFAULT") )
		   )
		   ccyIsObject = true;

		intRuleId = ARM_ConvIntRule (C_intRule);

		nxChange = ARM_NotionalExchange(C_nxChange);

		stubRuleId = ARM_ConvStubRule (C_stubRule);
		if ( spreadType == XL_TYPE_STRING )
		{
		   C_spread_double = (double) LocalGetNumObjectId(C_spread_str);

		   spreadType = 1L;
		}
		else
		{
		   spreadType = 0L;
		}

		if ((adjStartDateId = ARM_ConvYesOrNo(C_adjStartDate, C_result)) == ARM_DEFAULT_ERR)
		{
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}

		daycountId = ARM_ConvDayCount(C_daycount);

		long retCode;
		long objId;
		CCString prevClass;
		
		CCString curClass = LOCAL_SWAPLEG_CLASS;
		CCString stringId = GetLastCurCellEnvValue ();
		
		retCode = ARMLOCAL_LIVRETALEG (C_startDate,
									 C_endDate,
									 liborTypeId,
									 receiveOrPayId,
									 spreadType,
									 C_spread_double,
									 resetFreqId,
									 payFreqId,
									 resetTimingId,
									 payTimingId,
									 ccyIsObject,
									 C_ccy,
									 intRuleId,
									 (long) C_resetGap,
									 C_resetCal,
									 C_payCal,
									 (long) C_decompPricingFlag,
									 nxChange,
									 stubRuleId,
									 C_refDate,
									 adjStartDateId,
									 daycountId,
									 C_result);

		if(retCode == ARM_OK)
		{
			objId = C_result.getLong ();

			stringId = LocalMakeObjectId (objId, curClass);
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_LIVRETALEG" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////
						/* pour Livret A ( cap / floor hybride)  */


__declspec(dllexport) LPXLOPER WINAPI Local_AVERAGE_LIBOR_LEG(	LPXLOPER XL_irIndex,
																LPXLOPER XL_startDate,
																LPXLOPER XL_endDate,
																LPXLOPER XL_receiveOrPay,
																LPXLOPER XL_spread,
																LPXLOPER XL_ccy,
																LPXLOPER XL_dayCount,
																LPXLOPER XL_payCal,
																LPXLOPER XL_decompPricingFlag,
																LPXLOPER XL_nxChange,
																LPXLOPER XL_stubRule,
																LPXLOPER XL_refDate,
																LPXLOPER XL_adjStartDate,
																LPXLOPER XL_couru)
{
	ADD_LOG("Local_AVERAGE_LIBOR_LEG");
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
		ARM_NOCALCIFWIZ();

		// C variable
		CCString C_irIndex;

		double C_startDate;
		double C_endDate;
		
		CCString C_receiveOrPay;
		long receiveOrPayId;

		double C_spread_double;
		double C_couru;
		double C_couru_default(0.0);

		CCString C_spread_str;
		long   spreadType;
		
		CCString C_ccy;
		bool ccyIsObject = false;

		CCString C_dayCount;
		long dayCountId;

		CCString C_payCal;
		
		double C_decompPricingFlag;
		double C_decompPricingFlag_default(1.0);

		CCString C_nxChange;
		long     nxChange;

		CCString C_stubRule;
		long stubRuleId;

		double C_refDate;
		double C_refDate_default = -1.0;

		CCString C_adjStartDate;
		long     adjStartDateId;
;
		// error
		static int error;
		static char* reason = "";

		XL_readStrCell		(	XL_irIndex,				C_irIndex,			" ARM_ERR: interest rate index id: object expected",C_result);
		XL_readNumCell		(	XL_startDate,			C_startDate,		" ARM_ERR: start date: date expected",C_result);
		XL_readNumCell		(	XL_endDate,				C_endDate,			" ARM_ERR: end date: date expected",C_result);
		XL_readStrCell		(	XL_receiveOrPay,		C_receiveOrPay,		" ARM_ERR: receive or pay: string expected",C_result);
		XL_readStrOrNumCell	(	XL_spread,				C_spread_str,		C_spread_double, spreadType,  " ARM_ERR: spread: numeric or object ID string expected",C_result);
		XL_readStrCellWD	(	XL_ccy,					C_ccy,				"DEFAULT"," ARM_ERR: currency id: string expected",C_result);
		XL_readStrCellWD	(	XL_dayCount,			C_dayCount,			"-1"," ARM_ERR: day count: string expected",C_result);
		XL_readStrCellWD	(	XL_payCal,				C_payCal,			"NULL"," ARM_ERR: pay Calendar: string expected",C_result);
		XL_readNumCellWD	(	XL_decompPricingFlag,	C_decompPricingFlag,C_decompPricingFlag_default," ARM_ERR: decompPricingFlag: numeric expected",C_result);
		XL_readStrCellWD	(	XL_nxChange,			C_nxChange,			"NXNONE"," ARM_ERR: Notional xchge: string expected",C_result);
		XL_readStrCellWD	(	XL_stubRule,			C_stubRule,			"SS"," ARM_ERR: stub Rule: string expected",C_result);
		XL_readNumCellWD	(	XL_refDate,				C_refDate,			C_refDate_default," ARM_ERR: reference date: date expected",C_result);
		XL_readStrCellWD	(	XL_adjStartDate,		C_adjStartDate,		"YES"," ARM_ERR: adjStartDate: string expected",C_result);
		XL_readNumCellWD	(	XL_couru,				C_couru,			C_couru_default," ARM_ERR: start date: date expected",C_result);

		if (( C_ccy.GetLen() > 3 )
			&& 
			( !(C_ccy == "DEFAULT"))
		   )
		   ccyIsObject = true;

		if ((receiveOrPayId = ARM_ConvRecOrPay (C_receiveOrPay, C_result)) == ARM_DEFAULT_ERR)
		{
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}

		if ( spreadType == XL_TYPE_STRING )
		{
		   C_spread_double = (double) LocalGetNumObjectId(C_spread_str);

		   spreadType = 1L;
		}
		else
		{
		   spreadType = 0L;
		}

		dayCountId = ARM_ConvDayCount (C_dayCount);

		nxChange = ARM_NotionalExchange(C_nxChange);

		stubRuleId = ARM_ConvStubRule (C_stubRule);

		if ((adjStartDateId = ARM_ConvYesOrNo(C_adjStartDate, C_result)) == ARM_DEFAULT_ERR)
		{
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}

		long retCode;
		long objId;
		CCString prevClass;
		
		CCString curClass = LOCAL_SWAPLEG_CLASS;
		CCString stringId = GetLastCurCellEnvValue ();
		
		if (!stringId)
		{
			retCode = ARMLOCAL_AVERAGE_LIBOR_LEG (	LocalGetNumObjectId (C_irIndex),
													C_startDate,
													C_endDate,
													(long)	receiveOrPayId,
													spreadType,
													C_spread_double,
													ccyIsObject,
													C_ccy,
													(long)dayCountId,
													C_payCal,
													(long) C_decompPricingFlag,
													nxChange,
													stubRuleId,
													C_refDate,
													adjStartDateId,
													C_couru,
													C_result);

			if ( retCode == ARM_OK )
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
				retCode = ARMLOCAL_AVERAGE_LIBOR_LEG (	LocalGetNumObjectId (C_irIndex),
														C_startDate,
														C_endDate,
														(long) receiveOrPayId,
														spreadType,
														C_spread_double,
														ccyIsObject,
														C_ccy,
														(long) dayCountId,
														C_payCal,
														(long) C_decompPricingFlag,
														nxChange,
														stubRuleId,
														C_refDate,
														adjStartDateId,
														C_couru,
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
				
				retCode = ARMLOCAL_AVERAGE_LIBOR_LEG (	LocalGetNumObjectId (C_irIndex),
														C_startDate,
														C_endDate,
														(long) receiveOrPayId,
														spreadType,
														C_spread_double,
														ccyIsObject,
														C_ccy,
														(long) dayCountId,
														C_payCal,
														(long) C_decompPricingFlag,
														nxChange,
														stubRuleId,
														C_refDate,
														adjStartDateId,
														C_couru,
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in local_AVERAGE_LIBOR_LEG" )

	/// return the result as an LPXLOPER*/
	return (LPXLOPER)&XL_result;
}

__declspec(dllexport) LPXLOPER WINAPI Local_PXL_AVERAGE_LIBOR_LEG (	LPXLOPER XL_irIndex,
																	LPXLOPER XL_startDate,
																	LPXLOPER XL_endDate,
																	LPXLOPER XL_receiveOrPay,
																	LPXLOPER XL_spread,
																	LPXLOPER XL_ccy,
																	LPXLOPER XL_dayCount,
																	LPXLOPER XL_payCal,
																	LPXLOPER XL_decompPricingFlag,
																	LPXLOPER XL_nxChange,
																	LPXLOPER XL_stubRule,
																	LPXLOPER XL_refDate,
																	LPXLOPER XL_adjStartDate,
																	LPXLOPER XL_couru)
{
	ADD_LOG("Local_PXL_AVERAGE_LIBOR_LEG ");
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
		ARM_NOCALCIFWIZ();

		// C variable
		CCString C_irIndex;

		double C_startDate;
		double C_endDate;
		
		CCString C_receiveOrPay;
		long receiveOrPayId;

		double C_spread_double;
		double C_couru;
		double C_couru_default(0.0);

		CCString C_spread_str;
		long   spreadType;
		
		CCString C_ccy;
		bool ccyIsObject = false;

		CCString C_dayCount;
		long dayCountId;

		CCString C_payCal;
		
		double C_decompPricingFlag;
		double C_decompPricingFlag_default(1.0);

		CCString C_nxChange;
		long     nxChange;

		CCString C_stubRule;
		long stubRuleId;

		double C_refDate;
		double C_refDate_default = -1.0;

		CCString C_adjStartDate;
		long     adjStartDateId;

		// error
		static int error;
		static char* reason = "";

		XL_readStrCell		(	XL_irIndex,				C_irIndex,			" ARM_ERR: interest rate index id: object expected",C_result);
		XL_readNumCell		(	XL_startDate,			C_startDate,		" ARM_ERR: start date: date expected",C_result);
		XL_readNumCell		(	XL_endDate,				C_endDate,			" ARM_ERR: end date: date expected",C_result);
		XL_readStrCell		(	XL_receiveOrPay,		C_receiveOrPay,		" ARM_ERR: receive or pay: string expected",C_result);
		XL_readStrOrNumCell	(	XL_spread,				C_spread_str,		C_spread_double, spreadType,  " ARM_ERR: spread: numeric or object ID string expected",C_result);
		XL_readStrCellWD	(	XL_ccy,					C_ccy,				"DEFAULT"," ARM_ERR: currency id: string expected",C_result);
		XL_readStrCellWD	(	XL_dayCount,			C_dayCount,			"-1"," ARM_ERR: day count: string expected",C_result);
		XL_readStrCellWD	(	XL_payCal,				C_payCal,			"NULL"," ARM_ERR: pay Calendar: string expected",C_result);
		XL_readNumCellWD	(	XL_decompPricingFlag,	C_decompPricingFlag,C_decompPricingFlag_default," ARM_ERR: decompPricingFlag: numeric expected",C_result);
		XL_readStrCellWD	(	XL_nxChange,			C_nxChange,			"NXNONE"," ARM_ERR: Notional xchge: string expected",C_result);
		XL_readStrCellWD	(	XL_stubRule,			C_stubRule,			"SS"," ARM_ERR: stub Rule: string expected",C_result);
		XL_readNumCellWD	(	XL_refDate,				C_refDate,			C_refDate_default," ARM_ERR: reference date: date expected",C_result);
		XL_readStrCellWD	(	XL_adjStartDate,		C_adjStartDate,		"YES"," ARM_ERR: adjStartDate: string expected",C_result);
		XL_readNumCellWD	(	XL_couru,				C_couru,			C_couru_default," ARM_ERR: start date: date expected",C_result);

		if ((C_ccy.GetLen() > 3)
			&& !(C_ccy == "DEFAULT"))
			ccyIsObject = true;

		if((receiveOrPayId = ARM_ConvRecOrPay (C_receiveOrPay, C_result)) == ARM_DEFAULT_ERR)
		{
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}

		if ( spreadType == XL_TYPE_STRING )
		{
		   C_spread_double = (double) LocalGetNumObjectId(C_spread_str);

		   spreadType = 1L;
		}
		else
		{
		   spreadType = 0L;
		}

		dayCountId = ARM_ConvDayCount (C_dayCount);

		nxChange = ARM_NotionalExchange(C_nxChange);

		stubRuleId = ARM_ConvStubRule (C_stubRule);

		if ((adjStartDateId = ARM_ConvYesOrNo(C_adjStartDate, C_result)) == ARM_DEFAULT_ERR)
		{
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}

		long retCode;
		long objId;
		
		CCString curClass = LOCAL_SWAPLEG_CLASS;
		CCString stringId;


		retCode = ARMLOCAL_AVERAGE_LIBOR_LEG (	LocalGetNumObjectId	(	C_irIndex	),
												C_startDate,
												C_endDate,
												(long)	receiveOrPayId,
												spreadType,
												C_spread_double,
												ccyIsObject,
												C_ccy,
												(long)	dayCountId,
												C_payCal,
												(long)	C_decompPricingFlag,
												nxChange,
												stubRuleId,
												C_refDate,
												adjStartDateId,
												C_couru,
												C_result);

		if(retCode == ARM_OK)
		{
			objId = C_result.getLong ();

			stringId = LocalMakeObjectId (objId, curClass);
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

///	ARM_END();
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_PXL_AVERAGE_LIBOR_LEG" )

	/// return the result as an LPXLOPER*/
	return (LPXLOPER)&XL_result;
}

__declspec(dllexport) LPXLOPER WINAPI Local_UnderSwapFromBermuda(LPXLOPER XL_swaption,
																 LPXLOPER XL_NDIndex)
{
	ADD_LOG("Local_UnderSwapFromBermuda");
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
		ARM_NOCALCIFWIZ();

		// C variable
		CCString C_swaption;

		double   C_NDIndex;
		double	 C_NDIndex_default = -1.0;
		long	 NDIndex;

		// error
		static int error;
		static char* reason = "";

		XL_readStrCell	(XL_swaption, C_swaption, " ARM_ERR: swaption id: object expected", C_result);
		XL_readNumCellWD(XL_NDIndex, C_NDIndex, C_NDIndex_default, " ARM_ERR: Notice Date Index: integer expected", C_result);

		NDIndex	=	(long) C_NDIndex;

		long retCode;
		long objId;
		CCString prevClass;
		
		CCString curClass = LOCAL_SWAP_CLASS;
		CCString stringId = GetLastCurCellEnvValue ();
		
		if (!stringId)
		{
			retCode = ARMLOCAL_SWAP_FROM_SWAPTION (	LocalGetNumObjectId (C_swaption),
													NDIndex,
													C_result);

			if ( retCode == ARM_OK )
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
				
			if ( curClass == prevClass )
			{
				retCode = ARMLOCAL_SWAP_FROM_SWAPTION (	LocalGetNumObjectId (C_swaption),
														NDIndex,
														C_result,
														objId);

				if ( retCode == ARM_OK )
				{
					LocalSetCurCellEnvValue (curClass, objId); 

					stringId = LocalMakeObjectId (objId, curClass);
				}
			}
			else
			{
				FreeCurCellContent ();
				
				retCode = ARMLOCAL_SWAP_FROM_SWAPTION (	LocalGetNumObjectId (C_swaption),
														NDIndex,
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_SWAP_FROM_SWAPTION" )

	/// return the result as an LPXLOPER*/
	return (LPXLOPER)&XL_result;
}


__declspec(dllexport) LPXLOPER WINAPI Local_PXL_UnderSwapFromBermuda(LPXLOPER XL_swaption,
																	 LPXLOPER XL_NDIndex)
{
	ADD_LOG("Local_PXL_UnderSwapFromBermuda");
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
		ARM_NOCALCIFWIZ();

		// C variable
		CCString C_swaption;

		double   C_NDIndex;
		double	 C_NDIndex_default = -1.0;
		long	 NDIndex;

		// error
		static int error;
		static char* reason = "";

		XL_readStrCell	(XL_swaption, C_swaption, " ARM_ERR: swaption id: object expected", C_result);
		XL_readNumCellWD(XL_NDIndex, C_NDIndex, C_NDIndex_default, " ARM_ERR: Notice Date Index: integer expected", C_result);

		NDIndex	=	(long) C_NDIndex;

		long retCode;
		long objId;
		
		CCString curClass = LOCAL_SWAP_CLASS;
		CCString stringId;
		
		retCode = ARMLOCAL_SWAP_FROM_SWAPTION (	LocalGetNumObjectId (C_swaption),
												NDIndex,
												C_result);

		if (retCode == ARM_OK)
		{
			objId = C_result.getLong ();

			stringId = LocalMakeObjectId (objId, curClass);
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_PXL_SWAP_FROM_SWAPTION" )

	/// return the result as an LPXLOPER*/
	return (LPXLOPER)&XL_result;
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////


/*---- End Of File ----*/

// EOF %M% 
