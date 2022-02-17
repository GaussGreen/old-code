#pragma warning(disable : 4786)
#pragma warning(disable : 4005)

#include <libCCxll\CCxll.h>

#include <ARM\libarm_local\ARM_local_glob.h>
#include <ARM\libarm_local\ARM_local_class.h>
#include <ARM\libarm_local\ARM_local_swtion.h>
#include <ARM\libarm_local\ARM_local_swap.h>

#include "ARM_local_interface.h"
#include "ARM_local_interglob.h"
#include "ARM_xl_trycatch_local.h"

#include "util\tech_macro.h"

__declspec(dllexport) LPXLOPER WINAPI Local_SWAPTION (LPXLOPER XL_swapId,
													  LPXLOPER XL_receiveOrPay,
													  LPXLOPER XL_strike,
													  LPXLOPER XL_maturity,
													  LPXLOPER XL_exerciseType)
{
	ADD_LOG("Local_SWAPTION ");
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();


	// C variable
	CCString C_swapId;

	CCString C_receiveOrPay;
	long receiveOrPayId;

    CCString C_strike_str;
    double C_strike_double;
	long   strikeType;

	double C_maturity;
	
	CCString C_exerciseType;
	long exerciseTypeId;

	// error
	static int error;
	static char* reason = "";

	XL_readStrCell(XL_swapId,C_swapId," ARM_ERR: swap id: object expected",C_result);
	XL_readStrCell(XL_receiveOrPay,C_receiveOrPay," ARM_ERR: receive or pay: string expected",C_result);
	XL_readStrOrNumCell(XL_strike,C_strike_str,C_strike_double,strikeType," ARM_ERR: strike: numeric or object expected",C_result);
	XL_readNumCell(XL_maturity,C_maturity," ARM_ERR: maturity: date expected",C_result);
	XL_readStrCellWD(XL_exerciseType,C_exerciseType,"E"," ARM_ERR: exercise type: string expected",C_result);
	
	if((receiveOrPayId = ARM_ConvRecOrPay (C_receiveOrPay, C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	if((exerciseTypeId = ARM_ConvExerciseType (C_exerciseType, C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

    if ( strikeType == XL_TYPE_STRING )
	{
	   C_strike_double = (double) LocalGetNumObjectId(C_strike_str);

	   strikeType = 1L;
	}
	else
	{
	   strikeType = 0L;
	}

	long retCode;
	long objId;
	CCString prevClass;
	
	CCString curClass = LOCAL_SWAPTION_CLASS;
	CCString stringId = GetLastCurCellEnvValue ();

	if(!stringId)
	{
		retCode = ARMLOCAL_SWAPTION (LocalGetNumObjectId (C_swapId), receiveOrPayId,
									 strikeType, C_strike_double, C_maturity, exerciseTypeId, C_result);
							  
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
			retCode = ARMLOCAL_SWAPTION (LocalGetNumObjectId (C_swapId), receiveOrPayId,
										 strikeType, C_strike_double, C_maturity, exerciseTypeId, C_result, objId);
			if(retCode == ARM_OK)
			{			
				LocalSetCurCellEnvValue (curClass, objId); 

				stringId = LocalMakeObjectId (objId, curClass);
			}
		}
		else
		{
			FreeCurCellContent ();
			retCode = ARMLOCAL_SWAPTION (LocalGetNumObjectId (C_swapId), receiveOrPayId,
										 strikeType, C_strike_double, C_maturity, exerciseTypeId, C_result);
		
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_SWAPTION" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}



__declspec(dllexport) LPXLOPER WINAPI Local_PXL_SWAPTION (LPXLOPER XL_swapId,
														  LPXLOPER XL_receiveOrPay,
														  LPXLOPER XL_strike,
														  LPXLOPER XL_maturity,
														  LPXLOPER XL_exerciseType)
{
	ADD_LOG("Local_PXL_SWAPTION ");
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();


	// C variable
	CCString C_swapId;

	CCString C_receiveOrPay;
	long receiveOrPayId;

    CCString C_strike_str;
    double C_strike_double;
	long   strikeType;

	double C_maturity;
	
	CCString C_exerciseType;
	long exerciseTypeId;

	// error
	static int error;
	static char* reason = "";

	XL_readStrCell(XL_swapId,C_swapId," ARM_ERR: swap id: object expected",C_result);
	XL_readStrCell(XL_receiveOrPay,C_receiveOrPay," ARM_ERR: receive or pay: string expected",C_result);
	XL_readStrOrNumCell(XL_strike,C_strike_str,C_strike_double,strikeType," ARM_ERR: strike: numeric or object expected",C_result);
	XL_readNumCell(XL_maturity,C_maturity," ARM_ERR: maturity: date expected",C_result);
	XL_readStrCellWD(XL_exerciseType,C_exerciseType,"E"," ARM_ERR: exercise type: string expected",C_result);
	
	if((receiveOrPayId = ARM_ConvRecOrPay (C_receiveOrPay, C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	if((exerciseTypeId = ARM_ConvExerciseType (C_exerciseType, C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

    if ( strikeType == XL_TYPE_STRING )
	{
	   C_strike_double = (double) LocalGetNumObjectId(C_strike_str);

	   strikeType = 1L;
	}
	else
	{
	   strikeType = 0L;
	}

	long retCode;
	long objId;

	CCString curClass = LOCAL_SWAPTION_CLASS;
	CCString stringId;

	retCode = ARMLOCAL_SWAPTION (LocalGetNumObjectId (C_swapId), receiveOrPayId,
		                    strikeType, C_strike_double, C_maturity, exerciseTypeId, C_result);
							  
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_PXL_SWAPTION" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


__declspec(dllexport) LPXLOPER WINAPI Local_SwaptionFromExpiry(LPXLOPER XL_optionExpiry,
													           LPXLOPER XL_swapTerm,
													           LPXLOPER XL_liborType,
													           LPXLOPER XL_strike,
													           LPXLOPER XL_receiveOrPay,
													           LPXLOPER XL_spread,
													           LPXLOPER XL_ccy)
{
	ADD_LOG("Local_SwaptionFromExpiry");

//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();


	// C variable
	CCString C_optionExpiry;
	CCString C_swapTerm;

	CCString C_receiveOrPay;
	long receiveOrPayId;
	
	CCString C_liborType;
	long liborTypeId;

	double C_strike;

	CCString C_spread_str;
	double C_spread;
	double C_spread_default = 0.0;
	long spreadType;

	CCString C_ccy;
	bool ccyIsObject = false;
	
	// error
	static int error;
	static char* reason = "";
	
	XL_readStrCell(XL_optionExpiry,C_optionExpiry," ARM_ERR: maturity: date expected",C_result);
	XL_readStrCell(XL_swapTerm,C_swapTerm," ARM_ERR: start date: date expected",C_result);
	XL_readStrCell(XL_liborType,C_liborType," ARM_ERR: libor type: string expected",C_result);
	XL_readStrCell(XL_receiveOrPay,C_receiveOrPay," ARM_ERR: receive or pay: string expected",C_result);
	XL_readStrOrNumCellWD(XL_spread, C_spread_str, C_spread, C_spread_default, spreadType,
			   " ARM_ERR: spread: numeric or object ID string expected",C_result);
	XL_readNumCell(XL_strike,C_strike," ARM_ERR: strike: numeric expected",C_result);
	XL_readStrCellWD(XL_ccy,C_ccy,"DEFAULT"," ARM_ERR: currency: string expected",C_result);

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

	if ((C_ccy.GetLen() > 3)
		&&
		!(C_ccy == "DEFAULT"))
		ccyIsObject = true;

	if ( spreadType == XL_TYPE_STRING )
		{
		   C_spread = (double) LocalGetNumObjectId(C_spread_str);

		   spreadType = 1L;
		}
	else
	{
	   spreadType = 0L;
	}

	long retCode;
	long objId;
	CCString prevClass;
	
	CCString curClass = LOCAL_SWAPTION_CLASS;
	CCString stringId = GetLastCurCellEnvValue ();

	if (!stringId)
	{
		retCode = ARMLOCAL_SwaptionFromExpiry(C_optionExpiry, C_swapTerm, liborTypeId, 
										      receiveOrPayId, C_strike, spreadType,
										      C_spread, ccyIsObject, C_ccy, C_result);

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
			retCode = ARMLOCAL_SwaptionFromExpiry(C_optionExpiry, C_swapTerm, liborTypeId, 
										          receiveOrPayId, C_strike, spreadType,
										          C_spread, ccyIsObject, C_ccy, 
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

			retCode = ARMLOCAL_SwaptionFromExpiry (C_optionExpiry, C_swapTerm, liborTypeId, 
										           receiveOrPayId, C_strike, spreadType,
										           C_spread, ccyIsObject, C_ccy, C_result);
		
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_SWAPTIONFROMEXPIRYDATE" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


__declspec(dllexport) LPXLOPER WINAPI Local_PXL_SwaptionFromExpiry(LPXLOPER XL_optionExpiry,
													           LPXLOPER XL_swapTerm,
													           LPXLOPER XL_liborType,
													           LPXLOPER XL_strike,
													           LPXLOPER XL_receiveOrPay,
													           LPXLOPER XL_spread,
													           LPXLOPER XL_ccy)
{
	ADD_LOG("Local_PXL_SwaptionFromExpiry");

//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();


	// C variable
	CCString C_optionExpiry;
	CCString C_swapTerm;

	CCString C_receiveOrPay;
	long receiveOrPayId;
	
	CCString C_liborType;
	long liborTypeId;

	double C_strike;

	CCString C_spread_str;
	double C_spread;
	double C_spread_default = 0.0;
	long spreadType;

	CCString C_ccy;
	bool ccyIsObject = false;
	
	// error
	static int error;
	static char* reason = "";
	
	XL_readStrCell(XL_optionExpiry,C_optionExpiry," ARM_ERR: maturity: date expected",C_result);
	XL_readStrCell(XL_swapTerm,C_swapTerm," ARM_ERR: start date: date expected",C_result);
	XL_readStrCell(XL_liborType,C_liborType," ARM_ERR: libor type: string expected",C_result);
	XL_readStrCell(XL_receiveOrPay,C_receiveOrPay," ARM_ERR: receive or pay: string expected",C_result);
	XL_readStrOrNumCellWD(XL_spread, C_spread_str, C_spread, C_spread_default, spreadType,
			   " ARM_ERR: spread: numeric or object ID string expected",C_result);
	XL_readNumCell(XL_strike,C_strike," ARM_ERR: strike: numeric expected",C_result);
	XL_readStrCellWD(XL_ccy,C_ccy,"DEFAULT"," ARM_ERR: currency: string expected",C_result);

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

	if ((C_ccy.GetLen() > 3)
		&&
		!(C_ccy == "DEFAULT"))
		ccyIsObject = true;

	if ( spreadType == XL_TYPE_STRING )
		{
		   C_spread = (double) LocalGetNumObjectId(C_spread_str);

		   spreadType = 1L;
		}
		else
		{
		   spreadType = 0L;
		}

	long retCode;
	long objId;
	
	retCode = ARMLOCAL_SwaptionFromExpiry(C_optionExpiry, C_swapTerm, liborTypeId, 
										      receiveOrPayId, C_strike, spreadType,
										      C_spread, ccyIsObject, C_ccy, C_result);
	
	CCString stringId;
	CCString curClass = LOCAL_SWAPTION_CLASS;

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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_SWAPTIONFROMEXPIRYDATE" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


__declspec(dllexport) LPXLOPER WINAPI Local_LIBORSWAPTION (LPXLOPER XL_startDate,
														   LPXLOPER XL_endDate,
														   LPXLOPER XL_receiveOrPay,
														   LPXLOPER XL_strike,
														   LPXLOPER XL_maturity,
														   LPXLOPER XL_liborType,
														   LPXLOPER XL_spread,
														   LPXLOPER XL_exerciseType,
														   LPXLOPER XL_resetFreq,
														   LPXLOPER XL_payFreq,
														   LPXLOPER XL_ccy)
{
	ADD_LOG("Local_LIBORSWAPTION ");
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
	long receiveOrPayId;
	
	CCString C_liborType;
	long liborTypeId;

	double C_strike;
	double C_maturity;

	CCString C_spread_str;
	double C_spread;
	double C_spread_default = 0.0;
	long spreadType;

	CCString C_exerciseType;
	long exerciseTypeId;

	CCString C_resetFreq;
	long resetFreqId;

	CCString C_payFreq;
	long payFreqId;

	CCString C_ccy;
	bool ccyIsObject = false;
	
	// error
	static int error;
	static char* reason = "";

	XL_readNumCell(XL_startDate,C_startDate," ARM_ERR: start date: date expected",C_result);
	XL_readNumCell(XL_endDate,C_endDate," ARM_ERR: end date: date expected",C_result);
	XL_readStrCell(XL_receiveOrPay,C_receiveOrPay," ARM_ERR: receive or pay: string expected",C_result);
	XL_readStrCell(XL_liborType,C_liborType," ARM_ERR: libor type: string expected",C_result);
	XL_readNumCell(XL_strike,C_strike," ARM_ERR: strike: numeric expected",C_result);
	XL_readNumCell(XL_maturity,C_maturity," ARM_ERR: maturity: date expected",C_result);
	XL_readStrOrNumCellWD(XL_spread, C_spread_str, C_spread, C_spread_default, spreadType,
			   " ARM_ERR: spread: numeric or object ID string expected",C_result);
//	XL_readNumCellWD(XL_spread,C_spread,C_spread_default," ARM_ERR: spread: numeric expected",C_result);
	XL_readStrCellWD(XL_exerciseType,C_exerciseType,"E"," ARM_ERR: exercise type: string expected",C_result);
	XL_readStrCellWD(XL_resetFreq,C_resetFreq,"-1"," ARM_ERR: reset frequency: string expected",C_result);
	XL_readStrCellWD(XL_payFreq,C_payFreq,"-1"," ARM_ERR: pay frequency: string expected",C_result);
	XL_readStrCellWD(XL_ccy,C_ccy,"DEFAULT"," ARM_ERR: currency: string expected",C_result);

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

	if((exerciseTypeId = ARM_ConvExerciseType (C_exerciseType, C_result)) == ARM_DEFAULT_ERR)
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

	if ((C_ccy.GetLen() > 3)
		&&
		!(C_ccy == "DEFAULT"))
		ccyIsObject = true;

	if ( spreadType == XL_TYPE_STRING )
		{
		   C_spread = (double) LocalGetNumObjectId(C_spread_str);

		   spreadType = 1L;
		}
		else
		{
		   spreadType = 0L;
		}

	long retCode;
	long objId;
	CCString prevClass;
	
	CCString curClass = LOCAL_SWAPTION_CLASS;
	CCString stringId = GetLastCurCellEnvValue ();

	if(!stringId)
	{
		retCode = ARMLOCAL_LIBORSWAPTION (C_startDate, C_endDate, receiveOrPayId,
										 C_strike, C_maturity, liborTypeId,spreadType,
										 C_spread, exerciseTypeId, resetFreqId,
										 payFreqId, ccyIsObject, C_ccy, C_result);

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
				retCode = ARMLOCAL_LIBORSWAPTION (C_startDate, C_endDate, receiveOrPayId,
											 C_strike, C_maturity, liborTypeId,spreadType,
											 C_spread, exerciseTypeId, resetFreqId,
											 payFreqId, ccyIsObject, C_ccy, C_result, objId);

			if(retCode == ARM_OK)
			{			
				LocalSetCurCellEnvValue (curClass, objId); 

				stringId = LocalMakeObjectId (objId, curClass);
			}
		}
		else
		{
			FreeCurCellContent ();
			retCode = ARMLOCAL_LIBORSWAPTION (C_startDate, C_endDate, receiveOrPayId,
											 C_strike, C_maturity, liborTypeId,spreadType,
											 C_spread, exerciseTypeId, resetFreqId,
											 payFreqId, ccyIsObject, C_ccy, C_result);
		
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
//		XL_result.xltype = xltypeNum;
//		XL_result.val.num = C_result.getDouble ();
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_LIBORSWAPTION" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


__declspec(dllexport) LPXLOPER WINAPI Local_PXL_LIBORSWAPTION (LPXLOPER XL_startDate,
															   LPXLOPER XL_endDate,
															   LPXLOPER XL_receiveOrPay,
															   LPXLOPER XL_strike,
															   LPXLOPER XL_maturity,
															   LPXLOPER XL_liborType,
															   LPXLOPER XL_spread,
															   LPXLOPER XL_exerciseType,
															   LPXLOPER XL_resetFreq,
															   LPXLOPER XL_payFreq,
															   LPXLOPER XL_ccy)
{
	ADD_LOG("Local_PXL_LIBORSWAPTION ");
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
	long receiveOrPayId;
	
	CCString C_liborType;
	long liborTypeId;

	double C_strike;
	double C_maturity;

	CCString C_spread_str;
	double C_spread;
	double C_spread_default = 0.0;
	long spreadType;

	CCString C_exerciseType;
	long exerciseTypeId;

	CCString C_resetFreq;
	long resetFreqId;

	CCString C_payFreq;
	long payFreqId;

	CCString C_ccy;
	bool ccyIsObject = false;
	
	// error
	static int error;
	static char* reason = "";

	XL_readNumCell(XL_startDate,C_startDate," ARM_ERR: start date: date expected",C_result);
	XL_readNumCell(XL_endDate,C_endDate," ARM_ERR: end date: date expected",C_result);
	XL_readStrCell(XL_receiveOrPay,C_receiveOrPay," ARM_ERR: receive or pay: string expected",C_result);
	XL_readStrCell(XL_liborType,C_liborType," ARM_ERR: libor type: string expected",C_result);
	XL_readNumCell(XL_strike,C_strike," ARM_ERR: strike: numeric expected",C_result);
	XL_readNumCell(XL_maturity,C_maturity," ARM_ERR: maturity: date expected",C_result);
	XL_readStrOrNumCellWD(XL_spread, C_spread_str, C_spread, C_spread_default, spreadType,
			   " ARM_ERR: spread: numeric or object ID string expected",C_result);
//	XL_readNumCellWD(XL_spread,C_spread,C_spread_default," ARM_ERR: spread: numeric expected",C_result);
	XL_readStrCellWD(XL_exerciseType,C_exerciseType,"E"," ARM_ERR: exercise type: string expected",C_result);
	XL_readStrCellWD(XL_resetFreq,C_resetFreq,"-1"," ARM_ERR: reset frequency: string expected",C_result);
	XL_readStrCellWD(XL_payFreq,C_payFreq,"-1"," ARM_ERR: pay frequency: string expected",C_result);
	XL_readStrCellWD(XL_ccy,C_ccy,"DEFAULT"," ARM_ERR: currency: string expected",C_result);

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

	if((exerciseTypeId = ARM_ConvExerciseType (C_exerciseType, C_result)) == ARM_DEFAULT_ERR)
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

	if ((C_ccy.GetLen() > 3)
		&&
		!(C_ccy == "DEFAULT"))
		ccyIsObject = true;

	if ( spreadType == XL_TYPE_STRING )
	{
	   C_spread = (double) LocalGetNumObjectId(C_spread_str);

	   spreadType = 1L;
	}
	else
	{
	   spreadType = 0L;
	}

	long retCode;
	long objId;
	CCString prevClass;
	
	CCString curClass = LOCAL_SWAPTION_CLASS;
	CCString stringId;

	retCode = ARMLOCAL_LIBORSWAPTION (C_startDate, C_endDate, receiveOrPayId,
										 C_strike, C_maturity, liborTypeId,spreadType,
										 C_spread, exerciseTypeId, resetFreqId,
										 payFreqId, ccyIsObject, C_ccy, C_result);

	if(retCode == ARM_OK)
	{
		objId = C_result.getLong ();

		stringId = LocalMakeObjectId (objId, curClass);
	}

	if(retCode == ARM_OK)
	{
		FreeCurCellErr ();
		XL_result.xltype = xltypeStr;
//		XL_result.xltype = xltypeNum;
//		XL_result.val.num = C_result.getDouble ();
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_PXL_LIBORSWAPTION" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


__declspec(dllexport) LPXLOPER WINAPI Local_EXOSWAPTION (LPXLOPER XL_swapId,
														 LPXLOPER XL_receiveOrPay,
														 LPXLOPER XL_xStyleId,
														 LPXLOPER XL_kRefValId,
														 LPXLOPER XL_swapYearTerm)
{
	ADD_LOG("Local_EXOSWAPTION ");
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();


	// C variable
	CCString C_swapId;

	CCString C_receiveOrPay;
	long receiveOrPayId;

	CCString C_xStyleId;

	CCString C_kRefValId;

	double C_swapYearTerm;
	
	// error
	static int error;
	static char* reason = "";

	XL_readStrCell(XL_swapId,C_swapId," ARM_ERR: swap id: object expected",C_result);
	XL_readStrCell(XL_receiveOrPay,C_receiveOrPay," ARM_ERR: receive or pay: string expected",C_result);
	XL_readStrCell(XL_xStyleId,C_xStyleId," ARM_ERR: exercise style id: object expected",C_result);
	XL_readStrCell(XL_kRefValId,C_kRefValId," ARM_ERR: strike reference value id: object expected",C_result);
	XL_readNumCell(XL_swapYearTerm,C_swapYearTerm," ARM_ERR: swap year term: date expected",C_result);
	
	if((receiveOrPayId = ARM_ConvRecOrPay (C_receiveOrPay, C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	long retCode;
	long objId;
	CCString prevClass;
	
	CCString curClass = LOCAL_SWAPTION_CLASS;
	CCString stringId = GetLastCurCellEnvValue ();

	if(!stringId)
	{
		retCode = ARMLOCAL_EXOSWAPTION (LocalGetNumObjectId (C_swapId),
								   receiveOrPayId,
								   LocalGetNumObjectId (C_xStyleId),
								   LocalGetNumObjectId (C_kRefValId),
								   C_swapYearTerm,
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
			retCode = ARMLOCAL_EXOSWAPTION (LocalGetNumObjectId (C_swapId),
								       receiveOrPayId,
								       LocalGetNumObjectId (C_xStyleId),
								       LocalGetNumObjectId (C_kRefValId),
								       C_swapYearTerm,
								       C_result, objId);

			if(retCode == ARM_OK)
			{
				LocalSetCurCellEnvValue (curClass, objId); 

				stringId = LocalMakeObjectId (objId, curClass);
			}
		}
		else
		{
			FreeCurCellContent ();
			retCode = ARMLOCAL_EXOSWAPTION (LocalGetNumObjectId (C_swapId),
								       receiveOrPayId,
									   LocalGetNumObjectId (C_xStyleId),
									   LocalGetNumObjectId (C_kRefValId),
								       C_swapYearTerm,
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_EXOSWAPTION" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


__declspec(dllexport) LPXLOPER WINAPI Local_PXL_EXOSWAPTION (LPXLOPER XL_swapId,
															 LPXLOPER XL_receiveOrPay,
															 LPXLOPER XL_xStyleId,
															 LPXLOPER XL_kRefValId,
															 LPXLOPER XL_swapYearTerm)
{
	ADD_LOG("Local_PXL_EXOSWAPTION ");
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();


	// C variable
	CCString C_swapId;

	CCString C_receiveOrPay;
	long receiveOrPayId;

	CCString C_xStyleId;

	CCString C_kRefValId;

	double C_swapYearTerm;
	
	// error
	static int error;
	static char* reason = "";

	XL_readStrCell(XL_swapId,C_swapId," ARM_ERR: swap id: object expected",C_result);
	XL_readStrCell(XL_receiveOrPay,C_receiveOrPay," ARM_ERR: receive or pay: string expected",C_result);
	XL_readStrCell(XL_xStyleId,C_xStyleId," ARM_ERR: exercise style id: object expected",C_result);
	XL_readStrCell(XL_kRefValId,C_kRefValId," ARM_ERR: strike reference value id: object expected",C_result);
	XL_readNumCell(XL_swapYearTerm,C_swapYearTerm," ARM_ERR: swap year term: date expected",C_result);
	
	if((receiveOrPayId = ARM_ConvRecOrPay (C_receiveOrPay, C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	long retCode;
	long objId;
	CCString prevClass;
	
	CCString curClass = LOCAL_SWAPTION_CLASS;
	CCString stringId;

	retCode = ARMLOCAL_EXOSWAPTION (LocalGetNumObjectId (C_swapId),
								   receiveOrPayId,
								   LocalGetNumObjectId (C_xStyleId),
								   LocalGetNumObjectId (C_kRefValId),
								   C_swapYearTerm,
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_PXL_EXOSWAPTION" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


__declspec(dllexport) LPXLOPER WINAPI Local_VARFIXSWAPTION (LPXLOPER XL_startDate,
															LPXLOPER XL_endDate,
															LPXLOPER XL_fixSpreads,
															LPXLOPER XL_XStyle,
															LPXLOPER XL_receiveOrPay,
															LPXLOPER XL_strike,
															LPXLOPER XL_maturity,
															LPXLOPER XL_liborType,
															LPXLOPER XL_spread,
															LPXLOPER XL_resetFreq,
															LPXLOPER XL_payFreq,
															LPXLOPER XL_ccy)
{
	ADD_LOG("Local_VARFIXSWAPTION ");
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

	CCString C_fixSpreads;
	CCString C_XStyle;

	CCString C_receiveOrPay;
	long receiveOrPayId;

	CCString C_liborType;
	long liborTypeId;

	double C_strike;
	double C_maturity;

	double C_spread;
	double C_spread_default = 0.0;

	CCString C_resetFreq;
	long resetFreqId;

	CCString C_payFreq;
	long payFreqId;

	CCString C_ccy;
	long ccyId;

	// error
	static int error;
	static char* reason = "";

	XL_readNumCell(XL_startDate,C_startDate," ARM_ERR: start date: date expected",C_result);
	XL_readNumCell(XL_endDate,C_endDate," ARM_ERR: end date: date expected",C_result);
	XL_readStrCell(XL_fixSpreads,C_fixSpreads," ARM_ERR: spreads : object expected",C_result);
	XL_readStrCell(XL_XStyle,C_XStyle," ARM_ERR: exercise style : object expected",C_result);
	XL_readStrCell(XL_receiveOrPay,C_receiveOrPay," ARM_ERR: receive or pay: string expected",C_result);
	XL_readStrCell(XL_liborType,C_liborType," ARM_ERR: libor type: string expected",C_result);
	XL_readNumCell(XL_strike,C_strike," ARM_ERR: strike: numeric expected",C_result);
	XL_readNumCell(XL_maturity,C_maturity," ARM_ERR: maturity: date expected",C_result);
	XL_readNumCellWD(XL_spread,C_spread,C_spread_default," ARM_ERR: spread: numeric expected",C_result);
	XL_readStrCellWD(XL_resetFreq,C_resetFreq,"-1"," ARM_ERR: reset frequency: string expected",C_result);
	XL_readStrCellWD(XL_payFreq,C_payFreq,"-1"," ARM_ERR: pay frequency: string expected",C_result);
	XL_readStrCellWD(XL_ccy,C_ccy,"DEFAULT"," ARM_ERR: currency: string expected",C_result);

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

	if(C_ccy == "DEFAULT")
	{
		if((liborTypeId == K_PIBOR1M) ||
		   (liborTypeId == K_PIBOR3M) ||
		   (liborTypeId == K_PIBOR6M) ||
		   (liborTypeId == K_PIBOR1Y))
		{
			ccyId = ARM_FRF_CCY_OBJECT;
		}
		else
		{
			ccyId = ARM_NULL_OBJECT;
		}
	}
	else
	{
		ccyId = LocalGetNumObjectId (C_ccy);
	}

	long retCode;
	long objId;
	CCString prevClass;

	CCString curClass = LOCAL_SWAPTION_CLASS;
	CCString stringId = GetLastCurCellEnvValue ();

	if(!stringId)
	{
		retCode = ARMLOCAL_VARFIXSWAPTION (C_startDate, C_endDate, 
										 LocalGetNumObjectId(C_fixSpreads), LocalGetNumObjectId(C_XStyle), 
										 receiveOrPayId,
										 C_strike, C_maturity, liborTypeId,
										 C_spread, resetFreqId,
										 payFreqId, ccyId, C_result);

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
			retCode = ARMLOCAL_VARFIXSWAPTION (C_startDate, C_endDate, 
										 LocalGetNumObjectId(C_fixSpreads), LocalGetNumObjectId(C_XStyle), 
										 receiveOrPayId,
										 C_strike, C_maturity, liborTypeId,
										 C_spread, resetFreqId,
										 payFreqId, ccyId, C_result, objId);

			if(retCode == ARM_OK)
			{
				LocalSetCurCellEnvValue (curClass, objId); 

				stringId = LocalMakeObjectId (objId, curClass);
			}
		}
		else
		{
			FreeCurCellContent ();
			retCode = ARMLOCAL_VARFIXSWAPTION (C_startDate, C_endDate, 
										 LocalGetNumObjectId(C_fixSpreads), LocalGetNumObjectId(C_XStyle), 
										 receiveOrPayId,
										 C_strike, C_maturity, liborTypeId,
										 C_spread, resetFreqId,
										 payFreqId, ccyId, C_result);
		
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_VARFIXSWAPTION" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}



__declspec(dllexport) LPXLOPER WINAPI Local_PXL_VARFIXSWAPTION (LPXLOPER XL_startDate,
																LPXLOPER XL_endDate,
																LPXLOPER XL_fixSpreads,
																LPXLOPER XL_XStyle,
																LPXLOPER XL_receiveOrPay,
																LPXLOPER XL_strike,
																LPXLOPER XL_maturity,
																LPXLOPER XL_liborType,
																LPXLOPER XL_spread,
																LPXLOPER XL_resetFreq,
																LPXLOPER XL_payFreq,
																LPXLOPER XL_ccy)
{
	ADD_LOG("Local_PXL_VARFIXSWAPTION ");
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

	CCString C_fixSpreads;
	CCString C_XStyle;

	CCString C_receiveOrPay;
	long receiveOrPayId;

	CCString C_liborType;
	long liborTypeId;

	double C_strike;
	double C_maturity;

	double C_spread;
	double C_spread_default = 0.0;

	CCString C_resetFreq;
	long resetFreqId;

	CCString C_payFreq;
	long payFreqId;

	CCString C_ccy;
	long ccyId;

	// error
	static int error;
	static char* reason = "";

	XL_readNumCell(XL_startDate,C_startDate," ARM_ERR: start date: date expected",C_result);
	XL_readNumCell(XL_endDate,C_endDate," ARM_ERR: end date: date expected",C_result);
	XL_readStrCell(XL_fixSpreads,C_fixSpreads," ARM_ERR: spreads : object expected",C_result);
	XL_readStrCell(XL_XStyle,C_XStyle," ARM_ERR: exercise style : object expected",C_result);
	XL_readStrCell(XL_receiveOrPay,C_receiveOrPay," ARM_ERR: receive or pay: string expected",C_result);
	XL_readStrCell(XL_liborType,C_liborType," ARM_ERR: libor type: string expected",C_result);
	XL_readNumCell(XL_strike,C_strike," ARM_ERR: strike: numeric expected",C_result);
	XL_readNumCell(XL_maturity,C_maturity," ARM_ERR: maturity: date expected",C_result);
	XL_readNumCellWD(XL_spread,C_spread,C_spread_default," ARM_ERR: spread: numeric expected",C_result);
	XL_readStrCellWD(XL_resetFreq,C_resetFreq,"-1"," ARM_ERR: reset frequency: string expected",C_result);
	XL_readStrCellWD(XL_payFreq,C_payFreq,"-1"," ARM_ERR: pay frequency: string expected",C_result);
	XL_readStrCellWD(XL_ccy,C_ccy,"DEFAULT"," ARM_ERR: currency: string expected",C_result);

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

	if(C_ccy == "DEFAULT")
	{
		if((liborTypeId == K_PIBOR1M) ||
		   (liborTypeId == K_PIBOR3M) ||
		   (liborTypeId == K_PIBOR6M) ||
		   (liborTypeId == K_PIBOR1Y))
		{
			ccyId = ARM_FRF_CCY_OBJECT;
		}
		else
		{
			ccyId = ARM_NULL_OBJECT;
		}
	}
	else
	{
		ccyId = LocalGetNumObjectId (C_ccy);
	}

	long retCode;
	long objId;
	CCString prevClass;

	CCString curClass = LOCAL_SWAPTION_CLASS;
	CCString stringId;

	retCode = ARMLOCAL_VARFIXSWAPTION (C_startDate, C_endDate, 
									 LocalGetNumObjectId(C_fixSpreads), LocalGetNumObjectId(C_XStyle), 
									 receiveOrPayId,
									 C_strike, C_maturity, liborTypeId,
									 C_spread, resetFreqId,
									 payFreqId, ccyId, C_result);

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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_PXL_VARFIXSWAPTION" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


__declspec(dllexport) LPXLOPER WINAPI Local_ARM_OPTIONALACCRUALZCBOND(LPXLOPER XL_startDate,
																	  LPXLOPER XL_endDate,
																	  LPXLOPER XL_strike,
																	  LPXLOPER XL_payfreq,
																	  LPXLOPER XL_nbCurPerforAcc,
																	  LPXLOPER XL_ccy)
{
	ADD_LOG("Local_ARM_OPTIONALACCRUALZCBOND");
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
	double C_strike;

	double C_nbCurPerforAcc;
	double C_nbCurPerforAcc_default = 0;

	CCString C_payFreq;
	long payFreqId;

	CCString C_ccy;
	long ccyId;
	
	// error
	static int error;
	static char* reason = "";

	XL_readNumCell(XL_startDate,C_startDate," ARM_ERR: start date: date expected",C_result);
	XL_readNumCell(XL_endDate,C_endDate," ARM_ERR: end date: date expected",C_result);
    XL_readNumCell(XL_strike,C_strike," ARM_ERR: strike : numeric expected",C_result);
	XL_readStrCell(XL_payfreq,C_payFreq," ARM_ERR: pay frequency: string expected",C_result);
    XL_readNumCellWD(XL_nbCurPerforAcc,C_nbCurPerforAcc,C_nbCurPerforAcc_default," ARM_ERR: nb Current Periods for Accrued : numeric expected",C_result);
	XL_readStrCellWD(XL_ccy,C_ccy,"DEFAULT"," ARM_ERR: currency: string expected",C_result);

	if(C_ccy == "DEFAULT")
	{
		ccyId = ARM_NULL_OBJECT;
	}
	else
	{
		ccyId = LocalGetNumObjectId (C_ccy);
	}

	if((payFreqId = ARM_ConvFrequency (C_payFreq, C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	long retCode;
	long objId;
	CCString prevClass;
	
	CCString curClass = LOCAL_OPTIONALACCRUALZC_CLASS;
	CCString stringId = GetLastCurCellEnvValue ();
	
	if(!stringId)
	{
		retCode = ARMLOCAL_ARM_OPTIONALACCRUALZCBOND (C_startDate,
													  C_endDate,
													  C_strike,
													  (long)C_nbCurPerforAcc,
													  payFreqId,
													  ccyId,
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
			retCode = ARMLOCAL_ARM_OPTIONALACCRUALZCBOND (C_startDate,
														  C_endDate,
														  C_strike,
														  (long)C_nbCurPerforAcc,
														  payFreqId,
														  ccyId,
														  C_result,
														  objId);

			if(retCode == ARM_OK)
			{
				objId = C_result.getLong ();

				stringId = LocalMakeObjectId (objId, curClass);
			}
		}
		else
		{
			FreeCurCellContent ();
			retCode = ARMLOCAL_ARM_OPTIONALACCRUALZCBOND (C_startDate,
														  C_endDate,
														  C_strike,
														  (long)C_nbCurPerforAcc,
														  payFreqId,
														  ccyId,
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_ARM_OPTIONALACCRUALZCBOND" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}



__declspec(dllexport) LPXLOPER WINAPI Local_PXL_ARM_OPTIONALACCRUALZCBOND(LPXLOPER XL_startDate,
																		  LPXLOPER XL_endDate,
																		  LPXLOPER XL_strike,
																		  LPXLOPER XL_payfreq,
																		  LPXLOPER XL_nbCurPerforAcc,
																		  LPXLOPER XL_ccy)
{
	ADD_LOG("Local_PXL_ARM_OPTIONALACCRUALZCBOND");
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
	double C_strike;

	double C_nbCurPerforAcc;
	double C_nbCurPerforAcc_default = 0;

	CCString C_payFreq;
	long payFreqId;

	CCString C_ccy;
    long ccyId;

	
	// error
	static int error;
	static char* reason = "";

	XL_readNumCell(XL_startDate,C_startDate," ARM_ERR: start date: date expected",C_result);
	XL_readNumCell(XL_endDate,C_endDate," ARM_ERR: end date: date expected",C_result);
    XL_readNumCell(XL_strike,C_strike," ARM_ERR: strike : numeric expected",C_result);
    XL_readNumCellWD(XL_nbCurPerforAcc,C_nbCurPerforAcc,C_nbCurPerforAcc_default," ARM_ERR: nb Current Periods for Accrued : numeric expected",C_result);
	XL_readStrCell(XL_payfreq,C_payFreq," ARM_ERR: pay frequency: string expected",C_result);
	XL_readStrCellWD(XL_ccy,C_ccy,"DEFAULT"," ARM_ERR: currency: string expected",C_result);

	if(C_ccy == "DEFAULT")
	{
		ccyId = ARM_NULL_OBJECT;
	}
	else
	{
		ccyId = LocalGetNumObjectId (C_ccy);
	}

	if((payFreqId = ARM_ConvFrequency (C_payFreq, C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	long retCode;
	long objId;
	
	CCString curClass = LOCAL_OPTIONALACCRUALZC_CLASS;
	CCString stringId;
	
	retCode = ARMLOCAL_ARM_OPTIONALACCRUALZCBOND (C_startDate,
												  C_endDate,
												  C_strike,
												  (long)C_nbCurPerforAcc,
												  payFreqId,
												  ccyId,
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_PXL_ARM_OPTIONALACCRUALZCBOND" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


__declspec(dllexport) LPXLOPER WINAPI Local_EXOCFSWAPTION (LPXLOPER XL_swapId,
														   LPXLOPER XL_receiveOrPay,
														   LPXLOPER XL_capOrFloor,
														   LPXLOPER XL_xStyleId,
														   LPXLOPER XL_kSptionRefValId,
														   LPXLOPER XL_kCFloorRefValId,
														   LPXLOPER XL_floorPosition,
														   LPXLOPER XL_IsBarrierCF)
{
	ADD_LOG("Local_EXOCFSWAPTION ");
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();


	// C variable
	CCString C_swapId;

	CCString C_receiveOrPay;
	long receiveOrPayId;

	CCString C_capOrFloor;
	long capOrFloorId;

	CCString C_xStyleId;
	CCString C_kSptionRefValId;
	CCString C_kCFloorRefValId;

	double C_floorPosition;
	
    double C_IsBarrierCF;
    double defaultFlag = 0;

	// error
	static int error;
	static char* reason = "";

	XL_readStrCell(XL_swapId,C_swapId," ARM_ERR: swap id: object expected",C_result);
	XL_readStrCell(XL_receiveOrPay,C_receiveOrPay," ARM_ERR: receive or pay: string expected",C_result);
	XL_readStrCell(XL_capOrFloor,C_capOrFloor," ARM_ERR: cap or floor: string expected",C_result);
	XL_readStrCell(XL_xStyleId,C_xStyleId," ARM_ERR: exercise style id: object expected",C_result);
	XL_readStrCell(XL_kSptionRefValId,C_kSptionRefValId," ARM_ERR: strike swaption reference value id: object expected",C_result);
	XL_readStrCell(XL_kCFloorRefValId,C_kCFloorRefValId," ARM_ERR: strike capfloor reference value id: object expected",C_result);
	XL_readNumCell(XL_floorPosition,C_floorPosition," ARM_ERR: floor position: numeric expected",C_result);
	XL_readNumCellWD(XL_IsBarrierCF,C_IsBarrierCF, defaultFlag, " ARM_ERR: Barrier or not Flag: numeric expected",C_result);

	if((receiveOrPayId = ARM_ConvRecOrPay (C_receiveOrPay, C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	if((capOrFloorId = ARM_ConvCapOrFloor (C_capOrFloor, C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	long retCode;
	long objId;
	CCString prevClass;
	
	CCString curClass = LOCAL_SWAPTION_CAPFLOOR_CLASS;
	CCString stringId = GetLastCurCellEnvValue ();

	if(!stringId)
	{
		retCode = ARMLOCAL_EXOCFSWAPTION (LocalGetNumObjectId (C_swapId),
										  receiveOrPayId,
										  capOrFloorId,
										  LocalGetNumObjectId (C_xStyleId),
										  LocalGetNumObjectId (C_kSptionRefValId),
										  LocalGetNumObjectId (C_kCFloorRefValId),
										  C_floorPosition,
										  (long) C_IsBarrierCF,
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
			retCode = ARMLOCAL_EXOCFSWAPTION (LocalGetNumObjectId (C_swapId),
											  receiveOrPayId,
											  capOrFloorId,
											  LocalGetNumObjectId (C_xStyleId),
											  LocalGetNumObjectId (C_kSptionRefValId),
											  LocalGetNumObjectId (C_kCFloorRefValId),
											  C_floorPosition,
											  (long) C_IsBarrierCF,
											  C_result,
											  objId);
		}
		else
		{
			FreeCurCellContent ();
			retCode = ARMLOCAL_EXOCFSWAPTION (LocalGetNumObjectId (C_swapId),
											  receiveOrPayId,
											  capOrFloorId,
											  LocalGetNumObjectId (C_xStyleId),
											  LocalGetNumObjectId (C_kSptionRefValId),
											  LocalGetNumObjectId (C_kCFloorRefValId),
											  C_floorPosition,
											  (long) C_IsBarrierCF,
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_EXOCFSWAPTION" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}



__declspec(dllexport) LPXLOPER WINAPI Local_PXL_EXOCFSWAPTION (LPXLOPER XL_swapId,
															   LPXLOPER XL_receiveOrPay,
															   LPXLOPER XL_capOrFloor,
															   LPXLOPER XL_xStyleId,
															   LPXLOPER XL_kSptionRefValId,
															   LPXLOPER XL_kCFloorRefValId,
															   LPXLOPER XL_floorPosition,
															   LPXLOPER XL_IsBarrierCF)
{
	ADD_LOG("Local_PXL_EXOCFSWAPTION ");
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();


	// C variable
	CCString C_swapId;

	CCString C_receiveOrPay;
	long receiveOrPayId;

	CCString C_capOrFloor;
	long capOrFloorId;

	CCString C_xStyleId;
	CCString C_kSptionRefValId;
	CCString C_kCFloorRefValId;

	double C_floorPosition;

    double C_IsBarrierCF;
    double defaultFlag = 0;

	// error
	static int error;
	static char* reason = "";

	XL_readStrCell(XL_swapId,C_swapId," ARM_ERR: swap id: object expected",C_result);
	XL_readStrCell(XL_receiveOrPay,C_receiveOrPay," ARM_ERR: receive or pay: string expected",C_result);
	XL_readStrCell(XL_capOrFloor,C_capOrFloor," ARM_ERR: cap or floor: string expected",C_result);
	XL_readStrCell(XL_xStyleId,C_xStyleId," ARM_ERR: exercise style id: object expected",C_result);
	XL_readStrCell(XL_kSptionRefValId,C_kSptionRefValId," ARM_ERR: strike swaption reference value id: object expected",C_result);
	XL_readStrCell(XL_kCFloorRefValId,C_kCFloorRefValId," ARM_ERR: strike capfloor reference value id: object expected",C_result);
	XL_readNumCell(XL_floorPosition,C_floorPosition," ARM_ERR: floor position: numeric expected",C_result);
    XL_readNumCellWD(XL_IsBarrierCF,C_IsBarrierCF, defaultFlag, " ARM_ERR: Barrier or not Flag: numeric expected",C_result);	

    if((receiveOrPayId = ARM_ConvRecOrPay (C_receiveOrPay, C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	if((capOrFloorId = ARM_ConvCapOrFloor (C_capOrFloor, C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	long retCode;
	long objId;
	
	CCString curClass = LOCAL_SWAPTION_CAPFLOOR_CLASS;
	CCString stringId;

	retCode = ARMLOCAL_EXOCFSWAPTION (LocalGetNumObjectId (C_swapId),
									  receiveOrPayId,
									  capOrFloorId,
									  LocalGetNumObjectId (C_xStyleId),
									  LocalGetNumObjectId (C_kSptionRefValId),
									  LocalGetNumObjectId (C_kCFloorRefValId),
									  C_floorPosition,
									  (long) C_IsBarrierCF,
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_PXL_EXOCFSWAPTION" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}



__declspec(dllexport) LPXLOPER WINAPI Local_ARM_FlexAccretSwaption (LPXLOPER XL_startDate,
																	LPXLOPER XL_endDate,
																	LPXLOPER XL_fixedRate,
																	LPXLOPER XL_nbCurrentPeriodsForAccrued,
																	LPXLOPER XL_receiveOrPay,
																	LPXLOPER XL_freq,
																	LPXLOPER XL_liborType,
																	LPXLOPER XL_spread,
																	LPXLOPER XL_exerciseDates,
																	LPXLOPER XL_currency)
{
	ADD_LOG("Local_ARM_FlexAccretSwaption ");
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
	double C_fixedRate;

	double C_nbCurrentPeriodsForAccrued;
	double C_nbCurrentPeriodsForAccrued_default = 0.0;

	CCString C_receiveOrPay;
	long receiveOrPayId;

	CCString C_freq;
	long freqId;

	CCString C_liborType;
	long liborTypeId;

	double C_spread;
	double C_spread_default = 0.0;

	CCString C_exerciseDates;
	long exerciseDatesId;

	CCString C_currency;
	long ccyId;

	// error
	static int error;
	static char* reason = "";

	XL_readNumCell(XL_startDate,C_startDate," ARM_ERR: start date: numeric expected",C_result);
	XL_readNumCell(XL_endDate,C_endDate," ARM_ERR: end date: numeric expected",C_result);
	XL_readNumCell(XL_fixedRate,C_fixedRate," ARM_ERR: fixed rate: numeric expected",C_result);
	XL_readNumCellWD(XL_nbCurrentPeriodsForAccrued,C_nbCurrentPeriodsForAccrued,C_nbCurrentPeriodsForAccrued_default," ARM_ERR: nbCurrentPeriodsForAccrued: numeric expected",C_result);
	XL_readStrCellWD(XL_receiveOrPay,C_receiveOrPay,"P"," ARM_ERR: receive or pay: string expected",C_result);
	XL_readStrCellWD(XL_freq,C_freq,"A"," ARM_ERR: frequency: string expected",C_result);
	XL_readStrCellWD(XL_liborType,C_liborType,"LIBOR3M"," ARM_ERR: libor type: object expected",C_result);
	XL_readNumCellWD(XL_spread,C_spread,C_spread_default, " ARM_ERR: spread: numeric expected",C_result);
	XL_readStrCellWD(XL_exerciseDates,C_exerciseDates,"NULL"," ARM_ERR: exercise style id: object expected",C_result);
	XL_readStrCellWD(XL_currency,C_currency,"DEFAULT"," ARM_ERR: currency: object expected",C_result);

	if((receiveOrPayId = ARM_ConvRecOrPay (C_receiveOrPay, C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	if((freqId = ARM_ConvFrequency (C_freq, C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	if((liborTypeId = ARM_ConvIrIndName (C_liborType, C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

    if ( C_exerciseDates == "NULL" )
	{
		exerciseDatesId = ARM_NULL_OBJECT;
	}
	else
	{
		exerciseDatesId = LocalGetNumObjectId(C_exerciseDates);
	}

    if ( C_currency == "DEFAULT" )
	{
		ccyId = ARM_NULL_OBJECT;
	}
	else
	{
		ccyId = LocalGetNumObjectId(C_currency);
	}

	long retCode;
	long objId;
	CCString prevClass;
	
	CCString curClass = LOCAL_FLEX_SWAPTION_CLASS;
	CCString stringId = GetLastCurCellEnvValue ();

	if(!stringId)
	{
		retCode = ARMLOCAL_FlexAccretSwaption (C_startDate,
											   C_endDate,
											   C_fixedRate,
											   (long) C_nbCurrentPeriodsForAccrued,
											   receiveOrPayId,
											   freqId,
											   liborTypeId,
											   C_spread,
											   exerciseDatesId,
											   ccyId,
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
			retCode = ARMLOCAL_FlexAccretSwaption (C_startDate,
												   C_endDate,
												   C_fixedRate,
												   (long) C_nbCurrentPeriodsForAccrued,
												   receiveOrPayId,
												   freqId,
												   liborTypeId,
												   C_spread,
												   exerciseDatesId,
												   ccyId,
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
			retCode = ARMLOCAL_FlexAccretSwaption (C_startDate,
												   C_endDate,
												   C_fixedRate,
												   (long) C_nbCurrentPeriodsForAccrued,
												   receiveOrPayId,
												   freqId,
												   liborTypeId,
												   C_spread,
												   exerciseDatesId,
												   ccyId,
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_ARM_FlexAccretSwaption" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}



__declspec(dllexport) LPXLOPER WINAPI Local_PXL_ARM_FlexAccretSwaption (LPXLOPER XL_startDate,
																		LPXLOPER XL_endDate,
																		LPXLOPER XL_fixedRate,
																		LPXLOPER XL_nbCurrentPeriodsForAccrued,
																		LPXLOPER XL_receiveOrPay,
																		LPXLOPER XL_freq,
																		LPXLOPER XL_liborType,
																		LPXLOPER XL_spread,
																		LPXLOPER XL_exerciseDates,
																		LPXLOPER XL_currency)
{
	ADD_LOG("Local_PXL_ARM_FlexAccretSwaption ");
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();


	// C variable
	// C variable
	double C_startDate;
	double C_endDate;
	double C_fixedRate;

	double C_nbCurrentPeriodsForAccrued;
	double C_nbCurrentPeriodsForAccrued_default = 0.0;

	CCString C_receiveOrPay;
	long receiveOrPayId;

	CCString C_freq;
	long freqId;

	CCString C_liborType;
	long liborTypeId;

	double C_spread;
	double C_spread_default = 0.0;

	CCString C_exerciseDates;
	long exerciseDatesId;

	CCString C_currency;
	long ccyId;

	// error
	static int error;
	static char* reason = "";

	XL_readNumCell(XL_startDate,C_startDate," ARM_ERR: start date: numeric expected",C_result);
	XL_readNumCell(XL_endDate,C_endDate," ARM_ERR: end date: numeric expected",C_result);
	XL_readNumCell(XL_fixedRate,C_fixedRate," ARM_ERR: fixed rate: numeric expected",C_result);
	XL_readNumCellWD(XL_nbCurrentPeriodsForAccrued,C_nbCurrentPeriodsForAccrued,C_nbCurrentPeriodsForAccrued_default," ARM_ERR: nbCurrentPeriodsForAccrued: numeric expected",C_result);
	XL_readStrCellWD(XL_receiveOrPay,C_receiveOrPay,"P"," ARM_ERR: receive or pay: string expected",C_result);
	XL_readStrCellWD(XL_freq,C_freq,"A"," ARM_ERR: frequency: string expected",C_result);
	XL_readStrCellWD(XL_liborType,C_liborType,"LIBOR3M"," ARM_ERR: libor type: object expected",C_result);
	XL_readNumCellWD(XL_spread,C_spread,C_spread_default, " ARM_ERR: spread: numeric expected",C_result);
	XL_readStrCellWD(XL_exerciseDates,C_exerciseDates,"NULL"," ARM_ERR: exercise style id: object expected",C_result);
	XL_readStrCellWD(XL_currency,C_currency,"DEFAULT"," ARM_ERR: currency: object expected",C_result);

	if((receiveOrPayId = ARM_ConvRecOrPay (C_receiveOrPay, C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	if((freqId = ARM_ConvFrequency (C_freq, C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	if((liborTypeId = ARM_ConvIrIndName (C_liborType, C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

    if ( C_exerciseDates == "NULL" )
	{
		exerciseDatesId = ARM_NULL_OBJECT;
	}
	else
	{
		exerciseDatesId = LocalGetNumObjectId(C_exerciseDates);
	}

    if ( C_currency == "DEFAULT" )
	{
		ccyId = ARM_NULL_OBJECT;
	}
	else
	{
		ccyId = LocalGetNumObjectId(C_currency);
	}

	long retCode;
	long objId;
	
	CCString curClass = LOCAL_FLEX_SWAPTION_CLASS;
	CCString stringId;

	retCode = ARMLOCAL_FlexAccretSwaption (C_startDate,
										   C_endDate,
										   C_fixedRate,
										   (long) C_nbCurrentPeriodsForAccrued,
										   receiveOrPayId,
										   freqId,
										   liborTypeId,
										   C_spread,
										   exerciseDatesId,
										   ccyId,
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_PXL_ARM_FlexAccretSwaption" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}

//------------------------------------------------------------------------------------------------------
__declspec(dllexport) LPXLOPER WINAPI Local_SwaptionStickyDelta(LPXLOPER XL_swaptionId,
																LPXLOPER XL_modelId,
																LPXLOPER XL_perturbeDiscountCurvId)
{
	ADD_LOG("Local_SwaptionStickyDelta");
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();



	// C variable
	CCString C_swaptionId;
	CCString C_modelId;
	CCString C_perturbeDiscountCurvId;
	
	// error
	static int error;
	static char* reason = "";
	
	XL_readStrCell(XL_swaptionId,C_swaptionId," ARM_ERR: forecast swaption id: object expected",C_result);
	XL_readStrCell(XL_modelId,C_modelId," ARM_ERR: forecast model id: object expected",C_result);
	XL_readStrCell(XL_perturbeDiscountCurvId,C_perturbeDiscountCurvId," ARM_ERR: forecast perturbe discount curve id: object expected",C_result);
	
	long retCode;

	retCode = ARMLOCAL_SwaptionStickyDelta(LocalGetNumObjectId(C_swaptionId), 
										   LocalGetNumObjectId(C_modelId),
										   LocalGetNumObjectId(C_perturbeDiscountCurvId),
										   C_result);

	if(retCode == ARM_OK)
	{
		FreeCurCellErr ();
		XL_result.xltype = xltypeNum;
		XL_result.val.num = C_result.getDouble();
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_SwaptionStickyDelta" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}



/*---- End Of File ----*/

// EOF %M% 


