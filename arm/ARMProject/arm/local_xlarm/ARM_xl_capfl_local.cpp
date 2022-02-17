#pragma warning(disable : 4786)
#pragma warning(disable : 4005)

#include <libCCxll\CCxll.h>
#include "ARM_xl_trycatch_local.h"
#include <util\fromto.h>

#include "ARM_local_interglob.h"
#include "ARM_local_interface.h"

#include <ARM\libarm_local\ARM_local_class.h>
#include <ARM\libarm_local\ARM_local_capfl.h>

#include "util\tech_macro.h"

__declspec(dllexport) LPXLOPER WINAPI Local_LIBORCF (LPXLOPER XL_startDate,
													 LPXLOPER XL_endDate,
													 LPXLOPER XL_isItCapOrFloor,
													 LPXLOPER XL_strike,
													 LPXLOPER XL_liborType,
													 LPXLOPER XL_spread,
													 LPXLOPER XL_resetFreq,
													 LPXLOPER XL_payFreq,
													 LPXLOPER XL_ccy)
{
	ADD_LOG("Local_LIBORCF ");
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

	CCString C_isItCapOrFloor;
	double isItCapOrFloorId;
	
	double   C_strike_double;
	CCString C_strike_str;
	long     strikeType;

	CCString C_liborType;
	long liborTypeId;

	double C_spread;

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
	XL_readStrCell(XL_isItCapOrFloor,C_isItCapOrFloor," ARM_ERR: cap or floor: string expected",C_result);
	XL_readStrOrNumCell(XL_strike, C_strike_str, C_strike_double, strikeType,
		   " ARM_ERR: strike: numeric or object ID string expected",C_result);
	XL_readStrCell(XL_liborType,C_liborType," ARM_ERR: libor type: string expected",C_result);
	XL_readNumCell(XL_spread,C_spread," ARM_ERR: spread: numeric expected",C_result);
	XL_readStrCellWD(XL_resetFreq,C_resetFreq,"-1"," ARM_ERR: reset frequency: string expected",C_result);
	XL_readStrCellWD(XL_payFreq,C_payFreq,"-1"," ARM_ERR: pay frequency: string expected",C_result);
	XL_readStrCellWD(XL_ccy,C_ccy,"DEFAULT"," ARM_ERR: currency: string expected",C_result);

	if((liborTypeId = ARM_ConvIrIndName (C_liborType, C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}
		
	if ((C_ccy.GetLen() > 3)
		&&
		!(C_ccy == "DEFAULT"))
		ccyIsObject = true;

	if((isItCapOrFloorId = ARM_ConvCapOrFloor (C_isItCapOrFloor, C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	if((resetFreqId = ARM_ConvFrequency (C_resetFreq, C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	if((payFreqId = ARM_ConvFrequency (C_payFreq, C_result)) == ARM_DEFAULT_ERR)
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
	
	CCString curClass = LOCAL_CAPFLOOR_CLASS;
	CCString stringId = GetLastCurCellEnvValue ();
	
	if(!stringId)
	{
		retCode = ARMLOCAL_LIBORCF (C_startDate,
									C_endDate,
									isItCapOrFloorId,
									strikeType,
									C_strike_double,
									liborTypeId,
									C_spread,
									resetFreqId,
									payFreqId,
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
			retCode = ARMLOCAL_LIBORCF (C_startDate,
										C_endDate,
										isItCapOrFloorId,
										strikeType,
										C_strike_double,
										liborTypeId,
										C_spread,
										resetFreqId,
										payFreqId,
										ccyIsObject,
										C_ccy,
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
			retCode = ARMLOCAL_LIBORCF (C_startDate,
										C_endDate,
										isItCapOrFloorId,
										strikeType,
										C_strike_double,
										liborTypeId,
										C_spread,
										resetFreqId,
										payFreqId,
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_LIBORCF" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;

}


__declspec(dllexport) LPXLOPER WINAPI Local_PXL_LIBORCF (LPXLOPER XL_startDate,
														 LPXLOPER XL_endDate,
														 LPXLOPER XL_isItCapOrFloor,
														 LPXLOPER XL_strike,
														 LPXLOPER XL_liborType,
														 LPXLOPER XL_spread,
														 LPXLOPER XL_resetFreq,
														 LPXLOPER XL_payFreq,
														 LPXLOPER XL_ccy)
{
	ADD_LOG("Local_PXL_LIBORCF ");
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

	CCString C_isItCapOrFloor;
	long isItCapOrFloorId;
	
	double   C_strike_double;
	CCString C_strike_str;
	long     strikeType;

	CCString C_liborType;
	long liborTypeId;

	double C_spread;

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
	XL_readStrCell(XL_isItCapOrFloor,C_isItCapOrFloor," ARM_ERR: cap or floor: string expected",C_result);
	XL_readStrOrNumCell(XL_strike, C_strike_str, C_strike_double, strikeType,
		   " ARM_ERR: strike: numeric or object ID string expected",C_result);
	XL_readStrCell(XL_liborType,C_liborType," ARM_ERR: libor type: string expected",C_result);
	XL_readNumCell(XL_spread,C_spread," ARM_ERR: spread: numeric expected",C_result);
	XL_readStrCellWD(XL_resetFreq,C_resetFreq,"-1"," ARM_ERR: reset frequency: string expected",C_result);
	XL_readStrCellWD(XL_payFreq,C_payFreq,"-1"," ARM_ERR: pay frequency: string expected",C_result);
	XL_readStrCellWD(XL_ccy,C_ccy,"DEFAULT"," ARM_ERR: currency: string expected",C_result);

	if((liborTypeId = ARM_ConvIrIndName (C_liborType, C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}
		
	if ((C_ccy.GetLen() > 3)
		&&
		!(C_ccy == "DEFAULT"))
		ccyIsObject = true;

	if((isItCapOrFloorId = ARM_ConvCapOrFloor (C_isItCapOrFloor, C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	if((resetFreqId = ARM_ConvFrequency (C_resetFreq, C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	if((payFreqId = ARM_ConvFrequency (C_payFreq, C_result)) == ARM_DEFAULT_ERR)
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
	
	CCString curClass = LOCAL_CAPFLOOR_CLASS;
	CCString stringId;
	
	retCode = ARMLOCAL_LIBORCF (C_startDate,
								C_endDate,
								isItCapOrFloorId,
								strikeType,
								C_strike_double,
								liborTypeId,
								C_spread,
								resetFreqId,
								payFreqId,
								ccyIsObject,
								C_ccy,
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_PXL_LIBORCF" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


__declspec(dllexport) LPXLOPER WINAPI Local_LIBORFLEXCF (LPXLOPER XL_startDate,
														 LPXLOPER XL_endDate,
														 LPXLOPER XL_isItCapOrFloor,
														 LPXLOPER XL_strike,
														 LPXLOPER XL_nbEx,
														 LPXLOPER XL_exerciseType,
														 LPXLOPER XL_liborType,
														 LPXLOPER XL_spread,
														 LPXLOPER XL_resetFreq,
														 LPXLOPER XL_payFreq,
														 LPXLOPER XL_ccy)
{
	ADD_LOG("Local_LIBORFLEXCF ");
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

	CCString C_isItCapOrFloor;
	long isItCapOrFloorId;
	
	double C_strike;

	double C_nbEx;

	CCString C_exerciseType;
	long exerciseTypeId;

	CCString C_liborType;
	long liborTypeId;

	double C_spread;

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
	XL_readStrCell(XL_isItCapOrFloor,C_isItCapOrFloor," ARM_ERR: cap or floor: string expected",C_result);
	XL_readNumCell(XL_strike,C_strike," ARM_ERR: strike: numeric expected",C_result);
	XL_readNumCell(XL_nbEx,C_nbEx," ARM_ERR: number of exercises: numeric expected",C_result);
	XL_readStrCell(XL_exerciseType,C_exerciseType," ARM_ERR: reset frequency: string expected",C_result);
	XL_readStrCell(XL_liborType,C_liborType," ARM_ERR: libor type: string expected",C_result);
	XL_readNumCell(XL_spread,C_spread," ARM_ERR: spread: numeric expected",C_result);
	XL_readStrCellWD(XL_resetFreq,C_resetFreq,"-1"," ARM_ERR: reset frequency: string expected",C_result);
	XL_readStrCellWD(XL_payFreq,C_payFreq,"-1"," ARM_ERR: pay frequency: string expected",C_result);
	XL_readStrCellWD(XL_ccy,C_ccy,"DEFAULT"," ARM_ERR: currency: string expected",C_result);

	if((liborTypeId = ARM_ConvIrIndName (C_liborType, C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}
	
	if((exerciseTypeId = ARM_ConvExerciseType (C_exerciseType, C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}
		
	if((isItCapOrFloorId = ARM_ConvCapOrFloor (C_isItCapOrFloor, C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	if((resetFreqId = ARM_ConvFrequency (C_resetFreq, C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	if((payFreqId = ARM_ConvFrequency (C_payFreq, C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	if ((C_ccy.GetLen() > 3)
		&&
		!(C_ccy == "DEFAULT"))
		ccyIsObject = true;

	long retCode;
	long objId;
	CCString prevClass;
	
	CCString curClass = LOCAL_FLEXIBLECAPFLOOR_CLASS;
	CCString stringId = GetLastCurCellEnvValue ();
	
	if(!stringId)
	{
		retCode = ARMLOCAL_LIBORFLEXCF (C_startDate, C_endDate, isItCapOrFloorId,
							       C_strike, (long)C_nbEx, exerciseTypeId,
								   liborTypeId, C_spread, resetFreqId, payFreqId,
								   ccyIsObject, C_ccy, C_result);

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
			retCode = ARMLOCAL_LIBORFLEXCF (C_startDate, C_endDate, isItCapOrFloorId,
							           C_strike, (long)C_nbEx, exerciseTypeId,
								       liborTypeId, C_spread, resetFreqId, payFreqId,
								       ccyIsObject, C_ccy, C_result, objId);

			if(retCode == ARM_OK)
			{
				LocalSetCurCellEnvValue (curClass, objId); 

				stringId = LocalMakeObjectId (objId, curClass);
			}
		}
		else
		{
			FreeCurCellContent ();
			retCode = ARMLOCAL_LIBORFLEXCF (C_startDate, C_endDate, isItCapOrFloorId,
							           C_strike, (long)C_nbEx, exerciseTypeId,
								       liborTypeId, C_spread, resetFreqId, payFreqId,
								       ccyIsObject, C_ccy, C_result);
		
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_LIBORFLEXCF" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


__declspec(dllexport) LPXLOPER WINAPI Local_PXL_LIBORFLEXCF (LPXLOPER XL_startDate,
															 LPXLOPER XL_endDate,
															 LPXLOPER XL_isItCapOrFloor,
															 LPXLOPER XL_strike,
															 LPXLOPER XL_nbEx,
															 LPXLOPER XL_exerciseType,
															 LPXLOPER XL_liborType,
															 LPXLOPER XL_spread,
															 LPXLOPER XL_resetFreq,
															 LPXLOPER XL_payFreq,
															 LPXLOPER XL_ccy)
{
	ADD_LOG("Local_PXL_LIBORFLEXCF ");
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

	CCString C_isItCapOrFloor;
	long isItCapOrFloorId;
	
	double C_strike;

	double C_nbEx;

	CCString C_exerciseType;
	long exerciseTypeId;

	CCString C_liborType;
	long liborTypeId;

	double C_spread;

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
	XL_readStrCell(XL_isItCapOrFloor,C_isItCapOrFloor," ARM_ERR: cap or floor: string expected",C_result);
	XL_readNumCell(XL_strike,C_strike," ARM_ERR: strike: numeric expected",C_result);
	XL_readNumCell(XL_nbEx,C_nbEx," ARM_ERR: number of exercises: numeric expected",C_result);
	XL_readStrCell(XL_exerciseType,C_exerciseType," ARM_ERR: reset frequency: string expected",C_result);
	XL_readStrCell(XL_liborType,C_liborType," ARM_ERR: libor type: string expected",C_result);
	XL_readNumCell(XL_spread,C_spread," ARM_ERR: spread: numeric expected",C_result);
	XL_readStrCellWD(XL_resetFreq,C_resetFreq,"-1"," ARM_ERR: reset frequency: string expected",C_result);
	XL_readStrCellWD(XL_payFreq,C_payFreq,"-1"," ARM_ERR: pay frequency: string expected",C_result);
	XL_readStrCellWD(XL_ccy,C_ccy,"DEFAULT"," ARM_ERR: currency: string expected",C_result);

	if((liborTypeId = ARM_ConvIrIndName (C_liborType, C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}
	
	if((exerciseTypeId = ARM_ConvExerciseType (C_exerciseType, C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}
		
	if((isItCapOrFloorId = ARM_ConvCapOrFloor (C_isItCapOrFloor, C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	if((resetFreqId = ARM_ConvFrequency (C_resetFreq, C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	if((payFreqId = ARM_ConvFrequency (C_payFreq, C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	if ((C_ccy.GetLen() > 3)
		&&
		!(C_ccy == "DEFAULT"))
		ccyIsObject = true;

	long retCode;
	long objId;
	
	CCString curClass = LOCAL_FLEXIBLECAPFLOOR_CLASS;
	CCString stringId;
	
	retCode = ARMLOCAL_LIBORFLEXCF (C_startDate, C_endDate, isItCapOrFloorId,
							   C_strike, (long)C_nbEx, exerciseTypeId,
							   liborTypeId, C_spread, resetFreqId, payFreqId,
							   ccyIsObject, C_ccy, C_result);

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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_PXL_LIBORFLEXCF" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}




__declspec(dllexport) LPXLOPER WINAPI Local_CAPFLOOR (LPXLOPER XL_swapLeg,
													  LPXLOPER XL_isItCapOrFloor,
													  LPXLOPER XL_strike)
{
	ADD_LOG("Local_CAPFLOOR ");
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();

	// C variable
	CCString C_swapLeg;
	
	CCString C_isItCapOrFloor;

    long isItCapOrFloorId;

	double   C_strike_double;
	CCString C_strike_str;
	long     strikeType;

	// error
	static int error;
	static char* reason = "";

	XL_readStrCell(XL_swapLeg,C_swapLeg," ARM_ERR: swap leg id: object expected",C_result);
	XL_readStrCell(XL_isItCapOrFloor,C_isItCapOrFloor," ARM_ERR: cap or floor: string expected",C_result);
	XL_readStrOrNumCell(XL_strike, C_strike_str, C_strike_double, strikeType,
		   " ARM_ERR: strike: numeric or object ID string expected",C_result);
	
	if((isItCapOrFloorId = ARM_ConvCapOrFloor (C_isItCapOrFloor, C_result)) == ARM_DEFAULT_ERR)
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
	
	CCString curClass = CorrespondingCapFloorClasswId( LocalGetNumObjectId (C_swapLeg) );
	CCString stringId = GetLastCurCellEnvValue ();
	
	if(!stringId)
	{
		retCode = ARMLOCAL_CAPFLOOR (LocalGetNumObjectId (C_swapLeg),
									 isItCapOrFloorId,
									 strikeType,
									 C_strike_double,
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
			retCode = ARMLOCAL_CAPFLOOR (LocalGetNumObjectId (C_swapLeg),
										 isItCapOrFloorId,
										 strikeType,
										 C_strike_double,
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
			retCode = ARMLOCAL_CAPFLOOR (LocalGetNumObjectId (C_swapLeg),
										 isItCapOrFloorId,
										 strikeType,
										 C_strike_double,
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_CAPFLOOR" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


__declspec(dllexport) LPXLOPER WINAPI Local_PXL_CAPFLOOR (LPXLOPER XL_swapLeg,
														  LPXLOPER XL_isItCapOrFloor,
														  LPXLOPER XL_strike)
{
	ADD_LOG("Local_PXL_CAPFLOOR ");
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();

	// C variable
	CCString C_swapLeg;
	
	CCString C_isItCapOrFloor;

    long isItCapOrFloorId;

	double   C_strike_double;
	CCString C_strike_str;
	long     strikeType;

	// error
	static int error;
	static char* reason = "";

	XL_readStrCell(XL_swapLeg,C_swapLeg," ARM_ERR: swap leg id: object expected",C_result);
	XL_readStrCell(XL_isItCapOrFloor,C_isItCapOrFloor," ARM_ERR: cap or floor: string expected",C_result);
	XL_readStrOrNumCell(XL_strike, C_strike_str, C_strike_double, strikeType,
		   " ARM_ERR: strike: numeric or object ID string expected",C_result);
	
	if((isItCapOrFloorId = ARM_ConvCapOrFloor (C_isItCapOrFloor, C_result)) == ARM_DEFAULT_ERR)
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
	
	CCString curClass = LOCAL_CAPFLOOR_CLASS;
	CCString stringId;
	
	retCode = ARMLOCAL_CAPFLOOR (LocalGetNumObjectId (C_swapLeg),
								 isItCapOrFloorId,
								 strikeType,
								 C_strike_double,
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_PXL_CAPFLOOR" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


__declspec(dllexport) LPXLOPER WINAPI Local_FLEXCF (LPXLOPER XL_swapLegId,
													LPXLOPER XL_isItCapOrFloor,
													LPXLOPER XL_strike,
													LPXLOPER XL_nbEx,
													LPXLOPER XL_exerciseType)
{
	ADD_LOG("Local_FLEXCF ");
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();

	// C variable
	CCString C_swapLegId;

	CCString C_isItCapOrFloor;
	long isItCapOrFloorId;
	
	double C_strike;

	double C_nbEx;

	CCString C_exerciseType;
	long exerciseTypeId;

	// error
	static int error;
	static char* reason = "";

	XL_readStrCell(XL_swapLegId,C_swapLegId," ARM_ERR: swap leg id: object expected",C_result);
	XL_readStrCell(XL_isItCapOrFloor,C_isItCapOrFloor," ARM_ERR: cap or floor: string expected",C_result);
	XL_readNumCell(XL_strike,C_strike," ARM_ERR: strike: numeric expected",C_result);
	XL_readNumCell(XL_nbEx,C_nbEx," ARM_ERR: number of exercises: numeric expected",C_result);
	XL_readStrCell(XL_exerciseType,C_exerciseType," ARM_ERR: reset frequency: string expected",C_result);
	
	if((exerciseTypeId = ARM_ConvExerciseType (C_exerciseType, C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}
		
	if((isItCapOrFloorId = ARM_ConvCapOrFloor (C_isItCapOrFloor, C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	long retCode;
	long objId;
	CCString prevClass;
	
	CCString curClass = LOCAL_FLEXIBLECAPFLOOR_CLASS;
	CCString stringId = GetLastCurCellEnvValue ();
	
	if(!stringId)
	{
		retCode = ARMLOCAL_FLEXCF (LocalGetNumObjectId (C_swapLegId), isItCapOrFloorId,
							  C_strike, (long)C_nbEx, exerciseTypeId, 
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
			retCode = ARMLOCAL_FLEXCF (LocalGetNumObjectId (C_swapLegId), isItCapOrFloorId,
							      C_strike, (long)C_nbEx, exerciseTypeId,
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
			retCode = ARMLOCAL_FLEXCF (LocalGetNumObjectId (C_swapLegId), isItCapOrFloorId,
							      C_strike, (long)C_nbEx, exerciseTypeId, 
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_FLEXCF" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


__declspec(dllexport) LPXLOPER WINAPI Local_PXL_FLEXCF (LPXLOPER XL_swapLegId,
														LPXLOPER XL_isItCapOrFloor,
														LPXLOPER XL_strike,
														LPXLOPER XL_nbEx,
														LPXLOPER XL_exerciseType)
{
	ADD_LOG("Local_PXL_FLEXCF ");
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();

	// C variable
	CCString C_swapLegId;

	CCString C_isItCapOrFloor;
	long isItCapOrFloorId;
	
	double C_strike;

	double C_nbEx;

	CCString C_exerciseType;
	long exerciseTypeId;

	// error
	static int error;
	static char* reason = "";

	XL_readStrCell(XL_swapLegId,C_swapLegId," ARM_ERR: swap leg id: object expected",C_result);
	XL_readStrCell(XL_isItCapOrFloor,C_isItCapOrFloor," ARM_ERR: cap or floor: string expected",C_result);
	XL_readNumCell(XL_strike,C_strike," ARM_ERR: strike: numeric expected",C_result);
	XL_readNumCell(XL_nbEx,C_nbEx," ARM_ERR: number of exercises: numeric expected",C_result);
	XL_readStrCell(XL_exerciseType,C_exerciseType," ARM_ERR: reset frequency: string expected",C_result);
	
	if((exerciseTypeId = ARM_ConvExerciseType (C_exerciseType, C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}
		
	if((isItCapOrFloorId = ARM_ConvCapOrFloor (C_isItCapOrFloor, C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	long retCode;
	long objId;
	
	CCString curClass = LOCAL_FLEXIBLECAPFLOOR_CLASS;
	CCString stringId;
	
	retCode = ARMLOCAL_FLEXCF (LocalGetNumObjectId (C_swapLegId), isItCapOrFloorId,
						  C_strike, (long)C_nbEx, exerciseTypeId, 
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_PXL_FLEXCF" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


__declspec(dllexport) LPXLOPER WINAPI Local_MATCAPFLOOR(LPXLOPER XL_swapLegId,
														LPXLOPER XL_annuity,
														LPXLOPER XL_initNominal,
														LPXLOPER XL_isTRI,
														LPXLOPER XL_capOrFloor,
														LPXLOPER XL_coeff,
														LPXLOPER XL_firstTRIstrike,
														LPXLOPER XL_minStrikes,
														LPXLOPER XL_isDigitalPayoff,
														LPXLOPER XL_increasingCoef,
                                                        LPXLOPER XL_maxMaturityDate)
{
	ADD_LOG("Local_MATCAPFLOOR");
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();

	// C variable
	CCString C_swapLeg;
    
    double C_annuity;
    double C_initNominal;
    double C_isTRI; 
    double defTRI = 1.0;
    
    CCString C_capOrFloor;
    long isItCapOrFloorId;
    double C_coeff;
    double defCoeff = -1.0;

    double C_firstTRIstrike;
    double defStrike = -1.0;

    CCString C_minStrikes;
	long minStrikesId;

    double C_isDigitalPayoff;
    double C_isDigitalPayoff_default = 0.0;

    double C_increasingCoef;
    double C_increasingCoef_default = 1.0;

    double C_MaxMatDate_default = -1.0;
    double C_MaxMatDate;

	// error
	static int error;
	static char* reason = "";

	XL_readStrCell(XL_swapLegId,C_swapLeg," ARM_ERR: swap leg id: object expected",C_result);

	XL_readNumCell(XL_annuity,C_annuity," ARM_ERR: annuity: numeric expected",C_result);
	XL_readNumCell(XL_initNominal,C_initNominal," ARM_ERR: nominal: numeric expected",C_result);

	XL_readNumCellWD(XL_isTRI,C_isTRI,defTRI," ARM_ERR: isTRI: numeric expected",C_result);

	XL_readStrCellWD(XL_capOrFloor,C_capOrFloor,"C"," ARM_ERR: cap or floor: string expected",C_result);

	XL_readNumCellWD(XL_coeff,C_coeff, defCoeff," ARM_ERR: Coeff: numeric expected",C_result);
	XL_readNumCellWD(XL_firstTRIstrike,C_firstTRIstrike, defStrike," ARM_ERR: strike: numeric expected",C_result);

	XL_readStrCellWD(XL_minStrikes,C_minStrikes,"DEFAULT"," ARM_ERR: min Strikes: object expected",C_result);
	XL_readNumCellWD(XL_isDigitalPayoff,C_isDigitalPayoff, C_isDigitalPayoff_default," ARM_ERR: isDigital Payoff: numeric expected",C_result);
	XL_readNumCellWD(XL_increasingCoef,C_increasingCoef, C_increasingCoef_default," ARM_ERR: increasing coefficient: numeric expected",C_result);
	XL_readNumCellWD(XL_maxMaturityDate,C_MaxMatDate, C_MaxMatDate_default," ARM_ERR: max maturity date : date expected",C_result);

	if((isItCapOrFloorId = ARM_ConvCapOrFloor (C_capOrFloor, C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	if(C_minStrikes == "DEFAULT")
	{
		minStrikesId = ARM_NULL_OBJECT;
	}
	else
	{
		minStrikesId = LocalGetNumObjectId (C_minStrikes);
	}
	
	long retCode;
	long objId;
	CCString prevClass;
	
	CCString curClass = LOCAL_MATCAPFLOOR_CLASS;
	CCString stringId = GetLastCurCellEnvValue ();
	
	if(!stringId)
	{
		retCode = ARMLOCAL_MATCAPFLOOR (LocalGetNumObjectId(C_swapLeg),
										C_annuity,
										C_initNominal,
										(long) C_isTRI,
										isItCapOrFloorId,
										C_coeff,
										C_firstTRIstrike,
										minStrikesId,
										(long) C_isDigitalPayoff,
										C_increasingCoef,
                                        C_MaxMatDate,
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
			retCode = ARMLOCAL_MATCAPFLOOR (LocalGetNumObjectId(C_swapLeg),
											C_annuity,
											C_initNominal,
											(long) C_isTRI,
											isItCapOrFloorId,
											C_coeff,
											C_firstTRIstrike,
											minStrikesId,
											(long) C_isDigitalPayoff,
											C_increasingCoef,
                                            C_MaxMatDate,
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
            
			retCode = ARMLOCAL_MATCAPFLOOR (LocalGetNumObjectId(C_swapLeg),
											C_annuity,
											C_initNominal,
											(long) C_isTRI,
											isItCapOrFloorId,
											C_coeff,
											C_firstTRIstrike,
											minStrikesId,
											(long) C_isDigitalPayoff,
											C_increasingCoef,
                                            C_MaxMatDate,
											C_result);
		
			if ( retCode == ARM_OK )
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_MATCAPFLOOR" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}



__declspec(dllexport) LPXLOPER WINAPI Local_PXL_MATCAPFLOOR(LPXLOPER XL_swapLegId,
															LPXLOPER XL_annuity,
															LPXLOPER XL_initNominal,
															LPXLOPER XL_isTRI,
															LPXLOPER XL_capOrFloor,
															LPXLOPER XL_coeff,
															LPXLOPER XL_firstTRIstrike,
															LPXLOPER XL_minStrikes,
															LPXLOPER XL_isDigitalPayoff,
															LPXLOPER XL_increasingCoef,
                                                            LPXLOPER XL_maxMaturityDate)
{
	ADD_LOG("Local_PXL_MATCAPFLOOR");
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();

	// C variable

	CCString C_swapLeg;

    double C_annuity;
    double C_initNominal;
    double C_isTRI; 
    double defTRI = 1;

    CCString C_capOrFloor;
    long isItCapOrFloorId;
    double C_coeff;
    double defCoeff = -1;

    double C_firstTRIstrike;
    double defStrike = -1;

    CCString C_minStrikes;
	long minStrikesId;

    double C_isDigitalPayoff;
    double C_isDigitalPayoff_default = 0.0;

    double C_increasingCoef;
    double C_increasingCoef_default = 1.0;

    double C_MaxMatDate_default = -1.0;
    double C_MaxMatDate;

	// error
	static int error;
	static char* reason = "";

	XL_readStrCell(XL_swapLegId,C_swapLeg," ARM_ERR: swap leg id: object expected",C_result);
	
    XL_readNumCell(XL_annuity,C_annuity," ARM_ERR: annuity: numeric expected",C_result);
    XL_readNumCell(XL_initNominal,C_initNominal," ARM_ERR: nominal: numeric expected",C_result);

    XL_readNumCellWD(XL_isTRI,C_isTRI,defTRI," ARM_ERR: isTRI: numeric expected",C_result);

    XL_readStrCellWD(XL_capOrFloor,C_capOrFloor,"C"," ARM_ERR: cap or floor: string expected",C_result);

    XL_readNumCellWD(XL_coeff,C_coeff, defCoeff," ARM_ERR: Coeff: numeric expected",C_result);
    XL_readNumCellWD(XL_firstTRIstrike,C_firstTRIstrike, defStrike," ARM_ERR: strike: numeric expected",C_result);
	
    XL_readStrCellWD(XL_minStrikes,C_minStrikes,"DEFAULT"," ARM_ERR: min Strikes: object expected",C_result);
	XL_readNumCellWD(XL_isDigitalPayoff,C_isDigitalPayoff, C_isDigitalPayoff_default," ARM_ERR: isDigital Payoff: numeric expected",C_result);
	XL_readNumCellWD(XL_increasingCoef,C_increasingCoef, C_increasingCoef_default," ARM_ERR: increasing coefficient: numeric expected",C_result);
	XL_readNumCellWD(XL_maxMaturityDate,C_MaxMatDate, C_MaxMatDate_default," ARM_ERR: max maturity date : date expected",C_result);

	if((isItCapOrFloorId = ARM_ConvCapOrFloor (C_capOrFloor, C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	if(C_minStrikes == "DEFAULT")
	{
		minStrikesId = ARM_NULL_OBJECT;
	}
	else
	{
		minStrikesId = LocalGetNumObjectId (C_minStrikes);
	}

	long retCode;
	long objId;
	CCString prevClass;
	
	CCString curClass = LOCAL_MATCAPFLOOR_CLASS;
	CCString stringId;


    retCode = ARMLOCAL_MATCAPFLOOR (LocalGetNumObjectId(C_swapLeg),
									C_annuity,
									C_initNominal,
									(long) C_isTRI,
									isItCapOrFloorId,
									C_coeff,
									C_firstTRIstrike,
									minStrikesId,
									(long) C_isDigitalPayoff,
									C_increasingCoef,
                                    C_MaxMatDate,
    								C_result);

	if ( retCode == ARM_OK )
	{
	   objId = C_result.getLong (); 

	   stringId = LocalMakeObjectId (objId, curClass);
	}


	if ( retCode == ARM_OK )
	{
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_PXL_MATCAPFLOOR" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}



__declspec(dllexport) LPXLOPER WINAPI Local_STICKY (LPXLOPER XL_swapLeg,
													LPXLOPER XL_isItCapOrFloor,
													LPXLOPER XL_strike,
													LPXLOPER XL_spreadDates,
													LPXLOPER XL_spreadValues,
													LPXLOPER XL_kRefVal)
{
	ADD_LOG("Local_STICKY ");
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();

	// C variable
	CCString C_swapLeg;

	CCString C_isItCapOrFloor;
	long isItCapOrFloorId;

	double C_strike;

	VECTOR<double> C_spreadDates;
	VECTOR<double> C_spreadValues;

	CCString C_kRefVal;
	long kRefValId;

	// error
	static int error;
	static char* reason = "";

	XL_readStrCell(XL_swapLeg,C_swapLeg," ARM_ERR: swap leg id: object expected",C_result);
	XL_readStrCell(XL_isItCapOrFloor,C_isItCapOrFloor," ARM_ERR: cap or floor: string expected",C_result);
	XL_readNumCell(XL_strike,C_strike," ARM_ERR: strike: numeric expected",C_result);
	XL_readNumVector(XL_spreadDates,C_spreadDates," ARM_ERR: spread dates: array of date expected",C_result);
	XL_readNumVector(XL_spreadValues,C_spreadValues," ARM_ERR: spread values: array of numeric expected",C_result);
	XL_readStrCellWD(XL_kRefVal,C_kRefVal,"DEFAULT"," ARM_ERR: strike reference value id: object expected",C_result);
	
	if((isItCapOrFloorId = ARM_ConvCapOrFloor (C_isItCapOrFloor, C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	if(C_kRefVal == "DEFAULT")
	{
		kRefValId = ARM_NULL_OBJECT_ID;
	}
	else
	{
		kRefValId = LocalGetNumObjectId (C_kRefVal);
	}

	long retCode;
	long objId;
	CCString prevClass;
	
	CCString curClass = LOCAL_STICKY_CLASS;
	CCString stringId = GetLastCurCellEnvValue ();
	
	if(!stringId)
	{
		retCode = ARMLOCAL_STICKY (LocalGetNumObjectId (C_swapLeg),
							  isItCapOrFloorId,
							  C_strike,
							  C_spreadDates,
							  C_spreadValues,
							  kRefValId,
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
			retCode = ARMLOCAL_STICKY (LocalGetNumObjectId (C_swapLeg),
								  isItCapOrFloorId,
								  C_strike,
								  C_spreadDates,
								  C_spreadValues,
								  kRefValId,
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
			retCode = ARMLOCAL_STICKY (LocalGetNumObjectId (C_swapLeg),
								  isItCapOrFloorId,
								  C_strike,
								  C_spreadDates,
								  C_spreadValues,
								  kRefValId,
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_STICKY" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}



__declspec(dllexport) LPXLOPER WINAPI Local_PXL_STICKY (LPXLOPER XL_swapLeg,
														LPXLOPER XL_isItCapOrFloor,
														LPXLOPER XL_strike,
														LPXLOPER XL_spreadDates,
														LPXLOPER XL_spreadValues,
														LPXLOPER XL_kRefVal)
{
	ADD_LOG("Local_PXL_STICKY ");
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();

	// C variable
	CCString C_swapLeg;
	
	CCString C_isItCapOrFloor;
	long isItCapOrFloorId;

	double C_strike;

	VECTOR<double> C_spreadDates;
	VECTOR<double> C_spreadValues;

	CCString C_kRefVal;
	long kRefValId;

	// error
	static int error;
	static char* reason = "";

	XL_readStrCell(XL_swapLeg,C_swapLeg," ARM_ERR: swap leg id: object expected",C_result);
	XL_readStrCell(XL_isItCapOrFloor,C_isItCapOrFloor," ARM_ERR: cap or floor: string expected",C_result);
	XL_readNumCell(XL_strike,C_strike," ARM_ERR: strike: numeric expected",C_result);
	XL_readNumVector(XL_spreadDates,C_spreadDates," ARM_ERR: spread dates: array of date expected",C_result);
	XL_readNumVector(XL_spreadValues,C_spreadValues," ARM_ERR: spread values: array of numeric expected",C_result);
	XL_readStrCellWD(XL_kRefVal,C_kRefVal,"DEFAULT"," ARM_ERR: strike reference value id: object expected",C_result);
	
	if((isItCapOrFloorId = ARM_ConvCapOrFloor (C_isItCapOrFloor, C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	if(C_kRefVal == "DEFAULT")
	{
		kRefValId = ARM_NULL_OBJECT_ID;
	}
	else
	{
		kRefValId = LocalGetNumObjectId (C_kRefVal);
	}

	long retCode;
	long objId;
	
	CCString curClass = LOCAL_STICKY_CLASS;
	CCString stringId;

	retCode = ARMLOCAL_STICKY (LocalGetNumObjectId (C_swapLeg),
						  isItCapOrFloorId,
						  C_strike,
						  C_spreadDates,
						  C_spreadValues,
						  kRefValId,
						  C_result);

	if ( retCode == ARM_OK )
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_PXL_STICKY" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


__declspec(dllexport) LPXLOPER WINAPI Local_ARM_SPREADOPTION (LPXLOPER XL_startDate,
															  LPXLOPER XL_endDate,
															  LPXLOPER XL_capOrFloor,
															  LPXLOPER XL_strike,
															  LPXLOPER XL_liborType1,
															  LPXLOPER XL_liborType2,
															  LPXLOPER XL_weight1,
															  LPXLOPER XL_weight2,
															  LPXLOPER XL_daycount,
															  LPXLOPER XL_resetFreq,
															  LPXLOPER XL_payFreq,
															  LPXLOPER XL_resetTiming,
															  LPXLOPER XL_payTiming,
															  LPXLOPER XL_currency,
                                                              LPXLOPER XL_resetGap,
															  LPXLOPER XL_intRule,
															  LPXLOPER XL_stubRule,
															  LPXLOPER XL_fixing1,
															  LPXLOPER XL_fixing2,
															  LPXLOPER XL_cptStrikeMethod,
															  LPXLOPER XL_VCtcalibInfo)
{
	ADD_LOG("Local_ARM_SPREADOPTION ");
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

	CCString C_capOrFloor;
	long capOrFloorId;

	double C_strike_double;
	CCString C_strike_str;
	long strike_type;

	CCString C_liborType1;
	long liborType1Id;
	CCString C_liborType2;
	long liborType2Id;

	double C_weight1_double;
	double C_weight1_default = 1.0;
	CCString C_weight1_str;
	bool weight1IsReferenceValue = false;
	long weight1Id;
	double C_weight2_double;
	double C_weight2_default = 1.0;
	CCString C_weight2_str;
	bool weight2IsReferenceValue = false;
	long weight2Id;

	CCString C_daycount;
	long dayCountId;

	CCString C_resetFreq;
	long resetFreqId;
	CCString C_payFreq;
	long payFreqId;

	CCString C_resetTiming;
	long resetTimingId;
	CCString C_payTiming;
	long payTimingId;

	CCString C_ccy;
    long ccyId;

    double C_resetGap;
    double C_resetGap_default = 10000.0;

   	double C_fixing1_double;
   	VECTOR<double> C_fixing1_double_vect;
	CCString C_fixing1_str;
	long fixing1_type;
	VECTOR<double> C_fixing1_default;

   	double C_fixing2_double;
   	VECTOR<double> C_fixing2_double_vect;
	CCString C_fixing2_str;
	long fixing2_type;
	VECTOR<double> C_fixing2_default;
	
	CCString C_intRule;
	double intRuleId;

	CCString C_stubRule;
	double stubRuleId;

	VECTOR<CCString> C_VCtcptStrikeMethod;
	VECTOR<CCString> C_VCtcptStrikeMethod_default;

	double C_cptStrikeMethod;
	double C_cptStrikeMethod_default = 1.0;

	double C_computedFormula;
	double C_computedFormula_default = 1.0;

	VECTOR<double> C_VCtcalibInfo;
	VECTOR<double> C_VCtcalibInfo_default;

	// error
	static int error;
	static char* reason = "";

	XL_readNumCell(XL_startDate,C_startDate," ARM_ERR: start date: date expected",C_result);
	XL_readNumCell(XL_endDate,C_endDate," ARM_ERR: end date: date expected",C_result);
	XL_readStrCell(XL_capOrFloor,C_capOrFloor," ARM_ERR: cap or floor: string expected",C_result);
	XL_readStrOrNumCell(XL_strike,C_strike_str,C_strike_double,strike_type," ARM_ERR: strike: string or numeric expected",C_result);
	XL_readStrCellWD(XL_liborType1,C_liborType1,"LIBOR3M"," ARM_ERR: libor type 1 : string expected",C_result);
	XL_readStrCellWD(XL_liborType2,C_liborType2,"LIBOR3M"," ARM_ERR: libor type 2 : string expected",C_result);
    XL_readStrOrNumCellWD(XL_weight1, C_weight1_str, C_weight1_double, C_weight1_default, weight1Id,
							" ARM_ERR: weight 1 : numeric or ReferenceValue expected",C_result);
	XL_readStrOrNumCellWD(XL_weight2, C_weight2_str, C_weight2_double, C_weight2_default, weight2Id,
							" ARM_ERR: weight 2 : numeric or ReferenceValue expected",C_result);
	XL_readStrCellWD(XL_daycount,C_daycount,"30/360"," ARM_ERR: daycount : string expected",C_result);
	XL_readStrCellWD(XL_resetFreq,C_resetFreq,"A"," ARM_ERR: reset frequency : string expected",C_result);
	XL_readStrCellWD(XL_payFreq,C_payFreq,"A"," ARM_ERR: payment frequency : string expected",C_result);
	XL_readStrCellWD(XL_resetTiming,C_resetTiming,"ADV"," ARM_ERR: reset timing : string expected",C_result);
	XL_readStrCellWD(XL_payTiming,C_payTiming,"ARR"," ARM_ERR: payment timing : string expected",C_result);
	XL_readStrCellWD(XL_currency,C_ccy,"DEFAULT"," ARM_ERR: currency : object expected",C_result);
	XL_readNumCellWD(XL_resetGap,C_resetGap,C_resetGap_default," ARM_ERR: resetGap : numeric expected",C_result);
	XL_readStrCellWD(XL_intRule,C_intRule,"DEFAULT"," ARM_ERR: intRule : string expected",C_result);
	XL_readStrCellWD(XL_stubRule,C_stubRule,"DEFAULT"," ARM_ERR: stubRule : string expected",C_result);

	XL_readStrOrNumVectorWD(XL_fixing1,C_fixing1_str,C_fixing1_double_vect,C_fixing1_default,fixing1_type," ARM_ERR: fixing taux : array of numeric expected",C_result);
	XL_readStrOrNumVectorWD(XL_fixing2,C_fixing2_str,C_fixing2_double_vect,C_fixing2_default,fixing2_type," ARM_ERR: fixing taux : array of numeric expected",C_result);
	// XL_readNumVectorWD(XL_cptStrikeMethod,C_VCtcptStrikeMethod,C_VCtcptStrikeMethod_default," ARM_ERR: cptStrikeMethod : vector of integer(0/1) expected",C_result);

	XL_readStrVectorWD(XL_cptStrikeMethod,C_VCtcptStrikeMethod,C_VCtcptStrikeMethod_default," ARM_ERR: cptStrikeMethod : vector of integer(0/1) expected",XL_TYPE_STRING, C_result);

/*
	int real_size = C_VCtcptStrikeMethod.size ();
	double tmp;
	int j = 0;

	for(int i = 0; i < real_size; i++)
	{
		C_matu.push_back (C_matuRate[j]);
		if((sscanf ((const char*)C_matuRate[j+1], "%lf", &tmp) != 1))
		{
			C_result.setMsg ("ARM_ERR: check your maturities and rates array");
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}
		C_rate.push_back (tmp);
		j += 2;
	}
*/

	if(C_VCtcptStrikeMethod.size() == 0)
	{
		C_cptStrikeMethod = C_cptStrikeMethod_default;
		C_computedFormula = C_computedFormula_default;
	}
	else if(C_VCtcptStrikeMethod.size() == 1)
	{
		// C_cptStrikeMethod = C_VCtcptStrikeMethod[0];
		sscanf ((const char*)C_VCtcptStrikeMethod[0], "%lf", &C_cptStrikeMethod);
		C_computedFormula = C_computedFormula_default;
	}
	else
	{
		// C_cptStrikeMethod = C_VCtcptStrikeMethod[0];
		// C_computedFormula = C_VCtcptStrikeMethod[1];
		sscanf ((const char*)C_VCtcptStrikeMethod[0], "%lf", &C_cptStrikeMethod);
		sscanf ((const char*)C_VCtcptStrikeMethod[1], "%lf", &C_computedFormula);

		// get Sabr calib infos for each index

		double tmp;
		for (int j=0;j<C_VCtcptStrikeMethod.size()-2;j++)
		{
			sscanf ((const char*)C_VCtcptStrikeMethod[j+2], "%lf", &tmp);
			C_VCtcalibInfo.push_back (tmp);
		}
	}

	if ((XL_fixing1->xltype == xltypeMissing) || (XL_fixing1->xltype == xltypeNil))
	{
		fixing1_type = 1;
		C_fixing1_double = ARM_NULL_OBJECT;
    }
	else if ( fixing1_type == XL_TYPE_STRING )
	{
		C_fixing1_double = (double) LocalGetNumObjectId(C_fixing1_str);
		fixing1_type = 1;
	}
	else
	{
	    fixing1_type = 0;
	}

	if ((XL_fixing2->xltype == xltypeMissing) || (XL_fixing2->xltype == xltypeNil))
	{
		fixing2_type = 1;
		C_fixing2_double = ARM_NULL_OBJECT;
    }
	else if ( fixing2_type == XL_TYPE_STRING )
	{
		C_fixing2_double = (double) LocalGetNumObjectId(C_fixing2_str);
		fixing2_type = 1;
	}
	else
	{
	    fixing2_type = 0;
	}

	if((capOrFloorId = ARM_ConvCapOrFloor (C_capOrFloor, C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	liborType1Id = ARM_ConvIrType (C_liborType1);

	liborType2Id = ARM_ConvIrType (C_liborType2);

	dayCountId = ARM_ConvDayCount (C_daycount);
	
	if((resetFreqId = ARM_ConvFrequency (C_resetFreq, C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	if((payFreqId = ARM_ConvFrequency (C_payFreq, C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	resetTimingId = ARM_ConvPayResetRule (C_resetTiming);

	payTimingId = ARM_ConvPayResetRule (C_payTiming);

	if(C_ccy == "DEFAULT")
	{
		ccyId = ARM_NULL_OBJECT;
	}
	else
	{
		ccyId = LocalGetNumObjectId (C_ccy);
	}

	if ( strike_type == XL_TYPE_STRING )
	{
	   C_strike_double = (double) LocalGetNumObjectId(C_strike_str);

	   strike_type = 1;
	}
	else
	{
	   strike_type = 0;
	}

	if(C_intRule == "DEFAULT")
	{
		C_intRule = "ADJ";
	}
	intRuleId = ARM_ConvIntRule(C_intRule);

	if(C_stubRule == "DEFAULT")
	{
		C_stubRule = "SS";
	}
	stubRuleId = ARM_ConvStubRule(C_stubRule);

	if(weight1Id == XL_TYPE_STRING)
	{
		weight1IsReferenceValue = true;
		weight1Id = (long) LocalGetNumObjectId(C_weight1_str);
	}

	if(weight2Id == XL_TYPE_STRING)
	{
		weight2IsReferenceValue = true;
		weight2Id = (long) LocalGetNumObjectId(C_weight2_str);
	}

	long retCode;
	long objId;
	CCString prevClass;
	
	CCString curClass = LOCAL_SPREAD_OPTION_CLASS;
	CCString stringId = GetLastCurCellEnvValue ();

	if(!stringId)
	{
		retCode = ARMLOCAL_SPREADOPTION (C_startDate,
										 C_endDate,
										 capOrFloorId,
										 (long) strike_type,
										 C_strike_double,
										 (long)liborType1Id,
										 (long)liborType2Id,
										 C_weight1_double,
										 weight1IsReferenceValue,
										 weight1Id,
										 C_weight2_double,
										 weight1IsReferenceValue,
										 weight2Id,
										 dayCountId,
										 resetFreqId,
										 payFreqId,
										 resetTimingId,
										 payTimingId,
										 ccyId,
                                         (long) C_resetGap,
										 (long) fixing1_type,
										 C_fixing1_double,
										 C_fixing1_double_vect,
										 (long) fixing1_type,
										 C_fixing2_double,
										 C_fixing2_double_vect,
										 (long) intRuleId,
										 (long) stubRuleId,
										 (int) C_cptStrikeMethod,
										 (int) C_computedFormula,
										 C_VCtcalibInfo,
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
			retCode = ARMLOCAL_SPREADOPTION (C_startDate,
											 C_endDate,
											 capOrFloorId,
											 (long) strike_type,
											 C_strike_double,
											 (long)liborType1Id,
											 (long)liborType2Id,
											 C_weight1_double,
											 weight1IsReferenceValue,
											 weight1Id,
											 C_weight2_double,
											 weight1IsReferenceValue,
											 weight2Id,
											 dayCountId,
											 resetFreqId,
											 payFreqId,
											 resetTimingId,
											 payTimingId,
											 ccyId,
                                             long(C_resetGap),
											 (long) fixing1_type,
											 C_fixing1_double,
											 C_fixing1_double_vect,
											 (long) fixing1_type,
											 C_fixing2_double,
											 C_fixing2_double_vect,
											 (long) intRuleId,
											 (long) stubRuleId,
											 (int) C_cptStrikeMethod,
											 (int) C_computedFormula,
											 C_VCtcalibInfo,
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
			retCode = ARMLOCAL_SPREADOPTION (C_startDate,
											 C_endDate,
											 capOrFloorId,
											 (long) strike_type,
											 C_strike_double,
											 (long)liborType1Id,
											 (long)liborType2Id,
											 C_weight1_double,
											 weight1IsReferenceValue,
											 weight1Id,
											 C_weight2_double,
											 weight2IsReferenceValue,
											 weight2Id,
											 dayCountId,
											 resetFreqId,
											 payFreqId,
											 resetTimingId,
											 payTimingId,
											 ccyId,
                                             long(C_resetGap),
											 (long) fixing1_type,
											 C_fixing1_double,
											 C_fixing1_double_vect,
											 (long) fixing1_type,
											 C_fixing2_double,
											 C_fixing2_double_vect,
											 (long) intRuleId,
											 (long) stubRuleId,
											 (int) C_cptStrikeMethod,
											 (int) C_computedFormula,
											 C_VCtcalibInfo,
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_ARM_SPREADOPTION" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


__declspec(dllexport) LPXLOPER WINAPI Local_PXL_ARM_SPREADOPTION (LPXLOPER XL_startDate,
																  LPXLOPER XL_endDate,
																  LPXLOPER XL_capOrFloor,
																  LPXLOPER XL_strike,
																  LPXLOPER XL_liborType1,
																  LPXLOPER XL_liborType2,
																  LPXLOPER XL_weight1,
																  LPXLOPER XL_weight2,
																  LPXLOPER XL_daycount,
																  LPXLOPER XL_resetFreq,
																  LPXLOPER XL_payFreq,
																  LPXLOPER XL_resetTiming,
																  LPXLOPER XL_payTiming,
																  LPXLOPER XL_currency,
                                                                  LPXLOPER XL_resetGap,
																  LPXLOPER XL_intRule,
																  LPXLOPER XL_stubRule,
																  LPXLOPER XL_fixing1,
																  LPXLOPER XL_fixing2,
																  LPXLOPER XL_cptStrikeMethod)
{
	ADD_LOG("Local_PXL_ARM_SPREADOPTION ");
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

	CCString C_capOrFloor;
	long capOrFloorId;

	double C_strike_double;
	CCString C_strike_str;
	long strike_type;


	CCString C_liborType1;
	long liborType1Id;
	CCString C_liborType2;
	long liborType2Id;

	double C_weight1_double;
	double C_weight1_default = 1.0;
	CCString C_weight1_str;
	bool weight1IsReferenceValue = false;
	long weight1Id;
	double C_weight2_double;
	double C_weight2_default = 1.0;
	CCString C_weight2_str;
	bool weight2IsReferenceValue = false;
	long weight2Id;

	CCString C_daycount;
	long dayCountId;

	CCString C_resetFreq;
	long resetFreqId;
	CCString C_payFreq;
	long payFreqId;

	CCString C_resetTiming;
	long resetTimingId;
	CCString C_payTiming;
	long payTimingId;

	CCString C_ccy;
    long ccyId;

    double C_resetGap;
    double C_resetGap_default = 10000.0;

   	double C_fixing1_double;
   	VECTOR<double> C_fixing1_double_vect;
	CCString C_fixing1_str;
	long fixing1_type;
	VECTOR<double> C_fixing1_default;

   	double C_fixing2_double;
   	VECTOR<double> C_fixing2_double_vect;
	CCString C_fixing2_str;
	long fixing2_type;
	VECTOR<double> C_fixing2_default;

	CCString C_intRule;
    long intRuleId;

	CCString C_stubRule;
    long stubRuleId;

	VECTOR<double> C_VCtcptStrikeMethod;
	VECTOR<double> C_VCtcptStrikeMethod_default;

	double C_cptStrikeMethod;
	double C_cptStrikeMethod_default = 1.0;

	double C_computedFormula;
	double C_computedFormula_default = 1.0;

	VECTOR<double> C_VCtcalibInfo;
	VECTOR<double> C_VCtcalibInfo_default;
	// error
	static int error;
	static char* reason = "";

	XL_readNumCell(XL_startDate,C_startDate," ARM_ERR: start date: date expected",C_result);
	XL_readNumCell(XL_endDate,C_endDate," ARM_ERR: end date: date expected",C_result);
	XL_readStrCell(XL_capOrFloor,C_capOrFloor," ARM_ERR: cap or floor: string expected",C_result);
	XL_readStrOrNumCell(XL_strike,C_strike_str,C_strike_double,strike_type," ARM_ERR: strike: string or numeric expected",C_result);
	XL_readStrCellWD(XL_liborType1,C_liborType1,"LIBOR3M"," ARM_ERR: libor type 1 : string expected",C_result);
	XL_readStrCellWD(XL_liborType2,C_liborType2,"LIBOR3M"," ARM_ERR: libor type 2 : string expected",C_result);
	XL_readStrOrNumCellWD(XL_weight1, C_weight1_str, C_weight1_double, C_weight1_default, weight1Id,
							" ARM_ERR: weight 1 : numeric or ReferenceValue expected",C_result);
	XL_readStrOrNumCellWD(XL_weight2, C_weight2_str, C_weight2_double, C_weight2_default, weight2Id,
							" ARM_ERR: weight 2 : numeric or ReferenceValue expected",C_result);	XL_readStrCellWD(XL_daycount,C_daycount,"30/360"," ARM_ERR: daycount : string expected",C_result);
	XL_readStrCellWD(XL_resetFreq,C_resetFreq,"A"," ARM_ERR: reset frequency : string expected",C_result);
	XL_readStrCellWD(XL_payFreq,C_payFreq,"A"," ARM_ERR: payment frequency : string expected",C_result);
	XL_readStrCellWD(XL_resetTiming,C_resetTiming,"ADV"," ARM_ERR: reset timing : string expected",C_result);
	XL_readStrCellWD(XL_payTiming,C_payTiming,"ARR"," ARM_ERR: payment timing : string expected",C_result);
	XL_readStrCellWD(XL_currency,C_ccy,"DEFAULT"," ARM_ERR: currency : object expected",C_result);
	XL_readNumCellWD(XL_resetGap,C_resetGap,C_resetGap_default," ARM_ERR: resetGap : numeric expected",C_result);
	XL_readStrCellWD(XL_intRule,C_intRule,"DEFAULT"," ARM_ERR: intRule : string expected",C_result);
	XL_readStrCellWD(XL_stubRule,C_stubRule,"DEFAULT"," ARM_ERR: stubRule : string expected",C_result);
	
	XL_readStrOrNumVectorWD(XL_fixing1,C_fixing1_str,C_fixing1_double_vect,C_fixing1_default,fixing1_type," ARM_ERR: fixing taux : array of numeric expected",C_result);
	XL_readStrOrNumVectorWD(XL_fixing2,C_fixing2_str,C_fixing2_double_vect,C_fixing2_default,fixing2_type," ARM_ERR: fixing taux : array of numeric expected",C_result);
	XL_readNumVectorWD(XL_cptStrikeMethod,C_VCtcptStrikeMethod,C_VCtcptStrikeMethod_default," ARM_ERR: cptStrikeMethod : vector of integer(0/1) expected",C_result);

	if(C_VCtcptStrikeMethod.size() == 0)
	{
		C_cptStrikeMethod = C_cptStrikeMethod_default;
		C_computedFormula = C_computedFormula_default;
	}
	else if(C_VCtcptStrikeMethod.size() == 1)
	{
		C_cptStrikeMethod = C_VCtcptStrikeMethod[0];
		C_computedFormula = C_computedFormula_default;
	}
	else
	{
		C_cptStrikeMethod = C_VCtcptStrikeMethod[0];
		C_computedFormula = C_VCtcptStrikeMethod[1];
	}
	
	if ((XL_fixing1->xltype == xltypeMissing) || (XL_fixing1->xltype == xltypeNil))
	{
		fixing1_type = 1;
		C_fixing1_double = ARM_NULL_OBJECT;
    }
	else if ( fixing1_type == XL_TYPE_STRING )
	{
		C_fixing1_double = (double) LocalGetNumObjectId(C_fixing1_str);
		fixing1_type = 1;
	}
	else
	{
	    fixing1_type = 0;
	}

	if ((XL_fixing2->xltype == xltypeMissing) || (XL_fixing2->xltype == xltypeNil))
	{
		fixing2_type = 1;
		C_fixing2_double = ARM_NULL_OBJECT;
    }
	else if ( fixing2_type == XL_TYPE_STRING )
	{
		C_fixing2_double = (double) LocalGetNumObjectId(C_fixing2_str);
		fixing2_type = 1;
	}
	else
	{
	    fixing2_type = 0;
	}

	if((capOrFloorId = ARM_ConvCapOrFloor (C_capOrFloor, C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	liborType1Id = ARM_ConvIrType (C_liborType1);

	liborType2Id = ARM_ConvIrType (C_liborType2);

	dayCountId = ARM_ConvDayCount (C_daycount);

	if((resetFreqId = ARM_ConvFrequency (C_resetFreq, C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	if((payFreqId = ARM_ConvFrequency (C_payFreq, C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	resetTimingId = ARM_ConvPayResetRule (C_resetTiming);

	payTimingId = ARM_ConvPayResetRule (C_payTiming);

	if(C_ccy == "DEFAULT")
	{
		ccyId = ARM_NULL_OBJECT;
	}
	else
	{
		ccyId = LocalGetNumObjectId (C_ccy);
	}

	if ( strike_type == XL_TYPE_STRING )
	{
	   C_strike_double = (double) LocalGetNumObjectId(C_strike_str);

	   strike_type = 1;
	}
	else
	{
	   strike_type = 0;
	}

	if(C_intRule == "DEFAULT")
	{
		C_intRule = "ADJ";
	}
	intRuleId = ARM_ConvIntRule(C_intRule);
	
	if(C_stubRule == "DEFAULT")
	{
		C_stubRule = "SS";
	}
	stubRuleId = ARM_ConvStubRule(C_stubRule);

	if(weight1Id == XL_TYPE_STRING)
	{
		weight1IsReferenceValue = true;
		weight1Id = (long) LocalGetNumObjectId(C_weight1_str);
	}

	if(weight2Id == XL_TYPE_STRING)
	{
		weight2IsReferenceValue = true;
		weight2Id = (long) LocalGetNumObjectId(C_weight2_str);
	}

	long retCode;
	long objId;
	
	CCString curClass = LOCAL_SPREAD_OPTION_CLASS;
	CCString stringId;

	retCode = ARMLOCAL_SPREADOPTION (C_startDate,
									 C_endDate,
									 capOrFloorId,
									 (long) strike_type,
									 C_strike_double,
									 (long)liborType1Id,
									 (long)liborType2Id,
									 C_weight1_double,
									 weight1IsReferenceValue,
									 weight1Id,
									 C_weight2_double,
									 weight2IsReferenceValue,
									 weight2Id,
									 dayCountId,
									 resetFreqId,
									 payFreqId,
									 resetTimingId,
									 payTimingId,
									 ccyId,
                                     long(C_resetGap),
									 (long) fixing1_type,
									 C_fixing1_double,
									 C_fixing1_double_vect,
									 (long) fixing1_type,
									 C_fixing2_double,
									 C_fixing2_double_vect,
									 intRuleId,
									 stubRuleId,
									 (int) C_cptStrikeMethod,
									 (int) C_computedFormula,
									 C_VCtcalibInfo,
									 C_result);

	if ( retCode == ARM_OK )
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_PXL_ARM_SPREADOPTION" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}

__declspec(dllexport) LPXLOPER WINAPI Local_ARM_QUANTOSPREADOPTION(		LPXLOPER XL_startDate,
																		LPXLOPER XL_endDate,
 																		LPXLOPER XL_capOrFloor,
																		LPXLOPER XL_strikes,
																		LPXLOPER XL_liborIdx1,
																		LPXLOPER XL_liborIdx2,
																		LPXLOPER XL_leg1Weight,
																		LPXLOPER XL_leg2Weight,
																		LPXLOPER XL_leg1Fixings,
																		LPXLOPER XL_leg2Fixings,
																		LPXLOPER XL_leg1Spread,
																		LPXLOPER XL_leg2Spread,
																		LPXLOPER XL_modelParams,
																		LPXLOPER XL_prodParams,
																		LPXLOPER XL_notional)
{
	ADD_LOG("Local_ARM_QUANTOSPREADOPTION");
	//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();

	double C_startDate;
	double C_endDate;

	CCString C_capOrFloor;
	long capOrFloorId;

	CCString C_strikes_str;
	double C_strike_double;
	long strike_type;

	CCString C_liborIdx1;
	long liborIdx1Id;
	CCString C_liborIdx2;
	long liborIdx2Id;
	
	double C_Idx1weight;
	double C_Idx1weight_default = 1.0;

	double C_Idx2weight;
	double C_Idx2weight_default = 1.0;

	CCString C_Idx1fixings_str;
	CCString C_Idx1fixings_str_default("DEFAULT");
	long Idx1fixingsId = ARM_NULL_OBJECT;
	CCString C_Idx2fixings_str;
	CCString C_Idx2fixings_str_default("DEFAULT");
	long Idx2fixingsId = ARM_NULL_OBJECT;
	
	double C_Idx1spread;
	double C_Idx1spread_default = 0.0;

	double C_Idx2spread;
	double C_Idx2spread_default = 0.0;

	VECTOR<double> C_modelParams;
	VECTOR<double> C_modelParams_default(0);
	VECTOR<double> modelParams(2);
	
	modelParams[0] = 1.0;
	modelParams[1] = 1.0;
	
	double C_cptStrikeMethod;
	double C_computedFormula;

	// Product parameters
	VECTOR<CCString> C_prodParams;
	VECTOR<CCString> C_prodParams_default(0);

	VECTOR<CCString> prodParams(13);

	prodParams[0] = "DEFAULT";
	prodParams[1] = "30/360";
	prodParams[2] = "A";
	prodParams[3] = "A";
	prodParams[4] = "ADV";
	prodParams[5] = "ARR";
	prodParams[6] = "ADJ";
	prodParams[7] = "SS";
	prodParams[8] = "2";
	prodParams[9] = ARM_NULL_OBJECT;
	prodParams[10] = ARM_NULL_OBJECT;
	prodParams[11] = "MF";
	prodParams[12] = "NULL";

	CCString C_ccy;
	long ccyId;
	
	CCString C_daycount;
	long dayCountId;

	CCString C_resetFreq;
	long resetFreqId;

	CCString C_payFreq;
	long payFreqId;

	CCString C_resetTiming;
	long resetTimingId;
	
	CCString C_payTiming;
	long payTimingId;

	CCString C_intRule;
    long intRuleId;
	
	CCString C_stubRule;
    long stubRuleId;

    double C_resetGap;
	
	CCString C_resetCal;
	CCString C_payCal;
	CCString C_fwdRule;
	double fwdRuleId;

	CCString C_refDate;
	
	long notionalId = ARM_NULL_OBJECT;
	CCString C_notional_str;
	CCString C_notional_default("DEFAULT");

	// error
	static int error;
	static char* reason = "";

	XL_readNumCell(XL_startDate,C_startDate," ARM_ERR: start date: date expected",C_result);
	XL_readNumCell(XL_endDate,C_endDate," ARM_ERR: end date: date expected",C_result);
	XL_readStrCell(XL_capOrFloor,C_capOrFloor," ARM_ERR: cap or floor: string expected",C_result);
	XL_readStrOrNumCell(XL_strikes,C_strikes_str,C_strike_double,strike_type," ARM_ERR: strikes vector : array of numeric expected",C_result);
	XL_readStrCellWD(XL_liborIdx1,C_liborIdx1,"DEFAULT"," ARM_ERR: leg1 libor index : string expected",C_result);
	XL_readStrCellWD(XL_liborIdx2,C_liborIdx2,"DEFAULT"," ARM_ERR: leg2 libor index : string expected",C_result);
	XL_readNumCellWD(XL_leg1Weight,C_Idx1weight,C_Idx1weight_default," ARM_ERR: leg1 weight : numeric expected",C_result);
	XL_readNumCellWD(XL_leg2Weight,C_Idx2weight,C_Idx2weight_default," ARM_ERR: leg2 weight : numeric expected",C_result);
	XL_readStrCellWD(XL_leg1Fixings,C_Idx1fixings_str,C_Idx1fixings_str_default," ARM_ERR: leg1 fixings vector : object expected",C_result);
	XL_readStrCellWD(XL_leg2Fixings,C_Idx2fixings_str,C_Idx2fixings_str_default," ARM_ERR: leg2 fixings vector : object expected",C_result);
	XL_readNumCellWD(XL_leg1Spread,C_Idx1spread, C_Idx1spread_default, " ARM_ERR: leg1 spread : numeric expected",C_result);
	XL_readNumCellWD(XL_leg2Spread,C_Idx2spread,C_Idx2spread_default, " ARM_ERR: leg2 spread : numeric expected",C_result);
	XL_readNumVectorWD(XL_modelParams,C_modelParams,C_modelParams_default," ARM_ERR: model parameters vector : array of numeric expected",C_result);
	XL_readStrVectorWD(XL_prodParams,C_prodParams,C_prodParams_default," ARM_ERR: product parameters vector : array of numeric expected",DOUBLE_TYPE,C_result);
	XL_readStrCellWD(XL_notional,C_notional_str,C_notional_default," ARM_ERR: notional vector : object expected",C_result);

	if((capOrFloorId = ARM_ConvCapOrFloor (C_capOrFloor, C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}
		
	if ( strike_type == XL_TYPE_STRING )
	{
	   C_strike_double = (double) LocalGetNumObjectId(C_strikes_str);

	   strike_type = 1;
	}
	else
	{
	   strike_type = 0;
	}

	liborIdx1Id = strcmp(C_liborIdx1,"DEFAULT") ? (double) LocalGetNumObjectId(C_liborIdx1) : -1.0;	
	liborIdx2Id = strcmp(C_liborIdx2,"DEFAULT") ? (double) LocalGetNumObjectId(C_liborIdx2) : -1.0;	
	
	Idx1fixingsId = strcmp(C_Idx1fixings_str,"DEFAULT") ? (long) LocalGetNumObjectId(C_Idx1fixings_str) : -1.0;
	Idx2fixingsId = strcmp(C_Idx2fixings_str,"DEFAULT") ? (long) LocalGetNumObjectId(C_Idx2fixings_str) : -1.0;
	
	int i;
	for (i = 0; i < C_modelParams.size(); i++)
	{
		modelParams[i]  = C_modelParams[i];
	}
	
	C_cptStrikeMethod = modelParams[0];
	C_computedFormula = modelParams[1];

	for (i = 0; i < C_prodParams.size(); i++)
	{
		if (strcmp(C_prodParams[i],"DEFAULT"))
			prodParams[i] = C_prodParams[i];
	}

	C_ccy = prodParams[0];
	
	if(C_ccy == "DEFAULT")
	{
		ccyId = ARM_NULL_OBJECT;
	}
	else
	{
		ccyId = LocalGetNumObjectId (C_ccy);
	}
	
	C_daycount = prodParams[1];
	dayCountId = ARM_ConvDayCount (C_daycount);
	
	C_resetFreq = prodParams[2];
	if((resetFreqId = ARM_ConvFrequency (C_resetFreq, C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}
	
	C_payFreq = prodParams[3];
	if((payFreqId = ARM_ConvFrequency (C_payFreq, C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	C_resetTiming = prodParams[4];
	resetTimingId = ARM_ConvPayResetRule (C_resetTiming);
	
	C_payTiming = prodParams[5];
	payTimingId = ARM_ConvPayResetRule (C_payTiming);

	C_intRule  = prodParams[6];
    intRuleId = ARM_ConvIntRule(C_intRule);
	
	C_stubRule = prodParams[7];
	stubRuleId = ARM_ConvStubRule(C_stubRule);

    C_resetGap = atof(prodParams[8]);
	
	C_resetCal = prodParams[9];

	C_payCal = prodParams[10];

	C_fwdRule = prodParams[11];
	fwdRuleId = ARM_ConvFwdRule(C_fwdRule);

	C_refDate = prodParams[12];

	notionalId = strcmp(C_notional_str,"DEFAULT") ? (long) LocalGetNumObjectId(C_notional_str) : -1.0;

	long retCode;
	long objId;
	
	CCString prevClass;
	CCString curClass = LOCAL_SPREAD_OPTION_CLASS;
	CCString stringId = GetLastCurCellEnvValue();
	
	if(!stringId)
	{
		retCode = ARMLOCAL_QUANTOSPREADOPTION (		C_startDate,
													C_endDate,
													capOrFloorId,
													(long) strike_type,
													C_strike_double,
													(long) liborIdx1Id,
													(long) liborIdx2Id,
													C_Idx1weight,
													C_Idx2weight,
													Idx1fixingsId,
													Idx2fixingsId,
													C_Idx1spread,
													C_Idx2spread,
													C_cptStrikeMethod,
													(int) C_computedFormula,
													ccyId,
													dayCountId,
													resetFreqId,
													payFreqId,
													resetTimingId,
													payTimingId,
													intRuleId,
													stubRuleId,
													C_resetGap,
													C_resetCal,
													C_payCal,
													fwdRuleId,
													C_refDate,
													notionalId,
													C_result);
		if (retCode == ARM_OK)
		{
			objId = C_result.getLong();

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
			retCode = ARMLOCAL_QUANTOSPREADOPTION(C_startDate,
													C_endDate,
													capOrFloorId,
													(long) strike_type,
													C_strike_double,
													(long) liborIdx1Id,
													(long) liborIdx2Id,
													C_Idx1weight,
													C_Idx2weight,
													Idx1fixingsId,
													Idx2fixingsId,
													C_Idx1spread,
													C_Idx2spread,
													C_cptStrikeMethod,
													(int) C_computedFormula,
													ccyId,
													dayCountId,
													resetFreqId,
													payFreqId,
													resetTimingId,
													payTimingId,
													intRuleId,
													stubRuleId,
													C_resetGap,
													C_resetCal,
													C_payCal,
													fwdRuleId,
													C_refDate,
													notionalId,
													C_result,
													objId);
			if (retCode == ARM_OK)
			{
				objId = C_result.getLong();

				stringId = LocalMakeObjectId (objId, curClass);
			}
		}
		else
		{
			FreeCurCellContent ();
			retCode = ARMLOCAL_QUANTOSPREADOPTION (	C_startDate,
													C_endDate,
													capOrFloorId,
													(long) strike_type,
													C_strike_double,
													(long) liborIdx1Id,
													(long) liborIdx2Id,
													C_Idx1weight,
													C_Idx2weight,
													Idx1fixingsId,
													Idx2fixingsId,
													C_Idx1spread,
													C_Idx2spread,
													C_cptStrikeMethod,
													(int) C_computedFormula,
													ccyId,
													dayCountId,
													resetFreqId,
													payFreqId,
													resetTimingId,
													payTimingId,
													intRuleId,
													stubRuleId,
													C_resetGap,
													C_resetCal,
													C_payCal,
													fwdRuleId,
													C_refDate,
													notionalId,
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_ARM_QUANTOSPREADOPTION"  )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}

__declspec(dllexport) LPXLOPER WINAPI Local_ARM_SPREADOPTIONWithLegs (LPXLOPER XL_firstLeg,
																	  LPXLOPER XL_secondLeg,
																	  LPXLOPER XL_capOrFloor,
																	  LPXLOPER XL_strike,
																	  LPXLOPER XL_weight1,
																	  LPXLOPER XL_weight2,
																	  LPXLOPER XL_fixing1,
																	  LPXLOPER XL_fixing2,
																	  LPXLOPER XL_cptStrikeMethod)
{
	ADD_LOG("Local_ARM_SPREADOPTIONWithLegs ");
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();

	CCString C_firstLeg;
	long firstLegId;

	CCString C_secondLeg;
	long secondLegId;

	CCString C_capOrFloor;
	long capOrFloorId;

	double C_strike_double;
	CCString C_strike_str;
	long strike_type;

	double C_weight1_double;
	double C_weight1_default = 1.0;
	CCString C_weight1_str;
	bool weight1IsReferenceValue = false;
	long weight1Id;
	double C_weight2_double;
	double C_weight2_default = 1.0;
	CCString C_weight2_str;
	bool weight2IsReferenceValue = false;
	long weight2Id;

   	double C_fixing1_double;
   	VECTOR<double> C_fixing1_double_vect;
	CCString C_fixing1_str;
	long fixing1_type;
	VECTOR<double> C_fixing1_default;

   	double C_fixing2_double;
   	VECTOR<double> C_fixing2_double_vect;
	CCString C_fixing2_str;
	long fixing2_type;
	VECTOR<double> C_fixing2_default;

	VECTOR<double> C_VCtcptStrikeMethod;
	VECTOR<double> C_VCtcptStrikeMethod_default;

	double C_cptStrikeMethod;
	double C_cptStrikeMethod_default = 1.0;

	double C_computedFormula;
	double C_computedFormula_default = 1.0;

	// error
	static int error;
	static char* reason = "";

	XL_readStrCell(XL_firstLeg,C_firstLeg," ARM_ERR: first leg id: object expected",C_result);
	XL_readStrCell(XL_secondLeg,C_secondLeg," ARM_ERR: secondleg id: object expected",C_result);

	XL_readStrCell(XL_capOrFloor,C_capOrFloor," ARM_ERR: cap or floor: string expected",C_result);
	XL_readStrOrNumCell(XL_strike,C_strike_str,C_strike_double,strike_type," ARM_ERR: strike: string or numeric expected",C_result);
	XL_readStrOrNumCellWD(XL_weight1, C_weight1_str, C_weight1_double, C_weight1_default, weight1Id,
							" ARM_ERR: weight 1 : numeric or ReferenceValue expected",C_result);
	XL_readStrOrNumCellWD(XL_weight2, C_weight2_str, C_weight2_double, C_weight2_default, weight2Id,
							" ARM_ERR: weight 2 : numeric or ReferenceValue expected",C_result);

	XL_readStrOrNumVectorWD(XL_fixing1,C_fixing1_str,C_fixing1_double_vect,C_fixing1_default,fixing1_type," ARM_ERR: fixing taux : array of numeric expected",C_result);
	XL_readStrOrNumVectorWD(XL_fixing2,C_fixing2_str,C_fixing2_double_vect,C_fixing2_default,fixing2_type," ARM_ERR: fixing taux : array of numeric expected",C_result);

	XL_readNumVectorWD(XL_cptStrikeMethod,C_VCtcptStrikeMethod,C_VCtcptStrikeMethod_default," ARM_ERR: cptStrikeMethod : vector of integer(0/1) expected",C_result);

	if(C_VCtcptStrikeMethod.size() == 0)
	{
		C_cptStrikeMethod = C_cptStrikeMethod_default;
		C_computedFormula = C_computedFormula_default;
	}
	else if(C_VCtcptStrikeMethod.size() == 1)
	{
		C_cptStrikeMethod = C_VCtcptStrikeMethod[0];
		C_computedFormula = C_computedFormula_default;
	}
	else
	{
		C_cptStrikeMethod = C_VCtcptStrikeMethod[0];
		C_computedFormula = C_VCtcptStrikeMethod[1];
	}

	firstLegId = LocalGetNumObjectId (C_firstLeg);
	secondLegId = LocalGetNumObjectId (C_secondLeg);

	if ((XL_fixing1->xltype == xltypeMissing) || (XL_fixing1->xltype == xltypeNil))
	{
		fixing1_type = 1;
		C_fixing1_double = ARM_NULL_OBJECT;
    }
	else if ( fixing1_type == XL_TYPE_STRING )
	{
		C_fixing1_double = (double) LocalGetNumObjectId(C_fixing1_str);
		fixing1_type = 1;
	}
	else
	{
	    fixing1_type = 0;
	}

	if ((XL_fixing2->xltype == xltypeMissing) || (XL_fixing2->xltype == xltypeNil))
	{
		fixing2_type = 1;
		C_fixing2_double = ARM_NULL_OBJECT;
    }
	else if ( fixing2_type == XL_TYPE_STRING )
	{
		C_fixing2_double = (double) LocalGetNumObjectId(C_fixing2_str);
		fixing2_type = 1;
	}
	else
	{
	    fixing2_type = 0;
	}

	if((capOrFloorId = ARM_ConvCapOrFloor (C_capOrFloor, C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	if ( strike_type == XL_TYPE_STRING )
	{
	   C_strike_double = (double) LocalGetNumObjectId(C_strike_str);

	   strike_type = 1;
	}
	else
	{
	   strike_type = 0;
	}

	if(weight1Id == XL_TYPE_STRING)
	{
		weight1IsReferenceValue = true;
		weight1Id = (long) LocalGetNumObjectId(C_weight1_str);
	}

	if(weight2Id == XL_TYPE_STRING)
	{
		weight2IsReferenceValue = true;
		weight2Id = (long) LocalGetNumObjectId(C_weight2_str);
	}

	long retCode;
	long objId;
	CCString prevClass;
	
	CCString curClass = LOCAL_SPREAD_OPTION_CLASS;
	CCString stringId = GetLastCurCellEnvValue ();

	if(!stringId)
	{
		retCode = ARMLOCAL_SPREADOPTIONWithLegs(firstLegId,
												secondLegId,
												capOrFloorId,
												(long) strike_type,
												C_strike_double,
												C_weight1_double,
												weight1IsReferenceValue,
												weight1Id,
												C_weight2_double,
												weight1IsReferenceValue,
												weight2Id,
												(long) fixing1_type,
												C_fixing1_double,
												C_fixing1_double_vect,
												(long) fixing1_type,
												C_fixing2_double,
												C_fixing2_double_vect,
												(int) C_cptStrikeMethod,
												(int) C_computedFormula,
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
			retCode = ARMLOCAL_SPREADOPTIONWithLegs (firstLegId,
													 secondLegId,
													 capOrFloorId,
													 (long) strike_type,
													 C_strike_double,
													 C_weight1_double,
													 weight1IsReferenceValue,
													 weight1Id,
													 C_weight2_double,
													 weight1IsReferenceValue,
													 weight2Id,
													 (long) fixing1_type,
													 C_fixing1_double,
													 C_fixing1_double_vect,
													 (long) fixing1_type,
													 C_fixing2_double,
													 C_fixing2_double_vect,
													 (int) C_cptStrikeMethod,
													 (int) C_computedFormula,
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
			retCode = ARMLOCAL_SPREADOPTIONWithLegs (firstLegId,
													 secondLegId,
													 capOrFloorId,
													 (long) strike_type,
													 C_strike_double,
													 C_weight1_double,
													 weight1IsReferenceValue,
													 weight1Id,
													 C_weight2_double,
													 weight1IsReferenceValue,
													 weight2Id,													 
													 (long) fixing1_type,
													 C_fixing1_double,
													 C_fixing1_double_vect,
													 (long) fixing1_type,
													 C_fixing2_double,
													 C_fixing2_double_vect,
													 (int) C_cptStrikeMethod,
													 (int) C_computedFormula,
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_ARM_SPREADOPTIONWithLegs" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


__declspec(dllexport) LPXLOPER WINAPI Local_PXL_ARM_SPREADOPTIONWithLegs (LPXLOPER XL_firstLeg,
																		  LPXLOPER XL_secondLeg,
																		  LPXLOPER XL_capOrFloor,
																		  LPXLOPER XL_strike,
																		  LPXLOPER XL_weight1,
																		  LPXLOPER XL_weight2,
																		  LPXLOPER XL_fixing1,
																		  LPXLOPER XL_fixing2,
																		  LPXLOPER XL_cptStrikeMethod)
{
	ADD_LOG("Local_PXL_ARM_SPREADOPTIONWithLegs ");
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();

	CCString C_firstLeg;
	long firstLegId;

	CCString C_secondLeg;
	long secondLegId;

	CCString C_capOrFloor;
	long capOrFloorId;

	double C_strike_double;
	CCString C_strike_str;
	long strike_type;

	double C_weight1_double;
	double C_weight1_default = 1.0;
	CCString C_weight1_str;
	bool weight1IsReferenceValue = false;
	long weight1Id;
	double C_weight2_double;
	double C_weight2_default = 1.0;
	CCString C_weight2_str;
	bool weight2IsReferenceValue = false;
	long weight2Id;

   	double C_fixing1_double;
   	VECTOR<double> C_fixing1_double_vect;
	CCString C_fixing1_str;
	long fixing1_type;
	VECTOR<double> C_fixing1_default;

   	double C_fixing2_double;
   	VECTOR<double> C_fixing2_double_vect;
	CCString C_fixing2_str;
	long fixing2_type;
	VECTOR<double> C_fixing2_default;
	
	VECTOR<double> C_VCtcptStrikeMethod;
	VECTOR<double> C_VCtcptStrikeMethod_default;

	double C_cptStrikeMethod;
	double C_cptStrikeMethod_default = 1.0;

	double C_computedFormula;
	double C_computedFormula_default = 1.0;

	// error
	static int error;
	static char* reason = "";

	XL_readStrCell(XL_firstLeg,C_firstLeg," ARM_ERR: first leg id: object expected",C_result);
	XL_readStrCell(XL_secondLeg,C_secondLeg," ARM_ERR: secondleg id: object expected",C_result);

	XL_readStrCell(XL_capOrFloor,C_capOrFloor," ARM_ERR: cap or floor: string expected",C_result);
	XL_readStrOrNumCell(XL_strike,C_strike_str,C_strike_double,strike_type," ARM_ERR: strike: string or numeric expected",C_result);
	XL_readStrOrNumCellWD(XL_weight1, C_weight1_str, C_weight1_double, C_weight1_default, weight1Id,
							" ARM_ERR: weight 1 : numeric or ReferenceValue expected",C_result);
	XL_readStrOrNumCellWD(XL_weight2, C_weight2_str, C_weight2_double, C_weight2_default, weight2Id,
							" ARM_ERR: weight 2 : numeric or ReferenceValue expected",C_result);
	XL_readStrOrNumVectorWD(XL_fixing1,C_fixing1_str,C_fixing1_double_vect,C_fixing1_default,fixing1_type," ARM_ERR: fixing taux : array of numeric expected",C_result);
	XL_readStrOrNumVectorWD(XL_fixing2,C_fixing2_str,C_fixing2_double_vect,C_fixing2_default,fixing2_type," ARM_ERR: fixing taux : array of numeric expected",C_result);

	XL_readNumVectorWD(XL_cptStrikeMethod,C_VCtcptStrikeMethod,C_VCtcptStrikeMethod_default," ARM_ERR: cptStrikeMethod : vector of integer(0/1) expected",C_result);

	if(C_VCtcptStrikeMethod.size() == 0)
	{
		C_cptStrikeMethod = C_cptStrikeMethod_default;
		C_computedFormula = C_computedFormula_default;
	}
	else if(C_VCtcptStrikeMethod.size() == 1)
	{
		C_cptStrikeMethod = C_VCtcptStrikeMethod[0];
		C_computedFormula = C_computedFormula_default;
	}
	else
	{
		C_cptStrikeMethod = C_VCtcptStrikeMethod[0];
		C_computedFormula = C_VCtcptStrikeMethod[1];
	}

	firstLegId = LocalGetNumObjectId (C_firstLeg);
	secondLegId = LocalGetNumObjectId (C_secondLeg);

	if ((XL_fixing1->xltype == xltypeMissing) || (XL_fixing1->xltype == xltypeNil))
	{
		fixing1_type = 1;
		C_fixing1_double = ARM_NULL_OBJECT;
    }
	else if ( fixing1_type == XL_TYPE_STRING )
	{
		C_fixing1_double = (double) LocalGetNumObjectId(C_fixing1_str);
		fixing1_type = 1;
	}
	else
	{
	    fixing1_type = 0;
	}

	if ((XL_fixing2->xltype == xltypeMissing) || (XL_fixing2->xltype == xltypeNil))
	{
		fixing2_type = 1;
		C_fixing2_double = ARM_NULL_OBJECT;
    }
	else if ( fixing2_type == XL_TYPE_STRING )
	{
		C_fixing2_double = (double) LocalGetNumObjectId(C_fixing2_str);
		fixing2_type = 1;
	}
	else
	{
	    fixing2_type = 0;
	}

	if((capOrFloorId = ARM_ConvCapOrFloor (C_capOrFloor, C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	if ( strike_type == XL_TYPE_STRING )
	{
	   C_strike_double = (double) LocalGetNumObjectId(C_strike_str);

	   strike_type = 1;
	}
	else
	{
	   strike_type = 0;
	}

	if(weight1Id == XL_TYPE_STRING)
	{
		weight1IsReferenceValue = true;
		weight1Id = (long) LocalGetNumObjectId(C_weight1_str);
	}

	if(weight2Id == XL_TYPE_STRING)
	{
		weight2IsReferenceValue = true;
		weight2Id = (long) LocalGetNumObjectId(C_weight2_str);
	}

	long retCode;
	long objId;
	CCString prevClass;
	
	CCString curClass = LOCAL_SPREAD_OPTION_CLASS;
	CCString stringId;

	retCode = ARMLOCAL_SPREADOPTIONWithLegs(firstLegId,
											secondLegId,
											capOrFloorId,
											(long) strike_type,
											C_strike_double,
											C_weight1_double,
											weight1IsReferenceValue,
											weight1Id,
											C_weight2_double,
											weight1IsReferenceValue,
											weight2Id,
											(long) fixing1_type,
											C_fixing1_double,
											C_fixing1_double_vect,
											(long) fixing1_type,
											C_fixing2_double,
											C_fixing2_double_vect,
											(int) C_cptStrikeMethod,
											(int) C_computedFormula,
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_PXL_ARM_SPREADOPTIONWithLegs" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


__declspec(dllexport) LPXLOPER WINAPI Local_ARM_SPREADDIGITAL(LPXLOPER XL_startDate,
															  LPXLOPER XL_endDate,
															  LPXLOPER XL_capOrFloor,
															  LPXLOPER XL_strike,
															  LPXLOPER XL_spread,//vector
															  LPXLOPER XL_payoff,
															  LPXLOPER XL_liborType1,
															  LPXLOPER XL_liborType2,
															  LPXLOPER XL_weight,//vector
															  LPXLOPER XL_daycount,
															  LPXLOPER XL_resetFreq,
															  LPXLOPER XL_payFreq,
															  LPXLOPER XL_resetTiming,
															  LPXLOPER XL_payTiming,
															  LPXLOPER XL_currency,
                                                              LPXLOPER XL_resetGap,
															  LPXLOPER XL_intRule,
															  LPXLOPER XL_stubRule,
															  LPXLOPER XL_fixing1,
															  LPXLOPER XL_fixing2)

{
	ADD_LOG("Local_ARM_SPREADDIGITAL");
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

	CCString C_capOrFloor;
	long capOrFloorId;

	double C_strike_double;
	CCString C_strike_str;
	long strike_type;
	
	double C_payoff_double;
	CCString C_payoff_str;
	long payoff_type;
	double C_payoff_default = 0.01;

	CCString C_liborType1;
	long liborType1Id;
	CCString C_liborType2;
	long liborType2Id;

	double C_weight1,C_weight2;
	VECTOR<double> C_weight;
	VECTOR<double> C_weight_default(2,1.0);

	long C_slopeFlag;
	long C_slopeFlagDefault = 1;

	double C_cptStrikeMethod;
	double C_cptStrikeMethod_default = 1.0;

	double C_computedFormula;
	double C_computedFormula_default = 1.0;

	VECTOR<CCString> C_VCt_weight;
	VECTOR<CCString> C_VCt_weight_default;

	VECTOR<double> C_VCtcalibInfo;
	VECTOR<double> C_VCtcalibInfo_default;

	CCString C_daycount;
	long dayCountId;

	CCString C_resetFreq;
	long resetFreqId;
	CCString C_payFreq;
	long payFreqId;

	CCString C_resetTiming;
	long resetTimingId;
	CCString C_payTiming;
	long payTimingId;

	CCString C_ccy;
    long ccyId;

    double C_resetGap;
    double C_resetGap_default = 10000.0;

	double C_spread1;
	double C_spread2;
	VECTOR<double> C_spread;

	double C_fixing1_double;
   	VECTOR<double> C_fixing1_double_vect;
	CCString C_fixing1_str;
	long fixing1_type;
	VECTOR<double> C_fixing1_default;

   	double C_fixing2_double;
   	VECTOR<double> C_fixing2_double_vect;
	CCString C_fixing2_str;
	long fixing2_type;
	VECTOR<double> C_fixing2_default;

	CCString C_intRule;
    long intRuleId;

	CCString C_stubRule;
    long stubRuleId;

	// error
	static int error;
	static char* reason = "";

	XL_readNumCell(XL_startDate,C_startDate," ARM_ERR: start date: date expected",C_result);
	XL_readNumCell(XL_endDate,C_endDate," ARM_ERR: end date: date expected",C_result);
	XL_readStrCell(XL_capOrFloor,C_capOrFloor," ARM_ERR: cap or floor: string expected",C_result);
	XL_readStrOrNumCell(XL_strike,C_strike_str,C_strike_double,strike_type," ARM_ERR: strike: string or numeric expected",C_result);
	
	XL_readNumVector(XL_spread,C_spread," ARM_ERR: spread vector: array of numeric expected", C_result);

	XL_readStrOrNumCell(XL_payoff,C_payoff_str,C_payoff_double,payoff_type," ARM_ERR: payoff: string or numeric expected",C_result);
	
	XL_readStrCellWD(XL_liborType1,C_liborType1,"LIBOR3M"," ARM_ERR: libor type 1 : string expected",C_result);
	XL_readStrCellWD(XL_liborType2,C_liborType2,"LIBOR3M"," ARM_ERR: libor type 2 : string expected",C_result);
	// XL_readNumVectorWD(XL_weight,C_weight,C_weight_default," ARM_ERR: weight_slopeflag_cptStrikeMethod: array of numeric expected", C_result);
	XL_readStrVectorWD(XL_weight,C_VCt_weight,C_VCt_weight_default," ARM_ERR: weight_slopeflag_cptStrikeMethod",XL_TYPE_STRING, C_result);
	XL_readStrCellWD(XL_daycount,C_daycount,"30/360"," ARM_ERR: daycount : string expected",C_result);
	XL_readStrCellWD(XL_resetFreq,C_resetFreq,"A"," ARM_ERR: reset frequency : string expected",C_result);
	XL_readStrCellWD(XL_payFreq,C_payFreq,"A"," ARM_ERR: payment frequency : string expected",C_result);
	XL_readStrCellWD(XL_resetTiming,C_resetTiming,"ADV"," ARM_ERR: reset timing : string expected",C_result);
	XL_readStrCellWD(XL_payTiming,C_payTiming,"ARR"," ARM_ERR: payment timing : string expected",C_result);
	XL_readStrCellWD(XL_currency,C_ccy,"DEFAULT"," ARM_ERR: currency : object expected",C_result);
	XL_readNumCellWD(XL_resetGap,C_resetGap,C_resetGap_default," ARM_ERR: resetGap : numeric expected",C_result);
	XL_readStrCellWD(XL_intRule,C_intRule,"DEFAULT"," ARM_ERR: intRule : string expected",C_result);
	XL_readStrCellWD(XL_stubRule,C_stubRule,"DEFAULT"," ARM_ERR: stubRule : string expected",C_result);

	XL_readStrOrNumVectorWD(XL_fixing1,C_fixing1_str,C_fixing1_double_vect,C_fixing1_default,fixing1_type," ARM_ERR: fixing taux : array of numeric expected",C_result);
	XL_readStrOrNumVectorWD(XL_fixing2,C_fixing2_str,C_fixing2_double_vect,C_fixing2_default,fixing2_type," ARM_ERR: fixing taux : array of numeric expected",C_result);

	/*
	if(C_weight.size()<2)
	{
		C_result.setMsg ("ARM_ERR: Verify the Vector of Weight");
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}
	
	C_weight1 = C_weight[0];
	C_weight2 = C_weight[1];

	if(C_weight.size() == 2)
	{
		C_slopeFlag = C_slopeFlagDefault;
		C_cptStrikeMethod = C_cptStrikeMethod_default;
		C_computedFormula = C_computedFormula_default;
	}
	else if(C_weight.size() == 3)
	{
		C_slopeFlag = C_weight[2];
		C_cptStrikeMethod = C_cptStrikeMethod_default;
		C_computedFormula = C_computedFormula_default;
	}
	else if(C_weight.size() == 4)
	{
		C_slopeFlag = C_weight[2];
		C_cptStrikeMethod = C_weight[3];
		C_computedFormula = C_computedFormula_default;
	}
	else
	{
		C_slopeFlag = C_weight[2];
		C_cptStrikeMethod = C_weight[3];
		C_computedFormula = C_weight[4];
	}
	*/

	if(C_VCt_weight.size()<2)
	{
		C_result.setMsg ("ARM_ERR: Verify the Vector of Weight");
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	sscanf ((const char*)C_VCt_weight[0], "%lf", &C_weight1);
	sscanf ((const char*)C_VCt_weight[1], "%lf", &C_weight2);

	if(C_VCt_weight.size() == 2)
	{
		C_slopeFlag = C_slopeFlagDefault;
		C_cptStrikeMethod = C_cptStrikeMethod_default;
		C_computedFormula = C_computedFormula_default;
	}
	else if(C_VCt_weight.size() == 3)
	{
		sscanf ((const char*)C_VCt_weight[2], "%ld", &C_slopeFlag);
		C_cptStrikeMethod = C_cptStrikeMethod_default;
		C_computedFormula = C_computedFormula_default;
	}
	else if(C_VCt_weight.size() == 4)
	{
		sscanf ((const char*)C_VCt_weight[2], "%ld", &C_slopeFlag);
		sscanf ((const char*)C_VCt_weight[3], "%lf", &C_cptStrikeMethod);
		C_computedFormula = C_computedFormula_default;
	}
	else
	{
		sscanf ((const char*)C_VCt_weight[2], "%ld", &C_slopeFlag);
		sscanf ((const char*)C_VCt_weight[3], "%lf", &C_cptStrikeMethod);
		sscanf ((const char*)C_VCt_weight[4], "%lf", &C_computedFormula);

		// get Sabr calib infos for each index

		double tmp;
		for (int j=0;j<C_VCt_weight.size()-5;j++)
		{
			sscanf ((const char*)C_VCt_weight[j+5], "%lf", &tmp);
			C_VCtcalibInfo.push_back (tmp);
		}
	}
		
	if(C_spread.size()<2)
	{
		C_result.setMsg ("ARM_ERR: Verify the Vector of Spread");
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}
	C_spread1 = C_spread[0];
	C_spread2 = C_spread[1];

	if ((XL_fixing1->xltype == xltypeMissing) || (XL_fixing1->xltype == xltypeNil))
	{
		fixing1_type = 1;
		C_fixing1_double = ARM_NULL_OBJECT;
    }
	else if ( fixing1_type == XL_TYPE_STRING )
	{
		C_fixing1_double = (double) LocalGetNumObjectId(C_fixing1_str);
		fixing1_type = 1;
	}
	else
	{
	    fixing1_type = 0;
	}

	if ((XL_fixing2->xltype == xltypeMissing) || (XL_fixing2->xltype == xltypeNil))
	{
		fixing2_type = 1;
		C_fixing2_double = ARM_NULL_OBJECT;
    }
	else if ( fixing2_type == XL_TYPE_STRING )
	{
		C_fixing2_double = (double) LocalGetNumObjectId(C_fixing2_str);
		fixing2_type = 1;
	}
	else
	{
	    fixing2_type = 0;
	}
	
	if((capOrFloorId = ARM_ConvCapOrFloor (C_capOrFloor, C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	liborType1Id = ARM_ConvIrType (C_liborType1);

	liborType2Id = ARM_ConvIrType (C_liborType2);

	dayCountId = ARM_ConvDayCount (C_daycount);
	
	if((resetFreqId = ARM_ConvFrequency (C_resetFreq, C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	if((payFreqId = ARM_ConvFrequency (C_payFreq, C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	resetTimingId = ARM_ConvPayResetRule (C_resetTiming);

	payTimingId = ARM_ConvPayResetRule (C_payTiming);

	if(C_ccy == "DEFAULT")
	{
		ccyId = ARM_NULL_OBJECT;
	}
	else
	{
		ccyId = LocalGetNumObjectId (C_ccy);
	}

	if ( strike_type == XL_TYPE_STRING )
	{
	   C_strike_double = (double) LocalGetNumObjectId(C_strike_str);

	   strike_type = 1;
	}
	else
	{
	   strike_type = 0;
	}

	if ( payoff_type == XL_TYPE_STRING )
	{
	   C_payoff_double = (double) LocalGetNumObjectId(C_payoff_str);

	   payoff_type = 1;
	}
	else
	{
	   payoff_type = 0;
	}

	if(C_intRule == "DEFAULT")
	{
		C_intRule = "ADJ";
	}
	intRuleId = ARM_ConvIntRule(C_intRule);

	if(C_stubRule == "DEFAULT")
	{
		C_stubRule = "SS";
	}
	stubRuleId = ARM_ConvStubRule(C_stubRule);

	long retCode;
	long objId;
	CCString prevClass;
	
	CCString curClass = LOCAL_SPREAD_OPTION_CLASS;
	CCString stringId = GetLastCurCellEnvValue ();

	if(!stringId)
	{
		retCode = ARMLOCAL_SPREADDIGITAL(C_startDate,
										 C_endDate,
										 capOrFloorId,
										 (long) strike_type,
										 C_strike_double,
										 (long) payoff_type,
										 C_payoff_double,
										 (long)liborType1Id,
										 (long)liborType2Id,
										 C_weight1,
										 C_weight2,
										 dayCountId,
										 resetFreqId,
										 payFreqId,
										 resetTimingId,
										 payTimingId,
										 ccyId,
                                         long(C_resetGap),
										 C_spread1,
										 C_spread2,
										 (long) fixing1_type,
										 C_fixing1_double,
										 C_fixing1_double_vect,
										 (long) fixing1_type,
										 C_fixing2_double,
										 C_fixing2_double_vect,
										 intRuleId,
										 stubRuleId,
										 C_slopeFlag,
										 (int) C_cptStrikeMethod,
										 (int) C_computedFormula,
										 C_VCtcalibInfo,
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
			retCode = ARMLOCAL_SPREADDIGITAL(C_startDate,
											 C_endDate,
											 capOrFloorId,
											 (long) strike_type,
											 C_strike_double,
											 (long) payoff_type,
											 C_payoff_double,
											 (long)liborType1Id,
											 (long)liborType2Id,
											 C_weight1,
											 C_weight2,
											 dayCountId,
											 resetFreqId,
											 payFreqId,
											 resetTimingId,
											 payTimingId,
											 ccyId,
                                             long(C_resetGap),
											 C_spread1,
											 C_spread2,
											 (long) fixing1_type,
											 C_fixing1_double,
											 C_fixing1_double_vect,
											 (long) fixing1_type,
											 C_fixing2_double,
											 C_fixing2_double_vect,
											 intRuleId,
											 stubRuleId,
											 C_slopeFlag,
											 (int) C_cptStrikeMethod,
											 (int) C_computedFormula,
											 C_VCtcalibInfo,
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
			retCode = ARMLOCAL_SPREADDIGITAL(C_startDate,
											 C_endDate,
											 capOrFloorId,
											 (long) strike_type,
											 C_strike_double,
										     (long) payoff_type,
											 C_payoff_double,
											 (long)liborType1Id,
											 (long)liborType2Id,
											 C_weight1,
											 C_weight2,
											 dayCountId,
											 resetFreqId,
											 payFreqId,
											 resetTimingId,
											 payTimingId,
											 ccyId,
                                             long(C_resetGap),
											 C_spread1,
											 C_spread2,
											 (long) fixing1_type,
											 C_fixing1_double,
											 C_fixing1_double_vect,
											 (long) fixing1_type,
											 C_fixing2_double,
											 C_fixing2_double_vect,
											 intRuleId,
											 stubRuleId,
											 C_slopeFlag,
											 (int) C_cptStrikeMethod,
											 (int) C_computedFormula,
											 C_VCtcalibInfo,
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_ARM_SPREADDIGITAL" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


__declspec(dllexport) LPXLOPER WINAPI Local_PXL_ARM_SPREADDIGITAL(LPXLOPER XL_startDate,
																  LPXLOPER XL_endDate,
																  LPXLOPER XL_capOrFloor,
																  LPXLOPER XL_strike,
																  LPXLOPER XL_spread,//vector
																  LPXLOPER XL_payoff,
																  LPXLOPER XL_liborType1,
																  LPXLOPER XL_liborType2,
																  LPXLOPER XL_weight,//vector
																  LPXLOPER XL_daycount,
																  LPXLOPER XL_resetFreq,
																  LPXLOPER XL_payFreq,
																  LPXLOPER XL_resetTiming,
																  LPXLOPER XL_payTiming,
																  LPXLOPER XL_currency,
																  LPXLOPER XL_resetGap,
																  LPXLOPER XL_intRule,
																  LPXLOPER XL_stubRule,
																  LPXLOPER XL_fixing1,
																  LPXLOPER XL_fixing2)

{
	ADD_LOG("Local_PXL_ARM_SPREADDIGITAL");
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

	CCString C_capOrFloor;
	long capOrFloorId;

	double C_strike_double;
	CCString C_strike_str;
	long strike_type;
	
	double C_payoff_double;
	CCString C_payoff_str;
	long payoff_type;
	double C_payoff_default = 0.01;

	CCString C_liborType1;
	long liborType1Id;
	CCString C_liborType2;
	long liborType2Id;

	double C_weight1,C_weight2;
	VECTOR<double> C_weight;
	VECTOR<double> C_weight_default(2,1.0);

	long C_slopeFlag;
	long C_slopeFlagDefault = 1;

	CCString C_daycount;
	long dayCountId;

	CCString C_resetFreq;
	long resetFreqId;
	CCString C_payFreq;
	long payFreqId;

	CCString C_resetTiming;
	long resetTimingId;
	CCString C_payTiming;
	long payTimingId;

	CCString C_ccy;
    long ccyId;

    double C_resetGap;
    double C_resetGap_default = 10000.0;

	double C_spread1;
	double C_spread2;
	VECTOR<double> C_spread;

	double C_fixing1_double;
   	VECTOR<double> C_fixing1_double_vect;
	CCString C_fixing1_str;
	long fixing1_type;
	VECTOR<double> C_fixing1_default;

   	double C_fixing2_double;
   	VECTOR<double> C_fixing2_double_vect;
	CCString C_fixing2_str;
	long fixing2_type;
	VECTOR<double> C_fixing2_default;

	CCString C_intRule;
    long intRuleId;

	CCString C_stubRule;
    long stubRuleId;

	double C_cptStrikeMethod;
	double C_cptStrikeMethod_default = 1.0;

	double C_computedFormula;
	double C_computedFormula_default = 1.0;

	VECTOR<double> C_VCtcalibInfo;
	VECTOR<double> C_VCtcalibInfo_default;

	// error
	static int error;
	static char* reason = "";

	XL_readNumCell(XL_startDate,C_startDate," ARM_ERR: start date: date expected",C_result);
	XL_readNumCell(XL_endDate,C_endDate," ARM_ERR: end date: date expected",C_result);
	XL_readStrCell(XL_capOrFloor,C_capOrFloor," ARM_ERR: cap or floor: string expected",C_result);
	XL_readStrOrNumCell(XL_strike,C_strike_str,C_strike_double,strike_type," ARM_ERR: strike: string or numeric expected",C_result);
	
	XL_readNumVector(XL_spread,C_spread," ARM_ERR: spread vector: array of numeric expected", C_result);

	XL_readStrOrNumCell(XL_payoff,C_payoff_str,C_payoff_double,payoff_type," ARM_ERR: payoff: string or numeric expected",C_result);
	
	XL_readStrCellWD(XL_liborType1,C_liborType1,"LIBOR3M"," ARM_ERR: libor type 1 : string expected",C_result);
	XL_readStrCellWD(XL_liborType2,C_liborType2,"LIBOR3M"," ARM_ERR: libor type 2 : string expected",C_result);
	XL_readNumVectorWD(XL_weight,C_weight,C_weight_default," ARM_ERR: weight_slopeflag: array of numeric expected", C_result);
	XL_readStrCellWD(XL_daycount,C_daycount,"30/360"," ARM_ERR: daycount : string expected",C_result);
	XL_readStrCellWD(XL_resetFreq,C_resetFreq,"A"," ARM_ERR: reset frequency : string expected",C_result);
	XL_readStrCellWD(XL_payFreq,C_payFreq,"A"," ARM_ERR: payment frequency : string expected",C_result);
	XL_readStrCellWD(XL_resetTiming,C_resetTiming,"ADV"," ARM_ERR: reset timing : string expected",C_result);
	XL_readStrCellWD(XL_payTiming,C_payTiming,"ARR"," ARM_ERR: payment timing : string expected",C_result);
	XL_readStrCellWD(XL_currency,C_ccy,"DEFAULT"," ARM_ERR: currency : object expected",C_result);
	XL_readNumCellWD(XL_resetGap,C_resetGap,C_resetGap_default," ARM_ERR: resetGap : numeric expected",C_result);
	XL_readStrCellWD(XL_intRule,C_intRule,"DEFAULT"," ARM_ERR: intRule : string expected",C_result);
	XL_readStrCellWD(XL_stubRule,C_stubRule,"DEFAULT"," ARM_ERR: stubRule : string expected",C_result);

	XL_readStrOrNumVectorWD(XL_fixing1,C_fixing1_str,C_fixing1_double_vect,C_fixing1_default,fixing1_type," ARM_ERR: fixing taux : array of numeric expected",C_result);
	XL_readStrOrNumVectorWD(XL_fixing2,C_fixing2_str,C_fixing2_double_vect,C_fixing2_default,fixing2_type," ARM_ERR: fixing taux : array of numeric expected",C_result);

	if(C_weight.size()<2)
	{
		C_result.setMsg ("ARM_ERR: Verify the Vector of Weight");
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}
	
	C_weight1 = C_weight[0];
	C_weight2 = C_weight[1];

	if(C_weight.size() == 2)
	{
		C_slopeFlag = C_slopeFlagDefault;
		C_cptStrikeMethod = C_cptStrikeMethod_default;
		C_computedFormula = C_computedFormula_default;
	}
	else if(C_weight.size() == 3)
	{
		C_slopeFlag = C_weight[2];
		C_cptStrikeMethod = C_cptStrikeMethod_default;
		C_computedFormula = C_computedFormula_default;
	}
	else if(C_weight.size() == 4)
	{
		C_slopeFlag = C_weight[2];
		C_cptStrikeMethod = C_weight[3];
		C_computedFormula = C_computedFormula_default;
	}
	else
	{
		C_slopeFlag = C_weight[2];
		C_cptStrikeMethod = C_weight[3];
		C_computedFormula = C_weight[4];
	}

	if(C_spread.size()<2)
	{
		C_result.setMsg ("ARM_ERR: Verify the Vector of Spread");
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}
	C_spread1 = C_spread[0];
	C_spread2 = C_spread[1];

	if ((XL_fixing1->xltype == xltypeMissing) || (XL_fixing1->xltype == xltypeNil))
	{
		fixing1_type = 1;
		C_fixing1_double = ARM_NULL_OBJECT;
    }
	else if ( fixing1_type == XL_TYPE_STRING )
	{
		C_fixing1_double = (double) LocalGetNumObjectId(C_fixing1_str);
		fixing1_type = 1;
	}
	else
	{
	    fixing1_type = 0;
	}

	if ((XL_fixing2->xltype == xltypeMissing) || (XL_fixing2->xltype == xltypeNil))
	{
		fixing2_type = 1;
		C_fixing2_double = ARM_NULL_OBJECT;
    }
	else if ( fixing2_type == XL_TYPE_STRING )
	{
		C_fixing2_double = (double) LocalGetNumObjectId(C_fixing2_str);
		fixing2_type = 1;
	}
	else
	{
	    fixing2_type = 0;
	}
	
	if((capOrFloorId = ARM_ConvCapOrFloor (C_capOrFloor, C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	liborType1Id = ARM_ConvIrType (C_liborType1);

	liborType2Id = ARM_ConvIrType (C_liborType2);

	dayCountId = ARM_ConvDayCount (C_daycount);
	
	if((resetFreqId = ARM_ConvFrequency (C_resetFreq, C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	if((payFreqId = ARM_ConvFrequency (C_payFreq, C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	resetTimingId = ARM_ConvPayResetRule (C_resetTiming);

	payTimingId = ARM_ConvPayResetRule (C_payTiming);

	if(C_ccy == "DEFAULT")
	{
		ccyId = ARM_NULL_OBJECT;
	}
	else
	{
		ccyId = LocalGetNumObjectId (C_ccy);
	}

	if ( strike_type == XL_TYPE_STRING )
	{
	   C_strike_double = (double) LocalGetNumObjectId(C_strike_str);

	   strike_type = 1;
	}
	else
	{
	   strike_type = 0;
	}

	if ( payoff_type == XL_TYPE_STRING )
	{
	   C_payoff_double = (double) LocalGetNumObjectId(C_payoff_str);

	   payoff_type = 1;
	}
	else
	{
	   payoff_type = 0;
	}

	if(C_intRule == "DEFAULT")
	{
		C_intRule = "ADJ";
	}
	intRuleId = ARM_ConvIntRule(C_intRule);

	if(C_stubRule == "DEFAULT")
	{
		C_stubRule = "SS";
	}
	stubRuleId = ARM_ConvStubRule(C_stubRule);

	long retCode;
	long objId;
	CCString prevClass;
	
	CCString curClass = LOCAL_SPREAD_OPTION_CLASS;
	CCString stringId;

	retCode = ARMLOCAL_SPREADDIGITAL(C_startDate,
									 C_endDate,
									 capOrFloorId,
									 (long) strike_type,
									 C_strike_double,
									 (long) payoff_type,
									 C_payoff_double,
									 (long)liborType1Id,
									 (long)liborType2Id,
									 C_weight1,
									 C_weight2,
									 dayCountId,
									 resetFreqId,
									 payFreqId,
									 resetTimingId,
									 payTimingId,
									 ccyId,
									 long(C_resetGap),
									 C_spread1,
									 C_spread2,
									 (long) fixing1_type,
									 C_fixing1_double,
									 C_fixing1_double_vect,
									 (long) fixing1_type,
									 C_fixing2_double,
									 C_fixing2_double_vect,
									 intRuleId,
									 stubRuleId,
									 C_slopeFlag,
									 (int) C_cptStrikeMethod,
									 (int) C_computedFormula,
									 C_VCtcalibInfo,
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_PXL_ARM_SPREADDIGITAL" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}

__declspec(dllexport) LPXLOPER WINAPI Local_ARM_QUANTOSPREADDIGITAL(	LPXLOPER XL_startDate,
																		LPXLOPER XL_endDate,
 																		LPXLOPER XL_capOrFloor,
																		LPXLOPER XL_strikes,
																		LPXLOPER XL_payoff,
																		LPXLOPER XL_liborIdx1,
																		LPXLOPER XL_liborIdx2,
																		LPXLOPER XL_leg1Weight,
																		LPXLOPER XL_leg2Weight,
																		LPXLOPER XL_leg1Fixings,
																		LPXLOPER XL_leg2Fixings,
																		LPXLOPER XL_leg1Spread,
																		LPXLOPER XL_leg2Spread,
																		LPXLOPER XL_modelParams,
																		LPXLOPER XL_prodParams,
																		LPXLOPER XL_notional)
{
	ADD_LOG("Local_ARM_QUANTOSPREADDIGITAL");
	//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();

	double C_startDate;
	double C_endDate;

	CCString C_capOrFloor;
	long capOrFloorId;

	CCString C_strikes_str;
	double C_strike_double;
	long strike_type;

	CCString C_payoff_str;
	double C_payoff_double;
	long payoff_type;

	CCString C_liborIdx1;
	long liborIdx1Id;
	CCString C_liborIdx2;
	long liborIdx2Id;
	
	double C_Idx1weight;
	double C_Idx1weight_default = 1.0;

	double C_Idx2weight;
	double C_Idx2weight_default = 1.0;

	CCString C_Idx1fixings_str;
	CCString C_Idx1fixings_str_default("DEFAULT");
	long Idx1fixingsId = ARM_NULL_OBJECT;
	CCString C_Idx2fixings_str;
	CCString C_Idx2fixings_str_default("DEFAULT");
	long Idx2fixingsId = ARM_NULL_OBJECT;
	
	double C_Idx1spread;
	double C_Idx1spread_default = 0.0;

	double C_Idx2spread;
	double C_Idx2spread_default = 0.0;

	VECTOR<double> C_modelParams;
	VECTOR<double> C_modelParams_default(0);
	VECTOR<double> modelParams(3);
	
	modelParams[0] = 1.0;
	modelParams[1] = 1.0;
	modelParams[2] = 1.0;

	double C_slopeFlag;
	double C_cptStrikeMethod;
	double C_computedFormula;

	// Product parameters
	VECTOR<CCString> C_prodParams;
	VECTOR<CCString> C_prodParams_default(0);

	VECTOR<CCString> prodParams(13);

	prodParams[0] = "DEFAULT";
	prodParams[1] = "30/360";
	prodParams[2] = "A";
	prodParams[3] = "A";
	prodParams[4] = "ADV";
	prodParams[5] = "ARR";
	prodParams[6] = "ADJ";
	prodParams[7] = "SS";
	prodParams[8] = "2";
	prodParams[9] = "EUR";
	prodParams[10] = "EUR";
	prodParams[11] = "MF";
	prodParams[12] = "NULL";

	CCString C_ccy;
	long ccyId;
	
	CCString C_daycount;
	long dayCountId;

	CCString C_resetFreq;
	long resetFreqId;

	CCString C_payFreq;
	long payFreqId;

	CCString C_resetTiming;
	long resetTimingId;
	
	CCString C_payTiming;
	long payTimingId;

	CCString C_intRule;
    long intRuleId;
	
	CCString C_stubRule;
    long stubRuleId;

    double C_resetGap;
	
	CCString C_resetCal;
	CCString C_payCal;
	CCString C_fwdRule;
	double fwdRuleId;

	CCString C_refDate;
	
	long notionalId = ARM_NULL_OBJECT;
	CCString C_notional_str;
	CCString C_notional_default("DEFAULT");

	// error
	static int error;
	static char* reason = "";

	XL_readNumCell(XL_startDate,C_startDate," ARM_ERR: start date: date expected",C_result);
	XL_readNumCell(XL_endDate,C_endDate," ARM_ERR: end date: date expected",C_result);
	XL_readStrCell(XL_capOrFloor,C_capOrFloor," ARM_ERR: cap or floor: string expected",C_result);
	XL_readStrOrNumCell(XL_strikes,C_strikes_str,C_strike_double,strike_type," ARM_ERR: strikes vector : double or array of numeric expected",C_result);
	XL_readStrOrNumCell(XL_payoff,C_payoff_str,C_payoff_double,payoff_type," ARM_ERR: Pay Off vector : double or array of numeric expected",C_result);
	XL_readStrCellWD(XL_liborIdx1,C_liborIdx1,"DEFAULT"," ARM_ERR: leg1 libor index : string expected",C_result);
	XL_readStrCellWD(XL_liborIdx2,C_liborIdx2,"DEFAULT"," ARM_ERR: leg2 libor index : string expected",C_result);
	XL_readNumCellWD(XL_leg1Weight,C_Idx1weight,C_Idx1weight_default," ARM_ERR: leg1 weight : numeric expected",C_result);
	XL_readNumCellWD(XL_leg2Weight,C_Idx2weight,C_Idx2weight_default," ARM_ERR: leg2 weight : numeric expected",C_result);
	XL_readStrCellWD(XL_leg1Fixings,C_Idx1fixings_str,C_Idx1fixings_str_default," ARM_ERR: leg1 fixings vector : object expected",C_result);
	XL_readStrCellWD(XL_leg2Fixings,C_Idx2fixings_str,C_Idx2fixings_str_default," ARM_ERR: leg2 fixings vector : object expected",C_result);
	XL_readNumCellWD(XL_leg1Spread,C_Idx1spread, C_Idx1spread_default, " ARM_ERR: leg1 spread : numeric expected",C_result);
	XL_readNumCellWD(XL_leg2Spread,C_Idx2spread,C_Idx2spread_default, " ARM_ERR: leg2 spread : numeric expected",C_result);
	XL_readNumVectorWD(XL_modelParams,C_modelParams,C_modelParams_default," ARM_ERR: model parameters vector : array of numeric expected",C_result);
	XL_readStrVectorWD(XL_prodParams,C_prodParams,C_prodParams_default," ARM_ERR: product parameters vector : array of numeric expected",DOUBLE_TYPE,C_result);
	XL_readStrCellWD(XL_notional,C_notional_str,C_notional_default," ARM_ERR: notional vector : object expected",C_result);

	if((capOrFloorId = ARM_ConvCapOrFloor (C_capOrFloor, C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}
		
	if ( strike_type == XL_TYPE_STRING )
	{
	   C_strike_double = (double) LocalGetNumObjectId(C_strikes_str);

	   strike_type = 1;
	}
	else
	{
	   strike_type = 0;
	}

	if ( payoff_type == XL_TYPE_STRING )
	{
	   C_payoff_double = (double) LocalGetNumObjectId(C_payoff_str);

	   payoff_type = 1;
	}
	else
	{
	   payoff_type = 0;
	}

	liborIdx1Id = strcmp(C_liborIdx1,"DEFAULT") ? (double) LocalGetNumObjectId(C_liborIdx1) : -1.0;	
	liborIdx2Id = strcmp(C_liborIdx2,"DEFAULT") ? (double) LocalGetNumObjectId(C_liborIdx2) : -1.0;	
	
	Idx1fixingsId = strcmp(C_Idx1fixings_str,"DEFAULT") ? (long) LocalGetNumObjectId(C_Idx1fixings_str) : -1.0;
	Idx2fixingsId = strcmp(C_Idx2fixings_str,"DEFAULT") ? (long) LocalGetNumObjectId(C_Idx2fixings_str) : -1.0;
	
	int i = 0;
	for (i; i < C_modelParams.size(); i++)
	{
		modelParams[i]  = C_modelParams[i];
	}
	
	C_slopeFlag = modelParams[0];
	C_cptStrikeMethod = modelParams[1];
	C_computedFormula = modelParams[2];

	for (i = 0; i<C_prodParams.size(); i++)
	{
		if (strcmp(C_prodParams[i],"DEFAULT"))
			prodParams[i] = C_prodParams[i];
	}

	C_ccy = prodParams[0];
	
	if(C_ccy == "DEFAULT")
	{
		ccyId = ARM_NULL_OBJECT;
	}
	else
	{
		ccyId = LocalGetNumObjectId (C_ccy);
	}
	
	C_daycount = prodParams[1];
	dayCountId = ARM_ConvDayCount (C_daycount);
	
	C_resetFreq = prodParams[2];
	if((resetFreqId = ARM_ConvFrequency (C_resetFreq, C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}
	
	C_payFreq = prodParams[3];
	if((payFreqId = ARM_ConvFrequency (C_payFreq, C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	C_resetTiming = prodParams[4];
	resetTimingId = ARM_ConvPayResetRule (C_resetTiming);
	
	C_payTiming = prodParams[5];
	payTimingId = ARM_ConvPayResetRule (C_payTiming);

	C_intRule  = prodParams[6];
    intRuleId = ARM_ConvIntRule(C_intRule);
	
	C_stubRule = prodParams[7];
	stubRuleId = ARM_ConvStubRule(C_stubRule);

    C_resetGap = atof(prodParams[8]);
	
	C_resetCal = prodParams[9];

	C_payCal = prodParams[10];

	C_fwdRule = prodParams[11];
	fwdRuleId = ARM_ConvFwdRule(C_fwdRule);

	C_refDate = prodParams[12];

	notionalId = strcmp(C_notional_str,"DEFAULT") ? (long) LocalGetNumObjectId(C_notional_str) : -1.0;

	long retCode;
	long objId;
	
	CCString prevClass;
	CCString curClass = LOCAL_SPREAD_OPTION_CLASS;
	CCString stringId = GetLastCurCellEnvValue();
	
	if(!stringId)
	{
		retCode = ARMLOCAL_QUANTOSPREADDIGITAL (	C_startDate,
													C_endDate,
													capOrFloorId,
													(long) strike_type,
													C_strike_double,
													(long) payoff_type,
													C_payoff_double,
													(long) liborIdx1Id,
													(long) liborIdx2Id,
													C_Idx1weight,
													C_Idx2weight,
													Idx1fixingsId,
													Idx2fixingsId,
													C_Idx1spread,
													C_Idx2spread,
													C_slopeFlag,
													C_cptStrikeMethod,
													(int) C_computedFormula,
													ccyId,
													dayCountId,
													resetFreqId,
													payFreqId,
													resetTimingId,
													payTimingId,
													intRuleId,
													stubRuleId,
													C_resetGap,
													C_resetCal,
													C_payCal,
													fwdRuleId,
													C_refDate,
													notionalId,
													C_result);
		if (retCode == ARM_OK)
		{
			objId = C_result.getLong();

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
			retCode = ARMLOCAL_QUANTOSPREADDIGITAL( C_startDate,
													C_endDate,
													capOrFloorId,
													(long) strike_type,
													C_strike_double,
													(long) payoff_type,
													C_payoff_double,
													(long) liborIdx1Id,
													(long) liborIdx2Id,
													C_Idx1weight,
													C_Idx2weight,
													Idx1fixingsId,
													Idx2fixingsId,
													C_Idx1spread,
													C_Idx2spread,
													C_slopeFlag,
													C_cptStrikeMethod,
													(int) C_computedFormula,
													ccyId,
													dayCountId,
													resetFreqId,
													payFreqId,
													resetTimingId,
													payTimingId,
													intRuleId,
													stubRuleId,
													C_resetGap,
													C_resetCal,
													C_payCal,
													fwdRuleId,
													C_refDate,
													notionalId,
													C_result,
													objId);
			if (retCode == ARM_OK)
			{
				objId = C_result.getLong();

				stringId = LocalMakeObjectId (objId, curClass);
			}
		}
		else
		{
			FreeCurCellContent ();
			retCode = ARMLOCAL_QUANTOSPREADDIGITAL (C_startDate,
													C_endDate,
													capOrFloorId,
													(long) strike_type,
													C_strike_double,
													(long) payoff_type,
													C_payoff_double,
													(long) liborIdx1Id,
													(long) liborIdx2Id,
													C_Idx1weight,
													C_Idx2weight,
													Idx1fixingsId,
													Idx2fixingsId,
													C_Idx1spread,
													C_Idx2spread,
													C_slopeFlag,
													C_cptStrikeMethod,
													(int) C_computedFormula,
													ccyId,
													dayCountId,
													resetFreqId,
													payFreqId,
													resetTimingId,
													payTimingId,
													intRuleId,
													stubRuleId,
													C_resetGap,
													C_resetCal,
													C_payCal,
													fwdRuleId,
													C_refDate,
													notionalId,
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_ARM_QUANTOSPREADDIGITAL  "  )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}

__declspec(dllexport) LPXLOPER WINAPI Local_ARM_SPREADDIGITALWithLegs(LPXLOPER XL_firstLeg,
																	  LPXLOPER XL_secondLeg,
																	  LPXLOPER XL_capOrFloor,
																	  LPXLOPER XL_strike,
																	  LPXLOPER XL_spread1,
																	  LPXLOPER XL_spread2,
																	  LPXLOPER XL_payoff,
																	  LPXLOPER XL_weight1,
																	  LPXLOPER XL_weight2,
																	  LPXLOPER XL_fixing1,
																	  LPXLOPER XL_fixing2,
																	  LPXLOPER XL_slopeFlag,
																	  LPXLOPER XL_cptStrikeMethod)

{
	ADD_LOG("Local_ARM_SPREADDIGITALWithLegs");
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();

	// C variable
	CCString C_firstLeg;
	long firstLegId;

	CCString C_secondLeg;
	long secondLegId;

	CCString C_capOrFloor;
	long capOrFloorId;

	double C_strike_double;
	CCString C_strike_str;
	long strike_type;
	
	double C_payoff_double;
	CCString C_payoff_str;
	long payoff_type;
	double C_payoff_default = 0.01;

	double C_weight1;
	double C_weight1_default = 1.0;
	double C_weight2;
	double C_weight2_default = 1.0;

	double C_spread1;
	double C_spread2;

	double C_fixing1_double;
   	VECTOR<double> C_fixing1_double_vect;
	CCString C_fixing1_str;
	long fixing1_type;
	VECTOR<double> C_fixing1_default;

   	double C_fixing2_double;
   	VECTOR<double> C_fixing2_double_vect;
	CCString C_fixing2_str;
	long fixing2_type;
	VECTOR<double> C_fixing2_default;

	double C_slopeFlag;
	double C_slopeFlagDefault = 1.0; 

	VECTOR<double> C_VCtcptStrikeMethod;
	VECTOR<double> C_VCtcptStrikeMethod_default;

	double C_cptStrikeMethod;
	double C_cptStrikeMethod_default = 1.0;

	double C_computedFormula;
	double C_computedFormula_default = 1.0;

	// error
	static int error;
	static char* reason = "";

	XL_readStrCell(XL_firstLeg,C_firstLeg," ARM_ERR: first leg id: object expected",C_result);
	XL_readStrCell(XL_secondLeg,C_secondLeg," ARM_ERR: secondleg id: object expected",C_result);

	XL_readStrCell(XL_capOrFloor,C_capOrFloor," ARM_ERR: cap or floor: string expected",C_result);
	XL_readStrOrNumCell(XL_strike,C_strike_str,C_strike_double,strike_type," ARM_ERR: strike: string or numeric expected",C_result);
	
	XL_readNumCell(XL_spread1,C_spread1," ARM_ERR: spread vector: numeric expected", C_result);
	XL_readNumCell(XL_spread2,C_spread2," ARM_ERR: spread vector: numeric expected", C_result);

	XL_readStrOrNumCell(XL_payoff,C_payoff_str,C_payoff_double,payoff_type," ARM_ERR: payoff: string or numeric expected",C_result);
	
	XL_readNumCellWD(XL_weight1,C_weight1,C_weight1_default," ARM_ERR: weight 1 : numeric expected",C_result);
	XL_readNumCellWD(XL_weight2,C_weight2,C_weight2_default," ARM_ERR: weight 2 : numeric expected",C_result);

	XL_readStrOrNumVectorWD(XL_fixing1,C_fixing1_str,C_fixing1_double_vect,C_fixing1_default,fixing1_type," ARM_ERR: fixing taux : array of numeric expected",C_result);
	XL_readStrOrNumVectorWD(XL_fixing2,C_fixing2_str,C_fixing2_double_vect,C_fixing2_default,fixing2_type," ARM_ERR: fixing taux : array of numeric expected",C_result);

	XL_readNumCellWD(XL_slopeFlag,C_slopeFlag,C_slopeFlagDefault," ARM_ERR: slope flag : numeric expected",C_result);

	XL_readNumVectorWD(XL_cptStrikeMethod,C_VCtcptStrikeMethod,C_VCtcptStrikeMethod_default," ARM_ERR: cptStrikeMethod : vector of integer(0/1) expected",C_result);

	if(C_VCtcptStrikeMethod.size() == 0)
	{
		C_cptStrikeMethod = C_cptStrikeMethod_default;
		C_computedFormula = C_computedFormula_default;
	}
	else if(C_VCtcptStrikeMethod.size() == 1)
	{
		C_cptStrikeMethod = C_VCtcptStrikeMethod[0];
		C_computedFormula = C_computedFormula_default;
	}
	else
	{
		C_cptStrikeMethod = C_VCtcptStrikeMethod[0];
		C_computedFormula = C_VCtcptStrikeMethod[1];
	}

	firstLegId = LocalGetNumObjectId (C_firstLeg);
	secondLegId = LocalGetNumObjectId (C_secondLeg);

	if ((XL_fixing1->xltype == xltypeMissing) || (XL_fixing1->xltype == xltypeNil))
	{
		fixing1_type = 1;
		C_fixing1_double = ARM_NULL_OBJECT;
    }
	else if ( fixing1_type == XL_TYPE_STRING )
	{
		C_fixing1_double = (double) LocalGetNumObjectId(C_fixing1_str);
		fixing1_type = 1;
	}
	else
	{
	    fixing1_type = 0;
	}

	if ((XL_fixing2->xltype == xltypeMissing) || (XL_fixing2->xltype == xltypeNil))
	{
		fixing2_type = 1;
		C_fixing2_double = ARM_NULL_OBJECT;
    }
	else if ( fixing2_type == XL_TYPE_STRING )
	{
		C_fixing2_double = (double) LocalGetNumObjectId(C_fixing2_str);
		fixing2_type = 1;
	}
	else
	{
	    fixing2_type = 0;
	}
	
	if((capOrFloorId = ARM_ConvCapOrFloor (C_capOrFloor, C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	if ( strike_type == XL_TYPE_STRING )
	{
	   C_strike_double = (double) LocalGetNumObjectId(C_strike_str);

	   strike_type = 1;
	}
	else
	{
	   strike_type = 0;
	}

	if ( payoff_type == XL_TYPE_STRING )
	{
	   C_payoff_double = (double) LocalGetNumObjectId(C_payoff_str);

	   payoff_type = 1;
	}
	else
	{
	   payoff_type = 0;
	}

	long retCode;
	long objId;
	CCString prevClass;
	
	CCString curClass = LOCAL_SPREAD_OPTION_CLASS;
	CCString stringId = GetLastCurCellEnvValue ();

	if(!stringId)
	{
		retCode = ARMLOCAL_SPREADDIGITALWithLegs(firstLegId,
												 secondLegId,
												 capOrFloorId,
												 (long) strike_type,
												 C_strike_double,
												 (long) payoff_type,
												 C_payoff_double,
												 C_weight1,
												 C_weight2,
												 C_spread1,
												 C_spread2,
												 (long) fixing1_type,
												 C_fixing1_double,
												 C_fixing1_double_vect,
												 (long) fixing1_type,
												 C_fixing2_double,
												 C_fixing2_double_vect,
												 (long) C_slopeFlag,
												 (int) C_cptStrikeMethod,
												 (int) C_computedFormula,
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
			retCode = ARMLOCAL_SPREADDIGITALWithLegs(firstLegId,
													 secondLegId,
													 capOrFloorId,
													 (long) strike_type,
													 C_strike_double,
													 (long) payoff_type,
													 C_payoff_double,
													 C_weight1,
													 C_weight2,
													 C_spread1,
													 C_spread2,
													 (long) fixing1_type,
													 C_fixing1_double,
													 C_fixing1_double_vect,
													 (long) fixing1_type,
													 C_fixing2_double,
													 C_fixing2_double_vect,
													 (long) C_slopeFlag,
													 (int) C_cptStrikeMethod,
													 (int) C_computedFormula,
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
			retCode = ARMLOCAL_SPREADDIGITALWithLegs(firstLegId,
													 secondLegId,
													 capOrFloorId,
													 (long) strike_type,
													 C_strike_double,
													 (long) payoff_type,
													 C_payoff_double,
													 C_weight1,
													 C_weight2,
													 C_spread1,
													 C_spread2,
													 (long) fixing1_type,
													 C_fixing1_double,
													 C_fixing1_double_vect,
													 (long) fixing1_type,
													 C_fixing2_double,
													 C_fixing2_double_vect,
													 (long) C_slopeFlag,
													 (int) C_cptStrikeMethod,
													 (int) C_computedFormula,
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_ARM_SPREADDIGITALWithLegs" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


__declspec(dllexport) LPXLOPER WINAPI Local_PXL_ARM_SPREADDIGITALWithLegs(LPXLOPER XL_firstLeg,
																		  LPXLOPER XL_secondLeg,
																		  LPXLOPER XL_capOrFloor,
																		  LPXLOPER XL_strike,
																		  LPXLOPER XL_spread1,
																		  LPXLOPER XL_spread2,
																		  LPXLOPER XL_payoff,
																		  LPXLOPER XL_weight1,
																		  LPXLOPER XL_weight2,
																		  LPXLOPER XL_fixing1,
																		  LPXLOPER XL_fixing2,
																		  LPXLOPER XL_slopeFlag,
																		  LPXLOPER XL_cptStrikeMethod)
{
	ADD_LOG("Local_PXL_ARM_SPREADDIGITALWithLegs");
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();

	// C variable
	CCString C_firstLeg;
	long firstLegId;

	CCString C_secondLeg;
	long secondLegId;

	CCString C_capOrFloor;
	long capOrFloorId;

	double C_strike_double;
	CCString C_strike_str;
	long strike_type;
	
	double C_payoff_double;
	CCString C_payoff_str;
	long payoff_type;
	double C_payoff_default = 0.01;

	double C_weight1;
	double C_weight1_default = 1.0;
	double C_weight2;
	double C_weight2_default = 1.0;

	double C_spread1;
	double C_spread2;

	double C_fixing1_double;
   	VECTOR<double> C_fixing1_double_vect;
	CCString C_fixing1_str;
	long fixing1_type;
	VECTOR<double> C_fixing1_default;

   	double C_fixing2_double;
   	VECTOR<double> C_fixing2_double_vect;
	CCString C_fixing2_str;
	long fixing2_type;
	VECTOR<double> C_fixing2_default;

	double C_slopeFlag;
	double C_slopeFlagDefault = 1.0; 

	VECTOR<double> C_VCtcptStrikeMethod;
	VECTOR<double> C_VCtcptStrikeMethod_default;

	double C_cptStrikeMethod;
	double C_cptStrikeMethod_default = 1.0;

	double C_computedFormula;
	double C_computedFormula_default = 1.0;

	// error
	static int error;
	static char* reason = "";

	XL_readStrCell(XL_firstLeg,C_firstLeg," ARM_ERR: first leg id: object expected",C_result);
	XL_readStrCell(XL_secondLeg,C_secondLeg," ARM_ERR: secondleg id: object expected",C_result);

	XL_readStrCell(XL_capOrFloor,C_capOrFloor," ARM_ERR: cap or floor: string expected",C_result);
	XL_readStrOrNumCell(XL_strike,C_strike_str,C_strike_double,strike_type," ARM_ERR: strike: string or numeric expected",C_result);
	
	XL_readNumCell(XL_spread1,C_spread1," ARM_ERR: spread vector: numeric expected", C_result);
	XL_readNumCell(XL_spread2,C_spread2," ARM_ERR: spread vector: numeric expected", C_result);

	XL_readStrOrNumCell(XL_payoff,C_payoff_str,C_payoff_double,payoff_type," ARM_ERR: payoff: string or numeric expected",C_result);
	
	XL_readNumCellWD(XL_weight1,C_weight1,C_weight1_default," ARM_ERR: weight 1 : numeric expected",C_result);
	XL_readNumCellWD(XL_weight2,C_weight2,C_weight2_default," ARM_ERR: weight 2 : numeric expected",C_result);

	XL_readStrOrNumVectorWD(XL_fixing1,C_fixing1_str,C_fixing1_double_vect,C_fixing1_default,fixing1_type," ARM_ERR: fixing taux : array of numeric expected",C_result);
	XL_readStrOrNumVectorWD(XL_fixing2,C_fixing2_str,C_fixing2_double_vect,C_fixing2_default,fixing2_type," ARM_ERR: fixing taux : array of numeric expected",C_result);

	XL_readNumCellWD(XL_slopeFlag,C_slopeFlag,C_slopeFlagDefault," ARM_ERR: slope flag : numeric expected",C_result);

	XL_readNumVectorWD(XL_cptStrikeMethod,C_VCtcptStrikeMethod,C_VCtcptStrikeMethod_default," ARM_ERR: cptStrikeMethod : vector of integer(0/1) expected",C_result);

	if(C_VCtcptStrikeMethod.size() == 0)
	{
		C_cptStrikeMethod = C_cptStrikeMethod_default;
		C_computedFormula = C_computedFormula_default;
	}
	else if(C_VCtcptStrikeMethod.size() == 1)
	{
		C_cptStrikeMethod = C_VCtcptStrikeMethod[0];
		C_computedFormula = C_computedFormula_default;
	}
	else
	{
		C_cptStrikeMethod = C_VCtcptStrikeMethod[0];
		C_computedFormula = C_VCtcptStrikeMethod[1];
	}

	firstLegId = LocalGetNumObjectId (C_firstLeg);
	secondLegId = LocalGetNumObjectId (C_secondLeg);

	if ((XL_fixing1->xltype == xltypeMissing) || (XL_fixing1->xltype == xltypeNil))
	{
		fixing1_type = 1;
		C_fixing1_double = ARM_NULL_OBJECT;
    }
	else if ( fixing1_type == XL_TYPE_STRING )
	{
		C_fixing1_double = (double) LocalGetNumObjectId(C_fixing1_str);
		fixing1_type = 1;
	}
	else
	{
	    fixing1_type = 0;
	}

	if ((XL_fixing2->xltype == xltypeMissing) || (XL_fixing2->xltype == xltypeNil))
	{
		fixing2_type = 1;
		C_fixing2_double = ARM_NULL_OBJECT;
    }
	else if ( fixing2_type == XL_TYPE_STRING )
	{
		C_fixing2_double = (double) LocalGetNumObjectId(C_fixing2_str);
		fixing2_type = 1;
	}
	else
	{
	    fixing2_type = 0;
	}
	
	if((capOrFloorId = ARM_ConvCapOrFloor (C_capOrFloor, C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	if ( strike_type == XL_TYPE_STRING )
	{
	   C_strike_double = (double) LocalGetNumObjectId(C_strike_str);

	   strike_type = 1;
	}
	else
	{
	   strike_type = 0;
	}

	if ( payoff_type == XL_TYPE_STRING )
	{
	   C_payoff_double = (double) LocalGetNumObjectId(C_payoff_str);

	   payoff_type = 1;
	}
	else
	{
	   payoff_type = 0;
	}

	long retCode;
	long objId;
	CCString prevClass;
	
	CCString curClass = LOCAL_SPREAD_OPTION_CLASS;
	CCString stringId;

	retCode = ARMLOCAL_SPREADDIGITALWithLegs(firstLegId,
											 secondLegId,
											 capOrFloorId,
											 (long) strike_type,
											 C_strike_double,
											 (long) payoff_type,
											 C_payoff_double,
											 C_weight1,
											 C_weight2,
											 C_spread1,
											 C_spread2,
											 (long) fixing1_type,
											 C_fixing1_double,
											 C_fixing1_double_vect,
											 (long) fixing1_type,
											 C_fixing2_double,
											 C_fixing2_double_vect,
											 (long) C_slopeFlag,
											 (int) C_cptStrikeMethod,
											 (int) C_computedFormula,
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_PXL_ARM_SPREADDIGITALWithLegs" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}

__declspec(dllexport) LPXLOPER WINAPI Local_ARM_SPREADCORRIDOR(LPXLOPER XL_startDate,
															   LPXLOPER XL_endDate,
 															   LPXLOPER XL_capOrFloor,
															   LPXLOPER XL_strike,
															   LPXLOPER XL_spreads,//Vect: spread1, spread2
															   LPXLOPER XL_payIndexParams,//Vect:payIndexObj, itsMargin,itsWeight
															   LPXLOPER XL_spreadIdxTypes,//Vect: IndexType1, IndexType2
															   LPXLOPER XL_sprdWeights_slope,//Vect: Idx Weight1,Weight2, slopeflag and cptStrikeMethod
															   LPXLOPER XL_prodParamDatas,//Vect: dayCount, resetFreq, payFreq, resetTiming, payTiming, intRule, stubRule,resetGap 
															   LPXLOPER XL_currency,
															   LPXLOPER XL_fixing1,
															   LPXLOPER XL_fixing2,
															   LPXLOPER XL_fixing3,
															   LPXLOPER XL_freezeFixing)
{
	ADD_LOG("Local_ARM_SPREADCORRIDOR");
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

	CCString C_capOrFloor;
	long capOrFloorId;

	double C_strike_double;
	CCString C_strike_str;
	long strike_type;
	
	double C_spread1;
	double C_spread2;
	VECTOR<double> C_spreads;

	CCString C_payIndex;
	long payIndexId;

	CCString C_payMargin_str;
	long payMargin_type;
	double payMarginId;

	double C_payIdxWeight;

	CCString C_payIdxFixedRate_str;
		
	VECTOR<CCString> payIndexParams(3);
   	VECTOR<CCString> C_payIndexParams;
	VECTOR<CCString> C_payIndexParamsDF(3);
	C_payIndexParamsDF[0] = payIndexParams[0] = "DEFAULT";// pay index objet
	C_payIndexParamsDF[1] = payIndexParams[1] = "DEFAULT";      // Margin
	C_payIndexParamsDF[2] = payIndexParams[2] = "1";      // weight
//	C_payIndexParamsDF[3] = payIndexParams[3] = "DEFAULT";      // fixed rate
	
	CCString C_liborType1;
	long liborType1Id;
	CCString C_liborType2;
	long liborType2Id;
	VECTOR<CCString> liborTypes(2);
   	VECTOR<CCString> C_liborTypes;
	VECTOR<CCString> C_liborTypesDF(2);
	C_liborTypesDF[0] = liborTypes[0] = "LIBOR3M";
	C_liborTypesDF[1] = liborTypes[1] = "LIBOR3M";

	double C_weight1, C_weight2, C_cptStrikeMethod;
	VECTOR<double> weights_slope(2,1.0);
	VECTOR<double> C_weights_slope;
	VECTOR<double> C_weights_slope_df(2, 1.0);

	long C_slopeFlag;
	long C_slopeFlagDefault = 1;

	double C_cptStrikeMethod_default = 1.0;

	double C_computedFormula;
	double C_computedFormula_default = 1.0;

	// prodParamDatas
	VECTOR<CCString> prodParamDatas(8);
	VECTOR<CCString> C_prodParamDatas;
	VECTOR<CCString> C_prodParamDatasDF(8);
	C_prodParamDatasDF[0] = prodParamDatas[0] = "30/360"; // dayCount
	C_prodParamDatasDF[1] = prodParamDatas[1] = "A"; // resetFreq
	C_prodParamDatasDF[2] = prodParamDatas[2] = "A"; // payFreq
	C_prodParamDatasDF[3] = prodParamDatas[3] = "ADV"; // resetTiming
	C_prodParamDatasDF[4] = prodParamDatas[4] = "ARR"; // payTiming
	C_prodParamDatasDF[5] = prodParamDatas[5] = "ADJ"; // intRule
	C_prodParamDatasDF[6] = prodParamDatas[6] = "SS";  // stubRule
	C_prodParamDatasDF[7] = prodParamDatas[7] = "10000.0";// resetGap

	CCString C_daycount;
	long dayCountId;

	CCString C_resetFreq;
	long resetFreqId;

	CCString C_payFreq;
	long payFreqId;

	CCString C_resetTiming;
	long resetTimingId;
	
	CCString C_payTiming;
	long payTimingId;

	CCString C_intRule;
    long intRuleId;
	
	CCString C_stubRule;
    long stubRuleId;

    double C_resetGap;

	CCString C_ccy;
    long ccyId;

   	double C_fixing1_double;
   	VECTOR<double> C_fixing1_double_vect;
	CCString C_fixing1_str;
	long fixing1_type;
	VECTOR<double> C_fixing1_default(0.0);

   	double C_fixing2_double;
   	VECTOR<double> C_fixing2_double_vect;
	CCString C_fixing2_str;
	long fixing2_type;
	VECTOR<double> C_fixing2_default(0.0);

   	double C_fixing3_double;
   	VECTOR<double> C_fixing3_double_vect;
	CCString C_fixing3_str;
	long fixing3_type;
	VECTOR<double> C_fixing3_default(0.0);

	double C_freezeFixing;
	double C_freezeFixing_default = 0.0;

	// error
	static int error;
	static char* reason = "";

	XL_readNumCell(XL_startDate,C_startDate," ARM_ERR: start date: date expected",C_result);
	XL_readNumCell(XL_endDate,C_endDate," ARM_ERR: end date: date expected",C_result);
	XL_readStrCell(XL_capOrFloor,C_capOrFloor," ARM_ERR: cap or floor: string expected",C_result);
	XL_readStrOrNumCell(XL_strike,C_strike_str,C_strike_double,strike_type," ARM_ERR: strike: string or numeric expected",C_result);
	XL_readNumVector(XL_spreads,C_spreads," ARM_ERR: spread vector: array of numeric expected", C_result);

	XL_readStrVectorWD(XL_payIndexParams,C_payIndexParams,C_payIndexParamsDF," ARM_ERR: payIndex Params: array expected",DOUBLE_TYPE,C_result);

	XL_readStrVectorWD(XL_spreadIdxTypes,C_liborTypes,C_liborTypesDF," ARM_ERR: spread indexes : array of string expected",DOUBLE_TYPE,C_result);

	XL_readNumVectorWD(XL_sprdWeights_slope,C_weights_slope,C_weights_slope_df," ARM_ERR: weights_slopeflag: array of numeric expected", C_result);

	// paramDatas
	XL_readStrVectorWD(XL_prodParamDatas,C_prodParamDatas,C_prodParamDatasDF," ARM_ERR: daycount,resetFreq,payFreq,resetTiming,payTiming,resetGap,intRule,stubRule: array of string expected",DOUBLE_TYPE,C_result);
	
	XL_readStrCellWD(XL_currency,C_ccy,"DEFAULT"," ARM_ERR: currency : object expected",C_result);
	XL_readStrOrNumVectorWD(XL_fixing1,C_fixing1_str,C_fixing1_double_vect,C_fixing1_default,fixing1_type," ARM_ERR: fixing taux : array or string expected",C_result);
	XL_readStrOrNumVectorWD(XL_fixing2,C_fixing2_str,C_fixing2_double_vect,C_fixing2_default,fixing2_type," ARM_ERR: fixing taux : array or string expected",C_result);
	XL_readStrOrNumVectorWD(XL_fixing3,C_fixing3_str,C_fixing3_double_vect,C_fixing3_default,fixing3_type," ARM_ERR: fixing taux : array or string expected",C_result);
	XL_readNumCellWD(XL_freezeFixing,C_freezeFixing,C_freezeFixing_default," ARM_ERR: freeze fixing : numeric expected",C_result);

	//traitement de spreads
	if(C_spreads.size()<2)
	{
		C_result.setMsg ("ARM_ERR: Verify the Vector of Spread");
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}
	C_spread1 = C_spreads[0];
	C_spread2 = C_spreads[1];

	// traitement de payIndexType_itsmargin_itsWeight
	int paySize = C_payIndexParams.size() < 3 ? C_payIndexParams.size() : 3;
	for(int i=0; i< paySize; i++)
	{
		if( strcmp(C_payIndexParams[i], "DEFAULT") != 0)// fill in the value not by default
			payIndexParams[i] = C_payIndexParams[i];
	}

	C_payIndex = payIndexParams[0];
	C_payMargin_str = payIndexParams[1];
	C_payIdxWeight = atof(payIndexParams[2]);
//	C_payIdxFixedRate_str = payIndexParams[3];
	
	if(C_payIndex == "DEFAULT")
	{
		payIndexId = ARM_NULL_OBJECT;
	}
	else
	{
		payIndexId = LocalGetNumObjectId(C_payIndex);
	}

	if(C_payMargin_str == "DEFAULT")
	{
		payMarginId = 0.0;
		payMargin_type = 0;//nombre
	}
	else
	{
		payMarginId = LocalGetNumObjectId(C_payMargin_str);

		if( (int)payMarginId == -1 ) // non objet
		{
			payMarginId = atof(C_payMargin_str);
			payMargin_type = 0;
		}
		else
			payMargin_type = 1;//objet
	}

	//traitement de indexes
	int spreadIdxSize = C_liborTypes.size() < 2 ? C_liborTypes.size() : 2;
	for (i =0; i< spreadIdxSize; i++)
	{
		if(strcmp(C_liborTypes[i],"") != 0)
			liborTypes[i] = C_liborTypes[i];
	}
	C_liborType1 = liborTypes[0];
	C_liborType2 = liborTypes[1];

	// traitement de weights de spread Indexes
	if(C_weights_slope.size()<2)
	{
		C_result.setMsg ("ARM_ERR: Verify the Vector of Weight");
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}
	
	C_weight1 = C_weights_slope[0];
	C_weight2 = C_weights_slope[1];

	if(C_weights_slope.size() == 2)
	{
		C_slopeFlag = C_slopeFlagDefault;
		C_cptStrikeMethod = C_cptStrikeMethod_default;
		C_computedFormula = C_computedFormula_default;
	}
	else if(C_weights_slope.size() == 3)
	{
		C_slopeFlag = C_weights_slope[2];
		C_cptStrikeMethod = C_cptStrikeMethod_default;
		C_computedFormula = C_computedFormula_default;
	}
	else if(C_weights_slope.size() == 4)
	{
		C_slopeFlag = C_weights_slope[2];
		C_cptStrikeMethod = C_weights_slope[3];
		C_computedFormula = C_computedFormula_default;
	}
	else
	{
		C_slopeFlag = C_weights_slope[2];
		C_cptStrikeMethod = C_weights_slope[3];
		C_computedFormula = C_weights_slope[4];
	}

	// traitement de prodParamDatas
	int prodParamSize = C_prodParamDatas.size() < 8 ? C_prodParamDatas.size() : 8;
	for (i =0; i<prodParamSize; i++)
	{
		if( strcmp(C_prodParamDatas[i],"DEFAULT")!= 0 )
			prodParamDatas[i] = C_prodParamDatas[i];
	}
	// dayCount
	if (prodParamDatas.size() >= 1)
		C_daycount    = prodParamDatas[0];
	else
		C_daycount    = C_prodParamDatasDF[0];
	// resetFreq
	if (prodParamDatas.size() >= 2)
		C_resetFreq   = prodParamDatas[1];
	else
		C_resetFreq   = C_prodParamDatasDF[1];
	// payFreq
	if (prodParamDatas.size() >= 3)
		C_payFreq     = prodParamDatas[2];
	else
		C_payFreq     = C_prodParamDatasDF[2];
	// resetTiming
	if (prodParamDatas.size() >= 4)
		C_resetTiming = prodParamDatas[3];
	else
		C_resetTiming = C_prodParamDatasDF[3];
	// payTiming
	if (prodParamDatas.size() >= 5)
		C_payTiming   = prodParamDatas[4];
	else
		C_payTiming   = C_prodParamDatasDF[4];
	// intRule
	if (prodParamDatas.size() >= 6)
		C_intRule     = prodParamDatas[5];
	else
		C_intRule     = C_prodParamDatasDF[5];
	// stubRule
	if (prodParamDatas.size() >= 7)
		C_stubRule    = prodParamDatas[6];
	else
		C_stubRule    = C_prodParamDatasDF[6];
	// resetGap
	if (prodParamDatas.size() >= 8)
		C_resetGap    = atof(prodParamDatas[7]);
	else
		C_resetGap    = atof(C_prodParamDatasDF[7]);

	// traitement de Fixing 
	if ((XL_fixing1->xltype == xltypeMissing) || (XL_fixing1->xltype == xltypeNil))
	{
		fixing1_type = 1;
		C_fixing1_double = ARM_NULL_OBJECT;
    }
	else if ( fixing1_type == XL_TYPE_STRING )
	{
		C_fixing1_double = (double) LocalGetNumObjectId(C_fixing1_str);
		fixing1_type = 1;
	}
	else
	{
	    fixing1_type = 0;
	}

	if ((XL_fixing2->xltype == xltypeMissing) || (XL_fixing2->xltype == xltypeNil))
	{
		fixing2_type = 1;
		C_fixing2_double = ARM_NULL_OBJECT;
    }
	else if ( fixing2_type == XL_TYPE_STRING )
	{
		C_fixing2_double = (double) LocalGetNumObjectId(C_fixing2_str);
		fixing2_type = 1;
	}
	else
	{
	    fixing2_type = 0;
	}

	if ((XL_fixing3->xltype == xltypeMissing) || (XL_fixing3->xltype == xltypeNil))
	{
		fixing3_type = 1;
		C_fixing3_double = ARM_NULL_OBJECT;
    }
	else if ( fixing3_type == XL_TYPE_STRING )
	{
		C_fixing3_double = (double) LocalGetNumObjectId(C_fixing3_str);
		fixing3_type = 1;
	}
	else
	{
	    fixing3_type = 0;
	}

	if((capOrFloorId = ARM_ConvCapOrFloor (C_capOrFloor, C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	//payoffLiborTypeId = ARM_ConvIrType (C_payoffLiborType);

	liborType1Id = ARM_ConvIrType (C_liborType1);

	liborType2Id = ARM_ConvIrType (C_liborType2);

	dayCountId = ARM_ConvDayCount (C_daycount);
	
	if((resetFreqId = ARM_ConvFrequency (C_resetFreq, C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	if((payFreqId = ARM_ConvFrequency (C_payFreq, C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	resetTimingId = ARM_ConvPayResetRule (C_resetTiming);

	payTimingId = ARM_ConvPayResetRule (C_payTiming);

	if(C_ccy == "DEFAULT")
	{
		ccyId = ARM_NULL_OBJECT;
	}
	else
	{
		ccyId = LocalGetNumObjectId (C_ccy);
	}

	if ( strike_type == XL_TYPE_STRING )
	{
	   C_strike_double = (double) LocalGetNumObjectId(C_strike_str);

	   strike_type = 1;
	}
	else
	{
	   strike_type = 0;
	}
	
	intRuleId = ARM_ConvIntRule(C_intRule);

	stubRuleId = ARM_ConvStubRule(C_stubRule);

	long retCode;
	long objId;
	CCString prevClass;
	
	CCString curClass = LOCAL_SPREAD_OPTION_CLASS;
	CCString stringId = GetLastCurCellEnvValue ();

	if(!stringId)
	{
		retCode = ARMLOCAL_SPREADCORRIDOR(C_startDate,
										  C_endDate,
										  capOrFloorId,
										  (long) strike_type,
										  C_strike_double,
										  payIndexId,//(long)payoffLiborTypeId,
 										  (long)liborType1Id,
 										  (long)liborType2Id,
										  C_weight1,
										  C_weight2,
										  dayCountId,
										  resetFreqId,
										  payFreqId,
										  resetTimingId,
										  payTimingId,
										  ccyId,
										  long(C_resetGap),
										  C_spread1,
										  C_spread2,
										  (long) fixing1_type,
										  C_fixing1_double,
										  C_fixing1_double_vect,
										  (long) fixing2_type,
										  C_fixing2_double,
										  C_fixing2_double_vect,
										  intRuleId,
										  stubRuleId,
										  C_slopeFlag,
										  (int) C_cptStrikeMethod,
										  payMargin_type,
										  payMarginId,
										  C_payIdxWeight,
										  (long) fixing3_type,
										  C_fixing3_double,
										  C_fixing3_double_vect,
										  (int)C_freezeFixing,
										  (int) C_computedFormula,
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
			retCode = ARMLOCAL_SPREADCORRIDOR(C_startDate,
											  C_endDate,
											  capOrFloorId,
											  (long) strike_type,
											  C_strike_double,
											  payIndexId,//(long)payoffLiborTypeId,
											  (long)liborType1Id,
											  (long)liborType2Id,
											  C_weight1,
											  C_weight2,
											  dayCountId,
											  resetFreqId,
											  payFreqId,
											  resetTimingId,
											  payTimingId,
											  ccyId,
											  long(C_resetGap),
											  C_spread1,
											  C_spread2,
											  (long) fixing1_type,
											  C_fixing1_double,
											  C_fixing1_double_vect,
											  (long) fixing2_type,
											  C_fixing2_double,
											  C_fixing2_double_vect,
											  intRuleId,
											  stubRuleId,
											  C_slopeFlag,
											  (int) C_cptStrikeMethod,
											  payMargin_type,
											  payMarginId,
											  C_payIdxWeight,
	 										  (long) fixing3_type,
											  C_fixing3_double,
											  C_fixing3_double_vect,
											  (int)C_freezeFixing,
											  (int) C_computedFormula,
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
			retCode = ARMLOCAL_SPREADCORRIDOR(C_startDate,
											 C_endDate,
											 capOrFloorId,
											 (long) strike_type,
											 C_strike_double,
										     payIndexId,//(long)payoffLiborTypeId,
											 (long)liborType1Id,
											 (long)liborType2Id,
											 C_weight1,
											 C_weight2,
											 dayCountId,
											 resetFreqId,
											 payFreqId,
											 resetTimingId,
											 payTimingId,
											 ccyId,
                                             long(C_resetGap),
											 C_spread1,
											 C_spread2,
											 (long) fixing1_type,
											 C_fixing1_double,
											 C_fixing1_double_vect,
											 (long) fixing2_type,
											 C_fixing2_double,
											 C_fixing2_double_vect,
											 intRuleId,
											 stubRuleId,
											 C_slopeFlag,
											 (int) C_cptStrikeMethod,
										     payMargin_type,
										     payMarginId,
											 C_payIdxWeight,
											 (long) fixing3_type,
											 C_fixing3_double,
											 C_fixing3_double_vect,
											 (int)C_freezeFixing,
											 (int) C_computedFormula,
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_ARM_SPREADCORRIDOR" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}

__declspec(dllexport) LPXLOPER WINAPI Local_PXL_ARM_SPREADCORRIDOR( LPXLOPER XL_startDate,
																	LPXLOPER XL_endDate,
 																	LPXLOPER XL_capOrFloor,
																	LPXLOPER XL_strike,
																	LPXLOPER XL_spreads,//Vect: spread1, spread2
																	LPXLOPER XL_payIndexParams,//Vect:payIndexObj, itsWeight, itsMargin
																	LPXLOPER XL_spreadIdxTypes,//Vect: IndexType1, IndexType2
																	LPXLOPER XL_sprdWeights_slope,//Vect: Idx Weight1,Weight2, slopeflag and cptStrikeMethod
																	LPXLOPER XL_prodParamDatas,//Vect: dayCount, resetFreq, payFreq, resetTiming, payTiming, intRule, stubRule,resetGap 
																	LPXLOPER XL_currency,
																	LPXLOPER XL_fixing1,
																	LPXLOPER XL_fixing2,
																	LPXLOPER XL_fixing3,
																	LPXLOPER XL_freezeFixing)
{
	ADD_LOG("Local_PXL_ARM_SPREADCORRIDOR");
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

	CCString C_capOrFloor;
	long capOrFloorId;

	double C_strike_double;
	CCString C_strike_str;
	long strike_type;
	
	double C_spread1;
	double C_spread2;
	VECTOR<double> C_spreads;

	CCString C_payIndex;
	long payIndexId;
	
	CCString C_payMargin_str;
	long payMargin_type;
	double payMarginId;

	double C_payIdxWeight;
	
	VECTOR<CCString> payIndexParams(3);
   	VECTOR<CCString> C_payIndexParams;
	VECTOR<CCString> C_payIndexParamsDF(3);
	C_payIndexParamsDF[0] = payIndexParams[0] = "DEFAULT";// pay index
	C_payIndexParamsDF[1] = payIndexParams[1] = "DEFAULT";// pay margin referenceValue
	C_payIndexParamsDF[2] = payIndexParams[2] = "1";// pay weight
	
	CCString C_liborType1;
	long liborType1Id;
	CCString C_liborType2;
	long liborType2Id;
	VECTOR<CCString> liborTypes(2);
   	VECTOR<CCString> C_liborTypes;
	VECTOR<CCString> C_liborTypesDF(2);
	C_liborTypesDF[0] = liborTypes[0] = "LIBOR3M";
	C_liborTypesDF[1] = liborTypes[1] = "LIBOR3M";

	double C_weight1, C_weight2, C_cptStrikeMethod;
	VECTOR<double> weights_slope(2,1.0);
	VECTOR<double> C_weights_slope;
	VECTOR<double> C_weights_slope_df(2, 1.0);

	long C_slopeFlag;
	long C_slopeFlagDefault = 1;
	double C_cptStrikeMethod_default = 1.0;

	double C_computedFormula;
	double C_computedFormula_default = 1.0;

	// prodParamDatas
	VECTOR<CCString> prodParamDatas(8);
	VECTOR<CCString> C_prodParamDatas;
	VECTOR<CCString> C_prodParamDatasDF(8);
	C_prodParamDatasDF[0] = prodParamDatas[0] = "30/360"; // dayCount
	C_prodParamDatasDF[1] = prodParamDatas[1] = "A"; // resetFreq
	C_prodParamDatasDF[2] = prodParamDatas[2] = "A"; // payFreq
	C_prodParamDatasDF[3] = prodParamDatas[3] = "ADV"; // resetTiming
	C_prodParamDatasDF[4] = prodParamDatas[4] = "ARR"; // payTiming
	C_prodParamDatasDF[5] = prodParamDatas[5] = "ADJ"; // intRule
	C_prodParamDatasDF[6] = prodParamDatas[6] = "SS";  // stubRule
	C_prodParamDatasDF[7] = prodParamDatas[7] = "10000.0";// resetGap

	CCString C_daycount;
	long dayCountId;

	CCString C_resetFreq;
	long resetFreqId;

	CCString C_payFreq;
	long payFreqId;

	CCString C_resetTiming;
	long resetTimingId;
	
	CCString C_payTiming;
	long payTimingId;

	CCString C_intRule;
    long intRuleId;
	
	CCString C_stubRule;
    long stubRuleId;

    double C_resetGap;

	CCString C_ccy;
    long ccyId;

   	double C_fixing1_double;
   	VECTOR<double> C_fixing1_double_vect;
	CCString C_fixing1_str;
	long fixing1_type;
	VECTOR<double> C_fixing1_default(0.0);

   	double C_fixing2_double;
   	VECTOR<double> C_fixing2_double_vect;
	CCString C_fixing2_str;
	long fixing2_type;
	VECTOR<double> C_fixing2_default(0.0);

   	double C_fixing3_double;
   	VECTOR<double> C_fixing3_double_vect;
	CCString C_fixing3_str;
	long fixing3_type;
	VECTOR<double> C_fixing3_default(0.0);

	double C_freezeFixing;
	double C_freezeFixing_default = 0.0;

	// error
	static int error;
	static char* reason = "";

	XL_readNumCell(XL_startDate,C_startDate," ARM_ERR: start date: date expected",C_result);
	XL_readNumCell(XL_endDate,C_endDate," ARM_ERR: end date: date expected",C_result);
	XL_readStrCell(XL_capOrFloor,C_capOrFloor," ARM_ERR: cap or floor: string expected",C_result);
	XL_readStrOrNumCell(XL_strike,C_strike_str,C_strike_double,strike_type," ARM_ERR: strike: string or numeric expected",C_result);
	XL_readNumVector(XL_spreads,C_spreads," ARM_ERR: spread vector: array of numeric expected", C_result);

	XL_readStrVectorWD(XL_payIndexParams,C_payIndexParams,C_payIndexParamsDF," ARM_ERR: payIndex Params: array expected",DOUBLE_TYPE,C_result);

	XL_readStrVectorWD(XL_spreadIdxTypes,C_liborTypes,C_liborTypesDF," ARM_ERR: spread indexes : array of string expected",DOUBLE_TYPE,C_result);

	XL_readNumVectorWD(XL_sprdWeights_slope,C_weights_slope,C_weights_slope_df," ARM_ERR: weights_slopeflag: array of numeric expected", C_result);

	// paramDatas
	XL_readStrVectorWD(XL_prodParamDatas,C_prodParamDatas,C_prodParamDatasDF," ARM_ERR: daycount,resetFreq,payFreq,resetTiming,payTiming,resetGap,intRule,stubRule: array of string expected",DOUBLE_TYPE,C_result);
	
	XL_readStrCellWD(XL_currency,C_ccy,"DEFAULT"," ARM_ERR: currency : object expected",C_result);
	XL_readStrOrNumVectorWD(XL_fixing1,C_fixing1_str,C_fixing1_double_vect,C_fixing1_default,fixing1_type," ARM_ERR: fixing taux : array or string expected",C_result);
	XL_readStrOrNumVectorWD(XL_fixing2,C_fixing2_str,C_fixing2_double_vect,C_fixing2_default,fixing2_type," ARM_ERR: fixing taux : array or string expected",C_result);
	XL_readStrOrNumVectorWD(XL_fixing3,C_fixing3_str,C_fixing3_double_vect,C_fixing3_default,fixing3_type," ARM_ERR: fixing taux : array or string expected",C_result);
	XL_readNumCellWD(XL_freezeFixing,C_freezeFixing,C_freezeFixing_default," ARM_ERR: freeze fixing : numeric expected",C_result);

	//traitement de spreads
	if(C_spreads.size()<2)
	{
		C_result.setMsg ("ARM_ERR: Verify the Vector of Spread");
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}
	C_spread1 = C_spreads[0];
	C_spread2 = C_spreads[1];

	// traitement de payIndexType_itsWeight
	int paySize = C_payIndexParams.size() < 3 ? C_payIndexParams.size() : 3;
	for(int i=0; i<paySize; i++)
	{
		if( strcmp(C_payIndexParams[i], "DEFAULT") != 0)// fill in the value not by default
			payIndexParams[i] = C_payIndexParams[i];
	}

	C_payIndex = payIndexParams[0];
	C_payMargin_str = payIndexParams[1];
	C_payIdxWeight = atof(payIndexParams[2]);

	if(C_payIndex == "DEFAULT")
	{
		payIndexId = ARM_NULL_OBJECT;
	}
	else
	{
		payIndexId = LocalGetNumObjectId(C_payIndex);
	}

	if(C_payMargin_str == "DEFAULT")
	{
		payMarginId = 0.0;
		payMargin_type = 0;//nombre
	}
	else
	{
		payMarginId = LocalGetNumObjectId(C_payMargin_str);

		if( (int)payMarginId == -1 ) // non objet
		{
			payMarginId = atof(C_payMargin_str);
			payMargin_type = 0;
		}
		else
			payMargin_type = 1;//objet
	}


	//traitement de indexes
	int spreadIdxSize = C_liborTypes.size() < 2 ? C_liborTypes.size() : 2 ;
	for (i =0; i< spreadIdxSize; i++)
	{
		if(strcmp(C_liborTypes[i],"") != 0)
			liborTypes[i] = C_liborTypes[i];
	}
	C_liborType1 = liborTypes[0];
	C_liborType2 = liborTypes[1];

	// traitement de weights de spread Indexes
	if(C_weights_slope.size()<2)
	{
		C_result.setMsg ("ARM_ERR: Verify the Vector of Weight");
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}
	
	C_weight1 = C_weights_slope[0];
	C_weight2 = C_weights_slope[1];

	if(C_weights_slope.size() == 2)
	{
		C_slopeFlag = C_slopeFlagDefault;
		C_cptStrikeMethod = C_cptStrikeMethod_default;
		C_computedFormula = C_computedFormula_default;
	}
	else if(C_weights_slope.size() == 3)
	{
		C_slopeFlag = C_weights_slope[2];
		C_cptStrikeMethod = C_cptStrikeMethod_default;
		C_computedFormula = C_computedFormula_default;
	}
	else if(C_weights_slope.size() == 4)
	{
		C_slopeFlag = C_weights_slope[2];
		C_cptStrikeMethod = C_weights_slope[3];
		C_computedFormula = C_computedFormula_default;
	}
	else
	{
		C_slopeFlag = C_weights_slope[2];
		C_cptStrikeMethod = C_weights_slope[3];
		C_computedFormula = C_weights_slope[4];
	}

	// traitement de prodParamDatas
	int prodSize = C_prodParamDatas.size() < 8 ? C_prodParamDatas.size() : 8;
	for (i =0; i<prodSize; i++)
	{
		if( strcmp(C_prodParamDatas[i],"DEFAULT")!= 0 )
			prodParamDatas[i] = C_prodParamDatas[i];
	}
	// dayCount
	if (prodParamDatas.size() >= 1)
		C_daycount    = prodParamDatas[0];
	else
		C_daycount    = C_prodParamDatasDF[0];
	// resetFreq
	if (prodParamDatas.size() >= 2)
		C_resetFreq   = prodParamDatas[1];
	else
		C_resetFreq   = C_prodParamDatasDF[1];
	// payFreq
	if (prodParamDatas.size() >= 3)
		C_payFreq     = prodParamDatas[2];
	else
		C_payFreq     = C_prodParamDatasDF[2];
	// resetTiming
	if (prodParamDatas.size() >= 4)
		C_resetTiming = prodParamDatas[3];
	else
		C_resetTiming = C_prodParamDatasDF[3];
	// payTiming
	if (prodParamDatas.size() >= 5)
		C_payTiming   = prodParamDatas[4];
	else
		C_payTiming   = C_prodParamDatasDF[4];
	// intRule
	if (prodParamDatas.size() >= 6)
		C_intRule     = prodParamDatas[5];
	else
		C_intRule     = C_prodParamDatasDF[5];
	// stubRule
	if (prodParamDatas.size() >= 7)
		C_stubRule    = prodParamDatas[6];
	else
		C_stubRule    = C_prodParamDatasDF[6];
	// resetGap
	if (prodParamDatas.size() >= 8)
		C_resetGap    = atof(prodParamDatas[7]);
	else
		C_resetGap    = atof(C_prodParamDatasDF[7]);

	// traitement de Fixing 
	if ((XL_fixing1->xltype == xltypeMissing) || (XL_fixing1->xltype == xltypeNil))
	{
		fixing1_type = 1;
		C_fixing1_double = ARM_NULL_OBJECT;
    }
	else if ( fixing1_type == XL_TYPE_STRING )
	{
		C_fixing1_double = (double) LocalGetNumObjectId(C_fixing1_str);
		fixing1_type = 1;
	}
	else
	{
	    fixing1_type = 0;
	}

	if ((XL_fixing2->xltype == xltypeMissing) || (XL_fixing2->xltype == xltypeNil))
	{
		fixing2_type = 1;
		C_fixing2_double = ARM_NULL_OBJECT;
    }
	else if ( fixing2_type == XL_TYPE_STRING )
	{
		C_fixing2_double = (double) LocalGetNumObjectId(C_fixing2_str);
		fixing2_type = 1;
	}
	else
	{
	    fixing2_type = 0;
	}

	if ((XL_fixing3->xltype == xltypeMissing) || (XL_fixing3->xltype == xltypeNil))
	{
		fixing3_type = 1;
		C_fixing3_double = ARM_NULL_OBJECT;
    }
	else if ( fixing3_type == XL_TYPE_STRING )
	{
		C_fixing3_double = (double) LocalGetNumObjectId(C_fixing3_str);
		fixing3_type = 1;
	}
	else
	{
	    fixing3_type = 0;
	}

	if((capOrFloorId = ARM_ConvCapOrFloor (C_capOrFloor, C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	liborType1Id = ARM_ConvIrType (C_liborType1);

	liborType2Id = ARM_ConvIrType (C_liborType2);

	dayCountId = ARM_ConvDayCount (C_daycount);
	
	if((resetFreqId = ARM_ConvFrequency (C_resetFreq, C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	if((payFreqId = ARM_ConvFrequency (C_payFreq, C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	resetTimingId = ARM_ConvPayResetRule (C_resetTiming);

	payTimingId = ARM_ConvPayResetRule (C_payTiming);

	if(C_ccy == "DEFAULT")
	{
		ccyId = ARM_NULL_OBJECT;
	}
	else
	{
		ccyId = LocalGetNumObjectId (C_ccy);
	}

	if ( strike_type == XL_TYPE_STRING )
	{
	   C_strike_double = (double) LocalGetNumObjectId(C_strike_str);

	   strike_type = 1;
	}
	else
	{
	   strike_type = 0;
	}
	
	intRuleId = ARM_ConvIntRule(C_intRule);

	stubRuleId = ARM_ConvStubRule(C_stubRule);

	long retCode;
	long objId;
	
	CCString curClass = LOCAL_SPREAD_OPTION_CLASS;
	CCString stringId = GetLastCurCellEnvValue ();

	retCode = ARMLOCAL_SPREADCORRIDOR(C_startDate,
									  C_endDate,
									  capOrFloorId,
									  (long) strike_type,
									  C_strike_double,
									  payIndexId,
									  (long)liborType1Id,
									  (long)liborType2Id,
									  C_weight1,
									  C_weight2,
									  dayCountId,
									  resetFreqId,
									  payFreqId,
									  resetTimingId,
									  payTimingId,
									  ccyId,
									  long(C_resetGap),
									  C_spread1,
									  C_spread2,
									  (long) fixing1_type,
									  C_fixing1_double,
									  C_fixing1_double_vect,
									  (long) fixing2_type,
									  C_fixing2_double,
									  C_fixing2_double_vect,
									  intRuleId,
									  stubRuleId,
									  C_slopeFlag,
									  (int) C_cptStrikeMethod,
									  (long) payMargin_type,
									  payMarginId,
									  C_payIdxWeight,
									  (long) fixing3_type,
									  C_fixing3_double,
									  C_fixing3_double_vect,
									  (int)C_freezeFixing,
									  (int) C_computedFormula,
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_PXL_ARM_SPREADCORRIDOR"  )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}

//------------------------------------------------------------------------------------------//
// Interface function for VMS spread options corridor										//
//------------------------------------------------------------------------------------------//
__declspec(dllexport) LPXLOPER WINAPI Local_ARM_SPREADCORRIDORVMS(LPXLOPER XL_startDate,
																  LPXLOPER XL_endDate,
 																  LPXLOPER XL_capOrFloor,
																  LPXLOPER XL_strike,
																  LPXLOPER XL_spreads,				//Vect: spread1, spread2
																  LPXLOPER XL_payIndexParams,		//Vect: payIndexObj, itsMargin,itsWeight
															      LPXLOPER XL_CMSIndexes1,		
																  LPXLOPER XL_CMSIndexes2,		
															      LPXLOPER XL_sprdWeights_slope,	//Vect: Idx Weight1,Weight2, slopeflag and cptStrikeMethod
															      LPXLOPER XL_prodParamDatas,		//Vect: dayCount, resetFreq, payFreq, resetTiming, payTiming, intRule, stubRule,resetGap 
															      LPXLOPER XL_currency,
															      LPXLOPER XL_fixing1,
															      LPXLOPER XL_fixing2,
															      LPXLOPER XL_fixingPay,
															      LPXLOPER XL_freezeFixing)
{
	ADD_LOG("Local_ARM_SPREADCORRIDORVMS");
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

	CCString C_capOrFloor;
	long capOrFloorId;

	double C_strike_double;
	CCString C_strike_str;
	long strike_type;
	
	double C_spread1;
	double C_spread2;
	VECTOR<double> C_spreads;

	CCString C_payIndex;
	long payIndexId;

	CCString C_payMargin_str;
	long payMargin_type;
	double payMarginId;

	double C_payIdxWeight;

	CCString C_payIdxFixedRate_str;
		
	VECTOR<CCString> payIndexParams(3);
   	VECTOR<CCString> C_payIndexParams;
	VECTOR<CCString> C_payIndexParamsDF(3);
	C_payIndexParamsDF[0] = payIndexParams[0] = "DEFAULT";	// pay index objet
	C_payIndexParamsDF[1] = payIndexParams[1] = "DEFAULT";  // Margin
	C_payIndexParamsDF[2] = payIndexParams[2] = "1";		// weight

	CCString C_CMSIndexes1;
	double   C_CMSIndexes1_double;

	CCString C_CMSIndexes2;
	double   C_CMSIndexes2_double;

	double C_weight1, C_weight2, C_cptStrikeMethod;
	VECTOR<double> weights_slope(2,1.0);
	VECTOR<double> C_weights_slope;
	VECTOR<double> C_weights_slope_df(2, 1.0);

	long C_slopeFlag;
	long C_slopeFlagDefault = 1;

	double C_cptStrikeMethod_default = 1.0;

	double C_computedFormula;
	double C_computedFormula_default = 1.0;

	// prodParamDatas
	VECTOR<CCString> prodParamDatas(8);
	VECTOR<CCString> C_prodParamDatas;
	VECTOR<CCString> C_prodParamDatasDF(8);
	C_prodParamDatasDF[0] = prodParamDatas[0] = "30/360";	// dayCount
	C_prodParamDatasDF[1] = prodParamDatas[1] = "A";		// resetFreq
	C_prodParamDatasDF[2] = prodParamDatas[2] = "A";		// payFreq
	C_prodParamDatasDF[3] = prodParamDatas[3] = "ADV";		// resetTiming
	C_prodParamDatasDF[4] = prodParamDatas[4] = "ARR";		// payTiming
	C_prodParamDatasDF[5] = prodParamDatas[5] = "ADJ";		// intRule
	C_prodParamDatasDF[6] = prodParamDatas[6] = "SS";		// stubRule
	C_prodParamDatasDF[7] = prodParamDatas[7] = "10000.0";	// resetGap

	CCString C_daycount;
	long dayCountId;

	CCString C_resetFreq;
	long resetFreqId;

	CCString C_payFreq;
	long payFreqId;

	CCString C_resetTiming;
	long resetTimingId;
	
	CCString C_payTiming;
	long payTimingId;

	CCString C_intRule;
    long intRuleId;
	
	CCString C_stubRule;
    long stubRuleId;

    double C_resetGap;

	CCString C_ccy;
    long ccyId;

	CCString C_fixing1;
	double   C_fixing1_double;

	CCString C_fixing2;
	double   C_fixing2_double;

	CCString C_fixingPay;
	double   C_fixingPay_double;

	double C_freezeFixing;
	double C_freezeFixing_default = 0.0;

	// error
	static int error;
	static char* reason = "";


	// Parameters reading
	//----------------------------------------------

	XL_readNumCell(XL_startDate,C_startDate," ARM_ERR: start date: date expected",C_result);

	XL_readNumCell(XL_endDate,C_endDate," ARM_ERR: end date: date expected",C_result);
	
	XL_readStrCell(XL_capOrFloor,C_capOrFloor," ARM_ERR: cap or floor: string expected",C_result);
	
	XL_readStrOrNumCell(XL_strike,C_strike_str,C_strike_double,strike_type," ARM_ERR: strike: string or numeric expected",C_result);
	
	XL_readNumVector(XL_spreads,C_spreads," ARM_ERR: spread vector: array of numeric expected", C_result);
	
	XL_readStrVectorWD(XL_payIndexParams,C_payIndexParams,C_payIndexParamsDF," ARM_ERR: payIndex Params: array expected",DOUBLE_TYPE,C_result);	
	
	XL_readStrCellWD(XL_CMSIndexes1, C_CMSIndexes1, "", "ARM_ERR: CMS indexes : string expected", C_result);
	
	XL_readStrCellWD(XL_CMSIndexes2, C_CMSIndexes2, "", "ARM_ERR: CMS indexes : string expected", C_result);
	
	XL_readNumVectorWD(XL_sprdWeights_slope,C_weights_slope,C_weights_slope_df," ARM_ERR: weights_slopeflag: array of numeric expected", C_result);
	
	XL_readStrVectorWD(XL_prodParamDatas,C_prodParamDatas,C_prodParamDatasDF," ARM_ERR: daycount,resetFreq,payFreq,resetTiming,payTiming,resetGap,intRule,stubRule: array of string expected",DOUBLE_TYPE,C_result);
	
	XL_readStrCellWD(XL_currency,C_ccy,"DEFAULT"," ARM_ERR: currency : object expected",C_result);
	
	XL_readStrCellWD(XL_fixing1,C_fixing1,""," ARM_ERR: CMS fixing : string expected",C_result);
	
	XL_readStrCellWD(XL_fixing2,C_fixing2,""," ARM_ERR: CMS fixing : string expected",C_result);
	
	XL_readStrCellWD(XL_fixingPay,C_fixingPay,""," ARM_ERR: Pay fixing : string expected",C_result);
	
	XL_readNumCellWD(XL_freezeFixing,C_freezeFixing,C_freezeFixing_default," ARM_ERR: freeze fixing : numeric expected",C_result);

	
	// Spreads treatment
	//----------------------------------------------

	if(C_spreads.size()<2)
	{
		C_result.setMsg ("ARM_ERR: Verify the Vector of Spread");
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	C_spread1 = C_spreads[0];
	C_spread2 = C_spreads[1];


	// Pay index type treatment
	//----------------------------------------------

	int paySize = C_payIndexParams.size() < 3 ? C_payIndexParams.size() : 3;

	for(int i=0; i< paySize; i++)
	{
		if( strcmp(C_payIndexParams[i], "DEFAULT") != 0)	// fill in the value not by default
			payIndexParams[i] = C_payIndexParams[i];
	}

	C_payIndex		= payIndexParams[0];
	C_payMargin_str = payIndexParams[1];
	C_payIdxWeight	= atof(payIndexParams[2]);
	
	if(C_payIndex == "DEFAULT")
	{
		payIndexId = ARM_NULL_OBJECT;
	}
	else
	{
		payIndexId = LocalGetNumObjectId(C_payIndex);
	}

	if(C_payMargin_str == "DEFAULT")
	{
		payMarginId = 0.0;
		payMargin_type = 0;				// number
	}
	else
	{
		payMarginId = LocalGetNumObjectId(C_payMargin_str);

		if( (int)payMarginId == -1 )	// no object
		{
			payMarginId = atof(C_payMargin_str);
			payMargin_type = 0;
		}
		else
			payMargin_type = 1;			// object
	}


	// CMS Indexes treatment
	//----------------------------------------------

	if ( (XL_CMSIndexes1->xltype == xltypeMissing) || (XL_CMSIndexes1->xltype == xltypeNil) )
	{		
		C_CMSIndexes1_double = ARM_NULL_OBJECT;
    }
	else
	{
		C_CMSIndexes1_double = (double) LocalGetNumObjectId(C_CMSIndexes1);
	}



	if ( (XL_CMSIndexes2->xltype == xltypeMissing) || (XL_CMSIndexes2->xltype == xltypeNil) )
	{		
		C_CMSIndexes2_double = ARM_NULL_OBJECT;
    }
	else
	{
		C_CMSIndexes2_double = (double) LocalGetNumObjectId(C_CMSIndexes2);
	}


	// Spread weights treatment
	//----------------------------------------------

	if(C_weights_slope.size()<2)
	{
		C_result.setMsg ("ARM_ERR: Verify the Vector of Weight");
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}
	
	C_weight1 = C_weights_slope[0];
	C_weight2 = C_weights_slope[1];

	if(C_weights_slope.size() == 2)
	{
		C_slopeFlag			= C_slopeFlagDefault;
		C_cptStrikeMethod	= C_cptStrikeMethod_default;
		C_computedFormula	= C_computedFormula_default;
	}
	else if(C_weights_slope.size() == 3)
	{
		C_slopeFlag			= C_weights_slope[2];
		C_cptStrikeMethod	= C_cptStrikeMethod_default;
		C_computedFormula	= C_computedFormula_default;
	}
	else if(C_weights_slope.size() == 4)
	{
		C_slopeFlag			= C_weights_slope[2];
		C_cptStrikeMethod	= C_weights_slope[3];
		C_computedFormula	= C_computedFormula_default;
	}
	else
	{
		C_slopeFlag			= C_weights_slope[2];
		C_cptStrikeMethod	= C_weights_slope[3];
		C_computedFormula	= C_weights_slope[4];
	}


	// prodParamDatas treatment
	//----------------------------------------------

	int prodParamSize = C_prodParamDatas.size() < 8 ? C_prodParamDatas.size() : 8;
	for (i =0; i<prodParamSize; i++)
	{
		if( strcmp(C_prodParamDatas[i],"DEFAULT")!= 0 )
			prodParamDatas[i] = C_prodParamDatas[i];
	}
	// dayCount
	if (prodParamDatas.size() >= 1)
		C_daycount    = prodParamDatas[0];
	else
		C_daycount    = C_prodParamDatasDF[0];
	// resetFreq
	if (prodParamDatas.size() >= 2)
		C_resetFreq   = prodParamDatas[1];
	else
		C_resetFreq   = C_prodParamDatasDF[1];
	// payFreq
	if (prodParamDatas.size() >= 3)
		C_payFreq     = prodParamDatas[2];
	else
		C_payFreq     = C_prodParamDatasDF[2];
	// resetTiming
	if (prodParamDatas.size() >= 4)
		C_resetTiming = prodParamDatas[3];
	else
		C_resetTiming = C_prodParamDatasDF[3];
	// payTiming
	if (prodParamDatas.size() >= 5)
		C_payTiming   = prodParamDatas[4];
	else
		C_payTiming   = C_prodParamDatasDF[4];
	// intRule
	if (prodParamDatas.size() >= 6)
		C_intRule     = prodParamDatas[5];
	else
		C_intRule     = C_prodParamDatasDF[5];
	// stubRule
	if (prodParamDatas.size() >= 7)
		C_stubRule    = prodParamDatas[6];
	else
		C_stubRule    = C_prodParamDatasDF[6];
	// resetGap
	if (prodParamDatas.size() >= 8)
		C_resetGap    = atof(prodParamDatas[7]);
	else
		C_resetGap    = atof(C_prodParamDatasDF[7]);


	// Fixing1 treatment
	//----------------------------------------------

	if ((XL_fixing1->xltype == xltypeMissing) || (XL_fixing1->xltype == xltypeNil))
	{
		C_fixing1_double = ARM_NULL_OBJECT;
    }
	else
	{
		C_fixing1_double = (double) LocalGetNumObjectId(C_fixing1);
	}

	
	// Fixing2 treatment
	//----------------------------------------------
	
	if ((XL_fixing2->xltype == xltypeMissing) || (XL_fixing2->xltype == xltypeNil))
	{
		C_fixing2_double = ARM_NULL_OBJECT;
    }
	else
	{
		C_fixing2_double = (double) LocalGetNumObjectId(C_fixing2);
	}


	// Fixing3 treatment
	//----------------------------------------------
	
	if ((XL_fixingPay->xltype == xltypeMissing) || (XL_fixingPay->xltype == xltypeNil))
	{
		C_fixingPay_double = ARM_NULL_OBJECT;
    }
	else
	{
		C_fixingPay_double = (double) LocalGetNumObjectId(C_fixingPay);

	}


	// Ids recovery
	//----------------------------------------------

	if((capOrFloorId = ARM_ConvCapOrFloor (C_capOrFloor, C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	dayCountId		= ARM_ConvDayCount (C_daycount);
	
	if((resetFreqId = ARM_ConvFrequency (C_resetFreq, C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	if((payFreqId = ARM_ConvFrequency (C_payFreq, C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	resetTimingId	= ARM_ConvPayResetRule (C_resetTiming);

	payTimingId		= ARM_ConvPayResetRule (C_payTiming  );

	if(C_ccy == "DEFAULT")
	{
		ccyId = ARM_NULL_OBJECT;
	}
	else
	{
		ccyId = LocalGetNumObjectId (C_ccy);
	}


	// Strike treatment
	//----------------------------------------------

	if ( strike_type == XL_TYPE_STRING )
	{
	   C_strike_double = (double) LocalGetNumObjectId(C_strike_str);

	   strike_type = 1;
	}
	else
	{
	   strike_type = 0;
	}


	// Ids recovery
	//----------------------------------------------
	
	intRuleId  = ARM_ConvIntRule (C_intRule );

	stubRuleId = ARM_ConvStubRule(C_stubRule);

	long retCode;
	long objId;
	CCString prevClass;
	
	CCString curClass = LOCAL_SPREAD_OPTION_CLASS;
	CCString stringId = GetLastCurCellEnvValue ();

	if(!stringId)
	{
		retCode = ARMLOCAL_SPREADCORRIDORVMS(C_startDate,
										     C_endDate,
											 capOrFloorId,
											 (long) strike_type,
											 C_strike_double,
										     payIndexId,
											 (long) C_CMSIndexes1_double,
											 (long) C_CMSIndexes2_double,
										     C_weight1,
										     C_weight2,
										     dayCountId,
										     resetFreqId,
										     payFreqId,
										     resetTimingId,
										     payTimingId,
										     ccyId,
										     (long) C_resetGap,
										     C_spread1,
										     C_spread2,
										     (long) C_fixing1_double,
										     (long) C_fixing2_double,
										     intRuleId,
										     stubRuleId,
										     C_slopeFlag,
										     (int) C_cptStrikeMethod,
										     payMargin_type,
										     payMarginId,
										     C_payIdxWeight,
										     (long) C_fixingPay_double,
										     (int) C_freezeFixing,
										     (int) C_computedFormula,
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
			retCode = ARMLOCAL_SPREADCORRIDORVMS(C_startDate,
											     C_endDate,
												 capOrFloorId,
												 (long) strike_type,
												 C_strike_double,
												 payIndexId,
												 (long) C_CMSIndexes1_double,
												 (long) C_CMSIndexes2_double,
												 C_weight1,
												 C_weight2,
												 dayCountId,
												 resetFreqId,
												 payFreqId,
												 resetTimingId,
												 payTimingId,
												 ccyId,
												 (long) C_resetGap,
												 C_spread1,
												 C_spread2,
												 (long) C_fixing1_double,
												 (long) C_fixing2_double,
												 intRuleId,
												 stubRuleId,
												 C_slopeFlag,
												 (int) C_cptStrikeMethod,
												 payMargin_type,
												 payMarginId,
												 C_payIdxWeight,
												 (long) C_fixingPay_double,
												 (int) C_freezeFixing,
												 (int) C_computedFormula,
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

			retCode = ARMLOCAL_SPREADCORRIDORVMS(C_startDate,
												 C_endDate,
												 capOrFloorId,
												 (long) strike_type,
												 C_strike_double,
												 payIndexId,
												 (long) C_CMSIndexes1_double,
												 (long) C_CMSIndexes2_double,
												 C_weight1,
												 C_weight2,
												 dayCountId,
												 resetFreqId,
												 payFreqId,
												 resetTimingId,
												 payTimingId,
												 ccyId,
												 (long) C_resetGap,
												 C_spread1,
												 C_spread2,
												 (long) C_fixing1_double,
												 (long) C_fixing2_double,
												 intRuleId,
												 stubRuleId,
												 C_slopeFlag,
												 (int) C_cptStrikeMethod,
												 payMargin_type,
												 payMarginId,
												 C_payIdxWeight,
												 (long) C_fixingPay_double,
												 (int) C_freezeFixing,
												 (int) C_computedFormula,
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_ARM_SPREADCORRIDORVMS" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


//------------------------------------------------------------------------------------------//
// Interface function for VMS spread options corridor -- Interface to use with VBA			//
//------------------------------------------------------------------------------------------//
__declspec(dllexport) LPXLOPER WINAPI Local_PXL_ARM_SPREADCORRIDORVMS(LPXLOPER XL_startDate,
																	  LPXLOPER XL_endDate,
 																	  LPXLOPER XL_capOrFloor,
																	  LPXLOPER XL_strike,
																	  LPXLOPER XL_spreads,				//Vect: spread1, spread2
																	  LPXLOPER XL_payIndexParams,		//Vect: payIndexObj, itsMargin,itsWeight
																	  LPXLOPER XL_CMSIndexes1,		
																	  LPXLOPER XL_CMSIndexes2,		
																	  LPXLOPER XL_sprdWeights_slope,	//Vect: Idx Weight1,Weight2, slopeflag and cptStrikeMethod
																	  LPXLOPER XL_prodParamDatas,		//Vect: dayCount, resetFreq, payFreq, resetTiming, payTiming, intRule, stubRule,resetGap 
																	  LPXLOPER XL_currency,
																	  LPXLOPER XL_fixing1,
																	  LPXLOPER XL_fixing2,
																	  LPXLOPER XL_fixingPay,
																	  LPXLOPER XL_freezeFixing)
{
	ADD_LOG("Local_PXL_ARM_SPREADCORRIDORVMS");
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

	CCString C_capOrFloor;
	long capOrFloorId;

	double C_strike_double;
	CCString C_strike_str;
	long strike_type;
	
	double C_spread1;
	double C_spread2;
	VECTOR<double> C_spreads;

	CCString C_payIndex;
	long payIndexId;

	CCString C_payMargin_str;
	long payMargin_type;
	double payMarginId;

	double C_payIdxWeight;

	CCString C_payIdxFixedRate_str;
		
	VECTOR<CCString> payIndexParams(3);
   	VECTOR<CCString> C_payIndexParams;
	VECTOR<CCString> C_payIndexParamsDF(3);
	C_payIndexParamsDF[0] = payIndexParams[0] = "DEFAULT";	// pay index objet
	C_payIndexParamsDF[1] = payIndexParams[1] = "DEFAULT";  // Margin
	C_payIndexParamsDF[2] = payIndexParams[2] = "1";		// weight

	CCString C_CMSIndexes1;
	double   C_CMSIndexes1_double;

	CCString C_CMSIndexes2;
	double   C_CMSIndexes2_double;

	double C_weight1, C_weight2, C_cptStrikeMethod;
	VECTOR<double> weights_slope(2,1.0);
	VECTOR<double> C_weights_slope;
	VECTOR<double> C_weights_slope_df(2, 1.0);

	long C_slopeFlag;
	long C_slopeFlagDefault = 1;

	double C_cptStrikeMethod_default = 1.0;

	double C_computedFormula;
	double C_computedFormula_default = 1.0;

	// prodParamDatas
	VECTOR<CCString> prodParamDatas(8);
	VECTOR<CCString> C_prodParamDatas;
	VECTOR<CCString> C_prodParamDatasDF(8);
	C_prodParamDatasDF[0] = prodParamDatas[0] = "30/360";	// dayCount
	C_prodParamDatasDF[1] = prodParamDatas[1] = "A";		// resetFreq
	C_prodParamDatasDF[2] = prodParamDatas[2] = "A";		// payFreq
	C_prodParamDatasDF[3] = prodParamDatas[3] = "ADV";		// resetTiming
	C_prodParamDatasDF[4] = prodParamDatas[4] = "ARR";		// payTiming
	C_prodParamDatasDF[5] = prodParamDatas[5] = "ADJ";		// intRule
	C_prodParamDatasDF[6] = prodParamDatas[6] = "SS";		// stubRule
	C_prodParamDatasDF[7] = prodParamDatas[7] = "10000.0";	// resetGap

	CCString C_daycount;
	long dayCountId;

	CCString C_resetFreq;
	long resetFreqId;

	CCString C_payFreq;
	long payFreqId;

	CCString C_resetTiming;
	long resetTimingId;
	
	CCString C_payTiming;
	long payTimingId;

	CCString C_intRule;
    long intRuleId;
	
	CCString C_stubRule;
    long stubRuleId;

    double C_resetGap;

	CCString C_ccy;
    long ccyId;

	CCString C_fixing1;
	double   C_fixing1_double;

	CCString C_fixing2;
	double   C_fixing2_double;

	CCString C_fixingPay;
	double   C_fixingPay_double;

	double C_freezeFixing;
	double C_freezeFixing_default = 0.0;

	// error
	static int error;
	static char* reason = "";


	// Parameters reading
	//----------------------------------------------

	XL_readNumCell(XL_startDate,C_startDate," ARM_ERR: start date: date expected",C_result);

	XL_readNumCell(XL_endDate,C_endDate," ARM_ERR: end date: date expected",C_result);
	
	XL_readStrCell(XL_capOrFloor,C_capOrFloor," ARM_ERR: cap or floor: string expected",C_result);
	
	XL_readStrOrNumCell(XL_strike,C_strike_str,C_strike_double,strike_type," ARM_ERR: strike: string or numeric expected",C_result);
	
	XL_readNumVector(XL_spreads,C_spreads," ARM_ERR: spread vector: array of numeric expected", C_result);
	
	XL_readStrVectorWD(XL_payIndexParams,C_payIndexParams,C_payIndexParamsDF," ARM_ERR: payIndex Params: array expected",DOUBLE_TYPE,C_result);	
	
	XL_readStrCellWD(XL_CMSIndexes1, C_CMSIndexes1, "", "ARM_ERR: CMS indexes : string expected", C_result);
	
	XL_readStrCellWD(XL_CMSIndexes2, C_CMSIndexes2, "", "ARM_ERR: CMS indexes : string expected", C_result);
	
	XL_readNumVectorWD(XL_sprdWeights_slope,C_weights_slope,C_weights_slope_df," ARM_ERR: weights_slopeflag: array of numeric expected", C_result);
	
	XL_readStrVectorWD(XL_prodParamDatas,C_prodParamDatas,C_prodParamDatasDF," ARM_ERR: daycount,resetFreq,payFreq,resetTiming,payTiming,resetGap,intRule,stubRule: array of string expected",DOUBLE_TYPE,C_result);
	
	XL_readStrCellWD(XL_currency,C_ccy,"DEFAULT"," ARM_ERR: currency : object expected",C_result);
	
	XL_readStrCellWD(XL_fixing1,C_fixing1,""," ARM_ERR: CMS fixing : string expected",C_result);
	
	XL_readStrCellWD(XL_fixing2,C_fixing2,""," ARM_ERR: CMS fixing : string expected",C_result);
	
	XL_readStrCellWD(XL_fixingPay,C_fixingPay,""," ARM_ERR: Pay fixing : string expected",C_result);
	
	XL_readNumCellWD(XL_freezeFixing,C_freezeFixing,C_freezeFixing_default," ARM_ERR: freeze fixing : numeric expected",C_result);

	
	// Spreads treatment
	//----------------------------------------------

	if(C_spreads.size()<2)
	{
		C_result.setMsg ("ARM_ERR: Verify the Vector of Spread");
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	C_spread1 = C_spreads[0];
	C_spread2 = C_spreads[1];


	// Pay index type treatment
	//----------------------------------------------

	int paySize = C_payIndexParams.size() < 3 ? C_payIndexParams.size() : 3;

	for(int i=0; i< paySize; i++)
	{
		if( strcmp(C_payIndexParams[i], "DEFAULT") != 0)	// fill in the value not by default
			payIndexParams[i] = C_payIndexParams[i];
	}

	C_payIndex		= payIndexParams[0];
	C_payMargin_str = payIndexParams[1];
	C_payIdxWeight	= atof(payIndexParams[2]);
	
	if(C_payIndex == "DEFAULT")
	{
		payIndexId = ARM_NULL_OBJECT;
	}
	else
	{
		payIndexId = LocalGetNumObjectId(C_payIndex);
	}

	if(C_payMargin_str == "DEFAULT")
	{
		payMarginId = 0.0;
		payMargin_type = 0;				// number
	}
	else
	{
		payMarginId = LocalGetNumObjectId(C_payMargin_str);

		if( (int)payMarginId == -1 )	// no object
		{
			payMarginId = atof(C_payMargin_str);
			payMargin_type = 0;
		}
		else
			payMargin_type = 1;			// object
	}


	// CMS Indexes treatment
	//----------------------------------------------

	if ( (XL_CMSIndexes1->xltype == xltypeMissing) || (XL_CMSIndexes1->xltype == xltypeNil) )
	{		
		C_CMSIndexes1_double = ARM_NULL_OBJECT;
    }
	else
	{
		C_CMSIndexes1_double = (double) LocalGetNumObjectId(C_CMSIndexes1);
	}



	if ( (XL_CMSIndexes2->xltype == xltypeMissing) || (XL_CMSIndexes2->xltype == xltypeNil) )
	{		
		C_CMSIndexes2_double = ARM_NULL_OBJECT;
    }
	else
	{
		C_CMSIndexes2_double = (double) LocalGetNumObjectId(C_CMSIndexes2);
	}


	// Spread weights treatment
	//----------------------------------------------

	if(C_weights_slope.size()<2)
	{
		C_result.setMsg ("ARM_ERR: Verify the Vector of Weight");
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}
	
	C_weight1 = C_weights_slope[0];
	C_weight2 = C_weights_slope[1];

	if(C_weights_slope.size() == 2)
	{
		C_slopeFlag			= C_slopeFlagDefault;
		C_cptStrikeMethod	= C_cptStrikeMethod_default;
		C_computedFormula	= C_computedFormula_default;
	}
	else if(C_weights_slope.size() == 3)
	{
		C_slopeFlag			= C_weights_slope[2];
		C_cptStrikeMethod	= C_cptStrikeMethod_default;
		C_computedFormula	= C_computedFormula_default;
	}
	else if(C_weights_slope.size() == 4)
	{
		C_slopeFlag			= C_weights_slope[2];
		C_cptStrikeMethod	= C_weights_slope[3];
		C_computedFormula	= C_computedFormula_default;
	}
	else
	{
		C_slopeFlag			= C_weights_slope[2];
		C_cptStrikeMethod	= C_weights_slope[3];
		C_computedFormula	= C_weights_slope[4];
	}


	// prodParamDatas treatment
	//----------------------------------------------

	int prodParamSize = C_prodParamDatas.size() < 8 ? C_prodParamDatas.size() : 8;
	for (i =0; i<prodParamSize; i++)
	{
		if( strcmp(C_prodParamDatas[i],"DEFAULT")!= 0 )
			prodParamDatas[i] = C_prodParamDatas[i];
	}
	// dayCount
	if (prodParamDatas.size() >= 1)
		C_daycount    = prodParamDatas[0];
	else
		C_daycount    = C_prodParamDatasDF[0];
	// resetFreq
	if (prodParamDatas.size() >= 2)
		C_resetFreq   = prodParamDatas[1];
	else
		C_resetFreq   = C_prodParamDatasDF[1];
	// payFreq
	if (prodParamDatas.size() >= 3)
		C_payFreq     = prodParamDatas[2];
	else
		C_payFreq     = C_prodParamDatasDF[2];
	// resetTiming
	if (prodParamDatas.size() >= 4)
		C_resetTiming = prodParamDatas[3];
	else
		C_resetTiming = C_prodParamDatasDF[3];
	// payTiming
	if (prodParamDatas.size() >= 5)
		C_payTiming   = prodParamDatas[4];
	else
		C_payTiming   = C_prodParamDatasDF[4];
	// intRule
	if (prodParamDatas.size() >= 6)
		C_intRule     = prodParamDatas[5];
	else
		C_intRule     = C_prodParamDatasDF[5];
	// stubRule
	if (prodParamDatas.size() >= 7)
		C_stubRule    = prodParamDatas[6];
	else
		C_stubRule    = C_prodParamDatasDF[6];
	// resetGap
	if (prodParamDatas.size() >= 8)
		C_resetGap    = atof(prodParamDatas[7]);
	else
		C_resetGap    = atof(C_prodParamDatasDF[7]);


	// Fixing1 treatment
	//----------------------------------------------

	if ((XL_fixing1->xltype == xltypeMissing) || (XL_fixing1->xltype == xltypeNil))
	{
		C_fixing1_double = ARM_NULL_OBJECT;
    }
	else
	{
		C_fixing1_double = (double) LocalGetNumObjectId(C_fixing1);
	}

	
	// Fixing2 treatment
	//----------------------------------------------
	
	if ((XL_fixing2->xltype == xltypeMissing) || (XL_fixing2->xltype == xltypeNil))
	{
		C_fixing2_double = ARM_NULL_OBJECT;
    }
	else
	{
		C_fixing2_double = (double) LocalGetNumObjectId(C_fixing2);
	}


	// Fixing3 treatment
	//----------------------------------------------
	
	if ((XL_fixingPay->xltype == xltypeMissing) || (XL_fixingPay->xltype == xltypeNil))
	{
		C_fixingPay_double = ARM_NULL_OBJECT;
    }
	else
	{
		C_fixingPay_double = (double) LocalGetNumObjectId(C_fixingPay);

	}


	// Ids recovery
	//----------------------------------------------

	if((capOrFloorId = ARM_ConvCapOrFloor (C_capOrFloor, C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	dayCountId		= ARM_ConvDayCount (C_daycount);
	
	if((resetFreqId = ARM_ConvFrequency (C_resetFreq, C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	if((payFreqId = ARM_ConvFrequency (C_payFreq, C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	resetTimingId	= ARM_ConvPayResetRule (C_resetTiming);

	payTimingId		= ARM_ConvPayResetRule (C_payTiming  );

	if(C_ccy == "DEFAULT")
	{
		ccyId = ARM_NULL_OBJECT;
	}
	else
	{
		ccyId = LocalGetNumObjectId (C_ccy);
	}


	// Strike treatment
	//----------------------------------------------

	if ( strike_type == XL_TYPE_STRING )
	{
	   C_strike_double = (double) LocalGetNumObjectId(C_strike_str);

	   strike_type = 1;
	}
	else
	{
	   strike_type = 0;
	}


	// Ids recovery
	//----------------------------------------------
	
	intRuleId  = ARM_ConvIntRule (C_intRule );

	stubRuleId = ARM_ConvStubRule(C_stubRule);

	long retCode;
	long objId;
	
	CCString curClass = LOCAL_SPREAD_OPTION_CLASS;
	CCString stringId;

	retCode = ARMLOCAL_SPREADCORRIDORVMS(C_startDate,
										 C_endDate,
										 capOrFloorId,
										 (long) strike_type,
										 C_strike_double,
										 payIndexId,
										 (long) C_CMSIndexes1_double,
										 (long) C_CMSIndexes2_double,
										 C_weight1,
										 C_weight2,
										 dayCountId,
										 resetFreqId,
										 payFreqId,
										 resetTimingId,
										 payTimingId,
										 ccyId,
										 (long) C_resetGap,
										 C_spread1,
										 C_spread2,
										 (long) C_fixing1_double,
										 (long) C_fixing2_double,
										 intRuleId,
										 stubRuleId,
										 C_slopeFlag,
										 (int) C_cptStrikeMethod,
										 payMargin_type,
										 payMarginId,
										 C_payIdxWeight,
										 (long) C_fixingPay_double,
										 (int) C_freezeFixing,
										 (int) C_computedFormula,
										 C_result);


	if(retCode == ARM_OK)
	{
		objId	 = C_result.getLong ();

		stringId = LocalMakeObjectId (objId, curClass);
	}

	if(retCode == ARM_OK)
	{
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_ARM_SPREADCORRIDORVMS" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


__declspec(dllexport) LPXLOPER WINAPI Local_ARM_QUANTOSPREADCORRIDOR(	LPXLOPER XL_startDate,
																		LPXLOPER XL_endDate,
 																		LPXLOPER XL_capOrFloor,
																		LPXLOPER XL_strikes,
																		LPXLOPER XL_liborIdx1,
																		LPXLOPER XL_liborIdx2,
																		LPXLOPER XL_liborIdxP,
																		LPXLOPER XL_fixedRateP,
																		LPXLOPER XL_leg1Weight,
																		LPXLOPER XL_leg2Weight,
																		LPXLOPER XL_legPWeight,
																		LPXLOPER XL_leg1Fixings,
																		LPXLOPER XL_leg2Fixings,
																		LPXLOPER XL_legPFixings,
																		LPXLOPER XL_leg1Spread,
																		LPXLOPER XL_leg2Spread,
																		LPXLOPER XL_modelParams,
																		LPXLOPER XL_prodParams,
																		LPXLOPER XL_notional)
{
	ADD_LOG("Local_ARM_QUANTOSPREADCORRIDOR");
	//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();

	double C_startDate;
	double C_endDate;

	CCString C_capOrFloor;
	long capOrFloorId;

	CCString C_strikes_str;
	double C_strike_double;
	long strike_type;

	CCString C_liborIdx1;
	long liborIdx1Id;
	CCString C_liborIdx2;
	long liborIdx2Id;
	CCString C_liborIdxP;
	long liborIdxPId;

	CCString C_fixedRateP_str;
	CCString C_fixedRateP_str_default("DEFAULT");
	long fixedRatePId = ARM_NULL_OBJECT;

	double C_Idx1weight;
	double C_Idx1weight_default = 1.0;

	double C_Idx2weight;
	double C_Idx2weight_default = 1.0;

	double C_IdxPweight;
	double C_IdxPweight_default = 1.0;
	
	CCString C_Idx1fixings_str;
	CCString C_Idx1fixings_str_default("DEFAULT");
	long Idx1fixingsId = ARM_NULL_OBJECT;
	CCString C_Idx2fixings_str;
	CCString C_Idx2fixings_str_default("DEFAULT");
	long Idx2fixingsId = ARM_NULL_OBJECT;
	CCString C_IdxPfixings_str;
	CCString C_IdxPfixings_str_default("DEFAULT");
	long IdxPfixingsId = ARM_NULL_OBJECT;

	double C_Idx1spread;
	double C_Idx1spread_default = 0.0;

	double C_Idx2spread;
	double C_Idx2spread_default = 0.0;

	VECTOR<double> C_modelParams;
	VECTOR<double> C_modelParams_default(0);
	VECTOR<double> modelParams(3);
	
	modelParams[0] = 1.0;
	modelParams[1] = 1.0;
	modelParams[2] = 1.0;

	double C_slopeFlag;
	double C_cptStrikeMethod;
	double C_computedFormula;

	// Product parameters
	VECTOR<CCString> C_prodParams;
	VECTOR<CCString> C_prodParams_default(0);

	VECTOR<CCString> prodParams(15);

	prodParams[0] = "DEFAULT";
	prodParams[1] = "30/360";
	prodParams[2] = "A";
	prodParams[3] = "A";
	prodParams[4] = "ADV";
	prodParams[5] = "ARR";
	prodParams[6] = "ADJ";
	prodParams[7] = "SS";
	prodParams[8] = "-2";
	prodParams[9] = "EUR";
	prodParams[10] = "EUR";
	prodParams[11] = "EUR";
	prodParams[12] = "MF";
	prodParams[13] = "NULL";
	prodParams[14] = "0";

	CCString C_ccy;
	long ccyId;
	
	CCString C_daycount;
	long dayCountId;

	CCString C_resetFreq;
	long resetFreqId;

	CCString C_payFreq;
	long payFreqId;

	CCString C_resetTiming;
	long resetTimingId;
	
	CCString C_payTiming;
	long payTimingId;

	CCString C_intRule;
    long intRuleId;
	
	CCString C_stubRule;
    long stubRuleId;

    double C_resetGap;
	
	CCString C_resetCal;
	CCString C_payCal;
	CCString C_payIndexResetCal;

	CCString C_fwdRule;
	double fwdRuleId;

	CCString C_refDate;

	double C_freezeFixing;
	
	long notionalId = ARM_NULL_OBJECT;
	CCString C_notional_str;
	CCString C_notional_default("DEFAULT");

	// error
	static int error;
	static char* reason = "";

	XL_readNumCell(XL_startDate,C_startDate," ARM_ERR: start date: date expected",C_result);
	XL_readNumCell(XL_endDate,C_endDate," ARM_ERR: end date: date expected",C_result);
	XL_readStrCell(XL_capOrFloor,C_capOrFloor," ARM_ERR: cap or floor: string expected",C_result);
	XL_readStrOrNumCell(XL_strikes,C_strikes_str,C_strike_double,strike_type," ARM_ERR: strikes vector : array of numeric expected",C_result);
	XL_readStrCellWD(XL_liborIdx1,C_liborIdx1,"DEFAULT"," ARM_ERR: leg1 libor index : string expected",C_result);
	XL_readStrCellWD(XL_liborIdx2,C_liborIdx2,"DEFAULT"," ARM_ERR: leg2 libor index : string expected",C_result);
	XL_readStrCellWD(XL_liborIdxP,C_liborIdxP,"DEFAULT"," ARM_ERR: paiement leg libor index : string expected",C_result);
	XL_readStrCellWD(XL_fixedRateP,C_fixedRateP_str,C_fixedRateP_str_default," ARM_ERR: fixed rate vector for paiement leg : object expected",C_result);
	XL_readNumCellWD(XL_leg1Weight,C_Idx1weight,C_Idx1weight_default," ARM_ERR: leg1 weight : numeric expected",C_result);
	XL_readNumCellWD(XL_leg2Weight,C_Idx2weight,C_Idx2weight_default," ARM_ERR: leg2 weight : numeric expected",C_result);
	XL_readNumCellWD(XL_legPWeight,C_IdxPweight,C_IdxPweight_default," ARM_ERR: legP weight : numeric expected",C_result);
	XL_readStrCellWD(XL_leg1Fixings,C_Idx1fixings_str,C_Idx1fixings_str_default," ARM_ERR: leg1 fixings vector : object expected",C_result);
	XL_readStrCellWD(XL_leg2Fixings,C_Idx2fixings_str,C_Idx2fixings_str_default," ARM_ERR: leg2 fixings vector : object expected",C_result);
	XL_readStrCellWD(XL_legPFixings,C_IdxPfixings_str,C_IdxPfixings_str_default," ARM_ERR: legP fixings vector : object expected",C_result);
	XL_readNumCellWD(XL_leg1Spread,C_Idx1spread, C_Idx1spread_default, " ARM_ERR: leg1 spread : numeric expected",C_result);
	XL_readNumCellWD(XL_leg2Spread,C_Idx2spread,C_Idx2spread_default, " ARM_ERR: leg2 spread : numeric expected",C_result);
	XL_readNumVectorWD(XL_modelParams,C_modelParams,C_modelParams_default," ARM_ERR: model parameters vector : array of numeric expected",C_result);
	XL_readStrVectorWD(XL_prodParams,C_prodParams,C_prodParams_default," ARM_ERR: product parameters vector : array of numeric expected",DOUBLE_TYPE,C_result);
	XL_readStrCellWD(XL_notional,C_notional_str,C_notional_default," ARM_ERR: notional vector : object expected",C_result);

	if((capOrFloorId = ARM_ConvCapOrFloor (C_capOrFloor, C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}
		
	if ( strike_type == XL_TYPE_STRING )
	{
	   C_strike_double = (double) LocalGetNumObjectId(C_strikes_str);

	   strike_type = 1;
	}
	else
	{
	   strike_type = 0;
	}

	liborIdx1Id = strcmp(C_liborIdx1,"DEFAULT") ? (double) LocalGetNumObjectId(C_liborIdx1) : -1.0;	
	liborIdx2Id = strcmp(C_liborIdx2,"DEFAULT") ? (double) LocalGetNumObjectId(C_liborIdx2) : -1.0;	
	liborIdxPId = strcmp(C_liborIdxP,"DEFAULT") ? (double) LocalGetNumObjectId(C_liborIdxP) : -1.0;

	fixedRatePId = strcmp(C_Idx1fixings_str,"DEFAULT") ? (long) LocalGetNumObjectId(C_fixedRateP_str) : -1.0;

	Idx1fixingsId = strcmp(C_Idx1fixings_str,"DEFAULT") ? (long) LocalGetNumObjectId(C_Idx1fixings_str) : -1.0;
	Idx2fixingsId = strcmp(C_Idx2fixings_str,"DEFAULT") ? (long) LocalGetNumObjectId(C_Idx2fixings_str) : -1.0;
	IdxPfixingsId = strcmp(C_IdxPfixings_str,"DEFAULT") ? (long) LocalGetNumObjectId(C_IdxPfixings_str) : -1.0;
	
	int i = 0;
	for (i; i < C_modelParams.size(); i++)
	{
		modelParams[i]  = C_modelParams[i];
	}
	
	C_slopeFlag = modelParams[0];
	C_cptStrikeMethod = modelParams[1];
	C_computedFormula = modelParams[2];

	for (i = 0; i<C_prodParams.size(); i++)
	{
		if (strcmp(C_prodParams[i],"DEFAULT"))
			prodParams[i] = C_prodParams[i];
	}

	C_ccy = prodParams[0];
	
	if(C_ccy == "DEFAULT")
	{
		ccyId = ARM_NULL_OBJECT;
	}
	else
	{
		ccyId = LocalGetNumObjectId (C_ccy);
	}
	
	C_daycount = prodParams[1];
	dayCountId = ARM_ConvDayCount (C_daycount);
	
	C_resetFreq = prodParams[2];
	if((resetFreqId = ARM_ConvFrequency (C_resetFreq, C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}
	
	C_payFreq = prodParams[3];
	if((payFreqId = ARM_ConvFrequency (C_payFreq, C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	C_resetTiming = prodParams[4];
	resetTimingId = ARM_ConvPayResetRule (C_resetTiming);
	
	C_payTiming = prodParams[5];
	payTimingId = ARM_ConvPayResetRule (C_payTiming);

	C_intRule  = prodParams[6];
    intRuleId = ARM_ConvIntRule(C_intRule);
	
	C_stubRule = prodParams[7];
	stubRuleId = ARM_ConvStubRule(C_stubRule);

    C_resetGap = atof(prodParams[8]);
	
	C_resetCal = prodParams[9];

	C_payCal = prodParams[10];

	C_payIndexResetCal = prodParams[11];

	C_fwdRule = prodParams[12];
	fwdRuleId = ARM_ConvFwdRule(C_fwdRule);

	C_refDate = prodParams[13];

	C_freezeFixing = atof(prodParams[14]);

	notionalId = strcmp(C_notional_str,"DEFAULT") ? (long) LocalGetNumObjectId(C_notional_str) : -1.0;

	long retCode;
	long objId;
	
	CCString prevClass;
	CCString curClass = LOCAL_SPREAD_OPTION_CLASS;
	CCString stringId = GetLastCurCellEnvValue();
	
	if(!stringId)
	{
		retCode = ARMLOCAL_QUANTOSPREADCORRIDOR (	C_startDate,
													C_endDate,
													capOrFloorId,
													(long) strike_type,
													C_strike_double,
													(long) liborIdx1Id,
													(long) liborIdx2Id,
													(long) liborIdxPId,
													fixedRatePId,
													C_Idx1weight,
													C_Idx2weight,
													C_IdxPweight,
													Idx1fixingsId,
													Idx2fixingsId,
													IdxPfixingsId,
													C_Idx1spread,
													C_Idx2spread,
													C_slopeFlag,
													C_cptStrikeMethod,
													(int) C_computedFormula,
													ccyId,
													dayCountId,
													resetFreqId,
													payFreqId,
													resetTimingId,
													payTimingId,
													intRuleId,
													stubRuleId,
													C_resetGap,
													C_resetCal,
													C_payIndexResetCal,
													C_payCal,
													fwdRuleId,
													C_refDate,
													notionalId,
													C_freezeFixing,
													C_result);
		if (retCode == ARM_OK)
		{
			objId = C_result.getLong();

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
			retCode = ARMLOCAL_QUANTOSPREADCORRIDOR(C_startDate,
													C_endDate,
													capOrFloorId,
													(long) strike_type,
													C_strike_double,
													(long) liborIdx1Id,
													(long) liborIdx2Id,
													(long) liborIdxPId,
													fixedRatePId,
													C_Idx1weight,
													C_Idx2weight,
													C_IdxPweight,
													Idx1fixingsId,
													Idx2fixingsId,
													IdxPfixingsId,
													C_Idx1spread,
													C_Idx2spread,
													C_slopeFlag,
													C_cptStrikeMethod,
													(int) C_computedFormula,
													ccyId,
													dayCountId,
													resetFreqId,
													payFreqId,
													resetTimingId,
													payTimingId,
													intRuleId,
													stubRuleId,
													C_resetGap,
													C_resetCal,
													C_payCal,
													C_payIndexResetCal,
													fwdRuleId,
													C_refDate,
													notionalId,
													C_freezeFixing,
													C_result,
													objId);
			if (retCode == ARM_OK)
			{
				objId = C_result.getLong();

				stringId = LocalMakeObjectId (objId, curClass);
			}
		}
		else
		{
			FreeCurCellContent ();
			retCode = ARMLOCAL_QUANTOSPREADCORRIDOR (	C_startDate,
													C_endDate,
													capOrFloorId,
													(long) strike_type,
													C_strike_double,
													(long) liborIdx1Id,
													(long) liborIdx2Id,
													(long) liborIdxPId,
													fixedRatePId,
													C_Idx1weight,
													C_Idx2weight,
													C_IdxPweight,
													Idx1fixingsId,
													Idx2fixingsId,
													IdxPfixingsId,
													C_Idx1spread,
													C_Idx2spread,
													C_slopeFlag,
													C_cptStrikeMethod,
													(int) C_computedFormula,
													ccyId,
													dayCountId,
													resetFreqId,
													payFreqId,
													resetTimingId,
													payTimingId,
													intRuleId,
													stubRuleId,
													C_resetGap,
													C_resetCal,
													C_payCal,
													C_payIndexResetCal,
													fwdRuleId,
													C_refDate,
													notionalId,
													C_freezeFixing,
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_ARM_QUANTOSPREADCORRIDOR"  )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


__declspec(dllexport) LPXLOPER WINAPI Local_ARM_CORRIDORDBLCONDITION(LPXLOPER XL_startDate,
																	 LPXLOPER XL_endDate,
																	 LPXLOPER XL_DigitalcapOrFloor,
																	 LPXLOPER XL_SpreadcapOrFloor,
																	 LPXLOPER XL_digitalBarrier, 
																	 LPXLOPER XL_spreadBarrier,
																	 LPXLOPER XL_spreads,  //vect: spread1, spread2
																	 LPXLOPER XL_payIndexParams,//vect:payIndexObj,  itsMargin, itsWeight, itsFixedRate
																	 LPXLOPER XL_spreadIdxTypes,//vect: IndexType1, IndexType2, IndexType3
																	 LPXLOPER XL_Weights,//vect: Weight1,Weight2, Weight3
																	 LPXLOPER XL_prodParamDatas,//vect: dayCount, resetFreq, payFreq, resetTiming, payTiming, intRule, stubRule,resetGap 
																	 LPXLOPER XL_currency,
																	 LPXLOPER XL_fixing1,
																	 LPXLOPER XL_fixing2,
																	 LPXLOPER XL_fixing3,
																	 LPXLOPER XL_fixingPay,
																	 LPXLOPER XL_freezeFixing,
																	 LPXLOPER XL_DigitalSpreadCorrel,
																	 LPXLOPER XL_ShiftVolCorrel,
																	 LPXLOPER XL_fixing4)
{
	ADD_LOG("Local_ARM_CORRIDORDBLCONDITION");
	static XLOPER XL_result;
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();

	double C_startDate;
	double C_endDate;

	CCString C_DigitalcapOrFloor;
	long DigitalcapOrFloorId;
	CCString C_SpreadcapOrFloor;
	long SpreadcapOrFloorId;

	double C_digitalBarrier_double;
	CCString C_digitalBarrier_str;
	long digitalBarrier_type;
	double C_spreadBarrier_double;
	CCString C_spreadBarrier_str;
	long spreadBarrier_type;

	double C_ShiftVolCorrel_double;
	double C_Default_double = 0.0;
	CCString C_ShiftVolCorrel_str;
	long shiftVolCorrel_type;

	VECTOR<double> C_spreads;

	CCString C_payIndex;
	long payIndexId;

	CCString C_resetCal;
	CCString C_payCal;
	char* resetCal = NULL;
	char* payCal = NULL;
	
	CCString C_payMargin_str;
	long payMargin_type;
	double payMarginId;

	double C_payIdxWeight;

	VECTOR<CCString> payIndexParams(3);
   	VECTOR<CCString> C_payIndexParams;
	VECTOR<CCString> C_payIndexParamsDF(3);
	C_payIndexParamsDF[0] = payIndexParams[0] = "DEFAULT";// pay index
	C_payIndexParamsDF[1] = payIndexParams[1] = "DEFAULT";// pay margin referenceValue
	C_payIndexParamsDF[2] = payIndexParams[2] = "1";// pay weight
	
	CCString C_liborType1;
	long liborType1Id;
	CCString C_liborType2;
	long liborType2Id;
	CCString C_liborType3;
	long liborType3Id;
	CCString C_liborType4;
	long liborType4Id;
	VECTOR<CCString> liborTypes(4);
   	VECTOR<CCString> C_liborTypes;
	VECTOR<CCString> C_liborTypesDF(4);
	C_liborTypesDF[0] = liborTypes[0] = "LIBOR3M";
	C_liborTypesDF[1] = liborTypes[1] = "LIBOR3M";
	C_liborTypesDF[2] = liborTypes[2] = "LIBOR3M";
	C_liborTypesDF[3] = liborTypes[3] = "LIBOR3M";

	double C_weight1, C_weight2, C_weight3, C_weight4;
	VECTOR<double> weights(4,1.0);
	VECTOR<double> C_weights;
	VECTOR<double> C_weights_df(4, 1.0);

	VECTOR<CCString> prodParamDatas(11);
	VECTOR<CCString> C_prodParamDatas;
	VECTOR<CCString> C_prodParamDatasDF(11);
	C_prodParamDatasDF[0] = prodParamDatas[0] = "30/360"; // dayCount
	C_prodParamDatasDF[1] = prodParamDatas[1] = "A"; // resetFreq
	C_prodParamDatasDF[2] = prodParamDatas[2] = "A"; // payFreq
	C_prodParamDatasDF[3] = prodParamDatas[3] = "ADV"; // resetTiming
	C_prodParamDatasDF[4] = prodParamDatas[4] = "ARR"; // payTiming
	C_prodParamDatasDF[5] = prodParamDatas[5] = "ADJ"; // intRule
	C_prodParamDatasDF[6] = prodParamDatas[6] = "SS";  // stubRule
	C_prodParamDatasDF[7] = prodParamDatas[7] = "10000.0";// resetGap
	C_prodParamDatasDF[8] = prodParamDatas[8] = "N";	// no correlation corrector
	C_prodParamDatasDF[9] = prodParamDatas[9] = "";	// resetCal
	C_prodParamDatasDF[10] = prodParamDatas[10] = "";	// payCal

	CCString C_daycount;
	long dayCountId;

	CCString C_resetFreq;
	long resetFreqId;

	CCString C_payFreq;
	long payFreqId;

	CCString C_resetTiming;
	long resetTimingId;
	
	CCString C_payTiming;
	long payTimingId;

	CCString C_intRule;
    long intRuleId;
	
	CCString C_stubRule;
    long stubRuleId;

	CCString C_DigitalSpreadCorrel;
	long digitalSpreadcorrelId;

    double C_resetGap;
	bool C_correlCorrectFlag;

	CCString C_ccy;
    long ccyId;

   	double C_fixing1_double;
   	VECTOR<double> C_fixing1_double_vect;
	CCString C_fixing1_str;
	long fixing1_type;
	VECTOR<double> C_fixing1_default(0.0);

   	double C_fixing2_double;
   	VECTOR<double> C_fixing2_double_vect;
	CCString C_fixing2_str;
	long fixing2_type;
	VECTOR<double> C_fixing2_default(0.0);

   	double C_fixing3_double;
   	VECTOR<double> C_fixing3_double_vect;
	CCString C_fixing3_str;
	long fixing3_type;
	VECTOR<double> C_fixing3_default(0.0);

   	double C_fixing4_double;
   	VECTOR<double> C_fixing4_double_vect;
	CCString C_fixing4_str;
	long fixing4_type;
	VECTOR<double> C_fixing4_default(0.0);

	double C_fixingPay_double;
   	VECTOR<double> C_fixingPay_double_vect;
	CCString C_fixingPay_str;
	long fixingPay_type;
	VECTOR<double> C_fixingPay_default(0.0);

	double C_freezeFixing;
	double C_freezeFixing_default = 0.0;

	// error
	static int error;
	static char* reason = "";



	XL_readNumCell(XL_startDate,C_startDate," ARM_ERR: start date: date expected",C_result);
	XL_readNumCell(XL_endDate,C_endDate," ARM_ERR: end date: date expected",C_result);
	XL_readStrCell(XL_DigitalcapOrFloor,C_DigitalcapOrFloor," ARM_ERR: cap or floor: string expected",C_result);
	XL_readStrCell(XL_SpreadcapOrFloor,C_SpreadcapOrFloor," ARM_ERR: cap or floor: string expected",C_result);
	XL_readStrOrNumCell(XL_digitalBarrier,C_digitalBarrier_str,C_digitalBarrier_double,digitalBarrier_type," ARM_ERR: Rate/Spread#1 Barrier: string or numeric expected",C_result);
	XL_readStrOrNumCell(XL_spreadBarrier,C_spreadBarrier_str,C_spreadBarrier_double,spreadBarrier_type," ARM_ERR: Spread#2 Barrier: string or numeric expected",C_result);
	XL_readNumVector(XL_spreads,C_spreads," ARM_ERR: spread vector: array of numeric expected", C_result);
	XL_readStrVectorWD(XL_payIndexParams,C_payIndexParams,C_payIndexParamsDF," ARM_ERR: payIndex Params: array expected",DOUBLE_TYPE,C_result);
	XL_readStrVectorWD(XL_spreadIdxTypes,C_liborTypes,C_liborTypesDF," ARM_ERR: spread indexes : array of string expected",DOUBLE_TYPE,C_result);
	XL_readNumVectorWD(XL_Weights,C_weights,C_weights_df," ARM_ERR: weights: array of numeric expected", C_result);
	XL_readStrVectorWD(XL_prodParamDatas,C_prodParamDatas,C_prodParamDatasDF," ARM_ERR: daycount,resetFreq,payFreq,resetTiming,payTiming,intRule,stubRule,resetGap,correlCorrector: array of string expected",DOUBLE_TYPE,C_result);
	XL_readStrCellWD(XL_currency,C_ccy,"DEFAULT"," ARM_ERR: currency : object expected",C_result);
	XL_readStrOrNumVectorWD(XL_fixing1,C_fixing1_str,C_fixing1_double_vect,C_fixing1_default,fixing1_type," ARM_ERR: fixing1 taux : array or string expected",C_result);
	XL_readStrOrNumVectorWD(XL_fixing2,C_fixing2_str,C_fixing2_double_vect,C_fixing2_default,fixing2_type," ARM_ERR: fixing2 taux : array or string expected",C_result);
	XL_readStrOrNumVectorWD(XL_fixing3,C_fixing3_str,C_fixing3_double_vect,C_fixing3_default,fixing3_type," ARM_ERR: fixing3 taux : array or string expected",C_result);
	XL_readStrOrNumVectorWD(XL_fixing4,C_fixing4_str,C_fixing4_double_vect,C_fixing4_default,fixing4_type," ARM_ERR: fixing4 taux : array or string expected",C_result);
	XL_readStrOrNumVectorWD(XL_fixingPay,C_fixingPay_str,C_fixingPay_double_vect,C_fixingPay_default,fixingPay_type," ARM_ERR: fixingPay taux : array or string expected",C_result);
	XL_readNumCellWD(XL_freezeFixing,C_freezeFixing,C_freezeFixing_default," ARM_ERR: freeze fixing : numeric expected",C_result);
	XL_readStrCellWD(XL_DigitalSpreadCorrel,C_DigitalSpreadCorrel,"DEFAULT"," ARM_ERR: VolCorrel : object expected",C_result);
	XL_readStrOrNumCellWD(XL_ShiftVolCorrel,C_ShiftVolCorrel_str, C_ShiftVolCorrel_double, C_Default_double, shiftVolCorrel_type,
			   " ARM_ERR: ShiftCorrel : object expected",C_result);

	if(C_spreads.size()<2)
	{
		C_result.setMsg ("ARM_ERR: Verify the Vector of Spread");
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}
	
	// traitement de payIndexType_itsWeight
	int paySize = C_payIndexParams.size() < 3 ? C_payIndexParams.size() : 3;
	for(int i=0; i<paySize; i++)
	{
		if( strcmp(C_payIndexParams[i], "DEFAULT") != 0)// fill in the value not by default
			payIndexParams[i] = C_payIndexParams[i];
	}

	C_payIndex = payIndexParams[0];
	C_payMargin_str = payIndexParams[1];
	C_payIdxWeight = atof(payIndexParams[2]);

	if(C_payIndex == "DEFAULT")
	{
		payIndexId = ARM_NULL_OBJECT;
	}
	else
	{
		payIndexId = LocalGetNumObjectId(C_payIndex);
	}

	if(C_payMargin_str == "DEFAULT")
	{
		payMarginId = 0.0;
		payMargin_type = 0;//nombre
	}
	else
	{
		payMarginId = LocalGetNumObjectId(C_payMargin_str);

		if( (int)payMarginId == -1 ) // non objet
		{
			payMarginId = atof(C_payMargin_str);
			payMargin_type = 0;
		}
		else
			payMargin_type = 1;//objet
	}

	//traitement de indexes et des weights
	int IdxSize = C_liborTypes.size();
	bool isRate=true;
	if(IdxSize>3)
	{
		isRate=false; // 2 spreads
		IdxSize=4;
	}
	for (i =0; i< IdxSize; i++)
	{
		if(strcmp(C_liborTypes[i],"") != 0)
			liborTypes[i] = C_liborTypes[i];
	}
	if(isRate)
	{
		C_liborType1	= liborTypes[3]; // 2nd rate of 1st spread = default value
		C_liborType2	= liborTypes[0]; // Rate or 2st rate of 1st spread
		C_liborType3	= liborTypes[1]; // 1st rate of 2nd spread
		C_liborType4	= liborTypes[2]; // 2nd rate of 2nd spread
		C_weight1		= 0.0;			 // no use of 2nd rate of 1st spread
		C_weight2		= C_weights[0];
		C_weight3		= C_weights[1];
		C_weight4		= C_weights[2];
	}
	else
	{
		C_liborType1	= liborTypes[0]; // 1st rate of 1st spread
		C_liborType2	= liborTypes[1]; // 2st rate of 1st spread
		C_liborType3	= liborTypes[2]; // 1st rate of 2nd spread
		C_liborType4	= liborTypes[3]; // 2nd rate of 2nd spread
		C_weight1		= C_weights[0];
		C_weight2		= C_weights[1];
		C_weight3		= C_weights[2];
		C_weight4		= C_weights[3];
	}

	// traitement de prodParamDatas
	int prodSize = C_prodParamDatas.size() < 11 ? C_prodParamDatas.size() : 11;
	for (i =0; i<prodSize; i++)
	{
		if( strcmp(C_prodParamDatas[i],"DEFAULT")!= 0 )
			prodParamDatas[i] = C_prodParamDatas[i];
	}
	// dayCount
	if (prodParamDatas.size() >= 1)
		C_daycount    = prodParamDatas[0];
	else
		C_daycount    = C_prodParamDatasDF[0];
	// resetFreq
	if (prodParamDatas.size() >= 2)
		C_resetFreq   = prodParamDatas[1];
	else
		C_resetFreq   = C_prodParamDatasDF[1];
	// payFreq
	if (prodParamDatas.size() >= 3)
		C_payFreq     = prodParamDatas[2];
	else
		C_payFreq     = C_prodParamDatasDF[2];
	// resetTiming
	if (prodParamDatas.size() >= 4)
		C_resetTiming = prodParamDatas[3];
	else
		C_resetTiming = C_prodParamDatasDF[3];
	// payTiming
	if (prodParamDatas.size() >= 5)
		C_payTiming   = prodParamDatas[4];
	else
		C_payTiming   = C_prodParamDatasDF[4];
	// intRule
	if (prodParamDatas.size() >= 6)
		C_intRule     = prodParamDatas[5];
	else
		C_intRule     = C_prodParamDatasDF[5];
	// stubRule
	if (prodParamDatas.size() >= 7)
		C_stubRule    = prodParamDatas[6];
	else
		C_stubRule    = C_prodParamDatasDF[6];
	// resetGap
	if (prodParamDatas.size() >= 8)
		C_resetGap    = atof(prodParamDatas[7]);
	else
		C_resetGap    = atof(C_prodParamDatasDF[7]);
	// correl corrector
	if (prodParamDatas.size() >= 9)
		C_correlCorrectFlag    = (prodParamDatas[8]=="Y" || prodParamDatas[8]=="y");
	else
		C_correlCorrectFlag    = (C_prodParamDatasDF[8]=="Y" || C_prodParamDatasDF[8]=="y");
	// resetCal
	if (prodParamDatas.size() >= 10)
		C_resetCal = prodParamDatas[9];
	else
		C_resetCal   = C_prodParamDatasDF[9];
	// payCal
	if (prodParamDatas.size() >= 11)
		C_payCal   = prodParamDatas[10];
	else
		C_payCal   = C_prodParamDatasDF[10];
	
	if(!(C_resetCal == ""))
		resetCal = C_resetCal.c_str();
	if(!(C_payCal == ""))
		payCal = C_payCal.c_str();

	
	if(fixing1_type == XL_TYPE_STRING)
		C_fixing1_double = (double) LocalGetNumObjectId(C_fixing1_str);
	else
	    C_fixing1_double = ARM_NULL_OBJECT;

	if(fixing2_type == XL_TYPE_STRING)
		C_fixing2_double = (double) LocalGetNumObjectId(C_fixing2_str);
	else
	    C_fixing2_double = ARM_NULL_OBJECT;

	if(fixingPay_type == XL_TYPE_STRING)
		C_fixingPay_double = (double) LocalGetNumObjectId(C_fixingPay_str);
	else
	    C_fixingPay_double = ARM_NULL_OBJECT;

	if(fixing3_type == XL_TYPE_STRING)
		C_fixing3_double = (double) LocalGetNumObjectId(C_fixing3_str);
	else
	    C_fixing3_double = ARM_NULL_OBJECT;

	if(fixing4_type == XL_TYPE_STRING)
		C_fixing4_double = (double) LocalGetNumObjectId(C_fixing4_str);
	else
	    C_fixing4_double = ARM_NULL_OBJECT;

	if((DigitalcapOrFloorId = ARM_ConvCapOrFloor (C_DigitalcapOrFloor, C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}
	if((SpreadcapOrFloorId = ARM_ConvCapOrFloor (C_SpreadcapOrFloor, C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	liborType1Id = ARM_ConvIrType (C_liborType1);
	liborType2Id = ARM_ConvIrType (C_liborType2);
	liborType3Id = ARM_ConvIrType (C_liborType3);
	liborType4Id = ARM_ConvIrType (C_liborType4);

	dayCountId = ARM_ConvDayCount (C_daycount);
	
	if((resetFreqId = ARM_ConvFrequency (C_resetFreq, C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	if((payFreqId = ARM_ConvFrequency (C_payFreq, C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	resetTimingId = ARM_ConvPayResetRule (C_resetTiming);

	payTimingId = ARM_ConvPayResetRule (C_payTiming);

	if(C_ccy == "DEFAULT")
	{
		ccyId = ARM_NULL_OBJECT;
	}
	else
	{
		ccyId = LocalGetNumObjectId (C_ccy);
	}

	if(C_DigitalSpreadCorrel == "DEFAULT")
	{
		digitalSpreadcorrelId  = ARM_NULL_OBJECT;
	}
	else
	{
		digitalSpreadcorrelId = LocalGetNumObjectId (C_DigitalSpreadCorrel);
	}

	if ( digitalBarrier_type == XL_TYPE_STRING )
	{
	   C_digitalBarrier_double = (double) LocalGetNumObjectId(C_digitalBarrier_str);
	   digitalBarrier_type = 1;
	}
	else
	   digitalBarrier_type = 0;

	if ( spreadBarrier_type == XL_TYPE_STRING )
	{
	   C_spreadBarrier_double = (double) LocalGetNumObjectId(C_spreadBarrier_str);
	   spreadBarrier_type = 1;
	}
	else
	   spreadBarrier_type = 0;
	
	if ( shiftVolCorrel_type == XL_TYPE_STRING )
	{
	   C_ShiftVolCorrel_double = (double) LocalGetNumObjectId(C_ShiftVolCorrel_str);
	   shiftVolCorrel_type = 1;
	}
	else
	   shiftVolCorrel_type = 0;

	intRuleId = ARM_ConvIntRule(C_intRule);

	stubRuleId = ARM_ConvStubRule(C_stubRule);

	long retCode;
	long objId;
	CCString prevClass;

	CCString curClass = LOCAL_SPREAD_OPTION_CLASS;
	CCString stringId = GetLastCurCellEnvValue ();

	if(!stringId)
	{
		retCode = ARMLOCAL_CORRIDORDBLCONDITION(C_startDate,
												C_endDate,
												DigitalcapOrFloorId,
												SpreadcapOrFloorId,
											   (long)digitalBarrier_type,
											   C_digitalBarrier_double,
											   (long)spreadBarrier_type,
											   C_spreadBarrier_double,
											   payIndexId,
											   payMargin_type,
											   payMarginId,
											   (long)liborType1Id,
											   (long)liborType2Id,
											   (long)liborType3Id,
											   (long)liborType4Id,
											   C_weight1,
											   C_weight2,
											   C_weight3,
											   C_weight4,
											   (long)C_fixing1_double,
											   (long)C_fixing2_double,
											   (long)C_fixing3_double,
											   (long)C_fixing4_double,
											   (long)C_fixingPay_double,
											   dayCountId,
											   resetFreqId,
											   payFreqId,
											   resetTimingId,
											   payTimingId,
											   ccyId,
											   long(C_resetGap),
											   C_spreads[0],
											   C_spreads[1],
											   intRuleId,
											   stubRuleId,
											   resetCal,
											   payCal,
											   C_payIdxWeight,
											   (int)C_freezeFixing,
											   digitalSpreadcorrelId,
											   (long)shiftVolCorrel_type,
											   C_ShiftVolCorrel_double,
											   C_correlCorrectFlag,
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
			retCode = ARMLOCAL_CORRIDORDBLCONDITION(C_startDate,
												C_endDate,
												DigitalcapOrFloorId,
												SpreadcapOrFloorId,
											   (long)digitalBarrier_type,
											   C_digitalBarrier_double,
											   (long)spreadBarrier_type,
											   C_spreadBarrier_double,
											   payIndexId,
											   payMargin_type,
											   payMarginId,
											   (long)liborType1Id,
											   (long)liborType2Id,
											   (long)liborType3Id,
											   (long)liborType4Id,
											   C_weight1,
											   C_weight2,
											   C_weight3,
											   C_weight4,
											   (long)C_fixing1_double,
											   (long)C_fixing2_double,
											   (long)C_fixing3_double,
											   (long)C_fixing4_double,
											   (long)C_fixingPay_double,
											   dayCountId,
											   resetFreqId,
											   payFreqId,
											   resetTimingId,
											   payTimingId,
											   ccyId,
											   long(C_resetGap),
											   C_spreads[0],
											   C_spreads[1],
											   intRuleId,
											   stubRuleId,
											   resetCal,
											   payCal,
											   C_payIdxWeight,
											   (int)C_freezeFixing,
											   digitalSpreadcorrelId,
											   (long)shiftVolCorrel_type,
											   C_ShiftVolCorrel_double,
											   C_correlCorrectFlag,
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
			retCode = ARMLOCAL_CORRIDORDBLCONDITION(C_startDate,
												C_endDate,
												DigitalcapOrFloorId,
												SpreadcapOrFloorId,
											   (long)digitalBarrier_type,
											   C_digitalBarrier_double,
											   (long)spreadBarrier_type,
											   C_spreadBarrier_double,
											   payIndexId,
											   payMargin_type,
											   payMarginId,
											   (long)liborType1Id,
											   (long)liborType2Id,
											   (long)liborType3Id,
											   (long)liborType4Id,
											   C_weight1,
											   C_weight2,
											   C_weight3,
											   C_weight4,
											   (long)C_fixing1_double,
											   (long)C_fixing2_double,
											   (long)C_fixing3_double,
											   (long)C_fixing4_double,
											   (long)C_fixingPay_double,
											   dayCountId,
											   resetFreqId,
											   payFreqId,
											   resetTimingId,
											   payTimingId,
											   ccyId,
											   long(C_resetGap),
											   C_spreads[0],
											   C_spreads[1],
											   intRuleId,
											   stubRuleId,
											   resetCal,
											   payCal,
											   C_payIdxWeight,
											   (int)C_freezeFixing,
											   digitalSpreadcorrelId,
											   (long)shiftVolCorrel_type,
											   C_ShiftVolCorrel_double,
											   C_correlCorrectFlag,
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_ARM_CORRIDORDBLCONDITION" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}

__declspec(dllexport) LPXLOPER WINAPI Local_PXL_ARM_CORRIDORDBLCONDITION(LPXLOPER XL_startDate,
																		 LPXLOPER XL_endDate,
																		 LPXLOPER XL_DigitalcapOrFloor,
																		 LPXLOPER XL_SpreadcapOrFloor,
																		 LPXLOPER XL_digitalBarrier, 
																		 LPXLOPER XL_spreadBarrier,
																		 LPXLOPER XL_spreads,  //vect: spread1, spread2
																		 LPXLOPER XL_payIndexParams,//vect:payIndexObj,  itsMargin, itsWeight, itsFixedRate
																		 LPXLOPER XL_spreadIdxTypes,//vect: IndexType1, IndexType2, IndexType3
																		 LPXLOPER XL_Weights,//vect: Weight1,Weight2, Weight3
																		 LPXLOPER XL_prodParamDatas,//vect: dayCount, resetFreq, payFreq, resetTiming, payTiming, intRule, stubRule,resetGap 
																		 LPXLOPER XL_currency,
																		 LPXLOPER XL_fixing1,
																		 LPXLOPER XL_fixing2,
																		 LPXLOPER XL_fixing3,
																		 LPXLOPER XL_fixingPay,
																		 LPXLOPER XL_freezeFixing,
																		 LPXLOPER XL_DigitalSpreadCorrel,
																		 LPXLOPER XL_ShiftVolCorrel,
																		 LPXLOPER XL_fixing4)
{
	ADD_LOG("Local_PXL_ARM_CORRIDORDBLCONDITION");
	static XLOPER XL_result;
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();

	double C_startDate;
	double C_endDate;

	CCString C_DigitalcapOrFloor;
	long DigitalcapOrFloorId;

	CCString C_SpreadcapOrFloor;
	long SpreadcapOrFloorId;

	double C_digitalBarrier_double;
	CCString C_digitalBarrier_str;
	long digitalBarrier_type;
	double C_spreadBarrier_double;
	CCString C_spreadBarrier_str;
	long spreadBarrier_type;

	double C_ShiftVolCorrel_double;
	double C_Default_double = 0.0;
	CCString C_ShiftVolCorrel_str;
	long shiftVolCorrel_type;

	VECTOR<double> C_spreads;

	CCString C_payIndex;
	long payIndexId;

	CCString C_resetCal;
	CCString C_payCal;
	char* resetCal = NULL;
	char* payCal = NULL;
	
	CCString C_DigitalSpreadCorrel;
	long digitalSpreadcorrelId;

	CCString C_payMargin_str;
	long payMargin_type;
	double payMarginId;

	double C_payIdxWeight;
	
//	CCString C_payIdxFixedRate_str;
//	long payIdxFixedRate_type;
//	double payIdxFixedRateId;

	VECTOR<CCString> payIndexParams(3);
   	VECTOR<CCString> C_payIndexParams;
	VECTOR<CCString> C_payIndexParamsDF(3);
	C_payIndexParamsDF[0] = payIndexParams[0] = "DEFAULT";// pay index
	C_payIndexParamsDF[1] = payIndexParams[1] = "DEFAULT";// pay margin referenceValue
	C_payIndexParamsDF[2] = payIndexParams[2] = "1";// pay weight
//	C_payIndexParamsDF[3] = payIndexParams[3] = "DEFAULT"; // pay fixed rate

	
	CCString C_liborType1;
	long liborType1Id;
	CCString C_liborType2;
	long liborType2Id;
	CCString C_liborType3;
	long liborType3Id;
	CCString C_liborType4;
	long liborType4Id;
	VECTOR<CCString> liborTypes(4);
   	VECTOR<CCString> C_liborTypes;
	VECTOR<CCString> C_liborTypesDF(4);
	C_liborTypesDF[0] = liborTypes[0] = "LIBOR3M";
	C_liborTypesDF[1] = liborTypes[1] = "LIBOR3M";
	C_liborTypesDF[2] = liborTypes[2] = "LIBOR3M";
	C_liborTypesDF[3] = liborTypes[3] = "LIBOR3M";

	double C_weight1, C_weight2, C_weight3, C_weight4;
	VECTOR<double> weights(4,1.0);
	VECTOR<double> C_weights;
	VECTOR<double> C_weights_df(4, 1.0);

	VECTOR<CCString> prodParamDatas(11);
	VECTOR<CCString> C_prodParamDatas;
	VECTOR<CCString> C_prodParamDatasDF(11);
	C_prodParamDatasDF[0] = prodParamDatas[0] = "30/360"; // dayCount
	C_prodParamDatasDF[1] = prodParamDatas[1] = "A"; // resetFreq
	C_prodParamDatasDF[2] = prodParamDatas[2] = "A"; // payFreq
	C_prodParamDatasDF[3] = prodParamDatas[3] = "ADV"; // resetTiming
	C_prodParamDatasDF[4] = prodParamDatas[4] = "ARR"; // payTiming
	C_prodParamDatasDF[5] = prodParamDatas[5] = "ADJ"; // intRule
	C_prodParamDatasDF[6] = prodParamDatas[6] = "SS";  // stubRule
	C_prodParamDatasDF[7] = prodParamDatas[7] = "10000.0";// resetGap
	C_prodParamDatasDF[8] = prodParamDatas[8] = "N";	// no correlation corrector
	C_prodParamDatasDF[9] = prodParamDatas[9] = "";	// resetCal
	C_prodParamDatasDF[10] = prodParamDatas[10] = "";	// payCal

	CCString C_daycount;
	long dayCountId;

	CCString C_resetFreq;
	long resetFreqId;

	CCString C_payFreq;
	long payFreqId;

	CCString C_resetTiming;
	long resetTimingId;
	
	CCString C_payTiming;
	long payTimingId;

	CCString C_intRule;
    long intRuleId;
	
	CCString C_stubRule;
    long stubRuleId;

    double C_resetGap;
	bool C_correlCorrectFlag;

	CCString C_ccy;
    long ccyId;

   	double C_fixing1_double;
   	VECTOR<double> C_fixing1_double_vect;
	CCString C_fixing1_str;
	long fixing1_type;
	VECTOR<double> C_fixing1_default(0.0);

   	double C_fixing2_double;
   	VECTOR<double> C_fixing2_double_vect;
	CCString C_fixing2_str;
	long fixing2_type;
	VECTOR<double> C_fixing2_default(0.0);

   	double C_fixing3_double;
   	VECTOR<double> C_fixing3_double_vect;
	CCString C_fixing3_str;
	long fixing3_type;
	VECTOR<double> C_fixing3_default(0.0);

   	double C_fixing4_double;
   	VECTOR<double> C_fixing4_double_vect;
	CCString C_fixing4_str;
	long fixing4_type;
	VECTOR<double> C_fixing4_default(0.0);

	double C_fixingPay_double;
   	VECTOR<double> C_fixingPay_double_vect;
	CCString C_fixingPay_str;
	long fixingPay_type;
	VECTOR<double> C_fixingPay_default(0.0);

	double C_freezeFixing;
	double C_freezeFixing_default = 0.0;

	// error
	static int error;
	static char* reason = "";


	XL_readNumCell(XL_startDate,C_startDate," ARM_ERR: start date: date expected",C_result);
	XL_readNumCell(XL_endDate,C_endDate," ARM_ERR: end date: date expected",C_result);
	XL_readStrCell(XL_DigitalcapOrFloor,C_DigitalcapOrFloor," ARM_ERR: cap or floor: string expected",C_result);
	XL_readStrCell(XL_SpreadcapOrFloor,C_SpreadcapOrFloor," ARM_ERR: cap or floor: string expected",C_result);
	XL_readStrOrNumCell(XL_digitalBarrier,C_digitalBarrier_str,C_digitalBarrier_double,digitalBarrier_type," ARM_ERR: Rate/Spread#1 Barrier: string or numeric expected",C_result);
	XL_readStrOrNumCell(XL_spreadBarrier,C_spreadBarrier_str,C_spreadBarrier_double,spreadBarrier_type," ARM_ERR: Spread#2 Barrier: string or numeric expected",C_result);
	XL_readNumVector(XL_spreads,C_spreads," ARM_ERR: spread vector: array of numeric expected", C_result);
	XL_readStrVectorWD(XL_payIndexParams,C_payIndexParams,C_payIndexParamsDF," ARM_ERR: payIndex Params: array expected",DOUBLE_TYPE,C_result);
	XL_readStrVectorWD(XL_spreadIdxTypes,C_liborTypes,C_liborTypesDF," ARM_ERR: spread indexes : array of string expected",DOUBLE_TYPE,C_result);
	XL_readNumVectorWD(XL_Weights,C_weights,C_weights_df," ARM_ERR: weights: array of numeric expected", C_result);
	XL_readStrVectorWD(XL_prodParamDatas,C_prodParamDatas,C_prodParamDatasDF," ARM_ERR: daycount,resetFreq,payFreq,resetTiming,payTiming,intRule,stubRule,resetGap,correlCorrector: array of string expected",DOUBLE_TYPE,C_result);
	XL_readStrCellWD(XL_currency,C_ccy,"DEFAULT"," ARM_ERR: currency : object expected",C_result);
	XL_readStrOrNumVectorWD(XL_fixing1,C_fixing1_str,C_fixing1_double_vect,C_fixing1_default,fixing1_type," ARM_ERR: fixing1 taux : array or string expected",C_result);
	XL_readStrOrNumVectorWD(XL_fixing2,C_fixing2_str,C_fixing2_double_vect,C_fixing2_default,fixing2_type," ARM_ERR: fixing2 taux : array or string expected",C_result);
	XL_readStrOrNumVectorWD(XL_fixing3,C_fixing3_str,C_fixing3_double_vect,C_fixing3_default,fixing3_type," ARM_ERR: fixing3 taux : array or string expected",C_result);
	XL_readStrOrNumVectorWD(XL_fixing4,C_fixing4_str,C_fixing4_double_vect,C_fixing4_default,fixing4_type," ARM_ERR: fixing4 taux : array or string expected",C_result);
	XL_readStrOrNumVectorWD(XL_fixingPay,C_fixingPay_str,C_fixingPay_double_vect,C_fixingPay_default,fixingPay_type," ARM_ERR: fixingPay taux : array or string expected",C_result);
	XL_readNumCellWD(XL_freezeFixing,C_freezeFixing,C_freezeFixing_default," ARM_ERR: freeze fixing : numeric expected",C_result);
	XL_readStrCellWD(XL_DigitalSpreadCorrel,C_DigitalSpreadCorrel,"DEFAULT"," ARM_ERR: VolCorrel : object expected",C_result);
	XL_readStrOrNumCellWD(XL_ShiftVolCorrel,C_ShiftVolCorrel_str, C_ShiftVolCorrel_double, C_Default_double, shiftVolCorrel_type,
			   " ARM_ERR: ShiftCorrel : object expected",C_result);


	if(C_spreads.size()<2)
	{
		C_result.setMsg ("ARM_ERR: Verify the Vector of Spread");
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}
	
	// traitement de payIndexType_itsWeight
	int paySize = C_payIndexParams.size() < 3 ? C_payIndexParams.size() : 3;
	for(int i=0; i<paySize; i++)
	{
		if( strcmp(C_payIndexParams[i], "DEFAULT") != 0)// fill in the value not by default
			payIndexParams[i] = C_payIndexParams[i];
	}

	C_payIndex = payIndexParams[0];
	C_payMargin_str = payIndexParams[1];
	C_payIdxWeight = atof(payIndexParams[2]);
//	C_payIdxFixedRate_str = payIndexParams[3];

	if(C_payIndex == "DEFAULT")
	{
		payIndexId = ARM_NULL_OBJECT;
	}
	else
	{
		payIndexId = LocalGetNumObjectId(C_payIndex);
	}

	if(C_payMargin_str == "DEFAULT")
	{
		payMarginId = 0.0;
		payMargin_type = 0;//nombre
	}
	else
	{
		payMarginId = LocalGetNumObjectId(C_payMargin_str);

		if( (int)payMarginId == -1 ) // non objet
		{
			payMarginId = atof(C_payMargin_str);
			payMargin_type = 0;
		}
		else
			payMargin_type = 1;//objet
	}

/*	if(C_payIdxFixedRate_str == "DEFAULT")
	{
		payIdxFixedRateId = 0.0;
		payIdxFixedRate_type = 0;//nombre
	}
	else
	{
		payIdxFixedRateId = LocalGetNumObjectId(C_payIdxFixedRate_str);

		if( (int)payIdxFixedRateId == -1 ) // non objet
		{
			payIdxFixedRateId = atof(C_payIdxFixedRate_str);
			payIdxFixedRate_type = 0;
		}
		else
			payIdxFixedRate_type = 1;//objet
	}
*/
	//traitement de indexes et des weights
	int IdxSize = C_liborTypes.size();
	bool isRate=true;
	if(IdxSize>3)
	{
		isRate=false; // 2 spreads
		IdxSize=4;
	}
	for (i =0; i< IdxSize; i++)
	{
		if(strcmp(C_liborTypes[i],"") != 0)
			liborTypes[i] = C_liborTypes[i];
	}
	if(isRate)
	{
		C_liborType1	= liborTypes[3]; // 2nd rate of 1st spread = default value
		C_liborType2	= liborTypes[0]; // Rate or 2st rate of 1st spread
		C_liborType3	= liborTypes[1]; // 1st rate of 2nd spread
		C_liborType4	= liborTypes[2]; // 2nd rate of 2nd spread
		C_weight1		= 0.0;			 // no use of 2nd rate of 1st spread
		C_weight2		= C_weights[0];
		C_weight3		= C_weights[1];
		C_weight4		= C_weights[2];
	}
	else
	{
		C_liborType1	= liborTypes[0]; // 1st rate of 1st spread
		C_liborType2	= liborTypes[1]; // 2st rate of 1st spread
		C_liborType3	= liborTypes[2]; // 1st rate of 2nd spread
		C_liborType4	= liborTypes[3]; // 2nd rate of 2nd spread
		C_weight1		= C_weights[0];
		C_weight2		= C_weights[1];
		C_weight3		= C_weights[2];
		C_weight4		= C_weights[3];
	}

	// traitement de prodParamDatas
	int prodSize = C_prodParamDatas.size() < 11 ? C_prodParamDatas.size() : 11;
	for (i =0; i<prodSize; i++)
	{
		if( strcmp(C_prodParamDatas[i],"DEFAULT")!= 0 )
			prodParamDatas[i] = C_prodParamDatas[i];
	}
	// dayCount
	if (prodParamDatas.size() >= 1)
		C_daycount    = prodParamDatas[0];
	else
		C_daycount    = C_prodParamDatasDF[0];
	// resetFreq
	if (prodParamDatas.size() >= 2)
		C_resetFreq   = prodParamDatas[1];
	else
		C_resetFreq   = C_prodParamDatasDF[1];
	// payFreq
	if (prodParamDatas.size() >= 3)
		C_payFreq     = prodParamDatas[2];
	else
		C_payFreq     = C_prodParamDatasDF[2];
	// resetTiming
	if (prodParamDatas.size() >= 4)
		C_resetTiming = prodParamDatas[3];
	else
		C_resetTiming = C_prodParamDatasDF[3];
	// payTiming
	if (prodParamDatas.size() >= 5)
		C_payTiming   = prodParamDatas[4];
	else
		C_payTiming   = C_prodParamDatasDF[4];
	// intRule
	if (prodParamDatas.size() >= 6)
		C_intRule     = prodParamDatas[5];
	else
		C_intRule     = C_prodParamDatasDF[5];
	// stubRule
	if (prodParamDatas.size() >= 7)
		C_stubRule    = prodParamDatas[6];
	else
		C_stubRule    = C_prodParamDatasDF[6];
	// resetGap
	if (prodParamDatas.size() >= 8)
		C_resetGap    = atof(prodParamDatas[7]);
	else
		C_resetGap    = atof(C_prodParamDatasDF[7]);
	// correl corrector
	if (prodParamDatas.size() >= 9)
		C_correlCorrectFlag    = (prodParamDatas[8]=="Y" || prodParamDatas[8]=="y");
	else
		C_correlCorrectFlag    = (C_prodParamDatasDF[8]=="Y" || C_prodParamDatasDF[8]=="y");
	// resetCal
	if (prodParamDatas.size() >= 10)
		C_resetCal = prodParamDatas[9];
	else
		C_resetCal   = C_prodParamDatasDF[9];
	// payCal
	if (prodParamDatas.size() >= 11)
		C_payCal   = prodParamDatas[10];
	else
		C_payCal   = C_prodParamDatasDF[10];
	// resetCal
	if (prodParamDatas.size() >= 10)
		C_resetCal = prodParamDatas[9];
	else
		C_resetCal   = C_prodParamDatasDF[9];
	// payCal
	if (prodParamDatas.size() >= 11)
		C_payCal   = prodParamDatas[10];
	else
		C_payCal   = C_prodParamDatasDF[10];

	if(!(C_resetCal == ""))
		resetCal = C_resetCal.c_str();
	if(!(C_payCal == ""))
		payCal = C_payCal.c_str();
	
	if(fixing1_type == XL_TYPE_STRING)
		C_fixing1_double = (double) LocalGetNumObjectId(C_fixing1_str);
	else
	    C_fixing1_double = ARM_NULL_OBJECT;

	if(fixing2_type == XL_TYPE_STRING)
		C_fixing2_double = (double) LocalGetNumObjectId(C_fixing2_str);
	else
	    C_fixing2_double = ARM_NULL_OBJECT;

	if(fixingPay_type == XL_TYPE_STRING)
		C_fixingPay_double = (double) LocalGetNumObjectId(C_fixingPay_str);
	else
	    C_fixingPay_double = ARM_NULL_OBJECT;

	if(fixing3_type == XL_TYPE_STRING)
		C_fixing3_double = (double) LocalGetNumObjectId(C_fixing3_str);
	else
	    C_fixing3_double = ARM_NULL_OBJECT;

	if(fixing4_type == XL_TYPE_STRING)
		C_fixing4_double = (double) LocalGetNumObjectId(C_fixing4_str);
	else
	    C_fixing4_double = ARM_NULL_OBJECT;

	if((DigitalcapOrFloorId = ARM_ConvCapOrFloor (C_DigitalcapOrFloor, C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	if((SpreadcapOrFloorId = ARM_ConvCapOrFloor (C_SpreadcapOrFloor, C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	liborType1Id = ARM_ConvIrType (C_liborType1);
	liborType2Id = ARM_ConvIrType (C_liborType2);
	liborType3Id = ARM_ConvIrType (C_liborType3);
	liborType4Id = ARM_ConvIrType (C_liborType4);

	dayCountId = ARM_ConvDayCount (C_daycount);
	
	if((resetFreqId = ARM_ConvFrequency (C_resetFreq, C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	if((payFreqId = ARM_ConvFrequency (C_payFreq, C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	resetTimingId = ARM_ConvPayResetRule (C_resetTiming);

	payTimingId = ARM_ConvPayResetRule (C_payTiming);

	if(C_ccy == "DEFAULT")
	{
		ccyId = ARM_NULL_OBJECT;
	}
	else
	{
		ccyId = LocalGetNumObjectId (C_ccy);
	}

	if(C_DigitalSpreadCorrel == "DEFAULT")
	{
		digitalSpreadcorrelId  = ARM_NULL_OBJECT;
	}
	else
	{
		digitalSpreadcorrelId = LocalGetNumObjectId (C_DigitalSpreadCorrel);
	}

	if ( digitalBarrier_type == XL_TYPE_STRING )
	{
	   C_digitalBarrier_double = (double) LocalGetNumObjectId(C_digitalBarrier_str);
	   digitalBarrier_type = 1;
	}
	else
	   digitalBarrier_type = 0;

	if ( spreadBarrier_type == XL_TYPE_STRING )
	{
	   C_spreadBarrier_double = (double) LocalGetNumObjectId(C_spreadBarrier_str);
	   spreadBarrier_type = 1;
	}
	else
	   spreadBarrier_type = 0;

	if ( shiftVolCorrel_type == XL_TYPE_STRING )
	{
	   C_ShiftVolCorrel_double = (double) LocalGetNumObjectId(C_ShiftVolCorrel_str);
	   shiftVolCorrel_type = 1;
	}
	else
	   shiftVolCorrel_type = 0;
	
	intRuleId = ARM_ConvIntRule(C_intRule);

	stubRuleId = ARM_ConvStubRule(C_stubRule);

	long retCode;
	long objId;

	CCString curClass = LOCAL_SPREAD_OPTION_CLASS;
	CCString stringId = GetLastCurCellEnvValue ();

	retCode = ARMLOCAL_CORRIDORDBLCONDITION(C_startDate,
											C_endDate,
											DigitalcapOrFloorId,
											SpreadcapOrFloorId,
											(long)digitalBarrier_type,
											C_digitalBarrier_double,
											(long)spreadBarrier_type,
											C_spreadBarrier_double,
											payIndexId,
											payMargin_type,
											payMarginId,
										   (long)liborType1Id,
										   (long)liborType2Id,
										   (long)liborType3Id,
										   (long)liborType4Id,
										   C_weight1,
										   C_weight2,
										   C_weight3,
										   C_weight4,
										   (long)C_fixing1_double,
										   (long)C_fixing2_double,
										   (long)C_fixing3_double,
										   (long)C_fixing4_double,
										   (long)C_fixingPay_double,
											dayCountId,
											resetFreqId,
											payFreqId,
											resetTimingId,
											payTimingId,
											ccyId,
											long(C_resetGap),
											C_spreads[0],
											C_spreads[1],
											intRuleId,
											stubRuleId,
											resetCal,
											payCal,
											C_payIdxWeight,
											(int)C_freezeFixing,
											digitalSpreadcorrelId,
											(long)shiftVolCorrel_type,
											C_ShiftVolCorrel_double,
											C_correlCorrectFlag,
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_PXL_ARM_CORRIDORDBLCONDITION" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


__declspec(dllexport) LPXLOPER WINAPI Local_ARM_CHECK_CORRELS(LPXLOPER XL_expiries,
															  LPXLOPER XL_tenors,
															  LPXLOPER XL_model)
{
	ADD_LOG("Local_ARM_CHECK_CORRELS");
	static XLOPER XL_result;
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();

	// error
	static int error;
	static char* reason = "";


	VECTOR<double> C_expiries;
	XL_readNumVector(XL_expiries,C_expiries," ARM_ERR: expiries vector: array of numeric expected", C_result);
   	VECTOR<CCString> C_tenors;
   	VECTOR<CCString> C_tenorsDef(4,"CMS1Y");
	XL_readStrVectorWD(XL_tenors,C_tenors,C_tenorsDef," ARM_ERR: tenors vector: array of numeric expected",XL_TYPE_STRING,C_result);
	VECTOR<long> C_tenorTypes;
	for(size_t i=0;i<C_tenors.size();++i)
		C_tenorTypes.push_back(ARM_ConvIrType (C_tenors[i]));
	CCString C_model;
	XL_readStrCell(XL_model,C_model," ARM_ERR: model : object expected",C_result);
	long C_modelId = LocalGetNumObjectId(C_model);

	vector<double> C_DataResult;
	long C_rows,C_cols;
	if(ARMLOCAL_CHECK_CORRELS(C_expiries,C_tenorTypes,C_modelId,C_DataResult,C_rows,C_cols,C_result) == ARM_OK)
	{
		const int additionalLinesNb = 100;
		bool fillWithBlank = true;
		FreeCurCellContent ();
		XL_writeNumMatrixSizeWithOptions( XL_result, C_DataResult, C_rows, C_cols, " ARM_ERR: Could not get result data", C_result, additionalLinesNb, fillWithBlank );
	}
	else
	{
		ARM_ERR();
	}

	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_ARM_CHECK_CORRELS" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}

__declspec(dllexport) LPXLOPER WINAPI Local_ARM_CHECK_RA2_CORRELS_STATUS(LPXLOPER XL_RA2Id,
																		 LPXLOPER XL_ModelId)
{
	ADD_LOG("Local_ARM_CHECK_RA2_CORRELS_STATUS");
	static XLOPER XL_result;
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();

	// error
	static int error;
	static char* reason = "";


	CCString C_RA2;
	XL_readStrCell(XL_RA2Id,C_RA2," ARM_ERR: security : object expected",C_result);
	long C_RA2Id = LocalGetNumObjectId(C_RA2);

	CCString C_model;
	CCString C_modelD="DEFAULT";
	XL_readStrCellWD(XL_ModelId,C_model,C_modelD," ARM_ERR: model : object expected",C_result);
	long C_modelId=ARM_NULL_OBJECT;
	if(!(C_model == "DEFAULT"))
		C_modelId = LocalGetNumObjectId(C_model);

	vector<double> C_DataResult;
	long C_rows,C_cols;
	if(ARMLOCAL_CHECK_RA2_CORRELS_STATUTS(C_RA2Id,C_modelId,C_DataResult,C_rows,C_cols,C_result) == ARM_OK)
	{
		if(C_modelId==ARM_NULL_OBJECT)
		{
			FreeCurCellErr ();
			XL_result.xltype = xltypeStr;
			XL_result.val.str = XL_StrC2StrPascal (C_result.getString() );
			XL_result.xltype |= xlbitDLLFree;
		}
		else
		{
			const int additionalLinesNb = 100;
			bool fillWithBlank = true;
			FreeCurCellContent ();
			XL_writeNumMatrixSizeWithOptions( XL_result, C_DataResult, C_rows, C_cols, " ARM_ERR: Could not get result data", C_result, additionalLinesNb, fillWithBlank );
		}
	}
	else
	{
		ARM_ERR();
	}

	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_ARM_CHECK_RA2_CORRELS" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}



__declspec(dllexport) LPXLOPER WINAPI Local_ARM_SPREADDIGITALFLT(LPXLOPER XL_startDate,
															     LPXLOPER XL_endDate,
															     LPXLOPER XL_capOrFloor,
															     LPXLOPER XL_strike,
															     LPXLOPER XL_spread,//Vector
															     LPXLOPER XL_payoffLiborType,
															     LPXLOPER XL_liborType1,
															     LPXLOPER XL_liborType2,
																 LPXLOPER XL_weight,//Vector Weights, slopeflag ans cptStrikeMethod
															     LPXLOPER XL_daycount,
															     LPXLOPER XL_resetFreq,
															     LPXLOPER XL_payFreq,
															     LPXLOPER XL_resetTiming,
															     LPXLOPER XL_payTiming,
															     LPXLOPER XL_currency,
                                                                 LPXLOPER XL_resetGap,
																 LPXLOPER XL_intRule,
																 LPXLOPER XL_stubRule,
																 LPXLOPER XL_fixing1,
																 LPXLOPER XL_fixing2)

{
	ADD_LOG("Local_ARM_SPREADDIGITALFLT");
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

	CCString C_capOrFloor;
	long capOrFloorId;

	double C_strike_double;
	CCString C_strike_str;
	long strike_type;
	
	CCString C_payoffLiborType;
	long payoffLiborTypeId;
	
	CCString C_liborType1;
	long liborType1Id;
	CCString C_liborType2;
	long liborType2Id;

	double C_weight1, C_weight2, C_cptStrikeMethod;
	VECTOR<double> C_weight;
	VECTOR<double> C_weight_default(1.0);

	long C_slopeFlag;
	long C_slopeFlagDefault = 1;
	double C_cptStrikeMethod_default = 1.0;

	double C_computedFormula;
	double C_computedFormula_default = 1.0;

	CCString C_daycount;
	long dayCountId;

	CCString C_resetFreq;
	long resetFreqId;
	CCString C_payFreq;
	long payFreqId;

	CCString C_resetTiming;
	long resetTimingId;
	CCString C_payTiming;
	long payTimingId;

	CCString C_ccy;
    long ccyId;

    double C_resetGap;
    double C_resetGap_default = 10000.0;

	double C_spread1;
	double C_spread2;
	VECTOR<double> C_spread;

   	double C_fixing1_double;
   	VECTOR<double> C_fixing1_double_vect;
	CCString C_fixing1_str;
	long fixing1_type;
	VECTOR<double> C_fixing1_default;

   	double C_fixing2_double;
   	VECTOR<double> C_fixing2_double_vect;
	CCString C_fixing2_str;
	long fixing2_type;
	VECTOR<double> C_fixing2_default;

	CCString C_intRule;
    long intRuleId;

	CCString C_stubRule;
    long stubRuleId;

	// error
	static int error;
	static char* reason = "";

	XL_readNumCell(XL_startDate,C_startDate," ARM_ERR: start date: date expected",C_result);
	XL_readNumCell(XL_endDate,C_endDate," ARM_ERR: end date: date expected",C_result);
	XL_readStrCell(XL_capOrFloor,C_capOrFloor," ARM_ERR: cap or floor: string expected",C_result);
	XL_readStrOrNumCell(XL_strike,C_strike_str,C_strike_double,strike_type," ARM_ERR: strike: string or numeric expected",C_result);
	
	XL_readNumVector(XL_spread,C_spread," ARM_ERR: spread vector: array of numeric expected", C_result);

	XL_readStrCellWD(XL_payoffLiborType,C_payoffLiborType,"LIBOR3M"," ARM_ERR: payoff libor type  : string expected",C_result);

	XL_readStrCellWD(XL_liborType1,C_liborType1,"LIBOR3M"," ARM_ERR: libor type 1 : string expected",C_result);
	XL_readStrCellWD(XL_liborType2,C_liborType2,"LIBOR3M"," ARM_ERR: libor type 2 : string expected",C_result);
	XL_readNumVectorWD(XL_weight,C_weight,C_weight_default," ARM_ERR: weights_slopeflag: array of numeric expected", C_result);
	XL_readStrCellWD(XL_daycount,C_daycount,"30/360"," ARM_ERR: daycount : string expected",C_result);
	XL_readStrCellWD(XL_resetFreq,C_resetFreq,"A"," ARM_ERR: reset frequency : string expected",C_result);
	XL_readStrCellWD(XL_payFreq,C_payFreq,"A"," ARM_ERR: payment frequency : string expected",C_result);
	XL_readStrCellWD(XL_resetTiming,C_resetTiming,"ADV"," ARM_ERR: reset timing : string expected",C_result);
	XL_readStrCellWD(XL_payTiming,C_payTiming,"ARR"," ARM_ERR: payment timing : string expected",C_result);
	XL_readStrCellWD(XL_currency,C_ccy,"DEFAULT"," ARM_ERR: currency : object expected",C_result);
	XL_readNumCellWD(XL_resetGap,C_resetGap,C_resetGap_default," ARM_ERR: resetGap : numeric expected",C_result);
	XL_readStrCellWD(XL_intRule,C_intRule,"DEFAULT"," ARM_ERR: intRule : string expected",C_result);
	XL_readStrCellWD(XL_stubRule,C_stubRule,"DEFAULT"," ARM_ERR: stubRule : string expected",C_result);

	XL_readStrOrNumVectorWD(XL_fixing1,C_fixing1_str,C_fixing1_double_vect,C_fixing1_default,fixing1_type," ARM_ERR: fixing taux : array of numeric expected",C_result);
	XL_readStrOrNumVectorWD(XL_fixing2,C_fixing2_str,C_fixing2_double_vect,C_fixing2_default,fixing2_type," ARM_ERR: fixing taux : array of numeric expected",C_result);

	if(C_weight.size()<2)
	{
		C_result.setMsg ("ARM_ERR: Verify the Vector of Weight");
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}
	
	C_weight1 = C_weight[0];
	C_weight2 = C_weight[1];

	if(C_weight.size() == 2)
	{
		C_slopeFlag = C_slopeFlagDefault;
		C_cptStrikeMethod = C_cptStrikeMethod_default;
		C_computedFormula = C_computedFormula_default;
	}
	else if(C_weight.size() == 3)
	{
		C_slopeFlag = C_weight[2];
		C_cptStrikeMethod = C_cptStrikeMethod_default;
		C_computedFormula = C_computedFormula_default;
	}
	else if(C_weight.size() == 4)
	{
		C_slopeFlag = C_weight[2];
		C_cptStrikeMethod = C_weight[3];
		C_computedFormula = C_computedFormula_default;
	}
	else
	{
		C_slopeFlag = C_weight[2];
		C_cptStrikeMethod = C_weight[3];
		C_computedFormula = C_weight[4];
	}
	
	if(C_spread.size()<2)
	{
		C_result.setMsg ("ARM_ERR: Verify the Vector of Spread");
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}
	C_spread1 = C_spread[0];
	C_spread2 = C_spread[1];
	
	if ((XL_fixing1->xltype == xltypeMissing) || (XL_fixing1->xltype == xltypeNil))
	{
		fixing1_type = 1;
		C_fixing1_double = ARM_NULL_OBJECT;
    }
	else if ( fixing1_type == XL_TYPE_STRING )
	{
		C_fixing1_double = (double) LocalGetNumObjectId(C_fixing1_str);
		fixing1_type = 1;
	}
	else
	{
	    fixing1_type = 0;
	}

	if ((XL_fixing2->xltype == xltypeMissing) || (XL_fixing2->xltype == xltypeNil))
	{
		fixing2_type = 1;
		C_fixing2_double = ARM_NULL_OBJECT;
    }
	else if ( fixing2_type == XL_TYPE_STRING )
	{
		C_fixing2_double = (double) LocalGetNumObjectId(C_fixing2_str);
		fixing2_type = 1;
	}
	else
	{
	    fixing2_type = 0;
	}

	if((capOrFloorId = ARM_ConvCapOrFloor (C_capOrFloor, C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	payoffLiborTypeId = ARM_ConvIrType (C_payoffLiborType);

	liborType1Id = ARM_ConvIrType (C_liborType1);

	liborType2Id = ARM_ConvIrType (C_liborType2);

	dayCountId = ARM_ConvDayCount (C_daycount);
	
	if((resetFreqId = ARM_ConvFrequency (C_resetFreq, C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	if((payFreqId = ARM_ConvFrequency (C_payFreq, C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	resetTimingId = ARM_ConvPayResetRule (C_resetTiming);

	payTimingId = ARM_ConvPayResetRule (C_payTiming);

	if(C_ccy == "DEFAULT")
	{
		ccyId = ARM_NULL_OBJECT;
	}
	else
	{
		ccyId = LocalGetNumObjectId (C_ccy);
	}

	if ( strike_type == XL_TYPE_STRING )
	{
	   C_strike_double = (double) LocalGetNumObjectId(C_strike_str);

	   strike_type = 1;
	}
	else
	{
	   strike_type = 0;
	}
	
	if(C_intRule == "DEFAULT")
	{
		C_intRule = "ADJ";
	}
	intRuleId = ARM_ConvIntRule(C_intRule);

	if(C_stubRule == "DEFAULT")
	{
		C_stubRule = "SS";
	}
	stubRuleId = ARM_ConvStubRule(C_stubRule);

	long retCode;
	long objId;
	CCString prevClass;
	
	CCString curClass = LOCAL_SPREAD_OPTION_CLASS;
	CCString stringId = GetLastCurCellEnvValue ();

	if(!stringId)
	{
		retCode = ARMLOCAL_SPREADDIGITALFLT( C_startDate,
											 C_endDate,
											 capOrFloorId,
											 (long) strike_type,
											 C_strike_double,
											 (long)payoffLiborTypeId,
											 (long)liborType1Id,
											 (long)liborType2Id,
											 C_weight1,
											 C_weight2,
											 dayCountId,
											 resetFreqId,
											 payFreqId,
											 resetTimingId,
											 payTimingId,
											 ccyId,
											 long(C_resetGap),
											 C_spread1,
											 C_spread2,
											 (long) fixing1_type,
											 C_fixing1_double,
											 C_fixing1_double_vect,
											 (long) fixing2_type,
											 C_fixing2_double,
											 C_fixing2_double_vect,
											 intRuleId,
											 stubRuleId,
											 C_slopeFlag,
											 (int) C_cptStrikeMethod,
											 (int) C_computedFormula,
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
			retCode = ARMLOCAL_SPREADDIGITALFLT( C_startDate,
												 C_endDate,
												 capOrFloorId,
												 (long) strike_type,
												 C_strike_double,
												 (long)payoffLiborTypeId,
												 (long)liborType1Id,
												 (long)liborType2Id,
												 C_weight1,
												 C_weight2,
												 dayCountId,
												 resetFreqId,
												 payFreqId,
												 resetTimingId,
												 payTimingId,
												 ccyId,
												 long(C_resetGap),
												 C_spread1,
												 C_spread2,
												 (long) fixing1_type,
												 C_fixing1_double,
												 C_fixing1_double_vect,
												 (long) fixing2_type,
												 C_fixing2_double,
												 C_fixing2_double_vect,
												 intRuleId,
												 stubRuleId,
												 C_slopeFlag,
												 (int) C_cptStrikeMethod,
												 (int) C_computedFormula,
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
			retCode = ARMLOCAL_SPREADDIGITALFLT(C_startDate,
											 C_endDate,
											 capOrFloorId,
											 (long) strike_type,
											 C_strike_double,
										     (long)payoffLiborTypeId,
											 (long)liborType1Id,
											 (long)liborType2Id,
											 C_weight1,
											 C_weight2,
											 dayCountId,
											 resetFreqId,
											 payFreqId,
											 resetTimingId,
											 payTimingId,
											 ccyId,
                                             long(C_resetGap),
											 C_spread1,
											 C_spread2,
											 (long) fixing1_type,
											 C_fixing1_double,
											 C_fixing1_double_vect,
											 (long) fixing2_type,
											 C_fixing2_double,
											 C_fixing2_double_vect,
											 intRuleId,
											 stubRuleId,
											 C_slopeFlag,
											 (int) C_cptStrikeMethod,
											 (int) C_computedFormula,
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_ARM_SPREADDIGITALFLT" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


__declspec(dllexport) LPXLOPER WINAPI Local_PXL_ARM_SPREADDIGITALFLT(LPXLOPER XL_startDate,
																	 LPXLOPER XL_endDate,
																	 LPXLOPER XL_capOrFloor,
																	 LPXLOPER XL_strike,
																	 LPXLOPER XL_spread,//Vector
																	 LPXLOPER XL_payoffLiborType,
																	 LPXLOPER XL_liborType1,
																	 LPXLOPER XL_liborType2,
																	 LPXLOPER XL_weight,//Vector Weights, slopeflag and cptStrikeMethod
																	 LPXLOPER XL_daycount,
																	 LPXLOPER XL_resetFreq,
																	 LPXLOPER XL_payFreq,
																	 LPXLOPER XL_resetTiming,
																	 LPXLOPER XL_payTiming,
																	 LPXLOPER XL_currency,
																	 LPXLOPER XL_resetGap,
																	 LPXLOPER XL_intRule,
																	 LPXLOPER XL_stubRule,
																	 LPXLOPER XL_fixing1,
																	 LPXLOPER XL_fixing2)
{
	ADD_LOG("Local_PXL_ARM_SPREADDIGITALFLT");
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

	CCString C_capOrFloor;
	long capOrFloorId;

	double C_strike_double;
	CCString C_strike_str;
	long strike_type;
	
	CCString C_payoffLiborType;
	long payoffLiborTypeId;

	CCString C_liborType1;
	long liborType1Id;
	CCString C_liborType2;
	long liborType2Id;

	double C_weight1, C_weight2, C_cptStrikeMethod;
	VECTOR<double> C_weight;
	VECTOR<double> C_weight_default(1.0);

	long C_slopeFlag;
	long C_slopeFlagDefault = 1;
	double C_cptStrikeMethod_default = 1.0;

	double C_computedFormula;
	double C_computedFormula_default = 1.0;

	CCString C_daycount;
	long dayCountId;

	CCString C_resetFreq;
	long resetFreqId;
	CCString C_payFreq;
	long payFreqId;

	CCString C_resetTiming;
	long resetTimingId;
	CCString C_payTiming;
	long payTimingId;

	CCString C_ccy;
    long ccyId;

    double C_resetGap;
    double C_resetGap_default = 10000.0;

	double C_spread1;
	double C_spread2;
	VECTOR<double> C_spread;

   	double C_fixing1_double;
   	VECTOR<double> C_fixing1_double_vect;
	CCString C_fixing1_str;
	long fixing1_type;
	VECTOR<double> C_fixing1_default;

   	double C_fixing2_double;
   	VECTOR<double> C_fixing2_double_vect;
	CCString C_fixing2_str;
	long fixing2_type;
	VECTOR<double> C_fixing2_default;

	CCString C_intRule;
    long intRuleId;

	CCString C_stubRule;
    long stubRuleId;

	// error
	static int error;
	static char* reason = "";

	XL_readNumCell(XL_startDate,C_startDate," ARM_ERR: start date: date expected",C_result);
	XL_readNumCell(XL_endDate,C_endDate," ARM_ERR: end date: date expected",C_result);
	XL_readStrCell(XL_capOrFloor,C_capOrFloor," ARM_ERR: cap or floor: string expected",C_result);
	XL_readStrOrNumCell(XL_strike,C_strike_str,C_strike_double,strike_type," ARM_ERR: strike: string or numeric expected",C_result);
	
	XL_readNumVector(XL_spread,C_spread," ARM_ERR: spread vector: array of numeric expected", C_result);

	XL_readStrCellWD(XL_payoffLiborType,C_payoffLiborType,"LIBOR3M"," ARM_ERR: payoff libor type  : string expected",C_result);
	
	XL_readStrCellWD(XL_liborType1,C_liborType1,"LIBOR3M"," ARM_ERR: libor type 1 : string expected",C_result);
	XL_readStrCellWD(XL_liborType2,C_liborType2,"LIBOR3M"," ARM_ERR: libor type 2 : string expected",C_result);
	XL_readNumVectorWD(XL_weight,C_weight,C_weight_default," ARM_ERR: weights_slopeflag: array of numeric expected", C_result);
	XL_readStrCellWD(XL_daycount,C_daycount,"30/360"," ARM_ERR: daycount : string expected",C_result);
	XL_readStrCellWD(XL_resetFreq,C_resetFreq,"A"," ARM_ERR: reset frequency : string expected",C_result);
	XL_readStrCellWD(XL_payFreq,C_payFreq,"A"," ARM_ERR: payment frequency : string expected",C_result);
	XL_readStrCellWD(XL_resetTiming,C_resetTiming,"ADV"," ARM_ERR: reset timing : string expected",C_result);
	XL_readStrCellWD(XL_payTiming,C_payTiming,"ARR"," ARM_ERR: payment timing : string expected",C_result);
	XL_readStrCellWD(XL_currency,C_ccy,"DEFAULT"," ARM_ERR: currency : object expected",C_result);
	XL_readNumCellWD(XL_resetGap,C_resetGap,C_resetGap_default," ARM_ERR: resetGap : numeric expected",C_result);
	XL_readStrCellWD(XL_intRule,C_intRule,"DEFAULT"," ARM_ERR: intRule : string expected",C_result);
	XL_readStrCellWD(XL_stubRule,C_stubRule,"DEFAULT"," ARM_ERR: stubRule : string expected",C_result);

	XL_readStrOrNumVectorWD(XL_fixing1,C_fixing1_str,C_fixing1_double_vect,C_fixing1_default,fixing1_type," ARM_ERR: fixing taux : array of numeric expected",C_result);
	XL_readStrOrNumVectorWD(XL_fixing2,C_fixing2_str,C_fixing2_double_vect,C_fixing2_default,fixing2_type," ARM_ERR: fixing taux : array of numeric expected",C_result);

	if(C_weight.size()<2)
	{
		C_result.setMsg ("ARM_ERR: Verify the Vector of Weight");
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}
	
	C_weight1 = C_weight[0];
	C_weight2 = C_weight[1];

	if(C_weight.size() == 2)
	{
		C_slopeFlag = C_slopeFlagDefault;
		C_cptStrikeMethod = C_cptStrikeMethod_default;
		C_computedFormula = C_computedFormula_default;
	}
	else if(C_weight.size() == 3)
	{
		C_slopeFlag = C_weight[2];
		C_cptStrikeMethod = C_cptStrikeMethod_default;
		C_computedFormula = C_computedFormula_default;
	}
	else if(C_weight.size() == 4)
	{
		C_slopeFlag = C_weight[2];
		C_cptStrikeMethod = C_weight[3];
		C_computedFormula = C_computedFormula_default;
	}
	else
	{
		C_slopeFlag = C_weight[2];
		C_cptStrikeMethod = C_weight[3];
		C_computedFormula = C_weight[4];
	}

	if(C_spread.size()<2)
	{
		C_result.setMsg ("ARM_ERR: Verify the Vector of Spread");
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}
	C_spread1 = C_spread[0];
	C_spread2 = C_spread[1];
	
	if ((XL_fixing1->xltype == xltypeMissing) || (XL_fixing1->xltype == xltypeNil))
	{
		fixing1_type = 1;
		C_fixing1_double = ARM_NULL_OBJECT;
    }
	else if ( fixing1_type == XL_TYPE_STRING )
	{
		C_fixing1_double = (double) LocalGetNumObjectId(C_fixing1_str);
		fixing1_type = 1;
	}
	else
	{
	    fixing1_type = 0;
	}

	if ((XL_fixing2->xltype == xltypeMissing) || (XL_fixing2->xltype == xltypeNil))
	{
		fixing2_type = 1;
		C_fixing2_double = ARM_NULL_OBJECT;
    }
	else if ( fixing2_type == XL_TYPE_STRING )
	{
		C_fixing2_double = (double) LocalGetNumObjectId(C_fixing2_str);
		fixing2_type = 1;
	}
	else
	{
	    fixing2_type = 0;
	}

	if((capOrFloorId = ARM_ConvCapOrFloor (C_capOrFloor, C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	payoffLiborTypeId = ARM_ConvIrType (C_payoffLiborType);

	liborType1Id = ARM_ConvIrType (C_liborType1);

	liborType2Id = ARM_ConvIrType (C_liborType2);

	dayCountId = ARM_ConvDayCount (C_daycount);
	
	if((resetFreqId = ARM_ConvFrequency (C_resetFreq, C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	if((payFreqId = ARM_ConvFrequency (C_payFreq, C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	resetTimingId = ARM_ConvPayResetRule (C_resetTiming);

	payTimingId = ARM_ConvPayResetRule (C_payTiming);

	if(C_ccy == "DEFAULT")
	{
		ccyId = ARM_NULL_OBJECT;
	}
	else
	{
		ccyId = LocalGetNumObjectId (C_ccy);
	}

	if ( strike_type == XL_TYPE_STRING )
	{
	   C_strike_double = (double) LocalGetNumObjectId(C_strike_str);

	   strike_type = 1;
	}
	else
	{
	   strike_type = 0;
	}
	
	if(C_intRule == "DEFAULT")
	{
		C_intRule = "ADJ";
	}
	intRuleId = ARM_ConvIntRule(C_intRule);

	if(C_stubRule == "DEFAULT")
	{
		C_stubRule = "SS";
	}
	stubRuleId = ARM_ConvStubRule(C_stubRule);

	long retCode;
	long objId;
	CCString prevClass;
	
	CCString curClass = LOCAL_SPREAD_OPTION_CLASS;
	CCString stringId;

	retCode = ARMLOCAL_SPREADDIGITALFLT( C_startDate,
										 C_endDate,
										 capOrFloorId,
										 (long) strike_type,
										 C_strike_double,
										 (long)payoffLiborTypeId,
										 (long)liborType1Id,
										 (long)liborType2Id,
										 C_weight1,
										 C_weight2,
										 dayCountId,
										 resetFreqId,
										 payFreqId,
										 resetTimingId,
										 payTimingId,
										 ccyId,
										 long(C_resetGap),
										 C_spread1,
										 C_spread2,
										 (long) fixing1_type,
										 C_fixing1_double,
										 C_fixing1_double_vect,
										 (long) fixing2_type,
										 C_fixing2_double,
										 C_fixing2_double_vect,
										 intRuleId,
										 stubRuleId,
										 C_slopeFlag,
									     (int) C_cptStrikeMethod,
										 (int) C_computedFormula,
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_PXL_ARM_SPREADDIGITALFLT" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


__declspec(dllexport) LPXLOPER WINAPI Local_ARM_SPREADDIGITALFLTWithLegs(LPXLOPER XL_firstLeg,
																		 LPXLOPER XL_secondLeg,
																		 LPXLOPER XL_capOrFloor,
																		 LPXLOPER XL_strike,
																		 LPXLOPER XL_spread1,
																		 LPXLOPER XL_spread2,
																		 LPXLOPER XL_payoffLiborType,
																		 LPXLOPER XL_weightVec,
																		 LPXLOPER XL_fixing1,
																		 LPXLOPER XL_fixing2,
																		 LPXLOPER XL_slopeFlag,
																		 LPXLOPER XL_cptStrikeMethod)

{
	ADD_LOG("Local_ARM_SPREADDIGITALFLTWithLegs");
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();

	// C variable
	CCString C_firstLeg;
	long firstLegId;

	CCString C_secondLeg;
	long secondLegId;

	CCString C_capOrFloor;
	long capOrFloorId;

	double C_strike_double;
	CCString C_strike_str;
	long strike_type;

	CCString C_payoffLiborType;
	long payoffLiborTypeId;

	double C_weight1, C_weight2;
	VECTOR<double> C_weight;
	VECTOR<double> C_weight_default(1.0);

	double C_slopeFlag;
	double C_slopeFlagDefault = 1.0;

	VECTOR<double> C_VCtcptStrikeMethod;
	VECTOR<double> C_VCtcptStrikeMethod_default;

	double C_cptStrikeMethod;
	double C_cptStrikeMethod_default = 1.0;

	double C_computedFormula;
	double C_computedFormula_default = 1.0;

	double C_spread1;
	double C_spread2;

   	double C_fixing1_double;
   	VECTOR<double> C_fixing1_double_vect;
	CCString C_fixing1_str;
	long fixing1_type;
	VECTOR<double> C_fixing1_default;

   	double C_fixing2_double;
   	VECTOR<double> C_fixing2_double_vect;
	CCString C_fixing2_str;
	long fixing2_type;
	VECTOR<double> C_fixing2_default;

	// error
	static int error;
	static char* reason = "";
	
	XL_readStrCell(XL_firstLeg,C_firstLeg," ARM_ERR: first leg id: object expected",C_result);
	XL_readStrCell(XL_secondLeg,C_secondLeg," ARM_ERR: secondleg id: object expected",C_result);

	XL_readStrCell(XL_capOrFloor,C_capOrFloor," ARM_ERR: cap or floor: string expected",C_result);
	XL_readStrOrNumCell(XL_strike,C_strike_str,C_strike_double,strike_type," ARM_ERR: strike: string or numeric expected",C_result);
	
	XL_readNumCell(XL_spread1,C_spread1," ARM_ERR: spread1 : numeric expected",C_result);
	XL_readNumCell(XL_spread2,C_spread2," ARM_ERR: spread2 : numeric expected",C_result);

	XL_readStrCellWD(XL_payoffLiborType,C_payoffLiborType,"LIBOR3M"," ARM_ERR: payoff libor type  : string expected",C_result);

	XL_readNumVectorWD(XL_weightVec,C_weight,C_weight_default," ARM_ERR: weights_slopeflag: array of numeric expected", C_result);
	
	XL_readStrOrNumVectorWD(XL_fixing1,C_fixing1_str,C_fixing1_double_vect,C_fixing1_default,fixing1_type," ARM_ERR: fixing taux : array of numeric expected",C_result);
	XL_readStrOrNumVectorWD(XL_fixing2,C_fixing2_str,C_fixing2_double_vect,C_fixing2_default,fixing2_type," ARM_ERR: fixing taux : array of numeric expected",C_result);

	XL_readNumCellWD(XL_slopeFlag,C_slopeFlag,C_slopeFlagDefault," ARM_ERR: slope flag : numeric expected",C_result);

	XL_readNumVectorWD(XL_cptStrikeMethod,C_VCtcptStrikeMethod,C_VCtcptStrikeMethod_default," ARM_ERR: cptStrikeMethod : vector of integer(0/1) expected",C_result);

	if(C_VCtcptStrikeMethod.size() == 0)
	{
		C_cptStrikeMethod = C_cptStrikeMethod_default;
		C_computedFormula = C_computedFormula_default;
	}
	else if(C_VCtcptStrikeMethod.size() == 1)
	{
		C_cptStrikeMethod = C_VCtcptStrikeMethod[0];
		C_computedFormula = C_computedFormula_default;
	}
	else
	{
		C_cptStrikeMethod = C_VCtcptStrikeMethod[0];
		C_computedFormula = C_VCtcptStrikeMethod[1];
	}

	if(C_weight.size()<2)
	{
		C_result.setMsg ("ARM_ERR: Verify the Vector of Weight");
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}
	C_weight1 = C_weight[0];
	C_weight2 = C_weight[1];

	firstLegId = LocalGetNumObjectId (C_firstLeg);
	secondLegId = LocalGetNumObjectId (C_secondLeg);
	
	if ((XL_fixing1->xltype == xltypeMissing) || (XL_fixing1->xltype == xltypeNil))
	{
		fixing1_type = 1;
		C_fixing1_double = ARM_NULL_OBJECT;
    }
	else if ( fixing1_type == XL_TYPE_STRING )
	{
		C_fixing1_double = (double) LocalGetNumObjectId(C_fixing1_str);
		fixing1_type = 1;
	}
	else
	{
	    fixing1_type = 0;
	}

	if ((XL_fixing2->xltype == xltypeMissing) || (XL_fixing2->xltype == xltypeNil))
	{
		fixing2_type = 1;
		C_fixing2_double = ARM_NULL_OBJECT;
    }
	else if ( fixing2_type == XL_TYPE_STRING )
	{
		C_fixing2_double = (double) LocalGetNumObjectId(C_fixing2_str);
		fixing2_type = 1;
	}
	else
	{
	    fixing2_type = 0;
	}

	if((capOrFloorId = ARM_ConvCapOrFloor (C_capOrFloor, C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	payoffLiborTypeId = ARM_ConvIrType (C_payoffLiborType);

	if ( strike_type == XL_TYPE_STRING )
	{
	   C_strike_double = (double) LocalGetNumObjectId(C_strike_str);

	   strike_type = 1;
	}
	else
	{
	   strike_type = 0;
	}

	long retCode;
	long objId;
	CCString prevClass;
	
	CCString curClass = LOCAL_SPREAD_OPTION_CLASS;
	CCString stringId = GetLastCurCellEnvValue ();

	if(!stringId)
	{
		retCode = ARMLOCAL_SPREADDIGITALFLTWithLegs( firstLegId,
													 secondLegId,
													 capOrFloorId,
													 (long) strike_type,
													 C_strike_double,
													 (long)payoffLiborTypeId,
													 C_weight1,
													 C_weight2,
													 C_spread1,
													 C_spread2,
													 (long) fixing1_type,
													 C_fixing1_double,
													 C_fixing1_double_vect,
													 (long) fixing2_type,
													 C_fixing2_double,
													 C_fixing2_double_vect,
													 (long) C_slopeFlag,
													 (int) C_cptStrikeMethod,
													 (int) C_computedFormula,
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
			retCode = ARMLOCAL_SPREADDIGITALFLTWithLegs( firstLegId,
														 secondLegId,
														 capOrFloorId,
														 (long) strike_type,
														 C_strike_double,
														 (long)payoffLiborTypeId,
														 C_weight1,
														 C_weight2,
														 C_spread1,
														 C_spread2,
														 (long) fixing1_type,
														 C_fixing1_double,
														 C_fixing1_double_vect,
														 (long) fixing2_type,
														 C_fixing2_double,
														 C_fixing2_double_vect,
														 (long) C_slopeFlag,
														 (int) C_cptStrikeMethod,
														 (int) C_computedFormula,
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
			retCode = ARMLOCAL_SPREADDIGITALFLTWithLegs( firstLegId,
														 secondLegId,
														 capOrFloorId,
														 (long) strike_type,
														 C_strike_double,
														 (long)payoffLiborTypeId,
														 C_weight1,
														 C_weight2,
														 C_spread1,
														 C_spread2,
														 (long) fixing1_type,
														 C_fixing1_double,
														 C_fixing1_double_vect,
														 (long) fixing2_type,
														 C_fixing2_double,
														 C_fixing2_double_vect,
														 (long) C_slopeFlag,
														 (int) C_cptStrikeMethod,
														 (int) C_computedFormula,
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_ARM_SPREADDIGITALFLTWithLegs" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


__declspec(dllexport) LPXLOPER WINAPI Local_PXL_ARM_SPREADDIGITALFLTWithLegs(LPXLOPER XL_firstLeg,
																			  LPXLOPER XL_secondLeg,
																			  LPXLOPER XL_capOrFloor,
																			  LPXLOPER XL_strike,
																			  LPXLOPER XL_spread1,
																			  LPXLOPER XL_spread2,
																			  LPXLOPER XL_payoffLiborType,
																			  LPXLOPER XL_weightVec,
																			  LPXLOPER XL_fixing1,
																			  LPXLOPER XL_fixing2,
																			  LPXLOPER XL_slopeFlag,
																			  LPXLOPER XL_cptStrikeMethod)
{
	ADD_LOG("Local_PXL_ARM_SPREADDIGITALFLTWithLegs");
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();

	// C variable
	CCString C_firstLeg;
	long firstLegId;

	CCString C_secondLeg;
	long secondLegId;

	CCString C_capOrFloor;
	long capOrFloorId;

	double C_strike_double;
	CCString C_strike_str;
	long strike_type;
	
	CCString C_payoffLiborType;
	long payoffLiborTypeId;

	double C_weight1, C_weight2;
	VECTOR<double> C_weight;
	VECTOR<double> C_weight_default(1.0);

	double C_slopeFlag;
	double C_slopeFlagDefault = 1.0;

	VECTOR<double> C_VCtcptStrikeMethod;
	VECTOR<double> C_VCtcptStrikeMethod_default;

	double C_cptStrikeMethod;
	double C_cptStrikeMethod_default = 1.0;

	double C_computedFormula;
	double C_computedFormula_default = 1.0;

	double C_spread1;
	double C_spread2;

   	double C_fixing1_double;
   	VECTOR<double> C_fixing1_double_vect;
	CCString C_fixing1_str;
	long fixing1_type;
	VECTOR<double> C_fixing1_default;

   	double C_fixing2_double;
   	VECTOR<double> C_fixing2_double_vect;
	CCString C_fixing2_str;
	long fixing2_type;
	VECTOR<double> C_fixing2_default;

	// error
	static int error;
	static char* reason = "";
	
	XL_readStrCell(XL_firstLeg,C_firstLeg," ARM_ERR: first leg id: object expected",C_result);
	XL_readStrCell(XL_secondLeg,C_secondLeg," ARM_ERR: secondleg id: object expected",C_result);

	XL_readStrCell(XL_capOrFloor,C_capOrFloor," ARM_ERR: cap or floor: string expected",C_result);
	XL_readStrOrNumCell(XL_strike,C_strike_str,C_strike_double,strike_type," ARM_ERR: strike: string or numeric expected",C_result);
	
	XL_readNumCell(XL_spread1,C_spread1," ARM_ERR: spread1 : numeric expected",C_result);
	XL_readNumCell(XL_spread2,C_spread2," ARM_ERR: spread2 : numeric expected",C_result);

	XL_readStrCellWD(XL_payoffLiborType,C_payoffLiborType,"LIBOR3M"," ARM_ERR: payoff libor type  : string expected",C_result);

	XL_readNumVectorWD(XL_weightVec,C_weight,C_weight_default," ARM_ERR: weights_slopeflag: array of numeric expected", C_result);
	
	XL_readStrOrNumVectorWD(XL_fixing1,C_fixing1_str,C_fixing1_double_vect,C_fixing1_default,fixing1_type," ARM_ERR: fixing taux : array of numeric expected",C_result);
	XL_readStrOrNumVectorWD(XL_fixing2,C_fixing2_str,C_fixing2_double_vect,C_fixing2_default,fixing2_type," ARM_ERR: fixing taux : array of numeric expected",C_result);

	XL_readNumCellWD(XL_slopeFlag,C_slopeFlag,C_slopeFlagDefault," ARM_ERR: slope flag : numeric expected",C_result);

	XL_readNumVectorWD(XL_cptStrikeMethod,C_VCtcptStrikeMethod,C_VCtcptStrikeMethod_default," ARM_ERR: cptStrikeMethod : vector of integer(0/1) expected",C_result);

	if(C_VCtcptStrikeMethod.size() == 0)
	{
		C_cptStrikeMethod = C_cptStrikeMethod_default;
		C_computedFormula = C_computedFormula_default;
	}
	else if(C_VCtcptStrikeMethod.size() == 1)
	{
		C_cptStrikeMethod = C_VCtcptStrikeMethod[0];
		C_computedFormula = C_computedFormula_default;
	}
	else
	{
		C_cptStrikeMethod = C_VCtcptStrikeMethod[0];
		C_computedFormula = C_VCtcptStrikeMethod[1];
	}

	if(C_weight.size()<2)
	{
		C_result.setMsg ("ARM_ERR: Verify the Vector of Weight");
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}
	C_weight1 = C_weight[0];
	C_weight2 = C_weight[1];

	firstLegId = LocalGetNumObjectId (C_firstLeg);
	secondLegId = LocalGetNumObjectId (C_secondLeg);
	
	if ((XL_fixing1->xltype == xltypeMissing) || (XL_fixing1->xltype == xltypeNil))
	{
		fixing1_type = 1;
		C_fixing1_double = ARM_NULL_OBJECT;
    }
	else if ( fixing1_type == XL_TYPE_STRING )
	{
		C_fixing1_double = (double) LocalGetNumObjectId(C_fixing1_str);
		fixing1_type = 1;
	}
	else
	{
	    fixing1_type = 0;
	}

	if ((XL_fixing2->xltype == xltypeMissing) || (XL_fixing2->xltype == xltypeNil))
	{
		fixing2_type = 1;
		C_fixing2_double = ARM_NULL_OBJECT;
    }
	else if ( fixing2_type == XL_TYPE_STRING )
	{
		C_fixing2_double = (double) LocalGetNumObjectId(C_fixing2_str);
		fixing2_type = 1;
	}
	else
	{
	    fixing2_type = 0;
	}

	if((capOrFloorId = ARM_ConvCapOrFloor (C_capOrFloor, C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	payoffLiborTypeId = ARM_ConvIrType (C_payoffLiborType);

	if ( strike_type == XL_TYPE_STRING )
	{
	   C_strike_double = (double) LocalGetNumObjectId(C_strike_str);

	   strike_type = 1;
	}
	else
	{
	   strike_type = 0;
	}

	long retCode;
	long objId;
	CCString prevClass;
	
	CCString curClass = LOCAL_SPREAD_OPTION_CLASS;
	CCString stringId;

	retCode = ARMLOCAL_SPREADDIGITALFLTWithLegs( firstLegId,
												 secondLegId,
												 capOrFloorId,
												 (long) strike_type,
												 C_strike_double,
												 (long)payoffLiborTypeId,
												 C_weight1,
												 C_weight2,
												 C_spread1,
												 C_spread2,
												 (long) fixing1_type,
												 C_fixing1_double,
												 C_fixing1_double_vect,
												 (long) fixing2_type,
												 C_fixing2_double,
												 C_fixing2_double_vect,
												 (long) C_slopeFlag,
												 (int) C_cptStrikeMethod,
												 (int) C_computedFormula,
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_PXL_ARM_SPREADDIGITALFLTWithLegs" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}

__declspec(dllexport) LPXLOPER WINAPI Local_ARM_QUANTOSPREADDIGITALFLT(	LPXLOPER XL_startDate,
																		LPXLOPER XL_endDate,
 																		LPXLOPER XL_capOrFloor,
																		LPXLOPER XL_strikes,
																		LPXLOPER XL_liborIdx1,
																		LPXLOPER XL_liborIdx2,
																		LPXLOPER XL_liborIdxP,
																		LPXLOPER XL_leg1Weight,
																		LPXLOPER XL_leg2Weight,
																		LPXLOPER XL_leg1Fixings,
																		LPXLOPER XL_leg2Fixings,
																		LPXLOPER XL_leg1Spread,
																		LPXLOPER XL_leg2Spread,
																		LPXLOPER XL_modelParams,
																		LPXLOPER XL_prodParams,
																		LPXLOPER XL_notional)
{
	ADD_LOG("Local_ARM_QUANTOSPREADDIGITALFLT");
	//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();

	double C_startDate;
	double C_endDate;

	CCString C_capOrFloor;
	long capOrFloorId;

	CCString C_strikes_str;
	double C_strike_double;
	long strike_type;

	CCString C_liborIdx1;
	long liborIdx1Id;
	CCString C_liborIdx2;
	long liborIdx2Id;
	CCString C_liborIdxP;
	long liborIdxPId;
	
	double C_Idx1weight;
	double C_Idx1weight_default = 1.0;

	double C_Idx2weight;
	double C_Idx2weight_default = 1.0;

	CCString C_Idx1fixings_str;
	CCString C_Idx1fixings_str_default("DEFAULT");
	long Idx1fixingsId = ARM_NULL_OBJECT;
	CCString C_Idx2fixings_str;
	CCString C_Idx2fixings_str_default("DEFAULT");
	long Idx2fixingsId = ARM_NULL_OBJECT;
	
	double C_Idx1spread;
	double C_Idx1spread_default = 0.0;

	double C_Idx2spread;
	double C_Idx2spread_default = 0.0;

	VECTOR<double> C_modelParams;
	VECTOR<double> C_modelParams_default(0);
	VECTOR<double> modelParams(3);
	
	modelParams[0] = 1.0;
	modelParams[1] = 1.0;
	modelParams[2] = 1.0;

	double C_slopeFlag;
	double C_cptStrikeMethod;
	double C_computedFormula;

	// Product parameters
	VECTOR<CCString> C_prodParams;
	VECTOR<CCString> C_prodParams_default(0);

	VECTOR<CCString> prodParams(13);

	prodParams[0] = "DEFAULT";
	prodParams[1] = "30/360";
	prodParams[2] = "A";
	prodParams[3] = "A";
	prodParams[4] = "ADV";
	prodParams[5] = "ARR";
	prodParams[6] = "ADJ";
	prodParams[7] = "SS";
	prodParams[8] = "2";
	prodParams[9] = "EUR";
	prodParams[10] = "EUR";
	prodParams[11] = "MF";
	prodParams[12] = "NULL";

	CCString C_ccy;
	long ccyId;
	
	CCString C_daycount;
	long dayCountId;

	CCString C_resetFreq;
	long resetFreqId;

	CCString C_payFreq;
	long payFreqId;

	CCString C_resetTiming;
	long resetTimingId;
	
	CCString C_payTiming;
	long payTimingId;

	CCString C_intRule;
    long intRuleId;
	
	CCString C_stubRule;
    long stubRuleId;

    double C_resetGap;
	
	CCString C_resetCal;
	CCString C_payCal;
	CCString C_fwdRule;
	double fwdRuleId;

	CCString C_refDate;
	
	long notionalId = ARM_NULL_OBJECT;
	CCString C_notional_str;
	CCString C_notional_default("DEFAULT");

	// error
	static int error;
	static char* reason = "";

	XL_readNumCell(XL_startDate,C_startDate," ARM_ERR: start date: date expected",C_result);
	XL_readNumCell(XL_endDate,C_endDate," ARM_ERR: end date: date expected",C_result);
	XL_readStrCell(XL_capOrFloor,C_capOrFloor," ARM_ERR: cap or floor: string expected",C_result);
	XL_readStrOrNumCell(XL_strikes,C_strikes_str,C_strike_double,strike_type," ARM_ERR: strikes vector : double or array of numeric expected",C_result);
	XL_readStrCellWD(XL_liborIdx1,C_liborIdx1,"DEFAULT"," ARM_ERR: leg1 libor index : string expected",C_result);
	XL_readStrCellWD(XL_liborIdx2,C_liborIdx2,"DEFAULT"," ARM_ERR: leg2 libor index : string expected",C_result);
	XL_readStrCellWD(XL_liborIdxP,C_liborIdxP,"DEFAULT"," ARM_ERR: legP libor index : string expected",C_result);
	XL_readNumCellWD(XL_leg1Weight,C_Idx1weight,C_Idx1weight_default," ARM_ERR: leg1 weight : numeric expected",C_result);
	XL_readNumCellWD(XL_leg2Weight,C_Idx2weight,C_Idx2weight_default," ARM_ERR: leg2 weight : numeric expected",C_result);
	XL_readStrCellWD(XL_leg1Fixings,C_Idx1fixings_str,C_Idx1fixings_str_default," ARM_ERR: leg1 fixings vector : object expected",C_result);
	XL_readStrCellWD(XL_leg2Fixings,C_Idx2fixings_str,C_Idx2fixings_str_default," ARM_ERR: leg2 fixings vector : object expected",C_result);
	XL_readNumCellWD(XL_leg1Spread,C_Idx1spread, C_Idx1spread_default, " ARM_ERR: leg1 spread : numeric expected",C_result);
	XL_readNumCellWD(XL_leg2Spread,C_Idx2spread,C_Idx2spread_default, " ARM_ERR: leg2 spread : numeric expected",C_result);
	XL_readNumVectorWD(XL_modelParams,C_modelParams,C_modelParams_default," ARM_ERR: model parameters vector : array of numeric expected",C_result);
	XL_readStrVectorWD(XL_prodParams,C_prodParams,C_prodParams_default," ARM_ERR: product parameters vector : array of numeric expected",DOUBLE_TYPE,C_result);
	XL_readStrCellWD(XL_notional,C_notional_str,C_notional_default," ARM_ERR: notional vector : object expected",C_result);

	if((capOrFloorId = ARM_ConvCapOrFloor (C_capOrFloor, C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}
		
	if ( strike_type == XL_TYPE_STRING )
	{
	   C_strike_double = (double) LocalGetNumObjectId(C_strikes_str);

	   strike_type = 1;
	}
	else
	{
	   strike_type = 0;
	}

	liborIdx1Id = strcmp(C_liborIdx1,"DEFAULT") ? (double) LocalGetNumObjectId(C_liborIdx1) : -1.0;	
	liborIdx2Id = strcmp(C_liborIdx2,"DEFAULT") ? (double) LocalGetNumObjectId(C_liborIdx2) : -1.0;
	liborIdxPId = strcmp(C_liborIdxP,"DEFAULT") ? (double) LocalGetNumObjectId(C_liborIdxP) : -1.0;	
	
	Idx1fixingsId = strcmp(C_Idx1fixings_str,"DEFAULT") ? (long) LocalGetNumObjectId(C_Idx1fixings_str) : -1.0;
	Idx2fixingsId = strcmp(C_Idx2fixings_str,"DEFAULT") ? (long) LocalGetNumObjectId(C_Idx2fixings_str) : -1.0;
	
	int i = 0;
	for (i; i < C_modelParams.size(); i++)
	{
		modelParams[i]  = C_modelParams[i];
	}
	
	C_slopeFlag = modelParams[0];
	C_cptStrikeMethod = modelParams[1];
	C_computedFormula = modelParams[2];

	for (i = 0; i<C_prodParams.size(); i++)
	{
		if (strcmp(C_prodParams[i],"DEFAULT"))
			prodParams[i] = C_prodParams[i];
	}

	C_ccy = prodParams[0];
	
	if(C_ccy == "DEFAULT")
	{
		ccyId = ARM_NULL_OBJECT;
	}
	else
	{
		ccyId = LocalGetNumObjectId (C_ccy);
	}
	
	C_daycount = prodParams[1];
	dayCountId = ARM_ConvDayCount (C_daycount);
	
	C_resetFreq = prodParams[2];
	if((resetFreqId = ARM_ConvFrequency (C_resetFreq, C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}
	
	C_payFreq = prodParams[3];
	if((payFreqId = ARM_ConvFrequency (C_payFreq, C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	C_resetTiming = prodParams[4];
	resetTimingId = ARM_ConvPayResetRule (C_resetTiming);
	
	C_payTiming = prodParams[5];
	payTimingId = ARM_ConvPayResetRule (C_payTiming);

	C_intRule  = prodParams[6];
    intRuleId = ARM_ConvIntRule(C_intRule);
	
	C_stubRule = prodParams[7];
	stubRuleId = ARM_ConvStubRule(C_stubRule);

    C_resetGap = atof(prodParams[8]);
	
	C_resetCal = prodParams[9];

	C_payCal = prodParams[10];

	C_fwdRule = prodParams[11];
	fwdRuleId = ARM_ConvFwdRule(C_fwdRule);

	C_refDate = prodParams[12];

	notionalId = strcmp(C_notional_str,"DEFAULT") ? (long) LocalGetNumObjectId(C_notional_str) : -1.0;

	long retCode;
	long objId;
	
	CCString prevClass;
	CCString curClass = LOCAL_SPREAD_OPTION_CLASS;
	CCString stringId = GetLastCurCellEnvValue();
	
	if(!stringId)
	{
		retCode = ARMLOCAL_QUANTOSPREADDIGITALFLT (	C_startDate,
													C_endDate,
													capOrFloorId,
													(long) strike_type,
													C_strike_double,
													(long) liborIdx1Id,
													(long) liborIdx2Id,
													(long) liborIdxPId,
													C_Idx1weight,
													C_Idx2weight,
													Idx1fixingsId,
													Idx2fixingsId,
													C_Idx1spread,
													C_Idx2spread,
													C_slopeFlag,
													C_cptStrikeMethod,
													(int) C_computedFormula,
													ccyId,
													dayCountId,
													resetFreqId,
													payFreqId,
													resetTimingId,
													payTimingId,
													intRuleId,
													stubRuleId,
													C_resetGap,
													C_resetCal,
													C_payCal,
													fwdRuleId,
													C_refDate,
													notionalId,
													C_result);
		if (retCode == ARM_OK)
		{
			objId = C_result.getLong();

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
			retCode = ARMLOCAL_QUANTOSPREADDIGITALFLT( C_startDate,
													C_endDate,
													capOrFloorId,
													(long) strike_type,
													C_strike_double,
													(long) liborIdx1Id,
													(long) liborIdx2Id,
													(long) liborIdxPId,
													C_Idx1weight,
													C_Idx2weight,
													Idx1fixingsId,
													Idx2fixingsId,
													C_Idx1spread,
													C_Idx2spread,
													C_slopeFlag,
													C_cptStrikeMethod,
													(int) C_computedFormula,
													ccyId,
													dayCountId,
													resetFreqId,
													payFreqId,
													resetTimingId,
													payTimingId,
													intRuleId,
													stubRuleId,
													C_resetGap,
													C_resetCal,
													C_payCal,
													fwdRuleId,
													C_refDate,
													notionalId,
													C_result,
													objId);
			if (retCode == ARM_OK)
			{
				objId = C_result.getLong();

				stringId = LocalMakeObjectId (objId, curClass);
			}
		}
		else
		{
			FreeCurCellContent ();
			retCode = ARMLOCAL_QUANTOSPREADDIGITALFLT (C_startDate,
													C_endDate,
													capOrFloorId,
													(long) strike_type,
													C_strike_double,
													(long) liborIdx1Id,
													(long) liborIdx2Id,
													(long) liborIdxPId,
													C_Idx1weight,
													C_Idx2weight,
													Idx1fixingsId,
													Idx2fixingsId,
													C_Idx1spread,
													C_Idx2spread,
													C_slopeFlag,
													C_cptStrikeMethod,
													(int) C_computedFormula,
													ccyId,
													dayCountId,
													resetFreqId,
													payFreqId,
													resetTimingId,
													payTimingId,
													intRuleId,
													stubRuleId,
													C_resetGap,
													C_resetCal,
													C_payCal,
													fwdRuleId,
													C_refDate,
													notionalId,
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_ARM_QUANTOSPREADDIGITAL  "  )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}

__declspec(dllexport) LPXLOPER WINAPI Local_CapLetPrice (LPXLOPER XL_secId,
														 LPXLOPER XL_modId,
														 LPXLOPER XL_numEx)
{
	ADD_LOG("Local_CapLetPrice ");
	//ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();

	// C variable
	CCString C_secId;
	CCString C_modId;
	double C_numEx;
	
	// error
	static int error;
	static char* reason = "";

	XL_readStrCell(XL_secId,C_secId," ARM_ERR: security id: object expected",C_result);
	XL_readStrCell(XL_modId,C_modId," ARM_ERR: model id: object expected",C_result);
	XL_readNumCell(XL_numEx,C_numEx," ARM_ERR: number of exercises: numeric expected",C_result);
	
	long retCode = ARMLOCAL_CapLetPrice (LocalGetNumObjectId (C_secId), LocalGetNumObjectId (C_modId),
									(long)C_numEx, C_result);
		                                   
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

	//ARM_END();
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_CapLetPrice" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


__declspec(dllexport) LPXLOPER WINAPI Local_RATCHET (LPXLOPER XL_swapLeg,
													 LPXLOPER XL_isItCapOrFloor,
													 LPXLOPER XL_strike,
													 LPXLOPER XL_spreadDates,
													 LPXLOPER XL_spreadValues,
													 LPXLOPER XL_correlDates,
													 LPXLOPER XL_correlValues,
													 LPXLOPER XL_fwdVolsDates,
													 LPXLOPER XL_fwdVolsValues)
{
	ADD_LOG("Local_RATCHET ");
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

    // to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();

	// C variable
	CCString C_swapLeg;
	
	CCString C_isItCapOrFloor;
	long isItCapOrFloorId;

	double C_strike;

	VECTOR<double> C_spreadDates;
	VECTOR<double> C_spreadValues;

	VECTOR<double> C_correlDates;
	VECTOR<double> C_correlValues;

	VECTOR<double> C_fwdVolsDates;
	VECTOR<double> C_fwdVolsValues;

	// error
	static int error;
	static char* reason = "";

	XL_readStrCell(XL_swapLeg,C_swapLeg," ARM_ERR: swap leg id: object expected",C_result);
	XL_readStrCell(XL_isItCapOrFloor,C_isItCapOrFloor," ARM_ERR: cap or floor: string expected",C_result);
	XL_readNumCell(XL_strike,C_strike," ARM_ERR: strike: numeric expected",C_result);
	XL_readNumVector(XL_spreadDates,C_spreadDates," ARM_ERR: spread dates: array of date expected",C_result);
	XL_readNumVector(XL_spreadValues,C_spreadValues," ARM_ERR: spread values: array of numeric expected",C_result);
	XL_readNumVector(XL_correlDates,C_correlDates," ARM_ERR: correlation dates: array of date expected",C_result);
	XL_readNumVector(XL_correlValues,C_correlValues," ARM_ERR: correlation values: array of numeric expected",C_result);
	XL_readNumVector(XL_fwdVolsDates,C_fwdVolsDates," ARM_ERR: forward volatility dates: array of date expected",C_result);
	XL_readNumVector(XL_fwdVolsValues,C_fwdVolsValues," ARM_ERR: forward volatility values: array of numeric expected",C_result);
		
	if((isItCapOrFloorId = ARM_ConvCapOrFloor (C_isItCapOrFloor, C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	long retCode;
	long objId;
	CCString prevClass;
	
	CCString curClass = LOCAL_RATCHET_CLASS;
	CCString stringId = GetLastCurCellEnvValue ();
	
	if(!stringId)
	{
		retCode = ARMLOCAL_RATCHET (LocalGetNumObjectId (C_swapLeg),
									isItCapOrFloorId,
									C_strike,
									C_spreadDates,
									C_spreadValues,
									C_correlDates,
									C_correlValues,
									C_fwdVolsDates,
									C_fwdVolsValues,
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
			retCode = ARMLOCAL_RATCHET (LocalGetNumObjectId (C_swapLeg),
										isItCapOrFloorId,
										C_strike,
										C_spreadDates,
										C_spreadValues,
										C_correlDates,
										C_correlValues,
										C_fwdVolsDates,
										C_fwdVolsValues,
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
			retCode = ARMLOCAL_RATCHET (LocalGetNumObjectId (C_swapLeg),
										isItCapOrFloorId,
										C_strike,
										C_spreadDates,
										C_spreadValues,
										C_correlDates,
										C_correlValues,
										C_fwdVolsDates,
										C_fwdVolsValues,
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_RATCHET" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}



__declspec(dllexport) LPXLOPER WINAPI Local_PXL_RATCHET (LPXLOPER XL_swapLeg,
														 LPXLOPER XL_isItCapOrFloor,
														 LPXLOPER XL_strike,
														 LPXLOPER XL_spreadDates,
														 LPXLOPER XL_spreadValues,
														 LPXLOPER XL_correlDates,
														 LPXLOPER XL_correlValues,
														 LPXLOPER XL_fwdVolsDates,
														 LPXLOPER XL_fwdVolsValues)
{
	ADD_LOG("Local_PXL_RATCHET ");
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();

	// C variable
	CCString C_swapLeg;
	
	CCString C_isItCapOrFloor;
	long isItCapOrFloorId;

	double C_strike;

	VECTOR<double> C_spreadDates;
	VECTOR<double> C_spreadValues;

	VECTOR<double> C_correlDates;
	VECTOR<double> C_correlValues;

	VECTOR<double> C_fwdVolsDates;
	VECTOR<double> C_fwdVolsValues;

	// error
	static int error;
	static char* reason = "";

	XL_readStrCell(XL_swapLeg,C_swapLeg," ARM_ERR: swap leg id: object expected",C_result);
	XL_readStrCell(XL_isItCapOrFloor,C_isItCapOrFloor," ARM_ERR: cap or floor: string expected",C_result);
	XL_readNumCell(XL_strike,C_strike," ARM_ERR: strike: numeric expected",C_result);
	XL_readNumVector(XL_spreadDates,C_spreadDates," ARM_ERR: spread dates: array of date expected",C_result);
	XL_readNumVector(XL_spreadValues,C_spreadValues," ARM_ERR: spread values: array of numeric expected",C_result);
	XL_readNumVector(XL_correlDates,C_correlDates," ARM_ERR: correlation dates: array of date expected",C_result);
	XL_readNumVector(XL_correlValues,C_correlValues," ARM_ERR: correlation values: array of numeric expected",C_result);
	XL_readNumVector(XL_fwdVolsDates,C_fwdVolsDates," ARM_ERR: forward volatility dates: array of date expected",C_result);
	XL_readNumVector(XL_fwdVolsValues,C_fwdVolsValues," ARM_ERR: forward volatility values: array of numeric expected",C_result);
		
	if((isItCapOrFloorId = ARM_ConvCapOrFloor (C_isItCapOrFloor, C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	long retCode;
	long objId;
	
	CCString curClass = LOCAL_RATCHET_CLASS;
	CCString stringId;
	
	retCode = ARMLOCAL_RATCHET (LocalGetNumObjectId (C_swapLeg),
								isItCapOrFloorId,
								C_strike,
								C_spreadDates,
								C_spreadValues,
								C_correlDates,
								C_correlValues,
								C_fwdVolsDates,
								C_fwdVolsValues,
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_PXL_RATCHET" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


__declspec(dllexport) LPXLOPER WINAPI Local_ARM_DIGITAL (LPXLOPER XL_swapLeg,
														 LPXLOPER XL_isItCapOrFloor,
														 LPXLOPER XL_strike,
														 LPXLOPER XL_spread1,
														 LPXLOPER XL_spread2,
														 LPXLOPER XL_payoff)
{
	ADD_LOG("Local_ARM_DIGITAL ");
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

    // to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();

	// C variable
	CCString C_swapLeg;
	
	CCString C_isItCapOrFloor;
	long isItCapOrFloorId;

	double   C_strike_double;
	CCString C_strike_str;
	long     strikeType;

	double C_spread1;
	double C_spread2;

	double   C_payoff_double;
	CCString C_payoff_str;
	long     payoffType;

	// error
	static int error;
	static char* reason = "";

	XL_readStrCell(XL_swapLeg,C_swapLeg," ARM_ERR: swap leg id: object expected",C_result);
	XL_readStrCell(XL_isItCapOrFloor,C_isItCapOrFloor," ARM_ERR: cap or floor: string expected",C_result);
	XL_readStrOrNumCell(XL_strike, C_strike_str, C_strike_double, strikeType,
		   " ARM_ERR: strike: numeric or object ID string expected",C_result);
	XL_readNumCell(XL_spread1, C_spread1," ARM_ERR: spread1 : numeric expected",C_result);
	XL_readNumCell(XL_spread2, C_spread2," ARM_ERR: spread2 : numeric expected",C_result);
	XL_readStrOrNumCell(XL_payoff, C_payoff_str,C_payoff_double, payoffType,
		   " ARM_ERR: payoff: numeric or object ID string expected",C_result);

	if((isItCapOrFloorId = ARM_ConvCapOrFloor (C_isItCapOrFloor, C_result)) == ARM_DEFAULT_ERR)
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

	if ( payoffType == XL_TYPE_STRING )
	{
		C_payoff_double = (double) LocalGetNumObjectId(C_payoff_str);

		payoffType = 1L;
	}
	else
	{
		payoffType = 0L;
	}

	long retCode;
	long objId;
	CCString prevClass;

	CCString curClass = LOCAL_DIGITAL_CLASS;
	CCString stringId = GetLastCurCellEnvValue ();

	if(!stringId)
	{
		retCode = ARMLOCAL_DIGITAL (LocalGetNumObjectId (C_swapLeg),
									isItCapOrFloorId,
									strikeType,
									C_strike_double,
									C_spread1,
									C_spread2,
									payoffType,
									C_payoff_double,
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
			retCode = ARMLOCAL_DIGITAL (LocalGetNumObjectId (C_swapLeg),
										isItCapOrFloorId,
										strikeType,
										C_strike_double,
										C_spread1,
										C_spread2,
										payoffType,
										C_payoff_double,
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
			retCode = ARMLOCAL_DIGITAL (LocalGetNumObjectId (C_swapLeg),
										isItCapOrFloorId,
										strikeType,
										C_strike_double,
										C_spread1,
										C_spread2,
										payoffType,
										C_payoff_double,
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_ARM_DIGITAL" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}

__declspec(dllexport) LPXLOPER WINAPI Local_PXL_ARM_DIGITAL (LPXLOPER XL_swapLeg,
															 LPXLOPER XL_isItCapOrFloor,
															 LPXLOPER XL_strike,
															 LPXLOPER XL_spread1,
															 LPXLOPER XL_spread2,
															 LPXLOPER XL_payoff)
{
	ADD_LOG("Local_PXL_ARM_DIGITAL ");
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

    // to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();

	// C variable
	CCString C_swapLeg;
	
	CCString C_isItCapOrFloor;
	long isItCapOrFloorId;

	double   C_strike_double;
	CCString C_strike_str;
	long     strikeType;

	double C_spread1;
	double C_spread2;

	double   C_payoff_double;
	CCString C_payoff_str;
	long     payoffType;

	// error
	static int error;
	static char* reason = "";

	XL_readStrCell(XL_swapLeg,C_swapLeg," ARM_ERR: swap leg id: object expected",C_result);
	XL_readStrCell(XL_isItCapOrFloor,C_isItCapOrFloor," ARM_ERR: cap or floor: string expected",C_result);
	XL_readStrOrNumCell(XL_strike, C_strike_str, C_strike_double, strikeType,
		   " ARM_ERR: strike: numeric or object ID string expected",C_result);
	XL_readNumCell(XL_spread1, C_spread1," ARM_ERR: spread1 : numeric expected",C_result);
	XL_readNumCell(XL_spread2, C_spread2," ARM_ERR: spread2 : numeric expected",C_result);
	XL_readStrOrNumCell(XL_payoff, C_payoff_str,C_payoff_double, payoffType,
		   " ARM_ERR: payoff: numeric or object ID string expected",C_result);

	if((isItCapOrFloorId = ARM_ConvCapOrFloor (C_isItCapOrFloor, C_result)) == ARM_DEFAULT_ERR)
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

	if ( payoffType == XL_TYPE_STRING )
	{
		C_payoff_double = (double) LocalGetNumObjectId(C_payoff_str);

		payoffType = 1L;
	}
	else
	{
		payoffType = 0L;
	}

	long retCode;
	long objId;

	CCString curClass = LOCAL_DIGITAL_CLASS;
	CCString stringId;

	retCode = ARMLOCAL_DIGITAL (LocalGetNumObjectId (C_swapLeg),
								isItCapOrFloorId,
								strikeType,
								C_strike_double,
								C_spread1,
								C_spread2,
								payoffType,
								C_payoff_double,
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_PXL_ARM_DIGITAL" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}



__declspec(dllexport) LPXLOPER WINAPI Local_ARM_DUALCAP (LPXLOPER XL_startDate,
														 LPXLOPER XL_endDate,
														 LPXLOPER XL_isItCapOrFloor,
														 LPXLOPER XL_strike,
														 LPXLOPER XL_liborType1,
														 LPXLOPER XL_liborType2,
														 LPXLOPER XL_daycount,
														 LPXLOPER XL_Freq,
														 LPXLOPER XL_resetTiming1,
														 LPXLOPER XL_resetTiming2,
														 LPXLOPER XL_resetGap1,
														 LPXLOPER XL_resetGap2,
														 LPXLOPER XL_indexCcyId1,
														 LPXLOPER XL_indexCcyId2,
														 LPXLOPER XL_discountCcyId,
														 LPXLOPER XL_dresetCal,
														 LPXLOPER XL_fresetCal)
{
	ADD_LOG("Local_ARM_DUALCAP ");
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
	
	CCString C_isItCapOrFloor;
	long isItCapOrFloorId;

	double   C_strike_double;
	CCString C_strike_str;
	long     strikeType;

	CCString C_liborType1;
	long liborType1Id;

	CCString C_liborType2;
	long liborType2Id;

	CCString C_daycount;
	long dayCountId;

	CCString C_Freq;
	long FreqId;

	CCString C_resetTiming1;
	long resetTiming1Id;

	CCString C_resetTiming2;
	long resetTiming2Id;

	double C_resetGap1;
	double C_resetGap1_default = 10000.;

	double C_resetGap2;
	double C_resetGap2_default = 10000.;

	CCString C_indexCcy1Id;
	long ccy1Id;
	CCString C_indexCcy2Id;
	long ccy2Id;
	CCString C_discountCcyId;
	long discountCcyId;

	CCString C_dresetCal;
	
	CCString C_fresetCal;
	

	// error
	static int error;
	static char* reason = "";

	XL_readNumCell(XL_startDate, C_startDate," ARM_ERR: start Date : numeric expected",C_result);
	XL_readNumCell(XL_endDate, C_endDate," ARM_ERR: end Date : numeric expected",C_result);
	XL_readStrCell(XL_isItCapOrFloor,C_isItCapOrFloor," ARM_ERR: cap or floor: string expected",C_result);
	XL_readStrOrNumCell(XL_strike, C_strike_str, C_strike_double, strikeType,
		   " ARM_ERR: strike: numeric or object ID string expected",C_result);
	XL_readStrCellWD(XL_liborType1, C_liborType1, "LIBOR3M"," ARM_ERR: libor Type 1 : string expected",C_result);
	XL_readStrCellWD(XL_liborType2, C_liborType2, "LIBOR3M"," ARM_ERR: libor Type 2 : string expected",C_result);
	XL_readStrCellWD(XL_daycount, C_daycount, "30/360"," ARM_ERR: daycount : string expected",C_result);
	XL_readStrCellWD(XL_Freq, C_Freq, "Q"," ARM_ERR: indexes and payment Frequency : string expected",C_result);
	XL_readStrCellWD(XL_resetTiming1, C_resetTiming1, "ADV"," ARM_ERR: reset Timing 1 : string expected",C_result);
	XL_readStrCellWD(XL_resetTiming2, C_resetTiming2, "ADV"," ARM_ERR: reset Timing 2 : string expected",C_result);
	XL_readNumCellWD(XL_resetGap1,C_resetGap1,C_resetGap1_default," ARM_ERR: reset Gap of index 1: numeric expected",C_result);
	XL_readNumCellWD(XL_resetGap2,C_resetGap2,C_resetGap2_default," ARM_ERR: reset Gap of index 2: numeric expected",C_result);
	XL_readStrCellWD(XL_indexCcyId1, C_indexCcy1Id, "DEFAULT"," ARM_ERR: index currency Id 1 : string expected",C_result);
	XL_readStrCellWD(XL_indexCcyId2, C_indexCcy2Id, "DEFAULT"," ARM_ERR: index currency Id 2 : string expected",C_result);
	XL_readStrCellWD(XL_discountCcyId, C_discountCcyId, "DEFAULT"," ARM_ERR: discount currency Id  : string expected",C_result);
	XL_readStrCellWD(XL_dresetCal,C_dresetCal,"NULL"," ARM_ERR: first index reset Calendar: string expected",C_result);
	XL_readStrCellWD(XL_fresetCal,C_fresetCal,"NULL"," ARM_ERR: second index reset Calendar: string expected",C_result);
	
	if((isItCapOrFloorId = ARM_ConvCapOrFloor (C_isItCapOrFloor, C_result)) == ARM_DEFAULT_ERR)
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

	if((FreqId = ARM_ConvFrequency (C_Freq, C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	resetTiming1Id = ARM_ConvPayResetRule (C_resetTiming1);

	resetTiming2Id = ARM_ConvPayResetRule (C_resetTiming2);

	if((liborType1Id = ARM_ConvIrIndName (C_liborType1, C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	dayCountId = ARM_ConvDayCount (C_daycount);

	if((liborType2Id = ARM_ConvIrIndName (C_liborType2, C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	if(C_indexCcy1Id == "DEFAULT")
	{
		if((liborType1Id == K_PIBOR1M) ||
		   (liborType1Id == K_PIBOR3M) ||
		   (liborType1Id == K_PIBOR6M) ||
		   (liborType1Id == K_PIBOR1Y))
		{
			ccy1Id = ARM_FRF_CCY_OBJECT;
		}
		else
		{
			ccy1Id = ARM_NULL_OBJECT;
		}
	}
	else
	{
		ccy1Id = LocalGetNumObjectId (C_indexCcy1Id);
	}

	if(C_indexCcy2Id == "DEFAULT")
	{
		if((liborType2Id == K_PIBOR1M) ||
		   (liborType2Id == K_PIBOR3M) ||
		   (liborType2Id == K_PIBOR6M) ||
		   (liborType2Id == K_PIBOR1Y))
		{
			ccy2Id = ARM_FRF_CCY_OBJECT;
		}
		else
		{
			ccy2Id = ARM_NULL_OBJECT;
		}
	}
	else
	{
		ccy2Id = LocalGetNumObjectId (C_indexCcy2Id);
	}

	if(C_discountCcyId == "DEFAULT")
	{
		if((liborType1Id == K_PIBOR1M) ||
		   (liborType1Id == K_PIBOR3M) ||
		   (liborType1Id == K_PIBOR6M) ||
		   (liborType1Id == K_PIBOR1Y))
		{
			discountCcyId = ARM_FRF_CCY_OBJECT;
		}
		else
		{
			discountCcyId = ARM_NULL_OBJECT;
		}
	}
	else
	{
		discountCcyId = LocalGetNumObjectId (C_discountCcyId);
	}


	long retCode;
	long objId;
	CCString prevClass;

	CCString curClass = LOCAL_DUALCAP_CLASS;
	CCString stringId = GetLastCurCellEnvValue ();

	if(!stringId)
	{
		retCode = ARMLOCAL_DUALCAP (C_startDate,
									C_endDate,
									isItCapOrFloorId,
									strikeType,
									C_strike_double,
									liborType1Id,
									liborType2Id,
									dayCountId,
									FreqId,
									resetTiming1Id,
									resetTiming2Id,
									C_resetGap1,
									C_resetGap2,	
									ccy1Id,
									ccy2Id,
									discountCcyId,
									C_dresetCal,
									C_fresetCal,
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
			retCode = ARMLOCAL_DUALCAP (C_startDate,
										C_endDate,
										isItCapOrFloorId,
										strikeType,
										C_strike_double,
										liborType1Id,
										liborType2Id,
										dayCountId,
										FreqId,
										resetTiming1Id,
										resetTiming2Id,
										C_resetGap1,
										C_resetGap2,
										ccy1Id,
										ccy2Id,
										discountCcyId,
										C_dresetCal,
										C_fresetCal,
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
			retCode = ARMLOCAL_DUALCAP (C_startDate,
										C_endDate,
										isItCapOrFloorId,
										strikeType,
										C_strike_double,
										liborType1Id,
										liborType2Id,
										dayCountId,
										FreqId,
										resetTiming1Id,
										resetTiming2Id,
										C_resetGap1,
										C_resetGap2,
										ccy1Id,
										ccy2Id,
										discountCcyId,
										C_dresetCal,
										C_fresetCal,
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_ARM_DUALCAP" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


__declspec(dllexport) LPXLOPER WINAPI Local_PXL_ARM_DUALCAP (LPXLOPER XL_startDate,
															 LPXLOPER XL_endDate,
															 LPXLOPER XL_isItCapOrFloor,
															 LPXLOPER XL_strike,
															 LPXLOPER XL_liborType1,
															 LPXLOPER XL_liborType2,
															 LPXLOPER XL_daycount,
															 LPXLOPER XL_Freq,
															 LPXLOPER XL_resetTiming1,
															 LPXLOPER XL_resetTiming2,
															 LPXLOPER XL_resetGap1,
															 LPXLOPER XL_resetGap2,
															 LPXLOPER XL_indexCcyId1,
															 LPXLOPER XL_indexCcyId2,
															 LPXLOPER XL_discountCcyId,
															 LPXLOPER XL_dresetCal,
															 LPXLOPER XL_fresetCal)
{
	ADD_LOG("Local_PXL_ARM_DUALCAP ");
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
	
	CCString C_isItCapOrFloor;
	long isItCapOrFloorId;

	double   C_strike_double;
	CCString C_strike_str;
	long     strikeType;

	CCString C_liborType1;
	long liborType1Id;

	CCString C_liborType2;
	long liborType2Id;

	CCString C_daycount;
	long dayCountId;

	CCString C_Freq;
	long FreqId;

	CCString C_resetTiming1;
	long resetTiming1Id;

	CCString C_resetTiming2;
	long resetTiming2Id;

	double C_resetGap1;
	double C_resetGap1_default = 10000.;

	double C_resetGap2;
	double C_resetGap2_default = 10000.;

	CCString C_indexCcy1Id;
	long ccy1Id;
	CCString C_indexCcy2Id;
	long ccy2Id;
	CCString C_discountCcyId;
	long discountCcyId;

	CCString C_dresetCal;

	CCString C_fresetCal;

	// error
	static int error;
	static char* reason = "";

	XL_readNumCell(XL_startDate, C_startDate," ARM_ERR: start Date : numeric expected",C_result);
	XL_readNumCell(XL_endDate, C_endDate," ARM_ERR: end Date : numeric expected",C_result);
	XL_readStrCell(XL_isItCapOrFloor,C_isItCapOrFloor," ARM_ERR: cap or floor: string expected",C_result);
	XL_readStrOrNumCell(XL_strike, C_strike_str, C_strike_double, strikeType,
		   " ARM_ERR: strike: numeric or object ID string expected",C_result);
	XL_readStrCellWD(XL_liborType1, C_liborType1, "LIBOR3M"," ARM_ERR: libor Type 1 : string expected",C_result);
	XL_readStrCellWD(XL_liborType2, C_liborType2, "LIBOR3M"," ARM_ERR: libor Type 2 : string expected",C_result);
	XL_readStrCellWD(XL_daycount, C_daycount, "30/360"," ARM_ERR: daycount : string expected",C_result);
	XL_readStrCellWD(XL_Freq, C_Freq, "Q"," ARM_ERR: indexes and payment Frequency : string expected",C_result);
	XL_readStrCellWD(XL_resetTiming1, C_resetTiming1, "ADV"," ARM_ERR: reset Timing 1 : string expected",C_result);
	XL_readStrCellWD(XL_resetTiming2, C_resetTiming2, "ADV"," ARM_ERR: reset Timing 2 : string expected",C_result);
	XL_readNumCellWD(XL_resetGap1,C_resetGap1,C_resetGap1_default," ARM_ERR: reset Gap of index 1: numeric expected",C_result);
	XL_readNumCellWD(XL_resetGap2,C_resetGap2,C_resetGap2_default," ARM_ERR: reset Gap of index 2: numeric expected",C_result);
	XL_readStrCellWD(XL_indexCcyId1, C_indexCcy1Id, "DEFAULT"," ARM_ERR: index currency Id 1 : string expected",C_result);
	XL_readStrCellWD(XL_indexCcyId2, C_indexCcy2Id, "DEFAULT"," ARM_ERR: index currency Id 2 : string expected",C_result);
	XL_readStrCellWD(XL_discountCcyId, C_discountCcyId, "DEFAULT"," ARM_ERR: discount currency Id  : string expected",C_result);
	XL_readStrCellWD(XL_dresetCal,C_dresetCal,"NULL"," ARM_ERR: first index reset Calendar: string expected",C_result);
	XL_readStrCellWD(XL_fresetCal,C_fresetCal,"NULL"," ARM_ERR: second index reset Calendar: string expected",C_result);
	
	if((isItCapOrFloorId = ARM_ConvCapOrFloor (C_isItCapOrFloor, C_result)) == ARM_DEFAULT_ERR)
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

	if((FreqId = ARM_ConvFrequency (C_Freq, C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	resetTiming1Id = ARM_ConvPayResetRule (C_resetTiming1);

	resetTiming2Id = ARM_ConvPayResetRule (C_resetTiming2);

	dayCountId = ARM_ConvDayCount (C_daycount);

	if((liborType1Id = ARM_ConvIrIndName (C_liborType1, C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	if((liborType2Id = ARM_ConvIrIndName (C_liborType2, C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	if(C_indexCcy1Id == "DEFAULT")
	{
		if((liborType1Id == K_PIBOR1M) ||
		   (liborType1Id == K_PIBOR3M) ||
		   (liborType1Id == K_PIBOR6M) ||
		   (liborType1Id == K_PIBOR1Y))
		{
			ccy1Id = ARM_FRF_CCY_OBJECT;
		}
		else
		{
			ccy1Id = ARM_NULL_OBJECT;
		}
	}
	else
	{
		ccy1Id = LocalGetNumObjectId (C_indexCcy1Id);
	}

	if(C_indexCcy2Id == "DEFAULT")
	{
		if((liborType2Id == K_PIBOR1M) ||
		   (liborType2Id == K_PIBOR3M) ||
		   (liborType2Id == K_PIBOR6M) ||
		   (liborType2Id == K_PIBOR1Y))
		{
			ccy2Id = ARM_FRF_CCY_OBJECT;
		}
		else
		{
			ccy2Id = ARM_NULL_OBJECT;
		}
	}
	else
	{
		ccy2Id = LocalGetNumObjectId (C_indexCcy2Id);
	}

	if(C_discountCcyId == "DEFAULT")
	{
		if((liborType1Id == K_PIBOR1M) ||
		   (liborType1Id == K_PIBOR3M) ||
		   (liborType1Id == K_PIBOR6M) ||
		   (liborType1Id == K_PIBOR1Y))
		{
			discountCcyId = ARM_FRF_CCY_OBJECT;
		}
		else
		{
			discountCcyId = ARM_NULL_OBJECT;
		}
	}
	else
	{
		discountCcyId = LocalGetNumObjectId (C_discountCcyId);
	}

	long retCode;
	long objId;

	CCString curClass = LOCAL_DUALCAP_CLASS;
	CCString stringId;

	retCode = ARMLOCAL_DUALCAP (C_startDate,
								C_endDate,
								isItCapOrFloorId,
								strikeType,
								C_strike_double,
								liborType1Id,
								liborType2Id,
								dayCountId,
								FreqId,
								resetTiming1Id,
								resetTiming2Id,
								C_resetGap1,
								C_resetGap2,		
								ccy1Id,
								ccy2Id,
								discountCcyId,
								C_dresetCal,
								C_fresetCal,
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_PXL_ARM_DUALCAP" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}

__declspec(dllexport) LPXLOPER WINAPI Local_GlobalCap ( LPXLOPER XL_swapLeg,
														LPXLOPER XL_isItCapOrFloor,
														LPXLOPER XL_globalCapValue,
														LPXLOPER XL_globalCapSpreads,
														LPXLOPER XL_globalCapFixedRates,
														LPXLOPER XL_globalCapBarriers,
														LPXLOPER XL_finalRatio,
														LPXLOPER XL_MCNbIter,
														LPXLOPER XL_globalCapPastFixings)
{
	ADD_LOG("Local_GlobalCap ");
	static XLOPER XL_result;
	ARM_result C_result;

	ARM_XL_TRY_BLOCK_BEGIN
	{
		ARM_NOCALCIFWIZ();

		CCString C_swapLeg;
		
		CCString C_isItCapOrFloor;
		long isItCapOrFloorId;

		double C_globalCapValue;
		double C_finalRatio;
		double C_nbIter;

		CCString C_globalCapSpreadsStr;
		long C_globalCapSpreadsId;

		CCString C_globalCapFixedRatesStr;
		long C_globalCapFixedRatesId;

		CCString C_globalCapBarriersStr;
		long C_globalCapBarriersId;

		CCString C_globalCapPastFixingsStr;
		long C_globalCapPastFixingsId;

		// error
		static int error;
		static char* reason = "";

		XL_readStrCell(XL_swapLeg,C_swapLeg," ARM_ERR: swap leg id: object expected",C_result);
		XL_readStrCell(XL_isItCapOrFloor,C_isItCapOrFloor," ARM_ERR: cap or floor: string expected",C_result);
		
		if((isItCapOrFloorId = ARM_ConvCapOrFloor (C_isItCapOrFloor, C_result)) == ARM_DEFAULT_ERR)
		{
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}

		XL_readNumCell(XL_globalCapValue, C_globalCapValue," ARM_ERR: global cap value: number expected",C_result);
		XL_readNumCell(XL_finalRatio, C_finalRatio," ARM_ERR: global cap value: number expected",C_result);
		XL_readNumCell(XL_MCNbIter, C_nbIter," ARM_ERR: global cap value: number expected",C_result);

		XL_readStrCell(XL_globalCapSpreads,C_globalCapSpreadsStr," ARM_ERR: global cap spreads: reference value expected",C_result);
		XL_readStrCell(XL_globalCapFixedRates,C_globalCapFixedRatesStr," ARM_ERR: global cap fixed rates: reference value expected",C_result);
		XL_readStrCell(XL_globalCapBarriers,C_globalCapBarriersStr," ARM_ERR: global cap barriers: reference value expected",C_result);
		XL_readStrCellWD(XL_globalCapPastFixings,C_globalCapPastFixingsStr,"NULL"," ARM_ERR: global cap past fixings: reference value expected",C_result);

		C_globalCapSpreadsId = LocalGetNumObjectId(C_globalCapSpreadsStr);
		C_globalCapFixedRatesId = LocalGetNumObjectId(C_globalCapFixedRatesStr);
		C_globalCapBarriersId = LocalGetNumObjectId(C_globalCapBarriersStr);
		C_globalCapPastFixingsId = (C_globalCapPastFixingsStr == "NULL") ? ARM_NULL_OBJECT : LocalGetNumObjectId(C_globalCapPastFixingsStr);

		long retCode;
		long objId;
		CCString prevClass;
		
		CCString curClass = LOCAL_GLOBALCAP_CLASS;
		CCString stringId = GetLastCurCellEnvValue();
		
		if(!stringId)
		{
			retCode = ARMLOCAL_GlobalCap(LocalGetNumObjectId (C_swapLeg),
										 isItCapOrFloorId,
										 C_globalCapValue,
										 C_globalCapSpreadsId,
										 C_globalCapFixedRatesId,
										 C_globalCapBarriersId,
										 C_finalRatio,
										 C_nbIter,
										 C_globalCapPastFixingsId,
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
				retCode = ARMLOCAL_GlobalCap(LocalGetNumObjectId (C_swapLeg),
											 isItCapOrFloorId,
											 C_globalCapValue,
											 C_globalCapSpreadsId,
											 C_globalCapFixedRatesId,
											 C_globalCapBarriersId,
											 C_finalRatio,
											 C_nbIter,
											 C_globalCapPastFixingsId,
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
				retCode = ARMLOCAL_GlobalCap(LocalGetNumObjectId (C_swapLeg),
											 isItCapOrFloorId,
											 C_globalCapValue,
											 C_globalCapSpreadsId,
											 C_globalCapFixedRatesId,
											 C_globalCapBarriersId,
											 C_finalRatio,
											 C_nbIter,
											 C_globalCapPastFixingsId,
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
	}
	ARM_XL_TRY_BLOCK_END

	ARM_XL_CATCH_ARM_EXPT

	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_GlobalCap" )

	return (LPXLOPER)&XL_result;
}

__declspec(dllexport) LPXLOPER WINAPI Local_PXL_GlobalCap ( LPXLOPER XL_swapLeg,
															LPXLOPER XL_isItCapOrFloor,
															LPXLOPER XL_globalCapValue,
															LPXLOPER XL_globalCapSpreads,
															LPXLOPER XL_globalCapFixedRates,
															LPXLOPER XL_globalCapBarriers,
															LPXLOPER XL_finalRatio,
															LPXLOPER XL_MCNbIter,
															LPXLOPER XL_globalCapPastFixings)
{
	ADD_LOG("Local_PXL_GlobalCap ");
	static XLOPER XL_result;
	ARM_result C_result;

	ARM_XL_TRY_BLOCK_BEGIN
	{
		ARM_NOCALCIFWIZ();

		CCString C_swapLeg;
		
		CCString C_isItCapOrFloor;
		long isItCapOrFloorId;

		double C_globalCapValue;
		double C_finalRatio;
		double C_nbIter;

		CCString C_globalCapSpreadsStr;
		long C_globalCapSpreadsId;

		CCString C_globalCapFixedRatesStr;
		long C_globalCapFixedRatesId;

		CCString C_globalCapBarriersStr;
		long C_globalCapBarriersId;

		CCString C_globalCapPastFixingsStr;
		long C_globalCapPastFixingsId;

		// error
		static int error;
		static char* reason = "";

		XL_readStrCell(XL_swapLeg,C_swapLeg," ARM_ERR: swap leg id: object expected",C_result);
		XL_readStrCell(XL_isItCapOrFloor,C_isItCapOrFloor," ARM_ERR: cap or floor: string expected",C_result);
		
		if((isItCapOrFloorId = ARM_ConvCapOrFloor (C_isItCapOrFloor, C_result)) == ARM_DEFAULT_ERR)
		{
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}

		XL_readNumCell(XL_globalCapValue, C_globalCapValue," ARM_ERR: global cap value: number expected",C_result);
		XL_readNumCell(XL_finalRatio, C_finalRatio," ARM_ERR: global cap value: number expected",C_result);
		XL_readNumCell(XL_MCNbIter, C_nbIter," ARM_ERR: global cap value: number expected",C_result);

		XL_readStrCell(XL_globalCapSpreads,C_globalCapSpreadsStr," ARM_ERR: global cap spreads: reference value expected",C_result);
		XL_readStrCell(XL_globalCapFixedRates,C_globalCapFixedRatesStr," ARM_ERR: global cap fixed rates: reference value expected",C_result);
		XL_readStrCell(XL_globalCapBarriers,C_globalCapBarriersStr," ARM_ERR: global cap barriers: reference value expected",C_result);
		XL_readStrCellWD(XL_globalCapPastFixings,C_globalCapPastFixingsStr,"NULL"," ARM_ERR: global cap past fixings: reference value expected",C_result);

		C_globalCapSpreadsId = LocalGetNumObjectId(C_globalCapSpreadsStr);
		C_globalCapFixedRatesId = LocalGetNumObjectId(C_globalCapFixedRatesStr);
		C_globalCapBarriersId = LocalGetNumObjectId(C_globalCapBarriersStr);
		C_globalCapPastFixingsId = (C_globalCapPastFixingsStr == "NULL") ? ARM_NULL_OBJECT : LocalGetNumObjectId(C_globalCapPastFixingsStr);

		long retCode;
		long objId;
		CCString prevClass;
		
		CCString curClass = LOCAL_GLOBALCAP_CLASS;
		CCString stringId = GetLastCurCellEnvValue();
		
		retCode = ARMLOCAL_GlobalCap(LocalGetNumObjectId (C_swapLeg),
									 isItCapOrFloorId,
									 C_globalCapValue,
									 C_globalCapSpreadsId,
									 C_globalCapFixedRatesId,
									 C_globalCapBarriersId,
									 C_finalRatio,
									 C_nbIter,
									 C_globalCapPastFixingsId,
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
	}
	ARM_XL_TRY_BLOCK_END

	ARM_XL_CATCH_ARM_EXPT

	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_PXL_GlobalCap" )

	return (LPXLOPER)&XL_result;
}
