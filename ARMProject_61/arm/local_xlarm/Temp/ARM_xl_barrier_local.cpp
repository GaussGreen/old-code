#pragma warning(disable : 4786)
#pragma warning(disable : 4005)

#include "ARM_xl_barrier_local.h"
#include "ARM_xl_trycatch_local.h"
#include "ARM_local_interglob.h"
#include "ARM_local_interface.h"
#include <ARM\libarm_local\ARM_local_class.h>
#include <ARM\libarm_local\ARM_local_barrier.h>
#include <util\fromto.h>

#include "util\tech_macro.h"

__declspec(dllexport) LPXLOPER WINAPI Local_CONSTBARRIER (LPXLOPER XL_underlying,
														  LPXLOPER XL_tAsset,
														  LPXLOPER XL_maturity,
														  LPXLOPER XL_barrier1,
														  LPXLOPER XL_upDownDouble,
														  LPXLOPER XL_inOut,
														  LPXLOPER XL_triggerVar,
														  LPXLOPER XL_rebate,
														  LPXLOPER XL_firstX,
														  LPXLOPER XL_exerciseTiming,
														  LPXLOPER XL_barrier2)
{
	ADD_LOG("Local_CONSTBARRIER ");
	//ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{

	ARM_NOCALCIFWIZ();

	// C variable
	CCString C_underlying;
	CCString C_tAsset;

	double C_maturity;

	double C_barrier1;
	double C_barrier2;
	double C_barrier2_default = 0.0;

	CCString C_upDownDouble;
	long upDownDoubleId;

	CCString C_inOut;
	long inOutId;

	double C_triggerVar;
	double C_triggerVar_default = 0.0;

	double C_rebate;
	double C_rebate_default = 0.0;

	double C_firstX;
	double C_firstX_default = -1.0;

	CCString C_exerciseTiming;
	CCString C_exerciseTiming_default = "ADVANCE";
	long exerciseTimingId = K_ADVANCE;

	// error
	static int error;
	static char* reason = "";
	XL_readStrCell(XL_underlying,C_underlying," ARM_ERR: underlying id: object expected",C_result);
	XL_readStrCell(XL_tAsset,C_tAsset," ARM_ERR: trigger asset id: object expected",C_result);
	XL_readNumCell(XL_maturity,C_maturity," ARM_ERR: maturity: date expected",C_result);
	XL_readNumCell(XL_barrier1,C_barrier1," ARM_ERR: barrier 1: numeric expected",C_result);
	XL_readStrCell(XL_upDownDouble,C_upDownDouble," ARM_ERR: up or down: string expected",C_result);
	XL_readStrCell(XL_inOut,C_inOut," ARM_ERR: in or out: string expected",C_result);
	XL_readNumCellWD(XL_triggerVar,C_triggerVar,C_triggerVar_default," ARM_ERR: trigger var: numeric expected",C_result);
	XL_readNumCellWD(XL_rebate,C_rebate,C_rebate_default," ARM_ERR: rebate: numeric expected",C_result);
	XL_readNumCellWD(XL_firstX,C_firstX,C_firstX_default," ARM_ERR: first exercise: date expected",C_result);
	XL_readStrCellWD(XL_exerciseTiming,C_exerciseTiming,C_exerciseTiming_default," ARM_ERR: ExerciseTiming: string expected",C_result);
	XL_readNumCellWD(XL_barrier2,C_barrier2,C_barrier2_default," ARM_ERR: barrier 2: numeric expected",C_result);
		
	if((upDownDoubleId = ARM_ConvUpDownDouble (C_upDownDouble, C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	if((inOutId = ARM_ConvInOut (C_inOut, C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	exerciseTimingId = ARM_ConvPayResetRule(C_exerciseTiming);

	long retCode;
	long objId;
	CCString prevClass;
	
	CCString curClass = LOCAL_BARRIER_CLASS;
	CCString stringId = GetLastCurCellEnvValue ();

	if(!stringId)
	{
		retCode = ARMLOCAL_CONSTBARRIER (LocalGetNumObjectId (C_underlying),
										 LocalGetNumObjectId (C_tAsset),
										 C_maturity,
										 C_barrier1,
										 C_barrier2,
										 upDownDoubleId,
										 inOutId,
										 C_triggerVar,
										 C_rebate,
										 C_firstX,
										 exerciseTimingId,
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
			retCode = ARMLOCAL_CONSTBARRIER (LocalGetNumObjectId (C_underlying),
											 LocalGetNumObjectId (C_tAsset),
											 C_maturity,
											 C_barrier1,
											 C_barrier2,
											 upDownDoubleId,
											 inOutId,
											 C_triggerVar,
											 C_rebate,
											 C_firstX,
											 exerciseTimingId,
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
			retCode = ARMLOCAL_CONSTBARRIER (LocalGetNumObjectId (C_underlying),
											 LocalGetNumObjectId (C_tAsset),
											 C_maturity,
											 C_barrier1,
											 C_barrier2,
											 upDownDoubleId,
											 inOutId,
											 C_triggerVar,
											 C_rebate,
											 C_firstX,
											 exerciseTimingId,
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
//	ARM_END();
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_CONSTBARRIER" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


__declspec(dllexport) LPXLOPER WINAPI Local_PXL_CONSTBARRIER (LPXLOPER XL_underlying,
															  LPXLOPER XL_tAsset,
															  LPXLOPER XL_maturity,
															  LPXLOPER XL_barrier1,
															  LPXLOPER XL_upDownDouble,
															  LPXLOPER XL_inOut,
															  LPXLOPER XL_triggerVar,
															  LPXLOPER XL_rebate,
															  LPXLOPER XL_firstX,
															  LPXLOPER XL_exerciseTiming,
															  LPXLOPER XL_barrier2)
{
	ADD_LOG("Local_PXL_CONSTBARRIER ");
	//ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;
	
	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{

	ARM_NOCALCIFWIZ();
	
	// C variable
	CCString C_underlying;
	CCString C_tAsset;

	double C_maturity;

	double C_barrier1;
	double C_barrier2;
	double C_barrier2_default = 0.0;

	CCString C_upDownDouble;
	long upDownDoubleId;

	CCString C_inOut;
	long inOutId;

	double C_triggerVar;
	double C_triggerVar_default = 0.0;

	double C_rebate;
	double C_rebate_default = 0.0;

	double C_firstX;
	double C_firstX_default = -1.0;

	CCString C_exerciseTiming;
	CCString C_exerciseTiming_default = "ADVANCE";
	long exerciseTimingId = K_ADVANCE;

	// error
	static int error;
	static char* reason = "";

	XL_readStrCell(XL_underlying,C_underlying," ARM_ERR: underlying id: object expected",C_result);
	XL_readStrCell(XL_tAsset,C_tAsset," ARM_ERR: trigger asset id: object expected",C_result);
	XL_readNumCell(XL_maturity,C_maturity," ARM_ERR: maturity: date expected",C_result);
	XL_readNumCell(XL_barrier1,C_barrier1," ARM_ERR: barrier 1: numeric expected",C_result);
	XL_readStrCell(XL_upDownDouble,C_upDownDouble," ARM_ERR: up or down: string expected",C_result);
	XL_readStrCell(XL_inOut,C_inOut," ARM_ERR: in or out: string expected",C_result);
	XL_readNumCellWD(XL_triggerVar,C_triggerVar,C_triggerVar_default," ARM_ERR: trigger var: numeric expected",C_result);
	XL_readNumCellWD(XL_rebate,C_rebate,C_rebate_default," ARM_ERR: rebate: numeric expected",C_result);
	XL_readNumCellWD(XL_firstX,C_firstX,C_firstX_default," ARM_ERR: first exercise: date expected",C_result);
	XL_readStrCellWD(XL_exerciseTiming,C_exerciseTiming,C_exerciseTiming_default," ARM_ERR: ExerciseTiming: string expected",C_result);
	XL_readNumCellWD(XL_barrier2,C_barrier2,C_barrier2_default," ARM_ERR: barrier 2: numeric expected",C_result);
		
	if((upDownDoubleId = ARM_ConvUpDownDouble (C_upDownDouble, C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	if((inOutId = ARM_ConvInOut (C_inOut, C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	exerciseTimingId = ARM_ConvPayResetRule(C_exerciseTiming);

	long retCode;
	long objId;
	
	CCString curClass = LOCAL_BARRIER_CLASS;
	CCString stringId;

	retCode = ARMLOCAL_CONSTBARRIER (LocalGetNumObjectId (C_underlying),
									 LocalGetNumObjectId (C_tAsset),
									 C_maturity,
									 C_barrier1,
									 C_barrier2,
									 upDownDoubleId,
									 inOutId,
									 C_triggerVar,
									 C_rebate,
									 C_firstX,
									 exerciseTimingId,
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_PXL_CONSTBARRIER" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


__declspec(dllexport) LPXLOPER WINAPI Local_BARRIER (LPXLOPER XL_underlying,
													 LPXLOPER XL_tAsset,
													 LPXLOPER XL_xStyle,
													 LPXLOPER XL_refVal1,
													 LPXLOPER XL_upDownDouble,
													 LPXLOPER XL_inOut,
													 LPXLOPER XL_triggerVar,
													 LPXLOPER XL_rebate,
													 LPXLOPER XL_exerciseTiming,
													 LPXLOPER XL_refVal2)
{
	ADD_LOG("Local_BARRIER ");
	//ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;
	
	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{

	ARM_NOCALCIFWIZ();
	
	// C variable
	CCString C_underlying;
	CCString C_tAsset;
	CCString C_xStyle;
	CCString C_refVal1;
	CCString C_refVal2;
	CCString C_refVal2_default = "";

	CCString C_upDownDouble;
	long upDownDoubleId;

	CCString C_inOut;
	long inOutId;

	double C_triggerVar;
	double C_triggerVar_default = 0.0;

	double C_rebate;
	double C_rebate_default = 0.0;

	CCString C_exerciseTiming;
	CCString C_exerciseTiming_default = "ADVANCE";
	long exerciseTimingId = K_ADVANCE;

	// error
	static int error;
	static char* reason = "";

	XL_readStrCell(XL_underlying,C_underlying," ARM_ERR: underlying id: object expected",C_result);
	XL_readStrCell(XL_tAsset,C_tAsset," ARM_ERR: trigger asset id: object expected",C_result);
	XL_readStrCell(XL_xStyle,C_xStyle," ARM_ERR: exercise style id: object expected",C_result);
	XL_readStrCell(XL_refVal1,C_refVal1," ARM_ERR: reference value id 1: object expected",C_result);
	XL_readStrCell(XL_upDownDouble,C_upDownDouble," ARM_ERR: up or down: string expected",C_result);
	XL_readStrCell(XL_inOut,C_inOut," ARM_ERR: in or out: string expected",C_result);
	XL_readNumCellWD(XL_triggerVar,C_triggerVar,C_triggerVar_default," ARM_ERR: trigger var: numeric expected",C_result);
	XL_readNumCellWD(XL_rebate,C_rebate,C_rebate_default," ARM_ERR: rebate: numeric expected",C_result);
	XL_readStrCellWD(XL_exerciseTiming,C_exerciseTiming,C_exerciseTiming_default," ARM_ERR: ExerciseTiming: string expected",C_result);
	XL_readStrCellWD(XL_refVal2,C_refVal2,C_refVal2_default," ARM_ERR: reference value id 2: object expected",C_result);
		
	if((upDownDoubleId = ARM_ConvUpDownDouble (C_upDownDouble, C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	if((inOutId = ARM_ConvInOut (C_inOut, C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	exerciseTimingId = ARM_ConvPayResetRule(C_exerciseTiming);

	long retCode;
	long objId;
	CCString prevClass;
	
	CCString curClass = LOCAL_BARRIER_CLASS;
	CCString stringId = GetLastCurCellEnvValue ();

	if(!stringId)
	{
		retCode = ARMLOCAL_BARRIER (LocalGetNumObjectId (C_underlying),
									LocalGetNumObjectId (C_tAsset),
									LocalGetNumObjectId (C_xStyle),
									LocalGetNumObjectId (C_refVal1),
									LocalGetNumObjectId (C_refVal2),
									upDownDoubleId,
									inOutId,
									C_triggerVar,
									C_rebate,
									exerciseTimingId,
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
			retCode = ARMLOCAL_BARRIER (LocalGetNumObjectId (C_underlying),
										LocalGetNumObjectId (C_tAsset),
										LocalGetNumObjectId (C_xStyle),
										LocalGetNumObjectId (C_refVal1),
										LocalGetNumObjectId (C_refVal2),
										upDownDoubleId,
										inOutId,
										C_triggerVar,
										C_rebate,
										exerciseTimingId,
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
			retCode = ARMLOCAL_BARRIER (LocalGetNumObjectId (C_underlying),
										LocalGetNumObjectId (C_tAsset),
										LocalGetNumObjectId (C_xStyle),
										LocalGetNumObjectId (C_refVal1),
										LocalGetNumObjectId (C_refVal2),
										upDownDoubleId,
										inOutId,
										C_triggerVar,
										C_rebate,
										exerciseTimingId,
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_BARRIER" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


__declspec(dllexport) LPXLOPER WINAPI Local_PXL_BARRIER (LPXLOPER XL_underlying,
														 LPXLOPER XL_tAsset,
														 LPXLOPER XL_xStyle,
														 LPXLOPER XL_refVal1,
														 LPXLOPER XL_upDownDouble,
														 LPXLOPER XL_inOut,
														 LPXLOPER XL_triggerVar,
														 LPXLOPER XL_rebate,
														 LPXLOPER XL_exerciseTiming,
														 LPXLOPER XL_refVal2)
{
	ADD_LOG("Local_PXL_BARRIER ");
	//ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;
	
	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{

	ARM_NOCALCIFWIZ();
	
	// C variable
	CCString C_underlying;
	CCString C_tAsset;
	CCString C_xStyle;
	CCString C_refVal1;
	CCString C_refVal2;
	CCString C_refVal2_default = "";

	CCString C_upDownDouble;
	long upDownDoubleId;

	CCString C_inOut;
	long inOutId;

	double C_triggerVar;
	double C_triggerVar_default = 0.0;

	double C_rebate;
	double C_rebate_default = 0.0;

	CCString C_exerciseTiming;
	CCString C_exerciseTiming_default = "ADVANCE";
	long exerciseTimingId = K_ADVANCE;

	// error
	static int error;
	static char* reason = "";

	XL_readStrCell(XL_underlying,C_underlying," ARM_ERR: underlying id: object expected",C_result);
	XL_readStrCell(XL_tAsset,C_tAsset," ARM_ERR: trigger asset id: object expected",C_result);
	XL_readStrCell(XL_xStyle,C_xStyle," ARM_ERR: exercise style id: object expected",C_result);
	XL_readStrCell(XL_refVal1,C_refVal1," ARM_ERR: reference value id 1: object expected",C_result);
	XL_readStrCell(XL_upDownDouble,C_upDownDouble," ARM_ERR: up or down: string expected",C_result);
	XL_readStrCell(XL_inOut,C_inOut," ARM_ERR: in or out: string expected",C_result);
	XL_readNumCellWD(XL_triggerVar,C_triggerVar,C_triggerVar_default," ARM_ERR: trigger var: numeric expected",C_result);
	XL_readNumCellWD(XL_rebate,C_rebate,C_rebate_default," ARM_ERR: rebate: numeric expected",C_result);
	XL_readStrCellWD(XL_exerciseTiming,C_exerciseTiming,C_exerciseTiming_default," ARM_ERR: ExerciseTiming: string expected",C_result);
	XL_readStrCellWD(XL_refVal2,C_refVal2,C_refVal2_default," ARM_ERR: reference value id 2: object expected",C_result);
		
	if((upDownDoubleId = ARM_ConvUpDownDouble (C_upDownDouble, C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	if((inOutId = ARM_ConvInOut (C_inOut, C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	exerciseTimingId = ARM_ConvPayResetRule(C_exerciseTiming);

	long retCode;
	long objId;
	
	CCString curClass = LOCAL_BARRIER_CLASS;
	CCString stringId;

	retCode = ARMLOCAL_BARRIER (LocalGetNumObjectId (C_underlying),
								LocalGetNumObjectId (C_tAsset),
								LocalGetNumObjectId (C_xStyle),
								LocalGetNumObjectId (C_refVal1),
								LocalGetNumObjectId (C_refVal2),
								upDownDoubleId,
								inOutId,
								C_triggerVar,
								C_rebate,
								exerciseTimingId,
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

	//ARM_END();
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_PXL_BARRIER" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


// EOF %M% 
