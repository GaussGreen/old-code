#pragma warning(disable : 4786)
#pragma warning(disable : 4005)

#include <libCCxll\CCxll.h>

#include "ARM_local_interglob.h"
#include "ARM_local_interface.h"

#include "ARM_xl_trycatch_local.h"
#include <util\fromto.h>

#include <ARM\libarm_local\ARM_local_class.h>
#include <ARM\libarm_local\ARM_local_irfut.h>

#include "util\tech_macro.h"


__declspec(dllexport) LPXLOPER WINAPI Local_IRFUT (LPXLOPER XL_delivery,
												   LPXLOPER XL_underlying)
{
	ADD_LOG("Local_IRFUT ");
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();


	// C variable
	double C_delivery;
	CCString C_underlying;
	
	// error
	static int error;
	static char* reason = "";

	XL_readNumCell(XL_delivery,C_delivery," ARM_ERR: delivery: date expected",C_result);
	XL_readStrCell(XL_underlying,C_underlying," ARM_ERR: underlying id: object expected",C_result);
	
	long retCode;
	long objId;
	CCString prevClass;
	
	CCString curClass = LOCAL_IRFUT_CLASS;
	CCString stringId = GetLastCurCellEnvValue ();
	
	if(!stringId)
	{
		retCode = ARMLOCAL_IRFUT (C_delivery, LocalGetNumObjectId (C_underlying), C_result);

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
			retCode = ARMLOCAL_IRFUT (C_delivery, LocalGetNumObjectId (C_underlying), C_result, objId);

			if(retCode == ARM_OK)
			{
				LocalSetCurCellEnvValue (curClass, objId); 

				stringId = LocalMakeObjectId (objId, curClass);
			}
		}
		else
		{
			FreeCurCellContent ();
			retCode = ARMLOCAL_IRFUT (C_delivery, LocalGetNumObjectId (C_underlying), C_result);
		
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_IRFUT" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


__declspec(dllexport) LPXLOPER WINAPI Local_PXL_IRFUT (LPXLOPER XL_delivery,
													   LPXLOPER XL_underlying)
{
	ADD_LOG("Local_PXL_IRFUT ");
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();


	// C variable
	double C_delivery;
	CCString C_underlying;
	
	// error
	static int error;
	static char* reason = "";

	XL_readNumCell(XL_delivery,C_delivery," ARM_ERR: delivery: date expected",C_result);
	XL_readStrCell(XL_underlying,C_underlying," ARM_ERR: underlying id: object expected",C_result);
	
	long retCode;
	long objId;
	
	CCString curClass = LOCAL_IRFUT_CLASS;
	CCString stringId;
	
	retCode = ARMLOCAL_IRFUT (C_delivery, LocalGetNumObjectId (C_underlying), C_result);

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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_PXL_IRFUT" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


__declspec(dllexport) LPXLOPER WINAPI Local_THREE_MONTH_FUT (LPXLOPER XL_delivery,
															 LPXLOPER XL_market,
															 LPXLOPER XL_ccy)
{
	ADD_LOG("Local_THREE_MONTH_FUT ");
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();


	// C variable
	CCString C_delivery;
	double C_market;
	CCString C_ccy;
	
	// error
	static int error;
	static char* reason = "";

	XL_readStrCell(XL_delivery,C_delivery," ARM_ERR: delivery: string expected",C_result);
	XL_readNumCell(XL_market,C_market," ARM_ERR: market: numeric expected",C_result);
	XL_readStrCell(XL_ccy,C_ccy," ARM_ERR: currency: string expected",C_result);
		
	long retCode;
	long objId;
	CCString prevClass;
	
	CCString curClass = LOCAL_IRFUT_CLASS;
	CCString stringId = GetLastCurCellEnvValue ();
	
	if(!stringId)
	{
		retCode = ARMLOCAL_THREE_MONTH_FUT (C_delivery, (long)C_market, C_ccy, C_result);

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
			retCode = ARMLOCAL_THREE_MONTH_FUT (C_delivery, (long)C_market, C_ccy, C_result, objId);

			if(retCode == ARM_OK)
			{
				LocalSetCurCellEnvValue (curClass, objId); 

				stringId = LocalMakeObjectId (objId, curClass);
			}
		}
		else
		{
			FreeCurCellContent ();
			retCode = ARMLOCAL_THREE_MONTH_FUT (C_delivery, (long)C_market, C_ccy, C_result);
		
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_THREE_MONTH_FUT" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


__declspec(dllexport) LPXLOPER WINAPI Local_PXL_THREE_MONTH_FUT (LPXLOPER XL_delivery,
																 LPXLOPER XL_market,
																 LPXLOPER XL_ccy)
{
	ADD_LOG("Local_PXL_THREE_MONTH_FUT ");
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();


	// C variable
	CCString C_delivery;
	double C_market;
	CCString C_ccy;
	
	// error
	static int error;
	static char* reason = "";

	XL_readStrCell(XL_delivery,C_delivery," ARM_ERR: delivery: string expected",C_result);
	XL_readNumCell(XL_market,C_market," ARM_ERR: market: numeric expected",C_result);
	XL_readStrCell(XL_ccy,C_ccy," ARM_ERR: currency: string expected",C_result);
		
	long retCode;
	long objId;
	
	CCString curClass = LOCAL_IRFUT_CLASS;
	CCString stringId;
	
	retCode = ARMLOCAL_THREE_MONTH_FUT (C_delivery, (long)C_market, C_ccy, C_result);

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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_PXL_THREE_MONTH_FUT" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}
