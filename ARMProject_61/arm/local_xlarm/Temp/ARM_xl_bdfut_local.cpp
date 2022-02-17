#pragma warning(disable : 4786)
#pragma warning(disable : 4005)

#include "ARM_xl_bdfut_local.h"

#include "ARM_local_interglob.h"
#include "ARM_local_interface.h"

#include "ARM_xl_trycatch_local.h"
#include <util\fromto.h>

#include <ARM\libarm_local\ARM_local_class.h>
#include <ARM\libarm_local\ARM_local_bdfut.h>

#include "util\tech_macro.h"


__declspec(dllexport) LPXLOPER WINAPI Local_BDFUT (LPXLOPER XL_delivery,
												   LPXLOPER XL_underlyingId,
												   LPXLOPER XL_coupon,
												   LPXLOPER XL_cfVal)
{
	ADD_LOG("Local_BDFUT ");
	//ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();

	// C variable
	double C_delivery;

	CCString C_underlyingId;

	double C_coupon;

	double C_cfVal;
	double C_cfVal_default = 1000000;

	// error
	static int error;
	static char* reason = "";

	XL_readNumCell(XL_delivery,C_delivery," ARM_ERR: delivery: date expected",C_result);
	XL_readStrCell(XL_underlyingId,C_underlyingId," ARM_ERR: underlying id: object expected",C_result);
	XL_readNumCell(XL_coupon,C_coupon," ARM_ERR: coupon: numeric expected",C_result);
	XL_readNumCellWD(XL_cfVal,C_cfVal,C_cfVal_default," ARM_ERR: conversion factor: numeric expected",C_result);
	
	long retCode;
	long objId;
	CCString prevClass;
	
	CCString curClass = LOCAL_BDFUT_CLASS;
	CCString stringId = GetLastCurCellEnvValue ();

	if(!stringId)
	{
		retCode = ARMLOCAL_BDFUT (C_delivery, 1, LocalGetNumObjectId (C_underlyingId), C_coupon, C_cfVal, C_result);

		if(retCode == ARM_OK)
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
			
		if(curClass == prevClass)
		{
			retCode = ARMLOCAL_BDFUT (C_delivery, 1, LocalGetNumObjectId (C_underlyingId), 
								 C_coupon, C_cfVal, C_result, objId);
		}
		else
		{
			FreeCurCellContent ();
			retCode = ARMLOCAL_BDFUT (C_delivery, 1, LocalGetNumObjectId (C_underlyingId),
				                 C_coupon, C_cfVal, C_result);

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

	//ARM_END();
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_BDFUT" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}

__declspec(dllexport) LPXLOPER WINAPI Local_GetConversionFactor (LPXLOPER XL_bdFutId,
																 LPXLOPER XL_factId)
{
	ADD_LOG("Local_GetConversionFactor ");
	//ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();

	// C variable
	CCString C_bdFutId;
	double C_factId;
	
	// error
	static int error;
	static char* reason = "";

	XL_readStrCell(XL_bdFutId,C_bdFutId," ARM_ERR: bond future id: object expected",C_result);
	XL_readNumCell(XL_factId,C_factId," ARM_ERR: factor id: numeric expected",C_result);
	
	long retCode = ARMLOCAL_GetConversionFactor (LocalGetNumObjectId (C_bdFutId), (long)C_factId, C_result);

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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_GetConversionFactor" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}

__declspec(dllexport) LPXLOPER WINAPI Local_PXL_BDFUT (LPXLOPER XL_delivery,
													   LPXLOPER XL_underlyingId,
													   LPXLOPER XL_coupon,
													   LPXLOPER XL_cfVal)
{
	ADD_LOG("Local_PXL_BDFUT ");
	//ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();

	// C variable
	double C_delivery;

	CCString C_underlyingId;

	double C_coupon;

	double C_cfVal;
	double C_cfVal_default = 1000000;

	// error
	static int error;
	static char* reason = "";

	XL_readNumCell(XL_delivery,C_delivery," ARM_ERR: delivery: date expected",C_result);
	XL_readStrCell(XL_underlyingId,C_underlyingId," ARM_ERR: underlying id: object expected",C_result);
	XL_readNumCell(XL_coupon,C_coupon," ARM_ERR: coupon: numeric expected",C_result);
	XL_readNumCellWD(XL_cfVal,C_cfVal,C_cfVal_default," ARM_ERR: conversion factor: numeric expected",C_result);
	
	long retCode;
	long objId;
	
	CCString curClass =LOCAL_BDFUT_CLASS;
	CCString stringId;

	retCode = ARMLOCAL_BDFUT (C_delivery, 1, LocalGetNumObjectId (C_underlyingId), C_coupon, C_cfVal, C_result);

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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_PXL_BDFUT" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}

__declspec(dllexport) LPXLOPER WINAPI Local_GetCheapest (LPXLOPER XL_bdFutId)
{
	ADD_LOG("Local_GetCheapest ");
	//ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();

	// C variable
	CCString C_bdFutId;
	
	// error
	static int error;
	static char* reason = "";

	XL_readStrCell(XL_bdFutId,C_bdFutId," ARM_ERR: bond future id: object expected",C_result);
	
	long retCode = ARMLOCAL_GetCheapest (LocalGetNumObjectId (C_bdFutId), C_result);

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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_GetCheapest" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}

__declspec(dllexport) LPXLOPER WINAPI Local_GILT_NOTIONNAL_BUND (LPXLOPER XL_delivery,
																 LPXLOPER XL_underlyingId,
																 LPXLOPER XL_bdFutType,
																 LPXLOPER XL_market)
{
	ADD_LOG("Local_GILT_NOTIONNAL_BUND ");
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

	CCString C_underlyingId;

	double C_bdFutType;
	double C_market;
	
	// error
	static int error;
	static char* reason = "";

	XL_readStrCell(XL_delivery,C_delivery," ARM_ERR: delivery: string expected",C_result);
	XL_readStrCell(XL_underlyingId,C_underlyingId," ARM_ERR: underlying id: object expected",C_result);
	XL_readNumCell(XL_bdFutType,C_bdFutType," ARM_ERR: bond future type: numeric expected",C_result);
	XL_readNumCell(XL_market,C_market," ARM_ERR: market: numeric expected",C_result);
	
	long retCode;
	long objId;
	CCString prevClass;
	
	CCString curClass =LOCAL_BDFUT_CLASS;
	CCString stringId = GetLastCurCellEnvValue ();

	if(!stringId)
	{
		retCode = ARMLOCAL_GILT_NOTIONNAL_BUND (C_delivery, LocalGetNumObjectId (C_underlyingId),
										   (long)C_bdFutType, (long)C_market, C_result);

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
			retCode = ARMLOCAL_GILT_NOTIONNAL_BUND (C_delivery, LocalGetNumObjectId (C_underlyingId),
										       (long)C_bdFutType, (long)C_market, C_result, objId);
		}
		else
		{
			FreeCurCellContent ();
			retCode = ARMLOCAL_GILT_NOTIONNAL_BUND (C_delivery, LocalGetNumObjectId (C_underlyingId),
										       (long)C_bdFutType, (long)C_market, C_result);

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

	//ARM_END();
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_GILT_NOTIONNAL_BUND" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}

__declspec(dllexport) LPXLOPER WINAPI Local_PXL_GILT_NOTIONNAL_BUND (LPXLOPER XL_delivery,
																	 LPXLOPER XL_underlyingId,
																	 LPXLOPER XL_bdFutType,
																	 LPXLOPER XL_market)
{
	ADD_LOG("Local_PXL_GILT_NOTIONNAL_BUND ");
	//ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();

	// C variable
	CCString C_delivery;

	CCString C_underlyingId;

	double C_bdFutType;
	double C_market;
	
	// error
	static int error;
	static char* reason = "";

	XL_readStrCell(XL_delivery,C_delivery," ARM_ERR: delivery: string expected",C_result);
	XL_readStrCell(XL_underlyingId,C_underlyingId," ARM_ERR: underlying id: object expected",C_result);
	XL_readNumCell(XL_bdFutType,C_bdFutType," ARM_ERR: bond future type: numeric expected",C_result);
	XL_readNumCell(XL_market,C_market," ARM_ERR: market: numeric expected",C_result);
	
	long retCode;
	long objId;
	
	CCString curClass =LOCAL_BDFUT_CLASS;
	CCString stringId;

	retCode = ARMLOCAL_GILT_NOTIONNAL_BUND (C_delivery, LocalGetNumObjectId (C_underlyingId),
									   (long)C_bdFutType, (long)C_market, C_result);

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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_PXL_GILT_NOTIONNAL_BUND" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}

__declspec(dllexport) LPXLOPER WINAPI Local_GILT (LPXLOPER XL_delivery,
												  LPXLOPER XL_underlyingId)
{
	ADD_LOG("Local_GILT ");
	//ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;
	
	// C variable
	CCString C_delivery;
	CCString C_underlyingId;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();

	// error
	static int error;
	static char* reason = "";

	XL_readStrCell(XL_delivery,C_delivery," ARM_ERR: delivery: string expected",C_result);
	XL_readStrCell(XL_underlyingId,C_underlyingId," ARM_ERR: underlying id: object expected",C_result);
	
	long retCode;
	long objId;
	CCString prevClass;
	
	CCString curClass =LOCAL_BDFUT_CLASS;
	CCString stringId = GetLastCurCellEnvValue ();

	if(!stringId)
	{
		retCode = ARMLOCAL_GILT (C_delivery, LocalGetNumObjectId (C_underlyingId), C_result);
						    
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
			retCode = ARMLOCAL_GILT (C_delivery, LocalGetNumObjectId (C_underlyingId), C_result, objId);
		}
		else
		{
			FreeCurCellContent ();
			retCode = ARMLOCAL_GILT (C_delivery, LocalGetNumObjectId (C_underlyingId), C_result);

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

	//ARM_END();
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_GILT" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}

__declspec(dllexport) LPXLOPER WINAPI Local_PXL_GILT (LPXLOPER XL_delivery,
													  LPXLOPER XL_underlyingId)
{
	ADD_LOG("Local_PXL_GILT ");
	//ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();

	// C variable
	CCString C_delivery;
	CCString C_underlyingId;

	// error
	static int error;
	static char* reason = "";

	XL_readStrCell(XL_delivery,C_delivery," ARM_ERR: delivery: string expected",C_result);
	XL_readStrCell(XL_underlyingId,C_underlyingId," ARM_ERR: underlying id: object expected",C_result);
	
	long retCode;
	long objId;
	
	CCString curClass =LOCAL_BDFUT_CLASS;
	CCString stringId;

	retCode = ARMLOCAL_GILT (C_delivery, LocalGetNumObjectId (C_underlyingId), C_result);
						    
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_PXL_GILT" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}

__declspec(dllexport) LPXLOPER WINAPI Local_NOTIONNAL (LPXLOPER XL_delivery,
													   LPXLOPER XL_underlyingId)
{
	ADD_LOG("Local_NOTIONNAL ");
	//ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();

	// C variable
	CCString C_delivery;
	CCString C_underlyingId;

	// error
	static int error;
	static char* reason = "";

	XL_readStrCell(XL_delivery,C_delivery," ARM_ERR: delivery: string expected",C_result);
	XL_readStrCell(XL_underlyingId,C_underlyingId," ARM_ERR: underlying id: object expected",C_result);
	
	long retCode;
	long objId;
	CCString prevClass;
	
	CCString curClass =LOCAL_BDFUT_CLASS;
	CCString stringId = GetLastCurCellEnvValue ();

	if(!stringId)
	{
		retCode = ARMLOCAL_NOTIONNAL (C_delivery, LocalGetNumObjectId (C_underlyingId), C_result);
						    
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
			retCode = ARMLOCAL_NOTIONNAL (C_delivery, LocalGetNumObjectId (C_underlyingId), C_result, objId);
		}
		else
		{
			FreeCurCellContent ();
			retCode = ARMLOCAL_NOTIONNAL (C_delivery, LocalGetNumObjectId (C_underlyingId), C_result);

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

	//ARM_END();
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_NOTIONNAL" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}

__declspec(dllexport) LPXLOPER WINAPI Local_PXL_NOTIONNAL (LPXLOPER XL_delivery,
														   LPXLOPER XL_underlyingId)
{
	ADD_LOG("Local_PXL_NOTIONNAL ");
	//ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();

	// C variable
	CCString C_delivery;
	CCString C_underlyingId;

	// error
	static int error;
	static char* reason = "";

	XL_readStrCell(XL_delivery,C_delivery," ARM_ERR: delivery: string expected",C_result);
	XL_readStrCell(XL_underlyingId,C_underlyingId," ARM_ERR: underlying id: object expected",C_result);
	
	long retCode;
	long objId;
	
	CCString curClass =LOCAL_BDFUT_CLASS;
	CCString stringId;

	retCode = ARMLOCAL_NOTIONNAL (C_delivery, LocalGetNumObjectId (C_underlyingId), C_result);
						    
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_PXL_NOTIONNAL" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}

__declspec(dllexport) LPXLOPER WINAPI Local_BUND_LIFFE (LPXLOPER XL_delivery,
														LPXLOPER XL_underlyingId)
{
	ADD_LOG("Local_BUND_LIFFE ");
	//ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();

	// C variable
	CCString C_delivery;
	CCString C_underlyingId;

	// error
	static int error;
	static char* reason = "";

	XL_readStrCell(XL_delivery,C_delivery," ARM_ERR: delivery: string expected",C_result);
	XL_readStrCell(XL_underlyingId,C_underlyingId," ARM_ERR: underlying id: object expected",C_result);

	long retCode;
	long objId;
	CCString prevClass;
	
	CCString curClass =LOCAL_BDFUT_CLASS;
	CCString stringId = GetLastCurCellEnvValue ();

	if(!stringId)
	{
		retCode = ARMLOCAL_BUND_LIFFE (C_delivery, LocalGetNumObjectId (C_underlyingId), C_result);
						    
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
			retCode = ARMLOCAL_BUND_LIFFE (C_delivery, LocalGetNumObjectId (C_underlyingId), C_result, objId);
		}
		else
		{
			FreeCurCellContent ();
			retCode = ARMLOCAL_BUND_LIFFE (C_delivery, LocalGetNumObjectId (C_underlyingId), C_result);

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

	//ARM_END();
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_BUND_LIFFE" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}

__declspec(dllexport) LPXLOPER WINAPI Local_PXL_BUND_LIFFE (LPXLOPER XL_delivery,
															LPXLOPER XL_underlyingId)
{
	ADD_LOG("Local_PXL_BUND_LIFFE ");
	//ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();

	// C variable
	double C_delivery;
	CCString C_underlyingId;

	// error
	static int error;
	static char* reason = "";

	XL_readNumCell(XL_delivery,C_delivery," ARM_ERR: delivery: date expected",C_result);
	XL_readStrCell(XL_underlyingId,C_underlyingId," ARM_ERR: underlying id: object expected",C_result);
	
	long retCode;
	long objId;
	
	CCString curClass =LOCAL_BDFUT_CLASS;
	CCString stringId;

	retCode = ARMLOCAL_BUND_LIFFE (C_delivery, LocalGetNumObjectId (C_underlyingId), C_result);
						    
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_PXL_BUND_LIFFE" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}

__declspec(dllexport) LPXLOPER WINAPI Local_BUND_DTB (LPXLOPER XL_delivery,
													  LPXLOPER XL_underlyingId)
{
	ADD_LOG("Local_BUND_DTB ");
	//ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();

	// C variable
	CCString C_delivery;
	CCString C_underlyingId;

	// error
	static int error;
	static char* reason = "";

	XL_readStrCell(XL_delivery,C_delivery," ARM_ERR: delivery: date expected",C_result);
	XL_readStrCell(XL_underlyingId,C_underlyingId," ARM_ERR: underlying id: object expected",C_result);
	
	long retCode;
	long objId;
	CCString prevClass;
	
	CCString curClass =LOCAL_BDFUT_CLASS;
	CCString stringId = GetLastCurCellEnvValue ();

	if(!stringId)
	{
		retCode = ARMLOCAL_BUND_DTB (C_delivery, LocalGetNumObjectId (C_underlyingId), C_result);
						    
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
			retCode = ARMLOCAL_BUND_DTB (C_delivery, LocalGetNumObjectId (C_underlyingId), C_result, objId);
		}
		else
		{
			FreeCurCellContent ();
			retCode = ARMLOCAL_BUND_DTB (C_delivery, LocalGetNumObjectId (C_underlyingId), C_result);

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

	//ARM_END();
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_BUND_DTB" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}

__declspec(dllexport) LPXLOPER WINAPI Local_PXL_BUND_DTB (LPXLOPER XL_delivery,
														  LPXLOPER XL_underlyingId)
{
	ADD_LOG("Local_PXL_BUND_DTB ");
	//ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();

	// C variable
	double C_delivery;
	CCString C_underlyingId;

	// error
	static int error;
	static char* reason = "";

	XL_readNumCell(XL_delivery,C_delivery," ARM_ERR: delivery: date expected",C_result);
	XL_readStrCell(XL_underlyingId,C_underlyingId," ARM_ERR: underlying id: object expected",C_result);
	
	long retCode;
	long objId;
	
	CCString curClass =LOCAL_BDFUT_CLASS;
	CCString stringId;

	retCode = ARMLOCAL_BUND_DTB (C_delivery, LocalGetNumObjectId (C_underlyingId), C_result);
						    
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_PXL_BUND_DTB" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}

__declspec(dllexport) LPXLOPER WINAPI Local_BDFUTBASKET (LPXLOPER XL_delivery,
														 LPXLOPER XL_underlyingId,
														 LPXLOPER XL_coupon,
														 LPXLOPER XL_futurePrice)
{
	ADD_LOG("Local_BDFUTBASKET ");
	//ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();

	// C variable
	double C_delivery;
	CCString C_underlyingId;

	double C_coupon;
	double C_futurePrice;

	// error
	static int error;
	static char* reason = "";

	XL_readNumCell(XL_delivery,C_delivery," ARM_ERR: delivery: date expected",C_result);
	XL_readStrCell(XL_underlyingId,C_underlyingId," ARM_ERR: underlying id: object expected",C_result);
	XL_readNumCell(XL_coupon,C_coupon," ARM_ERR: coupon: numeric expected",C_result);
	XL_readNumCell(XL_futurePrice,C_futurePrice," ARM_ERR: future price: numeric expected",C_result);
	
	long retCode;
	long objId;
	CCString prevClass;
	
	CCString curClass =LOCAL_BDFUT_CLASS;
	CCString stringId = GetLastCurCellEnvValue ();

	if(!stringId)
	{
		retCode = ARMLOCAL_BDFUTBASKET (C_delivery, 0, LocalGetNumObjectId (C_underlyingId), 
								   C_coupon, C_futurePrice, C_result);
						    
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
			retCode = ARMLOCAL_BDFUTBASKET (C_delivery, 0, LocalGetNumObjectId (C_underlyingId), 
								       C_coupon, C_futurePrice, C_result, objId);
		}
		else
		{
			FreeCurCellContent ();
			retCode = ARMLOCAL_BDFUTBASKET (C_delivery, 0, LocalGetNumObjectId (C_underlyingId), 
								       C_coupon, C_futurePrice, C_result);
		
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

	//ARM_END();
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_BDFUTBASKET" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}

__declspec(dllexport) LPXLOPER WINAPI Local_PXL_BDFUTBASKET (LPXLOPER XL_delivery,
															 LPXLOPER XL_underlyingId,
															 LPXLOPER XL_coupon,
															 LPXLOPER XL_futurePrice)
{
	ADD_LOG("Local_PXL_BDFUTBASKET ");
	//ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();

	// C variable
	double C_delivery;
	CCString C_underlyingId;

	double C_coupon;
	double C_futurePrice;

	// error
	static int error;
	static char* reason = "";

	XL_readNumCell(XL_delivery,C_delivery," ARM_ERR: delivery: date expected",C_result);
	XL_readStrCell(XL_underlyingId,C_underlyingId," ARM_ERR: underlying id: object expected",C_result);
	XL_readNumCell(XL_coupon,C_coupon," ARM_ERR: coupon: numeric expected",C_result);
	XL_readNumCell(XL_futurePrice,C_futurePrice," ARM_ERR: future price: numeric expected",C_result);
	
	long retCode;
	long objId;
	
	CCString curClass =LOCAL_BDFUT_CLASS;
	CCString stringId;

	retCode = ARMLOCAL_BDFUTBASKET (C_delivery, 0, LocalGetNumObjectId (C_underlyingId), 
							   C_coupon, C_futurePrice, C_result);
						    
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_PXL_BDFUTBASKET" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


					  
// EOF %M% 
