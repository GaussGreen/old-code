#pragma warning(disable : 4786)
#pragma warning(disable : 4005)

#include <libCCxll\CCxll.h>

#include <ARM\libarm_local\ARM_local_glob.h>
#include <ARM\libarm_local\ARM_local_class.h>
#include <ARM\libarm_local\ARM_local_ccy.h>

#include "ARM_local_interglob.h"
#include "ARM_local_interface.h"

#include "ARM_xl_trycatch_local.h"
#include <util\fromto.h>

#include "util\tech_macro.h"

__declspec(dllexport) LPXLOPER WINAPI Local_ISOCCY (LPXLOPER XL_name)
{
	ADD_LOG("Local_ISOCCY ");
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();

	// C variable
	CCString C_name;
	
	// error
	static int error;
	static char* reason = "";

	XL_readStrCell(XL_name,C_name," ARM_ERR: name: string expected",C_result);
	
	long retCode;
	long objId;
	CCString prevClass;
	
	CCString curClass = LOCAL_CCY_CLASS;
	CCString stringId = GetLastCurCellEnvValue ();

	if(!stringId)
	{
		retCode = ARMLOCAL_ISOCCY (C_name,C_result);

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
			
		if(curClass == prevClass)
		{
			retCode = ARMLOCAL_ISOCCY (C_name,C_result,objId);

			if (retCode == ARM_OK)
			{
				LocalSetCurCellEnvValue (curClass, objId); 
				
				stringId = LocalMakeObjectId (objId, curClass);
			}
		}
		else
		{
			FreeCurCellContent ();

			retCode = ARMLOCAL_ISOCCY (C_name,C_result);
			
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_ISOCCY" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


__declspec(dllexport) LPXLOPER WINAPI Local_PXL_ISOCCY (LPXLOPER XL_name)
{
	ADD_LOG("Local_PXL_ISOCCY ");
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();

	// C variable
	CCString C_name;
	
	// error
	static int error;
	static char* reason = "";

	XL_readStrCell(XL_name,C_name," ARM_ERR: name: string expected",C_result);
	
	long retCode;
	long objId;
	
	CCString curClass (LOCAL_CCY_CLASS);
	CCString stringId;
	
	retCode = ARMLOCAL_ISOCCY (C_name, C_result);

	if(retCode == ARM_OK)
	{
		objId = C_result.getLong ();

		stringId = LocalMakeObjectId (objId, curClass);
	}
	
	if ( retCode == ARM_OK )
	{
		// No need: FreeCurCellErr();
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_PXL_ISOCCY" )

	/// return the result as an LPXLOPER	
	return (LPXLOPER)&XL_result;
}



__declspec(dllexport) LPXLOPER WINAPI Local_CCY (LPXLOPER XL_name,
												 LPXLOPER XL_idCurve,
												 LPXLOPER XL_crossValue,
												 LPXLOPER XL_dayCount)
{
	ADD_LOG("Local_CCY ");
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{

	ARM_NOCALCIFWIZ();

	// C variable
	CCString C_name;
	CCString C_idCurve;
	double C_crossValue;

	CCString C_dayCount;
    long dayCountId;

	// error
	static int error;
	static char* reason = "";

	XL_readStrCell(XL_name,C_name," ARM_ERR: name: string expected",C_result);
	XL_readStrCell(XL_idCurve,C_idCurve," ARM_ERR: curve id: string expected",C_result);
	XL_readNumCell(XL_crossValue,C_crossValue," ARM_ERR: cross value: numeric expected",C_result);
	XL_readStrCellWD(XL_dayCount,C_dayCount,"-1"," ARM_ERR: day count: String expected",C_result);

    dayCountId = ARM_ConvDayCount (C_dayCount);

	long retCode;
	long objId;
	CCString prevClass;
	
	CCString curClass = LOCAL_CCY_CLASS;
	CCString stringId = GetLastCurCellEnvValue ();

	if (!stringId)
	{
		retCode = ARMLOCAL_CCY (C_name, LocalGetNumObjectId (C_idCurve), C_crossValue, 
						   dayCountId, C_result);

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
			retCode = ARMLOCAL_CCY (C_name, LocalGetNumObjectId (C_idCurve), C_crossValue, 
						       dayCountId, C_result, objId);

			if(retCode == ARM_OK)
			{			
				LocalSetCurCellEnvValue (curClass, objId); 

				stringId = LocalMakeObjectId (objId, curClass);
			}
		}
		else
		{
			FreeCurCellContent ();
			retCode = ARMLOCAL_CCY (C_name, LocalGetNumObjectId (C_idCurve), C_crossValue, 
						       dayCountId, C_result);

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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_CCY" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}

__declspec(dllexport) LPXLOPER WINAPI Local_PXL_CCY (LPXLOPER XL_name,
													 LPXLOPER XL_idCurve,
													 LPXLOPER XL_crossValue,
													 LPXLOPER XL_dayCount)
{
	ADD_LOG("Local_PXL_CCY ");
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();

	// C variable
	CCString C_name;
	CCString C_idCurve;
	double C_crossValue;
	CCString C_dayCount;

    long dayCountId;

	// error
	static int error;
	static char* reason = "";

	XL_readStrCell(XL_name,C_name," ARM_ERR: name: string expected",C_result);
	XL_readStrCell(XL_idCurve,C_idCurve," ARM_ERR: curve id: string expected",C_result);
	XL_readNumCell(XL_crossValue,C_crossValue," ARM_ERR: cross value: numeric expected",C_result);
	XL_readStrCellWD(XL_dayCount,C_dayCount,"-1"," ARM_ERR: day count: String expected",C_result);

    dayCountId = ARM_ConvDayCount(C_dayCount);

	long retCode;
	long objId;
	
	CCString curClass = LOCAL_CCY_CLASS;
	CCString stringId;
	
	retCode = ARMLOCAL_CCY(C_name, LocalGetNumObjectId(C_idCurve), C_crossValue, 
					  dayCountId, C_result);

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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_PXL_CCY" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


__declspec(dllexport) LPXLOPER WINAPI Local_ARM_GetInfoFromCcy (LPXLOPER XL_ccyId,
																LPXLOPER XL_type)
{
	ADD_LOG("Local_ARM_GetInfoFromCcy ");
	long retCode;

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{

	ARM_NOCALCIFWIZ();

	// C variable
	CCString C_ccyId;
	CCString C_type;

	// error
	static int error;
	static char* reason = "";

	XL_readStrCell(XL_ccyId,C_ccyId," ARM_ERR: ccyId: object expected",C_result);
	XL_readStrCell(XL_type,C_type," ARM_ERR: type: string expected",C_result);

	C_type.toUpper ();

	retCode = ARMLOCAL_GetInfoFromCcy(LocalGetNumObjectId(C_ccyId),C_type,C_result);

	if(retCode == ARM_OK)
	{
		FreeCurCellErr ();
		XL_result.xltype = xltypeStr;
		XL_result.val.str = XL_StrC2StrPascal (C_result.getString());
		XL_result.xltype |= xlbitDLLFree;
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_ARM_GETINFOFROMPRCS" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}
