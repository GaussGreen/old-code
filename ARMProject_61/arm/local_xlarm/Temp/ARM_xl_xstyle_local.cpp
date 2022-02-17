#pragma warning(disable : 4786)
#pragma warning(disable : 4005)

#include <libCCxll\CCxll.h>

#include <ARM\libarm_local\ARM_local_class.h>
#include <ARM\libarm_local\ARM_local_glob.h>
#include <ARM\libarm_local\ARM_local_xstyle.h>

#include "ARM_local_interglob.h"
#include "ARM_local_interface.h"
#include "ARM_xl_trycatch_local.h"

#include "util\tech_macro.h"

__declspec(dllexport) LPXLOPER WINAPI Local_BERMUDANXSTYLE (LPXLOPER XL_xDates, LPXLOPER XL_xExpiryDates)
{
	ADD_LOG("Local_BERMUDANXSTYLE ");
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();

	// C variable
	VECTOR<double> C_xDates;
	VECTOR<double> C_xExpiryDates;

	
	// error
	static int error;
	static char* reason = "";

	XL_readNumVector(XL_xDates,C_xDates," ARM_ERR: exercise dates: array of date expected",C_result);
	XL_readNumVectorWD(XL_xExpiryDates, C_xExpiryDates, C_xExpiryDates,
					   " ARM_ERR: expiry dates: array of date expected", C_result);

	if( C_xExpiryDates.size() != 0 )
	{
		if( C_xExpiryDates.size() != C_xDates.size() )
		{
			C_result.setMsg("Array of notice dates and array of expiry dates should have the same size");
			ARM_ERR();
		}
	}

	long retCode;
	long objId;
	CCString prevClass;
	
	CCString curClass = LOCAL_XSTYLE_CLASS;
	CCString stringId = GetLastCurCellEnvValue ();

	if(!stringId)
	{
		retCode = ARMLOCAL_BERMUDANXSTYLE (C_xDates, C_xExpiryDates, C_result);

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
			retCode = ARMLOCAL_BERMUDANXSTYLE (C_xDates, C_xExpiryDates, C_result, objId);

			if(retCode == ARM_OK)
			{
				LocalSetCurCellEnvValue (curClass, objId); 

				stringId = LocalMakeObjectId (objId, curClass);
			}
		}
		else
		{
			FreeCurCellContent ();
			retCode = ARMLOCAL_BERMUDANXSTYLE (C_xDates, C_xExpiryDates, C_result);

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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_BERMUDANXSTYLE" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


__declspec(dllexport) LPXLOPER WINAPI Local_PXL_BERMUDANXSTYLE (LPXLOPER XL_xDates, LPXLOPER XL_xExpiryDates)
{
	ADD_LOG("Local_PXL_BERMUDANXSTYLE ");
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();

	// C variable
	VECTOR<double> C_xDates;
	VECTOR<double> C_xExpiryDates;
	
	// error
	static int error;
	static char* reason = "";

	XL_readNumVector(XL_xDates,C_xDates," ARM_ERR: exercise dates: array of date expected",C_result);
	XL_readNumVectorWD(XL_xExpiryDates, C_xExpiryDates, C_xExpiryDates,
					   " ARM_ERR: expiry dates: array of date expected", C_result);

	if( C_xExpiryDates.size() != 0 )
	{
		if( C_xExpiryDates.size() != C_xDates.size() )
		{
			C_result.setMsg("Array of notice dates and array of expiry dates should have the same size");
			ARM_ERR();
		}
	}

	long retCode;
	long objId;
	CCString prevClass;
	
	CCString curClass = LOCAL_XSTYLE_CLASS;
	CCString stringId;

	retCode = ARMLOCAL_BERMUDANXSTYLE (C_xDates, C_xExpiryDates, C_result);

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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_PXL_BERMUDANXSTYLE" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


__declspec(dllexport) LPXLOPER WINAPI Local_EUROPEANXSTYLE (LPXLOPER XL_xdate)
{
	ADD_LOG("Local_EUROPEANXSTYLE ");
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();

	// C variable
	double C_xdate;

	// error
	static int error;
	static char* reason = "";

	XL_readNumCell(XL_xdate,C_xdate," ARM_ERR: exercise date: date expected",C_result);
	
	long retCode;
	long objId;
	CCString prevClass;
	
	CCString curClass = LOCAL_XSTYLE_CLASS;
	CCString stringId = GetLastCurCellEnvValue ();

	if(!stringId)
	{
		retCode = ARMLOCAL_EUROPEANXSTYLE (C_xdate, C_result);

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
			retCode = ARMLOCAL_EUROPEANXSTYLE (C_xdate, C_result, objId);

			if(retCode == ARM_OK)
			{
				LocalSetCurCellEnvValue (curClass, objId); 

				stringId = LocalMakeObjectId (objId, curClass);
			}
		}
		else
		{
			FreeCurCellContent ();
			retCode = ARMLOCAL_EUROPEANXSTYLE (C_xdate, C_result);

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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_EUROPEANXSTYLE" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


__declspec(dllexport) LPXLOPER WINAPI Local_PXL_EUROPEANXSTYLE (LPXLOPER XL_xdate)
{
	ADD_LOG("Local_PXL_EUROPEANXSTYLE ");
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();

	// C variable
	double C_xdate;

	// error
	static int error;
	static char* reason = "";

	XL_readNumCell(XL_xdate,C_xdate," ARM_ERR: exercise date: date expected",C_result);
	
	long retCode;
	long objId;
	CCString prevClass;
	
	CCString curClass = LOCAL_XSTYLE_CLASS;
	CCString stringId;

	retCode = ARMLOCAL_EUROPEANXSTYLE (C_xdate, C_result);

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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_PXL_EUROPEANXSTYLE" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


__declspec(dllexport) LPXLOPER WINAPI Local_AMERICANXSTYLE (LPXLOPER XL_xStartDate,
															LPXLOPER XL_xEndDate)
{
	ADD_LOG("Local_AMERICANXSTYLE ");
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();

	// C variable
	double C_xStartDate;
	double C_xEndDate;

	// error
	static int error;
	static char* reason = "";

	XL_readNumCell(XL_xStartDate,C_xStartDate," ARM_ERR: exercise start date: date expected",C_result);
	XL_readNumCell(XL_xEndDate,C_xEndDate," ARM_ERR: exercise end date: date expected",C_result);
		
	long retCode;
	long objId;
	CCString prevClass;
	
	CCString curClass = LOCAL_XSTYLE_CLASS;
	CCString stringId = GetLastCurCellEnvValue ();

	if(!stringId)
	{
		retCode = ARMLOCAL_AMERICANXSTYLE (C_xStartDate, C_xEndDate, C_result);

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
			retCode = ARMLOCAL_AMERICANXSTYLE (C_xStartDate, C_xEndDate, C_result, objId);

			if(retCode == ARM_OK)
			{
				objId = C_result.getLong ();

				stringId = LocalMakeObjectId (objId, curClass);
			}
		}
		else
		{
			FreeCurCellContent ();

			retCode = ARMLOCAL_AMERICANXSTYLE (C_xStartDate, C_xEndDate, C_result);

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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_AMERICANXSTYLE" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}

__declspec(dllexport) LPXLOPER WINAPI Local_PXL_AMERICANXSTYLE (LPXLOPER XL_xStartDate,
																LPXLOPER XL_xEndDate)
{
	ADD_LOG("Local_PXL_AMERICANXSTYLE ");
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();

	// C variable
	double C_xStartDate;
	double C_xEndDate;

	// error
	static int error;
	static char* reason = "";

	XL_readNumCell(XL_xStartDate,C_xStartDate," ARM_ERR: exercise start date: date expected",C_result);
	XL_readNumCell(XL_xEndDate,C_xEndDate," ARM_ERR: exercise end date: date expected",C_result);
		
	long retCode;
	long objId;
	
	CCString curClass = LOCAL_XSTYLE_CLASS;
	CCString stringId;

	retCode = ARMLOCAL_AMERICANXSTYLE (C_xStartDate, C_xEndDate, C_result);

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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_PXL_AMERICANXSTYLE" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}
