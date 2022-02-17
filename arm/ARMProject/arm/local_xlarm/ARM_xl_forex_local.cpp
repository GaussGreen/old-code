#pragma warning(disable : 4786)
#pragma warning(disable : 4005)

#include <libCCxll\CCxll.h>

#include <ARM\libarm_local\ARM_local_glob.h>
#include <ARM\libarm_local\ARM_local_class.h>
#include <ARM\libarm_local\ARM_local_forex.h>

#include "ARM_local_interglob.h"
#include "ARM_local_interface.h"
#include "ARM_xl_trycatch_local.h"

#include "util\tech_macro.h"

__declspec(dllexport) LPXLOPER WINAPI Local_FOREX (LPXLOPER XL_LeftCcy,
												   LPXLOPER XL_RightCcy,
                                                   LPXLOPER XL_SpotValue)
{
	ADD_LOG("Local_FOREX ");
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();

	// C variable
	CCString C_LeftCcy;
	CCString C_RightCcy;

	// error
	static int error;
	static char* reason = "";

	XL_readStrCell(XL_LeftCcy,C_LeftCcy," ARM_ERR: left currency: string expected",C_result);
	XL_readStrCell(XL_RightCcy,C_RightCcy," ARM_ERR: right currency: string expected",C_result);

	int sizeLeft = C_LeftCcy.GetLen();
	int sizeRight = C_RightCcy.GetLen();

	double spotValue;
	double spotValueDef=1.0;
	XL_readNumCellWD(XL_SpotValue,spotValue,spotValueDef," ARM_ERR: spot value: numerical expected",C_result);
		
	long retCode;
	long objId;
	CCString prevClass;
	
	CCString curClass = LOCAL_FOREX_CLASS;
	CCString stringId = GetLastCurCellEnvValue ();

	if(!stringId)
	{
		if(sizeLeft == 3 && sizeRight == 3)
		{
			retCode = ARMLOCAL_FOREX (C_LeftCcy,
									  C_RightCcy,
									  spotValue,
									  C_result);
			
		}
		else
		{
			retCode = ARMLOCAL_FOREX (LocalGetNumObjectId (C_LeftCcy),
									  LocalGetNumObjectId (C_RightCcy),
									  spotValue,
									  C_result);
		}

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
			if(sizeLeft == 3 && sizeRight == 3)
			{
				retCode = ARMLOCAL_FOREX (C_LeftCcy,
										  C_RightCcy,
										  spotValue,
										  C_result,
										  objId);
				
			}
			else
			{
				retCode = ARMLOCAL_FOREX (LocalGetNumObjectId (C_LeftCcy),
										  LocalGetNumObjectId (C_RightCcy),
										  spotValue,
										  C_result,
										  objId);
			}

			if(retCode == ARM_OK)
			{
				LocalSetCurCellEnvValue (curClass, objId); 

				stringId = LocalMakeObjectId (objId, curClass);
			}
		}
		else
		{
			FreeCurCellContent ();
			if(sizeLeft == 3 && sizeRight == 3)
			{
				retCode = ARMLOCAL_FOREX (C_LeftCcy,
										  C_RightCcy,
										  spotValue,
										  C_result);
				
			}
			else
			{
				retCode = ARMLOCAL_FOREX (LocalGetNumObjectId (C_LeftCcy),
										  LocalGetNumObjectId (C_RightCcy),
										  spotValue,
										  C_result);
			}

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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_FOREX" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}

__declspec(dllexport) LPXLOPER WINAPI Local_PXL_FOREX (LPXLOPER XL_LeftCcy,
													   LPXLOPER XL_RightCcy,
                                                       LPXLOPER XL_SpotValue)
{
	ADD_LOG("Local_PXL_FOREX ");
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();

	// C variable
	CCString C_LeftCcy;
	CCString C_RightCcy;

	// error
	static int error;
	static char* reason = "";

	XL_readStrCell(XL_LeftCcy,C_LeftCcy," ARM_ERR: left currency: string expected",C_result);
	XL_readStrCell(XL_RightCcy,C_RightCcy," ARM_ERR: right currency: string expected",C_result);

    int sizeLeft = C_LeftCcy.GetLen();
	int sizeRight = C_RightCcy.GetLen();
		
	double spotValue;
	double spotValueDef=1.0;
	XL_readNumCellWD(XL_SpotValue,spotValue,spotValueDef," ARM_ERR: spot value: numerical expected",C_result);

	long retCode;
	long objId;
	
	CCString curClass = LOCAL_FOREX_CLASS;
	CCString stringId;

	retCode = ARMLOCAL_FOREX (LocalGetNumObjectId (C_LeftCcy),
							  LocalGetNumObjectId (C_RightCcy),
                              spotValue,
							  C_result);

	if(sizeLeft == 3 && sizeRight == 3)
	{
		retCode = ARMLOCAL_FOREX (C_LeftCcy,
								  C_RightCcy,
								  spotValue,
								  C_result);
		
	}
	else
	{
		retCode = ARMLOCAL_FOREX (LocalGetNumObjectId (C_LeftCcy),
								  LocalGetNumObjectId (C_RightCcy),
								  spotValue,
								  C_result);
	}

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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_PXL_FOREX" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}

					  
// EOF %M% 
