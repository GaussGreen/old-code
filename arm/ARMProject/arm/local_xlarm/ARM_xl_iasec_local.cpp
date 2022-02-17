#pragma warning(disable : 4786)
#pragma warning(disable : 4005)

#include <libCCxll\CCxll.h>

#include "ARM_local_interglob.h"
#include "ARM_local_interface.h"

#include <ARM\libarm_local\ARM_local_class.h>
#include <ARM\libarm_local\ARM_local_iasec.h>

#include "ARM_xl_trycatch_local.h"
#include <util\fromto.h>

#include "util\tech_macro.h"

__declspec(dllexport) LPXLOPER WINAPI Local_IASEC (LPXLOPER XL_underlying,
												   LPXLOPER XL_iaCtrl,
												   LPXLOPER XL_refVal,
												   LPXLOPER XL_iaCtrlType)
{
	ADD_LOG("Local_IASEC ");
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();

	// C variable
	CCString C_underlying;
	CCString C_iaCtrl;
	CCString C_refVal;
	CCString C_iaCtrlType;
	long iaCtrlTypeId;

	// error
	static int error;
	static char* reason = "";

	XL_readStrCell(XL_underlying,C_underlying," ARM_ERR: underlying id: object expected",C_result);
	XL_readStrCell(XL_iaCtrl,C_iaCtrl," ARM_ERR: index amortizing control id: object expected",C_result);
	XL_readStrCell(XL_refVal,C_refVal," ARM_ERR: reference value id: object expected",C_result);
	XL_readStrCellWD(XL_iaCtrlType,C_iaCtrlType,"P"," ARM_ERR: index amortizing control type: string expected",C_result);
		
	if((iaCtrlTypeId = ARM_ConvPriceYield (C_iaCtrlType, C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	long retCode;
	long objId;
	CCString prevClass;
	
	CCString curClass = LOCAL_IASEC_CLASS;
	CCString stringId = GetLastCurCellEnvValue ();

	if(!stringId)
	{
		retCode = ARMLOCAL_IASEC (LocalGetNumObjectId (C_underlying),
								  LocalGetNumObjectId (C_iaCtrl),
								  LocalGetNumObjectId (C_refVal),
								  iaCtrlTypeId,
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
			retCode = ARMLOCAL_IASEC (LocalGetNumObjectId (C_underlying),
									  LocalGetNumObjectId (C_iaCtrl),
									  LocalGetNumObjectId (C_refVal),
									  iaCtrlTypeId,
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
			retCode = ARMLOCAL_IASEC (LocalGetNumObjectId (C_underlying),
									  LocalGetNumObjectId (C_iaCtrl),
									  LocalGetNumObjectId (C_refVal),
									  iaCtrlTypeId,
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_IASEC" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}

__declspec(dllexport) LPXLOPER WINAPI Local_PXL_IASEC (LPXLOPER XL_underlying,
													   LPXLOPER XL_iaCtrl,
													   LPXLOPER XL_refVal,
													   LPXLOPER XL_iaCtrlType)
{
	ADD_LOG("Local_PXL_IASEC ");
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();

	// C variable
	CCString C_underlying;
	CCString C_iaCtrl;
	CCString C_refVal;
	CCString C_iaCtrlType;
	long iaCtrlTypeId;

	// error
	static int error;
	static char* reason = "";

	XL_readStrCell(XL_underlying,C_underlying," ARM_ERR: underlying id: object expected",C_result);
	XL_readStrCell(XL_iaCtrl,C_iaCtrl," ARM_ERR: index amortizing control id: object expected",C_result);
	XL_readStrCell(XL_refVal,C_refVal," ARM_ERR: reference value id: object expected",C_result);
	XL_readStrCellWD(XL_iaCtrlType,C_iaCtrlType,"P"," ARM_ERR: index amortizing control type: string expected",C_result);
		
	if((iaCtrlTypeId = ARM_ConvPriceYield (C_iaCtrlType, C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	long retCode;
	long objId;
	
	CCString curClass = LOCAL_IASEC_CLASS;
	CCString stringId;

	retCode = ARMLOCAL_IASEC (LocalGetNumObjectId (C_underlying),
							  LocalGetNumObjectId (C_iaCtrl),
							  LocalGetNumObjectId (C_refVal),
							  iaCtrlTypeId,
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_PXL_IASEC" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}



// EOF %M% 
