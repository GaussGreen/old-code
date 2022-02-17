#pragma warning(disable : 4786)
#pragma warning(disable : 4541)

#include "XL_local_api.h"
#include <libCCxll\CCxll.h>

#include "ARM_local_interglob.h"
#include "ARM_local_interface.h"

#include <ARM\libarm_frometk\arm_local_parsexml.h>
#include <ARM\libarm_frometk\ARM_local_etoolkit.h>

#include <ARM\libarm_local\ARM_local_class.h>
#include <ARM\libarm_local\ARM_local_glob.h>
#include <ARM\libicm_local\ICM_local_summit.h>

#include <ARM\libarm_local\ARM_local_persistent.h>

#include <ARM\libarm\ARM_result.h>
#include <ARM\local_xlarm\ARM_local_interglob.h>
#include <ARM\local_xlarm\ARM_local_interface.h>
#include <ARM\local_xlarm\XL_local_xlarm_common.h>
#include <ARM\libarm_local\arm_local_init.h>
#include <ARM\libarm_local\ARM_local_volcrv.h>
#include <ARM\libarm_local\ARM_local_class.h>
#include "ExcelTools.h"
#include <stdlib.h>

/** 
__declspec(dllexport) LPXLOPER WINAPI ARM_Credit_GetModelFromSummit(LPXLOPER XL_idSummit,
																	LPXLOPER XL_type,
																	LPXLOPER XL_ircurve,
																	LPXLOPER XL_CurveId,
																	LPXLOPER XL_CorrCurveId)
{
	long retCode;
	ARM_result C_result;
static XLOPER XL_result;
	// return
	try {

	

	//	ARM_BEGIN();

	ARM_NOCALCIFWIZ();

	// C variable
	CCString C_idSummit;
	CCString C_type;
	CCString C_ircurve;

	CCString C_CurveId;
	CCString C_CurveId_default = "MO";

	CCString C_CorrCurveId;
	CCString C_CorrCurveId_default = "NONE";

	CCString typeId;

	// error
	static int error;
	static char* reason = "";

	XL_readStrCell(XL_idSummit,C_idSummit," ARM_ERR: summit id: string expected",C_result);
	XL_readStrCell(XL_type,C_type," ARM_ERR: type: string expected",C_result);
	XL_readStrCell(XL_ircurve,C_ircurve," ARM_ERR: ircurve: string expected",C_result);
	XL_readStrCellWD(XL_CurveId,C_CurveId,C_CurveId_default," ARM_ERR: CurveId: string expected",C_result);
	XL_readStrCellWD(XL_CorrCurveId,C_CorrCurveId,C_CorrCurveId_default," ARM_ERR: CorrCurveId: string expected",C_result);

	if((typeId = ARM_ConvTypeModel (C_type, C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	CCString prevClass;
	
	CCString curClass = typeId;
	CCString stringId = GetLastCurCellEnvValue ();
	
	long objId;

	if(!stringId)
	{
		retCode = ICMLOCAL_GetModelFromSummit (LocalGetNumObjectId(C_ircurve),
											   C_idSummit,
											   C_type,
											   C_CurveId,
											   C_CorrCurveId,
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
			retCode = ICMLOCAL_GetModelFromSummit (LocalGetNumObjectId(C_ircurve), 
													C_idSummit,
													C_type,
													C_CurveId,
													C_CorrCurveId,
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
			retCode = ICMLOCAL_GetModelFromSummit(LocalGetNumObjectId(C_ircurve), 
													C_idSummit,
													C_type,
													C_CurveId,
													C_CorrCurveId,
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
	catch (Exception&e)
	{
		ExcelCaller::get().setError(e.GetErrorString()); 
		ExcelTools::convert("ARM_ERR",&XL_result) ;
	}
	catch (std::exception&e)
	{
		ExcelCaller::get().setError(e.what()); 
		ExcelTools::convert("ARM_ERR",&XL_result) ;
	}
	catch(...)
	{
		ExcelCaller::get().setError("Unknown Exception");
		ExcelTools::convert("ARM_ERR",&XL_result) ;
	}
	
	return (LPXLOPER)&XL_result;
}

**/ 