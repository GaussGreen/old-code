#include <ARM\libarm_local\firstToBeIncluded.h>

#include <libCCxll\CCxll.h>
#include <ARM\libarm_local\ARM_local_glob.h>
#include <ARM\libarm_local\ARM_local_class.h>
#include <ARM\libarm_frometk\ARM_local_etoolkit.h>
#include <ARM\libarm_frometk\ARM_local_etoolkitX.h>
#include <ARM\libarm_frometk\arm_local_paesexml_calypso.h> 
//#include <ARM\libarm_frometk\ARM_local_wsetk.h>

#include <ARM\libarm\ARM_result.h>

#include "ARM_local_interglob.h"
#include "ARM_xl_gp_fctorhelper.h"
#include "ARM_local_interface.h"
#include "XL_local_xlarm_common.h"
#include "ARM_xl_glob_local.h"
#include "ARM_local_help.h"
#include "ARM_xl_wrapper_local.h"

/// general macro for try catch
#include "ARM_xl_trycatch_local.h"
#include <util\fromto.h>
#include "ExcelTools.h"

#include "util\tech_macro.h"

static HGLOBAL Mmglob = NULL;

_declspec(dllexport) LPXLOPER WINAPI Local_bsflexible (LPXLOPER XL_F,LPXLOPER XL_V,
														LPXLOPER XL_B,LPXLOPER XL_K,
													   LPXLOPER XL_CallPutFlag)
{
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();

	// C variable
	double C_F;
	double C_V;
	double C_B;
	double C_K;
	double C_CallPutFlag;

	// error
	static int error;
	static char* reason = "";

	XL_readNumCell(XL_F,C_F," ARM_ERR: forward: numeric expected",C_result);
	XL_readNumCell(XL_V,C_V," ARM_ERR: volatility: numeric expected",C_result);
	XL_readNumCell(XL_B,C_B," ARM_ERR: bond price: numeric expected",C_result);
	XL_readNumCell(XL_K,C_K," ARM_ERR: strike: numeric expected",C_result);
	XL_readNumCell(XL_CallPutFlag,C_CallPutFlag," ARM_ERR: CallPutFlag: numeric expected",C_result);
	
	long retCode = ARMLOCAL_bsflexible (
		C_F,C_V,C_B,C_K,C_CallPutFlag,C_result
			);

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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_bsflexible" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}



/// we had a final try catch block to make Local_ARM_Price robust!
__declspec(dllexport) LPXLOPER WINAPI Local_ARM_Price (LPXLOPER XL_secId,
													   LPXLOPER XL_modId)
{
	ADD_LOG("Local_ARM_Price ");
	///	ARM_BEGIN();
	// return
	static XLOPER XL_result;
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
		ARM_NOCALCIFWIZ();

		// C variable
		CCString C_secId;
		CCString C_modId;
		long modId;
		
		// error
		static int error;
		static char* reason = "";

		XL_readStrCell(XL_secId,C_secId," ARM_ERR: security id: object expected",C_result);
		XL_readStrCellWD(XL_modId,C_modId,"DEFAULT"," ARM_ERR: model id: object expected",C_result);

		if(C_modId == "DEFAULT")
		{
			modId = ARM_NULL_OBJECT;
		}
		else
		{
			modId = LocalGetNumObjectId (C_modId);
		}

		long retCode = ARMLOCAL_ARM_Price (LocalGetNumObjectId (C_secId), modId, C_result);

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

	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_ARM_Price" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}



__declspec(dllexport) LPXLOPER WINAPI Local_ARM_Cover (LPXLOPER XL_secId,
													   LPXLOPER XL_modId)
{
	ADD_LOG("Local_ARM_Cover ");
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	LPXLOPER pxArray;
	ARM_result C_result;

	// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();

	// C variable
	CCString C_secId;
	CCString C_modId;
	
	// error
	static int error;
	static char* reason = "";

	XL_readStrCell(XL_secId,C_secId," ARM_ERR: security id: object expected",C_result);
	XL_readStrCell(XL_modId,C_modId," ARM_ERR: model id: object expected",C_result);
	
	long retCode = ARMLOCAL_ARM_Cover (LocalGetNumObjectId (C_secId), LocalGetNumObjectId (C_modId) , C_result);

	if(retCode == ARM_OK)
	{
			long nbrows = C_result.getLong();
			int nbcolumns = 1;
		
			FreeCurCellErr ();

			XL_result.xltype = xltypeMulti;
			XL_result.val.array.columns = nbcolumns;
			XL_result.val.array.rows = nbrows; 

			if (Mmglob)
				GlobalFree(Mmglob);

			Mmglob = GlobalAlloc (GMEM_ZEROINIT, nbrows * nbcolumns * sizeof (XLOPER));
			XL_result.val.array.lparray = pxArray = (LPXLOPER)(Mmglob);
           
            
			for (int i = 0; i < nbrows; i++)
			{	
                pxArray[XL_Coordonnate2Rank (i, 0, nbcolumns)].xltype = xltypeNum;
				pxArray[XL_Coordonnate2Rank (i, 0, nbcolumns)].val.num = C_result.getArray(i);          
			}
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_ARM_Cover" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


__declspec(dllexport) LPXLOPER WINAPI Local_FreeObject (LPXLOPER XL_secId)
{
	ADD_LOG("Local_FreeObject ");
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	// to make the infrastructure very robust we put a try catch!
    ARM_XL_TRY_BLOCK_BEGIN
	{
	    ARM_NOCALCIFWIZ();

	    // C variable
	    CCString C_secId;

	    // error
	    static int error;
	    static char* reason = "";

	    XL_readStrCell(XL_secId,C_secId," ARM_ERR: security id: object expected",C_result);

	    long retCode = ARMLOCAL_FreeObject (LocalGetNumObjectId (C_secId), C_result);

	    if(retCode == ARM_OK)
	    {
		    FreeCurCellErr ();
		    XL_result.xltype = xltypeNum;
		    XL_result.val.num = C_result.getLong ();
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_FreeObject" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


__declspec(dllexport) LPXLOPER WINAPI Local_FreeAllObjects ()
{
	ADD_LOG("Local_FreeAllObjects ");
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();

	// error
	static int error;
	static char* reason = "";

	long retCode = ARMLOCAL_FreeAllObjects (C_result);

	LOCALARM_PersistentListsClear ();

	if(retCode == ARM_OK)
	{
		FreeCurCellErr ();
		XL_result.xltype = xltypeNum;
		XL_result.val.num = C_result.getLong ();
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_FreeAllObjects" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}



__declspec(dllexport) LPXLOPER WINAPI Local_NextBusinessDay (LPXLOPER XL_date,
															 LPXLOPER XL_currency,
															 LPXLOPER XL_days)
{
	ADD_LOG("Local_NextBusinessDay ");
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();

	// C variable
	double C_date;
	CCString C_currency;
	double C_days;
	double C_days_default = 1;
	
	// error
	static int error;
	static char* reason = "";

	XL_readNumCell(XL_date,C_date," ARM_ERR: date: date expected",C_result);
	XL_readStrCellWD(XL_currency,C_currency,"DEFAULT"," ARM_ERR: currency: string expected",C_result);
	XL_readNumCellWD(XL_days,C_days,C_days_default," ARM_ERR: number of days: numeric expected",C_result);

	if(C_currency == "DEFAULT")
	{
		ARM_result currencyres;
		ARMLOCAL_ARM_GetDefaultCurrency (currencyres);
		if(currencyres.getRetCode () != ARM_OK)
		{
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}
		else
		{
			C_currency = currencyres.getString ();
		}
	}

	long retCode = ARMLOCAL_NextBusinessDay (C_date, C_currency, (long)C_days, C_result);

	if(retCode == ARM_OK)
	{
		FreeCurCellErr ();
		XL_result.xltype = xltypeNum;
		XL_result.val.num = Local_ARMDATE2XLDATE(C_result.getString ());
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_NextBusinessDay" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}



__declspec(dllexport) LPXLOPER WINAPI Local_IsBusinessDay (LPXLOPER XL_date,
														   LPXLOPER XL_currency)
{
	ADD_LOG("Local_IsBusinessDay ");
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();

	// C variable
	double C_date;
	CCString C_currency;
	
	// error
	static int error;
	static char* reason = "";

	XL_readNumCell(XL_date,C_date," ARM_ERR: date: date expected",C_result);
	XL_readStrCellWD(XL_currency,C_currency,"DEFAULT"," ARM_ERR: currency: string expected",C_result);
	
	if(C_currency == "DEFAULT")
	{
		ARM_result currencyres;
		ARMLOCAL_ARM_GetDefaultCurrency (currencyres);
		if(currencyres.getRetCode () != ARM_OK)
		{
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}
		else
		{
			C_currency = currencyres.getString ();
		}
	}
		
	long retCode = ARMLOCAL_IsBusinessDay (C_date, C_currency, C_result);

	if(retCode == ARM_OK)
	{
		FreeCurCellErr ();
		XL_result.xltype = xltypeNum;
		XL_result.val.num = C_result.getLong ();
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_IsBusinessDay" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}



__declspec(dllexport) LPXLOPER WINAPI Local_ARM_FutDelivery (LPXLOPER XL_Fut,
															 LPXLOPER XL_currency)
{
	ADD_LOG("Local_ARM_FutDelivery ");
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();

	// C variable
	CCString C_fut;
	CCString C_currency;
	
	// error
	static int error;
	static char* reason = "";

	XL_readStrCell(XL_Fut,C_fut," ARM_ERR: Future: string expected",C_result);
	XL_readStrCell(XL_currency,C_currency," ARM_ERR: currency: string expected",C_result);

	long retCode = ARMLOCAL_FutDelivery (C_fut, C_currency, C_result);

	if(retCode == ARM_OK)
	{
		FreeCurCellErr ();
		XL_result.xltype = xltypeNum;
		XL_result.val.num = Local_ARMDATE2XLDATE(C_result.getString ());
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_ARM_FutDelivery" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


__declspec(dllexport) LPXLOPER WINAPI Local_ARM_Accrued (LPXLOPER XL_secId,
														 LPXLOPER XL_fwdDate,
														 LPXLOPER XL_modId)
{
	ADD_LOG("Local_ARM_Accrued ");
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();

	// C variable
	CCString C_secId;

	double C_fwdDate;

	CCString C_modId;
	long modId;
	
	// error
	static int error;
	static char* reason = "";

	XL_readStrCell(XL_secId,C_secId," ARM_ERR: security id: object expected",C_result);
	XL_readNumCell(XL_fwdDate,C_fwdDate," ARM_ERR: forward date: date expected",C_result);
	XL_readStrCellWD(XL_modId,C_modId,"DEFAULT"," ARM_ERR: model id: object expected",C_result);
	
	if(C_modId == "DEFAULT")
	{
		modId = ARM_NULL_OBJECT;
	}
	else
	{
		modId = LocalGetNumObjectId (C_modId);
	}

	long retCode = ARMLOCAL_ARM_Accrued (LocalGetNumObjectId (C_secId), C_fwdDate,
							        modId, C_result);

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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_ARM_Accrued" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


__declspec(dllexport) LPXLOPER WINAPI Local_IMPLIEDVOL (LPXLOPER XL_instId,
														LPXLOPER XL_modId,
														LPXLOPER XL_price,
														LPXLOPER XL_LnOrNorVol)
{
	ADD_LOG("Local_IMPLIEDVOL ");
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();

	// C variable
	CCString C_instId;
	CCString C_modId;
	double C_price;
	
	// error
	static int error;
	static char* reason = "";

	XL_readStrCell(XL_instId,C_instId," ARM_ERR: instrument id: object expected",C_result);
	XL_readStrCell(XL_modId,C_modId," ARM_ERR: model id: object expected",C_result);
	XL_readNumCell(XL_price,C_price," ARM_ERR: price: numeric expected",C_result);

	CCString lnOrNorVol;
	XL_readStrCellWD(XL_LnOrNorVol,lnOrNorVol,"Y"," ARM_ERR: LogNor or Nor Vol: array of string expected (Y/N)",C_result);
	bool isLnVol=true;
    lnOrNorVol.toUpper();
    if(CCSTringToSTLString(lnOrNorVol)!="Y")
            isLnVol=false;
	
	long retCode = ARMLOCAL_IMPLIEDVOL (LocalGetNumObjectId (C_instId), 
									LocalGetNumObjectId (C_modId),
									C_price,isLnVol, C_result);

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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_IMPLIEDVOL" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


long LocalDisplayArmFile (const CCString& fileName)
{
	CCString	vEditor(EDITOR);
	int	vLen = fileName.GetLen();
	if(vLen > 3)
	{
		char*	vExtension = fileName.c_str() + vLen - 3;
		CCString	cExtension(vExtension);
		if(cExtension == "xml")
		{
			vEditor = XML_EDITOR;	
		}
	}

	/// the start command helps to start on a new thread the editor!
	CCString command = CCString("start ") + vEditor + " " + fileName;

	//_flushall();

	system((const char*) command);

	return(ARM_OK);
}


long LocalGetArmViewFile (const CCString& sockId)
{
    CCString clientViewFileName = CCString(VIEW_FILE_CLIENT_LOCATION)+CCString(VIEW_FILE_PREFIX)+sockId;

	LocalDisplayArmFile (clientViewFileName);

	/// Watch out! because we start in a new thread
	/// we cannot delete the file anymore!!!!
	/// otherwise we will delete while it is processed

	return ARM_OK;
}

long Local_ViewFile (const CCString& C_instId)
{
	ARM_result C_result;
	
	long retCode = ARMLOCAL_ARM_View (LocalGetNumObjectId (C_instId), C_result);

	if(retCode == ARM_OK)
	{
		char username[100];
		DWORD nbChar = sizeof(username);

		GetUserName(username,&nbChar);

		retCode = LocalGetArmViewFile ((CCString)"123" + (CCString)username);
	}

	return retCode;
}



__declspec(dllexport) LPXLOPER WINAPI Local_ARM_View (LPXLOPER XL_instId)
{
	ADD_LOG("Local_ARM_View ");
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();

	// C variable
	CCString C_instId;

	// error
	static int error;
	static char* reason = "";

	XL_readStrCell(XL_instId,C_instId," ARM_ERR: instrument id: object expected",C_result);

	long retCode = Local_ViewFile (C_instId);

	if(retCode == ARM_OK)
	{
		FreeCurCellErr ();
		XL_result.xltype = xltypeStr;
		XL_result.val.str = "\006ARM_OK";
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_ARM_View" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}

__declspec(dllexport) long WINAPI Local_View ()
{
	ARM_result C_result;
	
	static XLOPER XL_instId;

	// error
	static int error;
	static char* reason = "";

	CCString C_instId;
	
	XL_getActiveCellContent (&XL_instId);
	XL_getStrCell (&XL_instId, &error, &reason, XL_ERROR_VALUE_MISSING, XL_ERROR_VALUE_INCONSISTENCY, XL_ERROR_VALUE_IS_ERR, C_instId);

	long retCode = ARM_KO;

	if(error == XL_NO_ERROR)
	{
		retCode = ARMLOCAL_ARM_View (LocalGetNumObjectId (C_instId), C_result);

		if(retCode == ARM_OK)
		{
			char username[100];
			DWORD nbChar = sizeof(username);

			GetUserName(username,&nbChar);

			CCString	vViewFile( (CCString)"123" + (CCString)username );


//			if( LocalGetStringObjectClassAndError(C_instId) == "LMRCR" )
//				vViewFile += ".xml";

			retCode = LocalGetArmViewFile(vViewFile);
			if(retCode != ARM_OK)
			{
				Excel(xlcAlert, 0, 2, TempStr (" ARM_ERR: Can't retrieve the view file"), TempInt(2));
			}
		}
		else
		{
			CCString msg = CCString ("  ") + C_result.getMsg ();
			Excel(xlcAlert, 0, 2, TempStr ((char*)(const char*)msg), TempInt(2));
		}
	}
	else
	{
		CCString msg = CCString (" ARM_ERR: object id: string expected");
		Excel(xlcAlert, 0, 2, TempStr ((char*)(const char*)msg), TempInt(2));
	}
	
	return retCode;
}


__declspec(dllexport) long WINAPI Local_View_XML()
{
	ARM_result C_result;
	
	static XLOPER XL_instId;

	// error
	static int error;
	static char* reason = "";

	CCString C_instId;
	
	XL_getActiveCellContent (&XL_instId);
	XL_getStrCell (&XL_instId, &error, &reason, XL_ERROR_VALUE_MISSING, XL_ERROR_VALUE_INCONSISTENCY, XL_ERROR_VALUE_IS_ERR, C_instId);

	long retCode = ARM_KO;

	if( LocalGetStringObjectClassAndError(C_instId) == "LMRCR" )
	{
		if(error == XL_NO_ERROR)
		{
			retCode = ARMLOCAL_ARM_View (LocalGetNumObjectId (C_instId), C_result, true);

			if(retCode == ARM_OK)
			{
				char username[100];
				DWORD nbChar = sizeof(username);

				GetUserName(username,&nbChar);

				CCString	vViewFile( (CCString)"123" + (CCString)username + (CCString)".xml" );

				retCode = LocalGetArmViewFile(vViewFile);
				if(retCode != ARM_OK)
				{
					Excel(xlcAlert, 0, 2, TempStr (" ARM_ERR: Can't retrieve the view file"), TempInt(2));
				}
			}
			else
			{
				CCString msg = CCString ("  ") + C_result.getMsg ();
				Excel(xlcAlert, 0, 2, TempStr ((char*)(const char*)msg), TempInt(2));
			}
		}
		else
		{
			CCString msg = CCString (" ARM_ERR: object id: string expected");
			Excel(xlcAlert, 0, 2, TempStr ((char*)(const char*)msg), TempInt(2));
		}
	}
	else
	{
		CCString msg = CCString (" ARM_ERR: XML View available only for Mercure results");
		Excel(xlcAlert, 0, 2, TempStr ((char*)(const char*)msg), TempInt(2));
	}
	
	return retCode;
}


__declspec(dllexport) int WINAPI Local_DisplayErrorMessage ()
{
	CCString caller = XL_getActiveCell ();
	CCString errorMessage = GetErrValue (caller);

	/// Warning we need the two blank before the error message!

    char errMsgBuf[1024];

    strcpy(errMsgBuf, "  ");

    strncat(errMsgBuf, (const char *) errorMessage,1024-2);
	errMsgBuf[1024-1]=0; 

	// Excel(xlcAlert, 0, 2, TempStr ((char *) errMsgBuf), TempInt(2));
	::MessageBox(NULL,errorMessage.c_str(),"ARM Error",MB_OK|MB_ICONINFORMATION) ; 

	return(0);
}


__declspec(dllexport) LPXLOPER WINAPI Local_SetNotional (LPXLOPER XL_secId,
														 LPXLOPER XL_refValId,
														 LPXLOPER XL_percentRemainder,
														 LPXLOPER XL_datesForInterpol)
{
	ADD_LOG("Local_SetNotional ");
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();

	// C variable
	CCString C_secId;
	CCString C_refValId;

	double C_percentRemainder;
	double C_percentRemainder_default = 100.0;

	double C_datesForInterpol;
	double C_datesForInterpol_default = 3.0;

	// error
	static int error;
	static char* reason = "";

	XL_readStrCell(XL_secId,C_secId," ARM_ERR: security id: object expected",C_result);
	XL_readStrCell(XL_refValId,C_refValId," ARM_ERR: reference value id: object expected",C_result);
	XL_readNumCellWD(XL_percentRemainder,C_percentRemainder,C_percentRemainder_default," ARM_ERR: percentage of remainder: numeric expected",C_result);
	XL_readNumCellWD(XL_datesForInterpol,C_datesForInterpol,C_datesForInterpol_default," ARM_ERR: dates on which to interpolate: numeric expected",C_result);

	long retCode = ARMLOCAL_SetNotional (LocalGetNumObjectId (C_secId),
										 LocalGetNumObjectId (C_refValId),
										 C_percentRemainder,
										 C_result,
										 (long)C_datesForInterpol);

	if ( retCode == ARM_OK )
	{
		FreeCurCellErr();
	
        XL_result.xltype = xltypeStr;
		XL_result.val.str = XL_StrC2StrPascal(C_secId);
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_SetNotional" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


__declspec(dllexport) LPXLOPER WINAPI Local_ARM_today ()
{
	ADD_LOG("Local_ARM_today ");
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();

	// error
	static int error;
	static char* reason = "";

	int y, m, d, hh, mm, ss;

	DAT_gmt_to_struct (DAT_now, &y, &m, &d, &hh, &mm, &ss);
	
	XL_result.xltype = xltypeNum;
	XL_result.val.num =  DAT_struct_to_ssdate (y, m ,d);
	
//	ARM_END();
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_ARM_today" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


__declspec(dllexport) LPXLOPER WINAPI Local_GetExpiry (LPXLOPER XL_secId)
{
	ADD_LOG("Local_GetExpiry ");
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();

	// C variable
	CCString C_secId;
	
	// error
	static int error;
	static char* reason = "";

	XL_readStrCell(XL_secId,C_secId," ARM_ERR: security id: object expected",C_result);
	
	long retCode = ARMLOCAL_GetExpiry (LocalGetNumObjectId (C_secId), C_result);

	if(retCode == ARM_OK)
	{
		FreeCurCellErr ();
		XL_result.xltype = xltypeNum;
		XL_result.val.num = Local_ARMDATE2XLDATE(C_result.getString ());
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_GetExpiry" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


__declspec(dllexport) LPXLOPER WINAPI Local_Sensitivity (LPXLOPER XL_secId,
														 LPXLOPER XL_modId,
														 LPXLOPER XL_param)
{
	ADD_LOG("Local_Sensitivity ");
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
		ARM_NOCALCIFWIZ();

		// C variable
		CCString C_secId;
		CCString C_modId;

		CCString C_param;
		long paramId;

		// error
		static int error;
		static char* reason = "";

		XL_readStrCell(XL_secId,C_secId," ARM_ERR: security id: object expected",C_result);
		XL_readStrCell(XL_modId,C_modId," ARM_ERR: model id: object expected",C_result);
		XL_readStrCell(XL_param,C_param," ARM_ERR: param: string expected",C_result);
		
		if((paramId = ARM_ConvSvtyParam (C_param, C_result)) == ARM_DEFAULT_ERR)
		{
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}
			
		long retCode = ARMLOCAL_Sensitivity (LocalGetNumObjectId (C_secId), 
										LocalGetNumObjectId (C_modId),
										paramId, C_result);

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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_Sensitivity" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


__declspec(dllexport) LPXLOPER WINAPI Local_GetFRMShortRateVols(LPXLOPER XL_modelId)
{
	ADD_LOG("Local_GetFRMShortRateVols");
													      
//	ARM_BEGIN();
	
	// return
	static XLOPER XL_result;
	LPXLOPER pxArray;
	ARM_result C_result;

	// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();

	// C variable	
	CCString C_modelId;
	
	VECTOR<double> matu;
	VECTOR<double> rate;
	
	// error
	static int error;
	static char* reason = "";

	XL_readStrCell(XL_modelId,C_modelId," ARM_ERR: Model ID: String expected",C_result);
	
	long retCode;
	
	retCode = ARMLOCAL_GetFRMShortRateVols(LocalGetNumObjectId(C_modelId),&matu,&rate,
								      C_result);
	if ( retCode == ARM_OK )
	{
		int nbrows = matu.size ();
		int nbcolumns = 2;
	
		FreeCurCellErr ();

		XL_result.xltype = xltypeMulti;
		XL_result.val.array.columns = nbcolumns;
		XL_result.val.array.rows = nbrows;
		XL_result.val.array.lparray = pxArray = (LPXLOPER)GlobalAlloc (GMEM_ZEROINIT, nbrows * nbcolumns * sizeof (XLOPER));

		for (int i = 0; i < nbrows; i++)
		{
			pxArray[XL_Coordonnate2Rank(i, 0, nbcolumns)].xltype = xltypeNum;
			pxArray[XL_Coordonnate2Rank(i, 0, nbcolumns)].val.num = matu[i];
			pxArray[XL_Coordonnate2Rank(i, 1, nbcolumns)].xltype = xltypeNum;
			pxArray[XL_Coordonnate2Rank(i, 1, nbcolumns)].val.num = rate[i]; 
		}
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_GetFRMShortRateVols" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}					  


__declspec(dllexport) LPXLOPER WINAPI Local_XCcyAdjust (LPXLOPER XL_startDate,
														LPXLOPER XL_endDate, 
														LPXLOPER XL_payFreq,
														LPXLOPER XL_domCcy,
														LPXLOPER XL_forIndexType,
														LPXLOPER XL_forCcy,
														LPXLOPER XL_spreadsId,
														LPXLOPER XL_zcDomId,
														LPXLOPER XL_discDomId,
														LPXLOPER XL_zcForId,
														LPXLOPER XL_discForId,
														LPXLOPER XL_FX,
														LPXLOPER XL_CouponId,
														LPXLOPER XL_domDc,
														LPXLOPER XL_forDc)
{
	ADD_LOG("Local_XCcyAdjust ");
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();

	// C variable
	CCString C_spreadsId;
    CCString C_zcDomId;
    CCString C_discDomId;
    CCString C_zcForId;
    CCString C_discForId;
    CCString C_domCcy;
    CCString C_forCcy;
    CCString C_forIndexType;
    CCString C_couponId;
    CCString C_domDc;
    CCString C_forDc;

	double C_FX;
    double C_startDate;
    double C_endDate;
    double C_payFreq;

    int domDc;
    int forDc;
    int couponId;
	long forIndexType;
    
	// error
	static int error;
	static char* reason = "";

    XL_readNumCell(XL_startDate,C_startDate," ARM_ERR: Start date: date expected",C_result);
    XL_readNumCell(XL_endDate,C_endDate," ARM_ERR: End date: date expected",C_result);
    XL_readNumCell(XL_payFreq,C_payFreq," ARM_ERR: PayFreq: numeric expected",C_result);
	XL_readStrCell(XL_domCcy,C_domCcy," ARM_ERR: Domestic Currency: string expected",C_result);
	XL_readStrCell(XL_forIndexType,C_forIndexType," ARM_ERR: Foreign index type : string expected",C_result);
	XL_readStrCell(XL_forCcy,C_forCcy," ARM_ERR: Foreign Currency: string expected",C_result);
   	XL_readStrCell(XL_spreadsId, C_spreadsId," ARM_ERR: Spread Curve : object expected",C_result);
   	XL_readStrCell(XL_zcDomId, C_zcDomId," ARM_ERR: Domestic zero curve: object expected",C_result);
   	XL_readStrCell(XL_discDomId, C_discDomId," ARM_ERR: Domestic discount Curve : object expected",C_result);
   	XL_readStrCell(XL_zcForId, C_zcForId," ARM_ERR: Foreign zero Curve : object expected",C_result);
   	XL_readStrCell(XL_discForId, C_discForId," ARM_ERR: Foreign discount Curve : object expected",C_result);
    XL_readNumCell(XL_FX,C_FX," ARM_ERR: FX : numeric expected",C_result);
    XL_readStrCellWD(XL_CouponId, C_couponId,"DEFAULT"," ARM_ERR: coupon id: object expected",C_result);
    XL_readStrCellWD(XL_domDc, C_domDc, "A360"," ARM_ERR: domestic day count: string expected",C_result);
    XL_readStrCellWD(XL_forDc, C_forDc, "A360"," ARM_ERR: foreign day count: string expected",C_result);

	if(C_couponId == "DEFAULT")
    {
        couponId = ARM_NULL_OBJECT;
	}
	else
	{
		couponId = LocalGetNumObjectId (C_couponId);
	}

    domDc = ARM_ConvDayCount(C_domDc);

    forDc = ARM_ConvDayCount(C_forDc);

	if ((forIndexType = ARM_ConvIrIndName (C_forIndexType, C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	long retCode;
	long objId;
	CCString prevClass;
	
	CCString curClass = LOCAL_REFVAL_CLASS;             
	CCString stringId = GetLastCurCellEnvValue ();
	
	if(!stringId)
	{
		retCode = ARMLOCAL_XCccyAdjustment(C_startDate, C_endDate, C_payFreq,C_domCcy, forIndexType,
                                      C_forCcy, LocalGetNumObjectId(C_spreadsId), LocalGetNumObjectId(C_zcDomId),
                                      LocalGetNumObjectId(C_discDomId), LocalGetNumObjectId(C_zcForId), LocalGetNumObjectId(C_discForId),
                                      C_FX, couponId, domDc, forDc, C_result);

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
			retCode = ARMLOCAL_XCccyAdjustment(C_startDate, C_endDate, C_payFreq,C_domCcy, forIndexType,
                                      C_forCcy, LocalGetNumObjectId(C_spreadsId), LocalGetNumObjectId(C_zcDomId),
                                      LocalGetNumObjectId(C_discDomId), LocalGetNumObjectId(C_zcForId), LocalGetNumObjectId(C_discForId),
                                      C_FX, couponId, domDc, forDc,C_result, objId);

			if(retCode == ARM_OK)
			{
				LocalSetCurCellEnvValue (curClass, objId); 

				stringId = LocalMakeObjectId (objId, curClass);
			}
		}
		else
		{
			FreeCurCellContent ();
			retCode = ARMLOCAL_XCccyAdjustment(C_startDate, C_endDate, C_payFreq,C_domCcy, forIndexType,
                                      C_forCcy, LocalGetNumObjectId(C_spreadsId), LocalGetNumObjectId(C_zcDomId),
                                      LocalGetNumObjectId(C_discDomId), LocalGetNumObjectId(C_zcForId), LocalGetNumObjectId(C_discForId),
                                      C_FX, couponId, domDc, forDc,C_result);
		
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_XCcyAdjust" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


__declspec(dllexport) LPXLOPER WINAPI Local_PXL_XCcyAdjust (LPXLOPER XL_startDate,
															LPXLOPER XL_endDate, 
															LPXLOPER XL_payFreq,
															LPXLOPER XL_domCcy,
															LPXLOPER XL_forIndexType,
															LPXLOPER XL_forCcy,
															LPXLOPER XL_spreadsId,
															LPXLOPER XL_zcDomId,
															LPXLOPER XL_discDomId,
															LPXLOPER XL_zcForId,
															LPXLOPER XL_discForId,
															LPXLOPER XL_FX,
															LPXLOPER XL_CouponId,
															LPXLOPER XL_domDc,
															LPXLOPER XL_forDc)
{
	ADD_LOG("Local_PXL_XCcyAdjust ");
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();

	// C variable
	CCString C_spreadsId;
    CCString C_zcDomId;
    CCString C_discDomId;
    CCString C_zcForId;
    CCString C_discForId;
    CCString C_domCcy;
    CCString C_forCcy;
    CCString C_forIndexType;
    CCString C_couponId;
    CCString C_domDc;
    CCString C_forDc;

	double C_FX;
    double C_startDate;
    double C_endDate;
    double C_payFreq;

    int domDc;
    int forDc;
    int couponId;
	long forIndexType;
    
	// error
	static int error;
	static char* reason = "";

    XL_readNumCell(XL_startDate,C_startDate," ARM_ERR: Start date: date expected",C_result);
    XL_readNumCell(XL_endDate,C_endDate," ARM_ERR: End date: date expected",C_result);
    XL_readNumCell(XL_payFreq,C_payFreq," ARM_ERR: PayFreq: numeric expected",C_result);
	XL_readStrCell(XL_domCcy,C_domCcy," ARM_ERR: Domestic Currency: string expected",C_result);
	XL_readStrCell(XL_forIndexType,C_forIndexType," ARM_ERR: Foreign index type : string expected",C_result);
	XL_readStrCell(XL_forCcy,C_forCcy," ARM_ERR: Foreign Currency: string expected",C_result);
   	XL_readStrCell(XL_spreadsId, C_spreadsId," ARM_ERR: Spread Curve : object expected",C_result);
   	XL_readStrCell(XL_zcDomId, C_zcDomId," ARM_ERR: Domestic zero curve: object expected",C_result);
   	XL_readStrCell(XL_discDomId, C_discDomId," ARM_ERR: Domestic discount Curve : object expected",C_result);
   	XL_readStrCell(XL_zcForId, C_zcForId," ARM_ERR: Foreign zero Curve : object expected",C_result);
   	XL_readStrCell(XL_discForId, C_discForId," ARM_ERR: Foreign discount Curve : object expected",C_result);
    XL_readNumCell(XL_FX,C_FX," ARM_ERR: FX : numeric expected",C_result);
    XL_readStrCellWD(XL_CouponId, C_couponId,"DEFAULT"," ARM_ERR: coupon id: object expected",C_result);
    XL_readStrCellWD(XL_domDc, C_domDc, "A360"," ARM_ERR: domestic day count: string expected",C_result);
    XL_readStrCellWD(XL_forDc, C_forDc, "A360"," ARM_ERR: foreign day count: string expected",C_result);

	if(C_couponId == "DEFAULT")
    {
        couponId = ARM_NULL_OBJECT;
	}
	else
	{
		couponId = LocalGetNumObjectId (C_couponId);
	}


    domDc = ARM_ConvDayCount(C_domDc);

    forDc = ARM_ConvDayCount(C_forDc);

	if ((forIndexType = ARM_ConvIrIndName (C_forIndexType, C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	long retCode;
	long objId;
	
	CCString curClass = LOCAL_REFVAL_CLASS;             
	CCString stringId;

	retCode = ARMLOCAL_XCccyAdjustment(C_startDate, C_endDate, C_payFreq,C_domCcy, forIndexType,
                                  C_forCcy, LocalGetNumObjectId(C_spreadsId), LocalGetNumObjectId(C_zcDomId),
                                  LocalGetNumObjectId(C_discDomId), LocalGetNumObjectId(C_zcForId), LocalGetNumObjectId(C_discForId),
                                  C_FX, couponId, domDc, forDc, C_result);

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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_PXL_XCcyAdjust" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}



__declspec(dllexport) LPXLOPER WINAPI Local_ADJUSTTOBUSDATE (LPXLOPER XL_date,
															 LPXLOPER XL_currency,
															 LPXLOPER XL_rule)
{
	ADD_LOG("Local_ADJUSTTOBUSDATE ");
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();

	// C variable
	double C_date;
	CCString C_currency;

	CCString C_rule;
	long ruleId;

	// error
	static int error;
	static char* reason = "";

	XL_readNumCell(XL_date,C_date," ARM_ERR: date: date expected",C_result);
	XL_readStrCell(XL_currency,C_currency," ARM_ERR: currency: string expected",C_result);
	XL_readStrCell(XL_rule,C_rule," ARM_ERR: rule: string expected",C_result);

	if((ruleId = ARM_ConvRule (C_rule, C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}
		
	long retCode = ARMLOCAL_ADJUSTTOBUSDATE (C_date, C_currency, ruleId, C_result);

	if(retCode == ARM_OK)
	{
		FreeCurCellErr ();
		XL_result.xltype = xltypeNum;
		XL_result.val.num = Local_ARMDATE2XLDATE(C_result.getString ());
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_ADJUSTTOBUSDATE" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}



__declspec(dllexport) LPXLOPER WINAPI Local_FwdPrice (LPXLOPER XL_secId,
													  LPXLOPER XL_modId,
													  LPXLOPER XL_fwdDate)
{
	ADD_LOG("Local_FwdPrice ");
//	ARM_BEGIN();

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
	double C_fwdDate;
	
	// error
	static int error;
	static char* reason = "";

	XL_readStrCell(XL_secId,C_secId," ARM_ERR: security id: object expected",C_result);
	XL_readStrCell(XL_modId,C_modId," ARM_ERR: model id: object expected",C_result);
	XL_readNumCell(XL_fwdDate,C_fwdDate," ARM_ERR: forward date: date expected",C_result);
	
	long retCode = ARMLOCAL_FwdPrice (LocalGetNumObjectId (C_secId), LocalGetNumObjectId (C_modId), C_fwdDate, C_result);

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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_FwdPrice" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


__declspec(dllexport) LPXLOPER WINAPI Local_CvSensitivity (LPXLOPER XL_secId,
														   LPXLOPER XL_modId,
														   LPXLOPER XL_param)
{
	ADD_LOG("Local_CvSensitivity ");
//	ARM_BEGIN();

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

	CCString C_param;
	long paramId;

	CCString id ("");

	// error
	static int error;
	static char* reason = "";

	XL_readStrCell(XL_secId,C_secId," ARM_ERR: security id: object expected",C_result);
	XL_readStrCell(XL_modId,C_modId," ARM_ERR: model id: object expected",C_result);
	XL_readStrCell(XL_param,C_param," ARM_ERR: param: string expected",C_result);

	if((paramId = ARM_ConvSvtyParam (C_param, C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	long retCode = ARMLOCAL_CvSensitivity (LocalGetNumObjectId (C_secId),
										   LocalGetNumObjectId (C_modId),
										   paramId,
										   C_result);

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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_CvSensitivity" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


__declspec(dllexport) LPXLOPER WINAPI Local_ARM_BetweenDates (LPXLOPER XL_date1,
															  LPXLOPER XL_date2,
															  LPXLOPER XL_daycount,
															  LPXLOPER XL_isYearFrac)
{
	ADD_LOG("Local_ARM_BetweenDates ");
	// return
	static XLOPER XL_result;
	ARM_result C_result;

	// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();

	// C variable
	double C_date1;
	double C_date2;
	
	CCString C_daycount;
	long daycountId;

	double C_isYearFrac;
	double C_isYearFrac_default(1.0);


	static int error;
	static char* reason = "";

	XL_readNumCell(XL_date1,C_date1," ARM_ERR: date1 : numeric expected",C_result);
	XL_readNumCell(XL_date2,C_date2," ARM_ERR: date2 : numeric expected",C_result);
	XL_readStrCell(XL_daycount,C_daycount," ARM_ERR: Daycount Basis: string expected",C_result);
	XL_readNumCellWD(XL_isYearFrac,C_isYearFrac,C_isYearFrac_default," ARM_ERR: isYearFrac: numeric expected",C_result);

	daycountId = ARM_ConvDayCount(C_daycount);

	long retCode = ARMLOCAL_ARM_BetweenDates(C_date1, C_date2, daycountId, C_isYearFrac, C_result);

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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_ARM_BetweenDates" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


__declspec(dllexport) LPXLOPER WINAPI Local_ARM_CountBusinessDays ( LPXLOPER XL_date1,
																	LPXLOPER XL_date2,
																	LPXLOPER XL_calendar)
{
	ADD_LOG("Local_ARM_CountBusinessDays ");
	// return
	static XLOPER XL_result;
	ARM_result C_result;

	// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();

	// C variable
	double C_date1;
	double C_date2;
	
	CCString C_calendar;
	
	static int error;
	static char* reason = "";

	XL_readNumCell(XL_date1,C_date1," ARM_ERR: date1 : numeric expected",C_result);
	XL_readNumCell(XL_date2,C_date2," ARM_ERR: date2 : numeric expected",C_result);
	XL_readStrCell(XL_calendar,C_calendar," ARM_ERR: Calendar : string expected",C_result);

	long retCode = ARMLOCAL_ARM_CountBusinessDays(C_date1, C_date2, C_calendar, C_result);

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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_ARM_CountBusinessDays" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


__declspec(dllexport) LPXLOPER WINAPI Local_ARM_ADDMONTHS (LPXLOPER XL_date,
														   LPXLOPER XL_nb,
														   LPXLOPER XL_rule,
														   LPXLOPER XL_currency)
{
	ADD_LOG("Local_ARM_ADDMONTHS ");
	// return
	static XLOPER XL_result;
	ARM_result C_result;

	// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();

	// C variable
	double C_date;
	double C_nb;
	
	CCString C_rule;
	long ruleId;

	CCString C_currency;

	static int error;
	static char* reason = "";

	XL_readNumCell(XL_date,C_date," ARM_ERR: date : numeric expected",C_result);
	XL_readNumCell(XL_nb,C_nb," ARM_ERR: nb : numeric expected",C_result);
	XL_readStrCellWD(XL_rule,C_rule,"NONE"," ARM_ERR: rule: string expected",C_result);
	XL_readStrCellWD(XL_currency,C_currency,"DEFAULT"," ARM_ERR: currency: string expected",C_result);

	if((ruleId = ARM_ConvRule (C_rule, C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	if(C_currency == "DEFAULT")
	{
		ARM_result currencyres;
		ARMLOCAL_ARM_GetDefaultCurrency (currencyres);
		if(currencyres.getRetCode () != ARM_OK)
		{
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}
		else
		{
			C_currency = currencyres.getString ();
		}
	}

	long retCode = ARMLOCAL_ARM_ADDMONTHS(C_date, C_nb, ruleId, C_currency, C_result);

	if(retCode == ARM_OK)
	{
		FreeCurCellErr ();
		XL_result.xltype = xltypeNum;
		XL_result.val.num = Local_ARMDATE2XLDATE(C_result.getString ());
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_ARM_ADDMONTHS" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


__declspec(dllexport) LPXLOPER WINAPI Local_ARM_ADDYEARS (LPXLOPER XL_date,
														  LPXLOPER XL_nb,
														  LPXLOPER XL_rule,
														  LPXLOPER XL_currency)
{
	ADD_LOG("Local_ARM_ADDYEARS ");
	// return
	static XLOPER XL_result;
	ARM_result C_result;

	// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();

	// C variable
	double C_date;
	double C_nb;

	CCString C_rule;
	long ruleId;

	CCString C_currency;

	static int error;
	static char* reason = "";

	XL_readNumCell(XL_date,C_date," ARM_ERR: date : numeric expected",C_result);
	XL_readNumCell(XL_nb,C_nb," ARM_ERR: nb : numeric expected",C_result);
	XL_readStrCellWD(XL_rule,C_rule,"NONE"," ARM_ERR: rule: string expected",C_result);
	XL_readStrCellWD(XL_currency,C_currency,"DEFAULT"," ARM_ERR: currency: string expected",C_result);

	if((ruleId = ARM_ConvRule (C_rule, C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	if(C_currency == "DEFAULT")
	{
		ARM_result currencyres;
		ARMLOCAL_ARM_GetDefaultCurrency (currencyres);
		if(currencyres.getRetCode () != ARM_OK)
		{
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}
		else
		{
			C_currency = currencyres.getString ();
		}
	}

	long retCode = ARMLOCAL_ARM_ADDYEARS(C_date, C_nb, ruleId, C_currency, C_result);

	if(retCode == ARM_OK)
	{
		FreeCurCellErr ();
		XL_result.xltype = xltypeNum;
		XL_result.val.num = Local_ARMDATE2XLDATE(C_result.getString ());
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_ARM_ADDYEARS" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


__declspec(dllexport) LPXLOPER WINAPI Local_FxConvert (LPXLOPER XL_ccy1,
													   LPXLOPER XL_ccy2,
													   LPXLOPER XL_asOfDate,
													   LPXLOPER XL_amount,
													   LPXLOPER XL_cvname)
{
	ADD_LOG("Local_FxConvert ");
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();

	// C variable
	CCString C_ccy1;
	CCString C_ccy2;

	double C_asOfDate;
	double C_amount;
	double C_amount_default = 1.0;
	
	CCString C_cvname;

	// error
	static int error;
	static char* reason = "";

	XL_readStrCell(XL_ccy1,C_ccy1," ARM_ERR: currency 1: string expected",C_result);
	XL_readStrCell(XL_ccy2,C_ccy2," ARM_ERR: currency 2: string expected",C_result);
	XL_readNumCell(XL_asOfDate,C_asOfDate," ARM_ERR: as of date: date expected",C_result);
	XL_readNumCellWD(XL_amount,C_amount,C_amount_default," ARM_ERR: amount: numeric expected",C_result);
	XL_readStrCellWD(XL_cvname,C_cvname,"MO"," ARM_ERR: cvname: string expected",C_result);
	
	C_ccy1.toUpper ();
	C_ccy2.toUpper ();
	C_cvname.toUpper ();

	long retCode = ARMLOCAL_FxConvert (C_ccy1, C_ccy2, C_asOfDate, C_amount, C_cvname, C_result);

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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_FxConvert" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


__declspec(dllexport) LPXLOPER WINAPI Local_FxConvertFromCalypso (LPXLOPER XL_ccy1,
													   LPXLOPER XL_ccy2,
													   LPXLOPER XL_asOfDate,
													   LPXLOPER XL_cvname)
{
	ADD_LOG("Local_FxConvertFromCalypso ");
	static XLOPER XL_result;
	try 
	{
		std::string ccy1; ExcelTools::convert(XL_ccy1,ccy1); 
		std::string ccy2; ExcelTools::convert(XL_ccy2,ccy2); 
		std::string cvname; ExcelTools::convert(XL_cvname,"MO",cvname); 
		// if (cvname=="") cvname="MO"; 
		ARM_Date date; ExcelTools::convert(XL_asOfDate,date); 
		double fxRate ;
		ARM_CalypsoToolkit::GetFXRate(ccy1,ccy2,cvname,date,fxRate); 
		ExcelTools::convert(fxRate,&XL_result); 
		return &XL_result;		
	}
	catch (Exception&e)
	{
		ExcelCaller::get().get().setError(e.GetErrorString()); 
		ExcelTools::convert("ARM_ERR",&XL_result) ;
	}
	catch (std::exception&e)
	{
		ExcelCaller::get().get().setError(e.what()); 
		ExcelTools::convert("ARM_ERR",&XL_result) ;
	}
	catch(...)
	{
		ExcelCaller::get().setError("Unknown Exception");
		ExcelTools::convert("ARM_ERR",&XL_result) ;
	}
	return &XL_result;
}
 
__declspec(dllexport) LPXLOPER WINAPI Local_ARM_GetCurrency(LPXLOPER XL_Security)
{
	ADD_LOG("Local_ARM_GetCurrency");
	// return
	static XLOPER XL_result;
	ARM_result C_result;

	// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();

	// error
	static int error;
	static char* reason = "";

	CCString C_Security;
	XL_readStrCell(XL_Security,C_Security," ARM_ERR: Security: string expected",C_result);

	long retCode = ARMLOCAL_GetCurrency(LocalGetNumObjectId(C_Security), C_result);;

	if(retCode == ARM_OK)
	{
		FreeCurCellErr ();
		XL_result.xltype = xltypeStr;
		XL_result.val.str = XL_StrC2StrPascal (C_result.getString ());
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_ARM_GetCurrency" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}
__declspec(dllexport) LPXLOPER WINAPI Local_ARM_GetDefaultCurrency ()
{
	ADD_LOG("Local_ARM_GetDefaultCurrency ");
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();

	// error
	static int error;
	static char* reason = "";

	long retCode = ARMLOCAL_ARM_GetDefaultCurrency (C_result);

	if(retCode == ARM_OK)
	{
		FreeCurCellErr ();
		XL_result.xltype = xltypeStr;
		XL_result.val.str = XL_StrC2StrPascal (C_result.getString ());
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_ARM_GetDefaultCurrency" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}

__declspec(dllexport) LPXLOPER WINAPI Local_ARM_SetDefaultCurrency (LPXLOPER XL_currency)
{
	ADD_LOG("Local_ARM_SetDefaultCurrency ");
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();

	// C variable
	CCString C_currency;
	
	// error
	static int error;
	static char* reason = "";

	XL_readStrCell(XL_currency,C_currency," ARM_ERR: currency: string expected",C_result);
	
	long retCode = ARMLOCAL_ARM_SetDefaultCurrency (C_currency, C_result);

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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_ARM_SetDefaultCurrency" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


__declspec(dllexport) LPXLOPER WINAPI Local_ARM_ADDPERIOD (LPXLOPER XL_date,
														   LPXLOPER XL_freq,
														   LPXLOPER XL_ccy,
														   LPXLOPER XL_nbPeriods,
														   LPXLOPER XL_adjRule,
														   LPXLOPER XL_goToEndOfMonth)
{
	ADD_LOG("Local_ARM_ADDPERIOD ");
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();

	// C variable
	double C_date;

	CCString C_freq;
	long freqId;

	CCString C_ccy;

	double C_nbPeriods;
	double C_nbPeriods_default = 1;

	CCString C_adjRule;
	long adjRuleId;

	double C_goToEndOfMonth;
	double C_goToEndOfMonth_default = 0.0;
	long   goToEndOfMonth;

	static int error;
	static char* reason = "";

	XL_readNumCell(XL_date,C_date," ARM_ERR: date : numeric expected",C_result);
	XL_readStrCell(XL_freq,C_freq," ARM_ERR: frequency : string expected",C_result);
	XL_readStrCellWD(XL_ccy,C_ccy,"DEFAULT"," ARM_ERR: currency: string expected",C_result);
	XL_readNumCellWD(XL_nbPeriods,C_nbPeriods,C_nbPeriods_default," ARM_ERR: number of periods : numeric expected",C_result);
	XL_readStrCellWD(XL_adjRule,C_adjRule,"NONE"," ARM_ERR: adjusting Rule: string expected",C_result);
	XL_readNumCellWD(XL_goToEndOfMonth,C_goToEndOfMonth,C_goToEndOfMonth_default," ARM_ERR: Go to End of Month : 0 or 1 expected",C_result);

	if(C_ccy == "DEFAULT")
	{
		ARM_result currencyres;
		ARMLOCAL_ARM_GetDefaultCurrency (currencyres);
		if(currencyres.getRetCode () != ARM_OK)
		{
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}
		else
		{
			C_ccy = currencyres.getString ();
		}
	}

	long retCode;

	goToEndOfMonth = (long) C_goToEndOfMonth;

	if((adjRuleId = ARM_ConvRule (C_adjRule, C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	if ( (freqId = ARM_ConvFrequency(C_freq,C_result)) == ARM_DEFAULT_ERR)
	{
		int Nb;
		char matu;

		sscanf(C_freq, "%d%c", &Nb, &matu);

        matu = toupper(matu);

        if ( matu == 'D' ) // Ex : "1D"
        {    
			retCode = ARMLOCAL_ARM_ADDPERIOD(C_date, K_DAILY, C_ccy, (long) (C_nbPeriods*Nb), adjRuleId, goToEndOfMonth, C_result);
        }
        else if ( matu == 'W' )  
        {   //  Ex : "1W"    

			retCode = ARMLOCAL_ARM_ADDPERIOD(C_date, K_WEEKLY, C_ccy, (long) (C_nbPeriods*Nb), adjRuleId, goToEndOfMonth, C_result);
        }
        else if ( matu == 'M' ) 
        {   //  Ex : "9M"
			retCode = ARMLOCAL_ARM_ADDPERIOD(C_date, K_MONTHLY, C_ccy, (long) (C_nbPeriods*Nb), adjRuleId, goToEndOfMonth, C_result);
        }
        else if ( matu == 'Y')  // ->implicitement ce sont des taux de swap
        {   
			retCode = ARMLOCAL_ARM_ADDPERIOD(C_date, K_ANNUAL, C_ccy, (long) (C_nbPeriods*Nb), adjRuleId, goToEndOfMonth, C_result);
		}
		else
		{
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}
	}
	else
	{
		retCode = ARMLOCAL_ARM_ADDPERIOD(C_date, freqId, C_ccy, (long) C_nbPeriods, adjRuleId, goToEndOfMonth, C_result);
	}

	if(retCode == ARM_OK)
	{
		FreeCurCellErr ();
		XL_result.xltype = xltypeNum;
		XL_result.val.num = Local_ARMDATE2XLDATE(C_result.getString ());
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_ARM_ADDPERIOD" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}



__declspec(dllexport) LPXLOPER WINAPI Local_ARM_Price_OptUnder (LPXLOPER XL_secId,
																LPXLOPER XL_modId)
{
	ADD_LOG("Local_ARM_Price_OptUnder ");
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	LPXLOPER pxArray;
	ARM_result C_result;

	// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();

	// C variable
	CCString C_secId;
	CCString C_modId;
	
	// error
	static int error;
	static char* reason = "";

	XL_readStrCell(XL_secId,C_secId," ARM_ERR: security id: object expected",C_result);
	XL_readStrCell(XL_modId,C_modId," ARM_ERR: model id: object expected",C_result);
	
	long retCode = ARMLOCAL_ARM_Price_OptUnder (LocalGetNumObjectId (C_secId),
												LocalGetNumObjectId (C_modId),
												C_result);

	if ( retCode == ARM_OK )
	{
		int nbrows = (int) (C_result.getDouble());
		int nbcolumns = 1;
	
		FreeCurCellErr ();

		XL_result.xltype = xltypeMulti;
		XL_result.val.array.columns = nbcolumns;
		XL_result.val.array.rows = nbrows;
		XL_result.val.array.lparray = pxArray = (LPXLOPER)GlobalAlloc (GMEM_ZEROINIT, nbrows * nbcolumns * sizeof (XLOPER));

		for (int i = 0; i < nbrows; i++)
		{
			pxArray[XL_Coordonnate2Rank(i, 0, nbcolumns)].xltype = xltypeNum;
			pxArray[XL_Coordonnate2Rank(i, 0, nbcolumns)].val.num = C_result.getArray(i);
		}
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_ARM_Price_OptUnder" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


__declspec(dllexport) LPXLOPER WINAPI Local_ARM_GetPID(void)
{
	ADD_LOG("Local_ARM_GetPID");
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;
	// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();
    
    char* user = getenv ("USERNAME");
	char buf[100];
#ifdef _DEBUG
	sprintf(buf,"[ARMLOCAL_CLIENT(Debug Version:%s:%s,compiled by %s.)]",__DATE__, __TIME__, user);
#else
	sprintf(buf,"[ARMLOCAL_CLIENT(Release Version:%s:%s)]",__DATE__, __TIME__);
#endif

	FreeCurCellErr ();

	XL_result.xltype = xltypeStr;
	XL_result.val.str = XL_StrC2StrPascal((const CCString&) buf);
	XL_result.xltype |= xlbitDLLFree;

//	ARM_END();
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_ARM_GetPID" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


__declspec(dllexport) LPXLOPER WINAPI Local_ParallelShift (LPXLOPER XL_curveId,
														   LPXLOPER XL_value)
{
	ADD_LOG("Local_ParallelShift ");
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();

	// C variable
	CCString C_curveId;
	double C_value;

	// error
	static int error;
	static char* reason = "";

	XL_readStrCell(XL_curveId,C_curveId," ARM_ERR: curve id: object expected",C_result);
	XL_readNumCell(XL_value,C_value," ARM_ERR: value: numeric expected",C_result);
	
	long retCode;
	long objId;
	CCString prevClass;
	
	CCString curClass = LocalGetStringObjectClass (C_curveId);
	CCString stringId = GetLastCurCellEnvValue ();
	
	if(!stringId)
	{
		retCode = ARMLOCAL_ParallelShift (LocalGetNumObjectId (C_curveId), C_value, C_result);

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
			retCode = ARMLOCAL_ParallelShift (LocalGetNumObjectId (C_curveId), C_value, C_result, objId);
		
			if(retCode == ARM_OK)
			{
				LocalSetCurCellEnvValue (curClass, objId); 

				stringId = LocalMakeObjectId (objId, curClass);
			}
		}
		else
		{
			FreeCurCellContent ();
			retCode = ARMLOCAL_ParallelShift (LocalGetNumObjectId (C_curveId), C_value, C_result);
		
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_ParallelShift" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}

__declspec(dllexport) LPXLOPER WINAPI Local_PXL_ParallelShift (LPXLOPER XL_curveId,
															   LPXLOPER XL_value)
{
	ADD_LOG("Local_PXL_ParallelShift ");
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();

	// C variable
	CCString C_curveId;
	double C_value;

	// error
	static int error;
	static char* reason = "";

	XL_readStrCell(XL_curveId,C_curveId," ARM_ERR: curve id: object expected",C_result);
	XL_readNumCell(XL_value,C_value," ARM_ERR: value: numeric expected",C_result);
	
	long retCode;
	long objId;
	
	CCString curClass = LocalGetStringObjectClass (C_curveId);
	CCString stringId;

	retCode = ARMLOCAL_ParallelShift (LocalGetNumObjectId (C_curveId), C_value, C_result);

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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_PXL_ParallelShift" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


__declspec(dllexport) LPXLOPER WINAPI Local_ARM_ClonedAndSetNotional (LPXLOPER XL_secId,
																	  LPXLOPER XL_refValId,
																	  LPXLOPER XL_percentRemainder,
																	  LPXLOPER XL_datesForInterpol)
{
	ADD_LOG("Local_ARM_ClonedAndSetNotional ");
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();

	// C variable
	CCString C_secId;
	CCString C_refValId;

	double C_percentRemainder;
	double C_percentRemainder_default = 100.0;

	double C_datesForInterpol;
	double C_datesForInterpol_default = 3.0;

	// error
	static int error;
	static char* reason = "";

	XL_readStrCell(XL_secId,C_secId," ARM_ERR: security id: object expected",C_result);
	XL_readStrCell(XL_refValId,C_refValId," ARM_ERR: reference value id: object expected",C_result);
	XL_readNumCellWD(XL_percentRemainder,C_percentRemainder,C_percentRemainder_default," ARM_ERR: percentage of remainder: numeric expected",C_result);
	XL_readNumCellWD(XL_datesForInterpol,C_datesForInterpol,C_datesForInterpol_default," ARM_ERR: dates on which to interpolate: numeric expected",C_result);

	long retCode;
	long objId;
	CCString prevClass;
	
	CCString curClass = LocalGetStringObjectClass(C_secId);
	CCString stringId = GetLastCurCellEnvValue ();
	
	if(!stringId)
	{
		retCode = ARMLOCAL_ClonedAndSetNotional(LocalGetNumObjectId(C_secId),
												LocalGetNumObjectId(C_refValId),
												C_percentRemainder,
												C_result,
												(long)C_datesForInterpol);

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
			retCode = ARMLOCAL_ClonedAndSetNotional(LocalGetNumObjectId(C_secId),
													LocalGetNumObjectId(C_refValId),
													C_percentRemainder,
													C_result,
													(long)C_datesForInterpol,
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

			retCode = ARMLOCAL_ClonedAndSetNotional(LocalGetNumObjectId(C_secId),
													LocalGetNumObjectId(C_refValId),
													C_percentRemainder,
													C_result,
													(long)C_datesForInterpol);
		
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_ARM_ClonedAndSetNotional" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


__declspec(dllexport) LPXLOPER WINAPI Local_PXL_ARM_ClonedAndSetNotional (LPXLOPER XL_secId,
																		  LPXLOPER XL_refValId,
																		  LPXLOPER XL_percentRemainder,
																		  LPXLOPER XL_datesForInterpol)
{
	ADD_LOG("Local_PXL_ARM_ClonedAndSetNotional ");
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();

	// C variable
	CCString C_secId;
	CCString C_refValId;

	double C_percentRemainder;
	double C_percentRemainder_default = 100.;

	double C_datesForInterpol;
	double C_datesForInterpol_default = 3.0;

	// error
	static int error;
	static char* reason = "";

	XL_readStrCell(XL_secId,C_secId," ARM_ERR: security id: object expected",C_result);
	XL_readStrCell(XL_refValId,C_refValId," ARM_ERR: reference value id: object expected",C_result);
	XL_readNumCellWD(XL_percentRemainder,C_percentRemainder,C_percentRemainder_default," ARM_ERR: percentage of remainder: numeric expected",C_result);
	XL_readNumCellWD(XL_datesForInterpol,C_datesForInterpol,C_datesForInterpol_default," ARM_ERR: dates on which to interpolate: numeric expected",C_result);

	long retCode;
	long objId;

	CCString curClass = LocalGetStringObjectClass (C_secId);
	CCString stringId;

	retCode = ARMLOCAL_ClonedAndSetNotional(LocalGetNumObjectId(C_secId),
											LocalGetNumObjectId(C_refValId),
											C_percentRemainder,
											C_result,
											(long)C_datesForInterpol);

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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_PXL_ARM_ClonedAndSetNotional" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}





__declspec(dllexport) LPXLOPER WINAPI Local_ARM_ClonedAndSet (LPXLOPER XL_spdoptId,
															  LPXLOPER XL_valToSet,
															  LPXLOPER XL_typeToSet)
{
	ADD_LOG("Local_ARM_ClonedAndSet ");
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();

	// C variable
	CCString C_spdoptId;

	double   C_valtoset_double;
	CCString C_valtoset_str;
	long     valtosetType;

	CCString C_typeToSet;

	// error
	static int error;
	static char* reason = "";

	XL_readStrCell(XL_spdoptId,C_spdoptId," ARM_ERR: security id: object expected",C_result);
	XL_readStrOrNumCell(XL_valToSet, C_valtoset_str, C_valtoset_double, valtosetType,
			   " ARM_ERR: Value to set: numeric or object ID string expected",C_result);
	XL_readStrCellWD(XL_typeToSet,C_typeToSet,"STRIKE"," ARM_ERR: type: string expected",C_result);

	C_typeToSet.toUpper();

	if ( valtosetType == XL_TYPE_STRING )
	{
		C_valtoset_double = (double) LocalGetNumObjectId(C_valtoset_str);

		valtosetType = 1L;
	}
	else
	{
		valtosetType = 0L;
	}

	long retCode;
	long objId;
	CCString prevClass;
	
	CCString curClass = LocalGetStringObjectClass(C_spdoptId);
	CCString stringId = GetLastCurCellEnvValue ();

	if(!stringId)
	{
		retCode = ARMLOCAL_ClonedAndSet(LocalGetNumObjectId(C_spdoptId),
										valtosetType,
										C_valtoset_double,
										C_typeToSet,
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
			retCode = ARMLOCAL_ClonedAndSet(LocalGetNumObjectId(C_spdoptId),
											valtosetType,
											C_valtoset_double,
											C_typeToSet,
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

			retCode = ARMLOCAL_ClonedAndSet(LocalGetNumObjectId(C_spdoptId),
											valtosetType,
											C_valtoset_double,
											C_typeToSet,
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_ARM_ClonedAndSetSlopeFlag" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


__declspec(dllexport) LPXLOPER WINAPI Local_PXL_ARM_ClonedAndSet (LPXLOPER XL_spdoptId,
																  LPXLOPER XL_valToSet,
																  LPXLOPER XL_typeToSet)
{
	ADD_LOG("Local_PXL_ARM_ClonedAndSet ");
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();

	// C variable
	CCString C_spdoptId;

	double   C_valtoset_double;
	CCString C_valtoset_str;
	long     valtosetType;

	CCString C_typeToSet;

	// error
	static int error;
	static char* reason = "";

	XL_readStrCell(XL_spdoptId,C_spdoptId," ARM_ERR: security id: object expected",C_result);
	XL_readStrOrNumCell(XL_valToSet, C_valtoset_str, C_valtoset_double, valtosetType,
			   " ARM_ERR: Value to set: numeric or object ID string expected",C_result);
	XL_readStrCellWD(XL_typeToSet,C_typeToSet,"STRIKE"," ARM_ERR: type: string expected",C_result);

	C_typeToSet.toUpper();

	if ( valtosetType == XL_TYPE_STRING )
	{
		C_valtoset_double = (double) LocalGetNumObjectId(C_valtoset_str);

		valtosetType = 1L;
	}
	else
	{
		valtosetType = 0L;
	}

	long retCode;
	long objId;
	CCString prevClass;
	
	CCString curClass = LocalGetStringObjectClass(C_spdoptId);
	CCString stringId = GetLastCurCellEnvValue ();

	retCode = ARMLOCAL_ClonedAndSet(LocalGetNumObjectId(C_spdoptId),
									valtosetType,
									C_valtoset_double,
									C_typeToSet,
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_ARM_ClonedAndSet" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


__declspec(dllexport) LPXLOPER WINAPI Local_PXL_ARM_Clone (LPXLOPER XL_objectId)
{
	ADD_LOG("Local_PXL_ARM_Clone ");
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();

	// C variable
	CCString C_objectId;

	// error
	static int error;
	static char* reason = "";

	XL_readStrCell(XL_objectId,C_objectId," ARM_ERR: security id: object expected",C_result);

	long retCode;
	long objId;

	CCString curClass = LocalGetStringObjectClass (C_objectId);
	CCString stringId;

	retCode = ARMLOCAL_Clone(LocalGetNumObjectId(C_objectId),
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_PXL_ARM_Clone" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


__declspec(dllexport) LPXLOPER WINAPI Local_ARM_DoPastReset (LPXLOPER XL_secId,
															 LPXLOPER XL_resetMgrId,
															 LPXLOPER XL_AsOf)
{
	ADD_LOG("Local_ARM_DoPastReset ");
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();

	// C variable
	CCString C_secId;
	VECTOR<CCString> C_resetMgrId;
	double C_AsOf;

	// error
	static int error;
	static char* reason = "";

	XL_readStrCell(XL_secId,C_secId," ARM_ERR: security id: object expected",C_result);
	XL_readStrVector(XL_resetMgrId, C_resetMgrId, " ARM_ERR: resetMgrIds: list of object ID string expected",XL_TYPE_STRING,C_result);
	XL_readNumCell(XL_AsOf,C_AsOf," ARM_ERR: AsOf: numeric expected",C_result);

	long retCode;
	long objId;
	CCString prevClass;
	
	CCString curClass = LocalGetStringObjectClass(C_secId);
	CCString stringId = GetLastCurCellEnvValue ();

	if(!stringId)
	{
		retCode = ARMLOCAL_DoPastReset(LocalGetNumObjectId(C_secId),
									   C_resetMgrId,
									   (long)C_AsOf,
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
			retCode = ARMLOCAL_DoPastReset(LocalGetNumObjectId(C_secId),
										   C_resetMgrId,
										   (long)C_AsOf,
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

			retCode = ARMLOCAL_DoPastReset(LocalGetNumObjectId(C_secId),
										   C_resetMgrId,
										   (long)C_AsOf,
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_ARM_DoPastReset" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


__declspec(dllexport) LPXLOPER WINAPI Local_PXL_ARM_DoPastReset (LPXLOPER XL_secId,
																 LPXLOPER XL_resetMgrId,
																 LPXLOPER XL_AsOf)
{
	ADD_LOG("Local_PXL_ARM_DoPastReset ");
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();

	// C variable
	CCString C_secId;
	VECTOR<CCString> C_resetMgrId;
	double C_AsOf;

	// error
	static int error;
	static char* reason = "";

	XL_readStrCell(XL_secId,C_secId," ARM_ERR: security id: object expected",C_result);
	XL_readStrVector(XL_resetMgrId, C_resetMgrId, " ARM_ERR: resetMgrIds: list of object ID string expected",XL_TYPE_STRING,C_result);
	XL_readNumCell(XL_AsOf,C_AsOf," ARM_ERR: AsOf: numeric expected",C_result);

	long retCode;
	long objId;
	CCString prevClass;
	
	CCString curClass = LocalGetStringObjectClass(C_secId);
	CCString stringId;

	retCode = ARMLOCAL_DoPastReset(LocalGetNumObjectId(C_secId),
								   C_resetMgrId,
								   (long)C_AsOf,
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_PXL_ARM_DoPastReset" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


__declspec(dllexport) LPXLOPER WINAPI Local_FIXRATES (LPXLOPER XL_secId,
													  LPXLOPER XL_rate)
{
	ADD_LOG("Local_FIXRATES ");
	// return
	static XLOPER XL_result;
	ARM_result C_result;

	// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();

	// C variable
	CCString C_secId;
	VECTOR<double> C_rate;
	
	// error
	static int error;
	static char* reason = "";

	XL_readStrCell(XL_secId,C_secId," ARM_ERR: seurity id: object expected",C_result);
	XL_readNumVector(XL_rate,C_rate," ARM_ERR: rate: array of numeric expected",C_result);
	
	long retCode;
	long objId;
	CCString prevClass;

	CCString curClass = LocalGetStringObjectClass (C_secId);
	CCString stringId = GetLastCurCellEnvValue ();

	if(!stringId)
	{
		retCode = ARMLOCAL_FIXRATES(LocalGetNumObjectId(C_secId),
									C_rate,
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
			retCode = ARMLOCAL_FIXRATES(LocalGetNumObjectId(C_secId),
										C_rate,
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

			retCode = ARMLOCAL_FIXRATES(LocalGetNumObjectId(C_secId),
										C_rate,
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_FIXRATES" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}



__declspec(dllexport) LPXLOPER WINAPI Local_PXL_FIXRATES (LPXLOPER XL_secId,
														  LPXLOPER XL_rate)
{
	ADD_LOG("Local_PXL_FIXRATES ");
	// return
	static XLOPER XL_result;
	ARM_result C_result;

	// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();

	// C variable
	CCString C_secId;
	VECTOR<double> C_rate;
	
	// error
	static int error;
	static char* reason = "";

	XL_readStrCell(XL_secId,C_secId," ARM_ERR: seurity id: object expected",C_result);
	XL_readNumVector(XL_rate,C_rate," ARM_ERR: rate: array of numeric expected",C_result);
	
	long retCode;
	long objId;

	CCString curClass = LocalGetStringObjectClass (C_secId);
	CCString stringId;
	
	retCode = ARMLOCAL_FIXRATES(LocalGetNumObjectId(C_secId),
								C_rate,
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_PXL_FIXRATES" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


__declspec(dllexport) LPXLOPER WINAPI Local_ARM_DisplayScheduleValues (LPXLOPER XL_instId,
																	   LPXLOPER XL_typeValues,
																	   LPXLOPER XL_RecOrPay,
																	   LPXLOPER XL_modelId)
{
	ADD_LOG("Local_ARM_DisplayScheduleValues ");
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();

	// C variable
	CCString C_instId;
	
	CCString C_typeValues;
	long typeValuesId;

	CCString C_modelId;
	long modelId;

	CCString C_RecOrPay;
	long recId;

	// error
	static int error;
	static char* reason = "";

	XL_readStrCell(XL_instId,C_instId," ARM_ERR: instrument id: object expected",C_result);
	XL_readStrCell(XL_typeValues,C_typeValues," ARM_ERR: type Values: string expected",C_result);
	XL_readStrCellWD(XL_RecOrPay,C_RecOrPay,"R"," ARM_ERR: Rec or Pay: string expected",C_result);
	XL_readStrCellWD(XL_modelId,C_modelId,"DEFAULT"," ARM_ERR: model Id: string expected",C_result);

	if((typeValuesId = ARM_ConvTypeValues (C_typeValues, C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	if((recId = ARM_ConvRecOrPay (C_RecOrPay, C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	if(C_modelId == "DEFAULT")
	{
		modelId = ARM_NULL_OBJECT;
	}
	else
	{
		modelId = LocalGetNumObjectId (C_modelId);
	}

	long retCode = ARMLOCAL_ARM_DisplayScheduleValues(LocalGetNumObjectId(C_instId),
													  typeValuesId,
													  recId,
													  modelId,
													  C_result);

	if(retCode == ARM_OK)
	{
		VECTOR<double> dVal;
		long vecSize;
		retCode = ExtractVectorDoubleFromFile("123",dVal,vecSize);
		/// AAAAAARRRRRRRGGGGGGG the extract routine extracts one more data!
		dVal.resize(vecSize);

		if ( (retCode == ARM_OK) )
		{
			XL_writeNumVector( XL_result, dVal, " ARM_ERR: Could not set value", C_result );
		}
		else
		{
			ARM_ERR();
		}
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_ARM_DisplayScheduleValues" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


__declspec(dllexport) LPXLOPER WINAPI Local_ARM_DisplayScheduleDates (LPXLOPER XL_instId,
																	  LPXLOPER XL_typeDates,
																	  LPXLOPER XL_RecOrPay,
																	  LPXLOPER XL_viewInitExch)
{
	ADD_LOG("Local_ARM_DisplayScheduleDates ");
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();

	// C variable
	CCString C_instId;
	
	CCString C_typeDates;
	long typeDatesId;

	CCString C_RecOrPay;
	long recId;

	CCString C_viewInitExch;
	long viewInitExchId;

	// error
	static int error;
	static char* reason = "";

	XL_readStrCell(XL_instId,C_instId," ARM_ERR: instrument id: object expected",C_result);
	XL_readStrCell(XL_typeDates,C_typeDates," ARM_ERR: type Values: string expected",C_result);
	XL_readStrCellWD(XL_RecOrPay,C_RecOrPay,"R"," ARM_ERR: Rec or Pay: string expected",C_result);
	XL_readStrCellWD(XL_viewInitExch,C_viewInitExch,"Y"," ARM_ERR: view Initial Exch: string expected",C_result);

	if((typeDatesId = ARM_ConvTypeDates (C_typeDates, C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	if((recId = ARM_ConvRecOrPay (C_RecOrPay, C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	if((viewInitExchId = ARM_ConvYesOrNo (C_viewInitExch, C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	long retCode = ARMLOCAL_ARM_DisplayScheduleDates(LocalGetNumObjectId(C_instId),
													 typeDatesId,
													 recId,
													 viewInitExchId,
													 C_result);

	if(retCode == ARM_OK)
	{
		VECTOR<CCString> dDate;
		long vecSize;
		retCode = ExtractVectorDateFromFile("123",dDate,vecSize);

		/// AAAAAARRRRRRRGGGGGGG the extract routine extracts one more data!
		dDate.resize(vecSize);

		if ( (retCode == ARM_OK) )
		{
			VECTOR<double> dXLDateResult( dDate.size() );
			int i;
			for(i=0;i<dDate.size(); ++i)
				dXLDateResult[i] = Local_ARMDATE2XLDATE( dDate[i] );

			XL_writeNumVector( XL_result, dXLDateResult, " ARM_ERR: Could not set value", C_result);
		}
		else
		{
			ARM_ERR();
		}
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_ARM_DisplayScheduleDates" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}

__declspec(dllexport) LPXLOPER WINAPI Local_ARM_DisplayReplicPort(
LPXLOPER XL_instId,
LPXLOPER XL_WeightOrStrike,
LPXLOPER XL_PayoffOrSensi,
LPXLOPER XL_RecOrPay,
LPXLOPER XL_ModelId)
{
	ADD_LOG("Local_ARM_DisplayReplicPort");
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;
    VECTOR<double> C_DataResult;

	// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();

	// C variable
	CCString C_instId;
	
	CCString C_WeightOrStrike;
	long WeightOrStrike;

    CCString C_PayoffOrSensi;
	long PayoffOrSensi;

	CCString C_RecOrPay;
	long recId;

    CCString C_modelId;
	long modelId;

	// error
	static int error;
	static char* reason = "";

	XL_readStrCell(XL_instId,C_instId," ARM_ERR: instrument id: object expected",C_result);
	XL_readStrCell(XL_WeightOrStrike,C_WeightOrStrike," ARM_ERR: Weight or Strike: string expected",C_result);
    XL_readStrCell(XL_PayoffOrSensi,C_PayoffOrSensi," ARM_ERR: PayoffOrSensi or Strike: string expected",C_result);
	XL_readStrCellWD(XL_RecOrPay,C_RecOrPay,"R"," ARM_ERR: Rec or Pay: string expected",C_result);
    XL_readStrCellWD(XL_ModelId,C_modelId,"DEFAULT"," ARM_ERR: model Id: string expected",C_result);

	if(C_WeightOrStrike == "WEIGHT")
	{
		WeightOrStrike = 1;
	}
    else if(C_WeightOrStrike == "STRIKE")
    {
        WeightOrStrike = 0;
    }
    else
    {
        ARM_ARG_ERR();
    }

    if(C_PayoffOrSensi == "PAYOFF")
	{
		PayoffOrSensi = 1;
	}
    else if(C_PayoffOrSensi == "SENSI")
    {
        PayoffOrSensi = 0;
    }
    else
    {
        ARM_ARG_ERR();
    }

	if((recId = ARM_ConvRecOrPay (C_RecOrPay, C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

    if(C_modelId == "DEFAULT")
	{
		modelId = ARM_NULL_OBJECT;
	}
	else
	{
		modelId = LocalGetNumObjectId (C_modelId);
	}

    long retCode = ARMLOCAL_ARM_DisplayReplicPortfolio(
        LocalGetNumObjectId(C_instId),
        WeightOrStrike,
        PayoffOrSensi,
        recId,
        modelId,
        C_DataResult,
        C_result);

    if (retCode == ARM_OK)
	{
        XL_writeNumVector( XL_result, C_DataResult, " ARM_ERR: Could not set value", C_result );
    }

//	ARM_END();
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_ARM_DisplayReplicPort" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}



__declspec(dllexport) LPXLOPER WINAPI Local_ARM_INTERPOL (LPXLOPER XL_vexX,
														  LPXLOPER XL_vecY,
														  LPXLOPER XL_X,
														  LPXLOPER XL_typeInterpol)
{
	ADD_LOG("Local_ARM_INTERPOL ");
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();

	// C variable
	VECTOR<double> C_vecX;
	VECTOR<double> C_vecY;

	double C_X;
	
	CCString C_interp;
	long interpId;

	static int error;
	static char* reason = "";

	XL_readNumVector(XL_vexX,C_vecX," ARM_ERR: X vector: array of numeric expected",C_result);
	XL_readNumVector(XL_vecY,C_vecY," ARM_ERR: Y vector: array of numeric expected",C_result);
	XL_readNumCell(XL_X,C_X," ARM_ERR: x: numeric expected",C_result);
	XL_readStrCellWD(XL_typeInterpol,C_interp,"LINEAR"," ARM_ERR: interpolation type: string expected",C_result);

	if ( (interpId = ARM_ConvInterpMethod(C_interp,C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	long retCode = ARMLOCAL_INTERPOL(C_vecX, C_vecY, C_X, interpId, C_result);

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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_ARM_INTERPOL" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


__declspec(dllexport) LPXLOPER WINAPI Local_ARM_TRIANGULARINTERPOL (LPXLOPER XL_vexX,
																	LPXLOPER XL_vecY,
																	LPXLOPER XL_matZ,
																	LPXLOPER XL_X,
																	LPXLOPER XL_Y)
{
	ADD_LOG("Local_ARM_TRIANGULARINTERPOL ");
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();

	// C variable
	VECTOR<double> C_vecX;
	VECTOR<double> C_vecY;
	VECTOR<double> C_matZ;

	double C_X;
	double C_Y;
	
	static int error;
	static char* reason = "";

	XL_readNumVector(XL_vexX,C_vecX," ARM_ERR: X vector: array of numeric expected",C_result);
	XL_readNumVector(XL_vecY,C_vecY," ARM_ERR: Y vector: array of numeric expected",C_result);
	XL_readNumVector(XL_matZ,C_matZ," ARM_ERR: Z matrix: array of numeric expected",C_result);
	XL_readNumCell(XL_X,C_X," ARM_ERR: x: numeric expected",C_result);
	XL_readNumCell(XL_Y,C_Y," ARM_ERR: x: numeric expected",C_result);

	long retCode = ARMLOCAL_TRIANGULARINTERPOL(C_vecX, C_vecY, C_matZ, C_X, C_Y, C_result);

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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_ARM_TRIANGULARINTERPOL" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


__declspec(dllexport) LPXLOPER WINAPI Local_KImp (LPXLOPER XL_sec,
												  LPXLOPER XL_model,
												  LPXLOPER XL_price,
												  LPXLOPER XL_param)
{
	ADD_LOG("Local_KImp ");
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();

	// C variable
	CCString C_sec;
	CCString C_model;
	double C_price;

	CCString C_param;
	long paramId;
		
	// error
	static int error;
	static char* reason = "";

	XL_readStrCell(XL_sec,C_sec," ARM_ERR: security id: object expected",C_result);
	XL_readStrCell(XL_model,C_model," ARM_ERR: model id: object expected",C_result);
	XL_readNumCell(XL_price,C_price," ARM_ERR: price: numeric expected",C_result);
	XL_readStrCellWD(XL_param,C_param,"PRICE"," ARM_ERR: param: string expected",C_result);

	if((paramId = ARM_ConvParam (C_param, C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}
	
	long retCode = ARMLOCAL_KImp (LocalGetNumObjectId (C_sec),
								  LocalGetNumObjectId (C_model),
								  C_price,
								  paramId,
								  C_result);
		
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_KImp" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


__declspec(dllexport) LPXLOPER WINAPI Local_ARM_DisplayZC (LPXLOPER XL_zcId)
{
	ADD_LOG("Local_ARM_DisplayZC ");
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	LPXLOPER pxArray;
	ARM_result C_result;

	// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();

	// C variable
	CCString C_zcId;
	
	// error
	static int error;
	static char* reason = "";

	XL_readStrCell(XL_zcId,C_zcId," ARM_ERR: Zc id: object expected",C_result);

	long retCode = ARMLOCAL_ARM_DisplayZC(LocalGetNumObjectId(C_zcId),
										  C_result);

	if(retCode == ARM_OK)
	{
		if(retCode == ARM_OK)
		{
			int nbrows = C_result.getLong();
			int nbcolumns = 3;
		
			FreeCurCellErr ();
			XL_result.xltype = xltypeMulti;
			XL_result.val.array.columns = nbcolumns;
			XL_result.val.array.rows = nbrows; 
			XL_result.val.array.lparray = pxArray = (LPXLOPER)GlobalAlloc (GMEM_ZEROINIT, nbrows * nbcolumns * sizeof (XLOPER));

			for(int i = 0; i < nbrows; i++)
			{
				pxArray[XL_Coordonnate2Rank (i, 0, nbcolumns)].xltype = xltypeNum;
				pxArray[XL_Coordonnate2Rank (i, 0, nbcolumns)].val.num = C_result.getArray(i);

				pxArray[XL_Coordonnate2Rank (i, 1, nbcolumns)].xltype = xltypeNum;
				pxArray[XL_Coordonnate2Rank (i, 1, nbcolumns)].val.num = C_result.getArray(nbrows+i);

				pxArray[XL_Coordonnate2Rank (i, 2, nbcolumns)].xltype = xltypeNum;
				pxArray[XL_Coordonnate2Rank (i, 2, nbcolumns)].val.num = C_result.getArray(nbrows*2+i);
			}
		}
		else
		{
			ARM_ERR();
		}
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_ARM_DisplayZC" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


__declspec(dllexport) LPXLOPER WINAPI Local_BSSpot (LPXLOPER XL_secId,
													LPXLOPER XL_modId,
													LPXLOPER XL_date)
{
	ADD_LOG("Local_BSSpot ");
//	ARM_BEGIN();

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
	double C_date;
	
	// error
	static int error;
	static char* reason = "";

	XL_readStrCell(XL_secId,C_secId," ARM_ERR: security id: object expected",C_result);
	XL_readStrCell(XL_modId,C_modId," ARM_ERR: model id: object expected",C_result);
	XL_readNumCell(XL_date,C_date," ARM_ERR: date: date expected",C_result);
	
	long retCode = ARMLOCAL_BSSpot (LocalGetNumObjectId (C_secId), LocalGetNumObjectId (C_modId), C_date, C_result);

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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_BSSpot" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}

__declspec(dllexport) int WINAPI Local_ARM_Help (void)
{
    CCString msg = CCString (IDH_Advanced_Risk_Management);
	Excel(xlfHelp,0, 1,  (LPXLOPER) TempStr(((char*)(const char*)msg)));
	return 0;
}


__declspec(dllexport) void WINAPI Local_ARM_ProdConnect(void)
{
	if (GetDataRetrieverVersion() == WSETKRETRIEVER)
	{
		connection_wsetoolkit("PROD","OTC","wsotc");
	}
	else
	{
		switchToETK();

		connection_etoolkit(SUMMIT_PROD_CONNEXION_USERNAME,
							SUMMIT_PROD_CONNEXION_PASSWD,
							SUMMIT_PROD_CONNEXION_CONTEXT,
							SUMMIT_PROD_IT_CONFIG_DOMAINSDIR,
							SUMMIT_PROD_IT_DOMAIN_NAME);
	}
}

__declspec(dllexport) void WINAPI Local_ARM_CalypsoDevConnect(void)
{
	try {
	
	char * user = getenv("USERNAME") ;
	std::stringstream sstr ;
	sstr<<"-user "<<user<<" -password "<<user<<" -env calypso_dev" ; 
	ARM_CalypsoToolkit::init(sstr.str()); 
	}
	catch(...)
	{}
}
__declspec(dllexport) void WINAPI Local_ARM_CalypsoProdConnect(void)
{
	try {
	
	char * user = getenv("USERNAME") ;
	std::stringstream sstr ;
	sstr<<"-user "<<user<<" -password "<<user<<" -env calypso_prod" ; 
	ARM_CalypsoToolkit::init(sstr.str()); 
	}
	catch(...)
	{}
}
__declspec(dllexport) void WINAPI Local_ARM_CalypsoRecConnect(void)
{
	try {
	
	char * user = getenv("USERNAME") ;
	std::stringstream sstr ;
	sstr<<"-user "<<user<<" -password "<<user<<" -env calypso_rec" ; 
	ARM_CalypsoToolkit::init(sstr.str()); 
	}
	catch(...)
	{}
}
__declspec(dllexport) void WINAPI Local_ARM_ProdConnect_WithFallBack(void)
{
	if (GetDataRetrieverVersion() == WSETKRETRIEVER)
	{
		connection_wsetoolkit("PROD","OTC","wsotc");
	}
	else
	{
		switchToETK(1);

		connection_etoolkit(SUMMIT_PROD_CONNEXION_USERNAME,
							SUMMIT_PROD_CONNEXION_PASSWD,
							SUMMIT_PROD_CONNEXION_CONTEXT,
							SUMMIT_PROD_IT_CONFIG_DOMAINSDIR,
							SUMMIT_PROD_IT_DOMAIN_NAME);
	}
}


__declspec(dllexport) void WINAPI Local_ARM_RepliConnect(void)
{
	switchToETK();

	connection_etoolkit(SUMMIT_REPLI_CONNEXION_USERNAME,
						SUMMIT_REPLI_CONNEXION_PASSWD,
						SUMMIT_REPLI_CONNEXION_CONTEXT,
						SUMMIT_REPLI_IT_CONFIG_DOMAINSDIR,
						SUMMIT_REPLI_IT_DOMAIN_NAME);
}


__declspec(dllexport) void WINAPI Local_ARM_InfocConnect(void)
{
	switchToETK();

	connection_etoolkit(SUMMIT_INFOC_CONNEXION_USERNAME,
						SUMMIT_INFOC_CONNEXION_PASSWD,
						SUMMIT_INFOC_CONNEXION_CONTEXT,
						SUMMIT_INFOC_IT_CONFIG_DOMAINSDIR,
						SUMMIT_INFOC_IT_DOMAIN_NAME);
}

__declspec(dllexport) void WINAPI Local_ARM_RecConnect(void)
{
	if (GetDataRetrieverVersion() == WSETKRETRIEVER)
	{
		connection_wsetoolkit("REC","OTC","wsotc");
	}
	else
	{
		switchToETK();

		connection_etoolkit(SUMMIT_REC342_CONNEXION_USERNAME,
							SUMMIT_REC342_CONNEXION_PASSWD,
							SUMMIT_REC342_CONNEXION_CONTEXT,
							SUMMIT_REC342_IT_CONFIG_DOMAINSDIR,
							SUMMIT_REC342_IT_DOMAIN_NAME);
	}
}

__declspec(dllexport) void WINAPI Local_ARM_ShutDownETK(void)
{
	shutdown_etoolkit();
}

__declspec(dllexport) void WINAPI Local_ARM_SwitchToETK(void)
{
	switchToETK();
}

__declspec(dllexport) void WINAPI Local_ARM_SwitchToETK_WithFallBack(void)
{
	switchToETK(1);
}

__declspec(dllexport) void WINAPI Local_ARM_SwitchToFLATFILE(void)
{
	switchToFLATFILE();
}

__declspec(dllexport) void WINAPI Local_ARM_SwitchToWSETK(void)
{
	switchToWSETK();
}

__declspec(dllexport) LPXLOPER WINAPI Local_ARM_GetMeanRevFromSummit(LPXLOPER XL_ccy,
																	 LPXLOPER XL_index,
																	 LPXLOPER XL_cvname,
																	 LPXLOPER XL_date,
																	 LPXLOPER XL_2or3Factor)
{
	ADD_LOG("Local_ARM_GetMeanRevFromSummit");
	long retCode;
	ARM_result C_result;

	// return
	static XLOPER XL_result;

	//	ARM_BEGIN();

	// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();

	// C variable
	CCString C_ccy;
	CCString C_index;
	CCString C_cvname;
	double C_date;
	CCString C_2or3Factor;

	// error
	static int error;
	static char* reason = "";

	XL_readStrCell(XL_ccy,C_ccy," ARM_ERR: currency: string expected",C_result);
	XL_readStrCell(XL_index,C_index," ARM_ERR: index: string expected",C_result);
	XL_readStrCell(XL_cvname,C_cvname," ARM_ERR: curve name: string expected",C_result);
	XL_readNumCell(XL_date,C_date," ARM_ERR: date: numeric expected",C_result);
	XL_readStrCellWD(XL_2or3Factor,C_2or3Factor,"3F"," ARM_ERR: currency: string expected",C_result);

	C_ccy.toUpper ();
	C_index.toUpper ();
	C_cvname.toUpper ();
	C_2or3Factor.toUpper ();

	retCode = ARMLOCAL_ARM_GetMeanRevFromSummit (C_ccy, C_index, C_cvname, C_date, C_2or3Factor, C_result);

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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_ARM_GetMeanRevFromSummit" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


__declspec(dllexport) LPXLOPER WINAPI Local_ARM_GetCutOffFromSummit(LPXLOPER XL_ccy,
																	LPXLOPER XL_index,
																	LPXLOPER XL_cvname,
																	LPXLOPER XL_numfac,
																	LPXLOPER XL_date)
{
	ADD_LOG("Local_ARM_GetCutOffFromSummit");
	long retCode;
	ARM_result C_result;

	// return
	static XLOPER XL_result;

	//	ARM_BEGIN();

	// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();

	// C variable
	CCString C_ccy;
	CCString C_index;
	CCString C_cvname;
	CCString C_numfac;
	double C_date;

	// error
	static int error;
	static char* reason = "";

	XL_readStrCell(XL_ccy,C_ccy," ARM_ERR: currency: string expected",C_result);
	XL_readStrCell(XL_index,C_index," ARM_ERR: index: string expected",C_result);
	XL_readStrCell(XL_cvname,C_cvname," ARM_ERR: curve name: string expected",C_result);
	XL_readStrCell(XL_numfac,C_numfac," ARM_ERR: factor number: string expected",C_result);
	XL_readNumCell(XL_date,C_date," ARM_ERR: date: numeric expected",C_result);

	C_ccy.toUpper ();
	C_index.toUpper ();
	C_cvname.toUpper ();

	retCode = ARMLOCAL_ARM_GetCutOffFromSummit (C_ccy, C_index, C_cvname, C_numfac, C_date, C_result);

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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_ARM_GetCutOffFromSummit" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}



__declspec(dllexport) LPXLOPER WINAPI Local_ARM_Hedge(LPXLOPER XL_sec,
													  LPXLOPER XL_typeHedge)
{
	ADD_LOG("Local_ARM_Hedge");
	long retCode;
	ARM_result C_result;
	LPXLOPER pxArray;

	// return
	static XLOPER XL_result;

	//	ARM_BEGIN();

	// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();

	// C variable
	CCString C_sec;
	CCString C_typeHedge;
	long hedgeId;

	// error
	static int error;
	static char* reason = "";

	XL_readStrCell(XL_sec,C_sec," ARM_ERR: security id: object expected",C_result);
	XL_readStrCell(XL_typeHedge,C_typeHedge," ARM_ERR: hedge type: string expected",C_result);

	if((hedgeId = ARM_Convhedge (C_typeHedge, C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	retCode = ARMLOCAL_HEDGE (LocalGetNumObjectId (C_sec),
							  hedgeId,
							  C_result);

	if(retCode == ARM_OK)
	{
		int nbrows = C_result.getLong();
		int nbcolumns = (int) (C_result.getDouble());

		XL_result.xltype = xltypeMulti;
		XL_result.val.array.columns = nbcolumns;
		XL_result.val.array.rows = nbrows; 
		XL_result.val.array.lparray = pxArray = (LPXLOPER)GlobalAlloc (GMEM_ZEROINIT, nbrows * nbcolumns * sizeof (XLOPER));

		for (int i = 0; i < nbrows; i++)
		{
			for (int j = 0; j < nbcolumns; j++)
			{
				pxArray[XL_Coordonnate2Rank (i, j, nbcolumns)].xltype = xltypeNum;
				pxArray[XL_Coordonnate2Rank (i, j, nbcolumns)].val.num = C_result.getArray(i+j*nbrows);
			}
		}

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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_ARM_Hedge" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}




__declspec(dllexport) LPXLOPER WINAPI Local_ARM_GETINFOFROMPRCS(LPXLOPER XL_prcs,
																LPXLOPER XL_datatype)
{
	ADD_LOG("Local_ARM_GETINFOFROMPRCS");
	long retCode;
	ARM_result C_result;

	// return
	static XLOPER XL_result;

	//	ARM_BEGIN();

	// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	
        ARM_NOCALCIFWIZ();

	    // C variables
	    CCString C_prcs;
	    CCString C_datatype;

	    // error
	    static int error;
	    static char* reason = "";

	    XL_readStrCell(XL_prcs,C_prcs," ARM_ERR: prcs id: object expected",C_result);
	    XL_readStrCell(XL_datatype,C_datatype," ARM_ERR: data type: string expected",C_result);

	    retCode = ARMLOCAL_GETINFOFROMPRCS (LocalGetNumObjectId (C_prcs),
										    C_datatype,
										    C_result);

	    if ( retCode == ARM_OK )
	    {
		    if (( strcmp((const char *) C_datatype, "FXNUMNOT") == 0 ) 
                || 
                ( strcmp((const char *) C_datatype, "FUNDINGNOT") == 0 )
                ||
                ( strcmp((const char *) C_datatype, "FUNDINGPV") == 0 )
                ||
                ( strcmp((const char *) C_datatype, "FUNDPV") == 0 ) 
               )
		    {
			    FreeCurCellErr();

			    XL_result.xltype  = xltypeNum;
			    XL_result.val.num = C_result.getDouble ();
		    }
		    else
		    {
			    FreeCurCellErr();

			    XL_result.xltype  = xltypeStr;
			    XL_result.val.str = XL_StrC2StrPascal(C_result.getString());
			    XL_result.xltype |= xlbitDLLFree;
		    }
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_ARM_GETINFOFROMPRCS" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}



__declspec(dllexport) LPXLOPER WINAPI Local_ARM_GETOBJINFOFROMPRCS(LPXLOPER XL_prcs,
																   LPXLOPER XL_datatype)
{
	ADD_LOG("Local_ARM_GETOBJINFOFROMPRCS");
	long retCode;
	ARM_result C_result;

	// return
	static XLOPER XL_result;

	//	ARM_BEGIN();

	// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();

	// C variable
	CCString C_prcs;
	CCString C_datatype;

	// error
	static int error;
	static char* reason = "";

	XL_readStrCell(XL_prcs,C_prcs," ARM_ERR: prcs id: object expected",C_result);
	XL_readStrCell(XL_datatype,C_datatype," ARM_ERR: data type: string expected",C_result);

	C_datatype.toUpper();

	long objId;
	CCString prevClass;

	CCString curClass;

	if ( (strcmp(C_datatype,"FXNUMCPN") == 0)
		|| (strcmp(C_datatype,"FXUNDCPN") == 0)
		|| (strcmp(C_datatype,"FX0") == 0)
		|| (strcmp(C_datatype,"CAP") == 0)
		|| (strcmp(C_datatype,"FLOOR") == 0)
		|| (strcmp(C_datatype,"FUNDSPREAD") == 0)
		)
	{
		curClass = LOCAL_REFVAL_CLASS;
	}
	else
	{
		curClass = LOCAL_SWAPLEG_CLASS;
	}

	CCString stringId = GetLastCurCellEnvValue ();

	if(!stringId)
	{
		retCode = ARMLOCAL_GETOBJINFOFROMPRCS (LocalGetNumObjectId (C_prcs),
											   C_datatype,
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
			retCode = ARMLOCAL_GETOBJINFOFROMPRCS (LocalGetNumObjectId (C_prcs),
												   C_datatype,
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

			retCode = ARMLOCAL_GETOBJINFOFROMPRCS (LocalGetNumObjectId (C_prcs),
												   C_datatype,
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_ARM_GETOBJINFOFROMPRCS" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


__declspec(dllexport) LPXLOPER WINAPI Local_ARM_GetFixing(LPXLOPER XL_source,
														  LPXLOPER XL_index,
														  LPXLOPER XL_term,
														  LPXLOPER XL_ccy,
														  LPXLOPER XL_date)
{
	ADD_LOG("Local_ARM_GetFixing");
	long retCode;
	ARM_result C_result;

	// return
	static XLOPER XL_result;

	//	ARM_BEGIN();

	// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();

	// C variable
	CCString C_source;
	CCString C_index;
	CCString C_term;
	CCString C_ccy;

	double C_date;

	// error
	static int error;
	static char* reason = "";

	XL_readStrCell(XL_source,C_source," ARM_ERR: source: string expected",C_result);
	XL_readStrCell(XL_index,C_index," ARM_ERR: index: string expected",C_result);
	XL_readStrCell(XL_term,C_term," ARM_ERR: term: string expected",C_result);
	XL_readStrCell(XL_ccy,C_ccy," ARM_ERR: currency: string expected",C_result);
	XL_readNumCell(XL_date,C_date," ARM_ERR: date: numeric expected",C_result);
	
	retCode = ARMLOCAL_ARM_GetFixing (C_source,
									  C_index,
									  C_term,
									  C_ccy,
									  C_date,
									  C_result);

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

	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_ARM_GetFixing" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


__declspec(dllexport) LPXLOPER WINAPI Local_ARM_GetFixingFromCalypso(
														  LPXLOPER XL_index,
														  LPXLOPER XL_term,
														  LPXLOPER XL_ccy,
														  LPXLOPER XL_src,
														  LPXLOPER XL_cvName,
														  LPXLOPER XL_date)
{
	ADD_LOG("Local_ARM_GetFixingFromCalypso");
	static XLOPER XL_result ;
	try
	{
		std::string index; ExcelTools::convert(XL_index,index); 
		std::string term; ExcelTools::convert(XL_term,term); 
		std::string ccy; ExcelTools::convert(XL_ccy,ccy); 
		std::string src; ExcelTools::convert(XL_src,src); 
		std::string cvname; ExcelTools::convert(XL_cvName,"MO",cvname); 
		ARM_Date date; ExcelTools::convert(XL_date,date); 
		double fixing;
		ARM_CalypsoToolkit::GetFixing(index,term,ccy,src,cvname,date,fixing);
		ExcelTools::convert(fixing,&XL_result); 
		return &XL_result; 
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
	return &XL_result;
}



class setSecDataFunc : public ARMResultLong2LongFunc
{
public:
	setSecDataFunc(long secId,const VECTOR<double>& data) : C_secId(secId),C_data(data) {};
	
	long operator()( ARM_result& result, long objId )
        { return ARMLOCAL_SetSecurityData(C_secId,C_data,result,objId); }

private:
	long                C_secId;
	VECTOR< double >    C_data;
};

LPXLOPER Local_SetSecurityData_Common(
	LPXLOPER XL_secId,
	LPXLOPER XL_data,
	bool PersistentInXL )
{
	
	/// this is defined first because it is used in XL macros
	static XLOPER XL_result;
	
	/// Get the variables from the XLOper variables
	ARM_result C_result;
	
	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
		/// to avoid computation if called by the wizard
		ARM_NOCALCIFWIZ();
		
		/// this is used by macros 
		/// and therefore this has to be defined
		static int error;
		static char* reason = "";

		CCString secId;
		XL_readStrCell(XL_secId,secId," ARM_ERR: Security : Object expected",C_result);
	    CCString prevSecClass = LocalGetStringObjectClass (secId);

		VECTOR<double> data;
	    XL_readNumVector(XL_data,data," ARM_ERR: datas : array of numeric expected",C_result);

		setSecDataFunc ourFunc(LocalGetNumObjectId (secId),data);
		fillXL_Result( prevSecClass, ourFunc, C_result, XL_result, PersistentInXL );
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_ARM_GETOBJINFOFROMPRCS" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}

__declspec(dllexport) LPXLOPER WINAPI Local_ARM_SetSecurityData(
	LPXLOPER XL_secId,
	LPXLOPER XL_data)
{
	ADD_LOG("Local_ARM_SetSecurityData");
	bool PersistentInXL = true;
	return Local_SetSecurityData_Common(
	    XL_secId,
	    XL_data,
        PersistentInXL );
}


__declspec(dllexport) LPXLOPER WINAPI Local_PXL_ARM_SetSecurityData(
	LPXLOPER XL_secId,
	LPXLOPER XL_data)
{
	ADD_LOG("Local_PXL_ARM_SetSecurityData");
	bool PersistentInXL = false;
	return Local_SetSecurityData_Common(
	    XL_secId,
	    XL_data,
        PersistentInXL );
}


class getSecDataFunc : public ARMResultLong2LongFunc
{
public:
	getSecDataFunc(long secId,const CCString& data) : C_secId(secId),C_data(data) {};
	
	long operator()( ARM_result& result, long Idobj)
        { return ARMLOCAL_GetSecurityData(C_secId,C_data,result); }

private:
	long                C_secId;
	 const CCString&    C_data;
};


LPXLOPER Local_GetSecurityData_Common(
	LPXLOPER XL_secId,
	LPXLOPER XL_data,
	bool PersistentInXL )
{
	
	/// this is defined first because it is used in XL macros
	static XLOPER XL_result;
	
	/// Get the variables from the XLOper variables
	ARM_result C_result;
	
	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
		/// to avoid computation if called by the wizard
		ARM_NOCALCIFWIZ();
		
		/// this is used by macros 
		/// and therefore this has to be defined
		static int error;
		static char* reason = "";

		CCString secId;
		XL_readStrCell(XL_secId,secId," ARM_ERR: Security : Object expected",C_result);
	    CCString prevSecClass = LocalGetStringObjectClass (secId);

		CCString strdata;
	    XL_readStrCell(XL_data,strdata," ARM_ERR: data : string expected",C_result);

		getSecDataFunc ourFunc(LocalGetNumObjectId (secId),strdata);
        long retCode = ourFunc(C_result, LocalGetNumObjectId (secId));

        /// return the result as an LPXLOPER
	    if (retCode == ARM_OK)
	    {
		    FreeCurCellErr ();
		    XL_result.xltype = xltypeNum;
		    XL_result.val.num = C_result.getDouble();
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_GetSecurityData_Common" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}

__declspec(dllexport) LPXLOPER WINAPI Local_ARM_GetSecurityData(
	LPXLOPER XL_secId,
	LPXLOPER XL_data)
{
	ADD_LOG("Local_ARM_GetSecurityData");
	bool PersistentInXL = false;
	return Local_GetSecurityData_Common(
	    XL_secId,
	    XL_data,
        PersistentInXL );
}


__declspec(dllexport) LPXLOPER WINAPI Local_ARM_SecurityFlows(LPXLOPER XL_labels,
															  LPXLOPER XL_values)
{
	ADD_LOG("Local_ARM_SecurityFlows");
	static XLOPER	XL_result;
	ARM_result		C_result;

	ARM_XL_TRY_BLOCK_BEGIN
	{
	    ARM_NOCALCIFWIZ();
	    
		vector<CCString>	C_labels;
		vector<double>		C_values;

	    static int		error;
	    static char*	reason = "";

	    XL_readStrVector(XL_labels, C_labels, " ARM_ERR: labels list: vector of string expected", XL_TYPE_STRING, C_result);
		XL_readNumVector(XL_values, C_values, " ARM_ERR: values : array of double expected", C_result);

		if( C_values.size() % C_labels.size() )
		{
			C_result.setMsg("size of values not a multiple of size of labels");
			ARM_ERR();
		}
		else
		{
			long	retCode;
			long	objId;

			CCString	prevClass;
			CCString	curClass = LOCAL_SECURITY_FLOWS_CLASS;
			CCString	stringId = GetLastCurCellEnvValue ();
			
			if(!stringId)
			{
				retCode = ARMLOCAL_SecurityFlows(C_labels, C_values, C_result);

				if( retCode == ARM_OK )
				{
					objId = C_result.getLong();

					LocalSetCurCellEnvValue(curClass, objId); 

					stringId = LocalMakeObjectId(objId, curClass);
				}
			}
			else
			{
				prevClass = LocalGetStringObjectClass(stringId);
				
				objId = LocalGetNumObjectId(stringId);
					
				if( curClass == prevClass )
				{
					retCode = ARMLOCAL_SecurityFlows(C_labels, C_values, C_result, objId);

					if ( retCode == ARM_OK )
					{
						LocalSetCurCellEnvValue (curClass, objId); 

						stringId = LocalMakeObjectId (objId, curClass);
					}
				}
				else
				{
					FreeCurCellContent ();
					retCode = ARMLOCAL_SecurityFlows(C_labels, C_values, C_result);

					if( retCode == ARM_OK )
					{
						objId = C_result.getLong ();

						LocalSetCurCellEnvValue (curClass, objId); 

						stringId = LocalMakeObjectId (objId, curClass);
					}
				}
			}

			if( retCode == ARM_OK )
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
	}

	ARM_XL_TRY_BLOCK_END

	ARM_XL_CATCH_ARM_EXPT

	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_ARM_SecurityFlows" )

	return	(LPXLOPER)&XL_result;
}

///////////////////////////////////////
// Get Warning when variance squeeze for local model
///////////////////////////////////////

LPXLOPER Local_ARM_ViewCell_Common(
	LPXLOPER XL_ObjectId)
{
	
	/// this is defined first because it is used in XL macros
	static XLOPER XL_result;
	
	/// Get the variables from the XLOper variables
	ARM_result C_result;
	
	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
		/// to avoid computation if called by the wizard
		ARM_NOCALCIFWIZ();
		
		/// this is used by macros 
		/// and therefore this has to be defined
		static int error;
		static char* reason = "";

		long C_ObjectId;
		XL_GETOBJID( XL_ObjectId, C_ObjectId, " ARM_ERR: Object expected",C_result );

		exportFunc1Arg<long> ourFunc(C_ObjectId,ARMLOCAL_ViewCell);

		/// Simple updating of the calculator
		long retCode = ourFunc(C_result,ARM_NULL_OBJECT_ID);

		/// feed the LPXLOPER object result 
		if (retCode == ARM_OK)
		{
			FreeCurCellErr ();
			XL_result.xltype = xltypeStr;
			XL_result.val.str = XL_StrC2StrPascal (C_result.getString() );
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_ViewCell_Common" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}
							 


///////////////////////////////////
/// version that takes into account 
/// previous creation of object
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_ARM_ViewCell(
	LPXLOPER XL_ObjectId)
{
	ADD_LOG("Local_ARM_ViewCell");
	return Local_ARM_ViewCell_Common(
	    XL_ObjectId);
}


///////////////////////////////////////
// Get Warning when variance squeeze for local model
///////////////////////////////////////

LPXLOPER Local_ARM_MatrixVectorViewer_Common(
	LPXLOPER XL_ObjectId)
{
	
	/// this is defined first because it is used in XL macros
	static XLOPER XL_result;
	
	/// Get the variables from the XLOper variables
	ARM_result C_result;
	
	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
		/// to avoid computation if called by the wizard
		ARM_NOCALCIFWIZ();
		
		/// this is used by macros 
		/// and therefore this has to be defined
		static int error;
		static char* reason = "";

		long C_ObjectId;
		XL_GETOBJID( XL_ObjectId, C_ObjectId, " ARM_ERR: Object expected",C_result );

		/// a function with a context
		vector<double> C_DataResult;
		long C_rows,C_cols;

		long retCode = ARMLOCAL_MatrixVectorViewer(C_ObjectId, C_rows, C_cols, C_DataResult,C_result);


		/// feed the LPXLOPER object result 
		if (retCode == ARM_OK)
		{
			/// add these additional lines 
			/// to display blank lines
			const int additionalLinesNb = 100;
			bool fillWithBlank = true;
			FreeCurCellContent ();
			XL_writeNumMatrixSizeWithOptions( XL_result, C_DataResult, C_rows, C_cols, " ARM_ERR: Could not get result data for Vector/Matrix info ", C_result, additionalLinesNb, fillWithBlank );
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_MatrixVectorViewer_Common" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}
							 


///////////////////////////////////
/// version that takes into account 
/// previous creation of object
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_ARM_MatrixVectorViewer(
	LPXLOPER XL_ObjectId)
{
	ADD_LOG("Local_ARM_MatrixVectorViewer");
	return Local_ARM_MatrixVectorViewer_Common(
	    XL_ObjectId);
}


__declspec(dllexport) LPXLOPER WINAPI Local_GetFixingsFromInstrument(LPXLOPER XL_tradeId)
{
	ADD_LOG("Local_GetFixingsFromInstrument");

	static XLOPER XL_result;
	ARM_result C_result;

	ARM_XL_TRY_BLOCK_BEGIN
	{
		ARM_NOCALCIFWIZ();

		CCString C_tradeId;
		// error
		static int error;
		static char* reason = "";

		XL_readStrCell(XL_tradeId,C_tradeId," ARM_ERR: Summit trade ID expected",C_result);

		long retCode;
		long objId;
		CCString prevClass;

		CCString curClass = LOCAL_REFVAL_CLASS;
		CCString stringId = GetLastCurCellEnvValue ();

		if(!stringId)
		{
			retCode = ARMLOCAL_GetFixingsFromInstrument(C_tradeId,
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
				retCode = ARMLOCAL_GetFixingsFromInstrument(C_tradeId,
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
				retCode = ARMLOCAL_GetFixingsFromInstrument(C_tradeId,
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

	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_GetFixingsFromInstrument" )

	return (LPXLOPER)&XL_result;
}

__declspec(dllexport) LPXLOPER WINAPI Local_PXL_GetFixingsFromInstrument(LPXLOPER XL_tradeId)
{
	ADD_LOG("Local_PXL_GetFixingsFromInstrument");

	static XLOPER XL_result;
	ARM_result C_result;

	ARM_XL_TRY_BLOCK_BEGIN
	{
		ARM_NOCALCIFWIZ();

		CCString C_tradeId;
		// error
		static int error;
		static char* reason = "";

		XL_readStrCell(XL_tradeId,C_tradeId," ARM_ERR: Summit trade ID expected",C_result);

		CCString curClass = LOCAL_REFVAL_CLASS;
		long retCode;

		retCode = ARMLOCAL_GetFixingsFromInstrument(C_tradeId,
													C_result);

		if(retCode == ARM_OK)
		{
			long objId = C_result.getLong ();
			CCString stringId = LocalMakeObjectId (objId, curClass);
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

	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_PXL_GetFixingsFromInstrument" )

	return (LPXLOPER)&XL_result;
}