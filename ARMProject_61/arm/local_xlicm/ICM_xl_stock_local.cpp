
#pragma warning(disable :4005 4786)

#include <libCCxll\CCxll.h>

#include <ARM\libarm_local\ARM_local_glob.h>
#include <ARM\libarm_local\ARM_local_class.h>

#include <ARM\libicm_local\ICM_local_stock.h>

#include "ARM_local_interface.h"
#include "ARM_local_interglob.h"




__declspec(dllexport) LPXLOPER WINAPI ARM_Credit_Stock ( LPXLOPER XL_DefLessVol,
														LPXLOPER XL_DivYield,
														LPXLOPER XL_RefModel,
														LPXLOPER XL_Spot)

{
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	ARM_NOCALCIFWIZ();

	// C variable
	double FirmVol;
	double DivYield;
	double Spot;

	CCString C_Model;

	// error
	static int error;
	static char* reason = "";


	XL_readNumCell(XL_DefLessVol,FirmVol," ARM_ERR: Defaultless volatility : numeric expected",C_result);
	XL_readNumCell(XL_DivYield,DivYield," ARM_ERR: Dividend yield: numeric expected",C_result);
	XL_readNumCell(XL_Spot,Spot," ARM_ERR: Spot: numeric expected",C_result);
	XL_readStrCell(XL_RefModel,C_Model," ARM_ERR: tree model id: Object expected",C_result);

	long retCode;

	long objId;
	CCString prevClass;
	
	CCString curClass = LOCAL_DEFAULTSTOCK_CLASS;
	CCString stringId = GetLastCurCellEnvValue ();
	
	if(!stringId)
	{	
		retCode = ICMLOCAL_STOCK(FirmVol,DivYield,
								LocalGetNumObjectId (C_Model),
								Spot,
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
		retCode = ICMLOCAL_STOCK(FirmVol,DivYield,
									LocalGetNumObjectId (C_Model),
									Spot,
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

		retCode = ICMLOCAL_STOCK(FirmVol,DivYield,
									LocalGetNumObjectId (C_Model),
									Spot,
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
	
	return (LPXLOPER)&XL_result;
}

__declspec(dllexport) LPXLOPER WINAPI ARM_Credit_StockCallOption ( LPXLOPER XL_RefStock,
														LPXLOPER XL_Strike,
														LPXLOPER XL_FirstDateIn,
														LPXLOPER XL_LastDateIn)
{
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	ARM_NOCALCIFWIZ();

	// C variable
	double StrikePrice;
	double FirstDate;
	double LastDate;

	CCString C_Stock;

	// error
	static int error;
	static char* reason = "";


	XL_readStrCell(XL_RefStock,C_Stock," ARM_ERR: underlying stock id: Object expected",C_result);
	XL_readNumCell(XL_Strike,StrikePrice," ARM_ERR: Strike price : numeric expected",C_result);
	XL_readNumCell(XL_FirstDateIn,FirstDate,"ARM_ERR: First date: numeric expected",C_result);
	XL_readNumCell(XL_LastDateIn,LastDate,"ARM_ERR: Last date: numeric expected",C_result);

	long retCode;

	long objId;
	CCString prevClass;
	
	CCString curClass = LOCAL_STOCKCALLOPTION_CLASS;
	CCString stringId = GetLastCurCellEnvValue ();
	
	if(!stringId)
	{	
		retCode = ICMLOCAL_STOCKCALLOPTION(LocalGetNumObjectId (C_Stock),
								StrikePrice,
								FirstDate,LastDate,
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
		retCode = ICMLOCAL_STOCKCALLOPTION(LocalGetNumObjectId (C_Stock),
								StrikePrice,
								FirstDate,LastDate,
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

		retCode = ICMLOCAL_STOCKCALLOPTION(LocalGetNumObjectId (C_Stock),
								StrikePrice,
								FirstDate,LastDate,
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
	
	return (LPXLOPER)&XL_result;
}
/*---- End Of File ----*/

// EOF %M% 
