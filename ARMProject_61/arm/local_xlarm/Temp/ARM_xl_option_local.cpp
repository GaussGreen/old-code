#pragma warning(disable : 4786)
#pragma warning(disable : 4005)

#include <libCCxll\CCxll.h>

#include <ARM\libarm_local\ARM_local_option.h>
#include <ARM\libarm_local\ARM_local_class.h>
#include <ARM\libarm_local\ARM_local_glob.h>

#include "ARM_local_interglob.h"
#include "ARM_local_interface.h"

#include "ARM_xl_gp_fctorhelper.h"

#include "ARM_xl_trycatch_local.h"
#include <util\fromto.h>

#include "util\tech_macro.h"

__declspec(dllexport) LPXLOPER WINAPI Local_EXOPTION (LPXLOPER XL_underlyingId,
													  LPXLOPER XL_optionType,
													  LPXLOPER XL_styleId,
													  LPXLOPER XL_KRefValId)
{
	ADD_LOG("Local_EXOPTION ");
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();

	// C variable
	CCString C_underlyingId;

	CCString C_optionType;
	long optionTypeId;

	CCString C_styleId;
	
	CCString C_KRefValId;
	
	// error
	static int error;
	static char* reason = "";

	XL_readStrCell(XL_underlyingId,C_underlyingId," ARM_ERR: underlying id: object expected",C_result);
	XL_readStrCell(XL_optionType,C_optionType," ARM_ERR: option type: string expected",C_result);
	XL_readStrCell(XL_styleId,C_styleId," ARM_ERR: style id: object expected",C_result);
	XL_readStrCell(XL_KRefValId,C_KRefValId," ARM_ERR: strike reference value id: object expected",C_result);
		
	if((optionTypeId = ARM_ConvCallOrPut (C_optionType, C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	long retCode;
	long objId;
	CCString prevClass;
	
	CCString curClass = LOCAL_OPTION_CLASS;
	CCString stringId = GetLastCurCellEnvValue ();

	if(!stringId)
	{
		retCode = ARMLOCAL_EXOPTION (LocalGetNumObjectId (C_underlyingId), optionTypeId,
								LocalGetNumObjectId (C_styleId),
								LocalGetNumObjectId (C_KRefValId), C_result);
							  
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
			retCode = ARMLOCAL_EXOPTION (LocalGetNumObjectId (C_underlyingId), optionTypeId,
								    LocalGetNumObjectId (C_styleId),
								    LocalGetNumObjectId (C_KRefValId), C_result, objId);

			if(retCode == ARM_OK)
			{
				LocalSetCurCellEnvValue (curClass, objId); 

				stringId = LocalMakeObjectId (objId, curClass);
			}
		}
		else
		{
			FreeCurCellContent ();
			retCode = ARMLOCAL_EXOPTION (LocalGetNumObjectId (C_underlyingId), optionTypeId,
								    LocalGetNumObjectId (C_styleId),
								    LocalGetNumObjectId (C_KRefValId), C_result);
		
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_EXOPTION" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


__declspec(dllexport) LPXLOPER WINAPI Local_PXL_EXOPTION (LPXLOPER XL_underlyingId,
														  LPXLOPER XL_optionType,
														  LPXLOPER XL_styleId,
														  LPXLOPER XL_KRefValId)
{
	ADD_LOG("Local_PXL_EXOPTION ");
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();

	// C variable
	CCString C_underlyingId;

	CCString C_optionType;
	long optionTypeId;

	CCString C_styleId;
	
	CCString C_KRefValId;
	
	// error
	static int error;
	static char* reason = "";

	XL_readStrCell(XL_underlyingId,C_underlyingId," ARM_ERR: underlying id: object expected",C_result);
	XL_readStrCell(XL_optionType,C_optionType," ARM_ERR: option type: string expected",C_result);
	XL_readStrCell(XL_styleId,C_styleId," ARM_ERR: style id: object expected",C_result);
	XL_readStrCell(XL_KRefValId,C_KRefValId," ARM_ERR: strike reference value id: object expected",C_result);
		
	if((optionTypeId = ARM_ConvCallOrPut (C_optionType, C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	long retCode;
	long objId;
	
	CCString curClass = LOCAL_OPTION_CLASS;
	CCString stringId;

	retCode = ARMLOCAL_EXOPTION (LocalGetNumObjectId (C_underlyingId), optionTypeId,
							LocalGetNumObjectId (C_styleId),
							LocalGetNumObjectId (C_KRefValId), C_result);
						  
	if(retCode == ARM_OK)
	{
		objId = C_result.getLong ();
		
		stringId = LocalMakeObjectId (objId, curClass);
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_PXL_EXOPTION" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


__declspec(dllexport) LPXLOPER WINAPI Local_OPTION (LPXLOPER XL_underlyingId,
													LPXLOPER XL_maturity,
													LPXLOPER XL_strike,
													LPXLOPER XL_optionType,
													LPXLOPER XL_exerciseType,
													LPXLOPER XL_strikeType,
													LPXLOPER XL_FstXDate,
													LPXLOPER XL_PayDate)
{
	ADD_LOG("Local_OPTION ");
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();

	// C variable
	CCString C_underlyingId;

	double C_maturity;
	double C_strike;

	CCString C_optionType;
	long optionTypeId;

	CCString C_exerciseType;
	long exerciseTypeId;

	CCString C_strikeType;
	long strikeTypeId;

	double C_FstXDate;
	double C_FstXDate_default = -1;

	double C_PayDate;
	double C_PayDate_default = -1;

	// error
	static int error;
	static char* reason = "";

	XL_readStrCell(XL_underlyingId,C_underlyingId," ARM_ERR: underlying id: object expected",C_result);
	XL_readNumCell(XL_maturity,C_maturity," ARM_ERR: maturity: date expected",C_result);
	XL_readNumCell(XL_strike,C_strike," ARM_ERR: strike: numeric expected",C_result);
	XL_readStrCell(XL_optionType,C_optionType," ARM_ERR: option type: string expected",C_result);
	XL_readStrCell(XL_exerciseType,C_exerciseType," ARM_ERR: exercise type: string expected",C_result);
	XL_readStrCellWD(XL_strikeType,C_strikeType,"Y"," ARM_ERR: strike type: string expected",C_result);
	XL_readNumCellWD(XL_FstXDate,C_FstXDate,C_FstXDate_default," ARM_ERR: first exercise date: date expected",C_result);
	XL_readNumCellWD(XL_PayDate,C_PayDate,C_PayDate_default," ARM_ERR: payment date: date expected",C_result);

	if((optionTypeId = ARM_ConvCallOrPut (C_optionType, C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	if((exerciseTypeId = ARM_ConvExerciseType (C_exerciseType, C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	if((strikeTypeId = StrikeCode (C_strikeType, C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	long retCode;
	long objId;
	CCString prevClass;
	
	CCString curClass = LOCAL_OPTION_CLASS;
	CCString stringId = GetLastCurCellEnvValue ();

	if(!stringId)
	{
		retCode = ARMLOCAL_OPTION (LocalGetNumObjectId (C_underlyingId), C_maturity,
			                  C_strike, optionTypeId, exerciseTypeId, strikeTypeId,
							  C_FstXDate, C_PayDate, C_result);

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
			retCode = ARMLOCAL_OPTION (LocalGetNumObjectId (C_underlyingId), C_maturity,
			                      C_strike, optionTypeId, exerciseTypeId, strikeTypeId,
							      C_FstXDate, C_PayDate, C_result, objId);

			if(retCode == ARM_OK)
			{
				LocalSetCurCellEnvValue (curClass, objId); 

				stringId = LocalMakeObjectId (objId, curClass);
			}
		}
		else
		{
			FreeCurCellContent ();
			retCode = ARMLOCAL_OPTION (LocalGetNumObjectId (C_underlyingId), C_maturity,
			                      C_strike, optionTypeId, exerciseTypeId, strikeTypeId,
							      C_FstXDate, C_PayDate, C_result);

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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_OPTION" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}



__declspec(dllexport) LPXLOPER WINAPI Local_PXL_OPTION (LPXLOPER XL_underlyingId,
														LPXLOPER XL_maturity,
														LPXLOPER XL_strike,
														LPXLOPER XL_optionType,
														LPXLOPER XL_exerciseType,
														LPXLOPER XL_strikeType,
														LPXLOPER XL_FstXDate,
														LPXLOPER XL_PayDate)
{
	ADD_LOG("Local_PXL_OPTION ");
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();

	// C variable
	CCString C_underlyingId;

	double C_maturity;
	double C_strike;

	CCString C_optionType;
	long optionTypeId;

	CCString C_exerciseType;
	long exerciseTypeId;

	CCString C_strikeType;
	long strikeTypeId;

	double C_FstXDate;
	double C_FstXDate_default = -1;

	double C_PayDate;
	double C_PayDate_default = -1;

	// error
	static int error;
	static char* reason = "";

	XL_readStrCell(XL_underlyingId,C_underlyingId," ARM_ERR: underlying id: object expected",C_result);
	XL_readNumCell(XL_maturity,C_maturity," ARM_ERR: maturity: date expected",C_result);
	XL_readNumCell(XL_strike,C_strike," ARM_ERR: strike: numeric expected",C_result);
	XL_readStrCell(XL_optionType,C_optionType," ARM_ERR: option type: string expected",C_result);
	XL_readStrCell(XL_exerciseType,C_exerciseType," ARM_ERR: exercise type: string expected",C_result);
	XL_readStrCellWD(XL_strikeType,C_strikeType,"Y"," ARM_ERR: strike type: string expected",C_result);
	XL_readNumCellWD(XL_FstXDate,C_FstXDate,C_FstXDate_default," ARM_ERR: first exercise date: date expected",C_result);
	XL_readNumCellWD(XL_PayDate,C_PayDate,C_PayDate_default," ARM_ERR: payment date: date expected",C_result);

	if((optionTypeId = ARM_ConvCallOrPut (C_optionType, C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	if((exerciseTypeId = ARM_ConvExerciseType (C_exerciseType, C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	if((strikeTypeId = StrikeCode (C_strikeType, C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	long retCode;
	long objId;
	
	CCString curClass = LOCAL_OPTION_CLASS;
	CCString stringId;

	retCode = ARMLOCAL_OPTION (LocalGetNumObjectId (C_underlyingId), C_maturity,
			              C_strike, optionTypeId, exerciseTypeId, strikeTypeId,
						  C_FstXDate, C_PayDate, C_result);

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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_PXL_OPTION" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


__declspec(dllexport) LPXLOPER WINAPI Local_VolImp (LPXLOPER XL_secId,
													LPXLOPER XL_modId,
													LPXLOPER XL_price)
{
	ADD_LOG("Local_VolImp ");
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
	
	double C_price;
	
	// error
	static int error;
	static char* reason = "";

	XL_readStrCell(XL_secId,C_secId," ARM_ERR: security id: object expected",C_result);
	XL_readStrCell(XL_modId,C_modId," ARM_ERR: model id: object expected",C_result);
	XL_readNumCell(XL_price,C_price," ARM_ERR: price: numeric expected",C_result);

	long retCode = ARMLOCAL_VolImp (LocalGetNumObjectId (C_secId),
									LocalGetNumObjectId (C_modId),
									C_price,
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_VolImp" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


__declspec(dllexport) LPXLOPER WINAPI Local_bsOption (LPXLOPER XL_spot,
													  LPXLOPER XL_strike,
													  LPXLOPER XL_volatility,
													  LPXLOPER XL_dividend,
													  LPXLOPER XL_discountRate,
													  LPXLOPER XL_maturity,
													  LPXLOPER XL_CallPut)
{
	ADD_LOG("Local_bsOption ");
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();

	// C variable
	double C_spot;
    double C_strike;
    double C_volatility;
    double C_dividend;
    double C_discountRate;
    double C_maturity;
 	CCString C_CallPut;
	long l_CallPut;

	// error
	static int error;
	static char* reason = "";

	XL_readNumCell(XL_spot,C_spot," ARM_ERR: spot : Numeric expected ",C_result);
	XL_readNumCell(XL_strike,C_strike," ARM_ERR: strike : Numeric expected ",C_result);
	XL_readNumCell(XL_volatility,C_volatility," ARM_ERR: volatility : Numeric expected ",C_result);
	XL_readNumCell(XL_dividend,C_dividend," ARM_ERR: dividend : Numeric expected ",C_result);
	XL_readNumCell(XL_discountRate,C_discountRate," ARM_ERR: discountrate : Numeric expected ",C_result);
	XL_readNumCell(XL_maturity,C_maturity," ARM_ERR: maturity : Numeric expected ",C_result);
	XL_readStrCell(XL_CallPut,C_CallPut," ARM_ERR: callput : String expected ",C_result);

	if((l_CallPut = ARM_ConvCallOrPut (C_CallPut, C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	long retCode = ARMLOCAL_bsOption( C_spot,C_strike,C_volatility,C_dividend,C_discountRate,C_maturity,(double)l_CallPut, C_result);

	if(retCode == ARM_OK)
	{
		FreeCurCellErr ();
		XL_result.xltype = xltypeNum;
		XL_result.val.num = C_result.getDouble () * 100.;
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_bsOption" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


__declspec(dllexport) LPXLOPER WINAPI Local_bsVega (LPXLOPER XL_spot,
													LPXLOPER XL_strike,
													LPXLOPER XL_volatility,
													LPXLOPER XL_dividend,
													LPXLOPER XL_discountRate,
													LPXLOPER XL_maturity,
													LPXLOPER XL_CallPut)
{
	ADD_LOG("Local_bsVega ");
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();

	// C variable
	double C_spot;
    double C_strike;
    double C_volatility;
    double C_dividend;
    double C_discountRate;
    double C_maturity;
 	CCString C_CallPut;
	long l_CallPut;

	// error
	static int error;
	static char* reason = "";

	XL_readNumCell(XL_spot,C_spot," ARM_ERR: spot : Numeric expected ",C_result);
	XL_readNumCell(XL_strike,C_strike," ARM_ERR: strike : Numeric expected ",C_result);
	XL_readNumCell(XL_volatility,C_volatility," ARM_ERR: volatility : Numeric expected ",C_result);
	XL_readNumCell(XL_dividend,C_dividend," ARM_ERR: dividend : Numeric expected ",C_result);
	XL_readNumCell(XL_discountRate,C_discountRate," ARM_ERR: discountrate : Numeric expected ",C_result);
	XL_readNumCell(XL_maturity,C_maturity," ARM_ERR: maturity : Numeric expected ",C_result);
	XL_readStrCell(XL_CallPut,C_CallPut," ARM_ERR: callput : String expected ",C_result);

	if((l_CallPut = ARM_ConvCallOrPut (C_CallPut, C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	long retCode = ARMLOCAL_bsVega( C_spot,C_strike,C_volatility,C_dividend,C_discountRate,C_maturity,(double)l_CallPut, C_result);

	if(retCode == ARM_OK)
	{
		FreeCurCellErr ();
		XL_result.xltype = xltypeNum;
		XL_result.val.num = C_result.getDouble () * 100.;
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_bsVega" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}



__declspec(dllexport) LPXLOPER WINAPI Local_bsDelta (LPXLOPER XL_spot,
													 LPXLOPER XL_strike,
													 LPXLOPER XL_volatility,
													 LPXLOPER XL_dividend,
													 LPXLOPER XL_discountRate,
													 LPXLOPER XL_maturity,
													 LPXLOPER XL_CallPut)
{
	ADD_LOG("Local_bsDelta ");
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();

	// C variable
	double C_spot;
    double C_strike;
    double C_volatility;
    double C_dividend;
    double C_discountRate;
    double C_maturity;
 	CCString C_CallPut;
	long l_CallPut;

	// error
	static int error;
	static char* reason = "";

	XL_readNumCell(XL_spot,C_spot," ARM_ERR: spot : Numeric expected ",C_result);
	XL_readNumCell(XL_strike,C_strike," ARM_ERR: strike : Numeric expected ",C_result);
	XL_readNumCell(XL_volatility,C_volatility," ARM_ERR: volatility : Numeric expected ",C_result);
	XL_readNumCell(XL_dividend,C_dividend," ARM_ERR: dividend : Numeric expected ",C_result);
	XL_readNumCell(XL_discountRate,C_discountRate," ARM_ERR: discountrate : Numeric expected ",C_result);
	XL_readNumCell(XL_maturity,C_maturity," ARM_ERR: maturity : Numeric expected ",C_result);
	XL_readStrCell(XL_CallPut,C_CallPut," ARM_ERR: callput : String expected ",C_result);

	if((l_CallPut = ARM_ConvCallOrPut (C_CallPut, C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	long retCode = ARMLOCAL_bsDelta( C_spot,C_strike,C_volatility,C_dividend,C_discountRate,C_maturity,(double)l_CallPut, C_result);

	if(retCode == ARM_OK)
	{
		FreeCurCellErr ();
		XL_result.xltype = xltypeNum;
		XL_result.val.num = C_result.getDouble () * 100.;
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_bsDelta" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


__declspec(dllexport) LPXLOPER WINAPI Local_bsGamma (LPXLOPER XL_spot,
													 LPXLOPER XL_strike,
													 LPXLOPER XL_volatility,
													 LPXLOPER XL_dividend,
													 LPXLOPER XL_discountRate,
													 LPXLOPER XL_maturity,
													 LPXLOPER XL_CallPut)
{
	ADD_LOG("Local_bsGamma ");
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();

	// C variable
	double C_spot;
    double C_strike;
    double C_volatility;
    double C_dividend;
    double C_discountRate;
    double C_maturity;
 	CCString C_CallPut;
	long l_CallPut;

	// error
	static int error;
	static char* reason = "";

	XL_readNumCell(XL_spot,C_spot," ARM_ERR: spot : Numeric expected ",C_result);
	XL_readNumCell(XL_strike,C_strike," ARM_ERR: strike : Numeric expected ",C_result);
	XL_readNumCell(XL_volatility,C_volatility," ARM_ERR: volatility : Numeric expected ",C_result);
	XL_readNumCell(XL_dividend,C_dividend," ARM_ERR: dividend : Numeric expected ",C_result);
	XL_readNumCell(XL_discountRate,C_discountRate," ARM_ERR: discountrate : Numeric expected ",C_result);
	XL_readNumCell(XL_maturity,C_maturity," ARM_ERR: maturity : Numeric expected ",C_result);
	XL_readStrCell(XL_CallPut,C_CallPut," ARM_ERR: callput : String expected ",C_result);

	if((l_CallPut = ARM_ConvCallOrPut (C_CallPut, C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	long retCode = ARMLOCAL_bsGamma( C_spot,C_strike,C_volatility,C_dividend,C_discountRate,C_maturity,(double)l_CallPut, C_result);

	if(retCode == ARM_OK)
	{
		FreeCurCellErr ();
		XL_result.xltype = xltypeNum;
		XL_result.val.num = C_result.getDouble () * 100.;
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_bsGamma" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


__declspec(dllexport) LPXLOPER WINAPI Local_bsTheta (LPXLOPER XL_spot,
													 LPXLOPER XL_strike,
													 LPXLOPER XL_volatility,
													 LPXLOPER XL_dividend,
													 LPXLOPER XL_discountRate,
													 LPXLOPER XL_maturity,
													 LPXLOPER XL_CallPut)
{
	ADD_LOG("Local_bsTheta ");
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();

	// C variable
	double C_spot;
    double C_strike;
    double C_volatility;
    double C_dividend;
    double C_discountRate;
    double C_maturity;
 	CCString C_CallPut;
	long l_CallPut;

	// error
	static int error;
	static char* reason = "";

	XL_readNumCell(XL_spot,C_spot," ARM_ERR: spot : Numeric expected ",C_result);
	XL_readNumCell(XL_strike,C_strike," ARM_ERR: strike : Numeric expected ",C_result);
	XL_readNumCell(XL_volatility,C_volatility," ARM_ERR: volatility : Numeric expected ",C_result);
	XL_readNumCell(XL_dividend,C_dividend," ARM_ERR: dividend : Numeric expected ",C_result);
	XL_readNumCell(XL_discountRate,C_discountRate," ARM_ERR: discountrate : Numeric expected ",C_result);
	XL_readNumCell(XL_maturity,C_maturity," ARM_ERR: maturity : Numeric expected ",C_result);
	XL_readStrCell(XL_CallPut,C_CallPut," ARM_ERR: callput : String expected ",C_result);

	if((l_CallPut = ARM_ConvCallOrPut (C_CallPut, C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	long retCode = ARMLOCAL_bsTheta( C_spot,C_strike,C_volatility,C_dividend,C_discountRate,C_maturity,(double)l_CallPut, C_result);

	if(retCode == ARM_OK)
	{
		FreeCurCellErr ();
		XL_result.xltype = xltypeNum;
		XL_result.val.num = C_result.getDouble () * 100.;
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_bsTheta" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


__declspec(dllexport) LPXLOPER WINAPI Local_ARM_OPTIONONPORTFOLIO (LPXLOPER XL_portfolioId,
																   LPXLOPER XL_styleId,
																   LPXLOPER XL_strikesId,
																   LPXLOPER XL_optionType)
{
	ADD_LOG("Local_ARM_OPTIONONPORTFOLIO ");
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();

	// C variable
	CCString C_portfolioId;
	CCString C_styleId;

	CCString C_strikes;
	long C_strikesId;

	CCString C_optionType;
	long optionTypeId;

	// error
	static int error;
	static char* reason = "";

	XL_readStrCell(XL_portfolioId,C_portfolioId," ARM_ERR: portfolio id: object expected",C_result);
	XL_readStrCell(XL_styleId,C_styleId," ARM_ERR: style Id: object expected",C_result);
	XL_readStrCellWD(XL_strikesId,C_strikes,"DEFAULT"," ARM_ERR: strikes id: object expected",C_result);
	XL_readStrCellWD(XL_optionType,C_optionType,"C"," ARM_ERR: option type: string expected",C_result);
		
	if((optionTypeId = ARM_ConvCallOrPut (C_optionType, C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	if(C_strikes == "DEFAULT")
	{
		C_strikesId = ARM_NULL_OBJECT;
	}
	else
	{
		C_strikesId = LocalGetNumObjectId (C_strikes);
	}

	long retCode;
	long objId;
	CCString prevClass;
	
	CCString curClass = LOCAL_OPTION_PORTFOLIO_CLASS;
	CCString stringId = GetLastCurCellEnvValue ();

	if(!stringId)
	{
		retCode = ARMLOCAL_ARM_OPTIONPORTFOLIO (LocalGetNumObjectId (C_portfolioId),
												LocalGetNumObjectId (C_styleId),
												C_strikesId,
												optionTypeId,
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
			retCode = ARMLOCAL_ARM_OPTIONPORTFOLIO (LocalGetNumObjectId (C_portfolioId),
													LocalGetNumObjectId (C_styleId),
													C_strikesId,
													optionTypeId,
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

			retCode = ARMLOCAL_ARM_OPTIONPORTFOLIO (LocalGetNumObjectId (C_portfolioId),
													LocalGetNumObjectId (C_styleId),
													C_strikesId,
													optionTypeId,
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_ARM_OPTIONONPORTFOLIO" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}

__declspec(dllexport) LPXLOPER WINAPI Local_PXL_ARM_OPTIONONPORTFOLIO (LPXLOPER XL_portfolioId,
																	   LPXLOPER XL_styleId,
																	   LPXLOPER XL_strikesId,
																	   LPXLOPER XL_optionType)
{
	ADD_LOG("Local_PXL_ARM_OPTIONONPORTFOLIO ");
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();

	// C variable
	CCString C_portfolioId;
	CCString C_styleId;

	CCString C_strikes;
	long C_strikesId;

	CCString C_optionType;
	long optionTypeId;

	// error
	static int error;
	static char* reason = "";

	XL_readStrCell(XL_portfolioId,C_portfolioId," ARM_ERR: portfolio id: object expected",C_result);
	XL_readStrCell(XL_styleId,C_styleId," ARM_ERR: style Id: object expected",C_result);
	XL_readStrCellWD(XL_strikesId,C_strikes,"DEFAULT"," ARM_ERR: strikes id: object expected",C_result);
	XL_readStrCellWD(XL_optionType,C_optionType,"C"," ARM_ERR: option type: string expected",C_result);
		
	if((optionTypeId = ARM_ConvCallOrPut (C_optionType, C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	if(C_strikes == "DEFAULT")
	{
		C_strikesId = ARM_NULL_OBJECT;
	}
	else
	{
		C_strikesId = LocalGetNumObjectId (C_strikes);
	}

	long retCode;
	long objId;

	CCString curClass = LOCAL_OPTION_PORTFOLIO_CLASS;
	CCString stringId;

	retCode = ARMLOCAL_ARM_OPTIONPORTFOLIO (LocalGetNumObjectId (C_portfolioId),
											LocalGetNumObjectId (C_styleId),
											C_strikesId,
											optionTypeId,
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_PXL_ARM_OPTIONONPORTFOLIO" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


__declspec(dllexport) LPXLOPER WINAPI Local_ARM_GETPFFROMSFRMCALIBRATOROFCRA (LPXLOPER XL_CalibratorId,
															                  LPXLOPER XL_PortfolioType)
{
	ADD_LOG("Local_ARM_GETPFFROMSFRMCALIBRATOROFCRA ");
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	    ARM_NOCALCIFWIZ();

	    // error
	    static int error;
	    static char* reason = "";

        // C variable
        CCString C_CalibratorId;
	    XL_readStrCell(XL_CalibratorId,C_CalibratorId," ARM_ERR: SFRM Calibrator of CRA id: object expected",C_result);
        
        // C variable
	    CCString C_PortfolioType;
        XL_readStrCellWD(XL_PortfolioType,C_PortfolioType,"Volatility"," ARM_ERR: portfolio Type: string expected",C_result);

        long retCode;
	    long objId;

	    CCString prevClass;
	    CCString curClass = LOCAL_PF_CLASS;
	    CCString stringId = GetLastCurCellEnvValue ();

	    if(!stringId)
	    {
		    retCode = ARMLOCAL_ARM_GETPFFROMSFRMCALIBRATOROFCRA (LocalGetNumObjectId (C_CalibratorId),
                                               C_PortfolioType,
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
			    retCode = ARMLOCAL_ARM_GETPFFROMSFRMCALIBRATOROFCRA (LocalGetNumObjectId (C_CalibratorId),
                                                    C_PortfolioType,
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

			    retCode = ARMLOCAL_ARM_GETPFFROMSFRMCALIBRATOROFCRA (LocalGetNumObjectId (C_CalibratorId),
                                                    C_PortfolioType,
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_ARM_GETPFFROMSFRMCALIBRATOROFCRA" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}

__declspec(dllexport) LPXLOPER WINAPI Local_PXL_ARM_GETPFFROMSFRMCALIBRATOROFCRA (LPXLOPER XL_CalibratorId,
															                      LPXLOPER XL_PortfolioType)
{
	ADD_LOG("Local_PXL_ARM_GETPFFROMSFRMCALIBRATOROFCRA ");
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	    ARM_NOCALCIFWIZ();
        // error
	    static int error;
	    static char* reason = "";

	    // C variable
        CCString C_CalibratorId;
	    XL_readStrCell(XL_CalibratorId,C_CalibratorId," ARM_ERR: SFRM Calibrator of CRA id: object expected",C_result);
    
        // C variable
	    CCString C_PortfolioType;
        XL_readStrCellWD(XL_PortfolioType,C_PortfolioType,"Volatility"," ARM_ERR: portfolio Type: string expected",C_result);
	    
        long retCode;
	    long objId;

	    CCString curClass = LOCAL_PF_CLASS;
	    CCString stringId;

	    retCode = ARMLOCAL_ARM_GETPFFROMSFRMCALIBRATOROFCRA (LocalGetNumObjectId (C_CalibratorId),
                                                            C_PortfolioType,
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_PXL_ARM_GETPFFROMSFRMCALIBRATOROFCRA" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


_declspec(dllexport) LPXLOPER WINAPI Local_SpreadOptionFormula(
	LPXLOPER XL_fwd1, 
	LPXLOPER XL_fwd2, 
	LPXLOPER XL_vol1, 
	LPXLOPER XL_vol2, 
	LPXLOPER XL_Correl, 
	LPXLOPER XL_strike, 
	LPXLOPER XL_optMat, 
	LPXLOPER XL_optType, 
	LPXLOPER XL_modelType, 
	LPXLOPER XL_spreadVol )
{
	static XLOPER XL_result;
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();

	
	// error
	static int error;
	static char* reason = "";


	double C_fwd1;
	XL_readNumCell(XL_fwd1,  C_fwd1,  " ARM_ERR: fwd1: numeric expected",C_result);
	double C_fwd2;
	XL_readNumCell(XL_fwd2,  C_fwd2,  " ARM_ERR: fwd2: numeric expected",C_result);
	double C_vol1;
	XL_readNumCell(XL_vol1,  C_vol1,  " ARM_ERR: vol1: numeric expected",C_result);
	double C_vol2;
	XL_readNumCell(XL_vol2,	 C_vol2,  " ARM_ERR: vol2: numeric expected",C_result);
	double C_Correl;
	XL_readNumCell(XL_Correl,C_Correl," ARM_ERR: correl: numeric expected",C_result);
	double C_strike;
	XL_readNumCell(XL_strike,C_strike," ARM_ERR: strike: numeric expected",C_result);
	double C_optMat;
	XL_readNumCell(XL_optMat,C_optMat," ARM_ERR: option maturity: numeric expected",C_result);

	int C_optTypeId;
	CCString C_optType;
	XL_readStrCell(XL_optType,C_optType," ARM_ERR: call or put: string expected",C_result);

	if((C_optTypeId = ARM_ConvCallOrPut(C_optType, C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	CCString C_modelType;
	long C_modelTypeId;
	XL_readStrCellWD(XL_modelType, C_modelType,"2LOG"," ARM_ERR: model type: string expected",C_result);
	if ((C_modelTypeId = ARM_ConvModelType (C_modelType, C_result)) == ARM_DEFAULT_ERR)
	{
	   ARM_ARG_ERR();
	   return (LPXLOPER)&XL_result;
	}
	double C_spreadVol, spreadVolDefault=-1;
	XL_readNumCellWD(XL_spreadVol,C_spreadVol,spreadVolDefault," ARM_ERR: spread vol numeric expected",C_result);

	long retCode = ARMLOCAL_SpreadOptionFormula(
		C_fwd1, 
		C_fwd2, 
		C_vol1, 
		C_vol2, 
		C_Correl, 
		C_strike, 
		C_optMat, 
		C_optTypeId, 
		C_modelTypeId, 
		C_spreadVol,
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_SpreadOptionFormula" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}

_declspec(dllexport) LPXLOPER WINAPI Local_SumOptionCommon(
	LPXLOPER XL_dateStrip,
	LPXLOPER XL_strike,
	LPXLOPER XL_capFloor,
	LPXLOPER XL_coeff,
	LPXLOPER XL_payDate,
	LPXLOPER XL_indexDayCount,
	bool PersistentInXL)
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

		CCString C_dateStripIdStr;
		XL_readStrCell( XL_dateStrip, C_dateStripIdStr, " ARM_ERR: DateStrip Id: Object expected",C_result);
		long C_dateStripId = LocalGetNumObjectId(C_dateStripIdStr);


		double C_strike;
		XL_readNumCell( XL_strike, C_strike, " ARM_ERR: strike: numeric expected",	C_result);

		CCString capFloorStr;
		long C_capFloor;
		XL_readStrCell(XL_capFloor,capFloorStr," ARM_ERR: cap/floor: string expected",C_result);
		if((C_capFloor = ARM_ConvCapOrFloor (capFloorStr, C_result)) == ARM_DEFAULT_ERR)
		{
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}

		CCString C_coeffIdStr, C_coeffIdDefStr = "1.0";
		XL_readStrCellWD( XL_coeff, C_coeffIdStr, C_coeffIdDefStr, " ARM_ERR: DateStrip Id: Object expected",C_result);
		double C_coeff;
		long     C_coeffId;
		if(C_coeffIdStr[0] == 'L')
			C_coeffId = LocalGetNumObjectId(C_coeffIdStr);
		else
		{
			C_coeffId = ARM_NULL_OBJECT;
			C_coeff = atof(C_coeffIdStr);
		}

		double C_payDate, C_payDateDef = -1;
		XL_readNumCellWD( XL_payDate, C_payDate, C_payDateDef, " ARM_ERR: pay date: numeric expected",	C_result);

		// Except in the default case
		if (C_payDate != -1)
		{
			C_payDate = XLDateToJulian(C_payDate);
		}

		CCString C_indexDayCount, C_indexDayCountDefStr = "-1";
		long indexDayCount;
		XL_readStrCellWD( XL_indexDayCount, C_indexDayCount, C_indexDayCountDefStr, "ARM_ERR: index day count: string expected", C_result);
		indexDayCount = ARM_ConvDayCount (C_indexDayCount);

		/// use the concept of Functor to transfer the knowledge of
		/// a function with a context
		exportFunc7Args<long, double, long, double, long, double, long> 
			ourFunc(C_dateStripId, C_strike, C_capFloor, C_coeff, C_coeffId, C_payDate, indexDayCount, ARMLOCAL_SumOption);
		
		/// call the general function
		fillXL_Result( LOCAL_SUMOPT_CLASS, ourFunc, C_result, XL_result, PersistentInXL );
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_GenericCurve_Insert_Common" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}

_declspec(dllexport) LPXLOPER WINAPI Local_SumOption(
	LPXLOPER XL_dateStrip,
	LPXLOPER XL_strike,
	LPXLOPER XL_capFloor,
	LPXLOPER XL_coeff,
	LPXLOPER XL_payDate,
	LPXLOPER XL_indexDayCount)
{
	bool PersistentInXL = true;
	return Local_SumOptionCommon(
		XL_dateStrip,
		XL_strike,
		XL_capFloor,
		XL_coeff,
		XL_payDate,
		XL_indexDayCount,
		PersistentInXL);
}

_declspec(dllexport) LPXLOPER WINAPI Local_PXL_SumOption(
	LPXLOPER XL_dateStrip,
	LPXLOPER XL_strike,
	LPXLOPER XL_capFloor,
	LPXLOPER XL_coeff,
	LPXLOPER XL_payDate,
	LPXLOPER XL_indexDayCount)
{
		bool PersistentInXL = false;
		return Local_SumOptionCommon(
		XL_dateStrip,
		XL_strike,
		XL_capFloor,
		XL_coeff,
		XL_payDate,
		XL_indexDayCount,
		PersistentInXL);
}



_declspec(dllexport) LPXLOPER WINAPI Local_SmiledSwaptionCommon(
	LPXLOPER XL_swaption,
	LPXLOPER XL_data,
	bool PersistentInXL)
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

		CCString C_swaptionIdStr;
		XL_readStrCell( XL_swaption, C_swaptionIdStr, " ARM_ERR: Swaption Id: Object expected",C_result);
		long C_swaptionId = LocalGetNumObjectId(C_swaptionIdStr);

		VECTOR < double > C_data;
		XL_readNumVector(XL_data,C_data," ARM_ERR: params vector is not of a good type",C_result);

		
		/// use the concept of Functor to transfer the knowledge of
		/// a function with a context
		exportFunc2Args<long, VECTOR < double > > 
			ourFunc(C_swaptionId,C_data, ARMLOCAL_SmiledSwaption);
		
		/// call the general function
		fillXL_Result( LOCAL_SUMOPT_CLASS, ourFunc, C_result, XL_result, PersistentInXL );
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_GenericCurve_Insert_Common" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}

_declspec(dllexport) LPXLOPER WINAPI Local_SmiledSwaption(
	LPXLOPER XL_swaption,
	LPXLOPER XL_data)
{
	bool PersistentInXL = true;
	return Local_SmiledSwaptionCommon(
		XL_swaption,
		XL_data,
		PersistentInXL);
}

_declspec(dllexport) LPXLOPER WINAPI Local_PXL_SmiledSwaption(
	LPXLOPER XL_swaption,
	LPXLOPER XL_data)
{
		bool PersistentInXL = false;
		return Local_SmiledSwaptionCommon(
		XL_swaption,
		XL_data,
		PersistentInXL);
}



__declspec(dllexport) LPXLOPER WINAPI Local_STRIPOPTION (LPXLOPER XL_asOfDate,
														 LPXLOPER XL_underlyingId,
														 LPXLOPER XL_optionType,
														 LPXLOPER XL_strikesId,
														 LPXLOPER XL_scheduleId,
														 LPXLOPER XL_PorS,
														 LPXLOPER XL_fxFixingsId,
														 LPXLOPER XL_leverageId)
{
	ADD_LOG("Local_STRIPOPTION ");
	static XLOPER XL_result;
	ARM_result C_result;

	ARM_XL_TRY_BLOCK_BEGIN
	{
		ARM_NOCALCIFWIZ();

		// C variable
		double C_asOfDate;
		CCString C_underlyingId;
		CCString C_optionType;
		CCString C_strikesId;		
		CCString C_scheduleId;
		long optionTypeId;
		
		CCString C_PorS;
		CCString C_fxFixings;
		long fxFixingsId;
		CCString leverageStr;
		long leverageId;
		double leverageValue;
		double leverageValueDef = 1.0;

		// error
		static int error;
		static char* reason = "";

		XL_readNumCell( XL_asOfDate, C_asOfDate,	" ARM_ERR: as of date: date expected",	C_result );
		XL_readStrCell(XL_underlyingId,C_underlyingId," ARM_ERR: underlying id: object expected",C_result);
		XL_readStrCell(XL_optionType,C_optionType," ARM_ERR: option type: string expected",C_result);
		XL_readStrCell(XL_strikesId,C_strikesId," ARM_ERR: strikes curve id: object expected",C_result);
		XL_readStrCell(XL_scheduleId,C_scheduleId," ARM_ERR: schedule id: object expected",C_result);
		XL_readStrCellWD( XL_PorS, C_PorS, "RCV", " ARM_ERR: receive or pay : string expected",	C_result );
		XL_readStrCellWD( XL_fxFixingsId, C_fxFixings, "NULL", "ARR: fx fixings : vector expected",C_result );
		XL_readStrOrNumCellWD( XL_leverageId, leverageStr, leverageValue, leverageValueDef, leverageId, "ARR: leverage : reference value or number expected",C_result );

		C_asOfDate = XLDateToJulian(C_asOfDate);

		if((optionTypeId = ARM_ConvCallOrPut (C_optionType, C_result)) == ARM_DEFAULT_ERR)
		{
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}

		if (C_fxFixings == "NULL")
			fxFixingsId = ARM_NULL_OBJECT;
		else
			fxFixingsId = LocalGetNumObjectId(C_fxFixings);

		if (leverageId == XL_TYPE_STRING)
			leverageId = LocalGetNumObjectId(leverageStr);
		else
			leverageId = ARM_NULL_OBJECT;
		
		long retCode;
		long objId;
		CCString prevClass;
		
		CCString curClass = LOCAL_OPTION_CLASS;
		CCString stringId = GetLastCurCellEnvValue ();

		if(!stringId)
		{
			retCode = ARMLOCAL_STRIPOPTION (C_asOfDate,
											LocalGetNumObjectId (C_underlyingId),
											optionTypeId,
											LocalGetNumObjectId (C_strikesId), 
											LocalGetNumObjectId (C_scheduleId), 
											C_PorS,
											fxFixingsId,
											leverageId, 
											leverageValue,
											C_result,
											-1);
								  
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
				retCode = ARMLOCAL_STRIPOPTION (C_asOfDate,
												LocalGetNumObjectId (C_underlyingId),
												optionTypeId,
												LocalGetNumObjectId (C_strikesId), 
												LocalGetNumObjectId (C_scheduleId), 
												C_PorS,
												fxFixingsId,
												leverageId, 
												leverageValue,
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
				retCode = ARMLOCAL_STRIPOPTION (C_asOfDate,
												LocalGetNumObjectId (C_underlyingId),
												optionTypeId,
												LocalGetNumObjectId (C_strikesId), 
												LocalGetNumObjectId (C_scheduleId), 
												C_PorS,
												fxFixingsId,
												leverageId, 
												leverageValue,
												C_result,
												-1);
			
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

	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_STRIPOPTION" )

	return (LPXLOPER)&XL_result;
}



__declspec(dllexport) LPXLOPER WINAPI Local_PXL_STRIPOPTION (LPXLOPER XL_asOfDate,
															 LPXLOPER XL_underlyingId,
															 LPXLOPER XL_optionType,
															 LPXLOPER XL_strikesId,
															 LPXLOPER XL_scheduleId,
															 LPXLOPER XL_PorS,
															 LPXLOPER XL_fxFixingsId,
															 LPXLOPER XL_leverageId)
{
	ADD_LOG("Local_PXL_STRIPOPTION ");
	static XLOPER XL_result;
	ARM_result C_result;

	ARM_XL_TRY_BLOCK_BEGIN
	{
		ARM_NOCALCIFWIZ();

		// C variable
		double C_asOfDate;
		CCString C_underlyingId;
		CCString C_optionType;
		CCString C_strikesId;		
		CCString C_scheduleId;
		long optionTypeId;
		
		CCString C_PorS;
		CCString C_fxFixings;
		long fxFixingsId;

		CCString leverageStr;
		long leverageId;
		double leverageValue;
		double leverageValueDef = 1.0;

		// error
		static int error;
		static char* reason = "";

		XL_readNumCell( XL_asOfDate, C_asOfDate,	" ARM_ERR: as of date: date expected",	C_result );
		XL_readStrCell(XL_underlyingId,C_underlyingId," ARM_ERR: underlying id: object expected",C_result);
		XL_readStrCell(XL_optionType,C_optionType," ARM_ERR: option type: string expected",C_result);
		XL_readStrCell(XL_strikesId,C_strikesId," ARM_ERR: strikes curve id: object expected",C_result);
		XL_readStrCell(XL_scheduleId,C_scheduleId," ARM_ERR: schedule id: object expected",C_result);
		XL_readStrCellWD( XL_PorS, C_PorS, "RCV", " ARM_ERR: receive or pay : string expected",	C_result );
		XL_readStrCellWD( XL_fxFixingsId, C_fxFixings, "NULL", "ARR: fx fixings : vector expected",C_result );
		XL_readStrOrNumCellWD( XL_leverageId, leverageStr, leverageValue, leverageValueDef, leverageId, "ARR: leverage : reference value or number expected",C_result );

		C_asOfDate = XLDateToJulian(C_asOfDate);

		if((optionTypeId = ARM_ConvCallOrPut (C_optionType, C_result)) == ARM_DEFAULT_ERR)
		{
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}

		if (C_fxFixings == "NULL")
			fxFixingsId = ARM_NULL_OBJECT;
		else
			fxFixingsId = LocalGetNumObjectId(C_fxFixings);

		if (leverageId == XL_TYPE_STRING)
			leverageId = LocalGetNumObjectId(leverageStr);
		else
			leverageId = ARM_NULL_OBJECT;
		
		long retCode;
		CCString prevClass;
		
		CCString curClass = LOCAL_OPTION_CLASS;
		CCString stringId = GetLastCurCellEnvValue ();

		retCode = ARMLOCAL_STRIPOPTION (C_asOfDate,
										LocalGetNumObjectId (C_underlyingId),
										optionTypeId,
										LocalGetNumObjectId (C_strikesId), 
										LocalGetNumObjectId (C_scheduleId), 
										C_PorS,
										fxFixingsId,
										leverageId, 
										leverageValue,
										C_result,
										-1);
								  
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

	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_PXL_STRIPOPTION" )

	return (LPXLOPER)&XL_result;
}



__declspec(dllexport) LPXLOPER WINAPI Local_STRIPDIGITALOPTION ( LPXLOPER XL_asOfDate,
																 LPXLOPER XL_underlyingId,
																 LPXLOPER XL_optionType,
																 LPXLOPER XL_strikesId,
																 LPXLOPER XL_scheduleId,
																 LPXLOPER XL_correl,
																 LPXLOPER XL_PorS,
																 LPXLOPER XL_fxFixingsId,
																 LPXLOPER XL_digitalOptionData,
																 LPXLOPER XL_leverageId)
{
	ADD_LOG("Local_STRIPDIGITALOPTION ");
	static XLOPER XL_result;
	ARM_result C_result;

	ARM_XL_TRY_BLOCK_BEGIN
	{
		ARM_NOCALCIFWIZ();

		// C variable
		double C_asOfDate;
		CCString C_underlyingId;
		CCString C_optionType;
		CCString C_strikesId;
		CCString C_scheduleId;
		double C_correl;
		long optionTypeId;
		
		CCString C_PorS;
		CCString C_fxFixings;
		long fxFixingsId;
		vector<CCString> C_digitalOptionData;
		vector<CCString> C_digitalOptionData_def(0);
		int C_callSpreadFlag = 0;
		double C_epsilon = 0.1;
		long C_payoffCurveId = ARM_NULL_OBJECT;

		CCString leverageStr;
		long leverageId;
		double leverageValue;
		double leverageValueDef = 1.0;

		// error
		static int error;
		static char* reason = "";

		XL_readNumCell( XL_asOfDate, C_asOfDate,	" ARM_ERR: as of date: date expected",	C_result );
		XL_readStrCell(XL_underlyingId,C_underlyingId," ARM_ERR: underlying id: object expected",C_result);
		XL_readStrCell(XL_optionType,C_optionType," ARM_ERR: option type: string expected",C_result);
		XL_readStrCell(XL_strikesId,C_strikesId," ARM_ERR: strikes curve id: object expected",C_result);
		XL_readStrCell(XL_scheduleId,C_scheduleId," ARM_ERR: schedule id: object expected",C_result);
		XL_readNumCell(XL_correl,C_correl," ARM_ERR: correlation: numeric expected",C_result);
		XL_readStrCellWD( XL_PorS, C_PorS, "RCV", " ARM_ERR: receive or pay : string expected",	C_result );
		XL_readStrCellWD( XL_fxFixingsId, C_fxFixings, "NULL", "ARR: fx fixings : vector expected",C_result );
		XL_readStrVectorWD( XL_digitalOptionData, C_digitalOptionData, C_digitalOptionData_def, "ARR: digital Data : vector expected", DOUBLE_TYPE, C_result );
		XL_readStrOrNumCellWD( XL_leverageId, leverageStr, leverageValue, leverageValueDef, leverageId, "ARR: leverage : reference value or number expected",C_result );

		C_asOfDate = XLDateToJulian(C_asOfDate);

		if((optionTypeId = ARM_ConvCallOrPut (C_optionType, C_result)) == ARM_DEFAULT_ERR)
		{
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}

		if (C_fxFixings == "NULL")
			fxFixingsId = ARM_NULL_OBJECT;
		else
			fxFixingsId = LocalGetNumObjectId(C_fxFixings);

		if (leverageId == XL_TYPE_STRING)
			leverageId = LocalGetNumObjectId(leverageStr);
		else
			leverageId = ARM_NULL_OBJECT;
		
		if (C_digitalOptionData.size() > 0)
		{
			C_callSpreadFlag = atoi(C_digitalOptionData[0]);
			C_epsilon = atof(C_digitalOptionData[1]);
			C_payoffCurveId = LocalGetNumObjectId(C_digitalOptionData[2]);
		}

		long retCode;
		long objId;
		CCString prevClass;
		
		CCString curClass = LOCAL_OPTION_CLASS;
		CCString stringId = GetLastCurCellEnvValue ();

		if(!stringId)
		{
			retCode = ARMLOCAL_STRIPDIGITALOPTION (	C_asOfDate,
													LocalGetNumObjectId (C_underlyingId),
													optionTypeId,
													LocalGetNumObjectId (C_strikesId), 
													LocalGetNumObjectId (C_scheduleId),
													C_correl,
													C_PorS,
													fxFixingsId,
													C_callSpreadFlag,
													C_epsilon,
													C_payoffCurveId,
													leverageId, 
													leverageValue,
													C_result,
													-1);
								  
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
				retCode = ARMLOCAL_STRIPDIGITALOPTION (	C_asOfDate,
														LocalGetNumObjectId (C_underlyingId),
														optionTypeId,
														LocalGetNumObjectId (C_strikesId), 
														LocalGetNumObjectId (C_scheduleId), 
														C_correl,
														C_PorS,
														fxFixingsId,
														C_callSpreadFlag,
														C_epsilon,
														C_payoffCurveId,
														leverageId, 
														leverageValue,
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
				retCode = ARMLOCAL_STRIPDIGITALOPTION (	C_asOfDate,
														LocalGetNumObjectId (C_underlyingId),
														optionTypeId,
														LocalGetNumObjectId (C_strikesId), 
														LocalGetNumObjectId (C_scheduleId), 
														C_correl,
														C_PorS,
														fxFixingsId,
														C_callSpreadFlag,
														C_epsilon,
														C_payoffCurveId,
														leverageId, 
														leverageValue,
														C_result,
														-1);
			
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

	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_STRIPDIGITALOPTION" )

	return (LPXLOPER)&XL_result;
}



__declspec(dllexport) LPXLOPER WINAPI Local_PXL_STRIPDIGITALOPTION ( LPXLOPER XL_asOfDate,
																	 LPXLOPER XL_underlyingId,
																	 LPXLOPER XL_optionType,
																	 LPXLOPER XL_strikesId,
																	 LPXLOPER XL_scheduleId,
																	 LPXLOPER XL_correl,
																	 LPXLOPER XL_PorS,
																	 LPXLOPER XL_fxFixingsId,
																	 LPXLOPER XL_digitalOptionData,
																	 LPXLOPER XL_leverageId)
{
	ADD_LOG("Local_PXL_STRIPDIGITALOPTION ");
	static XLOPER XL_result;
	ARM_result C_result;

	ARM_XL_TRY_BLOCK_BEGIN
	{
		ARM_NOCALCIFWIZ();

		// C variable
		double C_asOfDate;
		CCString C_underlyingId;
		CCString C_optionType;
		CCString C_strikesId;
		CCString C_scheduleId;
		double C_correl;
		long optionTypeId;
		
		CCString C_PorS;
		CCString C_fxFixings;
		long fxFixingsId;
		vector<CCString> C_digitalOptionData;
		vector<CCString> C_digitalOptionData_def(0);
		int C_callSpreadFlag = 1;
		double C_epsilon = 0.1;
		long C_payoffCurveId = ARM_NULL_OBJECT;

		CCString leverageStr;
		long leverageId;
		double leverageValue;
		double leverageValueDef = 1.0;

		// error
		static int error;
		static char* reason = "";

		XL_readNumCell( XL_asOfDate, C_asOfDate,	" ARM_ERR: as of date: date expected",	C_result );
		XL_readStrCell(XL_underlyingId,C_underlyingId," ARM_ERR: underlying id: object expected",C_result);
		XL_readStrCell(XL_optionType,C_optionType," ARM_ERR: option type: string expected",C_result);
		XL_readStrCell(XL_strikesId,C_strikesId," ARM_ERR: strikes curve id: object expected",C_result);
		XL_readStrCell(XL_scheduleId,C_scheduleId," ARM_ERR: schedule id: object expected",C_result);
		XL_readNumCell(XL_correl,C_correl," ARM_ERR: correlation: numeric expected",C_result);
		XL_readStrCellWD( XL_PorS, C_PorS, "RCV", " ARM_ERR: receive or pay : string expected",	C_result );
		XL_readStrCellWD( XL_fxFixingsId, C_fxFixings, "NULL", "ARR: fx fixings : vector expected",C_result );
		XL_readStrVectorWD( XL_digitalOptionData, C_digitalOptionData, C_digitalOptionData_def, "ARR: digital Data : vector expected", DOUBLE_TYPE, C_result );
		XL_readStrOrNumCellWD( XL_leverageId, leverageStr, leverageValue, leverageValueDef, leverageId, "ARR: leverage : reference value or number expected",C_result );

		C_asOfDate = XLDateToJulian(C_asOfDate);
		
		if((optionTypeId = ARM_ConvCallOrPut (C_optionType, C_result)) == ARM_DEFAULT_ERR)
		{
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}

		if (C_fxFixings == "NULL")
			fxFixingsId = ARM_NULL_OBJECT;
		else
			fxFixingsId = LocalGetNumObjectId(C_fxFixings);

		if (C_digitalOptionData.size() > 0)
		{
			C_callSpreadFlag = atoi(C_digitalOptionData[0]);
			C_epsilon = atof(C_digitalOptionData[1]);
			C_payoffCurveId = LocalGetNumObjectId(C_digitalOptionData[2]);
		}

		if (leverageId == XL_TYPE_STRING)
			leverageId = LocalGetNumObjectId(leverageStr);
		else
			leverageId = ARM_NULL_OBJECT;

		long retCode;
		CCString prevClass;
		
		CCString curClass = LOCAL_OPTION_CLASS;
		CCString stringId = GetLastCurCellEnvValue ();

		retCode = ARMLOCAL_STRIPDIGITALOPTION (	C_asOfDate,
												LocalGetNumObjectId (C_underlyingId),
												optionTypeId,
												LocalGetNumObjectId (C_strikesId), 
												LocalGetNumObjectId (C_scheduleId),
												C_correl,
												C_PorS,
												fxFixingsId,
												C_callSpreadFlag,
												C_epsilon,
												C_payoffCurveId,
												leverageId, 
												leverageValue,
												C_result,
												-1);
	
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

	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_PXL_STRIPDIGITALOPTION" )

	return (LPXLOPER)&XL_result;
}




__declspec(dllexport) LPXLOPER WINAPI Local_FxOptionStrip(LPXLOPER XL_asOfDate,
														  LPXLOPER XL_underlyingId, 
														  LPXLOPER XL_strikesCurveId, 
														  LPXLOPER XL_optionType, 
														  LPXLOPER XL_startDate, 
														  LPXLOPER XL_endDate, 
														  LPXLOPER XL_notionalId,
														  LPXLOPER XL_resetData,
														  LPXLOPER XL_payData,
														  LPXLOPER XL_dayCount, 
														  LPXLOPER XL_fwdRule, 
														  LPXLOPER XL_intRule, 
														  LPXLOPER XL_stubRule,
														  LPXLOPER XL_PorS,
														  LPXLOPER XL_fxFixingsId,
														  LPXLOPER XL_digitalOptionData,
														  LPXLOPER XL_leverageId)
{
	ADD_LOG("Local_FxOptionStrip");
	static XLOPER XL_result;
	ARM_result C_result;

	ARM_XL_TRY_BLOCK_BEGIN
	{
		ARM_NOCALCIFWIZ();

		double C_asOfDate;
		CCString C_underlyingId;
		CCString C_strikesCurveId;
		CCString C_optionType;
		double C_startDate;
		double C_endDate;
		CCString C_notionalId;

		vector<CCString> C_resetData;
		vector<CCString> C_resetData_def(0);
		long C_resetFreq = K_ANNUAL;
		CCString C_resetCalendar = GETDEFAULTVALUESTR;
		double C_resetGap = GETDEFAULTVALUE;
		long C_resetTiming = K_ADVANCE; 

		vector<CCString> C_payData;
		vector<CCString> C_payData_def(0);
		CCString C_paymentCcy = GETDEFAULTVALUESTR;
		long C_payFreq = K_ANNUAL; 
		CCString C_payCalendar = GETDEFAULTVALUESTR; 
		double C_payGap = GETDEFAULTVALUE;
		long C_payTiming = K_ARREARS;

		long C_dayCount;
		long C_fwdRule;
		long C_intRule;
		long C_stubRule;
		long optionType;

		CCString C_PorS;
		CCString C_fxFixings;
		long fxFixingsId;
		vector<CCString> C_digitalOptionData;
		vector<CCString> C_digitalOptionData_def(0);
		bool isDigital = false;
		int C_callSpreadFlag = 1;
		double C_epsilon = 0.1;
		long C_payoffCurveId = ARM_NULL_OBJECT;

		CCString leverageStr;
		long leverageId;
		double leverageValue;
		double leverageValueDef = 1.0;

		static int error;
		static char* reason = "";

		XL_readNumCell( XL_asOfDate, C_asOfDate,	" ARM_ERR: as of date: date expected",	C_result );
		XL_readStrCell( XL_underlyingId, C_underlyingId,	" ARM_ERR: underlying id: object expected",	C_result );
		XL_readStrCell( XL_strikesCurveId, C_strikesCurveId," ARM_ERR: strikes curve id: object expected",C_result );
		XL_readStrCell( XL_optionType, C_optionType,		" ARM_ERR: option type: string expected",	C_result );
		XL_readNumCell( XL_startDate, C_startDate,			" ARM_ERR: start Date: date expected",		C_result );
		XL_readNumCell( XL_endDate, C_endDate,				" ARM_ERR: end Date: date expected",		C_result );
		XL_readStrCell( XL_notionalId, C_notionalId,		" ARM_ERR: notional id: object expected",	C_result );
		XL_readStrVectorWD( XL_resetData, C_resetData, C_resetData_def, "ARR: resetData : vector expected", DOUBLE_TYPE, C_result );
		XL_readStrVectorWD( XL_payData, C_payData, C_payData_def, "ARR: payData : vector expected",		DOUBLE_TYPE, C_result );
		XL_GETDAYCOUNTWD( XL_dayCount, C_dayCount,"ACTUAL",	" ARM_ERR: dayCount : string expected",		C_result );
		XL_GETFWDRULEWD( XL_fwdRule, C_fwdRule, "MF",		" ARM_ERR: fwdRule : string expected",		C_result );
		XL_GETINTRULEWD( XL_intRule, C_intRule, "ADJ",		" ARM_ERR: fwdRule : string expected",		C_result );
		XL_GETCONVRULEWD( XL_stubRule, C_stubRule, "SS",	" ARM_ERR: stub Rule : string expected",	C_result );
		XL_readStrCellWD( XL_PorS, C_PorS, "RCV",			" ARM_ERR: receive or pay : string expected",	C_result );
		XL_readStrCellWD( XL_fxFixingsId, C_fxFixings, "NULL", "ARR: fx fixings : vector expected",		C_result );
		XL_readStrVectorWD( XL_digitalOptionData, C_digitalOptionData, C_digitalOptionData_def, "ARR: digital Data : vector expected", DOUBLE_TYPE, C_result );
		XL_readStrOrNumCellWD( XL_leverageId, leverageStr, leverageValue, leverageValueDef, leverageId, "ARR: leverage : reference value or number expected",C_result );

		C_asOfDate = XLDateToJulian(C_asOfDate);

		if((optionType = ARM_ConvCallOrPut (C_optionType, C_result)) == ARM_DEFAULT_ERR)
		{
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}
		
		if (C_resetData.size() > 0)
		{
			C_resetFreq = ARM_ConvFrequency(C_resetData[0], C_result);
			
			if (C_resetData.size() > 1)
				C_resetCalendar = C_resetData[1];
			
			if (C_resetData.size() > 2)
				C_resetTiming = ARM_ConvPayResetRule(C_resetData[2]);
			
			if (C_resetData.size() > 3)
				C_resetGap = atof(C_resetData[3]);
			
			if (C_resetData.size() > 4)
				throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"reset data: 4 data max (freq, cal, timing, gap");
		}

		if (C_payData.size() > 0)
		{
			C_paymentCcy = C_payData[0];

			if (C_payData.size() > 1)
				C_payFreq = ARM_ConvFrequency(C_payData[1], C_result);
			
			if (C_payData.size() > 2)
				C_payCalendar = C_payData[2];
			
			if (C_payData.size() > 3)
				C_payTiming = ARM_ConvPayResetRule(C_payData[3]);
			
			if (C_payData.size() > 4)
				C_payGap = atof(C_payData[4]);
			
			if (C_payData.size() > 5)
				throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"pay data: 5 data max (ccy, freq, cal, timing, gap");
		}
		else
		{
			C_paymentCcy = GETDEFAULTVALUESTR;
		}

		if (C_fxFixings == "NULL")
			fxFixingsId = ARM_NULL_OBJECT;
		else
			fxFixingsId = LocalGetNumObjectId(C_fxFixings);

		if (C_digitalOptionData.size() > 0)
		{
			isDigital = (C_digitalOptionData[0] == "Y");
			if (isDigital && (C_digitalOptionData.size() > 1))
			{
				C_callSpreadFlag = atoi(C_digitalOptionData[1]);
				C_epsilon = atof(C_digitalOptionData[2]);
				C_payoffCurveId = LocalGetNumObjectId(C_digitalOptionData[3]);
			}
		}

		if (leverageId == XL_TYPE_STRING)
			leverageId = LocalGetNumObjectId(leverageStr);
		else
			leverageId = ARM_NULL_OBJECT;

		long retCode;
		long objId;
		CCString prevClass;
		
		CCString curClass = LOCAL_OPTION_CLASS;
		CCString stringId = GetLastCurCellEnvValue();

		if(!stringId)
		{
			retCode = ARMLOCAL_FxOptionStrip(C_asOfDate,
											 LocalGetNumObjectId(C_underlyingId),
											 LocalGetNumObjectId(C_strikesCurveId),
											 optionType,
											 C_startDate,
											 C_endDate,
											 LocalGetNumObjectId(C_notionalId),
											 C_paymentCcy,
											 C_resetFreq,
											 C_dayCount,
											 C_resetCalendar,
											 C_fwdRule,
											 C_intRule,
											 C_stubRule,
											 C_resetGap,
											 C_payFreq, 
											 C_payGap, 
											 C_payCalendar, 
											 C_resetTiming, 
											 C_payTiming,
											 C_PorS,
											 fxFixingsId,
											 isDigital,
											 C_callSpreadFlag,
											 C_epsilon,
											 C_payoffCurveId,
											 leverageId, 
											 leverageValue,
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
				retCode = ARMLOCAL_FxOptionStrip(C_asOfDate,
												 LocalGetNumObjectId(C_underlyingId),
												 LocalGetNumObjectId(C_strikesCurveId),
												 optionType,
												 C_startDate,
												 C_endDate,
												 LocalGetNumObjectId(C_notionalId),
												 C_paymentCcy,
												 C_resetFreq,
												 C_dayCount,
												 C_resetCalendar,
												 C_fwdRule,
												 C_intRule,
												 C_stubRule,
												 C_resetGap,
												 C_payFreq, 
												 C_payGap, 
												 C_payCalendar, 
												 C_resetTiming, 
												 C_payTiming,
												 C_PorS,
												 fxFixingsId,
												 isDigital,
												 C_callSpreadFlag,
												 C_epsilon,
												 C_payoffCurveId,
												 leverageId, 
												 leverageValue,
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
				retCode = ARMLOCAL_FxOptionStrip(C_asOfDate,
												 LocalGetNumObjectId(C_underlyingId),
												 LocalGetNumObjectId(C_strikesCurveId),
												 optionType,
												 C_startDate,
												 C_endDate,
												 LocalGetNumObjectId(C_notionalId),
												 C_paymentCcy,
												 C_resetFreq,
												 C_dayCount,
												 C_resetCalendar,
												 C_fwdRule,
												 C_intRule,
												 C_stubRule,
												 C_resetGap,
												 C_payFreq, 
												 C_payGap, 
												 C_payCalendar, 
												 C_resetTiming, 
												 C_payTiming,
												 C_PorS,
												 fxFixingsId,
												 isDigital,
												 C_callSpreadFlag,
												 C_epsilon,
												 C_payoffCurveId,
												 leverageId, 
												 leverageValue,
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

	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_FxOptionStrip" )

	return (LPXLOPER)&XL_result;
}


__declspec(dllexport) LPXLOPER WINAPI Local_PXL_FxOptionStrip(LPXLOPER XL_asOfDate,
															  LPXLOPER XL_underlyingId, 
															  LPXLOPER XL_strikesCurveId, 
															  LPXLOPER XL_optionType, 
															  LPXLOPER XL_startDate, 
															  LPXLOPER XL_endDate, 
															  LPXLOPER XL_notionalId,
															  LPXLOPER XL_resetData,
															  LPXLOPER XL_payData,
															  LPXLOPER XL_dayCount, 
															  LPXLOPER XL_fwdRule, 
															  LPXLOPER XL_intRule, 
															  LPXLOPER XL_stubRule,
															  LPXLOPER XL_PorS,
															  LPXLOPER XL_fxFixingsId,
															  LPXLOPER XL_digitalOptionData,
															  LPXLOPER XL_leverageId)
{
	ADD_LOG("Local_PXL_FxOptionStrip");
	static XLOPER XL_result;
	ARM_result C_result;

	ARM_XL_TRY_BLOCK_BEGIN
	{
		ARM_NOCALCIFWIZ();

		double C_asOfDate;
		CCString C_underlyingId;
		CCString C_strikesCurveId;
		CCString C_optionType;
		double C_startDate;
		double C_endDate;
		CCString C_notionalId;

		vector<CCString> C_resetData;
		vector<CCString> C_resetData_def(0);
		long C_resetFreq = K_ANNUAL;
		CCString C_resetCalendar = GETDEFAULTVALUESTR;
		double C_resetGap = GETDEFAULTVALUE;
		long C_resetTiming = K_ADVANCE; 

		vector<CCString> C_payData;
		vector<CCString> C_payData_def(0);
		CCString C_paymentCcy = GETDEFAULTVALUESTR;
		long C_payFreq = K_ANNUAL; 
		CCString C_payCalendar = GETDEFAULTVALUESTR; 
		double C_payGap = GETDEFAULTVALUE;
		long C_payTiming = K_ARREARS;

		long C_dayCount;
		long C_fwdRule;
		long C_intRule;
		long C_stubRule;
		long optionType;

		CCString C_PorS;
		CCString C_fxFixings;
		long fxFixingsId;
		vector<CCString> C_digitalOptionData;
		vector<CCString> C_digitalOptionData_def(0);
		bool isDigital = false;
		int C_callSpreadFlag = 1;
		double C_epsilon = 0.1;
		long C_payoffCurveId = ARM_NULL_OBJECT;

		CCString leverageStr;
		long leverageId;
		double leverageValue;
		double leverageValueDef = 1.0;

		static int error;
		static char* reason = "";

		XL_readNumCell( XL_asOfDate, C_asOfDate,	" ARM_ERR: as of date: date expected",	C_result );
		XL_readStrCell( XL_underlyingId, C_underlyingId,	" ARM_ERR: underlying id: object expected",	C_result );
		XL_readStrCell( XL_strikesCurveId, C_strikesCurveId," ARM_ERR: strikes curve id: object expected",C_result );
		XL_readStrCell( XL_optionType, C_optionType,		" ARM_ERR: option type: string expected",	C_result );
		XL_readNumCell( XL_startDate, C_startDate,			" ARM_ERR: start Date: date expected",		C_result );
		XL_readNumCell( XL_endDate, C_endDate,				" ARM_ERR: end Date: date expected",		C_result );
		XL_readStrCell( XL_notionalId, C_notionalId,		" ARM_ERR: notional id: object expected",	C_result );
		XL_readStrVectorWD( XL_resetData, C_resetData, C_resetData_def, "ARR: resetData : vector expected", DOUBLE_TYPE, C_result );
		XL_readStrVectorWD( XL_payData, C_payData, C_payData_def, "ARR: payData : vector expected",		DOUBLE_TYPE, C_result );
		XL_GETDAYCOUNTWD( XL_dayCount, C_dayCount,"ACTUAL",	" ARM_ERR: dayCount : string expected",		C_result );
		XL_GETFWDRULEWD( XL_fwdRule, C_fwdRule, "MF",		" ARM_ERR: fwdRule : string expected",		C_result );
		XL_GETINTRULEWD( XL_intRule, C_intRule, "ADJ",		" ARM_ERR: fwdRule : string expected",		C_result );
		XL_GETCONVRULEWD( XL_stubRule, C_stubRule, "SS",	" ARM_ERR: stub Rule : string expected",	C_result );
		XL_readStrCellWD( XL_PorS, C_PorS, "RCV",			" ARM_ERR: receive or pay : string expected",	C_result );
		XL_readStrCellWD( XL_fxFixingsId, C_fxFixings, "NULL", "ARR: fx fixings : vector expected",		C_result );
		XL_readStrVectorWD( XL_digitalOptionData, C_digitalOptionData, C_digitalOptionData_def, "ARR: digital Data : vector expected", DOUBLE_TYPE, C_result );
		XL_readStrOrNumCellWD( XL_leverageId, leverageStr, leverageValue, leverageValueDef, leverageId, "ARR: leverage : reference value or number expected",C_result );

		C_asOfDate = XLDateToJulian(C_asOfDate);

		if((optionType = ARM_ConvCallOrPut (C_optionType, C_result)) == ARM_DEFAULT_ERR)
		{
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}
		
		if (C_resetData.size() > 0)
		{
			C_resetFreq = ARM_ConvFrequency(C_resetData[0], C_result);
			
			if (C_resetData.size() > 1)
				C_resetCalendar = C_resetData[1];
			
			if (C_resetData.size() > 2)
				C_resetTiming = ARM_ConvPayResetRule(C_resetData[2]);
			
			if (C_resetData.size() > 3)
				C_resetGap = atof(C_resetData[3]);
			
			if (C_resetData.size() > 4)
				throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"reset data: 4 data max (freq, cal, timing, gap");
		}

		if (C_payData.size() > 0)
		{
			C_paymentCcy = C_payData[0];

			if (C_payData.size() > 1)
				C_payFreq = ARM_ConvFrequency(C_payData[1], C_result);
			
			if (C_payData.size() > 2)
				C_payCalendar = C_payData[2];
			
			if (C_payData.size() > 3)
				C_payTiming = ARM_ConvPayResetRule(C_payData[3]);
			
			if (C_payData.size() > 4)
				C_payGap = atof(C_payData[4]);
			
			if (C_payData.size() > 5)
				throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"pay data: 5 data max (ccy, freq, cal, timing, gap");
		}
		else
		{
			C_paymentCcy = GETDEFAULTVALUESTR;
		}

		if (C_fxFixings == "NULL")
			fxFixingsId = ARM_NULL_OBJECT;
		else
			fxFixingsId = LocalGetNumObjectId(C_fxFixings);

		if (C_digitalOptionData.size() > 0)
		{
			isDigital = (C_digitalOptionData[0] == "Y");
			if (isDigital && (C_digitalOptionData.size() > 1))
			{
				C_callSpreadFlag = atoi(C_digitalOptionData[1]);
				C_epsilon = atof(C_digitalOptionData[2]);
				C_payoffCurveId = LocalGetNumObjectId(C_digitalOptionData[3]);
			}
		}

		if (leverageId == XL_TYPE_STRING)
			leverageId = LocalGetNumObjectId(leverageStr);
		else
			leverageId = ARM_NULL_OBJECT;

		long retCode;
		long objId = -1;
		CCString prevClass;
		
		CCString curClass = LOCAL_OPTION_CLASS;
		CCString stringId = GetLastCurCellEnvValue();

		retCode = ARMLOCAL_FxOptionStrip(C_asOfDate,
										 LocalGetNumObjectId(C_underlyingId),
										 LocalGetNumObjectId(C_strikesCurveId),
										 optionType,
										 C_startDate,
										 C_endDate,
										 LocalGetNumObjectId(C_notionalId),
										 C_paymentCcy,
										 C_resetFreq,
										 C_dayCount,
										 C_resetCalendar,
										 C_fwdRule,
										 C_intRule,
										 C_stubRule,
										 C_resetGap,
										 C_payFreq, 
										 C_payGap, 
										 C_payCalendar, 
										 C_resetTiming, 
										 C_payTiming,
										 C_PorS,
										 fxFixingsId,
										 isDigital,
										 C_callSpreadFlag,
										 C_epsilon,
										 C_payoffCurveId,
										 leverageId, 
										 leverageValue,
										 C_result,
										 objId);
		
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

	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_PXL_FxOptionStrip" )

	return (LPXLOPER)&XL_result;
}
