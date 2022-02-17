#pragma warning(disable : 4786)
#pragma warning(disable : 4005)

#include <libCCxll\CCxll.h>

#include "ARM_local_interglob.h"
#include "ARM_local_interface.h"

#include <ARM\libarm_local\ARM_local_class.h>
#include <ARM\libarm_local\ARM_local_refval.h>

#include "ARM_xl_trycatch_local.h"

#include "util\tech_macro.h"

__declspec(dllexport) LPXLOPER WINAPI Local_REFVALUE (LPXLOPER XL_dates,
													  LPXLOPER XL_values,
													  LPXLOPER XL_valueType,
													  LPXLOPER XL_conversion,
													  LPXLOPER XL_calcMethod,
													  LPXLOPER XL_values2)
{
	ADD_LOG("Local_REFVALUE ");
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();

	// C variable
	VECTOR<double> C_dates;
	VECTOR<double> C_values;
	VECTOR<double> C_values2;
	VECTOR<double> C_values2_default;

	CCString C_valueType;
	long valueTypeId;

	double C_conversion;
	double C_conversion_default = 1;

    CCString C_calcMethod;
    long calcModId;

	// error
	static int error;
	static char* reason = "";

	XL_readNumVector(XL_dates,C_dates," ARM_ERR: dates: array of date expected",C_result);
	XL_readNumVector(XL_values,C_values," ARM_ERR: values: array of numeric expected",C_result);
	XL_readStrCellWD(XL_valueType,C_valueType,"1"," ARM_ERR: value type: string expected",C_result);
	XL_readNumCellWD(XL_conversion,C_conversion,C_conversion_default," ARM_ERR: conversion: numeric expected",C_result);
	XL_readStrCellWD(XL_calcMethod,C_calcMethod,"LIN"," ARM_ERR: calculation method: string expected",C_result);
	XL_readNumVectorWD(XL_values2,C_values2,C_values2_default," ARM_ERR: values2: array of numeric expected",C_result);

	if((valueTypeId = ARM_ConvPriceYield (C_valueType, C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	if((calcModId = ARM_ConvCalculationMethod(C_calcMethod, C_result)) == ARM_DEFAULT_ERR)
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
		retCode = ARMLOCAL_REFVALUE (C_dates, C_values, C_values2, valueTypeId, (long)C_conversion, (long) calcModId, C_result);

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
			retCode = ARMLOCAL_REFVALUE (C_dates, C_values, C_values2, valueTypeId, (long)C_conversion, (long) calcModId, C_result, objId);

			if(retCode == ARM_OK)
			{
				LocalSetCurCellEnvValue (curClass, objId); 

				stringId = LocalMakeObjectId (objId, curClass);
			}
		}
		else
		{
			FreeCurCellContent ();
			retCode = ARMLOCAL_REFVALUE (C_dates, C_values, C_values2, valueTypeId, (long)C_conversion, (long) calcModId, C_result);

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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_REFVALUE" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


__declspec(dllexport) LPXLOPER WINAPI Local_PXL_REFVALUE (LPXLOPER XL_dates,
														  LPXLOPER XL_values,
														  LPXLOPER XL_valueType,
														  LPXLOPER XL_conversion,
														  LPXLOPER XL_calcMethod,
														  LPXLOPER XL_values2)
{
	ADD_LOG("Local_PXL_REFVALUE ");
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();

	// C variable
	VECTOR<double> C_dates;
	VECTOR<double> C_values;
	VECTOR<double> C_values2;
	VECTOR<double> C_values2_default;

	CCString C_valueType;
	long valueTypeId;

	double C_conversion;
	double C_conversion_default = 1;

    CCString C_calcMethod;
    long calcModId;

	// error
	static int error;
	static char* reason = "";

	XL_readNumVector(XL_dates,C_dates," ARM_ERR: dates: array of date expected",C_result);
	XL_readNumVector(XL_values,C_values," ARM_ERR: values: array of numeric expected",C_result);
	XL_readStrCellWD(XL_valueType,C_valueType,"1"," ARM_ERR: value type: string expected",C_result);
	XL_readNumCellWD(XL_conversion,C_conversion,C_conversion_default," ARM_ERR: conversion: numeric expected",C_result);
	XL_readStrCellWD(XL_calcMethod,C_calcMethod,"LIN"," ARM_ERR: calculation method: string expected",C_result);
	XL_readNumVectorWD(XL_values2,C_values2,C_values2_default," ARM_ERR: values2: array of numeric expected",C_result);

	if((valueTypeId = ARM_ConvPriceYield (C_valueType, C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	if((calcModId = ARM_ConvCalculationMethod(C_calcMethod, C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	long retCode;
	long objId;

	CCString curClass = LOCAL_REFVAL_CLASS;
	CCString stringId;

	retCode = ARMLOCAL_REFVALUE (C_dates, C_values, C_values2, valueTypeId, (long)C_conversion, (long) calcModId, C_result);

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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_PXL_REFVALUE" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}



__declspec(dllexport) LPXLOPER WINAPI Local_CONSTREFVALUE (LPXLOPER XL_value)
{
	ADD_LOG("Local_CONSTREFVALUE ");
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();

	// C variable
	double C_value;
	
	// error
	static int error;
	static char* reason = "";

	XL_readNumCell(XL_value,C_value," ARM_ERR: value: numeric expected",C_result);

	long retCode;
	long objId;
	CCString prevClass;
	
	CCString curClass = LOCAL_REFVAL_CLASS;
	CCString stringId = GetLastCurCellEnvValue ();

	if(!stringId)
	{
		retCode = ARMLOCAL_CONSTREFVALUE (C_value, C_result);

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
			retCode = ARMLOCAL_CONSTREFVALUE (C_value, C_result, objId);

			if(retCode == ARM_OK)
			{
				LocalSetCurCellEnvValue (curClass, objId); 

				stringId = LocalMakeObjectId (objId, curClass);
			}
		}
		else
		{
			FreeCurCellContent ();
			retCode = ARMLOCAL_CONSTREFVALUE (C_value, C_result);

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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_CONSTREFVALUE" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}



__declspec(dllexport) LPXLOPER WINAPI Local_PXL_CONSTREFVALUE (LPXLOPER XL_value)
{
	ADD_LOG("Local_PXL_CONSTREFVALUE ");
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();

	// C variable
	double C_value;

	// error
	static int error;
	static char* reason = "";

	XL_readNumCell(XL_value,C_value," ARM_ERR: value: numeric expected",C_result);

	long retCode;
	long objId;

	CCString curClass = LOCAL_REFVAL_CLASS;
	CCString stringId;

	retCode = ARMLOCAL_CONSTREFVALUE (C_value, C_result);

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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_PXL_CONSTREFVALUE" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}



__declspec(dllexport) LPXLOPER WINAPI Local_IATHREELEVREFVAL (LPXLOPER XL_level1,
															  LPXLOPER XL_amort1,
															  LPXLOPER XL_level2,
															  LPXLOPER XL_amort2,
															  LPXLOPER XL_level3,
															  LPXLOPER XL_amort3,
															  LPXLOPER XL_notionnal)
{
	ADD_LOG("Local_IATHREELEVREFVAL ");
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();

	// C variable
	double C_level1;
	double C_amort1;

	double C_level2;
	double C_amort2;
	
	double C_level3;
	double C_amort3;

	double C_notionnal;
	double C_notionnal_default = 100;
			
	// error
	static int error;
	static char* reason = "";

	XL_readNumCell(XL_level1,C_level1," ARM_ERR: level1: numeric expected",C_result);
	XL_readNumCell(XL_amort1,C_amort1," ARM_ERR: amort1: numeric expected",C_result);
	XL_readNumCell(XL_level2,C_level2," ARM_ERR: level2: numeric expected",C_result);
	XL_readNumCell(XL_amort2,C_amort2," ARM_ERR: amort2: numeric expected",C_result);
	XL_readNumCell(XL_level3,C_level3," ARM_ERR: level3: numeric expected",C_result);
	XL_readNumCell(XL_amort3,C_amort3," ARM_ERR: amort3: numeric expected",C_result);
	XL_readNumCellWD(XL_notionnal,C_notionnal,C_notionnal_default," ARM_ERR: notionnal: numeric expected",C_result);

	long retCode;
	long objId;
	CCString prevClass;
	
	CCString curClass = LOCAL_IAREFVAL_CLASS;
	CCString stringId = GetLastCurCellEnvValue ();

	if(!stringId)
	{
		retCode = ARMLOCAL_IATHREELEVREFVAL (C_notionnal,
											 C_level1,
											 C_amort1,
											 C_level2,
											 C_amort2,
											 C_level3,
											 C_amort3,
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
			retCode = ARMLOCAL_IATHREELEVREFVAL (C_notionnal,
												 C_level1,
												 C_amort1,
												 C_level2,
												 C_amort2,
												 C_level3,
												 C_amort3,
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
			retCode = ARMLOCAL_IATHREELEVREFVAL (C_notionnal,
												 C_level1,
												 C_amort1,
												 C_level2,
												 C_amort2,
												 C_level3,
												 C_amort3,
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_IATHREELEVREFVAL" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}

__declspec(dllexport) LPXLOPER WINAPI Local_PXL_IATHREELEVREFVAL (LPXLOPER XL_level1,
																  LPXLOPER XL_amort1,
																  LPXLOPER XL_level2,
																  LPXLOPER XL_amort2,
																  LPXLOPER XL_level3,
																  LPXLOPER XL_amort3,
																  LPXLOPER XL_notionnal)
{
	ADD_LOG("Local_PXL_IATHREELEVREFVAL ");
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();

	// C variable
	double C_level1;
	double C_amort1;

	double C_level2;
	double C_amort2;
	
	double C_level3;
	double C_amort3;

	double C_notionnal;
	double C_notionnal_default = 100;
			
	// error
	static int error;
	static char* reason = "";

	XL_readNumCell(XL_level1,C_level1," ARM_ERR: level1: numeric expected",C_result);
	XL_readNumCell(XL_amort1,C_amort1," ARM_ERR: amort1: numeric expected",C_result);
	XL_readNumCell(XL_level2,C_level2," ARM_ERR: level2: numeric expected",C_result);
	XL_readNumCell(XL_amort2,C_amort2," ARM_ERR: amort2: numeric expected",C_result);
	XL_readNumCell(XL_level3,C_level3," ARM_ERR: level3: numeric expected",C_result);
	XL_readNumCell(XL_amort3,C_amort3," ARM_ERR: amort3: numeric expected",C_result);
	XL_readNumCellWD(XL_notionnal,C_notionnal,C_notionnal_default," ARM_ERR: notionnal: numeric expected",C_result);

	long retCode;
	long objId;
	
	CCString curClass = LOCAL_IAREFVAL_CLASS;
	CCString stringId;

	retCode = ARMLOCAL_IATHREELEVREFVAL (C_notionnal,
										 C_level1,
										 C_amort1,
										 C_level2,
										 C_amort2,
										 C_level3,
										 C_amort3,
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_PXL_IATHREELEVREFVAL" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


__declspec(dllexport) LPXLOPER WINAPI Local_ARM_CptRefValue (LPXLOPER XL_refval,
															 LPXLOPER XL_date)
{
	ADD_LOG("Local_ARM_CptRefValue ");
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();

	// C variable
	CCString C_refval;
	double C_date;

		// error
	static int error;
	static char* reason = "";

	XL_readStrCell(XL_refval,C_refval," ARM_ERR: refvalue id: object expected",C_result);
	XL_readNumCell(XL_date,C_date," ARM_ERR: maturity: numeric expected",C_result);

	long retCode;

	retCode = ARMLOCAL_CptRefValue (LocalGetNumObjectId (C_refval),
									C_date,
									C_result);

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

//	ARM_END();
}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_ARM_CptRefValue" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}



__declspec(dllexport) LPXLOPER WINAPI Local_ARM_SumRefValue (LPXLOPER XL_refval1,
															 LPXLOPER XL_refval2,
															 LPXLOPER XL_coef)
{
	ADD_LOG("Local_ARM_SumRefValue ");
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();

	// C variable
	CCString C_refval1;
	CCString C_refval2;
	long refVal2Id;

	double C_coef;
	double C_coef_default = 1.0;

		// error
	static int error;
	static char* reason = "";

	XL_readStrCell(XL_refval1,C_refval1," ARM_ERR: refvalue1 id: object expected",C_result);
	XL_readStrCellWD(XL_refval2,C_refval2,"DEFAULT"," ARM_ERR: refvalue2 id: object expected",C_result);
	XL_readNumCellWD(XL_coef,C_coef,C_coef_default," ARM_ERR: coef: numeric expected",C_result);

	if (C_refval2 == "DEFAULT")
		refVal2Id = ARM_NULL_OBJECT;
	else
		refVal2Id = LocalGetNumObjectId(C_refval2);

	long retCode;
	long objId;
	CCString prevClass;

	CCString curClass = LOCAL_REFVAL_CLASS;
	CCString stringId = GetLastCurCellEnvValue ();

	if(!stringId)
	{
		retCode = ARMLOCAL_SumRefValue (LocalGetNumObjectId(C_refval1),
										refVal2Id,
										C_coef,
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
			retCode = ARMLOCAL_SumRefValue (LocalGetNumObjectId(C_refval1),
											refVal2Id,
											C_coef,
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
			retCode = ARMLOCAL_SumRefValue (LocalGetNumObjectId(C_refval1),
											refVal2Id,
											C_coef,
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_IATHREELEVREFVAL" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}



__declspec(dllexport) LPXLOPER WINAPI Local_PXL_ARM_SumRefValue (LPXLOPER XL_refval1,
																 LPXLOPER XL_refval2,
																 LPXLOPER XL_coef)
{
	ADD_LOG("Local_PXL_ARM_SumRefValue ");
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();

	// C variable
	CCString C_refval1;
	CCString C_refval2;
	long refVal2Id;

	double C_coef;
	double C_coef_default = 1.0;

		// error
	static int error;
	static char* reason = "";

	XL_readStrCell(XL_refval1,C_refval1," ARM_ERR: refvalue1 id: object expected",C_result);
	XL_readStrCellWD(XL_refval2,C_refval2,"DEFAULT"," ARM_ERR: refvalue2 id: object expected",C_result);
	XL_readNumCellWD(XL_coef,C_coef,C_coef_default," ARM_ERR: coef: numeric expected",C_result);

	if (C_refval2 == "DEFAULT")
		refVal2Id = ARM_NULL_OBJECT;
	else
		refVal2Id = LocalGetNumObjectId(C_refval2);

	long retCode;
	long objId;

	CCString curClass = LOCAL_REFVAL_CLASS;
	CCString stringId;

	retCode = ARMLOCAL_SumRefValue (LocalGetNumObjectId(C_refval1),
									refVal2Id,
									C_coef,
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_IATHREELEVREFVAL" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


__declspec(dllexport) LPXLOPER WINAPI Local_ARM_DisplayRefValue (LPXLOPER XL_refvalue, LPXLOPER XL_isDate)
{
	ADD_LOG("Local_ARM_DisplayRefValue ");
//	ARM_BEGIN();
	
	// return
	static XLOPER XL_result;
	LPXLOPER pxArray;
	ARM_result C_result;

    /// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();

	// C variable
	CCString C_refvalue;
	double C_isDate = 1;
	double C_isDate_default = 1;

	// error
	static int error;
	static char* reason = "";

	XL_readStrCell(XL_refvalue,C_refvalue," ARM_ERR: refvalue: object expected",C_result);
	XL_readNumCellWD(XL_isDate,C_isDate,C_isDate_default," ARM_ERR: isDate: 0 or 1 expected",C_result);

	long retCode;

	retCode = ARMLOCAL_DisplayRefValue(LocalGetNumObjectId(C_refvalue), (bool) C_isDate, C_result);

	if(retCode == ARM_OK)
	{
		long nbrows = C_result.getLong();
		long nbcolumns = 2;

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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_ARM_DisplayRefValue" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}
