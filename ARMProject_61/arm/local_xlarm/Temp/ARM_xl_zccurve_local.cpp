#pragma warning(disable : 4786)
#pragma warning(disable : 4005)

#include <libCCxll\CCxll.h>

#include <ARM\libarm_local\ARM_local_zccurve.h>
#include <ARM\libarm_local\ARM_local_glob.h>
#include <ARM\libarm_local\ARM_local_ccy.h>
#include <ARM\libarm_local\ARM_local_class.h>

#include <ARM\libarm_frometk\ARM_local_etoolkit.h>

#include "ARM_local_interglob.h"
#include "ARM_local_interface.h"
#include "ARM_xl_trycatch_local.h"
#include "ARM_xl_gp_fctorhelper.h"



#include <util\fromto.h>
#include "ExcelTools.h"

#include "util\tech_macro.h"

#include "gpinflation\resetmanager.h"
using namespace ARM;




__declspec(dllexport) LPXLOPER WINAPI Local_CreateZCSwapInt (LPXLOPER XL_date,
															 LPXLOPER XL_matuRate,
															 LPXLOPER XL_mmVsFut,
															 LPXLOPER XL_swapVsFut,
															 LPXLOPER XL_raw,
															 LPXLOPER XL_interp,
															 LPXLOPER XL_ccy,
															 LPXLOPER XL_swapFrq,
															 LPXLOPER XL_fixDayCount)
{
	ADD_LOG("Local_CreateZCSwapInt ");
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();

	// C variable
	double C_date;
	VECTOR<CCString> C_matuRate;
	VECTOR<CCString> C_matu;
	VECTOR<double> C_rate; 

	CCString C_mmVsFut;
	long mmVsFutId;

	CCString C_swapVsFut;
	long swapVsFutId;

	CCString C_raw;
	long rawId;

	CCString C_interp;
	long interpId;

	CCString C_ccy;

	CCString C_swapFrq;
	long swapFrqId;

	CCString C_fixDayCount;
	long fixDayCount;

	// error
	static int error;
	static char* reason = "";

	XL_readNumCell(XL_date,C_date," ARM_ERR: as of date: date expected",C_result);
	
	XL_readStrVector(XL_matuRate,C_matuRate," ARM_ERR: maturities and rates: array of numeric expected",DOUBLE_TYPE,C_result);

	double tmp_rate = -1;
	CCString tmp_matu ("NULL");

	if(C_matuRate.size () % 2 != 0)
	{
		C_result.setMsg ("ARM_ERR: check your maturities and rates array");
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	int real_size = C_matuRate.size () / 2;
	double tmp;
	int j = 0;

	for(int i = 0; i < real_size; i++)
	{
		C_matu.push_back (C_matuRate[j]);
		if((sscanf ((const char*)C_matuRate[j+1], "%lf", &tmp) != 1))
		{
			C_result.setMsg ("ARM_ERR: check your maturities and rates array");
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}
		C_rate.push_back (tmp);
		j += 2;
	}

	CCString last = C_matu[C_matu.size () - 1];
	if((!last) || (last[0] == ' '))
	{
		C_result.setMsg ("ARM_ERR: check your maturities and rates array");
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	XL_readStrCell(XL_mmVsFut,C_mmVsFut," ARM_ERR: money market versus future: string expected",C_result);
	XL_readStrCell(XL_swapVsFut,C_swapVsFut," ARM_ERR: swap versus future: string expected",C_result);
	XL_readStrCellWD(XL_raw,C_raw,"P"," ARM_ERR: raw: string expected",C_result);
	XL_readStrCellWD(XL_interp,C_interp,"L"," ARM_ERR: interpolation mode: string expected",C_result);
	XL_readStrCellWD(XL_ccy,C_ccy,"DEFAULT"," ARM_ERR: currency: string expected",C_result);
	XL_readStrCellWD(XL_swapFrq,C_swapFrq,"-1"," ARM_ERR: swap frequency: string expected",C_result);
	XL_readStrCellWD(XL_fixDayCount,C_fixDayCount,"-1"," ARM_ERR: swap Fixed Leg DayCount: string expected",C_result); // Attn : 8 correspond à KNOBASE

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

	if((mmVsFutId = ARM_ConvMktType (C_mmVsFut, C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	if((swapVsFutId = ARM_ConvMktType (C_swapVsFut, C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	rawId = ARM_ConvCvMethod (C_raw);

	if((interpId = ARM_ConvInterpMethod (C_interp, C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	if ((swapFrqId = ARM_ConvFrequency(C_swapFrq, C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	if ((fixDayCount = ARM_ConvDayCount(C_fixDayCount)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	long retCode;
	long objId;
	CCString prevClass;
	
	CCString curClass = LOCAL_ZERO_CURVE_LIN_CLASS;
	CCString stringId = GetLastCurCellEnvValue ();
	
	if(!stringId)
	{
		retCode = ARMLOCAL_CreateZCSwapInt(C_date, C_matu, C_rate, (long)mmVsFutId, 
									   (long)swapVsFutId, (long)rawId,
									   (long)interpId, C_ccy, swapFrqId, fixDayCount,
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
			retCode = ARMLOCAL_CreateZCSwapInt(C_date, C_matu, C_rate, (long)mmVsFutId, 
									   (long)swapVsFutId, (long)rawId,
									   (long)interpId, C_ccy, swapFrqId, fixDayCount,
									   C_result,objId);

			if(retCode == ARM_OK)
			{			
				LocalSetCurCellEnvValue (curClass, objId); 

				stringId = LocalMakeObjectId (objId, curClass);
			}
		}
		else
		{
			FreeCurCellContent ();
			retCode = ARMLOCAL_CreateZCSwapInt(C_date, C_matu, C_rate, (long)mmVsFutId, 
									   (long)swapVsFutId, (long)rawId,
									   (long)interpId, C_ccy, swapFrqId, fixDayCount,
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_CreateZCSwapInt" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}



__declspec(dllexport) LPXLOPER WINAPI Local_PXL_CreateZCSwapInt (LPXLOPER XL_date,
																 LPXLOPER XL_matuRate,
																 LPXLOPER XL_mmVsFut,
																 LPXLOPER XL_swapVsFut,
																 LPXLOPER XL_raw,
																 LPXLOPER XL_interp,
																 LPXLOPER XL_ccy,
																 LPXLOPER XL_swapFrq,
															     LPXLOPER XL_fixDayCount)
{
	ADD_LOG("Local_PXL_CreateZCSwapInt ");
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();

	// C variable
	double C_date;
	VECTOR<CCString> C_matuRate;
	VECTOR<CCString> C_matu;
	VECTOR<double> C_rate; 

	CCString C_mmVsFut;
	long mmVsFutId;

	CCString C_swapVsFut;
	long swapVsFutId;

	CCString C_raw;
	long rawId;

	CCString C_interp;
	long interpId;

	CCString C_ccy;
	
	CCString C_swapFrq;
	long swapFrqId;

	CCString C_fixDayCount;
	long fixDayCount;

	// error
	static int error;
	static char* reason = "";

	XL_readNumCell(XL_date,C_date," ARM_ERR: as of date: date expected",C_result);
	
	XL_readStrVector(XL_matuRate,C_matuRate," ARM_ERR: maturities and rates: array of numeric expected",DOUBLE_TYPE,C_result);

	double tmp_rate = -1;
	CCString tmp_matu ("NULL");

	if(C_matuRate.size () % 2 != 0)
	{
		C_result.setMsg ("ARM_ERR: check your maturities and rates array");
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	int real_size = C_matuRate.size () / 2;
	double tmp;
	int j = 0;

	for(int i = 0; i < real_size; i++)
	{
		C_matu.push_back (C_matuRate[j]);
		if((sscanf ((const char*)C_matuRate[j+1], "%lf", &tmp) != 1))
		{
			C_result.setMsg ("ARM_ERR: check your maturities and rates array");
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}
		C_rate.push_back (tmp);
		j += 2;
	}

	CCString last = C_matu[C_matu.size () - 1];
	if((!last) || (last[0] == ' '))
	{
		C_result.setMsg ("ARM_ERR: check your maturities and rates array");
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	XL_readStrCell(XL_mmVsFut,C_mmVsFut," ARM_ERR: money market versus future: string expected",C_result);
	XL_readStrCell(XL_swapVsFut,C_swapVsFut," ARM_ERR: swap versus future: string expected",C_result);
	XL_readStrCellWD(XL_raw,C_raw,"P"," ARM_ERR: raw: string expected",C_result);
	XL_readStrCellWD(XL_interp,C_interp,"L"," ARM_ERR: interpolation mode: string expected",C_result);
	XL_readStrCellWD(XL_ccy,C_ccy,"DEFAULT"," ARM_ERR: currency: string expected",C_result);
	XL_readStrCellWD(XL_swapFrq,C_swapFrq,"-1"," ARM_ERR: swap frequency: string expected",C_result);
	XL_readStrCellWD(XL_fixDayCount,C_fixDayCount,"-1"," ARM_ERR: swap Fixed Leg DayCount: string expected",C_result);

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
	
	if((mmVsFutId = ARM_ConvMktType (C_mmVsFut, C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	if((swapVsFutId = ARM_ConvMktType (C_swapVsFut, C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	rawId = ARM_ConvCvMethod (C_raw);

	if((interpId = ARM_ConvInterpMethod (C_interp, C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	if ((swapFrqId = ARM_ConvFrequency(C_swapFrq, C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	if ((fixDayCount = ARM_ConvDayCount(C_fixDayCount)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	long retCode;
	long objId;
	
	CCString curClass = LOCAL_ZERO_CURVE_LIN_CLASS;
	CCString stringId;
	
	retCode = ARMLOCAL_CreateZCSwapInt (C_date, C_matu, C_rate, (long)mmVsFutId, 
								   (long)swapVsFutId, (long)rawId,
								   (long)interpId, C_ccy, swapFrqId, fixDayCount, C_result);

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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_PXL_CreateZCSwapInt" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}



__declspec(dllexport) LPXLOPER WINAPI Local_ZCFLAT(LPXLOPER XL_zeroFlat,
												   LPXLOPER XL_date,
                                                   LPXLOPER XL_ccy)
{
	ADD_LOG("Local_ZCFLAT");
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();

	// C variable
	double C_zeroFlat;
	double C_date;

    CCString C_discountCcy;

	// error
	static int error;
	static char* reason = "";

	XL_readNumCell(XL_zeroFlat,C_zeroFlat," ARM_ERR: zero flat: numeric expected",C_result);
	XL_readNumCell(XL_date,C_date," ARM_ERR: as of date: date expected",C_result);

    XL_readStrCellWD(XL_ccy, C_discountCcy, "DEFAULT",
                       " ARM_ERR: currency: string expected", C_result);

    if ( C_discountCcy == "DEFAULT" )
	{
		ARM_result currencyres;
		ARMLOCAL_ARM_GetDefaultCurrency (currencyres);
		
        if ( currencyres.getRetCode () != ARM_OK )
		{
		   ARM_ARG_ERR();

		   return (LPXLOPER)&XL_result;
		}
		else
		{
		   C_discountCcy = currencyres.getString();
		}
	}

	long retCode;
	long objId;
	CCString prevClass;
	
	CCString curClass = LOCAL_ZERO_CURVE_FLAT_CLASS;
	CCString stringId = GetLastCurCellEnvValue ();

	if(!stringId)
	{
		retCode = ARMLOCAL_zcflat(C_zeroFlat, C_date, C_discountCcy,
                                  C_result);
		
        if( retCode == ARM_OK )
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
			
		if ( curClass == prevClass )
		{
			retCode = ARMLOCAL_zcflat(C_zeroFlat, C_date,
                                      C_discountCcy,
                                      C_result, objId);

			if ( retCode == ARM_OK )
			{			
				LocalSetCurCellEnvValue (curClass, objId); 

				stringId = LocalMakeObjectId (objId, curClass);
			}
		}
		else
		{
			FreeCurCellContent();

			retCode = ARMLOCAL_zcflat(C_zeroFlat, C_date,
                                      C_discountCcy,
                                      C_result);
		
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
		FreeCurCellErr();

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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_ZCFLAT" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}



__declspec(dllexport) LPXLOPER WINAPI Local_PXL_ZCFLAT(LPXLOPER XL_zeroFlat,
													   LPXLOPER XL_date,
                                                       LPXLOPER XL_ccy)
{
	ADD_LOG("Local_PXL_ZCFLAT");
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();

	// C variable
	double C_zeroFlat;
	double C_date;

    CCString C_discountCcy;


	// error
	static int error;
	static char* reason = "";

	XL_readNumCell(XL_zeroFlat,C_zeroFlat," ARM_ERR: zero flat: numeric expected",C_result);
	XL_readNumCell(XL_date,C_date," ARM_ERR: as of date: date expected",C_result);

       XL_readStrCellWD(XL_ccy, C_discountCcy, "DEFAULT",
                       " ARM_ERR: currency: string expected", C_result);

    if ( C_discountCcy == "DEFAULT" )
	{
		ARM_result currencyres;
		ARMLOCAL_ARM_GetDefaultCurrency(currencyres);
		
        if ( currencyres.getRetCode () != ARM_OK )
		{
		   ARM_ARG_ERR();

		   return (LPXLOPER)&XL_result;
		}
		else
		{
		   C_discountCcy = currencyres.getString();
		}
	}


	long retCode;
	long objId;
	CCString prevClass;
	
	CCString curClass = LOCAL_ZERO_CURVE_FLAT_CLASS;
	CCString stringId;

	retCode = ARMLOCAL_zcflat(C_zeroFlat, C_date, C_discountCcy, C_result);

	if ( retCode == ARM_OK )
	{
	   objId = C_result.getLong ();

	   stringId = LocalMakeObjectId(objId, curClass);
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_PXL_ZCFLAT" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}




__declspec(dllexport) LPXLOPER WINAPI Local_DiscountPrice (LPXLOPER XL_curve,
														   LPXLOPER XL_matu)
{
	ADD_LOG("Local_DiscountPrice ");
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();

	// C variable
	CCString C_curve;
	double C_matu;

		// error
	static int error;
	static char* reason = "";

	XL_readStrCell(XL_curve,C_curve," ARM_ERR: curve id: object expected",C_result);
	XL_readNumCell(XL_matu,C_matu," ARM_ERR: maturity: numeric expected",C_result);

	long retCode;

	if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 1)
	{
		retCode = ARMLOCAL_DiscountPrice (LocalGetNumObjectId (C_curve), C_matu,
							  C_result);
	}
	else
		retCode = ARM_KO;

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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_DiscountPrice" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


__declspec(dllexport) LPXLOPER WINAPI Local_DiscountYield (LPXLOPER XL_curve,
														   LPXLOPER XL_matu,
														   LPXLOPER XL_meth)
{
	ADD_LOG("Local_DiscountYield ");
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();

	// C variable
	CCString C_curve;
	double C_matu;
	CCString C_strMeth;
	double C_meth;
	double C_meth_default = 0.;
	long typeMeth;

	// error
	static int error;
	static char* reason = "";

	XL_readStrCell(XL_curve,C_curve," ARM_ERR: curve id: object expected",C_result);
	XL_readNumCell(XL_matu,C_matu," ARM_ERR: maturity: numeric expected",C_result);
	XL_readStrOrNumCellWD(XL_meth,C_strMeth,C_meth,C_meth_default,typeMeth," ARM_ERR: decompounding method: string or numeric expected",C_result);

	if ( typeMeth == XL_TYPE_STRING )
	{
		if((C_meth = ARM_ConvForwardYieldMethod (C_strMeth, C_result)) == ARM_DEFAULT_ERR)
		{
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}
	}

	long retCode;

	if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 1)
	{
		retCode = ARMLOCAL_DiscountYield (LocalGetNumObjectId (C_curve),C_matu,(long)C_meth,
							  C_result);
	}
	else
		retCode = ARM_KO;

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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_DiscountYield" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


__declspec(dllexport) LPXLOPER WINAPI Local_ForwardPrice (LPXLOPER XL_curve,
														  LPXLOPER XL_matu1,
														  LPXLOPER XL_matu2)
{
	ADD_LOG("Local_ForwardPrice ");
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();

	// C variable
	CCString C_curve;
	double C_matu1;
	double C_matu2;
		
	// error
	static int error;
	static char* reason = "";

	XL_readStrCell(XL_curve,C_curve," ARM_ERR: curve id: object expected",C_result);
	XL_readNumCell(XL_matu1,C_matu1," ARM_ERR: maturity 1: numeric expected",C_result);
	XL_readNumCell(XL_matu2,C_matu2," ARM_ERR: maturity 2: numeric expected",C_result);

	long retCode;

	if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 1)
	{
		retCode = ARMLOCAL_ForwardPrice (LocalGetNumObjectId (C_curve),C_matu1,C_matu2,
							 C_result);
	}
	else
		retCode = ARM_KO;

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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_ForwardPrice" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


__declspec(dllexport) LPXLOPER WINAPI Local_ForwardYield (LPXLOPER XL_curve,
														  LPXLOPER XL_matu1,
														  LPXLOPER XL_matu2,
														  LPXLOPER XL_meth,
														  LPXLOPER XL_adjDaycount,
														  LPXLOPER XL_decompFreq,
														  LPXLOPER XL_daycount)
{
	ADD_LOG("Local_ForwardYield ");
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();

	// C variable
	CCString C_curve;
	double C_matu1;
	double C_matu2;

	CCString C_meth;
	long methId;

	CCString C_adjDaycount;
	long adjDaycountId;

	CCString C_decompFreq;
	long decompFreqId;

	CCString C_daycount;
	long daycountId;

	// error
	static int error;
	static char* reason = "";

	XL_readStrCell(XL_curve,C_curve," ARM_ERR: curve id: object expected",C_result);
	XL_readNumCell(XL_matu1,C_matu1," ARM_ERR: maturity 1: numeric expected",C_result);
	XL_readNumCell(XL_matu2,C_matu2," ARM_ERR: maturity 2: numeric expected",C_result);
	XL_readStrCellWD(XL_meth,C_meth,"CONT"," ARM_ERR: method: string expected",C_result);
	XL_readStrCellWD(XL_adjDaycount,C_adjDaycount,"NO"," ARM_ERR: adjustment: string expected",C_result);
	XL_readStrCellWD(XL_decompFreq,C_decompFreq,"-9999"," ARM_ERR: adjustment: string expected",C_result);
	XL_readStrCellWD(XL_daycount,C_daycount,"-1"," ARM_ERR: adjustment: string expected",C_result);

	if((methId = ARM_ConvForwardYieldMethod (C_meth, C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	if((adjDaycountId = ARM_ConvYesOrNo (C_adjDaycount, C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	if (strcmp(C_decompFreq,"-9999") == 0)
	{
		decompFreqId = -9999;
	}
	else
	{
		if((decompFreqId = ARM_ConvCompMeth (C_decompFreq, C_result)) == ARM_DEFAULT_ERR)
		{
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}
	}

	daycountId = ARM_ConvDayCount(C_daycount);

	long retCode;

	retCode = ARMLOCAL_ForwardYield (LocalGetNumObjectId (C_curve),
									 C_matu1,
									 C_matu2,
									 methId,
									 adjDaycountId,
									 decompFreqId,
									 daycountId,
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_ForwardYield" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}



__declspec(dllexport) LPXLOPER WINAPI Local_ZCLINT(LPXLOPER XL_matu,
												   LPXLOPER XL_rate,
												   LPXLOPER XL_meth,
												   LPXLOPER XL_aDate,
												   LPXLOPER XL_ccy,
												   LPXLOPER XL_interpMeth)
{
	ADD_LOG("Local_ZCLINT");
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	    ARM_NOCALCIFWIZ();

	    // C variable
	    VECTOR<CCString> C_matu;
	    VECTOR<double> matu;
	    VECTOR<double> C_rate;
	    CCString C_ccy;

        long isDouble;

	    CCString C_strMeth;
	    double C_meth;
	    long typeMeth;
	    
	    double C_aDate;

	    CCString C_interpMeth;
	    long interpMethId;

	    // error
	    static int error;
	    static char* reason = "";

	    XL_readStrVector(XL_matu,C_matu," ARM_ERR: maturities: array of numeric or string expected",DOUBLE_TYPE,C_result);
	    XL_readNumVector(XL_rate,C_rate," ARM_ERR: rates: array of numeric expected",C_result);
	    XL_readStrOrNumCell(XL_meth,C_strMeth,C_meth,typeMeth," ARM_ERR: decompounding method: string or numeric expected",C_result);
	    XL_readNumCell(XL_aDate,C_aDate," ARM_ERR: as of date: date expected",C_result);
	    XL_readStrCellWD(XL_ccy,C_ccy,"DEFAULT"," ARM_ERR: currency: string expected",C_result);
	    XL_readStrCellWD(XL_interpMeth,C_interpMeth,"LINEAR"," ARM_ERR: interpolation method: string expected",C_result);

	    if(C_ccy == "DEFAULT")
	    {
		    ARM_result currencyres;
		    ARMLOCAL_ARM_GetDefaultCurrency(currencyres);
		    if ( currencyres.getRetCode () != ARM_OK )
		    {
			    ARM_ARG_ERR();
			    return (LPXLOPER)&XL_result;
		    }
		    else
		    {
			    C_ccy = currencyres.getString ();
		    }
	    }

	    if ( typeMeth == XL_TYPE_STRING )
	    {
		    if((C_meth = ARM_ConvForwardYieldMethod(C_strMeth, C_result)) == ARM_DEFAULT_ERR)
		    {
			    ARM_ARG_ERR();
			    return (LPXLOPER)&XL_result;
		    }
	    }

	    if ((interpMethId = ARM_ConvInterpMethod(C_interpMeth, C_result)) == ARM_DEFAULT_ERR)
	    {
		    ARM_ARG_ERR();
		    return (LPXLOPER)&XL_result;
	    }

	    long retCode;
	    double dSettleDate;

	    retCode = ARMLOCAL_GetSpotDays(C_ccy,C_result);

	    if ( retCode == ARM_OK )
	    {
		    retCode = ARMLOCAL_NextBusinessDay(C_aDate,C_ccy,C_result.getDouble(),C_result);
		    dSettleDate = Local_ARMDATE2XLDATE(C_result.getString());
	    }
	    else
	    {
		    ARM_ERR();
		    return (LPXLOPER)&XL_result;
	    }

        ARM_CRV_TERMS sTerms;

	    for (int i = 0; i < C_matu.size(); i++)
	    {
		    isDouble = 0;

		    int Nb;
		    char cMatu;
		    long freqId;

		    sscanf(C_matu[i], "%d%c", &Nb, &cMatu);

            strcpy(sTerms[i], C_matu[i]);

		    cMatu = toupper(cMatu);

		    if ( cMatu == 'D' ) // Ex : "1D"
			    freqId = K_DAILY;
		    else if ( cMatu == 'W' )  
			    freqId = K_WEEKLY;
		    else if ( cMatu == 'M' ) 
			    freqId = K_MONTHLY;
		    else if ( cMatu == 'Y')  // ->implicitement ce sont des taux de swap
			    freqId = K_ANNUAL;
		    else
		    {
			    isDouble = 1;
			    matu.push_back (atof(C_matu[i]));
		    }
		    
		    if ( isDouble == 0 )
		    {
			    if (freqId == K_DAILY)
				    retCode = ARMLOCAL_NextBusinessDay(C_aDate,C_ccy,Nb,C_result);
			    else
				    retCode = ARMLOCAL_ARM_ADDPERIOD(dSettleDate, freqId, C_ccy, (long) Nb, K_FOLLOWING, 0, C_result);

			    if (retCode != ARM_OK)
			    {
				    ARM_ERR();
				    return (LPXLOPER)&XL_result;
			    }

			    double myDate = Local_ARMDATE2XLDATE(C_result.getString());

			    matu.push_back ((myDate - C_aDate) /365.);
		    }
	    }

        int           matuAreDoubles = isDouble; // 1: YES, O: NO
       


	    long objId;
	    CCString prevClass;

	    CCString curClass = LOCAL_ZERO_CURVE_LIN_CLASS;
	    CCString stringId = GetLastCurCellEnvValue ();

	    if (!stringId)
	    {
		    retCode = ARMLOCAL_zclint(matu, C_rate, (long) C_meth, C_aDate, C_ccy, interpMethId,
                                      matuAreDoubles,
                                      sTerms,
                                      C_result);

		    if ( retCode == ARM_OK )
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

		    if ( curClass == prevClass )
		    {
			    retCode = ARMLOCAL_zclint(matu, C_rate, (long)C_meth, C_aDate, C_ccy, interpMethId, 
                                          matuAreDoubles,
                                          sTerms,
                                          C_result, objId);

			    if ( retCode == ARM_OK )
			    {			
				    LocalSetCurCellEnvValue (curClass, objId); 

				    stringId = LocalMakeObjectId (objId, curClass);
			    }
		    }
		    else
		    {
			    FreeCurCellContent ();

			    retCode = ARMLOCAL_zclint(matu, C_rate, (long)C_meth, C_aDate, C_ccy, interpMethId,
                                          matuAreDoubles,
                                          sTerms,
                                          C_result);

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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_ZCLINT" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}



__declspec(dllexport) LPXLOPER WINAPI Local_PXL_ZCLINT(LPXLOPER XL_matu,
													   LPXLOPER XL_rate,
													   LPXLOPER XL_meth,
													   LPXLOPER XL_aDate,
													   LPXLOPER XL_ccy,
													   LPXLOPER XL_interpMeth)
{
	ADD_LOG("Local_PXL_ZCLINT");
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	    ARM_NOCALCIFWIZ();

	    // C variable
	    VECTOR<CCString> C_matu;
	    VECTOR<double> matu;
	    VECTOR<double> C_rate; 
	    CCString C_ccy;

        long isDouble;

	    CCString C_strMeth;
	    double C_meth;
	    long typeMeth;
	    
	    double C_aDate;

	    CCString C_interpMeth;
	    long interpMethId;

	    // error
	    static int error;
	    static char* reason = "";

	    XL_readStrVector(XL_matu,C_matu," ARM_ERR: maturities: array of string or numeric expected",DOUBLE_TYPE,C_result);
	    XL_readNumVector(XL_rate,C_rate," ARM_ERR: rates: array of numeric expected",C_result);
	    XL_readStrOrNumCell(XL_meth,C_strMeth,C_meth,typeMeth," ARM_ERR: decompounding method: string or numeric expected",C_result);
	    XL_readNumCell(XL_aDate,C_aDate," ARM_ERR: as of date: date expected",C_result);
	    XL_readStrCellWD(XL_ccy,C_ccy,"DEFAULT"," ARM_ERR: currency: string expected",C_result);
	    XL_readStrCellWD(XL_interpMeth,C_interpMeth,"LINEAR"," ARM_ERR: interpolation method: string expected",C_result);

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

	    if ( typeMeth == XL_TYPE_STRING )
	    {
		    if((C_meth = ARM_ConvForwardYieldMethod (C_strMeth, C_result)) == ARM_DEFAULT_ERR)
		    {
			    ARM_ARG_ERR();
			    return (LPXLOPER)&XL_result;
		    }
	    }

	    if((interpMethId = ARM_ConvInterpMethod (C_interpMeth, C_result)) == ARM_DEFAULT_ERR)
	    {
		    ARM_ARG_ERR();
		    return (LPXLOPER)&XL_result;
	    }

	    long retCode;
	    double dSettleDate;

	    retCode = ARMLOCAL_GetSpotDays(C_ccy,C_result);

	    if (retCode == ARM_OK)
	    {
		    retCode = ARMLOCAL_NextBusinessDay(C_aDate,C_ccy,C_result.getDouble(),C_result);
		    dSettleDate = Local_ARMDATE2XLDATE(C_result.getString());
	    }
	    else
	    {
		    ARM_ERR();
		    return (LPXLOPER)&XL_result;
	    }

        ARM_CRV_TERMS sTerms;

	    for (int i = 0; i < C_matu.size(); i++)
	    {
		    isDouble = 0;

		    int Nb;
		    char cMatu;
		    long freqId;

		    sscanf(C_matu[i], "%d%c", &Nb, &cMatu);

            strcpy(sTerms[i], C_matu[i]);

		    cMatu = toupper(cMatu);

		    if ( cMatu == 'D' ) // Ex : "1D"
			    freqId = K_DAILY;
		    else if ( cMatu == 'W' )  
			    freqId = K_WEEKLY;
		    else if ( cMatu == 'M' ) 
			    freqId = K_MONTHLY;
		    else if ( cMatu == 'Y')  // ->implicitement ce sont des taux de swap
			    freqId = K_ANNUAL;
		    else
		    {
			    isDouble = 1;
			    matu.push_back (atof(C_matu[i]));
		    }
		    
		    if ( isDouble == 0 )
		    {
			    if (freqId == K_DAILY)
				    retCode = ARMLOCAL_NextBusinessDay(C_aDate,C_ccy,Nb,C_result);
			    else
				    retCode = ARMLOCAL_ARM_ADDPERIOD(dSettleDate, freqId, C_ccy, (long) Nb, K_FOLLOWING, 0, C_result);

			    if (retCode != ARM_OK)
			    {
				    ARM_ERR();
				    return (LPXLOPER)&XL_result;
			    }

			    double myDate = Local_ARMDATE2XLDATE(C_result.getString());

			    matu.push_back ((myDate - C_aDate) /365.);
		    }
	    }

        int           matuAreDoubles = isDouble; // 1: YES, O: NO
        

	    long objId;
	    
	    CCString curClass = LOCAL_ZERO_CURVE_LIN_CLASS;
	    CCString stringId;
	    
	    retCode = ARMLOCAL_zclint(matu, C_rate, (long) C_meth, C_aDate, C_ccy, interpMethId, 
                                  matuAreDoubles,
                                  sTerms,
                                  C_result);

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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_PXL_ZCLINT" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}



__declspec(dllexport) LPXLOPER WINAPI Local_ZCSPREADED(LPXLOPER XL_zcSprId,
													   LPXLOPER XL_zcInitId,
													   LPXLOPER XL_date,
													   LPXLOPER XL_mmFreq,
													   LPXLOPER XL_swapFreq,
													   LPXLOPER XL_ccy)
{
	ADD_LOG("Local_ZCSPREADED");
	// return
	static XLOPER XL_result;
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	    ARM_NOCALCIFWIZ();

	    // C variable
	    CCString C_zcSprId;
	    CCString C_zcInitId;
	    double C_date;

	    CCString C_mmFreq;
	    long mmFreqId;

	    CCString C_swapFreq;
	    long swapFreqId;
	    
	    CCString C_ccy;
	    bool ccyIsObject = false;
	    
	    // error
	    static int error;
	    static char* reason = "";

	    XL_readStrCell(XL_zcSprId,C_zcSprId," ARM_ERR: zero-coupon curve spread id: object expected",C_result);
	    XL_readStrCell(XL_zcInitId,C_zcInitId," ARM_ERR: zero-coupon curve initialisation id: object expected",C_result);
	    XL_readNumCell(XL_date,C_date," ARM_ERR: as of date: date expected",C_result);
	    XL_readStrCell(XL_mmFreq,C_mmFreq," ARM_ERR: money market frequency: string expected",C_result);
	    XL_readStrCell(XL_swapFreq,C_swapFreq," ARM_ERR: swap frequency: string expected",C_result);
        XL_readStrCellWD(XL_ccy, C_ccy, "DEFAULT"," ARM_ERR: currency id: string expected",C_result);

	    long retCode;
	    long objId;
	    CCString prevClass;

		if ((C_ccy.GetLen() > 3)
			&&
			!(C_ccy == "DEFAULT"))
			ccyIsObject = true;

	    if ((mmFreqId = ARM_ConvFrequency(C_mmFreq, C_result)) == ARM_DEFAULT_ERR)
	    {
		    ARM_ARG_ERR();
		    return (LPXLOPER)&XL_result;
	    }

	    if ((swapFreqId = ARM_ConvFrequency(C_swapFreq, C_result)) == ARM_DEFAULT_ERR)
	    {
		    ARM_ARG_ERR();

		    return (LPXLOPER)&XL_result;
	    }
	    
	    CCString curClass = LOCAL_ZERO_CURVE_LIN_CLASS;
	    CCString stringId = GetLastCurCellEnvValue ();
	    
	    if (!stringId)
	    {
		    retCode = ARMLOCAL_zcspreaded(LocalGetNumObjectId(C_zcSprId), 
                                          LocalGetNumObjectId(C_zcInitId),
								          C_date, mmFreqId, swapFreqId, 
                                          ccyIsObject, C_ccy, C_result);

		    if ( retCode == ARM_OK )
		    {
			    objId = C_result.getLong ();

			    LocalSetCurCellEnvValue(curClass, objId); 

			    stringId = LocalMakeObjectId (objId, curClass);
		    }
	    }
	    else
	    {
		    prevClass = LocalGetStringObjectClass (stringId);
		    
		    objId = LocalGetNumObjectId (stringId);
			    
		    if ( curClass == prevClass )
		    {
			    retCode = ARMLOCAL_zcspreaded(LocalGetNumObjectId(C_zcSprId), 
                                              LocalGetNumObjectId(C_zcInitId),
								              C_date, mmFreqId, swapFreqId, 
                                              ccyIsObject, C_ccy, C_result, objId);

			    if ( retCode == ARM_OK )
			    {			
				    LocalSetCurCellEnvValue (curClass, objId); 

				    stringId = LocalMakeObjectId (objId, curClass);
			    }
		    }
		    else
		    {
			    FreeCurCellContent();

			    retCode = ARMLOCAL_zcspreaded(LocalGetNumObjectId(C_zcSprId), 
                                              LocalGetNumObjectId (C_zcInitId),
								              C_date, mmFreqId, swapFreqId, 
											  ccyIsObject, C_ccy, C_result);

			    if(retCode == ARM_OK)
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_ZCSPREADED" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}



__declspec(dllexport) LPXLOPER WINAPI Local_PXL_ZCSPREADED (LPXLOPER XL_zcSprId,
															LPXLOPER XL_zcInitId,
															LPXLOPER XL_date,
															LPXLOPER XL_mmFreq,
															LPXLOPER XL_swapFreq,
															LPXLOPER XL_ccy)
{
	ADD_LOG("Local_PXL_ZCSPREADED ");
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	    ARM_NOCALCIFWIZ();

	    // C variable
	    CCString C_zcSprId;
	    CCString C_zcInitId;
	    double C_date;

	    CCString C_mmFreq;
	    long mmFreqId;

	    CCString C_swapFreq;
	    long swapFreqId;
	    
	    CCString C_ccy;
	    bool ccyIsObject = false;
	    
	    // error
	    static int error;
	    static char* reason = "";

	    XL_readStrCell(XL_zcSprId,C_zcSprId," ARM_ERR: zero-coupon curve spread id: object expected",C_result);
	    XL_readStrCell(XL_zcInitId,C_zcInitId," ARM_ERR: zero-coupon curve initialisation id: object expected",C_result);
	    XL_readNumCell(XL_date,C_date," ARM_ERR: as of date: date expected",C_result);
	    XL_readStrCell(XL_mmFreq,C_mmFreq," ARM_ERR: money market frequency: string expected",C_result);
	    XL_readStrCell(XL_swapFreq,C_swapFreq," ARM_ERR: swap frequency: string expected",C_result);
	    XL_readStrCellWD(XL_ccy,C_ccy,"DEFAULT"," ARM_ERR: currency id: object expected",C_result);

		if ((C_ccy.GetLen() > 3)
			&&
			!(C_ccy == "DEFAULT"))
			ccyIsObject = true;

	    if((mmFreqId = ARM_ConvFrequency (C_mmFreq, C_result)) == ARM_DEFAULT_ERR)
	    {
		    ARM_ARG_ERR();
		    return (LPXLOPER)&XL_result;
	    }

	    if((swapFreqId = ARM_ConvFrequency (C_swapFreq, C_result)) == ARM_DEFAULT_ERR)
	    {
		    ARM_ARG_ERR();
		    return (LPXLOPER)&XL_result;
	    }
	    
	    long retCode;
	    long objId;
	    CCString curClass = LOCAL_ZERO_CURVE_LIN_CLASS;
	    CCString stringId;
	    
	    retCode = ARMLOCAL_zcspreaded(LocalGetNumObjectId(C_zcSprId),
                                      LocalGetNumObjectId(C_zcInitId),
							          C_date, mmFreqId, swapFreqId, 
                                      ccyIsObject, C_ccy, C_result);

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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_PXL_ZCSPREADED" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


__declspec(dllexport) LPXLOPER WINAPI Local_CreateZCSpreadedFromSummit(LPXLOPER XL_index,
																	   LPXLOPER XL_currency,
																	   LPXLOPER XL_cvName,
																	   LPXLOPER XL_aSdate,
																	   LPXLOPER XL_convAdj,
																	   LPXLOPER XL_raw,
																	   LPXLOPER XL_mmFreq,
																	   LPXLOPER XL_swapFreq,
																	   LPXLOPER XL_interp,
																	   LPXLOPER XL_zcinit)
{	
	ADD_LOG("Local_CreateZCSpreadedFromSummit");
	static XLOPER XL_result;
	ARM_result C_result;

	ARM_XL_TRY_BLOCK_BEGIN
	{
		ARM_NOCALCIFWIZ();

		// C variable
		CCString C_index;
		CCString C_currency;
		CCString C_cvName;
		double C_aSdate;
		
		CCString C_convAdj;
		long adjOrNotId;

		CCString C_raw;

		CCString C_mmFreq;
	    long mmFreq;

		CCString C_swapFreq;
		long swapFreq;

		CCString C_interp;
		long interpId;

		CCString C_zcInit;
		long zcInitId;

		// error
		static int error;
		static char* reason = "";

		XL_readStrCell(XL_index,C_index," ARM_ERR: index: string expected",C_result);
		XL_readStrCellWD(XL_currency,C_currency,"DEFAULT"," ARM_ERR: currency: string expected",C_result);
		XL_readStrCell(XL_cvName,C_cvName," ARM_ERR: curve name: string expected",C_result);
		XL_readNumCell(XL_aSdate,C_aSdate," ARM_ERR: as of date: date expected",C_result);
		XL_readStrCellWD(XL_convAdj,C_convAdj,"YES"," ARM_ERR: convexity adjustment: string expected",C_result);
		XL_readStrCellWD(XL_raw,C_raw,"DEFAULT"," ARM_ERR: FWD method: string expected",C_result);
		XL_readStrCellWD(XL_mmFreq,C_mmFreq,"-1"," ARM_ERR: money market frequency: string expected",C_result);
		XL_readStrCellWD(XL_swapFreq,C_swapFreq,"-1"," ARM_ERR: swap frequency: string expected",C_result);
		XL_readStrCellWD(XL_interp,C_interp,"L"," ARM_ERR: interpolation mode: string expected",C_result);
		XL_readStrCellWD(XL_zcinit,C_zcInit,"DEFAULT"," ARM_ERR: zc: object expected",C_result);

		C_index.toUpper ();
		C_currency.toUpper ();
		C_cvName.toUpper ();

		if((adjOrNotId = ARM_ConvYesOrNo (C_convAdj, C_result)) == ARM_DEFAULT_ERR)
		{
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}

		if ((mmFreq = ARM_ConvFrequency(C_mmFreq, C_result)) == ARM_DEFAULT_ERR)
	    {
		    ARM_ARG_ERR();
		    return (LPXLOPER)&XL_result;
	    }

		if ((swapFreq = ARM_ConvFrequency(C_swapFreq, C_result)) == ARM_DEFAULT_ERR)
		{
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}

		if((interpId = ARM_ConvInterpMethod (C_interp, C_result)) == ARM_DEFAULT_ERR)
		{
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}
	
		if(C_zcInit == "DEFAULT")
		{
			zcInitId = -1;
		}
		else
		{
			zcInitId = LocalGetNumObjectId(C_zcInit);
		}

	    long retCode;
	    long objId;
	    CCString prevClass;

	    
	    CCString curClass = LOCAL_ZERO_CURVE_LIN_CLASS;
	    CCString stringId = GetLastCurCellEnvValue ();

		if(!stringId)
		{
			retCode = ARMLOCAL_CreateZCSpreadedFromSummit(C_index, C_currency, C_cvName, C_aSdate, adjOrNotId, C_raw, 
														  swapFreq, mmFreq, interpId, zcInitId, C_result);

			if(retCode == ARM_OK)
				objId = C_result.getLong ();
		}
		else
		{
			prevClass = LocalGetStringObjectClass (stringId);
			objId = LocalGetNumObjectId (stringId);
				
			if(curClass == prevClass)
			{
				retCode = ARMLOCAL_CreateZCSpreadedFromSummit(C_index, C_currency, C_cvName, C_aSdate, adjOrNotId, C_raw, 
															  swapFreq, mmFreq, interpId, zcInitId, C_result, objId);
			}
			else
			{
				FreeCurCellContent ();
				retCode = ARMLOCAL_CreateZCSpreadedFromSummit(C_index, C_currency, C_cvName, C_aSdate, adjOrNotId, C_raw, 
															  swapFreq, mmFreq, interpId, zcInitId, C_result);
			
				if(retCode == ARM_OK)
					objId = C_result.getLong ();
			}
		}
		
		if(retCode == ARM_OK)
		{
			LocalSetCurCellEnvValue (curClass, objId); 
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

	}
	ARM_XL_TRY_BLOCK_END

	ARM_XL_CATCH_ARM_EXPT

	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_CreateZCSpreadedFromSummit" )

	return (LPXLOPER)&XL_result;
}


__declspec(dllexport) LPXLOPER WINAPI Local_PXL_CreateZCSpreadedFromSummit(LPXLOPER XL_index,
																		   LPXLOPER XL_currency,
																		   LPXLOPER XL_cvName,
																		   LPXLOPER XL_aSdate,
																		   LPXLOPER XL_convAdj,
																		   LPXLOPER XL_raw,
																		   LPXLOPER XL_mmFreq,
																		   LPXLOPER XL_swapFreq,
																		   LPXLOPER XL_interp,
																		   LPXLOPER XL_zcinit)
{	
	ADD_LOG("Local_PXL_CreateZCSpreadedFromSummit");
	static XLOPER XL_result;
	ARM_result C_result;

	ARM_XL_TRY_BLOCK_BEGIN
	{
		ARM_NOCALCIFWIZ();

		// C variable
		CCString C_index;
		CCString C_currency;
		CCString C_cvName;
		double C_aSdate;
		
		CCString C_convAdj;
		long adjOrNotId;

		CCString C_raw;

		CCString C_mmFreq;
	    long mmFreq;

		CCString C_swapFreq;
		long swapFreq;

		CCString C_interp;
		long interpId;

		CCString C_zcInit;
		long zcInitId;

		// error
		static int error;
		static char* reason = "";

		XL_readStrCell(XL_index,C_index," ARM_ERR: index: string expected",C_result);
		XL_readStrCellWD(XL_currency,C_currency,"DEFAULT"," ARM_ERR: currency: string expected",C_result);
		XL_readStrCell(XL_cvName,C_cvName," ARM_ERR: curve name: string expected",C_result);
		XL_readNumCell(XL_aSdate,C_aSdate," ARM_ERR: as of date: date expected",C_result);
		XL_readStrCellWD(XL_convAdj,C_convAdj,"YES"," ARM_ERR: convexity adjustment: string expected",C_result);
		XL_readStrCellWD(XL_raw,C_raw,"DEFAULT"," ARM_ERR: FWD method: string expected",C_result);
		XL_readStrCellWD(XL_mmFreq,C_mmFreq,"-1"," ARM_ERR: money market frequency: string expected",C_result);
		XL_readStrCellWD(XL_swapFreq,C_swapFreq,"-1"," ARM_ERR: swap frequency: string expected",C_result);
		XL_readStrCellWD(XL_interp,C_interp,"L"," ARM_ERR: interpolation mode: string expected",C_result);
		XL_readStrCellWD(XL_zcinit,C_zcInit,"DEFAULT"," ARM_ERR: zc: object expected",C_result);

		C_index.toUpper ();
		C_currency.toUpper ();
		C_cvName.toUpper ();

		if((adjOrNotId = ARM_ConvYesOrNo (C_convAdj, C_result)) == ARM_DEFAULT_ERR)
		{
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}

		if ((mmFreq = ARM_ConvFrequency(C_mmFreq, C_result)) == ARM_DEFAULT_ERR)
	    {
		    ARM_ARG_ERR();
		    return (LPXLOPER)&XL_result;
	    }

		if ((swapFreq = ARM_ConvFrequency(C_swapFreq, C_result)) == ARM_DEFAULT_ERR)
		{
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}

		if((interpId = ARM_ConvInterpMethod (C_interp, C_result)) == ARM_DEFAULT_ERR)
		{
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}
	
		if(C_zcInit == "DEFAULT")
		{
			zcInitId = -1;
		}
		else
		{
			zcInitId = LocalGetNumObjectId(C_zcInit);
		}

	    long retCode;
	    long objId;
	    CCString prevClass;
	    
	    CCString curClass = LOCAL_ZERO_CURVE_LIN_CLASS;
	    CCString stringId = GetLastCurCellEnvValue ();

		retCode = ARMLOCAL_CreateZCSpreadedFromSummit(C_index, C_currency, C_cvName, C_aSdate, adjOrNotId, C_raw, 
													  swapFreq, mmFreq, interpId, zcInitId, C_result);

		if(retCode == ARM_OK)
		{
			objId = C_result.getLong ();
			LocalSetCurCellEnvValue (curClass, objId); 
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

	}
	ARM_XL_TRY_BLOCK_END

	ARM_XL_CATCH_ARM_EXPT

	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_PXL_CreateZCSpreadedFromSummit" )

	return (LPXLOPER)&XL_result;
}


__declspec(dllexport) LPXLOPER WINAPI Local_CreateZCSwapIntSmooth (LPXLOPER XL_date,
																   LPXLOPER XL_matuRate,
																   LPXLOPER XL_mmVsFut,
																   LPXLOPER XL_swapVsFut,
																   LPXLOPER XL_raw,
																   LPXLOPER XL_interp,
																   LPXLOPER XL_ccy,
																   LPXLOPER XL_lambda,
																   LPXLOPER XL_prec)
{
	ADD_LOG("Local_CreateZCSwapIntSmooth ");
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();

	// C variable
	double C_date;
	VECTOR<CCString> C_matuRate;
	VECTOR<CCString> C_matu;
	VECTOR<double> C_rate; 

	CCString C_mmVsFut;
	long mmVsFutId;

	CCString C_swapVsFut;
	long swapVsFutId;

	CCString C_raw;
	long rawId;

	CCString C_interp;
	long interpId;

	CCString C_ccy;

	double C_lambda;
	double C_lambda_default = 0.0;

	double C_prec;
	double C_prec_default = 0;
	
	// error
	static int error;
	static char* reason = "";

	XL_readNumCell(XL_date,C_date," ARM_ERR: as of date: date expected",C_result);
	
	XL_readStrVector(XL_matuRate,C_matuRate," ARM_ERR: maturities and rates: array of numeric expected",DOUBLE_TYPE,C_result);

	double tmp_rate = -1;
	CCString tmp_matu ("NULL");

	if(C_matuRate.size () % 2 != 0)
	{
		C_result.setMsg ("ARM_ERR: check your maturities and rates array");
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	int real_size = C_matuRate.size () / 2;
	double tmp;
	int j = 0;

	for(int i = 0; i < real_size; i++)
	{
		C_matu.push_back (C_matuRate[j]);
		if((sscanf ((const char*)C_matuRate[j+1], "%lf", &tmp) != 1))
		{
			C_result.setMsg ("ARM_ERR: check your maturities and rates array");
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}
		C_rate.push_back (tmp);
		j += 2;
	}

	CCString last = C_matu[C_matu.size () - 1];
	if((!last) || (last[0] == ' '))
	{
		C_result.setMsg ("ARM_ERR: check your maturities and rates array");
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	XL_readStrCell(XL_mmVsFut,C_mmVsFut," ARM_ERR: money market versus future: string expected",C_result);
	XL_readStrCell(XL_swapVsFut,C_swapVsFut," ARM_ERR: swap versus future: string expected",C_result);
	XL_readStrCellWD(XL_raw,C_raw,"P"," ARM_ERR: raw: string expected",C_result);
	XL_readStrCellWD(XL_interp,C_interp,"L"," ARM_ERR: interpolation mode: string expected",C_result);
	XL_readStrCellWD(XL_ccy,C_ccy,"DEFAULT"," ARM_ERR: currency: string expected",C_result);
	XL_readNumCellWD(XL_lambda,C_lambda,C_lambda_default," ARM_ERR: lambda: numeric expected",C_result);
	XL_readNumCellWD(XL_prec,C_prec,C_prec_default," ARM_ERR: precision: numeric expected",C_result);

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

	if((mmVsFutId = ARM_ConvMktType (C_mmVsFut, C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	if((swapVsFutId = ARM_ConvMktType (C_swapVsFut, C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	rawId = ARM_ConvCvMethod (C_raw);

	if((interpId = ARM_ConvInterpMethod (C_interp, C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	long retCode;
	long objId;
	CCString prevClass;
	
	CCString curClass = LOCAL_ZERO_CURVE_INTSMO_CLASS;
	CCString stringId = GetLastCurCellEnvValue ();
	
	if(!stringId)
	{
		retCode = ARMLOCAL_CreateZCSwapIntSmooth (C_date, C_matu, C_rate, (long)mmVsFutId, 
											 (long)swapVsFutId, (long)rawId,
											 (long)interpId, C_ccy,
											 C_lambda,
											 (long)C_prec,
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
			retCode = ARMLOCAL_CreateZCSwapIntSmooth (C_date, C_matu, C_rate, (long)mmVsFutId, 
												 (long)swapVsFutId, (long)rawId,
												 (long)interpId, C_ccy, C_lambda, (long)C_prec, C_result, objId);

			if(retCode == ARM_OK)
			{			
				LocalSetCurCellEnvValue (curClass, objId); 

				stringId = LocalMakeObjectId (objId, curClass);
			}
		}
		else
		{
			FreeCurCellContent ();
			retCode = ARMLOCAL_CreateZCSwapIntSmooth (C_date, C_matu, C_rate, (long)mmVsFutId, 
												 (long)swapVsFutId, (long)rawId,
												 (long)interpId, C_ccy, C_lambda, (long)C_prec, C_result);
		
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

/*	if (cCcy)
		delete [] cCcy;
*/
//	ARM_END();
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_CreateZCSwapIntSmooth" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


__declspec(dllexport) LPXLOPER WINAPI Local_PXL_CreateZCSwapIntSmooth (LPXLOPER XL_date,
																	   LPXLOPER XL_matuRate,
																	   LPXLOPER XL_mmVsFut,
																	   LPXLOPER XL_swapVsFut,
																	   LPXLOPER XL_raw,
																	   LPXLOPER XL_interp,
																	   LPXLOPER XL_ccy,
																	   LPXLOPER XL_lambda,
																	   LPXLOPER XL_prec)
{
	ADD_LOG("Local_PXL_CreateZCSwapIntSmooth ");
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();

	// C variable
	double C_date;
	VECTOR<CCString> C_matuRate;
	VECTOR<CCString> C_matu;
	VECTOR<double> C_rate; 

	CCString C_mmVsFut;
	long mmVsFutId;

	CCString C_swapVsFut;
	long swapVsFutId;

	CCString C_raw;
	long rawId;

	CCString C_interp;
	long interpId;

	CCString C_ccy;

	double C_lambda;
	double C_lambda_default = 0.0;

	double C_prec;
	double C_prec_default = 0.0;
	
	// error
	static int error;
	static char* reason = "";

	XL_readNumCell(XL_date,C_date," ARM_ERR: as of date: date expected",C_result);
	
	XL_readStrVector(XL_matuRate,C_matuRate," ARM_ERR: maturities and rates: array of numeric expected",DOUBLE_TYPE,C_result);

	double tmp_rate = -1;
	CCString tmp_matu ("NULL");

	if(C_matuRate.size () % 2 != 0)
	{
		C_result.setMsg ("ARM_ERR: check your maturities and rates array");
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	int real_size = C_matuRate.size () / 2;
	double tmp;
	int j = 0;

	for(int i = 0; i < real_size; i++)
	{
		C_matu.push_back (C_matuRate[j]);
		if((sscanf ((const char*)C_matuRate[j+1], "%lf", &tmp) != 1))
		{
			C_result.setMsg ("ARM_ERR: check your maturities and rates array");
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}
		C_rate.push_back (tmp);
		j += 2;
	}

	CCString last = C_matu[C_matu.size () - 1];
	if((!last) || (last[0] == ' '))
	{
		C_result.setMsg ("ARM_ERR: check your maturities and rates array");
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	XL_readStrCell(XL_mmVsFut,C_mmVsFut," ARM_ERR: money market versus future: string expected",C_result);
	XL_readStrCell(XL_swapVsFut,C_swapVsFut," ARM_ERR: swap versus future: string expected",C_result);
	XL_readStrCellWD(XL_raw,C_raw,"P"," ARM_ERR: raw: string expected",C_result);
	XL_readStrCellWD(XL_interp,C_interp,"L"," ARM_ERR: interpolation mode: string expected",C_result);
	XL_readStrCellWD(XL_ccy,C_ccy,"DEFAULT"," ARM_ERR: currency: string expected",C_result);
	XL_readNumCellWD(XL_lambda,C_lambda,C_lambda_default," ARM_ERR: lambda: numeric expected",C_result);
	XL_readNumCellWD(XL_prec,C_prec,C_prec_default," ARM_ERR: precision: numeric expected",C_result);

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

	if((mmVsFutId = ARM_ConvMktType (C_mmVsFut, C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	if((swapVsFutId = ARM_ConvMktType (C_swapVsFut, C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	rawId = ARM_ConvCvMethod (C_raw);

	if((interpId = ARM_ConvInterpMethod (C_interp, C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	long retCode;
	long objId;
	
	CCString curClass = LOCAL_ZERO_CURVE_INTSMO_CLASS;
	CCString stringId;
	
	retCode = ARMLOCAL_CreateZCSwapIntSmooth (C_date, C_matu, C_rate, (long)mmVsFutId, 
										 (long)swapVsFutId, (long)rawId,
										 (long)interpId, C_ccy,
										 C_lambda, (long)C_prec,
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

/*	if (cCcy)
		delete [] cCcy;
*/
//	ARM_END();
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_PXL_CreateZCSwapIntSmooth" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


__declspec(dllexport) LPXLOPER WINAPI Local_CreateZCSwapFutInt (LPXLOPER XL_date,
																LPXLOPER XL_matuRate,
																LPXLOPER XL_mmVsFut,
																LPXLOPER XL_swapVsFut,
																LPXLOPER XL_raw,
																LPXLOPER XL_interp,
																LPXLOPER XL_ccy)
{
	ADD_LOG("Local_CreateZCSwapFutInt ");
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();

	// C variable
	double C_date;
	VECTOR<CCString> C_matuRate;
	VECTOR<CCString> C_matu;
	VECTOR<double> C_rate; 

	CCString C_mmVsFut;
	double mmVsFutId;

	CCString C_swapVsFut;
	double swapVsFutId;

	CCString C_raw;
	double rawId;

	CCString C_interp;
	double interpId;

	CCString C_ccy;
	
	// error
	static int error;
	static char* reason = "";

	XL_readNumCell(XL_date,C_date," ARM_ERR: as of date: date expected",C_result);
	
	XL_readStrVector(XL_matuRate,C_matuRate," ARM_ERR: maturities and rates: array of numeric expected",DOUBLE_TYPE,C_result);

	double tmp_rate = -1;
	CCString tmp_matu ("NULL");

	if(C_matuRate.size () % 2 != 0)
	{
		C_result.setMsg ("ARM_ERR: check your maturities and rates array");
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	int real_size = C_matuRate.size () / 2;
	double tmp;
	int j = 0;

	for(int i = 0; i < real_size; i++)
	{
		C_matu.push_back (C_matuRate[j]);
		if((sscanf ((const char*)C_matuRate[j+1], "%lf", &tmp) != 1))
		{
			C_result.setMsg ("ARM_ERR: check your maturities and rates array");
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}
		C_rate.push_back (tmp);
		j += 2;
	}

	CCString last = C_matu[C_matu.size () - 1];
	if((!last) || (last[0] == ' '))
	{
		C_result.setMsg ("ARM_ERR: check your maturities and rates array");
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	XL_readStrCell(XL_mmVsFut,C_mmVsFut," ARM_ERR: money market versus future: string expected",C_result);
	XL_readStrCell(XL_swapVsFut,C_swapVsFut," ARM_ERR: swap versus future: string expected",C_result);
	XL_readStrCellWD(XL_raw,C_raw,"P"," ARM_ERR: raw: string expected",C_result);
	XL_readStrCellWD(XL_interp,C_interp,"L"," ARM_ERR: interpolation mode: string expected",C_result);
	XL_readStrCellWD(XL_ccy,C_ccy,"DEFAULT"," ARM_ERR: currency: string expected",C_result);

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

	if((mmVsFutId = ARM_ConvMktType (C_mmVsFut, C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	if((swapVsFutId = ARM_ConvMktType (C_swapVsFut, C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	rawId = ARM_ConvCvMethod (C_raw);

	if((interpId = ARM_ConvInterpMethod (C_interp, C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	long retCode;
	long objId;
	CCString prevClass;
	
	CCString curClass = LOCAL_ZERO_CURVE_LIN_CLASS;
	CCString stringId = GetLastCurCellEnvValue ();
	
	if(!stringId)
	{
		retCode = ARMLOCAL_CreateZCSwapFutInt (C_date, C_matu, C_rate, (long)mmVsFutId, 
									      (long)swapVsFutId, (long)rawId,
										  (long)interpId, C_ccy, C_result);

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
			retCode = ARMLOCAL_CreateZCSwapFutInt (C_date, C_matu, C_rate, (long)mmVsFutId, 
											  (long)swapVsFutId, (long)rawId,
											  (long)interpId, C_ccy, C_result, objId);

			if(retCode == ARM_OK)
			{			
				LocalSetCurCellEnvValue (curClass, objId); 

				stringId = LocalMakeObjectId (objId, curClass);
			}
		}
		else
		{
			FreeCurCellContent ();
			retCode = ARMLOCAL_CreateZCSwapFutInt (C_date, C_matu, C_rate, (long)mmVsFutId, 
				                              (long)swapVsFutId, (long)rawId,
							                  (long)interpId, C_ccy, C_result);
		
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_CreateZCSwapFutInt" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}



__declspec(dllexport) LPXLOPER WINAPI Local_PXL_CreateZCSwapFutInt (LPXLOPER XL_date,
																	LPXLOPER XL_matuRate,
																	LPXLOPER XL_mmVsFut,
																	LPXLOPER XL_swapVsFut,
																	LPXLOPER XL_raw,
																	LPXLOPER XL_interp,
																	LPXLOPER XL_ccy)
{
	ADD_LOG("Local_PXL_CreateZCSwapFutInt ");
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();

	// C variable
	double C_date;
	VECTOR<CCString> C_matuRate;
	VECTOR<CCString> C_matu;
	VECTOR<double> C_rate; 

	CCString C_mmVsFut;
	double mmVsFutId;

	CCString C_swapVsFut;
	double swapVsFutId;

	CCString C_raw;
	double rawId;

	CCString C_interp;
	double interpId;

	CCString C_ccy;
	
	// error
	static int error;
	static char* reason = "";

	XL_readNumCell(XL_date,C_date," ARM_ERR: as of date: date expected",C_result);
	
	XL_readStrVector(XL_matuRate,C_matuRate," ARM_ERR: maturities and rates: array of numeric expected",DOUBLE_TYPE,C_result);

	double tmp_rate = -1;
	CCString tmp_matu ("NULL");

	if(C_matuRate.size () % 2 != 0)
	{
		C_result.setMsg ("ARM_ERR: check your maturities and rates array");
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	int real_size = C_matuRate.size () / 2;
	double tmp;
	int j = 0;

	for(int i = 0; i < real_size; i++)
	{
		C_matu.push_back (C_matuRate[j]);
		if((sscanf ((const char*)C_matuRate[j+1], "%lf", &tmp) != 1))
		{
			C_result.setMsg ("ARM_ERR: check your maturities and rates array");
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}
		C_rate.push_back (tmp);
		j += 2;
	}

	CCString last = C_matu[C_matu.size () - 1];
	if((!last) || (last[0] == ' '))
	{
		C_result.setMsg ("ARM_ERR: check your maturities and rates array");
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	XL_readStrCell(XL_mmVsFut,C_mmVsFut," ARM_ERR: money market versus future: string expected",C_result);
	XL_readStrCell(XL_swapVsFut,C_swapVsFut," ARM_ERR: swap versus future: string expected",C_result);
	XL_readStrCellWD(XL_raw,C_raw,"P"," ARM_ERR: raw: string expected",C_result);
	XL_readStrCellWD(XL_interp,C_interp,"L"," ARM_ERR: interpolation mode: string expected",C_result);
	XL_readStrCellWD(XL_ccy,C_ccy,"DEFAULT"," ARM_ERR: currency: string expected",C_result);

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

	if((mmVsFutId = ARM_ConvMktType (C_mmVsFut, C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	if((swapVsFutId = ARM_ConvMktType (C_swapVsFut, C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	rawId = ARM_ConvCvMethod (C_raw);

	if((interpId = ARM_ConvInterpMethod (C_interp, C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	long retCode;
	long objId;
	CCString prevClass;
	
	CCString curClass = LOCAL_ZERO_CURVE_LIN_CLASS;
	CCString stringId;
	
	retCode = ARMLOCAL_CreateZCSwapFutInt (C_date, C_matu, C_rate, (long)mmVsFutId, 
									      (long)swapVsFutId, (long)rawId,
										  (long)interpId, C_ccy, C_result);

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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_PXL_CreateZCSwapFutInt" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}



__declspec(dllexport) LPXLOPER WINAPI Local_CreateZCSwapFutIntSmooth (LPXLOPER XL_date,
																	  LPXLOPER XL_matuRate,
																	  LPXLOPER XL_mmVsFut,
																	  LPXLOPER XL_swapVsFut,
																	  LPXLOPER XL_raw,
																	  LPXLOPER XL_interp,
																	  LPXLOPER XL_ccy,
																	  LPXLOPER XL_lambda,
																	  LPXLOPER XL_prec)
{
	ADD_LOG("Local_CreateZCSwapFutIntSmooth ");
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();

	// C variable
	double C_date;
	VECTOR<CCString> C_matuRate;
	VECTOR<CCString> C_matu;
	VECTOR<double> C_rate; 

	CCString C_mmVsFut;
	long mmVsFutId;

	CCString C_swapVsFut;
	long swapVsFutId;

	CCString C_raw;
	long rawId;

	CCString C_interp;
	double interpId;

	CCString C_ccy;

	double C_lambda;
	double C_lambda_default = 0.0;

	double C_prec;
	double C_prec_default = 0.0;
	
	// error
	static int error;
	static char* reason = "";

	XL_readNumCell(XL_date,C_date," ARM_ERR: as of date: date expected",C_result);
	
	XL_readStrVector(XL_matuRate,C_matuRate," ARM_ERR: maturities and rates: array of numeric expected",DOUBLE_TYPE,C_result);

	double tmp_rate = -1;
	CCString tmp_matu ("NULL");

	if(C_matuRate.size () % 2 != 0)
	{
		C_result.setMsg ("ARM_ERR: check your maturities and rates array");
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	int real_size = C_matuRate.size () / 2;
	double tmp;
	int j = 0;

	for(int i = 0; i < real_size; i++)
	{
		C_matu.push_back (C_matuRate[j]);
		if((sscanf ((const char*)C_matuRate[j+1], "%lf", &tmp) != 1))
		{
			C_result.setMsg ("ARM_ERR: check your maturities and rates array");
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}
		C_rate.push_back (tmp);
		j += 2;
	}

	CCString last = C_matu[C_matu.size () - 1];
	if((!last) || (last[0] == ' '))
	{
		C_result.setMsg ("ARM_ERR: check your maturities and rates array");
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	XL_readStrCell(XL_mmVsFut,C_mmVsFut," ARM_ERR: money market versus future: string expected",C_result);
	XL_readStrCell(XL_swapVsFut,C_swapVsFut," ARM_ERR: swap versus future: string expected",C_result);
	XL_readStrCellWD(XL_raw,C_raw,"P"," ARM_ERR: raw: string expected",C_result);
	XL_readStrCellWD(XL_interp,C_interp,"L"," ARM_ERR: interpolation mode: string expected",C_result);
	XL_readStrCellWD(XL_ccy,C_ccy,"DEFAULT"," ARM_ERR: currency: string expected",C_result);
	XL_readNumCellWD(XL_lambda,C_lambda,C_lambda_default," ARM_ERR: lambda: numeric expected",C_result);
	XL_readNumCellWD(XL_prec,C_prec,C_prec_default," ARM_ERR: precision: numeric expected",C_result);

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

	if((mmVsFutId = ARM_ConvMktType (C_mmVsFut, C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	if((swapVsFutId = ARM_ConvMktType (C_swapVsFut, C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	rawId = ARM_ConvCvMethod (C_raw);

	if((interpId = ARM_ConvInterpMethod (C_interp, C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	long retCode;
	long objId;
	CCString prevClass;
	
	CCString curClass = LOCAL_ZERO_CURVE_INTSMO_CLASS;
	CCString stringId = GetLastCurCellEnvValue ();
	
	if(!stringId)
	{
		retCode = ARMLOCAL_CreateZCSwapFutIntSmooth (C_date, C_matu, C_rate, (long)mmVsFutId, 
												(long)swapVsFutId, (long)rawId,
											    (long)interpId, C_ccy,
												C_lambda, (long)C_prec, C_result);

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
			retCode = ARMLOCAL_CreateZCSwapFutIntSmooth (C_date, C_matu, C_rate, (long)mmVsFutId, 
													(long)swapVsFutId, (long)rawId,
													(long)interpId, C_ccy,
													C_lambda, (long)C_prec,
													C_result, objId);

			if(retCode == ARM_OK)
			{
				LocalSetCurCellEnvValue (curClass, objId); 

				stringId = LocalMakeObjectId (objId, curClass);
			}
		}
		else
		{
			FreeCurCellContent ();
			retCode = ARMLOCAL_CreateZCSwapFutIntSmooth (C_date, C_matu, C_rate, (long)mmVsFutId, 
													(long)swapVsFutId, (long)rawId,
													(long)interpId, C_ccy,
													C_lambda, (long)C_prec,
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_CreateZCSwapFutIntSmooth" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


__declspec(dllexport) LPXLOPER WINAPI Local_PXL_CreateZCSwapFutIntSmooth (LPXLOPER XL_date,
																		  LPXLOPER XL_matuRate,
																		  LPXLOPER XL_mmVsFut,
																		  LPXLOPER XL_swapVsFut,
																		  LPXLOPER XL_raw,
																		  LPXLOPER XL_interp,
																		  LPXLOPER XL_ccy,
																		  LPXLOPER XL_lambda,
																		  LPXLOPER XL_prec)
{
	ADD_LOG("Local_PXL_CreateZCSwapFutIntSmooth ");
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();

	// C variable
	double C_date;
	VECTOR<CCString> C_matuRate;
	VECTOR<CCString> C_matu;
	VECTOR<double> C_rate; 

	CCString C_mmVsFut;
	long mmVsFutId;

	CCString C_swapVsFut;
	long swapVsFutId;

	CCString C_raw;
	long rawId;

	CCString C_interp;
	double interpId;

	CCString C_ccy;

	double C_lambda;
	double C_lambda_default = 0.0;

	double C_prec;
	double C_prec_default = 0.0;
	
	// error
	static int error;
	static char* reason = "";

	XL_readNumCell(XL_date,C_date," ARM_ERR: as of date: date expected",C_result);
	
	XL_readStrVector(XL_matuRate,C_matuRate," ARM_ERR: maturities and rates: array of numeric expected",DOUBLE_TYPE,C_result);

	double tmp_rate = -1;
	CCString tmp_matu ("NULL");

	if(C_matuRate.size () % 2 != 0)
	{
		C_result.setMsg ("ARM_ERR: check your maturities and rates array");
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	int real_size = C_matuRate.size () / 2;
	double tmp;
	int j = 0;

	for(int i = 0; i < real_size; i++)
	{
		C_matu.push_back (C_matuRate[j]);
		if((sscanf ((const char*)C_matuRate[j+1], "%lf", &tmp) != 1))
		{
			C_result.setMsg ("ARM_ERR: check your maturities and rates array");
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}
		C_rate.push_back (tmp);
		j += 2;
	}

	CCString last = C_matu[C_matu.size () - 1];
	if((!last) || (last[0] == ' '))
	{
		C_result.setMsg ("ARM_ERR: check your maturities and rates array");
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	XL_readStrCell(XL_mmVsFut,C_mmVsFut," ARM_ERR: money market versus future: string expected",C_result);
	XL_readStrCell(XL_swapVsFut,C_swapVsFut," ARM_ERR: swap versus future: string expected",C_result);
	XL_readStrCellWD(XL_raw,C_raw,"P"," ARM_ERR: raw: string expected",C_result);
	XL_readStrCellWD(XL_interp,C_interp,"L"," ARM_ERR: interpolation mode: string expected",C_result);
	XL_readStrCellWD(XL_ccy,C_ccy,"DEFAULT"," ARM_ERR: currency: string expected",C_result);
	XL_readNumCellWD(XL_lambda,C_lambda,C_lambda_default," ARM_ERR: lambda: numeric expected",C_result);
	XL_readNumCellWD(XL_prec,C_prec,C_prec_default," ARM_ERR: precision: numeric expected",C_result);

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

	if((mmVsFutId = ARM_ConvMktType (C_mmVsFut, C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	if((swapVsFutId = ARM_ConvMktType (C_swapVsFut, C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	rawId = ARM_ConvCvMethod (C_raw);

	if((interpId = ARM_ConvInterpMethod (C_interp, C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	long retCode;
	long objId;
	CCString prevClass;
	
	CCString curClass = LOCAL_ZERO_CURVE_INTSMO_CLASS;
	CCString stringId;
	
	retCode = ARMLOCAL_CreateZCSwapFutIntSmooth (C_date, C_matu, C_rate, (long)mmVsFutId, 
												(long)swapVsFutId, (long)rawId,
											    (long)interpId, C_ccy,
												C_lambda, (long)C_prec, C_result);

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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_PXL_CreateZCSwapFutIntSmooth" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


__declspec(dllexport) LPXLOPER WINAPI Local_ShiftedZCSWAPINT (LPXLOPER XL_shiftValue,
															  LPXLOPER XL_nbPlot,
															  LPXLOPER XL_date,
															  LPXLOPER XL_matuRate,
															  LPXLOPER XL_mmVsFut,
															  LPXLOPER XL_swapVsFut,
															  LPXLOPER XL_raw,
															  LPXLOPER XL_interp,
															  LPXLOPER XL_ccy,
															  LPXLOPER XL_swapFrq,
															  LPXLOPER XL_fixDayCount)
{
	ADD_LOG("Local_ShiftedZCSWAPINT ");
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();

	// C variable
	double C_shiftValue;
	double C_nbPlot;

	double C_date;
	VECTOR<CCString> C_matuRate;
	VECTOR<CCString> C_matu;
	VECTOR<double> C_rate; 

	CCString C_mmVsFut;
	long mmVsFutId;

	CCString C_swapVsFut;
	long swapVsFutId;

	CCString C_raw;
	long rawId;

	CCString C_interp;
	long interpId;

	CCString C_ccy;

	CCString C_swapFrq;
	long swapFrqId;

	CCString C_fixDayCount;
	long fixDayCount;
	
	// error
	static int error;
	static char* reason = "";

	XL_readNumCell(XL_shiftValue,C_shiftValue," ARM_ERR: shift value: numeric expected",C_result);
	XL_readNumCell(XL_nbPlot,C_nbPlot," ARM_ERR: number of plots: numeric expected",C_result);
	XL_readNumCell(XL_date,C_date," ARM_ERR: as of date: date expected",C_result);	
	XL_readStrVector(XL_matuRate,C_matuRate," ARM_ERR: maturities and rates: array expected",DOUBLE_TYPE,C_result);

	double tmp_rate = -1;
	CCString tmp_matu ("NULL");

	if(C_matuRate.size () % 2 != 0)
	{
		C_result.setMsg ("ARM_ERR: check your maturities and rates array");
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	int real_size = C_matuRate.size () / 2;
	double tmp;
	int j = 0;

	for(int i = 0; i < real_size; i++)
	{
		C_matu.push_back (C_matuRate[j]);
		if((sscanf ((const char*)C_matuRate[j+1], "%lf", &tmp) != 1))
		{
			C_result.setMsg ("ARM_ERR: check your maturities and rates array");
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}
		C_rate.push_back (tmp);
		j += 2;
	}

	CCString last = C_matu[C_matu.size () - 1];
	if((!last) || (last[0] == ' '))
	{
		C_result.setMsg ("ARM_ERR: check your maturities and rates array");
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	if((long)C_nbPlot < 0)
		C_nbPlot = 0;
	else if((long)C_nbPlot > C_rate.size ())
		C_nbPlot = C_rate.size ();
	
	for(i = 0; i < (long)C_nbPlot; i++)
	{
		if (C_rate[i] > 50.0)
			C_rate[i] -= C_shiftValue;
		else
			C_rate[i] += C_shiftValue;
	}

	XL_readStrCell(XL_mmVsFut,C_mmVsFut," ARM_ERR: money market versus: string expected",C_result);
	XL_readStrCell(XL_swapVsFut,C_swapVsFut," ARM_ERR: swap versus future: string expected",C_result);
	XL_readStrCellWD(XL_raw,C_raw,"P"," ARM_ERR: raw: string expected",C_result);
	XL_readStrCellWD(XL_interp,C_interp,"L"," ARM_ERR: interpolation mode: string expected",C_result);
	XL_readStrCellWD(XL_ccy,C_ccy,"DEFAULT"," ARM_ERR: currency: string expected",C_result);
	XL_readStrCellWD(XL_swapFrq,C_swapFrq,"-1"," ARM_ERR: swap frequency: string expected",C_result);
	XL_readStrCellWD(XL_fixDayCount,C_fixDayCount,"-1"," ARM_ERR: swap Fixed Leg DayCount: string expected",C_result);

	if(C_ccy == "DEFAULT")
	{
		ARM_result ccyres;
		ARMLOCAL_ARM_GetDefaultCurrency (ccyres);
		if(ccyres.getRetCode () != ARM_OK)
		{
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}
		else
		{
			C_ccy = ccyres.getString ();
		}
	}
	
	if((mmVsFutId = ARM_ConvMktType (C_mmVsFut, C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	if((swapVsFutId = ARM_ConvMktType (C_swapVsFut, C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	rawId = ARM_ConvCvMethod (C_raw);

	if((interpId = ARM_ConvInterpMethod (C_interp, C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	if ((swapFrqId = ARM_ConvFrequency(C_swapFrq, C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	if ((fixDayCount = ARM_ConvDayCount(C_fixDayCount)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	long retCode;
	long objId;
	CCString prevClass;
	
	CCString curClass = LOCAL_ZERO_CURVE_LIN_CLASS;
	CCString stringId = GetLastCurCellEnvValue ();
	
	if(!stringId)
	{
		retCode = ARMLOCAL_CreateZCSwapInt (C_date, C_matu, C_rate, (long)mmVsFutId, 
									   (long)swapVsFutId, (long)rawId,
									   (long)interpId, C_ccy, swapFrqId, fixDayCount, C_result);

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
			retCode = ARMLOCAL_CreateZCSwapInt (C_date, C_matu, C_rate, (long)mmVsFutId, 
										   (long)swapVsFutId, (long)rawId,
										   (long)interpId, C_ccy, swapFrqId, fixDayCount, C_result, objId);

			if(retCode == ARM_OK)
			{
				LocalSetCurCellEnvValue (curClass, objId); 

				stringId = LocalMakeObjectId (objId, curClass);
			}
		}
		else
		{
			FreeCurCellContent ();
			retCode = ARMLOCAL_CreateZCSwapInt (C_date, C_matu, C_rate, (long)mmVsFutId, 
				                           (long)swapVsFutId, (long)rawId,
							               (long)interpId, C_ccy, swapFrqId, fixDayCount, C_result);

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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_ShiftedZCSWAPINT" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


__declspec(dllexport) LPXLOPER WINAPI Local_PXL_ShiftedZCSWAPINT (LPXLOPER XL_shiftValue,
																  LPXLOPER XL_nbPlot,
																  LPXLOPER XL_date,
																  LPXLOPER XL_matuRate,
																  LPXLOPER XL_mmVsFut,
																  LPXLOPER XL_swapVsFut,
																  LPXLOPER XL_raw,
																  LPXLOPER XL_interp,
																  LPXLOPER XL_ccy,
															      LPXLOPER XL_swapFrq,
															      LPXLOPER XL_fixDayCount)
{
	ADD_LOG("Local_PXL_ShiftedZCSWAPINT ");
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();

	// C variable
	double C_shiftValue;
	double C_nbPlot;

	double C_date;
	VECTOR<CCString> C_matuRate;
	VECTOR<CCString> C_matu;
	VECTOR<double> C_rate; 

	CCString C_mmVsFut;
	long mmVsFutId;

	CCString C_swapVsFut;
	long swapVsFutId;

	CCString C_raw;
	long rawId;

	CCString C_interp;
	long interpId;

	CCString C_ccy;
	
	CCString C_swapFrq;
	long swapFrqId;

	CCString C_fixDayCount;
	long fixDayCount;

	// error
	static int error;
	static char* reason = "";

	XL_readNumCell(XL_shiftValue,C_shiftValue," ARM_ERR: shift value: numeric expected",C_result);
	XL_readNumCell(XL_nbPlot,C_nbPlot," ARM_ERR: number of plots: numeric expected",C_result);
	XL_readNumCell(XL_date,C_date," ARM_ERR: as of date: date expected",C_result);	
	XL_readStrVector(XL_matuRate,C_matuRate," ARM_ERR: maturities and rates: array expected",DOUBLE_TYPE,C_result);

	double tmp_rate = -1;
	CCString tmp_matu ("NULL");

	if(C_matuRate.size () % 2 != 0)
	{
		C_result.setMsg ("ARM_ERR: check your maturities and rates array");
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	int real_size = C_matuRate.size () / 2;
	double tmp;
	int j = 0;

	for(int i = 0; i < real_size; i++)
	{
		C_matu.push_back (C_matuRate[j]);
		if((sscanf ((const char*)C_matuRate[j+1], "%lf", &tmp) != 1))
		{
			C_result.setMsg ("ARM_ERR: check your maturities and rates array");
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}
		C_rate.push_back (tmp);
		j += 2;
	}

	CCString last = C_matu[C_matu.size () - 1];
	if((!last) || (last[0] == ' '))
	{
		C_result.setMsg ("ARM_ERR: check your maturities and rates array");
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	if((long)C_nbPlot < 0)
		C_nbPlot = 0;
	else if((long)C_nbPlot > C_rate.size ())
		C_nbPlot = C_rate.size ();
	
	for(i = 0; i < (long)C_nbPlot; i++)
	{
		if (C_rate[i] > 50.0)
			C_rate[i] -= C_shiftValue;
		else
			C_rate[i] += C_shiftValue;
	}

	XL_readStrCell(XL_mmVsFut,C_mmVsFut," ARM_ERR: money market versus: string expected",C_result);
	XL_readStrCell(XL_swapVsFut,C_swapVsFut," ARM_ERR: swap versus future: string expected",C_result);
	XL_readStrCellWD(XL_raw,C_raw,"P"," ARM_ERR: raw: string expected",C_result);
	XL_readStrCellWD(XL_interp,C_interp,"L"," ARM_ERR: interpolation mode: string expected",C_result);
	XL_readStrCellWD(XL_ccy,C_ccy,"DEFAULT"," ARM_ERR: currency: string expected",C_result);
	XL_readStrCellWD(XL_swapFrq,C_swapFrq,"-1"," ARM_ERR: swap frequency: string expected",C_result);
	XL_readStrCellWD(XL_fixDayCount,C_fixDayCount,"-1"," ARM_ERR: swap Fixed Leg DayCount: string expected",C_result);


	if(C_ccy == "DEFAULT")
	{
		ARM_result ccyres;
		ARMLOCAL_ARM_GetDefaultCurrency (ccyres);
		if(ccyres.getRetCode () != ARM_OK)
		{
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}
		else
		{
			C_ccy = ccyres.getString ();
		}
	}
	
	if((mmVsFutId = ARM_ConvMktType (C_mmVsFut, C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	if((swapVsFutId = ARM_ConvMktType (C_swapVsFut, C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	rawId = ARM_ConvCvMethod (C_raw);

	if((interpId = ARM_ConvInterpMethod (C_interp, C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	if ((swapFrqId = ARM_ConvFrequency(C_swapFrq, C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	if ((fixDayCount = ARM_ConvDayCount(C_fixDayCount)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	long retCode;
	long objId;
	CCString prevClass;
	
	CCString curClass = LOCAL_ZERO_CURVE_LIN_CLASS;
	CCString stringId;
	
	retCode = ARMLOCAL_CreateZCSwapInt (C_date, C_matu, C_rate, (long)mmVsFutId, 
									   (long)swapVsFutId, (long)rawId,
									   (long)interpId, C_ccy, swapFrqId, fixDayCount, C_result);

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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_PXL_ShiftedZCSWAPINT" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}



__declspec(dllexport) LPXLOPER WINAPI Local_ShiftedZCSWAPINTSmooth (LPXLOPER XL_shiftValue,
																	LPXLOPER XL_nbPlot,
																	LPXLOPER XL_date,
																	LPXLOPER XL_matuRate,
																	LPXLOPER XL_mmVsFut,
																	LPXLOPER XL_swapVsFut,
																	LPXLOPER XL_raw,
																	LPXLOPER XL_interp,
																	LPXLOPER XL_lambda,
																	LPXLOPER XL_prec,
																	LPXLOPER XL_ccy)
{
	ADD_LOG("Local_ShiftedZCSWAPINTSmooth ");
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();

	// C variable
	double C_shiftValue;
	double C_nbPlot;

	double C_date;
	VECTOR<CCString> C_matuRate;
	VECTOR<CCString> C_matu;
	VECTOR<double> C_rate; 

	CCString C_mmVsFut;
	long mmVsFutId;

	CCString C_swapVsFut;
	long swapVsFutId;

	CCString C_raw;
	long rawId;

	CCString C_interp;
	long interpId;

	CCString C_ccy;

	double C_lambda;
	double C_lambda_default = 0.0;

	double C_prec;
	double C_prec_default = 0.0;
	
	// error
	static int error;
	static char* reason = "";

	XL_readNumCell(XL_shiftValue,C_shiftValue," ARM_ERR: shift value: numeric expected",C_result);
	XL_readNumCell(XL_nbPlot,C_nbPlot," ARM_ERR: number of plots: numeric expected",C_result);
	XL_readNumCell(XL_date,C_date," ARM_ERR: as of date: date expected",C_result);
	
	XL_readStrVector(XL_matuRate,C_matuRate," ARM_ERR: maturities and rates: array expected",DOUBLE_TYPE,C_result);

	double tmp_rate = -1;
	CCString tmp_matu ("NULL");

	if(C_matuRate.size () % 2 != 0)
	{
		C_result.setMsg ("ARM_ERR: check your maturities and rates array");
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	int real_size = C_matuRate.size () / 2;
	double tmp;
	int j = 0;

	for(int i = 0; i < real_size; i++)
	{
		C_matu.push_back (C_matuRate[j]);
		if((sscanf ((const char*)C_matuRate[j+1], "%lf", &tmp) != 1))
		{
			C_result.setMsg ("ARM_ERR: check your maturities and rates array");
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}
		C_rate.push_back (tmp);
		j += 2;
	}

	CCString last = C_matu[C_matu.size () - 1];
	if((!last) || (last[0] == ' '))
	{
		C_result.setMsg ("ARM_ERR: check your maturities and rates array");
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	if((long)C_nbPlot < 0)
		C_nbPlot = 0;
	else if((long)C_nbPlot > C_rate.size ())
		C_nbPlot = C_rate.size ();
	
	for(i = 0; i < (long)C_nbPlot; i++)
	{
		if (C_rate[i] > 50.0)
			C_rate[i] -= C_shiftValue;
		else
			C_rate[i] += C_shiftValue;
	}

	XL_readStrCell(XL_mmVsFut,C_mmVsFut," ARM_ERR: money market versus: string expected",C_result);
	XL_readStrCell(XL_swapVsFut,C_swapVsFut," ARM_ERR: swap versus future: string expected",C_result);
	XL_readStrCellWD(XL_raw,C_raw,"P"," ARM_ERR: raw: string expected",C_result);
	XL_readStrCellWD(XL_interp,C_interp,"L"," ARM_ERR: interpolation mode: string expected",C_result);
	XL_readStrCellWD(XL_ccy,C_ccy,"DEFAULT"," ARM_ERR: currency: string expected",C_result);
	XL_readNumCellWD(XL_lambda,C_lambda,C_lambda_default," ARM_ERR: lambda: numeric expected",C_result);
	XL_readNumCellWD(XL_prec,C_prec,C_prec_default," ARM_ERR: precision: numeric expected",C_result);

	if(C_ccy == "DEFAULT")
	{
		ARM_result ccyres;
		ARMLOCAL_ARM_GetDefaultCurrency (ccyres);
		if(ccyres.getRetCode () != ARM_OK)
		{
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}
		else
		{
			C_ccy = ccyres.getString ();
		}
	}
	
	if((mmVsFutId = ARM_ConvMktType (C_mmVsFut, C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	if((swapVsFutId = ARM_ConvMktType (C_swapVsFut, C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	rawId = ARM_ConvCvMethod (C_raw);

	if((interpId = ARM_ConvInterpMethod (C_interp, C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	long retCode;
	long objId;
	CCString prevClass;
	
	CCString curClass = LOCAL_ZERO_CURVE_INTSMO_CLASS;
	CCString stringId = GetLastCurCellEnvValue ();
	
	if(!stringId)
	{
		retCode = ARMLOCAL_CreateZCSwapIntSmooth (C_date, C_matu, C_rate, (long)mmVsFutId, 
											 (long)swapVsFutId, (long)rawId,
											 (long)interpId, C_ccy,
											 C_lambda, (long)C_prec, C_result);

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
			retCode = ARMLOCAL_CreateZCSwapIntSmooth (C_date, C_matu, C_rate, (long)mmVsFutId, 
												 (long)swapVsFutId, (long)rawId,
												 (long)interpId, C_ccy,
												 C_lambda, (long)C_prec, C_result, objId);

			if(retCode == ARM_OK)
			{
				LocalSetCurCellEnvValue (curClass, objId); 

				stringId = LocalMakeObjectId (objId, curClass);
			}
		}
		else
		{
			FreeCurCellContent ();
			retCode = ARMLOCAL_CreateZCSwapIntSmooth (C_date, C_matu, C_rate, (long)mmVsFutId, 
												 (long)swapVsFutId, (long)rawId,
												 (long)interpId, C_ccy,
												 C_lambda, (long)C_prec, C_result);
		
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_ShiftedZCSWAPINTSmooth" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}



__declspec(dllexport) LPXLOPER WINAPI Local_PXL_ShiftedZCSWAPINTSmooth (LPXLOPER XL_shiftValue,
																		LPXLOPER XL_nbPlot,
																		LPXLOPER XL_date,
																		LPXLOPER XL_matuRate,
																		LPXLOPER XL_mmVsFut,
																		LPXLOPER XL_swapVsFut,
																		LPXLOPER XL_raw,
																		LPXLOPER XL_interp,
																		LPXLOPER XL_lambda,
																		LPXLOPER XL_prec,
																		LPXLOPER XL_ccy)
{
	ADD_LOG("Local_PXL_ShiftedZCSWAPINTSmooth ");
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();

	// C variable
	double C_shiftValue;
	double C_nbPlot;

	double C_date;
	VECTOR<CCString> C_matuRate;
	VECTOR<CCString> C_matu;
	VECTOR<double> C_rate; 

	CCString C_mmVsFut;
	long mmVsFutId;

	CCString C_swapVsFut;
	long swapVsFutId;

	CCString C_raw;
	long rawId;

	CCString C_interp;
	long interpId;

	CCString C_ccy;

	double C_lambda;
	double C_lambda_default = 0.0;

	double C_prec;
	double C_prec_default = 0.0;
	
	// error
	static int error;
	static char* reason = "";

	XL_readNumCell(XL_shiftValue,C_shiftValue," ARM_ERR: shift value: numeric expected",C_result);
	XL_readNumCell(XL_nbPlot,C_nbPlot," ARM_ERR: number of plots: numeric expected",C_result);
	XL_readNumCell(XL_date,C_date," ARM_ERR: as of date: date expected",C_result);
	
	XL_readStrVector(XL_matuRate,C_matuRate," ARM_ERR: maturities and rates: array expected",DOUBLE_TYPE,C_result);

	double tmp_rate = -1;
	CCString tmp_matu ("NULL");

	if(C_matuRate.size () % 2 != 0)
	{
		C_result.setMsg ("ARM_ERR: check your maturities and rates array");
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	int real_size = C_matuRate.size () / 2;
	double tmp;
	int j = 0;

	for(int i = 0; i < real_size; i++)
	{
		C_matu.push_back (C_matuRate[j]);
		if((sscanf ((const char*)C_matuRate[j+1], "%lf", &tmp) != 1))
		{
			C_result.setMsg ("ARM_ERR: check your maturities and rates array");
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}
		C_rate.push_back (tmp);
		j += 2;
	}

	CCString last = C_matu[C_matu.size () - 1];
	if((!last) || (last[0] == ' '))
	{
		C_result.setMsg ("ARM_ERR: check your maturities and rates array");
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	if((long)C_nbPlot < 0)
		C_nbPlot = 0;
	else if((long)C_nbPlot > C_rate.size ())
		C_nbPlot = C_rate.size ();
	
	for(i = 0; i < (long)C_nbPlot; i++)
	{
		if (C_rate[i] > 50.0)
			C_rate[i] -= C_shiftValue;
		else
			C_rate[i] += C_shiftValue;
	}

	XL_readStrCell(XL_mmVsFut,C_mmVsFut," ARM_ERR: money market versus: string expected",C_result);
	XL_readStrCell(XL_swapVsFut,C_swapVsFut," ARM_ERR: swap versus future: string expected",C_result);
	XL_readStrCellWD(XL_raw,C_raw,"P"," ARM_ERR: raw: string expected",C_result);
	XL_readStrCellWD(XL_interp,C_interp,"L"," ARM_ERR: interpolation mode: string expected",C_result);
	XL_readStrCellWD(XL_ccy,C_ccy,"DEFAULT"," ARM_ERR: currency: string expected",C_result);
	XL_readNumCellWD(XL_lambda,C_lambda,C_lambda_default," ARM_ERR: lambda: numeric expected",C_result);
	XL_readNumCellWD(XL_prec,C_prec,C_prec_default," ARM_ERR: precision: numeric expected",C_result);

	if(C_ccy == "DEFAULT")
	{
		ARM_result ccyres;
		ARMLOCAL_ARM_GetDefaultCurrency (ccyres);
		if(ccyres.getRetCode () != ARM_OK)
		{
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}
		else
		{
			C_ccy = ccyres.getString ();
		}
	}
	
	if((mmVsFutId = ARM_ConvMktType (C_mmVsFut, C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	if((swapVsFutId = ARM_ConvMktType (C_swapVsFut, C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	rawId = ARM_ConvCvMethod (C_raw);

	if((interpId = ARM_ConvInterpMethod (C_interp, C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	long retCode;
	long objId;
	CCString prevClass;
	
	CCString curClass = LOCAL_ZERO_CURVE_INTSMO_CLASS;
	CCString stringId;
	
	retCode = ARMLOCAL_CreateZCSwapIntSmooth (C_date, C_matu, C_rate, (long)mmVsFutId, 
											 (long)swapVsFutId, (long)rawId,
											 (long)interpId, C_ccy,
											 C_lambda, (long)C_prec, C_result);

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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_PXL_ShiftedZCSWAPINTSmooth" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


__declspec(dllexport) LPXLOPER WINAPI Local_TRANS2SMOOTH (LPXLOPER XL_inCv,
														  LPXLOPER XL_lambda,
														  LPXLOPER XL_prec)
{
	ADD_LOG("Local_TRANS2SMOOTH ");
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();

	// C variable
	CCString C_cvString;
	long C_inCvId;

	double C_lambda;
	double C_lambda_default = 0.0;

	double C_prec;
	double C_prec_default = 0.0;

	// error
	static int error;
	static char* reason = "";

	XL_readStrCell(XL_inCv,C_cvString," ARM_ERR: curve id: object expected",C_result);
	XL_readNumCellWD(XL_lambda,C_lambda,C_lambda_default," ARM_ERR: lambda: numeric expected",C_result);
	XL_readNumCellWD(XL_prec,C_prec,C_prec_default," ARM_ERR: precision: numeric expected",C_result);

	C_inCvId =  LocalGetNumObjectId (C_cvString);

	long retCode;
	long objId;
	CCString prevClass;

	CCString curClass = LOCAL_ZERO_CURVE_INTSMO_CLASS;
	CCString stringId = GetLastCurCellEnvValue ();

	if(!stringId)
	{
		retCode = ARMLOCAL_ZCINTSMOOTH (C_inCvId, C_lambda, (long)C_prec, C_result);

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
			retCode = ARMLOCAL_ZCINTSMOOTH (C_inCvId, C_lambda, (long)C_prec, C_result, objId);

			if(retCode == ARM_OK)
			{
				LocalSetCurCellEnvValue (curClass, objId); 

				stringId = LocalMakeObjectId (objId, curClass);
			}
		}
		else
		{
			FreeCurCellContent ();
			retCode = ARMLOCAL_ZCINTSMOOTH (C_inCvId, C_lambda, (long)C_prec, C_result);
		
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_TRANS2SMOOTH" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}

__declspec(dllexport) LPXLOPER WINAPI Local_PXL_TRANS2SMOOTH (LPXLOPER XL_inCv,
															  LPXLOPER XL_lambda,
															  LPXLOPER XL_prec)
{
	ADD_LOG("Local_PXL_TRANS2SMOOTH ");
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();

	// C variable
    CCString C_cvString;     
	long C_inCvId;

	double C_lambda;
	double C_lambda_default = 0.0;

	double C_prec;
	double C_prec_default = 0.0;

	// error
	static int error;
	static char* reason = "";

    XL_readStrCell(XL_inCv,C_cvString," ARM_ERR: curve id: object expected",C_result);
	XL_readNumCellWD(XL_lambda,C_lambda,C_lambda_default," ARM_ERR: lambda: numeric expected",C_result);
	XL_readNumCellWD(XL_prec,C_prec,C_prec_default," ARM_ERR: precision: numeric expected",C_result);

	C_inCvId =  LocalGetNumObjectId (C_cvString);

	long retCode;
	long objId;
	
	CCString curClass = LOCAL_ZERO_CURVE_INTSMO_CLASS;
	CCString stringId;
	
	retCode = ARMLOCAL_ZCINTSMOOTH (C_inCvId, C_lambda, (long)C_prec, C_result);

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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_PXL_TRANS2SMOOTH" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}



__declspec(dllexport) LPXLOPER WINAPI Local_ZCSWAPCUBDIFF (LPXLOPER XL_date,
														   LPXLOPER XL_matuRate,
														   LPXLOPER XL_mmVsFut,
														   LPXLOPER XL_swapVsFut,
														   LPXLOPER XL_raw,
														   LPXLOPER XL_ccy)
{
	ADD_LOG("Local_ZCSWAPCUBDIFF ");
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();

	// C variable
	double C_date;
	VECTOR<CCString> C_matuRate;
	VECTOR<CCString> C_matu;
	VECTOR<double> C_rate; 

	CCString C_mmVsFut;
	long mmVsFutId;

	CCString C_swapVsFut;
	long swapVsFutId;

	CCString C_raw;
	long rawId;

	CCString C_ccy;

	long interpId;
	
	if((interpId = ARM_ConvInterpMethod ("L", C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}
	
	// error
	static int error;
	static char* reason = "";

	XL_readNumCell(XL_date,C_date," ARM_ERR: as of date: date expected",C_result);
	XL_readStrVector(XL_matuRate,C_matuRate," ARM_ERR: maturities and rates: array expected",DOUBLE_TYPE,C_result);

	double tmp_rate = -1;
	CCString tmp_matu ("NULL");

	if(C_matuRate.size () % 2 != 0)
	{
		C_result.setMsg ("ARM_ERR: check your maturities and rates array");
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	int real_size = C_matuRate.size () / 2;
	double tmp;
	int j = 0;

	for(int i = 0; i < real_size; i++)
	{
		C_matu.push_back (C_matuRate[j]);
		if((sscanf ((const char*)C_matuRate[j+1], "%lf", &tmp) != 1))
		{
			C_result.setMsg ("ARM_ERR: check your maturities and rates array");
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}
		C_rate.push_back (tmp);
		j += 2;
	}

	CCString last = C_matu[C_matu.size () - 1];
	if((!last) || (last[0] == ' '))
	{
		C_result.setMsg ("ARM_ERR: check your maturities and rates array");
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	XL_readStrCell(XL_mmVsFut,C_mmVsFut," ARM_ERR: money market versus future: string expected",C_result);
	XL_readStrCell(XL_swapVsFut,C_swapVsFut," ARM_ERR: swap versus future: string expected",C_result);
	XL_readStrCellWD(XL_raw,C_raw,"P"," ARM_ERR: raw: string expected",C_result);
	XL_readStrCellWD(XL_ccy,C_ccy,"DEFAULT"," ARM_ERR: currency: string expected",C_result);

	if(C_ccy == "DEFAULT")
	{
		ARM_result ccyres;
		ARMLOCAL_ARM_GetDefaultCurrency (ccyres);
		if(ccyres.getRetCode () != ARM_OK)
		{
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}
		else
		{
			C_ccy = ccyres.getString ();
		}
	}

	if((mmVsFutId = ARM_ConvMktType (C_mmVsFut, C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	if((swapVsFutId = ARM_ConvMktType (C_swapVsFut, C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	rawId = ARM_ConvCvMethod (C_raw);

	long retCode;
	long objId;
	CCString prevClass;
	
	CCString curClass = LOCAL_ZERO_CURVE_CUBDIFF_CLASS;
	CCString stringId = GetLastCurCellEnvValue ();
	
	if(!stringId)
	{
		retCode = ARMLOCAL_zcswapcubdiff (C_date, C_matu, C_rate, (long)mmVsFutId, 
									 (long)swapVsFutId, (long)rawId, interpId,
									 C_ccy, C_result);

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
			retCode = ARMLOCAL_zcswapcubdiff (C_date, C_matu, C_rate, (long)mmVsFutId, 
									     (long)swapVsFutId, (long)rawId, interpId,
								         C_ccy, C_result, objId);

			if(retCode == ARM_OK)
			{
				LocalSetCurCellEnvValue (curClass, objId); 

				stringId = LocalMakeObjectId (objId, curClass);
			}
		}
		else
		{
			FreeCurCellContent ();
			retCode = ARMLOCAL_zcswapcubdiff (C_date, C_matu, C_rate, (long)mmVsFutId, 
									     (long)swapVsFutId, (long)rawId, interpId,
									     C_ccy, C_result);
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_ZCSWAPCUBDIFF" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}



__declspec(dllexport) LPXLOPER WINAPI Local_PXL_ZCSWAPCUBDIFF (LPXLOPER XL_date,
															   LPXLOPER XL_matuRate,
															   LPXLOPER XL_mmVsFut,
															   LPXLOPER XL_swapVsFut,
															   LPXLOPER XL_raw,
															   LPXLOPER XL_ccy)
{
	ADD_LOG("Local_PXL_ZCSWAPCUBDIFF ");
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();

	// C variable
	double C_date;
	VECTOR<CCString> C_matuRate;
	VECTOR<CCString> C_matu;
	VECTOR<double> C_rate; 

	CCString C_mmVsFut;
	long mmVsFutId;

	CCString C_swapVsFut;
	long swapVsFutId;

	CCString C_raw;
	long rawId;

	CCString C_ccy;

	long interpId;
	
	if((interpId = ARM_ConvInterpMethod ("L", C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}
	
	// error
	static int error;
	static char* reason = "";

	XL_readNumCell(XL_date,C_date," ARM_ERR: as of date: date expected",C_result);	
	XL_readStrVector(XL_matuRate,C_matuRate," ARM_ERR: maturities and rates: array expected",DOUBLE_TYPE,C_result);

	double tmp_rate = -1;
	CCString tmp_matu ("NULL");

	if(C_matuRate.size () % 2 != 0)
	{
		C_result.setMsg ("ARM_ERR: check your maturities and rates array");
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	int real_size = C_matuRate.size () / 2;
	double tmp;
	int j = 0;

	for(int i = 0; i < real_size; i++)
	{
		C_matu.push_back (C_matuRate[j]);
		if((sscanf ((const char*)C_matuRate[j+1], "%lf", &tmp) != 1))
		{
			C_result.setMsg ("ARM_ERR: check your maturities and rates array");
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}
		C_rate.push_back (tmp);
		j += 2;
	}

	CCString last = C_matu[C_matu.size () - 1];
	if((!last) || (last[0] == ' '))
	{
		C_result.setMsg ("ARM_ERR: check your maturities and rates array");
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	XL_readStrCell(XL_mmVsFut,C_mmVsFut," ARM_ERR: money market versus future: string expected",C_result);
	XL_readStrCell(XL_swapVsFut,C_swapVsFut," ARM_ERR: swap versus future: string expected",C_result);
	XL_readStrCellWD(XL_raw,C_raw,"P"," ARM_ERR: raw: string expected",C_result);
	XL_readStrCellWD(XL_ccy,C_ccy,"DEFAULT"," ARM_ERR: currency: string expected",C_result);

	if(C_ccy == "DEFAULT")
	{
		ARM_result ccyres;
		ARMLOCAL_ARM_GetDefaultCurrency (ccyres);
		if(ccyres.getRetCode () != ARM_OK)
		{
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}
		else
		{
			C_ccy = ccyres.getString ();
		}
	}

	if((mmVsFutId = ARM_ConvMktType (C_mmVsFut, C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	if((swapVsFutId = ARM_ConvMktType (C_swapVsFut, C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	rawId = ARM_ConvCvMethod (C_raw);

	long retCode;
	long objId;

	CCString curClass = LOCAL_ZERO_CURVE_CUBDIFF_CLASS;
	CCString stringId;

	retCode = ARMLOCAL_zcswapcubdiff (C_date, C_matu, C_rate, (long)mmVsFutId, 
								 (long)swapVsFutId, (long)rawId, interpId,
								 C_ccy, C_result);

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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_PXL_ZCSWAPCUBDIFF" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}



__declspec(dllexport) LPXLOPER WINAPI Local_ZCSWAPSPLSUM (LPXLOPER XL_date,
														   LPXLOPER XL_matuRate,
														   LPXLOPER XL_mmVsFut,
														   LPXLOPER XL_swapVsFut,
														   LPXLOPER XL_raw,
														   LPXLOPER XL_ccy)
{
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();

	// C variable
	double C_date;
	VECTOR<CCString> C_matuRate;
	VECTOR<CCString> C_matu;
	VECTOR<double> C_rate; 

	CCString C_mmVsFut;
	long mmVsFutId;

	CCString C_swapVsFut;
	long swapVsFutId;

	CCString C_raw;
	long rawId;

	CCString C_ccy;

	long interpId;
	
	if((interpId = ARM_ConvInterpMethod ("L", C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}
	
	// error
	static int error;
	static char* reason = "";

	XL_readNumCell(XL_date,C_date," ARM_ERR: as of date: date expected",C_result);
	XL_readStrVector(XL_matuRate,C_matuRate," ARM_ERR: maturities and rates: array expected",DOUBLE_TYPE,C_result);

	double tmp_rate = -1;
	CCString tmp_matu ("NULL");

	if(C_matuRate.size () % 2 != 0)
	{
		C_result.setMsg ("ARM_ERR: check your maturities and rates array");
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	int real_size = C_matuRate.size () / 2;
	double tmp;
	int j = 0;

	for(int i = 0; i < real_size; i++)
	{
		C_matu.push_back (C_matuRate[j]);
		if((sscanf ((const char*)C_matuRate[j+1], "%lf", &tmp) != 1))
		{
			C_result.setMsg ("ARM_ERR: check your maturities and rates array");
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}
		C_rate.push_back (tmp);
		j += 2;
	}

	CCString last = C_matu[C_matu.size () - 1];
	if((!last) || (last[0] == ' '))
	{
		C_result.setMsg ("ARM_ERR: check your maturities and rates array");
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	XL_readStrCell(XL_mmVsFut,C_mmVsFut," ARM_ERR: money market versus future: string expected",C_result);
	XL_readStrCell(XL_swapVsFut,C_swapVsFut," ARM_ERR: swap versus future: string expected",C_result);
	XL_readStrCellWD(XL_raw,C_raw,"P"," ARM_ERR: raw: string expected",C_result);
	XL_readStrCellWD(XL_ccy,C_ccy,"DEFAULT"," ARM_ERR: currency: string expected",C_result);

	if(C_ccy == "DEFAULT")
	{
		ARM_result ccyres;
		ARMLOCAL_ARM_GetDefaultCurrency (ccyres);
		if(ccyres.getRetCode () != ARM_OK)
		{
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}
		else
		{
			C_ccy = ccyres.getString ();
		}
	}

	if((mmVsFutId = ARM_ConvMktType (C_mmVsFut, C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	if((swapVsFutId = ARM_ConvMktType (C_swapVsFut, C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	rawId = ARM_ConvCvMethod (C_raw);

	long retCode;
	long objId;
	CCString prevClass;
	
	CCString curClass = LOCAL_ZERO_CURVE_SPLSUM_CLASS;
	CCString stringId = GetLastCurCellEnvValue ();
	
	if(!stringId)
	{
		retCode = ARMLOCAL_zcswapsplsum (C_date, C_matu, C_rate, (long)mmVsFutId, 
									 (long)swapVsFutId, (long)rawId, interpId,
									 C_ccy, C_result);

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
			retCode = ARMLOCAL_zcswapsplsum (C_date, C_matu, C_rate, (long)mmVsFutId, 
									     (long)swapVsFutId, (long)rawId, interpId,
								         C_ccy, C_result, objId);

			if(retCode == ARM_OK)
			{
				LocalSetCurCellEnvValue (curClass, objId); 

				stringId = LocalMakeObjectId (objId, curClass);
			}
		}
		else
		{
			FreeCurCellContent ();
			retCode = ARMLOCAL_zcswapsplsum (C_date, C_matu, C_rate, (long)mmVsFutId, 
									     (long)swapVsFutId, (long)rawId, interpId,
									     C_ccy, C_result);
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_ZCSWAPSPLSUM" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}



__declspec(dllexport) LPXLOPER WINAPI Local_PXL_ZCSWAPSPLSUM (LPXLOPER XL_date,
															   LPXLOPER XL_matuRate,
															   LPXLOPER XL_mmVsFut,
															   LPXLOPER XL_swapVsFut,
															   LPXLOPER XL_raw,
															   LPXLOPER XL_ccy)
{
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();

	// C variable
	double C_date;
	VECTOR<CCString> C_matuRate;
	VECTOR<CCString> C_matu;
	VECTOR<double> C_rate; 

	CCString C_mmVsFut;
	long mmVsFutId;

	CCString C_swapVsFut;
	long swapVsFutId;

	CCString C_raw;
	long rawId;

	CCString C_ccy;

	long interpId;
	
	if((interpId = ARM_ConvInterpMethod ("L", C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}
	
	// error
	static int error;
	static char* reason = "";

	XL_readNumCell(XL_date,C_date," ARM_ERR: as of date: date expected",C_result);	
	XL_readStrVector(XL_matuRate,C_matuRate," ARM_ERR: maturities and rates: array expected",DOUBLE_TYPE,C_result);

	double tmp_rate = -1;
	CCString tmp_matu ("NULL");

	if(C_matuRate.size () % 2 != 0)
	{
		C_result.setMsg ("ARM_ERR: check your maturities and rates array");
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	int real_size = C_matuRate.size () / 2;
	double tmp;
	int j = 0;

	for(int i = 0; i < real_size; i++)
	{
		C_matu.push_back (C_matuRate[j]);
		if((sscanf ((const char*)C_matuRate[j+1], "%lf", &tmp) != 1))
		{
			C_result.setMsg ("ARM_ERR: check your maturities and rates array");
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}
		C_rate.push_back (tmp);
		j += 2;
	}

	CCString last = C_matu[C_matu.size () - 1];
	if((!last) || (last[0] == ' '))
	{
		C_result.setMsg ("ARM_ERR: check your maturities and rates array");
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	XL_readStrCell(XL_mmVsFut,C_mmVsFut," ARM_ERR: money market versus future: string expected",C_result);
	XL_readStrCell(XL_swapVsFut,C_swapVsFut," ARM_ERR: swap versus future: string expected",C_result);
	XL_readStrCellWD(XL_raw,C_raw,"P"," ARM_ERR: raw: string expected",C_result);
	XL_readStrCellWD(XL_ccy,C_ccy,"DEFAULT"," ARM_ERR: currency: string expected",C_result);

	if(C_ccy == "DEFAULT")
	{
		ARM_result ccyres;
		ARMLOCAL_ARM_GetDefaultCurrency (ccyres);
		if(ccyres.getRetCode () != ARM_OK)
		{
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}
		else
		{
			C_ccy = ccyres.getString ();
		}
	}

	if((mmVsFutId = ARM_ConvMktType (C_mmVsFut, C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	if((swapVsFutId = ARM_ConvMktType (C_swapVsFut, C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	rawId = ARM_ConvCvMethod (C_raw);

	long retCode;
	long objId;

	CCString curClass = LOCAL_ZERO_CURVE_SPLSUM_CLASS;
	CCString stringId;

	retCode = ARMLOCAL_zcswapsplsum (C_date, C_matu, C_rate, (long)mmVsFutId, 
								 (long)swapVsFutId, (long)rawId, interpId,
								 C_ccy, C_result);

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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_PXL_ZCSWAPSPLSUM" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}



__declspec(dllexport) LPXLOPER WINAPI Local_ZCINTSMOOTH (LPXLOPER XL_matu,
														 LPXLOPER XL_rate,
														 LPXLOPER XL_aDate,
														 LPXLOPER XL_meth,
														 LPXLOPER XL_lambda,
														 LPXLOPER XL_prec)
{
	ADD_LOG("Local_ZCINTSMOOTH ");
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();

	// C variable
	VECTOR<double> C_matu;
	VECTOR<double> C_rate; 
	double C_meth;
	double C_aDate;
	
	double C_lambda;
	double C_lambda_default = 0.0;

	double C_prec;
	double C_prec_default = 0.0;

	// error
	static int error;
	static char* reason = "";

	XL_readNumVector(XL_matu,C_matu," ARM_ERR: maturities: array of numeric expected",C_result);
	XL_readNumVector(XL_rate,C_rate," ARM_ERR: rates: array of numeric expected",C_result);
	XL_readNumCell(XL_aDate,C_aDate," ARM_ERR: as of date: date expected",C_result);
	XL_readNumCell(XL_meth,C_meth," ARM_ERR: method: numeric expected",C_result);
	XL_readNumCellWD(XL_lambda,C_lambda,C_lambda_default," ARM_ERR: lambda: numeric expected",C_result);
	XL_readNumCellWD(XL_prec,C_prec,C_prec_default," ARM_ERR: precision: numeric expected",C_result);

	long retCode;
	long objId;
	CCString prevClass;
	
	CCString curClass = LOCAL_ZERO_CURVE_INTSMO_CLASS;
	CCString stringId = GetLastCurCellEnvValue ();
	
	if(!stringId)
	{
		retCode = ARMLOCAL_ZCINTSMOOTH (C_matu, C_rate, C_aDate, (long)C_meth, C_lambda, (long)C_prec, C_result);

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
			retCode = ARMLOCAL_ZCINTSMOOTH (C_matu, C_rate, C_aDate, (long)C_meth, C_lambda, (long)C_prec, C_result, objId);

			if(retCode == ARM_OK)
			{
				LocalSetCurCellEnvValue (curClass, objId); 

				stringId = LocalMakeObjectId (objId, curClass);
			}
		}
		else
		{
			FreeCurCellContent ();
			retCode = ARMLOCAL_ZCINTSMOOTH (C_matu, C_rate, C_aDate, (long)C_meth, C_lambda, (long)C_prec, C_result);
		
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_ZCINTSMOOTH" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}



__declspec(dllexport) LPXLOPER WINAPI Local_PXL_ZCINTSMOOTH (LPXLOPER XL_matu,
															 LPXLOPER XL_rate,
															 LPXLOPER XL_aDate,
															 LPXLOPER XL_meth,
															 LPXLOPER XL_lambda,
															 LPXLOPER XL_prec)
{
	ADD_LOG("Local_PXL_ZCINTSMOOTH ");
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();

	// C variable
	VECTOR<double> C_matu;
	VECTOR<double> C_rate; 
	double C_meth;
	double C_aDate;

	double C_lambda;
	double C_lambda_default = 0.0;

	double C_prec;
	double C_prec_default = 0.0;

	// error
	static int error;
	static char* reason = "";

	XL_readNumVector(XL_matu,C_matu," ARM_ERR: maturities: array of numeric expected",C_result);
	XL_readNumVector(XL_rate,C_rate," ARM_ERR: rates: array of numeric expected",C_result);
	XL_readNumCell(XL_meth,C_meth," ARM_ERR: method: numeric expected",C_result);
	XL_readNumCell(XL_aDate,C_aDate," ARM_ERR: as of date: date expected",C_result);
	XL_readNumCellWD(XL_lambda,C_lambda,C_lambda_default," ARM_ERR: lambda: numeric expected",C_result);
	XL_readNumCellWD(XL_prec,C_prec,C_prec_default," ARM_ERR: precision: numeric expected",C_result);

	long retCode;
	long objId;
	
	CCString curClass = LOCAL_ZERO_CURVE_INTSMO_CLASS;
	CCString stringId;
	
	retCode = ARMLOCAL_ZCINTSMOOTH (C_matu, C_rate, C_aDate, (long)C_meth, C_lambda, (long)C_prec, C_result);

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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_PXL_ZCINTSMOOTH" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}

__declspec(dllexport) LPXLOPER WINAPI Local_GetZCFromSummit (LPXLOPER XL_index,
															 LPXLOPER XL_currency,
															 LPXLOPER XL_cvName,
															 LPXLOPER XL_aSdate,
															 LPXLOPER XL_interp)
{
	ADD_LOG("Local_GetZCFromSummit ");
//	ARM_BEGIN();
	
	// return
	static XLOPER XL_result;
	ARM_result C_result;

    /// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();

    // C variable
	CCString C_index;
	CCString C_currency;
	CCString C_cvName;
	double C_aSdate;
	
	CCString C_interp;
	long interpId;

	// error
	static int error;
	static char* reason = "";

	XL_readStrCell(XL_index,C_index," ARM_ERR: index: string expected",C_result);
	XL_readStrCell(XL_currency,C_currency," ARM_ERR: currency: string expected",C_result);
	XL_readStrCell(XL_cvName,C_cvName," ARM_ERR: curve name: string expected",C_result);
	XL_readNumCell(XL_aSdate,C_aSdate," ARM_ERR: as of date: date expected",C_result);
	XL_readStrCellWD(XL_interp,C_interp,"L"," ARM_ERR: interpolation mode: string expected",C_result);

	C_index.toUpper ();
	C_currency.toUpper ();
	C_cvName.toUpper ();

	if((interpId = ARM_ConvInterpMethod (C_interp, C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	long retCode;
	long objId;
	CCString prevClass;

	CCString curClass = LOCAL_ZERO_CURVE_LIN_CLASS;
	CCString stringId = GetLastCurCellEnvValue ();
	
	if(!stringId)
	{
		retCode = ARMLOCAL_GetZCFromSummit (C_index, C_currency, C_cvName, C_aSdate, interpId, C_result);

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
			retCode = ARMLOCAL_GetZCFromSummit (C_index, C_currency, C_cvName, C_aSdate, interpId, C_result, objId);

			if(retCode == ARM_OK)
			{
				LocalSetCurCellEnvValue (curClass, objId); 

				stringId = LocalMakeObjectId (objId, curClass);
			}
		}
		else
		{
			FreeCurCellContent ();
			retCode = ARMLOCAL_GetZCFromSummit (C_index, C_currency, C_cvName, C_aSdate, interpId, C_result);
		
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_GetZCFromSummit" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}



__declspec(dllexport) LPXLOPER WINAPI Local_GetZCFromCalypso (LPXLOPER XL_index,
															 LPXLOPER XL_currency,
															 LPXLOPER XL_term,
															 LPXLOPER XL_pricingEnv,
															 LPXLOPER XL_aSdate,
															 LPXLOPER XL_interp,
															 LPXLOPER XL_forceCurveName,
															 LPXLOPER XL_xmlFileName
															 )
{
	ADD_LOG("Local_GetZCFromCalypso ");
	static XLOPER XL_result;
	ARM_result C_result;

    /// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	FreeCurCellErr ();
	ARM_NOCALCIFWIZ();

   
	long interpId;
	static int error;
	static char* reason = "";

	std::string index ; ExcelTools::convert(XL_index,"",index); 
	std::string currency ; ExcelTools::convert(XL_currency,"",currency); 
	std::string term; ExcelTools::convert(XL_term,"",term); 
	std::string pricingEnv; ExcelTools::convert(XL_pricingEnv,"MO",pricingEnv); 
	std::string interp ; ExcelTools::convert(XL_interp,"L",interp); 
	std::string forceCurveName; ExcelTools::convert(XL_forceCurveName,"",forceCurveName); 
	std::string xmlFileName; ExcelTools::convert(XL_xmlFileName,"",xmlFileName); 

//	if (interp=="") interp="L"; 

	if((interpId = ARM_ConvInterpMethod (interp.c_str(), C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	
	long objId;
	CCString prevClass;

	CCString curClass = LOCAL_ZERO_CURVE_LIN_CLASS;
	CCString stringId = GetLastCurCellEnvValue ();
	

	ARM_Date AsOf ;
	ExcelTools::convert(XL_aSdate,AsOf); 

	// case 1 :: cell was empty 
	if(!stringId)
	{
		long objId = ARMLOCAL_GetZCFromCalypso(AsOf,index, currency, term,pricingEnv, forceCurveName,interpId, xmlFileName /** ,-1 **/ );
		LocalSetCurCellEnvValue (curClass, objId); 
		stringId = LocalMakeObjectId (objId, curClass);
		ExcelTools::convert(CCSTringToSTLString(stringId),&XL_result); 
		return &XL_result ;
	}
	
	//  cell not empty: 
	prevClass = LocalGetStringObjectClass (stringId);
	objId = LocalGetNumObjectId (stringId);
		
	//	cell was containing the same kind of object 
	if(curClass == prevClass)
	{
		objId= ARMLOCAL_GetZCFromCalypso(AsOf,index, currency, term,pricingEnv, forceCurveName,interpId, xmlFileName, objId );
		LocalSetCurCellEnvValue (curClass, objId); 
		stringId = LocalMakeObjectId (objId, curClass);
		std::string toto = CCSTringToSTLString(stringId); 
		ExcelTools::convert(toto,&XL_result); 
		return &XL_result ;
	}
	
	// cell was containing other kind of object
	FreeCurCellContent ();
	objId= ARMLOCAL_GetZCFromCalypso(AsOf,index, currency, term,pricingEnv, forceCurveName,interpId, xmlFileName /** ,-1 **/ );
	LocalSetCurCellEnvValue (curClass, objId); 
	stringId = LocalMakeObjectId (objId, curClass);
	ExcelTools::convert(CCSTringToSTLString(stringId),&XL_result); 

	return &XL_result ;
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_GetZCFromCalypso" )

	/// return the result as an LPXLOPER
	return &XL_result;
}


__declspec(dllexport) LPXLOPER WINAPI Local_PXL_GetZCFromSummit (LPXLOPER XL_index,
																 LPXLOPER XL_currency,
																 LPXLOPER XL_cvName,
																 LPXLOPER XL_aSdate,
																 LPXLOPER XL_interp)
{
	ADD_LOG("Local_PXL_GetZCFromSummit ");
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();

    // C variable
	CCString C_index;
	CCString C_currency;
	CCString C_cvName;
	double C_aSdate;
	
	CCString C_interp;
	long interpId;

	// error
	static int error;
	static char* reason = "";

	XL_readStrCell(XL_index,C_index," ARM_ERR: index: string expected",C_result);
	XL_readStrCell(XL_currency,C_currency," ARM_ERR: currency: string expected",C_result);
	XL_readStrCell(XL_cvName,C_cvName," ARM_ERR: curve name: string expected",C_result);
	XL_readNumCell(XL_aSdate,C_aSdate," ARM_ERR: as of date: date expected",C_result);
	XL_readStrCellWD(XL_interp,C_interp,"L"," ARM_ERR: interpolation mode: string expected",C_result);

	C_index.toUpper ();
	C_currency.toUpper ();
	C_cvName.toUpper ();

	if((interpId = ARM_ConvInterpMethod (C_interp, C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	long retCode;
	long objId;
	CCString prevClass;

	CCString curClass = LOCAL_ZERO_CURVE_LIN_CLASS;
	CCString stringId;
	
	retCode = ARMLOCAL_GetZCFromSummit (C_index, C_currency, C_cvName, C_aSdate, interpId, C_result);

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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_PXL_GetZCFromSummit" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}

__declspec(dllexport) LPXLOPER WINAPI Local_GetMaturitiesFromZC(LPXLOPER XL_cvName)
{
	ADD_LOG("Local_GetMaturitiesFromZC");
	// return
	static XLOPER XL_result;
	ARM_result C_result;
	LPXLOPER pxArray;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
		ARM_NOCALCIFWIZ();
		
		CCString C_cvName;
		long ObjId;
		long retCode;

		// error
		static int error;
		static char* reason = "";

		XL_readStrCell(XL_cvName,C_cvName," ARM_ERR: curve name: string expected",C_result);
		
		ObjId = LocalGetNumObjectId(C_cvName);
		retCode = ARMLOCAL_GetMaturitiesFromZC(C_result,ObjId);

		if(retCode == ARM_OK)
		{
			const VECTOR<CCString> vMaturities = C_result.getStringVector();
			int nbcolumns = 1;
			int nbrows = vMaturities.size();
			FreeCurCellErr ();

			XL_result.xltype= xltypeMulti ;
			XL_result.val.array.columns = nbcolumns;
			XL_result.val.array.rows = nbrows;
			XL_result.val.array.lparray = pxArray = (LPXLOPER)GlobalAlloc(GMEM_ZEROINIT, nbrows * nbcolumns * sizeof(XLOPER));
		
			for (int j = 0 ; j < nbrows ; j++)
			{	
				pxArray[XL_Coordonnate2Rank (j, 0, nbcolumns)].xltype = xltypeStr;
				pxArray[XL_Coordonnate2Rank (j, 0, nbcolumns)].val.str= XL_StrC2StrPascal(vMaturities[j]);
				pxArray[XL_Coordonnate2Rank (j, 0, nbcolumns)].xltype |= xlbitDLLFree;
			}
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_GetMaturitiesFromZC" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;


}

__declspec(dllexport) LPXLOPER WINAPI Local_GetInitialCurveFromCalypso(LPXLOPER XL_index,
																	   LPXLOPER XL_currency,
																	   LPXLOPER XL_term,
																	   LPXLOPER XL_pricingEnv,
																	   LPXLOPER XL_Date,
																	   LPXLOPER XL_forceCurveName,
																	   LPXLOPER XL_xmlFile,
																	   LPXLOPER XL_convAdj)
{
	ADD_LOG("Local_GetInitialCurveFromCalypso");

	static XLOPER XL_result;
	LPXLOPER pxArray ;

	try
	{
		ARM_NOCALCIFWIZ();
		
		std::string index ; ExcelTools::convert(XL_index ,"",index);
		std::string currency ; ExcelTools::convert(XL_currency,"",currency);
		std::string term ; ExcelTools::convert(XL_term,"",term);
		std::string pricingEnv ; ExcelTools::convert(XL_pricingEnv,"",pricingEnv);
		std::string forceCurveName ; ExcelTools::convert(XL_forceCurveName,"",forceCurveName);
		std::string xmlFile ; ExcelTools::convert(XL_xmlFile,"",xmlFile);
		std::string convAdj ; ExcelTools::convert(XL_convAdj,"",convAdj);
		ARM_Date AsOf ; ExcelTools::convert(XL_Date,AsOf);

		std::vector<std::string> matus;
		std::vector<double> yields;

		bool doAdj;
		doAdj = ARM_ConvYesOrNo(convAdj);
		ARMLOCAL_GetInitialCurveFromCalypso(AsOf,pricingEnv,index,currency,term,forceCurveName,xmlFile,doAdj,matus,yields); 

		int nbcolumns = 2;
		int nbrows = matus.size();
		FreeCurCellErr ();

		XL_result.xltype= xltypeMulti ;
		XL_result.val.array.columns = nbcolumns;
		XL_result.val.array.rows = nbrows;
		XL_result.val.array.lparray = pxArray = (LPXLOPER)GlobalAlloc(GMEM_ZEROINIT, nbrows * nbcolumns * sizeof(XLOPER));
		
		char * tmp_str ;
		CCString tmpCCstring ;
	
		for (int j = 0 ; j < nbrows ; j++)
		{	
			tmp_str  =(char *) matus[j].c_str();
			tmpCCstring = CCString(tmp_str);

			pxArray[XL_Coordonnate2Rank (j, 0, nbcolumns)].xltype = xltypeStr;
			pxArray[XL_Coordonnate2Rank (j, 0, nbcolumns)].val.str= XL_StrC2StrPascal(tmpCCstring);
			pxArray[XL_Coordonnate2Rank (j, 0, nbcolumns)].xltype |= xlbitDLLFree;

			pxArray[XL_Coordonnate2Rank (j, 1, nbcolumns)].xltype = xltypeNum;
			pxArray[XL_Coordonnate2Rank (j, 1, nbcolumns)].val.num = yields[j];
		}
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

__declspec(dllexport) LPXLOPER WINAPI Local_GetInitialCurveFromSummit (LPXLOPER XL_index,
																	   LPXLOPER XL_currency,
																	   LPXLOPER XL_cvName,
																	   LPXLOPER XL_aSdate,
																	   LPXLOPER XL_convAdj)
{
	ADD_LOG("Local_GetInitialCurveFromSummit ");
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
	CCString C_index;
	CCString C_currency;
	CCString C_cvName;
	double C_aSdate;
	VECTOR<CCString> matu;
	VECTOR<double> rate;

	CCString C_convAdj;
	long adjOrNotId;

	// error
	static int error;
	static char* reason = "";

	XL_readStrCell(XL_index,C_index," ARM_ERR: index: string expected",C_result);
	XL_readStrCell(XL_currency,C_currency," ARM_ERR: currency: string expected",C_result);
	XL_readStrCell(XL_cvName,C_cvName," ARM_ERR: curve name: string expected",C_result);
	XL_readNumCell(XL_aSdate,C_aSdate," ARM_ERR: as of date: date expected",C_result);
	XL_readStrCellWD(XL_convAdj,C_convAdj,"YES"," ARM_ERR: convexity adjustment: string expected",C_result);

	C_index.toUpper ();
	C_currency.toUpper ();
	C_cvName.toUpper ();

	if((adjOrNotId = ARM_ConvYesOrNo (C_convAdj, C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	long retCode;

	retCode = ARMLOCAL_GetInitialCurveFromSummit (C_index, C_currency, C_cvName, C_aSdate, adjOrNotId, &matu, &rate, C_result);
//	retCode = ARMLOCAL_GetInitialCurveFromSummit (C_index, C_currency, C_cvName, C_aSdate, C_result);

	// Dans le long du C_result, on sette le isETK
	if(retCode == ARM_OK)
	{
		long result;

		if (C_result.getLong() >= ETKRETRIEVER)
		{
			result = ARM_OK;
		}
		else
		{
			if (strncmp((const char*)C_index,"BS",2) == 0)
				result = LocalExtractCurveFromFileBS (C_result.getString(), matu, rate);
			else
				result = LocalExtractCurveFromFileMO (C_result.getString(), matu, rate,adjOrNotId);
		}
		
		if(result == ARM_OK)
		{
			int isVBA = IsVBA();

			int nbrows;
			if (C_result.getLong() >= 1)
				nbrows= matu.size ();
			else
				nbrows = matu.size ()-1;
			int nbcolumns = 2;

            /// add these additional lines 
			/// to display blank lines
			const int additionalLinesNb = 100;
			bool fillWithBlank = true;

            int nbRowswithblank;
			
			if (isVBA == 1)
				nbRowswithblank = nbrows;
			else
				nbRowswithblank = nbrows + additionalLinesNb;
		
			FreeCurCellErr ();
			XL_result.xltype = xltypeMulti;
			XL_result.val.array.columns = nbcolumns;
			XL_result.val.array.rows = nbRowswithblank; 
			XL_result.val.array.lparray = pxArray = (LPXLOPER)GlobalAlloc (GMEM_ZEROINIT, nbRowswithblank * nbcolumns * sizeof (XLOPER));

			for(int i = 0; i < nbrows; i++)
			{
				if (((C_index == "IFRF") || (C_index == "EMU"))
					&& 
					((string(matu[i]).find(".") != -1) || (string(matu[i]).find("/") != -1)))
				{
					pxArray[XL_Coordonnate2Rank (i, 0, nbcolumns)].xltype = xltypeNum;
					pxArray[XL_Coordonnate2Rank (i, 0, nbcolumns)].val.num = Local_ARMDATE2XLDATE (matu[i]);
				}
				else
				{
					pxArray[XL_Coordonnate2Rank (i, 0, nbcolumns)].xltype = xltypeStr;
					pxArray[XL_Coordonnate2Rank (i, 0, nbcolumns)].val.str = XL_StrC2StrPascal (matu[i]);
				}
				pxArray[XL_Coordonnate2Rank (i, 0, nbcolumns)].xltype |= xlbitDLLFree;
				pxArray[XL_Coordonnate2Rank (i, 1, nbcolumns)].xltype = xltypeNum;
				pxArray[XL_Coordonnate2Rank (i, 1, nbcolumns)].val.num = rate[i]; 
			}

            for(; i < nbRowswithblank; i++)
			{

                if( fillWithBlank )
			    {
				    pxArray[XL_Coordonnate2Rank (i, 0, nbcolumns)].xltype = xltypeStr;
				    pxArray[XL_Coordonnate2Rank (i, 0, nbcolumns)].val.str= XL_StrC2StrPascal("");
				    pxArray[XL_Coordonnate2Rank (i, 0, nbcolumns)].xltype |= xlbitDLLFree;
                    pxArray[XL_Coordonnate2Rank (i, 1, nbcolumns)].xltype = xltypeStr;
				    pxArray[XL_Coordonnate2Rank (i, 1, nbcolumns)].val.str= XL_StrC2StrPascal("");
				    pxArray[XL_Coordonnate2Rank (i, 1, nbcolumns)].xltype |= xlbitDLLFree;
			    }
			    else
			    {
				    pxArray[XL_Coordonnate2Rank (i, 0, nbcolumns)].xltype = xltypeNil;
                    pxArray[XL_Coordonnate2Rank (i, 1, nbcolumns)].xltype = xltypeNil;
			    }
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_GetInitialCurveFromSummit" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


__declspec(dllexport) LPXLOPER WINAPI Local_ARM_CreateZCFromSummit (LPXLOPER XL_index,
																	LPXLOPER XL_currency,
																	LPXLOPER XL_cvName,
																	LPXLOPER XL_aSdate,
																	LPXLOPER XL_convAdj,
																	LPXLOPER XL_raw,
																	LPXLOPER XL_swapFrq)
{
	ADD_LOG("Local_ARM_CreateZCFromSummit ");
//	ARM_BEGIN();
	
	// return
	static XLOPER XL_result;
	ARM_result C_result;

    /// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();

    // C variable
	CCString C_index;
	CCString C_currency;
	CCString C_cvName;
	double C_aSdate;
	
	CCString C_convAdj;
	long adjOrNotId;

	CCString C_raw;

	CCString C_swapFrq;
	long swapFrqId;

	// error
	static int error;
	static char* reason = "";

	XL_readStrCell(XL_index,C_index," ARM_ERR: index: string expected",C_result);
	XL_readStrCell(XL_currency,C_currency," ARM_ERR: currency: string expected",C_result);
	XL_readStrCell(XL_cvName,C_cvName," ARM_ERR: curve name: string expected",C_result);
	XL_readNumCell(XL_aSdate,C_aSdate," ARM_ERR: as of date: date expected",C_result);
	XL_readStrCellWD(XL_convAdj,C_convAdj,"YES"," ARM_ERR: convexity adjustment: string expected",C_result);
	XL_readStrCellWD(XL_raw,C_raw,"DEFAULT"," ARM_ERR: FWD method: string expected",C_result);
	XL_readStrCellWD(XL_swapFrq,C_swapFrq,"-1"," ARM_ERR: swap frequency: string expected",C_result);

	C_index.toUpper ();
	C_currency.toUpper ();
	C_cvName.toUpper ();

	if((adjOrNotId = ARM_ConvYesOrNo (C_convAdj, C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	if ((swapFrqId = ARM_ConvFrequency(C_swapFrq, C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	long retCode;
	long objId;
	CCString prevClass;

	CCString curClass = LOCAL_ZERO_CURVE_LIN_CLASS;
	CCString stringId = GetLastCurCellEnvValue ();
	
	switchToETK();

	if(!stringId)
	{
		retCode = ARMLOCAL_CreateZCFromSummit (C_index, C_currency, C_cvName, C_aSdate, adjOrNotId, C_raw, swapFrqId, C_result);

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
			retCode = ARMLOCAL_CreateZCFromSummit (C_index, C_currency, C_cvName, C_aSdate, adjOrNotId, C_raw, swapFrqId, C_result, objId);

			if(retCode == ARM_OK)
			{
				LocalSetCurCellEnvValue (curClass, objId); 

				stringId = LocalMakeObjectId (objId, curClass);
			}
		}
		else
		{
			FreeCurCellContent ();
			retCode = ARMLOCAL_CreateZCFromSummit (C_index, C_currency, C_cvName, C_aSdate, adjOrNotId, C_raw, swapFrqId, C_result);
		
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_ARM_CreateZCFromSummit" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


__declspec(dllexport) LPXLOPER WINAPI Local_PXL_ARM_CreateZCFromSummit (LPXLOPER XL_index,
																		LPXLOPER XL_currency,
																		LPXLOPER XL_cvName,
																		LPXLOPER XL_aSdate,
																		LPXLOPER XL_convAdj,
																		LPXLOPER XL_raw,
																		LPXLOPER XL_swapFrq)
{
	ADD_LOG("Local_PXL_ARM_CreateZCFromSummit ");
//	ARM_BEGIN();
	
	// return
	static XLOPER XL_result;
	ARM_result C_result;

    /// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();

    // C variable
	CCString C_index;
	CCString C_currency;
	CCString C_cvName;
	double C_aSdate;
	
	CCString C_convAdj;
	long adjOrNotId;

	CCString C_raw;

	CCString C_swapFrq;
	long swapFrqId;

	// error
	static int error;
	static char* reason = "";

	XL_readStrCell(XL_index,C_index," ARM_ERR: index: string expected",C_result);
	XL_readStrCell(XL_currency,C_currency," ARM_ERR: currency: string expected",C_result);
	XL_readStrCell(XL_cvName,C_cvName," ARM_ERR: curve name: string expected",C_result);
	XL_readNumCell(XL_aSdate,C_aSdate," ARM_ERR: as of date: date expected",C_result);
	XL_readStrCellWD(XL_convAdj,C_convAdj,"YES"," ARM_ERR: convexity adjustment: string expected",C_result);
	XL_readStrCellWD(XL_raw,C_raw,"DEFAULT"," ARM_ERR: FWD method: string expected",C_result);
	XL_readStrCellWD(XL_swapFrq,C_swapFrq,"-1"," ARM_ERR: swap frequency: string expected",C_result);

	C_index.toUpper ();
	C_currency.toUpper ();
	C_cvName.toUpper ();

	if((adjOrNotId = ARM_ConvYesOrNo (C_convAdj, C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	if ((swapFrqId = ARM_ConvFrequency(C_swapFrq, C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	long retCode;
	long objId;
	CCString prevClass;

	CCString curClass = LOCAL_ZERO_CURVE_LIN_CLASS;
	CCString stringId;
	
	retCode = ARMLOCAL_CreateZCFromSummit (C_index, C_currency, C_cvName, C_aSdate, adjOrNotId, C_raw, swapFrqId, C_result);

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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_PXL_ARM_CreateZCFromSummit" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}



__declspec(dllexport) LPXLOPER WINAPI Local_ZCVSK (LPXLOPER XL_param,
												   LPXLOPER XL_date)
{
	ADD_LOG("Local_ZCVSK ");
	//ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();

	// C variable
	VECTOR<double> C_param;
	double C_date;

	// error
	static int error;
	static char* reason = "";

	XL_readNumVector(XL_param,C_param," ARM_ERR: parameters: array of numeric expected",C_result);
	XL_readNumCell(XL_date,C_date," ARM_ERR: as of date: date expected",C_result);

	long retCode;
	long objId;
	CCString prevClass;
	
	CCString curClass = LOCAL_ZERO_CURVE_VSK_CLASS;
	CCString stringId = GetLastCurCellEnvValue ();
	
	if(!stringId)
	{
		retCode = ARMLOCAL_zcvsk (C_param, C_date, C_result);

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
			retCode = ARMLOCAL_zcvsk (C_param, C_date, C_result, objId);

			if(retCode == ARM_OK)
			{
				LocalSetCurCellEnvValue (curClass, objId); 

				stringId = LocalMakeObjectId (objId, curClass);
			}
		}
		else
		{
			FreeCurCellContent ();
			retCode = ARMLOCAL_zcvsk (C_param, C_date, C_result);

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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_ZCVSK" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}

__declspec(dllexport) LPXLOPER WINAPI Local_PXL_ZCVSK (LPXLOPER XL_param,
													   LPXLOPER XL_date)
{
	ADD_LOG("Local_PXL_ZCVSK ");
	//ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();

	// C variable
	VECTOR<double> C_param;
	double C_date;

	// error
	static int error;
	static char* reason = "";

	XL_readNumVector(XL_param,C_param," ARM_ERR: parameters: array of numeric expected",C_result);
	XL_readNumCell(XL_date,C_date," ARM_ERR: as of date: date expected",C_result);

	long retCode;
	long objId;
	
	CCString curClass = LOCAL_ZERO_CURVE_VSK_CLASS;
	CCString stringId;
	
	retCode = ARMLOCAL_zcvsk (C_param, C_date, C_result);

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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_PXL_ZCVSK" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}



__declspec(dllexport) LPXLOPER WINAPI Local_ZCSPLICUB (LPXLOPER XL_matu,
													   LPXLOPER XL_rate,
													   LPXLOPER XL_meth,
													   LPXLOPER XL_date,
													   LPXLOPER XL_lastBucket)
{
	ADD_LOG("Local_ZCSPLICUB ");
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();

	// C variable
	VECTOR<double> C_matu;
	VECTOR<double> C_rate; 
	double C_meth;
	double C_meth_default = 0;
	double C_date;
	double C_lastBucket;
	double C_lastBucket_default = 1;

	// error
	static int error;
	static char* reason = "";

	XL_readNumVector(XL_matu,C_matu," ARM_ERR: maturities: array of numeric expected",C_result);
	XL_readNumVector(XL_rate,C_rate," ARM_ERR: rates: array of numeric expected",C_result);
	XL_readNumCellWD(XL_meth,C_meth,C_meth_default," ARM_ERR: method: numeric expected",C_result);
	XL_readNumCell(XL_date,C_date," ARM_ERR: as of date: date expected",C_result);
	XL_readNumCellWD(XL_lastBucket,C_lastBucket,C_lastBucket_default," ARM_ERR: last bucket: numeric expected",C_result);

	long retCode;
	long objId;
	CCString prevClass;
	
	CCString curClass = LOCAL_ZERO_CURVE_SPLCUB_CLASS;
	CCString stringId = GetLastCurCellEnvValue ();
	
	if(!stringId)
	{
		retCode = ARMLOCAL_zcsplicub (C_matu, C_rate, (long)C_meth, C_date, (long)C_lastBucket, C_result);

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
			retCode = ARMLOCAL_zcsplicub (C_matu, C_rate, (long)C_meth, C_date, (long)C_lastBucket, C_result, objId);

			if(retCode == ARM_OK)
			{
				LocalSetCurCellEnvValue (curClass, objId); 

				stringId = LocalMakeObjectId (objId, curClass);
			}
		}
		else
		{
			FreeCurCellContent ();
			retCode = ARMLOCAL_zcsplicub (C_matu, C_rate, (long)C_meth, C_date, (long)C_lastBucket, C_result);

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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_ZCSPLICUB" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}

__declspec(dllexport) LPXLOPER WINAPI Local_PXL_ZCSPLICUB (LPXLOPER XL_matu,
														   LPXLOPER XL_rate,
														   LPXLOPER XL_meth,
														   LPXLOPER XL_date,
														   LPXLOPER XL_lastBucket)
{
	ADD_LOG("Local_PXL_ZCSPLICUB ");
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();

	// C variable
	VECTOR<double> C_matu;
	VECTOR<double> C_rate; 
	double C_meth;
	double C_meth_default = 0;
	double C_date;
	double C_lastBucket;
	double C_lastBucket_default = 1;

	// error
	static int error;
	static char* reason = "";

	XL_readNumVector(XL_matu,C_matu," ARM_ERR: maturities: array of numeric expected",C_result);
	XL_readNumVector(XL_rate,C_rate," ARM_ERR: rates: array of numeric expected",C_result);
	XL_readNumCellWD(XL_meth,C_meth,C_meth_default," ARM_ERR: method: numeric expected",C_result);
	XL_readNumCell(XL_date,C_date," ARM_ERR: as of date: date expected",C_result);
	XL_readNumCellWD(XL_lastBucket,C_lastBucket,C_lastBucket_default," ARM_ERR: last bucket: numeric expected",C_result);

	long retCode;
	long objId;
	
	CCString curClass = LOCAL_ZERO_CURVE_SPLCUB_CLASS;
	CCString stringId;
	
	retCode = ARMLOCAL_zcsplicub (C_matu, C_rate, (long)C_meth, C_date, (long)C_lastBucket, C_result);

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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_PXL_ZCSPLICUB" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


__declspec(dllexport) LPXLOPER WINAPI Local_ZCSPL (LPXLOPER XL_param,
												   LPXLOPER XL_date)
{
	ADD_LOG("Local_ZCSPL ");
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

    /// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();

	// C variable
	VECTOR<double> C_param;
	double C_date;

	// error
	static int error;
	static char* reason = "";

	XL_readNumVector(XL_param,C_param," ARM_ERR: parameters: array of numeric expected",C_result);
	XL_readNumCell(XL_date,C_date," ARM_ERR: as of date: date expected",C_result);

	long retCode;
	long objId;
	CCString prevClass;
	
	CCString curClass = LOCAL_ZERO_CURVE_SPL_CLASS;
	CCString stringId = GetLastCurCellEnvValue ();
	
	if(!stringId)
	{
		retCode = ARMLOCAL_zcspl (C_param, C_date, C_result);

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
			retCode = ARMLOCAL_zcspl (C_param, C_date, C_result, objId);

			if(retCode == ARM_OK)
			{
				LocalSetCurCellEnvValue (curClass, objId); 

				stringId = LocalMakeObjectId (objId, curClass);
			}
		}
		else
		{
			FreeCurCellContent ();
			retCode = ARMLOCAL_zcspl (C_param, C_date, C_result);
		
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_ZCSPL" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


__declspec(dllexport) LPXLOPER WINAPI Local_PXL_ZCSPL (LPXLOPER XL_param,
													   LPXLOPER XL_date)
{
	ADD_LOG("Local_PXL_ZCSPL ");
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();

	// C variable
	VECTOR<double> C_param;
	double C_date;

	// error
	static int error;
	static char* reason = "";

	XL_readNumVector(XL_param,C_param," ARM_ERR: parameters: array of numeric expected",C_result);
	XL_readNumCell(XL_date,C_date," ARM_ERR: as of date: date expected",C_result);

	long retCode;
	long objId;

	CCString curClass = LOCAL_ZERO_CURVE_SPL_CLASS;
	CCString stringId;

	retCode = ARMLOCAL_zcspl (C_param, C_date, C_result);

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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_PXL_ZCSPL" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}



__declspec(dllexport) LPXLOPER WINAPI Local_ARM_BUMPCURVE (LPXLOPER XL_inCv,
														   LPXLOPER XL_epsilon,
														   LPXLOPER XL_meth)
{
	ADD_LOG("Local_ARM_BUMPCURVE ");
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();

	// C variable
	CCString C_cvString;
	long C_inCvId;
	double C_meth;
	double C_meth_default = 0;

	VECTOR<CCString> C_epsilon;
	double C_epsilon_dbl;
	long epsilonType;

	VECTOR<CCString> C_matu;
	VECTOR<double> C_rate; 

	// error
	static int error;
	static char* reason = "";

	XL_readStrCell(XL_inCv,C_cvString," ARM_ERR: curve id: object expected",C_result);
	XL_readNumOrStrVector(XL_epsilon,C_epsilon_dbl,C_epsilon,epsilonType," ARM_ERR: epsilon: array of string expected",C_result);
	XL_readNumCellWD(XL_meth,C_meth,C_meth_default," ARM_ERR: method: numeric expected (0 - 'Classic'; 1 - 'Forward')",C_result);

	if (epsilonType == XL_TYPE_DOUBLE)
	{
		C_rate.push_back (C_epsilon_dbl);
	}
	else
	{
		double tmp_rate = -1;
		CCString tmp_matu ("NULL");

		if(C_epsilon.size () % 2 != 0)
		{
			C_result.setMsg ("ARM_ERR: check your maturities and rates array");
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}

		int real_size = C_epsilon.size () / 2;
		double tmp;
		int j = 0;

		for(int i = 0; i < real_size; i++)
		{
			C_epsilon[j].toUpper();
			C_matu.push_back (C_epsilon[j]);
			if((sscanf ((const char*)C_epsilon[j+1], "%lf", &tmp) != 1))
			{
				C_result.setMsg ("ARM_ERR: check your maturities and rates array");
				ARM_ARG_ERR();
				return (LPXLOPER)&XL_result;
			}
			C_rate.push_back (tmp);
			j += 2;
		}

		CCString last = C_matu[C_matu.size () - 1];
		if((!last) || (last[0] == ' '))
		{
			C_result.setMsg ("ARM_ERR: check your maturities and rates array");
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}
	}

	C_inCvId =  LocalGetNumObjectId (C_cvString);

	long retCode;
	long objId;
	CCString prevClass;

	CCString curClass = LocalGetStringObjectClass(C_cvString);
	CCString stringId = GetLastCurCellEnvValue ();

	if(!stringId)
	{
		retCode = ARMLOCAL_BumpCurve (C_inCvId, C_matu, C_rate, C_meth, C_result);

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
			retCode = ARMLOCAL_BumpCurve (C_inCvId, C_matu, C_rate, C_meth, C_result, objId);

			if(retCode == ARM_OK)
			{
				LocalSetCurCellEnvValue (curClass, objId); 

				stringId = LocalMakeObjectId (objId, curClass);
			}
		}
		else
		{
			FreeCurCellContent ();
			retCode = ARMLOCAL_BumpCurve (C_inCvId, C_matu, C_rate, C_meth, C_result);
		
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_ARM_BUMPCURVE" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


__declspec(dllexport) LPXLOPER WINAPI Local_PXL_ARM_BUMPCURVE (LPXLOPER XL_inCv,
															   LPXLOPER XL_epsilon,
															   LPXLOPER XL_meth)
{
	ADD_LOG("Local_PXL_ARM_BUMPCURVE ");
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();

	// C variable
	CCString C_cvString;
	long C_inCvId;
	double C_meth;
	double C_meth_default = 0;

	VECTOR<CCString> C_epsilon;
	double C_epsilon_dbl;
	long epsilonType;

	VECTOR<CCString> C_matu;
	VECTOR<double> C_rate; 

	// error
	static int error;
	static char* reason = "";

	XL_readStrCell(XL_inCv,C_cvString," ARM_ERR: curve id: object expected",C_result);
	XL_readNumOrStrVector(XL_epsilon,C_epsilon_dbl,C_epsilon,epsilonType," ARM_ERR: epsilon: array of string expected",C_result);
	XL_readNumCellWD(XL_meth,C_meth,C_meth_default," ARM_ERR: method: numeric expected (0 - 'Classic'; 1 - 'Forward')",C_result);

	if (epsilonType == XL_TYPE_DOUBLE)
	{
		C_rate.push_back (C_epsilon_dbl);
	}
	else
	{
		double tmp_rate = -1;
		CCString tmp_matu ("NULL");

		if(C_epsilon.size () % 2 != 0)
		{
			C_result.setMsg ("ARM_ERR: check your maturities and rates array");
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}

		int real_size = C_epsilon.size () / 2;
		double tmp;
		int j = 0;

		for(int i = 0; i < real_size; i++)
		{
			C_epsilon[j].toUpper();
			C_matu.push_back (C_epsilon[j]);
			if((sscanf ((const char*)C_epsilon[j+1], "%lf", &tmp) != 1))
			{
				C_result.setMsg ("ARM_ERR: check your maturities and rates array");
				ARM_ARG_ERR();
				return (LPXLOPER)&XL_result;
			}
			C_rate.push_back (tmp);
			j += 2;
		}

		CCString last = C_matu[C_matu.size () - 1];
		if((!last) || (last[0] == ' '))
		{
			C_result.setMsg ("ARM_ERR: check your maturities and rates array");
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}
	}

	C_inCvId =  LocalGetNumObjectId (C_cvString);

	long retCode;
	long objId;

	CCString curClass = LocalGetStringObjectClass(C_cvString);
	CCString stringId;

	retCode = ARMLOCAL_BumpCurve (C_inCvId, C_matu, C_rate, C_meth, C_result);

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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_PXL_ARM_BUMPCURVE" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


__declspec(dllexport) LPXLOPER WINAPI Local_ARM_BumpSpreadedCurve(	LPXLOPER XL_inCv,
																	LPXLOPER XL_epsilon,
																	LPXLOPER XL_curveToBump,
																	LPXLOPER XL_meth)
{
	ADD_LOG("Local_ARM_BumpSpreadedCurve");
	static XLOPER XL_result;
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
		ARM_NOCALCIFWIZ();

		// C variable
		CCString C_cvString;
		long C_inCvId;
		double C_meth;
		double C_meth_default = 0;

		VECTOR<CCString> C_epsilon;
		double C_epsilon_dbl;
		long epsilonType;

		CCString	C_CurveToBump;

		VECTOR<CCString> C_matu;
		VECTOR<double> C_rate; 

		// error
		static int error;
		static char* reason = "";

		XL_readStrCell(XL_inCv,C_cvString," ARM_ERR: curve id: object expected",C_result);
		XL_readNumOrStrVector(XL_epsilon,C_epsilon_dbl,C_epsilon,epsilonType," ARM_ERR: epsilon: array of string expected",C_result);
		XL_readStrCellWD(XL_curveToBump,C_CurveToBump,"ZC"," ARM_ERR: curve to bump: string expected",C_result);
		XL_readNumCellWD(XL_meth,C_meth,C_meth_default," ARM_ERR: method: numeric expected (0 - 'Classic'; 1 - 'Forward')",C_result);

		if (epsilonType == XL_TYPE_DOUBLE)
		{
			C_rate.push_back (C_epsilon_dbl);
		}
		else
		{
			double tmp_rate = -1;
			CCString tmp_matu ("NULL");

			if(C_epsilon.size () % 2 != 0)
			{
				C_result.setMsg ("ARM_ERR: check your maturities and rates array");
				ARM_ARG_ERR();
				return (LPXLOPER)&XL_result;
			}

			int real_size = C_epsilon.size () / 2;
			double tmp;
			int j = 0;

			for(int i = 0; i < real_size; i++)
			{
				C_epsilon[j].toUpper();
				C_matu.push_back (C_epsilon[j]);
				if((sscanf ((const char*)C_epsilon[j+1], "%lf", &tmp) != 1))
				{
					C_result.setMsg ("ARM_ERR: check your maturities and rates array");
					ARM_ARG_ERR();
					return (LPXLOPER)&XL_result;
				}
				C_rate.push_back (tmp);
				j += 2;
			}

			CCString last = C_matu[C_matu.size () - 1];
			if((!last) || (last[0] == ' '))
			{
				C_result.setMsg ("ARM_ERR: check your maturities and rates array");
				ARM_ARG_ERR();
				return (LPXLOPER)&XL_result;
			}
		}

		C_inCvId =  LocalGetNumObjectId (C_cvString);

		long retCode;
		long objId;
		CCString prevClass;

		CCString curClass = LocalGetStringObjectClass(C_cvString);
		CCString stringId = GetLastCurCellEnvValue ();

		if(!stringId)
		{
			retCode = ARMLOCAL_BumpSpreadedCurve(C_inCvId, C_matu, C_rate, C_CurveToBump, C_meth, C_result);

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
				retCode = ARMLOCAL_BumpSpreadedCurve (C_inCvId, C_matu, C_rate, C_CurveToBump, C_meth, C_result, objId);

				if(retCode == ARM_OK)
				{
					LocalSetCurCellEnvValue (curClass, objId); 

					stringId = LocalMakeObjectId (objId, curClass);
				}
			}
			else
			{
				FreeCurCellContent ();
				retCode = ARMLOCAL_BumpSpreadedCurve (C_inCvId, C_matu, C_rate, C_CurveToBump, C_meth, C_result);
			
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
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_ARM_BumpSpreadedCurve" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


__declspec(dllexport) LPXLOPER WINAPI Local_PXL_ARM_BumpSpreadedCurve(	LPXLOPER XL_inCv,
																		LPXLOPER XL_epsilon,
																		LPXLOPER XL_curveToBump,
																		LPXLOPER XL_meth)
{
	ADD_LOG("Local_PXL_ARM_BumpSpreadedCurve");
	static XLOPER XL_result;
	ARM_result C_result;

	ARM_XL_TRY_BLOCK_BEGIN
	{
		ARM_NOCALCIFWIZ();

		CCString C_cvString;
		long C_inCvId;
		double C_meth;
		double C_meth_default = 0;

		VECTOR<CCString> C_epsilon;
		double C_epsilon_dbl;
		long epsilonType;

		CCString	C_CurveToBump;

		VECTOR<CCString> C_matu;
		VECTOR<double> C_rate; 

		// error
		static int error;
		static char* reason = "";

		XL_readStrCell(XL_inCv,C_cvString," ARM_ERR: curve id: object expected",C_result);
		XL_readNumOrStrVector(XL_epsilon,C_epsilon_dbl,C_epsilon,epsilonType," ARM_ERR: epsilon: array of string expected",C_result);
		XL_readStrCellWD(XL_curveToBump,C_CurveToBump,"ZC"," ARM_ERR: curve to bump: string expected",C_result);
		XL_readNumCellWD(XL_meth,C_meth,C_meth_default," ARM_ERR: method: numeric expected (0 - 'Classic'; 1 - 'Forward')",C_result);

		if (epsilonType == XL_TYPE_DOUBLE)
		{
			C_rate.push_back (C_epsilon_dbl);
		}
		else
		{
			double tmp_rate = -1;
			CCString tmp_matu ("NULL");

			if(C_epsilon.size () % 2 != 0)
			{
				C_result.setMsg ("ARM_ERR: check your maturities and rates array");
				ARM_ARG_ERR();
				return (LPXLOPER)&XL_result;
			}

			int real_size = C_epsilon.size () / 2;
			double tmp;
			int j = 0;

			for(int i = 0; i < real_size; i++)
			{
				C_epsilon[j].toUpper();
				C_matu.push_back (C_epsilon[j]);
				if((sscanf ((const char*)C_epsilon[j+1], "%lf", &tmp) != 1))
				{
					C_result.setMsg ("ARM_ERR: check your maturities and rates array");
					ARM_ARG_ERR();
					return (LPXLOPER)&XL_result;
				}
				C_rate.push_back (tmp);
				j += 2;
			}

			CCString last = C_matu[C_matu.size () - 1];
			if((!last) || (last[0] == ' '))
			{
				C_result.setMsg ("ARM_ERR: check your maturities and rates array");
				ARM_ARG_ERR();
				return (LPXLOPER)&XL_result;
			}
		}

		C_inCvId =  LocalGetNumObjectId (C_cvString);

		long retCode;
		long objId;

		CCString curClass = LocalGetStringObjectClass(C_cvString);
		CCString stringId;

		retCode = ARMLOCAL_BumpSpreadedCurve(C_inCvId, C_matu, C_rate, C_CurveToBump, C_meth, C_result);

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

	}
	ARM_XL_TRY_BLOCK_END

	ARM_XL_CATCH_ARM_EXPT

	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_PXL_ARM_BumpSpreadedCurve" )

	return (LPXLOPER)&XL_result;
}



__declspec(dllexport) LPXLOPER WINAPI Local_CreateZeroCurveLin (LPXLOPER XL_matuRate,
																LPXLOPER XL_meth,
																LPXLOPER XL_date,
																LPXLOPER XL_ccy,
																LPXLOPER XL_interpMeth)
{
	ADD_LOG("Local_CreateZeroCurveLin ");
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

    /// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	    ARM_NOCALCIFWIZ();

	    // C variable
	    VECTOR<CCString> C_matuRate;
	    VECTOR<double> C_matu;
	    VECTOR<double> C_rate;

	    CCString C_ccy;
	    CCString C_strMeth;
	    double C_meth;
	    long typeMeth;
	    
        long isDouble;

	    double C_date;

	    CCString C_interpMeth;
	    long interpMethId;

	    // error
	    static int error;
	    static char* reason = "";

	    XL_readStrVector(XL_matuRate,C_matuRate," ARM_ERR: maturities and rates: array of numeric expected",DOUBLE_TYPE,C_result);	
	    XL_readStrOrNumCell(XL_meth,C_strMeth,C_meth,typeMeth," ARM_ERR: decompounding method: string or numeric expected",C_result);
	    XL_readNumCell(XL_date,C_date," ARM_ERR: as of date: date expected",C_result);
	    XL_readStrCellWD(XL_ccy,C_ccy,"DEFAULT"," ARM_ERR: currency: string expected",C_result);
	    XL_readStrCellWD(XL_interpMeth,C_interpMeth,"LINEAR"," ARM_ERR: interpolation method: string expected",C_result);

	    double tmp_rate = -1;

	    if(C_matuRate.size () % 2 != 0)
	    {
		    C_result.setMsg ("ARM_ERR: check your maturities and rates array");
		    ARM_ARG_ERR();
		    return (LPXLOPER)&XL_result;
	    }

	    int real_size = C_matuRate.size () / 2;
	    int j = 0;

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

	    if ( typeMeth == XL_TYPE_STRING )
	    {
		    if((C_meth = ARM_ConvForwardYieldMethod (C_strMeth, C_result)) == ARM_DEFAULT_ERR)
		    {
			    ARM_ARG_ERR();
			    return (LPXLOPER)&XL_result;
		    }
	    }

	    if((interpMethId = ARM_ConvInterpMethod (C_interpMeth, C_result)) == ARM_DEFAULT_ERR)
	    {
		    ARM_ARG_ERR();
		    return (LPXLOPER)&XL_result;
	    }

	    long retCode;
	    double dSettleDate;

	    retCode = ARMLOCAL_GetSpotDays(C_ccy,C_result);

	    if (retCode == ARM_OK)
	    {
		    retCode = ARMLOCAL_NextBusinessDay(C_date,C_ccy,C_result.getDouble(),C_result);
		    dSettleDate = Local_ARMDATE2XLDATE(C_result.getString());
	    }
	    else
	    {
		    ARM_ERR();
		    return (LPXLOPER)&XL_result;
	    }

         ARM_CRV_TERMS sTerms;

	    for (int i = 0; i < real_size; i++)
	    {
		    isDouble = 0;

		    int Nb;
		    char cMatu;
		    long freqId;


		    sscanf(C_matuRate[j], "%d%c", &Nb, &cMatu);

            strcpy(sTerms[i], C_matuRate[j]);

		    cMatu = toupper(cMatu);

		    if ( cMatu == 'D' ) // Ex : "1D"
			    freqId = K_DAILY;
		    else if ( cMatu == 'W' )  
			    freqId = K_WEEKLY;
		    else if ( cMatu == 'M' ) 
			    freqId = K_MONTHLY;
		    else if ( cMatu == 'Y')  // ->implicitement ce sont des taux de swap
			    freqId = K_ANNUAL;
		    else
		    {
			    isDouble = 1;
			    C_matu.push_back (atof(C_matuRate[j]));
		    }
		    
		    if (isDouble == 0)
		    {
			    if (freqId == K_DAILY)
				    retCode = ARMLOCAL_NextBusinessDay(C_date,C_ccy,Nb,C_result);
			    else
				    retCode = ARMLOCAL_ARM_ADDPERIOD(dSettleDate, freqId, C_ccy, (long) Nb, K_FOLLOWING, 0, C_result);

			    if (retCode != ARM_OK)
			    {
				    ARM_ERR();
				    return (LPXLOPER)&XL_result;
			    }

			    double myDate = Local_ARMDATE2XLDATE(C_result.getString());

			    C_matu.push_back ((myDate - C_date) /365.);
		    }

		    C_rate.push_back(atof(C_matuRate[j+1]));

		    j += 2;
	    }

        int           matuAreDoubles = isDouble; // 1: YES, O: NO
       

	    long objId;
	    CCString prevClass;
	    
	    CCString curClass = LOCAL_ZERO_CURVE_LIN_CLASS;
	    CCString stringId = GetLastCurCellEnvValue ();
	    
	    if (!stringId)
	    {
		    retCode = ARMLOCAL_zclint(C_matu, C_rate, (long)C_meth, C_date, C_ccy, interpMethId,
                                      matuAreDoubles,
                                      sTerms,
                                      C_result);

		    if ( retCode == ARM_OK )
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
			    
		    if ( curClass == prevClass )
		    {
			    retCode = ARMLOCAL_zclint(C_matu, C_rate, (long) C_meth, C_date, C_ccy, interpMethId,
                                          matuAreDoubles,
                                          sTerms,
                                          C_result, objId);

			    if ( retCode == ARM_OK )
			    {
				    LocalSetCurCellEnvValue (curClass, objId); 

				    stringId = LocalMakeObjectId (objId, curClass);
			    }
		    }
		    else
		    {
			    FreeCurCellContent ();
			    
                retCode = ARMLOCAL_zclint(C_matu, C_rate, (long)C_meth, C_date, C_ccy, interpMethId,
                                          matuAreDoubles,
                                          sTerms,
                                          C_result);

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
		    FreeCurCellErr();

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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_CreateZeroCurveLin" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}



__declspec(dllexport) LPXLOPER WINAPI Local_PXL_CreateZeroCurveLin(LPXLOPER XL_matuRate,
																   LPXLOPER XL_meth,
																   LPXLOPER XL_date,
																   LPXLOPER XL_ccy,
																   LPXLOPER XL_interpMeth)
{
	ADD_LOG("Local_PXL_CreateZeroCurveLin");
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	    ARM_NOCALCIFWIZ();

	    // C variable
	    VECTOR<CCString> C_matuRate;
	    VECTOR<double> C_matu;
	    VECTOR<double> C_rate;

        long isDouble;

	    CCString C_ccy;
	    CCString C_strMeth;
	    double C_meth;
	    long typeMeth;
	    
	    double C_date;

	    CCString C_interpMeth;
	    long interpMethId;

	    // error
	    static int error;
	    static char* reason = "";

	    XL_readStrVector(XL_matuRate,C_matuRate," ARM_ERR: maturities and rates: array of numeric expected",DOUBLE_TYPE,C_result);	
	    XL_readStrOrNumCell(XL_meth,C_strMeth,C_meth,typeMeth," ARM_ERR: decompounding method: string or numeric expected",C_result);
	    XL_readNumCell(XL_date,C_date," ARM_ERR: as of date: date expected",C_result);
	    XL_readStrCellWD(XL_ccy,C_ccy,"DEFAULT"," ARM_ERR: currency: string expected",C_result);
	    XL_readStrCellWD(XL_interpMeth,C_interpMeth,"LINEAR"," ARM_ERR: interpolation method: string expected",C_result);

	    double tmp_rate = -1;

	    if(C_matuRate.size () % 2 != 0)
	    {
		    C_result.setMsg ("ARM_ERR: check your maturities and rates array");
		    ARM_ARG_ERR();
		    return (LPXLOPER)&XL_result;
	    }

	    int real_size = C_matuRate.size () / 2;
	    int j = 0;

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

	    if ( typeMeth == XL_TYPE_STRING )
	    {
		    if((C_meth = ARM_ConvForwardYieldMethod (C_strMeth, C_result)) == ARM_DEFAULT_ERR)
		    {
			    ARM_ARG_ERR();
			    return (LPXLOPER)&XL_result;
		    }
	    }

	    if((interpMethId = ARM_ConvInterpMethod (C_interpMeth, C_result)) == ARM_DEFAULT_ERR)
	    {
		    ARM_ARG_ERR();
		    return (LPXLOPER)&XL_result;
	    }

	    long retCode;
	    double dSettleDate;

	    retCode = ARMLOCAL_GetSpotDays(C_ccy,C_result);

	    if (retCode == ARM_OK)
	    {
		    retCode = ARMLOCAL_NextBusinessDay(C_date,C_ccy,C_result.getDouble(),C_result);
		    dSettleDate = Local_ARMDATE2XLDATE(C_result.getString());
	    }
	    else
	    {
		    ARM_ERR();
		    return (LPXLOPER)&XL_result;
	    }

        ARM_CRV_TERMS sTerms;

	    for (int i = 0; i < real_size; i++)
	    {
		    isDouble = 0;

		    int Nb;
		    char cMatu;
		    long freqId;

		    sscanf(C_matuRate[j], "%d%c", &Nb, &cMatu);

            strcpy(sTerms[i], C_matuRate[j]);

		    cMatu = toupper(cMatu);

		    if ( cMatu == 'D' ) // Ex : "1D"
			    freqId = K_DAILY;
		    else if ( cMatu == 'W' )  
			    freqId = K_WEEKLY;
		    else if ( cMatu == 'M' ) 
			    freqId = K_MONTHLY;
		    else if ( cMatu == 'Y')  // ->implicitement ce sont des taux de swap
			    freqId = K_ANNUAL;
		    else
		    {
			    isDouble = 1;
			    C_matu.push_back (atof(C_matuRate[j]));
		    }
		    
		    if ( isDouble == 0 )
		    {
			    if (freqId == K_DAILY)
				    retCode = ARMLOCAL_NextBusinessDay(C_date,C_ccy,Nb,C_result);
			    else
				    retCode = ARMLOCAL_ARM_ADDPERIOD(dSettleDate, freqId, C_ccy, (long) Nb, K_FOLLOWING, 0, C_result);

			    if (retCode != ARM_OK)
			    {
				    ARM_ERR();
				    return (LPXLOPER)&XL_result;
			    }

			    double myDate = Local_ARMDATE2XLDATE(C_result.getString());

			    C_matu.push_back ((myDate - C_date) /365.);
		    }

		    C_rate.push_back (atof(C_matuRate[j+1]));

		    j += 2;
	    }

        int           matuAreDoubles = isDouble; // 1: YES, O: NO
        
	    
        long objId;
	    
	    CCString curClass = LOCAL_ZERO_CURVE_LIN_CLASS;
	    CCString stringId;
	    
	    retCode = ARMLOCAL_zclint(C_matu, C_rate, (long) C_meth, C_date, C_ccy, interpMethId,
                                  matuAreDoubles,
                                  sTerms,
                                  C_result);

	    if ( retCode == ARM_OK )
	    {
		    objId = C_result.getLong ();

		    stringId = LocalMakeObjectId (objId, curClass);
	    }
	    
	    if ( retCode == ARM_OK )
	    {
		    FreeCurCellErr();

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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_PXL_CreateZeroCurveLin" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}





__declspec(dllexport) LPXLOPER WINAPI Local_ShiftZeroCurveLin (LPXLOPER XL_value,
                                                               LPXLOPER XL_nbPlot,
                                                               LPXLOPER XL_matuRate,
																LPXLOPER XL_meth,
																LPXLOPER XL_date,
																LPXLOPER XL_ccy,
																LPXLOPER XL_interpMeth)
{
	ADD_LOG("Local_ShiftZeroCurveLin ");
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

    /// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();

	// C variable
	VECTOR<CCString> C_matuRate;
	VECTOR<double> C_matu;
	VECTOR<double> C_rate;

	CCString C_ccy;
	CCString C_strMeth;
	double C_meth;
	long typeMeth;
	
	double C_date;

	CCString C_interpMeth;
	long interpMethId;

    double C_value;
    double C_nbPlot;

	// error
	static int error;
	static char* reason = "";
    XL_readNumCell(XL_value,C_value," ARM_ERR: value: numeric expected",C_result);
    XL_readNumCell(XL_nbPlot,C_nbPlot," ARM_ERR: nbPlot: numeric expected",C_result);

	XL_readStrVector(XL_matuRate,C_matuRate," ARM_ERR: maturities and rates: array of numeric expected",DOUBLE_TYPE,C_result);	
	XL_readStrOrNumCell(XL_meth,C_strMeth,C_meth,typeMeth," ARM_ERR: decompounding method: string or numeric expected",C_result);
	XL_readNumCell(XL_date,C_date," ARM_ERR: as of date: date expected",C_result);
	XL_readStrCellWD(XL_ccy,C_ccy,"DEFAULT"," ARM_ERR: currency: string expected",C_result);
	XL_readStrCellWD(XL_interpMeth,C_interpMeth,"LINEAR"," ARM_ERR: interpolation method: string expected",C_result);

	double tmp_rate = -1;

	if(C_matuRate.size () % 2 != 0)
	{
		C_result.setMsg ("ARM_ERR: check your maturities and rates array");
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	int real_size = C_matuRate.size () / 2;
	int j = 0;

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

	if ( typeMeth == XL_TYPE_STRING )
	{
		if((C_meth = ARM_ConvForwardYieldMethod (C_strMeth, C_result)) == ARM_DEFAULT_ERR)
		{
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}
	}

	if((interpMethId = ARM_ConvInterpMethod (C_interpMeth, C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	long retCode;
	double dSettleDate;

	retCode = ARMLOCAL_GetSpotDays(C_ccy,C_result);

	if (retCode == ARM_OK)
	{
		retCode = ARMLOCAL_NextBusinessDay(C_date,C_ccy,C_result.getDouble(),C_result);
		dSettleDate = Local_ARMDATE2XLDATE(C_result.getString());
	}
	else
	{
		ARM_ERR();
		return (LPXLOPER)&XL_result;
	}

	for (int i=0;i<real_size;i++)
	{
		long isDouble = 0;

		int Nb;
		char cMatu;
		long freqId;

		sscanf(C_matuRate[j], "%d%c", &Nb, &cMatu);

		cMatu = toupper(cMatu);

		if ( cMatu == 'D' ) // Ex : "1D"
			freqId = K_DAILY;
		else if ( cMatu == 'W' )  
			freqId = K_WEEKLY;
		else if ( cMatu == 'M' ) 
			freqId = K_MONTHLY;
		else if ( cMatu == 'Y')  // ->implicitement ce sont des taux de swap
			freqId = K_ANNUAL;
		else
		{
			isDouble = 1;
			C_matu.push_back (atof(C_matuRate[j]));
		}
		
		if (isDouble == 0)
		{
			if (freqId == K_DAILY)
				retCode = ARMLOCAL_NextBusinessDay(C_date,C_ccy,Nb,C_result);
			else
				retCode = ARMLOCAL_ARM_ADDPERIOD(dSettleDate, freqId, C_ccy, (long) Nb, K_FOLLOWING, 0, C_result);

			if (retCode != ARM_OK)
			{
				ARM_ERR();
				return (LPXLOPER)&XL_result;
			}

			double myDate = Local_ARMDATE2XLDATE(C_result.getString());

			C_matu.push_back ((myDate - C_date) /365.);
		}

		C_rate.push_back (atof(C_matuRate[j+1]));
		j += 2;
	}

	long objId;
	CCString prevClass;
	
	CCString curClass = LOCAL_ZERO_CURVE_LIN_CLASS;
	CCString stringId = GetLastCurCellEnvValue ();
	
	if(!stringId)
	{
		retCode = ARMLOCAL_shiftzclint (C_value,C_nbPlot,C_matu, C_rate, (long)C_meth, C_date, C_ccy, interpMethId, C_result);

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
			retCode = ARMLOCAL_shiftzclint (C_value,C_nbPlot,C_matu, C_rate, (long)C_meth, C_date, C_ccy, interpMethId, C_result, objId);

			if(retCode == ARM_OK)
			{
				LocalSetCurCellEnvValue (curClass, objId); 

				stringId = LocalMakeObjectId (objId, curClass);
			}
		}
		else
		{
			FreeCurCellContent ();
			retCode = ARMLOCAL_shiftzclint (C_value,C_nbPlot,C_matu, C_rate, (long)C_meth, C_date, C_ccy, interpMethId, C_result);

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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_CreateZeroCurveLin" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}

__declspec(dllexport) LPXLOPER WINAPI Local_PXL_ShiftZeroCurveLin (LPXLOPER XL_value,
                                                                    LPXLOPER XL_nbPlot,
                                                                    LPXLOPER XL_matuRate,
																	LPXLOPER XL_meth,
																	LPXLOPER XL_date,
																	LPXLOPER XL_ccy,
																	LPXLOPER XL_interpMeth)
{
	ADD_LOG("Local_PXL_ShiftZeroCurveLin ");
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();

	// C variable
	VECTOR<CCString> C_matuRate;
	VECTOR<double> C_matu;
	VECTOR<double> C_rate;

	CCString C_ccy;
	CCString C_strMeth;
	double C_meth;
	long typeMeth;
	
	double C_date;

	CCString C_interpMeth;
	long interpMethId;

    double C_value;
    double C_nbPlot;

	// error
	static int error;
	static char* reason = "";
    XL_readNumCell(XL_value,C_value," ARM_ERR: value: numeric expected",C_result);
    XL_readNumCell(XL_nbPlot,C_nbPlot," ARM_ERR: nbPlot: numeric expected",C_result);

	XL_readStrVector(XL_matuRate,C_matuRate," ARM_ERR: maturities and rates: array of numeric expected",DOUBLE_TYPE,C_result);	
	XL_readStrOrNumCell(XL_meth,C_strMeth,C_meth,typeMeth," ARM_ERR: decompounding method: string or numeric expected",C_result);
	XL_readNumCell(XL_date,C_date," ARM_ERR: as of date: date expected",C_result);
	XL_readStrCellWD(XL_ccy,C_ccy,"DEFAULT"," ARM_ERR: currency: string expected",C_result);
	XL_readStrCellWD(XL_interpMeth,C_interpMeth,"LINEAR"," ARM_ERR: interpolation method: string expected",C_result);

	double tmp_rate = -1;

	if(C_matuRate.size () % 2 != 0)
	{
		C_result.setMsg ("ARM_ERR: check your maturities and rates array");
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	int real_size = C_matuRate.size () / 2;
	int j = 0;

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

	if ( typeMeth == XL_TYPE_STRING )
	{
		if((C_meth = ARM_ConvForwardYieldMethod (C_strMeth, C_result)) == ARM_DEFAULT_ERR)
		{
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}
	}

	if((interpMethId = ARM_ConvInterpMethod (C_interpMeth, C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	long retCode;
	double dSettleDate;

	retCode = ARMLOCAL_GetSpotDays(C_ccy,C_result);

	if (retCode == ARM_OK)
	{
		retCode = ARMLOCAL_NextBusinessDay(C_date,C_ccy,C_result.getDouble(),C_result);
		dSettleDate = Local_ARMDATE2XLDATE(C_result.getString());
	}
	else
	{
		ARM_ERR();
		return (LPXLOPER)&XL_result;
	}

	for (int i=0;i<real_size;i++)
	{
		long isDouble = 0;

		int Nb;
		char cMatu;
		long freqId;

		sscanf(C_matuRate[j], "%d%c", &Nb, &cMatu);

		cMatu = toupper(cMatu);

		if ( cMatu == 'D' ) // Ex : "1D"
			freqId = K_DAILY;
		else if ( cMatu == 'W' )  
			freqId = K_WEEKLY;
		else if ( cMatu == 'M' ) 
			freqId = K_MONTHLY;
		else if ( cMatu == 'Y')  // ->implicitement ce sont des taux de swap
			freqId = K_ANNUAL;
		else
		{
			isDouble = 1;
			C_matu.push_back (atof(C_matuRate[j]));
		}
		
		if (isDouble == 0)
		{
			if (freqId == K_DAILY)
				retCode = ARMLOCAL_NextBusinessDay(C_date,C_ccy,Nb,C_result);
			else
				retCode = ARMLOCAL_ARM_ADDPERIOD(dSettleDate, freqId, C_ccy, (long) Nb, K_FOLLOWING, 0, C_result);

			if (retCode != ARM_OK)
			{
				ARM_ERR();
				return (LPXLOPER)&XL_result;
			}

			double myDate = Local_ARMDATE2XLDATE(C_result.getString());

			C_matu.push_back ((myDate - C_date) /365.);
		}

		C_rate.push_back (atof(C_matuRate[j+1]));
		j += 2;
	}

	long objId;
	
	CCString curClass = LOCAL_ZERO_CURVE_LIN_CLASS;
	CCString stringId;
	
	retCode = ARMLOCAL_shiftzclint (C_value,C_nbPlot,C_matu, C_rate, (long)C_meth, C_date, C_ccy, interpMethId, C_result);

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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_PXL_CreateZeroCurveLin" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}



















__declspec(dllexport) LPXLOPER WINAPI Local_CreateTOYNYZCSwapInt (LPXLOPER XL_date,
																  LPXLOPER XL_matuRate,
																  LPXLOPER XL_mmVsFut,
																  LPXLOPER XL_swapVsFut,
																  LPXLOPER XL_raw,
																  LPXLOPER XL_interp,
																  LPXLOPER XL_ccy,
																  LPXLOPER XL_frq)
{
	ADD_LOG("Local_CreateTOYNYZCSwapInt ");
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();

	// C variable
	double C_date;
	VECTOR<CCString> C_matuRate;
	VECTOR<CCString> C_matu;
	VECTOR<double> C_rate; 

	CCString C_mmVsFut;
	long mmVsFutId;

	CCString C_swapVsFut;
	long swapVsFutId;

	CCString C_raw;
	long rawId;

	CCString C_interp;
	long interpId;

	CCString C_ccy;

	double C_frq;
	double C_frq_default = 0.;
	
	// error
	static int error;
	static char* reason = "";

	XL_readNumCell(XL_date,C_date," ARM_ERR: as of date: date expected",C_result);
	
	XL_readStrVector(XL_matuRate,C_matuRate," ARM_ERR: maturities and rates: array of numeric expected",DOUBLE_TYPE,C_result);

	double tmp_rate = -1;
	CCString tmp_matu ("NULL");

	if(C_matuRate.size () % 2 != 0)
	{
		C_result.setMsg ("ARM_ERR: check your maturities and rates array");
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	int real_size = C_matuRate.size () / 2;
	double tmp;
	int j = 0;

	for(int i = 0; i < real_size; i++)
	{
		C_matu.push_back (C_matuRate[j]);
		if((sscanf ((const char*)C_matuRate[j+1], "%lf", &tmp) != 1))
		{
			C_result.setMsg ("ARM_ERR: check your maturities and rates array");
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}
		C_rate.push_back (tmp);
		j += 2;
	}

	CCString last = C_matu[C_matu.size () - 1];
	if((!last) || (last[0] == ' '))
	{
		C_result.setMsg ("ARM_ERR: check your maturities and rates array");
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	XL_readStrCellWD(XL_mmVsFut,C_mmVsFut,"MM"," ARM_ERR: money market versus future: string expected",C_result);
	XL_readStrCellWD(XL_swapVsFut,C_swapVsFut,"SWAP"," ARM_ERR: swap versus future: string expected",C_result);
	XL_readStrCellWD(XL_raw,C_raw,"P"," ARM_ERR: raw: string expected",C_result);
	XL_readStrCellWD(XL_interp,C_interp,"C"," ARM_ERR: interpolation mode: string expected",C_result);
	XL_readStrCellWD(XL_ccy,C_ccy,"DEFAULT"," ARM_ERR: currency: string expected",C_result);
	XL_readNumCellWD(XL_frq,C_frq,C_frq_default," ARM_ERR: frequency: numeric expected",C_result);

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

	if((mmVsFutId = ARM_ConvMktType (C_mmVsFut, C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	if((swapVsFutId = ARM_ConvMktType (C_swapVsFut, C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	rawId = ARM_ConvCvMethod (C_raw);

	if((interpId = ARM_ConvInterpMethod (C_interp, C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	long retCode;
	long objId;
	CCString prevClass;
	
	CCString curClass = LOCAL_ZERO_CURVE_LIN_CLASS;
	CCString stringId = GetLastCurCellEnvValue ();
	
	if(!stringId)
	{
		retCode = ARMLOCAL_CreateTOYNYZCSwapInt(C_date, C_matu, C_rate, (long)mmVsFutId,
												(long)swapVsFutId, (long)rawId,
												(long)interpId, C_ccy, (long)C_frq,
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
			retCode = ARMLOCAL_CreateTOYNYZCSwapInt(C_date, C_matu, C_rate, (long)mmVsFutId,
													(long)swapVsFutId, (long)rawId,
													(long)interpId, C_ccy, (long)C_frq,
													C_result,objId);

			if(retCode == ARM_OK)
			{			
				LocalSetCurCellEnvValue (curClass, objId); 

				stringId = LocalMakeObjectId (objId, curClass);
			}
		}
		else
		{
			FreeCurCellContent ();
			retCode = ARMLOCAL_CreateTOYNYZCSwapInt(C_date, C_matu, C_rate, (long)mmVsFutId,
													(long)swapVsFutId, (long)rawId,
													(long)interpId, C_ccy, (long)C_frq,
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_CreateTOYNYZCSwapInt" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}

__declspec(dllexport) LPXLOPER WINAPI Local_PXL_CreateTOYNYZCSwapInt (LPXLOPER XL_date,
																	  LPXLOPER XL_matuRate,
																	  LPXLOPER XL_mmVsFut,
																	  LPXLOPER XL_swapVsFut,
																	  LPXLOPER XL_raw,
																	  LPXLOPER XL_interp,
																	  LPXLOPER XL_ccy,
																	  LPXLOPER XL_frq)
{
	ADD_LOG("Local_PXL_CreateTOYNYZCSwapInt ");
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();

	// C variable
	double C_date;
	VECTOR<CCString> C_matuRate;
	VECTOR<CCString> C_matu;
	VECTOR<double> C_rate; 

	CCString C_mmVsFut;
	long mmVsFutId;

	CCString C_swapVsFut;
	long swapVsFutId;

	CCString C_raw;
	long rawId;

	CCString C_interp;
	long interpId;

	CCString C_ccy;

	double C_frq;
	double C_frq_default = 0;

	// error
	static int error;
	static char* reason = "";

	XL_readNumCell(XL_date,C_date," ARM_ERR: as of date: date expected",C_result);
	
	XL_readStrVector(XL_matuRate,C_matuRate," ARM_ERR: maturities and rates: array of numeric expected",DOUBLE_TYPE,C_result);

	double tmp_rate = -1;
	CCString tmp_matu ("NULL");

	if(C_matuRate.size () % 2 != 0)
	{
		C_result.setMsg ("ARM_ERR: check your maturities and rates array");
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	int real_size = C_matuRate.size () / 2;
	double tmp;
	int j = 0;

	for(int i = 0; i < real_size; i++)
	{
		C_matu.push_back (C_matuRate[j]);
		if((sscanf ((const char*)C_matuRate[j+1], "%lf", &tmp) != 1))
		{
			C_result.setMsg ("ARM_ERR: check your maturities and rates array");
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}
		C_rate.push_back (tmp);
		j += 2;
	}

	CCString last = C_matu[C_matu.size () - 1];
	if((!last) || (last[0] == ' '))
	{
		C_result.setMsg ("ARM_ERR: check your maturities and rates array");
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	XL_readStrCellWD(XL_mmVsFut,C_mmVsFut,"MM"," ARM_ERR: money market versus future: string expected",C_result);
	XL_readStrCellWD(XL_swapVsFut,C_swapVsFut,"SWAP"," ARM_ERR: swap versus future: string expected",C_result);
	XL_readStrCellWD(XL_raw,C_raw,"P"," ARM_ERR: raw: string expected",C_result);
	XL_readStrCellWD(XL_interp,C_interp,"C"," ARM_ERR: interpolation mode: string expected",C_result);
	XL_readStrCellWD(XL_ccy,C_ccy,"USD"," ARM_ERR: currency: string expected",C_result);
	XL_readNumCellWD(XL_frq,C_frq,C_frq_default," ARM_ERR: frequency: numeric expected",C_result);

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
	
	if((mmVsFutId = ARM_ConvMktType (C_mmVsFut, C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	if((swapVsFutId = ARM_ConvMktType (C_swapVsFut, C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	rawId = ARM_ConvCvMethod (C_raw);

	if((interpId = ARM_ConvInterpMethod (C_interp, C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	long retCode;
	long objId;
	
	CCString curClass = LOCAL_ZERO_CURVE_LIN_CLASS;
	CCString stringId;

	retCode = ARMLOCAL_CreateTOYNYZCSwapInt (C_date, C_matu, C_rate, (long)mmVsFutId,
											(long)swapVsFutId, (long)rawId,
											(long)interpId, C_ccy, (long)C_frq, C_result);

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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_PXL_CreateTOYNYZCSwapInt" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}



__declspec(dllexport) LPXLOPER WINAPI Local_CreateZCCashInt (LPXLOPER XL_date,
															 LPXLOPER XL_matuRate,
															 LPXLOPER XL_bonds,
															 LPXLOPER XL_mmVsFut,
															 LPXLOPER XL_ccy)
{
	ADD_LOG("Local_CreateZCCashInt ");
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

    /// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();

	// C variable
	double C_date;

	VECTOR<CCString> C_matuRate;
	VECTOR<CCString> C_matu;
	VECTOR<double> C_rate; 

	VECTOR<CCString> C_bonds;
	VECTOR<long> C_bondsId;
	VECTOR<double> C_yields; 

	CCString C_mmVsFut;
	long mmVsFutId;

	CCString C_ccy;
	
	// error
	static int error;
	static char* reason = "";

	XL_readNumCell(XL_date,C_date," ARM_ERR: as of date: date expected",C_result);
	
	XL_readStrVector(XL_matuRate,C_matuRate," ARM_ERR: maturities and rates: array of numeric expected",DOUBLE_TYPE,C_result);

	double tmp_rate = -1;
	CCString tmp_matu ("NULL");

	if(C_matuRate.size () % 2 != 0)
	{
		C_result.setMsg ("ARM_ERR: check your maturities and rates array");
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	int real_size = C_matuRate.size () / 2;
	double tmp;
	int j = 0;

	for(int i = 0; i < real_size; i++)
	{
		C_matu.push_back (C_matuRate[j]);
		if((sscanf ((const char*)C_matuRate[j+1], "%lf", &tmp) != 1))
		{
			C_result.setMsg ("ARM_ERR: check your maturities and rates array");
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}
		C_rate.push_back (tmp);
		j += 2;
	}

	CCString last = C_matu[C_matu.size () - 1];
	if((!last) || (last[0] == ' '))
	{
		C_result.setMsg ("ARM_ERR: check your maturities and rates array");
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	XL_readStrVector(XL_bonds,C_bonds," ARM_ERR: bonds: array expected",DOUBLE_TYPE,C_result);

	if(C_bonds.size () % 2 != 0)
	{
		C_result.setMsg ("ARM_ERR: check your bonds and yields array");
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	real_size = C_bonds.size () / 2;

	j = 0;

	for(i = 0; i < real_size; i++)
	{
		C_bondsId.push_back (LocalGetNumObjectId (C_bonds[j]));
		if((sscanf ((const char*)C_bonds[j+1], "%lf", &tmp) != 1))
		{
			C_result.setMsg ("ARM_ERR: check your bonds and yields array");
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}
		C_yields.push_back (tmp);
		j += 2;
	}

	last = C_bonds[C_bonds.size () - 1];
	if((!last) || (last[0] == ' '))
	{
		C_result.setMsg ("ARM_ERR: check your bonds and yields array");
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	XL_readStrCell(XL_mmVsFut,C_mmVsFut," ARM_ERR: money market versus future: string expected",C_result);
	XL_readStrCellWD(XL_ccy,C_ccy,"DEFAULT"," ARM_ERR: currency: string expected",C_result);

	if(C_ccy == "DEFAULT")
	{
		ARM_result ccyres;
		ARMLOCAL_ARM_GetDefaultCurrency (ccyres);
		if(ccyres.getRetCode () != ARM_OK)
		{
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}
		else
		{
			C_ccy = ccyres.getString ();
		}
	}
	
	if((mmVsFutId = ARM_ConvMktType (C_mmVsFut, C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	long retCode;
	long objId;
	CCString prevClass;
	
	CCString curClass = LOCAL_ZERO_CURVE_LIN_CLASS;
	CCString stringId = GetLastCurCellEnvValue ();
	
	if(!stringId)
	{
		retCode = ARMLOCAL_CreateZCCashInt (C_date,
											C_matu,
											C_rate, 
											C_bondsId,
											C_yields,
											(long)mmVsFutId,
											C_ccy,
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
			retCode = ARMLOCAL_CreateZCCashInt (C_date,
												C_matu,
												C_rate,
												C_bondsId,
												C_yields,
												(long)mmVsFutId,
												C_ccy,
												C_result,
												objId);
		}
		else
		{
			FreeCurCellContent ();
			retCode = ARMLOCAL_CreateZCCashInt (C_date,
												C_matu,
												C_rate,
												C_bondsId,
												C_yields,
												(long)mmVsFutId,
												C_ccy,
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_CreateZCCashInt" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}

__declspec(dllexport) LPXLOPER WINAPI Local_PXL_CreateZCCashInt (LPXLOPER XL_date,
																 LPXLOPER XL_matuRate,
																 LPXLOPER XL_bonds,
																 LPXLOPER XL_mmVsFut,
																 LPXLOPER XL_ccy)
{
	ADD_LOG("Local_PXL_CreateZCCashInt ");
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();

	// C variable
	double C_date;

	VECTOR<CCString> C_matuRate;
	VECTOR<CCString> C_matu;
	VECTOR<double> C_rate; 

	VECTOR<CCString> C_bonds;
	VECTOR<long> C_bondsId;
	VECTOR<double> C_yields; 

	CCString C_mmVsFut;
	long mmVsFutId;

	CCString C_ccy;
	
	// error
	static int error;
	static char* reason = "";

	XL_readNumCell(XL_date,C_date," ARM_ERR: as of date: date expected",C_result);
	
	XL_readStrVector(XL_matuRate,C_matuRate," ARM_ERR: maturities and rates: array of numeric expected",DOUBLE_TYPE,C_result);

	double tmp_rate = -1;
	CCString tmp_matu ("NULL");

	if(C_matuRate.size () % 2 != 0)
	{
		C_result.setMsg ("ARM_ERR: check your maturities and rates array");
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	int real_size = C_matuRate.size () / 2;
	double tmp;
	int j = 0;

	for(int i = 0; i < real_size; i++)
	{
		C_matu.push_back (C_matuRate[j]);
		if((sscanf ((const char*)C_matuRate[j+1], "%lf", &tmp) != 1))
		{
			C_result.setMsg ("ARM_ERR: check your maturities and rates array");
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}
		C_rate.push_back (tmp);
		j += 2;
	}

	CCString last = C_matu[C_matu.size () - 1];
	if((!last) || (last[0] == ' '))
	{
		C_result.setMsg ("ARM_ERR: check your maturities and rates array");
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	XL_readStrVector(XL_bonds,C_bonds," ARM_ERR: bonds: array expected",DOUBLE_TYPE,C_result);

	if(C_bonds.size () % 2 != 0)
	{
		C_result.setMsg ("ARM_ERR: check your bonds and yields array");
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	real_size = C_bonds.size () / 2;

	j = 0;

	for(i = 0; i < real_size; i++)
	{
		C_bondsId.push_back (LocalGetNumObjectId (C_bonds[j]));
		if((sscanf ((const char*)C_bonds[j+1], "%lf", &tmp) != 1))
		{
			C_result.setMsg ("ARM_ERR: check your bonds and yields array");
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}
		C_yields.push_back (tmp);
		j += 2;
	}

	last = C_bonds[C_bonds.size () - 1];
	if((!last) || (last[0] == ' '))
	{
		C_result.setMsg ("ARM_ERR: check your bonds and yields array");
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	XL_readStrCell(XL_mmVsFut,C_mmVsFut," ARM_ERR: money market versus future: string expected",C_result);
	XL_readStrCellWD(XL_ccy,C_ccy,"DEFAULT"," ARM_ERR: currency: string expected",C_result);

	if(C_ccy == "DEFAULT")
	{
		ARM_result ccyres;
		ARMLOCAL_ARM_GetDefaultCurrency (ccyres);
		if(ccyres.getRetCode () != ARM_OK)
		{
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}
		else
		{
			C_ccy = ccyres.getString ();
		}
	}
	
	if((mmVsFutId = ARM_ConvMktType (C_mmVsFut, C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	long retCode;
	long objId;
	
	CCString curClass = LOCAL_ZERO_CURVE_LIN_CLASS;
	CCString stringId;

	retCode = ARMLOCAL_CreateZCCashInt (C_date,
										C_matu,
										C_rate,
										C_bondsId,
										C_yields,
										(long)mmVsFutId,
										C_ccy,
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_PXL_CreateZCCashInt" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}



__declspec(dllexport) LPXLOPER WINAPI Local_GenerateBasisAdjCurve(LPXLOPER XL_Crv1,
																  LPXLOPER XL_BSCrv1,
																  LPXLOPER XL_Crv2,
																  LPXLOPER XL_BSCrv2,
																  LPXLOPER XL_RetSprds,
																  LPXLOPER XL_matuVec,
																  LPXLOPER XL_BSasSprds)
{
	ADD_LOG("Local_GenerateBasisAdjCurve");
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
		ARM_NOCALCIFWIZ();

		// C variable
		CCString C_Crv1;
		CCString C_BSCrv1;
		CCString C_Crv2;
		CCString C_BSCrv2;
		CCString C_BSasSprds;
		CCString C_RetSprds;

		long flagInputAsSprds;	
		long flagRetSprds;	
		VECTOR<CCString> C_matu;
		VECTOR<CCString> C_matuDef(0);
		// error
		static int error;
		static char* reason = "";

		XL_readStrCell(XL_Crv1,C_Crv1," ARM_ERR:  curve1 id: object expected",C_result);
		XL_readStrCell(XL_BSCrv1,C_BSCrv1," ARM_ERR: BS curve1 id: object expected",C_result);
		XL_readStrCell(XL_Crv2,C_Crv2," ARM_ERR: curve2 id: object expected",C_result);
		XL_readStrCell(XL_BSCrv2,C_BSCrv2," ARM_ERR: BS curve2 id: object expected",C_result);
		XL_readStrCellWD(XL_RetSprds,C_RetSprds,"N"," ARM_ERR: RetSprds: string expected",C_result);
		XL_readStrVectorWD(XL_matuVec,C_matu,C_matuDef," ARM_ERR: maturities: array of string expected",DOUBLE_TYPE,C_result);
		XL_readStrCellWD(XL_BSasSprds,C_BSasSprds,"N"," ARM_ERR: BSasSprds: string expected",C_result);
		
		if(C_BSasSprds == "Y")
		{
			flagInputAsSprds = 1;
		}
		else if(C_BSasSprds == "N")
		{
			flagInputAsSprds = 0;
		}
		else
		{
			throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
							   "String \"Y\" or \"N\" Expected  BS_AsSpreads");
		}

		if(C_RetSprds == "Y")
		{
			flagRetSprds = 1;
		}
		else if(C_RetSprds == "N")
		{
			flagRetSprds = 0;
		}
		else
		{
			throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
							   "String \"Y\" or \"N\" Expected  ReturnSprds");
		}

		long retCode;
		long objId;
		CCString prevClass;
		
		CCString curClass = LOCAL_ZERO_CURVE_LIN_CLASS;
		CCString stringId = GetLastCurCellEnvValue ();

		if(!stringId)
		{
			retCode = ARMLOCAL_GenerateBasisAdjCurve(LocalGetNumObjectId(C_Crv1),
													 LocalGetNumObjectId(C_BSCrv1),
													 LocalGetNumObjectId(C_Crv2),
													 LocalGetNumObjectId(C_BSCrv2),
													 flagInputAsSprds,
													 flagRetSprds,
													 C_matu,
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
				retCode = ARMLOCAL_GenerateBasisAdjCurve(LocalGetNumObjectId(C_Crv1),
														 LocalGetNumObjectId(C_BSCrv1),
														 LocalGetNumObjectId(C_Crv2),
														 LocalGetNumObjectId(C_BSCrv2),
														 flagInputAsSprds,
														 flagRetSprds,
														 C_matu,
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
				retCode = ARMLOCAL_GenerateBasisAdjCurve(LocalGetNumObjectId(C_Crv1),
														 LocalGetNumObjectId(C_BSCrv1),
														 LocalGetNumObjectId(C_Crv2),
														 LocalGetNumObjectId(C_BSCrv2),
														 flagInputAsSprds,
														 flagRetSprds,
														 C_matu,
														 C_result);
					
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure! in Local_GenerateBasisAdjCurve" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}

__declspec(dllexport) LPXLOPER WINAPI Local_PXL_GenerateBasisAdjCurve(LPXLOPER XL_Crv1,
																  LPXLOPER XL_BSCrv1,
																  LPXLOPER XL_Crv2,
																  LPXLOPER XL_BSCrv2,
																  LPXLOPER XL_RetSprds,
																  LPXLOPER XL_matuVec,
																  LPXLOPER XL_BSasSprds)
{
	ADD_LOG("Local_PXL_GenerateBasisAdjCurve");
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
		ARM_NOCALCIFWIZ();

		// C variable
		CCString C_Crv1;
		CCString C_BSCrv1;
		CCString C_Crv2;
		CCString C_BSasSprds;
		CCString C_BSCrv2;
		CCString C_RetSprds;

		long flagInputAsSprds;
		long flagRetSprds;	
		VECTOR<CCString> C_matu;
		VECTOR<CCString> C_matuDef(0);
		// error
		static int error;
		static char* reason = "";

		XL_readStrCell(XL_Crv1,C_Crv1," ARM_ERR:  curve1 id: object expected",C_result);
		XL_readStrCell(XL_BSCrv1,C_BSCrv1," ARM_ERR: BS curve1 id: object expected",C_result);
		XL_readStrCell(XL_Crv2,C_Crv2," ARM_ERR: curve2 id: object expected",C_result);
		XL_readStrCell(XL_BSCrv2,C_BSCrv2," ARM_ERR: BS curve2 id: object expected",C_result);
		XL_readStrCellWD(XL_RetSprds,C_RetSprds,"N"," ARM_ERR: RetSprds: string expected",C_result);
		XL_readStrVectorWD(XL_matuVec,C_matu,C_matuDef," ARM_ERR: maturities: array of string expected",DOUBLE_TYPE,C_result);
		XL_readStrCellWD(XL_BSasSprds,C_BSasSprds,"N"," ARM_ERR: BSasSprds: string expected",C_result);
		
		if(C_BSasSprds == "Y")
		{
			flagInputAsSprds = 1;
		}
		else if(C_BSasSprds == "N")
		{
			flagInputAsSprds = 0;
		}
		else
		{
			throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
							   "String \"Y\" or \"N\" Expected  BS_AsSpreads");
		}

		if(C_RetSprds == "Y")
		{
			flagRetSprds = 1;
		}
		else if(C_RetSprds == "N")
		{
			flagRetSprds = 0;
		}
		else
		{
			throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
							   "String \"Y\" or \"N\" Expected  ReturnSprds");
		}

		long retCode;
		long objId;
		CCString prevClass;
		
		CCString curClass = LOCAL_ZERO_CURVE_LIN_CLASS;
		CCString stringId;
		
		retCode = ARMLOCAL_GenerateBasisAdjCurve(LocalGetNumObjectId(C_Crv1),
													 LocalGetNumObjectId(C_BSCrv1),
													 LocalGetNumObjectId(C_Crv2),
													 LocalGetNumObjectId(C_BSCrv2),
													 flagInputAsSprds,
													 flagRetSprds,
													 C_matu,
													 C_result);			
		if(retCode == ARM_OK)
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure! in Local_GenerateBasisAdjCurve" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}



__declspec(dllexport) LPXLOPER WINAPI Local_GenForwardYield (LPXLOPER XL_curve,
															 LPXLOPER XL_mat,
															 LPXLOPER XL_tenor,
															 LPXLOPER XL_isSwapRate,
															 LPXLOPER XL_decompFreq,
															 LPXLOPER XL_daycount)
{
	ADD_LOG("Local_GenForwardYield ");
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();

	// C variable
	CCString C_curve;
	double C_matu1;
	double C_matu2;

	CCString C_isSwapRate;
	long isSwapRateId;

	CCString C_decompFreq;
	long decompFreqId;
	
	CCString C_daycount;
	long daycountId;

	// error
	static int error;
	static char* reason = "";

	XL_readStrCell(XL_curve,C_curve," ARM_ERR: curve id: object expected",C_result);
	XL_readNumCell(XL_mat,C_matu1," ARM_ERR: maturity : numeric expected",C_result);
	XL_readNumCell(XL_tenor,C_matu2," ARM_ERR: maturity 2: numeric expected",C_result);
	XL_readStrCellWD(XL_isSwapRate,C_isSwapRate,"YES"," ARM_ERR: isSwapRate: string expected",C_result);
	XL_readStrCellWD(XL_decompFreq,C_decompFreq,"-9999"," ARM_ERR: adjustment: string expected",C_result);
	XL_readStrCellWD(XL_daycount,C_daycount,"-1"," ARM_ERR: adjustment: string expected",C_result);

	if((isSwapRateId = ARM_ConvYesOrNo(C_isSwapRate, C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	if (strcmp(C_decompFreq,"-9999") == 0)
	{
		decompFreqId = -1;
	}
	else
	{
		if((decompFreqId = ARM_ConvCompMeth (C_decompFreq, C_result)) == ARM_DEFAULT_ERR)
		{
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}
	}
	
	daycountId = ARM_ConvDayCount(C_daycount);

	long retCode;

	retCode = ARMLOCAL_GenForwardYield (LocalGetNumObjectId (C_curve),
										C_matu1,
										C_matu2,
										isSwapRateId,
										decompFreqId,
										daycountId,
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_GenForwardYield" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


__declspec(dllexport) LPXLOPER WINAPI Local_OldZcCurve_Common(
		LPXLOPER XL_zcCurve,
		LPXLOPER XL_asOfDate,
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

		double C_asOfDate;
		XL_readNumCell(XL_asOfDate,C_asOfDate," ARM_ERR: as of date: date expected",C_result);
	
		long C_zcCurveId;
		XL_GETOBJID( XL_zcCurve,	C_zcCurveId,	" ARM_ERR: zc Curve Id: Object expected",	C_result);

		exportFunc2Args< long, double >
			ourFunc(C_zcCurveId, C_asOfDate, ARMLOCAL_OldZcCurve);			

		bool PersistentInXL = true;

		/// call the general function
		fillXL_Result( LOCAL_ZERO_CURVE_LIN_CLASS, ourFunc, C_result, XL_result, PersistentInXL );
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_CRAQuantoCalculator_Create" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


__declspec(dllexport) LPXLOPER WINAPI Local_OldZcCurve(
	LPXLOPER XL_zcCurve,
	LPXLOPER XL_asOfDate)
{
	bool PersistentInXL = true;
	
	return Local_OldZcCurve_Common(
			XL_zcCurve,
			XL_asOfDate,
			PersistentInXL);
}

__declspec(dllexport) LPXLOPER WINAPI Local_PXL_OldZcCurve(
	LPXLOPER XL_zcCurve,
	LPXLOPER XL_asOfDate)
{
	bool PersistentInXL = false;
			
	return Local_OldZcCurve_Common(
			XL_zcCurve,
			XL_asOfDate,
			PersistentInXL);
}


__declspec(dllexport) LPXLOPER WINAPI Local_ARM_FixingSched(LPXLOPER XL_asOfDate,
															LPXLOPER XL_LiborFixing,
															LPXLOPER XL_FXFixing)
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
		

		// error
		static int error;
		static char* reason = "";

		/// this is used by macros 
		/// and therefore this has to be defined
		double C_asOfDate;
		XL_readNumCell(XL_asOfDate,C_asOfDate," ARM_ERR: as of date: date expected",C_result);


		vector<CCString> C_LiborFixing;
		XL_readStrVector(XL_LiborFixing,C_LiborFixing," ARM_ERR: LiborFixing : string expected key+curve",DOUBLE_TYPE,C_result);


		vector<string> C_LiborFixingKeys(C_LiborFixing.size() / 2);
		vector<long>   C_LiborFixingCurve(C_LiborFixing.size() / 2);

		for (int i = 0; i < C_LiborFixing.size()/2 ; i++)
		{
			C_LiborFixingKeys[i] = C_LiborFixing[2*i];
			C_LiborFixingCurve[i]  = LocalGetNumObjectId(C_LiborFixing[2*i+1]);
		}

		vector<CCString> C_FXFixing;
		XL_readStrVector(XL_FXFixing,C_FXFixing," ARM_ERR: FXFixing : string expected key+curve",DOUBLE_TYPE,C_result);


		vector<string> C_FXFixingKeys(C_FXFixing.size() / 2);
		vector<long>   C_FXFixingCurve(C_FXFixing.size() / 2);

		for (i = 0; i < C_FXFixing.size()/2 ; i++)
		{
			C_FXFixingKeys[i] = C_FXFixing[2*i];
			C_FXFixingCurve[i]  = LocalGetNumObjectId(C_FXFixing[2*i+1]);
		}

		long objId;
	    CCString prevClass;

	    CCString curClass = LOCAL_FIXING_SCHED_CLASS;
	    CCString stringId = GetLastCurCellEnvValue ();

		long retCode;

		if(!stringId)
		{
			retCode = ARMLOCAL_FixingSched (C_asOfDate,
											C_LiborFixingKeys,
											C_FXFixingKeys,
											C_LiborFixingCurve,
											C_FXFixingCurve,
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
				retCode = ARMLOCAL_FixingSched (C_asOfDate,
										C_LiborFixingKeys,
										C_FXFixingKeys,
										C_LiborFixingCurve,
										C_FXFixingCurve,
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
				retCode = ARMLOCAL_FixingSched (C_asOfDate,
										C_LiborFixingKeys,
										C_FXFixingKeys,
										C_LiborFixingCurve,
										C_FXFixingCurve,
										C_result);
					
				if ( retCode == ARM_OK )
				{
					objId = C_result.getLong ();

					LocalSetCurCellEnvValue (curClass, objId); 

					stringId = LocalMakeObjectId (objId, curClass);
				}
			}
		}



		if (retCode == ARM_OK)
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


	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_ARM_FixingSched" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}

__declspec(dllexport) LPXLOPER WINAPI Local_PXL_ARM_FixingSched(LPXLOPER XL_asOfDate,
																LPXLOPER XL_LiborFixing,
																LPXLOPER XL_FXFixing)
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
		

		// error
		static int error;
		static char* reason = "";

		/// this is used by macros 
		/// and therefore this has to be defined
		double C_asOfDate;
		XL_readNumCell(XL_asOfDate,C_asOfDate," ARM_ERR: as of date: date expected",C_result);


		vector<CCString> C_LiborFixing;
		XL_readStrVector(XL_LiborFixing,C_LiborFixing," ARM_ERR: LiborFixing : string expected key+curve",DOUBLE_TYPE,C_result);


		vector<string> C_LiborFixingKeys(C_LiborFixing.size() / 2);
		vector<long>   C_LiborFixingCurve(C_LiborFixing.size() / 2);

		for (int i = 0; i < C_LiborFixing.size()/2 ; i++)
		{
			C_LiborFixingKeys[i] = C_LiborFixing[2*i];
			C_LiborFixingCurve[i]  = LocalGetNumObjectId(C_LiborFixing[2*i+1]);
		}

		vector<CCString> C_FXFixing;
		XL_readStrVector(XL_FXFixing,C_FXFixing," ARM_ERR: FXFixing : string expected key+curve",DOUBLE_TYPE,C_result);


		vector<string> C_FXFixingKeys(C_FXFixing.size() / 2);
		vector<long>   C_FXFixingCurve(C_FXFixing.size() / 2);

		for (i = 0; i < C_FXFixing.size()/2 ; i++)
		{
			C_FXFixingKeys[i] = C_FXFixing[2*i];
			C_FXFixingCurve[i]  = LocalGetNumObjectId(C_FXFixing[2*i+1]);
		}


		long retCode;

		long objId;
		CCString curClass = LOCAL_FIXING_SCHED_CLASS;

		retCode = ARMLOCAL_FixingSched (C_asOfDate,
										C_LiborFixingKeys,
										C_FXFixingKeys,
										C_LiborFixingCurve,
										C_FXFixingCurve,
										C_result);

		if(retCode == ARM_OK)
		{
			objId = C_result.getLong ();
		}

		if ( retCode == ARM_OK )
		{
			FreeCurCellErr ();
			XL_result.xltype = xltypeStr;
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_PXL_ARM_FixingSched" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


__declspec(dllexport) LPXLOPER WINAPI Local_ARM_GetFixingSchedFromSummit(LPXLOPER XL_ListOfKeys,
																		 LPXLOPER XL_AsOf,
																		 LPXLOPER XL_Source,
																		 LPXLOPER XL_dateStrip)
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
		
		// error
		static int error;
		static char* reason = "";

		/// this is used by macros 
		/// and therefore this has to be defined
		double C_aSdate;

		XL_readNumCell(XL_AsOf,C_aSdate," ARM_ERR: as of date: date expected",C_result);

		VECTOR<CCString> C_ListOfKeys;
		XL_readStrVector(XL_ListOfKeys,C_ListOfKeys," ARM_ERR: mkt Tag: array expected",DOUBLE_TYPE,C_result);

		vector<string> Libor_and_FX_Keys( C_ListOfKeys.size() );

		Libor_and_FX_Keys[0] = CCSTringToSTLString(C_ListOfKeys[0]);
		Libor_and_FX_Keys[1] = CCSTringToSTLString(C_ListOfKeys[1]);

		CCString C_Source;
		XL_readStrCell(XL_Source,C_Source," ARM_ERR: Source: string expected",C_result);

		CCString C_dateStripIdStr;
		XL_readStrCell( XL_dateStrip, C_dateStripIdStr, " ARM_ERR: DateStrip Id: Object expected",C_result);
		long C_dateStripId = LocalGetNumObjectId(C_dateStripIdStr);

		long objId;
	    CCString prevClass;

	    CCString curClass = LOCAL_FIXING_SCHED_CLASS;
	    CCString stringId = GetLastCurCellEnvValue ();

		long retCode;

		if(!stringId)
		{
			retCode = ARMLOCAL_GetFixingSchedFromSummit(Libor_and_FX_Keys,	
														C_aSdate,
														C_Source,
														C_dateStripId,
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
				retCode = ARMLOCAL_GetFixingSchedFromSummit(Libor_and_FX_Keys,	
															C_aSdate,
															C_Source,
															C_dateStripId,
															C_result);
					
				if(retCode == ARM_OK)
				{
					LocalSetCurCellEnvValue (curClass, objId); 

					stringId = LocalMakeObjectId (objId, curClass);
				}
			}
			else
			{
				FreeCurCellContent ();
				retCode = ARMLOCAL_GetFixingSchedFromSummit(Libor_and_FX_Keys,	
															C_aSdate,
															C_Source,
															C_dateStripId,
															C_result);

					
				if ( retCode == ARM_OK )
				{
					objId = C_result.getLong ();

					LocalSetCurCellEnvValue (curClass, objId); 

					stringId = LocalMakeObjectId (objId, curClass);
				}
			}
		}



		if (retCode == ARM_OK)
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



	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_ARM_GetFixingSchedFromSummit" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}

