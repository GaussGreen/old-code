#pragma warning(disable : 4786)
#pragma warning(disable : 4005)

#include <libCCxll\CCxll.h>


#include "ARM_local_interglob.h"
#include "ARM_local_interface.h"
#include <ARM\libarm_local\ARM_local_class.h>
#include <ARM\libarm_local\ARM_local_util.h>
#include "ARM_xl_trycatch_local.h"

#include "util\tech_macro.h"

__declspec(dllexport) LPXLOPER WINAPI Local_GetCorrelInst(LPXLOPER XL_date1,
														  LPXLOPER XL_date2,
														  LPXLOPER XL_ccy1,
														  LPXLOPER XL_index1,
														  LPXLOPER XL_expiry1,
														  LPXLOPER XL_tenor1,
														  LPXLOPER XL_curve1_ccy1,	
														  LPXLOPER XL_curve2_ccy1,	
														  LPXLOPER XL_nbmonths_curve1_ccy1,	
														  LPXLOPER XL_ccy2,
														  LPXLOPER XL_index2,
														  LPXLOPER XL_expiry2,
														  LPXLOPER XL_tenor2,
														  LPXLOPER XL_curve1_ccy2,	
														  LPXLOPER XL_curve2_ccy2,	
														  LPXLOPER XL_nbmonths_curve1_ccy2,	
														  LPXLOPER XL_type,
														  LPXLOPER XL_lambda,
														  LPXLOPER XL_precision,
														  LPXLOPER XL_ccy)
{
	ADD_LOG("Local_GetCorrelInst");
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
	double C_date1;
	double C_date2;

	CCString C_ccy1;
	CCString C_index1;
	CCString C_expiry1;
	CCString C_tenor1;

	CCString C_ccy2;
	CCString C_index2;
	CCString C_expiry2;
	CCString C_tenor2;

	CCString C_curve1_ccy1;
    CCString C_curve2_ccy1;	
    CCString C_nbmonths_curve1_ccy1;
	CCString C_nbmonths_curve1_ccy1_default="0M";

	CCString C_curve1_ccy2;
    CCString C_curve2_ccy2;	
    CCString C_nbmonths_curve1_ccy2;
	CCString C_nbmonths_curve1_ccy2_default="0M";

	CCString C_type;
	long typeId;

	double C_lambda;
	double C_lambda_default = 0.0;

	double C_precision;
	double C_precision_default = 0;

	CCString C_ccy;

	// error
	static int error;
	static char* reason = "";

	XL_readNumCell(XL_date1,C_date1," ARM_ERR: date 1: date expected",C_result);
	XL_readNumCell(XL_date2,C_date2," ARM_ERR: date 2: date expected",C_result);
	XL_readStrCell(XL_ccy1,C_ccy1," ARM_ERR: currency1: string expected",C_result);
	XL_readStrCell(XL_index1,C_index1," ARM_ERR: index1: string expected",C_result);
	XL_readStrCell(XL_expiry1,C_expiry1," ARM_ERR: expiry 1: string expected (1M, 2M, ... 20Y)",C_result);
	XL_readStrCell(XL_tenor1,C_tenor1," ARM_ERR: tenor1: string expected",C_result);
	XL_readStrCell(XL_ccy2,C_ccy2," ARM_ERR: currency2: string expected",C_result);
	XL_readStrCell(XL_index2,C_index2," ARM_ERR: index2: string expected",C_result);
	XL_readStrCell(XL_expiry2,C_expiry2," ARM_ERR: expiry 2: string expected (1M, 2M, ... 20Y)",C_result);
	XL_readStrCell(XL_tenor2,C_tenor2," ARM_ERR: tenor2: string expected",C_result);
	XL_readStrCell(XL_type,C_type," ARM_ERR: type: string expected (NOR or LOGNOR)",C_result);

	XL_readStrCell(XL_curve1_ccy1,C_curve1_ccy1," ARM_ERR: tenor2: string expected",C_result);
	XL_readStrCellWD(XL_curve2_ccy1,C_curve2_ccy1, C_curve1_ccy1, " ARM_ERR: tenor2: string expected",C_result);
	XL_readStrCellWD(XL_nbmonths_curve1_ccy1,C_nbmonths_curve1_ccy1, C_nbmonths_curve1_ccy1_default," ARM_ERR: tenor2: string expected",C_result);
	XL_readStrCell(XL_curve1_ccy2,C_curve1_ccy2," ARM_ERR: tenor2: string expected",C_result);
	XL_readStrCellWD(XL_curve2_ccy2,C_curve2_ccy2,C_curve1_ccy2," ARM_ERR: tenor2: string expected",C_result);
	XL_readStrCellWD(XL_nbmonths_curve1_ccy2,C_nbmonths_curve1_ccy2, C_nbmonths_curve1_ccy2_default," ARM_ERR: tenor2: string expected",C_result);
	XL_readNumCellWD(XL_lambda,C_lambda,C_lambda_default," ARM_ERR: lambda: numeric expected",C_result);
	XL_readNumCellWD(XL_precision,C_precision,C_precision_default," ARM_ERR: precision: numeric expected",C_result);
	XL_readStrCellWD(XL_ccy,C_ccy,"EUR"," ARM_ERR: currency: string expected",C_result);

	if((typeId = ARM_ConvCalcMod (C_type, C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

    long retCode;

	retCode = ARMLOCAL_GetCorrelInst(C_date1,
									 C_date2,
									 C_ccy1,
									 C_index1,
									 C_expiry1,
									 C_tenor1,
									 C_ccy2,
									 C_index2,
									 C_expiry2,
									 C_tenor2,
									 C_curve1_ccy1,
									 C_curve2_ccy1,	
									 C_nbmonths_curve1_ccy1,
									 C_curve1_ccy2,
									 C_curve2_ccy2,	
									 C_nbmonths_curve1_ccy2,
									 typeId,
									 C_lambda,
									 C_precision,
									 C_ccy,
									 C_result);


	if ( retCode == ARM_OK )
	{
		int nbrows = 4;
		int nbcolumns = 1;
	
		FreeCurCellErr ();
		XL_result.xltype = xltypeMulti;
		XL_result.val.array.columns = nbcolumns;
		XL_result.val.array.rows = nbrows; 
		XL_result.val.array.lparray = pxArray = (LPXLOPER)GlobalAlloc (GMEM_ZEROINIT, nbrows * nbcolumns * sizeof (XLOPER));

		for(int i = 0; i < nbrows; i++)
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_GetCorrelInst" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


__declspec(dllexport) LPXLOPER WINAPI Local_GetMoyCorrel(LPXLOPER XL_date1,
														 LPXLOPER XL_date2,
														 LPXLOPER XL_ccy1,
														 LPXLOPER XL_index1,
														 LPXLOPER XL_expiry1,
														 LPXLOPER XL_tenor1,
													     LPXLOPER XL_curve1_ccy1,	
													     LPXLOPER XL_curve2_ccy1,	
													     LPXLOPER XL_nbmonths_curve1_ccy1,	
														 LPXLOPER XL_ccy2,
														 LPXLOPER XL_index2,
														 LPXLOPER XL_expiry2,
														 LPXLOPER XL_tenor2,
													     LPXLOPER XL_curve1_ccy2,	
													     LPXLOPER XL_curve2_ccy2,	
													     LPXLOPER XL_nbmonths_curve1_ccy2,	
														 LPXLOPER XL_type,
														 LPXLOPER XL_lambda,
														 LPXLOPER XL_precision,
														 LPXLOPER XL_ccy)
{
	ADD_LOG("Local_GetMoyCorrel");
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
	double C_date1;
	double C_date2;

	CCString C_ccy1;
	CCString C_index1;
	CCString C_expiry1;
	CCString C_tenor1;

	CCString C_ccy2;
	CCString C_index2;
	CCString C_expiry2;
	CCString C_tenor2;

	CCString C_type;
	long typeId;

	CCString C_curve1_ccy1;
    CCString C_curve2_ccy1;	
    CCString C_nbmonths_curve1_ccy1;
	CCString C_nbmonths_curve1_ccy1_default="0M";

	CCString C_curve1_ccy2;
    CCString C_curve2_ccy2;	
    CCString C_nbmonths_curve1_ccy2;
	CCString C_nbmonths_curve1_ccy2_default="0M";

	double C_lambda;
	double C_lambda_default = 0.0;

	double C_precision;
	double C_precision_default = 0;

	CCString C_ccy;

	// error
	static int error;
	static char* reason = "";

	XL_readNumCell(XL_date1,C_date1," ARM_ERR: date 1: date expected",C_result);
	XL_readNumCell(XL_date2,C_date2," ARM_ERR: date 2: date expected",C_result);
	XL_readStrCell(XL_ccy1,C_ccy1," ARM_ERR: currency1: string expected",C_result);
	XL_readStrCell(XL_index1,C_index1," ARM_ERR: index1: string expected",C_result);
	XL_readStrCell(XL_expiry1,C_expiry1," ARM_ERR: expiry 1: string expected (1M, 2M, ... 20Y)",C_result);
	XL_readStrCell(XL_tenor1,C_tenor1," ARM_ERR: tenor1: string expected",C_result);
	XL_readStrCell(XL_ccy2,C_ccy2," ARM_ERR: currency2: string expected",C_result);
	XL_readStrCell(XL_index2,C_index2," ARM_ERR: index2: string expected",C_result);
	XL_readStrCell(XL_expiry2,C_expiry2," ARM_ERR: expiry 2: string expected (1M, 2M, ... 20Y)",C_result);
	XL_readStrCell(XL_tenor2,C_tenor2," ARM_ERR: tenor2: string expected",C_result);

	XL_readStrCell(XL_curve1_ccy1,C_curve1_ccy1," ARM_ERR: tenor2: string expected",C_result);
	XL_readStrCell(XL_curve2_ccy1,C_curve2_ccy1," ARM_ERR: tenor2: string expected",C_result);
	XL_readStrCellWD(XL_nbmonths_curve1_ccy1,C_nbmonths_curve1_ccy1, C_nbmonths_curve1_ccy1_default," ARM_ERR: tenor2: string expected",C_result);
	XL_readStrCell(XL_curve1_ccy2,C_curve1_ccy2," ARM_ERR: tenor2: string expected",C_result);
	XL_readStrCell(XL_curve2_ccy2,C_curve2_ccy2," ARM_ERR: tenor2: string expected",C_result);
	XL_readStrCellWD(XL_nbmonths_curve1_ccy2,C_nbmonths_curve1_ccy2, C_nbmonths_curve1_ccy2_default," ARM_ERR: tenor2: string expected",C_result);

	XL_readStrCell(XL_type,C_type," ARM_ERR: type: string expected (NOR or LOGNOR)",C_result);
	XL_readNumCellWD(XL_lambda,C_lambda,C_lambda_default," ARM_ERR: lambda: numeric expected",C_result);
	XL_readNumCellWD(XL_precision,C_precision,C_precision_default," ARM_ERR: precision: numeric expected",C_result);
	XL_readStrCellWD(XL_ccy,C_ccy,"EUR"," ARM_ERR: currency: string expected",C_result);

	if((typeId = ARM_ConvCalcMod (C_type, C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

    long retCode;

	retCode = ARMLOCAL_GetMoyCorrel(C_date1,
									C_date2,
									C_ccy1,
									C_index1,
									C_expiry1,
									C_tenor1,
									C_ccy2,
									C_index2,
									C_expiry2,
									C_tenor2,
									C_curve1_ccy1,
									C_curve2_ccy1,	
									C_nbmonths_curve1_ccy1,
									C_curve1_ccy2,
									C_curve2_ccy2,	
									C_nbmonths_curve1_ccy2,
									typeId,
									C_lambda,
									C_precision,
									C_ccy,
									C_result);

	if ( retCode == ARM_OK )
	{
		int nbrows = 4;
		int nbcolumns = 1;
	
		FreeCurCellErr ();
		XL_result.xltype = xltypeMulti;
		XL_result.val.array.columns = nbcolumns;
		XL_result.val.array.rows = nbrows; 
		XL_result.val.array.lparray = pxArray = (LPXLOPER)GlobalAlloc (GMEM_ZEROINIT, nbrows * nbcolumns * sizeof (XLOPER));

		for(int i = 0; i < nbrows; i++)
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_GetMoyCorrel" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}



__declspec(dllexport) LPXLOPER WINAPI Local_GetCorrelQuanto(LPXLOPER XL_date1,
															LPXLOPER XL_date2,
															LPXLOPER XL_ccy,
															LPXLOPER XL_index,
															LPXLOPER XL_expiry,
															LPXLOPER XL_tenor,
															LPXLOPER XL_cvname1,
															LPXLOPER XL_cvname2,
															LPXLOPER XL_switchinmonth,
															LPXLOPER XL_domccy,
															LPXLOPER XL_domindex,
															LPXLOPER XL_forccy,
															LPXLOPER XL_forindex,
															LPXLOPER XL_type,
															LPXLOPER XL_lambda,
															LPXLOPER XL_precision,
															LPXLOPER XL_calccy,
															LPXLOPER XL_fwdOrNot)
{
	ADD_LOG("Local_GetCorrelQuanto");
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
	double C_date1;
	double C_date2;

	CCString C_ccy;
	CCString C_index;
	CCString C_expiry;
	CCString C_tenor;

	CCString C_cvname1;
    CCString C_cvname2;	
    CCString C_switchinmonth;
	CCString C_switchinmonth_default="0M";

	CCString C_domccy;
	CCString C_domindex;
	CCString C_forccy;
	CCString C_forindex;

	CCString C_type;
	long typeId;

	double C_lambda;
	double C_lambda_default = 0.0;

	double C_precision;
	double C_precision_default = 0;

	CCString C_calccy;

	double C_fwdOrNot;
	double C_fwdOrNot_default = 1;

	// error
	static int error;
	static char* reason = "";

	XL_readNumCell(XL_date1,C_date1," ARM_ERR: date 1: date expected",C_result);
	XL_readNumCell(XL_date2,C_date2," ARM_ERR: date 2: date expected",C_result);
	XL_readStrCell(XL_ccy,C_ccy," ARM_ERR: Yield currency: string expected",C_result);
	XL_readStrCell(XL_index,C_index," ARM_ERR: Yield index: string expected",C_result);
	XL_readStrCell(XL_expiry,C_expiry," ARM_ERR: Yield fixing: string expected (1M, 2M, ... 20Y)",C_result);
	XL_readStrCell(XL_tenor,C_tenor," ARM_ERR: Yield tenor: string expected",C_result);

	XL_readStrCell(XL_cvname1,C_cvname1," ARM_ERR: Yield cvname1: string expected",C_result);
	XL_readStrCellWD(XL_cvname2,C_cvname2,C_cvname1," ARM_ERR: Yield cvname2: string expected",C_result);
	XL_readStrCellWD(XL_switchinmonth,C_switchinmonth, C_switchinmonth_default," ARM_ERR: switch in month: string expected (1M,1Y...)",C_result);

	XL_readStrCell(XL_domccy,C_domccy," ARM_ERR: dom currency: string expected",C_result);
	XL_readStrCell(XL_domindex,C_domindex," ARM_ERR: dom index: string expected",C_result);
	XL_readStrCell(XL_forccy,C_forccy," ARM_ERR: for currency: string expected",C_result);
	XL_readStrCell(XL_forindex,C_forindex," ARM_ERR: for index: string expected",C_result);

	XL_readStrCell(XL_type,C_type," ARM_ERR: type: string expected (NOR or LOGNOR)",C_result);

	XL_readNumCellWD(XL_lambda,C_lambda,C_lambda_default," ARM_ERR: lambda: numeric expected",C_result);
	XL_readNumCellWD(XL_precision,C_precision,C_precision_default," ARM_ERR: precision: numeric expected",C_result);
	XL_readStrCellWD(XL_calccy,C_calccy,"EUR"," ARM_ERR: currency: string expected",C_result);
	XL_readNumCellWD(XL_fwdOrNot,C_fwdOrNot,C_fwdOrNot_default," ARM_ERR: fwd or not: numeric expected",C_result);

	if((typeId = ARM_ConvCalcMod (C_type, C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

    long retCode;

	retCode = ARMLOCAL_GetCorrelQuanto(C_date1,
									   C_date2,
									   C_ccy,
									   C_index,
									   C_expiry,
									   C_tenor,
									   C_cvname1,
									   C_cvname2,
									   C_switchinmonth,
									   C_domccy,
									   C_domindex,
									   C_forccy,
									   C_forindex,
									   typeId,
									   C_lambda,
									   C_precision,
									   C_calccy,
									   (long)C_fwdOrNot,
									   C_result);


	if ( retCode == ARM_OK )
	{
		int nbrows = 4;
		int nbcolumns = 1;
	
		FreeCurCellErr ();
		XL_result.xltype = xltypeMulti;
		XL_result.val.array.columns = nbcolumns;
		XL_result.val.array.rows = nbrows; 
		XL_result.val.array.lparray = pxArray = (LPXLOPER)GlobalAlloc (GMEM_ZEROINIT, nbrows * nbcolumns * sizeof (XLOPER));

		for(int i = 0; i < nbrows; i++)
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_GetCorrelQuanto" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


__declspec(dllexport) LPXLOPER WINAPI Local_GetMoyCorrelQuanto(LPXLOPER XL_date1,
															   LPXLOPER XL_date2,
															   LPXLOPER XL_ccy,
															   LPXLOPER XL_index,
															   LPXLOPER XL_expiry,
															   LPXLOPER XL_tenor,
															   LPXLOPER XL_cvname1,
															   LPXLOPER XL_cvname2,
															   LPXLOPER XL_switchinmonth,
															   LPXLOPER XL_domccy,
															   LPXLOPER XL_domindex,
															   LPXLOPER XL_forccy,
															   LPXLOPER XL_forindex,
															   LPXLOPER XL_type,
															   LPXLOPER XL_lambda,
															   LPXLOPER XL_precision,
															   LPXLOPER XL_calccy,
															   LPXLOPER XL_fwdOrNot)
{
	ADD_LOG("Local_GetMoyCorrelQuanto");
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
	double C_date1;
	double C_date2;

	CCString C_ccy;
	CCString C_index;
	CCString C_expiry;
	CCString C_tenor;

	CCString C_cvname1;
    CCString C_cvname2;	
    CCString C_switchinmonth;
	CCString C_switchinmonth_default="0M";

	CCString C_domccy;
	CCString C_domindex;
	CCString C_forccy;
	CCString C_forindex;

	CCString C_type;
	long typeId;

	double C_lambda;
	double C_lambda_default = 0.0;

	double C_precision;
	double C_precision_default = 0;

	CCString C_calccy;

	double C_fwdOrNot;
	double C_fwdOrNot_default = 1;

	// error
	static int error;
	static char* reason = "";

	XL_readNumCell(XL_date1,C_date1," ARM_ERR: date 1: date expected",C_result);
	XL_readNumCell(XL_date2,C_date2," ARM_ERR: date 2: date expected",C_result);
	XL_readStrCell(XL_ccy,C_ccy," ARM_ERR: Yield currency: string expected",C_result);
	XL_readStrCell(XL_index,C_index," ARM_ERR: Yield index: string expected",C_result);
	XL_readStrCell(XL_expiry,C_expiry," ARM_ERR: Yield fixing: string expected (1M, 2M, ... 20Y)",C_result);
	XL_readStrCell(XL_tenor,C_tenor," ARM_ERR: Yield tenor: string expected",C_result);

	XL_readStrCell(XL_cvname1,C_cvname1," ARM_ERR: Yield cvname1: string expected",C_result);
	XL_readStrCellWD(XL_cvname2,C_cvname2,C_cvname1," ARM_ERR: Yield cvname2: string expected",C_result);
	XL_readStrCellWD(XL_switchinmonth,C_switchinmonth, C_switchinmonth_default," ARM_ERR: switch in month: string expected (1M,1Y...)",C_result);

	XL_readStrCell(XL_domccy,C_domccy," ARM_ERR: dom currency: string expected",C_result);
	XL_readStrCell(XL_domindex,C_domindex," ARM_ERR: dom index: string expected",C_result);
	XL_readStrCell(XL_forccy,C_forccy," ARM_ERR: for currency: string expected",C_result);
	XL_readStrCell(XL_forindex,C_forindex," ARM_ERR: for index: string expected",C_result);

	XL_readStrCell(XL_type,C_type," ARM_ERR: type: string expected (NOR or LOGNOR)",C_result);

	XL_readNumCellWD(XL_lambda,C_lambda,C_lambda_default," ARM_ERR: lambda: numeric expected",C_result);
	XL_readNumCellWD(XL_precision,C_precision,C_precision_default," ARM_ERR: precision: numeric expected",C_result);
	XL_readStrCellWD(XL_calccy,C_calccy,"EUR"," ARM_ERR: currency: string expected",C_result);
	XL_readNumCellWD(XL_fwdOrNot,C_fwdOrNot,C_fwdOrNot_default," ARM_ERR: fwd or not: numeric expected",C_result);

	if((typeId = ARM_ConvCalcMod (C_type, C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

    long retCode;

	retCode = ARMLOCAL_GetMoyCorrelQuanto(C_date1,
										  C_date2,
										  C_ccy,
										  C_index,
										  C_expiry,
										  C_tenor,
										  C_cvname1,
										  C_cvname2,
										  C_switchinmonth,
										  C_domccy,
										  C_domindex,
										  C_forccy,
										  C_forindex,
										  typeId,
										  C_lambda,
										  C_precision,
										  C_calccy,
										  (long)C_fwdOrNot,
										  C_result);


	if ( retCode == ARM_OK )
	{
		int nbrows = 4;
		int nbcolumns = 1;
	
		FreeCurCellErr ();
		XL_result.xltype = xltypeMulti;
		XL_result.val.array.columns = nbcolumns;
		XL_result.val.array.rows = nbrows; 
		XL_result.val.array.lparray = pxArray = (LPXLOPER)GlobalAlloc (GMEM_ZEROINIT, nbrows * nbcolumns * sizeof (XLOPER));

		for(int i = 0; i < nbrows; i++)
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_GetMoyCorrelQuanto" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}



__declspec(dllexport) LPXLOPER WINAPI Local_ARM_GetDealsFromSummitFilter (LPXLOPER XL_filter)
{
	ADD_LOG("Local_ARM_GetDealsFromSummitFilter ");
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
	CCString C_filter;

	VECTOR<CCString> listDeals;
	VECTOR<CCString> listTypes;

	// error
	static int error;
	static char* reason = "";

	XL_readStrCell(XL_filter,C_filter," ARM_ERR: filter: string expected",C_result);

	long retCode;

	retCode = ARMLOCAL_GetDealsFromSummitFilter (C_filter, listDeals, listTypes, C_result);

	if(retCode == ARM_OK)
	{
		int nbrows = listDeals.size();

		int nbcolumns = 2;
		
		FreeCurCellErr ();
		XL_result.xltype = xltypeMulti;
		XL_result.val.array.columns = nbcolumns;
		XL_result.val.array.rows = nbrows; 
		XL_result.val.array.lparray = pxArray = (LPXLOPER)GlobalAlloc (GMEM_ZEROINIT, nbrows * nbcolumns * sizeof (XLOPER));

		for(int i = 0; i < nbrows; i++)
		{
			pxArray[XL_Coordonnate2Rank (i, 0, nbcolumns)].xltype = xltypeStr;
			pxArray[XL_Coordonnate2Rank (i, 0, nbcolumns)].val.str = XL_StrC2StrPascal (listDeals[i]);
			pxArray[XL_Coordonnate2Rank (i, 0, nbcolumns)].xltype |= xlbitDLLFree;
			pxArray[XL_Coordonnate2Rank (i, 1, nbcolumns)].xltype = xltypeStr;
			pxArray[XL_Coordonnate2Rank (i, 1, nbcolumns)].val.str = XL_StrC2StrPascal (listTypes[i]);
			pxArray[XL_Coordonnate2Rank (i, 1, nbcolumns)].xltype |= xlbitDLLFree;
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_ARM_GetDealsFromSummitFilter" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


__declspec(dllexport) LPXLOPER WINAPI Local_ARM_GetAsOfVolOrRate(LPXLOPER XL_AsOfDate,
																 LPXLOPER XL_ccy,
																 LPXLOPER XL_index,
																 LPXLOPER XL_cvName,
																 LPXLOPER XL_expiry,
																 LPXLOPER XL_matu,
																 LPXLOPER XL_yieldOrVol,
																 LPXLOPER XL_calcMod,
																 LPXLOPER XL_volType)
{
	ADD_LOG("Local_ARM_GetAsOfVolOrRate");
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;

	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();

	// C variable
	double C_asOfDate;

	CCString C_ccy;
	CCString C_index;
	CCString C_cvName;

	CCString C_expiry;
	CCString C_matu;

	CCString C_yieldOrVol;
	long yieldOrVolId;

	CCString C_calcMod;
	long calcModId;

	CCString C_volType;
	
	VECTOR<CCString> matu;
	VECTOR<double> rate;
	
	// error
	static int error;
	static char* reason = "";

	XL_readNumCell(XL_AsOfDate,C_asOfDate," ARM_ERR: As of Date: date expected",C_result);

	XL_readStrCell(XL_ccy,C_ccy," ARM_ERR: currency: string expected",C_result);
	XL_readStrCell(XL_index,C_index," ARM_ERR: index: string expected",C_result);
	XL_readStrCell(XL_cvName,C_cvName," ARM_ERR: Summit curve ID: string expected",C_result);

    XL_readStrCell(XL_expiry,C_expiry," ARM_ERR: expiry: string expected (2M, 5Y, ..)",
                   C_result);
	XL_readStrCell(XL_matu,C_matu," ARM_ERR: maturity: string expected (3M, 10Y, ..)",C_result);

    XL_readStrCellWD(XL_yieldOrVol, C_yieldOrVol, "Y",
       " ARM_ERR: Yield or Volatily type expected: string (Y, V) expected",C_result);
    
    XL_readStrCellWD(XL_calcMod, C_calcMod, "LOGNOR",
       " ARM_ERR: Calculation mode: string expected (LOGNOR, NOR, SQR)",C_result);

	XL_readStrCellWD(XL_volType, C_volType, "IRG",
       " ARM_ERR: Volatility type: string expected(IRG, SWOPT, ..)",C_result);	

	if((yieldOrVolId = ARM_ConvYieldOrVol (C_yieldOrVol, C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	if((calcModId = ARM_ConvCalcMod (C_calcMod, C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	long retCode;

	retCode = ARMLOCAL_GetAsOfVolOrRate(C_asOfDate,
										C_ccy,
										C_index,
										C_cvName,
										C_expiry,
										C_matu,
										yieldOrVolId,
										calcModId,
										C_volType,
										C_result);

	if ( retCode == ARM_OK )
	{
	   FreeCurCellErr();

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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_GetCorrelQuanto" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}
					        

__declspec(dllexport) LPXLOPER WINAPI Local_ARM_GetFwdRatesMatrix(LPXLOPER XL_AsOfDate,
																  LPXLOPER XL_ccy,
																  LPXLOPER XL_index,
																  LPXLOPER XL_cvName)
{
	ADD_LOG("Local_ARM_GetFwdRatesMatrix");
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
	double C_asOfDate;

	CCString C_ccy;
	CCString C_index;
	CCString C_cvName;

	// error
	static int error;
	static char* reason = "";

	XL_readNumCell(XL_AsOfDate,C_asOfDate," ARM_ERR: As of Date: date expected",C_result);

	XL_readStrCell(XL_ccy,C_ccy," ARM_ERR: currency: string expected",C_result);
	XL_readStrCell(XL_index,C_index," ARM_ERR: index: string expected",C_result);
	XL_readStrCell(XL_cvName,C_cvName," ARM_ERR: Summit curve ID: string expected",C_result);

	long retCode = ARMLOCAL_GetFwdRatesMatrix(C_asOfDate,
											  C_ccy,
											  C_index,
											  C_cvName,
											  C_result);

	if ( retCode == ARM_OK )
	{
		int nbrows = MATU_LINES_SIZE+1;
		int nbcolumns = MATU_COLS_SIZE+1;
	
		FreeCurCellErr ();
		XL_result.xltype = xltypeMulti;
		XL_result.val.array.columns = nbcolumns;
		XL_result.val.array.rows = nbrows; 
		XL_result.val.array.lparray = pxArray = (LPXLOPER)GlobalAlloc (GMEM_ZEROINIT, nbrows * nbcolumns * sizeof (XLOPER));

		pxArray[XL_Coordonnate2Rank (0, 0, nbcolumns)].xltype = xltypeStr;
		pxArray[XL_Coordonnate2Rank (0, 0, nbcolumns)].val.str = XL_StrC2StrPascal ("");
		pxArray[XL_Coordonnate2Rank (0, 0, nbcolumns)].xltype |= xlbitDLLFree;

		for (int i = 0; i < nbrows; i++)
		{
			for (int j = 0; j < nbcolumns; j++)
			{
				if ( (i != 0) || (j != 0) )
				{
					pxArray[XL_Coordonnate2Rank (i, j, nbcolumns)].xltype = xltypeNum;
					pxArray[XL_Coordonnate2Rank (i, j, nbcolumns)].val.num = C_result.getArray(i*nbcolumns+j); 
				}
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_ARM_GetFwdRatesMatrix" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


__declspec(dllexport) LPXLOPER WINAPI Local_ARM_GetInfo(LPXLOPER XL_secId,
														LPXLOPER XL_type)
{
	ADD_LOG("Local_ARM_GetInfo");
	// return
	static XLOPER XL_result;
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();

	// C variable
	CCString C_secId;
	CCString C_type;

	// error
	static int error;
	static char* reason = "";

	XL_readStrCell(XL_secId,C_secId," ARM_ERR: security id: object expected",C_result);
	XL_readStrCell(XL_type,C_type," ARM_ERR: type: string expected",C_result);

	C_type.toUpper ();

	long retCode;
	long objId;
	CCString prevClass;
	CCString curClass;

	if ( (C_type == "LEG1") || (C_type == "LEG2") || (C_type == "PAYLEG") )
		curClass = LOCAL_SWAPLEG_CLASS;
	else
		curClass = LOCAL_REFVAL_CLASS;

	CCString secClass = LocalGetStringObjectClass (C_secId);
	CCString stringId = GetLastCurCellEnvValue ();

	if(!stringId)
	{
		retCode = ARMLOCAL_ARM_GetInfo (LocalGetNumObjectId (C_secId), secClass, C_type, C_result);

		if (retCode == ARM_OK)
		{
			if (C_result.getLong() != -9999)
			{
				objId = C_result.getLong();

				LocalSetCurCellEnvValue (curClass, objId);

				stringId = LocalMakeObjectId (objId, curClass);
			}
		}
	}
	else
	{
		prevClass = LocalGetStringObjectClass (stringId);

		objId = LocalGetNumObjectId (stringId);
			
		if(curClass == prevClass)
		{
			retCode = ARMLOCAL_ARM_GetInfo (LocalGetNumObjectId (C_secId), secClass, C_type, C_result, objId);

			if (retCode == ARM_OK)
			{
				if (C_result.getLong() != -9999)
				{
					LocalSetCurCellEnvValue (curClass, objId); 
					
					stringId = LocalMakeObjectId (objId, curClass);
				}
			}
		}
		else
		{
			FreeCurCellContent ();

			retCode = ARMLOCAL_ARM_GetInfo (LocalGetNumObjectId (C_secId), secClass, C_type, C_result);
			
			if(retCode == ARM_OK)
			{
				if (C_result.getLong() != -9999)
				{
					objId = C_result.getLong ();
				
					LocalSetCurCellEnvValue (curClass, objId); 

					stringId = LocalMakeObjectId (objId, curClass);
				}
			}
		}
	}

	if(retCode == ARM_OK)
	{
		FreeCurCellErr ();
		if (C_result.getLong() != -9999)
		{
			XL_result.xltype = xltypeStr;
			XL_result.val.str = XL_StrC2StrPascal (stringId);
			XL_result.xltype |= xlbitDLLFree;
		}
		else if (C_result.getDouble() != -9999.0)
		{
			XL_result.xltype = xltypeNum;
			XL_result.val.num = C_result.getDouble();
		}
		else
		{
			XL_result.xltype = xltypeStr;
			XL_result.val.str = XL_StrC2StrPascal (C_result.getString());
			XL_result.xltype |= xlbitDLLFree;
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_ARM_GetInfo" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


__declspec(dllexport) LPXLOPER WINAPI Local_PXL_ARM_GetInfo(LPXLOPER XL_secId,
															LPXLOPER XL_type)
{
	ADD_LOG("Local_PXL_ARM_GetInfo");
	// return
	static XLOPER XL_result;
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();

	// C variable
	CCString C_secId;
	CCString C_type;

	// error
	static int error;
	static char* reason = "";

	XL_readStrCell(XL_secId,C_secId," ARM_ERR: security id: object expected",C_result);
	XL_readStrCell(XL_type,C_type," ARM_ERR: type: string expected",C_result);

	C_type.toUpper ();

	long retCode;
	long objId;
	CCString prevClass;

	CCString curClass;

	if ( (C_type == "LEG1") || (C_type == "LEG2") || (C_type == "PAYLEG") )
		curClass = LOCAL_SWAPLEG_CLASS;
	else
		curClass = LOCAL_REFVAL_CLASS;

	CCString secClass = LocalGetStringObjectClass (C_secId);
	CCString stringId;

	retCode = ARMLOCAL_ARM_GetInfo (LocalGetNumObjectId (C_secId), secClass, C_type, C_result);

	if (retCode == ARM_OK)
	{
		if (C_result.getLong() != -9999)
		{
			objId = C_result.getLong();

			stringId = LocalMakeObjectId (objId, curClass);
		}
	}

	if(retCode == ARM_OK)
	{
		FreeCurCellErr ();
		if (C_result.getLong() != -9999)
		{
			XL_result.xltype = xltypeStr;
			XL_result.val.str = XL_StrC2StrPascal (stringId);
			XL_result.xltype |= xlbitDLLFree;
		}
		else if (C_result.getDouble() != -9999.0)
		{
			XL_result.xltype = xltypeNum;
			XL_result.val.num = C_result.getDouble();
		}
		else
		{
			XL_result.xltype = xltypeStr;
			XL_result.val.str = XL_StrC2StrPascal (C_result.getString());
			XL_result.xltype |= xlbitDLLFree;
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_PXL_ARM_GetInfo" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}

__declspec(dllexport) LPXLOPER WINAPI Local_ARM_GetModelFactorFromSummit(LPXLOPER XL_date,
																		 LPXLOPER XL_model,
																		 LPXLOPER XL_type,
																		 LPXLOPER XL_factorName,
																		 LPXLOPER XL_ccy,
																		 LPXLOPER XL_index,
																		 LPXLOPER XL_cvName,
																		 LPXLOPER XL_calcMethod)
{
	ADD_LOG("Local_ARM_GetModelFactorFromSummit");
	static XLOPER XL_result;

	ARM_result C_result;

	ARM_XL_TRY_BLOCK_BEGIN
	{
		ARM_NOCALCIFWIZ();

		// C variable
		double C_date;
		CCString C_model;
		CCString C_type;
		CCString C_factorName;
		CCString C_ccy;
		CCString C_index;
		CCString C_cvName;
		CCString C_calcMethod;
		long calcModId;

		// error
		static int error;
		static char* reason = "";

		XL_readNumCell(XL_date,C_date," ARM_ERR: date : date expected",C_result);
		XL_readStrCell(XL_model,C_model," ARM_ERR: model: string expected",C_result);
		XL_readStrCell(XL_type,C_type," ARM_ERR: inst type: string expected",C_result);
		XL_readStrCell(XL_factorName,C_factorName," ARM_ERR: factor name: string expected",C_result);
		XL_readStrCell(XL_ccy,C_ccy," ARM_ERR: currency: string expected",C_result);
		XL_readStrCell(XL_index,C_index," ARM_ERR: index1: string expected",C_result);
		XL_readStrCell(XL_cvName,C_cvName," ARM_ERR: CV name: string expected",C_result);
		XL_readStrCellWD(XL_calcMethod,C_calcMethod,"LIN"," ARM_ERR: calculation method: string expected",C_result);

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
			retCode = ARMLOCAL_ARM_GetModelFactor ( C_date,
													C_model,
													C_type,
													C_factorName,
													C_ccy,
													C_index,
													C_cvName,
													calcModId,
													C_result );

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
				retCode = ARMLOCAL_ARM_GetModelFactor ( C_date,
														C_model,
														C_type,
														C_factorName,
														C_ccy,
														C_index,
														C_cvName,
														calcModId,
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
				retCode = ARMLOCAL_ARM_GetModelFactor ( C_date,
														C_model,
														C_type,
														C_factorName,
														C_ccy,
														C_index,
														C_cvName,
														calcModId,
														C_result );
			
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

	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_ARM_GetModelFactorFromSummit" )

	return (LPXLOPER)&XL_result;
}
