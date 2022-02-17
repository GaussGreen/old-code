#pragma warning(disable : 4786)
#pragma warning(disable : 4005)

#include <libCCxll\CCxll.h>

#include "ARM_local_interglob.h"
#include "ARM_local_interface.h"

#include "ARM_xl_trycatch_local.h"
#include <util\fromto.h>

#include <ARM\libarm_local\ARM_local_assetswap.h>
#include <ARM\libarm_local\ARM_local_glob.h>

#include "util\tech_macro.h"

__declspec(dllexport) LPXLOPER WINAPI Local_BasisSwap (LPXLOPER XL_asOfDate,
													   LPXLOPER XL_delivery,
													   LPXLOPER XL_maturity,
													   LPXLOPER XL_margin1,
													   LPXLOPER XL_ccy1,
													   LPXLOPER XL_index1,
													   LPXLOPER XL_forwardCurve1,
													   LPXLOPER XL_discountCurve1,
													   LPXLOPER XL_ccy2,
													   LPXLOPER XL_index2,
													   LPXLOPER XL_forwardCurve2,
													   LPXLOPER XL_discountCurve2,
													   LPXLOPER XL_outMode,
													   LPXLOPER XL_solve,
													   LPXLOPER XL_amortizationId)
{
	ADD_LOG("Local_BasisSwap ");
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
	double C_asOfDate;
	double C_delivery;
	double C_maturity;

	CCString C_ccy1;

	CCString C_index1;
	long liborType1Id;

	CCString C_discountCurve1;
	CCString C_forwardCurve1;
	
	double C_margin1;
	
	CCString C_ccy2;

	CCString C_index2;
	long liborType2Id;

	CCString C_discountCurve2;
	CCString C_forwardCurve2;

	CCString C_amortizationId;

	double C_solve;
	double C_solve_default = BasisSwap_Method_Alg;

	double C_outMode;
	double C_outMode_default = 0.0;

	// error
	static int error;
	static char* reason = "";

	XL_readNumCell(XL_asOfDate,C_asOfDate," ARM_ERR: as of date: date expected",C_result);
	XL_readNumCell(XL_delivery,C_delivery," ARM_ERR: delivery: date expected",C_result);
	XL_readNumCell(XL_maturity,C_maturity," ARM_ERR: maturity: date expected",C_result);
	XL_readStrCell(XL_ccy1,C_ccy1," ARM_ERR: currency 1: string expected",C_result);
	XL_readStrCell(XL_index1,C_index1," ARM_ERR: index 1: string expected",C_result);
	XL_readStrCell(XL_discountCurve1,C_discountCurve1," ARM_ERR: discount curve 1: string expected",C_result);
	XL_readStrCell(XL_forwardCurve1,C_forwardCurve1," ARM_ERR: forward curve 1: string expected",C_result);
	XL_readNumCell(XL_margin1,C_margin1," ARM_ERR: margin 1: numeric expected",C_result);
	XL_readStrCell(XL_ccy2,C_ccy2," ARM_ERR: currency 2: string expected",C_result);
	XL_readStrCell(XL_index2,C_index2," ARM_ERR: index 2: string expected",C_result);
	XL_readStrCell(XL_discountCurve2,C_discountCurve2," ARM_ERR: discount curve 2: string expected",C_result);
	XL_readStrCell(XL_forwardCurve2,C_forwardCurve2," ARM_ERR: forward curve 2: string expected",C_result);
	XL_readStrCellWD(XL_amortizationId,C_amortizationId,"DEFAULT"," ARM_ERR: amortization: object expected",C_result);
	XL_readNumCellWD(XL_outMode,C_outMode,C_outMode_default," ARM_ERR: out mode: numeric expected",C_result);
	XL_readNumCellWD(XL_solve,C_solve,C_solve_default," ARM_ERR: solve: numeric expected",C_result);
	
	if((liborType1Id = ARM_ConvIrIndName (C_index1, C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	if((liborType2Id = ARM_ConvIrIndName (C_index2, C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}
	
	long retCode = ARMLOCAL_BasisSwap (C_asOfDate,
								  C_delivery,
								  C_maturity,
								  C_margin1,
								  C_ccy1,
								  liborType1Id,
								  LocalGetNumObjectId (C_forwardCurve1),
								  LocalGetNumObjectId (C_discountCurve1),
								  C_ccy2,
								  liborType2Id,
								  LocalGetNumObjectId (C_forwardCurve2),
								  LocalGetNumObjectId (C_discountCurve2),
								  LocalGetNumObjectId (C_amortizationId),
								  (long)C_solve,
								  C_result);

	if(retCode == ARM_OK)
	{
		if((C_outMode != 0.0) && ((C_solve == BasisSwap_Method_Num) || (C_solve == BasisSwap_Method_Num_Spr)))
		{
			int nbrows = 9;
			int nbcolumns = 2;
		
			FreeCurCellErr ();
			XL_result.xltype = xltypeMulti;
			XL_result.val.array.columns = nbcolumns;
			XL_result.val.array.rows = nbrows; 
			XL_result.val.array.lparray = pxArray = (LPXLOPER)GlobalAlloc (GMEM_ZEROINIT, nbrows * nbcolumns * sizeof (XLOPER));

			pxArray[XL_Coordonnate2Rank (0, 0, nbcolumns)].xltype = xltypeNum;
			pxArray[XL_Coordonnate2Rank (0, 0, nbcolumns)].val.num = C_result.getDouble (); 
			pxArray[XL_Coordonnate2Rank (0, 1, nbcolumns)].xltype = xltypeStr;
			pxArray[XL_Coordonnate2Rank (0, 1, nbcolumns)].val.str = XL_StrC2StrPascal ("");
			pxArray[XL_Coordonnate2Rank (0, 1, nbcolumns)].xltype |= xlbitDLLFree;
			
			pxArray[XL_Coordonnate2Rank (1, 0, nbcolumns)].xltype = xltypeStr;
			pxArray[XL_Coordonnate2Rank (1, 0, nbcolumns)].val.str = XL_StrC2StrPascal ("StartXNL1");
			pxArray[XL_Coordonnate2Rank (1, 0, nbcolumns)].xltype |= xlbitDLLFree;
			pxArray[XL_Coordonnate2Rank (1, 1, nbcolumns)].xltype = xltypeNum;
			pxArray[XL_Coordonnate2Rank (1, 1, nbcolumns)].val.num = ASP_getStartXNL1 (); 

			pxArray[XL_Coordonnate2Rank (2, 0, nbcolumns)].xltype = xltypeStr;
			pxArray[XL_Coordonnate2Rank (2, 0, nbcolumns)].val.str = XL_StrC2StrPascal ("EndXNL1");
			pxArray[XL_Coordonnate2Rank (2, 0, nbcolumns)].xltype |= xlbitDLLFree;
			pxArray[XL_Coordonnate2Rank (2, 1, nbcolumns)].xltype = xltypeNum;
			pxArray[XL_Coordonnate2Rank (2, 1, nbcolumns)].val.num = ASP_getEndXNL1 (); 

			pxArray[XL_Coordonnate2Rank (3, 0, nbcolumns)].xltype = xltypeStr;
			pxArray[XL_Coordonnate2Rank (3, 0, nbcolumns)].val.str = XL_StrC2StrPascal ("StartXNL2");
			pxArray[XL_Coordonnate2Rank (3, 0, nbcolumns)].xltype |= xlbitDLLFree;
			pxArray[XL_Coordonnate2Rank (3, 1, nbcolumns)].xltype = xltypeNum;
			pxArray[XL_Coordonnate2Rank (3, 1, nbcolumns)].val.num = ASP_getStartXNL2 (); 

			pxArray[XL_Coordonnate2Rank (4, 0, nbcolumns)].xltype = xltypeStr;
			pxArray[XL_Coordonnate2Rank (4, 0, nbcolumns)].val.str = XL_StrC2StrPascal ("EndXNL2");
			pxArray[XL_Coordonnate2Rank (4, 0, nbcolumns)].xltype |= xlbitDLLFree;
			pxArray[XL_Coordonnate2Rank (4, 1, nbcolumns)].xltype = xltypeNum;
			pxArray[XL_Coordonnate2Rank (4, 1, nbcolumns)].val.num = ASP_getEndXNL2 (); 

			pxArray[XL_Coordonnate2Rank (5, 0, nbcolumns)].xltype = xltypeStr;
			pxArray[XL_Coordonnate2Rank (5, 0, nbcolumns)].val.str = XL_StrC2StrPascal ("Price1");
			pxArray[XL_Coordonnate2Rank (5, 0, nbcolumns)].xltype |= xlbitDLLFree;
			pxArray[XL_Coordonnate2Rank (5, 1, nbcolumns)].xltype = xltypeNum;
			pxArray[XL_Coordonnate2Rank (5, 1, nbcolumns)].val.num = ASP_getPrice1 (); 

			pxArray[XL_Coordonnate2Rank (6, 0, nbcolumns)].xltype = xltypeStr;
			pxArray[XL_Coordonnate2Rank (6, 0, nbcolumns)].val.str = XL_StrC2StrPascal ("Price2");
			pxArray[XL_Coordonnate2Rank (6, 0, nbcolumns)].xltype |= xlbitDLLFree;
			pxArray[XL_Coordonnate2Rank (6, 1, nbcolumns)].xltype = xltypeNum;
			pxArray[XL_Coordonnate2Rank (6, 1, nbcolumns)].val.num = ASP_getPrice2 ();

			pxArray[XL_Coordonnate2Rank (7, 0, nbcolumns)].xltype = xltypeStr;
			pxArray[XL_Coordonnate2Rank (7, 0, nbcolumns)].val.str = XL_StrC2StrPascal ("Minimum");
			pxArray[XL_Coordonnate2Rank (7, 0, nbcolumns)].xltype |= xlbitDLLFree;
			pxArray[XL_Coordonnate2Rank (7, 1, nbcolumns)].xltype = xltypeNum;
			pxArray[XL_Coordonnate2Rank (7, 1, nbcolumns)].val.num = ASP_getMinimum ();

			pxArray[XL_Coordonnate2Rank (8, 0, nbcolumns)].xltype = xltypeStr;
			pxArray[XL_Coordonnate2Rank (8, 0, nbcolumns)].val.str = XL_StrC2StrPascal ("Iterations");
			pxArray[XL_Coordonnate2Rank (8, 0, nbcolumns)].xltype |= xlbitDLLFree;
			pxArray[XL_Coordonnate2Rank (8, 1, nbcolumns)].xltype = xltypeNum;
			pxArray[XL_Coordonnate2Rank (8, 1, nbcolumns)].val.num = ASP_getNbIter ();
		}
		else
		{
			FreeCurCellErr ();
			XL_result.xltype = xltypeNum;
			XL_result.val.num = C_result.getDouble ();
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_BasisSwap" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


__declspec(dllexport) LPXLOPER WINAPI Local_ASWMargin (LPXLOPER XL_bondMaturity,
													   LPXLOPER XL_bondCoupon,
													   LPXLOPER XL_bondFrequency,
													   LPXLOPER XL_bondBase,
													   LPXLOPER XL_bondPrice,
													   LPXLOPER XL_bondRedemptionPrice,
													   LPXLOPER XL_asOfDate,
													   LPXLOPER XL_delivery,
													   LPXLOPER XL_fixDecompFrequency,
													   LPXLOPER XL_ccy1,
													   LPXLOPER XL_index1,
													   LPXLOPER XL_forwardCurve1,
													   LPXLOPER XL_discountCurve1,
													   LPXLOPER XL_ccy2,
													   LPXLOPER XL_index2,
													   LPXLOPER XL_forwardCurve2,
													   LPXLOPER XL_discountCurve2,
													   LPXLOPER XL_solve,
													   LPXLOPER XL_amortizationId)
{
	ADD_LOG("Local_ASWMargin ");
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;
// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();

	// C variable
	double C_bondMaturity;

    CCString C_bondcpn_str;
    double C_bondcpn_double;
	long   cpnType;

	CCString C_bondFrequency;
	long bondFrequencyId;

	CCString C_bondBase;
	long bondBaseId;

	double C_bondPrice;

	double C_bondRedemptionPrice;
	double C_bondRedemptionPrice_default = 100;

	double C_asOfDate;
	double C_delivery;

	CCString C_fixDecompFrequency;
	long fixDecompFrequencyId;
	
	CCString C_ccy1;

	CCString C_index1;
	long liborType1Id;

	CCString C_discountCurve1;
	CCString C_forwardCurve1;
	
	CCString C_ccy2;

	CCString C_index2;
	long liborType2Id;

	CCString C_discountCurve2 ("");
	CCString C_forwardCurve2 ("");

	double C_viewFlag;
	double C_viewFlag_default = 0;

	CCString C_amortizationId;

	double C_solve;
	double C_solve_default = BasisSwap_Method_Alg;
	
	double C_minValue = -1000.0;
	double C_maxValue = 1000.0;

	// error
	static int error;
	static char* reason = "";

	XL_readNumCell(XL_bondMaturity,C_bondMaturity," ARM_ERR: bond maturity: date expected",C_result);
	XL_readStrOrNumCell(XL_bondCoupon, C_bondcpn_str, C_bondcpn_double, cpnType,
		   " ARM_ERR: cpn: numeric or object ID string expected",C_result);
	XL_readStrCell(XL_bondFrequency,C_bondFrequency," ARM_ERR: bond frequency: string expected",C_result);
	XL_readStrCell(XL_bondBase,C_bondBase," ARM_ERR: bond base: string expected",C_result);
	XL_readNumCell(XL_bondPrice,C_bondPrice," ARM_ERR: bond price: numeric expected",C_result);
	XL_readNumCellWD(XL_bondRedemptionPrice,C_bondRedemptionPrice,C_bondRedemptionPrice_default," ARM_ERR: bond redemption price: numeric expected",C_result);
	XL_readNumCell(XL_asOfDate,C_asOfDate," ARM_ERR: as of date: date expected",C_result);
	XL_readNumCell(XL_delivery,C_delivery," ARM_ERR: delivery: date expected",C_result);
	XL_readStrCellWD(XL_fixDecompFrequency,C_fixDecompFrequency,"1"," ARM_ERR: fix decomp frequency: string expected",C_result);
	XL_readStrCell(XL_ccy1,C_ccy1," ARM_ERR: currency 1: string expected",C_result);
	XL_readStrCell(XL_index1,C_index1," ARM_ERR: index 1: string expected",C_result);
	XL_readStrCell(XL_forwardCurve1,C_forwardCurve1," ARM_ERR: forward curve 1: string expected",C_result);
	XL_readStrCellWD(XL_ccy2,C_ccy2,"NONE"," ARM_ERR: currency 2: string expected",C_result);
	XL_readStrCellWD(XL_amortizationId,C_amortizationId,"DEFAULT"," ARM_ERR: amortization: object expected",C_result);
	XL_readNumCellWD(XL_solve,C_solve,C_solve_default," ARM_ERR: solve: numeric expected",C_result);
	
	// Modif 21/06/2001 : prise en compte de l'index2 même si la devise 2 n'est pas définie
	XL_readStrCellWD(XL_index2,C_index2,C_index1," ARM_ERR: index 2: string expected",C_result);
	/*************/

	C_viewFlag = C_viewFlag_default;

	long floatResetFreqId;
	long floatPayFreqId;
	long discountCurve1Id = ARM_NULL_OBJECT_ID; 

	if(!(C_ccy2 == "NONE") && (C_ccy1 != C_ccy2))
	{
		XL_readStrCell(XL_index2,C_index2," ARM_ERR: index 2: string expected",C_result);
		XL_readStrCell(XL_discountCurve2,C_discountCurve2," ARM_ERR: discount curve 2: string expected",C_result);
		XL_readStrCell(XL_forwardCurve2,C_forwardCurve2," ARM_ERR: forward curve 2: string expected",C_result);
		XL_readStrCell(XL_discountCurve1,C_discountCurve1," ARM_ERR: discount curve 1: string expected",C_result);
		if((liborType2Id = ARM_ConvIrIndName (C_index2, C_result)) == ARM_DEFAULT_ERR)
		{
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}
		discountCurve1Id = LocalGetNumObjectId (C_discountCurve1);
	}
	else
	{
		if((liborType2Id = ARM_ConvIrIndName (C_index1, C_result)) == ARM_DEFAULT_ERR)
		{
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}
	}

	if((liborType1Id = ARM_ConvIrIndName (C_index1, C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	if((fixDecompFrequencyId = ARM_ConvFrequency (C_fixDecompFrequency, C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	if((bondFrequencyId = ARM_ConvFrequency (C_bondFrequency, C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	bondBaseId = ARM_ConvDayCount (C_bondBase);

	if((floatResetFreqId = ARM_ConvIrIndNameToFreq (C_index1, C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	if((floatPayFreqId = ARM_ConvIrIndNameToFreq (C_index2, C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

    if ( cpnType == XL_TYPE_STRING )
	{
	   C_bondcpn_double = (double) LocalGetNumObjectId(C_bondcpn_str);

	   cpnType = 1L;
	}
	else
	{
	   cpnType = 0L;
	}

	CCString C_sockId ("");

/*	if(C_viewFlag == 1)
	{
		C_sockId = getSockId ();
	}
*/	
	long retCode = ARMLOCAL_ASWMarginNew (C_bondMaturity,
										  cpnType,
										  C_bondcpn_double,
										  bondFrequencyId,
										  bondBaseId,
										  C_bondPrice,
										  C_bondRedemptionPrice,
										  C_asOfDate,
										  C_delivery,
										  fixDecompFrequencyId,
										  floatResetFreqId,
										  floatPayFreqId,
										  C_ccy1,
										  liborType1Id,
										  LocalGetNumObjectId (C_forwardCurve1),
										  discountCurve1Id,
										  C_ccy2,
										  liborType2Id,
										  LocalGetNumObjectId (C_forwardCurve2),
										  LocalGetNumObjectId (C_discountCurve2),
										  LocalGetNumObjectId (C_amortizationId),
										  (long)C_solve,
										  C_minValue,
										  C_maxValue,
										  C_result);

/*	if(C_viewFlag == 1)
	{
		getArmViewFile (C_sockId);
	}
*/												
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_ASWMargin" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}



__declspec(dllexport) LPXLOPER WINAPI Local_BondASWMargin (LPXLOPER XL_bond,
														   LPXLOPER XL_bondPrice,
														   LPXLOPER XL_asOfDate,
														   LPXLOPER XL_delivery,
														   LPXLOPER XL_ccy1,
														   LPXLOPER XL_ccy2,
														   LPXLOPER XL_index2,
														   LPXLOPER XL_forwardCurve1,
														   LPXLOPER XL_discountCurve1,
														   LPXLOPER XL_forwardCurve2,
														   LPXLOPER XL_discountCurve2,
														   LPXLOPER XL_solve,
														   LPXLOPER XL_amortizationId)
{
	ADD_LOG("Local_BondASWMargin ");
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();

	// C variable
	CCString C_bond;

	double C_bondPrice;

	double C_asOfDate;
	double C_delivery;
	
	CCString C_ccy1;

	CCString C_discountCurve1;
	CCString C_forwardCurve1;
	
	CCString C_ccy2;

	CCString C_index2;
	long liborType2Id;

	CCString C_discountCurve2 ("");
	CCString C_forwardCurve2 ("");

	double C_viewFlag;
	double C_viewFlag_default = 0;

	CCString C_amortizationId;

	double C_solve;
	double C_solve_default = BasisSwap_Method_Alg;
	
	double C_minValue = 0.0;
	double C_maxValue = 1000.0;

	// error
	static int error;
	static char* reason = "";

	XL_readStrCell(XL_bond,C_bond," ARM_ERR: bond object: string expected",C_result);
	XL_readNumCell(XL_bondPrice,C_bondPrice," ARM_ERR: bond price: numeric expected",C_result);
	XL_readNumCell(XL_asOfDate,C_asOfDate," ARM_ERR: as of date: date expected",C_result);
	XL_readNumCell(XL_delivery,C_delivery," ARM_ERR: delivery: date expected",C_result);
	XL_readStrCell(XL_ccy1,C_ccy1," ARM_ERR: currency 1: string expected",C_result);
	XL_readStrCell(XL_forwardCurve1,C_forwardCurve1," ARM_ERR: forward curve 1: string expected",C_result);
	XL_readStrCellWD(XL_ccy2,C_ccy2,C_ccy1," ARM_ERR: currency 2: string expected",C_result);
	XL_readStrCellWD(XL_amortizationId,C_amortizationId,"DEFAULT"," ARM_ERR: amortization: object expected",C_result);
	XL_readNumCellWD(XL_solve,C_solve,C_solve_default," ARM_ERR: solve: numeric expected",C_result);
	
	// Modif 21/06/2001 : prise en compte de l'index2 même si la devise 2 n'est pas définie
	XL_readStrCell(XL_index2,C_index2," ARM_ERR: index 2: string expected",C_result);
	/*************/

	C_viewFlag = C_viewFlag_default;

	long floatResetFreqId;
	long floatPayFreqId;
	long discountCurve1Id = ARM_NULL_OBJECT_ID; 

	if(C_ccy1 != C_ccy2)
	{
		XL_readStrCell(XL_discountCurve2,C_discountCurve2," ARM_ERR: discount curve 2: string expected",C_result);
		XL_readStrCell(XL_forwardCurve2,C_forwardCurve2," ARM_ERR: forward curve 2: string expected",C_result);
		XL_readStrCell(XL_discountCurve1,C_discountCurve1," ARM_ERR: discount curve 1: string expected",C_result);
		discountCurve1Id = LocalGetNumObjectId (C_discountCurve1);
	}
	else
	{
	}

	if((liborType2Id = ARM_ConvIrIndName (C_index2, C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	if((floatResetFreqId = ARM_ConvIrIndNameToFreq (C_index2, C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	if((floatPayFreqId = ARM_ConvIrIndNameToFreq (C_index2, C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}
	
	long retCode = ARMLOCAL_BondASWMargin (LocalGetNumObjectId(C_bond),
										   C_bondPrice,
										   C_asOfDate,
										   C_delivery,
										   floatResetFreqId,
										   floatPayFreqId,
										   C_ccy1,
										   LocalGetNumObjectId (C_forwardCurve1),
										   discountCurve1Id,
										   C_ccy2,
										   liborType2Id,
										   LocalGetNumObjectId (C_forwardCurve2),
										   LocalGetNumObjectId (C_discountCurve2),
										   LocalGetNumObjectId (C_amortizationId),
										   (long)C_solve,
										   C_minValue,
										   C_maxValue,
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_BondASWMargin" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}



__declspec(dllexport) LPXLOPER WINAPI Local_ASWPrice (LPXLOPER XL_bondMaturity,
													  LPXLOPER XL_bondCoupon,
													  LPXLOPER XL_bondFrequency,
													  LPXLOPER XL_bondBase,
													  LPXLOPER XL_bondMargin,
													  LPXLOPER XL_bondRedemptionPrice,
													  LPXLOPER XL_asOfDate,
													  LPXLOPER XL_delivery,
													  LPXLOPER XL_fixDecompFrequency,
													  LPXLOPER XL_ccy1,
													  LPXLOPER XL_index1,
													  LPXLOPER XL_forwardCurve1,
													  LPXLOPER XL_discountCurve1,
													  LPXLOPER XL_ccy2,
													  LPXLOPER XL_index2,
													  LPXLOPER XL_forwardCurve2,
													  LPXLOPER XL_discountCurve2,
													  LPXLOPER XL_solve,
													  LPXLOPER XL_amortizationId)
{
	ADD_LOG("Local_ASWPrice ");
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();

	// C variable
	double C_bondMaturity;

    CCString C_bondcpn_str;
    double C_bondcpn_double;
	long   cpnType;

	CCString C_bondFrequency;
	long bondFrequencyId;

	CCString C_bondBase;
	long bondBaseId;

	double C_bondMargin;

	double C_bondRedemptionPrice;
	double C_bondRedemptionPrice_default = 100.0;

	double C_asOfDate;
	double C_delivery;

	CCString C_fixDecompFrequency;
	long fixDecompFrequencyId;
	
	CCString C_ccy1;

	CCString C_index1;
	long liborType1Id;

	CCString C_discountCurve1;
	CCString C_forwardCurve1;
	
	CCString C_ccy2;

	CCString C_index2;
	long liborType2Id;

	CCString C_discountCurve2 ("");
	CCString C_forwardCurve2 ("");

	CCString C_amortizationId;

	double C_solve;
	double C_solve_default = BasisSwap_Method_Alg;

	double C_minValue = -1000.0;
	double C_maxValue = 2000.0;
	
	// error
	static int error;
	static char* reason = "";

	XL_readNumCell(XL_bondMaturity,C_bondMaturity," ARM_ERR: bond maturity: date expected",C_result);
	XL_readStrOrNumCell(XL_bondCoupon, C_bondcpn_str, C_bondcpn_double, cpnType,
		   " ARM_ERR: cpn: numeric or object ID string expected",C_result);
	XL_readStrCell(XL_bondFrequency,C_bondFrequency," ARM_ERR: bond frequency: string expected",C_result);
	XL_readStrCell(XL_bondBase,C_bondBase," ARM_ERR: bond base: string expected",C_result);
	XL_readNumCell(XL_bondMargin,C_bondMargin," ARM_ERR: bond margin: numeric expected",C_result);
	XL_readNumCellWD(XL_bondRedemptionPrice,C_bondRedemptionPrice,C_bondRedemptionPrice_default," ARM_ERR: bond redemption price: numeric expected",C_result);
	XL_readNumCell(XL_asOfDate,C_asOfDate," ARM_ERR: as of date: date expected",C_result);
	XL_readNumCell(XL_delivery,C_delivery," ARM_ERR: delivery: date expected",C_result);
	XL_readStrCellWD(XL_fixDecompFrequency,C_fixDecompFrequency,"1"," ARM_ERR: fix decomp frequency: string expected",C_result);
	XL_readStrCell(XL_ccy1,C_ccy1," ARM_ERR: currency 1: string expected",C_result);
	XL_readStrCell(XL_index1,C_index1," ARM_ERR: index 1: string expected",C_result);
	XL_readStrCell(XL_forwardCurve1,C_forwardCurve1," ARM_ERR: forward curve 1: string expected",C_result);
	XL_readStrCellWD(XL_ccy2,C_ccy2,"NONE"," ARM_ERR: currency 2: string expected",C_result);
	XL_readStrCellWD(XL_amortizationId,C_amortizationId,"DEFAULT"," ARM_ERR: amortization: object expected",C_result);

	// Modif 21/06/2001 : prise en compte de l'index2 même si la devise 2 n'est pas définie
	XL_readStrCellWD(XL_index2,C_index2,C_index1," ARM_ERR: index 2: string expected",C_result);
	/*************/
	
	XL_readNumCellWD(XL_solve,C_solve,C_solve_default," ARM_ERR: solve: numeric expected",C_result);

	long floatResetFreqId;
	long floatPayFreqId;
	long discountCurve1Id = ARM_NULL_OBJECT_ID;

	if(!(C_ccy2 == "NONE") && (C_ccy1 != C_ccy2))
	{
		XL_readStrCell(XL_index2,C_index2," ARM_ERR: index 2: string expected",C_result);
		XL_readStrCell(XL_discountCurve2,C_discountCurve2," ARM_ERR: discount curve 2: string expected",C_result);
		XL_readStrCell(XL_forwardCurve2,C_forwardCurve2," ARM_ERR: forward curve 2: string expected",C_result);
		XL_readStrCell(XL_discountCurve1,C_discountCurve1," ARM_ERR: discount curve 1: string expected",C_result);
		if((liborType2Id = ARM_ConvIrIndName (C_index2, C_result)) == ARM_DEFAULT_ERR)
		{
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}
		discountCurve1Id = LocalGetNumObjectId (C_discountCurve1);
	}
	else
	{		
		if((liborType2Id = ARM_ConvIrIndName (C_index1, C_result)) == ARM_DEFAULT_ERR)
		{
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}
	}

	if((liborType1Id = ARM_ConvIrIndName (C_index1, C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	if((fixDecompFrequencyId = ARM_ConvFrequency (C_fixDecompFrequency, C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	if((bondFrequencyId = ARM_ConvFrequency (C_bondFrequency, C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	bondBaseId = ARM_ConvDayCount (C_bondBase);

	if((floatResetFreqId = ARM_ConvIrIndNameToFreq (C_index1, C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	if((floatPayFreqId = ARM_ConvIrIndNameToFreq (C_index2, C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	long retCode = ARMLOCAL_ASWPriceNew (C_bondMaturity,
										 cpnType,
										 C_bondcpn_double,
										 bondFrequencyId,
										 bondBaseId,
										 C_bondMargin,
										 C_bondRedemptionPrice,
										 C_asOfDate,
										 C_delivery,
										 fixDecompFrequencyId,
										 floatResetFreqId,
										 floatPayFreqId,
										 C_ccy1,
										 liborType1Id,
										 LocalGetNumObjectId (C_forwardCurve1),
										 discountCurve1Id,
										 C_ccy2,
										 liborType2Id,
										 LocalGetNumObjectId (C_forwardCurve2),
										 LocalGetNumObjectId (C_discountCurve2),
										 LocalGetNumObjectId (C_amortizationId),
										 (long)C_solve,
										 C_minValue,
										 C_maxValue,
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_ASWPrice" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


__declspec(dllexport) LPXLOPER WINAPI Local_BondASWPrice (LPXLOPER XL_bond,
														  LPXLOPER XL_bondMargin,
														  LPXLOPER XL_asOfDate,
														  LPXLOPER XL_delivery,
														  LPXLOPER XL_ccy1,
														  LPXLOPER XL_ccy2,
														  LPXLOPER XL_index2,
														  LPXLOPER XL_forwardCurve1,
														  LPXLOPER XL_discountCurve1,
														  LPXLOPER XL_forwardCurve2,
														  LPXLOPER XL_discountCurve2,
														  LPXLOPER XL_solve,
														  LPXLOPER XL_amortizationId)
{
	ADD_LOG("Local_BondASWPrice ");
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();

	// C variable
	CCString C_bond;

	double C_bondMargin;

	double C_asOfDate;
	double C_delivery;

	CCString C_ccy1;

	CCString C_discountCurve1;
	CCString C_forwardCurve1;
	
	CCString C_ccy2;

	CCString C_index2;
	long liborType2Id;

	CCString C_discountCurve2 ("");
	CCString C_forwardCurve2 ("");

	CCString C_amortizationId;

	double C_solve;
	double C_solve_default = BasisSwap_Method_Alg;

	double C_minValue = -1000.0;
	double C_maxValue = 2000.0;
	
	// error
	static int error;
	static char* reason = "";

	XL_readStrCell(XL_bond,C_bond," ARM_ERR: bond object: string expected",C_result);
	XL_readNumCell(XL_bondMargin,C_bondMargin," ARM_ERR: bond margin: numeric expected",C_result);
	XL_readNumCell(XL_asOfDate,C_asOfDate," ARM_ERR: as of date: date expected",C_result);
	XL_readNumCell(XL_delivery,C_delivery," ARM_ERR: delivery: date expected",C_result);
	XL_readStrCell(XL_ccy1,C_ccy1," ARM_ERR: currency 1: string expected",C_result);
	XL_readStrCell(XL_forwardCurve1,C_forwardCurve1," ARM_ERR: forward curve 1: string expected",C_result);
	XL_readStrCellWD(XL_ccy2,C_ccy2,C_ccy1," ARM_ERR: currency 2: string expected",C_result);
	XL_readStrCellWD(XL_amortizationId,C_amortizationId,"DEFAULT"," ARM_ERR: amortization: object expected",C_result);
	XL_readNumCellWD(XL_solve,C_solve,C_solve_default," ARM_ERR: solve: numeric expected",C_result);
	
	// Modif 21/06/2001 : prise en compte de l'index2 même si la devise 2 n'est pas définie
	XL_readStrCell(XL_index2,C_index2," ARM_ERR: index 2: string expected",C_result);
	/*************/

	long floatResetFreqId;
	long floatPayFreqId;
	long discountCurve1Id = ARM_NULL_OBJECT_ID; 

	if(C_ccy1 != C_ccy2)
	{
		XL_readStrCell(XL_discountCurve2,C_discountCurve2," ARM_ERR: discount curve 2: string expected",C_result);
		XL_readStrCell(XL_forwardCurve2,C_forwardCurve2," ARM_ERR: forward curve 2: string expected",C_result);
		XL_readStrCell(XL_discountCurve1,C_discountCurve1," ARM_ERR: discount curve 1: string expected",C_result);
		discountCurve1Id = LocalGetNumObjectId (C_discountCurve1);
	}
	else
	{
	}

	if((liborType2Id = ARM_ConvIrIndName (C_index2, C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	if((floatResetFreqId = ARM_ConvIrIndNameToFreq (C_index2, C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	if((floatPayFreqId = ARM_ConvIrIndNameToFreq (C_index2, C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	long retCode = ARMLOCAL_BondASWPrice (LocalGetNumObjectId(C_bond),
										  C_bondMargin,
										  C_asOfDate,
										  C_delivery,
										  floatResetFreqId,
										  floatPayFreqId,
										  C_ccy1,
										  LocalGetNumObjectId (C_forwardCurve1),
										  discountCurve1Id,
										  C_ccy2,
										  liborType2Id,
										  LocalGetNumObjectId (C_forwardCurve2),
										  LocalGetNumObjectId (C_discountCurve2),
										  LocalGetNumObjectId (C_amortizationId),
										  (long)C_solve,
										  C_minValue,
										  C_maxValue,
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_BondASWPrice" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}





__declspec(dllexport) LPXLOPER WINAPI Local_FRNMargin (LPXLOPER XL_asOfDate,
													   LPXLOPER XL_delivery,
													   LPXLOPER XL_maturity,
													   LPXLOPER XL_ccy1,
													   LPXLOPER XL_index1,
													   LPXLOPER XL_forwardCurve1,
													   LPXLOPER XL_discountCurve1,
													   LPXLOPER XL_facialMargin,
													   LPXLOPER XL_price,
													   LPXLOPER XL_ccy2,
													   LPXLOPER XL_index2,
													   LPXLOPER XL_forwardCurve2,
													   LPXLOPER XL_discountCurve2,
													   LPXLOPER XL_fixing,
													   LPXLOPER XL_spread,
													   LPXLOPER XL_outMode,
													   LPXLOPER XL_solve,
													   LPXLOPER XL_amortizationId)
{
	ADD_LOG("Local_FRNMargin ");
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
	double C_asOfDate;
	double C_delivery;
	double C_maturity;

	CCString C_ccy1;

	CCString C_index1;
	long liborType1Id;

	double C_facialMargin;

	CCString C_ccy2;

	CCString C_index2;
	long liborType2Id = ARM_NULL_OBJECT_ID;

	double C_price;
	
	CCString C_discountCurve1;
	CCString C_forwardCurve1;

	long frequencyId;
	long frequencyId2;

	CCString C_discountCurve2;
	CCString C_forwardCurve2;

	double C_fixing;
	double C_fixing_default = 0.0;

	CCString C_amortizationId;

	double C_spread;
	double C_spread_default = 0.0;

	double C_solve;
	double C_solve_default = BasisSwap_Method_Alg;

	double C_outMode;
	double C_outMode_default = 0.0;

	// error
	static int error;
	static char* reason = "";

	XL_readNumCell(XL_asOfDate,C_asOfDate," ARM_ERR: as of date: date expected",C_result);
	XL_readNumCell(XL_delivery,C_delivery," ARM_ERR: delivery: date expected",C_result);
	XL_readNumCell(XL_maturity,C_maturity," ARM_ERR: maturity: date expected",C_result);
	XL_readStrCell(XL_ccy1,C_ccy1," ARM_ERR: currency 1: string expected",C_result);
	XL_readStrCell(XL_index1,C_index1," ARM_ERR: index 1: string expected",C_result);
	XL_readNumCell(XL_facialMargin,C_facialMargin," ARM_ERR: facial margin: numeric expected",C_result);
	XL_readStrCellWD(XL_ccy2,C_ccy2,"NONE"," ARM_ERR: currency 2: string expected",C_result);
	// Modif du 21/06/2001 : prise en compte du 2eme index en monocurrency
	XL_readStrCellWD(XL_index2,C_index2,C_index1," ARM_ERR: index 2: string expected",C_result);
	XL_readNumCell(XL_price,C_price," ARM_ERR: price: numeric expected",C_result);
	XL_readStrCell(XL_discountCurve1,C_discountCurve1," ARM_ERR: discount curve 1: string expected",C_result);
	XL_readStrCell(XL_forwardCurve1,C_forwardCurve1," ARM_ERR: forward curve 1: string expected",C_result);
	XL_readNumCellWD(XL_fixing,C_fixing,C_fixing_default," ARM_ERR: fixing: numeric expected",C_result);
	XL_readNumCellWD(XL_spread,C_spread,C_spread_default," ARM_ERR: spread: numeric expected",C_result);
	XL_readNumCellWD(XL_outMode,C_outMode,C_outMode_default," ARM_ERR: out mode: numeric expected",C_result);
	XL_readNumCellWD(XL_solve,C_solve,C_solve_default," ARM_ERR: solve: numeric expected",C_result);
	XL_readStrCellWD(XL_amortizationId,C_amortizationId,"DEFAULT"," ARM_ERR: amortization: object expected",C_result);
	
	if(!(C_ccy2 == "NONE"))
	{
		XL_readStrCell(XL_index2,C_index2," ARM_ERR: index 2: string expected",C_result);
		XL_readStrCell(XL_discountCurve2,C_discountCurve2," ARM_ERR: discount curve 2: string expected",C_result);
		XL_readStrCell(XL_forwardCurve2,C_forwardCurve2," ARM_ERR: forward curve 2: string expected",C_result);
		if((liborType2Id = ARM_ConvIrIndName (C_index2, C_result)) == ARM_DEFAULT_ERR)
		{
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}
	}
	
	if((liborType1Id = ARM_ConvIrIndName (C_index1, C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	if((frequencyId = ARM_ConvIrIndNameToFreq (C_index1, C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	if((frequencyId2 = ARM_ConvIrIndNameToFreq (C_index2, C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	long retCode = ARMLOCAL_FRNMarginNew (C_asOfDate,
							      C_delivery,
								  C_maturity,
								  C_ccy1,
								  liborType1Id,
								  LocalGetNumObjectId (C_forwardCurve1),
								  LocalGetNumObjectId (C_discountCurve1),
								  C_facialMargin,
								  C_price,
								  frequencyId,
								  C_ccy2,
								  liborType2Id,
								  LocalGetNumObjectId (C_forwardCurve2),
								  LocalGetNumObjectId (C_discountCurve2),
								  LocalGetNumObjectId (C_amortizationId),
								  frequencyId2,
								  C_fixing,
								  C_spread,
								  (long)C_solve,
								  C_result);
											
	if(retCode == ARM_OK)
	{
		if((C_outMode != 0.0) && ((C_solve == BasisSwap_Method_Num) || (C_solve == BasisSwap_Method_Num_Spr)))
		{
			int nbrows = 9;
			int nbcolumns = 2;
		
			FreeCurCellErr ();
			XL_result.xltype = xltypeMulti;
			XL_result.val.array.columns = nbcolumns;
			XL_result.val.array.rows = nbrows; 
			XL_result.val.array.lparray = pxArray = (LPXLOPER)GlobalAlloc (GMEM_ZEROINIT, nbrows * nbcolumns * sizeof (XLOPER));

			pxArray[XL_Coordonnate2Rank (0, 0, nbcolumns)].xltype = xltypeNum;
			pxArray[XL_Coordonnate2Rank (0, 0, nbcolumns)].val.num = C_result.getDouble (); 
			pxArray[XL_Coordonnate2Rank (0, 1, nbcolumns)].xltype = xltypeStr;
			pxArray[XL_Coordonnate2Rank (0, 1, nbcolumns)].val.str = XL_StrC2StrPascal ("");
			pxArray[XL_Coordonnate2Rank (0, 1, nbcolumns)].xltype |= xlbitDLLFree;
			
			pxArray[XL_Coordonnate2Rank (1, 0, nbcolumns)].xltype = xltypeStr;
			pxArray[XL_Coordonnate2Rank (1, 0, nbcolumns)].val.str = XL_StrC2StrPascal ("StartXNL1");
			pxArray[XL_Coordonnate2Rank (1, 0, nbcolumns)].xltype |= xlbitDLLFree;
			pxArray[XL_Coordonnate2Rank (1, 1, nbcolumns)].xltype = xltypeNum;
			pxArray[XL_Coordonnate2Rank (1, 1, nbcolumns)].val.num = ASP_getStartXNL1 (); 

			pxArray[XL_Coordonnate2Rank (2, 0, nbcolumns)].xltype = xltypeStr;
			pxArray[XL_Coordonnate2Rank (2, 0, nbcolumns)].val.str = XL_StrC2StrPascal ("EndXNL1");
			pxArray[XL_Coordonnate2Rank (2, 0, nbcolumns)].xltype |= xlbitDLLFree;
			pxArray[XL_Coordonnate2Rank (2, 1, nbcolumns)].xltype = xltypeNum;
			pxArray[XL_Coordonnate2Rank (2, 1, nbcolumns)].val.num = ASP_getEndXNL1 (); 

			pxArray[XL_Coordonnate2Rank (3, 0, nbcolumns)].xltype = xltypeStr;
			pxArray[XL_Coordonnate2Rank (3, 0, nbcolumns)].val.str = XL_StrC2StrPascal ("StartXNL2");
			pxArray[XL_Coordonnate2Rank (3, 0, nbcolumns)].xltype |= xlbitDLLFree;
			pxArray[XL_Coordonnate2Rank (3, 1, nbcolumns)].xltype = xltypeNum;
			pxArray[XL_Coordonnate2Rank (3, 1, nbcolumns)].val.num = ASP_getStartXNL2 (); 

			pxArray[XL_Coordonnate2Rank (4, 0, nbcolumns)].xltype = xltypeStr;
			pxArray[XL_Coordonnate2Rank (4, 0, nbcolumns)].val.str = XL_StrC2StrPascal ("EndXNL2");
			pxArray[XL_Coordonnate2Rank (4, 0, nbcolumns)].xltype |= xlbitDLLFree;
			pxArray[XL_Coordonnate2Rank (4, 1, nbcolumns)].xltype = xltypeNum;
			pxArray[XL_Coordonnate2Rank (4, 1, nbcolumns)].val.num = ASP_getEndXNL2 (); 

			pxArray[XL_Coordonnate2Rank (5, 0, nbcolumns)].xltype = xltypeStr;
			pxArray[XL_Coordonnate2Rank (5, 0, nbcolumns)].val.str = XL_StrC2StrPascal ("Price1");
			pxArray[XL_Coordonnate2Rank (5, 0, nbcolumns)].xltype |= xlbitDLLFree;
			pxArray[XL_Coordonnate2Rank (5, 1, nbcolumns)].xltype = xltypeNum;
			pxArray[XL_Coordonnate2Rank (5, 1, nbcolumns)].val.num = ASP_getPrice1 (); 

			pxArray[XL_Coordonnate2Rank (6, 0, nbcolumns)].xltype = xltypeStr;
			pxArray[XL_Coordonnate2Rank (6, 0, nbcolumns)].val.str = XL_StrC2StrPascal ("Price2");
			pxArray[XL_Coordonnate2Rank (6, 0, nbcolumns)].xltype |= xlbitDLLFree;
			pxArray[XL_Coordonnate2Rank (6, 1, nbcolumns)].xltype = xltypeNum;
			pxArray[XL_Coordonnate2Rank (6, 1, nbcolumns)].val.num = ASP_getPrice2 ();

			pxArray[XL_Coordonnate2Rank (7, 0, nbcolumns)].xltype = xltypeStr;
			pxArray[XL_Coordonnate2Rank (7, 0, nbcolumns)].val.str = XL_StrC2StrPascal ("Minimum");
			pxArray[XL_Coordonnate2Rank (7, 0, nbcolumns)].xltype |= xlbitDLLFree;
			pxArray[XL_Coordonnate2Rank (7, 1, nbcolumns)].xltype = xltypeNum;
			pxArray[XL_Coordonnate2Rank (7, 1, nbcolumns)].val.num = ASP_getMinimum ();

			pxArray[XL_Coordonnate2Rank (8, 0, nbcolumns)].xltype = xltypeStr;
			pxArray[XL_Coordonnate2Rank (8, 0, nbcolumns)].val.str = XL_StrC2StrPascal ("Iterations");
			pxArray[XL_Coordonnate2Rank (8, 0, nbcolumns)].xltype |= xlbitDLLFree;
			pxArray[XL_Coordonnate2Rank (8, 1, nbcolumns)].xltype = xltypeNum;
			pxArray[XL_Coordonnate2Rank (8, 1, nbcolumns)].val.num = ASP_getNbIter ();
		}
		else
		{
			FreeCurCellErr ();
			XL_result.xltype = xltypeNum;
			XL_result.val.num = C_result.getDouble ();
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_FRNMargin" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


__declspec(dllexport) LPXLOPER WINAPI Local_FRNPrice (LPXLOPER XL_asOfDate,
													  LPXLOPER XL_delivery,
													  LPXLOPER XL_maturity,
													  LPXLOPER XL_ccy1,
													  LPXLOPER XL_index1,
													  LPXLOPER XL_forwardCurve1,
													  LPXLOPER XL_discountCurve1,
													  LPXLOPER XL_facialMargin,
													  LPXLOPER XL_valoMargin,
													  LPXLOPER XL_ccy2,
													  LPXLOPER XL_index2,
													  LPXLOPER XL_forwardCurve2,
													  LPXLOPER XL_discountCurve2,
													  LPXLOPER XL_fixing,
													  LPXLOPER XL_spread,
													  LPXLOPER XL_outMode,
													  LPXLOPER XL_solve,
													  LPXLOPER XL_amortizationId)
{
	ADD_LOG("Local_FRNPrice ");
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
	double C_asOfDate;
	double C_delivery;
	double C_maturity;

	CCString C_ccy1;

	CCString C_index1;
	long liborType1Id;

	long frequencyId;
	long frequencyId2;
	
	double C_facialMargin;

	CCString C_discountCurve1;
	CCString C_forwardCurve1;
	
	CCString C_ccy2;

	CCString C_index2;
	long liborType2Id = ARM_NULL_OBJECT_ID;

	CCString C_discountCurve2;
	CCString C_forwardCurve2;

	double C_valoMargin;

	double C_fixing;
	double C_fixing_default = 0.0;

	double C_spread;
	double C_spread_default = 0.0;

	CCString C_amortizationId;

	double C_solve;
	double C_solve_default = BasisSwap_Method_Alg;

	double C_outMode;
	double C_outMode_default = 0.0;

	// error
	static int error;
	static char* reason = "";

	XL_readNumCell(XL_asOfDate,C_asOfDate," ARM_ERR: as of date: date expected",C_result);
	XL_readNumCell(XL_delivery,C_delivery," ARM_ERR: delivery: date expected",C_result);
	XL_readNumCell(XL_maturity,C_maturity," ARM_ERR: maturity: date expected",C_result);
	XL_readStrCell(XL_ccy1,C_ccy1," ARM_ERR: currency 1: string expected",C_result);
	XL_readStrCell(XL_index1,C_index1," ARM_ERR: index 1: string expected",C_result);
	XL_readNumCell(XL_facialMargin,C_facialMargin," ARM_ERR: facial margin: numeric expected",C_result);
	XL_readStrCell(XL_discountCurve1,C_discountCurve1," ARM_ERR: discount curve 1: string expected",C_result);
	XL_readStrCell(XL_forwardCurve1,C_forwardCurve1," ARM_ERR: forward curve 1: string expected",C_result);
	XL_readStrCellWD(XL_ccy2,C_ccy2,"NONE"," ARM_ERR: currency 2: string expected",C_result);
	// Modif du 21/06/2001 : prise en compte du 2eme index en monocurrency
	XL_readStrCellWD(XL_index2,C_index2,C_index1," ARM_ERR: index 2: string expected",C_result);
	XL_readNumCell(XL_valoMargin,C_valoMargin," ARM_ERR: valo margin: numeric expected",C_result);
	XL_readNumCellWD(XL_fixing,C_fixing,C_fixing_default," ARM_ERR: fixing: numeric expected",C_result);
	XL_readNumCellWD(XL_spread,C_spread,C_spread_default," ARM_ERR: spread: numeric expected",C_result);
	XL_readStrCellWD(XL_amortizationId,C_amortizationId,"DEFAULT"," ARM_ERR: amortization: object expected",C_result);
	XL_readNumCellWD(XL_outMode,C_outMode,C_outMode_default," ARM_ERR: out mode: numeric expected",C_result);
	XL_readNumCellWD(XL_solve,C_solve,C_solve_default," ARM_ERR: solve: numeric expected",C_result);
	
	if(!(C_ccy2 == "NONE"))
	{
		XL_readStrCell(XL_index2,C_index2," ARM_ERR: index 2: string expected",C_result);
		XL_readStrCell(XL_discountCurve2,C_discountCurve2," ARM_ERR: discount curve 2: string expected",C_result);
		XL_readStrCell(XL_forwardCurve2,C_forwardCurve2," ARM_ERR: forward curve 2: string expected",C_result);
		if((liborType2Id = ARM_ConvIrIndName (C_index2, C_result)) == ARM_DEFAULT_ERR)
		{
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}
	}

	if((liborType1Id = ARM_ConvIrIndName (C_index1, C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	if((frequencyId = ARM_ConvIrIndNameToFreq (C_index1, C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	if((frequencyId2 = ARM_ConvIrIndNameToFreq (C_index2, C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}
	
	long retCode = ARMLOCAL_FRNPriceNew (C_asOfDate,
								 C_delivery,
								 C_maturity,
								 C_ccy1,
								 liborType1Id,
								 LocalGetNumObjectId (C_forwardCurve1),
								 LocalGetNumObjectId (C_discountCurve1),
								 C_facialMargin,
								 C_valoMargin,
								 frequencyId,
								 C_ccy2,
								 liborType2Id,
								 LocalGetNumObjectId (C_forwardCurve2),
								 LocalGetNumObjectId (C_discountCurve2),
								 LocalGetNumObjectId (C_amortizationId),
								 frequencyId2,
								 C_fixing,
								 C_spread,
								 (long)C_solve,
								 C_result);
												
	if(retCode == ARM_OK)
	{
		if((C_outMode != 0.0) && ((C_solve == BasisSwap_Method_Num) || (C_solve == BasisSwap_Method_Num_Spr)))
		{
			int nbrows = 9;
			int nbcolumns = 2;
		
			FreeCurCellErr ();
			XL_result.xltype = xltypeMulti;
			XL_result.val.array.columns = nbcolumns;
			XL_result.val.array.rows = nbrows; 
			XL_result.val.array.lparray = pxArray = (LPXLOPER)GlobalAlloc (GMEM_ZEROINIT, nbrows * nbcolumns * sizeof (XLOPER));

			pxArray[XL_Coordonnate2Rank (0, 0, nbcolumns)].xltype = xltypeNum;
			pxArray[XL_Coordonnate2Rank (0, 0, nbcolumns)].val.num = C_result.getDouble (); 
			pxArray[XL_Coordonnate2Rank (0, 1, nbcolumns)].xltype = xltypeStr;
			pxArray[XL_Coordonnate2Rank (0, 1, nbcolumns)].val.str = XL_StrC2StrPascal ("");
			pxArray[XL_Coordonnate2Rank (0, 1, nbcolumns)].xltype |= xlbitDLLFree;
			
			pxArray[XL_Coordonnate2Rank (1, 0, nbcolumns)].xltype = xltypeStr;
			pxArray[XL_Coordonnate2Rank (1, 0, nbcolumns)].val.str = XL_StrC2StrPascal ("StartXNL1");
			pxArray[XL_Coordonnate2Rank (1, 0, nbcolumns)].xltype |= xlbitDLLFree;
			pxArray[XL_Coordonnate2Rank (1, 1, nbcolumns)].xltype = xltypeNum;
			pxArray[XL_Coordonnate2Rank (1, 1, nbcolumns)].val.num = ASP_getStartXNL1 (); 

			pxArray[XL_Coordonnate2Rank (2, 0, nbcolumns)].xltype = xltypeStr;
			pxArray[XL_Coordonnate2Rank (2, 0, nbcolumns)].val.str = XL_StrC2StrPascal ("EndXNL1");
			pxArray[XL_Coordonnate2Rank (2, 0, nbcolumns)].xltype |= xlbitDLLFree;
			pxArray[XL_Coordonnate2Rank (2, 1, nbcolumns)].xltype = xltypeNum;
			pxArray[XL_Coordonnate2Rank (2, 1, nbcolumns)].val.num = ASP_getEndXNL1 (); 

			pxArray[XL_Coordonnate2Rank (3, 0, nbcolumns)].xltype = xltypeStr;
			pxArray[XL_Coordonnate2Rank (3, 0, nbcolumns)].val.str = XL_StrC2StrPascal ("StartXNL2");
			pxArray[XL_Coordonnate2Rank (3, 0, nbcolumns)].xltype |= xlbitDLLFree;
			pxArray[XL_Coordonnate2Rank (3, 1, nbcolumns)].xltype = xltypeNum;
			pxArray[XL_Coordonnate2Rank (3, 1, nbcolumns)].val.num = ASP_getStartXNL2 (); 

			pxArray[XL_Coordonnate2Rank (4, 0, nbcolumns)].xltype = xltypeStr;
			pxArray[XL_Coordonnate2Rank (4, 0, nbcolumns)].val.str = XL_StrC2StrPascal ("EndXNL2");
			pxArray[XL_Coordonnate2Rank (4, 0, nbcolumns)].xltype |= xlbitDLLFree;
			pxArray[XL_Coordonnate2Rank (4, 1, nbcolumns)].xltype = xltypeNum;
			pxArray[XL_Coordonnate2Rank (4, 1, nbcolumns)].val.num = ASP_getEndXNL2 (); 

			pxArray[XL_Coordonnate2Rank (5, 0, nbcolumns)].xltype = xltypeStr;
			pxArray[XL_Coordonnate2Rank (5, 0, nbcolumns)].val.str = XL_StrC2StrPascal ("Price1");
			pxArray[XL_Coordonnate2Rank (5, 0, nbcolumns)].xltype |= xlbitDLLFree;
			pxArray[XL_Coordonnate2Rank (5, 1, nbcolumns)].xltype = xltypeNum;
			pxArray[XL_Coordonnate2Rank (5, 1, nbcolumns)].val.num = ASP_getPrice1 (); 

			pxArray[XL_Coordonnate2Rank (6, 0, nbcolumns)].xltype = xltypeStr;
			pxArray[XL_Coordonnate2Rank (6, 0, nbcolumns)].val.str = XL_StrC2StrPascal ("Price2");
			pxArray[XL_Coordonnate2Rank (6, 0, nbcolumns)].xltype |= xlbitDLLFree;
			pxArray[XL_Coordonnate2Rank (6, 1, nbcolumns)].xltype = xltypeNum;
			pxArray[XL_Coordonnate2Rank (6, 1, nbcolumns)].val.num = ASP_getPrice2 ();

			pxArray[XL_Coordonnate2Rank (7, 0, nbcolumns)].xltype = xltypeStr;
			pxArray[XL_Coordonnate2Rank (7, 0, nbcolumns)].val.str = XL_StrC2StrPascal ("Minimum");
			pxArray[XL_Coordonnate2Rank (7, 0, nbcolumns)].xltype |= xlbitDLLFree;
			pxArray[XL_Coordonnate2Rank (7, 1, nbcolumns)].xltype = xltypeNum;
			pxArray[XL_Coordonnate2Rank (7, 1, nbcolumns)].val.num = ASP_getMinimum ();

			pxArray[XL_Coordonnate2Rank (8, 0, nbcolumns)].xltype = xltypeStr;
			pxArray[XL_Coordonnate2Rank (8, 0, nbcolumns)].val.str = XL_StrC2StrPascal ("Iterations");
			pxArray[XL_Coordonnate2Rank (8, 0, nbcolumns)].xltype |= xlbitDLLFree;
			pxArray[XL_Coordonnate2Rank (8, 1, nbcolumns)].xltype = xltypeNum;
			pxArray[XL_Coordonnate2Rank (8, 1, nbcolumns)].val.num = ASP_getNbIter ();
		}
		else
		{
			FreeCurCellErr ();
			XL_result.xltype = xltypeNum;
			XL_result.val.num = C_result.getDouble ();
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_FRNPrice" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


__declspec(dllexport) LPXLOPER WINAPI Local_CptBPV (LPXLOPER XL_asOfDate,
													LPXLOPER XL_delivery,
													LPXLOPER XL_maturity,
													LPXLOPER XL_zc,
													LPXLOPER XL_frequency,
													LPXLOPER XL_dayCount,
													LPXLOPER XL_ccy,
													LPXLOPER XL_amortizationId)
{
	ADD_LOG("Local_CptBPV ");
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();

	// C variable
	double C_asOfDate;
	double C_delivery;
	double C_maturity;

	CCString C_ccy;
	
	CCString C_amortizationId;
	CCString C_frequency;
	long frequencyId;

	CCString C_dayCount;
	long dayCountId;

	CCString C_zc;
	
	// error
	static int error;
	static char* reason = "";

	XL_readNumCell(XL_asOfDate,C_asOfDate," ARM_ERR: as of date: date expected",C_result);
	XL_readNumCell(XL_delivery,C_delivery," ARM_ERR: delivery: date expected",C_result);
	XL_readNumCell(XL_maturity,C_maturity," ARM_ERR: maturity: date expected",C_result);
	XL_readStrCell(XL_zc,C_zc," ARM_ERR: zc: string expected",C_result);
	XL_readStrCellWD(XL_frequency,C_frequency,"1"," ARM_ERR: frequency: string expected",C_result);
	XL_readStrCellWD(XL_dayCount,C_dayCount,"1"," ARM_ERR: day count: string expected",C_result);
	XL_readStrCellWD(XL_ccy,C_ccy,"DEFAULT"," ARM_ERR: currency: string expected",C_result);
	XL_readStrCellWD(XL_amortizationId,C_amortizationId,"DEFAULT"," ARM_ERR: amortization: string expected",C_result);
		
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

	if((frequencyId = ARM_ConvFrequency (C_frequency, C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	dayCountId = ARM_ConvDayCount (C_dayCount);
	
	long retCode = ARMLOCAL_CptBPV (C_asOfDate,
							   C_delivery,
							   C_maturity,
							   LocalGetNumObjectId (C_zc),
							   frequencyId,
							   dayCountId,
							   C_ccy,
							   LocalGetNumObjectId (C_amortizationId),
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_CptBPV" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}




__declspec(dllexport) LPXLOPER WINAPI Local_ARM_NextCpnDate (LPXLOPER XL_asOfDate,
															 LPXLOPER XL_maturity,
															 LPXLOPER XL_frequency,
															 LPXLOPER XL_rule,
															 LPXLOPER XL_ccy,
															 LPXLOPER XL_intrule)
{
	ADD_LOG("Local_ARM_NextCpnDate ");
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();

	// C variable
	double C_asOfDate;

	double C_maturity;
	
	CCString C_frequency;
	long frequencyId;

	CCString C_rule;
	long ruleId;
	
	CCString C_ccy;
	
	CCString C_intrule;
	long intruleId;

	// error
	static int error;
	static char* reason = "";

	XL_readNumCell(XL_asOfDate,C_asOfDate," ARM_ERR: as of date: date expected",C_result);
	XL_readNumCell(XL_maturity,C_maturity," ARM_ERR: maturity: date expected",C_result);
	XL_readStrCell(XL_frequency,C_frequency," ARM_ERR: frequency: string expected",C_result);
	XL_readStrCellWD(XL_rule,C_rule,"F"," ARM_ERR: rule: string expected",C_result);
	XL_readStrCellWD(XL_ccy,C_ccy,"DEFAULT"," ARM_ERR: currency: string expected",C_result);
	XL_readStrCellWD(XL_intrule,C_intrule,"ADJ"," ARM_ERR: intrule: string expected",C_result);

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

	if((ruleId = ARM_ConvRule (C_rule, C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}
	
	if((frequencyId = ARM_ConvFrequency (C_frequency, C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	intruleId = ARM_ConvIntRule (C_intrule);

	long retCode = ARMLOCAL_ARM_NextCpnDate (C_asOfDate,
											 C_maturity,
											 frequencyId,
											 ruleId,
											 C_ccy,
											 C_result,
											 intruleId);
											 	
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_ARM_NextCpnDate" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


__declspec(dllexport) LPXLOPER WINAPI Local_ARM_PrevCpnDate (LPXLOPER XL_asOfDate,
															 LPXLOPER XL_maturity,
															 LPXLOPER XL_frequency,
															 LPXLOPER XL_rule,
															 LPXLOPER XL_ccy,
															 LPXLOPER XL_intrule)
{
	ADD_LOG("Local_ARM_PrevCpnDate ");
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();

	// C variable
	double C_asOfDate;

	double C_maturity;
	
	CCString C_frequency;
	long frequencyId;

	CCString C_rule;
	long ruleId;

	CCString C_intrule;
	long intruleId;

	CCString C_ccy;
	
	// error
	static int error;
	static char* reason = "";

	XL_readNumCell(XL_asOfDate,C_asOfDate," ARM_ERR: as of date: date expected",C_result);
	XL_readNumCell(XL_maturity,C_maturity," ARM_ERR: maturity: date expected",C_result);
	XL_readStrCell(XL_frequency,C_frequency," ARM_ERR: frequency: string expected",C_result);
	XL_readStrCellWD(XL_rule,C_rule,"F"," ARM_ERR: rule: string expected",C_result);
	XL_readStrCellWD(XL_ccy,C_ccy,"DEFAULT"," ARM_ERR: currency: string expected",C_result);
	XL_readStrCellWD(XL_intrule,C_intrule,"ADJ"," ARM_ERR: intrule: string expected",C_result);

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

	if((ruleId = ARM_ConvRule (C_rule, C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}
	
	if((frequencyId = ARM_ConvFrequency (C_frequency, C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	intruleId = ARM_ConvIntRule (C_intrule);

	long retCode = ARMLOCAL_ARM_PrevCpnDate (C_asOfDate,
											 C_maturity,
											 frequencyId,
											 ruleId,
											 C_ccy,
											 C_result,
											 intruleId);

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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_ARM_PrevCpnDate" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}
