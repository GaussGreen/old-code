#pragma warning(disable : 4786)
#pragma warning(disable : 4005)

#include <libCCxll\CCxll.h>

#include <ARM\libarm_local\ARM_local_mod.h>
#include <ARM\libarm_local\ARM_local_class.h>
#include <ARM\libarm_local\ARM_local_glob.h>
#include <ARM\libarm_local\ARM_local_zccurve.h>
#include <ARM\libarm_local\ARM_local_volcrv.h>
#include <ARM\libarm_local\ARM_local_dkprcs.h>

#include <ARM\libarm\ARM_result.h>

#include "ARM_local_interglob.h"
#include "ARM_local_interface.h"

#include "ARM_xl_trycatch_local.h"
#include <util\fromto.h>


#include "XL_local_xlarm_common.h"

#include "ARM_xl_dkprcs_local.h"

#include "util\tech_macro.h"

static HGLOBAL Mmglob = NULL;







__declspec(dllexport) LPXLOPER WINAPI Local_PRCS3F_LatticePricing(LPXLOPER XL_dLatticeGeometryDataIn,
		 								    LPXLOPER XL_dNumTimeLinesBeforeFirstNotice,
										    LPXLOPER XL_dNumTimeLinesPerYearAfterFirstNoticeAndBeforeOptimal,
										    LPXLOPER XL_dNumTimeLines, 
										    LPXLOPER XL_evalDate,
                                            LPXLOPER XL_curves,
											/*
										    LPXLOPER XL_dBaseYieldCurveId,
										    LPXLOPER XL_dForeignYieldCurveId,
										    LPXLOPER XL_volSwopBaseId,
                                            LPXLOPER XL_volSwopForeignId,
                                            LPXLOPER XL_volFxId,
											*/
										    LPXLOPER XL_dNoticeDatesIn,
										    LPXLOPER XL_dStrike,
										    LPXLOPER XL_dType,
										    LPXLOPER XL_dOptionExpiry, 
										    LPXLOPER XL_dMeanReversionBase,
										    LPXLOPER XL_dMeanReversionForeign,  
										    LPXLOPER XL_dSpotFX,
										    LPXLOPER XL_dBaseForeignCorrelationId,
										    LPXLOPER XL_dBaseSpotFXCorrelationId,
										    LPXLOPER XL_dForeignSpotFXCorrelationId,
										    LPXLOPER XL_dProductModelCodeId,
                                            LPXLOPER XL_paraModel,
										    LPXLOPER XL_dBoosterDataIn)
{
	ADD_LOG("Local_PRCS3F_LatticePricing");
	//	ARM_BEGIN();

	
	static XLOPER XL_result;
    LPXLOPER pxArray;

	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();

	// C variables 
    VECTOR<CCString> C_curves;
    VECTOR<double> C_dLatticeGeometryDataIn;
	VECTOR<double> C_dNoticeDatesIn;
	VECTOR<double> C_dBoosterDataIn;
	VECTOR<double> C_paraModel;

	long C_dBaseYieldCurveId;
	long C_dForeignYieldCurveId;
	long C_dBaseRatesNoBasisCurveId;
    long C_dForeignRatesNoBasisCurveId;
    long C_volSwopBaseId;
    long C_volSwopForeignId;
    long C_volFxId;
     
	CCString C_dBaseForeignCorrelationId;
	CCString C_dBaseSpotFXCorrelationId;
	CCString C_dForeignSpotFXCorrelationId;

	double C_dNumTimeLinesBeforeFirstNotice;
	double C_dNumTimeLinesPerYearAfterFirstNoticeAndBeforeOptimal;
    double C_dNumTimeLines; 
    double C_evalDate; 
	double C_dStrike;
	double C_dType;
	double C_dOptionExpiry; 
	double C_dMeanReversionBase;
	double C_dMeanReversionForeign;
	double C_dSpotFX;
	double C_dProductModelCodeId;
	double C_dX1Limit;
	double C_dX2Limit;
	double C_dX3Limit;
	double C_dI1Limit;
	double C_dI2Limit;
	double C_dI3Limit;
	double C_dADPLimit;
	double C_dOptimal;
	double C_dTimeBoost;
	double C_dDeltaFlag;
	double C_dStringModel;

	long nbBoosterRows;


	// error
	static int error;
	static char* reason = "";


  

	// a Vector
	XL_readNumVector(XL_dLatticeGeometryDataIn, 
		             C_dLatticeGeometryDataIn," ARM_ERR: LatticeGeometry, array of numeric expected", C_result);

	// Doubles
	XL_readNumCell(XL_dNumTimeLinesBeforeFirstNotice,
		           C_dNumTimeLinesBeforeFirstNotice," ARM_ERR: NumTimeLinesBeforeFirstNotice, numeric expected", C_result);

    XL_readNumCell(XL_dNumTimeLinesPerYearAfterFirstNoticeAndBeforeOptimal,
		           C_dNumTimeLinesPerYearAfterFirstNoticeAndBeforeOptimal,
				   " ARM_ERR: NumTimeLinesPerYearAfterFirstNoticeAndBeforeOptimal, numeric expected", C_result);

    XL_readNumCell(XL_dNumTimeLines,
		           C_dNumTimeLines," ARM_ERR: NumTimeLines, numeric expected", C_result);

	XL_readNumCell(XL_evalDate,
		           C_evalDate," ARM_ERR: (Excel Date) numeric expected", C_result);


    // Objects Curves
	/*
	XL_readStrCell(XL_dBaseYieldCurveId, C_dBaseYieldCurveId," ARM_ERR: BaseYieldCurveId, Curve Id: object expected",C_result);
    XL_readStrCell(XL_dForeignYieldCurveId, C_dForeignYieldCurveId," ARM_ERR: ForeignYieldCurveId, Curve Id: object expected",C_result);
    XL_dBaseRatesNoBasisCurveId,
    XL_dForeignRatesNoBasisCurveId
	XL_readStrCell(XL_volSwopBaseId, C_volSwopBaseId," ARM_ERR: volSwopBaseId, Curve Id: object expected",C_result);
    XL_readStrCell(XL_volSwopForeignId, C_volSwopForeignId," ARM_ERR: volSwopForeignId, Curve Id: object expected",C_result);
    XL_readStrCell(XL_volFxId, C_volFxId," ARM_ERR: volFxId, Curve Id: object expected",C_result);
    */
	XL_readStrVector(XL_curves, C_curves," ARM_ERR: Yield and vols. curves: array of 5 objects expected",DOUBLE_TYPE,C_result);


    if ( C_curves.size() != 7 )
	{
	   // Error

	   XL_result.xltype = xltypeStr;
	   XL_result.val.str = "\007ARM_ERR";

	   SetCurCellErrValue("Curves Vector size must be 7");

	   return((LPXLOPER) &XL_result);
	}
	else
	{
	   int i = 0;

       C_dBaseYieldCurveId    = LocalGetNumObjectId(C_curves[i++]);

       C_dForeignYieldCurveId = LocalGetNumObjectId(C_curves[i++]);

	   C_dBaseRatesNoBasisCurveId = LocalGetNumObjectId(C_curves[i++]);
       
	   C_dForeignRatesNoBasisCurveId = LocalGetNumObjectId(C_curves[i++]);
       
	   C_volSwopBaseId        = LocalGetNumObjectId(C_curves[i++]);

       C_volSwopForeignId     = LocalGetNumObjectId(C_curves[i++]);
 
       C_volFxId              = LocalGetNumObjectId(C_curves[i++]);	   
	}


    // a Vector
	XL_readNumVector(XL_dNoticeDatesIn, 
		             C_dNoticeDatesIn," ARM_ERR: NoticeDatesIn, array of numeric expected", C_result);

	// Doubles
	XL_readNumCell(XL_dStrike,
		           C_dStrike," ARM_ERR: Strike, numeric expected", C_result);
     
	XL_readNumCell(XL_dType,
		           C_dType," ARM_ERR: Type, numeric expected", C_result);

	XL_readNumCell(XL_dOptionExpiry,
		           C_dOptionExpiry," ARM_ERR: OptionExpiry, numeric expected", C_result);

    XL_readNumCell(XL_dMeanReversionBase,
		           C_dMeanReversionBase," ARM_ERR: MeanReversionBase, numeric expected", C_result);

	XL_readNumCell(XL_dMeanReversionForeign,
		           C_dMeanReversionForeign, " ARM_ERR: MeanReversionForeign, numeric expected", C_result);
	
	
	XL_readNumCell(XL_dSpotFX,
		           C_dSpotFX, " ARM_ERR: SpotFX, numeric expected", C_result);

	// Correlation objects
    XL_readStrCell(XL_dBaseForeignCorrelationId, C_dBaseForeignCorrelationId,
		           " ARM_ERR: BaseForeignCorrelationId, Vol. Curve Id: object expected",C_result);

    XL_readStrCell(XL_dBaseSpotFXCorrelationId, C_dBaseSpotFXCorrelationId,
		           " ARM_ERR: BaseSpotFXCorrelationId, Vol. Curve Id: object expected", C_result);

    XL_readStrCell(XL_dForeignSpotFXCorrelationId, C_dForeignSpotFXCorrelationId,
		           " ARM_ERR: ForeignSpotFXCorrelationId, Curve Id: object expected", C_result);

	// a Double 
	XL_readNumCell(XL_dProductModelCodeId,
		           C_dProductModelCodeId, " ARM_ERR: ProductModelCodeId, numeric expected", C_result);


	// model parameters vector
	XL_readNumVector(XL_paraModel, 
		             C_paraModel," ARM_ERR: paraModel, array of numeric expected", C_result);
    if ( C_paraModel.size() != 11 )
	{
	   // Error

	   XL_result.xltype = xltypeStr;
	   XL_result.val.str = "\007ARM_ERR";
	   SetCurCellErrValue("paraModel size must be 11");

	   return((LPXLOPER) &XL_result);
	}
	else
	{
	   int i = 0;

	   C_dX1Limit     = C_paraModel[i++]; 
       C_dX2Limit     = C_paraModel[i++];
	   C_dX3Limit     = C_paraModel[i++];
	   C_dI1Limit     = C_paraModel[i++];
	   C_dI2Limit     = C_paraModel[i++];
	   C_dI3Limit     = C_paraModel[i++];
	   C_dADPLimit    = C_paraModel[i++];
	   C_dOptimal     = C_paraModel[i++];
	   C_dTimeBoost   = C_paraModel[i++];
	   C_dDeltaFlag   = C_paraModel[i++];
	   C_dStringModel = C_paraModel[i++];
	}

	// Matrix

    long nbLin;
    long nbCol;

    if (( error = XL_getNumVectorAndSize(XL_dBoosterDataIn, nbLin, nbCol,C_dBoosterDataIn)) != XL_NO_ERROR)
	{
	   if ( error != XL_NO_ERROR )
       {
		  CCString local_msg(" ARM_ERR: BoosterDataIn, array of numeric expected");

          C_result.setMsg(local_msg);

          ARM_ARG_ERR();
		
          return((LPXLOPER) &XL_result);
       }
	}

    nbBoosterRows = nbCol;
    
    // nbBoosterRows = XL_dBoosterDataIn->val.array.rows;

	long retCode = ARMLOCAL_PRCS3F_Lattice_HWVFDK_Pricing(C_dLatticeGeometryDataIn,
		 						                          C_dNumTimeLinesBeforeFirstNotice,
								                          C_dNumTimeLinesPerYearAfterFirstNoticeAndBeforeOptimal,
								                          C_dNumTimeLines, 
								                          C_evalDate, 
								                          C_dBaseYieldCurveId,
								                          C_dForeignYieldCurveId,
														  C_dBaseRatesNoBasisCurveId,
                                                          C_dForeignRatesNoBasisCurveId,
								                          C_volSwopBaseId,
                                                          C_volSwopForeignId,
                                                          C_volFxId,
								                          C_dNoticeDatesIn,
							 	                          C_dStrike,
								                          C_dType,
								                          C_dOptionExpiry, 
								                          C_dMeanReversionBase,
								                          C_dMeanReversionForeign,  
								                          C_dSpotFX,
								                          LocalGetNumObjectId(C_dBaseForeignCorrelationId),
								                          LocalGetNumObjectId(C_dBaseSpotFXCorrelationId),
								                          LocalGetNumObjectId(C_dForeignSpotFXCorrelationId),
								                          C_dProductModelCodeId,
								                          C_dX1Limit,
								                          C_dX2Limit,
								                          C_dX3Limit,
								                          C_dI1Limit,
								                          C_dI2Limit,
								                          C_dI3Limit,
								                          C_dADPLimit,
								                          C_dOptimal,
								                          C_dTimeBoost,
								                          C_dDeltaFlag,
								                          C_dStringModel,
								                          nbBoosterRows,
								                          C_dBoosterDataIn,
								                          C_result);

	if ( retCode == ARM_OK )
	{
		int nbrows    = C_result.getLong(); // in fact the size
		int nbcolumns = 1;
	
		FreeCurCellErr ();

		XL_result.xltype = xltypeMulti;
		XL_result.val.array.columns = nbcolumns;
		XL_result.val.array.rows = nbrows; 
		XL_result.val.array.lparray = pxArray = (LPXLOPER) GlobalAlloc(GMEM_ZEROINIT, 
			                                             nbrows*nbcolumns*sizeof(XLOPER));

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
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_PRCS3F_LatticePricing" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}



__declspec(dllexport) LPXLOPER WINAPI Local_ARM_TREE3FACT(LPXLOPER XL_asof,
														  LPXLOPER XL_xbsfx,
														  LPXLOPER XL_volSwopBaseId,
														  LPXLOPER XL_volSwopForeignId,
														  LPXLOPER XL_dBaseForeignCorrelationId,
														  LPXLOPER XL_dLatticeGeometryDataIn,
														  LPXLOPER XL_dNumTimeLinesBeforeFirstNotice,
														  LPXLOPER XL_dNumTimeLinesPerYearAfterFirstNoticeAndBeforeOptimal,
														  LPXLOPER XL_dNumTimeLines,
														  LPXLOPER XL_dMeanReversionBase,
														  LPXLOPER XL_dMeanReversionForeign,
														  LPXLOPER XL_limitVector,
														  LPXLOPER XL_dOptimal,
														  LPXLOPER XL_dTimeBoost,
														  LPXLOPER XL_dDeltaFlag,
														  LPXLOPER XL_dSmoothingValue,
														  LPXLOPER XL_sCalcProbaSurvOrNot,
														  LPXLOPER XL_QBaseSmile,
				                                          LPXLOPER XL_QForeignSmile,
														  LPXLOPER XL_convInVolFwd)
{
	ADD_LOG("Local_ARM_TREE3FACT");
	// return
	static XLOPER XL_result;
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();

	// C variable 
	double C_asof;

	CCString C_xbxfx;
	CCString C_volSwopBaseId;
	CCString C_volSwopForeignId;
	CCString C_dBaseForeignCorrelationId;

	VECTOR<double> C_dLatticeGeometryDataIn;

	double C_dNumTimeLinesBeforeFirstNotice;
	double C_dNumTimeLinesPerYearAfterFirstNoticeAndBeforeOptimal;
    double C_dNumTimeLines; 
	double C_dMeanReversionBase;
	double C_dMeanReversionForeign;

	VECTOR<double> C_limitVector;
	double C_dX1Limit;
	double C_dX2Limit;
	double C_dX3Limit;
	double C_dI1Limit;
	double C_dI2Limit;
	double C_dI3Limit;
	double C_dADPLimit;
	long   C_spotVolOrFwdVol;
	double C_histoVolLongTerm;

	double C_dOptimal;
	double C_dTimeBoost;
	double C_dDeltaFlag;
	double C_dSmoothingValue;

	CCString C_sCalcProbaSurvOrNot;
	long calcProbaSurvOrNotId;

	double C_QBaseSmile;
	double C_QForeignSmile;

    double zeroVal = 0.0; 

	CCString C_convInVolFwd;
	long convInVolFwdId;

	long C_calibBasisIncluded = 1;

	// error
	static int error;
	static char* reason = "";

	XL_readNumCell(XL_asof,C_asof," ARM_ERR: as of date: numeric expected",C_result);

	XL_readStrCell(XL_xbsfx,C_xbxfx," ARM_ERR: xbsfx model id: object expected",C_result);
	XL_readStrCell(XL_volSwopBaseId,C_volSwopBaseId," ARM_ERR: dom swopt vol id: object expected",C_result);
	XL_readStrCell(XL_volSwopForeignId,C_volSwopForeignId," ARM_ERR: for swopt vol id: object expected",C_result);
	XL_readStrCell(XL_dBaseForeignCorrelationId,C_dBaseForeignCorrelationId," ARM_ERR: correlation id: object expected",C_result);

	XL_readNumVector(XL_dLatticeGeometryDataIn, C_dLatticeGeometryDataIn," ARM_ERR: LatticeGeometry :  array of numeric expected", C_result);

	XL_readNumCell(XL_dNumTimeLinesBeforeFirstNotice, C_dNumTimeLinesBeforeFirstNotice," ARM_ERR: NumTimeLinesBeforeFirstNotice : numeric expected", C_result);
	XL_readNumCell(XL_dNumTimeLinesPerYearAfterFirstNoticeAndBeforeOptimal, C_dNumTimeLinesPerYearAfterFirstNoticeAndBeforeOptimal,
				   " ARM_ERR: NumTimeLinesPerYearAfterFirstNoticeAndBeforeOptimal : numeric expected", C_result);
	XL_readNumCell(XL_dNumTimeLines, C_dNumTimeLines," ARM_ERR: NumTimeLines : numeric expected", C_result);
	XL_readNumCell(XL_dMeanReversionBase, C_dMeanReversionBase," ARM_ERR: MeanReversionBase : numeric expected", C_result);
	XL_readNumCell(XL_dMeanReversionForeign, C_dMeanReversionForeign, " ARM_ERR: MeanReversionForeign : numeric expected", C_result);

	// model parameters vector
	XL_readNumVector(XL_limitVector, C_limitVector," ARM_ERR: limit Vector : array of numeric expected", C_result);

	if ( ( C_limitVector.size() != 9 ) && ( C_limitVector.size() != 10 ) )
	{
	   // Error

	   XL_result.xltype = xltypeStr;
	   XL_result.val.str = "\007ARM_ERR";
	   SetCurCellErrValue("limit vector size must be 9");

	   return((LPXLOPER) &XL_result);
	}
	else
	{
		int i = 0;

		C_dX1Limit			= C_limitVector[i++]; 
		C_dX2Limit			= C_limitVector[i++];
		C_dX3Limit			= C_limitVector[i++];
		C_dI1Limit			= C_limitVector[i++];
		C_dI2Limit			= C_limitVector[i++];
		C_dI3Limit			= C_limitVector[i++];
		C_dADPLimit			= C_limitVector[i++];
		C_spotVolOrFwdVol	= (long) C_limitVector[i++];
		C_histoVolLongTerm	= C_limitVector[i++];

		if (C_limitVector.size() == 10)
		{
			C_calibBasisIncluded = C_limitVector[i++];
		}
	}

	XL_readNumCell(XL_dOptimal, C_dOptimal," ARM_ERR: optimal : numeric expected", C_result);
	XL_readNumCell(XL_dTimeBoost, C_dTimeBoost," ARM_ERR: time boost : numeric expected", C_result);
	XL_readNumCell(XL_dDeltaFlag, C_dDeltaFlag, " ARM_ERR: delta Flag : numeric expected", C_result);
	XL_readNumCell(XL_dSmoothingValue, C_dSmoothingValue, " ARM_ERR: smoothing value : numeric expected", C_result);

	XL_readStrCellWD(XL_sCalcProbaSurvOrNot,C_sCalcProbaSurvOrNot,"YES"," ARM_ERR: Calc proba survive or not: string expected",C_result);

	XL_readNumCellWD(XL_QBaseSmile, C_QBaseSmile, zeroVal,
		             " ARM_ERR: Q-Value for Base Smile : numeric expected", C_result);

    XL_readNumCellWD(XL_QForeignSmile, C_QForeignSmile, zeroVal,
		             " ARM_ERR: Q-Value for Foreign Smile : numeric expected", C_result);

	XL_readStrCellWD(XL_convInVolFwd,C_convInVolFwd,"N"," ARM_ERR: conversion in Vol Forward ?: string expected",C_result);

	if ((calcProbaSurvOrNotId = ARM_ConvYesOrNo (C_sCalcProbaSurvOrNot, C_result)) == ARM_DEFAULT_ERR)
	{
	   ARM_ARG_ERR();
		
	   return (LPXLOPER)&XL_result;
	}

	if ((convInVolFwdId = ARM_ConvYesOrNo (C_convInVolFwd, C_result)) == ARM_DEFAULT_ERR)
	{
	   ARM_ARG_ERR();
		
	   return (LPXLOPER)&XL_result;
	}

	long retCode;
	long objId;
	CCString prevClass;
	
	CCString curClass = LOCAL_TREE_3FACT_CLASS;
	CCString stringId = GetLastCurCellEnvValue ();
	
	if (!stringId)
	{
		retCode = ARMLOCAL_TREE3FACT (C_asof,
									  LocalGetNumObjectId (C_xbxfx),
									  LocalGetNumObjectId (C_volSwopBaseId),
									  LocalGetNumObjectId (C_volSwopForeignId),
									  LocalGetNumObjectId (C_dBaseForeignCorrelationId),
									  C_dLatticeGeometryDataIn,
									  C_dNumTimeLinesBeforeFirstNotice,
									  C_dNumTimeLinesPerYearAfterFirstNoticeAndBeforeOptimal,
									  C_dNumTimeLines,
									  C_dMeanReversionBase,
									  C_dMeanReversionForeign,
									  C_dX1Limit,
									  C_dX2Limit,
									  C_dX3Limit,
									  C_dI1Limit,
									  C_dI2Limit,
									  C_dI3Limit,
									  C_dADPLimit,
									  C_dOptimal,
									  C_dTimeBoost,
									  C_dDeltaFlag,
									  C_dSmoothingValue,
									  calcProbaSurvOrNotId,
                                      C_QBaseSmile,
                                      C_QForeignSmile,
									  C_spotVolOrFwdVol,
									  C_histoVolLongTerm,
									  convInVolFwdId,
									  C_calibBasisIncluded,
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
			retCode = ARMLOCAL_TREE3FACT (C_asof,
										  LocalGetNumObjectId (C_xbxfx),
										  LocalGetNumObjectId (C_volSwopBaseId),
										  LocalGetNumObjectId (C_volSwopForeignId),
										  LocalGetNumObjectId (C_dBaseForeignCorrelationId),
										  C_dLatticeGeometryDataIn,
										  C_dNumTimeLinesBeforeFirstNotice,
										  C_dNumTimeLinesPerYearAfterFirstNoticeAndBeforeOptimal,
										  C_dNumTimeLines,
										  C_dMeanReversionBase,
										  C_dMeanReversionForeign,
										  C_dX1Limit,
										  C_dX2Limit,
										  C_dX3Limit,
										  C_dI1Limit,
										  C_dI2Limit,
										  C_dI3Limit,
										  C_dADPLimit,
										  C_dOptimal,
										  C_dTimeBoost,
										  C_dDeltaFlag,
										  C_dSmoothingValue,
										  calcProbaSurvOrNotId,
										  C_QBaseSmile,
                                          C_QForeignSmile,
										  C_spotVolOrFwdVol,
										  C_histoVolLongTerm,
										  convInVolFwdId,
										  C_calibBasisIncluded,
										  C_result,
										  objId);

			if ( retCode == ARM_OK )
			{			
				LocalSetCurCellEnvValue (curClass, objId); 

				stringId = LocalMakeObjectId (objId, curClass);
			}
		}
		else
		{
			FreeCurCellContent ();

			retCode = ARMLOCAL_TREE3FACT (C_asof,
										  LocalGetNumObjectId (C_xbxfx),
										  LocalGetNumObjectId (C_volSwopBaseId),
										  LocalGetNumObjectId (C_volSwopForeignId),
										  LocalGetNumObjectId (C_dBaseForeignCorrelationId),
										  C_dLatticeGeometryDataIn,
										  C_dNumTimeLinesBeforeFirstNotice,
										  C_dNumTimeLinesPerYearAfterFirstNoticeAndBeforeOptimal,
										  C_dNumTimeLines,
										  C_dMeanReversionBase,
										  C_dMeanReversionForeign,
										  C_dX1Limit,
										  C_dX2Limit,
										  C_dX3Limit,
										  C_dI1Limit,
										  C_dI2Limit,
										  C_dI3Limit,
										  C_dADPLimit,
										  C_dOptimal,
										  C_dTimeBoost,
										  C_dDeltaFlag,
										  C_dSmoothingValue,
										  calcProbaSurvOrNotId,
										  C_QBaseSmile,
                                          C_QForeignSmile,
										  C_spotVolOrFwdVol,
										  C_histoVolLongTerm,
										  convInVolFwdId,
										  C_calibBasisIncluded,
										  C_result);

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
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_ARM_TREE3FACT" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


__declspec(dllexport) LPXLOPER WINAPI Local_PXL_ARM_TREE3FACT(LPXLOPER XL_asof,
															  LPXLOPER XL_xbsfx,
															  LPXLOPER XL_volSwopBaseId,
															  LPXLOPER XL_volSwopForeignId,
															  LPXLOPER XL_dBaseForeignCorrelationId,
															  LPXLOPER XL_dLatticeGeometryDataIn,
															  LPXLOPER XL_dNumTimeLinesBeforeFirstNotice,
															  LPXLOPER XL_dNumTimeLinesPerYearAfterFirstNoticeAndBeforeOptimal,
															  LPXLOPER XL_dNumTimeLines,
															  LPXLOPER XL_dMeanReversionBase,
															  LPXLOPER XL_dMeanReversionForeign,
															  LPXLOPER XL_limitVector,
															  LPXLOPER XL_dOptimal,
															  LPXLOPER XL_dTimeBoost,
															  LPXLOPER XL_dDeltaFlag,
															  LPXLOPER XL_dSmoothingValue,
															  LPXLOPER XL_sCalcProbaSurvOrNot,
															  LPXLOPER XL_QBaseSmile,
				                                              LPXLOPER XL_QForeignSmile,
															  LPXLOPER XL_convInVolFwd)
{
	ADD_LOG("Local_PXL_ARM_TREE3FACT");
	// return
	static XLOPER XL_result;
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();

	// C variable 
	double C_asof;

	CCString C_xbxfx;
	CCString C_volSwopBaseId;
	CCString C_volSwopForeignId;
	CCString C_dBaseForeignCorrelationId;

	VECTOR<double> C_dLatticeGeometryDataIn;

	double C_dNumTimeLinesBeforeFirstNotice;
	double C_dNumTimeLinesPerYearAfterFirstNoticeAndBeforeOptimal;
    double C_dNumTimeLines; 
	double C_dMeanReversionBase;
	double C_dMeanReversionForeign;

	VECTOR<double> C_limitVector;
	double C_dX1Limit;
	double C_dX2Limit;
	double C_dX3Limit;
	double C_dI1Limit;
	double C_dI2Limit;
	double C_dI3Limit;
	double C_dADPLimit;
	long   C_spotVolOrFwdVol;
	double C_histoVolLongTerm;

	double C_dOptimal;
	double C_dTimeBoost;
	double C_dDeltaFlag;
	double C_dSmoothingValue;

	CCString C_sCalcProbaSurvOrNot;
	long calcProbaSurvOrNotId;
	
	double C_QBaseSmile;

	double C_QForeignSmile;

    double zeroVal = 0.0;

	CCString C_convInVolFwd;
	long convInVolFwdId;

	long C_calibBasisIncluded = 1;

	// error
	static int error;
	static char* reason = "";

	XL_readNumCell(XL_asof,C_asof," ARM_ERR: as of date: numeric expected",C_result);

	XL_readStrCell(XL_xbsfx,C_xbxfx," ARM_ERR: xbsfx model id: object expected",C_result);
	XL_readStrCell(XL_volSwopBaseId,C_volSwopBaseId," ARM_ERR: dom swopt vol id: object expected",C_result);
	XL_readStrCell(XL_volSwopForeignId,C_volSwopForeignId," ARM_ERR: for swopt vol id: object expected",C_result);
	XL_readStrCell(XL_dBaseForeignCorrelationId,C_dBaseForeignCorrelationId," ARM_ERR: correlation id: object expected",C_result);

	XL_readNumVector(XL_dLatticeGeometryDataIn, C_dLatticeGeometryDataIn," ARM_ERR: LatticeGeometry :  array of numeric expected", C_result);

	XL_readNumCell(XL_dNumTimeLinesBeforeFirstNotice, C_dNumTimeLinesBeforeFirstNotice," ARM_ERR: NumTimeLinesBeforeFirstNotice : numeric expected", C_result);
	XL_readNumCell(XL_dNumTimeLinesPerYearAfterFirstNoticeAndBeforeOptimal, C_dNumTimeLinesPerYearAfterFirstNoticeAndBeforeOptimal,
				   " ARM_ERR: NumTimeLinesPerYearAfterFirstNoticeAndBeforeOptimal : numeric expected", C_result);
	XL_readNumCell(XL_dNumTimeLines, C_dNumTimeLines," ARM_ERR: NumTimeLines : numeric expected", C_result);
	XL_readNumCell(XL_dMeanReversionBase, C_dMeanReversionBase," ARM_ERR: MeanReversionBase : numeric expected", C_result);
	XL_readNumCell(XL_dMeanReversionForeign, C_dMeanReversionForeign, " ARM_ERR: MeanReversionForeign : numeric expected", C_result);

	// model parameters vector
	XL_readNumVector(XL_limitVector, C_limitVector," ARM_ERR: limit Vector : array of numeric expected", C_result);

	if ( ( C_limitVector.size() != 9 ) && ( C_limitVector.size() != 10 ) )
	{
	   // Error

	   XL_result.xltype = xltypeStr;
	   XL_result.val.str = "\007ARM_ERR";
	   SetCurCellErrValue("limit vector size must be 9");

	   return((LPXLOPER) &XL_result);
	}
	else
	{
		int i = 0;

		C_dX1Limit			= C_limitVector[i++]; 
		C_dX2Limit			= C_limitVector[i++];
		C_dX3Limit			= C_limitVector[i++];
		C_dI1Limit			= C_limitVector[i++];
		C_dI2Limit			= C_limitVector[i++];
		C_dI3Limit			= C_limitVector[i++];
		C_dADPLimit			= C_limitVector[i++];
		C_spotVolOrFwdVol	= (long) C_limitVector[i++];
		C_histoVolLongTerm	= C_limitVector[i++];

		if (C_limitVector.size() == 10)
		{
			C_calibBasisIncluded = C_limitVector[i++];
		}
	}

	XL_readNumCell(XL_dOptimal, C_dOptimal," ARM_ERR: optimal : numeric expected", C_result);
	XL_readNumCell(XL_dTimeBoost, C_dTimeBoost," ARM_ERR: time boost : numeric expected", C_result);
	XL_readNumCell(XL_dDeltaFlag, C_dDeltaFlag, " ARM_ERR: delta Flag : numeric expected", C_result);
	XL_readNumCell(XL_dSmoothingValue, C_dSmoothingValue, " ARM_ERR: smoothing value : numeric expected", C_result);

	XL_readStrCellWD(XL_sCalcProbaSurvOrNot,C_sCalcProbaSurvOrNot,"YES"," ARM_ERR: Calc proba survive or not: string expected",C_result);

	XL_readNumCellWD(XL_QBaseSmile, C_QBaseSmile, zeroVal,
		             " ARM_ERR: Q-Value for Base Smile : numeric expected", C_result);

    XL_readNumCellWD(XL_QForeignSmile, C_QForeignSmile, zeroVal,
		             " ARM_ERR: Q-Value for Foreign Smile : numeric expected", C_result);

	XL_readStrCellWD(XL_convInVolFwd,C_convInVolFwd,"N"," ARM_ERR: conversion in Vol Forward ?: string expected",C_result);

	if ((calcProbaSurvOrNotId = ARM_ConvYesOrNo (C_sCalcProbaSurvOrNot, C_result)) == ARM_DEFAULT_ERR)
	{
	   ARM_ARG_ERR();
		
	   return (LPXLOPER)&XL_result;
	}

	if ((convInVolFwdId = ARM_ConvYesOrNo (C_convInVolFwd, C_result)) == ARM_DEFAULT_ERR)
	{
	   ARM_ARG_ERR();
		
	   return (LPXLOPER)&XL_result;
	}

	long retCode;
	long objId;
	
	CCString curClass = LOCAL_TREE_3FACT_CLASS;

	CCString stringId;
	
	retCode = ARMLOCAL_TREE3FACT (C_asof,
								  LocalGetNumObjectId (C_xbxfx),
								  LocalGetNumObjectId (C_volSwopBaseId),
								  LocalGetNumObjectId (C_volSwopForeignId),
								  LocalGetNumObjectId (C_dBaseForeignCorrelationId),
								  C_dLatticeGeometryDataIn,
								  C_dNumTimeLinesBeforeFirstNotice,
								  C_dNumTimeLinesPerYearAfterFirstNoticeAndBeforeOptimal,
								  C_dNumTimeLines,
								  C_dMeanReversionBase,
								  C_dMeanReversionForeign,
								  C_dX1Limit,
								  C_dX2Limit,
								  C_dX3Limit,
								  C_dI1Limit,
								  C_dI2Limit,
								  C_dI3Limit,
								  C_dADPLimit,
								  C_dOptimal,
								  C_dTimeBoost,
								  C_dDeltaFlag,
								  C_dSmoothingValue,
								  calcProbaSurvOrNotId,
								  C_QBaseSmile,
                                  C_QForeignSmile,
								  C_spotVolOrFwdVol,
								  C_histoVolLongTerm,
								  convInVolFwdId,
								  C_calibBasisIncluded,
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

	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_PXL_ARM_TREE3FACT" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}



__declspec(dllexport) LPXLOPER WINAPI Local_ARM_Bootstrapping_VFDK_HW1To3F(LPXLOPER XL_volCurve,
                                                                           LPXLOPER XL_zc,
                                                                           LPXLOPER XL_noticeDates,
                                                                           LPXLOPER XL_swapStartDates,
                                                                           LPXLOPER XL_swapEndDates,
                                                                           LPXLOPER XL_HW3FParameters,
                                                                           LPXLOPER XL_observationDate)
{
	ADD_LOG("Local_ARM_Bootstrapping_VFDK_HW1To3F");
	// return
    LPXLOPER pxArray;

	static XLOPER XL_result;
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();

	// C variable 
	double C_observationDate;
	CCString C_volId;
    CCString C_zcId;
    VECTOR<double> C_noticeDates;
    VECTOR<double> C_swapStartDates;
    VECTOR<double> C_swapEndDates;
    VECTOR<double> C_HW3FParameters;


	// error
	static int error;
	static char* reason = "";

	

	XL_readStrCell(XL_volCurve,C_volId," ARM_ERR: volcurve id: object expected",C_result);
	XL_readStrCell(XL_zc,C_zcId," ARM_ERR: zc id: object expected",C_result);
	XL_readNumVector(XL_noticeDates,C_noticeDates," ARM_ERR: noticeDates: array of numeric expected",C_result);
    XL_readNumVector(XL_swapStartDates,C_swapStartDates," ARM_ERR: swapStartDates: array of numeric expected",C_result);
    XL_readNumVector(XL_swapEndDates,C_swapEndDates," ARM_ERR: swapEndDates: array of numeric expected",C_result);
    XL_readNumVector(XL_HW3FParameters,C_HW3FParameters," ARM_ERR: HW3FParameters: array of numeric expected",C_result);
    XL_readNumCell(XL_observationDate,C_observationDate," ARM_ERR: observation date: numeric expected",C_result);
	
	long retCode;
    VECTOR<double> matu;
	VECTOR<double> sigma;

	retCode = ARMLOCAL_Bootstrapping_VFDK_HW1To3F (LocalGetNumObjectId (C_volId),
								                  LocalGetNumObjectId (C_zcId),
								                  C_noticeDates,
                                                  C_swapStartDates,
                                                  C_swapEndDates,
                                                  C_HW3FParameters,
                                                  C_observationDate,
                                                  C_result);

	if(retCode == ARM_OK)
	{
		
		int nbrows;
		nbrows = C_noticeDates.size ();

		int nbcolumns = 2;
	
		FreeCurCellErr ();
		XL_result.xltype = xltypeMulti;
		XL_result.val.array.columns = nbcolumns;
		XL_result.val.array.rows = nbrows; 
		XL_result.val.array.lparray = pxArray = (LPXLOPER)GlobalAlloc (GMEM_ZEROINIT, nbrows * nbcolumns * sizeof (XLOPER));

		for(int i = 0; i < nbrows; i++)
		{
		
            //pxArray[XL_Coordonnate2Rank (i, 0, nbcolumns)].xltype |= xlbitDLLFree;
            pxArray[XL_Coordonnate2Rank (i, 0, nbcolumns)].xltype = xltypeNum;
			pxArray[XL_Coordonnate2Rank (i, 0, nbcolumns)].val.num = C_noticeDates[i];
			pxArray[XL_Coordonnate2Rank (i, 1, nbcolumns)].xltype = xltypeNum;
			pxArray[XL_Coordonnate2Rank (i, 1, nbcolumns)].val.num = C_result.getArray(i); 
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_PRCS3F_LatticePricing" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}

__declspec(dllexport) LPXLOPER WINAPI Local_ARM_SwaptionPrice_VFDK_HW1To3F(LPXLOPER XL_dSwaptionExpiryInYears,
                                                                           LPXLOPER XL_dSwaptionTenorInYears,
                                                                           LPXLOPER XL_dNoticePeriodInDays,
                                                                           LPXLOPER XL_dStrike,
                                                                           LPXLOPER XL_dCallPut,
                                                                           LPXLOPER XL_zc,
                                                                           LPXLOPER XL_noticeDates,
                                                                           LPXLOPER XL_sigma,
                                                                           LPXLOPER XL_HW3FParameters,
                                                                           LPXLOPER XL_observationDate)
{
	ADD_LOG("Local_ARM_SwaptionPrice_VFDK_HW1To3F");
	// return
	static XLOPER XL_result;
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();

	// C variable 
    double C_swaptionExpiryInYears;
    double C_swaptionTenorInYears;
    double C_noticePeriodInDays;
    double C_strike;
  	CCString C_optionType;
	long optionTypeId;
    double C_callput;
	double C_observationDate;
	CCString C_zcId;
    VECTOR<double> C_noticeDates;
    VECTOR<double> C_sigma;
    VECTOR<double> C_HW3FParameters;


	// error
	static int error;
	static char* reason = "";

	XL_readNumCell(XL_dSwaptionExpiryInYears,C_swaptionExpiryInYears," ARM_ERR: Swaption Expiry: numeric expected",C_result);
    XL_readNumCell(XL_dSwaptionTenorInYears,C_swaptionTenorInYears," ARM_ERR: Swaption Tenor: numeric expected",C_result);
    XL_readNumCell(XL_dNoticePeriodInDays,C_noticePeriodInDays," ARM_ERR: Notice Period: numeric expected",C_result);
    XL_readNumCell(XL_dStrike,C_strike," ARM_ERR: Strike: numeric expected",C_result);
    XL_readStrCell(XL_dCallPut,C_optionType," ARM_ERR: option type: string expected",C_result);
	XL_readStrCell(XL_zc,C_zcId," ARM_ERR: zc id: object expected",C_result);
	XL_readNumVector(XL_noticeDates,C_noticeDates," ARM_ERR: noticeDates: array of numeric expected",C_result);
    XL_readNumVector(XL_sigma,C_sigma," ARM_ERR: Sigma: array of numeric expected",C_result);
    XL_readNumVector(XL_HW3FParameters,C_HW3FParameters," ARM_ERR: HW3FParameters: array of numeric expected",C_result);
    XL_readNumCell(XL_observationDate,C_observationDate," ARM_ERR: observation date: numeric expected",C_result);
	
    if((optionTypeId = ARM_ConvCallOrPut (C_optionType, C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}
    if (optionTypeId==1)
        C_callput=0.0;
    else
        C_callput=1.0;

    long retCode;
    
    retCode = ARMLOCAL_SwaptionPrice_VFDK_HW1To3F ( C_swaptionExpiryInYears,
                                                    C_swaptionTenorInYears,
                                                    C_noticePeriodInDays,
                                                    C_strike,
                                                    C_callput,
                                                    LocalGetNumObjectId (C_zcId),
								                    C_noticeDates,
                                                    C_sigma,
                                                    C_HW3FParameters,
                                                    C_observationDate,
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_ARM_SwaptionPrice_VFDK_HW1To3F" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}
	

__declspec(dllexport) LPXLOPER WINAPI Local_ARM_ImpliedFwdCorrelation_VFDK_HW1To3F(LPXLOPER XL_dSwaptionExpiryInYears,
                                                                                   LPXLOPER XL_dSwaptionTenorInYears,
                                                                                   LPXLOPER XL_dSwaptionTenor2InYears,
                                                                                   LPXLOPER XL_dNoticePeriodInDays,
                                                                                   LPXLOPER XL_zc,
                                                                                   LPXLOPER XL_noticeDates,
                                                                                   LPXLOPER XL_sigma,
                                                                                   LPXLOPER XL_HW3FParameters,
                                                                                   LPXLOPER XL_observationDate)
{
	ADD_LOG("Local_ARM_ImpliedFwdCorrelation_VFDK_HW1To3F");
	// return
	static XLOPER XL_result;
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();

	// C variable 
    double C_swaptionExpiryInYears;
    double C_swaptionTenorInYears;
    double C_swaptionTenor2InYears;
    double C_noticePeriodInDays;
  	double C_observationDate;
	CCString C_zcId;
    VECTOR<double> C_noticeDates;
    VECTOR<double> C_sigma;
    VECTOR<double> C_HW3FParameters;


	// error
	static int error;
	static char* reason = "";

	XL_readNumCell(XL_dSwaptionExpiryInYears,C_swaptionExpiryInYears," ARM_ERR: Swaption Expiry: numeric expected",C_result);
    XL_readNumCell(XL_dSwaptionTenorInYears,C_swaptionTenorInYears," ARM_ERR: Swaption Tenor: numeric expected",C_result);
    XL_readNumCell(XL_dSwaptionTenor2InYears,C_swaptionTenor2InYears," ARM_ERR: Swaption Tenor: numeric expected",C_result);
    XL_readNumCell(XL_dNoticePeriodInDays,C_noticePeriodInDays," ARM_ERR: Notice Period: numeric expected",C_result);
	XL_readStrCell(XL_zc,C_zcId," ARM_ERR: zc id: object expected",C_result);
	XL_readNumVector(XL_noticeDates,C_noticeDates," ARM_ERR: noticeDates: array of numeric expected",C_result);
    XL_readNumVector(XL_sigma,C_sigma," ARM_ERR: Sigma: array of numeric expected",C_result);
    XL_readNumVector(XL_HW3FParameters,C_HW3FParameters," ARM_ERR: HW3FParameters: array of numeric expected",C_result);
    XL_readNumCell(XL_observationDate,C_observationDate," ARM_ERR: observation date: numeric expected",C_result);
	

    long retCode;
    
    retCode = ARMLOCAL_ImpliedFwdCorrelation_VFDK_HW1To3F ( C_swaptionExpiryInYears,
                                                            C_swaptionTenorInYears,
                                                            C_swaptionTenor2InYears,
                                                            C_noticePeriodInDays,
                                                            LocalGetNumObjectId (C_zcId),
                                                            C_noticeDates,
                                                            C_sigma,
                                                            C_HW3FParameters,
                                                            C_observationDate,
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_ARM_ImpliedFwdCorrelation_VFDK_HW1To3F" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}




__declspec(dllexport) LPXLOPER WINAPI Local_ARM_HW3F_CALIBRATION(LPXLOPER XL_dSpotDate,
                                                                 LPXLOPER XL_zc,
                                                                 LPXLOPER XL_HW3FParametersIn,
                                                                 LPXLOPER XL_volCurve,
                                                                 LPXLOPER XL_correlationCurve,
                                                                 LPXLOPER XL_volWeightsCurve,
                                                                 LPXLOPER XL_correlationWeightsCurve,
                                                                 LPXLOPER XL_dNoticeDates,
                                                                 LPXLOPER XL_dSwapStartDates,
                                                                 LPXLOPER XL_dSwapEndDates)
{
	ADD_LOG("Local_ARM_HW3F_CALIBRATION");
    // return
    LPXLOPER pxArray;

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();

	// C variable 
  	double C_spotDate;
	CCString C_zcId;
    CCString C_volId;
    CCString C_volWeightId;
    CCString C_correlationId;
    CCString C_correlationWeightId;
    VECTOR<double> C_noticeDates;
    VECTOR<double> C_swapStartDates;
    VECTOR<double> C_swapEndDates;;
    VECTOR<double> C_HW3FParametersIn;
    

	// error
	static int error;
	static char* reason = "";

	XL_readStrCell(XL_zc,C_zcId," ARM_ERR: zc id: object expected",C_result);
    XL_readStrCell(XL_volCurve,C_volId," ARM_ERR: volcurve id: object expected",C_result);
    
	XL_readStrCellWD(XL_correlationCurve,C_correlationId, "DEFAULT", " ARM_ERR: volcurve id: object expected",C_result);
    XL_readStrCellWD(XL_volWeightsCurve,C_volWeightId, "DEFAULT", " ARM_ERR: volcurve id: object expected",C_result);
    
	XL_readStrCellWD(XL_correlationWeightsCurve,C_correlationWeightId, "DEFAULT", " ARM_ERR: volcurve id: object expected",C_result);
	
	XL_readNumVector(XL_dNoticeDates,C_noticeDates," ARM_ERR: noticeDates: array of numeric expected",C_result);
    XL_readNumVector(XL_dSwapStartDates,C_swapStartDates," ARM_ERR: noticeDates: array of numeric expected",C_result);
    XL_readNumVector(XL_dSwapEndDates,C_swapEndDates," ARM_ERR: noticeDates: array of numeric expected",C_result);
    XL_readNumVector(XL_HW3FParametersIn,C_HW3FParametersIn," ARM_ERR: HW3FParameters: array of numeric expected",C_result);
    XL_readNumCell(XL_dSpotDate,C_spotDate," ARM_ERR: observation date: numeric expected",C_result);
	

    long retCode;
    
	long correlObjId;

    if ( C_correlationId == "DEFAULT" )
	{
	   correlObjId = ARM_NULL_OBJECT;
	}
	else
	{
	   correlObjId = LocalGetNumObjectId(C_correlationId);
	}

	long volWeightObjId;

    if ( C_volWeightId == "DEFAULT" )
	{
	   volWeightObjId = ARM_NULL_OBJECT;
	}
	else
	{
	   volWeightObjId = LocalGetNumObjectId(C_volWeightId);
	}

    long correlationWeightObjId;

    if ( C_correlationWeightId == "DEFAULT" )

	{
	   correlationWeightObjId = ARM_NULL_OBJECT;
	}
	else
	{
	   correlationWeightObjId = LocalGetNumObjectId(C_correlationWeightId);
	}

    retCode = ARMLOCAL_HW3F_CALIBRATION(C_spotDate,
                                        LocalGetNumObjectId(C_zcId),
                                        C_HW3FParametersIn,
                                        LocalGetNumObjectId(C_volId),
                                        correlObjId,
                                        volWeightObjId,
                                        correlationWeightObjId,
                                        C_noticeDates,
                                        C_swapStartDates,
                                        C_swapEndDates,
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


    if ( retCode == ARM_OK )
	{		
	   int nbrows;
	   nbrows = C_result.getLong();

	   int nbcolumns = 1;
	
	   FreeCurCellErr ();
	   XL_result.xltype = xltypeMulti;
	   XL_result.val.array.columns = nbcolumns;
	   XL_result.val.array.rows = nbrows; 
	   XL_result.val.array.lparray = pxArray = (LPXLOPER) GlobalAlloc(GMEM_ZEROINIT, nbrows*nbcolumns*sizeof(XLOPER));

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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_ARM_HW3F_CALIBRATION" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


/*--------------------------------------------------------------------------------------------*/
/*---- End Of File ----*/