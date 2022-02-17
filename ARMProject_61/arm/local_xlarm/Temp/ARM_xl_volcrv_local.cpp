#pragma warning(disable : 4786)
#pragma warning(disable : 4005)

#include <libCCxll\CCxll.h>

#include <ARM\libarm_local\ARM_local_class.h>
#include <ARM\libarm_local\ARM_local_glob.h>
#include <ARM\libarm_local\ARM_local_volcrv.h>

#include <ARM\libarm_frometk\ARM_local_parsexml.h>
#include <ARM\libarm_frometk\ARM_local_etoolkit.h>

#include "ARM_xl_wrapper_local.h"
#include "ARM_local_interglob.h"
#include "ARM_local_interface.h"
#include "ARM_xl_trycatch_local.h"
#include "ARM_xl_gp_fctorhelper.h"


#include "ExcelTools.h"

#include "util\tech_macro.h"

 

__declspec(dllexport) LPXLOPER WINAPI Local_VolCurv(LPXLOPER XL_matu,
													LPXLOPER XL_strikes,
													LPXLOPER XL_vols,
													LPXLOPER XL_date,
													LPXLOPER XL_volType,
													LPXLOPER XL_strikeType,
                                                    LPXLOPER XL_InterpolType,
                                                    LPXLOPER XL_ccy,
													LPXLOPER XL_indexName,
													LPXLOPER XL_indexId)
{
	ADD_LOG("Local_VolCurv");
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
	    VECTOR<double> C_strikes;
	    VECTOR<double> C_vols;
	    double C_date;

	    CCString C_volType;
	    long C_volTypeId;

	    CCString C_strikeType;
	    long C_strikeTypeId;

        CCString C_discountCcy;

        CCString C_StrInterpolType;

		CCString C_indexName;

        long C_InterpolType = K_LINEAR; // OR K_DIAG_INTERPOL

	    // error
	    static int error;
	    static char* reason = "";

	    XL_readNumVector(XL_matu,C_matu," ARM_ERR: maturities: array of numeric expected",C_result);
	    XL_readNumVector(XL_strikes,C_strikes," ARM_ERR: strikes: array of numeric expected",C_result);
	    XL_readNumVector(XL_vols,C_vols," ARM_ERR: volatilities: array of numeric expected",C_result);
	    XL_readNumCell(XL_date,C_date," ARM_ERR: date: date expected",C_result);
	    XL_readStrCellWD(XL_volType,C_volType,"A"," ARM_ERR: volatility type: string expected",C_result);
	    XL_readStrCellWD(XL_strikeType,C_strikeType,"Y"," ARM_ERR: strike type: string expected",C_result);
    
        XL_readStrCellWD(XL_InterpolType, C_StrInterpolType, "DEFAULT",
                           " ARM_ERR: Interpol. Type: (LINEAR, DIAG, SPLINE) expected", C_result);

		XL_readStrCellWD(XL_indexName, C_indexName, "",
                         " ARM_ERR: Index Name: string expected", C_result);

		std::string indexId; ExcelTools::convert(XL_indexId,"",indexId); 

        if ( C_StrInterpolType == "DIAG" )
        {
           C_InterpolType = K_DIAG_INTERPOL;       
        }
		else
		{
			if( C_StrInterpolType == "SPLINE" )
			{	C_InterpolType = K_SPLINE;  }
			else if( C_StrInterpolType == "FIX_RIGHT_MATU" ) //credit only
			{	C_InterpolType = 1001;  }					 //credit only			
			else if( C_StrInterpolType == "FIX_LEFT_MATU" )  //credit only
			{	C_InterpolType = 1002;  }					 //credit only			 
			else
			{	C_InterpolType = K_LINEAR;}
		}

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

	    if ( ( C_volTypeId = ARM_ConvVolType(C_volType, C_result) ) == ARM_DEFAULT_ERR )
	    {
	       ARM_ARG_ERR();

	       return (LPXLOPER)&XL_result;
	    }

	    if((C_strikeTypeId = StrikeCode (C_strikeType, C_result)) == ARM_DEFAULT_ERR)
	    {
		    ARM_ARG_ERR();
		    return (LPXLOPER)&XL_result;
	    }

	    long retCode;
	    long objId;
	    CCString prevClass;
	    long indexId_ = LocalPersistent::get().getObjectId(indexId);
	    CCString curClass = LOCAL_VOL_CURVE_LIN_CLASS;
	    CCString stringId = GetLastCurCellEnvValue ();
	    
	    if(!stringId)
	    {
		    retCode = ARMLOCAL_volcurv(C_matu, C_strikes, C_vols, C_date, C_strikeTypeId, 
                                       C_volTypeId,
                                       C_InterpolType,
                                       C_discountCcy,
									   C_indexName,
									   indexId_,
                                       C_result);

		    if ( retCode == ARM_OK )
		    {
			    objId = C_result.getLong();

			    LocalSetCurCellEnvValue(curClass, objId); 

			    stringId = LocalMakeObjectId(objId, curClass);
		    }
	    }
	    else
	    {
		    prevClass = LocalGetStringObjectClass (stringId);
		    
		    objId = LocalGetNumObjectId (stringId);
			    
		    if ( curClass == prevClass )
		    {
			    retCode = ARMLOCAL_volcurv(C_matu, C_strikes, C_vols, C_date, 
                                           C_strikeTypeId, C_volTypeId,
                                           C_InterpolType,
                                           C_discountCcy,
										   C_indexName,
										   indexId_,
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
			    
                retCode = ARMLOCAL_volcurv(C_matu, C_strikes, C_vols, C_date, 
                                           C_strikeTypeId, C_volTypeId,
                                           C_InterpolType,
                                           C_discountCcy,
										   C_indexName,
										   indexId_,
                                           C_result);
		    
			    if ( retCode == ARM_OK )
			    {
				    objId = C_result.getLong ();
			    
				    LocalSetCurCellEnvValue(curClass, objId); 

				    stringId = LocalMakeObjectId(objId, curClass);
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_VolCurv" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}



__declspec(dllexport) LPXLOPER WINAPI Local_ARM_FxVolCurv(LPXLOPER XL_date,
                                                          LPXLOPER XL_matus,
                                                          LPXLOPER XL_fxFwds,
													      LPXLOPER XL_pivotVols,
                                                          LPXLOPER XL_pivotTypes,
													      LPXLOPER XL_deltaPuts,
													      LPXLOPER XL_deltaCalls,
                                                          LPXLOPER XL_volsPuts,
                                                          LPXLOPER XL_volsCalls,
                                                          LPXLOPER XL_InterpolTypes,
													      LPXLOPER XL_WhatIsInterpolated,
													      LPXLOPER XL_correctSplineWithLinear,
                                                          LPXLOPER XL_fxSpot,
                                                          LPXLOPER XL_domZcCurve,
                                                          LPXLOPER XL_forZcCurve,
                                                          LPXLOPER XL_inRRSTR,
                                                          LPXLOPER XL_isATM)
{
	ADD_LOG("Local_ARM_FxVolCurv");
	static XLOPER XL_result;
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	    ARM_NOCALCIFWIZ();

	    // C variables

	    VECTOR<double> C_matu;
	    VECTOR<double> C_fxFwds;
        
        VECTOR<double> C_fxFwdsDefault;

        VECTOR<double> C_pivotVols;
        VECTOR<double> C_pivotTypes;
        VECTOR<double> C_deltaPuts;
        VECTOR<double> C_deltaCalls;
        VECTOR<double> C_volsPuts;
	    VECTOR<double> C_volsCalls;

	    double C_date;
	    double C_correctSplineWithLinear;
    
        CCString C_domZcCurve;
        CCString C_forZcCurve;

        long C_domZcCurveId;
        long C_forZcCurveId;

        long C_inRRSTR;
        CCString C_inRRSTRstring;
 
        long C_isATM;
        CCString C_isATMstring;

	    CCString C_StrWhatIsInterpolated;


        long C_InterpolType = K_LINEAR; // OR K_DPLINE
	    long C_WhatIsInterpolated = FXINTERP_STRIKE;

	    // error
	    static int error;
	    static char* reason = "";

        XL_readNumCell(XL_date, C_date," ARM_ERR: date: date expected",C_result);
	    XL_readNumVector(XL_matus, C_matu," ARM_ERR: maturities: array of numeric expected",C_result);
	    
        XL_readNumVectorWD(XL_fxFwds, C_fxFwds, C_fxFwdsDefault,
                           " ARM_ERR: FX fwds: array of numeric expected",C_result);
	    
        XL_readNumVector(XL_pivotVols, C_pivotVols," ARM_ERR: volatilities: array of numeric expected",C_result);
        XL_readNumVector(XL_pivotTypes, C_pivotTypes," ARM_ERR: : (0|1) array of numeric expected",C_result);
        XL_readNumVector(XL_deltaPuts, C_deltaPuts," ARM_ERR: deltas PUT: array of numeric expected",C_result);
        XL_readNumVector(XL_deltaCalls, C_deltaCalls," ARM_ERR: deltas CALL: array of numeric expected",C_result);
        XL_readNumVector(XL_volsPuts, C_volsPuts," ARM_ERR: PUT volatilities: array of numeric expected",C_result);
        XL_readNumVector(XL_volsCalls, C_volsCalls," ARM_ERR: CALL volatilities: array of numeric expected",C_result);
    
        VECTOR<CCString> C_StrInterpolTypes(C_matu.size());

	    VECTOR<CCString> C_InterpolDef(C_matu.size()); // Spline by default

        for (int h = 0; h < C_matu.size(); h++)
        {
            C_InterpolDef[h] = "S";
        }
    
        XL_readStrVectorWD(XL_InterpolTypes, C_StrInterpolTypes, C_InterpolDef,
              " ARM_ERR: flags: array of string expected", DOUBLE_TYPE, C_result);    
  
        VECTOR<double> C_InterpolTypes(C_matu.size());

        int i;
        for (i = 0; i < C_StrInterpolTypes.size(); i++)
        {
            long interpType; 

            if (( interpType = ARM_ConvInterpMethod(C_StrInterpolTypes[i], C_result) ) == ARM_DEFAULT_ERR )
            {
	           ARM_ARG_ERR();

	           return((LPXLOPER) &XL_result);
            }

            C_InterpolTypes[i] = interpType;
        }

        if ( C_StrInterpolTypes.size() != C_matu.size() )
        {
           for (int h  = i; h < C_matu.size(); h++)
           {
               C_InterpolTypes[h] = C_InterpolTypes[i-1];
           }
        }

	    XL_readStrCellWD(XL_WhatIsInterpolated, C_StrWhatIsInterpolated, "DEFAULT",
               " ARM_ERR: whatIsInterpoled. Type: S)TRIKE, D)ELTA ot E) Smooth Delta expected", C_result);
    
        if (( C_WhatIsInterpolated = ARM_ConvWhatIsInterp(C_StrWhatIsInterpolated, C_result)) == ARM_DEFAULT_ERR )
	    {
	       ARM_ARG_ERR();

	       return((LPXLOPER) &XL_result);
	    }

        CCString C_correctSplineWithLinearStr;

	    XL_readStrCell(XL_correctSplineWithLinear, C_correctSplineWithLinearStr,
                       " ARM_ERR: correctSplineWithLinear: Y/N or 1/0 expected",C_result);

        if (( C_correctSplineWithLinear = ARM_ConvYesOrNo(C_correctSplineWithLinearStr, C_result) ) == ARM_DEFAULT_ERR )
	    {
	       ARM_ARG_ERR();

	       return (LPXLOPER)&XL_result;
	    }

        double C_FxSpot    = 0.0;
        double C_FxSpotDef = 0.0;

        XL_readNumCellWD(XL_fxSpot, C_FxSpot, C_FxSpotDef, 
                       " ARM_ERR: FX spot: numeric expected",C_result);

        XL_readStrCellWD(XL_domZcCurve, C_domZcCurve, "DEFAULT", " ARM_ERR: Domestic ZC curve ID expected", C_result); 
        XL_readStrCellWD(XL_forZcCurve, C_forZcCurve, "DEFAULT", " ARM_ERR: Foreign ZC curve ID expected", C_result);

        if ( C_domZcCurve == "DEFAULT" )
	       C_domZcCurveId = ARM_NULL_OBJECT;
	    else
	       C_domZcCurveId = LocalGetNumObjectId(C_domZcCurve);

        if ( C_forZcCurve == "DEFAULT" )
	       C_forZcCurveId = ARM_NULL_OBJECT;
	    else
	       C_forZcCurveId = LocalGetNumObjectId(C_forZcCurve);

        XL_readStrCellWD(XL_inRRSTR, C_inRRSTRstring,
               "N"," ARM_ERR: (Y/N)Inputs are RR or STR : string expected", C_result);
       
        if (( C_inRRSTR = ARM_ConvYesOrNo(C_inRRSTRstring, C_result) ) == ARM_DEFAULT_ERR )
	    {
	       ARM_ARG_ERR();

	       return (LPXLOPER)&XL_result;
	    }

        XL_readStrCellWD(XL_isATM, C_isATMstring,
               "N"," ARM_ERR: (Y/N)Inputs are RR or STR : string expected", C_result);
       
        if (( C_isATM = ARM_ConvYesOrNo(C_isATMstring, C_result) ) == ARM_DEFAULT_ERR )
	    {
	       ARM_ARG_ERR();

	       return (LPXLOPER)&XL_result;
	    }

	    long retCode; 
	    long objId;
	    CCString prevClass;
	    
	    CCString curClass = LOCAL_VOLAT_FX_CLASS;
	    CCString stringId = GetLastCurCellEnvValue ();

    
	    if (!stringId)
	    {
           retCode = ARMLOCAL_FxVolcurv(C_date,
                                        C_matu,
					                    C_fxFwds,
                                        C_pivotVols,
                                        C_pivotTypes,
                                        C_deltaPuts,
                                        C_deltaCalls,
					                    C_volsPuts,
                                        C_volsCalls,
                                        C_InterpolTypes,
									    C_WhatIsInterpolated,
									    C_correctSplineWithLinear,
                                        C_FxSpot,
                                        C_domZcCurveId,
                                        C_forZcCurveId,
                                        C_inRRSTR,
                                        C_isATM,
					                    C_result);

		    if ( retCode == ARM_OK )
		    {
			    objId = C_result.getLong();

			    LocalSetCurCellEnvValue(curClass, objId); 

			    stringId = LocalMakeObjectId(objId, curClass);
		    }
	    }
	    else
	    {
		    prevClass = LocalGetStringObjectClass (stringId);
		    
		    objId = LocalGetNumObjectId (stringId);
			    
		    if ( curClass == prevClass )
		    {
			    retCode = ARMLOCAL_FxVolcurv(C_date,
                                             C_matu,
					                         C_fxFwds,
                                             C_pivotVols,
                                             C_pivotTypes,
                                             C_deltaPuts,
                                             C_deltaCalls,
					                         C_volsPuts,
                                             C_volsCalls,
                                             C_InterpolTypes,
										     C_WhatIsInterpolated,
										     C_correctSplineWithLinear,
                                             C_FxSpot,
                                             C_domZcCurveId,
                                             C_forZcCurveId,
                                             C_inRRSTR,
                                             C_isATM,
                                             C_result,
                                             objId);

			    if ( retCode == ARM_OK )
			    {
			       LocalSetCurCellEnvValue(curClass, objId); 

			       stringId = LocalMakeObjectId(objId, curClass);
			    }
		    }
		    else
		    {
			    FreeCurCellContent ();
			    
                retCode = ARMLOCAL_FxVolcurv(C_date,
                                             C_matu,
					                         C_fxFwds,
                                             C_pivotVols,
                                             C_pivotTypes,
                                             C_deltaPuts,
                                             C_deltaCalls,
					                         C_volsPuts,
                                             C_volsCalls,
                                             C_InterpolTypes,
										     C_WhatIsInterpolated,
										     C_correctSplineWithLinear,
                                             C_FxSpot,
                                             C_domZcCurveId,
                                             C_forZcCurveId,
                                             C_inRRSTR,
                                             C_isATM,
					                         C_result);
		    
			    if ( retCode == ARM_OK )
			    {
				    objId = C_result.getLong ();
			    
				    LocalSetCurCellEnvValue(curClass, objId); 

				    stringId = LocalMakeObjectId(objId, curClass);
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG("Unrecognized failure in Local_VolCurv")

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}



__declspec(dllexport) LPXLOPER WINAPI Local_PXL_ARM_FxVolCurv(LPXLOPER XL_date,
	                                                          LPXLOPER XL_matus,
	                                                          LPXLOPER XL_fxFwds,
														      LPXLOPER XL_pivotVols,
	                                                          LPXLOPER XL_pivotTypes,
														      LPXLOPER XL_deltaPuts,
														      LPXLOPER XL_deltaCalls,
	                                                          LPXLOPER XL_volsPuts,
	                                                          LPXLOPER XL_volsCalls,
	                                                          LPXLOPER XL_InterpolTypes,
														      LPXLOPER XL_WhatIsInterpolated,
														      LPXLOPER XL_correctSplineWithLinear,
	                                                          LPXLOPER XL_fxSpot,
	                                                          LPXLOPER XL_domZcCurve,
	                                                          LPXLOPER XL_forZcCurve,
	                                                          LPXLOPER XL_inRRSTR,
	                                                          LPXLOPER XL_isATM)
{
	ADD_LOG("Local_PXL_ARM_FxVolCurv");
	static XLOPER XL_result;
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	    ARM_NOCALCIFWIZ();

	    // C variables

	    VECTOR<double> C_matu;
	    VECTOR<double> C_fxFwds;
        
        VECTOR<double> C_fxFwdsDefault;

        VECTOR<double> C_pivotVols;
        VECTOR<double> C_pivotTypes;
        VECTOR<double> C_deltaPuts;
        VECTOR<double> C_deltaCalls;
        VECTOR<double> C_volsPuts;
	    VECTOR<double> C_volsCalls;

	    double C_date;
	    double C_correctSplineWithLinear;
    
        CCString C_domZcCurve;
        CCString C_forZcCurve;

        long C_domZcCurveId;
        long C_forZcCurveId;

        long C_inRRSTR;
        CCString C_inRRSTRstring;
 
        long C_isATM;
        CCString C_isATMstring;

	    CCString C_StrWhatIsInterpolated;


        long C_InterpolType = K_LINEAR; // OR K_DPLINE
	    long C_WhatIsInterpolated = FXINTERP_STRIKE;

	    // error
	    static int error;
	    static char* reason = "";

        XL_readNumCell(XL_date, C_date," ARM_ERR: date: date expected",C_result);
	    XL_readNumVector(XL_matus, C_matu," ARM_ERR: maturities: array of numeric expected",C_result);
	    
        XL_readNumVectorWD(XL_fxFwds, C_fxFwds, C_fxFwdsDefault,
                           " ARM_ERR: FX fwds: array of numeric expected",C_result);
	    
        XL_readNumVector(XL_pivotVols, C_pivotVols," ARM_ERR: volatilities: array of numeric expected",C_result);
        XL_readNumVector(XL_pivotTypes, C_pivotTypes," ARM_ERR: : (0|1) array of numeric expected",C_result);
        XL_readNumVector(XL_deltaPuts, C_deltaPuts," ARM_ERR: deltas PUT: array of numeric expected",C_result);
        XL_readNumVector(XL_deltaCalls, C_deltaCalls," ARM_ERR: deltas CALL: array of numeric expected",C_result);
        XL_readNumVector(XL_volsPuts, C_volsPuts," ARM_ERR: PUT volatilities: array of numeric expected",C_result);
        XL_readNumVector(XL_volsCalls, C_volsCalls," ARM_ERR: CALL volatilities: array of numeric expected",C_result);
    
        VECTOR<CCString> C_StrInterpolTypes(C_matu.size());

	    VECTOR<CCString> C_InterpolDef(C_matu.size()); // Spline by default

        for (int h = 0; h < C_matu.size(); h++)
        {
            C_InterpolDef[h] = "S";
        }
    
        XL_readStrVectorWD(XL_InterpolTypes, C_StrInterpolTypes, C_InterpolDef,
              " ARM_ERR: flags: array of string expected", DOUBLE_TYPE, C_result);    
  
        VECTOR<double> C_InterpolTypes(C_matu.size());

        int i;
        for (i = 0; i < C_StrInterpolTypes.size(); i++)
        {
            long interpType; 

            if (( interpType = ARM_ConvInterpMethod(C_StrInterpolTypes[i], C_result) ) == ARM_DEFAULT_ERR )
            {
	           ARM_ARG_ERR();

	           return((LPXLOPER) &XL_result);
            }

            C_InterpolTypes[i] = interpType;
        }

        if ( C_StrInterpolTypes.size() != C_matu.size() )
        {
           for (int h  = i; h < C_matu.size(); h++)
           {
               C_InterpolTypes[h] = C_InterpolTypes[i-1];
           }
        }

	    XL_readStrCellWD(XL_WhatIsInterpolated, C_StrWhatIsInterpolated, "DEFAULT",
               " ARM_ERR: whatIsInterpoled. Type: S)TRIKE, D)ELTA ot E) Smooth Delta expected", C_result);
    
        if (( C_WhatIsInterpolated = ARM_ConvWhatIsInterp(C_StrWhatIsInterpolated, C_result)) == ARM_DEFAULT_ERR )
	    {
	       ARM_ARG_ERR();

	       return((LPXLOPER) &XL_result);
	    }

        CCString C_correctSplineWithLinearStr;

	    XL_readStrCell(XL_correctSplineWithLinear, C_correctSplineWithLinearStr,
                       " ARM_ERR: correctSplineWithLinear: Y/N or 1/0 expected",C_result);

        if (( C_correctSplineWithLinear = ARM_ConvYesOrNo(C_correctSplineWithLinearStr, C_result) ) == ARM_DEFAULT_ERR )
	    {
	       ARM_ARG_ERR();

	       return (LPXLOPER)&XL_result;
	    }

        double C_FxSpot    = 0.0;
        double C_FxSpotDef = 0.0;

        XL_readNumCellWD(XL_fxSpot, C_FxSpot, C_FxSpotDef, 
                       " ARM_ERR: FX spot: numeric expected",C_result);

        XL_readStrCellWD(XL_domZcCurve, C_domZcCurve, "DEFAULT", " ARM_ERR: Domestic ZC curve ID expected", C_result); 
        XL_readStrCellWD(XL_forZcCurve, C_forZcCurve, "DEFAULT", " ARM_ERR: Foreign ZC curve ID expected", C_result);

        if ( C_domZcCurve == "DEFAULT" )
	       C_domZcCurveId = ARM_NULL_OBJECT;
	    else
	       C_domZcCurveId = LocalGetNumObjectId(C_domZcCurve);

        if ( C_forZcCurve == "DEFAULT" )
	       C_forZcCurveId = ARM_NULL_OBJECT;
	    else
	       C_forZcCurveId = LocalGetNumObjectId(C_forZcCurve);

        XL_readStrCellWD(XL_inRRSTR, C_inRRSTRstring,
               "N"," ARM_ERR: (Y/N)Inputs are RR or STR : string expected", C_result);
       
        if (( C_inRRSTR = ARM_ConvYesOrNo(C_inRRSTRstring, C_result) ) == ARM_DEFAULT_ERR )
	    {
	       ARM_ARG_ERR();

	       return (LPXLOPER)&XL_result;
	    }

        XL_readStrCellWD(XL_isATM, C_isATMstring,
               "N"," ARM_ERR: (Y/N)Inputs are RR or STR : string expected", C_result);
       
        if (( C_isATM = ARM_ConvYesOrNo(C_isATMstring, C_result) ) == ARM_DEFAULT_ERR )
	    {
	       ARM_ARG_ERR();

	       return (LPXLOPER)&XL_result;
	    }

	    long retCode; 
	    long objId;
	    
	    CCString curClass = LOCAL_VOLAT_FX_CLASS;
	    CCString stringId = GetLastCurCellEnvValue ();

    
        retCode = ARMLOCAL_FxVolcurv(C_date,
                                     C_matu,
				                     C_fxFwds,
                                     C_pivotVols,
                                     C_pivotTypes,
                                     C_deltaPuts,
                                     C_deltaCalls,
				                     C_volsPuts,
                                     C_volsCalls,
                                     C_InterpolTypes,
								     C_WhatIsInterpolated,
								     C_correctSplineWithLinear,
                                     C_FxSpot,
                                     C_domZcCurveId,
                                     C_forZcCurveId,
                                     C_inRRSTR,
                                     C_isATM,
				                     C_result);

	    if ( retCode == ARM_OK )
	    {
		    objId = C_result.getLong();

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
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG("Unrecognized failure in Local_VolCurv")

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}



__declspec(dllexport) LPXLOPER WINAPI Local_ARM_NewComputeFxVolatility(LPXLOPER XL_curve,
															           LPXLOPER XL_matu,
															           LPXLOPER XL_moneyness)
{
	ADD_LOG("Local_ARM_NewComputeFxVolatility");
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
		double C_moneyness;
	

		// error
		static int error;
		static char* reason = "";

		XL_readStrCell(XL_curve, C_curve," ARM_ERR: string expected", C_result);
		XL_readNumCell(XL_matu, C_matu," ARM_ERR: maturity: numeric expected", C_result);
		XL_readNumCell(XL_moneyness, C_moneyness," ARM_ERR: moneyness: numeric expected", C_result);


		long retCode = ARMLOCAL_NewComputeFxVolatility(LocalGetNumObjectId (C_curve), 
                                                       C_moneyness,
												       C_matu, C_result);

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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_NewComputeFXVolatility" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}



int MIC_setNumMatrixSize(XLOPER& XL_result,
                         double* values,
                         long nbRows, long nbCols )
{
	/// function to free an xloper
	if ( XLOPER_Holder::FreeXLOPER_Multi(XL_result) == XL_ERROR )
	   return XL_ERROR;

	LPXLOPER pxArray;
	

    if (!(XL_result.val.array.lparray = 
             pxArray = 
             (LPXLOPER) GlobalAlloc(GMEM_ZEROINIT, nbRows*nbCols*sizeof(XLOPER))))
	{
		/// problem with memory allocation!
		return XL_ERROR; 
	}

	XL_result.xltype			= xltypeMulti;
	XL_result.val.array.columns = nbCols;
	XL_result.val.array.rows	= nbRows; 

	int rowIdx, colIdx;
    int idx = 0;

	for (rowIdx = 0; rowIdx < nbRows; ++rowIdx)
	{
	    for (colIdx = 0; colIdx < nbCols; ++colIdx)
	    {
		    pxArray[XL_Coordonnate2Rank(rowIdx, colIdx, nbCols)].val.num = values[idx];
		    pxArray[XL_Coordonnate2Rank(rowIdx, colIdx, nbCols)].xltype = xltypeNum;
			
			++idx;
        }
	}

	return(XL_NO_ERROR);
}



__declspec(dllexport) LPXLOPER WINAPI Local_ARM_GetInfoFromFxVolatility(LPXLOPER XL_curve)
{
	ADD_LOG("Local_ARM_GetInfoFromFxVolatility");
	static XLOPER XL_result;
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
		ARM_NOCALCIFWIZ();
		
        // C variable
		CCString C_curve;
		int nbColumns, nbRows;
		double* res = NULL;
	

		// error
		static int error;
		static char* reason = "";

		XL_readStrCell(XL_curve, C_curve," ARM_ERR: string expected",C_result);


		long retCode = ARMLOCAL_GetInfoFromFxVolatility(LocalGetNumObjectId(C_curve), 
                                                        res, 
                                                        nbColumns, nbRows, 
														C_result);

		
		if ( retCode == ARM_OK )
		{
           FreeCurCellErr();

           MIC_setNumMatrixSize(XL_result, res, nbRows, nbColumns);

           if (res)
              delete [] res;
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_GetInfoFromFXVolatility" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}



__declspec(dllexport) LPXLOPER WINAPI Local_PXL_VolCurv(LPXLOPER XL_matu,
														LPXLOPER XL_strikes,
														LPXLOPER XL_vols,
														LPXLOPER XL_date,
														LPXLOPER XL_volType,
														LPXLOPER XL_strikeType,
                                                        LPXLOPER XL_InterpolType,
                                                        LPXLOPER XL_ccy,
														LPXLOPER XL_indexName,
														LPXLOPER XL_indexId)
{
	ADD_LOG("Local_PXL_VolCurv");
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
	    VECTOR<double> C_strikes;
	    VECTOR<double> C_vols;
	    double C_date;

	    CCString C_volType;
	    long C_volTypeId;

	    CCString C_strikeType;
	    long C_strikeTypeId;

        CCString C_discountCcy;

        CCString C_StrInterpolType;

		CCString C_indexName;

        long C_InterpolType = K_LINEAR; // OR K_DIAG_INTERPOL

	    // error
	    static int error;
	    static char* reason = "";

	    XL_readNumVector(XL_matu,C_matu," ARM_ERR: maturities: array of numeric expected",C_result);
	    XL_readNumVector(XL_strikes,C_strikes," ARM_ERR: strikes: array of numeric expected",C_result);
	    XL_readNumVector(XL_vols,C_vols," ARM_ERR: volatilities: array of numeric expected",C_result);
	    XL_readNumCell(XL_date,C_date," ARM_ERR: date: date expected",C_result);
	    XL_readStrCellWD(XL_volType,C_volType,"A"," ARM_ERR: volatility type: string expected",C_result);
	    XL_readStrCellWD(XL_strikeType,C_strikeType,"Y"," ARM_ERR: strike type: string expected",C_result);

        XL_readStrCellWD(XL_InterpolType, C_StrInterpolType, "DEFAULT",
                           " ARM_ERR: Interpol. Type: (LINEAR, DIAG) expected", C_result);

		XL_readStrCellWD(XL_indexName, C_indexName, "",
                         " ARM_ERR: Index Name: string expected", C_result);
		std::string indexId; ExcelTools::convert(XL_indexId,"",indexId); 

        if ( C_StrInterpolType == "DIAG" )
        {
           C_InterpolType = K_DIAG_INTERPOL;       
        }
        else
        {
           C_InterpolType = K_LINEAR;
        }

        XL_readStrCellWD(XL_ccy, C_discountCcy,"DEFAULT",
                             " ARM_ERR: currency: string expected",C_result);

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

	    if((C_volTypeId = ARM_ConvVolType (C_volType, C_result)) == ARM_DEFAULT_ERR)
	    {
		    ARM_ARG_ERR();
		    return (LPXLOPER)&XL_result;
	    }

	    if((C_strikeTypeId = StrikeCode (C_strikeType, C_result)) == ARM_DEFAULT_ERR)
	    {
		    ARM_ARG_ERR();
		    return (LPXLOPER)&XL_result;
	    }

	    long retCode;
	    long objId;
	    
	    CCString curClass = LOCAL_VOL_CURVE_LIN_CLASS;
	    CCString stringId;
		long indexId_ = LocalPersistent::get().getObjectId(indexId);
	    retCode = ARMLOCAL_volcurv(C_matu, C_strikes, C_vols, C_date, C_volTypeId, 
                                   C_strikeTypeId,
                                   C_InterpolType,
                                   C_discountCcy,
								   C_indexName,
								   indexId_,
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_PXL_VolCurv" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}



__declspec(dllexport) LPXLOPER WINAPI Local_VolFlat(LPXLOPER XL_vol,
													LPXLOPER XL_date,
                                                    LPXLOPER XL_ccy)
{
	ADD_LOG("Local_VolFlat");
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	    ARM_NOCALCIFWIZ();


	    // C variable
	    double C_vol;
	    double C_date;

        CCString C_discountCcy;

	    // error
	    static int error;
	    static char* reason = "";

	    XL_readNumCell(XL_vol,C_vol," ARM_ERR: volatility: numeric expected",C_result);
	    XL_readNumCell(XL_date,C_date," ARM_ERR: date: date expected",C_result);
	    
        XL_readStrCellWD(XL_ccy, C_discountCcy, "DEFAULT",
                             " ARM_ERR: currency: string expected",C_result);

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
	    
	    CCString curClass = LOCAL_VOL_FLAT_CLASS;
	    CCString stringId = GetLastCurCellEnvValue ();
	    
	    if (!stringId)
	    {
		    retCode = ARMLOCAL_volflat(C_vol, C_date, 
                                       C_discountCcy,
                                       C_result);

		    if ( retCode == ARM_OK )
		    {
			    objId = C_result.getLong ();

			    LocalSetCurCellEnvValue(curClass, objId); 

			    stringId = LocalMakeObjectId(objId, curClass);
		    }
	    }
	    else
	    {
		    prevClass = LocalGetStringObjectClass (stringId);
		    
		    objId = LocalGetNumObjectId(stringId);
			    
		    if ( curClass == prevClass )
		    {
			    retCode = ARMLOCAL_volflat(C_vol, C_date, 
                                           C_discountCcy,
                                           C_result, objId);

			    if (retCode == ARM_OK)
			    {
				    LocalSetCurCellEnvValue (curClass, objId); 

				    stringId = LocalMakeObjectId (objId, curClass);
			    }
		    }
		    else
		    {
			    FreeCurCellContent();
			    
                retCode = ARMLOCAL_volflat(C_vol, C_date, 
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_VolFlat" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}



__declspec(dllexport) LPXLOPER WINAPI Local_PXL_VolFlat(LPXLOPER XL_vol,
														LPXLOPER XL_date,
                                                        LPXLOPER XL_ccy)
{
	ADD_LOG("Local_PXL_VolFlat");
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	    ARM_NOCALCIFWIZ();


	    // C variable
	    double C_vol;
	    double C_date;

        CCString C_discountCcy;

	    // error
	    static int error;
	    static char* reason = "";

	    XL_readNumCell(XL_vol,C_vol," ARM_ERR: volatility: numeric expected",C_result);
	    XL_readNumCell(XL_date,C_date," ARM_ERR: date: date expected",C_result);
	    
        XL_readStrCellWD(XL_ccy, C_discountCcy, "DEFAULT",
                             " ARM_ERR: currency: string expected",C_result);

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
	    
	    CCString curClass = LOCAL_VOL_FLAT_CLASS;
	    CCString stringId;
	    
	    retCode = ARMLOCAL_volflat(C_vol, C_date, 
                                   C_discountCcy,
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_PXL_VolFlat" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}



__declspec(dllexport) LPXLOPER WINAPI Local_VolCube(LPXLOPER XL_ATMVolId,
													LPXLOPER XL_volCurves,
													LPXLOPER XL_tenors,
                                                    LPXLOPER XL_VolType,
													LPXLOPER XL_CheckCcy)
{
	ADD_LOG("Local_VolCube");
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

    /// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	    ARM_NOCALCIFWIZ();


	    // C variable
	    CCString C_ATMVolId;
	    
	    VECTOR<CCString>	C_volCurves;
	    VECTOR<long>		C_volCurveIds;

	    VECTOR<double>		C_tenors;

		double C_checkCcy;
		double default_CheckCcy = 1.0;


        CCString C_volType;
        long C_volTypeId;

	    // error
	    static int error;
	    static char* reason = "";

	    XL_readStrCell(XL_ATMVolId,C_ATMVolId," ARM_ERR: ATM volatility curve : string expected",C_result);
	    XL_readStrVector(XL_volCurves,C_volCurves," ARM_ERR: vol curves : array of string expected",XL_TYPE_STRING,C_result);
	    XL_readNumVector(XL_tenors,C_tenors," ARM_ERR: tenors : array of numeric expected",C_result);
	    XL_readStrCellWD(XL_VolType,C_volType,"A"," ARM_ERR: volatility type: string expected",C_result);
		XL_readNumCellWD(XL_CheckCcy,C_checkCcy,default_CheckCcy," ARM_ERR: CheckCcy: num expected",C_result);

	    if((C_volTypeId = ARM_ConvVolType(C_volType, C_result)) == ARM_DEFAULT_ERR)
	    {
		    ARM_ARG_ERR();
		    return (LPXLOPER)&XL_result;
	    }

		if (C_volCurves.size() != C_tenors.size())
			throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"number of tenors must be equal to number of curves");

	    for (int i = 0; i < C_volCurves.size(); i++)
		    C_volCurveIds.push_back(LocalGetNumObjectId(C_volCurves[i]));

	    CCString prevClass;
	    long retCode;
	    long objId;

	    CCString curClass = LOCAL_VOL_CUBE_CLASS;
	    CCString stringId = GetLastCurCellEnvValue ();
	    
	    if (!stringId)
	    {
		    retCode = ARMLOCAL_VolCube(LocalGetNumObjectId(C_ATMVolId), C_volCurveIds, 
                                       C_tenors, C_volTypeId, (long)C_checkCcy, C_result);

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
			    retCode = ARMLOCAL_VolCube(LocalGetNumObjectId(C_ATMVolId), C_volCurveIds,
                                           C_tenors, C_volTypeId, (long)C_checkCcy,
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
			    retCode = ARMLOCAL_VolCube(LocalGetNumObjectId(C_ATMVolId), C_volCurveIds,
                                           C_tenors,
                                           C_volTypeId, (long)C_checkCcy, C_result);
		    
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_VolCube" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}



__declspec(dllexport) LPXLOPER WINAPI Local_PXL_VolCube (LPXLOPER XL_ATMVolId,
														 LPXLOPER XL_volCurves,
														 LPXLOPER XL_tenors,
                                                         LPXLOPER XL_VolType,
														 LPXLOPER XL_CheckCcy)
{
	ADD_LOG("Local_PXL_VolCube ");
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	    ARM_NOCALCIFWIZ();


	    // C variable
	    CCString C_ATMVolId;
	    
	    VECTOR<CCString>	C_volCurves;
	    VECTOR<long>		C_volCurveIds;

	    VECTOR<double>		C_tenors;
        CCString            C_volType;
        long                C_volTypeId;
		double				C_checkCcy;
		double				default_checkCcy = 1.0;

	    // error
	    static int error;
	    static char* reason = "";

	    XL_readStrCell(XL_ATMVolId,C_ATMVolId," ARM_ERR: ATM volatility curve : string expected",C_result);
	    XL_readStrVector(XL_volCurves,C_volCurves," ARM_ERR: vol curves : array of string expected",XL_TYPE_STRING,C_result);
	    XL_readNumVector(XL_tenors,C_tenors," ARM_ERR: tenors : array of numeric expected",C_result);
	    XL_readStrCellWD(XL_VolType,C_volType,"A"," ARM_ERR: volatility type: string expected",C_result);
		XL_readNumCellWD(XL_CheckCcy,C_checkCcy,default_checkCcy," ARM_ERR: CheckCcy: num expected",C_result);

	    if((C_volTypeId = ARM_ConvVolType (C_volType, C_result)) == ARM_DEFAULT_ERR)
	    {
		    ARM_ARG_ERR();
		    return (LPXLOPER)&XL_result;
	    }

	    for (int i=0; i<C_volCurves.size();i++)
		    C_volCurveIds.push_back(LocalGetNumObjectId(C_volCurves[i]));

	    long retCode;
	    long objId;

	    CCString curClass = LOCAL_VOL_CUBE_CLASS;
	    CCString stringId;

	    retCode = ARMLOCAL_VolCube(LocalGetNumObjectId(C_ATMVolId), C_volCurveIds, 
                                   C_tenors, C_volTypeId, (long)C_checkCcy, C_result);

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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_PXL_VolCube" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}



__declspec(dllexport) LPXLOPER WINAPI Local_ComputeVolatility (LPXLOPER XL_curve,
															   LPXLOPER XL_matu,
															   LPXLOPER XL_strike,
															   LPXLOPER XL_tenor)
{
	ADD_LOG("Local_ComputeVolatility ");
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
		double C_strike;
		double C_tenor;
		double default_tenor = 0.0;

		// error
		static int error;
		static char* reason = "";

		XL_readStrCell(XL_curve,C_curve," ARM_ERR: string expected",C_result);
		XL_readNumCell(XL_matu,C_matu," ARM_ERR: maturity: numeric expected",C_result);
		XL_readNumCell(XL_strike,C_strike," ARM_ERR: strike: numeric expected",C_result);
		XL_readNumCellWD(XL_tenor, C_tenor, default_tenor,
							" ARM_ERR: tenor: numeric expected",C_result);

		long retCode = ARMLOCAL_ComputeVolatility (LocalGetNumObjectId (C_curve), C_matu,
													C_strike, C_tenor, C_result);

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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_ComputeVolatility" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}



__declspec(dllexport) LPXLOPER WINAPI Local_ComputeModelVolatility (LPXLOPER XL_model,
																	LPXLOPER XL_matu,
																	LPXLOPER XL_tenor,
																	LPXLOPER XL_fwd,
																	LPXLOPER XL_strike,
																	LPXLOPER XL_useSabr)
{
	ADD_LOG("Local_ComputeModelVolatility ");
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	    ARM_NOCALCIFWIZ();


	    // C variable
	    CCString C_model;
	    double C_matu;
	    double C_strike;
	    double C_tenor;
	    double C_fwd;
		CCString C_useSabr;

	    // error
	    static int error;
	    static char* reason = "";

	    XL_readStrCell(XL_model,C_model," ARM_ERR: string expected",C_result);
	    XL_readNumCell(XL_matu,C_matu," ARM_ERR: maturity: numeric expected",C_result);
	    XL_readNumCell(XL_strike,C_strike," ARM_ERR: strike: numeric expected",C_result);
	    XL_readNumCell(XL_tenor, C_tenor," ARM_ERR: tenor: numeric expected",C_result);
	    XL_readNumCell(XL_fwd, C_fwd," ARM_ERR: tenor: numeric expected",C_result);
	    XL_readStrCellWD(XL_useSabr, C_useSabr, "N", " ARM_ERR: use SABR: string expected (Y/N)",C_result);

		int useSabr = (C_useSabr == "Y" || C_useSabr == "YES") ? K_YES : K_NO;

	    long retCode = ARMLOCAL_ComputeModelVolatility (LocalGetNumObjectId (C_model),
													    C_matu,
													    C_tenor,
													    C_fwd,
													    C_strike,
													    C_result,
														useSabr);

	    if ( retCode == ARM_OK )
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_ComputeModelVolatility" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}




__declspec(dllexport) LPXLOPER WINAPI Local_ComputeCorrel(LPXLOPER XL_Correl,
	                                                      LPXLOPER XL_X,
	                                                      LPXLOPER XL_Y)
{
	ADD_LOG("Local_ComputeCorrel");
	/// use a dummy XL_X since not used
	return Local_ComputeVolatility( XL_Correl, XL_X, XL_Y, XL_X );
}


__declspec(dllexport) LPXLOPER WINAPI Local_GetVolFromCalypso(LPXLOPER XL_index,
															 LPXLOPER XL_currency,
															 LPXLOPER XL_cvName,
															 LPXLOPER XL_date,
															 LPXLOPER XL_vtype,
															 LPXLOPER XL_matIndex,
															 LPXLOPER XL_impOrHist,
															 LPXLOPER XL_indexId)
{
	ADD_LOG("Local_GetVolFromCalypso");
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	    ARM_NOCALCIFWIZ();


	    // C variable
	    std::string C_index;	ExcelTools::convert(XL_index,"",C_index);
	    std::string C_currency; ExcelTools::convert(XL_currency,"",C_currency);
	    std::string C_cvName;	ExcelTools::convert(XL_cvName,"",C_cvName);
	    ARM_Date C_date;		ExcelTools::convert(XL_date,C_date);
	    std::string C_vtype;	ExcelTools::convert(XL_vtype,"",C_vtype);
	    std::string C_matIndex;	ExcelTools::convert(XL_matIndex,"ATM",C_matIndex); 
	    std::string C_impOrHist;ExcelTools::convert(XL_impOrHist,"IRFWDVOL",C_impOrHist);
	    
		std::string indexId; ExcelTools::convert(XL_indexId,"",indexId); 
		long indexId_ = LocalPersistent::get().getObjectId(indexId); 

	    
		_strupr((char *) C_index.c_str());
		_strupr((char *) C_currency.c_str());
		_strupr((char *) C_cvName.c_str());
		_strupr((char *) C_matIndex.c_str());
		_strupr((char *) C_impOrHist.c_str());

	    
	    
	    CCString prevClass;
	    
	    CCString curClass = LOCAL_VOL_CURVE_LIN_CLASS;
	    CCString stringId = GetLastCurCellEnvValue ();

        long objId=-1;
		if (!(!stringId)){
			prevClass = LocalGetStringObjectClass (stringId);
			if ( curClass == prevClass ){
				objId = LocalGetNumObjectId (stringId);
			}else{
				FreeCurCellContent();
			}
		}
		long retCode = ARMLOCAL_GetVolFromCalypso (C_index,
											     C_currency,
											     C_cvName,
											     C_date,
											     C_vtype,
											     C_matIndex,
											     C_impOrHist,
												 indexId_,
											     C_result,
												 objId);


	    
        if ( retCode == ARM_OK )
	    {
		   	std::string objectLabel =ExcelCaller::get().setObject(C_result.getLong(), LOCAL_VOL_CURVE_LIN_CLASS) ;
		    // ExcelCaller::get().setObject(objectLabel) ;
		    ExcelTools::convert(objectLabel,&XL_result); 
		    return &XL_result;		
	    }
	    else
	    {
		    ARM_ERR();
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
	return &XL_result;
}

__declspec(dllexport) LPXLOPER WINAPI Local_GetVolFromSummit(LPXLOPER XL_index,
															 LPXLOPER XL_currency,
															 LPXLOPER XL_cvName,
															 LPXLOPER XL_date,
															 LPXLOPER XL_vtype,
															 LPXLOPER XL_matIndex,
															 LPXLOPER XL_impOrHist,
															 LPXLOPER XL_indexId)
{
	ADD_LOG("Local_GetVolFromSummit");
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
	    double C_date;
	    CCString C_vtype;
	    CCString C_matIndex;
	    CCString C_impOrHist;

	    // error
	    static int error;
	    static char* reason = "";

	    XL_readStrCell(XL_index,C_index," ARM_ERR: index: string expected",C_result);
	    XL_readStrCell(XL_currency,C_currency," ARM_ERR: currency: string expected",C_result);
	    XL_readStrCell(XL_cvName,C_cvName," ARM_ERR: curve name: string expected",C_result);
	    XL_readNumCell(XL_date,C_date," ARM_ERR: as of date: date expected",C_result);
	    XL_readStrCell(XL_vtype,C_vtype," ARM_ERR: volatility type: string expected",C_result);
	    XL_readStrCellWD(XL_matIndex,C_matIndex,"ATM"," ARM_ERR: maturity index: string expected",C_result);	
	    XL_readStrCellWD(XL_impOrHist,C_impOrHist,"IRFWDVOL"," ARM_ERR: implied or Historical?: string expected",C_result);

		std::string indexId; ExcelTools::convert(XL_indexId,"",indexId); 
		long indexId_ = LocalPersistent::get().getObjectId(indexId); 

	    C_index.toUpper ();
	    C_currency.toUpper ();
	    C_cvName.toUpper ();
	    C_matIndex.toUpper ();
	    C_impOrHist.toUpper ();

	    long retCode;
	    long objId;
	    CCString prevClass;
	    
	    CCString curClass = LOCAL_VOL_CURVE_LIN_CLASS;
	    CCString stringId = GetLastCurCellEnvValue ();
	    
	    if (!stringId)
	    {
		    retCode = ARMLOCAL_GetVolFromSummit (C_index,
											     C_currency,
											     C_cvName,
											     C_date,
											     C_vtype,
											     C_matIndex,
											     C_impOrHist,
												 indexId_,
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
			    retCode = ARMLOCAL_GetVolFromSummit (C_index,
												     C_currency,
												     C_cvName,
												     C_date,
												     C_vtype,
												     C_matIndex,
												     C_impOrHist,
												     indexId_,
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
			    retCode = ARMLOCAL_GetVolFromSummit (C_index,
												     C_currency,
												     C_cvName,
												     C_date,
												     C_vtype,
												     C_matIndex,
												     C_impOrHist,
													 indexId_,
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in GetVolFromSummit" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}



__declspec(dllexport) LPXLOPER WINAPI Local_PXL_GetVolFromSummit(LPXLOPER XL_index,
																 LPXLOPER XL_currency,
																 LPXLOPER XL_cvName,
																 LPXLOPER XL_date,
																 LPXLOPER XL_vtype,
																 LPXLOPER XL_matIndex,
																 LPXLOPER XL_impOrHist,
																 LPXLOPER XL_IndexId)
{
	ADD_LOG("Local_PXL_GetVolFromSummit");
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
	    double C_date;
	    CCString C_vtype;
	    CCString C_matIndex;
	    CCString C_impOrHist;

	    // error
	    static int error;
	    static char* reason = "";

	    XL_readStrCell(XL_index,C_index," ARM_ERR: index: string expected",C_result);
	    XL_readStrCell(XL_currency,C_currency," ARM_ERR: currency: string expected",C_result);
	    XL_readStrCell(XL_cvName,C_cvName," ARM_ERR: curve name: string expected",C_result);
	    XL_readNumCell(XL_date,C_date," ARM_ERR: as of date: date expected",C_result);
	    XL_readStrCell(XL_vtype,C_vtype," ARM_ERR: volatility type: string expected",C_result);
	    XL_readStrCellWD(XL_matIndex,C_matIndex,"ATM"," ARM_ERR: maturity index: string expected",C_result);	
	    XL_readStrCellWD(XL_impOrHist,C_impOrHist,"IRFWDVOL"," ARM_ERR: implied or Historical?: string expected",C_result);

		std::string indexId; ExcelTools::convert(XL_IndexId,"",indexId); 
		long indexId_=LocalPersistent::get().getObjectId(indexId); 

	    C_index.toUpper ();
	    C_currency.toUpper ();
	    C_cvName.toUpper ();
	    C_matIndex.toUpper ();
	    C_impOrHist.toUpper ();

	    long retCode;
	    long objId;
	    
	    CCString curClass = LOCAL_VOL_CURVE_LIN_CLASS;
	    CCString stringId;
	    
	    retCode = ARMLOCAL_GetVolFromSummit(C_index,
										    C_currency,
										    C_cvName,
										    C_date,
										    C_vtype,
										    C_matIndex,
										    C_impOrHist,
											indexId_,
										    C_result);

	    if ( retCode == ARM_OK )
	    {
		    objId = C_result.getLong ();

		    stringId = LocalMakeObjectId (objId, curClass);
	    }
	    
	    if ( retCode == ARM_OK )
	    {
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_PXL_GetVolFromSummit" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}



__declspec(dllexport) LPXLOPER WINAPI Local_GetVolCubeFromCalypso(LPXLOPER XL_index,
																 LPXLOPER XL_currency,
																 LPXLOPER XL_cvName,
																 LPXLOPER XL_date,
																 LPXLOPER XL_vtype,
                                                                 LPXLOPER XL_suffix,
                                                                 LPXLOPER XL_type,
																 LPXLOPER XL_tenors,
																 LPXLOPER XL_SmileOrNot,
																 LPXLOPER XL_indexid)
{
	ADD_LOG("Local_GetVolCubeFromCalypso");
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

    /// to make the infrastructure very robust we put a try catch!
	try
	{
	    ARM_NOCALCIFWIZ();


	    // C variable
		std::string C_index;     ExcelTools::convert(XL_index,"",C_index);
	    std::string C_currency;	 ExcelTools::convert(XL_currency,"",C_currency);
	    std::string C_cvName;	 ExcelTools::convert(XL_cvName,"",C_cvName);
	    ARM_Date C_date;		 ExcelTools::convert(XL_date,C_date);
	    std::string C_vtype;	 ExcelTools::convert(XL_vtype,"",C_vtype);
        std::string C_type;	 ExcelTools::convert(XL_type,"",C_type);
        std::string C_suffix;	 ExcelTools::convert(XL_suffix,"",C_suffix);


	    VECTOR<string> C_tenors;		 ; // ExcelTools::convert(XL_vtype,C_tenors);
	    VECTOR<string> C_tenors_default; ExcelTools::convert(XL_vtype,C_tenors_default,C_tenors);

		std::string C_SmileOrNot;	 ExcelTools::convert(XL_SmileOrNot,"",C_SmileOrNot) ;
	    // to remove if not necessary
		//long C_SmileOrNotId;	    
	    
		std::string indexId;	 ExcelTools::convert(XL_indexid,"",indexId); 
		long indexId_= LocalPersistent::get().getObjectId(indexId) ; 
	    

		
		_strupr((char *) C_index.c_str());
		_strupr((char *) C_currency.c_str());
		_strupr((char *) C_cvName.c_str());

	    for (int i = 0; i < C_tenors.size(); i++)
		   _strupr((char *) C_tenors[i].c_str());

		
    // to remove if not necessary
	/*	
		 if((C_SmileOrNotId = ARM_ConvSmileNoSmile ((new CCString(C_SmileOrNot.c_str())), C_result)) == ARM_DEFAULT_ERR)
	    {
		    ARM_ARG_ERR();
		    return (LPXLOPER)&XL_result;
	    }
    */

        	CCString curClass = LOCAL_VOL_CUBE_CLASS;
	    CCString stringId = GetLastCurCellEnvValue ();
		CCString prevClass;
		long objId=-1;
		if (!(!stringId)){
			prevClass = LocalGetStringObjectClass (stringId);
			if ( curClass == prevClass ){
				objId = LocalGetNumObjectId (stringId);
			}else{
				FreeCurCellContent();
			}
		}

		//long prevId = ExcelCaller::get().getObjectId();  
		long retCode =  ARMLOCAL_GetVolCubeFromCalypso (C_index,
												     C_currency,
												     C_cvName,
												     C_date,
												     C_vtype,
                                                     C_type,
                                                     C_suffix,
												     C_tenors,
												     C_SmileOrNot,
													 indexId_,
												     C_result,
													 objId);

	    	    
	   
			
		
     if ( retCode == ARM_OK )
	    {
		   	std::string objectLabel =ExcelCaller::get().setObject(C_result.getLong(), LOCAL_VOL_CUBE_CLASS) ;
	    	// ExcelCaller::get().setObject(objectLabel) ;
	    	ExcelTools::convert(objectLabel,&XL_result); 
		    return &XL_result;		
	    }
	    else
	    {
		    ARM_ERR();
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
	return &XL_result;
}

__declspec(dllexport) LPXLOPER WINAPI Local_GetVolCubeFromSummit(LPXLOPER XL_index,
																 LPXLOPER XL_currency,
																 LPXLOPER XL_cvName,
																 LPXLOPER XL_date,
																 LPXLOPER XL_vtype,
																 LPXLOPER XL_tenors,
																 LPXLOPER XL_SmileOrNot,
																 LPXLOPER XL_indexid,
																 LPXLOPER XL_SmileType)
{
	ADD_LOG("Local_GetVolCubeFromSummit");
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
	    double C_date;
	    CCString C_vtype;

	    VECTOR<CCString> C_tenors;
	    VECTOR<CCString> C_tenors_default;

		CCString C_SmileOrNot;
	    long C_SmileOrNotId;
	    CCString C_StrSmileType, defST("LINEAR");
		int C_SmileType;

	    // error
	    static int error;
	    static char* reason = "";

	    XL_readStrCell(XL_index,C_index," ARM_ERR: index: string expected",C_result);
	    XL_readStrCell(XL_currency,C_currency," ARM_ERR: currency: string expected",C_result);
	    XL_readStrCell(XL_cvName,C_cvName," ARM_ERR: curve name: string expected",C_result);
	    XL_readNumCell(XL_date,C_date," ARM_ERR: as of date: date expected",C_result);
	    XL_readStrCell(XL_vtype,C_vtype," ARM_ERR: volatility type: string expected",C_result);
        XL_readStrVectorWD(XL_tenors,C_tenors,C_tenors_default," ARM_ERR: tenors (3M, 1Y, ...: array of numeric expected",XL_TYPE_STRING,C_result);
	    XL_readStrCellWD(XL_SmileOrNot,C_SmileOrNot,"NOSMILE"," ARM_ERR: Smile, NoSmile: string expected",C_result);
		XL_readStrCellWD(XL_SmileType, C_StrSmileType, defST, "ARM_ERR : SmileType : string expected",C_result);

		std::string indexId ; ExcelTools::convert(XL_indexid,"",indexId); 
		long indexId_=LocalPersistent::get().getObjectId(indexId) ; 
	    C_index.toUpper ();
	    C_currency.toUpper ();
	    C_cvName.toUpper ();

	    for (int i = 0; i < C_tenors.size(); i++)
		    C_tenors[i].toUpper();

	    if((C_SmileOrNotId = ARM_ConvSmileNoSmile (&C_SmileOrNot, C_result)) == ARM_DEFAULT_ERR)
	    {
		    ARM_ARG_ERR();
		    return (LPXLOPER)&XL_result;
	    }

	    long retCode;
	    long objId;
	    CCString prevClass;
	    
	    CCString curClass = LOCAL_VOL_CUBE_CLASS;
	    CCString stringId = GetLastCurCellEnvValue ();
	    
		C_StrSmileType.toUpper();
		C_SmileType = C_StrSmileType == "SPLINE" ? 1 : 0;

	    if (!stringId)
	    {
		    retCode = ARMLOCAL_GetVolCubeFromSummit (C_index,
												     C_currency,
												     C_cvName,
												     C_date,
												     C_vtype,
												     C_tenors,
												     C_SmileOrNot,
													 indexId_,
													 C_SmileType,
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
			    retCode = ARMLOCAL_GetVolCubeFromSummit(C_index,
													    C_currency,
													    C_cvName,
													    C_date,
													    C_vtype,
													    C_tenors,
													    C_SmileOrNot,
														indexId_,
														C_SmileType,
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
			    FreeCurCellContent();

			    retCode = ARMLOCAL_GetVolCubeFromSummit(C_index,
													    C_currency,
													    C_cvName,
													    C_date,
													    C_vtype,
													    C_tenors,
													    C_SmileOrNot,
														indexId_,
														C_SmileType,
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_GetVolCubeFromSummit" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}



_declspec(dllexport) LPXLOPER WINAPI Local_PXL_GetVolCubeFromSummit(LPXLOPER XL_index,
																	LPXLOPER XL_currency,
																	LPXLOPER XL_cvName,
																	LPXLOPER XL_date,
																	LPXLOPER XL_vtype,
																	LPXLOPER XL_tenors,
																	LPXLOPER XL_SmileOrNot,
																	LPXLOPER XL_IndexId,
																	LPXLOPER XL_SmileType)
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
	    CCString C_index;
	    CCString C_currency;
	    CCString C_cvName;
	    double C_date;
	    CCString C_vtype;

		VECTOR<CCString> C_tenors;
	    VECTOR<CCString> C_tenors_default;

		CCString C_SmileOrNot;
	    long C_SmileOrNotId;
	    
	    CCString C_StrSmileType;
		int C_SmileType;

		// error
	    static int error;
	    static char* reason = "";

	    XL_readStrCell(XL_index,C_index," ARM_ERR: index: string expected",C_result);
	    XL_readStrCell(XL_currency,C_currency," ARM_ERR: currency: string expected",C_result);
	    XL_readStrCell(XL_cvName,C_cvName," ARM_ERR: curve name: string expected",C_result);
	    XL_readNumCell(XL_date,C_date," ARM_ERR: as of date: date expected",C_result);
	    XL_readStrCell(XL_vtype,C_vtype," ARM_ERR: volatility type: string expected",C_result);
        XL_readStrVectorWD(XL_tenors,C_tenors,C_tenors_default," ARM_ERR: tenors (3M, 1Y, ...: array of numeric expected",XL_TYPE_STRING,C_result);
	    XL_readStrCellWD(XL_SmileOrNot,C_SmileOrNot,"NOSMILE"," ARM_ERR: Smile, NoSmile: string expected",C_result);
		XL_readStrCellWD(XL_SmileType, C_StrSmileType, "LINEAR", "ARM_ERR : SmileType : string expected",C_result);

		std::string indexId ; ExcelTools::convert(XL_IndexId,"",indexId) ;
		long indexId_ = LocalPersistent::get().getObjectId(indexId); 

	    C_index.toUpper ();
	    C_currency.toUpper ();
	    C_cvName.toUpper ();

	    for (int i = 0; i < C_tenors.size(); i++)
		    C_tenors[i].toUpper();

	    if((C_SmileOrNotId = ARM_ConvSmileNoSmile (&C_SmileOrNot, C_result)) == ARM_DEFAULT_ERR)
	    {
		    ARM_ARG_ERR();
		    return (LPXLOPER)&XL_result;
	    }

	    long retCode;
	    long objId;
	    
	    CCString curClass = LOCAL_VOL_CUBE_CLASS;
	    CCString stringId;
	    
		C_StrSmileType.toUpper();
		C_SmileType = C_StrSmileType == "SPLINE" ? 1 : 0;

		retCode = ARMLOCAL_GetVolCubeFromSummit(C_index,
											    C_currency,
											    C_cvName,
											    C_date,
											    C_vtype,
											    C_tenors,
											    C_SmileOrNot,
											    indexId_,
												C_SmileType,
												C_result);

	    if ( retCode == ARM_OK )
	    {
	       objId = C_result.getLong ();

	       stringId = LocalMakeObjectId (objId, curClass);
	    }

	    if ( retCode == ARM_OK )
	    {
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_PXL_GetVolCubeFromSummit" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}



__declspec(dllexport) LPXLOPER WINAPI Local_ARM_BumpVolatility (LPXLOPER XL_vol,
																LPXLOPER XL_value,
																LPXLOPER XL_nthLine,
																LPXLOPER XL_nthCol,
																LPXLOPER XL_isCumul,
																LPXLOPER XL_isAbsolute)
{
	ADD_LOG("Local_ARM_BumpVolatility ");
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	    ARM_NOCALCIFWIZ();


	    // C variable
	    CCString C_vol;

	    double C_value;

	    double C_nthLine;
	    double C_nthLine_default = 0.;

	    double C_nthCol;
	    double C_nthCol_default = 0.;

	    CCString C_isCumul;
	    long cumulId;

	    CCString C_isAbsolute;
	    long absoluteId;

		// error
	    static int error;
	    static char* reason = "";

	    XL_readStrCell(XL_vol,C_vol," ARM_ERR: volatility: object expected",C_result);
	    XL_readNumCell(XL_value,C_value," ARM_ERR: value to bump: numeric expected",C_result);
	    XL_readNumCellWD(XL_nthLine,C_nthLine,C_nthLine_default," ARM_ERR: nth Line: numeric expected",C_result);
	    XL_readNumCellWD(XL_nthCol,C_nthCol,C_nthCol_default," ARM_ERR: nth Column: numeric expected",C_result);
	    XL_readStrCellWD(XL_isCumul,C_isCumul,"NO"," ARM_ERR: is cumulative?: string expected",C_result);
	    XL_readStrCellWD(XL_isAbsolute,C_isAbsolute,"YES"," ARM_ERR: is absolute?: string expected",C_result);

		if(C_isCumul[0]=='Y' || C_isCumul[0]=='1')
			cumulId = 1;
		else if(C_isCumul[0]=='N' || C_isCumul[0]=='0')
			cumulId = 0;
		else if( C_isCumul[0] == 'C')
			cumulId = -1;
		else
			ARM_THROW(ERR_INVALID_ARGUMENT, "isCumulative should be Y[ES], N[O], or C[ustom]"); 


	    if((absoluteId = ARM_ConvYesOrNo (C_isAbsolute, C_result)) == ARM_DEFAULT_ERR)
	    {
		    ARM_ARG_ERR();
		    return (LPXLOPER)&XL_result;
	    }

		long retCode;
	    long objId;
	    CCString prevClass;

	    CCString curClass = LocalGetStringObjectClass (C_vol);
	    CCString stringId = GetLastCurCellEnvValue ();

	    if(!stringId)
	    {
		    retCode = ARMLOCAL_ARM_BumpVolatility (LocalGetNumObjectId(C_vol),
											       C_value,
											       (long)C_nthLine,
											       (long)C_nthCol,
											       cumulId,
												   absoluteId,
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
			    retCode = ARMLOCAL_ARM_BumpVolatility (LocalGetNumObjectId(C_vol),
												       C_value,
												       (long)C_nthLine,
												       (long)C_nthCol,
												       cumulId,
													   absoluteId,
												       C_result,
												       objId);

			    if (retCode == ARM_OK)
			    {
				    LocalSetCurCellEnvValue (curClass, objId); 

				    stringId = LocalMakeObjectId (objId, curClass);
			    }
		    }
		    else
		    {
			    FreeCurCellContent ();

			    retCode = ARMLOCAL_ARM_BumpVolatility (LocalGetNumObjectId(C_vol),
												       C_value,
												       (long)C_nthLine,
												       (long)C_nthCol,
												       cumulId,
													   absoluteId,
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_ARM_BumpVolatility" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}



__declspec(dllexport) LPXLOPER WINAPI Local_PXL_ARM_BumpVolatility(LPXLOPER XL_vol,
																   LPXLOPER XL_value,
																   LPXLOPER XL_nthLine,
																   LPXLOPER XL_nthCol,
																   LPXLOPER XL_isCumul,
																   LPXLOPER XL_isAbsolute)
{
	ADD_LOG("Local_PXL_ARM_BumpVolatility");
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	    ARM_NOCALCIFWIZ();


	    // C variable
	    CCString C_vol;

	    double C_value;

	    double C_nthLine;
	    double C_nthLine_default = 0.;

	    double C_nthCol;
	    double C_nthCol_default = 0.;

	    CCString C_isCumul;
	    long cumulId;

	    CCString C_isAbsolute;
	    long absoluteId;

		// error
	    static int error;
	    static char* reason = "";

	    XL_readStrCell(XL_vol,C_vol," ARM_ERR: volatility: object expected",C_result);
	    XL_readNumCell(XL_value,C_value," ARM_ERR: value to bump: numeric expected",C_result);
	    XL_readNumCellWD(XL_nthLine,C_nthLine,C_nthLine_default," ARM_ERR: nth Line: numeric expected",C_result);
	    XL_readNumCellWD(XL_nthCol,C_nthCol,C_nthCol_default," ARM_ERR: nth Column: numeric expected",C_result);
	    XL_readStrCellWD(XL_isCumul,C_isCumul,"NO"," ARM_ERR: is cumulative?: string expected",C_result);
	    XL_readStrCellWD(XL_isAbsolute,C_isAbsolute,"YES"," ARM_ERR: is absolute?: string expected",C_result);

	    if((cumulId = ARM_ConvYesOrNo (C_isCumul, C_result)) == ARM_DEFAULT_ERR)
	    {
		    ARM_ARG_ERR();
		    return (LPXLOPER)&XL_result;
	    }

	    if((absoluteId = ARM_ConvYesOrNo (C_isAbsolute, C_result)) == ARM_DEFAULT_ERR)
	    {
		    ARM_ARG_ERR();
		    return (LPXLOPER)&XL_result;
	    }

	    long retCode;
	    long objId;

	    CCString curClass = LocalGetStringObjectClass (C_vol);
	    CCString stringId;
	    
	    retCode = ARMLOCAL_ARM_BumpVolatility (LocalGetNumObjectId(C_vol),
										       C_value,
										       (long)C_nthLine,
										       (long)C_nthCol,
										       cumulId,
											   absoluteId,
										       C_result);

	    if ( retCode == ARM_OK )
	    {
		    objId = C_result.getLong ();

		    stringId = LocalMakeObjectId (objId, curClass);
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_PXL_ARM_BumpVolatility" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}

__declspec(dllexport) LPXLOPER WINAPI Local_ARM_BumpSmile(LPXLOPER XL_vol,
														  LPXLOPER XL_value,
														  LPXLOPER XL_tenor,
														  LPXLOPER XL_nthLine,
														  LPXLOPER XL_nthCol,
														  LPXLOPER XL_isCumul,
														  LPXLOPER XL_isAbsolute)
{
	ADD_LOG("Local_ARM_BumpSmile");
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	    ARM_NOCALCIFWIZ();


	    // C variable
	    CCString C_vol;

	    double C_value;
		double C_tenor;
		double C_tenor_default = 0.;

	    double C_nthLine;
	    double C_nthLine_default = 0.;

	    double C_nthCol;
	    double C_nthCol_default = 0.;

	    CCString C_isCumul;
	    long cumulId;

	    CCString C_isAbsolute;
	    long absoluteId;

		// error
	    static int error;
	    static char* reason = "";

	    XL_readStrCell(XL_vol,C_vol," ARM_ERR: volatility: object expected",C_result);
	    XL_readNumCell(XL_value,C_value," ARM_ERR: value to bump: numeric expected",C_result);
		XL_readNumCellWD(XL_tenor,C_tenor,C_tenor_default," ARM_ERR: tenor: numeric expected",C_result);
	    XL_readNumCellWD(XL_nthLine,C_nthLine,C_nthLine_default," ARM_ERR: nth Line: numeric expected",C_result);
	    XL_readNumCellWD(XL_nthCol,C_nthCol,C_nthCol_default," ARM_ERR: nth Column: numeric expected",C_result);
	    XL_readStrCellWD(XL_isCumul,C_isCumul,"NO"," ARM_ERR: is cumulative?: string expected",C_result);
	    XL_readStrCellWD(XL_isAbsolute,C_isAbsolute,"YES"," ARM_ERR: is absolute?: string expected",C_result);

	    if((cumulId = ARM_ConvYesOrNo (C_isCumul, C_result)) == ARM_DEFAULT_ERR)
	    {
		    ARM_ARG_ERR();
		    return (LPXLOPER)&XL_result;
	    }

	    if((absoluteId = ARM_ConvYesOrNo (C_isAbsolute, C_result)) == ARM_DEFAULT_ERR)
	    {
		    ARM_ARG_ERR();
		    return (LPXLOPER)&XL_result;
	    }

		long retCode;
	    long objId;
	    CCString prevClass;

	    CCString curClass = LocalGetStringObjectClass (C_vol);
	    CCString stringId = GetLastCurCellEnvValue ();

	    if(!stringId)
	    {
		    retCode = ARMLOCAL_ARM_BumpSmile(LocalGetNumObjectId(C_vol),
											 C_value,
											 C_tenor,
											 (long)C_nthLine,
											 (long)C_nthCol,
											 cumulId,
											 absoluteId,
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
			    retCode = ARMLOCAL_ARM_BumpSmile(LocalGetNumObjectId(C_vol),
												 C_value,
												 C_tenor,
												 (long)C_nthLine,
												 (long)C_nthCol,
												 cumulId,
												 absoluteId,
												 C_result,
												 objId);

			    if (retCode == ARM_OK)
			    {
				    LocalSetCurCellEnvValue (curClass, objId); 

				    stringId = LocalMakeObjectId (objId, curClass);
			    }
		    }
		    else
		    {
			    FreeCurCellContent ();

			    retCode = ARMLOCAL_ARM_BumpSmile(LocalGetNumObjectId(C_vol),
												 C_value,
												 C_tenor,
												 (long)C_nthLine,
												 (long)C_nthCol,
												 cumulId,
												 absoluteId,
												 C_result);

			    if(retCode == ARM_OK)
			    {
				    objId = C_result.getLong();
			    
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_ARM_BumpSmile" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}



__declspec(dllexport) LPXLOPER WINAPI Local_PXL_ARM_BumpSmile(LPXLOPER XL_vol,
															  LPXLOPER XL_value,
															  LPXLOPER XL_tenor,
															  LPXLOPER XL_nthLine,
															  LPXLOPER XL_nthCol,
															  LPXLOPER XL_isCumul,
															  LPXLOPER XL_isAbsolute)
{
	ADD_LOG("Local_PXL_ARM_BumpSmile");
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	    ARM_NOCALCIFWIZ();


	    // C variable
	    CCString C_vol;

	    double C_value;
		double C_tenor;
		double C_tenor_default = 0.;

	    double C_nthLine;
	    double C_nthLine_default = 0.;

	    double C_nthCol;
	    double C_nthCol_default = 0.;

	    CCString C_isCumul;
	    long cumulId;

	    CCString C_isAbsolute;
	    long absoluteId;

		// error
	    static int error;
	    static char* reason = "";

	    XL_readStrCell(XL_vol,C_vol," ARM_ERR: volatility: object expected",C_result);
	    XL_readNumCell(XL_value,C_value," ARM_ERR: value to bump: numeric expected",C_result);
		XL_readNumCellWD(XL_tenor,C_tenor,C_tenor_default," ARM_ERR: tenor: numeric expected",C_result);
	    XL_readNumCellWD(XL_nthLine,C_nthLine,C_nthLine_default," ARM_ERR: nth Line: numeric expected",C_result);
	    XL_readNumCellWD(XL_nthCol,C_nthCol,C_nthCol_default," ARM_ERR: nth Column: numeric expected",C_result);
	    XL_readStrCellWD(XL_isCumul,C_isCumul,"NO"," ARM_ERR: is cumulative?: string expected",C_result);
	    XL_readStrCellWD(XL_isAbsolute,C_isAbsolute,"YES"," ARM_ERR: is absolute?: string expected",C_result);

	    if((cumulId = ARM_ConvYesOrNo (C_isCumul, C_result)) == ARM_DEFAULT_ERR)
	    {
		    ARM_ARG_ERR();
		    return (LPXLOPER)&XL_result;
	    }

	    if((absoluteId = ARM_ConvYesOrNo (C_isAbsolute, C_result)) == ARM_DEFAULT_ERR)
	    {
		    ARM_ARG_ERR();
		    return (LPXLOPER)&XL_result;
	    }

	    long retCode;
	    long objId;

	    CCString curClass = LocalGetStringObjectClass (C_vol);
	    CCString stringId;
	    
	    retCode = ARMLOCAL_ARM_BumpSmile(LocalGetNumObjectId(C_vol),
										 C_value,
										 C_tenor,
										 (long)C_nthLine,
										 (long)C_nthCol,
										 cumulId,
										 absoluteId,
										 C_result);

	    if ( retCode == ARM_OK )
	    {
		    objId = C_result.getLong ();

		    stringId = LocalMakeObjectId (objId, curClass);
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_PXL_ARM_BumpSmile" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


__declspec(dllexport) LPXLOPER WINAPI Local_ARM_BumpHyperCubeSmile(LPXLOPER XL_vol,
																  LPXLOPER XL_value,
																  LPXLOPER XL_cubeTenor,
																  LPXLOPER XL_smileTenor,
																  LPXLOPER XL_nthLine,
																  LPXLOPER XL_nthCol,
																  LPXLOPER XL_isCumul,
																  LPXLOPER XL_isAbsolute)
{
	ADD_LOG("Local_ARM_BumpHyperCubeSmile");
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	    ARM_NOCALCIFWIZ();


	    // C variable
	    CCString C_vol;

	    double C_value;

		double C_cubeTenor;
		double C_cubeTenor_default = 0.;

		double C_smileTenor;
		double C_smileTenor_default = 0.;

	    double C_nthLine;
	    double C_nthLine_default = 0.;

	    double C_nthCol;
	    double C_nthCol_default = 0.;

	    CCString C_isCumul;
	    long cumulId;

	    CCString C_isAbsolute;
	    long absoluteId;

		// error
	    static int error;
	    static char* reason = "";

	    XL_readStrCell(XL_vol,C_vol," ARM_ERR: volatility: object expected",C_result);
	    XL_readNumCell(XL_value,C_value," ARM_ERR: value to bump: numeric expected",C_result);
		XL_readNumCellWD(XL_cubeTenor,C_cubeTenor,C_cubeTenor_default," ARM_ERR: cube tenor: numeric expected",C_result);
		XL_readNumCellWD(XL_smileTenor,C_smileTenor,C_smileTenor_default," ARM_ERR: smile tenor: numeric expected",C_result);
	    XL_readNumCellWD(XL_nthLine,C_nthLine,C_nthLine_default," ARM_ERR: nth Line: numeric expected",C_result);
	    XL_readNumCellWD(XL_nthCol,C_nthCol,C_nthCol_default," ARM_ERR: nth Column: numeric expected",C_result);
	    XL_readStrCellWD(XL_isCumul,C_isCumul,"NO"," ARM_ERR: is cumulative?: string expected",C_result);
	    XL_readStrCellWD(XL_isAbsolute,C_isAbsolute,"YES"," ARM_ERR: is absolute?: string expected",C_result);

	    if((cumulId = ARM_ConvYesOrNo (C_isCumul, C_result)) == ARM_DEFAULT_ERR)
	    {
		    ARM_ARG_ERR();
		    return (LPXLOPER)&XL_result;
	    }

	    if((absoluteId = ARM_ConvYesOrNo (C_isAbsolute, C_result)) == ARM_DEFAULT_ERR)
	    {
		    ARM_ARG_ERR();
		    return (LPXLOPER)&XL_result;
	    }

		long retCode;
	    long objId;
	    CCString prevClass;

	    CCString curClass = LocalGetStringObjectClass (C_vol);
	    CCString stringId = GetLastCurCellEnvValue ();

	    if(!stringId)
	    {
		    retCode = ARMLOCAL_ARM_BumpHyperCubeSmile(LocalGetNumObjectId(C_vol),
													 C_value,
													 C_cubeTenor,
													 C_smileTenor,
													 (long)C_nthLine,
													 (long)C_nthCol,
													 cumulId,
													 absoluteId,
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
				retCode = ARMLOCAL_ARM_BumpHyperCubeSmile(LocalGetNumObjectId(C_vol),
														 C_value,
														 C_cubeTenor,
														 C_smileTenor,
														 (long)C_nthLine,
														 (long)C_nthCol,
														 cumulId,
														 absoluteId,
														 C_result,
														 objId);

			    if (retCode == ARM_OK)
			    {
				    LocalSetCurCellEnvValue (curClass, objId); 
				    stringId = LocalMakeObjectId (objId, curClass);
			    }
		    }
		    else
		    {
			    FreeCurCellContent ();
				retCode = ARMLOCAL_ARM_BumpHyperCubeSmile(LocalGetNumObjectId(C_vol),
														 C_value,
														 C_cubeTenor,
														 C_smileTenor,
														 (long)C_nthLine,
														 (long)C_nthCol,
														 cumulId,
														 absoluteId,
														 C_result);

			    if(retCode == ARM_OK)
			    {
				    objId = C_result.getLong();			    
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_ARM_BumpHyperCubeSmile" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}



__declspec(dllexport) LPXLOPER WINAPI Local_PXL_ARM_BumpHyperCubeSmile(LPXLOPER XL_vol,
																	  LPXLOPER XL_value,
																	  LPXLOPER XL_cubeTenor,
																	  LPXLOPER XL_smileTenor,
																	  LPXLOPER XL_nthLine,
																	  LPXLOPER XL_nthCol,
																	  LPXLOPER XL_isCumul,
																	  LPXLOPER XL_isAbsolute)
{
	ADD_LOG("Local_PXL_ARM_BumpHyperCubeSmile");
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	    ARM_NOCALCIFWIZ();


	    // C variable
	    CCString C_vol;

	    double C_value;

		double C_cubeTenor;
		double C_cubeTenor_default = 0.;

		double C_smileTenor;
		double C_smileTenor_default = 0.;
	    
		double C_nthLine;
	    double C_nthLine_default = 0.;

	    double C_nthCol;
	    double C_nthCol_default = 0.;

	    CCString C_isCumul;
	    long cumulId;

	    CCString C_isAbsolute;
	    long absoluteId;

		// error
	    static int error;
	    static char* reason = "";

	    XL_readStrCell(XL_vol,C_vol," ARM_ERR: volatility: object expected",C_result);
	    XL_readNumCell(XL_value,C_value," ARM_ERR: value to bump: numeric expected",C_result);
		XL_readNumCellWD(XL_cubeTenor,C_cubeTenor,C_cubeTenor_default," ARM_ERR: cube tenor: numeric expected",C_result);
		XL_readNumCellWD(XL_smileTenor,C_smileTenor,C_smileTenor_default," ARM_ERR: smile tenor: numeric expected",C_result);
	    XL_readNumCellWD(XL_nthLine,C_nthLine,C_nthLine_default," ARM_ERR: nth Line: numeric expected",C_result);
	    XL_readNumCellWD(XL_nthCol,C_nthCol,C_nthCol_default," ARM_ERR: nth Column: numeric expected",C_result);
	    XL_readStrCellWD(XL_isCumul,C_isCumul,"NO"," ARM_ERR: is cumulative?: string expected",C_result);
	    XL_readStrCellWD(XL_isAbsolute,C_isAbsolute,"YES"," ARM_ERR: is absolute?: string expected",C_result);

	    if((cumulId = ARM_ConvYesOrNo (C_isCumul, C_result)) == ARM_DEFAULT_ERR)
	    {
		    ARM_ARG_ERR();
		    return (LPXLOPER)&XL_result;
	    }

	    if((absoluteId = ARM_ConvYesOrNo (C_isAbsolute, C_result)) == ARM_DEFAULT_ERR)
	    {
		    ARM_ARG_ERR();
		    return (LPXLOPER)&XL_result;
	    }

	    long retCode;
	    long objId;

	    CCString curClass = LocalGetStringObjectClass (C_vol);
	    CCString stringId;
	    
		retCode = ARMLOCAL_ARM_BumpHyperCubeSmile(LocalGetNumObjectId(C_vol),
												 C_value,
												 C_cubeTenor,
												 C_smileTenor,
												 (long)C_nthLine,
												 (long)C_nthCol,
												 cumulId,
												 absoluteId,
												 C_result);
	    if ( retCode == ARM_OK )
	    {
		    objId = C_result.getLong ();

		    stringId = LocalMakeObjectId (objId, curClass);
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_PXL_ARM_BumpHyperCubeSmile" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}



__declspec(dllexport) LPXLOPER WINAPI Local_ARM_FXBumpRRorSTR ( LPXLOPER XL_vol,
																LPXLOPER XL_value,
																LPXLOPER XL_nthLine,
																LPXLOPER XL_nthCol,
																LPXLOPER XL_spotFX,
																LPXLOPER XL_isCumul,
																LPXLOPER XL_isAbsolute,
																LPXLOPER XL_isRR)
{
	ADD_LOG("Local_ARM_FXBumpRRorSTR ");
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	    ARM_NOCALCIFWIZ();


	    // C variable
	    CCString C_vol;

	    double C_value;

	    double C_nthLine;
	    double C_nthLine_default = 0.;

	    double C_nthCol;
	    double C_nthCol_default = 0.;

		double C_spotFX;
	    double C_spotFX_default = 0.;

	    CCString C_isCumul;
	    long cumulId;

	    CCString C_isAbsolute;
	    long absoluteId;

		CCString C_isRR;
	    long isRRId;

		// error
	    static int error;
	    static char* reason = "";

	    XL_readStrCell(XL_vol,C_vol," ARM_ERR: volatility: object expected",C_result);
	    XL_readNumCell(XL_value,C_value," ARM_ERR: value to bump: numeric expected",C_result);
	    XL_readNumCellWD(XL_nthLine,C_nthLine,C_nthLine_default," ARM_ERR: nth Line: numeric expected",C_result);
	    XL_readNumCellWD(XL_nthCol,C_nthCol,C_nthCol_default," ARM_ERR: nth Column: numeric expected",C_result);
		XL_readNumCellWD(XL_spotFX,C_spotFX,C_spotFX_default," ARM_ERR: spot FX: numeric expected",C_result);
	    XL_readStrCellWD(XL_isCumul,C_isCumul,"NO"," ARM_ERR: is cumulative?: string expected",C_result);
	    XL_readStrCellWD(XL_isAbsolute,C_isAbsolute,"YES"," ARM_ERR: is absolute?: string expected",C_result);
		XL_readStrCellWD(XL_isRR,C_isRR,"YES"," ARM_ERR: is RR?: string expected",C_result);

	    if((cumulId = ARM_ConvYesOrNo (C_isCumul, C_result)) == ARM_DEFAULT_ERR)
	    {
		    ARM_ARG_ERR();
		    return (LPXLOPER)&XL_result;
	    }

	    if((absoluteId = ARM_ConvYesOrNo (C_isAbsolute, C_result)) == ARM_DEFAULT_ERR)
	    {
		    ARM_ARG_ERR();
		    return (LPXLOPER)&XL_result;
	    }

		if((isRRId = ARM_ConvYesOrNo (C_isRR, C_result)) == ARM_DEFAULT_ERR)
	    {
		    ARM_ARG_ERR();
		    return (LPXLOPER)&XL_result;
	    }

		long retCode;
	    long objId;
	    CCString prevClass;

	    CCString curClass = LocalGetStringObjectClass (C_vol);
	    CCString stringId = GetLastCurCellEnvValue ();

	    if(!stringId)
	    {
		    retCode = ARMLOCAL_ARM_FXBumpRRorSTR ( LocalGetNumObjectId(C_vol),
											       C_value,
											       (long)C_nthLine,
											       (long)C_nthCol,
												   C_spotFX,
											       cumulId,
												   absoluteId,
												   isRRId,
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
			    retCode = ARMLOCAL_ARM_FXBumpRRorSTR(	LocalGetNumObjectId(C_vol),
														C_value,
														(long)C_nthLine,
														(long)C_nthCol,
														C_spotFX,
														cumulId,
														absoluteId,
														isRRId,
														C_result,
														objId);

			    if (retCode == ARM_OK)
			    {
				    LocalSetCurCellEnvValue (curClass, objId); 

				    stringId = LocalMakeObjectId (objId, curClass);
			    }
		    }
		    else
		    {
			    FreeCurCellContent ();

			    retCode = ARMLOCAL_ARM_FXBumpRRorSTR (	LocalGetNumObjectId(C_vol),
														C_value,
														(long)C_nthLine,
														(long)C_nthCol,
														C_spotFX,
														cumulId,
														absoluteId,
														isRRId,
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_ARM_BumpVolatility" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}



_declspec(dllexport) LPXLOPER WINAPI Local_PXL_ARM_FXBumpRRorSTR ( LPXLOPER XL_vol,
																LPXLOPER XL_value,
																LPXLOPER XL_nthLine,
																LPXLOPER XL_nthCol,
																LPXLOPER XL_spotFX,
																LPXLOPER XL_isCumul,
																LPXLOPER XL_isAbsolute,
																LPXLOPER XL_isRR)
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
	    CCString C_vol;

	    double C_value;

	    double C_nthLine;
	    double C_nthLine_default = 0.;

	    double C_nthCol;
	    double C_nthCol_default = 0.;

		double C_spotFX;
	    double C_spotFX_default = 0.;

	    CCString C_isCumul;
	    long cumulId;

	    CCString C_isAbsolute;
	    long absoluteId;

		CCString C_isRR;
	    long isRRId;

		// error
	    static int error;
	    static char* reason = "";

	    XL_readStrCell(XL_vol,C_vol," ARM_ERR: volatility: object expected",C_result);
	    XL_readNumCell(XL_value,C_value," ARM_ERR: value to bump: numeric expected",C_result);
	    XL_readNumCellWD(XL_nthLine,C_nthLine,C_nthLine_default," ARM_ERR: nth Line: numeric expected",C_result);
	    XL_readNumCellWD(XL_nthCol,C_nthCol,C_nthCol_default," ARM_ERR: nth Column: numeric expected",C_result);
		XL_readNumCellWD(XL_spotFX,C_spotFX,C_spotFX_default," ARM_ERR: spot FX: numeric expected",C_result);
	    XL_readStrCellWD(XL_isCumul,C_isCumul,"NO"," ARM_ERR: is cumulative?: string expected",C_result);
	    XL_readStrCellWD(XL_isAbsolute,C_isAbsolute,"YES"," ARM_ERR: is absolute?: string expected",C_result);
		XL_readStrCellWD(XL_isRR,C_isRR,"YES"," ARM_ERR: is RR?: string expected",C_result);

	    if((cumulId = ARM_ConvYesOrNo (C_isCumul, C_result)) == ARM_DEFAULT_ERR)
	    {
		    ARM_ARG_ERR();
		    return (LPXLOPER)&XL_result;
	    }

	    if((absoluteId = ARM_ConvYesOrNo (C_isAbsolute, C_result)) == ARM_DEFAULT_ERR)
	    {
		    ARM_ARG_ERR();
		    return (LPXLOPER)&XL_result;
	    }

		if((isRRId = ARM_ConvYesOrNo (C_isRR, C_result)) == ARM_DEFAULT_ERR)
	    {
		    ARM_ARG_ERR();
		    return (LPXLOPER)&XL_result;
	    }

		long retCode;
	    long objId;
	    CCString prevClass;

	    CCString curClass = LocalGetStringObjectClass (C_vol);
	    CCString stringId;

		retCode = ARMLOCAL_ARM_FXBumpRRorSTR ( LocalGetNumObjectId(C_vol),
											   C_value,
											   (long)C_nthLine,
											   (long)C_nthCol,
											   C_spotFX,
											   cumulId,
											   absoluteId,
											   isRRId,
											   C_result);

		if(retCode == ARM_OK)
		{
			objId = C_result.getLong ();

			LocalSetCurCellEnvValue (curClass, objId); 
			stringId = LocalMakeObjectId (objId, curClass);

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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_ARM_BumpVolatility" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


__declspec(dllexport) LPXLOPER WINAPI Local_ARM_GetFXVolFromSummit(LPXLOPER XL_ccy1,
																   LPXLOPER XL_ccy2,
																   LPXLOPER XL_date,
																   LPXLOPER XL_cvName,
																   LPXLOPER XL_impOrHist,
																   LPXLOPER XL_volType)
{
	ADD_LOG("Local_ARM_GetFXVolFromSummit");
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

    /// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	    ARM_NOCALCIFWIZ();


	    // C variable
	    CCString C_ccy1;
	    CCString C_ccy2;
	    CCString C_cvName;
	    double C_date;
	    CCString C_impOrHist;
	    CCString C_volType;
	    
	    // error
	    static int error;
	    static char* reason = "";

	    XL_readStrCell(XL_ccy1,C_ccy1," ARM_ERR: ccy1: string expected",C_result);
	    XL_readStrCell(XL_ccy2,C_ccy2," ARM_ERR: ccy2: string expected",C_result);
	    XL_readNumCell(XL_date,C_date," ARM_ERR: as of date: date expected",C_result);
	    XL_readStrCellWD(XL_cvName,C_cvName,"MO"," ARM_ERR: curve name: string expected",C_result);
	    XL_readStrCellWD(XL_impOrHist,C_impOrHist,"FXVOL"," ARM_ERR: implied or historical vol: string expected",C_result);
	    XL_readStrCellWD(XL_volType,C_volType,"ATM"," ARM_ERR: Volatility type: string expected",C_result);

	    C_ccy1.toUpper ();
	    C_ccy2.toUpper ();
	    C_cvName.toUpper ();
	    C_impOrHist.toUpper ();
	    C_volType.toUpper ();

	    long retCode;
	    long objId;
	    CCString prevClass;

	    CCString curClass;
	    CCString stringId = GetLastCurCellEnvValue ();
	    
	    if ((strcmp((const char*)C_volType,"CSMILE")==0) || (strcmp((const char*)C_volType,"CFXSPI")==0))
		    curClass = LOCAL_VOL_CUBE_CLASS;
	    else
		    curClass = LOCAL_VOL_CURVE_LIN_CLASS;

	    if (!stringId)
	    {
		    retCode = ARMLOCAL_GetFXVolFromSummit (C_ccy1,
											       C_ccy2,
											       C_date,
											       C_cvName,
											       C_impOrHist,
											       C_volType,
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

		    if ( (prevClass == LOCAL_VOL_CUBE_CLASS) || (prevClass == LOCAL_VOL_CURVE_LIN_CLASS) )
		    {
			    retCode = ARMLOCAL_GetFXVolFromSummit (C_ccy1,
												       C_ccy2,
												       C_date,
												       C_cvName,
												       C_impOrHist,
												       C_volType,
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
			    retCode = ARMLOCAL_GetFXVolFromSummit (C_ccy1,
												       C_ccy2,
												       C_date,
												       C_cvName,
												       C_impOrHist,
												       C_volType,
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_ARM_GetFXVolFromSummit" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}




__declspec(dllexport) LPXLOPER WINAPI Local_PXL_ARM_GetFXVolFromSummit(LPXLOPER XL_ccy1,
																	   LPXLOPER XL_ccy2,
																	   LPXLOPER XL_date,
																	   LPXLOPER XL_cvName,
																	   LPXLOPER XL_impOrHist,
																	   LPXLOPER XL_volType)
{
	ADD_LOG("Local_PXL_ARM_GetFXVolFromSummit");
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	    ARM_NOCALCIFWIZ();


        // C variable
	    CCString C_ccy1;
	    CCString C_ccy2;
	    CCString C_cvName;
	    double C_date;
	    CCString C_impOrHist;
	    CCString C_volType;
	    
	    // error
	    static int error;
	    static char* reason = "";

	    XL_readStrCell(XL_ccy1,C_ccy1," ARM_ERR: ccy1: string expected",C_result);
	    XL_readStrCell(XL_ccy2,C_ccy2," ARM_ERR: ccy2: string expected",C_result);
	    XL_readNumCell(XL_date,C_date," ARM_ERR: as of date: date expected",C_result);
	    XL_readStrCellWD(XL_cvName,C_cvName,"MO"," ARM_ERR: curve name: string expected",C_result);
	    XL_readStrCellWD(XL_impOrHist,C_impOrHist,"FXVOL"," ARM_ERR: implied or historical vol: string expected",C_result);
	    XL_readStrCellWD(XL_volType,C_volType,"ATM"," ARM_ERR: Volatility type: string expected",C_result);

	    C_ccy1.toUpper ();
	    C_ccy2.toUpper ();
	    C_cvName.toUpper ();
	    C_impOrHist.toUpper ();
	    C_volType.toUpper ();
	    
	    long retCode;
	    long objId;
	    
	    CCString curClass = LOCAL_VOL_CURVE_LIN_CLASS;
	    CCString stringId;
	    
	    retCode = ARMLOCAL_GetFXVolFromSummit (C_ccy1,
										       C_ccy2,
										       C_date,
										       C_cvName,
										       C_impOrHist,
										       C_volType,
										       C_result);

	    if ( retCode == ARM_OK )
	    {
		    objId = C_result.getLong ();

		    stringId = LocalMakeObjectId (objId, curClass);
	    }

	    if ( retCode == ARM_OK )
	    {
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_PXL_ARM_GetFXVolFromSummit" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}



__declspec(dllexport) LPXLOPER WINAPI Local_ARM_GetNewFXVolFromSummit(LPXLOPER XL_ccy1,
																	  LPXLOPER XL_ccy2,
																	  LPXLOPER XL_date,
																	  LPXLOPER XL_cvName,
																	  LPXLOPER XL_domZc,
																	  LPXLOPER XL_forZc,
																	  LPXLOPER XL_fxSpot,
																	  LPXLOPER XL_forwards,
																	  LPXLOPER XL_WhatIsInterpolated,
																	  LPXLOPER XL_correctSplineWithLinear,
																	  LPXLOPER XL_isATM)
{
	ADD_LOG("Local_ARM_GetNewFXVolFromSummit");
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

    /// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	    ARM_NOCALCIFWIZ();


	    // C variable
	    CCString C_ccy1;
	    CCString C_ccy2;
	    CCString C_cvName;
	    double C_date;

	    CCString C_domZc;
	    long domZcId;
	    CCString C_forZc;
	    long forZcId;
	    double C_fxSpot;
	    
	    CCString C_StrWhatIsInterpolated;
	    long C_WhatIsInterpolated;

	    double C_correctSplineWithLinear;
	    double C_correctSplineWithLinear_default(1.0);

        VECTOR<double> C_Forwards;

	    CCString C_isATM;

	    // error
	    static int error;
	    static char* reason = "";

	    XL_readStrCell(XL_ccy1,C_ccy1," ARM_ERR: ccy1: string expected",C_result);
	    XL_readStrCell(XL_ccy2,C_ccy2," ARM_ERR: ccy2: string expected",C_result);
	    XL_readNumCell(XL_date,C_date," ARM_ERR: as of date: date expected",C_result);
	    XL_readStrCellWD(XL_cvName,C_cvName,"MO"," ARM_ERR: curve name: string expected",C_result);
	    XL_readStrCellWD(XL_domZc,C_domZc,"DEFAULT"," ARM_ERR: Domestic Zc Curve: object expected",C_result);
	    XL_readStrCellWD(XL_WhatIsInterpolated, C_StrWhatIsInterpolated, "DEFAULT",
                           " ARM_ERR: whatIsInterpoled. Type: (STRIKE, DELTA) expected", C_result);
	    XL_readNumCellWD(XL_correctSplineWithLinear, C_correctSplineWithLinear,C_correctSplineWithLinear_default,
                       " ARM_ERR: correctSplineWithLinear: number expected",C_result);

        if (( C_WhatIsInterpolated = ARM_ConvWhatIsInterp(C_StrWhatIsInterpolated, C_result)) == ARM_DEFAULT_ERR )
	    {
	       ARM_ARG_ERR();

	       return((LPXLOPER) &XL_result);
	    }

	    if (C_domZc == "DEFAULT")
	    {
	        XL_readNumVector(XL_forwards,C_Forwards," ARM_ERR: Forwards Values: array of numeric expected",C_result);
	    }
	    else
	    {
		    XL_readStrCell(XL_forZc,C_forZc," ARM_ERR: Foreign Zc Curve: object expected",C_result);
		    XL_readNumCell(XL_fxSpot,C_fxSpot," ARM_ERR: Forex Spot: double expected",C_result);
	    }

	    XL_readStrCellWD(XL_isATM, C_isATM, "N",
                           " ARM_ERR: isATM : string expected", C_result);

	    long isATMId;

        if (( isATMId = ARM_ConvYesOrNo(C_isATM, C_result)) == ARM_DEFAULT_ERR )
	    {
	       ARM_ARG_ERR();

	       return((LPXLOPER) &XL_result);
	    }

	    C_ccy1.toUpper ();
	    C_ccy2.toUpper ();
	    C_cvName.toUpper ();

	    if(C_domZc == "DEFAULT")
	    {
		    domZcId = ARM_NULL_OBJECT;
		    forZcId = ARM_NULL_OBJECT;
	    }
	    else
	    {
		    domZcId = LocalGetNumObjectId (C_domZc);
		    forZcId = LocalGetNumObjectId (C_forZc);
	    }

	    long retCode;
	    long objId;
	    CCString prevClass;

    //	CCString curClass = LOCAL_VOLAT_FX_CLASS;
	    CCString curClass;
	    CCString stringId = GetLastCurCellEnvValue ();
	    
	    if (!stringId)
	    {
		    retCode = ARMLOCAL_GetNewFXVolFromSummit (C_ccy1,
												      C_ccy2,
												      C_date,
												      C_cvName,
												      domZcId,
												      forZcId,
												      C_fxSpot,
												      C_Forwards,
												      (long)C_WhatIsInterpolated,
												      (long)C_correctSplineWithLinear,
												      isATMId,
												      curClass,
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

		    if (prevClass == LOCAL_VOLAT_FX_CLASS)
		    {
			    retCode = ARMLOCAL_GetNewFXVolFromSummit (C_ccy1,
													      C_ccy2,
													      C_date,
													      C_cvName,
													      domZcId,
													      forZcId,
													      C_fxSpot,
													      C_Forwards,
													      (long)C_WhatIsInterpolated,
													      (long)C_correctSplineWithLinear,
													      isATMId,
													      curClass,
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

			    retCode = ARMLOCAL_GetNewFXVolFromSummit (C_ccy1,
													      C_ccy2,
													      C_date,
													      C_cvName,
													      domZcId,
													      forZcId,
													      C_fxSpot,
													      C_Forwards,
													      (long)C_WhatIsInterpolated,
													      (long)C_correctSplineWithLinear,													  
													      isATMId,
													      curClass,
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_ARM_GetNewFXVolFromSummit" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}



__declspec(dllexport) LPXLOPER WINAPI Local_ARM_GetInitialFXVolFromSummit(LPXLOPER XL_ccy1,
																		  LPXLOPER XL_ccy2,
																		  LPXLOPER XL_date,
																		  LPXLOPER XL_cvName,
																		  LPXLOPER XL_impOrHist,
																		  LPXLOPER XL_volType)
{
	ADD_LOG("Local_ARM_GetInitialFXVolFromSummit");
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
	    CCString C_ccy1;
	    CCString C_ccy2;
	    CCString C_cvName;
	    double C_date;
	    CCString C_impOrHist;
	    CCString C_volType;

	    VECTOR<CCString> yearterm;
	    VECTOR<double> strike;
	    VECTOR<double> vol;
	    
	    // error
	    static int error;
	    static char* reason = "";

	    XL_readStrCell(XL_ccy1,C_ccy1," ARM_ERR: ccy1: string expected",C_result);
	    XL_readStrCell(XL_ccy2,C_ccy2," ARM_ERR: ccy2: string expected",C_result);
	    XL_readNumCell(XL_date,C_date," ARM_ERR: as of date: date expected",C_result);
	    XL_readStrCellWD(XL_cvName,C_cvName,"MO"," ARM_ERR: curve name: string expected",C_result);
	    XL_readStrCellWD(XL_impOrHist,C_impOrHist,"FXVOL"," ARM_ERR: implied or historical vol: string expected",C_result);
	    XL_readStrCellWD(XL_volType,C_volType,"ATM"," ARM_ERR: Volatility type: string expected",C_result);

	    C_ccy1.toUpper ();
	    C_ccy2.toUpper ();
	    C_cvName.toUpper ();
	    C_impOrHist.toUpper ();
	    C_volType.toUpper ();

	    VECTOR<CCString>* maturities = new VECTOR<CCString>;
	    VECTOR<double>* tenors = new VECTOR<double>;
	    VECTOR<double>* volatilities = new VECTOR<double>;

	    long retCode;

	    retCode = ARMLOCAL_GetInitialFXVolFromSummit (C_ccy1,
												      C_ccy2,
												      C_date,
												      C_cvName,
												      C_impOrHist,
												      C_volType,
												      maturities,
												      tenors,
												      volatilities,
												      C_result);

	    if ( retCode == ARM_OK )
	    {
		    int nbcolumns;
		    
		    if ( strcmp(C_volType,"ATM") == 0 )
		       nbcolumns = 2;
		    else
		       nbcolumns = tenors->size()+1;
		    
		    int nbrows = maturities->size()+1;

		    FreeCurCellErr ();

            XL_result.xltype = xltypeMulti;
		    XL_result.val.array.columns = nbcolumns;
		    XL_result.val.array.rows = nbrows; 
		    XL_result.val.array.lparray = pxArray 
                = (LPXLOPER) GlobalAlloc(GMEM_ZEROINIT, nbrows*nbcolumns*sizeof(XLOPER));

		    // 1ere ligne : intitulé des tenors

		    pxArray[XL_Coordonnate2Rank (0, 0, nbcolumns)].xltype = xltypeStr;
		    pxArray[XL_Coordonnate2Rank (0, 0, nbcolumns)].val.str = XL_StrC2StrPascal ("");
		    pxArray[XL_Coordonnate2Rank (0, 0, nbcolumns)].xltype |= xlbitDLLFree;

		    int i;

		    if ( strcmp(C_volType,"ATM") != 0 )
		    {
			    for (i = 1; i < nbcolumns; i++)
			    {
				    pxArray[XL_Coordonnate2Rank (0, i, nbcolumns)].xltype = xltypeNum;
				    pxArray[XL_Coordonnate2Rank (0, i, nbcolumns)].val.num = (*tenors)[i-1];
			    }
		    }
		    else
		    {
			    pxArray[XL_Coordonnate2Rank (0, 1, nbcolumns)].xltype = xltypeNum;
			    pxArray[XL_Coordonnate2Rank (0, 1, nbcolumns)].val.num = 0.0;
		    }

		    for(i = 1; i < nbrows; i++)
		    {
			    // 1ere colonne : yearTerm
			    pxArray[XL_Coordonnate2Rank (i, 0, nbcolumns)].xltype = xltypeStr;
			    pxArray[XL_Coordonnate2Rank (i, 0, nbcolumns)].val.str = XL_StrC2StrPascal((*maturities)[i-1]);
			    pxArray[XL_Coordonnate2Rank (i, 0, nbcolumns)].xltype |= xlbitDLLFree;

			    for (int j = 1; j < nbcolumns; j++)
			    {
				    pxArray[XL_Coordonnate2Rank (i, j, nbcolumns)].xltype = xltypeNum;
				    pxArray[XL_Coordonnate2Rank (i, j, nbcolumns)].val.num = (*volatilities)[(j-1)*(nbrows-1)+(i-1)];
			    }
		    }
	    }
	    else
	    {
		    ARM_ERR();
	    }

	    if (maturities)
		    delete maturities;
	    maturities = NULL;

	    if (tenors)
		    delete tenors;
	    tenors = NULL;

	    if (volatilities)
		    delete volatilities;
	    volatilities = NULL;
    }
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_ARM_GetInitialFXVolFromSummit" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}



__declspec(dllexport) LPXLOPER WINAPI Local_ARM_GetInitialVolFromSummit(LPXLOPER XL_index,
																		LPXLOPER XL_currency,
																		LPXLOPER XL_cvName,
																		LPXLOPER XL_date,
																		LPXLOPER XL_vtype,
																		LPXLOPER XL_matIndex)
{
	ADD_LOG("Local_ARM_GetInitialVolFromSummit");
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
	    double C_date;
	    CCString C_vtype;
	    CCString C_matIndex;

	    VECTOR<CCString> yearterm;
	    VECTOR<CCString> tenor;
	    VECTOR<double> strike;
	    VECTOR<double> vol;
	    
	    // error
	    static int error;
	    static char* reason = "";

	    XL_readStrCell(XL_index,C_index," ARM_ERR: index: string expected",C_result);
	    XL_readStrCell(XL_currency,C_currency," ARM_ERR: currency: string expected",C_result);
	    XL_readStrCell(XL_cvName,C_cvName," ARM_ERR: curve name: string expected",C_result);
	    XL_readNumCell(XL_date,C_date," ARM_ERR: as of date: date expected",C_result);
	    XL_readStrCell(XL_vtype,C_vtype," ARM_ERR: volatility type: string expected",C_result);
	    XL_readStrCellWD(XL_matIndex,C_matIndex,"ATM"," ARM_ERR: maturity index: string expected",C_result);

	    C_index.toUpper ();
	    C_currency.toUpper ();
	    C_cvName.toUpper ();
	    C_vtype.toUpper ();
	    C_matIndex.toUpper ();
	    
	    long retCode;

	    VECTOR<CCString>* maturities = new VECTOR<CCString>;
	    VECTOR<CCString>* tenors = new VECTOR<CCString>;
	    VECTOR<double>* volatilities = new VECTOR<double>;

	    retCode = ARMLOCAL_GetInitialVolFromSummit (C_index,
												    C_currency,
												    C_cvName,
												    C_date,
												    C_vtype,
												    C_matIndex,
												    maturities,
												    tenors,
												    volatilities,
												    C_result);

	    if(retCode == ARM_OK)
	    {
			// on sette le IsEtk dans C_result
		    if (C_result.getLong() >= ETKRETRIEVER)
		    {
			    int nbcolumns = tenors->size () + 1;
			    int nbrows = maturities->size () + 1;

			    FreeCurCellErr ();
			    XL_result.xltype = xltypeMulti;
			    XL_result.val.array.columns = nbcolumns;
			    XL_result.val.array.rows = nbrows; 
			    XL_result.val.array.lparray = pxArray = (LPXLOPER)GlobalAlloc (GMEM_ZEROINIT, nbrows * nbcolumns * sizeof (XLOPER));

			    // 1ere ligne : intitulé des tenors

			    pxArray[XL_Coordonnate2Rank (0, 0, nbcolumns)].xltype = xltypeStr;
			    pxArray[XL_Coordonnate2Rank (0, 0, nbcolumns)].val.str = XL_StrC2StrPascal ("");
			    pxArray[XL_Coordonnate2Rank (0, 0, nbcolumns)].xltype |= xlbitDLLFree;

			    int i;

			    if (strcmp(C_matIndex,"ATM") == 0)
			    {
				    for (i = 1; i < nbcolumns; i++)
				    {
					    pxArray[XL_Coordonnate2Rank (0, i, nbcolumns)].xltype = xltypeStr;
					    pxArray[XL_Coordonnate2Rank (0, i, nbcolumns)].val.str = XL_StrC2StrPascal ((*tenors)[i-1]);
					    pxArray[XL_Coordonnate2Rank (0, i, nbcolumns)].xltype |= xlbitDLLFree;
				    }
			    }
			    else
			    {
				    for (i = 1; i < nbcolumns; i++)
				    {
					    pxArray[XL_Coordonnate2Rank (0, i, nbcolumns)].xltype = xltypeNum;
					    pxArray[XL_Coordonnate2Rank (0, i, nbcolumns)].val.num = atof ((*tenors)[i-1]) / 100.;
				    }
			    }

			    for(i = 1; i < nbrows; i++)
			    {
				    // 1ere colonne : yearTerm
				    pxArray[XL_Coordonnate2Rank (i, 0, nbcolumns)].xltype = xltypeStr;
				    pxArray[XL_Coordonnate2Rank (i, 0, nbcolumns)].val.str = XL_StrC2StrPascal ((*maturities)[i-1]);
				    pxArray[XL_Coordonnate2Rank (i, 0, nbcolumns)].xltype |= xlbitDLLFree;

				    for (int j = 1; j < nbcolumns; j++)
				    {
					    pxArray[XL_Coordonnate2Rank (i, j, nbcolumns)].xltype = xltypeNum;
					    pxArray[XL_Coordonnate2Rank (i, j, nbcolumns)].val.num = (*volatilities)[(j-1)+(i-1)*(nbcolumns-1)];
				    }
			    }
		    }
		    else
		    {
			    long result;
			    if (strcmp((const char*)C_matIndex,"ATM") == 0)
				    result = LocalExtractVolFromFile (C_result.getString(), yearterm, tenor, vol);
			    else
				    result = LocalExtractSmileFromFile (C_result.getString(), yearterm, strike, vol);

			    if(result == ARM_OK)
			    {
				    int nbcolumns;
				    int nbrows = yearterm.size () + 1;
				    if (strcmp((const char*)C_matIndex,"ATM") == 0)
					    nbcolumns = tenor.size () + 1;
				    else
					    nbcolumns = strike.size () + 1;

				    FreeCurCellErr ();
				    XL_result.xltype = xltypeMulti;
				    XL_result.val.array.columns = nbcolumns;
				    XL_result.val.array.rows = nbrows; 
				    XL_result.val.array.lparray = pxArray = (LPXLOPER)GlobalAlloc (GMEM_ZEROINIT, nbrows * nbcolumns * sizeof (XLOPER));

				    // 1ere ligne : intitulé des tenors

				    pxArray[XL_Coordonnate2Rank (0, 0, nbcolumns)].xltype = xltypeStr;
				    pxArray[XL_Coordonnate2Rank (0, 0, nbcolumns)].val.str = XL_StrC2StrPascal ("");
				    pxArray[XL_Coordonnate2Rank (0, 0, nbcolumns)].xltype |= xlbitDLLFree;

				    for (int i = 1; i < nbcolumns; i++)
				    {
					    if (strcmp((const char*)C_matIndex,"ATM") == 0)
					    {
						    pxArray[XL_Coordonnate2Rank (0, i, nbcolumns)].xltype = xltypeStr;
						    pxArray[XL_Coordonnate2Rank (0, i, nbcolumns)].val.str = XL_StrC2StrPascal (tenor[i-1]);
						    pxArray[XL_Coordonnate2Rank (0, i, nbcolumns)].xltype |= xlbitDLLFree;
					    }
					    else
					    {
						    pxArray[XL_Coordonnate2Rank (0, i, nbcolumns)].xltype = xltypeNum;
						    pxArray[XL_Coordonnate2Rank (0, i, nbcolumns)].val.num = strike[i-1];
					    }
				    }

				    for(i = 1; i < nbrows; i++)
				    {
					    // 1ere colonne : yearTerm
					    pxArray[XL_Coordonnate2Rank (i, 0, nbcolumns)].xltype = xltypeStr;
					    pxArray[XL_Coordonnate2Rank (i, 0, nbcolumns)].val.str = XL_StrC2StrPascal (yearterm[i-1]);
					    pxArray[XL_Coordonnate2Rank (i, 0, nbcolumns)].xltype |= xlbitDLLFree;

					    for (int j = 1; j < nbcolumns; j++)
					    {
						    pxArray[XL_Coordonnate2Rank (i, j, nbcolumns)].xltype = xltypeNum;
						    if (strcmp((const char*)C_matIndex,"ATM") == 0)
							    pxArray[XL_Coordonnate2Rank (i, j, nbcolumns)].val.num = vol[(j-1)*(nbrows-1)+(i-1)];
						    else
							    pxArray[XL_Coordonnate2Rank (i, j, nbcolumns)].val.num = vol[(j-1)+(i-1)*(nbcolumns-1)];
					    }
				    }
			    }
		    }
	    }
	    else
	    {
		    ARM_ERR();
	    }

	    if (maturities)
		    delete maturities;
	    maturities = NULL;

	    if (tenors)
		    delete tenors;
	    tenors = NULL;

	    if (volatilities)
		    delete volatilities;
	    volatilities = NULL;
    }
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_ARM_GetInitialVolFromSummit" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}



__declspec(dllexport) LPXLOPER WINAPI Local_ARM_GetFXCorrelFromSummit(LPXLOPER XL_ccy1,
																	  LPXLOPER XL_index,
																	  LPXLOPER XL_ccy2,
																	  LPXLOPER XL_asof,
																	  LPXLOPER XL_cvName,
																	  LPXLOPER XL_tenors)
{
	ADD_LOG("Local_ARM_GetFXCorrelFromSummit");
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

    /// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	    ARM_NOCALCIFWIZ();


	    // C variable
	    CCString C_ccy1;
	    CCString C_index;
	    CCString C_ccy2;
	    CCString C_cvName;
	    double C_date;
	    VECTOR<CCString> C_tenors;
	    
	    // error
	    static int error;
	    static char* reason = "";

	    XL_readStrCell(XL_ccy1,C_ccy1," ARM_ERR: ccy1: string expected",C_result);
	    XL_readStrCell(XL_index,C_index," ARM_ERR: index: string expected",C_result);
	    XL_readStrCell(XL_ccy2,C_ccy2," ARM_ERR: ccy2: string expected",C_result);
	    XL_readStrCell(XL_cvName,C_cvName," ARM_ERR: curve name: string expected",C_result);
	    XL_readNumCell(XL_asof,C_date," ARM_ERR: as of date: date expected",C_result);
	    XL_readStrVector(XL_tenors,C_tenors," ARM_ERR: tenors (3M, 1Y, ...: array of numeric expected",XL_TYPE_STRING,C_result);

	    C_ccy1.toUpper ();
	    C_index.toUpper ();
	    C_ccy2.toUpper ();
	    C_cvName.toUpper ();

	    for (int i = 0; i < C_tenors.size(); i++)
		    C_tenors[i].toUpper();

	    long retCode;
	    long objId;
	    CCString prevClass;
	    
	    CCString curClass = LOCAL_VOL_CURVE_LIN_CLASS;
	    CCString stringId = GetLastCurCellEnvValue ();
	    
	    if (!stringId)
	    {
		    retCode = ARMLOCAL_GetFXCorrelFromSummit (C_ccy1,
												      C_index,
												      C_ccy2,
												      C_date,
												      C_cvName,
												      C_tenors,
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
			    retCode = ARMLOCAL_GetFXCorrelFromSummit (C_ccy1,
													      C_index,
													      C_ccy2,
													      C_date,
													      C_cvName,
													      C_tenors,
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

			    retCode = ARMLOCAL_GetFXCorrelFromSummit (C_ccy1,
													      C_index,
													      C_ccy2,
													      C_date,
													      C_cvName,
													      C_tenors,
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_ARM_GetFXCorrelFromSummit" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}



__declspec(dllexport) LPXLOPER WINAPI Local_PXL_ARM_GetFXCorrelFromSummit(LPXLOPER XL_ccy1,
																		  LPXLOPER XL_index,
																		  LPXLOPER XL_ccy2,
																		  LPXLOPER XL_asof,
																		  LPXLOPER XL_cvName,
																		  LPXLOPER XL_tenors)
{
	ADD_LOG("Local_PXL_ARM_GetFXCorrelFromSummit");
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

    /// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	    ARM_NOCALCIFWIZ();

	    // C variable
	    CCString C_ccy1;
	    CCString C_index;
	    CCString C_ccy2;
	    CCString C_cvName;
	    double C_date;
	    VECTOR<CCString> C_tenors;
	    
	    // error
	    static int error;
	    static char* reason = "";

	    XL_readStrCell(XL_ccy1,C_ccy1," ARM_ERR: ccy1: string expected",C_result);
	    XL_readStrCell(XL_index,C_index," ARM_ERR: index: string expected",C_result);
	    XL_readStrCell(XL_ccy2,C_ccy2," ARM_ERR: ccy2: string expected",C_result);
	    XL_readStrCell(XL_cvName,C_cvName," ARM_ERR: curve name: string expected",C_result);
	    XL_readNumCell(XL_asof,C_date," ARM_ERR: as of date: date expected",C_result);
	    XL_readStrVector(XL_tenors,C_tenors," ARM_ERR: tenors (3M, 1Y, ...: array of numeric expected",XL_TYPE_STRING,C_result);

	    C_ccy1.toUpper ();
	    C_index.toUpper ();
	    C_ccy2.toUpper ();
	    C_cvName.toUpper ();

	    for (int i = 0; i < C_tenors.size(); i++)
		    C_tenors[i].toUpper();
	    
	    long retCode;
	    long objId;
	    CCString prevClass;
	    
	    CCString curClass = LOCAL_VOL_CURVE_LIN_CLASS;
	    CCString stringId;
	    
	    retCode = ARMLOCAL_GetFXCorrelFromSummit (C_ccy1,
											      C_index,
											      C_ccy2,
											      C_date,
											      C_cvName,
											      C_tenors,
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_PXL_ARM_GetFXCorrelFromSummit" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}



__declspec(dllexport) LPXLOPER WINAPI Local_ARM_GetCorrelFromSummit(LPXLOPER XL_ccy1,
																	LPXLOPER XL_index1,
																	LPXLOPER XL_ccy2,
																	LPXLOPER XL_index2,
																	LPXLOPER XL_asof,
																	LPXLOPER XL_cvName)
{
	ADD_LOG("Local_ARM_GetCorrelFromSummit");
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

    /// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	    ARM_NOCALCIFWIZ();


	    // C variable
	    CCString C_ccy1;
	    CCString C_index1;
	    CCString C_ccy2;
	    CCString C_index2;
	    CCString C_cvName;
	    double C_date;
	    
	    // error
	    static int error;
	    static char* reason = "";

	    XL_readStrCell(XL_ccy1,C_ccy1," ARM_ERR: ccy1: string expected",C_result);
	    XL_readStrCell(XL_index1,C_index1," ARM_ERR: index1: string expected",C_result);
	    XL_readStrCell(XL_ccy2,C_ccy2," ARM_ERR: ccy2: string expected",C_result);
	    XL_readStrCell(XL_index2,C_index2," ARM_ERR: index2: string expected",C_result);
	    XL_readStrCell(XL_cvName,C_cvName," ARM_ERR: curve name: string expected",C_result);
	    XL_readNumCell(XL_asof,C_date," ARM_ERR: as of date: date expected",C_result);

	    C_ccy1.toUpper ();
	    C_index1.toUpper ();
	    C_ccy2.toUpper ();
	    C_index2.toUpper ();
	    C_cvName.toUpper ();

	    long retCode;
	    long objId;
	    CCString prevClass;
	    
	    CCString curClass = LOCAL_VOL_CURVE_LIN_CLASS;
	    CCString stringId = GetLastCurCellEnvValue ();
	    
	    if (!stringId)
	    {
		    retCode = ARMLOCAL_GetCorrelFromSummit (C_ccy1,
												    C_index1,
												    C_ccy2,
												    C_index2,
												    C_date,
												    C_cvName,
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
			    retCode = ARMLOCAL_GetCorrelFromSummit (C_ccy1,
													    C_index1,
													    C_ccy2,
													    C_index2,
													    C_date,
													    C_cvName,
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

			    retCode = ARMLOCAL_GetCorrelFromSummit (C_ccy1,
													    C_index1,
													    C_ccy2,
													    C_index2,
													    C_date,
													    C_cvName,
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_ARM_GetCorrelFromSummit" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}



__declspec(dllexport) LPXLOPER WINAPI Local_PXL_ARM_GetCorrelFromSummit(LPXLOPER XL_ccy1,
																		LPXLOPER XL_index1,
																		LPXLOPER XL_ccy2,
																		LPXLOPER XL_index2,
																		LPXLOPER XL_asof,
																		LPXLOPER XL_cvName)
{
	ADD_LOG("Local_PXL_ARM_GetCorrelFromSummit");
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

    /// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	    ARM_NOCALCIFWIZ();


	    // C variable
	    CCString C_ccy1;
	    CCString C_index1;
	    CCString C_ccy2;
	    CCString C_index2;
	    CCString C_cvName;
	    double C_date;
	    
	    // error
	    static int error;
	    static char* reason = "";

	    XL_readStrCell(XL_ccy1,C_ccy1," ARM_ERR: ccy1: string expected",C_result);
	    XL_readStrCell(XL_index1,C_index1," ARM_ERR: index1: string expected",C_result);
	    XL_readStrCell(XL_ccy2,C_ccy2," ARM_ERR: ccy2: string expected",C_result);
	    XL_readStrCell(XL_index2,C_index2," ARM_ERR: index2: string expected",C_result);
	    XL_readStrCell(XL_cvName,C_cvName," ARM_ERR: curve name: string expected",C_result);
	    XL_readNumCell(XL_asof,C_date," ARM_ERR: as of date: date expected",C_result);
	    
	    C_ccy1.toUpper ();
	    C_index1.toUpper ();
	    C_ccy2.toUpper ();
	    C_index2.toUpper ();
	    C_cvName.toUpper ();

	    long retCode;
	    long objId;
	    CCString prevClass;
	    
	    CCString curClass = LOCAL_VOL_CURVE_LIN_CLASS;
	    CCString stringId;
	    
	    retCode = ARMLOCAL_GetCorrelFromSummit (C_ccy1,
											    C_index1,
											    C_ccy2,
											    C_index2,
											    C_date,
											    C_cvName,
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_PXL_ARM_GetCorrelFromSummit" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}




__declspec(dllexport) LPXLOPER WINAPI Local_ARM_GetNthMaturity(LPXLOPER XL_curve,
															   LPXLOPER XL_nLine)
{
	ADD_LOG("Local_ARM_GetNthMaturity");
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
	    double C_nLine;
	    double C_nLine_default = -1.0;

	    // error
	    static int error;
	    static char* reason = "";

	    XL_readStrCell(XL_curve,C_curve," ARM_ERR: string expected",C_result);
	    XL_readNumCellWD(XL_nLine,C_nLine,C_nLine_default," ARM_ERR: nth Line: numeric expected",C_result);

	    long retCode = ARMLOCAL_GetNthMaturity (LocalGetNumObjectId (C_curve), (long)C_nLine, C_result);

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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_ARM_GetNthMaturity" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}



__declspec(dllexport) LPXLOPER WINAPI Local_ARM_CONVERTFROMBSTONORMALVOL (LPXLOPER XL_vol,
																		  LPXLOPER XL_zc,
																		  LPXLOPER XL_isSwoptVol,
																		  LPXLOPER XL_inPct)
{
	ADD_LOG("Local_ARM_CONVERTFROMBSTONORMALVOL ");
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	    ARM_NOCALCIFWIZ();


	    // C variable
	    CCString C_vol;
	    CCString C_zc;

	    double C_isSwoptVol;
	    double C_isSwoptVol_default = 1.;

	    double C_inPct;
	    double C_inPct_default = 1.;

	    // error
	    static int error;
	    static char* reason = "";

	    XL_readStrCell(XL_vol,C_vol," ARM_ERR: volatility: object expected",C_result);
	    XL_readStrCell(XL_zc,C_zc," ARM_ERR: zc curve: object expected",C_result);
	    XL_readNumCellWD(XL_isSwoptVol,C_isSwoptVol,C_isSwoptVol_default," ARM_ERR: is Swopt Vol?: numeric expected",C_result);
	    XL_readNumCellWD(XL_inPct,C_inPct,C_inPct_default," ARM_ERR: in percentage: numeric expected",C_result);

	    long retCode;
	    long objId;
	    CCString prevClass;
	    
	    CCString curClass = LOCAL_VOL_CURVE_LIN_CLASS;
	    CCString stringId = GetLastCurCellEnvValue ();

	    if(!stringId)
	    {
		    retCode = ARMLOCAL_ARM_CONVERTFROMBSTONORMALVOL (LocalGetNumObjectId(C_vol),
														     LocalGetNumObjectId(C_zc),
														     (long)C_isSwoptVol,
														     (long)C_inPct,
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
			    retCode = ARMLOCAL_ARM_CONVERTFROMBSTONORMALVOL (LocalGetNumObjectId(C_vol),
															     LocalGetNumObjectId(C_zc),
															     (long)C_isSwoptVol,
															     (long)C_inPct,
															     C_result,
															     objId);

			    if (retCode == ARM_OK)
			    {
				    LocalSetCurCellEnvValue (curClass, objId); 

				    stringId = LocalMakeObjectId (objId, curClass);
			    }
		    }
		    else
		    {
			    FreeCurCellContent ();

			    retCode = ARMLOCAL_ARM_CONVERTFROMBSTONORMALVOL (LocalGetNumObjectId(C_vol),
															     LocalGetNumObjectId(C_zc),
															     (long)C_isSwoptVol,
															     (long)C_inPct,
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_ARM_CONVERTFROMBSTONORMALVOL" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}



__declspec(dllexport) LPXLOPER WINAPI Local_PXL_ARM_CONVERTFROMBSTONORMALVOL (LPXLOPER XL_vol,
																			  LPXLOPER XL_zc,
																			  LPXLOPER XL_isSwoptVol,
																			  LPXLOPER XL_inPct)
{
	ADD_LOG("Local_PXL_ARM_CONVERTFROMBSTONORMALVOL ");
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	    ARM_NOCALCIFWIZ();


	    // C variable
	    CCString C_vol;
	    CCString C_zc;

	    double C_isSwoptVol;
	    double C_isSwoptVol_default = 1.;

	    double C_inPct;
	    double C_inPct_default = 1.;

	    // error
	    static int error;
	    static char* reason = "";

	    XL_readStrCell(XL_vol,C_vol," ARM_ERR: volatility: object expected",C_result);
	    XL_readStrCell(XL_zc,C_zc," ARM_ERR: zc curve: object expected",C_result);
	    XL_readNumCellWD(XL_isSwoptVol,C_isSwoptVol,C_isSwoptVol_default," ARM_ERR: is Swopt Vol?: numeric expected",C_result);
	    XL_readNumCellWD(XL_inPct,C_inPct,C_inPct_default," ARM_ERR: in percentage: numeric expected",C_result);

	    long retCode;
	    long objId;
	    CCString prevClass;
	    
	    CCString curClass = LOCAL_VOL_CURVE_LIN_CLASS;
	    CCString stringId;

	    retCode = ARMLOCAL_ARM_CONVERTFROMBSTONORMALVOL (LocalGetNumObjectId(C_vol),
													     LocalGetNumObjectId(C_zc),
													     (long)C_isSwoptVol,
													     (long)C_inPct,
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_PXL_ARM_CONVERTFROMBSTONORMALVOL" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}



__declspec(dllexport) LPXLOPER WINAPI Local_ARM_CONVERTFROMNORMALTOBSVOL (LPXLOPER XL_vol,
																		  LPXLOPER XL_zc,
																		  LPXLOPER XL_isSwoptVol,
																		  LPXLOPER XL_inPct,
																		  LPXLOPER XL_post100)
{
	ADD_LOG("Local_ARM_CONVERTFROMNORMALTOBSVOL ");
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	    ARM_NOCALCIFWIZ();


	    // C variable
	    CCString C_vol;
	    CCString C_zc;

	    double C_isSwoptVol;
	    double C_isSwoptVol_default = 1.;

	    double C_inPct;
	    double C_inPct_default = 1.;

	    double C_outPct;
	    double C_outPct_default = 1.;

	    // error
	    static int error;
	    static char* reason = "";

	    XL_readStrCell(XL_vol,C_vol," ARM_ERR: volatility: object expected",C_result);
	    XL_readStrCell(XL_zc,C_zc," ARM_ERR: zc curve: object expected",C_result);
	    XL_readNumCellWD(XL_isSwoptVol,C_isSwoptVol,C_isSwoptVol_default," ARM_ERR: is Swopt Vol?: numeric expected",C_result);
	    XL_readNumCellWD(XL_inPct,C_inPct,C_inPct_default," ARM_ERR: in percentage: numeric expected",C_result);
	    XL_readNumCellWD(XL_post100,C_outPct,C_outPct_default," ARM_ERR: post 100: numeric expected",C_result);

	    long retCode;
	    long objId;
	    CCString prevClass;
	    
	    CCString curClass = LOCAL_VOL_CURVE_LIN_CLASS;
	    CCString stringId = GetLastCurCellEnvValue ();

	    if(!stringId)
	    {
		    retCode = ARMLOCAL_ARM_CONVERTFROMNORMALTOBSVOL (LocalGetNumObjectId(C_vol),
														     LocalGetNumObjectId(C_zc),
														     (long)C_isSwoptVol,
														     (long)C_inPct,
														     (long)C_outPct,
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
			    retCode = ARMLOCAL_ARM_CONVERTFROMNORMALTOBSVOL (LocalGetNumObjectId(C_vol),
															     LocalGetNumObjectId(C_zc),
															     (long)C_isSwoptVol,
															     (long)C_inPct,
															     (long)C_outPct,
															     C_result,
															     objId);

			    if (retCode == ARM_OK)
			    {
				    LocalSetCurCellEnvValue (curClass, objId); 

				    stringId = LocalMakeObjectId (objId, curClass);
			    }
		    }
		    else
		    {
			    FreeCurCellContent ();

			    retCode = ARMLOCAL_ARM_CONVERTFROMNORMALTOBSVOL (LocalGetNumObjectId(C_vol),
															     LocalGetNumObjectId(C_zc),
															     (long)C_isSwoptVol,
															     (long)C_inPct,
															     (long)C_outPct,
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_ARM_CONVERTFROMNORMALTOBSVOL" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}



__declspec(dllexport) LPXLOPER WINAPI Local_PXL_ARM_CONVERTFROMNORMALTOBSVOL (LPXLOPER XL_vol,
																			  LPXLOPER XL_zc,
																			  LPXLOPER XL_isSwoptVol,
																			  LPXLOPER XL_inPct,
																			  LPXLOPER XL_post100)
{
	ADD_LOG("Local_PXL_ARM_CONVERTFROMNORMALTOBSVOL ");
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	    ARM_NOCALCIFWIZ();


	    // C variable
	    CCString C_vol;
	    CCString C_zc;

	    double C_isSwoptVol;
	    double C_isSwoptVol_default = 1.;

	    double C_inPct;
	    double C_inPct_default = 1.;

	    double C_outPct;
	    double C_outPct_default = 1.;

	    // error
	    static int error;
	    static char* reason = "";

	    XL_readStrCell(XL_vol,C_vol," ARM_ERR: volatility: object expected",C_result);
	    XL_readStrCell(XL_zc,C_zc," ARM_ERR: zc curve: object expected",C_result);
	    XL_readNumCellWD(XL_isSwoptVol,C_isSwoptVol,C_isSwoptVol_default," ARM_ERR: is Swopt Vol?: numeric expected",C_result);
	    XL_readNumCellWD(XL_inPct,C_inPct,C_inPct_default," ARM_ERR: in percentage: numeric expected",C_result);
	    XL_readNumCellWD(XL_post100,C_outPct,C_outPct_default," ARM_ERR: post 100: numeric expected",C_result);

	    long retCode;
	    long objId;
	    CCString prevClass;
	    
	    CCString curClass = LOCAL_VOL_CURVE_LIN_CLASS;
	    CCString stringId;

	    retCode = ARMLOCAL_ARM_CONVERTFROMNORMALTOBSVOL (LocalGetNumObjectId(C_vol),
													     LocalGetNumObjectId(C_zc),
													     (long)C_isSwoptVol,
													     (long)C_inPct,
													     (long)C_outPct,
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_PXL_ARM_CONVERTFROMNORMALTOBSVOL" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}



__declspec(dllexport) LPXLOPER WINAPI Local_ARM_CONV3FFROMSPOTTOFWDVOL_IN_PRCS(LPXLOPER XL_tree3f,
																	           LPXLOPER XL_PrcsObject,
                                                                               LPXLOPER XL_ForwardVolDates)
                                                                      
{
	ADD_LOG("Local_ARM_CONV3FFROMSPOTTOFWDVOL_IN_PRCS");
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	    ARM_NOCALCIFWIZ();


	    // C variable
	    CCString C_tree3f;
	    CCString C_PrcsObject;
        VECTOR<double> C_ForwardVolDates;
        VECTOR<double> C_ForwardVolDates_Default;
    
	    
	    // error
	    static int error;
	    static char* reason = "";

	    XL_readStrCell(XL_tree3f,C_tree3f," ARM_ERR: tree 3f model: object expected",C_result);
        XL_readStrCell(XL_PrcsObject,C_PrcsObject," ARM_ERR: PRCS Object: object expected",C_result);
        XL_readNumVectorWD(XL_ForwardVolDates,C_ForwardVolDates,C_ForwardVolDates_Default," ARM_ERR: Forward Vol dates: array of numeric expected",C_result);
   
	    long retCode;
	    long objId;
	    CCString prevClass;
	    
	    CCString curClass = LOCAL_VOL_CURVE_LIN_CLASS;
	    CCString stringId = GetLastCurCellEnvValue ();

	    if(!stringId)
	    {
		    retCode = ARMLOCAL_ARM_CONV3FFROMSPOTTOFWDVOL(  LocalGetNumObjectId(C_tree3f),
													        LocalGetNumObjectId(C_PrcsObject),
                                                            C_ForwardVolDates,
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
			    retCode = ARMLOCAL_ARM_CONV3FFROMSPOTTOFWDVOL(  LocalGetNumObjectId(C_tree3f),
													            LocalGetNumObjectId(C_PrcsObject),
                                                                C_ForwardVolDates,
													            C_result,
                                                                objId);

			    if (retCode == ARM_OK)
			    {
				    LocalSetCurCellEnvValue (curClass, objId); 

				    stringId = LocalMakeObjectId (objId, curClass);
			    }
		    }
		    else
		    {
			    FreeCurCellContent ();

			    retCode = ARMLOCAL_ARM_CONV3FFROMSPOTTOFWDVOL(  LocalGetNumObjectId(C_tree3f),
													            LocalGetNumObjectId(C_PrcsObject),
                                                                C_ForwardVolDates,
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_ARM_CONV3FFROMSPOTTOFWDVOL_IN_PRCS" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}



__declspec(dllexport) LPXLOPER WINAPI Local_PXL_ARM_CONV3FFROMSPOTTOFWDVOL_IN_PRCS(LPXLOPER XL_tree3f,
																	               LPXLOPER XL_PrcsObject,
                                                                                   LPXLOPER XL_ForwardVolDates)
                                                                            
{
	ADD_LOG("Local_PXL_ARM_CONV3FFROMSPOTTOFWDVOL_IN_PRCS");
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	    ARM_NOCALCIFWIZ();


	    // C variable
	    CCString C_tree3f;
        CCString C_PrcsObject;
        VECTOR<double> C_ForwardVolDates;
        VECTOR<double> C_ForwardVolDates_Default;
	    

	    // error
	    static int error;
	    static char* reason = "";

	    XL_readStrCell(XL_tree3f,C_tree3f," ARM_ERR: tree 3f model: object expected",C_result);
        XL_readStrCell(XL_PrcsObject,C_PrcsObject," ARM_ERR: PRCS Object: object expected",C_result);
        XL_readNumVectorWD(XL_ForwardVolDates,C_ForwardVolDates,C_ForwardVolDates_Default," ARM_ERR: Forward Vol dates: array of numeric expected",C_result);
   

	    long retCode;
	    long objId;
	    CCString prevClass;
	    
	    CCString curClass = LOCAL_VOL_CURVE_LIN_CLASS;
	    CCString stringId;

	    retCode = ARMLOCAL_ARM_CONV3FFROMSPOTTOFWDVOL(LocalGetNumObjectId(C_tree3f),
													  LocalGetNumObjectId(C_PrcsObject),
                                                      C_ForwardVolDates,
													  C_result);

	    if ( retCode == ARM_OK)
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_PXL_ARM_CONV3FFROMSPOTTOFWDVOL_IN_PRCS" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}



__declspec(dllexport) LPXLOPER WINAPI Local_InterpolInStrikeFwdTime (LPXLOPER XL_curve,
																	 LPXLOPER XL_forward,
																	 LPXLOPER XL_strike,
																	 LPXLOPER XL_matu,
																	 LPXLOPER XL_precision,
																	 LPXLOPER XL_sigmaATM,
																	 LPXLOPER XL_y2NULL)
{
	ADD_LOG("Local_InterpolInStrikeFwdTime ");
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
	    double C_forward;
	    double C_strike;
	    double C_matu;
	    double C_precision;
	    double C_sigmaATM;
	    double C_y2NULL;
	    double default_y2NULL = 0.0;

	    // error
	    static int error;
	    static char* reason = "";

	    XL_readStrCell(XL_curve,C_curve," ARM_ERR: string expected",C_result);
	    XL_readNumCell(XL_forward,C_forward," ARM_ERR: forward: numeric expected",C_result);
	    XL_readNumCell(XL_strike,C_strike," ARM_ERR: strike: numeric expected",C_result);
	    XL_readNumCell(XL_matu,C_matu," ARM_ERR: maturity: numeric expected",C_result);
	    XL_readNumCell(XL_precision,C_precision," ARM_ERR: precision: numeric expected",C_result);
	    XL_readNumCell(XL_sigmaATM,C_sigmaATM," ARM_ERR: sigmaATM: numeric expected",C_result);
	    XL_readNumCellWD(XL_y2NULL, C_y2NULL, default_y2NULL," ARM_ERR: y2NULL: numeric expected",C_result);

	    long retCode = ARMLOCAL_InterpolInStrikeFwdTime(LocalGetNumObjectId (C_curve),
													    C_forward,
													    C_strike,		
													    C_matu,
													    C_precision,
													    C_sigmaATM,
													    (long)C_y2NULL,
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_InterpolInStrikeFwdTime" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


__declspec(dllexport) LPXLOPER WINAPI Local_ComputeFxVol (LPXLOPER XL_fxvol,
														  LPXLOPER XL_asof,
														  LPXLOPER XL_matu,
														  LPXLOPER XL_calcmatu,
														  LPXLOPER XL_fxspot,
														  LPXLOPER XL_strike,
														  LPXLOPER XL_discCrv,
														  LPXLOPER XL_divCrv)
{
	ADD_LOG("Local_ComputeFxVol ");
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();


	// C variable
	CCString C_fxvol;
	double C_asof;
	double C_matu;
	double C_calcmatu;
	double C_fxspot;
	double C_strike;
	CCString C_discCrv;
	CCString C_divCrv;

	// error
	static int error;
	static char* reason = "";

	XL_readStrCell(XL_fxvol,C_fxvol," ARM_ERR: string expected",C_result);
	XL_readNumCell(XL_asof,C_asof," ARM_ERR: As of: numeric expected",C_result);
	XL_readNumCell(XL_matu,C_matu," ARM_ERR: maturity: numeric expected",C_result);
	XL_readNumCell(XL_calcmatu,C_calcmatu," ARM_ERR: Calc maturity: numeric expected",C_result);
	XL_readNumCell(XL_fxspot,C_fxspot," ARM_ERR: Forex Spot: numeric expected",C_result);
	XL_readNumCell(XL_strike,C_strike," ARM_ERR: strike: numeric expected",C_result);
	XL_readStrCell(XL_discCrv,C_discCrv," ARM_ERR: string expected",C_result);
	XL_readStrCell(XL_divCrv,C_divCrv," ARM_ERR: string expected",C_result);

	long retCode = ARMLOCAL_ComputeFxVol (	LocalGetNumObjectId (C_fxvol),
											C_asof,
											C_matu,
											C_calcmatu,
											C_fxspot,
											C_strike,
											LocalGetNumObjectId (C_discCrv),
											LocalGetNumObjectId (C_divCrv),
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_ComputeFxVol" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


__declspec(dllexport) LPXLOPER WINAPI Local_ARM_CONV3FFROMSPOTTOFWDVOL( LPXLOPER XL_asOf,
                                                                        LPXLOPER XL_dZcId,
                                                                        LPXLOPER XL_fZcId,
                                                                        LPXLOPER XL_dBSZcId,
                                                                        LPXLOPER XL_fBSZcId,
                                                                        LPXLOPER XL_volSwopBaseId,
                                                                        LPXLOPER XL_volSwopForeignId,
                                                                        LPXLOPER XL_fxVolId,
                                                                        LPXLOPER XL_dMeanReversionBase,
                                                                        LPXLOPER XL_dMeanReversionForeign,
                                                                        LPXLOPER XL_dFxRdCorrId,
                                                                        LPXLOPER XL_dFxRfCorrId,
                                                                        LPXLOPER XL_dRdRfCorrId,
                                                                        LPXLOPER XL_dCutOff,
                                                                        LPXLOPER XL_dVolLongTerm,
                                                                        LPXLOPER XL_calibBasisIncluded,
                                                                        LPXLOPER XL_ForwardVolDates)                                                                    
{
	ADD_LOG("Local_ARM_CONV3FFROMSPOTTOFWDVOL");
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();


	// C variable
    double C_asOf;
    CCString C_dZcId;
    CCString C_fZcId;
    CCString C_dBSZcId;
    CCString C_fBSZcId;
    CCString C_volSwopBaseId;
    CCString C_volSwopForeignId;
    CCString C_fxVolId;
    double C_dMeanReversionBase;
    double C_dMeanReversionForeign;
    CCString C_dFxRdCorrId;
    CCString C_dFxRfCorrId;
    CCString C_dRdRfCorrId;
    double C_dCutOff;
    double C_dVolLongTerm;
    double C_calibBasisIncluded;
	VECTOR<double> C_ForwardVolDates;
    VECTOR<double> C_ForwardVolDates_Default;
    
	
	// error
	static int error;
	static char* reason = "";

    XL_readNumCell(XL_asOf,C_asOf," ARM_ERR: asOf Date: numeric expected",C_result);
	XL_readStrCell(XL_dZcId,C_dZcId," ARM_ERR: Domestic ZC: object expected",C_result);
    XL_readStrCell(XL_fZcId,C_fZcId," ARM_ERR: Foreign ZC: object expected",C_result);
    XL_readStrCell(XL_dBSZcId,C_dBSZcId," ARM_ERR: Domestic BS ZC: object expected",C_result);
    XL_readStrCell(XL_fBSZcId,C_fBSZcId," ARM_ERR: Foreign BS ZC: object expected",C_result);
    XL_readStrCell(XL_volSwopBaseId,C_volSwopBaseId," ARM_ERR: Domestic Swaption Vol: object expected",C_result);
    XL_readStrCell(XL_volSwopForeignId,C_volSwopForeignId," ARM_ERR: Foreign Swaption Vol: object expected",C_result);
    XL_readStrCell(XL_fxVolId,C_fxVolId," ARM_ERR: FX Vol: object expected",C_result);
    XL_readNumCell(XL_dMeanReversionBase,C_dMeanReversionBase," ARM_ERR: Domestic Mean Rev: numeric expected",C_result);
    XL_readNumCell(XL_dMeanReversionForeign,C_dMeanReversionForeign," ARM_ERR: Foreifn Mean Rev: numeric expected",C_result);
    XL_readStrCell(XL_dFxRdCorrId,C_dFxRdCorrId," ARM_ERR: Fx Rd Correlation : object expected",C_result);
    XL_readStrCell(XL_dFxRfCorrId,C_dFxRfCorrId," ARM_ERR: Fx Rf Correlation : object expected",C_result);
    XL_readStrCell(XL_dRdRfCorrId,C_dRdRfCorrId," ARM_ERR: Rd Rf Correlation : object expected",C_result);
    XL_readNumCell(XL_dCutOff,C_dCutOff," ARM_ERR: CutOff: numeric expected",C_result);
    XL_readNumCell(XL_dVolLongTerm,C_dVolLongTerm," ARM_ERR: Vol Long Term: numeric expected",C_result);
    XL_readNumCell(XL_calibBasisIncluded,C_calibBasisIncluded," ARM_ERR: Calibration with Basis ? : numeric expected",C_result);
    XL_readNumVectorWD(XL_ForwardVolDates,C_ForwardVolDates,C_ForwardVolDates_Default," ARM_ERR: Forward Vol dates: array of numeric expected",C_result);
   
	long retCode;
	long objId;
	CCString prevClass;
	
	CCString curClass = LOCAL_VOL_CURVE_LIN_CLASS;
	CCString stringId = GetLastCurCellEnvValue ();

	if(!stringId)
	{
		retCode = ARMLOCAL_ARM_CONV3FFROMSPOTTOFWDVOL(  C_asOf,
                                                        LocalGetNumObjectId(C_dZcId),
                                                        LocalGetNumObjectId(C_fZcId),
                                                        LocalGetNumObjectId(C_dBSZcId),
                                                        LocalGetNumObjectId(C_fBSZcId),
                                                        LocalGetNumObjectId(C_volSwopBaseId),
                                                        LocalGetNumObjectId(C_volSwopForeignId),
                                                        LocalGetNumObjectId(C_fxVolId),
                                                        C_dMeanReversionBase,
                                                        C_dMeanReversionForeign,
                                                        LocalGetNumObjectId(C_dFxRdCorrId),
                                                        LocalGetNumObjectId(C_dFxRfCorrId),
                                                        LocalGetNumObjectId(C_dRdRfCorrId),
                                                        C_dCutOff,
                                                        C_dVolLongTerm,
                                                        (long) C_calibBasisIncluded,
                                                        C_ForwardVolDates,
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
			retCode = ARMLOCAL_ARM_CONV3FFROMSPOTTOFWDVOL(  C_asOf,
                                                            LocalGetNumObjectId(C_dZcId),
                                                            LocalGetNumObjectId(C_fZcId),
                                                            LocalGetNumObjectId(C_dBSZcId),
                                                            LocalGetNumObjectId(C_fBSZcId),
                                                            LocalGetNumObjectId(C_volSwopBaseId),
                                                            LocalGetNumObjectId(C_volSwopForeignId),
                                                            LocalGetNumObjectId(C_fxVolId),
                                                            C_dMeanReversionBase,
                                                            C_dMeanReversionForeign,
                                                            LocalGetNumObjectId(C_dFxRdCorrId),
                                                            LocalGetNumObjectId(C_dFxRfCorrId),
                                                            LocalGetNumObjectId(C_dRdRfCorrId),
                                                            C_dCutOff,
                                                            C_dVolLongTerm,
                                                            C_calibBasisIncluded,
                                                            C_ForwardVolDates,
                                                            C_result,
                                                            objId);

			if (retCode == ARM_OK)
			{
				LocalSetCurCellEnvValue (curClass, objId); 

				stringId = LocalMakeObjectId (objId, curClass);
			}
		}
		else
		{
			FreeCurCellContent ();

			retCode = ARMLOCAL_ARM_CONV3FFROMSPOTTOFWDVOL(  C_asOf,
                                                            LocalGetNumObjectId(C_dZcId),
                                                            LocalGetNumObjectId(C_fZcId),
                                                            LocalGetNumObjectId(C_dBSZcId),
                                                            LocalGetNumObjectId(C_fBSZcId),
                                                            LocalGetNumObjectId(C_volSwopBaseId),
                                                            LocalGetNumObjectId(C_volSwopForeignId),
                                                            LocalGetNumObjectId(C_fxVolId),
                                                            C_dMeanReversionBase,
                                                            C_dMeanReversionForeign,
                                                            LocalGetNumObjectId(C_dFxRdCorrId),
                                                            LocalGetNumObjectId(C_dFxRfCorrId),
                                                            LocalGetNumObjectId(C_dRdRfCorrId),
                                                            C_dCutOff,
                                                            C_dVolLongTerm,
                                                            C_calibBasisIncluded,
                                                            C_ForwardVolDates,
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_ARM_CONV3FFROMSPOTTOFWDVOL" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}



__declspec(dllexport) LPXLOPER WINAPI Local_PXL_ARM_CONV3FFROMSPOTTOFWDVOL( LPXLOPER XL_asOf,
                                                                            LPXLOPER XL_dZcId,
                                                                            LPXLOPER XL_fZcId,
                                                                            LPXLOPER XL_dBSZcId,
                                                                            LPXLOPER XL_fBSZcId,
                                                                            LPXLOPER XL_volSwopBaseId,
                                                                            LPXLOPER XL_volSwopForeignId,
                                                                            LPXLOPER XL_fxVolId,
                                                                            LPXLOPER XL_dMeanReversionBase,
                                                                            LPXLOPER XL_dMeanReversionForeign,
                                                                            LPXLOPER XL_dFxRdCorrId,
                                                                            LPXLOPER XL_dFxRfCorrId,
                                                                            LPXLOPER XL_dRdRfCorrId,
                                                                            LPXLOPER XL_dCutOff,
                                                                            LPXLOPER XL_dVolLongTerm,
                                                                            LPXLOPER XL_calibBasisIncluded,
                                                                            LPXLOPER XL_ForwardVolDates)                                                                    
{
	ADD_LOG("Local_PXL_ARM_CONV3FFROMSPOTTOFWDVOL");
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();


	// C variable
    double C_asOf;
    CCString C_dZcId;
    CCString C_fZcId;
    CCString C_dBSZcId;
    CCString C_fBSZcId;
    CCString C_volSwopBaseId;
    CCString C_volSwopForeignId;
    CCString C_fxVolId;
    double C_dMeanReversionBase;
    double C_dMeanReversionForeign;
    CCString C_dFxRdCorrId;
    CCString C_dFxRfCorrId;
    CCString C_dRdRfCorrId;
    double C_dCutOff;
    double C_dVolLongTerm;
    double C_calibBasisIncluded;
	VECTOR<double> C_ForwardVolDates;
    VECTOR<double> C_ForwardVolDates_Default;
    
	
	// error
	static int error;
	static char* reason = "";

    XL_readNumCell(XL_asOf,C_asOf," ARM_ERR: asOf Date: numeric expected",C_result);
	XL_readStrCell(XL_dZcId,C_dZcId," ARM_ERR: Domestic ZC: object expected",C_result);
    XL_readStrCell(XL_fZcId,C_fZcId," ARM_ERR: Foreign ZC: object expected",C_result);
    XL_readStrCell(XL_dBSZcId,C_dBSZcId," ARM_ERR: Domestic BS ZC: object expected",C_result);
    XL_readStrCell(XL_fBSZcId,C_fBSZcId," ARM_ERR: Foreign BS ZC: object expected",C_result);
    XL_readStrCell(XL_volSwopBaseId,C_volSwopBaseId," ARM_ERR: Domestic Swaption Vol: object expected",C_result);
    XL_readStrCell(XL_volSwopForeignId,C_volSwopForeignId," ARM_ERR: Foreign Swaption Vol: object expected",C_result);
    XL_readStrCell(XL_fxVolId,C_fxVolId," ARM_ERR: FX Vol: object expected",C_result);
    XL_readNumCell(XL_dMeanReversionBase,C_dMeanReversionBase," ARM_ERR: Domestic Mean Rev: numeric expected",C_result);
    XL_readNumCell(XL_dMeanReversionForeign,C_dMeanReversionForeign," ARM_ERR: Foreifn Mean Rev: numeric expected",C_result);
    XL_readStrCell(XL_dFxRdCorrId,C_dFxRdCorrId," ARM_ERR: Fx Rd Correlation : object expected",C_result);
    XL_readStrCell(XL_dFxRfCorrId,C_dFxRfCorrId," ARM_ERR: Fx Rf Correlation : object expected",C_result);
    XL_readStrCell(XL_dRdRfCorrId,C_dRdRfCorrId," ARM_ERR: Rd Rf Correlation : object expected",C_result);
    XL_readNumCell(XL_dCutOff,C_dCutOff," ARM_ERR: CutOff: numeric expected",C_result);
    XL_readNumCell(XL_dVolLongTerm,C_dVolLongTerm," ARM_ERR: Vol Long Term: numeric expected",C_result);
    XL_readNumCell(XL_calibBasisIncluded,C_calibBasisIncluded," ARM_ERR: Calibration with Basis ? : numeric expected",C_result);
    XL_readNumVectorWD(XL_ForwardVolDates,C_ForwardVolDates,C_ForwardVolDates_Default," ARM_ERR: Forward Vol dates: array of numeric expected",C_result);
   
	long retCode;
	long objId;
	CCString prevClass;
	
	CCString curClass = LOCAL_VOL_CURVE_LIN_CLASS;
	CCString stringId = GetLastCurCellEnvValue ();

	retCode = ARMLOCAL_ARM_CONV3FFROMSPOTTOFWDVOL(  C_asOf,
                                                    LocalGetNumObjectId(C_dZcId),
                                                    LocalGetNumObjectId(C_fZcId),
                                                    LocalGetNumObjectId(C_dBSZcId),
                                                    LocalGetNumObjectId(C_fBSZcId),
                                                    LocalGetNumObjectId(C_volSwopBaseId),
                                                    LocalGetNumObjectId(C_volSwopForeignId),
                                                    LocalGetNumObjectId(C_fxVolId),
                                                    C_dMeanReversionBase,
                                                    C_dMeanReversionForeign,
                                                    LocalGetNumObjectId(C_dFxRdCorrId),
                                                    LocalGetNumObjectId(C_dFxRfCorrId),
                                                    LocalGetNumObjectId(C_dRdRfCorrId),
                                                    C_dCutOff,
                                                    C_dVolLongTerm,
                                                    (long) C_calibBasisIncluded,
                                                    C_ForwardVolDates,
													C_result);

	if(retCode == ARM_OK)
	{
		objId = C_result.getLong ();

		LocalSetCurCellEnvValue (curClass, objId); 

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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_PXL_ARM_CONV3FFROMSPOTTOFWDVOL" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


__declspec(dllexport) LPXLOPER WINAPI Local_ARM_CONV3FFROMFWDTOSPOTVOL_IN_PRCS( LPXLOPER XL_tree3f,
																	            LPXLOPER XL_PrcsObject)
                                                                      
{
	ADD_LOG("Local_ARM_CONV3FFROMFWDTOSPOTVOL_IN_PRCS");
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();


	// C variable
	CCString C_tree3f;
	CCString C_PrcsObject;
    
	// error
	static int error;
	static char* reason = "";

	XL_readStrCell(XL_tree3f,C_tree3f," ARM_ERR: tree 3f model: object expected",C_result);
    XL_readStrCell(XL_PrcsObject,C_PrcsObject," ARM_ERR: PRCS Object: object expected",C_result);
    
	long retCode;
	long objId;
	CCString prevClass;
	
	CCString curClass = LOCAL_VOL_CURVE_LIN_CLASS;
	CCString stringId = GetLastCurCellEnvValue ();

	if(!stringId)
	{
		retCode = ARMLOCAL_ARM_CONV3FFROMFWDTOSPOTVOL(  LocalGetNumObjectId(C_tree3f),
													    LocalGetNumObjectId(C_PrcsObject),
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
			retCode = ARMLOCAL_ARM_CONV3FFROMFWDTOSPOTVOL(  LocalGetNumObjectId(C_tree3f),
													        LocalGetNumObjectId(C_PrcsObject),
                                                            C_result,
                                                            objId);

			if (retCode == ARM_OK)
			{
				LocalSetCurCellEnvValue (curClass, objId); 

				stringId = LocalMakeObjectId (objId, curClass);
			}
		}
		else
		{
			FreeCurCellContent ();

			retCode = ARMLOCAL_ARM_CONV3FFROMFWDTOSPOTVOL(  LocalGetNumObjectId(C_tree3f),
													        LocalGetNumObjectId(C_PrcsObject),
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_ARM_CONV3FFROMFWDTOSPOTVOL_IN_PRCS" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


__declspec(dllexport) LPXLOPER WINAPI Local_PXL_ARM_CONV3FFROMFWDTOSPOTVOL_IN_PRCS (LPXLOPER XL_tree3f,
																	                LPXLOPER XL_PrcsObject)
                                                                            
{
	ADD_LOG("Local_PXL_ARM_CONV3FFROMFWDTOSPOTVOL_IN_PRCS ");
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();


	// C variable
	CCString C_tree3f;
    CCString C_PrcsObject;
    

	// error
	static int error;
	static char* reason = "";

	XL_readStrCell(XL_tree3f,C_tree3f," ARM_ERR: tree 3f model: object expected",C_result);
    XL_readStrCell(XL_PrcsObject,C_PrcsObject," ARM_ERR: PRCS Object: object expected",C_result);
    

	long retCode;
	long objId;
	CCString prevClass;
	
	CCString curClass = LOCAL_VOL_CURVE_LIN_CLASS;
	CCString stringId;

	retCode = ARMLOCAL_ARM_CONV3FFROMFWDTOSPOTVOL(  LocalGetNumObjectId(C_tree3f),
													LocalGetNumObjectId(C_PrcsObject),
                                                    C_result);

	if ( retCode == ARM_OK)
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_PXL_ARM_CONV3FFROMFWDTOSPOTVOL_IN_PRCS" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


__declspec(dllexport) LPXLOPER WINAPI Local_ARM_CONV3FFROMFWDTOSPOTVOL( LPXLOPER XL_asOf,
                                                                        LPXLOPER XL_dZcId,
                                                                        LPXLOPER XL_fZcId,
                                                                        LPXLOPER XL_dBSZcId,
                                                                        LPXLOPER XL_fBSZcId,
                                                                        LPXLOPER XL_volSwopBaseId,
                                                                        LPXLOPER XL_volSwopForeignId,
                                                                        LPXLOPER XL_fxVolId,
                                                                        LPXLOPER XL_dMeanReversionBase,
                                                                        LPXLOPER XL_dMeanReversionForeign,
                                                                        LPXLOPER XL_dFxRdCorrId,
                                                                        LPXLOPER XL_dFxRfCorrId,
                                                                        LPXLOPER XL_dRdRfCorrId,
                                                                        LPXLOPER XL_dCutOff,
                                                                        LPXLOPER XL_dVolLongTerm,
                                                                        LPXLOPER XL_calibBasisIncluded)                                                                    
{
	ADD_LOG("Local_ARM_CONV3FFROMFWDTOSPOTVOL");
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();


	// C variable
    double C_asOf;
    CCString C_dZcId;
    CCString C_fZcId;
    CCString C_dBSZcId;
    CCString C_fBSZcId;
    CCString C_volSwopBaseId;
    CCString C_volSwopForeignId;
    CCString C_fxVolId;
    double C_dMeanReversionBase;
    double C_dMeanReversionForeign;
    CCString C_dFxRdCorrId;
    CCString C_dFxRfCorrId;
    CCString C_dRdRfCorrId;
    double C_dCutOff;
    double C_dVolLongTerm;
    double C_calibBasisIncluded;
	
	
	// error
	static int error;
	static char* reason = "";

    XL_readNumCell(XL_asOf,C_asOf," ARM_ERR: asOf Date: numeric expected",C_result);
	XL_readStrCell(XL_dZcId,C_dZcId," ARM_ERR: Domestic ZC: object expected",C_result);
    XL_readStrCell(XL_fZcId,C_fZcId," ARM_ERR: Foreign ZC: object expected",C_result);
    XL_readStrCell(XL_dBSZcId,C_dBSZcId," ARM_ERR: Domestic BS ZC: object expected",C_result);
    XL_readStrCell(XL_fBSZcId,C_fBSZcId," ARM_ERR: Foreign BS ZC: object expected",C_result);
    XL_readStrCell(XL_volSwopBaseId,C_volSwopBaseId," ARM_ERR: Domestic Swaption Vol: object expected",C_result);
    XL_readStrCell(XL_volSwopForeignId,C_volSwopForeignId," ARM_ERR: Foreign Swaption Vol: object expected",C_result);
    XL_readStrCell(XL_fxVolId,C_fxVolId," ARM_ERR: FX Vol: object expected",C_result);
    XL_readNumCell(XL_dMeanReversionBase,C_dMeanReversionBase," ARM_ERR: Domestic Mean Rev: numeric expected",C_result);
    XL_readNumCell(XL_dMeanReversionForeign,C_dMeanReversionForeign," ARM_ERR: Foreifn Mean Rev: numeric expected",C_result);
    XL_readStrCell(XL_dFxRdCorrId,C_dFxRdCorrId," ARM_ERR: Fx Rd Correlation : object expected",C_result);
    XL_readStrCell(XL_dFxRfCorrId,C_dFxRfCorrId," ARM_ERR: Fx Rf Correlation : object expected",C_result);
    XL_readStrCell(XL_dRdRfCorrId,C_dRdRfCorrId," ARM_ERR: Rd Rf Correlation : object expected",C_result);
    XL_readNumCell(XL_dCutOff,C_dCutOff," ARM_ERR: CutOff: numeric expected",C_result);
    XL_readNumCell(XL_dVolLongTerm,C_dVolLongTerm," ARM_ERR: Vol Long Term: numeric expected",C_result);
    XL_readNumCell(XL_calibBasisIncluded,C_calibBasisIncluded," ARM_ERR: Calibration with Basis ? : numeric expected",C_result);
    
	long retCode;
	long objId;
	CCString prevClass;
	
	CCString curClass = LOCAL_VOL_CURVE_LIN_CLASS;
	CCString stringId = GetLastCurCellEnvValue ();

	if(!stringId)
	{
		retCode = ARMLOCAL_ARM_CONV3FFROMFWDTOSPOTVOL(  C_asOf,
                                                        LocalGetNumObjectId(C_dZcId),
                                                        LocalGetNumObjectId(C_fZcId),
                                                        LocalGetNumObjectId(C_dBSZcId),
                                                        LocalGetNumObjectId(C_fBSZcId),
                                                        LocalGetNumObjectId(C_volSwopBaseId),
                                                        LocalGetNumObjectId(C_volSwopForeignId),
                                                        LocalGetNumObjectId(C_fxVolId),
                                                        C_dMeanReversionBase,
                                                        C_dMeanReversionForeign,
                                                        LocalGetNumObjectId(C_dFxRdCorrId),
                                                        LocalGetNumObjectId(C_dFxRfCorrId),
                                                        LocalGetNumObjectId(C_dRdRfCorrId),
                                                        C_dCutOff,
                                                        C_dVolLongTerm,
                                                        (long) C_calibBasisIncluded,
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
			retCode = ARMLOCAL_ARM_CONV3FFROMFWDTOSPOTVOL(  C_asOf,
                                                            LocalGetNumObjectId(C_dZcId),
                                                            LocalGetNumObjectId(C_fZcId),
                                                            LocalGetNumObjectId(C_dBSZcId),
                                                            LocalGetNumObjectId(C_fBSZcId),
                                                            LocalGetNumObjectId(C_volSwopBaseId),
                                                            LocalGetNumObjectId(C_volSwopForeignId),
                                                            LocalGetNumObjectId(C_fxVolId),
                                                            C_dMeanReversionBase,
                                                            C_dMeanReversionForeign,
                                                            LocalGetNumObjectId(C_dFxRdCorrId),
                                                            LocalGetNumObjectId(C_dFxRfCorrId),
                                                            LocalGetNumObjectId(C_dRdRfCorrId),
                                                            C_dCutOff,
                                                            C_dVolLongTerm,
                                                            C_calibBasisIncluded,
                                                            C_result,
                                                            objId);

			if (retCode == ARM_OK)
			{
				LocalSetCurCellEnvValue (curClass, objId); 

				stringId = LocalMakeObjectId (objId, curClass);
			}
		}
		else
		{
			FreeCurCellContent ();

			retCode = ARMLOCAL_ARM_CONV3FFROMFWDTOSPOTVOL(  C_asOf,
                                                            LocalGetNumObjectId(C_dZcId),
                                                            LocalGetNumObjectId(C_fZcId),
                                                            LocalGetNumObjectId(C_dBSZcId),
                                                            LocalGetNumObjectId(C_fBSZcId),
                                                            LocalGetNumObjectId(C_volSwopBaseId),
                                                            LocalGetNumObjectId(C_volSwopForeignId),
                                                            LocalGetNumObjectId(C_fxVolId),
                                                            C_dMeanReversionBase,
                                                            C_dMeanReversionForeign,
                                                            LocalGetNumObjectId(C_dFxRdCorrId),
                                                            LocalGetNumObjectId(C_dFxRfCorrId),
                                                            LocalGetNumObjectId(C_dRdRfCorrId),
                                                            C_dCutOff,
                                                            C_dVolLongTerm,
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_ARM_CONV3FFROMFWDTOSPOTVOL" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


__declspec(dllexport) LPXLOPER WINAPI Local_PXL_ARM_CONV3FFROMFWDTOSPOTVOL( LPXLOPER XL_asOf,
                                                                            LPXLOPER XL_dZcId,
                                                                            LPXLOPER XL_fZcId,
                                                                            LPXLOPER XL_dBSZcId,
                                                                            LPXLOPER XL_fBSZcId,
                                                                            LPXLOPER XL_volSwopBaseId,
                                                                            LPXLOPER XL_volSwopForeignId,
                                                                            LPXLOPER XL_fxVolId,
                                                                            LPXLOPER XL_dMeanReversionBase,
                                                                            LPXLOPER XL_dMeanReversionForeign,
                                                                            LPXLOPER XL_dFxRdCorrId,
                                                                            LPXLOPER XL_dFxRfCorrId,
                                                                            LPXLOPER XL_dRdRfCorrId,
                                                                            LPXLOPER XL_dCutOff,
                                                                            LPXLOPER XL_dVolLongTerm,
                                                                            LPXLOPER XL_calibBasisIncluded)
{
	ADD_LOG("Local_PXL_ARM_CONV3FFROMFWDTOSPOTVOL");
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();


	// C variable
    double C_asOf;
    CCString C_dZcId;
    CCString C_fZcId;
    CCString C_dBSZcId;
    CCString C_fBSZcId;
    CCString C_volSwopBaseId;
    CCString C_volSwopForeignId;
    CCString C_fxVolId;
    double C_dMeanReversionBase;
    double C_dMeanReversionForeign;
    CCString C_dFxRdCorrId;
    CCString C_dFxRfCorrId;
    CCString C_dRdRfCorrId;
    double C_dCutOff;
    double C_dVolLongTerm;
    double C_calibBasisIncluded;

	// error
	static int error;
	static char* reason = "";

    XL_readNumCell(XL_asOf,C_asOf," ARM_ERR: asOf Date: numeric expected",C_result);
	XL_readStrCell(XL_dZcId,C_dZcId," ARM_ERR: Domestic ZC: object expected",C_result);
    XL_readStrCell(XL_fZcId,C_fZcId," ARM_ERR: Foreign ZC: object expected",C_result);
    XL_readStrCell(XL_dBSZcId,C_dBSZcId," ARM_ERR: Domestic BS ZC: object expected",C_result);
    XL_readStrCell(XL_fBSZcId,C_fBSZcId," ARM_ERR: Foreign BS ZC: object expected",C_result);
    XL_readStrCell(XL_volSwopBaseId,C_volSwopBaseId," ARM_ERR: Domestic Swaption Vol: object expected",C_result);
    XL_readStrCell(XL_volSwopForeignId,C_volSwopForeignId," ARM_ERR: Foreign Swaption Vol: object expected",C_result);
    XL_readStrCell(XL_fxVolId,C_fxVolId," ARM_ERR: FX Vol: object expected",C_result);
    XL_readNumCell(XL_dMeanReversionBase,C_dMeanReversionBase," ARM_ERR: Domestic Mean Rev: numeric expected",C_result);
    XL_readNumCell(XL_dMeanReversionForeign,C_dMeanReversionForeign," ARM_ERR: Foreifn Mean Rev: numeric expected",C_result);
    XL_readStrCell(XL_dFxRdCorrId,C_dFxRdCorrId," ARM_ERR: Fx Rd Correlation : object expected",C_result);
    XL_readStrCell(XL_dFxRfCorrId,C_dFxRfCorrId," ARM_ERR: Fx Rf Correlation : object expected",C_result);
    XL_readStrCell(XL_dRdRfCorrId,C_dRdRfCorrId," ARM_ERR: Rd Rf Correlation : object expected",C_result);
    XL_readNumCell(XL_dCutOff,C_dCutOff," ARM_ERR: CutOff: numeric expected",C_result);
    XL_readNumCell(XL_dVolLongTerm,C_dVolLongTerm," ARM_ERR: Vol Long Term: numeric expected",C_result);
    XL_readNumCell(XL_calibBasisIncluded,C_calibBasisIncluded," ARM_ERR: Calibration with Basis ? : numeric expected",C_result);
    
	long retCode;
	long objId;
	CCString prevClass;
	
	CCString curClass = LOCAL_VOL_CURVE_LIN_CLASS;
	CCString stringId = GetLastCurCellEnvValue ();

	retCode = ARMLOCAL_ARM_CONV3FFROMFWDTOSPOTVOL(  C_asOf,
                                                    LocalGetNumObjectId(C_dZcId),
                                                    LocalGetNumObjectId(C_fZcId),
                                                    LocalGetNumObjectId(C_dBSZcId),
                                                    LocalGetNumObjectId(C_fBSZcId),
                                                    LocalGetNumObjectId(C_volSwopBaseId),
                                                    LocalGetNumObjectId(C_volSwopForeignId),
                                                    LocalGetNumObjectId(C_fxVolId),
                                                    C_dMeanReversionBase,
                                                    C_dMeanReversionForeign,
                                                    LocalGetNumObjectId(C_dFxRdCorrId),
                                                    LocalGetNumObjectId(C_dFxRfCorrId),
                                                    LocalGetNumObjectId(C_dRdRfCorrId),
                                                    C_dCutOff,
                                                    C_dVolLongTerm,
                                                    (long) C_calibBasisIncluded,
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_PXL_ARM_CONV3FFROMFWDTOSPOTVOL" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}

__declspec(dllexport) LPXLOPER WINAPI Local_SetVolCurveName(LPXLOPER XL_curve,
															LPXLOPER XL_name)
{
	ADD_LOG("Local_SetVolCurveName");
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
		CCString C_name;
		// error
		static int error;
		static char* reason = "";

		XL_readStrCell(XL_curve,C_curve," ARM_ERR: string expected",C_result);
		XL_readStrCell(XL_name,C_name," ARM_ERR: string expected",C_result);

		long retCode = ARMLOCAL_SetVolCurveName(LocalGetNumObjectId (C_curve), C_name, C_result);

		if(retCode == ARM_OK)
		{
			ExcelTools::convert(CCSTringToSTLString(C_curve),&XL_result); 
			FreeCurCellErr ();
			// .xltype = xltypeNum;
			// XL_result.val.num = C_result.getDouble ();
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_SetVolCurveName" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


// Hypercube
__declspec(dllexport) LPXLOPER WINAPI Local_HyperCube(LPXLOPER XL_volCubes,
													  LPXLOPER XL_keys)
{
	ADD_LOG("Local_HyperCube");
	static XLOPER	XL_result;
	ARM_result		C_result;

	ARM_XL_TRY_BLOCK_BEGIN
	{
	    ARM_NOCALCIFWIZ();
	    
		vector<CCString>	C_keys;
	    vector<CCString>	C_volCubes;
	    vector<long>		C_volCubeIds;

	    static int		error;
	    static char*	reason = "";

	    XL_readStrVector(XL_volCubes, C_volCubes, " ARM_ERR: vol cubes : array of string expected", XL_TYPE_STRING, C_result);
		XL_readStrVector(XL_keys, C_keys, " ARM_ERR: keys : array of string expected", XL_TYPE_STRING, C_result);
	    
		
		if( C_volCubes.size() != C_keys.size() )
		{
			C_result.setMsg("Array of keys and array of volCubes should have the same size");
			ARM_ERR();
		}
		else
		{
			for(int i=0; i < C_volCubes.size(); i++)
				C_volCubeIds.push_back(LocalGetNumObjectId(C_volCubes[i]));

			long	retCode;
			long	objId;

			CCString	prevClass;
			CCString	curClass = LOCAL_HYPER_CUBE_CLASS;
			CCString	stringId = GetLastCurCellEnvValue ();
			
			if(!stringId)
			{
				retCode = ARMLOCAL_HyperCube(C_volCubeIds, C_keys, C_result);

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
					retCode = ARMLOCAL_HyperCube(C_volCubeIds, C_keys, C_result, objId);

					if ( retCode == ARM_OK )
					{
						LocalSetCurCellEnvValue (curClass, objId); 

						stringId = LocalMakeObjectId (objId, curClass);
					}
				}
				else
				{
					FreeCurCellContent ();
					retCode = ARMLOCAL_HyperCube(C_volCubeIds, C_keys, C_result);
					
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

	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_HyperCube" )

	return	(LPXLOPER)&XL_result;
}


__declspec(dllexport) LPXLOPER WINAPI Local_PXL_HyperCube(LPXLOPER XL_volCubes,
														  LPXLOPER XL_keys)
{
	ADD_LOG("Local_PXL_HyperCube");
	static XLOPER	XL_result;
	ARM_result		C_result;

	ARM_XL_TRY_BLOCK_BEGIN
	{
	    ARM_NOCALCIFWIZ();
	    
		vector<CCString>	C_keys;
	    vector<CCString>	C_volCubes;
	    vector<long>		C_volCubeIds;

	    static int		error;
	    static char*	reason = "";

	    XL_readStrVector(XL_volCubes, C_volCubes, " ARM_ERR: vol cubes : array of string expected", XL_TYPE_STRING, C_result);
		XL_readStrVector(XL_keys, C_keys, " ARM_ERR: keys : array of string expected", XL_TYPE_STRING, C_result);
	    
		
		if( C_volCubes.size() != C_keys.size() )
		{
			C_result.setMsg("Array of keys and array of volCubes should have the same size");
			ARM_ERR();
		}
		else
		{
			for(int i=0; i < C_volCubes.size(); i++)
				C_volCubeIds.push_back(LocalGetNumObjectId(C_volCubes[i]));

			long	retCode;
			long	objId;

			CCString	curClass = LOCAL_HYPER_CUBE_CLASS;
			CCString	stringId = GetLastCurCellEnvValue ();
			
			retCode = ARMLOCAL_HyperCube(C_volCubeIds, C_keys, C_result);

			if( retCode == ARM_OK )
			{
				objId = C_result.getLong();

				stringId = LocalMakeObjectId(objId, curClass);
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

	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_HyperCube" )

	return	(LPXLOPER)&XL_result;
}


__declspec(dllexport) LPXLOPER WINAPI Local_CreateCorrelCubeByExpiry(LPXLOPER XL_hyperCube,
																	 LPXLOPER XL_TenorList,
																	 LPXLOPER XL_ExpiryList,
																	 LPXLOPER XL_IntersurfaceInterpol)
{
	ADD_LOG("Local_CreateCorrelCubeByExpiry");
	static XLOPER	XL_result;
	ARM_result		C_result;

	ARM_XL_TRY_BLOCK_BEGIN
	{
	    ARM_NOCALCIFWIZ();

		CCString			C_hyperCube;
		vector<CCString>	C_tenorList;
		vector<CCString>	C_expiryList;
		CCString			C_intersurfaceInterpol;

	    static int		error;
	    static char*	reason = "";

	    XL_readStrCell(XL_hyperCube, C_hyperCube, " ARM_ERR: HyperCube : string expected", C_result);
		XL_readStrVectorWD(XL_TenorList, C_tenorList, C_tenorList, 
						   " ARM_ERR: tenor list : array of string expected", XL_TYPE_STRING, C_result);
		XL_readStrVectorWD(XL_ExpiryList, C_expiryList, C_expiryList,
						 " ARM_ERR: expiry list : array of string expected", XL_TYPE_STRING, C_result);
		XL_readStrCellWD(XL_IntersurfaceInterpol, C_intersurfaceInterpol, "NO",
						 " ARM_ERR: Intersurface Interpol : YES/NO expectded", C_result);
	    		
		long	retCode;
		long	objId;

		CCString	prevClass;
		CCString	curClass = LOCAL_VOL_CUBE_CLASS;
		CCString	stringId = GetLastCurCellEnvValue ();
		
		if(!stringId)
		{
			retCode = ARMLOCAL_CreateCorrelCubeByExpiry(LocalGetNumObjectId(C_hyperCube), 
														C_tenorList, 
														C_expiryList,
														C_intersurfaceInterpol,
														C_result);

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
				retCode = ARMLOCAL_CreateCorrelCubeByExpiry(LocalGetNumObjectId(C_hyperCube), 
															C_tenorList, 
															C_expiryList,
															C_intersurfaceInterpol,
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
				retCode = ARMLOCAL_CreateCorrelCubeByExpiry(LocalGetNumObjectId(C_hyperCube), 
															C_tenorList, 
															C_expiryList, 
															C_intersurfaceInterpol,
															C_result);				
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

	ARM_XL_TRY_BLOCK_END

	ARM_XL_CATCH_ARM_EXPT

	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_CreateCorrelCubeByExpiry" )

	return	(LPXLOPER)&XL_result;
}

__declspec(dllexport) LPXLOPER WINAPI Local_PXL_CreateCorrelCubeByExpiry(LPXLOPER XL_hyperCube,
																	     LPXLOPER XL_TenorList,
																	     LPXLOPER XL_ExpiryList,
																	     LPXLOPER XL_IntersurfaceInterpol)
{
	ADD_LOG("Local_PXL_CreateCorrelCubeByExpiry");
	static XLOPER	XL_result;
	ARM_result		C_result;

	ARM_XL_TRY_BLOCK_BEGIN
	{
	    ARM_NOCALCIFWIZ();

		CCString			C_hyperCube;
		vector<CCString>	C_tenorList;
		vector<CCString>	C_expiryList;
		CCString			C_intersurfaceInterpol;

	    static int		error;
	    static char*	reason = "";

	    XL_readStrCell(XL_hyperCube, C_hyperCube, " ARM_ERR: HyperCube : string expected", C_result);
		XL_readStrVectorWD(XL_TenorList, C_tenorList, C_tenorList, 
						   " ARM_ERR: tenor list : array of string expected", XL_TYPE_STRING, C_result);
		XL_readStrVectorWD(XL_ExpiryList, C_expiryList, C_expiryList,
						 " ARM_ERR: expiry list : array of string expected", XL_TYPE_STRING, C_result);
		XL_readStrCellWD(XL_IntersurfaceInterpol, C_intersurfaceInterpol, "NO",
						 " ARM_ERR: Intersurface Interpol : YES/NO expectded", C_result);
	    		
		long	retCode;
		long	objId;

		CCString	curClass = LOCAL_VOL_CUBE_CLASS;
		CCString stringId;

		retCode = ARMLOCAL_CreateCorrelCubeByExpiry(LocalGetNumObjectId(C_hyperCube), 
														C_tenorList, 
														C_expiryList,
														C_intersurfaceInterpol,
														C_result);

		if( retCode == ARM_OK )
		{
			objId = C_result.getLong();

			stringId = LocalMakeObjectId(objId, curClass);
	
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

	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_PXL_CreateCorrelCubeByExpiry" )

	return	(LPXLOPER)&XL_result;
}






__declspec(dllexport) LPXLOPER WINAPI Local_ComputeCorrelByExpiry(LPXLOPER XL_correlCube,
																  LPXLOPER XL_Expiry,
																  LPXLOPER XL_Tenor1,
																  LPXLOPER XL_Tenor2)
{
	ADD_LOG("Local_ComputeCorrelByExpiry");
	static XLOPER	XL_result;
	ARM_result		C_result;

	ARM_XL_TRY_BLOCK_BEGIN
	{
	    ARM_NOCALCIFWIZ();
	    
		CCString	C_correlCube;
		double		C_Expiry;
		double		C_Tenor1;
		double		C_Tenor2;

	    static int		error;
	    static char*	reason = "";

	    XL_readStrCell(XL_correlCube, C_correlCube, " ARM_ERR: correl cube : string expected", C_result);
		XL_readNumCell(XL_Expiry, C_Expiry, " ARM_ERR: Expiry : double expected", C_result);
		XL_readNumCell(XL_Tenor1, C_Tenor1, " ARM_ERR: Tenor1 : double expected", C_result);
		XL_readNumCell(XL_Tenor2, C_Tenor2, " ARM_ERR: Tenor2 : double expected", C_result);		
	    
		long retCode = ARMLOCAL_ComputeCorrelByExpiry(LocalGetNumObjectId(C_correlCube), 
													  C_Expiry, C_Tenor1, C_Tenor2,
													  C_result);
		if( retCode == ARM_OK )
		{
			FreeCurCellErr();

			XL_result.xltype = xltypeNum;
			XL_result.val.num = C_result.getDouble ();
		}
		else
		{
			ARM_ERR();
		}
	}

	ARM_XL_TRY_BLOCK_END

	ARM_XL_CATCH_ARM_EXPT

	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_ComputeCorrelByExpiry" )

	return	(LPXLOPER)&XL_result;
}

__declspec(dllexport) LPXLOPER WINAPI Local_ComputeHyperCorrel(LPXLOPER XL_correlCube,
														       LPXLOPER XL_Tenor1,
															   LPXLOPER XL_Tenor2,
															   LPXLOPER XL_Expiry,
															   LPXLOPER XL_Moneyness)
{
	ADD_LOG("Local_ComputeHyperCorrel");
	static XLOPER	XL_result;
	ARM_result		C_result;

	ARM_XL_TRY_BLOCK_BEGIN
	{
	    ARM_NOCALCIFWIZ();
	    
		CCString	C_correlCube;
		double		C_Expiry;
		double		C_Tenor1;
		double		C_Tenor2;
		double		C_Moneyness;
		double		C_defaultMoneyness = 0.0;

	    static int		error;
	    static char*	reason = "";

	    XL_readStrCell(XL_correlCube, C_correlCube, " ARM_ERR: correl cube : string expected", C_result);
		XL_readNumCell(XL_Expiry, C_Expiry, " ARM_ERR: Expiry : double expected", C_result);
		XL_readNumCell(XL_Tenor1, C_Tenor1, " ARM_ERR: Tenor1 : double expected", C_result);
		XL_readNumCell(XL_Tenor2, C_Tenor2, " ARM_ERR: Tenor2 : double expected", C_result);
		XL_readNumCellWD(XL_Moneyness, C_Moneyness, C_defaultMoneyness, " ARM_ERR: Moneyness : double expected", C_result);
		
		long retCode = ARMLOCAL_ComputeHyperCorrel(LocalGetNumObjectId(C_correlCube), 
													   C_Tenor1, C_Tenor2,C_Expiry,
													   C_Moneyness, C_result);
		if( retCode == ARM_OK )
		{
			FreeCurCellErr();

			XL_result.xltype = xltypeNum;
			XL_result.val.num = C_result.getDouble ();
		}
		else
		{
			ARM_ERR();
		}
	}

	ARM_XL_TRY_BLOCK_END

	ARM_XL_CATCH_ARM_EXPT

	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_ComputeHyperCorrel" )

	return	(LPXLOPER)&XL_result;
}
	  


// IndexIndexCorrelCube
__declspec(dllexport) LPXLOPER WINAPI Local_IndexIndexCorrelCube(LPXLOPER XL_correlCurves,
																 LPXLOPER XL_Tenors1List,
																 LPXLOPER XL_Tenors2List,
																 LPXLOPER XL_IntersurfaceInterpol)
{
	ADD_LOG("Local_IndexIndexCorrelCube");
	static XLOPER	XL_result;
	ARM_result		C_result;

	ARM_XL_TRY_BLOCK_BEGIN
	{
	    ARM_NOCALCIFWIZ();
	    
		vector<CCString>	C_tenors1List;
		vector<CCString>	C_tenors2List;
	    vector<CCString>	C_correlCurves;
	    vector<long>		C_correlCurveIds;
		CCString			C_intersurfaceInterpol;

	    static int		error;
	    static char*	reason = "";

	    XL_readStrVector(XL_correlCurves, C_correlCurves, " ARM_ERR: correl curves : array of string expected", XL_TYPE_STRING, C_result);
		XL_readStrVector(XL_Tenors1List, C_tenors1List, " ARM_ERR: tenors1 list : array of string expected", XL_TYPE_STRING, C_result);
		XL_readStrVector(XL_Tenors2List, C_tenors2List, " ARM_ERR: tenors2 list : array of string expected", XL_TYPE_STRING, C_result);
		XL_readStrCellWD(XL_IntersurfaceInterpol, C_intersurfaceInterpol, "YES",
						 " ARM_ERR: Intersurface Interpol : YES/NO expectded", C_result);
	    
		
		if( C_correlCurves.size() != C_tenors1List.size() )
		{
			C_result.setMsg("Array of tenors1 and array of correlCurves should have the same size");
			ARM_ERR();
		}
		else if( C_tenors2List.size() != C_tenors1List.size() )
		{
			C_result.setMsg("Array of tenors2 and array of correlCurves should have the same size");
			ARM_ERR();
		}
		else
		{
			for(int i=0; i < C_correlCurves.size(); i++)
				C_correlCurveIds.push_back(LocalGetNumObjectId(C_correlCurves[i]));

			long	retCode;
			long	objId;

			CCString	prevClass;
			CCString	curClass = LOCAL_INDEX_INDEX_CORREL_CUBE_CLASS;
			CCString	stringId = GetLastCurCellEnvValue ();
			
			if(!stringId)
			{
				retCode = ARMLOCAL_IndexIndexCorrelCube(C_correlCurveIds, 
														C_tenors1List, C_tenors2List, 
														C_intersurfaceInterpol, C_result);

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
					retCode = ARMLOCAL_IndexIndexCorrelCube(C_correlCurveIds, 
															C_tenors1List, C_tenors2List, 
															C_intersurfaceInterpol, C_result, objId);

					if ( retCode == ARM_OK )
					{
						LocalSetCurCellEnvValue (curClass, objId); 

						stringId = LocalMakeObjectId (objId, curClass);
					}
				}
				else
				{
					FreeCurCellContent ();
					retCode = ARMLOCAL_IndexIndexCorrelCube(C_correlCurveIds, 
															C_tenors1List, C_tenors2List, 
															C_intersurfaceInterpol, C_result);
					
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

	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_IndexIndexCorrelCube" )

	return	(LPXLOPER)&XL_result;
}


__declspec(dllexport) LPXLOPER WINAPI Local_PXL_IndexIndexCorrelCube(LPXLOPER XL_correlCurves,
																	 LPXLOPER XL_Tenors1List,
																	 LPXLOPER XL_Tenors2List,
																	 LPXLOPER XL_IntersurfaceInterpol)
{
	ADD_LOG("Local_PXL_IndexIndexCorrelCube");
	static XLOPER	XL_result;
	ARM_result		C_result;

	ARM_XL_TRY_BLOCK_BEGIN
	{
	    ARM_NOCALCIFWIZ();
	    
		vector<CCString>	C_tenors1List;
		vector<CCString>	C_tenors2List;
	    vector<CCString>	C_correlCurves;
	    vector<long>		C_correlCurveIds;
		CCString			C_intersurfaceInterpol;

	    static int		error;
	    static char*	reason = "";

	    XL_readStrVector(XL_correlCurves, C_correlCurves, " ARM_ERR: correl curves : array of string expected", XL_TYPE_STRING, C_result);
		XL_readStrVector(XL_Tenors1List, C_tenors1List, " ARM_ERR: tenors1 list : array of string expected", XL_TYPE_STRING, C_result);
		XL_readStrVector(XL_Tenors2List, C_tenors2List, " ARM_ERR: tenors2 list : array of string expected", XL_TYPE_STRING, C_result);
		XL_readStrCellWD(XL_IntersurfaceInterpol, C_intersurfaceInterpol, "YES",
						 " ARM_ERR: Intersurface Interpol : YES/NO expectded", C_result);
	    
		
		if( C_correlCurves.size() != C_tenors1List.size() )
		{
			C_result.setMsg("Array of tenors1 and array of correlCurves should have the same size");
			ARM_ERR();
		}
		else if( C_tenors2List.size() != C_tenors1List.size() )
		{
			C_result.setMsg("Array of tenors2 and array of correlCurves should have the same size");
			ARM_ERR();
		}
		else
		{
			for(int i=0; i < C_correlCurves.size(); i++)
				C_correlCurveIds.push_back(LocalGetNumObjectId(C_correlCurves[i]));

			long	retCode;
			long	objId;

			CCString	curClass = LOCAL_INDEX_INDEX_CORREL_CUBE_CLASS;
			CCString	stringId = GetLastCurCellEnvValue ();
			
			retCode = ARMLOCAL_IndexIndexCorrelCube(C_correlCurveIds, 
													C_tenors1List, C_tenors2List, 
													C_intersurfaceInterpol, C_result);

			if( retCode == ARM_OK )
			{
				objId = C_result.getLong();

				stringId = LocalMakeObjectId(objId, curClass);
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

	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_IndexIndexCorrelCube" )

	return	(LPXLOPER)&XL_result;
}


__declspec(dllexport) LPXLOPER WINAPI Local_ComputeIndexIndexCorrel(LPXLOPER XL_correlCube,
																	LPXLOPER XL_Tenor1,
																	LPXLOPER XL_Tenor2,
																	LPXLOPER XL_Expiry1,
																	LPXLOPER XL_Expiry2)
{
	ADD_LOG("Local_ComputeIndexIndexCorrel");
	static XLOPER	XL_result;
	ARM_result		C_result;

	ARM_XL_TRY_BLOCK_BEGIN
	{
	    ARM_NOCALCIFWIZ();
	    
		CCString	C_correlCube;
		double		C_Tenor1;
		double		C_Tenor2;
		double		C_Expiry1;
	    double		C_Expiry2;

	    static int		error;
	    static char*	reason = "";

	    XL_readStrCell(XL_correlCube, C_correlCube, " ARM_ERR: correl cube : string expected", C_result);
		XL_readNumCell(XL_Tenor1, C_Tenor1, " ARM_ERR: Tenor1 : double expected", C_result);
		XL_readNumCell(XL_Tenor2, C_Tenor2, " ARM_ERR: Tenor2 : double expected", C_result);
		XL_readNumCell(XL_Expiry1, C_Expiry1, " ARM_ERR: Expiry1 : double expected", C_result);
		XL_readNumCell(XL_Expiry2, C_Expiry2, " ARM_ERR: Expiry2 : double expected", C_result);
	    

		long retCode = ARMLOCAL_ComputeIndexIndexCorrel(LocalGetNumObjectId(C_correlCube), 
														C_Tenor1, C_Tenor2, 
														C_Expiry1, C_Expiry2,
														C_result);
		if( retCode == ARM_OK )
		{
			FreeCurCellErr();

			XL_result.xltype = xltypeNum;
			XL_result.val.num = C_result.getDouble ();
		}
		else
		{
			ARM_ERR();
		}
	}

	ARM_XL_TRY_BLOCK_END

	ARM_XL_CATCH_ARM_EXPT

	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_ComputeIndexIndexCorrel" )

	return	(LPXLOPER)&XL_result;
}

__declspec(dllexport) LPXLOPER WINAPI Local_BumpVolatilityCorrelManager(LPXLOPER XL_CorrelManager,
																		LPXLOPER XL_TypeCorr,
																		LPXLOPER XL_Value,
																		LPXLOPER XL_nthLine,
																		LPXLOPER XL_nthCol,
																		LPXLOPER XL_isCumul,
																		LPXLOPER XL_isAbsolute,
																		LPXLOPER XL_Currency)
{
	ADD_LOG("Local_BumpVolatilityCorrelManager");
	static XLOPER XL_result;
	ARM_result C_result;

	ARM_XL_TRY_BLOCK_BEGIN
	{
	    ARM_NOCALCIFWIZ();
	
		CCString C_CorrelManager;
		
		double C_Value;

		double C_nthLine;
	    double C_nthLine_default = 0.;

	    double C_nthCol;
	    double C_nthCol_default = 0.;

	    CCString C_isCumul;
	    long cumulId;

	    CCString C_isAbsolute;
	    long absoluteId;

		CCString C_isToClone;
	    long isToCloneId = 1;

		CCString C_Currency;
		vector<CCString> C_TypeCorr;
		long TypeCorrId = 0;
		vector<string> C_mktTag;
		vector<string> C_intraMktTag;

		static int error;
	    static char* reason = "";

	    XL_readStrCell(XL_CorrelManager,C_CorrelManager," ARM_ERR: CorrelManager: object expected",C_result);
		XL_readStrVector(XL_TypeCorr,C_TypeCorr," ARM_ERR: correl type: string expected",DOUBLE_TYPE,C_result);
	    XL_readNumCell(XL_Value,C_Value," ARM_ERR: value to bump: numeric expected",C_result);
	    XL_readNumCellWD(XL_nthLine,C_nthLine,C_nthLine_default," ARM_ERR: nth Line: numeric expected",C_result);
	    XL_readNumCellWD(XL_nthCol,C_nthCol,C_nthCol_default," ARM_ERR: nth Column: numeric expected",C_result);
	    XL_readStrCellWD(XL_isCumul,C_isCumul,"NO"," ARM_ERR: is cumulative?: string expected",C_result);
	    XL_readStrCellWD(XL_isAbsolute,C_isAbsolute,"YES"," ARM_ERR: is absolute?: string expected",C_result);
		XL_readStrCellWD(XL_Currency,C_Currency,ARM_DEFAULT_COUNTRY," ARM_ERR: Currency: string expected",C_result);

		if((cumulId = ARM_ConvYesOrNo (C_isCumul, C_result)) == ARM_DEFAULT_ERR)
	    {
		    ARM_ARG_ERR();
		    return (LPXLOPER)&XL_result;
	    }

	    if((absoluteId = ARM_ConvYesOrNo (C_isAbsolute, C_result)) == ARM_DEFAULT_ERR)
	    {
		    ARM_ARG_ERR();
		    return (LPXLOPER)&XL_result;
	    }

		C_TypeCorr[0].toUpper();
		bool isNewCorrelManager = true;

		if (C_TypeCorr.size() > 1)
		{
			// ARM_CorrelManager instead of CorrelatorManager
			// not the same method : pass tags instead of correlation type
			// can be an array of {mktTag;intraMktTag}
			isNewCorrelManager = false;

			for (int i=0; i<C_TypeCorr.size()/2; i++)
			{
				C_mktTag.push_back(CCSTringToSTLString(C_TypeCorr[2*i]));
				C_intraMktTag.push_back(CCSTringToSTLString(C_TypeCorr[2*i+1]));
			}
		}
		else
		{
			if( C_TypeCorr[0] == "D" || C_TypeCorr[0] == "DIAG" )
			{
				TypeCorrId = 0;
			}
			else if( C_TypeCorr[0] == "I" || C_TypeCorr[0] == "INDEX" )
			{
				TypeCorrId = 1;
			}
			else if(C_TypeCorr[0] == "C"|| C_TypeCorr[0] == "CORR" )
			{
				TypeCorrId = 2;
			}
			else if(C_TypeCorr[0] == "S"|| C_TypeCorr[0] == "SIMPLE" )
			{
				TypeCorrId = 3;
			}
			else if(C_TypeCorr[0] == "F" || C_TypeCorr[0] == "FX" )
			{
				TypeCorrId = 6;
			}
			else
			{
				C_result.setMsg("ARM_ERR: typeCorrel should be one of [D or DIAG,I or INDEX,C or CORR, S or SIMPLE]");
				ARM_ARG_ERR();
				return (LPXLOPER)&XL_result;
			}
		}

		long retCode;
	    long objId;
	    CCString prevClass;

	    CCString curClass = LocalGetStringObjectClass (C_CorrelManager);
	    CCString stringId = GetLastCurCellEnvValue ();
		if(!stringId)
		{
			if (isNewCorrelManager)
			{
				retCode = ARMLOCAL_BumpVolatilityCorrelManager( LocalGetNumObjectId(C_CorrelManager), 
																CCSTringToSTLString(C_Currency), TypeCorrId, 
																C_Value, C_nthLine, C_nthCol, cumulId, 
																absoluteId, isToCloneId, C_result);
			}
			else
			{
				retCode = ARMLOCAL_BumpVolatilityCorrelManager( LocalGetNumObjectId(C_CorrelManager), 
																C_mktTag, C_intraMktTag, 
																C_Value, C_nthLine, C_nthCol, cumulId, 
																absoluteId, isToCloneId, C_result);
			}

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
				if (isNewCorrelManager)
				{
					retCode = ARMLOCAL_BumpVolatilityCorrelManager( LocalGetNumObjectId(C_CorrelManager), 
																	CCSTringToSTLString(C_Currency), TypeCorrId, 
																	C_Value, C_nthLine, C_nthCol, cumulId, 
																	absoluteId, isToCloneId, C_result, objId);
				}
				else
				{
					retCode = ARMLOCAL_BumpVolatilityCorrelManager( LocalGetNumObjectId(C_CorrelManager), 
																	C_mktTag, C_intraMktTag, 
																	C_Value, C_nthLine, C_nthCol, cumulId, 
																	absoluteId, isToCloneId, C_result, objId);
				}

				if (retCode == ARM_OK)
			    {
				    LocalSetCurCellEnvValue (curClass, objId); 

				    stringId = LocalMakeObjectId (objId, curClass);
			    }
		    }
			else
		    {
			    FreeCurCellContent ();

				if (isNewCorrelManager)
				{
				    retCode = ARMLOCAL_BumpVolatilityCorrelManager( LocalGetNumObjectId(C_CorrelManager), 
																	CCSTringToSTLString(C_Currency), TypeCorrId, 
																	C_Value, C_nthLine, C_nthCol, cumulId, 
																	absoluteId, isToCloneId, C_result);
				}
				else
				{
					retCode = ARMLOCAL_BumpVolatilityCorrelManager( LocalGetNumObjectId(C_CorrelManager), 
																	C_mktTag, C_intraMktTag, 
																	C_Value, C_nthLine, C_nthCol, cumulId, 
																	absoluteId, isToCloneId, C_result);
				}

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
		    FreeCurCellErr();

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

	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_BumpVolatilityCorrelManager" )

	return (LPXLOPER)&XL_result;
}

__declspec(dllexport) LPXLOPER WINAPI Local_PXL_BumpVolatilityCorrelManager(LPXLOPER XL_CorrelManager,
																			LPXLOPER XL_TypeCorr,
																		    LPXLOPER XL_Value,
																		    LPXLOPER XL_nthLine,
																		    LPXLOPER XL_nthCol,
																		    LPXLOPER XL_isCumul,
																		    LPXLOPER XL_isAbsolute,
																		    LPXLOPER XL_Currency/*,
																		    LPXLOPER XL_isToClone*/)

{
	ADD_LOG("Local_PXL_BumpVolatilityCorrelManager");
	static XLOPER XL_result;
	ARM_result C_result;

	ARM_XL_TRY_BLOCK_BEGIN
	{
	    ARM_NOCALCIFWIZ();
	
		CCString C_CorrelManager;
		
		double C_Value;

		double C_nthLine;
	    double C_nthLine_default = 0.;

	    double C_nthCol;
	    double C_nthCol_default = 0.;

	    CCString C_isCumul;
	    long cumulId;

	    CCString C_isAbsolute;
	    long absoluteId;

		CCString C_isToClone;
	    long isToCloneId = 1;

		CCString C_Currency;

		vector<CCString> C_TypeCorr;
		long TypeCorrId = 0;
		vector<string> C_mktTag;
		vector<string> C_intraMktTag;

		static int error;
	    static char* reason = "";

	    XL_readStrCell(XL_CorrelManager,C_CorrelManager," ARM_ERR: CorrelManager: object expected",C_result);
		XL_readStrVector(XL_TypeCorr,C_TypeCorr," ARM_ERR: correl type: string expected",DOUBLE_TYPE,C_result);
	    XL_readNumCell(XL_Value,C_Value," ARM_ERR: value to bump: numeric expected",C_result);
	    XL_readNumCellWD(XL_nthLine,C_nthLine,C_nthLine_default," ARM_ERR: nth Line: numeric expected",C_result);
	    XL_readNumCellWD(XL_nthCol,C_nthCol,C_nthCol_default," ARM_ERR: nth Column: numeric expected",C_result);
	    XL_readStrCellWD(XL_isCumul,C_isCumul,"NO"," ARM_ERR: is cumulative?: string expected",C_result);
	    XL_readStrCellWD(XL_isAbsolute,C_isAbsolute,"YES"," ARM_ERR: is absolute?: string expected",C_result);
		XL_readStrCellWD(XL_Currency,C_Currency,ARM_DEFAULT_COUNTRY," ARM_ERR: Currency: string expected",C_result);

		if((cumulId = ARM_ConvYesOrNo (C_isCumul, C_result)) == ARM_DEFAULT_ERR)
	    {
		    ARM_ARG_ERR();
		    return (LPXLOPER)&XL_result;
	    }

	    if((absoluteId = ARM_ConvYesOrNo (C_isAbsolute, C_result)) == ARM_DEFAULT_ERR)
	    {
		    ARM_ARG_ERR();
		    return (LPXLOPER)&XL_result;
	    }

		C_TypeCorr[0].toUpper();
		bool isNewCorrelManager = true;

		if (C_TypeCorr.size() > 1)
		{
			// ARM_CorrelManager instead of CorrelatorManager
			// not the same method : pass tags instead of correlation type
			// can be an array of {mktTag;intraMktTag}
			isNewCorrelManager = false;

			for (int i=0; i<C_TypeCorr.size()/2; i++)
			{
				C_mktTag.push_back(CCSTringToSTLString(C_TypeCorr[2*i]));
				C_intraMktTag.push_back(CCSTringToSTLString(C_TypeCorr[2*i+1]));
			}
		}
		else
		{
			if( C_TypeCorr[0] == "D" || C_TypeCorr[0] == "DIAG" )
			{
				TypeCorrId = 0;
			}
			else if( C_TypeCorr[0] == "I" || C_TypeCorr[0] == "INDEX" )
			{
				TypeCorrId = 1;
			}
			else if(C_TypeCorr[0] == "C"|| C_TypeCorr[0] == "CORR" )
			{
				TypeCorrId = 2;
			}
			else if(C_TypeCorr[0] == "S"|| C_TypeCorr[0] == "SIMPLE" )
			{
				TypeCorrId = 3;
			}
			else if(C_TypeCorr[0] == "F"|| C_TypeCorr[0] == "FX" )
			{
				TypeCorrId = 6;
			}
			else
			{
				C_result.setMsg("ARM_ERR: typeCorrel should be one of [D or DIAG,I or INDEX,C or CORR, S or SIMPLE]");
				ARM_ARG_ERR();
				return (LPXLOPER)&XL_result;
			}
		}

		long retCode;
	    long objId;
	    CCString prevClass;

	    CCString curClass = LocalGetStringObjectClass (C_CorrelManager);
	    CCString stringId;

		if (isNewCorrelManager)
		{
			retCode = ARMLOCAL_BumpVolatilityCorrelManager( LocalGetNumObjectId(C_CorrelManager), 
															CCSTringToSTLString(C_Currency), TypeCorrId, 
															C_Value, C_nthLine, C_nthCol, cumulId, 
															absoluteId, isToCloneId, C_result);
		}
		else
		{
			retCode = ARMLOCAL_BumpVolatilityCorrelManager( LocalGetNumObjectId(C_CorrelManager), 
															C_mktTag, C_intraMktTag, 
															C_Value, C_nthLine, C_nthCol, cumulId, 
															absoluteId, isToCloneId, C_result);
		}

		if(retCode == ARM_OK)
		{
			objId = C_result.getLong ();

			stringId = LocalMakeObjectId (objId, curClass);

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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_PXL_BumpVolatilityCorrelManager" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}



__declspec(dllexport) LPXLOPER WINAPI Local_GetCorrelDiagFromIndex(LPXLOPER XL_IndexIndexCubeId,
																   LPXLOPER XL_Tenor,
														      	   LPXLOPER XL_ListTenor)
{
	ADD_LOG("Local_GetCorrelDiagFromIndex");
	static XLOPER XL_result;
	ARM_result C_result;

	ARM_XL_TRY_BLOCK_BEGIN
	{
	    ARM_NOCALCIFWIZ();
	
		CCString C_IndexIndexCubeId;
		CCString C_Tenor;
		vector<CCString> C_ListTenor;
		vector<CCString> C_ListTenorDefault(0);

		static int error;
	    static char* reason = "";

		XL_readStrCell(XL_IndexIndexCubeId,C_IndexIndexCubeId," ARM_ERR: IndexIndexCubeId: string expected",C_result);
		XL_readStrCell(XL_Tenor,C_Tenor," ARM_ERR: Tenor: string expected",C_result);
		XL_readStrVectorWD(XL_ListTenor,C_ListTenor,C_ListTenorDefault,"ARM_ERR : List Tenor : Vecteur expcected",DOUBLE_TYPE, C_result);

		long retCode;
	    long objId;
	    CCString prevClass;

	    CCString curClass = LOCAL_VOL_CURVE_LIN_CLASS;
	    CCString stringId = GetLastCurCellEnvValue ();
		if(!stringId)
		{
			retCode = ARMLOCAL_GetCorrelDiag(LocalGetNumObjectId(C_IndexIndexCubeId), CCSTringToSTLString(C_Tenor),
											 C_ListTenor,C_result);
							
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
				retCode = ARMLOCAL_GetCorrelDiag(LocalGetNumObjectId(C_IndexIndexCubeId), CCSTringToSTLString(C_Tenor),
											     C_ListTenor,C_result,objId);

				if (retCode == ARM_OK)
			    {
				    LocalSetCurCellEnvValue (curClass, objId); 

				    stringId = LocalMakeObjectId (objId, curClass);
			    }
		    }
			else
		    {
			    FreeCurCellContent ();

			    retCode = ARMLOCAL_GetCorrelDiag(LocalGetNumObjectId(C_IndexIndexCubeId), CCSTringToSTLString(C_Tenor),
											     C_ListTenor,C_result);

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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_GetCorrelDiagFromIndex" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}

__declspec(dllexport) LPXLOPER WINAPI Local_CreateCorrelator( LPXLOPER XL_mktTags,
															  LPXLOPER XL_HyperDiagCurveIds,
															  LPXLOPER XL_IndexIndexCurveIds,
															  LPXLOPER XL_CorrelCurveIds,
															  LPXLOPER XL_IndexCurveIds,
															  LPXLOPER XL_IRVolHyperCubeIds,
															  LPXLOPER XL_VolVolHyperCubeIds,
															  LPXLOPER XL_FXVolHyperCubeIds)
{
	ADD_LOG("Local_CreateCorrelator");
	/// this is defined first because it is used in XL macros
	static XLOPER XL_result;
	/// Get the variables from the XLOper variables
	ARM_result C_result;

	// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
		/// to avoid computation if called by the wizard
		ARM_NOCALCIFWIZ();

		// error
		static int error;
		static char* reason = "";

		VECTOR<CCString> C_mktTags;
		VECTOR<CCString> C_HyperDiagVolIds;
		VECTOR<CCString> C_IndexIndexVolIds;
		VECTOR<CCString> C_CorrelVolIds;
		VECTOR<CCString> C_IndexVolIds;
		VECTOR<CCString> C_IRVolHyperCubeIds;
		VECTOR<CCString> C_VolVolHyperCubeIds;
		VECTOR<CCString> C_FXVolHyperCubeIds;
		VECTOR<CCString> C_DefaultVolIds(0);
		//C_DefaultVolIds[0] = "NULL";

		XL_readStrVector( XL_mktTags, C_mktTags," ARM_ERR: mkt Tag: array expected",DOUBLE_TYPE, C_result);
		XL_readStrVectorWD( XL_HyperDiagCurveIds,C_HyperDiagVolIds,C_DefaultVolIds," ARM_ERR: Hyper Diag id: array expected",DOUBLE_TYPE,	C_result);
		XL_readStrVectorWD( XL_IndexIndexCurveIds,C_IndexIndexVolIds,C_DefaultVolIds," ARM_ERR: IndexIndex id: string or array of string expected",DOUBLE_TYPE, C_result);
		XL_readStrVectorWD( XL_CorrelCurveIds,C_CorrelVolIds,C_DefaultVolIds," ARM_ERR: CORR id: string or array of string expected",DOUBLE_TYPE, C_result);
		XL_readStrVectorWD( XL_IndexCurveIds,C_IndexVolIds,C_DefaultVolIds," ARM_ERR: IndexCorr id: string or array of string expected",DOUBLE_TYPE, C_result);
		XL_readStrVectorWD( XL_IRVolHyperCubeIds,C_IRVolHyperCubeIds,C_DefaultVolIds," ARM_ERR: IR/VOL id: string or array of string expected",DOUBLE_TYPE, C_result);
		XL_readStrVectorWD( XL_VolVolHyperCubeIds,C_VolVolHyperCubeIds,C_DefaultVolIds," ARM_ERR: VOL/VOL id: string or array of string expected",DOUBLE_TYPE, C_result);
		XL_readStrVectorWD( XL_FXVolHyperCubeIds,C_FXVolHyperCubeIds,C_DefaultVolIds," ARM_ERR: FX/VOL id: string or array of string expected",DOUBLE_TYPE, C_result);


		int i;
		long sizeCurve = C_HyperDiagVolIds.size();
		vector<long> v_HyperDiagVolIds(sizeCurve);
		for(i = 0; i < sizeCurve; i++ )
			v_HyperDiagVolIds[i] = ( strcmp(C_HyperDiagVolIds[i], "NULL") == 0 ) ?
						           ARM_NULL_OBJECT : LocalGetNumObjectId(C_HyperDiagVolIds[i] );

		sizeCurve = C_IndexIndexVolIds.size();
		vector<long> v_IndexIndexVolIds(sizeCurve);
		for(i = 0; i < sizeCurve; i++ )
			v_IndexIndexVolIds[i] = ( strcmp(C_IndexIndexVolIds[i], "NULL") == 0 ) ?
						           ARM_NULL_OBJECT : LocalGetNumObjectId(C_IndexIndexVolIds[i] );

		sizeCurve = C_CorrelVolIds.size();
		vector<long> v_CorrelVolIds(sizeCurve);
		for(i = 0; i < sizeCurve; i++ )
			v_CorrelVolIds[i] = ( strcmp(C_CorrelVolIds[i], "NULL") == 0 ) ?
						           ARM_NULL_OBJECT : LocalGetNumObjectId(C_CorrelVolIds[i] );

		sizeCurve = C_IndexVolIds.size();
		vector<long> v_IndexVolIds(sizeCurve);
		for(i = 0; i < sizeCurve; i++ )
			v_IndexVolIds[i] = ( strcmp(C_IndexVolIds[i], "NULL") == 0 ) ?
						           ARM_NULL_OBJECT : LocalGetNumObjectId(C_IndexVolIds[i] );

		sizeCurve = C_IRVolHyperCubeIds.size();
		vector<long> vIRVolHyperCubeIds(sizeCurve);
		for(i = 0; i < sizeCurve; i++ )
			vIRVolHyperCubeIds[i] = ( strcmp(C_IRVolHyperCubeIds[i], "NULL") == 0 ) ?
						           ARM_NULL_OBJECT : LocalGetNumObjectId(C_IRVolHyperCubeIds[i] );

		sizeCurve = C_VolVolHyperCubeIds.size();
		vector<long> vVolVolHyperCubeIds(sizeCurve);
		for(i = 0; i < sizeCurve; i++ )
			vVolVolHyperCubeIds[i] = ( strcmp(C_VolVolHyperCubeIds[i], "NULL") == 0 ) ?
						           ARM_NULL_OBJECT : LocalGetNumObjectId(C_VolVolHyperCubeIds[i] );

		sizeCurve = C_FXVolHyperCubeIds.size();
		vector<long> vFXVolHyperCubeIds(sizeCurve);
		for(i = 0; i < sizeCurve; i++ )
			vFXVolHyperCubeIds[i] = ( strcmp(C_FXVolHyperCubeIds[i], "NULL") == 0 ) ?
						           ARM_NULL_OBJECT : LocalGetNumObjectId(C_FXVolHyperCubeIds[i] );


		long retCode;
		long objId;
		CCString prevClass;
		
		CCString curClass = LOCAL_CORRELMANAGER_CLASS;
		CCString stringId = GetLastCurCellEnvValue ();

		if(!stringId)
		{
			retCode =  ARMLOCAL_CreateGenCorrelatorManager(	C_mktTags,
															v_HyperDiagVolIds,
															v_IndexIndexVolIds,
															v_CorrelVolIds,
															v_IndexVolIds,
															vIRVolHyperCubeIds,
															vVolVolHyperCubeIds,
															vFXVolHyperCubeIds,
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
				retCode = ARMLOCAL_CreateGenCorrelatorManager (	C_mktTags,
																v_HyperDiagVolIds,
																v_IndexIndexVolIds,
																v_CorrelVolIds,
																v_IndexVolIds,
																vIRVolHyperCubeIds,
																vVolVolHyperCubeIds,
																vFXVolHyperCubeIds,
																C_result, 
																objId);

				if (retCode == ARM_OK)
				{
					LocalSetCurCellEnvValue (curClass, objId); 

					stringId = LocalMakeObjectId (objId, curClass);
				}
			}
			else
			{
				FreeCurCellContent ();
				retCode = ARMLOCAL_CreateGenCorrelatorManager (	C_mktTags,
																v_HyperDiagVolIds,
																v_IndexIndexVolIds,
																v_CorrelVolIds,
																v_IndexVolIds,
																vIRVolHyperCubeIds,
																vVolVolHyperCubeIds,
																vFXVolHyperCubeIds,
																C_result, 
																objId);
			
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_CreateCorrelator" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
	
}

__declspec(dllexport) LPXLOPER WINAPI Local_PXL_CreateCorrelator( LPXLOPER XL_mktTags,
															      LPXLOPER XL_HyperDiagCurveIds,
															      LPXLOPER XL_IndexIndexCurveIds,
															      LPXLOPER XL_CorrelCurveIds,
															      LPXLOPER XL_IndexCurveIds,
																  LPXLOPER XL_IRVolHyperCubeIds,
																  LPXLOPER XL_VolVolHyperCubeIds,
																  LPXLOPER XL_FXVolHyperCubeIds)
{
	ADD_LOG("Local_PXL_CreateCorrelator");
	/// this is defined first because it is used in XL macros
	static XLOPER XL_result;
	/// Get the variables from the XLOper variables
	ARM_result C_result;

	// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
		/// to avoid computation if called by the wizard
		ARM_NOCALCIFWIZ();

		// error
		static int error;
		static char* reason = "";

		VECTOR<CCString> C_mktTags;
		VECTOR<CCString> C_HyperDiagVolIds;
		VECTOR<CCString> C_IndexIndexVolIds;
		VECTOR<CCString> C_CorrelVolIds;
		VECTOR<CCString> C_IndexVolIds;
		VECTOR<CCString> C_IRVolHyperCubeIds;
		VECTOR<CCString> C_VolVolHyperCubeIds;
		VECTOR<CCString> C_FXVolHyperCubeIds;
		VECTOR<CCString> C_DefaultVolIds(0);
		//C_DefaultVolIds[0] = "NULL";

		XL_readStrVector( XL_mktTags, C_mktTags," ARM_ERR: mkt Tag: array expected",DOUBLE_TYPE, C_result);
		XL_readStrVector( XL_HyperDiagCurveIds,C_HyperDiagVolIds," ARM_ERR: Hyper Diag id: array expected",DOUBLE_TYPE,	C_result);
		XL_readStrVectorWD( XL_IndexIndexCurveIds,C_IndexIndexVolIds,C_DefaultVolIds," ARM_ERR: IndexIndex id: string or array of string expected",DOUBLE_TYPE, C_result);
		XL_readStrVectorWD( XL_CorrelCurveIds,C_CorrelVolIds,C_DefaultVolIds," ARM_ERR: CORR id: string or array of string expected",DOUBLE_TYPE, C_result);
		XL_readStrVectorWD( XL_IndexCurveIds,C_IndexVolIds,C_DefaultVolIds," ARM_ERR: IndexCorr id: string or array of string expected",DOUBLE_TYPE, C_result);
		XL_readStrVectorWD( XL_IRVolHyperCubeIds,C_IRVolHyperCubeIds,C_DefaultVolIds," ARM_ERR: IR/VOL id: string or array of string expected",DOUBLE_TYPE, C_result);
		XL_readStrVectorWD( XL_VolVolHyperCubeIds,C_VolVolHyperCubeIds,C_DefaultVolIds," ARM_ERR: VOL/VOL id: string or array of string expected",DOUBLE_TYPE, C_result);
		XL_readStrVectorWD( XL_FXVolHyperCubeIds,C_FXVolHyperCubeIds,C_DefaultVolIds," ARM_ERR: FX/VOL id: string or array of string expected",DOUBLE_TYPE, C_result);


		int i;
		long sizeCurve = C_HyperDiagVolIds.size();
		vector<long> v_HyperDiagVolIds(sizeCurve);
		for(i = 0; i < sizeCurve; i++ )
			v_HyperDiagVolIds[i] = ( strcmp(C_HyperDiagVolIds[i], "NULL") == 0 ) ?
						           ARM_NULL_OBJECT : LocalGetNumObjectId(C_HyperDiagVolIds[i] );

		sizeCurve = C_IndexIndexVolIds.size();
		vector<long> v_IndexIndexVolIds(sizeCurve);
		for(i = 0; i < sizeCurve; i++ )
			v_IndexIndexVolIds[i] = ( strcmp(C_IndexIndexVolIds[i], "NULL") == 0 ) ?
						           ARM_NULL_OBJECT : LocalGetNumObjectId(C_IndexIndexVolIds[i] );

		sizeCurve = C_CorrelVolIds.size();
		vector<long> v_CorrelVolIds(sizeCurve);
		for(i = 0; i < sizeCurve; i++ )
			v_CorrelVolIds[i] = ( strcmp(C_CorrelVolIds[i], "NULL") == 0 ) ?
						           ARM_NULL_OBJECT : LocalGetNumObjectId(C_CorrelVolIds[i] );

		sizeCurve = C_IndexVolIds.size();
		vector<long> v_IndexVolIds(sizeCurve);
		for(i = 0; i < sizeCurve; i++ )
			v_IndexVolIds[i] = ( strcmp(C_IndexVolIds[i], "NULL") == 0 ) ?
						           ARM_NULL_OBJECT : LocalGetNumObjectId(C_IndexVolIds[i] );

		sizeCurve = C_IRVolHyperCubeIds.size();
		vector<long> vIRVolHyperCubeIds(sizeCurve);
		for(i = 0; i < sizeCurve; i++ )
			vIRVolHyperCubeIds[i] = ( strcmp(C_IRVolHyperCubeIds[i], "NULL") == 0 ) ?
						           ARM_NULL_OBJECT : LocalGetNumObjectId(C_IRVolHyperCubeIds[i] );

		sizeCurve = C_VolVolHyperCubeIds.size();
		vector<long> vVolVolHyperCubeIds(sizeCurve);
		for(i = 0; i < sizeCurve; i++ )
			vVolVolHyperCubeIds[i] = ( strcmp(C_VolVolHyperCubeIds[i], "NULL") == 0 ) ?
						           ARM_NULL_OBJECT : LocalGetNumObjectId(C_VolVolHyperCubeIds[i] );

		sizeCurve = C_FXVolHyperCubeIds.size();
		vector<long> vFXVolHyperCubeIds(sizeCurve);
		for(i = 0; i < sizeCurve; i++ )
			vFXVolHyperCubeIds[i] = ( strcmp(C_FXVolHyperCubeIds[i], "NULL") == 0 ) ?
						           ARM_NULL_OBJECT : LocalGetNumObjectId(C_FXVolHyperCubeIds[i] );


		long retCode;
		long objId;


		CCString curClass = LOCAL_CORRELMANAGER_CLASS;
		CCString stringId;

		retCode =  ARMLOCAL_CreateGenCorrelatorManager(	C_mktTags,
														v_HyperDiagVolIds,
														v_IndexIndexVolIds,
														v_CorrelVolIds,
														v_IndexVolIds,
														vIRVolHyperCubeIds,
														vVolVolHyperCubeIds,
														vFXVolHyperCubeIds,
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
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_PXL_CreateCorrelator" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}

__declspec(dllexport) LPXLOPER WINAPI Local_ComputeIdIdCorrelFromManager(
    LPXLOPER XL_correlManagerId,
	LPXLOPER XL_tenor1,
	LPXLOPER XL_tenor2,
	LPXLOPER XL_expiry1,
	LPXLOPER XL_expiry2,
	LPXLOPER XL_currecy)
{
	ADD_LOG("Local_ComputeIdIdCorrelFromManager");
	static XLOPER XL_result;

	ARM_result C_result;

	ARM_XL_TRY_BLOCK_BEGIN
	{
		ARM_NOCALCIFWIZ();

		static int error;
		static char* reason = "";

		CCString C_Currency;
		CCString C_Tenor1;
		CCString C_Tenor2;
		double C_Expiry1;
		double C_Expiry2;
		long C_CorrelManagerId;

		XL_readStrCell( XL_tenor1, C_Tenor1, " ARM_ERR: tenor1 Tag: string expected", C_result );
		XL_readStrCell( XL_tenor2, C_Tenor2, " ARM_ERR: tenor2 Tag: string expected", C_result );
		XL_readNumCell( XL_expiry1, C_Expiry1, " ARM_ERR: Expiry1: double expected", C_result );
		XL_readNumCell( XL_expiry2, C_Expiry2, " ARM_ERR: Expiry2: double expected", C_result );
		XL_GETOBJID( XL_correlManagerId, C_CorrelManagerId,	" ARM_ERR: correl manager: Object expected",	C_result);
		XL_readStrCellWD( XL_currecy, C_Currency, ARM_DEFAULT_COUNTRY, " ARM_ERR: currency Tag: string expected", C_result );

		long retCode;
		if( CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 1 )
			retCode = ARMLOCAL_ComputeIdIdCorrelFromCorrelatorManager( CCSTringToSTLString( C_Currency ), CCSTringToSTLString( C_Tenor1 ), CCSTringToSTLString( C_Tenor2 ),
																	   C_Expiry1, C_Expiry2, C_CorrelManagerId, C_result);
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

	}

	ARM_XL_TRY_BLOCK_END

	ARM_XL_CATCH_ARM_EXPT

	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_ComputeIdIdCorrelFromManager" )

	return (LPXLOPER)&XL_result;
}

__declspec(dllexport) LPXLOPER WINAPI Local_ComputeHyperCorrelFromManager(
    LPXLOPER XL_correlManagerId,
	LPXLOPER XL_tenor1,
	LPXLOPER XL_tenor2,
	LPXLOPER XL_expiry,
	LPXLOPER XL_currecy,
	LPXLOPER XL_byExpiry)
{
	ADD_LOG("Local_ComputeHyperCorrelFromManager");
	static XLOPER XL_result;

	ARM_result C_result;

	ARM_XL_TRY_BLOCK_BEGIN
	{
		ARM_NOCALCIFWIZ();

		static int error;
		static char* reason = "";

		CCString C_Currency;
		CCString C_Tenor1;
		CCString C_Tenor2;
		CCString C_byExpiry;
		double C_Expiry;
		long C_CorrelManagerId;
		bool isByExpiry = false;

		XL_readStrCell( XL_tenor1, C_Tenor1, " ARM_ERR: tenor1 Tag: string expected", C_result );
		XL_readStrCell( XL_tenor2, C_Tenor2, " ARM_ERR: tenor2 Tag: string expected", C_result );
		XL_readNumCell( XL_expiry, C_Expiry, " ARM_ERR: Expiry: double expected", C_result );
		XL_GETOBJID( XL_correlManagerId, C_CorrelManagerId,	" ARM_ERR: correl manager: Object expected",	C_result);
		XL_readStrCellWD( XL_currecy, C_Currency, ARM_DEFAULT_COUNTRY, " ARM_ERR: currency Tag: string expected", C_result );
		XL_readStrCellWD( XL_byExpiry, C_byExpiry, "N", " ARM_ERR: byExpiry flag Y/N: string expected", C_result );

		C_byExpiry.toUpper();
		if(CCSTringToSTLString(C_byExpiry)!="N")
            isByExpiry=true;

		long retCode;
		if( CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 1 )
			retCode = ARMLOCAL_ComputeHyperCorrelFromCorrelatorManager( CCSTringToSTLString( C_Currency ), 
																		CCSTringToSTLString( C_Tenor1 ),
																		CCSTringToSTLString( C_Tenor2 ),
																	    C_Expiry, C_CorrelManagerId, C_result, isByExpiry );
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

	}

	ARM_XL_TRY_BLOCK_END

	ARM_XL_CATCH_ARM_EXPT

	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_ComputeHyperCorrelFromManager" )

	return (LPXLOPER)&XL_result;
}

__declspec(dllexport) LPXLOPER WINAPI Local_ComputeCorrSimpleFromManager(LPXLOPER XL_correlManagerId,
																		 LPXLOPER XL_tenor,
																		 LPXLOPER XL_expiry,
																		 LPXLOPER XL_currency,
																		 LPXLOPER XL_Corr)
{
	ADD_LOG("Local_ComputeCorrSimpleFromManager");
	static XLOPER XL_result;

	ARM_result C_result;

	ARM_XL_TRY_BLOCK_BEGIN
	{
		ARM_NOCALCIFWIZ();

		static int error;
		static char* reason = "";

		CCString C_Currency;
		CCString C_Tenor;
		CCString C_Corr;
		double C_Expiry;
		long C_CorrelManagerId;
		bool isCorr = true;

		XL_readStrCell( XL_tenor, C_Tenor, " ARM_ERR: tenor1 Tag: string expected", C_result );
		XL_readNumCell( XL_expiry, C_Expiry, " ARM_ERR: Expiry: double expected", C_result );
		XL_GETOBJID( XL_correlManagerId, C_CorrelManagerId,	" ARM_ERR: correl manager: Object expected",	C_result);
		XL_readStrCellWD( XL_currency, C_Currency, ARM_DEFAULT_COUNTRY, " ARM_ERR: currency Tag: string expected", C_result );
		XL_readStrCellWD( XL_Corr, C_Corr, "Y", " ARM_ERR: Corr flag Y/N: string expected", C_result );

		C_Corr.toUpper();
		if(CCSTringToSTLString(C_Corr)!="Y")
            isCorr=false;

		long retCode;
		if( CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 1 )
			retCode = ARMLOCAL_ComputeCorrSimplFromCorrelatorManager( CCSTringToSTLString( C_Currency ), 
																CCSTringToSTLString( C_Tenor ),
																C_Expiry, C_CorrelManagerId, C_result, isCorr );
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

	}

	ARM_XL_TRY_BLOCK_END

	ARM_XL_CATCH_ARM_EXPT

	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_ComputeCorrSimpleFromManager" )

	return (LPXLOPER)&XL_result;
}

__declspec(dllexport) LPXLOPER WINAPI Local_ComputeFXCorrelFromManager(LPXLOPER XL_correlManagerId,
																	   LPXLOPER XL_currency1,
																	   LPXLOPER XL_currency2,
																	   LPXLOPER XL_tenor,
																	   LPXLOPER XL_expiry)
{
	ADD_LOG("Local_ComputeFXCorrelFromManager");
	static XLOPER XL_result;

	ARM_result C_result;

	ARM_XL_TRY_BLOCK_BEGIN
	{
		ARM_NOCALCIFWIZ();

		static int error;
		static char* reason = "";

		CCString C_Currency1;
		CCString C_Currency2;
		CCString C_Tenor;
		double C_Expiry;
		long C_CorrelManagerId;

		XL_readStrCell( XL_tenor, C_Tenor, " ARM_ERR: tenor1 Tag: string expected", C_result );
		XL_readNumCell( XL_expiry, C_Expiry, " ARM_ERR: Expiry: double expected", C_result );
		XL_GETOBJID( XL_correlManagerId, C_CorrelManagerId,	" ARM_ERR: correl manager: Object expected",	C_result);
		XL_readStrCell( XL_currency1, C_Currency1, " ARM_ERR: currency1 Tag: string expected", C_result );
		XL_readStrCell( XL_currency2, C_Currency2, " ARM_ERR: currency2 Tag: string expected", C_result );

		long retCode;
		if( CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 1 )
			retCode = ARMLOCAL_ComputeFXCorrelFromCorrelatorManager(CCSTringToSTLString( C_Currency1 ),
																	CCSTringToSTLString( C_Currency2 ),
																	CCSTringToSTLString( C_Tenor ),
																	C_Expiry, C_CorrelManagerId, C_result);
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

	}

	ARM_XL_TRY_BLOCK_END

	ARM_XL_CATCH_ARM_EXPT

	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_ComputeFXCorrelFromManager" )

	return (LPXLOPER)&XL_result;
}

__declspec(dllexport) LPXLOPER WINAPI Local_PXL_GetCorrelDiagFromIndex(LPXLOPER XL_IndexIndexCubeId,
																       LPXLOPER XL_Tenor,
														      	       LPXLOPER XL_ListTenor)
{
	ADD_LOG("Local_PXL_GetCorrelDiagFromIndex");
	static XLOPER XL_result;
	ARM_result C_result;

	ARM_XL_TRY_BLOCK_BEGIN
	{
	    ARM_NOCALCIFWIZ();
	
		CCString C_IndexIndexCubeId;
		CCString C_Tenor;
		vector<CCString> C_ListTenor;
		vector<CCString> C_ListTenorDefault(0);

		static int error;
	    static char* reason = "";

		XL_readStrCell(XL_IndexIndexCubeId,C_IndexIndexCubeId," ARM_ERR: IndexIndexCubeId: string expected",C_result);
		XL_readStrCell(XL_Tenor,C_Tenor," ARM_ERR: Tenor: string expected",C_result);
		XL_readStrVectorWD(XL_ListTenor,C_ListTenor,C_ListTenorDefault,"ARM_ERR : List Tenor : Vecteur expcected",DOUBLE_TYPE, C_result);

		long retCode;
	    long objId;
	    CCString prevClass;

	    CCString curClass = LOCAL_VOL_CURVE_LIN_CLASS;
	    CCString stringId;

		retCode = ARMLOCAL_GetCorrelDiag(LocalGetNumObjectId(C_IndexIndexCubeId), CCSTringToSTLString(C_Tenor),
											 C_ListTenor,C_result);
							
		if(retCode == ARM_OK)
		{
			objId = C_result.getLong ();

			stringId = LocalMakeObjectId (objId, curClass);
		
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_PXL_GetCorrelDiagFromIndex" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


__declspec(dllexport) LPXLOPER WINAPI Local_GetMixtureParamsFromSummit ( LPXLOPER XL_index,
																		 LPXLOPER XL_currency,
																		 LPXLOPER XL_cvName,
																		 LPXLOPER XL_date,
																		 LPXLOPER XL_interpolMeth )
{
	ADD_LOG("Local_GetMixtureParamsFromSummit ");
	static XLOPER XL_result;
	ARM_result C_result;

	ARM_XL_TRY_BLOCK_BEGIN
	{
	    ARM_NOCALCIFWIZ();

	    // C variable
	    CCString C_index;
	    CCString C_currency;
	    CCString C_cvName;
	    double C_date;
		CCString C_interpolMeth;

	    // error
	    static int error;
	    static char* reason = "";

	    XL_readStrCell(XL_index,C_index," ARM_ERR: index: string expected",C_result);
	    XL_readStrCell(XL_currency,C_currency," ARM_ERR: currency: string expected",C_result);
	    XL_readStrCell(XL_cvName,C_cvName," ARM_ERR: curve name: string expected",C_result);
	    XL_readNumCell(XL_date,C_date," ARM_ERR: as of date: date expected",C_result);
	    XL_readStrCellWD(XL_interpolMeth,C_interpolMeth,"LINEAR"," ARM_ERR: curve name: string expected",C_result);

	    C_index.toUpper();
	    C_currency.toUpper();
	    C_cvName.toUpper();

        if ( (CCSTringToSTLString(C_interpolMeth) != "LINEAR") 
			&& (CCSTringToSTLString(C_interpolMeth) != "STEPUPRIGHT") 
			&& (CCSTringToSTLString(C_interpolMeth) != "STEPUPLEFT") )
			throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT, "interpol method should be 'LINEAR', 'STEPUPRIGHT', STEPUPLEFT'...");

		long retCode;
	    long objId;
	    CCString prevClass;
	    
	    CCString stringId = GetLastCurCellEnvValue ();
		CCString curClass = LOCAL_MIXTURE_PARAMS_CLASS;
	    
	    if (!stringId)
	    {
		    retCode = ARMLOCAL_GetMixtureParamsFromSummit(C_index,
														  C_currency,
														  C_cvName,
														  C_date,
														  C_interpolMeth,
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
			    retCode = ARMLOCAL_GetMixtureParamsFromSummit(C_index,
															  C_currency,
															  C_cvName,
															  C_date,
															  C_interpolMeth,
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
			    retCode = ARMLOCAL_GetMixtureParamsFromSummit(C_index,
															  C_currency,
															  C_cvName,
															  C_date,
															  C_interpolMeth,
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
	}
	ARM_XL_TRY_BLOCK_END

	ARM_XL_CATCH_ARM_EXPT

	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in GetMixtureParamsFromSummit" )

	return (LPXLOPER)&XL_result;
}


__declspec(dllexport) LPXLOPER WINAPI Local_PXL_GetMixtureParamsFromSummit( LPXLOPER XL_index,
																			LPXLOPER XL_currency,
																			LPXLOPER XL_cvName,
																			LPXLOPER XL_date,
																			LPXLOPER XL_interpolMeth )
{
	ADD_LOG("Local_PXL_GetMixtureParamsFromSummit");
	static XLOPER XL_result;
	ARM_result C_result;

	ARM_XL_TRY_BLOCK_BEGIN
	{
	    ARM_NOCALCIFWIZ();

	    // C variable
	    CCString C_index;
	    CCString C_currency;
	    CCString C_cvName;
	    double C_date;
		CCString C_interpolMeth;

	    // error
	    static int error;
	    static char* reason = "";

	    XL_readStrCell(XL_index,C_index," ARM_ERR: index: string expected",C_result);
	    XL_readStrCell(XL_currency,C_currency," ARM_ERR: currency: string expected",C_result);
	    XL_readStrCell(XL_cvName,C_cvName," ARM_ERR: curve name: string expected",C_result);
	    XL_readNumCell(XL_date,C_date," ARM_ERR: as of date: date expected",C_result);
	    XL_readStrCellWD(XL_interpolMeth,C_interpolMeth,"LINEAR"," ARM_ERR: curve name: string expected",C_result);

	    C_index.toUpper();
	    C_currency.toUpper();
	    C_cvName.toUpper();

        if ( (CCSTringToSTLString(C_interpolMeth) != "LINEAR") 
			&& (CCSTringToSTLString(C_interpolMeth) != "STEPUPRIGHT") 
			&& (CCSTringToSTLString(C_interpolMeth) != "STEPUPLEFT") )
			throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT, "interpol method should be 'LINEAR', 'STEPUPRIGHT', STEPUPLEFT'...");

		long retCode;
	    long objId;
	    CCString prevClass;
	    
	    CCString stringId = GetLastCurCellEnvValue ();
		CCString curClass = LOCAL_MIXTURE_PARAMS_CLASS;
	    
		retCode = ARMLOCAL_GetMixtureParamsFromSummit(C_index,
													  C_currency,
													  C_cvName,
													  C_date,
													  C_interpolMeth,
													  C_result);

	    if ( retCode == ARM_OK )
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

	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in PXL_GetMixtureParamsFromSummit" )

	return (LPXLOPER)&XL_result;
}


//--------------------------------------------------------------------------//
// Building of a SABR volatility thanks to a 4 volatility curves :			//
// SigmaOrAlpha, Rho, Beta, Nu												//
//--------------------------------------------------------------------------//
__declspec(dllexport) LPXLOPER WINAPI Local_SABRVol(LPXLOPER XL_SigmaOrAlphaId,
													LPXLOPER XL_RhoId,
													LPXLOPER XL_NuId,
													LPXLOPER XL_BetaId,
													LPXLOPER XL_SOrAFlag,
													LPXLOPER XL_ModelType,
													LPXLOPER XL_Weight)
{
	ADD_LOG("Local_SABRVol");
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	    ARM_NOCALCIFWIZ();


	    // C variables
		CCString C_SigmaOrAlpha;
		CCString C_Rho;
		CCString C_Beta;
		CCString C_Nu;

		CCString C_ModType;
		long     ModelType;

		long   C_Beta_Id;

        double C_SOrA;
		double C_SOrA_default = 1.0;		// Default : Sigma
		long   SOrA_Flag;

		double C_Weight;
		double C_Weight_default = 0.5;

	    // error
	    static int error;
	    static char* reason = "";

		XL_readStrCellWD(XL_SigmaOrAlphaId, C_SigmaOrAlpha, "", " ARM_ERR: SigmaOrAlpha: string expected", C_result);
		XL_readStrCellWD(XL_RhoId, C_Rho, "", " ARM_ERR: Rho: string expected", C_result);
		XL_readStrCellWD(XL_BetaId, C_Beta, "", " ARM_ERR: Beta: string expected", C_result);
		XL_readStrCellWD(XL_NuId, C_Nu, "", " ARM_ERR: Nu: string expected", C_result);	    
		XL_readNumCellWD(XL_SOrAFlag, C_SOrA, C_SOrA_default, " ARM_ERR: Sigma or Alpha Flag: numeric expected",C_result);
		XL_readStrCellWD(XL_ModelType, C_ModType, "LD", " ARM_ERR: smiled model type (LD, SABR_G, SABR_A): string expected", C_result);
		XL_readNumCellWD(XL_Weight, C_Weight, C_Weight_default, " ARM_ERR: Weight: numeric expected", C_result);

		
		if ( (XL_BetaId->xltype == xltypeMissing) || (XL_BetaId->xltype == xltypeNil) )
		{
			C_Beta_Id	=	ARM_NULL_OBJECT;
		}
		else
		{
			C_Beta_Id	=	(long) LocalGetNumObjectId(C_Beta);
		}

		SOrA_Flag	=	(long) C_SOrA;

	    if ( (ModelType = ARM_ConvSmiledModelFlag(C_ModType, C_result)) == ARM_DEFAULT_ERR )
	    {
	       ARM_ARG_ERR();

	       return (LPXLOPER)&XL_result;
	    }

	    long retCode;
	    long objId;
	    CCString prevClass;

	    CCString curClass = LOCAL_SABR_VOL_CLASS;
	    CCString stringId = GetLastCurCellEnvValue ();
	    
	    if(!stringId)
	    {
		    retCode = ARMLOCAL_SABRVol(LocalGetNumObjectId(C_SigmaOrAlpha),
									   LocalGetNumObjectId(C_Rho),
									   C_Beta_Id,
									   LocalGetNumObjectId(C_Nu),
									   SOrA_Flag,
									   ModelType,
									   C_Weight,
									   C_result);

		    if ( retCode == ARM_OK )
		    {
			    objId = C_result.getLong();

			    LocalSetCurCellEnvValue(curClass, objId); 

			    stringId = LocalMakeObjectId(objId, curClass);
		    }
	    }
	    else
	    {
		    prevClass = LocalGetStringObjectClass (stringId);
		    
		    objId = LocalGetNumObjectId (stringId);
			    
		    if ( curClass == prevClass )
		    {
				retCode = ARMLOCAL_SABRVol(LocalGetNumObjectId(C_SigmaOrAlpha),
										   LocalGetNumObjectId(C_Rho),
										   C_Beta_Id,
										   LocalGetNumObjectId(C_Nu),
										   SOrA_Flag,
										   ModelType,
										   C_Weight,
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
			    
				retCode = ARMLOCAL_SABRVol(LocalGetNumObjectId(C_SigmaOrAlpha),
										   LocalGetNumObjectId(C_Rho),
										   C_Beta_Id,
										   LocalGetNumObjectId(C_Nu),
										   SOrA_Flag,
										   ModelType,
										   C_Weight,
										   C_result);
		    
			    if ( retCode == ARM_OK )
			    {
				    objId = C_result.getLong ();
			    
				    LocalSetCurCellEnvValue(curClass, objId); 

				    stringId = LocalMakeObjectId(objId, curClass);
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_SABRVol" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


//--------------------------------------------------------------------------//
// Building of a SABR volatility thanks to a 4 volatility curves :			//
// SigmaOrAlpha, Rho, Beta, Nu	-> for Excel VBA							//
//--------------------------------------------------------------------------//
__declspec(dllexport) LPXLOPER WINAPI Local_PXL_SABRVol(LPXLOPER XL_SigmaOrAlphaId,
														LPXLOPER XL_RhoId,
														LPXLOPER XL_NuId,
														LPXLOPER XL_BetaId,
														LPXLOPER XL_SOrAFlag,
														LPXLOPER XL_ModelType,
														LPXLOPER XL_Weight)
{
	ADD_LOG("Local_PXL_SABRVol");
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	    ARM_NOCALCIFWIZ();


	    // C variables
		CCString C_SigmaOrAlpha;
		CCString C_Rho;
		CCString C_Beta;
		CCString C_Nu;

		CCString C_ModType;
		long     ModelType;

		long   C_Beta_Id;

        double C_SOrA;
		double C_SOrA_default = 1.0; // Default : Sigma
		long   SOrA_Flag;

		double C_Weight;
		double C_Weight_default = 0.5;

	    // error
	    static int error;
	    static char* reason = "";

		XL_readStrCellWD(XL_SigmaOrAlphaId, C_SigmaOrAlpha, "", " ARM_ERR: SigmaOrAlpha: string expected", C_result);
		XL_readStrCellWD(XL_RhoId, C_Rho, "", " ARM_ERR: Rho: string expected", C_result);
		XL_readStrCellWD(XL_BetaId, C_Beta, "", " ARM_ERR: Beta: string expected", C_result);
		XL_readStrCellWD(XL_NuId, C_Nu, "", " ARM_ERR: Nu: string expected", C_result);	    
		XL_readNumCellWD(XL_SOrAFlag, C_SOrA, C_SOrA_default, " ARM_ERR: Sigma or Alpha Flag: numeric expected",C_result);
		XL_readStrCellWD(XL_ModelType, C_ModType, "LD", " ARM_ERR: smiled model type (LD, SABR_G, SABR_A): string expected", C_result);
		XL_readNumCellWD(XL_Weight, C_Weight, C_Weight_default, " ARM_ERR: Weight: numeric expected", C_result);


		if ( (XL_BetaId->xltype == xltypeMissing) || (XL_BetaId->xltype == xltypeNil) )
		{
			C_Beta_Id	=	ARM_NULL_OBJECT;
		}
		else
		{
			C_Beta_Id	=	(long) LocalGetNumObjectId(C_Beta);
		}

		SOrA_Flag	=	(long) C_SOrA;

	    if ( (ModelType = ARM_ConvSmiledModelFlag(C_ModType, C_result)) == ARM_DEFAULT_ERR )
	    {
	       ARM_ARG_ERR();

	       return (LPXLOPER)&XL_result;
	    }

	    long retCode;
	    long objId;
	    CCString prevClass;

	    CCString curClass = LOCAL_SABR_VOL_CLASS;
	    CCString stringId;
	    
		retCode = ARMLOCAL_SABRVol(LocalGetNumObjectId(C_SigmaOrAlpha),
								   LocalGetNumObjectId(C_Rho),
								   C_Beta_Id,
								   LocalGetNumObjectId(C_Nu),
								   SOrA_Flag,
								   ModelType,
								   C_Weight,
								   C_result);

		if ( retCode == ARM_OK )
		{
			objId = C_result.getLong(); 

			stringId = LocalMakeObjectId(objId, curClass);
		}

	    if ( retCode == ARM_OK )
	    {
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_PXL_SABRVol" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


//--------------------------------------------------------------------------//
// Building of a SABR volatility structure thanks to Summit					//
//--------------------------------------------------------------------------//
__declspec(dllexport) LPXLOPER WINAPI Local_GetSABRVolFromSummit(LPXLOPER XL_index,
																 LPXLOPER XL_currency,
																 LPXLOPER XL_cvName,
																 LPXLOPER XL_date,
																 LPXLOPER XL_vtype,
																 LPXLOPER XL_matIndex,
																 LPXLOPER XL_impOrHist,
																 LPXLOPER XL_indexId,
																 LPXLOPER XL_sigmaOrAlpha,
																 LPXLOPER XL_modelType,
																 LPXLOPER XL_weight)
{
	ADD_LOG("Local_GetSABRVolFromSummit");
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
	    double C_date;
	    CCString C_vtype;
	    CCString C_matIndex;
	    CCString C_impOrHist;

		CCString C_ModType;
		long     ModelType;

		double   C_Flag;
		double   C_DefaultFlag = 1.0;
		long	 SigmaOrAlphaFlag;

		double   C_weight;
		double   C_weight_default = 0.5;

	    // error
	    static int error;
	    static char* reason = "";

	    XL_readStrCell(XL_index,C_index," ARM_ERR: index: string expected",C_result);
	    XL_readStrCell(XL_currency,C_currency," ARM_ERR: currency: string expected",C_result);
	    XL_readStrCell(XL_cvName,C_cvName," ARM_ERR: curve name: string expected",C_result);
	    XL_readNumCell(XL_date,C_date," ARM_ERR: as of date: date expected",C_result);
	    XL_readStrCell(XL_vtype,C_vtype," ARM_ERR: volatility type: string expected",C_result);
	    XL_readStrCellWD(XL_matIndex,C_matIndex,"ATM"," ARM_ERR: maturity index: string expected",C_result);	
	    XL_readStrCellWD(XL_impOrHist,C_impOrHist,"IRFWDVOL"," ARM_ERR: implied or Historical?: string expected",C_result);
		XL_readNumCellWD(XL_sigmaOrAlpha, C_Flag, C_DefaultFlag, " ARM_ERR: boolean expected", C_result);
		XL_readStrCellWD(XL_modelType, C_ModType, "LD", " ARM_ERR: smiled model type (LD, SABR_G, SABR_A): string expected", C_result);
		XL_readNumCellWD(XL_weight, C_weight, C_weight_default, " ARM_ERR: weight: numeric expected", C_result);

		std::string indexId; ExcelTools::convert(XL_indexId,"",indexId); 
		long indexId_ = LocalPersistent::get().getObjectId(indexId); 

	    C_index.toUpper ();
	    C_currency.toUpper ();
	    C_cvName.toUpper ();
	    C_matIndex.toUpper ();
	    C_impOrHist.toUpper ();

	    if ( (ModelType = ARM_ConvSmiledModelFlag(C_ModType, C_result)) == ARM_DEFAULT_ERR )
	    {
	       ARM_ARG_ERR();

	       return (LPXLOPER)&XL_result;
	    }

	    long retCode;
	    long objId;
	    CCString prevClass;
	    
	    CCString curClass = LOCAL_SABR_VOL_CLASS;
	    CCString stringId = GetLastCurCellEnvValue ();

		SigmaOrAlphaFlag  = (long) C_Flag;
	    
	    if (!stringId)
	    {
		    retCode = ARMLOCAL_GetSABRVolFromSummit (C_index,
												     C_currency,
												     C_cvName,
												     C_date,
													 C_vtype,
													 C_matIndex,
													 C_impOrHist,
													 indexId_,
													 SigmaOrAlphaFlag,
													 ModelType,
													 C_weight,
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
			    retCode = ARMLOCAL_GetSABRVolFromSummit (C_index,
														 C_currency,
														 C_cvName,
													     C_date,
												         C_vtype,
												         C_matIndex,
												         C_impOrHist,
												         indexId_,
														 SigmaOrAlphaFlag,
														 ModelType,
														 C_weight,
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
			    retCode = ARMLOCAL_GetSABRVolFromSummit (C_index,
													     C_currency,
														 C_cvName,
														 C_date,
														 C_vtype,
														 C_matIndex,
														 C_impOrHist,
														 indexId_,
														 SigmaOrAlphaFlag,
														 ModelType,
														 C_weight,
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_GetSABRVolFromSummit" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}



//--------------------------------------------------------------------------//
// Building of a SABR volatility structure thanks to Summit					//
// -> for Excel VBA															//
//--------------------------------------------------------------------------//
__declspec(dllexport) LPXLOPER WINAPI Local_PXL_GetSABRVolFromSummit(LPXLOPER XL_index,
																	 LPXLOPER XL_currency,
																	 LPXLOPER XL_cvName,
																	 LPXLOPER XL_date,
																	 LPXLOPER XL_vtype,
																	 LPXLOPER XL_matIndex,
																	 LPXLOPER XL_impOrHist,
																	 LPXLOPER XL_indexId,
																	 LPXLOPER XL_sigmaOrAlpha,
																	 LPXLOPER XL_modelType,
																	 LPXLOPER XL_weight)
{
	ADD_LOG("Local_PXL_GetSABRVolFromSummit");
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
	    double C_date;
	    CCString C_vtype;
	    CCString C_matIndex;
	    CCString C_impOrHist;

		CCString C_ModType;
		long     ModelType;

		double   C_Flag;
		double   C_DefaultFlag = 1.0;
		long     SigmaOrAlphaFlag;

		double   C_weight;
		double   C_weight_default = 0.5;

	    // error
	    static int error;
	    static char* reason = "";

	    XL_readStrCell(XL_index,C_index," ARM_ERR: index: string expected",C_result);
	    XL_readStrCell(XL_currency,C_currency," ARM_ERR: currency: string expected",C_result);
	    XL_readStrCell(XL_cvName,C_cvName," ARM_ERR: curve name: string expected",C_result);
	    XL_readNumCell(XL_date,C_date," ARM_ERR: as of date: date expected",C_result);
	    XL_readStrCell(XL_vtype,C_vtype," ARM_ERR: volatility type: string expected",C_result);
	    XL_readStrCellWD(XL_matIndex,C_matIndex,"ATM"," ARM_ERR: maturity index: string expected",C_result);	
	    XL_readStrCellWD(XL_impOrHist,C_impOrHist,"IRFWDVOL"," ARM_ERR: implied or Historical?: string expected",C_result);
		XL_readNumCellWD(XL_sigmaOrAlpha, C_Flag, C_DefaultFlag, " ARM_ERR: boolean expected", C_result);
		XL_readStrCellWD(XL_modelType, C_ModType, "LD", " ARM_ERR: smiled model type (LD, SABR_G, SABR_A): string expected", C_result);
		XL_readNumCellWD(XL_weight, C_weight, C_weight_default, " ARM_ERR: weight: numeric expected", C_result);

		std::string indexId; ExcelTools::convert(XL_indexId,"",indexId); 
		long indexId_ = LocalPersistent::get().getObjectId(indexId); 

	    C_index.toUpper ();
	    C_currency.toUpper ();
	    C_cvName.toUpper ();
	    C_matIndex.toUpper ();
	    C_impOrHist.toUpper ();

	    if ( (ModelType = ARM_ConvSmiledModelFlag(C_ModType, C_result)) == ARM_DEFAULT_ERR )
	    {
	       ARM_ARG_ERR();

	       return (LPXLOPER)&XL_result;
	    }

	    long retCode;
	    long objId;
	    
	    CCString curClass = LOCAL_SABR_VOL_CLASS;
		CCString stringId;

		SigmaOrAlphaFlag  = (long) C_Flag;
	
		retCode = ARMLOCAL_GetSABRVolFromSummit (C_index,
												 C_currency,
												 C_cvName,
												 C_date,
												 C_vtype,
												 C_matIndex,
												 C_impOrHist,
												 indexId_,
												 C_Flag,
												 ModelType,
												 C_weight,
												 C_result);

		if ( retCode == ARM_OK )
		{
			objId = C_result.getLong (); 

			stringId = LocalMakeObjectId (objId, curClass);
		}
	    
	    if ( retCode == ARM_OK )
	    {
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_PXL_GetSABRVolFromSummit" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


__declspec(dllexport) LPXLOPER WINAPI Local_OldVolCurve_Common(
		LPXLOPER XL_volCurve,
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
	
		long C_volCurveId;
		XL_GETOBJID( XL_volCurve,	C_volCurveId,	" ARM_ERR: vol Curve Id: Object expected",	C_result);

		exportFunc2Args< long, double >
			ourFunc(C_volCurveId, C_asOfDate, ARMLOCAL_OldVolCurve);

		bool PersistentInXL = true;

		/// call the general function
		fillXL_Result( LOCAL_VOL_CURVE_LIN_CLASS, ourFunc, C_result, XL_result, PersistentInXL );
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


__declspec(dllexport) LPXLOPER WINAPI Local_OldVolCurve(
	LPXLOPER XL_volCurve,
	LPXLOPER XL_asOfDate)
{
	bool PersistentInXL = true;
	
	return Local_OldVolCurve_Common(
			XL_volCurve,
			XL_asOfDate,
			PersistentInXL);
}

__declspec(dllexport) LPXLOPER WINAPI Local_PXL_OldVolCurve(
	LPXLOPER XL_volCurve,
	LPXLOPER XL_asOfDate)
{
	bool PersistentInXL = false;
			
	return Local_OldVolCurve_Common(
			XL_volCurve,
			XL_asOfDate,
			PersistentInXL);
}

__declspec(dllexport) LPXLOPER WINAPI Local_ARM_ConvIndexInYearTerm(LPXLOPER XL_index,
															           LPXLOPER XL_asof,
															           LPXLOPER XL_ccy)
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
	double C_asof;
	CCString C_currency;
	CCString C_index;

	// error
	static int error;
	static char* reason = "";
	
	XL_readStrCellWD(XL_index,C_index,"DEFAULT"," ARM_ERR: index: string expected",C_result);
	XL_readNumCell(XL_asof,C_asof," ARM_ERR: asofdate: numeric expected",C_result);
	XL_readStrCellWD(XL_ccy,C_currency,"DEFAULT"," ARM_ERR: currency: string expected",C_result);
	
	long retCode = ARMLOCAL_ConvIndexInYearTerm(C_index,C_asof,C_currency,C_result);

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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_ARM_ConvIndexInYearTerm" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}