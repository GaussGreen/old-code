#pragma warning(disable : 4786)
#pragma warning(disable : 4541)

#include "XL_local_api.h"
#include <libCCxll\CCxll.h>

#include "ARM_local_interglob.h"
#include "ARM_local_interface.h"
#include "ARM_gp_local_interglob.h"
#include "ARM_xl_trycatch_local.h"

#include <ARM\libarm_frometk\arm_local_parsexml.h>
#include <ARM\libarm_frometk\ARM_local_etoolkit.h>

#include <ARM\libarm_local\ARM_local_class.h>
#include <ARM\libarm_local\ARM_local_glob.h>
#include <ARM\libicm_local\ICM_local_summit.h>

#include <ARM\libarm_local\ARM_local_persistent.h>

#include <ARM\libarm\ARM_result.h>
#include <ARM\local_xlarm\XL_local_xlarm_common.h>
#include <ARM\libarm_local\arm_local_init.h>
#include <ARM\libarm_local\ARM_local_volcrv.h>
#include <GP_Base\gpbase\stringmanip.h>


#include <deque>
#include <stdlib.h>



//to do mettre dans Arm_glob.h
long ARMLOCAL_SetMDFromXML(const CCString& XMLString,
						   const CCString& type,
						   const CCString& ccy,
						   double asof,
						   ARM_result& result,
						   long objId);


__declspec(dllexport) void WINAPI Local_ARM_deconnection_etoolkit()
{
	deconnection_etoolkit();
}

__declspec(dllexport) void WINAPI Local_ARM_etoolkit_getasof()
{
	etoolkit_getasof();
}

__declspec(dllexport) int WINAPI Local_ARM_etoolkit_connecte()
{
	return etoolkit_connecte();
}


long FreeCurCellContent (const CCString & envVal);


__declspec(dllexport) LPXLOPER WINAPI Local_ARM_GetLastDateWarm()
{
	long retCode;
	ARM_result C_result;

	// return
	static XLOPER XL_result;

	//	ARM_BEGIN();
	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();

	// C variable
	// error
	static int error;
	static char* reason = "";

	retCode = ARMLOCAL_GetLastDateWarm (C_result);

	if(retCode == ARM_OK)
	{
		FreeCurCellErr ();
		XL_result.xltype = xltypeStr;
		XL_result.val.str = XL_StrC2StrPascal (C_result.getString());
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_ARM_GetLastDateWarm" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;

}


__declspec(dllexport) LPXLOPER WINAPI Local_ARM_GETINSTRUMENTFROMSUMMIT(LPXLOPER XL_idSummit,
																		LPXLOPER XL_type,
																		LPXLOPER XL_asOf,
																		LPXLOPER XL_ExoticFilter)
{
	long retCode;
	ARM_result C_result;

	// return
	static XLOPER XL_result;

	//	ARM_BEGIN();
	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();

	// C variable
	CCString C_idSummit;

	CCString C_type;
	CCString typeId;

	double C_asOf;
	double C_asOf_default = -1.0;
	
	CCString C_filter;

	// error
	static int error;
	static char* reason = "";

	XL_readStrCell(XL_idSummit,C_idSummit," ARM_ERR: summit id: string expected",C_result);
	XL_readStrCell(XL_type,C_type," ARM_ERR: type: string expected",C_result);
	XL_readNumCellWD(XL_asOf,C_asOf,C_asOf_default," ARM_ERR: asOf: date expected",C_result);
	XL_readStrCellWD(XL_ExoticFilter,C_filter,"ALL"," ARM_ERR: Exotic filter: string expected ie ALL, SPDOPT,...",C_result);

	C_idSummit.toUpper ();
	C_type.toUpper ();

	if((typeId = ARM_ConvTypeDeal (C_type, C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	CCString prevClass;
	
	CCString curClass = typeId;
	// NB :
	// In summit, caption legs may be included in an exotic portfolio.
	// the type ID should then be "EXOTIC"
	// because GP Calculators don't inherit from Security, 
	// if we specify the caption leg in filter, we will build only the caption leg
	if (C_type == "EXOTIC" && C_filter == "ALMCAPTION")
		curClass = LOCAL_GC_CAPTION_CLASS;

	if (C_type == "GLOBALCAP")
	{
		C_type = "IRG";
	}

	CCString stringId = GetLastCurCellEnvValue ();
	
	long objId;

	if(!stringId)
	{
		retCode = ARMLOCAL_GETINSTRUMENTFROMSUMMIT (C_idSummit,
													C_type,
													C_asOf,
													C_filter,
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
			retCode = ARMLOCAL_GETINSTRUMENTFROMSUMMIT (C_idSummit,
														C_type,
														C_asOf,
														C_filter,
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
			retCode = ARMLOCAL_GETINSTRUMENTFROMSUMMIT (C_idSummit,
														C_type,
														C_asOf,
														C_filter,
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
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_ARM_GETINSTRUMENTFROMSUMMIT" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}



__declspec(dllexport) LPXLOPER WINAPI Local_PXL_ARM_GETINSTRUMENTFROMSUMMIT(LPXLOPER XL_idSummit,
																			LPXLOPER XL_type,
																			LPXLOPER XL_asOf,
																			LPXLOPER XL_ExoticFilter)
{
	long retCode;
	ARM_result C_result;

	// return
	static XLOPER XL_result;

	//	ARM_BEGIN();

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();

	// C variable
	CCString C_idSummit;

	CCString C_type;
	CCString typeId;

	double C_asOf;
	double C_asOf_default = -1.0;
	CCString C_filter;

	// error
	static int error;
	static char* reason = "";

	XL_readStrCell(XL_idSummit,C_idSummit," ARM_ERR: summit id: string expected",C_result);
	XL_readStrCell(XL_type,C_type," ARM_ERR: type: string expected",C_result);
	XL_readNumCellWD(XL_asOf,C_asOf,C_asOf_default," ARM_ERR: asOf: date expected",C_result);
	XL_readStrCellWD(XL_ExoticFilter,C_filter,"ALL"," ARM_ERR: Exotic filter: string expected ie ALL, SPDOPT,...",C_result);


	C_idSummit.toUpper ();
	C_type.toUpper ();

	if((typeId = ARM_ConvTypeDeal (C_type, C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	CCString prevClass;
	
	CCString curClass = typeId;
	CCString stringId;
	
	long objId;

	retCode = ARMLOCAL_GETINSTRUMENTFROMSUMMIT (C_idSummit,
												C_type,
												C_asOf,
												C_filter,
												C_result);

	if ( retCode == ARM_OK )
	{
		objId = C_result.getLong ();

		stringId = LocalMakeObjectId (objId, curClass);
	}

	if ( retCode == ARM_OK )
	{
		// FreeCurCellErr ();
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_PXL_ARM_GETINSTRUMENTFROMSUMMIT" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


__declspec(dllexport) LPXLOPER WINAPI Local_ARM_GETINSTRUMENTFROMCALYPSO(LPXLOPER XL_idCalypso,
																		 LPXLOPER XL_type,
																		 LPXLOPER XL_asOf,
																		 LPXLOPER XL_modelType)
{
	long retCode;
	ARM_result C_result;

	// return
	static XLOPER XL_result;

	//	ARM_BEGIN();
	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();

	// C variable
	CCString C_idCalypso;
	CCString C_modelType;
	CCString C_type;
	CCString typeId;

	double C_asOf;
	double C_asOf_default = -1.0;
	
	CCString C_filter;

	// error
	static int error;
	static char* reason = "";


	XL_readStrCell(XL_idCalypso,C_idCalypso," ARM_ERR: calypso id: numeric expected",C_result);
	XL_readStrCell(XL_type,C_type," ARM_ERR: type: string expected",C_result);
	XL_readNumCellWD(XL_asOf,C_asOf,C_asOf_default," ARM_ERR: asOf: date expected",C_result);
	XL_readStrCellWD(XL_modelType,C_modelType,"CRF_HW","ARM_ERR: type: string expected",C_result);
	
	C_type.toUpper ();

	if((typeId = ARM_ConvTypeDeal (C_type, C_result)) == ARM_DEFAULT_ERR)
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
		retCode = ARMLOCAL_GETINSTRUMENTFROMCALYPSO (C_idCalypso,
													 C_type,
													 C_modelType,
													 C_asOf,
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
			retCode = ARMLOCAL_GETINSTRUMENTFROMCALYPSO (C_idCalypso,
														 C_type,
														 C_modelType,
														 C_asOf,
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
			retCode = ARMLOCAL_GETINSTRUMENTFROMCALYPSO (C_idCalypso,
														 C_type,
														 C_modelType,
														 C_asOf,
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
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_ARM_GETINSTRUMENTFROMCALYPSO" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}



__declspec(dllexport) LPXLOPER WINAPI Local_PXL_ARM_GETINSTRUMENTFROMCALYPSO(LPXLOPER XL_idCalypso,
																			 LPXLOPER XL_type,
																			 LPXLOPER XL_asOf,
																			 LPXLOPER XL_modelType)
{
	long retCode;
	ARM_result C_result;

	// return
	static XLOPER XL_result;

	//	ARM_BEGIN();
	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();

	// C variable
	CCString C_idCalypso;

	CCString C_type;
	CCString C_modelType;
	CCString typeId;

	double C_asOf;
	double C_asOf_default = -1.0;
	
	CCString C_filter;

	// error
	static int error;
	static char* reason = "";

	XL_readStrCell(XL_idCalypso,C_idCalypso," ARM_ERR: calypso id: string expected",C_result);
	XL_readStrCell(XL_type,C_type," ARM_ERR: type: string expected",C_result);
	XL_readNumCellWD(XL_asOf,C_asOf,C_asOf_default," ARM_ERR: asOf: date expected",C_result);
	XL_readStrCellWD(XL_modelType,C_modelType,"CRF_HW","ARM_ERR: type: string expected",C_result);

	C_type.toUpper ();

	if((typeId = ARM_ConvTypeDeal (C_type, C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	CCString prevClass;
	
	CCString curClass = typeId;
	CCString stringId;
	
	long objId;

	retCode = ARMLOCAL_GETINSTRUMENTFROMCALYPSO(C_idCalypso,
												C_type,
												C_modelType,
												C_asOf,
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

	
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_ARM_GETINSTRUMENTFROMCALYPSO" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}



__declspec(dllexport) LPXLOPER WINAPI Local_ARM_INITCRF(LPXLOPER XL_crfId,
														LPXLOPER XL_zcCpnId,
														LPXLOPER XL_swoptVol,
														LPXLOPER XL_capVol,
														LPXLOPER XL_rhoCap,
														LPXLOPER XL_nuCap,
														LPXLOPER XL_zcFund,
														LPXLOPER XL_zcCpnBasis,
														LPXLOPER XL_zcFundBasis,
														LPXLOPER XL_fx,
														LPXLOPER XL_isUpdate,
														LPXLOPER XL_betaCap,
														LPXLOPER XL_rhoSwopt,
														LPXLOPER XL_nuSwopt,
														LPXLOPER XL_betaSwopt,
														LPXLOPER XL_modelType,
														LPXLOPER XL_meanReversionId,
                                                        LPXLOPER XL_skewRecalibFlag,
                                                        LPXLOPER XL_SigmaOrAlpha)
{
	long retCode;
	ARM_result C_result;

	static XLOPER XL_result;

	ARM_XL_TRY_BLOCK_BEGIN
	{
		ARM_NOCALCIFWIZ();

		// C variable
		CCString C_crfId;
		CCString C_zcCpnId;
		CCString C_swVolId;
		CCString C_capVolId;

		CCString C_rhoCap;
		long rhoCapId;
		CCString C_nuCap;
		long nuCapId;
		CCString C_betaCap;
		long betaCapId;

		CCString C_rhoSwopt;
		long rhoSwoptId;
		CCString C_nuSwopt;
		long nuSwoptId;
		CCString C_betaSwopt;
		long betaSwoptId;

		CCString C_zcFund;
		long zcFundId;
		CCString C_zcCpnBasis;
		long zcCpnBasisId;
		CCString C_zcFundBasis;
		long zcFundBasisId;
		double C_fx;
		double C_fxDef = 0.0;

		CCString C_isUpdate;
		long isUpdateId;

		CCString C_modelType;

		CCString C_meanReversionId;
		long meanReversionId;

        CCString C_skewRecalibFlag;
		long skewRecalibFlag;

        CCString C_SigmaOrAlpha;


		// error
		static int error;
		static char* reason = "";

		XL_readStrCell(	  XL_crfId,		 C_crfId,		" ARM_ERR: crf id: object expected",		C_result);
		XL_readStrCell(   XL_zcCpnId,	 C_zcCpnId,		" ARM_ERR: zc id: object expected",			C_result);
		XL_readStrCell(   XL_swoptVol,	 C_swVolId,		" ARM_ERR: swopt Vol id: object expected",	C_result);
		XL_readStrCell(   XL_capVol,	 C_capVolId,	" ARM_ERR: crf id: object expected",		C_result);
		XL_readStrCellWD( XL_rhoCap,	 C_rhoCap,		"DEFAULT"," ARM_ERR: rhoCap(Rho): object expected",		C_result);	
		XL_readStrCellWD( XL_nuCap,		 C_nuCap,		"DEFAULT"," ARM_ERR: nuCap(Nu): object expected",		C_result);
		XL_readStrCellWD( XL_betaCap,	 C_betaCap,		"DEFAULT"," ARM_ERR: betaCap(beta): object expected",	C_result);
		XL_readStrCellWD( XL_rhoSwopt,	 C_rhoSwopt,	"DEFAULT"," ARM_ERR: rhoSwopt(Rho): object expected",	C_result);	
		XL_readStrCellWD( XL_nuSwopt,	 C_nuSwopt,		"DEFAULT"," ARM_ERR: nuSwopt(Nu): object expected",		C_result);
		XL_readStrCellWD( XL_betaSwopt,	 C_betaSwopt,	"DEFAULT"," ARM_ERR: betaSwopt(beta): object expected",	C_result);
		XL_readStrCellWD( XL_zcFund,	 C_zcFund,		"DEFAULT"," ARM_ERR: zc funding: object expected",	C_result);
		XL_readStrCellWD( XL_zcCpnBasis, C_zcCpnBasis,	"DEFAULT"," ARM_ERR: zc basis: object expected",	C_result);
		XL_readStrCellWD( XL_zcFundBasis,C_zcFundBasis,	"DEFAULT"," ARM_ERR: zc basis: object expected",	C_result);
		XL_readNumCellWD( XL_fx,		 C_fx,			C_fxDef,  " ARM_ERR: forex: object expected",		C_result);
		XL_readStrCellWD( XL_isUpdate,	 C_isUpdate,	"NO",	  " ARM_ERR: is Update: string expected",	C_result);
		XL_readStrCellWD( XL_modelType,	 C_modelType,	"HWM1F",  " ARM_ERR: Model type: string expected",	C_result);
		XL_readStrCellWD( XL_meanReversionId, C_meanReversionId, "DEFAULT",  
							" ARM_ERR: Mean reversion Id: object expected", C_result);

        XL_readStrCellWD(XL_skewRecalibFlag, C_skewRecalibFlag, "DEFAULT",  
					     " ARM_ERR: Skew Recalib Flag: Y|N|YES|NO expected", C_result);

        XL_readStrCellWD(XL_SigmaOrAlpha,
                         C_SigmaOrAlpha, "S",
                         " ARM_ERR: SABR vol mode: S, SIGMA, A, ALPHA expected", C_result);

        if ( C_skewRecalibFlag == "DEFAULT" )
        {
           skewRecalibFlag = -1; 
        }
        else
        {
           const char* YesOrNo = (const char *) C_skewRecalibFlag;

           if (( C_skewRecalibFlag[0] == 'Y' ) || ( C_skewRecalibFlag[0] == 'y' ))
           {
              skewRecalibFlag = 1L;
           }
           else if (( C_skewRecalibFlag[0] == 'N' ) || ( C_skewRecalibFlag[0] == 'n' ))
           {
              skewRecalibFlag = 0L;
           }
           else
           {
              skewRecalibFlag = -1;
           }
        }

		if ( C_rhoCap == "DEFAULT" )
		{
		   rhoCapId = ARM_NULL_OBJECT;
		}
		else
		{
		   rhoCapId = LocalGetNumObjectId(C_rhoCap);
		}

		if ( C_nuCap == "DEFAULT" )
		{
		   nuCapId = ARM_NULL_OBJECT;
		}
		else
		{
		   nuCapId = LocalGetNumObjectId(C_nuCap);
		}

		if ( C_betaCap == "DEFAULT" )
		{
		   betaCapId = ARM_NULL_OBJECT;
		}
		else
		{
		   betaCapId = LocalGetNumObjectId(C_betaCap);
		}

		if ( C_rhoSwopt == "DEFAULT" )
		{
		   rhoSwoptId = ARM_NULL_OBJECT;
		}
		else
		{
		   rhoSwoptId = LocalGetNumObjectId(C_rhoSwopt);
		}

		if ( C_nuSwopt == "DEFAULT" )
		{
		   nuSwoptId = ARM_NULL_OBJECT;
		}
		else
		{
		   nuSwoptId = LocalGetNumObjectId(C_nuSwopt);
		}

		if ( C_betaSwopt == "DEFAULT" )
		{
		   betaSwoptId = ARM_NULL_OBJECT;
		}
		else
		{
		   betaSwoptId = LocalGetNumObjectId(C_betaSwopt);
		}

		if ( C_zcFund == "DEFAULT" )
		{
		   zcFundId = ARM_NULL_OBJECT;
		}
		else
		{
		   zcFundId = LocalGetNumObjectId(C_zcFund);
		}

		if ( C_zcCpnBasis == "DEFAULT" )
		{
		   zcCpnBasisId = ARM_NULL_OBJECT;
		}
		else
		{
		   zcCpnBasisId = LocalGetNumObjectId(C_zcCpnBasis);
		}

		if ( C_zcFundBasis == "DEFAULT" )
		{
		   zcFundBasisId = ARM_NULL_OBJECT;
		}
		else
		{
		   zcFundBasisId = LocalGetNumObjectId(C_zcFundBasis);
		}

		if((isUpdateId = ARM_ConvYesOrNo (C_isUpdate, C_result)) == ARM_DEFAULT_ERR)
		{
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}

		if ( C_meanReversionId == "DEFAULT" )
		{
			meanReversionId = ARM_NULL_OBJECT;
		}
		else
		{
			meanReversionId = LocalGetNumObjectId(C_meanReversionId);
		}

        long SigmaOrAlpha;

        if (( C_SigmaOrAlpha == "S" )
            ||
            ( C_SigmaOrAlpha == "SIGMA" )
            ||
            ( C_SigmaOrAlpha == "sigma" )
            ||
            ( C_SigmaOrAlpha == "Sigma" )
           )
        {
           SigmaOrAlpha = 1;
        }
        else
        {
           SigmaOrAlpha = 0;
        }

		CCString curClass = LOCAL_GC_CRF_CLASS;

		long objId;

		CCString stringId = GetLastCurCellEnvValue ();
		
		if (!stringId)
		   objId = -1;
		else
		{
			CCString prevClass = LocalGetStringObjectClass (stringId);
			
			objId = LocalGetNumObjectId (stringId);
				
			if ( curClass != prevClass )
			{
				FreeCurCellContent();

				objId = -1;
			}
		}

		retCode = ARMLOCAL_INITCRF (LocalGetNumObjectId(C_crfId),
									LocalGetNumObjectId(C_zcCpnId),
									LocalGetNumObjectId(C_swVolId),
									LocalGetNumObjectId(C_capVolId),
									rhoCapId,
									nuCapId,
									betaCapId,
									rhoSwoptId,
									nuSwoptId,
									betaSwoptId,
									zcFundId,
									zcCpnBasisId,
									zcFundBasisId,
									C_fx,
									isUpdateId,
									C_modelType,
									meanReversionId,
                                    skewRecalibFlag,
                                    SigmaOrAlpha,
									C_result,
									objId);

		if ( retCode == ARM_OK )
		{
			objId = C_result.getLong ();

			LocalSetCurCellEnvValue (curClass, objId); 

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

	ARM_XL_TRY_BLOCK_END

	ARM_XL_CATCH_ARM_EXPT

	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_ARM_INITCRF" )

	return (LPXLOPER)&XL_result;
}




__declspec(dllexport) LPXLOPER WINAPI Local_PXL_ARM_INITCRF(LPXLOPER XL_crfId,
															LPXLOPER XL_zcCpnId,
															LPXLOPER XL_swoptVol,
															LPXLOPER XL_capVol,
															LPXLOPER XL_rhoCap,
															LPXLOPER XL_nuCap,
															LPXLOPER XL_zcFund,
															LPXLOPER XL_zcCpnBasis,
															LPXLOPER XL_zcFundBasis,
															LPXLOPER XL_fx,
															LPXLOPER XL_isUpdate,
															LPXLOPER XL_betaCap,
															LPXLOPER XL_rhoSwopt,
															LPXLOPER XL_nuSwopt,
															LPXLOPER XL_betaSwopt,
															LPXLOPER XL_modelType,
															LPXLOPER XL_meanReversionId,
                                                            LPXLOPER XL_skewRecalibFlag,
                                                            LPXLOPER XL_SigmaOrAlpha)
{
	long retCode;
	ARM_result C_result;

	static XLOPER XL_result;

	ARM_XL_TRY_BLOCK_BEGIN
	{
		ARM_NOCALCIFWIZ();

		// C variable
		CCString C_crfId;
		CCString C_zcCpnId;
		CCString C_swVolId;
		CCString C_capVolId;

		CCString C_rhoCap;
		long rhoCapId;
		CCString C_nuCap;
		long nuCapId;
		CCString C_betaCap;
		long betaCapId;

		CCString C_rhoSwopt;
		long rhoSwoptId;
		CCString C_nuSwopt;
		long nuSwoptId;
		CCString C_betaSwopt;
		long betaSwoptId;

		CCString C_zcFund;
		long zcFundId;
		CCString C_zcCpnBasis;
		long zcCpnBasisId;
		CCString C_zcFundBasis;
		long zcFundBasisId;
		double C_fx;
		double C_fxDef = 0.0;

		CCString C_isUpdate;
		long isUpdateId;

		CCString C_modelType;

		CCString C_meanReversionId;
		long meanReversionId;

        CCString C_skewRecalibFlag;
		long skewRecalibFlag;

        CCString C_SigmaOrAlpha;

		// error
		static int error;
		static char* reason = "";

		XL_readStrCell(	  XL_crfId,		 C_crfId,		" ARM_ERR: crf id: object expected",		C_result);
		XL_readStrCell(   XL_zcCpnId,	 C_zcCpnId,		" ARM_ERR: zc id: object expected",			C_result);
		XL_readStrCell(   XL_swoptVol,	 C_swVolId,		" ARM_ERR: swopt Vol id: object expected",	C_result);
		XL_readStrCell(   XL_capVol,	 C_capVolId,	" ARM_ERR: crf id: object expected",		C_result);
		XL_readStrCellWD( XL_rhoCap,	 C_rhoCap,		"DEFAULT"," ARM_ERR: rhoCap(Rho): object expected",		C_result);	
		XL_readStrCellWD( XL_nuCap,		 C_nuCap,		"DEFAULT"," ARM_ERR: nuCap(Nu): object expected",		C_result);
		XL_readStrCellWD( XL_betaCap,	 C_betaCap,		"DEFAULT"," ARM_ERR: betaCap(beta): object expected",	C_result);
		XL_readStrCellWD( XL_rhoSwopt,	 C_rhoSwopt,	"DEFAULT"," ARM_ERR: rhoSwopt(Rho): object expected",	C_result);	
		XL_readStrCellWD( XL_nuSwopt,	 C_nuSwopt,		"DEFAULT"," ARM_ERR: nuSwopt(Nu): object expected",		C_result);
		XL_readStrCellWD( XL_betaSwopt,	 C_betaSwopt,	"DEFAULT"," ARM_ERR: betaSwopt(beta): object expected",	C_result);
		XL_readStrCellWD( XL_zcFund,	 C_zcFund,		"DEFAULT"," ARM_ERR: zc funding: object expected",	C_result);
		XL_readStrCellWD( XL_zcCpnBasis, C_zcCpnBasis,	"DEFAULT"," ARM_ERR: zc basis: object expected",	C_result);
		XL_readStrCellWD( XL_zcFundBasis,C_zcFundBasis,	"DEFAULT"," ARM_ERR: zc basis: object expected",	C_result);
		XL_readNumCellWD( XL_fx,		 C_fx,			C_fxDef,  " ARM_ERR: forex: object expected",		C_result);
		XL_readStrCellWD( XL_isUpdate,	 C_isUpdate,	"NO",	  " ARM_ERR: is Update: string expected",	C_result);
		XL_readStrCellWD( XL_modelType,	 C_modelType,	"HWM1F",  " ARM_ERR: Model type: string expected",	C_result);
		XL_readStrCellWD( XL_meanReversionId, C_meanReversionId, "DEFAULT",  
							" ARM_ERR: Mean reversion Id: object expected", C_result);

        XL_readStrCellWD(XL_skewRecalibFlag, C_skewRecalibFlag, "DEFAULT",  
					     " ARM_ERR: Skew Recalib Flag: Y|N|YES|NO expected", C_result);

        XL_readStrCellWD(XL_SigmaOrAlpha,
                         C_SigmaOrAlpha, "S",
                         " ARM_ERR: SABR vol mode: S, SIGMA, A, ALPHA expected", C_result);

        if ( C_skewRecalibFlag == "DEFAULT" )
        {
           skewRecalibFlag = -1; 
        }
        else
        {
           const char* YesOrNo = (const char *) C_skewRecalibFlag;

           if (( C_skewRecalibFlag[0] == 'Y' ) || ( C_skewRecalibFlag[0] == 'y' ))
           {
              skewRecalibFlag = 1L;
           }
           else if (( C_skewRecalibFlag[0] == 'N' ) || ( C_skewRecalibFlag[0] == 'n' ))
           {
              skewRecalibFlag = 0L;
           }
           else
           {
              skewRecalibFlag = -1;
           }
        }

		if ( C_rhoCap == "DEFAULT" )
		{
		   rhoCapId = ARM_NULL_OBJECT;
		}
		else
		{
		   rhoCapId = LocalGetNumObjectId(C_rhoCap);
		}

		if ( C_nuCap == "DEFAULT" )
		{
		   nuCapId = ARM_NULL_OBJECT;
		}
		else
		{
		   nuCapId = LocalGetNumObjectId(C_nuCap);
		}

		if ( C_betaCap == "DEFAULT" )
		{
		   betaCapId = ARM_NULL_OBJECT;
		}
		else
		{
		   betaCapId = LocalGetNumObjectId(C_betaCap);
		}

		if ( C_rhoSwopt == "DEFAULT" )
		{
		   rhoSwoptId = ARM_NULL_OBJECT;
		}
		else
		{
		   rhoSwoptId = LocalGetNumObjectId(C_rhoSwopt);
		}

		if ( C_nuSwopt == "DEFAULT" )
		{
		   nuSwoptId = ARM_NULL_OBJECT;
		}
		else
		{
		   nuSwoptId = LocalGetNumObjectId(C_nuSwopt);
		}

		if ( C_betaSwopt == "DEFAULT" )
		{
		   betaSwoptId = ARM_NULL_OBJECT;
		}
		else
		{
		   betaSwoptId = LocalGetNumObjectId(C_betaSwopt);
		}

		if ( C_zcFund == "DEFAULT" )
		{
		   zcFundId = ARM_NULL_OBJECT;
		}
		else
		{
		   zcFundId = LocalGetNumObjectId(C_zcFund);
		}

		if ( C_zcCpnBasis == "DEFAULT" )
		{
		   zcCpnBasisId = ARM_NULL_OBJECT;
		}
		else
		{
		   zcCpnBasisId = LocalGetNumObjectId(C_zcCpnBasis);
		}

		if ( C_zcFundBasis == "DEFAULT" )
		{
		   zcFundBasisId = ARM_NULL_OBJECT;
		}
		else
		{
		   zcFundBasisId = LocalGetNumObjectId(C_zcFundBasis);
		}

		if ((isUpdateId = ARM_ConvYesOrNo (C_isUpdate, C_result)) == ARM_DEFAULT_ERR)
		{
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}

		if ( C_meanReversionId == "DEFAULT" )
		{
			meanReversionId = ARM_NULL_OBJECT;
		}
		else
		{
			meanReversionId = LocalGetNumObjectId(C_meanReversionId);
		}

        long SigmaOrAlpha;

        if (( C_SigmaOrAlpha == "S" )
            ||
            ( C_SigmaOrAlpha == "SIGMA" )
            ||
            ( C_SigmaOrAlpha == "sigma" )
            ||
            ( C_SigmaOrAlpha == "Sigma" )
           )
        {
           SigmaOrAlpha = 1;
        }
        else
        {
           SigmaOrAlpha = 0;
        }

		CCString curClass = LOCAL_GC_CRF_CLASS;

		long objId;

		CCString stringId;

		retCode = ARMLOCAL_INITCRF(LocalGetNumObjectId(C_crfId),
								   LocalGetNumObjectId(C_zcCpnId),
								   LocalGetNumObjectId(C_swVolId),
								   LocalGetNumObjectId(C_capVolId),
								   rhoCapId,
								   nuCapId,
								   betaCapId,
								   rhoSwoptId,
								   nuSwoptId,
								   betaSwoptId,
								   zcFundId,
								   zcCpnBasisId,
								   zcFundBasisId,
								   C_fx,
								   isUpdateId,
								   C_modelType,
								   meanReversionId,
                                   skewRecalibFlag,
                                   SigmaOrAlpha,
								   C_result);

		if ( retCode == ARM_OK )
		{
			objId = C_result.getLong ();

			stringId = LocalMakeObjectId (objId, curClass);
		}

		if ( retCode == ARM_OK )
		{
			// No need FreeCurCellErr();

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

	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_PXL_ARM_INITCRF" )

	return (LPXLOPER)&XL_result;
}



__declspec(dllexport) LPXLOPER WINAPI Local_ARM_INITMATCAP(LPXLOPER XL_matcapId,
														   LPXLOPER XL_zcId,
														   LPXLOPER XL_capVol,
														   LPXLOPER XL_rho,
														   LPXLOPER XL_nu,
														   LPXLOPER XL_beta,
														   LPXLOPER XL_nbpas,
                                                           LPXLOPER XL_SigmaOrAlpha)
{
	long retCode;
	ARM_result C_result;

	static XLOPER XL_result;

	ARM_XL_TRY_BLOCK_BEGIN
	{
		ARM_NOCALCIFWIZ();

		// C variable
		CCString C_matcapId;
		CCString C_zcId;
		CCString C_capVolId;
		CCString C_rho;
		long rhoId;
		CCString C_nu;
		long nuId;
		CCString C_beta;
		long betaId;
		double nbpas;
		double nbpas_default = 10000.0;
		
        CCString C_SigmaOrAlpha;

        // error
		static int error;
		static char* reason = "";

		XL_readStrCell(XL_matcapId,C_matcapId," ARM_ERR: maturity cap id: object expected",C_result);
		XL_readStrCell(XL_zcId,C_zcId," ARM_ERR: zc id: object expected",C_result);
		XL_readStrCell(XL_capVol,C_capVolId," ARM_ERR: crf id: object expected",C_result);
		XL_readStrCellWD(XL_rho,C_rho,"DEFAULT"," ARM_ERR: rho(Rho): object expected",C_result);	
		XL_readStrCellWD(XL_nu,C_nu,"DEFAULT"," ARM_ERR: nu(Nu): object expected",C_result);
		XL_readStrCellWD(XL_beta,C_beta,"DEFAULT"," ARM_ERR: beta(Beta): object expected",C_result);
		XL_readNumCellWD(XL_nbpas,nbpas,nbpas_default," ARM_ERR: nb de pas: numeric expected",C_result);

        XL_readStrCellWD(XL_SigmaOrAlpha,
                         C_SigmaOrAlpha, "S",
                         " ARM_ERR: SABR vol mode: S, SIGMA, A, ALPHA expected", C_result);

		if ( C_rho == "DEFAULT" )
		{
		   rhoId = ARM_NULL_OBJECT;
		}
		else
		{
		   rhoId = LocalGetNumObjectId(C_rho);
		}

		if ( C_nu == "DEFAULT" )
		{
		   nuId = ARM_NULL_OBJECT;
		}
		else
		{
		   nuId = LocalGetNumObjectId(C_nu);
		}

		if ( C_beta == "DEFAULT" )
		{
		   betaId = ARM_NULL_OBJECT;
		}
		else
		{
		   betaId = LocalGetNumObjectId(C_beta);
		}

		CCString prevClass;

        long SigmaOrAlpha;

        if (( C_SigmaOrAlpha == "S" )
            ||
            ( C_SigmaOrAlpha == "SIGMA" )
            ||
            ( C_SigmaOrAlpha == "sigma" )
            ||
            ( C_SigmaOrAlpha == "Sigma" )
           )
        {
           SigmaOrAlpha = 1;
        }
        else
        {
           SigmaOrAlpha = 0;
        }
		
		CCString curClass = LOCAL_GC_MATURITYCAP_CLASS;

		long objId;

        CCString stringId = GetLastCurCellEnvValue();
		
		if (!stringId)
		   objId = -1;
		else
		{
			CCString prevClass = LocalGetStringObjectClass (stringId);
			
			objId = LocalGetNumObjectId (stringId);
				
			if ( curClass != prevClass )
			{
				FreeCurCellContent();

				objId = -1;
			}
		}

		retCode = ARMLOCAL_INITMATCAP(LocalGetNumObjectId(C_matcapId),
									  LocalGetNumObjectId(C_zcId),
									  LocalGetNumObjectId(C_capVolId),
									  rhoId,
									  nuId,
									  betaId,
									  (long) nbpas,
                                      SigmaOrAlpha,
									  C_result,
                                      objId);

		if ( retCode == ARM_OK )
		{
		   objId = C_result.getLong ();

           LocalSetCurCellEnvValue(curClass, objId); 

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

	ARM_XL_TRY_BLOCK_END

	ARM_XL_CATCH_ARM_EXPT

	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_ARM_INITMATCAP" )

	return (LPXLOPER)&XL_result;
}



__declspec(dllexport) LPXLOPER WINAPI Local_PXL_ARM_INITMATCAP(LPXLOPER XL_matcapId,
															   LPXLOPER XL_zcId,
															   LPXLOPER XL_capVol,
															   LPXLOPER XL_rho,
															   LPXLOPER XL_nu,
															   LPXLOPER XL_beta,
															   LPXLOPER XL_nbpas,
                                                               LPXLOPER XL_SigmaOrAlpha)
{
	long retCode;
	ARM_result C_result;

	static XLOPER XL_result;

	ARM_XL_TRY_BLOCK_BEGIN
	{
		ARM_NOCALCIFWIZ();

		// C variable
		CCString C_matcapId;
		CCString C_zcId;
		CCString C_capVolId;
		CCString C_rho;
		long rhoId;
		CCString C_nu;
		long nuId;
		CCString C_beta;
		long betaId;
		double nbpas;
		double nbpas_default = 10000.0;

        CCString C_SigmaOrAlpha;

        // error
		static int error;
		static char* reason = "";

		XL_readStrCell(XL_matcapId,C_matcapId," ARM_ERR: maturity cap id: object expected",C_result);
		XL_readStrCell(XL_zcId,C_zcId," ARM_ERR: zc id: object expected",C_result);
		XL_readStrCell(XL_capVol,C_capVolId," ARM_ERR: crf id: object expected",C_result);
		XL_readStrCellWD(XL_rho,C_rho,"DEFAULT"," ARM_ERR: rho(Rho): object expected",C_result);	
		XL_readStrCellWD(XL_nu,C_nu,"DEFAULT"," ARM_ERR: nu(Nu): object expected",C_result);
		XL_readStrCellWD(XL_beta,C_beta,"DEFAULT"," ARM_ERR: beta(Beta): object expected",C_result);
		XL_readNumCellWD(XL_nbpas,nbpas,nbpas_default," ARM_ERR: nb de pas: numeric expected",C_result);

        XL_readStrCellWD(XL_SigmaOrAlpha,
                         C_SigmaOrAlpha, "S",
                         " ARM_ERR: SABR vol mode: S, SIGMA, A, ALPHA expected", C_result);

		if ( C_rho == "DEFAULT" )
		{
		   rhoId = ARM_NULL_OBJECT;
		}
		else
		{
		   rhoId = LocalGetNumObjectId(C_rho);
		}

		if ( C_nu == "DEFAULT" )
		{
		   nuId = ARM_NULL_OBJECT;
		}
		else
		{
		   nuId = LocalGetNumObjectId(C_nu);
		}

		if ( C_beta == "DEFAULT" )
		{
		   betaId = ARM_NULL_OBJECT;
		}
		else
		{
		   betaId = LocalGetNumObjectId(C_beta);
		}

		CCString prevClass;

        long SigmaOrAlpha;

        if (( C_SigmaOrAlpha == "S" )
            ||
            ( C_SigmaOrAlpha == "SIGMA" )
            ||
            ( C_SigmaOrAlpha == "sigma" )
            ||
            ( C_SigmaOrAlpha == "Sigma" )
           )
        {
           SigmaOrAlpha = 1;
        }
        else
        {
           SigmaOrAlpha = 0;
        }
		
		CCString curClass = LOCAL_GC_MATURITYCAP_CLASS;
		CCString stringId;

		long objId;

		retCode = ARMLOCAL_INITMATCAP(LocalGetNumObjectId(C_matcapId),
									  LocalGetNumObjectId(C_zcId),
									  LocalGetNumObjectId(C_capVolId),
									  rhoId,
									  nuId,
									  betaId,
									  (long) nbpas,
                                      SigmaOrAlpha,
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
	ARM_XL_TRY_BLOCK_END

	ARM_XL_CATCH_ARM_EXPT

	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_PXL_ARM_INITMATCAP" )

	return (LPXLOPER)&XL_result;
}



__declspec(dllexport) LPXLOPER WINAPI Local_ARM_INITTARN(LPXLOPER XL_tarnId,
														 LPXLOPER XL_zcCpnId,
														 LPXLOPER XL_swoptVol,
														 LPXLOPER XL_capVol,
														 LPXLOPER XL_rhoCap,
														 LPXLOPER XL_nuCap,
														 LPXLOPER XL_betaCap,
														 LPXLOPER XL_rhoSwopt,
														 LPXLOPER XL_nuSwopt,
														 LPXLOPER XL_betaSwopt,
														 LPXLOPER XL_zcFund,
														 LPXLOPER XL_zcCpnBasis,
														 LPXLOPER XL_zcFundBasis,
														 LPXLOPER XL_fx,
                                                         LPXLOPER XL_modelType,
                                                         LPXLOPER XL_betaCorrel,
                                                         LPXLOPER XL_hump,
                                                         LPXLOPER XL_SigmaOrAlpha)
{
	long retCode;
	ARM_result C_result;

	static XLOPER XL_result;

	ARM_XL_TRY_BLOCK_BEGIN
	{
		ARM_NOCALCIFWIZ();

		// C variable
		CCString C_tarnId;
		CCString C_zcCpnId;

		CCString C_swVolId;
		CCString C_capVolId;

		CCString C_rhoCap;
		long rhoCapId;
		CCString C_nuCap;
		long nuCapId;
		CCString C_betaCap;
		long betaCapId;

		CCString C_rhoSwopt;
		long rhoSwoptId;
		CCString C_nuSwopt;
		long nuSwoptId;
		CCString C_betaSwopt;
		long betaSwoptId;

		CCString C_zcFund;
		long zcFundId;
		CCString C_zcCpnBasis;
		long zcCpnBasisId;
		CCString C_zcFundBasis;
		long zcFundBasisId;

		double C_fx;
		double C_fxDef = 0.0;

        double C_betaCorrel;
		double C_betaCorrelDef = -10000.0;

        double C_hump;
		double C_humpDef = -10000.0;
		
        CCString C_modelType;

        CCString C_SigmaOrAlpha;

        // error
		static int error;
		static char* reason = "";

		XL_readStrCell(	  XL_tarnId,	 C_tarnId,					" ARM_ERR: Tarn id: object expected",		C_result);
		XL_readStrCell(   XL_zcCpnId,	 C_zcCpnId,					" ARM_ERR: zc id: object expected",		C_result);
		
		XL_readStrCell(   XL_swoptVol,	 C_swVolId,					" ARM_ERR: swopt Vol id: object expected", C_result);
		XL_readStrCell(   XL_capVol,	 C_capVolId,				" ARM_ERR: crf id: object expected",		C_result);

		XL_readStrCellWD( XL_rhoCap,	 C_rhoCap,				    "DEFAULT"," ARM_ERR: rhoCap(Rho): object expected",		C_result);	
		XL_readStrCellWD( XL_nuCap,		 C_nuCap,					"DEFAULT"," ARM_ERR: nuCap(Nu): object expected",		C_result);
		XL_readStrCellWD( XL_betaCap,	 C_betaCap,					"DEFAULT"," ARM_ERR: betaCap(Nu): object expected",		C_result);

		XL_readStrCellWD( XL_rhoSwopt,	 C_rhoSwopt,	  "DEFAULT"," ARM_ERR: rhoSwopt(Rho): object expected",		C_result);	
		XL_readStrCellWD( XL_nuSwopt,	 C_nuSwopt,		  "DEFAULT"," ARM_ERR: nuSwopt(Nu): object expected",		C_result);
		XL_readStrCellWD( XL_betaSwopt,	 C_betaSwopt,	  "DEFAULT"," ARM_ERR: betaSwopt(Nu): object expected",		C_result);
		
		XL_readStrCellWD( XL_zcFund,	 C_zcFund,		"DEFAULT"," ARM_ERR: zc funding: object expected",	C_result);
		XL_readStrCellWD( XL_zcCpnBasis, C_zcCpnBasis,	"DEFAULT"," ARM_ERR: zc basis: object expected",		C_result);
		XL_readStrCellWD( XL_zcFundBasis,C_zcFundBasis,	"DEFAULT"," ARM_ERR: zc basis: object expected",		C_result);
		XL_readNumCellWD( XL_fx,		 C_fx,		   C_fxDef,  " ARM_ERR: forex: object expected", C_result);

        XL_readStrCellWD(XL_modelType, C_modelType, "SFRM2F",
                         " ARM_ERR: SFRM2F or SBGM string expected", C_result);	

        XL_readNumCellWD(XL_betaCorrel,		
                         C_betaCorrel, C_betaCorrelDef,  " ARM_ERR: Beta correl: numeric expected", C_result);

        XL_readNumCellWD(XL_hump,		
                         C_hump, C_humpDef,  " ARM_ERR: hump: numeric expected", C_result);

        XL_readStrCellWD(XL_SigmaOrAlpha,
                         C_SigmaOrAlpha, "S",
                         " ARM_ERR: SABR vol mode: S, SIGMA, A, ALPHA expected", C_result);

		if ( C_rhoCap == "DEFAULT" )
		{
		   rhoCapId = ARM_NULL_OBJECT;
		}
		else
		{
		   rhoCapId = LocalGetNumObjectId(C_rhoCap);
		}

		if ( C_nuCap == "DEFAULT" )
		{
		   nuCapId = ARM_NULL_OBJECT;
		}
		else
		{
		   nuCapId = LocalGetNumObjectId(C_nuCap);
		}

		if ( C_betaCap == "DEFAULT" )
		{
		   betaCapId = ARM_NULL_OBJECT;
		}
		else
		{
		   betaCapId = LocalGetNumObjectId(C_betaCap);
		}

		if ( C_rhoSwopt == "DEFAULT" )
		{
		   rhoSwoptId = ARM_NULL_OBJECT;
		}
		else
		{
		   rhoSwoptId = LocalGetNumObjectId(C_rhoSwopt);
		}

		if ( C_nuSwopt == "DEFAULT" )
		{
		   nuSwoptId = ARM_NULL_OBJECT;
		}
		else
		{
		   nuSwoptId = LocalGetNumObjectId(C_nuSwopt);
		}

		if ( C_betaSwopt == "DEFAULT" )
		{
		   betaSwoptId = ARM_NULL_OBJECT;
		}
		else
		{
		   betaSwoptId = LocalGetNumObjectId(C_betaSwopt);
		}

		if ( C_zcFund == "DEFAULT" )
		{
		   zcFundId = ARM_NULL_OBJECT;
		}
		else
		{
		   zcFundId = LocalGetNumObjectId(C_zcFund);
		}

		if ( C_zcCpnBasis == "DEFAULT" )
		{
		   zcCpnBasisId = ARM_NULL_OBJECT;
		}
		else
		{
		   zcCpnBasisId = LocalGetNumObjectId(C_zcCpnBasis);
		}

		if ( C_zcFundBasis == "DEFAULT" )
		{
		   zcFundBasisId = ARM_NULL_OBJECT;
		}
		else
		{
		   zcFundBasisId = LocalGetNumObjectId(C_zcFundBasis);
		}

        long SigmaOrAlpha;

        if (( C_SigmaOrAlpha == "S" )
            ||
            ( C_SigmaOrAlpha == "SIGMA" )
            ||
            ( C_SigmaOrAlpha == "sigma" )
            ||
            ( C_SigmaOrAlpha == "Sigma" )
           )
        {
           SigmaOrAlpha = 1;
        }
        else
        {
           SigmaOrAlpha = 0;
        }
		
        CCString curClass = LOCAL_GC_TARN_CLASS;

		long objId;

		CCString stringId = GetLastCurCellEnvValue ();
		
		if (!stringId)
		   objId = -1;
		else
		{
		   CCString prevClass = LocalGetStringObjectClass (stringId);
			
		   objId = LocalGetNumObjectId (stringId);
				
		   if ( curClass != prevClass)
           {
			  FreeCurCellContent();

			  objId = -1;
           }
		}
 
		retCode = ARMLOCAL_INITTARN(LocalGetNumObjectId(C_tarnId),
									LocalGetNumObjectId(C_zcCpnId),
									LocalGetNumObjectId(C_swVolId),
									LocalGetNumObjectId(C_capVolId),
									rhoCapId,
									nuCapId,
									betaCapId,
									rhoSwoptId,
									nuSwoptId,
									betaSwoptId,
									zcFundId,
									zcCpnBasisId,
									zcFundBasisId,
									C_fx,
                                    C_modelType,
                                    C_betaCorrel,
                                    C_hump,
                                    SigmaOrAlpha,
									C_result,
									objId);

		if ( retCode == ARM_OK )
		{
			objId = C_result.getLong ();

			LocalSetCurCellEnvValue (curClass, objId); 

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
	}
	ARM_XL_TRY_BLOCK_END

	ARM_XL_CATCH_ARM_EXPT

	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_ARM_INITTARN" )

	return (LPXLOPER)&XL_result;
}



__declspec(dllexport) LPXLOPER WINAPI Local_PXL_ARM_INITTARN(LPXLOPER XL_tarnId,
															 LPXLOPER XL_zcCpnId,
															 LPXLOPER XL_swoptVol,
															 LPXLOPER XL_capVol,
															 LPXLOPER XL_rhoCap,
															 LPXLOPER XL_nuCap,
															 LPXLOPER XL_betaCap,
															 LPXLOPER XL_rhoSwopt,
															 LPXLOPER XL_nuSwopt,
															 LPXLOPER XL_betaSwopt,
															 LPXLOPER XL_zcFund,
															 LPXLOPER XL_zcCpnBasis,
															 LPXLOPER XL_zcFundBasis,
															 LPXLOPER XL_fx,
                                                             LPXLOPER XL_modelType,
                                                             LPXLOPER XL_betaCorrel,
                                                             LPXLOPER XL_hump,
                                                             LPXLOPER XL_SigmaOrAlpha)
{
	long retCode;
	ARM_result C_result;

	static XLOPER XL_result;

	ARM_XL_TRY_BLOCK_BEGIN
	{
		ARM_NOCALCIFWIZ();

		// C variable
		CCString C_tarnId;
		CCString C_zcCpnId;

		CCString C_swVolId;
		CCString C_capVolId;

		CCString C_rhoCap;
		long rhoCapId;
		CCString C_nuCap;
		long nuCapId;
		CCString C_betaCap;
		long betaCapId;

		CCString C_rhoSwopt;
		long rhoSwoptId;
		CCString C_nuSwopt;
		long nuSwoptId;
		CCString C_betaSwopt;
		long betaSwoptId;

		CCString C_zcFund;
		long zcFundId;
		CCString C_zcCpnBasis;
		long zcCpnBasisId;
		CCString C_zcFundBasis;
		long zcFundBasisId;

		double C_fx;
		double C_fxDef = 0.0;

        double C_betaCorrel;
		double C_betaCorrelDef = -10000.0;

        double C_hump;
		double C_humpDef = -10000.0;

        CCString C_modelType;

        CCString C_SigmaOrAlpha;


		// error
		static int error;
		static char* reason = "";

		XL_readStrCell(	  XL_tarnId,	 C_tarnId,					" ARM_ERR: Tarn id: object expected",		C_result);
		XL_readStrCell(   XL_zcCpnId,	 C_zcCpnId,					" ARM_ERR: zc id: object expected",		C_result);
		
		XL_readStrCell(   XL_swoptVol,	 C_swVolId,					" ARM_ERR: swopt Vol id: object expected", C_result);
		XL_readStrCell(   XL_capVol,	 C_capVolId,				" ARM_ERR: cap Vol id: object expected",		C_result);

		XL_readStrCellWD( XL_rhoCap,	 C_rhoCap,				    "DEFAULT"," ARM_ERR: rhoCap(Rho): object expected",		C_result);	
		XL_readStrCellWD( XL_nuCap,		 C_nuCap,					"DEFAULT"," ARM_ERR: nuCap(Nu): object expected",		C_result);
		XL_readStrCellWD( XL_betaCap,	 C_betaCap,					"DEFAULT"," ARM_ERR: betaCap(Nu): object expected",		C_result);

		XL_readStrCellWD( XL_rhoSwopt,	 C_rhoSwopt,	  "DEFAULT"," ARM_ERR: rhoSwopt(Rho): object expected",		C_result);	
		XL_readStrCellWD( XL_nuSwopt,	 C_nuSwopt,		  "DEFAULT"," ARM_ERR: nuSwopt(Nu): object expected",		C_result);
		XL_readStrCellWD( XL_betaSwopt,	 C_betaSwopt,	  "DEFAULT"," ARM_ERR: betaSwopt(Nu): object expected",		C_result);
		
		XL_readStrCellWD( XL_zcFund,	 C_zcFund,		"DEFAULT"," ARM_ERR: zc funding: object expected",	C_result);
		XL_readStrCellWD( XL_zcCpnBasis, C_zcCpnBasis,	"DEFAULT"," ARM_ERR: zc basis: object expected",		C_result);
		XL_readStrCellWD( XL_zcFundBasis,C_zcFundBasis,	"DEFAULT"," ARM_ERR: zc basis: object expected",		C_result);
		XL_readNumCellWD( XL_fx,		 C_fx,		   C_fxDef,  " ARM_ERR: forex: object expected",		C_result);

        XL_readStrCellWD(XL_modelType, C_modelType, "SFRM2F",
                         " ARM_ERR: SFRM2F or SBGM string expected", C_result);

        XL_readNumCellWD(XL_betaCorrel,		
                         C_betaCorrel, C_betaCorrelDef,  " ARM_ERR: Beta correl: numeric expected", C_result);

        XL_readNumCellWD(XL_hump,		
                         C_hump, C_humpDef,  " ARM_ERR: hump: numeric expected", C_result);

        XL_readStrCellWD(XL_SigmaOrAlpha,
                         C_SigmaOrAlpha, "S",
                         " ARM_ERR: SABR vol mode: S, SIGMA, A, ALPHA expected", C_result);

		if ( C_rhoCap == "DEFAULT" )
		{
		   rhoCapId = ARM_NULL_OBJECT;
		}
		else
		{
		   rhoCapId = LocalGetNumObjectId(C_rhoCap);
		}

		if ( C_nuCap == "DEFAULT" )
		{
		   nuCapId = ARM_NULL_OBJECT;
		}
		else
		{
		   nuCapId = LocalGetNumObjectId(C_nuCap);
		}

		if ( C_betaCap == "DEFAULT" )
		{
		   betaCapId = ARM_NULL_OBJECT;
		}
		else
		{
		   betaCapId = LocalGetNumObjectId(C_betaCap);
		}

		if ( C_rhoSwopt == "DEFAULT" )
		{
		   rhoSwoptId = ARM_NULL_OBJECT;
		}
		else
		{
		   rhoSwoptId = LocalGetNumObjectId(C_rhoSwopt);
		}

		if ( C_nuSwopt == "DEFAULT" )
		{
		   nuSwoptId = ARM_NULL_OBJECT;
		}
		else
		{
		   nuSwoptId = LocalGetNumObjectId(C_nuSwopt);
		}

		if ( C_betaSwopt == "DEFAULT" )
		{
		   betaSwoptId = ARM_NULL_OBJECT;
		}
		else
		{
		   betaSwoptId = LocalGetNumObjectId(C_betaSwopt);
		}

		if ( C_zcFund == "DEFAULT" )
		{
		   zcFundId = ARM_NULL_OBJECT;
		}
		else
		{
		   zcFundId = LocalGetNumObjectId(C_zcFund);
		}

		if ( C_zcCpnBasis == "DEFAULT" )
		{
		   zcCpnBasisId = ARM_NULL_OBJECT;
		}
		else
		{
		   zcCpnBasisId = LocalGetNumObjectId(C_zcCpnBasis);
		}

		if ( C_zcFundBasis == "DEFAULT" )
		{
		   zcFundBasisId = ARM_NULL_OBJECT;
		}
		else
		{
		   zcFundBasisId = LocalGetNumObjectId(C_zcFundBasis);
		}

        long SigmaOrAlpha;

        if (( C_SigmaOrAlpha == "S" )
            ||
            ( C_SigmaOrAlpha == "SIGMA" )
            ||
            ( C_SigmaOrAlpha == "sigma" )
            ||
            ( C_SigmaOrAlpha == "Sigma" )
           )
        {
           SigmaOrAlpha = 1;
        }
        else
        {
           SigmaOrAlpha = 0;
        }

		CCString curClass = LOCAL_GC_TARN_CLASS;

		long objId;

		CCString stringId;

		retCode = ARMLOCAL_INITTARN(LocalGetNumObjectId(C_tarnId),
									LocalGetNumObjectId(C_zcCpnId),
									LocalGetNumObjectId(C_swVolId),
									LocalGetNumObjectId(C_capVolId),
									rhoCapId,
									nuCapId,
									betaCapId,
									rhoSwoptId,
									nuSwoptId,
									betaSwoptId,
									zcFundId,
									zcCpnBasisId,
									zcFundBasisId,
									C_fx,
                                    C_modelType,
                                    C_betaCorrel,
                                    C_hump,
                                    SigmaOrAlpha,
									C_result);

		if ( retCode == ARM_OK )
		{
			objId = C_result.getLong ();

			stringId = LocalMakeObjectId (objId, curClass);
		}

		if ( retCode == ARM_OK )
		{
			// FreeCurCellErr();

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

	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_PXL_ARM_INITTARN" )

	return (LPXLOPER)&XL_result;
}



__declspec(dllexport) LPXLOPER WINAPI Local_ARM_INITCSB(LPXLOPER XL_csbId,
														LPXLOPER XL_zcCpnId,
														LPXLOPER XL_swoptVol,
														LPXLOPER XL_capVol,
														LPXLOPER XL_rhoCap,
														LPXLOPER XL_nuCap,
														LPXLOPER XL_betaCap,
														LPXLOPER XL_rhoSwopt,
														LPXLOPER XL_nuSwopt,
														LPXLOPER XL_betaSwopt,
														LPXLOPER XL_hump,
														LPXLOPER XL_betaCorrel,
														LPXLOPER XL_reCorrel,
														LPXLOPER XL_SigmaOrAlpha)
{
	long retCode;
	ARM_result C_result;

	static XLOPER XL_result;

	ARM_XL_TRY_BLOCK_BEGIN
	{
		ARM_NOCALCIFWIZ();

		// C variable
		CCString C_csbId;
		CCString C_zcCpnId;

		CCString C_swVolId;
		CCString C_capVolId;

		CCString C_rhoCap;
		long rhoCapId;
		CCString C_nuCap;
		long nuCapId;
		CCString C_betaCap;
		long betaCapId;

		CCString C_rhoSwopt;
		long rhoSwoptId;
		CCString C_nuSwopt;
		long nuSwoptId;
		CCString C_betaSwopt;
		long betaSwoptId;

        double C_hump;
		double C_humpDef = -10000.0;
        double C_betaCorrel;
		double C_betaCorrelDef = -10000.0;
        double C_reCorrel;
		double C_reCorrelDef = -10000.0;
		
        CCString C_modelType;

        CCString C_SigmaOrAlpha;

        // error
		static int error;
		static char* reason = "";

		XL_readStrCell(	  XL_csbId,		 C_csbId,					" ARM_ERR: CSB id: object expected",		C_result);
		XL_readStrCell(   XL_zcCpnId,	 C_zcCpnId,					" ARM_ERR: zc id: object expected",		C_result);
		
		XL_readStrCell(   XL_swoptVol,	 C_swVolId,					" ARM_ERR: swopt Vol id: object expected", C_result);
		XL_readStrCell(   XL_capVol,	 C_capVolId,				" ARM_ERR: cap Vol id: object expected",		C_result);

		XL_readStrCellWD( XL_rhoCap,	 C_rhoCap,				    "DEFAULT"," ARM_ERR: rhoCap(Rho): object expected",		C_result);	
		XL_readStrCellWD( XL_nuCap,		 C_nuCap,					"DEFAULT"," ARM_ERR: nuCap(Nu): object expected",		C_result);
		XL_readStrCellWD( XL_betaCap,	 C_betaCap,					"DEFAULT"," ARM_ERR: betaCap(Nu): object expected",		C_result);

		XL_readStrCellWD( XL_rhoSwopt,	 C_rhoSwopt,	  "DEFAULT"," ARM_ERR: rhoSwopt(Rho): object expected",		C_result);	
		XL_readStrCellWD( XL_nuSwopt,	 C_nuSwopt,		  "DEFAULT"," ARM_ERR: nuSwopt(Nu): object expected",		C_result);
		XL_readStrCellWD( XL_betaSwopt,	 C_betaSwopt,	  "DEFAULT"," ARM_ERR: betaSwopt(Nu): object expected",		C_result);
		
        XL_readNumCellWD(XL_hump,		
                         C_hump, C_humpDef,  " ARM_ERR: hump: numeric expected", C_result);
        XL_readNumCellWD(XL_betaCorrel,		
                         C_betaCorrel, C_betaCorrelDef,  " ARM_ERR: Beta correl: numeric expected", C_result);
        XL_readNumCellWD(XL_reCorrel,		
                         C_reCorrel, C_reCorrelDef,  " ARM_ERR: reCorrel: numeric expected", C_result);
        XL_readStrCellWD(XL_SigmaOrAlpha,
                         C_SigmaOrAlpha, "S",
                         " ARM_ERR: SABR vol mode: S, SIGMA, A, ALPHA expected", C_result);


		if ( C_rhoCap == "DEFAULT" )
		{
		   rhoCapId = ARM_NULL_OBJECT;
		}
		else
		{
		   rhoCapId = LocalGetNumObjectId(C_rhoCap);
		}

		if ( C_nuCap == "DEFAULT" )
		{
		   nuCapId = ARM_NULL_OBJECT;
		}
		else
		{
		   nuCapId = LocalGetNumObjectId(C_nuCap);
		}

		if ( C_betaCap == "DEFAULT" )
		{
		   betaCapId = ARM_NULL_OBJECT;
		}
		else
		{
		   betaCapId = LocalGetNumObjectId(C_betaCap);
		}

		if ( C_rhoSwopt == "DEFAULT" )
		{
		   rhoSwoptId = ARM_NULL_OBJECT;
		}
		else
		{
		   rhoSwoptId = LocalGetNumObjectId(C_rhoSwopt);
		}

		if ( C_nuSwopt == "DEFAULT" )
		{
		   nuSwoptId = ARM_NULL_OBJECT;
		}
		else
		{
		   nuSwoptId = LocalGetNumObjectId(C_nuSwopt);
		}

		if ( C_betaSwopt == "DEFAULT" )
		{
		   betaSwoptId = ARM_NULL_OBJECT;
		}
		else
		{
		   betaSwoptId = LocalGetNumObjectId(C_betaSwopt);
		}

        long SigmaOrAlpha;

        if (( C_SigmaOrAlpha == "S" )
            ||
            ( C_SigmaOrAlpha == "SIGMA" )
            ||
            ( C_SigmaOrAlpha == "sigma" )
            ||
            ( C_SigmaOrAlpha == "Sigma" )
           )
        {
           SigmaOrAlpha = 1;
        }
        else
        {
           SigmaOrAlpha = 0;
        }

		CCString curClass = LOCAL_GC_CALLABLE_SB_CLASS;

		long objId;
		CCString stringId = GetLastCurCellEnvValue ();
		
		if(!stringId)
			objId = -1;
		else
		{
			CCString prevClass = LocalGetStringObjectClass (stringId);
			
			objId = LocalGetNumObjectId (stringId);
				
			if(curClass != prevClass)
			{
				FreeCurCellContent ();
				objId = -1;
			}
		}
 
		retCode = ARMLOCAL_INITCSB(LocalGetNumObjectId(C_csbId),
								   LocalGetNumObjectId(C_zcCpnId),
								   LocalGetNumObjectId(C_swVolId),
								   LocalGetNumObjectId(C_capVolId),
								   rhoCapId,
								   nuCapId,
								   betaCapId,
								   rhoSwoptId,
								   nuSwoptId,
								   betaSwoptId,
								   C_hump,
								   C_betaCorrel,
								   C_reCorrel,
								   SigmaOrAlpha,
								   C_result,
								   objId);

		if ( retCode == ARM_OK )
		{
			objId = C_result.getLong ();

			LocalSetCurCellEnvValue (curClass, objId); 

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
	}
	ARM_XL_TRY_BLOCK_END

	ARM_XL_CATCH_ARM_EXPT

	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_ARM_INITCSB" )

	return (LPXLOPER)&XL_result;
}



__declspec(dllexport) LPXLOPER WINAPI Local_PXL_ARM_INITCSB(LPXLOPER XL_csbId,
															LPXLOPER XL_zcCpnId,
															LPXLOPER XL_swoptVol,
															LPXLOPER XL_capVol,
															LPXLOPER XL_rhoCap,
															LPXLOPER XL_nuCap,
															LPXLOPER XL_betaCap,
															LPXLOPER XL_rhoSwopt,
															LPXLOPER XL_nuSwopt,
															LPXLOPER XL_betaSwopt,
															LPXLOPER XL_hump,
															LPXLOPER XL_betaCorrel,
															LPXLOPER XL_reCorrel,
															LPXLOPER XL_SigmaOrAlpha)
{
	long retCode;
	ARM_result C_result;

	static XLOPER XL_result;

	ARM_XL_TRY_BLOCK_BEGIN
	{
		ARM_NOCALCIFWIZ();

		// C variable
		CCString C_csbId;
		CCString C_zcCpnId;

		CCString C_swVolId;
		CCString C_capVolId;

		CCString C_rhoCap;
		long rhoCapId;
		CCString C_nuCap;
		long nuCapId;
		CCString C_betaCap;
		long betaCapId;

		CCString C_rhoSwopt;
		long rhoSwoptId;
		CCString C_nuSwopt;
		long nuSwoptId;
		CCString C_betaSwopt;
		long betaSwoptId;

        double C_hump;
		double C_humpDef = -10000.0;
        double C_betaCorrel;
		double C_betaCorrelDef = -10000.0;
        double C_reCorrel;
		double C_reCorrelDef = -10000.0;
		
        CCString C_modelType;

        CCString C_SigmaOrAlpha;

        // error
		static int error;
		static char* reason = "";

		XL_readStrCell(	  XL_csbId,		 C_csbId,					" ARM_ERR: CSB id: object expected",		C_result);
		XL_readStrCell(   XL_zcCpnId,	 C_zcCpnId,					" ARM_ERR: zc id: object expected",		C_result);
		
		XL_readStrCell(   XL_swoptVol,	 C_swVolId,					" ARM_ERR: swopt Vol id: object expected", C_result);
		XL_readStrCell(   XL_capVol,	 C_capVolId,				" ARM_ERR: cap Vol id: object expected",		C_result);

		XL_readStrCellWD( XL_rhoCap,	 C_rhoCap,				    "DEFAULT"," ARM_ERR: rhoCap(Rho): object expected",		C_result);	
		XL_readStrCellWD( XL_nuCap,		 C_nuCap,					"DEFAULT"," ARM_ERR: nuCap(Nu): object expected",		C_result);
		XL_readStrCellWD( XL_betaCap,	 C_betaCap,					"DEFAULT"," ARM_ERR: betaCap(Nu): object expected",		C_result);

		XL_readStrCellWD( XL_rhoSwopt,	 C_rhoSwopt,	  "DEFAULT"," ARM_ERR: rhoSwopt(Rho): object expected",		C_result);	
		XL_readStrCellWD( XL_nuSwopt,	 C_nuSwopt,		  "DEFAULT"," ARM_ERR: nuSwopt(Nu): object expected",		C_result);
		XL_readStrCellWD( XL_betaSwopt,	 C_betaSwopt,	  "DEFAULT"," ARM_ERR: betaSwopt(Nu): object expected",		C_result);
		
        XL_readNumCellWD(XL_hump,		
                         C_hump, C_humpDef,  " ARM_ERR: hump: numeric expected", C_result);
        XL_readNumCellWD(XL_betaCorrel,		
                         C_betaCorrel, C_betaCorrelDef,  " ARM_ERR: Beta correl: numeric expected", C_result);
        XL_readNumCellWD(XL_reCorrel,		
                         C_reCorrel, C_reCorrelDef,  " ARM_ERR: reCorrel: numeric expected", C_result);
        XL_readStrCellWD(XL_SigmaOrAlpha,
                         C_SigmaOrAlpha, "S",
                         " ARM_ERR: SABR vol mode: S, SIGMA, A, ALPHA expected", C_result);


		if ( C_rhoCap == "DEFAULT" )
		{
		   rhoCapId = ARM_NULL_OBJECT;
		}
		else
		{
		   rhoCapId = LocalGetNumObjectId(C_rhoCap);
		}

		if ( C_nuCap == "DEFAULT" )
		{
		   nuCapId = ARM_NULL_OBJECT;
		}
		else
		{
		   nuCapId = LocalGetNumObjectId(C_nuCap);
		}

		if ( C_betaCap == "DEFAULT" )
		{
		   betaCapId = ARM_NULL_OBJECT;
		}
		else
		{
		   betaCapId = LocalGetNumObjectId(C_betaCap);
		}

		if ( C_rhoSwopt == "DEFAULT" )
		{
		   rhoSwoptId = ARM_NULL_OBJECT;
		}
		else
		{
		   rhoSwoptId = LocalGetNumObjectId(C_rhoSwopt);
		}

		if ( C_nuSwopt == "DEFAULT" )
		{
		   nuSwoptId = ARM_NULL_OBJECT;
		}
		else
		{
		   nuSwoptId = LocalGetNumObjectId(C_nuSwopt);
		}

		if ( C_betaSwopt == "DEFAULT" )
		{
		   betaSwoptId = ARM_NULL_OBJECT;
		}
		else
		{
		   betaSwoptId = LocalGetNumObjectId(C_betaSwopt);
		}

        long SigmaOrAlpha;

        if (( C_SigmaOrAlpha == "S" )
            ||
            ( C_SigmaOrAlpha == "SIGMA" )
            ||
            ( C_SigmaOrAlpha == "sigma" )
            ||
            ( C_SigmaOrAlpha == "Sigma" )
           )
        {
           SigmaOrAlpha = 1;
        }
        else
        {
           SigmaOrAlpha = 0;
        }

		CCString curClass = LOCAL_GC_CALLABLE_SB_CLASS;

		long objId;
		CCString stringId;
 
		retCode = ARMLOCAL_INITCSB(LocalGetNumObjectId(C_csbId),
								   LocalGetNumObjectId(C_zcCpnId),
								   LocalGetNumObjectId(C_swVolId),
								   LocalGetNumObjectId(C_capVolId),
								   rhoCapId,
								   nuCapId,
								   betaCapId,
								   rhoSwoptId,
								   nuSwoptId,
								   betaSwoptId,
								   C_hump,
								   C_betaCorrel,
								   C_reCorrel,
								   SigmaOrAlpha,
								   C_result);

		if ( retCode == ARM_OK )
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
	ARM_XL_TRY_BLOCK_END

	ARM_XL_CATCH_ARM_EXPT

	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_PXL_ARM_INITCSB" )

	return (LPXLOPER)&XL_result;
}


__declspec(dllexport) LPXLOPER WINAPI Local_ARM_INITBERMUDASWAPTION(LPXLOPER XL_bsId,
																	LPXLOPER XL_mktDataManager,
																	LPXLOPER XL_calibParams,
																	LPXLOPER XL_controlVariates)
{
	long retCode;
	ARM_result C_result;

	static XLOPER XL_result;

	ARM_XL_TRY_BLOCK_BEGIN
	{
		ARM_NOCALCIFWIZ();

		// C variable
		CCString C_bsId;

		// error
		static int error;
		static char* reason = "";

		//Bermuda Swaption Calculator to update
		XL_readStrCell(XL_bsId, C_bsId, " ARM_ERR: Bermuda Swaption id: object expected", C_result);
		
		//MKTDATAMANAGER: persistent object
		VECTOR<CCString> mktDataManager;
		XL_readStrVector(XL_mktDataManager,mktDataManager," ARM_ERR: Market datas: array of string expected",DOUBLE_TYPE,C_result);
		if(mktDataManager.size()<4)
		{
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}

		long mktDataManagerId = LocalGetNumObjectId (mktDataManager[0]);
		vector< string > mdmKeys(mktDataManager.size()-1);
		for(int i=1;i<mktDataManager.size();++i)
		mdmKeys[i-1]=CCSTringToSTLString(mktDataManager[i]);
		
		//CONTROLVARIATES
		VECTOR<CCString> controlVariates;
		VECTOR<CCString> controlVariatesDef(0);
		XL_readStrVectorWD(XL_controlVariates, controlVariates,controlVariatesDef, " ARM_ERR: control variates: array of string expected",DOUBLE_TYPE,C_result);
		int ctrlSize = controlVariates.size();
		vector<int>* ctrlVariateVec = new vector<int>(ctrlSize);
		vector<double>* ctrlVariatePrices = new vector<double>(ctrlSize);
		CCString cCurrCtrl;
		string sCurrCtrl;
		for (i=0; i<ctrlSize;i++)
		{
			cCurrCtrl = controlVariates[i];
			sCurrCtrl = CCSTringToSTLString(cCurrCtrl);
			(*ctrlVariateVec)[i] = atoi(sCurrCtrl.c_str());
			(*ctrlVariatePrices)[i] = 0;
		}

		//CALIB PARAMS
		//Model Parameters
		VECTOR<CCString> calibParamsXl;
		XL_readStrVector(XL_calibParams, calibParamsXl," ARM_ERR: Calib Parameters: array expected",DOUBLE_TYPE,C_result);
		if(calibParamsXl.size()<25)
		{
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}
		
		//Model Type
		CCString modelTypeXl = calibParamsXl[0];
		string smodelTypeXl = CCSTringToSTLString(modelTypeXl);
		string smodelType(ARM::stringGetUpper(smodelTypeXl));
		int modelType;
		if(smodelType == "SFRM1F")
		{
			modelType = 0;
		}
		else if (smodelType == "SFRM2F")
		{
			modelType = 1;
		}
		else if (smodelType == "QGM1F")
		{
			modelType = 2;
		}
		else if (smodelType == "HW1F")
		{
			modelType = 3;
		}
		else
			throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"only SFRM1F, SFRM2F, QGM1F or H&W are Valid ");


		vector<string> modelParams(7);

		CCString initVolXl = calibParamsXl[1];
		string sinitVol = CCSTringToSTLString(initVolXl);
		modelParams[0] = sinitVol;

		CCString initMrsXl = calibParamsXl[2];
		string sinitMrs = CCSTringToSTLString(initMrsXl);
		modelParams[1] = sinitMrs;

		CCString initBetaXl = calibParamsXl[3];
		string sinitBeta = CCSTringToSTLString(initBetaXl);
		modelParams[2] = sinitBeta;

		CCString initThetaXl = calibParamsXl[4];
		string sinitTheta = CCSTringToSTLString(initThetaXl);
		modelParams[3] = sinitTheta;

		CCString initSkewXl = calibParamsXl[5];
		string sinitSkew = CCSTringToSTLString(initSkewXl);
		modelParams[4] = sinitSkew;

		CCString initMrsLXl = calibParamsXl[6];
		string sinitMrsL = CCSTringToSTLString(initMrsLXl);
		modelParams[5] = sinitMrsL;

		CCString initMrsUXl = calibParamsXl[7];
		string sinitMrsU = CCSTringToSTLString(initMrsUXl);
		modelParams[6] = sinitMrsU;

		//Calib MRS/Beta ? (Y/N
		CCString calibMrsXl  = calibParamsXl[8];
		CCString calibBetaXl = calibParamsXl[9];
		string calibMrs  = CCSTringToSTLString(calibMrsXl);
		string calibBeta = CCSTringToSTLString(calibBetaXl);
		
		bool mrsCalibFlag = false;
		bool betaCalibFlag = false;
		if (calibMrs == "Y")
		{
			mrsCalibFlag = true;
		}
		if (calibBeta == "Y")
		{
			betaCalibFlag = true;
		}

		//Num Method Type
		CCString methodTypeXl = calibParamsXl[10];
		string smethodTypeXl = CCSTringToSTLString(methodTypeXl);
		string smethodType(ARM::stringGetUpper(smethodTypeXl));
		int numMethodType;
		if(smethodType == "TREE")
		{
			numMethodType = 1;
		}
		else if (smethodType == "MONTECARLO")
		{
			numMethodType = 0;
		}
		else if (smethodType == "PDE")
		{
			numMethodType = 2;
		}
		else
			throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"only TREE or MONTECARLO are valid ");

		//Num Amc Iter
		CCString amcIterXl = calibParamsXl[11];
		string samcIter = CCSTringToSTLString(amcIterXl);
		int amcIter = atoi(samcIter.c_str());
		
		//Num Mc Iter
		CCString mcIterXl = calibParamsXl[12];
		string smcIter = CCSTringToSTLString(mcIterXl);
		int mcIter = atoi(smcIter.c_str());

		//Num Max Bucket Size
		CCString maxBucketSizeXl = calibParamsXl[13];
		string smaxBucket = CCSTringToSTLString(maxBucketSizeXl );
		int maxBucketSize = atoi(smaxBucket.c_str());

		//Num Tree Steps
		CCString treeStepsXl = calibParamsXl[14];
		string streeSteps = CCSTringToSTLString(treeStepsXl);
		int treeSteps = atoi(streeSteps.c_str());

		//Portfolio Mode
		vector<int> portfolioMode(5);

		//Calib Mode
		CCString ptflModeXl = calibParamsXl[15];
		string sptflModeXl = CCSTringToSTLString(ptflModeXl);
		string sptflMode(ARM::stringGetUpper(sptflModeXl));
		int ptflMode;
		if(sptflMode == "SUMMIT")
		{
			ptflMode = 0;
		}
		else if (sptflMode == "MANUAL")
		{
			ptflMode = 1;
		}
		else if(sptflMode == "FASTER")
		{
			ptflMode = 2;
		}
		else if (sptflMode == "NEWMODE")
		{
			ptflMode = 3;
		}
		else if (sptflMode == "BASKET_ATM")
		{
			ptflMode = 4;
		}
		else if (sptflMode == "BASKET_EQUIV")
		{
			ptflMode = 5;
		}
		else if (sptflMode == "FRONTIER")
		{
			ptflMode = 6;
		}
		else
			throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"Portfolio Mode should be SUMMIT, BASKET_ATM, BASKET_EQUIV, NEWMODE, MANUAL or FASTER");
			
		portfolioMode[0] = ptflMode;
		
		//Calib Freq
		CCString ptflFreqXl = calibParamsXl[16];
		long lPtflFreq;
		if((lPtflFreq = ARM_ConvFrequency (ptflFreqXl, C_result)) == ARM_DEFAULT_ERR)
		{
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}
		int ptflFreq = (int)lPtflFreq;
		portfolioMode[1] = ptflFreq;

		//Calib Strike
		CCString ptflStrikeXl = calibParamsXl[17];
		string sptflStrikeXl = CCSTringToSTLString(ptflStrikeXl);
		string sptflStrike(ARM::stringGetUpper(sptflStrikeXl));
		int ptflStrike; 
		if(sptflStrike == "CALIB")
		{
			ptflStrike = 0;
		}
		else if (sptflStrike == "ATM")
		{
			ptflStrike = 1;
		}
		else if (sptflStrike == "MONEYNESS")
		{
			ptflStrike = 2;
		}
		else 
			throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"only ATM or CALIB mode are valid");
			
		portfolioMode[2] = ptflStrike;

		//Fix Boundary Flag
		CCString boundaryFlagXl = calibParamsXl[18];
		string sboundaryFlagXl = CCSTringToSTLString(boundaryFlagXl);
		string sboundaryFlag(ARM::stringGetUpper(sboundaryFlagXl));
		bool boundaryFlag; 
		if(sboundaryFlag == "Y")
		{
			boundaryFlag = true;
		}
		else if (sboundaryFlag == "N")
		{
			boundaryFlag = false;
		}
		else
			throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"only Y or N are valid ");

		//Approx Margin Flag
		CCString approxMarginFlagXl = calibParamsXl[19];
		string sapproxMarginFlagXl = CCSTringToSTLString(approxMarginFlagXl);
		string sapproxMarginFlag(ARM::stringGetUpper(sapproxMarginFlagXl));
		bool approxMarginFlag; 
		if(sapproxMarginFlag == "Y")
		{
			approxMarginFlag = true;
		}
		else if (sapproxMarginFlag == "N")
		{
			approxMarginFlag = false;
		}
		else
			throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"only Y or N are valid ");

		//Freeze Betas Flag
		CCString freezeBetasFlagXl = calibParamsXl[20];
		string sfreezeBetasFlagXl = CCSTringToSTLString(freezeBetasFlagXl);
		string sfreezeBetasFlag(ARM::stringGetUpper(sfreezeBetasFlagXl));
		bool freezeBetasFlag; 
		if(sfreezeBetasFlag == "Y")
		{
			freezeBetasFlag = true;
		}
		else if (sfreezeBetasFlag == "N")
		{
			freezeBetasFlag = false;
		}
		else
			throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"only Y or N are valid for FreezeBetasFlag ");
			
		//Generator Type
		CCString genTypeXl = calibParamsXl[21];
		string 	genType = CCSTringToSTLString(genTypeXl);
		genType = ARM::stringGetUpper(genType);

		//STM Calib Mode
		CCString mrsptflModeXl = calibParamsXl[22];
		string smrsptflModeXl = CCSTringToSTLString(mrsptflModeXl);
		string smrsptflMode(ARM::stringGetUpper(smrsptflModeXl));
		int mrsptflMode;
		if(smrsptflMode == "COLUMN")
		{
			mrsptflMode = 0;
		}
		else if (smrsptflMode == "CORREL")
		{
			mrsptflMode = 1;
		}
		else if (smrsptflMode == "SURDIAG")
		{
			mrsptflMode = 2;
		}
		else if (smrsptflMode == "ANTIDIAG")
		{
			mrsptflMode = 3;
		}
		else
			throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"STM Portfolio Mode should be COLUMN or ANTIDIAG");
		
		portfolioMode[3] = mrsptflMode;

		//Calculate Proba Flag
		CCString calculateProbaFlagXl = calibParamsXl[23];
		string scalculateProbaFlagXl = CCSTringToSTLString(calculateProbaFlagXl);
		string scalculateProbaFlag(ARM::stringGetUpper(scalculateProbaFlagXl));
		bool calculateProbaFlag; 
		if(scalculateProbaFlag == "Y")
		{
			calculateProbaFlag = true;
		}
		else if (scalculateProbaFlag == "N")
		{
			calculateProbaFlag = false;
		}
		else
			throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"only Y or N are valid for FreezeBetasFlag ");
		//End of calibParams

		// Lag for volatility calibration
		CCString volLagXl = calibParamsXl[24];
		string svolLag = CCSTringToSTLString(volLagXl );
		int volLag = atoi(svolLag .c_str());
		portfolioMode[4] = volLag;

		CCString curClass = LOCAL_GC_BERMUDASWAPTION_CLASS;

		long objId;

		CCString stringId = GetLastCurCellEnvValue ();
		
		if(!stringId)
			objId = -1;
		else
		{
			CCString prevClass = LocalGetStringObjectClass (stringId);
			
			objId = LocalGetNumObjectId (stringId);
				
			if(curClass != prevClass)
			{
				FreeCurCellContent ();
				objId = -1;
			}
		}

		retCode = ARMLOCAL_INITBERMUDASWAPTION (LocalGetNumObjectId(C_bsId),
												mktDataManagerId,
												mdmKeys,
												ctrlVariateVec,
												ctrlVariatePrices,
												modelParams,
												mrsCalibFlag,
												betaCalibFlag,
												numMethodType,
												amcIter,
												mcIter,
												maxBucketSize,
												genType,
												treeSteps,
												portfolioMode,
												boundaryFlag,
												approxMarginFlag,
												freezeBetasFlag,
												modelType,
												calculateProbaFlag,
												C_result,
												objId);

		if(retCode == ARM_OK)
		{
			objId = C_result.getLong ();

			LocalSetCurCellEnvValue (curClass, objId); 

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
	ARM_XL_TRY_BLOCK_END

	ARM_XL_CATCH_ARM_EXPT

	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_ARM_INITBERMUDASWAPTION" )

	return (LPXLOPER)&XL_result;
}



__declspec(dllexport) LPXLOPER WINAPI Local_PXL_ARM_INITBERMUDASWAPTION(LPXLOPER XL_bsId,
																		LPXLOPER XL_mktDataManager,
																		LPXLOPER XL_calibParams,
																		LPXLOPER XL_controlVariates)
{
	long retCode;
	ARM_result C_result;

	static XLOPER XL_result;

	ARM_XL_TRY_BLOCK_BEGIN
	{
		ARM_NOCALCIFWIZ();

		// C variable
		CCString C_bsId;

		// error
		static int error;
		static char* reason = "";

		//Bermuda Swaption Calculator to update
		XL_readStrCell(XL_bsId, C_bsId, " ARM_ERR: Bermuda Swaption id: object expected", C_result);
		
		//MKTDATAMANAGER: persistent object
		VECTOR<CCString> mktDataManager;
		XL_readStrVector(XL_mktDataManager,mktDataManager," ARM_ERR: Market datas: array of string expected",DOUBLE_TYPE,C_result);
		if(mktDataManager.size()<4)
		{
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}

		long mktDataManagerId = LocalGetNumObjectId (mktDataManager[0]);
		vector< string > mdmKeys(mktDataManager.size()-1);
		for(int i=1;i<mktDataManager.size();++i)
		mdmKeys[i-1]=CCSTringToSTLString(mktDataManager[i]);
		
		//CONTROLVARIATES
		VECTOR<CCString> controlVariates;
		VECTOR<CCString> controlVariatesDef(0);
		XL_readStrVectorWD(XL_controlVariates, controlVariates,controlVariatesDef, " ARM_ERR: control variates: array of string expected",DOUBLE_TYPE,C_result);
		int ctrlSize = controlVariates.size();
		vector<int>* ctrlVariateVec = new vector<int>(ctrlSize);
		vector<double>* ctrlVariatePrices = new vector<double>(ctrlSize);
		CCString cCurrCtrl;
		string sCurrCtrl;
		for (i=0; i<ctrlSize;i++)
		{
			cCurrCtrl = controlVariates[i];
			sCurrCtrl = CCSTringToSTLString(cCurrCtrl);
			(*ctrlVariateVec)[i] = atoi(sCurrCtrl.c_str());
			(*ctrlVariatePrices)[i] = 0;
		}

		//CALIB PARAMS
		//Model Parameters
		VECTOR<CCString> calibParamsXl;
		XL_readStrVector(XL_calibParams, calibParamsXl," ARM_ERR: Calib Parameters: array expected",DOUBLE_TYPE,C_result);
		if(calibParamsXl.size()<25)
		{
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}
		
		//Model Type
		CCString modelTypeXl = calibParamsXl[0];
		string smodelTypeXl = CCSTringToSTLString(modelTypeXl);
		string smodelType(ARM::stringGetUpper(smodelTypeXl));
		int modelType;
		if(smodelType == "SFRM1F")
		{
			modelType = 0;
		}
		else if (smodelType == "SFRM2F")
		{
			modelType = 1;
		}
		else if (smodelType == "QGM1F")
		{
			modelType = 2;
		}
		else if (smodelType == "HW1F")
		{
			modelType = 3;
		}
		else
			throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"only SFRM1F, SFRM2F, QGM1F or H&W are Valid ");


		vector<string> modelParams(7);

		CCString initVolXl = calibParamsXl[1];
		string sinitVol = CCSTringToSTLString(initVolXl);
		modelParams[0] = sinitVol;

		CCString initMrsXl = calibParamsXl[2];
		string sinitMrs = CCSTringToSTLString(initMrsXl);
		modelParams[1] = sinitMrs;

		CCString initBetaXl = calibParamsXl[3];
		string sinitBeta = CCSTringToSTLString(initBetaXl);
		modelParams[2] = sinitBeta;

		CCString initThetaXl = calibParamsXl[4];
		string sinitTheta = CCSTringToSTLString(initThetaXl);
		modelParams[3] = sinitTheta;

		CCString initSkewXl = calibParamsXl[5];
		string sinitSkew = CCSTringToSTLString(initSkewXl);
		modelParams[4] = sinitSkew;

		CCString initMrsLXl = calibParamsXl[6];
		string sinitMrsL = CCSTringToSTLString(initMrsLXl);
		modelParams[5] = sinitMrsL;

		CCString initMrsUXl = calibParamsXl[7];
		string sinitMrsU = CCSTringToSTLString(initMrsUXl);
		modelParams[6] = sinitMrsU;

		//Calib MRS/Beta ? (Y/N
		CCString calibMrsXl  = calibParamsXl[8];
		CCString calibBetaXl = calibParamsXl[9];
		string calibMrs  = CCSTringToSTLString(calibMrsXl);
		string calibBeta = CCSTringToSTLString(calibBetaXl);
		
		bool mrsCalibFlag = false;
		bool betaCalibFlag = false;
		if (calibMrs == "Y")
		{
			mrsCalibFlag = true;
		}
		if (calibBeta == "Y")
		{
			betaCalibFlag = true;
		}

		//Num Method Type
		CCString methodTypeXl = calibParamsXl[10];
		string smethodTypeXl = CCSTringToSTLString(methodTypeXl);
		string smethodType(ARM::stringGetUpper(smethodTypeXl));
		int numMethodType;
		if(smethodType == "TREE")
		{
			numMethodType = 1;
		}
		else if (smethodType == "MONTECARLO")
		{
			numMethodType = 0;
		}
		else if (smethodType == "PDE")
		{
			numMethodType = 2;
		}
		else
			throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"only TREE or MONTECARLO are valid ");

		//Num Amc Iter
		CCString amcIterXl = calibParamsXl[11];
		string samcIter = CCSTringToSTLString(amcIterXl);
		int amcIter = atoi(samcIter.c_str());
		
		//Num Mc Iter
		CCString mcIterXl = calibParamsXl[12];
		string smcIter = CCSTringToSTLString(mcIterXl);
		int mcIter = atoi(smcIter.c_str());

		//Num Max Bucket Size
		CCString maxBucketSizeXl = calibParamsXl[13];
		string smaxBucket = CCSTringToSTLString(maxBucketSizeXl );
		int maxBucketSize = atoi(smaxBucket.c_str());

		//Num Tree Steps
		CCString treeStepsXl = calibParamsXl[14];
		string streeSteps = CCSTringToSTLString(treeStepsXl);
		int treeSteps = atoi(streeSteps.c_str());

		//Portfolio Mode
		vector<int> portfolioMode(5);

		//Calib Mode
		CCString ptflModeXl = calibParamsXl[15];
		string sptflModeXl = CCSTringToSTLString(ptflModeXl);
		string sptflMode(ARM::stringGetUpper(sptflModeXl));
		int ptflMode;
		if(sptflMode == "SUMMIT")
		{
			ptflMode = 0;
		}
		else if (sptflMode == "MANUAL")
		{
			ptflMode = 1;
		}
		else if(sptflMode == "FASTER")
		{
			ptflMode = 2;
		}
		else if (sptflMode == "NEWMODE")
		{
			ptflMode = 3;
		}
		else if (sptflMode == "BASKET_ATM")
		{
			ptflMode = 4;
		}
		else if (sptflMode == "BASKET_EQUIV")
		{
			ptflMode = 5;
		}
		else if (sptflMode == "FRONTIER")
		{
			ptflMode = 6;
		}
		else
			throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"Portfolio Mode should be SUMMIT, BASKET_ATM, BASKET_EQUIV, NEWMODE, MANUAL or FASTER");
			
		portfolioMode[0] = ptflMode;
		
		//Calib Freq
		CCString ptflFreqXl = calibParamsXl[16];
		long lPtflFreq;
		if((lPtflFreq = ARM_ConvFrequency (ptflFreqXl, C_result)) == ARM_DEFAULT_ERR)
		{
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}
		int ptflFreq = (int)lPtflFreq;
		portfolioMode[1] = ptflFreq;

		//Calib Strike
		CCString ptflStrikeXl = calibParamsXl[17];
		string sptflStrikeXl = CCSTringToSTLString(ptflStrikeXl);
		string sptflStrike(ARM::stringGetUpper(sptflStrikeXl));
		int ptflStrike; 
		if(sptflStrike == "CALIB")
		{
			ptflStrike = 0;
		}
		else if (sptflStrike == "ATM")
		{
			ptflStrike = 1;
		}
		else if (sptflStrike == "MONEYNESS")
		{
			ptflStrike = 2;
		}
		else 
			throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"only ATM or CALIB mode are valid");
			
		portfolioMode[2] = ptflStrike;

		//Fix Boundary Flag
		CCString boundaryFlagXl = calibParamsXl[18];
		string sboundaryFlagXl = CCSTringToSTLString(boundaryFlagXl);
		string sboundaryFlag(ARM::stringGetUpper(sboundaryFlagXl));
		bool boundaryFlag; 
		if(sboundaryFlag == "Y")
		{
			boundaryFlag = true;
		}
		else if (sboundaryFlag == "N")
		{
			boundaryFlag = false;
		}
		else
			throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"only Y or N are valid ");

		//Approx Margin Flag
		CCString approxMarginFlagXl = calibParamsXl[19];
		string sapproxMarginFlagXl = CCSTringToSTLString(approxMarginFlagXl);
		string sapproxMarginFlag(ARM::stringGetUpper(sapproxMarginFlagXl));
		bool approxMarginFlag; 
		if(sapproxMarginFlag == "Y")
		{
			approxMarginFlag = true;
		}
		else if (sapproxMarginFlag == "N")
		{
			approxMarginFlag = false;
		}
		else
			throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"only Y or N are valid ");

		//Freeze Betas Flag
		CCString freezeBetasFlagXl = calibParamsXl[20];
		string sfreezeBetasFlagXl = CCSTringToSTLString(freezeBetasFlagXl);
		string sfreezeBetasFlag(ARM::stringGetUpper(sfreezeBetasFlagXl));
		bool freezeBetasFlag; 
		if(sfreezeBetasFlag == "Y")
		{
			freezeBetasFlag = true;
		}
		else if (sfreezeBetasFlag == "N")
		{
			freezeBetasFlag = false;
		}
		else
			throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"only Y or N are valid for FreezeBetasFlag ");
			
		//Generator Type
		CCString genTypeXl = calibParamsXl[21];
		string 	genType = CCSTringToSTLString(genTypeXl);
		genType = ARM::stringGetUpper(genType);

		//STM Calib Mode
		CCString mrsptflModeXl = calibParamsXl[22];
		string smrsptflModeXl = CCSTringToSTLString(mrsptflModeXl);
		string smrsptflMode(ARM::stringGetUpper(smrsptflModeXl));
		int mrsptflMode;
		if(smrsptflMode == "COLUMN")
		{
			mrsptflMode = 0;
		}
		else if (smrsptflMode == "CORREL")
		{
			mrsptflMode = 1;
		}
		else if (smrsptflMode == "SURDIAG")
		{
			mrsptflMode = 2;
		}
		else if (smrsptflMode == "ANTIDIAG")
		{
			mrsptflMode = 3;
		}
		else
			throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"STM Portfolio Mode should be COLUMN or ANTIDIAG");
		
		portfolioMode[3] = mrsptflMode;

		//Calculate Proba Flag
		CCString calculateProbaFlagXl = calibParamsXl[23];
		string scalculateProbaFlagXl = CCSTringToSTLString(calculateProbaFlagXl);
		string scalculateProbaFlag(ARM::stringGetUpper(scalculateProbaFlagXl));
		bool calculateProbaFlag; 
		if(scalculateProbaFlag == "Y")
		{
			calculateProbaFlag = true;
		}
		else if (scalculateProbaFlag == "N")
		{
			calculateProbaFlag = false;
		}
		else
			throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"only Y or N are valid for FreezeBetasFlag ");
		//End of calibParams

		// Lag for volatility calibration
		CCString volLagXl = calibParamsXl[24];
		string svolLag = CCSTringToSTLString(volLagXl );
		int volLag = atoi(svolLag .c_str());
		portfolioMode[4] = volLag;

		CCString curClass = LOCAL_GC_BERMUDASWAPTION_CLASS;

		long objId;

		CCString stringId;
		
		retCode = ARMLOCAL_INITBERMUDASWAPTION (LocalGetNumObjectId(C_bsId),
												mktDataManagerId,
												mdmKeys,
												ctrlVariateVec,
												ctrlVariatePrices,
												modelParams,
												mrsCalibFlag,
												betaCalibFlag,
												numMethodType,
												amcIter,
												mcIter,
												maxBucketSize,
												genType,
												treeSteps,
												portfolioMode,
												boundaryFlag,
												approxMarginFlag,
												freezeBetasFlag,
												modelType,
												calculateProbaFlag,
												C_result);

		
        if ( retCode == ARM_OK )
		{
			objId = C_result.getLong ();

			stringId = LocalMakeObjectId (objId, curClass);
		}

		if ( retCode == ARM_OK )
		{
			// FreeCurCellErr ();
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

	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_PXL_ARM_INITBERMUDASWAPTION" )

	return (LPXLOPER)&XL_result;
}



__declspec(dllexport) LPXLOPER WINAPI Local_ARM_INITSWAPTIONBERMUDA(LPXLOPER XL_bsId,
																	LPXLOPER XL_calibParams,
																	LPXLOPER XL_zcCpn,
																	LPXLOPER XL_swoptVc,
																	LPXLOPER XL_capVc,
																	LPXLOPER XL_capRo,
																	LPXLOPER XL_capNu,
																	LPXLOPER XL_capBeta,
																	LPXLOPER XL_swoptRo,
																	LPXLOPER XL_swoptNu,
																	LPXLOPER XL_swoptBeta,
																	LPXLOPER XL_normalModel,
																	LPXLOPER XL_controlVariates,
																	LPXLOPER XL_SigmaOrAlpha)
{
	long retCode;
	ARM_result C_result;

	static XLOPER XL_result;

	ARM_XL_TRY_BLOCK_BEGIN
	{
		ARM_NOCALCIFWIZ();
		
		// C variable
		CCString C_bsId;

		// error
		static int error;
		static char* reason = "";

		//Bermuda Swaption Calculator to update
		XL_readStrCell(XL_bsId, C_bsId, " ARM_ERR: Bermuda Swaption id: object expected", C_result);
		
		//CALIB PARAMS
		//Model Parameters
		VECTOR<CCString> calibParamsXl;
		XL_readStrVector(XL_calibParams, calibParamsXl," ARM_ERR: Calib Parameters: array expected",DOUBLE_TYPE,C_result);
		if(calibParamsXl.size()<25)
		{
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}
		
		//Model Type
		CCString modelTypeXl = calibParamsXl[0];
		string smodelTypeXl = CCSTringToSTLString(modelTypeXl);
		string smodelType(ARM::stringGetUpper(smodelTypeXl));
		int modelType;
		if(smodelType == "SFRM1F")
		{
			modelType = 0;
		}
		else if (smodelType == "SFRM2F")
		{
			modelType = 1;
		}
		else if (smodelType == "QGM1F")
		{
			modelType = 2;
		}
		else if (smodelType == "HW1F")
		{
			modelType = 3;
		}
		else
			throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"only SFRM1F, SFRM2F, QGM1F or H&W are Valid ");

		vector<string> modelParams(7);

		CCString initVolXl = calibParamsXl[1];
		string sinitVol = CCSTringToSTLString(initVolXl);
		double initVol = atof(sinitVol.c_str());
		modelParams[0] = sinitVol;

		CCString initMrsXl = calibParamsXl[2];
		string sinitMrs = CCSTringToSTLString(initMrsXl);
		modelParams[1] = sinitMrs;

		CCString initBetaXl = calibParamsXl[3];
		string sinitBeta = CCSTringToSTLString(initBetaXl);
		modelParams[2] = sinitBeta;

		CCString initThetaXl = calibParamsXl[4];
		string sinitTheta = CCSTringToSTLString(initThetaXl);
		modelParams[3] = sinitTheta;

		CCString initSkewXl = calibParamsXl[5];
		string sinitSkew = CCSTringToSTLString(initSkewXl);
		modelParams[4] = sinitSkew;

		CCString initMrsLXl = calibParamsXl[6];
		string sinitMrsL = CCSTringToSTLString(initMrsLXl);
		modelParams[5] = sinitMrsL;

		CCString initMrsUXl = calibParamsXl[7];
		string sinitMrsU = CCSTringToSTLString(initMrsUXl);
		modelParams[6] = sinitMrsU;

		//Calib MRS/Beta ? (Y/N
		CCString calibMrsXl  = calibParamsXl[8];
		CCString calibBetaXl = calibParamsXl[9];
		string calibMrs  = CCSTringToSTLString(calibMrsXl);
		string calibBeta = CCSTringToSTLString(calibBetaXl);
		
		bool mrsCalibFlag = false;
		bool betaCalibFlag = false;
		if (calibMrs == "Y")
		{
			mrsCalibFlag = true;
		}
		if (calibBeta == "Y")
		{
			betaCalibFlag = true;
		}

		//Num Method Type
		CCString methodTypeXl = calibParamsXl[10];
		string smethodTypeXl = CCSTringToSTLString(methodTypeXl);
		string smethodType(ARM::stringGetUpper(smethodTypeXl));
		int numMethodType;
		if(smethodType == "TREE")
		{
			numMethodType = 1;
		}
		else if (smethodType == "MONTECARLO")
		{
			numMethodType = 0;
		}
		else if (smethodType == "PDE")
		{
			numMethodType = 2;
		}
		else
			throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"only TREE or MONTECARLO are valid ");

		//Num Amc Iter
		CCString amcIterXl = calibParamsXl[11];
		string samcIter = CCSTringToSTLString(amcIterXl);
		int amcIter = atoi(samcIter.c_str());
		
		//Num Mc Iter
		CCString mcIterXl = calibParamsXl[12];
		string smcIter = CCSTringToSTLString(mcIterXl);
		int mcIter = atoi(smcIter.c_str());

		//Num Max Bucket Size
		CCString maxBucketSizeXl = calibParamsXl[13];
		string smaxBucket = CCSTringToSTLString(maxBucketSizeXl );
		int maxBucketSize = atoi(smaxBucket.c_str());

		//Num Tree Steps
		CCString treeStepsXl = calibParamsXl[14];
		string streeSteps = CCSTringToSTLString(treeStepsXl);
		int treeSteps = atoi(streeSteps.c_str());

		//Portfolio Mode
		vector<int> portfolioMode(5);

		//Calib Mode
		CCString ptflModeXl = calibParamsXl[15];
		string sptflModeXl = CCSTringToSTLString(ptflModeXl);
		string sptflMode(ARM::stringGetUpper(sptflModeXl));
		int ptflMode;
		if(sptflMode == "SUMMIT")
		{
			ptflMode = 0;
		}
		else if (sptflMode == "MANUAL")
		{
			ptflMode = 1;
		}
		else if(sptflMode == "FASTER")
		{
			ptflMode = 2;
		}
		else if (sptflMode == "NEWMODE")
		{
			ptflMode = 3;
		}
		else if (sptflMode == "BASKET_ATM")
		{
			ptflMode = 4;
		}
		else if (sptflMode == "BASKET_EQUIV")
		{
			ptflMode = 5;
		}
		else if (sptflMode == "FRONTIER")
		{
			ptflMode = 6;
		}
		else
			throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"Portfolio Mode should be SUMMIT, BASKET_ATM, BASKET_EQUIV, NEWMODE, MANUAL or FASTER");
			
		portfolioMode[0] = ptflMode;
		
		//Calib Freq
		CCString ptflFreqXl = calibParamsXl[16];
		long lPtflFreq;
		if((lPtflFreq = ARM_ConvFrequency (ptflFreqXl, C_result)) == ARM_DEFAULT_ERR)
		{
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}
		int ptflFreq = (int)lPtflFreq;
		portfolioMode[1] = ptflFreq;

		//Calib Strike
		CCString ptflStrikeXl = calibParamsXl[17];
		string sptflStrikeXl = CCSTringToSTLString(ptflStrikeXl);
		string sptflStrike(ARM::stringGetUpper(sptflStrikeXl));
		int ptflStrike; 
		if(sptflStrike == "CALIB")
		{
			ptflStrike = 0;
		}
		else if (sptflStrike == "ATM")
		{
			ptflStrike = 1;
		}
		else if (sptflStrike == "MONEYNESS")
		{
			ptflStrike = 2;
		}
		else 
			throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"only ATM or CALIB mode are valid");
			
		portfolioMode[2] = ptflStrike;

		//Fix Boundary Flag
		CCString boundaryFlagXl = calibParamsXl[18];
		string sboundaryFlagXl = CCSTringToSTLString(boundaryFlagXl);
		string sboundaryFlag(ARM::stringGetUpper(sboundaryFlagXl));
		bool boundaryFlag; 
		if(sboundaryFlag == "Y")
		{
			boundaryFlag = true;
		}
		else if (sboundaryFlag == "N")
		{
			boundaryFlag = false;
		}
		else
			throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"only Y or N are valid ");

		//Approx Margin Flag
		CCString approxMarginFlagXl = calibParamsXl[19];
		string sapproxMarginFlagXl = CCSTringToSTLString(approxMarginFlagXl);
		string sapproxMarginFlag(ARM::stringGetUpper(sapproxMarginFlagXl));
		bool approxMarginFlag; 
		if(sapproxMarginFlag == "Y")
		{
			approxMarginFlag = true;
		}
		else if (sapproxMarginFlag == "N")
		{
			approxMarginFlag = false;
		}
		else
			throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"only Y or N are valid ");

		//Freeze Betas Flag
		CCString freezeBetasFlagXl = calibParamsXl[20];
		string sfreezeBetasFlagXl = CCSTringToSTLString(freezeBetasFlagXl);
		string sfreezeBetasFlag(ARM::stringGetUpper(sfreezeBetasFlagXl));
		bool freezeBetasFlag; 
		if(sfreezeBetasFlag == "Y")
		{
			freezeBetasFlag = true;
		}
		else if (sfreezeBetasFlag == "N")
		{
			freezeBetasFlag = false;
		}
		else
			throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"only Y or N are valid for FreezeBetasFlag ");
			
		//Generator Type
		CCString genTypeXl = calibParamsXl[21];
		string 	genType = CCSTringToSTLString(genTypeXl);
		genType = ARM::stringGetUpper(genType);

		//STM Calib Mode
		CCString mrsptflModeXl = calibParamsXl[22];
		string smrsptflModeXl = CCSTringToSTLString(mrsptflModeXl);
		string smrsptflMode(ARM::stringGetUpper(smrsptflModeXl));
		int mrsptflMode;
		if(smrsptflMode == "COLUMN")
		{
			mrsptflMode = 0;
		}
		else if (smrsptflMode == "CORREL")
		{
			mrsptflMode = 1;
		}
		else if (smrsptflMode == "SURDIAG")
		{
			mrsptflMode = 2;
		}
		else if (smrsptflMode == "ANTIDIAG")
		{
			mrsptflMode = 3;
		}
		else
			throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"STM Portfolio Mode should be COLUMN or ANTIDIAG");
		
		portfolioMode[3] = mrsptflMode;

		//Calculate Proba Flag
		CCString calculateProbaFlagXl = calibParamsXl[23];
		string scalculateProbaFlagXl = CCSTringToSTLString(calculateProbaFlagXl);
		string scalculateProbaFlag(ARM::stringGetUpper(scalculateProbaFlagXl));
		bool calculateProbaFlag; 
		if(scalculateProbaFlag == "Y")
		{
			calculateProbaFlag = true;
		}
		else if (scalculateProbaFlag == "N")
		{
			calculateProbaFlag = false;
		}
		else
			throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"only Y or N are valid for FreezeBetasFlag ");
		//End of calibParams

		// Lag for volatility calibration
		CCString volLagXl = calibParamsXl[24];
		string svolLag = CCSTringToSTLString(volLagXl );
		int volLag = atoi(svolLag .c_str());
		portfolioMode[4] = volLag;

		//zcCpn
		CCString C_zcCpnId;
		long zcCpnId = ARM_NULL_OBJECT;
		XL_readStrCell( XL_zcCpn, C_zcCpnId, "ARM_ERR: zc cpn id: object expected", C_result);
		if ( string(C_zcCpnId) != "DEFAULT" )
			zcCpnId = LocalGetNumObjectId(C_zcCpnId);

		//swoptVc
		CCString C_swoptVcId;
		long swoptVcId = ARM_NULL_OBJECT;
		XL_readStrCell( XL_swoptVc, C_swoptVcId, "ARM_ERR: swopt vol: object expected", C_result);
		if ( string(C_swoptVcId) != "DEFAULT" )
			swoptVcId = LocalGetNumObjectId(C_swoptVcId);

		//capVc
		CCString C_capVcId;
		long capVcId = ARM_NULL_OBJECT;
		XL_readStrCell( XL_capVc, C_capVcId, "ARM_ERR: cap vol id: object expected", C_result);
		if ( string(C_capVcId) != "DEFAULT" )
			capVcId = LocalGetNumObjectId(C_capVcId);

		//capRo
		CCString C_capRoId;
		long capRoId = ARM_NULL_OBJECT;
		XL_readStrCellWD( XL_capRo,	C_capRoId, "DEFAULT"," ARM_ERR: capRo: object expected", C_result);	
		if ( string(C_capRoId) != "DEFAULT" )
			capRoId = LocalGetNumObjectId(C_capRoId);

		//capNu
		CCString C_capNuId;
		long capNuId = ARM_NULL_OBJECT;
		XL_readStrCellWD( XL_capNu,	C_capNuId, "DEFAULT"," ARM_ERR: capNu: object expected", C_result);	
		if ( string(C_capNuId) != "DEFAULT" )
			capNuId = LocalGetNumObjectId(C_capNuId);
		
		//capBeta
		CCString C_capBetaId;
		long capBetaId = ARM_NULL_OBJECT;
		XL_readStrCellWD( XL_capBeta,	C_capBetaId, "DEFAULT"," ARM_ERR: capBeta: object expected", C_result);	
		if ( string(C_capBetaId) != "DEFAULT" )
			capBetaId = LocalGetNumObjectId(C_capBetaId);
		
		//swoptRo
		CCString C_swoptRoId;
		long swoptRoId = ARM_NULL_OBJECT;
		XL_readStrCellWD( XL_swoptRo,	C_swoptRoId, "DEFAULT"," ARM_ERR: swopt ro: object expected", C_result);	
		if ( string(C_swoptRoId) != "DEFAULT" )
			swoptRoId = LocalGetNumObjectId(C_swoptRoId);
		
		//swoptNu	
		CCString C_swoptNuId;
		long swoptNuId = ARM_NULL_OBJECT;
		XL_readStrCellWD( XL_swoptNu,	C_swoptNuId, "DEFAULT"," ARM_ERR: swopt nu: object expected", C_result);	
		if ( string(C_swoptNuId) != "DEFAULT" )
			swoptNuId = LocalGetNumObjectId(C_swoptNuId);

		//swoptBeta
		CCString C_swoptBetaId;
		long swoptBetaId = ARM_NULL_OBJECT;
		XL_readStrCellWD( XL_swoptBeta,	C_swoptBetaId, "DEFAULT"," ARM_ERR: swopt beta: object expected", C_result);	
		if ( string(C_swoptBetaId) != "DEFAULT" )
			swoptBetaId = LocalGetNumObjectId(C_swoptBetaId);

		//sigmaOrAlpha
		CCString C_SigmaOrAlpha;
		XL_readStrCellWD(XL_SigmaOrAlpha, C_SigmaOrAlpha, "S", " ARM_ERR: SABR vol mode: S, SIGMA, A, ALPHA expected", C_result);
		C_SigmaOrAlpha.toUpper();
        
		long SigmaOrAlpha = 1; // sigma
		if (( C_SigmaOrAlpha == "A" ) || ( C_SigmaOrAlpha == "ALPHA" ))
			SigmaOrAlpha = 0;

		//normalModel
		CCString C_normalModelId;
		long normalModelId = ARM_NULL_OBJECT;
		XL_readStrCellWD( XL_normalModel, C_normalModelId, "DEFAULT"," ARM_ERR: normal model: object expected", C_result);	
		if ( string(C_normalModelId) != "DEFAULT" )
			normalModelId = LocalGetNumObjectId(C_normalModelId);

		//CONTROLVARIATES
		VECTOR<CCString> controlVariates;
		VECTOR<CCString> controlVariatesDef(0);
		XL_readStrVectorWD(XL_controlVariates, controlVariates,controlVariatesDef, " ARM_ERR: control variates: array of string expected",DOUBLE_TYPE,C_result);
		int ctrlSize						= controlVariates.size();
		vector<int>* ctrlVariateVec			= new vector<int>(ctrlSize);
		vector<double>* ctrlVariatePrices	= new vector<double>(ctrlSize);
		CCString cCurrCtrl;
		string sCurrCtrl;
		for (int i = 0 ; i < ctrlSize ; i++)
		{
			cCurrCtrl				= controlVariates[i];
			sCurrCtrl				= CCSTringToSTLString(cCurrCtrl);
			(*ctrlVariateVec)[i]	= atoi(sCurrCtrl.c_str());
			(*ctrlVariatePrices)[i] = 0;
		}

		CCString curClass = LOCAL_GC_BERMUDASWAPTION_CLASS;

		long objId;

		CCString stringId = GetLastCurCellEnvValue ();
		
		if(!stringId)
			objId = -1;
		else
		{
			CCString prevClass = LocalGetStringObjectClass (stringId);
			
			objId = LocalGetNumObjectId (stringId);
				
			if(curClass != prevClass)
			{
				FreeCurCellContent ();
				objId = -1;
			}
		}

		retCode = ARMLOCAL_INITSWAPTIONBERMUDA (//BermudaSwaption partial calculator
												LocalGetNumObjectId(C_bsId),
												//CalibParameters
												ctrlVariateVec,
												ctrlVariatePrices,
												modelParams,
												mrsCalibFlag,
												betaCalibFlag,
												numMethodType,
												amcIter,
												mcIter,
												maxBucketSize,
												genType,
												treeSteps,
												portfolioMode,
												boundaryFlag,
												approxMarginFlag,
												freezeBetasFlag,
												modelType,
												calculateProbaFlag,
												//MarketDatas
												zcCpnId,
												swoptVcId,
												capVcId,
												capRoId,
												capNuId,
												capBetaId,
												swoptRoId,
												swoptNuId,
												swoptBetaId,
												normalModelId,
												SigmaOrAlpha,
												C_result,
												objId);

		if(retCode == ARM_OK)
		{
			objId = C_result.getLong ();
			LocalSetCurCellEnvValue (curClass, objId); 
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
	ARM_XL_TRY_BLOCK_END

	ARM_XL_CATCH_ARM_EXPT

	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_ARM_INITSWAPTIONBERMUDA" )

	return (LPXLOPER)&XL_result;
}

__declspec(dllexport) LPXLOPER WINAPI Local_PXL_ARM_INITSWAPTIONBERMUDA(LPXLOPER XL_bsId,
																		LPXLOPER XL_calibParams,
																		LPXLOPER XL_zcCpn,
																		LPXLOPER XL_swoptVc,
																		LPXLOPER XL_capVc,
																		LPXLOPER XL_capRo,
																		LPXLOPER XL_capNu,
																		LPXLOPER XL_capBeta,
																		LPXLOPER XL_swoptRo,
																		LPXLOPER XL_swoptNu,
																		LPXLOPER XL_swoptBeta,
																		LPXLOPER XL_normalModel,
																		LPXLOPER XL_controlVariates,
																		LPXLOPER XL_SigmaOrAlpha)
{
	long retCode;
	ARM_result C_result;

	static XLOPER XL_result;

	ARM_XL_TRY_BLOCK_BEGIN
	{
		ARM_NOCALCIFWIZ();
		
		// C variable
		CCString C_bsId;

		// error
		static int error;
		static char* reason = "";

		//Bermuda Swaption Calculator to update
		XL_readStrCell(XL_bsId, C_bsId, " ARM_ERR: Bermuda Swaption id: object expected", C_result);
		
		//CALIB PARAMS
		//Model Parameters
		VECTOR<CCString> calibParamsXl;
		XL_readStrVector(XL_calibParams, calibParamsXl," ARM_ERR: Calib Parameters: array expected",DOUBLE_TYPE,C_result);
		if(calibParamsXl.size()<25)
		{
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}
		
		//Model Type
		CCString modelTypeXl = calibParamsXl[0];
		string smodelTypeXl = CCSTringToSTLString(modelTypeXl);
		string smodelType(ARM::stringGetUpper(smodelTypeXl));
		int modelType;
		if(smodelType == "SFRM1F")
		{
			modelType = 0;
		}
		else if (smodelType == "SFRM2F")
		{
			modelType = 1;
		}
		else if (smodelType == "QGM1F")
		{
			modelType = 2;
		}
		else if (smodelType == "HW1F")
		{
			modelType = 3;
		}
		else
			throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"only SFRM1F, SFRM2F, QGM1F or H&W are Valid ");

		vector<string> modelParams(7);

		CCString initVolXl = calibParamsXl[1];
		string sinitVol = CCSTringToSTLString(initVolXl);
		double initVol = atof(sinitVol.c_str());
		modelParams[0] = sinitVol;

		CCString initMrsXl = calibParamsXl[2];
		string sinitMrs = CCSTringToSTLString(initMrsXl);
		modelParams[1] = sinitMrs;

		CCString initBetaXl = calibParamsXl[3];
		string sinitBeta = CCSTringToSTLString(initBetaXl);
		modelParams[2] = sinitBeta;

		CCString initThetaXl = calibParamsXl[4];
		string sinitTheta = CCSTringToSTLString(initThetaXl);
		modelParams[3] = sinitTheta;

		CCString initSkewXl = calibParamsXl[5];
		string sinitSkew = CCSTringToSTLString(initSkewXl);
		modelParams[4] = sinitSkew;

		CCString initMrsLXl = calibParamsXl[6];
		string sinitMrsL = CCSTringToSTLString(initMrsLXl);
		modelParams[5] = sinitMrsL;

		CCString initMrsUXl = calibParamsXl[7];
		string sinitMrsU = CCSTringToSTLString(initMrsUXl);
		modelParams[6] = sinitMrsU;

		//Calib MRS/Beta ? (Y/N
		CCString calibMrsXl  = calibParamsXl[8];
		CCString calibBetaXl = calibParamsXl[9];
		string calibMrs  = CCSTringToSTLString(calibMrsXl);
		string calibBeta = CCSTringToSTLString(calibBetaXl);
		
		bool mrsCalibFlag = false;
		bool betaCalibFlag = false;
		if (calibMrs == "Y")
		{
			mrsCalibFlag = true;
		}
		if (calibBeta == "Y")
		{
			betaCalibFlag = true;
		}

		//Num Method Type
		CCString methodTypeXl = calibParamsXl[10];
		string smethodTypeXl = CCSTringToSTLString(methodTypeXl);
		string smethodType(ARM::stringGetUpper(smethodTypeXl));
		int numMethodType;
		if(smethodType == "TREE")
		{
			numMethodType = 1;
		}
		else if (smethodType == "MONTECARLO")
		{
			numMethodType = 0;
		}
		else if (smethodType == "PDE")
		{
			numMethodType = 2;
		}
		else
			throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"only TREE or MONTECARLO are valid ");

		//Num Amc Iter
		CCString amcIterXl = calibParamsXl[11];
		string samcIter = CCSTringToSTLString(amcIterXl);
		int amcIter = atoi(samcIter.c_str());
		
		//Num Mc Iter
		CCString mcIterXl = calibParamsXl[12];
		string smcIter = CCSTringToSTLString(mcIterXl);
		int mcIter = atoi(smcIter.c_str());

		//Num Max Bucket Size
		CCString maxBucketSizeXl = calibParamsXl[13];
		string smaxBucket = CCSTringToSTLString(maxBucketSizeXl );
		int maxBucketSize = atoi(smaxBucket.c_str());

		//Num Tree Steps
		CCString treeStepsXl = calibParamsXl[14];
		string streeSteps = CCSTringToSTLString(treeStepsXl);
		int treeSteps = atoi(streeSteps.c_str());

		//Portfolio Mode
		vector<int> portfolioMode(5);

		//Calib Mode
		CCString ptflModeXl = calibParamsXl[15];
		string sptflModeXl = CCSTringToSTLString(ptflModeXl);
		string sptflMode(ARM::stringGetUpper(sptflModeXl));
		int ptflMode;
		if(sptflMode == "SUMMIT")
		{
			ptflMode = 0;
		}
		else if (sptflMode == "MANUAL")
		{
			ptflMode = 1;
		}
		else if(sptflMode == "FASTER")
		{
			ptflMode = 2;
		}
		else if (sptflMode == "NEWMODE")
		{
			ptflMode = 3;
		}
		else if (sptflMode == "BASKET_ATM")
		{
			ptflMode = 4;
		}
		else if (sptflMode == "BASKET_EQUIV")
		{
			ptflMode = 5;
		}
		else if (sptflMode == "FRONTIER")
		{
			ptflMode = 6;
		}
		else
			throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"Portfolio Mode should be SUMMIT, BASKET_ATM, BASKET_EQUIV, NEWMODE, MANUAL or FASTER");
			
		portfolioMode[0] = ptflMode;
		
		//Calib Freq
		CCString ptflFreqXl = calibParamsXl[16];
		long lPtflFreq;
		if((lPtflFreq = ARM_ConvFrequency (ptflFreqXl, C_result)) == ARM_DEFAULT_ERR)
		{
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}
		int ptflFreq = (int)lPtflFreq;
		portfolioMode[1] = ptflFreq;

		//Calib Strike
		CCString ptflStrikeXl = calibParamsXl[17];
		string sptflStrikeXl = CCSTringToSTLString(ptflStrikeXl);
		string sptflStrike(ARM::stringGetUpper(sptflStrikeXl));
		int ptflStrike; 
		if(sptflStrike == "CALIB")
		{
			ptflStrike = 0;
		}
		else if (sptflStrike == "ATM")
		{
			ptflStrike = 1;
		}
		else if (sptflStrike == "MONEYNESS")
		{
			ptflStrike = 2;
		}
		else 
			throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"only ATM or CALIB mode are valid");
			
		portfolioMode[2] = ptflStrike;

		//Fix Boundary Flag
		CCString boundaryFlagXl = calibParamsXl[18];
		string sboundaryFlagXl = CCSTringToSTLString(boundaryFlagXl);
		string sboundaryFlag(ARM::stringGetUpper(sboundaryFlagXl));
		bool boundaryFlag; 
		if(sboundaryFlag == "Y")
		{
			boundaryFlag = true;
		}
		else if (sboundaryFlag == "N")
		{
			boundaryFlag = false;
		}
		else
			throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"only Y or N are valid ");

		//Approx Margin Flag
		CCString approxMarginFlagXl = calibParamsXl[19];
		string sapproxMarginFlagXl = CCSTringToSTLString(approxMarginFlagXl);
		string sapproxMarginFlag(ARM::stringGetUpper(sapproxMarginFlagXl));
		bool approxMarginFlag; 
		if(sapproxMarginFlag == "Y")
		{
			approxMarginFlag = true;
		}
		else if (sapproxMarginFlag == "N")
		{
			approxMarginFlag = false;
		}
		else
			throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"only Y or N are valid ");

		//Freeze Betas Flag
		CCString freezeBetasFlagXl = calibParamsXl[20];
		string sfreezeBetasFlagXl = CCSTringToSTLString(freezeBetasFlagXl);
		string sfreezeBetasFlag(ARM::stringGetUpper(sfreezeBetasFlagXl));
		bool freezeBetasFlag; 
		if(sfreezeBetasFlag == "Y")
		{
			freezeBetasFlag = true;
		}
		else if (sfreezeBetasFlag == "N")
		{
			freezeBetasFlag = false;
		}
		else
			throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"only Y or N are valid for FreezeBetasFlag ");
			
		//Generator Type
		CCString genTypeXl = calibParamsXl[21];
		string 	genType = CCSTringToSTLString(genTypeXl);
		genType = ARM::stringGetUpper(genType);

		//STM Calib Mode
		CCString mrsptflModeXl = calibParamsXl[22];
		string smrsptflModeXl = CCSTringToSTLString(mrsptflModeXl);
		string smrsptflMode(ARM::stringGetUpper(smrsptflModeXl));
		int mrsptflMode;
		if(smrsptflMode == "COLUMN")
		{
			mrsptflMode = 0;
		}
		else if (smrsptflMode == "CORREL")
		{
			mrsptflMode = 1;
		}
		else if (smrsptflMode == "SURDIAG")
		{
			mrsptflMode = 2;
		}
		else if (smrsptflMode == "ANTIDIAG")
		{
			mrsptflMode = 3;
		}
		else
			throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"STM Portfolio Mode should be COLUMN or ANTIDIAG");
		
		portfolioMode[3] = mrsptflMode;

		//Calculate Proba Flag
		CCString calculateProbaFlagXl = calibParamsXl[23];
		string scalculateProbaFlagXl = CCSTringToSTLString(calculateProbaFlagXl);
		string scalculateProbaFlag(ARM::stringGetUpper(scalculateProbaFlagXl));
		bool calculateProbaFlag; 
		if(scalculateProbaFlag == "Y")
		{
			calculateProbaFlag = true;
		}
		else if (scalculateProbaFlag == "N")
		{
			calculateProbaFlag = false;
		}
		else
			throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"only Y or N are valid for FreezeBetasFlag ");
		//End of calibParams

		// Lag for volatility calibration
		CCString volLagXl = calibParamsXl[24];
		string svolLag = CCSTringToSTLString(volLagXl );
		int volLag = atoi(svolLag .c_str());
		portfolioMode[4] = volLag;

		//zcCpn
		CCString C_zcCpnId;
		long zcCpnId = ARM_NULL_OBJECT;
		XL_readStrCell( XL_zcCpn, C_zcCpnId, "ARM_ERR: zc cpn id: object expected", C_result);
		if ( string(C_zcCpnId) != "DEFAULT" )
			zcCpnId = LocalGetNumObjectId(C_zcCpnId);

		//swoptVc
		CCString C_swoptVcId;
		long swoptVcId = ARM_NULL_OBJECT;
		XL_readStrCell( XL_swoptVc, C_swoptVcId, "ARM_ERR: swopt vol: object expected", C_result);
		if ( string(C_swoptVcId) != "DEFAULT" )
			swoptVcId = LocalGetNumObjectId(C_swoptVcId);

		//capVc
		CCString C_capVcId;
		long capVcId = ARM_NULL_OBJECT;
		XL_readStrCell( XL_capVc, C_capVcId, "ARM_ERR: cap vol id: object expected", C_result);
		if ( string(C_capVcId) != "DEFAULT" )
			capVcId = LocalGetNumObjectId(C_capVcId);

		//capRo
		CCString C_capRoId;
		long capRoId = ARM_NULL_OBJECT;
		XL_readStrCellWD( XL_capRo,	C_capRoId, "DEFAULT"," ARM_ERR: capRo: object expected", C_result);	
		if ( string(C_capRoId) != "DEFAULT" )
			capRoId = LocalGetNumObjectId(C_capRoId);

		//capNu
		CCString C_capNuId;
		long capNuId = ARM_NULL_OBJECT;
		XL_readStrCellWD( XL_capNu,	C_capNuId, "DEFAULT"," ARM_ERR: capNu: object expected", C_result);	
		if ( string(C_capNuId) != "DEFAULT" )
			capNuId = LocalGetNumObjectId(C_capNuId);
		
		//capBeta
		CCString C_capBetaId;
		long capBetaId = ARM_NULL_OBJECT;
		XL_readStrCellWD( XL_capBeta,	C_capBetaId, "DEFAULT"," ARM_ERR: capBeta: object expected", C_result);	
		if ( string(C_capBetaId) != "DEFAULT" )
			capBetaId = LocalGetNumObjectId(C_capBetaId);
		
		//swoptRo
		CCString C_swoptRoId;
		long swoptRoId = ARM_NULL_OBJECT;
		XL_readStrCellWD( XL_swoptRo,	C_swoptRoId, "DEFAULT"," ARM_ERR: swopt ro: object expected", C_result);	
		if ( string(C_swoptRoId) != "DEFAULT" )
			swoptRoId = LocalGetNumObjectId(C_swoptRoId);
		
		//swoptNu	
		CCString C_swoptNuId;
		long swoptNuId = ARM_NULL_OBJECT;
		XL_readStrCellWD( XL_swoptNu,	C_swoptNuId, "DEFAULT"," ARM_ERR: swopt nu: object expected", C_result);	
		if ( string(C_swoptNuId) != "DEFAULT" )
			swoptNuId = LocalGetNumObjectId(C_swoptNuId);

		//swoptBeta
		CCString C_swoptBetaId;
		long swoptBetaId = ARM_NULL_OBJECT;
		XL_readStrCellWD( XL_swoptBeta,	C_swoptBetaId, "DEFAULT"," ARM_ERR: swopt beta: object expected", C_result);	
		if ( string(C_swoptBetaId) != "DEFAULT" )
			swoptBetaId = LocalGetNumObjectId(C_swoptBetaId);

		//sigma or alpha
		CCString C_SigmaOrAlpha;
		XL_readStrCellWD(XL_SigmaOrAlpha, C_SigmaOrAlpha, "S", " ARM_ERR: SABR vol mode: S, SIGMA, A, ALPHA expected", C_result);
		C_SigmaOrAlpha.toUpper();
        
		long SigmaOrAlpha = 1; // sigma
		if (( C_SigmaOrAlpha == "A" ) || ( C_SigmaOrAlpha == "ALPHA" ))
			SigmaOrAlpha = 0;

		//normalModel
		CCString C_normalModelId;
		long normalModelId = ARM_NULL_OBJECT;
		XL_readStrCellWD( XL_normalModel, C_normalModelId, "DEFAULT"," ARM_ERR: normal model: object expected", C_result);	
		if ( string(C_normalModelId) != "DEFAULT" )
			normalModelId = LocalGetNumObjectId(C_normalModelId);

		//CONTROLVARIATES
		VECTOR<CCString> controlVariates;
		VECTOR<CCString> controlVariatesDef(0);
		XL_readStrVectorWD(XL_controlVariates, controlVariates,controlVariatesDef, " ARM_ERR: control variates: array of string expected",DOUBLE_TYPE,C_result);
		int ctrlSize						= controlVariates.size();
		vector<int>* ctrlVariateVec			= new vector<int>(ctrlSize);
		vector<double>* ctrlVariatePrices	= new vector<double>(ctrlSize);
		CCString cCurrCtrl;
		string sCurrCtrl;
		for (int i = 0 ; i < ctrlSize ; i++)
		{
			cCurrCtrl				= controlVariates[i];
			sCurrCtrl				= CCSTringToSTLString(cCurrCtrl);
			(*ctrlVariateVec)[i]	= atoi(sCurrCtrl.c_str());
			(*ctrlVariatePrices)[i] = 0;
		}

		CCString curClass = LOCAL_GC_BERMUDASWAPTION_CLASS;

		long objId;
		CCString stringId;

		retCode = ARMLOCAL_INITSWAPTIONBERMUDA (//BermudaSwaption partial calculator
												LocalGetNumObjectId(C_bsId),
												//CalibParameters
												ctrlVariateVec,
												ctrlVariatePrices,
												modelParams,
												mrsCalibFlag,
												betaCalibFlag,
												numMethodType,
												amcIter,
												mcIter,
												maxBucketSize,
												genType,
												treeSteps,
												portfolioMode,
												boundaryFlag,
												approxMarginFlag,
												freezeBetasFlag,
												modelType,
												calculateProbaFlag,
												//MarketDatas
												zcCpnId,
												swoptVcId,
												capVcId,
												capRoId,
												capNuId,
												capBetaId,
												swoptRoId,
												swoptNuId,
												swoptBetaId,
												normalModelId,
												SigmaOrAlpha,
												C_result);

		if(retCode == ARM_OK)
		{
			objId = C_result.getLong ();
			stringId = LocalMakeObjectId (objId, curClass);
		}

		if(retCode == ARM_OK)
		{
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

	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_PXL_ARM_INITSWAPTIONBERMUDA" )

	return (LPXLOPER)&XL_result;
}



__declspec(dllexport) LPXLOPER WINAPI Local_ARM_INITCAPTION(LPXLOPER XL_captionId,
															LPXLOPER XL_mktDataManager,
															LPXLOPER XL_calibParams,
															LPXLOPER XL_productsToPrice)
{
	long retCode;
	ARM_result C_result;

	static XLOPER XL_result;

	ARM_XL_TRY_BLOCK_BEGIN
	{
		ARM_NOCALCIFWIZ();

		// C variable
		CCString C_captionId;

		// error
		static int error;
		static char* reason = "";

		// Caption Calculator to initialise
		XL_readStrCell(XL_captionId, C_captionId, " ARM_ERR: Caption id: object expected", C_result);
		
		// MKTDATAMANAGER: persistent object
		VECTOR<CCString> mktDataManager;
		XL_readStrVector(XL_mktDataManager,mktDataManager," ARM_ERR: Market datas: array of string expected",DOUBLE_TYPE,C_result);
		if(mktDataManager.size()<10)
		{
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}

		long mktDataManagerId = LocalGetNumObjectId (mktDataManager[0]);

		vector<string> mdmKeys(mktDataManager.size()-1);
		for (int i = 1; i < mktDataManager.size(); ++i)
			mdmKeys[i-1] = CCSTringToSTLString(mktDataManager[i]);
		
		// CALIB PARAMS
		VECTOR<CCString> C_calibParams;
		XL_readStrVector(XL_calibParams, C_calibParams," ARM_ERR: Calib Parameters: array expected",DOUBLE_TYPE,C_result);
		if(C_calibParams.size()<4)
		{
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}

		int factorNb = atoi(C_calibParams[0].c_str());
		int SFRMVolType = K_DIAG;
		if (C_calibParams[1] == "ROW")
			SFRMVolType = K_ROW;
		string SwoptCalibMode = CCSTringToSTLString(C_calibParams[2]);
		string BetaCalibMode = CCSTringToSTLString(C_calibParams[3]);

		// Products to price
		VECTOR<CCString> C_productsToPrice;
		VECTOR<CCString> C_productsToPriceDef(0);
		XL_readStrVectorWD(XL_productsToPrice, C_productsToPrice, C_productsToPriceDef, " ARM_ERR: control variates: array of string expected",DOUBLE_TYPE,C_result);
		int size = C_productsToPrice.size();
		std::deque<bool> productsToPrice(size);
		for (i = 0; i < size; i++)
		{
			productsToPrice[i] = (C_productsToPrice[i] == "Y" || C_productsToPrice[i] == "YES");
		}

		/// WARNING !!
		/// we have extracted an EXOTIC from Summit
		/// but we will return only the CAPTION leg
		/// so, instead of a Portfolio object, we will get a CAPTION CALCULATOR object
		CCString curClass = LOCAL_GC_CAPTION_CLASS;

		long objId;

		CCString stringId = GetLastCurCellEnvValue ();
		
		if(!stringId)
			objId = -1;
		else
		{
			CCString prevClass = LocalGetStringObjectClass (stringId);
			
			objId = LocalGetNumObjectId (stringId);
				
			if(curClass != prevClass)
			{
				FreeCurCellContent ();
				objId = -1;
			}
		}

		retCode = ARMLOCAL_INITCAPTION (LocalGetNumObjectId(C_captionId),
										mktDataManagerId,
										mdmKeys,
										factorNb,
										SFRMVolType,
										SwoptCalibMode,
										BetaCalibMode,
										productsToPrice,
										C_result,
										objId);

		if(retCode == ARM_OK)
		{
			objId = C_result.getLong ();

			LocalSetCurCellEnvValue (curClass, objId); 

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
	ARM_XL_TRY_BLOCK_END

	ARM_XL_CATCH_ARM_EXPT

	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_ARM_INITCAPTION" )

	return (LPXLOPER)&XL_result;
}



__declspec(dllexport) LPXLOPER WINAPI Local_PXL_ARM_INITCAPTION(LPXLOPER XL_captionId,
																LPXLOPER XL_mktDataManager,
																LPXLOPER XL_calibParams,
																LPXLOPER XL_productsToPrice)
{
	long retCode;
	ARM_result C_result;

	static XLOPER XL_result;

	ARM_XL_TRY_BLOCK_BEGIN
	{
		ARM_NOCALCIFWIZ();

		// C variable
		CCString C_captionId;

		// error
		static int error;
		static char* reason = "";

		// Caption Calculator to initialise
		XL_readStrCell(XL_captionId, C_captionId, " ARM_ERR: Caption id: object expected", C_result);
		
		// MKTDATAMANAGER: persistent object
		VECTOR<CCString> mktDataManager;
		XL_readStrVector(XL_mktDataManager,mktDataManager," ARM_ERR: Market datas: array of string expected",DOUBLE_TYPE,C_result);
		if(mktDataManager.size()<10)
		{
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}

		long mktDataManagerId = LocalGetNumObjectId (mktDataManager[0]);

		vector<string> mdmKeys(mktDataManager.size()-1);
		for (int i = 1; i < mktDataManager.size(); ++i)
			mdmKeys[i-1] = CCSTringToSTLString(mktDataManager[i]);
		
		// CALIB PARAMS
		VECTOR<CCString> C_calibParams;
		XL_readStrVector(XL_calibParams, C_calibParams," ARM_ERR: Calib Parameters: array expected",DOUBLE_TYPE,C_result);
		if(C_calibParams.size()<4)
		{
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}

		int factorNb = atoi(C_calibParams[0].c_str());
		int SFRMVolType = K_DIAG;
		if (C_calibParams[1] == "ROW")
			SFRMVolType = K_ROW;
		string SwoptCalibMode = CCSTringToSTLString(C_calibParams[2]);
		string BetaCalibMode = CCSTringToSTLString(C_calibParams[3]);

		// Products to price
		VECTOR<CCString> C_productsToPrice;
		VECTOR<CCString> C_productsToPriceDef(0);
		XL_readStrVectorWD(XL_productsToPrice, C_productsToPrice, C_productsToPriceDef, " ARM_ERR: control variates: array of string expected",DOUBLE_TYPE,C_result);
		int size = C_productsToPrice.size();
		std::deque<bool> productsToPrice(size);
		for (i = 0; i < size; i++)
		{
			productsToPrice[i] = (C_productsToPrice[i] == "Y" || C_productsToPrice[i] == "YES");
		}

		/// WARNING !!
		/// we have extracted an EXOTIC from Summit
		/// but we will return only the CAPTION leg
		/// so, instead of a PORTFOLIO object, we will get a CAPTION CALCULATOR object
		CCString curClass = LOCAL_GC_CAPTION_CLASS;

		long objId;

		CCString stringId;

		retCode = ARMLOCAL_INITCAPTION (LocalGetNumObjectId(C_captionId),
										mktDataManagerId,
										mdmKeys,
										factorNb,
										SFRMVolType,
										SwoptCalibMode,
										BetaCalibMode,
										productsToPrice,
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
	}
	ARM_XL_TRY_BLOCK_END

	ARM_XL_CATCH_ARM_EXPT

	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_PXL_ARM_INITCAPTION" )

	return (LPXLOPER)&XL_result;
}



__declspec(dllexport) LPXLOPER WINAPI Local_ARM_INITCSO_withMDM(LPXLOPER XL_csoId,
																LPXLOPER XL_mktDataManager,
																LPXLOPER XL_modelParams,
																LPXLOPER XL_calibParams,
																LPXLOPER XL_productsToPrice)
{
	long retCode;
	ARM_result C_result;

	static XLOPER XL_result;

	ARM_XL_TRY_BLOCK_BEGIN
	{
		ARM_NOCALCIFWIZ();

		// C variable
		CCString C_csoId;

		// error
		static int error;
		static char* reason = "";

		// CSO Calculator to initialise
		XL_readStrCell(XL_csoId, C_csoId, " ARM_ERR: cso id: object expected", C_result);
		
		// MKTDATAMANAGER: persistent object
		VECTOR<CCString> mktDataManager;
		XL_readStrVector(XL_mktDataManager,mktDataManager," ARM_ERR: Market data: array of string expected",DOUBLE_TYPE,C_result);
        if(mktDataManager.size()<1)
		{
			throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"Mkt Data size should be greater than 1");
		}
		long mktDataManagerId = LocalGetNumObjectId (mktDataManager[0]);

		vector<string> mdmKeys(mktDataManager.size()-1);
		for (int i = 1; i < mktDataManager.size(); ++i)
			mdmKeys[i-1] = CCSTringToSTLString(mktDataManager[i]);
		
		// MODEL PARAMS
		VECTOR<double> C_modelParams;
		VECTOR<double> C_modelParams_def(1, 0.0);
		XL_readNumVectorWD(XL_modelParams, C_modelParams,C_modelParams_def," ARM_ERR: Model Parameters: array expected",C_result);

		// CALIB PARAMS
		VECTOR<CCString> C_calibParams;
		VECTOR<CCString> C_calibParams_def(3);
		C_calibParams_def[0]= CCString("DIAG"); // diagonal calibration
		C_calibParams_def[1]= CCString("ATM");  // atm swaptions
		C_calibParams_def[2]= CCString("HWM1F");  // hw1f model
		XL_readStrVectorWD(XL_calibParams, C_calibParams,C_calibParams_def," ARM_ERR: Calib Parameters: array expected",DOUBLE_TYPE,C_result);
		XL_readStrVector(XL_calibParams, C_calibParams," ARM_ERR: Calib Parameters: array expected",DOUBLE_TYPE,C_result);
		if(C_calibParams.size()>8)
		{
			throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"Calib Parameters: calibType,calibStrikeType,modelType,pricingMethod ");
		}
		vector<string> calibParams(C_calibParams.size());
		for (i = 0; i < C_calibParams.size(); i++)
			calibParams[i] = CCSTringToSTLString(C_calibParams[i]);

		// Products to price
		VECTOR<CCString> C_productsToPrice;
		VECTOR<CCString> C_productsToPriceDef(5, "N");
		XL_readStrVectorWD(XL_productsToPrice, C_productsToPrice, C_productsToPriceDef, " ARM_ERR: control variates: array of string expected",DOUBLE_TYPE,C_result);
		int size = C_productsToPrice.size();
		std::deque<bool> productsToPrice(6,false);
		for (i = 0; i < size; i++)
		{
			productsToPrice[i] = (C_productsToPrice[i] == "Y" || C_productsToPrice[i] == "YES");
		}

		CCString curClass = LOCAL_GC_CALLABLE_SO_CLASS;

		long objId;

		CCString stringId = GetLastCurCellEnvValue ();
		
		if (!stringId)
		   objId = -1;
		else
		{
			CCString prevClass = LocalGetStringObjectClass (stringId);
			
			objId = LocalGetNumObjectId (stringId);
				
			if ( curClass != prevClass )
			{
			   FreeCurCellContent();
				
               objId = -1;
			}
		}

		retCode = ARMLOCAL_INITCSO (LocalGetNumObjectId(C_csoId),
									mktDataManagerId,
									// mdmKeys, : NOW irrelevant
									C_modelParams,
									calibParams,
									productsToPrice,
									C_result,
									objId);

		if ( retCode == ARM_OK )
		{
			objId = C_result.getLong ();

			LocalSetCurCellEnvValue (curClass, objId); 

			stringId = LocalMakeObjectId (objId, curClass);
		}

		if ( retCode == ARM_OK )
		{
			FreeCurCellErr ();
		
            XL_result.xltype  = xltypeStr;
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

	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_ARM_INITCSO_withMDM" )

	return (LPXLOPER)&XL_result;
}



__declspec(dllexport) LPXLOPER WINAPI Local_PXL_ARM_INITCSO_withMDM(LPXLOPER XL_csoId,
																	LPXLOPER XL_mktDataManager,
																	LPXLOPER XL_modelParams,
																	LPXLOPER XL_calibParams,
																	LPXLOPER XL_productsToPrice)
{
	long retCode;
	ARM_result C_result;

	static XLOPER XL_result;

	ARM_XL_TRY_BLOCK_BEGIN
	{
		ARM_NOCALCIFWIZ();

		// C variable
		CCString C_csoId;

		// error
		static int error;
		static char* reason = "";

		// CSO Calculator to initialise
		XL_readStrCell(XL_csoId, C_csoId, " ARM_ERR: cso id: object expected", C_result);
		
		// MKTDATAMANAGER: persistent object
		VECTOR<CCString> mktDataManager;
		XL_readStrVector(XL_mktDataManager,mktDataManager," ARM_ERR: Market data: array of string expected",DOUBLE_TYPE,C_result);
		if(mktDataManager.size()!=9)
		{
			throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"MDM + 8 model parameters expected");
		}
		long mktDataManagerId = LocalGetNumObjectId (mktDataManager[0]);

		vector<string> mdmKeys(mktDataManager.size()-1);
		for (int i = 1; i < mktDataManager.size(); ++i)
			mdmKeys[i-1] = CCSTringToSTLString(mktDataManager[i]);
		
		// MODEL PARAMS
		VECTOR<double> C_modelParams;
		VECTOR<double> C_modelParams_def(1, 0.0);
		XL_readNumVectorWD(XL_modelParams, C_modelParams,C_modelParams_def," ARM_ERR: Model Parameters: array expected",C_result);

		// CALIB PARAMS
		VECTOR<CCString> C_calibParams;
		VECTOR<CCString> C_calibParams_def(3);
		C_calibParams_def[0]= CCString("DIAG"); // diagonal calibration
		C_calibParams_def[1]= CCString("ATM");  // atm swaptions
		C_calibParams_def[2]= CCString("HWM1F");  // hw1f model
		XL_readStrVectorWD(XL_calibParams, C_calibParams,C_calibParams_def," ARM_ERR: Calib Parameters: array expected",DOUBLE_TYPE,C_result);
		XL_readStrVector(XL_calibParams, C_calibParams," ARM_ERR: Calib Parameters: array expected",DOUBLE_TYPE,C_result);
		if(C_calibParams.size()>8)
		{
			throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"Calib Parameters: calibType,calibStrikeType,modelType,pricingMethod ");
		}
		vector<string> calibParams(C_calibParams.size());
		for (i = 0; i < C_calibParams.size(); i++)
			calibParams[i] = CCSTringToSTLString(C_calibParams[i]);

		// Products to price
		VECTOR<CCString> C_productsToPrice;
		VECTOR<CCString> C_productsToPriceDef(5, "N");
		XL_readStrVectorWD(XL_productsToPrice, C_productsToPrice, C_productsToPriceDef, " ARM_ERR: control variates: array of string expected",DOUBLE_TYPE,C_result);
		int size = C_productsToPrice.size();
		std::deque<bool> productsToPrice(6,false);
		for (i = 0; i < size; i++)
		{
			productsToPrice[i] = (C_productsToPrice[i] == "Y" || C_productsToPrice[i] == "YES");
		}

		CCString curClass = LOCAL_GC_CALLABLE_SO_CLASS;

		long objId;

		CCString stringId;

		retCode = ARMLOCAL_INITCSO (LocalGetNumObjectId(C_csoId),
									mktDataManagerId,
									// mdmKeys, : NOW irrelevant
									C_modelParams,
									calibParams,
									productsToPrice,
									C_result);

		if ( retCode == ARM_OK )
		{
			objId = C_result.getLong ();

			stringId = LocalMakeObjectId (objId, curClass);
		}

		if ( retCode == ARM_OK )
		{
			FreeCurCellErr ();
		
            XL_result.xltype  = xltypeStr;
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

	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_PXL_ARM_INITCSO_withMDM" )

	return (LPXLOPER)&XL_result;
}



__declspec(dllexport) LPXLOPER WINAPI Local_ARM_INITCSO(LPXLOPER XL_csoId,
														LPXLOPER XL_zcId,
														LPXLOPER XL_capVolId,
														LPXLOPER XL_swoptVolId,
														LPXLOPER XL_sabrCapVolId,
														LPXLOPER XL_sabrSwoptVolId,
														LPXLOPER XL_sigmaOrAlpha,
														LPXLOPER XL_flatVolId,
														LPXLOPER XL_convAdjustVolId,
														LPXLOPER XL_convAdjustManagerId,
														LPXLOPER XL_correlDiagCapId,
														LPXLOPER XL_correlDiagSwoptId,
														LPXLOPER XL_correlCorrId,
														LPXLOPER XL_curveModelParamsId,
														LPXLOPER XL_modelData,
														LPXLOPER XL_calibParams,
														LPXLOPER XL_productsToPrice,
														LPXLOPER XL_forexId,
														LPXLOPER XL_fund_basis_curvesId)
{
	long retCode;
	ARM_result C_result;

	static XLOPER XL_result;

	ARM_XL_TRY_BLOCK_BEGIN
	{
		ARM_NOCALCIFWIZ();

		CCString C_csoId;
		CCString C_zcId;
		//vol
		CCString C_capVolId, C_swoptVolId, C_flatVolId;

		CCString C_convAdjustVolId, C_convAdjustManagerId;

		//correl
		CCString C_correlDiagCapId, C_correlDiagSwoptId, C_correlCorrId;

		//basis case:
		CCString C_forexId;
		VECTOR<CCString> C_fund_basis_curvesId(3, "NULL");

		// error
		static int error;
		static char* reason = "";

		XL_readStrCell(	  XL_csoId,	C_csoId, " ARM_ERR: cso id: object expected", C_result);
		XL_readStrCell(   XL_zcId, C_zcId, " ARM_ERR: zc id: object expected", C_result);
		XL_readStrCell(   XL_capVolId, C_capVolId, " ARM_ERR: cap vol id: object expected", C_result);
		XL_readStrCell(   XL_swoptVolId, C_swoptVolId, " ARM_ERR: swopt vol id: object expected",	C_result);
		XL_readStrCellWD(   XL_convAdjustVolId, C_convAdjustVolId, "NULL", " ARM_ERR: convexity adjustment vol id: object expected", C_result);
		XL_readStrCellWD(   XL_convAdjustManagerId, C_convAdjustManagerId, "NULL", " ARM_ERR: convecity adjustment manager id: object expected",	C_result);
		XL_readStrCellWD(   XL_flatVolId, C_flatVolId, "NULL", " ARM_ERR: flat vol id: object expected", C_result);
		XL_readStrCell(   XL_correlDiagCapId, C_correlDiagCapId, " ARM_ERR: correl diag cap id: object expected", C_result);
		XL_readStrCellWD(   XL_correlDiagSwoptId, C_correlDiagSwoptId, "NULL", " ARM_ERR: correl diag cap id: object expected", C_result);
		XL_readStrCellWD(   XL_correlCorrId, C_correlCorrId, "NULL", " ARM_ERR: correl diag corr id: object expected", C_result);
		XL_readStrCellWD(   XL_forexId, C_forexId, "NULL", " ARM_ERR: forex id: object expected", C_result);
		if (C_forexId != CCString("NULL"))
			XL_readStrVector(   XL_fund_basis_curvesId, C_fund_basis_curvesId, " ARM_ERR: funding, domestic basis and funding basis objects expected", DOUBLE_TYPE, C_result);

		//sigma or alpha
		CCString C_SigmaOrAlpha;
		XL_readStrCellWD(XL_sigmaOrAlpha, C_SigmaOrAlpha, "S", " ARM_ERR: SABR vol mode: S, SIGMA, A, ALPHA expected", C_result);
		C_SigmaOrAlpha.toUpper();
        
		long sigmaOrAlpha = 1; // sigma
		if (( C_SigmaOrAlpha == "A" ) || ( C_SigmaOrAlpha == "ALPHA" ))
			sigmaOrAlpha = 0;

		// SABR Cap : Rho, Nu, Beta
		vector<CCString> C_sabrCapVolId;
		vector<CCString> C_sabrCapVolId_def(0);
		XL_readStrVectorWD( XL_sabrCapVolId, C_sabrCapVolId, C_sabrCapVolId_def," ARM_ERR: sabr cap volatilities: vector of objects expected", DOUBLE_TYPE, C_result);
		long rhoCapId = ARM_NULL_OBJECT;
		long nuCapId  = ARM_NULL_OBJECT;
		long betaCapId = ARM_NULL_OBJECT;

		if (C_sabrCapVolId.size() == 2)
		{
			rhoCapId = LocalGetNumObjectId(C_sabrCapVolId[0]);
			nuCapId = LocalGetNumObjectId(C_sabrCapVolId[1]);
		}
		else if (C_sabrCapVolId.size() == 3)
		{
			rhoCapId = LocalGetNumObjectId(C_sabrCapVolId[0]);
			nuCapId = LocalGetNumObjectId(C_sabrCapVolId[1]);
			betaCapId = LocalGetNumObjectId(C_sabrCapVolId[2]);
		}

		// SABR Swopt : Rho, Nu, Beta
		vector<CCString> C_sabrSwoptVolId;
		vector<CCString> C_sabrSwoptVolId_def(0);
		XL_readStrVectorWD( XL_sabrSwoptVolId, C_sabrSwoptVolId, C_sabrSwoptVolId_def," ARM_ERR: sabr swaption volatilities: vector of objects expected", DOUBLE_TYPE, C_result);
		long rhoSwoptId = ARM_NULL_OBJECT;
		long nuSwoptId  = ARM_NULL_OBJECT;
		long betaSwoptId = ARM_NULL_OBJECT;

		if (C_sabrSwoptVolId.size() == 2)
		{
			rhoSwoptId = LocalGetNumObjectId(C_sabrSwoptVolId[0]);
			nuSwoptId = LocalGetNumObjectId(C_sabrSwoptVolId[1]);
		}
		else if (C_sabrSwoptVolId.size() == 3)
		{
			rhoSwoptId = LocalGetNumObjectId(C_sabrSwoptVolId[0]);
			nuSwoptId = LocalGetNumObjectId(C_sabrSwoptVolId[1]);
			betaSwoptId = LocalGetNumObjectId(C_sabrSwoptVolId[2]);
		}

        long convAdjustVolId = (C_convAdjustVolId == "NULL")? ARM_NULL_OBJECT : LocalGetNumObjectId (C_convAdjustVolId);
        long convAdjustManagerId = (C_convAdjustManagerId == "NULL")? ARM_NULL_OBJECT : LocalGetNumObjectId (C_convAdjustManagerId);
        long flatVolId = (C_flatVolId == "NULL")? ARM_NULL_OBJECT : LocalGetNumObjectId (C_flatVolId);
        long correlDiagSwoptId = (C_correlDiagSwoptId == "NULL")? ARM_NULL_OBJECT : LocalGetNumObjectId (C_correlDiagSwoptId);
        long correlCorrId = (C_correlCorrId == "NULL")? ARM_NULL_OBJECT : LocalGetNumObjectId (C_correlCorrId);

		// MODEL PARAMS : in Basis case, only MRS
		vector<CCString> C_curveModelParamsId;
		XL_readStrVector(XL_curveModelParamsId, C_curveModelParamsId, " ARM_ERR: model params id: vector of objects expected", DOUBLE_TYPE, C_result);
		if ((C_curveModelParamsId.size() != 1) && (C_curveModelParamsId.size() != 4))
			throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"model params : 1 or 4 curves expected");

		long mrsId = LocalGetNumObjectId(C_curveModelParamsId[0]);
		long correlId = ARM_NULL_OBJECT;
		long volRatioId = ARM_NULL_OBJECT;
		long mrsSpreadId = ARM_NULL_OBJECT;
		if (C_curveModelParamsId.size() == 4)
		{
			correlId = LocalGetNumObjectId(C_curveModelParamsId[1]);
			volRatioId = LocalGetNumObjectId(C_curveModelParamsId[2]);
			mrsSpreadId = LocalGetNumObjectId(C_curveModelParamsId[3]);
		}

		// MODEL DATA
		VECTOR<double> C_modelData;
		VECTOR<double> C_modelData_def(1, 0.0);
		XL_readNumVectorWD(XL_modelData, C_modelData,C_modelData_def," ARM_ERR: Model Data: array expected",C_result);

		// CALIB PARAMS
		VECTOR<CCString> C_calibParams;
		VECTOR<CCString> C_calibParams_def(3);
		C_calibParams_def[0]= CCString("DIAG"); // diagonal calibration
		C_calibParams_def[1]= CCString("ATM");  // atm swaptions
		C_calibParams_def[2]= CCString("HWM1F");  // hw1f model
		XL_readStrVectorWD(XL_calibParams, C_calibParams,C_calibParams_def," ARM_ERR: Calib Parameters: array expected",DOUBLE_TYPE,C_result);
		if(C_calibParams.size()>8)
		{
			throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"Calib Parameters: calibType,calibStrikeType,modelType,pricingMethod ");
		}
		vector<string> calibParams(C_calibParams.size());
		for (int i = 0; i < C_calibParams.size(); i++)
			calibParams[i] = CCSTringToSTLString(C_calibParams[i]);

		// Products to price
		VECTOR<CCString> C_productsToPrice;
		VECTOR<CCString> C_productsToPriceDef(5, "N");
		XL_readStrVectorWD(XL_productsToPrice, C_productsToPrice, C_productsToPriceDef, " ARM_ERR: control variates: array of string expected",DOUBLE_TYPE,C_result);
		int size = C_productsToPrice.size();
		if (size > 6)
			throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"Products to price flags : invalid size");

		std::deque<bool> productsToPrice(6,false);
		for (i = 0; i < size; i++)
		{
			productsToPrice[i] = (C_productsToPrice[i] == "Y" || C_productsToPrice[i] == "YES");
		}

		// CSO BASIS :
		long forexId = (C_forexId == "NULL") ? ARM_NULL_OBJECT : LocalGetNumObjectId (C_forexId);
        long fundZcId = (C_fund_basis_curvesId[0] == "NULL") ? ARM_NULL_OBJECT : LocalGetNumObjectId (C_fund_basis_curvesId[0]);
        long domBasisZcId = (C_fund_basis_curvesId[1] == "NULL") ? ARM_NULL_OBJECT : LocalGetNumObjectId (C_fund_basis_curvesId[1]);
        long fundBasisZcId = (C_fund_basis_curvesId[2] == "NULL") ? ARM_NULL_OBJECT : LocalGetNumObjectId (C_fund_basis_curvesId[2]);

		CCString curClass = LOCAL_GC_CALLABLE_SO_CLASS;

		long objId;

		CCString stringId = GetLastCurCellEnvValue ();
		
		if(!stringId)
			objId = -1;
		else
		{
			CCString prevClass = LocalGetStringObjectClass (stringId);
			
			objId = LocalGetNumObjectId (stringId);
				
			if(curClass != prevClass)
			{
				FreeCurCellContent ();
				objId = -1;
			}
		}

		retCode = ARMLOCAL_INITCSO (LocalGetNumObjectId(C_csoId),
									LocalGetNumObjectId(C_zcId),
									LocalGetNumObjectId(C_capVolId),
									LocalGetNumObjectId(C_swoptVolId),
									sigmaOrAlpha,
									rhoCapId,
									nuCapId,
									betaCapId,
									rhoSwoptId,
									nuSwoptId,
									betaSwoptId,
									flatVolId,
									convAdjustVolId,
									convAdjustManagerId,
									LocalGetNumObjectId(C_correlDiagCapId),
									correlDiagSwoptId,
									correlCorrId,
									mrsId,
									correlId,
									volRatioId,
									mrsSpreadId,
									C_modelData,
									calibParams,
									productsToPrice,
									forexId,
									fundZcId,
									domBasisZcId,
									fundBasisZcId,
									C_result,
									objId);

		if(retCode == ARM_OK)
		{
			objId = C_result.getLong ();
			LocalSetCurCellEnvValue (curClass, objId); 
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
	ARM_XL_TRY_BLOCK_END

	ARM_XL_CATCH_ARM_EXPT

	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_ARM_INITCSO" )

	return (LPXLOPER)&XL_result;
}



__declspec(dllexport) LPXLOPER WINAPI Local_PXL_ARM_INITCSO(LPXLOPER XL_csoId,
															LPXLOPER XL_zcId,
															LPXLOPER XL_capVolId,
															LPXLOPER XL_swoptVolId,
															LPXLOPER XL_sabrCapVolId,
															LPXLOPER XL_sabrSwoptVolId,
															LPXLOPER XL_sigmaOrAlpha,
															LPXLOPER XL_flatVolId,
															LPXLOPER XL_convAdjustVolId,
															LPXLOPER XL_convAdjustManagerId,
															LPXLOPER XL_correlDiagCapId,
															LPXLOPER XL_correlDiagSwoptId,
															LPXLOPER XL_correlCorrId,
															LPXLOPER XL_curveModelParamsId,
															LPXLOPER XL_modelData,
															LPXLOPER XL_calibParams,
															LPXLOPER XL_productsToPrice,
															LPXLOPER XL_forexId,
															LPXLOPER XL_fund_basis_curvesId)
{
	long retCode;
	ARM_result C_result;

	static XLOPER XL_result;

	ARM_XL_TRY_BLOCK_BEGIN
	{
		ARM_NOCALCIFWIZ();

		CCString C_csoId;
		CCString C_zcId;
		//vol
		CCString C_capVolId, C_swoptVolId, C_flatVolId;

		CCString C_convAdjustVolId, C_convAdjustManagerId;

		//correl
		CCString C_correlDiagCapId, C_correlDiagSwoptId, C_correlCorrId;

		//basis case:
		CCString C_forexId;
		VECTOR<CCString> C_fund_basis_curvesId(3, "NULL");

		// error
		static int error;
		static char* reason = "";

		XL_readStrCell(	  XL_csoId,	C_csoId, " ARM_ERR: cso id: object expected", C_result);
		XL_readStrCell(   XL_zcId, C_zcId, " ARM_ERR: zc id: object expected", C_result);
		XL_readStrCell(   XL_capVolId, C_capVolId, " ARM_ERR: cap vol id: object expected", C_result);
		XL_readStrCell(   XL_swoptVolId, C_swoptVolId, " ARM_ERR: swopt vol id: object expected",	C_result);
		XL_readStrCellWD(   XL_convAdjustVolId, C_convAdjustVolId, "NULL", " ARM_ERR: convexity adjustment vol id: object expected", C_result);
		XL_readStrCellWD(   XL_convAdjustManagerId, C_convAdjustManagerId, "NULL", " ARM_ERR: convecity adjustment manager id: object expected",	C_result);
		XL_readStrCellWD(   XL_flatVolId, C_flatVolId, "NULL", " ARM_ERR: flat vol id: object expected", C_result);
		XL_readStrCell(   XL_correlDiagCapId, C_correlDiagCapId, " ARM_ERR: correl diag cap id: object expected", C_result);
		XL_readStrCellWD(   XL_correlDiagSwoptId, C_correlDiagSwoptId, "NULL", " ARM_ERR: correl diag cap id: object expected", C_result);
		XL_readStrCellWD(   XL_correlCorrId, C_correlCorrId, "NULL", " ARM_ERR: correl diag corr id: object expected", C_result);
		XL_readStrCellWD(   XL_forexId, C_forexId, "NULL", " ARM_ERR: forex id: object expected", C_result);
		if (C_forexId != CCString("NULL"))
			XL_readStrVector(   XL_fund_basis_curvesId, C_fund_basis_curvesId, " ARM_ERR: funding, domestic basis and funding basis objects expected", DOUBLE_TYPE, C_result);

		//sigma or alpha
		CCString C_SigmaOrAlpha;
		XL_readStrCellWD(XL_sigmaOrAlpha, C_SigmaOrAlpha, "S", " ARM_ERR: SABR vol mode: S, SIGMA, A, ALPHA expected", C_result);
		C_SigmaOrAlpha.toUpper();
        
		long sigmaOrAlpha = 1; // sigma
		if (( C_SigmaOrAlpha == "A" ) || ( C_SigmaOrAlpha == "ALPHA" ))
			sigmaOrAlpha = 0;

		// SABR Cap : Rho, Nu, Beta
		vector<CCString> C_sabrCapVolId;
		vector<CCString> C_sabrCapVolId_def(0);
		XL_readStrVectorWD( XL_sabrCapVolId, C_sabrCapVolId, C_sabrCapVolId_def," ARM_ERR: sabr cap volatilities: vector of objects expected", DOUBLE_TYPE, C_result);
		long rhoCapId = ARM_NULL_OBJECT;
		long nuCapId  = ARM_NULL_OBJECT;
		long betaCapId = ARM_NULL_OBJECT;

		if (C_sabrCapVolId.size() == 2)
		{
			rhoCapId = LocalGetNumObjectId(C_sabrCapVolId[0]);
			nuCapId = LocalGetNumObjectId(C_sabrCapVolId[1]);
		}
		else if (C_sabrCapVolId.size() == 3)
		{
			rhoCapId = LocalGetNumObjectId(C_sabrCapVolId[0]);
			nuCapId = LocalGetNumObjectId(C_sabrCapVolId[1]);
			betaCapId = LocalGetNumObjectId(C_sabrCapVolId[2]);
		}

		// SABR Swopt : Rho, Nu, Beta
		vector<CCString> C_sabrSwoptVolId;
		vector<CCString> C_sabrSwoptVolId_def(0);
		XL_readStrVectorWD( XL_sabrSwoptVolId, C_sabrSwoptVolId, C_sabrSwoptVolId_def," ARM_ERR: sabr swaption volatilities: vector of objects expected", DOUBLE_TYPE, C_result);
		long rhoSwoptId = ARM_NULL_OBJECT;
		long nuSwoptId  = ARM_NULL_OBJECT;
		long betaSwoptId = ARM_NULL_OBJECT;

		if (C_sabrSwoptVolId.size() == 2)
		{
			rhoSwoptId = LocalGetNumObjectId(C_sabrSwoptVolId[0]);
			nuSwoptId = LocalGetNumObjectId(C_sabrSwoptVolId[1]);
		}
		else if (C_sabrSwoptVolId.size() == 3)
		{
			rhoSwoptId = LocalGetNumObjectId(C_sabrSwoptVolId[0]);
			nuSwoptId = LocalGetNumObjectId(C_sabrSwoptVolId[1]);
			betaSwoptId = LocalGetNumObjectId(C_sabrSwoptVolId[2]);
		}

        long convAdjustVolId = (C_convAdjustVolId == "NULL")? ARM_NULL_OBJECT : LocalGetNumObjectId (C_convAdjustVolId);
        long convAdjustManagerId = (C_convAdjustManagerId == "NULL")? ARM_NULL_OBJECT : LocalGetNumObjectId (C_convAdjustManagerId);
        long flatVolId = (C_flatVolId == "NULL")? ARM_NULL_OBJECT : LocalGetNumObjectId (C_flatVolId);
        long correlDiagSwoptId = (C_correlDiagSwoptId == "NULL")? ARM_NULL_OBJECT : LocalGetNumObjectId (C_correlDiagSwoptId);
        long correlCorrId = (C_correlCorrId == "NULL")? ARM_NULL_OBJECT : LocalGetNumObjectId (C_correlCorrId);

		// MODEL PARAMS : in Basis case, only MRS
		vector<CCString> C_curveModelParamsId;
		XL_readStrVector(XL_curveModelParamsId, C_curveModelParamsId, " ARM_ERR: model params id: vector of objects expected", DOUBLE_TYPE, C_result);
		if ((C_curveModelParamsId.size() != 1) && (C_curveModelParamsId.size() != 4))
			throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"model params : 1 or 4 curves expected");

		long mrsId = LocalGetNumObjectId(C_curveModelParamsId[0]);
		long correlId = ARM_NULL_OBJECT;
		long volRatioId = ARM_NULL_OBJECT;
		long mrsSpreadId = ARM_NULL_OBJECT;
		if (C_curveModelParamsId.size() == 4)
		{
			correlId = LocalGetNumObjectId(C_curveModelParamsId[1]);
			volRatioId = LocalGetNumObjectId(C_curveModelParamsId[2]);
			mrsSpreadId = LocalGetNumObjectId(C_curveModelParamsId[3]);
		}

		// MODEL DATA
		VECTOR<double> C_modelData;
		VECTOR<double> C_modelData_def(1, 0.0);
		XL_readNumVectorWD(XL_modelData, C_modelData,C_modelData_def," ARM_ERR: Model Data: array expected",C_result);

		// CALIB PARAMS
		VECTOR<CCString> C_calibParams;
		VECTOR<CCString> C_calibParams_def(3);
		C_calibParams_def[0]= CCString("DIAG"); // diagonal calibration
		C_calibParams_def[1]= CCString("ATM");  // atm swaptions
		C_calibParams_def[2]= CCString("HWM1F");  // hw1f model
		XL_readStrVectorWD(XL_calibParams, C_calibParams,C_calibParams_def," ARM_ERR: Calib Parameters: array expected",DOUBLE_TYPE,C_result);
		if(C_calibParams.size()>8)
		{
			throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"Calib Parameters: calibType,calibStrikeType,modelType,pricingMethod ");
		}
		vector<string> calibParams(C_calibParams.size());
		for (int i = 0; i < C_calibParams.size(); i++)
			calibParams[i] = CCSTringToSTLString(C_calibParams[i]);

		// Products to price
		VECTOR<CCString> C_productsToPrice;
		VECTOR<CCString> C_productsToPriceDef(5, "N");
		XL_readStrVectorWD(XL_productsToPrice, C_productsToPrice, C_productsToPriceDef, " ARM_ERR: control variates: array of string expected",DOUBLE_TYPE,C_result);
		int size = C_productsToPrice.size();
		if (size > 6)
			throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"Products to price flags : invalid size");

		std::deque<bool> productsToPrice(6,false);
		for (i = 0; i < size; i++)
		{
			productsToPrice[i] = (C_productsToPrice[i] == "Y" || C_productsToPrice[i] == "YES");
		}

		// CSO BASIS :
		long forexId = (C_forexId == "NULL") ? ARM_NULL_OBJECT : LocalGetNumObjectId (C_forexId);
        long fundZcId = (C_fund_basis_curvesId[0] == "NULL") ? ARM_NULL_OBJECT : LocalGetNumObjectId (C_fund_basis_curvesId[0]);
        long domBasisZcId = (C_fund_basis_curvesId[1] == "NULL") ? ARM_NULL_OBJECT : LocalGetNumObjectId (C_fund_basis_curvesId[1]);
        long fundBasisZcId = (C_fund_basis_curvesId[2] == "NULL") ? ARM_NULL_OBJECT : LocalGetNumObjectId (C_fund_basis_curvesId[2]);

		CCString curClass = LOCAL_GC_CALLABLE_SO_CLASS;

		long objId;
		CCString stringId;

		retCode = ARMLOCAL_INITCSO(LocalGetNumObjectId(C_csoId),
								   LocalGetNumObjectId(C_zcId),
								   LocalGetNumObjectId(C_capVolId),
								   LocalGetNumObjectId(C_swoptVolId),
								   sigmaOrAlpha,
								   rhoCapId,
								   nuCapId,
								   betaCapId,
								   rhoSwoptId,
								   nuSwoptId,
								   betaSwoptId,
								   flatVolId,
								   convAdjustVolId,
								   convAdjustManagerId,
								   LocalGetNumObjectId(C_correlDiagCapId),
								   correlDiagSwoptId,
								   correlCorrId,
								   mrsId,
								   correlId,
								   volRatioId,
								   mrsSpreadId,
								   C_modelData,
								   calibParams,
								   productsToPrice,
								   forexId,
								   fundZcId,
								   domBasisZcId,
								   fundBasisZcId,
								   C_result);

		if ( retCode == ARM_OK )
		{
			objId = C_result.getLong ();
			stringId = LocalMakeObjectId (objId, curClass);
		}

		if ( retCode == ARM_OK )
		{
			XL_result.xltype  = xltypeStr;
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

	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_PXL_ARM_INITCSO" )

	return (LPXLOPER)&XL_result;
}



__declspec(dllexport) LPXLOPER WINAPI Local_ARM_INITCRASPREAD ( LPXLOPER XL_ccsoId,
																LPXLOPER XL_zcId,
																LPXLOPER XL_capVolId,
																LPXLOPER XL_capSabrId,
																LPXLOPER XL_swoptVolId,
																LPXLOPER XL_swoptSabrId,
																LPXLOPER XL_sigmaOrAlpha,
																LPXLOPER XL_convAdjustVolCapId,
																LPXLOPER XL_convAdjustVolSwoptId,
																LPXLOPER XL_convAdjustType,
																LPXLOPER XL_correlCorrId,
																LPXLOPER XL_correlDiagCapId,
																LPXLOPER XL_correlDiagSwoptId,
																LPXLOPER XL_mrs,
																LPXLOPER XL_volRatio,
																LPXLOPER XL_mrsSpread,
																LPXLOPER XL_correl,
																LPXLOPER XL_modelParams,
																LPXLOPER XL_productsToPrice,
																LPXLOPER XL_localCalibFlags)
{
	long retCode;
	ARM_result C_result;

	static XLOPER XL_result;

	ARM_XL_TRY_BLOCK_BEGIN
	{
		ARM_NOCALCIFWIZ();

		int i=0;
		CCString C_ccsoId;
		CCString C_zcId;
		//vol
		CCString C_capVolId, C_swoptVolId, C_flatVolId;
		vector<CCString> C_capSabrId;
		vector<CCString> C_capSabrId_def(0);

		vector<CCString> C_swoptSabrId;
		vector<CCString> C_swoptSabrId_def(0);

		CCString C_SigmaOrAlpha;

		CCString C_convAdjustVolCapId, C_convAdjustVolSwoptId;
		long convAdjustVolCapId, convAdjustVolSwoptId;

		CCString C_convAdjustTypeStr;
		long convAdjustType;

		//correl
		CCString C_correlDiagCapId, C_correlDiagSwoptId, C_correlCorrId;
		long correlDiagCapId, correlDiagSwoptId, correlCorrId;

		//model parameters
		double C_mrs, C_correl, C_volRatio, C_mrsSpread;

		// error
		static int error;
		static char* reason = "";

		XL_readStrCell(	  XL_ccsoId, C_ccsoId, " ARM_ERR: cra spread id: object expected", C_result);
		XL_readStrCell(   XL_zcId, C_zcId, " ARM_ERR: zc id: object expected", C_result);
		XL_readStrCell(   XL_capVolId, C_capVolId, " ARM_ERR: cap vol id: object expected", C_result);
		XL_readStrVectorWD( XL_capSabrId, C_capSabrId, C_capSabrId_def," ARM_ERR: cap sabr volatilities: vector of objects expected", DOUBLE_TYPE, C_result);
		XL_readStrCell(   XL_swoptVolId, C_swoptVolId, " ARM_ERR: swopt vol id: object expected",	C_result);
		XL_readStrVectorWD( XL_swoptSabrId, C_swoptSabrId, C_swoptSabrId_def," ARM_ERR: sabr volatilities: vector of objects expected", DOUBLE_TYPE, C_result);
		XL_readStrCellWD(   XL_convAdjustVolCapId, C_convAdjustVolCapId, "NULL", " ARM_ERR: convexity adjustment vol (cap) id: object expected", C_result);
		XL_readStrCellWD(   XL_convAdjustVolSwoptId, C_convAdjustVolSwoptId, "NULL", " ARM_ERR: convexity adjustment vol (swaption) id: object expected", C_result);
		XL_readStrCellWD(	XL_convAdjustType,C_convAdjustTypeStr,"SUMEXP"," ARM_ERR: ConvAdjustType: string expected", C_result);
		XL_readStrCellWD(   XL_correlDiagCapId, C_correlDiagCapId, "NULL", " ARM_ERR: correl diag cap id: object expected", C_result);
		XL_readStrCellWD(   XL_correlDiagSwoptId, C_correlDiagSwoptId, "NULL", " ARM_ERR: correl diag swaption id: object expected", C_result);
		XL_readStrCellWD(   XL_correlCorrId, C_correlCorrId, "NULL", " ARM_ERR: correl corr id: object expected", C_result);
		XL_readNumCell(   XL_mrs, C_mrs, " ARM_ERR: mrs: double expected", C_result);
		XL_readNumCell(   XL_volRatio, C_volRatio, " ARM_ERR: vol ratio: double expected", C_result);
		XL_readNumCell(   XL_mrsSpread, C_mrsSpread, " ARM_ERR: mrs spread: double expected", C_result);
		XL_readNumCell(   XL_correl, C_correl, " ARM_ERR: correl: double expected", C_result);
		XL_readStrCellWD(	XL_sigmaOrAlpha, C_SigmaOrAlpha, "S", " ARM_ERR: SABR vol mode: S, SIGMA, A, ALPHA expected", C_result);

		long rhoCapId = ARM_NULL_OBJECT;
		long nuCapId = ARM_NULL_OBJECT;
		long betaCapId = ARM_NULL_OBJECT;
		
		if (C_capSabrId.size() == 2)
		{
			rhoCapId = LocalGetNumObjectId(C_capSabrId[0]);
			nuCapId = LocalGetNumObjectId(C_capSabrId[1]);
			betaCapId = ARM_NULL_OBJECT;
		}
		else if (C_capSabrId.size() == 3)
		{
			rhoCapId = LocalGetNumObjectId(C_capSabrId[0]);
			nuCapId = LocalGetNumObjectId(C_capSabrId[1]);
			betaCapId = LocalGetNumObjectId(C_capSabrId[2]);
		}

		long rhoSwoptId = ARM_NULL_OBJECT;
		long nuSwoptId = ARM_NULL_OBJECT;
		long betaSwoptId = ARM_NULL_OBJECT;

		if (C_swoptSabrId.size() == 0)
		{
			rhoSwoptId = ARM_NULL_OBJECT;
			nuSwoptId = ARM_NULL_OBJECT;
			betaSwoptId = ARM_NULL_OBJECT;
		}
		else if (C_swoptSabrId.size() == 2)
		{
			rhoSwoptId = LocalGetNumObjectId(C_swoptSabrId[0]);
			nuSwoptId = LocalGetNumObjectId(C_swoptSabrId[1]);
			betaSwoptId = ARM_NULL_OBJECT;
		}
		else if (C_swoptSabrId.size() == 3)
		{
			rhoSwoptId = LocalGetNumObjectId(C_swoptSabrId[0]);
			nuSwoptId = LocalGetNumObjectId(C_swoptSabrId[1]);
			betaSwoptId = LocalGetNumObjectId(C_swoptSabrId[2]);
		}

		C_SigmaOrAlpha.toUpper();
		long sigmaOrAlpha = 1; // sigma
		if (( C_SigmaOrAlpha == "A" ) || ( C_SigmaOrAlpha == "ALPHA" ))
			sigmaOrAlpha = 0;

        convAdjustVolCapId = (C_convAdjustVolCapId == "NULL")? ARM_NULL_OBJECT : LocalGetNumObjectId (C_convAdjustVolCapId);
		convAdjustVolSwoptId = (C_convAdjustVolSwoptId == "NULL")? ARM_NULL_OBJECT : LocalGetNumObjectId (C_convAdjustVolSwoptId);
		if ( (convAdjustType = ARM_ConvSummitFormulae (C_convAdjustTypeStr, C_result)) == ARM_DEFAULT_ERR)
		{
			throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"invalid convexity adjustment type");
		}
        correlDiagCapId = (C_correlDiagCapId == "NULL")? ARM_NULL_OBJECT : LocalGetNumObjectId (C_correlDiagCapId);
        correlDiagSwoptId = (C_correlDiagSwoptId == "NULL")? ARM_NULL_OBJECT : LocalGetNumObjectId (C_correlDiagSwoptId);
        correlCorrId = (C_correlCorrId == "NULL")? ARM_NULL_OBJECT : LocalGetNumObjectId (C_correlCorrId);

		// MODEL PARAMS
		VECTOR<CCString> C_modelParams;
		XL_readStrVector(XL_modelParams, C_modelParams," ARM_ERR: Model Parameters: array expected",DOUBLE_TYPE,C_result);
		vector<string> modelParams(C_modelParams.size());
		for (i = 0; i < C_modelParams.size(); i++)
			modelParams[i] = CCSTringToSTLString(C_modelParams[i]);

		// Products to price
		VECTOR<CCString> C_productsToPrice;
		VECTOR<CCString> C_productsToPriceDef(5, "N");
		XL_readStrVectorWD(XL_productsToPrice, C_productsToPrice, C_productsToPriceDef, " ARM_ERR: control variates: array of string expected",DOUBLE_TYPE,C_result);
		int size = C_productsToPrice.size();
		if (size > 5)
			throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"Products to price flags : invalid size");

		std::deque<bool> productsToPrice(size);
		for (i = 0; i < size; i++)
			productsToPrice[i] = (C_productsToPrice[i] == "Y" || C_productsToPrice[i] == "YES");

		CCString curClass = LOCAL_GC_CCSO_CLASS;

		long objId;

		CCString stringId = GetLastCurCellEnvValue ();
		
		if (!stringId)
		   objId = -1;
		else
		{
			CCString prevClass = LocalGetStringObjectClass (stringId);
			
			objId = LocalGetNumObjectId (stringId);
				
			if ( curClass != prevClass )
			{
				FreeCurCellContent ();
				objId = -1;
			}
		}

		/// Local calib flags
		VECTOR<CCString> LocalCalibFlagsXL;
		VECTOR<CCString> LocalCalibFlagsDefault(4,"N");
		XL_readStrVectorWD(XL_localCalibFlags,LocalCalibFlagsXL,LocalCalibFlagsDefault," ARM_ERR: Local Calib Flag: array of string expected",DOUBLE_TYPE,C_result);
		vector< string > localCalibFlags(5,"N");

		/// Vol Unsqueezer
		localCalibFlags[1]=CCSTringToSTLString(LocalCalibFlagsXL[0]);

		/// PV adjuster
		localCalibFlags[2]=CCSTringToSTLString(LocalCalibFlagsXL[1]);

		/// Bootstrap Optimizer
		localCalibFlags[3]=CCSTringToSTLString(LocalCalibFlagsXL[2]);

		/// Switch to old calib
		localCalibFlags[4]=CCSTringToSTLString(LocalCalibFlagsXL[3]);
		
		retCode = ARMLOCAL_INITCRASPREAD(LocalGetNumObjectId(C_ccsoId),
										LocalGetNumObjectId(C_zcId),
										LocalGetNumObjectId(C_capVolId),
										rhoCapId,
										nuCapId,
										betaCapId,
										LocalGetNumObjectId(C_swoptVolId),
										rhoSwoptId,
										nuSwoptId,
										betaSwoptId,
										sigmaOrAlpha,
										convAdjustVolCapId,
										convAdjustVolSwoptId,
										convAdjustType,
										correlCorrId,
										correlDiagCapId,
										correlDiagSwoptId,
										C_mrs,
										C_volRatio,
										C_mrsSpread,
										C_correl,
										modelParams,
										productsToPrice,
										localCalibFlags,
										C_result,
										objId);

		if(retCode == ARM_OK)
		{
			objId = C_result.getLong ();
			LocalSetCurCellEnvValue (curClass, objId); 
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
	ARM_XL_TRY_BLOCK_END

	ARM_XL_CATCH_ARM_EXPT

	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_ARM_INITCRASPREAD" )

	return (LPXLOPER)&XL_result;
}



__declspec(dllexport) LPXLOPER WINAPI Local_PXL_ARM_INITCRASPREAD ( LPXLOPER XL_ccsoId,
																	LPXLOPER XL_zcId,
																	LPXLOPER XL_capVolId,
																	LPXLOPER XL_capSabrId,
																	LPXLOPER XL_swoptVolId,
																	LPXLOPER XL_swoptSabrId,
																	LPXLOPER XL_sigmaOrAlpha,
																	LPXLOPER XL_convAdjustVolCapId,
																	LPXLOPER XL_convAdjustVolSwoptId,
																	LPXLOPER XL_convAdjustType,
																	LPXLOPER XL_correlCorrId,
																	LPXLOPER XL_correlDiagCapId,
																	LPXLOPER XL_correlDiagSwoptId,
																	LPXLOPER XL_mrs,
																	LPXLOPER XL_correl,
																	LPXLOPER XL_volRatio,
																	LPXLOPER XL_mrsSpread,
																	LPXLOPER XL_modelParams,
																	LPXLOPER XL_productsToPrice,
																	LPXLOPER XL_localCalibFlags)
{
	long retCode;
	ARM_result C_result;

	static XLOPER XL_result;

	ARM_XL_TRY_BLOCK_BEGIN
	{
		ARM_NOCALCIFWIZ();

		int i=0;
		CCString C_ccsoId;
		CCString C_zcId;
		//vol
		CCString C_capVolId, C_swoptVolId, C_flatVolId;
		vector<CCString> C_capSabrId;
		vector<CCString> C_capSabrId_def(0);

		vector<CCString> C_swoptSabrId;
		vector<CCString> C_swoptSabrId_def(0);

		CCString C_SigmaOrAlpha;

		CCString C_convAdjustVolCapId, C_convAdjustVolSwoptId;
		long convAdjustVolCapId, convAdjustVolSwoptId;

		CCString C_convAdjustTypeStr;
		long convAdjustType;

		//correl
		CCString C_correlDiagCapId, C_correlDiagSwoptId, C_correlCorrId;
		long correlDiagCapId, correlDiagSwoptId, correlCorrId;

		//model parameters
		double C_mrs, C_correl, C_volRatio, C_mrsSpread;

		// error
		static int error;
		static char* reason = "";

		XL_readStrCell(	  XL_ccsoId, C_ccsoId, " ARM_ERR: cra spread id: object expected", C_result);
		XL_readStrCell(   XL_zcId, C_zcId, " ARM_ERR: zc id: object expected", C_result);
		XL_readStrCell(   XL_capVolId, C_capVolId, " ARM_ERR: cap vol id: object expected", C_result);
		XL_readStrVectorWD( XL_capSabrId, C_capSabrId, C_capSabrId_def," ARM_ERR: cap sabr volatilities: vector of objects expected", DOUBLE_TYPE, C_result);
		XL_readStrCell(   XL_swoptVolId, C_swoptVolId, " ARM_ERR: swopt vol id: object expected",	C_result);
		XL_readStrVectorWD( XL_swoptSabrId, C_swoptSabrId, C_swoptSabrId_def," ARM_ERR: sabr volatilities: vector of objects expected", DOUBLE_TYPE, C_result);
		XL_readStrCellWD(   XL_convAdjustVolCapId, C_convAdjustVolCapId, "NULL", " ARM_ERR: convexity adjustment vol (cap) id: object expected", C_result);
		XL_readStrCellWD(   XL_convAdjustVolSwoptId, C_convAdjustVolSwoptId, "NULL", " ARM_ERR: convexity adjustment vol (swaption) id: object expected", C_result);
		XL_readStrCellWD(	XL_convAdjustType,C_convAdjustTypeStr,"SUMEXP"," ARM_ERR: ConvAdjustType: string expected", C_result);
		XL_readStrCellWD(   XL_correlDiagCapId, C_correlDiagCapId, "NULL", " ARM_ERR: correl diag cap id: object expected", C_result);
		XL_readStrCellWD(   XL_correlDiagSwoptId, C_correlDiagSwoptId, "NULL", " ARM_ERR: correl diag cap id: object expected", C_result);
		XL_readStrCellWD(   XL_correlCorrId, C_correlCorrId, "NULL", " ARM_ERR: correl diag corr id: object expected", C_result);
		XL_readNumCell(   XL_mrs, C_mrs, " ARM_ERR: mrs: double expected", C_result);
		XL_readNumCell(   XL_correl, C_correl, " ARM_ERR: correl: double expected", C_result);
		XL_readNumCell(   XL_volRatio, C_volRatio, " ARM_ERR: vol ratio: double expected", C_result);
		XL_readNumCell(   XL_mrsSpread, C_mrsSpread, " ARM_ERR: mrs spread: double expected", C_result);
		XL_readStrCellWD(XL_sigmaOrAlpha, C_SigmaOrAlpha, "S", " ARM_ERR: SABR vol mode: S, SIGMA, A, ALPHA expected", C_result);

		long rhoCapId = ARM_NULL_OBJECT;
		long nuCapId = ARM_NULL_OBJECT;
		long betaCapId = ARM_NULL_OBJECT;
		
		if (C_capSabrId.size() == 2)
		{
			rhoCapId = LocalGetNumObjectId(C_capSabrId[0]);
			nuCapId = LocalGetNumObjectId(C_capSabrId[1]);
			betaCapId = ARM_NULL_OBJECT;
		}
		else if (C_capSabrId.size() == 3)
		{
			rhoCapId = LocalGetNumObjectId(C_capSabrId[0]);
			nuCapId = LocalGetNumObjectId(C_capSabrId[1]);
			betaCapId = LocalGetNumObjectId(C_capSabrId[2]);
		}

		long rhoSwoptId = ARM_NULL_OBJECT;
		long nuSwoptId = ARM_NULL_OBJECT;
		long betaSwoptId = ARM_NULL_OBJECT;
		
		if (C_swoptSabrId.size() == 2)
		{
			rhoSwoptId = LocalGetNumObjectId(C_swoptSabrId[0]);
			nuSwoptId = LocalGetNumObjectId(C_swoptSabrId[1]);
			betaSwoptId = ARM_NULL_OBJECT;
		}
		else if (C_swoptSabrId.size() == 3)
		{
			rhoSwoptId = LocalGetNumObjectId(C_swoptSabrId[0]);
			nuSwoptId = LocalGetNumObjectId(C_swoptSabrId[1]);
			betaSwoptId = LocalGetNumObjectId(C_swoptSabrId[2]);
		}

		C_SigmaOrAlpha.toUpper();        
		long sigmaOrAlpha = 1; // sigma
		if (( C_SigmaOrAlpha == "A" ) || ( C_SigmaOrAlpha == "ALPHA" ))
			sigmaOrAlpha = 0;

        convAdjustVolCapId = (C_convAdjustVolCapId == "NULL")? ARM_NULL_OBJECT : LocalGetNumObjectId (C_convAdjustVolCapId);
		convAdjustVolSwoptId = (C_convAdjustVolSwoptId == "NULL")? ARM_NULL_OBJECT : LocalGetNumObjectId (C_convAdjustVolSwoptId);
		if ( (convAdjustType = ARM_ConvSummitFormulae (C_convAdjustTypeStr, C_result)) == ARM_DEFAULT_ERR)
		{
			throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"invalid convexity adjustment type");
		}
        correlDiagCapId = (C_correlDiagCapId == "NULL")? ARM_NULL_OBJECT : LocalGetNumObjectId (C_correlDiagCapId);
        correlDiagSwoptId = (C_correlDiagSwoptId == "NULL")? ARM_NULL_OBJECT : LocalGetNumObjectId (C_correlDiagSwoptId);
        correlCorrId = (C_correlCorrId == "NULL")? ARM_NULL_OBJECT : LocalGetNumObjectId (C_correlCorrId);

		// MODEL PARAMS
		VECTOR<CCString> C_modelParams;
		XL_readStrVector(XL_modelParams, C_modelParams," ARM_ERR: Model Parameters: array expected",DOUBLE_TYPE,C_result);
		vector<string> modelParams(C_modelParams.size());
		for (i = 0; i < C_modelParams.size(); i++)
			modelParams[i] = CCSTringToSTLString(C_modelParams[i]);

		// Products to price
		VECTOR<CCString> C_productsToPrice;
		VECTOR<CCString> C_productsToPriceDef(5, "N");
		XL_readStrVectorWD(XL_productsToPrice, C_productsToPrice, C_productsToPriceDef, " ARM_ERR: control variates: array of string expected",DOUBLE_TYPE,C_result);
		int size = C_productsToPrice.size();
		if (size > 5)
			throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"Products to price flags : invalid size");

		std::deque<bool> productsToPrice(size);
		for (i = 0; i < size; i++)
			productsToPrice[i] = (C_productsToPrice[i] == "Y" || C_productsToPrice[i] == "YES");

		CCString curClass = LOCAL_GC_CCSO_CLASS;

		long objId;

		CCString stringId;

		/// Local calib flags
		VECTOR<CCString> LocalCalibFlagsXL;
		VECTOR<CCString> LocalCalibFlagsDefault(4,"N");
		XL_readStrVectorWD(XL_localCalibFlags,LocalCalibFlagsXL,LocalCalibFlagsDefault," ARM_ERR: Local Calib Flag: array of string expected",DOUBLE_TYPE,C_result);
		vector< string > localCalibFlags(5,"N");

		/// Vol Unsqueezer
		localCalibFlags[1]=CCSTringToSTLString(LocalCalibFlagsXL[0]);

		/// PV adjuster
		localCalibFlags[2]=CCSTringToSTLString(LocalCalibFlagsXL[1]);

		/// Bootstrap Optimizer
		localCalibFlags[3]=CCSTringToSTLString(LocalCalibFlagsXL[2]);

		/// Switch to old calib
		localCalibFlags[4]=CCSTringToSTLString(LocalCalibFlagsXL[3]);

		retCode = ARMLOCAL_INITCRASPREAD(LocalGetNumObjectId(C_ccsoId),
										LocalGetNumObjectId(C_zcId),
										LocalGetNumObjectId(C_capVolId),
										rhoCapId,
										nuCapId,
										betaCapId,
										LocalGetNumObjectId(C_swoptVolId),
										rhoSwoptId,
										nuSwoptId,
										betaSwoptId,
										sigmaOrAlpha,
										convAdjustVolCapId,
										convAdjustVolSwoptId,
										convAdjustType,
										correlCorrId,
										correlDiagCapId,
										correlDiagSwoptId,
										C_mrs,
										C_volRatio,
										C_mrsSpread,
										C_correl,
										modelParams,
										productsToPrice,
										localCalibFlags,
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
	}
	ARM_XL_TRY_BLOCK_END

	ARM_XL_CATCH_ARM_EXPT

	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_PXL_ARM_INITCRASPREAD" )

	return (LPXLOPER)&XL_result;
}


__declspec(dllexport) LPXLOPER WINAPI Local_InitGlobalCap(  LPXLOPER XL_GCCalculator,
															LPXLOPER XL_zc,
															LPXLOPER XL_capVol,
															LPXLOPER XL_rhoCap,
															LPXLOPER XL_nuCap,
															LPXLOPER XL_betaCap,
															LPXLOPER XL_mrs,
															LPXLOPER XL_calibParams,
															LPXLOPER XL_modelParams,
															LPXLOPER XL_productsToPrice )
{
	long retCode;
	ARM_result C_result;

	static XLOPER XL_result;

	ARM_XL_TRY_BLOCK_BEGIN
	{
		ARM_NOCALCIFWIZ();

		int i=0;
		CCString C_calculatorId;
		long calculatorId;

		CCString C_zcId;
		long zcId;

		//vol
		CCString C_capVolId, C_rhoId, C_nuId, C_betaId;
		long capVolId, rhoId, nuId, betaId;
		double C_mrs;

		// error
		static int error;
		static char* reason = "";

		XL_readStrCell(	  XL_GCCalculator, C_calculatorId, " ARM_ERR: global cap id: object expected", C_result);
		XL_readStrCell(   XL_zc, C_zcId, " ARM_ERR: zc id: object expected", C_result);
		XL_readStrCell(   XL_capVol, C_capVolId, " ARM_ERR: cap vol id: object expected", C_result);
		XL_readStrCellWD(   XL_rhoCap, C_rhoId, "NULL", " ARM_ERR: rho vol id: object expected", C_result);
		XL_readStrCellWD(   XL_nuCap, C_nuId, "NULL", " ARM_ERR: nu vol id: object expected", C_result);
		XL_readStrCellWD(   XL_betaCap, C_betaId, "NULL", " ARM_ERR: beta vol id: object expected", C_result);
		XL_readNumCell(   XL_mrs, C_mrs, " ARM_ERR: mrs: double expected", C_result);

		calculatorId = LocalGetNumObjectId(C_calculatorId);
		zcId = LocalGetNumObjectId(C_zcId);
		capVolId = LocalGetNumObjectId(C_capVolId);
		
		rhoId = (C_rhoId == "NULL") ? ARM_NULL_OBJECT : LocalGetNumObjectId(C_rhoId);
		nuId = (C_nuId == "NULL") ? ARM_NULL_OBJECT : LocalGetNumObjectId(C_nuId);
		betaId = (C_betaId == "NULL") ? ARM_NULL_OBJECT : LocalGetNumObjectId(C_betaId);

		//Calib Params
		//-------------
		VECTOR<double> calibParamsVect;
		XL_readNumVector (XL_calibParams, calibParamsVect, "ARM_ERR: Calib Params: array of double expected", C_result);
		if (calibParamsVect.size() != 5)
		{
			throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"Calib Params: 5 parameters expected");
		}
		ARM_Vector* calibParams = CreateARMVectorFromVECTOR(calibParamsVect);

		//Model Params
		//-------------
		VECTOR<CCString> modelParams;
		XL_readStrVector (XL_modelParams, modelParams, "ARM_ERR: Model Params: array of string expected", DOUBLE_TYPE, C_result);
		if ((modelParams.size() != 5) && (modelParams.size() != 9))
		{
			throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"Model Params: 5 or 9 parameters expected");
		}

		//NbSteps
		CCString nbStepsStr = modelParams[0];
		string sNbStepsStr = CCSTringToSTLString(nbStepsStr);
		int nbSteps = atoi(sNbStepsStr.c_str());

		// Random Generator
		vector<string> randGenerator;

		//Generator Type
		CCString generatorTypeStr = modelParams[1];
		string generatorType = CCSTringToSTLString(generatorTypeStr);
		randGenerator.push_back(generatorType);

		//Inversion Method
		CCString inversionMethodStr = modelParams[2];
		string inversionMethod = CCSTringToSTLString(inversionMethodStr);
		randGenerator.push_back(inversionMethod);

		//Transform Algo
		CCString algoStr = modelParams[3];
		string algo = CCSTringToSTLString(algoStr);
		randGenerator.push_back(algo);

		int samplerType;
		if (modelParams.size() == 9)
		{
			//Generator Type
			CCString generatorTypeStr = modelParams[4];
			string generatorType = CCSTringToSTLString(generatorTypeStr);
			randGenerator.push_back(generatorType);

			//Inversion Method
			CCString inversionMethodStr = modelParams[5];
			string inversionMethod = CCSTringToSTLString(inversionMethodStr);
			randGenerator.push_back(inversionMethod);

			//Transform Algo
			CCString algoStr = modelParams[6];
			string algo = CCSTringToSTLString(algoStr);
			randGenerator.push_back(algo);

			// Nb First Dim
			CCString NbFirstDimStr = modelParams[7];
			string nbFirstDim = CCSTringToSTLString(NbFirstDimStr);
			randGenerator.push_back(nbFirstDim);

			//Sampler Type
			CCString samplerTypeStr = modelParams[8];
			if( (samplerType = ARM_ConvGPSamplerType( samplerTypeStr, C_result)) == ARM_DEFAULT_ERR )
			{
				throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"Model Params: invalid sampler type");
			}
		}
		else
		{
			//Sampler Type
			CCString samplerTypeStr = modelParams[4];
			if( (samplerType = ARM_ConvGPSamplerType( samplerTypeStr, C_result)) == ARM_DEFAULT_ERR )
			{
				throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"Model Params: invalid sampler type");
			}
		}

		// Products to price
		VECTOR<CCString> C_productsToPrice;
		VECTOR<CCString> C_productsToPriceDef(5, "N");
		XL_readStrVectorWD(XL_productsToPrice, C_productsToPrice, C_productsToPriceDef, " ARM_ERR: control variates: array of string expected",DOUBLE_TYPE,C_result);
		int size = C_productsToPrice.size();
		if (size > 5)
			throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"Products to price flags : invalid size");

		std::deque<bool> productsToPrice(size);
		for (i = 0; i < size; i++)
			productsToPrice[i] = (C_productsToPrice[i] == "Y" || C_productsToPrice[i] == "YES");

		CCString curClass = LOCAL_GC_GLOBALCAP_CLASS;

		long objId;

		CCString stringId = GetLastCurCellEnvValue ();
		
		if(!stringId)
			objId = -1;
		else
		{
			CCString prevClass = LocalGetStringObjectClass (stringId);
			
			objId = LocalGetNumObjectId (stringId);
				
			if(curClass != prevClass)
			{
				FreeCurCellContent ();
				objId = -1;
			}
		}

		retCode = ARMLOCAL_INITGLOBALCAP(	calculatorId,
											zcId,
											capVolId,
											rhoId,
											nuId,
											betaId,
											C_mrs,
											calibParams,
											nbSteps,
											randGenerator,
											samplerType,
											productsToPrice,
											C_result,
											objId);

		if(retCode == ARM_OK)
		{
			objId = C_result.getLong ();
			LocalSetCurCellEnvValue (curClass, objId); 
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
	ARM_XL_TRY_BLOCK_END

	ARM_XL_CATCH_ARM_EXPT

	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_ARM_INITGLOBALCAP" )

	return (LPXLOPER)&XL_result;
}


__declspec(dllexport) LPXLOPER WINAPI Local_ARM_EtkConnect(LPXLOPER XL_username,
														   LPXLOPER XL_passwd,
														   LPXLOPER XL_context,
														   LPXLOPER XL_itconfigdomainsdir,
														   LPXLOPER XL_itdomainname)
{
	ARM_result C_result;

	// return
	static XLOPER XL_result;

	//	ARM_BEGIN();
	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();

	// C variable
	CCString C_username;
	CCString C_passwd;
	CCString C_context;
	CCString C_itconfigdomainsdir;
	CCString C_itdomainname;

	// error
	static int error;
	static char* reason = "";

	XL_readStrCell(XL_username,C_username," ARM_ERR: username: string expected",C_result);
	XL_readStrCell(XL_passwd,C_passwd," ARM_ERR: passwd: string expected",C_result);
	XL_readStrCell(XL_context,C_context," ARM_ERR: context: string expected",C_result);
	XL_readStrCell(XL_itconfigdomainsdir,C_itconfigdomainsdir," ARM_ERR: itconfigdomainsdir: string expected",C_result);
	XL_readStrCell(XL_itdomainname,C_itdomainname," ARM_ERR: itdomainname: string expected",C_result);

	switchToETK();

	connection_etoolkit(C_username,
						C_passwd,
						C_context,
						C_itconfigdomainsdir,
						C_itdomainname);


	FreeCurCellErr ();
	XL_result.xltype = xltypeNum;
	XL_result.val.num = 1;
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_ARM_EtkConnect" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


__declspec(dllexport) LPXLOPER WINAPI Local_ARM_INITTARNFX( LPXLOPER XL_tarnId,
															LPXLOPER XL_zcCurvesId,
															LPXLOPER XL_basisCurvesId,
															LPXLOPER XL_forexId,
															LPXLOPER XL_ATMVolId,
															LPXLOPER XL_bsFxVolId,
															LPXLOPER XL_mixtureParamsId,
															LPXLOPER XL_mrsParamsId,
															LPXLOPER XL_QParamsId,
															LPXLOPER XL_correlMatrixId,
															LPXLOPER XL_MCParams,
															LPXLOPER XL_modelType,
															LPXLOPER XL_nbFactor,
															LPXLOPER XL_PDEParams,
															LPXLOPER XL_rescalling,
															LPXLOPER XL_smileFlag,
															LPXLOPER XL_mixCalib,
															LPXLOPER XL_oneFactorFlag,
															LPXLOPER XL_correlType )
{
	long retCode;
	ARM_result C_result;

	static XLOPER XL_result;

	ARM_XL_TRY_BLOCK_BEGIN
	{
		ARM_NOCALCIFWIZ();

		// C variable
		CCString C_tarnId;
		vector<CCString> C_zcCurvesId;
		vector<CCString> C_basisCurvesId;
		vector<CCString> C_forexId;
		vector<CCString> C_ATMVolId;
		vector<CCString> C_bsFxVolId;
		vector<CCString> C_bsFxVolId_default(0);
		vector<CCString> C_mixtureParamsId;
		vector<CCString> C_mixtureParamsId_default(0);
		vector<CCString> C_mrsParamsId;
		vector<CCString> C_QParamsId;
		vector<CCString> C_MCParams;
		CCString C_correlMatrixId;
		vector<CCString> C_PDEParams;
		vector<CCString> C_PDEParams_default(4);
		C_PDEParams_default[0] = "1000.0";
		C_PDEParams_default[1] = "901.0";
		C_PDEParams_default[2] = "6.0";
		C_PDEParams_default[3] = "N";
		double C_nbFactor;
		double C_nbFactor_default = 6.0;
		CCString C_rescalling, C_modelType, C_correlType, C_smileFlag, C_mixCalib, C_oneFactorFlag;

        // error
		static int error;
		static char* reason = "";

		XL_readStrCell(	  XL_tarnId,	 C_tarnId,					" ARM_ERR: Tarn id: object expected",		C_result);
		XL_readStrVector( XL_zcCurvesId, C_zcCurvesId,				" ARM_ERR: zc id: object expected",	DOUBLE_TYPE,	C_result);
		XL_readStrVector( XL_basisCurvesId, C_basisCurvesId,		" ARM_ERR: basis id: object expected", DOUBLE_TYPE,		C_result);
		XL_readStrVector( XL_forexId, C_forexId,					" ARM_ERR: fx id: object expected",	DOUBLE_TYPE,	C_result);
		XL_readStrVector( XL_ATMVolId, C_ATMVolId,					" ARM_ERR: ATM vol id: object expected", DOUBLE_TYPE,		C_result);
		XL_readStrVectorWD( XL_bsFxVolId, C_bsFxVolId, C_bsFxVolId_default,	" ARM_ERR: fx vol id: object expected",	DOUBLE_TYPE,	C_result);
		XL_readStrVectorWD( XL_mixtureParamsId, C_mixtureParamsId, C_mixtureParamsId_default, " ARM_ERR: mixture vol params: object expected", DOUBLE_TYPE,		C_result);

		if ((C_bsFxVolId.size() == 0) && (C_mixtureParamsId.size() == 0))
			throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"must provide BS fxVol or Mixture");

		XL_readStrVector( XL_mrsParamsId, C_mrsParamsId,			" ARM_ERR: mrs params id: object expected",	DOUBLE_TYPE,	C_result);
		XL_readStrVector( XL_QParamsId, C_QParamsId,				" ARM_ERR: Q params id: object expected", DOUBLE_TYPE,		C_result);
		XL_readStrVector( XL_MCParams, C_MCParams,					" ARM_ERR: MC params: vector of string expected", DOUBLE_TYPE,		C_result);
		XL_readStrCell(	  XL_correlMatrixId, C_correlMatrixId,		" ARM_ERR: correl matrix id: object expected",		C_result);
		XL_readNumCellWD( XL_nbFactor, C_nbFactor, C_nbFactor_default, " ARM_ERR: nb factors: integer expected",		C_result);
		XL_readStrVectorWD( XL_PDEParams, C_PDEParams, C_PDEParams_default, " ARM_ERR: PDE parameters: vector expected", DOUBLE_TYPE,		C_result);
		XL_readStrCellWD( XL_rescalling, C_rescalling, "N",		" ARM_ERR: rescalling: Y/N expected",		C_result);
		XL_readStrCell(	  XL_modelType, C_modelType,				" ARM_ERR: model type: string expected",		C_result);
		XL_readStrCellWD( XL_smileFlag, C_smileFlag, "Y",			" ARM_ERR: smile: Y/N expected",		C_result);
		XL_readStrCellWD( XL_mixCalib, C_mixCalib, "N",			" ARM_ERR: mix calib: Y/N expected",		C_result);
		XL_readStrCellWD( XL_oneFactorFlag, C_oneFactorFlag, "N",	" ARM_ERR: one factor flag: Y/N expected",		C_result);
		XL_readStrCellWD( XL_correlType, C_correlType, "CorrelMatrix",	" ARM_ERR: correl type: string expected",		C_result);

		long objId;

		CCString stringId = GetLastCurCellEnvValue ();
		
		CCString curClass = LOCAL_GC_TARNFX_CLASS;

		if (!stringId)
		   objId = -1;
		else
		{
		   CCString prevClass = LocalGetStringObjectClass (stringId);
			
		   objId = LocalGetNumObjectId (stringId);
				
		   if ( curClass != prevClass)
           {
			  FreeCurCellContent();

			  objId = -1;
           }
		}
 
		retCode = ARMLOCAL_INITTARNFX ( LocalGetNumObjectId(C_tarnId),
										C_MCParams,
										(int)C_nbFactor,
										atoi(C_PDEParams[0].c_str()),
										atoi(C_PDEParams[1].c_str()),
										atof(C_PDEParams[2].c_str()),
										(C_PDEParams[3] == "Y"),
										(C_rescalling == "Y"),
										string(C_modelType),
										(C_smileFlag == "Y"),
										(C_mixCalib == "Y"),
										(C_oneFactorFlag == "Y"),
										string(C_correlType),
										C_zcCurvesId,
										C_basisCurvesId,
										C_forexId,
										C_ATMVolId, //for swopt BSGen
										C_bsFxVolId, //for BS fx models
										C_mixtureParamsId, //for mixture fx models
										C_mrsParamsId,
										C_QParamsId,
										LocalGetNumObjectId(C_correlMatrixId),
										C_result,
										objId);

		if ( retCode == ARM_OK )
		{
			objId = C_result.getLong ();

			LocalSetCurCellEnvValue (curClass, objId); 

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
	}
	ARM_XL_TRY_BLOCK_END

	ARM_XL_CATCH_ARM_EXPT

	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_ARM_INITTARNFX" )

	return (LPXLOPER)&XL_result;
}



__declspec(dllexport) LPXLOPER WINAPI Local_PXL_ARM_INITTARNFX( LPXLOPER XL_tarnId,
																LPXLOPER XL_zcCurvesId,
																LPXLOPER XL_basisCurvesId,
																LPXLOPER XL_forexId,
																LPXLOPER XL_ATMVolId,
																LPXLOPER XL_bsFxVolId,
																LPXLOPER XL_mixtureParamsId,
																LPXLOPER XL_mrsParamsId,
																LPXLOPER XL_QParamsId,
																LPXLOPER XL_correlMatrixId,
																LPXLOPER XL_MCParams,
																LPXLOPER XL_modelType,
																LPXLOPER XL_nbFactor,
																LPXLOPER XL_PDEParams,
																LPXLOPER XL_rescalling,
																LPXLOPER XL_smileFlag,
																LPXLOPER XL_mixCalib,
																LPXLOPER XL_oneFactorFlag,
																LPXLOPER XL_correlType )
{
	long retCode;
	ARM_result C_result;

	static XLOPER XL_result;

	ARM_XL_TRY_BLOCK_BEGIN
	{
		ARM_NOCALCIFWIZ();

		// C variable
		CCString C_tarnId;
		vector<CCString> C_zcCurvesId;
		vector<CCString> C_basisCurvesId;
		vector<CCString> C_forexId;
		vector<CCString> C_ATMVolId;
		vector<CCString> C_bsFxVolId;
		vector<CCString> C_bsFxVolId_default(0);
		vector<CCString> C_mixtureParamsId;
		vector<CCString> C_mixtureParamsId_default(0);
		vector<CCString> C_mrsParamsId;
		vector<CCString> C_QParamsId;
		vector<CCString> C_MCParams;
		CCString C_correlMatrixId;
		vector<CCString> C_PDEParams;
		vector<CCString> C_PDEParams_default(4);
		C_PDEParams_default[0] = "1000.0";
		C_PDEParams_default[1] = "901.0";
		C_PDEParams_default[2] = "6.0";
		C_PDEParams_default[3] = "N";
		double C_nbFactor;
		double C_nbFactor_default = 6.0;
		CCString C_rescalling, C_modelType, C_correlType, C_smileFlag, C_mixCalib, C_oneFactorFlag;

        // error
		static int error;
		static char* reason = "";

		XL_readStrCell(	  XL_tarnId,	 C_tarnId,					" ARM_ERR: Tarn id: object expected",		C_result);
		XL_readStrVector( XL_zcCurvesId, C_zcCurvesId,				" ARM_ERR: zc id: object expected",	DOUBLE_TYPE,	C_result);
		XL_readStrVector( XL_basisCurvesId, C_basisCurvesId,		" ARM_ERR: basis id: object expected", DOUBLE_TYPE,		C_result);
		XL_readStrVector( XL_forexId, C_forexId,					" ARM_ERR: fx id: object expected",	DOUBLE_TYPE,	C_result);
		XL_readStrVector( XL_ATMVolId, C_ATMVolId,					" ARM_ERR: ATM vol id: object expected", DOUBLE_TYPE,		C_result);
		XL_readStrVectorWD( XL_bsFxVolId, C_bsFxVolId, C_bsFxVolId_default,	" ARM_ERR: fx vol id: object expected",	DOUBLE_TYPE,	C_result);
		XL_readStrVectorWD( XL_mixtureParamsId, C_mixtureParamsId, C_mixtureParamsId_default, " ARM_ERR: mixture vol params: object expected", DOUBLE_TYPE,		C_result);

		if ((C_bsFxVolId.size() == 0) && (C_mixtureParamsId.size() == 0))
			throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"must provide BS fxVol or Mixture");

		XL_readStrVector( XL_mrsParamsId, C_mrsParamsId,			" ARM_ERR: mrs params id: object expected",	DOUBLE_TYPE,	C_result);
		XL_readStrVector( XL_QParamsId, C_QParamsId,				" ARM_ERR: Q params id: object expected", DOUBLE_TYPE,		C_result);
		XL_readStrVector( XL_MCParams, C_MCParams,					" ARM_ERR: MC params: vector of string expected", DOUBLE_TYPE,		C_result);
		XL_readStrCell(	  XL_correlMatrixId, C_correlMatrixId,		" ARM_ERR: correl matrix id: object expected",		C_result);
		XL_readNumCellWD( XL_nbFactor, C_nbFactor, C_nbFactor_default, " ARM_ERR: nb factors: integer expected",		C_result);
		XL_readStrVectorWD( XL_PDEParams, C_PDEParams, C_PDEParams_default, " ARM_ERR: PDE parameters: vector expected", DOUBLE_TYPE,		C_result);
		XL_readStrCellWD( XL_rescalling, C_rescalling, "N",		" ARM_ERR: rescalling: Y/N expected",		C_result);
		XL_readStrCell(	  XL_modelType, C_modelType,				" ARM_ERR: model type: string expected",		C_result);
		XL_readStrCellWD( XL_smileFlag, C_smileFlag, "Y",			" ARM_ERR: smile: Y/N expected",		C_result);
		XL_readStrCellWD( XL_mixCalib, C_mixCalib, "N",			" ARM_ERR: mix calib: Y/N expected",		C_result);
		XL_readStrCellWD( XL_oneFactorFlag, C_oneFactorFlag, "N",	" ARM_ERR: one factor flag: Y/N expected",		C_result);
		XL_readStrCellWD( XL_correlType, C_correlType, "CorrelMatrix",	" ARM_ERR: correl type: string expected",		C_result);

		long objId;
		CCString stringId;
		CCString curClass = LOCAL_GC_TARNFX_CLASS;
		 
		retCode = ARMLOCAL_INITTARNFX ( LocalGetNumObjectId(C_tarnId),
										C_MCParams,
										(int)C_nbFactor,
										atoi(C_PDEParams[0].c_str()),
										atoi(C_PDEParams[1].c_str()),
										atof(C_PDEParams[2].c_str()),
										(C_PDEParams[3] == "Y"),
										(C_rescalling == "Y"),
										string(C_modelType),
										(C_smileFlag == "Y"),
										(C_mixCalib == "Y"),
										(C_oneFactorFlag == "Y"),
										string(C_correlType),
										C_zcCurvesId,
										C_basisCurvesId,
										C_forexId,
										C_ATMVolId, //for swopt BSGen
										C_bsFxVolId, //for BS fx models
										C_mixtureParamsId, //for mixture fx models
										C_mrsParamsId,
										C_QParamsId,
										LocalGetNumObjectId(C_correlMatrixId),
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
	}
	ARM_XL_TRY_BLOCK_END

	ARM_XL_CATCH_ARM_EXPT

	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_PXL_ARM_INITTARNFX" )

	return (LPXLOPER)&XL_result;
}


/*------------------------------------------------------------------------------------------*/
/*---- End Of File ----*/
