#pragma warning(disable : 4786)
#pragma warning(disable : 4005)

#include <libCCxll\CCxll.h>

#include <ARM\libarm_local\ARM_local_mod.h>
#include <ARM\libarm_local\ARM_local_class.h>
#include <ARM\libarm_local\ARM_local_glob.h>
#include "ARM_local_interglob.h"
#include "ARM_local_interface.h"
#include "ARM_xl_wrapper_local.h"
#include <ARM\libarm_local\ARM_local_gp_calculators.h>

#include <mod3f\hw_vfdk_LDHD_lattice.h>
#include <GP_Base\gpbase\gpmatrix.h>


#include <util\fromto.h>

/// general macro for try catch
#include "ARM_xl_trycatch_local.h"

#include "util\tech_macro.h"



__declspec(dllexport) LPXLOPER WINAPI Local_YCMOD(LPXLOPER XL_zeroCurve,
												  LPXLOPER XL_discCurve)
{
	ADD_LOG("Local_YCMOD");
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

	    CCString C_zeroCurve;
	    CCString C_discCurve;
	    long C_discCurveId ;
	    
	    XL_readStrCell(XL_zeroCurve,C_zeroCurve," ARM_ERR: forecast zero-coupon curve id: object expected",C_result);
	    XL_readStrCellWD(XL_discCurve,C_discCurve,"DEFAULT"," ARM_ERR: discount zero-coupon curve id: object expected",C_result);

	    if (C_discCurve == "DEFAULT")
		    C_discCurveId = ARM_NULL_OBJECT;
	    else
		    C_discCurveId = LocalGetNumObjectId(C_discCurve);

	    long retCode;
	    long objId;
	    CCString prevClass;
	    
	    CCString curClass = LOCAL_YIELD_CURVE_BASIC_CLASS;
	    CCString stringId = GetLastCurCellEnvValue ();
	    
	    if (!stringId)
	    {
		    retCode = ARMLOCAL_ycmod( LocalGetNumObjectId(C_zeroCurve),
								      C_discCurveId,
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
			    
		    if (curClass == prevClass)
		    {
			    retCode = ARMLOCAL_ycmod( LocalGetNumObjectId(C_zeroCurve),
									      C_discCurveId,
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

			    retCode = ARMLOCAL_ycmod( LocalGetNumObjectId(C_zeroCurve),
									      C_discCurveId,
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_YCMOD" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


__declspec(dllexport) LPXLOPER WINAPI Local_PXL_YCMOD(LPXLOPER XL_zeroCurve,
													  LPXLOPER XL_discCurve
)
{
	ADD_LOG("Local_PXL_YCMOD");
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

	    CCString C_zeroCurve;
	    CCString C_discCurve;
	    long C_discCurveId ;
	    
	    XL_readStrCell(XL_zeroCurve,C_zeroCurve," ARM_ERR: forecast zero-coupon curve id: object expected",C_result);
	    XL_readStrCellWD(XL_discCurve,C_discCurve,"DEFAULT"," ARM_ERR: discount zero-coupon curve id: object expected",C_result);

	    if (C_discCurve == "DEFAULT")
		    C_discCurveId = ARM_NULL_OBJECT;
	    else
		    C_discCurveId = LocalGetNumObjectId(C_discCurve);

	    long retCode;
	    long objId;
	    
	    CCString curClass = LOCAL_YIELD_CURVE_BASIC_CLASS;
	    CCString stringId;

	    retCode = ARMLOCAL_ycmod( LocalGetNumObjectId(C_zeroCurve),
							      C_discCurveId,
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_PXL_YCMOD" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


__declspec(dllexport) LPXLOPER WINAPI Local_BSMODEL(LPXLOPER XL_date,
													LPXLOPER XL_spot,
													LPXLOPER XL_dividend,
													LPXLOPER XL_discrate,
													LPXLOPER XL_volat,
													LPXLOPER XL_typstk)
{
	ADD_LOG("Local_BSMODEL");

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
        double C_spot;

        double C_dividend_double;
        CCString C_dividend_str;
        long dividend_type;

        double C_discrate_double;
        CCString C_discrate_str;
        long discrate_type;

        double C_volat_double;
        CCString C_volat_str;
        long volat_type;

        CCString C_typstk;
        long typstk = K_YIELD;

        // error
        static int error;
        static char* reason = "";

        XL_readNumCell(XL_date,C_date," ARM_ERR: date: date expected",C_result);
        XL_readNumCell(XL_spot,C_spot," ARM_ERR: spot: numeric expected",C_result);

        XL_readStrOrNumCell(XL_dividend,C_dividend_str,C_dividend_double,dividend_type," ARM_ERR: dividend: string or numeric expected",C_result);
        XL_readStrOrNumCell(XL_discrate,C_discrate_str,C_discrate_double,discrate_type," ARM_ERR: discrate: string or numeric expected",C_result);
        XL_readStrOrNumCell(XL_volat,C_volat_str,C_volat_double,volat_type," ARM_ERR: volatility: string or numeric expected",C_result);

        XL_readStrCellWD(XL_typstk,C_typstk,"YIELD"," ARM_ERR: strike type: string expected",C_result);

        if ((typstk = ARM_ConvPriceYield (C_typstk, C_result)) == ARM_DEFAULT_ERR)
        {
           ARM_ARG_ERR();
           return (LPXLOPER)&XL_result;
        }

        if ( dividend_type == XL_TYPE_STRING )
        {
           C_dividend_double = (double) LocalGetNumObjectId(C_dividend_str);

           dividend_type = 1;
        }
        else
        {
           dividend_type = 0;
        }

        if ( discrate_type == XL_TYPE_STRING )
        {
           C_discrate_double = (double) LocalGetNumObjectId (C_discrate_str);
           discrate_type = 1;
        }
        else
        {
           discrate_type = 0;
        }

        if ( volat_type == XL_TYPE_STRING )
        {
           C_volat_double = (double) LocalGetNumObjectId (C_volat_str);

           volat_type = 1;
        }
        else
        {
           volat_type = 0;
        }

        long retCode;
        long objId;
        CCString prevClass;

        CCString curClass = LOCAL_BSMODEL_CLASS;
        CCString stringId = GetLastCurCellEnvValue ();

        if (!stringId)
        {
	        retCode = ARMLOCAL_bsmodel (C_date,C_spot,(long)dividend_type,C_dividend_double,
							        (long)discrate_type,C_discrate_double,(long)volat_type,C_volat_double,
							        typstk, C_result);

	        if(retCode == ARM_OK)
	        {
		        objId = C_result.getLong ();

		        LocalSetCurCellEnvValue (curClass, objId); 

		        stringId = LocalMakeObjectId (objId, curClass);
	        }
	        else
		        retCode = ARM_KO;
        }
        else
        {
	        prevClass = LocalGetStringObjectClass (stringId);
	        
	        objId = LocalGetNumObjectId (stringId);
		        
	        if (curClass == prevClass)
	        {
		        retCode = ARMLOCAL_bsmodel (C_date,C_spot,(long)dividend_type,C_dividend_double,
							        (long)discrate_type,C_discrate_double,(long)volat_type,C_volat_double,
							        typstk, C_result,objId);


		        if(retCode == ARM_OK)
		        {			
			        LocalSetCurCellEnvValue (curClass, objId); 

			        stringId = LocalMakeObjectId (objId, curClass);
		        }
	        }
	        else
	        {
		        FreeCurCellContent ();

		        retCode = ARMLOCAL_bsmodel (C_date,C_spot,(long)dividend_type,C_dividend_double,
							        (long)discrate_type,C_discrate_double,(long)volat_type,C_volat_double,
							        typstk, C_result);

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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_BSModel" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}
/*
__declspec(dllexport) LPXLOPER WINAPI Local_BSMODEL(LPXLOPER XL_date,
													LPXLOPER XL_spot,
													LPXLOPER XL_dividend,
													LPXLOPER XL_discrate,
													LPXLOPER XL_volat,
													LPXLOPER XL_typstk)
{
ADD_LOG("Local_BSMODEL");

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

		double C_spot_double;
        CCString C_spot_str;
        long spot_type;

        double C_dividend_double;
        CCString C_dividend_str;
        long dividend_type;

        double C_discrate_double;
        CCString C_discrate_str;
        long discrate_type;

        double C_volat_double;
        CCString C_volat_str;
        long volat_type;

        CCString C_typstk;
        long typstk = K_YIELD;

        // error
        static int error;
        static char* reason = "";

        XL_readNumCell(XL_date,C_date," ARM_ERR: date: date expected",C_result);
        XL_readNumCell(XL_spot,C_spot," ARM_ERR: spot: numeric expected",C_result);

        XL_readStrOrNumCell(XL_spot,	C_spot,			C_spot_double,		spot_type," ARM_ERR: spot: string or numeric expected",C_result);
        XL_readStrOrNumCell(XL_dividend,C_dividend_str,C_dividend_double,dividend_type," ARM_ERR: dividend: string or numeric expected",C_result);
        XL_readStrOrNumCell(XL_discrate,C_discrate_str,C_discrate_double,discrate_type," ARM_ERR: discrate: string or numeric expected",C_result);
        XL_readStrOrNumCell(XL_volat,C_volat_str,C_volat_double,volat_type," ARM_ERR: volatility: string or numeric expected",C_result);

        XL_readStrCellWD(XL_typstk,C_typstk,"YIELD"," ARM_ERR: strike type: string expected",C_result);

        if ((typstk = ARM_ConvPriceYield (C_typstk, C_result)) == ARM_DEFAULT_ERR)
        {
           ARM_ARG_ERR();
           return (LPXLOPER)&XL_result;
        }

        if ( dividend_type == XL_TYPE_STRING )
        {
           C_dividend_double = (double) LocalGetNumObjectId(C_dividend_str);

           dividend_type = 1;
        }
        else
        {
           dividend_type = 0;
        }

        if ( discrate_type == XL_TYPE_STRING )
        {
           C_discrate_double = (double) LocalGetNumObjectId (C_discrate_str);
           discrate_type = 1;
        }
        else
        {
           discrate_type = 0;
        }

        if ( volat_type == XL_TYPE_STRING )
        {
           C_volat_double = (double) LocalGetNumObjectId (C_volat_str);

           volat_type = 1;
        }
        else
        {
           volat_type = 0;
        }

        long retCode;
        long objId;
        CCString prevClass;

        CCString curClass = LOCAL_BSMODEL_CLASS;
        CCString stringId = GetLastCurCellEnvValue ();

        if (!stringId)
        {
	        retCode = ARMLOCAL_bsmodel (	C_date,
											(long)C_spot_type,	C_spot_double,
											(long)dividend_type,C_dividend_double,
											(long)discrate_type,C_discrate_double,
											(long)volat_type,	C_volat_double,
											typstk, 
											C_result);

	        if(retCode == ARM_OK)
	        {
		        objId = C_result.getLong ();

		        LocalSetCurCellEnvValue (curClass, objId); 

		        stringId = LocalMakeObjectId (objId, curClass);
	        }
	        else
		        retCode = ARM_KO;
        }
        else
        {
	        prevClass = LocalGetStringObjectClass (stringId);
	        
	        objId = LocalGetNumObjectId (stringId);
		        
	        if (curClass == prevClass)
	        {
		        retCode = ARMLOCAL_bsmodel (	C_date,
												(long)C_spot_type,	C_spot_double,
												(long)dividend_type,C_dividend_double,
												(long)discrate_type,C_discrate_double,
												(long)volat_type,	C_volat_double,
												typstk, 
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

		        retCode = ARMLOCAL_bsmodel (	C_date,
												(long)C_spot_type,	C_spot_double,
												(long)dividend_type,C_dividend_double,
												(long)discrate_type,C_discrate_double,
												(long)volat_type,	C_volat_double,
												typstk, 
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_BSModel" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}
*/

__declspec(dllexport) LPXLOPER WINAPI Local_PXL_BSMODEL(LPXLOPER XL_date,
														LPXLOPER XL_spot,
														LPXLOPER XL_dividend,
														LPXLOPER XL_discrate,
														LPXLOPER XL_volat,
														LPXLOPER XL_typstk)
{
	ADD_LOG("Local_PXL_BSMODEL");
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
	    double C_spot;
	    
	    double C_dividend_double;
	    CCString C_dividend_str;
	    long dividend_type;
	    
	    double C_discrate_double;
	    CCString C_discrate_str;
	    long discrate_type;
	    
	    double C_volat_double;
	    CCString C_volat_str;
	    long volat_type;
	    
	    CCString C_typstk;
	    long typstk = K_YIELD;
	    
	    // error
	    static int error;
	    static char* reason = "";

	    XL_readNumCell(XL_date,C_date," ARM_ERR: date: date expected",C_result);
	    XL_readNumCell(XL_spot,C_spot," ARM_ERR: spot: numeric expected",C_result);

	    XL_readStrOrNumCell(XL_dividend,C_dividend_str,C_dividend_double,dividend_type," ARM_ERR: dividend: string or numeric expected",C_result);
	    XL_readStrOrNumCell(XL_discrate,C_discrate_str,C_discrate_double,discrate_type," ARM_ERR: discrate: string or numeric expected",C_result);
	    XL_readStrOrNumCell(XL_volat,C_volat_str,C_volat_double,volat_type," ARM_ERR: volatility: string or numeric expected",C_result);

	    XL_readStrCellWD(XL_typstk,C_typstk,"YIELD"," ARM_ERR: strike type: string expected",C_result);

	    if ((typstk = ARM_ConvPriceYield (C_typstk, C_result)) == ARM_DEFAULT_ERR)
	    {
	       ARM_ARG_ERR();

	       return (LPXLOPER)&XL_result;
	    }
	    
	    if ( dividend_type == XL_TYPE_STRING )
	    {
	       C_dividend_double = (double) LocalGetNumObjectId(C_dividend_str);

	       dividend_type = 1;
	    }
	    else
	    {
	       dividend_type = 0;
	    }

	    if ( discrate_type == XL_TYPE_STRING )
	    {
	       C_discrate_double = (double)LocalGetNumObjectId (C_discrate_str);
	       discrate_type = 1;
	    }
	    else
	    {
	       discrate_type = 0;
	    }

	    if ( volat_type == XL_TYPE_STRING )
	    {
	       C_volat_double = (double) LocalGetNumObjectId (C_volat_str);

	       volat_type = 1;
	    }
	    else
	    {
	       volat_type = 0;
	    }

	    long retCode;
	    long objId;
	    
	    CCString curClass = LOCAL_BSMODEL_CLASS;
	    CCString stringId;

	    retCode = ARMLOCAL_bsmodel(C_date, C_spot, (long)dividend_type, C_dividend_double, 
						     (long) discrate_type, C_discrate_double, 
						     (long) volat_type, C_volat_double,
						     typstk, C_result);

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
	       PXL_ARM_ERR();
	    }

    //	ARM_END();	
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_PXL_BSMODEL" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}

/////////////////////////////////////////
/// Functor to create a seasonality manager
/////////////////////////////////////////
class NewBSFunc : public ARMResultLong2LongFunc
{
public:
	NewBSFunc(long YieldCurveId,
					  long volatilityCurveId,
					  long ConvAdjuVolatilityCurveId,
					  long CorrelManagerId,
					  long CnvxManagerId,
					  long SpreadLockCurveId,
                      long DiscountCurveId)
					 : 
					C_YieldCurveId(YieldCurveId),
					C_volatilityCurveId(volatilityCurveId),
					C_ConvAdjVolatilityCurveId(ConvAdjuVolatilityCurveId),
					C_CorrelManagerId(CorrelManagerId),
					C_CnvxManagerId(CnvxManagerId),
					C_SpreadLockCurveId(SpreadLockCurveId),
                    C_DiscountCurveId(DiscountCurveId)
	{};

	long operator()( ARM_result& result, long objId )
	{
		return ARMLOCAL_NewBSModel(
			C_YieldCurveId,
			C_volatilityCurveId,
			C_ConvAdjVolatilityCurveId,
			C_CorrelManagerId,
			C_CnvxManagerId,
			C_SpreadLockCurveId,
            C_DiscountCurveId,
			result, 
			objId);
	};

private:
	long C_YieldCurveId;
    long C_DiscountCurveId;
	long C_volatilityCurveId;
	long C_ConvAdjVolatilityCurveId;
	long C_CorrelManagerId;
	long C_CnvxManagerId;
	long C_SpreadLockCurveId;
};


__declspec(dllexport) LPXLOPER WINAPI Local_ARM_NEWBSMODELCommon(LPXLOPER XL_yieldcurve,
													        LPXLOPER XL_volatility,
															LPXLOPER XL_capletVol,
                                                            LPXLOPER XL_correlmanger,
															LPXLOPER XL_cnvxManager,
															LPXLOPER XL_spreadLock,
                                                            LPXLOPER XL_discountcurve,
															bool PersistentInXL)
{
	ADD_LOG("Local_ARM_NEWBSMODELCommon");
	static XLOPER XL_result;
	/// Get the variables from the XLOper variables
	ARM_result C_result;
	
	
	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{/// to avoid computation if called by the wizard
	    ARM_NOCALCIFWIZ();

	    // error
	    static int error;
	    static char* reason = "";

	    // C variable		
	    CCString C_yieldcurve;
	    long C_yieldcurveId;		
	    XL_readStrCell(XL_yieldcurve,C_yieldcurve," ARM_ERR: yieldcurve: string expected",C_result);
        C_yieldcurveId = LocalGetNumObjectId(C_yieldcurve);

        CCString C_volatility_str;
        long C_volatilityId;
	    XL_readStrCell(XL_volatility,C_volatility_str," ARM_ERR: volatility id : object expected",C_result);
        C_volatilityId = LocalGetNumObjectId (C_volatility_str);

        CCString C_capletVol_str;
        long C_capletVolId;
	    XL_readStrCellWD(XL_capletVol,C_capletVol_str,"NULL"," ARM_ERR: caplet volatility id : object expected",C_result);
        C_capletVolId = (C_capletVol_str == "NULL")? ARM_NULL_OBJECT : LocalGetNumObjectId (C_capletVol_str);

        CCString C_correlmanger;
        long C_correlmangerId;
        XL_readStrCellWD(XL_correlmanger, C_correlmanger,"NULL"," ARM_ERR: CorrelManager id: object expected",C_result);
        C_correlmangerId = (C_correlmanger == "NULL")? ARM_NULL_OBJECT : LocalGetNumObjectId (C_correlmanger);

        CCString C_cnvxManager;
	    long C_cnvxManagerId;
	    XL_readStrCellWD(XL_cnvxManager, C_cnvxManager, "NULL", "Convexity Manager Id: string expected", C_result);
        C_cnvxManagerId = (C_cnvxManager == "NULL")? ARM_NULL_OBJECT : LocalGetNumObjectId (C_cnvxManager);

	    CCString C_spreadLock_str;
	    long C_spreadLockId;
	    XL_readStrCellWD(XL_spreadLock, C_spreadLock_str, "NULL", " ARM_ERR: spread lock id : string expected", C_result);
        C_spreadLockId = (C_spreadLock_str == "NULL")? ARM_NULL_OBJECT : LocalGetNumObjectId (C_spreadLock_str);

        CCString C_discountcurve_str;
	    long C_discountcurveId;		
	    XL_readStrCellWD(XL_discountcurve,C_discountcurve_str, "NULL", " ARM_ERR: discountcurve: string expected",C_result);
        C_discountcurveId = (C_discountcurve_str == "NULL")? ARM_NULL_OBJECT : LocalGetNumObjectId (C_discountcurve_str);
	    

	    NewBSFunc ourFunc(
		    C_yieldcurveId,
		    C_volatilityId,
		    C_capletVolId,
		    C_correlmangerId,
		    C_cnvxManagerId,
		    C_spreadLockId,
            C_discountcurveId
		    );

	    /// call the general function
	    fillXL_Result( LOCAL_BSMODEL_CLASS, ourFunc, C_result, XL_result, PersistentInXL );

	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_ARM_NEWBSMODELCommon" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}



__declspec(dllexport) LPXLOPER WINAPI Local_ARM_NEWBSMODEL(LPXLOPER XL_yieldcurve,													
													        LPXLOPER XL_volatility,
                                                            LPXLOPER XL_correlmanager,
															LPXLOPER XL_cnvxManager,
                                                            LPXLOPER XL_convadjVolat,
															LPXLOPER XL_spreadLock,
                                                            LPXLOPER XL_discountcurve)
{
	ADD_LOG("Local_ARM_NEWBSMODEL");
	bool PersistentInXL = true;
	return Local_ARM_NEWBSMODELCommon(
		XL_yieldcurve,
		XL_volatility,
		XL_convadjVolat,
		XL_correlmanager,
		XL_cnvxManager,
		XL_spreadLock,
        XL_discountcurve,
		PersistentInXL);

}


__declspec(dllexport) LPXLOPER WINAPI Local_PXL_NEWBSMODEL(LPXLOPER XL_yieldcurve,													
													        LPXLOPER XL_volatility,
                                                            LPXLOPER XL_correlmanager,
															LPXLOPER XL_cnvxManager,
                                                            LPXLOPER  XL_convadjVolat,
															LPXLOPER XL_spreadLock,
                                                            LPXLOPER XL_discountcurve)
{
	ADD_LOG("Local_PXL_NEWBSMODEL");
	bool PersistentInXL = false;
	return Local_ARM_NEWBSMODELCommon(
		XL_yieldcurve,
		XL_volatility,
		XL_convadjVolat,
		XL_correlmanager,
		XL_cnvxManager,
		XL_spreadLock,
        XL_discountcurve,
		PersistentInXL);
}


__declspec(dllexport) LPXLOPER WINAPI Local_ARM_BSMODELGEN(	LPXLOPER XL_yieldcurve,
															LPXLOPER XL_volatility,
															LPXLOPER XL_correlmanger,
															LPXLOPER XL_cnvxManager,
															LPXLOPER XL_capletVol,//ConvAdjVolatility
															LPXLOPER XL_spreadLock,
															LPXLOPER XL_discountcurve,
															LPXLOPER XL_correl,
															LPXLOPER XL_cashVol,
															LPXLOPER XL_spreadVol,
															LPXLOPER XL_modelType,
															LPXLOPER XL_spreadVolType,
															LPXLOPER XL_sabrMod,
                                                            LPXLOPER XL_LnOrNorVol,
															LPXLOPER XL_numSteps)
{
	ADD_LOG("Local_ARM_BSMODELGEN");
	static XLOPER XL_result;
	/// Get the variables from the XLOper variables
	ARM_result C_result;
	
	
	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{/// to avoid computation if called by the wizard
	    ARM_NOCALCIFWIZ();

	    // error
	    static int error;
	    static char* reason = "";

	    // C variable		
	    CCString C_yieldcurve;
	    long C_yieldcurveId;		
	    XL_readStrCell(XL_yieldcurve,C_yieldcurve," ARM_ERR: yieldcurve: string expected",C_result);
		C_yieldcurveId = LocalGetNumObjectId(C_yieldcurve);

        CCString C_volatility_str;
        long C_volatilityId;
	    XL_readStrCell(XL_volatility,C_volatility_str," ARM_ERR: volatility id : object expected",C_result);
        C_volatilityId = LocalGetNumObjectId (C_volatility_str);

        CCString C_correlmanger;
        long C_correlmangerId;
        XL_readStrCellWD(XL_correlmanger, C_correlmanger,"NULL"," ARM_ERR: CorrelManager id: object expected",C_result);
        C_correlmangerId = (C_correlmanger == "NULL")? ARM_NULL_OBJECT : LocalGetNumObjectId (C_correlmanger);

        CCString C_cnvxManager;
	    long C_cnvxManagerId;
	    XL_readStrCellWD(XL_cnvxManager, C_cnvxManager, "NULL", "Convexity Manager Id: string expected", C_result);
        C_cnvxManagerId = (C_cnvxManager == "NULL")? ARM_NULL_OBJECT : LocalGetNumObjectId (C_cnvxManager);

        CCString C_capletVol_str;
        long C_capletVolId;
	    XL_readStrCellWD(XL_capletVol,C_capletVol_str,"NULL"," ARM_ERR: caplet volatility id : object expected",C_result);
        C_capletVolId = (C_capletVol_str == "NULL")? ARM_NULL_OBJECT : LocalGetNumObjectId (C_capletVol_str);

	    CCString C_spreadLock_str;
	    long C_spreadLockId;
	    XL_readStrCellWD(XL_spreadLock, C_spreadLock_str, "NULL", " ARM_ERR: spread lock id : string expected", C_result);
        C_spreadLockId = (C_spreadLock_str == "NULL")? ARM_NULL_OBJECT : LocalGetNumObjectId (C_spreadLock_str);

        CCString C_discountcurve_str;
	    long C_discountcurveId;		
	    XL_readStrCellWD(XL_discountcurve,C_discountcurve_str, "NULL", " ARM_ERR: discountcurve: string expected",C_result);
        C_discountcurveId = (C_discountcurve_str == "NULL")? ARM_NULL_OBJECT : LocalGetNumObjectId (C_discountcurve_str);
	    
		CCString C_correl_str;
		long C_correlId;
		XL_readStrCellWD(XL_correl,C_correl_str,"NULL"," ARM_ERR: correlation id : object expected",C_result);
		C_correlId = (C_correl_str == "NULL")? ARM_NULL_OBJECT : LocalGetNumObjectId (C_correl_str);

		CCString C_cashVol_str;
		long C_cashVolId;
		XL_readStrCellWD(XL_cashVol,C_cashVol_str,"NULL"," ARM_ERR: cash volatility id : object expected",C_result);
		C_cashVolId = (C_cashVol_str == "NULL")? ARM_NULL_OBJECT : LocalGetNumObjectId (C_cashVol_str);

		CCString C_spreadVol_str;
		long C_spreadVolId;
		XL_readStrCellWD(XL_spreadVol,C_spreadVol_str,"NULL"," ARM_ERR: spread volatility id : object expected",C_result);
		C_spreadVolId = (C_spreadVol_str == "NULL")? ARM_NULL_OBJECT : LocalGetNumObjectId (C_spreadVol_str);

		CCString C_modelType;
		long modelTypeId;
		XL_readStrCellWD(XL_modelType, C_modelType,"2LOG"," ARM_ERR: model type: string expected",C_result);

		CCString C_spreadVolType;
		long spreadVolTypeId;
		XL_readStrCellWD(XL_spreadVolType, C_spreadVolType,"COMPUTE"," ARM_ERR: spread type: string expected",C_result);

		CCString C_sabrModId;
		long sabrModId;
		XL_readStrCellWD(XL_sabrMod, C_sabrModId,"DEFAULT"," ARM_ERR: sabr model id: string expected",C_result);
		
		if(C_sabrModId == "DEFAULT")
		{
			sabrModId = ARM_NULL_OBJECT;
		}
		else
		{
			sabrModId = LocalGetNumObjectId(C_sabrModId);
		}

	    CCString lnOrNorVol;
	    XL_readStrCellWD(XL_LnOrNorVol,lnOrNorVol,"Y"," ARM_ERR: LogNor or Nor Vol: array of string expected (Y/N)",C_result);
	    bool isLnVol=true;
        lnOrNorVol.toUpper();
        if(CCSTringToSTLString(lnOrNorVol)!="Y")
                isLnVol=false;

		if ((modelTypeId = ARM_ConvModelType (C_modelType, C_result)) == ARM_DEFAULT_ERR)
		{
		   ARM_ARG_ERR();
		   return (LPXLOPER)&XL_result;
		}

		if ((spreadVolTypeId = ARM_ConvVolType2 (C_spreadVolType, C_result)) == ARM_DEFAULT_ERR)
		{
		   ARM_ARG_ERR();
		   return (LPXLOPER)&XL_result;
		}

		double hundred = 100.;
		double C_numSteps;
		XL_readNumCellWD(XL_numSteps, C_numSteps, hundred, " ARM_ERR: Numerical Num Steps : numeric expected",C_result);
	
		long retCode;
		long objId;
		CCString prevClass;
		
		CCString curClass = LOCAL_BSMODELGEN_CLASS;
		CCString stringId = GetLastCurCellEnvValue ();

		if (!stringId)
		{
			retCode = ARMLOCAL_BSMODELGEN(C_yieldcurveId,
										   C_spreadLockId,
										   C_capletVolId,
										   C_volatilityId,
										   C_correlmangerId,
										   C_cnvxManagerId,
										   C_discountcurveId,
										   C_correlId,
										   C_cashVolId,
										   C_spreadVolId,
										   long(modelTypeId),
										   long(spreadVolTypeId),
										   long(sabrModId),
                                           isLnVol,
										   C_numSteps,
										   C_result);

			if ( retCode == ARM_OK )
			{
				objId = C_result.getLong();
				
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
				retCode = ARMLOCAL_BSMODELGEN(C_yieldcurveId,
											     C_spreadLockId,
											     C_capletVolId,
											     C_volatilityId,
											     C_correlmangerId,
											     C_cnvxManagerId,
											     C_discountcurveId,
											     C_correlId,
											     C_cashVolId,
												 C_spreadVolId,
											     long(modelTypeId),
												 long(spreadVolTypeId),
												 long(sabrModId),
                                                 isLnVol,
												 long(C_numSteps),
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

				retCode = ARMLOCAL_BSMODELGEN(C_yieldcurveId,
												 C_spreadLockId,
											     C_capletVolId,
											     C_volatilityId,
											     C_correlmangerId,
											     C_cnvxManagerId,
											     C_discountcurveId,
											     C_correlId,
											     C_cashVolId,
												 C_spreadVolId,
											     long(modelTypeId),
												 long(spreadVolTypeId),
												 long(sabrModId),
                                                 isLnVol,
												 long(C_numSteps),
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

	//ARM_END();

	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_ARM_NEWBSMODELGELCommon" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}

__declspec(dllexport) LPXLOPER WINAPI Local_PXL_ARM_BSMODELGEN(	LPXLOPER XL_yieldcurve,
																LPXLOPER XL_volatility,
																LPXLOPER XL_correlmanger,
																LPXLOPER XL_cnvxManager,
																LPXLOPER XL_capletVol,//ConvAdjVolatility
																LPXLOPER XL_spreadLock,
																LPXLOPER XL_discountcurve,
																LPXLOPER XL_correl,
																LPXLOPER XL_cashVol,
																LPXLOPER XL_spreadVol,
																LPXLOPER XL_modelType,
																LPXLOPER XL_spreadVolType,
																LPXLOPER XL_sabrMod,
                                                                LPXLOPER XL_LnOrNorVol,
																LPXLOPER XL_numSteps)
{
	ADD_LOG("Local_PXL_ARM_BSMODELGEN");
	static XLOPER XL_result;
	/// Get the variables from the XLOper variables
	ARM_result C_result;
	
	
	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{/// to avoid computation if called by the wizard
	    ARM_NOCALCIFWIZ();

	    // error
	    static int error;
	    static char* reason = "";

	    // C variable		
	    CCString C_yieldcurve;
	    long C_yieldcurveId;		
	    XL_readStrCell(XL_yieldcurve,C_yieldcurve," ARM_ERR: yieldcurve: string expected",C_result);
		C_yieldcurveId = LocalGetNumObjectId(C_yieldcurve);

        CCString C_volatility_str;
        long C_volatilityId;
	    XL_readStrCell(XL_volatility,C_volatility_str," ARM_ERR: volatility id : object expected",C_result);
        C_volatilityId = LocalGetNumObjectId (C_volatility_str);

        CCString C_correlmanger;
        long C_correlmangerId;
        XL_readStrCellWD(XL_correlmanger, C_correlmanger,"NULL"," ARM_ERR: CorrelManager id: object expected",C_result);
        C_correlmangerId = (C_correlmanger == "NULL")? ARM_NULL_OBJECT : LocalGetNumObjectId (C_correlmanger);

        CCString C_cnvxManager;
	    long C_cnvxManagerId;
	    XL_readStrCellWD(XL_cnvxManager, C_cnvxManager, "NULL", "Convexity Manager Id: string expected", C_result);
        C_cnvxManagerId = (C_cnvxManager == "NULL")? ARM_NULL_OBJECT : LocalGetNumObjectId (C_cnvxManager);

        CCString C_capletVol_str;
        long C_capletVolId;
	    XL_readStrCellWD(XL_capletVol,C_capletVol_str,"NULL"," ARM_ERR: caplet volatility id : object expected",C_result);
        C_capletVolId = (C_capletVol_str == "NULL")? ARM_NULL_OBJECT : LocalGetNumObjectId (C_capletVol_str);

	    CCString C_spreadLock_str;
	    long C_spreadLockId;
	    XL_readStrCellWD(XL_spreadLock, C_spreadLock_str, "NULL", " ARM_ERR: spread lock id : string expected", C_result);
        C_spreadLockId = (C_spreadLock_str == "NULL")? ARM_NULL_OBJECT : LocalGetNumObjectId (C_spreadLock_str);

        CCString C_discountcurve_str;
	    long C_discountcurveId;		
	    XL_readStrCellWD(XL_discountcurve,C_discountcurve_str, "NULL", " ARM_ERR: discountcurve: string expected",C_result);
        C_discountcurveId = (C_discountcurve_str == "NULL")? ARM_NULL_OBJECT : LocalGetNumObjectId (C_discountcurve_str);
	    
		CCString C_correl_str;
		long C_correlId;
		XL_readStrCellWD(XL_correl,C_correl_str,"NULL"," ARM_ERR: correlation id : object expected",C_result);
		C_correlId = (C_correl_str == "NULL")? ARM_NULL_OBJECT : LocalGetNumObjectId (C_correl_str);

		CCString C_cashVol_str;
		long C_cashVolId;
		XL_readStrCellWD(XL_cashVol,C_cashVol_str,"NULL"," ARM_ERR: cash volatility id : object expected",C_result);
		C_cashVolId = (C_cashVol_str == "NULL")? ARM_NULL_OBJECT : LocalGetNumObjectId (C_cashVol_str);

		CCString C_spreadVol_str;
		long C_spreadVolId;
		XL_readStrCellWD(XL_spreadVol,C_spreadVol_str,"NULL"," ARM_ERR: spread volatility id : object expected",C_result);
		C_spreadVolId = (C_spreadVol_str == "NULL")? ARM_NULL_OBJECT : LocalGetNumObjectId (C_spreadVol_str);

		CCString C_modelType;
		long modelTypeId;
		XL_readStrCellWD(XL_modelType, C_modelType,"2LOG"," ARM_ERR: model type: string expected",C_result);

		CCString C_spreadVolType;
		long spreadVolTypeId;
		XL_readStrCellWD(XL_spreadVolType, C_spreadVolType,"COMPUTE"," ARM_ERR: spread type: string expected",C_result);

		CCString C_sabrModId;
		long sabrModId;
		XL_readStrCellWD(XL_sabrMod, C_sabrModId,"DEFAULT"," ARM_ERR: sabr model id: string expected",C_result);
		
		if(C_sabrModId == "DEFAULT")
		{
			sabrModId = ARM_NULL_OBJECT;
		}
		else
		{
			sabrModId = LocalGetNumObjectId(C_sabrModId);
		}

	    CCString lnOrNorVol;
	    XL_readStrCellWD(XL_LnOrNorVol,lnOrNorVol,"Y"," ARM_ERR: LogNor or Nor Vol: array of string expected (Y/N)",C_result);
	    bool isLnVol=true;
        lnOrNorVol.toUpper();
        if(CCSTringToSTLString(lnOrNorVol)!="Y")
                isLnVol=false;

		if ((modelTypeId = ARM_ConvModelType (C_modelType, C_result)) == ARM_DEFAULT_ERR)
		{
		   ARM_ARG_ERR();
		   return (LPXLOPER)&XL_result;
		}

		if ((spreadVolTypeId = ARM_ConvVolType2 (C_spreadVolType, C_result)) == ARM_DEFAULT_ERR)
		{
		   ARM_ARG_ERR();
		   return (LPXLOPER)&XL_result;
		}

		double hundred = 100.;
		double C_numSteps;
		XL_readNumCellWD(XL_numSteps, C_numSteps, hundred, " ARM_ERR: Numerical Num Steps : numeric expected",C_result);
	
		long retCode;
		long objId;
		CCString curClass = LOCAL_BSMODELGEN_CLASS;
		CCString stringId ; 

		retCode = ARMLOCAL_BSMODELGEN(C_yieldcurveId,
										   C_spreadLockId,
										   C_capletVolId,
										   C_volatilityId,
										   C_correlmangerId,
										   C_cnvxManagerId,
										   C_discountcurveId,
										   C_correlId,
										   C_cashVolId,
										   C_spreadVolId,
										   long(modelTypeId),
										   long(spreadVolTypeId),
										   long(sabrModId),
                                           isLnVol,
										   long (C_numSteps),
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

	//ARM_END();

	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_ARM_NEWBSMODELGELCommon" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}



__declspec(dllexport) LPXLOPER WINAPI Local_ARM_BSPricingModel(LPXLOPER XL_yieldCurve,
															   LPXLOPER XL_convexityManager,
															   LPXLOPER XL_capModel,
															   LPXLOPER XL_swoptModel,
															   LPXLOPER XL_correlManager,
															   LPXLOPER XL_modelType,
															   LPXLOPER XL_spreadVolCurve,
															   LPXLOPER XL_discountCurve,
															   LPXLOPER XL_adjConvVolCurve,
															   LPXLOPER XL_calibInfos,
															   LPXLOPER XL_numInfos)
{
	ADD_LOG("Local_ARM_BSPricingModel");
	static XLOPER XL_result;
	ARM_result C_result;
	
	ARM_XL_TRY_BLOCK_BEGIN
	{
	    ARM_NOCALCIFWIZ();

	    // error
	    static int error;
	    static char* reason = "";

	    // C variable		
	    CCString C_yieldCurve;
	    long C_yieldCurveId;		
	    XL_readStrCell(XL_yieldCurve,C_yieldCurve," ARM_ERR: yield curve: object expected",C_result);
		C_yieldCurveId = LocalGetNumObjectId(C_yieldCurve);

        CCString C_convexityManager;
	    long C_convexityManagerId;
	    XL_readStrCellWD(XL_convexityManager, C_convexityManager, "NULL", "Convexity Manager: object expected", C_result);
        C_convexityManagerId = (C_convexityManager == "NULL")? ARM_NULL_OBJECT : LocalGetNumObjectId (C_convexityManager);

        CCString C_capModel;
	    long C_capModelId;
	    XL_readStrCellWD(XL_capModel, C_capModel, "NULL", "cap model: object expected", C_result);
        C_capModelId = (C_capModel == "NULL")? ARM_NULL_OBJECT : LocalGetNumObjectId (C_capModel);

        CCString C_swoptModel;
	    long C_swoptModelId;
	    XL_readStrCellWD(XL_swoptModel, C_swoptModel, "NULL", "Swopt model: object expected", C_result);
        C_swoptModelId = (C_swoptModel == "NULL")? ARM_NULL_OBJECT : LocalGetNumObjectId (C_swoptModel);

        CCString C_correlManager;
        long C_correlManagerId;
        XL_readStrCellWD(XL_correlManager, C_correlManager,"NULL"," ARM_ERR: CorrelManager: object expected",C_result);
        C_correlManagerId = (C_correlManager == "NULL")? ARM_NULL_OBJECT : LocalGetNumObjectId (C_correlManager);

		CCString C_modelType;
		long modelTypeId;
		XL_readStrCellWD(XL_modelType, C_modelType,"2LOG"," ARM_ERR: model type: string expected",C_result);

		if ((modelTypeId = ARM_ConvModelType (C_modelType, C_result)) == ARM_DEFAULT_ERR)
		{
		   ARM_ARG_ERR();
		   return (LPXLOPER)&XL_result;
		}

		CCString C_spreadVolCurve;
		long C_spreadVolCurveId;
		XL_readStrCellWD(XL_spreadVolCurve, C_spreadVolCurve,"NULL"," ARM_ERR: spread volatility: object expected",C_result);
        C_spreadVolCurveId = (C_spreadVolCurve == "NULL")? ARM_NULL_OBJECT : LocalGetNumObjectId (C_spreadVolCurve);

        CCString C_discountCurve;
	    long C_discountCurveId;		
	    XL_readStrCellWD(XL_discountCurve,C_discountCurve, "NULL", " ARM_ERR: discount curve: object expected",C_result);
        C_discountCurveId = (C_discountCurve == "NULL")? ARM_NULL_OBJECT : LocalGetNumObjectId (C_discountCurve);
	    
	    CCString C_adjConvVolCurve;
	    long C_adjConvVolCurveId;
	    XL_readStrCellWD(XL_adjConvVolCurve, C_adjConvVolCurve, "NULL", " ARM_ERR: adj vol Curve: object expected", C_result);
        C_adjConvVolCurveId = (C_adjConvVolCurve == "NULL")? ARM_NULL_OBJECT : LocalGetNumObjectId(C_adjConvVolCurve);

		VECTOR<double> C_calibInfos;
		VECTOR<double> C_calibInfosDef(0);
		XL_readNumVectorWD(XL_calibInfos, C_calibInfos, C_calibInfosDef, " ARM_ERR: Calib infos: array of number expected", C_result);		

		VECTOR<double> C_numInfos;
		VECTOR<double> C_numInfosDef(1, 100.0);
		XL_readNumVectorWD(XL_numInfos, C_numInfos, C_numInfosDef, " ARM_ERR: Num infos: array of number expected", C_result);		
	
		long retCode;
		long objId;
		CCString prevClass;
		
		CCString curClass = LOCAL_BSMODEL_CLASS;
		CCString stringId = GetLastCurCellEnvValue ();

		if (!stringId)
		{
			retCode = ARMLOCAL_BSPricingModel(C_yieldCurveId,
											  C_convexityManagerId,
											  C_capModelId,
											  C_swoptModelId,
											  C_correlManagerId,
											  modelTypeId,
											  C_spreadVolCurveId,
											  C_discountCurveId,
											  C_adjConvVolCurveId,
											  C_calibInfos,
											  C_numInfos,
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
			prevClass = LocalGetStringObjectClass(stringId);

			objId = LocalGetNumObjectId(stringId);
				
			if ( curClass == prevClass )
			{
				retCode = ARMLOCAL_BSPricingModel(C_yieldCurveId,
												  C_convexityManagerId,
												  C_capModelId,
												  C_swoptModelId,
												  C_correlManagerId,
												  modelTypeId,
												  C_spreadVolCurveId,
												  C_discountCurveId,
												  C_adjConvVolCurveId,
												  C_calibInfos,
												  C_numInfos,													C_result,
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

				retCode = ARMLOCAL_BSPricingModel(C_yieldCurveId,
												  C_convexityManagerId,
												  C_capModelId,
												  C_swoptModelId,
												  C_correlManagerId,
												  modelTypeId,
												  C_spreadVolCurveId,
												  C_discountCurveId,
												  C_adjConvVolCurveId,
												  C_calibInfos,
												  C_numInfos,
												  C_result );

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

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_ARM_BSPricingModel" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}



__declspec(dllexport) LPXLOPER WINAPI Local_PXL_ARM_BSPricingModel(LPXLOPER XL_yieldCurve,
																   LPXLOPER XL_convexityManager,
																   LPXLOPER XL_capModel,
																   LPXLOPER XL_swoptModel,
																   LPXLOPER XL_correlManager,
																   LPXLOPER XL_modelType,
																   LPXLOPER XL_spreadVolCurve,
																   LPXLOPER XL_discountCurve,
																   LPXLOPER XL_adjConvVolCurve,
																   LPXLOPER XL_calibInfos,
																   LPXLOPER XL_numInfos)
{
	ADD_LOG("Local_PXL_ARM_BSPricingModel");
	static XLOPER XL_result;
	ARM_result C_result;
	
	ARM_XL_TRY_BLOCK_BEGIN
	{
	    ARM_NOCALCIFWIZ();

	    // error
	    static int error;
	    static char* reason = "";

	    // C variable		
	    CCString C_yieldCurve;
	    long C_yieldCurveId;		
	    XL_readStrCell(XL_yieldCurve,C_yieldCurve," ARM_ERR: yield curve: object expected",C_result);
		C_yieldCurveId = LocalGetNumObjectId(C_yieldCurve);

        CCString C_convexityManager;
	    long C_convexityManagerId;
	    XL_readStrCellWD(XL_convexityManager, C_convexityManager, "NULL", "Convexity Manager: object expected", C_result);
        C_convexityManagerId = (C_convexityManager == "NULL")? ARM_NULL_OBJECT : LocalGetNumObjectId (C_convexityManager);

        CCString C_capModel;
	    long C_capModelId;
	    XL_readStrCellWD(XL_capModel, C_capModel, "NULL", "cap model: object expected", C_result);
        C_capModelId = (C_capModel == "NULL")? ARM_NULL_OBJECT : LocalGetNumObjectId (C_capModel);

        CCString C_swoptModel;
	    long C_swoptModelId;
	    XL_readStrCellWD(XL_swoptModel, C_swoptModel, "NULL", "Swopt model: object expected", C_result);
        C_swoptModelId = (C_swoptModel == "NULL")? ARM_NULL_OBJECT : LocalGetNumObjectId (C_swoptModel);

        CCString C_correlManager;
        long C_correlManagerId;
        XL_readStrCellWD(XL_correlManager, C_correlManager,"NULL"," ARM_ERR: CorrelManager: object expected",C_result);
        C_correlManagerId = (C_correlManager == "NULL")? ARM_NULL_OBJECT : LocalGetNumObjectId (C_correlManager);

		CCString C_modelType;
		long modelTypeId;
		XL_readStrCellWD(XL_modelType, C_modelType,"2LOG"," ARM_ERR: model type: string expected",C_result);

		if ((modelTypeId = ARM_ConvModelType (C_modelType, C_result)) == ARM_DEFAULT_ERR)
		{
		   ARM_ARG_ERR();
		   return (LPXLOPER)&XL_result;
		}

		CCString C_spreadVolCurve;
		long C_spreadVolCurveId;
		XL_readStrCellWD(XL_spreadVolCurve,C_spreadVolCurve,"NULL"," ARM_ERR: spread volatility: object expected",C_result);
        C_spreadVolCurveId = (C_spreadVolCurve == "NULL")? ARM_NULL_OBJECT : LocalGetNumObjectId (C_spreadVolCurve);


        CCString C_discountCurve;
	    long C_discountCurveId;		
	    XL_readStrCellWD(XL_discountCurve,C_discountCurve, "NULL", " ARM_ERR: discount curve: object expected",C_result);
        C_discountCurveId = (C_discountCurve == "NULL")? ARM_NULL_OBJECT : LocalGetNumObjectId (C_discountCurve);
	    
	    CCString C_adjConvVolCurve;
	    long C_adjConvVolCurveId;
	    XL_readStrCellWD(XL_adjConvVolCurve, C_adjConvVolCurve, "NULL", " ARM_ERR: adj vol Curve: object expected", C_result);
        C_adjConvVolCurveId = (C_adjConvVolCurve == "NULL")? ARM_NULL_OBJECT : LocalGetNumObjectId(C_adjConvVolCurve);

		VECTOR<double> C_calibInfos;
		VECTOR<double> C_calibInfosDef(0);
		XL_readNumVectorWD(XL_calibInfos, C_calibInfos, C_calibInfosDef, " ARM_ERR: Calib infos: array of number expected", C_result);		

		VECTOR<double> C_numInfos;
		VECTOR<double> C_numInfosDef(1, 100.0);
		XL_readNumVectorWD(XL_numInfos, C_numInfos, C_numInfosDef, " ARM_ERR: Num infos: array of number expected", C_result);		
	
		long retCode;
		long objId;
		CCString prevClass;
		
		CCString curClass = LOCAL_BSMODEL_CLASS;
		CCString stringId = GetLastCurCellEnvValue ();

		retCode = ARMLOCAL_BSPricingModel(C_yieldCurveId,
										  C_convexityManagerId,
										  C_capModelId,
										  C_swoptModelId,
										  C_correlManagerId,
										  modelTypeId,
										  C_spreadVolCurveId,
										  C_discountCurveId,
										  C_adjConvVolCurveId,
										  C_calibInfos,
										  C_numInfos,
										  C_result );

		if ( retCode == ARM_OK )
		{
			objId = C_result.getLong();
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

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_PXL_ARM_BSPricingModel" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


__declspec(dllexport) LPXLOPER WINAPI Local_GTWOYC (LPXLOPER XL_dMeanRevSpeed,
													LPXLOPER XL_dSigma,
													LPXLOPER XL_dZcId,
													LPXLOPER XL_fZcId,
													LPXLOPER XL_ratesCorr,
													LPXLOPER XL_fMeanRevSpeed,
													LPXLOPER XL_fSigma)
{
	ADD_LOG("Local_GTWOYC ");
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	    ARM_NOCALCIFWIZ();

	    // C variable
	    double C_dMeanRevSpeed;
	    double C_dSigma;

	    CCString C_dZcId;
	    CCString C_fZcId;

	    double C_ratesCorr;
	    double C_fMeanRevSpeed;
	    double C_fSigma;
	    
	    // error
	    static int error;
	    static char* reason = "";

	    XL_readNumCell(XL_dMeanRevSpeed,C_dMeanRevSpeed," ARM_ERR: domestic mean reversion speed: numeric expected",C_result);
	    XL_readNumCell(XL_dSigma,C_dSigma," ARM_ERR: domestic sigma: numeric expected",C_result);
	    XL_readStrCell(XL_dZcId,C_dZcId," ARM_ERR: domestic zero-coupon curve id: object expected",C_result);
	    XL_readStrCell(XL_fZcId,C_fZcId," ARM_ERR: foreign zero-coupon curve id: object expected",C_result);
	    XL_readNumCell(XL_ratesCorr,C_ratesCorr," ARM_ERR: rates correlation: numeric expected",C_result);
	    XL_readNumCell(XL_fMeanRevSpeed,C_fMeanRevSpeed," ARM_ERR: foreign mean reversion speed: numeric expected",C_result);
	    XL_readNumCell(XL_fSigma,C_fSigma," ARM_ERR: foreign sigma: numeric expected",C_result);
	    
	    long retCode;
	    long objId;
	    CCString prevClass;
	    
	    CCString curClass = LOCAL_GAUSSIAN_2YC_MODEL_CLASS;
	    CCString stringId = GetLastCurCellEnvValue ();

	    if (!stringId)
	    {
		    retCode = ARMLOCAL_GTWOYC (C_dMeanRevSpeed, C_dSigma,
									    LocalGetNumObjectId (C_dZcId),
									    LocalGetNumObjectId (C_fZcId),
									    C_ratesCorr,
									    C_fMeanRevSpeed,
									    C_fSigma,
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
			    
		    if(curClass == prevClass)
		    {
			    retCode = ARMLOCAL_GTWOYC (C_dMeanRevSpeed, C_dSigma,
					                  LocalGetNumObjectId (C_dZcId),
							          LocalGetNumObjectId (C_fZcId),
							          C_ratesCorr,
							          C_fMeanRevSpeed,
							          C_fSigma,
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
			    retCode = ARMLOCAL_GTWOYC (C_dMeanRevSpeed, C_dSigma,
					                  LocalGetNumObjectId (C_dZcId),
							          LocalGetNumObjectId (C_fZcId),
							          C_ratesCorr,
							          C_fMeanRevSpeed,
							          C_fSigma,
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

    //	ARM_END();
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_GTWOYC" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}




__declspec(dllexport) LPXLOPER WINAPI Local_PXL_GTWOYC (LPXLOPER XL_dMeanRevSpeed,
														LPXLOPER XL_dSigma,
														LPXLOPER XL_dZcId,
														LPXLOPER XL_fZcId,
														LPXLOPER XL_ratesCorr,
														LPXLOPER XL_fMeanRevSpeed,
														LPXLOPER XL_fSigma)
{
	ADD_LOG("Local_PXL_GTWOYC ");
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	    ARM_NOCALCIFWIZ();

	    // C variable
	    double C_dMeanRevSpeed;
	    double C_dSigma;

	    CCString C_dZcId;
	    CCString C_fZcId;

	    double C_ratesCorr;
	    double C_fMeanRevSpeed;
	    double C_fSigma;
	    
	    // error
	    static int error;
	    static char* reason = "";

	    XL_readNumCell(XL_dMeanRevSpeed,C_dMeanRevSpeed," ARM_ERR: domestic mean reversion speed: numeric expected",C_result);
	    XL_readNumCell(XL_dSigma,C_dSigma," ARM_ERR: domestic sigma: numeric expected",C_result);
	    XL_readStrCell(XL_dZcId,C_dZcId," ARM_ERR: domestic zero-coupon curve id: object expected",C_result);
	    XL_readStrCell(XL_fZcId,C_fZcId," ARM_ERR: foreign zero-coupon curve id: object expected",C_result);
	    XL_readNumCell(XL_ratesCorr,C_ratesCorr," ARM_ERR: rates correlation: numeric expected",C_result);
	    XL_readNumCell(XL_fMeanRevSpeed,C_fMeanRevSpeed," ARM_ERR: foreign mean reversion speed: numeric expected",C_result);
	    XL_readNumCell(XL_fSigma,C_fSigma," ARM_ERR: foreign sigma: numeric expected",C_result);
	    
	    long retCode;
	    long objId;
	    
	    CCString curClass = LOCAL_GAUSSIAN_2YC_MODEL_CLASS;
	    CCString stringId;

	    retCode = ARMLOCAL_GTWOYC (C_dMeanRevSpeed, C_dSigma,
								    LocalGetNumObjectId (C_dZcId),
								    LocalGetNumObjectId (C_fZcId),
								    C_ratesCorr,
								    C_fMeanRevSpeed,
								    C_fSigma,
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_PXL_GTWOYC" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}




__declspec(dllexport) LPXLOPER WINAPI Local_GYCMODEL(LPXLOPER XL_zcId,
													 LPXLOPER XL_a,
													 LPXLOPER XL_sigma)
{
	ADD_LOG("Local_GYCMODEL");
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	    ARM_NOCALCIFWIZ();

	    // C variable
	    CCString C_zcId;

	    double C_a;
	    double C_sigma;

	    // error
	    static int error;
	    static char* reason = "";

	    XL_readStrCell(XL_zcId, C_zcId, " ARM_ERR: curve id: object expected",C_result);
	    XL_readNumCell(XL_a, C_a," ARM_ERR: a: numeric expected",C_result);
	    XL_readNumCell(XL_sigma,C_sigma," ARM_ERR: sigma: numeric expected", C_result);
   

	    long retCode;
	    long objId;
	    CCString prevClass;

	    CCString curClass = LOCAL_GYCMODEL_CLASS;
	    CCString stringId = GetLastCurCellEnvValue ();

	    if (!stringId)
	    {
            retCode = ARMLOCAL_GYCMODEL(LocalGetNumObjectId(C_zcId),
                                   C_a, C_sigma, C_result);

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
               retCode = ARMLOCAL_GYCMODEL(LocalGetNumObjectId (C_zcId),
                                      C_a, C_sigma, C_result, objId);

			    if ( retCode == ARM_OK )
			    {
			       LocalSetCurCellEnvValue (curClass, objId); 

			       stringId = LocalMakeObjectId (objId, curClass);
			    }
		    }
		    else
		    {
			    FreeCurCellContent ();

                retCode = ARMLOCAL_GYCMODEL(LocalGetNumObjectId (C_zcId),
                                       C_a, C_sigma, C_result);
	    
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_GYCMODEL" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}



__declspec(dllexport) LPXLOPER WINAPI Local_PXL_GYCMODEL(LPXLOPER XL_zcId,
														 LPXLOPER XL_a,
														 LPXLOPER XL_sigma)
{
	ADD_LOG("Local_PXL_GYCMODEL");
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	    ARM_NOCALCIFWIZ();

	    // C variable
	    CCString C_zcId;

	    double C_a;
	    double C_sigma;

	    // error
	    static int error;
	    static char* reason = "";

	    XL_readStrCell(XL_zcId, C_zcId, " ARM_ERR: curve id: object expected",C_result);
	    XL_readNumCell(XL_a, C_a," ARM_ERR: a: numeric expected",C_result);
	    XL_readNumCell(XL_sigma,C_sigma," ARM_ERR: sigma: numeric expected", C_result);

	    long retCode;
	    long objId;

	    CCString curClass = LOCAL_GYCMODEL_CLASS;
	    CCString stringId;

        retCode = ARMLOCAL_GYCMODEL(LocalGetNumObjectId(C_zcId),
                               C_a, C_sigma, C_result);

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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_PXL_GYCMODEL" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}



__declspec(dllexport) LPXLOPER WINAPI Local_HWFNMONTECARLOSV (LPXLOPER XL_zc,
															  LPXLOPER XL_horizon,
															  LPXLOPER XL_a,
															  LPXLOPER XL_sigmaDate,
															  LPXLOPER XL_sigmaVal,
															  LPXLOPER XL_nbTraj)
{
	ADD_LOG("Local_HWFNMONTECARLOSV ");
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	    ARM_NOCALCIFWIZ();

	    // C variable
	    CCString C_zc;
	    double C_horizon;
	    double C_a;
	    VECTOR<double> C_sigmaDate;
	    VECTOR<double> C_sigmaVal;
	    double C_nbTraj;
	    
	    // error
	    static int error;
	    static char* reason = "";

	    XL_readStrCell(XL_zc,C_zc," ARM_ERR: zero-coupon curve id: object expected",C_result);
	    XL_readNumCell(XL_horizon,C_horizon," ARM_ERR: horizon date: date expected",C_result);
	    XL_readNumCell(XL_a,C_a," ARM_ERR: a: numeric expected",C_result);
	    XL_readNumVector(XL_sigmaDate,C_sigmaDate," ARM_ERR: sigma dates: array of date expected",C_result);
	    XL_readNumVector(XL_sigmaVal,C_sigmaVal," ARM_ERR: sigma values: array of numeric expected",C_result);
	    XL_readNumCell(XL_nbTraj,C_nbTraj," ARM_ERR: number of trajects: numeric expected",C_result);
	    
	    long retCode;
	    long objId;
	    CCString prevClass;
	    
	    CCString curClass = LOCAL_HWFNMONTECARLOSV_CLASS;
	    CCString stringId = GetLastCurCellEnvValue ();

	    if (!stringId)
	    {
	       retCode = ARMLOCAL_HWFNMONTECARLOSV (LocalGetNumObjectId (C_zc),
									       C_horizon,
									       C_a,
									       C_sigmaDate,
									       C_sigmaVal,
									       (long)C_nbTraj,
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
			    
		    if(curClass == prevClass)
		    {
			    retCode = ARMLOCAL_HWFNMONTECARLOSV (LocalGetNumObjectId (C_zc),
									            C_horizon,
									            C_a,
											    C_sigmaDate,
											    C_sigmaVal,
											    (long)C_nbTraj,
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
			    retCode = ARMLOCAL_HWFNMONTECARLOSV (LocalGetNumObjectId (C_zc),
											    C_horizon,
											    C_a,
											    C_sigmaDate,
											    C_sigmaVal,
											    (long)C_nbTraj,
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

   //	ARM_END();
    }
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_HWFNMONTECARLOSV" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}



__declspec(dllexport) LPXLOPER WINAPI Local_PXL_HWFNMONTECARLOSV (LPXLOPER XL_zc,
																  LPXLOPER XL_horizon,
																  LPXLOPER XL_a,
																  LPXLOPER XL_sigmaDate,
																  LPXLOPER XL_sigmaVal,
																  LPXLOPER XL_nbTraj)
{
	ADD_LOG("Local_PXL_HWFNMONTECARLOSV ");
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	    ARM_NOCALCIFWIZ();

	    // C variable
	    CCString C_zc;
	    double C_horizon;
	    double C_a;
	    VECTOR<double> C_sigmaDate;
	    VECTOR<double> C_sigmaVal;
	    double C_nbTraj;
	    
	    // error
	    static int error;
	    static char* reason = "";

	    XL_readStrCell(XL_zc,C_zc," ARM_ERR: zero-coupon curve id: object expected",C_result);
	    XL_readNumCell(XL_horizon,C_horizon," ARM_ERR: horizon date: date expected",C_result);
	    XL_readNumCell(XL_a,C_a," ARM_ERR: a: numeric expected",C_result);
	    XL_readNumVector(XL_sigmaDate,C_sigmaDate," ARM_ERR: sigma dates: array of date expected",C_result);
	    XL_readNumVector(XL_sigmaVal,C_sigmaVal," ARM_ERR: sigma values: array of numeric expected",C_result);
	    XL_readNumCell(XL_nbTraj,C_nbTraj," ARM_ERR: number of trajects: numeric expected",C_result);
	    
	    long retCode;
	    long objId;
	    
	    CCString curClass = LOCAL_HWFNMONTECARLOSV_CLASS;
	    CCString stringId;

	    retCode = ARMLOCAL_HWFNMONTECARLOSV (LocalGetNumObjectId (C_zc),
								        C_horizon,
									    C_a,
									    C_sigmaDate,
									    C_sigmaVal,
									    (long)C_nbTraj,
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
	       PXL_ARM_ERR();
	    }

    //	ARM_END();
    }
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_PXL_HWFNMONTECARLOSV" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}



__declspec(dllexport) LPXLOPER WINAPI Local_HWFNMONTECARLO (LPXLOPER XL_zc,
															LPXLOPER XL_horizon,
															LPXLOPER XL_a,
															LPXLOPER XL_sigma,
															LPXLOPER XL_nbTraj)
{
	ADD_LOG("Local_HWFNMONTECARLO ");
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	    ARM_NOCALCIFWIZ();

	    // C variable
	    CCString C_zc;
	    double C_horizon;
	    double C_a;
	    double C_sigma;
	    double C_nbTraj;
	    
	    // error
	    static int error;
	    static char* reason = "";

	    XL_readStrCell(XL_zc,C_zc," ARM_ERR: zero-coupon curve id: object expected",C_result);
	    XL_readNumCell(XL_horizon,C_horizon," ARM_ERR: horizon date: date expected",C_result);
	    XL_readNumCell(XL_a,C_a," ARM_ERR: a: numeric expected",C_result);
	    XL_readNumCell(XL_sigma,C_sigma," ARM_ERR: sigma: numeric expected",C_result);
	    XL_readNumCell(XL_nbTraj,C_nbTraj," ARM_ERR: number of trajects: numeric expected",C_result);
	    
	    long retCode;
	    long objId;
	    CCString prevClass;
	    
	    CCString curClass = LOCAL_HWFNMONTECARLO_CLASS;
	    CCString stringId = GetLastCurCellEnvValue ();

	    if (!stringId)
	    {
		    retCode = ARMLOCAL_HWFNMONTECARLO (LocalGetNumObjectId (C_zc),
										       C_horizon,
										       C_a,
										       C_sigma,
										       (long)C_nbTraj,
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
			    
		    if(curClass == prevClass)
		    {
			    retCode = ARMLOCAL_HWFNMONTECARLO (LocalGetNumObjectId (C_zc),
											       C_horizon,
											       C_a,
											       C_sigma,
											       (long)C_nbTraj,
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
			    retCode = ARMLOCAL_HWFNMONTECARLO (LocalGetNumObjectId (C_zc),
											       C_horizon,
											       C_a,
											       C_sigma,
											       (long)C_nbTraj,
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

    //	ARM_END();	
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_HWFNMONTECARLO" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}



__declspec(dllexport) LPXLOPER WINAPI Local_PXL_HWFNMONTECARLO (LPXLOPER XL_zc,
																LPXLOPER XL_horizon,
																LPXLOPER XL_a,
																LPXLOPER XL_sigma,
																LPXLOPER XL_nbTraj)
{
	ADD_LOG("Local_PXL_HWFNMONTECARLO ");
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	    ARM_NOCALCIFWIZ();

	    // C variable
	    CCString C_zc;
	    double C_horizon;
	    double C_a;
	    double C_sigma;
	    double C_nbTraj;

	    // error
	    static int error;
	    static char* reason = "";

	    XL_readStrCell(XL_zc,C_zc," ARM_ERR: zero-coupon curve id: object expected",C_result);
	    XL_readNumCell(XL_horizon,C_horizon," ARM_ERR: horizon date: date expected",C_result);
	    XL_readNumCell(XL_a,C_a," ARM_ERR: a: numeric expected",C_result);
	    XL_readNumCell(XL_sigma,C_sigma," ARM_ERR: sigma: numeric expected",C_result);
	    XL_readNumCell(XL_nbTraj,C_nbTraj," ARM_ERR: number of trajects: numeric expected",C_result);

	    long retCode;
	    long objId;

	    CCString curClass = LOCAL_HWFNMONTECARLO_CLASS;
	    CCString stringId;

	    retCode = ARMLOCAL_HWFNMONTECARLO (LocalGetNumObjectId (C_zc),
									       C_horizon,
									       C_a,
									       C_sigma,
									       (long)C_nbTraj,
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
	       PXL_ARM_ERR();
	    }

        //	ARM_END();    
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_PXL_HWFNMONTECARLO" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}



__declspec(dllexport) LPXLOPER WINAPI Local_FRMTREE_AUTO(LPXLOPER XL_zcId,
														 LPXLOPER XL_volId,
														 LPXLOPER XL_smileId,
														 LPXLOPER XL_autoMode,
														 LPXLOPER XL_horizon,
														 LPXLOPER XL_fineMonth,
														 LPXLOPER XL_shapeDecay,
														 LPXLOPER XL_shapeSlope,
														 LPXLOPER XL_shapeAsymptote,
														 LPXLOPER XL_nbFactor,
														 LPXLOPER XL_corrMatu,
														 LPXLOPER XL_corrMatrix,
														 LPXLOPER XL_corr2Matu)
{
	ADD_LOG("Local_FRMTREE_AUTO");
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	    ARM_NOCALCIFWIZ();

	    // C variable
	    CCString C_zcId;
        CCString C_volId;
        CCString C_smileId;
        CCString CS_autoMode;
        CCString CS_fineMonth;
        double C_horizon;
        double C_shapeDecay;
        double C_shapeSlope;
        double C_shapeAsymptote;
        double C_nbFactor;
   
        VECTOR<double> C_corrMatu;
        VECTOR<double> C_corrMatrix;
        VECTOR<double> C_corr2Matu;	

        double mFineDefault = K_L_DEFAULT;
        double epsDefault = 1e-8;
        double huegeDefault = 1e8;
        double zero = 0;
        double one = 1;
	    
	    // error
	    static int error;
	    static char* reason = "";

        XL_readStrCell(XL_zcId, C_zcId, " ARM_ERR: Curve id: string expected",C_result);
        XL_readStrCell(XL_volId, C_volId, " ARM_ERR: Vol id: string expected",C_result);
        XL_readStrCell(XL_smileId, C_smileId, " ARM_ERR: Smile id: string expected",C_result);
        XL_readStrCell(XL_autoMode, CS_autoMode, " ARM_ERR: AutoMode: string expected",C_result);
	    XL_readNumCell(XL_horizon,C_horizon," ARM_ERR: Horizon date: date expected",C_result);
        XL_readStrCellWD(XL_fineMonth, CS_fineMonth, "L_DEFAULT", " ARM_ERR: fine month. : string expected",C_result);
        XL_readNumCellWD(XL_shapeDecay, C_shapeDecay, zero, " ARM_ERR: shape decay. : numeric expected",C_result);
        XL_readNumCellWD(XL_shapeSlope, C_shapeSlope, zero, " ARM_ERR: shape slope. : numeric expected",C_result);
        XL_readNumCellWD(XL_shapeAsymptote, C_shapeAsymptote, zero, " ARM_ERR: shape Asymptote. : numeric expected",C_result);
        XL_readNumCellWD(XL_nbFactor, C_nbFactor, one, " ARM_ERR: Nb factor. : numeric expected",C_result);

   	    if((XL_corrMatu->xltype == xltypeMissing) || (XL_corrMatu->xltype == xltypeNil))
	    {
        }
	    else
	    {
            XL_readNumVector(XL_corrMatu, C_corrMatu, " ARM_ERR: Correlation maturities: array of numeric expected",C_result);
        }

   	    if((XL_corrMatrix->xltype == xltypeMissing) || (XL_corrMatrix->xltype == xltypeNil))
	    {
        }
	    else
	    {
            XL_readNumVector(XL_corrMatrix, C_corrMatrix, " ARM_ERR: Correlation matrix: array of numeric expected",C_result);
        }

   	    if((XL_corr2Matu->xltype == xltypeMissing) || (XL_corr2Matu->xltype == xltypeNil))
	    {
        }
	    else
	    {
            XL_readNumVector(XL_corr2Matu, C_corr2Matu, " ARM_ERR: Correlation correlated maturities : array of numeric expected",C_result);
        }

	    long autoModeId;
	    long fineMonthId;

	    if((autoModeId = ARM_ConvAutoMode (CS_autoMode, C_result)) == ARM_DEFAULT_ERR)
	    {
		    ARM_ARG_ERR();
		    return (LPXLOPER)&XL_result;
	    }

	    if((fineMonthId = ARM_ConvFineMode (CS_fineMonth, C_result)) == ARM_DEFAULT_ERR)
	    {
		    ARM_ARG_ERR();
		    return (LPXLOPER)&XL_result;
	    }

	    long retCode;
	    long objId;
	    CCString prevClass;

	    CCString curClass = LOCAL_FRMTREE_CLASS;
	    CCString stringId = GetLastCurCellEnvValue ();

	    if (!stringId)
	    {
            retCode = ARMLOCAL_FRMTREE_AUTO(LocalGetNumObjectId(C_zcId), LocalGetNumObjectId(C_volId), LocalGetNumObjectId(C_smileId), 
                                       autoModeId,  C_horizon, fineMonthId, C_shapeDecay, C_shapeSlope, C_shapeAsymptote,
                                       (long) C_nbFactor, C_corrMatu, C_corrMatrix, C_corr2Matu, C_result);

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
               retCode = ARMLOCAL_FRMTREE_AUTO(LocalGetNumObjectId(C_zcId), LocalGetNumObjectId(C_volId), LocalGetNumObjectId(C_smileId),
                                       autoModeId, C_horizon, fineMonthId, C_shapeDecay, C_shapeSlope, C_shapeAsymptote,
                                       (long) C_nbFactor, C_corrMatu, C_corrMatrix, C_corr2Matu, C_result, objId);

			    if ( retCode == ARM_OK )
			    {
			       LocalSetCurCellEnvValue (curClass, objId); 

			       stringId = LocalMakeObjectId (objId, curClass);
			    }
		    }
		    else
		    {
			    FreeCurCellContent ();

                retCode = ARMLOCAL_FRMTREE_AUTO(LocalGetNumObjectId(C_zcId), LocalGetNumObjectId(C_volId), LocalGetNumObjectId(C_smileId),
                                       autoModeId, C_horizon, fineMonthId, C_shapeDecay, C_shapeSlope, C_shapeAsymptote,
                                       (long) C_nbFactor, C_corrMatu, C_corrMatrix, C_corr2Matu, C_result);
	    
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_FRMTREE_AUTO" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}



__declspec(dllexport) LPXLOPER WINAPI Local_PXL_FRMTREE_AUTO(LPXLOPER XL_zcId,
															 LPXLOPER XL_volId,
															 LPXLOPER XL_smileId,
															 LPXLOPER XL_autoMode,
															 LPXLOPER XL_horizon,
															 LPXLOPER XL_fineMonth,
															 LPXLOPER XL_shapeDecay,
															 LPXLOPER XL_shapeSlope,
															 LPXLOPER XL_shapeAsymptote,
															 LPXLOPER XL_nbFactor,
															 LPXLOPER XL_corrMatu,
															 LPXLOPER XL_corrMatrix,
															 LPXLOPER XL_corr2Matu)
{
	ADD_LOG("Local_PXL_FRMTREE_AUTO");
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();

	// C variable
	CCString C_zcId;
    CCString C_volId;
    CCString C_smileId;
    CCString CS_autoMode;
    CCString CS_fineMonth;
    double C_horizon;
    double C_shapeDecay;
    double C_shapeSlope;
    double C_shapeAsymptote;
    double C_nbFactor;
   
    VECTOR<double> C_corrMatu;
    VECTOR<double> C_corrMatrix;
    VECTOR<double> C_corr2Matu;	

    double mFineDefault = K_L_DEFAULT;
    double epsDefault = 1e-8;
    double huegeDefault = 1e8;
    double zero = 0;
    double one = 1;
	
	// error
	static int error;
	static char* reason = "";

    XL_readStrCell(XL_zcId, C_zcId, " ARM_ERR: Curve id: string expected",C_result);
    XL_readStrCell(XL_volId, C_volId, " ARM_ERR: Vol id: string expected",C_result);
    XL_readStrCell(XL_smileId, C_smileId, " ARM_ERR: Smile id: string expected",C_result);
    XL_readStrCell(XL_autoMode, CS_autoMode, " ARM_ERR: AutoMode: string expected",C_result);
	XL_readNumCell(XL_horizon,C_horizon," ARM_ERR: Horizon date: date expected",C_result);
    XL_readStrCellWD(XL_fineMonth, CS_fineMonth, "L_DEFAULT", " ARM_ERR: fine month. : string expected",C_result);
    XL_readNumCellWD(XL_shapeDecay, C_shapeDecay, zero, " ARM_ERR: shape decay. : numeric expected",C_result);
    XL_readNumCellWD(XL_shapeSlope, C_shapeSlope, zero, " ARM_ERR: shape slope. : numeric expected",C_result);
    XL_readNumCellWD(XL_shapeAsymptote, C_shapeAsymptote, zero, " ARM_ERR: shape Asymptote. : numeric expected",C_result);
    XL_readNumCellWD(XL_nbFactor, C_nbFactor, one, " ARM_ERR: Nb factor. : numeric expected",C_result);

   	if((XL_corrMatu->xltype == xltypeMissing) || (XL_corrMatu->xltype == xltypeNil))
	{
    }
	else
	{
        XL_readNumVector(XL_corrMatu, C_corrMatu, " ARM_ERR: Correlation maturities: array of numeric expected",C_result);
    }

   	if((XL_corrMatrix->xltype == xltypeMissing) || (XL_corrMatrix->xltype == xltypeNil))
	{
    }
	else
	{
        XL_readNumVector(XL_corrMatrix, C_corrMatrix, " ARM_ERR: Correlation matrix: array of numeric expected",C_result);
    }

   	if((XL_corr2Matu->xltype == xltypeMissing) || (XL_corr2Matu->xltype == xltypeNil))
	{
    }
	else
	{
        XL_readNumVector(XL_corr2Matu, C_corr2Matu, " ARM_ERR: Correlation correlated maturities : array of numeric expected",C_result);
    }

	long autoModeId;
	long fineMonthId;

	if((autoModeId = ARM_ConvAutoMode (CS_autoMode, C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	if((fineMonthId = ARM_ConvFineMode (CS_fineMonth, C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	long retCode;
	long objId;

	CCString curClass = LOCAL_FRMTREE_CLASS;
	CCString stringId;

    retCode = ARMLOCAL_FRMTREE_AUTO(LocalGetNumObjectId(C_zcId), LocalGetNumObjectId(C_volId), LocalGetNumObjectId(C_smileId), 
                               autoModeId,  C_horizon, fineMonthId, C_shapeDecay, C_shapeSlope, C_shapeAsymptote,
                               (long) C_nbFactor, C_corrMatu, C_corrMatrix, C_corr2Matu, C_result);

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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_PXL_FRMTREE_AUTO" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}



__declspec(dllexport) LPXLOPER WINAPI Local_FRMTREE_AUTO_G(LPXLOPER XL_zcId,
														   LPXLOPER XL_zdId,
														   LPXLOPER XL_volId,
														   LPXLOPER XL_smileId,
														   LPXLOPER XL_irgvolId,
														   LPXLOPER XL_irgsmileId,
														   LPXLOPER XL_autoMode,
														   LPXLOPER XL_horizon,
														   LPXLOPER XL_fineMonth,
														   LPXLOPER XL_shapeDecay,
														   LPXLOPER XL_shapeSlope,
														   LPXLOPER XL_shapeAsymptote,
														   LPXLOPER XL_nbFactor,
														   LPXLOPER XL_corrMatu,
														   LPXLOPER XL_corrMatrix,
														   LPXLOPER XL_corr2Matu)
{
	ADD_LOG("Local_FRMTREE_AUTO_G");
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();

	// C variable
	CCString C_zcId;
	CCString C_zdId;

    CCString C_volId;
    CCString C_smileId;
    CCString C_irgvolId;
    CCString C_irgsmileId;
    CCString CS_autoMode;
    CCString CS_fineMonth;
    double C_horizon;
    double C_shapeDecay;
    double C_shapeSlope;
    double C_shapeAsymptote;
    double C_nbFactor;
   
    VECTOR<double> C_corrMatu;
    VECTOR<double> C_corrMatrix;
    VECTOR<double> C_corr2Matu;	

    double mFineDefault = K_L_DEFAULT;
    double epsDefault = 1e-8;
    double huegeDefault = 1e8;
    double zero = 0;
    double one = 1;
	
	// error
	static int error;
	static char* reason = "";

    XL_readStrCell(XL_zcId, C_zcId, " ARM_ERR: Curve id: string expected",C_result);
	XL_readStrCell(XL_zdId, C_zdId, " ARM_ERR: Curve id: string expected",C_result);
    XL_readStrCell(XL_volId, C_volId, " ARM_ERR: Vol id: string expected",C_result);
    XL_readStrCell(XL_smileId, C_smileId, " ARM_ERR: Smile id: string expected",C_result);
	XL_readStrCell(XL_irgvolId, C_irgvolId, " ARM_ERR: Irg Vol id: string expected",C_result);
    XL_readStrCell(XL_irgsmileId, C_irgsmileId, " ARM_ERR: Irg Smile id: string expected",C_result);
    XL_readStrCell(XL_autoMode, CS_autoMode, " ARM_ERR: AutoMode: string expected",C_result);
	XL_readNumCell(XL_horizon,C_horizon," ARM_ERR: Horizon date: date expected",C_result);
    XL_readStrCellWD(XL_fineMonth, CS_fineMonth, "L_DEFAULT", " ARM_ERR: fine month. : string expected",C_result);
    XL_readNumCellWD(XL_shapeDecay, C_shapeDecay, zero, " ARM_ERR: shape decay. : numeric expected",C_result);
    XL_readNumCellWD(XL_shapeSlope, C_shapeSlope, zero, " ARM_ERR: shape slope. : numeric expected",C_result);
    XL_readNumCellWD(XL_shapeAsymptote, C_shapeAsymptote, zero, " ARM_ERR: shape Asymptote. : numeric expected",C_result);
    XL_readNumCellWD(XL_nbFactor, C_nbFactor, one, " ARM_ERR: Nb factor. : numeric expected",C_result);

   	if((XL_corrMatu->xltype == xltypeMissing) || (XL_corrMatu->xltype == xltypeNil))
	{
    }
	else
	{
        XL_readNumVector(XL_corrMatu, C_corrMatu, " ARM_ERR: Correlation maturities: array of numeric expected",C_result);
    }

   	if((XL_corrMatrix->xltype == xltypeMissing) || (XL_corrMatrix->xltype == xltypeNil))
	{
    }
	else
	{
        XL_readNumVector(XL_corrMatrix, C_corrMatrix, " ARM_ERR: Correlation matrix: array of numeric expected",C_result);
    }

   	if((XL_corr2Matu->xltype == xltypeMissing) || (XL_corr2Matu->xltype == xltypeNil))
	{
    }
	else
	{
        XL_readNumVector(XL_corr2Matu, C_corr2Matu, " ARM_ERR: Correlation correlated maturities : array of numeric expected",C_result);
    }

	long autoModeId;
	long fineMonthId;

	if((autoModeId = ARM_ConvAutoMode2 (CS_autoMode, C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	if((fineMonthId = ARM_ConvFineMode (CS_fineMonth, C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	long retCode;
	long objId;
	CCString prevClass;

	CCString curClass = LOCAL_FRMTREE_CLASS;
	CCString stringId = GetLastCurCellEnvValue ();

	if (!stringId)
	{
        retCode = ARMLOCAL_FRMTREE_AUTO_G(LocalGetNumObjectId(C_zcId), LocalGetNumObjectId(C_zdId),
									LocalGetNumObjectId(C_volId), LocalGetNumObjectId(C_smileId),
									LocalGetNumObjectId(C_irgvolId), LocalGetNumObjectId(C_irgsmileId),
                                   autoModeId,  C_horizon, fineMonthId, C_shapeDecay, C_shapeSlope, C_shapeAsymptote,
                                   (long) C_nbFactor, C_corrMatu, C_corrMatrix, C_corr2Matu, C_result);

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
           retCode = ARMLOCAL_FRMTREE_AUTO_G(LocalGetNumObjectId(C_zcId), LocalGetNumObjectId(C_zdId),
									LocalGetNumObjectId(C_volId), LocalGetNumObjectId(C_smileId),
									LocalGetNumObjectId(C_irgvolId), LocalGetNumObjectId(C_irgsmileId),
                                   autoModeId, C_horizon, fineMonthId, C_shapeDecay, C_shapeSlope, C_shapeAsymptote,
                                   (long) C_nbFactor, C_corrMatu, C_corrMatrix, C_corr2Matu, C_result, objId);

			if ( retCode == ARM_OK )
			{
			   LocalSetCurCellEnvValue (curClass, objId); 

			   stringId = LocalMakeObjectId (objId, curClass);
			}
		}
		else
		{
			FreeCurCellContent ();

            retCode = ARMLOCAL_FRMTREE_AUTO_G(LocalGetNumObjectId(C_zcId), LocalGetNumObjectId(C_zdId),
									LocalGetNumObjectId(C_volId), LocalGetNumObjectId(C_smileId),
									LocalGetNumObjectId(C_irgvolId), LocalGetNumObjectId(C_irgsmileId),
                                   autoModeId, C_horizon, fineMonthId, C_shapeDecay, C_shapeSlope, C_shapeAsymptote,
                                   (long) C_nbFactor, C_corrMatu, C_corrMatrix, C_corr2Matu, C_result);
	
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_FRMTREE_AUTO_G" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


__declspec(dllexport) LPXLOPER WINAPI Local_PXL_FRMTREE_AUTO_G(LPXLOPER XL_zcId,
															   LPXLOPER XL_zdId,
															   LPXLOPER XL_volId,
															   LPXLOPER XL_smileId,
															   LPXLOPER XL_irgvolId,
															   LPXLOPER XL_irgsmileId,
															   LPXLOPER XL_autoMode,
															   LPXLOPER XL_horizon,
															   LPXLOPER XL_fineMonth,
															   LPXLOPER XL_shapeDecay,
															   LPXLOPER XL_shapeSlope,
															   LPXLOPER XL_shapeAsymptote,
															   LPXLOPER XL_nbFactor,
															   LPXLOPER XL_corrMatu,
															   LPXLOPER XL_corrMatrix,
															   LPXLOPER XL_corr2Matu)
{
	ADD_LOG("Local_PXL_FRMTREE_AUTO_G");
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();

	// C variable
	CCString C_zcId;
	CCString C_zdId;

    CCString C_volId;
    CCString C_smileId;
    CCString C_irgvolId;
    CCString C_irgsmileId;
    CCString CS_autoMode;
    CCString CS_fineMonth;
    double C_horizon;
    double C_shapeDecay;
    double C_shapeSlope;
    double C_shapeAsymptote;
    double C_nbFactor;
   
    VECTOR<double> C_corrMatu;
    VECTOR<double> C_corrMatrix;
    VECTOR<double> C_corr2Matu;	

    double mFineDefault = K_L_DEFAULT;
    double epsDefault = 1e-8;
    double huegeDefault = 1e8;
    double zero = 0;
    double one = 1;
	
	// error
	static int error;
	static char* reason = "";

    XL_readStrCell(XL_zcId, C_zcId, " ARM_ERR: Curve id: string expected",C_result);
	XL_readStrCell(XL_zdId, C_zdId, " ARM_ERR: Curve id: string expected",C_result);
    XL_readStrCell(XL_volId, C_volId, " ARM_ERR: Vol id: string expected",C_result);
    XL_readStrCell(XL_smileId, C_smileId, " ARM_ERR: Smile id: string expected",C_result);
	XL_readStrCell(XL_irgvolId, C_irgvolId, " ARM_ERR: Irg Vol id: string expected",C_result);
    XL_readStrCell(XL_irgsmileId, C_irgsmileId, " ARM_ERR: Irg Smile id: string expected",C_result);
    XL_readStrCell(XL_autoMode, CS_autoMode, " ARM_ERR: AutoMode: string expected",C_result);
	XL_readNumCell(XL_horizon,C_horizon," ARM_ERR: Horizon date: date expected",C_result);
    XL_readStrCellWD(XL_fineMonth, CS_fineMonth, "L_DEFAULT", " ARM_ERR: fine month. : string expected",C_result);
    XL_readNumCellWD(XL_shapeDecay, C_shapeDecay, zero, " ARM_ERR: shape decay. : numeric expected",C_result);
    XL_readNumCellWD(XL_shapeSlope, C_shapeSlope, zero, " ARM_ERR: shape slope. : numeric expected",C_result);
    XL_readNumCellWD(XL_shapeAsymptote, C_shapeAsymptote, zero, " ARM_ERR: shape Asymptote. : numeric expected",C_result);
    XL_readNumCellWD(XL_nbFactor, C_nbFactor, one, " ARM_ERR: Nb factor. : numeric expected",C_result);

   	if((XL_corrMatu->xltype == xltypeMissing) || (XL_corrMatu->xltype == xltypeNil))
	{
    }
	else
	{
        XL_readNumVector(XL_corrMatu, C_corrMatu, " ARM_ERR: Correlation maturities: array of numeric expected",C_result);
    }

   	if((XL_corrMatrix->xltype == xltypeMissing) || (XL_corrMatrix->xltype == xltypeNil))
	{
    }
	else
	{
        XL_readNumVector(XL_corrMatrix, C_corrMatrix, " ARM_ERR: Correlation matrix: array of numeric expected",C_result);
    }

   	if((XL_corr2Matu->xltype == xltypeMissing) || (XL_corr2Matu->xltype == xltypeNil))
	{
    }
	else
	{
        XL_readNumVector(XL_corr2Matu, C_corr2Matu, " ARM_ERR: Correlation correlated maturities : array of numeric expected",C_result);
    }

	long autoModeId;
	long fineMonthId;

	if((autoModeId = ARM_ConvAutoMode2 (CS_autoMode, C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	if((fineMonthId = ARM_ConvFineMode (CS_fineMonth, C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	long retCode;
	long objId;

	CCString curClass = LOCAL_FRMTREE_CLASS;
	CCString stringId;

    retCode = ARMLOCAL_FRMTREE_AUTO_G(LocalGetNumObjectId(C_zcId), LocalGetNumObjectId(C_zdId),
								LocalGetNumObjectId(C_volId), LocalGetNumObjectId(C_smileId),
								LocalGetNumObjectId(C_irgvolId), LocalGetNumObjectId(C_irgsmileId),
                               autoModeId,  C_horizon, fineMonthId, C_shapeDecay, C_shapeSlope, C_shapeAsymptote,
                               (long) C_nbFactor, C_corrMatu, C_corrMatrix, C_corr2Matu, C_result);

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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_PXL_FRMTREE_AUTO_G" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}



__declspec(dllexport) LPXLOPER WINAPI Local_DFHWSIGVARTREE (LPXLOPER XL_startDate,
															LPXLOPER XL_horizon,
															LPXLOPER XL_numSteps,
															LPXLOPER XL_dZc,
															LPXLOPER XL_fZc,
															LPXLOPER XL_dMeanRevSpeed,
															LPXLOPER XL_dDate,
															LPXLOPER XL_dSigma,
															LPXLOPER XL_fMeanRevSpeed,
															LPXLOPER XL_fDate,
															LPXLOPER XL_fSigma,
															LPXLOPER XL_dFxCorr,
															LPXLOPER XL_fFxCorr,
															LPXLOPER XL_fxVol,
															LPXLOPER XL_ratesCorr,
															LPXLOPER XL_fxSpotRate)
{
	ADD_LOG("Local_DFHWSIGVARTREE ");
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();

	// C variable
	double C_startDate;
	double C_horizon;
	double C_numSteps;
	CCString C_dZc;
	CCString C_fZc;
	double C_dMeanRevSpeed;
	VECTOR<double> C_dDate;
	VECTOR<double> C_dSigma;
	double C_fMeanRevSpeed;
	VECTOR<double> C_fDate;
	VECTOR<double> C_fSigma;
	CCString C_dFxCorr;
	CCString C_fFxCorr;
	CCString C_fxVol;
	double C_ratesCorr;
	double C_fxSpotRate;

	// error
	static int error;
	static char* reason = "";

	XL_readNumCell(XL_startDate,C_startDate," ARM_ERR: start date: date expected",C_result);
	XL_readNumCell(XL_horizon,C_horizon," ARM_ERR: horizon date: date expected",C_result);
	XL_readNumCell(XL_numSteps,C_numSteps," ARM_ERR: number of steps: numeric expected",C_result);
	XL_readStrCell(XL_dZc,C_dZc," ARM_ERR: domestic zero-coupon curve id: object expected",C_result);
	XL_readStrCell(XL_fZc,C_fZc," ARM_ERR: foreign zero-coupon curve id: object expected",C_result);
	XL_readNumCell(XL_dMeanRevSpeed,C_dMeanRevSpeed," ARM_ERR: domestic mean reversion speed: numeric expected",C_result);
	XL_readNumVector (XL_dDate, C_dDate, " ARM_ERR: domestic date: array of date expected",C_result);
	XL_readNumVector (XL_dSigma, C_dSigma, " ARM_ERR: domestic sigma: array of numeric expected",C_result);
	XL_readNumCell(XL_fMeanRevSpeed,C_fMeanRevSpeed," ARM_ERR: foreign mean reversion speed: numeric expected",C_result);
	XL_readNumVector (XL_fDate, C_fDate, " ARM_ERR: foreign date: array of date expected",C_result);
	XL_readNumVector (XL_fSigma, C_fSigma, " ARM_ERR: foreign sigma: array of numeric expected",C_result);
	XL_readStrCell(XL_dFxCorr,C_dFxCorr," ARM_ERR: domestic forex correlation: object expected",C_result);
	XL_readStrCell(XL_fFxCorr,C_fFxCorr," ARM_ERR: foreign forex correlation: object expected",C_result);
	XL_readStrCell(XL_fxVol,C_fxVol," ARM_ERR: forex volatility: object expected",C_result);
	XL_readNumCell(XL_ratesCorr,C_ratesCorr," ARM_ERR: rates correlation: numeric expected",C_result);
	XL_readNumCell(XL_fxSpotRate,C_fxSpotRate," ARM_ERR: forex spot rates: numeric expected",C_result);
		
	long retCode;
	long objId;
	CCString prevClass;
	
	CCString curClass = LOCAL_DFHWSIGVARTREE_CLASS;
	CCString stringId = GetLastCurCellEnvValue ();

	if (!stringId)
	{
		retCode = ARMLOCAL_DFHWSIGVARTREE (C_startDate,
										   C_horizon,
										   (long)C_numSteps,
										   LocalGetNumObjectId (C_dZc),
										   LocalGetNumObjectId (C_fZc),
										   C_dMeanRevSpeed,
										   C_dDate,
										   C_dSigma,
										   C_fMeanRevSpeed,
										   C_fDate,
										   C_fSigma,
										   LocalGetNumObjectId (C_dFxCorr),
										   LocalGetNumObjectId (C_fFxCorr),
										   LocalGetNumObjectId (C_fxVol),
										   C_ratesCorr,
										   C_fxSpotRate,
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
			
		if(curClass == prevClass)
		{
			retCode = ARMLOCAL_DFHWSIGVARTREE (C_startDate,
										  C_horizon,
									      (long)C_numSteps,
									      LocalGetNumObjectId (C_dZc),
										  LocalGetNumObjectId (C_fZc),
										  C_dMeanRevSpeed,
										  C_dDate,
										  C_dSigma,
										  C_fMeanRevSpeed,
										  C_fDate,
										  C_fSigma,
										  LocalGetNumObjectId (C_dFxCorr),
										  LocalGetNumObjectId (C_fFxCorr),
										  LocalGetNumObjectId (C_fxVol),
										  C_ratesCorr,
										  C_fxSpotRate,
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
			retCode = ARMLOCAL_DFHWSIGVARTREE (C_startDate,
										  C_horizon,
										  (long)C_numSteps,
										  LocalGetNumObjectId (C_dZc),
										  LocalGetNumObjectId (C_fZc),
										  C_dMeanRevSpeed,
										  C_dDate,
										  C_dSigma,
										  C_fMeanRevSpeed,
										  C_fDate,
										  C_fSigma,
										  LocalGetNumObjectId (C_dFxCorr),
										  LocalGetNumObjectId (C_fFxCorr),
										  LocalGetNumObjectId (C_fxVol),
										  C_ratesCorr,
										  C_fxSpotRate,
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

//	ARM_END();
	
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_DFHWSIGVARTREE" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


__declspec(dllexport) LPXLOPER WINAPI Local_PXL_DFHWSIGVARTREE (LPXLOPER XL_startDate,
																LPXLOPER XL_horizon,
																LPXLOPER XL_numSteps,
																LPXLOPER XL_dZc,
																LPXLOPER XL_fZc,
																LPXLOPER XL_dMeanRevSpeed,
																LPXLOPER XL_dDate,
																LPXLOPER XL_dSigma,
																LPXLOPER XL_fMeanRevSpeed,
																LPXLOPER XL_fDate,
																LPXLOPER XL_fSigma,
																LPXLOPER XL_dFxCorr,
																LPXLOPER XL_fFxCorr,
																LPXLOPER XL_fxVol,
																LPXLOPER XL_ratesCorr,
																LPXLOPER XL_fxSpotRate)
{
	ADD_LOG("Local_PXL_DFHWSIGVARTREE ");
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();

	// C variable
	double C_startDate;
	double C_horizon;
	double C_numSteps;
	CCString C_dZc;
	CCString C_fZc;
	double C_dMeanRevSpeed;
	VECTOR<double> C_dDate;
	VECTOR<double> C_dSigma;
	double C_fMeanRevSpeed;
	VECTOR<double> C_fDate;
	VECTOR<double> C_fSigma;
	CCString C_dFxCorr;
	CCString C_fFxCorr;
	CCString C_fxVol;
	double C_ratesCorr;
	double C_fxSpotRate;

	// error
	static int error;
	static char* reason = "";

	XL_readNumCell(XL_startDate,C_startDate," ARM_ERR: start date: date expected",C_result);
	XL_readNumCell(XL_horizon,C_horizon," ARM_ERR: horizon date: date expected",C_result);
	XL_readNumCell(XL_numSteps,C_numSteps," ARM_ERR: number of steps: numeric expected",C_result);
	XL_readStrCell(XL_dZc,C_dZc," ARM_ERR: domestic zero-coupon curve id: object expected",C_result);
	XL_readStrCell(XL_fZc,C_fZc," ARM_ERR: foreign zero-coupon curve id: object expected",C_result);
	XL_readNumCell(XL_dMeanRevSpeed,C_dMeanRevSpeed," ARM_ERR: domestic mean reversion speed: numeric expected",C_result);
	XL_readNumVector (XL_dDate, C_dDate, " ARM_ERR: domestic date: array of date expected",C_result);
	XL_readNumVector (XL_dSigma, C_dSigma, " ARM_ERR: domestic sigma: array of numeric expected",C_result);
	XL_readNumCell(XL_fMeanRevSpeed,C_fMeanRevSpeed," ARM_ERR: foreign mean reversion speed: numeric expected",C_result);
	XL_readNumVector (XL_fDate, C_fDate, " ARM_ERR: foreign date: array of date expected",C_result);
	XL_readNumVector (XL_fSigma, C_fSigma, " ARM_ERR: foreign sigma: array of numeric expected",C_result);
	XL_readStrCell(XL_dFxCorr,C_dFxCorr," ARM_ERR: domestic forex correlation: object expected",C_result);
	XL_readStrCell(XL_fFxCorr,C_fFxCorr," ARM_ERR: foreign forex correlation: object expected",C_result);
	XL_readStrCell(XL_fxVol,C_fxVol," ARM_ERR: forex volatility: object expected",C_result);
	XL_readNumCell(XL_ratesCorr,C_ratesCorr," ARM_ERR: rates correlation: numeric expected",C_result);
	XL_readNumCell(XL_fxSpotRate,C_fxSpotRate," ARM_ERR: forex spot rates: numeric expected",C_result);
		
	long retCode;
	long objId;
	
	CCString curClass = LOCAL_DFHWSIGVARTREE_CLASS;
	CCString stringId;

	retCode = ARMLOCAL_DFHWSIGVARTREE (C_startDate,
								  C_horizon,
								  (long)C_numSteps,
								  LocalGetNumObjectId (C_dZc),
								  LocalGetNumObjectId (C_fZc),
								  C_dMeanRevSpeed,
								  C_dDate,
								  C_dSigma,
								  C_fMeanRevSpeed,
								  C_fDate,
								  C_fSigma,
								  LocalGetNumObjectId (C_dFxCorr),
								  LocalGetNumObjectId (C_fFxCorr),
								  LocalGetNumObjectId (C_fxVol),
								  C_ratesCorr,
								  C_fxSpotRate,
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
	   PXL_ARM_ERR();
	}

//	ARM_END();
	
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_PXL_DFHWSIGVARTREE" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


__declspec(dllexport) LPXLOPER WINAPI Local_DFHWBASISTREE (LPXLOPER XL_startDate,
														   LPXLOPER XL_horizon,
														   LPXLOPER XL_numSteps,
														   LPXLOPER XL_dZc,
														   LPXLOPER XL_fZc,
														   LPXLOPER XL_dMeanRevSpeed,
														   LPXLOPER XL_dDate,
														   LPXLOPER XL_dSigma,
														   LPXLOPER XL_fDate,
														   LPXLOPER XL_fSigma,
														   LPXLOPER XL_FxCorr,
														   LPXLOPER XL_fxVol,
														   LPXLOPER XL_fxSpotRate,
														   LPXLOPER XL_dNonBasisZeroCurve,
														   LPXLOPER XL_fNonBasisZeroCurve,
														   LPXLOPER XL_dSwaptionBSVol,
														   LPXLOPER XL_dSwaptionBSSmile,
														   LPXLOPER XL_calagetype,
														   LPXLOPER XL_pftype,
														   LPXLOPER XL_data)
{
	ADD_LOG("Local_DFHWBASISTREE ");
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();

	// C variable
	double C_startDate;
	double C_horizon;
	double C_numSteps;
	CCString C_dZc;
	CCString C_fZc;
	double C_dMeanRevSpeed;
	VECTOR<double> C_dDate;
	VECTOR<double> C_dSigma;
	double C_fMeanRevSpeed;
	VECTOR<double> C_fDate;
	VECTOR<double> C_fSigma;

	VECTOR<double> C_data;

	VECTOR<CCString> C_FxCorr;
	CCString C_dFxCorr;
	CCString C_fFxCorr;
	CCString C_fxVol;
	double C_ratesCorr;
	double C_fxSpotRate;

	CCString C_dNonBasisZeroCurve;
	CCString C_fNonBasisZeroCurve;
	CCString C_dSwaptionBSVol;
	CCString C_dSwaptionBSSmile;

	double C_amin;
	double C_amax;
	double C_volmin;
	double C_volmax;

	double C_calagetype;
	double C_pftype;

	// error
	static int error;
	static char* reason = "";

	XL_readNumCell(XL_startDate,C_startDate," ARM_ERR: start date: date expected",C_result);
	XL_readNumCell(XL_horizon,C_horizon," ARM_ERR: horizon date: date expected",C_result);
	XL_readNumCell(XL_numSteps,C_numSteps," ARM_ERR: number of steps: numeric expected",C_result);
	XL_readStrCell(XL_dZc,C_dZc," ARM_ERR: domestic zero-coupon curve id: object expected",C_result);
	XL_readStrCell(XL_fZc,C_fZc," ARM_ERR: foreign zero-coupon curve id: object expected",C_result);
	XL_readNumCell(XL_dMeanRevSpeed,C_dMeanRevSpeed," ARM_ERR: domestic mean reversion speed: numeric expected",C_result);
	XL_readNumVector (XL_dDate, C_dDate, " ARM_ERR: domestic date: array of date expected",C_result);
	XL_readNumVector (XL_dSigma, C_dSigma, " ARM_ERR: domestic sigma: array of numeric expected",C_result);
	XL_readNumVector (XL_fDate, C_fDate, " ARM_ERR: foreign date: array of date expected",C_result);
	XL_readNumVector (XL_fSigma, C_fSigma, " ARM_ERR: foreign sigma: array of numeric expected",C_result);
	XL_readStrVector (XL_FxCorr,C_FxCorr," ARM_ERR: forex correlation: array of object expected",DOUBLE_TYPE,C_result);
	XL_readStrCell(XL_fxVol,C_fxVol," ARM_ERR: forex volatility: object expected",C_result);
	XL_readNumCell(XL_fxSpotRate,C_fxSpotRate," ARM_ERR: forex spot rates: numeric expected",C_result);
	XL_readStrCell(XL_dNonBasisZeroCurve,C_dNonBasisZeroCurve," ARM_ERR: domestic non basis zc: object expected",C_result);
	XL_readStrCell(XL_fNonBasisZeroCurve,C_fNonBasisZeroCurve," ARM_ERR: foreign non basis zc: object expected",C_result);
	XL_readStrCell(XL_dSwaptionBSVol,C_dSwaptionBSVol," ARM_ERR: Swaption BS Vol: object expected",C_result);
	XL_readStrCell(XL_dSwaptionBSSmile,C_dSwaptionBSSmile," ARM_ERR: Swaption BS Smile: object expected",C_result);
	XL_readNumVector (XL_data, C_data, " ARM_ERR: data: array of numeric expected",C_result);
	XL_readNumCell(XL_calagetype,C_calagetype," ARM_ERR: calage type: numeric expected",C_result);
	XL_readNumCell(XL_pftype,C_pftype," ARM_ERR: pf type: numeric expected",C_result);

	C_amin = C_data[0];
	C_amax = C_data[1];
	C_volmin = C_data[2];
	C_volmax = C_data[3];
	C_ratesCorr = C_data[4];
	C_fMeanRevSpeed = C_data[5];
		

	C_dFxCorr = C_FxCorr[0];
	C_fFxCorr = C_FxCorr[1];

	long retCode;
	long objId;
	CCString prevClass;

	CCString curClass = LOCAL_DFHWSIGVARTREE_CLASS;
	CCString stringId = GetLastCurCellEnvValue ();

	if (!stringId)
	{
		retCode = ARMLOCAL_DFHWBASISTREE (C_startDate,
										  C_horizon,
										  (long)C_numSteps,
										  LocalGetNumObjectId (C_dZc),
										  LocalGetNumObjectId (C_fZc),
										  C_dMeanRevSpeed,
										  C_dDate,
										  C_dSigma,
										  C_fMeanRevSpeed,
										  C_fDate,
										  C_fSigma,
										  LocalGetNumObjectId (C_dFxCorr),
										  LocalGetNumObjectId (C_fFxCorr),
										  LocalGetNumObjectId (C_fxVol),
										  C_ratesCorr,
										  C_fxSpotRate,
										  LocalGetNumObjectId (C_dNonBasisZeroCurve),
										  LocalGetNumObjectId (C_fNonBasisZeroCurve),
										  LocalGetNumObjectId (C_dSwaptionBSVol),
										  LocalGetNumObjectId (C_dSwaptionBSSmile),
										  (long) C_calagetype,
										  (long) C_pftype,
										  C_amin,
										  C_amax,
										  C_volmin,
										  C_volmax,
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
			
		if(curClass == prevClass)
		{
			retCode = ARMLOCAL_DFHWBASISTREE (C_startDate,
											  C_horizon,
											  (long)C_numSteps,
											  LocalGetNumObjectId (C_dZc),
											  LocalGetNumObjectId (C_fZc),
											  C_dMeanRevSpeed,
											  C_dDate,
											  C_dSigma,
											  C_fMeanRevSpeed,
											  C_fDate,
											  C_fSigma,
											  LocalGetNumObjectId (C_dFxCorr),
											  LocalGetNumObjectId (C_fFxCorr),
											  LocalGetNumObjectId (C_fxVol),
											  C_ratesCorr,
											  C_fxSpotRate,
											  LocalGetNumObjectId (C_dNonBasisZeroCurve),
											  LocalGetNumObjectId (C_fNonBasisZeroCurve),
											  LocalGetNumObjectId (C_dSwaptionBSVol),
											  LocalGetNumObjectId (C_dSwaptionBSSmile),
											  (long) C_calagetype,
											  (long) C_pftype,
											  C_amin,
											  C_amax,
											  C_volmin,
											  C_volmax,
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
			retCode = ARMLOCAL_DFHWBASISTREE (C_startDate,
											  C_horizon,
											  (long)C_numSteps,
											  LocalGetNumObjectId (C_dZc),
											  LocalGetNumObjectId (C_fZc),
											  C_dMeanRevSpeed,
											  C_dDate,
											  C_dSigma,
											  C_fMeanRevSpeed,
											  C_fDate,
											  C_fSigma,
											  LocalGetNumObjectId (C_dFxCorr),
											  LocalGetNumObjectId (C_fFxCorr),
											  LocalGetNumObjectId (C_fxVol),
											  C_ratesCorr,
											  C_fxSpotRate,
											  LocalGetNumObjectId (C_dNonBasisZeroCurve),
											  LocalGetNumObjectId (C_fNonBasisZeroCurve),
											  LocalGetNumObjectId (C_dSwaptionBSVol),
											  LocalGetNumObjectId (C_dSwaptionBSSmile),
											  (long) C_calagetype,
											  (long) C_pftype,
											  C_amin,
											  C_amax,
											  C_volmin,
											  C_volmax,
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

//	ARM_END();
	
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_DFHWBASISTREE" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}



__declspec(dllexport) LPXLOPER WINAPI Local_PXL_DFHWBASISTREE (LPXLOPER XL_startDate,
															   LPXLOPER XL_horizon,
															   LPXLOPER XL_numSteps,
															   LPXLOPER XL_dZc,
															   LPXLOPER XL_fZc,
															   LPXLOPER XL_dMeanRevSpeed,
															   LPXLOPER XL_dDate,
															   LPXLOPER XL_dSigma,
															   LPXLOPER XL_fDate,
															   LPXLOPER XL_fSigma,
															   LPXLOPER XL_FxCorr,
															   LPXLOPER XL_fxVol,
															   LPXLOPER XL_fxSpotRate,
															   LPXLOPER XL_dNonBasisZeroCurve,
															   LPXLOPER XL_fNonBasisZeroCurve,
															   LPXLOPER XL_dSwaptionBSVol,
															   LPXLOPER XL_dSwaptionBSSmile,
															   LPXLOPER XL_calagetype,
															   LPXLOPER XL_pftype,
															   LPXLOPER XL_data)
{
	ADD_LOG("Local_PXL_DFHWBASISTREE ");
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();

	// C variable
	double C_startDate;
	double C_horizon;
	double C_numSteps;
	CCString C_dZc;
	CCString C_fZc;
	double C_dMeanRevSpeed;
	VECTOR<double> C_dDate;
	VECTOR<double> C_dSigma;
	double C_fMeanRevSpeed;
	VECTOR<double> C_fDate;
	VECTOR<double> C_fSigma;

	VECTOR<double> C_data;

	VECTOR<CCString> C_FxCorr;
	CCString C_dFxCorr;
	CCString C_fFxCorr;
	CCString C_fxVol;
	double C_ratesCorr;
	double C_fxSpotRate;

	CCString C_dNonBasisZeroCurve;
	CCString C_fNonBasisZeroCurve;
	CCString C_dSwaptionBSVol;
	CCString C_dSwaptionBSSmile;

	double C_amin;
	double C_amax;
	double C_volmin;
	double C_volmax;

	double C_calagetype;
	double C_pftype;

	// error
	static int error;
	static char* reason = "";

	XL_readNumCell(XL_startDate,C_startDate," ARM_ERR: start date: date expected",C_result);
	XL_readNumCell(XL_horizon,C_horizon," ARM_ERR: horizon date: date expected",C_result);
	XL_readNumCell(XL_numSteps,C_numSteps," ARM_ERR: number of steps: numeric expected",C_result);
	XL_readStrCell(XL_dZc,C_dZc," ARM_ERR: domestic zero-coupon curve id: object expected",C_result);
	XL_readStrCell(XL_fZc,C_fZc," ARM_ERR: foreign zero-coupon curve id: object expected",C_result);
	XL_readNumCell(XL_dMeanRevSpeed,C_dMeanRevSpeed," ARM_ERR: domestic mean reversion speed: numeric expected",C_result);
	XL_readNumVector (XL_dDate, C_dDate, " ARM_ERR: domestic date: array of date expected",C_result);
	XL_readNumVector (XL_dSigma, C_dSigma, " ARM_ERR: domestic sigma: array of numeric expected",C_result);
	XL_readNumVector (XL_fDate, C_fDate, " ARM_ERR: foreign date: array of date expected",C_result);
	XL_readNumVector (XL_fSigma, C_fSigma, " ARM_ERR: foreign sigma: array of numeric expected",C_result);
	XL_readStrVector (XL_FxCorr,C_FxCorr," ARM_ERR: forex correlation: array of object expected",DOUBLE_TYPE,C_result);
	XL_readStrCell(XL_fxVol,C_fxVol," ARM_ERR: forex volatility: object expected",C_result);
	XL_readNumCell(XL_fxSpotRate,C_fxSpotRate," ARM_ERR: forex spot rates: numeric expected",C_result);
	XL_readStrCell(XL_dNonBasisZeroCurve,C_dNonBasisZeroCurve," ARM_ERR: domestic non basis zc: object expected",C_result);
	XL_readStrCell(XL_fNonBasisZeroCurve,C_fNonBasisZeroCurve," ARM_ERR: foreign non basis zc: object expected",C_result);
	XL_readStrCell(XL_dSwaptionBSVol,C_dSwaptionBSVol," ARM_ERR: Swaption BS Vol: object expected",C_result);
	XL_readStrCell(XL_dSwaptionBSSmile,C_dSwaptionBSSmile," ARM_ERR: Swaption BS Smile: object expected",C_result);
	XL_readNumVector (XL_data, C_data, " ARM_ERR: data: array of numeric expected",C_result);
	XL_readNumCell(XL_calagetype,C_calagetype," ARM_ERR: calage type: numeric expected",C_result);
	XL_readNumCell(XL_pftype,C_pftype," ARM_ERR: pf type: numeric expected",C_result);

	C_amin = C_data[0];
	C_amax = C_data[1];
	C_volmin = C_data[2];
	C_volmax = C_data[3];
	C_ratesCorr = C_data[4];
	C_fMeanRevSpeed = C_data[5];
		

	C_dFxCorr = C_FxCorr[0];
	C_fFxCorr = C_FxCorr[1];
		
	long retCode;
	long objId;
	
	CCString curClass = LOCAL_DFHWSIGVARTREE_CLASS;
	CCString stringId;

			retCode = ARMLOCAL_DFHWBASISTREE (C_startDate,
											  C_horizon,
											  (long)C_numSteps,
											  LocalGetNumObjectId (C_dZc),
											  LocalGetNumObjectId (C_fZc),
											  C_dMeanRevSpeed,
											  C_dDate,
											  C_dSigma,
											  C_fMeanRevSpeed,
											  C_fDate,
											  C_fSigma,
											  LocalGetNumObjectId (C_dFxCorr),
											  LocalGetNumObjectId (C_fFxCorr),
											  LocalGetNumObjectId (C_fxVol),
											  C_ratesCorr,
											  C_fxSpotRate,
											  LocalGetNumObjectId (C_dNonBasisZeroCurve),
											  LocalGetNumObjectId (C_fNonBasisZeroCurve),
											  LocalGetNumObjectId (C_dSwaptionBSVol),
											  LocalGetNumObjectId (C_dSwaptionBSSmile),
											  (long) C_calagetype,
											  (long) C_pftype,
											  C_amin,
											  C_amax,
											  C_volmin,
											  C_volmax,
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
	   PXL_ARM_ERR();
	}

//	ARM_END();
	
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_PXL_DFHWBASISTREE" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}



__declspec(dllexport) LPXLOPER WINAPI Local_GetParameter(LPXLOPER XL_modId,
														 LPXLOPER XL_paramId)
{
	ADD_LOG("Local_GetParameter");
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();

	// C variable
	CCString C_modId;
    double   C_paramId;
	
	// error
	static int error;
	static char* reason = "";
	
    XL_readStrCell(XL_modId, C_modId, " ARM_ERR: model id: object expected",C_result);
    XL_readNumCell(XL_paramId,C_paramId, " ARM_ERR: rank factor: numeric expected",C_result);
    
	long retCode = ARMLOCAL_GetParameter(LocalGetNumObjectId (C_modId), 
                                    (long) C_paramId, C_result);

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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_GetParameter" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}



__declspec(dllexport) LPXLOPER WINAPI Local_FRMTREE(LPXLOPER XL_FrmAnaId, 
													LPXLOPER XL_horizon,
													LPXLOPER XL_fineMonth,
													LPXLOPER XL_corrMatu,
													LPXLOPER XL_corrMatrix,
													LPXLOPER XL_corr2Matu)
{
	ADD_LOG("Local_FRMTREE");
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();

	// C variable
   	CCString C_anaModId;
    double C_horizon;
    CCString CS_fineMonth;
    VECTOR<double> C_corrMatu;
    VECTOR<double> C_corrMatrix;
    VECTOR<double> C_corr2Matu;	

	long fineMonth;

	// error
	static int error;
	static char* reason = "";

	XL_readStrCell(XL_FrmAnaId, C_anaModId, " ARM_ERR: Model id: string expected",C_result);
	XL_readNumCell(XL_horizon,C_horizon," ARM_ERR: Horizon date: date expected",C_result);
	XL_readStrCellWD(XL_fineMonth, CS_fineMonth, "L_DEFAULT", " ARM_ERR: FineMonth: string expected",C_result);

   	if((XL_corrMatu->xltype == xltypeMissing) || (XL_corrMatu->xltype == xltypeNil))
	{
    }
	else
	{
        XL_readNumVector(XL_corrMatu,C_corrMatu, " ARM_ERR: Correlation maturities: array of numeric expected",C_result);
    }

   	if((XL_corrMatrix->xltype == xltypeMissing) || (XL_corrMatrix->xltype == xltypeNil))
	{
    }
	else
	{
        XL_readNumVector(XL_corrMatrix,C_corrMatrix, " ARM_ERR: Correlation matrix: array of numeric expected",C_result);
    }

   	if((XL_corr2Matu->xltype == xltypeMissing) || (XL_corr2Matu->xltype == xltypeNil))
	{
    }
	else
	{
        XL_readNumVector(XL_corr2Matu,C_corr2Matu, " ARM_ERR: Correlation correlated maturities : array of numeric expected",C_result);
    }

	if ((fineMonth = ARM_ConvFineMode(CS_fineMonth, C_result)) == ARM_DEFAULT_ERR)
	{
	   ARM_ARG_ERR();
	   return (LPXLOPER)&XL_result;
	}

	long retCode;
	long objId;
	CCString prevClass;
	
	CCString curClass = LOCAL_FRMTREE_CLASS;
	CCString stringId = GetLastCurCellEnvValue ();

	if (!stringId)
	{
        retCode = ARMLOCAL_FRMTREE(LocalGetNumObjectId (C_anaModId), C_horizon,
                                     fineMonth, C_corrMatu, C_corrMatrix, C_corr2Matu, C_result);

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
			retCode = ARMLOCAL_FRMTREE(LocalGetNumObjectId (C_anaModId), C_horizon,
                                     fineMonth, C_corrMatu, C_corrMatrix, C_corr2Matu, C_result, objId);

			if ( retCode == ARM_OK )
			{
			   LocalSetCurCellEnvValue (curClass, objId); 

			   stringId = LocalMakeObjectId (objId, curClass);
			}
		}
		else
		{
			FreeCurCellContent ();

            retCode = ARMLOCAL_FRMTREE(LocalGetNumObjectId (C_anaModId), C_horizon,
                                     fineMonth, C_corrMatu, C_corrMatrix, C_corr2Matu, C_result);
	
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_FRMTREE" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}



__declspec(dllexport) LPXLOPER WINAPI Local_PXL_FRMTREE(LPXLOPER XL_FrmAnaId, 
														LPXLOPER XL_horizon, 
														LPXLOPER XL_fineMonth,
														LPXLOPER XL_corrMatu,
														LPXLOPER XL_corrMatrix,
														LPXLOPER XL_corr2Matu)                                          
{
	ADD_LOG("Local_PXL_FRMTREE");
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();

	// C variable
   	CCString C_anaModId;
    double C_horizon;
    CCString CS_fineMonth;
    VECTOR<double> C_corrMatu;
    VECTOR<double> C_corrMatrix;
    VECTOR<double> C_corr2Matu;	

	long fineMonth;

	// error
	static int error;
	static char* reason = "";

	XL_readStrCell(XL_FrmAnaId, C_anaModId, " ARM_ERR: Model id: string expected",C_result);
	XL_readNumCell(XL_horizon,C_horizon," ARM_ERR: Horizon date: date expected",C_result);
	XL_readStrCellWD(XL_fineMonth, CS_fineMonth, "L_DEFAULT", " ARM_ERR: FineMonth: string expected",C_result);

   	if((XL_corrMatu->xltype == xltypeMissing) || (XL_corrMatu->xltype == xltypeNil))
	{
    }
	else
	{
        XL_readNumVector(XL_corrMatu,C_corrMatu, " ARM_ERR: Correlation maturities: array of numeric expected",C_result);
    }

   	if((XL_corrMatrix->xltype == xltypeMissing) || (XL_corrMatrix->xltype == xltypeNil))
	{
    }
	else
	{
        XL_readNumVector(XL_corrMatrix,C_corrMatrix, " ARM_ERR: Correlation matrix: array of numeric expected",C_result);
    }

   	if((XL_corr2Matu->xltype == xltypeMissing) || (XL_corr2Matu->xltype == xltypeNil))
	{
    }
	else
	{
        XL_readNumVector(XL_corr2Matu,C_corr2Matu, " ARM_ERR: Correlation correlated maturities : array of numeric expected",C_result);
    }

	if ((fineMonth = ARM_ConvFineMode(CS_fineMonth, C_result)) == ARM_DEFAULT_ERR)
	{
	   ARM_ARG_ERR();
	   return (LPXLOPER)&XL_result;
	}

	long retCode;
	long objId;
	
	CCString curClass = LOCAL_FRMTREE_CLASS;
	CCString stringId;

    retCode = ARMLOCAL_FRMTREE(LocalGetNumObjectId (C_anaModId), C_horizon,
                                 fineMonth, C_corrMatu, C_corrMatrix, C_corr2Matu, C_result);

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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_PXL_FRMTREE" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}



__declspec(dllexport) LPXLOPER WINAPI Local_FRMANA(LPXLOPER XL_zcId,
												   LPXLOPER XL_resetDates,
												   LPXLOPER XL_spotVols,
												   LPXLOPER XL_shapeType,
												   LPXLOPER XL_shapeDecay,
												   LPXLOPER XL_shapeSlope,
												   LPXLOPER XL_shapeAsymptote,
												   LPXLOPER XL_nbFactor,
												   LPXLOPER XL_corrMatu,
												   LPXLOPER XL_corrMatrix,
												   LPXLOPER XL_corr2Matu)
{
	ADD_LOG("Local_FRMANA");
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();

	// C variable
	CCString C_zcId;
    CCString C_portId;
    VECTOR<CCString> C_resetDates;
    VECTOR<double> C_resetDatesDouble;	
    VECTOR<double> C_spotVols;	

    double C_shapeDecay;
    double C_shapeSlope;
    double C_shapeAsymptote;
    double C_nbFactor;
    CCString CS_shapeType; 

   
    VECTOR<double> C_corrMatu;
    VECTOR<double> C_corrMatrix;
    VECTOR<double> C_corr2Matu;	

    double zero = 0;
    double one = 1;
    double row = K_ROW;

	long shapeType;

	// error
	static int error;
	static char* reason = "";

    XL_readStrCell(XL_zcId, C_zcId, " ARM_ERR: Curve id: string expected",C_result);

    XL_readNumVector(XL_resetDates,C_resetDatesDouble, " ARM_ERR: Reset Dates: array of dates expected",C_result);

    XL_readNumVector(XL_spotVols,C_spotVols, " ARM_ERR: Spot Vols: array of numeric expected",C_result);

    XL_readStrCellWD(XL_shapeType, CS_shapeType, "K_ROW", " ARM_ERR: shape type. : string expected",C_result);
  
    XL_readNumCellWD(XL_shapeDecay, C_shapeDecay, zero, " ARM_ERR: shape decay. : numeric expected",C_result);

    XL_readNumCellWD(XL_shapeSlope, C_shapeSlope, zero, " ARM_ERR: shape slope. : numeric expected",C_result);

    XL_readNumCellWD(XL_shapeAsymptote, C_shapeAsymptote, zero, " ARM_ERR: shape Asymptote. : numeric expected",C_result);

    XL_readNumCellWD(XL_nbFactor, C_nbFactor, one, " ARM_ERR: Nb factor. : numeric expected",C_result);

   	if((XL_corrMatu->xltype == xltypeMissing) || (XL_corrMatu->xltype == xltypeNil))
	{
    }
	else
	{
        XL_readNumVector(XL_corrMatu, C_corrMatu, " ARM_ERR: Correlation maturities: array of numeric expected",C_result);
    }

   	if((XL_corrMatrix->xltype == xltypeMissing) || (XL_corrMatrix->xltype == xltypeNil))
	{
    }
	else
	{
        XL_readNumVector(XL_corrMatrix, C_corrMatrix, " ARM_ERR: Correlation matrix: array of numeric expected",C_result);
    }

   	if((XL_corr2Matu->xltype == xltypeMissing) || (XL_corr2Matu->xltype == xltypeNil))
	{
    }
	else
	{
        XL_readNumVector(XL_corr2Matu, C_corr2Matu, " ARM_ERR: Correlation correlated maturities : array of numeric expected",C_result);
    }

	if ((shapeType = ARM_ConvShapeType(CS_shapeType, C_result)) == ARM_DEFAULT_ERR)
	{
	   ARM_ARG_ERR();
	   return (LPXLOPER)&XL_result;
	}


    long real_size = C_resetDatesDouble.size();

	char* sDate = new char[11];
	for (int i = 0; i < real_size; i++)
	{
		Local_XLDATE2ARMDATE(C_resetDatesDouble[i],sDate);
		C_resetDates.push_back((const char*)sDate);
	}
	delete [] sDate;

	long retCode;
	long objId;
	CCString prevClass;
	
	CCString curClass = LOCAL_FRMANA_CLASS;
	CCString stringId = GetLastCurCellEnvValue ();

	if (!stringId)
	{
		retCode = ARMLOCAL_FRMANA(LocalGetNumObjectId(C_zcId), C_resetDates, C_spotVols,
                                   shapeType, C_shapeDecay, C_shapeSlope, C_shapeAsymptote,
                                   (long) C_nbFactor, C_corrMatu, C_corrMatrix, C_corr2Matu, C_result);

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
           retCode = ARMLOCAL_FRMANA(LocalGetNumObjectId(C_zcId), C_resetDates, C_spotVols,
                                   shapeType, C_shapeDecay, C_shapeSlope, C_shapeAsymptote,
                                   (long) C_nbFactor, C_corrMatu, C_corrMatrix, C_corr2Matu, C_result, objId);
		}
		else
		{
			FreeCurCellContent ();

            retCode = ARMLOCAL_FRMANA(LocalGetNumObjectId(C_zcId), C_resetDates, C_spotVols,
                                   shapeType, C_shapeDecay, C_shapeSlope, C_shapeAsymptote,
                                   (long) C_nbFactor, C_corrMatu, C_corrMatrix, C_corr2Matu, C_result);
	
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_FRMANA" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}



__declspec(dllexport) LPXLOPER WINAPI Local_PXL_FRMANA(LPXLOPER XL_zcId,
													   LPXLOPER XL_resetDates,
													   LPXLOPER XL_spotVols,
													   LPXLOPER XL_shapeType,
													   LPXLOPER XL_shapeDecay,
													   LPXLOPER XL_shapeSlope,
													   LPXLOPER XL_shapeAsymptote,
													   LPXLOPER XL_nbFactor,
													   LPXLOPER XL_corrMatu,
													   LPXLOPER XL_corrMatrix,
													   LPXLOPER XL_corr2Matu)
{
	ADD_LOG("Local_PXL_FRMANA");
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();

	// C variable
	CCString C_zcId;
    CCString C_portId;
    VECTOR<CCString> C_resetDates;
    VECTOR<double> C_resetDatesDouble;	
    VECTOR<double> C_spotVols;	

    double C_shapeDecay;
    double C_shapeSlope;
    double C_shapeAsymptote;
    double C_nbFactor;
    CCString CS_shapeType; 

   
    VECTOR<double> C_corrMatu;
    VECTOR<double> C_corrMatrix;
    VECTOR<double> C_corr2Matu;	

    double zero = 0;
    double one = 1;
    double row = K_ROW;

	long shapeType;

	// error
	static int error;
	static char* reason = "";

    XL_readStrCell(XL_zcId, C_zcId, " ARM_ERR: Curve id: string expected",C_result);

    XL_readNumVector(XL_resetDates,C_resetDatesDouble, " ARM_ERR: Reset Dates: array of dates expected",C_result);

    XL_readNumVector(XL_spotVols,C_spotVols, " ARM_ERR: Spot Vols: array of numeric expected",C_result);

    XL_readStrCellWD(XL_shapeType, CS_shapeType, "K_ROW", " ARM_ERR: shape type. : string expected",C_result);
  
    XL_readNumCellWD(XL_shapeDecay, C_shapeDecay, zero, " ARM_ERR: shape decay. : numeric expected",C_result);

    XL_readNumCellWD(XL_shapeSlope, C_shapeSlope, zero, " ARM_ERR: shape slope. : numeric expected",C_result);

    XL_readNumCellWD(XL_shapeAsymptote, C_shapeAsymptote, zero, " ARM_ERR: shape Asymptote. : numeric expected",C_result);

    XL_readNumCellWD(XL_nbFactor, C_nbFactor, one, " ARM_ERR: Nb factor. : numeric expected",C_result);

   	if((XL_corrMatu->xltype == xltypeMissing) || (XL_corrMatu->xltype == xltypeNil))
	{
    }
	else
	{
        XL_readNumVector(XL_corrMatu, C_corrMatu, " ARM_ERR: Correlation maturities: array of numeric expected",C_result);
    }

   	if((XL_corrMatrix->xltype == xltypeMissing) || (XL_corrMatrix->xltype == xltypeNil))
	{
    }
	else
	{
        XL_readNumVector(XL_corrMatrix, C_corrMatrix, " ARM_ERR: Correlation matrix: array of numeric expected",C_result);
    }

   	if((XL_corr2Matu->xltype == xltypeMissing) || (XL_corr2Matu->xltype == xltypeNil))
	{
    }
	else
	{
        XL_readNumVector(XL_corr2Matu, C_corr2Matu, " ARM_ERR: Correlation correlated maturities : array of numeric expected",C_result);
    }

	if ((shapeType = ARM_ConvShapeType(CS_shapeType, C_result)) == ARM_DEFAULT_ERR)
	{
	   ARM_ARG_ERR();
	   return (LPXLOPER)&XL_result;
	}


    long real_size = C_resetDatesDouble.size();

	char* sDate = new char[11];
	for (int i = 0; i < real_size; i++)
	{
		Local_XLDATE2ARMDATE(C_resetDatesDouble[i],sDate);
		C_resetDates.push_back((const char*)sDate);
	}
	delete [] sDate;

	long retCode;
	long objId;
	
	CCString curClass = LOCAL_FRMANA_CLASS;
	CCString stringId;

	retCode = ARMLOCAL_FRMANA(LocalGetNumObjectId(C_zcId), C_resetDates, C_spotVols,
                               shapeType, C_shapeDecay, C_shapeSlope, C_shapeAsymptote,
                               (long) C_nbFactor, C_corrMatu, C_corrMatrix, C_corr2Matu, C_result);

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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_PXL_FRMANA" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


__declspec(dllexport) LPXLOPER WINAPI Local_FRMANA_PORT(LPXLOPER XL_zcId,
														LPXLOPER XL_portId,
														LPXLOPER XL_resetDates,
														LPXLOPER XL_precision,
														LPXLOPER XL_min_paras,
														LPXLOPER XL_max_paras,
														LPXLOPER XL_max_iters,
														LPXLOPER XL_shapeType,
														LPXLOPER XL_shapeDecay,
														LPXLOPER XL_shapeSlope,
														LPXLOPER XL_shapeAsymptote,
														LPXLOPER XL_nbFactor,
														LPXLOPER XL_corrMatu,
														LPXLOPER XL_corrMatrix,
														LPXLOPER XL_corr2Matu) 
                                           
{
	ADD_LOG("Local_FRMANA_PORT");
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();

	// C variable
	CCString C_zcId;
    CCString C_portId;
    VECTOR<CCString> C_resetDates;
    VECTOR<double> C_resetDatesDouble;	
    double C_precision;
    double C_min_paras;
    double C_max_paras;
    double C_max_iters;
    double C_shapeDecay;
    double C_shapeSlope;
    double C_shapeAsymptote;
    double C_nbFactor;
    CCString CS_shapeType; 

    VECTOR<double> C_corrMatu;
    VECTOR<double> C_corrMatrix;
    VECTOR<double> C_corr2Matu;	

    double mFineDefault = K_L_DEFAULT;
    double epsDefault = 1e-8;
    double minusEpsDefault = -1e-8;
    double hugeDefault = 1e8;
    double zero = 0;
    double one = 1;
    double row = K_ROW;

	long shapeType;

	// error
	static int error;
	static char* reason = "";

	XL_readStrCell(XL_zcId, C_zcId, " ARM_ERR: Curve id: string expected",C_result);
	XL_readStrCell(XL_portId, C_portId, " ARM_ERR: Port id: string expected",C_result);
	XL_readNumVector(XL_resetDates, C_resetDatesDouble, " ARM_ERR: Reset Dates: array of dates expected",C_result);
	XL_readNumCellWD(XL_precision,C_precision, epsDefault, " ARM_ERR: precision: numeric expected",C_result);
	XL_readNumCellWD(XL_min_paras, C_min_paras, minusEpsDefault, " ARM_ERR: min paras. : numeric expected",C_result);
	XL_readNumCellWD(XL_max_paras, C_max_paras,hugeDefault," ARM_ERR: max paras. : numeric expected",C_result);
	XL_readNumCellWD(XL_max_iters, C_max_iters,hugeDefault, " ARM_ERR: max iters. : numeric expected",C_result);
	XL_readStrCellWD(XL_shapeType, CS_shapeType, "ROW", " ARM_ERR: shape type. : string expected",C_result);
	XL_readNumCellWD(XL_shapeDecay, C_shapeDecay, zero, " ARM_ERR: shape decay. : numeric expected",C_result);
	XL_readNumCellWD(XL_shapeSlope, C_shapeSlope, zero, " ARM_ERR: shape slope. : numeric expected",C_result);
	XL_readNumCellWD(XL_shapeAsymptote, C_shapeAsymptote, zero, " ARM_ERR: shape Asymptote. : numeric expected",C_result);
	XL_readNumCellWD(XL_nbFactor, C_nbFactor, one, " ARM_ERR: Nb factor. : numeric expected",C_result);

   	if((XL_corrMatu->xltype == xltypeMissing) || (XL_corrMatu->xltype == xltypeNil))
	{
    }
	else
	{
        XL_readNumVector(XL_corrMatu, C_corrMatu, " ARM_ERR: Correlation maturities: array of numeric expected",C_result);
    }

   	if((XL_corrMatrix->xltype == xltypeMissing) || (XL_corrMatrix->xltype == xltypeNil))
	{
    }
	else
	{
        XL_readNumVector(XL_corrMatrix, C_corrMatrix, " ARM_ERR: Correlation matrix: array of numeric expected",C_result);
    }

   	if((XL_corr2Matu->xltype == xltypeMissing) || (XL_corr2Matu->xltype == xltypeNil))
	{
    }
	else
	{
        XL_readNumVector(XL_corr2Matu, C_corr2Matu, " ARM_ERR: Correlation correlated maturities : array of numeric expected",C_result);
    }

	if ((shapeType = ARM_ConvShapeType(CS_shapeType, C_result)) == ARM_DEFAULT_ERR)
	{
	   ARM_ARG_ERR();
	   return (LPXLOPER)&XL_result;
	}

    long real_size = C_resetDatesDouble.size();

	char* sDate = new char[11];
	for (int i = 0; i < real_size; i++)
	{
		Local_XLDATE2ARMDATE(C_resetDatesDouble[i],sDate);
		C_resetDates.push_back((const char*)sDate);
	}

	delete [] sDate;
  
	long retCode;
	long objId;
	CCString prevClass;
	
	CCString curClass = LOCAL_FRMANA_CLASS;
	CCString stringId = GetLastCurCellEnvValue ();

	if (!stringId)
	{

        retCode = ARMLOCAL_FRMANA_PORT(LocalGetNumObjectId(C_zcId), LocalGetNumObjectId(C_portId), C_resetDates, 
                                   C_precision, C_min_paras, C_max_paras, (long) C_max_iters, shapeType,
                                   C_shapeDecay, C_shapeSlope, C_shapeAsymptote,
                                   (long) C_nbFactor, C_corrMatu, C_corrMatrix, C_corr2Matu, C_result);

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
			retCode = ARMLOCAL_FRMANA_PORT(LocalGetNumObjectId(C_zcId), LocalGetNumObjectId(C_portId), C_resetDates, 
								   C_precision, C_min_paras, C_max_paras, (long) C_max_iters, shapeType,
								   C_shapeDecay, C_shapeSlope, C_shapeAsymptote,
								   (long) C_nbFactor, C_corrMatu, C_corrMatrix, C_corr2Matu, C_result, objId);

			if ( retCode == ARM_OK )
			{
			   LocalSetCurCellEnvValue (curClass, objId); 

			   stringId = LocalMakeObjectId (objId, curClass);
			}
		}
		else
		{
			FreeCurCellContent ();

            retCode = ARMLOCAL_FRMANA_PORT(LocalGetNumObjectId(C_zcId), LocalGetNumObjectId(C_portId), C_resetDates, 
                                   C_precision, C_min_paras, C_max_paras, (long) C_max_iters, shapeType,
                                   C_shapeDecay, C_shapeSlope, C_shapeAsymptote,
                                   (long) C_nbFactor, C_corrMatu, C_corrMatrix, C_corr2Matu, C_result);
	
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_FRMANA_PORT" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}



__declspec(dllexport) LPXLOPER WINAPI Local_PXL_FRMANA_PORT(LPXLOPER XL_zcId,
															LPXLOPER XL_portId,
															LPXLOPER XL_resetDates,
															LPXLOPER XL_precision,
															LPXLOPER XL_min_paras,
															LPXLOPER XL_max_paras,
															LPXLOPER XL_max_iters,
															LPXLOPER XL_shapeType,
															LPXLOPER XL_shapeDecay,
															LPXLOPER XL_shapeSlope,
															LPXLOPER XL_shapeAsymptote,
															LPXLOPER XL_nbFactor,
															LPXLOPER XL_corrMatu,
															LPXLOPER XL_corrMatrix,
															LPXLOPER XL_corr2Matu)
{
	ADD_LOG("Local_PXL_FRMANA_PORT");
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();

	// C variable
	CCString C_zcId;
    CCString C_portId;
    VECTOR<CCString> C_resetDates;
    VECTOR<double> C_resetDatesDouble;
    double C_precision;
    double C_min_paras;
    double C_max_paras;
    double C_max_iters;
    double C_shapeDecay;
    double C_shapeSlope;
    double C_shapeAsymptote;
    double C_nbFactor;

    CCString CS_shapeType;
	long shapeType;
   
    VECTOR<double> C_corrMatu;
    VECTOR<double> C_corrMatrix;
    VECTOR<double> C_corr2Matu;	

    double epsDefault = 1e-8;
    double minusEpsDefault = -1e-8;
    double hugeDefault = 1e8;
    double zero = 0;
    double one = 1;

	// error
	static int error;
	static char* reason = "";

    XL_readStrCell(XL_zcId, C_zcId," ARM_ERR: Curve id: string expected",C_result);
    XL_readStrCell(XL_portId, C_portId," ARM_ERR: Port id: string expected",C_result);
    XL_readNumVector(XL_resetDates,C_resetDatesDouble," ARM_ERR: Reset Dates: array of dates expected",C_result);
	XL_readNumCellWD(XL_precision,C_precision, epsDefault, " ARM_ERR: precision: numeric expected",C_result);
    XL_readNumCellWD(XL_min_paras, C_min_paras, minusEpsDefault," ARM_ERR: min paras. : numeric expected",C_result);
    XL_readNumCellWD(XL_max_paras, C_max_paras, hugeDefault,  " ARM_ERR: max paras. : numeric expected",C_result);
    XL_readNumCellWD(XL_max_iters, C_max_iters, hugeDefault, " ARM_ERR: max iters. : numeric expected",C_result);
    XL_readStrCellWD(XL_shapeType, CS_shapeType, "ROW",  " ARM_ERR: shape type. : string expected",C_result);
    XL_readNumCellWD(XL_shapeDecay, C_shapeDecay, zero, " ARM_ERR: shape decay. : numeric expected",C_result);
    XL_readNumCellWD(XL_shapeSlope, C_shapeSlope, zero, " ARM_ERR: shape slope. : numeric expected",C_result);
    XL_readNumCellWD(XL_shapeAsymptote, C_shapeAsymptote, zero, " ARM_ERR: shape Asymptote. : numeric expected",C_result);
    XL_readNumCellWD(XL_nbFactor, C_nbFactor,one, " ARM_ERR: Nb factor. : numeric expected",C_result);

   	if((XL_corrMatu->xltype == xltypeMissing) || (XL_corrMatu->xltype == xltypeNil))
	{
    }
	else
	{
        XL_readNumVector(XL_corrMatu,C_corrMatu," ARM_ERR: Correlation maturities: array of numeric expected",C_result);
    }

   	if((XL_corrMatrix->xltype == xltypeMissing) || (XL_corrMatrix->xltype == xltypeNil))
	{
    }
	else
	{
        XL_readNumVector(XL_corrMatrix,C_corrMatrix," ARM_ERR: Correlation matrix: array of numeric expected",C_result);
    }

   	if((XL_corr2Matu->xltype == xltypeMissing) || (XL_corr2Matu->xltype == xltypeNil))
	{
    }
	else
	{
        XL_readNumVector(XL_corr2Matu,C_corr2Matu," ARM_ERR: Correlation correlated maturities : array of numeric expected",C_result);
    }

	if ((shapeType = ARM_ConvShapeType(CS_shapeType, C_result)) == ARM_DEFAULT_ERR)
	{
	   ARM_ARG_ERR();
	   return (LPXLOPER)&XL_result;
	}

    long real_size = C_resetDatesDouble.size();

	char* sDate = new char[11];
	for (int i = 0; i < real_size; i++)
	{
		Local_XLDATE2ARMDATE(C_resetDatesDouble[i],sDate);
		C_resetDates.push_back((const char*)sDate);
	}

	delete [] sDate;

	long retCode;
	long objId;

	CCString curClass = LOCAL_FRMANA_CLASS;
	CCString stringId;

    retCode = ARMLOCAL_FRMANA_PORT(LocalGetNumObjectId(C_zcId), LocalGetNumObjectId(C_portId), C_resetDates, 
                                   C_precision, C_min_paras, C_max_paras, (long) C_max_iters,  shapeType,
                                   C_shapeDecay, C_shapeSlope, C_shapeAsymptote,
                                   (long) C_nbFactor, C_corrMatu, C_corrMatrix, C_corr2Matu, C_result);

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
		PXL_ARM_ERR();
	}

//	ARM_END();
	
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_PXL_FRMANA_PORT" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


__declspec(dllexport) LPXLOPER WINAPI Local_FRMLSMC(LPXLOPER XL_FrmAnaId,
													LPXLOPER XL_horizon,
													LPXLOPER XL_nbTraj,
													LPXLOPER XL_fineMonth,
													LPXLOPER XL_mcMethod,
													LPXLOPER XL_corrMatu,
													LPXLOPER XL_corrMatrix,
													LPXLOPER XL_corr2Matu)
{
	ADD_LOG("Local_FRMLSMC");
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();

	// C variable
   	CCString C_anaModId;
    double C_horizon;
    double C_fineMonth;
    double C_nbTraj;
    CCString CS_mcMethod;
    VECTOR<double> C_corrMatu;
    VECTOR<double> C_corrMatrix;
    VECTOR<double> C_corr2Matu;	 
	
	long mcMethod;

	// error
	static int error;
	static char* reason = "";

    double nbTrajDef = 1000;
    double trois = 3;

    XL_readStrCell(XL_FrmAnaId, C_anaModId, " ARM_ERR: Model id: string expected",C_result);
	XL_readNumCell(XL_horizon, C_horizon, " ARM_ERR: Horizon date: date expected",C_result);
    XL_readNumCellWD(XL_nbTraj, C_nbTraj, nbTrajDef, " ARM_ERR: Number of trajectories . : numeric expected",C_result);
    XL_readNumCellWD(XL_fineMonth, C_fineMonth, trois, " ARM_ERR: fine month. : numeric expected",C_result);
    XL_readStrCellWD(XL_mcMethod, CS_mcMethod, "FAURE", " ARM_ERR: MC method : string expected",C_result);

   	if((XL_corrMatu->xltype == xltypeMissing) || (XL_corrMatu->xltype == xltypeNil))
	{
    }
	else
	{
        XL_readNumVector(XL_corrMatu, C_corrMatu, " ARM_ERR: Correlation maturities: array of numeric expected",C_result);
    }

   	if((XL_corrMatrix->xltype == xltypeMissing) || (XL_corrMatrix->xltype == xltypeNil))
	{
    }
	else
	{
        XL_readNumVector(XL_corrMatrix, C_corrMatrix, " ARM_ERR: Correlation matrix: array of numeric expected",C_result);
    }

   	if((XL_corr2Matu->xltype == xltypeMissing) || (XL_corr2Matu->xltype == xltypeNil))
	{
    }
	else
	{
        XL_readNumVector(XL_corr2Matu, C_corr2Matu, " ARM_ERR: Correlation correlated maturities : array of numeric expected",C_result);
    }

	if ((mcMethod = ARM_ConvMcMethod(CS_mcMethod, C_result)) == ARM_DEFAULT_ERR)
	{
	   ARM_ARG_ERR();
	   return (LPXLOPER)&XL_result;
	}

	long retCode;
	long objId;
	CCString prevClass;
	
	CCString curClass = LOCAL_FRMLSMONTECARLO_CLASS;
	CCString stringId = GetLastCurCellEnvValue ();

	if (!stringId)
	{
        retCode = ARMLOCAL_FRMLSMC(LocalGetNumObjectId (C_anaModId), C_horizon, (long) C_nbTraj, 
                                     C_fineMonth, mcMethod, C_corrMatu, C_corrMatrix, C_corr2Matu, C_result);

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
           retCode = ARMLOCAL_FRMLSMC(LocalGetNumObjectId (C_anaModId), C_horizon, C_nbTraj, 
                                     C_fineMonth, mcMethod, C_corrMatu, C_corrMatrix, C_corr2Matu, C_result, objId);

			if ( retCode == ARM_OK )
			{
			   LocalSetCurCellEnvValue (curClass, objId); 

			   stringId = LocalMakeObjectId (objId, curClass);
			}
		}
		else
		{
			FreeCurCellContent ();

            retCode = ARMLOCAL_FRMLSMC(LocalGetNumObjectId (C_anaModId), C_horizon, C_nbTraj, 
                                     C_fineMonth, mcMethod, C_corrMatu, C_corrMatrix, C_corr2Matu, C_result);
	
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_FRMLSMC" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}



__declspec(dllexport) LPXLOPER WINAPI Local_PXL_FRMLSMC(LPXLOPER XL_FrmAnaId,
														LPXLOPER XL_horizon,
														LPXLOPER XL_nbTraj,
														LPXLOPER XL_fineMonth,
														LPXLOPER XL_mcMethod,
														LPXLOPER XL_corrMatu,
														LPXLOPER XL_corrMatrix,
														LPXLOPER XL_corr2Matu)
{
	ADD_LOG("Local_PXL_FRMLSMC");
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();

	// C variable
   	CCString C_anaModId;
    CCString CS_fineMonth;
    CCString CS_mcMethod;
	long mcMethod;

    double C_nbTraj;
    double C_horizon;
    double C_fineMonth;

    double mFineDefault = K_L_DEFAULT;
    VECTOR<double> C_corrMatu;
    VECTOR<double> C_corrMatrix;
    VECTOR<double> C_corr2Matu;	
	
	// error
	static int error;
	static char* reason = "";

    double nbTrajDef = 1000;

    XL_readStrCell(XL_FrmAnaId,C_anaModId," ARM_ERR: Model id: string expected",C_result);
	XL_readNumCell(XL_horizon,C_horizon," ARM_ERR: Horizon date: date expected",C_result);
    XL_readNumCellWD(XL_nbTraj, C_nbTraj, nbTrajDef, " ARM_ERR: Number of trajectories . : numeric expected",C_result);
    XL_readNumCellWD(XL_fineMonth, C_fineMonth, mFineDefault, " ARM_ERR: fine month. : numeric expected",C_result);
    XL_readStrCellWD(XL_mcMethod, CS_mcMethod, "FAURE"," ARM_ERR: MC method : string expected",C_result);

   	if((XL_corrMatu->xltype == xltypeMissing) || (XL_corrMatu->xltype == xltypeNil))
	{
    }
	else
	{
        XL_readNumVector(XL_corrMatu,C_corrMatu,
          " ARM_ERR: Correlation maturities: array of numeric expected",C_result);
    }

   	if((XL_corrMatrix->xltype == xltypeMissing) || (XL_corrMatrix->xltype == xltypeNil))
	{
    }
	else
	{
        XL_readNumVector(XL_corrMatrix,C_corrMatrix,
          " ARM_ERR: Correlation matrix: array of numeric expected",C_result);
    }

   	if((XL_corr2Matu->xltype == xltypeMissing) || (XL_corr2Matu->xltype == xltypeNil))
	{
    }
	else
	{
        XL_readNumVector(XL_corr2Matu,C_corr2Matu,
          " ARM_ERR: Correlation correlated maturities : array of numeric expected",C_result);
    }

	if ((mcMethod = ARM_ConvMcMethod(CS_mcMethod, C_result)) == ARM_DEFAULT_ERR)
	{
	   ARM_ARG_ERR();
	   return (LPXLOPER)&XL_result;
	}

	long retCode;
	long objId;
	
	CCString curClass = LOCAL_FRMLSMONTECARLO_CLASS;
	CCString stringId;

    retCode = ARMLOCAL_FRMLSMC(LocalGetNumObjectId (C_anaModId), C_horizon, (long) C_nbTraj, 
                                     C_fineMonth, mcMethod, C_corrMatu, C_corrMatrix, C_corr2Matu, C_result);

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
	   PXL_ARM_ERR();
	}

//	ARM_END();
	
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_PXL_FRMLSMC" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}



__declspec(dllexport) LPXLOPER WINAPI Local_FRMLSMC_AUTO(LPXLOPER XL_zcId,
														 LPXLOPER XL_volId,
														 LPXLOPER XL_smileId,
														 LPXLOPER XL_autoMode,
														 LPXLOPER XL_horizon,
														 LPXLOPER XL_nbTraj,
														 LPXLOPER XL_fineMonth,
														 LPXLOPER XL_mcMethod,
														 LPXLOPER XL_noControl,
														 LPXLOPER XL_shapeDecay,
														 LPXLOPER XL_shapeSlope,
														 LPXLOPER XL_shapeAsymptote,
														 LPXLOPER XL_nbFactor,
														 LPXLOPER XL_corrMatu,
														 LPXLOPER XL_corrMatrix,
														 LPXLOPER XL_corr2Matu) 
                                           
{
	ADD_LOG("Local_FRMLSMC_AUTO");
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();

	// C variable
	CCString C_zcId;
    CCString C_volId;
    CCString C_smileId;
    CCString CS_autoMode;
    CCString CS_fineMonth;
    CCString CS_mcMethod;
    CCString CS_noControl;
    double C_horizon;
    double C_nbTraj;
    double C_shapeDecay;
    double C_shapeSlope;
    double C_shapeAsymptote;
    double C_nbFactor;
    
   
    VECTOR<double> C_corrMatu;
    VECTOR<double> C_corrMatrix;
    VECTOR<double> C_corr2Matu;	

    double mFineDefault = K_L_DEFAULT;
    double epsDefault = 1e-8;
    double huegeDefault = 1e8;
    double zero = 0;
    double one = 1;
    double thousand = 1000;

	long intAutoMode;
	long fineMonth; 
	long mcMethod;
	long noControl;

	// error
	static int error;
	static char* reason = "";

    XL_readStrCell(XL_zcId, C_zcId, " ARM_ERR: Curve id: string expected",C_result);

    XL_readStrCell(XL_volId, C_volId, " ARM_ERR: Vol id: string expected",C_result);

    XL_readStrCell(XL_smileId, C_smileId, " ARM_ERR: Smile id: string expected",C_result);

    XL_readStrCell(XL_autoMode, CS_autoMode, " ARM_ERR: AutoMode: string expected",C_result);

	XL_readNumCell(XL_horizon, C_horizon," ARM_ERR: Horizon date: date expected",C_result);

    XL_readNumCellWD(XL_nbTraj, C_nbTraj, thousand," ARM_ERR: Number of trajectories : numeric expected",C_result);

    XL_readStrCellWD(XL_fineMonth, CS_fineMonth, "L_DEFAULT", " ARM_ERR: fine month. : string expected",C_result);

    XL_readStrCellWD(XL_mcMethod, CS_mcMethod, "FAURE", " ARM_ERR: MC method : string expected",C_result);

    XL_readStrCellWD(XL_noControl, CS_noControl, "IN", " ARM_ERR: Control variate: string expected : IN or OUT",C_result);

    XL_readNumCellWD(XL_shapeDecay, C_shapeDecay, zero, " ARM_ERR: shape decay. : numeric expected",C_result);

    XL_readNumCellWD(XL_shapeSlope, C_shapeSlope, zero, " ARM_ERR: shape slope. : numeric expected",C_result);

    XL_readNumCellWD(XL_shapeAsymptote, C_shapeAsymptote, zero, " ARM_ERR: shape Asymptote. : numeric expected",C_result);

    XL_readNumCellWD(XL_nbFactor, C_nbFactor, one, " ARM_ERR: Nb factor. : numeric expected",C_result);

   	if((XL_corrMatu->xltype == xltypeMissing) || (XL_corrMatu->xltype == xltypeNil))
	{
    }
	else
	{
        XL_readNumVector(XL_corrMatu, C_corrMatu, " ARM_ERR: Correlation maturities: array of numeric expected",C_result);
    }

   	if((XL_corrMatrix->xltype == xltypeMissing) || (XL_corrMatrix->xltype == xltypeNil))
	{
    }
	else
	{
        XL_readNumVector(XL_corrMatrix, C_corrMatrix, " ARM_ERR: Correlation matrix: array of numeric expected",C_result);
    }

   	if((XL_corr2Matu->xltype == xltypeMissing) || (XL_corr2Matu->xltype == xltypeNil))
	{
    }
	else
	{
        XL_readNumVector(XL_corr2Matu, C_corr2Matu, " ARM_ERR: Correlation correlated maturities : array of numeric expected",C_result);
    }

	if ((intAutoMode = ARM_ConvAutoMode(CS_autoMode, C_result)) == ARM_DEFAULT_ERR)
	{
	   ARM_ARG_ERR();
	   return (LPXLOPER)&XL_result;
	}

	if ((fineMonth = ARM_ConvFineMode(CS_fineMonth, C_result)) == ARM_DEFAULT_ERR)
	{
	   ARM_ARG_ERR();
	   return (LPXLOPER)&XL_result;
	}

	if ((mcMethod = ARM_ConvMcMethod(CS_mcMethod, C_result)) == ARM_DEFAULT_ERR)
	{
	   ARM_ARG_ERR();
	   return (LPXLOPER)&XL_result;
	}

	if ((noControl = ARM_ConvInOut(CS_noControl, C_result)) == ARM_DEFAULT_ERR)
	{
	   ARM_ARG_ERR();
	   return (LPXLOPER)&XL_result;
	}

	long retCode;
	long objId;
	CCString prevClass;

	CCString curClass = LOCAL_FRMLSMONTECARLO_CLASS;
	CCString stringId = GetLastCurCellEnvValue ();

	if (!stringId)
	{
        retCode = ARMLOCAL_FRMLSMC_AUTO(LocalGetNumObjectId(C_zcId), LocalGetNumObjectId(C_volId), LocalGetNumObjectId(C_smileId), 
                                   intAutoMode,  C_horizon, fineMonth, C_nbTraj, mcMethod, noControl,
                                   C_shapeDecay, C_shapeSlope, C_shapeAsymptote,
                                   (long) C_nbFactor, C_corrMatu, C_corrMatrix, C_corr2Matu, C_result);

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
           retCode = ARMLOCAL_FRMLSMC_AUTO(LocalGetNumObjectId(C_zcId), LocalGetNumObjectId(C_volId), LocalGetNumObjectId(C_smileId), 
                                   intAutoMode,  C_horizon, fineMonth, C_nbTraj, mcMethod, noControl,
                                   C_shapeDecay, C_shapeSlope, C_shapeAsymptote,
                                   (long) C_nbFactor, C_corrMatu, C_corrMatrix, C_corr2Matu, C_result, objId);

			if ( retCode == ARM_OK )
			{
			   LocalSetCurCellEnvValue (curClass, objId); 

			   stringId = LocalMakeObjectId (objId, curClass);
			}
		}
		else
		{
			FreeCurCellContent ();

            retCode = ARMLOCAL_FRMLSMC_AUTO(LocalGetNumObjectId(C_zcId), LocalGetNumObjectId(C_volId), LocalGetNumObjectId(C_smileId), 
                                   intAutoMode,  C_horizon, fineMonth, C_nbTraj, mcMethod, noControl,
                                   C_shapeDecay, C_shapeSlope, C_shapeAsymptote,
                                   (long) C_nbFactor, C_corrMatu, C_corrMatrix, C_corr2Matu, C_result, objId);
	
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_FRMLSMC_AUTO" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


__declspec(dllexport) LPXLOPER WINAPI Local_PXL_FRMLSMC_AUTO(LPXLOPER XL_zcId,
															 LPXLOPER XL_volId,
															 LPXLOPER XL_smileId,
															 LPXLOPER XL_autoMode,
															 LPXLOPER XL_horizon,
															 LPXLOPER XL_nbTraj,
															 LPXLOPER XL_fineMonth,
															 LPXLOPER XL_mcMethod,
															 LPXLOPER XL_noControl,
															 LPXLOPER XL_shapeDecay,
															 LPXLOPER XL_shapeSlope,
															 LPXLOPER XL_shapeAsymptote,
															 LPXLOPER XL_nbFactor,
															 LPXLOPER XL_corrMatu,
															 LPXLOPER XL_corrMatrix,
															 LPXLOPER XL_corr2Matu)
{
	ADD_LOG("Local_PXL_FRMLSMC_AUTO");
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;
	
	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	// C variable
	CCString C_zcId;
	CCString C_volId;
	CCString C_smileId;
	CCString CS_autoMode;
	CCString CS_fineMonth;
	CCString CS_mcMethod;
	CCString CS_noControl;
	double C_horizon;
	double C_nbTraj;
	double C_shapeDecay;
	double C_shapeSlope;
	double C_shapeAsymptote;
	double C_nbFactor;


	VECTOR<double> C_corrMatu;
	VECTOR<double> C_corrMatrix;
	VECTOR<double> C_corr2Matu;	

	double mFineDefault = K_L_DEFAULT;
	double epsDefault = 1e-8;
	double huegeDefault = 1e8;
	double zero = 0;
	double one = 1;
	double thousand = 1000;

	long intAutoMode;
	long fineMonth; 
	long mcMethod;
	long noControl;

	// error
	static int error;
	static char* reason = "";

    XL_readStrCell(XL_zcId, C_zcId, " ARM_ERR: Curve id: string expected",C_result);

    XL_readStrCell(XL_volId, C_volId, " ARM_ERR: Vol id: string expected",C_result);

    XL_readStrCell(XL_smileId, C_smileId, " ARM_ERR: Smile id: string expected",C_result);

    XL_readStrCell(XL_autoMode, CS_autoMode, " ARM_ERR: AutoMode: string expected",C_result);

	XL_readNumCell(XL_horizon, C_horizon," ARM_ERR: Horizon date: date expected",C_result);

    XL_readNumCellWD(XL_nbTraj, C_nbTraj, thousand," ARM_ERR: Number of trajectories : numeric expected",C_result);

    XL_readStrCellWD(XL_fineMonth, CS_fineMonth, "L_DEFAULT", " ARM_ERR: fine month. : string expected",C_result);

    XL_readStrCellWD(XL_mcMethod, CS_mcMethod, "FAURE", " ARM_ERR: MC method : string expected",C_result);

    XL_readStrCellWD(XL_noControl, CS_noControl, "IN", " ARM_ERR: Control variate: string expected : IN or OUT",C_result);

    XL_readNumCellWD(XL_shapeDecay, C_shapeDecay, zero, " ARM_ERR: shape decay. : numeric expected",C_result);

    XL_readNumCellWD(XL_shapeSlope, C_shapeSlope, zero, " ARM_ERR: shape slope. : numeric expected",C_result);

    XL_readNumCellWD(XL_shapeAsymptote, C_shapeAsymptote, zero, " ARM_ERR: shape Asymptote. : numeric expected",C_result);

    XL_readNumCellWD(XL_nbFactor, C_nbFactor, one, " ARM_ERR: Nb factor. : numeric expected",C_result);

   	if((XL_corrMatu->xltype == xltypeMissing) || (XL_corrMatu->xltype == xltypeNil))
	{
    }
	else
	{
        XL_readNumVector(XL_corrMatu, C_corrMatu, " ARM_ERR: Correlation maturities: array of numeric expected",C_result);
    }

   	if((XL_corrMatrix->xltype == xltypeMissing) || (XL_corrMatrix->xltype == xltypeNil))
	{
    }
	else
	{
        XL_readNumVector(XL_corrMatrix, C_corrMatrix, " ARM_ERR: Correlation matrix: array of numeric expected",C_result);
    }

   	if((XL_corr2Matu->xltype == xltypeMissing) || (XL_corr2Matu->xltype == xltypeNil))
	{
    }
	else
	{
        XL_readNumVector(XL_corr2Matu, C_corr2Matu, " ARM_ERR: Correlation correlated maturities : array of numeric expected",C_result);
    }

	if ((intAutoMode = ARM_ConvAutoMode(CS_autoMode, C_result)) == ARM_DEFAULT_ERR)
	{
	   ARM_ARG_ERR();
	   return (LPXLOPER)&XL_result;
	}

	if ((fineMonth = ARM_ConvFineMode(CS_fineMonth, C_result)) == ARM_DEFAULT_ERR)
	{
	   ARM_ARG_ERR();
	   return (LPXLOPER)&XL_result;
	}

	if ((mcMethod = ARM_ConvMcMethod(CS_mcMethod, C_result)) == ARM_DEFAULT_ERR)
	{
	   ARM_ARG_ERR();
	   return (LPXLOPER)&XL_result;
	}

	if ((noControl = ARM_ConvInOut(CS_noControl, C_result)) == ARM_DEFAULT_ERR)
	{
	   ARM_ARG_ERR();
	   return (LPXLOPER)&XL_result;
	}

	long retCode;
	long objId;

	CCString curClass = LOCAL_FRMLSMONTECARLO_CLASS;
	CCString stringId;

    retCode = ARMLOCAL_FRMLSMC_AUTO(LocalGetNumObjectId(C_zcId), LocalGetNumObjectId(C_volId), LocalGetNumObjectId(C_smileId), 
                               intAutoMode,  C_horizon, fineMonth, C_nbTraj, mcMethod, noControl,
                               C_shapeDecay, C_shapeSlope, C_shapeAsymptote,
                               (long) C_nbFactor, C_corrMatu, C_corrMatrix, C_corr2Matu, C_result);

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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_PXL_FRMLSMC_AUTO" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}




__declspec(dllexport) LPXLOPER WINAPI Local_FRMLSMC_AUTO_2(LPXLOPER XL_zcId,
														   LPXLOPER XL_swopVolId,
														   LPXLOPER XL_swopSmileId,
														   LPXLOPER XL_irgVolId,
														   LPXLOPER XL_irgSmileId,
														   LPXLOPER XL_autoMode,
														   LPXLOPER XL_horizon,
														   LPXLOPER XL_nbTraj,
														   LPXLOPER XL_fineMonth,
														   LPXLOPER XL_mcMethod,
														   LPXLOPER XL_noControl,
														   LPXLOPER XL_shapeDecay,
														   LPXLOPER XL_shapeSlope,
														   LPXLOPER XL_shapeAsymptote,
														   LPXLOPER XL_nbFactor,
														   LPXLOPER XL_corrMatu,
														   LPXLOPER XL_corrMatrix,
														   LPXLOPER XL_corr2Matu) 
                                           
{
	ADD_LOG("Local_FRMLSMC_AUTO_2");
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();

	// C variable
	CCString C_zcId;
    CCString C_swopVolId;
    CCString C_swopSmileId;
    CCString C_irgVolId;
    CCString C_irgSmileId;
    CCString CS_autoMode;
	long C_autoMode;
    CCString CS_fineMonth;
	long C_fineMonth;
    CCString CS_mcMethod;
	long C_mcMethod;
    CCString CS_noControl;
	long C_noControl;
    double C_horizon;
    double C_nbTraj;
    double C_shapeDecay;
    double C_shapeSlope;
    double C_shapeAsymptote;
    double C_nbFactor;

    VECTOR<double> C_corrMatu;
    VECTOR<double> C_corrMatrix;
    VECTOR<double> C_corr2Matu;	

    double mFineDefault = K_L_DEFAULT;
    double epsDefault = 1e-8;
    double huegeDefault = 1e8;
    double zero = 0;
    double one = 1;
    double thousand = 1000;

	// error
	static int error;
	static char* reason = "";

    XL_readStrCell(XL_zcId, C_zcId, " ARM_ERR: Curve id: string expected",C_result);

    XL_readStrCell(XL_swopVolId, C_swopVolId, " ARM_ERR: Swaption Vol id: string expected",C_result);

    XL_readStrCell(XL_swopSmileId, C_swopSmileId, " ARM_ERR: Swaption Smile id: string expected",C_result);

    XL_readStrCell(XL_irgVolId, C_irgVolId, " ARM_ERR: Cap Vol id: string expected",C_result);

    XL_readStrCell(XL_irgSmileId, C_irgSmileId, " ARM_ERR: Cap Smile id: string expected",C_result);

    XL_readStrCell(XL_autoMode, CS_autoMode, " ARM_ERR: AutoMode: string expected",C_result);

	XL_readNumCell(XL_horizon,C_horizon," ARM_ERR: Horizon date: date expected",C_result);

    XL_readNumCellWD(XL_nbTraj,C_nbTraj, thousand," ARM_ERR: Number of trajectories : numeric expected",C_result);

    XL_readStrCellWD(XL_fineMonth, CS_fineMonth, "L_DEFAULT", " ARM_ERR: fine month. : string expected",C_result);

    XL_readStrCellWD(XL_mcMethod, CS_mcMethod, "FAURE", " ARM_ERR: MC method : string expected",C_result);

    XL_readStrCellWD(XL_noControl, CS_noControl, "IN", " ARM_ERR: Control variate: string expected : IN or OUT",C_result);

    XL_readNumCellWD(XL_shapeDecay, C_shapeDecay, zero, " ARM_ERR: shape decay. : numeric expected",C_result);

    XL_readNumCellWD(XL_shapeSlope, C_shapeSlope, zero, " ARM_ERR: shape slope. : numeric expected",C_result);

    XL_readNumCellWD(XL_shapeAsymptote, C_shapeAsymptote, zero, " ARM_ERR: shape Asymptote. : numeric expected",C_result);

    XL_readNumCellWD(XL_nbFactor, C_nbFactor, one, " ARM_ERR: Nb factor. : numeric expected",C_result);

   	if((XL_corrMatu->xltype == xltypeMissing) || (XL_corrMatu->xltype == xltypeNil))
	{
    }
	else
	{
        XL_readNumVector(XL_corrMatu,C_corrMatu,
          " ARM_ERR: Correlation maturities: array of numeric expected",C_result);
    }

   	if((XL_corrMatrix->xltype == xltypeMissing) || (XL_corrMatrix->xltype == xltypeNil))
	{
    }
	else
	{
        XL_readNumVector(XL_corrMatrix,C_corrMatrix,
          " ARM_ERR: Correlation matrix: array of numeric expected",C_result);
    }

   	if (( XL_corr2Matu->xltype == xltypeMissing ) 
		|| (XL_corr2Matu->xltype == xltypeNil))
	{
    }
	else
	{
        XL_readNumVector(XL_corr2Matu,C_corr2Matu,
          " ARM_ERR: Correlation correlated maturities : array of numeric expected",C_result);
    }

	if ( (C_autoMode = ARM_ConvAutoMode2(CS_autoMode, C_result)) == ARM_DEFAULT_ERR)
	{
		 ARM_ARG_ERR();

		 return(LPXLOPER)&XL_result;
	}

	if ( (C_fineMonth = ARM_ConvFineMode(CS_fineMonth, C_result)) == ARM_DEFAULT_ERR)
	{
		 ARM_ARG_ERR();

		 return(LPXLOPER)&XL_result;
	}

	if ( (C_mcMethod = ARM_ConvMcMethod(CS_mcMethod, C_result)) == ARM_DEFAULT_ERR)
	{
		 ARM_ARG_ERR();

		 return(LPXLOPER)&XL_result;
	}

	if ( (C_noControl = ARM_ConvInOut(CS_noControl, C_result)) == ARM_DEFAULT_ERR)
	{
		 ARM_ARG_ERR();

		 return(LPXLOPER)&XL_result;
	}

	long retCode;
	long objId;
	CCString prevClass;
	
	CCString curClass = LOCAL_FRMLSMONTECARLO_CLASS;
	CCString stringId = GetLastCurCellEnvValue ();

	if (!stringId)
	{
		retCode = ARMLOCAL_FRMLSMC_AUTO_2(LocalGetNumObjectId(C_zcId), LocalGetNumObjectId(C_swopVolId), LocalGetNumObjectId(C_swopSmileId),
									   LocalGetNumObjectId(C_irgVolId), LocalGetNumObjectId(C_irgSmileId),
									   C_autoMode,  C_horizon, C_fineMonth, C_nbTraj, C_mcMethod, C_noControl,
									   C_shapeDecay, C_shapeSlope, C_shapeAsymptote,
									   (long) C_nbFactor, C_corrMatu, C_corrMatrix, C_corr2Matu, C_result);

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
			retCode = ARMLOCAL_FRMLSMC_AUTO_2(LocalGetNumObjectId(C_zcId), LocalGetNumObjectId(C_swopVolId), LocalGetNumObjectId(C_swopSmileId),
										   LocalGetNumObjectId(C_irgVolId), LocalGetNumObjectId(C_irgSmileId),
										   C_autoMode,  C_horizon, C_fineMonth, C_nbTraj, C_mcMethod, C_noControl,
										   C_shapeDecay, C_shapeSlope, C_shapeAsymptote,
										   (long) C_nbFactor, C_corrMatu, C_corrMatrix, C_corr2Matu, C_result, objId);

			if ( retCode == ARM_OK )
			{
				LocalSetCurCellEnvValue (curClass, objId); 

				stringId = LocalMakeObjectId (objId, curClass);
			}
		}
		else
		{
			FreeCurCellContent ();

			retCode = ARMLOCAL_FRMLSMC_AUTO_2(LocalGetNumObjectId(C_zcId), LocalGetNumObjectId(C_swopVolId), LocalGetNumObjectId(C_swopSmileId),
										   LocalGetNumObjectId(C_irgVolId), LocalGetNumObjectId(C_irgSmileId),
										   C_autoMode,  C_horizon, C_fineMonth, C_nbTraj, C_mcMethod, C_noControl,
										   C_shapeDecay, C_shapeSlope, C_shapeAsymptote,
										   (long) C_nbFactor, C_corrMatu, C_corrMatrix, C_corr2Matu, C_result, objId);
	
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_FRMLSMC_AUTO_2" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}




__declspec(dllexport) LPXLOPER WINAPI Local_PXL_FRMLSMC_AUTO_2(LPXLOPER XL_zcId,
															   LPXLOPER XL_swopVolId,
															   LPXLOPER XL_swopSmileId,
															   LPXLOPER XL_irgVolId,
															   LPXLOPER XL_irgSmileId,
															   LPXLOPER XL_autoMode,
															   LPXLOPER XL_horizon,
															   LPXLOPER XL_nbTraj,
															   LPXLOPER XL_fineMonth,
															   LPXLOPER XL_mcMethod,
															   LPXLOPER XL_noControl,
															   LPXLOPER XL_shapeDecay,
															   LPXLOPER XL_shapeSlope,
															   LPXLOPER XL_shapeAsymptote,
															   LPXLOPER XL_nbFactor,
															   LPXLOPER XL_corrMatu,
															   LPXLOPER XL_corrMatrix,
															   LPXLOPER XL_corr2Matu)
{
	ADD_LOG("Local_PXL_FRMLSMC_AUTO_2");
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();

	// C variable
	CCString C_zcId;
    CCString C_swopVolId;
    CCString C_swopSmileId;
    CCString C_irgVolId;
    CCString C_irgSmileId;
    CCString CS_autoMode;
	long C_autoMode;
    CCString CS_fineMonth;
	long C_fineMonth;
    CCString CS_mcMethod;
	long C_mcMethod;
    CCString CS_noControl;
	long C_noControl;
    double C_horizon;
    double C_nbTraj;
    double C_shapeDecay;
    double C_shapeSlope;
    double C_shapeAsymptote;
    double C_nbFactor;
    
   
    VECTOR<double> C_corrMatu;
    VECTOR<double> C_corrMatrix;
    VECTOR<double> C_corr2Matu;	

    double mFineDefault = K_L_DEFAULT;
    double epsDefault = 1e-8;
    double huegeDefault = 1e8;
    double zero = 0;
    double one = 1;
    double thousand = 1000;


	
	// error
	static int error;
	static char* reason = "";


    XL_readStrCell(XL_zcId, C_zcId,
                   " ARM_ERR: Curve id: string expected",C_result);

    XL_readStrCell(XL_swopVolId, C_swopVolId,
                   " ARM_ERR: Swaption Vol id: string expected",C_result);

    XL_readStrCell(XL_swopSmileId, C_swopSmileId,
                   " ARM_ERR: Swaption Smile id: string expected",C_result);

    XL_readStrCell(XL_irgVolId, C_irgVolId,
                   " ARM_ERR: Cap Vol id: string expected",C_result);

    XL_readStrCell(XL_irgSmileId, C_irgSmileId,
                   " ARM_ERR: Cap Smile id: string expected",C_result);

    XL_readStrCell(XL_autoMode, CS_autoMode,
                   " ARM_ERR: AutoMode: string expected",C_result);

	XL_readNumCell(XL_horizon,C_horizon," ARM_ERR: Horizon date: date expected",C_result);

    XL_readNumCellWD(XL_nbTraj,C_nbTraj, thousand," ARM_ERR: Number of trajectories : numeric expected",C_result);

    XL_readStrCellWD(XL_fineMonth, CS_fineMonth, "L_DEFAULT", " ARM_ERR: fine month. : string expected",C_result);

    XL_readStrCellWD(XL_mcMethod, CS_mcMethod, "FAURE", 
                   " ARM_ERR: MC method : string expected",C_result);

    XL_readStrCellWD(XL_noControl, CS_noControl, "IN", 
                   " ARM_ERR: Control variate: string expected : IN or OUT",C_result);

    XL_readNumCellWD(XL_shapeDecay, C_shapeDecay, zero, " ARM_ERR: shape decay. : numeric expected",C_result);

    XL_readNumCellWD(XL_shapeSlope, C_shapeSlope, zero, " ARM_ERR: shape slope. : numeric expected",C_result);

    XL_readNumCellWD(XL_shapeAsymptote, C_shapeAsymptote, zero, " ARM_ERR: shape Asymptote. : numeric expected",C_result);

    XL_readNumCellWD(XL_nbFactor, C_nbFactor, one, " ARM_ERR: Nb factor. : numeric expected",C_result);

   	if((XL_corrMatu->xltype == xltypeMissing) || (XL_corrMatu->xltype == xltypeNil))
	{
    }
	else
	{
        XL_readNumVector(XL_corrMatu,C_corrMatu,
          " ARM_ERR: Correlation maturities: array of numeric expected",C_result);
    }

   	if((XL_corrMatrix->xltype == xltypeMissing) || (XL_corrMatrix->xltype == xltypeNil))
	{
    }
	else
	{
        XL_readNumVector(XL_corrMatrix,C_corrMatrix,
          " ARM_ERR: Correlation matrix: array of numeric expected",C_result);
    }

   	if((XL_corr2Matu->xltype == xltypeMissing) || (XL_corr2Matu->xltype == xltypeNil))
	{
    }
	else
	{
        XL_readNumVector(XL_corr2Matu,C_corr2Matu,
          " ARM_ERR: Correlation correlated maturities : array of numeric expected",C_result);
    }

	if ( (C_autoMode = ARM_ConvAutoMode(CS_autoMode, C_result)) == ARM_DEFAULT_ERR)
	{
		 ARM_ARG_ERR();

		 return(LPXLOPER)&XL_result;
	}

	if ( (C_fineMonth = ARM_ConvFineMode(CS_fineMonth, C_result)) == ARM_DEFAULT_ERR)
	{
		 ARM_ARG_ERR();

		 return(LPXLOPER)&XL_result;
	}

	if ( (C_mcMethod = ARM_ConvMcMethod(CS_mcMethod, C_result)) == ARM_DEFAULT_ERR)
	{
		 ARM_ARG_ERR();

		 return(LPXLOPER)&XL_result;
	}

	if ( (C_noControl = ARM_ConvInOut(CS_noControl, C_result)) == ARM_DEFAULT_ERR)
	{
		 ARM_ARG_ERR();

		 return(LPXLOPER)&XL_result;
	}

	long retCode;
	long objId;
	CCString prevClass;
	
	CCString curClass = LOCAL_FRMLSMONTECARLO_CLASS;
	CCString stringId;



    retCode = ARMLOCAL_FRMLSMC_AUTO_2(LocalGetNumObjectId(C_zcId), LocalGetNumObjectId(C_swopVolId), LocalGetNumObjectId(C_swopSmileId),
                                   LocalGetNumObjectId(C_irgVolId), LocalGetNumObjectId(C_irgSmileId),
                                   C_autoMode,  C_horizon, C_fineMonth, C_nbTraj, C_mcMethod, C_noControl,
                                   C_shapeDecay, C_shapeSlope, C_shapeAsymptote,
                                   (long) C_nbFactor, C_corrMatu, C_corrMatrix, C_corr2Matu, C_result);

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
		PXL_ARM_ERR();
	}

//	ARM_END();
	
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_PXL_FRMLSMC_AUTO_2" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}



__declspec(dllexport) LPXLOPER WINAPI Local_FRMLSMC_AUTO_G(LPXLOPER XL_zcId,
														   LPXLOPER XL_zdId,
														   LPXLOPER XL_swopVolId,
														   LPXLOPER XL_swopSmileId,
														   LPXLOPER XL_irgVolId,
														   LPXLOPER XL_irgSmileId,
														   LPXLOPER XL_autoMode,
														   LPXLOPER XL_horizon,
														   LPXLOPER XL_nbTraj,
														   LPXLOPER XL_fineMonth,
														   LPXLOPER XL_mcMethod,
														   LPXLOPER XL_noControl,
														   LPXLOPER XL_shapeDecay,
														   LPXLOPER XL_shapeSlope,
														   LPXLOPER XL_shapeAsymptote,
														   LPXLOPER XL_nbFactor,
														   LPXLOPER XL_corrMatu,
														   LPXLOPER XL_corrMatrix,
														   LPXLOPER XL_corr2Matu)
{
	ADD_LOG("Local_FRMLSMC_AUTO_G");
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();

	// C variable
	CCString C_zcId;
	CCString C_zdId;

    CCString C_swopVolId;
    CCString C_swopSmileId;
    CCString C_irgVolId;
    CCString C_irgSmileId;
    CCString CS_autoMode;
	long C_autoMode;
    CCString CS_fineMonth;
	long C_fineMonth;
    CCString CS_mcMethod;
	long C_mcMethod;
    CCString CS_noControl;
	long C_noControl;
    double C_horizon;
    double C_nbTraj;
    double C_shapeDecay;
    double C_shapeSlope;
    double C_shapeAsymptote;
    double C_nbFactor;

    VECTOR<double> C_corrMatu;
    VECTOR<double> C_corrMatrix;
    VECTOR<double> C_corr2Matu;	

    double mFineDefault = K_L_DEFAULT;
    double epsDefault = 1e-8;
    double huegeDefault = 1e8;
    double zero = 0;
    double one = 1;
    double thousand = 1000;

	// error
	static int error;
	static char* reason = "";

    XL_readStrCell(XL_zcId, C_zcId, " ARM_ERR: Domestic Curve id: string expected",C_result);

    XL_readStrCell(XL_zdId, C_zdId, " ARM_ERR: Domestic BS Curve id: string expected",C_result);

    XL_readStrCell(XL_swopVolId, C_swopVolId, " ARM_ERR: Swaption Vol id: string expected",C_result);

    XL_readStrCell(XL_swopSmileId, C_swopSmileId, " ARM_ERR: Swaption Smile id: string expected",C_result);

    XL_readStrCell(XL_irgVolId, C_irgVolId, " ARM_ERR: Cap Vol id: string expected",C_result);

    XL_readStrCell(XL_irgSmileId, C_irgSmileId, " ARM_ERR: Cap Smile id: string expected",C_result);

    XL_readStrCell(XL_autoMode, CS_autoMode, " ARM_ERR: AutoMode: string expected",C_result);

	XL_readNumCell(XL_horizon,C_horizon," ARM_ERR: Horizon date: date expected",C_result);

    XL_readNumCellWD(XL_nbTraj,C_nbTraj, thousand," ARM_ERR: Number of trajectories : numeric expected",C_result);

    XL_readStrCellWD(XL_fineMonth, CS_fineMonth, "L_DEFAULT", " ARM_ERR: fine month. : string expected",C_result);

    XL_readStrCellWD(XL_mcMethod, CS_mcMethod, "FAURE", " ARM_ERR: MC method : string expected",C_result);

    XL_readStrCellWD(XL_noControl, CS_noControl, "IN", " ARM_ERR: Control variate: string expected : IN or OUT",C_result);

    XL_readNumCellWD(XL_shapeDecay, C_shapeDecay, zero, " ARM_ERR: shape decay. : numeric expected",C_result);

    XL_readNumCellWD(XL_shapeSlope, C_shapeSlope, zero, " ARM_ERR: shape slope. : numeric expected",C_result);

    XL_readNumCellWD(XL_shapeAsymptote, C_shapeAsymptote, zero, " ARM_ERR: shape Asymptote. : numeric expected",C_result);

    XL_readNumCellWD(XL_nbFactor, C_nbFactor, one, " ARM_ERR: Nb factor. : numeric expected",C_result);

   	if((XL_corrMatu->xltype == xltypeMissing) || (XL_corrMatu->xltype == xltypeNil))
	{
    }
	else
	{
        XL_readNumVector(XL_corrMatu,C_corrMatu,
          " ARM_ERR: Correlation maturities: array of numeric expected",C_result);
    }

   	if((XL_corrMatrix->xltype == xltypeMissing) || (XL_corrMatrix->xltype == xltypeNil))
	{
    }
	else
	{
        XL_readNumVector(XL_corrMatrix,C_corrMatrix,
          " ARM_ERR: Correlation matrix: array of numeric expected",C_result);
    }

   	if (( XL_corr2Matu->xltype == xltypeMissing ) 
		|| (XL_corr2Matu->xltype == xltypeNil))
	{
    }
	else
	{
        XL_readNumVector(XL_corr2Matu,C_corr2Matu,
          " ARM_ERR: Correlation correlated maturities : array of numeric expected",C_result);
    }

	if ( (C_autoMode = ARM_ConvAutoMode2(CS_autoMode, C_result)) == ARM_DEFAULT_ERR)
	{
		 ARM_ARG_ERR();

		 return(LPXLOPER)&XL_result;
	}

	if ( (C_fineMonth = ARM_ConvFineMode(CS_fineMonth, C_result)) == ARM_DEFAULT_ERR)
	{
		 ARM_ARG_ERR();

		 return(LPXLOPER)&XL_result;
	}

	if ( (C_mcMethod = ARM_ConvMcMethod(CS_mcMethod, C_result)) == ARM_DEFAULT_ERR)
	{
		 ARM_ARG_ERR();

		 return(LPXLOPER)&XL_result;
	}

	if ( (C_noControl = ARM_ConvInOut(CS_noControl, C_result)) == ARM_DEFAULT_ERR)
	{
		 ARM_ARG_ERR();

		 return(LPXLOPER)&XL_result;
	}

	long retCode;
	long objId;
	CCString prevClass;
	
	CCString curClass = LOCAL_FRMLSMONTECARLO_CLASS;
	CCString stringId = GetLastCurCellEnvValue ();

	if (!stringId)
	{
		retCode = ARMLOCAL_FRMLSMC_AUTO_G(LocalGetNumObjectId(C_zcId), LocalGetNumObjectId(C_zdId), 
										LocalGetNumObjectId(C_swopVolId), LocalGetNumObjectId(C_swopSmileId),
									   LocalGetNumObjectId(C_irgVolId), LocalGetNumObjectId(C_irgSmileId),
									   C_autoMode,  C_horizon, C_fineMonth, C_nbTraj, C_mcMethod, C_noControl,
									   C_shapeDecay, C_shapeSlope, C_shapeAsymptote,
									   (long) C_nbFactor, C_corrMatu, C_corrMatrix, C_corr2Matu, C_result);

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
			retCode = ARMLOCAL_FRMLSMC_AUTO_G(LocalGetNumObjectId(C_zcId), LocalGetNumObjectId(C_zdId), 											LocalGetNumObjectId(C_swopVolId), LocalGetNumObjectId(C_swopSmileId),
										   LocalGetNumObjectId(C_irgVolId), LocalGetNumObjectId(C_irgSmileId),
										   C_autoMode,  C_horizon, C_fineMonth, C_nbTraj, C_mcMethod, C_noControl,
										   C_shapeDecay, C_shapeSlope, C_shapeAsymptote,
										   (long) C_nbFactor, C_corrMatu, C_corrMatrix, C_corr2Matu, C_result, objId);

			if ( retCode == ARM_OK )
			{
				LocalSetCurCellEnvValue (curClass, objId); 

				stringId = LocalMakeObjectId (objId, curClass);
			}
		}
		else
		{
			FreeCurCellContent ();

			retCode = ARMLOCAL_FRMLSMC_AUTO_G(LocalGetNumObjectId(C_zcId), LocalGetNumObjectId(C_zdId), 
											LocalGetNumObjectId(C_swopVolId), LocalGetNumObjectId(C_swopSmileId),
										   LocalGetNumObjectId(C_irgVolId), LocalGetNumObjectId(C_irgSmileId),
										   C_autoMode,  C_horizon, C_fineMonth, C_nbTraj, C_mcMethod, C_noControl,
										   C_shapeDecay, C_shapeSlope, C_shapeAsymptote,
										   (long) C_nbFactor, C_corrMatu, C_corrMatrix, C_corr2Matu, C_result, objId);
	
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_FRMLSMC_AUTO_G" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


__declspec(dllexport) LPXLOPER WINAPI Local_PXL_FRMLSMC_AUTO_G(LPXLOPER XL_zcId,
															   LPXLOPER XL_zdId,
															   LPXLOPER XL_swopVolId,
															   LPXLOPER XL_swopSmileId,
															   LPXLOPER XL_irgVolId,
															   LPXLOPER XL_irgSmileId,
															   LPXLOPER XL_autoMode,
															   LPXLOPER XL_horizon,
															   LPXLOPER XL_nbTraj,
															   LPXLOPER XL_fineMonth,
															   LPXLOPER XL_mcMethod,
															   LPXLOPER XL_noControl,
															   LPXLOPER XL_shapeDecay,
															   LPXLOPER XL_shapeSlope,
															   LPXLOPER XL_shapeAsymptote,
															   LPXLOPER XL_nbFactor,
															   LPXLOPER XL_corrMatu,
															   LPXLOPER XL_corrMatrix,
															   LPXLOPER XL_corr2Matu)                                            
{
	ADD_LOG("Local_PXL_FRMLSMC_AUTO_G");
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();

	// C variable
	CCString C_zcId;
	CCString C_zdId;

    CCString C_swopVolId;
    CCString C_swopSmileId;
    CCString C_irgVolId;
    CCString C_irgSmileId;
    CCString CS_autoMode;
	long C_autoMode;
    CCString CS_fineMonth;
	long C_fineMonth;
    CCString CS_mcMethod;
	long C_mcMethod;
    CCString CS_noControl;
	long C_noControl;
    double C_horizon;
    double C_nbTraj;
    double C_shapeDecay;
    double C_shapeSlope;
    double C_shapeAsymptote;
    double C_nbFactor;

    VECTOR<double> C_corrMatu;
    VECTOR<double> C_corrMatrix;
    VECTOR<double> C_corr2Matu;	

    double mFineDefault = K_L_DEFAULT;
    double epsDefault = 1e-8;
    double huegeDefault = 1e8;
    double zero = 0;
    double one = 1;
    double thousand = 1000;

	// error
	static int error;
	static char* reason = "";

    XL_readStrCell(XL_zcId, C_zcId, " ARM_ERR: Domestic Curve id: string expected",C_result);

    XL_readStrCell(XL_zdId, C_zdId, " ARM_ERR: Domestic BS Curve id: string expected",C_result);

    XL_readStrCell(XL_swopVolId, C_swopVolId, " ARM_ERR: Swaption Vol id: string expected",C_result);

    XL_readStrCell(XL_swopSmileId, C_swopSmileId, " ARM_ERR: Swaption Smile id: string expected",C_result);

    XL_readStrCell(XL_irgVolId, C_irgVolId, " ARM_ERR: Cap Vol id: string expected",C_result);

    XL_readStrCell(XL_irgSmileId, C_irgSmileId, " ARM_ERR: Cap Smile id: string expected",C_result);

    XL_readStrCell(XL_autoMode, CS_autoMode, " ARM_ERR: AutoMode: string expected",C_result);

	XL_readNumCell(XL_horizon,C_horizon," ARM_ERR: Horizon date: date expected",C_result);

    XL_readNumCellWD(XL_nbTraj,C_nbTraj, thousand," ARM_ERR: Number of trajectories : numeric expected",C_result);

    XL_readStrCellWD(XL_fineMonth, CS_fineMonth, "L_DEFAULT", " ARM_ERR: fine month. : string expected",C_result);

    XL_readStrCellWD(XL_mcMethod, CS_mcMethod, "FAURE", " ARM_ERR: MC method : string expected",C_result);

    XL_readStrCellWD(XL_noControl, CS_noControl, "IN", " ARM_ERR: Control variate: string expected : IN or OUT",C_result);

    XL_readNumCellWD(XL_shapeDecay, C_shapeDecay, zero, " ARM_ERR: shape decay. : numeric expected",C_result);

    XL_readNumCellWD(XL_shapeSlope, C_shapeSlope, zero, " ARM_ERR: shape slope. : numeric expected",C_result);

    XL_readNumCellWD(XL_shapeAsymptote, C_shapeAsymptote, zero, " ARM_ERR: shape Asymptote. : numeric expected",C_result);

    XL_readNumCellWD(XL_nbFactor, C_nbFactor, one, " ARM_ERR: Nb factor. : numeric expected",C_result);

   	if((XL_corrMatu->xltype == xltypeMissing) || (XL_corrMatu->xltype == xltypeNil))
	{
    }
	else
	{
        XL_readNumVector(XL_corrMatu,C_corrMatu,
          " ARM_ERR: Correlation maturities: array of numeric expected",C_result);
    }

   	if((XL_corrMatrix->xltype == xltypeMissing) || (XL_corrMatrix->xltype == xltypeNil))
	{
    }
	else
	{
        XL_readNumVector(XL_corrMatrix,C_corrMatrix,
          " ARM_ERR: Correlation matrix: array of numeric expected",C_result);
    }

   	if (( XL_corr2Matu->xltype == xltypeMissing ) 
		|| (XL_corr2Matu->xltype == xltypeNil))
	{
    }
	else
	{
        XL_readNumVector(XL_corr2Matu,C_corr2Matu,
          " ARM_ERR: Correlation correlated maturities : array of numeric expected",C_result);
    }

	if ( (C_autoMode = ARM_ConvAutoMode2(CS_autoMode, C_result)) == ARM_DEFAULT_ERR)
	{
		 ARM_ARG_ERR();

		 return(LPXLOPER)&XL_result;
	}

	if ( (C_fineMonth = ARM_ConvFineMode(CS_fineMonth, C_result)) == ARM_DEFAULT_ERR)
	{
		 ARM_ARG_ERR();

		 return(LPXLOPER)&XL_result;
	}

	if ( (C_mcMethod = ARM_ConvMcMethod(CS_mcMethod, C_result)) == ARM_DEFAULT_ERR)
	{
		 ARM_ARG_ERR();

		 return(LPXLOPER)&XL_result;
	}

	if ( (C_noControl = ARM_ConvInOut(CS_noControl, C_result)) == ARM_DEFAULT_ERR)
	{
		 ARM_ARG_ERR();

		 return(LPXLOPER)&XL_result;
	}

	long retCode;
	long objId;
	
	CCString curClass = LOCAL_FRMLSMONTECARLO_CLASS;
	CCString stringId;

	retCode = ARMLOCAL_FRMLSMC_AUTO_G(LocalGetNumObjectId(C_zcId), LocalGetNumObjectId(C_zdId), 
									LocalGetNumObjectId(C_swopVolId), LocalGetNumObjectId(C_swopSmileId),
								   LocalGetNumObjectId(C_irgVolId), LocalGetNumObjectId(C_irgSmileId),
								   C_autoMode,  C_horizon, C_fineMonth, C_nbTraj, C_mcMethod, C_noControl,
								   C_shapeDecay, C_shapeSlope, C_shapeAsymptote,
								   (long) C_nbFactor, C_corrMatu, C_corrMatrix, C_corr2Matu, C_result);

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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_PXL_FRMLSMC_AUTO_G" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}



__declspec(dllexport) LPXLOPER WINAPI Local_BSSLMODEL(LPXLOPER XL_date,
													  LPXLOPER XL_zerocurve,
													  LPXLOPER XL_volSpreadLock,
													  LPXLOPER XL_cvCapVol,
													  LPXLOPER XL_cvIndexVol)
{
	ADD_LOG("Local_BSSLMODEL");
//	ARM_BEGIN();

    static XLOPER XL_result;
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	    ARM_NOCALCIFWIZ();

	    // C variable
	    double C_date;
	    
	    double C_zerocurve_double;
	    CCString C_zerocurve_str;
	    long zerocurve_type;
	    
	    double C_volSpreadLock_double;
	    CCString C_volSpreadLock_str;
	    long volSpreadLock_type;
	    
	    double C_cvCapVol_double;
	    CCString C_cvCapVol_str;
	    long cvCapVol_type;

        double C_cvIndexVol_double;
	    double C_cvIndexVol_double_default = 0.0;
	    CCString C_cvIndexVol_str;
	    long cvIndexVol_type;
	    
	    // error
	    static int error;
	    static char* reason = "";

	    XL_readNumCell(XL_date,C_date," ARM_ERR: date: date expected",C_result);
	    XL_readStrOrNumCell(XL_zerocurve,C_zerocurve_str,C_zerocurve_double,zerocurve_type," ARM_ERR: curve id or rate: string or numeric expected",C_result);
	    XL_readStrOrNumCell(XL_volSpreadLock,C_volSpreadLock_str,C_volSpreadLock_double,volSpreadLock_type," ARM_ERR: volatility SPLK curve or volatility: string or numeric expected",C_result);
	    XL_readStrOrNumCell(XL_cvCapVol,C_cvCapVol_str,C_cvCapVol_double,cvCapVol_type," ARM_ERR: volatility curve or volatility: string or numeric expected",C_result);
        XL_readStrOrNumCellWD(XL_cvIndexVol,C_cvIndexVol_str,C_cvIndexVol_double,C_cvIndexVol_double_default,cvIndexVol_type," ARM_ERR: volatility curve or volatility: string or numeric expected",C_result);
	    
	    if ( zerocurve_type == XL_TYPE_STRING )
	    {
	       C_zerocurve_double = (double) LocalGetNumObjectId(C_zerocurve_str);

	       zerocurve_type = 1;
	    }
	    else
	    {
	       zerocurve_type = 0;
	    }


	    if ( volSpreadLock_type == XL_TYPE_STRING )
	    {
	       C_volSpreadLock_double = (double) LocalGetNumObjectId (C_volSpreadLock_str);

	       volSpreadLock_type = 1;
	    }
	    else
	    {
	       volSpreadLock_type = 0;
	    }

	    if ( cvCapVol_type == XL_TYPE_STRING )
	    {
	       C_cvCapVol_double = (double)LocalGetNumObjectId (C_cvCapVol_str);

           cvCapVol_type = 1;
	    }
	    else
	    {
	       cvCapVol_type = 0;
	    }

	    if((XL_cvIndexVol->xltype == xltypeMissing) || (XL_cvIndexVol->xltype == xltypeNil))
	    {
		    cvIndexVol_type = 0;
        }
	    else
        if ( cvIndexVol_type == XL_TYPE_STRING )
	    {
	       C_cvIndexVol_double = (double) LocalGetNumObjectId (C_cvIndexVol_str);

	       cvIndexVol_type = 1;
	    }
	    else
	    {
	       cvIndexVol_type = 0;
	    }

	    long retCode;
	    long objId;
	    CCString prevClass;
	    
	    CCString curClass = LOCAL_BSMODEL_CLASS;
	    CCString stringId = GetLastCurCellEnvValue ();

	    if (!stringId)
	    {
	       retCode = ARMLOCAL_bsslmodel(C_date, (long)zerocurve_type, C_zerocurve_double, 
							     (long) volSpreadLock_type, C_volSpreadLock_double, 
							     (long) cvCapVol_type, C_cvCapVol_double,
							     (long) cvIndexVol_type, C_cvIndexVol_double,
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
			    retCode = ARMLOCAL_bsslmodel(C_date, (long)zerocurve_type, C_zerocurve_double, 
							            (long) volSpreadLock_type, C_volSpreadLock_double, 
							            (long) cvCapVol_type, C_cvCapVol_double,
						    	        (long) cvIndexVol_type, C_cvIndexVol_double,
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

			    retCode = ARMLOCAL_bsslmodel(C_date, (long)zerocurve_type, C_zerocurve_double, 
							            (long) volSpreadLock_type, C_volSpreadLock_double, 
							            (long) cvCapVol_type, C_cvCapVol_double,
							            (long) cvIndexVol_type, C_cvIndexVol_double, C_result);
			    
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_BSSLMODEL" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}



__declspec(dllexport) LPXLOPER WINAPI Local_PXL_BSSLMODEL(LPXLOPER XL_date,
														  LPXLOPER XL_zerocurve,
														  LPXLOPER XL_volSpreadLock,
														  LPXLOPER XL_cvCapVol,
														  LPXLOPER XL_cvIndexVol)
{
	ADD_LOG("Local_PXL_BSSLMODEL");
//	ARM_BEGIN();

    static XLOPER XL_result;
	ARM_result C_result;
	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	    // C variable
	    double C_date;

	    double C_zerocurve_double;
	    CCString C_zerocurve_str;
	    long zerocurve_type;

	    double C_volSpreadLock_double;
	    CCString C_volSpreadLock_str;
	    long volSpreadLock_type;

	    double C_cvCapVol_double;
	    CCString C_cvCapVol_str;
	    long cvCapVol_type;

        double C_cvIndexVol_double;
	    double C_cvIndexVol_double_default = 0.0;
	    CCString C_cvIndexVol_str;
	    long cvIndexVol_type;

	    // error
	    static int error;
	    static char* reason = "";

	    XL_readNumCell(XL_date,C_date," ARM_ERR: date: date expected",C_result);
	    XL_readStrOrNumCell(XL_zerocurve,C_zerocurve_str,C_zerocurve_double,zerocurve_type," ARM_ERR: curve id or rate: string or numeric expected",C_result);
	    XL_readStrOrNumCell(XL_volSpreadLock,C_volSpreadLock_str,C_volSpreadLock_double,volSpreadLock_type," ARM_ERR: volatility SPLK curve or volatility: string or numeric expected",C_result);
	    XL_readStrOrNumCell(XL_cvCapVol,C_cvCapVol_str,C_cvCapVol_double,cvCapVol_type," ARM_ERR: volatility curve or volatility: string or numeric expected",C_result);
        XL_readStrOrNumCellWD(XL_cvIndexVol,C_cvIndexVol_str,C_cvIndexVol_double,C_cvIndexVol_double_default,cvIndexVol_type," ARM_ERR: volatility curve or volatility: string or numeric expected",C_result);
	    
	    if ( zerocurve_type == XL_TYPE_STRING )
	    {
	       C_zerocurve_double = (double) LocalGetNumObjectId(C_zerocurve_str);

	       zerocurve_type = 1;
	    }
	    else
	    {
	       zerocurve_type = 0;
	    }


	    if ( volSpreadLock_type == XL_TYPE_STRING )
	    {
	       C_volSpreadLock_double = (double) LocalGetNumObjectId (C_volSpreadLock_str);

	       volSpreadLock_type = 1;
	    }
	    else
	    {
	       volSpreadLock_type = 0;
	    }

	    if ( cvCapVol_type == XL_TYPE_STRING )
	    {
	       C_cvCapVol_double = (double)LocalGetNumObjectId (C_cvCapVol_str);

           cvCapVol_type = 1;
	    }
	    else
	    {
	       cvCapVol_type = 0;
	    }

	    if((XL_cvIndexVol->xltype == xltypeMissing) || (XL_cvIndexVol->xltype == xltypeNil))
	    {
		    cvIndexVol_type = 0;
        }
	    else
        if ( cvIndexVol_type == XL_TYPE_STRING )
	    {
	       C_cvIndexVol_double = (double) LocalGetNumObjectId (C_cvIndexVol_str);

	       cvIndexVol_type = 1;
	    }
	    else
	    {
	       cvIndexVol_type = 0;
	    }

	    long retCode;
	    long objId;
	    
	    CCString curClass = LOCAL_BSMODEL_CLASS;
	    CCString stringId;

       retCode = ARMLOCAL_bsslmodel(C_date, (long)zerocurve_type, C_zerocurve_double, 
						     (long) volSpreadLock_type, C_volSpreadLock_double, 
						     (long) cvCapVol_type, C_cvCapVol_double,
						     (long) cvIndexVol_type, C_cvIndexVol_double,
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_PXL_BSSLMODEL" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


__declspec(dllexport) LPXLOPER WINAPI Local_BKIRTREE (LPXLOPER XL_zc,
													  LPXLOPER XL_startDate,
													  LPXLOPER XL_endDate,
													  LPXLOPER XL_pas,
													  LPXLOPER XL_a,
													  LPXLOPER XL_sigma)
{
	ADD_LOG("Local_BKIRTREE ");
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();

	// C variable
	CCString C_zc;
	double C_startDate;
	double C_endDate;
	double C_pas;
	double C_a;
	double C_sigma;
	
	// error
	static int error;
	static char* reason = "";

	XL_readStrCell(XL_zc,C_zc," ARM_ERR: zero-coupon curve id: object expected",C_result);
	XL_readNumCell(XL_startDate,C_startDate," ARM_ERR: start date: date expected",C_result);
	XL_readNumCell(XL_endDate,C_endDate," ARM_ERR: end date: date expected",C_result);
	XL_readNumCell(XL_pas,C_pas," ARM_ERR: number of steps: numeric expected",C_result);
	XL_readNumCell(XL_a,C_a," ARM_ERR: a: numeric expected",C_result);
	XL_readNumCell(XL_sigma,C_sigma," ARM_ERR: sigma: numeric expected",C_result);
	
	long retCode;
	long objId;
	CCString prevClass;
	
	CCString curClass = LOCAL_BKIRTREE_CLASS;
	CCString stringId = GetLastCurCellEnvValue ();

	if (!stringId)
	{
		retCode = ARMLOCAL_BKIRTREE (LocalGetNumObjectId (C_zc),
							    C_startDate,
							    C_endDate,
							    (long)C_pas,
							    C_a,
							    C_sigma,
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
			
		if(curClass == prevClass)
		{
			retCode = ARMLOCAL_BKIRTREE (LocalGetNumObjectId (C_zc),
									C_startDate,
									C_endDate,
									(long)C_pas,
									C_a,
									C_sigma,
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
			retCode = ARMLOCAL_BKIRTREE (LocalGetNumObjectId (C_zc),
									C_startDate,
									C_endDate,
									(long)C_pas,
									C_a,
									C_sigma,
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

//	ARM_END();
	
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_BKIRTREE" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


__declspec(dllexport) LPXLOPER WINAPI Local_PXL_BKIRTREE (LPXLOPER XL_zc,
														  LPXLOPER XL_startDate,
														  LPXLOPER XL_endDate,
														  LPXLOPER XL_pas,
														  LPXLOPER XL_a,
														  LPXLOPER XL_sigma)
{
	ADD_LOG("Local_PXL_BKIRTREE ");
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();

	// C variable
	CCString C_zc;
	double C_startDate;
	double C_endDate;
	double C_pas;
	double C_a;
	double C_sigma;

	// error
	static int error;
	static char* reason = "";

	XL_readStrCell(XL_zc,C_zc," ARM_ERR: zero-coupon curve id: object expected",C_result);
	XL_readNumCell(XL_startDate,C_startDate," ARM_ERR: start date: date expected",C_result);
	XL_readNumCell(XL_endDate,C_endDate," ARM_ERR: end date: date expected",C_result);
	XL_readNumCell(XL_pas,C_pas," ARM_ERR: number of steps: numeric expected",C_result);
	XL_readNumCell(XL_a,C_a," ARM_ERR: a: numeric expected",C_result);
	XL_readNumCell(XL_sigma,C_sigma," ARM_ERR: sigma: numeric expected",C_result);
	
	long retCode;
	long objId;

	CCString curClass = LOCAL_BKIRTREE_CLASS;
	CCString stringId;

	retCode = ARMLOCAL_BKIRTREE (LocalGetNumObjectId (C_zc),
							C_startDate,
							C_endDate,
							(long)C_pas,
							C_a,
							C_sigma,
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_PXL_BKIRTREE" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}



__declspec(dllexport) LPXLOPER WINAPI Local_LDCMC_FROM_ANA(LPXLOPER XL_anaModId,
														   LPXLOPER XL_horizon,
														   LPXLOPER XL_nbTraj,
														   LPXLOPER XL_mcmethod,
														   LPXLOPER XL_pricerType) 
                                           
{
	ADD_LOG("Local_LDCMC_FROM_ANA");
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();

	// C variable
	CCString C_anaModId;

	double C_horizon;

	double C_nbTraj;

	double C_mcmethod;
	double C_mcmethod_default=0.0;

	CCString C_pricerType;
	long pricerTypeId;

	// error
	static int error;
	static char* reason = "";

	XL_readStrCell(XL_anaModId,C_anaModId," ARM_ERR: Model id: string expected",C_result);
	XL_readNumCell(XL_horizon,C_horizon," ARM_ERR: Horizon date: date expected",C_result);
	XL_readNumCell(XL_nbTraj, C_nbTraj," ARM_ERR: number of traj. : numeric expected",C_result);
	XL_readNumCellWD(XL_mcmethod, C_mcmethod, C_mcmethod_default, " ARM_ERR: Monte Carlo Method. : numeric expected",C_result);
	XL_readStrCellWD(XL_pricerType,C_pricerType,"RNMC"," ARM_ERR: Pricer type: string expected",C_result);

	if ( (pricerTypeId = ARM_ConvPricerType(C_pricerType, C_result)) == ARM_DEFAULT_ERR)
	{
		 ARM_ARG_ERR();

		 return(LPXLOPER)&XL_result;
	}

	long retCode;
	long objId;
	CCString prevClass;

	CCString curClass = LOCAL_LOGDECMONTECARLO_CLASS;
	CCString stringId = GetLastCurCellEnvValue ();

	if (!stringId)
	{

        retCode = ARMLOCAL_LDCMC_FROM_ANA(LocalGetNumObjectId (C_anaModId),
										  C_horizon,
										  (long) C_nbTraj,
										  (long) C_mcmethod,
										  pricerTypeId,
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
			retCode = ARMLOCAL_LDCMC_FROM_ANA(LocalGetNumObjectId (C_anaModId),
											  C_horizon,
											  (long) C_nbTraj,
											  (long) C_mcmethod,
											  pricerTypeId,
											  C_result,
											  objId);

			if ( retCode == ARM_OK )
			{
			   objId = C_result.getLong ();

			   stringId = LocalMakeObjectId (objId, curClass);
			}
		}
		else
		{
			FreeCurCellContent ();

            retCode = ARMLOCAL_LDCMC_FROM_ANA(LocalGetNumObjectId (C_anaModId),
											  C_horizon,
											  (long) C_nbTraj,
											  (long) C_mcmethod,
											  pricerTypeId,
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_LDCMC_FROM_ANA" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


__declspec(dllexport) LPXLOPER WINAPI Local_PXL_LDCMC_FROM_ANA(LPXLOPER XL_anaModId,
															   LPXLOPER XL_horizon,
															   LPXLOPER XL_nbTraj,
															   LPXLOPER XL_mcmethod,
															   LPXLOPER XL_pricerType) 
{
	ADD_LOG("Local_PXL_LDCMC_FROM_ANA");
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();

	// C variable
	CCString C_anaModId;

	double C_horizon;

	double C_nbTraj;

	double C_mcmethod;
	double C_mcmethod_default=0.0;

	CCString C_pricerType;
	long pricerTypeId;

	// error
	static int error;
	static char* reason = "";

	XL_readStrCell(XL_anaModId,C_anaModId," ARM_ERR: Model id: string expected",C_result);
	XL_readNumCell(XL_horizon,C_horizon," ARM_ERR: Horizon date: date expected",C_result);
	XL_readNumCell(XL_nbTraj, C_nbTraj," ARM_ERR: number of traj. : numeric expected",C_result);
	XL_readNumCellWD(XL_mcmethod, C_mcmethod, C_mcmethod_default, " ARM_ERR: Monte Carlo Method. : numeric expected",C_result);
	XL_readStrCellWD(XL_pricerType,C_pricerType,"RNMC"," ARM_ERR: Pricer type: string expected",C_result);

	if ( (pricerTypeId = ARM_ConvPricerType(C_pricerType, C_result)) == ARM_DEFAULT_ERR)
	{
		 ARM_ARG_ERR();

		 return(LPXLOPER)&XL_result;
	}

	long retCode;
	long objId;

	CCString curClass = LOCAL_LOGDECMONTECARLO_CLASS;
	CCString stringId;

    retCode = ARMLOCAL_LDCMC_FROM_ANA(LocalGetNumObjectId (C_anaModId), C_horizon,
                                 (long) C_nbTraj, (long) C_mcmethod, pricerTypeId, C_result);

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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_PXL_LDCMC_FROM_ANA" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}



__declspec(dllexport) LPXLOPER WINAPI Local_HWTREE(LPXLOPER XL_zcId,
												   LPXLOPER XL_begDate,
												   LPXLOPER XL_endDate,
												   LPXLOPER XL_nbSteps,
												   LPXLOPER XL_a,
												   LPXLOPER XL_sigma)
{
	ADD_LOG("Local_HWTREE");
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();

	// C variable
	CCString C_zcId;
   
    double C_begDate;
	double C_endDate;	

	double C_nbSteps;
	double C_a;
	double C_sigma;
	
	// error
	static int error;
	static char* reason = "";

    XL_readStrCell(XL_zcId, C_zcId," ARM_ERR: Curve id: string expected",C_result);
	XL_readNumCell(XL_begDate,C_begDate," ARM_ERR: start date: date expected",C_result);
    XL_readNumCell(XL_endDate,C_endDate," ARM_ERR: end date: date expected",C_result);
    XL_readNumCell(XL_nbSteps, C_nbSteps," ARM_ERR: number of steps: numeric expected",C_result);
	XL_readNumCell(XL_sigma,C_sigma," ARM_ERR: sigma: numeric expected",C_result);
    XL_readNumCell(XL_a, C_a," ARM_ERR: a: numeric expected",C_result);

	long retCode;
	long objId;
	CCString prevClass;

	CCString curClass = LOCAL_YIELD_CURVE_HWTREE_CLASS;
	CCString stringId = GetLastCurCellEnvValue ();

	if (!stringId)
	{
        retCode = ARMLOCAL_HWTREE(LocalGetNumObjectId (C_zcId), C_begDate, C_endDate,
                             (long) C_nbSteps, C_a, C_sigma, C_result);

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
           retCode = ARMLOCAL_HWTREE(LocalGetNumObjectId (C_zcId), C_begDate, C_endDate,
                                (long) C_nbSteps, C_a, C_sigma, C_result, objId);

			if ( retCode == ARM_OK )
			{
			   LocalSetCurCellEnvValue (curClass, objId); 

			   stringId = LocalMakeObjectId (objId, curClass);
			}
		}
		else
		{
			FreeCurCellContent ();

            retCode = ARMLOCAL_HWTREE(LocalGetNumObjectId (C_zcId), C_begDate, C_endDate,
                                 (long) C_nbSteps, C_a, C_sigma, C_result);
	
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_HWTREE" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


__declspec(dllexport) LPXLOPER WINAPI Local_PXL_HWTREE(LPXLOPER XL_zcId,
													   LPXLOPER XL_begDate,
													   LPXLOPER XL_endDate,
													   LPXLOPER XL_nbSteps,
													   LPXLOPER XL_a,
													   LPXLOPER XL_sigma)
{
	ADD_LOG("Local_PXL_HWTREE");
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();

	// C variable
	CCString C_zcId;
   
    double C_begDate;
	double C_endDate;	

	double C_nbSteps;
	double C_a;
	double C_sigma;
	
	// error
	static int error;
	static char* reason = "";

    XL_readStrCell(XL_zcId, C_zcId," ARM_ERR: Curve id: string expected",C_result);
	XL_readNumCell(XL_begDate,C_begDate," ARM_ERR: start date: date expected",C_result);
    XL_readNumCell(XL_endDate,C_endDate," ARM_ERR: end date: date expected",C_result);
    XL_readNumCell(XL_nbSteps, C_nbSteps," ARM_ERR: number of steps: numeric expected",C_result);
	XL_readNumCell(XL_sigma,C_sigma," ARM_ERR: sigma: numeric expected",C_result);
    XL_readNumCell(XL_a, C_a," ARM_ERR: a: numeric expected",C_result);

	long retCode;
	long objId;

	CCString curClass = LOCAL_YIELD_CURVE_HWTREE_CLASS;
	CCString stringId;

    retCode = ARMLOCAL_HWTREE(LocalGetNumObjectId (C_zcId), C_begDate, C_endDate,
                         (long) C_nbSteps, C_a, C_sigma, C_result);

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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_PXL_HWTREE" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}



__declspec(dllexport) LPXLOPER WINAPI Local_HWSIGCONST(LPXLOPER XL_zcId,
													   LPXLOPER XL_begDate,
													   LPXLOPER XL_endDate,
													   LPXLOPER XL_nbSteps, 
													   LPXLOPER XL_a,
													   LPXLOPER XL_sigma)
{
	ADD_LOG("Local_HWSIGCONST");
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();

	// C variable
	CCString C_zcId;
   
    double C_begDate;
	double C_endDate;
	
	double C_nbSteps;
	double C_a;
	double C_sigma;

	// error
	static int error;
	static char* reason = "";


    XL_readStrCell(XL_zcId, C_zcId," ARM_ERR: Curve id: string expected",C_result);
	XL_readNumCell(XL_begDate,C_begDate," ARM_ERR: start date: date expected",C_result);
    XL_readNumCell(XL_endDate,C_endDate," ARM_ERR: end date: date expected",C_result);
    XL_readNumCell(XL_nbSteps, C_nbSteps," ARM_ERR: number of steps: numeric expected",C_result);
	XL_readNumCell(XL_sigma,C_sigma," ARM_ERR: sigma: numeric expected",C_result);
    XL_readNumCell(XL_a, C_a," ARM_ERR: a: numeric expected",C_result);

	long retCode;
	long objId;
	CCString prevClass;
	
	CCString curClass = LOCAL_HWSIGCONST_TREE_CLASS;
	CCString stringId = GetLastCurCellEnvValue ();

	if (!stringId)
	{
        retCode = ARMLOCAL_HWSIGCONST(LocalGetNumObjectId (C_zcId), C_begDate, C_endDate,
                                 (long) C_nbSteps, C_a, C_sigma, C_result);

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
			retCode = ARMLOCAL_HWSIGCONST(LocalGetNumObjectId (C_zcId), C_begDate, C_endDate,
									(long) C_nbSteps, C_a, C_sigma, C_result, objId);

			if ( retCode == ARM_OK )
			{
				LocalSetCurCellEnvValue (curClass, objId); 

				stringId = LocalMakeObjectId (objId, curClass);
			}
		}
		else
		{
			FreeCurCellContent ();

            retCode = ARMLOCAL_HWSIGCONST(LocalGetNumObjectId (C_zcId), C_begDate, C_endDate,
                                     (long) C_nbSteps, C_a, C_sigma, C_result);
	
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_HWSIGCONST" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


__declspec(dllexport) LPXLOPER WINAPI Local_PXL_HWSIGCONST(LPXLOPER XL_zcId,
														   LPXLOPER XL_begDate,
														   LPXLOPER XL_endDate,
														   LPXLOPER XL_nbSteps,
														   LPXLOPER XL_a,
														   LPXLOPER XL_sigma)
{
	ADD_LOG("Local_PXL_HWSIGCONST");
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();

	// C variable
	CCString C_zcId;

	double C_begDate;
	double C_endDate;

	double C_nbSteps;
	double C_a;
	double C_sigma;

	// error
	static int error;
	static char* reason = "";


    XL_readStrCell(XL_zcId, C_zcId," ARM_ERR: Curve id: string expected",C_result);
	XL_readNumCell(XL_begDate,C_begDate," ARM_ERR: start date: date expected",C_result);
    XL_readNumCell(XL_endDate,C_endDate," ARM_ERR: end date: date expected",C_result);
    XL_readNumCell(XL_nbSteps, C_nbSteps," ARM_ERR: number of steps: numeric expected",C_result);
	XL_readNumCell(XL_sigma,C_sigma," ARM_ERR: sigma: numeric expected",C_result);
    XL_readNumCell(XL_a, C_a," ARM_ERR: a: numeric expected",C_result);

	long retCode;
	long objId;
	
	CCString curClass = LOCAL_HWSIGCONST_TREE_CLASS;
	CCString stringId;

    retCode = ARMLOCAL_HWSIGCONST(LocalGetNumObjectId (C_zcId), C_begDate, C_endDate,
                             (long) C_nbSteps, C_a, C_sigma, C_result);

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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_PXL_HWSIGCONST" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


__declspec(dllexport) LPXLOPER WINAPI Local_HWSIGVAR(LPXLOPER XL_zcId,
													 LPXLOPER XL_begDate,
													 LPXLOPER XL_endDate,
													 LPXLOPER XL_nbSteps,
													 LPXLOPER XL_a,
													 LPXLOPER XL_datesSigma)
{
	ADD_LOG("Local_HWSIGVAR");
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();

	// C variable
	CCString C_zcId;

	double C_begDate;
	double C_endDate;	

	double C_nbSteps;
	double C_a;
    VECTOR<double> C_datesSigma;

    VECTOR<CCString> C_dates;
    VECTOR<double>   C_sigmas;

	// error
	static int error;
	static char* reason = "";

    XL_readStrCell(XL_zcId, C_zcId," ARM_ERR: Curve id: string expected",C_result);
	XL_readNumCell(XL_begDate,C_begDate," ARM_ERR: start date: date expected",C_result);
    XL_readNumCell(XL_endDate,C_endDate," ARM_ERR: end date: date expected",C_result);
    XL_readNumCell(XL_nbSteps, C_nbSteps," ARM_ERR: number of steps: numeric expected",C_result);
    XL_readNumCell(XL_a, C_a," ARM_ERR: a: numeric expected",C_result);
    XL_readNumVector(XL_datesSigma, C_datesSigma," ARM_ERR: [date, sigma]: array expected",C_result);

	double tmp_rate = -1;
	
	if ( C_datesSigma.size () % 2 != 0 )
	{
		C_result.setMsg (" ARM_ERR: check your [Dates, Sigma(s)] array");

		ARM_ARG_ERR();

		return(LPXLOPER)&XL_result;
	}

	long real_size = C_datesSigma.size () / 2;
	long j = 0;

	char* sDate = new char[11];

	for (int i = 0; i < real_size; i++)
	{
		Local_XLDATE2ARMDATE(C_datesSigma[j],sDate);

		C_dates.push_back(sDate);
		C_sigmas.push_back(C_datesSigma[j+1]);
		j += 2;
	}
	delete [] sDate;

    if ( real_size == 1 )
    {
       C_dates.push_back(C_dates[0]);
       C_sigmas.push_back(C_sigmas[0]);
       
       real_size = 2;
    }

	long retCode;
	long objId;
	CCString prevClass;
	
	CCString curClass = LOCAL_HWSIGVAR_TREE_CLASS;
	CCString stringId = GetLastCurCellEnvValue();

	if (!stringId)
	{
        retCode = ARMLOCAL_HWSIGVAR(LocalGetNumObjectId (C_zcId), C_begDate, C_endDate,
                               (long) C_nbSteps, C_a, real_size,
                               C_dates, C_sigmas, C_result);

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
           retCode = ARMLOCAL_HWSIGVAR(LocalGetNumObjectId (C_zcId), 
                                  C_begDate, C_endDate,
                                  (long) C_nbSteps, C_a, real_size,
                                  C_dates, C_sigmas,
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

            retCode = ARMLOCAL_HWSIGVAR(LocalGetNumObjectId (C_zcId), C_begDate, C_endDate,
                                   (long) C_nbSteps, C_a, real_size,
                                   C_dates, C_sigmas, C_result);

	
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_HWSIGVAR" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


__declspec(dllexport) LPXLOPER WINAPI Local_PXL_HWSIGVAR(LPXLOPER XL_zcId,
														 LPXLOPER XL_begDate,
														 LPXLOPER XL_endDate,
														 LPXLOPER XL_nbSteps,
														 LPXLOPER XL_a,
														 LPXLOPER XL_datesSigma)
{
	ADD_LOG("Local_PXL_HWSIGVAR");
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();

	// C variable
	CCString C_zcId;

	double C_begDate;
	double C_endDate;	

	double C_nbSteps;
	double C_a;
    VECTOR<double> C_datesSigma;

    VECTOR<CCString> C_dates;
    VECTOR<double>   C_sigmas;

	// error
	static int error;
	static char* reason = "";

    XL_readStrCell(XL_zcId, C_zcId," ARM_ERR: Curve id: string expected",C_result);
	XL_readNumCell(XL_begDate,C_begDate," ARM_ERR: start date: date expected",C_result);
    XL_readNumCell(XL_endDate,C_endDate," ARM_ERR: end date: date expected",C_result);
    XL_readNumCell(XL_nbSteps, C_nbSteps," ARM_ERR: number of steps: numeric expected",C_result);
    XL_readNumCell(XL_a, C_a," ARM_ERR: a: numeric expected",C_result);
    XL_readNumVector(XL_datesSigma, C_datesSigma," ARM_ERR: [date, sigma]: array expected",C_result);

	double tmp_rate = -1;
	
	if ( C_datesSigma.size () % 2 != 0 )
	{
		C_result.setMsg (" ARM_ERR: check your [Dates, Sigma(s)] array");

		ARM_ARG_ERR();

		return(LPXLOPER)&XL_result;
	}

	long real_size = C_datesSigma.size () / 2;
	long j = 0;

	char* sDate = new char[11];

	for (int i = 0; i < real_size; i++)
	{
		Local_XLDATE2ARMDATE(C_datesSigma[j],sDate);

		C_dates.push_back(sDate);
		C_sigmas.push_back(C_datesSigma[j+1]);
		j += 2;
	}
	delete [] sDate;

    if ( real_size == 1 )
    {
       C_dates.push_back(C_dates[0]);
       C_sigmas.push_back(C_sigmas[0]);
       
       real_size = 2;
    }

	long retCode;
	long objId;
	
	CCString curClass = LOCAL_HWSIGVAR_TREE_CLASS;
	CCString stringId;

    retCode = ARMLOCAL_HWSIGVAR(LocalGetNumObjectId (C_zcId), C_begDate, C_endDate,
                           (long) C_nbSteps, C_a, real_size,
                           C_dates, C_sigmas, C_result);

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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_PXL_HWSIGVAR" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


__declspec(dllexport) LPXLOPER WINAPI Local_DFGYC (LPXLOPER XL_dMeanRevSpeed,
												   LPXLOPER XL_fMeanRevSpeed,
												   LPXLOPER XL_dSigma,
												   LPXLOPER XL_fSigma,
												   LPXLOPER XL_fxCorr,
												   LPXLOPER XL_fxVol,
												   LPXLOPER XL_ratesCorr,
												   LPXLOPER XL_dZc,
												   LPXLOPER XL_fZc)
{
	ADD_LOG("Local_DFGYC ");
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();

	// C variable
	double C_dMeanRevSpeed;
	double C_fMeanRevSpeed;
	double C_dSigma;
	double C_fSigma;
	double C_fxCorr;
	double C_fxVol;
	double C_ratesCorr;
	CCString C_dZc;
	CCString C_fZc;
	
	// error
	static int error;
	static char* reason = "";

	XL_readNumCell(XL_dMeanRevSpeed,C_dMeanRevSpeed," ARM_ERR: domestic mean reversion speed: numeric expected",C_result);
	XL_readNumCell(XL_fMeanRevSpeed,C_fMeanRevSpeed," ARM_ERR: foreign mean reversion speed: numeric expected",C_result);
	XL_readNumCell(XL_dSigma,C_dSigma," ARM_ERR: domestic sigma: numeric expected",C_result);
	XL_readNumCell(XL_fSigma,C_fSigma," ARM_ERR: foreign sigma: numeric expected",C_result);
	XL_readNumCell(XL_fxCorr,C_fxCorr," ARM_ERR: forex correlation: numeric expected",C_result);
	XL_readNumCell(XL_fxVol,C_fxVol," ARM_ERR: forex volatility: numeric expected",C_result);
	XL_readNumCell(XL_ratesCorr,C_ratesCorr," ARM_ERR: rates correlation: numeric expected",C_result);
	XL_readStrCell(XL_dZc,C_dZc," ARM_ERR: domestic zero-coupon curve id: object expected",C_result);
	XL_readStrCell(XL_fZc,C_fZc," ARM_ERR: foreign zero-coupon curve id: object expected",C_result);
	
	long retCode;
	long objId;
	CCString prevClass;
	
	CCString curClass = LOCAL_DFG_YIELD_CURVE_HWTREE_CLASS;
	CCString stringId = GetLastCurCellEnvValue ();

	if (!stringId)
	{
		retCode = ARMLOCAL_DFGYC (C_dMeanRevSpeed,
							 C_fMeanRevSpeed,
							 C_dSigma,
							 C_fSigma,
							 C_fxCorr,
							 C_fxVol,
							 C_ratesCorr,
							 LocalGetNumObjectId (C_dZc),
							 LocalGetNumObjectId (C_fZc),
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
			
		if(curClass == prevClass)
		{
			retCode = ARMLOCAL_DFGYC (C_dMeanRevSpeed,
								 C_fMeanRevSpeed,
								 C_dSigma,
								 C_fSigma,
								 C_fxCorr,
								 C_fxVol,
								 C_ratesCorr,
								 LocalGetNumObjectId (C_dZc),
								 LocalGetNumObjectId (C_fZc),
								 C_result,
								 objId);
		}
		else
		{
			FreeCurCellContent ();
			retCode = ARMLOCAL_DFGYC (C_dMeanRevSpeed,
								 C_fMeanRevSpeed,
								 C_dSigma,
								 C_fSigma,
								 C_fxCorr,
								 C_fxVol,
								 C_ratesCorr,
								 LocalGetNumObjectId (C_dZc),
								 LocalGetNumObjectId (C_fZc),
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

//	ARM_END();
	
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_DFGYC" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


__declspec(dllexport) LPXLOPER WINAPI Local_PXL_DFGYC (LPXLOPER XL_dMeanRevSpeed,
													   LPXLOPER XL_fMeanRevSpeed,
													   LPXLOPER XL_dSigma,
													   LPXLOPER XL_fSigma,
													   LPXLOPER XL_fxCorr,
													   LPXLOPER XL_fxVol,
													   LPXLOPER XL_ratesCorr,
													   LPXLOPER XL_dZc,
													   LPXLOPER XL_fZc)
{
	ADD_LOG("Local_PXL_DFGYC ");
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();

	// C variable
	double C_dMeanRevSpeed;
	double C_fMeanRevSpeed;
	double C_dSigma;
	double C_fSigma;
	double C_fxCorr;
	double C_fxVol;
	double C_ratesCorr;
	CCString C_dZc;
	CCString C_fZc;
	
	// error
	static int error;
	static char* reason = "";

	XL_readNumCell(XL_dMeanRevSpeed,C_dMeanRevSpeed," ARM_ERR: domestic mean reversion speed: numeric expected",C_result);
	XL_readNumCell(XL_fMeanRevSpeed,C_fMeanRevSpeed," ARM_ERR: foreign mean reversion speed: numeric expected",C_result);
	XL_readNumCell(XL_dSigma,C_dSigma," ARM_ERR: domestic sigma: numeric expected",C_result);
	XL_readNumCell(XL_fSigma,C_fSigma," ARM_ERR: foreign sigma: numeric expected",C_result);
	XL_readNumCell(XL_fxCorr,C_fxCorr," ARM_ERR: forex correlation: numeric expected",C_result);
	XL_readNumCell(XL_fxVol,C_fxVol," ARM_ERR: forex volatility: numeric expected",C_result);
	XL_readNumCell(XL_ratesCorr,C_ratesCorr," ARM_ERR: rates correlation: numeric expected",C_result);
	XL_readStrCell(XL_dZc,C_dZc," ARM_ERR: domestic zero-coupon curve id: object expected",C_result);
	XL_readStrCell(XL_fZc,C_fZc," ARM_ERR: foreign zero-coupon curve id: object expected",C_result);
	
	long retCode;
	long objId;
	
	CCString curClass = LOCAL_DFG_YIELD_CURVE_HWTREE_CLASS;
	CCString stringId;

	retCode = ARMLOCAL_DFGYC (C_dMeanRevSpeed,
						 C_fMeanRevSpeed,
						 C_dSigma,
						 C_fSigma,
						 C_fxCorr,
						 C_fxVol,
						 C_ratesCorr,
						 LocalGetNumObjectId (C_dZc),
						 LocalGetNumObjectId (C_fZc),
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_PXL_DFGYC" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


__declspec(dllexport) LPXLOPER WINAPI Local_IDTREEHW(LPXLOPER XL_startDate,
													 LPXLOPER XL_horizon,
													 LPXLOPER XL_nbSteps,
													 LPXLOPER XL_dMeanRevSpeed,
													 LPXLOPER XL_fMeanRevSpeed,
													 LPXLOPER XL_dSigma,
													 LPXLOPER XL_fSigma,
													 LPXLOPER XL_prtyCorr,
													 LPXLOPER XL_prtyVol,
													 LPXLOPER XL_ratesCorr,
													 LPXLOPER XL_dZcId,
													 LPXLOPER XL_fZcId)
{
	ADD_LOG("Local_IDTREEHW");
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();

	// C variable

	double C_startDate;
	double C_horizon;

	double C_nbSteps;

	double C_dMeanRevSpeed;
	double C_fMeanRevSpeed;

	double C_dSigma;
	double C_fSigma;

	double C_prtyCorr;
	double C_prtyVol;
	double C_ratesCorr;

	CCString C_dZcId;
	CCString C_fZcId;

	// error
	static int error;
	static char* reason = "";

	XL_readNumCell(XL_startDate,C_startDate," ARM_ERR: start date: date expected",C_result);
	XL_readNumCell(XL_horizon, C_horizon," ARM_ERR: horizon date: date expected",C_result);
	XL_readNumCell(XL_nbSteps, C_nbSteps," ARM_ERR: number of steps: numeric expected",C_result);
	XL_readNumCell(XL_dMeanRevSpeed,C_dMeanRevSpeed," ARM_ERR: domestic mean reversion: numeric expected",C_result);
	XL_readNumCell(XL_fMeanRevSpeed,C_fMeanRevSpeed," ARM_ERR: foreign mean reversion: numeric expected",C_result);
	XL_readNumCell(XL_dSigma, C_dSigma," ARM_ERR: domestic sigma: numeric expected",C_result);
	XL_readNumCell(XL_fSigma, C_fSigma," ARM_ERR: foreign sigma: numeric expected",C_result);
	XL_readNumCell(XL_prtyCorr, C_prtyCorr," ARM_ERR: correlation: numeric expected",C_result);
	XL_readNumCell(XL_prtyVol, C_prtyVol," ARM_ERR: volatility: numeric expected",C_result);
	XL_readNumCell(XL_ratesCorr, C_ratesCorr," ARM_ERR: rates correlation: numeric expected",C_result);
	XL_readStrCell(XL_dZcId, C_dZcId," ARM_ERR: domestic curve id: object expected",C_result);
	XL_readStrCell(XL_fZcId, C_fZcId," ARM_ERR: foreign curve id: object expected",C_result);

	long retCode;
	long objId;
	CCString prevClass;
	
	CCString curClass = LOCAL_I3DTREE_HWTREE_CLASS;
	CCString stringId = GetLastCurCellEnvValue ();

	if (!stringId)
	{
        retCode = ARMLOCAL_IDTREEHW(C_startDate, C_horizon, (long) C_nbSteps,
                               C_dMeanRevSpeed, C_fMeanRevSpeed,
                               C_dSigma, C_fSigma,
                               C_prtyCorr,  C_prtyVol,
                               C_ratesCorr, 
                               LocalGetNumObjectId(C_dZcId),
                               LocalGetNumObjectId(C_fZcId), 
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
           retCode = ARMLOCAL_IDTREEHW(C_startDate, C_horizon, (long) C_nbSteps,
                                  C_dMeanRevSpeed, C_fMeanRevSpeed,
                                  C_dSigma, C_fSigma,
                                  C_prtyCorr, C_prtyVol,
                                  C_ratesCorr, 
                                  LocalGetNumObjectId(C_dZcId),
                                  LocalGetNumObjectId(C_fZcId), 
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

            retCode = ARMLOCAL_IDTREEHW(C_startDate, C_horizon, (long) C_nbSteps,
                                   C_dMeanRevSpeed, C_fMeanRevSpeed,
                                   C_dSigma, C_fSigma,
                                   C_prtyCorr, C_prtyVol,
                                   C_ratesCorr, 
                                   LocalGetNumObjectId(C_dZcId),
                                   LocalGetNumObjectId(C_fZcId), C_result);

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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_IDTREEHW" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}



__declspec(dllexport) LPXLOPER WINAPI Local_PXL_IDTREEHW(LPXLOPER XL_startDate,
														 LPXLOPER XL_horizon,
														 LPXLOPER XL_nbSteps,
														 LPXLOPER XL_dMeanRevSpeed,
														 LPXLOPER XL_fMeanRevSpeed,
														 LPXLOPER XL_dSigma,
														 LPXLOPER XL_fSigma,
														 LPXLOPER XL_prtyCorr,
														 LPXLOPER XL_prtyVol,
														 LPXLOPER XL_ratesCorr,
														 LPXLOPER XL_dZcId,
														 LPXLOPER XL_fZcId)
{
	ADD_LOG("Local_PXL_IDTREEHW");
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();

	// C variable
	double C_startDate;
	double C_horizon;

	double C_nbSteps;

	double C_dMeanRevSpeed;
	double C_fMeanRevSpeed;

	double C_dSigma;
	double C_fSigma;

	double C_prtyCorr;
	double C_prtyVol;
	double C_ratesCorr;

	CCString C_dZcId;
	CCString C_fZcId;

	// error
	static int error;
	static char* reason = "";

	XL_readNumCell(XL_startDate,C_startDate," ARM_ERR: start date: date expected",C_result);
	XL_readNumCell(XL_horizon, C_horizon," ARM_ERR: horizon date: date expected",C_result);
	XL_readNumCell(XL_nbSteps, C_nbSteps," ARM_ERR: number of steps: numeric expected",C_result);
	XL_readNumCell(XL_dMeanRevSpeed,C_dMeanRevSpeed," ARM_ERR: domestic mean reversion: numeric expected",C_result);
	XL_readNumCell(XL_fMeanRevSpeed,C_fMeanRevSpeed," ARM_ERR: foreign mean reversion: numeric expected",C_result);
	XL_readNumCell(XL_dSigma, C_dSigma," ARM_ERR: domestic sigma: numeric expected",C_result);
	XL_readNumCell(XL_fSigma, C_fSigma," ARM_ERR: foreign sigma: numeric expected",C_result);
	XL_readNumCell(XL_prtyCorr, C_prtyCorr," ARM_ERR: correlation: numeric expected",C_result);
	XL_readNumCell(XL_prtyVol, C_prtyVol," ARM_ERR: volatility: numeric expected",C_result);
	XL_readNumCell(XL_ratesCorr, C_ratesCorr," ARM_ERR: rates correlation: numeric expected",C_result);
	XL_readStrCell(XL_dZcId, C_dZcId," ARM_ERR: domestic curve id: object expected",C_result);
	XL_readStrCell(XL_fZcId, C_fZcId," ARM_ERR: foreign curve id: object expected",C_result);

	long retCode;
	long objId;
	
	CCString curClass = LOCAL_I3DTREE_HWTREE_CLASS;
	CCString stringId;

    retCode = ARMLOCAL_IDTREEHW(C_startDate, C_horizon, (long) C_nbSteps,
                           C_dMeanRevSpeed, C_fMeanRevSpeed,
                           C_dSigma, C_fSigma,
                           C_prtyCorr,  C_prtyVol,
                           C_ratesCorr, 
                           LocalGetNumObjectId(C_dZcId),
                           LocalGetNumObjectId(C_fZcId), 
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_PXL_IDTREEHW" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}



__declspec(dllexport) LPXLOPER WINAPI Local_HWANALYTICSIGVAR(LPXLOPER XL_zcId,
															 LPXLOPER XL_a,
															 LPXLOPER XL_datesSigma)
{
	ADD_LOG("Local_HWANALYTICSIGVAR");
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();

	// C variable
	CCString C_zcId;

	double C_a;
	VECTOR<double> C_datesSigma;

	VECTOR<CCString> C_dates;
	VECTOR<double>   C_sigmas;	

	// error
	static int error;
	static char* reason = "";

    XL_readStrCell(XL_zcId, C_zcId," ARM_ERR: curve id: object expected",C_result);
    XL_readNumCell(XL_a, C_a," ARM_ERR: a: numeric expected",C_result);
    XL_readNumVector(XL_datesSigma, C_datesSigma," ARM_ERR: [date, sigma]: array expected",C_result);

	double tmp_rate = -1;
	
	if ( C_datesSigma.size () % 2 != 0 )
	{
		C_result.setMsg (" ARM_ERR: check your [Dates, Sigma(s)] array");

		ARM_ARG_ERR();

		return(LPXLOPER)&XL_result;
	}

	long real_size = C_datesSigma.size () / 2;
	long j = 0;

	char* sDate = new char[11];

	for (int i = 0; i < real_size; i++)
	{
		Local_XLDATE2ARMDATE(C_datesSigma[j],sDate);
		C_dates.push_back(sDate);
		C_sigmas.push_back(C_datesSigma[j+1]);
		j += 2;
	}
	if (sDate)
		delete [] sDate;
	sDate = NULL;

    if ( real_size == 1 )
    {
       C_dates.push_back(C_dates[0]);
       C_sigmas.push_back(C_sigmas[0]);
       
       real_size = 2;
    }

	long retCode;
	long objId;
	CCString prevClass;
	
	CCString curClass = LOCAL_HWSIGVAR_ANALYTIC_CLASS;
	CCString stringId = GetLastCurCellEnvValue();

	if (!stringId)
	{
        retCode = ARMLOCAL_HWANALYTICSIGVAR(LocalGetNumObjectId(C_zcId),
                                       C_a, real_size,
                                       C_dates, C_sigmas, C_result);

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
           retCode = ARMLOCAL_HWANALYTICSIGVAR(LocalGetNumObjectId(C_zcId), 
                                          C_a, real_size,
                                          C_dates, C_sigmas,
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

            retCode = ARMLOCAL_HWANALYTICSIGVAR(LocalGetNumObjectId(C_zcId),
                                           C_a, real_size,
                                           C_dates, C_sigmas, C_result);

	
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_HWANALYTICSIGVAR" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


__declspec(dllexport) LPXLOPER WINAPI Local_PXL_HWANALYTICSIGVAR(LPXLOPER XL_zcId,
																 LPXLOPER XL_a,
																 LPXLOPER XL_datesSigma)
{
	ADD_LOG("Local_PXL_HWANALYTICSIGVAR");
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();

	// C variable
	CCString C_zcId;

	double C_a;
	VECTOR<double> C_datesSigma;

	VECTOR<CCString> C_dates;
	VECTOR<double>   C_sigmas;	

	// error
	static int error;
	static char* reason = "";

    XL_readStrCell(XL_zcId, C_zcId," ARM_ERR: curve id: object expected",C_result);
    XL_readNumCell(XL_a, C_a," ARM_ERR: a: numeric expected",C_result);
    XL_readNumVector(XL_datesSigma, C_datesSigma," ARM_ERR: [date, sigma]: array expected",C_result);

	double tmp_rate = -1;
	
	if ( C_datesSigma.size () % 2 != 0 )
	{
		C_result.setMsg (" ARM_ERR: check your [Dates, Sigma(s)] array");

		ARM_ARG_ERR();

		return(LPXLOPER)&XL_result;
	}

	long real_size = C_datesSigma.size () / 2;
	long j = 0;

	char* sDate = new char[11];

	for (int i = 0; i < real_size; i++)
	{
		Local_XLDATE2ARMDATE(C_datesSigma[j],sDate);
		C_dates.push_back(sDate);
		C_sigmas.push_back(C_datesSigma[j+1]);
		j += 2;
	}
	if (sDate)
		delete [] sDate;
	sDate = NULL;

    if ( real_size == 1 )
    {
       C_dates.push_back(C_dates[0]);
       C_sigmas.push_back(C_sigmas[0]);
       
       real_size = 2;
    }

	long retCode;
	long objId;

	CCString curClass = LOCAL_HWSIGVAR_ANALYTIC_CLASS;
	CCString stringId;

    retCode = ARMLOCAL_HWANALYTICSIGVAR(LocalGetNumObjectId(C_zcId),
                                   C_a, real_size,
                                   C_dates, C_sigmas, C_result);

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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_PXL_HWANALYTICSIGVAR" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


__declspec(dllexport) LPXLOPER WINAPI Local_BASISMCFRM(LPXLOPER XL_zcId,
													   LPXLOPER XL_baZcId,
													   LPXLOPER XL_volId,
													   LPXLOPER XL_smileId,
													   LPXLOPER XL_productType,
													   LPXLOPER XL_horizon,
													   LPXLOPER XL_nbTraj,
													   LPXLOPER XL_MCGeneratorType,
													   LPXLOPER XL_shapeDecay,
													   LPXLOPER XL_shapeSlope,
													   LPXLOPER XL_shapeAsymptote,
													   LPXLOPER XL_nbFactor,
													   LPXLOPER XL_indexes,
													   LPXLOPER XL_correlatedIndex,
													   LPXLOPER XL_corrMatrix,
													   LPXLOPER XL_control,
													   LPXLOPER XL_seed)
{
	ADD_LOG("Local_BASISMCFRM");
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();

	// C variable
	CCString C_zcId;
	CCString C_baZcId;
	CCString C_volId;
	CCString C_smileId;
	CCString CS_productType;

	long C_productType;

	double C_horizon;
	double C_nbTraj;
	CCString CS_MCGeneratorType;
	double C_shapeDecay;
	double C_shapeSlope;
	double C_shapeAsymptote;
	double C_nbFactor;

	long   C_MCGeneratorType;

    VECTOR<double> C_indexes;

    VECTOR<double> C_correlatedIndex;	

	VECTOR<double> C_corrMatrix;

	CCString CS_control;
	long   C_control;
	double C_seed;
	double C_seedDefault = 10000.0;

	double zero = 0.0;
	double one = 1.0;
	double thousand = 1000.0;

	// error
	static int error;
	static char* reason = "";


	XL_readStrCell(XL_zcId, C_zcId, " ARM_ERR: Curve id: string expected",C_result);
	XL_readStrCell(XL_baZcId, C_baZcId, " ARM_ERR: Curve id: string expected",C_result);
	XL_readStrCell(XL_volId, C_volId, " ARM_ERR: Swaption Vol id: string expected",C_result);
	XL_readStrCell(XL_smileId, C_smileId, " ARM_ERR: Smile id: string expected",C_result);
	XL_readStrCell(XL_productType, CS_productType, " ARM_ERR: Product Type: string expected",C_result);

	if ( (C_productType = ARM_ConvAutoMode(CS_productType, C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();

		return(LPXLOPER)&XL_result;
	}

	XL_readNumCell(XL_horizon, C_horizon," ARM_ERR: Horizon date: date expected",C_result);
    XL_readNumCellWD(XL_nbTraj,C_nbTraj, thousand," ARM_ERR: Number of trajectories : numeric expected",C_result);
    XL_readStrCellWD(XL_MCGeneratorType, CS_MCGeneratorType, "DEFAULT", " ARM_ERR: Generator type. : string expected",C_result);

    if ( CS_MCGeneratorType == (CCString) "DEFAULT" )
		C_MCGeneratorType = K_MC_FAURE;
    else
	{
		C_MCGeneratorType = ARM_ConvMcMethod(CS_MCGeneratorType, C_result);
		if ( C_MCGeneratorType == ARM_DEFAULT_ERR )
		{
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}
	}

	XL_readNumCellWD(XL_shapeDecay, C_shapeDecay, zero, " ARM_ERR: shape decay. : numeric expected",C_result);
	XL_readNumCellWD(XL_shapeSlope, C_shapeSlope, zero, " ARM_ERR: shape slope. : numeric expected",C_result);
	XL_readNumCellWD(XL_shapeAsymptote, C_shapeAsymptote, zero, " ARM_ERR: shape Asymptote. : numeric expected",C_result);
	XL_readNumCellWD(XL_nbFactor, C_nbFactor, one, " ARM_ERR: Nb factor. : numeric expected",C_result);

   	if(( XL_indexes->xltype == xltypeMissing ) || ( XL_indexes->xltype == xltypeNil ))
	{
    }
	else
	{
		XL_readNumVector(XL_indexes, C_indexes, " ARM_ERR: indexes: array of numeric expected",C_result);
	}

  

   	if (( XL_correlatedIndex->xltype == xltypeMissing ) || ( XL_correlatedIndex->xltype == xltypeNil ))
	{
    }
	else
	{
		XL_readNumVector(XL_correlatedIndex, C_correlatedIndex, " ARM_ERR: Correlation correlated maturities : array of numeric expected",C_result);
    }

	if(( XL_corrMatrix->xltype == xltypeMissing ) || ( XL_corrMatrix->xltype == xltypeNil ))
	{
	}
	else
	{
		XL_readNumVector(XL_corrMatrix, C_corrMatrix, " ARM_ERR: Correlation matrix: array of numeric expected",C_result);
    }


	XL_readStrCellWD(XL_control, CS_control, "DEFAULT", " ARM_ERR: Control variate: string expected : IN or OUT",C_result);  

	if ( CS_control == (CCString) "DEFAULT" )
	{
		C_control = 0;
	}
	else
	{
		C_control = ARM_ConvInOut(CS_control, C_result);

		if ( C_control == ARM_DEFAULT_ERR )
		{
			ARM_ARG_ERR();

			return(LPXLOPER)&XL_result;
		}
	}

	XL_readNumCellWD(XL_seed, C_seed, C_seedDefault," ARM_ERR: seed : numeric expected",C_result);

	long retCode;
	long objId;
	CCString prevClass;
	
	CCString curClass = LOCAL_BMCFRM_CLASS;
	CCString stringId = GetLastCurCellEnvValue ();

	if (!stringId)
	{
        retCode = ARMLOCAL_BASISMCFRM(LocalGetNumObjectId(C_zcId),
									  LocalGetNumObjectId(C_baZcId),
									  LocalGetNumObjectId(C_volId), 
									  LocalGetNumObjectId(C_smileId),
									  C_productType,
									  C_horizon,
									  (long) C_nbTraj,
									  C_MCGeneratorType,
									  C_shapeDecay,
									  C_shapeSlope,
									  (long) C_shapeAsymptote,
									  (long) C_nbFactor,
									  C_indexes,
									  C_correlatedIndex,
									  C_corrMatrix,
									  C_control,
									  (long) C_seed,
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
			retCode = ARMLOCAL_BASISMCFRM(LocalGetNumObjectId(C_zcId),
										  LocalGetNumObjectId(C_baZcId),
										  LocalGetNumObjectId(C_volId), 
										  LocalGetNumObjectId(C_smileId),
										  C_productType,
										  C_horizon,
										  (long) C_nbTraj,
										  C_MCGeneratorType,
										  C_shapeDecay,
										  C_shapeSlope,
										  (long) C_shapeAsymptote,
										  (long) C_nbFactor,
										  C_indexes,
										  C_correlatedIndex,
										  C_corrMatrix,
										  C_control,
										  (long) C_seed, 
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

            retCode = ARMLOCAL_BASISMCFRM(LocalGetNumObjectId(C_zcId),
										  LocalGetNumObjectId(C_baZcId),
										  LocalGetNumObjectId(C_volId), 
										  LocalGetNumObjectId(C_smileId),
										  C_productType,
										  C_horizon,
										  (long) C_nbTraj,
										  C_MCGeneratorType,
										  C_shapeDecay,
										  C_shapeSlope,
										  (long) C_shapeAsymptote,
										  (long) C_nbFactor,
										  C_indexes,
										  C_correlatedIndex,
										  C_corrMatrix,
										  C_control,
										  (long) C_seed, 
										  C_result, 
										  objId);
	
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_BASISMCFRM" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}





__declspec(dllexport) LPXLOPER WINAPI Local_PXL_BASISMCFRM(LPXLOPER XL_zcId,
														   LPXLOPER XL_baZcId,
														   LPXLOPER XL_volId,
														   LPXLOPER XL_smileId,
														   LPXLOPER XL_productType,
														   LPXLOPER XL_horizon,
														   LPXLOPER XL_nbTraj,
														   LPXLOPER XL_MCGeneratorType,
														   LPXLOPER XL_shapeDecay,
														   LPXLOPER XL_shapeSlope,
														   LPXLOPER XL_shapeAsymptote,
														   LPXLOPER XL_nbFactor,
														   LPXLOPER XL_indexes,
														   LPXLOPER XL_correlatedIndex,
														   LPXLOPER XL_corrMatrix,
														   LPXLOPER XL_control,
														   LPXLOPER XL_seed)
{
	ADD_LOG("Local_PXL_BASISMCFRM");
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();

	// C variable
	CCString C_zcId;
	CCString C_baZcId;
    CCString C_volId;
    CCString C_smileId;
    CCString CS_productType;

	long C_productType;
    
	
    double C_horizon;
    double C_nbTraj;
	CCString CS_MCGeneratorType;
    double C_shapeDecay;
    double C_shapeSlope;
    double C_shapeAsymptote;
    double C_nbFactor;

	long   C_MCGeneratorType;
       
    VECTOR<double> C_indexes;
   
    VECTOR<double> C_correlatedIndex;	

	VECTOR<double> C_corrMatrix;

    CCString CS_control;
	long   C_control;
	double C_seed;
    double C_seedDefault = 10000.0;
      
    double zero = 0.0;
    double one = 1.0;
    double thousand = 1000.0;

	// error
	static int error;
	static char* reason = "";

    XL_readStrCell(XL_zcId, C_zcId, " ARM_ERR: Curve id: string expected",C_result);
	XL_readStrCell(XL_baZcId, C_baZcId, " ARM_ERR: Curve id: string expected",C_result);
    XL_readStrCell(XL_volId, C_volId, " ARM_ERR: Swaption Vol id: string expected",C_result);
    XL_readStrCell(XL_smileId, C_smileId, " ARM_ERR: Smile id: string expected",C_result);
	XL_readStrCell(XL_productType, CS_productType, " ARM_ERR: Product Type : string expected",C_result);

	C_productType = ARM_ConvAutoMode(CS_productType, C_result);

	if ( C_productType == ARM_DEFAULT_ERR )
	{
		 ARM_ARG_ERR();

		 return(LPXLOPER)&XL_result;
	}

	XL_readNumCell(XL_horizon, C_horizon," ARM_ERR: Horizon date: date expected",C_result);
    XL_readNumCellWD(XL_nbTraj,C_nbTraj, thousand," ARM_ERR: Number of trajectories : numeric expected",C_result);
    XL_readStrCellWD(XL_MCGeneratorType, CS_MCGeneratorType, "DEFAULT", " ARM_ERR: Generator type. : string expected",C_result);

    if ( CS_MCGeneratorType == (CCString) "DEFAULT" )
	{
	   C_MCGeneratorType = K_MC_FAURE;
	}
    else
	{
       C_MCGeneratorType = ARM_ConvMcMethod(CS_MCGeneratorType, C_result);

	   if ( C_MCGeneratorType == ARM_DEFAULT_ERR )
	   {
		  ARM_ARG_ERR();

		  return (LPXLOPER)&XL_result;
	   }
	}

    XL_readNumCellWD(XL_shapeDecay, C_shapeDecay, zero, " ARM_ERR: shape decay. : numeric expected",C_result);
    XL_readNumCellWD(XL_shapeSlope, C_shapeSlope, zero, " ARM_ERR: shape slope. : numeric expected",C_result);
    XL_readNumCellWD(XL_shapeAsymptote, C_shapeAsymptote, zero, " ARM_ERR: shape Asymptote. : numeric expected",C_result);
    XL_readNumCellWD(XL_nbFactor, C_nbFactor, one, " ARM_ERR: Nb factor. : numeric expected",C_result);

   	if(( XL_indexes->xltype == xltypeMissing ) || ( XL_indexes->xltype == xltypeNil ))
	{
    }
	else
	{
		XL_readNumVector(XL_indexes, C_indexes, " ARM_ERR: indexes: array of numeric expected",C_result);
    }

   	if (( XL_correlatedIndex->xltype == xltypeMissing ) 
		|| 
		( XL_correlatedIndex->xltype == xltypeNil ) 
	   )
	{
    }
	else
	{
       XL_readNumVector(XL_correlatedIndex, C_correlatedIndex, " ARM_ERR: Correlation correlated maturities : array of numeric expected",C_result);
    }

	if(( XL_corrMatrix->xltype == xltypeMissing ) 
	   || 
	   ( XL_corrMatrix->xltype == xltypeNil )
	  )
	{
    }
	else
	{
        XL_readNumVector(XL_corrMatrix,C_corrMatrix, " ARM_ERR: Correlation matrix: array of numeric expected",C_result);
    }


    XL_readStrCellWD(XL_control, CS_control, "DEFAULT", " ARM_ERR: Control variate: string expected : IN or OUT",C_result);  

	if ( CS_control == (CCString) "DEFAULT" )
	{
	   C_control = 0;
	}
	else
	{
	   C_control = ARM_ConvInOut(CS_control, C_result);

	   if ( C_control == ARM_DEFAULT_ERR )
	   {
		  ARM_ARG_ERR();

		  return(LPXLOPER)&XL_result;
	   }
	}

	XL_readNumCellWD(XL_seed, C_seed, C_seedDefault," ARM_ERR: seed : numeric expected",C_result);

	long retCode;
	long objId;
	CCString prevClass;
	
	CCString curClass = LOCAL_BMCFRM_CLASS;
	CCString stringId;

    retCode = ARMLOCAL_BASISMCFRM(LocalGetNumObjectId(C_zcId),
								  LocalGetNumObjectId(C_baZcId),
								  LocalGetNumObjectId(C_volId), 
								  LocalGetNumObjectId(C_smileId),
								  C_productType,
								  C_horizon,
								  (long) C_nbTraj,
								  C_MCGeneratorType,
								  C_shapeDecay,
								  C_shapeSlope,
								  (long) C_shapeAsymptote,
								  (long) C_nbFactor,
								  C_indexes,
								  C_correlatedIndex,
								  C_corrMatrix,
								  C_control,
								  (long) C_seed,
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
		PXL_ARM_ERR();
	}

//	ARM_END();
	
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_PXL_BASISMCFRM" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}



__declspec(dllexport) LPXLOPER WINAPI Local_BASISMCFRM2CR(LPXLOPER XL_zcId,
														  LPXLOPER XL_baZcId,
														  LPXLOPER XL_swoptvolId,
														  LPXLOPER XL_swoptsmileId,
														  LPXLOPER XL_irgvolId,
														  LPXLOPER XL_irgsmileId,
														  LPXLOPER XL_productType,
														  LPXLOPER XL_horizon,
														  LPXLOPER XL_nbTraj,
														  LPXLOPER XL_MCGeneratorType,
														  LPXLOPER XL_shapeDecay,
														  LPXLOPER XL_shapeSlope,
														  LPXLOPER XL_shapeAsymptote,
														  LPXLOPER XL_nbFactor,
														  LPXLOPER XL_indexes,
														  LPXLOPER XL_correlatedIndex,
														  LPXLOPER XL_corrMatrix,
														  LPXLOPER XL_control,
														  LPXLOPER XL_seed)
{
	ADD_LOG("Local_BASISMCFRM2CR");
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();

	// C variable
	CCString C_zcId;
	CCString C_baZcId;
    CCString C_swoptvolId;
    CCString C_swoptsmileId;
    CCString C_irgvolId;
    CCString C_irgsmileId;
    CCString CS_productType;

	long C_productType;
    
	
    double C_horizon;
    double C_nbTraj;
	CCString CS_MCGeneratorType;
    double C_shapeDecay;
    double C_shapeSlope;
    double C_shapeAsymptote;
    double C_nbFactor;

	long   C_MCGeneratorType;
    
   
    VECTOR<double> C_indexes;
   
    VECTOR<double> C_correlatedIndex;	

	VECTOR<double> C_corrMatrix;

    CCString CS_control;
	long   C_control;
	double C_seed;
    double C_seedDefault = 10000.0;
      
    double zero = 0.0;
    double one = 1.0;
    double thousand = 1000.0;
	
	// error
	static int error;
	static char* reason = "";


    XL_readStrCell(XL_zcId, C_zcId, " ARM_ERR: Curve id: string expected",C_result);
	XL_readStrCell(XL_baZcId, C_baZcId, " ARM_ERR: Curve id: string expected",C_result);
    XL_readStrCell(XL_swoptvolId, C_swoptvolId, " ARM_ERR: Swopt Vol id: string expected",C_result);
    XL_readStrCell(XL_swoptsmileId, C_swoptsmileId, " ARM_ERR: Swopt Smile id: string expected",C_result);
    XL_readStrCell(XL_irgvolId, C_irgvolId, " ARM_ERR: IRG Vol id: string expected",C_result);
    XL_readStrCell(XL_irgsmileId, C_irgsmileId, " ARM_ERR: IRG Smile id: string expected",C_result);
	XL_readStrCell(XL_productType, CS_productType, " ARM_ERR: Product Type: string expected",C_result);

	C_productType = ARM_ConvAutoMode(CS_productType, C_result);

	if ( C_productType == ARM_DEFAULT_ERR )
	{
		 ARM_ARG_ERR();

		 return(LPXLOPER)&XL_result;
	}

	XL_readNumCell(XL_horizon, C_horizon," ARM_ERR: Horizon date: date expected",C_result);
    XL_readNumCellWD(XL_nbTraj,C_nbTraj, thousand," ARM_ERR: Number of trajectories : numeric expected",C_result);
    XL_readStrCellWD(XL_MCGeneratorType, CS_MCGeneratorType, "DEFAULT", " ARM_ERR: Generator type. : string expected",C_result);

    if ( CS_MCGeneratorType == (CCString) "DEFAULT" )
	{
	   C_MCGeneratorType = K_MC_FAURE;
	}
    else
	{
       C_MCGeneratorType = ARM_ConvMcMethod(CS_MCGeneratorType, C_result);

	   if ( C_MCGeneratorType == ARM_DEFAULT_ERR )
	   {
		  ARM_ARG_ERR();

		  return (LPXLOPER)&XL_result;
	   }
	}

    XL_readNumCellWD(XL_shapeDecay, C_shapeDecay, zero, " ARM_ERR: shape decay. : numeric expected",C_result);
    XL_readNumCellWD(XL_shapeSlope, C_shapeSlope, zero, " ARM_ERR: shape slope. : numeric expected",C_result);
    XL_readNumCellWD(XL_shapeAsymptote, C_shapeAsymptote, zero, " ARM_ERR: shape Asymptote. : numeric expected",C_result);
    XL_readNumCellWD(XL_nbFactor, C_nbFactor, one, " ARM_ERR: Nb factor. : numeric expected",C_result);

   	if(( XL_indexes->xltype == xltypeMissing ) || ( XL_indexes->xltype == xltypeNil ))
	{
    }
	else
	{
       XL_readNumVector(XL_indexes, C_indexes, " ARM_ERR: indexes: array of numeric expected",C_result);
    }

  

   	if (( XL_correlatedIndex->xltype == xltypeMissing ) 
		|| 
		( XL_correlatedIndex->xltype == xltypeNil ) 
	   )
	{
    }
	else
	{
       XL_readNumVector(XL_correlatedIndex, C_correlatedIndex, " ARM_ERR: Correlation correlated maturities : array of numeric expected",C_result);
    }

	if(( XL_corrMatrix->xltype == xltypeMissing ) 
	   || 
	   ( XL_corrMatrix->xltype == xltypeNil )
	  )
	{
    }
	else
	{
        XL_readNumVector(XL_corrMatrix,C_corrMatrix, " ARM_ERR: Correlation matrix: array of numeric expected",C_result);
    }


    XL_readStrCellWD(XL_control, CS_control, "DEFAULT", " ARM_ERR: Control variate: string expected : IN or OUT",C_result);  

	if ( CS_control == (CCString) "DEFAULT" )
	{
	   C_control = 0;
	}
	else
	{
	   C_control = ARM_ConvInOut(CS_control, C_result);

	   if ( C_control == ARM_DEFAULT_ERR )
	   {
		  ARM_ARG_ERR();

		  return(LPXLOPER)&XL_result;
	   }
	}

	XL_readNumCellWD(XL_seed, C_seed, C_seedDefault," ARM_ERR: seed : numeric expected",C_result);

	long retCode;
	long objId;
	CCString prevClass;
	
	CCString curClass = LOCAL_BMCFRM_CLASS;
	CCString stringId = GetLastCurCellEnvValue ();

	if (!stringId)
	{
        retCode = ARMLOCAL_BASISMCFRM2CR(LocalGetNumObjectId(C_zcId),
										 LocalGetNumObjectId(C_baZcId),
										 LocalGetNumObjectId(C_swoptvolId),
										 LocalGetNumObjectId(C_swoptsmileId),
										 LocalGetNumObjectId(C_irgvolId),
										 LocalGetNumObjectId(C_irgsmileId),
										 C_productType,
										 C_horizon,
										 (long) C_nbTraj,
										 C_MCGeneratorType,
										 C_shapeDecay,
										 C_shapeSlope,
										 (long) C_shapeAsymptote,
										 (long) C_nbFactor,
										 C_indexes,
										 C_correlatedIndex,
										 C_corrMatrix,
										 C_control,
										 (long) C_seed,
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
			retCode = ARMLOCAL_BASISMCFRM2CR(LocalGetNumObjectId(C_zcId),
											 LocalGetNumObjectId(C_baZcId),
											 LocalGetNumObjectId(C_swoptvolId),
											 LocalGetNumObjectId(C_swoptsmileId),
											 LocalGetNumObjectId(C_irgvolId),
											 LocalGetNumObjectId(C_irgsmileId),
											 C_productType,
											 C_horizon,
											 (long) C_nbTraj,
											 C_MCGeneratorType,
											 C_shapeDecay,
											 C_shapeSlope,
											 (long) C_shapeAsymptote,
											 (long) C_nbFactor,
											 C_indexes,
											 C_correlatedIndex,
											 C_corrMatrix,
											 C_control,
											 (long) C_seed,
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

            retCode = ARMLOCAL_BASISMCFRM2CR(LocalGetNumObjectId(C_zcId),
											 LocalGetNumObjectId(C_baZcId),
											 LocalGetNumObjectId(C_swoptvolId),
											 LocalGetNumObjectId(C_swoptsmileId),
											 LocalGetNumObjectId(C_irgvolId),
											 LocalGetNumObjectId(C_irgsmileId),
											 C_productType,
											 C_horizon,
											 (long) C_nbTraj,
											 C_MCGeneratorType,
											 C_shapeDecay,
											 C_shapeSlope,
											 (long) C_shapeAsymptote,
											 (long) C_nbFactor,
											 C_indexes,
											 C_correlatedIndex,
											 C_corrMatrix,
											 C_control,
											 (long) C_seed,
											 C_result,
											 objId);
	
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_BASISMCFRM2CR" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


__declspec(dllexport) LPXLOPER WINAPI Local_PXL_BASISMCFRM2CR(LPXLOPER XL_zcId,
															  LPXLOPER XL_baZcId,
															  LPXLOPER XL_swoptvolId,
															  LPXLOPER XL_swoptsmileId,
															  LPXLOPER XL_irgvolId,
															  LPXLOPER XL_irgsmileId,
															  LPXLOPER XL_productType,
															  LPXLOPER XL_horizon,
															  LPXLOPER XL_nbTraj,
															  LPXLOPER XL_MCGeneratorType,
															  LPXLOPER XL_shapeDecay,
															  LPXLOPER XL_shapeSlope,
															  LPXLOPER XL_shapeAsymptote,
															  LPXLOPER XL_nbFactor,
															  LPXLOPER XL_indexes,
															  LPXLOPER XL_correlatedIndex,
															  LPXLOPER XL_corrMatrix,
															  LPXLOPER XL_control,
															  LPXLOPER XL_seed)
{
	ADD_LOG("Local_PXL_BASISMCFRM2CR");
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();

	// C variable
	CCString C_zcId;
	CCString C_baZcId;
    CCString C_swoptvolId;
    CCString C_swoptsmileId;
    CCString C_irgvolId;
    CCString C_irgsmileId;
    CCString CS_productType;

	long C_productType;

    double C_horizon;
    double C_nbTraj;
	CCString CS_MCGeneratorType;
    double C_shapeDecay;
    double C_shapeSlope;
    double C_shapeAsymptote;
    double C_nbFactor;

	long   C_MCGeneratorType;
       
    VECTOR<double> C_indexes;
   
    VECTOR<double> C_correlatedIndex;	

	VECTOR<double> C_corrMatrix;

    CCString CS_control;
	long   C_control;
	double C_seed;
    double C_seedDefault = 10000.0;
    
  
    double zero = 0.0;
    double one = 1.0;
    double thousand = 1000.0;
	
	// error
	static int error;
	static char* reason = "";

    XL_readStrCell(XL_zcId, C_zcId, " ARM_ERR: Curve id: string expected",C_result);
	XL_readStrCell(XL_baZcId, C_baZcId, " ARM_ERR: Curve id: string expected",C_result);
    XL_readStrCell(XL_swoptvolId, C_swoptvolId, " ARM_ERR: Swopt Vol id: string expected",C_result);
    XL_readStrCell(XL_swoptsmileId, C_swoptsmileId, " ARM_ERR: Swopt Smile id: string expected",C_result);
    XL_readStrCell(XL_irgvolId, C_irgvolId, " ARM_ERR: IRG Vol id: string expected",C_result);
    XL_readStrCell(XL_irgsmileId, C_irgsmileId, " ARM_ERR: IRG Smile id: string expected",C_result);
	XL_readStrCell(XL_productType, CS_productType, " ARM_ERR: Product Type: string expected",C_result);

	C_productType = ARM_ConvAutoMode(CS_productType, C_result);

	if ( C_productType == ARM_DEFAULT_ERR )
	{
		 ARM_ARG_ERR();

		 return(LPXLOPER)&XL_result;
	}

	XL_readNumCell(XL_horizon, C_horizon," ARM_ERR: Horizon date: date expected",C_result);
    XL_readNumCellWD(XL_nbTraj,C_nbTraj, thousand," ARM_ERR: Number of trajectories : numeric expected",C_result);
    XL_readStrCellWD(XL_MCGeneratorType, CS_MCGeneratorType, "DEFAULT", " ARM_ERR: Generator type. : string expected",C_result);

    if ( CS_MCGeneratorType == (CCString) "DEFAULT" )
	{
	   C_MCGeneratorType = K_MC_FAURE;
	}
    else
	{
       C_MCGeneratorType = ARM_ConvMcMethod(CS_MCGeneratorType, C_result);

	   if ( C_MCGeneratorType == ARM_DEFAULT_ERR )
	   {
		  ARM_ARG_ERR();

		  return (LPXLOPER)&XL_result;
	   }
	}

    XL_readNumCellWD(XL_shapeDecay, C_shapeDecay, zero, " ARM_ERR: shape decay. : numeric expected",C_result);
    XL_readNumCellWD(XL_shapeSlope, C_shapeSlope, zero, " ARM_ERR: shape slope. : numeric expected",C_result);
    XL_readNumCellWD(XL_shapeAsymptote, C_shapeAsymptote, zero, " ARM_ERR: shape Asymptote. : numeric expected",C_result);
    XL_readNumCellWD(XL_nbFactor, C_nbFactor, one, " ARM_ERR: Nb factor. : numeric expected",C_result);

   	if(( XL_indexes->xltype == xltypeMissing ) || ( XL_indexes->xltype == xltypeNil ))
	{
    }
	else
	{
       XL_readNumVector(XL_indexes, C_indexes, " ARM_ERR: indexes: array of numeric expected",C_result);
    }

  

   	if (( XL_correlatedIndex->xltype == xltypeMissing ) 
		|| 
		( XL_correlatedIndex->xltype == xltypeNil ) 
	   )
	{
    }
	else
	{
       XL_readNumVector(XL_correlatedIndex, C_correlatedIndex, " ARM_ERR: Correlation correlated maturities : array of numeric expected",C_result);
    }

	if(( XL_corrMatrix->xltype == xltypeMissing ) 
	   || 
	   ( XL_corrMatrix->xltype == xltypeNil )
	  )
	{
    }
	else
	{
        XL_readNumVector(XL_corrMatrix,C_corrMatrix, " ARM_ERR: Correlation matrix: array of numeric expected",C_result);
    }


    XL_readStrCellWD(XL_control, CS_control, "DEFAULT", " ARM_ERR: Control variate: string expected : IN or OUT",C_result);  

	if ( CS_control == (CCString) "DEFAULT" )
	{
	   C_control = 0;
	}
	else
	{
	   C_control = ARM_ConvInOut(CS_control, C_result);

	   if ( C_control == ARM_DEFAULT_ERR )
	   {
		  ARM_ARG_ERR();

		  return(LPXLOPER)&XL_result;
	   }
	}

	XL_readNumCellWD(XL_seed, C_seed, C_seedDefault," ARM_ERR: seed : numeric expected",C_result);

	long retCode;
	long objId;
	CCString prevClass;
	
	CCString curClass = LOCAL_BMCFRM_CLASS;
	CCString stringId;

    retCode = ARMLOCAL_BASISMCFRM2CR(LocalGetNumObjectId(C_zcId),
									 LocalGetNumObjectId(C_baZcId),
									 LocalGetNumObjectId(C_swoptvolId),
									 LocalGetNumObjectId(C_swoptsmileId),
									 LocalGetNumObjectId(C_irgvolId),
									 LocalGetNumObjectId(C_irgsmileId),
									 C_productType,
									 C_horizon,
									 (long) C_nbTraj,
									 C_MCGeneratorType,
									 C_shapeDecay,
									 C_shapeSlope,
									 (long) C_shapeAsymptote,
									 (long) C_nbFactor,
									 C_indexes,
									 C_correlatedIndex,
									 C_corrMatrix,
									 C_control,
									 (long) C_seed,
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
		PXL_ARM_ERR();
	}

//	ARM_END();
	
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_PXL_BASISMCFRM2CR" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


__declspec(dllexport) LPXLOPER WINAPI Local_GetCalibrationOutput(LPXLOPER XL_modId)
{
	ADD_LOG("Local_GetCalibrationOutput");
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
	CCString C_modId;
	
	long pfSize, hwVolSize;

	long i;

	// error
	static int error;
	static char* reason = "";

    XL_readStrCell(XL_modId, C_modId, " ARM_ERR: model id: object expected",C_result);

	long retCode = ARMLOCAL_GetCalibrationOutputPFSize(LocalGetNumObjectId (C_modId), C_result);

	if ( retCode == ARM_OK )
	{
//		FreeCurCellErr ();
//		XL_result.xltype = xltypeNum;
//		XL_result.val.num = C_result.getDouble ();
		if ( (C_result.getDouble () > 0) && (C_result.getDouble () < 50) )
			pfSize = C_result.getDouble ();
		else
		{
			ARM_ERR();
			return (LPXLOPER)&XL_result;
		}
	}
	else
	{
		ARM_ERR();
		return (LPXLOPER)&XL_result;
	}

	retCode = ARMLOCAL_GetCalibrationOutputHWVolSize(LocalGetNumObjectId (C_modId), C_result);

	if ( retCode == ARM_OK )
	{
//		FreeCurCellErr ();
//		XL_result.xltype = xltypeNum;
		hwVolSize = C_result.getDouble ();
	}
	else
	{
		ARM_ERR();
		return (LPXLOPER)&XL_result;
	}

	int nbrows = hwVolSize + pfSize + 1;
	int nbcolumns = 10;

	FreeCurCellErr ();
	XL_result.xltype = xltypeMulti;
	XL_result.val.array.columns = nbcolumns;
	XL_result.val.array.rows = nbrows; 
	XL_result.val.array.lparray = pxArray = (LPXLOPER)GlobalAlloc (GMEM_ZEROINIT, nbrows * nbcolumns * sizeof (XLOPER));

	for (i=0; i<nbrows; i++)
	{
		for (int j=0; j<nbcolumns; j++)
		{
			pxArray[XL_Coordonnate2Rank (i, j, nbcolumns)].xltype = xltypeStr;
			pxArray[XL_Coordonnate2Rank (i, j, nbcolumns)].val.str = XL_StrC2StrPascal ("-");
			pxArray[XL_Coordonnate2Rank (i, j, nbcolumns)].xltype |= xlbitDLLFree;
		}
	}

	retCode = ARMLOCAL_GetCalibrationOutputMeanRev(LocalGetNumObjectId (C_modId), C_result);
	pxArray[XL_Coordonnate2Rank (0, 0, nbcolumns)].xltype = xltypeNum;
	pxArray[XL_Coordonnate2Rank (0, 0, nbcolumns)].val.num = C_result.getDouble ();

	pxArray[XL_Coordonnate2Rank (1, 0, nbcolumns)].xltype = xltypeNum;
	pxArray[XL_Coordonnate2Rank (1, 0, nbcolumns)].val.num = hwVolSize;

	pxArray[XL_Coordonnate2Rank (2, 0, nbcolumns)].xltype = xltypeNum;
	pxArray[XL_Coordonnate2Rank (2, 0, nbcolumns)].val.num = pfSize;

/*	for (i = 3; i <= hwVolSize; i++)
	{
		pxArray[XL_Coordonnate2Rank (i, 0, nbcolumns)].xltype = xltypeStr;
		pxArray[XL_Coordonnate2Rank (i, 0, nbcolumns)].val.str = XL_StrC2StrPascal ("-");
		pxArray[XL_Coordonnate2Rank (i, 0, nbcolumns)].xltype |= xlbitDLLFree;
	}
*/
	for (i = 0;i < hwVolSize; i++)
	{
		retCode = ARMLOCAL_GetCalibrationOutputDate(LocalGetNumObjectId (C_modId), i, C_result);
		pxArray[XL_Coordonnate2Rank (i, 1, nbcolumns)].xltype = xltypeNum;
		pxArray[XL_Coordonnate2Rank (i, 1, nbcolumns)].val.num = Local_ARMDATE2XLDATE(C_result.getString ());
	}
/*
	pxArray[XL_Coordonnate2Rank (i, 1, nbcolumns)].xltype = xltypeStr;
	pxArray[XL_Coordonnate2Rank (i, 1, nbcolumns)].val.str = XL_StrC2StrPascal ("-");
	pxArray[XL_Coordonnate2Rank (i, 1, nbcolumns)].xltype |= xlbitDLLFree;
*/
	for (i = 0;i < hwVolSize; i++)
	{
		retCode = ARMLOCAL_GetCalibrationVolSched(LocalGetNumObjectId (C_modId), i, C_result);
		pxArray[XL_Coordonnate2Rank (i, 2, nbcolumns)].xltype = xltypeNum;
		pxArray[XL_Coordonnate2Rank (i, 2, nbcolumns)].val.num = C_result.getDouble ();
	}

/*	pxArray[XL_Coordonnate2Rank (i, 2, nbcolumns)].xltype = xltypeStr;
	pxArray[XL_Coordonnate2Rank (i, 2, nbcolumns)].val.str = XL_StrC2StrPascal ("-");
	pxArray[XL_Coordonnate2Rank (i, 2, nbcolumns)].xltype |= xlbitDLLFree;

	for (int j = 3;j <= 10;j++)
	{
		for (i = 0; i <= hwVolSize;i++)
		{
			pxArray[XL_Coordonnate2Rank (i, j, nbcolumns)].xltype = xltypeStr;
			pxArray[XL_Coordonnate2Rank (i, j, nbcolumns)].val.str = XL_StrC2StrPascal ("-");
			pxArray[XL_Coordonnate2Rank (i, j, nbcolumns)].xltype |= xlbitDLLFree;
		}
	}
*/
	for (i = 0;i < pfSize; i++)
	{
		retCode = ARMLOCAL_GetCalibrationOutputSecExeDate(LocalGetNumObjectId (C_modId), i, C_result);
		pxArray[XL_Coordonnate2Rank (i+hwVolSize+1, 0, nbcolumns)].xltype = xltypeNum;
		pxArray[XL_Coordonnate2Rank (i+hwVolSize+1, 0, nbcolumns)].val.num = Local_ARMDATE2XLDATE(C_result.getString ());
	}

	for (i = 0;i < pfSize; i++)
	{
		retCode = ARMLOCAL_GetCalibrationOutputSecStartDate(LocalGetNumObjectId (C_modId), i, C_result);
		pxArray[XL_Coordonnate2Rank (i+hwVolSize+1, 1, nbcolumns)].xltype = xltypeNum;
		pxArray[XL_Coordonnate2Rank (i+hwVolSize+1, 1, nbcolumns)].val.num = Local_ARMDATE2XLDATE(C_result.getString ());
	}

	for (i = 0;i < pfSize; i++)
	{
		retCode = ARMLOCAL_GetCalibrationOutputSecEndDate(LocalGetNumObjectId (C_modId), i, C_result);
		pxArray[XL_Coordonnate2Rank (i+hwVolSize+1, 2, nbcolumns)].xltype = xltypeNum;
		pxArray[XL_Coordonnate2Rank (i+hwVolSize+1, 2, nbcolumns)].val.num = Local_ARMDATE2XLDATE(C_result.getString ());
	}

	for (i = 0;i < pfSize; i++)
	{
		retCode = ARMLOCAL_GetCalibrationOutputSecStrike(LocalGetNumObjectId (C_modId), i, C_result);
		pxArray[XL_Coordonnate2Rank (i+hwVolSize+1, 3, nbcolumns)].xltype = xltypeNum;
		pxArray[XL_Coordonnate2Rank (i+hwVolSize+1, 3, nbcolumns)].val.num = C_result.getDouble ();
	}

	for (i = 0;i < pfSize; i++)
	{
		retCode = ARMLOCAL_GetCalibrationOutputInputVol(LocalGetNumObjectId (C_modId), i, C_result);
		pxArray[XL_Coordonnate2Rank (i+hwVolSize+1, 4, nbcolumns)].xltype = xltypeNum;
		pxArray[XL_Coordonnate2Rank (i+hwVolSize+1, 4, nbcolumns)].val.num = C_result.getDouble ();
	}

	for (i = 0;i < pfSize; i++)
	{
		retCode = ARMLOCAL_GetCalibrationOutputOutputVol(LocalGetNumObjectId (C_modId), i, C_result);
		pxArray[XL_Coordonnate2Rank (i+hwVolSize+1, 5, nbcolumns)].xltype = xltypeNum;
		pxArray[XL_Coordonnate2Rank (i+hwVolSize+1, 5, nbcolumns)].val.num = C_result.getDouble ();
	}

	for (i = 0;i < pfSize; i++)
	{
		retCode = ARMLOCAL_GetCalibrationOutputErrorVol(LocalGetNumObjectId (C_modId), i, C_result);
		pxArray[XL_Coordonnate2Rank (i+hwVolSize+1, 6, nbcolumns)].xltype = xltypeNum;
		pxArray[XL_Coordonnate2Rank (i+hwVolSize+1, 6, nbcolumns)].val.num = C_result.getDouble ();
	}

	for (i = 0;i < pfSize; i++)
	{
		retCode = ARMLOCAL_GetCalibrationOutputInputPrice(LocalGetNumObjectId (C_modId), i, C_result);
		pxArray[XL_Coordonnate2Rank (i+hwVolSize+1, 7, nbcolumns)].xltype = xltypeNum;
		pxArray[XL_Coordonnate2Rank (i+hwVolSize+1, 7, nbcolumns)].val.num = C_result.getDouble ();
	}

	for (i = 0;i < pfSize; i++)
	{
		retCode = ARMLOCAL_GetCalibrationOutputOutputPrice(LocalGetNumObjectId (C_modId), i, C_result);
		pxArray[XL_Coordonnate2Rank (i+hwVolSize+1, 8, nbcolumns)].xltype = xltypeNum;
		pxArray[XL_Coordonnate2Rank (i+hwVolSize+1, 8, nbcolumns)].val.num = C_result.getDouble ();
	}

	for (i = 0;i < pfSize; i++)
	{
		retCode = ARMLOCAL_GetCalibrationOutputErrorPrice(LocalGetNumObjectId (C_modId), i, C_result);
		pxArray[XL_Coordonnate2Rank (i+hwVolSize+1, 9, nbcolumns)].xltype = xltypeNum;
		pxArray[XL_Coordonnate2Rank (i+hwVolSize+1, 9, nbcolumns)].val.num = C_result.getDouble ();
	}

	//	ARM_END();
	
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_GetCalibrationOutput" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


__declspec(dllexport) LPXLOPER WINAPI Local_SMILEDLDCANA(LPXLOPER XL_anaModId,
														 LPXLOPER XL_dVolSpotUp,
														 LPXLOPER XL_Asymp_RowUp,
														 LPXLOPER XL_Asymp_ColUp,
														 LPXLOPER XL_pUp,
														 LPXLOPER XL_dVolSpotDo,
														 LPXLOPER XL_Asymp_RowDo,
														 LPXLOPER XL_Asymp_ColDo,
														 LPXLOPER XL_pDo)
{
	ADD_LOG("Local_SMILEDLDCANA");
	//ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();

	// C variable
	CCString C_anaModId;
   
    double C_dVolSpotUp;
	double C_Asymp_RowUp;
	double C_Asymp_ColUp;
	double C_pUp;
	double C_dVolSpotDo;
	double C_Asymp_RowDo;
	double C_Asymp_ColDo;
	double C_pDo;
	
    double C_dVolSpotUp_default = 0.;
	double C_Asymp_RowUp_default = 0.;
	double C_Asymp_ColUp_default = 0.;
	double C_pUp_default = 0.;
	double C_dVolSpotDo_default = 0.;
	double C_Asymp_RowDo_default = 0.;
	double C_Asymp_ColDo_default = 0.;
	double C_pDo_default = 0.;

	// error
	static int error;
	static char* reason = "";

	XL_readStrCell(XL_anaModId, C_anaModId, " ARM_ERR: Model id: string expected",C_result);
	XL_readNumCellWD(XL_dVolSpotUp,C_dVolSpotUp,C_dVolSpotUp_default, " ARM_ERR: Vol Spot Up: numeric expected",C_result);
	XL_readNumCellWD(XL_Asymp_RowUp,C_Asymp_RowUp,C_Asymp_RowUp_default, " ARM_ERR: Asymp_RowUp: numeric expected",C_result);
	XL_readNumCellWD(XL_Asymp_ColUp,C_Asymp_ColUp,C_Asymp_ColUp_default, " ARM_ERR: Asymp_ColUp: numeric expected",C_result);
	XL_readNumCellWD(XL_pUp,C_pUp,C_pUp_default, " ARM_ERR: Probability Up: numeric expected",C_result);
	XL_readNumCellWD(XL_dVolSpotDo,C_dVolSpotDo,C_dVolSpotDo_default, " ARM_ERR: Vol Spot Down: numeric expected",C_result);
	XL_readNumCellWD(XL_Asymp_RowDo,C_Asymp_RowDo,C_Asymp_RowDo_default, " ARM_ERR: Asymp_RowDo: numeric expected",C_result);
	XL_readNumCellWD(XL_Asymp_ColDo,C_Asymp_ColDo,C_Asymp_ColDo_default, " ARM_ERR: Asymp_ColDo: numeric expected",C_result);
	XL_readNumCellWD(XL_pDo,C_pDo,C_pDo_default, " ARM_ERR: Probability Down: numeric expected",C_result);

	if (
		(C_dVolSpotUp < 0)		||
		(C_dVolSpotDo > 0)		||
		(C_pUp < 0)				||
		(C_pDo < 0)				||
		((C_pDo + C_pUp) > 1)
		)
	{
		C_result.setMsg ("ARM_ERR: check your parameters");
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}


	long retCode;
	long objId;
	CCString prevClass;
	
	CCString curClass = LOCAL_LOGDECANA_CLASS;
	CCString stringId = GetLastCurCellEnvValue ();

	if (!stringId)
	{
        retCode = ARMLOCAL_SMILEDLDCANA(LocalGetNumObjectId (C_anaModId),
										C_dVolSpotUp,
										C_Asymp_RowUp,
										C_Asymp_ColUp,
										C_pUp,
										C_dVolSpotDo,
										C_Asymp_RowDo,
										C_Asymp_ColDo,
										C_pDo,
										C_result);

		if ( retCode == ARM_OK )
		{
			objId = C_result.getLong();
			
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
		    retCode = ARMLOCAL_SMILEDLDCANA(LocalGetNumObjectId (C_anaModId),
											C_dVolSpotUp,
											C_Asymp_RowUp,
											C_Asymp_ColUp,
											C_pUp,
											C_dVolSpotDo,
											C_Asymp_RowDo,
											C_Asymp_ColDo,
											C_pDo,
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

	        retCode = ARMLOCAL_SMILEDLDCANA(LocalGetNumObjectId (C_anaModId),
											C_dVolSpotUp,
											C_Asymp_RowUp,
											C_Asymp_ColUp,
											C_pUp,
											C_dVolSpotDo,
											C_Asymp_RowDo,
											C_Asymp_ColDo,
											C_pDo,
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

	//ARM_END();
	
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_SMILEDLDCANA" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


__declspec(dllexport) LPXLOPER WINAPI Local_PXL_SMILEDLDCANA(LPXLOPER XL_anaModId,
															 LPXLOPER XL_dVolSpotUp,
															 LPXLOPER XL_Asymp_RowUp,
															 LPXLOPER XL_Asymp_ColUp,
															 LPXLOPER XL_pUp,
															 LPXLOPER XL_dVolSpotDo,
															 LPXLOPER XL_Asymp_RowDo,
															 LPXLOPER XL_Asymp_ColDo,
															 LPXLOPER XL_pDo)
{
	ADD_LOG("Local_PXL_SMILEDLDCANA");
	//ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();

	// C variable
	CCString C_anaModId;
   
    double C_dVolSpotUp;
	double C_Asymp_RowUp;
	double C_Asymp_ColUp;
	double C_pUp;
	double C_dVolSpotDo;
	double C_Asymp_RowDo;
	double C_Asymp_ColDo;
	double C_pDo;
	
    double C_dVolSpotUp_default = 0.;
	double C_Asymp_RowUp_default = 0.;
	double C_Asymp_ColUp_default = 0.;
	double C_pUp_default = 0.;
	double C_dVolSpotDo_default = 0.;
	double C_Asymp_RowDo_default = 0.;
	double C_Asymp_ColDo_default = 0.;
	double C_pDo_default = 0.;
	
	// error
	static int error;
	static char* reason = "";


	XL_readStrCell(XL_anaModId, C_anaModId, " ARM_ERR: Model id: string expected",C_result);
	XL_readNumCellWD(XL_dVolSpotUp,C_dVolSpotUp,C_dVolSpotUp_default, " ARM_ERR: Vol Spot Up: numeric expected",C_result);
	XL_readNumCellWD(XL_Asymp_RowUp,C_Asymp_RowUp,C_Asymp_RowUp_default, " ARM_ERR: Asymp_RowUp: numeric expected",C_result);
	XL_readNumCellWD(XL_Asymp_ColUp,C_Asymp_ColUp,C_Asymp_ColUp_default, " ARM_ERR: Asymp_ColUp: numeric expected",C_result);
	XL_readNumCellWD(XL_pUp,C_pUp,C_pUp_default, " ARM_ERR: Probability Up: numeric expected",C_result);
	XL_readNumCellWD(XL_dVolSpotDo,C_dVolSpotDo,C_dVolSpotDo_default, " ARM_ERR: Vol Spot Down: numeric expected",C_result);
	XL_readNumCellWD(XL_Asymp_RowDo,C_Asymp_RowDo,C_Asymp_RowDo_default, " ARM_ERR: Asymp_RowDo: numeric expected",C_result);
	XL_readNumCellWD(XL_Asymp_ColDo,C_Asymp_ColDo,C_Asymp_ColDo_default, " ARM_ERR: Asymp_ColDo: numeric expected",C_result);
	XL_readNumCellWD(XL_pDo,C_pDo,C_pDo_default, " ARM_ERR: Probability Down: numeric expected",C_result);

	if (
		(C_dVolSpotUp < 0)		||
		(C_dVolSpotDo > 0)		||
		(C_pUp < 0)				||
		(C_pDo < 0)				||
		((C_pDo + C_pUp) > 1)
		)
	{
		C_result.setMsg ("ARM_ERR: check your parameters");
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	long retCode;
	long objId;
	CCString prevClass;
	
	CCString curClass = LOCAL_LOGDECANA_CLASS;
	CCString stringId;

	retCode = ARMLOCAL_SMILEDLDCANA(LocalGetNumObjectId (C_anaModId),
														 C_dVolSpotUp,
														 C_Asymp_RowUp,
														 C_Asymp_ColUp,
														 C_pUp,
														 C_dVolSpotDo,
														 C_Asymp_RowDo,
														 C_Asymp_ColDo,
														 C_pDo,
														 C_result);

    if ( retCode == ARM_OK )
    {
	   objId = C_result.getLong();

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
		PXL_ARM_ERR();
	}

	//ARM_END();
	
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_PXL_SMILEDLDCANA" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


__declspec(dllexport) LPXLOPER WINAPI Local_SMILEDLDCFROMANA(LPXLOPER XL_anaModId,
															 LPXLOPER XL_horizon,
															 LPXLOPER XL_nbTraj,
															 LPXLOPER XL_mcmethod)
{
	ADD_LOG("Local_SMILEDLDCFROMANA");
	//ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();

	// C variable
	CCString C_anaModId;
   
    double C_horizon;
	double C_nbTraj;
	
	double C_mcmethod;
	double C_mcmethod_default=0.0;

	
	// error
	static int error;
	static char* reason = "";


    XL_readStrCell(XL_anaModId, C_anaModId, " ARM_ERR: Model id: string expected",C_result);
	XL_readNumCell(XL_horizon,C_horizon," ARM_ERR: Horizon date: date expected",C_result);
    XL_readNumCell(XL_nbTraj, C_nbTraj," ARM_ERR: number of traj. : numeric expected",C_result);
	XL_readNumCellWD(XL_mcmethod, C_mcmethod, C_mcmethod_default, " ARM_ERR: Monte Carlo Method. : numeric expected",C_result);

	long retCode;
	long objId;
	CCString prevClass;
	
	CCString curClass = LOCAL_LOGDECMONTECARLO_CLASS;
	CCString stringId = GetLastCurCellEnvValue ();

	if (!stringId)
	{
        retCode = ARMLOCAL_SMILEDLDCFROMANA(LocalGetNumObjectId (C_anaModId),
											C_horizon,
											(long) C_nbTraj,
											(long) C_mcmethod,
											C_result);

		if ( retCode == ARM_OK )
		{
			objId = C_result.getLong();
			
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
	        retCode = ARMLOCAL_SMILEDLDCFROMANA(LocalGetNumObjectId (C_anaModId),
												C_horizon,
												(long) C_nbTraj,
												(long) C_mcmethod,
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

	        retCode = ARMLOCAL_SMILEDLDCFROMANA(LocalGetNumObjectId (C_anaModId),
												C_horizon,
												(long) C_nbTraj,
												(long) C_mcmethod,
												C_result);
	
			if ( retCode == ARM_OK )
			{
				objId = C_result.getLong();

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

	//ARM_END();
	
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_SMILEDLDCFROMANA" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


__declspec(dllexport) LPXLOPER WINAPI Local_PXL_SMILEDLDCFROMANA(LPXLOPER XL_anaModId,
																 LPXLOPER XL_horizon,
																 LPXLOPER XL_nbTraj,
																 LPXLOPER XL_mcmethod)
{
	ADD_LOG("Local_PXL_SMILEDLDCFROMANA");
	//ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();

	// C variable
	CCString C_anaModId;
   
    double C_horizon;
	double C_nbTraj;
	
	double C_mcmethod;
	double C_mcmethod_default=0.0;
	
	// error
	static int error;
	static char* reason = "";

    XL_readStrCell(XL_anaModId, C_anaModId, " ARM_ERR: Model id: string expected",C_result);
	XL_readNumCell(XL_horizon,C_horizon," ARM_ERR: Horizon date: date expected",C_result);
    XL_readNumCell(XL_nbTraj, C_nbTraj," ARM_ERR: number of traj. : numeric expected",C_result);
	XL_readNumCellWD(XL_mcmethod, C_mcmethod, C_mcmethod_default, " ARM_ERR: Monte Carlo Method. : numeric expected",C_result);

	long retCode;
	long objId;
	CCString prevClass;
	
	CCString curClass = LOCAL_LOGDECMONTECARLO_CLASS;
	CCString stringId;

	retCode = ARMLOCAL_SMILEDLDCFROMANA(LocalGetNumObjectId (C_anaModId),
										C_horizon,
										(long) C_nbTraj,
										(long) C_mcmethod,
										C_result);

    if ( retCode == ARM_OK )
    {
	   objId = C_result.getLong();

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
		PXL_ARM_ERR();
	}

	//ARM_END();
	
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_PXL_SMILEDLDCFROMANA" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


__declspec(dllexport) LPXLOPER WINAPI Local_BSSMILEDMODEL(LPXLOPER XL_date,
														  LPXLOPER XL_spot,
														  LPXLOPER XL_dividend,
														  LPXLOPER XL_discrate,
														  LPXLOPER XL_volat,
														  LPXLOPER XL_typstk,
														  LPXLOPER XL_maturities,
														  LPXLOPER XL_rho,
														  LPXLOPER XL_nu,
														  LPXLOPER XL_isSABR,
														  LPXLOPER XL_beta,
														  LPXLOPER XL_weight,
                                                          LPXLOPER XL_sigmaOrAlphaFlag,
                                                          LPXLOPER XL_correlManager,
														  LPXLOPER XL_rdiscrate,
														  LPXLOPER XL_adjCvxVol,
														  LPXLOPER XL_convToAdjVolWhenAlpha)
{
	ADD_LOG("Local_BSSMILEDMODEL");
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
	    double C_spot;
	    
	    double C_dividend_double;
	    CCString C_dividend_str;
	    long dividend_type;
	    
	    double C_discrate_double;
	    CCString C_discrate_str;
	    long discrate_type;
	    
	    double C_volat_double;
	    CCString C_volat_str;
	    long volat_type;
	    
	    CCString C_typstk;
	    long typstk = K_YIELD;

	    VECTOR<double> C_maturities;
        VECTOR<double> C_maturitiesDef;

        long rhoType;
	    VECTOR<double> C_rho;
        long           rhoObj;
        CCString       C_rhoStr;

        long nuType;
	    VECTOR<double> C_nu;
        long           nuObj;
        CCString       C_nuStr;   


	    CCString C_isSABR;
	    long isSabrId;

   	    double C_beta_double;
   	    VECTOR<double> C_beta_double_vect;
	    CCString C_beta_str;
	    long beta_type;
	    VECTOR<double> C_beta_default;

	    double C_weight;
	    double C_weight_def=0.0;

        double C_sigmaOrAlphaFlag;
        double C_sigmaOrAlphaFlag_def = 1;

		double C_rdiscrate_double;
	    CCString C_rdiscrate_str;
	    long rdiscrate_type;
		
		double C_rdiscrateDef = -1;

        CCString C_correlManager; 
        long C_correlManagerId;

		CCString C_adjCvxVol; 
        long	 C_adjCvxVolId;

		double C_ConversionFlag;
		double C_ConversionFlag_def = 1.0;

	    // error
	    static int error;
	    static char* reason = "";

	    XL_readNumCell(XL_date,C_date," ARM_ERR: date: date expected",C_result);
	    XL_readNumCell(XL_spot,C_spot," ARM_ERR: spot: numeric expected",C_result);

	    XL_readStrOrNumCell(XL_dividend,C_dividend_str,C_dividend_double,dividend_type," ARM_ERR: dividend: string or numeric expected",C_result);
	    XL_readStrOrNumCell(XL_discrate,C_discrate_str,C_discrate_double,discrate_type," ARM_ERR: discrate: string or numeric expected",C_result);
	    XL_readStrOrNumCell(XL_volat,C_volat_str,C_volat_double,volat_type," ARM_ERR: volatility: string or numeric expected",C_result);

	    XL_readStrCellWD(XL_typstk,C_typstk,"YIELD"," ARM_ERR: strike type: string expected",C_result);
	    XL_readNumVectorWD(XL_maturities, C_maturities, C_maturitiesDef, " ARM_ERR: maturities: array of numeric expected",C_result);
	    
        XL_readStrOrNumVector(XL_rho, C_rhoStr, C_rho, rhoType, " ARM_ERR: rho(Rho): array of numeric expected", C_result)
        XL_readStrOrNumVector(XL_nu, C_nuStr, C_nu, nuType," ARM_ERR: nu(Nu): array of numeric expected",C_result);
	    
        XL_readStrCellWD(XL_isSABR,C_isSABR,"LD"," ARM_ERR: smiled model type (LD, SABR_G, SABR_A): string expected",C_result);

		XL_readStrOrNumVectorWD(XL_beta,C_beta_str,C_beta_double_vect,C_beta_default,beta_type," ARM_ERR: beta: string or numeric expected",C_result);

		XL_readNumCellWD(XL_weight,C_weight,C_weight_def," ARM_ERR: weight: numeric expected",C_result);
	    XL_readNumCellWD(XL_sigmaOrAlphaFlag,C_sigmaOrAlphaFlag,C_sigmaOrAlphaFlag_def," ARM_ERR: sigmaOrAlphaFlag: numeric expected(1 or 0)",C_result);
	    
        XL_readStrCellWD(XL_correlManager, C_correlManager,"NULL"," ARM_ERR: CorrelManager id: object expected",C_result);
		XL_readStrOrNumCellWD(XL_rdiscrate,C_rdiscrate_str,C_rdiscrate_double,C_rdiscrateDef,rdiscrate_type," ARM_ERR: rdiscrate: string or numeric expected",C_result);

		XL_readStrCellWD(XL_adjCvxVol, C_adjCvxVol, "NULL", " ARM_ERR: Convexity adjustments: object expected", C_result);
		XL_readNumCellWD(XL_convToAdjVolWhenAlpha, C_ConversionFlag, C_ConversionFlag_def, " ARM_ERR: Conversion when alpha flag: numeric expected", C_result);

		
		C_correlManagerId = ( C_correlManager == "NULL" ? ARM_NULL_OBJECT : LocalGetNumObjectId (C_correlManager) );
		C_adjCvxVolId     = ( C_adjCvxVol     == "NULL" ? ARM_NULL_OBJECT : LocalGetNumObjectId (C_adjCvxVol    ) );
	    
	    if (( C_maturities.size() != C_rho.size() )
		    || ( C_maturities.size() != C_nu.size() ) )
	    {
	       ARM_ARG_ERR();

	       return((LPXLOPER) &XL_result);
	    }

	    if ( ( C_maturities.size() != C_beta_double_vect.size() )
		    && (C_beta_double_vect.size() != 0)
		    && (C_beta_double_vect.size() != 1) )
	    {
	       ARM_ARG_ERR();

	       return((LPXLOPER) &XL_result);
	    }

	    int size = C_maturities.size();

	    if (( typstk = ARM_ConvPriceYield (C_typstk, C_result)) == ARM_DEFAULT_ERR )
	    {
	       ARM_ARG_ERR();

	       return (LPXLOPER)&XL_result;
	    }

        if ( rhoType == XL_TYPE_STRING )
	    {
	       rhoObj = (long) LocalGetNumObjectId(C_rhoStr);

	       rhoType = 1;
	    }
	    else // A vector
	    {
	       rhoType = 0;
	    }

        if ( nuType == XL_TYPE_STRING )
	    {
	       nuObj = (long) LocalGetNumObjectId(C_nuStr);

	       nuType = 1;
	    }
	    else // A vector
	    {
	       nuType = 0;
	    }

	    
	    if ( dividend_type == XL_TYPE_STRING )
	    {
	       C_dividend_double = (double) LocalGetNumObjectId(C_dividend_str);

	       dividend_type = 1;
	    }
	    else
	    {
	       dividend_type = 0;
	    }

	    if ( discrate_type == XL_TYPE_STRING )
	    {
	       C_discrate_double = (double)LocalGetNumObjectId(C_discrate_str);
	       discrate_type = 1;
	    }
	    else
	    {
	       discrate_type = 0;
	    }
		
	    if ( volat_type == XL_TYPE_STRING )
	    {
	       C_volat_double = (double) LocalGetNumObjectId(C_volat_str);

	       volat_type = 1;
	    }
	    else
	    {
	       volat_type = 0;
	    }

	    if ( (isSabrId = ARM_ConvSmiledModelFlag(C_isSABR, C_result)) == ARM_DEFAULT_ERR )
	    {
	       ARM_ARG_ERR();

	       return (LPXLOPER)&XL_result;
	    }

	    if ((XL_beta->xltype == xltypeMissing) || (XL_beta->xltype == xltypeNil))
	    {
		    beta_type = 1;
		    C_beta_double = ARM_NULL_OBJECT;
        }
	    else if ( beta_type == XL_TYPE_STRING )
	    {
		    C_beta_double = (double) LocalGetNumObjectId(C_beta_str);
		    beta_type = 1;
	    }
	    else
	    {
	       beta_type = 0;
	    }
		
		if ( rdiscrate_type == XL_TYPE_STRING )
	    {
	       C_rdiscrate_double = (double)LocalGetNumObjectId(C_rdiscrate_str);
	       rdiscrate_type = 1;
	    }
	    else
	    {
	       rdiscrate_type = 0;
	    }

	    long retCode;
	    long objId;
	    CCString prevClass;
	    
	    CCString curClass = LOCAL_BSMODEL_CLASS;
	    CCString stringId = GetLastCurCellEnvValue ();

	    if (!stringId)
	    {
	       retCode = ARMLOCAL_bssmiledmodel(C_date,
										    C_spot,
										    (long)dividend_type,
										    C_dividend_double,
										    (long) discrate_type,
										    C_discrate_double,
										    (long) volat_type,
										    C_volat_double,
										    typstk,
										    C_maturities,
                                            (long)rhoType,
                                            rhoObj,
										    C_rho,
                                            (long)nuType,
                                            nuObj,
										    C_nu,
										    isSabrId,
										    (long) beta_type,
										    C_beta_double,
										    C_beta_double_vect,
										    C_weight,
                                            (long) C_sigmaOrAlphaFlag,
                                            C_correlManagerId,
											(long) rdiscrate_type,
										    C_rdiscrate_double,
											C_adjCvxVolId,
											(long) C_ConversionFlag,
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

		    objId = LocalGetNumObjectId(stringId);

		    if ( curClass == prevClass )
		    {
			    retCode = ARMLOCAL_bssmiledmodel(C_date,
											    C_spot,
											    (long)dividend_type,
											    C_dividend_double,
											    (long) discrate_type,
											    C_discrate_double,
											    (long) volat_type,
											    C_volat_double,
											    typstk,
											    C_maturities,
											    (long)rhoType,
											    rhoObj,
											    C_rho,
											    (long)nuType,
											    nuObj,
											    C_nu,
											    isSabrId,
											    (long) beta_type,
											    C_beta_double,
											    C_beta_double_vect,
											    C_weight,
                                                (long) C_sigmaOrAlphaFlag,
                                                C_correlManagerId,
												(long) rdiscrate_type,
												C_rdiscrate_double,
												C_adjCvxVolId,
												(long) C_ConversionFlag,
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

			    retCode = ARMLOCAL_bssmiledmodel(C_date,
											    C_spot,
											    (long)dividend_type,
											    C_dividend_double,
											    (long) discrate_type,
											    C_discrate_double,
											    (long) volat_type,
											    C_volat_double,
											    typstk,
											    C_maturities,
											    (long)rhoType,
											    rhoObj,
											    C_rho,
											    (long)nuType,
											    nuObj,
											    C_nu,
											    isSabrId,
											    (long) beta_type,
											    C_beta_double,
											    C_beta_double_vect,
											    C_weight,
                                                (long) C_sigmaOrAlphaFlag,
                                                C_correlManagerId,
												(long) rdiscrate_type,
												C_rdiscrate_double,
												C_adjCvxVolId,
												(long) C_ConversionFlag,
											    C_result);

			    if ( retCode == ARM_OK )
			    {
				    objId = C_result.getLong();

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
        //ARM_END();
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_BSSMILEDMODEL" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}



__declspec(dllexport) LPXLOPER WINAPI Local_PXL_BSSMILEDMODEL(LPXLOPER XL_date,
														  LPXLOPER XL_spot,
														  LPXLOPER XL_dividend,
														  LPXLOPER XL_discrate,
														  LPXLOPER XL_volat,
														  LPXLOPER XL_typstk,
														  LPXLOPER XL_maturities,
														  LPXLOPER XL_rho,
														  LPXLOPER XL_nu,
														  LPXLOPER XL_isSABR,
														  LPXLOPER XL_beta,
														  LPXLOPER XL_weight,
                                                          LPXLOPER XL_sigmaOrAlphaFlag,
                                                          LPXLOPER XL_correlManager,
														  LPXLOPER XL_rdiscrate,
														  LPXLOPER XL_adjCvxVol,
														  LPXLOPER XL_convToAdjVolWhenAlpha)
{
	ADD_LOG("Local_PXL_BSSMILEDMODEL");
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
	    double C_spot;
	    
	    double C_dividend_double;
	    CCString C_dividend_str;
	    long dividend_type;
	    
	    double C_discrate_double;
	    CCString C_discrate_str;
	    long discrate_type;
	    
	    double C_volat_double;
	    CCString C_volat_str;
	    long volat_type;
	    
	    CCString C_typstk;
	    long typstk = K_YIELD;
	    
	    VECTOR<double> C_maturities;
        VECTOR<double>C_maturitiesDef;

        long rhoType;
	    VECTOR<double> C_rho;
        long           rhoObj;
        CCString       C_rhoStr;

        long nuType;
	    VECTOR<double> C_nu;
        long           nuObj;
        CCString       C_nuStr;   

	    CCString C_isSABR;
	    long isSabrId;

   	    double C_beta_double;
   	    VECTOR<double> C_beta_double_vect;
	    CCString C_beta_str;
	    long beta_type;
	    VECTOR<double> C_beta_default;

	    double C_weight;
	    double C_weight_def=0.0;

        double C_sigmaOrAlphaFlag;
        double C_sigmaOrAlphaFlag_def = 1;
		
		double C_rdiscrate_double;
	    CCString C_rdiscrate_str;
	    long rdiscrate_type;
		
		double C_rdiscrateDef = -1;

        CCString C_correlManager; 
        long C_correlManagerId;

		CCString C_adjCvxVol;
		long     C_adjCvxVolId;

		double C_ConversionFlag;
		double C_ConversionFlag_def = 1.0;

	    // error
	    static int error;
	    static char* reason = "";

	    XL_readNumCell(XL_date,C_date," ARM_ERR: date: date expected",C_result);
	    XL_readNumCell(XL_spot,C_spot," ARM_ERR: spot: numeric expected",C_result);

	    XL_readStrOrNumCell(XL_dividend,C_dividend_str,C_dividend_double,dividend_type," ARM_ERR: dividend: string or numeric expected",C_result);
	    XL_readStrOrNumCell(XL_discrate,C_discrate_str,C_discrate_double,discrate_type," ARM_ERR: discrate: string or numeric expected",C_result);
	    XL_readStrOrNumCell(XL_volat,C_volat_str,C_volat_double,volat_type," ARM_ERR: volatility: string or numeric expected",C_result);

	    XL_readStrCellWD(XL_typstk,C_typstk,"YIELD"," ARM_ERR: strike type: string expected",C_result);
	    XL_readNumVectorWD(XL_maturities, C_maturities, C_maturitiesDef," ARM_ERR: maturities: array of numeric expected",C_result);
	    
        XL_readStrOrNumVector(XL_rho, C_rhoStr, C_rho, rhoType, " ARM_ERR: rho(Rho): array of numeric expected", C_result)
        XL_readStrOrNumVector(XL_nu, C_nuStr, C_nu, nuType," ARM_ERR: nu(Nu): array of numeric expected",C_result);
	    
        XL_readStrCellWD(XL_isSABR,C_isSABR,"LD"," ARM_ERR: smiled model type (LD, SABR_G, SABR_A): string expected",C_result);
	    XL_readStrOrNumVectorWD(XL_beta,C_beta_str,C_beta_double_vect,C_beta_default,beta_type," ARM_ERR: beta: string or numeric expected",C_result);
	    
		XL_readNumCellWD(XL_weight,C_weight,C_weight_def," ARM_ERR: weight: numeric expected",C_result);
        XL_readNumCellWD(XL_sigmaOrAlphaFlag,C_sigmaOrAlphaFlag,C_sigmaOrAlphaFlag_def," ARM_ERR: sigmaOrAlphaFlag: numeric expected(O or 1)",C_result);
        XL_readStrCellWD(XL_correlManager, C_correlManager,"NULL"," ARM_ERR: CorrelManager id: object expected",C_result);

		XL_readStrOrNumCellWD(XL_rdiscrate,C_rdiscrate_str,C_rdiscrate_double,C_rdiscrateDef,rdiscrate_type," ARM_ERR: rdiscrate: string or numeric expected",C_result);

		XL_readStrCellWD(XL_adjCvxVol, C_adjCvxVol, "NULL", " ARM_ERR: Convexity adjustments: object expected", C_result);
		XL_readNumCellWD(XL_convToAdjVolWhenAlpha, C_ConversionFlag, C_ConversionFlag_def, " ARM_ERR: conversion when alpha flag: numeric expected", C_result);
        

		C_correlManagerId = ( C_correlManager == "NULL" ? ARM_NULL_OBJECT : LocalGetNumObjectId (C_correlManager) );
		C_adjCvxVolId     = ( C_adjCvxVol     == "NULL" ? ARM_NULL_OBJECT : LocalGetNumObjectId (C_adjCvxVol    ) );

	    if ( (C_maturities.size() != C_rho.size())
		    || (C_maturities.size() != C_nu.size())
		    || (C_rho.size() != C_nu.size()) )
	    {
		    ARM_ARG_ERR();

		    return (LPXLOPER)&XL_result;
	    }

	    int size = C_maturities.size();
	    if (size == 1)
	    {
	    }

	    typstk = ARM_ConvPriceYield (C_typstk, C_result);

	    if ( typstk == ARM_DEFAULT_ERR )
	    {
	       ARM_ARG_ERR();

	       return (LPXLOPER)&XL_result;
	    }

        if ( rhoType == XL_TYPE_STRING )
	    {
	       rhoObj = (long) LocalGetNumObjectId(C_rhoStr);

	       rhoType = 1;
	    }
	    else // A vector
	    {
	       rhoType = 0;
	    }

        if ( nuType == XL_TYPE_STRING )
	    {
	       nuObj = (long) LocalGetNumObjectId(C_nuStr);

	       nuType = 1;
	    }
	    else // A vector
	    {
	       nuType = 0;
	    }

	    if ( dividend_type == XL_TYPE_STRING )
	    {
	       C_dividend_double = (double) LocalGetNumObjectId(C_dividend_str);

	       dividend_type = 1;
	    }
	    else
	    {
	       dividend_type = 0;
	    }

	    if ( discrate_type == XL_TYPE_STRING )
	    {
	       C_discrate_double = (double)LocalGetNumObjectId (C_discrate_str);
	       discrate_type = 1;
	    }
	    else
	    {
	       discrate_type = 0;
	    }

	    if ( volat_type == XL_TYPE_STRING )
	    {
	       C_volat_double = (double) LocalGetNumObjectId (C_volat_str);

	       volat_type = 1;
	    }
	    else
	    {
	       volat_type = 0;
	    }

	    if ( (isSabrId = ARM_ConvSmiledModelFlag (C_isSABR, C_result)) == ARM_DEFAULT_ERR)
	    {
	       ARM_ARG_ERR();

	       return (LPXLOPER)&XL_result;
	    }

	    if ((XL_beta->xltype == xltypeMissing) || (XL_beta->xltype == xltypeNil))
	    {
		    beta_type = 1;
		    C_beta_double = ARM_NULL_OBJECT;
        }
	    else if ( beta_type == XL_TYPE_STRING )
	    {
		    C_beta_double = (double) LocalGetNumObjectId(C_beta_str);
		    beta_type = 1;
	    }
	    else
	    {
		    beta_type = 0;
	    }
	    
		if ( rdiscrate_type == XL_TYPE_STRING )
	    {
	       C_rdiscrate_double = (double)LocalGetNumObjectId (C_rdiscrate_str);
	       rdiscrate_type = 1;
	    }
	    else
	    {
	       rdiscrate_type = 0;
	    }

	    long retCode;
	    long objId;
	    
	    CCString curClass = LOCAL_BSMODEL_CLASS;
	    CCString stringId;

	    retCode = ARMLOCAL_bssmiledmodel(C_date,
									    C_spot,
									    (long)dividend_type,
									    C_dividend_double,
									    (long) discrate_type,
									    C_discrate_double,
									    (long) volat_type,
									    C_volat_double,
									    typstk,
									    C_maturities,
									    (long)rhoType,
									    rhoObj,
									    C_rho,
									    (long)nuType,
									    nuObj,
									    C_nu,
									    isSabrId,
									    (long) beta_type,
									    C_beta_double,
									    C_beta_double_vect,
									    C_weight,
                                        (long) C_sigmaOrAlphaFlag,
                                        C_correlManagerId,
										(long) rdiscrate_type,
									    C_rdiscrate_double,
										(long) C_adjCvxVolId,
										(long) C_ConversionFlag,
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
	       PXL_ARM_ERR();
	    }

    //	ARM_END();	
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_PXL_BSSMILEDMODEL" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


//------------------------------------------------------------------------------------------//
// BS Smiled Model with SABR structure (class ARM_SABRVol)									//
//------------------------------------------------------------------------------------------//																							
__declspec(dllexport) LPXLOPER WINAPI Local_SABR_BSSMILEDMODEL(LPXLOPER XL_date,
															   LPXLOPER XL_spot,
															   LPXLOPER XL_dividend,
															   LPXLOPER XL_discrate,
															   LPXLOPER XL_SABRVol,
															   LPXLOPER XL_typstk,
															   LPXLOPER XL_correlManager,
															   LPXLOPER XL_rdiscrate,
															   LPXLOPER XL_adjCvxVol,
															   LPXLOPER XL_convToAdjVolWhenAlpha)
{
	ADD_LOG("Local_SABR_BSSMILEDMODEL");
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;


     /// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
        ARM_NOCALCIFWIZ();

	    // C variable
	    double   C_date;
	    double   C_spot;
	    
	    double   C_dividend_double;
		CCString C_dividend_str;
	    long	 dividend_type;
	    
	    double	 C_discrate_double;
	    CCString C_discrate_str;
	    long	 discrate_type;
	    
	    CCString C_typstk;
	    long	 typstk = K_YIELD;

		double	 C_rdiscrate_double;
	    CCString C_rdiscrate_str;
	    long	 rdiscrate_type;
		
		double   C_rdiscrateDef = -1.0;

		CCString C_correlManager; 
        long	 C_correlManagerId;

		CCString C_adjCvxVol; 
        long	 C_adjCvxVolId;

		CCString C_SABRVol;
		long	 C_SABRVolId;

		double   C_ConversionFlag;
		double   C_ConversionFlag_def = 1.0;

	    // error
	    static int error;
	    static char* reason = "";

		XL_readNumCell(XL_date,C_date," ARM_ERR: date: date expected",C_result);
		XL_readNumCell(XL_spot,C_spot," ARM_ERR: spot: numeric expected",C_result);

		XL_readStrOrNumCell(XL_dividend,C_dividend_str,C_dividend_double,dividend_type," ARM_ERR: dividend: string or numeric expected",C_result);
		XL_readStrOrNumCell(XL_discrate,C_discrate_str,C_discrate_double,discrate_type," ARM_ERR: discrate: string or numeric expected",C_result);
		XL_readStrOrNumCellWD(XL_rdiscrate,C_rdiscrate_str,C_rdiscrate_double,C_rdiscrateDef,rdiscrate_type," ARM_ERR: rdiscrate: string or numeric expected",C_result);
		
		XL_readStrCellWD(XL_SABRVol, C_SABRVol, "NULL", " ARM_ERR: SABR volatility: object expected", C_result);
		XL_readStrCellWD(XL_typstk,C_typstk,"YIELD"," ARM_ERR: strike type: string expected",C_result);    
		XL_readStrCellWD(XL_correlManager, C_correlManager,"NULL"," ARM_ERR: CorrelManager id: object expected",C_result);
		XL_readStrCellWD(XL_adjCvxVol, C_adjCvxVol, "NULL", " ARM_ERR: Convexity adjustments: object expected", C_result);
		XL_readNumCellWD(XL_convToAdjVolWhenAlpha, C_ConversionFlag, C_ConversionFlag_def, " ARM_ERR: conversion when alpha flag: numeric expected", C_result);
	
	    if (( typstk = ARM_ConvPriceYield (C_typstk, C_result)) == ARM_DEFAULT_ERR )
	    {
	       ARM_ARG_ERR();

	       return (LPXLOPER)&XL_result;
	    }
	    
	    if ( dividend_type == XL_TYPE_STRING )
	    {
			C_dividend_double = (double) LocalGetNumObjectId (C_dividend_str);
			dividend_type = 1;
	    }
	    else
	    {
			dividend_type = 0;
	    }

	    if ( discrate_type == XL_TYPE_STRING )
	    {
			C_discrate_double = (double) LocalGetNumObjectId (C_discrate_str);
	        discrate_type = 1;
	    }
	    else
	    {
	        discrate_type = 0;
	    }

		if ( rdiscrate_type == XL_TYPE_STRING )
	    {
			C_rdiscrate_double = (double) LocalGetNumObjectId (C_rdiscrate_str);
	        rdiscrate_type = 1;
	    }
	    else
	    {
	        rdiscrate_type = 0;
	    }

        C_correlManagerId = ( C_correlManager == "NULL" ? ARM_NULL_OBJECT : LocalGetNumObjectId (C_correlManager) );
		C_adjCvxVolId	  = ( C_adjCvxVol     == "NULL" ? ARM_NULL_OBJECT : LocalGetNumObjectId (C_adjCvxVol    ) );
		C_SABRVolId		  = ( C_SABRVol       == "NULL" ? ARM_NULL_OBJECT : LocalGetNumObjectId (C_SABRVol      ) );
		
	    
		long retCode;
	    long objId;
	    CCString prevClass;
	    
	    CCString curClass = LOCAL_BSMODEL_CLASS;
	    CCString stringId = GetLastCurCellEnvValue ();

	    if (!stringId)
	    {
	       retCode = ARMLOCAL_bssmiledmodel_SABR(C_date,
												 C_spot,
												 (long)dividend_type,
												 C_dividend_double,
												 (long) discrate_type,
												 C_discrate_double,
												 C_SABRVolId,
												 typstk,
												 C_correlManagerId,
												 (long) rdiscrate_type,
												 C_rdiscrate_double,
												 C_adjCvxVolId,
												 (long) C_ConversionFlag,
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

		    objId = LocalGetNumObjectId(stringId);

		    if ( curClass == prevClass )
		    {
			    retCode = ARMLOCAL_bssmiledmodel_SABR(C_date,
													  C_spot,
													  (long)dividend_type,
													  C_dividend_double,
													  (long) discrate_type,
													  C_discrate_double,
													  C_SABRVolId,
													  typstk,
													  C_correlManagerId,
													  (long) rdiscrate_type,
													  C_rdiscrate_double,
													  C_adjCvxVolId,
													  (long) C_ConversionFlag,
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

			    retCode = ARMLOCAL_bssmiledmodel_SABR(C_date,
													  C_spot,
													  (long)dividend_type,
													  C_dividend_double,
													  (long) discrate_type,
													  C_discrate_double,
													  C_SABRVolId,
													  typstk,
													  C_correlManagerId,
													  (long) rdiscrate_type,
													  C_rdiscrate_double,
													  C_adjCvxVolId,
													  (long) C_ConversionFlag,
													  C_result);

			    if ( retCode == ARM_OK )
			    {
				    objId = C_result.getLong();

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
        //ARM_END();
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_SABR_BSSMILEDMODEL" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


//------------------------------------------------------------------------------------------//
// BS Smiled Model with SABR structure (class ARM_SABRVol)									//
// Function to use with Excel VBA															//
//------------------------------------------------------------------------------------------//																							
__declspec(dllexport) LPXLOPER WINAPI Local_PXL_SABR_BSSMILEDMODEL(LPXLOPER XL_date,
																   LPXLOPER XL_spot,
																   LPXLOPER XL_dividend,
																   LPXLOPER XL_discrate,
																   LPXLOPER XL_SABRVol,
																   LPXLOPER XL_typstk,
																   LPXLOPER XL_correlManager,
																   LPXLOPER XL_rdiscrate,
																   LPXLOPER XL_adjCvxVol,
																   LPXLOPER XL_convToAdjVolWhenAlpha)
{
	ADD_LOG("Local_PXL_SABR_BSSMILEDMODEL");
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;


     /// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
        ARM_NOCALCIFWIZ();

	    // C variable
	    double   C_date;
	    double   C_spot;
	    
	    double   C_dividend_double;
		CCString C_dividend_str;
	    long	 dividend_type;
	    
	    double	 C_discrate_double;
	    CCString C_discrate_str;
	    long	 discrate_type;
	    
	    CCString C_typstk;
	    long	 typstk = K_YIELD;

		double	 C_rdiscrate_double;
	    CCString C_rdiscrate_str;
	    long	 rdiscrate_type;
		
		double   C_rdiscrateDef = -1.0;

		CCString C_correlManager; 
        long	 C_correlManagerId;

		CCString C_adjCvxVol; 
        long	 C_adjCvxVolId;

		CCString C_SABRVol;
		long	 C_SABRVolId;

		double   C_ConversionFlag;
		double   C_ConversionFlag_def = 1.0;

	    // error
	    static int error;
	    static char* reason = "";

		XL_readNumCell(XL_date,C_date," ARM_ERR: date: date expected",C_result);
		XL_readNumCell(XL_spot,C_spot," ARM_ERR: spot: numeric expected",C_result);

		XL_readStrOrNumCell(XL_dividend,C_dividend_str,C_dividend_double,dividend_type," ARM_ERR: dividend: string or numeric expected",C_result);
		XL_readStrOrNumCell(XL_discrate,C_discrate_str,C_discrate_double,discrate_type," ARM_ERR: discrate: string or numeric expected",C_result);
		XL_readStrOrNumCellWD(XL_rdiscrate,C_rdiscrate_str,C_rdiscrate_double,C_rdiscrateDef,rdiscrate_type," ARM_ERR: rdiscrate: string or numeric expected",C_result);
		
		XL_readStrCellWD(XL_SABRVol, C_SABRVol, "NULL", " ARM_ERR: SABR volatility: object expected", C_result);
		XL_readStrCellWD(XL_typstk,C_typstk,"YIELD"," ARM_ERR: strike type: string expected",C_result);    
		XL_readStrCellWD(XL_correlManager, C_correlManager,"NULL"," ARM_ERR: CorrelManager id: object expected",C_result);
		XL_readStrCellWD(XL_adjCvxVol, C_adjCvxVol, "NULL", " ARM_ERR: Convexity adjustments: object expected", C_result);
		XL_readNumCellWD(XL_convToAdjVolWhenAlpha, C_ConversionFlag, C_ConversionFlag_def, " ARM_ERR: conversion when alpha flag: numeric expected", C_result);
	
	    if (( typstk = ARM_ConvPriceYield (C_typstk, C_result)) == ARM_DEFAULT_ERR )
	    {
	       ARM_ARG_ERR();

	       return (LPXLOPER)&XL_result;
	    }
	    
	    if ( dividend_type == XL_TYPE_STRING )
	    {
			C_dividend_double = (double) LocalGetNumObjectId (C_dividend_str);
			dividend_type = 1;
	    }
	    else
	    {
			dividend_type = 0;
	    }

	    if ( discrate_type == XL_TYPE_STRING )
	    {
			C_discrate_double = (double) LocalGetNumObjectId (C_discrate_str);
	        discrate_type = 1;
	    }
	    else
	    {
	        discrate_type = 0;
	    }

		if ( rdiscrate_type == XL_TYPE_STRING )
	    {
			C_rdiscrate_double = (double) LocalGetNumObjectId (C_rdiscrate_str);
	        rdiscrate_type = 1;
	    }
	    else
	    {
	        rdiscrate_type = 0;
	    }

        C_correlManagerId = ( C_correlManager == "NULL" ? ARM_NULL_OBJECT : LocalGetNumObjectId (C_correlManager) );
		C_adjCvxVolId	  = ( C_adjCvxVol     == "NULL" ? ARM_NULL_OBJECT : LocalGetNumObjectId (C_adjCvxVol    ) );
		C_SABRVolId		  = ( C_SABRVol       == "NULL" ? ARM_NULL_OBJECT : LocalGetNumObjectId (C_SABRVol      ) );
		
	    
		long retCode;
	    long objId;
	    CCString prevClass;
	    
	    CCString curClass = LOCAL_BSMODEL_CLASS;
	    CCString stringId = GetLastCurCellEnvValue ();


	   retCode = ARMLOCAL_bssmiledmodel_SABR(C_date,
											 C_spot,
											 (long)dividend_type,
											 C_dividend_double,
											 (long) discrate_type,
											 C_discrate_double,
											 C_SABRVolId,
											 typstk,
											 C_correlManagerId,
											 (long) rdiscrate_type,
											 C_rdiscrate_double,
											 C_adjCvxVolId,
											 (long) C_ConversionFlag,
											 C_result);

		if ( retCode == ARM_OK )
		{
			objId    = C_result.getLong (); 

			stringId = LocalMakeObjectId (objId, curClass);
		}

	    if ( retCode == ARM_OK )
	    {			
	       XL_result.xltype   = xltypeStr;
	       XL_result.val.str  = XL_StrC2StrPascal (stringId);
	       XL_result.xltype  |= xlbitDLLFree;
	    }
	    else
	    {
	       PXL_ARM_ERR();
	    }
        //ARM_END();
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_PXL_SABR_BSSMILEDMODEL" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


__declspec(dllexport) LPXLOPER WINAPI Local_CRRTREE(LPXLOPER XL_begDate,
													LPXLOPER XL_endDate,
													LPXLOPER XL_nbSteps,
													LPXLOPER XL_spot,
													LPXLOPER XL_dividend,
													LPXLOPER XL_discrate,
													LPXLOPER XL_volat)
{
	ADD_LOG("Local_CRRTREE");
	//ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();

	// C variable   
    double C_begDate;
	double C_endDate;
	
		
	double C_nbSteps;
	double C_spot;
	double C_dividend;
    double C_discrate;
	double C_volat;
	
	
	// error
	static int error;
	static char* reason = "";


	XL_readNumCell(XL_begDate,C_begDate," ARM_ERR: start date: date expected",C_result);
    XL_readNumCell(XL_endDate,C_endDate," ARM_ERR: end date: date expected",C_result);
    XL_readNumCell(XL_nbSteps, C_nbSteps," ARM_ERR: number of steps: numeric expected",C_result);
	XL_readNumCell(XL_spot,C_spot," ARM_ERR: spot rate: numeric expected",C_result);
    XL_readNumCell(XL_dividend, C_dividend," ARM_ERR: dividend: numeric expected",C_result);
    XL_readNumCell(XL_discrate, C_discrate," ARM_ERR: discount rate: numeric expected",C_result);
    XL_readNumCell(XL_volat, C_volat," ARM_ERR: volatility: numeric expected",C_result);

	long retCode;
	long objId;
	CCString prevClass;
	
	CCString curClass = LOCAL_CRR_TREE_CLASS;
	CCString stringId = GetLastCurCellEnvValue ();

	if (!stringId)
	{
        retCode = ARMLOCAL_CRRTREE(C_begDate,
								   C_endDate,
								   (long) C_nbSteps,
								   C_spot,
								   C_dividend,
								   C_discrate,
								   C_volat,
								   C_result);

		if ( retCode == ARM_OK )
		{
			objId = C_result.getLong();
			
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
           retCode = ARMLOCAL_CRRTREE(C_begDate,
									  C_endDate,
									  (long) C_nbSteps,
									  C_spot,
									  C_dividend,
									  C_discrate,
									  C_volat,
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

            retCode = ARMLOCAL_CRRTREE(C_begDate,
									   C_endDate,
									   (long) C_nbSteps,
									   C_spot,
									   C_dividend,
									   C_discrate,
									   C_volat,
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

	//ARM_END();
	
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_CRRTREE" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}



__declspec(dllexport) LPXLOPER WINAPI Local_PXL_CRRTREE(LPXLOPER XL_begDate,
														LPXLOPER XL_endDate,
														LPXLOPER XL_nbSteps,
														LPXLOPER XL_spot,
														LPXLOPER XL_dividend,
														LPXLOPER XL_discrate,
														LPXLOPER XL_volat)
{
	ADD_LOG("Local_PXL_CRRTREE");
	//ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();

	// C variable   
    double C_begDate;
	double C_endDate;
	
		
	double C_nbSteps;
	double C_spot;
	double C_dividend;
    double C_discrate;
	double C_volat;
	
	
	// error
	static int error;
	static char* reason = "";


	XL_readNumCell(XL_begDate,C_begDate," ARM_ERR: start date: date expected",C_result);
    XL_readNumCell(XL_endDate,C_endDate," ARM_ERR: end date: date expected",C_result);
    XL_readNumCell(XL_nbSteps, C_nbSteps," ARM_ERR: number of steps: numeric expected",C_result);
	XL_readNumCell(XL_spot,C_spot," ARM_ERR: spot rate: numeric expected",C_result);
    XL_readNumCell(XL_dividend, C_dividend," ARM_ERR: dividend: numeric expected",C_result);
    XL_readNumCell(XL_discrate, C_discrate," ARM_ERR: discount rate: numeric expected",C_result);
    XL_readNumCell(XL_volat, C_volat," ARM_ERR: volatility: numeric expected",C_result);

	long retCode;
	long objId;
	
	CCString curClass = LOCAL_CRR_TREE_CLASS;
	CCString stringId;

	retCode = ARMLOCAL_CRRTREE(C_begDate,
							   C_endDate,
							   (long) C_nbSteps,
							   C_spot,
							   C_dividend,
							   C_discrate,
							   C_volat,
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
		PXL_ARM_ERR();
	}

	//ARM_END();
	
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_PXL_CRRTREE" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


__declspec(dllexport) LPXLOPER WINAPI Local_BSCORRMODEL(LPXLOPER XL_startDate,
														LPXLOPER XL_zeroCurve,
														LPXLOPER XL_spreadLock,
														LPXLOPER XL_capIrgVol,
														LPXLOPER XL_capCashVol,
														LPXLOPER XL_indexVolAdj,
														LPXLOPER XL_spreadVol,
														LPXLOPER XL_correlations,
														LPXLOPER XL_modelType,
														LPXLOPER XL_volType)
{
	ADD_LOG("Local_BSCORRMODEL");
	//ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();

	// C variable   
	double C_startDate;

	CCString C_zc;
	CCString C_spreadlock;
	CCString C_capirgvol;

	CCString C_capcashvol;
	long capcashvolId;
	CCString C_indexadjvol;
	long indexadjvolId;
	CCString C_spreadvol;
	long spreadvolId;
	CCString C_correlations;
	long correlationsId;

	CCString C_modelType;
	long modelTypeId;
	CCString C_volType;
	long volTypeId;

	// error
	static int error;
	static char* reason = "";

	XL_readNumCell(XL_startDate,C_startDate," ARM_ERR: start date: date expected",C_result);
	XL_readStrCell(XL_zeroCurve, C_zc," ARM_ERR: zero curve: object expected",C_result);
	XL_readStrCell(XL_spreadLock, C_spreadlock," ARM_ERR: Spread Lock: object expected",C_result);
	XL_readStrCell(XL_capIrgVol, C_capirgvol," ARM_ERR: Cap IRG Vol: object expected",C_result);
	XL_readStrCellWD(XL_capCashVol, C_capcashvol,"DEFAULT"," ARM_ERR: Cap Cash Vol: object expected",C_result);
	XL_readStrCellWD(XL_indexVolAdj, C_indexadjvol,"DEFAULT"," ARM_ERR: Index Adj Vol: object expected",C_result);
	XL_readStrCellWD(XL_spreadVol, C_spreadvol,"DEFAULT"," ARM_ERR: Spread Vol: object expected",C_result);
	XL_readStrCellWD(XL_correlations, C_correlations,"DEFAULT"," ARM_ERR: Correlations: object expected",C_result);
	XL_readStrCellWD(XL_modelType, C_modelType,"2LOG"," ARM_ERR: model type: string expected",C_result);
	XL_readStrCellWD(XL_volType, C_volType,"INPUT"," ARM_ERR: vol type: string expected",C_result);

	if ((modelTypeId = ARM_ConvModelType (C_modelType, C_result)) == ARM_DEFAULT_ERR)
	{
	   ARM_ARG_ERR();
	   return (LPXLOPER)&XL_result;
	}

	if ((volTypeId = ARM_ConvVolType2 (C_volType, C_result)) == ARM_DEFAULT_ERR)
	{
	   ARM_ARG_ERR();
	   return (LPXLOPER)&XL_result;
	}

	if (C_capcashvol == "DEFAULT")
		capcashvolId = ARM_NULL_OBJECT;
	else
		capcashvolId = LocalGetNumObjectId(C_capcashvol);

	if (C_indexadjvol == "DEFAULT")
		indexadjvolId = ARM_NULL_OBJECT;
	else
		indexadjvolId = LocalGetNumObjectId(C_indexadjvol);

	if (C_spreadvol == "DEFAULT")
		spreadvolId = ARM_NULL_OBJECT;
	else
		spreadvolId = LocalGetNumObjectId(C_spreadvol);

	if (C_correlations == "DEFAULT")
		correlationsId = ARM_NULL_OBJECT;
	else
		correlationsId = LocalGetNumObjectId(C_correlations);

	long retCode;
	long objId;
	CCString prevClass;
	
	CCString curClass = LOCAL_BSCORRMODEL_CLASS;
	CCString stringId = GetLastCurCellEnvValue ();

	if (!stringId)
	{
        retCode = ARMLOCAL_BSCORRMODEL(C_startDate,
									   LocalGetNumObjectId(C_zc),
									   LocalGetNumObjectId(C_spreadlock),
									   LocalGetNumObjectId(C_capirgvol),
									   capcashvolId,
									   indexadjvolId,
									   spreadvolId,
									   correlationsId,
									   modelTypeId,
									   volTypeId,
									   C_result);

		if ( retCode == ARM_OK )
		{
			objId = C_result.getLong();
			
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
			retCode = ARMLOCAL_BSCORRMODEL(C_startDate,
										   LocalGetNumObjectId(C_zc),
										   LocalGetNumObjectId(C_spreadlock),
										   LocalGetNumObjectId(C_capirgvol),
										   capcashvolId,
										   indexadjvolId,
										   spreadvolId,
										   correlationsId,
										   modelTypeId,
										   volTypeId,
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

			retCode = ARMLOCAL_BSCORRMODEL(C_startDate,
										   LocalGetNumObjectId(C_zc),
										   LocalGetNumObjectId(C_spreadlock),
										   LocalGetNumObjectId(C_capirgvol),
										   capcashvolId,
										   indexadjvolId,
										   spreadvolId,
										   correlationsId,
										   modelTypeId,
										   volTypeId,
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

	//ARM_END();
	
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_BSCORRMODEL" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}

__declspec(dllexport) LPXLOPER WINAPI Local_PXL_BSCORRMODEL(LPXLOPER XL_startDate,
															LPXLOPER XL_zeroCurve,
															LPXLOPER XL_spreadLock,
															LPXLOPER XL_capIrgVol,
															LPXLOPER XL_capCashVol,
															LPXLOPER XL_indexVolAdj,
															LPXLOPER XL_spreadVol,
															LPXLOPER XL_correlations,
															LPXLOPER XL_modelType,
															LPXLOPER XL_volType)
{
	ADD_LOG("Local_PXL_BSCORRMODEL");
	//ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();

	// C variable   
	double C_startDate;

	CCString C_zc;
	CCString C_spreadlock;
	CCString C_capirgvol;

	CCString C_capcashvol;
	long capcashvolId;
	CCString C_indexadjvol;
	long indexadjvolId;
	CCString C_spreadvol;
	long spreadvolId;
	CCString C_correlations;
	long correlationsId;

	CCString C_modelType;
	long modelTypeId;
	CCString C_volType;
	long volTypeId;

	// error
	static int error;
	static char* reason = "";

	XL_readNumCell(XL_startDate,C_startDate," ARM_ERR: start date: date expected",C_result);
	XL_readStrCell(XL_zeroCurve, C_zc," ARM_ERR: zero curve: object expected",C_result);
	XL_readStrCell(XL_spreadLock, C_spreadlock," ARM_ERR: Spread Lock: object expected",C_result);
	XL_readStrCell(XL_capIrgVol, C_capirgvol," ARM_ERR: Cap IRG Vol: object expected",C_result);
	XL_readStrCellWD(XL_capCashVol, C_capcashvol,"DEFAULT"," ARM_ERR: Cap Cash Vol: object expected",C_result);
	XL_readStrCellWD(XL_indexVolAdj, C_indexadjvol,"DEFAULT"," ARM_ERR: Index Adj Vol: object expected",C_result);
	XL_readStrCellWD(XL_spreadVol, C_spreadvol,"DEFAULT"," ARM_ERR: Spread Vol: object expected",C_result);
	XL_readStrCellWD(XL_correlations, C_correlations,"DEFAULT"," ARM_ERR: Correlations: object expected",C_result);
	XL_readStrCellWD(XL_modelType, C_modelType,"2LOG"," ARM_ERR: model type: string expected",C_result);
	XL_readStrCellWD(XL_volType, C_volType,"INPUT"," ARM_ERR: vol type: string expected",C_result);

	if ((modelTypeId = ARM_ConvModelType (C_modelType, C_result)) == ARM_DEFAULT_ERR)
	{
	   ARM_ARG_ERR();
	   return (LPXLOPER)&XL_result;
	}

	if ((volTypeId = ARM_ConvVolType2 (C_volType, C_result)) == ARM_DEFAULT_ERR)
	{
	   ARM_ARG_ERR();
	   return (LPXLOPER)&XL_result;
	}

	if (C_capcashvol == "DEFAULT")
		capcashvolId = ARM_NULL_OBJECT;
	else
		capcashvolId = LocalGetNumObjectId(C_capcashvol);

	if (C_indexadjvol == "DEFAULT")
		indexadjvolId = ARM_NULL_OBJECT;
	else
		indexadjvolId = LocalGetNumObjectId(C_indexadjvol);

	if (C_spreadvol == "DEFAULT")
		spreadvolId = ARM_NULL_OBJECT;
	else
		spreadvolId = LocalGetNumObjectId(C_spreadvol);

	if (C_correlations == "DEFAULT")
		correlationsId = ARM_NULL_OBJECT;
	else
		correlationsId = LocalGetNumObjectId(C_correlations);

	long retCode;
	long objId;
	
	CCString curClass = LOCAL_BSCORRMODEL_CLASS;
	CCString stringId;

	retCode = ARMLOCAL_BSCORRMODEL(C_startDate,
								   LocalGetNumObjectId(C_zc),
								   LocalGetNumObjectId(C_spreadlock),
								   LocalGetNumObjectId(C_capirgvol),
								   capcashvolId,
								   indexadjvolId,
								   spreadvolId,
								   correlationsId,
								   modelTypeId,
								   volTypeId,
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
		PXL_ARM_ERR();
	}

	//ARM_END();
	
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_PXL_BSCORRMODEL" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


__declspec(dllexport) LPXLOPER WINAPI Local_FRMTREE_AUTO_B(LPXLOPER XL_zcId,
														   LPXLOPER XL_zdId,
														   LPXLOPER XL_volId,
														   LPXLOPER XL_smileId,
														   LPXLOPER XL_autoMode,
														   LPXLOPER XL_horizon,
														   LPXLOPER XL_fineMonth,
														   LPXLOPER XL_shapeDecay,
														   LPXLOPER XL_shapeSlope,
														   LPXLOPER XL_shapeAsymptote,
														   LPXLOPER XL_nbFactor,
														   LPXLOPER XL_corrMatu,
														   LPXLOPER XL_corrMatrix,
														   LPXLOPER XL_corr2Matu)
{
	ADD_LOG("Local_FRMTREE_AUTO_B");
	//ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();

	// C variable
	CCString C_zcId;
	CCString C_zdId;

    CCString C_volId;
    CCString C_smileId;
    CCString CS_autoMode;
    CCString CS_fineMonth;
    double C_horizon;
    double C_shapeDecay;
    double C_shapeSlope;
    double C_shapeAsymptote;
    double C_nbFactor;
   
    VECTOR<double> C_corrMatu;
    VECTOR<double> C_corrMatrix;
    VECTOR<double> C_corr2Matu;	

    double mFineDefault = K_L_DEFAULT;
    double epsDefault = 1e-8;
    double huegeDefault = 1e8;
    double zero = 0;
    double one = 1;
	
	// error
	static int error;
	static char* reason = "";

    XL_readStrCell(XL_zcId, C_zcId, " ARM_ERR: Curve id: string expected",C_result);
	XL_readStrCell(XL_zdId, C_zdId, " ARM_ERR: Curve id: string expected",C_result);
    XL_readStrCell(XL_volId, C_volId, " ARM_ERR: Vol id: string expected",C_result);
    XL_readStrCell(XL_smileId, C_smileId, " ARM_ERR: Smile id: string expected",C_result);
    XL_readStrCell(XL_autoMode, CS_autoMode, " ARM_ERR: AutoMode: string expected",C_result);
	XL_readNumCell(XL_horizon,C_horizon," ARM_ERR: Horizon date: date expected",C_result);
    XL_readStrCellWD(XL_fineMonth, CS_fineMonth, "L_DEFAULT", " ARM_ERR: fine month. : string expected",C_result);
    XL_readNumCellWD(XL_shapeDecay, C_shapeDecay, zero, " ARM_ERR: shape decay. : numeric expected",C_result);
    XL_readNumCellWD(XL_shapeSlope, C_shapeSlope, zero, " ARM_ERR: shape slope. : numeric expected",C_result);
    XL_readNumCellWD(XL_shapeAsymptote, C_shapeAsymptote, zero, " ARM_ERR: shape Asymptote. : numeric expected",C_result);
    XL_readNumCellWD(XL_nbFactor, C_nbFactor, one, " ARM_ERR: Nb factor. : numeric expected",C_result);

   	if((XL_corrMatu->xltype == xltypeMissing) || (XL_corrMatu->xltype == xltypeNil))
	{
    }
	else
	{
        XL_readNumVector(XL_corrMatu, C_corrMatu, " ARM_ERR: Correlation maturities: array of numeric expected",C_result);
    }

   	if((XL_corrMatrix->xltype == xltypeMissing) || (XL_corrMatrix->xltype == xltypeNil))
	{
    }
	else
	{
        XL_readNumVector(XL_corrMatrix, C_corrMatrix, " ARM_ERR: Correlation matrix: array of numeric expected",C_result);
    }

   	if((XL_corr2Matu->xltype == xltypeMissing) || (XL_corr2Matu->xltype == xltypeNil))
	{
    }
	else
	{
        XL_readNumVector(XL_corr2Matu, C_corr2Matu, " ARM_ERR: Correlation correlated maturities : array of numeric expected",C_result);
    }

	long autoModeId;
	long fineMonthId;

	if((autoModeId = ARM_ConvAutoMode2 (CS_autoMode, C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	if((fineMonthId = ARM_ConvFineMode (CS_fineMonth, C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	long retCode;
	long objId;
	CCString prevClass;

	CCString curClass = LOCAL_FRMTREE_CLASS;
	CCString stringId = GetLastCurCellEnvValue ();

	if (!stringId)
	{
        retCode = ARMLOCAL_FRMTREE_AUTO_B(LocalGetNumObjectId(C_zcId),
									 LocalGetNumObjectId(C_zdId),
									 LocalGetNumObjectId(C_volId),
									 LocalGetNumObjectId(C_smileId),
									 autoModeId,
									 C_horizon,
									 fineMonthId,
									 C_shapeDecay,
									 C_shapeSlope,
									 C_shapeAsymptote,
									 (long) C_nbFactor,
									 C_corrMatu,
									 C_corrMatrix,
									 C_corr2Matu,
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
			retCode = ARMLOCAL_FRMTREE_AUTO_B(LocalGetNumObjectId(C_zcId),
										 LocalGetNumObjectId(C_zdId),
										 LocalGetNumObjectId(C_volId),
										 LocalGetNumObjectId(C_smileId),
										 autoModeId,
										 C_horizon,
										 fineMonthId,
										 C_shapeDecay,
										 C_shapeSlope,
										 C_shapeAsymptote,
										 (long) C_nbFactor,
										 C_corrMatu,
										 C_corrMatrix,
										 C_corr2Matu,
										 C_result,
										 objId);
		}
		else
		{
			FreeCurCellContent ();

			retCode = ARMLOCAL_FRMTREE_AUTO_B(LocalGetNumObjectId(C_zcId),
										 LocalGetNumObjectId(C_zdId),
										 LocalGetNumObjectId(C_volId),
										 LocalGetNumObjectId(C_smileId),
										 autoModeId,
										 C_horizon,
										 fineMonthId,
										 C_shapeDecay,
										 C_shapeSlope,
										 C_shapeAsymptote,
										 (long) C_nbFactor,
										 C_corrMatu,
										 C_corrMatrix,
										 C_corr2Matu,
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

	//ARM_END();
	
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_FRMTREE_AUTO_B" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


__declspec(dllexport) LPXLOPER WINAPI Local_PXL_FRMTREE_AUTO_B(LPXLOPER XL_zcId,
															   LPXLOPER XL_zdId,
															   LPXLOPER XL_volId,
															   LPXLOPER XL_smileId,
															   LPXLOPER XL_autoMode,
															   LPXLOPER XL_horizon,
															   LPXLOPER XL_fineMonth,
															   LPXLOPER XL_shapeDecay,
															   LPXLOPER XL_shapeSlope,
															   LPXLOPER XL_shapeAsymptote,
															   LPXLOPER XL_nbFactor,
															   LPXLOPER XL_corrMatu,
															   LPXLOPER XL_corrMatrix,
															   LPXLOPER XL_corr2Matu)
{
	ADD_LOG("Local_PXL_FRMTREE_AUTO_B");
	//ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();

	// C variable
	CCString C_zcId;
	CCString C_zdId;

    CCString C_volId;
    CCString C_smileId;
    CCString CS_autoMode;
    CCString CS_fineMonth;
    double C_horizon;
    double C_shapeDecay;
    double C_shapeSlope;
    double C_shapeAsymptote;
    double C_nbFactor;
   
    VECTOR<double> C_corrMatu;
    VECTOR<double> C_corrMatrix;
    VECTOR<double> C_corr2Matu;	

    double mFineDefault = K_L_DEFAULT;
    double epsDefault = 1e-8;
    double huegeDefault = 1e8;
    double zero = 0;
    double one = 1;
	
	// error
	static int error;
	static char* reason = "";

    XL_readStrCell(XL_zcId, C_zcId, " ARM_ERR: Curve id: string expected",C_result);
	XL_readStrCell(XL_zdId, C_zdId, " ARM_ERR: Curve id: string expected",C_result);
    XL_readStrCell(XL_volId, C_volId, " ARM_ERR: Vol id: string expected",C_result);
    XL_readStrCell(XL_smileId, C_smileId, " ARM_ERR: Smile id: string expected",C_result);
    XL_readStrCell(XL_autoMode, CS_autoMode, " ARM_ERR: AutoMode: string expected",C_result);
	XL_readNumCell(XL_horizon,C_horizon," ARM_ERR: Horizon date: date expected",C_result);
    XL_readStrCellWD(XL_fineMonth, CS_fineMonth, "L_DEFAULT", " ARM_ERR: fine month. : string expected",C_result);
    XL_readNumCellWD(XL_shapeDecay, C_shapeDecay, zero, " ARM_ERR: shape decay. : numeric expected",C_result);
    XL_readNumCellWD(XL_shapeSlope, C_shapeSlope, zero, " ARM_ERR: shape slope. : numeric expected",C_result);
    XL_readNumCellWD(XL_shapeAsymptote, C_shapeAsymptote, zero, " ARM_ERR: shape Asymptote. : numeric expected",C_result);
    XL_readNumCellWD(XL_nbFactor, C_nbFactor, one, " ARM_ERR: Nb factor. : numeric expected",C_result);

   	if((XL_corrMatu->xltype == xltypeMissing) || (XL_corrMatu->xltype == xltypeNil))
	{
    }
	else
	{
        XL_readNumVector(XL_corrMatu, C_corrMatu, " ARM_ERR: Correlation maturities: array of numeric expected",C_result);
    }

   	if((XL_corrMatrix->xltype == xltypeMissing) || (XL_corrMatrix->xltype == xltypeNil))
	{
    }
	else
	{
        XL_readNumVector(XL_corrMatrix, C_corrMatrix, " ARM_ERR: Correlation matrix: array of numeric expected",C_result);
    }

   	if((XL_corr2Matu->xltype == xltypeMissing) || (XL_corr2Matu->xltype == xltypeNil))
	{
    }
	else
	{
        XL_readNumVector(XL_corr2Matu, C_corr2Matu, " ARM_ERR: Correlation correlated maturities : array of numeric expected",C_result);
    }

	long autoModeId;
	long fineMonthId;

	if((autoModeId = ARM_ConvAutoMode2 (CS_autoMode, C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	if((fineMonthId = ARM_ConvFineMode (CS_fineMonth, C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	long retCode;
	long objId;

	CCString curClass = LOCAL_FRMTREE_CLASS;
	CCString stringId;

    retCode = ARMLOCAL_FRMTREE_AUTO_B(LocalGetNumObjectId(C_zcId),
								 LocalGetNumObjectId(C_zdId),
								 LocalGetNumObjectId(C_volId),
								 LocalGetNumObjectId(C_smileId),
								 autoModeId,
								 C_horizon,
								 fineMonthId,
								 C_shapeDecay,
								 C_shapeSlope,
								 C_shapeAsymptote,
								 (long) C_nbFactor,
								 C_corrMatu,
								 C_corrMatrix,
								 C_corr2Matu,
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

	//ARM_END();
	
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_PXL_FRMTREE_AUTO_B" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}



__declspec(dllexport) LPXLOPER WINAPI Local_DFFXBS (LPXLOPER XL_dVol,
													LPXLOPER XL_fVol,
													LPXLOPER XL_dZc,
													LPXLOPER XL_fZc,
													LPXLOPER XL_dFxCorr,
													LPXLOPER XL_fFxCorr,
													LPXLOPER XL_fxVol,
													LPXLOPER XL_ratesCorr,
													LPXLOPER XL_dBSZc,
													LPXLOPER XL_fBSZc,
													LPXLOPER XL_spot,
                                                    LPXLOPER XL_FundZc,
													LPXLOPER XL_FundBSZc,
													LPXLOPER XL_fundFxSpot,
													LPXLOPER XL_fxSmile,
                                                    LPXLOPER XL_LnOrNorVols,
													LPXLOPER XL_sabrDom,
													LPXLOPER XL_sabrFor,
													LPXLOPER XL_fxVolModel)
{
	ADD_LOG("Local_DFFXBS ");
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
		ARM_NOCALCIFWIZ();

		// C variable
		CCString C_dVol;
		CCString C_fVol;
		CCString C_dZc;
		CCString C_fZc;
		CCString C_dFxCorr;
		CCString C_fFxCorr;
		CCString C_fxVol;
		
		double C_ratesCorr;

		CCString C_dBSZc;
		long dBSZcId;

		CCString C_fBSZc;
		long fBSZcId;

		double C_spot;
		double C_spot_default = 0.0;

		CCString C_FundZC;
		long FundZCId;

		CCString C_FundBSZC;
		long FundBSZCId;

		CCString C_fxSmile;
		long FxSmileId;

		double C_fundFxSpot;
		double C_fundFxSpot_default = 0.0;

		vector<CCString> C_sabrDom;
		vector<CCString> C_sabrDom_default(0);
		long rhoDomId = ARM_NULL_OBJECT;
		long nuDomId = ARM_NULL_OBJECT;
		long betaDomId = ARM_NULL_OBJECT;

		vector<CCString> C_sabrFor;
		vector<CCString> C_sabrFor_default(0);
		long rhoForId = ARM_NULL_OBJECT;
		long nuForId = ARM_NULL_OBJECT;
		long betaForId = ARM_NULL_OBJECT;

		CCString C_fxVolModel;
		long fxVolModelId;

		// error
		static int error;
		static char* reason = "";

		XL_readStrCell(XL_dVol,C_dVol," ARM_ERR: domestic volatility curve id: object expected",C_result);
		XL_readStrCell(XL_fVol,C_fVol," ARM_ERR: foreign volatility curve id: object expected",C_result);
		XL_readStrCell(XL_dZc,C_dZc," ARM_ERR: domestic zero-coupon curve id: object expected",C_result);
		XL_readStrCell(XL_fZc,C_fZc," ARM_ERR: foreign zero-coupon curve id: object expected",C_result);
		XL_readStrCell(XL_dFxCorr,C_dFxCorr," ARM_ERR: domestic forex correlation: object expected",C_result);
		XL_readStrCell(XL_fFxCorr,C_fFxCorr," ARM_ERR: foreign forex correlation: object expected",C_result);
		XL_readStrCell(XL_fxVol,C_fxVol," ARM_ERR: forex volatility: object expected",C_result);
		XL_readNumCell(XL_ratesCorr,C_ratesCorr," ARM_ERR: rates correlation: numeric expected",C_result);
		
		XL_readStrCellWD(XL_dBSZc,C_dBSZc,"DEFAULT"," ARM_ERR: domestic BS zero-coupon id: object expected",C_result);
		XL_readStrCellWD(XL_fBSZc,C_fBSZc,"DEFAULT"," ARM_ERR: foreign BS zero-coupon curve id: object expected",C_result);
		XL_readNumCellWD(XL_spot,C_spot,C_spot_default," ARM_ERR: FX spot: numeric expected",C_result);

		XL_readStrCellWD(XL_FundZc,C_FundZC,"DEFAULT"," ARM_ERR: Funding zero-coupon id: object expected",C_result);
		XL_readStrCellWD(XL_FundBSZc,C_FundBSZC,"DEFAULT"," ARM_ERR: Funding BS zero-coupon curve id: object expected",C_result);
		XL_readNumCellWD(XL_fundFxSpot,C_fundFxSpot,C_fundFxSpot_default," ARM_ERR: Funding FX spot: numeric expected",C_result);
		XL_readStrCellWD(XL_fxSmile,C_fxSmile,"DEFAULT"," ARM_ERR: forex smile: object expected",C_result);

		CCString lnOrNorVols;
		XL_readStrCellWD(XL_LnOrNorVols,lnOrNorVols,"Y"," ARM_ERR: LogNor or Nor Vol: array of string expected (Y/N)",C_result);
		bool isLnVols=true;
		lnOrNorVols.toUpper();
		if(CCSTringToSTLString(lnOrNorVols)!="Y")
				isLnVols=false;

		XL_readStrVectorWD(XL_sabrDom,C_sabrDom,C_sabrDom_default," ARM_ERR: Sabr domestic: vector of objects expected",DOUBLE_TYPE,C_result);	
		XL_readStrVectorWD(XL_sabrFor,C_sabrFor,C_sabrFor_default," ARM_ERR: Sabr foreign: vector of objects expected",DOUBLE_TYPE,C_result);	

		XL_readStrCellWD(XL_fxVolModel,C_fxVolModel,"DEFAULT"," ARM_ERR: fx vol model: object expected",C_result);

		if ( C_sabrDom.size() == 1 )
		{
			rhoDomId = LocalGetNumObjectId(C_sabrDom[0]);
		}
		else if ( C_sabrDom.size() == 2 )
		{
			rhoDomId = LocalGetNumObjectId(C_sabrDom[0]);
			nuDomId = LocalGetNumObjectId(C_sabrDom[1]);
		}
		else if ( C_sabrDom.size() == 3 )
		{
			rhoDomId = LocalGetNumObjectId(C_sabrDom[0]);
			nuDomId = LocalGetNumObjectId(C_sabrDom[1]);
			betaDomId = LocalGetNumObjectId(C_sabrDom[2]);
		}
		else if ( C_sabrDom.size() > 3 )
			throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT, "vector of Sabr : Rho, Nu, Beta ");

		if ( C_sabrFor.size() == 1 )
		{
			rhoForId = LocalGetNumObjectId(C_sabrFor[0]);
		}
		else if ( C_sabrFor.size() == 2 )
		{
			rhoForId = LocalGetNumObjectId(C_sabrFor[0]);
			nuForId = LocalGetNumObjectId(C_sabrFor[1]);
		}
		else if ( C_sabrFor.size() == 3 )
		{
			rhoForId = LocalGetNumObjectId(C_sabrFor[0]);
			nuForId = LocalGetNumObjectId(C_sabrFor[1]);
			betaForId = LocalGetNumObjectId(C_sabrFor[2]);
		}
		else if ( C_sabrFor.size() > 3 )
			throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT, "vector of Sabr : Rho, Nu, Beta ");

		if(C_dBSZc == "DEFAULT")
		{
			dBSZcId = ARM_NULL_OBJECT;
		}
		else
		{
			dBSZcId = LocalGetNumObjectId (C_dBSZc);
		}

		if(C_fBSZc == "DEFAULT")
		{
			fBSZcId = ARM_NULL_OBJECT;
		}
		else
		{
			fBSZcId = LocalGetNumObjectId (C_fBSZc);
		}

		if (C_FundZC == "DEFAULT")
		{
		   FundZCId = ARM_NULL_OBJECT;
		}
		else
		{
		   FundZCId = LocalGetNumObjectId (C_FundZC);
		}

		if ( C_FundBSZC == "DEFAULT" )
		{
		   FundBSZCId = ARM_NULL_OBJECT;
		}
		else
		{
		   FundBSZCId = LocalGetNumObjectId(C_FundBSZC);
		}

		if(C_fxSmile == "DEFAULT")
		{
			FxSmileId = ARM_NULL_OBJECT;
		}
		else
		{
			FxSmileId = LocalGetNumObjectId (C_fxSmile);
		}

		if(C_fxVolModel == "DEFAULT")
		{
			fxVolModelId = ARM_NULL_OBJECT;
		}
		else
		{
			fxVolModelId = LocalGetNumObjectId (C_fxVolModel);
		}
		
		long retCode;
		long objId;
		CCString prevClass;

		CCString curClass = LOCAL_DFFXBS_CLASS;
		CCString stringId = GetLastCurCellEnvValue ();

		if (!stringId)
		{
			retCode = ARMLOCAL_DFFXBS (LocalGetNumObjectId (C_dVol),
									   LocalGetNumObjectId (C_fVol),
									   LocalGetNumObjectId (C_dZc),
									   LocalGetNumObjectId (C_fZc),
									   LocalGetNumObjectId (C_dFxCorr),
									   LocalGetNumObjectId (C_fFxCorr),
									   LocalGetNumObjectId (C_fxVol),
									   C_ratesCorr,
									   dBSZcId,
									   fBSZcId,
									   C_spot,
									   FundZCId,
									   FundBSZCId,
									   C_fundFxSpot,
									   FxSmileId,
									   isLnVols,
									   rhoDomId,
									   nuDomId,
									   betaDomId,
									   rhoForId,
									   nuForId,
									   betaForId,
									   fxVolModelId,
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
				
			if(curClass == prevClass)
			{
				retCode = ARMLOCAL_DFFXBS (LocalGetNumObjectId (C_dVol),
										   LocalGetNumObjectId (C_fVol),
										   LocalGetNumObjectId (C_dZc),
										   LocalGetNumObjectId (C_fZc),
										   LocalGetNumObjectId (C_dFxCorr),
										   LocalGetNumObjectId (C_fFxCorr),
										   LocalGetNumObjectId (C_fxVol),
										   C_ratesCorr,
										   dBSZcId,
										   fBSZcId,
										   C_spot,
										   FundZCId,
										   FundBSZCId,
										   C_fundFxSpot,
										   FxSmileId,
										   isLnVols,
										   rhoDomId,
										   nuDomId,
										   betaDomId,
										   rhoForId,
										   nuForId,
										   betaForId,
										   fxVolModelId,
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
				retCode = ARMLOCAL_DFFXBS (LocalGetNumObjectId (C_dVol),
										   LocalGetNumObjectId (C_fVol),
										   LocalGetNumObjectId (C_dZc),
										   LocalGetNumObjectId (C_fZc),
										   LocalGetNumObjectId (C_dFxCorr),
										   LocalGetNumObjectId (C_fFxCorr),
										   LocalGetNumObjectId (C_fxVol),
										   C_ratesCorr,
										   dBSZcId,
										   fBSZcId,
										   C_spot,
										   FundZCId,
										   FundBSZCId,
										   C_fundFxSpot,
										   FxSmileId,
										   isLnVols,
										   rhoDomId,
										   nuDomId,
										   betaDomId,
										   rhoForId,
										   nuForId,
										   betaForId,
										   fxVolModelId,
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

	//	ARM_END();
	
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_DFFXBS" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}



__declspec(dllexport) LPXLOPER WINAPI Local_PXL_DFFXBS (LPXLOPER XL_dVol,
														LPXLOPER XL_fVol,
														LPXLOPER XL_dZc,
														LPXLOPER XL_fZc,
														LPXLOPER XL_dFxCorr,
														LPXLOPER XL_fFxCorr,
														LPXLOPER XL_fxVol,
														LPXLOPER XL_ratesCorr,
														LPXLOPER XL_dBSZc,
														LPXLOPER XL_fBSZc,
														LPXLOPER XL_spot,
                                                        LPXLOPER XL_FundZC,
													    LPXLOPER XL_FundBSZC,
													    LPXLOPER XL_fundFxSpot,
														LPXLOPER XL_fxSmile,
                                                        LPXLOPER XL_LnOrNorVols,
														LPXLOPER XL_sabrDom,
														LPXLOPER XL_sabrFor,
														LPXLOPER XL_fxVolModel)
{
	ADD_LOG("Local_PXL_DFFXBS ");
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
		ARM_NOCALCIFWIZ();

		// C variable
		CCString C_dVol;
		CCString C_fVol;
		CCString C_dZc;
		CCString C_fZc;
		CCString C_dFxCorr;
		CCString C_fFxCorr;
		CCString C_fxVol;
		double C_ratesCorr;
		
		CCString C_dBSZc;
		long dBSZcId;

		CCString C_fBSZc;
		long fBSZcId;

		double C_spot;
		double C_spot_default;
		
		CCString C_FundZC;
		long FundZCId;

		CCString C_FundBSZC;
		long FundBSZCId;

		CCString C_fxSmile;
		long FxSmileId;

		double C_fundFxSpot;
		double C_fundFxSpot_default = 0.0;  
    
		vector<CCString> C_sabrDom;
		vector<CCString> C_sabrDom_default(0);
		long rhoDomId = ARM_NULL_OBJECT;
		long nuDomId = ARM_NULL_OBJECT;
		long betaDomId = ARM_NULL_OBJECT;

		vector<CCString> C_sabrFor;
		vector<CCString> C_sabrFor_default(0);
		long rhoForId = ARM_NULL_OBJECT;
		long nuForId = ARM_NULL_OBJECT;
		long betaForId = ARM_NULL_OBJECT;

		CCString C_fxVolModel;
		long fxVolModelId;

		// error
		static int error;
		static char* reason = "";

		XL_readStrCell(XL_dVol,C_dVol," ARM_ERR: domestic volatility curve id: object expected",C_result);
		XL_readStrCell(XL_fVol,C_fVol," ARM_ERR: foreign volatility curve id: object expected",C_result);
		XL_readStrCell(XL_dZc,C_dZc," ARM_ERR: domestic zero-coupon curve id: object expected",C_result);
		XL_readStrCell(XL_fZc,C_fZc," ARM_ERR: foreign zero-coupon curve id: object expected",C_result);
		XL_readStrCell(XL_dFxCorr,C_dFxCorr," ARM_ERR: domestic forex correlation: object expected",C_result);
		XL_readStrCell(XL_fFxCorr,C_fFxCorr," ARM_ERR: foreign forex correlation: object expected",C_result);
		XL_readStrCell(XL_fxVol,C_fxVol," ARM_ERR: forex volatility: object expected",C_result);
		XL_readNumCell(XL_ratesCorr,C_ratesCorr," ARM_ERR: rates correlation: numeric expected",C_result);
		
		XL_readStrCellWD(XL_dBSZc,C_dBSZc,"DEFAULT"," ARM_ERR: domestic BS zero-coupon id: object expected",C_result);
		XL_readStrCellWD(XL_fBSZc,C_fBSZc,"DEFAULT"," ARM_ERR: foreign BS zero-coupon curve id: object expected",C_result);
		XL_readNumCellWD(XL_spot,C_spot,C_spot_default," ARM_ERR: FX spot: numeric expected",C_result);
		
		XL_readStrCellWD(XL_FundZC,C_FundZC,"DEFAULT"," ARM_ERR: Funding zero-coupon id: object expected",C_result);
		XL_readStrCellWD(XL_FundBSZC,C_FundBSZC,"DEFAULT"," ARM_ERR: BS Funding zero-coupon curve id: object expected",C_result);
		XL_readNumCellWD(XL_fundFxSpot,C_fundFxSpot,C_fundFxSpot_default," ARM_ERR: funding FX spot: numeric expected",C_result);
		XL_readStrCellWD(XL_fxSmile,C_fxSmile,"DEFAULT"," ARM_ERR: forex smile: object expected",C_result);

		XL_readNumCellWD(XL_fundFxSpot,C_fundFxSpot,C_fundFxSpot_default," ARM_ERR: Funding FX spot: numeric expected",C_result);
		XL_readStrCellWD(XL_fxSmile,C_fxSmile,"DEFAULT"," ARM_ERR: forex smile: object expected",C_result);

		CCString lnOrNorVols;
		XL_readStrCellWD(XL_LnOrNorVols,lnOrNorVols,"Y"," ARM_ERR: LogNor or Nor Vol: array of string expected (Y/N)",C_result);
		bool isLnVols=true;
		lnOrNorVols.toUpper();
		if(CCSTringToSTLString(lnOrNorVols)!="Y")
				isLnVols=false;

		XL_readStrVectorWD(XL_sabrDom,C_sabrDom,C_sabrDom_default," ARM_ERR: Sabr domestic: vector of objects expected",DOUBLE_TYPE,C_result);	
		XL_readStrVectorWD(XL_sabrFor,C_sabrFor,C_sabrFor_default," ARM_ERR: Sabr foreign: vector of objects expected",DOUBLE_TYPE,C_result);	

		XL_readStrCellWD(XL_fxVolModel,C_fxVolModel,"DEFAULT"," ARM_ERR: fx vol model: object expected",C_result);

		if ( C_sabrDom.size() == 1 )
		{
			rhoDomId = LocalGetNumObjectId(C_sabrDom[0]);
		}
		else if ( C_sabrDom.size() == 2 )
		{
			rhoDomId = LocalGetNumObjectId(C_sabrDom[0]);
			nuDomId = LocalGetNumObjectId(C_sabrDom[1]);
		}
		else if ( C_sabrDom.size() == 3 )
		{
			rhoDomId = LocalGetNumObjectId(C_sabrDom[0]);
			nuDomId = LocalGetNumObjectId(C_sabrDom[1]);
			betaDomId = LocalGetNumObjectId(C_sabrDom[2]);
		}
		else if ( C_sabrDom.size() > 3 )
			throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT, "vector of Sabr : Rho, Nu, Beta ");

		if ( C_sabrFor.size() == 1 )
		{
			rhoForId = LocalGetNumObjectId(C_sabrFor[0]);
		}
		else if ( C_sabrFor.size() == 2 )
		{
			rhoForId = LocalGetNumObjectId(C_sabrFor[0]);
			nuForId = LocalGetNumObjectId(C_sabrFor[1]);
		}
		else if ( C_sabrFor.size() == 3 )
		{
			rhoForId = LocalGetNumObjectId(C_sabrFor[0]);
			nuForId = LocalGetNumObjectId(C_sabrFor[1]);
			betaForId = LocalGetNumObjectId(C_sabrFor[2]);
		}
		else if ( C_sabrFor.size() > 3 )
			throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT, "vector of Sabr : Rho, Nu, Beta ");

		if(C_dBSZc == "DEFAULT")
		{
			dBSZcId = ARM_NULL_OBJECT;
		}
		else
		{
			dBSZcId = LocalGetNumObjectId (C_dBSZc);
		}

		if(C_fBSZc == "DEFAULT")
		{
			fBSZcId = ARM_NULL_OBJECT;
		}
		else
		{
			fBSZcId = LocalGetNumObjectId (C_fBSZc);
		}

		if (C_FundZC == "DEFAULT")
		{
		   FundZCId = ARM_NULL_OBJECT;
		}
		else
		{
		   FundZCId = LocalGetNumObjectId (C_FundZC);
		}

		if ( C_FundBSZC == "DEFAULT" )
		{
		   FundBSZCId = ARM_NULL_OBJECT;
		}
		else
		{
		   FundBSZCId = LocalGetNumObjectId(C_FundBSZC);
		}

		if(C_fxSmile == "DEFAULT")
		{
			FxSmileId = ARM_NULL_OBJECT;
		}
		else
		{
			FxSmileId = LocalGetNumObjectId (C_fxSmile);
		}

		if(C_fxVolModel == "DEFAULT")
		{
			fxVolModelId = ARM_NULL_OBJECT;
		}
		else
		{
			fxVolModelId = LocalGetNumObjectId (C_fxVolModel);
		}
		
		long retCode;
		long objId;
		
		CCString curClass = LOCAL_DFFXBS_CLASS;
		CCString stringId;

		retCode = ARMLOCAL_DFFXBS (LocalGetNumObjectId (C_dVol),
								   LocalGetNumObjectId (C_fVol),
								   LocalGetNumObjectId (C_dZc),
								   LocalGetNumObjectId (C_fZc),
								   LocalGetNumObjectId (C_dFxCorr),
								   LocalGetNumObjectId (C_fFxCorr),
								   LocalGetNumObjectId (C_fxVol),
								   C_ratesCorr,
								   dBSZcId,
								   fBSZcId,
								   C_spot,
								   FundZCId,
								   FundBSZCId,
								   C_fundFxSpot,
								   FxSmileId,
								   isLnVols,
								   rhoDomId,
								   nuDomId,
								   betaDomId,
								   rhoForId,
								   nuForId,
								   betaForId,
								   fxVolModelId,
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
		   PXL_ARM_ERR();
		}

	//	ARM_END();
	
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_PXL_DFFXBS" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


__declspec(dllexport) LPXLOPER WINAPI Local_DFFXBS_from2BS (LPXLOPER XL_dBS,
															LPXLOPER XL_fBS,
															LPXLOPER XL_dFxCorr,
															LPXLOPER XL_fFxCorr,
															LPXLOPER XL_fxVol,
															LPXLOPER XL_ratesCorr,
															LPXLOPER XL_dBSZc,
															LPXLOPER XL_fBSZc,
															LPXLOPER XL_spot,
															LPXLOPER XL_FundZc,
															LPXLOPER XL_FundBSZc,
															LPXLOPER XL_fundFxSpot,
															LPXLOPER XL_fxSmile,
															LPXLOPER XL_fxVolModel)
{
	ADD_LOG("Local_DFFXBS_from2BS ");
	static XLOPER XL_result;
	ARM_result C_result;

	ARM_XL_TRY_BLOCK_BEGIN
	{
		ARM_NOCALCIFWIZ();

		// C variable
		CCString C_dBS;
		CCString C_fBS;
		CCString C_dFxCorr;
		CCString C_fFxCorr;
		CCString C_fxVol;
		
		double C_ratesCorr;

		CCString C_dBSZc;
		long dBSZcId;

		CCString C_fBSZc;
		long fBSZcId;

		double C_spot;
		double C_spot_default = 0.0;

		CCString C_FundZC;
		long FundZCId;

		CCString C_FundBSZC;
		long FundBSZCId;

		CCString C_fxSmile;
		long FxSmileId;

		double C_fundFxSpot;
		double C_fundFxSpot_default = 0.0;

		CCString C_fxVolModel;
		long fxVolModelId;

		// error
		static int error;
		static char* reason = "";

		XL_readStrCell(XL_dBS,C_dBS," ARM_ERR: domestic BS model id: object expected",C_result);
		XL_readStrCell(XL_fBS,C_fBS," ARM_ERR: foreign BS model id: object expected",C_result);
		XL_readStrCell(XL_dFxCorr,C_dFxCorr," ARM_ERR: domestic forex correlation: object expected",C_result);
		XL_readStrCell(XL_fFxCorr,C_fFxCorr," ARM_ERR: foreign forex correlation: object expected",C_result);
		XL_readStrCell(XL_fxVol,C_fxVol," ARM_ERR: forex volatility: object expected",C_result);
		XL_readNumCell(XL_ratesCorr,C_ratesCorr," ARM_ERR: rates correlation: numeric expected",C_result);
		
		XL_readStrCellWD(XL_dBSZc,C_dBSZc,"DEFAULT"," ARM_ERR: domestic BS zero-coupon id: object expected",C_result);
		XL_readStrCellWD(XL_fBSZc,C_fBSZc,"DEFAULT"," ARM_ERR: foreign BS zero-coupon curve id: object expected",C_result);
		XL_readNumCellWD(XL_spot,C_spot,C_spot_default," ARM_ERR: FX spot: numeric expected",C_result);

		XL_readStrCellWD(XL_FundZc,C_FundZC,"DEFAULT"," ARM_ERR: Funding zero-coupon id: object expected",C_result);
		XL_readStrCellWD(XL_FundBSZc,C_FundBSZC,"DEFAULT"," ARM_ERR: Funding BS zero-coupon curve id: object expected",C_result);
		XL_readNumCellWD(XL_fundFxSpot,C_fundFxSpot,C_fundFxSpot_default," ARM_ERR: Funding FX spot: numeric expected",C_result);
		XL_readStrCellWD(XL_fxSmile,C_fxSmile,"DEFAULT"," ARM_ERR: forex smile: object expected",C_result);

		XL_readStrCellWD(XL_fxVolModel,C_fxVolModel,"DEFAULT"," ARM_ERR: fx vol model: object expected",C_result);

		if(C_dBSZc == "DEFAULT")
		{
			dBSZcId = ARM_NULL_OBJECT;
		}
		else
		{
			dBSZcId = LocalGetNumObjectId (C_dBSZc);
		}

		if(C_fBSZc == "DEFAULT")
		{
			fBSZcId = ARM_NULL_OBJECT;
		}
		else
		{
			fBSZcId = LocalGetNumObjectId (C_fBSZc);
		}

		if (C_FundZC == "DEFAULT")
		{
		   FundZCId = ARM_NULL_OBJECT;
		}
		else
		{
		   FundZCId = LocalGetNumObjectId (C_FundZC);
		}

		if ( C_FundBSZC == "DEFAULT" )
		{
		   FundBSZCId = ARM_NULL_OBJECT;
		}
		else
		{
		   FundBSZCId = LocalGetNumObjectId(C_FundBSZC);
		}

		if(C_fxSmile == "DEFAULT")
		{
			FxSmileId = ARM_NULL_OBJECT;
		}
		else
		{
			FxSmileId = LocalGetNumObjectId (C_fxSmile);
		}

		if(C_fxVolModel == "DEFAULT")
		{
			fxVolModelId = ARM_NULL_OBJECT;
		}
		else
		{
			fxVolModelId = LocalGetNumObjectId (C_fxVolModel);
		}
		
		long retCode;
		long objId;
		CCString prevClass;

		CCString curClass = LOCAL_DFFXBS_CLASS;
		CCString stringId = GetLastCurCellEnvValue ();

		if (!stringId)
		{
			retCode = ARMLOCAL_DFFXBS (LocalGetNumObjectId (C_dBS),
									   LocalGetNumObjectId (C_fBS),
									   LocalGetNumObjectId (C_dFxCorr),
									   LocalGetNumObjectId (C_fFxCorr),
									   LocalGetNumObjectId (C_fxVol),
									   C_ratesCorr,
									   dBSZcId,
									   fBSZcId,
									   C_spot,
									   FundZCId,
									   FundBSZCId,
									   C_fundFxSpot,
									   FxSmileId,
									   fxVolModelId,
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
				
			if(curClass == prevClass)
			{
				retCode = ARMLOCAL_DFFXBS (LocalGetNumObjectId (C_dBS),
										   LocalGetNumObjectId (C_fBS),
										   LocalGetNumObjectId (C_dFxCorr),
										   LocalGetNumObjectId (C_fFxCorr),
										   LocalGetNumObjectId (C_fxVol),
										   C_ratesCorr,
										   dBSZcId,
										   fBSZcId,
										   C_spot,
										   FundZCId,
										   FundBSZCId,
										   C_fundFxSpot,
										   FxSmileId,
										   fxVolModelId,
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
				retCode = ARMLOCAL_DFFXBS (LocalGetNumObjectId (C_dBS),
										   LocalGetNumObjectId (C_fBS),
										   LocalGetNumObjectId (C_dFxCorr),
										   LocalGetNumObjectId (C_fFxCorr),
										   LocalGetNumObjectId (C_fxVol),
										   C_ratesCorr,
										   dBSZcId,
										   fBSZcId,
										   C_spot,
										   FundZCId,
										   FundBSZCId,
										   C_fundFxSpot,
										   FxSmileId,
										   fxVolModelId,
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
	ARM_XL_TRY_BLOCK_END

	ARM_XL_CATCH_ARM_EXPT

	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_DFFXBS_from2BS" )

	return (LPXLOPER)&XL_result;
}

__declspec(dllexport) LPXLOPER WINAPI Local_PXL_DFFXBS_from2BS (LPXLOPER XL_dBS,
																LPXLOPER XL_fBS,
																LPXLOPER XL_dFxCorr,
																LPXLOPER XL_fFxCorr,
																LPXLOPER XL_fxVol,
																LPXLOPER XL_ratesCorr,
																LPXLOPER XL_dBSZc,
																LPXLOPER XL_fBSZc,
																LPXLOPER XL_spot,
																LPXLOPER XL_FundZc,
																LPXLOPER XL_FundBSZc,
																LPXLOPER XL_fundFxSpot,
																LPXLOPER XL_fxSmile,
																LPXLOPER XL_fxVolModel)
{
	ADD_LOG("Local_PXL_DFFXBS_from2BS ");
	static XLOPER XL_result;
	ARM_result C_result;

	ARM_XL_TRY_BLOCK_BEGIN
	{
		ARM_NOCALCIFWIZ();

		// C variable
		CCString C_dBS;
		CCString C_fBS;
		CCString C_dFxCorr;
		CCString C_fFxCorr;
		CCString C_fxVol;
		double C_ratesCorr;
		
		CCString C_dBSZc;
		long dBSZcId;

		CCString C_fBSZc;
		long fBSZcId;

		double C_spot;
		double C_spot_default;
		
		CCString C_FundZC;
		long FundZCId;

		CCString C_FundBSZC;
		long FundBSZCId;

		CCString C_fxSmile;
		long FxSmileId;

		double C_fundFxSpot;
		double C_fundFxSpot_default = 0.0;  

		CCString C_fxVolModel;
		long fxVolModelId;

		// error
		static int error;
		static char* reason = "";

		XL_readStrCell(XL_dBS,C_dBS," ARM_ERR: domestic BS model id: object expected",C_result);
		XL_readStrCell(XL_fBS,C_fBS," ARM_ERR: foreign BS model id: object expected",C_result);
		XL_readStrCell(XL_dFxCorr,C_dFxCorr," ARM_ERR: domestic forex correlation: object expected",C_result);
		XL_readStrCell(XL_fFxCorr,C_fFxCorr," ARM_ERR: foreign forex correlation: object expected",C_result);
		XL_readStrCell(XL_fxVol,C_fxVol," ARM_ERR: forex volatility: object expected",C_result);
		XL_readNumCell(XL_ratesCorr,C_ratesCorr," ARM_ERR: rates correlation: numeric expected",C_result);
		
		XL_readStrCellWD(XL_dBSZc,C_dBSZc,"DEFAULT"," ARM_ERR: domestic BS zero-coupon id: object expected",C_result);
		XL_readStrCellWD(XL_fBSZc,C_fBSZc,"DEFAULT"," ARM_ERR: foreign BS zero-coupon curve id: object expected",C_result);
		XL_readNumCellWD(XL_spot,C_spot,C_spot_default," ARM_ERR: FX spot: numeric expected",C_result);
		
		XL_readStrCellWD(XL_FundZc,C_FundZC,"DEFAULT"," ARM_ERR: Funding zero-coupon id: object expected",C_result);
		XL_readStrCellWD(XL_FundBSZc,C_FundBSZC,"DEFAULT"," ARM_ERR: BS Funding zero-coupon curve id: object expected",C_result);
		XL_readNumCellWD(XL_fundFxSpot,C_fundFxSpot,C_fundFxSpot_default," ARM_ERR: funding FX spot: numeric expected",C_result);
		XL_readStrCellWD(XL_fxSmile,C_fxSmile,"DEFAULT"," ARM_ERR: forex smile: object expected",C_result);

		XL_readNumCellWD(XL_fundFxSpot,C_fundFxSpot,C_fundFxSpot_default," ARM_ERR: Funding FX spot: numeric expected",C_result);
		XL_readStrCellWD(XL_fxSmile,C_fxSmile,"DEFAULT"," ARM_ERR: forex smile: object expected",C_result);

		XL_readStrCellWD(XL_fxVolModel,C_fxVolModel,"DEFAULT"," ARM_ERR: fx vol model: object expected",C_result);

		if(C_dBSZc == "DEFAULT")
		{
			dBSZcId = ARM_NULL_OBJECT;
		}
		else
		{
			dBSZcId = LocalGetNumObjectId (C_dBSZc);
		}

		if(C_fBSZc == "DEFAULT")
		{
			fBSZcId = ARM_NULL_OBJECT;
		}
		else
		{
			fBSZcId = LocalGetNumObjectId (C_fBSZc);
		}

		if (C_FundZC == "DEFAULT")
		{
		   FundZCId = ARM_NULL_OBJECT;
		}
		else
		{
		   FundZCId = LocalGetNumObjectId (C_FundZC);
		}

		if ( C_FundBSZC == "DEFAULT" )
		{
		   FundBSZCId = ARM_NULL_OBJECT;
		}
		else
		{
		   FundBSZCId = LocalGetNumObjectId(C_FundBSZC);
		}

		if(C_fxSmile == "DEFAULT")
		{
			FxSmileId = ARM_NULL_OBJECT;
		}
		else
		{
			FxSmileId = LocalGetNumObjectId (C_fxSmile);
		}

		if(C_fxVolModel == "DEFAULT")
		{
			fxVolModelId = ARM_NULL_OBJECT;
		}
		else
		{
			fxVolModelId = LocalGetNumObjectId (C_fxVolModel);
		}
		
		long retCode;
		long objId;
		
		CCString curClass = LOCAL_DFFXBS_CLASS;
		CCString stringId;

		retCode = ARMLOCAL_DFFXBS (LocalGetNumObjectId (C_dBS),
								   LocalGetNumObjectId (C_fBS),
								   LocalGetNumObjectId (C_dFxCorr),
								   LocalGetNumObjectId (C_fFxCorr),
								   LocalGetNumObjectId (C_fxVol),
								   C_ratesCorr,
								   dBSZcId,
								   fBSZcId,
								   C_spot,
								   FundZCId,
								   FundBSZCId,
								   C_fundFxSpot,
								   FxSmileId,
								   fxVolModelId,
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
		   PXL_ARM_ERR();
		}
	}
	ARM_XL_TRY_BLOCK_END

	ARM_XL_CATCH_ARM_EXPT

	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_PXL_DFFXBS_from2BS" )

	return (LPXLOPER)&XL_result;
}

__declspec(dllexport) LPXLOPER WINAPI Local_ARM_GetComputedSigmaSABR(LPXLOPER XL_modId)
{
	ADD_LOG("Local_ARM_GetComputedSigmaSABR");
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();

	// C variable
	CCString C_modId;

	// error
	static int error;
	static char* reason = "";

    XL_readStrCell(XL_modId, C_modId, " ARM_ERR: model id: object expected",C_result);
	    
	long retCode = ARMLOCAL_GetComputedSigmaSABR(LocalGetNumObjectId (C_modId),
												C_result);

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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_ARM_GetComputedSigmaSABR" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}



__declspec(dllexport) LPXLOPER WINAPI Local_ARM_CALIBRATIONHWSV(LPXLOPER XL_asof,
																LPXLOPER XL_curve,
																LPXLOPER XL_vol,
																LPXLOPER XL_sec,
																LPXLOPER XL_amin,
																LPXLOPER XL_amax,
																LPXLOPER XL_volmin,
																LPXLOPER XL_volmax,
																LPXLOPER XL_dates,
																LPXLOPER XL_pf,
																LPXLOPER XL_rhoswopt,
																LPXLOPER XL_nuswopt,
																LPXLOPER XL_betaswopt)
{
	ADD_LOG("Local_ARM_CALIBRATIONHWSV");
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();

	// C variable
	double C_asof;

	CCString C_Zc;
	CCString C_Vol;
	CCString C_sec;

	double C_amin;
	double C_amin_default = 0.001;

	double C_amax;
	double C_amax_default = 0.5;

	double C_volmin;
	double C_volmin_default = 0.01;

	double C_volmax;
	double C_volmax_default = 2;

	VECTOR<double> C_dates;

	CCString C_pf;
	long pfId;

	CCString C_RhoSwopt;
	long RhoSwoptId;
	CCString C_NuSwopt;
	long NuSwoptId;
	CCString C_BetaSwopt;
	long BetaSwoptId;


	// error
	static int error;
	static char* reason = "";

	XL_readNumCell(XL_asof,C_asof," ARM_ERR: as of date: numeric expected",C_result);
	XL_readStrCell(XL_curve,C_Zc," ARM_ERR: curve id: object expected",C_result);
	XL_readStrCell(XL_vol,C_Vol," ARM_ERR: volcurve id: object expected",C_result);
	XL_readStrCell(XL_sec,C_sec," ARM_ERR: security id: object expected",C_result);
	XL_readNumCellWD(XL_amin,C_amin,C_amin_default," ARM_ERR: a min: numeric expected",C_result);
	XL_readNumCellWD(XL_amax,C_amax,C_amax_default," ARM_ERR: a max: numeric expected",C_result);
	XL_readNumCellWD(XL_volmin,C_volmin,C_volmin_default," ARM_ERR: vol min: numeric expected",C_result);
	XL_readNumCellWD(XL_volmax,C_volmax,C_volmax_default," ARM_ERR: vol max: numeric expected",C_result);

   	if((XL_dates->xltype == xltypeMissing) || (XL_dates->xltype == xltypeNil))
	{
    }
	else
	{
        XL_readNumVector(XL_dates, C_dates, " ARM_ERR: dates: array of numeric expected",C_result);
    }

	XL_readStrCellWD(XL_pf,C_pf,"DEFAULT"," ARM_ERR: portfolio id: object expected",C_result);

	if(C_pf == "DEFAULT")
	{
		pfId = ARM_NULL_OBJECT;
	}
	else
	{
		pfId = LocalGetNumObjectId (C_pf);
	}

	XL_readStrCellWD(XL_rhoswopt,C_RhoSwopt,"DEFAULT"," ARM_ERR: rhoswopt id: object expected",C_result);

	if(C_RhoSwopt == "DEFAULT")
	{
		RhoSwoptId = ARM_NULL_OBJECT;
	}
	else
	{
		RhoSwoptId = LocalGetNumObjectId (C_RhoSwopt);
	}

	XL_readStrCellWD(XL_nuswopt,C_NuSwopt,"DEFAULT"," ARM_ERR: nuswopt id: object expected",C_result);

	if(C_NuSwopt == "DEFAULT")
	{
		NuSwoptId = ARM_NULL_OBJECT;
	}
	else
	{
		NuSwoptId = LocalGetNumObjectId (C_NuSwopt);
	}

	XL_readStrCellWD(XL_betaswopt,C_BetaSwopt,"DEFAULT"," ARM_ERR: betaswopt id: object expected",C_result);

	if(C_BetaSwopt == "DEFAULT")
	{
		BetaSwoptId = ARM_NULL_OBJECT;
	}
	else
	{
		BetaSwoptId = LocalGetNumObjectId (C_BetaSwopt);
	}

	long retCode;
	long objId;
	CCString prevClass;
	
	CCString curClass = LOCAL_CALIB_HWSV_CLASS;
	CCString stringId = GetLastCurCellEnvValue ();

	if (!stringId)
	{
		retCode =  ARMLOCAL_CALIBRATIONHWSV(C_asof,
											LocalGetNumObjectId (C_Zc),
											LocalGetNumObjectId (C_Vol),
											LocalGetNumObjectId (C_sec),
											C_amin,
											C_amax,
											C_volmin,
											C_volmax,
											C_dates,
											pfId,
											RhoSwoptId,
											NuSwoptId,
											BetaSwoptId,
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
			
		if(curClass == prevClass)
		{
			retCode = ARMLOCAL_CALIBRATIONHWSV (C_asof,
												LocalGetNumObjectId (C_Zc),
												LocalGetNumObjectId (C_Vol),
												LocalGetNumObjectId (C_sec),
												C_amin,
												C_amax,
												C_volmin,
												C_volmax,
												C_dates,
												pfId,
												RhoSwoptId,
												NuSwoptId,
												BetaSwoptId,
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

			retCode = ARMLOCAL_CALIBRATIONHWSV (C_asof,
												LocalGetNumObjectId (C_Zc),
												LocalGetNumObjectId (C_Vol),
												LocalGetNumObjectId (C_sec),
												C_amin,
												C_amax,
												C_volmin,
												C_volmax,
												C_dates,
												pfId,
												RhoSwoptId,
												NuSwoptId,
												BetaSwoptId,
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

//	ARM_END();
	
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_ARM_CALIBRATIONHWSV" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


__declspec(dllexport) LPXLOPER WINAPI Local_PXL_ARM_CALIBRATIONHWSV(LPXLOPER XL_asof,
																	LPXLOPER XL_curve,
																	LPXLOPER XL_vol,
																	LPXLOPER XL_sec,
																	LPXLOPER XL_amin,
																	LPXLOPER XL_amax,
																	LPXLOPER XL_volmin,
																	LPXLOPER XL_volmax,
																	LPXLOPER XL_dates,
																	LPXLOPER XL_pf,
																	LPXLOPER XL_rhoswopt,
																	LPXLOPER XL_nuswopt,
																	LPXLOPER XL_betaswopt)
{
	ADD_LOG("Local_PXL_ARM_CALIBRATIONHWSV");
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();

	// C variable
	double C_asof;

	CCString C_Zc;
	CCString C_Vol;
	CCString C_sec;

	double C_amin;
	double C_amin_default = 0.001;

	double C_amax;
	double C_amax_default = 0.5;

	double C_volmin;
	double C_volmin_default = 0.01;

	double C_volmax;
	double C_volmax_default = 2;

	VECTOR<double> C_dates;

	CCString C_pf;
	long pfId;

	CCString C_RhoSwopt;
	long RhoSwoptId;
	CCString C_NuSwopt;
	long NuSwoptId;
	CCString C_BetaSwopt;
	long BetaSwoptId;

	// error
	static int error;
	static char* reason = "";

	XL_readNumCell(XL_asof,C_asof," ARM_ERR: as of date: numeric expected",C_result);
	XL_readStrCell(XL_curve,C_Zc," ARM_ERR: curve id: object expected",C_result);
	XL_readStrCell(XL_vol,C_Vol," ARM_ERR: volcurve id: object expected",C_result);
	XL_readStrCell(XL_sec,C_sec," ARM_ERR: security id: object expected",C_result);
	XL_readNumCellWD(XL_amin,C_amin,C_amin_default," ARM_ERR: a min: numeric expected",C_result);
	XL_readNumCellWD(XL_amax,C_amax,C_amax_default," ARM_ERR: a max: numeric expected",C_result);
	XL_readNumCellWD(XL_volmin,C_volmin,C_volmin_default," ARM_ERR: vol min: numeric expected",C_result);
	XL_readNumCellWD(XL_volmax,C_volmax,C_volmax_default," ARM_ERR: vol max: numeric expected",C_result);

   	if((XL_dates->xltype == xltypeMissing) || (XL_dates->xltype == xltypeNil))
	{
    }
	else
	{
        XL_readNumVector(XL_dates, C_dates, " ARM_ERR: dates: array of numeric expected",C_result);
    }

	XL_readStrCellWD(XL_pf,C_pf,"DEFAULT"," ARM_ERR: portfolio id: object expected",C_result);

	if(C_pf == "DEFAULT")
	{
		pfId = ARM_NULL_OBJECT;
	}
	else
	{
		pfId = LocalGetNumObjectId (C_pf);
	}

	XL_readStrCellWD(XL_rhoswopt,C_RhoSwopt,"DEFAULT"," ARM_ERR: rhoswopt id: object expected",C_result);

	if(C_RhoSwopt == "DEFAULT")
	{
		RhoSwoptId = ARM_NULL_OBJECT;
	}
	else
	{
		RhoSwoptId = LocalGetNumObjectId (C_RhoSwopt);
	}

	XL_readStrCellWD(XL_nuswopt,C_NuSwopt,"DEFAULT"," ARM_ERR: nuswopt id: object expected",C_result);

	if(C_NuSwopt == "DEFAULT")
	{
		NuSwoptId = ARM_NULL_OBJECT;
	}
	else
	{
		NuSwoptId = LocalGetNumObjectId (C_NuSwopt);
	}

	XL_readStrCellWD(XL_betaswopt,C_BetaSwopt,"DEFAULT"," ARM_ERR: betaswopt id: object expected",C_result);

	if(C_BetaSwopt == "DEFAULT")
	{
		BetaSwoptId = ARM_NULL_OBJECT;
	}
	else
	{
		BetaSwoptId = LocalGetNumObjectId (C_BetaSwopt);
	}

	long retCode;
	long objId;
	
	CCString curClass = LOCAL_CALIB_HWSV_CLASS;
	CCString stringId;

	retCode = ARMLOCAL_CALIBRATIONHWSV(C_asof,
										LocalGetNumObjectId (C_Zc),
										LocalGetNumObjectId (C_Vol),
										LocalGetNumObjectId (C_sec),
										C_amin,
										C_amax,
										C_volmin,
										C_volmax,
										C_dates,
										pfId,
										RhoSwoptId,
										NuSwoptId,
										BetaSwoptId,
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
	   PXL_ARM_ERR();
	}

//	ARM_END();
	
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_PXL_ARM_CALIBRATIONHWSV" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}



__declspec(dllexport) LPXLOPER WINAPI Local_ARM_CALIBRATE(LPXLOPER XL_calibrator,
														  LPXLOPER XL_calibType)
{
	ADD_LOG("Local_ARM_CALIBRATE");
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();

	// C variable
	CCString C_calib;

	double C_calibType;
	double C_calibType_default = 2.;

	// error
	static int error;
	static char* reason = "";

	XL_readStrCell(XL_calibrator,C_calib," ARM_ERR: calibrator: object expected",C_result);
	XL_readNumCellWD(XL_calibType,C_calibType,C_calibType_default," ARM_ERR: calibration type: numeric expected",C_result);

	long retCode;
	long objId;
	CCString prevClass;
	
	CCString curClass = LOCAL_HWSIGVAR_ANALYTIC_CLASS;
	CCString stringId = GetLastCurCellEnvValue ();

	if (!stringId)
	{
		retCode = ARMLOCAL_CALIBRATE (LocalGetNumObjectId (C_calib),
									  (long) C_calibType,
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
			
		if(curClass == prevClass)
		{
			retCode = ARMLOCAL_CALIBRATE (LocalGetNumObjectId (C_calib),
										  (long) C_calibType,
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

			retCode = ARMLOCAL_CALIBRATE (LocalGetNumObjectId (C_calib),
										  (long) C_calibType,
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

//	ARM_END();
	
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_ARM_CALIBRATE" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


__declspec(dllexport) LPXLOPER WINAPI Local_PXL_ARM_CALIBRATE(LPXLOPER XL_calibrator,
															  LPXLOPER XL_calibType)
{
	ADD_LOG("Local_PXL_ARM_CALIBRATE");
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();

	// C variable
	CCString C_calib;

	double C_calibType;
	double C_calibType_default = 2.;

	// error
	static int error;
	static char* reason = "";

	XL_readStrCell(XL_calibrator,C_calib," ARM_ERR: calibrator: object expected",C_result);
	XL_readNumCellWD(XL_calibType,C_calibType,C_calibType_default," ARM_ERR: calibration type: numeric expected",C_result);

	long retCode;
	long objId;
	
	CCString curClass = LOCAL_HWSIGVAR_ANALYTIC_CLASS;
	CCString stringId;

	retCode = ARMLOCAL_CALIBRATE (LocalGetNumObjectId (C_calib),
								  (long) C_calibType,
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
	   PXL_ARM_ERR();
	}

//	ARM_END();
	
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_PXL_ARM_CALIBRATE" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


__declspec(dllexport) LPXLOPER WINAPI Local_ARM_HWSIGVARFromAna(LPXLOPER XL_modId,
																LPXLOPER XL_endDate,
																LPXLOPER XL_nbSteps)
{
	ADD_LOG("Local_ARM_HWSIGVARFromAna");
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();

	// C variable
	CCString C_modId;
	double C_endDate;	
	double C_nbSteps;

	// error
	static int error;
	static char* reason = "";

    XL_readStrCell(XL_modId, C_modId," ARM_ERR: Mod id: string expected",C_result);
    XL_readNumCell(XL_endDate,C_endDate," ARM_ERR: end date: date expected",C_result);
    XL_readNumCell(XL_nbSteps, C_nbSteps," ARM_ERR: number of steps: numeric expected",C_result);

	long j = 0;

	long retCode;
	long objId;
	CCString prevClass;
	CCString curClass = LOCAL_HWSIGVAR_TREE_CLASS;
	CCString stringId = GetLastCurCellEnvValue();

	if (!stringId)
	{
        retCode = ARMLOCAL_HWSIGVARFROMANA(LocalGetNumObjectId (C_modId), C_endDate,(long) C_nbSteps, C_result);

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
           retCode = ARMLOCAL_HWSIGVARFROMANA(	LocalGetNumObjectId (C_modId), 
												C_endDate,
												(long) C_nbSteps, C_result, objId);

			if ( retCode == ARM_OK )
			{
			   LocalSetCurCellEnvValue (curClass, objId); 

			   stringId = LocalMakeObjectId (objId, curClass);
			}
		}
		else
		{
			FreeCurCellContent ();

            retCode = ARMLOCAL_HWSIGVARFROMANA(LocalGetNumObjectId (C_modId), C_endDate,
										       (long) C_nbSteps, C_result);

	
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_ARM_HWSIGVARFromAna" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}



__declspec(dllexport) LPXLOPER WINAPI Local_PXL_ARM_HWSIGVARFromAna(LPXLOPER XL_modId,
																	LPXLOPER XL_endDate,
																	LPXLOPER XL_nbSteps)
{
	ADD_LOG("Local_PXL_ARM_HWSIGVARFromAna");
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	    ARM_NOCALCIFWIZ();

	    // C variable
	    CCString C_modId;
	    double C_endDate;	
	    double C_nbSteps;

	    // error
	    static int error;
	    static char* reason = "";

        XL_readStrCell(XL_modId, C_modId," ARM_ERR: Mod id: string expected",C_result);
        XL_readNumCell(XL_endDate,C_endDate," ARM_ERR: end date: date expected",C_result);
        XL_readNumCell(XL_nbSteps, C_nbSteps," ARM_ERR: number of steps: numeric expected",C_result);


	    long retCode;
	    long objId;
	    
	    CCString curClass = LOCAL_HWSIGVAR_TREE_CLASS;
	    CCString stringId;

        retCode = ARMLOCAL_HWSIGVARFROMANA(LocalGetNumObjectId (C_modId),C_endDate,
							      	      (long) C_nbSteps, C_result);

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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_PXL_ARM_HWSIGVARFromAna" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}



__declspec(dllexport) LPXLOPER WINAPI Local_ARM_BSSMILEDCALIBRATE(LPXLOPER XL_model,
																  LPXLOPER XL_pf,
																  LPXLOPER XL_calVolOrNot,
																  LPXLOPER XL_calRhoOrNot,
																  LPXLOPER XL_calNuOrNot,
																  LPXLOPER XL_volTenor,
																  LPXLOPER XL_timeStep,
																  LPXLOPER XL_minSig,
																  LPXLOPER XL_maxSig,
																  LPXLOPER XL_minRho,
																  LPXLOPER XL_maxRho,
																  LPXLOPER XL_minNu,
																  LPXLOPER XL_maxNu,
																  LPXLOPER XL_interpMethod,
																  LPXLOPER XL_tol,
																  LPXLOPER XL_maxIter,
																  LPXLOPER XL_gradCalc,
																  LPXLOPER XL_lambda,
																  LPXLOPER XL_globOrBootstrap,
																  LPXLOPER XL_SmoothMatuBetaParam)
{
	ADD_LOG("Local_ARM_BSSMILEDCALIBRATE");
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
	    CCString C_pf;

	    CCString C_calVolOrNot;
	    long calVolOrNotId;

	    CCString C_calRhoOrNot;
	    long calRhoOrNotId;

	    CCString C_calNuOrNot;
	    long calNuOrNotId;

        CCString C_calBetaOrNot;
	    long calBetaOrNotId;

	    double C_volTenor;
	    double C_volTenor_default = -1.;

	    double C_timeStep;
	    double C_timeStep_default = -1.;

	    double C_minSig;
	    double C_minSig_default = LOW_INFINITE_BOUND;

	    double C_maxSig;
	    double C_maxSig_default = UPPER_INFINITE_BOUND;

	    double C_minRho;
	    double C_minRho_default = LOW_INFINITE_BOUND;

	    double C_maxRho;
	    double C_maxRho_default = UPPER_INFINITE_BOUND;

	    double C_minNu;
	    double C_minNu_default = LOW_INFINITE_BOUND;

	    double C_maxNu;
	    double C_maxNu_default = UPPER_INFINITE_BOUND;

        double C_minBeta;
	    double C_minBeta_default = LOW_INFINITE_BOUND;

	    double C_maxBeta;
	    double C_maxBeta_default = UPPER_INFINITE_BOUND;


	    CCString C_interpMethod;
	    long interpMethodId;

	    double C_tol;
	    double C_tol_default = -1.;

	    double C_maxIter;
	    double C_maxIter_default = ARM_DEF_MAX_ITER;

	    CCString C_gradCalc;
	    long gradCalcId;

	    double C_lambda;
	    double C_lambda_default = 0.;

	    CCString C_globOrBootstrap;
	    long globOrBootstrapId;

	    double C_beginSmoothMatu;
	    double C_beginSmoothMatu_default = -1.;


	    // error
	    static int error;
	    static char* reason = "";

	    XL_readStrCell(XL_model,C_model," ARM_ERR: model: object expected",C_result);
	    XL_readStrCell(XL_pf,C_pf," ARM_ERR: pf: object expected",C_result);

	    XL_readStrCellWD(XL_calVolOrNot, C_calVolOrNot,"Y"," ARM_ERR: calVolOrNot: string expected",C_result);
	    XL_readStrCellWD(XL_calRhoOrNot, C_calRhoOrNot,"Y"," ARM_ERR: calRhoOrNot: string expected",C_result);
	    XL_readStrCellWD(XL_calNuOrNot,  C_calNuOrNot, "Y"," ARM_ERR: calNuOrNot: string expected",C_result);

	    XL_readNumCellWD(XL_volTenor,C_volTenor,C_volTenor_default," ARM_ERR: vol Tenor: numeric expected",C_result);
	    XL_readNumCellWD(XL_timeStep,C_timeStep,C_timeStep_default," ARM_ERR: time step: numeric expected",C_result);

	    XL_readNumCellWD(XL_minSig,C_minSig,C_minSig_default," ARM_ERR: minSig: numeric expected",C_result);
	    XL_readNumCellWD(XL_maxSig,C_maxSig,C_maxSig_default," ARM_ERR: maxSig: numeric expected",C_result);
	    XL_readNumCellWD(XL_minRho,C_minRho,C_minRho_default," ARM_ERR: minRho: numeric expected",C_result);
	    XL_readNumCellWD(XL_maxRho,C_maxRho,C_maxRho_default," ARM_ERR: maxRho: numeric expected",C_result);
	    XL_readNumCellWD(XL_minNu,C_minNu,C_minNu_default," ARM_ERR: minNu: numeric expected",C_result);
	    XL_readNumCellWD(XL_maxNu,C_maxNu,C_maxNu_default," ARM_ERR: maxNu: numeric expected",C_result);

	    XL_readStrCellWD(XL_interpMethod,C_interpMethod,"LINEAR"," ARM_ERR: interpolation method: string expected",C_result);
	    XL_readNumCellWD(XL_tol,C_tol,C_tol_default," ARM_ERR: tolerance: numeric expected",C_result);
	    XL_readNumCellWD(XL_maxIter,C_maxIter,C_maxIter_default," ARM_ERR: maximum of Iterations: numeric expected",C_result);
	    XL_readStrCellWD(XL_gradCalc,C_gradCalc,"Y"," ARM_ERR: calcGrad: string expected",C_result);

	    XL_readNumCellWD(XL_lambda,C_lambda,C_lambda_default," ARM_ERR: lambda : numeric expected",C_result);
	    XL_readStrCellWD(XL_globOrBootstrap,C_globOrBootstrap,"B"," ARM_ERR: calcGrad: string expected",C_result);

        // Beta

        VECTOR<double> C_SmoothMatuBetaParam;
        VECTOR<double> C_SmoothMatuBetaParamDef(4);
    
        C_SmoothMatuBetaParamDef[0] = C_beginSmoothMatu_default;
        C_SmoothMatuBetaParamDef[1] = 0;
        C_SmoothMatuBetaParamDef[2] = C_minBeta_default;
        C_SmoothMatuBetaParamDef[3] = C_maxBeta_default;

        XL_readNumVectorWD(XL_SmoothMatuBetaParam, C_SmoothMatuBetaParam, C_SmoothMatuBetaParamDef, 
           " ARM_ERR: Smooth matu. and beta calib. params.: numeric array expected: [(smooth matu(or -1)], [FitOrNot(0|1)], [Bmin], [Bmax]",
          C_result);

        switch(C_SmoothMatuBetaParam.size())
        {
            // case 0 : impossible

            case 1: 
            {
                 C_SmoothMatuBetaParam.push_back(0); // Dont fit beta
            };

            case 2 :
            {
                 C_SmoothMatuBetaParam.push_back(C_minBeta_default);
            };

            case 3:
            {
                 C_SmoothMatuBetaParam.push_back(C_maxBeta_default);
            };
            break;

            // case 4 : all parameters are updated
        }

        C_beginSmoothMatu = C_SmoothMatuBetaParam[0];
        calBetaOrNotId    = C_SmoothMatuBetaParam[1];
        C_minBeta         = C_SmoothMatuBetaParam[2];
        C_maxBeta         = C_SmoothMatuBetaParam[3]; 

	    if ( (calVolOrNotId = ARM_ConvYesOrNo (C_calVolOrNot, C_result)) == ARM_DEFAULT_ERR)
	    {
	       ARM_ARG_ERR();

	       return (LPXLOPER)&XL_result;
	    }

	    if ( (calRhoOrNotId = ARM_ConvYesOrNo (C_calRhoOrNot, C_result)) == ARM_DEFAULT_ERR)
	    {
	       ARM_ARG_ERR();

	       return (LPXLOPER)&XL_result;
	    }

	    if ( (calNuOrNotId = ARM_ConvYesOrNo (C_calNuOrNot, C_result)) == ARM_DEFAULT_ERR)
	    {
	       ARM_ARG_ERR();

	       return (LPXLOPER)&XL_result;
	    }

	    if ( (gradCalcId = ARM_ConvYesOrNo (C_gradCalc, C_result)) == ARM_DEFAULT_ERR)
	    {
	       ARM_ARG_ERR();

	       return (LPXLOPER)&XL_result;
	    }

	    if ( (interpMethodId = ARM_ConvInterpMethod (C_interpMethod, C_result)) == ARM_DEFAULT_ERR )
	    {
		    ARM_ARG_ERR();
		    return (LPXLOPER)&XL_result;
	    }

	    if ( (globOrBootstrapId = ARM_ConvBootstrap (C_globOrBootstrap, C_result)) == ARM_DEFAULT_ERR)
	    {
	       ARM_ARG_ERR();

	       return (LPXLOPER)&XL_result;
	    }

        C_beginSmoothMatu = C_SmoothMatuBetaParam[0];
        calBetaOrNotId    = C_SmoothMatuBetaParam[1];
        C_minBeta         = C_SmoothMatuBetaParam[2];
        C_maxBeta         = C_SmoothMatuBetaParam[3];

	    long retCode;
	    long objId;
	    CCString prevClass;
	    
	    CCString curClass = LOCAL_BSMODEL_CLASS;
	    CCString stringId = GetLastCurCellEnvValue ();

	    if (!stringId)
	    {
		    retCode = ARMLOCAL_BSSMILEDCALIBRATE(LocalGetNumObjectId(C_model),
											     LocalGetNumObjectId(C_pf),
											     calVolOrNotId,
											     calRhoOrNotId,
											     calNuOrNotId,
                                                 calBetaOrNotId,
											     C_volTenor,
											     C_timeStep,
											     C_minSig,
											     C_maxSig,
											     C_minRho,
											     C_maxRho,
											     C_minNu,
											     C_maxNu,
                                                 C_minBeta,
											     C_maxBeta,
											     interpMethodId,
											     C_tol,
											     (long)C_maxIter,
											     gradCalcId,
											     C_lambda,
											     globOrBootstrapId,
											     C_beginSmoothMatu,
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

		    if(curClass == prevClass)
		    {
			    retCode = ARMLOCAL_BSSMILEDCALIBRATE(LocalGetNumObjectId(C_model),
												     LocalGetNumObjectId(C_pf),
												     calVolOrNotId,
												     calRhoOrNotId,
												     calNuOrNotId,
                                                     calBetaOrNotId,
												     C_volTenor,
												     C_timeStep,
												     C_minSig,
												     C_maxSig,
												     C_minRho,
												     C_maxRho,
												     C_minNu,
												     C_maxNu,
                                                     C_minBeta,
											         C_maxBeta,
												     interpMethodId,
												     C_tol,
												     (long)C_maxIter,
												     gradCalcId,
												     C_lambda,
												     globOrBootstrapId,
												     C_beginSmoothMatu,
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

			    retCode = ARMLOCAL_BSSMILEDCALIBRATE(LocalGetNumObjectId(C_model),
												     LocalGetNumObjectId(C_pf),
												     calVolOrNotId,
												     calRhoOrNotId,
												     calNuOrNotId,
                                                     calBetaOrNotId,
												     C_volTenor,
												     C_timeStep,
												     C_minSig,
												     C_maxSig,
												     C_minRho,
												     C_maxRho,
												     C_minNu,
												     C_maxNu,
                                                     C_minBeta,
											         C_maxBeta,
												     interpMethodId,
												     C_tol,
												     (long)C_maxIter,
												     gradCalcId,
												     C_lambda,
												     globOrBootstrapId,
												     C_beginSmoothMatu,
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
    ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_ARM_BSSMILEDCALIBRATE");

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}



__declspec(dllexport) LPXLOPER WINAPI Local_PXL_ARM_BSSMILEDCALIBRATE(LPXLOPER XL_model,
																	  LPXLOPER XL_pf,
																	  LPXLOPER XL_calVolOrNot,
																	  LPXLOPER XL_calRhoOrNot,
																	  LPXLOPER XL_calNuOrNot,
																	  LPXLOPER XL_volTenor,
																	  LPXLOPER XL_timeStep,
																	  LPXLOPER XL_minSig,
																	  LPXLOPER XL_maxSig,
																	  LPXLOPER XL_minRho,
																	  LPXLOPER XL_maxRho,
																	  LPXLOPER XL_minNu,
																	  LPXLOPER XL_maxNu,
																	  LPXLOPER XL_interpMethod,
																	  LPXLOPER XL_tol,
																	  LPXLOPER XL_maxIter,
																	  LPXLOPER XL_gradCalc,
																	  LPXLOPER XL_lambda,
																	  LPXLOPER XL_globOrBootstrap,
																	  LPXLOPER XL_SmoothMatuBetaParam)
{
	ADD_LOG("Local_PXL_ARM_BSSMILEDCALIBRATE");
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
	    CCString C_pf;

	    CCString C_calVolOrNot;
	    long calVolOrNotId;

	    CCString C_calRhoOrNot;
	    long calRhoOrNotId;

	    CCString C_calNuOrNot;
	    long calNuOrNotId;

        CCString C_calBetaOrNot;
	    long calBetaOrNotId;

	    double C_volTenor;
	    double C_volTenor_default = -1.;

	    double C_timeStep;
	    double C_timeStep_default = -1.;

	    double C_minSig;
	    double C_minSig_default = LOW_INFINITE_BOUND;

	    double C_maxSig;
	    double C_maxSig_default = UPPER_INFINITE_BOUND;

	    double C_minRho;
	    double C_minRho_default = LOW_INFINITE_BOUND;

	    double C_maxRho;
	    double C_maxRho_default = UPPER_INFINITE_BOUND;

	    double C_minNu;
	    double C_minNu_default = LOW_INFINITE_BOUND;

	    double C_maxNu;
	    double C_maxNu_default = UPPER_INFINITE_BOUND;

        double C_minBeta;
	    double C_minBeta_default = LOW_INFINITE_BOUND;

        double C_maxBeta;
	    double C_maxBeta_default = UPPER_INFINITE_BOUND;

	    CCString C_interpMethod;
	    long interpMethodId;

	    double C_tol;
	    double C_tol_default = -1.;

	    double C_maxIter;
	    double C_maxIter_default = ARM_DEF_MAX_ITER;

	    CCString C_gradCalc;
	    long gradCalcId;

	    double C_lambda;
	    double C_lambda_default = 0.;

	    CCString C_globOrBootstrap;
	    long globOrBootstrapId;

	    double C_beginSmoothMatu;
	    double C_beginSmoothMatu_default = -1.;

	    // error
	    static int error;
	    static char* reason = "";

	    XL_readStrCell(XL_model,C_model," ARM_ERR: model: object expected",C_result);
	    XL_readStrCell(XL_pf,C_pf," ARM_ERR: pf: object expected",C_result);

	    XL_readStrCellWD(XL_calVolOrNot,C_calVolOrNot,"Y"," ARM_ERR: calVolOrNot: string expected",C_result);
	    XL_readStrCellWD(XL_calRhoOrNot,C_calRhoOrNot,"Y"," ARM_ERR: calRhoOrNot: string expected",C_result);
	    XL_readStrCellWD(XL_calNuOrNot,C_calNuOrNot,"Y"," ARM_ERR: calNuOrNot: string expected",C_result);

	    XL_readNumCellWD(XL_volTenor,C_volTenor,C_volTenor_default," ARM_ERR: vol Tenor: numeric expected",C_result);
	    XL_readNumCellWD(XL_timeStep,C_timeStep,C_timeStep_default," ARM_ERR: time step: numeric expected",C_result);

	    XL_readNumCellWD(XL_minSig,C_minSig,C_minSig_default," ARM_ERR: minSig: numeric expected",C_result);
	    XL_readNumCellWD(XL_maxSig,C_maxSig,C_maxSig_default," ARM_ERR: maxSig: numeric expected",C_result);
	    XL_readNumCellWD(XL_minRho,C_minRho,C_minRho_default," ARM_ERR: minRho: numeric expected",C_result);
	    XL_readNumCellWD(XL_maxRho,C_maxRho,C_maxRho_default," ARM_ERR: maxRho: numeric expected",C_result);
	    XL_readNumCellWD(XL_minNu,C_minNu,C_minNu_default," ARM_ERR: minNu: numeric expected",C_result);
	    XL_readNumCellWD(XL_maxNu,C_maxNu,C_maxNu_default," ARM_ERR: maxNu: numeric expected",C_result);

	    XL_readStrCellWD(XL_interpMethod,C_interpMethod,"LINEAR"," ARM_ERR: interpolation method: string expected",C_result);
	    XL_readNumCellWD(XL_tol,C_tol,C_tol_default," ARM_ERR: tolerance: numeric expected",C_result);
	    XL_readNumCellWD(XL_maxIter,C_maxIter,C_maxIter_default," ARM_ERR: maximum of Iterations: numeric expected",C_result);
	    XL_readStrCellWD(XL_gradCalc,C_gradCalc,"Y"," ARM_ERR: calcGrad: string expected",C_result);

	    XL_readNumCellWD(XL_lambda,C_lambda,C_lambda_default," ARM_ERR: lambda : numeric expected",C_result);
	    XL_readStrCellWD(XL_globOrBootstrap,C_globOrBootstrap,"B"," ARM_ERR: calcGrad: string expected",C_result);
	
      // Beta

        VECTOR<double> C_SmoothMatuBetaParam;
        VECTOR<double> C_SmoothMatuBetaParamDef(4);
    
        C_SmoothMatuBetaParamDef[0] = C_beginSmoothMatu_default;
        C_SmoothMatuBetaParamDef[1] = 0;
        C_SmoothMatuBetaParamDef[2] = C_minBeta_default;
        C_SmoothMatuBetaParamDef[3] = C_maxBeta_default;

        XL_readNumVectorWD(XL_SmoothMatuBetaParam, C_SmoothMatuBetaParam, C_SmoothMatuBetaParamDef, 
           " ARM_ERR: Smooth matu. and beta calib. params.: numeric array expected: [smooth matu(or -1)], [FitOrNot(0|1)], [Bmin], [Bmax]",
          C_result);

        switch(C_SmoothMatuBetaParam.size())
        {
            // case 0 : impossible

            case 1: 
            {
                 C_SmoothMatuBetaParam.push_back(0); // Dont fit beta
            };

            case 2 :
            {
                 C_SmoothMatuBetaParam.push_back(C_minBeta_default);
            };

            case 3:
            {
                 C_SmoothMatuBetaParam.push_back(C_maxBeta_default);
            };
            break;

            // case 4 : all parameters are updated
        }

        C_beginSmoothMatu = C_SmoothMatuBetaParam[0];
        calBetaOrNotId    = C_SmoothMatuBetaParam[1];
        C_minBeta         = C_SmoothMatuBetaParam[2];
        C_maxBeta         = C_SmoothMatuBetaParam[3]; 

	    if ( (calVolOrNotId = ARM_ConvYesOrNo (C_calVolOrNot, C_result)) == ARM_DEFAULT_ERR)
	    {
	       ARM_ARG_ERR();

	       return (LPXLOPER)&XL_result;
	    }

	    if ( (calRhoOrNotId = ARM_ConvYesOrNo (C_calRhoOrNot, C_result)) == ARM_DEFAULT_ERR)
	    {
	       ARM_ARG_ERR();

	       return (LPXLOPER)&XL_result;
	    }

	    if ( (calNuOrNotId = ARM_ConvYesOrNo (C_calNuOrNot, C_result)) == ARM_DEFAULT_ERR)
	    {
	       ARM_ARG_ERR();

	       return (LPXLOPER)&XL_result;
	    }

	    if ( (gradCalcId = ARM_ConvYesOrNo (C_gradCalc, C_result)) == ARM_DEFAULT_ERR)
	    {
	       ARM_ARG_ERR();

	       return (LPXLOPER)&XL_result;
	    }

	    if((interpMethodId = ARM_ConvInterpMethod (C_interpMethod, C_result)) == ARM_DEFAULT_ERR)
	    {
		    ARM_ARG_ERR();
		    return (LPXLOPER)&XL_result;
	    }

	    if ( (globOrBootstrapId = ARM_ConvBootstrap (C_globOrBootstrap, C_result)) == ARM_DEFAULT_ERR)
	    {
	       ARM_ARG_ERR();

	       return (LPXLOPER)&XL_result;
	    }

	    long retCode;
	    long objId;
	    
	    CCString curClass = LOCAL_BSMODEL_CLASS;
	    CCString stringId;

	    retCode = ARMLOCAL_BSSMILEDCALIBRATE(LocalGetNumObjectId(C_model),
										     LocalGetNumObjectId(C_pf),
										     calVolOrNotId,
										     calRhoOrNotId,
										     calNuOrNotId,
                                             calBetaOrNotId,
										     C_volTenor,
										     C_timeStep,
										     C_minSig,
										     C_maxSig,
										     C_minRho,
										     C_maxRho,
										     C_minNu,
										     C_maxNu,
                                             C_minBeta,
											 C_maxBeta,
										     interpMethodId,
										     C_tol,
										     (long)C_maxIter,
										     gradCalcId,
										     C_lambda,
										     globOrBootstrapId,
										     C_beginSmoothMatu,
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
	       PXL_ARM_ERR();
	    }

    //	ARM_END();
	
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_PXL_ARM_BSSMILEDCALIBRATE" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}




__declspec(dllexport) LPXLOPER WINAPI Local_ARM_GETCALIBRATED_SIGRHONU(LPXLOPER XL_model,
																	   LPXLOPER XL_param)
{
	ADD_LOG("Local_ARM_GETCALIBRATED_SIGRHONU");
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

	CCString C_param;
	long paramId;

	// error
	static int error;
	static char* reason = "";

	XL_readStrCell(XL_model,C_model," ARM_ERR: model: object expected",C_result);
	XL_readStrCell(XL_param,C_param," ARM_ERR: param: string expected",C_result);

	if((paramId = ARM_ConvSigRhoNu (C_param, C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	long retCode;
	long objId;
	CCString prevClass;
	
	CCString curClass = LOCAL_VOL_CURVE_LIN_CLASS;
	CCString stringId = GetLastCurCellEnvValue ();

	if (!stringId)
	{
		retCode = ARMLOCAL_GETCALIBRATED_SIGRHONU(LocalGetNumObjectId(C_model),
												  paramId,
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

		if(curClass == prevClass)
		{
			retCode = ARMLOCAL_GETCALIBRATED_SIGRHONU(LocalGetNumObjectId(C_model),
													  paramId,
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

			retCode = ARMLOCAL_GETCALIBRATED_SIGRHONU(LocalGetNumObjectId(C_model),
													  paramId,
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

//	ARM_END();
	
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_ARM_GETCALIBRATED_SIGRHONU" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}



__declspec(dllexport) LPXLOPER WINAPI Local_PXL_ARM_GETCALIBRATED_SIGRHONU(LPXLOPER XL_model,
																		   LPXLOPER XL_param)
{
	ADD_LOG("Local_PXL_ARM_GETCALIBRATED_SIGRHONU");
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

	CCString C_param;
	long paramId;

	// error
	static int error;
	static char* reason = "";

	XL_readStrCell(XL_model,C_model," ARM_ERR: model: object expected",C_result);
	XL_readStrCell(XL_param,C_param," ARM_ERR: param: string expected",C_result);

	if((paramId = ARM_ConvSigRhoNu (C_param, C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	long retCode;
	long objId;
	
	CCString curClass = LOCAL_VOL_CURVE_LIN_CLASS;
	CCString stringId;

	retCode = ARMLOCAL_GETCALIBRATED_SIGRHONU(LocalGetNumObjectId(C_model),
											  paramId,
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
	   PXL_ARM_ERR();
	}

//	ARM_END();
	
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_PXL_ARM_GETCALIBRATED_SIGRHONU" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}

_declspec(dllexport) LPXLOPER WINAPI Local_ARM_CalibrationFRMModel(LPXLOPER XL_zcId,
																	LPXLOPER XL_liborType,
                                                                    LPXLOPER XL_nbFactors,
                                                                    LPXLOPER XL_initA,
																	LPXLOPER XL_initK,
                                                                    LPXLOPER XL_meanRevA,
																	LPXLOPER XL_meanRevK,                                                                    
                                                                    LPXLOPER XL_mdec,
                                                                    LPXLOPER XL_schedPrice,
                                                                    LPXLOPER XL_pf1Id,
                                                                    LPXLOPER XL_pf2Id,
                                                                    LPXLOPER XL_ACalibrationSchedule,
																	LPXLOPER XL_KCalibrationSchedule,                                                                    
                                                                    LPXLOPER XL_powers,
																	LPXLOPER XL_smoothParams,
																	LPXLOPER XL_bounds,
																	LPXLOPER XL_optimizerParams,
                                                                    LPXLOPER XL_calib)
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
	CCString C_zcId;
	CCString C_pf1Id;
    CCString C_pf2Id;
    long pf2Id;

	CCString C_liborType;
	long liborTypeId;

	CCString C_calib;
    long calibId;

	VECTOR<double> C_Powers;
	VECTOR<double> C_smoothParams;

	double C_nbFactors;

	VECTOR<double> C_ACalibrationSchedule;
	VECTOR<double> C_KCalibrationSchedule;

    VECTOR<double> C_schedPrice;

	CCString C_initA;
	CCString C_initK;

	VECTOR<double> C_meanRevA;
	VECTOR<double> C_meanRevK;
    
    VECTOR<double> C_mdec;
	
    VECTOR<double> C_bounds;

	VECTOR<double> C_optimizerParams;

	// error
	static int error;
	static char* reason = "";

    XL_readStrCell(XL_zcId, C_zcId," ARM_ERR: Zc id: object expected",C_result);
    XL_readNumCell(XL_nbFactors, C_nbFactors," ARM_ERR: NbFactors id: object expected",C_result);
    XL_readStrCell(XL_liborType, C_liborType," ARM_ERR: libor type: string expected",C_result);
    XL_readStrCell(XL_initA, C_initA," ARM_ERR: init A: object expected",C_result);
    XL_readStrCell(XL_initK, C_initK," ARM_ERR: init K: object expected",C_result);
    XL_readNumVector(XL_meanRevA, C_meanRevA," ARM_ERR: Mean Reversion A: vector expected",C_result);
    XL_readNumVector(XL_meanRevK, C_meanRevK," ARM_ERR: Mean Reversion K: vector expected",C_result);
    XL_readNumVector(XL_schedPrice, C_schedPrice," ARM_ERR: Sched Price: vector expected",C_result);
    XL_readNumVector(XL_mdec, C_mdec," ARM_ERR: mdec: matrix expected",C_result);
    XL_readStrCell(XL_pf1Id, C_pf1Id," ARM_ERR: Pf1 id: object expected",C_result);

    if((liborTypeId = ARM_ConvIrIndName (C_liborType, C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

    XL_readStrCellWD(XL_calib, C_calib,"NO"," ARM_ERR: Calib: string expected",C_result);

    if((calibId = ARM_ConvYesOrNo (C_calib, C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

    if (calibId == K_YES)
    {
        XL_readStrCellWD(XL_pf2Id, C_pf2Id,"NULL"," ARM_ERR: Pf2 id: object expected",C_result);
        XL_readNumVector(XL_ACalibrationSchedule, C_ACalibrationSchedule," ARM_ERR: Calibration Schedule A: matrix expected",C_result);
        XL_readNumVector(XL_KCalibrationSchedule, C_KCalibrationSchedule," ARM_ERR: Calibration Schedule K: matrix expected",C_result);
        XL_readNumVector(XL_powers,C_Powers," ARM_ERR: market and smooth powers: vector expected",C_result);
        XL_readNumVector(XL_smoothParams,C_smoothParams," ARM_ERR: smooth parameters: vector expected",C_result);    
        XL_readNumVector(XL_bounds, C_bounds," ARM_ERR: Bounds: matrix expected",C_result);
        XL_readNumVector(XL_optimizerParams, C_optimizerParams," ARM_ERR: Optimizer Parameters: vector expected",C_result);
    }
    else
    {
        C_pf2Id = "NULL";
       /* C_ACalibrationSchedule = NULL;
        C_KCalibrationSchedule= NULL;
        C_Powers= NULL;
        C_smoothParams = = NULL;
        C_bounds = NULL;
        C_optimizerParams = NULL;  */     

    }   

    if (C_pf2Id == "NULL")
        pf2Id = ARM_NULL_OBJECT;
    else
        pf2Id = LocalGetNumObjectId (C_pf2Id);

	long retCode;
	long objId;
	CCString prevClass;

	CCString curClass = LOCAL_FRMMODEL_CLASS;
	CCString stringId = GetLastCurCellEnvValue();
	
	if (!stringId)
	{
		retCode = ARMLOCAL_CalibrationFRMModel(LocalGetNumObjectId (C_zcId),
											   LocalGetNumObjectId (C_pf1Id),
                                               pf2Id,
											   C_nbFactors,
											   liborTypeId,
											   C_Powers,
											   C_smoothParams,
											   C_ACalibrationSchedule,
											   C_KCalibrationSchedule,
											   LocalGetNumObjectId(C_initA),
											   LocalGetNumObjectId(C_initK),
											   C_meanRevA,
											   C_meanRevK,
                                               C_schedPrice,
                                               C_mdec,
											   C_bounds,
											   C_optimizerParams,
                                               calibId,
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
			retCode = ARMLOCAL_CalibrationFRMModel(LocalGetNumObjectId (C_zcId),
												   LocalGetNumObjectId (C_pf1Id),
                                                   pf2Id,
												   C_nbFactors,
												   liborTypeId,
												   C_Powers,
												   C_smoothParams,
												   C_ACalibrationSchedule,
												   C_KCalibrationSchedule,
												   LocalGetNumObjectId(C_initA),
												   LocalGetNumObjectId(C_initK),
												   C_meanRevA,
												   C_meanRevK,
                                                   C_schedPrice,
                                                   C_mdec,
												   C_bounds,
												   C_optimizerParams,
                                                   calibId,
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

			retCode = ARMLOCAL_CalibrationFRMModel(LocalGetNumObjectId (C_zcId),
												   LocalGetNumObjectId (C_pf1Id),
                                                   pf2Id,
												   C_nbFactors,
												   liborTypeId,
												   C_Powers,
											       C_smoothParams,
												   C_ACalibrationSchedule,
												   C_KCalibrationSchedule,
												   LocalGetNumObjectId(C_initA),
												   LocalGetNumObjectId(C_initK),
												   C_meanRevA,
												   C_meanRevK,
                                                   C_schedPrice,
                                                   C_mdec,
												   C_bounds,
												   C_optimizerParams,
                                                   calibId,
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_ARM_CalibrationFRMModel" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}

_declspec(dllexport) LPXLOPER WINAPI Local_PXL_ARM_CalibrationFRMModel(LPXLOPER XL_zcId,
																	LPXLOPER XL_liborType,
                                                                    LPXLOPER XL_nbFactors,
                                                                    LPXLOPER XL_initA,
																	LPXLOPER XL_initK,
                                                                    LPXLOPER XL_meanRevA,
																	LPXLOPER XL_meanRevK,
                                                                    LPXLOPER XL_mdec,
                                                                    LPXLOPER XL_schedPrice,
                                                                    LPXLOPER XL_pf1Id,
                                                                    LPXLOPER XL_pf2Id,
                                                                    LPXLOPER XL_ACalibrationSchedule,
																	LPXLOPER XL_KCalibrationSchedule,                                                                    
                                                                    LPXLOPER XL_powers,
																	LPXLOPER XL_smoothParams,
																	LPXLOPER XL_bounds,
																	LPXLOPER XL_optimizerParams,
                                                                    LPXLOPER XL_calib)
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
	CCString C_zcId;
	CCString C_pf1Id;
    CCString C_pf2Id;
    long pf2Id;

	CCString C_liborType;
	long liborTypeId;

	CCString C_calib;
    long calibId;

	VECTOR<double> C_Powers;
	VECTOR<double> C_smoothParams;

	double C_nbFactors;

	VECTOR<double> C_ACalibrationSchedule;
	VECTOR<double> C_KCalibrationSchedule;

    VECTOR<double> C_schedPrice;

	CCString C_initA;
	CCString C_initK;

	VECTOR<double> C_meanRevA;
	VECTOR<double> C_meanRevK;
    
    VECTOR<double> C_mdec;
	
    VECTOR<double> C_bounds;

	VECTOR<double> C_optimizerParams;

	// error
	static int error;
	static char* reason = "";

    XL_readStrCell(XL_zcId, C_zcId," ARM_ERR: Zc id: object expected",C_result);
    XL_readNumCell(XL_nbFactors, C_nbFactors," ARM_ERR: NbFactors id: object expected",C_result);
    XL_readStrCell(XL_liborType, C_liborType," ARM_ERR: libor type: string expected",C_result);
    XL_readStrCell(XL_initA, C_initA," ARM_ERR: init A: object expected",C_result);
    XL_readStrCell(XL_initK, C_initK," ARM_ERR: init K: object expected",C_result);
    XL_readNumVector(XL_meanRevA, C_meanRevA," ARM_ERR: Mean Reversion A: vector expected",C_result);
    XL_readNumVector(XL_meanRevK, C_meanRevK," ARM_ERR: Mean Reversion K: vector expected",C_result);
    XL_readNumVector(XL_schedPrice, C_schedPrice," ARM_ERR: Sched Price: vector expected",C_result);
    XL_readNumVector(XL_mdec, C_mdec," ARM_ERR: init K: matrix expected",C_result);
    XL_readStrCell(XL_pf1Id, C_pf1Id," ARM_ERR: Pf1 id: object expected",C_result);

    if((liborTypeId = ARM_ConvIrIndName (C_liborType, C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

    XL_readStrCellWD(XL_calib, C_calib,"NO"," ARM_ERR: Calib: string expected",C_result);

    if((calibId = ARM_ConvYesOrNo (C_calib, C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

    if (calibId == K_YES)
    {
        XL_readStrCellWD(XL_pf2Id, C_pf2Id,"NULL"," ARM_ERR: Pf2 id: object expected",C_result);
        XL_readNumVector(XL_ACalibrationSchedule, C_ACalibrationSchedule," ARM_ERR: Calibration Schedule A: matrix expected",C_result);
        XL_readNumVector(XL_KCalibrationSchedule, C_KCalibrationSchedule," ARM_ERR: Calibration Schedule K: matrix expected",C_result);
        XL_readNumVector(XL_powers,C_Powers," ARM_ERR: market and smooth powers: vector expected",C_result);
        XL_readNumVector(XL_smoothParams,C_smoothParams," ARM_ERR: smooth parameters: vector expected",C_result);    
        XL_readNumVector(XL_bounds, C_bounds," ARM_ERR: Bounds: matrix expected",C_result);
        XL_readNumVector(XL_optimizerParams, C_optimizerParams," ARM_ERR: Optimizer Parameters: vector expected",C_result);
    }
    else
    {
        C_pf2Id = "NULL";    
    }    

    if (C_pf2Id == "NULL")
        pf2Id = ARM_NULL_OBJECT;
    else
        pf2Id = LocalGetNumObjectId (C_pf2Id);

	long retCode;
	long objId;
	CCString prevClass;

	CCString curClass = LOCAL_FRMMODEL_CLASS;
	CCString stringId = GetLastCurCellEnvValue();	
	
	retCode = ARMLOCAL_CalibrationFRMModel(LocalGetNumObjectId (C_zcId),
										   LocalGetNumObjectId (C_pf1Id),
                                           pf2Id,
										   C_nbFactors,
										   liborTypeId,
										   C_Powers,
										   C_smoothParams,
										   C_ACalibrationSchedule,
										   C_KCalibrationSchedule,
										   LocalGetNumObjectId(C_initA),
										   LocalGetNumObjectId(C_initK),
										   C_meanRevA,
										   C_meanRevK,
                                           C_schedPrice,
                                           C_mdec,
										   C_bounds,
										   C_optimizerParams,
                                           calibId,
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
	   PXL_ARM_ERR();
	}

//	ARM_END();
	
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_PXL_ARM_CalibrationFRMModel" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}



__declspec(dllexport) LPXLOPER WINAPI Local_ARM_BootstCalibFRMModel(LPXLOPER XL_zcId,                                                 
																	LPXLOPER XL_liborType,                                                                    
																	LPXLOPER XL_CalibParams,
																	LPXLOPER XL_initCurve,
																	LPXLOPER XL_mdec,										
																	LPXLOPER XL_meanRev,
                                                                    LPXLOPER XL_pf1Id,
                                                                    LPXLOPER XL_nbfactor,										
																	LPXLOPER XL_correlmatrix,
                                                                    LPXLOPER XL_pf2Id,
                                                                    LPXLOPER XL_VolType,
                                                                    LPXLOPER XL_Calib,
                                                                    LPXLOPER XL_pf3Id,
                                                                    LPXLOPER XL_marketmodelId,
                                                                    LPXLOPER XL_flagPreInit,
                                                                    LPXLOPER XL_precision,
                                                                    LPXLOPER XL_vegaLevel)
{
	ADD_LOG("Local_ARM_BootstCalibFRMModel");
	
	//	ARM_BEGIN();
	
	// return
	static XLOPER XL_result;
	ARM_result C_result;
	
	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();
		
		// C variable
		CCString C_zcId;
		CCString C_volId;
		CCString C_pf1Id;
		CCString C_pf2Id;
		CCString C_pf3Id;
		CCString C_marketmodelId;
		
		long pf2Id;
		long pf3Id;
		long marketmodelId;
		
		CCString C_liborType;
		
		CCString C_VolType;
		long volTypeId;
		
		CCString C_Calib;
		long calibId;
		
		long liborTypeId;
		
		VECTOR<double> C_CalibParams;
		
		VECTOR<double> C_initCurve;
		VECTOR<double> C_mdec;
		
		
		double C_meanRev;
		
		double C_nbfactor;
		CCString C_flagPreInit;
		
		VECTOR<double> C_correlmatrix;
		long C_nbrows;
		long C_nbcolumns;
		
		// error
		static int error;
		static char* reason = "";
		
		XL_readStrCell(XL_zcId, C_zcId," ARM_ERR: Zc id: object expected",C_result);
		XL_readStrCell(XL_liborType, C_liborType," ARM_ERR: libor type: string expected",C_result);
		XL_readNumVector(XL_CalibParams, C_CalibParams," ARM_ERR: Calib Params: vector expected",C_result);
		XL_readNumVector(XL_initCurve, C_initCurve," ARM_ERR: init Curve: vector expected",C_result);
		XL_readNumVector(XL_mdec, C_mdec," ARM_ERR: init shift or LogNorProbability Curve: vector expected",C_result);
		XL_readNumCell(XL_meanRev,C_meanRev," ARM_ERR: Mean Reversion: numeric expected",C_result);	   
		XL_readStrCell(XL_pf1Id, C_pf1Id," ARM_ERR: Pf1 id: object expected",C_result);
		XL_readNumCell(XL_nbfactor,C_nbfactor," ARM_ERR:Nb of factors: numeric expected",C_result);
		
		if(C_nbfactor > 1)
			XL_readNumVectorAndSize(XL_correlmatrix,C_nbrows,C_nbcolumns,C_correlmatrix," ARM_ERR: init Curve: matrix expected",C_result);
		
		XL_readStrCellWD(XL_pf2Id, C_pf2Id,"NULL"," ARM_ERR: Pf2 id: object expected",C_result);
		XL_readStrCellWD(XL_VolType,C_VolType,"ROW"," ARM_ERR: Vol Type: string expected",C_result);
		XL_readStrCellWD(XL_Calib,C_Calib,"YES"," ARM_ERR: Calib: string expected",C_result);
		XL_readStrCellWD(XL_pf3Id, C_pf3Id,"NULL"," ARM_ERR: PfToCalibrateShift id: object expected",C_result);
		XL_readStrCellWD(XL_marketmodelId, C_marketmodelId,"NULL"," ARM_ERR: Market Model id: object expected",C_result);
		XL_readStrCellWD(XL_flagPreInit,C_flagPreInit,"NO"," ARM_ERR: FlagPreInit: string expected",C_result);
		
		double C_precision;
		double C_precision_default = 1.0e-3;
		XL_readNumCellWD(XL_precision,C_precision,C_precision_default," ARM_ERR: precision: numeric expected",C_result);

        double C_vegaLevel;
		double C_vegaLevel_default = 1.0e-10;
		XL_readNumCellWD(XL_vegaLevel,C_vegaLevel,C_vegaLevel_default," ARM_ERR: precision: numeric expected",C_result);
        
		if((liborTypeId = ARM_ConvIrIndName (C_liborType, C_result)) == ARM_DEFAULT_ERR)
		{
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}
		
		if((volTypeId = ARM_ConvShapeType (C_VolType, C_result)) == ARM_DEFAULT_ERR)
		{
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}
		
		if((calibId = ARM_ConvYesOrNo (C_Calib, C_result)) == ARM_DEFAULT_ERR)
		{
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}
		
		if (C_pf2Id == "NULL")
			pf2Id = ARM_NULL_OBJECT;
		else
			pf2Id = LocalGetNumObjectId (C_pf2Id);
		if (C_pf3Id == "NULL")
			pf3Id = ARM_NULL_OBJECT;
		else
			pf3Id = LocalGetNumObjectId (C_pf3Id);
		
		if (C_marketmodelId == "NULL")
			marketmodelId = ARM_NULL_OBJECT;
		else
			marketmodelId = LocalGetNumObjectId (C_marketmodelId);
		
		long flagPreInit ;
		if((flagPreInit = ARM_ConvYesOrNo (C_flagPreInit, C_result)) == ARM_DEFAULT_ERR)
		{
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}
		
		long retCode;
		long objId;
		CCString prevClass;
		
		CCString curClass = LOCAL_FRMMODEL_CLASS;
		CCString stringId = GetLastCurCellEnvValue();
		
		if (!stringId)
		{
			retCode = ARMLOCAL_BootstCalibFRMModel(LocalGetNumObjectId (C_zcId),
				LocalGetNumObjectId (C_pf1Id),
				pf2Id,
				pf3Id,
				marketmodelId,
				liborTypeId,
				C_CalibParams,
				C_initCurve,
				C_mdec,
				C_meanRev,
				(long) C_nbfactor,
				C_nbrows,
				C_nbcolumns,
				C_correlmatrix,
				volTypeId,
				calibId,
				flagPreInit,
				C_precision,
                C_vegaLevel,
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
				retCode = ARMLOCAL_BootstCalibFRMModel(LocalGetNumObjectId (C_zcId),
					LocalGetNumObjectId (C_pf1Id),
					pf2Id,
					pf3Id,
					marketmodelId,
					liborTypeId,
					C_CalibParams,
					C_initCurve,
					C_mdec,
					C_meanRev,
					(long) C_nbfactor,
					C_nbrows,
					C_nbcolumns,
					C_correlmatrix,
					volTypeId,
					calibId,
					flagPreInit,
					C_precision,
                    C_vegaLevel,
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
				
				retCode = ARMLOCAL_BootstCalibFRMModel(LocalGetNumObjectId (C_zcId),
					LocalGetNumObjectId (C_pf1Id),
					pf2Id,
					pf3Id,
					marketmodelId,
					liborTypeId,
					C_CalibParams,
					C_initCurve,
					C_mdec,
					C_meanRev,
					(long) C_nbfactor,
					C_nbrows,
					C_nbcolumns,
					C_correlmatrix,
					volTypeId,
					calibId,
					flagPreInit,
					C_precision,
                    C_vegaLevel,
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_ARM_BootstCalibFRMModel" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


//___________________________________________________________________________________________

__declspec(dllexport) LPXLOPER WINAPI Local_PXL_ARM_BootstCalibFRMModel(LPXLOPER XL_zcId,																		
																		LPXLOPER XL_liborType,
																		LPXLOPER XL_CalibParams,
																		LPXLOPER XL_initCurve,
																		LPXLOPER XL_mdec,
																		LPXLOPER XL_meanRev,
                                                                        LPXLOPER XL_pf1Id,
																		LPXLOPER XL_nbfactor,										
																	    LPXLOPER XL_correlmatrix,
                                                                        LPXLOPER XL_pf2Id,
                                                                        LPXLOPER XL_VolType,
																		LPXLOPER XL_Calib,
                                                                        LPXLOPER XL_pf3Id,
                                                                        LPXLOPER XL_marketmodelId,
                                                                        LPXLOPER XL_flagPreInit,
                                                                        LPXLOPER XL_precision,
                                                                        LPXLOPER XL_vegaLevel)
{
	ADD_LOG("Local_PXL_ARM_BootstCalibFRMModel");
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();

	// C variable
	CCString C_zcId;
	CCString C_volId;
	CCString C_pf1Id;
	CCString C_pf2Id;
    CCString C_pf3Id;
    CCString C_marketmodelId;

    long pf2Id;
    long pf3Id;
    long marketmodelId;

	CCString C_liborType;

    CCString C_VolType;
    long volTypeId;

    CCString C_Calib;
    long calibId;

	long liborTypeId;

	VECTOR<double> C_CalibParams;

	VECTOR<double> C_initCurve;
	VECTOR<double> C_mdec;
	
	double C_meanRev;
    	
    double C_nbfactor;
    CCString C_flagPreInit;
      
    VECTOR<double> C_correlmatrix;
	long C_nbrows;
	long C_nbcolumns;

	// error
	static int error;
	static char* reason = "";

    XL_readStrCell(XL_zcId, C_zcId," ARM_ERR: Zc id: object expected",C_result);
    XL_readStrCell(XL_liborType, C_liborType," ARM_ERR: libor type: string expected",C_result);
    XL_readNumVector(XL_CalibParams, C_CalibParams," ARM_ERR: Calib Params: vector expected",C_result);
    XL_readNumVector(XL_initCurve, C_initCurve," ARM_ERR: init Curve: vector expected",C_result);
	XL_readNumVector(XL_mdec, C_mdec," ARM_ERR: init shift or LogNorProbability Curve: vector expected",C_result);
	XL_readNumCell(XL_meanRev,C_meanRev," ARM_ERR: Mean Reversion: numeric expected",C_result);	   
    XL_readStrCell(XL_pf1Id, C_pf1Id," ARM_ERR: Pf1 id: object expected",C_result);
	XL_readNumCell(XL_nbfactor,C_nbfactor," ARM_ERR:Nb of factors: numeric expected",C_result);
    
    if(C_nbfactor > 1)
        XL_readNumVectorAndSize(XL_correlmatrix,C_nbrows,C_nbcolumns,C_correlmatrix," ARM_ERR: init Curve: matrix expected",C_result);

	XL_readStrCellWD(XL_pf2Id, C_pf2Id,"NULL"," ARM_ERR: Pf2 id: object expected",C_result);
    XL_readStrCellWD(XL_VolType,C_VolType,"ROW"," ARM_ERR: Vol Type: string expected",C_result);
    XL_readStrCellWD(XL_Calib,C_Calib,"YES"," ARM_ERR: Calib: string expected",C_result);
    XL_readStrCellWD(XL_pf3Id, C_pf3Id,"NULL"," ARM_ERR: PfToCalibrateShift id: object expected",C_result);
    XL_readStrCellWD(XL_marketmodelId, C_marketmodelId,"NULL"," ARM_ERR: Market Model id: object expected",C_result);
    XL_readStrCellWD(XL_flagPreInit,C_flagPreInit,"NO"," ARM_ERR: FlagPreInit: string expected",C_result);

    double C_precision;
    double C_precision_default = 1.0e-3;
    XL_readNumCellWD(XL_precision,C_precision,C_precision_default," ARM_ERR: precision: numeric expected",C_result);

    double C_vegaLevel;
	double C_vegaLevel_default = 1.0e-10;
	XL_readNumCellWD(XL_vegaLevel,C_vegaLevel,C_vegaLevel_default," ARM_ERR: precision: numeric expected",C_result);

    if((liborTypeId = ARM_ConvIrIndName (C_liborType, C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}
   
    if((volTypeId = ARM_ConvShapeType (C_VolType, C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

    if((calibId = ARM_ConvYesOrNo (C_Calib, C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

    if (C_pf2Id == "NULL")
        pf2Id = ARM_NULL_OBJECT;
    else
        pf2Id = LocalGetNumObjectId (C_pf2Id);

    if (C_pf3Id == "NULL")
        pf3Id = ARM_NULL_OBJECT;
    else
        pf3Id = LocalGetNumObjectId (C_pf3Id);

    if (C_marketmodelId == "NULL")
        marketmodelId = ARM_NULL_OBJECT;
    else
        marketmodelId = LocalGetNumObjectId (C_marketmodelId);

    long flagPreInit ;
    if((flagPreInit = ARM_ConvYesOrNo (C_flagPreInit, C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}


    long retCode;
	long objId;

	CCString curClass = LOCAL_FRMMODEL_CLASS;
	CCString stringId;
	
	retCode = ARMLOCAL_BootstCalibFRMModel(LocalGetNumObjectId (C_zcId),
										   LocalGetNumObjectId (C_pf1Id),
										   pf2Id,
                                           pf3Id,
                                           marketmodelId,
										   liborTypeId,
										   C_CalibParams,
										   C_initCurve,
										   C_mdec,
										   C_meanRev,
                                           (long) C_nbfactor,
										   C_nbrows,
										   C_nbcolumns,
                                           C_correlmatrix,
                                           volTypeId,
        								   calibId,
                                           flagPreInit,
                                           C_precision,
                                           C_vegaLevel,
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
	   PXL_ARM_ERR();
	}

//	ARM_END();
	
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_PXL_ARM_BootstCalibFRMModel" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}



__declspec(dllexport) LPXLOPER WINAPI Local_ARM_GETFROMFRMMODEL(LPXLOPER XL_model,
																LPXLOPER XL_param)
{
	ADD_LOG("Local_ARM_GETFROMFRMMODEL");
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
	CCString C_model;

	CCString C_param;

	// error
	static int error;
	static char* reason = "";

	XL_readStrCell(XL_model,C_model," ARM_ERR: model: object expected",C_result);
	XL_readStrCell(XL_param,C_param," ARM_ERR: param: string expected",C_result);

	long retCode;

	retCode = ARMLOCAL_GETFROMFRMMODEL(LocalGetNumObjectId(C_model),
									   C_param,
									   C_result);

	if(retCode == ARM_OK)
	{
		long nbrows = C_result.getLong();
		int nbcolumns;
		
		if ( (!(strcmp((const char*)C_param,"KCURVE"))) || (!(strcmp((const char*)C_param,"ACURVE")))
                                                        || (!(strcmp((const char*)C_param,"SHIFTCURVE")))
                                                        || (!(strcmp((const char*)C_param,"PROBA_LOG"))))
			nbcolumns = 3;
		else if ( (!(strcmp((const char*)C_param,"KMEAN"))) || (!(strcmp((const char*)C_param,"AMEAN")))
            || (!(strcmp((const char*)C_param,"ERROR"))) || (!(strcmp((const char*)C_param,"MEANREV"))))
			nbcolumns = 1;
	
		FreeCurCellErr ();
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_ARM_GETFROMFRMMODEL" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}

//___________________________________________________________________________________________

__declspec(dllexport) LPXLOPER WINAPI Local_ARM_BootstCalibFRMModelMixture(
																	LPXLOPER XL_N,
																	LPXLOPER XL_zcId,                                                 
																	LPXLOPER XL_liborType,
																	LPXLOPER XL_lambdas,
																	LPXLOPER XL_meanRev,
                                                                    LPXLOPER XL_pf1Id,
																	LPXLOPER XL_pf2Id,	
																	LPXLOPER XL_nbproducts,
                                                                    LPXLOPER XL_nbfactor,
																	LPXLOPER XL_correlmatrix,
                                                                    LPXLOPER XL_VolType,
																	LPXLOPER XL_LowerBound,
																	LPXLOPER XL_UpperBound,
                                                                    LPXLOPER XL_Calib)
{
	ADD_LOG("Local_ARM_BootstCalibFRMModelMixture");
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();

	// C variable
	double C_N;
	CCString C_zcId;
	CCString C_volId;
	CCString C_pf1Id;
	CCString C_pf2Id;
	long pf2Id;
	double C_nbproducts;

	CCString C_liborType;

    CCString C_VolType;
    long volTypeId;

    CCString C_Calib;
    long calibId;

	long liborTypeId;
    
	VECTOR<double> C_lambdas;
	double C_meanRev;
    
    double C_nbfactor;

	VECTOR<double> C_lowerBound;
	VECTOR<double> C_upperBound;
       
    VECTOR<double> C_correlmatrix;
	long C_nbrows;
	long C_nbcolumns;

	// error
	static int error;
	static char* reason = "";

	XL_readNumCell(XL_N, C_N," ARM_ERR: N: numeric expected",C_result);
    XL_readStrCell(XL_zcId, C_zcId," ARM_ERR: Zc id: object expected",C_result);
    XL_readStrCell(XL_liborType, C_liborType," ARM_ERR: libor type: string expected",C_result);
	XL_readNumVector(XL_lambdas, C_lambdas," ARM_ERR: lambda: vector expected",C_result);
	XL_readNumCell(XL_meanRev,C_meanRev," ARM_ERR: Mean Reversion: numeric expected",C_result);	   
    XL_readStrCell(XL_pf1Id, C_pf1Id," ARM_ERR: Pf1 id: object expected",C_result);
	XL_readStrCellWD(XL_pf2Id, C_pf2Id,"NULL"," ARM_ERR: Pf2 id: object expected",C_result);
	XL_readNumCell(XL_nbproducts,C_nbproducts," ARM_ERR:Nb of products: numeric expected",C_result);
	XL_readNumCell(XL_nbfactor,C_nbfactor," ARM_ERR:Nb of factors: numeric expected",C_result);
    XL_readNumVectorAndSize(XL_correlmatrix,C_nbrows,C_nbcolumns,C_correlmatrix," ARM_ERR: init Curve: matrix expected",C_result);

    XL_readStrCellWD(XL_VolType,C_VolType,"ROW"," ARM_ERR: Vol Type: string expected",C_result);
	VECTOR<double> defaultLowerBound;
	VECTOR<double> defaultUpperBound;
	// Bornes de minimisation par défaut
	for (int i = 0; i < C_N+1; ++i)
	{
		defaultLowerBound.push_back(0.0);
		defaultUpperBound.push_back(1.0);
	}
	XL_readNumVectorWD(XL_LowerBound, C_lowerBound, defaultLowerBound, " ARM_ERR: lower bound: vector expected",C_result);
	XL_readNumVectorWD(XL_UpperBound, C_upperBound, defaultUpperBound, " ARM_ERR: upper bound: vector expected",C_result);
    XL_readStrCellWD(XL_Calib,C_Calib,"YES"," ARM_ERR: Calib: string expected",C_result);

    if((liborTypeId = ARM_ConvIrIndName (C_liborType, C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}
   
    if((volTypeId = ARM_ConvShapeType (C_VolType, C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	if((calibId = ARM_ConvYesOrNo (C_Calib, C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	if (C_pf2Id == "NULL")
        pf2Id = ARM_NULL_OBJECT;
    else
        pf2Id = LocalGetNumObjectId (C_pf2Id);

    long retCode;
	long objId;
	CCString prevClass;

	CCString curClass = LOCAL_FRMMODELMIXTURE_CLASS;
	CCString stringId = GetLastCurCellEnvValue();
	
	if (!stringId)
	{
		retCode = ARMLOCAL_BootstCalibFRMModelMixture(
											   C_N,
											   LocalGetNumObjectId (C_zcId),
											   LocalGetNumObjectId (C_pf1Id),
											   pf2Id,
											   C_nbproducts,
											   liborTypeId,
											   C_lambdas,
											   C_meanRev,
											   (long) C_nbfactor,
                                               C_nbrows,
											   C_nbcolumns,
                                               C_correlmatrix,
                                               volTypeId,
											   C_lowerBound,
											   C_upperBound,
                                               calibId,
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
			retCode = ARMLOCAL_BootstCalibFRMModelMixture(
											   C_N,
											   LocalGetNumObjectId (C_zcId),
											   LocalGetNumObjectId (C_pf1Id),
											   pf2Id,
											   C_nbproducts,
											   liborTypeId,
											   C_lambdas,
											   C_meanRev,
											   (long) C_nbfactor,
                                               C_nbrows,
											   C_nbcolumns,
                                               C_correlmatrix,
                                               volTypeId,
											   C_lowerBound,
											   C_upperBound,
                                               calibId,
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

			retCode = ARMLOCAL_BootstCalibFRMModelMixture(
											   C_N,
											   LocalGetNumObjectId (C_zcId),
											   LocalGetNumObjectId (C_pf1Id),
											   pf2Id,
											   C_nbproducts,
											   liborTypeId,
											   C_lambdas,
											   C_meanRev,
											   (long) C_nbfactor,
                                               C_nbrows,
											   C_nbcolumns,
                                               C_correlmatrix,
                                               volTypeId,
											   C_lowerBound,
											   C_upperBound,
                                               calibId,
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_ARM_BootstCalibFRMModelMixture" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}

//___________________________________________________________________________________________


__declspec(dllexport) LPXLOPER WINAPI Local_ARM_GETFROMFRMMODELMIXTURE(LPXLOPER XL_model,
																	   LPXLOPER XL_modelNumber,	
																	   LPXLOPER XL_param)
{
	ADD_LOG("Local_ARM_GETFROMFRMMODELMIXTURE");
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
	CCString C_model;
	double C_modelNumber;
	CCString C_param;

	// error
	static int error;
	static char* reason = "";

	XL_readStrCell(XL_model,C_model," ARM_ERR: model: object expected",C_result);
	XL_readNumCell(XL_modelNumber, C_modelNumber," ARM_ERR: modelNumber: numeric expected",C_result);
	XL_readStrCell(XL_param,C_param," ARM_ERR: param: string expected",C_result);

	long retCode;

	retCode = ARMLOCAL_GETFROMFRMMODELMIXTURE(LocalGetNumObjectId(C_model),
									   (int) C_modelNumber,
									   C_param,
									   C_result);

	if(retCode == ARM_OK)
	{
		long nbrows = C_result.getLong();
		int nbcolumns;
		
		if (((!(strcmp((const char*)C_param,"ACURVE"))) || (!(strcmp((const char*)C_param,"SPREADCURVE")))))
			nbcolumns = 3;
		else if ((!(strcmp((const char*)C_param,"MEANREV"))))
			nbcolumns = 1;
	
		FreeCurCellErr ();
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_ARM_GETFROMFRMMODELMIXTURE" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}

//___________________________________________________________________________________________

__declspec(dllexport) LPXLOPER WINAPI Local_ARM_MCFRMMODEL(LPXLOPER XL_FRMModId,
													       LPXLOPER XL_nbTraj,
													       LPXLOPER XL_nbStepIn,
                                                           LPXLOPER XL_pricerType,
                                                           LPXLOPER XL_mcMethod) 
                                           
{
	ADD_LOG("Local_ARM_MCFRMMODEL");
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();

	// C variable
	CCString C_FRMModId;
    CCString CS_mcMethod;
	
    double C_nbTraj;
	double C_nbStepIn;

    CCString C_pricerType;
	long pricerTypeId;
    long mcMethod;

	// error
	static int error;
	static char* reason = "";

	XL_readStrCell(XL_FRMModId,C_FRMModId," ARM_ERR: Model id: string expected",C_result);
	XL_readNumCell(XL_nbTraj, C_nbTraj," ARM_ERR: number of traj. : numeric expected",C_result);
	XL_readNumCell(XL_nbStepIn, C_nbStepIn," ARM_ERR: Monte Carlo Method. : numeric expected",C_result);
    XL_readStrCellWD(XL_pricerType,C_pricerType,"FNMC"," ARM_ERR: Pricer type: string expected",C_result);
    XL_readStrCellWD(XL_mcMethod, CS_mcMethod, "SIMPLE", " ARM_ERR: MC method : string expected",C_result);


    if ( (pricerTypeId = ARM_ConvPricerType(C_pricerType, C_result)) == ARM_DEFAULT_ERR)
	{
		 ARM_ARG_ERR();

		 return(LPXLOPER)&XL_result;
	}

    if ((mcMethod = ARM_ConvMcMethod(CS_mcMethod, C_result)) == ARM_DEFAULT_ERR)
	{
	   ARM_ARG_ERR();
	   return (LPXLOPER)&XL_result;
	}


	long retCode;
	long objId;
	CCString prevClass;

	CCString curClass = LOCAL_FRMMCMODEL_CLASS;
	CCString stringId = GetLastCurCellEnvValue ();

	if (!stringId)
	{
        retCode = ARMLOCAL_MCFRMMODEL(LocalGetNumObjectId (C_FRMModId),
									  (long) C_nbTraj,
									  (long) C_nbStepIn,
                                      (long) pricerTypeId, 
                                      mcMethod,
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
			retCode = ARMLOCAL_MCFRMMODEL(LocalGetNumObjectId (C_FRMModId),
									        (long) C_nbTraj,
									        (long) C_nbStepIn,
                                            (long) pricerTypeId,
                                            mcMethod,
									        C_result,
								            objId);

			if ( retCode == ARM_OK )
			{
			   objId = C_result.getLong ();

			   stringId = LocalMakeObjectId (objId, curClass);
			}
		}
		else
		{
			FreeCurCellContent ();

            retCode = ARMLOCAL_MCFRMMODEL(LocalGetNumObjectId (C_FRMModId),
									      (long) C_nbTraj,
									      (long) C_nbStepIn,
                                          (long) pricerTypeId,
                                          mcMethod,
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_ARM_MCFRMMODEL" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}



__declspec(dllexport) LPXLOPER WINAPI Local_PXL_ARM_MCFRMMODEL(LPXLOPER XL_FRMModId,
															   LPXLOPER XL_nbTraj,
															   LPXLOPER XL_nbStepIn,
															   LPXLOPER XL_pricerType,
                                                               LPXLOPER XL_mcMethod)
{
	ADD_LOG("Local_PXL_ARM_MCFRMMODEL");
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();

	// C variable
	CCString C_FRMModId;
    CCString CS_mcMethod;

    double C_nbTraj;
	double C_nbStepIn;

    CCString C_pricerType;
	long pricerTypeId;
    long mcMethod;

	// error
	static int error;
	static char* reason = "";

	XL_readStrCell(XL_FRMModId,C_FRMModId," ARM_ERR: Model id: string expected",C_result);
	XL_readNumCell(XL_nbTraj, C_nbTraj," ARM_ERR: number of traj. : numeric expected",C_result);
	XL_readNumCell(XL_nbStepIn, C_nbStepIn," ARM_ERR: Monte Carlo Method. : numeric expected",C_result);
    XL_readStrCellWD(XL_pricerType,C_pricerType,"FNMC"," ARM_ERR: Pricer type: string expected",C_result);
    XL_readStrCellWD(XL_mcMethod, CS_mcMethod, "SIMPLE", " ARM_ERR: MC method : string expected",C_result);


    if ( (pricerTypeId = ARM_ConvPricerType(C_pricerType, C_result)) == ARM_DEFAULT_ERR)
	{
		 ARM_ARG_ERR();

		 return(LPXLOPER)&XL_result;
	}
    if ((mcMethod = ARM_ConvMcMethod(CS_mcMethod, C_result)) == ARM_DEFAULT_ERR)
	{
	   ARM_ARG_ERR();
	   return (LPXLOPER)&XL_result;
	}

	long retCode;
	long objId;

	CCString curClass = LOCAL_FRMMCMODEL_CLASS;
	CCString stringId;

    retCode = ARMLOCAL_MCFRMMODEL(LocalGetNumObjectId (C_FRMModId),
								  (long) C_nbTraj,
								  (long) C_nbStepIn,
                                  (long) pricerTypeId,
                                  mcMethod,
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
	   PXL_ARM_ERR();
	}

//	ARM_END();
	
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_PXL_ARM_MCFRMMODEL" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}

//___________________________________________________________________________________________

__declspec(dllexport) LPXLOPER WINAPI Local_ARM_MCFRMMODELMIXTURE(LPXLOPER XL_FRMModMixId,
													       LPXLOPER XL_nbTraj,
													       LPXLOPER XL_nbStepIn,
                                                           LPXLOPER XL_pricerType) 
                                           
{
	ADD_LOG("Local_ARM_MCFRMMODELMIXTURE");
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();

	// C variable
	CCString C_FRMModMixId;
	
	double C_nbTraj;

	double C_nbStepIn;

    CCString C_pricerType;
	long pricerTypeId;

	// error
	static int error;
	static char* reason = "";

	XL_readStrCell(XL_FRMModMixId,C_FRMModMixId," ARM_ERR: Model Mixture id: string expected",C_result);
	XL_readNumCell(XL_nbTraj, C_nbTraj," ARM_ERR: number of traj. : numeric expected",C_result);
	XL_readNumCell(XL_nbStepIn, C_nbStepIn," ARM_ERR: Monte Carlo Method. : numeric expected",C_result);
    XL_readStrCellWD(XL_pricerType,C_pricerType,"FNMC"," ARM_ERR: Pricer type: string expected",C_result);

    if ( (pricerTypeId = ARM_ConvPricerType(C_pricerType, C_result)) == ARM_DEFAULT_ERR)
	{
		 ARM_ARG_ERR();

		 return(LPXLOPER)&XL_result;
	}


	long retCode;
	long objId;
	CCString prevClass;

	CCString curClass = LOCAL_FRMMCMODELMIXTURE_CLASS;
	CCString stringId = GetLastCurCellEnvValue ();

	if (!stringId)
	{
        retCode = ARMLOCAL_MCFRMMODELMIXTURE(LocalGetNumObjectId (C_FRMModMixId),
									  (long) C_nbTraj,
									  (long) C_nbStepIn,
                                      (long) pricerTypeId,                                      
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
			retCode = ARMLOCAL_MCFRMMODELMIXTURE(LocalGetNumObjectId (C_FRMModMixId),
									        (long) C_nbTraj,
									        (long) C_nbStepIn,
                                            (long) pricerTypeId,
									        C_result,
								            objId);

			if ( retCode == ARM_OK )
			{
			   objId = C_result.getLong ();

			   stringId = LocalMakeObjectId (objId, curClass);
			}
		}
		else
		{
			FreeCurCellContent ();

            retCode = ARMLOCAL_MCFRMMODELMIXTURE(LocalGetNumObjectId (C_FRMModMixId),
									      (long) C_nbTraj,
									      (long) C_nbStepIn,
                                          (long) pricerTypeId,
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_ARM_MCFRMMODELMIXTURE" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}



__declspec(dllexport) LPXLOPER WINAPI Local_PXL_ARM_MCFRMMODELMIXTURE(LPXLOPER XL_FRMModMixId,
															   LPXLOPER XL_nbTraj,
															   LPXLOPER XL_nbStepIn,
															   LPXLOPER XL_pricerType)
{
	ADD_LOG("Local_PXL_ARM_MCFRMMODELMIXTURE");
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();

	// C variable
	CCString C_FRMModMixId;

    double C_nbTraj;

	double C_nbStepIn;

    CCString C_pricerType;
	long pricerTypeId;

	// error
	static int error;
	static char* reason = "";

	XL_readStrCell(XL_FRMModMixId,C_FRMModMixId," ARM_ERR: Model id: string expected",C_result);
	XL_readNumCell(XL_nbTraj, C_nbTraj," ARM_ERR: number of traj. : numeric expected",C_result);
	XL_readNumCell(XL_nbStepIn, C_nbStepIn," ARM_ERR: Monte Carlo Method. : numeric expected",C_result);
    XL_readStrCellWD(XL_pricerType,C_pricerType,"RNMC"," ARM_ERR: Pricer type: string expected",C_result);

    if ( (pricerTypeId = ARM_ConvPricerType(C_pricerType, C_result)) == ARM_DEFAULT_ERR)
	{
		 ARM_ARG_ERR();

		 return(LPXLOPER)&XL_result;
	}

	long retCode;
	long objId;

	CCString curClass = LOCAL_FRMMCMODEL_CLASS;
	CCString stringId;

    retCode = ARMLOCAL_MCFRMMODELMIXTURE(LocalGetNumObjectId (C_FRMModMixId),
								  (long) C_nbTraj,
								  (long) C_nbStepIn,
                                  (long) pricerTypeId,                                      
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
	   PXL_ARM_ERR();
	}

//	ARM_END();
	
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_PXL_ARM_MCFRMMODELMIXTURE" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


//______________________________________________________________________________________________

__declspec(dllexport) LPXLOPER WINAPI Local_ARM_FRMMARKOVTREE(LPXLOPER XL_startDate,
                                                              LPXLOPER XL_horizon,
                                                              LPXLOPER XL_zc,                                                              
                                                              LPXLOPER XL_pathNumber,
                                                              LPXLOPER XL_Params,
                                                              LPXLOPER XL_FRMModel)                                                             
{
	ADD_LOG("Local_ARM_FRMMARKOVTREE");
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();

	// C variable
	double C_startDate;
	double C_horizon;

	CCString C_zc;
   	VECTOR<double> C_Params; 

    CCString C_FRMModel;
    double C_pathNumber;

	// error
	static int error;
	static char* reason = "";

	XL_readNumCell(XL_startDate,C_startDate," ARM_ERR: start date: numeric expected",C_result);
	XL_readNumCell(XL_horizon,C_horizon," ARM_ERR: start date: numeric expected",C_result);
	XL_readStrCell(XL_zc,C_zc," ARM_ERR: zc: object expected",C_result);
	XL_readNumCell(XL_pathNumber,C_pathNumber," ARM_ERR: pathNumber: numeric expected",C_result);
	XL_readNumVector(XL_Params,C_Params," ARM_ERR: Params: array of numeric expected",C_result);
    XL_readStrCellWD(XL_FRMModel,C_FRMModel,"DEFAULT"," ARM_ERR: FRMModel: object expected",C_result);
        
	long retCode;
	long objId;
	CCString prevClass;
	
	CCString curClass = LOCAL_FRM_MARKOVTREE_CLASS;
	CCString stringId = GetLastCurCellEnvValue ();

	if (!stringId)
	{
		retCode = ARMLOCAL_FRMMARKOVTREE(C_startDate,
                                         C_horizon,
                                         LocalGetNumObjectId(C_zc),                                         
                                         (long)C_pathNumber,
                                         C_Params,
                                         LocalGetNumObjectId(C_FRMModel),                                                                                         
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

		if(curClass == prevClass)
		{
		    retCode = ARMLOCAL_FRMMARKOVTREE(C_startDate,
                                             C_horizon,
                                             LocalGetNumObjectId(C_zc),                                         
                                             (long)C_pathNumber,
                                             C_Params,
                                             LocalGetNumObjectId(C_FRMModel),
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

		    retCode = ARMLOCAL_FRMMARKOVTREE(C_startDate,
                                             C_horizon,
                                             LocalGetNumObjectId(C_zc),                                         
                                             (long)C_pathNumber,
                                             C_Params,
                                             LocalGetNumObjectId(C_FRMModel),
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

//	ARM_END();
	
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_ARM_FRMMARKOVTREE" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


__declspec(dllexport) LPXLOPER WINAPI Local_PXL_ARM_FRMMARKOVTREE(LPXLOPER XL_startDate,
                                                              LPXLOPER XL_horizon,
                                                              LPXLOPER XL_zc,                                                              
                                                              LPXLOPER XL_pathNumber,
                                                              LPXLOPER XL_Params,
                                                              LPXLOPER XL_FRMModel)
{
	ADD_LOG("Local_PXL_ARM_FRMMARKOVTREE");
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();

	// C variable
	double C_startDate;
	double C_horizon;

	CCString C_zc;   	
   	VECTOR<double> C_Params; 
    CCString C_FRMModel;    
    double C_pathNumber;   
   
	// error
	static int error;
	static char* reason = "";

	XL_readNumCell(XL_startDate,C_startDate," ARM_ERR: start date: numeric expected",C_result);
	XL_readNumCell(XL_horizon,C_horizon," ARM_ERR: start date: numeric expected",C_result);
	XL_readStrCell(XL_zc,C_zc," ARM_ERR: zc: object expected",C_result);
	XL_readNumCell(XL_pathNumber,C_pathNumber," ARM_ERR: pathNumber: numeric expected",C_result);
	XL_readNumVector(XL_Params,C_Params," ARM_ERR: Params: array of numeric expected",C_result);
    XL_readStrCell(XL_FRMModel,C_FRMModel," ARM_ERR: FRMModel: object expected",C_result);
    
	long retCode;
	long objId;

	CCString curClass = LOCAL_FRM_MARKOVTREE_CLASS;
	CCString stringId;

	retCode = ARMLOCAL_FRMMARKOVTREE(C_startDate,
                                     C_horizon,
                                     LocalGetNumObjectId(C_zc),                                         
                                     (long)C_pathNumber,
                                     C_Params,
                                     LocalGetNumObjectId(C_FRMModel),
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
	   PXL_ARM_ERR();
	}

//	ARM_END();
	
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_PXL_ARM_FRMMARKOVTREE" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


__declspec(dllexport) LPXLOPER WINAPI Local_ARM_LOGDECALANA(LPXLOPER XL_zc,
															LPXLOPER XL_freq,
															LPXLOPER XL_resetMatu,
															LPXLOPER XL_shifts,
															LPXLOPER XL_fwdVols,
															LPXLOPER XL_isWeighted)
{
	ADD_LOG("Local_ARM_LOGDECALANA");
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();

	// C variable
	CCString C_zc;

	CCString C_freq;
	long freqId;

   	VECTOR<double> C_resetMatu; 
   	VECTOR<double> C_shifts; 
   	VECTOR<double> C_fwdVols; 

	CCString C_isWeighted;
	long isWeightedId;

	// error
	static int error;
	static char* reason = "";

	XL_readStrCell(XL_zc,C_zc," ARM_ERR: zc: object expected",C_result);
	XL_readStrCell(XL_freq,C_freq," ARM_ERR: frequency: string expected",C_result);
	XL_readNumVector(XL_resetMatu,C_resetMatu," ARM_ERR: reset Maturities: array of numeric expected",C_result);
	XL_readNumVector(XL_shifts,C_shifts," ARM_ERR: shifts: array of numeric expected",C_result);
	XL_readNumVector(XL_fwdVols,C_fwdVols," ARM_ERR: fwd Vols: array of numeric expected",C_result);
	XL_readStrCellWD(XL_isWeighted,C_isWeighted,"Y"," ARM_ERR: is Weighted ?: string expected",C_result);

	if((freqId = ARM_ConvFrequency (C_freq, C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	if((isWeightedId = ARM_ConvYesOrNo (C_isWeighted, C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	long retCode;
	long objId;
	CCString prevClass;
	
	CCString curClass = LOCAL_LOGDECANA_CLASS;
	CCString stringId = GetLastCurCellEnvValue ();

	if (!stringId)
	{
		retCode = ARMLOCAL_LOGDECANA(LocalGetNumObjectId(C_zc),
									 freqId,
									 C_resetMatu,
									 C_shifts,
									 C_fwdVols,
									 isWeightedId,
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

		if(curClass == prevClass)
		{
			retCode = ARMLOCAL_LOGDECANA(LocalGetNumObjectId(C_zc),
										 freqId,
										 C_resetMatu,
										 C_shifts,
										 C_fwdVols,
										 isWeightedId,
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

			retCode = ARMLOCAL_LOGDECANA(LocalGetNumObjectId(C_zc),
										 freqId,
										 C_resetMatu,
										 C_shifts,
										 C_fwdVols,
										 isWeightedId,
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

//	ARM_END();
	
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_ARM_LOGDECALANA" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


__declspec(dllexport) LPXLOPER WINAPI Local_PXL_ARM_LOGDECALANA(LPXLOPER XL_zc,
																LPXLOPER XL_freq,
																LPXLOPER XL_resetMatu,
																LPXLOPER XL_shifts,
																LPXLOPER XL_fwdVols,
																LPXLOPER XL_isWeighted)
{
	ADD_LOG("Local_PXL_ARM_LOGDECALANA");
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();

	// C variable
	CCString C_zc;

	CCString C_freq;
	long freqId;

   	VECTOR<double> C_resetMatu; 
   	VECTOR<double> C_shifts; 
   	VECTOR<double> C_fwdVols; 

	CCString C_isWeighted;
	long isWeightedId;

	// error
	static int error;
	static char* reason = "";

	XL_readStrCell(XL_zc,C_zc," ARM_ERR: zc: object expected",C_result);
	XL_readStrCell(XL_freq,C_freq," ARM_ERR: frequency: string expected",C_result);
	XL_readNumVector(XL_resetMatu,C_resetMatu," ARM_ERR: reset Maturities: array of numeric expected",C_result);
	XL_readNumVector(XL_shifts,C_shifts," ARM_ERR: shifts: array of numeric expected",C_result);
	XL_readNumVector(XL_fwdVols,C_fwdVols," ARM_ERR: fwd Vols: array of numeric expected",C_result);
	XL_readStrCellWD(XL_isWeighted,C_isWeighted,"Y"," ARM_ERR: is Weighted ?: string expected",C_result);

	if((freqId = ARM_ConvFrequency (C_freq, C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	if((isWeightedId = ARM_ConvYesOrNo (C_isWeighted, C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	long retCode;
	long objId;

	CCString curClass = LOCAL_LOGDECANA_CLASS;
	CCString stringId;

	retCode = ARMLOCAL_LOGDECANA(LocalGetNumObjectId(C_zc),
								 freqId,
								 C_resetMatu,
								 C_shifts,
								 C_fwdVols,
								 isWeightedId,
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
	   PXL_ARM_ERR();
	}

//	ARM_END();
	
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_PXL_ARM_LOGDECALANA" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


__declspec(dllexport) LPXLOPER WINAPI Local_ARM_QMODEL(LPXLOPER XL_date,
													   LPXLOPER XL_spot,
													   LPXLOPER XL_dividend,
													   LPXLOPER XL_discrate,
													   LPXLOPER XL_volat,
													   LPXLOPER XL_q0,
													   LPXLOPER XL_q1,
													   LPXLOPER XL_typstk)
{
	ADD_LOG("Local_ARM_QMODEL");


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
		double C_spot;
		
		double C_dividend_double;
		CCString C_dividend_str;
		long dividend_type;
		
		double C_discrate_double;
		CCString C_discrate_str;
		long discrate_type;
		
		double C_volat_double;
		CCString C_volat_str;
		long volat_type;
		
		CCString C_typstk;
		long typstk = K_YIELD;

		double C_q0;
		double C_q0_default = 0.0;

		double C_q1;
		double C_q1_default = 0.0;

		// error
		static int error;
		static char* reason = "";

		XL_readNumCell(XL_date,C_date," ARM_ERR: date: date expected",C_result);
		XL_readNumCell(XL_spot,C_spot," ARM_ERR: spot: numeric expected",C_result);

		XL_readStrOrNumCell(XL_dividend,C_dividend_str,C_dividend_double,dividend_type," ARM_ERR: dividend: string or numeric expected",C_result);
		XL_readStrOrNumCell(XL_discrate,C_discrate_str,C_discrate_double,discrate_type," ARM_ERR: discrate: string or numeric expected",C_result);
		XL_readStrOrNumCell(XL_volat,C_volat_str,C_volat_double,volat_type," ARM_ERR: volatility: string or numeric expected",C_result);

		XL_readNumCellWD(XL_q0,C_q0,C_q0_default," ARM_ERR: q0: numeric expected",C_result);
		XL_readNumCellWD(XL_q1,C_q1,C_q1_default," ARM_ERR: q1: numeric expected",C_result);
		XL_readStrCellWD(XL_typstk,C_typstk,"YIELD"," ARM_ERR: strike type: string expected",C_result);

		if ((typstk = ARM_ConvPriceYield (C_typstk, C_result)) == ARM_DEFAULT_ERR)
		{
		   ARM_ARG_ERR();
		   return (LPXLOPER)&XL_result;
		}
		
		if ( dividend_type == XL_TYPE_STRING )
		{
		   C_dividend_double = (double) LocalGetNumObjectId(C_dividend_str);

		   dividend_type = 1;
		}
		else
		{
		   dividend_type = 0;
		}

		if ( discrate_type == XL_TYPE_STRING )
		{
		   C_discrate_double = (double) LocalGetNumObjectId (C_discrate_str);
		   discrate_type = 1;
		}
		else
		{
		   discrate_type = 0;
		}

		if ( volat_type == XL_TYPE_STRING )
		{
		   C_volat_double = (double) LocalGetNumObjectId (C_volat_str);

		   volat_type = 1;
		}
		else
		{
		   volat_type = 0;
		}

		long retCode;
		long objId;
		CCString prevClass;
		
		CCString curClass = LOCAL_QMODEL_CLASS;
		CCString stringId = GetLastCurCellEnvValue ();

		if (!stringId)
		{
			retCode = ARMLOCAL_QMODEL (C_date,
									   C_spot,
									   dividend_type,
									   C_dividend_double,
									   discrate_type,
									   C_discrate_double,
									   volat_type,
									   C_volat_double,
									   C_q0,
									   C_q1,
									   typstk,
									   C_result);

			if(retCode == ARM_OK)
			{
				objId = C_result.getLong ();

				LocalSetCurCellEnvValue (curClass, objId); 

				stringId = LocalMakeObjectId (objId, curClass);
			}
			else
				retCode = ARM_KO;
		}
		else
		{
			prevClass = LocalGetStringObjectClass (stringId);
			
			objId = LocalGetNumObjectId (stringId);
				
			if (curClass == prevClass)
			{
				retCode = ARMLOCAL_QMODEL (C_date,
										   C_spot,
										   dividend_type,
										   C_dividend_double,
										   discrate_type,
										   C_discrate_double,
										   volat_type,
										   C_volat_double,
										   C_q0,
										   C_q1,
										   typstk,
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

				retCode = ARMLOCAL_QMODEL (C_date,
										   C_spot,
										   dividend_type,
										   C_dividend_double,
										   discrate_type,
										   C_discrate_double,
										   volat_type,
										   C_volat_double,
										   C_q0,
										   C_q1,
										   typstk,
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_ARM_QMODEL" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


__declspec(dllexport) LPXLOPER WINAPI Local_PXL_ARM_QMODEL(LPXLOPER XL_date,
														   LPXLOPER XL_spot,
														   LPXLOPER XL_dividend,
														   LPXLOPER XL_discrate,
														   LPXLOPER XL_volat,
														   LPXLOPER XL_q0,
														   LPXLOPER XL_q1,
														   LPXLOPER XL_typstk)
{
	ADD_LOG("Local_PXL_ARM_QMODEL");
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
	double C_spot;
	
	double C_dividend_double;
	CCString C_dividend_str;
	long dividend_type;
	
	double C_discrate_double;
	CCString C_discrate_str;
	long discrate_type;
	
	double C_volat_double;
	CCString C_volat_str;
	long volat_type;
	
	CCString C_typstk;
	long typstk = K_YIELD;
	
	double C_q0;
	double C_q0_default = 0.0;

	double C_q1;
	double C_q1_default = 0.0;

	// error
	static int error;
	static char* reason = "";

	XL_readNumCell(XL_date,C_date," ARM_ERR: date: date expected",C_result);
	XL_readNumCell(XL_spot,C_spot," ARM_ERR: spot: numeric expected",C_result);

	XL_readStrOrNumCell(XL_dividend,C_dividend_str,C_dividend_double,dividend_type," ARM_ERR: dividend: string or numeric expected",C_result);
	XL_readStrOrNumCell(XL_discrate,C_discrate_str,C_discrate_double,discrate_type," ARM_ERR: discrate: string or numeric expected",C_result);
	XL_readStrOrNumCell(XL_volat,C_volat_str,C_volat_double,volat_type," ARM_ERR: volatility: string or numeric expected",C_result);

	XL_readNumCellWD(XL_q0,C_q0,C_q0_default," ARM_ERR: q0: numeric expected",C_result);
	XL_readNumCellWD(XL_q1,C_q1,C_q1_default," ARM_ERR: q1: numeric expected",C_result);
	XL_readStrCellWD(XL_typstk,C_typstk,"YIELD"," ARM_ERR: strike type: string expected",C_result);

	if ((typstk = ARM_ConvPriceYield (C_typstk, C_result)) == ARM_DEFAULT_ERR)
	{
	   ARM_ARG_ERR();

	   return (LPXLOPER)&XL_result;
	}
	
	if ( dividend_type == XL_TYPE_STRING )
	{
	   C_dividend_double = (double) LocalGetNumObjectId(C_dividend_str);

	   dividend_type = 1;
	}
	else
	{
	   dividend_type = 0;
	}

	if ( discrate_type == XL_TYPE_STRING )
	{
	   C_discrate_double = (double)LocalGetNumObjectId (C_discrate_str);
	   discrate_type = 1;
	}
	else
	{
	   discrate_type = 0;
	}

	if ( volat_type == XL_TYPE_STRING )
	{
	   C_volat_double = (double) LocalGetNumObjectId (C_volat_str);

	   volat_type = 1;
	}
	else
	{
	   volat_type = 0;
	}

	long retCode;
	long objId;
	
	CCString curClass = LOCAL_QMODEL_CLASS;
	CCString stringId;

	retCode = ARMLOCAL_QMODEL (C_date,
							   C_spot,
							   dividend_type,
							   C_dividend_double,
							   discrate_type,
							   C_discrate_double,
							   volat_type,
							   C_volat_double,
							   C_q0,
							   C_q1,
							   typstk,
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
	   PXL_ARM_ERR();
	}

//	ARM_END();
	
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_PXL_ARM_QMODEL" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


__declspec(dllexport) LPXLOPER WINAPI Local_GLOBDFBS(LPXLOPER XL_DomBSId,
													 LPXLOPER XL_DomCurrId,
													 LPXLOPER XL_FrgBSId,
													 LPXLOPER XL_FrgCurrId,
													 LPXLOPER XL_fxVolCrvId,
													 LPXLOPER XL_FFxCorrId,
													 LPXLOPER XL_RatesCorrId,
													 LPXLOPER XL_fxVolModelId)
{
	ADD_LOG("Local_GLOBDFBS");
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
		ARM_NOCALCIFWIZ();

		// C variable
		CCString C_DomBSId;
		CCString C_DomCurrId;
		CCString C_FrgBSId;
		CCString C_FrgCurrId;
		CCString C_fxVolCrvId;
		CCString C_FFxCorrId;
		CCString C_RatesCorrId;
		CCString C_fxVolModel;
		long fxVolModelId;
		
		// error
		static int error;
		static char* reason = "";

		XL_readStrCell(XL_DomBSId, C_DomBSId, " ARM_ERR: Modele BS Domestique id: string expected",C_result);
		XL_readStrCell(XL_DomCurrId, C_DomCurrId, " ARM_ERR:  Currency Domestique id: string expected",C_result);
		XL_readStrCell(XL_FrgBSId, C_FrgBSId, " ARM_ERR:  Modele BS Foreign id: string expected",C_result);
		XL_readStrCell(XL_FrgCurrId, C_FrgCurrId, " ARM_ERR:  Currency Foreign id: string expected",C_result);
		XL_readStrCell(XL_fxVolCrvId, C_fxVolCrvId, " ARM_ERR:  Vol Curve FX id: string expected",C_result);
		XL_readStrCell(XL_FFxCorrId, C_FFxCorrId, " ARM_ERR:  Correl Foreign-FX id: string expected",C_result);
		XL_readStrCell(XL_RatesCorrId, C_RatesCorrId, " ARM_ERR:  Correl Foreign-Domestique id: string expected",C_result);
		XL_readStrCellWD(XL_fxVolModelId,C_fxVolModel,"DEFAULT"," ARM_ERR: fx vol model: object expected",C_result);

		if(C_fxVolModel == "DEFAULT")
		{
			fxVolModelId = ARM_NULL_OBJECT;
		}
		else
		{
			fxVolModelId = LocalGetNumObjectId (C_fxVolModel);
		}
			
		long retCode;
		long objId;
		CCString prevClass;
		
		CCString curClass = LOCAL_DFFXBS_CLASS;
		CCString stringId = GetLastCurCellEnvValue ();

		if (!stringId)
		{
			retCode = ARMLOCAL_GLOBDFBS(LocalGetNumObjectId (C_DomBSId),
										LocalGetNumObjectId (C_DomCurrId),
										LocalGetNumObjectId (C_FrgBSId),
										LocalGetNumObjectId (C_FrgCurrId),
										LocalGetNumObjectId (C_fxVolCrvId),
										LocalGetNumObjectId (C_FFxCorrId),
										LocalGetNumObjectId (C_RatesCorrId),
										fxVolModelId,
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
				retCode = ARMLOCAL_GLOBDFBS(LocalGetNumObjectId (C_DomBSId),
											LocalGetNumObjectId (C_DomCurrId),
											LocalGetNumObjectId (C_FrgBSId),
											LocalGetNumObjectId (C_FrgCurrId),
											LocalGetNumObjectId (C_fxVolCrvId),
											LocalGetNumObjectId (C_FFxCorrId),
											LocalGetNumObjectId (C_RatesCorrId),
											fxVolModelId,
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

				retCode = ARMLOCAL_GLOBDFBS(LocalGetNumObjectId (C_DomBSId),
											LocalGetNumObjectId (C_DomCurrId),
											LocalGetNumObjectId (C_FrgBSId),
											LocalGetNumObjectId (C_FrgCurrId),
											LocalGetNumObjectId (C_fxVolCrvId),
											LocalGetNumObjectId (C_FFxCorrId),
											LocalGetNumObjectId (C_RatesCorrId),
											fxVolModelId,
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_GLOBDFBS" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}



__declspec(dllexport) LPXLOPER WINAPI Local_PXL_GLOBDFBS(LPXLOPER XL_DomBSId,
														 LPXLOPER XL_DomCurrId,
														 LPXLOPER XL_FrgBSId,
														 LPXLOPER XL_FrgCurrId,
														 LPXLOPER XL_fxVolCrvId,
														 LPXLOPER XL_FFxCorrId,
														 LPXLOPER XL_RatesCorrId,
														 LPXLOPER XL_fxVolModelId)
{
	ADD_LOG("Local_PXL_GLOBDFBS");
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;
	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
		// C variable
		CCString C_DomBSId;
		CCString C_DomCurrId;
		CCString C_FrgBSId;
		CCString C_FrgCurrId;
		CCString C_fxVolCrvId;
		CCString C_FFxCorrId;
		CCString C_RatesCorrId;
		CCString C_fxVolModel;
		long fxVolModelId;
		
		// error
		static int error;
		static char* reason = "";

		XL_readStrCell(XL_DomBSId, C_DomBSId, " ARM_ERR: Modele BS Domestique id: string expected",C_result);
		XL_readStrCell(XL_DomCurrId, C_DomCurrId, " ARM_ERR:  Currency Domestique id: string expected",C_result);
		XL_readStrCell(XL_FrgBSId, C_FrgBSId, " ARM_ERR:  Modele BS Foreign id: string expected",C_result);
		XL_readStrCell(XL_FrgCurrId, C_FrgCurrId, " ARM_ERR:  Currency Foreign id: string expected",C_result);
		XL_readStrCell(XL_fxVolCrvId, C_fxVolCrvId, " ARM_ERR:  Vol Curve FX id: string expected",C_result);
		XL_readStrCell(XL_FFxCorrId, C_FFxCorrId, " ARM_ERR:  Correl Foreign-FX id: string expected",C_result);
		XL_readStrCell(XL_RatesCorrId, C_RatesCorrId, " ARM_ERR:  Correl Foreign-Domestique id: string expected",C_result);
		XL_readStrCellWD(XL_fxVolModelId,C_fxVolModel,"DEFAULT"," ARM_ERR: fx vol model: object expected",C_result);

		if(C_fxVolModel == "DEFAULT")
		{
			fxVolModelId = ARM_NULL_OBJECT;
		}
		else
		{
			fxVolModelId = LocalGetNumObjectId (C_fxVolModel);
		}

		long retCode;
		long objId;
		
		CCString curClass = LOCAL_DFFXBS_CLASS;
		CCString stringId;

		retCode = ARMLOCAL_GLOBDFBS(LocalGetNumObjectId (C_DomBSId),
									LocalGetNumObjectId (C_DomCurrId),
									LocalGetNumObjectId (C_FrgBSId),
									LocalGetNumObjectId (C_FrgCurrId),
									LocalGetNumObjectId (C_fxVolCrvId),
									LocalGetNumObjectId (C_FFxCorrId),
									LocalGetNumObjectId (C_RatesCorrId),
									fxVolModelId,
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
			PXL_ARM_ERR();
		}

//	ARM_END();
	
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_PXL_GLOBDFBS" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}



__declspec(dllexport) LPXLOPER WINAPI Local_ARM_CROSSMODEL (LPXLOPER XL_date,
															LPXLOPER XL_dBSmodel,
															LPXLOPER XL_fBSmodel,
															LPXLOPER XL_fxVol,
															LPXLOPER XL_dFxCorr,
															LPXLOPER XL_fFxCorr,
															LPXLOPER XL_domfgnCorr,
															LPXLOPER XL_discountZC,
															LPXLOPER XL_correlforAdj,
															LPXLOPER XL_adjustflag,
															LPXLOPER XL_slopeflag)
{
	ADD_LOG("Local_ARM_CROSSMODEL ");
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
	CCString C_dBSmodel;
	CCString C_fBSmodel;
	CCString C_domfgnCorr;
	CCString C_dFxCorr;
	CCString C_fFxCorr;
	CCString C_fxVol;
	CCString C_discountZc;
	long C_discCurveId;

	double C_correlforAdj;
	double C_correlforAdj_default = 0.0;
  
	double C_adjflag;
	double C_adjflag_default = 1.0;
  
    double C_slopeflag;
	double C_slopeflag_default = 0.0;  

	// error
	static int error;
	static char* reason = "";

	XL_readNumCell(XL_date,C_date," ARM_ERR: date: date expected",C_result);
	XL_readStrCell(XL_dBSmodel,C_dBSmodel," ARM_ERR: Pay Index B&S model id: object expected",C_result);
	XL_readStrCell(XL_fBSmodel,C_fBSmodel," ARM_ERR: RangeIndex B&S model id: object expected",C_result);
	XL_readStrCell(XL_fxVol,C_fxVol," ARM_ERR: forex volatility: object expected",C_result);
	XL_readStrCell(XL_dFxCorr,C_dFxCorr," ARM_ERR: pay idx forex correlation: object expected",C_result);
	XL_readStrCell(XL_fFxCorr,C_fFxCorr," ARM_ERR: range idx forex correlation: object expected",C_result);
	XL_readStrCell(XL_domfgnCorr,C_domfgnCorr," ARM_ERR: pay/range index correlation: object expected",C_result);
	XL_readStrCellWD(XL_discountZC,C_discountZc,"DEFAULT"," ARM_ERR: discount curve id: object expected",C_result);
	XL_readNumCellWD(XL_correlforAdj,C_correlforAdj,C_correlforAdj_default," ARM_ERR: correlation for convexity adjustment: numeric expected",C_result);
	XL_readNumCellWD(XL_adjustflag,C_adjflag,C_adjflag_default," ARM_ERR: adjustFlag: numeric expected",C_result);
	XL_readNumCellWD(XL_slopeflag,C_slopeflag,C_slopeflag_default," ARM_ERR: slopeFlag: numeric expected",C_result);


	if (C_discountZc == "DEFAULT")
		C_discCurveId = ARM_NULL_OBJECT;
	else
		C_discCurveId = LocalGetNumObjectId(C_discountZc);
	
	long retCode;
	long objId;
	CCString prevClass;

	CCString curClass = LOCAL_CROSSMOD_CLASS;
	CCString stringId = GetLastCurCellEnvValue ();

	if (!stringId)
	{
		retCode = ARMLOCAL_CROSSMODEL (	C_date,
										LocalGetNumObjectId (C_dBSmodel),
										LocalGetNumObjectId (C_fBSmodel),
										LocalGetNumObjectId (C_fxVol),
										LocalGetNumObjectId (C_dFxCorr),
										LocalGetNumObjectId (C_fFxCorr),
										LocalGetNumObjectId (C_domfgnCorr),
										C_discCurveId,
										C_correlforAdj,
										C_adjflag,
										C_slopeflag,
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
			
		if(curClass == prevClass)
		{
		
			retCode = ARMLOCAL_CROSSMODEL (	C_date,
											LocalGetNumObjectId (C_dBSmodel),
											LocalGetNumObjectId (C_fBSmodel),
											LocalGetNumObjectId (C_fxVol),
											LocalGetNumObjectId (C_dFxCorr),
											LocalGetNumObjectId (C_fFxCorr),
											LocalGetNumObjectId (C_domfgnCorr),
											C_discCurveId,
											C_correlforAdj,
											C_adjflag,
											C_slopeflag,
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
			
			retCode = ARMLOCAL_CROSSMODEL (	C_date,
											LocalGetNumObjectId (C_dBSmodel),
											LocalGetNumObjectId (C_fBSmodel),
											LocalGetNumObjectId (C_fxVol),
											LocalGetNumObjectId (C_dFxCorr),
											LocalGetNumObjectId (C_fFxCorr),
											LocalGetNumObjectId (C_domfgnCorr),
											C_discCurveId,
											C_correlforAdj,
											C_adjflag,
											C_slopeflag,
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

//	ARM_END();
	
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_ARM_CROSSMODEL" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}



__declspec(dllexport) LPXLOPER WINAPI Local_PXL_ARM_CROSSMODEL (LPXLOPER XL_date,
																LPXLOPER XL_dBSmodel,
																LPXLOPER XL_fBSmodel,
																LPXLOPER XL_fxVol,
																LPXLOPER XL_dFxCorr,
																LPXLOPER XL_fFxCorr,
																LPXLOPER XL_domfgnCorr,
																LPXLOPER XL_discountZC,
																LPXLOPER XL_correlforAdj,
																LPXLOPER XL_adjustflag,
																LPXLOPER XL_slopeflag)
{
	ADD_LOG("Local_PXL_ARM_CROSSMODEL ");
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
	CCString C_dBSmodel;
	CCString C_fBSmodel;
	CCString C_domfgnCorr;
	CCString C_dFxCorr;
	CCString C_fFxCorr;
	CCString C_fxVol;
	CCString C_discountZc;
	long C_discCurveId;

	double C_correlforAdj;
	double C_correlforAdj_default = 0.0;

	double C_adjflag;
	double C_adjflag_default = 1.0;
  
    double C_slopeflag;
	double C_slopeflag_default = 0.0;  

	// error
	static int error;
	static char* reason = "";

	XL_readNumCell(XL_date,C_date," ARM_ERR: date: date expected",C_result);
	XL_readStrCell(XL_dBSmodel,C_dBSmodel," ARM_ERR: Pay Index B&S model id: object expected",C_result);
	XL_readStrCell(XL_fBSmodel,C_fBSmodel," ARM_ERR: RangeIndex B&S model id: object expected",C_result);
	XL_readStrCell(XL_fxVol,C_fxVol," ARM_ERR: forex volatility: object expected",C_result);
	XL_readStrCell(XL_dFxCorr,C_dFxCorr," ARM_ERR: pay idx forex correlation: object expected",C_result);
	XL_readStrCell(XL_fFxCorr,C_fFxCorr," ARM_ERR: range idx forex correlation: object expected",C_result);
	XL_readStrCell(XL_domfgnCorr,C_domfgnCorr," ARM_ERR: pay/range index correlation: object expected",C_result);
	XL_readStrCellWD(XL_discountZC,C_discountZc,"DEFAULT"," ARM_ERR: discount curve id: object expected",C_result);
	XL_readNumCellWD(XL_correlforAdj,C_correlforAdj,C_correlforAdj_default," ARM_ERR: correlation for convexity adjustment: numeric expected",C_result);
	XL_readNumCellWD(XL_adjustflag,C_adjflag,C_adjflag_default," ARM_ERR: adjustFlag: numeric expected",C_result);
	XL_readNumCellWD(XL_slopeflag,C_slopeflag,C_slopeflag_default," ARM_ERR: slopeFlag: numeric expected",C_result);

	if (C_discountZc == "DEFAULT")
		C_discCurveId = ARM_NULL_OBJECT;
	else
		C_discCurveId = LocalGetNumObjectId(C_discountZc);

	long retCode;
	long objId;
	
	CCString curClass = LOCAL_CROSSMOD_CLASS;
	CCString stringId;

retCode = ARMLOCAL_CROSSMODEL (	C_date,
									LocalGetNumObjectId (C_dBSmodel),
									LocalGetNumObjectId (C_fBSmodel),
									LocalGetNumObjectId (C_fxVol),
									LocalGetNumObjectId (C_dFxCorr),
									LocalGetNumObjectId (C_fFxCorr),
									LocalGetNumObjectId (C_domfgnCorr),
									C_discCurveId,
									C_correlforAdj,
									C_adjflag,
									C_slopeflag,
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
	   PXL_ARM_ERR();
	}

//	ARM_END();
	
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_PXL_ARM_CROSSMODEL" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


/////////////////////////////////////////
/// Functor to create a bs convexity adjustment manager
/////////////////////////////////////////
class BSConvAdjFunc : public ARMResultLong2LongFunc
{
public:
	BSConvAdjFunc(int SUMMITFormulaeUsed, int UseSabrCMS)
    : itsSUMMITFormulaeUsed(SUMMITFormulaeUsed),
	  itsUseSabrCMS(UseSabrCMS)
	{};

	long operator()( ARM_result& result, long objId )
	{
		return ARMLOCAL_BSConvAdjust_Create(itsSUMMITFormulaeUsed, itsUseSabrCMS, result, objId);
	};
private:
    int itsSUMMITFormulaeUsed;
	int itsUseSabrCMS;
};


/////////////////////////////////////////
/// BS Convexity Adjust Manager create common function
/////////////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_BSConvAdjust_CreateCommon(
    LPXLOPER XL_SUMMITFormulaeUsed, LPXLOPER XL_UseSabrCMS,   
	bool PersistentInXL )
{
	ADD_LOG("Local_BSConvAdjust_CreateCommon");
	/// this is defined first because it is used in XL macros
	static XLOPER XL_result;
	/// Get the variables from the XLOper variables
	ARM_result C_result;

	/// to avoid computation if called by the wizard
	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();
	
    CCString C_SUMMITFormulaeUsedStr;
    int C_SUMMITFormulaeUsed;

	CCString C_UseSabrCMSStr;
    int C_UseSabrCMS;

	// error
	static int error;
	static char* reason = "";

    //XL_readStrCellWD(XL_SUMMITFormulaeUsed,C_SUMMITFormulaeUsedStr,"N"," ARM_ERR: SUMMITFormulaeUsed: numeric expected", C_result);
	XL_readStrCellWD(XL_SUMMITFormulaeUsed,C_SUMMITFormulaeUsedStr,"EXP"," ARM_ERR: SUMMITFormulaeUsed: numeric expected", C_result);

    if ( (C_SUMMITFormulaeUsed = ARM_ConvSummitFormulae (C_SUMMITFormulaeUsedStr, C_result)) == ARM_DEFAULT_ERR)
	{
	   ARM_ARG_ERR();

	   return (LPXLOPER)&XL_result;
	}

	XL_readStrCellWD(XL_UseSabrCMS,C_UseSabrCMSStr,"NO"," ARM_ERR: UseSabrCMSStr: numeric expected", C_result);

    if ( (C_UseSabrCMS = ARM_ConvYesOrNo (C_UseSabrCMSStr, C_result)) == ARM_DEFAULT_ERR)
	{
	   ARM_ARG_ERR();

	   return (LPXLOPER)&XL_result;
	}

	/// use the concept of Functor to transfer the knowledge of
	/// a function with a context
	BSConvAdjFunc ourFunc(C_SUMMITFormulaeUsed, C_UseSabrCMS);
	
	/// call the general function
	fillXL_Result( LOCAL_BSCONVADJUST_CLASS, ourFunc, C_result, XL_result, PersistentInXL );

	/// return the result as an LPXLOPER
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_BSConvAdjust_CreateCommon" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


//////////////////////////////////////////////////////
/// Addin to create a Convexity Adjust Manager
//////////////////////////////////////////////////////

__declspec(dllexport) LPXLOPER WINAPI Local_BSConvAdjust_Create(LPXLOPER XL_SUMMITFormulaeUsed,
																LPXLOPER XL_UseSabrCMS)
{
	ADD_LOG("Local_BSConvAdjust_Create");
	bool PersistentInXL = true;
	return Local_BSConvAdjust_CreateCommon(
        XL_SUMMITFormulaeUsed,
		XL_UseSabrCMS,
		PersistentInXL);
}



//////////////////////////////////////////////////////
/// Addin to create a Convexity Adjust Manager
/// Version for VBA
//////////////////////////////////////////////////////

__declspec(dllexport) LPXLOPER WINAPI Local_PXL_BSConvAdjust_Create(
    LPXLOPER XL_SUMMITFormulaeUsed, LPXLOPER XL_UseSabrCMS)
{
	ADD_LOG("Local_PXL_BSConvAdjust_Create");
	bool PersistentInXL = false;
	return Local_BSConvAdjust_CreateCommon(
        XL_SUMMITFormulaeUsed, XL_UseSabrCMS,
		PersistentInXL);
}


/////////////////////////////////////////
/// Functor to create a replication convexity adjustment manager
/////////////////////////////////////////
class ReplicConvAdjFunc : public ARMResultLong2LongFunc
{
public:
	ReplicConvAdjFunc(int Payoff_ReplicMode,
					  double Payoff_ReplicPrecision,
					  int Payoff_StopMode,
					  double Payoff_StopThreshold,
					  int Sensi_ReplicMode,
					  double Sensi_ReplicPrecision,
					  int Sensi_StopMode,
					  double Sensi_StopThreshold,
					  long UsedModel,
					  double StrikeMinReplic,
					  double StrikeMaxReplic)
					 : C_Payoff_ReplicMode(Payoff_ReplicMode),
					   C_Payoff_ReplicPrecision(Payoff_ReplicPrecision),
					   C_Payoff_StopMode(Payoff_StopMode),
					   C_Payoff_StopThreshold(Payoff_StopThreshold),
					   C_Sensi_ReplicMode(Sensi_ReplicMode),
					   C_Sensi_ReplicPrecision(Sensi_ReplicPrecision),
					   C_Sensi_StopMode(Sensi_StopMode),
					   C_Sensi_StopThreshold(Sensi_StopThreshold),
					   C_UsedModel(UsedModel),
					   C_StrikeMinReplic(StrikeMinReplic),
					   C_StrikeMaxReplic(StrikeMaxReplic)
	{};

	long operator()( ARM_result& result, long objId )
	{
		return ARMLOCAL_ReplicConvAdjust_Create(
			C_Payoff_ReplicMode,
			C_Payoff_ReplicPrecision,
			C_Payoff_StopMode,
			C_Payoff_StopThreshold,
			C_Sensi_ReplicMode,
			C_Sensi_ReplicPrecision,
			C_Sensi_StopMode,
			C_Sensi_StopThreshold,
			C_UsedModel,
			C_StrikeMinReplic,
			C_StrikeMaxReplic,
			result, 
			objId);
	};

private:
	int C_Payoff_ReplicMode;
	double C_Payoff_ReplicPrecision;
	int C_Payoff_StopMode;
	double C_Payoff_StopThreshold;
	int C_Sensi_ReplicMode;
	double C_Sensi_ReplicPrecision;
	int C_Sensi_StopMode;
	double C_Sensi_StopThreshold;
	long C_UsedModel;
	double C_StrikeMinReplic;
	double C_StrikeMaxReplic;

};


/////////////////////////////////////////
/// ReplicationConvexity Adjust Manager create common function
/////////////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_ReplicConvAdjust_CreateCommon(
	LPXLOPER XL_Payoff_ReplicMode,
	LPXLOPER XL_Payoff_ReplicPrecision,
	LPXLOPER XL_Payoff_StopMode,
	LPXLOPER XL_Payoff_StopThreshold,
	LPXLOPER XL_Sensi_ReplicMode,
	LPXLOPER XL_Sensi_ReplicPrecision,
	LPXLOPER XL_Sensi_StopMode,
	LPXLOPER XL_Sensi_StopThreshold,
	LPXLOPER XL_UsedModel,
	LPXLOPER XL_StrikeMinReplic,
	LPXLOPER XL_StrikeMaxReplic,
	bool PersistentInXL )
{
	ADD_LOG("Local_ReplicConvAdjust_CreateCommon");
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

	CCString C_Payoff_ReplicModeStr;
	int C_Payoff_ReplicMode;
	double C_Payoff_ReplicPrecision;
	CCString C_Payoff_StopModeStr;
	int C_Payoff_StopMode;
	double C_Payoff_StopThreshold;
	CCString C_Sensi_ReplicModeStr;
	int C_Sensi_ReplicMode;
	double C_Sensi_ReplicPrecision;
	CCString C_Sensi_StopModeStr;
	int C_Sensi_StopMode;
	double C_Sensi_StopThreshold;
	CCString C_UsedModel;
	long C_UsedModelId;
	double C_StrikeMinReplic;
	double C_StrikeMinReplic_Default = (double) K_STRIKEMIN_REPLIC;
	double C_StrikeMaxReplic_Default = (double) K_STRIKEMAX_REPLIC;
	double C_StrikeMaxReplic;

	XL_readStrCell(	 XL_Payoff_ReplicMode,		C_Payoff_ReplicModeStr,	" ARM_ERR: replicMode: string expected",				C_result);
	XL_readNumCell(	 XL_Payoff_ReplicPrecision,		C_Payoff_ReplicPrecision,	" ARM_ERR: replicPrecision: numeric expected",				C_result);
	XL_readStrCell(	 XL_Payoff_StopMode,		C_Payoff_StopModeStr,	" ARM_ERR: stopMode: string expected",				C_result);
	XL_readNumCell(	 XL_Payoff_StopThreshold,		C_Payoff_StopThreshold,	" ARM_ERR: stopThreshold: numeric expected",				C_result);

	XL_readStrCell(	 XL_Sensi_ReplicMode,		C_Sensi_ReplicModeStr,	" ARM_ERR: replicMode: string expected",				C_result);
	XL_readNumCell(	 XL_Sensi_ReplicPrecision,		C_Sensi_ReplicPrecision,	" ARM_ERR: replicPrecision: numeric expected",				C_result);
	XL_readStrCell(	 XL_Sensi_StopMode,		C_Sensi_StopModeStr,	" ARM_ERR: stopMode: string expected",				C_result);
	XL_readNumCell(	 XL_Sensi_StopThreshold,		C_Sensi_StopThreshold,	" ARM_ERR: stopThreshold: numeric expected",				C_result);

	XL_readStrCellWD(	 XL_UsedModel,		C_UsedModel, "NULL",	" ARM_ERR: Model Id: string expected",				C_result);

	XL_readNumCellWD(	 XL_StrikeMinReplic,		C_StrikeMinReplic,	C_StrikeMinReplic_Default, " ARM_ERR: StrikeMinReplic: numeric expected",				C_result);
	XL_readNumCellWD(	 XL_StrikeMaxReplic,		C_StrikeMaxReplic,	C_StrikeMaxReplic_Default, " ARM_ERR: StrikeMaxReplic: numeric expected",				C_result);

	if( ((C_Payoff_ReplicMode	= ARM_ConvReplicMode( C_Payoff_ReplicModeStr, C_result)) == ARM_DEFAULT_ERR)
		|| ((C_Payoff_StopMode	= ARM_ConvStopMode( C_Payoff_StopModeStr, C_result)) == ARM_DEFAULT_ERR)
		|| ((C_Sensi_ReplicMode	= ARM_ConvReplicMode( C_Sensi_ReplicModeStr, C_result)) == ARM_DEFAULT_ERR)
		|| ((C_Sensi_StopMode	= ARM_ConvStopMode( C_Sensi_StopModeStr, C_result)) == ARM_DEFAULT_ERR))
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	if (C_UsedModel == "NULL")
		C_UsedModelId = ARM_NULL_OBJECT;
	else
		C_UsedModelId = LocalGetNumObjectId (C_UsedModel);

	/// use the concept of Functor to transfer the knowledge of
	/// a function with a context
	ReplicConvAdjFunc ourFunc(
		C_Payoff_ReplicMode,
		C_Payoff_ReplicPrecision,
		C_Payoff_StopMode,
		C_Payoff_StopThreshold,
		C_Sensi_ReplicMode,
		C_Sensi_ReplicPrecision,
		C_Sensi_StopMode,
		C_Sensi_StopThreshold,
		C_UsedModelId,
		C_StrikeMinReplic,
		C_StrikeMaxReplic);
	
	/// call the general function
	fillXL_Result( LOCAL_REPLICCONVADJUST_CLASS, ourFunc, C_result, XL_result, PersistentInXL );

	/// return the result as an LPXLOPER
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_ReplicConvAdjust_CreateCommon" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


//////////////////////////////////////////////////////
/// Addin to create a Replication Convexity Adjust Manager
//////////////////////////////////////////////////////

__declspec(dllexport) LPXLOPER WINAPI Local_ReplicConvAdjust_Create(
	LPXLOPER XL_Payoff_ReplicMode,
	LPXLOPER XL_Payoff_StepOrReplicPrecision,
	LPXLOPER XL_Payoff_StopMode,
	LPXLOPER XL_Payoff_StopThreshold,
	LPXLOPER XL_Sensi_ReplicMode,
	LPXLOPER XL_Sensi_StepOrReplicPrecision,
	LPXLOPER XL_Sensi_StopMode,
	LPXLOPER XL_Sensi_StopThreshold,
	LPXLOPER XL_UsedModel,
	LPXLOPER XL_StrikeMinReplic,
	LPXLOPER XL_StrikeMaxReplic)
{
	ADD_LOG("Local_ReplicConvAdjust_Create");
	bool PersistentInXL = true;
	return Local_ReplicConvAdjust_CreateCommon(
			XL_Payoff_ReplicMode,
			XL_Payoff_StepOrReplicPrecision,
			XL_Payoff_StopMode,
			XL_Payoff_StopThreshold,
			XL_Sensi_ReplicMode,
			XL_Sensi_StepOrReplicPrecision,
			XL_Sensi_StopMode,
			XL_Sensi_StopThreshold,
			XL_UsedModel,
			XL_StrikeMinReplic,
			XL_StrikeMaxReplic,
			PersistentInXL);
}



//////////////////////////////////////////////////////
/// Addin to create a Replication Convexity Adjust Manager
/// Version for VBA
//////////////////////////////////////////////////////

__declspec(dllexport) LPXLOPER WINAPI Local_PXL_ReplicConvAdjust_Create(
    LPXLOPER XL_Payoff_ReplicMode,
	LPXLOPER XL_Payoff_ReplicPrecision,
	LPXLOPER XL_Payoff_StopMode,
	LPXLOPER XL_Payoff_StopThreshold,
	LPXLOPER XL_Sensi_ReplicMode,
	LPXLOPER XL_Sensi_ReplicPrecision,
	LPXLOPER XL_Sensi_StopMode,
	LPXLOPER XL_Sensi_StopThreshold,
	LPXLOPER XL_UsedModel,
	LPXLOPER XL_StrikeMinReplic,
	LPXLOPER XL_StrikeMaxReplic)
{
	ADD_LOG("Local_PXL_ReplicConvAdjust_Create");
	bool PersistentInXL = false;
	return Local_ReplicConvAdjust_CreateCommon(
			XL_Payoff_ReplicMode,
			XL_Payoff_ReplicPrecision,
			XL_Payoff_StopMode,
			XL_Payoff_StopThreshold,
			XL_Sensi_ReplicMode,
			XL_Sensi_ReplicPrecision,
			XL_Sensi_StopMode,
			XL_Sensi_StopThreshold,
			XL_UsedModel,
			XL_StrikeMinReplic,
			XL_StrikeMaxReplic,
			PersistentInXL);
}


class BSConvAdjustRepFunc : public ARMResultLong2LongFunc
{
public:
	BSConvAdjustRepFunc(long UsedModelId,
						long SwoptVolId,
						const VECTOR<double>& Stddev,
						int NbPtsForRepliq,
						long MR,
						bool FullRepliq,
						double upperProba,
						double lowerProba):
						C_UsedModelId(UsedModelId), C_SwoptVolId(SwoptVolId),
						C_Stddev(Stddev), C_NbPtsForRepliq(NbPtsForRepliq),
						C_MR(MR), C_FullRepliq(FullRepliq), C_UpperProba(upperProba), C_LowerProba(lowerProba)
	{};

	long operator()( ARM_result& result, long objId )
	{
		return ARMLOCAL_BSConvAdjustRep_Create(
			C_UsedModelId,
			C_SwoptVolId,
			C_Stddev,
			C_NbPtsForRepliq,
			C_MR,
			C_FullRepliq,
			C_UpperProba,
			C_LowerProba,
			result, 
			objId);
	};

private:
	long C_UsedModelId;
	long C_SwoptVolId;
	VECTOR<double> C_Stddev;
	int C_NbPtsForRepliq;
	long C_MR;
	bool C_FullRepliq;
	double C_UpperProba;
	double C_LowerProba;
};

_declspec(dllexport) LPXLOPER WINAPI Local_BSConvAdjustRep_CreateCommon(
	LPXLOPER XL_UsedModelId,
	LPXLOPER XL_SwoptVolId,
	LPXLOPER XL_MR,
	LPXLOPER XL_Stddev,
	LPXLOPER XL_NbPtsForRepliq,
	LPXLOPER XL_FullRepliq, 
	LPXLOPER XL_UpperProba,
	LPXLOPER XL_LowerProba,
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
	
	// error
	static int error;
	static char* reason = "";

	CCString C_UsedModel, C_SwoptVol, C_MR, C_FullRepliq, defRep("N");
	long C_UsedModelId, C_SwoptVolId, C_MRId;
	VECTOR<double> C_Stddev, defVec(0);
	double C_NbPtsForRepliq, defNb(8), C_UpperProba, C_LowerProba, defUProba(0.999), defLProba(0.001);
	bool fullRepliq;

	XL_readStrCell( XL_UsedModelId,		C_UsedModel, " ARM_ERR: Model Id: string expected",				C_result);
	XL_readStrCell( XL_SwoptVolId, C_SwoptVol, "ARM_ERR : Swopt Vol Id : string expected", C_result);
	XL_readStrCell( XL_MR, C_MR, "ARM_ERR : MR Id : string expected", C_result);
	XL_readNumVectorWD(XL_Stddev, C_Stddev, defVec, "ARM_ERR : Stddev : array expected", C_result);
	XL_readNumCellWD(XL_NbPtsForRepliq, C_NbPtsForRepliq, defNb, "ARM_ERR : Nb Pts for Repliq : numeric expected", C_result);
	XL_readStrCellWD(XL_FullRepliq, C_FullRepliq, defRep, "ARM_ERR : Full Replication : string expected", C_result);
	XL_readNumCellWD(XL_UpperProba, C_UpperProba, defUProba, "ARM_ERR : Upper Proba : numeric expected", C_result);
	XL_readNumCellWD(XL_LowerProba, C_LowerProba, defLProba, "ARM_ERR : Upper Proba : numeric expected", C_result);

	if (C_UsedModel == "NULL")
		C_UsedModelId = ARM_NULL_OBJECT;
	else
		C_UsedModelId = LocalGetNumObjectId (C_UsedModel);

	if (C_SwoptVol == "NULL")
		C_SwoptVolId = ARM_NULL_OBJECT;
	else
		C_SwoptVolId = LocalGetNumObjectId (C_SwoptVol);

	if (C_MR == "NULL")
		C_MRId = ARM_NULL_OBJECT;
	else
		C_MRId = LocalGetNumObjectId (C_MR);

	C_FullRepliq.toUpper();
	fullRepliq = C_FullRepliq == "Y" ? true : false;

	BSConvAdjustRepFunc ourFunc(C_UsedModelId, C_SwoptVolId, C_Stddev, C_NbPtsForRepliq, C_MRId, fullRepliq, C_UpperProba, C_LowerProba);

	/// call the general function
	fillXL_Result( LOCAL_REPLICCONVADJUST_CLASS, ourFunc, C_result, XL_result, PersistentInXL );

	/// return the result as an LPXLOPER
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_ReplicConvAdjust_CreateCommon" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;

}

_declspec(dllexport) LPXLOPER WINAPI Local_BSConvAdjustRep_Create(
	LPXLOPER XL_UsedModelId,
	LPXLOPER XL_SwoptVolId,
	LPXLOPER XL_MR,
	LPXLOPER XL_Stddev,
	LPXLOPER XL_NbPtsForRepliq,
	LPXLOPER XL_FullRepliq, 
	LPXLOPER XL_UpperProba,
	LPXLOPER XL_LowerProba)
{
	bool PersistentInXL = true;

	return Local_BSConvAdjustRep_CreateCommon(
		XL_UsedModelId,
		XL_SwoptVolId,
		XL_MR,
		XL_Stddev,
		XL_NbPtsForRepliq,
		XL_FullRepliq,
		XL_UpperProba,
		XL_LowerProba,
		PersistentInXL);
}

_declspec(dllexport) LPXLOPER WINAPI Local_PXL_BSConvAdjustRep_Create(
	LPXLOPER XL_UsedModelId,
	LPXLOPER XL_SwoptVolId,
	LPXLOPER XL_MR,
	LPXLOPER XL_Stddev,
	LPXLOPER XL_NbPtsForRepliq,
	LPXLOPER XL_FullRepliq, 
	LPXLOPER XL_UpperProba,
	LPXLOPER XL_LowerProba)
{
	bool PersistentInXL = false;

	return Local_BSConvAdjustRep_CreateCommon(
		XL_UsedModelId,
		XL_SwoptVolId,
		XL_MR,
		XL_Stddev,
		XL_NbPtsForRepliq,
		XL_FullRepliq,
		XL_UpperProba,
		XL_LowerProba,
		PersistentInXL);
}
/////////////////////////////////////////
/// Functor to create a map convexity 
/// adjustment manager
/////////////////////////////////////////
class MapConvAdjFunc : public ARMResultLong2LongFunc
{
public:
	MapConvAdjFunc(long LiborArrearAdjId,
					  long NaturalCMSAdjId,
					  long PaymentLagAdjId)
					 : C_LiborArrearAdjId(LiborArrearAdjId),
					   C_NaturalCMSAdjId(NaturalCMSAdjId),
					   C_PaymentLagAdjId(PaymentLagAdjId)
	{};

	long operator()( ARM_result& result, long objId )
	{
		return ARMLOCAL_MapConvAdjust_Create(
			C_LiborArrearAdjId,
			C_NaturalCMSAdjId,
			C_PaymentLagAdjId,
			result, 
			objId);
	};

private:
	long C_LiborArrearAdjId;
	long C_NaturalCMSAdjId;
	long C_PaymentLagAdjId;

};


/////////////////////////////////////////
/// MapationConvexity Adjust Manager 
/// create common function
/////////////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_MapConvAdjust_CreateCommon(
	LPXLOPER XL_LiborArrearAdj,
	LPXLOPER XL_NaturalCMSAdj,
	LPXLOPER XL_PaymentLagAdj,
	bool PersistentInXL )
{
	ADD_LOG("Local_MapConvAdjust_CreateCommon");
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

	CCString C_LiborArrearAdj;
	long C_LiborArrearAdjId;

	CCString C_NaturalCMSAdj;
	long C_NaturalCMSAdjId;

	CCString C_PaymentLagAdj;
	long C_PaymentLagAdjId;

	XL_readStrCell(	 XL_LiborArrearAdj,		C_LiborArrearAdj,	" ARM_ERR: Cnvx Adj Id: string expected",				C_result);
	XL_readStrCell(	 XL_NaturalCMSAdj,		C_NaturalCMSAdj,	" ARM_ERR: Cnvx Adj Id: string expected",				C_result);
	XL_readStrCell(	 XL_PaymentLagAdj,		C_PaymentLagAdj,	" ARM_ERR: Cnvx Adj Id: string expected",				C_result);

	C_LiborArrearAdjId = LocalGetNumObjectId (C_LiborArrearAdj);
	C_NaturalCMSAdjId = LocalGetNumObjectId (C_NaturalCMSAdj);
	C_PaymentLagAdjId = LocalGetNumObjectId (C_PaymentLagAdj);

	/// use the concept of Functor to transfer the knowledge of
	/// a function with a context
	MapConvAdjFunc ourFunc(
		C_LiborArrearAdjId,
		C_NaturalCMSAdjId,
		C_PaymentLagAdjId);
	
	/// call the general function
	fillXL_Result( LOCAL_MAPCONVADJUST_CLASS, ourFunc, C_result, XL_result, PersistentInXL );

	/// return the result as an LPXLOPER
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_MapConvAdjust_CreateCommon" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


//////////////////////////////////////////////////////
/// Addin to create a Mapping Convexity Adjustment 
///	Manager
//////////////////////////////////////////////////////

__declspec(dllexport) LPXLOPER WINAPI Local_MapConvAdjust_Create(
	LPXLOPER XL_LiborArrearAdj,
	LPXLOPER XL_NaturalCMSAdj,
	LPXLOPER XL_PaymentLagAdj)
{
	ADD_LOG("Local_MapConvAdjust_Create");
	bool PersistentInXL = true;
	return Local_MapConvAdjust_CreateCommon(
			XL_LiborArrearAdj,
			XL_NaturalCMSAdj,
			XL_PaymentLagAdj,
			PersistentInXL);
}



//////////////////////////////////////////////////////
/// Addin to create a Mapping Convexity Adjust Manager
/// Version for VBA
//////////////////////////////////////////////////////

__declspec(dllexport) LPXLOPER WINAPI Local_PXL_MapConvAdjust_Create(
    LPXLOPER XL_LiborArrearAdj,
	LPXLOPER XL_NaturalCMSAdj,
	LPXLOPER XL_PaymentLagAdj)
{
	ADD_LOG("Local_PXL_MapConvAdjust_Create");
	bool PersistentInXL = false;
	return Local_MapConvAdjust_CreateCommon(
			XL_LiborArrearAdj,
			XL_NaturalCMSAdj,
			XL_PaymentLagAdj,
			PersistentInXL);
}


/////////////////////////////////////////
/// Functor to create a Replic Model
/////////////////////////////////////////
class ReplicModelFunc : public ARMResultLong2LongFunc
{
public:
	ReplicModelFunc(int ReplicMode,
					double ReplicPrecision,
					int StopMode,
					double StopThreshold,
                    int SensiReplicMode,
					double SensiReplicPrecision,
					int SensiStopMode,
					double SensiStopThreshold,
					long UsedModel,
                    int SUMMITFormulaeUsed,
					double StrikeMinReplic,
					double StrikeMaxReplic)
					 : C_ReplicMode(ReplicMode),
					   C_ReplicPrecision(ReplicPrecision),
					   C_StopMode(StopMode),
					   C_StopThreshold(StopThreshold),
                       C_SensiReplicMode(SensiReplicMode),
					   C_SensiReplicPrecision(SensiReplicPrecision),
					   C_SensiStopMode(SensiStopMode),
					   C_SensiStopThreshold(SensiStopThreshold),
					   C_UsedModel(UsedModel),
                       C_SUMMITFormulaeUsed(SUMMITFormulaeUsed),
					   C_StrikeMinReplic(StrikeMinReplic),
					   C_StrikeMaxReplic(StrikeMaxReplic)
	{};

	long operator()( ARM_result& result, long objId )
	{
		return ARMLOCAL_ReplicModel_Create(
			C_ReplicMode,
			C_ReplicPrecision,
			C_StopMode,
			C_StopThreshold,
            C_SensiReplicMode,
			C_SensiReplicPrecision,
			C_SensiStopMode,
			C_SensiStopThreshold,
			C_UsedModel,
            C_SUMMITFormulaeUsed,
			C_StrikeMinReplic,
			C_StrikeMaxReplic,
			result, 
			objId);
	};

private:
	int C_ReplicMode;
	double C_ReplicPrecision;
	int C_StopMode;
	double C_StopThreshold;
    int C_SensiReplicMode;
	double C_SensiReplicPrecision;
	int C_SensiStopMode;
	double C_SensiStopThreshold;
	long C_UsedModel;
    int C_SUMMITFormulaeUsed;
	double C_StrikeMinReplic;
	double C_StrikeMaxReplic;

};


/////////////////////////////////////////
/// Replic Model create common function
/////////////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_ReplicModel_CreateCommon(
	LPXLOPER XL_ReplicMode,
	LPXLOPER XL_ReplicPrecision,
	LPXLOPER XL_StopMode,
	LPXLOPER XL_StopThreshold,
    LPXLOPER XL_SensiReplicMode,
	LPXLOPER XL_SensiReplicPrecision,
	LPXLOPER XL_SensiStopMode,
	LPXLOPER XL_SensiStopThreshold,
    LPXLOPER XL_UsedModel,
    LPXLOPER XL_SUMMITFormulaeUsed,
	LPXLOPER XL_StrikeMinReplic,
	LPXLOPER XL_StrikeMaxReplic,
	bool PersistentInXL )
{
	ADD_LOG("Local_ReplicModel_CreateCommon");
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

	CCString C_ReplicModeStr;
	int C_ReplicMode;
	double C_ReplicPrecision;
	CCString C_StopModeStr;
	int C_StopMode;
	double C_StopThreshold;
    CCString C_SensiReplicModeStr;
	int C_SensiReplicMode;
	double C_SensiReplicPrecision;
	CCString C_SensiStopModeStr;
	int C_SensiStopMode;
	double C_SensiStopThreshold;
    CCString C_SUMMITFormulaeUsedStr;
    int C_SUMMITFormulaeUsed;
	CCString C_UsedModel;
	long C_UsedModelId;
	double C_StrikeMinReplic;
	double C_StrikeMaxReplic;
	double C_StrikeMinReplic_Default = (double) K_STRIKEMIN_REPLIC;
	double C_StrikeMaxReplic_Default = (double) K_STRIKEMAX_REPLIC;


	XL_readStrCell(	 XL_ReplicMode,		C_ReplicModeStr,	" ARM_ERR: replicMode: string expected",				C_result);
	XL_readNumCell(	 XL_ReplicPrecision,		C_ReplicPrecision,	" ARM_ERR: replicPrecision: numeric expected",				C_result);
	XL_readStrCell(	 XL_StopMode,		C_StopModeStr,	" ARM_ERR: stopMode: string expected",				C_result);
	XL_readNumCell(	 XL_StopThreshold,		C_StopThreshold,	" ARM_ERR: stopThreshold: numeric expected",				C_result);

    XL_readStrCell(	 XL_SensiReplicMode,		C_SensiReplicModeStr,	" ARM_ERR: replicMode: string expected",				C_result);
	XL_readNumCell(	 XL_SensiReplicPrecision,		C_SensiReplicPrecision,	" ARM_ERR: replicPrecision: numeric expected",				C_result);
	XL_readStrCell(	 XL_SensiStopMode,		C_SensiStopModeStr,	" ARM_ERR: stopMode: string expected",				C_result);
	XL_readNumCell(	 XL_SensiStopThreshold,		C_SensiStopThreshold,	" ARM_ERR: stopThreshold: numeric expected",				C_result);

	XL_readStrCell(	 XL_UsedModel,		C_UsedModel, " ARM_ERR: Model Id: string expected",				C_result);

    XL_readStrCellWD(	 XL_SUMMITFormulaeUsed,		C_SUMMITFormulaeUsedStr, "N",	" ARM_ERR: SUMMITFormulaeUsedStr: string expected",				C_result);

	XL_readNumCellWD(	 XL_StrikeMinReplic,		C_StrikeMinReplic, C_StrikeMinReplic_Default,	" ARM_ERR: StrikeMinReplic: numeric expected",				C_result);
	XL_readNumCellWD(	 XL_StrikeMaxReplic,		C_StrikeMaxReplic, C_StrikeMaxReplic_Default,	" ARM_ERR: StrikeMinReplic: numeric expected",				C_result);

	if( ((C_ReplicMode	= ARM_ConvReplicMode( C_ReplicModeStr, C_result)) == ARM_DEFAULT_ERR)
		|| ((C_StopMode	= ARM_ConvStopMode( C_StopModeStr, C_result)) == ARM_DEFAULT_ERR) 
        || ((C_SensiReplicMode	= ARM_ConvReplicMode( C_SensiReplicModeStr, C_result)) == ARM_DEFAULT_ERR)
		|| ((C_SensiStopMode	= ARM_ConvStopMode( C_SensiStopModeStr, C_result)) == ARM_DEFAULT_ERR))
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	if (C_UsedModel == "NULL")
		C_UsedModelId = ARM_NULL_OBJECT;
	else
		C_UsedModelId = LocalGetNumObjectId (C_UsedModel);

    if ( (C_SUMMITFormulaeUsed = ARM_ConvYesOrNo (C_SUMMITFormulaeUsedStr, C_result)) == ARM_DEFAULT_ERR)
	{
	   ARM_ARG_ERR();

	   return (LPXLOPER)&XL_result;
	}

	/// use the concept of Functor to transfer the knowledge of
	/// a function with a context
	ReplicModelFunc ourFunc(
		C_ReplicMode,
		C_ReplicPrecision,
		C_StopMode,
		C_StopThreshold,
        C_SensiReplicMode,
		C_SensiReplicPrecision,
		C_SensiStopMode,
		C_SensiStopThreshold,
		C_UsedModelId,
        C_SUMMITFormulaeUsed,
		C_StrikeMinReplic,
		C_StrikeMaxReplic);
	
	/// call the general function
	fillXL_Result( LOCAL_REPLICMODEL_CLASS, ourFunc, C_result, XL_result, PersistentInXL );

	/// return the result as an LPXLOPER
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_ReplicModel_CreateCommon" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


//////////////////////////////////////////////////////
/// Addin to create a Replic Model
//////////////////////////////////////////////////////

__declspec(dllexport) LPXLOPER WINAPI Local_ReplicModel_Create(
	LPXLOPER XL_ReplicMode,
	LPXLOPER XL_StepOrReplicPrecision,
	LPXLOPER XL_StopMode,
	LPXLOPER XL_StopThreshold,
    LPXLOPER XL_SensiReplicMode,
	LPXLOPER XL_SensiStepOrReplicPrecision,
	LPXLOPER XL_SensiStopMode,
	LPXLOPER XL_SensiStopThreshold,
	LPXLOPER XL_UsedModel,
    LPXLOPER XL_SUMMITFormulaeUsed,
	LPXLOPER XL_StrikeMinReplic,
	LPXLOPER XL_StrikeMaxReplic)
{
	ADD_LOG("Local_ReplicModel_Create");
	bool PersistentInXL = true;
	return Local_ReplicModel_CreateCommon(
			XL_ReplicMode,
			XL_StepOrReplicPrecision,
			XL_StopMode,
			XL_StopThreshold,
            XL_SensiReplicMode,
			XL_SensiStepOrReplicPrecision,
			XL_SensiStopMode,
			XL_SensiStopThreshold,
			XL_UsedModel,
            XL_SUMMITFormulaeUsed,
			XL_StrikeMinReplic,
			XL_StrikeMaxReplic,
			PersistentInXL);
}



//////////////////////////////////////////////////////
/// Addin to create a Replic Model
/// Version for VBA
//////////////////////////////////////////////////////

__declspec(dllexport) LPXLOPER WINAPI Local_PXL_ReplicModel_Create(
    LPXLOPER XL_ReplicMode,
	LPXLOPER XL_StepOrReplicPrecision,
	LPXLOPER XL_StopMode,
	LPXLOPER XL_StopThreshold,
    LPXLOPER XL_SensiReplicMode,
	LPXLOPER XL_SensiStepOrReplicPrecision,
	LPXLOPER XL_SensiStopMode,
	LPXLOPER XL_SensiStopThreshold,
	LPXLOPER XL_UsedModel,
    LPXLOPER XL_SUMMITFormulaeUsed,
	LPXLOPER XL_StrikeMinReplic,
	LPXLOPER XL_StrikeMaxReplic)
{
	ADD_LOG("Local_PXL_ReplicModel_Create");
	bool PersistentInXL = false;
	return Local_ReplicModel_CreateCommon(
			XL_ReplicMode,
			XL_StepOrReplicPrecision,
			XL_StopMode,
			XL_StopThreshold,
            XL_SensiReplicMode,
			XL_SensiStepOrReplicPrecision,
			XL_SensiStopMode,
			XL_SensiStopThreshold,
			XL_UsedModel,
            XL_SUMMITFormulaeUsed,
			XL_StrikeMinReplic,
			XL_StrikeMaxReplic,
			PersistentInXL);
}

/////////////////////////////////////////
/// Set the replic debug mode create common function
/////////////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_SetReplicDebugMode(
	LPXLOPER XL_DebugMode)
{
	ADD_LOG("Local_SetReplicDebugMode");
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

	double C_DebugMode;

	
	XL_readNumCell(	 XL_DebugMode,		C_DebugMode,	" ARM_ERR: replicPrecision: numeric expected",				C_result);

	/// use the concept of Functor to transfer the knowledge of
	/// a function with a context
    long retCode;
	retCode = ARMLOCAL_SetReplicDebugMode(C_DebugMode, C_result);

    if( retCode == ARM_OK)
	{
		FreeCurCellErr ();
		XL_result.xltype = xltypeStr;
		XL_result.val.str = XL_StrC2StrPascal (C_result.getString());
		XL_result.xltype |= xlbitDLLFree;
	}

	/// return the result as an LPXLOPER
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_SetReplicDebugMode" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


/////////////////////////////////////////
/// Functor to create a TriBsModel
/////////////////////////////////////////
class TriBSModelFunc : public ARMResultLong2LongFunc
{
public:
	TriBSModelFunc(long BS1,
                   long BS2,
 				   long DiscountBS,
				   long Fx1DiscVol,
                   long Fx2DiscVol,
				   long Idx1Idx2Corr,
				   long Idx1DiscIdxCorr,
				   long Idx2DiscIdxCorr,
				   long Idx1FxCorr,
				   long Idx2FxCorr,
				   int QuantoFlag)
					 : C_BS1(BS1),
					   C_BS2(BS2),
					   C_DiscountBS(DiscountBS),
					   C_Fx1DiscVol(Fx1DiscVol),
					   C_Fx2DiscVol(Fx2DiscVol),
					   C_Idx1Idx2Corr(Idx1Idx2Corr),
					   C_Idx1DiscIdxCorr(Idx1DiscIdxCorr),
					   C_Idx2DiscIdxCorr(Idx2DiscIdxCorr),
					   C_Idx1FxCorr(Idx1FxCorr),
					   C_Idx2FxCorr(Idx2FxCorr),
					   C_QuantoFlag(QuantoFlag)
	{};

	long operator()( ARM_result& result, long objId )
	{
		return ARMLOCAL_TriBSModel(
			C_BS1,
			C_BS2,
			C_DiscountBS,
			C_Fx1DiscVol,
			C_Fx2DiscVol,
			C_Idx1Idx2Corr,
			C_Idx1DiscIdxCorr,
			C_Idx2DiscIdxCorr,
            C_Idx1FxCorr,
			C_Idx2FxCorr,
			C_QuantoFlag,
			result, 
			objId);
	};

private:
	
	long C_BS1;
    long C_BS2;
    long C_DiscountBS;
    long C_Fx1DiscVol;
    long C_Fx2DiscVol;
    long C_Idx1Idx2Corr;
    long C_Idx1DiscIdxCorr;
    long C_Idx2DiscIdxCorr;
    long C_Idx1FxCorr;
    long C_Idx2FxCorr;
    int C_QuantoFlag;

};


/////////////////////////////////////////
/// ReplicationConvexity Adjust Manager create common function
/////////////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_ARM_TRIBSCommun (LPXLOPER XL_Model1,
                                                               LPXLOPER XL_Model2,
                                                               LPXLOPER XL_DiscountModel,
                                                               LPXLOPER XL_Fx1DiscVol,
                                                               LPXLOPER XL_Fx2DiscVol,
                                                               LPXLOPER XL_Idx1Idx2Corr,
                                                               LPXLOPER XL_Idx1DiscIdxCorr,
                                                               LPXLOPER XL_Idx2DiscIdxCorr,
                                                               LPXLOPER XL_Idx1FxCorr,
                                                               LPXLOPER XL_Idx2FxCorr,
                                                               LPXLOPER XL_quantoflag,
                                                               bool PersistentInXL )
{
	ADD_LOG("Local_ARM_TRIBSCommun ");
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

	CCString C_Model1;
	CCString C_Model2;
    CCString C_DiscountModel;
    long C_DiscountModelId;
	CCString C_Fx1DiscVol;
	CCString C_Fx2DiscVol;
	CCString C_Idx1Idx2Corr;
	CCString C_Idx1DiscIdxCorr;
    CCString C_Idx2DiscIdxCorr;
	CCString C_Idx1FxCorr;
	CCString C_Idx2FxCorr;
	double C_quantoflag;
	double C_adjflag_default = 1.0;
	
	
	XL_readStrCell(	 XL_Model1,		    C_Model1,	    " ARM_ERR: First BSModel Id: string expected",				C_result);
	XL_readStrCell(	 XL_Model2,		    C_Model2,	    " ARM_ERR: Second BSModel Id: string expected",				C_result);
	XL_readStrCellWD(	 XL_DiscountModel,	C_DiscountModel, "DEFAULT"," ARM_ERR: Discount BSModel Id: string expected",				C_result);
	XL_readStrCell(	 XL_Fx1DiscVol,		C_Fx1DiscVol,	" ARM_ERR: Fx1DiscVol Id: string expected",				C_result);
	XL_readStrCell(	 XL_Fx2DiscVol,		C_Fx2DiscVol,	" ARM_ERR: Fx2DiscVol Id: string expected",				C_result);
	XL_readStrCell(	 XL_Idx1Idx2Corr,	C_Idx1Idx2Corr,	" ARM_ERR: Idx1Idx2Corr Id: string expected",				C_result);
	XL_readStrCell(	 XL_Idx1DiscIdxCorr,C_Idx1DiscIdxCorr,	" ARM_ERR: Idx1DiscIdxCorr Id: string expected",				C_result);
	XL_readStrCell(	 XL_Idx2DiscIdxCorr,C_Idx2DiscIdxCorr,	" ARM_ERR: Idx2DiscIdxCorr Id: string expected",				C_result);
	XL_readStrCell(	 XL_Idx1FxCorr,		C_Idx1FxCorr,	" ARM_ERR: Idx1FxCorr Id: string expected",				C_result);
	XL_readStrCell(	 XL_Idx2FxCorr,		C_Idx2FxCorr,	" ARM_ERR: Idx2FxCorr Id: string expected",				C_result);
	XL_readNumCellWD(XL_quantoflag,C_quantoflag,C_adjflag_default," ARM_ERR: quanto convexity adjustment flag: numeric expected",C_result);

    if (C_DiscountModel == "DEFAULT")
	{
	   C_DiscountModelId = ARM_NULL_OBJECT;
	}
	else
	{
	   C_DiscountModelId = LocalGetNumObjectId (C_DiscountModel);
	}

	/// use the concept of Functor to transfer the knowledge of
	/// a function with a context
	TriBSModelFunc ourFunc(LocalGetNumObjectId (C_Model1),
			               LocalGetNumObjectId (C_Model2),
                           C_DiscountModelId,
                           LocalGetNumObjectId (C_Fx1DiscVol),
                           LocalGetNumObjectId (C_Fx2DiscVol),
                           LocalGetNumObjectId (C_Idx1Idx2Corr),
                           LocalGetNumObjectId (C_Idx1DiscIdxCorr),
                           LocalGetNumObjectId (C_Idx2DiscIdxCorr),
                           LocalGetNumObjectId (C_Idx1FxCorr),
                           LocalGetNumObjectId (C_Idx2FxCorr),
                		    C_quantoflag);
	
	/// call the general function
	fillXL_Result( LOCAL_TRIBSMOD_CLASS, ourFunc, C_result, XL_result, PersistentInXL );

	/// return the result as an LPXLOPER
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_ARM_TRIBSCommun" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


//////////////////////////////////////////////////////
/// Addin to create a Replication Convexity Adjust Manager
//////////////////////////////////////////////////////

__declspec(dllexport) LPXLOPER WINAPI Local_ARM_TRIBS (LPXLOPER XL_Model1,
                                                       LPXLOPER XL_Model2,
                                                       LPXLOPER XL_DiscountModel,
                                                       LPXLOPER XL_Fx1DiscVol,
                                                       LPXLOPER XL_Fx2DiscVol,
                                                       LPXLOPER XL_Idx1Idx2Corr,
                                                       LPXLOPER XL_Idx1DiscIdxCorr,
                                                       LPXLOPER XL_Idx2DiscIdxCorr,
                                                       LPXLOPER XL_Idx1FxCorr,
                                                       LPXLOPER XL_Idx2FxCorr,
                                                       LPXLOPER XL_quantoflag)
{
	ADD_LOG("Local_ARM_TRIBS ");
	bool PersistentInXL = true;
	return Local_ARM_TRIBSCommun(
               XL_Model1,
               XL_Model2,
               XL_DiscountModel,
               XL_Fx1DiscVol,
               XL_Fx2DiscVol,
               XL_Idx1Idx2Corr,
               XL_Idx1DiscIdxCorr,
               XL_Idx2DiscIdxCorr,
               XL_Idx1FxCorr,
               XL_Idx2FxCorr,
               XL_quantoflag,
			   PersistentInXL);
}



//////////////////////////////////////////////////////
/// Addin to create a Replication Convexity Adjust Manager
/// Version for VBA
//////////////////////////////////////////////////////

__declspec(dllexport) LPXLOPER WINAPI Local_PXL_TRIBS( LPXLOPER XL_Model1,
                                                       LPXLOPER XL_Model2,
                                                       LPXLOPER XL_DiscountModel,
                                                       LPXLOPER XL_Fx1DiscVol,
                                                       LPXLOPER XL_Fx2DiscVol,
                                                       LPXLOPER XL_Idx1Idx2Corr,
                                                       LPXLOPER XL_Idx1DiscIdxCorr,
                                                       LPXLOPER XL_Idx2DiscIdxCorr,
                                                       LPXLOPER XL_Idx1FxCorr,
                                                       LPXLOPER XL_Idx2FxCorr,
                                                       LPXLOPER XL_quantoflag)
{
	ADD_LOG("Local_PXL_TRIBS");
	bool PersistentInXL = false;
	return Local_ARM_TRIBSCommun(
			   XL_Model1,
               XL_Model2,
               XL_DiscountModel,
               XL_Fx1DiscVol,
               XL_Fx2DiscVol,
               XL_Idx1Idx2Corr,
               XL_Idx1DiscIdxCorr,
               XL_Idx2DiscIdxCorr,
               XL_Idx1FxCorr,
               XL_Idx2FxCorr,
               XL_quantoflag,
			   PersistentInXL);
}
	


__declspec(dllexport) LPXLOPER WINAPI Local_ARM_SetDiscountPricingMode(LPXLOPER XL_model,
																	   LPXLOPER XL_discPriceMode)
{
	ADD_LOG("Local_ARM_SetDiscountPricingMode");
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

	CCString C_model;
	double C_discPriceMode;
	
	XL_readStrCell(XL_model,C_model," ARM_ERR: model id: object expected",C_result);
	XL_readNumCell(XL_discPriceMode,C_discPriceMode," ARM_ERR: disc pricing method: numeric expected",C_result);

	long retCode;
	long objId;
	CCString prevClass;
	
	CCString curClass = LocalGetStringObjectClass(C_model);
	CCString stringId = GetLastCurCellEnvValue ();
	
	retCode = ARMLOCAL_SetDiscountPricingMode(LocalGetNumObjectId(C_model),
											  (long)C_discPriceMode,
											  C_result);

	if(retCode == ARM_OK)
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
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_ARM_SetDiscountPricingMode" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


__declspec(dllexport) LPXLOPER WINAPI Local_ARM_CALIBRATORSFRM(LPXLOPER XL_sec,
															   LPXLOPER XL_capModel,
															   LPXLOPER XL_swoptModel,
															   LPXLOPER XL_meanRev,
															   LPXLOPER XL_calibParams,
															   LPXLOPER XL_preInitFlag,
															   LPXLOPER XL_initSigmaCurve,
															   LPXLOPER XL_initBetaOrShift,
															   LPXLOPER XL_sigmaPF,
															   LPXLOPER XL_betaPF,
															   LPXLOPER XL_meanrevPF,
															   LPXLOPER XL_voltype,
															   LPXLOPER XL_correl,
                                                               LPXLOPER XL_SecurityParams,
                                                               LPXLOPER XL_tocalswaptatm)
{
	ADD_LOG("Local_ARM_CALIBRATORSFRM");
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();

	// C variable
	CCString C_sec;
	CCString C_capModel;
	CCString C_swoptModel;

	double C_meanRev;

	VECTOR<double> C_calibParams;

	VECTOR<double>C_initSigmaCurve;
	VECTOR<double>C_initSigmaCurve_default;
	VECTOR<double> C_initBetaOrShift;
	VECTOR<double> C_initBetaOrShift_default;

	CCString C_preInitFlag;
	long preInitFlagId;

	CCString C_sigmaPF;
	long sigmaPFId;
	CCString C_betaPF;
	long betaPFId;
	CCString C_meanrevPF;
	long meanrevPFId;

	CCString C_correl;
	long correlId;

	CCString C_voltype;
	long voltypeId;

	// error
	static int error;
	static char* reason = "";

	XL_readStrCell(XL_sec,C_sec," ARM_ERR: security id: object expected",C_result);
	XL_readStrCell(XL_capModel,C_capModel," ARM_ERR: cap model id: object expected",C_result);
	XL_readStrCell(XL_swoptModel,C_swoptModel," ARM_ERR: swopt model id: object expected",C_result);
	XL_readNumCell(XL_meanRev,C_meanRev," ARM_ERR: mean reversion: numeric expected",C_result);
	XL_readNumVector(XL_calibParams,C_calibParams," ARM_ERR: calib params : array of numeric expected",C_result);

	XL_readStrCellWD(XL_preInitFlag,C_preInitFlag,"NO"," ARM_ERR: pre Init Flag: string expected",C_result);
	XL_readNumVectorWD(XL_initSigmaCurve,C_initSigmaCurve,C_initSigmaCurve_default," ARM_ERR: init Sigma Curve : array of numeric expected",C_result);
	XL_readNumVectorWD(XL_initBetaOrShift,C_initBetaOrShift,C_initBetaOrShift_default," ARM_ERR: init Beta or Shift : array of numeric expected",C_result);

	XL_readStrCellWD(XL_voltype,C_voltype,"DIAG"," ARM_ERR: vol type: string expected",C_result);

	XL_readStrCellWD(XL_sigmaPF,C_sigmaPF,"DEFAULT"," ARM_ERR: sigma PF id: object expected",C_result);
	XL_readStrCellWD(XL_betaPF,C_betaPF,"DEFAULT"," ARM_ERR: beta PF id: object expected",C_result);
	XL_readStrCellWD(XL_meanrevPF,C_meanrevPF,"DEFAULT"," ARM_ERR: meanrev PF id: object expected",C_result);
	XL_readStrCellWD(XL_correl,C_correl,"DEFAULT"," ARM_ERR: correl matrix id: object expected",C_result);
    VECTOR<double> C_SecurityParams;
	VECTOR<double> C_SecurityParams_default;
	XL_readNumVectorWD(XL_SecurityParams,C_SecurityParams,C_SecurityParams_default," ARM_ERR: calib Params: array of numeric expected",C_result);

    // C variable
	CCString C_Tocalswaptatm;
    XL_readStrCellWD(XL_tocalswaptatm,C_Tocalswaptatm,"NO"," ARM_ERR: flag to generate swapt Pf ATM: string expected",C_result);

	if ((voltypeId = ARM_ConvShapeType (C_voltype, C_result)) == ARM_DEFAULT_ERR)
	{
	   ARM_ARG_ERR();
	   return (LPXLOPER)&XL_result;
	}

	if ((preInitFlagId = ARM_ConvYesOrNo (C_preInitFlag, C_result)) == ARM_DEFAULT_ERR)
	{
	   ARM_ARG_ERR();
	   return (LPXLOPER)&XL_result;
	}

    if (C_sigmaPF == "DEFAULT")
	{
		sigmaPFId = ARM_NULL_OBJECT;
	}
	else
	{
		sigmaPFId = LocalGetNumObjectId (C_sigmaPF);
	}

	if ( C_betaPF == "DEFAULT" )
	{
		betaPFId = ARM_NULL_OBJECT;
	}
	else
	{
		betaPFId = LocalGetNumObjectId(C_betaPF);
	}

	if ( C_meanrevPF == "DEFAULT" )
	{
	   meanrevPFId = ARM_NULL_OBJECT;
	}
	else
	{
	   meanrevPFId = LocalGetNumObjectId(C_meanrevPF);
	}

	if ( C_correl == "DEFAULT" )
	{
	   correlId = ARM_NULL_OBJECT;
	}
	else
	{
	   correlId = LocalGetNumObjectId(C_correl);
	}

	long retCode;
	long objId;
	CCString prevClass;

	CCString curClass = LOCAL_CALIBRATORSFRM_CLASS;
	CCString stringId = GetLastCurCellEnvValue ();

	if (!stringId)
	{
		retCode = ARMLOCAL_CALIBRATORSFRM (LocalGetNumObjectId (C_sec),
										   LocalGetNumObjectId (C_capModel),
										   LocalGetNumObjectId (C_swoptModel),
										   C_meanRev,
										   C_calibParams,
										   preInitFlagId,
										   C_initSigmaCurve,
										   C_initBetaOrShift,
										   sigmaPFId,
										   betaPFId,
										   meanrevPFId,
										   voltypeId,
										   correlId,
                                           C_SecurityParams,
                                           C_Tocalswaptatm,
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
			
		if(curClass == prevClass)
		{
			retCode = ARMLOCAL_CALIBRATORSFRM (LocalGetNumObjectId (C_sec),
											   LocalGetNumObjectId (C_capModel),
											   LocalGetNumObjectId (C_swoptModel),
											   C_meanRev,
											   C_calibParams,
											   preInitFlagId,
											   C_initSigmaCurve,
											   C_initBetaOrShift,
											   sigmaPFId,
											   betaPFId,
											   meanrevPFId,
											   voltypeId,
											   correlId,
                                               C_SecurityParams,
                                               C_Tocalswaptatm,
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

			retCode = ARMLOCAL_CALIBRATORSFRM (LocalGetNumObjectId (C_sec),
											   LocalGetNumObjectId (C_capModel),
											   LocalGetNumObjectId (C_swoptModel),
											   C_meanRev,
											   C_calibParams,
											   preInitFlagId,
											   C_initSigmaCurve,
											   C_initBetaOrShift,
											   sigmaPFId,
											   betaPFId,
											   meanrevPFId,
											   voltypeId,
											   correlId,
                                               C_SecurityParams,
                                               C_Tocalswaptatm,
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

//	ARM_END();
	
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_ARM_CALIBRATORSFRM" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


__declspec(dllexport) LPXLOPER WINAPI Local_PXL_ARM_CALIBRATORSFRM(LPXLOPER XL_sec,
																   LPXLOPER XL_capModel,
																   LPXLOPER XL_swoptModel,
																   LPXLOPER XL_meanRev,
																   LPXLOPER XL_calibParams,
																   LPXLOPER XL_preInitFlag,
																   LPXLOPER XL_initSigmaCurve,
																   LPXLOPER XL_initBetaOrShift,
																   LPXLOPER XL_sigmaPF,
																   LPXLOPER XL_betaPF,
																   LPXLOPER XL_meanrevPF,
																   LPXLOPER XL_voltype,
																   LPXLOPER XL_correl,
                                                                   LPXLOPER XL_SecurityParams,
                                                                   LPXLOPER XL_tocalswaptatm)
{
	ADD_LOG("Local_PXL_ARM_CALIBRATORSFRM");
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();

	// C variable
	CCString C_sec;
	CCString C_capModel;
	CCString C_swoptModel;

	double C_meanRev;

	VECTOR<double> C_calibParams;

	VECTOR<double>C_initSigmaCurve;
	VECTOR<double>C_initSigmaCurve_default;
	VECTOR<double> C_initBetaOrShift;
	VECTOR<double> C_initBetaOrShift_default;

	CCString C_preInitFlag;
	long preInitFlagId;

	CCString C_sigmaPF;
	long sigmaPFId;
	CCString C_betaPF;
	long betaPFId;
	CCString C_meanrevPF;
	long meanrevPFId;

	CCString C_correl;
	long correlId;

	CCString C_voltype;
	long voltypeId;

	// error
	static int error;
	static char* reason = "";

	XL_readStrCell(XL_sec,C_sec," ARM_ERR: security id: object expected",C_result);
	XL_readStrCell(XL_capModel,C_capModel," ARM_ERR: cap model id: object expected",C_result);
	XL_readStrCell(XL_swoptModel,C_swoptModel," ARM_ERR: swopt model id: object expected",C_result);
	XL_readNumCell(XL_meanRev,C_meanRev," ARM_ERR: mean reversion: numeric expected",C_result);
	XL_readNumVector(XL_calibParams,C_calibParams," ARM_ERR: calib params : array of numeric expected",C_result);

	XL_readStrCellWD(XL_preInitFlag,C_preInitFlag,"NO"," ARM_ERR: pre Init Flag: string expected",C_result);
	XL_readNumVectorWD(XL_initSigmaCurve,C_initSigmaCurve,C_initSigmaCurve_default," ARM_ERR: init Sigma Curve : array of numeric expected",C_result);
	XL_readNumVectorWD(XL_initBetaOrShift,C_initBetaOrShift,C_initBetaOrShift_default," ARM_ERR: init Beta or Shift : array of numeric expected",C_result);

	XL_readStrCellWD(XL_voltype,C_voltype,"DIAG"," ARM_ERR: vol type: string expected",C_result);

	XL_readStrCellWD(XL_sigmaPF,C_sigmaPF,"DEFAULT"," ARM_ERR: sigma PF id: object expected",C_result);
	XL_readStrCellWD(XL_betaPF,C_betaPF,"DEFAULT"," ARM_ERR: beta PF id: object expected",C_result);
	XL_readStrCellWD(XL_meanrevPF,C_meanrevPF,"DEFAULT"," ARM_ERR: meanrev PF id: object expected",C_result);
	XL_readStrCellWD(XL_correl,C_correl,"DEFAULT"," ARM_ERR: correl matrix id: object expected",C_result);

    VECTOR<double> C_SecurityParams;
	VECTOR<double> C_SecurityParams_default;
	XL_readNumVectorWD(XL_SecurityParams,C_SecurityParams,C_SecurityParams_default," ARM_ERR: calib Params: array of numeric expected",C_result);

    // C variable
	CCString C_Tocalswaptatm;
    XL_readStrCellWD(XL_tocalswaptatm,C_Tocalswaptatm,"NO"," ARM_ERR: flag to generate swapt Pf ATM: string expected",C_result);


	if ((voltypeId = ARM_ConvShapeType (C_voltype, C_result)) == ARM_DEFAULT_ERR)
	{
	   ARM_ARG_ERR();
	   return (LPXLOPER)&XL_result;
	}

	if ((preInitFlagId = ARM_ConvYesOrNo (C_preInitFlag, C_result)) == ARM_DEFAULT_ERR)
	{
	   ARM_ARG_ERR();
	   return (LPXLOPER)&XL_result;
	}

    if (C_sigmaPF == "DEFAULT")
	{
		sigmaPFId = ARM_NULL_OBJECT;
	}
	else
	{
		sigmaPFId = LocalGetNumObjectId (C_sigmaPF);
	}

	if ( C_betaPF == "DEFAULT" )
	{
		betaPFId = ARM_NULL_OBJECT;
	}
	else
	{
		betaPFId = LocalGetNumObjectId(C_betaPF);
	}

	if ( C_meanrevPF == "DEFAULT" )
	{
	   meanrevPFId = ARM_NULL_OBJECT;
	}
	else
	{
	   meanrevPFId = LocalGetNumObjectId(C_meanrevPF);
	}

	if ( C_correl == "DEFAULT" )
	{
	   correlId = ARM_NULL_OBJECT;
	}
	else
	{
	   correlId = LocalGetNumObjectId(C_correl);
	}

	long retCode;
	long objId;

	CCString curClass = LOCAL_CALIBRATORSFRM_CLASS;
	CCString stringId;

	retCode = ARMLOCAL_CALIBRATORSFRM (LocalGetNumObjectId (C_sec),
									   LocalGetNumObjectId (C_capModel),
									   LocalGetNumObjectId (C_swoptModel),
									   C_meanRev,
									   C_calibParams,
									   preInitFlagId,
									   C_initSigmaCurve,
									   C_initBetaOrShift,
									   sigmaPFId,
									   betaPFId,
									   meanrevPFId,
									   voltypeId,
									   correlId,
                                       C_SecurityParams,
                                       C_Tocalswaptatm,
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_PXL_ARM_CALIBRATORSFRM" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


/////////////////////////////////////////
/// Functor to create a TriBsDual
/////////////////////////////////////////
class TriBSDualModelFunc : public ARMResultLong2LongFunc
{
public:
	TriBSDualModelFunc(long BS1,
                       long BS2,
 				       long DiscountBS,
				       long Fx1DiscVol,
                       long Fx2DiscVol,
				       long Idx1Idx2Corr,
				       long Idx1DiscIdxCorr,
				       long Idx2DiscIdxCorr,
				       long Idx1FxCorr,
				       long Idx2FxCorr,
				       int QuantoFlag,
                       double CorrelForAdj,
                       int WhithSlopeFlag)
					 : C_BS1(BS1),
					   C_BS2(BS2),
					   C_DiscountBS(DiscountBS),
					   C_Fx1DiscVol(Fx1DiscVol),
					   C_Fx2DiscVol(Fx2DiscVol),
					   C_Idx1Idx2Corr(Idx1Idx2Corr),
					   C_Idx1DiscIdxCorr(Idx1DiscIdxCorr),
					   C_Idx2DiscIdxCorr(Idx2DiscIdxCorr),
					   C_Idx1FxCorr(Idx1FxCorr),
					   C_Idx2FxCorr(Idx2FxCorr),
					   C_QuantoFlag(QuantoFlag),
                       C_CorrelForAdj(CorrelForAdj),
                       C_WhithSlopeFlag(WhithSlopeFlag)

	{};

	long operator()( ARM_result& result, long objId )
	{
		return ARMLOCAL_TriBSDualModel(
			C_BS1,
			C_BS2,
			C_DiscountBS,
			C_Fx1DiscVol,
			C_Fx2DiscVol,
			C_Idx1Idx2Corr,
			C_Idx1DiscIdxCorr,
			C_Idx2DiscIdxCorr,
            C_Idx1FxCorr,
			C_Idx2FxCorr,
			C_QuantoFlag,
            C_CorrelForAdj,
            C_WhithSlopeFlag,
			result, 
			objId);
	};

private:
	
	long C_BS1;
    long C_BS2;
    long C_DiscountBS;
    long C_Fx1DiscVol;
    long C_Fx2DiscVol;
    long C_Idx1Idx2Corr;
    long C_Idx1DiscIdxCorr;
    long C_Idx2DiscIdxCorr;
    long C_Idx1FxCorr;
    long C_Idx2FxCorr;
    int C_QuantoFlag;
    double C_CorrelForAdj;
    int C_WhithSlopeFlag;

};


/////////////////////////////////////////
/// TriBsDual Model create common function
/////////////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_ARM_TRIBSDualCommun (LPXLOPER XL_Model1,
                                                               LPXLOPER XL_Model2,
                                                               LPXLOPER XL_DiscountModel,
                                                               LPXLOPER XL_Fx1DiscVol,
                                                               LPXLOPER XL_Fx2DiscVol,
                                                               LPXLOPER XL_Idx1Idx2Corr,
                                                               LPXLOPER XL_Idx1DiscIdxCorr,
                                                               LPXLOPER XL_Idx2DiscIdxCorr,
                                                               LPXLOPER XL_Idx1FxCorr,
                                                               LPXLOPER XL_Idx2FxCorr,
                                                               LPXLOPER XL_quantoflag,
                                                               LPXLOPER XL_correlforadj,
                                                               LPXLOPER XL_withslopeflag,
                                                               bool PersistentInXL )
{
	ADD_LOG("Local_ARM_TRIBSDualCommun ");
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

	CCString C_Model1;
	CCString C_Model2;
    CCString C_DiscountModel;
    long C_DiscountModelId;
	CCString C_Fx1DiscVol;
	CCString C_Fx2DiscVol;
	CCString C_Idx1Idx2Corr;
	CCString C_Idx1DiscIdxCorr;
    CCString C_Idx2DiscIdxCorr;
	CCString C_Idx1FxCorr;
	CCString C_Idx2FxCorr;
	double C_quantoflag;
	double C_adjflag_default = 1.0;
	double C_correlforadj;
    double C_withslopeflag;
	double C_withslopeflag_default = 1.0;
	
	XL_readStrCell(	 XL_Model1,		    C_Model1,	    " ARM_ERR: First BSModel Id: string expected",				C_result);
	XL_readStrCell(	 XL_Model2,		    C_Model2,	    " ARM_ERR: Second BSModel Id: string expected",				C_result);
	XL_readStrCellWD(	 XL_DiscountModel,	C_DiscountModel, "DEFAULT"," ARM_ERR: Discount BSModel Id: string expected",				C_result);
	XL_readStrCell(	 XL_Fx1DiscVol,		C_Fx1DiscVol,	" ARM_ERR: Fx1DiscVol Id: string expected",				C_result);
	XL_readStrCell(	 XL_Fx2DiscVol,		C_Fx2DiscVol,	" ARM_ERR: Fx2DiscVol Id: string expected",				C_result);
	XL_readStrCell(	 XL_Idx1Idx2Corr,	C_Idx1Idx2Corr,	" ARM_ERR: Idx1Idx2Corr Id: string expected",				C_result);
	XL_readStrCell(	 XL_Idx1DiscIdxCorr,C_Idx1DiscIdxCorr,	" ARM_ERR: Idx1DiscIdxCorr Id: string expected",				C_result);
	XL_readStrCell(	 XL_Idx2DiscIdxCorr,C_Idx2DiscIdxCorr,	" ARM_ERR: Idx2DiscIdxCorr Id: string expected",				C_result);
	XL_readStrCell(	 XL_Idx1FxCorr,		C_Idx1FxCorr,	" ARM_ERR: Idx1FxCorr Id: string expected",				C_result);
	XL_readStrCell(	 XL_Idx2FxCorr,		C_Idx2FxCorr,	" ARM_ERR: Idx2FxCorr Id: string expected",				C_result);
	XL_readNumCellWD(XL_quantoflag,C_quantoflag,C_adjflag_default," ARM_ERR: quanto convexity adjustment flag: numeric expected",C_result);
    XL_readNumCell(	 XL_correlforadj,		C_correlforadj,   " ARM_ERR: correlforadj: numeric expected",				C_result);
    XL_readNumCellWD(XL_withslopeflag,C_withslopeflag,C_withslopeflag_default," ARM_ERR: with slope  flag: numeric expected",C_result);

    if (C_DiscountModel == "DEFAULT")
	{
	   C_DiscountModelId = ARM_NULL_OBJECT;
	}
	else
	{
	   C_DiscountModelId = LocalGetNumObjectId (C_DiscountModel);
	}

	/// use the concept of Functor to transfer the knowledge of
	/// a function with a context
	TriBSDualModelFunc ourFunc(LocalGetNumObjectId (C_Model1),
			                   LocalGetNumObjectId (C_Model2),
                               C_DiscountModelId,
                               LocalGetNumObjectId (C_Fx1DiscVol),
                               LocalGetNumObjectId (C_Fx2DiscVol),
                               LocalGetNumObjectId (C_Idx1Idx2Corr),
                               LocalGetNumObjectId (C_Idx1DiscIdxCorr),
                               LocalGetNumObjectId (C_Idx2DiscIdxCorr),
                               LocalGetNumObjectId (C_Idx1FxCorr),
                               LocalGetNumObjectId (C_Idx2FxCorr),
                		        C_quantoflag,
                                C_correlforadj,
                                C_withslopeflag);
	
	/// call the general function
	fillXL_Result( LOCAL_TRIBSMOD_CLASS, ourFunc, C_result, XL_result, PersistentInXL );

	/// return the result as an LPXLOPER
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_ARM_TRIBSDualCommun" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


__declspec(dllexport) LPXLOPER WINAPI Local_ARM_TRIBSDUAL (LPXLOPER XL_Model1,
                                                           LPXLOPER XL_Model2,
                                                           LPXLOPER XL_DiscountModel,
                                                           LPXLOPER XL_Fx1DiscVol,
                                                           LPXLOPER XL_Fx2DiscVol,
                                                           LPXLOPER XL_Idx1Idx2Corr,
                                                           LPXLOPER XL_Idx1DiscIdxCorr,
                                                           LPXLOPER XL_Idx2DiscIdxCorr,
                                                           LPXLOPER XL_Idx1FxCorr,
                                                           LPXLOPER XL_Idx2FxCorr,
                                                           LPXLOPER XL_quantoflag,
                                                           LPXLOPER XL_correlforadj,
                                                           LPXLOPER XL_withslopeflag)
{
	ADD_LOG("Local_ARM_TRIBSDUAL ");
	bool PersistentInXL = true;
	return Local_ARM_TRIBSDualCommun(
               XL_Model1,
               XL_Model2,
               XL_DiscountModel,
               XL_Fx1DiscVol,
               XL_Fx2DiscVol,
               XL_Idx1Idx2Corr,
               XL_Idx1DiscIdxCorr,
               XL_Idx2DiscIdxCorr,
               XL_Idx1FxCorr,
               XL_Idx2FxCorr,
               XL_quantoflag,
               XL_correlforadj,
               XL_withslopeflag,
			   PersistentInXL);
}



//////////////////////////////////////////////////////
/// Addin to create a Replication Convexity Adjust Manager
/// Version for VBA
//////////////////////////////////////////////////////

__declspec(dllexport) LPXLOPER WINAPI Local_PXL_TRIBSDUAL (LPXLOPER XL_Model1,
                                                           LPXLOPER XL_Model2,
                                                           LPXLOPER XL_DiscountModel,
                                                           LPXLOPER XL_Fx1DiscVol,
                                                           LPXLOPER XL_Fx2DiscVol,
                                                           LPXLOPER XL_Idx1Idx2Corr,
                                                           LPXLOPER XL_Idx1DiscIdxCorr,
                                                           LPXLOPER XL_Idx2DiscIdxCorr,
                                                           LPXLOPER XL_Idx1FxCorr,
                                                           LPXLOPER XL_Idx2FxCorr,
                                                           LPXLOPER XL_quantoflag,
                                                           LPXLOPER XL_correlforadj,
                                                           LPXLOPER XL_withslopeflag)

{
	ADD_LOG("Local_PXL_TRIBSDUAL ");
	bool PersistentInXL = false;
	return Local_ARM_TRIBSDualCommun(
			   XL_Model1,
               XL_Model2,
               XL_DiscountModel,
               XL_Fx1DiscVol,
               XL_Fx2DiscVol,
               XL_Idx1Idx2Corr,
               XL_Idx1DiscIdxCorr,
               XL_Idx2DiscIdxCorr,
               XL_Idx1FxCorr,
               XL_Idx2FxCorr,
               XL_quantoflag,
               XL_correlforadj,
               XL_withslopeflag,
			   PersistentInXL);
}

__declspec(dllexport) LPXLOPER WINAPI Local_ARM_SFRMCALIBRATE(LPXLOPER XL_calibratorSFRM,
                                                              LPXLOPER XL_kerneltoGP,
                                                              LPXLOPER XL_MktCapModel,
                                                              LPXLOPER XL_MktSwaptModel,
                                                              LPXLOPER XL_ToCalibrateBeta,
                                                              LPXLOPER XL_ToCalibrateMR)
                                                              
{
	ADD_LOG("Local_ARM_SFRMCALIBRATE");
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
	CCString C_calibratorSFRM;
	XL_readStrCell(XL_calibratorSFRM,C_calibratorSFRM," ARM_ERR: calibrator SFRM id: object expected",C_result);

    // C variable
	CCString C_mktCapModel;
    XL_readStrCellWD(XL_MktCapModel,C_mktCapModel,"DEFAULT"," ARM_ERR: MktCapModel id: object expected",C_result);
    long C_mktCapModelid = C_mktCapModel == "DEFAULT" ? ARM_NULL_OBJECT : LocalGetNumObjectId(C_mktCapModel);

    // C variable
	CCString C_mktSwaptModel;
    XL_readStrCellWD(XL_MktSwaptModel,C_mktSwaptModel,"DEFAULT"," ARM_ERR: MktSwaptModel id: object expected",C_result);
    long C_mktSwapModelid = C_mktSwaptModel == "DEFAULT" ? ARM_NULL_OBJECT : LocalGetNumObjectId(C_mktSwaptModel);

    // C variable
	CCString C_ToCalibrateBeta;
    XL_readStrCellWD(XL_ToCalibrateBeta,C_ToCalibrateBeta,"YES"," ARM_ERR: Flag ToCalibrateBeta: string expected",C_result);

    // C variable
	CCString C_ToCalibrateMR;
    XL_readStrCellWD(XL_ToCalibrateMR,C_ToCalibrateMR,"YES"," ARM_ERR: Flag ToCalibrateMR: string expected",C_result);

    double default_kerneltoGP = 0;
    double C_kerneltoGP;
    XL_readNumCellWD(XL_kerneltoGP,C_kerneltoGP,default_kerneltoGP," ARM_ERR:Flag to choose KernelTo GP: numeric expected",C_result);

	long retCode;
	long objId;
	CCString prevClass;

	CCString curClass = LOCAL_FRMMODEL_CLASS;
	CCString stringId = GetLastCurCellEnvValue ();

	if (!stringId)
	{
		retCode = ARMLOCAL_SFRMCALIBRATE (LocalGetNumObjectId (C_calibratorSFRM),
                                            C_mktCapModelid,
                                            C_mktSwapModelid,
                                            C_ToCalibrateBeta,
                                            C_ToCalibrateMR,
                                            C_kerneltoGP,
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
			
		if(curClass == prevClass)
		{
			retCode = ARMLOCAL_SFRMCALIBRATE (LocalGetNumObjectId (C_calibratorSFRM),
                                            C_mktCapModelid,
                                            C_mktSwapModelid,
                                            C_ToCalibrateBeta,
                                            C_ToCalibrateMR,
											  C_kerneltoGP,
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

			retCode = ARMLOCAL_SFRMCALIBRATE (LocalGetNumObjectId (C_calibratorSFRM),
                                            C_mktCapModelid,
                                            C_mktSwapModelid,
                                            C_ToCalibrateBeta,
                                            C_ToCalibrateMR,
											  C_kerneltoGP,
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

//	ARM_END();
	
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_ARM_SFRMCALIBRATE" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


__declspec(dllexport) LPXLOPER WINAPI Local_PXL_ARM_SFRMCALIBRATE(LPXLOPER XL_calibratorSFRM,
                                                                  LPXLOPER XL_kerneltoGP,
                                                                  LPXLOPER XL_MktCapModel,
                                                                  LPXLOPER XL_MktSwaptModel,
                                                                  LPXLOPER XL_ToCalibrateBeta,
                                                                  LPXLOPER XL_ToCalibrateMR)
{
	ADD_LOG("Local_PXL_ARM_SFRMCALIBRATE");
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
	CCString C_calibratorSFRM;
	XL_readStrCell(XL_calibratorSFRM,C_calibratorSFRM," ARM_ERR: calibrator SFRM id: object expected",C_result);

    // C variable
	CCString C_mktCapModel;
    XL_readStrCellWD(XL_MktCapModel,C_mktCapModel,"DEFAULT"," ARM_ERR: MktCapModel id: object expected",C_result);
    long C_mktCapModelid = C_mktCapModel == "DEFAULT" ? ARM_NULL_OBJECT : LocalGetNumObjectId(C_mktCapModel);

    // C variable
	CCString C_mktSwaptModel;
    XL_readStrCellWD(XL_MktSwaptModel,C_mktSwaptModel,"DEFAULT"," ARM_ERR: MktSwaptModel id: object expected",C_result);
    long C_mktSwapModelid = C_mktSwaptModel == "DEFAULT" ? ARM_NULL_OBJECT : LocalGetNumObjectId(C_mktSwaptModel);

        // C variable
	CCString C_ToCalibrateBeta;
    XL_readStrCellWD(XL_ToCalibrateBeta,C_ToCalibrateBeta,"YES"," ARM_ERR: Flag ToCalibrateBeta: string expected",C_result);

    // C variable
	CCString C_ToCalibrateMR;
    XL_readStrCellWD(XL_ToCalibrateMR,C_ToCalibrateMR,"YES"," ARM_ERR: Flag ToCalibrateMR: string expected",C_result);

    
    double default_kerneltoGP = 0;
    double C_kerneltoGP;
    XL_readNumCellWD(XL_kerneltoGP,C_kerneltoGP,default_kerneltoGP," ARM_ERR:flag to choose KernelTo GP: numeric expected",C_result);

	long retCode;
	long objId;
	CCString prevClass;

	CCString curClass = LOCAL_FRMMODEL_CLASS;
	CCString stringId = GetLastCurCellEnvValue ();

	retCode = ARMLOCAL_SFRMCALIBRATE (LocalGetNumObjectId (C_calibratorSFRM),
                                        C_mktCapModelid,
                                        C_mktSwapModelid,
                                        C_ToCalibrateBeta,
                                        C_ToCalibrateMR,
                                        C_kerneltoGP,
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_PXL_ARM_SFRMCALIBRATE" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}



__declspec(dllexport) LPXLOPER WINAPI Local_ARM_CRACALCULATOR(LPXLOPER XL_sec,
															  LPXLOPER XL_zc,
															  LPXLOPER XL_capVolATM,
															  LPXLOPER XL_rho,
															  LPXLOPER XL_nu,
															  LPXLOPER XL_swptVol,
															  LPXLOPER XL_meanRev,
															  LPXLOPER XL_calibParams,
															  LPXLOPER XL_preInitFlag,
															  LPXLOPER XL_initSigmaCurve,
															  LPXLOPER XL_initBetaOrShift,
															  LPXLOPER XL_voltype,
                                                              LPXLOPER XL_SecurityParams,
															  LPXLOPER XL_horizon,
															  LPXLOPER XL_pathNumber,
															  LPXLOPER XL_TreeParams,
															  LPXLOPER XL_CalibVect,
															  LPXLOPER XL_CalswaptATM,
															  LPXLOPER XL_betaSABR,
															  LPXLOPER XL_SwoptSABRParams)
{
	ADD_LOG("Local_ARM_CRACALCULATOR");
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	    ARM_NOCALCIFWIZ();

	    // C variable
	    CCString C_sec;
	    CCString C_zc;
	    CCString C_capvolat;
        CCString C_rho;
	    long rhoId;
        CCString C_nu;
	    long nuId;
		CCString C_betaSABR;
		long betaSABRId;
	    CCString C_swpvolat;
	    
	    double C_meanRev;
	    
        VECTOR<double> C_calibParams;
	    VECTOR<double> C_initSigmaCurve;
	    VECTOR<double> C_initSigmaCurve_default;
	    VECTOR<double> C_initBetaOrShift;
	    VECTOR<double> C_initBetaOrShift_default;
        VECTOR<double> C_SecurityParams;
	    VECTOR<double> C_SecurityParams_default;
	    VECTOR<double> C_ReCalibParams;
	    VECTOR<double> C_ReCalibParams_default;
	    C_ReCalibParams_default.push_back(1);
	    C_ReCalibParams_default.push_back(1);
		VECTOR<CCString> C_SwoptSABRParams;
		VECTOR<CCString> C_SwoptSABRParams_default;
		
        long rhoSwoptId  = ARM_NULL_OBJECT;
		long nuSwoptId   = ARM_NULL_OBJECT;
		long betaSwoptId = ARM_NULL_OBJECT;
	    
        long SABRSigmaOrAlpha = 1;

	    double C_horizon;
	    double C_pathNumber;
	    VECTOR<double> C_TreeParams; 

	    CCString C_preInitFlag;
	    long preInitFlagId;

	    CCString C_voltype;
	    long voltypeId;
	    
	    CCString C_CalswaptATM;
	    
	    // error
	    static int error;
	    static char* reason = "";

	    XL_readStrCell(XL_sec,C_sec," ARM_ERR: security id: object expected",C_result);
	    XL_readStrCell(XL_zc,C_zc," ARM_ERR: zc: object expected",C_result);
	    XL_readStrCell(XL_capVolATM,C_capvolat," ARM_ERR: Cap Volatility: object expected",C_result);
	    XL_readStrCellWD(XL_rho,C_rho,"DEFAULT"," ARM_ERR: rho(Rho): object expected",C_result);	
	    XL_readStrCellWD(XL_nu,C_nu,"DEFAULT"," ARM_ERR: nu(Nu): object expected",C_result);
        XL_readStrCell(XL_swptVol,C_swpvolat," ARM_ERR: Swaption Volatility: object expected",C_result);
		XL_readStrCellWD(XL_betaSABR,C_betaSABR,"DEFAULT"," ARM_ERR: beta(beta): object expected",C_result);

	    XL_readNumCell(XL_meanRev,C_meanRev," ARM_ERR: mean reversion: numeric expected",C_result);
	    XL_readNumVector(XL_calibParams,C_calibParams," ARM_ERR: calib params : array of numeric expected",C_result);

	    XL_readStrCellWD(XL_preInitFlag,C_preInitFlag,"NO"," ARM_ERR: pre Init Flag: string expected",C_result);
	    XL_readNumVectorWD(XL_initSigmaCurve,C_initSigmaCurve,C_initSigmaCurve_default," ARM_ERR: init Sigma Curve : array of numeric expected",C_result);
	    XL_readNumVectorWD(XL_initBetaOrShift,C_initBetaOrShift,C_initBetaOrShift_default," ARM_ERR: init Beta or Shift : array of numeric expected",C_result);

	    XL_readStrCellWD(XL_voltype,C_voltype,"DIAG"," ARM_ERR: vol type: string expected",C_result);

	    XL_readNumVectorWD(XL_SecurityParams,C_SecurityParams,C_SecurityParams_default," ARM_ERR: calib Params: array of numeric expected",C_result);

	    XL_readNumCell(XL_horizon,C_horizon," ARM_ERR: start date: numeric expected",C_result);
  	    XL_readNumCell(XL_pathNumber,C_pathNumber," ARM_ERR: pathNumber: numeric expected",C_result);
	    XL_readNumVector(XL_TreeParams,C_TreeParams," ARM_ERR: Tree Params: array of numeric expected",C_result);
	    XL_readNumVectorWD(XL_CalibVect,C_ReCalibParams,C_ReCalibParams_default, " ARM_ERR: Re calib params flag: array of numeric expected",C_result);
	    XL_readStrCellWD(XL_CalswaptATM,C_CalswaptATM,"NO"," ARM_ERR: CalswaptATM: string expected",C_result);
		XL_readStrVectorWD(XL_SwoptSABRParams,C_SwoptSABRParams,C_SwoptSABRParams_default," ARM_ERR: swaption SABR parameters: objects expected",DOUBLE_TYPE,C_result);	

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

		if ( C_betaSABR == "DEFAULT" )
	    {
	       betaSABRId = ARM_NULL_OBJECT;
	    }
	    else
	    {
	       betaSABRId = LocalGetNumObjectId(C_betaSABR);
	    }

	    if ((voltypeId = ARM_ConvShapeType (C_voltype, C_result)) == ARM_DEFAULT_ERR)
	    {
	       ARM_ARG_ERR();
	       return (LPXLOPER)&XL_result;
	    }

	    if ((preInitFlagId = ARM_ConvYesOrNo (C_preInitFlag, C_result)) == ARM_DEFAULT_ERR)
	    {
	       ARM_ARG_ERR();
	       return (LPXLOPER)&XL_result;
	    }

		if (C_SwoptSABRParams.size() != 0)
		{
			switch (C_SwoptSABRParams.size())
			{
				case 1:
				{
					rhoSwoptId = LocalGetNumObjectId(C_SwoptSABRParams[0]);
					break;
				}

				case 2:
				{
					rhoSwoptId = LocalGetNumObjectId(C_SwoptSABRParams[0]);
					nuSwoptId  = LocalGetNumObjectId(C_SwoptSABRParams[1]);
					break;
				}

				case 3:
				{
					rhoSwoptId  = LocalGetNumObjectId(C_SwoptSABRParams[0]);
					nuSwoptId   = LocalGetNumObjectId(C_SwoptSABRParams[1]);
					betaSwoptId = LocalGetNumObjectId(C_SwoptSABRParams[2]);
					break;
				}

                case 4:
                {
                    rhoSwoptId = LocalGetNumObjectId(C_SwoptSABRParams[0]);
					nuSwoptId = LocalGetNumObjectId(C_SwoptSABRParams[1]);
					betaSwoptId = LocalGetNumObjectId(C_SwoptSABRParams[2]);

                    if (( C_SwoptSABRParams[3] == "s" )
                        ||
                        ( C_SwoptSABRParams[3] == "S" )
                        ||
                        ( C_SwoptSABRParams[3] == "SIGMA" )
                        ||
                        ( C_SwoptSABRParams[3] == "Sigma" )
                        ||
                        ( C_SwoptSABRParams[3] == "sigma" )
                       )
                    {
                       SABRSigmaOrAlpha = 1;
                    }
                    else
                    {
                       SABRSigmaOrAlpha = 0;
                    }
                }
			}
		}

	    long retCode;
	    long objId;
	    CCString prevClass;

	    CCString curClass = LOCAL_FRM_MARKOVTREE_CLASS;
	    CCString stringId = GetLastCurCellEnvValue ();

	    if (!stringId)
	    {
		    retCode = ARMLOCAL_CRACALCULATOR(LocalGetNumObjectId (C_sec),
										    LocalGetNumObjectId (C_zc),
										    LocalGetNumObjectId (C_capvolat),
										    rhoId,
										    nuId,
										    LocalGetNumObjectId (C_swpvolat),
											betaSABRId,
										    C_meanRev,
									        C_calibParams,
									        preInitFlagId,
									        C_initSigmaCurve,
									        C_initBetaOrShift,
									        voltypeId,
                                            C_SecurityParams,
										    C_horizon,
										    (long)C_pathNumber,
										    C_TreeParams,
										    C_ReCalibParams,
										    C_CalswaptATM,
											rhoSwoptId,
											nuSwoptId,
											betaSwoptId,
                                            SABRSigmaOrAlpha,
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
			    
		    if(curClass == prevClass)
		    {
			    retCode = ARMLOCAL_CRACALCULATOR(LocalGetNumObjectId (C_sec),
											    LocalGetNumObjectId (C_zc),
											    LocalGetNumObjectId (C_capvolat),
											    rhoId,
											    nuId,
											    LocalGetNumObjectId (C_swpvolat),
												betaSABRId,
											    C_meanRev,
											    C_calibParams,
											    preInitFlagId,
											    C_initSigmaCurve,
											    C_initBetaOrShift,
											    voltypeId,
											    C_SecurityParams,
											    C_horizon,
											    (long)C_pathNumber,
											    C_TreeParams,
											    C_ReCalibParams,
											    C_CalswaptATM,
												rhoSwoptId,
												nuSwoptId,
												betaSwoptId,
                                                SABRSigmaOrAlpha,
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

			    retCode = ARMLOCAL_CRACALCULATOR(LocalGetNumObjectId (C_sec),
											    LocalGetNumObjectId (C_zc),
											    LocalGetNumObjectId (C_capvolat),
											    rhoId,
											    nuId,
											    LocalGetNumObjectId (C_swpvolat),
												betaSABRId,
											    C_meanRev,
											    C_calibParams,
											    preInitFlagId,
											    C_initSigmaCurve,
											    C_initBetaOrShift,
											    voltypeId,
											    C_SecurityParams,
											    C_horizon,
											    (long)C_pathNumber,
											    C_TreeParams,
											    C_ReCalibParams,
											    C_CalswaptATM,
												rhoSwoptId,
												nuSwoptId,
												betaSwoptId,
                                                SABRSigmaOrAlpha,
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

    //	ARM_END();
	
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_ARM_CRACALCULATOR" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}



__declspec(dllexport) LPXLOPER WINAPI Local_PXL_ARM_CRACALCULATOR(LPXLOPER XL_sec,
																  LPXLOPER XL_zc,
																  LPXLOPER XL_capVolATM,
																  LPXLOPER XL_rho,
																  LPXLOPER XL_nu,
																  LPXLOPER XL_swptVol,
																  LPXLOPER XL_meanRev,
																  LPXLOPER XL_calibParams,
																  LPXLOPER XL_preInitFlag,
																  LPXLOPER XL_initSigmaCurve,
																  LPXLOPER XL_initBetaOrShift,
																  LPXLOPER XL_voltype,
																  LPXLOPER XL_SecurityParams,
																  LPXLOPER XL_horizon,
																  LPXLOPER XL_pathNumber,
																  LPXLOPER XL_TreeParams,
																  LPXLOPER XL_CalibVect,
																  LPXLOPER XL_CalswaptATM,
																  LPXLOPER XL_betaSABR,
																  LPXLOPER XL_SwoptSABRParams)
{
	ADD_LOG("Local_PXL_ARM_CRACALCULATOR");
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	    ARM_NOCALCIFWIZ();

	    // C variable
	    CCString C_sec;
	    CCString C_zc;
	    CCString C_capvolat;
        CCString C_rho;
        long rhoId;
        CCString C_nu;
	    long nuId;
		CCString C_betaSABR;
		long betaSABRId;
	    CCString C_swpvolat;
	    
	    double C_meanRev;

	    VECTOR<double> C_calibParams;
	    VECTOR<double> C_initSigmaCurve;
	    VECTOR<double> C_initSigmaCurve_default;
	    VECTOR<double> C_initBetaOrShift;
	    VECTOR<double> C_initBetaOrShift_default;
        VECTOR<double> C_SecurityParams;
	    VECTOR<double> C_SecurityParams_default;
	    VECTOR<double> C_ReCalibParams;
	    
        VECTOR<double> C_ReCalibParams_default;
	    C_ReCalibParams_default.push_back(1);
	    C_ReCalibParams_default.push_back(1);
	    
	    
	    double C_horizon;
	    double C_pathNumber;
	    VECTOR<double> C_TreeParams; 

	    CCString C_preInitFlag;
	    long preInitFlagId;

	    CCString C_voltype;
	    long voltypeId;

	    CCString C_CalswaptATM;
	    
		VECTOR<CCString> C_SwoptSABRParams;
		VECTOR<CCString> C_SwoptSABRParams_default;
		
        long rhoSwoptId  = ARM_NULL_OBJECT;
		long nuSwoptId   = ARM_NULL_OBJECT;
		long betaSwoptId = ARM_NULL_OBJECT;

        long SABRSigmaOrAlpha = 1;

		// error
	    static int error;
	    static char* reason = "";

	    XL_readStrCell(XL_sec,C_sec," ARM_ERR: security id: object expected",C_result);
	    XL_readStrCell(XL_zc,C_zc," ARM_ERR: zc: object expected",C_result);
	    XL_readStrCell(XL_capVolATM,C_capvolat," ARM_ERR: Cap Volatility: object expected",C_result);
	    XL_readStrCellWD(XL_rho,C_rho,"DEFAULT"," ARM_ERR: rho(Rho): object expected",C_result);	
	    XL_readStrCellWD(XL_nu,C_nu,"DEFAULT"," ARM_ERR: nu(Nu): object expected",C_result);
        XL_readStrCell(XL_swptVol,C_swpvolat," ARM_ERR: Swaption Volatility: object expected",C_result);
		XL_readStrCellWD(XL_betaSABR,C_betaSABR,"DEFAULT"," ARM_ERR: beta(beta): object expected",C_result);

	    XL_readNumCell(XL_meanRev,C_meanRev," ARM_ERR: mean reversion: numeric expected",C_result);
	    XL_readNumVector(XL_calibParams,C_calibParams," ARM_ERR: calib params : array of numeric expected",C_result);

	    XL_readStrCellWD(XL_preInitFlag,C_preInitFlag,"NO"," ARM_ERR: pre Init Flag: string expected",C_result);
	    XL_readNumVectorWD(XL_initSigmaCurve,C_initSigmaCurve,C_initSigmaCurve_default," ARM_ERR: init Sigma Curve : array of numeric expected",C_result);
	    XL_readNumVectorWD(XL_initBetaOrShift,C_initBetaOrShift,C_initBetaOrShift_default," ARM_ERR: init Beta or Shift : array of numeric expected",C_result);

	    XL_readStrCellWD(XL_voltype,C_voltype,"DIAG"," ARM_ERR: vol type: string expected",C_result);

	    XL_readNumVectorWD(XL_SecurityParams,C_SecurityParams,C_SecurityParams_default," ARM_ERR: calib Params: array of numeric expected",C_result);

	    XL_readNumCell(XL_horizon,C_horizon," ARM_ERR: start date: numeric expected",C_result);
  	    XL_readNumCell(XL_pathNumber,C_pathNumber," ARM_ERR: pathNumber: numeric expected",C_result);
	    XL_readNumVector(XL_TreeParams,C_TreeParams," ARM_ERR: Tree Params: array of numeric expected",C_result);
	    XL_readNumVectorWD(XL_CalibVect,C_ReCalibParams,C_ReCalibParams_default, " ARM_ERR: Re calib params flag: array of numeric expected",C_result);
	    XL_readStrCellWD(XL_CalswaptATM,C_CalswaptATM,"NO"," ARM_ERR: CalswaptATM: string expected",C_result);
		XL_readStrVectorWD(XL_SwoptSABRParams,C_SwoptSABRParams,C_SwoptSABRParams_default," ARM_ERR: swaption SABR parameters: objects expected",DOUBLE_TYPE,C_result);	

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

		if ( C_betaSABR == "DEFAULT" )
	    {
	       betaSABRId = ARM_NULL_OBJECT;
	    }
	    else
	    {
	       betaSABRId = LocalGetNumObjectId(C_betaSABR);
	    }

	    if ((voltypeId = ARM_ConvShapeType (C_voltype, C_result)) == ARM_DEFAULT_ERR)
	    {
	       ARM_ARG_ERR();
	       return (LPXLOPER)&XL_result;
	    }

	    if ((preInitFlagId = ARM_ConvYesOrNo (C_preInitFlag, C_result)) == ARM_DEFAULT_ERR)
	    {
	       ARM_ARG_ERR();
	       return (LPXLOPER)&XL_result;
	    }

		if (C_SwoptSABRParams.size() != 0)
		{
			switch (C_SwoptSABRParams.size())
			{
				case 1:
				{
					rhoSwoptId = LocalGetNumObjectId(C_SwoptSABRParams[0]);
				};
                break;

				case 2:
				{
					rhoSwoptId = LocalGetNumObjectId(C_SwoptSABRParams[0]);
					nuSwoptId  = LocalGetNumObjectId(C_SwoptSABRParams[1]);
				};
                break;

				case 3:
				{
					rhoSwoptId  = LocalGetNumObjectId(C_SwoptSABRParams[0]);
					nuSwoptId   = LocalGetNumObjectId(C_SwoptSABRParams[1]);
					betaSwoptId = LocalGetNumObjectId(C_SwoptSABRParams[2]);
				};
                break;

                case 4:
                {
                    rhoSwoptId  = LocalGetNumObjectId(C_SwoptSABRParams[0]);
					nuSwoptId   = LocalGetNumObjectId(C_SwoptSABRParams[1]);
					betaSwoptId = LocalGetNumObjectId(C_SwoptSABRParams[2]);

                    if (( C_SwoptSABRParams[3] == "s" )
                        ||
                        ( C_SwoptSABRParams[3] == "S" )
                        ||
                        ( C_SwoptSABRParams[3] == "SIGMA" )
                        ||
                        ( C_SwoptSABRParams[3] == "Sigma" )
                        ||
                        ( C_SwoptSABRParams[3] == "sigma" )
                       )
                    {
                       SABRSigmaOrAlpha = 1;
                    }
                    else
                    {
                       SABRSigmaOrAlpha = 0;
                    }
                }

				default:
				{
				   ARM_ARG_ERR();
				   return (LPXLOPER)&XL_result;
				}
			}
		}

	    long retCode;
	    long objId;
	    CCString prevClass;

	    CCString curClass = LOCAL_FRM_MARKOVTREE_CLASS;
	    CCString stringId = GetLastCurCellEnvValue ();

	    retCode = ARMLOCAL_CRACALCULATOR(LocalGetNumObjectId (C_sec),
									    LocalGetNumObjectId (C_zc),
									    LocalGetNumObjectId (C_capvolat),
									    rhoId,
									    nuId,
									    LocalGetNumObjectId (C_swpvolat),
										betaSABRId,
									    C_meanRev,
									    C_calibParams,
									    preInitFlagId,
									    C_initSigmaCurve,
									    C_initBetaOrShift,
									    voltypeId,
                                        C_SecurityParams,
									    C_horizon,
									    (long)C_pathNumber,
									    C_TreeParams,
									    C_ReCalibParams,
									    C_CalswaptATM,
										rhoSwoptId,
										nuSwoptId,
										betaSwoptId,
                                        SABRSigmaOrAlpha,
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_PXL_ARM_CRACALCULATOR" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}



class craGetModelFunc : public ARMResultLong2LongFunc
{
public:
	craGetModelFunc(long craId)
    :
    C_craId(craId)
    {};
	
	long operator()( ARM_result& result, long objId ){
		return ARMLOCAL_CRA_GetModel(
            C_craId,
            result,
            objId);
    }

private:
	long    C_craId;
};



LPXLOPER Local_CraCalculator_GetModel_Common(
	LPXLOPER XL_craId,
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

		long craId;
		XL_GETOBJID( XL_craId, craId, "ARM_ERR: o CRA Calculator: Object expected", C_result);

		craGetModelFunc ourFunc(craId);

		/// call the general function
		fillXL_Result( LOCAL_SFRMMODEL_CLASS, ourFunc, C_result, XL_result, PersistentInXL );
	}

	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_CraCalculator_GetModel_Common" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}

__declspec(dllexport) LPXLOPER WINAPI Local_ARM_CRACALCULATOR_GetModel(
	LPXLOPER XL_craId)
{
	ADD_LOG("Local_ARM_CRACALCULATOR_GetModel");
	bool PersistentInXL = true;
	return Local_CraCalculator_GetModel_Common(
	    XL_craId,
	    PersistentInXL );
}


__declspec(dllexport) LPXLOPER WINAPI Local_PXL_ARM_CRACALCULATOR_GetModel(
	LPXLOPER XL_craId)
{
	ADD_LOG("Local_PXL_ARM_CRACALCULATOR_GetModel");
	bool PersistentInXL = false;
	return Local_CraCalculator_GetModel_Common(
	    XL_craId,
        PersistentInXL );
}


class craGetCalibDataFunc : public ARMResultLong2LongFunc
{
public:
	craGetCalibDataFunc(long craId,
						string calibOrPortfolio,
						string dataType)
    :
    C_craId(craId),
	C_calibOrPortfolio(calibOrPortfolio),
	C_dataType(dataType)
    {};
	
	long operator()( ARM_result& result, long objId ){
		return ARMLOCAL_CRA_GetCalibData(
            C_craId,
			C_calibOrPortfolio,
			C_dataType,
            result,
            objId);
    }

private:
	long    C_craId;
	string	C_calibOrPortfolio;
	string	C_dataType;
};

LPXLOPER Local_CraCalculator_GetCalibData_Common(
	LPXLOPER XL_craId,
	LPXLOPER XL_calibOrPortfolio,
	LPXLOPER XL_dataType,
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

		//CRA ID
		long craId;
		XL_GETOBJID( XL_craId, craId, "ARM_ERR: CRA Calculator: Object expected", C_result);

		//CALIB OR PORTFOLIO
		CCString C_calibOrPortfolio;
		XL_readStrCell(XL_calibOrPortfolio,C_calibOrPortfolio," ARM_ERR: calib or portfolio: string expected",C_result);
		C_calibOrPortfolio.toUpper ();
		string calibOrPortfolio = CCSTringToSTLString(C_calibOrPortfolio);
	
		//DATA TYPE
		CCString C_dataType;
		XL_readStrCell(XL_dataType,C_dataType," ARM_ERR: data type: string expected",C_result);
		C_dataType.toUpper();
		string dataType = CCSTringToSTLString(C_dataType);

		if ((calibOrPortfolio == "CALIB") || (calibOrPortfolio == "PORTFOLIO"))
		{
			if ((dataType == "SIGMA") || (dataType == "BETA") || (dataType == "MRS"))
			{
				craGetCalibDataFunc ourFunc(craId, calibOrPortfolio, dataType);

				if (calibOrPortfolio == "CALIB")
					fillXL_Result( LOCAL_CALIBMETHOD_CLASS, ourFunc, C_result, XL_result, PersistentInXL );
				else if (calibOrPortfolio == "PORTFOLIO")
					fillXL_Result( LOCAL_PF_CLASS, ourFunc, C_result, XL_result, PersistentInXL );
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
	}

	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Local_CraCalculator_GetCalibData_Common: Did you enter a good calibOrPortfolioFlag and good dataType ?" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}

__declspec(dllexport) LPXLOPER WINAPI Local_ARM_CRACALCULATOR_GetCalibData(
	LPXLOPER XL_craId,
	LPXLOPER XL_calibOrPortfolio,
	LPXLOPER XL_dataType)
{
	ADD_LOG("Local_ARM_CRACALCULATOR_GetCalibData");
	bool PersistentInXL = true;
	return Local_CraCalculator_GetCalibData_Common(
	    XL_craId,
		XL_calibOrPortfolio,
		XL_dataType,
	    PersistentInXL );
}


__declspec(dllexport) LPXLOPER WINAPI Local_PXL_ARM_CRACALCULATOR_GetCalibData(
	LPXLOPER XL_craId,
	LPXLOPER XL_calibOrPortfolio,
	LPXLOPER XL_dataType)
{
	ADD_LOG("Local_PXL_ARM_CRACALCULATOR_GetCalibData");
	bool PersistentInXL = false;
	return Local_CraCalculator_GetCalibData_Common(
	    XL_craId,
		XL_calibOrPortfolio,
		XL_dataType,
        PersistentInXL );
}


class craGetUnderlyingFunc : public ARMResultLong2LongFunc
{
public:
	craGetUnderlyingFunc(long craId)
    :
    C_craId(craId)
    {};
	
	long operator()( ARM_result& result, long objId ){
		return ARMLOCAL_CRA_GetUnderlying(  C_craId,
											result,
											objId);
    }

private:
	long    C_craId;
};


LPXLOPER Local_CraCalculator_GetUnderlying_Common(
	LPXLOPER XL_craId,
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

		//CRA ID
		long craId;
		XL_GETOBJID( XL_craId, craId, "ARM_ERR: CRA Calculator: Object expected", C_result);

		craGetUnderlyingFunc ourFunc(craId);

		fillXL_Result( LOCAL_SWAP_CLASS, ourFunc, C_result, XL_result, PersistentInXL );
	}

	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception 
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Local_CraCalculator_GetUnderlying: error, please contact R & D desk" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


__declspec(dllexport) LPXLOPER WINAPI Local_ARM_CRACALCULATOR_GetUnderlying(
	LPXLOPER XL_craId)
{
	ADD_LOG("Local_ARM_CRACALCULATOR_GetUnderlying");
	bool PersistentInXL = true;
	return Local_CraCalculator_GetUnderlying_Common(
	    XL_craId,
	    PersistentInXL );
}

__declspec(dllexport) LPXLOPER WINAPI Local_PXL_ARM_CRACALCULATOR_GetUnderlying(
	LPXLOPER XL_craId)
{
	ADD_LOG("Local_PXL_ARM_CRACALCULATOR_GetUnderlying");
	bool PersistentInXL = false;
	return Local_CraCalculator_GetUnderlying_Common(
	    XL_craId,
        PersistentInXL );
}



__declspec(dllexport) LPXLOPER WINAPI Local_ARM_BUMPCRACALCULATOR(LPXLOPER XL_cracalc,
																  LPXLOPER XL_zc,
																  LPXLOPER XL_capVolATM,
																  LPXLOPER XL_rho,
																  LPXLOPER XL_nu,
																  LPXLOPER XL_swptVol,
																  LPXLOPER XL_betaSABR,
																  LPXLOPER XL_rhoSwopt,
																  LPXLOPER XL_nuSwopt,
																  LPXLOPER XL_betaSwopt,
                                                                  LPXLOPER XL_SigmaOrAlpha)
{
	ADD_LOG("Local_ARM_BUMPCRACALCULATOR");
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	    ARM_NOCALCIFWIZ();

	    // C variable
	    CCString C_cracalc;
	    CCString C_zc;
	    CCString C_capvolat;
        
        CCString C_rho;
	    long rhoId;
        
        CCString C_nu;
	    long nuId;

	    CCString C_betaSABR;
	    long betaSABRId;
	    
        CCString C_swpvolat;
        
        CCString C_rhoSwopt;
	    long rhoSwoptId;
        
        CCString C_nuSwopt;
	    long nuSwoptId;

	    CCString C_betaSwopt;
	    long betaSwoptId;

        CCString C_SigmaOrAlpha;

	    // error
	    static int error;
	    static char* reason = "";

	    XL_readStrCell(XL_cracalc,C_cracalc," ARM_ERR: Cra Calculator id: object expected",C_result);
	    XL_readStrCell(XL_zc,C_zc," ARM_ERR: zc: object expected",C_result);
	    XL_readStrCell(XL_capVolATM,C_capvolat," ARM_ERR: Cap Volatility: object expected",C_result);
	    XL_readStrCellWD(XL_rho,C_rho,"DEFAULT"," ARM_ERR: rho(Rho): object expected",C_result);	
	    XL_readStrCellWD(XL_nu,C_nu,"DEFAULT"," ARM_ERR: nu(Nu): object expected",C_result);
        XL_readStrCell(XL_swptVol,C_swpvolat," ARM_ERR: Swaption Volatility: object expected",C_result);
	    XL_readStrCellWD(XL_betaSABR,C_betaSABR,"DEFAULT"," ARM_ERR: beta(beta): object expected",C_result);
	    XL_readStrCellWD(XL_rhoSwopt,C_rhoSwopt,"DEFAULT"," ARM_ERR: rho swaption: object expected",C_result);	
	    XL_readStrCellWD(XL_nuSwopt,C_nuSwopt,"DEFAULT"," ARM_ERR: nu swaption: object expected",C_result);
	    XL_readStrCellWD(XL_betaSwopt,C_betaSwopt,"DEFAULT"," ARM_ERR: beta swaption: object expected",C_result);
   
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

	    if ( C_betaSABR == "DEFAULT" )
		   betaSABRId = ARM_NULL_OBJECT;
	    else
		   betaSABRId = LocalGetNumObjectId(C_betaSABR);

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
		   betaSwoptId = ARM_NULL_OBJECT;
	    else
		   betaSwoptId = LocalGetNumObjectId(C_betaSwopt);

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

	    long retCode = ARMLOCAL_BUMPCRACALCULATOR(LocalGetNumObjectId (C_cracalc),
											    LocalGetNumObjectId (C_zc),
											    LocalGetNumObjectId (C_capvolat),
											    rhoId,
											    nuId,
											    LocalGetNumObjectId (C_swpvolat),
											    betaSABRId,
											    rhoSwoptId,
											    nuSwoptId,
											    betaSwoptId,
                                                SigmaOrAlpha,
											    C_result);

	    if ( retCode == ARM_OK )
	    {
		    FreeCurCellErr();
	    
            XL_result.xltype = xltypeStr;
		    XL_result.val.str = XL_StrC2StrPascal(C_cracalc);
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_ARM_BUMPCRACALCULATOR" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}



__declspec(dllexport) LPXLOPER WINAPI Local_ARM_ComputeFwdBSDelta(	LPXLOPER XL_fwd,
																	LPXLOPER XL_strike,
																	LPXLOPER XL_vol,
																	LPXLOPER XL_T,
																	LPXLOPER XL_CallPut)
{
	ADD_LOG("Local_ARM_ComputeFwdBSDelta");
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
		ARM_NOCALCIFWIZ();


		// C variable
		double C_fwd;
		double C_strike;
		double C_vol;
		double C_T;
		double C_CallPut; 
		
		// error
		static int error;
		static char* reason = "";
		
		XL_readNumCell(XL_fwd,C_fwd," ARM_ERR: fwd: fwd expected",C_result);
		XL_readNumCell(XL_strike,C_strike," ARM_ERR: strike: strike expected",C_result);
		XL_readNumCell(XL_vol,C_vol," ARM_ERR: vol: volatility expected",C_result);
		XL_readNumCell(XL_T,C_T," ARM_ERR: T: date expected",C_result);
		XL_readNumCell(XL_CallPut,C_CallPut," ARM_ERR: CallPut: flag Call or Put expected",C_result);
		
		long retCode;

		retCode = ARMLOCAL_COMPUTEFWDBSDELTA(C_fwd, C_strike,C_vol,C_T, int(C_CallPut), C_result);

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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_ARM_ComputeFwdBSDelta" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


__declspec(dllexport) LPXLOPER WINAPI Local_ARM_ComputeSplinedSigmaATMF(LPXLOPER XL_deltas,
																		LPXLOPER XL_sigmas,
																		LPXLOPER XL_matu,
																		LPXLOPER XL_SigmaZDS,
																		LPXLOPER XL_Precision,
                                                                        LPXLOPER XL_FxSpot)
{
	ADD_LOG("Local_ARM_ComputeSplinedSigmaATMF");
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

    ARM_NOCALCIFWIZ();

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
		// C variable
		VECTOR<double> C_deltas; 
		VECTOR<double> C_sigmas; 

		double C_matu;
	
		double C_SigmaZDS;
		double C_SigmaZDS_Default = -1;

		double C_Precision;
		double C_Precision_Default = 1e-4;

        double C_FxSpot;
		double C_FxSpot_Default = -1.0;

		// error
		static int error;
		static char* reason = "";
		
		XL_readNumVector(XL_deltas,C_deltas," ARM_ERR: deltas : array of numeric expected",C_result);
		XL_readNumVector(XL_sigmas,C_sigmas," ARM_ERR: sigmas: array of numeric expected",C_result);
		
		XL_readNumCell(XL_matu,C_matu," ARM_ERR: matu: matu expected",C_result);
		XL_readNumCellWD(XL_SigmaZDS,C_SigmaZDS,C_SigmaZDS_Default, " ARM_ERR: SigmaZDS: SigmaZDS expected",C_result);
		XL_readNumCellWD(XL_Precision,C_Precision,C_Precision_Default," ARM_ERR: Precision: Precision expected",C_result);
		
        XL_readNumCellWD(XL_FxSpot, C_FxSpot, C_FxSpot_Default, " ARM_ERR: FX Spot expected",C_result);

		long retCode;

		retCode = ARMLOCAL_ComputeSplinedSigmaATMF(C_deltas, C_sigmas, C_matu, C_SigmaZDS, 
                                                   C_Precision, C_FxSpot, C_result);

		if ( retCode == ARM_OK )
		{
			FreeCurCellErr();
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_ARM_ComputeFwdBSDelta" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


__declspec(dllexport) LPXLOPER WINAPI Local_ARM_ComputeDeltaFwdFromDeltaWP
																	(	LPXLOPER XL_AsOf,
																		LPXLOPER XL_matu,
																		LPXLOPER XL_sigma,
																		LPXLOPER XL_fxSpot,
																		LPXLOPER XL_deltaWithPremium,
																		LPXLOPER XL_domCrvId,
																		LPXLOPER XL_foreignCrvId)
{
	ADD_LOG("Local_ARM_ComputeDeltaFwdFromDeltaWP");
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

    ARM_NOCALCIFWIZ();
	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
		// C variable
		double C_AsOf;
		double C_matu;
		double C_sigma;
		double C_fxSpot;
		double C_deltaWithPremium;
		
		CCString C_domCrvId;
		CCString C_foreignCrvId;


		// error
		static int error;
		static char* reason = "";
		
		
		XL_readNumCell(XL_AsOf,C_AsOf," ARM_ERR: As of date: date expected",C_result);

		XL_readNumCell(XL_matu,C_matu," ARM_ERR: matu: matu expected",C_result);
		
		XL_readNumCell(XL_sigma,C_sigma, " ARM_ERR: sigma: sigma expected",C_result);

		XL_readNumCell(XL_fxSpot,C_fxSpot," ARM_ERR: fxSpot: fxSpot expected",C_result);

		XL_readNumCell(XL_deltaWithPremium,C_deltaWithPremium," ARM_ERR: delta: delta expected",C_result);

		XL_readStrCell(XL_domCrvId,C_domCrvId," ARM_ERR: dom curve id: object expected",C_result);

		XL_readStrCell(XL_foreignCrvId,C_foreignCrvId," ARM_ERR:foreign curve id: object expected",C_result);
		
		long retCode;

		retCode = ARMLOCAL_ComputeDeltaFwdFromDeltaWP(C_AsOf, C_matu, C_sigma, C_fxSpot, 
													  C_deltaWithPremium, 
													  LocalGetNumObjectId(C_domCrvId),
													  LocalGetNumObjectId(C_foreignCrvId),
													  C_result);

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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_ARM_ComputeDeltaFwdFromDeltaWP" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}



_declspec(dllexport) LPXLOPER WINAPI Local_ARM_CalculateImpliedStrikeFromDeltaWithPremium  (LPXLOPER XL_AsOf,
																							LPXLOPER XL_matu,
																							LPXLOPER XL_sigma,
																							LPXLOPER XL_fxSpot,
																							LPXLOPER XL_deltaWithPremium,
																							LPXLOPER XL_domCrvId,
																							LPXLOPER XL_foreignCrvId)
{
	//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

    ARM_NOCALCIFWIZ();
	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
		// C variable
		double C_AsOf;
		double C_matu;
		double C_sigma;
		double C_fxSpot;
		double C_deltaWithPremium;
		
		CCString C_domCrvId;
		CCString C_foreignCrvId;


		// error
		static int error;
		static char* reason = "";
		
		
		XL_readNumCell(XL_AsOf,C_AsOf," ARM_ERR: As of date: date expected",C_result);

		XL_readNumCell(XL_matu,C_matu," ARM_ERR: matu: matu expected",C_result);
		
		XL_readNumCell(XL_sigma,C_sigma, " ARM_ERR: sigma: sigma expected",C_result);

		XL_readNumCell(XL_fxSpot,C_fxSpot," ARM_ERR: fxSpot: fxSpot expected",C_result);

		XL_readNumCell(XL_deltaWithPremium,C_deltaWithPremium," ARM_ERR: delta: delta expected",C_result);

		XL_readStrCell(XL_domCrvId,C_domCrvId," ARM_ERR: dom curve id: object expected",C_result);

		XL_readStrCell(XL_foreignCrvId,C_foreignCrvId," ARM_ERR:foreign curve id: object expected",C_result);
		
		long retCode;

		retCode = ARMLOCAL_CalculateImpliedStrikeFromDeltaWithPremium(C_AsOf, C_matu, C_sigma, C_fxSpot, 
																	  C_deltaWithPremium, 
																	  LocalGetNumObjectId(C_domCrvId),
																	  LocalGetNumObjectId(C_foreignCrvId),
																	  C_result);

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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_ARM_CalculateImpliedStrikeFromDeltaWithPremium" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


_declspec(dllexport) LPXLOPER WINAPI Local_ARM_CalcFwdFXSpot   (LPXLOPER XL_AsOfDate,
																LPXLOPER XL_Spot,
																LPXLOPER XL_aFwdDate,
																LPXLOPER XL_NumDiscountCurveID,
																LPXLOPER XL_UndDiscountCurveID)
{
	//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

    ARM_NOCALCIFWIZ();
	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
		// C variable
		double C_AsOfDate;
		double C_Spot;
		double C_aFwdDate;
		
		CCString C_NumDiscountCurveID;
		CCString C_UndDiscountCurveID;


		// error
		static int error;
		static char* reason = "";
		
		
		XL_readNumCell(XL_AsOfDate,C_AsOfDate," ARM_ERR: As of date: date expected",C_result);

		XL_readNumCell(XL_Spot,C_Spot," ARM_ERR: spot: value expected",C_result);
		
		XL_readNumCell(XL_aFwdDate,C_aFwdDate, " ARM_ERR: Forward Date: date expected",C_result);

		XL_readStrCell(XL_NumDiscountCurveID,C_NumDiscountCurveID," ARM_ERR: Domestic curve id: object expected",C_result);

		XL_readStrCell(XL_UndDiscountCurveID,C_UndDiscountCurveID," ARM_ERR:Foreign curve id: object expected",C_result);
		
		long retCode;

		retCode = ARMLOCAL_CalcFwdFXSpot(C_AsOfDate, C_Spot, C_aFwdDate,  
										  LocalGetNumObjectId(C_NumDiscountCurveID),
										  LocalGetNumObjectId(C_UndDiscountCurveID),
										  C_result);

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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_ARM_CalcFwdFXSpot" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}



__declspec(dllexport) LPXLOPER WINAPI Local_ARM_GetQModelQVolatility(	LPXLOPER XL_QModelId)
{
	ADD_LOG("Local_ARM_GetQModelQVolatility");
	// this is defined first because it is used in XL macros
	static XLOPER XL_result;
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

		long QModelId;
		XL_GETOBJID( XL_QModelId, QModelId,	" ARM_ERR: TARN Calculator: Object expected",C_result);

		/// call the function with nomore LPXLOPER objects
		/// this fction store the result into C_result
		long retCode;
		if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 1)
			retCode = ARMLOCAL_GetQModelQVolatility( QModelId,C_result );
		else
			retCode = ARM_KO;
	
		/// feed the LPXLOPER object result 
		if (retCode == ARM_OK)
		{
			FreeCurCellErr ();
			XL_result.xltype    = xltypeNum;
			XL_result.val.num	= C_result.getDouble();
		}
	
		/// be aware that ARM_ERR is a macro
		/// hence the bracket are necessary
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_ARM_GetQVolatility" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}



_declspec(dllexport) LPXLOPER WINAPI Local_ARM_FromSigmaToAlpha(LPXLOPER XL_bssmiledId,
																LPXLOPER XL_sigmaCurveId,
																LPXLOPER XL_strike)
{
	//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

    ARM_NOCALCIFWIZ();
	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
		// C variable
		CCString C_bssmiledId;
		long bssmiledId;

		CCString C_sigmaCurveId;
		long sigmaCurveId;

		double C_strike;
		double C_strike_default = -1.0;
		
		// error
		static int error;
		static char* reason = "";
		
		XL_readStrCell(XL_bssmiledId,C_bssmiledId," ARM_ERR: bssmiled model id: object expected",C_result);

		XL_readStrCellWD(XL_sigmaCurveId,C_sigmaCurveId,"DEFAULT"," ARM_ERR: sigma curve id: object expected",C_result);
		
		XL_readNumCellWD(XL_strike,C_strike,C_strike_default," ARM_ERR: strike: strike expected",C_result);
		
		if ( C_sigmaCurveId == "DEFAULT" )
		{
		   sigmaCurveId = ARM_NULL_OBJECT;
		}
		else
		{
		   sigmaCurveId = LocalGetNumObjectId(C_sigmaCurveId);
		}
		
		bssmiledId = LocalGetNumObjectId(C_bssmiledId);

		long retCode;
		long objId;
		CCString prevClass;

		CCString curClass = LOCAL_VOL_CURVE_LIN_CLASS;

		CCString stringId = GetLastCurCellEnvValue ();

		if (!stringId)
		{
		   retCode = ARMLOCAL_FromSigmaToAlpha(bssmiledId,
												sigmaCurveId,
												C_strike,
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
			prevClass = LocalGetStringObjectClass(stringId);

			objId = LocalGetNumObjectId(stringId);

			if ( curClass == prevClass )
			{
				retCode = ARMLOCAL_FromSigmaToAlpha(bssmiledId,
													sigmaCurveId,
													C_strike,
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
				
				retCode = ARMLOCAL_FromSigmaToAlpha(bssmiledId,
													sigmaCurveId,
													C_strike,
													C_result);
				if ( retCode == ARM_OK )
				{
					objId = C_result.getLong();

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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_ARM_FromSigmaToAlpha" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}



_declspec(dllexport) LPXLOPER WINAPI Local_ARM_FromAlphaToSigma (LPXLOPER XL_bssmiledId,
																 LPXLOPER XL_alphaCurveId,
																 LPXLOPER XL_strike)
{
	//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

    ARM_NOCALCIFWIZ();
	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
		// C variable
		CCString C_bssmiledId;
		long bssmiledId;

		CCString C_alphaCurveId;
		long alphaCurveId;

		double C_strike;
		double C_strike_default = -1.0;
		
		// error
		static int error;
		static char* reason = "";
		
		XL_readStrCell(XL_bssmiledId,C_bssmiledId," ARM_ERR: bssmile model id: object expected",C_result);

		XL_readStrCellWD(XL_alphaCurveId,C_alphaCurveId,"DEFAULT"," ARM_ERR: alpha curve id: object expected",C_result);
		
		XL_readNumCellWD(XL_strike,C_strike,C_strike_default," ARM_ERR: strike: strike expected",C_result);
		
		if ( C_alphaCurveId == "DEFAULT" )
		{
			alphaCurveId = ARM_NULL_OBJECT;
		}
		else
		{
			alphaCurveId = LocalGetNumObjectId(C_alphaCurveId);
		}
		
		bssmiledId = LocalGetNumObjectId(C_bssmiledId);

		long retCode;
		long objId;
		CCString prevClass;

		CCString curClass = LOCAL_VOL_CURVE_LIN_CLASS;
		CCString stringId = GetLastCurCellEnvValue ();

		if (!stringId)
		{
			retCode = ARMLOCAL_FromAlphaToSigma(bssmiledId,
												alphaCurveId,
												C_strike,
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
			prevClass = LocalGetStringObjectClass(stringId);

			objId = LocalGetNumObjectId(stringId);

			if ( curClass == prevClass )
			{
				retCode = ARMLOCAL_FromAlphaToSigma(bssmiledId,
													alphaCurveId,
													C_strike,
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
				
				retCode = ARMLOCAL_FromAlphaToSigma(bssmiledId,
													alphaCurveId,
													C_strike,
													C_result);
				if ( retCode == ARM_OK )
				{
					objId = C_result.getLong();

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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_ARM_FromAlphaToSigma" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}

_declspec(dllexport) LPXLOPER WINAPI Local_PXL_ARM_FromSigmaToAlpha(LPXLOPER XL_bssmiledId,
																	LPXLOPER XL_sigmaCurveId,
																	LPXLOPER XL_strike)
{
	//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

    ARM_NOCALCIFWIZ();
	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
		// C variable
		CCString C_bssmiledId;
		long bssmiledId;

		CCString C_sigmaCurveId;
		long sigmaCurveId;

		double C_strike;
		double C_strike_default = -1.0;
		
		// error
		static int error;
		static char* reason = "";
		
		XL_readStrCell(XL_bssmiledId,C_bssmiledId," ARM_ERR: bssmiled model id: object expected",C_result);

		XL_readStrCellWD(XL_sigmaCurveId,C_sigmaCurveId,"DEFAULT"," ARM_ERR: sigma curve id: object expected",C_result);
		
		XL_readNumCellWD(XL_strike,C_strike,C_strike_default," ARM_ERR: strike: strike expected",C_result);
		
		if ( C_sigmaCurveId == "DEFAULT" )
		{
		   sigmaCurveId = ARM_NULL_OBJECT;
		}
		else
		{
		   sigmaCurveId = LocalGetNumObjectId(C_sigmaCurveId);
		}
		
		bssmiledId = LocalGetNumObjectId(C_bssmiledId);

		long retCode;
		long objId;

		CCString curClass = LOCAL_VOL_CURVE_LIN_CLASS;

		CCString stringId = GetLastCurCellEnvValue ();

     	retCode = ARMLOCAL_FromSigmaToAlpha(bssmiledId,
											sigmaCurveId,
											C_strike,
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

	//	ARM_END();
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_PXL_ARM_FromSigmaToAlpha" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}



_declspec(dllexport) LPXLOPER WINAPI Local_PXL_ARM_FromAlphaToSigma (LPXLOPER XL_bssmiledId,
																	 LPXLOPER XL_alphaCurveId,
																	 LPXLOPER XL_strike)
{
	//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

    ARM_NOCALCIFWIZ();
	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
		// C variable
		CCString C_bssmiledId;
		long bssmiledId;

		CCString C_alphaCurveId;
		long alphaCurveId;

		double C_strike;
		double C_strike_default = -1.0;
		
		// error
		static int error;
		static char* reason = "";
		
		XL_readStrCell(XL_bssmiledId,C_bssmiledId," ARM_ERR: bssmile model id: object expected",C_result);

		XL_readStrCellWD(XL_alphaCurveId,C_alphaCurveId,"DEFAULT"," ARM_ERR: alpha curve id: object expected",C_result);
		
		XL_readNumCellWD(XL_strike,C_strike,C_strike_default," ARM_ERR: strike: strike expected",C_result);
		
		if ( C_alphaCurveId == "DEFAULT" )
		{
			alphaCurveId = ARM_NULL_OBJECT;
		}
		else
		{
			alphaCurveId = LocalGetNumObjectId(C_alphaCurveId);
		}
		
		bssmiledId = LocalGetNumObjectId(C_bssmiledId);

		long retCode;
		long objId;

		CCString curClass = LOCAL_VOL_CURVE_LIN_CLASS;
		CCString stringId = GetLastCurCellEnvValue ();

		retCode = ARMLOCAL_FromAlphaToSigma(bssmiledId,
											alphaCurveId,
											C_strike,
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

	//	ARM_END();
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_PXL_ARM_FromAlphaToSigma" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


_declspec(dllexport) LPXLOPER WINAPI Local_ExtractTree3FVol(LPXLOPER XL_FxOptIdx)
{
	//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

    ARM_NOCALCIFWIZ();
	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
		static int error;
		static char* reason = "";

        double C_FxOptIdx;
        double C_FxOptIdxDef=0;
		XL_readNumCellWD(XL_FxOptIdx,C_FxOptIdx,C_FxOptIdxDef," ARM_ERR: FxOptIdx: numerical expected",C_result);

        FxOptIdxGlobal = (int)C_FxOptIdx;

        int nbDomRows = DomVolGlobalExtraction.GetNumLines();
        int nbDomCols = DomVolGlobalExtraction.GetNumCols();
        int nbForRows = ForVolGlobalExtraction.GetNumLines();
        int nbForCols = ForVolGlobalExtraction.GetNumCols();
        int nbFxRows = FxVolGlobalExtraction.GetNumLines();
        int nbFxCols = FxVolGlobalExtraction.GetNumCols();
        int nbAnaRows = VolAnalyticsGlobalExtraction.GetNumLines();
        int nbAnaCols = VolAnalyticsGlobalExtraction.GetNumCols();
        int nbTreeVolRows = LatticeVolsGlobalExtraction.GetNumLines();
        int nbTreeVolCols = LatticeVolsGlobalExtraction.GetNumCols();
        int nbNodesRows = NbNodesPerSliceGlobalExtraction.GetSize();
        int nbFxOptionRows = AnalyticFxOptionGlobalExtraction.GetNumLines();
        int nbFxOptionCols = AnalyticFxOptionGlobalExtraction.GetNumCols();

        int nbRows = nbDomRows > nbForRows ? nbDomRows : nbForRows;
        nbRows = nbFxRows > nbRows ? nbFxRows : nbRows;
        nbRows = nbAnaRows > nbRows ? nbAnaRows : nbRows;
        nbRows = nbTreeVolRows > nbRows ? nbTreeVolRows : nbRows;
        nbRows = nbNodesRows > nbRows ? nbNodesRows : nbRows;
        nbRows = nbFxOptionRows > nbRows ? nbFxOptionRows : nbRows;

        if((nbDomRows>0 && nbDomCols!=2) || (nbForRows>0 && nbForCols!=2) ||
            (nbFxRows>0 && nbFxCols!=2) || (nbAnaRows>0 && nbAnaCols!=4) ||
            (nbTreeVolRows>0 && nbTreeVolCols!=7) || (nbFxOptionRows>0 && nbFxOptionCols!=3))
        {
			ARM_ERR();
        }

        int nbCols=nbDomCols+nbForCols+nbFxCols+nbAnaCols+nbTreeVolCols+1+nbFxOptionCols;

		VECTOR<double> vectorResult(nbRows*nbCols,0.0);

        size_t j,offset;
        if(nbRows>0)
        {
		    for ( size_t i=0; i < nbRows; ++i)
            {
                offset=0;
                if(i<nbDomRows)
                {
                    for (j=0; j < nbDomCols; ++j)
			            vectorResult[i*nbCols+offset+j] = DomVolGlobalExtraction.Elt(i,j);
                }
                offset+=nbDomCols;
                if(i<nbForRows)
                {
                    for ( j=0; j < nbForCols; ++j)
			            vectorResult[i*nbCols+offset+j] = ForVolGlobalExtraction.Elt(i,j);
                }
                offset+=nbForCols;
                if(i<nbFxRows)
                {
                    for ( j=0; j < nbFxCols; ++j)
			            vectorResult[i*nbCols+offset+j] = FxVolGlobalExtraction.Elt(i,j);
                }
                offset+=nbFxCols;
                if(i<nbAnaRows)
                {
                    for ( j=0; j < nbAnaCols; ++j)
			            vectorResult[i*nbCols+offset+j] = VolAnalyticsGlobalExtraction.Elt(i,j);
                }
                offset+=nbAnaCols;
                if(i<nbTreeVolRows)
                {
                    for ( j=0; j < nbTreeVolCols; ++j)
    			        vectorResult[i*nbCols+offset+j] = LatticeVolsGlobalExtraction.Elt(i,j);
                }

                offset+=nbTreeVolCols;
                if(i<nbNodesRows)
                {
    			    vectorResult[i*nbCols+offset] = NbNodesPerSliceGlobalExtraction[i];
                }

                offset+=1;
                if(i<nbFxOptionRows)
                {
                    for ( j=0; j < nbFxOptionCols; ++j)
    			        vectorResult[i*nbCols+offset+j] = AnalyticFxOptionGlobalExtraction.Elt(i,j);
                }
            }
            /// add these additional lines 
		    /// to display blank lines
		    const int additionalLinesNb = 100;
		    bool fillWithBlank = true;
		    XL_writeNumMatrixSizeWithOptions( XL_result, vectorResult, nbRows, nbCols, " ARM_ERR: Could not set the num matrix", C_result,additionalLinesNb,fillWithBlank );
        }
        else
        {
			FreeCurCellErr ();
			XL_result.xltype = xltypeStr;
			XL_result.val.str = XL_StrC2StrPascal ( "No Vols" );
			XL_result.xltype |= xlbitDLLFree;
        }



	//	ARM_END();
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_ExtractTree3FVol" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


/////////////////////////////////////////
/// Functor to create a TriXBsModel
/////////////////////////////////////////
class TriXBSModelFunc : public ARMResultLong2LongFunc
{
public:
	TriXBSModelFunc(long Model1,
					long Model2,
					long Model3,
					long Idx1Idx2Corr,
					long Idx1Idx3Corr,
					long Idx2Idx3Corr)
					: C_Model1(Model1),
					C_Model2(Model2),
					C_Model3(Model3),
					C_Idx1Idx2Corr(Idx1Idx2Corr),
					C_Idx1Idx3Corr(Idx1Idx3Corr),
					C_Idx2Idx3Corr(Idx2Idx3Corr)
	{};

	long operator()( ARM_result& result, long objId )
	{
		return ARMLOCAL_TriXBSModel(C_Model1,
									C_Model2,
									C_Model3,
									C_Idx1Idx2Corr,
									C_Idx1Idx3Corr,
									C_Idx2Idx3Corr,
									result, 
									objId);
	};

private:
	
	long C_Model1;
    long C_Model2;
    long C_Model3;
    long C_Idx1Idx2Corr;
    long C_Idx1Idx3Corr;
    long C_Idx2Idx3Corr;
};


/////////////////////////////////////////
/// ReplicationConvexity Adjust Manager create common function
/////////////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_ARM_TRIXBSCommun (	LPXLOPER XL_Model1,
																LPXLOPER XL_Model2,
																LPXLOPER XL_Model3,
																LPXLOPER XL_Idx1Idx2Corr,
																LPXLOPER XL_Idx1Idx3Corr,
																LPXLOPER XL_Idx2Idx3Corr,
																bool PersistentInXL )
{
	ADD_LOG("Local_ARM_TRIXBSCommun ");
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

	CCString C_Model1;
	CCString C_Model2;
    CCString C_Model3;
    CCString C_Idx1Idx2Corr;
	CCString C_Idx1Idx3Corr;
	CCString C_Idx2Idx3Corr;
	
	long C_Model3Id;
	long C_Idx1Idx2CorrId;
	long C_Idx1Idx3CorrId;
	long C_Idx2Idx3CorrId;

	XL_readStrCell	(	 XL_Model1,		    C_Model1,					" ARM_ERR: First BSModel Id: string expected",		C_result);
	XL_readStrCell	(	 XL_Model2,		    C_Model2,					" ARM_ERR: Second BSModel Id: string expected",		C_result);
	XL_readStrCellWD(	 XL_Model3,			C_Model3,		"DEFAULT",	" ARM_ERR: Discount BSModel Id: string expected",	C_result);
	XL_readStrCellWD(	 XL_Idx1Idx2Corr,	C_Idx1Idx2Corr,	"DEFAULT",	" ARM_ERR: Idx1Idx2Corr Id: string expected",		C_result);
	XL_readStrCellWD(	 XL_Idx1Idx3Corr,	C_Idx1Idx3Corr,	"DEFAULT",	" ARM_ERR: Idx1Idx3Corr Id: string expected",		C_result);
	XL_readStrCellWD(	 XL_Idx2Idx3Corr,	C_Idx2Idx3Corr,	"DEFAULT",	" ARM_ERR: Idx2Idx3Corr Id: string expected",		C_result);
	
    if (C_Model3 == "DEFAULT")
	{
	   C_Model3Id = ARM_NULL_OBJECT;
	}
	else
	{
	   C_Model3Id = LocalGetNumObjectId (C_Model3);
	}

	if (C_Idx1Idx2Corr == "DEFAULT")
	{
	   C_Idx1Idx2CorrId = ARM_NULL_OBJECT;
	}
	else
	{
	   C_Idx1Idx2CorrId = LocalGetNumObjectId (C_Idx1Idx2Corr);
	}

	if (C_Idx1Idx3Corr == "DEFAULT")
	{
	   C_Idx1Idx3CorrId = ARM_NULL_OBJECT;
	}
	else
	{
	   C_Idx1Idx3CorrId = LocalGetNumObjectId (C_Idx1Idx3Corr);
	}

	if (C_Idx2Idx3Corr == "DEFAULT")
	{
	   C_Idx2Idx3CorrId = ARM_NULL_OBJECT;
	}
	else
	{
	   C_Idx2Idx3CorrId = LocalGetNumObjectId (C_Idx2Idx3Corr);
	}

	/// use the concept of Functor to transfer the knowledge of
	/// a function with a context
	TriXBSModelFunc ourFunc(LocalGetNumObjectId (C_Model1),
							LocalGetNumObjectId (C_Model2),
							C_Model3Id,
							C_Idx1Idx2CorrId,
							C_Idx1Idx3CorrId,
							C_Idx2Idx3CorrId);
	
	/// call the general function
	fillXL_Result( LOCAL_TRIXBSMOD_CLASS, ourFunc, C_result, XL_result, PersistentInXL );

	/// return the result as an LPXLOPER
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_ARM_TRIXBSCommun" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


//////////////////////////////////////////////////////
/// Addin to create a Replication Convexity Adjust Manager
//////////////////////////////////////////////////////

__declspec(dllexport) LPXLOPER WINAPI Local_ARM_TRIXBS (LPXLOPER XL_Model1,
														LPXLOPER XL_Model2,
														LPXLOPER XL_Model3,
														LPXLOPER XL_Idx1Idx2Corr,
														LPXLOPER XL_Idx1Idx3Corr,
														LPXLOPER XL_Idx2Idx3Corr)
{
	ADD_LOG("Local_ARM_TRIXBS ");
	bool PersistentInXL = true;
	return Local_ARM_TRIXBSCommun(	XL_Model1,
									XL_Model2,
									XL_Model3,
									XL_Idx1Idx2Corr,
									XL_Idx1Idx3Corr,
									XL_Idx2Idx3Corr,
									PersistentInXL);
}



//////////////////////////////////////////////////////
/// Addin to create a Replication Convexity Adjust Manager
/// Version for VBA
//////////////////////////////////////////////////////

__declspec(dllexport) LPXLOPER WINAPI Local_PXL_TRIXBS(	LPXLOPER XL_Model1,
														LPXLOPER XL_Model2,
														LPXLOPER XL_Model3,
														LPXLOPER XL_Idx1Idx2Corr,
														LPXLOPER XL_Idx1Idx3Corr,
														LPXLOPER XL_Idx2Idx3Corr)
{
	ADD_LOG("Local_PXL_TRIXBS");
	bool PersistentInXL = false;
	return Local_ARM_TRIXBSCommun(	XL_Model1,
									XL_Model2,
									XL_Model3,
									XL_Idx1Idx2Corr,
									XL_Idx1Idx3Corr,
									XL_Idx2Idx3Corr,
									PersistentInXL);
}