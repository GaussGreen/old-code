#pragma warning(disable : 4786)
#pragma warning(disable : 4005)


#include "ARM_xl_pf_local.h"
#include "ARM_xl_wrapper_local.h"

#include "ARM_local_interglob.h"
#include "ARM_local_interface.h"

#include <ARM\libarm_local\ARM_local_class.h>
#include <ARM\libarm_local\ARM_local_pf.h>

#include "ARM_xl_trycatch_local.h"

#include "util\tech_macro.h"

__declspec(dllexport) LPXLOPER WINAPI Local_PF (LPXLOPER XL_insts,
												LPXLOPER XL_coeffs,
												LPXLOPER XL_marketPrices,
												LPXLOPER XL_precisions)
{
	ADD_LOG("Local_PF ");
//    ARM_BEGIN();

    // return
    static XLOPER XL_result;
    ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();

    // C variable
    VECTOR<CCString> C_insts;
    VECTOR<double> C_coeffs;
    VECTOR<double> C_marketPrices;
    VECTOR<double> C_precisions;

    // error
    static int error;
    static char* reason = "";

    XL_readStrVector(XL_insts, C_insts," ARM_ERR: instruments: array of object expected",DOUBLE_TYPE,C_result);
    XL_readNumVector(XL_coeffs, C_coeffs," ARM_ERR: coefficients: array of numeric expected",C_result);
    XL_readNumVector(XL_marketPrices, C_marketPrices," ARM_ERR: market prices: array of numeric expected",C_result);

   	if((XL_precisions->xltype == xltypeMissing) || (XL_precisions->xltype == xltypeNil))
	{
    }
	else
	{
        XL_readNumVector(XL_precisions, C_precisions, " ARM_ERR: precisions: array of numeric expected",C_result);
    }

    VECTOR<long> C_instsId;
    long sz = C_insts.size ();

    for (int i = 0; i < sz; i++)
    {
        C_instsId.push_back (LocalGetNumObjectId (C_insts[i]));
    }
    
    long retCode;
    long objId;
    CCString prevClass;
    
    CCString curClass = LOCAL_PF_CLASS;
    CCString stringId = GetLastCurCellEnvValue();
    
    if (!stringId)
    {
        retCode = ARMLOCAL_PF (C_instsId,
							   C_coeffs,
							   C_marketPrices,
							   C_precisions,
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
        prevClass = LocalGetStringObjectClass(stringId);
        
        objId = LocalGetNumObjectId(stringId);
            
        if ( curClass == prevClass )
        {
            retCode = ARMLOCAL_PF(C_instsId,
								  C_coeffs,
								  C_marketPrices,
								  C_precisions,
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

            retCode = ARMLOCAL_PF (C_instsId,
								   C_coeffs,
								   C_marketPrices,
								   C_precisions,
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

//    ARM_END();
    }
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_PF" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}



__declspec(dllexport) LPXLOPER WINAPI Local_PXL_PF (LPXLOPER XL_insts,
													LPXLOPER XL_coeffs,
													LPXLOPER XL_marketPrices,
													LPXLOPER XL_precisions)
{
	ADD_LOG("Local_PXL_PF ");
//    ARM_BEGIN();

    // return
    static XLOPER XL_result;
    ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();

    // C variable
    VECTOR<CCString> C_insts;
    VECTOR<double> C_coeffs;
    VECTOR<double> C_marketPrices;
    VECTOR<double> C_precisions;

    // error
    static int error;
    static char* reason = "";

    XL_readStrVector(XL_insts, C_insts," ARM_ERR: instruments: array of object expected",DOUBLE_TYPE,C_result);
    XL_readNumVector(XL_coeffs, C_coeffs," ARM_ERR: coefficients: array of numeric expected",C_result);
    XL_readNumVector(XL_marketPrices, C_marketPrices," ARM_ERR: market prices: array of numeric expected",C_result);

   	if((XL_precisions->xltype == xltypeMissing) || (XL_precisions->xltype == xltypeNil))
	{
    }
	else
	{
        XL_readNumVector(XL_precisions, C_precisions, " ARM_ERR: precisions: array of numeric expected",C_result);
    }

    VECTOR<long> C_instsId;
    long sz = C_insts.size ();

    for (int i = 0; i < sz; i++)
    {
        C_instsId.push_back (LocalGetNumObjectId (C_insts[i]));
    }
    
    long retCode;
    long objId;
    
    CCString curClass = LOCAL_PF_CLASS;
    CCString stringId;
    
    retCode = ARMLOCAL_PF (C_instsId,
						   C_coeffs,
						   C_marketPrices,
						   C_precisions,
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

//    ARM_END();
    }
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_PXL_PF" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}

///----------------------------------------------
///----------------------------------------------
/// Portfolio incremental construction
///----------------------------------------------
///----------------------------------------------
////////////////////////////////////////////////
/// very rapid definition of ourlonglongFunctor
////////////////////////////////////////////////
class PFFILLFunc : public ARMResultLong2LongFunc
{
public:
	PFFILLFunc(
		const VECTOR<long >& assetsIdVec,
		const VECTOR<double >& weightsVec,
        const VECTOR<double >& mktpricesVec,
        const VECTOR<double >& vegasVec,
		long portfolioId)
    :
		C_assetsIdVec(assetsIdVec),
		C_weightsVec(weightsVec),
		C_mktpricesVec(mktpricesVec),
		C_vegasVec(vegasVec),
		C_portfolioId(portfolioId)
    {};

	long operator()( ARM_result& result, long objId ){
		return ARMLOCAL_PFFILL_Create(
            C_assetsIdVec,
            C_weightsVec,
            C_mktpricesVec,
            C_vegasVec,
            C_portfolioId,
            result,
            objId);			
	}
private:
        VECTOR<long>    C_assetsIdVec;
        VECTOR<double > C_weightsVec;
        VECTOR<double > C_mktpricesVec;
        VECTOR<double > C_vegasVec;
        long            C_portfolioId;
};


/////////////////////////////////////////////////////////////
/// central function that does the creation of the XL function
/////////////////////////////////////////////////////////////
LPXLOPER Local_PFFILLCommon(LPXLOPER XL_Assets,
	LPXLOPER XL_Weights,
	LPXLOPER XL_MktPrices,
	LPXLOPER XL_VegaPrecisions,
    LPXLOPER XL_PortfolioId,
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

        VECTOR<CCString> C_assetsStr;
        XL_readStrVector(XL_Assets, C_assetsStr," ARM_ERR: Assets: array of object expected",DOUBLE_TYPE,C_result);

        long size = C_assetsStr.size ();
        VECTOR<long> C_assetsId(size);

        for (int i = 0; i < size; i++)
            C_assetsId[i] = LocalGetNumObjectId (C_assetsStr[i]);

        VECTOR<double> C_weights;
        XL_readNumVector(XL_Weights, C_weights," ARM_ERR: coefficients: array of numeric expected",C_result);

        VECTOR<double> C_mktPrices;
        XL_readNumVector(XL_MktPrices, C_mktPrices," ARM_ERR: market prices: array of numeric expected",C_result);

        VECTOR<double> defaultBound;
        VECTOR<double> C_vegaPrecisions;
        XL_readNumVectorWD(XL_VegaPrecisions,C_vegaPrecisions,defaultBound," ARM_ERR: lower bound: array of numeric expected",C_result);

        CCString C_pfStrId;
		XL_readStrCellWD( XL_PortfolioId, C_pfStrId,"NULL"," ARM_ERR: Portfolio Id: Object expected",		C_result);
        long C_PflId =  (C_pfStrId == "NULL") ? ARM_NULL_OBJECT : LocalGetNumObjectId(C_pfStrId);
		
		/// use the concept of Functor to transfer the knowledge of
		/// a function with a context
		PFFILLFunc ourFunc(
			C_assetsId,
			C_weights,
			C_mktPrices,
			C_vegaPrecisions,
			C_PflId);

		/// call the general function
		fillXL_Result( LOCAL_PF_CLASS, ourFunc, C_result, XL_result, PersistentInXL );
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_PFFILLCommon" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}



///////////////////////////////////
/// version that takes into account 
/// previous creation of object
///////////////////////////////////

__declspec(dllexport) LPXLOPER WINAPI Local_PF_FILL (LPXLOPER XL_Assets,
												LPXLOPER XL_Weights,
												LPXLOPER XL_MktPrices,
												LPXLOPER XL_VegaPrecisions,
                                                LPXLOPER XL_PortfolioId)
{
	ADD_LOG("Local_PF_FILL ");
	bool PersistentInXL = true;
	return Local_PFFILLCommon(
        XL_Assets,
		XL_Weights,
		XL_MktPrices,
		XL_VegaPrecisions,
		XL_PortfolioId,
		PersistentInXL );
}

///////////////////////////////////
/// version for non persistentInXl
///////////////////////////////////

__declspec(dllexport) LPXLOPER WINAPI Local_PXL_PF_FILL (LPXLOPER XL_Assets,
												LPXLOPER XL_Weights,
												LPXLOPER XL_MktPrices,
												LPXLOPER XL_VegaPrecisions,
                                                LPXLOPER XL_PortfolioId)
{
	ADD_LOG("Local_PXL_PF_FILL ");
	bool PersistentInXL = false;
	return Local_PFFILLCommon(
        XL_Assets,
		XL_Weights,
		XL_MktPrices,
		XL_VegaPrecisions,
		XL_PortfolioId,
		PersistentInXL );
}


///----------------------------------------------
///----------------------------------------------
/// Portfolio merging
///----------------------------------------------
///----------------------------------------------
class PFMergeFunc : public ARMResultLong2LongFunc
{
public:
	PFMergeFunc(const VECTOR<long> portfolioIds)
    : C_portfolioIds(portfolioIds) {};

	long operator()( ARM_result& result, long objId ){
		return ARMLOCAL_PF_Merge(C_portfolioIds,result,objId);}
private:
        VECTOR<long> C_portfolioIds;
};


/////////////////////////////////////////////////////////////
/// central function that does the creation of the XL function
/////////////////////////////////////////////////////////////
LPXLOPER Local_PFMergeCommon(LPXLOPER XL_PortfolioIds,bool PersistentInXL)
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

        VECTOR<CCString> C_PfStrIds;
		XL_readStrVector(XL_PortfolioIds,C_PfStrIds," ARM_ERR: Portfolio Id vector: array of object expected",DOUBLE_TYPE,C_result);
		VECTOR<long> C_PfIds(C_PfStrIds.size());
		for(size_t i=0;i<C_PfStrIds.size();++i)
			C_PfIds[i] = (C_PfStrIds[i] == "NULL") ? ARM_NULL_OBJECT : LocalGetNumObjectId(C_PfStrIds[i]);
		
		/// use the concept of Functor to transfer the knowledge of
		/// a function with a context
		PFMergeFunc ourFunc(C_PfIds);

		/// call the general function
		fillXL_Result( LOCAL_PF_CLASS, ourFunc, C_result, XL_result, PersistentInXL );
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_PFMergeCommon" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}



///////////////////////////////////
/// version that takes into account 
/// previous creation of object
///////////////////////////////////

__declspec(dllexport) LPXLOPER WINAPI Local_PF_Merge (LPXLOPER XL_PortfolioIds)
{
	ADD_LOG("Local_PF_Merge ");
	bool PersistentInXL = true;
	return Local_PFMergeCommon(XL_PortfolioIds,PersistentInXL);
}

///////////////////////////////////
/// version for non persistentInXl
///////////////////////////////////

__declspec(dllexport) LPXLOPER WINAPI Local_PXL_PF_Merge (LPXLOPER XL_PortfolioIds)
{
	ADD_LOG("Local_PXL_PF_Merge ");
	bool PersistentInXL = false;
	return Local_PFMergeCommon(XL_PortfolioIds,PersistentInXL);
}








__declspec(dllexport) LPXLOPER WINAPI Local_PFGYCSIGVARFIT (LPXLOPER XL_pf,
															LPXLOPER XL_curve,
															LPXLOPER XL_matCurve,
															LPXLOPER XL_in_min_meanRev,
															LPXLOPER XL_in_max_meanRev,
															LPXLOPER XL_in_min_vol,
															LPXLOPER XL_in_max_vol,
															LPXLOPER XL_in_precision_meanRev,
															LPXLOPER XL_in_precision_vol)
{
	ADD_LOG("Local_PFGYCSIGVARFIT ");
//    ARM_BEGIN();

    // return
    static XLOPER XL_result;
    ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();

    // C variable
    CCString C_pf;

    CCString C_curve;

    VECTOR<double> C_matCurve;

    double C_in_min_meanRev;
    double C_in_min_meanRev_default = -0.01;

    double C_in_max_meanRev;
    double C_in_max_meanRev_default = 0.2;

    double C_in_min_vol;
    double C_in_min_vol_default = 0.01;

    double C_in_max_vol;
    double C_in_max_vol_default = 2.0;

    double C_in_precision_meanRev;
    double C_in_precision_meanRev_default = 0.01;

    double C_in_precision_vol;
    double C_in_precision_vol_default = 0.0001;

    long C_nbMaxIter = 1000;
        
    // error
    static int error;
    static char* reason = "";

    XL_readStrCell(XL_pf,C_pf," ARM_ERR: portfolio id: object expected",C_result);
    XL_readStrCell(XL_curve,C_curve," ARM_ERR: curve id: object expected",C_result);
    XL_readNumVector(XL_matCurve,C_matCurve," ARM_ERR: maturity curve: array of dates expected",C_result);
    XL_readNumCellWD(XL_in_min_meanRev,C_in_min_meanRev,C_in_min_meanRev_default," ARM_ERR: in minimum mean reversion: numeric expected",C_result);
    XL_readNumCellWD(XL_in_max_meanRev,C_in_max_meanRev,C_in_max_meanRev_default," ARM_ERR: in maximum mean reversion: numeric expected",C_result);
    XL_readNumCellWD(XL_in_min_vol,C_in_min_vol,C_in_min_vol_default," ARM_ERR: in min volatility: numeric expected",C_result);
    XL_readNumCellWD(XL_in_max_vol,C_in_max_vol,C_in_max_vol_default," ARM_ERR: in max volatility: numeric expected",C_result);
    XL_readNumCellWD(XL_in_precision_meanRev,C_in_precision_meanRev,C_in_precision_meanRev_default," ARM_ERR: in precision mean reversion: numeric expected",C_result);
    XL_readNumCellWD(XL_in_precision_vol,C_in_precision_vol,C_in_precision_vol_default," ARM_ERR: in precision volatility: numeric expected",C_result);

    long retCode;
    long objId;
    CCString prevClass;
    
    CCString curClass = LOCAL_HWSIGVAR_ANALYTIC_CLASS;
    CCString stringId = GetLastCurCellEnvValue ();
    
    if(!stringId)
    {
        retCode = ARMLOCAL_PFGYCSIGVARFIT (LocalGetNumObjectId (C_pf),
                                      LocalGetNumObjectId (C_curve),
                                      C_matCurve,
                                      C_in_min_meanRev,
                                      C_in_max_meanRev,
                                      C_in_min_vol,
                                      C_in_max_vol,
                                      C_in_precision_meanRev,
                                      C_in_precision_vol,
                                      C_nbMaxIter,
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
            retCode = ARMLOCAL_PFGYCSIGVARFIT (LocalGetNumObjectId (C_pf),
                                          LocalGetNumObjectId (C_curve),
                                          C_matCurve,
                                          C_in_min_meanRev,
                                          C_in_max_meanRev,
                                          C_in_min_vol,
                                          C_in_max_vol,
                                          C_in_precision_meanRev,
                                          C_in_precision_vol,
                                          C_nbMaxIter,
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
            retCode = ARMLOCAL_PFGYCSIGVARFIT (LocalGetNumObjectId (C_pf),
                                          LocalGetNumObjectId (C_curve),
                                          C_matCurve,
                                          C_in_min_meanRev,
                                          C_in_max_meanRev,
                                          C_in_min_vol,
                                          C_in_max_vol,
                                          C_in_precision_meanRev,
                                          C_in_precision_vol,
                                          C_nbMaxIter,
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

//    ARM_END();
    }
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_PFGYCSIGVARFIT" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


__declspec(dllexport) LPXLOPER WINAPI Local_PXL_PFGYCSIGVARFIT (LPXLOPER XL_pf,
																LPXLOPER XL_curve,
																LPXLOPER XL_matCurve,
																LPXLOPER XL_in_min_meanRev,
																LPXLOPER XL_in_max_meanRev,
																LPXLOPER XL_in_min_vol,
																LPXLOPER XL_in_max_vol,
																LPXLOPER XL_in_precision_meanRev,
																LPXLOPER XL_in_precision_vol)
{
	ADD_LOG("Local_PXL_PFGYCSIGVARFIT ");
//    ARM_BEGIN();

    // return
    static XLOPER XL_result;
    ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();

    // C variable
	CCString C_pf;

	CCString C_curve;

	VECTOR<double> C_matCurve;

	double C_in_min_meanRev;
	double C_in_min_meanRev_default = -0.01;

	double C_in_max_meanRev;
	double C_in_max_meanRev_default = 0.2;

	double C_in_min_vol;
	double C_in_min_vol_default = 0.01;

	double C_in_max_vol;
	double C_in_max_vol_default = 2.0;

	double C_in_precision_meanRev;
	double C_in_precision_meanRev_default = 0.01;

	double C_in_precision_vol;
	double C_in_precision_vol_default = 0.0001;

	long C_nbMaxIter = 1000;
        
    // error
    static int error;
    static char* reason = "";

	XL_readStrCell(XL_pf,C_pf," ARM_ERR: portfolio id: object expected",C_result);
	XL_readStrCell(XL_curve,C_curve," ARM_ERR: curve id: object expected",C_result);
	XL_readNumVector(XL_matCurve,C_matCurve," ARM_ERR: maturity curve: array of dates expected",C_result);
	XL_readNumCellWD(XL_in_min_meanRev,C_in_min_meanRev,C_in_min_meanRev_default," ARM_ERR: in minimum mean reversion: numeric expected",C_result);
	XL_readNumCellWD(XL_in_max_meanRev,C_in_max_meanRev,C_in_max_meanRev_default," ARM_ERR: in maximum mean reversion: numeric expected",C_result);
	XL_readNumCellWD(XL_in_min_vol,C_in_min_vol,C_in_min_vol_default," ARM_ERR: in min volatility: numeric expected",C_result);
	XL_readNumCellWD(XL_in_max_vol,C_in_max_vol,C_in_max_vol_default," ARM_ERR: in max volatility: numeric expected",C_result);
	XL_readNumCellWD(XL_in_precision_meanRev,C_in_precision_meanRev,C_in_precision_meanRev_default," ARM_ERR: in precision mean reversion: numeric expected",C_result);
	XL_readNumCellWD(XL_in_precision_vol,C_in_precision_vol,C_in_precision_vol_default," ARM_ERR: in precision volatility: numeric expected",C_result);

	long retCode;
	long objId;

	CCString curClass = LOCAL_HWSIGVAR_ANALYTIC_CLASS;
	CCString stringId;
    
    retCode = ARMLOCAL_PFGYCSIGVARFIT (LocalGetNumObjectId (C_pf),
                                  LocalGetNumObjectId (C_curve),
                                  C_matCurve,
                                  C_in_min_meanRev,
                                  C_in_max_meanRev,
                                  C_in_min_vol,
                                  C_in_max_vol,
                                  C_in_precision_meanRev,
                                  C_in_precision_vol,
                                  C_nbMaxIter,
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

//    ARM_END();
    }
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_PXL_PFGYCSIGVARFIT" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}



__declspec(dllexport) LPXLOPER WINAPI Local_PFLOGDECVOLFIT(LPXLOPER XL_pfId,
														   LPXLOPER XL_curveId,
														   LPXLOPER XL_proba,
														   LPXLOPER XL_accuracy,
														   LPXLOPER XL_shapeType,
														   LPXLOPER XL_decay,
														   LPXLOPER XL_slope,
														   LPXLOPER XL_asymptote,
														   LPXLOPER XL_matCurve,
														   LPXLOPER XL_volinit_vect,
														   LPXLOPER XL_coeff_vect)
{
	ADD_LOG("Local_PFLOGDECVOLFIT");
	//    ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();

	// C variable
	CCString C_pfId;

	CCString C_curveId;

	VECTOR<double> C_proba;
	VECTOR<double> C_matCurve;

	VECTOR<double> C_volinit_vect;
	VECTOR<double> C_coeff_vect;

	double C_accuracy;
	double C_accuracy_default  = 0.00000000001;
	double C_shapeType;
	double C_shape_default = 0;
	double C_decay, C_slope, C_asymptote;
	double C_value_default = 0.0;

	// error
	static int error;
	static char* reason = "";

	XL_readStrCell(XL_pfId, C_pfId," ARM_ERR: portfolio id: object expected",C_result);
	XL_readStrCell(XL_curveId, C_curveId," ARM_ERR: curve id: object expected",C_result);
	XL_readNumVector(XL_proba, C_proba," ARM_ERR: Proba. Vector: array of numeric values expected",C_result);
	XL_readNumCellWD(XL_accuracy, C_accuracy,C_accuracy_default," ARM_ERR: the accuracy: numeric expected",C_result);
	XL_readNumCellWD(XL_shapeType, C_shapeType,C_shape_default," ARM_ERR: the Shape_Type: integer expected",C_result);
	XL_readNumCellWD(XL_decay, C_decay,C_value_default," ARM_ERR: the decay: numeric expected",C_result);
	XL_readNumCellWD(XL_slope, C_slope,C_value_default," ARM_ERR: the slope: numeric expected",C_result);
	XL_readNumCellWD(XL_asymptote, C_asymptote,C_value_default," ARM_ERR: the asymptote: numeric expected",C_result);

	if (( XL_matCurve->xltype == xltypeMissing )
		||
		( XL_matCurve->xltype == xltypeNil )
	   )
	{
		/*--- The size of C_matCurve must be == 0 ----*/
	}
	else
	{
		XL_readNumVector(XL_matCurve, C_matCurve," ARM_ERR: maturity curve: array of dates expected",C_result);
	}

	if (( XL_volinit_vect->xltype == xltypeMissing )
		||
		( XL_volinit_vect->xltype == xltypeNil )
	   )
	{
		if (C_shapeType)
		{
			for (int k = 0; k < C_proba.size (); k++)
			{
				C_volinit_vect.push_back(0.0);
			}      
		}
		else
		{
			for (int k = 0; k < C_matCurve.size (); k++)
			{
				C_volinit_vect.push_back(0.0);
			}      
		}
	}
	else
	{
		XL_readNumVector(XL_volinit_vect, C_volinit_vect," ARM_ERR: penalties VECTOR: array of numeric expected",C_result);
	}

	if (( XL_coeff_vect->xltype == xltypeMissing )
		||
		( XL_coeff_vect->xltype == xltypeNil )
	   )
	{
		for (int h = 0; h < C_volinit_vect.size (); h++)
		{
			C_coeff_vect.push_back(0.0);
		}  
	}
	else
	{
		XL_readNumVector(XL_coeff_vect, C_coeff_vect," ARM_ERR: coefficients VECTOR: array of numeric expected",C_result);
	}

	long retCode;
	long objId;
	CCString prevClass;

	CCString curClass = LOCAL_LOGDECANA_CLASS;
	CCString stringId = GetLastCurCellEnvValue ();
    
    if(!stringId)
    {
		retCode = ARMLOCAL_PFLOGDECVOLFIT(LocalGetNumObjectId(C_pfId),
                                     LocalGetNumObjectId(C_curveId),
                                     C_proba,
                                     C_accuracy,
                                     (long) C_shapeType,
                                     C_decay,
                                     C_slope, 
                                     C_asymptote,
                                     C_matCurve,
                                     C_volinit_vect,
                                     C_coeff_vect,
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
			retCode = ARMLOCAL_PFLOGDECVOLFIT(LocalGetNumObjectId(C_pfId),
                                            LocalGetNumObjectId(C_curveId),
                                            C_proba,
                                            C_accuracy,
                                            (long) C_shapeType,
                                            C_decay,
                                            C_slope, 
                                            C_asymptote,
                                            C_matCurve,
                                            C_volinit_vect,
                                            C_coeff_vect,
                                            C_result,
                                            objId);
		}
		else
		{
            FreeCurCellContent ();

            retCode = ARMLOCAL_PFLOGDECVOLFIT(LocalGetNumObjectId(C_pfId),
                                         LocalGetNumObjectId(C_curveId),
                                         C_proba,
                                         C_accuracy,
                                         (long) C_shapeType,
                                         C_decay,
                                         C_slope, 
                                         C_asymptote,
                                         C_matCurve,
                                         C_volinit_vect,
                                         C_coeff_vect,
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

//    ARM_END();
    }
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_PFLOGDECVOLFIT" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}



__declspec(dllexport) LPXLOPER WINAPI Local_PXL_PFLOGDECVOLFIT(LPXLOPER XL_pfId,
															   LPXLOPER XL_curveId,
															   LPXLOPER XL_proba,
															   LPXLOPER XL_accuracy,
															   LPXLOPER XL_shapeType,
															   LPXLOPER XL_decay,
															   LPXLOPER XL_slope,
															   LPXLOPER XL_asymptote,
															   LPXLOPER XL_matCurve,
															   LPXLOPER XL_volinit_vect,
															   LPXLOPER XL_coeff_vect)
{
	ADD_LOG("Local_PXL_PFLOGDECVOLFIT");
	//    ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();

	// C variable
	CCString C_pfId;

	CCString C_curveId;

	VECTOR<double> C_proba;
	VECTOR<double> C_matCurve;

	VECTOR<double> C_volinit_vect;
	VECTOR<double> C_coeff_vect;

	double C_accuracy;
	double C_accuracy_default  = 0.00000000001;
	double C_shapeType;
	double C_shape_default = 0;
	double C_decay, C_slope, C_asymptote;
	double C_value_default = 0.0;

	// error
	static int error;
	static char* reason = "";

	XL_readStrCell(XL_pfId, C_pfId," ARM_ERR: portfolio id: object expected",C_result);
	XL_readStrCell(XL_curveId, C_curveId," ARM_ERR: curve id: object expected",C_result);
	XL_readNumVector(XL_proba, C_proba," ARM_ERR: Proba. Vector: array of numeric values expected",C_result);
	XL_readNumCellWD(XL_accuracy, C_accuracy,C_accuracy_default," ARM_ERR: the accuracy: numeric expected",C_result);
	XL_readNumCellWD(XL_shapeType, C_shapeType,C_shape_default," ARM_ERR: the Shape_Type: integer expected",C_result);
	XL_readNumCellWD(XL_decay, C_decay,C_value_default," ARM_ERR: the decay: numeric expected",C_result);
	XL_readNumCellWD(XL_slope, C_slope,C_value_default," ARM_ERR: the slope: numeric expected",C_result);
	XL_readNumCellWD(XL_asymptote, C_asymptote,C_value_default," ARM_ERR: the asymptote: numeric expected",C_result);

	if (( XL_matCurve->xltype == xltypeMissing )
		||
		( XL_matCurve->xltype == xltypeNil )
	   )
	{
		/*--- The size of C_matCurve must be == 0 ----*/
	}
	else
	{
		XL_readNumVector(XL_matCurve, C_matCurve," ARM_ERR: maturity curve: array of dates expected",C_result);
	}

	if (( XL_volinit_vect->xltype == xltypeMissing )
		||
		( XL_volinit_vect->xltype == xltypeNil )
	   )
	{
		if (C_shapeType)
		{
			for (int k = 0; k < C_proba.size (); k++)
			{
				C_volinit_vect.push_back(0.0);
			}      
		}
		else
		{
			for (int k = 0; k < C_matCurve.size (); k++)
			{
				C_volinit_vect.push_back(0.0);
			}      
		}
	}
	else
	{
		XL_readNumVector(XL_volinit_vect, C_volinit_vect," ARM_ERR: penalties VECTOR: array of numeric expected",C_result);
	}

	if (( XL_coeff_vect->xltype == xltypeMissing )
		||
		( XL_coeff_vect->xltype == xltypeNil )
	   )
	{
		for (int h = 0; h < C_volinit_vect.size (); h++)
		{
			C_coeff_vect.push_back(0.0);
		}  
	}
	else
	{
		XL_readNumVector(XL_coeff_vect, C_coeff_vect," ARM_ERR: coefficients VECTOR: array of numeric expected",C_result);
	}

	long retCode;
	long objId;

	CCString curClass = LOCAL_LOGDECANA_CLASS;
	CCString stringId;
    
	retCode = ARMLOCAL_PFLOGDECVOLFIT(LocalGetNumObjectId(C_pfId),
                                 LocalGetNumObjectId(C_curveId),
                                 C_proba,
                                 C_accuracy,
                                 (long) C_shapeType,
                                 C_decay,
                                 C_slope, 
                                 C_asymptote,
                                 C_matCurve,
                                 C_volinit_vect,
                                 C_coeff_vect,
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

//    ARM_END();
    }
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_PXL_PFLOGDECVOLFIT" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


__declspec(dllexport) LPXLOPER WINAPI Local_PFGYCSIGVARPENALFIT(LPXLOPER XL_pfId,
																LPXLOPER XL_curveId,
																LPXLOPER XL_matCurve,
																LPXLOPER XL_start_meanRev,
																LPXLOPER XL_start_vol,
																LPXLOPER XL_penal_vect,
																LPXLOPER XL_coeff_vect,
																LPXLOPER XL_accuracy)
{
	ADD_LOG("Local_PFGYCSIGVARPENALFIT");
//    ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();

	// C variable
	CCString C_pfId;

	CCString C_curveId;

	VECTOR<double> C_matCurve;
	VECTOR<double> C_start_vol;

	VECTOR<double> C_penal_vect;
	VECTOR<double> C_coeff_vect;

	double C_start_meanRev;
	double C_start_meanRev_default  = 0.1;

	double C_start_vol_default      = 1.0;

	double C_coeff_vect_default     = 0.0;

	double C_accuracy;
	double C_accuracy_default  = 0.1;
    
	// error
	static int error;
	static char* reason = "";

	XL_readStrCell(XL_pfId, C_pfId," ARM_ERR: portfolio id: object expected",C_result);
	XL_readStrCell(XL_curveId, C_curveId," ARM_ERR: curve id: object expected",C_result);
	XL_readNumVector(XL_matCurve, C_matCurve," ARM_ERR: maturity curve: array of dates expected",C_result);
	XL_readNumCellWD(XL_start_meanRev, C_start_meanRev,C_start_meanRev_default," ARM_ERR: start mean reversion: numeric expected", C_result);
    
	if (( XL_start_vol->xltype == xltypeMissing )
		||
		( XL_start_vol->xltype == xltypeNil )
	   )
	{
		C_start_vol.push_back(C_start_vol_default);
	}
	else
	{
	   XL_readNumVector(XL_start_vol, C_start_vol," ARM_ERR: start volatility: array of numeric expected",C_result);
	}

	if (( XL_penal_vect->xltype == xltypeMissing )
		||
		( XL_penal_vect->xltype == xltypeNil )
	   )
	{
		C_penal_vect.push_back(C_start_meanRev_default);
		C_penal_vect.push_back(C_start_vol_default);
	}
	else
	{
		XL_readNumVector(XL_penal_vect, C_penal_vect," ARM_ERR: penalties VECTOR: array of numeric expected",C_result);
	}

	if (( XL_coeff_vect->xltype == xltypeMissing )
		||
		( XL_coeff_vect->xltype == xltypeNil )
	   )
	{
		C_coeff_vect.push_back(C_coeff_vect_default);
		C_coeff_vect.push_back(C_coeff_vect_default);
	}
	else
	{
		XL_readNumVector(XL_coeff_vect, C_coeff_vect," ARM_ERR: coefficients VECTOR: array of numeric expected",C_result);
	}

	XL_readNumCellWD(XL_accuracy, C_accuracy,C_accuracy_default," ARM_ERR: the accuracy: numeric expected",C_result);
    
	if ( C_start_vol.size() == 1 )
	{
		for (int i = 1; i < C_matCurve.size(); i++)
		{
			C_start_vol.push_back(C_start_vol[0]);
		}
	}

    long retCode;
    long objId;
    CCString prevClass;
    
    CCString curClass = LOCAL_HWSIGVAR_ANALYTIC_CLASS;
    CCString stringId = GetLastCurCellEnvValue ();

	if(!stringId)
	{
		retCode = ARMLOCAL_PFGYCSIGVARPENALFIT(LocalGetNumObjectId(C_pfId),
										  LocalGetNumObjectId(C_curveId),
										  C_matCurve,
										  C_start_meanRev,
										  C_start_vol,
										  C_penal_vect,
										  C_coeff_vect,
										  C_accuracy,
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
			retCode = ARMLOCAL_PFGYCSIGVARPENALFIT(LocalGetNumObjectId(C_pfId),
											 LocalGetNumObjectId(C_curveId),
											 C_matCurve,
											 C_start_meanRev,
											 C_start_vol,
											 C_penal_vect,
											 C_coeff_vect,
											 C_accuracy,
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

			retCode = ARMLOCAL_PFGYCSIGVARPENALFIT(LocalGetNumObjectId(C_pfId),
											  LocalGetNumObjectId(C_curveId),
											  C_matCurve,
											  C_start_meanRev,
											  C_start_vol,
											  C_penal_vect,
											  C_coeff_vect,
											  C_accuracy,
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

//    ARM_END();
    }
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_PFGYCSIGVARPENALFIT" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


__declspec(dllexport) LPXLOPER WINAPI Local_PXL_PFGYCSIGVARPENALFIT(LPXLOPER XL_pfId,
																	LPXLOPER XL_curveId,
																	LPXLOPER XL_matCurve,
																	LPXLOPER XL_start_meanRev,
																	LPXLOPER XL_start_vol,
																	LPXLOPER XL_penal_vect,
																	LPXLOPER XL_coeff_vect,
																	LPXLOPER XL_accuracy)
{
	ADD_LOG("Local_PXL_PFGYCSIGVARPENALFIT");
//    ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();

	// C variable
	CCString C_pfId;

	CCString C_curveId;

	VECTOR<double> C_matCurve;
	VECTOR<double> C_start_vol;

	VECTOR<double> C_penal_vect;
	VECTOR<double> C_coeff_vect;

	double C_start_meanRev;
	double C_start_meanRev_default  = 0.1;

	double C_start_vol_default      = 1.0;

	double C_coeff_vect_default     = 0.0;

	double C_accuracy;
	double C_accuracy_default  = 0.1;
    
	// error
	static int error;
	static char* reason = "";

	XL_readStrCell(XL_pfId, C_pfId," ARM_ERR: portfolio id: object expected",C_result);
	XL_readStrCell(XL_curveId, C_curveId," ARM_ERR: curve id: object expected",C_result);
	XL_readNumVector(XL_matCurve, C_matCurve," ARM_ERR: maturity curve: array of dates expected",C_result);
	XL_readNumCellWD(XL_start_meanRev, C_start_meanRev,C_start_meanRev_default," ARM_ERR: start mean reversion: numeric expected", C_result);
    
	if (( XL_start_vol->xltype == xltypeMissing )
		||
		( XL_start_vol->xltype == xltypeNil )
	   )
	{
		C_start_vol.push_back(C_start_vol_default);
	}
	else
	{
	   XL_readNumVector(XL_start_vol, C_start_vol," ARM_ERR: start volatility: array of numeric expected",C_result);
	}

	if (( XL_penal_vect->xltype == xltypeMissing )
		||
		( XL_penal_vect->xltype == xltypeNil )
	   )
	{
		C_penal_vect.push_back(C_start_meanRev_default);
		C_penal_vect.push_back(C_start_vol_default);
	}
	else
	{
		XL_readNumVector(XL_penal_vect, C_penal_vect," ARM_ERR: penalties VECTOR: array of numeric expected",C_result);
	}

	if (( XL_coeff_vect->xltype == xltypeMissing )
		||
		( XL_coeff_vect->xltype == xltypeNil )
	   )
	{
		C_coeff_vect.push_back(C_coeff_vect_default);
		C_coeff_vect.push_back(C_coeff_vect_default);
	}
	else
	{
		XL_readNumVector(XL_coeff_vect, C_coeff_vect," ARM_ERR: coefficients VECTOR: array of numeric expected",C_result);
	}

	XL_readNumCellWD(XL_accuracy, C_accuracy,C_accuracy_default," ARM_ERR: the accuracy: numeric expected",C_result);
    
	if ( C_start_vol.size() == 1 )
	{
		for (int i = 1; i < C_matCurve.size(); i++)
		{
			C_start_vol.push_back(C_start_vol[0]);
		}
	}

    long retCode;
    long objId;
    
    CCString curClass = LOCAL_HWSIGVAR_ANALYTIC_CLASS;
    CCString stringId;

	retCode = ARMLOCAL_PFGYCSIGVARPENALFIT(LocalGetNumObjectId(C_pfId),
									  LocalGetNumObjectId(C_curveId),
									  C_matCurve,
									  C_start_meanRev,
									  C_start_vol,
									  C_penal_vect,
									  C_coeff_vect,
									  C_accuracy,
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

//    ARM_END();
   }
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_PXL_PFGYCSIGVARPENALFIT" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;

}


__declspec(dllexport) LPXLOPER WINAPI Local_PFINSTLOGDECVOLFIT(LPXLOPER XL_pfId,
															   LPXLOPER XL_secId,
															   LPXLOPER XL_curveId,
															   LPXLOPER XL_proba,
															   LPXLOPER XL_UsePFResetDates,
															   LPXLOPER XL_accuracy,
															   LPXLOPER XL_shapeType,
															   LPXLOPER XL_decay,
															   LPXLOPER XL_slope,
															   LPXLOPER XL_asymptote,
															   LPXLOPER XL_VolBSId,
															   LPXLOPER XL_matCurve,
															   LPXLOPER XL_volinit_vect,
															   LPXLOPER XL_coeff_vect)
{
	ADD_LOG("Local_PFINSTLOGDECVOLFIT");
//    ARM_BEGIN();

    // return
    static XLOPER XL_result;
    ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();

    // C variable
	CCString C_secId;
	CCString C_curveId;
	CCString C_pfId;

	CCString C_VolBSId;

	VECTOR<double> C_proba;
	VECTOR<double> C_matCurve;

	VECTOR<double> C_volinit_vect;
	VECTOR<double> C_coeff_vect;

	double C_accuracy;
	double C_accuracy_default  = 0.00000000001;
	double C_shapeType;
	double C_shape_default = 0;
	double C_decay, C_slope, C_asymptote;
	double C_value_default = 0.0;

	double C_UsePFResetDates;
	double C_UsePFResetDates_default = 0.0;

	// error
	static int error;
	static char* reason = "";

	XL_readStrCell(XL_pfId, C_pfId," ARM_ERR: portfolio id: object expected", C_result);
	XL_readStrCell(XL_secId, C_secId," ARM_ERR: security id: object expected", C_result);
	XL_readStrCell(XL_curveId, C_curveId," ARM_ERR: curve id: object expected", C_result); 
	XL_readNumVector(XL_proba, C_proba, " ARM_ERR: Proba. Vector: array of numeric values expected", C_result);
	XL_readNumCellWD(XL_UsePFResetDates, C_UsePFResetDates, C_UsePFResetDates_default, " ARM_ERR: the UsePFResetDates: numeric expected",C_result);
	XL_readNumCellWD(XL_accuracy, C_accuracy, C_accuracy_default, " ARM_ERR: the accuracy: numeric expected",C_result);
	XL_readNumCellWD(XL_shapeType, C_shapeType, C_shape_default, " ARM_ERR: the Shape_Type: integer expected",C_result);
	XL_readNumCellWD(XL_decay, C_decay, C_value_default, " ARM_ERR: the decay: numeric expected",C_result);
	XL_readNumCellWD(XL_slope, C_slope, C_value_default, " ARM_ERR: the slope: numeric expected",C_result);
	XL_readNumCellWD(XL_asymptote, C_asymptote, C_value_default, " ARM_ERR: the asymptote: numeric expected",C_result);
	XL_readStrCellWD(XL_VolBSId, C_VolBSId,"DEFVOLBS"," ARM_ERR: BS Vol Curve ATM expected",C_result);

	long VolBSId = -1;
	if (!( C_VolBSId == "DEFVOLBS"))
	{
		VolBSId = LocalGetNumObjectId(C_VolBSId);
	}

    if (( XL_matCurve->xltype == xltypeMissing )
        ||
        ( XL_matCurve->xltype == xltypeNil )
       )
    {
       /*--- The size of C_matCurve must be == 0 ----*/
    }
    else
    {
        XL_readNumVector(XL_matCurve, C_matCurve, " ARM_ERR: maturity curve: array of dates expected", C_result);
    }

	if (( XL_volinit_vect->xltype == xltypeMissing )
		||
		( XL_volinit_vect->xltype == xltypeNil )
		)
	{
		if (C_shapeType)
		{
			for (int k = 0; k < C_proba.size (); k++)
				C_volinit_vect.push_back(0.0);
		}
		else
		{
			for (int k = 0; k < C_matCurve.size (); k++)
				C_volinit_vect.push_back(0.0);
		}
	}
	else
		XL_readNumVector(XL_volinit_vect, C_volinit_vect, " ARM_ERR: penalties VECTOR: array of numeric expected", C_result);

	if (( XL_coeff_vect->xltype == xltypeMissing )
		||
		( XL_coeff_vect->xltype == xltypeNil )
		)
	{
		for (int h = 0; h < C_volinit_vect.size (); h++)
			C_coeff_vect.push_back(0.0);
	}
	else
		XL_readNumVector(XL_coeff_vect, C_coeff_vect, " ARM_ERR: coefficients VECTOR: array of numeric expected", C_result);

	long retCode;
	long objId;
	CCString prevClass;
    
	CCString curClass = LOCAL_LOGDECANA_CLASS;
	CCString stringId = GetLastCurCellEnvValue ();

	if(!stringId)
	{
		retCode = ARMLOCAL_PFINSTLOGDECVOLFIT(LocalGetNumObjectId(C_pfId),
											  LocalGetNumObjectId(C_secId),
											  LocalGetNumObjectId(C_curveId),
											  C_proba,
											  (long) C_UsePFResetDates,
											  C_accuracy,
											  (long) C_shapeType,
											  C_decay,
											  C_slope,
											  C_asymptote,
											  VolBSId,
											  C_matCurve,
											  C_volinit_vect,
											  C_coeff_vect,
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
			retCode = ARMLOCAL_PFINSTLOGDECVOLFIT(LocalGetNumObjectId(C_pfId),
												  LocalGetNumObjectId(C_secId),
												  LocalGetNumObjectId(C_curveId),
												  C_proba,
												  (long) C_UsePFResetDates,
												  C_accuracy,
												  (long) C_shapeType,
												  C_decay,
												  C_slope,
												  C_asymptote,
												  VolBSId,
												  C_matCurve,
												  C_volinit_vect,
												  C_coeff_vect,
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

			retCode = ARMLOCAL_PFINSTLOGDECVOLFIT(LocalGetNumObjectId(C_pfId),
												  LocalGetNumObjectId(C_secId),
												  LocalGetNumObjectId(C_curveId),
												  C_proba,
												  (long) C_UsePFResetDates,
												  C_accuracy,
												  (long) C_shapeType,
												  C_decay,
												  C_slope,
												  C_asymptote,
												  VolBSId,
												  C_matCurve,
												  C_volinit_vect,
												  C_coeff_vect,
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

//    ARM_END();
    }
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_PFINSTLOGDECVOLFIT" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;

}


__declspec(dllexport) LPXLOPER WINAPI Local_PXL_PFINSTLOGDECVOLFIT(LPXLOPER XL_pfId,
																   LPXLOPER XL_secId,
																   LPXLOPER XL_curveId,
																   LPXLOPER XL_proba,
																   LPXLOPER XL_UsePFResetDates,
																   LPXLOPER XL_accuracy,
																   LPXLOPER XL_shapeType,
																   LPXLOPER XL_decay,
																   LPXLOPER XL_slope,
																   LPXLOPER XL_asymptote,
																   LPXLOPER XL_VolBSId,
																   LPXLOPER XL_matCurve,
																   LPXLOPER XL_volinit_vect,
																   LPXLOPER XL_coeff_vect)
{
	ADD_LOG("Local_PXL_PFINSTLOGDECVOLFIT");
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
	CCString C_curveId;
	CCString C_pfId;
	CCString C_VolBSId;

	VECTOR<double> C_proba;
	VECTOR<double> C_matCurve;

	VECTOR<double> C_volinit_vect;
	VECTOR<double> C_coeff_vect;

	double C_accuracy;
	double C_accuracy_default  = 0.00000000001;
	double C_shapeType;
	double C_shape_default = 0;
	double C_decay, C_slope, C_asymptote;
	double C_value_default = 0.0;

	double C_UsePFResetDates;
	double C_UsePFResetDates_default = 0.0;

	// error
	static int error;
	static char* reason = "";

	XL_readStrCell(XL_pfId, C_pfId," ARM_ERR: portfolio id: object expected", C_result);
	XL_readStrCell(XL_secId, C_secId," ARM_ERR: security id: object expected",C_result);
	XL_readStrCell(XL_curveId, C_curveId," ARM_ERR: curve id: object expected",C_result);
	XL_readNumVector(XL_proba, C_proba, " ARM_ERR: Proba. Vector: array of numeric values expected",C_result);
	XL_readNumCellWD(XL_UsePFResetDates, C_UsePFResetDates, C_UsePFResetDates_default," ARM_ERR: the UsePFResetDates: numeric expected",C_result);
	XL_readNumCellWD(XL_accuracy, C_accuracy, C_accuracy_default," ARM_ERR: the accuracy: numeric expected",C_result);
	XL_readNumCellWD(XL_shapeType, C_shapeType, C_shape_default," ARM_ERR: the Shape_Type: integer expected",C_result);
	XL_readNumCellWD(XL_decay, C_decay, C_value_default," ARM_ERR: the decay: numeric expected",C_result);
	XL_readNumCellWD(XL_slope, C_slope, C_value_default," ARM_ERR: the slope: numeric expected",C_result);
	XL_readNumCellWD(XL_asymptote, C_asymptote, C_value_default," ARM_ERR: the asymptote: numeric expected",C_result);
	XL_readStrCellWD(XL_VolBSId, C_VolBSId,"DEFVOLBS"," ARM_ERR: BS Vol Curve ATM expected",C_result);

    long VolBSId = -1;
    if (!( C_VolBSId == "DEFVOLBS"))
    {
        VolBSId = LocalGetNumObjectId(C_VolBSId);
    }

    if (( XL_matCurve->xltype == xltypeMissing )
        ||
        ( XL_matCurve->xltype == xltypeNil )
       )
    {
       /*--- The size of C_matCurve must be == 0 ----*/
    }
    else
        XL_readNumVector(XL_matCurve, C_matCurve," ARM_ERR: maturity curve: array of dates expected",C_result);

    if (( XL_volinit_vect->xltype == xltypeMissing )
        ||
        ( XL_volinit_vect->xltype == xltypeNil )
       )
    {
        if (C_shapeType)
        {
            for (int k = 0; k < C_proba.size (); k++)
                C_volinit_vect.push_back(0.0);
        }
        else
        {
            for (int k = 0; k < C_matCurve.size (); k++)
                C_volinit_vect.push_back(0.0);
        }
    }
    else
        XL_readNumVector(XL_volinit_vect, C_volinit_vect," ARM_ERR: penalties VECTOR: array of numeric expected",C_result);

    if (( XL_coeff_vect->xltype == xltypeMissing )
        ||
        ( XL_coeff_vect->xltype == xltypeNil )
       )
    {
        for (int h = 0; h < C_volinit_vect.size (); h++)
            C_coeff_vect.push_back(0.0);
    }
    else
        XL_readNumVector(XL_coeff_vect, C_coeff_vect," ARM_ERR: coefficients VECTOR: array of numeric expected",C_result);

    long retCode;
    long objId;
    CCString prevClass;
    
    CCString curClass = LOCAL_LOGDECANA_CLASS;
    CCString stringId = GetLastCurCellEnvValue ();
    
    retCode = ARMLOCAL_PFINSTLOGDECVOLFIT(LocalGetNumObjectId(C_pfId),
										  LocalGetNumObjectId(C_secId),
										  LocalGetNumObjectId(C_curveId),
										  C_proba,
										  (long) C_UsePFResetDates,
										  C_accuracy,
										  (long) C_shapeType,
										  C_decay,
										  C_slope, 
										  C_asymptote,
										  VolBSId,
										  C_matCurve,
										  C_volinit_vect,
										  C_coeff_vect,
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

//    ARM_END();
    }
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_PXL_PFINSTLOGDECVOLFIT" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;

}


__declspec(dllexport) LPXLOPER WINAPI Local_PFMODFIT (LPXLOPER XL_modName,
													  LPXLOPER XL_pf,
													  LPXLOPER XL_settlement,
													  LPXLOPER XL_curve,
													  LPXLOPER XL_vList,
													  LPXLOPER XL_fList,
													  LPXLOPER XL_nagAlgo,
													  LPXLOPER XL_step,
													  LPXLOPER XL_horizon)
{
	ADD_LOG("Local_PFMODFIT ");
    //ARM_BEGIN();

    // return
    static XLOPER XL_result;
    ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();

    // C variable
    CCString C_modName;
    CCString C_pf;
    double C_settlement;

    CCString C_curve;
    long curveId;

    VECTOR<double> C_vList;
    VECTOR<double> C_fList;

    double C_nagAlgo;
    double C_nagAlgo_default = 1.0;

    double C_step;
    double C_step_default = 100;

    double C_horizon;
    double C_horizon_default = 46022; // 31/12/2025 corresponds to : 46022

    // error
    static int error;
    static char* reason = "";



    XL_readStrCell(XL_modName,C_modName," ARM_ERR: model name: string expected",C_result);
    XL_readStrCell(XL_pf,C_pf," ARM_ERR: portfolio id: object expected",C_result);
    XL_readNumCell(XL_settlement,C_settlement," ARM_ERR: as of date: date expected",C_result);
    XL_readStrCellWD(XL_curve,C_curve,"DEFAULT"," ARM_ERR: curve id: object expected",C_result);

    if((XL_vList->xltype == xltypeMissing) || (XL_vList->xltype == xltypeNil))
    {
    }
    else
    {
        XL_readNumVector(XL_vList,C_vList," ARM_ERR: values list: array of numeric expected",C_result);
    }
    if((XL_fList->xltype == xltypeMissing) || (XL_fList->xltype == xltypeNil))
    {
    }
    else
    {
        XL_readNumVector(XL_fList,C_fList," ARM_ERR: flags list: array of numeric expected",C_result);
    }

    XL_readNumCellWD(XL_nagAlgo,C_nagAlgo,C_nagAlgo_default," ARM_ERR: NAG algorithm: numeric expected",C_result);
    XL_readNumCellWD(XL_step,C_step,C_step_default," ARM_ERR: number of steps: numeric expected",C_result);
    XL_readNumCellWD(XL_horizon,C_horizon,C_horizon_default," ARM_ERR: horizon date: date expected",C_result);

    if ( C_curve == "DEFAULT" )
    {
        curveId = ARM_NULL_OBJECT;
    }
    else
    {
        curveId = LocalGetNumObjectId (C_curve);
    }

    long retCode;
    long objId;
    CCString prevClass;
    CCString curClass;
    
    if (C_modName == "GYCM")
    {
        curClass = LOCAL_GYCMODEL_CLASS;
    }
    else
    if (C_modName == "HW2F")
    {
        curClass = LOCAL_HW2FMODEL_CLASS;
    }
    else
    if (C_modName == "ZCVS")
    {
        curClass = LOCAL_YIELD_CURVE_BASIC_CLASS;
    }
    else
    if (C_modName == "BKTR")
    {
        curClass = LOCAL_BKIRTREE_CLASS;
    }
    else
    if (C_modName == "ZCSP")
    {
        curClass = LOCAL_YIELD_CURVE_BASIC_CLASS;
    }
    else
    if (C_modName == "GYLS")
    {
        curClass = LOCAL_GYCLSMODEL_CLASS;
    }
    else
    {
        C_result.setMsg ("ARM_ERR: model name invalid");
        ARM_ARG_ERR();
        return (LPXLOPER)&XL_result;
    }

    CCString stringId = GetLastCurCellEnvValue ();
    
    if(!stringId)
    {
        retCode = ARMLOCAL_PFMODFIT (C_modName,
                                LocalGetNumObjectId (C_pf),
                                C_settlement,
                                curveId,
                                C_vList,
                                C_fList,
                                (long)C_step,
                                C_horizon,
                                (long)C_nagAlgo,
                                C_result);

        if ( retCode == ARM_OK) 
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
            retCode = ARMLOCAL_PFMODFIT (C_modName,
                                    LocalGetNumObjectId (C_pf),
                                    C_settlement,
                                    curveId,
                                    C_vList,
                                    C_fList,
                                    (long)C_step,
                                    C_horizon,
                                    (long)C_nagAlgo,
                                    C_result,
                                    objId);

			if ( retCode == ARM_OK) 
			{
				LocalSetCurCellEnvValue (curClass, objId); 

				stringId = LocalMakeObjectId (objId, curClass);
			}
		}
        else
        {
            FreeCurCellContent ();
            retCode = ARMLOCAL_PFMODFIT (C_modName,
                                    LocalGetNumObjectId (C_pf),
                                    C_settlement,
                                    curveId,
                                    C_vList,
                                    C_fList,
                                    (long)C_step,
                                    C_horizon,
                                    (long)C_nagAlgo,
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

//    ARM_END();
    }
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_PFMODFIT" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}

__declspec(dllexport) LPXLOPER WINAPI Local_PXL_PFMODFIT (LPXLOPER XL_modName,
														  LPXLOPER XL_pf,
														  LPXLOPER XL_settlement,
														  LPXLOPER XL_curve,
														  LPXLOPER XL_vList,
														  LPXLOPER XL_fList,
														  LPXLOPER XL_nagAlgo,
														  LPXLOPER XL_step,
														  LPXLOPER XL_horizon)
{
	ADD_LOG("Local_PXL_PFMODFIT ");
//    ARM_BEGIN();

    // return
    static XLOPER XL_result;
    ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();

    // C variable
    CCString C_modName;
    CCString C_pf;
    double C_settlement;

    CCString C_curve;
    long curveId;

    VECTOR<double> C_vList;
    VECTOR<double> C_fList;

    double C_nagAlgo;
    double C_nagAlgo_default = 1.0;

    double C_step;
    double C_step_default = 50;

    double C_horizon;
    double C_horizon_default = 40908; // 31/12/2011 corresponds to : 40908.0

    // error
    static int error;
    static char* reason = "";

    XL_readStrCell(XL_modName,C_modName," ARM_ERR: model name: string expected",C_result);
    XL_readStrCell(XL_pf,C_pf," ARM_ERR: portfolio id: object expected",C_result);
    XL_readNumCell(XL_settlement,C_settlement," ARM_ERR: as of date: date expected",C_result);
    XL_readStrCellWD(XL_curve,C_curve,"DEFAULT"," ARM_ERR: curve id: object expected",C_result);

    if((XL_vList->xltype == xltypeMissing) || (XL_vList->xltype == xltypeNil))
    {
    }
    else
    {
        XL_readNumVector(XL_vList,C_vList," ARM_ERR: values list: array of numeric expected",C_result);
    }
    if((XL_fList->xltype == xltypeMissing) || (XL_fList->xltype == xltypeNil))
    {
    }
    else
    {
        XL_readNumVector(XL_fList,C_fList," ARM_ERR: flags list: array of numeric expected",C_result);
    }

    XL_readNumCellWD(XL_nagAlgo,C_nagAlgo,C_nagAlgo_default," ARM_ERR: NAG algorithm: numeric expected",C_result);
    XL_readNumCellWD(XL_step,C_step,C_step_default," ARM_ERR: number of steps: numeric expected",C_result);
    XL_readNumCellWD(XL_horizon,C_horizon,C_horizon_default," ARM_ERR: horizon date: date expected",C_result);

    if ( C_curve == "DEFAULT" )
    {
        curveId = ARM_NULL_OBJECT;
    }
    else
    {
        curveId = LocalGetNumObjectId (C_curve);
    }

    long retCode;
    long objId;
    CCString curClass;
    
    if (C_modName == "GYCM")
    {
        curClass = LOCAL_GYCMODEL_CLASS;
    }
    else
    if (C_modName == "HW2F")
    {
        curClass = LOCAL_HW2FMODEL_CLASS;
    }
    else
    if (C_modName == "ZCVS")
    {
        curClass = LOCAL_YIELD_CURVE_BASIC_CLASS;
    }
    else
    if (C_modName == "BKTR")
    {
        curClass = LOCAL_BKIRTREE_CLASS;
    }
    else
    if (C_modName == "ZCSP")
    {
        curClass = LOCAL_YIELD_CURVE_BASIC_CLASS;
    }
    else
    if (C_modName == "GYLS")
    {
        curClass = LOCAL_GYCLSMODEL_CLASS;
    }
    else
    {
        C_result.setMsg ("ARM_ERR: model name invalid");
        ARM_ARG_ERR();
        return (LPXLOPER)&XL_result;
    }

    CCString stringId;
    
    retCode = ARMLOCAL_PFMODFIT (C_modName,
                            LocalGetNumObjectId (C_pf),
                            C_settlement,
                            curveId,
                            C_vList,
                            C_fList,
                            (long)C_step,
                            C_horizon,
                            (long)C_nagAlgo,
                            C_result);

    if ( retCode == ARM_OK) 
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

//    ARM_END();
   }
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_PXL_PFMODFIT" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


//YK CalibratePF////////////////////////////////////////////////////////////////////////////////////

__declspec(dllexport) LPXLOPER WINAPI Local_ARM_CalibrationPF(LPXLOPER XL_curve,
															  LPXLOPER XL_vol,
															  LPXLOPER XL_sec,
															  LPXLOPER XL_modeType,
															  LPXLOPER XL_portfolioType,
															  LPXLOPER XL_shift)
{
	ADD_LOG("Local_ARM_CalibrationPF");
	// return
	static XLOPER XL_result;
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();

	// C variable
	
	CCString C_Zc;
	CCString C_Vol;
	CCString C_sec;

	CCString C_modeType;
	CCString C_modeType_default = "SUMMIT";
	long modeTypeId;


	double C_portfolioType;
	double C_Shift;

	double C_portfolioTypeDefault = 0;
	double C_ShiftDefault = -0.5;
		
	// error
	static int error;
	static char* reason = "";

	XL_readStrCell(XL_curve,C_Zc," ARM_ERR: curve id: object expected",C_result);
	XL_readStrCell(XL_vol,C_Vol," ARM_ERR: volcurve id: object expected",C_result);
	XL_readStrCell(XL_sec,C_sec," ARM_ERR: security id: object expected",C_result);
	XL_readStrCellWD(XL_modeType,C_modeType,C_modeType_default," ARM_ERR: mode Type: numeric expected",C_result);

	XL_readNumCellWD(XL_portfolioType,C_portfolioType,C_portfolioTypeDefault," ARM_ERR: portfolio Type: numeric expected",C_result);
	XL_readNumCellWD(XL_shift,C_Shift,C_ShiftDefault," ARM_ERR: shift: numeric expected",C_result);

	
	
	if((modeTypeId = ARM_ConvSummitManual (C_modeType, C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	long retCode;
	long objId;
	CCString prevClass;
	
	CCString curClass = LOCAL_PF_CLASS;
	CCString stringId = GetLastCurCellEnvValue ();

	if (!stringId)
	{
		retCode = ARMLOCAL_CalibrationPF (LocalGetNumObjectId (C_Zc),
										  LocalGetNumObjectId (C_Vol),
										  LocalGetNumObjectId (C_sec),
										  (long) modeTypeId,
										  (long) C_portfolioType,
										  C_Shift,
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
			retCode = ARMLOCAL_CalibrationPF (LocalGetNumObjectId (C_Zc),
											  LocalGetNumObjectId (C_Vol),
											  LocalGetNumObjectId (C_sec),
											  (long) modeTypeId,
											  (long) C_portfolioType,
											  C_Shift,
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

			retCode = ARMLOCAL_CalibrationPF (LocalGetNumObjectId (C_Zc),
											  LocalGetNumObjectId (C_Vol),
											  LocalGetNumObjectId (C_sec),
											  (long) modeTypeId,
											  (long) C_portfolioType,
											  C_Shift,
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_ARM_CalibrationPF" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}



__declspec(dllexport) LPXLOPER WINAPI Local_PXL_ARM_CalibrationPF (LPXLOPER XL_curve,
																   LPXLOPER XL_vol,
																   LPXLOPER XL_sec,
																   LPXLOPER XL_modeType,
																   LPXLOPER XL_portfolioType,
																   LPXLOPER XL_shift)
{
	ADD_LOG("Local_PXL_ARM_CalibrationPF ");
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();

	// C variable
	CCString C_Zc;
	CCString C_Vol;
	CCString C_sec;

	CCString C_modeType;
	CCString C_modeType_default = "SUMMIT";
	long modeTypeId;

	double C_portfolioType;
	double C_Shift;

	double C_portfolioTypeDefault = 0;
	double C_ShiftDefault = -0.5;


	// error
	static int error;
	static char* reason = "";

	XL_readStrCell(XL_curve,C_Zc," ARM_ERR: curve id: object expected",C_result);
	XL_readStrCell(XL_vol,C_Vol," ARM_ERR: volcurve id: object expected",C_result);
	XL_readStrCell(XL_sec,C_sec," ARM_ERR: security id: object expected",C_result);
	XL_readStrCellWD(XL_modeType,C_modeType,C_modeType_default," ARM_ERR: mode Type: numeric expected",C_result);
	XL_readNumCellWD(XL_portfolioType,C_portfolioType,C_portfolioTypeDefault," ARM_ERR: portfolio Type: numeric expected",C_result);
	XL_readNumCellWD(XL_shift,C_Shift,C_ShiftDefault," ARM_ERR: shift: numeric expected",C_result);

	
	if((modeTypeId = ARM_ConvSummitManual (C_modeType, C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	long retCode;
	long objId;
	
	CCString curClass = LOCAL_CALIB_HWSV_CLASS;
	CCString stringId;

	retCode = ARMLOCAL_CalibrationPF (LocalGetNumObjectId (C_Zc),
									  LocalGetNumObjectId (C_Vol),
									  LocalGetNumObjectId (C_sec),
									  (long) modeTypeId,
									  (long) C_portfolioType,
									  C_Shift,
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_PXL_ARM_CalibrationPF" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}



//YK CalibratePF////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////Portfolio to calibrate Mean Reversion/////////////////////////////


__declspec(dllexport) LPXLOPER WINAPI Local_ARM_MRSCalibrationPF(LPXLOPER XL_curve,
															  LPXLOPER XL_vol,
															  LPXLOPER XL_sec,
															  LPXLOPER XL_firstPortfolio,
															  LPXLOPER XL_freq,
															  LPXLOPER XL_atmflag)
{
	ADD_LOG("Local_ARM_MRSCalibrationPF");
	// return
	static XLOPER XL_result;
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();

	// C variable
	
	CCString C_Zc;
	CCString C_Vol;
	CCString C_sec;
	CCString C_Portfolio;

	double C_Freq;
	double c_ATMFlag;

	// error
	static int error;
	static char* reason = "";

	XL_readStrCell(XL_curve,C_Zc," ARM_ERR: curve id: object expected",C_result);
	XL_readStrCell(XL_vol,C_Vol," ARM_ERR: volcurve id: object expected",C_result);
	XL_readStrCell(XL_sec,C_sec," ARM_ERR: security id: object expected",C_result);
	XL_readStrCell(XL_firstPortfolio,C_Portfolio," ARM_ERR: portfolio id: object expected",C_result);

	XL_readNumCell(XL_freq,C_Freq," ARM_ERR: portfolio Type: numeric expected",C_result);
	XL_readNumCell(XL_atmflag,c_ATMFlag," ARM_ERR: shift: numeric expected",C_result);

	bool atmFlag = false;
	if (c_ATMFlag!=0)
		atmFlag = true;


	long retCode;
	long objId;
	CCString prevClass;
	
	CCString curClass = LOCAL_PF_CLASS;
	CCString stringId = GetLastCurCellEnvValue ();

	if (!stringId)
	{
		retCode = ARMLOCAL_MRSCalibrationPF (LocalGetNumObjectId (C_Zc),
										  LocalGetNumObjectId (C_Vol),
										  LocalGetNumObjectId (C_sec),
										  LocalGetNumObjectId (C_Portfolio),
										  (long) C_Freq,
										  atmFlag,
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
			retCode = ARMLOCAL_MRSCalibrationPF (LocalGetNumObjectId (C_Zc),
											  LocalGetNumObjectId (C_Vol),
											  LocalGetNumObjectId (C_sec),
											  LocalGetNumObjectId (C_Portfolio),
											  (long) C_Freq,
											  atmFlag,
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

			retCode = ARMLOCAL_MRSCalibrationPF (LocalGetNumObjectId (C_Zc),
											  LocalGetNumObjectId (C_Vol),
											  LocalGetNumObjectId (C_sec),
											  LocalGetNumObjectId (C_Portfolio),
											  (long) C_Freq,
											  atmFlag,
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_ARM_MRSCalibrationPF" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}

__declspec(dllexport) LPXLOPER WINAPI Local_PXL_ARM_MRSCalibrationPF(LPXLOPER XL_curve,
															  LPXLOPER XL_vol,
															  LPXLOPER XL_sec,
															  LPXLOPER XL_firstPortfolio,
															  LPXLOPER XL_freq,
															  LPXLOPER XL_atmflag)
{
	ADD_LOG("Local_PXL_ARM_MRSCalibrationPF");
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();

	CCString C_Zc;
	CCString C_Vol;
	CCString C_sec;
	CCString C_Portfolio;

	double C_Freq;
	double c_ATMFlag;

	// error
	static int error;
	static char* reason = "";

	XL_readStrCell(XL_curve,C_Zc," ARM_ERR: curve id: object expected",C_result);
	XL_readStrCell(XL_vol,C_Vol," ARM_ERR: volcurve id: object expected",C_result);
	XL_readStrCell(XL_sec,C_sec," ARM_ERR: security id: object expected",C_result);
	XL_readStrCell(XL_firstPortfolio,C_Portfolio," ARM_ERR: portfolio id: object expected",C_result);

	XL_readNumCell(XL_freq,C_Freq," ARM_ERR: portfolio Type: numeric expected",C_result);
	XL_readNumCell(XL_atmflag,c_ATMFlag," ARM_ERR: shift: numeric expected",C_result);

	bool atmFlag = false;
	if (c_ATMFlag!=0)
		atmFlag = true;

	long retCode;
	long objId;
	
	CCString curClass = LOCAL_CALIB_HWSV_CLASS;
	CCString stringId;

	retCode = ARMLOCAL_MRSCalibrationPF (LocalGetNumObjectId (C_Zc),
											  LocalGetNumObjectId (C_Vol),
								 			  LocalGetNumObjectId (C_sec),
											  LocalGetNumObjectId (C_Portfolio),
											  (long) C_Freq,
											  atmFlag,
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_PXL_ARM_CalibrationPF" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}

///----------------------------------------------
///----------------------------------------------
///             Get from Id
/// Inputs :
///     Portfolio Id
///     Index of Asset Id
///----------------------------------------------
///----------------------------------------------
////////////////////////////////////////////////
/// very rapid definition of ourlonglongFunctor
////////////////////////////////////////////////
class GetAssetFromPFFunc : public ARMResultLong2LongFunc
{
public:
	GetAssetFromPFFunc(long portfolioId,long assetIdIndex)
                    :C_portfolioId(portfolioId), 
                     C_assetIdIndex(assetIdIndex)
    {};    
	
    long operator()( ARM_result& result, long objId ){
        return ARMLOCAL_GetAssetFromPF(C_portfolioId, C_assetIdIndex,result, objId);
    }
			

private:
	long C_portfolioId;
    long C_assetIdIndex;
};



/////////////////////////////////////////////////////////////
/// central function that does the creation of the XL function
/////////////////////////////////////////////////////////////
LPXLOPER Local_GetAssetFromPFCommon(
	LPXLOPER XL_PortfolioId,
    LPXLOPER XL_AssetIdIndex,
	bool PersistentInXL )
{
	/// this is defined first because it is used in XL macros
	static XLOPER XL_result;
	/// Get the variables from the XLOper variables
	ARM_result C_result;	

	/// to avoid computation if called by the wizard
	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();
		
	/// this is used by macros 
	/// and therefore this has to be defined
	static int error;
	static char* reason = "";	

    CCString C_portfolioStrId;
	XL_readStrCell( XL_PortfolioId, C_portfolioStrId,	" ARM_ERR: Portfolio Id: Object expected",C_result);
	long C_portfolioId = LocalGetNumObjectId(C_portfolioStrId);


    double assetid_idex = 0;
    double C_AssetIndex;
	XL_readNumCellWD(XL_AssetIdIndex, C_AssetIndex, assetid_idex, " ARM_ERR: Asset id Index: numeric expected",C_result);	  

	/// use the concept of Functor to transfer the knowledge of
	/// a function with a context
	GetAssetFromPFFunc ourFunc(C_portfolioId, (long)C_AssetIndex);

    /// given the name class appeared in excel
	string	vSecurityTypeToClass = SecurityTypeToClass(C_portfolioId, (long)C_AssetIndex);
    CCString AssetGetClass(vSecurityTypeToClass.c_str());


	/// call the general function
    fillXL_Result(AssetGetClass, ourFunc, C_result, XL_result, PersistentInXL );

	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_GetAssetFromPFCommon" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}	

///////////////////////////////////
/// version that takes into account 
/// previous creation of object
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_GetAssetFromPF(
	LPXLOPER XL_PortfolioId,
    LPXLOPER XL_AssetIdIndex)
{
	ADD_LOG("Local_GetAssetFromPF");
	bool PersistentInXL = true;
	return Local_GetAssetFromPFCommon( XL_PortfolioId,
                                    XL_AssetIdIndex,
                                    PersistentInXL );
}

///////////////////////////////////
/// version for non persistentInXl
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_PXL_GetAssetFromPF(
	LPXLOPER XL_PortfolioId,
    LPXLOPER XL_AssetIdIndex)
{
	ADD_LOG("Local_PXL_GetAssetFromPF");
	bool PersistentInXL = false;
	return Local_GetAssetFromPFCommon(XL_PortfolioId,
                                    XL_AssetIdIndex,
                                    PersistentInXL );
}

////////////////////////////////////////////////
/// very rapid definition of ourlonglongFunctor
////////////////////////////////////////////////
class GetSizeOfPFFunc : public ARMResultLong2LongFunc
{
public:
	GetSizeOfPFFunc(long portfolioId)
                    :C_portfolioId(portfolioId)
    {};             
	
    long operator()( ARM_result& result, long objId ){
        return ARMLOCAL_GetSizeOfPF(C_portfolioId,C_SizeofPf,result);
    }

    inline long  const GetSize() const { return C_SizeofPf; };			

private:
	long C_portfolioId;
    long C_SizeofPf;
};

////////////////////////////////////////////////////////////
/// central function that does the creation of the XL function
/////////////////////////////////////////////////////////////
LPXLOPER Local_GetSizeOfPFCommon(
	LPXLOPER XL_PortfolioId,
	bool PersistentInXL )
{
	/// this is defined first because it is used in XL macros
	static XLOPER XL_result;
	/// Get the variables from the XLOper variables
	ARM_result C_result;	

	/// to avoid computation if called by the wizard
	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();
		
	/// this is used by macros 
	/// and therefore this has to be defined
	static int error;
	static char* reason = "";	

    CCString C_portfolioStrId;
	XL_readStrCell( XL_PortfolioId, C_portfolioStrId,	" ARM_ERR: Portfolio Id: Object expected",C_result);
	long C_portfolioId = LocalGetNumObjectId(C_portfolioStrId);
    
	/// use the concept of Functor to transfer the knowledge of
	/// a function with a context
	GetSizeOfPFFunc ourFunc(C_portfolioId);  

    if(ourFunc(C_result,C_portfolioId) == ARM_OK)
    {
        XL_result.xltype  = xltypeNum;
        XL_result.val.num = ourFunc.GetSize();
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_GetSizeOfPFCommon" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
    


}						 					 


///////////////////////////////////
/// version that takes into account 
/// previous creation of object
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_GetSizeOfPF(
	LPXLOPER XL_PortfolioId)
{
	ADD_LOG("Local_GetSizeOfPF");
	bool PersistentInXL = true;
	return Local_GetSizeOfPFCommon( XL_PortfolioId,
                                    PersistentInXL );
}


__declspec(dllexport) LPXLOPER WINAPI Local_ComputeAll(	LPXLOPER XL_pfId,
													    LPXLOPER XL_modelId)
{
	ADD_LOG("Local_ComputeAll");
//    ARM_BEGIN();

    // return
    static XLOPER XL_result;
    ARM_result C_result;

	LPXLOPER pxArray;

		/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();

    // C variable
	CCString C_pfId;
	CCString C_modelId;
    
    // error
    static int error;
    static char* reason = "";

	XL_readStrCell(XL_pfId,C_pfId," ARM_ERR: pf id: object expected",C_result);
	XL_readStrCell(XL_modelId,C_modelId," ARM_ERR: model id: object expected",C_result);

    long retCode;

    retCode = ARMLOCAL_ComputeAll(LocalGetNumObjectId(C_pfId),LocalGetNumObjectId(C_modelId),C_result);

    if ( retCode == ARM_OK )
    {
        int nbrows = int (C_result.getDouble());
		int nbcolumns = 1;

		FreeCurCellErr ();
        
		XL_result.xltype = xltypeMulti;
        XL_result.val.array.columns = nbcolumns;
		XL_result.val.array.rows = nbrows;
		XL_result.val.array.lparray = pxArray = (LPXLOPER)GlobalAlloc (GMEM_ZEROINIT, nbrows * nbcolumns * sizeof (XLOPER));

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

//    ARM_END();
   	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_ComputeAll" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
    

}
