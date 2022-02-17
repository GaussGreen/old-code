/*
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file ARM_xl_gp_nummethod_local.cpp
 *
 *  \brief file for the model part in the generic pricer
 *
 *	\author  E Benhamou
 *	\version 1.0
 *	\date September 2003
 */

#include <ARM\libarm_local\firstToBeIncluded.h>
#include <functional>
#include <libCCxll\CCxll.h>

#include <ARM\libarm_local\ARM_local_gp_model.h>
#include "ARM_xl_gp_model_local.h"
#include "ARM_xl_gp_nummethod_local.h"
#include "ARM_xl_wrapper_local.h"
#include "ARM_gp_local_interglob.h"
#include "ARM_xl_trycatch_local.h"
#include "ARM_xl_gp_fctorhelper.h"

#include <GP_Base\gpbase\numericconstant.h>
#include <GP_Base\gpbase\stringmanip.h>
#include <GP_Models\gpmodels\argconvdefault.h>

#include <GP_Base\gpbase\gpvector.h>
using ARM::ARM_GP_Vector;
using ARM::ARM_NumericConstants;
using ARM::ARM_ArgConv_VnsPricingMethod;
#include <util\fromto.h>

#include "util\tech_macro.h"

CC_USING_NS(std,string)
using ARM::stringGetUpper;


#ifdef _PURIFY
	const  struct fpos_t std::_Fpz = {0, 0};
#endif

////////////////////////////////////////////////
/// 
///  THIS INTERFACE IS MORE FACTORISED THAN THE
///  REST OF ARM .... IN ORDER TO RECODE
///  ONLY THE PART TO DEFINE A FUNCTION CALL
///  WE PROVIDE FOR EACH A FUNCTOR
///  WE ALSO USE THE SAME API FOR PERSISTENT 
///  AND NON PERSISTENT ADDIN
///  BY COMMON DECISION.... WE KEEP THIS METHOD
///  TO THIS FILES AND DO NOT SPREAD ELSEWHERE
/// 
////////////////////////////////////////////////


////////////////////////////////////////////////
/// very rapid definition of ourlonglongFunctor
////////////////////////////////////////////////
class IRFwdModFunc : public ARMResultLong2LongFunc
{
public:
	IRFwdModFunc(long zeroCurveId)
	:	C_zeroCurveId( zeroCurveId){};
	
	long operator()( ARM_result& result, long objId ){
		return ARMLOCAL_IRFwd_Create(C_zeroCurveId, result, objId );			
	}

private:
	long C_zeroCurveId;
};


/////////////////////////////////////////////////////////////
/// central function that does the creation of the XL function
/////////////////////////////////////////////////////////////
LPXLOPER Local_IRFwdModCommon(
	LPXLOPER XL_ZeroCurveId,
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
		long C_ZeroCurveId;
		
		XL_GETOBJID( XL_ZeroCurveId, C_ZeroCurveId,	" ARM_ERR: Zero Curve: Object expected",		C_result);
		
		/// use the concept of Functor to transfer the knowledge of
		/// a function with a context
		IRFwdModFunc ourFunc(C_ZeroCurveId);

		/// call the general function
		fillXL_Result( LOCAL_IRFWDMOD_CLASS, ourFunc, C_result, XL_result, PersistentInXL );
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_IRFwdModCommon" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}
							 

///////////////////////////////////
/// version that takes into account 
/// previous creation of object
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_IRFwdMod_Create(
	LPXLOPER XL_ZeroCurveId )
{
	ADD_LOG("Local_IRFwdMod_Create");
	bool PersistentInXL = true;
	return Local_IRFwdModCommon( XL_ZeroCurveId, PersistentInXL );
}


///////////////////////////////////
/// version for non persistentInXl
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_PXL_IRFwdMod_Create(
	LPXLOPER XL_ZeroCurveId )
{
	ADD_LOG("Local_PXL_IRFwdMod_Create");
	bool PersistentInXL = false;
	return Local_IRFwdModCommon( XL_ZeroCurveId, PersistentInXL );
}


///----------------------------------------------
///----------------------------------------------
///             InfFwdModel
/// Inputs :
///     ZcCurve
///     InfCurve
///		IRModel
///----------------------------------------------
///----------------------------------------------

////////////////////////////////////////////////
/// very rapid definition of ourlonglongFunctor
////////////////////////////////////////////////
class InfFwdModFunc : public ARMResultLong2LongFunc
{
public:
	InfFwdModFunc( long infCurveId )
	:	C_infCurveId(infCurveId) {};
	
	long operator()( ARM_result& result, long objId ){
		return ARMLOCAL_InfFwd_Create(C_infCurveId,  result, objId );			
	}

private:
	long C_infCurveId;
};


/////////////////////////////////////////////////////////////
/// central function that does the creation of the XL function
/////////////////////////////////////////////////////////////
LPXLOPER Local_InfFwdModCommon(
	LPXLOPER XL_InfCurveId,
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
		long C_InfCurveId;
		
		XL_GETOBJID( XL_InfCurveId, C_InfCurveId,	" ARM_ERR: Inflation Curve: Object expected",		C_result);

		/// use the concept of Functor to transfer the knowledge of
		/// a function with a context
		InfFwdModFunc ourFunc(C_InfCurveId);

		/// call the general function
		fillXL_Result_withName( ourFunc, C_result, XL_result, PersistentInXL );
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_InfFwdModCommon" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}
							 

///////////////////////////////////
/// version that takes into account 
/// previous creation of object
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_InfFwdMod_Create(
	LPXLOPER XL_InfCurveId)
{
	ADD_LOG("Local_InfFwdMod_Create");
	bool PersistentInXL = true;
	return Local_InfFwdModCommon( XL_InfCurveId, PersistentInXL );
}


///////////////////////////////////
/// version for non persistentInXl
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_PXL_InfFwdMod_Create(
	LPXLOPER XL_InfCurveId )
{
	ADD_LOG("Local_PXL_InfFwdMod_Create");
	bool PersistentInXL = false;
	return Local_InfFwdModCommon( XL_InfCurveId, PersistentInXL );
}

///----------------------------------------------
///----------------------------------------------
///             Inflation Equity Model
/// Inputs :
///     ZcCurve
///     InfCurve
///		IRModel
///----------------------------------------------
///----------------------------------------------

////////////////////////////////////////////////
/// very rapid definition of ourlonglongFunctor
////////////////////////////////////////////////
class InflationEquityModelFunc : public ARMResultLong2LongFunc
{
public:
	InflationEquityModelFunc(long infCurveId, double PublicationLag, long Param1Id, long Param2Id )
	:	C_infCurveId(infCurveId),C_PublicationLag(PublicationLag),C_Param1Id(Param1Id),C_Param2Id(Param2Id) {};
	
	long operator()( ARM_result& result, long objId ){
		return ARMLOCAL_InflationEquityModel_Create(C_infCurveId, C_PublicationLag,C_Param1Id,C_Param2Id,result, objId );			
	}

private:
	long C_infCurveId,C_Param1Id,C_Param2Id;
	double C_PublicationLag;
};






/////////////////////////////////////////////////////////////
/// central function that does the creation of the XL function
/////////////////////////////////////////////////////////////
LPXLOPER Local_InflationEquityModelCommon(
	LPXLOPER XL_InfCurveId,
	LPXLOPER XL_PublicationLag,
	LPXLOPER XL_Param1Id,
	LPXLOPER XL_Param2Id,
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
		long C_InfCurveId, C_Param1Id, C_Param2Id;
		double C_PublicationLag;
		
		XL_GETOBJID( XL_InfCurveId, C_InfCurveId,	" ARM_ERR: Inflation Curve: Object expected",		C_result);
		XL_readNumCell( XL_PublicationLag, C_PublicationLag, " ARM_ERR: publication lag: double expected",C_result);		
		XL_GETOBJID( XL_Param1Id, C_Param1Id,	" ARM_ERR: ModelParam1 : Object expected",		C_result);
		XL_GETOBJID( XL_Param2Id, C_Param2Id,	" ARM_ERR: ModelParam2 : Object expected",		C_result);
		
		/// use the concept of Functor to transfer the knowledge of
		/// a function with a context
		InflationEquityModelFunc ourFunc(C_InfCurveId,C_PublicationLag, C_Param1Id, C_Param2Id);

		/// call the general function
		fillXL_Result_withName( ourFunc, C_result, XL_result, PersistentInXL );
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_InfFwdModCommon" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}
							 

///////////////////////////////////
/// version that takes into account 
/// previous creation of object
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_InflationEquityModel_Create(
	LPXLOPER XL_InfCurveId,LPXLOPER XL_PublicationLag,LPXLOPER XL_Param1Id,LPXLOPER XL_Param2Id)
{
	ADD_LOG("Local_InflationEquityModel_Create");
	bool PersistentInXL = true;
	return Local_InflationEquityModelCommon( XL_InfCurveId, XL_PublicationLag, XL_Param1Id,XL_Param2Id, PersistentInXL );
}


///////////////////////////////////
/// version for non persistentInXl
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_PXL_InflationEquityModel_Create(
	LPXLOPER XL_InfCurveId,LPXLOPER XL_PublicationLag,LPXLOPER XL_Param1Id,LPXLOPER XL_Param2Id)
{
	ADD_LOG("Local_PXL_InflationEquityModel_Create");
	bool PersistentInXL = false;
;	return Local_InflationEquityModelCommon( XL_InfCurveId, XL_PublicationLag, XL_Param1Id,XL_Param2Id, PersistentInXL );
}


///----------------------------------------------
///----------------------------------------------
///             SABR  Equity Model
/// Inputs :
///     ZcCurve
///		DvdCurve
///		Spot
///     beta (double)
///		alpha (double)
///		rho (double)
///		nu (double)
///----------------------------------------------
///----------------------------------------------
///----------------------------------------------



/////////////////////////////////////////////////////////////
/// central function that does the creation of the XL function
/////////////////////////////////////////////////////////////
LPXLOPER Local_SABREquityModelCommon(
	LPXLOPER XL_ZcCurveId,
	LPXLOPER XL_Spot,
	LPXLOPER XL_modelParamsIds,
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
		long C_ZcCurveId;
		XL_GETOBJID( XL_ZcCurveId, C_ZcCurveId,	" ARM_ERR: Zc Curve: Object expected",		C_result);

		double C_Spot;
	    XL_readNumCell(XL_Spot,C_Spot,		" ARM_ERR: Spot numeric expected",C_result);	

		vector<CCString> C_paramsIds;
		XL_readStrVector (XL_modelParamsIds,C_paramsIds," ARM_ERR: Model Parameters: array of object expected",DOUBLE_TYPE,C_result);
		size_t i, size = C_paramsIds.size();
		vector<long> C_paramsIdVec(size);
		for(i = 0; i < size; ++i )    
			C_paramsIdVec[i] = LocalGetNumObjectId(C_paramsIds[i]); 

		/// a function with a context
		exportFunc3Args< long, double, vector< long > >  ourFunc(C_ZcCurveId, C_Spot, C_paramsIdVec, ARMLOCAL_SABREquityModel_Create );
		
		/// call the general function
		fillXL_Result_withName( ourFunc, C_result, XL_result, PersistentInXL );
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_InfFwdModCommon" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}
		

///////////////////////////////////
/// version that takes into account 
/// previous creation of object
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_SABREquityModel_Create(
	LPXLOPER XL_ZcCurveId,
	LPXLOPER XL_Spot,
	LPXLOPER XL_modelParamsIds )
{
	ADD_LOG("Local_SABREquityModel_Create");
	bool PersistentInXL = true;
	return Local_SABREquityModelCommon( 
		XL_ZcCurveId,
		XL_Spot,		
		XL_modelParamsIds,
		PersistentInXL );
}


///////////////////////////////////
/// version for non persistentInXl
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_PXL_SABREquityModel_Create(
	LPXLOPER XL_ZcCurveId,
	LPXLOPER XL_Spot,
	LPXLOPER XL_modelParamsIds )
{
	ADD_LOG("Local_PXL_SABREquityModel_Create");
	bool PersistentInXL = false;
	return Local_SABREquityModelCommon( 
		XL_ZcCurveId,
		XL_Spot,		
		XL_modelParamsIds,
		PersistentInXL );
}


///----------------------------------------------
///----------------------------------------------
///             Model Parameter
/// Inputs :
///     Type
///     Time lags vector
///     Value vector
///     Name
///----------------------------------------------
///----------------------------------------------
////////////////////////////////////////////////
/// very rapid definition of ourlonglongFunctor
////////////////////////////////////////////////
class ModelParamFunc : public ARMResultLong2LongFunc
{
public:
	ModelParamFunc(
        long modelParamType,
        const vector<double>&  paramTimes,
        const vector<double>&  paramValues,
        const CCString&	modelParamName,
        const vector<double>& LowerBoundary,
        const vector<double>& UpperBoundary,
        const CCString& InterpolMethodName,
		bool adviseBreakPointTimes,
		const CCString& currency)
    :
    C_ModelParamType(modelParamType),
    C_ParamTimes(paramTimes),
    C_ParamValues(paramValues),
    C_ModelParamName(modelParamName),
    C_LowerBoundary(LowerBoundary),
    C_UpperBoundary(UpperBoundary),
    C_InterpolMethodName(InterpolMethodName),
	C_AdviseBreakPointTimes(adviseBreakPointTimes),
	C_Currency(currency)
    {};
	
	long operator()( ARM_result& result, long objId ){
		return ARMLOCAL_ModelParam_Create(
            C_ModelParamType,
            C_ParamTimes,
            C_ParamValues,
            C_ModelParamName,
            C_LowerBoundary,
            C_UpperBoundary,
            C_InterpolMethodName,
			C_AdviseBreakPointTimes,
			C_Currency,
            result, objId);			
	}

private:
    long            C_ModelParamType;
    vector<double>  C_ParamTimes;
    vector<double>  C_ParamValues;
    CCString        C_ModelParamName;
    vector<double>  C_LowerBoundary;
    vector<double>  C_UpperBoundary;
    const CCString& C_InterpolMethodName;
	bool			C_AdviseBreakPointTimes;
	CCString		C_Currency;
};


/////////////////////////////////////////////////////////////
/// central function that does the creation of the XL function
/////////////////////////////////////////////////////////////
LPXLOPER Local_ModelParamCommon(
	LPXLOPER XL_ModelParamType,
	LPXLOPER XL_ParamTimes,
	LPXLOPER XL_ParamValues,
	LPXLOPER XL_ModelParamName,
    LPXLOPER XL_LowerBoundary,
    LPXLOPER XL_UpperBoundary,
    LPXLOPER XL_InterpolMethod,
	LPXLOPER XL_AdviseBreakPointTimes,
	LPXLOPER XL_Currency,
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

		CCString C_ModelParamTypeStr;
		XL_readStrCell(XL_ModelParamType,C_ModelParamTypeStr," ARM_ERR: Param Type: String expected",C_result);
		long C_ModelParamType;
		if( (C_ModelParamType = ARM_ConvGPModelParam( C_ModelParamTypeStr, C_result)) == ARM_DEFAULT_ERR )
		{
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}

		vector<double> defaultVector(0);
		defaultVector.clear();

		vector<double> C_ParamTimes;
		XL_readNumVectorWD(XL_ParamTimes,C_ParamTimes,defaultVector," ARM_ERR: Param Times: array of numeric expected",C_result);

		vector<double> C_ParamValues;
		if(C_ParamTimes.size() > 0)
        {
		    XL_readNumVector(XL_ParamValues,C_ParamValues," ARM_ERR: Param Values: array of numeric expected",C_result);
        }

		CCString C_ModelParamName;
		CCString defaultModelParamName="";
		XL_readStrCellWD(XL_ModelParamName,C_ModelParamName,defaultModelParamName," ARM_ERR: Param Name: String expected",C_result);


		vector<double> C_LowerBoundary;
		vector<double> C_UpperBoundary;

		XL_readNumVectorWD(XL_LowerBoundary,C_LowerBoundary,defaultVector," ARM_ERR: lower bound: array of numeric expected",C_result);
		if(C_LowerBoundary.size() > 0)
        {
			XL_readNumVector(XL_UpperBoundary,C_UpperBoundary," ARM_ERR: Upper bound: array of numeric expected",C_result);
        }
		else 
			C_UpperBoundary.clear(); 
        
        CCString C_InterpolMethod;
        CCString default_InterpolMethod = "STEPUPRIGHT";
		XL_readStrCellWD(XL_InterpolMethod, C_InterpolMethod, default_InterpolMethod, " ARM_ERR: InterpolMethod : numeric expected",C_result);	  

		
		double C_AdviseBreakPointTimes;
		double C_AdviseBreakPointTimesDef=0;
		XL_readNumCellWD( XL_AdviseBreakPointTimes, C_AdviseBreakPointTimes, C_AdviseBreakPointTimesDef, " ARM_ERR: keep times: boolean expected",C_result);
		bool C_AdviseBreakPointTimesBool = C_AdviseBreakPointTimes != 0;



		CCString C_Currency;
		CCString defaultCurrency="";
		XL_readStrCellWD(XL_Currency,C_Currency,defaultCurrency," ARM_ERR: Currency Name: String expected",C_result);


		/// use the concept of Functor to transfer the knowledge of
		/// a function with a context
		ModelParamFunc ourFunc(
			C_ModelParamType,
			C_ParamTimes,
			C_ParamValues,
			C_ModelParamName,
			C_LowerBoundary,
			C_UpperBoundary,
            C_InterpolMethod,
			C_AdviseBreakPointTimesBool,
			C_Currency)	;

		/// call the general function
		fillXL_Result( LOCAL_MODELPARAM_CLASS, ourFunc, C_result, XL_result, PersistentInXL );
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_ModelParamCommon" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}
							 

///////////////////////////////////
/// version that takes into account 
/// previous creation of object
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_ModelParam_Create(
	LPXLOPER XL_ModelParamType,
	LPXLOPER XL_ParamTimes,
	LPXLOPER XL_ParamValues,
	LPXLOPER XL_ModelParamName,
    LPXLOPER XL_LowerBoundary,
    LPXLOPER XL_UpperBoundary,
    LPXLOPER XL_InterpolMethod,
	LPXLOPER XL_AdviseBreakPointTimes,
	LPXLOPER XL_Currency)
{
	ADD_LOG("Local_ModelParam_Create");
	bool PersistentInXL = true;
	return Local_ModelParamCommon(
	    XL_ModelParamType,
	    XL_ParamTimes,
	    XL_ParamValues,
	    XL_ModelParamName,
        XL_LowerBoundary,
        XL_UpperBoundary,
        XL_InterpolMethod,
		XL_AdviseBreakPointTimes,
		XL_Currency,
        PersistentInXL );
}


///////////////////////////////////
/// version for non persistentInXl
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_PXL_ModelParam_Create(
	LPXLOPER XL_ModelParamType,
	LPXLOPER XL_ParamTimes,
	LPXLOPER XL_ParamValues,
	LPXLOPER XL_ModelParamName,
    LPXLOPER XL_LowerBoundary,
    LPXLOPER XL_UpperBoundary,
    LPXLOPER XL_InterpolMethod,
	LPXLOPER XL_AdviseBreakPointTimes,
	LPXLOPER XL_Currency)
{
	ADD_LOG("Local_PXL_ModelParam_Create");
	bool PersistentInXL = false;
	return Local_ModelParamCommon(
	    XL_ModelParamType,
	    XL_ParamTimes,
	    XL_ParamValues,
	    XL_ModelParamName,
        XL_LowerBoundary,
        XL_UpperBoundary,
        XL_InterpolMethod,
		XL_AdviseBreakPointTimes,
		XL_Currency,
        PersistentInXL );
}




////////////////////////////////////////////////
/// very rapid definition of ourlonglongFunctor
////////////////////////////////////////////////
class CstModelParamFunc : public ARMResultLong2LongFunc
{
public:
	CstModelParamFunc(
		long modelParamType,
		double paramValue,
		bool adviseBreakPointTimes)
    :
		C_ModelParamType(modelParamType),
		C_ParamValue(paramValue),
		C_AdviseBreakPointTimes(adviseBreakPointTimes)
    {};

	long operator()( ARM_result& result, long objId )
	{
		return ARMLOCAL_CstModelParam_Create(
            C_ModelParamType,
            C_ParamValue,
			C_AdviseBreakPointTimes,
            result, objId);			
	}
private:
    long            C_ModelParamType;
    double			C_ParamValue;
	bool			C_AdviseBreakPointTimes;
};


/////////////////////////////////////////////////////////////
/// central function that does the creation of the XL function
/////////////////////////////////////////////////////////////
LPXLOPER Local_CstModelParamCommon(
	LPXLOPER XL_ModelParamType,
	LPXLOPER XL_ParamValue,
	LPXLOPER XL_AdviseBreakPointTimes,
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

		CCString C_ModelParamTypeStr;
		XL_readStrCell(XL_ModelParamType,C_ModelParamTypeStr," ARM_ERR: Param Type: String expected",C_result);
		long C_ModelParamType;
		if( (C_ModelParamType = ARM_ConvGPModelParam( C_ModelParamTypeStr, C_result)) == ARM_DEFAULT_ERR )
		{
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}

		double C_ParamValue;
		XL_readNumCell( XL_ParamValue, C_ParamValue, " ARM_ERR: value of the parameter: numeric expected",C_result);

		double C_AdviseBreakPointTimes;
		double C_AdviseBreakPointTimesDef=0;
		XL_readNumCellWD( XL_AdviseBreakPointTimes, C_AdviseBreakPointTimes, C_AdviseBreakPointTimesDef, " ARM_ERR: keep times: boolean expected",C_result);
		bool C_AdviseBreakPointTimesBool = C_AdviseBreakPointTimes != 0;


		/// use the concept of Functor to transfer the knowledge of
		/// a function with a context
		CstModelParamFunc ourFunc(
			C_ModelParamType,
			C_ParamValue,
			C_AdviseBreakPointTimesBool );

		/// call the general function
		fillXL_Result( LOCAL_MODELPARAM_CLASS, ourFunc, C_result, XL_result, PersistentInXL );
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_CstModelParamCommon" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}
							 

///////////////////////////////////
/// version that takes into account 
/// previous creation of object
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_CstModelParam_Create(
	LPXLOPER XL_ModelParamType,
	LPXLOPER XL_ParamValue,
	LPXLOPER XL_AdviseBreakPointTimes )
{
	ADD_LOG("Local_CstModelParam_Create");
	bool PersistentInXL = true;
	return Local_CstModelParamCommon(
	    XL_ModelParamType,
	    XL_ParamValue,
		XL_AdviseBreakPointTimes,
        PersistentInXL );
}


///////////////////////////////////
/// version for non persistentInXl
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_PXL_CstModelParam_Create(
	LPXLOPER XL_ModelParamType,
	LPXLOPER XL_ParamValue,
	LPXLOPER XL_AdviseBreakPointTimes )
{
	ADD_LOG("Local_PXL_CstModelParam_Create");
	bool PersistentInXL = false;
	return Local_CstModelParamCommon(
	    XL_ModelParamType,
	    XL_ParamValue,
		XL_AdviseBreakPointTimes,
        PersistentInXL );
}


///----------------------------------------------
///----------------------------------------------
///             Trigo Correl Param
/// Inputs :
///     DateStrip
///     Theta
///----------------------------------------------
///----------------------------------------------
////////////////////////////////////////////////
/// very rapid definition of ourlonglongFunctor
////////////////////////////////////////////////
class TrigoCorrelParamFunc : public ARMResultLong2LongFunc
{
public:
	TrigoCorrelParamFunc(
		double asOfDate,
		long dateStripId,
		double theta,
		const CCString& interpolatorName,
		const vector<double>& upperBound,
		const vector<double>& lowerBound )
    :
		C_asOfDate(asOfDate),
		C_dateStripId(dateStripId),
		C_theta(theta),
		C_interpolatorName(interpolatorName),
		C_upperBound(upperBound),
		C_lowerBound(lowerBound)
    {};

	long operator()( ARM_result& result, long objId )
	{
		return ARMLOCAL_TrigoCorrelParam_Create(
			C_asOfDate,
            C_dateStripId,
            C_theta,
			C_interpolatorName,
			C_upperBound,
			C_lowerBound,
            result, 
			objId);			
	}

private:
	double			C_asOfDate;
    long            C_dateStripId;
    double			C_theta;
	CCString		C_interpolatorName;
	vector<double>  C_upperBound;
	vector<double>  C_lowerBound;
};


/////////////////////////////////////////////////////////////
/// central function that does the creation of the XL function
/////////////////////////////////////////////////////////////
LPXLOPER Local_TrigoCorrelParamCommon(
	LPXLOPER XL_AsOfDate,
	LPXLOPER XL_DateStripId,
	LPXLOPER XL_Theta,
    LPXLOPER XL_lowerBound,
    LPXLOPER XL_upperBound,
	LPXLOPER XL_interpolatorName,
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

		double C_asOfDate;
		XL_readNumCell( XL_AsOfDate, C_asOfDate, " ARM_ERR: as Of Date: date expected",C_result);

		long C_dateStripId;
		XL_GETOBJID( XL_DateStripId, C_dateStripId,	" ARM_ERR: date strip: Object expected",	C_result);

		double C_theta;
		XL_readNumCell( XL_Theta, C_theta, " ARM_ERR: value of the parameter: numeric expected",C_result);

		vector<double> defaultVector(0);
		vector<double> C_LowerBoundary;
		vector<double> C_UpperBoundary;
		XL_readNumVectorWD(XL_lowerBound,C_LowerBoundary,defaultVector," ARM_ERR: lower bound: array of numeric expected",C_result);
		if(C_LowerBoundary.size() > 0)
		{
			XL_readNumVector(XL_upperBound,C_UpperBoundary," ARM_ERR: Upper bound: array of numeric expected",C_result);
		}
		else 
			C_UpperBoundary.clear(); 
        
        CCString C_InterpolMethod;
        CCString default_InterpolMethod = "STEPUPRIGHT";
		XL_readStrCellWD(XL_interpolatorName, C_InterpolMethod, default_InterpolMethod, " ARM_ERR: InterpolMethod : numeric expected",C_result);	  


		/// use the concept of Functor to transfer the knowledge of
		/// a function with a context
		TrigoCorrelParamFunc ourFunc(
			C_asOfDate,
			C_dateStripId,
			C_theta,
			C_InterpolMethod,
			C_LowerBoundary,
			C_UpperBoundary);

		/// call the general function
		fillXL_Result( LOCAL_MODELPARAM_CLASS, ourFunc, C_result, XL_result, PersistentInXL );
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_TrigoCorrelParamCommon" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}
							 

///////////////////////////////////
/// version that takes into account 
/// previous creation of object
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_TrigoCorrelParam_Create(
	LPXLOPER XL_AsOfDate,
	LPXLOPER XL_DateStripId,
	LPXLOPER XL_Theta,
    LPXLOPER XL_lowerBound,
    LPXLOPER XL_upperBound,
	LPXLOPER XL_interpolatorName )
{
	ADD_LOG("Local_TrigoCorrelParam_Create");
	bool PersistentInXL = true;
	return Local_TrigoCorrelParamCommon(
		XL_AsOfDate,
	    XL_DateStripId,
	    XL_Theta,
		XL_lowerBound,
		XL_upperBound,
		XL_interpolatorName,
        PersistentInXL );
}


///////////////////////////////////
/// version for non persistentInXl
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_PXL_TrigoCorrelParam_Create(
	LPXLOPER XL_AsOfDate,
	LPXLOPER XL_DateStripId,
	LPXLOPER XL_Theta,
    LPXLOPER XL_lowerBound,
    LPXLOPER XL_upperBound,
	LPXLOPER XL_interpolatorName )
{
	ADD_LOG("Local_PXL_TrigoCorrelParam_Create");
	bool PersistentInXL = false;
	return Local_TrigoCorrelParamCommon(
		XL_AsOfDate,
	    XL_DateStripId,
	    XL_Theta,
		XL_lowerBound,
		XL_upperBound,
		XL_interpolatorName,
        PersistentInXL );
}


///----------------------------------------------
///----------------------------------------------
///             HW1F Model
/// Inputs :
///     Zc curve Id
///     Sigma param Id
///     Mean Reversion param Id
///----------------------------------------------
///----------------------------------------------
////////////////////////////////////////////////
/// very rapid definition of ourlonglongFunctor
////////////////////////////////////////////////
class HW1FModelFunc : public ARMResultLong2LongFunc
{
public:
	HW1FModelFunc(
        long zeroCurveId,
        long sigmaParamId,
        long meanReversionParamId,
		vector<double> flags)
    :
    C_zeroCurveId(zeroCurveId),
    C_SigmaParamId(sigmaParamId),
    C_MeanReversionParamId(meanReversionParamId),
	C_Flags(flags)
    {};
	
	long operator()( ARM_result& result, long objId ){
		return ARMLOCAL_HW1FModel_Create(
            C_zeroCurveId,
            C_SigmaParamId,
            C_MeanReversionParamId,
			C_Flags,
            result,
            objId);			
	}

private:
	long C_zeroCurveId;
	long C_SigmaParamId;
	long C_MeanReversionParamId;
	vector< double > C_Flags;
};


/////////////////////////////////////////////////////////////
/// central function that does the creation of the XL function
/////////////////////////////////////////////////////////////
LPXLOPER Local_HW1FModelCommon(
	LPXLOPER XL_ZeroCurveId,
	LPXLOPER XL_SigmaParamId,
	LPXLOPER XL_MeanReversionParamId,
	LPXLOPER XL_Flags,
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

		long C_ZeroCurveId;
		XL_GETOBJID( XL_ZeroCurveId, C_ZeroCurveId,	" ARM_ERR: Zero Curve: Object expected",C_result);

		long C_SigmaParamId;
		XL_GETOBJID( XL_SigmaParamId, C_SigmaParamId,	" ARM_ERR: Model Parameter: Object expected",C_result);

		long C_MeanReversionParamId;
		XL_GETOBJID( XL_MeanReversionParamId, C_MeanReversionParamId, " ARM_ERR: Model Parameter: Object expected",C_result );
		
		/// Flags
		VECTOR<double> C_Flags;
		VECTOR<double> flagsDefault(2,1.0);
		XL_readNumVectorWD(XL_Flags,C_Flags,flagsDefault," ARM_ERR: flags: array of double expected",C_result);

		/// use the concept of Functor to transfer the knowledge of
		/// a function with a context
		HW1FModelFunc ourFunc(
			C_ZeroCurveId,
			C_SigmaParamId,
			C_MeanReversionParamId,
			C_Flags);

		/// call the general function
		fillXL_Result( LOCAL_HW1FMOD_CLASS, ourFunc, C_result, XL_result, PersistentInXL );
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_HW1FModelCommon" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}
							 

///////////////////////////////////
/// version that takes into account 
/// previous creation of object
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_HW1FModel_Create(
	LPXLOPER XL_ZeroCurveId,
	LPXLOPER XL_SigmaParamId,
	LPXLOPER XL_MeanReversionParamId,
	LPXLOPER XL_Flags)
{
	ADD_LOG("Local_HW1FModel_Create");
	bool PersistentInXL = true;
	return Local_HW1FModelCommon(
        XL_ZeroCurveId,
        XL_SigmaParamId,
        XL_MeanReversionParamId,
		XL_Flags,
        PersistentInXL );
}


///////////////////////////////////
/// version for non persistentInXl
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_PXL_HW1FModel_Create(
	LPXLOPER XL_ZeroCurveId,
	LPXLOPER XL_SigmaParamId,
	LPXLOPER XL_MeanReversionParamId,
	LPXLOPER XL_Flags)
{
	ADD_LOG("Local_PXL_HW1FModel_Create");
	bool PersistentInXL = false;
	return Local_HW1FModelCommon(
        XL_ZeroCurveId,
        XL_SigmaParamId,
        XL_MeanReversionParamId,
		XL_Flags,
        PersistentInXL );
}


///----------------------------------------------
///----------------------------------------------
///             MF1F Model
/// Inputs :
///     Zc curve Id
///     Sigma param Id
///----------------------------------------------
///----------------------------------------------
////////////////////////////////////////////////
/// very rapid definition of ourlonglongFunctor
////////////////////////////////////////////////
class MF1FModelFunc : public ARMResultLong2LongFunc
{
public:
	MF1FModelFunc(
        long zeroCurveId,
        long Param1Id,
		long Param2Id)
    :
    C_zeroCurveId(zeroCurveId),
    C_Param1Id(Param1Id),
	C_Param2Id(Param2Id)
    {};
	
	long operator()( ARM_result& result, long objId ){
		return ARMLOCAL_MarkovFunctionalModel_Create(
            C_zeroCurveId,
            C_Param1Id,
			C_Param2Id,
            result,
            objId);			
	}

private:
	long C_zeroCurveId;
	long C_Param1Id;
	long C_Param2Id;
};


/////////////////////////////////////////////////////////////
/// central function that does the creation of the XL function
/////////////////////////////////////////////////////////////
LPXLOPER Local_MF1FModelCommon(
	LPXLOPER XL_ZeroCurveId,
	LPXLOPER XL_Param1Id, /// vol
	LPXLOPER XL_Param2Id, /// MR
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

		long C_ZeroCurveId;
		XL_GETOBJID( XL_ZeroCurveId, C_ZeroCurveId,	" ARM_ERR: Zero Curve: Object expected",C_result);

		long C_Param1Id;
		XL_GETOBJID( XL_Param1Id, C_Param1Id,	" ARM_ERR: Model Parameter: Object expected",C_result);

		long C_Param2Id;
		XL_GETOBJIDWD( XL_Param2Id, C_Param2Id,	"NULL OBJECT"," ARM_ERR: Model Parameter: Object expected",C_result);

		/// use the concept of Functor to transfer the knowledge of
		/// a function with a context
		MF1FModelFunc ourFunc(
			C_ZeroCurveId,
			C_Param1Id,
			C_Param2Id);

		/// call the general function
		fillXL_Result( LOCAL_MF1FMOD_CLASS, ourFunc, C_result, XL_result, PersistentInXL );
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_MF1FModelCommon" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}

///////////////////////////////////
/// version that takes into account 
/// previous creation of object
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_MF1FModel_Create(
	LPXLOPER XL_ZeroCurveId,
	LPXLOPER XL_Param1Id,
	LPXLOPER XL_Param2Id)
{
	ADD_LOG("Local_MF1FModel_Create");
	bool PersistentInXL = true;
	return Local_MF1FModelCommon(
        XL_ZeroCurveId,
        XL_Param1Id,
		XL_Param2Id,
        PersistentInXL );
}


///////////////////////////////////
/// version for non persistentInXl
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_PXL_MF1FModel_Create(
	LPXLOPER XL_ZeroCurveId,
	LPXLOPER XL_Param1Id,
	LPXLOPER XL_Param2Id)
{
	ADD_LOG("Local_PXL_MF1FModel_Create");
	bool PersistentInXL = false;
	return Local_MF1FModelCommon(
        XL_ZeroCurveId,
        XL_Param1Id,
		XL_Param2Id,
        PersistentInXL );
}

///----------------------------------------------
///----------------------------------------------
///             HW2F Model
/// Inputs :
///     Zc curve Id
///     Sigma param Id
///     Mean Reversion param Id
///     Sigma Ratio param Id
///     Mean Reversion Spread param Id
///     Correlation param Id
///		Flags
///----------------------------------------------
///----------------------------------------------
////////////////////////////////////////////////
/// very rapid definition of ourlonglongFunctor
////////////////////////////////////////////////
class HW2FModelFunc : public ARMResultLong2LongFunc
{
public:
	HW2FModelFunc(
        long zeroCurveId,
        long sigmaParamId,
        long meanReversionParamId,
        long sigmaRatioParamId,
        long meanReversionSpreadParamId,
        long correlationParamId,
		vector<double> flags)
    :
    C_zeroCurveId(zeroCurveId),
    C_SigmaParamId(sigmaParamId),
    C_MeanReversionParamId(meanReversionParamId),
    C_SigmaRatioParamId(sigmaRatioParamId),
    C_MeanReversionSpreadParamId(meanReversionSpreadParamId),
    C_CorrelationParamId(correlationParamId),
	C_Flags(flags)
    {};
	
	long operator()( ARM_result& result, long objId ){
		return ARMLOCAL_HW2FModel_Create(
            C_zeroCurveId,
            C_SigmaParamId,
            C_MeanReversionParamId,
            C_SigmaRatioParamId,
            C_MeanReversionSpreadParamId,
            C_CorrelationParamId,
			C_Flags,
            result,
            objId);			
	}

private:
	long C_zeroCurveId;
	long C_SigmaParamId;
	long C_MeanReversionParamId;
	long C_SigmaRatioParamId;
	long C_MeanReversionSpreadParamId;
	long C_CorrelationParamId;
	vector< double > C_Flags;
};


/////////////////////////////////////////////////////////////
/// central function that does the creation of the XL function
/////////////////////////////////////////////////////////////
LPXLOPER Local_HW2FModelCommon(
	LPXLOPER XL_ZeroCurveId,
	LPXLOPER XL_SigmaParamId,
	LPXLOPER XL_MeanReversionParamId,
	LPXLOPER XL_SigmaRatioParamId,
	LPXLOPER XL_MeanReversionSpreadParamId,
	LPXLOPER XL_CorrelationParamId,
	LPXLOPER XL_Flags,
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

		long C_ZeroCurveId;
		XL_GETOBJID( XL_ZeroCurveId, C_ZeroCurveId,	" ARM_ERR: Zero Curve: Object expected",C_result);

		long C_SigmaParamId;
		XL_GETOBJID( XL_SigmaParamId, C_SigmaParamId,	" ARM_ERR: Model Parameter: Object expected",C_result);

		long C_MeanReversionParamId;
		XL_GETOBJID( XL_MeanReversionParamId, C_MeanReversionParamId, " ARM_ERR: Model Parameter: Object expected",C_result );
		
		long C_SigmaRatioParamId;
		XL_GETOBJID( XL_SigmaRatioParamId, C_SigmaRatioParamId, " ARM_ERR: Model Parameter: Object expected",C_result );
		
		long C_MeanReversionSpreadParamId;
		XL_GETOBJID( XL_MeanReversionSpreadParamId, C_MeanReversionSpreadParamId, " ARM_ERR: Model Parameter: Object expected",C_result );
		
		long C_CorrelationParamId;
		XL_GETOBJID( XL_CorrelationParamId, C_CorrelationParamId, " ARM_ERR: Model Parameter: Object expected",C_result );
		
		/// Flags
		VECTOR<double> C_Flags;
		VECTOR<double> flagsDefault(2,1.0);
		XL_readNumVectorWD(XL_Flags,C_Flags,flagsDefault," ARM_ERR: flags: array of double expected",C_result);

		/// use the concept of Functor to transfer the knowledge of
		/// a function with a context
		HW2FModelFunc ourFunc(
			C_ZeroCurveId,
			C_SigmaParamId,
			C_MeanReversionParamId,
			C_SigmaRatioParamId,
			C_MeanReversionSpreadParamId,
			C_CorrelationParamId,
			C_Flags);

		/// call the general function
		fillXL_Result( LOCAL_HW2FMOD_CLASS, ourFunc, C_result, XL_result, PersistentInXL );
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_HW2FModelCommon" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}
							 

///////////////////////////////////
/// version that takes into account 
/// previous creation of object
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_HW2FModel_Create(
	LPXLOPER XL_ZeroCurveId,
	LPXLOPER XL_SigmaParamId,
	LPXLOPER XL_MeanReversionParamId,
	LPXLOPER XL_SigmaRatioParamId,
	LPXLOPER XL_MeanReversionSpreadParamId,
	LPXLOPER XL_CorrelationParamId,
	LPXLOPER XL_Flags)
{
	ADD_LOG("Local_HW2FModel_Create");
	bool PersistentInXL = true;
	return Local_HW2FModelCommon(
        XL_ZeroCurveId,
        XL_SigmaParamId,
        XL_MeanReversionParamId,
	    XL_SigmaRatioParamId,
	    XL_MeanReversionSpreadParamId,
	    XL_CorrelationParamId,
		XL_Flags,
        PersistentInXL );
}


///////////////////////////////////
/// version for non persistentInXl
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_PXL_HW2FModel_Create(
	LPXLOPER XL_ZeroCurveId,
	LPXLOPER XL_SigmaParamId,
	LPXLOPER XL_MeanReversionParamId,
	LPXLOPER XL_SigmaRatioParamId,
	LPXLOPER XL_MeanReversionSpreadParamId,
	LPXLOPER XL_CorrelationParamId,
	LPXLOPER XL_Flags)
{
	ADD_LOG("Local_PXL_HW2FModel_Create");
	bool PersistentInXL = false;
	return Local_HW2FModelCommon(
        XL_ZeroCurveId,
        XL_SigmaParamId,
        XL_MeanReversionParamId,
	    XL_SigmaRatioParamId,
	    XL_MeanReversionSpreadParamId,
	    XL_CorrelationParamId,
		XL_Flags,
        PersistentInXL );
}

///----------------------------------------------
///----------------------------------------------
///             SFMR Model
/// Inputs :
///     Zc curve Id
///     Vector of param Id (volatility, mean reversion, shift, correlation)
///		volType
///		factorsNb
///		IRIndexId
///----------------------------------------------
///----------------------------------------------
////////////////////////////////////////////////
/// very rapid definition of ourlonglongFunctor
////////////////////////////////////////////////
class SFRMModelFunc : public ARMResultLong2LongFunc
{
public:
	SFRMModelFunc(
		long zeroCurveId,
		const vector<long >& paramsIdVec,
		long volType,
		long factorsNb,
		long IRIndexId,
		long shiftConvPortId,
		bool diffusionDrift)
    :
		C_zeroCurveId(zeroCurveId),
		C_paramsIdVec(paramsIdVec),
		C_volType(volType),
		C_factorsNb(factorsNb),
		C_IRIndexId(IRIndexId),
		C_ShiftConvPortId(shiftConvPortId),
		C_diffusionDrift(diffusionDrift)
    {};

	long operator()( ARM_result& result, long objId ){
		return ARMLOCAL_SFRMModel_Create(
            C_zeroCurveId,
            C_paramsIdVec,
            C_volType,
            C_factorsNb,
            C_IRIndexId,
			C_ShiftConvPortId,
			C_diffusionDrift,
            result,
            objId);			
	}
private:
	long C_zeroCurveId;
	vector<long> C_paramsIdVec;
	long C_volType;
	long C_factorsNb;
	long C_IRIndexId;
	long C_ShiftConvPortId;
	bool C_diffusionDrift;
};


/////////////////////////////////////////////////////////////
/// central function that does the creation of the XL function
/////////////////////////////////////////////////////////////
LPXLOPER Local_SFRMModelCommon(
	LPXLOPER XL_ZeroCurveId,
	LPXLOPER XL_ParamsIdVec,
	LPXLOPER XL_volType,
	LPXLOPER XL_factorsNb,
	LPXLOPER XL_IRIndexId,
	LPXLOPER XL_ShiftConvPortId,
	LPXLOPER XL_NonParamDrift,
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

		long C_ZeroCurveId;
		XL_GETOBJID( XL_ZeroCurveId, C_ZeroCurveId,	" ARM_ERR: Zero Curve: Object expected",C_result);
		
		vector<CCString> C_paramsIds;
		XL_readStrVector (XL_ParamsIdVec,C_paramsIds," ARM_ERR: Model Parameters: array of object expected",DOUBLE_TYPE,C_result);
		size_t i, size = C_paramsIds.size();
		vector<long> C_paramsIdVec(size);
		for(i = 0; i < size; ++i )    
			C_paramsIdVec[i] = LocalGetNumObjectId(C_paramsIds[i]); 

		long C_volType;
		XL_GETCONVSHAPETYPEWD( XL_volType, C_volType, "ROW"," ARM_ERR: Vol Type: string expected",C_result);
		
		double C_factorsNb;
		XL_readNumCell( XL_factorsNb, C_factorsNb, " ARM_ERR: Nb of factors: numeric expected",C_result);

		long C_IRIndexId;
		XL_GETOBJID( XL_IRIndexId, C_IRIndexId,	" ARM_ERR: IR Index: Object expected",C_result);
		
		long C_ShiftConvPortId;
		XL_GETOBJIDWD( XL_ShiftConvPortId, C_ShiftConvPortId, GETDEFAULTVALUESTR, " ARM_ERR: portfolio: Object expected",C_result);

		CCString NonParamDriftXl;
		XL_readStrCellWD( XL_NonParamDrift,NonParamDriftXl,"N"," ARM_ERR: Drift Diffusion: string expected",C_result);
		bool C_NonParamDrift;
		NonParamDriftXl.toUpper();
		if (NonParamDriftXl == "Y" )
			C_NonParamDrift = true;
		else if (NonParamDriftXl == "N")
			C_NonParamDrift = false;
		else
			throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"String \"Y\" or \"N\" Expected for Diffusion Drift");

		

		/// use the concept of Functor to transfer the knowledge of
		/// a function with a context
		exportFunc7Args< long, vector<long>, long, double, long, long, bool >  ourFunc(
			C_ZeroCurveId,
			C_paramsIdVec,
			C_volType,
			C_factorsNb,
			C_IRIndexId,
			C_ShiftConvPortId,
			C_NonParamDrift,
			ARMLOCAL_SFRMModel_Create);

		/// call the general function
		fillXL_Result( LOCAL_SFRMMODEL_CLASS, ourFunc, C_result, XL_result, PersistentInXL );
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_SFRMModelCommon" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


/////////////////////////////////////////////////////////////
/// central function that does the creation of the XL function
/////////////////////////////////////////////////////////////
__declspec(dllexport) LPXLOPER Local_SFRMModelVolSwapVolFRADump(
	LPXLOPER XL_SwaptionId,
	LPXLOPER XL_SFRMId)
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

		long C_SwaptionId;
		XL_GETOBJID( XL_SwaptionId, C_SwaptionId,	" ARM_ERR: Swaption: Object expected",C_result);
		long C_SFRMId;
		XL_GETOBJID( XL_SFRMId, C_SFRMId,	" ARM_ERR: SFRM: Object expected",C_result);

		VECTOR<double> C_OutputMatrix;
		long C_outputRows, C_outputCols;
		
		long retCode = ARMLOCAL_SFRMModel_VolSwapVolFRADump(
			C_SwaptionId, 
			C_SFRMId,
			C_OutputMatrix,
			C_outputRows,
			C_outputCols,
			C_result);

		if (retCode == ARM_OK)
		{
			XL_writeNumMatrixSize( XL_result, C_OutputMatrix, C_outputRows, C_outputCols, " ARM_ERR: Could not set the num matrix", C_result );
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_SFRMModelCommon" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


///////////////////////////////////
/// version that takes into account 
/// previous creation of object
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_SFRMModel_Create(
	LPXLOPER XL_ZeroCurveId,
	LPXLOPER XL_ParamsIdVec,
	LPXLOPER XL_volType,
	LPXLOPER XL_factorsNb,
	LPXLOPER XL_IRIndexId,
	LPXLOPER XL_ShiftConvPortId,
	LPXLOPER XL_NonParamDrift)
{
	ADD_LOG("Local_SFRMModel_Create");
	bool PersistentInXL = true;
	return Local_SFRMModelCommon(
        XL_ZeroCurveId,
		XL_ParamsIdVec,
		XL_volType,
		XL_factorsNb,
		XL_IRIndexId,
		XL_ShiftConvPortId,
		XL_NonParamDrift,
		PersistentInXL );
}

///////////////////////////////////
/// version for non persistentInXl
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_PXL_SFRMModel_Create(
	LPXLOPER XL_ZeroCurveId,
	LPXLOPER XL_ParamsIdVec,
	LPXLOPER XL_volType,
	LPXLOPER XL_factorsNb,
	LPXLOPER XL_IRIndexId,
	LPXLOPER XL_ShiftConvPortId,
	LPXLOPER XL_NonParamDrift)
{
	ADD_LOG("Local_PXL_SFRMModel_Create");
	bool PersistentInXL = false;
	return Local_SFRMModelCommon(
        XL_ZeroCurveId,
		XL_ParamsIdVec,
		XL_volType,
		XL_factorsNb,
		XL_IRIndexId,
		XL_ShiftConvPortId,
		XL_NonParamDrift,
		PersistentInXL );
}


///----------------------------------------------
///----------------------------------------------
///             Set SFRM Model Fix Scheduler
/// Inputs :
///     SFRM Model Id
///     Date Strip Id
///----------------------------------------------
///----------------------------------------------
////////////////////////////////////////////////
/// very rapid definition of ourlonglongFunctor
////////////////////////////////////////////////

/////////////////////////////////////////////////////////////
/// central function that does the creation of the XL function
/////////////////////////////////////////////////////////////
LPXLOPER Local_SetSFRMFixSchedulerCommon(
	LPXLOPER XL_SFRMModId,
	LPXLOPER XL_StartDate,
	LPXLOPER XL_EndDate,
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

		long C_SFRMModId;
		XL_GETOBJID( XL_SFRMModId, C_SFRMModId,	" ARM_ERR: SFRM Model: Object expected",C_result);

		double C_StartDate,C_StartDateDef = 0.0;
		XL_readNumCellWD( XL_StartDate, C_StartDate, C_StartDateDef,	" ARM_ERR: Start Date: Date expected",C_result);

		double C_EndDate,C_EndDateDef = 0.0;
		XL_readNumCellWD( XL_EndDate, C_EndDate, C_EndDateDef,	" ARM_ERR: End Date: Date expected",C_result);
		

		/// use the concept of Functor to transfer the knowledge of
		/// a function with a context
		exportFunc3Args< long, double, double >  ourFunc( C_SFRMModId, C_StartDate, C_EndDate, ARMLOCAL_SetSFRMFixScheduler );

		/// call the general function
		fillXL_Result( LOCAL_SFRMMODEL_CLASS, ourFunc, C_result, XL_result, PersistentInXL );
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_SetSFRMSchedulerCommon" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}



///////////////////////////////////
/// version that takes into account 
/// previous creation of object
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_SetSFRMFixScheduler(
		LPXLOPER XL_SFRMModId,
		LPXLOPER XL_StartDate,
		LPXLOPER XL_EndDate)
{
	ADD_LOG("Local_SetSFRMFixScheduler");
	bool PersistentInXL = true;
	return Local_SetSFRMFixSchedulerCommon(
        XL_SFRMModId,
		XL_StartDate,
		XL_EndDate,
		PersistentInXL );
}

///////////////////////////////////
/// version for non persistentInXl
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_PXL_SetSFRMFixScheduler(
		LPXLOPER XL_SFRMModId,
		LPXLOPER XL_StartDate,
		LPXLOPER XL_EndDate)
{
	ADD_LOG("Local_PXL_SetSFRMFixScheduler");
	bool PersistentInXL = false;
	return Local_SetSFRMFixSchedulerCommon(
        XL_SFRMModId,
		XL_StartDate,
		XL_EndDate,
		PersistentInXL);
}


/////////////////////////////////////////////////////////////
/// central function that does the creation of the XL function
/////////////////////////////////////////////////////////////
LPXLOPER Local_PricingModel_GetModelParam(
	LPXLOPER XL_ModelId,
	LPXLOPER XL_ModelParamType,
	LPXLOPER XL_ModelParamDataType,
	LPXLOPER XL_ModelParamIndex,
	LPXLOPER XL_FactorNb )
{
	/// to remove memory leak put a result holder!
	static XLOPER_Holder XL_resultHolder;
	XLOPER& XL_result = XL_resultHolder.GetResult();
	
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

		CCString C_modelStrId;
		XL_readStrCell( XL_ModelId, C_modelStrId,	" ARM_ERR: Model Id: Object expected",		C_result);
		long C_modelId = LocalGetNumObjectId(C_modelStrId);

		CCString C_ModelParamTypeStr;
		XL_readStrCell(XL_ModelParamType,C_ModelParamTypeStr," ARM_ERR: Param Type: String expected",C_result);
		long C_ModelParamType;
		if( (C_ModelParamType = ARM_ConvGPModelParam( C_ModelParamTypeStr, C_result)) == ARM_DEFAULT_ERR )
		{
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}

		CCString C_ModelParamDataTypeStr;
		XL_readStrCell(XL_ModelParamDataType,C_ModelParamDataTypeStr," ARM_ERR: Param Data Type: String expected",C_result);
		long C_ModelParamDataType;
		if( (C_ModelParamDataType = ARM_ConvGPModelParamDataType( C_ModelParamDataTypeStr, C_result)) == ARM_DEFAULT_ERR )
		{
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}

		double C_ModelParamIndex;
		double C_ModelParamIndexDef=0.;
		XL_readNumCellWD( XL_ModelParamIndex, C_ModelParamIndex, C_ModelParamIndexDef, " ARM_ERR: Param Index: numeric expected",C_result);

		double C_FactorNb;
		double C_FactorNbDefault=0.;
		XL_readNumCellWD( XL_FactorNb, C_FactorNb, C_FactorNbDefault, " ARM_ERR: Factor Nb: numeric expected",C_result);

		vector<double> C_DataResult;
		long C_rows,C_cols;
		if( ARMLOCAL_PricingModel_GetModelParam( C_modelId, C_ModelParamType, C_ModelParamDataType, (double) C_ModelParamIndex, (long) C_FactorNb, C_DataResult, C_rows, C_cols, C_result ) == ARM_OK )
		{
			if( C_cols == 1 )
			{
				/// add these additional lines 
				/// to display blank lines
				const int additionalLinesNb = 100;
				bool fillWithBlank = true;
				FreeCurCellContent ();
				XL_writeNumVectorWithOptions( XL_result, C_DataResult, " ARM_ERR: Could not get result data", C_result, additionalLinesNb, fillWithBlank );
			}
			else
			{
				/// add these additional lines 
				/// to display blank lines
				const int additionalLinesNb = 100;
				bool fillWithBlank = true;
				FreeCurCellContent ();
				XL_writeNumMatrixSizeWithOptions( XL_result, C_DataResult, C_rows, C_cols, " ARM_ERR: Could not get result data", C_result, additionalLinesNb, fillWithBlank );
			}
		}
		/// be ware that ARM_ERR is a macro
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_PricingModel_GetModelParam" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}
							 
///----------------------------------------------
///----------------------------------------------
///             Hybrid Basis Fwd & IR Model
/// Inputs :
///     Reference pricing model Id (=> 1st Zc curve)
///     2nd Zc curve Id
///     Basis Zc curve Id
///     Basis reference model (1st or 2nd Zc curve)
///     Forex
///----------------------------------------------
///----------------------------------------------
////////////////////////////////////////////////
/// very rapid definition of ourlonglongFunctor
////////////////////////////////////////////////
class HybridBasisFwdIRModelFunc : public ARMResultLong2LongFunc
{
public:
	HybridBasisFwdIRModelFunc(
        long refIRModelId,
        long zeroCurveId,
        long basisZcCurveId,
        long forexId,
        const vector< string >& modelNames)
    :
    C_refIRModelId(refIRModelId),
    C_zeroCurveId(zeroCurveId),
    C_basisZcCurveId(basisZcCurveId),
    C_forexId(forexId),
    C_modelNames(modelNames)
    {};
	
	long operator()( ARM_result& result, long objId ){
		return ARMLOCAL_HybridBasisFwdIRModel_Create(
            C_refIRModelId,
            C_zeroCurveId,
            C_basisZcCurveId,
            C_forexId,
            C_modelNames,
            result,
            objId);			
	}

private:
	long                C_refIRModelId;
	long                C_zeroCurveId;
	long                C_basisZcCurveId;
    long                C_forexId;
    vector< string >    C_modelNames;
};


/////////////////////////////////////////////////////////////
/// central function that does the creation of the XL function
/////////////////////////////////////////////////////////////
LPXLOPER Local_HybridBasisFwdIRModelCommon(
	LPXLOPER XL_refIRModelId,
	LPXLOPER XL_zeroCurveId,
	LPXLOPER XL_basisZcCurveId,
	LPXLOPER XL_forexId,
    LPXLOPER XL_modelNames,
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

		long C_refIRModelId;
		XL_GETOBJID( XL_refIRModelId, C_refIRModelId,	" ARM_ERR: Reference IR model : Object expected",C_result);

		long C_zeroCurveId;
		XL_GETOBJID( XL_zeroCurveId, C_zeroCurveId,	" ARM_ERR: Zc curve of the Forward margin IR model : Object expected",C_result);

		long C_basisZcCurveId;
		XL_GETOBJID( XL_basisZcCurveId, C_basisZcCurveId, " ARM_ERR: Basis swap Zc curve : Object expected",C_result );

		long C_forexId;
		XL_GETOBJID( XL_forexId, C_forexId, " ARM_ERR: Forex : Object expected",C_result );

		vector<CCString> modelNames;
		XL_readStrVector(XL_modelNames,modelNames," ARM_ERR: Model names : array of string expected",DOUBLE_TYPE,C_result);
		vector< string > C_modelNames(modelNames.size());
		for(size_t i=0;i<modelNames.size();++i)
			C_modelNames[i]=CCSTringToSTLString(modelNames[i]);

		/// use the concept of Functor to transfer the knowledge of
		/// a function with a context
		HybridBasisFwdIRModelFunc ourFunc(
			C_refIRModelId,
			C_zeroCurveId,
			C_basisZcCurveId,
            C_forexId,
            C_modelNames);

		/// call the general function
		fillXL_Result( LOCAL_BASISFWDIRMOD_CLASS, ourFunc, C_result, XL_result, PersistentInXL );
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_HybridBasisFwdIRModelCommon" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}
							 

///////////////////////////////////
/// version that takes into account 
/// previous creation of object
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_HybridBasisFwdIRModel_Create(
	LPXLOPER XL_refIRModelId,
	LPXLOPER XL_zeroCurveId,
	LPXLOPER XL_basisZcCurveId,
	LPXLOPER XL_forexId,
    LPXLOPER XL_modelNames)
{
	ADD_LOG("Local_HybridBasisFwdIRModel_Create");
	bool PersistentInXL = true;
	return Local_HybridBasisFwdIRModelCommon(
	    XL_refIRModelId,
	    XL_zeroCurveId,
	    XL_basisZcCurveId,
        XL_forexId,
        XL_modelNames,
        PersistentInXL );
}


///////////////////////////////////
/// version for non persistentInXl
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_PXL_HybridBasisFwdIRModel_Create(
	LPXLOPER XL_refIRModelId,
	LPXLOPER XL_zeroCurveId,
	LPXLOPER XL_basisZcCurveId,
	LPXLOPER XL_forexId,
    LPXLOPER XL_modelNames)
{
	ADD_LOG("Local_PXL_HybridBasisFwdIRModel_Create");
	bool PersistentInXL = false;
	return Local_HybridBasisFwdIRModelCommon(
	    XL_refIRModelId,
	    XL_zeroCurveId,
	    XL_basisZcCurveId,
        XL_forexId,
        XL_modelNames,
        PersistentInXL );
}


///----------------------------------------------
///----------------------------------------------
///             Set ZC curve in IR Model Fix Scheduler
/// Inputs :
///     Model Id
///     ZC curve
///----------------------------------------------
///----------------------------------------------


////////////////////////////////////////////////
/// ZC CURVE
/// functor to set zc curve in model
////////////////////////////////////////////////

class ZCCurveSetFunc : public ARMResultLong2LongFunc
{
public:

	ZCCurveSetFunc(
		long modelId,
		long zcCurveId)
	:
		C_modelId(modelId),
		C_zcCurveId(zcCurveId)
	{};

	long operator()( ARM_result& result, long objId = ARM_NULL_OBJECT_ID )
	{
		return ARMLOCAL_Model_SetZCCurve(
				C_modelId,
				C_zcCurveId,
				result,
				objId);
	}
			
private:
	long C_modelId;
	long C_zcCurveId;
};


/////////////////////////////////////////////////////////////
/// central function that does the creation of the XL function
/////////////////////////////////////////////////////////////
LPXLOPER Local_Model_ZCCurveSet_Common(
	LPXLOPER XL_modelId,
	LPXLOPER XL_zcCurveId,
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

		/*long C_modelId;
		XL_GETOBJID( XL_modelId, C_modelId,	" ARM_ERR: Model id: object expected",	C_result);*/

		CCString C_modelStrId;
		XL_readStrCell( XL_modelId, C_modelStrId,	" ARM_ERR: Model Id: Object expected",		C_result);
		long C_modelId = LocalGetNumObjectId(C_modelStrId);

		long C_zcCurveId;
		XL_GETOBJID( XL_zcCurveId, C_zcCurveId,	" ARM_ERR: zero curve id: object expected",	C_result);
	
		/// use the concept of Functor to transfer the knowledge of
		/// a function with a context
		ZCCurveSetFunc ourFunc(
			C_modelId,
			C_zcCurveId );

		CCString curClass = LocalGetStringObjectClassAndError(C_modelStrId);

		/// this is a market data manager
		if( curClass != CCString(" ") )
		{
			fillXL_Result( curClass, ourFunc, C_result, XL_result, PersistentInXL );
		}
		else
		{
			ARM_ERR_AND_EXIT( "Could not load the model ", C_result);
		}
		
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_Model_ZCCurveSet" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


///////////////////////////////////
/// version that takes into account 
/// previous creation of object
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_Model_ZCCurveSet(
	LPXLOPER XL_modelId,
	LPXLOPER XL_zcCurveId )
{
	ADD_LOG("Local_Model_ZCCurveSet");
	bool PersistentInXL = true;
	return Local_Model_ZCCurveSet_Common(
		XL_modelId,
		XL_zcCurveId ,
		PersistentInXL );
}


///////////////////////////////////
/// version for non persistentInXl
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_PXL_Model_ZCCurveSet(
	LPXLOPER XL_modelId,
	LPXLOPER XL_zcCurveId )
{
	ADD_LOG("Local_PXL_Model_ZCCurveSet");
	bool PersistentInXL = false;
	return Local_Model_ZCCurveSet_Common(
		XL_modelId,
		XL_zcCurveId ,
		PersistentInXL );
}

///----------------------------------------------
///----------------------------------------------
///             QGM1F Model
/// Inputs :
///     Zc curve Id
///     Sigma param Id
///     Mean Reversion param Id
///     Skew param Id
///----------------------------------------------
///----------------------------------------------
////////////////////////////////////////////////
/// very rapid definition of ourlonglongFunctor
////////////////////////////////////////////////
class QGM1FModelFunc : public ARMResultLong2LongFunc
{
public:
	QGM1FModelFunc(
        long zeroCurveId,
        long sigmaParamId,
        long meanReversionParamId,
        long skewParamId)
    :
    C_zeroCurveId(zeroCurveId),
    C_SigmaParamId(sigmaParamId),
    C_MeanReversionParamId(meanReversionParamId),
    C_SkewParamId(skewParamId)
    {};
	
	long operator()( ARM_result& result, long objId ){
		return ARMLOCAL_QGM1FModel_Create(
            C_zeroCurveId,
            C_SigmaParamId,
            C_MeanReversionParamId,
            C_SkewParamId,
            result,
            objId);			
	}

private:
	long C_zeroCurveId;
	long C_SigmaParamId;
	long C_MeanReversionParamId;
	long C_SkewParamId;
};


/////////////////////////////////////////////////////////////
/// central function that does the creation of the XL function
/////////////////////////////////////////////////////////////
LPXLOPER Local_QGM1FModelCommon(
	LPXLOPER XL_ZeroCurveId,
	LPXLOPER XL_SigmaParamId,
	LPXLOPER XL_MeanReversionParamId,
	LPXLOPER XL_SkewParamId,
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

		long C_ZeroCurveId;
		XL_GETOBJID( XL_ZeroCurveId, C_ZeroCurveId,	" ARM_ERR: Zero Curve: Object expected",C_result);

		long C_SigmaParamId;
		XL_GETOBJID( XL_SigmaParamId, C_SigmaParamId,	" ARM_ERR: Model Parameter: Object expected",C_result);

		long C_MeanReversionParamId;
		XL_GETOBJID( XL_MeanReversionParamId, C_MeanReversionParamId, " ARM_ERR: Model Parameter: Object expected",C_result );
		
		long C_SkewParamId;
		XL_GETOBJID( XL_SkewParamId, C_SkewParamId, " ARM_ERR: Model Parameter: Object expected",C_result );
		
		/// use the concept of Functor to transfer the knowledge of
		/// a function with a context
		QGM1FModelFunc ourFunc(
			C_ZeroCurveId,
			C_SigmaParamId,
			C_MeanReversionParamId,
            C_SkewParamId);

		/// call the general function
		fillXL_Result( LOCAL_QGM1FMOD_CLASS, ourFunc, C_result, XL_result, PersistentInXL );
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_QGM1FModelCommon" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}
							 

///////////////////////////////////
/// version that takes into account 
/// previous creation of object
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_QGM1FModel_Create(
	LPXLOPER XL_ZeroCurveId,
	LPXLOPER XL_SigmaParamId,
	LPXLOPER XL_MeanReversionParamId,
	LPXLOPER XL_SkewParamId)
{
	ADD_LOG("Local_QGM1FModel_Create");
	bool PersistentInXL = true;
	return Local_QGM1FModelCommon(
        XL_ZeroCurveId,
        XL_SigmaParamId,
        XL_MeanReversionParamId,
        XL_SkewParamId,
        PersistentInXL );
}


///////////////////////////////////
/// version for non persistentInXl
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_PXL_QGM1FModel_Create(
	LPXLOPER XL_ZeroCurveId,
	LPXLOPER XL_SigmaParamId,
	LPXLOPER XL_MeanReversionParamId,
	LPXLOPER XL_SkewParamId)
{
	ADD_LOG("Local_PXL_QGM1FModel_Create");
	bool PersistentInXL = false;
	return Local_QGM1FModelCommon(
        XL_ZeroCurveId,
        XL_SigmaParamId,
        XL_MeanReversionParamId,
        XL_SkewParamId,
        PersistentInXL );
}

__declspec(dllexport) LPXLOPER WINAPI Local_QGM1F_Test(
	LPXLOPER XL_QGM1FId,LPXLOPER XL_t,LPXLOPER XL_T,LPXLOPER XL_Tn,LPXLOPER XL_Xt,
    LPXLOPER XL_Ts,LPXLOPER XL_Te,LPXLOPER XL_K,LPXLOPER XL_CapPayFloorRec,
    LPXLOPER XL_Tp,LPXLOPER XL_YF)
{
	ADD_LOG("Local_QGM1F_Test");
	int i;

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

		long QGM1FId;
		XL_GETOBJID(XL_QGM1FId,QGM1FId," ARM_ERR: QGM1F model id: object expected",C_result);

		double t,T,Tn,Xt,Ts,Te,K,CapPayFloorRec;

		XL_readNumCell( XL_t, t, " ARM_ERR: time: numeric expected",C_result);
		XL_readNumCell( XL_T, T, " ARM_ERR: maturity or expiry: numeric expected",C_result);
		XL_readNumCell( XL_Tn, Tn, " ARM_ERR: numeraire or pay time: numeric expected",C_result);
		XL_readNumCell( XL_Xt, Xt, " ARM_ERR: X(time): numeric expected",C_result);
		XL_readNumCell( XL_Ts, Ts, " ARM_ERR: start: numeric expected",C_result);
		XL_readNumCell( XL_Te, Te, " ARM_ERR: end: numeric expected",C_result);
		XL_readNumCell( XL_K, K, " ARM_ERR: strike: numeric expected",C_result);
		XL_readNumCell( XL_CapPayFloorRec, CapPayFloorRec, " ARM_ERR: CapOrPay(1)/FloorOrRec(-1): numeric expected",C_result);

		vector<double> C_Tp,C_YF;
		XL_readNumVector (XL_Tp,C_Tp," ARM_ERR: fixed leg payment dates: array of dates expected",C_result);
        ARM_GP_Vector Tp(C_Tp.size());
        for(i=0;i<Tp.size();++i) Tp[i]=C_Tp[i];

		XL_readNumVector (XL_YF,C_YF," ARM_ERR: fixed flow interest term : array of double expected",C_result);
        ARM_GP_Vector YF(C_YF.size());
        for(i=0;i<YF.size();++i) YF[i]=C_YF[i];

		long retCode = ARMLOCAL_QGM1F_Test(
				QGM1FId,t,T,Tn,Xt,Ts,Te,K,static_cast<int>(CapPayFloorRec),
                Tp,YF,C_result );

        if ( retCode == ARM_OK )
        {
	        LPXLOPER pxArray;
	        int nbrows = (int) (C_result.getDouble());
	        int nbcolumns = 1;

	        FreeCurCellErr ();

	        XL_result.xltype = xltypeMulti;
	        XL_result.val.array.columns = nbcolumns;
	        XL_result.val.array.rows = nbrows;
	        XL_result.val.array.lparray = pxArray = (LPXLOPER)GlobalAlloc (GMEM_ZEROINIT, nbrows * nbcolumns * sizeof (XLOPER));

	        for (i = 0; i < nbrows; i++)
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_QGM1F_Test" )

	
	return (LPXLOPER)&XL_result;
}

/////YKKKKKKKKKKK
///----------------------------------------------
///----------------------------------------------
///             QGM2F Model
/// Inputs :
///     Zc curve Id
///     Sigma param Id
///     Mean Reversion param Id
///     Skew param Id
///----------------------------------------------
///----------------------------------------------
////////////////////////////////////////////////
/// very rapid definition of ourlonglongFunctor
////////////////////////////////////////////////
class QGM2FModelFunc : public ARMResultLong2LongFunc
{
public:
	QGM2FModelFunc(
        long zeroCurveId,
        const vector<long >& paramsVecFactor1Id,
		const vector<long >& paramsVecFactor2Id)
    :
    C_zeroCurveId(zeroCurveId),
    C_ParamsVecFactor1Id(paramsVecFactor1Id),
    C_ParamsVecFactor2Id(paramsVecFactor2Id)
    {};
	
	long operator()( ARM_result& result, long objId ){
		return ARMLOCAL_QGM2FModel_Create(
            C_zeroCurveId,
            C_ParamsVecFactor1Id,
            C_ParamsVecFactor2Id,
            result,
            objId);			
	}

private:
	long C_zeroCurveId;
	vector<long > C_ParamsVecFactor1Id;
	vector<long > C_ParamsVecFactor2Id;
};


/////////////////////////////////////////////////////////////
/// central function that does the creation of the XL function
/////////////////////////////////////////////////////////////
LPXLOPER Local_QGM2FModelCommon(
	LPXLOPER XL_ZeroCurveId,
	LPXLOPER XL_ParamsVecFactor1Id,
	LPXLOPER XL_ParamsVecFactor2Id,
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

		long C_ZeroCurveId;
		XL_GETOBJID( XL_ZeroCurveId, C_ZeroCurveId,	" ARM_ERR: Zero Curve: Object expected",C_result);

		
		vector<CCString> C_ParamsVecFactor1Id;
		XL_readStrVector (XL_ParamsVecFactor1Id,C_ParamsVecFactor1Id," ARM_ERR: First Model Parameters: array of object expected",DOUBLE_TYPE,C_result);
		size_t i, size1 = C_ParamsVecFactor1Id.size();
		vector<long> C_ParamsVecFactor1Id_Vec(size1);
		for(i = 0; i < size1; ++i )    
			C_ParamsVecFactor1Id_Vec[i] = LocalGetNumObjectId(C_ParamsVecFactor1Id[i]); 

		vector<CCString> C_ParamsVecFactor2Id;
		XL_readStrVector (XL_ParamsVecFactor2Id,C_ParamsVecFactor2Id," ARM_ERR: Second Model Parameters: array of object expected",DOUBLE_TYPE,C_result);
		size_t size2 = C_ParamsVecFactor2Id.size();
		vector<long> C_ParamsVecFactor2Id_Vec(size2);
		for(i = 0; i < size2; ++i )    
			C_ParamsVecFactor2Id_Vec[i] = LocalGetNumObjectId(C_ParamsVecFactor2Id[i]); 
		
				
		/// use the concept of Functor to transfer the knowledge of
		/// a function with a context
		QGM2FModelFunc ourFunc(
			C_ZeroCurveId,
			C_ParamsVecFactor1Id_Vec,
			C_ParamsVecFactor2Id_Vec);

		/// call the general function
		fillXL_Result( LOCAL_QGM2FMOD_CLASS, ourFunc, C_result, XL_result, PersistentInXL );
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_QGM2FModelCommon" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}
							 

///////////////////////////////////
/// version that takes into account 
/// previous creation of object
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_QGM2FModel_Create(
	LPXLOPER XL_ZeroCurveId,
	LPXLOPER XL_ParamsVecFactor1Id,
	LPXLOPER XL_ParamsVecFactor2Id)
{
	ADD_LOG("Local_QGM2FModel_Create");
	bool PersistentInXL = true;
	return Local_QGM2FModelCommon(
        XL_ZeroCurveId,
        XL_ParamsVecFactor1Id,
        XL_ParamsVecFactor2Id,
        PersistentInXL );
}


///////////////////////////////////
/// version for non persistentInXl
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_PXL_QGM2FModel_Create(
	LPXLOPER XL_ZeroCurveId,
	LPXLOPER XL_ParamsVecFactor1Id,
	LPXLOPER XL_ParamsVecFactor2Id)
{
	ADD_LOG("Local_PXL_QGM2FModel_Create");
	bool PersistentInXL = false;
	return Local_QGM2FModelCommon(
        XL_ZeroCurveId,
        XL_ParamsVecFactor1Id,
		XL_ParamsVecFactor2Id,
        PersistentInXL );
}


///----------------------------------------------
///----------------------------------------------
///             QModel 1F Analytic Model
/// Inputs :
///     Zc curve Id
///     Sigma param Id
///     Q Model param Id
///----------------------------------------------
///----------------------------------------------
////////////////////////////////////////////////
/// very rapid definition of ourlonglongFunctor
////////////////////////////////////////////////
class Q1FAnaModelFunc : public ARMResultLong2LongFunc
{
public:
	Q1FAnaModelFunc(
        long zeroCurveId,
        long sigmaParamId,
        long QParamId)
    :
    C_zeroCurveId(zeroCurveId),
    C_SigmaParamId(sigmaParamId),
    C_QParamId(QParamId)
    {};
	
	long operator()( ARM_result& result, long objId ){
		return ARMLOCAL_Q1FAnaModel_Create(
            C_zeroCurveId,
            C_SigmaParamId,
            C_QParamId,
            result,
            objId);			
	}

private:
	long C_zeroCurveId;
	long C_SigmaParamId;
	long C_QParamId;
};


/////////////////////////////////////////////////////////////
/// central function that does the creation of the XL function
/////////////////////////////////////////////////////////////
LPXLOPER Local_QModel1FAnaCommon(
	LPXLOPER XL_ZeroCurveId,
	LPXLOPER XL_SigmaParamId,
	LPXLOPER XL_QParamId,
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

		long C_ZeroCurveId;
		XL_GETOBJID( XL_ZeroCurveId, C_ZeroCurveId,	" ARM_ERR: Zero Curve: Object expected",C_result);

		long C_SigmaParamId;
		XL_GETOBJID( XL_SigmaParamId, C_SigmaParamId,	" ARM_ERR: Model Parameter: Object expected",C_result);

		long C_QParamId;
		XL_GETOBJID( XL_QParamId, C_QParamId, " ARM_ERR: Q Model Parameter: Object expected",C_result );
		
		
		/// use the concept of Functor to transfer the knowledge of
		/// a function with a context
		Q1FAnaModelFunc ourFunc(
			C_ZeroCurveId,
			C_SigmaParamId,
			C_QParamId);

		/// call the general function
		fillXL_Result( LOCAL_Q1FMOD_CLASS, ourFunc, C_result, XL_result, PersistentInXL );
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_QModel1FCommon" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}
							 

///////////////////////////////////
/// version that takes into account 
/// previous creation of object
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_QModel1FAna_Create(
	LPXLOPER XL_ZeroCurveId,
	LPXLOPER XL_SigmaParamId,
	LPXLOPER XL_QParamId)
{
	ADD_LOG("Local_QModel1FAna_Create");
	bool PersistentInXL = true;
	return Local_QModel1FAnaCommon(
        XL_ZeroCurveId,
        XL_SigmaParamId,
        XL_QParamId,
        PersistentInXL );
}


///////////////////////////////////
/// version for non persistentInXl
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_PXL_QModel1FAna_Create(
	LPXLOPER XL_ZeroCurveId,
	LPXLOPER XL_SigmaParamId,
	LPXLOPER XL_QParamId)
{
	ADD_LOG("Local_PXL_QModel1FAna_Create");
	bool PersistentInXL = false;
	return Local_QModel1FAnaCommon(
        XL_ZeroCurveId,
        XL_SigmaParamId,
        XL_QParamId,
        PersistentInXL );
}







///----------------------------------------------
///----------------------------------------------
///             QModel 1F Model
/// Inputs :
///     Zc curve Id
///     Sigma param Id
///     Q Model param Id
///----------------------------------------------
///----------------------------------------------
////////////////////////////////////////////////
/// very rapid definition of ourlonglongFunctor
////////////////////////////////////////////////
class Q1FModelFunc : public ARMResultLong2LongFunc
{
public:
	Q1FModelFunc(
        long zeroCurveId,
		const vector<long >& paramsIdVec,
		bool DegenerateInHW )
    :
    C_zeroCurveId(zeroCurveId),
    C_paramsIdVec(paramsIdVec),
	C_DegenerateInHW( DegenerateInHW )
    {};
	
	long operator()( ARM_result& result, long objId ){
		return ARMLOCAL_Q1FModel_Create(
            C_zeroCurveId,
            C_paramsIdVec,
			C_DegenerateInHW,
            result,
            objId);			
	}

private:
	long			C_zeroCurveId;
	vector<long >	C_paramsIdVec;
	bool			C_DegenerateInHW;
};


/////////////////////////////////////////////////////////////
/// central function that does the creation of the XL function
/////////////////////////////////////////////////////////////
LPXLOPER Local_QModel1FCommon(
	LPXLOPER XL_ZeroCurveId,
	LPXLOPER XL_ParamsIdVec,
	LPXLOPER XL_DegenerateInHW,
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

		long C_ZeroCurveId;
		XL_GETOBJID( XL_ZeroCurveId, C_ZeroCurveId,	" ARM_ERR: Zero Curve: Object expected",C_result);

		vector<CCString> C_paramsIds;
		XL_readStrVector (XL_ParamsIdVec,C_paramsIds," ARM_ERR: Model Parameters: array of object expected",DOUBLE_TYPE,C_result);
		size_t i, size = C_paramsIds.size();
		vector<long> C_paramsIdVec(size);
		for(i = 0; i < size; ++i )    
			C_paramsIdVec[i] = LocalGetNumObjectId(C_paramsIds[i]); 

		double DegenerateInHWDble;
		double DegenerateInHWDefault = false;
		XL_DegenerateInHW;
		XL_readNumCellWD( XL_DegenerateInHW, DegenerateInHWDble, DegenerateInHWDefault, " ARM_ERR: degenerated in HW: boolean compatible expected",C_result);
		bool C_DegenerateInHW = DegenerateInHWDble != 0;
		
		/// use the concept of Functor to transfer the knowledge of
		/// a function with a context
		Q1FModelFunc ourFunc(
			C_ZeroCurveId,
			C_paramsIdVec,
			C_DegenerateInHW );


		/// call the general function
		fillXL_Result( LOCAL_Q1FMOD_CLASS, ourFunc, C_result, XL_result, PersistentInXL );
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_QModel1FCommon" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}
							 

///////////////////////////////////
/// version that takes into account 
/// previous creation of object
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_QModel1F_Create(
	LPXLOPER XL_ZeroCurveId,
	LPXLOPER XL_ParamsIdVec,
	LPXLOPER XL_DegenerateInHW )
{
	ADD_LOG("Local_QModel1F_Create");
	bool PersistentInXL = true;
	return Local_QModel1FCommon(
        XL_ZeroCurveId,
        XL_ParamsIdVec,
		XL_DegenerateInHW,
        PersistentInXL );
}


///////////////////////////////////
/// version for non persistentInXl
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_PXL_QModel1F_Create(
	LPXLOPER XL_ZeroCurveId,
	LPXLOPER XL_ParamsIdVec,
	LPXLOPER XL_DegenerateInHW )
{
	ADD_LOG("Local_PXL_QModel1F_Create");
	bool PersistentInXL = false;
	return Local_QModel1FCommon(
        XL_ZeroCurveId,
        XL_ParamsIdVec,
		XL_DegenerateInHW,
        PersistentInXL );
}


///----------------------------------------------
///----------------------------------------------
///             Surface Model Param
/// Inputs :
///		modelParamType,
///		surfaceId,
///		modelParamName
///----------------------------------------------
///----------------------------------------------

LPXLOPER Local_SurfaceModelParam_Create_Common(
	LPXLOPER XL_ModelParamType,
	LPXLOPER XL_SurfaceId,
	LPXLOPER XL_ModelParamName,
	LPXLOPER XL_LowerBoundary,
    LPXLOPER XL_UpperBoundary,
	LPXLOPER XL_AdviseBreakPointTimes,
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

		CCString C_ModelParamTypeStr;
		XL_readStrCell(XL_ModelParamType,C_ModelParamTypeStr," ARM_ERR: Param Type: String expected",C_result);
		long C_ModelParamType;
		if( (C_ModelParamType = ARM_ConvGPModelParam( C_ModelParamTypeStr, C_result)) == ARM_DEFAULT_ERR )
		{
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}

		CCString XL_surfaceIdString;
		XL_readStrCell( XL_SurfaceId, XL_surfaceIdString, " ARM_ERR: Curve Id: Object expected",C_result);
		long C_surfaceId = LocalGetNumObjectId(XL_surfaceIdString);

		CCString C_ModelParamName;
		CCString defaultModelParamName="";
		XL_readStrCellWD(XL_ModelParamName,C_ModelParamName,defaultModelParamName," ARM_ERR: Param Name: String expected",C_result);

		double C_LowerBoundary;
		double defaultC_LowerBoundary=-ARM_NumericConstants::ARM_BIGGEST_POSITIVE_NUMBER;
		XL_readNumCellWD( XL_LowerBoundary, C_LowerBoundary,defaultC_LowerBoundary, " ARM_ERR: LowerBoundary: numeric expected",C_result);

		double C_UpperBoundary;
		double defaultC_UpperBoundary=ARM_NumericConstants::ARM_BIGGEST_POSITIVE_NUMBER;
		XL_readNumCellWD( XL_UpperBoundary, C_UpperBoundary,defaultC_UpperBoundary, " ARM_ERR: UpperBoundary: numeric expected",C_result);

		double C_AdviseBreakPointTimes;
		double C_AdviseBreakPointTimesDef=0;
		XL_readNumCellWD( XL_AdviseBreakPointTimes, C_AdviseBreakPointTimes, C_AdviseBreakPointTimesDef, " ARM_ERR: keep times: boolean expected",C_result);
		bool C_AdviseBreakPointTimesBool = C_AdviseBreakPointTimes != 0;


		/// a function with a context
		exportFunc6Args< long, long,double,double, CCString, bool >  ourFunc(C_ModelParamType, C_surfaceId,C_LowerBoundary,C_UpperBoundary,C_ModelParamName, C_AdviseBreakPointTimesBool, ARMLOCAL_SurfaceParam_Create );

		/// call the general function
		fillXL_Result( LOCAL_SURFACE_MODEL_PARAM_CLASS, ourFunc, C_result, XL_result, PersistentInXL );  /// fix fix fix
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_SurfaceModelParam_Create_Common" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}



///////////////////////////////////
/// version for persistentInXl
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_SurfaceModelParam_Create(
	LPXLOPER XL_ModelParamType,
	LPXLOPER XL_SurfaceId,
	LPXLOPER XL_ModelParamName,
	LPXLOPER XL_LowerBoundary,
    LPXLOPER XL_UpperBoundary,
	LPXLOPER XL_AdviseBreakPointTimes )
{
	ADD_LOG("Local_SurfaceModelParam_Create");
	bool PersistentInXL = true;
	return Local_SurfaceModelParam_Create_Common(
		XL_ModelParamType,
		XL_SurfaceId,
		XL_ModelParamName,
		XL_LowerBoundary,
		XL_UpperBoundary,
		XL_AdviseBreakPointTimes,
		PersistentInXL );
}

///////////////////////////////////
/// version for non persistentInXl
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_PXL_SurfaceModelParam_Create(
	LPXLOPER XL_ModelParamType,
	LPXLOPER XL_SurfaceId,
	LPXLOPER XL_ModelParamName,
	LPXLOPER XL_LowerBoundary,
    LPXLOPER XL_UpperBoundary,
	LPXLOPER XL_AdviseBreakPointTimes )
{
	ADD_LOG("Local_PXL_SurfaceModelParam_Create");
	bool PersistentInXL = false;
	return Local_SurfaceModelParam_Create_Common(
		XL_ModelParamType,
		XL_SurfaceId,
		XL_ModelParamName,
		XL_LowerBoundary,
		XL_UpperBoundary,
		XL_AdviseBreakPointTimes,
		PersistentInXL );
}



///----------------------------------------------
///----------------------------------------------
///             Surface List Model Param
/// Inputs :
///		modelParamType,
///		index,      // array of Index
///		surfaceIds, // array of SurfaceId
///		modelParamName
///----------------------------------------------
///----------------------------------------------

LPXLOPER Local_SurfaceListModelParam_Create_Common(
	LPXLOPER XL_ModelParamType,
	LPXLOPER XL_Index,
	LPXLOPER XL_SurfaceIdList,
	LPXLOPER XL_ModelParamName,
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

		CCString C_ModelParamTypeStr;
		XL_readStrCell(XL_ModelParamType,C_ModelParamTypeStr," ARM_ERR: Param Type: String expected",C_result);
		long C_ModelParamType;
		if( (C_ModelParamType = ARM_ConvGPModelParam( C_ModelParamTypeStr, C_result)) == ARM_DEFAULT_ERR )
		{
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}

		vector<double> C_Index;
		XL_readNumVector(XL_Index,C_Index," ARM_ERR: Index : array of numeric expected",C_result);

		// CCString XL_surfaceIdString;
		// XL_readStrCell( XL_SurfaceId, XL_surfaceIdString, " ARM_ERR: Curve Id: Object expected",C_result);
		// long C_surfaceId = LocalGetNumObjectId(XL_surfaceIdString);

		vector<CCString> C_SurfaceListStrId;
		XL_readStrVector (XL_SurfaceIdList,C_SurfaceListStrId," ARM_ERR: SurfaceParam List: array of object id expected",DOUBLE_TYPE,C_result);  // type to checked
    
		size_t size = C_SurfaceListStrId.size();
		vector<long> C_SurfaceListId;
		C_SurfaceListId.resize(size);
		size_t i;
		for(i = 0; i < size; ++i )    
			C_SurfaceListId[i] = LocalGetNumObjectId(C_SurfaceListStrId[i]); 


		CCString C_ModelParamName;
		CCString defaultModelParamName="";
		XL_readStrCellWD(XL_ModelParamName,C_ModelParamName,defaultModelParamName," ARM_ERR: Param Name: String expected",C_result);

		/// a function with a context
		exportFunc4Args< long, vector<double>, vector<long>, CCString >  ourFunc(C_ModelParamType, C_Index, C_SurfaceListId, C_ModelParamName, ARMLOCAL_SurfaceListParam_Create );

		/// call the general function
		fillXL_Result( LOCAL_SURFACE_MODEL_PARAM_CLASS, ourFunc, C_result, XL_result, PersistentInXL );  /// fix fix fix
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_SurfaceModelParam_Create_Common" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}



///////////////////////////////////
/// version for persistentInXl
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_SurfaceListModelParam_Create(
	LPXLOPER XL_ModelParamType,
	LPXLOPER XL_Index,
	LPXLOPER XL_SurfaceListId,
	LPXLOPER XL_ModelParamName )
{
	ADD_LOG("Local_SurfaceListModelParam_Create");
	bool PersistentInXL = true;
	return Local_SurfaceListModelParam_Create_Common(
		XL_ModelParamType,
		XL_Index,
		XL_SurfaceListId,
		XL_ModelParamName,
		PersistentInXL );
}

///////////////////////////////////
/// version for non persistentInXl
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_PXL_SurfaceListModelParam_Create(
	LPXLOPER XL_ModelParamType,
	LPXLOPER XL_Index,
	LPXLOPER XL_SurfaceListId,
	LPXLOPER XL_ModelParamName )
{
	ADD_LOG("Local_PXL_SurfaceListModelParam_Create");
	bool PersistentInXL = false;
	return Local_SurfaceListModelParam_Create_Common(
		XL_ModelParamType,
		XL_Index,
		XL_SurfaceListId,
		XL_ModelParamName,
		PersistentInXL );
}



///----------------------------------------------
///----------------------------------------------
///             GHeston Model
/// Inputs :
///     Zc curve Id
///     Vector of model params Id
///----------------------------------------------
///----------------------------------------------

LPXLOPER Local_Heston_Model_Create_Common(
	LPXLOPER XL_ZeroCurveId,
	LPXLOPER XL_ParamsIdVec,
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

		long C_ZeroCurveId;
		XL_GETOBJID( XL_ZeroCurveId, C_ZeroCurveId,	" ARM_ERR: Zero Curve: Object expected",		C_result);

		vector<CCString> C_paramsIds;
		XL_readStrVector (XL_ParamsIdVec,C_paramsIds," ARM_ERR: Model Parameters: array of object expected",DOUBLE_TYPE,C_result);
		size_t i, size = C_paramsIds.size();
		vector<long> C_paramsIdVec(size);
		for(i = 0; i < size; ++i )    
			C_paramsIdVec[i] = LocalGetNumObjectId(C_paramsIds[i]); 

		/// a function with a context
		exportFunc2Args< long, vector<long> >  ourFunc(C_ZeroCurveId, C_paramsIdVec, ARMLOCAL_Heston_Model_Create );

		/// call the general function
		fillXL_Result( LOCAL_HESTON_MODEL_CLASS, ourFunc, C_result, XL_result, PersistentInXL );
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_Heston_Model_Create_Common" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}

///////////////////////////////////
/// version for persistentInXl
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_Heston_Model_Create(
	LPXLOPER XL_ZeroCurveId,
	LPXLOPER XL_ModelParamsId )
{
	ADD_LOG("Local_Heston_Model_Create");
	bool PersistentInXL = true;
	return Local_Heston_Model_Create_Common(
        XL_ZeroCurveId,
        XL_ModelParamsId,
        PersistentInXL );
}

///////////////////////////////////
/// version for non persistentInXl
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_PXL_Heston_Model_Create(
	LPXLOPER XL_ZeroCurveId,
	LPXLOPER XL_ModelParamsId )
{
	ADD_LOG("Local_PXL_Heston_Model_Create");
	bool PersistentInXL = false;
	return Local_Heston_Model_Create_Common(
        XL_ZeroCurveId,
        XL_ModelParamsId,
        PersistentInXL );
}

///----------------------------------------------
///----------------------------------------------
///             ShiftedHeston Model
/// Inputs :
///     Zc curve Id
///     Vector of model params Id
///----------------------------------------------
///----------------------------------------------

LPXLOPER Local_ShiftedHeston_Model_Create_Common(
	LPXLOPER XL_ZeroCurveId,
	LPXLOPER XL_ParamsIdVec,
	LPXLOPER XL_IsMCPrice, 
	LPXLOPER XL_NbSteps,
	LPXLOPER XL_NbSimul,
	LPXLOPER XL_NbIntegrationSteps,
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

		long C_ZeroCurveId;
		XL_GETOBJID( XL_ZeroCurveId, C_ZeroCurveId,	" ARM_ERR: Zero Curve: Object expected",		C_result);

		vector<CCString> C_paramsIds;
		XL_readStrVector (XL_ParamsIdVec,C_paramsIds," ARM_ERR: Model Parameters: array of object expected",DOUBLE_TYPE,C_result);
		size_t i, size = C_paramsIds.size();
		vector<long> C_paramsIdVec(size);
		for(i = 0; i < size; ++i )    
			C_paramsIdVec[i] = LocalGetNumObjectId(C_paramsIds[i]); 

//////// For Monte carlo Tests
		double C_IsMCPriceDble;
		double DefaultIsMC = false;
		XL_readNumCellWD( XL_IsMCPrice, C_IsMCPriceDble, DefaultIsMC, " ARM_ERR: use MC: boolean compatible expected",C_result);
		bool C_IsMCPrice= C_IsMCPriceDble != 0;

		double C_nbSteps;
		double defaultC_nbSteps = 12;
		XL_readNumCellWD( XL_NbSteps, C_nbSteps,defaultC_nbSteps, " ARM_ERR: NbSteps: numeric expected",C_result);
		long C_nbStepsInt  = floor(C_nbSteps);

		double C_nbSimulations;
		double defaultC_nbSimulations = 100000;
		XL_readNumCellWD( XL_NbSimul, C_nbSimulations,defaultC_nbSimulations, " ARM_ERR: NbSimulations: numeric expected",C_result);
		long C_nbSimulationsInt  = floor(C_nbSimulations);

		double C_nbIntegrationSteps;
		double defaultC_nbIntegrationSteps = 140;
		XL_readNumCellWD( XL_NbIntegrationSteps, C_nbIntegrationSteps,defaultC_nbIntegrationSteps," ARM_ERR: Nb Integration Steps: numeric expected",C_result);
		long C_nbIntegrationStepsInt = floor(C_nbIntegrationSteps);
/////////////////////////////////

		/// a function with a context
		exportFunc6Args< long, vector<long>, bool, long, long, long >  ourFunc(C_ZeroCurveId, C_paramsIdVec,C_IsMCPrice, C_nbStepsInt,C_nbSimulationsInt, C_nbIntegrationStepsInt, ARMLOCAL_ShiftedHeston_Model_Create );

		/// call the general function
		fillXL_Result_withName(ourFunc, C_result, XL_result, PersistentInXL );
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_ShiftedHeston_Model_Create_Common" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}

///////////////////////////////////
/// version for persistentInXl
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_ShiftedHeston_Model_Create(
	LPXLOPER XL_ZeroCurveId,
	LPXLOPER XL_ModelParamsId,
	LPXLOPER XL_IsMCPrice, 
	LPXLOPER XL_NbSteps,
	LPXLOPER XL_NbSimul,
	LPXLOPER XL_NbIntegrationSteps)
{
	ADD_LOG("Local_ShiftedHeston_Model_Create");
	bool PersistentInXL = true;
	return Local_ShiftedHeston_Model_Create_Common(
        XL_ZeroCurveId,
        XL_ModelParamsId,
		XL_IsMCPrice,
		XL_NbSteps,
		XL_NbSimul,
		XL_NbIntegrationSteps,
        PersistentInXL );
}

///////////////////////////////////
/// version for non persistentInXl
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_PXL_ShiftedHeston_Model_Create(
	LPXLOPER XL_ZeroCurveId,
	LPXLOPER XL_ModelParamsId,
	LPXLOPER XL_IsMCPrice, 
	LPXLOPER XL_NbSteps,
	LPXLOPER XL_NbSimul,
	LPXLOPER XL_NbIntegrationSteps)
{
	ADD_LOG("Local_PXL_ShiftedHeston_Model_Create");
	bool PersistentInXL = false;
	return Local_ShiftedHeston_Model_Create_Common(
        XL_ZeroCurveId,
        XL_ModelParamsId,
		XL_IsMCPrice,
		XL_NbSteps,
		XL_NbSimul,
		XL_NbIntegrationSteps,
        PersistentInXL );
}

///----------------------------------------------
///----------------------------------------------
///             MSV Model
/// Inputs :
///     Zc curve Id
///     Vector of model params Id
///----------------------------------------------
///----------------------------------------------

LPXLOPER Local_MSV1FModel_Create_Common(
	LPXLOPER XL_ZeroCurveId,
	LPXLOPER XL_ParamsIdVec,
	LPXLOPER XL_fwdTerm,
	LPXLOPER XL_IRIndexId,
	LPXLOPER XL_FwdType,
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

		long C_ZeroCurveId;
		XL_GETOBJID( XL_ZeroCurveId, C_ZeroCurveId,	" ARM_ERR: Zero Curve: Object expected",		C_result);

		vector<CCString> C_paramsIds;
		XL_readStrVector (XL_ParamsIdVec,C_paramsIds," ARM_ERR: Model Parameters: array of object expected",DOUBLE_TYPE,C_result);
		size_t i, size = C_paramsIds.size();
		vector<long> C_paramsIdVec(size);
		for(i = 0; i < size; ++i )    
			C_paramsIdVec[i] = LocalGetNumObjectId(C_paramsIds[i]); 

		double C_fwdTerm;
		XL_readNumCell(XL_fwdTerm,C_fwdTerm," ARM_ERR: numeric expected",C_result);

		long C_IRIndexId;
		XL_GETOBJIDWD( XL_IRIndexId, C_IRIndexId,	"NULL OBJECT",	" ARM_ERR: IRIndex, Object Expected!",	C_result);

		CCString C_fwdType;
		CCString defaultFwdType = "Y";
		XL_readStrCellWD(XL_FwdType,C_fwdType,defaultFwdType," ARM_ERR: : String Y or N expected",C_result);
		string sfwdTypeXl = stringGetUpper(CCSTringToSTLString(C_fwdType));

		/// a function with a context
		exportFunc5Args< long, vector<long>, double, long , string>  ourFunc(C_ZeroCurveId, C_paramsIdVec, C_fwdTerm, C_IRIndexId, sfwdTypeXl, ARMLOCAL_MSV1FModel_Create );

		/// call the general function
		fillXL_Result_withName(ourFunc, C_result, XL_result, PersistentInXL );
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_MSV1FModel_Create_Common" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}

///////////////////////////////////
/// version for persistentInXl
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_MSV1FModel_Create(
	LPXLOPER XL_ZeroCurveId,
	LPXLOPER XL_ModelParamsId,
	LPXLOPER XL_fwdTerm,
	LPXLOPER XL_IRIndexId,
	LPXLOPER XL_FwdType)
{
	ADD_LOG("Local_MSV1FModel_Create");
	bool PersistentInXL = true;
	return Local_MSV1FModel_Create_Common(
        XL_ZeroCurveId,
        XL_ModelParamsId,
		XL_fwdTerm,
		XL_IRIndexId,
		XL_FwdType,
        PersistentInXL );
}

///////////////////////////////////
/// version for non persistentInXl
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_PXL_MSV1FModel_Create(
	LPXLOPER XL_ZeroCurveId,
	LPXLOPER XL_ModelParamsId, 
	LPXLOPER XL_fwdTerm,
	LPXLOPER XL_IRIndexId,
	LPXLOPER XL_FwdType)
{
	ADD_LOG("Local_PXL_MSV1FModel_Create");
	bool PersistentInXL = false;
	return Local_MSV1FModel_Create_Common(
        XL_ZeroCurveId,
        XL_ModelParamsId,
		XL_fwdTerm,
		XL_IRIndexId,
		XL_FwdType,
        PersistentInXL );
}

///----------------------------------------------
///----------------------------------------------
///             FRMSV Model
/// Inputs :
///     Zc curve Id
///     Vector of model params Id
///----------------------------------------------
///----------------------------------------------

LPXLOPER Local_FRMSVModel_Create_Common(
	LPXLOPER XL_ZeroCurveId,
	LPXLOPER XL_ParamsIdVec,
	LPXLOPER XL_ParamsIdVec2,
	LPXLOPER XL_IRIndexId,
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

		long C_ZeroCurveId;
		XL_GETOBJID( XL_ZeroCurveId, C_ZeroCurveId,	" ARM_ERR: Zero Curve: Object expected",		C_result);

		vector<CCString> C_paramsIds;
		XL_readStrVector (XL_ParamsIdVec,C_paramsIds," ARM_ERR: Model Parameters: array of object expected",DOUBLE_TYPE,C_result);
		size_t i, size = C_paramsIds.size();
		vector<long> C_paramsIdVec(size);
		for(i = 0; i < size; ++i )    
			C_paramsIdVec[i] = LocalGetNumObjectId(C_paramsIds[i]); 

		vector<CCString> C_paramsIds2;
		XL_readStrVector (XL_ParamsIdVec2,C_paramsIds2," ARM_ERR: Model Parameters: array of object expected",DOUBLE_TYPE,C_result);
		size = C_paramsIds2.size();
		vector<long> C_paramsIdVec2(size);
		for(i = 0; i < size; ++i )    
			C_paramsIdVec2[i] = LocalGetNumObjectId(C_paramsIds2[i]); 

		long C_IRIndexId;
		XL_GETOBJIDWD( XL_IRIndexId, C_IRIndexId,	"NULL OBJECT",	" ARM_ERR: IRIndex, Object Expected!",	C_result);


		/// a function with a context
		exportFunc4Args< long, vector<long>,vector<long>, long>  ourFunc(C_ZeroCurveId, C_paramsIdVec, C_paramsIdVec2, C_IRIndexId, ARMLOCAL_FRMSVModel_Create );

		/// call the general function
		fillXL_Result_withName(ourFunc, C_result, XL_result, PersistentInXL );
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_MSV1FModel_Create_Common" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}

///////////////////////////////////
/// version for persistentInXl
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_FRMSVModel_Create(
	LPXLOPER XL_ZeroCurveId,
	LPXLOPER XL_ParamsIdVec,
	LPXLOPER XL_ParamsIdVec2,
	LPXLOPER XL_IRIndexId)
{
	ADD_LOG("Local_FRMSVModel_Create");
	bool PersistentInXL = true;
	return Local_FRMSVModel_Create_Common(
        XL_ZeroCurveId,
        XL_ParamsIdVec,
		XL_ParamsIdVec2,
		XL_IRIndexId,
        PersistentInXL );
}

///////////////////////////////////
/// version for non persistentInXl
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_PXL_FRMSVModel_Create(
	LPXLOPER XL_ZeroCurveId,
	LPXLOPER XL_ParamsIdVec,
	LPXLOPER XL_ParamsIdVec2,
	LPXLOPER XL_IRIndexId)
{
	ADD_LOG("Local_PXL_FRMSVModel_Create");
	bool PersistentInXL = false;
	return Local_FRMSVModel_Create_Common(
        XL_ZeroCurveId,
        XL_ParamsIdVec,
		XL_ParamsIdVec2,
		XL_IRIndexId,
        PersistentInXL );
}


///----------------------------------------------
///----------------------------------------------
///             HWSV1F Model
/// Inputs :
///     Zc curve Id
///     Vector of model params Id
///----------------------------------------------
///----------------------------------------------

LPXLOPER Local_HWSV1FModel_Create_Common(
	LPXLOPER XL_ZeroCurveId,
	LPXLOPER XL_ParamsIdVec,
	LPXLOPER XL_SolverType,
	LPXLOPER XL_SolverParams,
	LPXLOPER XL_FormulaType,
	LPXLOPER XL_FormulaParams,
	LPXLOPER XL_FormulaTypeSO,
	LPXLOPER XL_FormulaParamsSO,
	LPXLOPER XL_MaxDecay,
	LPXLOPER XL_MaxDecaySO,
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

		long C_ZeroCurveId;
		XL_GETOBJID( XL_ZeroCurveId, C_ZeroCurveId,	" ARM_ERR: Zero Curve: Object expected",		C_result);

		vector<CCString> C_paramsIds;
		XL_readStrVector (XL_ParamsIdVec,C_paramsIds," ARM_ERR: Model Parameters: array of object expected",DOUBLE_TYPE,C_result);
		size_t i, size = C_paramsIds.size();
		vector<long> C_paramsIdVec(size);
		for(i = 0; i < size; ++i )    
			C_paramsIdVec[i] = LocalGetNumObjectId(C_paramsIds[i]); 

		CCString C_SolverTypeStr;
		CCString defaultSolverTypeStr = "RK5Adaptative";
		XL_readStrCellWD(XL_SolverType,C_SolverTypeStr,defaultSolverTypeStr," ARM_ERR: Solver type: String expected",C_result);
		long C_SolverType;
		if( (C_SolverType = ARM_ConvGPODESolverType( C_SolverTypeStr, C_result)) == ARM_DEFAULT_ERR )
		{
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}

		vector<double> defaultSolverParams(0);

		vector<double> C_SolverParams;
		XL_readNumVectorWD(XL_SolverParams,C_SolverParams,defaultSolverParams," ARM_ERR: Solver Params: array of numeric expected",C_result);

		CCString C_FormulaTypeStr;
		CCString defaultFormulaTypeStr = "Lewis";
		XL_readStrCellWD(XL_FormulaType,C_FormulaTypeStr,defaultFormulaTypeStr," ARM_ERR: Formula type: String expected",C_result);
		long C_FormulaType;
		if( (C_FormulaType = ARM_ConvHWSVFormulaType( C_FormulaTypeStr, C_result)) == ARM_DEFAULT_ERR )
		{
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}

		vector<double> defaultFormulaParams(0);

		vector<double> C_FormulaParams;
		XL_readNumVectorWD(XL_FormulaParams,C_FormulaParams,defaultFormulaParams," ARM_ERR: Formula Params: array of numeric expected",C_result);

		/// Formula parameters for spreadoption
		CCString C_FormulaTypeStrSO;
		CCString defaultFormulaTypeStrSO = "Heston";
		XL_readStrCellWD(XL_FormulaTypeSO,C_FormulaTypeStrSO,defaultFormulaTypeStrSO," ARM_ERR: Formula type: String expected",C_result);
		long C_FormulaTypeSO;
		if( (C_FormulaTypeSO = ARM_ConvHWSVFormulaType( C_FormulaTypeStrSO, C_result)) == ARM_DEFAULT_ERR )
		{
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}

		vector<double> defaultFormulaParamsSO(0);

		vector<double> C_FormulaParamsSO;
		XL_readNumVectorWD(XL_FormulaParamsSO,C_FormulaParamsSO,defaultFormulaParamsSO," ARM_ERR: Formula Params: array of numeric expected",C_result);

		double defaultMaxDecay=0;
		double C_MaxDecay;
		XL_readNumCellWD(XL_MaxDecay,C_MaxDecay,defaultMaxDecay," ARM_ERR: MaxDecay : numeric expected",C_result);

		double defaultMaxDecaySO=0;
		double C_MaxDecaySO;
		XL_readNumCellWD(XL_MaxDecaySO,C_MaxDecaySO,defaultMaxDecaySO," ARM_ERR: MaxDecay : numeric expected",C_result);

		/// a function with a context
		exportFunc10Args< long, vector<long>, long, vector<double>, long, vector<double>, long, vector<double>, double, double >  ourFunc(C_ZeroCurveId, C_paramsIdVec, C_SolverType, C_SolverParams, C_FormulaType, C_FormulaParams, C_FormulaTypeSO, C_FormulaParamsSO, C_MaxDecay, C_MaxDecaySO, ARMLOCAL_HWSV1FModel_Create );

		/// call the general function
		fillXL_Result_withName(ourFunc, C_result, XL_result, PersistentInXL );
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_HWSV1FModel_Create_Common" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}

///////////////////////////////////
/// version for persistentInXl
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_HWSV1FModel_Create(
	LPXLOPER XL_ZeroCurveId,
	LPXLOPER XL_ModelParamsId,
	LPXLOPER XL_SolverType,
	LPXLOPER XL_SolverParams,
	LPXLOPER XL_FormulaType,
	LPXLOPER XL_FormulaParams,
	LPXLOPER XL_MaxDecay,
	LPXLOPER XL_FormulaTypeSO,
	LPXLOPER XL_FormulaParamsSO,
	LPXLOPER XL_MaxDecaySO)
{
	ADD_LOG("Local_HWSV1FModel_Create");
	bool PersistentInXL = true;
	return Local_HWSV1FModel_Create_Common(
        XL_ZeroCurveId,
        XL_ModelParamsId,
		XL_SolverType,
		XL_SolverParams,
		XL_FormulaType,
		XL_FormulaParams,
		XL_FormulaTypeSO,
		XL_FormulaParamsSO,
		XL_MaxDecay,
		XL_MaxDecaySO,
        PersistentInXL );
}

///////////////////////////////////
/// version for non persistentInXl
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_PXL_HWSV1FModel_Create(
	LPXLOPER XL_ZeroCurveId,
	LPXLOPER XL_ModelParamsId,
	LPXLOPER XL_SolverType,
	LPXLOPER XL_SolverParams,
	LPXLOPER XL_FormulaType,
	LPXLOPER XL_FormulaParams,
	LPXLOPER XL_MaxDecay,
	LPXLOPER XL_FormulaTypeSO,
	LPXLOPER XL_FormulaParamsSO,
	LPXLOPER XL_MaxDecaySO)
{
	ADD_LOG("Local_PXL_HWSV1FModel_Create");
	bool PersistentInXL = false;
	return Local_HWSV1FModel_Create_Common(
        XL_ZeroCurveId,
        XL_ModelParamsId,
		XL_SolverType,
		XL_SolverParams,
		XL_FormulaType,
		XL_FormulaParams,
		XL_FormulaTypeSO,
		XL_FormulaParamsSO,
		XL_MaxDecay,
		XL_MaxDecaySO,
        PersistentInXL );
}


///----------------------------------------------
///----------------------------------------------
///             HWSV2F Model
/// Inputs :
///     Zc curve Id
///     Vector of model params Id
///		SolverParams
///----------------------------------------------
///----------------------------------------------

LPXLOPER Local_HWSV2FModel_Create_Common(
	LPXLOPER XL_ZeroCurveId,
	LPXLOPER XL_ParamsIdVec,
	LPXLOPER XL_SolverParams,
	LPXLOPER XL_FormulaType,
	LPXLOPER XL_FormulaParams,
	LPXLOPER XL_FormulaTypeSO,
	LPXLOPER XL_FormulaParamsSO,
	LPXLOPER XL_MaxDecay,
	LPXLOPER XL_MaxDecaySO,
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

		long C_ZeroCurveId;
		XL_GETOBJID( XL_ZeroCurveId, C_ZeroCurveId,	" ARM_ERR: Zero Curve: Object expected",		C_result);

		vector<CCString> C_paramsIds;
		XL_readStrVector (XL_ParamsIdVec,C_paramsIds," ARM_ERR: Model Parameters: array of object expected",DOUBLE_TYPE,C_result);
		size_t i, size = C_paramsIds.size();
		vector<long> C_paramsIdVec(size);
		for(i = 0; i < size; ++i )    
			C_paramsIdVec[i] = LocalGetNumObjectId(C_paramsIds[i]); 

		vector<double> defaultSolverParams(0);

		vector<double> C_SolverParams;
		XL_readNumVectorWD(XL_SolverParams,C_SolverParams,defaultSolverParams," ARM_ERR: Solver Params: array of numeric expected",C_result);

		/// Formula parameters for caplet/swaption
		CCString C_FormulaTypeStr;
		CCString defaultFormulaTypeStr = "Lewis";
		XL_readStrCellWD(XL_FormulaType,C_FormulaTypeStr,defaultFormulaTypeStr," ARM_ERR: Formula type: String expected",C_result);
		long C_FormulaType;
		if( (C_FormulaType = ARM_ConvHWSVFormulaType( C_FormulaTypeStr, C_result)) == ARM_DEFAULT_ERR )
		{
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}

		vector<double> defaultFormulaParams(0);

		vector<double> C_FormulaParams;
		XL_readNumVectorWD(XL_FormulaParams,C_FormulaParams,defaultFormulaParams," ARM_ERR: Formula Params: array of numeric expected",C_result);

		/// Formula parameters for spreadoption
		CCString C_FormulaTypeStrSO;
		CCString defaultFormulaTypeStrSO = "Heston";
		XL_readStrCellWD(XL_FormulaTypeSO,C_FormulaTypeStrSO,defaultFormulaTypeStrSO," ARM_ERR: Formula type: String expected",C_result);
		long C_FormulaTypeSO;
		if( (C_FormulaTypeSO = ARM_ConvHWSVFormulaType( C_FormulaTypeStrSO, C_result)) == ARM_DEFAULT_ERR )
		{
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}

		vector<double> defaultFormulaParamsSO(0);

		vector<double> C_FormulaParamsSO;
		XL_readNumVectorWD(XL_FormulaParamsSO,C_FormulaParamsSO,defaultFormulaParamsSO," ARM_ERR: Formula Params: array of numeric expected",C_result);

		double defaultMaxDecay=0;
		double C_MaxDecay;
		XL_readNumCellWD(XL_MaxDecay,C_MaxDecay,defaultMaxDecay," ARM_ERR: MaxDecay : numeric expected",C_result);

		double defaultMaxDecaySO=0;
		double C_MaxDecaySO;
		XL_readNumCellWD(XL_MaxDecaySO,C_MaxDecaySO,defaultMaxDecaySO," ARM_ERR: MaxDecay : numeric expected",C_result);

		/// a function with a context
		exportFunc9Args< long, vector<long>, vector<double>, long, vector<double>, long, vector<double>, double, double >  ourFunc(C_ZeroCurveId, C_paramsIdVec, C_SolverParams, C_FormulaType, C_FormulaParams, C_FormulaTypeSO, C_FormulaParamsSO, C_MaxDecay,  C_MaxDecaySO, ARMLOCAL_HWSV2FModel_Create );

		/// call the general function
		fillXL_Result_withName(ourFunc, C_result, XL_result, PersistentInXL );
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_HWSV2FModel_Create_Common" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


///////////////////////////////////
/// version for persistentInXl
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_HWSV2FModel_Create(
	LPXLOPER XL_ZeroCurveId,
	LPXLOPER XL_ModelParamsId,
	LPXLOPER XL_SolverParams,
	LPXLOPER XL_FormulaType,
	LPXLOPER XL_FormulaParams,
	LPXLOPER XL_MaxDecay,
	LPXLOPER XL_FormulaTypeSO,
	LPXLOPER XL_FormulaParamsSO,
	LPXLOPER XL_MaxDecaySO)
{
	ADD_LOG("Local_HWSV2FModel_Create");
	bool PersistentInXL = true;
	return Local_HWSV2FModel_Create_Common(
        XL_ZeroCurveId,
        XL_ModelParamsId,
		XL_SolverParams,
		XL_FormulaType,
		XL_FormulaParams,
		XL_FormulaTypeSO,
		XL_FormulaParamsSO,
		XL_MaxDecay,
		XL_MaxDecaySO,
        PersistentInXL );
}

///////////////////////////////////
/// version for non persistentInXl
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_PXL_HWSV2FModel_Create(
	LPXLOPER XL_ZeroCurveId,
	LPXLOPER XL_ModelParamsId,
	LPXLOPER XL_SolverParams,
	LPXLOPER XL_FormulaType,
	LPXLOPER XL_FormulaParams,
	LPXLOPER XL_MaxDecay,
	LPXLOPER XL_FormulaTypeSO,
	LPXLOPER XL_FormulaParamsSO,
	LPXLOPER XL_MaxDecaySO)
{
	ADD_LOG("Local_PXL_HWSV2FModel_Create");
	bool PersistentInXL = false;
	return Local_HWSV2FModel_Create_Common(
        XL_ZeroCurveId,
        XL_ModelParamsId,
		XL_SolverParams,
		XL_FormulaType,
		XL_FormulaParams,
		XL_FormulaTypeSO,
		XL_FormulaParamsSO,
		XL_MaxDecay,
		XL_MaxDecaySO,
        PersistentInXL );
}

/***********************************************

            EQHWVS Parameter Model

	Inputs :
		Vector of model params Id

***********************************************/

LPXLOPER Local_EQHWSV_ModelParamsCreate_Common(
	LPXLOPER XL_ParamsIdVec,
	bool PersistentInXL ){
	
	static XLOPER XL_result;
	ARM_result C_result;
	ARM_XL_TRY_BLOCK_BEGIN
	{
		ARM_NOCALCIFWIZ();
		static int error;
		static char* reason = "";

		vector<CCString> C_paramsIds;
		XL_readStrVector (XL_ParamsIdVec,C_paramsIds," ARM_ERR: Model Parameters: array of object expected",DOUBLE_TYPE,C_result);
		vector<long> C_paramsIdVec( C_paramsIds.size() );
		for(int i = 0; i < C_paramsIds.size(); ++i )  	C_paramsIdVec[i] = LocalGetNumObjectId(C_paramsIds[i]); 

		exportFunc1Arg< vector<long>  >  ourFunc( C_paramsIdVec,  ARMLOCAL_EQHWSV_ModelParamsCreate );
		fillXL_Result_withName(ourFunc, C_result, XL_result, PersistentInXL );
	}

	ARM_XL_TRY_BLOCK_END
	ARM_XL_CATCH_ARM_EXPT
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_EQHWVS_ParamModelCreate_Common" )
	return (LPXLOPER) &XL_result;
}

__declspec(dllexport) LPXLOPER WINAPI Local_EQHWSV_ModelParamsCreate(
	LPXLOPER XL_ModelParamsId)
{
	ADD_LOG("Local_EQHWSV_ModelParamsCreate");
	bool PersistentInXL = true;
	return Local_EQHWSV_ModelParamsCreate_Common(
        XL_ModelParamsId,
        PersistentInXL );
}

__declspec(dllexport) LPXLOPER WINAPI Local_PXL_EQHWSV_ModelParamsCreate(
	LPXLOPER XL_ModelParamsId)
{
	ADD_LOG("Local_PXL_EQHWSV_ModelParamsCreate");
	bool PersistentInXL = false;
	return Local_EQHWSV_ModelParamsCreate_Common(
        XL_ModelParamsId,
        PersistentInXL );
}


/***********************************************

            EQHWVS Numeric Method

	Inputs :
		IntStep
		ImAxis
		MaxDecay

***********************************************/

LPXLOPER Local_EQHWSV_NumMethodsCreate_Common(
	LPXLOPER XL_IntStep,
	LPXLOPER XL_ImAxis,
	LPXLOPER XL_MaxDecay,
	bool PersistentInXL )
{	
	static XLOPER XL_result;
	ARM_result C_result;
	ARM_XL_TRY_BLOCK_BEGIN
	{
		ARM_NOCALCIFWIZ();

		static int error;
		static char* reason = "";

		double defaultMaxDecay=0;
		double C_MaxDecay;
		XL_readNumCellWD(XL_MaxDecay,C_MaxDecay,defaultMaxDecay," ARM_ERR: MaxDecay : numeric expected",C_result);

		double defaultImAxis=0.5;
		double C_ImAxis;
		XL_readNumCellWD(XL_ImAxis,C_ImAxis,defaultImAxis," ARM_ERR: ImAxis : numeric expected",C_result);
		
		double defaultIntStep=20;
		double C_IntStep;
		XL_readNumCellWD(XL_IntStep,C_IntStep,defaultIntStep," ARM_ERR: IntStep : numeric expected",C_result);
		
		exportFunc3Args< long, double, double >  ourFunc((long)C_IntStep, C_ImAxis, C_MaxDecay, ARMLOCAL_EQHWSV_NumMethodsCreate );


		fillXL_Result_withName(ourFunc, C_result, XL_result, PersistentInXL );
	}

	ARM_XL_TRY_BLOCK_END
	ARM_XL_CATCH_ARM_EXPT
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_HWSV1FModel_Create_Common" )
	return (LPXLOPER)&XL_result;
}


__declspec(dllexport) LPXLOPER WINAPI Local_EQHWSV_NumMethodsCreate(
	LPXLOPER XL_IntStep,
	LPXLOPER XL_ImAxis,
	LPXLOPER XL_MaxDecay)
{
	ADD_LOG("Local_EQHWSV_NumMethodsCreate");
	bool PersistentInXL = true;
	return Local_EQHWSV_NumMethodsCreate_Common(
		XL_IntStep,
		XL_ImAxis,
		XL_MaxDecay,
        PersistentInXL );
}

__declspec(dllexport) LPXLOPER WINAPI Local_PXL_EQHWSV_NumMethodsCreate(
	LPXLOPER XL_IntStep,
	LPXLOPER XL_ImAxis,
	LPXLOPER XL_MaxDecay)
{
	ADD_LOG("Local_PXL_EQHWSV_NumMethodsCreate");
	bool PersistentInXL = false;
	return Local_EQHWSV_NumMethodsCreate_Common(
		XL_IntStep,
		XL_ImAxis,
		XL_MaxDecay,
        PersistentInXL );
}

/***********************************************

            EQHWVS

	Inputs :
		EQHWSV_ModelParams
		EQHWSV_NumMethods

***********************************************/

LPXLOPER Local_EQHWSV_Create_Common(
	LPXLOPER XL_ZeroCurveId,
	LPXLOPER XL_ModelParamsId,
	LPXLOPER XL_NumMethodsId,
	LPXLOPER XL_Dilatation,
	bool PersistentInXL ){
	
	static XLOPER XL_result;
	ARM_result C_result;
	ARM_XL_TRY_BLOCK_BEGIN
	{
		ARM_NOCALCIFWIZ();
		static int error;
		static char* reason = "";

		long C_ModelParamsId;
		XL_GETOBJID( XL_ModelParamsId, C_ModelParamsId, " ARM_ERR: EQHWSV Model Params: Object expected",	C_result);

		long C_NumMethodsId;
		XL_GETOBJID( XL_NumMethodsId, C_NumMethodsId,	" ARM_ERR: EQHWSV Num Methods: Object expected",	C_result);

		long C_ZeroCurveId;
		XL_GETOBJID( XL_ZeroCurveId, C_ZeroCurveId,	" ARM_ERR: ZeroCurve: Object expected",	C_result);

		double C_Dilatation;
		double D_Dilatation=1.0;
		XL_readNumCellWD(XL_Dilatation,C_Dilatation,D_Dilatation," ARM_ERR: numeric expected",C_result);


		exportFunc4Args<long, long, long, double >  ourFunc(C_ZeroCurveId, C_ModelParamsId, C_NumMethodsId,C_Dilatation, ARMLOCAL_EQHWSV_Create );
		fillXL_Result_withName(ourFunc, C_result, XL_result, PersistentInXL );
	}

	ARM_XL_TRY_BLOCK_END
	ARM_XL_CATCH_ARM_EXPT
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_EQHWVS_MktDatasCreate_Common" )
	return (LPXLOPER) &XL_result;
}

__declspec(dllexport) LPXLOPER WINAPI Local_EQHWSV_Create(
	LPXLOPER XL_ZeroCurveId,
	LPXLOPER XL_ModelParamsId,
	LPXLOPER XL_NumMethodsId,
	LPXLOPER XL_Dilatation)
{
	ADD_LOG("Local_EQHWSV_Create");
	bool PersistentInXL = true;
	return Local_EQHWSV_Create_Common(
		XL_ZeroCurveId,
        XL_ModelParamsId,
        XL_NumMethodsId,
		XL_Dilatation,
        PersistentInXL);
}

__declspec(dllexport) LPXLOPER WINAPI Local_PXL_EQHWSV_Create(
	LPXLOPER XL_ZeroCurveId,
	LPXLOPER XL_ModelParamsId,
	LPXLOPER XL_NumMethodsId,
	LPXLOPER XL_Dilatation	)
{
	ADD_LOG("Local_PXL_EQHWSV_Create");
	bool PersistentInXL = false;
	return Local_EQHWSV_Create_Common(
		XL_ZeroCurveId,
        XL_ModelParamsId,
        XL_NumMethodsId,
		XL_Dilatation,
        PersistentInXL);
}


///----------------------------------------------
///----------------------------------------------
///             BS Model
/// Inputs :
///     Zc curve Id
///     Vector of model params Id
///----------------------------------------------
///----------------------------------------------

LPXLOPER Local_BS_Model_Create_Common(
	LPXLOPER XL_ZeroCurveId,
	LPXLOPER XL_ParamsIdVec,
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

		long C_ZeroCurveId;
		XL_GETOBJID( XL_ZeroCurveId, C_ZeroCurveId,	" ARM_ERR: Zero Curve: Object expected",		C_result);

		vector<CCString> C_paramsIds;
		XL_readStrVector (XL_ParamsIdVec,C_paramsIds," ARM_ERR: Model Parameters: array of object expected",DOUBLE_TYPE,C_result);
		size_t i, size = C_paramsIds.size();
		vector<long> C_paramsIdVec(size);
		for(i = 0; i < size; ++i )    
			C_paramsIdVec[i] = LocalGetNumObjectId(C_paramsIds[i]); 

		/// a function with a context
		exportFunc2Args< long, vector<long> >  ourFunc(C_ZeroCurveId, C_paramsIdVec, ARMLOCAL_BS_Model_Create );

		/// call the general function
		fillXL_Result( LOCAL_BS_MODEL_CLASS, ourFunc, C_result, XL_result, PersistentInXL );
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_BS_Model_Create_Common" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


///////////////////////////////////
/// version for persistentInXl
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_BS_Model_Create(
	LPXLOPER XL_ZeroCurveId,
	LPXLOPER XL_ModelParamsId )
{
	ADD_LOG("Local_BS_Model_Create");
	bool PersistentInXL = true;
	return Local_BS_Model_Create_Common(
        XL_ZeroCurveId,
        XL_ModelParamsId,
        PersistentInXL );
}

///////////////////////////////////
/// version for non persistentInXl
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_PXL_BS_Model_Create(
	LPXLOPER XL_ZeroCurveId,
	LPXLOPER XL_ModelParamsId )
{
	ADD_LOG("Local_PXL_BS_Model_Create");
	bool PersistentInXL = false;
	return Local_BS_Model_Create_Common(
        XL_ZeroCurveId,
        XL_ModelParamsId,
        PersistentInXL );
}

///----------------------------------------------
///----------------------------------------------
///             CEV Model
/// Inputs :
///     Zc curve Id
///     Vector of model params Id
///----------------------------------------------
///----------------------------------------------

LPXLOPER Local_CEV_Model_Create_Common(
	LPXLOPER XL_ZeroCurveId,
	LPXLOPER XL_ParamsIdVec,
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

		long C_ZeroCurveId;
		XL_GETOBJID( XL_ZeroCurveId, C_ZeroCurveId,	" ARM_ERR: Zero Curve: Object expected",		C_result);

		vector<CCString> C_paramsIds;
		XL_readStrVector (XL_ParamsIdVec,C_paramsIds," ARM_ERR: Model Parameters: array of object expected",DOUBLE_TYPE,C_result);
		size_t i, size = C_paramsIds.size();
		vector<long> C_paramsIdVec(size);
		for(i = 0; i < size; ++i )    
			C_paramsIdVec[i] = LocalGetNumObjectId(C_paramsIds[i]); 

		/// a function with a context
		exportFunc2Args< long, vector<long> >  ourFunc(C_ZeroCurveId, C_paramsIdVec, ARMLOCAL_CEV_Model_Create );

		/// call the general function
		fillXL_Result( LOCAL_CEV_MODEL_CLASS, ourFunc, C_result, XL_result, PersistentInXL );
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_CEV_Model_Create_Common" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


///////////////////////////////////
/// version for persistentInXl
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_CEV_Model_Create(
	LPXLOPER XL_ZeroCurveId,
	LPXLOPER XL_ModelParamsId )
{
	ADD_LOG("Local_CEV_Model_Create");
	bool PersistentInXL = true;
	return Local_CEV_Model_Create_Common(
        XL_ZeroCurveId,
        XL_ModelParamsId,
        PersistentInXL );
}

///////////////////////////////////
/// version for non persistentInXl
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_PXL_CEV_Model_Create(
	LPXLOPER XL_ZeroCurveId,
	LPXLOPER XL_ModelParamsId )
{
	ADD_LOG("Local_PXL_CEV_Model_Create");
	bool PersistentInXL = false;
	return Local_CEV_Model_Create_Common(
        XL_ZeroCurveId,
        XL_ModelParamsId,
        PersistentInXL );
}


///----------------------------------------------
///----------------------------------------------
///             Merton Model
/// Inputs :
///     Zc curve Id
///     Vector of model params Id
///----------------------------------------------
///----------------------------------------------

LPXLOPER Local_Merton_Model_Create_Common(
	LPXLOPER XL_ZeroCurveId,
	LPXLOPER XL_ParamsIdVec,
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

		long C_ZeroCurveId;
		XL_GETOBJID( XL_ZeroCurveId, C_ZeroCurveId,	" ARM_ERR: Zero Curve: Object expected",		C_result);

		vector<CCString> C_paramsIds;
		XL_readStrVector (XL_ParamsIdVec,C_paramsIds," ARM_ERR: Model Parameters: array of object expected",DOUBLE_TYPE,C_result);
		size_t i, size = C_paramsIds.size();
		vector<long> C_paramsIdVec(size);
		for(i = 0; i < size; ++i )    
			C_paramsIdVec[i] = LocalGetNumObjectId(C_paramsIds[i]); 

		/// a function with a context
		exportFunc2Args< long, vector<long> >  ourFunc(C_ZeroCurveId, C_paramsIdVec, ARMLOCAL_Merton_Model_Create );

		/// call the general function
		fillXL_Result( LOCAL_MERTON_MODEL_CLASS, ourFunc, C_result, XL_result, PersistentInXL );
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_Merton_Model_Create_Common" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


///////////////////////////////////
/// version for persistentInXl
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_Merton_Model_Create(
	LPXLOPER XL_ZeroCurveId,
	LPXLOPER XL_ModelParamsId )
{
	ADD_LOG("Local_Merton_Model_Create");
	bool PersistentInXL = true;
	return Local_Merton_Model_Create_Common(
        XL_ZeroCurveId,
        XL_ModelParamsId,
        PersistentInXL );
}

///////////////////////////////////
/// version for non persistentInXl
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_PXL_Merton_Model_Create(
	LPXLOPER XL_ZeroCurveId,
	LPXLOPER XL_ModelParamsId )
{
	ADD_LOG("Local_PXL_Merton_Model_Create");
	bool PersistentInXL = false;
	return Local_Merton_Model_Create_Common(
        XL_ZeroCurveId,
        XL_ModelParamsId,
        PersistentInXL );
}


///----------------------------------------------
///----------------------------------------------
///             Normal Model
/// Inputs :
///     Zc curve Id
///     Vector of model params Id
///----------------------------------------------
///----------------------------------------------

LPXLOPER Local_Normal_Model_Create_Common(
	LPXLOPER XL_ZeroCurveId,
	LPXLOPER XL_ParamsIdVec,
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

		long C_ZeroCurveId;
		XL_GETOBJID( XL_ZeroCurveId, C_ZeroCurveId,	" ARM_ERR: Zero Curve: Object expected",		C_result);

		vector<CCString> C_paramsIds;
		XL_readStrVector (XL_ParamsIdVec,C_paramsIds," ARM_ERR: Model Parameters: array of object expected",DOUBLE_TYPE,C_result);
		size_t i, size = C_paramsIds.size();
		vector<long> C_paramsIdVec(size);
		for(i = 0; i < size; ++i )    
			C_paramsIdVec[i] = LocalGetNumObjectId(C_paramsIds[i]); 

		/// a function with a context
		exportFunc2Args< long, vector<long> >  ourFunc(C_ZeroCurveId, C_paramsIdVec, ARMLOCAL_Normal_Model_Create );

		/// call the general function
		fillXL_Result( LOCAL_NORMAL_MODEL_CLASS, ourFunc, C_result, XL_result, PersistentInXL );
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_Normal_Model_Create_Common" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


///////////////////////////////////
/// version for persistentInXl
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_Normal_Model_Create(
	LPXLOPER XL_ZeroCurveId,
	LPXLOPER XL_ModelParamsId )
{
	ADD_LOG("Local_Normal_Model_Create");
	bool PersistentInXL = true;
	return Local_Normal_Model_Create_Common(
        XL_ZeroCurveId,
        XL_ModelParamsId,
        PersistentInXL );
}

///////////////////////////////////
/// version for non persistentInXl
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_PXL_Normal_Model_Create(
	LPXLOPER XL_ZeroCurveId,
	LPXLOPER XL_ModelParamsId )
{
	ADD_LOG("Local_PXL_Normal_Model_Create");
	bool PersistentInXL = false;
	return Local_Normal_Model_Create_Common(
        XL_ZeroCurveId,
        XL_ModelParamsId,
        PersistentInXL );
}


///----------------------------------------------
///----------------------------------------------
///             SABR Model
/// Inputs :
///     Zc curve Id
///     Vector of model params Id
///----------------------------------------------
///----------------------------------------------

LPXLOPER Local_SABR_Model_Create_Common(
	LPXLOPER XL_ZeroCurveId,
	LPXLOPER XL_ParamsIdVec,
	LPXLOPER XL_ImpliedVolType,
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

		long C_ZeroCurveId;
		XL_GETOBJID( XL_ZeroCurveId, C_ZeroCurveId,	" ARM_ERR: Zero Curve: Object expected",		C_result);

		vector<CCString> C_paramsIds;
		XL_readStrVector (XL_ParamsIdVec,C_paramsIds," ARM_ERR: Model Parameters: array of object expected",DOUBLE_TYPE,C_result);
		size_t i, size = C_paramsIds.size();
		vector<long> C_paramsIdVec(size);
		for(i = 0; i < size; ++i )    
			C_paramsIdVec[i] = LocalGetNumObjectId(C_paramsIds[i]); 

		CCString C_impliedVoltr;
		CCString defaultImpliedVol = "DIRECTEXACT";
		XL_readStrCellWD(XL_ImpliedVolType,C_impliedVoltr,defaultImpliedVol," ARM_ERR: SABR Implied Vol: String expected",C_result);
		long impliedVolType = ARM_ConvGP_CFSABR_ImplicitVol_Formula_Extended_Flag(C_impliedVoltr,C_result);


		/// a function with a context
		exportFunc3Args< long, vector<long>, long >  ourFunc(C_ZeroCurveId, C_paramsIdVec, impliedVolType, ARMLOCAL_SABR_Model_Create );

		/// call the general function
		fillXL_Result( LOCAL_SABR_MODEL_CLASS, ourFunc, C_result, XL_result, PersistentInXL );
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_SABR_Model_Create_Common" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


///////////////////////////////////
/// version for persistentInXl
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_SABR_Model_Create(
	LPXLOPER XL_ZeroCurveId,
	LPXLOPER XL_ModelParamsId,
	LPXLOPER XL_ImpliedVolType)
{
	ADD_LOG("Local_SABR_Model_Create");
	bool PersistentInXL = true;
	return Local_SABR_Model_Create_Common(
        XL_ZeroCurveId,
        XL_ModelParamsId,
		XL_ImpliedVolType,
        PersistentInXL );
}

///////////////////////////////////
/// version for non persistentInXl
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_PXL_SABR_Model_Create(
	LPXLOPER XL_ZeroCurveId,
	LPXLOPER XL_ModelParamsId,
	LPXLOPER XL_ImpliedVolType)
{
	ADD_LOG("Local_PXL_SABR_Model_Create");
	bool PersistentInXL = false;
	return Local_SABR_Model_Create_Common(
        XL_ZeroCurveId,
        XL_ModelParamsId,
		XL_ImpliedVolType,
        PersistentInXL );
}

///----------------------------------------------
///----------------------------------------------
///             SLN Model
/// Inputs :
///     Zc curve Id
///     Vector of model params Id
///----------------------------------------------
///----------------------------------------------

LPXLOPER Local_SLN_Model_Create_Common(
	LPXLOPER XL_ZeroCurveId,
	LPXLOPER XL_ParamsIdVec,
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

		long C_ZeroCurveId;
		XL_GETOBJID( XL_ZeroCurveId, C_ZeroCurveId,	" ARM_ERR: Zero Curve: Object expected",		C_result);

		vector<CCString> C_paramsIds;
		XL_readStrVector (XL_ParamsIdVec,C_paramsIds," ARM_ERR: Model Parameters: array of object expected",DOUBLE_TYPE,C_result);
		size_t i, size = C_paramsIds.size();
		vector<long> C_paramsIdVec(size);
		for(i = 0; i < size; ++i )    
			C_paramsIdVec[i] = LocalGetNumObjectId(C_paramsIds[i]); 

		/// a function with a context
		exportFunc2Args< long, vector<long> >  ourFunc(C_ZeroCurveId, C_paramsIdVec, ARMLOCAL_SLN_Model_Create );

		/// call the general function
		fillXL_Result( LOCAL_SLN_MODEL_CLASS, ourFunc, C_result, XL_result, PersistentInXL );
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_SLN_Model_Create_Common" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


///////////////////////////////////
/// version for persistentInXl
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_SLN_Model_Create(
	LPXLOPER XL_ZeroCurveId,
	LPXLOPER XL_ModelParamsId )
{
	ADD_LOG("Local_SLN_Model_Create");
	bool PersistentInXL = true;
	return Local_SLN_Model_Create_Common(
        XL_ZeroCurveId,
        XL_ModelParamsId,
        PersistentInXL );
}

///////////////////////////////////
/// version for non persistentInXl
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_PXL_SLN_Model_Create(
	LPXLOPER XL_ZeroCurveId,
	LPXLOPER XL_ModelParamsId )
{
	ADD_LOG("Local_PXL_SLN_Model_Create");
	bool PersistentInXL = false;
	return Local_SLN_Model_Create_Common(
        XL_ZeroCurveId,
        XL_ModelParamsId,
        PersistentInXL );
}

///----------------------------------------------
///----------------------------------------------
///             Get Varaince squeeze info from a local model
/// Inputs :
///     model Id
///----------------------------------------------
///----------------------------------------------
LPXLOPER Local_GetVarianceSqueeze_Common(
	LPXLOPER XL_modelId,
	LPXLOPER XL_detail,
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

		long C_modelId;
		XL_GETOBJID( XL_modelId, C_modelId,	" ARM_ERR: Local model id : Object expected",		C_result);


		CCString C_detailStr;
		XL_readStrCellWD( XL_detail,C_detailStr,"N"," ARM_ERR: detail flag : boolean expected",C_result);
		bool detailFlag = false;
		C_detailStr.toUpper();
		if (C_detailStr == "Y" || C_detailStr == "YES")
			detailFlag = false;
		else if (C_detailStr == "N" || C_detailStr == "NO")
			detailFlag = false;
		else
			throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"String \"Y\" or \"N\" Expected for variance squeeze flag");


		/// a function with a context
		vector<double> C_DataResult;
		long C_rows,C_cols;
		bool  varSqueStatus;
		
		
		long retCode = ARMLOCAL_GetVarianceSqueeze(C_modelId,C_rows,C_cols,C_DataResult,varSqueStatus,C_result );

		/// feed the LPXLOPER object result 
		if (retCode == ARM_OK)
		{
			if (!detailFlag){
				FreeCurCellErr ();
				XL_result.xltype = xltypeStr;
				XL_result.val.str = XL_StrC2StrPascal (C_result.getString() );
				XL_result.xltype |= xlbitDLLFree;
			}

			else{
				/// add these additional lines 
				/// to display blank lines
				const int additionalLinesNb = 100;
				bool fillWithBlank = true;
				FreeCurCellContent ();
				XL_writeNumMatrixSizeWithOptions( XL_result, C_DataResult, C_rows, C_cols, " ARM_ERR: Could not get result data for variance squeeze info ", C_result, additionalLinesNb, fillWithBlank );
			}
		}
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_GetVarianceSqueeze_Common" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}

///////////////////////////////////
/// version for  persistentInXl
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_GetVarianceSqueeze(
	LPXLOPER XL_ModelId,
	LPXLOPER XL_detail)
{
	ADD_LOG("Local_GetVarianceSqueeze");
	bool PersistentInXL = true;
	return Local_GetVarianceSqueeze_Common(
        XL_ModelId,
		XL_detail,
        PersistentInXL );
}

///----------------------------------------------
///----------------------------------------------
///             HW1F Model Param Create
/// Inputs :
///     Vector of model params Id
///----------------------------------------------
///----------------------------------------------

LPXLOPER Local_HW1FModelParam_Create_Common(
	LPXLOPER XL_ParamsIdVec,
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

		vector<CCString> C_paramsIds;
		XL_readStrVector (XL_ParamsIdVec,C_paramsIds," ARM_ERR: Model Parameters: array of object expected",DOUBLE_TYPE,C_result);
		size_t i, size = C_paramsIds.size();
		vector<long> C_paramsIdVec(size);
		for(i = 0; i < size; ++i )    
			C_paramsIdVec[i] = LocalGetNumObjectId(C_paramsIds[i]); 

		/// a function with a context
		exportFunc1Arg< vector<long> >  ourFunc( C_paramsIdVec, ARMLOCAL_HW1FModelParam_Create );

		/// call the general function
		fillXL_Result( LOCAL_HW1F_MODELPARAM_CLASS, ourFunc, C_result, XL_result, PersistentInXL );
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_HW1FModelParam_Create_Common" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


///////////////////////////////////
/// version for non persistentInXl
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_HW1FModelParam_Create(
	LPXLOPER XL_ModelParamsId )
{
	ADD_LOG("Local_HW1FModelParam_Create");
	bool PersistentInXL = true;
	return Local_HW1FModelParam_Create_Common(
        XL_ModelParamsId,
        PersistentInXL );
}


///////////////////////////////////
/// version for non persistentInXl
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_PXL_HW1FModelParam_Create(
	LPXLOPER XL_ModelParamsId )
{
	ADD_LOG("Local_PXL_HW1FModelParam_Create");
	bool PersistentInXL = false;
	return Local_HW1FModelParam_Create_Common(
        XL_ModelParamsId,
        PersistentInXL );
}



///----------------------------------------------
///----------------------------------------------
///             QNF Model Param Create
/// Inputs :
///     QParam Id
///     Vector of model params Id
///     CorrelMat Id
///----------------------------------------------
///----------------------------------------------

LPXLOPER Local_QNFModelParam_Create_Common(
	LPXLOPER XL_QParamId,
	LPXLOPER XL_ParamsIdVec,
	LPXLOPER XL_CorrelMatId,
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

		CCString C_QParamIdStr;
		XL_readStrCell(XL_QParamId,C_QParamIdStr," ARM_ERR: Q Param: object expected",C_result);
		long C_QParamId=LocalGetNumObjectId(C_QParamIdStr);

		vector<CCString> C_paramsIds;
		XL_readStrVector (XL_ParamsIdVec,C_paramsIds," ARM_ERR: Model Parameters: array of object expected",DOUBLE_TYPE,C_result);
		size_t i, size = C_paramsIds.size();
		vector<long> C_paramsIdVec(size);
		for(i = 0; i < size; ++i )    
			C_paramsIdVec[i] = LocalGetNumObjectId(C_paramsIds[i]); 

		CCString C_CorrelMatIdStr;
		XL_readStrCell (XL_CorrelMatId,C_CorrelMatIdStr," ARM_ERR: Correlation Matrix: object expected", C_result);
		long C_CorrelMatId =LocalGetNumObjectId(C_CorrelMatIdStr);

		/// a function with a context
		exportFunc3Args< long, vector<long>, long >  ourFunc(C_QParamId, C_paramsIdVec, C_CorrelMatId, ARMLOCAL_QNFModelParam_Create );

		/// call the general function
		fillXL_Result( LOCAL_QNF_MODELPARAM_CLASS, ourFunc, C_result, XL_result, PersistentInXL );
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_QNFModelParam_Create_Common" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


///////////////////////////////////
/// version for non persistentInXl
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_QNFModelParam_Create(
	LPXLOPER XL_QParamId,
	LPXLOPER XL_ParamsIdVec,
	LPXLOPER XL_CorrelMatId )
{
	ADD_LOG("Local_QNFModelParam_Create");
	bool PersistentInXL = true;
	return Local_QNFModelParam_Create_Common(
		XL_QParamId,
        XL_ParamsIdVec,
		XL_CorrelMatId,
        PersistentInXL );
}


///////////////////////////////////
/// version for non persistentInXl
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_PXL_QNFModelParam_Create(
	LPXLOPER XL_QParamId,
	LPXLOPER XL_ParamsIdVec,
	LPXLOPER XL_CorrelMatId )
{
	ADD_LOG("Local_PXL_QNFModelParam_Create");
	bool PersistentInXL = false;
	return Local_QNFModelParam_Create_Common(
		XL_QParamId,
        XL_ParamsIdVec,
		XL_CorrelMatId,
        PersistentInXL );
}


///----------------------------------------------
///----------------------------------------------
///             QNF Model Create
/// Inputs :
///     Zc curve Id
///     QNF Model param Id
///----------------------------------------------
///----------------------------------------------

LPXLOPER Local_QNFModel_Create_Common(
	LPXLOPER XL_ZeroCurveId,
	LPXLOPER XL_QNFModelParamsId,
	LPXLOPER XL_DegenerateInHW,
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

		long C_ZeroCurveId;
		XL_GETOBJID( XL_ZeroCurveId, C_ZeroCurveId,	" ARM_ERR: Zero Curve: Object expected",		C_result);

		long C_QNFModelParamsId;
		XL_GETOBJID( XL_QNFModelParamsId, C_QNFModelParamsId,	" ARM_ERR: QNF Model Param: Object expected",		C_result);

		double DegenerateInHWDble;
		double DegenerateInHWDefault = false;
		XL_DegenerateInHW;
		XL_readNumCellWD( XL_DegenerateInHW, DegenerateInHWDble, DegenerateInHWDefault, " ARM_ERR: degenerated in HW: boolean compatible expected",C_result);
		bool C_DegenerateInHW= DegenerateInHWDble != 0;

		/// a function with a context
		exportFunc3Args< long, long, bool >  ourFunc(C_ZeroCurveId, C_QNFModelParamsId, C_DegenerateInHW, ARMLOCAL_QNFModel_Create );

		/// call the general function
		fillXL_Result( LOCAL_QNF_MODEL_CLASS, ourFunc, C_result, XL_result, PersistentInXL );
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_QNFModel_Create_Common" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


///////////////////////////////////
/// version for non persistentInXl
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_QNFModel_Create(
	LPXLOPER XL_ZeroCurveId,
	LPXLOPER XL_ModelParamsId,
	LPXLOPER XL_DegenerateInHW )
{
	ADD_LOG("Local_QNFModel_Create");
	bool PersistentInXL = true;
	return Local_QNFModel_Create_Common(
		XL_ZeroCurveId,
        XL_ModelParamsId,
		XL_DegenerateInHW,
        PersistentInXL );
}


///////////////////////////////////
/// version for non persistentInXl
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_PXL_QNFModel_Create(
	LPXLOPER XL_ZeroCurveId,
	LPXLOPER XL_ModelParamsId,
	LPXLOPER XL_DegenerateInHW )
{
	ADD_LOG("Local_PXL_QNFModel_Create");
	bool PersistentInXL = false;
	return Local_QNFModel_Create_Common(
		XL_ZeroCurveId,
        XL_ModelParamsId,
		XL_DegenerateInHW,
        PersistentInXL );
}



///----------------------------------------------
///----------------------------------------------
///             FX Model Create
/// Inputs :
///     Zc curve Id
///     Q1F Model param Id
///		Spot
///		Foreign Curve Id
///----------------------------------------------
///----------------------------------------------

LPXLOPER Local_Model_FX_Common(
	LPXLOPER XL_ZeroCurveId,
	LPXLOPER XL_ParamsIdVec,
	LPXLOPER XL_Spot,
	LPXLOPER XL_ForCurveId,
	LPXLOPER XL_MCScheme,
	long WhichModel,
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

		long C_ZeroCurveId;
		XL_GETOBJID( XL_ZeroCurveId, C_ZeroCurveId,	" ARM_ERR: Zero Curve: Object expected",		C_result);

		vector<CCString> C_paramsIds;
		XL_readStrVector (XL_ParamsIdVec,C_paramsIds," ARM_ERR: Model Parameters: array of object expected",DOUBLE_TYPE,C_result);
		size_t i, size = C_paramsIds.size();
		vector<long> C_paramsIdVec(size);
		for(i = 0; i < size; ++i )    
			C_paramsIdVec[i] = LocalGetNumObjectId(C_paramsIds[i]); 

		double C_Spot;
		XL_readNumCell( XL_Spot, C_Spot, " ARM_ERR: spot; double expected",C_result);

		long C_ForCurveId;
		XL_GETOBJID( XL_ForCurveId, C_ForCurveId, " ARM_ERR: foreign Curve Id: Object compatible expected",	C_result);

		CCString C_MCSchemeStr, C_MCSchemeDef="ANDREASEN";
		string C_MCScheme;
		XL_readStrCellWD(XL_MCScheme, C_MCSchemeStr, C_MCSchemeDef, "ARM_ERR: MC Scheme string expected.",C_result);

		C_MCScheme = CCSTringToSTLString(C_MCSchemeStr);

		exportFunc6Args< long, vector<long>, double, long, string, long >  ourFunc(C_ZeroCurveId, C_paramsIdVec, C_Spot, C_ForCurveId, C_MCScheme, WhichModel, ARMLOCAL_FXModel_Create );
				
		/// call the general function
		fillXL_Result_withName( ourFunc, C_result, XL_result, PersistentInXL );
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_Model_FX_Common" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


///////////////////////////////////
/// version for non persistentInXl
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_Q1FModel_FX_Create(
	LPXLOPER XL_ZeroCurveId,
	LPXLOPER XL_ParamsIdVec,
	LPXLOPER XL_Spot,
	LPXLOPER XL_ForCurveId )
{
	ADD_LOG("Local_Q1FModel_FX_Create");
	bool PersistentInXL = true;

	XLOPER XL_MCScheme;
	XL_MCScheme.xltype = xltypeMissing;

	return Local_Model_FX_Common(
		XL_ZeroCurveId,
        XL_ParamsIdVec,
		XL_Spot,
		XL_ForCurveId,
		&XL_MCScheme,
		1,			/// Q1F Model
        PersistentInXL );
}


///////////////////////////////////
/// version for non persistentInXl
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_PXL_Q1FModel_FX_Create(
	LPXLOPER XL_ZeroCurveId,
	LPXLOPER XL_ParamsIdVec,
	LPXLOPER XL_Spot,
	LPXLOPER XL_ForCurveId )
{
	ADD_LOG("Local_PXL_Q1FModel_FX_Create");
	bool PersistentInXL = false;

	XLOPER XL_MCScheme;
	XL_MCScheme.xltype = xltypeMissing;

	return Local_Model_FX_Common(
		XL_ZeroCurveId,
        XL_ParamsIdVec,
		XL_Spot,
		XL_ForCurveId,
		&XL_MCScheme,
		1,			/// Q1F Model
        PersistentInXL );
}

///////////////////////////////////
/// version for non persistentInXl
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_BSModel_FX_Create(
	LPXLOPER XL_ZeroCurveId,
	LPXLOPER XL_ParamsIdVec,
	LPXLOPER XL_Spot,
	LPXLOPER XL_ForCurveId )
{
	ADD_LOG("Local_BSModel_FX_Create");
	bool PersistentInXL = true;

	XLOPER XL_MCScheme;
	XL_MCScheme.xltype = xltypeMissing;

	return Local_Model_FX_Common(
		XL_ZeroCurveId,
        XL_ParamsIdVec,
		XL_Spot,
		XL_ForCurveId,
		&XL_MCScheme,
		0,			/// BS Model
        PersistentInXL );
}


///////////////////////////////////
/// version for non persistentInXl
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_PXL_BSModel_FX_Create(
	LPXLOPER XL_ZeroCurveId,
	LPXLOPER XL_ParamsIdVec,
	LPXLOPER XL_Spot,
	LPXLOPER XL_ForCurveId )
{
	ADD_LOG("Local_PXL_BSModel_FX_Create");
	bool PersistentInXL = false;

	XLOPER XL_MCScheme;
	XL_MCScheme.xltype = xltypeMissing;

	return Local_Model_FX_Common(
		XL_ZeroCurveId,
        XL_ParamsIdVec,
		XL_Spot,
		XL_ForCurveId,
		&XL_MCScheme,
		0,			/// BS Model
        PersistentInXL );
}

///////////////////////////////////
/// version for non persistentInXl
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_CEVModel_FX_Create(
	LPXLOPER XL_ZeroCurveId,
	LPXLOPER XL_ParamsIdVec,
	LPXLOPER XL_Spot,
	LPXLOPER XL_ForCurveId )
{
	ADD_LOG("Local_CEVModel_FX_Create");
	bool PersistentInXL = true;

	XLOPER XL_MCScheme;
	XL_MCScheme.xltype = xltypeMissing;

	return Local_Model_FX_Common(
		XL_ZeroCurveId,
        XL_ParamsIdVec,
		XL_Spot,
		XL_ForCurveId,
		&XL_MCScheme,
		4,			/// CEV Model
        PersistentInXL );
}


///////////////////////////////////
/// version for non persistentInXl
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_PXL_CEVModel_FX_Create(
	LPXLOPER XL_ZeroCurveId,
	LPXLOPER XL_ParamsIdVec,
	LPXLOPER XL_Spot,
	LPXLOPER XL_ForCurveId )
{
	ADD_LOG("Local_PXL_CEVModel_FX_Create");
	bool PersistentInXL = false;

	XLOPER XL_MCScheme;
	XL_MCScheme.xltype = xltypeMissing;

	return Local_Model_FX_Common(
		XL_ZeroCurveId,
        XL_ParamsIdVec,
		XL_Spot,
		XL_ForCurveId,
		&XL_MCScheme,
		4,			/// CEV Model
        PersistentInXL );
}
__declspec(dllexport) LPXLOPER WINAPI Local_HestonModel_FX_Create(
	LPXLOPER XL_ZeroCurveId,
	LPXLOPER XL_ParamsIdVec,
	LPXLOPER XL_Spot,
	LPXLOPER XL_ForCurveId,
	LPXLOPER XL_MCScheme)
{
	ADD_LOG("Local_HestonModel_FX_Create");
	bool PersistentInXL = true;
	return Local_Model_FX_Common(
		XL_ZeroCurveId,
        XL_ParamsIdVec,
		XL_Spot,
		XL_ForCurveId,
		XL_MCScheme,
		2,			/// Heston Model
        PersistentInXL );
}


///////////////////////////////////
/// version for non persistentInXl
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_PXL_HestonModel_FX_Create(
	LPXLOPER XL_ZeroCurveId,
	LPXLOPER XL_ParamsIdVec,
	LPXLOPER XL_Spot,
	LPXLOPER XL_ForCurveId,
	LPXLOPER XL_MCScheme )
{
	ADD_LOG("Local_PXL_HestonModel_FX_Create");
	bool PersistentInXL = false;

	return Local_Model_FX_Common(
		XL_ZeroCurveId,
        XL_ParamsIdVec,
		XL_Spot,
		XL_ForCurveId,
		XL_MCScheme,
		2,			/// Heston Model
        PersistentInXL );
}

///////////////////////////////////
/// version for non persistentInXl
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_SABRModel_FX_Create(
	LPXLOPER XL_ZeroCurveId,
	LPXLOPER XL_ParamsIdVec,
	LPXLOPER XL_Spot,
	LPXLOPER XL_ForCurveId )
{
	ADD_LOG("Local_SABRModel_FX_Create");
	bool PersistentInXL = true;

	XLOPER XL_MCScheme;
	XL_MCScheme.xltype = xltypeMissing;

	return Local_Model_FX_Common(
		XL_ZeroCurveId,
        XL_ParamsIdVec,
		XL_Spot,
		XL_ForCurveId,
		&XL_MCScheme,
		3,			/// SABR Model
        PersistentInXL );
}


///////////////////////////////////
/// version for non persistentInXl
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_PXL_SABRModel_FX_Create(
	LPXLOPER XL_ZeroCurveId,
	LPXLOPER XL_ParamsIdVec,
	LPXLOPER XL_Spot,
	LPXLOPER XL_ForCurveId )
{
	ADD_LOG("Local_PXL_SABRModel_FX_Create");
	bool PersistentInXL = false;

	XLOPER XL_MCScheme;
	XL_MCScheme.xltype = xltypeMissing;

	return Local_Model_FX_Common(
		XL_ZeroCurveId,
        XL_ParamsIdVec,
		XL_Spot,
		XL_ForCurveId,
		&XL_MCScheme,
		3,			/// SABR Model
        PersistentInXL );
}


///////////////////////////////////
/// version for non persistentInXl
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_MixtureModel_FX_Create(
	LPXLOPER XL_ZeroCurveId,
	LPXLOPER XL_ParamsIdVec,
	LPXLOPER XL_Spot,
	LPXLOPER XL_ForCurveId )
{
	ADD_LOG("Local_MixtureModel_FX_Create");
	bool PersistentInXL = true;

	XLOPER XL_MCScheme;
	XL_MCScheme.xltype = xltypeMissing;

	return Local_Model_FX_Common(
		XL_ZeroCurveId,
        XL_ParamsIdVec,
		XL_Spot,
		XL_ForCurveId,
		&XL_MCScheme,
		5,			/// Mixture Model
        PersistentInXL );
}


///////////////////////////////////
/// version for non persistentInXl
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_PXL_MixtureModel_FX_Create(
	LPXLOPER XL_ZeroCurveId,
	LPXLOPER XL_ParamsIdVec,
	LPXLOPER XL_Spot,
	LPXLOPER XL_ForCurveId )
{
	ADD_LOG("Local_PXL_MixtureModel_FX_Create");
	bool PersistentInXL = false;

	XLOPER XL_MCScheme;
	XL_MCScheme.xltype = xltypeMissing;

	return Local_Model_FX_Common(
		XL_ZeroCurveId,
        XL_ParamsIdVec,
		XL_Spot,
		XL_ForCurveId,
		&XL_MCScheme,
		5,			/// Mixture Model
        PersistentInXL );
}



///----------------------------------------------
///----------------------------------------------
///             Eq Model Create
/// Inputs :
///     Zc curve Id
///     Q1F Model param Id
///		Spot
///----------------------------------------------
///----------------------------------------------

LPXLOPER Local_Model_Eq_Common(
	LPXLOPER XL_ZeroCurveId,
	LPXLOPER XL_ParamsIdVec,
	LPXLOPER XL_Spot,
	long WhichModel,
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

		long C_ZeroCurveId;
		XL_GETOBJID( XL_ZeroCurveId, C_ZeroCurveId,	" ARM_ERR: Zero Curve: Object expected",		C_result);

		vector<CCString> C_paramsIds;
		XL_readStrVector (XL_ParamsIdVec,C_paramsIds," ARM_ERR: Model Parameters: array of object expected",DOUBLE_TYPE,C_result);
		size_t i, size = C_paramsIds.size();
		vector<long> C_paramsIdVec(size);
		for(i = 0; i < size; ++i )    
			C_paramsIdVec[i] = LocalGetNumObjectId(C_paramsIds[i]); 

		double C_Spot;
		XL_readNumCell( XL_Spot, C_Spot, " ARM_ERR: spot; double expected",C_result);

		/// a function with a context
		exportFunc4Args< long, vector<long>, double, long >  ourFunc(C_ZeroCurveId, C_paramsIdVec, C_Spot, WhichModel, ARMLOCAL_EqModel_Create );

		/// call the general function
		fillXL_Result( LOCAL_EQ_MODEL_CLASS, ourFunc, C_result, XL_result, PersistentInXL );
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_Model_Eq_Common" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}

///////////////////////////////////
/// version for non persistentInXl
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_HestonModel_Eq_Create(
	LPXLOPER XL_ZeroCurveId,
	LPXLOPER XL_ParamsIdVec,
	LPXLOPER XL_Spot )
{
	ADD_LOG("Local_HestonModel_Eq_Create");
	bool PersistentInXL = true;
	return Local_Model_Eq_Common(
		XL_ZeroCurveId,
        XL_ParamsIdVec,
		XL_Spot,
		2,		/// Heston Model
        PersistentInXL );
}


///////////////////////////////////
/// version for non persistentInXl
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_PXL_HestonModel_Eq_Create(
	LPXLOPER XL_ZeroCurveId,
	LPXLOPER XL_ParamsIdVec,
	LPXLOPER XL_Spot )
{
	ADD_LOG("Local_PXL_HestonModel_Eq_Create");
	bool PersistentInXL = false;

	return Local_Model_Eq_Common(
		XL_ZeroCurveId,
        XL_ParamsIdVec,
		XL_Spot,
		2,		/// Heston Model
        PersistentInXL );
}

///////////////////////////////////
/// version for non persistentInXl
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_Q1FModel_Eq_Create(
	LPXLOPER XL_ZeroCurveId,
	LPXLOPER XL_ParamsIdVec,
	LPXLOPER XL_Spot )
{
	ADD_LOG("Local_Q1FModel_Eq_Create");
	bool PersistentInXL = true;
	return Local_Model_Eq_Common(
		XL_ZeroCurveId,
        XL_ParamsIdVec,
		XL_Spot,
		1,		/// Q1F Model
        PersistentInXL );
}


///////////////////////////////////
/// version for non persistentInXl
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_PXL_Q1FModel_Eq_Create(
	LPXLOPER XL_ZeroCurveId,
	LPXLOPER XL_ParamsIdVec,
	LPXLOPER XL_Spot )
{
	ADD_LOG("Local_PXL_Q1FModel_Eq_Create");
	bool PersistentInXL = false;

	return Local_Model_Eq_Common(
		XL_ZeroCurveId,
        XL_ParamsIdVec,
		XL_Spot,
		1,		/// Q1F Model
        PersistentInXL );
}

///////////////////////////////////
/// version for non persistentInXl
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_BSModel_Eq_Create(
	LPXLOPER XL_ZeroCurveId,
	LPXLOPER XL_ParamsIdVec,
	LPXLOPER XL_Spot )
{
	ADD_LOG("Local_BSModel_Eq_Create");
	bool PersistentInXL = true;
	return Local_Model_Eq_Common(
		XL_ZeroCurveId,
        XL_ParamsIdVec,
		XL_Spot,
		0,		/// BS Model
        PersistentInXL );
}


///////////////////////////////////
/// version for non persistentInXl
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_PXL_BSModel_Eq_Create(
	LPXLOPER XL_ZeroCurveId,
	LPXLOPER XL_ParamsIdVec,
	LPXLOPER XL_Spot )
{
	ADD_LOG("Local_PXL_BSModel_Eq_Create");
	bool PersistentInXL = false;

	return Local_Model_Eq_Common(
		XL_ZeroCurveId,
        XL_ParamsIdVec,
		XL_Spot,
		0,		/// BS Model
        PersistentInXL );
}



///----------------------------------------------
///----------------------------------------------
///             Model Name Map Create
/// Inputs :
///     vector of names
///     vector of models id
///		vector of vector of string
///----------------------------------------------
///----------------------------------------------

LPXLOPER Local_ModelNameMap_Common(
	LPXLOPER XL_Names,
	LPXLOPER XL_ModelsIdsVec,
	LPXLOPER XL_OtherModelNames,
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

		vector<CCString> C_NamesCStr;
		XL_readStrVector (XL_Names,C_NamesCStr," ARM_ERR: names: array of string expected", DOUBLE_TYPE,C_result);
		vector<string> C_names(C_NamesCStr.size() );
		size_t i;
		for( i=0; i<C_names.size(); ++i )
			C_names[i] = CCSTringToSTLString(C_NamesCStr[i]);

		vector<CCString> C_ModelsIdsVec;
		XL_readStrVector(XL_ModelsIdsVec,C_ModelsIdsVec," ARM_ERR: Model: array of object expected",DOUBLE_TYPE,C_result);
		size_t size = C_ModelsIdsVec.size();
		vector<long> C_ModelsIdsVecLong(size);
		for(i = 0; i < size; ++i )    
			C_ModelsIdsVecLong[i] = LocalGetNumObjectId(C_ModelsIdsVec[i]); 

		vector<CCString> C_OtherModelNames;
		vector<CCString> C_OtherModelNamesDefault(size,CCString(""));
		long nbrows,nbcolumns;
		XL_readStrVectorAndSizeWD( XL_OtherModelNames,nbrows,nbcolumns,C_OtherModelNames,C_OtherModelNamesDefault," ARM_ERR: other names: array of string expected",DOUBLE_TYPE,C_result);
		vector< vector <string> > C_otherModelNames( nbrows );
		size_t j,k;
		for( i=0, k=0; i<nbrows; ++i )
		{
			for( j=0; j<nbcolumns; ++j )
			{
				if( C_OtherModelNames[k].GetLen() )
					C_otherModelNames[i].push_back( CCSTringToSTLString( C_OtherModelNames[k] ) );
				k++;
			}
		}

		/// a function with a context
		exportFunc3Args< vector<string>, vector<long>, vector< vector< string > > >  ourFunc(C_names, C_ModelsIdsVecLong, C_otherModelNames, ARMLOCAL_ModelNameMap_Create );

		/// call the general function
		fillXL_Result( LOCAL_MODELNAMEMAP_CLASS, ourFunc, C_result, XL_result, PersistentInXL );
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_ModelNameMap_Common" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


///////////////////////////////////
/// version for non persistentInXl
///////////////////////////////////
__declspec(dllexport) LPXLOPER Local_ModelNameMap_Create(
	LPXLOPER XL_Names,
	LPXLOPER XL_ModelsIdsVec,
	LPXLOPER XL_OtherModelNames )
{
	bool PersistentInXL = true;
	return Local_ModelNameMap_Common(
		XL_Names,
        XL_ModelsIdsVec,
		XL_OtherModelNames,
        PersistentInXL );
}


///////////////////////////////////
/// version for non persistentInXl
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_PXL_ModelNameMap_Create(
	LPXLOPER XL_Names,
	LPXLOPER XL_ModelsIdsVec,
	LPXLOPER XL_OtherModelNames )
{
	ADD_LOG("Local_PXL_ModelNameMap_Create");
	bool PersistentInXL = false;
	return Local_ModelNameMap_Common(
		XL_Names,
        XL_ModelsIdsVec,
		XL_OtherModelNames,
        PersistentInXL );
}


///----------------------------------------------
///----------------------------------------------
///             Multi-Assets Model Create
/// Inputs :
///     Model Name Map
///     Correlation matrix
///----------------------------------------------
///----------------------------------------------

LPXLOPER Local_MultiAssetsModel_Common(
	LPXLOPER XL_ModelNameMapId,
	LPXLOPER XL_CorrelationMatrixId,
	LPXLOPER XL_MultiAssetsModelName,
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

		long C_ModelNameMapId;
		XL_GETOBJID( XL_ModelNameMapId, C_ModelNameMapId,	" ARM_ERR: Model Name Map: Object expected",		C_result);

		long C_CorrelationMatrixId;
		XL_GETOBJIDWD( XL_CorrelationMatrixId, C_CorrelationMatrixId,"NULL OBJECT"," ARM_ERR: Correlation Matrix Id: Object expected",		C_result);

		CCString MultiAssetsModelNameXL;
		XL_readStrCellWD(XL_MultiAssetsModelName,MultiAssetsModelNameXL,"Unknown"," ARM_ERR: MultiAsset name: string expected",C_result);
		MultiAssetsModelNameXL.toUpper();
        string MultiAssetsModelName= CCSTringToSTLString(MultiAssetsModelNameXL);

		/// a function with a context
		exportFunc3Args< long, long, string >  ourFunc(C_ModelNameMapId, C_CorrelationMatrixId, MultiAssetsModelName, ARMLOCAL_MultiAssetsModel_Create );

		/// call the general function
		fillXL_Result( LOCAL_MULTIASSETSMODEL_CLASS, ourFunc, C_result, XL_result, PersistentInXL );
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_MultiAssetsModel_Common" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


///////////////////////////////////
/// version for non persistentInXl
///////////////////////////////////
__declspec(dllexport) LPXLOPER Local_MultiAssetsModel_Create(
	LPXLOPER XL_ModelNameMapId,
	LPXLOPER XL_CorrelationMatrixId,
	LPXLOPER XL_MultiAssetsModelName)
{
	bool PersistentInXL = true;
	return Local_MultiAssetsModel_Common(
		XL_ModelNameMapId,
		XL_CorrelationMatrixId,
		XL_MultiAssetsModelName,
		PersistentInXL );
}


///////////////////////////////////
/// version for non persistentInXl
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_PXL_MultiAssetsModel_Create(
	LPXLOPER XL_ModelNameMapId,
	LPXLOPER XL_CorrelationMatrixId,
	LPXLOPER XL_MultiAssetsModelName)
{
	ADD_LOG("Local_PXL_MultiAssetsModel_Create");
	bool PersistentInXL = false;
	return Local_MultiAssetsModel_Common(
		XL_ModelNameMapId,
		XL_CorrelationMatrixId,
		XL_MultiAssetsModelName,
		PersistentInXL );
}

///----------------------------------------------
///----------------------------------------------
///             Get Model Param Id from Model
/// Inputs :
///     ModelId
///----------------------------------------------
///----------------------------------------------

LPXLOPER Local_PricingModel_GetModelParamIdCommon(
	LPXLOPER XL_ModelId,
	LPXLOPER XL_ModelParamType,
	LPXLOPER XL_FactorNb,
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

		CCString C_modelStrId;
		XL_readStrCell( XL_ModelId, C_modelStrId,	" ARM_ERR: Model Id: Object expected",		C_result);
		long C_modelId = LocalGetNumObjectId(C_modelStrId);

		CCString C_ModelParamTypeStr;
		XL_readStrCell(XL_ModelParamType,C_ModelParamTypeStr," ARM_ERR: Param Type: String expected",C_result);
		long C_ModelParamType;
		if( (C_ModelParamType = ARM_ConvGPModelParam( C_ModelParamTypeStr, C_result)) == ARM_DEFAULT_ERR )
		{
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}

		double C_FactorNb;
		double C_FactorNbDefault=0.;
		long C_FactorNbLong;
		XL_readNumCellWD( XL_FactorNb, C_FactorNb, C_FactorNbDefault, " ARM_ERR: Factor Nb: numeric expected",C_result);

		C_FactorNbLong = C_FactorNb;

		/// a function with a context
		exportFunc3Args< long, long, long >  ourFunc(C_modelId, C_ModelParamType, C_FactorNbLong, ARMLOCAL_PricingModel_GetModelParamId );

		/// call the general function
		fillXL_Result( LOCAL_MODELPARAM_CLASS, ourFunc, C_result, XL_result, PersistentInXL );
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_PricingModel_GetModelParamIdCommon" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


///////////////////////////////////
/// version for non persistentInXl
///////////////////////////////////
__declspec(dllexport) LPXLOPER Local_PricingModel_GetModelParamId(
	LPXLOPER XL_ModelId,
	LPXLOPER XL_ModelParamType,
	LPXLOPER XL_FactorNb )
{
	bool PersistentInXL = true;
	return Local_PricingModel_GetModelParamIdCommon(
		XL_ModelId,
		XL_ModelParamType,
		XL_FactorNb,
		PersistentInXL );
}


///////////////////////////////////
/// version for non persistentInXl
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_PXL_PricingModel_GetModelParamId(
	LPXLOPER XL_ModelId,
	LPXLOPER XL_ModelParamType,
	LPXLOPER XL_FactorNb )
{
	ADD_LOG("Local_PXL_PricingModel_GetModelParamId");
	bool PersistentInXL = false;
	return Local_PricingModel_GetModelParamIdCommon(
		XL_ModelId,
		XL_ModelParamType,
		XL_FactorNb,
		PersistentInXL );
}

///----------------------------------------------
///----------------------------------------------
///             Get ModelNameMap from Model
/// Inputs :
///     Pricing Model
///----------------------------------------------
///----------------------------------------------

LPXLOPER Local_PricingModel_GetModelMapCommon(
	LPXLOPER XL_ModelId,
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

		CCString C_modelStrId;
		XL_readStrCell( XL_ModelId, C_modelStrId,	" ARM_ERR: Model Id: Object expected",		C_result);
		long C_modelId = LocalGetNumObjectId(C_modelStrId);

		/// a function with a context
		exportFunc1Arg< long>  ourFunc(C_modelId,ARMLOCAL_PricingModel_GetModelMap );

		/// call the general function
		fillXL_Result( LOCAL_MODELNAMEMAP_CLASS, ourFunc, C_result, XL_result, PersistentInXL );
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_PricingModel_GetModelMapCommon" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


///////////////////////////////////
/// version for non persistentInXl
///////////////////////////////////
__declspec(dllexport) LPXLOPER Local_PricingModel_GetModelMap(
	LPXLOPER XL_ModelId)
{
	bool PersistentInXL = true;
	return Local_PricingModel_GetModelMapCommon(
		XL_ModelId,
		PersistentInXL );
}


///////////////////////////////////
/// version for non persistentInXl
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_PXL_PricingModel_GetModelMap(
	LPXLOPER XL_ModelId)
{
	ADD_LOG("Local_PXL_PricingModel_GetModelMap");
	bool PersistentInXL = false;
	return Local_PricingModel_GetModelMapCommon(
		XL_ModelId,
		PersistentInXL );
}

///////////////////////////////////
/// version that takes into account 
/// previous creation of object
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_PricingModel_SetModelMap(
	LPXLOPER XL_modelId,
	LPXLOPER XL_ModelMapId )
{
	ADD_LOG("Local_PricingModel_SetModelMap");
	bool PersistentInXL = true;
	return Local_SetSomethingToModelCommon( XL_modelId, XL_ModelMapId, 
		"ModelMap", ARMLOCAL_PricingModel_SetModelMap, PersistentInXL );
}


///////////////////////////////////
/// version for non persistentInXl
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_PXL_PricingModel_SetModelMap(
	LPXLOPER XL_modelId,
	LPXLOPER XL_ModelMapId )
{
	ADD_LOG("Local_PXL_PricingModel_SetModelMap");
	bool PersistentInXL = false;
	return Local_SetSomethingToModelCommon( XL_modelId, XL_ModelMapId, 
		"ModelMap", ARMLOCAL_PricingModel_SetModelMap, PersistentInXL );
}

///----------------------------------------------
///----------------------------------------------
///             Get Model from ModelMap
/// Inputs :
///		ModelMap
///     Pricing Model
///		Name
///----------------------------------------------
///----------------------------------------------

LPXLOPER Local_GetModelFromModelMapCommon(
	LPXLOPER XL_ModelMapId,
	LPXLOPER XL_name,
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

		CCString C_modelMapStrId;
		XL_readStrCell( XL_ModelMapId, C_modelMapStrId,	" ARM_ERR: ModelMapId: Object expected",		C_result);
		long C_modelMapId = LocalGetNumObjectId(C_modelMapStrId);

		CCString C_nameStrId;
		XL_readStrCell( XL_name, C_nameStrId,	" ARM_ERR: ModelMapId: Object expected",		C_result);
		string C_name = CCSTringToSTLString(C_nameStrId);

		/// a function with a context
		exportFunc2Args< long, string>  ourFunc(C_modelMapId,C_name, ARMLOCAL_GetModelFromModelMap );

		/// call the general function
		fillXL_Result_withName( ourFunc, C_result, XL_result, PersistentInXL );
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_GetModelFromModelMapCommon" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


///////////////////////////////////
/// version for non persistentInXl
///////////////////////////////////
__declspec(dllexport) LPXLOPER Local_GetModelFromModelMap(
	LPXLOPER XL_ModelMapId,
	LPXLOPER XL_name)
{
	bool PersistentInXL = true;
	return Local_GetModelFromModelMapCommon(
		XL_ModelMapId,
		XL_name,
		PersistentInXL );
}


///////////////////////////////////
/// version for non persistentInXl
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_PXL_GetModelFromModelMap(
	LPXLOPER XL_ModelMapId,
	LPXLOPER XL_name)
{
	ADD_LOG("Local_PXL_GetModelFromModelMap");
	bool PersistentInXL = false;
	return Local_GetModelFromModelMapCommon(
		XL_ModelMapId,
		XL_name,
		PersistentInXL );
}

///----------------------------------------------
///----------------------------------------------
///             Set Model To Model Map
/// Inputs :
///     Pricing Model
///----------------------------------------------
///----------------------------------------------

LPXLOPER Local_SetModelToModelMapCommon(
	LPXLOPER XL_ModelMapId,
	LPXLOPER XL_name,
	LPXLOPER XL_ModelId,
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

		CCString C_modelMapStrId;
		XL_readStrCell( XL_ModelMapId, C_modelMapStrId,	" ARM_ERR: ModelMapId: Object expected",		C_result);
		long C_modelMapId = LocalGetNumObjectId(C_modelMapStrId);

		CCString C_nameStrId;
		XL_readStrCell( XL_name, C_nameStrId,	" ARM_ERR: ModelMapId: Object expected",		C_result);
		string C_name = CCSTringToSTLString(C_nameStrId);


		CCString C_modelStrId;
		XL_readStrCell( XL_ModelId, C_modelStrId,	" ARM_ERR: Model Id: Object expected",		C_result);
		long C_modelId = LocalGetNumObjectId(C_modelStrId);

		/// a function with a context
		exportFunc3Args< long, string, long>  ourFunc(C_modelMapId, C_name, C_modelId, ARMLOCAL_SetModelToModelMap );

		/// call the general function
		fillXL_Result( LOCAL_MODELNAMEMAP_CLASS, ourFunc, C_result, XL_result, PersistentInXL );
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_SetModelToModelMapCommon" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


///////////////////////////////////
/// version for non persistentInXl
///////////////////////////////////
__declspec(dllexport) LPXLOPER Local_SetModelToModelMap(
	LPXLOPER XL_ModelMapId,
	LPXLOPER XL_name,
	LPXLOPER XL_ModelId)
{
	bool PersistentInXL = true;
	return Local_SetModelToModelMapCommon(
		XL_ModelMapId,
		XL_name,
		XL_ModelId,
		PersistentInXL );
}


///////////////////////////////////
/// version for non persistentInXl
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_PXL_SetModelToModelMap(
		LPXLOPER XL_ModelMapId,
		LPXLOPER XL_name,
		LPXLOPER XL_ModelId)
{
	ADD_LOG("Local_PXL_SetModelToModelMap");
	bool PersistentInXL = false;
	return Local_SetModelToModelMapCommon(
		XL_ModelMapId,
		XL_name,
		XL_ModelId,
		PersistentInXL );
}






///----------------------------------------------
///----------------------------------------------
///             Local_Create2IRFXModelCommon
///		function to create a two factors interest rates and fx model
/// Inputs :
///     names
///		modelsIds
///		correlation matrix
///----------------------------------------------
///----------------------------------------------

LPXLOPER Local_Create2IRFXModelCommon(
	LPXLOPER XL_ModelNames,
	LPXLOPER XL_ModelsIdsVec,
	LPXLOPER XL_CorrelationId,
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

		vector<CCString> modelNames;
		XL_readStrVector(XL_ModelNames,modelNames," ARM_ERR: Model names : array of string expected",DOUBLE_TYPE,C_result);
		vector< string > C_modelNames(modelNames.size());
		for(size_t i=0;i<modelNames.size();++i)
			C_modelNames[i]=CCSTringToSTLString(modelNames[i]);

		vector<CCString> C_ModelsIdsVec;
		XL_readStrVector(XL_ModelsIdsVec,C_ModelsIdsVec," ARM_ERR: Model: array of object expected",DOUBLE_TYPE,C_result);
		size_t size = C_ModelsIdsVec.size();
		vector<long> C_ModelsIdsVecLong(size);
		for(i = 0; i < size; ++i )    
			C_ModelsIdsVecLong[i] = LocalGetNumObjectId(C_ModelsIdsVec[i]); 

		CCString C_CorrelationStrId;
		XL_readStrCell( XL_CorrelationId, C_CorrelationStrId,	" ARM_ERR: Correlation Id: Object expected",		C_result);
		long C_CorrelationId = LocalGetNumObjectId(C_CorrelationStrId);


		/// a function with a context
		exportFunc3Args< vector< string >, vector< long >, long>  ourFunc(C_modelNames, C_ModelsIdsVecLong, C_CorrelationId, ARMLOCAL_Create2IRFXModel );

		/// call the general function
		fillXL_Result_withName( ourFunc, C_result, XL_result, PersistentInXL );
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_Create2IRFXModelCommon" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


///////////////////////////////////
/// version for non persistentInXl
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_Create2IRFXModel(
	LPXLOPER XL_ModelNames,
	LPXLOPER XL_ModelsIdsVec,
	LPXLOPER XL_CorrelationId )
{
	ADD_LOG("Local_Create2IRFXModel");
	bool PersistentInXL = true;
	return Local_Create2IRFXModelCommon(
		XL_ModelNames,
		XL_ModelsIdsVec,
		XL_CorrelationId,
		PersistentInXL);
}


///////////////////////////////////
/// version for non persistentInXl
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_PXL_Create2IRFXModel(
	LPXLOPER XL_ModelNames,
	LPXLOPER XL_ModelsIdsVec,
	LPXLOPER XL_CorrelationId )
{
	ADD_LOG("Local_PXL_Create2IRFXModel");
	bool PersistentInXL = false;
	return Local_Create2IRFXModelCommon(
		XL_ModelNames,
		XL_ModelsIdsVec,
		XL_CorrelationId,
		PersistentInXL);
}



///----------------------------------------------
///----------------------------------------------
///             Local_Create1IRFXModelCommon
///		function to create a two factors interest rates and fx model
/// Inputs :
///     names
///		modelsIds
///		correlation matrix
///----------------------------------------------
///----------------------------------------------

LPXLOPER Local_Create1IRFXModelCommon(
	LPXLOPER XL_ModelNames,
	LPXLOPER XL_ModelsIdsVec,
	LPXLOPER XL_CorrelationId,
	LPXLOPER XL_Model2IRFXId,
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

		vector<CCString> modelNames;
		XL_readStrVector(XL_ModelNames,modelNames," ARM_ERR: Model names : array of string expected",DOUBLE_TYPE,C_result);
		vector< string > C_modelNames(modelNames.size());
		for(size_t i=0;i<modelNames.size();++i)
			C_modelNames[i]=CCSTringToSTLString(modelNames[i]);

		vector<CCString> C_ModelsIdsVec;
		XL_readStrVector(XL_ModelsIdsVec,C_ModelsIdsVec," ARM_ERR: Model: array of object expected",DOUBLE_TYPE,C_result);
		size_t size = C_ModelsIdsVec.size();
		vector<long> C_ModelsIdsVecLong(size);
		for(i = 0; i < size; ++i )    
			C_ModelsIdsVecLong[i] = LocalGetNumObjectId(C_ModelsIdsVec[i]); 

		CCString C_CorrelationStrId;
		XL_readStrCell( XL_CorrelationId, C_CorrelationStrId,	" ARM_ERR: Correlation Id: Object expected",		C_result);
		long C_CorrelationId = LocalGetNumObjectId(C_CorrelationStrId);

		long C_Model2IRFXId;
		XL_GETOBJIDWD( XL_Model2IRFXId, C_Model2IRFXId,	"NULL OBJECT",	" ARM_ERR: Model 2IRFX: object expected",	C_result);

		/// a function with a context
		exportFunc4Args< vector< string >, vector< long >, long, long>  ourFunc(C_modelNames, C_ModelsIdsVecLong, C_CorrelationId, C_Model2IRFXId, ARMLOCAL_Create1IRFXModel );

		/// call the general function
		fillXL_Result_withName( ourFunc, C_result, XL_result, PersistentInXL );
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_Create2IRFXModelCommon" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


///////////////////////////////////
/// version for non persistentInXl
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_Create1IRFXModel(
	LPXLOPER XL_ModelNames,
	LPXLOPER XL_ModelsIdsVec,
	LPXLOPER XL_CorrelationId,
	LPXLOPER XL_Model2IRFXId)
{
	ADD_LOG("Local_Create1IRFXModel");
	bool PersistentInXL = true;
	return Local_Create1IRFXModelCommon(
		XL_ModelNames,
		XL_ModelsIdsVec,
		XL_CorrelationId,
		XL_Model2IRFXId,
		PersistentInXL);
}


///////////////////////////////////
/// version for non persistentInXl
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_PXL_Create1IRFXModel(
	LPXLOPER XL_ModelNames,
	LPXLOPER XL_ModelsIdsVec,
	LPXLOPER XL_CorrelationId,
	LPXLOPER XL_Model2IRFXId)
{
	ADD_LOG("Local_PXL_Create1IRFXModel");
	bool PersistentInXL = false;
	return Local_Create1IRFXModelCommon(
		XL_ModelNames,
		XL_ModelsIdsVec,
		XL_CorrelationId,
		XL_Model2IRFXId,
		PersistentInXL);
}

///----------------------------------------------
///----------------------------------------------
///             Local_Create2IRFXModelCommon
///		function to create a two factors interest rates and fx model
/// Inputs :
///     names
///		modelsIds
///		correlation matrix
///----------------------------------------------
///----------------------------------------------

LPXLOPER Local_HWHWQtoModel_CreateCommon(
	LPXLOPER XL_ModelNames,
	LPXLOPER XL_ModelsIdsVec,
	LPXLOPER XL_CorrelationId,
	LPXLOPER XL_FxFlag,
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

		vector<CCString> modelNames;
		XL_readStrVector(XL_ModelNames,modelNames," ARM_ERR: Model names : array of string expected",DOUBLE_TYPE,C_result);
		vector< string > C_modelNames(modelNames.size());
		for(size_t i=0;i<modelNames.size();++i)
			C_modelNames[i]=CCSTringToSTLString(modelNames[i]);

		vector<CCString> C_ModelsIdsVec;
		XL_readStrVector(XL_ModelsIdsVec,C_ModelsIdsVec," ARM_ERR: Model: array of object expected",DOUBLE_TYPE,C_result);
		size_t size = C_ModelsIdsVec.size();
		vector<long> C_ModelsIdsVecLong(size);
		for(i = 0; i < size; ++i )    
			C_ModelsIdsVecLong[i] = LocalGetNumObjectId(C_ModelsIdsVec[i]); 

		CCString C_CorrelationStrId;
		XL_readStrCell( XL_CorrelationId, C_CorrelationStrId,	" ARM_ERR: Correlation Id: Object expected",		C_result);
		long C_CorrelationId = LocalGetNumObjectId(C_CorrelationStrId);


		CCString C_FxFlag;
		XL_readStrCellWD( XL_FxFlag,C_FxFlag,"N"," ARM_ERR: Fx Flag : string expected",C_result);
		bool C_FxFlagBool;
		C_FxFlag.toUpper();
		if (C_FxFlag == "Y" )
			C_FxFlagBool = true;
		else if (C_FxFlag == "N")
			C_FxFlagBool = false;
		else
			throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"String \"Y\" or \"N\" Expected for FX flag");

		/// a function with a context
		exportFunc4Args< vector< string >, vector< long >, long, bool>  ourFunc(C_modelNames, C_ModelsIdsVecLong, C_CorrelationId, C_FxFlagBool, ARMLOCAL_HWHWQtoModel_Create );

		/// call the general function
		fillXL_Result(LOCAL_HWHWQTOMODEL_CLASS, ourFunc, C_result, XL_result, PersistentInXL );
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_HWHWQtoModel_CreateCommon" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


///////////////////////////////////
/// version for non persistentInXl
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_HWHWQtoModel_Create(
	LPXLOPER XL_ModelNames,
	LPXLOPER XL_ModelsIdsVec,
	LPXLOPER XL_CorrelationId,
	LPXLOPER XL_FxFlag)
{
	ADD_LOG("Local_HWHWQtoModel_Create");
	bool PersistentInXL = true;
	return Local_HWHWQtoModel_CreateCommon(
		XL_ModelNames,
		XL_ModelsIdsVec,
		XL_CorrelationId,
		XL_FxFlag,
		PersistentInXL);
}


///////////////////////////////////
/// version for non persistentInXl
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_PXL_HWHWQtoModel_Create(
	LPXLOPER XL_ModelNames,
	LPXLOPER XL_ModelsIdsVec,
	LPXLOPER XL_CorrelationId,
	LPXLOPER XL_FxFlag )
{
	ADD_LOG("Local_PXL_HWHWQtoModel_Create");
	bool PersistentInXL = false;
	return Local_HWHWQtoModel_CreateCommon(
		XL_ModelNames,
		XL_ModelsIdsVec,
		XL_CorrelationId,
		XL_FxFlag,
		PersistentInXL);
}

///----------------------------------------------
///----------------------------------------------
///     Local_HWHW2FQtoModel_CreateCommon
///		function to create a domestic 1f HW, a foreign 2f HW and fx model
/// Inputs :
///     names
///		modelsIds
///		correlation matrix
///----------------------------------------------
///----------------------------------------------

LPXLOPER Local_HWHW2FQtoModel_CreateCommon(
	LPXLOPER XL_ModelNames,
	LPXLOPER XL_ModelsIdsVec,
	LPXLOPER XL_CorrelationId,
	LPXLOPER XL_FxFlag,
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

		vector<CCString> modelNames;
		XL_readStrVector(XL_ModelNames,modelNames," ARM_ERR: Model names : array of string expected",DOUBLE_TYPE,C_result);
		vector< string > C_modelNames(modelNames.size());
		for(size_t i=0;i<modelNames.size();++i)
			C_modelNames[i]=CCSTringToSTLString(modelNames[i]);

		vector<CCString> C_ModelsIdsVec;
		XL_readStrVector(XL_ModelsIdsVec,C_ModelsIdsVec," ARM_ERR: Model: array of object expected",DOUBLE_TYPE,C_result);
		size_t size = C_ModelsIdsVec.size();
		vector<long> C_ModelsIdsVecLong(size);
		for(i = 0; i < size; ++i )    
			C_ModelsIdsVecLong[i] = LocalGetNumObjectId(C_ModelsIdsVec[i]); 

		CCString C_CorrelationStrId;
		XL_readStrCell( XL_CorrelationId, C_CorrelationStrId,	" ARM_ERR: Correlation Id: Object expected",		C_result);
		long C_CorrelationId = LocalGetNumObjectId(C_CorrelationStrId);


		CCString C_FxFlag;
		XL_readStrCellWD( XL_FxFlag,C_FxFlag,"N"," ARM_ERR: Fx Flag : string expected",C_result);
		bool C_FxFlagBool;
		C_FxFlag.toUpper();
		if (C_FxFlag == "Y" )
			C_FxFlagBool = true;
		else if (C_FxFlag == "N")
			C_FxFlagBool = false;
		else
			throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"String \"Y\" or \"N\" Expected for FX flag");

		/// a function with a context
		exportFunc4Args< vector< string >, vector< long >, long, bool>  ourFunc(C_modelNames, C_ModelsIdsVecLong, C_CorrelationId, C_FxFlagBool, ARMLOCAL_HWHW2FQtoModel_Create );

		/// call the general function
		fillXL_Result(LOCAL_HWHW2FQTOMODEL_CLASS, ourFunc, C_result, XL_result, PersistentInXL );
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_HWHW2FQtoModel_CreateCommon" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


///////////////////////////////////
/// version for non persistentInXl
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_HWHW2FQtoModel_Create(
	LPXLOPER XL_ModelNames,
	LPXLOPER XL_ModelsIdsVec,
	LPXLOPER XL_CorrelationId,
	LPXLOPER XL_FxFlag)
{
	ADD_LOG("Local_HWHW2FQtoModel_Create");
	bool PersistentInXL = true;
	return Local_HWHW2FQtoModel_CreateCommon(
		XL_ModelNames,
		XL_ModelsIdsVec,
		XL_CorrelationId,
		XL_FxFlag,
		PersistentInXL);
}


///////////////////////////////////
/// version for non persistentInXl
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_PXL_HWHW2FQtoModel_Create(
	LPXLOPER XL_ModelNames,
	LPXLOPER XL_ModelsIdsVec,
	LPXLOPER XL_CorrelationId,
	LPXLOPER XL_FxFlag )
{
	ADD_LOG("Local_PXL_HWHW2FQtoModel_Create");
	bool PersistentInXL = false;
	return Local_HWHW2FQtoModel_CreateCommon(
		XL_ModelNames,
		XL_ModelsIdsVec,
		XL_CorrelationId,
		XL_FxFlag,
		PersistentInXL);
}

///----------------------------------------------
///----------------------------------------------
///             Local_CreateFwdMarginModelCommon
///		function to create a basis curve margin model
/// Inputs :
///     curve Id
///----------------------------------------------
///----------------------------------------------

LPXLOPER Local_CreateFwdMarginModelCommon(
	LPXLOPER XL_basisZcCurveId,
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

		long C_basisZcCurveId;
		XL_GETOBJID( XL_basisZcCurveId, C_basisZcCurveId, " ARM_ERR: Basis swap Zc curve : Object expected",C_result );

		/// a function with a context
		exportFunc1Arg< long>  ourFunc(C_basisZcCurveId, ARMLOCAL_CreateFwdMarginModel );

		/// call the general function
		fillXL_Result_withName( ourFunc, C_result, XL_result, PersistentInXL );
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_CreateFwdMarginModelCommon" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


///////////////////////////////////
/// version for non persistentInXl
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_CreateFwdMarginModel(
	LPXLOPER XL_basisZcCurveId )
{
	ADD_LOG("Local_CreateFwdMarginModel");
	bool PersistentInXL = true;
	return Local_CreateFwdMarginModelCommon(
		XL_basisZcCurveId,
		PersistentInXL);
}



///////////////////////////////////
/// version for non persistentInXl
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_PXL_CreateFwdMarginModel(
	LPXLOPER XL_basisZcCurveId )
{
	ADD_LOG("Local_PXL_CreateFwdMarginModel");
	bool PersistentInXL = false;
	return Local_CreateFwdMarginModelCommon(
		XL_basisZcCurveId,
		PersistentInXL);
}



///----------------------------------------------
///----------------------------------------------
///             Local_SetRefModelNameToMultiAsset
///		function to set the reference model to a multi-asset
///			using its name
/// Inputs :
///			Multi-asset Model
///			Name
///----------------------------------------------
///----------------------------------------------

LPXLOPER Local_SetRefModelNameToMultiAssetCommon(
	LPXLOPER XL_Name,
	LPXLOPER XL_MultiAssetModel,
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

		CCString C_NameStr;
		XL_readStrCell(XL_Name,C_NameStr," ARM_ERR: Name of The Ref Model: String expected",C_result);
		string C_Name = CCSTringToSTLString( C_NameStr );

		long C_MultiAssetModelId;
		XL_GETOBJID( XL_MultiAssetModel, C_MultiAssetModelId, " ARM_ERR: Multi Asset Model: Object expected",C_result );

		/// a function with a context
		exportFunc2Args<string,long>  ourFunc(C_Name, C_MultiAssetModelId, ARMLOCAL_SetRefModelNameToMultiAsset );

		/// call the general function
		fillXL_Result_withName( ourFunc, C_result, XL_result, PersistentInXL );
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_SetRefModelNameToMultiAssetCommon" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


///////////////////////////////////
/// version for non persistentInXl
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_SetRefModelNameToMultiAsset(
	LPXLOPER XL_Name,
	LPXLOPER XL_MultiAssetModel )
{
	ADD_LOG("Local_SetRefModelNameToMultiAsset");
	bool PersistentInXL = true;
	return Local_SetRefModelNameToMultiAssetCommon(
		XL_Name,
		XL_MultiAssetModel,
		PersistentInXL);
}



///////////////////////////////////
/// version for non persistentInXl
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_PXL_SetRefModelNameToMultiAsset(
	LPXLOPER XL_Name,
	LPXLOPER XL_MultiAssetModel )
{
	ADD_LOG("Local_PXL_SetRefModelNameToMultiAsset");
	bool PersistentInXL = false;
	return Local_SetRefModelNameToMultiAssetCommon(
		XL_Name,
		XL_MultiAssetModel,
		PersistentInXL);
}


///----------------------------------------------
///----------------------------------------------
///             Local Normal Model
/// Inputs :
///     Zc curve Id
///     Vector of model params Id
///----------------------------------------------
///----------------------------------------------

LPXLOPER Local_LocalNormal_Model_Create_Common(
	LPXLOPER XL_ZeroCurveId,
	LPXLOPER XL_ParamsIdVec,
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

		long C_ZeroCurveId;
		XL_GETOBJID( XL_ZeroCurveId, C_ZeroCurveId,	" ARM_ERR: Zero Curve: Object expected",		C_result);

		VECTOR<CCString> C_paramsIdsDef(0);
		vector<CCString> C_paramsIds;
		XL_readStrVectorWD (XL_ParamsIdVec,C_paramsIds,C_paramsIdsDef," ARM_ERR: Model Parameters: array of object expected",DOUBLE_TYPE,C_result);
		size_t i, size = C_paramsIds.size();
		vector<long> C_paramsIdVec(size);
		for(i = 0; i < size; ++i )    
			C_paramsIdVec[i] = LocalGetNumObjectId(C_paramsIds[i]); 


		/// a function with a context
		exportFunc2Args< long, vector<long> >  ourFunc(C_ZeroCurveId, C_paramsIdVec, ARMLOCAL_LocalNormal_Model_Create );

		/// call the general function
		fillXL_Result_withName(ourFunc, C_result, XL_result, PersistentInXL );
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_LocalNormal_Model_Create_Common" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


///////////////////////////////////
/// version for persistentInXl
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_LocalNormal_Model_Create(
	LPXLOPER XL_ZeroCurveId,
	LPXLOPER XL_ModelParamsId)
{
	ADD_LOG("Local_LocalNormal_Model_Create");
	bool PersistentInXL = true;
	return Local_LocalNormal_Model_Create_Common(
        XL_ZeroCurveId,
        XL_ModelParamsId,
		PersistentInXL );
}

///////////////////////////////////
/// version for non persistentInXl
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_PXL_LocalNormal_Model_Create(
	LPXLOPER XL_ZeroCurveId,
	LPXLOPER XL_ModelParamsId)
{
	ADD_LOG("Local_PXL_LocalNormal_Model_Create");
	bool PersistentInXL = false;
	return Local_LocalNormal_Model_Create_Common(
        XL_ZeroCurveId,
        XL_ModelParamsId,
		PersistentInXL );
}

///----------------------------------------------
///----------------------------------------------
///             Local Shifted LogNormal Model
/// Inputs :
///     Zc curve Id
///     Vector of model params Id
///----------------------------------------------
///----------------------------------------------

LPXLOPER Local_LocalSLN_Model_Create_Common(
	LPXLOPER XL_ParamsIdVec,
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

		vector<CCString> C_paramsIds;
		XL_readStrVector (XL_ParamsIdVec,C_paramsIds," ARM_ERR: Model Parameters: array of object expected",DOUBLE_TYPE,C_result);
		size_t i, size = C_paramsIds.size();
		vector<long> C_paramsIdVec(size);
		for(i = 0; i < size; ++i )    
			C_paramsIdVec[i] = LocalGetNumObjectId(C_paramsIds[i]); 


		/// a function with a context
		exportFunc1Arg< vector<long> >  ourFunc(C_paramsIdVec, ARMLOCAL_LocalSLN_Model_Create );

		/// call the general function
		fillXL_Result_withName(ourFunc, C_result, XL_result, PersistentInXL );
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_LocalSLN_Model_Create_Common" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


///////////////////////////////////
/// version for persistentInXl
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_LocalSLN_Model_Create(
	LPXLOPER XL_ModelParamsId)
{
	ADD_LOG("Local_LocalSLN_Model_Create");
	bool PersistentInXL = true;
	return Local_LocalSLN_Model_Create_Common(
        XL_ModelParamsId,
		PersistentInXL );
}

///////////////////////////////////
/// version for non persistentInXl
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_PXL_LocalSLN_Model_Create(
	LPXLOPER XL_ModelParamsId)
{
	ADD_LOG("Local_PXL_LocalSLN_Model_Create");
	bool PersistentInXL = false;
	return Local_LocalSLN_Model_Create_Common(
        XL_ModelParamsId,
		PersistentInXL );
}

///----------------------------------------------
///----------------------------------------------
///             Local Model Calibration
/// Inputs :
///     MultiAssets model Id
///     Name of Local Model embedded in MultiAssets
///		Portfolio Id
///		Eval Dates
///----------------------------------------------
///----------------------------------------------

LPXLOPER Local_Local_Model_Calibrate_Common(
	LPXLOPER XL_MultiAssetsModelId,
	LPXLOPER XL_LocalModelName,
	LPXLOPER XL_PortfolioId,
	LPXLOPER XL_EvalDates,
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

		long C_MultiAssetsModelId;
		XL_GETOBJID( XL_MultiAssetsModelId, C_MultiAssetsModelId, " ARM_ERR: Multi Asset Model: Object expected",C_result );

		CCString C_LocalModelNameStr;
		XL_readStrCell(XL_LocalModelName,C_LocalModelNameStr," ARM_ERR: Local Model Name: String expected",C_result);
		string C_LocalModelName = CCSTringToSTLString(C_LocalModelNameStr);

		long C_PortfolioId;
		XL_GETOBJID( XL_PortfolioId, C_PortfolioId, " ARM_ERR: Portfolio: Object expected",C_result );

		VECTOR<double> C_EvalDates;
		XL_readNumVector(XL_EvalDates,C_EvalDates," ARM_ERR: Eval dates: array of dates expected",C_result);


		/// a function with a context
		exportFunc4Args < long, string, long, vector<double> >  ourFunc (	C_MultiAssetsModelId, 
																			C_LocalModelName,
																			C_PortfolioId, 
																			C_EvalDates,
																			ARMLOCAL_Local_Model_Calibrate );

		/// call the general function
		fillXL_Result_withName(ourFunc, C_result, XL_result, PersistentInXL );
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_Local_Model_Calibrate_Common" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


///////////////////////////////////
/// version for persistentInXl
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_Local_Model_Calibrate(
	LPXLOPER XL_MultiAssetsModelId,
	LPXLOPER XL_LocalModelName,
	LPXLOPER XL_PortfolioId,
	LPXLOPER XL_EvalDates)
{
	ADD_LOG("Local_Local_Model_Calibrate");
	bool PersistentInXL = true;
	return Local_Local_Model_Calibrate_Common(
        XL_MultiAssetsModelId,
        XL_LocalModelName,
		XL_PortfolioId,
		XL_EvalDates,
		PersistentInXL );
}

///////////////////////////////////
/// version for non persistentInXl
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_PXL_Local_Model_Calibrate(
	LPXLOPER XL_MultiAssetsModelId,
	LPXLOPER XL_LocalModelName,
	LPXLOPER XL_PortfolioId,
	LPXLOPER XL_EvalDates)
{
	ADD_LOG("Local_PXL_Local_Model_Calibrate");
	bool PersistentInXL = false;
	return Local_Local_Model_Calibrate_Common(
        XL_MultiAssetsModelId,
        XL_LocalModelName,
		XL_PortfolioId,
		XL_EvalDates,
		PersistentInXL );
}

LPXLOPER Local_LocalModel_CalibrateFunctional_Common(
	LPXLOPER XL_MultiAssetsModelId,
	LPXLOPER XL_LocalModelName,
	LPXLOPER XL_SecuritiesId,
	LPXLOPER XL_DensitiesId,
	LPXLOPER XL_Rescaling,
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

		long C_MultiAssetsModelId;
		XL_GETOBJID( XL_MultiAssetsModelId, C_MultiAssetsModelId, " ARM_ERR: Multi Asset Model: Object expected",C_result );

		CCString C_LocalModelNameStr;
		XL_readStrCell(XL_LocalModelName,C_LocalModelNameStr," ARM_ERR: Local Model Name: String expected",C_result);
		string C_LocalModelName = CCSTringToSTLString(C_LocalModelNameStr);

		VECTOR<CCString> C_DensitiesDefault(0);
		VECTOR<CCString> C_DensitiesId;
		XL_readStrVectorWD (XL_DensitiesId,C_DensitiesId,C_DensitiesDefault," ARM_ERR: Securities: array of object expected",DOUBLE_TYPE,C_result);
		VECTOR<CCString> C_SecuritiesId;
		XL_readStrVector (XL_SecuritiesId,C_SecuritiesId," ARM_ERR: Securities: array of object expected",DOUBLE_TYPE,C_result);

		size_t sizeSec	= C_SecuritiesId.size();
		size_t sizeDens	= C_DensitiesId.size();
		if (sizeDens>0 && sizeDens != sizeSec)
			throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"Incoherent size for securities and densities");

		VECTOR<long> C_DensitiesIdLong (sizeDens);
		VECTOR<long> C_SecuritiesIdLong (sizeSec);
		size_t i;
		for(i = 0; i < sizeDens; ++i )
			C_DensitiesIdLong[i] = LocalGetNumObjectId(C_DensitiesId[i]); 

		for(i = 0; i < sizeSec; ++i )
			C_SecuritiesIdLong[i] = LocalGetNumObjectId(C_SecuritiesId[i]); 

		CCString C_RescalingStr;
		XL_readStrCellWD( XL_Rescaling,C_RescalingStr,"N"," ARM_ERR: Is Direct: string expected",C_result);
		bool C_RescalingBool;
		C_RescalingStr.toUpper();
		if (C_RescalingStr == "Y" )
			C_RescalingBool = true;
		else if (C_RescalingStr == "N")
			C_RescalingBool = false;
		else
			throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"String \"Y\" or \"N\" Expected for Is Direct");


		/// a function with a context
		exportFunc5Args < long, string, VECTOR<long>, VECTOR<long>, bool>  ourFunc (	C_MultiAssetsModelId, 
																			C_LocalModelName,
																			C_SecuritiesIdLong, 
																			C_DensitiesIdLong,
																			C_RescalingBool,
																			ARMLOCAL_LocalModel_CalibrateFunctional );

		/// call the general function
		fillXL_Result_withName(ourFunc, C_result, XL_result, PersistentInXL );
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_Local_Model_Calibrate_Common" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


///////////////////////////////////
/// version for persistentInXl
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_LocalModel_CalibrateFunctional(
	LPXLOPER XL_MultiAssetsModelId,
	LPXLOPER XL_LocalModelName,
	LPXLOPER XL_SecIds,
	LPXLOPER XL_DensIds,
	LPXLOPER XL_Rescaling)
{
	ADD_LOG("Local_LocalModel_CalibrateFunctional");
	bool PersistentInXL = true;
	return Local_LocalModel_CalibrateFunctional_Common(
        XL_MultiAssetsModelId,
        XL_LocalModelName,
		XL_SecIds,
		XL_DensIds,
		XL_Rescaling,
		PersistentInXL );
}

///////////////////////////////////
/// version for non persistentInXl
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_PXL_LocalModel_CalibrateFunctional(
	LPXLOPER XL_MultiAssetsModelId,
	LPXLOPER XL_LocalModelName,
	LPXLOPER XL_SecIds,
	LPXLOPER XL_DensIds,
	LPXLOPER XL_Rescaling)
{
	ADD_LOG("Local_PXL_LocalModel_CalibrateFunctional");
	bool PersistentInXL = false;
	return Local_LocalModel_CalibrateFunctional_Common(
        XL_MultiAssetsModelId,
        XL_LocalModelName,
		XL_SecIds,
		XL_DensIds,
		XL_Rescaling,
		PersistentInXL );
}

///////////////////////////////////////
// Get Warning when variance squeeze for local model
///////////////////////////////////////

LPXLOPER Local_Local_Model_VarSqueeze_Common(
	LPXLOPER XL_MultiAssetsModelId,
	LPXLOPER XL_LocalModelName)
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

		long C_MultiAssetsModelId;
		XL_GETOBJID( XL_MultiAssetsModelId, C_MultiAssetsModelId, " ARM_ERR: Multi Asset Model: Object expected",C_result );

		CCString C_LocalModelNameStr;
		XL_readStrCell(XL_LocalModelName,C_LocalModelNameStr," ARM_ERR: Local Model Name: String expected",C_result);
		string C_LocalModelName = CCSTringToSTLString(C_LocalModelNameStr);


		exportFunc2Args<long, string> ourFunc( 
								C_MultiAssetsModelId,
								C_LocalModelName,
								ARMLOCAL_Local_Model_VarSqueeze);

		/// Simple updating of the calculator
		long retCode = ourFunc(C_result,ARM_NULL_OBJECT_ID);

		/// feed the LPXLOPER object result 
		if (retCode == ARM_OK)
		{
			FreeCurCellErr ();
			XL_result.xltype = xltypeStr;
			XL_result.val.str = XL_StrC2StrPascal (C_result.getString() );
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_CRASet_Common" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}
							 


///////////////////////////////////
/// version that takes into account 
/// previous creation of object
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_Local_Model_VarSqueeze(
	LPXLOPER XL_craId,
	LPXLOPER XL_getType)
{
	ADD_LOG("Local_Local_Model_VarSqueeze");
	bool PersistentInXL = true;
	return Local_Local_Model_VarSqueeze_Common(
	    XL_craId,
	    XL_getType);
}


///----------------------------------------------
///----------------------------------------------
///             Market IR Model
/// Inputs :
///     MktDataManager Id
///     Keys
///----------------------------------------------
///----------------------------------------------

LPXLOPER Local_MarketIRModel_Create_Common(
	LPXLOPER XL_MktDatas,
	LPXLOPER XL_VnsPricingMethod,
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

		VECTOR<CCString> mktDatas;
		XL_readStrVector(XL_MktDatas,mktDatas," ARM_ERR: Market datas: array of string expected",DOUBLE_TYPE,C_result);
        
		long mktDataManagerId = LocalGetNumObjectId (mktDatas[0]);
		vector< string > keys(mktDatas.size()-1);
        size_t i;
		for(i=1;i<mktDatas.size();++i)
			keys[i-1] = CCSTringToSTLString(mktDatas[i]);


		CCString C_VnsPricingMethodStr;
		XL_readStrCellWD(XL_VnsPricingMethod,C_VnsPricingMethodStr, "MONEYNESS"," ARM_ERR: VNS vol type: String expected",C_result);
		string C_VnsPricingMethod = CCSTringToSTLString(C_VnsPricingMethodStr);

		int vnsmethod = ARM_ArgConv_VnsPricingMethod.GetNumber(C_VnsPricingMethod);

		/// a function with a context
		exportFunc3Args< long, vector<string>, int >  ourFunc(mktDataManagerId, keys, vnsmethod, ARMLOCAL_MarketIRModel_Create );

		/// call the general function
		fillXL_Result_withName(ourFunc, C_result, XL_result, PersistentInXL );
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_MarketIRModel_Create_Common" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


///////////////////////////////////
/// version for persistentInXl
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_MarketIRModel_Create(
	LPXLOPER XL_MktDatas, LPXLOPER XL_VnsPricingMethod)
{
	ADD_LOG("Local_MarketIRModel_Create");
	bool PersistentInXL = true;
	return Local_MarketIRModel_Create_Common(
        XL_MktDatas,
		XL_VnsPricingMethod,
        PersistentInXL );
}

///////////////////////////////////
/// version for non persistentInXl
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_PXL_MarketIRModel_Create(
	LPXLOPER XL_MktDatas, LPXLOPER XL_VnsPricingMethod)
{
	ADD_LOG("Local_PXL_MarketIRModel_Create");
	bool PersistentInXL = false;
	return Local_MarketIRModel_Create_Common(
        XL_MktDatas,
		XL_VnsPricingMethod,
        PersistentInXL );
}


///----------------------------------------------
///----------------------------------------------
///             SmiledFRM Model
/// Inputs :
///     Zc curve Id
///     Hump Param Id
///     Correlation Param Id
///		FactorsNb
///----------------------------------------------
///----------------------------------------------


/////////////////////////////////////////////////////////////
/// central function that does the creation of the XL function
/////////////////////////////////////////////////////////////
LPXLOPER Local_SmiledFRMModelCommon(
	LPXLOPER XL_ZeroCurveId,
	LPXLOPER XL_CorrelParamId,
	LPXLOPER XL_HumpId,
	LPXLOPER XL_FactorsNb,
	LPXLOPER XL_TimeStepsNb,	
	LPXLOPER XL_GridSize,		
	LPXLOPER XL_StdDevNb,
	LPXLOPER XL_SkipPDE,
	LPXLOPER XL_CorrelType,
	LPXLOPER XL_AllowInterpol,
	LPXLOPER XL_SwaptionApprox,
	LPXLOPER XL_Recorrel,
	LPXLOPER XL_Rescalling,
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

		long C_ZeroCurveId;
		XL_GETOBJID( XL_ZeroCurveId, C_ZeroCurveId,	" ARM_ERR: Zero Curve: Object expected",C_result);
		
		long C_CorrelParamId;
		XL_GETOBJID( XL_CorrelParamId, C_CorrelParamId,	" ARM_ERR: CorrelParam: Object expected",C_result);
		
		long C_HumpId;
		XL_GETOBJID( XL_HumpId, C_HumpId,	" ARM_ERR: Hump: Object expected",C_result);
		
		double C_FactorsNb;
		XL_readNumCell( XL_FactorsNb, C_FactorsNb, " ARM_ERR: Nb of factors: numeric expected",C_result);
		int C_FactorsNbInt = C_FactorsNb;

		double C_TimeStepsNb;
		double C_TimeStepsNbDefault=1000;
		XL_readNumCellWD(XL_TimeStepsNb, C_TimeStepsNb, C_TimeStepsNbDefault, " ARM_ERR: dim : numeric expected", C_result);
		int C_TimeStepsNbInt = C_TimeStepsNb;

		double C_GridSize;
		double C_GridSizeDefault=801;
		XL_readNumCellWD(XL_GridSize, C_GridSize, C_GridSizeDefault, " ARM_ERR: GridSize : numeric expected", C_result);
		int C_GridSizeInt = C_GridSize;

		double C_StdDevNb;
		double C_StdDevNbDefault=6;
		XL_readNumCellWD(XL_StdDevNb, C_StdDevNb, C_StdDevNbDefault, " ARM_ERR: StdDevNb : numeric expected", C_result);
		
		double C_Recorrel;
		double C_RecorrelDefault=0;
		XL_readNumCellWD(XL_Recorrel, C_Recorrel, C_RecorrelDefault, " ARM_ERR: Recorrel : numeric expected", C_result);
		
		CCString C_SkipPDE;
		XL_readStrCellWD( XL_SkipPDE,C_SkipPDE,"Y"," ARM_ERR: Skip PDE: string expected",C_result);
		bool C_SkipPDEBool;
		C_SkipPDE.toUpper();
		if (C_SkipPDE == "Y" )
			C_SkipPDEBool = true;
		else if (C_SkipPDE == "N")
			C_SkipPDEBool = false;
		else
			throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"String \"Y\" or \"N\" Expected for Skip PDE");

		CCString C_CorrelTypeStr;
		XL_readStrCellWD(XL_CorrelType,C_CorrelTypeStr,"Beta"," ARM_ERR: Param Type: String expected",C_result);
		string C_CorrelType = CCSTringToSTLString(C_CorrelTypeStr);

		double C_AllowInterpol;
		double C_AllowInterpolDefault=0;
		XL_readNumCellWD(XL_AllowInterpol, C_AllowInterpol, C_AllowInterpolDefault, " ARM_ERR: AllowInterpol : numeric expected", C_result);
		bool C_AllowInterpolBool = (C_AllowInterpol==1);
		
		CCString C_SwaptionApproxStr;
		XL_readStrCellWD(XL_SwaptionApprox,C_SwaptionApproxStr,"Local+"," ARM_ERR: Param Type: String expected",C_result);
		string C_SwaptionApprox = CCSTringToSTLString(C_SwaptionApproxStr);

		CCString C_Rescalling;
		XL_readStrCellWD( XL_Rescalling,C_Rescalling,"N"," ARM_ERR: Rescalling: string expected",C_result);
		bool C_RescallingBool;
		C_Rescalling.toUpper();
		if (C_Rescalling == "Y" )
			C_RescallingBool = true;
		else if (C_Rescalling == "N")
			C_RescallingBool = false;
		else
			throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"String \"Y\" or \"N\" Expected for Rescalling");

		exportFunc13Args< long, long ,long ,int ,int ,int ,double ,bool, string, bool, string , double, bool>  ourFunc(
			C_ZeroCurveId,
			C_CorrelParamId,
			C_HumpId,
			C_FactorsNbInt,
			C_TimeStepsNbInt,
			C_GridSizeInt,
			C_StdDevNb,
			C_SkipPDEBool,
			C_CorrelType,
			C_AllowInterpolBool,
			C_SwaptionApprox,
			C_Recorrel,
			C_RescallingBool,
			ARMLOCAL_SmiledFRMModel_Create);

		/// call the general function
		fillXL_Result( LOCAL_SMILEDFRMMODEL_CLASS, ourFunc, C_result, XL_result, PersistentInXL );
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_SmiledFRMModelCommon" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}

///////////////////////////////////
/// version that takes into account 
/// previous creation of object
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_SmiledFRMModel_Create(
	LPXLOPER XL_ZeroCurveId,
	LPXLOPER XL_CorrelParamId,
	LPXLOPER XL_HumpId,
	LPXLOPER XL_FactorsNb,
	LPXLOPER XL_TimeStepsNb,	
	LPXLOPER XL_GridSize,		
	LPXLOPER XL_StdDevNb,
	LPXLOPER XL_SkipPDE,
	LPXLOPER XL_SwitchTheta,
	LPXLOPER XL_AllowInterpol,
	LPXLOPER XL_SwaptionApprox,
	LPXLOPER XL_Recorrel,
	LPXLOPER XL_Rescalling)
{
	ADD_LOG("Local_SmiledFRMModel_Create");
	bool PersistentInXL = true;
	return Local_SmiledFRMModelCommon(
        XL_ZeroCurveId,
		XL_CorrelParamId,
		XL_HumpId,
		XL_FactorsNb,
		XL_TimeStepsNb,	
		XL_GridSize,		
		XL_StdDevNb,
		XL_SkipPDE,
		XL_SwitchTheta,
		XL_AllowInterpol,
		XL_SwaptionApprox,
		XL_Recorrel,
		XL_Rescalling,
		PersistentInXL );
}

///////////////////////////////////
/// version for non persistentInXl
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_PXL_SmiledFRMModel_Create(
	LPXLOPER XL_ZeroCurveId,
	LPXLOPER XL_CorrelParamId,
	LPXLOPER XL_HumpId,
	LPXLOPER XL_FactorsNb,
	LPXLOPER XL_TimeStepsNb,	
	LPXLOPER XL_GridSize,		
	LPXLOPER XL_StdDevNb,
	LPXLOPER XL_SkipPDE,
	LPXLOPER XL_SwitchTheta,
	LPXLOPER XL_AllowInterpol,
	LPXLOPER XL_SwaptionApprox,
	LPXLOPER XL_Recorrel,
	LPXLOPER XL_Rescalling)
{
	ADD_LOG("Local_PXL_SmiledFRMModel_Create");
	bool PersistentInXL = false;
	return Local_SmiledFRMModelCommon(
        XL_ZeroCurveId,
		XL_CorrelParamId,
		XL_HumpId,
		XL_FactorsNb,
		XL_TimeStepsNb,	
		XL_GridSize,		
		XL_StdDevNb,
		XL_SkipPDE,
		XL_SwitchTheta,
		XL_AllowInterpol,
		XL_SwaptionApprox,
		XL_Recorrel,
		XL_Rescalling,
		PersistentInXL );
}

///----------------------------------------------
///----------------------------------------------
///             SmiledFRM Model
/// Inputs :
///     Zc curve Id
///     Hump Param Id
///     Correlation Param Id
///		FactorsNb
///----------------------------------------------
///----------------------------------------------


/////////////////////////////////////////////////////////////
/// central function that does the creation of the XL function
/////////////////////////////////////////////////////////////
LPXLOPER Local_SmiledMarketModelCommon(
	LPXLOPER XL_CalibPattern,
	LPXLOPER XL_ZeroCurveId,
	LPXLOPER XL_CorrelParamId,
	LPXLOPER XL_HumpId,
	LPXLOPER XL_FactorsNb,
	LPXLOPER XL_TimeStepsNb,	
	LPXLOPER XL_GridSize,		
	LPXLOPER XL_StdDevNb,
	LPXLOPER XL_SkipPDE,
	LPXLOPER XL_CorrelType,
	LPXLOPER XL_AllowInterpol,
	LPXLOPER XL_CalibProxy,
	LPXLOPER XL_Recorrel,
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

		CCString C_CalibPatternStr;
		XL_readStrCellWD(XL_CalibPattern,C_CalibPatternStr,"LIBOR"," ARM_ERR: Param Type: String expected",C_result);
		string C_CalibPattern = CCSTringToSTLString(C_CalibPatternStr);

		long C_ZeroCurveId;
		XL_GETOBJID( XL_ZeroCurveId, C_ZeroCurveId,	" ARM_ERR: Zero Curve: Object expected",C_result);
		
		long C_CorrelParamId;
		XL_GETOBJID( XL_CorrelParamId, C_CorrelParamId,	" ARM_ERR: CorrelParam: Object expected",C_result);
		
		long C_HumpId;
		XL_GETOBJID( XL_HumpId, C_HumpId,	" ARM_ERR: Hump: Object expected",C_result);
		
		double C_FactorsNb;
		XL_readNumCell( XL_FactorsNb, C_FactorsNb, " ARM_ERR: Nb of factors: numeric expected",C_result);
		int C_FactorsNbInt = C_FactorsNb;

		double C_TimeStepsNb;
		double C_TimeStepsNbDefault=1000;
		XL_readNumCellWD(XL_TimeStepsNb, C_TimeStepsNb, C_TimeStepsNbDefault, " ARM_ERR: dim : numeric expected", C_result);
		int C_TimeStepsNbInt = C_TimeStepsNb;

		double C_GridSize;
		double C_GridSizeDefault=801;
		XL_readNumCellWD(XL_GridSize, C_GridSize, C_GridSizeDefault, " ARM_ERR: GridSize : numeric expected", C_result);
		int C_GridSizeInt = C_GridSize;

		double C_StdDevNb;
		double C_StdDevNbDefault=6;
		XL_readNumCellWD(XL_StdDevNb, C_StdDevNb, C_StdDevNbDefault, " ARM_ERR: StdDevNb : numeric expected", C_result);
		
		double C_Recorrel;
		double C_RecorrelDefault=0;
		XL_readNumCellWD(XL_Recorrel, C_Recorrel, C_RecorrelDefault, " ARM_ERR: Recorrel : numeric expected", C_result);
		
		CCString C_SkipPDE;
		XL_readStrCellWD( XL_SkipPDE,C_SkipPDE,"Y"," ARM_ERR: Skip PDE: string expected",C_result);
		bool C_SkipPDEBool;
		C_SkipPDE.toUpper();
		if (C_SkipPDE == "Y" )
			C_SkipPDEBool = true;
		else if (C_SkipPDE == "N")
			C_SkipPDEBool = false;
		else
			throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"String \"Y\" or \"N\" Expected for Skip PDE");

		CCString C_AllowInterpol;
		XL_readStrCellWD( XL_AllowInterpol,C_AllowInterpol,"N"," ARM_ERR: Skip PDE: string expected",C_result);
		bool C_AllowInterpolBool;
		C_AllowInterpol.toUpper();
		if (C_AllowInterpol == "Y" )
			C_AllowInterpolBool = true;
		else if (C_AllowInterpol == "N")
			C_AllowInterpolBool = false;
		else
			throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"String \"Y\" or \"N\" Expected for Allow Interpol");

		CCString C_CalibProxyStr;
		XL_readStrCellWD(XL_CalibProxy,C_CalibProxyStr,"Local"," ARM_ERR: Param Type: String expected",C_result);
		string C_CalibProxy = CCSTringToSTLString(C_CalibProxyStr);

		CCString C_CorrelTypeStr;
		XL_readStrCellWD(XL_CorrelType,C_CorrelTypeStr,"Beta"," ARM_ERR: Param Type: String expected",C_result);
		string C_CorrelType = CCSTringToSTLString(C_CorrelTypeStr);

		exportFunc13Args<string, long, long ,long ,int ,int ,int ,double ,bool, string, bool, string, double >  ourFunc(
			C_CalibPattern,
			C_ZeroCurveId,
			C_CorrelParamId,
			C_HumpId,
			C_FactorsNbInt,
			C_TimeStepsNbInt,
			C_GridSizeInt,
			C_StdDevNb,
			C_SkipPDEBool,
			C_CorrelType,
			C_AllowInterpolBool,
			C_CalibProxy,
			C_Recorrel,
			ARMLOCAL_SmiledMarketModel_Create);

		/// call the general function
		fillXL_Result( LOCAL_SMILEDFRMMODEL_CLASS, ourFunc, C_result, XL_result, PersistentInXL );
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_SmiledMarketModelCommon" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}

///////////////////////////////////
/// version that takes into account 
/// previous creation of object
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_SmiledMarketModel_Create(
	LPXLOPER XL_Pattern,
	LPXLOPER XL_ZeroCurveId,
	LPXLOPER XL_CorrelParamId,
	LPXLOPER XL_HumpId,
	LPXLOPER XL_FactorsNb,
	LPXLOPER XL_TimeStepsNb,	
	LPXLOPER XL_GridSize,		
	LPXLOPER XL_StdDevNb,
	LPXLOPER XL_SkipPDE,
	LPXLOPER XL_SwitchTheta,
	LPXLOPER XL_AllowInterpol,
	LPXLOPER XL_CalibProxy,
	LPXLOPER XL_Recorrel)
{
	ADD_LOG("Local_SmiledMarketModel_Create");
	bool PersistentInXL = true;
	return Local_SmiledMarketModelCommon(
		XL_Pattern,
        XL_ZeroCurveId,
		XL_CorrelParamId,
		XL_HumpId,
		XL_FactorsNb,
		XL_TimeStepsNb,	
		XL_GridSize,		
		XL_StdDevNb,
		XL_SkipPDE,
		XL_SwitchTheta,
		XL_AllowInterpol,
		XL_CalibProxy,
		XL_Recorrel,
		PersistentInXL );
}

///////////////////////////////////
/// version for non persistentInXl
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_PXL_SmiledMarketModel_Create(
	LPXLOPER XL_Pattern,
	LPXLOPER XL_ZeroCurveId,
	LPXLOPER XL_CorrelParamId,
	LPXLOPER XL_HumpId,
	LPXLOPER XL_FactorsNb,
	LPXLOPER XL_TimeStepsNb,	
	LPXLOPER XL_GridSize,		
	LPXLOPER XL_StdDevNb,
	LPXLOPER XL_SkipPDE,
	LPXLOPER XL_SwitchTheta,
	LPXLOPER XL_AllowInterpol,
	LPXLOPER XL_CalibProxy,
	LPXLOPER XL_Recorrel)
{
	ADD_LOG("Local_PXL_SmiledMarketModel_Create");
	bool PersistentInXL = false;
	return Local_SmiledMarketModelCommon(
		XL_Pattern,
        XL_ZeroCurveId,
		XL_CorrelParamId,
		XL_HumpId,
		XL_FactorsNb,
		XL_TimeStepsNb,	
		XL_GridSize,		
		XL_StdDevNb,
		XL_SkipPDE,
		XL_SwitchTheta,
		XL_AllowInterpol,
		XL_CalibProxy,
		XL_Recorrel,
		PersistentInXL );
}

///----------------------------------------------
///----------------------------------------------
///             SmiledFRM Model
/// Inputs :
///     Zc curve Id
///     Hump Param Id
///     Correlation Param Id
///		FactorsNb
///----------------------------------------------
///----------------------------------------------


/////////////////////////////////////////////////////////////
/// central function that does the creation of the XL function
/////////////////////////////////////////////////////////////
LPXLOPER Local_SmiledMarketModelDSCommon(
	LPXLOPER XL_CalibPattern,
	LPXLOPER XL_StartDate,
	LPXLOPER XL_EndDate,
	LPXLOPER XL_ResetFreq,
	LPXLOPER XL_IndexFreq,
	LPXLOPER XL_IndexType,	
	LPXLOPER XL_ResetTiming,		
	LPXLOPER XL_DayCount,
	LPXLOPER XL_ResetCalendar,
	LPXLOPER XL_FwdRule,
	LPXLOPER XL_IntRule,
	LPXLOPER XL_StubRule,
	LPXLOPER XL_ResetGap,
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

		CCString C_CalibPatternStr;
		XL_readStrCellWD(XL_CalibPattern,C_CalibPatternStr,"LIBOR"," ARM_ERR: Param Type: String expected",C_result);
		string C_CalibPattern = CCSTringToSTLString(C_CalibPatternStr);

		double C_indexTypeDbl;
		XL_readNumCell( XL_IndexType, C_indexTypeDbl, " ARM_ERR: Nb of factors: numeric expected",C_result);
		int C_indexType = C_indexTypeDbl;

		CCString C_resetCalendar;
		XL_readStrCellWD(XL_ResetCalendar, C_resetCalendar, GETDEFAULTVALUESTR," ARM_ERR: reset Calendar: string expected",C_result);
		
		double C_startDate;
		double C_endDate;
		long C_resetFreq;
		long C_indexFreq;
		long C_dayCount;
		long C_fwdRule;
		long C_intRule;
		long C_stubRule;
		double C_resetGapDbl;
		long C_resetTiming;

		double gapDefault = GETDEFAULTVALUE;
		XL_readNumCellWD(XL_ResetGap,C_resetGapDbl,gapDefault, " ARM_ERR: reset gap: numeric expected",	C_result);
		long C_resetGap = C_resetGapDbl;

		
		XL_readNumCell(XL_StartDate,C_startDate,			" ARM_ERR: start Date: date expected",		C_result );
		XL_readNumCell(XL_EndDate,C_endDate,				" ARM_ERR: end Date: date expected",		C_result );
		XL_GETFREQUENCYWD(XL_ResetFreq,C_resetFreq, "A",	" ARM_ERR: reset freq: string expected",	C_result );
		XL_GETFREQUENCYWD(XL_IndexFreq,C_indexFreq, "A",	" ARM_ERR: reset freq: string expected",	C_result );
		XL_GETDAYCOUNTWD( XL_DayCount, C_dayCount,"ACTUAL",	" ARM_ERR: dayCount : string expected",		C_result );
		XL_GETFWDRULEWD( XL_FwdRule, C_fwdRule, "MF",		" ARM_ERR: fwdRule : string expected",		C_result );
		XL_GETINTRULEWD( XL_IntRule, C_intRule, "ADJ",		" ARM_ERR: fwdRule : string expected",		C_result );
		XL_GETCONVRULEWD( XL_StubRule, C_stubRule, "SS",	" ARM_ERR: stub Rule : string expected",	C_result );
		XL_GETPAYRESETTIMINGWD(XL_ResetTiming,C_resetTiming,"ADV","ARM_ERR: reset timing: string expected",	C_result );
	
		exportFunc13Args<string, double, double ,long ,long ,int ,long ,long ,CCString, long ,long ,long ,long>  ourFunc(
			C_CalibPattern,
			C_startDate,
			C_endDate,
			C_resetFreq,
			C_indexFreq,
			C_indexType,
			C_resetTiming,
			C_dayCount,
			C_resetCalendar,
			C_fwdRule,
			C_intRule,
			C_stubRule,
			C_resetGap,
			ARMLOCAL_SmiledMarketModelDS_Create);

		/// call the general function
		fillXL_Result( LOCAL_DATESTRIP_CLASS, ourFunc, C_result, XL_result, PersistentInXL );
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_SmiledMarketModelDSCommon" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}

///////////////////////////////////
/// version that takes into account 
/// previous creation of object
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_SmiledMarketModelDS_Create(
	LPXLOPER XL_CalibPattern,
	LPXLOPER XL_StartDate,
	LPXLOPER XL_EndDate,
	LPXLOPER XL_ResetFreq,
	LPXLOPER XL_IndexFreq,
	LPXLOPER XL_IndexType,	
	LPXLOPER XL_ResetTiming,		
	LPXLOPER XL_DayCount,
	LPXLOPER XL_ResetCalendar,
	LPXLOPER XL_FwdRule,
	LPXLOPER XL_IntRule,
	LPXLOPER XL_StubRule,
	LPXLOPER XL_ResetGap)
{
	ADD_LOG("Local_SmiledMarketModelDS_Create");
	bool PersistentInXL = true;
	return Local_SmiledMarketModelDSCommon(
		XL_CalibPattern,
		XL_StartDate,
		XL_EndDate,
		XL_ResetFreq,
		XL_IndexFreq,
		XL_IndexType,	
		XL_ResetTiming,		
		XL_DayCount,
		XL_ResetCalendar,
		XL_FwdRule,
		XL_IntRule,
		XL_StubRule,
		XL_ResetGap,
		PersistentInXL );
}

///////////////////////////////////
/// version for non persistentInXl
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_PXL_SmiledMarketModelDS_Create(
	LPXLOPER XL_CalibPattern,
	LPXLOPER XL_StartDate,
	LPXLOPER XL_EndDate,
	LPXLOPER XL_ResetFreq,
	LPXLOPER XL_IndexFreq,
	LPXLOPER XL_IndexType,	
	LPXLOPER XL_ResetTiming,		
	LPXLOPER XL_DayCount,
	LPXLOPER XL_ResetCalendar,
	LPXLOPER XL_FwdRule,
	LPXLOPER XL_IntRule,
	LPXLOPER XL_StubRule,
	LPXLOPER XL_ResetGap)
{
	ADD_LOG("Local_PXL_SmiledMarketModelDS_Create");
	bool PersistentInXL = false;
	return Local_SmiledMarketModelDSCommon(
		XL_CalibPattern,
		XL_StartDate,
		XL_EndDate,
		XL_ResetFreq,
		XL_IndexFreq,
		XL_IndexType,	
		XL_ResetTiming,		
		XL_DayCount,
		XL_ResetCalendar,
		XL_FwdRule,
		XL_IntRule,
		XL_StubRule,
		XL_ResetGap,
		PersistentInXL );
}

///----------------------------------------------
///----------------------------------------------
///             SmiledFX Model
/// Inputs :
///     Zc curve Id
///     Hump Param Id
///     Correlation Param Id
///		FactorsNb
///----------------------------------------------
///----------------------------------------------


/////////////////////////////////////////////////////////////
/// central function that does the creation of the XL function
/////////////////////////////////////////////////////////////
LPXLOPER Local_SmiledFXModelCommon(
	LPXLOPER XL_DomZeroCurveId,
	LPXLOPER XL_ForZeroCurveId,
	LPXLOPER XL_FXSpot,
	LPXLOPER XL_CorrelParamId,
	LPXLOPER XL_HumpId,
	LPXLOPER XL_FactorsNb,
	LPXLOPER XL_TimeStepsNb,	
	LPXLOPER XL_GridSize,		
	LPXLOPER XL_StdDevNb,
	LPXLOPER XL_SkipPDE,
	LPXLOPER XL_CorrelType,
	LPXLOPER XL_Recorrel,
	LPXLOPER XL_Rescalling,
	LPXLOPER XL_Model2IRFXId,
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

		long C_DomZeroCurveId;
		XL_GETOBJID( XL_DomZeroCurveId, C_DomZeroCurveId,	" ARM_ERR: Domestic Zero Curve: Object expected",C_result);

		long C_ForZeroCurveId;
		XL_GETOBJID( XL_ForZeroCurveId, C_ForZeroCurveId,	" ARM_ERR: Foreign Zero Curve: Object expected",C_result);

		double C_FXSpot;
		XL_readNumCell( XL_FXSpot, C_FXSpot,	" ARM_ERR: FX Spot: numeric expected",C_result);
		
		long C_CorrelParamId;
		XL_GETOBJID( XL_CorrelParamId, C_CorrelParamId,	" ARM_ERR: CorrelParam: Object expected",C_result);
		
		long C_HumpId;
		XL_GETOBJID( XL_HumpId, C_HumpId,	" ARM_ERR: Hump: Object expected",C_result);

		double C_FactorsNb;
		XL_readNumCell( XL_FactorsNb, C_FactorsNb, " ARM_ERR: Nb of factors: numeric expected",C_result);
		int C_FactorsNbInt = C_FactorsNb;

		double C_TimeStepsNb;
		double C_TimeStepsNbDefault=1000;
		XL_readNumCellWD(XL_TimeStepsNb, C_TimeStepsNb, C_TimeStepsNbDefault, " ARM_ERR: dim : numeric expected", C_result);
		int C_TimeStepsNbInt = C_TimeStepsNb;

		double C_GridSize;
		double C_GridSizeDefault=801;
		XL_readNumCellWD(XL_GridSize, C_GridSize, C_GridSizeDefault, " ARM_ERR: GridSize : numeric expected", C_result);
		int C_GridSizeInt = C_GridSize;

		double C_StdDevNb;
		double C_StdDevNbDefault=6;
		XL_readNumCellWD(XL_StdDevNb, C_StdDevNb, C_StdDevNbDefault, " ARM_ERR: StdDevNb : numeric expected", C_result);
		
		double C_Recorrel;
		double C_RecorrelDefault=0;
		XL_readNumCellWD(XL_Recorrel, C_Recorrel, C_RecorrelDefault, " ARM_ERR: Recorrel : numeric expected", C_result);
		
		CCString C_SkipPDE;
		XL_readStrCellWD( XL_SkipPDE,C_SkipPDE,"Y"," ARM_ERR: Skip PDE: string expected",C_result);
		bool C_SkipPDEBool;
		C_SkipPDE.toUpper();
		if (C_SkipPDE == "Y" )
			C_SkipPDEBool = true;
		else if (C_SkipPDE == "N")
			C_SkipPDEBool = false;
		else
			throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"String \"Y\" or \"N\" Expected for Skip PDE");

		CCString C_CorrelTypeStr;
		XL_readStrCellWD(XL_CorrelType,C_CorrelTypeStr,"Beta"," ARM_ERR: Param Type: String expected",C_result);
		string C_CorrelType = CCSTringToSTLString(C_CorrelTypeStr);

		CCString C_Rescalling;
		XL_readStrCellWD( XL_Rescalling,C_Rescalling,"N"," ARM_ERR: Rescalling: string expected",C_result);
		bool C_RescallingBool;
		C_Rescalling.toUpper();
		if (C_Rescalling == "Y" )
			C_RescallingBool = true;
		else if (C_Rescalling == "N")
			C_RescallingBool = false;
		else
			throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"String \"Y\" or \"N\" Expected for Rescalling");

		long C_Model2IRFXId;
		XL_GETOBJIDWD( XL_Model2IRFXId, C_Model2IRFXId,	"NULL OBJECT",	" ARM_ERR: Model 2IRFX: object expected",	C_result);

		exportFunc14Args< long, long,  double, long ,long, int, int ,int ,double ,bool, string, double, bool, long>  ourFunc(
			C_DomZeroCurveId,
			C_ForZeroCurveId,
			C_FXSpot,
			C_CorrelParamId,
			C_HumpId,
			C_FactorsNbInt,
			C_TimeStepsNbInt,
			C_GridSizeInt,
			C_StdDevNb,
			C_SkipPDEBool,
			C_CorrelType,
			C_Recorrel,
			C_RescallingBool,
			C_Model2IRFXId,
			ARMLOCAL_SmiledFXModel_Create);

		/// call the general function
		fillXL_Result( LOCAL_SMILEDFRMMODEL_CLASS, ourFunc, C_result, XL_result, PersistentInXL );
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_SmiledFRMModelCommon" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}

///////////////////////////////////
/// version that takes into account 
/// previous creation of object
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_SmiledFXModel_Create(
	LPXLOPER XL_DomZeroCurveId,
	LPXLOPER XL_ForZeroCurveId,
	LPXLOPER XL_FXSpot,
	LPXLOPER XL_CorrelParamId,
	LPXLOPER XL_HumpId,
	LPXLOPER XL_FactorsNb,
	LPXLOPER XL_TimeStepsNb,	
	LPXLOPER XL_GridSize,		
	LPXLOPER XL_StdDevNb,
	LPXLOPER XL_SkipPDE,
	LPXLOPER XL_CorrelType,
	LPXLOPER XL_Recorrel,
	LPXLOPER XL_Rescalling,
	LPXLOPER XL_Model2IRFXId)
{
	ADD_LOG("Local_SmiledFXModel_Create");
	bool PersistentInXL = true;
	return Local_SmiledFXModelCommon(
        XL_DomZeroCurveId,
		XL_ForZeroCurveId,
		XL_FXSpot,
		XL_CorrelParamId,
		XL_HumpId,
		XL_FactorsNb,
		XL_TimeStepsNb,	
		XL_GridSize,		
		XL_StdDevNb,
		XL_SkipPDE,
		XL_CorrelType,
		XL_Recorrel,
		XL_Rescalling,
		XL_Model2IRFXId,
		PersistentInXL );
}

///////////////////////////////////
/// version for non persistentInXl
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_PXL_SmiledFXModel_Create(
	LPXLOPER XL_DomZeroCurveId,
	LPXLOPER XL_ForZeroCurveId,
	LPXLOPER XL_FXSpot,
	LPXLOPER XL_CorrelParamId,
	LPXLOPER XL_HumpId,
	LPXLOPER XL_FactorsNb,
	LPXLOPER XL_TimeStepsNb,	
	LPXLOPER XL_GridSize,		
	LPXLOPER XL_StdDevNb,
	LPXLOPER XL_SkipPDE,
	LPXLOPER XL_CorrelType,
	LPXLOPER XL_Recorrel,
	LPXLOPER XL_Rescalling,
	LPXLOPER XL_Model2IRFXId)
{
	ADD_LOG("Local_PXL_SmiledFXModel_Create");
	bool PersistentInXL = false;
	return Local_SmiledFXModelCommon(
        XL_DomZeroCurveId,
		XL_ForZeroCurveId,
		XL_FXSpot,
		XL_CorrelParamId,
		XL_HumpId,
		XL_FactorsNb,
		XL_TimeStepsNb,	
		XL_GridSize,		
		XL_StdDevNb,
		XL_SkipPDE,
		XL_CorrelType,
		XL_Recorrel,
		XL_Rescalling,
		XL_Model2IRFXId,
		PersistentInXL );
}

///----------------------------------------------
///----------------------------------------------
///             SmiledFX Model
/// Inputs :
///     Zc curve Id
///     Hump Param Id
///     Correlation Param Id
///		FactorsNb
///----------------------------------------------
///----------------------------------------------


/////////////////////////////////////////////////////////////
/// central function that does the creation of the XL function
/////////////////////////////////////////////////////////////
LPXLOPER Local_MixtureFXModel_Calibrate(
	LPXLOPER XL_Fwd,
	LPXLOPER XL_Expiry,
	LPXLOPER XL_CallPut,
	LPXLOPER XL_Strikes,
	LPXLOPER XL_Vols,
	LPXLOPER XL_DecVol,
	LPXLOPER XL_Alpha,	
	LPXLOPER XL_Lambda)
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

		double C_Fwd;
		XL_readNumCell(XL_Fwd, C_Fwd, " ARM_ERR: Fwd : numeric expected", C_result);

		double C_Expiry;
		XL_readNumCell(XL_Expiry, C_Expiry, " ARM_ERR: Expiry : numeric expected", C_result);

		CCString C_CallPutStr;
		long C_CallPut;
		XL_readStrCell(XL_CallPut, C_CallPutStr, " ARM_ERR: CallPut : string expected", C_result);
		C_CallPut = ARM_ConvCapOrFloor(C_CallPutStr,C_result);

		vector<double> C_Strikes;
		XL_readNumVector(XL_Strikes,C_Strikes," ARM_ERR: Strikes: array of numeric expected",C_result);

		vector<double> C_Vols;
		XL_readNumVector(XL_Vols,C_Vols," ARM_ERR: Vols: array of numeric expected",C_result);

		vector<double> C_DecVol;
		XL_readNumVector(XL_DecVol,C_DecVol," ARM_ERR: DecVol: array of numeric expected",C_result);

		vector<double> C_Alpha;
		XL_readNumVector(XL_Alpha,C_Alpha," ARM_ERR: Alpha: array of numeric expected",C_result);

		vector<double> C_Lambda;
		XL_readNumVector(XL_Lambda,C_Lambda," ARM_ERR: Lambda: array of numeric expected",C_result);

		vector<double> outParams;

		long retCode = ARMLOCAL_MixtureFXModel_Calibrate(
			C_Fwd,
			C_Expiry,
			C_CallPut,
			C_Strikes,
			C_Vols,
			C_DecVol,
			C_Alpha,
			C_Lambda,
			outParams,
			C_result);

		if (retCode == ARM_OK)
		{
			XL_writeNumVector( XL_result, outParams, " ARM_ERR: Could not get result data", C_result);
		}
	
		/// be ware that ARM_ERR is a macro
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_SmiledFRMModelCommon" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


///----------------------------------------------
///----------------------------------------------
///             SVBGM Model
/// Inputs :
///     Zc curve Id
///     Hump Param Id
///     Correlation Param Id
///		FactorsNb
///----------------------------------------------
///----------------------------------------------


/////////////////////////////////////////////////////////////
/// central function that does the creation of the XL function
/////////////////////////////////////////////////////////////
LPXLOPER Local_SVBGMModelCommon(
	LPXLOPER XL_ZeroCurveId,
	LPXLOPER XL_ShiftId,
	LPXLOPER XL_AlphaId,
	LPXLOPER XL_NuId,
	LPXLOPER XL_RhoId,
	LPXLOPER XL_RRCorrelParamId,
	LPXLOPER XL_RVCorrelParamId,
	LPXLOPER XL_VVCorrelParamId,
	LPXLOPER XL_Recorrel,
	LPXLOPER XL_FactorsNb,
	LPXLOPER XL_MinRatio,
	LPXLOPER XL_Proxy,
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

		long C_ZeroCurveId;
		XL_GETOBJID( XL_ZeroCurveId, C_ZeroCurveId,	" ARM_ERR: Zero Curve: Object expected",C_result);
		
		long C_ShiftId;
		XL_GETOBJID( XL_ShiftId, C_ShiftId,	" ARM_ERR: shift : Object expected",C_result);
			
		long C_AlphaId;
		XL_GETOBJID( XL_AlphaId, C_AlphaId,	" ARM_ERR: Alpha : Object expected",C_result);

		long C_RhoId;
		XL_GETOBJID( XL_RhoId, C_RhoId,	" ARM_ERR: Rho : Object expected",C_result);
		
		long C_NuId;
		XL_GETOBJID( XL_NuId, C_NuId,	" ARM_ERR: Nu : Object expected",C_result);

		long C_RRCorrelParamId;
		XL_GETOBJID( XL_RRCorrelParamId, C_RRCorrelParamId,	" ARM_ERR: RRCorrelParam : Object expected",C_result);

		long C_RVCorrelParamId;
		XL_GETOBJID( XL_RVCorrelParamId, C_RVCorrelParamId,	" ARM_ERR: RVCorrelParam : Object expected",C_result);

		long C_VVCorrelParamId;
		XL_GETOBJID( XL_VVCorrelParamId, C_VVCorrelParamId,	" ARM_ERR: VVCorrelParam : Object expected",C_result);
		
		double C_Recorrel;
		double C_RecorrelDefault=0.;
		XL_readNumCellWD(XL_Recorrel, C_Recorrel, C_RecorrelDefault, " ARM_ERR: Recorrel : numeric expected", C_result);
		
		double C_FactorsNb;
		double C_FactorsNbDefault = 0.;
		XL_readNumCellWD( XL_FactorsNb, C_FactorsNb, C_FactorsNbDefault, " ARM_ERR: Nb of factors: numeric expected",C_result);
		int C_FactorsNbInt = C_FactorsNb;

		double C_MinRatio;
		double C_MinRatioDefault=1.;
		XL_readNumCellWD(XL_MinRatio, C_MinRatio, C_MinRatioDefault, " ARM_ERR: MinRatio : numeric expected", C_result);

		CCString C_Proxy;
		XL_readStrCellWD( XL_Proxy,C_Proxy,"N"," ARM_ERR: Proxy: string expected",C_result);
		bool C_ProxyBool;
		C_Proxy.toUpper();
		if (C_Proxy == "Y" )
			C_ProxyBool = true;
		else
			C_ProxyBool = false;

		exportFunc12Args< long, long ,long ,long, long, long, long, long, double, int, double, bool> ourFunc(
			C_ZeroCurveId,
			C_ShiftId,
			C_AlphaId,
			C_NuId,
			C_RhoId,
			C_RRCorrelParamId,
			C_RVCorrelParamId,
			C_VVCorrelParamId,
			C_Recorrel,
			C_FactorsNbInt,
			C_MinRatio,
			C_ProxyBool,
			ARMLOCAL_SVBGMModel_Create);

		/// call the general function
		fillXL_Result( LOCAL_SVBGMMODEL_CLASS, ourFunc, C_result, XL_result, PersistentInXL );
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_SVBGMModelCommon" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}

///////////////////////////////////
/// version that takes into account 
/// previous creation of object
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_SVBGMModel_Create(
	LPXLOPER XL_zeroCurveId,
	LPXLOPER XL_shiftId,
	LPXLOPER XL_alphaId,
	LPXLOPER XL_nuId,
	LPXLOPER XL_rhoId,
	LPXLOPER XL_rrcorrelParamId,
	LPXLOPER XL_rvcorrelParamId,
	LPXLOPER XL_vvcorrelParamId,
	LPXLOPER XL_recorrel,
	LPXLOPER XL_factorsNb,
	LPXLOPER XL_minratio,
	LPXLOPER XL_Proxy)
{
	ADD_LOG("Local_SVBGMModel_Create");
	bool PersistentInXL = true;
	return Local_SVBGMModelCommon(
		XL_zeroCurveId,
		XL_shiftId,
		XL_alphaId,
		XL_nuId,
		XL_rhoId,
		XL_rrcorrelParamId,
		XL_rvcorrelParamId,
		XL_vvcorrelParamId,
		XL_recorrel,
		XL_factorsNb,
		XL_minratio,
		XL_Proxy,
		PersistentInXL );
}

///////////////////////////////////
/// version for non persistentInXl
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_PXL_SVBGMModel_Create(
	LPXLOPER XL_zeroCurveId,
	LPXLOPER XL_shiftId,
	LPXLOPER XL_alphaId,
	LPXLOPER XL_nuId,
	LPXLOPER XL_rhoId,
	LPXLOPER XL_rrcorrelParamId,
	LPXLOPER XL_rvcorrelParamId,
	LPXLOPER XL_vvcorrelParamId,
	LPXLOPER XL_recorrel,
	LPXLOPER XL_factorsNb,
	LPXLOPER XL_minratio,
	LPXLOPER XL_Proxy)
{
	ADD_LOG("Local_PXL_SVBGMModel_Create");
	bool PersistentInXL = false;
	return Local_SVBGMModelCommon(
		XL_zeroCurveId,
		XL_shiftId,
		XL_alphaId,
		XL_nuId,
		XL_rhoId,
		XL_rrcorrelParamId,
		XL_rvcorrelParamId,
		XL_vvcorrelParamId,
		XL_recorrel,
		XL_factorsNb,
		XL_minratio,
		XL_Proxy,
		PersistentInXL );
}

LPXLOPER Local_BGMSV1FModelCommon(
	LPXLOPER XL_ZeroCurveId,
	LPXLOPER XL_ShiftId,
	LPXLOPER XL_LevelId,
	LPXLOPER XL_InitVarId,
	LPXLOPER XL_LongTermVarId,
	LPXLOPER XL_VarVolId,
	LPXLOPER XL_VarMeanRevId,
	LPXLOPER XL_RhoId,
	LPXLOPER XL_RRCorrelParamId,
	LPXLOPER XL_LocalCalibration,
	LPXLOPER XL_Recorrel,
	LPXLOPER XL_FactorsNb,
	LPXLOPER XL_MinRatio,
	LPXLOPER XL_localRhoCalib,
	LPXLOPER XL_StdDevForCalib,
	LPXLOPER XL_Proxy,
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

		long C_ZeroCurveId;
		XL_GETOBJID( XL_ZeroCurveId, C_ZeroCurveId,	" ARM_ERR: Zero Curve: Object expected",C_result);
		
		long C_ShiftId;
		XL_GETOBJID( XL_ShiftId, C_ShiftId,	" ARM_ERR: shift : Object expected",C_result);
			
		long C_LevelId;
		XL_GETOBJID( XL_LevelId, C_LevelId,	" ARM_ERR: level : Object expected",C_result);

		long C_RhoId;
		XL_GETOBJID( XL_RhoId, C_RhoId,	" ARM_ERR: Rho : Object expected",C_result);
		
		long C_InitVarId;
		XL_GETOBJID( XL_InitVarId, C_InitVarId,	" ARM_ERR: initial var : Object expected",C_result);

		long C_VarMeanRevId;
		XL_GETOBJID( XL_VarMeanRevId, C_VarMeanRevId," ARM_ERR: var mean rev : Object expected",C_result);

		long C_LongTermVarId;
		XL_GETOBJID( XL_LongTermVarId, C_LongTermVarId,	" ARM_ERR: Long term var : Object expected",C_result);

		long C_VarVolId;
		XL_GETOBJID( XL_VarVolId, C_VarVolId,	" ARM_ERR: var volatilities : Object expected",C_result);

		long C_RRCorrelParamId;
		XL_GETOBJID( XL_RRCorrelParamId, C_RRCorrelParamId,	" ARM_ERR: RRCorrelParam : Object expected",C_result);

		double C_Recorrel;
		double C_RecorrelDefault=0.;
		XL_readNumCellWD(XL_Recorrel, C_Recorrel, C_RecorrelDefault, " ARM_ERR: Recorrel : numeric expected", C_result);
		
		double C_FactorsNb;
		double C_FactorsNbDefault = 0.;
		XL_readNumCellWD( XL_FactorsNb, C_FactorsNb, C_FactorsNbDefault, " ARM_ERR: Nb of factors: numeric expected",C_result);
		int C_FactorsNbInt = C_FactorsNb;

		double C_MinRatio;
		double C_MinRatioDefault=1.;
		XL_readNumCellWD(XL_MinRatio, C_MinRatio, C_MinRatioDefault, " ARM_ERR: MinRatio : numeric expected", C_result);

		CCString C_locrhocalib;
		XL_readStrCellWD( XL_localRhoCalib,C_locrhocalib,"Y"," ARM_ERR: Global Calib: string expected",C_result);
		C_locrhocalib.toUpper();
		bool C_locrhocalibBool = C_locrhocalib == "Y" ? true : false;

		CCString C_globalcalib;
		XL_readStrCell(XL_LocalCalibration, C_globalcalib, "ARM_ERR : Local Calib : string expected", C_result);
		C_globalcalib.toUpper();
		bool C_globalcalibBool = C_globalcalib == "Y" ? false : true;

		vector<double> C_Stddev(0);
		XL_getNumVector(XL_StdDevForCalib, C_Stddev);

		CCString C_Proxy;
		XL_readStrCellWD( XL_Proxy,C_Proxy,"N"," ARM_ERR: Proxy: string expected",C_result);
		bool C_ProxyBool;
		C_Proxy.toUpper();
		if (C_Proxy == "Y" )
			C_ProxyBool = true;
		else
			C_ProxyBool = false;

		exportFunc16Args< long, long ,long ,long, long, long, long, long, long, bool, double, int, double, bool, vector<double>, bool> ourFunc(
			C_ZeroCurveId,
			C_ShiftId,
			C_LevelId,
			C_InitVarId,
			C_LongTermVarId,
			C_VarVolId,
			C_VarMeanRevId,
			C_RhoId,
			C_RRCorrelParamId,
			C_globalcalibBool,
			C_Recorrel,
			C_FactorsNb,
			C_MinRatio,
			C_locrhocalibBool,
			C_Stddev,
			C_ProxyBool,
			ARMLOCAL_BGMSV1FModel_Create);

		/// call the general function
		fillXL_Result( LOCAL_BGMSV1FMODEL_CLASS, ourFunc, C_result, XL_result, PersistentInXL );
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_SVBGMModelCommon" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}

__declspec(dllexport) LPXLOPER WINAPI Local_BGMSV1FModel_Create(
	LPXLOPER XL_zeroCurveId,
	LPXLOPER XL_shiftId,
	LPXLOPER XL_levelId,
	LPXLOPER XL_initVar,
	LPXLOPER XL_LongTermVarId,
	LPXLOPER XL_VarVolId,
	LPXLOPER XL_varMeanRev,
	LPXLOPER XL_rhoId,
	LPXLOPER XL_rrcorrelParamId,
	LPXLOPER XL_LocalCalibration,
	LPXLOPER XL_recorrel,
	LPXLOPER XL_factorsNb,
	LPXLOPER XL_minratio,
	LPXLOPER XL_localRhoCalib,
	LPXLOPER XL_StdDevForCalib,
	LPXLOPER XL_Proxy)
{
	ADD_LOG("Local_BGMSV1FModel_Create");
	bool PersistentInXL = true;

	return Local_BGMSV1FModelCommon(
		XL_zeroCurveId,
		XL_shiftId,
		XL_levelId,
		XL_initVar,
		XL_LongTermVarId,
		XL_VarVolId,
		XL_varMeanRev,
		XL_rhoId,
		XL_rrcorrelParamId,
		XL_LocalCalibration,
		XL_recorrel,
		XL_factorsNb,
		XL_minratio,
		XL_localRhoCalib,
		XL_StdDevForCalib,
		XL_Proxy,
		PersistentInXL );
}


__declspec(dllexport) LPXLOPER WINAPI Local_PXL_BGMSV1FModel_Create(
	LPXLOPER XL_zeroCurveId,
	LPXLOPER XL_shiftId,
	LPXLOPER XL_levelId,
	LPXLOPER XL_initVar,
	LPXLOPER XL_LongTermVarId,
	LPXLOPER XL_VarVolId,
	LPXLOPER XL_varMeanRev,
	LPXLOPER XL_rhoId,
	LPXLOPER XL_rrcorrelParamId,
	LPXLOPER XL_LocalCalibration,
	LPXLOPER XL_recorrel,
	LPXLOPER XL_factorsNb,
	LPXLOPER XL_minratio,
	LPXLOPER XL_localRhoCalib,
	LPXLOPER XL_StdDevForCalib,
	LPXLOPER XL_Proxy)
{
	ADD_LOG("Local_PXL_BGMSV1FModel_Create");
	bool PersistentInXL = false;

	return Local_BGMSV1FModelCommon(
		XL_zeroCurveId,
		XL_shiftId,
		XL_levelId,
		XL_initVar,
		XL_LongTermVarId,
		XL_VarVolId,
		XL_varMeanRev,
		XL_rhoId,
		XL_rrcorrelParamId,
		XL_LocalCalibration,
		XL_recorrel,
		XL_factorsNb,
		XL_minratio,
		XL_localRhoCalib,
		XL_StdDevForCalib,
		XL_Proxy,
		PersistentInXL );
}

__declspec(dllexport) LPXLOPER WINAPI Local_BGMSV2FModel_CreateCommon(
	LPXLOPER XL_ZeroCurveId,
	LPXLOPER XL_BetaCorrelId,
	LPXLOPER XL_V01,
	LPXLOPER XL_Kappa1,
	LPXLOPER XL_V02,
	LPXLOPER XL_Kappa2,
	LPXLOPER XL_Rho1,
	LPXLOPER XL_Rho2,
	LPXLOPER XL_LocalRho1Calib,
	LPXLOPER XL_LocalRho2Calib,
	LPXLOPER XL_Shift,
	LPXLOPER XL_Recorrel,
	LPXLOPER XL_FactorsNb,
	LPXLOPER XL_MinRatio,
	LPXLOPER XL_StdDevForCalib,
	LPXLOPER XL_Proxy,
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

		long C_ZeroCurveId;
		XL_GETOBJID( XL_ZeroCurveId, C_ZeroCurveId,	" ARM_ERR: Zero Curve: Object expected",C_result);
		
		long C_BetaCorrelId;
		XL_GETOBJID( XL_BetaCorrelId, C_BetaCorrelId,	" ARM_ERR: BetaCorrel : Object expected",C_result);

		double v01, v02, kappa1, kappa2, rho1, rho2, shift, defRho(-2.), defShift(0.);

		XL_readNumCell(XL_V01, v01, " ARM_ERR : v01 : numeric expected", C_result);
		XL_readNumCell(XL_V02, v02, " ARM_ERR : v02 : numeric expected", C_result);
		XL_readNumCell(XL_Kappa1, kappa1, " ARM_ERR : kappa1 : numeric expected", C_result);
		XL_readNumCell(XL_Kappa2, kappa2, " ARM_ERR : kappa2 : numeric expected", C_result);

		XL_readNumCellWD(XL_Rho1, rho1, defRho, " ARM_ERR : rho1 : numeric expected", C_result);
		XL_readNumCellWD(XL_Rho2, rho2, defRho, " ARM_ERR : rho2 : numeric expected", C_result);
		XL_readNumCellWD(XL_Shift, shift, defShift, " ARM_ERR : shift: numeric expected", C_result);

		double C_Recorrel;
		double C_RecorrelDefault=0.;
		XL_readNumCellWD(XL_Recorrel, C_Recorrel, C_RecorrelDefault, " ARM_ERR: Recorrel : numeric expected", C_result);
		
		double C_FactorsNb;
		double C_FactorsNbDefault = 0.;
		XL_readNumCellWD( XL_FactorsNb, C_FactorsNb, C_FactorsNbDefault, " ARM_ERR: Nb of factors: numeric expected",C_result);
		int C_FactorsNbInt = C_FactorsNb;

		double C_MinRatio;
		double C_MinRatioDefault=1.;
		XL_readNumCellWD(XL_MinRatio, C_MinRatio, C_MinRatioDefault, " ARM_ERR: MinRatio : numeric expected", C_result);

		CCString C_locrho1calib;
		XL_readStrCellWD( XL_LocalRho1Calib,C_locrho1calib,"N"," ARM_ERR: Global Calib: string expected",C_result);
		C_locrho1calib.toUpper();
		bool C_locrho1calibBool = C_locrho1calib == "Y" ? true : false;

		CCString C_locrho2calib;
		XL_readStrCellWD( XL_LocalRho2Calib,C_locrho2calib,"N"," ARM_ERR: Global Calib: string expected",C_result);
		C_locrho2calib.toUpper();
		bool C_locrho2calibBool = C_locrho2calib == "Y" ? true : false;

		vector<double> C_Stddev(0);
		XL_getNumVector(XL_StdDevForCalib, C_Stddev);

		CCString C_Proxy;
		XL_readStrCellWD( XL_Proxy,C_Proxy,"N"," ARM_ERR: Proxy: string expected",C_result);
		bool C_ProxyBool;
		C_Proxy.toUpper();
		if (C_Proxy == "Y" )
			C_ProxyBool = true;
		else
			C_ProxyBool = false;

		exportFunc16Args< long, long ,double, int, double, double, double, double, double, double, double, double, bool, bool, vector<double>, bool> ourFunc(
			C_ZeroCurveId,
			C_BetaCorrelId,
			C_Recorrel,
			C_FactorsNb,
			C_MinRatio,
			v01, 
			kappa1,
			rho1,
			v02,
			kappa2,
			rho2,
			shift,
			C_locrho1calibBool,
			C_locrho2calibBool,
			C_Stddev,
			C_ProxyBool,
			ARMLOCAL_BGMSV2FModel_Create);

		/// call the general function
		fillXL_Result( LOCAL_BGMSV1FMODEL_CLASS, ourFunc, C_result, XL_result, PersistentInXL );
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_SVBGMModelCommon" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}

__declspec(dllexport) LPXLOPER WINAPI Local_BGMSV2FModel_Create(
	LPXLOPER XL_ZeroCurveId,
	LPXLOPER XL_BetaCorrelId,
	LPXLOPER XL_V01,
	LPXLOPER XL_Kappa1,
	LPXLOPER XL_V02,
	LPXLOPER XL_Kappa2,
	LPXLOPER XL_Rho1,
	LPXLOPER XL_Rho2,
	LPXLOPER XL_LocalRho1Calib,
	LPXLOPER XL_LocalRho2Calib,
	LPXLOPER XL_Shift,
	LPXLOPER XL_recorrel,
	LPXLOPER XL_factorsNb,
	LPXLOPER XL_minratio,
	LPXLOPER XL_StdDevForCalib,
	LPXLOPER XL_Proxy)
{
	ADD_LOG("Local_BGMSV2FModel_Create");

	bool PersistentInXL = false;

	return Local_BGMSV2FModel_CreateCommon(XL_ZeroCurveId, XL_BetaCorrelId,
		XL_V01, XL_Kappa1, 
		XL_V02, XL_Kappa2,
		XL_Rho1, XL_Rho2,
		XL_LocalRho1Calib, XL_LocalRho2Calib,
		XL_Shift, XL_recorrel,
		XL_factorsNb, XL_minratio, XL_StdDevForCalib, XL_Proxy,
		PersistentInXL);
}

__declspec(dllexport) LPXLOPER WINAPI Local_PXL_BGMSV2FModel_Create(
	LPXLOPER XL_ZeroCurveId,
	LPXLOPER XL_BetaCorrelId,
	LPXLOPER XL_V01,
	LPXLOPER XL_Kappa1,
	LPXLOPER XL_V02,
	LPXLOPER XL_Kappa2,
	LPXLOPER XL_Rho1,
	LPXLOPER XL_Rho2,
	LPXLOPER XL_LocalRho1Calib,
	LPXLOPER XL_LocalRho2Calib,
	LPXLOPER XL_Shift,
	LPXLOPER XL_recorrel,
	LPXLOPER XL_factorsNb,
	LPXLOPER XL_minratio,
	LPXLOPER XL_StdDevForCalib,
	LPXLOPER XL_Proxy)
{
	ADD_LOG("Local_PXL_BGMSV2FModel_Create");

	bool PersistentInXL = true;

	return Local_BGMSV2FModel_CreateCommon(XL_ZeroCurveId, XL_BetaCorrelId,
		XL_V01, XL_Kappa1, 
		XL_V02, XL_Kappa2,
		XL_Rho1, XL_Rho2,
		XL_LocalRho1Calib, XL_LocalRho2Calib,
		XL_Shift, XL_recorrel,
		XL_factorsNb, XL_minratio, XL_StdDevForCalib, XL_Proxy,
		PersistentInXL);
}

LPXLOPER Local_SVMMSpreadModelCommon(
	LPXLOPER XL_ZeroCurveId,
	LPXLOPER XL_LevelId,
	LPXLOPER XL_InitVarId,
	LPXLOPER XL_LongTermVarId,
	LPXLOPER XL_VarVolId,
	LPXLOPER XL_VarMeanRevId,
	LPXLOPER XL_RhoId,
	LPXLOPER XL_RRCorrelParamId,
	LPXLOPER XL_Recorrel,
	LPXLOPER XL_FactorsNb,
	LPXLOPER XL_MinRatio,
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

		long C_ZeroCurveId;
		XL_GETOBJID( XL_ZeroCurveId, C_ZeroCurveId,	" ARM_ERR: Zero Curve: Object expected",C_result);
		
		long C_LevelId;
		XL_GETOBJID( XL_LevelId, C_LevelId,	" ARM_ERR: level : Object expected",C_result);

		long C_RhoId;
		XL_GETOBJID( XL_RhoId, C_RhoId,	" ARM_ERR: Rho : Object expected",C_result);
		
		long C_InitVarId;
		XL_GETOBJID( XL_InitVarId, C_InitVarId,	" ARM_ERR: initial var : Object expected",C_result);

		long C_VarMeanRevId;
		XL_GETOBJID( XL_VarMeanRevId, C_VarMeanRevId," ARM_ERR: var mean rev : Object expected",C_result);

		long C_LongTermVarId;
		XL_GETOBJID( XL_LongTermVarId, C_LongTermVarId,	" ARM_ERR: Long term var : Object expected",C_result);

		long C_VarVolId;
		XL_GETOBJID( XL_VarVolId, C_VarVolId,	" ARM_ERR: var volatilities : Object expected",C_result);

		long C_RRCorrelParamId;
		XL_GETOBJID( XL_RRCorrelParamId, C_RRCorrelParamId,	" ARM_ERR: RRCorrelParam : Object expected",C_result);

		double C_Recorrel;
		double C_RecorrelDefault=0.;
		XL_readNumCellWD(XL_Recorrel, C_Recorrel, C_RecorrelDefault, " ARM_ERR: Recorrel : numeric expected", C_result);
		
		double C_FactorsNb;
		double C_FactorsNbDefault = 0.;
		XL_readNumCellWD( XL_FactorsNb, C_FactorsNb, C_FactorsNbDefault, " ARM_ERR: Nb of factors: numeric expected",C_result);
		int C_FactorsNbInt = C_FactorsNb;

		double C_MinRatio;
		double C_MinRatioDefault=1.;
		XL_readNumCellWD(XL_MinRatio, C_MinRatio, C_MinRatioDefault, " ARM_ERR: MinRatio : numeric expected", C_result);

		exportFunc11Args< long ,long ,long, long, long, long, long, long, double, int, double> ourFunc(
			C_ZeroCurveId,
			C_LevelId,
			C_InitVarId,
			C_LongTermVarId,
			C_VarVolId,
			C_VarMeanRevId,
			C_RhoId,
			C_RRCorrelParamId,
			C_Recorrel,
			C_FactorsNb,
			C_MinRatio,
			ARMLOCAL_SVMMSpreadModel_Create);

		/// call the general function
		fillXL_Result( LOCAL_BGMSV1FMODEL_CLASS, ourFunc, C_result, XL_result, PersistentInXL );
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_SVBGMModelCommon" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}

__declspec(dllexport) LPXLOPER WINAPI Local_SVMMSpreadModel_Create(
	LPXLOPER XL_zeroCurveId,
	LPXLOPER XL_levelId,
	LPXLOPER XL_initVar,
	LPXLOPER XL_LongTermVarId,
	LPXLOPER XL_VarVolId,
	LPXLOPER XL_varMeanRev,
	LPXLOPER XL_rhoId,
	LPXLOPER XL_rrcorrelParamId,
	LPXLOPER XL_recorrel,
	LPXLOPER XL_factorsNb,
	LPXLOPER XL_minratio)
{
	ADD_LOG("Local_SVMMSpreadModel_Create");
	bool PersistentInXL = true;

	return Local_SVMMSpreadModelCommon(
		XL_zeroCurveId,
		XL_levelId,
		XL_initVar,
		XL_LongTermVarId,
		XL_VarVolId,
		XL_varMeanRev,
		XL_rhoId,
		XL_rrcorrelParamId,
		XL_recorrel,
		XL_factorsNb,
		XL_minratio,
		PersistentInXL );
}


__declspec(dllexport) LPXLOPER WINAPI Local_PXL_SVMMSpreadModel_Create(
	LPXLOPER XL_zeroCurveId,
	LPXLOPER XL_levelId,
	LPXLOPER XL_initVar,
	LPXLOPER XL_LongTermVarId,
	LPXLOPER XL_VarVolId,
	LPXLOPER XL_varMeanRev,
	LPXLOPER XL_rhoId,
	LPXLOPER XL_rrcorrelParamId,
	LPXLOPER XL_recorrel,
	LPXLOPER XL_factorsNb,
	LPXLOPER XL_minratio)
{
	ADD_LOG("Local_PXL_SVMMSpreadModel_Create");
	bool PersistentInXL = false;

	return Local_SVMMSpreadModelCommon(
		XL_zeroCurveId,
		XL_levelId,
		XL_initVar,
		XL_LongTermVarId,
		XL_VarVolId,
		XL_varMeanRev,
		XL_rhoId,
		XL_rrcorrelParamId,
		XL_recorrel,
		XL_factorsNb,
		XL_minratio,
		PersistentInXL );
}

LPXLOPER Local_BiSVMMCommon(
	LPXLOPER XL_ModelNames,
	LPXLOPER XL_ModelIds,
	LPXLOPER XL_corrMatrixId,
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

		vector<CCString> modelNames;
		XL_readStrVector(XL_ModelNames,modelNames," ARM_ERR: Model names : array of string expected",DOUBLE_TYPE,C_result);
		vector< string > C_modelNames(modelNames.size());
		for(size_t i=0;i<modelNames.size();++i)
			C_modelNames[i]=CCSTringToSTLString(modelNames[i]);

		vector<CCString> C_ModelsIdsVec;
		XL_readStrVector(XL_ModelIds,C_ModelsIdsVec," ARM_ERR: Model: array of object expected",DOUBLE_TYPE,C_result);
		size_t size = C_ModelsIdsVec.size();
		vector<long> C_ModelsIdsVecLong(size);
		for(i = 0; i < size; ++i )    
			C_ModelsIdsVecLong[i] = LocalGetNumObjectId(C_ModelsIdsVec[i]); 

		CCString C_CorrelationStrId;
		XL_readStrCell( XL_corrMatrixId, C_CorrelationStrId,	" ARM_ERR: Correlation Id: Object expected",		C_result);
		long C_CorrelationId = LocalGetNumObjectId(C_CorrelationStrId);

		exportFunc3Args< vector<string>, vector<long>, long> ourFunc(
			C_modelNames,
			C_ModelsIdsVecLong,
			C_CorrelationId,
			ARMLOCAL_BiSVMM_Create);

		/// call the general function
		fillXL_Result( LOCAL_BiSVMM_CLASS, ourFunc, C_result, XL_result, PersistentInXL );
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_BiSVMMCommon" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}

__declspec(dllexport) LPXLOPER WINAPI Local_BiSVMM_Create(
	LPXLOPER XL_modelNames,
	LPXLOPER XL_modelIds,
	LPXLOPER XL_corrMatrixId)
{
	ADD_LOG("Local_BiSVMM_Create");
	bool PersistentInXL = true;

	return Local_BiSVMMCommon(XL_modelNames, XL_modelIds, XL_corrMatrixId, PersistentInXL);
}

__declspec(dllexport) LPXLOPER WINAPI Local_PXL_BiSVMM_Create(
	LPXLOPER XL_modelNames,
	LPXLOPER XL_modelIds,
	LPXLOPER XL_corrMatrixId)
{
	ADD_LOG("Local_PXL_BiSVMM_Create");
	bool PersistentInXL = false;

	return Local_BiSVMMCommon(XL_modelNames, XL_modelIds, XL_corrMatrixId, PersistentInXL);
}


///----------------------------------------------
///----------------------------------------------
///             Local_Create1IRFXModelCommon
///		function to create a two factors interest rates and fx model
/// Inputs :
///     names
///		modelsIds
///		correlation matrix
///----------------------------------------------
///----------------------------------------------

LPXLOPER Local_NP1IRNFXModel_CreateCommon(
	LPXLOPER XL_ModelNames,
	LPXLOPER XL_ModelsIdsVec,
	LPXLOPER XL_CorrelationId,
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

		vector<CCString> modelNames;
		XL_readStrVector(XL_ModelNames,modelNames," ARM_ERR: Model names : array of string expected",DOUBLE_TYPE,C_result);
		vector< string > C_modelNames(modelNames.size());
		for(size_t i=0;i<modelNames.size();++i)
			C_modelNames[i]=CCSTringToSTLString(modelNames[i]);

		vector<CCString> C_ModelsIdsVec;
		XL_readStrVector(XL_ModelsIdsVec,C_ModelsIdsVec," ARM_ERR: Model: array of object expected",DOUBLE_TYPE,C_result);
		size_t size = C_ModelsIdsVec.size();
		vector<long> C_ModelsIdsVecLong(size);
		for(i = 0; i < size; ++i )    
			C_ModelsIdsVecLong[i] = LocalGetNumObjectId(C_ModelsIdsVec[i]); 

		CCString C_CorrelationStrId;
		XL_readStrCell( XL_CorrelationId, C_CorrelationStrId,	" ARM_ERR: Correlation Id: Object expected",		C_result);
		long C_CorrelationId = LocalGetNumObjectId(C_CorrelationStrId);

		/// a function with a context
		exportFunc3Args< vector< string >, vector< long >, long>  ourFunc(C_modelNames, C_ModelsIdsVecLong, C_CorrelationId, ARMLOCAL_NP1IRNFXModel_Create );

		/// call the general function
		fillXL_Result_withName( ourFunc, C_result, XL_result, PersistentInXL );
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_Create2IRFXModelCommon" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


///////////////////////////////////
/// version for non persistentInXl
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_NP1IRNFXModel_Create(
	LPXLOPER XL_ModelNames,
	LPXLOPER XL_ModelsIdsVec,
	LPXLOPER XL_CorrelationId)
{
	ADD_LOG("Local_NP1IRNFXModel_Create");
	bool PersistentInXL = true;
	return Local_NP1IRNFXModel_CreateCommon(
		XL_ModelNames,
		XL_ModelsIdsVec,
		XL_CorrelationId,
		PersistentInXL);
}


///////////////////////////////////
/// version for non persistentInXl
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_PXL_NP1IRNFXModel_Create(
	LPXLOPER XL_ModelNames,
	LPXLOPER XL_ModelsIdsVec,
	LPXLOPER XL_CorrelationId)
{
	ADD_LOG("Local_PXL_NP1IRNFXModel_Create");
	bool PersistentInXL = false;
	return Local_NP1IRNFXModel_CreateCommon(
		XL_ModelNames,
		XL_ModelsIdsVec,
		XL_CorrelationId,
		PersistentInXL);
}


LPXLOPER Local_NP1IRNFX_CalibrateFunctional_Common(
	LPXLOPER XL_NP1IRNFXId,
	LPXLOPER XL_ResetDates,
	LPXLOPER XL_DensitiesId,
	LPXLOPER XL_GridSize,
	LPXLOPER XL_StdDevNb,
	LPXLOPER XL_Rescaling,
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

		long C_NP1IRNFXId;
		XL_GETOBJID( XL_NP1IRNFXId, C_NP1IRNFXId, " ARM_ERR: NP1IRNFXId: Object expected",C_result );

		VECTOR<double> C_ResetDates;
		XL_readNumVector(XL_ResetDates,C_ResetDates," ARM_ERR: Reset Dates: array of dates expected",C_result);

		VECTOR<CCString> C_DensitiesId;
		long C_NbRows, C_NbCols;
		XL_readStrVectorAndSize (XL_DensitiesId,C_NbRows,C_NbCols,C_DensitiesId," ARM_ERR: Densities: array of object expected",DOUBLE_TYPE,C_result);

		size_t size = C_DensitiesId.size();

		VECTOR<long> C_DensitiesIdLong (size);
		size_t i;
		for(i = 0; i < size; ++i )
			C_DensitiesIdLong[i] = LocalGetNumObjectId(C_DensitiesId[i]);
		
		double C_GridSize;
		double C_GridSizeDefault=501;
		XL_readNumCellWD(XL_GridSize, C_GridSize, C_GridSizeDefault, " ARM_ERR: GridSize : numeric expected", C_result);

		double C_StdDevNb;
		double C_StdDevNbDefault=6;
		XL_readNumCellWD(XL_StdDevNb, C_StdDevNb, C_StdDevNbDefault, " ARM_ERR: StdDevNb : numeric expected", C_result);

		CCString C_Rescaling;
		XL_readStrCellWD( XL_Rescaling,C_Rescaling,"N"," ARM_ERR: Rescalling: string expected",C_result);
		int C_RescalingBool;
		C_Rescaling.toUpper();
		if (C_Rescaling == "Y" )
			C_RescalingBool = true;
		else if (C_Rescaling == "N")
			C_RescalingBool = false;
		else
			throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"String \"Y\" or \"N\" Expected for Rescalling");

		/// a function with a context
		exportFunc8Args < long, VECTOR<double>, VECTOR<long>, long, long, double, double, int >  ourFunc (
				C_NP1IRNFXId,
				C_ResetDates,
				C_DensitiesIdLong,
				C_NbRows,
				C_NbCols,
				C_GridSize,
				C_StdDevNb,
				C_RescalingBool,
				ARMLOCAL_NP1IRNFX_CalibrateFunctional );

		/// call the general function
		fillXL_Result_withName(ourFunc, C_result, XL_result, PersistentInXL );
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_Local_Model_Calibrate_Common" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


///////////////////////////////////
/// version for persistentInXl
///////////////////////////////////
LPXLOPER Local_NP1IRNFX_CalibrateFunctional(
	LPXLOPER XL_NP1IRNFXId,
	LPXLOPER XL_ResetDates,
	LPXLOPER XL_DensitiesId,
	LPXLOPER XL_GridSize,
	LPXLOPER XL_StdDevNb,
	LPXLOPER XL_Rescaling)
{
	bool PersistentInXL = true;
	return Local_NP1IRNFX_CalibrateFunctional_Common(
        XL_NP1IRNFXId,
		XL_ResetDates,
		XL_DensitiesId,
		XL_GridSize,
		XL_StdDevNb,
		XL_Rescaling,
		PersistentInXL );
}

///////////////////////////////////
/// version for non persistentInXl
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_PXL_NP1IRNFX_CalibrateFunctional(
	LPXLOPER XL_NP1IRNFXId,
	LPXLOPER XL_ResetDates,
	LPXLOPER XL_DensitiesId,
	LPXLOPER XL_GridSize,
	LPXLOPER XL_StdDevNb,
	LPXLOPER XL_Rescaling)
{
	ADD_LOG("Local_PXL_NP1IRNFX_CalibrateFunctional");
	bool PersistentInXL = false;
	return Local_NP1IRNFX_CalibrateFunctional_Common(
        XL_NP1IRNFXId,
		XL_ResetDates,
		XL_DensitiesId,
		XL_GridSize,
		XL_StdDevNb,
		XL_Rescaling,
		PersistentInXL );
}


LPXLOPER Local_2IRFXSV_CalibrateFunctional_Common(
	LPXLOPER XL_2IRFXSVId,
	LPXLOPER XL_ResetDates,
	LPXLOPER XL_DensitiesId,
	LPXLOPER XL_GridSize,
	LPXLOPER XL_StdDevNb,
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

		long C_2IRFXSVId;
		XL_GETOBJID( XL_2IRFXSVId, C_2IRFXSVId, " ARM_ERR: 2IRFXSVId: Object expected",C_result );

		VECTOR<double> C_ResetDates;
		XL_readNumVector(XL_ResetDates,C_ResetDates," ARM_ERR: Reset Dates: array of dates expected",C_result);

		VECTOR<CCString> C_DensitiesId;
		long C_NbRows, C_NbCols;
		XL_readStrVectorAndSize (XL_DensitiesId,C_NbRows,C_NbCols,C_DensitiesId," ARM_ERR: Densities: array of object expected",DOUBLE_TYPE,C_result);

		size_t size = C_DensitiesId.size();

		VECTOR<long> C_DensitiesIdLong (size);
		size_t i;
		for(i = 0; i < size; ++i )
			C_DensitiesIdLong[i] = LocalGetNumObjectId(C_DensitiesId[i]);
		
		double C_GridSize;
		double C_GridSizeDefault=501;
		XL_readNumCellWD(XL_GridSize, C_GridSize, C_GridSizeDefault, " ARM_ERR: GridSize : numeric expected", C_result);

		double C_StdDevNb;
		double C_StdDevNbDefault=6;
		XL_readNumCellWD(XL_StdDevNb, C_StdDevNb, C_StdDevNbDefault, " ARM_ERR: StdDevNb : numeric expected", C_result);

		/// a function with a context
		exportFunc5Args < long, VECTOR<double>, VECTOR<long>, double, double>  ourFunc (
				C_2IRFXSVId,
				C_ResetDates,
				C_DensitiesIdLong,
				C_GridSize,
				C_StdDevNb,
				ARMLOCAL_2IRFXSV_CalibrateFunctional );

		/// call the general function
		fillXL_Result_withName(ourFunc, C_result, XL_result, PersistentInXL );
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_Local_Model_Calibrate_Common" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


///////////////////////////////////
/// version for persistentInXl
///////////////////////////////////
LPXLOPER Local_2IRFXSV_CalibrateFunctional(
	LPXLOPER XL_2IRFXSVId,
	LPXLOPER XL_ResetDates,
	LPXLOPER XL_DensitiesId,
	LPXLOPER XL_GridSize,
	LPXLOPER XL_StdDevNb)
{
	bool PersistentInXL = true;
	return Local_2IRFXSV_CalibrateFunctional_Common(
        XL_2IRFXSVId,
		XL_ResetDates,
		XL_DensitiesId,
		XL_GridSize,
		XL_StdDevNb,
		PersistentInXL );
}

///////////////////////////////////
/// version for non persistentInXl
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_PXL_2IRFXSV_CalibrateFunctional(
	LPXLOPER XL_2IRFXSVId,
	LPXLOPER XL_ResetDates,
	LPXLOPER XL_DensitiesId,
	LPXLOPER XL_GridSize,
	LPXLOPER XL_StdDevNb)
{
	ADD_LOG("Local_PXL_2IRFXSV_CalibrateFunctional");
	bool PersistentInXL = false;
	return Local_2IRFXSV_CalibrateFunctional_Common(
        XL_2IRFXSVId,
		XL_ResetDates,
		XL_DensitiesId,
		XL_GridSize,
		XL_StdDevNb,
		PersistentInXL );
}

LPXLOPER Local_HWxSVMMSpread_Common(
	LPXLOPER XL_modelNames,
	LPXLOPER XL_modelIds,
	LPXLOPER XL_hw2fId,
	LPXLOPER XL_corrIndexEndTimes,
	LPXLOPER XL_constantCrossCorrel,
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

		vector<CCString> modelNames;
		XL_readStrVector(XL_modelNames,modelNames," ARM_ERR: Model names : array of string expected",DOUBLE_TYPE,C_result);
		vector< string > C_modelNames(modelNames.size());
		for(size_t i=0;i<modelNames.size();++i)
			C_modelNames[i]=CCSTringToSTLString(modelNames[i]);

		vector<CCString> C_ModelsIdsVec;
		XL_readStrVector(XL_modelIds,C_ModelsIdsVec," ARM_ERR: Model: array of object expected",DOUBLE_TYPE,C_result);
		size_t size = C_ModelsIdsVec.size();
		vector<long> C_ModelsIdsVecLong(size);
		for(i = 0; i < size; ++i )    
			C_ModelsIdsVecLong[i] = LocalGetNumObjectId(C_ModelsIdsVec[i]); 

		CCString C_hw2fIdStr, defStr = "";
		XL_readStrCellWD( XL_hw2fId, C_hw2fIdStr, defStr,	" ARM_ERR: Correlation Id: Object expected",		C_result);
		long C_hw2fId = LocalGetNumObjectId(C_hw2fIdStr);

		vector<double> C_corrIdxEndTimes, defvec(0);
		XL_readNumVectorWD(XL_corrIndexEndTimes, C_corrIdxEndTimes, defvec, "ARM_ERR : CorrIndexEndTimes : array of numeric expected", C_result);

		double C_CstCorrel, defcorr = 0.;
		XL_readNumCellWD(XL_constantCrossCorrel, C_CstCorrel, defcorr, "ARM_ERR : constant correl : numeric expected", C_result);


		exportFunc5Args< vector<string>, vector<long>, long, vector<double>, double> ourFunc(
			C_modelNames,
			C_ModelsIdsVecLong,
			C_hw2fId,
			C_corrIdxEndTimes,
			C_CstCorrel,
			ARMLOCAL_HWxSVMMSpread_Create);

		/// call the general function
		fillXL_Result( LOCAL_BiSVMM_CLASS, ourFunc, C_result, XL_result, PersistentInXL );
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in ARMLOCAL_HWxSVMMSpread_Create" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}

__declspec(dllexport) LPXLOPER WINAPI Local_HWxSVMMSpread_Create(
	LPXLOPER XL_modelNames,
	LPXLOPER XL_modelIds,
	LPXLOPER XL_hw2fId,
	LPXLOPER XL_corrIndexEndTimes,
	LPXLOPER XL_constantCrossCorrel)
{
	ADD_LOG("Local_HWxSVMMSpread_Create");
	bool PersistentInXL = true;

	return Local_HWxSVMMSpread_Common(
		XL_modelNames,
		XL_modelIds,
		XL_hw2fId,
		XL_corrIndexEndTimes,
		XL_constantCrossCorrel,
		PersistentInXL);
}

__declspec(dllexport) LPXLOPER WINAPI Local_PXL_HWxSVMMSpread_Create(
	LPXLOPER XL_modelNames,
	LPXLOPER XL_modelIds,
	LPXLOPER XL_hw2fId,
	LPXLOPER XL_corrIndexEndTimes,
	LPXLOPER XL_constantCrossCorrel)
{
	ADD_LOG("Local_PXL_HWxSVMMSpread_Create");
	bool PersistentInXL = false;

	return Local_HWxSVMMSpread_Common(
		XL_modelNames,
		XL_modelIds,
		XL_hw2fId,
		XL_corrIndexEndTimes,
		XL_constantCrossCorrel,
		PersistentInXL);
}

LPXLOPER Local_HWSBGMQtoModel_Create_Common(
	LPXLOPER XL_ModelNames,
	LPXLOPER XL_ModelsIdsVec,
	LPXLOPER XL_CorrelationId,
	bool PersistentInXL)
{
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

		vector<CCString> modelNames;
		XL_readStrVector(XL_ModelNames,modelNames," ARM_ERR: Model names : array of string expected",DOUBLE_TYPE,C_result);
		vector< string > C_modelNames(modelNames.size());
		for(size_t i=0;i<modelNames.size();++i)
			C_modelNames[i]=CCSTringToSTLString(modelNames[i]);

		vector<CCString> C_ModelsIdsVec;
		XL_readStrVector(XL_ModelsIdsVec,C_ModelsIdsVec," ARM_ERR: Model: array of object expected",DOUBLE_TYPE,C_result);
		size_t size = C_ModelsIdsVec.size();
		vector<long> C_ModelsIdsVecLong(size);
		for(i = 0; i < size; ++i )    
			C_ModelsIdsVecLong[i] = LocalGetNumObjectId(C_ModelsIdsVec[i]); 

		CCString C_CorrelationStrId;
		XL_readStrCell( XL_CorrelationId, C_CorrelationStrId,	" ARM_ERR: Correlation Id: Object expected",		C_result);
		long C_CorrelationId = LocalGetNumObjectId(C_CorrelationStrId);

		exportFunc3Args< vector<string>, vector<long>, long> ourFunc(
			C_modelNames,
			C_ModelsIdsVecLong,
			C_CorrelationId,
			ARMLOCAL_HWSBGMQtoModel_Create);

		/// call the general function
		fillXL_Result( LOCAL_HWHWQTOMODEL_CLASS, ourFunc, C_result, XL_result, PersistentInXL );
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_HWSBGMQtoModel_Create" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


_declspec(dllexport) LPXLOPER WINAPI Local_HWSBGMQtoModel_Create(
	LPXLOPER XL_modelNames,
	LPXLOPER XL_modelIds,
	LPXLOPER XL_corrMatrixId)
{
	ADD_LOG("Local_HWSBGMQtoModel_Create");
	bool PersistentInXL = true;

	return Local_HWSBGMQtoModel_Create_Common(XL_modelNames, XL_modelIds, XL_corrMatrixId, PersistentInXL);
}

_declspec(dllexport) LPXLOPER WINAPI Local_PXL_HWSBGMQtoModel_Create(
	LPXLOPER XL_modelNames,
	LPXLOPER XL_modelIds,
	LPXLOPER XL_corrMatrixId)
{
	ADD_LOG("Local_PXL_HWSBGMQtoModel_Create");
	bool PersistentInXL = false;

	return Local_HWSBGMQtoModel_Create_Common(XL_modelNames, XL_modelIds, XL_corrMatrixId, PersistentInXL);
}

LPXLOPER Local_HWSVBGMQtoModel_Create_Common(
	LPXLOPER XL_ModelNames,
	LPXLOPER XL_ModelsIdsVec,
	LPXLOPER XL_CorrelationId,
	bool PersistentInXL)
{
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

		vector<CCString> modelNames;
		XL_readStrVector(XL_ModelNames,modelNames," ARM_ERR: Model names : array of string expected",DOUBLE_TYPE,C_result);
		vector< string > C_modelNames(modelNames.size());
		for(size_t i=0;i<modelNames.size();++i)
			C_modelNames[i]=CCSTringToSTLString(modelNames[i]);

		vector<CCString> C_ModelsIdsVec;
		XL_readStrVector(XL_ModelsIdsVec,C_ModelsIdsVec," ARM_ERR: Model: array of object expected",DOUBLE_TYPE,C_result);
		size_t size = C_ModelsIdsVec.size();
		vector<long> C_ModelsIdsVecLong(size);
		for(i = 0; i < size; ++i )    
			C_ModelsIdsVecLong[i] = LocalGetNumObjectId(C_ModelsIdsVec[i]); 

		CCString C_CorrelationStrId;
		XL_readStrCell( XL_CorrelationId, C_CorrelationStrId,	" ARM_ERR: Correlation Id: Object expected",		C_result);
		long C_CorrelationId = LocalGetNumObjectId(C_CorrelationStrId);

		exportFunc3Args< vector<string>, vector<long>, long> ourFunc(
			C_modelNames,
			C_ModelsIdsVecLong,
			C_CorrelationId,
			ARMLOCAL_HWSVBGMQtoModel_Create);

		/// call the general function
		fillXL_Result( LOCAL_HWHWQTOMODEL_CLASS, ourFunc, C_result, XL_result, PersistentInXL );
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_HWSVBGMQtoModel_Create" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


_declspec(dllexport) LPXLOPER WINAPI Local_HWSVBGMQtoModel_Create(
	LPXLOPER XL_modelNames,
	LPXLOPER XL_modelIds,
	LPXLOPER XL_corrMatrixId)
{
	ADD_LOG("Local_HWSVBGMQtoModel_Create");
	bool PersistentInXL = true;

	return Local_HWSVBGMQtoModel_Create_Common(XL_modelNames, XL_modelIds, XL_corrMatrixId, PersistentInXL);
}

_declspec(dllexport) LPXLOPER WINAPI Local_PXL_HWSVBGMQtoModel_Create(
	LPXLOPER XL_modelNames,
	LPXLOPER XL_modelIds,
	LPXLOPER XL_corrMatrixId)
{
	ADD_LOG("Local_PXL_HWSVBGMQtoModel_Create");
	bool PersistentInXL = false;

	return Local_HWSVBGMQtoModel_Create_Common(XL_modelNames, XL_modelIds, XL_corrMatrixId, PersistentInXL);
}


LPXLOPER Local_2IRFXSV_CreateCommon(
	LPXLOPER XL_ModelNames,
	LPXLOPER XL_ModelsIdsVec,
	LPXLOPER XL_CorrelationId,
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

		vector<CCString> modelNames;
		XL_readStrVector(XL_ModelNames,modelNames," ARM_ERR: Model names : array of string expected",DOUBLE_TYPE,C_result);
		vector< string > C_modelNames(modelNames.size());
		for(size_t i=0;i<modelNames.size();++i)
			C_modelNames[i]=CCSTringToSTLString(modelNames[i]);

		vector<CCString> C_ModelsIdsVec;
		XL_readStrVector(XL_ModelsIdsVec,C_ModelsIdsVec," ARM_ERR: Model: array of object expected",DOUBLE_TYPE,C_result);
		size_t size = C_ModelsIdsVec.size();
		vector<long> C_ModelsIdsVecLong(size);
		for(i = 0; i < size; ++i )    
			C_ModelsIdsVecLong[i] = LocalGetNumObjectId(C_ModelsIdsVec[i]); 

		CCString C_CorrelationStrId;
		XL_readStrCell( XL_CorrelationId, C_CorrelationStrId,	" ARM_ERR: Correlation Id: Object expected",		C_result);
		long C_CorrelationId = LocalGetNumObjectId(C_CorrelationStrId);

		/// a function with a context
		exportFunc3Args< vector< string >, vector< long >, long>  ourFunc(C_modelNames, C_ModelsIdsVecLong, C_CorrelationId, ARMLOCAL_2IRFXSV_Create );

		/// call the general function
		fillXL_Result_withName( ourFunc, C_result, XL_result, PersistentInXL );
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_Create2IRFXModelCommon" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


///////////////////////////////////
/// version for non persistentInXl
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_2IRFXSV_Create(
	LPXLOPER XL_ModelNames,
	LPXLOPER XL_ModelsIdsVec,
	LPXLOPER XL_CorrelationId)
{
	ADD_LOG("Local_2IRFXSV_Create");
	bool PersistentInXL = true;
	return Local_2IRFXSV_CreateCommon(
		XL_ModelNames,
		XL_ModelsIdsVec,
		XL_CorrelationId,
		PersistentInXL);
}


///////////////////////////////////
/// version for non persistentInXl
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_PXL_2IRFXSV_Create(
	LPXLOPER XL_ModelNames,
	LPXLOPER XL_ModelsIdsVec,
	LPXLOPER XL_CorrelationId)
{
	ADD_LOG("Local_PXL_2IRFXSV_Create");
	bool PersistentInXL = false;
	return Local_2IRFXSV_CreateCommon(
		XL_ModelNames,
		XL_ModelsIdsVec,
		XL_CorrelationId,
		PersistentInXL);
}
