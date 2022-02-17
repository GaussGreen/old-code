/*
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 * $Log: ARM_xl_gp_inflation_local.cpp,v $
 * Revision 1.1  2004/02/07 15:08:43  ebenhamou
 * Initial version
 *
 */

/*! \file ARM_xl_gp_inflation_local.cpp
 *
 *  \brief file for the inflation addins
 *
 *	\author  E Benhamou
 *	\version 1.0
 *	\date March 2004
 */


#include <ARM\libarm_local\firstToBeIncluded.h>
#include <functional>
#include <libCCxll\CCxll.h>
#include <ARM\libarm_local\ARM_local_gp_inflation.h>
//#include <ARM\libarm_local\ARM_local_gp_inflationCapFloor.h>
#include "ARM_xl_gp_inflation_local.h"
#include "ARM_xl_gp_fctorhelper.h"
#include "ARM_xl_wrapper_local.h"
#include "ARM_xl_trycatch_local.h"
#include <util\fromto.h>

#include "glob\paramview.h"
#include "ARM_gp_local_interglob.h"
#include "ARM_xl_gp_common.h"
#include <GP_Infra\gpinfra\gramfunctorarg.h>
#include <memory>
#include <string>
#include <sstream>
#include "util\tech_macro.h"

using namespace std;

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
class infCurvFctor : public ARMResultLong2LongFunc
{
public:
	infCurvFctor(	
		double					asOfDate,
		const CCString&			indexName,
		double					CPIIndexValue,
		double					CPIIndexDate,
		const VECTOR<CCString>&	maturities,
		const VECTOR<double>&	values,
		long					MonthlyInterpType,
		long					DailyInterpType,
		long					DCFMonthly,
		long					DCFDaily,
		long					ExtrapolType,
 		long					ResetManagerId,
		long					SeasonManagerId)
	:
		C_asOfDate(asOfDate ),
		C_indexName(indexName),
		C_CPIIndexValue(CPIIndexValue),
		C_CPIIndexDate(CPIIndexDate),
		C_maturities(maturities),
		C_values(values),
		C_MonthlyInterpType(MonthlyInterpType),
		C_DailyInterpType(DailyInterpType),
		C_DCFMonthly( DCFMonthly ),
		C_DCFDaily( DCFDaily ),
		C_ExtrapolType(ExtrapolType),
		C_ResetManagerId( ResetManagerId ),
		C_SeasonManagerId( SeasonManagerId )
	{};
	
	long operator()( ARM_result& result, long objId ) 
	{
		return ARMLOCAL_InfCurv_Create(
			C_asOfDate,
			C_indexName,			
			C_CPIIndexValue,
			C_CPIIndexDate,
			C_maturities,
			C_values,
			C_MonthlyInterpType,
			C_DailyInterpType,
			C_DCFMonthly,
			C_DCFDaily,
			C_ExtrapolType,
			C_ResetManagerId,
			C_SeasonManagerId,
			result,
			objId );			
	}

	
private:
	double				C_asOfDate;
	CCString			C_indexName;
	double				C_CPIIndexValue;
	double				C_CPIIndexDate;
	VECTOR<CCString>	C_maturities;
	VECTOR<double>		C_values;
	long				C_MonthlyInterpType;
	long				C_DailyInterpType;
	long				C_DCFMonthly;
	long				C_DCFDaily;
	long				C_ExtrapolType;
	long				C_ResetManagerId;
	long				C_SeasonManagerId;
};


/////////////////////////////////////////////////////////////
/// central function that does the creation of the XL function
/////////////////////////////////////////////////////////////
LPXLOPER Local_InfCurvCommon(
	LPXLOPER XL_asOfDate,
	LPXLOPER XL_indexName,
	LPXLOPER XL_CPIIndexValue,
	LPXLOPER XL_CPIIndexDate,
	LPXLOPER XL_maturities,
	LPXLOPER XL_values,
	LPXLOPER XL_MonthlyInterpType,
	LPXLOPER XL_DailyInterpType,
	LPXLOPER XL_DCFMonthly,
	LPXLOPER XL_DCFDaily,
	LPXLOPER XL_ExtrapolType,
	LPXLOPER XL_ResetManager,
	LPXLOPER XL_SeasonManager,
	bool PersistentInXL )
{
	
	/// this is defined first because it is used in XL macros
	static XLOPER XL_result;
	/// Get the variables from the XLOper variables
	ARM_result C_result;
	
	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
		ARM_NOCALCIFWIZ();
		
		double C_asOfDate;
		CCString C_indexName;
		double C_CPIIndexValue;
		double C_CPIIndexDate;
		VECTOR<CCString> C_maturities;
		VECTOR<double> C_values;
		
		/// this is used by macros 
		/// and therefore this has to be defined
		static int error;
		static char* reason = "";
		
		XL_readNumCell(	 XL_asOfDate,		C_asOfDate,		" ARM_ERR: asOfDate: date expected",				C_result);
		XL_readStrCell(	 XL_indexName,		C_indexName,	" ARM_ERR: index name string expected",				C_result);
		XL_readNumCell(	 XL_CPIIndexValue,	C_CPIIndexValue," ARM_ERR: index value: double expected",			C_result);
		XL_readNumCell(	 XL_CPIIndexDate,	C_CPIIndexDate,	" ARM_ERR: index Date: date expected",				C_result);
		XL_readStrVector(XL_maturities,		C_maturities,	" ARM_ERR: maturities and rates: array of numeric expected", DOUBLE_TYPE,C_result);
		XL_readNumVector(XL_values,			C_values,		" ARM_ERR: index values: array of numeric expected",C_result);
		
		/// reading of all the flags as string
		/// if missing defaults are the following
		CCString C_MonthlyInterpTypeStr;
		CCString C_DailyInterpTypeStr;
		CCString C_DCFMonthlyStr;
		CCString C_DCFDailyStr;
		CCString C_ExtrapolTypeStr;
		
		long C_MonthlyInterpType;
		long C_DailyInterpType;
		long C_DCFMonthly;
		long C_DCFDaily;
		long C_ExtrapolType;
		
		char defaultValue[3] = "-1";
		
		XL_readStrCellWD(XL_MonthlyInterpType,C_MonthlyInterpTypeStr, defaultValue, " ARM_ERR: monthly interpolation type: string expected",C_result);
		XL_readStrCellWD(XL_DailyInterpType,C_DailyInterpTypeStr, defaultValue ,	" ARM_ERR: daily interpolation type: string expected",	C_result);
		XL_readStrCellWD(XL_DCFMonthly,		C_DCFMonthlyStr,	defaultValue ,		" ARM_ERR: DCF Monthly type: string expected",	C_result);
		XL_readStrCellWD(XL_DCFDaily,		C_DCFDailyStr,		defaultValue ,		" ARM_ERR: DCD Daily type: string expected",	C_result);
		XL_readStrCellWD(XL_ExtrapolType,	C_ExtrapolTypeStr,	defaultValue ,		" ARM_ERR: extrapolaiton type: string expected",C_result);
		
		/// conversion to long
		if( ( C_MonthlyInterpType	= ARM_ConvCPIMonthlyInterpMethod( C_MonthlyInterpTypeStr, C_result)) == ARM_DEFAULT_ERR
			||	(C_DailyInterpType	= ARM_ConvCPIDailyInterpMethod( C_DailyInterpTypeStr, C_result)) == ARM_DEFAULT_ERR
			||	(C_ExtrapolType		= ARM_ConvCPIExtrapolMethod( C_ExtrapolTypeStr, C_result)) == ARM_DEFAULT_ERR
			)
		{
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}
		C_DCFMonthly		= ARM_ConvDayCount( C_DCFMonthlyStr);
		C_DCFDaily			= ARM_ConvDayCount( C_DCFDailyStr);

		long C_ResetManagerId;
		XL_GETOBJIDWD( XL_ResetManager, C_ResetManagerId, GETDEFAULTVALUESTR, " ARM_ERR: Reset Manager: Obj expected", C_result );

		long C_SeasonManagerId;

		XL_GETOBJIDWD( XL_SeasonManager, C_SeasonManagerId, GETDEFAULTVALUESTR, " ARM_ERR: Season Manager: Obj expected", C_result );


		/// use the concept of Functor to transfer the knowledge of
		/// a function with a context
		infCurvFctor ourFunc(
			C_asOfDate,
			C_indexName,
			C_CPIIndexValue,
			C_CPIIndexDate,
			C_maturities,
			C_values,
			C_MonthlyInterpType,
			C_DailyInterpType,
			C_DCFMonthly,
			C_DCFDaily,
			C_ExtrapolType,
			C_ResetManagerId,
			C_SeasonManagerId);
		
		/// call the general function
		fillXL_Result( LOCAL_INFCURV_CLASS, ourFunc, C_result, XL_result, PersistentInXL );
		
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_InfCurvCommon" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}
							 
							 
///////////////////////////////////
/// version that takes into account 
/// previous creation of object
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_CreateInfCurv(
	LPXLOPER XL_asOfDate,
	LPXLOPER XL_indexName,
	LPXLOPER XL_CPIIndexValue,
	LPXLOPER XL_CPIIndexDate,
	LPXLOPER XL_maturities,
	LPXLOPER XL_values,
	LPXLOPER XL_MonthlyInterpType,
	LPXLOPER XL_DailyInterpType,
	LPXLOPER XL_DCFMonthly,
	LPXLOPER XL_DCFDaily,
	LPXLOPER XL_ExtrapolType,
	LPXLOPER XL_ResetManager,
	LPXLOPER XL_SeasonManager)
{
	ADD_LOG("Local_CreateInfCurv");
	bool PersistentInXL = true;
	return Local_InfCurvCommon(
		XL_asOfDate,
		XL_indexName,
		XL_CPIIndexValue,
		XL_CPIIndexDate,
		XL_maturities,
		XL_values,
		XL_MonthlyInterpType,
		XL_DailyInterpType,
		XL_DCFMonthly,
		XL_DCFDaily,
		XL_ExtrapolType,
		XL_ResetManager,
		XL_SeasonManager,
		PersistentInXL );
}


///////////////////////////////////
/// version for non persistentInXl
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_PXL_CreateInfCurv(
	LPXLOPER XL_asOfDate,
	LPXLOPER XL_indexName,
	LPXLOPER XL_CPIIndexValue,
	LPXLOPER XL_CPIIndexDate,
	LPXLOPER XL_maturities,
	LPXLOPER XL_values,
	LPXLOPER XL_MonthlyInterpType,
	LPXLOPER XL_DailyInterpType,
	LPXLOPER XL_DCFMonthly,
	LPXLOPER XL_DCFDaily,
	LPXLOPER XL_ExtrapolType,
	LPXLOPER XL_ResetManager,
	LPXLOPER XL_SeasonManager)
{
	ADD_LOG("Local_PXL_CreateInfCurv");
	bool PersistentInXL = false;
	return Local_InfCurvCommon(
		XL_asOfDate,
		XL_indexName,
		XL_CPIIndexValue,
		XL_CPIIndexDate,
		XL_maturities,
		XL_values,
		XL_MonthlyInterpType,
		XL_DailyInterpType,
		XL_DCFMonthly,
		XL_DCFDaily,
		XL_ExtrapolType,
		XL_ResetManager,
		XL_SeasonManager,
		PersistentInXL );
}



////////////////////////////////////////////
/// addin to interpolate an inflation curve
////////////////////////////////////////////
LPXLOPER WINAPI Local_InterpInfCurv(	
	LPXLOPER XL_curv,
	LPXLOPER XL_CPIDate,
	LPXLOPER XL_DCFLag,
	LPXLOPER XL_DailyInterpType,
	LPXLOPER XL_CPILag,
	LPXLOPER XL_weight,
	instrument what  )
{
	// this is defined first because it is used in XL macros
	static XLOPER XL_result;
	ARM_result C_result;
	
	// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	/// to avoid computation if called by the wizard
	ARM_NOCALCIFWIZ();
	
	/// C variable to get the variables from LPXLOPER 
	/// Variables
	CCString C_curv;
	double	 C_CPIDate;
	CCString C_DCFLag;
	CCString C_DailyInterpTypeStr;
	CCString C_CPILag;
	double	 C_weight;
	
	// error
	static int error;
	static char* reason = "";
	
	XL_readStrCell(XL_curv, C_curv,				" ARM_ERR: curve id: object expected",				C_result);
	XL_readNumCell(XL_CPIDate,C_CPIDate,		" ARM_ERR: maturity: date expected",				C_result);
	
	/// reading of all the flags as string
	/// if missing defaults are the following
	char defaultValue[3] = "-1";
	
	XL_readStrCellWD(XL_DCFLag, C_DCFLag, defaultValue,	" ARM_ERR: DCF Flag: string expected",				C_result);
	XL_readStrCellWD(XL_DailyInterpType, C_DailyInterpTypeStr, defaultValue, 
														" ARM_ERR: Daily Interpolation Flag: string expected",				C_result);
	XL_readStrCellWD(XL_CPILag, C_CPILag, defaultValue,	" ARM_ERR: CPI Lag: string expected",				C_result);
	
	/// if weight is not provided, initialised to zero
	double defaultWeight = 0;
	XL_readNumCellWD(XL_weight,C_weight, defaultWeight,	" ARM_ERR: weight: weight expected",		C_result);
	
	/// checks that the weight is meaningful
	if( 0 > C_weight || 1 < C_weight )
	{
		C_result.setMsg ("ARM_ERR: weight should be between 0 and 1");
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}
	
	/// convert from string to long interpType 
	// using general conversion fctions
	long C_DailyInterpType;
	
	if( (C_DailyInterpType = ARM_ConvCPIDailyInterpMethod( C_DailyInterpTypeStr, C_result)) == ARM_DEFAULT_ERR )
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}
	
	/// call the function with nomore LPXLOPER objects
	// this fction store the result into C_result
	long retCode;
	if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 1)
		retCode = ARMLOCAL_InfCurv_Interp( 
			LocalGetNumObjectId (C_curv), 
			C_CPIDate,
			C_result,
			C_DCFLag,
			C_DailyInterpType,
			C_CPILag,
			C_weight,
			what );
	else
		retCode = ARM_KO;
	
	/// feed the LPXLOPER object result 
	if (retCode == ARM_OK)
	{
		FreeCurCellErr ();
		XL_result.xltype = xltypeNum;
		XL_result.val.num = C_result.getDouble();
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_InterpInfCurv" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


////////////////////////////
/// interpolation of the CPI
////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_InfCurv_CPIInterp(	
	LPXLOPER XL_curv,
	LPXLOPER XL_CPIDate,
	LPXLOPER XL_DCFLag,
	LPXLOPER XL_DailyInterpType,
	LPXLOPER XL_CPILag,
	LPXLOPER XL_weight )
{
	ADD_LOG("Local_InfCurv_CPIInterp");
	instrument inst = CPI;
	return Local_InterpInfCurv( XL_curv,		
		XL_CPIDate,	XL_DCFLag, XL_DailyInterpType,		
		XL_CPILag, XL_weight, inst );
}



////////////////////////////////
/// interpolation of the ZC Rate
////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_InfCurv_ZCRateInterp(	
	LPXLOPER XL_curv,
	LPXLOPER XL_CPIDate,
	LPXLOPER XL_DCFLag,
	LPXLOPER XL_DailyInterpType,
	LPXLOPER XL_CPILag,
	LPXLOPER XL_weight )
{
	ADD_LOG("Local_InfCurv_ZCRateInterp");
	instrument inst = ZCRate;
	return Local_InterpInfCurv( XL_curv,		
		XL_CPIDate,	XL_DCFLag, XL_DailyInterpType,		
		XL_CPILag, XL_weight, inst );
}


///////////////////////////////////////////////
/// very rapid definition of ourlonglongFunctor
///////////////////////////////////////////////
class infIdxFctor : public ARMResultLong2LongFunc
{
public:
	infIdxFctor(	
		const CCString&	indexName,
		const CCString&	resetLag,
		const CCString&	DCFLag,
 		long ccyId, long infCurveId )
	:
		C_indexName( indexName ),
		C_resetLag( resetLag ),
		C_DCFLag( DCFLag ),
 		C_ccyId( ccyId ),
		C_infCurveId( infCurveId)
	{};
	
	long operator()( ARM_result& result, long objId ) 
	{
		return ARMLOCAL_InfIdx_Create(
			/// context
			C_indexName,			
			C_resetLag,
			C_DCFLag,
			C_ccyId,C_infCurveId,
			/// arguments
			result,								
			objId);
	}

private:
	CCString C_indexName;
	CCString C_resetLag;
	CCString C_DCFLag;
 	long C_ccyId, C_infCurveId;
};



/////////////////////////////////////////////////////////////
/// central function that does the creation of the XL function
/////////////////////////////////////////////////////////////
LPXLOPER Local_InfIdxCreateCommon(
	LPXLOPER XL_indexName,
	LPXLOPER XL_resetLag,
	LPXLOPER XL_DCFLag,
	LPXLOPER XL_ccyId,
	LPXLOPER XL_infCurveId,
	bool PersistentInXL )
{
	// this is defined first because it is used in XL macros
	static XLOPER XL_result;
	ARM_result C_result;
	// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
		/// to avoid computation if called by the wizard
	ARM_NOCALCIFWIZ();
	
	/// Get the variables from the XLOper variables
	CCString C_indexName;
	CCString C_resetLag;
	CCString C_DCFLag;
	long C_ccyId, C_infCurveId ;

	/// this is used by macros 
	/// and therefore this has to be defined
	static int error;
	static char* reason = "";

	XL_readStrCell	(XL_indexName,	C_indexName,						" ARM_ERR: index name string expected",			C_result );
	XL_readStrCellWD(XL_resetLag,	C_resetLag,		GETDEFAULTVALUESTR,	" ARM_ERR: DCF Flag: string expected",			C_result );
	XL_readStrCellWD(XL_DCFLag,		C_DCFLag,		GETDEFAULTVALUESTR,	" ARM_ERR: DCF Flag: string expected",			C_result );
	XL_GETOBJIDWD	(XL_ccyId,		C_ccyId,		GETDEFAULTVALUESTR,	" ARM_ERR: Ccy Flag: string expected",			C_result );
	XL_GETOBJIDWD	(XL_infCurveId,	C_infCurveId,	GETDEFAULTVALUESTR,	" ARM_ERR: InfCurve Flag: string expected",		C_result );

	/// use the concept of Functor to transfer the knowledge of
	/// a function with a context
	infIdxFctor ourFunc(
		C_indexName,
		C_resetLag,
		C_DCFLag,
		C_ccyId,
		C_infCurveId);

	/// call the general function
	fillXL_Result( LOCAL_INFIDX_CLASS, ourFunc, C_result, XL_result, PersistentInXL );
	
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_InfIdxCreateCommon" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


//////////////////////////////////
/// version that takes into account 
/// previous creation of object
//////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_InfIdxCreate(
	LPXLOPER XL_indexName,
	LPXLOPER XL_resetLag,
	LPXLOPER XL_DCFLag,
	LPXLOPER XL_ccyId,
	LPXLOPER XL_infCurveId)
{
	ADD_LOG("Local_InfIdxCreate");
 bool PersistentInXL = true;
 return Local_InfIdxCreateCommon(
	XL_indexName,
	XL_resetLag,
	XL_DCFLag,
 	XL_ccyId,
	XL_infCurveId,
	PersistentInXL );
}

///////////////////////////////////
/// version for non persistentInXl
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_PXL_InfIdxCreate(
	LPXLOPER XL_indexName,
	LPXLOPER XL_resetLag,
	LPXLOPER XL_DCFLag,
	LPXLOPER XL_ccyId,
	LPXLOPER XL_infCurveId)
{
	ADD_LOG("Local_PXL_InfIdxCreate");
 bool PersistentInXL = false;
 return Local_InfIdxCreateCommon(
	XL_indexName,
	XL_resetLag,
	XL_DCFLag,
 	XL_ccyId,
	XL_infCurveId,
	PersistentInXL );
}



/// Functor for the inflation swap leg generic
/// uses datestrip as underlying object
class infLegwDateStripFunc: public ARMResultLong2LongFunc
{
public:
	infLegwDateStripFunc(
		const CCString& infIdxName,
		long rcvOrPay,
		int interpType,
		double multiple,
        double constant,
		long finalNotionalType,
		long numStripDateId,
		long denomStripDateId )
	:
		C_infIdxName( infIdxName ),
		C_rcvOrPay( rcvOrPay ),
		C_interpType( interpType ),
		C_multiple( multiple ),
		C_constant( constant ),
		C_finalNotionalType( finalNotionalType ),
		C_numStripDateId( numStripDateId ),
		C_denomStripDateId( denomStripDateId )
	{};

	long operator()( ARM_result& result, long objId ) 
	{
		return ARMLOCAL_InfLegwDateStrip_Create(
			/// context
			C_infIdxName,
			C_rcvOrPay,
			C_interpType,
			C_multiple,
			C_constant,
			C_finalNotionalType,
			C_numStripDateId,
			C_denomStripDateId,
			/// arguments
			result,								
			objId );
	}
private:
	CCString C_infIdxName;
	long C_rcvOrPay;
	int C_interpType;
	double C_multiple;
    double C_constant;
	long C_finalNotionalType;
	long C_numStripDateId;
	long C_denomStripDateId;
};


////////////////////////////////////////////////////////////
/// central function that does the creation of the XL function
//////////////////////////////////////////////////////////////
LPXLOPER Local_InfLegwDatesStripCreateCommon(
    LPXLOPER XL_infIdx, 
	LPXLOPER XL_rcvOrPay,
	LPXLOPER XL_interpType,
	LPXLOPER XL_multiple,
    LPXLOPER XL_constant,
	LPXLOPER XL_finalNotionalType,
	LPXLOPER XL_numStripDateId,
	LPXLOPER XL_denomStripDateId,
	bool PersistentInXL )
{
	/// this is defined first because it is used in XL macros
	static XLOPER XL_result;
	/// Get the variables from the XLOper variables
	ARM_result C_result;

	// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	/// to avoid computation if called by the wizard
	ARM_NOCALCIFWIZ();
	
	/// in the case of a null object reads an index name
	CCString C_infIdxName;
	long C_rcvOrPay;
	int C_interpType;
	double C_multiple;
    double C_constant;
	long C_finalNotionalType;
	long C_numStripDateId;
	long C_denomStripDateId;

	/// this is used by macros 
	/// and therefore this has to be defined
	static int error;
	static char* reason = "";
	double defaultinterp	= 0;
	double defaultConstant	= 0;
	double defaultMultiple	= 1;

	// read 
	XL_readStrCell(	 XL_infIdx,	C_infIdxName, " ARM_ERR: index name string expected", C_result );
	XL_GETRCVORPAY( XL_rcvOrPay,		C_rcvOrPay,			" ARM_ERR: receive or pay: string expected",			C_result);
	XL_GETCPIDLYINTERPWD( XL_interpType,C_interpType,"CPILINEAR",
														" ARM_ERR: interpType: numeric expected",			C_result);
	XL_readNumCellWD(XL_multiple,		C_multiple,			defaultMultiple," ARM_ERR: multiple: numeric expected",		C_result);
	XL_readNumCellWD(XL_constant,		C_constant,			defaultConstant," ARM_ERR: multiple: numeric expected",		C_result);
	XL_GETNXTYPEWD(	XL_finalNotionalType,C_finalNotionalType,"NONE", " ARM_ERR: notional type: string expected",	C_result);
	XL_GETOBJID(	XL_numStripDateId,	C_numStripDateId,	" ARM_ERR: Numerator Date Strip: Object expected",		C_result);
	XL_GETOBJID(	XL_denomStripDateId,C_denomStripDateId,	" ARM_ERR: Denominator Date Strip: Object expected",	C_result);
	
	/// use the concept of Functor to transfer the knowledge of
	/// a function with a context
	infLegwDateStripFunc ourFunc(
		C_infIdxName,
		C_rcvOrPay,
		C_interpType,
		C_multiple,
		C_constant,
		C_finalNotionalType,
		C_numStripDateId,
		C_denomStripDateId );

	/// call the general function
	fillXL_Result( LOCAL_INFSWAPLEG_CLASS, ourFunc, C_result, XL_result, PersistentInXL );
	
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_InfLegwDatesStripCreateCommon" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


//////////////////////////////////
/// version that takes into account 
/// previous creation of object
//////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_InfLegwDateStripCreate(
    LPXLOPER XL_infIdx, 
	LPXLOPER XL_rcvOrPay,
	LPXLOPER XL_interpType,
	LPXLOPER XL_multiple,
    LPXLOPER XL_constant,
	LPXLOPER XL_finalNotionalType,
	LPXLOPER XL_numStripDateId,
	LPXLOPER XL_denomStripDateId )
{
	ADD_LOG("Local_InfLegwDateStripCreate");
	bool PersistentInXL = true;
	return Local_InfLegwDatesStripCreateCommon(
		XL_infIdx, 
		XL_rcvOrPay,
		XL_interpType,
		XL_multiple,
		XL_constant,
		XL_finalNotionalType,
		XL_numStripDateId,
		XL_denomStripDateId,
		PersistentInXL );
}

///////////////////////////////////
/// version for non persistentInXl
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_PXL_InfLegwDateStripCreate(
    LPXLOPER XL_infIdx,
	LPXLOPER XL_rcvOrPay,
	LPXLOPER XL_interpType,
	LPXLOPER XL_multiple,
    LPXLOPER XL_constant,
	LPXLOPER XL_finalNotionalType,
	LPXLOPER XL_numStripDateId,
	LPXLOPER XL_denomStripDateId )
{
	ADD_LOG("Local_PXL_InfLegwDateStripCreate");
	bool PersistentInXL = false;
	return Local_InfLegwDatesStripCreateCommon(
		XL_infIdx, 
		XL_rcvOrPay,
		XL_interpType,
		XL_multiple,
		XL_constant,
		XL_finalNotionalType,
		XL_numStripDateId,
		XL_denomStripDateId,
		PersistentInXL );
}

/*****************************************************************************************************/

//////////////////////////////////////////////////
/// Version that specifies everything
/// and do not use any underlying datestrip object
//////////////////////////////////////////////////

class infLegAllInputsFunc: public ARMResultLong2LongFunc
{
public:
	infLegAllInputsFunc(
		double startDate, 
		double endDate,
		const CCString& infIdxName,
		int swapType,
		int rcvOrPay,
	    int interpType,
		double multiple,
		double constant,
		int resetFreq,
		int dayCount,
		const CCString& resetCalendar,
		int fwdRule,
		int intRule,
		int stubRule,
		int resetNumGap,
	    int resetDenomGap,
		int payFreq,
		int payGap,
		const CCString& payCalendar,
		int adjFirstDate,
		long finalNotionalType = K_NX_NONE,
		double firstReset = GETDEFAULTVALUE,
		double Comultiple = 1.0)
	:
		C_startDate( startDate ),  
		C_endDate( endDate ),
		C_infIdxName( infIdxName ),
		C_swapType( swapType ),
		C_rcvOrPay( rcvOrPay ),
		C_interpType( interpType ),
		C_multiple( multiple ),
		C_Comultiple( Comultiple ),
		C_constant( constant ),
		C_resetFreq( resetFreq ),
		C_dayCount( dayCount ),
		C_resetCalendar( resetCalendar ),
		C_fwdRule( fwdRule ),
		C_intRule( intRule ),
		C_stubRule( stubRule ),
		C_resetNumGap( resetNumGap ),
		C_resetDenomGap( resetDenomGap ),
		C_payFreq( payFreq ),
		C_payGap( payGap ),
		C_payCalendar( payCalendar ),
		C_adjFirstDate( adjFirstDate ),
		C_finalNotionalType ( finalNotionalType  ),
		C_firstReset( firstReset )
	{};


infLegAllInputsFunc(
		double startDate, 
		double endDate,
		const CCString& infIdxName,
		const CCString& irIdxName,
		int swapType,
		int rcvOrPay,
	    int interpType,
		double multiple,
		double constant,
		int resetFreq,
		int dayCount,
		const CCString& resetCalendar,
		int fwdRule,
		int intRule,
		int stubRule,
		int resetNumGap,
	    int resetDenomGap,
		int payFreq,
		int payGap,
		const CCString& payCalendar,
		int adjFirstDate,
		long finalNotionalType = K_NX_NONE,
		double firstReset = GETDEFAULTVALUE,
		double Comultiple = 1.0)
	:
		C_startDate( startDate ),  
		C_endDate( endDate ),
		C_infIdxName( infIdxName ),
		C_irIdxName( irIdxName ),
		C_swapType( swapType ),
		C_rcvOrPay( rcvOrPay ),
		C_interpType( interpType ),
		C_multiple( multiple ),
		C_Comultiple( Comultiple ),
		C_constant( constant ),
		C_resetFreq( resetFreq ),
		C_dayCount( dayCount ),
		C_resetCalendar( resetCalendar ),
		C_fwdRule( fwdRule ),
		C_intRule( intRule ),
		C_stubRule( stubRule ),
		C_resetNumGap( resetNumGap ),
		C_resetDenomGap( resetDenomGap ),
		C_payFreq( payFreq ),
		C_payGap( payGap ),
		C_payCalendar( payCalendar ),
		C_adjFirstDate( adjFirstDate ),
		C_finalNotionalType ( finalNotionalType  ),
		C_firstReset( firstReset )
	{};





	long operator()( ARM_result& result, long objId ) 
	{
		return ARMLOCAL_InfLegAllInputs_Create(
			/// context
			C_startDate, 
			C_endDate,
			C_infIdxName,
			C_swapType,
			C_rcvOrPay,
			C_interpType,
			C_multiple,
			C_Comultiple,
			C_constant,
			C_resetFreq,
			C_dayCount,
			C_resetCalendar,
			C_fwdRule,
			C_intRule,
			C_stubRule,
			C_resetNumGap,
			C_resetDenomGap,
			C_payFreq,
			C_payGap,
			C_payCalendar,
			C_adjFirstDate,
			C_finalNotionalType,
			C_firstReset,
			/// arguments
			result,								
			objId );
	}
private:
	double C_startDate; 
	double C_endDate;
	CCString C_infIdxName;
	CCString C_irIdxName;
	int C_swapType;
	int C_rcvOrPay;
	int C_interpType;
	double C_multiple;
	double C_Comultiple;
	double C_constant;
	int C_resetFreq;
	int C_dayCount;
	CCString C_resetCalendar;
	int C_fwdRule;
	int C_intRule;
	int C_stubRule;
	int C_resetNumGap;
	int C_resetDenomGap;
	int C_payFreq;
	int C_payGap;
	CCString C_payCalendar;
	int C_adjFirstDate;
	long C_finalNotionalType;
	double C_firstReset;
};



/////////////////////////////////////////////////////////////
/// central function that does the creation of the XL function
/////////////////////////////////////////////////////////////
LPXLOPER Local_InfLegAllInputsCreateCommon(
    LPXLOPER XL_startDate, 
	LPXLOPER XL_endDate,
	LPXLOPER XL_infIdx, 
    LPXLOPER XL_rcvOrPay,
    LPXLOPER XL_interpType,
	LPXLOPER XL_multiple,
	LPXLOPER XL_constant,
    LPXLOPER XL_resetFreq,
	LPXLOPER XL_dayCount,
    LPXLOPER XL_resetCalendar,
	LPXLOPER XL_fwdRule,
	LPXLOPER XL_intRule,
	LPXLOPER XL_stubRule,
	LPXLOPER XL_resetNumGap,
	LPXLOPER XL_resetDenomGap,
	LPXLOPER XL_payFreq,
	LPXLOPER XL_payGap,
	LPXLOPER XL_payCalendar,
	LPXLOPER XL_adjFirstDate,
	LPXLOPER XL_Comultiple,
	long finalNotionalType,
	double firstReset,
	int swapType,
	bool PersistentInXL )
{
	// this is defined first because it is used in XL macros
	static XLOPER XL_result;
	ARM_result C_result;

	
	// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	/// to avoid computation if called by the wizard
	ARM_NOCALCIFWIZ();
	
	/// Get the variables from the XLOper variables
	double C_startDate; 
	double C_endDate;
	CCString C_infIdxName;
	
	int C_rcvOrPay;
	int C_interpType;
	double C_multiple;
	double C_Comultiple;
	double C_constant;
	int C_stubRule;
	int C_dayCount;
	int C_resetFreq;
	int C_fwdRule;
	int C_intRule;
	double C_resetNumGap;
	double C_resetDenomGap;
	CCString C_resetCalendar;
	int C_payFreq;
	double C_payGap;
	CCString C_payCalendar;
	double C_adjFirstDate;

	/// this is used by macros 
	/// and therefore this has to be defined
	static int error;
	static char* reason = "";
	/// we adjust the first date 
	double adjFirstDateDefaultValue = K_UNADJUSTED;
	double gapDefault= 0;
	double defaultMultiple = 1;
	double defaultCoMultiple = 1;
	double defaultConstant = 0;

	// read 
	XL_readNumCell(	 XL_startDate,		C_startDate,	" ARM_ERR: start Date: date expected	",			C_result);
	XL_readNumCell(	 XL_endDate,		C_endDate,		" ARM_ERR: end Date: date expected",				C_result);

	XL_readStrCell(	 XL_infIdx,	C_infIdxName,					" ARM_ERR: index name string expected", C_result );
	XL_GETRCVORPAY( XL_rcvOrPay,		C_rcvOrPay,				" ARM_ERR: receive or pay: string expected",		C_result);
	XL_GETCPIDLYINTERPWD( XL_interpType,C_interpType,"CPILINEAR"," ARM_ERR: interpType: numeric expected",			C_result);
	XL_readNumCellWD( XL_multiple, C_multiple, defaultMultiple, " ARM_ERR: multiple: numeric expected",				C_result);
	XL_readNumCellWD( XL_Comultiple, C_Comultiple, defaultCoMultiple, " ARM_ERR: Comultiple: numeric expected",				C_result);
	XL_readNumCellWD( XL_constant, C_constant, defaultConstant, " ARM_ERR: constant: numeric expected",				C_result);
	XL_GETFREQUENCYWD(XL_resetFreq,C_resetFreq, "-1",			" ARM_ERR: reset freq: string expected",			C_result );
	XL_GETDAYCOUNTWD( XL_dayCount, C_dayCount, "30/360",		" ARM_ERR: dayCount : string expected",				C_result );
	XL_readStrCellWD(XL_resetCalendar,C_resetCalendar, "INF",	" ARM_ERR: reset Calendar: string expected",	C_result);
	XL_GETFWDRULEWD( XL_fwdRule, C_fwdRule, "MF",				" ARM_ERR: fwdRule : string expected",				C_result );
	XL_GETINTRULEWD( XL_intRule, C_intRule, "UNADJ",			" ARM_ERR: fwdRule : string expected",				C_result );
	XL_GETCONVRULEWD( XL_stubRule, C_stubRule, "SS",			" ARM_ERR: stub Rule : string expected",			C_result );
	XL_readNumCellWD(XL_resetNumGap,C_resetNumGap,gapDefault,	" ARM_ERR: numerator reset gap: numeric expected",			C_result);
	
	/// in the case of not specified input as reset DenomGap.. take the same input as numResetGap
	XL_readNumCellWD(XL_resetDenomGap,C_resetDenomGap,C_resetNumGap, " ARM_ERR: denominator reset gap: numeric expected",			C_result);
	XL_GETFREQUENCYWD( XL_payFreq, C_payFreq, "-1",				" ARM_ERR: pay frequency: string expected",			C_result );
	XL_readNumCellWD(XL_payGap, C_payGap, gapDefault,			" ARM_ERR: pay gap: numeric expected",				C_result);
	XL_readStrCellWD(XL_payCalendar, C_payCalendar,	GETDEFAULTVALUESTR," ARM_ERR: pay Calendar: string expected",		C_result);
	XL_readNumCellWD(XL_adjFirstDate,C_adjFirstDate,adjFirstDateDefaultValue, " ARM_ERR: adjust First Date: numeric expected",C_result);
	
	/// use the concept of Functor to transfer the knowledge of
	/// a function with a context
	infLegAllInputsFunc ourFunc(																							
		C_startDate, 
		C_endDate,
		C_infIdxName, 
		swapType,
		C_rcvOrPay,
		C_interpType,
		C_multiple,
		C_constant,
		C_resetFreq,
		C_dayCount,
		C_resetCalendar,
		C_fwdRule,
		C_intRule,
		C_stubRule,
		C_resetNumGap,
		C_resetDenomGap,
		C_payFreq,
		C_payGap,
		C_payCalendar,
		C_adjFirstDate,
		finalNotionalType,
		firstReset,
		C_Comultiple);

	/// call the general functionourFunc
	fillXL_Result( LOCAL_INFSWAPLEG_CLASS, ourFunc, C_result, XL_result, PersistentInXL );
	
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_InfLegAllInputsCreateCommon" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}

LPXLOPER Local_InfLegAllInputsCreateCommon(
    LPXLOPER XL_startDate, 
	LPXLOPER XL_endDate,
	LPXLOPER XL_infIdx, 
    LPXLOPER XL_rcvOrPay,
    LPXLOPER XL_interpType,
	LPXLOPER XL_multiple,
	LPXLOPER XL_constant,
    LPXLOPER XL_resetFreq,
	LPXLOPER XL_dayCount,
    LPXLOPER XL_resetCalendar,
	LPXLOPER XL_fwdRule,
	LPXLOPER XL_intRule,
	LPXLOPER XL_stubRule,
	LPXLOPER XL_resetNumGap,
	LPXLOPER XL_resetDenomGap,
	LPXLOPER XL_payFreq,
	LPXLOPER XL_payGap,
	LPXLOPER XL_payCalendar,
	LPXLOPER XL_adjFirstDate,

	long finalNotionalType,
	double firstReset,
	int swapType,
	bool PersistentInXL )
{
	// this is defined first because it is used in XL macros
	static XLOPER XL_result;
	ARM_result C_result;

	
	// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	/// to avoid computation if called by the wizard
	ARM_NOCALCIFWIZ();
	
	/// Get the variables from the XLOper variables
	double C_startDate; 
	double C_endDate;
	CCString C_infIdxName;
	
	int C_rcvOrPay;
	int C_interpType;
	double C_multiple;
	double C_constant;
	int C_stubRule;
	int C_dayCount;
	int C_resetFreq;
	int C_fwdRule;
	int C_intRule;
	double C_resetNumGap;
	double C_resetDenomGap;
	CCString C_resetCalendar;
	int C_payFreq;
	double C_payGap;
	CCString C_payCalendar;
	double C_adjFirstDate;

	/// this is used by macros 
	/// and therefore this has to be defined
	static int error;
	static char* reason = "";
	/// we adjust the first date 
	double adjFirstDateDefaultValue = K_UNADJUSTED;
	double gapDefault= 0;
	double defaultMultiple = 1;
	double defaultCoMultiple = 1;
	double defaultConstant = 0;

	// read 
	XL_readNumCell(	 XL_startDate,		C_startDate,	" ARM_ERR: start Date: date expected	",			C_result);
	XL_readNumCell(	 XL_endDate,		C_endDate,		" ARM_ERR: end Date: date expected",				C_result);

	XL_readStrCell(	 XL_infIdx,	C_infIdxName,					" ARM_ERR: index name string expected", C_result );
	XL_GETRCVORPAY( XL_rcvOrPay,		C_rcvOrPay,				" ARM_ERR: receive or pay: string expected",		C_result);
	XL_GETCPIDLYINTERPWD( XL_interpType,C_interpType,"CPILINEAR"," ARM_ERR: interpType: numeric expected",			C_result);
	XL_readNumCellWD( XL_multiple, C_multiple, defaultMultiple, " ARM_ERR: multiple: numeric expected",				C_result);
	XL_readNumCellWD( XL_constant, C_constant, defaultConstant, " ARM_ERR: constant: numeric expected",				C_result);
	XL_GETFREQUENCYWD(XL_resetFreq,C_resetFreq, "-1",			" ARM_ERR: reset freq: string expected",			C_result );
	XL_GETDAYCOUNTWD( XL_dayCount, C_dayCount, "30/360",		" ARM_ERR: dayCount : string expected",				C_result );
	XL_readStrCellWD(XL_resetCalendar,C_resetCalendar, "INF",	" ARM_ERR: reset Calendar: string expected",	C_result);
	XL_GETFWDRULEWD( XL_fwdRule, C_fwdRule, "MF",				" ARM_ERR: fwdRule : string expected",				C_result );
	XL_GETINTRULEWD( XL_intRule, C_intRule, "UNADJ",			" ARM_ERR: fwdRule : string expected",				C_result );
	XL_GETCONVRULEWD( XL_stubRule, C_stubRule, "SS",			" ARM_ERR: stub Rule : string expected",			C_result );
	XL_readNumCellWD(XL_resetNumGap,C_resetNumGap,gapDefault,	" ARM_ERR: numerator reset gap: numeric expected",			C_result);
	
	/// in the case of not specified input as reset DenomGap.. take the same input as numResetGap
	XL_readNumCellWD(XL_resetDenomGap,C_resetDenomGap,C_resetNumGap, " ARM_ERR: denominator reset gap: numeric expected",			C_result);
	XL_GETFREQUENCYWD( XL_payFreq, C_payFreq, "-1",				" ARM_ERR: pay frequency: string expected",			C_result );
	XL_readNumCellWD(XL_payGap, C_payGap, gapDefault,			" ARM_ERR: pay gap: numeric expected",				C_result);
	XL_readStrCellWD(XL_payCalendar, C_payCalendar,	GETDEFAULTVALUESTR," ARM_ERR: pay Calendar: string expected",		C_result);
	XL_readNumCellWD(XL_adjFirstDate,C_adjFirstDate,adjFirstDateDefaultValue, " ARM_ERR: adjust First Date: numeric expected",C_result);
	
	/// use the concept of Functor to transfer the knowledge of
	/// a function with a context
	infLegAllInputsFunc ourFunc(																							
		C_startDate, 
		C_endDate,
		C_infIdxName, 
		swapType,
		C_rcvOrPay,
		C_interpType,
		C_multiple,
		C_constant,
		C_resetFreq,
		C_dayCount,
		C_resetCalendar,
		C_fwdRule,
		C_intRule,
		C_stubRule,
		C_resetNumGap,
		C_resetDenomGap,
		C_payFreq,
		C_payGap,
		C_payCalendar,
		C_adjFirstDate,
		finalNotionalType,
		firstReset);

	/// call the general functionourFunc
	fillXL_Result( LOCAL_INFSWAPLEG_CLASS, ourFunc, C_result, XL_result, PersistentInXL );
	
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_InfLegAllInputsCreateCommon" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}

///////////////////////////////////
/// version that takes into account 
/// previous creation of object
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_InfLegYtYCreate(
    LPXLOPER XL_startDate, 
	LPXLOPER XL_endDate,
	LPXLOPER XL_infIdx, 
    LPXLOPER XL_rcvOrPay,
    LPXLOPER XL_interpType,
	LPXLOPER XL_leverage,
	LPXLOPER XL_spread,
    LPXLOPER XL_resetFreq,
	LPXLOPER XL_dayCount,
    LPXLOPER XL_resetCalendar,
	LPXLOPER XL_fwdRule,
	LPXLOPER XL_intRule,
	LPXLOPER XL_stubRule,
	LPXLOPER XL_resetNumGap,
	LPXLOPER XL_resetDenomGap,
	LPXLOPER XL_payFreq,
	LPXLOPER XL_payGap,
	LPXLOPER XL_payCalendar,
	LPXLOPER XL_adjFirstDate,
	LPXLOPER XL_Coleverage)
{
	ADD_LOG("Local_InfLegYtYCreate");
	bool PersistentInXL		= true;
	int swapType			= K_YEARTOYEAR_LEG;
	double finalNotionalType= K_NX_NONE;
	double firstReset		= GETDEFAULTVALUE;

	return Local_InfLegAllInputsCreateCommon(
		XL_startDate, 
		XL_endDate,
		XL_infIdx, 
		XL_rcvOrPay,
		XL_interpType,
		XL_leverage,
		XL_spread,
		XL_resetFreq,
		XL_dayCount,
		XL_resetCalendar,
		XL_fwdRule,
		XL_intRule,
		XL_stubRule,
		XL_resetNumGap,
		XL_resetDenomGap,
		XL_payFreq,
		XL_payGap,
		XL_payCalendar,
		XL_adjFirstDate,
		XL_Coleverage,
		finalNotionalType,
		firstReset,
		swapType,
		PersistentInXL );
}


///////////////////////////////////
/// version for non persistentInXl
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_PXL_InfLegYtYCreate(
    LPXLOPER XL_startDate, 
	LPXLOPER XL_endDate,
	LPXLOPER XL_infIdx, 
    LPXLOPER XL_rcvOrPay,
    LPXLOPER XL_interpType,
	LPXLOPER XL_leverage,
	LPXLOPER XL_spread,
    LPXLOPER XL_resetFreq,
	LPXLOPER XL_dayCount,
    LPXLOPER XL_resetCalendar,
	LPXLOPER XL_fwdRule,
	LPXLOPER XL_intRule,
	LPXLOPER XL_stubRule,
	LPXLOPER XL_resetNumGap,
	LPXLOPER XL_resetDenomGap,
	LPXLOPER XL_payFreq,
	LPXLOPER XL_payGap,
	LPXLOPER XL_payCalendar,
	LPXLOPER XL_adjFirstDate,
	LPXLOPER XL_Coleverage)
{
	ADD_LOG("Local_PXL_InfLegYtYCreate");


	bool PersistentInXL = false;
	int swapType			= K_YEARTOYEAR_LEG;
	double finalNotionalType= K_NX_NONE;
	double firstReset		= GETDEFAULTVALUE;

	return Local_InfLegAllInputsCreateCommon(
		XL_startDate, 
		XL_endDate,
		XL_infIdx, 
		XL_rcvOrPay,
		XL_interpType,
		XL_leverage,
		XL_spread,
		XL_resetFreq,
		XL_dayCount,
		XL_resetCalendar,
		XL_fwdRule,
		XL_intRule,
		XL_stubRule,
		XL_resetNumGap,
		XL_resetDenomGap,
		XL_payFreq,
		XL_payGap,
		XL_payCalendar,
		XL_adjFirstDate,
		XL_Coleverage,
		finalNotionalType,
		firstReset,
		swapType,
		PersistentInXL );
}
/**************************************************************************************************************/

/////////////////////////////////////////////////////////////
/// central function that does the creation of the XL function
/////////////////////////////////////////////////////////////

LPXLOPER Local_InfLegZCCreateCommon(
    LPXLOPER XL_startDate, 
	LPXLOPER XL_endDate,
	LPXLOPER XL_infIdx, 
    LPXLOPER XL_rcvOrPay,
    LPXLOPER XL_interpType,
	LPXLOPER XL_multiple,
	LPXLOPER XL_constant,
    LPXLOPER XL_resetCalendar,
	LPXLOPER XL_resetNumGap,
	LPXLOPER XL_resetDenomGap,
	LPXLOPER XL_payGap,
	LPXLOPER XL_payCalendar,
	LPXLOPER XL_adjFirstDate,
	LPXLOPER XL_firstReset,
	bool PersistentInXL )
{
	// this is defined first because it is used in XL macros
	static XLOPER XL_result;
	/// Get the variables from the XLOper variables
	ARM_result C_result;

	// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	/// to avoid computation if called by the wizard
	ARM_NOCALCIFWIZ();
	
	double C_startDate; 
	double C_endDate;
	CCString C_infIdxName;
	int C_rcvOrPay;
	int C_interpType;
	double C_multiple;
	double C_constant;
	double C_resetNumGap;
	double C_resetDenomGap;
	CCString C_resetCalendar;
	double C_payGap;
	CCString C_payCalendar;
	double C_adjFirstDate;
	double C_firstReset;

	/// this is used by macros 
	/// and therefore this has to be defined
	static int error;
	static char* reason = "";
	/// we adjust the first date 
	double adjFirstDateDefaultValue = 1;
	double gapDefault= 0;
	double defaultMultiple	= 1;
	double defaultConstant	= 0;
	double defaultReset		= GETDEFAULTVALUE;

	// read 
	XL_readNumCell(	XL_startDate, C_startDate,	" ARM_ERR: start Date: date expected	", C_result);
	XL_readNumCell(	XL_endDate,	C_endDate,		" ARM_ERR: end Date: date expected",	C_result);

	//// we want to allow to have either an object or a string
	XL_readStrCell( XL_infIdx, C_infIdxName,						" ARM_ERR: index name string expected", C_result );
	XL_GETRCVORPAY( XL_rcvOrPay, C_rcvOrPay,						" ARM_ERR: receive or pay: string expected", C_result);
	XL_GETCPIDLYINTERPWD( XL_interpType, C_interpType, "CPILINEAR", " ARM_ERR: interpType: numeric expected", C_result);
	XL_readNumCellWD( XL_multiple, C_multiple, defaultMultiple,		" ARM_ERR: multiple: numeric expected",	C_result);
	XL_readNumCellWD( XL_constant, C_constant, defaultConstant,		" ARM_ERR: constant: numeric expected",	C_result);
	XL_readStrCellWD(XL_resetCalendar, C_resetCalendar, "INF",		" ARM_ERR: reset Calendar: string expected", C_result);
	XL_readNumCellWD(XL_resetNumGap, C_resetNumGap,gapDefault,		" ARM_ERR: numerator reset gap: numeric expected", C_result);
	
	/// in the case of not specified input as reset DenomGap.. take the same input as numResetGap
	XL_readNumCellWD(XL_resetDenomGap,C_resetDenomGap,C_resetNumGap, " ARM_ERR: denominator reset gap: numeric expected",			C_result);
	XL_readNumCellWD(XL_payGap, C_payGap, gapDefault,				" ARM_ERR: pay gap: numeric expected", C_result);
	XL_readStrCellWD(XL_payCalendar, C_payCalendar,	GETDEFAULTVALUESTR," ARM_ERR: pay Calendar: string expected", C_result);
	XL_readNumCellWD(XL_adjFirstDate, C_adjFirstDate,adjFirstDateDefaultValue,	" ARM_ERR: adjust First Date: numeric expected", C_result);
	XL_readNumCellWD(XL_firstReset, C_firstReset,defaultReset,					" ARM_ERR: first Reset: numeric expected",C_result);
	
	long finalNotionalType = K_NX_NONE;

	/// use the concept of Functor to transfer the knowledge of
	/// a function with a context
	infLegAllInputsFunc ourFunc(																							
		C_startDate, 
		C_endDate,
		C_infIdxName, 
		K_ZEROCOUPON_LEG,	// hard coded type
		C_rcvOrPay,
		C_interpType,
		C_multiple,
		C_constant,
		K_ZEROCOUPON,		// hard coded type
		K30_360,
		C_resetCalendar,
		K_MOD_FOLLOWING,
		K_MOD_FOLLOWING,
		K_SHORTSTART,
		C_resetNumGap,
		C_resetDenomGap,
		K_ZEROCOUPON,		// hard coded type
		C_payGap,
		C_payCalendar,
		C_adjFirstDate,
		finalNotionalType,
		C_firstReset );

	/// call the general functionourFunc
	fillXL_Result( LOCAL_INFSWAPLEG_CLASS, ourFunc, C_result, XL_result, PersistentInXL );
	
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_InfLegZCCreateCommon" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


////////////////////////////////////////////////
/// Standard addins that handles persistence
////////////////////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_InfLegZCCreate(
    LPXLOPER XL_startDate, 
	LPXLOPER XL_endDate,
	LPXLOPER XL_infIdx, 
    LPXLOPER XL_rcvOrPay,
    LPXLOPER XL_interpType,
	LPXLOPER XL_multiple,
	LPXLOPER XL_constant,
    LPXLOPER XL_resetCalendar,
	LPXLOPER XL_resetNumGap,
	LPXLOPER XL_resetDenomGap,
	LPXLOPER XL_payGap,
	LPXLOPER XL_payCalendar,
	LPXLOPER XL_adjFirstDate,
	LPXLOPER XL_firstReset
)
{
	ADD_LOG("Local_InfLegZCCreate");
	bool PersistentInXL = true;
	return Local_InfLegZCCreateCommon(
		XL_startDate, 
		XL_endDate,
		XL_infIdx, 
		XL_rcvOrPay,
		XL_interpType,
		XL_multiple,
		XL_constant,
		XL_resetCalendar,
		XL_resetNumGap,
		XL_resetDenomGap,
		XL_payGap,
		XL_payCalendar,
		XL_adjFirstDate,
		XL_firstReset,
		PersistentInXL );
}

////////////////////////////////////////////////
/// Version for Visual Basic
////////////////////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_PXL_InfLegZCCreate(
    LPXLOPER XL_startDate, 
	LPXLOPER XL_endDate,
	LPXLOPER XL_infIdx, 
    LPXLOPER XL_rcvOrPay,
    LPXLOPER XL_interpType,
	LPXLOPER XL_multiple,
	LPXLOPER XL_constant,
    LPXLOPER XL_resetCalendar,
	LPXLOPER XL_resetNumGap,
	LPXLOPER XL_resetDenomGap,
	LPXLOPER XL_payGap,
	LPXLOPER XL_payCalendar,
	LPXLOPER XL_adjFirstDate,
	LPXLOPER XL_firstReset
)
{
	ADD_LOG("Local_PXL_InfLegZCCreate");
	bool PersistentInXL = false;
	return Local_InfLegZCCreateCommon(
		XL_startDate, 
		XL_endDate,
		XL_infIdx, 
		XL_rcvOrPay,
		XL_interpType,
		XL_multiple,
		XL_constant,
		XL_resetCalendar,
		XL_resetNumGap,
		XL_resetDenomGap,
		XL_payGap,
		XL_payCalendar,
		XL_adjFirstDate,
		XL_firstReset,
		PersistentInXL );
}


///////////////////////////////////
/// version for the OAT type inflation swap
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_InfLegOATCreate(
    LPXLOPER XL_startDate, 
	LPXLOPER XL_endDate,
	LPXLOPER XL_infIdx, 
    LPXLOPER XL_rcvOrPay,
    LPXLOPER XL_interpType,
	LPXLOPER XL_leverage,
	LPXLOPER XL_spread,
    LPXLOPER XL_resetFreq,
	LPXLOPER XL_dayCount,
    LPXLOPER XL_resetCalendar,
	LPXLOPER XL_fwdRule,
	LPXLOPER XL_intRule,
	LPXLOPER XL_stubRule,
	LPXLOPER XL_resetNumGap,
	LPXLOPER XL_resetDenomGap,
	LPXLOPER XL_payFreq,
	LPXLOPER XL_payGap,
	LPXLOPER XL_payCalendar,
	LPXLOPER XL_finalNotionalType,
	LPXLOPER XL_firstReset)
{
	ADD_LOG("Local_InfLegOATCreate");
	bool PersistentInXL = true;
	int swapType = K_OATTYPE_LEG;

	// this is defined first because it is used in XL macros
	static XLOPER XL_result;

	/// Get the variables from the XLOper variables
	ARM_result C_result;

	/// this is used by macros 
	/// and therefore this has to be defined
	static int error;
	static char* reason = "";

	int C_finalNotionalType;
	XL_GETNXTYPEWD(	XL_finalNotionalType, C_finalNotionalType, "NXINFEND", " ARM_ERR: notional type: string expected", C_result);

	double C_firstReset;
	double defaultReset	= GETDEFAULTVALUE;
	XL_readNumCellWD(XL_firstReset, C_firstReset,defaultReset, " ARM_ERR: first Reset: numeric expected",C_result);

	XLOPER XL_adjFirstDateRef;
	XL_adjFirstDateRef.xltype	= xltypeNum;
	XL_adjFirstDateRef.val.num	= true;

	return Local_InfLegAllInputsCreateCommon(
		XL_startDate, 
		XL_endDate,
		XL_infIdx, 
		XL_rcvOrPay,
		XL_interpType,
		XL_leverage,
		XL_spread,
		XL_resetFreq,
		XL_dayCount,
		XL_resetCalendar,
		XL_fwdRule,
		XL_intRule,
		XL_stubRule,
		XL_resetNumGap,
		XL_resetDenomGap,
		XL_payFreq,
		XL_payGap,
		XL_payCalendar,
		&XL_adjFirstDateRef,
		C_finalNotionalType,
		C_firstReset,
		swapType,
		PersistentInXL );
}

///////////////////////////////////
/// version for non persistentInXl
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_PXL_InfLegOATCreate(
    LPXLOPER XL_startDate, 
	LPXLOPER XL_endDate,
	LPXLOPER XL_infIdx, 
    LPXLOPER XL_rcvOrPay,
    LPXLOPER XL_interpType,
	LPXLOPER XL_leverage,
	LPXLOPER XL_spread,
	LPXLOPER XL_stubRule,
	LPXLOPER XL_dayCount,
    LPXLOPER XL_resetFreq,
	LPXLOPER XL_fwdRule,
	LPXLOPER XL_intRule,
	LPXLOPER XL_resetNumGap,
	LPXLOPER XL_resetDenomGap,
    LPXLOPER XL_resetCalendar,
	LPXLOPER XL_payFreq,
	LPXLOPER XL_payGap,
	LPXLOPER XL_payCalendar,
	LPXLOPER XL_finalNotionalType,
	LPXLOPER XL_firstReset)

{
	ADD_LOG("Local_PXL_InfLegOATCreate");
	bool PersistentInXL = false;
	int swapType = K_OATTYPE_LEG;

	// this is defined first because it is used in XL macros
	static XLOPER XL_result;

	/// Get the variables from the XLOper variables
	ARM_result C_result;

	/// this is used by macros 
	/// and therefore this has to be defined
	static int error;
	static char* reason = "";

	int C_finalNotionalType;
	XL_GETNXTYPEWD(	XL_finalNotionalType, C_finalNotionalType, "NXINFEND", " ARM_ERR: notional type: string expected", C_result);

	double C_firstReset;
	double defaultReset	= GETDEFAULTVALUE;
	XL_readNumCellWD(XL_firstReset, C_firstReset,defaultReset, " ARM_ERR: first Reset: numeric expected",C_result);

	XLOPER XL_adjFirstDateRef;
	XL_adjFirstDateRef.xltype	= xltypeNum;
	XL_adjFirstDateRef.val.num	= true;

	return Local_InfLegAllInputsCreateCommon(
		XL_startDate, 
		XL_endDate,
		XL_infIdx, 
		XL_rcvOrPay,
		XL_interpType,
		XL_leverage,
		XL_spread,
		XL_resetFreq,
		XL_dayCount,
		XL_resetCalendar,
		XL_fwdRule,
		XL_intRule,
		XL_stubRule,
		XL_resetNumGap,
		XL_resetDenomGap,
		XL_payFreq,
		XL_payGap,
		XL_payCalendar,
		&XL_adjFirstDateRef,
		C_finalNotionalType,
		C_firstReset,
		swapType,
		PersistentInXL );
}


/**************************************************************************************************************/


////////////////////////////////////////////////
/// functor for InfYCMod
////////////////////////////////////////////////
class InfYCFunc : public ARMResultLong2LongFunc
{
public:
	InfYCFunc( long FirstModId, long SecondModId ) 
	:	C_FirstModId( FirstModId ),	C_SecondModId( SecondModId )
	{};

	long operator()( ARM_result& result, long objId ) 
	{
		return ARMLOCAL_infYCmod( C_FirstModId, C_SecondModId,	result,	objId );
	}
private:
	long C_FirstModId;
	long C_SecondModId;

};


////////////////////////////////////////////////
/// general function to branch
////////////////////////////////////////////////

LPXLOPER Local_2ModCommon(
	LPXLOPER XL_zeroCurve,
	LPXLOPER XL_secondModel,
	bool PersistentInXL )
{
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

	long C_zeroCurveId;
	long C_secondModId ;

	XL_GETOBJIDWD(	XL_zeroCurve,	C_zeroCurveId,	"NULL OBJECT",	" ARM_ERR: Zero Curve: Object expected",		C_result);
	XL_GETOBJIDWD(	XL_secondModel,	C_secondModId ,	"NULL OBJECT",	" ARM_ERR: Second Model: Object expected",		C_result);

	InfYCFunc ourFunc( C_zeroCurveId, C_secondModId );

	/// call the general function
	fillXL_Result( LOCAL_INFCURVMODEL_CLASS, ourFunc, C_result, XL_result, PersistentInXL );

	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_2ModCommon" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


////////////////////////////////////////////////
/// Version for XL exportation of the InfYCMod
////////////////////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_InfYCMOD(
	LPXLOPER XL_zeroCurve,
	LPXLOPER XL_discCurve
)
{
	ADD_LOG("Local_InfYCMOD");
	bool PersistentInXL = true;
	return Local_2ModCommon(
		XL_zeroCurve,
		XL_discCurve,
		PersistentInXL );
}

////////////////////////////////////////////////
/// Version for XL exportation 
////////////////////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_PXL_InfYCMOD(
	LPXLOPER XL_zeroCurve,
	LPXLOPER XL_discCurve
)
{
	ADD_LOG("Local_PXL_InfYCMOD");
	bool PersistentInXL = false;
	return Local_2ModCommon(
		XL_zeroCurve,
		XL_discCurve,
		PersistentInXL );
}




////////////////////////////////////////////////
/// Fix Zero coupon
////////////////////////////////////////////////

class fixZCFunc: public ARMResultLong2LongFunc
{
public:
	fixZCFunc(
		double startDate, 
		double endDate,
		double fixRate,
		int rcvOrPay,
		int dayCount,
		int intRule,
		int stubRule,
		int payGap,
		const CCString& payCalendar,
		long ccyId )
	:
		C_startDate( startDate ),  
		C_endDate( endDate ),
		C_fixRate( fixRate ),
		C_rcvOrPay( rcvOrPay ),
		C_intRule( intRule ),
		C_stubRule( stubRule ),
		C_dayCount( dayCount ),
		C_payGap( payGap ),
		C_payCalendar( payCalendar ),
		C_ccyId( ccyId )
	{};

	long operator()( ARM_result& result, long objId ) 
	{
		return ARMLOCAL_fixZC(
			/// context
			C_startDate, 
			C_endDate,
			C_fixRate,
			C_rcvOrPay,
			C_dayCount,
			C_intRule,
			C_stubRule,
			C_payGap,
			C_payCalendar,
			C_ccyId,
			/// arguments
			result,								
			objId );
	}

private:
	double C_startDate;
	double C_endDate;
	double C_fixRate;
	int C_rcvOrPay;
	int C_dayCount;
	int C_intRule;
	int C_stubRule;
	int C_payGap;
	CCString C_payCalendar;
	long C_ccyId;
};



////////////////////////////////////////////////
/// general function to branch
////////////////////////////////////////////////

LPXLOPER Local_FixZCCommon(
    LPXLOPER XL_startDate, 
	LPXLOPER XL_endDate,
	LPXLOPER XL_fixRate, 
    LPXLOPER XL_rcvOrPay,
	LPXLOPER XL_dayCount,
	LPXLOPER XL_intRule,
	LPXLOPER XL_stubRule,
	LPXLOPER XL_payGap,
	LPXLOPER XL_payCalendar,
	LPXLOPER XL_ccyId,
	bool PersistentInXL )
{
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

	// read
	double C_startDate; 
	double C_endDate;
	double C_fixRate;
	int C_rcvOrPay;
	int C_dayCount;
	int C_intRule;
	int C_stubRule;
	double C_payGap;
	CCString C_payCalendar;
	long C_ccyId;

	double gapDefault = 0.0;


	XL_readNumCell(	XL_startDate, C_startDate,	" ARM_ERR: start Date: date expected",	C_result);
	XL_readNumCell(	XL_endDate,	C_endDate, " ARM_ERR: end Date: date expected",	C_result);
	XL_readNumCell(	XL_fixRate,	C_fixRate, " ARM_ERR: fix Rate: numeric expected",	C_result);
	XL_GETRCVORPAY( XL_rcvOrPay, C_rcvOrPay, " ARM_ERR: receive or pay: string expected", C_result);
	XL_GETDAYCOUNTWD( XL_dayCount, C_dayCount, "30/360"," ARM_ERR: dayCount : string expected",	C_result );
	XL_GETINTRULEWD( XL_intRule, C_intRule, "UNADJ",		" ARM_ERR: fwdRule : string expected",				C_result );
	XL_GETCONVRULEWD( XL_stubRule, C_stubRule, "SS",	" ARM_ERR: stub Rule : string expected",			C_result );
	XL_readNumCellWD( XL_payGap, C_payGap, gapDefault, " ARM_ERR: pay gap: numeric expected", C_result);
	XL_readStrCellWD( XL_payCalendar, C_payCalendar, GETDEFAULTVALUESTR," ARM_ERR: pay Calendar: string expected", C_result);
	XL_GETOBJIDWD( XL_ccyId, C_ccyId, GETDEFAULTVALUESTR, " ARM_ERR: Ccy Flag: string expected", C_result );

	fixZCFunc ourFunc(
		C_startDate,
		C_endDate,
		C_fixRate,
		C_rcvOrPay,
		C_dayCount,
		C_intRule,
		C_stubRule,
		C_payGap,
		C_payCalendar,
		C_ccyId );

	/// call the general function
	fillXL_Result( LOCAL_SWAPLEG_CLASS, ourFunc, C_result, XL_result, PersistentInXL );

	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_FixZCCommon" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


////////////////////////////////////////////////
/// Version for XL exportation of the InfYCMod
////////////////////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_FixZC(
    LPXLOPER XL_startDate, 
	LPXLOPER XL_endDate,
	LPXLOPER XL_fixRate, 
    LPXLOPER XL_rcvOrPay,
	LPXLOPER XL_dayCount,
	LPXLOPER XL_intRule,
	LPXLOPER XL_stubRule,
	LPXLOPER XL_payGap,
	LPXLOPER XL_payCalendar,
	LPXLOPER XL_ccyId
)
{
	ADD_LOG("Local_FixZC");
	bool PersistentInXL = true;
	return Local_FixZCCommon(
		XL_startDate, 
		XL_endDate,
		XL_fixRate, 
		XL_rcvOrPay,
		XL_dayCount,
		XL_intRule,
		XL_stubRule,
		XL_payGap,
		XL_payCalendar,
		XL_ccyId,
		PersistentInXL );
}

////////////////////////////////////////////////
/// Version for XL exportation 
////////////////////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_PXL_FixZC(
    LPXLOPER XL_startDate, 
	LPXLOPER XL_endDate,
	LPXLOPER XL_fixRate, 
    LPXLOPER XL_rcvOrPay,
	LPXLOPER XL_intRule,
	LPXLOPER XL_stubRule,
	LPXLOPER XL_dayCount,
	LPXLOPER XL_payGap,
	LPXLOPER XL_payCalendar,
	LPXLOPER XL_ccyId
)
{
	ADD_LOG("Local_PXL_FixZC");
	bool PersistentInXL = false;
	return Local_FixZCCommon(
		XL_startDate, 
		XL_endDate,
		XL_fixRate, 
		XL_rcvOrPay,
		XL_dayCount,
		XL_intRule,
		XL_stubRule,
		XL_payGap,
		XL_payCalendar,
		XL_ccyId,
		PersistentInXL );
}





class resetManagerFunc: public ARMResultLong2LongFunc
{
public:
	resetManagerFunc(const VECTOR<CCString>& data, int nbrows, int nbcolumns )
	:	C_data(data), C_nbrows( nbrows ), C_nbcolumns( nbcolumns )
	{};

	long operator()( ARM_result& result, long objId )
	{
		return ARMLOCAL_GetData( C_data, C_nbrows, C_nbcolumns, result, objId );
	}

private:
	VECTOR<CCString> C_data;
	int C_nbrows;
	int C_nbcolumns;
};



__declspec(dllexport) LPXLOPER WINAPI Local_ResetManagerCommon(
    LPXLOPER XL_data,
	bool PersistentInXL )
{
	ADD_LOG("Local_ResetManagerCommon");
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

	long nbrows, nbcolumns;
	VECTOR<CCString> C_data;
	VECTOR<CCString> C_dataDefault;
	XL_readStrVectorAndSizeWD( XL_data,nbrows,nbcolumns,C_data,C_dataDefault," ARM_ERR: reset data: array of data expected",DOUBLE_TYPE,C_result);

	/// use the concept of Functor to transfer the knowledge of
	/// a function with a context
	resetManagerFunc ourFunc( C_data, nbrows, nbcolumns  );
	
	/// call the general function
	fillXL_Result( LOCAL_RESETMANAGER_CLASS, ourFunc, C_result, XL_result, PersistentInXL );

	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_ResetManagerCommon" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}



__declspec(dllexport) LPXLOPER WINAPI Local_ResetManager(
    LPXLOPER XL_data )
{
	ADD_LOG("Local_ResetManager");
	bool PersistentInXL = true;
	return Local_ResetManagerCommon(
		XL_data, 
		PersistentInXL );
}



__declspec(dllexport) LPXLOPER WINAPI Local_PXL_ResetManager(
    LPXLOPER XL_data )
{
	ADD_LOG("Local_PXL_ResetManager");
	bool PersistentInXL = false;
	return Local_ResetManagerCommon(
		XL_data, 
		PersistentInXL );
}



/////////////////////////////////////////////
///// Sparse vol functor
/////////////////////////////////////////////
class sparseVolCubeCreateFunc: public ARMResultLong2LongFunc
{
public:
	sparseVolCubeCreateFunc(
		double					asOfDate,
		double					lastKnownDate,
		const CCString&			indexName,
		long					Dim1Type,
		const VECTOR<CCString>&	Dim1Value,
		long					Dim2Type,
		double					Dim2Value,
		const VECTOR<double>&	strikes,
		const VECTOR<double>&	vols,
		long					strikeType,
		long					volType )
	:	
		C_asOfDate( asOfDate ), 
		C_lastKnownDate( lastKnownDate ),
		C_indexName( indexName ),
		C_Dim1Type( Dim1Type ), 
		C_Dim1Value( Dim1Value ), 
		C_Dim2Type( Dim2Type ),
		C_Dim2Value( Dim2Value ), 
		C_strikes( strikes ), 
		C_vols( vols ), 
		C_strikeType( strikeType ),
		C_volType( volType )
	{};

	virtual long operator()( ARM_result& result, long objId )
	{
		return ARMLOCAL_SparseVolCube_CreateNFill(
			C_asOfDate,
			C_lastKnownDate,
			C_indexName,
			C_Dim1Type,
			C_Dim1Value,
			C_Dim2Type,
			C_Dim2Value,
			C_strikes,
			C_vols,		
			C_strikeType,
			C_volType,
			result,
			objId );
	}

private:
	double C_asOfDate;
	double C_lastKnownDate;
	CCString C_indexName;
	long C_Dim1Type;
	VECTOR<CCString> C_Dim1Value;
	long C_Dim2Type;
	double C_Dim2Value;
	VECTOR<double> C_strikes;
	VECTOR<double> C_vols;
	long C_strikeType;
	long C_volType;
};


///////////////////////////////////////////////////
///// General function
///////////////////////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_SparseVolCube_CreateCommon(
	LPXLOPER XL_asOfDate,
	LPXLOPER XL_lastKnonwnDate,
	LPXLOPER XL_indexName,
	LPXLOPER XL_Dim1Type,
	LPXLOPER XL_Dim1Value,
	LPXLOPER XL_strikes,
	LPXLOPER XL_vols,
	LPXLOPER XL_Dim2Type,
	LPXLOPER XL_Dim2Value,
	LPXLOPER XL_strikeType,
	LPXLOPER XL_volType,
	bool PersistentInXL )
{
	ADD_LOG("Local_SparseVolCube_CreateCommon");
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

	double C_asOfDate;
	double C_lastKnownDate;
	CCString C_indexName;
	VECTOR<CCString> C_Dim1Value;
	long C_Dim1Type;
	VECTOR<double> C_strikes;
	VECTOR<double> C_vols;
	long C_Dim2Type;
	double C_Dim2Value;
	long C_strikeType;
	long C_volType;

	XL_readNumCell(	 XL_asOfDate,		C_asOfDate,		" ARM_ERR: asOfDate: date expected",			C_result);
	XL_readNumCell(	 XL_lastKnonwnDate,	C_lastKnownDate," ARM_ERR: last Known Date: date expected",		C_result);
	XL_readStrCellWD( XL_indexName,		C_indexName, GETDEFAULTVALUESTR,	" ARM_ERR: index name string expected",				C_result);
	XL_GETDIMTYPE(	 XL_Dim1Type,		C_Dim1Type,		" ARM_ERR: dim type: string expected",	C_result);
	XL_readStrVector(XL_Dim1Value,		C_Dim1Value,	" ARM_ERR: Dim1: array of numeric expected",	DOUBLE_TYPE, C_result);
	XL_GETDIMTYPE(	 XL_Dim2Type,		C_Dim2Type,		" ARM_ERR: dim type: string expected",			C_result);
	XL_readNumCell(	 XL_Dim2Value,		C_Dim2Value,	" ARM_ERR: other Dim Value: double expected",	C_result);
	XL_readNumVector(XL_strikes,		C_strikes,		" ARM_ERR: strikes: array of numeric expected",	C_result);
	XL_readNumVector(XL_vols,			C_vols,			" ARM_ERR: volatilities: array of numeric expected",C_result);
	XL_GETSTRIKETYPEWD( XL_strikeType,	C_strikeType,	"P", " ARM_ERR: strike type: string expected",	C_result);
	XL_GETVOLTYPEWD( XL_volType,		C_volType,		"S", " ARM_ERR: volatility type: string expected",C_result);

	/// use the concept of Functor to transfer the knowledge of
	/// a function with a context
	sparseVolCubeCreateFunc ourFunc(
		C_asOfDate,
		C_lastKnownDate,
		C_indexName,
		C_Dim1Type,
		C_Dim1Value,
		C_Dim2Type,
		C_Dim2Value,
		C_strikes,
		C_vols,	
		C_strikeType,
		C_volType );

	/// call the general function
	fillXL_Result( LOCAL_SPARSE_VOL_CUBE_CLASS, ourFunc, C_result, XL_result, PersistentInXL );

	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_SparseVolCube_CreateCommon" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


//////////////////////////////////////////////////////
/// Version for XL exportation of the Sparse Vol Cube
/// in Excel
//////////////////////////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_SparseVolCube_CreateNFill(
	LPXLOPER XL_asOfDate,
	LPXLOPER XL_lastKnonwnDate,
	LPXLOPER XL_indexName,
	LPXLOPER XL_Dim1Type,
	LPXLOPER XL_Dim1Value,
	LPXLOPER XL_Dim2Type,
	LPXLOPER XL_Dim2Value,
	LPXLOPER XL_strikes,
	LPXLOPER XL_vols,
	LPXLOPER XL_strikeType,
	LPXLOPER XL_volType )
{
	ADD_LOG("Local_SparseVolCube_CreateNFill");
	bool PersistentInXL = true;
	return Local_SparseVolCube_CreateCommon(
		XL_asOfDate,
		XL_lastKnonwnDate,
		XL_indexName,
		XL_Dim1Type,
		XL_Dim1Value,
		XL_strikes,
		XL_vols,
		XL_Dim2Type,
		XL_Dim2Value,
		XL_strikeType,
		XL_volType,
		PersistentInXL );
}


//////////////////////////////////////////////////////
/// Version for XL exportation of the Sparse Vol Cube
/// in VB
//////////////////////////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_PXL_SparseVolCube_CreateNFill(
	LPXLOPER XL_asOfDate,
	LPXLOPER XL_lastKnonwnDate,
	LPXLOPER XL_indexName,
	LPXLOPER XL_Dim1Type,
	LPXLOPER XL_Dim1Value,
	LPXLOPER XL_Dim2Type,
	LPXLOPER XL_Dim2Value,
	LPXLOPER XL_strikes,
	LPXLOPER XL_vols,
	LPXLOPER XL_strikeType,
	LPXLOPER XL_volType )
{
	ADD_LOG("Local_PXL_SparseVolCube_CreateNFill");
	bool PersistentInXL = false;
	long C_sparseVolCubeId = ARM_NULL_OBJECT_ID;
	return Local_SparseVolCube_CreateCommon(
		XL_asOfDate,
		XL_lastKnonwnDate,
		XL_indexName,
		XL_Dim1Type,
		XL_Dim1Value,
		XL_strikes,
		XL_vols,
		XL_Dim2Type,
		XL_Dim2Value,
		XL_strikeType,
		XL_volType,
		PersistentInXL );
}



/////////////////////////////////////////////
///// Sparse vol functor to fill an already
///// existing object
/////////////////////////////////////////////
class sparseVolCubeFillFunc: public ARMResultLong2LongFunc
{
public:
	sparseVolCubeFillFunc(
		long					Dim1Type,
		const VECTOR<CCString>&	Dim1Value,
		long					Dim2Type,
		double					Dim2Value,
		const VECTOR<double>&	strikes,
		const VECTOR<double>&	vols,
		long					sparseVolCubeId)
	:	
		C_Dim1Type( Dim1Type ),
		C_Dim1Value( Dim1Value ),
		C_Dim2Type( Dim2Type ),
		C_Dim2Value( Dim2Value ), 
		C_strikes( strikes ), 
		C_vols( vols ), 
		C_sparseVolCubeId( sparseVolCubeId )
	{};

	virtual long operator()( ARM_result& result, long objId )
	{
		return ARMLOCAL_SparseVolCube_Fill(
			C_Dim1Type,
			C_Dim1Value,
			C_Dim2Type,
			C_Dim2Value,
			C_strikes,
			C_vols,		
			C_sparseVolCubeId,
			result,
			objId );
	}

private:
	long C_Dim1Type;
	VECTOR<CCString> C_Dim1Value;
	long C_Dim2Type;
	double C_Dim2Value;
	VECTOR<double> C_strikes;
	VECTOR<double> C_vols;
	long C_sparseVolCubeId;
};



///////////////////////////////////////////////////
///// General function
///////////////////////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_SparseVolCube_FillCommon(
	LPXLOPER XL_sparseVolCubeId,
	LPXLOPER XL_Dim1Type,
	LPXLOPER XL_Dim1Value,
	LPXLOPER XL_Dim2Type,
	LPXLOPER XL_Dim2Value,
	LPXLOPER XL_strikes,
	LPXLOPER XL_vols,
	bool PersistentInXL )
{
	ADD_LOG("Local_SparseVolCube_FillCommon");
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

	long C_Dim1Type;
	VECTOR<CCString> C_Dim1Value;
	VECTOR<double> C_strikes;
	VECTOR<double> C_vols;
	long C_Dim2Type;
	double C_Dim2Value;
	long C_sparseVolCubeId;


	XL_GETOBJID( XL_sparseVolCubeId,	C_sparseVolCubeId,	" ARM_ERR: sparse vol cube: Object expected",	C_result);
	XL_GETDIMTYPE(	 XL_Dim1Type,		C_Dim1Type,			" ARM_ERR: dim type: string expected",			C_result);
	XL_readStrVector(XL_Dim1Value,		C_Dim1Value,		" ARM_ERR: Dim1: array of numeric expected",	DOUBLE_TYPE, C_result);
	XL_GETDIMTYPE(	 XL_Dim2Type,		C_Dim2Type,			" ARM_ERR: dim type: string expected",			C_result);
	XL_readNumCell(	 XL_Dim2Value,		C_Dim2Value,		" ARM_ERR: other Dim Value: double expected",	C_result);
	XL_readNumVector(XL_strikes,		C_strikes,			" ARM_ERR: strikes: array of numeric expected",	C_result);
	XL_readNumVector(XL_vols,			C_vols,				" ARM_ERR: volatilities: array of numeric expected",C_result);

	/// use the concept of Functor to transfer the knowledge of
	/// a function with a context
	sparseVolCubeFillFunc ourFunc(
		C_Dim1Type,
		C_Dim1Value,
		C_Dim2Type,
		C_Dim2Value,
		C_strikes,
		C_vols,	
		C_sparseVolCubeId );

	/// call the general function
	fillXL_Result( LOCAL_SPARSE_VOL_CUBE_CLASS, ourFunc, C_result, XL_result, PersistentInXL );

	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_SparseVolCube_FillCommon" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}
	
	
	
//////////////////////////////////////////////////////
/// Addin to fill an already existing vol cube
/// in Excel
//////////////////////////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_SparseVolCube_Fill(
	LPXLOPER XL_sparseVolCubeId,
	LPXLOPER XL_Dim1Type,
	LPXLOPER XL_Dim1Value,
	LPXLOPER XL_Dim2Type,
	LPXLOPER XL_Dim2Value,
	LPXLOPER XL_strikes,
	LPXLOPER XL_vols )
{
	ADD_LOG("Local_SparseVolCube_Fill");
	bool PersistentInXL = true;
	return Local_SparseVolCube_FillCommon(
		XL_sparseVolCubeId,
		XL_Dim1Type,
		XL_Dim1Value,
		XL_Dim2Type,
		XL_Dim2Value,
		XL_strikes,
		XL_vols,
		PersistentInXL );
}



//////////////////////////////////////////////////////
/// Version for XL exportation of the Sparse Vol Cube
/// in Visual Basic
//////////////////////////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_PXL_SparseVolCube_Fill(
	LPXLOPER XL_sparseVolCubeId,
	LPXLOPER XL_Dim1Type,
	LPXLOPER XL_Dim1Value,
	LPXLOPER XL_Dim2Type,
	LPXLOPER XL_Dim2Value,
	LPXLOPER XL_strikes,
	LPXLOPER XL_vols )
{
	ADD_LOG("Local_PXL_SparseVolCube_Fill");
	bool PersistentInXL = false;
	return Local_SparseVolCube_FillCommon(
		XL_sparseVolCubeId,
		XL_Dim1Type,
		XL_Dim1Value,
		XL_Dim2Type,
		XL_Dim2Value,
		XL_strikes,
		XL_vols,
		PersistentInXL );
}




/////////////////////////////////////////////
///// Sparse vol functor to fill an already
///// existing object
/////////////////////////////////////////////
class volCubeFromSparseVolFunc: public ARMResultLong2LongFunc
{
public:
	volCubeFromSparseVolFunc( long sparseVolCubeId)
	:	C_sparseVolCubeId( sparseVolCubeId )
	{};

	virtual long operator()( ARM_result& result, long objId )	{
		return ARMLOCAL_VolCubeFromSparseVolCube( C_sparseVolCubeId, result, objId );
	}
private:
	long C_sparseVolCubeId;
};




///////////////////////////////////////////////////
///// General function
///////////////////////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_VolCubeFromSparseVolCubeCommon(
	LPXLOPER XL_sparseVolCubeId,
	bool PersistentInXL )
{
	ADD_LOG("Local_VolCubeFromSparseVolCubeCommon");
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

	long C_sparseVolCubeId;
	XL_GETOBJID( XL_sparseVolCubeId,	C_sparseVolCubeId,	" ARM_ERR: sparse vol cube: Object expected",	C_result);

	/// use the concept of Functor to transfer the knowledge of
	/// a function with a context
	volCubeFromSparseVolFunc ourFunc( C_sparseVolCubeId );

	/// call the general function
	fillXL_Result( LOCAL_VOL_CUBE_CLASS, ourFunc, C_result, XL_result, PersistentInXL );

	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_VolCubeFromSparseVolCubeCommon" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}
	
	
	
//////////////////////////////////////////////////////
/// Addin to fill an already existing vol cube
/// in Excel
//////////////////////////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_VolCubeFromSparseVolCube(
	LPXLOPER XL_sparseVolCubeId )
{
	ADD_LOG("Local_VolCubeFromSparseVolCube");
	bool PersistentInXL = true;
	return Local_VolCubeFromSparseVolCubeCommon(
		XL_sparseVolCubeId,
		PersistentInXL );
}



//////////////////////////////////////////////////////
/// Version for XL exportation of the Sparse Vol Cube
/// in Visual Basic
//////////////////////////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_PXL_VolCubeFromSparseVolCube(
	LPXLOPER XL_sparseVolCubeId )
{
	ADD_LOG("Local_PXL_VolCubeFromSparseVolCube");
	bool PersistentInXL = false;
	return Local_VolCubeFromSparseVolCubeCommon(
		XL_sparseVolCubeId,
		PersistentInXL );
}



////////////////////////////////////////////////
/// functor for InfYCMod
////////////////////////////////////////////////
class InfBSModFunc : public ARMResultLong2LongFunc
{
public:
	InfBSModFunc( double asOfDate, long discountCurvId, long infFwdCurvId, long volCurvId,
		long correlManagerId, long IRBSModelId, long infSwoptCurveId,
		long IRSwoptCurveId )
	:	C_AsOfDate( asOfDate ), C_DiscountCurvId( discountCurvId ),
		C_InfFwdCurvId( infFwdCurvId ), C_VolCurvId( volCurvId ),
		C_correlManagerId( correlManagerId ), 
		C_IRBSModelId( IRBSModelId ),
		C_infSwoptCurveId( infSwoptCurveId ),
		C_IRSwoptCurveId( IRSwoptCurveId )
	{};

	long operator()( ARM_result& result, long objId ) 
	{
		return ARMLOCAL_infBSMod( C_AsOfDate, C_DiscountCurvId,	C_InfFwdCurvId, C_VolCurvId, 
			C_correlManagerId, C_IRBSModelId, C_infSwoptCurveId, 
			C_IRSwoptCurveId, result, objId );
	}
private:
	double C_AsOfDate;
	long C_DiscountCurvId;
	long C_InfFwdCurvId;
	long C_VolCurvId;
	long C_correlManagerId;
	long C_IRBSModelId;
	long C_infSwoptCurveId;
	long C_IRSwoptCurveId;
};


////////////////////////////////////////////////
/// general function to branch
////////////////////////////////////////////////

LPXLOPER Local_InfBSMod_Common(
	LPXLOPER XL_AsOfDate,
	LPXLOPER XL_DiscountCurvId,
	LPXLOPER XL_InfFwdCurvId,
	LPXLOPER XL_VolCurvId,
	LPXLOPER XL_CorrelManagerId,
	LPXLOPER XL_IRBSModelId,
	LPXLOPER XL_InfSwoptCurveId,
	LPXLOPER XL_IRSwoptCurveId,
	bool PersistentInXL )
{
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

		double C_AsOfDate;
		long C_DiscountCurvId;
		long C_InfFwdCurvId;
		long C_VolCurvId;
		long C_correlManagerId;
		long C_IRBSModelId;
		long C_InfSwoptCurveId;
		long C_IRSwoptCurveId;
		
		XL_readNumCell(	 XL_AsOfDate,		C_AsOfDate,							" ARM_ERR: asOfDate: date expected",			C_result);
		XL_GETOBJIDWD(	XL_DiscountCurvId,	C_DiscountCurvId,	"NULL OBJECT",	" ARM_ERR: Discount Curv: Object expected",		C_result);
		XL_GETOBJIDWD(	XL_InfFwdCurvId,	C_InfFwdCurvId ,	"NULL OBJECT",	" ARM_ERR: Inflation Fwd Curv: Object expected",C_result);
		XL_GETOBJIDWD(	XL_VolCurvId,		C_VolCurvId ,		"NULL OBJECT",	" ARM_ERR: Volatility cube: Object expected",	C_result);
		XL_GETOBJIDWD(	XL_CorrelManagerId,	C_correlManagerId ,	"NULL OBJECT",	" ARM_ERR: Correl Manager: Object expected",	C_result);
		XL_GETOBJIDWD(	XL_IRBSModelId,		C_IRBSModelId ,		"NULL OBJECT",	" ARM_ERR: IR BS Model: Object expected",		C_result);
		XL_GETOBJIDWD(	XL_InfSwoptCurveId,	C_InfSwoptCurveId,	"NULL OBJECT",	" ARM_ERR: Inflation Volatility cube: Object expected",	C_result);
		XL_GETOBJIDWD(	XL_IRSwoptCurveId,	C_IRSwoptCurveId,	"NULL OBJECT",	" ARM_ERR: Interest Rate Volatility cube: Object expected",	C_result);

		InfBSModFunc ourFunc( C_AsOfDate, C_DiscountCurvId, C_InfFwdCurvId, C_VolCurvId,
			C_correlManagerId, C_IRBSModelId, C_InfSwoptCurveId, C_IRSwoptCurveId );

		/// call the general function
		fillXL_Result( LOCAL_INFCURVMODEL_CLASS, ourFunc, C_result, XL_result, PersistentInXL );

	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_InfBSMod_Common" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


////////////////////////////////////////////////
/// Version for XL exportation of the InfYCMod
////////////////////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_InfBSMOD(
	LPXLOPER XL_AsOfDate,
	LPXLOPER XL_DiscountCurvId,
	LPXLOPER XL_InfFwdCurvId,
	LPXLOPER XL_VolCurvId,
	LPXLOPER XL_CorrelManagerId,
	LPXLOPER XL_IRBSModelId,
	LPXLOPER XL_InfSwoptCurveId,
	LPXLOPER XL_IRSwoptCurveId )
{
	ADD_LOG("Local_InfBSMOD");
	bool PersistentInXL = true;
	return Local_InfBSMod_Common(
		XL_AsOfDate,
		XL_DiscountCurvId,
		XL_InfFwdCurvId,
		XL_VolCurvId,
		XL_CorrelManagerId,
		XL_IRBSModelId,
		XL_InfSwoptCurveId,
		XL_IRSwoptCurveId,
		PersistentInXL );
}



////////////////////////////////////////////////
/// Version for XL exportation 
////////////////////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_PXL_InfBSMOD(
	LPXLOPER XL_AsOfDate,
	LPXLOPER XL_DiscountCurvId,
	LPXLOPER XL_InfFwdCurvId,
	LPXLOPER XL_VolCurvId,
	LPXLOPER XL_CorrelManagerId,
	LPXLOPER XL_IRBSModelId,
	LPXLOPER XL_InfSwoptCurveId,
	LPXLOPER XL_IRSwoptCurveId )
{
	ADD_LOG("Local_PXL_InfBSMOD");
	bool PersistentInXL = false;
	return Local_InfBSMod_Common(
		XL_AsOfDate,
		XL_DiscountCurvId,
		XL_InfFwdCurvId,
		XL_VolCurvId,
		XL_CorrelManagerId,
		XL_IRBSModelId,
		XL_InfSwoptCurveId,
		XL_IRSwoptCurveId,
		PersistentInXL );
}


//////////////////////////////////////////////////
/// functor for the inflation cap and floor
//////////////////////////////////////////////////

class infCapFloorFunc: public ARMResultLong2LongFunc
{
public:
	infCapFloorFunc(
		double startDate,
		double endDate,
		const CCString& infIdxName,
		int capOrFloor,
		double strike,
		int swapType,
		int rcvOrPay,
	    int interpType,
		double leverage,
		double spread,
		int resetFreq,
		int dayCount,
		const CCString& resetCalendar,
		int fwdRule,
		int intRule,
		int stubRule,
		int resetGap,
		int payFreq,
		int payGap,
		const CCString& payCalendar,
		int adjFirstDate,
		double firstReset = GETDEFAULTVALUE )
	:
		C_startDate( startDate ),  
		C_endDate( endDate ),
		C_infIdxName( infIdxName ),
		C_capOrFloor( capOrFloor ),
		C_strike( strike ),
		C_swapType( swapType ),
		C_rcvOrPay( rcvOrPay ),
		C_interpType( interpType ),
		C_leverage( leverage ),
		C_spread( spread ),
		C_resetFreq( resetFreq ),
		C_dayCount( dayCount ),
		C_resetCalendar( resetCalendar ),
		C_fwdRule( fwdRule ),
		C_intRule( intRule ),
		C_stubRule( stubRule ),
		C_resetGap( resetGap ),
		C_payFreq( payFreq ),
		C_payGap( payGap ),
		C_payCalendar( payCalendar ),
		C_adjFirstDate( adjFirstDate ),
		C_firstReset( firstReset )
	{};

	long operator()( ARM_result& result, long objId ) 
	{
		return ARMLOCAL_InfCapFloor_Create(
			/// context
			C_startDate, 
			C_endDate,
			C_infIdxName,
			C_capOrFloor,
			C_strike,
			C_leverage,
			C_spread,
			C_swapType,
			C_rcvOrPay,
			C_interpType,
			C_resetFreq,
			C_dayCount,
			C_resetCalendar,
			C_fwdRule,
			C_intRule,
			C_stubRule,
			C_resetGap,
			C_payFreq,
			C_payGap,
			C_payCalendar,
			C_adjFirstDate,
			C_firstReset,
			ARM_NULL_OBJECT,
			/// arguments
			result,								
			objId );
	}
private:
	double C_startDate;									
	double C_endDate;
	CCString C_infIdxName;
	int C_capOrFloor;
	double C_strike;
	int C_swapType;
	int C_rcvOrPay;
	int C_interpType;
	double C_leverage;
	double C_spread;
	int C_resetFreq;
	int C_dayCount;
	CCString C_resetCalendar;
	int C_fwdRule;
	int C_intRule;
	int C_stubRule;
	int C_resetGap;
	int C_payFreq;
	int C_payGap;
	CCString C_payCalendar;
	int C_adjFirstDate;
	double C_firstReset;
};




/////////////////////////////////////////////////////////////
/// central function that does the creation of the XL function
/// Creates an inflation cap and/or floor 
/////////////////////////////////////////////////////////////
LPXLOPER Local_InfCapFloorCreateCommon(
    LPXLOPER XL_startDate, 
	LPXLOPER XL_endDate,
	LPXLOPER XL_infIdx, 
	LPXLOPER XL_capOrFloor,
	LPXLOPER XL_strike,
    LPXLOPER XL_swapType,
    LPXLOPER XL_rcvOrPay,
    LPXLOPER XL_interpType,
	LPXLOPER XL_leverage,
	LPXLOPER XL_spread,
    LPXLOPER XL_resetFreq,
	LPXLOPER XL_dayCount,
    LPXLOPER XL_resetCalendar,
	LPXLOPER XL_fwdRule,
	LPXLOPER XL_intRule,
	LPXLOPER XL_stubRule,
	LPXLOPER XL_resetGap,
	LPXLOPER XL_payFreq,
	LPXLOPER XL_payGap,
	LPXLOPER XL_payCalendar,
	bool PersistentInXL )
{
	// this is defined first because it is used in XL macros
	static XLOPER XL_result;
	/// Get the variables from the XLOper variables
	ARM_result C_result;

	// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	/// to avoid computation if called by the wizard
	ARM_NOCALCIFWIZ();
	
	double C_startDate; 
	double C_endDate;
	CCString C_infIdxName;
	int C_capOrFloor;
	double C_strike;
	double C_leverage;
	double C_spread;
	int C_swapType;
	int C_rcvOrPay;
	int C_interpType;
	int C_stubRule;
	int C_dayCount;
	int C_resetFreq;
	int C_fwdRule;
	int C_intRule;
	double C_resetGap;
	CCString C_resetCalendar;
	int C_payFreq;
	double C_payGap;
	CCString C_payCalendar;

	/// this is used by macros 
	/// and therefore this has to be defined
	static int error;
	static char* reason = "";
	/// we adjust the first date 
	double adjFirstDateDefaultValue = K_UNADJUSTED;
	double gapDefault= 0;
	double defaultReset		= GETDEFAULTVALUE;
	double defaultLeverage	= 1.0;
	double defaultSpread	= 0.0;

	// read 
	XL_readNumCell(	 XL_startDate,		C_startDate,				" ARM_ERR: start Date: date expected	",		C_result);
	XL_readNumCell(	 XL_endDate,		C_endDate,					" ARM_ERR: end Date: date expected",			C_result);
	XL_readStrCell(	 XL_infIdx,	C_infIdxName,						" ARM_ERR: index name string expected",			C_result );
	XL_GETCAPORFLOOR(XL_capOrFloor,	C_capOrFloor,					" ARM_ERR: cap or Floor: string expected",		C_result );
	XL_readNumCell(	 XL_strike,		C_strike,						" ARM_ERR: strike: numeric expected",			C_result);
	XL_GETINFSWAPTYPE(XL_swapType, C_swapType,						" ARM_ERR: swap type: numeric expected",			C_result);
	XL_GETRCVORPAY( XL_rcvOrPay,		C_rcvOrPay,					" ARM_ERR: receive or pay: string expected",		C_result);
	XL_GETCPIDLYINTERPWD( XL_interpType,C_interpType,"CPILINEAR",	" ARM_ERR: interpType: numeric expected",			C_result);
	XL_readNumCellWD(	 XL_leverage,	C_leverage,	defaultLeverage," ARM_ERR: leverage: numeric expected",			C_result);
	XL_readNumCellWD(	 XL_spread,		C_spread, defaultSpread,	" ARM_ERR: spread: numeric expected",			C_result);
	XL_GETFREQUENCYWD(XL_resetFreq,C_resetFreq, "-1",				" ARM_ERR: reset freq: string expected",			C_result );
	XL_GETDAYCOUNTWD( XL_dayCount, C_dayCount, "30/360",			" ARM_ERR: dayCount : string expected",				C_result );
	XL_readStrCellWD(XL_resetCalendar,C_resetCalendar, "INF",		" ARM_ERR: reset Calendar: string expected",	C_result);
	XL_GETFWDRULEWD( XL_fwdRule, C_fwdRule, "MF",					" ARM_ERR: fwdRule : string expected",				C_result );
	XL_GETINTRULEWD( XL_intRule, C_intRule, "UNADJ",				" ARM_ERR: fwdRule : string expected",				C_result );
	XL_GETCONVRULEWD( XL_stubRule, C_stubRule, "SS",				" ARM_ERR: stub Rule : string expected",			C_result );
	XL_readNumCellWD(XL_resetGap,C_resetGap,gapDefault,				" ARM_ERR: reset gap: numeric expected",			C_result);
	XL_GETFREQUENCYWD( XL_payFreq, C_payFreq, "-1",					" ARM_ERR: pay frequency: string expected",			C_result );
	XL_readNumCellWD(XL_payGap, C_payGap, gapDefault,				" ARM_ERR: pay gap: numeric expected",				C_result);
	XL_readStrCellWD(XL_payCalendar, C_payCalendar,	GETDEFAULTVALUESTR," ARM_ERR: pay Calendar: string expected",		C_result);

	/// use the concept of Functor to transfer the knowledge of
	/// a function with a context
	infCapFloorFunc ourFunc(																							
		C_startDate, 
		C_endDate,
		C_infIdxName,
		C_capOrFloor,
		C_strike,
		C_swapType,
		C_rcvOrPay,
		C_interpType,
		C_leverage,
		C_spread,
		C_resetFreq,
		C_dayCount,
		C_resetCalendar,
		C_fwdRule,
		C_intRule,
		C_stubRule,
		C_resetGap,
		C_payFreq,
		C_payGap,
		C_payCalendar,
		K_UNADJUSTED,
		GETDEFAULTVALUE );

	/// call the general functionourFunc
	fillXL_Result( LOCAL_INFCAPFLOOR_CLASS, ourFunc, C_result, XL_result, PersistentInXL );
	
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_InfCapFloorCreateCommon" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}





//////////////////////////////////////////////////////
/// Addin to create an inflation cap floor
//////////////////////////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_InfCapFloorCreate(
    LPXLOPER XL_startDate, 
	LPXLOPER XL_endDate,
	LPXLOPER XL_infIdx, 
	LPXLOPER XL_capOrFloor,
	LPXLOPER XL_strike,
    LPXLOPER XL_swapType,
    LPXLOPER XL_rcvOrPay,
    LPXLOPER XL_interpType,
	LPXLOPER XL_leverage,
	LPXLOPER XL_spread,
    LPXLOPER XL_resetFreq,
	LPXLOPER XL_dayCount,
    LPXLOPER XL_resetCalendar,
	LPXLOPER XL_fwdRule,
	LPXLOPER XL_intRule,
	LPXLOPER XL_stubRule,
	LPXLOPER XL_resetGap,
	LPXLOPER XL_payFreq,
	LPXLOPER XL_payGap,
	LPXLOPER XL_payCalendar )
{
	ADD_LOG("Local_InfCapFloorCreate");
	bool PersistentInXL = true;
	return Local_InfCapFloorCreateCommon(
		XL_startDate, 
		XL_endDate,
		XL_infIdx, 
		XL_capOrFloor,
		XL_strike,		
		XL_swapType,
		XL_rcvOrPay,		
		XL_interpType,		
		XL_leverage,		
		XL_spread,		
		XL_resetFreq,		
		XL_dayCount,		
		XL_resetCalendar,		
		XL_fwdRule,		
		XL_intRule,
		XL_stubRule,
		XL_resetGap,
		XL_payFreq,		
		XL_payGap,		
		XL_payCalendar,		
		PersistentInXL);
}



//////////////////////////////////////////////////////
/// Version for XL exportation of the function to create
/// an inflation cap and floor in Visual Basic
//////////////////////////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_PXL_InfCapFloorCreate(
    LPXLOPER XL_startDate, 
	LPXLOPER XL_endDate,
	LPXLOPER XL_infIdx, 
	LPXLOPER XL_capOrFloor,
	LPXLOPER XL_strike,
    LPXLOPER XL_swapType,
    LPXLOPER XL_rcvOrPay,
    LPXLOPER XL_interpType,
	LPXLOPER XL_leverage,
	LPXLOPER XL_spread,
    LPXLOPER XL_resetFreq,
	LPXLOPER XL_dayCount,
    LPXLOPER XL_resetCalendar,
	LPXLOPER XL_fwdRule,
	LPXLOPER XL_intRule,
	LPXLOPER XL_stubRule,
	LPXLOPER XL_resetGap,
	LPXLOPER XL_payFreq,
	LPXLOPER XL_payGap,
	LPXLOPER XL_payCalendar )
{
	ADD_LOG("Local_PXL_InfCapFloorCreate");
	bool PersistentInXL = false;
	return Local_InfCapFloorCreateCommon(
		XL_startDate, 
		XL_endDate,
		XL_infIdx, 
		XL_capOrFloor,
		XL_strike,		
		XL_swapType,
		XL_rcvOrPay,		
		XL_interpType,		
		XL_leverage,		
		XL_spread,		
		XL_resetFreq,		
		XL_dayCount,		
		XL_resetCalendar,		
		XL_fwdRule,		
		XL_intRule,		
		XL_stubRule,		
		XL_resetGap,
		XL_payFreq,		
		XL_payGap,		
		XL_payCalendar,		
		PersistentInXL);
}

/*
class infCapFloorFunctor: public ARMResultLong2LongFunc
{
public:
	infCapFloorFunctor(
		double startDate,
		double endDate,
		const CCString& infIdxName,
		int capOrFloor,
		double strike,
	    int interpType,
		int resetFreq,
		int dayCount,
		const CCString& resetCalendar,
		int fwdRule,
		int intRule,
		int stubRule,
		int resetGap,
		int payFreq,
		int payGap,
		const CCString& payCalendar,
		double firstReset = GETDEFAULTVALUE )
	:
		C_startDate( startDate ),  
		C_endDate( endDate ),
		C_infIdxName( infIdxName ),
		C_capOrFloor( capOrFloor ),
		C_strike( strike ),
		C_interpType( interpType ),
		C_resetFreq( resetFreq ),
		C_dayCount( dayCount ),
		C_resetCalendar( resetCalendar ),
		C_fwdRule( fwdRule ),
		C_intRule( intRule ),
		C_stubRule( stubRule ),
		C_resetGap( resetGap ),
		C_payFreq( payFreq ),
		C_payGap( payGap ),
		C_payCalendar( payCalendar ),
		C_firstReset( firstReset )
	{};

	long operator()( ARM_result& result, long objId ) 
	{
		return ARMLOCAL_InfCapFloor_Load(
			/// context
			C_startDate, 
			C_endDate,
			C_infIdxName,
			C_capOrFloor,
			C_strike,
			C_interpType,
			C_resetFreq,
			C_dayCount,
			C_resetCalendar,
			C_fwdRule,
			C_intRule,
			C_stubRule,
			C_resetGap,
			C_payFreq,
			C_payGap,
			C_payCalendar,
			C_firstReset,
			ARM_NULL_OBJECT,
			/// arguments
			result,								
			objId );
	}
private:
	double C_startDate;									
	double C_endDate;
	CCString C_infIdxName;
	int C_capOrFloor;
	double C_strike;
	int C_interpType;
	int C_resetFreq;
	int C_dayCount;
	CCString C_resetCalendar;
	int C_fwdRule;
	int C_intRule;
	int C_stubRule;
	int C_resetGap;
	int C_payFreq;
	int C_payGap;
	CCString C_payCalendar;
	double C_firstReset;
};



/////////////////////////////////////////////////////////////
/// central function that does the creation of the XL function
/// Creates an inflation cap and/or floor 
/////////////////////////////////////////////////////////////
LPXLOPER Local_InfCapFloorLoadCommon(
    LPXLOPER XL_startDate, 
	LPXLOPER XL_endDate,
	LPXLOPER XL_infIdx, 
	LPXLOPER XL_capOrFloor,
	LPXLOPER XL_strike,
    LPXLOPER XL_interpType,
    LPXLOPER XL_resetFreq,
	LPXLOPER XL_dayCount,
    LPXLOPER XL_resetCalendar,
	LPXLOPER XL_fwdRule,
	LPXLOPER XL_intRule,
	LPXLOPER XL_stubRule,
	LPXLOPER XL_resetGap,
	LPXLOPER XL_payFreq,
	LPXLOPER XL_payGap,
	LPXLOPER XL_payCalendar,
	bool PersistentInXL )
{
	// this is defined first because it is used in XL macros
	static XLOPER XL_result;
	/// Get the variables from the XLOper variables
	ARM_result C_result;

	// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	/// to avoid computation if called by the wizard
	ARM_NOCALCIFWIZ();
	
	double C_startDate; 
	double C_endDate;
	CCString C_infIdxName;
	int C_capOrFloor;
	double C_strike;
	int C_interpType;
	int C_stubRule;
	int C_dayCount;
	int C_resetFreq;
	int C_fwdRule;
	int C_intRule;
	double C_resetGap;
	CCString C_resetCalendar;
	int C_payFreq;
	double C_payGap;
	CCString C_payCalendar;

	/// this is used by macros 
	/// and therefore this has to be defined
	static int error;
	static char* reason = "";
	/// we adjust the first date 
	double adjFirstDateDefaultValue = K_UNADJUSTED;
	double gapDefault= 0;
	double defaultReset		= GETDEFAULTVALUE;
	double defaultLeverage	= 1.0;
	double defaultSpread	= 0.0;

	// read 
	XL_readNumCell(	 XL_startDate,		C_startDate,				" ARM_ERR: start Date: date expected	",		C_result);
	XL_readNumCell(	 XL_endDate,		C_endDate,					" ARM_ERR: end Date: date expected",			C_result);
	XL_readStrCell(	 XL_infIdx,	C_infIdxName,						" ARM_ERR: index name string expected",			C_result );
	XL_GETCAPORFLOOR(XL_capOrFloor,	C_capOrFloor,					" ARM_ERR: cap or Floor: string expected",		C_result );
	XL_readNumCell(	 XL_strike,		C_strike,						" ARM_ERR: strike: numeric expected",			C_result);
	XL_GETCPIDLYINTERPWD( XL_interpType,C_interpType,"CPILINEAR",	" ARM_ERR: interpType: numeric expected",			C_result);
	XL_GETFREQUENCYWD(XL_resetFreq,C_resetFreq, "-1",				" ARM_ERR: reset freq: string expected",			C_result );
	XL_GETDAYCOUNTWD( XL_dayCount, C_dayCount, "30/360",			" ARM_ERR: dayCount : string expected",				C_result );
	XL_readStrCellWD(XL_resetCalendar,C_resetCalendar, "INF",		" ARM_ERR: reset Calendar: string expected",	C_result);
	XL_GETFWDRULEWD( XL_fwdRule, C_fwdRule, "MF",					" ARM_ERR: fwdRule : string expected",				C_result );
	XL_GETINTRULEWD( XL_intRule, C_intRule, "UNADJ",				" ARM_ERR: fwdRule : string expected",				C_result );
	XL_GETCONVRULEWD( XL_stubRule, C_stubRule, "SS",				" ARM_ERR: stub Rule : string expected",			C_result );
	XL_readNumCellWD(XL_resetGap,C_resetGap,gapDefault,				" ARM_ERR: reset gap: numeric expected",			C_result);
	XL_GETFREQUENCYWD( XL_payFreq, C_payFreq, "-1",					" ARM_ERR: pay frequency: string expected",			C_result );
	XL_readNumCellWD(XL_payGap, C_payGap, gapDefault,				" ARM_ERR: pay gap: numeric expected",				C_result);
	XL_readStrCellWD(XL_payCalendar, C_payCalendar,	GETDEFAULTVALUESTR," ARM_ERR: pay Calendar: string expected",		C_result);

	/// use the concept of Functor to transfer the knowledge of
	/// a function with a context
	infCapFloorFunctor ourFunc(																							
		C_startDate, 
		C_endDate,
		C_infIdxName,
		C_capOrFloor,
		C_strike,
		C_interpType,
		C_resetFreq,
		C_dayCount,
		C_resetCalendar,
		C_fwdRule,
		C_intRule,
		C_stubRule,
		C_resetGap,
		C_payFreq,
		C_payGap,
		C_payCalendar,
		K_UNADJUSTED,
		GETDEFAULTVALUE );

	/// call the general functionourFunc
	fillXL_Result( LOCAL_INFCAPFLOOR_CLASS, ourFunc, C_result, XL_result, PersistentInXL );
	
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_InfCapFloorCreateCommon" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}





//////////////////////////////////////////////////////
/// Addin to create an inflation cap floor
//////////////////////////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_InfCapFloorLoad(
    LPXLOPER XL_startDate, 
	LPXLOPER XL_endDate,
	LPXLOPER XL_infIdx, 
	LPXLOPER XL_capOrFloor,
	LPXLOPER XL_strike,
    LPXLOPER XL_interpType,
    LPXLOPER XL_resetFreq,
	LPXLOPER XL_dayCount,
    LPXLOPER XL_resetCalendar,
	LPXLOPER XL_fwdRule,
	LPXLOPER XL_intRule,
	LPXLOPER XL_stubRule,
	LPXLOPER XL_resetGap,
	LPXLOPER XL_payFreq,
	LPXLOPER XL_payGap,
	LPXLOPER XL_payCalendar )
{
	ADD_LOG("Local_InfCapFloorLoad");
	bool PersistentInXL = true;
	return Local_InfCapFloorLoadCommon(
		XL_startDate, 
		XL_endDate,
		XL_infIdx, 
		XL_capOrFloor,
		XL_strike,		
		XL_interpType,		
		XL_resetFreq,		
		XL_dayCount,		
		XL_resetCalendar,		
		XL_fwdRule,		
		XL_intRule,
		XL_stubRule,
		XL_resetGap,
		XL_payFreq,		
		XL_payGap,		
		XL_payCalendar,		
		PersistentInXL);
}



//////////////////////////////////////////////////////
/// Version for XL exportation of the function to create
/// an inflation cap and floor in Visual Basic
//////////////////////////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_PXL_InfCapFloorLoad(
    LPXLOPER XL_startDate, 
	LPXLOPER XL_endDate,
	LPXLOPER XL_infIdx, 
	LPXLOPER XL_capOrFloor,
	LPXLOPER XL_strike,
    LPXLOPER XL_interpType,
    LPXLOPER XL_resetFreq,
	LPXLOPER XL_dayCount,
    LPXLOPER XL_resetCalendar,
	LPXLOPER XL_fwdRule,
	LPXLOPER XL_intRule,
	LPXLOPER XL_stubRule,
	LPXLOPER XL_resetGap,
	LPXLOPER XL_payFreq,
	LPXLOPER XL_payGap,
	LPXLOPER XL_payCalendar )
{
	ADD_LOG("Local_PXL_InfCapFloorCreate");
	bool PersistentInXL = false;
	return Local_InfCapFloorLoadCommon(
		XL_startDate, 
		XL_endDate,
		XL_infIdx, 
		XL_capOrFloor,
		XL_strike,		
		XL_interpType,		
		XL_resetFreq,		
		XL_dayCount,		
		XL_resetCalendar,		
		XL_fwdRule,		
		XL_intRule,		
		XL_stubRule,		
		XL_resetGap,
		XL_payFreq,		
		XL_payGap,		
		XL_payCalendar,		
		PersistentInXL);
}


*/
/////////////////////////////////////////
/// Functor to create a correlation matrix
/////////////////////////////////////////
class correlMatFunc: public ARMResultLong2LongFunc
{
public:
	correlMatFunc( double asOfDate, 
		const VECTOR<CCString>& X, 
		const VECTOR<CCString>& Y,
		const VECTOR<double>&   Z )
	:	C_asOfDate( asOfDate ), C_X( X ), C_Y( Y ), C_Z( Z )
	{};

	long operator()( ARM_result& result, long objId )
	{
		return ARMLOCAL_CorrelMat_Create( C_asOfDate, C_X, C_Y, C_Z, result, objId );
	}

private:
	double	C_asOfDate;
	VECTOR<CCString> C_X;
	VECTOR<CCString> C_Y;
	VECTOR<double>	 C_Z;
};


/////////////////////////////////////////
/// Correlation matrix common function
/////////////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_CorrelMatCommon(
    LPXLOPER XL_asOfDate,
    LPXLOPER XL_X,
    LPXLOPER XL_Y,
    LPXLOPER XL_Z,
	bool PersistentInXL )
{
	ADD_LOG("Local_CorrelMatCommon");
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

	double	C_asOfDate;
	VECTOR<CCString> C_X;
	VECTOR<CCString> C_Y;
	VECTOR<double> C_Z;

	XL_readNumCell(	 XL_asOfDate,		C_asOfDate,		" ARM_ERR: asOfDate: date expected",				C_result);
	XL_readStrVector(XL_X,				C_X,			" ARM_ERR: X: array of numeric/maturity expected",	DOUBLE_TYPE, C_result);
	XL_readStrVector(XL_Y,				C_Y,			" ARM_ERR: Y: array of numeric/maturity expected",	DOUBLE_TYPE, C_result);
	XL_readNumVector(XL_Z,				C_Z,			" ARM_ERR: Z: matrix of numeric expected",			C_result);

	/// use the concept of Functor to transfer the knowledge of
	/// a function with a context
	correlMatFunc ourFunc( C_asOfDate, C_X, C_Y, C_Z );
	
	/// call the general function
	fillXL_Result( LOCAL_CORRELMATRIX_CLASS, ourFunc, C_result, XL_result, PersistentInXL );

	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_CorrelMatCommon" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}



//////////////////////////////////////////////////////
/// Addin to create a correlation matrix
//////////////////////////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_CorrelMat_Create(
    LPXLOPER XL_asOfDate,
    LPXLOPER XL_X,
    LPXLOPER XL_Y,
    LPXLOPER XL_Z )
{
	ADD_LOG("Local_CorrelMat_Create");
	bool PersistentInXL = true;
	return Local_CorrelMatCommon(
		XL_asOfDate,		
		XL_X,		
		XL_Y,		
		XL_Z,
		PersistentInXL );
}



//////////////////////////////////////////////////////
/// Addin to create a correlation matrix
/// Version for VBA
//////////////////////////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_PXL_CorrelMat_Create(
    LPXLOPER XL_asOfDate,
    LPXLOPER XL_X,
    LPXLOPER XL_Y,
    LPXLOPER XL_Z )
{
	ADD_LOG("Local_PXL_CorrelMat_Create");
	bool PersistentInXL = false;
	return Local_CorrelMatCommon(
		XL_asOfDate,		
		XL_X,		
		XL_Y,		
		XL_Z,
		PersistentInXL );
}


/////////////////////////////////////////
/// Functor to create a correlation manager
/////////////////////////////////////////
class correlManagerFunc: public ARMResultLong2LongFunc
{
public:
	correlManagerFunc( const CCString& mktTag, const CCString& intraMktTag, 
		double asOfDate, const VECTOR<CCString>& X,
		const VECTOR<CCString>& Y, const VECTOR<double>& Z )
	:	C_mktTag( mktTag ), C_intraMktTag( intraMktTag ), 
		C_asOfDate( asOfDate ), C_X( X ), C_Y( Y ), C_Z( Z )
	{};
	long operator()( ARM_result& result, long objId )
	{
		return ARMLOCAL_CorrelManager_Create( C_mktTag, C_intraMktTag,
			C_asOfDate, C_X, C_Y, C_Z, result, objId );
	}

private:
	CCString C_mktTag;
	CCString C_intraMktTag;
	double	C_asOfDate;
	VECTOR<CCString> C_X;
	VECTOR<CCString> C_Y;
	VECTOR<double> C_Z;
};


/////////////////////////////////////////
/// Correlation manager create common function
/////////////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_CorrelManagerCommon(
    LPXLOPER XL_mktTag,
    LPXLOPER XL_intraMktTag,
    LPXLOPER XL_asOfDate,
    LPXLOPER XL_X,
    LPXLOPER XL_Y,
    LPXLOPER XL_Z,
	bool PersistentInXL )
{
	ADD_LOG("Local_CorrelManagerCommon");
	/// this is defined first because it is used in XL macros
	static XLOPER XL_result;
	/// Get the variables from the XLOper variables
	ARM_result C_result;
	
	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
		ARM_NOCALCIFWIZ();
		
		

		// error
		static int error;
		static char* reason = "";

		double	C_asOfDate;
		CCString C_mktTag;
		CCString C_intraMktTag;
		VECTOR<CCString> C_X;
		VECTOR<CCString> C_Y;
		VECTOR<double>	 C_Z;

		XL_readNumCell(	 XL_asOfDate,		C_asOfDate,		" ARM_ERR: asOfDate: date expected",			C_result);
		XL_readStrCell(  XL_mktTag,			C_mktTag,		" ARM_ERR: mkt Tag: string expected",			C_result);
		XL_readStrCell(	 XL_intraMktTag,	C_intraMktTag,	" ARM_ERR: intra Mkt Tag: string expected",		C_result);
		XL_readStrVector(XL_X,				C_X,			" ARM_ERR: X: array of numeric/maturity expected",	DOUBLE_TYPE, C_result);
		XL_readStrVector(XL_Y,				C_Y,			" ARM_ERR: Y: array of numeric/maturity expected",	DOUBLE_TYPE, C_result);
		XL_readNumVector(XL_Z,              C_Z,			" ARM_ERR: Z: matrix of numeric expected",C_result);


		/// use the concept of Functor to transfer the knowledge of
		/// a function with a context
		// reconfigure C_intraMktTag according to ARM conventions
		// ie: EURIBOR12M -> EURIBOR1Y
    
		VECTOR<CCString> listString;
		C_intraMktTag.Parser ('_', listString);
		string temp= "";
		long irType;
		VECTOR<CCString>::iterator iter= listString.begin();
		while( iter!=listString.end() ){        
			if (iter!= listString.begin())
				temp += string("_");

			irType=ARM_ConvIrTypeWithoutException(*iter);
			if ( irType==ARM_DEFAULT_ERR)
				temp += string(*iter);
			else
				temp += string(ARM_ParamView::GetMappingName(S_INDEX_TYPES, irType));

			iter++;
		}
 
    
		correlManagerFunc ourFunc( C_mktTag, CCString(temp.c_str()), C_asOfDate, C_X, C_Y, C_Z );
		
		/// call the general function
		fillXL_Result( LOCAL_CORRELMANAGER_CLASS, ourFunc, C_result, XL_result, PersistentInXL );

	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_CorrelManagerCommon" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}



//////////////////////////////////////////////////////
/// Addin to create a correlation manager
//////////////////////////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_CorrelManager_Create(
    LPXLOPER XL_mktTag,
    LPXLOPER XL_intraMktTag,
    LPXLOPER XL_asOfDate,
    LPXLOPER XL_X,
    LPXLOPER XL_Y,
    LPXLOPER XL_Z )
{
	ADD_LOG("Local_CorrelManager_Create");
	bool PersistentInXL = true;
	return Local_CorrelManagerCommon(
		XL_mktTag,
		XL_intraMktTag,		
		XL_asOfDate,		
		XL_X,		
		XL_Y,		
		XL_Z,		
		PersistentInXL );
}



//////////////////////////////////////////////////////
/// Addin to create a correlation manager
/// Version for VBA
//////////////////////////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_PXL_CorrelManager_Create(
    LPXLOPER XL_mktTag,
    LPXLOPER XL_intraMktTag,
    LPXLOPER XL_asOfDate,
    LPXLOPER XL_X,
    LPXLOPER XL_Y,
    LPXLOPER XL_Z )
{
	ADD_LOG("Local_PXL_CorrelManager_Create");
	bool PersistentInXL = false;
	return Local_CorrelManagerCommon(
		XL_mktTag,
		XL_intraMktTag,		
		XL_asOfDate,		
		XL_X,		
		XL_Y,		
		XL_Z,		
		PersistentInXL );
}








/////////////////////////////////////////
/// Functor to create a correlation manager from
/// correlMat
/////////////////////////////////////////
class correlManagerFromMatFunc: public ARMResultLong2LongFunc
{
public:
	correlManagerFromMatFunc( const CCString& mktTag, const CCString& intraMktTag, long correlMatId )
	:	C_mktTag( mktTag ), C_intraMktTag( intraMktTag ),  C_correlMatId( correlMatId )
	{};
	long operator()( ARM_result& result, long objId )
	{
		return ARMLOCAL_CorrelManager_CreateFromMat( C_mktTag, C_intraMktTag,
			C_correlMatId, result, objId );
	}

private:
	CCString C_mktTag;
	CCString C_intraMktTag;
	long C_correlMatId;
};


/////////////////////////////////////////
/// Correlation manager create common function
/// from a correlation matrix
/////////////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_CorrelManagerFromMatCommon(
    LPXLOPER XL_mktTag,
    LPXLOPER XL_intraMktTag,
    LPXLOPER XL_correlMatId,
	bool PersistentInXL )
{
	ADD_LOG("Local_CorrelManagerFromMatCommon");
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

		CCString C_mktTag;
		CCString C_intraMktTag;
		long C_correlMatId;

		XL_readStrCell(  XL_mktTag,			C_mktTag,		" ARM_ERR: mkt Tag: string expected",			C_result);
		XL_readStrCell(	 XL_intraMktTag,	C_intraMktTag,	" ARM_ERR: intra Mkt Tag: string expected",		C_result);
		XL_GETOBJID( XL_correlMatId,		C_correlMatId,		" ARM_ERR: correl Mat id: object expected",		C_result);

		/// use the concept of Functor to transfer the knowledge of
		/// a function with a context
		correlManagerFromMatFunc ourFunc( C_mktTag, C_intraMktTag, LocalGetNumObjectId( C_correlMatId ) );
		
		/// call the general function
		fillXL_Result( LOCAL_CORRELMANAGER_CLASS, ourFunc, C_result, XL_result, PersistentInXL );

	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_CorrelManagerFromMatCommon" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}



//////////////////////////////////////////////////////
/// Addin to create a correlation manager from a correl mat
//////////////////////////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_CorrelManagerFromMat_Create(
    LPXLOPER XL_mktTag,
    LPXLOPER XL_intraMktTag,
    LPXLOPER XL_correlMatId )
{
	ADD_LOG("Local_CorrelManagerFromMat_Create");
	bool PersistentInXL = true;
	return Local_CorrelManagerFromMatCommon(
		XL_mktTag,
		XL_intraMktTag,
		XL_correlMatId,
		PersistentInXL );
}


//////////////////////////////////////////////////////
/// Addin to create a correlation manager from a correl mat
/// Version for VBA 
//////////////////////////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_PXL_CorrelManagerFromMat_Create(
    LPXLOPER XL_mktTag,
    LPXLOPER XL_intraMktTag,
    LPXLOPER XL_correlMatId )
{
	ADD_LOG("Local_PXL_CorrelManagerFromMat_Create");
	bool PersistentInXL = false;
	return Local_CorrelManagerFromMatCommon(
		XL_mktTag,
		XL_intraMktTag,
		XL_correlMatId,
		PersistentInXL );
}


	
////////////////////////////////////////////////////////////////
/// Addin to create a correlation manager from a correl curve///
/// parameters are three Vectors for create the correlManager////
////////////////////////////////////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_CreateGenCorrelManager( LPXLOPER XL_mktTags,
																	LPXLOPER XL_intraMktTags,
																	LPXLOPER XL_correlCurveIds)
{
	ADD_LOG("Local_CreateGenCorrelManager");
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
		VECTOR<CCString> C_intraMktTags;
		VECTOR<CCString> C_correlVolIds;

		XL_readStrVector( XL_mktTags, C_mktTags," ARM_ERR: mkt Tag: array expected",DOUBLE_TYPE, C_result);
		XL_readStrVector( XL_intraMktTags,C_intraMktTags," ARM_ERR: intra Mkt Tag:  array expected",DOUBLE_TYPE,	C_result);
		XL_readStrVector( XL_correlCurveIds,C_correlVolIds," ARM_ERR: correl id: array expected",DOUBLE_TYPE,	C_result);

		int i;
		long sizeCurve = C_correlVolIds.size();
		vector<long> vCorrelVolIds(sizeCurve);
		for( i=0; i<sizeCurve; i++ )
		{
			vCorrelVolIds[i] = ( strcmp(C_correlVolIds[i], "NULL") == 0 )
								? ARM_NULL_OBJECT : LocalGetNumObjectId ( C_correlVolIds[i] );
		}

		long retCode;
		long objId;
		CCString prevClass;
		
		CCString curClass = LOCAL_CORRELMANAGER_CLASS;
		CCString stringId = GetLastCurCellEnvValue ();

		if(!stringId)
		{
			retCode =  ARMLOCAL_CreateGenCorrelManager(	C_mktTags,
														C_intraMktTags,
														vCorrelVolIds, 
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
				retCode = ARMLOCAL_CreateGenCorrelManager (	C_mktTags,
															C_intraMktTags,
															vCorrelVolIds, 
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
				retCode = ARMLOCAL_CreateGenCorrelManager (	C_mktTags,
															C_intraMktTags,
															vCorrelVolIds, 
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_CreateGenCorrelManager" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


////////////////////////////////////////////////////////////////
/// Addin to create a correlation manager from a correl curve///
/// parameters are three Vectors for create the correlManager////
////////////////////////////////////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_PXL_CreateGenCorrelManager( LPXLOPER XL_mktTags,
																	    LPXLOPER XL_intraMktTags,
																		LPXLOPER XL_correlCurveIds)
{
	ADD_LOG("Local_PXL_CreateGenCorrelManager");
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
		VECTOR<CCString> C_intraMktTags;
		VECTOR<CCString> C_correlVolIds;

		XL_readStrVector( XL_mktTags, C_mktTags," ARM_ERR: mkt Tag: array expected",DOUBLE_TYPE, C_result);
		XL_readStrVector( XL_intraMktTags,C_intraMktTags," ARM_ERR: intra Mkt Tag:  array expected",DOUBLE_TYPE,	C_result);
		XL_readStrVector( XL_correlCurveIds,C_correlVolIds," ARM_ERR: correl id: array expected",DOUBLE_TYPE,	C_result);

		int i;
		long sizeCurve = C_correlVolIds.size();
		vector<long> vCorrelVolIds(sizeCurve);
		for( i=0; i<sizeCurve; i++ )
		{
			vCorrelVolIds[i] = ( strcmp(C_correlVolIds[i], "NULL") == 0 )
								? ARM_NULL_OBJECT : LocalGetNumObjectId ( C_correlVolIds[i] );
		}

		long retCode;
		long objId;
		CCString prevClass;
		
		CCString curClass = LOCAL_CORRELMANAGER_CLASS;
		CCString stringId;

		retCode =  ARMLOCAL_CreateGenCorrelManager(	C_mktTags,
													C_intraMktTags,
													vCorrelVolIds, 
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_PXL_CreateGenCorrelManager" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}



/////////////////////////////////////////
/// Functor to fill a correlation manager
/////////////////////////////////////////
class correlManagerFillFunc: public ARMResultLong2LongFunc
{
public:
	correlManagerFillFunc( const CCString& mktTag, const CCString& intraMktTag, 
		const VECTOR<CCString>& X,
		const VECTOR<CCString>& Y, const VECTOR<double>& Z,
		long correlManagerId )
	:	C_mktTag( mktTag ), C_intraMktTag( intraMktTag ), 
		C_X( X ), C_Y( Y ), C_Z( Z ),
		C_correlManagerId( correlManagerId )
	{};
	long operator()( ARM_result& result, long objId )
	{
		return ARMLOCAL_CorrelManager_Fill( C_mktTag, C_intraMktTag,
			C_X, C_Y, C_Z, C_correlManagerId, result, objId );
	}

private:
	CCString C_mktTag;
	CCString C_intraMktTag;
	double	C_asOfDate;
	VECTOR<CCString> C_X;
	VECTOR<CCString> C_Y;
	VECTOR<double> C_Z;
	long C_correlManagerId;
};


/////////////////////////////////////////
/// Correlation manager fill common function
/////////////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_CorrelManagerFillCommon(
    LPXLOPER XL_mktTag,
    LPXLOPER XL_intraMktTag,
    LPXLOPER XL_X,
    LPXLOPER XL_Y,
    LPXLOPER XL_Z,
    LPXLOPER XL_correlManagerId,
	bool PersistentInXL )
{
	ADD_LOG("Local_CorrelManagerFillCommon");
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

		CCString C_mktTag;
		CCString C_intraMktTag;
		VECTOR<CCString> C_X;
		VECTOR<CCString> C_Y;
		VECTOR<double>	 C_Z;
		long C_correlManagerId;

		XL_readStrCell(  XL_mktTag,			C_mktTag,			" ARM_ERR: mkt Tag: string expected",			C_result);
		XL_readStrCell(	 XL_intraMktTag,	C_intraMktTag,		" ARM_ERR: intra Mkt Tag: string expected",		C_result);
		XL_readStrVector(XL_X,				C_X,				" ARM_ERR: X: array of numeric/maturity expected",	DOUBLE_TYPE, C_result);
		XL_readStrVector(XL_Y,				C_Y,				" ARM_ERR: Y: array of numeric/maturity expected",	DOUBLE_TYPE, C_result);
		XL_readNumVector(XL_Z,				C_Z,				" ARM_ERR: Z: matrix of numeric expected",C_result);

		XL_GETOBJID( XL_correlManagerId,	C_correlManagerId,	" ARM_ERR: correl manager: Object expected",	C_result);

		/// use the concept of Functor to transfer the knowledge of
		/// a function with a context
		correlManagerFillFunc ourFunc( C_mktTag, C_intraMktTag, C_X, C_Y, C_Z, C_correlManagerId );
		
		/// call the general function
		fillXL_Result( LOCAL_CORRELMANAGER_CLASS, ourFunc, C_result, XL_result, PersistentInXL );

	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_CorrelManagerFillCommon" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}



//////////////////////////////////////////////////////
/// Addin to fill a correlation manager
//////////////////////////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_CorrelManager_Fill(
    LPXLOPER XL_mktTag,
    LPXLOPER XL_intraMktTag,
    LPXLOPER XL_X,
    LPXLOPER XL_Y,
    LPXLOPER XL_Z,
    LPXLOPER XL_correlManagerId )
{
	ADD_LOG("Local_CorrelManager_Fill");
	bool PersistentInXL = true;
	return Local_CorrelManagerFillCommon(
		XL_mktTag,
		XL_intraMktTag,
		XL_X,
		XL_Y,
		XL_Z,
		XL_correlManagerId,
		PersistentInXL );
}



//////////////////////////////////////////////////////
/// Addin to fill a correlation manager
/// Version for VBA
//////////////////////////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_PXL_CorrelManager_Fill(
    LPXLOPER XL_mktTag,
    LPXLOPER XL_intraMktTag,
    LPXLOPER XL_X,
    LPXLOPER XL_Y,
    LPXLOPER XL_Z,
    LPXLOPER XL_correlManagerId )
{
	ADD_LOG("Local_PXL_CorrelManager_Fill");
	bool PersistentInXL = false;
	return Local_CorrelManagerFillCommon(
		XL_mktTag,
		XL_intraMktTag,
		XL_X,
		XL_Y,
		XL_Z,
		XL_correlManagerId,
		PersistentInXL );
}





/////////////////////////////////////////
/// Functor to fill a correlation manager
/////////////////////////////////////////
class correlManagerFillFromMatFunc: public ARMResultLong2LongFunc
{
public:
	correlManagerFillFromMatFunc( const CCString& mktTag, const CCString& intraMktTag, 
	long correlMatId, long correlManagerId )
	:	C_mktTag( mktTag ), C_intraMktTag( intraMktTag ),  C_correlMatId( correlMatId ),
		C_correlManagerId( correlManagerId )
	{};
	long operator()( ARM_result& result, long objId )
	{
		return ARMLOCAL_CorrelManager_FillFromMat( C_mktTag, C_intraMktTag,
			C_correlMatId, C_correlManagerId, result, objId );
	}

private:
	CCString C_mktTag;
	CCString C_intraMktTag;
	long C_correlMatId;
	long C_correlManagerId;
};


/////////////////////////////////////////
/// Correlation manager fill common function
/////////////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_CorrelManagerFillCommon(
    LPXLOPER XL_mktTag,
    LPXLOPER XL_intraMktTag,
    LPXLOPER XL_correlMatId,
	LPXLOPER XL_correlManagerId,
	bool PersistentInXL )
{
	ADD_LOG("Local_CorrelManagerFillCommon");
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

		CCString C_mktTag;
		CCString C_intraMktTag;
		VECTOR<CCString> C_X;
		VECTOR<CCString> C_Y;
		VECTOR<double>	 C_Z;
		long C_correlMatId;
		long C_correlManagerId;

		XL_readStrCell(  XL_mktTag,			C_mktTag,			" ARM_ERR: mkt Tag: string expected",			C_result);
		XL_readStrCell(	 XL_intraMktTag,	C_intraMktTag,		" ARM_ERR: intra Mkt Tag: string expected",		C_result);
		XL_GETOBJID( XL_correlMatId,		C_correlMatId,		" ARM_ERR: correl Mat id: object expected",		C_result);
		XL_GETOBJID( XL_correlManagerId,	C_correlManagerId,	" ARM_ERR: correl manager: Object expected",	C_result);


		/// use the concept of Functor to transfer the knowledge of
		/// a function with a context
		correlManagerFillFromMatFunc ourFunc( C_mktTag, C_intraMktTag, C_correlMatId, C_correlManagerId );
		
		/// call the general function
		fillXL_Result( LOCAL_CORRELMANAGER_CLASS, ourFunc, C_result, XL_result, PersistentInXL );

	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_CorrelManagerFillCommon" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}



//////////////////////////////////////////////////////
/// Addin to fill a correlation manager
//////////////////////////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_CorrelManagerFromMat_Fill(
    LPXLOPER XL_mktTag,
    LPXLOPER XL_intraMktTag,
    LPXLOPER XL_correlMatId,
    LPXLOPER XL_correlManagerId )
{
	ADD_LOG("Local_CorrelManagerFromMat_Fill");
	bool PersistentInXL = true;
	return Local_CorrelManagerFillCommon(
		XL_mktTag,
		XL_intraMktTag,
		XL_correlMatId,
		XL_correlManagerId,
		PersistentInXL );
}



//////////////////////////////////////////////////////
/// Addin to fill a correlation manager
/// Version for VBA
//////////////////////////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_PXL_CorrelManagerFromMat_Fill(
    LPXLOPER XL_mktTag,
    LPXLOPER XL_intraMktTag,
    LPXLOPER XL_correlMatId,
    LPXLOPER XL_correlManagerId )
{
	ADD_LOG("Local_PXL_CorrelManagerFromMat_Fill");
	bool PersistentInXL = false;
	return Local_CorrelManagerFillCommon(
		XL_mktTag,
		XL_intraMktTag,
		XL_correlMatId,
		XL_correlManagerId,
		PersistentInXL );
}





/////////////////////////////////////////
/// Correlation manager fill common function
/////////////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_ComputeCorrelFromManagerCommon(
    LPXLOPER XL_mktTag,
    LPXLOPER XL_intraMktTag,
    LPXLOPER XL_X,
    LPXLOPER XL_Y,
    LPXLOPER XL_correlManagerId,
	bool PersistentInXL )
{
	ADD_LOG("Local_ComputeCorrelFromManagerCommon");
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

		CCString C_mktTag;
		CCString C_intraMktTag;
		double C_X;
		double C_Y;
		long C_correlManagerId;

		XL_readStrCell(  XL_mktTag,			C_mktTag,			" ARM_ERR: mkt Tag: string expected",			C_result);
		XL_readStrCell(	 XL_intraMktTag,	C_intraMktTag,		" ARM_ERR: intra Mkt Tag: string expected",		C_result);
		XL_readNumCell(	 XL_X,				C_X,				" ARM_ERR: X: double expected",					C_result);
		XL_readNumCell(	 XL_Y,				C_Y,				" ARM_ERR: Y: double expected",					C_result);
		XL_GETOBJID( XL_correlManagerId,	C_correlManagerId,	" ARM_ERR: correl manager: Object expected",	C_result);

		/// use the concept of Functor to transfer the knowledge of
		/// a function with a context
		long retCode;
		if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 1)
			retCode = ARMLOCAL_ComputeCorrelFromCorrelManager( C_mktTag, C_intraMktTag,
				C_X, C_Y, C_correlManagerId, C_result );
		else
			retCode = ARM_KO;
		
		/// feed the LPXLOPER object result 
		if (retCode == ARM_OK)
		{
			FreeCurCellErr ();
			XL_result.xltype = xltypeNum;
			XL_result.val.num = C_result.getDouble();
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_ComputeCorrelFromManagerCommon" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}



//////////////////////////////////////////////////////
/// Addin to compute a correlation from a correlation manager
//////////////////////////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_ComputeCorrelFromManager(
    LPXLOPER XL_mktTag,
    LPXLOPER XL_intraMktTag,
    LPXLOPER XL_X,
    LPXLOPER XL_Y,
    LPXLOPER XL_correlManagerId )
{
	ADD_LOG("Local_ComputeCorrelFromManager");
	bool PersistentInXL = true;
	return Local_ComputeCorrelFromManagerCommon(
		XL_mktTag,
		XL_intraMktTag,
		XL_X,
		XL_Y,
		XL_correlManagerId,
		PersistentInXL );
}



////////////////////////////////////////////
/// addin to compute implicit correlations 
/// from volatilities
////////////////////////////////////////////
LPXLOPER WINAPI Local_Vol_to_Cor(
		LPXLOPER XL_ZCVolV, 
		LPXLOPER XL_YtYVolV, 
		LPXLOPER XL_MaturityV )
{
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
		/// C variable to get the variables from LPXLOPER 
		/// Variables
		VECTOR<double> C_ZCVolV;
		XL_readNumVector(XL_ZCVolV,		C_ZCVolV,		" ARM_ERR: Zero coupon volatilities: array of numeric expected",C_result);
		VECTOR<double> C_YtYVolV; 
		XL_readNumVector(XL_YtYVolV,	 C_YtYVolV,		" ARM_ERR: Year to year volatilities: array of numeric expected",C_result);
		VECTOR<double> C_MaturityV;
		XL_readNumVector(XL_MaturityV,	 C_MaturityV,	" ARM_ERR: Maturities: array of numeric expected",C_result);


		/// output 
		VECTOR<double> C_CorV;
		
		/// call the function with nomore LPXLOPER objects
		// this fction store the result into C_result
		long retCode = ARMLOCAL_Vol_to_Cor( 
				C_ZCVolV, 
				C_YtYVolV, 
				C_MaturityV,
				C_CorV,
				C_result );

		/// feed the LPXLOPER object result 
		if (retCode == ARM_OK)
		{
			FreeCurCellErr ();
			XL_writeNumVector( XL_result, C_CorV, "ARM ERR: could not write the output", C_result );
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_Vol_to_Cor" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


////////////////////////////////////////////
/// addin to compute zero coupon volatilities
/// from implicit correlations and year to year
/// volatilities
////////////////////////////////////////////
LPXLOPER WINAPI Local_YtYCor_to_ZC(
		LPXLOPER XL_YtYVolV, 
		LPXLOPER XL_CorV,
		LPXLOPER XL_MaturityV )
{
	// this is defined first because it is used in XL macros
	static XLOPER XL_result;
	ARM_result C_result;
	
	// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
		/// to avoid computation if called by the wizard
		ARM_NOCALCIFWIZ();
		
		// error
		static int error;
		static char* reason = "";

		/// C variable to get the variables from LPXLOPER 
		/// Variables
		VECTOR<double> C_CorV;
		XL_readNumVector(XL_CorV,		C_CorV,			" ARM_ERR: Zero coupon volatilities: array of numeric expected",C_result);
		VECTOR<double> C_YtYVolV; 
		XL_readNumVector(XL_YtYVolV,	C_YtYVolV,		" ARM_ERR: Year to year volatilities: array of numeric expected",C_result);
		VECTOR<double> C_MaturityV;
		XL_readNumVector(XL_MaturityV,	C_MaturityV,	" ARM_ERR: Maturities: array of numeric expected",C_result);

		/// output 
		VECTOR<double> C_ZCVolV;
		
		/// call the function with nomore LPXLOPER objects
		// this fction store the result into C_result
		long retCode = ARMLOCAL_YtYCor_to_ZC( 
				C_YtYVolV,
				C_CorV,  
				C_MaturityV,
				C_ZCVolV,
				C_result );

		/// feed the LPXLOPER object result 
		if (retCode == ARM_OK)
		{
			FreeCurCellErr ();
			XL_writeNumVector( XL_result, C_ZCVolV, "ARM ERR: could not write the output", C_result );
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_YtYCor_to_ZC" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


////////////////////////////////////////////
/// addin to compute year to year volatilities
/// from implicit correlations and zero coupon
/// volatilities
////////////////////////////////////////////
LPXLOPER WINAPI Local_ZCCor_to_YtY(
		LPXLOPER XL_ZCVolV, 
		LPXLOPER XL_CorV,
		LPXLOPER XL_MaturityV )
{
	// this is defined first because it is used in XL macros
	static XLOPER XL_result;
	ARM_result C_result;
	
	// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
		/// to avoid computation if called by the wizard
		ARM_NOCALCIFWIZ();
		
		// error
		static int error;
		static char* reason = "";

		/// C variable to get the variables from LPXLOPER 
		/// Variables
		VECTOR<double> C_CorV;
		XL_readNumVector(XL_CorV,		C_CorV,			" ARM_ERR: Zero coupon volatilities: array of numeric expected",	C_result);
		VECTOR<double> C_ZCVolV; 
		XL_readNumVector(XL_ZCVolV,		C_ZCVolV,		" ARM_ERR: Year to year volatilities: array of numeric expected",	C_result);
		VECTOR<double> C_MaturityV;
		XL_readNumVector(XL_MaturityV,	C_MaturityV,	" ARM_ERR: Maturities: array of numeric expected",					C_result);


		/// output 
		VECTOR<double> C_YtYVolV;
		
		/// call the function with nomore LPXLOPER objects
		// this fction store the result into C_result
		long retCode = ARMLOCAL_ZCCor_to_YtY( 
				C_ZCVolV,
				C_CorV,  
				C_MaturityV,
				C_YtYVolV,
				C_result );

		/// feed the LPXLOPER object result 
		if (retCode == ARM_OK)
		{
			FreeCurCellErr ();
			XL_writeNumVector( XL_result, C_YtYVolV, "ARM ERR: could not write the output", C_result );
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_ZCCor_to_YtY" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}

////////////////////////////////////////////
/// addin to compute year to year volatilities
/// from implicit correlations and zero coupon
/// volatilities
////////////////////////////////////////////
LPXLOPER WINAPI Local_Bounds(
		LPXLOPER XL_ZCVolV, 
		LPXLOPER XL_YtYVolV, 
		LPXLOPER XL_CorV,
		LPXLOPER XL_MaturityV,
		LPXLOPER XL_Choice )

{
	// this is defined first because it is used in XL macros
	static XLOPER XL_result;
	LPXLOPER pxArray;
	ARM_result C_result;
	
	// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
		/// to avoid computation if called by the wizard
		ARM_NOCALCIFWIZ();
		
		// error
		static int error;
		static char* reason = "";

		/// C variable to get the variables from LPXLOPER 
		/// Variables
		VECTOR<double> C_CorV;
		XL_readNumVector(XL_CorV,		C_CorV,			" ARM_ERR: Zero coupon volatilities: array of numeric expected",	C_result);
		VECTOR<double> C_ZCVolV; 
		XL_readNumVector(XL_ZCVolV,		C_ZCVolV,		" ARM_ERR: Year to year volatilities: array of numeric expected",	C_result);
		VECTOR<double> C_YtYVolV; 
		XL_readNumVector(XL_YtYVolV,	C_YtYVolV,		" ARM_ERR: Year to year volatilities: array of numeric expected",	C_result);
		VECTOR<double> C_MaturityV;
		XL_readNumVector(XL_MaturityV,	C_MaturityV,	" ARM_ERR: Maturities: array of numeric expected",					C_result);
		CCString C_Choice;
		XL_readStrCell(XL_Choice,		C_Choice,		" ARM_ERR: Maturities: array of numeric expected",					C_result);

		/// output 
		VECTOR<double> C_UBound;
		VECTOR<double> C_LBound;
		CCString C_TBound;
		string C_internTBound( " " );
		string C_internChoice;

		
		C_internChoice = CCSTringToSTLString( C_Choice );
		
		
		/// call the function with nomore LPXLOPER objects
		// this fction store the result into C_result
		long retCode = ARMLOCAL_Bounds( 
				C_ZCVolV,
				C_YtYVolV,
				C_CorV,  
				C_MaturityV,
				C_UBound,
				C_LBound,
				C_internChoice,
				C_internTBound,
				C_result );
		
		if (C_internTBound=="internal") {C_TBound=(char*)" internal ";};	
		if (C_internTBound=="external") {C_TBound=(char*)" external ";};	
		

		/// feed the LPXLOPER object result 
		if (retCode == ARM_OK)
		{
			/// free previous allocation!
			/// if any
			if( XL_result.xltype == xltypeMulti )
				HGLOBAL ReturnPtr = GlobalFree( (HGLOBAL) XL_result.val.array.lparray );

			long nbrows = C_result.getLong();
			int nbcolumns;
			nbcolumns = 3;
			
			FreeCurCellErr ();
			XL_result.xltype = xltypeMulti;
			XL_result.val.array.columns = nbcolumns;
			XL_result.val.array.rows = nbrows; 


			/// allocation
			XL_result.val.array.lparray = pxArray = (LPXLOPER)GlobalAlloc (GMEM_ZEROINIT, nbrows * nbcolumns * sizeof (XLOPER));

			for (int i = 0; i < nbrows/3; i++)
			{	
				for (int j = 0; j < nbcolumns; j++)
				{
					if (j == 2)
					{
						CCString internalStr = "internal";
						CCString externalStr = "external";

						pxArray[XL_Coordonnate2Rank (i, j, nbcolumns)].xltype = xltypeStr;
						if (C_result.getArray(i+j*nbrows/3)==1)
						{

							pxArray[XL_Coordonnate2Rank (i, j, nbcolumns)].val.str = XL_StrC2StrPascal (internalStr);
						}
						else
						{
							pxArray[XL_Coordonnate2Rank (i, j, nbcolumns)].val.str = XL_StrC2StrPascal (externalStr);
						}
					}
					
					else
					{
						pxArray[XL_Coordonnate2Rank (i, j, nbcolumns)].xltype = xltypeNum;
						pxArray[XL_Coordonnate2Rank (i, j, nbcolumns)].val.num = C_result.getArray(i+j*nbrows/3);
					}

				}
			}
		
		
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_Bounds" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}

////////////////////////////////////////////
/// addin to compute implicit correlations 
/// from volatilities in the homogenuous case
////////////////////////////////////////////
LPXLOPER WINAPI Local_HmgVol_to_Cor(
		LPXLOPER XL_ZCVolV, 
		LPXLOPER XL_YtYVolV, 
		LPXLOPER XL_MaturityV,
		LPXLOPER XL_length)
{
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
		/// C variable to get the variables from LPXLOPER 
		/// Variables
		VECTOR<double> C_ZCVolV;
		XL_readNumVector(XL_ZCVolV,		C_ZCVolV,		" ARM_ERR: Zero coupon volatilities: array of numeric expected",	C_result);
		VECTOR<double> C_YtYVolV; 
		XL_readNumVector(XL_YtYVolV,	C_YtYVolV,		" ARM_ERR: Year to year volatilities: array of numeric expected",	C_result);
		VECTOR<double> C_MaturityV;
		XL_readNumVector(XL_MaturityV,	C_MaturityV,	" ARM_ERR: Maturities: array of numeric expected",					C_result);
		double C_length;
		XL_readNumCell(XL_length,		C_length,		" ARM_ERR: lenght: for fixing step of interpolation ",				C_result);


		/// output 
		VECTOR<double> C_CorV;
		
		/// call the function with nomore LPXLOPER objects
		// this fction store the result into C_result
		long retCode = ARMLOCAL_HmgVol_to_Cor( 
				C_ZCVolV, 
				C_YtYVolV, 
				C_MaturityV,
				C_CorV,
				C_length,
				C_result );

		/// feed the LPXLOPER object result 
		if (retCode == ARM_OK)
		{
			FreeCurCellErr ();
			XL_writeNumVector( XL_result, C_CorV, "ARM ERR: could not write the output", C_result );
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_HmgVol_to_Cor" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


LPXLOPER WINAPI Local_HmgZCCor_to_YtY(
		LPXLOPER XL_ZCVolV, 
		LPXLOPER XL_CorV,
		LPXLOPER XL_MaturityV,
		LPXLOPER XL_length )
{
	// this is defined first because it is used in XL macros
	static XLOPER XL_result;
	ARM_result C_result;
	
	// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
		/// to avoid computation if called by the wizard
		ARM_NOCALCIFWIZ();
		
		// error
		static int error;
		static char* reason = "";

		/// C variable to get the variables from LPXLOPER 
		/// Variables
		VECTOR<double> C_CorV;
		XL_readNumVector(XL_CorV,		C_CorV,			" ARM_ERR: Zero coupon volatilities: array of numeric expected",	C_result);
		VECTOR<double> C_ZCVolV; 
		XL_readNumVector(XL_ZCVolV,		C_ZCVolV,		" ARM_ERR: Year to year volatilities: array of numeric expected",	C_result);
		VECTOR<double> C_MaturityV;
		XL_readNumVector(XL_MaturityV,	C_MaturityV,	" ARM_ERR: Maturities: array of numeric expected",					C_result);
		double C_length;
		XL_readNumCell(XL_length,		C_length,		" ARM_ERR: lenght: for fixing step of interpolation ",				C_result);


		/// output 
		VECTOR<double> C_YtYVolV;
		
		/// call the function with nomore LPXLOPER objects
		// this fction store the result into C_result
		long retCode = ARMLOCAL_HmgZCCor_to_YtY( 
				C_ZCVolV,
				C_CorV,  
				C_MaturityV,
				C_length,
				C_YtYVolV,
				C_result );

		/// feed the LPXLOPER object result 
		if (retCode == ARM_OK)
		{
			FreeCurCellErr ();
			XL_writeNumVector( XL_result, C_YtYVolV, "ARM ERR: could not write the output", C_result );
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_HmgZCCor_to_YtY" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


////////////////////////////////////////////
/// addin to compute zero coupon volatilities
/// from implicit correlations and year to year
/// volatilities in the homogeneous case
////////////////////////////////////////////
LPXLOPER WINAPI Local_HmgYtYCor_to_ZC(
		LPXLOPER XL_YtYVolV, 
		LPXLOPER XL_CorV,
		LPXLOPER XL_MaturityV,
		LPXLOPER XL_length )
{
	// this is defined first because it is used in XL macros
	static XLOPER XL_result;
	ARM_result C_result;
	
	// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
		/// to avoid computation if called by the wizard
		ARM_NOCALCIFWIZ();
		
		// error
		static int error;
		static char* reason = "";

		/// C variable to get the variables from LPXLOPER 
		/// Variables
		VECTOR<double> C_CorV;
		XL_readNumVector(XL_CorV,		C_CorV,			" ARM_ERR: Zero coupon volatilities: array of numeric expected",C_result);
		VECTOR<double> C_YtYVolV; 
		XL_readNumVector(XL_YtYVolV,	C_YtYVolV,		" ARM_ERR: Year to year volatilities: array of numeric expected",C_result);
		VECTOR<double> C_MaturityV;
		XL_readNumVector(XL_MaturityV,	C_MaturityV,	" ARM_ERR: Maturities: array of numeric expected",				C_result);
		double C_length;
		XL_readNumCell(XL_length,		C_length,		" ARM_ERR: lenght: for fixing step of interpolation ",			C_result);

		/// output 
		VECTOR<double> C_ZCVolV;
		
		/// call the function with nomore LPXLOPER objects
		// this fction store the result into C_result
		long retCode = ARMLOCAL_HmgYtYCor_to_ZC( 
				C_YtYVolV,
				C_CorV,  
				C_MaturityV,
				C_length,
				C_ZCVolV,
				C_result );

		/// feed the LPXLOPER object result 
		if (retCode == ARM_OK)
		{
			FreeCurCellErr ();
			XL_writeNumVector( XL_result, C_ZCVolV, "ARM ERR: could not write the output", C_result );
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_HmgYtYCor_to_ZC" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}



/////////////////////////////////////////
/// Function to compute the inflation swap's 
/// implied volatility from the year-on-year ones
/////////////////////////////////////////

LPXLOPER WINAPI Local_VolYoY_to_VolSwp(
		LPXLOPER XL_DFactor, 
		LPXLOPER XL_FwdCPI, 
		LPXLOPER XL_Vol_DF,
		LPXLOPER XL_Vol_YoY, 
		LPXLOPER XL_AvgCor, 
		LPXLOPER XL_Dates, 
		LPXLOPER XL_Tenors, 
		LPXLOPER XL_SwpRate)
{
	// this is defined first because it is used in XL macros
	static XLOPER XL_result;
	ARM_result C_result;
	
	// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
		/// to avoid computation if called by the wizard
		ARM_NOCALCIFWIZ();
		
		// error
		static int error;
		static char* reason = "";
		/// C variable to get the variables from LPXLOPER 
		/// Variables

		VECTOR<double> C_DFactor;
		XL_readNumVector(XL_DFactor,	C_DFactor,	" ARM_ERR: DFactor: array of numeric expected",	C_result);

		VECTOR<double> C_FwdCPI; 
		XL_readNumVector(XL_FwdCPI,		C_FwdCPI,	" ARM_ERR: FwdCPI: array of numeric expected",	C_result);

		VECTOR<double> C_Vol_DF;
		XL_readNumVector(XL_Vol_DF,		C_Vol_DF,	" ARM_ERR: Vol_DF: array of numeric expected",	C_result);

		VECTOR<double> C_Vol_YoY;
		XL_readNumVector(XL_Vol_YoY,	C_Vol_YoY,	" ARM_ERR: Vol_YoY: array of numeric expected",	C_result);

		VECTOR<double> C_AvgCor;
		XL_readNumVector(XL_AvgCor,		C_AvgCor,	" ARM_ERR: AvgCor: array of numeric expected",	C_result);

		VECTOR<double> C_Dates;
		XL_readNumVector(XL_Dates,		C_Dates,	" ARM_ERR: Dates: array of numeric expected",	C_result);

		VECTOR<double> C_Tenors;
		XL_readNumVector(XL_Tenors,		C_Tenors,	" ARM_ERR: Tenors: array of numeric expected",	C_result);

		double C_SwpRate;
		XL_readNumCell(XL_SwpRate,		C_SwpRate,	" ARM_ERR: SwpRate: the swap rate",				C_result);

		/// output 
		double C_Vol_Swp = 0.;
		
		/// call the function with nomore LPXLOPER objects
		// this fction store the result into C_result
		long retCode = ARMLOCAL_VolYoY_to_VolSwp(
			C_DFactor, 
			C_FwdCPI, 
			C_Vol_DF,
			C_Vol_YoY, 
			C_AvgCor, 
			C_Dates, 
			C_Tenors, 
			C_SwpRate, 
			C_Vol_Swp,
			C_result );

		/// feed the LPXLOPER object result 
		if (retCode == ARM_OK)
		{
			FreeCurCellErr ();
			XL_result.xltype = xltypeNum;
			XL_result.val.num = C_Vol_Swp;//C_result.getDouble();
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_VolYoY_to_VolSwp" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;



}


/////////////////////////////////////////
/// Functor to create a seasonality manager
/////////////////////////////////////////
class seasonalityManagerFunc: public ARMResultLong2LongFunc
{
public:
	seasonalityManagerFunc( const VECTOR<CCString>& monthList, 
		const VECTOR<double>& seasonSpreadList, 
		const VECTOR<double>& seasonHorizonList,
		long seasonAdjMode )
		:	C_monthList( monthList ), C_seasonSpreadList( seasonSpreadList ), C_seasonHorizonList(seasonHorizonList), C_seasonAdjMode( seasonAdjMode )
		
		
	{};
	long operator()( ARM_result& result, long objId )
	{
		return ARMLOCAL_SeasonalityManager_Create(
			C_monthList,
			C_seasonSpreadList, 
			C_seasonHorizonList,
			C_seasonAdjMode,
			result,
			objId);
	}
	
private:
	VECTOR<CCString> C_monthList;
	VECTOR<double> C_seasonSpreadList;
	VECTOR<double> C_seasonHorizonList;
	long C_seasonAdjMode;
};


/////////////////////////////////////////
/// Correlation manager create common function
/////////////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_SeasonalityManagerCommon(
																	 LPXLOPER XL_monthList,
																	 LPXLOPER XL_seasonSpreadList,
																	 LPXLOPER XL_seasonHorizonList,
																	 LPXLOPER XL_seasonAdjMode,
																	 bool PersistentInXL )
{
	ADD_LOG("Local_SeasonalityManagerCommon");
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
		
		VECTOR<CCString> C_monthList;
		XL_readStrVector(XL_monthList,			C_monthList, " ARM_ERR: month list: string expected", DOUBLE_TYPE, C_result);
		
		VECTOR<double>	 C_seasonSpreadList;
		XL_readNumVector(XL_seasonSpreadList,	C_seasonSpreadList,	" ARM_ERR: X: array of numeric", C_result);
		
		//LPXLOPER XL_seasonHorizonList,

		VECTOR<double>	 C_seasonHorizonList;
		VECTOR<double>	 C_seasonHorizonListDef;
		XL_readNumVectorWD(XL_seasonHorizonList, C_seasonHorizonList, C_seasonHorizonListDef, " ARM_ERR: Y: array of numeric", C_result);

		CCString C_seasonAdjModeStr;
		char defaultValue[] = "PLUS";
		XL_readStrCellWD(XL_seasonAdjMode, C_seasonAdjModeStr, defaultValue, " ARM_ERR: seasonality correction mode : string expected", C_result);
		long C_seasonAdjMode = ARM_ConvSeasonalityMode( C_seasonAdjModeStr, C_result );
		
		/// use the concept of Functor to transfer the knowledge of
		/// a function with a context
		seasonalityManagerFunc ourFunc( C_monthList, C_seasonSpreadList, C_seasonHorizonList, C_seasonAdjMode );
		
		/// call the general function
		fillXL_Result( LOCAL_SEASONMANAGER_CLASS, ourFunc, C_result, XL_result, PersistentInXL );
		
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END
		
		/// to catch arm exception
		ARM_XL_CATCH_ARM_EXPT
		
		/// to cath all the other exceptions
		ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_SeasonalityManagerCommon" )
		
		/// return the result as an LPXLOPER
		return (LPXLOPER)&XL_result;
}



//////////////////////////////////////////////////////
/// Addin to create a seasonality manager
//////////////////////////////////////////////////////

__declspec(dllexport) LPXLOPER WINAPI Local_SeasonalityManager_Create(
	LPXLOPER XL_monthList,
	LPXLOPER XL_seasonSpreadList,
	LPXLOPER XL_seasonHorizonList,
	LPXLOPER XL_seasonAdjMode )
{
	ADD_LOG("Local_SeasonalityManager_Create");
	bool PersistentInXL = true;
	return Local_SeasonalityManagerCommon(
		XL_monthList,
		XL_seasonSpreadList,
		XL_seasonHorizonList,
		XL_seasonAdjMode,
		PersistentInXL);
}



//////////////////////////////////////////////////////
/// Addin to create a seasonality manager
/// Version for VBA
//////////////////////////////////////////////////////

__declspec(dllexport) LPXLOPER WINAPI Local_PXL_SeasonalityManager_Create(
	LPXLOPER XL_monthList,
	LPXLOPER XL_seasonSpreadList,
	LPXLOPER XL_seasonHorizonList,
	LPXLOPER XL_seasonAdjMode )
{
	ADD_LOG("Local_PXL_SeasonalityManager_Create");
	bool PersistentInXL = false;
	return Local_SeasonalityManagerCommon(
		XL_monthList,
		XL_seasonSpreadList,
		XL_seasonHorizonList,
		XL_seasonAdjMode,
		PersistentInXL);
}


/////////////////////////////////////////
/// Functor to create an inflation curve from Summit
/////////////////////////////////////////
class getInfZcFromSummitFunc: public ARMResultLong2LongFunc
{
public:
	getInfZcFromSummitFunc( const CCString& index,
		const CCString& ccy,
		const CCString& cvname,
		double date,
		long seasonAdjId,
		long seasonAdjModeId)
	:	C_index( index ), C_ccy( ccy ), C_cvname( cvname ), C_date( date ), C_seasonAdjId(seasonAdjId), C_seasonAdjModeId(seasonAdjModeId)
	{};
	long operator()( ARM_result& result, long objId )
	{
		return ARMLOCAL_GetInfZcFromSummit_Create( C_index, C_ccy, C_cvname, C_date, C_seasonAdjId, C_seasonAdjModeId, result, objId);
	}

private:
	CCString C_index;
	CCString C_ccy;
	CCString C_cvname;
	double   C_date;
	long     C_seasonAdjId;
	long	 C_seasonAdjModeId;
};


/////////////////////////////////////////
/// Inflation curve from Summit create common function
/////////////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_GetInfZCFromSummit(
    LPXLOPER XL_index,
    LPXLOPER XL_ccy,
    LPXLOPER XL_cvname,
    LPXLOPER XL_date,
	LPXLOPER XL_seasonAdj,
	LPXLOPER XL_seasonAdjMode,
	bool PersistentInXL )
{
	ADD_LOG("Local_GetInfZCFromSummit");
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

		CCString C_index;
		CCString C_ccy;
		CCString C_cvname;
		double	 C_date;
		CCString C_seasonAdj;
		long seasonAdjId;
		CCString C_seasonAdjMode;
		long seasonAdjModeId;

		XL_readStrCell(  XL_index,			C_index,			" ARM_ERR: index: string expected",			C_result);
		XL_readStrCell(  XL_ccy,			C_ccy,				" ARM_ERR: currency: string expected",		C_result);
		XL_readStrCell(  XL_cvname,			C_cvname,			" ARM_ERR: cvname: string expected",		C_result);
		XL_readNumCell(	 XL_date,			C_date,				" ARM_ERR: as of date: double expected",	C_result);
		XL_readStrCellWD(XL_seasonAdj,		C_seasonAdj,"N",	" ARM_ERR: seasonAdj: string expected",		C_result);
		XL_readStrCellWD(XL_seasonAdjMode, C_seasonAdjMode, "PLUS", " ARM_ERR: seasonality correction mode : string expected", C_result);

		if((seasonAdjId = ARM_ConvYesOrNo (C_seasonAdj, C_result)) == ARM_DEFAULT_ERR)
		{
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}

		if((seasonAdjModeId = ARM_ConvSeasonalityMode (C_seasonAdjMode, C_result)) == ARM_DEFAULT_ERR)
		{
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}

		/// use the concept of Functor to transfer the knowledge of
		/// a function with a context
		getInfZcFromSummitFunc ourFunc( C_index, C_ccy, C_cvname, C_date, seasonAdjId, seasonAdjModeId);
		
		/// call the general function
		fillXL_Result( LOCAL_INFCURV_CLASS, ourFunc, C_result, XL_result, PersistentInXL );
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_GetInfZCFromSummit" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}



//////////////////////////////////////////////////////
/// Addin to create an inflation curve from Summit
//////////////////////////////////////////////////////

__declspec(dllexport) LPXLOPER WINAPI Local_GetInfZcFromSummit_Create(
    LPXLOPER XL_index,
    LPXLOPER XL_ccy,
    LPXLOPER XL_cvname,
    LPXLOPER XL_date,
	LPXLOPER XL_seasonAdj,
	LPXLOPER XL_seasonAdjMode)
{
	ADD_LOG("Local_GetInfZcFromSummit_Create");
	bool PersistentInXL = true;
	return Local_GetInfZCFromSummit(
		XL_index,
		XL_ccy,
		XL_cvname,
		XL_date,
		XL_seasonAdj,
		XL_seasonAdjMode,
		PersistentInXL);
}



//////////////////////////////////////////////////////
/// Addin to create an inflation curve from Summit
/// Version for VBA
//////////////////////////////////////////////////////

__declspec(dllexport) LPXLOPER WINAPI Local_PXL_GetInfZcFromSummit_Create(
    LPXLOPER XL_index,
    LPXLOPER XL_ccy,
    LPXLOPER XL_cvname,
    LPXLOPER XL_date,
	LPXLOPER XL_seasonAdj,
	LPXLOPER XL_seasonAdjMode)
{
	ADD_LOG("Local_PXL_GetInfZcFromSummit_Create");
	bool PersistentInXL = false;
	return Local_GetInfZCFromSummit(
		XL_index,
		XL_ccy,
		XL_cvname,
		XL_date,
		XL_seasonAdj,
		XL_seasonAdjMode,
		PersistentInXL);
}




class LivretACurveFctor : public ARMResultLong2LongFunc
{
public:
	LivretACurveFctor(	
		double asOfDateDble ,
		const CCString& infCurvId,
		const CCString& euribCurvId,
		long flagRounding,
		const CCString& infresetManagerId,
		const CCString& fixingLivretAId,
		const CCString& fixingEuribId,
		long monthForAugustId,
		long monthForFebruaryId)
	:
		C_asOfDateDble( asOfDateDble ),
		C_infCurvId( infCurvId ),
		C_euribCurvId( euribCurvId ),
		C_flagRounding( flagRounding ),
		C_infresetManagerId( infresetManagerId ),
		C_fixingLivretAId( fixingLivretAId ),
		C_fixingEuribId( fixingEuribId ),
		C_monthForAugustId (monthForAugustId),
		C_monthForFebruaryId (monthForFebruaryId)
	{};
	
	long operator()( ARM_result& result, long objId ) 
	{
		return ARMLOCAL_LIVRETACURVE (
			C_asOfDateDble,
			LocalGetNumObjectId (C_infCurvId),
			LocalGetNumObjectId (C_euribCurvId),
			C_flagRounding,
			LocalGetNumObjectId(C_infresetManagerId),
			LocalGetNumObjectId(C_fixingLivretAId),
			LocalGetNumObjectId(C_fixingEuribId), 
			C_monthForAugustId,
			C_monthForFebruaryId,
			result,
			objId);
	}

private:
	double C_asOfDateDble;
	CCString C_infCurvId;
	CCString C_euribCurvId;
	long C_flagRounding;
	CCString C_infresetManagerId;
	CCString C_fixingLivretAId;
	CCString C_fixingEuribId;
	long C_monthForAugustId;
	long C_monthForFebruaryId;

};



/////// addins for livretA
__declspec(dllexport) LPXLOPER WINAPI Local_ARM_LIVRETACURVE_COMMON(
	LPXLOPER XL_asOfDate,
	LPXLOPER XL_infCurvId,
	LPXLOPER XL_euribCurvId,
	LPXLOPER XL_flagRounding,
	LPXLOPER XL_infresetManagerId,
	LPXLOPER XL_fixingLivretAId,
	LPXLOPER XL_fixingEuriborId, 
	LPXLOPER XL_monthForAugust, 
	LPXLOPER XL_monthForFebruary, 
	bool PersistentInXL )
{
	ADD_LOG("Local_ARM_LIVRETACURVE_COMMON");
	static XLOPER XL_result;
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
		ARM_NOCALCIFWIZ();

		// C variable
		double C_flagRounding_default = 0.0;

		// error
		static int error;
		static char* reason = "";
		
		double C_asOfDate;
		XL_readNumCell(XL_asOfDate,C_asOfDate," ARM_ERR: as of date: date expected",C_result);

		CCString C_euribCurvId;
		XL_readStrCell(XL_euribCurvId,C_euribCurvId," ARM_ERR: forecast euribor curve id: object expected",C_result);

		CCString C_infCurveId;
		XL_readStrCell(XL_infCurvId,C_infCurveId," ARM_ERR: forecast inflation curve id: object expected",C_result);

		double C_flagRounding;
		XL_readNumCellWD(XL_flagRounding,C_flagRounding,C_flagRounding_default," ARM_ERR: flagRounding : 0 or 1 expected",C_result);

		CCString C_infresetManagerId;
		XL_readStrCellWD(XL_infresetManagerId,C_infresetManagerId,""," ARM_ERR: inflation resetManager id: object expected",C_result);

		CCString C_fixingLivretAId;
		XL_readStrCellWD(XL_fixingLivretAId,C_fixingLivretAId,""," ARM_ERR: fixingLivretA id: object expected",C_result);

		CCString C_fixingEuriborId;
		XL_readStrCellWD(XL_fixingEuriborId, C_fixingEuriborId,"", " ARM_ERR: fixingEuribor id: object expected",C_result);

		CCString C_monthForAugust;
		long monthForAugustId;
		XL_readStrCellWD(XL_monthForAugust, C_monthForAugust,"DEFAULT", " ARM_ERR: month for august: string expected",C_result);
		monthForAugustId = ARM_ConvMonth(C_monthForAugust);

		CCString C_monthForFebruary;
		long monthForFebruaryId;
		XL_readStrCellWD(XL_monthForFebruary, C_monthForFebruary,"DEFAULT", " ARM_ERR: month for february: string expected",C_result);
		monthForFebruaryId = ARM_ConvMonth(C_monthForFebruary);

		/// use the concept of Functor to transfer the knowledge of
		/// a function with a context
		LivretACurveFctor ourFunc( 
			C_asOfDate,
			C_infCurveId,
			C_euribCurvId,
			C_flagRounding,
			C_infresetManagerId,
			C_fixingLivretAId,			
			C_fixingEuriborId,
			monthForAugustId,
			monthForFebruaryId);
		
		/// call the general function
		fillXL_Result( LOCAL_LIVRET_A_CURVE_CLASS, ourFunc, C_result, XL_result, PersistentInXL );
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_ARM_LIVRETACURVE" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}



__declspec(dllexport) LPXLOPER WINAPI Local_ARM_LIVRETACURVE(
	LPXLOPER XL_asOfDate,
	LPXLOPER XL_infCurvId,
	LPXLOPER XL_euribCurvId,
	LPXLOPER XL_flagRounding,
	LPXLOPER XL_infresetManagerId,
	LPXLOPER XL_fixingLivretAId,
	LPXLOPER XL_fixingEuriborId,
	LPXLOPER XL_monthForAugust, 
	LPXLOPER XL_monthForFebruary) 
{
	ADD_LOG("Local_ARM_LIVRETACURVE");
	bool PersistentInXL = true;
	return Local_ARM_LIVRETACURVE_COMMON(
		XL_asOfDate,
		XL_infCurvId,
		XL_euribCurvId,
		XL_flagRounding,
		XL_infresetManagerId,
		XL_fixingLivretAId,
		XL_fixingEuriborId,
		XL_monthForAugust,
		XL_monthForFebruary,
		PersistentInXL);
}


__declspec(dllexport) LPXLOPER WINAPI Local_PXL_ARM_LIVRETACURVE(
	LPXLOPER XL_asOfDate,
	LPXLOPER XL_infCurvId,
	LPXLOPER XL_euribCurvId,
	LPXLOPER XL_flagRounding,
	LPXLOPER XL_infresetManagerId,
	LPXLOPER XL_fixingLivretAId,
	LPXLOPER XL_fixingEuriborId,
	LPXLOPER XL_monthForAugust,
	LPXLOPER XL_monthForFebruary) 
{
	ADD_LOG("Local_PXL_ARM_LIVRETACURVE");
	bool PersistentInXL = false;
	return Local_ARM_LIVRETACURVE_COMMON(
		XL_asOfDate,
		XL_infCurvId,
		XL_euribCurvId,
		XL_flagRounding,
		XL_infresetManagerId,
		XL_fixingLivretAId,
		XL_fixingEuriborId,
		XL_monthForAugust,
		XL_monthForFebruary,
		PersistentInXL);
}



__declspec(dllexport) LPXLOPER WINAPI Local_ARM_LivreACurveGetRateDate(LPXLOPER XL_livretACurvId,
																	   LPXLOPER XL_dateIn)
{
	ADD_LOG("Local_ARM_LivreACurveGetRateDate");
	static XLOPER XL_result;
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
		ARM_NOCALCIFWIZ();

		// C variable
		CCString C_livretACurvId;
		double C_dateIn;
		
		// error
		static int error;
		static char* reason = "";
		
		XL_readNumCell(XL_dateIn,C_dateIn," ARM_ERR: as of date: date expected",C_result);
		XL_readStrCell(XL_livretACurvId,C_livretACurvId," ARM_ERR: forecast livret A curve id: object expected",C_result);
		long retCode;

		retCode = ARMLOCAL_LIVRETACURVEGETRATEDATE (LocalGetNumObjectId(C_livretACurvId), C_dateIn, C_result);

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
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_ARM_LivreACurveGetRateDate" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}



////////////////////////////////////////////////
/// very rapid definition of ourlonglongFunctor
////////////////////////////////////////////////
class InfSwoVolCurvFctor : public ARMResultLong2LongFunc
{
public:
	InfSwoVolCurvFctor(
		const ARM_Date& asOfDate,
		long InfIRModelId,
		const VECTOR<double>& tenors,
		const VECTOR<double>& expiries,
		long computationMethod)
	:
		C_asOfDate(asOfDate),
		C_InfIRModelId(InfIRModelId ),
		C_tenors(tenors),
		C_expiries(expiries),
		C_ComputationMethod(computationMethod)
	{};
	
	long operator()( ARM_result& result, long objId ) 
	{
		return ARMLOCAL_InfSwoVolCurveFromModel_Create(
			C_asOfDate,
			C_InfIRModelId,
			C_tenors,
			C_expiries,
			C_ComputationMethod,
			result,
			objId );			
	}
private:
	ARM_Date C_asOfDate;
	long C_InfIRModelId;
	VECTOR<double> C_tenors;
	VECTOR<double> C_expiries;
	long C_ComputationMethod;
};



/////////////////////////////////////////////////////////////
/// central function that does the creation of the XL function
/////////////////////////////////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_InfSwoVolCurveFromModel_Common(
	LPXLOPER XL_asOfDate,
	LPXLOPER XL_InfIRModelId,
	LPXLOPER XL_tenors,
	LPXLOPER XL_expiries,
	LPXLOPER XL_ComputationMethod,
	bool PersistentInXL )
{
	ADD_LOG("Local_InfSwoVolCurveFromModel_Common");
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

		double C_asOfDble;
		XL_readNumCell(	 XL_asOfDate,	C_asOfDble,		" ARM_ERR: asOfDate: date expected", C_result);
		ARM_Date C_asOfDate = ConvertToARMDATE(C_asOfDble);

		long C_InfIRModelId;
		XL_GETOBJID( XL_InfIRModelId, C_InfIRModelId,	" ARM_ERR: Inflation model id: object expected", C_result);

		VECTOR<double> C_tenors;
		VECTOR<double> C_tenorsDefault;
		XL_readNumVectorWD(XL_tenors,C_tenors,C_tenorsDefault," ARM_ERR: tenors: array of numeric expected", C_result);

		VECTOR<double> C_expiries;
		VECTOR<double> C_expiriesDefault;
		XL_readNumVectorWD(XL_expiries,C_expiries,C_expiriesDefault," ARM_ERR: tenors: array of numeric expected",C_result);

		CCString C_ComputationMethodStr;
		CCString defaultComputationMethodStr = "Std";
		XL_readStrCellWD(XL_ComputationMethod,C_ComputationMethodStr,defaultComputationMethodStr," ARM_ERR: computation method: String expected",C_result);

		long C_ComputationMethod;
		if( (C_ComputationMethod = ARM_ConvGPInfSwoptComputationMethod( C_ComputationMethodStr, C_result)) == ARM_DEFAULT_ERR )
		{
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}


		/// use the concept of Functor to transfer the knowledge of
		/// a function with a context
		InfSwoVolCurvFctor ourFunc(
			C_asOfDate,
			C_InfIRModelId,
			C_tenors,
			C_expiries,
			C_ComputationMethod);

		/// call the general function
		fillXL_Result( LOCAL_VOL_CURVE_LIN_CLASS, ourFunc, C_result, XL_result, PersistentInXL );
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_InfSwoVolCurveFromModel_Common" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}



//////////////////////////////////////////////////////
/// Addin to create an inflation volatility swaption from a model
//////////////////////////////////////////////////////

__declspec(dllexport) LPXLOPER WINAPI Local_InfSwoVolCurveFromModel_Create(
	LPXLOPER XL_asOfDate,
	LPXLOPER XL_InfIRModelId,
	LPXLOPER XL_tenors,
	LPXLOPER XL_expiries,
	LPXLOPER XL_ComputationMethod )
{
	ADD_LOG("Local_InfSwoVolCurveFromModel_Create");
	bool PersistentInXL = true;
	return 	Local_InfSwoVolCurveFromModel_Common(
		XL_asOfDate,
		XL_InfIRModelId,
		XL_tenors,
		XL_expiries,
		XL_ComputationMethod,
		PersistentInXL );
}


//////////////////////////////////////////////////////
/// Addin to create an inflation volatility swaption from a model
/// Version for VBA
//////////////////////////////////////////////////////

__declspec(dllexport) LPXLOPER WINAPI Local_PXL_InfSwoVolCurveFromModel_Create(
	LPXLOPER XL_asOfDate,
	LPXLOPER XL_InfIRModelId,
	LPXLOPER XL_tenors,
	LPXLOPER XL_expiries,
	LPXLOPER XL_ComputationMethod )
{
	ADD_LOG("Local_PXL_InfSwoVolCurveFromModel_Create");
	bool PersistentInXL = false;
	return 	Local_InfSwoVolCurveFromModel_Common(
		XL_asOfDate,
		XL_InfIRModelId,
		XL_tenors,
		XL_expiries,
		XL_ComputationMethod,
		PersistentInXL );
}

//YK
////////////////////////////////////////////////
/// very rapid definition of ourlonglongFunctor
////////////////////////////////////////////////
class InfSwoVolCubeFctor : public ARMResultLong2LongFunc
{
public:
	InfSwoVolCubeFctor(
		const ARM_Date& asOfDate,
		long InfIRModelId,
		const VECTOR<double>& tenors,
		const VECTOR<double>& expiries,
		const VECTOR<double>& smiledTenors,
		const VECTOR<double>& strikes,
		long computationMethod)
	:
		C_asOfDate(asOfDate),
		C_InfIRModelId(InfIRModelId ),
		C_tenors(tenors),
		C_expiries(expiries),
		C_smiledTenors(smiledTenors),
		C_strikes(strikes),
		C_ComputationMethod(computationMethod)
	{};
	
	long operator()( ARM_result& result, long objId ) 
	{
		return ARMLOCAL_InfSwoVolCubeFromModel_Create(
			C_asOfDate,
			C_InfIRModelId,
			C_tenors,
			C_expiries,
			C_smiledTenors,
			C_strikes,
			C_ComputationMethod,
			result,
			objId );			
	}
private:
	ARM_Date C_asOfDate;
	long C_InfIRModelId;
	VECTOR<double> C_tenors;
	VECTOR<double> C_expiries;
	VECTOR<double> C_smiledTenors;
	VECTOR<double> C_strikes;
	long C_ComputationMethod;
};



/////////////////////////////////////////////////////////////
/// central function that does the creation of the XL function
/////////////////////////////////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_InfSwoVolCubeFromModel_Common(
	LPXLOPER XL_asOfDate,
	LPXLOPER XL_InfIRModelId,
	LPXLOPER XL_tenors,
	LPXLOPER XL_expiries,
	LPXLOPER XL_smiledTenors,
	LPXLOPER XL_strikes,
	LPXLOPER XL_ComputationMethod,
	bool PersistentInXL )
{
	ADD_LOG("Local_InfSwoVolCubeFromModel_Common");
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

		double C_asOfDble;
		XL_readNumCell(	 XL_asOfDate,	C_asOfDble,		" ARM_ERR: asOfDate: date expected", C_result);
		ARM_Date C_asOfDate = ConvertToARMDATE(C_asOfDble);

		long C_InfIRModelId;
		XL_GETOBJID( XL_InfIRModelId, C_InfIRModelId,	" ARM_ERR: Inflation model id: object expected", C_result);

		VECTOR<double> C_tenors;
		VECTOR<double> C_tenorsDefault;
		XL_readNumVectorWD(XL_tenors,C_tenors,C_tenorsDefault," ARM_ERR: tenors: array of numeric expected", C_result);

		VECTOR<double> C_expiries;
		VECTOR<double> C_expiriesDefault;
		XL_readNumVectorWD(XL_expiries,C_expiries,C_expiriesDefault," ARM_ERR: expiries: array of numeric expected",C_result);

		
		VECTOR<double> C_smiledTenors;
		VECTOR<double> C_smiledTenorsDefault;
		XL_readNumVectorWD(XL_smiledTenors,C_smiledTenors,C_smiledTenorsDefault," ARM_ERR: expiries: array of numeric expected",C_result);

		VECTOR<double> C_strikes;
		VECTOR<double> C_strikesDefault;
		XL_readNumVectorWD(XL_strikes,C_strikes,C_strikesDefault," ARM_ERR: strikes: array of numeric expected",C_result);

		CCString C_ComputationMethodStr;
		CCString defaultComputationMethodStr = "Std";
		XL_readStrCellWD(XL_ComputationMethod,C_ComputationMethodStr,defaultComputationMethodStr," ARM_ERR: computation method: String expected",C_result);

		long C_ComputationMethod;
		if( (C_ComputationMethod = ARM_ConvGPInfSwoptComputationMethod( C_ComputationMethodStr, C_result)) == ARM_DEFAULT_ERR )
		{
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}


		/// use the concept of Functor to transfer the knowledge of
		/// a function with a context
		InfSwoVolCubeFctor ourFunc(
			C_asOfDate,
			C_InfIRModelId,
			C_tenors,
			C_expiries,
			C_smiledTenors,
			C_strikes,
			C_ComputationMethod);

		/// call the general function
		fillXL_Result( LOCAL_VOL_CURVE_LIN_CLASS, ourFunc, C_result, XL_result, PersistentInXL );
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_InfSwoVolCurveFromModel_Common" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}



//////////////////////////////////////////////////////
/// Addin to create an inflation volatility swaption from a model
//////////////////////////////////////////////////////

__declspec(dllexport) LPXLOPER WINAPI Local_InfSwoVolCubeFromModel_Create(
	LPXLOPER XL_asOfDate,
	LPXLOPER XL_InfIRModelId,
	LPXLOPER XL_tenors,
	LPXLOPER XL_expiries,
	LPXLOPER XL_smiledTenors,
	LPXLOPER XL_strikes,
	LPXLOPER XL_ComputationMethod )
{
	ADD_LOG("Local_InfSwoVolCubeFromModel_Create");
	bool PersistentInXL = true;
	return 	Local_InfSwoVolCubeFromModel_Common(
		XL_asOfDate,
		XL_InfIRModelId,
		XL_tenors,
		XL_expiries,
		XL_smiledTenors,
		XL_strikes,
		XL_ComputationMethod,
		PersistentInXL );
}


//////////////////////////////////////////////////////
/// Addin to create an inflation volatility swaption from a model
/// Version for VBA
//////////////////////////////////////////////////////

__declspec(dllexport) LPXLOPER WINAPI Local_PXL_InfSwoVolCubeFromModel_Create(
	LPXLOPER XL_asOfDate,
	LPXLOPER XL_InfIRModelId,
	LPXLOPER XL_tenors,
	LPXLOPER XL_expiries,
	LPXLOPER XL_smiledTenors,
	LPXLOPER XL_strikes,
	LPXLOPER XL_ComputationMethod )
{
	ADD_LOG("Local_PXL_InfSwoVolCubeFromModel_Create");
	bool PersistentInXL = false;
	return 	Local_InfSwoVolCubeFromModel_Common(
		XL_asOfDate,
		XL_InfIRModelId,
		XL_tenors,
		XL_expiries,
		XL_smiledTenors,
		XL_strikes,		
		XL_ComputationMethod,
		PersistentInXL );
}

//YK


////////////////////////////////////////////////
/// very rapid definition of ourlonglongFunctor
////////////////////////////////////////////////
class InfOATSwoVolCurvFctor : public ARMResultLong2LongFunc
{
public:
	InfOATSwoVolCurvFctor(
		const	ARM_Date& asOfDate,
		long	InfIRModelId,
		const	VECTOR<double>& tenors,
		const	VECTOR<double>& expiries,
		double	coupon,
		long	choice,
		long	ComputationMethod)
	:
		C_asOfDate(asOfDate),
		C_InfIRModelId(InfIRModelId ),
		C_tenors(tenors),
		C_expiries(expiries),
		C_coupon(coupon),
		C_choice(choice),
		C_ComputationMethod(ComputationMethod)
		{};
	
	long operator()( ARM_result& result, long objId )
	{
		return ARMLOCAL_InfOATSwoVolCurveFromModel_Create(
			C_asOfDate,
			C_InfIRModelId,
			C_tenors,
			C_expiries,
			C_coupon,
			C_choice,
			C_ComputationMethod,
			result,
			objId );			
	}
private:
	ARM_Date		C_asOfDate;
	long			C_InfIRModelId;
	VECTOR<double>	C_tenors;
	VECTOR<double>	C_expiries;
	const double	C_coupon;
	long			C_choice;
	long			C_ComputationMethod;
};



/////////////////////////////////////////////////////////////
/// central function that does the creation of the XL function
/////////////////////////////////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_InfOATSwoVolCurveFromModel_Common(
	LPXLOPER XL_asOfDate,
	LPXLOPER XL_InfIRModelId,
	LPXLOPER XL_tenors,
	LPXLOPER XL_expiries,
	LPXLOPER XL_ComputationMethod,
	LPXLOPER XL_coupon,
	LPXLOPER XL_choice,
	bool PersistentInXL )
{
	ADD_LOG("Local_InfOATSwoVolCurveFromModel_Common");
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

		double C_asOfDble;
		XL_readNumCell(	 XL_asOfDate,	C_asOfDble,		" ARM_ERR: asOfDate: date expected", C_result);
		ARM_Date C_asOfDate = ConvertToARMDATE(C_asOfDble);

		long C_InfIRModelId;
		XL_GETOBJID( XL_InfIRModelId, C_InfIRModelId,	" ARM_ERR: Inflation model id: object expected", C_result);

		VECTOR<double> C_tenors;
		VECTOR<double> C_tenorsDefault;
		XL_readNumVectorWD(XL_tenors, C_tenors, C_tenorsDefault, " ARM_ERR: tenors: array of numeric expected", C_result);

		VECTOR<double> C_expiries;
		VECTOR<double> C_expiriesDefault;
		XL_readNumVectorWD(XL_expiries, C_expiries, C_expiriesDefault, " ARM_ERR: tenors: array of numeric expected", C_result);

		CCString C_ComputationMethodStr;
		CCString defaultComputationMethodStr = "Std";
		XL_readStrCellWD(XL_ComputationMethod, C_ComputationMethodStr, defaultComputationMethodStr, " ARM_ERR: computation method: String expected", C_result);

		double C_coupon;
		double C_couponDefalut = 1.;
		//CCString C_choiceStr;
		//CCString C_choiceDefaultStr = "Yes";
		//XL_readStrCellWD(XL_choice, C_choiceDefaultStr, C_choiceDefaultStr, " ARM_ERR: String Yes or No expected", C_result);
		
		long C_ComputationMethod;
		if( (C_ComputationMethod = ARM_ConvGPInfSwoptComputationMethod( C_ComputationMethodStr, C_result)) == ARM_DEFAULT_ERR )
		{
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}

		double C_choice;
		double C_choiceDefalut = 0;
		XL_readNumCellWD( XL_choice, C_choice, C_choiceDefalut, " ARM_ERR: coupon: double expected", C_result);
		
		if ( C_choice == 0 )
		{
			XL_readNumCellWD( XL_coupon, C_coupon, C_couponDefalut, " ARM_ERR: coupon: double expected", C_result);
		}
		else
		{
			C_choice = 1;
			XL_readNumCell( XL_coupon, C_coupon, " ARM_ERR: coupon: double expected", C_result);
			
		}

				
		/// use the concept of Functor to transfer the knowledge of
		/// a function with a context
		InfOATSwoVolCurvFctor ourFunc(
			C_asOfDate,
			C_InfIRModelId,
			C_tenors,
			C_expiries,
			C_coupon,
			C_choice,
			C_ComputationMethod);

		

		/// call the general function
		fillXL_Result( LOCAL_VOL_CURVE_LIN_CLASS, ourFunc, C_result, XL_result, PersistentInXL );
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_InfSwoVolCurveFromModel_Common" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}

//////////////////////////////////////////////////////
/// Addin to create an inflation volatility swaption from a model
//////////////////////////////////////////////////////

__declspec(dllexport) LPXLOPER WINAPI Local_InfOATSwoVolCurveFromModel_Create(
									  									  
	LPXLOPER XL_asOfDate,
	LPXLOPER XL_InfIRModelId,
	LPXLOPER XL_tenors,
	LPXLOPER XL_expiries,
	LPXLOPER XL_coupon,
	LPXLOPER XL_choice,
	LPXLOPER XL_ComputationMethod)
{
	ADD_LOG("Local_InfOATSwoVolCurveFromModel_Create");
	bool PersistentInXL = true;
	return 	Local_InfOATSwoVolCurveFromModel_Common(
		XL_asOfDate,
		XL_InfIRModelId,
		XL_tenors,
		XL_expiries,
		XL_coupon,
		XL_choice,
		XL_ComputationMethod,

		PersistentInXL );
}





////////////////////////////////////////////////
/// functor for InfYCMod
////////////////////////////////////////////////
class InfMultiBSModFunc : public ARMResultLong2LongFunc
{
public:
	InfMultiBSModFunc( const VECTOR<long>& infMultiBSModIdVec )
	:	C_infMultiBSModIdVec( infMultiBSModIdVec)
	{};

	long operator()( ARM_result& result, long objId ) 
	{
		return ARMLOCAL_infMultiBSMod_Create( C_infMultiBSModIdVec, result, objId );
	}
private:
	VECTOR<long> C_infMultiBSModIdVec;
};


////////////////////////////////////////////////
/// general function to branch
////////////////////////////////////////////////

LPXLOPER Local_InfMultiBSMod_Common(
	LPXLOPER XL_infMultiBSModIds,
	bool PersistentInXL )
{
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

		VECTOR<CCString>	C_infMultiBSModStrVec;
		VECTOR<long> C_infMultiBSModIds;
		
		XL_readStrVector(XL_infMultiBSModIds,C_infMultiBSModStrVec," ARM_ERR: inflation bs model: array of string expected",XL_TYPE_STRING,C_result);
		for (int i=0; i<C_infMultiBSModStrVec.size();i++)
			C_infMultiBSModIds.push_back(LocalGetNumObjectId(C_infMultiBSModStrVec[i]));
		InfMultiBSModFunc ourFunc( C_infMultiBSModIds );

		/// call the general function
		fillXL_Result( LOCAL_INFMULTIBSMODEL_CLASS, ourFunc, C_result, XL_result, PersistentInXL );
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_InfBSMod_Common" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


////////////////////////////////////////////////
/// Version for XL exportation of the InfYCMod
////////////////////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_InfMultiBSMOD(
	LPXLOPER XL_infMultiBSModIdVec )
{
	ADD_LOG("Local_InfMultiBSMOD");
	bool PersistentInXL = true;
	return Local_InfMultiBSMod_Common(
		XL_infMultiBSModIdVec,
		PersistentInXL );
}



////////////////////////////////////////////////
/// Version for XL exportation 
////////////////////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_PXL_InfMultiBSMOD(
	LPXLOPER XL_infMultiBSModIdVec )
{
	ADD_LOG("Local_PXL_InfMultiBSMOD");
	bool PersistentInXL = false;
	return Local_InfMultiBSMod_Common(
		XL_infMultiBSModIdVec,
		PersistentInXL );
}


///----------------------------------------------
///----------------------------------------------
///             set reset manager
/// Inputs :
///     resetManager id
///     inf curve id
///----------------------------------------------
///----------------------------------------------


/////////////////////////////////////////////////////////////
/// central function that does the creation of the XL function
/////////////////////////////////////////////////////////////
LPXLOPER Local_InfCurv_SetResetManager_Common(
	LPXLOPER XL_infCurvId,
	LPXLOPER XL_ResetManagerId,
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

		long C_infCurvId;
		XL_GETOBJID( XL_infCurvId, C_infCurvId, " ARM_ERR: inf curve: Obj expected", C_result );
		
		long C_ResetManagerId;
		XL_GETOBJID( XL_ResetManagerId, C_ResetManagerId, " ARM_ERR: Reset Manager: Obj expected", C_result );

		/// use the concept of Functor to transfer the knowledge of
		/// a function with a context
		exportFunc2Args< long, long >  ourFunc( C_infCurvId, C_ResetManagerId, ARMLOCAL_InfCurv_SetResetManager );

		/// call the general function
		fillXL_Result( LOCAL_INFCURV_CLASS, ourFunc, C_result, XL_result, PersistentInXL );
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_InfCurv_SetResetManager_Common" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}

							 


///////////////////////////////////
/// version that takes into account 
/// previous creation of object
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_InfCurv_SetResetManager(
	LPXLOPER XL_infCurvId,
	LPXLOPER XL_ResetManagerId )
{
	ADD_LOG("Local_InfCurv_SetResetManager");
	bool PersistentInXL = true;
	return Local_InfCurv_SetResetManager_Common(
		XL_infCurvId,
		XL_ResetManagerId,
		PersistentInXL );
}


///////////////////////////////////
/// version for non persistentInXl
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_PXL_InfCurv_SetResetManager(
	LPXLOPER XL_infCurvId,
	LPXLOPER XL_ResetManagerId )
{
	ADD_LOG("Local_PXL_InfCurv_SetResetManager");
	bool PersistentInXL = false;
	return Local_InfCurv_SetResetManager_Common(
		XL_infCurvId,
		XL_ResetManagerId,
		PersistentInXL );
}




///----------------------------------------------
///----------------------------------------------
///             set season manager
/// Inputs :
///     seasonality manager id
///     inf curve id
///----------------------------------------------
///----------------------------------------------


/////////////////////////////////////////////////////////////
/// central function that does the creation of the XL function
/////////////////////////////////////////////////////////////
LPXLOPER Local_InfCurv_SetSeasonalitytManager_Common(
	LPXLOPER XL_infCurvId,
	LPXLOPER XL_SeasonalitytManagerId,
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

		long C_infCurvId;
		XL_GETOBJID( XL_infCurvId, C_infCurvId, " ARM_ERR: inf curve: Obj expected", C_result );
		
		long C_SeasonalitytManagerId;
		XL_GETOBJID( XL_SeasonalitytManagerId, C_SeasonalitytManagerId, " ARM_ERR: SeasonalitytManager: Obj expected", C_result );

		/// use the concept of Functor to transfer the knowledge of
		/// a function with a context
		exportFunc2Args< long, long >  ourFunc(C_infCurvId, C_SeasonalitytManagerId, ARMLOCAL_InfCurv_SetSeasonalityManager );

		/// call the general function
		fillXL_Result( LOCAL_INFCURV_CLASS, ourFunc, C_result, XL_result, PersistentInXL );
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_InfCurv_SetSeasonalitytManager_Common" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}

							 


///////////////////////////////////
/// version that takes into account 
/// previous creation of object
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_InfCurv_SetSeasonalityManager(
	LPXLOPER XL_infCurvId,
	LPXLOPER XL_SeasonalitytManagerId )
{
	ADD_LOG("Local_InfCurv_SetSeasonalityManager");
	bool PersistentInXL = true;
	return Local_InfCurv_SetSeasonalitytManager_Common(
		XL_infCurvId,
		XL_SeasonalitytManagerId,
		PersistentInXL );
}


///////////////////////////////////
/// version for non persistentInXl
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_PXL_InfCurv_SetSeasonalityManager(
	LPXLOPER XL_infCurvId,
	LPXLOPER XL_SeasonalitytManagerId )
{
	ADD_LOG("Local_PXL_InfCurv_SetSeasonalityManager");
	bool PersistentInXL = false;
	return Local_InfCurv_SetSeasonalitytManager_Common(
		XL_infCurvId,
		XL_SeasonalitytManagerId,
		PersistentInXL );
}




__declspec(dllexport) LPXLOPER WINAPI Local_ARM_GetSeasonMgrFromSummit(LPXLOPER XL_index,
																	   LPXLOPER XL_ccy,
																	   LPXLOPER XL_cvname,
																	   LPXLOPER XL_asof,
																	   LPXLOPER XL_mode)
{
	ADD_LOG("Local_ARM_GetSeasonMgrFromSummit");
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
	CCString C_index;
	CCString C_ccy;
	CCString C_cvname;

	double C_asOf;

	CCString C_mode;
	long modeId;
	
	// error
	static int error;
	static char* reason = "";

	XL_readStrCell(XL_index,C_index," ARM_ERR: index: string expected",C_result);
	XL_readStrCell(XL_ccy,C_ccy," ARM_ERR: ccy: string expected",C_result);
	XL_readStrCell(XL_cvname,C_cvname," ARM_ERR: cvname: string expected",C_result);
	XL_readNumCell(XL_asof,C_asOf," ARM_ERR: asofdate: numeric expected",C_result);
	XL_readStrCellWD(XL_mode, C_mode, "PLUS", " ARM_ERR: seasonality correction mode : string expected", C_result);

	if((modeId = ARM_ConvSeasonalityMode (C_mode, C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	CCString prevClass;
	
	CCString curClass = LOCAL_SEASONMANAGER_CLASS;
	CCString stringId = GetLastCurCellEnvValue ();
	
	long objId;

	if(!stringId)
	{
		retCode = ARMLOCAL_GetSeasonMgrFromSummit (C_index,
												   C_ccy,
												   C_cvname,
												   C_asOf,
												   modeId,
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
			retCode = ARMLOCAL_GetSeasonMgrFromSummit (C_index,
													   C_ccy,
													   C_cvname,
													   C_asOf,
													   modeId,
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
			retCode = ARMLOCAL_GetSeasonMgrFromSummit (C_index,
													   C_ccy,
													   C_cvname,
													   C_asOf,
													   modeId,
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_ARM_GetSeasonMgrFromSummit" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}



__declspec(dllexport) LPXLOPER WINAPI Local_PXL_ARM_GetSeasonMgrFromSummit(LPXLOPER XL_index,
																		   LPXLOPER XL_ccy,
																		   LPXLOPER XL_cvname,
																		   LPXLOPER XL_asof,
																		   LPXLOPER XL_mode)
{
	ADD_LOG("Local_PXL_ARM_GetSeasonMgrFromSummit");
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
	CCString C_index;
	CCString C_ccy;
	CCString C_cvname;

	double C_asOf;

	CCString C_mode;
	long modeId;
	
	// error
	static int error;
	static char* reason = "";

	XL_readStrCell(XL_index,C_index," ARM_ERR: index: string expected",C_result);
	XL_readStrCell(XL_ccy,C_ccy," ARM_ERR: ccy: string expected",C_result);
	XL_readStrCell(XL_cvname,C_cvname," ARM_ERR: cvname: string expected",C_result);
	XL_readNumCell(XL_asof,C_asOf," ARM_ERR: asofdate: numeric expected",C_result);
	XL_readStrCellWD(XL_mode, C_mode, "PLUS", " ARM_ERR: seasonality correction mode : string expected", C_result);

	if((modeId = ARM_ConvSeasonalityMode (C_mode, C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	CCString prevClass;
	
	CCString curClass = LOCAL_SEASONMANAGER_CLASS;
	CCString stringId;
	long objId;
	
	retCode = ARMLOCAL_GetSeasonMgrFromSummit (C_index,
											   C_ccy,
											   C_cvname,
											   C_asOf,
											   modeId,
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_PXL_ARM_GetSeasonMgrFromSummit" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}



__declspec(dllexport) LPXLOPER WINAPI Local_ARM_GetResetMgrFromSummit(LPXLOPER XL_asof,
																	  LPXLOPER XL_index,
																	  LPXLOPER XL_source,
																	  LPXLOPER XL_ccy,
																	  LPXLOPER XL_isInflatIndex,
																	  LPXLOPER XL_term)
{
	ADD_LOG("Local_ARM_GetResetMgrFromSummit");
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
	double C_asof;
	CCString C_index;
	CCString C_source;
	CCString C_ccy;
	CCString C_isInflatIndex;
	CCString C_term;
	
	// error
	static int error;
	static char* reason = "";

	XL_readNumCell(XL_asof,C_asof," ARM_ERR: asof: numeric expected",C_result);
	XL_readStrCell(XL_index,C_index," ARM_ERR: Index: string expected",C_result);
	XL_readStrCellWD(XL_source,C_source,"MO"," ARM_ERR: Source: string expected",C_result);
	XL_readStrCellWD(XL_ccy,C_ccy,"EUR"," ARM_ERR: Currency: string expected",C_result);
	XL_readStrCellWD(XL_isInflatIndex,C_isInflatIndex,"Y"," ARM_ERR: Currency: string expected",C_result);

	long isInflatIndexId;

	if((isInflatIndexId = ARM_ConvYesOrNo (C_isInflatIndex, C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	if (isInflatIndexId == K_YES)
	{
		XL_readStrCellWD(XL_term,C_term,"1D"," ARM_ERR: Term: string expected",C_result);
	}
	else
	{
		XL_readStrCell(XL_term,C_term," ARM_ERR: Term: string expected",C_result);
	}

	CCString prevClass;
	
	CCString curClass = LOCAL_RESETMANAGER_CLASS;
	CCString stringId = GetLastCurCellEnvValue ();	

	long objId;

	if(!stringId)
	{
		retCode = ARMLOCAL_GetResetMgrFromSummit (C_asof,
												  C_index,
												  C_source,
												  C_ccy,
												  isInflatIndexId,
												  C_term,
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
			retCode = ARMLOCAL_GetResetMgrFromSummit (C_asof,
													  C_index,
													  C_source,
													  C_ccy,
													  isInflatIndexId,
													  C_term,
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
			retCode = ARMLOCAL_GetResetMgrFromSummit (C_asof,
													  C_index,
													  C_source,
													  C_ccy,
													  isInflatIndexId,
													  C_term,
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_ARM_GetResetMgrFromSummit" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


__declspec(dllexport) LPXLOPER WINAPI Local_PXL_ARM_GetResetMgrFromSummit(LPXLOPER XL_asof,
																		  LPXLOPER XL_index,
																		  LPXLOPER XL_source,
																		  LPXLOPER XL_ccy,
																		  LPXLOPER XL_isInflatIndex,
																		  LPXLOPER XL_term)
{
	ADD_LOG("Local_PXL_ARM_GetResetMgrFromSummit");
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
	double C_asof;
	CCString C_index;
	CCString C_source;
	CCString C_ccy;
	CCString C_isInflatIndex;
	CCString C_term;
	
	// error
	static int error;
	static char* reason = "";

	XL_readNumCell(XL_asof,C_asof," ARM_ERR: asof: numeric expected",C_result);
	XL_readStrCell(XL_index,C_index," ARM_ERR: Index: string expected",C_result);
	XL_readStrCell(XL_source,C_source," ARM_ERR: Source: string expected",C_result);
	XL_readStrCellWD(XL_ccy,C_ccy,"EUR"," ARM_ERR: Currency: string expected",C_result);
	XL_readStrCellWD(XL_isInflatIndex,C_isInflatIndex,"Y"," ARM_ERR: Currency: string expected",C_result);

	long isInflatIndexId;

	if((isInflatIndexId = ARM_ConvYesOrNo (C_isInflatIndex, C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	if (isInflatIndexId == K_YES)
	{
		XL_readStrCellWD(XL_term,C_term,"1D"," ARM_ERR: Term: string expected",C_result);
	}
	else
	{
		XL_readStrCell(XL_term,C_term," ARM_ERR: Term: string expected",C_result);
	}

	CCString prevClass;
	
	CCString curClass = LOCAL_RESETMANAGER_CLASS;
	CCString stringId;
	long objId;
	
	retCode = ARMLOCAL_GetResetMgrFromSummit (C_asof,
											  C_index,
											  C_source,
											  C_ccy,
											  isInflatIndexId,
											  C_term,
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_PXL_ARM_GetResetMgrFromSummit" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}

////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////
/// very rapid definition of ourlonglongFunctor
///////////////////////////////////////////////
class GPInfCapFloorFctor : public ARMResultLong2LongFunc
{
public:
	GPInfCapFloorFctor(	
		long swapId,
		long CF,
		double strike,
        long strikeId)
	:
		C_swapId( swapId ),
		C_CF( CF ),
		C_strike( strike ),
 		C_strikeId( strikeId )
	{};
	
	long operator()( ARM_result& result, long objId ) 
	{
		return ARMLOCAL_GP_INFCAPFLOOR(
			C_swapId,
			C_CF,
			C_strike,
 			C_strikeId,
			result,								
			objId);
	}

private:
	long C_swapId;
	long C_CF;
	double C_strike;
    long C_strikeId;
};



/////////////////////////////////////////////////////////////
/// central function that does the creation of the XL function
/////////////////////////////////////////////////////////////
LPXLOPER Local_ARM_GP_INFCAPFLOORCommon(
	LPXLOPER XL_swapId,
	LPXLOPER XL_CAPOrFLOOR,
	LPXLOPER XL_strike,
	bool PersistentInXL )
{
	// this is defined first because it is used in XL macros
	static XLOPER XL_result;
	ARM_result C_result;
	// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
		/// to avoid computation if called by the wizard
		ARM_NOCALCIFWIZ();
		
		/// Get the variables from the XLOper variables
		CCString swapStr;
		long swapId;

		CCString CF;
		long CFId;

		//strike
		double strike;
		CCString strikeStr;
		long     strikeId;

		// error
		static int error;
		static char* reason = "";
		
		XL_readStrOrNumCell(XL_strike, strikeStr, strike, strikeId,
			   " ARM_ERR: strike: numeric or srike Id expected",C_result);

		XL_readStrCell(XL_swapId,swapStr," ARM_ERR: swap id: object expected",C_result);
		XL_readStrCell(XL_CAPOrFLOOR,CF," ARM_ERR: Cap or Floor: string expected",C_result);

		swapId = LocalGetNumObjectId (swapStr);
		if((CFId = ARM_ConvCapOrFloor (CF, C_result)) == ARM_DEFAULT_ERR)
		{
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}
		
		if(strikeId == XL_TYPE_STRING)
			strikeId = LocalGetNumObjectId(strikeStr);
		else
			strikeId = ARM_NULL_OBJECT;

		
		GPInfCapFloorFctor ourFunc(
				swapId,
				CFId,
				strike,
 				strikeId);

		/// call the general function
		fillXL_Result( LOCAL_GP_INFCAPFLOOR_CLASS, ourFunc, C_result, XL_result, PersistentInXL );
	
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_ARM_GP_INFCAPFLOORCommon" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


//////////////////////////////////
/// version that takes into account 
/// previous creation of object
//////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_ARM_GP_INFCAPFLOOR(
	 LPXLOPER XL_swapId,
	 LPXLOPER XL_CAPOrFLOOR,
	 LPXLOPER XL_strike)
{
	ADD_LOG("Local_ARM_GP_INFCAPFLOOR");
	 bool PersistentInXL = true;
	 return Local_ARM_GP_INFCAPFLOORCommon(
		XL_swapId,
		XL_CAPOrFLOOR,
		XL_strike,
 		PersistentInXL);
}

///////////////////////////////////
/// version for non persistentInXl
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_PXL_ARM_GP_INFCAPFLOOR(
	 LPXLOPER XL_swapId,
	 LPXLOPER XL_CAPOrFLOOR,
	 LPXLOPER XL_strike)
{
	ADD_LOG("Local_PXL_ARM_GP_INFCAPFLOOR");
	 bool PersistentInXL = false;
	 return Local_ARM_GP_INFCAPFLOORCommon(
		XL_swapId,
		XL_CAPOrFLOOR,
		XL_strike,
 		PersistentInXL);
}
			/************************************************************************/
			/*																		*/
			/*				Begin Digital ( call spread ) Inf Function				*/
			/*																		*/
			/************************************************************************/
class GPInfCallSpreadFctor : public ARMResultLong2LongFunc
{
public:
	GPInfCallSpreadFctor(	long	swapId,
							long	CF,
							double	strike,
							long	strikeId):	C_swapId	( swapId	),
												C_CF		( CF		),
												C_strike	( strike	),
 												C_strikeId	( strikeId	){};
	
	long operator()( ARM_result& result, long objId ) {
		return ARMLOCAL_GP_INFCALLSPREAD(	C_swapId,	C_CF,	C_strike,	C_strikeId,	result,	objId	);
	}

private:
	long	C_swapId;
	long	C_CF;
	double	C_strike;
    long	C_strikeId;
};


/////////////////////////////////////////////////////////////
/// central function that does the creation of the XL function
/////////////////////////////////////////////////////////////
LPXLOPER Local_ARM_GP_INFCALLSPREADCommon(
	LPXLOPER XL_swapId,
	LPXLOPER XL_CAPOrFLOOR,
	LPXLOPER XL_strike,
	bool PersistentInXL )
{
	// this is defined first because it is used in XL macros
	static XLOPER XL_result;
	ARM_result C_result;
	// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
		/// to avoid computation if called by the wizard
		ARM_NOCALCIFWIZ();
		
		/// Get the variables from the XLOper variables
		CCString swapStr;
		long swapId;

		CCString CF;
		long CFId;

		//strike
		double strike;
		CCString strikeStr;
		long     strikeId;

		// error
		static int error;
		static char* reason = "";
		
		XL_readStrOrNumCell(XL_strike, strikeStr, strike, strikeId,
			   " ARM_ERR: strike: numeric or srike Id expected",C_result);

		XL_readStrCell(XL_swapId,swapStr," ARM_ERR: swap id: object expected",C_result);
		XL_readStrCell(XL_CAPOrFLOOR,CF," ARM_ERR: Cap or Floor: string expected",C_result);

		swapId = LocalGetNumObjectId (swapStr);
		if((CFId = ARM_ConvCapOrFloor (CF, C_result)) == ARM_DEFAULT_ERR)
		{
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}
		
		if(strikeId == XL_TYPE_STRING)
			strikeId = LocalGetNumObjectId(strikeStr);
		else
			strikeId = ARM_NULL_OBJECT;

		
		GPInfCallSpreadFctor ourFunc(	swapId,
										CFId,
										strike,
 										strikeId);

		/// call the general function
		fillXL_Result( LOCAL_GP_INFCALLSPREAD_CLASS, ourFunc, C_result, XL_result, PersistentInXL );
	
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_ARM_GP_INFCALLSPREADCommon" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


//////////////////////////////////
/// version that takes into account 
/// previous creation of object
//////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_ARM_GP_INFCALLSPREAD(
	 LPXLOPER XL_swapId,
	 LPXLOPER XL_CAPOrFLOOR,
	 LPXLOPER XL_strike)
{
	ADD_LOG("Local_ARM_GP_INFCALLSPREAD");
	 bool PersistentInXL = true;
	 return Local_ARM_GP_INFCALLSPREADCommon(
		XL_swapId,
		XL_CAPOrFLOOR,
		XL_strike,
 		PersistentInXL);
}

///////////////////////////////////
/// version for non persistentInXl
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_PXL_ARM_GP_INFCALLSPREAD(
	 LPXLOPER XL_swapId,
	 LPXLOPER XL_CAPOrFLOOR,
	 LPXLOPER XL_strike)
{
	ADD_LOG("Local_PXL_ARM_GP_INFCALLSPREAD");
	 bool PersistentInXL = false;
	 return Local_ARM_GP_INFCALLSPREADCommon(
		XL_swapId,
		XL_CAPOrFLOOR,
		XL_strike,
 		PersistentInXL);
}

			/************************************************************************/
			/*						End Digital Inf Function						*/
			/************************************************************************/

__declspec(dllexport) LPXLOPER WINAPI Local_ARM_GetReset(LPXLOPER XL_resetManager,
															 LPXLOPER XL_date)
{
	ADD_LOG("Local_ARM_GetReset");
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();

	// C variable
	CCString C_resetMgr;
	double C_date;

	// error
	static int error;
	static char* reason = "";

	XL_readStrCell(XL_resetManager,C_resetMgr," ARM_ERR: resetMgr id: object expected",C_result);
	XL_readNumCell(XL_date,C_date," ARM_ERR: date: numeric expected",C_result);

	long retCode;

	if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 1)
	{
		retCode = ARMLOCAL_GetReset (LocalGetNumObjectId (C_resetMgr), C_date, C_result);
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_ARM_GetReset" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}

///////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////GP_INF_DIGITAL/////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////
/// very rapid definition of ourlonglongFunctor
///////////////////////////////////////////////
class GPInfDigitalFctor : public ARMResultLong2LongFunc
{
	
public:
	GPInfDigitalFctor(	
		long PaySwapId,
		long DigitSwapId,
		long PayOffType,
		double Barrier,
        long BarrierId,
		long CFId,
		long RecOrPayId)
	:
		C_PaySwapId( PaySwapId ),
		C_DigitSwapId( DigitSwapId ),
		C_PayOffType( PayOffType ),
 		C_Barrier( Barrier ),
		C_BarrierId( BarrierId ),
		C_CFId(CFId),
		C_RecOrPayId(RecOrPayId)
	{};
	
	long operator()( ARM_result& result, long objId ) 
	{
		return ARMLOCAL_GP_INFDIGITAL(
			C_PaySwapId,
			C_DigitSwapId,
			C_PayOffType,
 			C_Barrier,
			C_BarrierId,
			C_CFId,
			C_RecOrPayId,
			result,								
			objId);
	}

private:
	long C_PaySwapId;
	long C_DigitSwapId;
	long C_PayOffType;
	double C_Barrier;
    long C_BarrierId;
	long C_CFId;
	long C_RecOrPayId;
};



/////////////////////////////////////////////////////////////
/// central function that does the creation of the XL function
/////////////////////////////////////////////////////////////
LPXLOPER Local_ARM_GP_INFDIGITALCommon(
	LPXLOPER XL_PaySwapId,
    LPXLOPER XL_DigitSwapId,
    LPXLOPER XL_PayOffType,
    LPXLOPER XL_Barrier,
	LPXLOPER XL_CapOrFloor,
	LPXLOPER XL_PayOrRec,
	bool PersistentInXL )
{
	// this is defined first because it is used in XL macros
	static XLOPER XL_result;
	ARM_result C_result;
	// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
		/// to avoid computation if called by the wizard
		ARM_NOCALCIFWIZ();
		
		/// Get the variables from the XLOper variables
		CCString paylegStr;
		long paylegId;

		CCString digitlegStr;
		long digitlegId;

		CCString payoffType;
		long payoffTypeId;

		//strike
		double barrier;
		CCString barrierStr;
		long     barrierId;

		CCString CF;
		long CFId;

		CCString RecOrPay;
		long RecOrPayId;

		// error
		static int error;
		static char* reason = "";
		
		XL_readStrOrNumCell(XL_Barrier, barrierStr, barrier, barrierId,
			   " ARM_ERR: barrier: numeric or barrier Id expected",C_result);

		XL_readStrCell(XL_PaySwapId,paylegStr," ARM_ERR: payLeg id: object expected",C_result);
		XL_readStrCell(XL_DigitSwapId,digitlegStr," ARM_ERR: digitleg id: object expected",C_result);

		XL_readStrCell(XL_PayOffType,payoffType," ARM_ERR: PayOff type: string expected",C_result);

		paylegId = LocalGetNumObjectId (paylegStr);
		digitlegId = LocalGetNumObjectId (digitlegStr);

		if((payoffTypeId = ARM_ConvINFDigitalPayoffType (payoffType, C_result)) == ARM_DEFAULT_ERR)
		{
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}
		
		if(barrierId == XL_TYPE_STRING)
			barrierId = LocalGetNumObjectId(barrierStr);
		else
			barrierId = ARM_NULL_OBJECT;


		XL_readStrCell(XL_CapOrFloor,CF," ARM_ERR: Cap or Floor: string expected",C_result);

		if((CFId = ARM_ConvCapOrFloor (CF, C_result)) == ARM_DEFAULT_ERR)
		{
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}

		XL_readStrCell(XL_PayOrRec,RecOrPay," ARM_ERR: Rec or Pay: string expected",C_result);

		if((RecOrPayId = ARM_ConvRecOrPay (RecOrPay, C_result)) == ARM_DEFAULT_ERR)
		{
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}

		GPInfDigitalFctor ourFunc(
				paylegId,
				digitlegId,
				payoffTypeId,
				barrier,
 				barrierId,
				CFId,
				RecOrPayId);

		/// call the general function
		fillXL_Result( LOCAL_GP_INFCALLSPREAD_CLASS, ourFunc, C_result, XL_result, PersistentInXL );
	
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_ARM_GP_INFDIGITALCommon" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


//////////////////////////////////
/// version that takes into account 
/// previous creation of object
//////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_ARM_GP_INFDIGITAL(
	LPXLOPER XL_PaySwapId,
    LPXLOPER XL_DigitSwapId,
    LPXLOPER XL_PayOffType,
    LPXLOPER XL_Barrier,
	LPXLOPER XL_CapOrFloor,
	LPXLOPER XL_PayOrRec)
{
	ADD_LOG("Local_ARM_GP_INFDIGITAL");
	 bool PersistentInXL = true;
	 return Local_ARM_GP_INFDIGITALCommon(
		XL_PaySwapId,
		XL_DigitSwapId,
		XL_PayOffType,
		XL_Barrier,
		XL_CapOrFloor,
		XL_PayOrRec,
 		PersistentInXL);
}

///////////////////////////////////
/// version for non persistentInXl
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_PXL_ARM_GP_INFDIGITAL(
	 LPXLOPER XL_PaySwapId,
     LPXLOPER XL_DigitSwapId,
     LPXLOPER XL_PayOffType,
     LPXLOPER XL_Barrier,
	 LPXLOPER XL_CapOrFloor,
	 LPXLOPER XL_PayOrRec)
{
	ADD_LOG("Local_PXL_ARM_GP_INFDIGITAL");
	 bool PersistentInXL = false;
	 return Local_ARM_GP_INFDIGITALCommon(
		XL_PaySwapId,
		XL_DigitSwapId,
		XL_PayOffType,
		XL_Barrier,
		XL_CapOrFloor,
		XL_PayOrRec,
 		PersistentInXL);
}


//////////////////////////////////////////////////////////////////////
/// central function that does the creation of the Hybrid Inf Ir Model
//////////////////////////////////////////////////////////////////////
LPXLOPER Local_ARM_INF_HybridInfIrMkt_CreateCommon(
	LPXLOPER XL_AsOf,
	LPXLOPER XL_Keys,
	LPXLOPER XL_ModId,
	bool PersistentInXL ){

	static XLOPER XL_result;
	ARM_result C_result;
	ARM_XL_TRY_BLOCK_BEGIN	{
	
		ARM_NOCALCIFWIZ();
		static int error;
		static char* reason = "";

		double C_asOfDble;
		XL_readNumCell(	 XL_AsOf,	C_asOfDble,	" ARM_ERR: asOfDate: date expected",	C_result);
		ARM_Date C_asOf = ConvertToARMDATE(C_asOfDble);
  
		VECTOR<CCString> Keys;
		vector<string>	 KeysSTL;
		XL_readStrVector (XL_Keys,Keys," ARM_ERR: Market datas keys : array of string expected",DOUBLE_TYPE,C_result);

		VECTOR<CCString> ModIdStr;
		vector<long>	 ModId; 
		XL_readStrVector (XL_ModId,ModIdStr," ARM_ERR: Market datas keys : array of string expected",DOUBLE_TYPE,C_result);

		KeysSTL.resize(Keys.size());
		ModId.resize(ModIdStr.size());
		for (int i=0; i < ModIdStr.size(); ++i)	{	
			KeysSTL[i] = CCSTringToSTLString(Keys[i]);
			ModId[i] = LocalGetNumObjectId(ModIdStr[i]);
		}

		exportFunc3Args<ARM_Date, vector<string>,vector<long> > ourFunc(
			C_asOf,
			KeysSTL,
			ModId,
			ARMLOCAL_HybridInfIrMkt_Create);

		fillXL_Result( LOCAL_INFHYBRIDMKT_CLASS, ourFunc, C_result, XL_result, PersistentInXL );

	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_MktDataManager_MktDataGetSet_Common" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


//////////////////////////////////
/// version that takes into account 
/// previous creation of object
//////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_ARM_INF_HybridInfIrMkt_Create(
	LPXLOPER XL_Asof,
	LPXLOPER XL_Keys,
	LPXLOPER XL_ModId)
{
	ADD_LOG("Local_ARM_INF_HybridInfIrMkt_Create");
	 bool PersistentInXL = true;
	 return Local_ARM_INF_HybridInfIrMkt_CreateCommon(
		XL_Asof,
		XL_Keys,
		XL_ModId,
 		PersistentInXL);
}

///////////////////////////////////
/// version for non persistentInXl
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_PXL_ARM_INF_HybridInfIrMkt_Create(
	LPXLOPER XL_Asof,
	LPXLOPER XL_Keys,
	LPXLOPER XL_ModId)
{
	ADD_LOG("Local_PXL_ARM_INF_HybridInfIrMkt_Create");
	 bool PersistentInXL = true;
	 return Local_ARM_INF_HybridInfIrMkt_CreateCommon(
		XL_Asof,
		XL_Keys,
		XL_ModId,
 		PersistentInXL);
}

////////////////////////////////////////////////
/// general function to branch
////////////////////////////////////////////////


class InfBSSmiledModelFunc : public ARMResultLong2LongFunc
{
public:
	InfBSSmiledModelFunc(	double	asOfDate, 
							long	discountCurvId, 
							long	infFwdCurvId, 
							long	volSigmaId,
							long	volNuId,
							long	volRhoId,
							long	volBetaId,
							long	volAtmIrId,
							long	correlId,
							long	correlAdjId):
		C_AsOfDate			( asOfDate			), 
		C_DiscountCurvId	( discountCurvId	),
		C_InfFwdCurvId		( infFwdCurvId		), 
		C_VolSigmaId		( volSigmaId		),
		C_VolNuId			( volNuId			),
		C_VolRhoId			( volRhoId			),
		C_VolBetaId			( volBetaId			),
		C_VolAtmIrId		( volAtmIrId		),
		C_CorrelId			( correlId			),
		C_CorrelAdjId		( correlAdjId		){};

	long operator()( ARM_result& result, long objId ) {

		return ARMLOCAL_infBSSmiledModel( 	C_AsOfDate, 
											C_DiscountCurvId,
											C_InfFwdCurvId, 
											C_VolSigmaId,
											C_VolNuId,
											C_VolRhoId,
											C_VolBetaId,
											C_VolAtmIrId,
											C_CorrelId,
											C_CorrelAdjId,
											result, 
											objId );
	}
private:
	double	C_AsOfDate;
	long	C_DiscountCurvId;
	long	C_InfFwdCurvId;
	long	C_VolSigmaId;
	long	C_VolNuId;
	long	C_VolRhoId;
	long	C_VolBetaId;
	long	C_VolAtmIrId;
	long	C_CorrelId;
	long	C_CorrelAdjId;
};

LPXLOPER Local_InfBSSmiledMod_Common(
	LPXLOPER XL_AsOfDate,
	LPXLOPER XL_DiscountCurvId,
	LPXLOPER XL_InfFwdCurvId,
	LPXLOPER XL_VolSigmaId,
	LPXLOPER XL_VolNuId,
	LPXLOPER XL_VolRhoId,
	LPXLOPER XL_VolBetaId,
	LPXLOPER XL_VolAtmIrId,
	LPXLOPER XL_CorrelId,
	LPXLOPER XL_CorrelAdjId,
	bool PersistentInXL )
{
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

		double	C_AsOfDate;
		long	C_DiscountCurvId;
		long	C_InfFwdCurvId;
		long	C_VolSigmaId;
		long	C_VolNuId;
		long	C_VolRhoId;
		long	C_VolBetaId;
		long	C_VolAtmIrId;
		long	C_CorrelId;
		long	C_CorrelAdjId;
		
		XL_readNumCell(	XL_AsOfDate,		C_AsOfDate,							" ARM_ERR: asOfDate: date expected",			C_result);
		XL_GETOBJIDWD(	XL_DiscountCurvId,	C_DiscountCurvId,	"NULL OBJECT",	" ARM_ERR: Discount Curv: Object expected",		C_result);
		XL_GETOBJIDWD(	XL_InfFwdCurvId,	C_InfFwdCurvId ,	"NULL OBJECT",	" ARM_ERR: Inflation Fwd Curv: Object expected",C_result);
		XL_GETOBJIDWD(	XL_VolSigmaId,		C_VolSigmaId ,		"NULL OBJECT",	" ARM_ERR: Volatility curve: Object expected",	C_result);
		XL_GETOBJIDWD(	XL_VolNuId,			C_VolNuId ,			"NULL OBJECT",	" ARM_ERR: Volatility curve: Object expected",	C_result);
		XL_GETOBJIDWD(	XL_VolRhoId,		C_VolRhoId ,		"NULL OBJECT",	" ARM_ERR: Volatility curve: Object expected",	C_result);
		XL_GETOBJIDWD(	XL_VolBetaId,		C_VolBetaId ,		"NULL OBJECT",	" ARM_ERR: Volatility curve: Object expected",	C_result);
		XL_GETOBJIDWD(	XL_VolAtmIrId,		C_VolAtmIrId ,		"NULL OBJECT",	" ARM_ERR: Volatility curve: Object expected",	C_result);
		XL_GETOBJIDWD(	XL_CorrelId,		C_CorrelId ,		"NULL OBJECT",	" ARM_ERR: Volatility curve: Object expected",	C_result);
		XL_GETOBJIDWD(	XL_CorrelAdjId,		C_CorrelAdjId ,		"NULL OBJECT",	" ARM_ERR: Volatility curve: Object expected",	C_result);

		InfBSSmiledModelFunc ourFunc( C_AsOfDate, C_DiscountCurvId, C_InfFwdCurvId, C_VolSigmaId,	C_VolNuId, C_VolRhoId,	C_VolBetaId, C_VolAtmIrId, C_CorrelId, C_CorrelAdjId );

		/// call the general function
		fillXL_Result( LOCAL_INFCURVMODEL_CLASS, ourFunc, C_result, XL_result, PersistentInXL );

	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_InfBSMod_Common" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


////////////////////////////////////////////////
/// Version for XL exportation of the InfYCMod
////////////////////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_InfBSSmiledModel(
	LPXLOPER XL_AsOfDate,
	LPXLOPER XL_DiscountCurvId,
	LPXLOPER XL_InfFwdCurvId,
	LPXLOPER XL_VolSigmaId,
	LPXLOPER XL_VolNuId,
	LPXLOPER XL_VolRhoId,
	LPXLOPER XL_VolBetaId,
	LPXLOPER XL_VolAtmIrId,
	LPXLOPER XL_CorrelId,
	LPXLOPER XL_CorrelAdjId	)
{
	ADD_LOG("Local_InfBSSmiledModel");
	bool PersistentInXL = true;
	return Local_InfBSSmiledMod_Common(
		XL_AsOfDate,
		XL_DiscountCurvId,
		XL_InfFwdCurvId,
		XL_VolSigmaId,
		XL_VolNuId,
		XL_VolRhoId,
		XL_VolBetaId,
		XL_VolAtmIrId,
		XL_CorrelId,
		XL_CorrelAdjId,
		PersistentInXL );
}



////////////////////////////////////////////////
/// Version for XL exportation 
////////////////////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_PXL_InfBSSmiledModel(
	LPXLOPER XL_AsOfDate,
	LPXLOPER XL_DiscountCurvId,
	LPXLOPER XL_InfFwdCurvId,
	LPXLOPER XL_VolSigmaId,
	LPXLOPER XL_VolNuId,
	LPXLOPER XL_VolRhoId,
	LPXLOPER XL_VolBetaId,
	LPXLOPER XL_VolAtmIrId,
	LPXLOPER XL_CorrelId,
	LPXLOPER XL_CorrelAdjId	)
{
	ADD_LOG("Local_PXL_InfBSSmiledModel");
	bool PersistentInXL = false;
	return Local_InfBSSmiledMod_Common(
		XL_AsOfDate,
		XL_DiscountCurvId,
		XL_InfFwdCurvId,
		XL_VolSigmaId,
		XL_VolNuId,
		XL_VolRhoId,
		XL_VolBetaId,
		XL_VolAtmIrId,
		XL_CorrelId,
		XL_CorrelAdjId,
		PersistentInXL );
}

////////////////////////////////////////////////
/// Version for Hybrid Inf Ir Common create
////////////////////////////////////////////////



LPXLOPER Local_ARM_INF_HybridInfIr_LoadCommon(
	LPXLOPER XL_InsId,		
	LPXLOPER XL_MktId,
	LPXLOPER XL_ModId,
	LPXLOPER XL_PayId,
	bool PersistentInXL ){

	static XLOPER XL_result;
	ARM_result C_result;
	ARM_XL_TRY_BLOCK_BEGIN	{
	
	ARM_NOCALCIFWIZ();
	static int error;
	static char* reason = "";

	long C_InsId;			// instrument
	long C_MktId;			// market Datas
	long C_ModId;			// model
	long C_PayId;			// Payoff
			
	XL_GETOBJIDWD(	XL_InsId,	C_InsId,	"NULL OBJECT",	" ARM_ERR: instrument: Object expected",	C_result);
	XL_GETOBJIDWD(	XL_MktId,	C_MktId ,	"NULL OBJECT",	" ARM_ERR: market datas: Object expected",	C_result);
	XL_GETOBJIDWD(	XL_ModId,	C_ModId ,	"NULL OBJECT",	" ARM_ERR: model: Object expected",			C_result);
	XL_GETOBJIDWD(	XL_PayId,	C_PayId ,	"NULL OBJECT",	" ARM_ERR: payoff: Object expected",		C_result);

	exportFunc4Args<long, long, long, long > ourFunc(
			C_InsId,
			C_MktId,
			C_ModId,
			C_PayId,
			ARMLOCAL_HybridInfIr_Load);

		fillXL_Result( LOCAL_INFHYBRIDMKT_CLASS, ourFunc, C_result, XL_result, PersistentInXL );

	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_MktDataManager_MktDataGetSet_Common" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


//////////////////////////////////
/// version that takes into account 
/// previous creation of object
//////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_ARM_INF_HybridInfIr_Load(
	LPXLOPER XL_Ins,		
	LPXLOPER XL_Mkt,
	LPXLOPER XL_Mod,
	LPXLOPER XL_Pay)
{
	ADD_LOG("Local_ARM_INF_HybridInfIr_Load");
	 bool PersistentInXL = true;
	 return Local_ARM_INF_HybridInfIr_LoadCommon(
		XL_Ins,		
		XL_Mkt,
		XL_Mod,
		XL_Pay,
 		PersistentInXL);
}

__declspec(dllexport) LPXLOPER WINAPI Local_PXL_ARM_INF_HybridInfIr_Load(
	LPXLOPER XL_Ins,		
	LPXLOPER XL_Mkt,
	LPXLOPER XL_Mod,
	LPXLOPER XL_Pay)
{
	ADD_LOG("Local_PXL_ARM_INF_HybridInfIr_Load");
	 bool PersistentInXL = false;
	 return Local_ARM_INF_HybridInfIr_LoadCommon(
		XL_Ins,		
		XL_Mkt,
		XL_Mod,
		XL_Pay,
 		PersistentInXL);
}

////////////*********************************//////////////

//////////////////////////////////////////////////////////////////////
/// central function that does the creation of the Hybrid Inf Ir PayOff
//////////////////////////////////////////////////////////////////////
LPXLOPER Local_ARM_INF_HybridInfIrPayOff_CreateCommon(
		LPXLOPER XL_CstCpnCoef,
		LPXLOPER XL_CstOptCoef,

		LPXLOPER XL_MainCpnName,
		LPXLOPER XL_MainCpnCoef,
		LPXLOPER XL_MainOptName,
		LPXLOPER XL_MainOptCoef,

		LPXLOPER XL_SubCpnName,
		LPXLOPER XL_SubCpnCoef,
		LPXLOPER XL_SubOptName,
		LPXLOPER XL_SubOptCoef,

		LPXLOPER XL_SupCpnName,
		LPXLOPER XL_SupCpnCoef,
		LPXLOPER XL_SupOptName,
		LPXLOPER XL_SupOptCoef,

		bool PersistentInXL ){

	static XLOPER XL_result;
	ARM_result C_result;
	ARM_XL_TRY_BLOCK_BEGIN	{
	
		ARM_NOCALCIFWIZ();
		static int error;
		static char* reason = "";

		CCString C_MainCpnName;
		CCString C_SubCpnName;
		CCString C_SupCpnName;

		CCString C_MainOptName;
		CCString C_SubOptName;
		CCString C_SupOptName;

		long C_MainCpnCoef;
		long C_SubCpnCoef;
		long C_SupCpnCoef;
		long C_CstCpnCoef;

		long C_MainOptCoef;
		long C_SubOptCoef;
		long C_SupOptCoef;
		long C_CstOptCoef;

		XL_readStrCellWD(	XL_MainCpnName,		C_MainCpnName,		"NO",	" ARM_ERR: Main Cpn: string expected",			C_result);
		XL_readStrCellWD(	XL_SubCpnName,		C_SubCpnName,		"NO",	" ARM_ERR: Sub  Cpn: string expected",			C_result);
		XL_readStrCellWD(	XL_SupCpnName,		C_SupCpnName,		"NO",	" ARM_ERR: Sup  Cpn: string expected",			C_result);

		XL_readStrCellWD(	XL_MainOptName,		C_MainOptName,		"NO",	" ARM_ERR: Main Opt: string expected",			C_result);
		XL_readStrCellWD(	XL_SubOptName,		C_SubOptName,		"NO",	" ARM_ERR: Sub  Opt: string expected",			C_result);
		XL_readStrCellWD(	XL_SupOptName,		C_SupOptName,		"NO",	" ARM_ERR: Sup  Opt: string expected",			C_result);

		string	MainCpnName	=	CCSTringToSTLString( C_MainCpnName);
		string	SubCpnName	=	CCSTringToSTLString( C_SubCpnName );
		string	SupCpnName	=	CCSTringToSTLString( C_SupCpnName );
		string	MainOptName	=	CCSTringToSTLString( C_MainOptName);
		string	SubOptName	=	CCSTringToSTLString( C_SubOptName );
		string	SupOptName	=	CCSTringToSTLString( C_SupOptName );
		

		XL_GETOBJIDWD	(	XL_MainCpnCoef,		C_MainCpnCoef,		"NULL OBJECT",	" ARM_ERR: Main Cpn: Object expected",	C_result);
		XL_GETOBJIDWD	(	XL_SubCpnCoef,		C_SubCpnCoef ,		"NULL OBJECT",	" ARM_ERR: Sub  Cpn: Object expected",	C_result);
		XL_GETOBJIDWD	(	XL_SupCpnCoef,		C_SupCpnCoef ,		"NULL OBJECT",	" ARM_ERR: Sup  Cpn: Object expected",	C_result);
		XL_GETOBJIDWD	(	XL_CstCpnCoef,		C_CstCpnCoef ,		"NULL OBJECT",	" ARM_ERR: Cst  Cpn: Object expected",	C_result);

		XL_GETOBJIDWD	(	XL_MainOptCoef,		C_MainOptCoef,		"NULL OBJECT",	" ARM_ERR: Main Opt: Object expected",	C_result);
		XL_GETOBJIDWD	(	XL_SubOptCoef,		C_SubOptCoef ,		"NULL OBJECT",	" ARM_ERR: Sub  Opt: Object expected",	C_result);
		XL_GETOBJIDWD	(	XL_SupOptCoef,		C_SupOptCoef ,		"NULL OBJECT",	" ARM_ERR: Sup  Opt: Object expected",	C_result);
		XL_GETOBJIDWD	(	XL_CstOptCoef,		C_CstOptCoef ,		"NULL OBJECT",	" ARM_ERR: Cst  Opt: Object expected",	C_result);


		exportFunc14Args<	long,
							long,
							string, long, 			 
							string, long,
							string, long,
							string, long,
							string, long,
							string, long	> ourFunc(
			C_CstCpnCoef,
			C_CstOptCoef,
			MainCpnName,
			C_MainCpnCoef,
			MainOptName,
			C_MainOptCoef,

			SubCpnName,
			C_SubCpnCoef,
			SubOptName,
			C_SubOptCoef,

			SupCpnName,
			C_SupCpnCoef,
			SupOptName,
			C_SupOptCoef,

			ARMLOCAL_HybridInfIrPayOff_Create);

		fillXL_Result( LOCAL_INFHYBRIDPAYOFF_CLASS, ourFunc, C_result, XL_result, PersistentInXL );

	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_MktDataManager_MktDataGetSet_Common" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


//////////////////////////////////
/// version that takes into account 
/// previous creation of object
//////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_ARM_INF_HybridInfIrPayOff_Create(
		LPXLOPER XL_CstCpnCoef,
		LPXLOPER XL_CstOptCoef,

		LPXLOPER XL_MainCpnName,
		LPXLOPER XL_MainCpnCoef,
		LPXLOPER XL_MainOptName,
		LPXLOPER XL_MainOptCoef,

		LPXLOPER XL_SubCpnName,
		LPXLOPER XL_SubCpnCoef,
		LPXLOPER XL_SubOptName,
		LPXLOPER XL_SubOptCoef,

		LPXLOPER XL_SupCpnName,
		LPXLOPER XL_SupCpnCoef,
		LPXLOPER XL_SupOptName,
		LPXLOPER XL_SupOptCoef	)
{
	bool PersistentInXL = true;
	return Local_ARM_INF_HybridInfIrPayOff_CreateCommon(
		XL_CstCpnCoef,
		XL_CstOptCoef,

		XL_MainCpnName,
		XL_MainCpnCoef,
		XL_MainOptName,
		XL_MainOptCoef,

		XL_SubCpnName,
		XL_SubCpnCoef,
		XL_SubOptName,
		XL_SubOptCoef,

		XL_SupCpnName,
		XL_SupCpnCoef,
		XL_SupOptName,
		XL_SupOptCoef,
		
		PersistentInXL);
}

///////////////////////////////////
/// version for non persistentInXl
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_PXL_ARM_INF_HybridInfIrPayOff_Create(
		LPXLOPER XL_CstCpnCoef,
		LPXLOPER XL_CstOptCoef,

		LPXLOPER XL_MainCpnName,
		LPXLOPER XL_MainCpnCoef,
		LPXLOPER XL_MainOptName,
		LPXLOPER XL_MainOptCoef,

		LPXLOPER XL_SubCpnName,
		LPXLOPER XL_SubCpnCoef,
		LPXLOPER XL_SubOptName,
		LPXLOPER XL_SubOptCoef,

		LPXLOPER XL_SupCpnName,
		LPXLOPER XL_SupCpnCoef,
		LPXLOPER XL_SupOptName,
		LPXLOPER XL_SupOptCoef)
{
	 bool PersistentInXL = true;
	 return Local_ARM_INF_HybridInfIrPayOff_CreateCommon(
		XL_CstCpnCoef,
		XL_CstOptCoef,

		XL_MainCpnName,
		XL_MainCpnCoef,
		XL_MainOptName,
		XL_MainOptCoef,

		XL_SubCpnName,
		XL_SubCpnCoef,
		XL_SubOptName,
		XL_SubOptCoef,

		XL_SupCpnName,
		XL_SupCpnCoef,
		XL_SupOptName,
		XL_SupOptCoef,

		PersistentInXL);
}


//////////////////////////////////////////////////////////////////////
/// central function that does the creation of the Hybrid Inf Ir PayOff
//////////////////////////////////////////////////////////////////////
LPXLOPER Local_ARM_INF_SpreadCap_CreateCommon(
		LPXLOPER XL_MainIndex,
		LPXLOPER XL_MainType,
		LPXLOPER XL_MainLeverage,
		LPXLOPER XL_SubIndex,
		LPXLOPER XL_SubType,
		LPXLOPER XL_SubLeverage,
		LPXLOPER XL_Strike,
		LPXLOPER XL_Notional,
		bool PersistentInXL ){

	static XLOPER XL_result;
	ARM_result C_result;
	ARM_XL_TRY_BLOCK_BEGIN	{
	
		ARM_NOCALCIFWIZ();
		static int error;
		static char* reason = "";

		CCString C_MainIndex;
		CCString C_SubIndex;
		CCString C_MainType;
		CCString C_SubType;

		long MainLeverage;
		long SubLeverage;
		long Strike;
		long Notional;
	
		XL_readStrCellWD(	XL_MainIndex,		C_MainIndex,		"NO",	" ARM_ERR: Main Index:	string expected",			C_result);
		XL_readStrCellWD(	XL_SubIndex,		C_SubIndex,			"NO",	" ARM_ERR: Sub  Index:	string expected",			C_result);
		XL_readStrCellWD(	XL_MainType,		C_MainType,			"IR",	" ARM_ERR: Main Type:	string expected",			C_result);
		XL_readStrCellWD(	XL_SubType,			C_SubType,			"IR",	" ARM_ERR: Sub  Type:	string expected",			C_result);

		string	MainIndex	=	CCSTringToSTLString( C_MainIndex);
		string	SubIndex	=	CCSTringToSTLString( C_SubIndex );
		string	MainType	=	CCSTringToSTLString( C_MainType);
		string	SubType		=	CCSTringToSTLString( C_SubType );

		XL_GETOBJIDWD	(	XL_MainLeverage,	MainLeverage,		"NULL OBJECT",	" ARM_ERR: Main Leverage: Object expected",	C_result);
		XL_GETOBJIDWD	(	XL_SubLeverage,		SubLeverage ,		"NULL OBJECT",	" ARM_ERR: Sub  Leverage: Object expected",	C_result);
		XL_GETOBJIDWD	(	XL_Strike,			Strike ,			"NULL OBJECT",	" ARM_ERR: Strike		: Object expected",	C_result);
		XL_GETOBJIDWD	(	XL_Notional,		Notional ,			"NULL OBJECT",	" ARM_ERR: Notional		: Object expected",	C_result);

		exportFunc8Args<	string, string, long, 			 
							string, string, long,
							long,
							long> ourFunc(
			MainIndex,
			MainType,
			MainLeverage,
			SubIndex,
			SubType,
			SubLeverage,
			Strike,
			Notional,
			ARMLOCAL_InfSpreadCap_Create);

		fillXL_Result( LOCAL_INFHYBRIDPAYOFF_CLASS, ourFunc, C_result, XL_result, PersistentInXL );

	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_MktDataManager_MktDataGetSet_Common" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


//////////////////////////////////
/// version that takes into account 
/// previous creation of object
//////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_ARM_INF_SpreadCap_Create(
		LPXLOPER XL_MainIndex,
		LPXLOPER XL_MainType,
		LPXLOPER XL_MainLeverage,
		LPXLOPER XL_SubIndex,
		LPXLOPER XL_SubType,
		LPXLOPER XL_SubLeverage,
		LPXLOPER XL_Strike,
		LPXLOPER XL_Notional )
{
	bool PersistentInXL = true;
	return Local_ARM_INF_SpreadCap_CreateCommon(
		XL_MainIndex,
		XL_MainType,
		XL_MainLeverage,
		XL_SubIndex,
		XL_SubType,
		XL_SubLeverage,
		XL_Strike,
		XL_Notional,
		PersistentInXL);
}

///////////////////////////////////
/// version for non persistentInXl
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_PXL_ARM_INF_SpreadCap_Create(
		LPXLOPER XL_MainIndex,
		LPXLOPER XL_MainType,
		LPXLOPER XL_MainLeverage,
		LPXLOPER XL_SubIndex,
		LPXLOPER XL_SubType,
		LPXLOPER XL_SubLeverage,
		LPXLOPER XL_Strike,
		LPXLOPER XL_Notional)
{
	 bool PersistentInXL = true;
	 return Local_ARM_INF_SpreadCap_CreateCommon(
		XL_MainIndex,
		XL_MainType,
		XL_MainLeverage,
		XL_SubIndex,
		XL_SubType,
		XL_SubLeverage,
		XL_Strike,
		XL_Notional,
		PersistentInXL);
}


//
//////////////////////////////////////////////////////////////////////
/// central function that does the creation of the Hybrid Inf Ir PayOff
//////////////////////////////////////////////////////////////////////
LPXLOPER Local_ARM_INF_SpreadDigital_CreateCommon(
		LPXLOPER XL_MainIndex,
		LPXLOPER XL_MainType,
		LPXLOPER XL_MainLeverage,
		LPXLOPER XL_SubIndex,
		LPXLOPER XL_SubType,
		LPXLOPER XL_SubLeverage,
		LPXLOPER XL_Strike,
		LPXLOPER XL_Notional,
		bool PersistentInXL ){

	static XLOPER XL_result;
	ARM_result C_result;
	ARM_XL_TRY_BLOCK_BEGIN	{
	
		ARM_NOCALCIFWIZ();
		static int error;
		static char* reason = "";

		CCString C_MainIndex;
		CCString C_SubIndex;
		CCString C_MainType;
		CCString C_SubType;

		long MainLeverage;
		long SubLeverage;
		long Strike;
		long Notional;
	
		XL_readStrCellWD(	XL_MainIndex,		C_MainIndex,		"NO",	" ARM_ERR: Main Index:	string expected",			C_result);
		XL_readStrCellWD(	XL_SubIndex,		C_SubIndex,			"NO",	" ARM_ERR: Sub  Index:	string expected",			C_result);
		XL_readStrCellWD(	XL_MainType,		C_MainType,			"IR",	" ARM_ERR: Main Type:	string expected",			C_result);
		XL_readStrCellWD(	XL_SubType,			C_SubType,			"IR",	" ARM_ERR: Sub  Type:	string expected",			C_result);

		string	MainIndex	=	CCSTringToSTLString( C_MainIndex);
		string	SubIndex	=	CCSTringToSTLString( C_SubIndex );
		string	MainType	=	CCSTringToSTLString( C_MainType);
		string	SubType		=	CCSTringToSTLString( C_SubType );

		XL_GETOBJIDWD	(	XL_MainLeverage,	MainLeverage,		"NULL OBJECT",	" ARM_ERR: Main Leverage: Object expected",	C_result);
		XL_GETOBJIDWD	(	XL_SubLeverage,		SubLeverage ,		"NULL OBJECT",	" ARM_ERR: Sub  Leverage: Object expected",	C_result);
		XL_GETOBJIDWD	(	XL_Strike,			Strike ,			"NULL OBJECT",	" ARM_ERR: Strike		: Object expected",	C_result);
		XL_GETOBJIDWD	(	XL_Notional,		Notional ,			"NULL OBJECT",	" ARM_ERR: Notional		: Object expected",	C_result);

		exportFunc8Args<	string, string, long, 			 
							string, string, long,
							long,
							long> ourFunc(
			MainIndex,
			MainType,
			MainLeverage,
			SubIndex,
			SubType,
			SubLeverage,
			Strike,
			Notional,
			ARMLOCAL_InfSpreadDigital_Create);

		fillXL_Result( LOCAL_INFHYBRIDPAYOFF_CLASS, ourFunc, C_result, XL_result, PersistentInXL );

	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_MktDataManager_MktDataGetSet_Common" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


//////////////////////////////////
/// version that takes into account 
/// previous creation of object
//////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_ARM_INF_SpreadDigital_Create(
		LPXLOPER XL_MainIndex,
		LPXLOPER XL_MainType,
		LPXLOPER XL_MainLeverage,
		LPXLOPER XL_SubIndex,
		LPXLOPER XL_SubType,
		LPXLOPER XL_SubLeverage,
		LPXLOPER XL_Strike,
		LPXLOPER XL_Notional )
{
	bool PersistentInXL = true;
	return Local_ARM_INF_SpreadDigital_CreateCommon(
		XL_MainIndex,
		XL_MainType,
		XL_MainLeverage,
		XL_SubIndex,
		XL_SubType,
		XL_SubLeverage,
		XL_Strike,
		XL_Notional,
		PersistentInXL);
}

///////////////////////////////////
/// version for non persistentInXl
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_PXL_ARM_INF_SpreadDigital_Create(
		LPXLOPER XL_MainIndex,
		LPXLOPER XL_MainType,
		LPXLOPER XL_MainLeverage,
		LPXLOPER XL_SubIndex,
		LPXLOPER XL_SubType,
		LPXLOPER XL_SubLeverage,
		LPXLOPER XL_Strike,
		LPXLOPER XL_Notional)
{
	 bool PersistentInXL = true;
	 return Local_ARM_INF_SpreadDigital_CreateCommon(
		XL_MainIndex,
		XL_MainType,
		XL_MainLeverage,
		XL_SubIndex,
		XL_SubType,
		XL_SubLeverage,
		XL_Strike,
		XL_Notional,
		PersistentInXL);
}

//
//////////////////////////////////////////////////////////////////////
/// central function that does the creation of the Hybrid Inf Ir PayOff
//////////////////////////////////////////////////////////////////////
LPXLOPER Local_ARM_INF_DoubleDigital_CreateCommon(
		LPXLOPER XL_MainIndex,
		LPXLOPER XL_MainType,
		LPXLOPER XL_MainLeverage,
		LPXLOPER XL_MainStrike,
		LPXLOPER XL_SubIndex,
		LPXLOPER XL_SubType,
		LPXLOPER XL_SubLeverage,
		LPXLOPER XL_SubStrike,
		LPXLOPER XL_Spread,
		LPXLOPER XL_Notional,
		bool PersistentInXL ){

	static XLOPER XL_result;
	ARM_result C_result;
	ARM_XL_TRY_BLOCK_BEGIN	{
	
		ARM_NOCALCIFWIZ();
		static int error;
		static char* reason = "";

		CCString C_MainIndex;
		CCString C_SubIndex;
		CCString C_MainType;
		CCString C_SubType;

		long MainLeverage;
		long SubLeverage;
		long MainStrike;
		long SubStrike;
		long Spread;
		long Notional;
	
		XL_readStrCellWD(	XL_MainIndex,		C_MainIndex,		"NO",	" ARM_ERR: Main Index:	string expected",			C_result);
		XL_readStrCellWD(	XL_SubIndex,		C_SubIndex,			"NO",	" ARM_ERR: Sub  Index:	string expected",			C_result);
		XL_readStrCellWD(	XL_MainType,		C_MainType,			"IR",	" ARM_ERR: Main Type:	string expected",			C_result);
		XL_readStrCellWD(	XL_SubType,			C_SubType,			"IR",	" ARM_ERR: Sub  Type:	string expected",			C_result);

		string	MainIndex	=	CCSTringToSTLString( C_MainIndex);
		string	SubIndex	=	CCSTringToSTLString( C_SubIndex );
		string	MainType	=	CCSTringToSTLString( C_MainType);
		string	SubType		=	CCSTringToSTLString( C_SubType );

		XL_GETOBJIDWD	(	XL_MainLeverage,	MainLeverage,		"NULL OBJECT",	" ARM_ERR: Main Leverage: Object expected",	C_result);
		XL_GETOBJIDWD	(	XL_SubLeverage,		SubLeverage ,		"NULL OBJECT",	" ARM_ERR: Sub  Leverage: Object expected",	C_result);
		XL_GETOBJIDWD	(	XL_MainStrike,		MainStrike ,		"NULL OBJECT",	" ARM_ERR: MainStrike	: Object expected",	C_result);
		XL_GETOBJIDWD	(	XL_SubStrike,		SubStrike ,			"NULL OBJECT",	" ARM_ERR: SubStrike	: Object expected",	C_result);
		XL_GETOBJIDWD	(	XL_Spread,			Spread ,			"NULL OBJECT",	" ARM_ERR: Spread		: Object expected",	C_result);
		XL_GETOBJIDWD	(	XL_Notional,		Notional ,			"NULL OBJECT",	" ARM_ERR: Notional		: Object expected",	C_result);

		exportFunc10Args<	string, string, long, long,			 
							string, string, long, long,
							long, long> ourFunc(
			MainIndex,
			MainType,
			MainLeverage,
			MainStrike,
			SubIndex,
			SubType,
			SubLeverage,
			SubStrike,
			Spread,
			Notional,
			ARMLOCAL_InfDoubleDigital_Create);

		fillXL_Result( LOCAL_INFHYBRIDPAYOFF_CLASS, ourFunc, C_result, XL_result, PersistentInXL );

	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_MktDataManager_MktDataGetSet_Common" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


//////////////////////////////////
/// version that takes into account 
/// previous creation of object
//////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_ARM_INF_DoubleDigital_Create(
		LPXLOPER XL_MainIndex,
		LPXLOPER XL_MainType,
		LPXLOPER XL_MainLeverage,
		LPXLOPER XL_MainStrike,
		LPXLOPER XL_SubIndex,
		LPXLOPER XL_SubType,
		LPXLOPER XL_SubLeverage,
		LPXLOPER XL_SubStrike,
		LPXLOPER XL_Spread,
		LPXLOPER XL_Notional )
{
	bool PersistentInXL = true;
	return Local_ARM_INF_DoubleDigital_CreateCommon(
		XL_MainIndex,
		XL_MainType,
		XL_MainLeverage,
		XL_MainStrike,
		XL_SubIndex,
		XL_SubType,
		XL_SubLeverage,
		XL_SubStrike,
		XL_Notional,
		XL_Spread,
		PersistentInXL);
}

///////////////////////////////////
/// version for non persistentInXl
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_PXL_ARM_INF_DoubleDigital_Create(
		LPXLOPER XL_MainIndex,
		LPXLOPER XL_MainType,
		LPXLOPER XL_MainLeverage,
		LPXLOPER XL_MainStrike,
		LPXLOPER XL_SubIndex,
		LPXLOPER XL_SubType,
		LPXLOPER XL_SubLeverage,
		LPXLOPER XL_SubStrike,
		LPXLOPER XL_Spread,
		LPXLOPER XL_Notional)
{
	 bool PersistentInXL = true;
	 return Local_ARM_INF_DoubleDigital_CreateCommon(
		XL_MainIndex,
		XL_MainType,
		XL_MainLeverage,
		XL_MainStrike,
		XL_SubIndex,
		XL_SubType,
		XL_SubLeverage,
		XL_SubStrike,
		XL_Spread,
		XL_Notional,
		PersistentInXL);
}

LPXLOPER Local_ARM_INF_Corridor_CreateCommon(
		LPXLOPER XL_MainIndex,
		LPXLOPER XL_MainType,
		LPXLOPER XL_MainLeverage,
		LPXLOPER XL_SubIndex,
		LPXLOPER XL_SubType,
		LPXLOPER XL_SubLeverage,
		LPXLOPER XL_StrikeInf,
		LPXLOPER XL_StrikeSup,
		LPXLOPER XL_Notional,
		bool PersistentInXL ){

	static XLOPER XL_result;
	ARM_result C_result;
	ARM_XL_TRY_BLOCK_BEGIN	{
	
		ARM_NOCALCIFWIZ();
		static int error;
		static char* reason = "";

		CCString C_MainIndex;
		CCString C_SubIndex;
		CCString C_MainType;
		CCString C_SubType;

		long MainLeverage;
		long SubLeverage;
		long StrikeInf;
		long StrikeSup;
		long Notional;
	
		XL_readStrCellWD(	XL_MainIndex,		C_MainIndex,		"NO",	" ARM_ERR: Main Index:	string expected",			C_result);
		XL_readStrCellWD(	XL_SubIndex,		C_SubIndex,			"NO",	" ARM_ERR: Sub  Index:	string expected",			C_result);
		XL_readStrCellWD(	XL_MainType,		C_MainType,			"IR",	" ARM_ERR: Main Type:	string expected",			C_result);
		XL_readStrCellWD(	XL_SubType,			C_SubType,			"IR",	" ARM_ERR: Sub  Type:	string expected",			C_result);

		string	MainIndex	=	CCSTringToSTLString( C_MainIndex);
		string	SubIndex	=	CCSTringToSTLString( C_SubIndex );
		string	MainType	=	CCSTringToSTLString( C_MainType);
		string	SubType		=	CCSTringToSTLString( C_SubType );

		XL_GETOBJIDWD	(	XL_MainLeverage,	MainLeverage,		"NULL OBJECT",	" ARM_ERR: Main Leverage: Object expected",	C_result);
		XL_GETOBJIDWD	(	XL_SubLeverage,		SubLeverage ,		"NULL OBJECT",	" ARM_ERR: Sub  Leverage: Object expected",	C_result);
		XL_GETOBJIDWD	(	XL_StrikeInf,		StrikeInf ,			"NULL OBJECT",	" ARM_ERR: Strike Inf	: Object expected",	C_result);
		XL_GETOBJIDWD	(	XL_StrikeSup,		StrikeSup ,			"NULL OBJECT",	" ARM_ERR: Strike Sup	: Object expected",	C_result);
		XL_GETOBJIDWD	(	XL_Notional,		Notional ,			"NULL OBJECT",	" ARM_ERR: Notional		: Object expected",	C_result);

		exportFunc9Args<string, string, long, 			 
						 string, string, long,
						 long,	 long,   long> ourFunc(
			MainIndex,
			MainType,
			MainLeverage,
			SubIndex,
			SubType,
			SubLeverage,
			StrikeInf,
			StrikeSup,
			Notional,
			ARMLOCAL_InfCorridor_Create);

		fillXL_Result( LOCAL_INFHYBRIDPAYOFF_CLASS, ourFunc, C_result, XL_result, PersistentInXL );

	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_MktDataManager_MktDataGetSet_Common" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


//////////////////////////////////
/// version that takes into account 
/// previous creation of object
//////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_ARM_INF_Corridor_Create(
		LPXLOPER XL_MainIndex,
		LPXLOPER XL_MainType,
		LPXLOPER XL_MainLeverage,
		LPXLOPER XL_SubIndex,
		LPXLOPER XL_SubType,
		LPXLOPER XL_SubLeverage,
		LPXLOPER XL_StrikeInf,
		LPXLOPER XL_StrikeSup,
		LPXLOPER XL_Notional )
{
	bool PersistentInXL = true;
	return Local_ARM_INF_Corridor_CreateCommon(
		XL_MainIndex,
		XL_MainType,
		XL_MainLeverage,
		XL_SubIndex,
		XL_SubType,
		XL_SubLeverage,
		XL_StrikeInf,
		XL_StrikeSup,
		XL_Notional,
		PersistentInXL);
}

///////////////////////////////////
/// version for non persistentInXl
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_PXL_ARM_INF_Corridor_Create(
		LPXLOPER XL_MainIndex,
		LPXLOPER XL_MainType,
		LPXLOPER XL_MainLeverage,
		LPXLOPER XL_SubIndex,
		LPXLOPER XL_SubType,
		LPXLOPER XL_SubLeverage,
		LPXLOPER XL_StrikeInf,
		LPXLOPER XL_StrikeSup,
		LPXLOPER XL_Notional)
{
	 bool PersistentInXL = true;
	 return Local_ARM_INF_Corridor_CreateCommon(
		XL_MainIndex,
		XL_MainType,
		XL_MainLeverage,
		XL_SubIndex,
		XL_SubType,
		XL_SubLeverage,
		XL_StrikeInf,
		XL_StrikeSup,
		XL_Notional,
		PersistentInXL);
}


LPXLOPER Local_ARM_INF_GetPriceCommon(
		LPXLOPER XL_InfPricer,
		LPXLOPER XL_Key,
		bool PersistentInXL ){

		static XLOPER XL_result;
		ARM_result C_result;
		ARM_XL_TRY_BLOCK_BEGIN	{

		ARM_NOCALCIFWIZ();
		static int error;
		static char* reason = "";
		
		CCString C_InfPricer;
		XL_readStrCellWD(XL_InfPricer, C_InfPricer,"NULL OBJECT",	" ARM_ERR: Inflation pricer is empty",	C_result);

		CCString	CC_Key;
		XL_readStrCellWD (XL_Key, CC_Key,"Price"," ARM_XXX_ERR: string expected",C_result);

		ARM::ARM_GramFctorArg argResult;
		long retCode= ARMLOCAL_Inf_GetPrice (LocalGetNumObjectId (C_InfPricer), CCSTringToSTLString( CC_Key ), argResult, C_result);

		if( retCode == ARM_KO )
		{
			ARM_ERR();
		}
		else
		{
			ARM_GramFunctorToXLOPER( argResult, XL_result, C_result, true);
			FreeCurCellErr ();
		}
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_MktDataManager_MktDataGetSet_Common" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


//////////////////////////////////
/// version that takes into account 
/// previous creation of object
//////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_ARM_INF_GetPrice(
		LPXLOPER XL_InfPricer,
		LPXLOPER XL_Key)
{
	bool PersistentInXL = true;
	return Local_ARM_INF_GetPriceCommon(
		XL_InfPricer,
		XL_Key,
		PersistentInXL);
}

///////////////////////////////////
/// version for non persistentInXl
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_PXL_ARM_INF_GetPrice(
		LPXLOPER XL_InfPricer,
		LPXLOPER XL_Key)
{
	bool PersistentInXL = false;
	return Local_ARM_INF_GetPriceCommon(
		XL_InfPricer,
		XL_Key,
		PersistentInXL);
}

LPXLOPER Local_ARM_INF_GetScheduleCommon(
		LPXLOPER XL_InfLeg,
		LPXLOPER XL_Key,
		bool PersistentInXL ){

		static XLOPER XL_result;
		ARM_result C_result;
		ARM_XL_TRY_BLOCK_BEGIN	{

		ARM_NOCALCIFWIZ();
		static int error;
		static char* reason = "";
		
		CCString C_InfLeg;
		XL_readStrCellWD(XL_InfLeg, C_InfLeg,"NULL OBJECT",	" ARM_ERR: Inflation Leg is not built",	C_result);

		CCString	CC_Key;
		XL_readStrCellWD (XL_Key, CC_Key,"Price"," ARM_XXX_ERR: string expected",C_result);

		ARM::ARM_GramFctorArg argResult;
		long retCode= ARMLOCAL_Inf_GetSchedule (LocalGetNumObjectId (C_InfLeg), CCSTringToSTLString( CC_Key ), argResult, C_result);

		if( retCode == ARM_KO )
		{
			ARM_ERR();
		}
		else
		{
			ARM_GramFunctorToXLOPER( argResult, XL_result, C_result, true);
			FreeCurCellErr ();
		}
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_MktDataManager_MktDataGetSet_Common" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}
//////////////////////////////////
/// version that takes into account 
/// previous creation of object
//////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_ARM_INF_GetSchedule(
		LPXLOPER XL_InfLeg,
		LPXLOPER XL_Key)
{
	bool PersistentInXL = true;
	return Local_ARM_INF_GetScheduleCommon(
		XL_InfLeg,
		XL_Key,
		PersistentInXL);
}

///////////////////////////////////
/// version for non persistentInXl
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_PXL_ARM_INF_GetSchedule(
		LPXLOPER XL_InfLeg,
		LPXLOPER XL_Key)
{
	bool PersistentInXL = false;
	return Local_ARM_INF_GetScheduleCommon(
		XL_InfLeg,
		XL_Key,
		PersistentInXL);
}
/***********************************************

            ARM_INF_EQHWSV_Laplace

	Inputs :
		return the laplace transform of an underlying distributed as EQHWSV

***********************************************/

LPXLOPER Local_ARM_INF_EQHWSV_Laplace(
	LPXLOPER XL_ModelId,
	LPXLOPER XL_evalTime,
	LPXLOPER XL_startTime,
	LPXLOPER XL_endTime,
	LPXLOPER XL_xt,
	LPXLOPER XL_vt,
	LPXLOPER XL_k_real,
	LPXLOPER XL_k_imag,
	LPXLOPER XL_isReal){

	ADD_LOG("Local_ARM_INF_EQHWSV_Laplace");

	static XLOPER XL_result;
	ARM_result C_result;

	ARM_XL_TRY_BLOCK_BEGIN
	{
		ARM_NOCALCIFWIZ();
		static int error;
		static char* reason = "";

		CCString C_ModelId;
		XL_readStrCell(XL_ModelId,C_ModelId,	" ARM_ERR: model id: object expected",C_result);

		double C_evalTime;
		XL_readNumCell(XL_evalTime,C_evalTime,	" ARM_ERR: evalTime : numeric expected",C_result);

		double C_startTime;
		XL_readNumCell(XL_startTime,C_startTime," ARM_ERR: startTime : numeric expected",C_result);

		double C_endTime;
		XL_readNumCell(XL_endTime,C_endTime,	" ARM_ERR: startTime : numeric expected",C_result);

		double D_xt=0.0;
		double C_xt;
		XL_readNumCellWD(XL_xt,C_xt,D_xt,		" ARM_ERR: xt : numeric expected",C_result);

		double D_vt=1.0;
		double C_vt;
		XL_readNumCellWD(XL_vt,C_vt,D_vt,		" ARM_ERR: vt : numeric expected",C_result);

		double D_k_real=0.0;
		double C_k_real;
		XL_readNumCellWD(XL_k_real,C_k_real,D_k_real,	" ARM_ERR: k_real : numeric expected",C_result);

		double D_k_imag=0.0;
		double C_k_imag;
		XL_readNumCellWD(XL_k_imag,C_k_imag,D_k_imag,	" ARM_ERR: k_imag : numeric expected",C_result);


		CCString C_isReal;
		bool isReal = false;
		XL_readStrCell(XL_isReal,C_isReal,	" ARM_ERR: model id: object expected",C_result);
		if (  CCSTringToSTLString( C_isReal ) =="Y"  || CCSTringToSTLString( C_isReal ).size() ==0)
			isReal = true;

		long retCode = ARMLOCAL_InfEqHwSV_Laplace (LocalGetNumObjectId (C_ModelId),C_evalTime, C_startTime, C_endTime, C_xt, C_vt, C_k_real, C_k_imag, isReal, C_result);

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

	}
	ARM_XL_TRY_BLOCK_END
	ARM_XL_CATCH_ARM_EXPT
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_ARM_Price" )
	return (LPXLOPER)&XL_result;
}

LPXLOPER Local_ARM_INF_EQHWSV_Density(
	LPXLOPER XL_ModelId,
	LPXLOPER XL_evalTime,
	LPXLOPER XL_startTime,
	LPXLOPER XL_endTime,
	LPXLOPER XL_xt,
	LPXLOPER XL_vt,
	LPXLOPER XL_x,
	LPXLOPER XL_period,	
	LPXLOPER XL_frequency){

	ADD_LOG("Local_ARM_INF_EQHWSV_Density");

	static XLOPER XL_result;
	ARM_result C_result;

	ARM_XL_TRY_BLOCK_BEGIN
	{
		ARM_NOCALCIFWIZ();
		static int error;
		static char* reason = "";

		CCString C_ModelId;
		XL_readStrCell(XL_ModelId,C_ModelId,	" ARM_ERR: model id: object expected",C_result);

		double C_evalTime;
		XL_readNumCell(XL_evalTime,C_evalTime,	" ARM_ERR: evalTime : numeric expected",C_result);

		double C_startTime;
		XL_readNumCell(XL_startTime,C_startTime," ARM_ERR: startTime : numeric expected",C_result);

		double C_endTime;
		XL_readNumCell(XL_endTime,C_endTime,	" ARM_ERR: startTime : numeric expected",C_result);

		double D_xt=0.0;
		double C_xt;
		XL_readNumCellWD(XL_xt,C_xt,D_xt,		" ARM_ERR: xt : numeric expected",C_result);

		double D_vt=1.0;
		double C_vt;
		XL_readNumCellWD(XL_vt,C_vt,D_vt,		" ARM_ERR: vt : numeric expected",C_result);

		double D_x=1.0;
		double C_x;
		XL_readNumCellWD(XL_x,C_x,D_x,			" ARM_ERR: x : numeric expected",C_result);

		double D_period=1.0;
		double C_period;
		XL_readNumCellWD(XL_period,C_period,D_period,			" ARM_ERR: period : numeric expected",C_result);
		
		double D_frequency=1.0;
		double C_frequency;
		XL_readNumCellWD(XL_frequency,C_frequency,D_frequency,	" ARM_ERR: frequency : numeric expected",C_result);


		long retCode = ARMLOCAL_InfEqHwSV_Density (LocalGetNumObjectId (C_ModelId),C_evalTime, C_startTime, C_endTime, C_xt, C_vt, C_x, C_period, C_frequency, C_result);

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

	}
	ARM_XL_TRY_BLOCK_END
	ARM_XL_CATCH_ARM_EXPT
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_ARM_Price" )
	return (LPXLOPER)&XL_result;
}

LPXLOPER Local_ARM_INF_GetAdjCorrelCommon(
		LPXLOPER XL_InfBsSmiledModel,
		bool PersistentInXL ){

		static XLOPER XL_result;
		ARM_result C_result;
		ARM_XL_TRY_BLOCK_BEGIN	{

		ARM_NOCALCIFWIZ();
		static int error;
		static char* reason = "";
		
		CCString C_InfBsSmiledModel;
		XL_readStrCellWD(XL_InfBsSmiledModel,C_InfBsSmiledModel,"NULL OBJECT",	" ARM_ERR: Inflation model is empty",	C_result);

		exportFunc1Arg<	long > ourFunc(
			LocalGetNumObjectId( C_InfBsSmiledModel ),
			ARMLOCAL_Inf_GetAdjCorrel);

		fillXL_Result( LOCAL_INFHYBRIDPAYOFF_CLASS, ourFunc, C_result, XL_result, PersistentInXL );
  
		}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_MktDataManager_MktDataGetSet_Common" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}

//////////////////////////////////
/// version that takes into account 
/// previous creation of object
//////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_ARM_INF_GetAdjCorrel(
		LPXLOPER XL_InfBsSmiledModel)
{
	bool PersistentInXL = true;
	return Local_ARM_INF_GetAdjCorrelCommon(
		XL_InfBsSmiledModel,
		PersistentInXL);
}

///////////////////////////////////
/// version for non persistentInXl
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_PXL_ARM_INF_GetAdjCorrel(
		LPXLOPER XL_InfBsSmiledModel)
{
	bool PersistentInXL = false;
	return Local_ARM_INF_GetAdjCorrelCommon(
		XL_InfBsSmiledModel,
		PersistentInXL);
}
