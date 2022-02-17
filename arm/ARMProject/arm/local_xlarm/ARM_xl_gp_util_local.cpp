/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file ARM_xl_gp_util_local.cpp
 *
 *  \brief file for the various utilities in the generic pricer
 *
 *	\author  J M Prie
 *	\version 1.0
 *	\date January 2004
 */

#include <ARM\libarm_local\firstToBeIncluded.h>
#include <libCCxll\CCxll.h>

#include <ARM\libarm_local\ARM_local_gp_util.h>
#include "ARM_xl_gp_util_local.h"
#include "ARM_xl_wrapper_local.h"
#include "ARM_xl_trycatch_local.h"
#include "ARM_xl_gp_fctorhelper.h"

#include "ARM_xl_gp_curve_local.h"
#include <ARM\libarm_local\ARM_local_gp_curve.h>
#include <libCCxll\CCxll.h>
#include "ARM_xl_wrapper_local.h"
#include "ARM_xl_trycatch_local.h"
#include "ARM_xl_gp_fctorhelper.h"
#include "ARM_local_interface.h"
#include <GP_Base\gpbase\vectormanip.h>
#include <GP_Base\gpbase\stringmanip.h>


#include <gpbase\gpmatrix.h>

#include "util\tech_macro.h"

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

///----------------------------------------------
///----------------------------------------------
///    Trigo Matrix Function
/// Inputs :
///     alpha : constant to compute coefficient
///     n : size of the square matrix
///----------------------------------------------
///----------------------------------------------
__declspec(dllexport) LPXLOPER WINAPI Local_TrigoMatrix(
    LPXLOPER XL_N,
	LPXLOPER XL_Alpha)
{
	ADD_LOG("Local_TrigoMatrix");
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

		double C_N;
		XL_readNumCell(XL_N,C_N," ARM_ERR: n: double expected",C_result);
		double C_Alpha;
		XL_readNumCell(XL_Alpha,C_Alpha," ARM_ERR: alpha: double expected",C_result);

		VECTOR<double> vectorResult;
		long nbRows;
		long nbCols;

		/// call the function
		long retCode=ARMLOCAL_TrigoMatrix(
				C_N,
				C_Alpha,
				vectorResult,
				nbRows,
				nbCols,
				C_result);

		/// return the result as an LPXLOPER
		if (retCode == ARM_OK)
		{
			FreeCurCellErr ();
			XL_writeNumMatrixSize( XL_result, vectorResult, nbRows, nbCols, " ARM_ERR: Could not set the num matrix", C_result );
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_TrigoMatrix" )

	return (LPXLOPER)&XL_result;
}


///----------------------------------------------
///----------------------------------------------
///    Local Covariance & Correlation Functions
/// Inputs :
///     Model Id
///     fromDate, toDate
///     startDate1, endDate1
///     startDate2, endDate2
///----------------------------------------------
///----------------------------------------------
LPXLOPER Local_LocalCovarianceCommon(
	LPXLOPER XL_ModelId,
	LPXLOPER XL_UnderlyingType,
    LPXLOPER XL_FromDate,
    LPXLOPER XL_ToDate,
    LPXLOPER XL_StartDate1,
    LPXLOPER XL_EndDate1,
    LPXLOPER XL_StartDate2,
    LPXLOPER XL_EndDate2,
    LPXLOPER XL_StartDate3,
    LPXLOPER XL_EndDate3,
	LPXLOPER XL_StartDate4,
    LPXLOPER XL_EndDate4,
    long (*Function)(long,string,double,double,double,string,double,long,double,string,double,long,double,string,double,long,double,string,double,long,ARM_result&))
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
		XL_readStrCell( XL_ModelId, C_modelStrId,	" ARM_ERR: Model Id: Object expected",C_result);
		long C_ModelId = LocalGetNumObjectId(C_modelStrId);

		CCString underlyingTypeStr;
		XL_readStrCell(XL_UnderlyingType,underlyingTypeStr," ARM_ERR: Funding Index Term: string expected",C_result);
		string C_UnderlyingType=CCSTringToSTLString(underlyingTypeStr);		

		double C_FromDate;
		XL_readNumCell(XL_FromDate,C_FromDate," ARM_ERR: from date: date expected",C_result);
		double C_ToDate;
		XL_readNumCell(XL_ToDate,C_ToDate," ARM_ERR: to date: date expected",C_result);

		double C_StartDate1;
		XL_readNumCell(XL_StartDate1,C_StartDate1," ARM_ERR: first start date: date expected",C_result);

		double C_EndDate1;
		CCString C_EndDate1Str;
		long     EndDate1Type;
		XL_readStrOrNumCell(XL_EndDate1, C_EndDate1Str, C_EndDate1, EndDate1Type,
			   " ARM_ERR: first tenor or first end date expected",C_result);	
		
		string EndDate1Tenor("") ;
		if(EndDate1Type == XL_TYPE_STRING)
		{
			EndDate1Tenor = CCSTringToSTLString(C_EndDate1Str);
		}
		
		double C_StartDate2;
		XL_readNumCell(XL_StartDate2,C_StartDate2," ARM_ERR: second start date: date expected",C_result);

		double C_EndDate2;
		CCString C_EndDate2Str;
		long     EndDate2Type;
		XL_readStrOrNumCell(XL_EndDate2, C_EndDate2Str, C_EndDate2, EndDate2Type,
			   " ARM_ERR: second tenor or second end date expected",C_result);	

		string EndDate2Tenor("") ;
		if(EndDate2Type == XL_TYPE_STRING)
		{
			EndDate2Tenor = CCSTringToSTLString(C_EndDate2Str);
		}

		/// Default is StartDate1
		double C_StartDate3;
		XL_readNumCellWD(XL_StartDate3,C_StartDate3,C_StartDate1," ARM_ERR: third start date: date expected",C_result);

		/// Default is EndDate1
		double C_EndDate3;
		CCString C_EndDate3Str;
		long     EndDate3Type;
		XL_readStrOrNumCellWD(XL_EndDate3, C_EndDate3Str, C_EndDate3, C_EndDate1, EndDate3Type,
			   " ARM_ERR: third tenor or third end date expected",C_result);	

		string EndDate3Tenor("") ;
		if(EndDate3Type == XL_TYPE_STRING)
		{
			EndDate3Tenor = CCSTringToSTLString(C_EndDate3Str);
		}
		else if(C_EndDate3 == C_EndDate1)
		{
			/// Default
			if(EndDate1Type == XL_TYPE_STRING)
			{
				EndDate3Type	= EndDate1Type;
				EndDate3Tenor	= EndDate1Tenor;
			}
		}

		/// Default is StartDate1
		double C_StartDate4;
		XL_readNumCellWD(XL_StartDate4,C_StartDate4,C_StartDate1," ARM_ERR: third start date: date expected",C_result);

		/// Default is EndDate1
		double C_EndDate4;
		CCString C_EndDate4Str;
		long     EndDate4Type;
		XL_readStrOrNumCellWD(XL_EndDate4, C_EndDate4Str, C_EndDate4, C_EndDate1, EndDate4Type,
			   " ARM_ERR: fourth tenor or fourth end date expected",C_result);	

		string EndDate4Tenor("") ;
		if(EndDate4Type == XL_TYPE_STRING)
		{
			EndDate4Tenor = CCSTringToSTLString(C_EndDate4Str);
		}
		else if(C_EndDate4 == C_EndDate1)
		{
			/// Default
			if(EndDate1Type == XL_TYPE_STRING)
			{
				EndDate4Type	= EndDate1Type;
				EndDate4Tenor	= EndDate1Tenor;
			}
		}


		/// call the function
		long retCode=(*Function)(
				C_ModelId,
				C_UnderlyingType,
				C_FromDate,
				C_ToDate,
				C_StartDate1,
				EndDate1Tenor, 
				C_EndDate1,
				EndDate1Type,
				C_StartDate2,
				EndDate2Tenor, 
				C_EndDate2, 
				EndDate2Type,
				C_StartDate3,
				EndDate3Tenor, 
				C_EndDate3, 
				EndDate3Type,
				C_StartDate4,
				EndDate4Tenor, 
				C_EndDate4, 
				EndDate4Type,
				C_result);

		/// return the result as an LPXLOPER
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
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_LocalCovarianceCommon" )

	return (LPXLOPER)&XL_result;
}						 



/*__declspec(dllexport) LPXLOPER WINAPI Local_LocalCovariance(
	LPXLOPER XL_ModelId,
	LPXLOPER XL_UnderlyingType,
    LPXLOPER XL_FromDate,
    LPXLOPER XL_ToDate,
    LPXLOPER XL_StartDate1,
    LPXLOPER XL_EndDate1,
    LPXLOPER XL_StartDate2,
    LPXLOPER XL_EndDate2,
    LPXLOPER XL_StartDate3,
    LPXLOPER XL_EndDate3)
{
	ADD_LOG("Local_LocalCovariance");
    return Local_LocalCovarianceCommon(XL_ModelId,
		XL_UnderlyingType,
        XL_FromDate,XL_ToDate,
        XL_StartDate1,XL_EndDate1,XL_StartDate2,XL_EndDate2,XL_StartDate3,XL_EndDate3,
        ARMLOCAL_LocalCovariance);
}

__declspec(dllexport) LPXLOPER WINAPI Local_LocalCorrelation(
	LPXLOPER XL_ModelId,
	LPXLOPER XL_UnderlyingType,
    LPXLOPER XL_FromDate,
    LPXLOPER XL_ToDate,
    LPXLOPER XL_StartDate1,
    LPXLOPER XL_EndDate1,
    LPXLOPER XL_StartDate2,
    LPXLOPER XL_EndDate2,
    LPXLOPER XL_StartDate3,
    LPXLOPER XL_EndDate3)
{
	ADD_LOG("Local_LocalCorrelation");
    return Local_LocalCovarianceCommon(XL_ModelId,
		XL_UnderlyingType,
        XL_FromDate,XL_ToDate,
        XL_StartDate1,XL_EndDate1,XL_StartDate2,XL_EndDate2,XL_StartDate3,XL_EndDate3,
        ARMLOCAL_LocalCorrelation);
}*/

__declspec(dllexport) LPXLOPER WINAPI Local_LocalCovariance(
	LPXLOPER XL_ModelId,
	LPXLOPER XL_UnderlyingType,
    LPXLOPER XL_FromDate,
    LPXLOPER XL_ToDate,
    LPXLOPER XL_StartDate1,
    LPXLOPER XL_EndDate1,
    LPXLOPER XL_StartDate2,
    LPXLOPER XL_EndDate2,
    LPXLOPER XL_StartDate3,
    LPXLOPER XL_EndDate3,
    LPXLOPER XL_StartDate4,
    LPXLOPER XL_EndDate4)
{
	ADD_LOG("Local_LocalCovariance");
    return Local_LocalCovarianceCommon(XL_ModelId,
		XL_UnderlyingType,
        XL_FromDate,XL_ToDate,
        XL_StartDate1,XL_EndDate1,XL_StartDate2,XL_EndDate2,XL_StartDate3,XL_EndDate3,XL_StartDate4,XL_EndDate4,
        ARMLOCAL_LocalCovariance);
}

__declspec(dllexport) LPXLOPER WINAPI Local_LocalCorrelation(
	LPXLOPER XL_ModelId,
	LPXLOPER XL_UnderlyingType,
    LPXLOPER XL_FromDate,
    LPXLOPER XL_ToDate,
    LPXLOPER XL_StartDate1,
    LPXLOPER XL_EndDate1,
    LPXLOPER XL_StartDate2,
    LPXLOPER XL_EndDate2,
    LPXLOPER XL_StartDate3,
    LPXLOPER XL_EndDate3,
    LPXLOPER XL_StartDate4,
    LPXLOPER XL_EndDate4)
{
	ADD_LOG("Local_LocalCorrelation");
    return Local_LocalCovarianceCommon(XL_ModelId,
		XL_UnderlyingType,
        XL_FromDate,XL_ToDate,
        XL_StartDate1,XL_EndDate1,XL_StartDate2,XL_EndDate2,XL_StartDate3,XL_EndDate3,XL_StartDate4,XL_EndDate4,
        ARMLOCAL_LocalCorrelation);
}


__declspec(dllexport) LPXLOPER WINAPI Local_IntegratedCorrelation_Common(
	LPXLOPER XL_ModelId,
	LPXLOPER XL_Tenor1,
	LPXLOPER XL_Tenors,
	LPXLOPER XL_Expiries,
	bool PersistentInXL )
{
	ADD_LOG("Local_IntegratedCorrelation_Common");
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
		XL_readStrCell( XL_ModelId, C_modelStrId, " ARM_ERR: Model Id: Object expected", C_result);
		long C_ModelId = LocalGetNumObjectId(C_modelStrId);

		double C_Tenor1;
		XL_readNumCell( XL_Tenor1, C_Tenor1, " ARM_ERR: Tenor1: ", C_result);
		
		VECTOR<double> C_Tenors;
		XL_readNumVector( XL_Tenors, C_Tenors,	" ARM_ERR: Tenors vector: array of numeric expected", C_result);

		VECTOR<double> C_Expiries;
		XL_readNumVector( XL_Expiries, C_Expiries,	" ARM_ERR: Expiries vector: array of numeric expected", C_result);

		/// a function with a context
		exportFunc4Args<long, double, vector< double>, vector< double> >  ourFunc(C_ModelId, C_Tenor1, C_Tenors, C_Expiries, ARMLOCAL_IntegratedCorrelation_Compute );

		/// call the general function
		fillXL_Result_withName( ourFunc, C_result, XL_result, PersistentInXL );

	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_IntegratedCorrelation_Common" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}



__declspec(dllexport) LPXLOPER WINAPI Local_IntegratedCorrelation(
	LPXLOPER XL_ModelId,
	LPXLOPER XL_Tenor1,
	LPXLOPER XL_Tenors,
	LPXLOPER XL_Expiries )
{
	ADD_LOG("Local_IntegratedCorrelation");
	return Local_IntegratedCorrelation_Common(
		XL_ModelId,
		XL_Tenor1,
		XL_Tenors,
		XL_Expiries,
		true );
}

__declspec(dllexport) LPXLOPER WINAPI Local_PXL_IntegratedCorrelation(
	LPXLOPER XL_ModelId,
	LPXLOPER XL_Tenor1,
	LPXLOPER XL_Tenors,
	LPXLOPER XL_Expiries )
{
	ADD_LOG("Local_PXL_IntegratedCorrelation");
	return Local_IntegratedCorrelation_Common(
		XL_ModelId,
		XL_Tenor1,
		XL_Tenors,
		XL_Expiries,
		false );
}



///----------------------------------------------
///----------------------------------------------
///    To convert Fx Mkt Data to Totem Format
/// Inputs :
///     Portfolios Id
///     fx ATM Vol
///     Delta calls
/// 	Delta puts
///     bs Model
///----------------------------------------------
///----------------------------------------------
__declspec(dllexport) LPXLOPER WINAPI Local_FxMktToTotemCalibrate_Common(
	LPXLOPER XL_portfolioId,
	LPXLOPER XL_ATMVol,
	LPXLOPER XL_DeltaCalls,
	LPXLOPER XL_DeltaPuts,
	LPXLOPER XL_Model,
	LPXLOPER XL_InitPoint,
	LPXLOPER XL_LowerBound,
	LPXLOPER XL_UpperBound,
	LPXLOPER XL_MaxIter,
	LPXLOPER XL_AlgoType,
	LPXLOPER XL_Tolerance,
	LPXLOPER XL_OutPutId,
	bool PersistentInXL )
{
	ADD_LOG("Local_FxMktToTotemCalibrate_Common");
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

		/// portfolio Id 
		CCString C_portfolioId;
		XL_readStrCell( XL_portfolioId, C_portfolioId, " ARM_ERR: Portfolio Id: Object expected",C_result);
		long pfId = LocalGetNumObjectId(C_portfolioId);

		double C_atmVol;
		XL_readNumCell( XL_ATMVol, C_atmVol, " ARM_ERR: ATM Vol: ", C_result);
		
		VECTOR<double> C_deltaCalls;
		XL_readNumVector( XL_DeltaCalls, C_deltaCalls,	" ARM_ERR: Delta Calls vector: array of numeric expected", C_result);

		VECTOR<double> C_deltaPuts;
		XL_readNumVector( XL_DeltaPuts, C_deltaPuts,	" ARM_ERR: Delta Puts vector: array of numeric expected", C_result);

		/// Model Id g
		CCString C_modelStr;
		XL_readStrCell( XL_Model, C_modelStr, " ARM_ERR: Model Id: Object expected",C_result);
		long modelId = LocalGetNumObjectId(C_modelStr);

		/// Initial point
		CCString InitPointStr;
		XL_readStrCellWD(XL_InitPoint, InitPointStr,"NULL"," ARM_ERR: Initial variable: object expected",C_result);
		long initPointId = InitPointStr == "NULL"? ARM_NULL_OBJECT: LocalGetNumObjectId (InitPointStr);

		/// lower bound
		CCString LowerBoundStr;
		XL_readStrCellWD(XL_LowerBound, LowerBoundStr,"NULL"," ARM_ERR: Lower Bound: object expected",C_result);
		long lowerBoundId = LowerBoundStr == "NULL"? ARM_NULL_OBJECT: LocalGetNumObjectId (LowerBoundStr);

		/// upper bound
		CCString UpperBoundStr;
		XL_readStrCellWD(XL_UpperBound, UpperBoundStr,"NULL"," ARM_ERR: Upper Bound: object expected",C_result);
		long upperBoundId = UpperBoundStr == "NULL"? ARM_NULL_OBJECT: LocalGetNumObjectId (UpperBoundStr);

		/// MAx Iteraion
		double C_MaxIter;
		double maxTer_default = 50.0;
		XL_readNumCellWD( XL_MaxIter, C_MaxIter, maxTer_default, " ARM_ERR: Max Iteration : numeric expected",	C_result);

		/// Algorithem type
		CCString C_AlgoTypeStr;
		XL_readStrCellWD(XL_AlgoType, C_AlgoTypeStr ,"NLIN_LSQ"," ARM_ERR: Upper Bound: object expected",C_result);
		string AlgoTypeStr = CCSTringToSTLString(C_AlgoTypeStr);

		/// Tolerance
		double C_Tolerance;
		double tolerance_default = 1.0e-6;
		XL_readNumCellWD( XL_Tolerance, C_Tolerance, tolerance_default, " ARM_ERR: Tolerance : numeric expected",	C_result);

		/// Algorithem type
		CCString C_OutPutStr;
		XL_readStrCellWD(XL_OutPutId, C_OutPutStr ,"GENCURVE"," ARM_ERR: outPut Id: object expected",C_result);
		string OutPutStr = CCSTringToSTLString(C_OutPutStr);


		/// a function with a context
		exportFunc12Args< long, double, vector< double>, vector< double>, long, long, long, long, double, string, double, string > 
			ourFunc(pfId, C_atmVol, C_deltaCalls, C_deltaPuts,modelId, initPointId,lowerBoundId,upperBoundId,C_MaxIter,AlgoTypeStr,C_Tolerance, OutPutStr,ARMLOCAL_FxMktToTotemCalibrate );

		/// call the general function
		fillXL_Result_withName( ourFunc, C_result, XL_result, PersistentInXL );

	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_FxMktToTotemCalibrate_Common" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}



__declspec(dllexport) LPXLOPER WINAPI Local_FxMktToTotemCalibrate(
	LPXLOPER XL_portfolioId,
	LPXLOPER XL_ATMVol,
	LPXLOPER XL_DeltaCalls,
	LPXLOPER XL_DeltaPuts,
	LPXLOPER XL_Model,
	LPXLOPER XL_InitPoint,
	LPXLOPER XL_LowerBound,
	LPXLOPER XL_UpperBound,
	LPXLOPER XL_MaxIter,
	LPXLOPER XL_AlgoType,
	LPXLOPER XL_Tolerance,
	LPXLOPER XL_OutPutId)
{
	ADD_LOG("Local_FxMktToTotemCalibrate");
	return Local_FxMktToTotemCalibrate_Common(
	XL_portfolioId,
	XL_ATMVol,
	XL_DeltaCalls,
	XL_DeltaPuts,
	XL_Model,
	XL_InitPoint,
	XL_LowerBound,
	XL_UpperBound,
	XL_MaxIter,
	XL_AlgoType,
	XL_Tolerance,
	XL_OutPutId,
	true );
}

__declspec(dllexport) LPXLOPER WINAPI Local_PXL_FxMktToTotemCalibrate(
	LPXLOPER XL_portfolioId,
	LPXLOPER XL_ATMVol,
	LPXLOPER XL_DeltaCalls,
	LPXLOPER XL_DeltaPuts,
	LPXLOPER XL_Model,
	LPXLOPER XL_InitPoint,
	LPXLOPER XL_LowerBound,
	LPXLOPER XL_UpperBound,
	LPXLOPER XL_MaxIter,
	LPXLOPER XL_AlgoType,
	LPXLOPER XL_Tolerance,
	LPXLOPER XL_OutPutId)
{
	ADD_LOG("Local_PXL_FxMktToTotemCalibrate");
	return Local_FxMktToTotemCalibrate_Common(
	XL_portfolioId,
	XL_ATMVol,
	XL_DeltaCalls,
	XL_DeltaPuts,
	XL_Model,
	XL_InitPoint,
	XL_LowerBound,
	XL_UpperBound,
	XL_MaxIter,
	XL_AlgoType,
	XL_Tolerance,
	XL_OutPutId,
	false );
}




///----------------------------------------------
///----------------------------------------------
///             event Viewer
/// Inputs : None
///----------------------------------------------
///----------------------------------------------

////////////////////////////////////////////////
/// very rapid definition of ourlonglongFunctor
////////////////////////////////////////////////
struct EventViewerFunc : public ARMResultLong2LongFunc
{
	EventViewerFunc(){};
	
	long operator()( ARM_result& result, long objId ){
		return ARMLOCAL_GetEventViewerRep( result, objId);			
	}
};


/////////////////////////////////////////////////////////////
/// central function that does the creation of the XL function
/////////////////////////////////////////////////////////////
LPXLOPER Local_Event_Create_Common( bool PersistentInXL )
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

 		/// use the concept of Functor to transfer the knowledge of
		/// a function with a context
		EventViewerFunc ourFunc;

		/// call the general function
		fillXL_Result( LOCAL_EVENT_VIEWER_CLASS, ourFunc, C_result, XL_result, PersistentInXL );
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_DateStripCombiner_Common" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}
					

////////////////////////////////////////////
//// function to create a representant of the event viewer
////////////////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_EventViewer_Create()
{
	ADD_LOG("Local_EventViewer_Create");
	bool PersistentInXL = true;
	return Local_Event_Create_Common( PersistentInXL );
}

////////////////////////////////////////////
//// function to create a representant of the event viewer
////////////////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_PXL_EventViewer_Create()
{
	ADD_LOG("Local_PXL_EventViewer_Create");
	bool PersistentInXL = false;
	return Local_Event_Create_Common( PersistentInXL );
}

////////////////////////////////////////////
//// function to reset message in the event viewer
////////////////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_EventViewer_ResetMssg(
	LPXLOPER XL_EventViewerId )
{
	ADD_LOG("Local_EventViewer_ResetMssg");
	
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

		long C_eventViewerId;
		XL_GETOBJID( XL_EventViewerId, C_eventViewerId,	" ARM_ERR: event viewer id: object expected",	C_result);

		long retCode = ARMLOCAL_GetEventViewer_ResetMssg(
				C_eventViewerId,
				C_result );

		/// feed the LPXLOPER object result 
		if (retCode == ARM_OK)
		{
			FreeCurCellErr ();
			XL_result.xltype  = xltypeStr;
			XL_result.val.str = XL_StrC2StrPascal (C_result.getString() );
			XL_result.xltype |= xlbitDLLFree;
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_DateStripCombiner_Common" )

	return (LPXLOPER)&XL_result;
}

////////////////////////////////////////////
//// function to set verbose mode in the event viewer
////////////////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_EventViewer_SetVerboseMode(
	LPXLOPER XL_verboseMode )
{
	ADD_LOG("Local_EventViewer_SetVerboseMode");
	
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

		double C_verboseMode;
		XL_readNumCell( XL_verboseMode, C_verboseMode,	" ARM_ERR: verbose mode Flag: boolean compatible expected",	C_result);
		bool C_verboseModeFlag = C_verboseMode != 0;

		long retCode = ARMLOCAL_EventViewer_SetVerboseMode(
				C_verboseModeFlag,
				C_result );

		/// feed the LPXLOPER object result 
		if (retCode == ARM_OK)
		{
			FreeCurCellErr ();
			XL_result.xltype  = xltypeStr;
			XL_result.val.str = XL_StrC2StrPascal (C_result.getString() );
			XL_result.xltype |= xlbitDLLFree;
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_DateStripCombiner_Common" )

	return (LPXLOPER)&XL_result;
}



///----------------------------------------------
///----------------------------------------------
///             error Viewer
/// Inputs : None
///----------------------------------------------
///----------------------------------------------

////////////////////////////////////////////////
/// very rapid definition of ourlonglongFunctor
////////////////////////////////////////////////
struct ErrViewerFunc : public ARMResultLong2LongFunc
{
	ErrViewerFunc (bool resetFlag): itsResetFlag(resetFlag){};
	
	long operator()( ARM_result& result, long objId ){
		return ARMLOCAL_ErrViewer_Create( itsResetFlag, result, objId);			
	}
private:
	bool itsResetFlag;
};


/////////////////////////////////////////////////////////////
/// central function that does the creation of the XL function
/////////////////////////////////////////////////////////////
LPXLOPER Local_ErrViewer_Create_Common( 
	LPXLOPER XL_resetFlag,
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

		double C_resetFlagDble;
		double C_resetFlagDefault=0;
		XL_readNumCellWD( XL_resetFlag, C_resetFlagDble,C_resetFlagDefault, " ARM_ERR: reset Flag: boolean compatible expected",	C_result);
		bool C_resetFlag = C_resetFlagDble != 0;

 		/// use the concept of Functor to transfer the knowledge of
		/// a function with a context
		ErrViewerFunc ourFunc(C_resetFlag);

		/// call the general function
		fillXL_Result( LOCAL_ERR_VIEWER_CLASS, ourFunc, C_result, XL_result, PersistentInXL );
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_DateStripCombiner_Common" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


////////////////////////////////////////////
//// function to create a representant of the event viewer
////////////////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_ErrViewer_Create(
	LPXLOPER XL_resetFlag )
{
	ADD_LOG("Local_ErrViewer_Create");
	bool PersistentInXL = true;
	return Local_ErrViewer_Create_Common( XL_resetFlag, PersistentInXL );
}


////////////////////////////////////////////
//// function to create a representant of the event viewer
////////////////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_PXL_ErrViewer_Create(
	LPXLOPER XL_resetFlag )
{
	ADD_LOG("Local_PXL_ErrViewer_Create");
	bool PersistentInXL = false;
	return Local_ErrViewer_Create_Common( XL_resetFlag, PersistentInXL );
}

///----------------------------------------------
///----------------------------------------------
///             GP Vector creation
/// Inputs : Values
///----------------------------------------------
///----------------------------------------------

/////////////////////////////////////////////////////////////
/// central function that does the creation of the XL function
/////////////////////////////////////////////////////////////
LPXLOPER Local_GPVector_Create_Common( 
	LPXLOPER XL_values,
	LPXLOPER XL_type,
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

		CCString C_Type;
		XL_readStrCellWD( XL_type,	C_Type,	"Real", " ARM_ERR: interpolator type string expected",	C_result);
		string typeName = CCSTringToSTLString( C_Type );
		ARM::stringToUpper( typeName );

		if(typeName == string("REAL")){
			VECTOR<double> C_Values;
			XL_readNumVector( XL_values, C_Values, " ARM_ERR: values vector expected", C_result);
			/// a function with a context
			exportFunc1Arg< vector<double> >  ourFunc(C_Values, ARMLOCAL_GP_Vector_Create );
			
			/// call the general function
			fillXL_Result_withName( ourFunc, C_result, XL_result, PersistentInXL );
		}
		else if(typeName == string("STRING")){
			VECTOR<CCString> C_Values;
			XL_readStrVector( XL_values, C_Values, " ARM_ERR: string vector expected", DOUBLE_TYPE,C_result);
			/// a function with a context
			size_t size = C_Values.size();
			VECTOR<string> C_ValueStr(size);
			for(size_t i(0); i<size; ++i)
				C_ValueStr[i] = CCSTringToSTLString( C_Values[i] );

			exportFunc1Arg< VECTOR<string> >  ourFunc(C_ValueStr, ARMLOCAL_GP_StrVector_Create );
			
			/// call the general function
			fillXL_Result_withName( ourFunc, C_result, XL_result, PersistentInXL );
		}
		else
			throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"Tyep : invalid type, please try Real or String");

	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_GPVector_Create_Common" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


////////////////////////////////////////////
//// function to create a gp vector
////////////////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_GPVector_Create(
	LPXLOPER XL_values ,
	LPXLOPER Xl_type)
{
	ADD_LOG("Local_GPVector_Create");
	bool PersistentInXL = true;
	return Local_GPVector_Create_Common( XL_values, Xl_type, PersistentInXL );
}


////////////////////////////////////////////
//// function to create a gp vector
////////////////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_PXL_GPVector_Create(
	LPXLOPER XL_values ,
	LPXLOPER Xl_type)
{
	ADD_LOG("Local_PXL_GPVector_Create");
	bool PersistentInXL = false;
	return Local_GPVector_Create_Common( XL_values, Xl_type, PersistentInXL );
}






///----------------------------------------------
///----------------------------------------------
///             GP Matrix creation
/// Inputs : Values
///----------------------------------------------
///----------------------------------------------

/////////////////////////////////////////////////////////////
/// central function that does the creation of the XL function
/////////////////////////////////////////////////////////////
LPXLOPER Local_GPMatrix_Create_Common( 
	LPXLOPER XL_values,
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

		VECTOR<double> C_Values;
		long C_rows, C_Cols;
		XL_readNumVectorAndSize( XL_values, C_rows, C_Cols, C_Values, " ARM_ERR: values vector or matrix expected", C_result);
		
		/// a function with a context
		exportFunc3Args< vector<double>, long, long >  ourFunc(C_Values, C_rows, C_Cols, ARMLOCAL_GP_Matrix_Create );

		/// call the general function
		fillXL_Result( LOCAL_GP_MATRIX_CLASS, ourFunc, C_result, XL_result, PersistentInXL );
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_GPMatrix_Create_Common" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


////////////////////////////////////////////
//// function to create a gp matrix
////////////////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_GPMatrix_Create(
	LPXLOPER XL_values )
{
	ADD_LOG("Local_GPMatrix_Create");
	bool PersistentInXL = true;
	return Local_GPMatrix_Create_Common( XL_values, PersistentInXL );
}


////////////////////////////////////////////
//// function to create a gp matrix
////////////////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_PXL_GPMatrix_Create(
	LPXLOPER XL_values )
{
	ADD_LOG("Local_PXL_GPMatrix_Create");
	bool PersistentInXL = false;
	return Local_GPMatrix_Create_Common( XL_values, PersistentInXL );
}

///----------------------------------------------
///----------------------------------------------
///             Regression Calculation
/// Inputs : X matrix, Y values
///----------------------------------------------
///----------------------------------------------

__declspec(dllexport) LPXLOPER WINAPI Local_GP_LeastSquareRegression( 
	LPXLOPER XL_X,
	LPXLOPER XL_Y)
{
	ADD_LOG("Local_GP_LeastSquareRegression");
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

		VECTOR<double> C_X;
		long C_rows, C_cols;
		XL_readNumVectorAndSize( XL_X, C_rows, C_cols, C_X, " ARM_ERR: X: matrix expected", C_result);

		VECTOR<double> C_Y;
		XL_readNumVector(XL_Y, C_Y,	" ARM_ERR: Y: array of numeric expected",	C_result);


		VECTOR<double> C_Coeffs;


		long retCode = ARMLOCAL_GP_LeastSquareRegression(
			C_X,
			C_rows,
			C_cols,
			C_Y,
			C_Coeffs,
			C_result);

		if (retCode == ARM_OK)
		{
			XL_writeNumVector( XL_result, C_Coeffs, "ARM ERR: could not write the output", C_result );
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_GP_LeastSquareRegression" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}

///----------------------------------------------
///----------------------------------------------
///             ACP Calculation
/// Inputs : matrix
///----------------------------------------------
///----------------------------------------------

__declspec(dllexport) LPXLOPER WINAPI Local_GP_ACP( 
	LPXLOPER XL_Matrix)
{
	ADD_LOG("Local_GP_ACP");
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

		VECTOR<double> C_Matrix;
		long C_rows, C_cols;
		XL_readNumVectorAndSize( XL_Matrix, C_rows, C_cols, C_Matrix, " ARM_ERR: Matrix: matrix expected", C_result);


		VECTOR<double> C_OutputMatrix;
		long C_outputRows, C_outputCols;

		long retCode = ARMLOCAL_GP_ACP(
			C_Matrix,
			C_rows,
			C_cols,
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_GP_ACP" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}

///----------------------------------------------
///----------------------------------------------
///             ACP Calculation
/// Inputs : 
///			_ vector
///			_ matrix
///			_ sring
///			_ double
///----------------------------------------------
///----------------------------------------------

__declspec(dllexport) LPXLOPER WINAPI Local_GP_Regression( 
	LPXLOPER XL_Y,
	LPXLOPER XL_X,
	LPXLOPER XL_XInter,
	LPXLOPER XL_RegMode,
	LPXLOPER XL_Span)
{
	ADD_LOG("Local_GP_Regression");
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

		VECTOR<double> C_Y;
		XL_readNumVector( XL_Y, C_Y, " ARM_ERR: vector of numeric expected", C_result);

		VECTOR<double> C_X;
		long C_XRows, C_XCols;
		XL_readNumVectorAndSize( XL_X, C_XRows, C_XCols, C_X, " ARM_ERR: Matrix: matrix expected", C_result);

		VECTOR<double> C_XInter;
		long C_XInterRows, C_XCInterCols;
		XL_readNumVectorAndSize( XL_XInter, C_XInterRows, C_XCInterCols, C_XInter, " ARM_ERR: Matrix: matrix expected", C_result);

		CCString C_RegModeCCStr;
		XL_readStrCell(XL_RegMode,C_RegModeCCStr," ARM_ERR: RegMode: String expected",C_result);
		string C_RegMode = CCSTringToSTLString(C_RegModeCCStr);

		double spanDef=0.8;
		double C_Span;
		XL_readNumCellWD( XL_Span, C_Span, spanDef, " ARM_ERR: Span: numeric expected", C_result);

		VECTOR<double> C_OutputVector;

		long retCode = ARMLOCAL_GP_Regression(
			C_Y,
			C_X,
			C_XRows,
			C_XCols,
			C_XInter,
			C_XInterRows,
			C_XCInterCols,
			C_RegMode,
			C_Span,
			C_OutputVector,
			C_result);

		if (retCode == ARM_OK)
		{
			XL_writeNumVector( XL_result, C_OutputVector, "ARM ERR: could not write the output", C_result );
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_GP_Regression" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}