/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file ARM_xl_gp_pricer_local.cpp
 *
 *  \brief file for the pricer part in the generic pricer
 *
 *	\author  E Benhamou
 *	\version 1.0
 *	\date September 2003
 */

#include <ARM\libarm_local\firstToBeIncluded.h>
#include <functional>
#include <libCCxll\CCxll.h>

#include <ARM\libarm_local\ARM_local_gp_pricer.h>
#include "ARM_xl_gp_gensecurity_local.h"
#include "ARM_xl_gp_common.h"
#include "ARM_xl_wrapper_local.h"
#include <GP_Infra\gpinfra\gramfunctorarg.h>
#include "ARM_xl_trycatch_local.h"
#include "ARM_xl_gp_fctorhelper.h"
#include <gpbase\gpmatrix.h>

#include "util\tech_macro.h"

using ARM::ARM_GramFctorArg;
using ARM::GFAT_DOUBLE_TYPE;
using ARM::GFAT_VECTOR_TYPE;
using ARM::GFAT_MATRIX_TYPE;
using ARM::GFAT_STRING_TYPE;
using ARM::GFAT_DATE_TYPE;

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
class GenPricerFctor : public ARMResultLong2LongFunc
{
public:
	GenPricerFctor(
		long genSecurityId,
		long pricingModelId,
		const vector<string>& columnNames,
		const vector<double>& columnPrices,
		const string& refColumn,
		const vector<double>& betas )
	:
		C_genSecurityId( genSecurityId ),
		C_pricingModelId( pricingModelId ),
		C_ColumnNames( columnNames ),
		C_ColumnPrices( columnPrices ),
		C_RefColumn( refColumn ),
		C_Betas(betas)
	{};
	
	long operator()( ARM_result& result, long objId ) 
	{
		return ARMLOCAL_GenPricer_Create(
			C_genSecurityId,
			C_pricingModelId,
			C_ColumnNames, 
			C_ColumnPrices,
			C_RefColumn,
			C_Betas,
			result,
			objId );			
	}
private:
	long C_genSecurityId;
	long C_pricingModelId;
	vector<string> C_ColumnNames;
	vector<double> C_ColumnPrices;
	string C_RefColumn;
	vector<double> C_Betas;
};



/////////////////////////////////////////////////////////////
/// central function that does the creation of the XL function
/////////////////////////////////////////////////////////////
LPXLOPER Local_GenPricer_Common(
	LPXLOPER XL_GenSecurityId,
	LPXLOPER XL_PricingModelId,
	LPXLOPER XL_columnNames,
	LPXLOPER XL_columnPrices,
	LPXLOPER XL_refColumn,
	LPXLOPER XL_betas,
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

		long retCode = ARM_OK;

		long C_genSecurityId;
		XL_GETOBJID( XL_GenSecurityId,	C_genSecurityId,	" ARM_ERR: Generic Security id: object expected",			C_result);

		long C_pricingModelId;
		XL_GETOBJID( XL_PricingModelId,	C_pricingModelId,	" ARM_ERR: Pricing Model  id: object expected",			C_result);

		VECTOR<CCString> C_columnNamesVecCCStr;
		VECTOR<CCString> C_columnNamesVecCCStrDef;
		XL_readStrVectorWD(XL_columnNames,C_columnNamesVecCCStr,C_columnNamesVecCCStrDef," ARM_ERR: column anmes: array of string expected", DOUBLE_TYPE, C_result);
		vector<string> C_columnNames(C_columnNamesVecCCStr.size());
		for( size_t i=0; i<C_columnNamesVecCCStr.size(); ++i )
			C_columnNames[i] = CCSTringToSTLString(C_columnNamesVecCCStr[i]);

		vector<double> C_columnPrices;
		vector<double> C_columnPricesDefault;
		XL_readNumVectorWD(XL_columnPrices,C_columnPrices,C_columnPricesDefault, " ARM_ERR: columnPrices: numeric vector expected", C_result);

		CCString C_refColumnCCStr;
		CCString  C_refColumnCCStrDef("");
		XL_readStrCellWD( XL_refColumn, C_refColumnCCStr, C_refColumnCCStrDef, " ARM_ERR: refColumn: String expected",C_result);
		string C_refColumn = CCSTringToSTLString( C_refColumnCCStr );

		vector<double> C_Betas;
		vector<double> C_BetasDefault;
		XL_readNumVectorWD(XL_betas,C_Betas,C_BetasDefault, " ARM_ERR: Betas: numeric vector expected", C_result);

		GenPricerFctor ourFunc( C_genSecurityId, C_pricingModelId, C_columnNames, C_columnPrices, C_refColumn, C_Betas );

		fillXL_Result( LOCAL_GENPRICER_CLASS, ourFunc, C_result, XL_result, PersistentInXL );
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_GenPricer_Common" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}
							 
							 
///////////////////////////////////
/// version that takes into account 
/// previous creation of object
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_GenPricer_Create(
	LPXLOPER XL_GenSecurityId,
	LPXLOPER XL_PricingModelId,
	LPXLOPER XL_columnNames,
	LPXLOPER XL_columnPrices,
	LPXLOPER XL_refColumn,
	LPXLOPER XL_betas )
{
	ADD_LOG("Local_GenPricer_Create");
	bool PersistentInXL = true;
	return Local_GenPricer_Common(
		XL_GenSecurityId,
		XL_PricingModelId,
		XL_columnNames,
		XL_columnPrices,
		XL_refColumn,
		XL_betas,
		PersistentInXL );
}


///////////////////////////////////
/// version for non persistentInXl
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_PXL_GenPricer_Create(
	LPXLOPER XL_GenSecurityId,
	LPXLOPER XL_PricingModelId,
	LPXLOPER XL_columnNames,
	LPXLOPER XL_columnPrices,
	LPXLOPER XL_refColumn,
	LPXLOPER XL_betas )
{
	ADD_LOG("Local_PXL_GenPricer_Create");
	bool PersistentInXL = false;
	return Local_GenPricer_Common(
		XL_GenSecurityId,
		XL_PricingModelId,
		XL_columnNames,
		XL_columnPrices,
		XL_refColumn,
		XL_betas,
		PersistentInXL );
}


///----------------------------------------------
///----------------------------------------------
///             Local Variance Function
/// Inputs :
///     Model Id
///     fromDate, toDate
///     startDate, endDate
///----------------------------------------------
///----------------------------------------------

__declspec(dllexport) LPXLOPER WINAPI Local_GetData_FromGenPricer_Common(
	LPXLOPER XL_PricerId,
    LPXLOPER XL_KeyId, 
	LPXLOPER XL_ColumnName,
	LPXLOPER XL_Index,
	bool fromVB )
{
	ADD_LOG("Local_GetData_FromGenPricer_Common");
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
		
		CCString C_PricerString;
		XL_readStrCell( XL_PricerId, C_PricerString, " ARM_ERR: Pricer Id: Object expected",C_result);
		long C_PricerId = LocalGetNumObjectId(C_PricerString);
		
		CCString C_KeyString;
		CCString C_KeyStringDef;
		XL_readStrCellWD( XL_KeyId, C_KeyString, C_KeyStringDef, " ARM_ERR: Key Id: String expected",C_result);

		CCString C_ColumnName;
		XL_readStrCellWD( XL_ColumnName, C_ColumnName, "", " ARM_ERR: Column Name: String expected",C_result);

		double C_IndexDble;
		double C_IndexDefault = 0;
		XL_readNumCellWD( XL_Index, C_IndexDble, C_IndexDefault, " ARM_ERR: Index expected",C_result);
		int  C_Index = floor(C_IndexDble);


		
		/// call the function
		ARM_GramFctorArg argResult;
		
		long retCode = ARMLOCAL_GenPricer_GetData(
			C_PricerId,
			CCSTringToSTLString( C_KeyString ),
			CCSTringToSTLString( C_ColumnName ),
			argResult,
			C_result );
		
		if( retCode == ARM_KO )
		{
			ARM_ERR();
		}
		else
		{
			ARM_GramFunctorToXLOPER(argResult, XL_result, C_result, fromVB, C_Index);
			FreeCurCellErr ();
		}
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_GetData_FromGenPricer" )

	return (LPXLOPER)&XL_result;
}


__declspec(dllexport) LPXLOPER WINAPI Local_GetData_FromGenPricer(
	LPXLOPER XL_PricerId,
    LPXLOPER XL_KeyId,
	LPXLOPER XL_ColumnName,
	LPXLOPER XL_Index)
{
	ADD_LOG("Local_GetData_FromGenPricer");
	return Local_GetData_FromGenPricer_Common( XL_PricerId, XL_KeyId, XL_ColumnName, XL_Index, true );
}


__declspec(dllexport) LPXLOPER WINAPI Local_PXL_GetData_FromGenPricer(
	LPXLOPER XL_PricerId,
    LPXLOPER XL_KeyId,
	LPXLOPER XL_ColumnName,
	LPXLOPER XL_Index)
{
	ADD_LOG("Local_PXL_GetData_FromGenPricer");
	return Local_GetData_FromGenPricer_Common( XL_PricerId, XL_KeyId, XL_ColumnName, XL_Index, false );
}



/////////////////////////////////////////////////////////////
/// central function that does the creation of the XL function
/////////////////////////////////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_GenPricer_SetDetailMode(
	LPXLOPER XL_PricerId,
	LPXLOPER XL_detailMode )
{
	ADD_LOG("Local_GenPricer_SetDetailMode");
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

		long C_pricerId;
		XL_GETOBJID( XL_PricerId, C_pricerId,	" ARM_ERR: Pricer id: object expected",	C_result);

		double C_detailMode;
		XL_readNumCell( XL_detailMode, C_detailMode,	" ARM_ERR: detail mode Flag: boolean compatible expected",	C_result);
		bool C_detailModeFlag = C_detailMode != 0;
			
		long retCode = ARMLOCAL_GenPricer_SetDetailMode(
				C_pricerId,
				C_detailModeFlag,
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_GenPricer_SetDetailMode" )

	return (LPXLOPER)&XL_result;
}


/////////////////////////////////////////////////////////////
/// central function that does the creation of a GenSecurity
/// object from the std ARM security. 
/////////////////////////////////////////////////////////////
LPXLOPER Local_ChangeSecurityIntoGenSec_Common(
	LPXLOPER XL_SecurityId, LPXLOPER XL_AsOfDate, LPXLOPER XL_ModelName, 
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

		long retCode = ARM_OK;

		long C_SecurityId;
		double AsOfDate;
		XL_readNumCell(XL_AsOfDate,AsOfDate," ARM_ERR: AsOfDate : date expected",C_result);


		CCString C_ModelName;
		XL_readStrCell(XL_ModelName,C_ModelName," ARM_ERR: ModelName Not of a good Type! ",C_result);
		string modelName = CCSTringToSTLString(C_ModelName);

		XL_GETOBJID( XL_SecurityId,	C_SecurityId,	" ARM_ERR: Security id: object expected",			C_result);
		exportFunc3Args< long, double, string >  ourFunc( C_SecurityId, AsOfDate, modelName, ARMLOCAL_ChangeSecurityIntoGenSec);
		fillXL_Result( LOCAL_GENSEC_CLASS, ourFunc, C_result, XL_result, PersistentInXL );
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_DealDesAndGenSecCommon" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


///////////////////////////////////
/// Creates a generic security object
/// from Another Generic security
///	(changes a bermudan option into a trigger
/// using its exercise boundary)
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_ChangeSecurityIntoGenSec(
	LPXLOPER XL_SecurityId, LPXLOPER XL_AsOfDate, LPXLOPER XL_ModelName)
{
	ADD_LOG("Local_ChangeSecurityIntoGenSec");
	bool PersistentInXL = true;
	return Local_ChangeSecurityIntoGenSec_Common( XL_SecurityId, XL_AsOfDate,XL_ModelName, PersistentInXL );
}

///////////////////////////////////
/// Creates a generic security object
/// from another generic security
///	used for call in VBA
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_PXL_ChangeSecurityIntoGenSec(
	LPXLOPER XL_SecurityId, LPXLOPER XL_AsOfDate, LPXLOPER XL_ModelName)
{
	ADD_LOG("Local_PXL_ChangeSecurityIntoGenSec");
	bool PersistentInXL = false;
	return Local_ChangeSecurityIntoGenSec_Common( XL_SecurityId, XL_AsOfDate,XL_ModelName, PersistentInXL );
}