/*
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 * $Log: ARM_xl_gp_gensecurity_local.cpp,v $
 * Revision 1.1  2003/13/07 15:08:43  ebenhamou
 * Initial version
 *
 */

/*! \file ARM_xl_gp_gensecurity_local.cpp
 *
 *  \brief file for the generic security part in the generic pricer
 *
 *	\author  E Benhamou
 *	\version 1.0
 *	\date September 2003
 */


#include <ARM\libarm_local\firstToBeIncluded.h>
#include <functional>
#include <libCCxll\CCxll.h>

#include <ARM\libarm_local\ARM_local_gp_gensecurity.h>
#include "ARM_xl_gp_gensecurity_local.h"
#include "ARM_xl_gp_fctorhelper.h"
#include "ARM_xl_wrapper_local.h"
#include "ARM_xl_trycatch_local.h"
#include "ARM_gp_local_interglob.h"
#include "ARM_gp_local_excelManip.h"

/// to debug release only crash
#include "gpbase/eventviewerfwd.h"

#include "util\tech_macro.h"

/// to debug release only crash
using ARM::ARM_EventViewerImp;
using ARM::ARM_TheEventViewer;

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

typedef long (*DealDesAndGenSecCtor_Function)(
	const VECTOR<string>&,
	const VECTOR<ARM_GP_VALUE_TYPE>&,
	long, 
	long,
    const string&,
	long,
	bool,
    bool,
	const VECTOR<string>&,
	bool,
	ARM_result&, 
	long);


////////////////////////////////////////////////
/// very rapid definition of ourlonglongFunctor
////////////////////////////////////////////////
class DealDesAndGenSecCtor : public ARMResultLong2LongFunc
{
public:
	DealDesAndGenSecCtor(
		const VECTOR<string>&  vecString ,
		const VECTOR<ARM_GP_VALUE_TYPE>& types,
		long rowsNb,
		long colsNb,
        const string& payCurveName,
		long cstmanagerId,
		bool ExercBoundaryResetFlag,
		bool OtherPayoffsFlag,
		const VECTOR<string>&  pricedColumns,
		bool IVFlag,
		const DealDesAndGenSecCtor_Function& f)
	:
		C_vecString( vecString ),
		C_types( types ),
		C_rowsNb( rowsNb ),
		C_colsNb( colsNb ),
        C_payCurveName(payCurveName),
		C_cstmanagerId(cstmanagerId),
		C_ExercBoundaryResetFlag(ExercBoundaryResetFlag),
        C_OtherPayoffsFlag(OtherPayoffsFlag),
		C_PricedColumns(pricedColumns),
		C_IVFlag(IVFlag),
		itsFunction(f)
	{};
	
	long operator()( ARM_result& result, long objId ) 
	{
		return (*itsFunction)(
			C_vecString,
			C_types,			
			C_rowsNb,
			C_colsNb,
            C_payCurveName,
			C_cstmanagerId,
			C_ExercBoundaryResetFlag,
            C_OtherPayoffsFlag,
			C_PricedColumns,
			C_IVFlag,
			result,
			objId );			
	}
private:
	VECTOR<string> C_vecString;
	VECTOR<ARM_GP_VALUE_TYPE> C_types;
	long C_rowsNb;
	long C_colsNb;
    string C_payCurveName;
	long C_cstmanagerId;
	bool C_ExercBoundaryResetFlag;
    bool C_OtherPayoffsFlag;
	bool C_IVFlag;
	VECTOR<string> C_PricedColumns;
	DealDesAndGenSecCtor_Function itsFunction;
};


/////////////////////////////////////////////////////////////
/// central function that does the creation of the XL function
/////////////////////////////////////////////////////////////
LPXLOPER Local_DealDesAndGenSecCommon(
	LPXLOPER XL_Data1,
	LPXLOPER XL_Data2,
	LPXLOPER XL_PayCurveName,
	LPXLOPER XL_cstManagerId,
	LPXLOPER XL_ExercBoundaryResetFlag,
    LPXLOPER XL_OtherPayoffsFlag,
	LPXLOPER XL_PricedColumns,
	LPXLOPER XL_IVFlag,
	const DealDesAndGenSecCtor_Function& function,
	const CCString& curClass,
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

		VECTOR<CCString> C_vecCCString1, C_vecCCString2,C_vecCCStringDefault;
		VECTOR<long> C_types1, C_types2,C_typesDefault;
		long C_rowsNb1 = 0, C_colsNb1 = 0,
			 C_rowsNb2 = 0, C_colsNb2 = 0;
		
		/// this is used by macros 
		/// and therefore this has to be defined
		static int error;
		static char* reason = "";

		XL_readStrVectorSizeAndType( XL_Data1,C_rowsNb1,C_colsNb1,C_vecCCString1,C_types1, " ARM_ERR: deal description: matrix expected", DOUBLE_TYPE,C_result);
		XL_readStrVectorSizeAndTypeWD( XL_Data2,C_rowsNb2,C_colsNb2,C_vecCCString2,&C_vecCCStringDefault,C_types2,&C_typesDefault, " ARM_ERR: deal description: matrix expected", DOUBLE_TYPE,C_result);

		/// FIX FIX discounting curve ignored at this stage!
		CCString payCurveNameXl;
		XL_readStrCellWD(XL_PayCurveName,payCurveNameXl,""," ARM_ERR: discount curve name: string expected",C_result);
        string payCurveName = CCSTringToSTLString(payCurveNameXl);

		long C_cstManagerId;
		XL_GETOBJIDWD( XL_cstManagerId,	C_cstManagerId,	"NULL OBJECT", " ARM_ERR: Cst MAnager id: object expected",C_result);

		/// ExercBoundaryResetFlag
		CCString ExercBoundaryResetFlagXl;
		XL_readStrCellWD( XL_ExercBoundaryResetFlag,ExercBoundaryResetFlagXl,"Y"," ARM_ERR: ExerciseBoundaryResetFlag: string expected",C_result);
		bool ExercBoundaryResetFlag;
		ExercBoundaryResetFlagXl.toUpper();
		if (ExercBoundaryResetFlagXl == "Y" )
			ExercBoundaryResetFlag = true;
		else if (ExercBoundaryResetFlagXl == "N")
			ExercBoundaryResetFlag = false;
		else
			throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"String \"Y\" or \"N\" Expected for ExercBoundaryResetFlag");


		/// OtherPayoffFlag (to avoid IntermediatePayoffs & Snapshots computing)
		CCString OtherPayoffsFlagXl;
		XL_readStrCellWD( XL_OtherPayoffsFlag,OtherPayoffsFlagXl,"Y"," ARM_ERR: OtherPayoffsFlag: string expected",C_result);
		bool OtherPayoffsFlag;
		OtherPayoffsFlagXl.toUpper();
		if (OtherPayoffsFlagXl == "Y" )
			OtherPayoffsFlag = true;
		else if (OtherPayoffsFlagXl == "N")
			OtherPayoffsFlag = false;
		else
			throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"String \"Y\" or \"N\" Expected for OtherPayoffsFlag");

		/// IVFlag (to avoid IntermediateValues)
		CCString IVFlagXl;
		XL_readStrCellWD( XL_IVFlag,IVFlagXl,"Y"," ARM_ERR: IVFlag: string expected",C_result);
		bool IVFlag;
		IVFlagXl.toUpper();
		if (IVFlagXl == "Y" )
			IVFlag = true;
		else if (IVFlagXl == "N")
			IVFlag = false;
		else
			throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"String \"Y\" or \"N\" Expected for IVFlag");

		VECTOR<CCString> C_CCpricedColumns;
		VECTOR<CCString> pricedColumnsDef;
		XL_readStrVectorWD(XL_PricedColumns,C_CCpricedColumns,pricedColumnsDef," ARM_ERR: priced columns: array of string expected",DOUBLE_TYPE,C_result);

		long retCode = ARM_OK;

		// If the second cell contains a new array with the same number of columns, we add it
		// on the bottom of the fist one.
		if( C_rowsNb2 && C_colsNb2 )
		{
			if (C_colsNb1 == C_colsNb2)
			{
				int i = 0;

				for (i = 0; i < C_rowsNb2*C_colsNb2; ++i)
				{
					C_vecCCString1.push_back(C_vecCCString2[i]);
					C_types1.push_back(C_types2[i]);
				}
				C_rowsNb1 += C_rowsNb2;
			}
			else
			{
				CCString local_msg ("The two data table must have the same number of columns.");
				C_result.setMsg (local_msg);
				retCode = ARM_KO;
			}
		}

		if(retCode==ARM_OK)
		{
			/// conversion from CCString to string
			vector<string> vecString( C_vecCCString1.size() );
			int i;
			for( i=0; i<C_vecCCString1.size(); ++i )
				vecString[i] = C_vecCCString1[i];

			/// conversion to ARM_ValueType
			vector<ARM_GP_VALUE_TYPE> vecTypes( C_types1.size() );
			for( i=0; i<vecTypes.size(); ++i )
				vecTypes[i] = ConvertToGenPricerType( ConvertXLtoARMType( C_types1[i] ), vecString[i] );

			vector<string> C_pricedColumns( C_CCpricedColumns.size() );

			for( i=0; i<C_CCpricedColumns.size(); ++i )
				C_pricedColumns[i] = C_CCpricedColumns[i];


			/// use the concept of Functor to transfer the knowledge of
			/// a function with a context
			DealDesAndGenSecCtor ourFunc(
				vecString,
				vecTypes,
				C_rowsNb1,
				C_colsNb1,
                payCurveName,
				C_cstManagerId,
				ExercBoundaryResetFlag,
                OtherPayoffsFlag,
				C_pricedColumns,
				IVFlag,
				function );

			/// call the general function
			fillXL_Result( curClass, ourFunc, C_result, XL_result, PersistentInXL );
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_DealDesAndGenSecCommon" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


/////////////////////////////////////////////////////////////
/// central function that does the creation of the new XL 
/// function from the old one. 
/////////////////////////////////////////////////////////////
LPXLOPER Local_GenSec_ChangeAmericanIntoTrigger_Common(
	LPXLOPER XL_GenSecurityId, LPXLOPER XL_PricingModelId, 
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
		long C_PricingModelId;

		XL_GETOBJID( XL_GenSecurityId,	C_genSecurityId,	" ARM_ERR: Generic Security id: object expected",			C_result);
		XL_GETOBJID( XL_PricingModelId,	C_PricingModelId,	" ARM_ERR: Pricing Model id: object expected",			C_result);

		exportFunc2Args< long, long >  ourFunc( C_genSecurityId, C_PricingModelId, ARMLOCAL_GenSec_ChangeAmericanIntoTrigger_Common );

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

/////////////////////////////////////////////////////////////
/// central function that does the creation of the new XL 
/// function from the old one. 
/////////////////////////////////////////////////////////////
LPXLOPER Local_GenSec_ChangeMAXPVIntoExercise_Common(
	LPXLOPER XL_GenSecurityId, 
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

		XL_GETOBJID( XL_GenSecurityId, C_genSecurityId,	" ARM_ERR: Generic Security id: object expected",			C_result);

		exportFunc1Arg< long >  ourFunc( C_genSecurityId, ARMLOCAL_GenSec_ChangeMAXPVIntoExercise_Common );

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
/// version that takes into account 
/// previous creation of object
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_DealDes_Create(
	LPXLOPER XL_Data1,
	LPXLOPER XL_Data2,
	LPXLOPER XL_PricedColumns)
{
	ADD_LOG("Local_DealDes_Create");
	bool PersistentInXL = true;
	return Local_DealDesAndGenSecCommon(
		XL_Data1,
		XL_Data2,
		NULL,
		NULL,
		NULL,
        NULL,
		XL_PricedColumns,
		NULL,
		&ARMLOCAL_DealDes_Create,
		LOCAL_DEALDES_CLASS,
		PersistentInXL );
}


///////////////////////////////////
/// version for non persistentInXl
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_PXL_DealDes_Create(
	LPXLOPER XL_Data1,
	LPXLOPER XL_Data2,
	LPXLOPER XL_PricedColumns)
{
	ADD_LOG("Local_PXL_DealDes_Create");
	bool PersistentInXL = false;
	return Local_DealDesAndGenSecCommon(
		XL_Data1,
		XL_Data2,
		NULL,
		NULL,
		NULL,
        NULL,
		XL_PricedColumns,
		NULL,
		&ARMLOCAL_DealDes_Create,
		LOCAL_DEALDES_CLASS,
		PersistentInXL );
}





///////////////////////////////////
/// version that takes into account 
/// previous creation of object
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_GenSec_Create(
	LPXLOPER XL_Data1,
	LPXLOPER XL_Data2,
	LPXLOPER XL_PayCurveName,
	LPXLOPER XL_CstManagerId,
	LPXLOPER XL_ExercBoundaryResetFlag,
    LPXLOPER XL_OtherPayoffsFlag,
	LPXLOPER XL_PricedColumns,
	LPXLOPER XL_IVFlag)
{
	ADD_LOG("Local_GenSec_Create");
	bool PersistentInXL = true;
	return Local_DealDesAndGenSecCommon(
		XL_Data1,
		XL_Data2,
		XL_PayCurveName,
		XL_CstManagerId,
		XL_ExercBoundaryResetFlag,
        XL_OtherPayoffsFlag,
		XL_PricedColumns,
		XL_IVFlag,
		&ARMLOCAL_GenSec_Create,
		LOCAL_GENSEC_CLASS,
		PersistentInXL );
}


///////////////////////////////////
/// version for non persistentInXl
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_PXL_GenSec_Create(
	LPXLOPER XL_Data1,
	LPXLOPER XL_Data2,
	LPXLOPER XL_PayCurveName,
	LPXLOPER XL_CstManagerId,
	LPXLOPER XL_ExercBoundaryResetFlag,
    LPXLOPER XL_OtherPayoffsFlag,
	LPXLOPER XL_PricedColumns,
	LPXLOPER XL_IVFlag)
{
	ADD_LOG("Local_PXL_GenSec_Create");
	bool PersistentInXL = false;
	return Local_DealDesAndGenSecCommon(
		XL_Data1,
		XL_Data2,
		XL_PayCurveName,
		XL_CstManagerId,
		XL_ExercBoundaryResetFlag,
        XL_OtherPayoffsFlag,
		XL_PricedColumns,
		XL_IVFlag,
		&ARMLOCAL_GenSec_Create,
		LOCAL_GENSEC_CLASS,
		PersistentInXL );
}

///////////////////////////////////
/// Creates a generic security object
/// from Another Generic security
///	(changes a bermudan option into a trigger
/// using its exercise boundary)
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_GenSec_ChangeAmericanIntoTrigger(
	LPXLOPER XL_GenSecurityId, LPXLOPER XL_PricingModelId )
{
	ADD_LOG("Local_GenSec_ChangeAmericanIntoTrigger");
	bool PersistentInXL = true;
	return Local_GenSec_ChangeAmericanIntoTrigger_Common( XL_GenSecurityId, XL_PricingModelId, PersistentInXL );
}

///////////////////////////////////
/// Creates a generic security object
/// from another generic security
///	used for call in VBA
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_PXL_GenSec_ChangeAmericanIntoTrigger(
	LPXLOPER XL_GenSecurityId, LPXLOPER XL_PricingModelId )
{
	ADD_LOG("Local_PXL_GenSec_ChangeAmericanIntoTrigger");
	bool PersistentInXL = false;
	return Local_GenSec_ChangeAmericanIntoTrigger_Common( XL_GenSecurityId, XL_PricingModelId, PersistentInXL );
}

///////////////////////////////////
/// Creates a generic security object
/// from Another Generic security
///	(changes MAX(PV) keywords into 
/// exericse)
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_GenSec_ChangeMAXPVIntoExercise(
	LPXLOPER XL_GenSecurityId )
{
	ADD_LOG("Local_GenSec_ChangeMAXPVIntoExercise");
	bool PersistentInXL = true;
	return Local_GenSec_ChangeMAXPVIntoExercise_Common( XL_GenSecurityId, PersistentInXL );
}

///////////////////////////////////
/// Creates a generic security object
/// from another generic security
///	used for call in VBA
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_PXL_GenSec_ChangeMAXPVIntoExercise(
	LPXLOPER XL_GenSecurityId )
{
	ADD_LOG("Local_PXL_GenSec_ChangeMAXPVIntoExercise");
	bool PersistentInXL = false;
	return Local_GenSec_ChangeMAXPVIntoExercise_Common( XL_GenSecurityId, PersistentInXL );
}


/////////////////////////////////////////////////////////////
/// central function that does the creation of the XL function
/////////////////////////////////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_GenSec_SetPTFlag(
	LPXLOPER XL_GenSec,
	LPXLOPER XL_PTFlag )
{
	ADD_LOG("Local_GenSec_SetPTFlag");
	
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

		double C_Flag;
		long C_GenSecId;
		XL_GETOBJID( XL_GenSec,	C_GenSecId,	" ARM_ERR: Generic Security id: object expected",			C_result);
		XL_readNumCell( XL_PTFlag, C_Flag,	" ARM_ERR: Parse Tree Flag: boolean compatible expected",	C_result);
		bool ParseTreeFlag = C_Flag != 0;
			
		long retCode = ARMLOCAL_GenSec_SetParseTreeFlag(
				C_GenSecId,
				ParseTreeFlag,
				C_result );

		/// feed the LPXLOPER object result 
		if (retCode == ARM_OK)
		{
			FreeCurCellErr ();
			XL_result.xltype = xltypeStr;
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG("Unrecognized failure in Local_GenSec_SetPTFlag" )


	return (LPXLOPER)&XL_result;
}
							 



////////////////////////////////////////////////
/// very rapid definition of ourlonglongFunctor
////////////////////////////////////////////////
class GramHelperFctor: public ARMResultLong2LongFunc
{
public:
	GramHelperFctor(const string& FuncName)
	:	C_FuncName(FuncName){};

	long operator()( ARM_result& result, long objId ){
		return ARMLOCAL_GramHelper(
			C_FuncName,
			result,
			objId );			
	}
private:
	string C_FuncName;
};


/////////////////////////////////////////////////////////////
/// central function that does the creation of the XL function
/////////////////////////////////////////////////////////////
LPXLOPER Local_GramHelperCommon(
	LPXLOPER XL_FuncName,
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

		CCString C_FuncName;
		const char FuncNameDefault[8]("ALLFUNC");
		
		/// this is used by macros 
		/// and therefore this has to be defined
		static int error;
		static char* reason = "";
		
		XL_readStrCellWD(	 XL_FuncName,	C_FuncName,	FuncNameDefault, " ARM_ERR: Function Name string expected",	C_result);

		/// use the concept of Functor to transfer the knowledge of
		/// a function with a context
		char* funcNameChar = C_FuncName.c_str();
		string FuncNameString( funcNameChar );
		delete funcNameChar;
		GramHelperFctor ourFunc( FuncNameString );

		/// call the general function
		fillXL_Result( LOCAL_GRAMHELPER_CLASS, ourFunc, C_result, XL_result, PersistentInXL );

	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG("Unrecognized failure in Local_GramHelperCommon" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}
							 
							 
///////////////////////////////////
/// version that takes into account 
/// previous creation of object
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_GramHelper_Create(
	LPXLOPER XL_FuncName )
{
	ADD_LOG("Local_GramHelper_Create");
	bool PersistentInXL = true;
	return Local_GramHelperCommon(
		XL_FuncName,
		PersistentInXL );
}


///////////////////////////////////
/// version for non persistentInXl
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_PXL_GramHelper_Create(
	LPXLOPER XL_FuncName )
{
	ADD_LOG("Local_PXL_GramHelper_Create");
	bool PersistentInXL = false;
	return Local_GramHelperCommon(
		XL_FuncName,
		PersistentInXL );
}



/////////////////////////////////////////////////////////////
/// central function that does the creation of the XL function
/////////////////////////////////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_GenSec_SetCRFlag(
	LPXLOPER XL_CRFlag )
{
	ADD_LOG("Local_GenSec_SetCRFlag");
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

		double C_Flag;
		XL_readNumCell( XL_CRFlag, C_Flag,	" ARM_ERR: Parse Tree Flag: boolean compatible expected",	C_result);
		bool CircularRefFlag = C_Flag != 0;
		
		long retCode;
		if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 1)
			retCode = ARMLOCAL_GenSec_SetCircularRefFlag(
				CircularRefFlag,
				C_result );
		else
			retCode = ARM_KO;
		
		/// feed the LPXLOPER object result 
		if (retCode == ARM_OK)
		{
			FreeCurCellErr ();
			XL_result.xltype = xltypeStr;
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG("Unrecognized failure in Local_GenSec_SetCRFlag" )

	return (LPXLOPER)&XL_result;
}

/////////////////////////////////////////////////////////////
/// central function that does the creation of the XL function
/////////////////////////////////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_GenSec_GetDealDesTable(
	LPXLOPER XL_GenSec)
{
	ADD_LOG("Local_GenSec_GetDealDesTable");
	
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

		long C_GenSecId;
		XL_GETOBJID( XL_GenSec,	C_GenSecId,	" ARM_ERR: Generic Security id: object expected",C_result);

		vector<string> dealDesValues;
		vector<ARM_GP_VALUE_TYPE> dealDesTypes;
		long nbRows,nbCols;
		long retCode = ARMLOCAL_GenSec_GetDealDesTable(
				C_GenSecId,
				dealDesValues,
				dealDesTypes,
				nbRows, nbCols,
				C_result );

		/// conversion to XL types
		VECTOR<CCString> C_DealDesValues(dealDesValues.size());
		int i;
		for( i=0; i<dealDesValues.size(); ++i )
		{
			if(dealDesTypes[i] == ARM_DATE_TYPE)
			{
				/// Convert julian to Xl date
				double xlDate = JulianToXLDate(atof(dealDesValues[i].c_str()));
				char xlDateStr[30];
				sprintf(xlDateStr,"%lf",xlDate);
				C_DealDesValues[i] = xlDateStr;
			}
			else
				C_DealDesValues[i] = dealDesValues[i].c_str();
		}

		VECTOR<long> C_DealDesTypes(dealDesTypes.size());
		for( i=0; i<dealDesTypes.size(); ++i )
			C_DealDesTypes[i] = ConvertARMToXLType(ConvertFromGenPricerType(dealDesTypes[i]));

		/// feed the LPXLOPER object result 
		if (retCode == ARM_OK)
		{	const int additionalLinesNb =100;
			bool fillWithBlank = true;
			FreeCurCellContent ();
			XL_writeStrMatrixSizeAndTypeWithOptions( XL_result, C_DealDesValues, C_DealDesTypes, nbRows, nbCols, " ARM_ERR: Could not set the string matrix", C_result,additionalLinesNb,fillWithBlank );
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG("Unrecognized failure in Local_GenSec_GetDealDesTable" )
	
	return (LPXLOPER)&XL_result;
}




///----------------------------------------------
///----------------------------------------------
///             cst manager
/// Inputs :
///     vector of names
///     vector of values
///----------------------------------------------
///----------------------------------------------

/////////////////////////////////////////////////////////////
/// central function that does the creation of the XL function
/////////////////////////////////////////////////////////////
LPXLOPER Local_CstManager_Create_Common(
	LPXLOPER XL_names,
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

		VECTOR<CCString> C_namesCCStr;
	    XL_readStrVector(XL_names,C_namesCCStr," ARM_ERR: names: array of string expected", XL_TYPE_STRING, C_result);

		VECTOR<double> C_values;
		XL_readNumVector(XL_values,C_values," ARM_ERR: values: array of numeric expected",C_result);

		/// conversion to vector<string>
		vector<string> C_Names( C_namesCCStr.size() );
		for(size_t i=0; i<C_namesCCStr.size() ; ++i )
			C_Names[i] = C_namesCCStr[i];

		/// use the concept of Functor to transfer the knowledge of
		/// a function with a context
		exportFunc2Args< VECTOR<string>, VECTOR<double> >  ourFunc(C_Names,C_values, ARMLOCAL_CstManager_Create );

		/// call the general function
		fillXL_Result( LOCAL_CSTMANAGER_CLASS, ourFunc, C_result, XL_result, PersistentInXL );
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_CstManager_Create_Common" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}
							 

///////////////////////////////////
/// version that takes into account 
/// previous creation of object
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_CstManager_Create(
	LPXLOPER XL_names,
	LPXLOPER XL_values )
{
	ADD_LOG("Local_CstManager_Create");
	bool PersistentInXL = true;
	return Local_CstManager_Create_Common(
		XL_names,
		XL_values,
		PersistentInXL );
}


///////////////////////////////////
/// version for non persistentInXl
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_PXL_CstManager_Create(
	LPXLOPER XL_names,
	LPXLOPER XL_values)
{
	ADD_LOG("Local_PXL_CstManager_Create");
	bool PersistentInXL = false;
	return Local_CstManager_Create_Common(
		XL_names,
		XL_values,
		PersistentInXL );
}

/////////////////////////////////////////////////////////////
/// central function that gets the cst manager from a generic security
/////////////////////////////////////////////////////////////
LPXLOPER Local_GenSec_GetCstManager_Common(
	LPXLOPER XL_GenSec,
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

		long C_GenSecId;
		XL_GETOBJID( XL_GenSec,	C_GenSecId,	" ARM_ERR: Generic Security id: object expected",C_result);

		/// use the concept of Functor to transfer the knowledge of
		/// a function with a context
		exportFunc1Arg< long >  ourFunc(C_GenSecId, ARMLOCAL_GenSec_GetCstManager );

		/// call the general function
		fillXL_Result( LOCAL_CSTMANAGER_CLASS, ourFunc, C_result, XL_result, PersistentInXL );
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_CstManager_Create_Common" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}
							 

///////////////////////////////////
/// version that takes into account 
/// previous creation of object
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_GenSec_GetCstManager(
	LPXLOPER XL_GenSec)
{
	ADD_LOG("Local_GenSec_GetCstManager");
	bool PersistentInXL = true;
	return Local_GenSec_GetCstManager_Common(
		XL_GenSec,
		PersistentInXL );
}


///////////////////////////////////
/// version for non persistentInXl
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_PXL_GenSec_GetCstManager(
	LPXLOPER XL_GenSec)
{
	ADD_LOG("Local_PXL_GenSec_GetCstManager");
	bool PersistentInXL = false;
	return Local_GenSec_GetCstManager_Common(
		XL_GenSec,
		PersistentInXL );
}

///----------------------------------------------
///----------------------------------------------
///             obj manager
/// Inputs :
///     vector of names
///     vector of objects
///----------------------------------------------
///----------------------------------------------

/////////////////////////////////////////////////////////////
/// central function that does the creation of the XL function
/////////////////////////////////////////////////////////////
LPXLOPER Local_ObjManager_Create_Common(
	LPXLOPER XL_names,
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

		VECTOR<CCString> C_namesCCStr;
	    XL_readStrVector(XL_names,C_namesCCStr," ARM_ERR: names: array of string expected", XL_TYPE_STRING, C_result);

		/// conversion to vector<string>
		vector<string> C_Names( C_namesCCStr.size() );
		for(size_t i=0; i<C_namesCCStr.size() ; ++i )
			C_Names[i] = C_namesCCStr[i];

		VECTOR<CCString> C_Objects;
		VECTOR<CCString> C_ObjectsDef(0);
		XL_readStrVectorWD (XL_values,C_Objects,C_ObjectsDef," ARM_ERR: Array of objects expected",DOUBLE_TYPE,C_result);
		size_t size = C_Objects.size();

		VECTOR<long> CObjectsId;
		CObjectsId.resize(size);
		for(int i = 0; i < size; ++i )    
			CObjectsId[i] = LocalGetNumObjectId(C_Objects[i]);


		/// use the concept of Functor to transfer the knowledge of
		/// a function with a context
		exportFunc2Args< VECTOR<string>, VECTOR<long> >  ourFunc(C_Names,CObjectsId, ARMLOCAL_ObjManager_Create );

		/// call the general function
		fillXL_Result( LOCAL_CSTMANAGER_CLASS, ourFunc, C_result, XL_result, PersistentInXL );
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_ObjManager_Create_Common" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}
							 

///////////////////////////////////
/// version that takes into account 
/// previous creation of object
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_ObjManager_Create(
	LPXLOPER XL_names,
	LPXLOPER XL_values )
{
	ADD_LOG("Local_ObjManager_Create");
	bool PersistentInXL = true;
	return Local_ObjManager_Create_Common(
		XL_names,
		XL_values,
		PersistentInXL );
}


///////////////////////////////////
/// version for non persistentInXl
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_PXL_ObjManager_Create(
	LPXLOPER XL_names,
	LPXLOPER XL_values)
{
	ADD_LOG("Local_PXL_ObjManager_Create");
	bool PersistentInXL = false;
	return Local_ObjManager_Create_Common(
		XL_names,
		XL_values,
		PersistentInXL );
}




/////////////////////////////////////////////////////////////
/// central function that does the creation of the XL function
/////////////////////////////////////////////////////////////
LPXLOPER Local_GenSec_ExtractSubDealDesCommon(
	LPXLOPER XL_genSecId,
	LPXLOPER XL_cutoffColName,
	LPXLOPER XL_columnNames,
	LPXLOPER XL_discounting,
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

		long C_genSecId;
		XL_GETOBJID( XL_genSecId,	C_genSecId, " ARM_ERR: Cst MAnager id: object expected",C_result);

		CCString colNameXl;
		XL_readStrCell(XL_cutoffColName,colNameXl," ARM_ERR: column name: string expected",C_result);
        string C_colName = CCSTringToSTLString(colNameXl);

		VECTOR<CCString> C_columnNamesVecCCStr;
		VECTOR<CCString> C_columnNamesVecCCStrDef;
		XL_readStrVectorWD(XL_columnNames,C_columnNamesVecCCStr,C_columnNamesVecCCStrDef," ARM_ERR: column anmes: array of string expected", DOUBLE_TYPE, C_result);
		vector<string> C_columnNames(C_columnNamesVecCCStr.size());
		for( size_t i=0; i<C_columnNamesVecCCStr.size(); ++i )
			C_columnNames[i] = CCSTringToSTLString(C_columnNamesVecCCStr[i]);

		CCString payCurveNameXl;
		XL_readStrCellWD(XL_discounting, payCurveNameXl,""," ARM_ERR: discount curve name: string expected",C_result);
        string payCurveName = CCSTringToSTLString(payCurveNameXl);

		exportFunc4Args< long, string, vector< string >, string  >  ourFunc( C_genSecId, C_colName, C_columnNames, payCurveName, ARMLOCAL_GenSec_ExtractSubDealDes );

		/// call the general function
		fillXL_Result( LOCAL_GENSEC_CLASS, ourFunc, C_result, XL_result, PersistentInXL );
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_GenSec_ExtractSubDealDesCommon" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}

///////////////////////////////////
/// version that takes into account 
/// previous creation of object
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_DealDes_ExtractSubDealDes(
	LPXLOPER XL_genSecId,
	LPXLOPER XL_cutoffColName,
	LPXLOPER XL_columnNames,
	LPXLOPER XL_discounting )
{
	ADD_LOG("Local_DealDes_ExtractSubDealDes");
	bool PersistentInXL = true;
	return Local_GenSec_ExtractSubDealDesCommon(
		XL_genSecId,
		XL_cutoffColName,
		XL_columnNames,
		XL_discounting,
		PersistentInXL );
}


///////////////////////////////////
/// version for non persistentInXl
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_PXL_DealDes_ExtractSubDealDes(
	LPXLOPER XL_genSecId,
	LPXLOPER XL_cutoffColName,
	LPXLOPER XL_columnNames,
	LPXLOPER XL_discounting )
{
	ADD_LOG("Local_PXL_DealDes_ExtractSubDealDes");
	bool PersistentInXL = false;
	return Local_GenSec_ExtractSubDealDesCommon(
		XL_genSecId,
		XL_cutoffColName,
		XL_columnNames,
		XL_discounting,
		PersistentInXL );
}

