/*
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 * $Log: ARM_xl_gp_calib_local.cpp,v $
 * Revision 1.1  2003/13/07 15:08:43  ebenhamou
 * Initial version
 *
 */


#include <ARM\libarm_local\firstToBeIncluded.h>
#include <functional>
#include <libCCxll\CCxll.h>

#include <ARM\libarm_local\ARM_local_gp_calib.h>
#include "ARM_xl_gp_calib_local.h"
#include "ARM_xl_wrapper_local.h"
#include "ARM_xl_trycatch_local.h"
#include "ARM_gp_local_interglob.h"
#include "ARM_xl_gp_fctorhelper.h"
#include <util\fromto.h>

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
///             Calibrate Model
/// Inputs :
///     Model Id
///     Portfolio Id
///----------------------------------------------
///----------------------------------------------
////////////////////////////////////////////////
/// very rapid definition of ourlonglongFunctor
////////////////////////////////////////////////
class CalibratePFFunc : public ARMResultLong2LongFunc
{
public:
	CalibratePFFunc(long modelId,long CalibMethodId)
                    :C_modelId(modelId), 
                     C_CalibMethodId(CalibMethodId)
    {};    
	
    long operator()( ARM_result& result, long objId ){
        return ARMLOCAL_Calibrate(C_modelId, C_CalibMethodId,result, objId );
    }
			

private:
	long C_modelId;
    long C_CalibMethodId;
};

/////////////////////////////////////////////////////////////
/// central function that does the creation of the XL function
/////////////////////////////////////////////////////////////
LPXLOPER Local_CalibratePFCommon(
	LPXLOPER XL_ModelId,
    LPXLOPER XL_CalibMethodId,
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
		
		/// call the general function
		/// this is used by macros 
		/// and therefore this has to be defined
		static int error;
		static char* reason = "";	

		CCString C_modelStrId;
		XL_readStrCell( XL_ModelId, C_modelStrId,	" ARM_ERR: Model Id: Object expected",		C_result);
		
		long C_ModelId = LocalGetNumObjectId(C_modelStrId);

		CCString C_methodcalibStrId;
		XL_readStrCell(XL_CalibMethodId,C_methodcalibStrId," ARM_ERR: CalibMethod Id: object expected",C_result);
		long C_methodcalibId = LocalGetNumObjectId(C_methodcalibStrId); 

		/// use the concept of Functor to transfer the knowledge of
		/// a function with a context
		CalibratePFFunc ourFunc(C_ModelId,C_methodcalibId);
		CCString curClass = LocalGetStringObjectClassAndError(C_modelStrId);

		/// call the general function

		if( curClass != CCString(" ") )
			fillXL_Result( curClass, ourFunc, C_result, XL_result, PersistentInXL );
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_CalibratePFCommon" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}

				 

///////////////////////////////////
/// version that takes into account 
/// previous creation of object
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_Calibrate(
	LPXLOPER XL_ModelId,
    LPXLOPER XL_CalibMethodId)
{
	ADD_LOG("Local_Calibrate");
	bool PersistentInXL = true;
	return Local_CalibratePFCommon( XL_ModelId,
                                    XL_CalibMethodId,
                                    PersistentInXL );
}

///////////////////////////////////
/// version for non persistentInXl
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_PXL_Calibrate(
	LPXLOPER XL_ModelId,
	LPXLOPER XL_CalibMethodId)
{
	ADD_LOG("Local_PXL_Calibrate");
	bool PersistentInXL = false;
	return Local_CalibratePFCommon(XL_ModelId,
                                   XL_CalibMethodId,
                                   PersistentInXL );
}


///----------------------------------------------
///----------------------------------------------
///             Calib Method
/// Inputs :
///     Portfolio Id
///     Vector of CalibParms Ids
///     Type of calibration
///     CalibMethod (optional = NULL)
///----------------------------------------------
///----------------------------------------------
////////////////////////////////////////////////


/////////////////////////////////////////////////////////////
/// central function that does the creation of the XL function
/////////////////////////////////////////////////////////////
LPXLOPER Local_CalibMethodCommon(
	LPXLOPER XL_PortfolioId,
	LPXLOPER XL_CalibParamsIds,
    LPXLOPER XL_Type,
    LPXLOPER XL_Max_iter,
    LPXLOPER XL_TargetFuncType,
	LPXLOPER XL_LinkedCalibMethod,
    LPXLOPER XL_PreviousCalibMethod,
    LPXLOPER XL_FactorNb,
	LPXLOPER XL_Validate,
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

		CCString C_pfStrId;
		XL_readStrCell( XL_PortfolioId, C_pfStrId,	" ARM_ERR: Portfolio Id: Object expected",		C_result);
		long C_PflId = LocalGetNumObjectId(C_pfStrId);    
    
		VECTOR<CCString> C_CalibParms;
		XL_readStrVector (XL_CalibParamsIds,C_CalibParms," ARM_ERR: CalibParmam: array of object expected",DOUBLE_TYPE,C_result);
    
		size_t size = C_CalibParms.size();
		VECTOR<long> C_CalibParmsId;
		C_CalibParmsId.resize(size);
		size_t i;
		for(i = 0; i < size; ++i )    
			C_CalibParmsId[i] = LocalGetNumObjectId(C_CalibParms[i]); 

		CCString C_CalibNamestr;
		CCString defaultCalibName = "Optimize";
		XL_readStrCellWD(XL_Type,C_CalibNamestr,defaultCalibName," ARM_ERR: CalibType Name: String expected",C_result);

		double C_max_iter;
		double max_iterDefault = 100;
		XL_readNumCellWD(XL_Max_iter, C_max_iter, max_iterDefault, " ARM_ERR: Max_iter : numeric expected",C_result);	  

		CCString C_TargetFuncTypeNamestr;
		CCString defaultTargetFuncTypeName = "PRICE_TAR";
		XL_readStrCellWD(XL_TargetFuncType,C_TargetFuncTypeNamestr,defaultTargetFuncTypeName," ARM_ERR: TargetFuncType Name: String expected",C_result);

		CCString C_LinkedCalibMethodStrId;
		XL_readStrCellWD(XL_LinkedCalibMethod,C_LinkedCalibMethodStrId,"NULL"," ARM_ERR: Linked CalibMethod Name: object expected",C_result);
        long C_LinkedCalibMethodId =  (C_LinkedCalibMethodStrId == "NULL") ? ARM_NULL_OBJECT : LocalGetNumObjectId(C_LinkedCalibMethodStrId);

        CCString C_PreviousCalibMethodStrId;
		XL_readStrCellWD(XL_PreviousCalibMethod,C_PreviousCalibMethodStrId,"NULL"," ARM_ERR: Linked CalibMethod Name: object expected",C_result);
        long C_PreviousCalibMethodId =  (C_PreviousCalibMethodStrId == "NULL") ? ARM_NULL_OBJECT : LocalGetNumObjectId(C_PreviousCalibMethodStrId);

		double C_FactorNb;
		double C_FactorNbDbleDefault=0;
		XL_readNumCellWD(XL_FactorNb,C_FactorNb,C_FactorNbDbleDefault," ARM_ERR: Factor Nb: object expected",C_result);

		double C_validateDble;
		double validateDefault = 1;
		XL_readNumCellWD( XL_Validate, C_validateDble, validateDefault, " ARM_ERR: Validate: boolean expected",	C_result);
		bool C_validate = C_validateDble != 0;

		/// use the concept of Functor to transfer the knowledge of
		/// a function with a context
		exportFunc9Args<long, VECTOR<long>, CCString, double, CCString, long, long, double, bool > 
			ourFunc(C_PflId,C_CalibParmsId,C_CalibNamestr,C_max_iter,C_TargetFuncTypeNamestr,C_LinkedCalibMethodId, 
				C_PreviousCalibMethodId,C_FactorNb,C_validate, ARMLOCAL_CalibMethod_Create);

		//CalibMethodFunc ourFunc(C_PflId,C_CalibParmsId,C_CalibType,C_Max_iter,C_TargetFuncType,C_LinkedCalibMethodId, C_PreviousCalibMethodId,C_FactorNb);

		/// call the general function
		fillXL_Result( LOCAL_CALIBMETHOD_CLASS, ourFunc, C_result, XL_result, PersistentInXL );
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_CalibMethodCommon" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}
							 

///////////////////////////////////
/// version that takes into account 
/// previous creation of object
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_CalibMethod_Create(
    LPXLOPER XL_PortfolioId,
	LPXLOPER XL_CalibParamsIds,
    LPXLOPER XL_Type,
    LPXLOPER XL_Max_iter,
    LPXLOPER XL_TargetFuncType,
	LPXLOPER XL_LinkedCalibMethod,
    LPXLOPER XL_PreviousCalibMethod,
	LPXLOPER XL_FactorNb,
	LPXLOPER XL_Validate)
{
	ADD_LOG("Local_CalibMethod_Create");
	bool PersistentInXL = true;
	return Local_CalibMethodCommon(
         XL_PortfolioId,
	     XL_CalibParamsIds,
         XL_Type,
         XL_Max_iter,
         XL_TargetFuncType,
	     XL_LinkedCalibMethod,
         XL_PreviousCalibMethod,
		 XL_FactorNb,
		 XL_Validate,
         PersistentInXL );
}


///////////////////////////////////
/// version for non persistentInXl
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_PXL_CalibMethod_Create(
    LPXLOPER XL_PortfolioId,
	LPXLOPER XL_CalibParamsIds,
    LPXLOPER XL_Type,
    LPXLOPER XL_Max_iter,
    LPXLOPER XL_TargetFuncType,
	LPXLOPER XL_LinkedCalibMethod,
    LPXLOPER XL_PreviousCalibMethod,
	LPXLOPER XL_FactorNb,
	LPXLOPER XL_Validate)
{
	ADD_LOG("Local_PXL_CalibMethod_Create");
	bool PersistentInXL = false;
	return Local_CalibMethodCommon(
         XL_PortfolioId,
	     XL_CalibParamsIds,
         XL_Type,
         XL_Max_iter,
         XL_TargetFuncType,
	     XL_LinkedCalibMethod,
         XL_PreviousCalibMethod,
		 XL_FactorNb,
		 XL_Validate,
         PersistentInXL );
}

//////////////////////////////////////////////////////////////////////////////
/// central function that does the creation of the XL function with Description
/////////////////////////////////////////////////////////////////////////////
LPXLOPER Local_CalibMethodWithDescriptionCommon(
	LPXLOPER XL_MethodType,
	LPXLOPER XL_PortfolioId,
	LPXLOPER XL_CalibParamsIds,
    LPXLOPER XL_ModelFitterDesId,
    LPXLOPER XL_TargetFuncType,
	LPXLOPER XL_LinkedCalibMethod,
    LPXLOPER XL_PreviousCalibMethod,
	LPXLOPER XL_FactorNb,
	LPXLOPER XL_NbIteration,
	LPXLOPER XL_Validate,
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

		CCString C_MethodType;
		XL_readStrCell( XL_MethodType, C_MethodType, " ARM_ERR: Method Type: string expected",	C_result);

		CCString C_pfStrId;
		XL_readStrCell( XL_PortfolioId, C_pfStrId,	" ARM_ERR: Portfolio Id: Object expected",		C_result);
		long C_PflId = LocalGetNumObjectId(C_pfStrId);    
    
		VECTOR<CCString> C_CalibParms;
		XL_readStrVector (XL_CalibParamsIds,C_CalibParms," ARM_ERR: CalibParmam: array of object expected",DOUBLE_TYPE,C_result);
    
		size_t size = C_CalibParms.size();
		VECTOR<long> C_CalibParmsId;
		C_CalibParmsId.resize(size);
		size_t i;
		for(i = 0; i < size; ++i )    
			C_CalibParmsId[i] = LocalGetNumObjectId(C_CalibParms[i]); 

		CCString C_mfdesstrId;
		XL_readStrCell( XL_ModelFitterDesId, C_mfdesstrId,	" ARM_ERR: ModelFitterDes Id: Object expected",		C_result);
		long C_MFDesId = LocalGetNumObjectId(C_mfdesstrId);

		CCString C_TargetFuncTypeNamestr;
		CCString defaultTargetFuncTypeName = "PRICE_TAR";
		XL_readStrCellWD(XL_TargetFuncType,C_TargetFuncTypeNamestr,defaultTargetFuncTypeName," ARM_ERR: TargetFuncType Name: String expected",C_result);

		CCString C_LinkedCalibMethodStrId;
		XL_readStrCellWD(XL_LinkedCalibMethod,C_LinkedCalibMethodStrId,"NULL"," ARM_ERR: Linked CalibMethod Name: object expected",C_result);
        long C_LinkedCalibMethodId =  (C_LinkedCalibMethodStrId == "NULL") ? ARM_NULL_OBJECT : LocalGetNumObjectId(C_LinkedCalibMethodStrId);

        CCString C_PreviousCalibMethodStrId;
		XL_readStrCellWD(XL_PreviousCalibMethod,C_PreviousCalibMethodStrId,"NULL"," ARM_ERR: Linked CalibMethod Name: object expected",C_result);
        long C_PreviousCalibMethodId =  (C_PreviousCalibMethodStrId == "NULL") ? ARM_NULL_OBJECT : LocalGetNumObjectId(C_PreviousCalibMethodStrId);

		double C_FactorNb;
		double C_FactorNbDbleDefault=0;
		XL_readNumCellWD(XL_FactorNb,C_FactorNb,C_FactorNbDbleDefault," ARM_ERR: Factor Nb: object expected",C_result);

		double C_NbIteration;
		double C_NbIterDbleDefault=1;
		XL_readNumCellWD(XL_NbIteration,C_NbIteration,C_NbIterDbleDefault," ARM_ERR: Nb Iteration: object expected",C_result);

		double C_validateDble;
		double validateDefault = 1;
		XL_readNumCellWD( XL_Validate, C_validateDble, validateDefault, " ARM_ERR: Validate: boolean expected",	C_result);
		bool C_validate = C_validateDble != 0;

		/// use the concept of Functor to transfer the knowledge of
		/// a function with a context
		exportFunc10Args<CCString,long, VECTOR<long>, long, CCString, long, long, double,double,bool > 
			ourFunc(C_MethodType,C_PflId,C_CalibParmsId,C_MFDesId,C_TargetFuncTypeNamestr,C_LinkedCalibMethodId, 
			C_PreviousCalibMethodId,C_FactorNb,C_NbIteration,C_validate,ARMLOCAL_CalibMethodWithDescription_Create);

		/// call the general function
		fillXL_Result( LOCAL_CALIBMETHOD_CLASS, ourFunc, C_result, XL_result, PersistentInXL );
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_CalibMethodWithDescriptionCommon" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}
							 

///////////////////////////////////
/// version that takes into account 
/// previous creation of object
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_CalibMethodWithDesc_Create(
	LPXLOPER XL_CalibMethod,
    LPXLOPER XL_PortfolioId,
	LPXLOPER XL_CalibParamsIds,
    LPXLOPER XL_ModelFitterDesId,
    LPXLOPER XL_TargetFuncType,
	LPXLOPER XL_LinkedCalibMethod,
    LPXLOPER XL_PreviousCalibMethod,
	LPXLOPER XL_FactorNb,
	LPXLOPER XL_NbIteration,
	LPXLOPER XL_Validate )
{
	ADD_LOG("Local_CalibMethodWithDesc_Create");
	bool PersistentInXL = true;
	return Local_CalibMethodWithDescriptionCommon(
		 XL_CalibMethod,
         XL_PortfolioId,
	     XL_CalibParamsIds,
         XL_ModelFitterDesId,
         XL_TargetFuncType,
	     XL_LinkedCalibMethod,
         XL_PreviousCalibMethod,
		 XL_FactorNb,
		 XL_NbIteration,
		 XL_Validate,
         PersistentInXL );
}


///////////////////////////////////
/// version for non persistentInXl
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_PXL_CalibMethodWithDesc_Create(
	LPXLOPER XL_CalibMethod,
    LPXLOPER XL_PortfolioId,
	LPXLOPER XL_CalibParamsIds,
    LPXLOPER XL_ModelFitterDesId,
    LPXLOPER XL_TargetFuncType,
	LPXLOPER XL_LinkedCalibMethod,
    LPXLOPER XL_PreviousCalibMethod,
	LPXLOPER XL_FactorNb,
	LPXLOPER XL_NbIteration,
	LPXLOPER XL_Validate )
{
	ADD_LOG("Local_PXL_CalibMethodWithDesc_Create");
	bool PersistentInXL = false;
	return Local_CalibMethodWithDescriptionCommon(
		 XL_CalibMethod,
         XL_PortfolioId,
	     XL_CalibParamsIds,
         XL_ModelFitterDesId,
         XL_TargetFuncType,
	     XL_LinkedCalibMethod,
         XL_PreviousCalibMethod,
		 XL_FactorNb,
		 XL_NbIteration,
		 XL_Validate,
         PersistentInXL );
}


////////////////////////////////////////////////////////////
/// function to create a Nag Optimizer
/////////////////////////////////////////////////////////////
LPXLOPER Local_OptimizerCommon(
	LPXLOPER XL_algoType,
    LPXLOPER XL_maxIter,
	LPXLOPER XL_tol,
    LPXLOPER XL_stepMax,
	LPXLOPER XL_localSearch,
	LPXLOPER XL_printLevel,
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

		CCString C_algoType;
		XL_readStrCell( XL_algoType, C_algoType, " ARM_ERR: Algorithm Type: string expected",	C_result);
		
		double C_Max_iterDble;
		double maxIterDefault = 100.0;
		XL_readNumCellWD( XL_maxIter, C_Max_iterDble, maxIterDefault, " ARM_ERR: max iterations: numeric expected",	C_result);

		double C_tol;
		double toleranceDefault = 3.26e-12;
		XL_readNumCellWD(XL_tol, C_tol, toleranceDefault, " ARM_ERR: Tolerance: numeric expected",	C_result);

		double C_stepMax;
		double stepMaxDefault = 2.0;
		XL_readNumCellWD( XL_stepMax, C_stepMax, stepMaxDefault, " ARM_ERR: Step Max: numeric expected",	C_result);

		double C_localSearchDble;
		double localSearchDefault = 0;
		XL_readNumCellWD( XL_localSearch, C_localSearchDble, localSearchDefault, " ARM_ERR: Local Search: boolean expected",	C_result);
		bool C_LocalSearch = C_localSearchDble != 0;

		double C_printLevelDble;
		double printLevelDefault = 0;
		XL_readNumCellWD( XL_printLevel, C_printLevelDble, printLevelDefault, " ARM_ERR: Print Level: boolean expected",	C_result);
		bool C_printLevel = C_printLevelDble != 0;

		/// a function with a context
		exportFunc6Args< CCString, double, double, double, bool, bool >
			ourFunc(C_algoType, C_Max_iterDble, C_tol, C_stepMax, C_LocalSearch, C_printLevel, ARMLOCAL_Optimizer_Create );

		/// call the general function
		fillXL_Result( LOCAL_MODELFITTERDES_CLASS, ourFunc, C_result, XL_result, PersistentInXL );
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_OptimizerCommon" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


///////////////////////////////////
////To Create a Nag Optimizer
/// version that takes into account 
/// previous creation of object
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_ARM_Optimizer_Create(
    LPXLOPER XL_algoType,
	LPXLOPER XL_maxIter,
    LPXLOPER XL_tol,
    LPXLOPER XL_stepMax,
	LPXLOPER XL_localSearch,
    LPXLOPER XL_printLevel)
{
	ADD_LOG("Local_ARM_Optimizer_Create");
	bool PersistentInXL = true;
	return Local_OptimizerCommon(
         XL_algoType,
	     XL_maxIter,
         XL_tol,
         XL_stepMax,
	     XL_localSearch,
         XL_printLevel,
         PersistentInXL );
}


///////////////////////////////////
/// version for non persistentInXl
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_PXL_Optimizer_Create(
    LPXLOPER XL_algoType,
	LPXLOPER XL_maxIter,
    LPXLOPER XL_tol,
    LPXLOPER XL_stepMax,
	LPXLOPER XL_localSearch,
    LPXLOPER XL_printLevel)
{
	ADD_LOG("Local_PXL_Optimizer_Create");
	bool PersistentInXL = false;
	return Local_OptimizerCommon(
         XL_algoType,
	     XL_maxIter,
         XL_tol,
         XL_stepMax,
	     XL_localSearch,
         XL_printLevel,
         PersistentInXL );
}

////////////////////////////////////////////////////////////
/// function to create a Nag Optimizer
/////////////////////////////////////////////////////////////
LPXLOPER Local_SolverCommon(
	LPXLOPER XL_algoType,
	LPXLOPER XL_maxIter,
    LPXLOPER XL_xTol,
	LPXLOPER XL_fxTol,
	LPXLOPER XL_gradTol,
    LPXLOPER XL_printLevel,
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

		CCString C_algoType;
		XL_readStrCell( XL_algoType, C_algoType, " ARM_ERR: Solver Type: string expected",	C_result);
		
		double C_Max_iterDble;
		double maxIterDefault = 100.0;
		XL_readNumCellWD( XL_maxIter, C_Max_iterDble, maxIterDefault, " ARM_ERR: max iterations: numeric expected",	C_result);

		double C_xTol;
		double xToleranceDefault = 3.26e-12;
		XL_readNumCellWD(XL_xTol, C_xTol, xToleranceDefault, " ARM_ERR: xTolerance: numeric expected",	C_result);

		double C_fxTol;
		double fxToleranceDefault = 1.0e-14;
		XL_readNumCellWD(XL_fxTol, C_fxTol, fxToleranceDefault, " ARM_ERR: fxTolerance: numeric expected",	C_result);

		double C_gradTol;
		double gradToleranceDefault = 1.0e-14;
		XL_readNumCellWD(XL_gradTol, C_gradTol, gradToleranceDefault, " ARM_ERR: gradTolerance: numeric expected",	C_result);

		double C_printLevelDble;
		double printLevelDefault = 0;
		XL_readNumCellWD( XL_printLevel, C_printLevelDble, printLevelDefault, " ARM_ERR: Print Level: boolean expected",	C_result);
		bool C_printLevel = C_printLevelDble != 0;

		/// a function with a context
		exportFunc6Args< CCString, double, double, double, double, bool > 
			ourFunc(C_algoType, C_Max_iterDble, C_xTol, C_fxTol, C_gradTol, C_printLevel, ARMLOCAL_Solver_Create );

		/// call the general function
		fillXL_Result( LOCAL_MODELFITTERDES_CLASS, ourFunc, C_result, XL_result, PersistentInXL );
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_SolverCommon" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


///////////////////////////////////
////To Create a Solver
/// version that takes into account 
/// previous creation of object
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_ARM_Solver_Create(
    LPXLOPER XL_algoType,
	LPXLOPER XL_maxIter,
    LPXLOPER XL_xTol,
	LPXLOPER XL_fxTol,
	LPXLOPER XL_gradTol,
    LPXLOPER XL_printLevel)
{
	ADD_LOG("Local_ARM_Solver_Create");
	bool PersistentInXL = false;
	return Local_SolverCommon(
         XL_algoType,
	     XL_maxIter,
         XL_xTol,
         XL_fxTol,
	     XL_gradTol,
         XL_printLevel,
         PersistentInXL );
}

///////////////////////////////////
/// version for non persistentInXl
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_PXL_ARM_Solver_Create(
    LPXLOPER XL_algoType,
	LPXLOPER XL_maxIter,
    LPXLOPER XL_xTol,
	LPXLOPER XL_fxTol,
	LPXLOPER XL_gradTol,
    LPXLOPER XL_printLevel)
{
	ADD_LOG("Local_PXL_ARM_Solver_Create");
	bool PersistentInXL = false;
	return Local_SolverCommon(
         XL_algoType,
	     XL_maxIter,
         XL_xTol,
         XL_fxTol,
	     XL_gradTol,
         XL_printLevel,
         PersistentInXL );
}

//////////////////////////////////////////////////////////////////////////////
/// central function that does the creation of the XL function with Description
/////////////////////////////////////////////////////////////////////////////
LPXLOPER Local_CalibMethod2DCommon(
	LPXLOPER XL_Portfolio1Id,
	LPXLOPER XL_Portfolio2Id,
	LPXLOPER XL_CalibParams1Ids,
	LPXLOPER XL_CalibParams2Ids,
    LPXLOPER XL_ModelFitterDes1Id,
	LPXLOPER XL_ModelFitterDes2Id,
    LPXLOPER XL_TargetFuncType,
	LPXLOPER XL_LinkedCalibMethod,
    LPXLOPER XL_PreviousCalibMethod,
	LPXLOPER XL_Calib2DDirection,
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

		CCString C_pf1StrId;
		XL_readStrCell( XL_Portfolio1Id, C_pf1StrId,	" ARM_ERR: Portfolio Id: Object expected",		C_result);
		long C_PflId = LocalGetNumObjectId(C_pf1StrId);   
		
		
		CCString C_pf2StrId;
		XL_readStrCell( XL_Portfolio2Id, C_pf2StrId,	" ARM_ERR: Portfolio Id: Object expected",		C_result);
		long C_Pf2Id = LocalGetNumObjectId(C_pf2StrId);   
    
		VECTOR<CCString> C_CalibParms1;
		XL_readStrVector (XL_CalibParams1Ids,C_CalibParms1," ARM_ERR: CalibParmam: array of object expected",DOUBLE_TYPE,C_result);

		VECTOR<CCString> C_CalibParms2;
		XL_readStrVector (XL_CalibParams2Ids,C_CalibParms2," ARM_ERR: CalibParmam: array of object expected",DOUBLE_TYPE,C_result);
    
		size_t size1 = C_CalibParms1.size();
		VECTOR<long> C_CalibParms1Id;
		C_CalibParms1Id.resize(size1);
		size_t i;
		for(i = 0; i < size1; ++i )    
			C_CalibParms1Id[i] = LocalGetNumObjectId(C_CalibParms1[i]); 

		size_t size2 = C_CalibParms2.size();
		if(size2!=size1)
		{
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}

		VECTOR<long> C_CalibParms2Id;
		C_CalibParms2Id.resize(size2);
		for(i = 0; i < size2; ++i )    
			C_CalibParms2Id[i] = LocalGetNumObjectId(C_CalibParms2[i]); 


		CCString C_mfdes1strId;
		XL_readStrCell( XL_ModelFitterDes1Id, C_mfdes1strId,	" ARM_ERR: ModelFitterDes Id: Object expected",		C_result);
		long C_MFDes1Id = LocalGetNumObjectId(C_mfdes1strId);

		CCString C_mfdes2strId;
		XL_readStrCell( XL_ModelFitterDes2Id, C_mfdes2strId,	" ARM_ERR: ModelFitterDes Id: Object expected",		C_result);
		long C_MFDes2Id = LocalGetNumObjectId(C_mfdes2strId);

		CCString C_TargetFuncTypeNamestr;
		CCString defaultTargetFuncTypeName = "PRICE_TAR";
		XL_readStrCellWD(XL_TargetFuncType,C_TargetFuncTypeNamestr,defaultTargetFuncTypeName," ARM_ERR: TargetFuncType Name: String expected",C_result);

		long C_TargetFuncType;
		if( (C_TargetFuncType = ARM_ConvGPTargetFuncType( C_TargetFuncTypeNamestr, C_result)) == ARM_DEFAULT_ERR )
		{
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}

		CCString C_LinkedCalibMethodStrId;
		XL_readStrCellWD(XL_LinkedCalibMethod,C_LinkedCalibMethodStrId,"NULL"," ARM_ERR: Linked CalibMethod Name: object expected",C_result);
        long C_LinkedCalibMethodId =  (C_LinkedCalibMethodStrId == "NULL") ? ARM_NULL_OBJECT : LocalGetNumObjectId(C_LinkedCalibMethodStrId);

        CCString C_PreviousCalibMethodStrId;
		XL_readStrCellWD(XL_PreviousCalibMethod,C_PreviousCalibMethodStrId,"NULL"," ARM_ERR: Linked CalibMethod Name: object expected",C_result);
        long C_PreviousCalibMethodId =  (C_PreviousCalibMethodStrId == "NULL") ? ARM_NULL_OBJECT : LocalGetNumObjectId(C_PreviousCalibMethodStrId);


		CCString C_Calib2DDirectionstr;
		CCString defaultCalib2DDirection = "Forward";
		XL_readStrCellWD(XL_Calib2DDirection,C_Calib2DDirectionstr,defaultCalib2DDirection," ARM_ERR: Calib 2D Direction: String expected",C_result);

		long C_Calib2DDirection;
		if( (C_Calib2DDirection = ARM_ConvGPCalib2DDirection( C_Calib2DDirectionstr, C_result)) == ARM_DEFAULT_ERR )
		{
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}
		/// use the concept of Functor to transfer the knowledge of
		/// a function with a context
		exportFunc10Args<long, long,VECTOR<long>, VECTOR<long>, long, long, long, long, long, long> ourFunc(C_PflId,C_Pf2Id,C_CalibParms1Id,C_CalibParms2Id,C_MFDes1Id,C_MFDes2Id,C_TargetFuncType,C_LinkedCalibMethodId, C_PreviousCalibMethodId,C_Calib2DDirection,ARMLOCAL_CalibMethod2D_Create);
		/// call the general function
		fillXL_Result( LOCAL_CALIBMETHOD_CLASS, ourFunc, C_result, XL_result, PersistentInXL );
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_CalibMethodWithDescriptionCommon" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}
							 

///////////////////////////////////
/// version that takes into account 
/// previous creation of object
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_CalibMethod2D_Create(
    LPXLOPER XL_Portfolio1Id,
	LPXLOPER XL_Portfolio2Id,
	LPXLOPER XL_CalibParams1Ids,
	LPXLOPER XL_CalibParams2Ids,
    LPXLOPER XL_ModelFitterDes1Id,
	LPXLOPER XL_ModelFitterDes2Id,
    LPXLOPER XL_TargetFuncType,
	LPXLOPER XL_LinkedCalibMethod,
    LPXLOPER XL_PreviousCalibMethod,
	LPXLOPER XL_Calib2DDirection)
{
	ADD_LOG("Local_CalibMethod2D_Create");
	bool PersistentInXL = true;
	return Local_CalibMethod2DCommon(
         XL_Portfolio1Id,
		 XL_Portfolio2Id,
	     XL_CalibParams1Ids,
		 XL_CalibParams2Ids,
         XL_ModelFitterDes1Id,
		 XL_ModelFitterDes2Id,
         XL_TargetFuncType,
	     XL_LinkedCalibMethod,
         XL_PreviousCalibMethod,
		 XL_Calib2DDirection,
         PersistentInXL );
}


///////////////////////////////////
/// version for non persistentInXl
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_PXL_CalibMethod2D_Create(
    LPXLOPER XL_Portfolio1Id,
	LPXLOPER XL_Portfolio2Id,
	LPXLOPER XL_CalibParams1Ids,
	LPXLOPER XL_CalibParams2Ids,
    LPXLOPER XL_ModelFitterDes1Id,
	LPXLOPER XL_ModelFitterDes2Id,
    LPXLOPER XL_TargetFuncType,
	LPXLOPER XL_LinkedCalibMethod,
    LPXLOPER XL_PreviousCalibMethod,
	LPXLOPER XL_Calib2DDirection)
{
	ADD_LOG("Local_PXL_CalibMethod2D_Create");
	bool PersistentInXL = false;
	return Local_CalibMethod2DCommon(
         XL_Portfolio1Id,
		 XL_Portfolio2Id,
	     XL_CalibParams1Ids,
		 XL_CalibParams2Ids,
         XL_ModelFitterDes1Id,
		 XL_ModelFitterDes2Id,
         XL_TargetFuncType,
	     XL_LinkedCalibMethod,
         XL_PreviousCalibMethod,
		 XL_Calib2DDirection,
         PersistentInXL );
}

///----------------------------------------------
///----------------------------------------------
///             SetDetailFlagToCalibMethod 
///----------------------------------------------
///----------------------------------------------

/////////////////////////////////////////////////////////////
/// central function that does the creation of the XL function
/////////////////////////////////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_SetDetailFlagToCalibMethod(
	LPXLOPER XL_CalibMethodId,
	LPXLOPER XL_DetailFlag )
{
	ADD_LOG("Local_SetDetailFlagToCalibMethod");
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

		double C_DetailFlag;
		long C_CalibMethodId;
		XL_GETOBJID( XL_CalibMethodId,	C_CalibMethodId,	" ARM_ERR: Calib Method id: object expected",			C_result);
		XL_readNumCell( XL_DetailFlag, C_DetailFlag,	    " ARM_ERR: Detail Flag: boolean compatible expected",	C_result);
		bool detailFlag = C_DetailFlag != 0;
			
		long retCode;
		if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 1)
			retCode = ARMLOCAL_SetDetailFlagToCalibMethod(
				C_CalibMethodId,
				detailFlag,
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_SetDetailFlagToCalibMethod" )

	return (LPXLOPER)&XL_result;
}



///----------------------------------------------
///----------------------------------------------
///             Local_GetDurationFromCalibMethod
///----------------------------------------------
///----------------------------------------------

/////////////////////////////////////////////////////////////
/// central function that does the creation of the XL function
/////////////////////////////////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_GetDurationFromCalibMethod(
	LPXLOPER XL_CalibMethodId )
{
	ADD_LOG("Local_GetDurationFromCalibMethod");
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

		long C_CalibMethodId;
		XL_GETOBJID( XL_CalibMethodId,	C_CalibMethodId,	" ARM_ERR: Calib Method id: object expected",			C_result);
		
		long retCode;
		if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 1)
			retCode = ARMLOCAL_GetDurationFromCalibMethod(
				C_CalibMethodId,
				C_result );
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_GetDurationFromCalibMethod" )

	return (LPXLOPER)&XL_result;
}
							 
///----------------------------------------------
///             Get from Object Id
/// Inputs :
///     CalibMethid Id
///     string to localise data object
///----------------------------------------------
///----------------------------------------------
////////////////////////////////////////////////
/// very rapid definition of ourlonglongFunctor
////////////////////////////////////////////////
class GetDataFromCalibMethodFunc : public ARMResultLong2LongFunc
{
public:
	GetDataFromCalibMethodFunc(long calibMethodId,const CCString& dataType)
                    :C_calibMethodId(calibMethodId), 
                     C_dataType(dataType)
    {};    
	
    long operator()( ARM_result& result, long objId ){
        return ARMLOCAL_DataFromCalibMethod(C_calibMethodId, C_dataType,result, objId);
    }
			

private:
	long C_calibMethodId;
    const CCString& C_dataType;
};



/////////////////////////////////////////////////////////////
/// central function that does the creation of the XL function
/////////////////////////////////////////////////////////////
LPXLOPER Local_GetDataFromCalibMethodCommon(
	LPXLOPER XL_CalibMethodId,
    LPXLOPER XL_DataType,
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

    CCString C_CalibMethodStrId;
	XL_readStrCell( XL_CalibMethodId, C_CalibMethodStrId,	" ARM_ERR: CalibMethod Id: Object expected",C_result);
	long C_CalibMethodId = LocalGetNumObjectId(C_CalibMethodStrId);


    CCString C_DataType;
    XL_readStrCell(XL_DataType, C_DataType, " ARM_ERR: DataType : string expected",C_result);	  

	/// use the concept of Functor to transfer the knowledge of
	/// a function with a context
	GetDataFromCalibMethodFunc ourFunc(C_CalibMethodId,C_DataType);

    /// given the name class appeared in excel
    CCString DataGetClass(DataCalibMethodGetClass(C_DataType).c_str());

	/// call the general function
    fillXL_Result(DataGetClass, ourFunc, C_result, XL_result, PersistentInXL );

	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_GetDataFromCalibMethodCommon" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}	

///////////////////////////////////
/// version that takes into account 
/// previous creation of object
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_GetDataFromCalibMethod(
	LPXLOPER XL_CalibMethodId,
    LPXLOPER XL_DataType)
{
	ADD_LOG("Local_GetDataFromCalibMethod");
	bool PersistentInXL = true;
	return Local_GetDataFromCalibMethodCommon( XL_CalibMethodId,
                                    XL_DataType,
                                    PersistentInXL );
}

///////////////////////////////////
/// version for non persistentInXl
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_PXL_GetDataFromCalibMethod(
	LPXLOPER XL_CalibMethodId,
    LPXLOPER XL_DataType)
{
	ADD_LOG("Local_PXL_GetDataFromCalibMethod");
	bool PersistentInXL = false;
	return Local_GetDataFromCalibMethodCommon(XL_CalibMethodId,
                                    XL_DataType,
                                    PersistentInXL );
}


/////////////////////////////////////////////////////////////////////////////::
/////////////////////////////////////////////////////////////////////////////::
/////////   Hunt Kennedy Part
/////////////////////////////////////////////////////////////////////////////::
/////////////////////////////////////////////////////////////////////////////::

/////////////////////////////////////////////////////////////
/// central function that does the creation of the XL function
/////////////////////////////////////////////////////////////


LPXLOPER Local_NumericalCalibMethod_CreateCommon(
	LPXLOPER XL_CalibDateStripId,
    LPXLOPER XL_VanillaSecDensitiesId,
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

		CCString C_CalibDateStripStrId;
		XL_readStrCell( XL_CalibDateStripId, C_CalibDateStripStrId,	" ARM_ERR: Portfolio Id: Object expected",		C_result);
		long C_CalibDateStripId = LocalGetNumObjectId(C_CalibDateStripStrId);    
    
		VECTOR<CCString> C_VanillaSecDensitiesId;
		XL_readStrVector (XL_VanillaSecDensitiesId,C_VanillaSecDensitiesId," ARM_ERR: Securities: array of object expected",DOUBLE_TYPE,C_result);

		size_t size = C_VanillaSecDensitiesId.size();
		VECTOR<long> C_VanillaSecDensitiesIdLong (size);
		size_t i;
		for(i = 0; i < size; ++i )
		{
			C_VanillaSecDensitiesIdLong[i] = LocalGetNumObjectId(C_VanillaSecDensitiesId[i]); 
		}

		/// portfolio : optional
		CCString C_PortfolioStrId;
		XL_readStrCellWD(XL_PortfolioId, C_PortfolioStrId, GETDEFAULTVALUESTR," ARM_ERR: porfolio: string expected", C_result);
		long C_PortfolioId ;
		if (C_PortfolioStrId == GETDEFAULTVALUESTR)
			C_PortfolioId = ARM_NULL_OBJECT;    
		else
			C_PortfolioId = LocalGetNumObjectId(C_PortfolioStrId);
		

		/// use the concept of Functor to transfer the knowledge of
		/// a function with a context
		exportFunc3Args< long, VECTOR<long>, long > 
			ourFunc(C_CalibDateStripId,C_VanillaSecDensitiesIdLong,C_PortfolioId,ARMLOCAL_NumericalCalibMethod_Create);

		/// call the general function
		fillXL_Result( LOCAL_CALIBMETHOD_CLASS, ourFunc, C_result, XL_result, PersistentInXL );
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_NumericalCalibMethodCommon" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}



///////////////////////////////////
/// Create a Numerical Calib Method
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_NumericalCalibMethod_Create(
    LPXLOPER XL_CalibDateStripId,
    LPXLOPER XL_VanillaSecDensitiesId,
	LPXLOPER XL_PortfolioId)
{
	ADD_LOG("Local_NumericalCalibMethod_Create");
	bool PersistentInXL = true;
	return Local_NumericalCalibMethod_CreateCommon(
                                    XL_CalibDateStripId,
									XL_VanillaSecDensitiesId, 
									XL_PortfolioId,
									PersistentInXL );
}

__declspec(dllexport) LPXLOPER WINAPI Local_PXL_NumericalCalibMethod_Create(
	LPXLOPER XL_CalibDateStripId,
	LPXLOPER XL_VanillaSecDensitiesId,
	LPXLOPER XL_PortfolioId)
{
	ADD_LOG("Local_PXL_NumericalCalibMethod_Create");
	bool PersistentInXL = false;
	return Local_NumericalCalibMethod_CreateCommon(
                                    XL_CalibDateStripId,
									XL_VanillaSecDensitiesId, 
									XL_PortfolioId,
									PersistentInXL );
}

////////////////////////////////////////////////////////
/// Create Shifted Lognormal Density Functors (HK)
////////////////////////////////////////////////////////

LPXLOPER Local_SLNDensityFunctorCommon(
    LPXLOPER XL_Volatility,
    LPXLOPER XL_Shift,
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

	    double C_Volatility;
	    XL_readNumCell(XL_Volatility,C_Volatility,		" ARM_ERR: volatility numeric expected",C_result);	

	    double C_Shift;
	    XL_readNumCell(XL_Shift,C_Shift,		" ARM_ERR: shift numeric expected",C_result);	

		exportFunc2Args< double, double > 
			ourFunc(C_Volatility,C_Shift,ARMLOCAL_SLNDensityFunctor_Create);

		/// call the general function
		fillXL_Result_withName(ourFunc, C_result, XL_result, PersistentInXL );

	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_NumericalCalibMethodCommon" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}

__declspec(dllexport) LPXLOPER WINAPI Local_PXL_SLNDensityFunctor(
    LPXLOPER XL_Volatility,
    LPXLOPER XL_Shift )
{
	ADD_LOG("Local_PXL_SLNDensityFunctor");
	bool PersistentInXL = false;
	return Local_SLNDensityFunctorCommon(XL_Volatility,XL_Shift,PersistentInXL);
}

__declspec(dllexport) LPXLOPER WINAPI Local_SLNDensityFunctor(
    LPXLOPER XL_Volatility,
    LPXLOPER XL_Shift )
{
	ADD_LOG("Local_SLNDensityFunctor");
	bool PersistentInXL = true;
	return Local_SLNDensityFunctorCommon(XL_Volatility,XL_Shift,PersistentInXL);
}



////////////////////////////////////////////////////////
/// Create Mixture Density Functors (HK)
////////////////////////////////////////////////////////
		
LPXLOPER Local_MixtureDensityFunctorCommon(
    LPXLOPER XL_Volatility1,
	LPXLOPER XL_Volatility2,
	LPXLOPER XL_Alpha,
	LPXLOPER XL_Lambda,
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

	    double C_Volatility1;
	    XL_readNumCell(XL_Volatility1,C_Volatility1,		" ARM_ERR: volatility 1 numeric expected",C_result);	

		double C_Volatility2;
	    XL_readNumCell(XL_Volatility2,C_Volatility2,		" ARM_ERR: volatility 2 numeric expected",C_result);	

		double C_Alpha;
	    XL_readNumCell(XL_Alpha,C_Alpha,		" ARM_ERR: alpha numeric expected",C_result);	

	    double C_Lambda;
	    XL_readNumCell(XL_Lambda,C_Lambda,		" ARM_ERR: lambda numeric expected",C_result);	

		exportFunc4Args< double, double, double, double > 
			ourFunc(C_Volatility1,C_Volatility2,C_Alpha,C_Lambda,ARMLOCAL_MixtureDensityFunctor_Create);

		/// call the general function
		fillXL_Result_withName(ourFunc, C_result, XL_result, PersistentInXL );

	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_NumericalCalibMethodCommon" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}

__declspec(dllexport) LPXLOPER WINAPI Local_PXL_MixtureDensityFunctor(
    LPXLOPER XL_Volatility1,
	LPXLOPER XL_Volatility2,
	LPXLOPER XL_Alpha,
	LPXLOPER XL_Lambda,
    LPXLOPER XL_Shift )
{
	ADD_LOG("Local_PXL_MixtureDensityFunctor");
	bool PersistentInXL = false;
	return Local_MixtureDensityFunctorCommon(XL_Volatility1,XL_Volatility2,XL_Alpha,XL_Lambda,PersistentInXL);
}

									  
__declspec(dllexport) LPXLOPER WINAPI Local_MixtureDensityFunctor(
    LPXLOPER XL_Volatility1,
	LPXLOPER XL_Volatility2,
	LPXLOPER XL_Alpha,
	LPXLOPER XL_Lambda)
{
	ADD_LOG("Local_MixtureDensityFunctor");
	bool PersistentInXL = true;
	return Local_MixtureDensityFunctorCommon(XL_Volatility1,XL_Volatility2,XL_Alpha,XL_Lambda,PersistentInXL);
}


////////////////////////////////////////////////////////
/// Create Mixture Density Functors (HK)
////////////////////////////////////////////////////////
		
LPXLOPER Local_MixtureDensityFunctorWithATMVolCommon(
	LPXLOPER XL_Fwd,
	LPXLOPER XL_Maturity,
    LPXLOPER XL_VolATM,
	LPXLOPER XL_DecVol,
	LPXLOPER XL_Alpha,
	LPXLOPER XL_Lambda,
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

		double C_Fwd;
	    XL_readNumCell(XL_Fwd,C_Fwd,		" ARM_ERR: Fwd numeric expected",C_result);

		double C_Maturity;
	    XL_readNumCell(XL_Maturity,C_Maturity,	" ARM_ERR: Maturity numeric expected",C_result);

	    double C_VolATM;
	    XL_readNumCell(XL_VolATM,C_VolATM,		" ARM_ERR: volatility ATM numeric expected",C_result);	

		double C_DecVol;
	    XL_readNumCell(XL_DecVol,C_DecVol,		" ARM_ERR: dev vol numeric expected",C_result);	

		double C_Alpha;
	    XL_readNumCell(XL_Alpha,C_Alpha,		" ARM_ERR: alpha numeric expected",C_result);	

	    double C_Lambda;
	    XL_readNumCell(XL_Lambda,C_Lambda,		" ARM_ERR: lambda numeric expected",C_result);	

		exportFunc6Args< double, double, double, double, double, double > 
			ourFunc(C_Fwd,C_Maturity,C_VolATM,C_DecVol,C_Alpha,C_Lambda,ARMLOCAL_MixtureDensityFunctor_CreateWithATMVol);

		/// call the general function
		fillXL_Result_withName(ourFunc, C_result, XL_result, PersistentInXL );

	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_NumericalCalibMethodCommon" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}

__declspec(dllexport) LPXLOPER WINAPI Local_PXL_MixtureDensityFunctorWithATMVol(
    LPXLOPER XL_Fwd,
	LPXLOPER XL_Maturity,
    LPXLOPER XL_VolATM,
	LPXLOPER XL_DecVol,
	LPXLOPER XL_Alpha,
	LPXLOPER XL_Lambda)
{
	ADD_LOG("Local_PXL_MixtureDensityFunctorWithATMVol");
	bool PersistentInXL = false;
	return Local_MixtureDensityFunctorWithATMVolCommon(XL_Fwd,XL_Maturity,XL_VolATM,XL_DecVol,XL_Alpha,XL_Lambda,PersistentInXL);
}

									  
__declspec(dllexport) LPXLOPER WINAPI Local_MixtureDensityFunctorWithATMVol(
    LPXLOPER XL_Fwd,
	LPXLOPER XL_Maturity,
    LPXLOPER XL_VolATM,
	LPXLOPER XL_DecVol,
	LPXLOPER XL_Alpha,
	LPXLOPER XL_Lambda)
{
	ADD_LOG("Local_MixtureDensityFunctorWithATMVol");
	bool PersistentInXL = true;
	return Local_MixtureDensityFunctorWithATMVolCommon(XL_Fwd,XL_Maturity,XL_VolATM,XL_DecVol,XL_Alpha,XL_Lambda,PersistentInXL);
}


////////////////////////////////////////////////////////
/// Create Heston Density Functors (HK)
////////////////////////////////////////////////////////
		
LPXLOPER Local_HestonDensityFunctorCommon(
	LPXLOPER XL_V0,
	LPXLOPER XL_Kappa,
	LPXLOPER XL_Theta,
	LPXLOPER XL_VVol,
	LPXLOPER XL_Rho,
	LPXLOPER XL_Shift,
	LPXLOPER XL_Level,
	LPXLOPER XL_Sigma,
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

	    double C_V0;
	    XL_readNumCell(XL_V0,C_V0,				" ARM_ERR: V0 numeric expected",C_result);	

		double C_Kappa;
	    XL_readNumCell(XL_Kappa,C_Kappa,		" ARM_ERR: kappa numeric expected",C_result);	

		double C_Theta;
	    XL_readNumCell(XL_Theta,C_Theta,		" ARM_ERR: theta numeric expected",C_result);	

	    double C_VVol;
	    XL_readNumCell(XL_VVol,C_VVol,			" ARM_ERR: vvol numeric expected",C_result);	

		double C_Rho;
	    XL_readNumCell(XL_VVol,C_Rho,			" ARM_ERR: rho numeric expected",C_result);	

		double C_Shift;
	    XL_readNumCell(XL_Shift,C_Shift,		" ARM_ERR: shift numeric expected",C_result);

		double C_Level;
	    XL_readNumCell(XL_Level,C_Level,		" ARM_ERR: Level numeric expected",C_result);	

		double C_Sigma, C_SigmaDef = 0.0;
	    XL_readNumCellWD(XL_Sigma,C_Sigma, C_SigmaDef, " ARM_ERR: sigma numeric expected",C_result);	

		exportFunc8Args< double, double, double, double, double, double, double, double > 
			ourFunc(C_V0,C_Kappa,C_Theta,C_VVol, C_Rho, C_Shift, C_Level, C_Sigma, ARMLOCAL_HestonDensityFunctor_Create);

		/// call the general function
		fillXL_Result_withName(ourFunc, C_result, XL_result, PersistentInXL );

	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_NumericalCalibMethodCommon" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}

__declspec(dllexport) LPXLOPER WINAPI Local_PXL_HestonDensityFunctor(
	LPXLOPER XL_V0,
	LPXLOPER XL_Kappa,
	LPXLOPER XL_Theta,
	LPXLOPER XL_VVol,
	LPXLOPER XL_Rho,
	LPXLOPER XL_Shift,
	LPXLOPER XL_Level,
	LPXLOPER XL_Sigma)
{
	ADD_LOG("Local_PXL_HestonDensityFunctor");
	bool PersistentInXL = false;
	return Local_HestonDensityFunctorCommon(XL_V0,XL_Kappa,XL_Theta,XL_VVol,XL_Rho,XL_Shift,XL_Level,XL_Sigma,PersistentInXL);
}

									  
__declspec(dllexport) LPXLOPER WINAPI Local_HestonDensityFunctor(
	LPXLOPER XL_V0,
	LPXLOPER XL_Kappa,
	LPXLOPER XL_Theta,
	LPXLOPER XL_VVol,
	LPXLOPER XL_Rho,
	LPXLOPER XL_Shift,
	LPXLOPER XL_Level,
	LPXLOPER XL_Sigma)
{
	ADD_LOG("Local_HestonDensityFunctor");
	bool PersistentInXL = true;
	return Local_HestonDensityFunctorCommon(XL_V0,XL_Kappa,XL_Theta,XL_VVol,XL_Rho,XL_Shift,XL_Level,XL_Sigma,PersistentInXL);
}


////////////////////////////////////////////////////////
/// Create No vol Density Functor (HK)
////////////////////////////////////////////////////////

LPXLOPER Local_IrFwdDensityFunctorCommon(bool PersistentInXL )
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

		exportFunc0Arg ourFunc(ARMLOCAL_IrFwdDensityFunctor_Create);

		/// call the general function
		fillXL_Result_withName(ourFunc, C_result, XL_result, PersistentInXL );

	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_NumericalCalibMethodCommon" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}

__declspec(dllexport) LPXLOPER WINAPI Local_PXL_IrFwdDensityFunctor()
{
	ADD_LOG("Local_PXL_IrFwdDensityFunctor");
	bool PersistentInXL = false;
	return Local_IrFwdDensityFunctorCommon(PersistentInXL);
}

__declspec(dllexport) LPXLOPER WINAPI Local_IrFwdDensityFunctor()
{
	ADD_LOG("Local_IrFwdDensityFunctor");
	bool PersistentInXL = true;
	return Local_IrFwdDensityFunctorCommon(PersistentInXL);
}

 

////////////////////////////////////////////////////////
/// Create SABR Density Functor (HK)
////////////////////////////////////////////////////////

/// Sabr
LPXLOPER Local_SABRDensityFunctorCommon(
    LPXLOPER XL_Alpha,
	LPXLOPER XL_Beta,
	LPXLOPER XL_Rho,
	LPXLOPER XL_Nu,
	LPXLOPER XL_SabrType,
	LPXLOPER XL_GridSize,
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

	    double C_Alpha;
	    XL_readNumCell(XL_Alpha,C_Alpha,	" ARM_ERR: numeric expected",C_result);	

	    double C_Beta;
	    XL_readNumCell(XL_Beta,C_Beta,		" ARM_ERR: numeric expected",C_result);	

		double C_Rho;
	    XL_readNumCell(XL_Rho,C_Rho,		" ARM_ERR: numeric expected",C_result);	

		double C_Nu;
	    XL_readNumCell(XL_Nu,C_Nu,			" ARM_ERR: numeric expected",C_result);	

		long C_SabrType;
		CCString DefaultSabrType = "SABR_IMPLNVOL";
		XL_GETCONVSABRFLAGWD(XL_SabrType,C_SabrType, DefaultSabrType, " ARM_ERR: SabrType flag string expected",C_result);

		double C_GridSize;
		double C_DefaultGridSize = 251;
	    XL_readNumCellWD(XL_GridSize, C_GridSize, C_DefaultGridSize, " ARM_ERR: numeric expected",C_result);	

		exportFunc6Args< double, double, double, double, long, double > 
			ourFunc(C_Alpha,C_Beta, C_Rho, C_Nu, C_SabrType, C_GridSize, ARMLOCAL_SABRDensityFunctor_Create);

		/// call the general function
		fillXL_Result_withName(ourFunc, C_result, XL_result, PersistentInXL );

	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_NumericalCalibMethodCommon" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;

}

__declspec(dllexport) LPXLOPER WINAPI Local_PXL_SABRDensityFunctor(
    LPXLOPER XL_Alpha,
	LPXLOPER XL_Beta,
	LPXLOPER XL_Rho,
	LPXLOPER XL_Nu,
	LPXLOPER XL_SabrType,
	LPXLOPER XL_GridSize)
{
	ADD_LOG("Local_PXL_SABRDensityFunctor");
	bool PersistentInXL = false;
	return Local_SABRDensityFunctorCommon(XL_Alpha, XL_Beta, XL_Rho, XL_Nu, XL_SabrType, XL_GridSize, PersistentInXL);
}

__declspec(dllexport) LPXLOPER WINAPI Local_SABRDensityFunctor(
    LPXLOPER XL_Alpha,
	LPXLOPER XL_Beta,
	LPXLOPER XL_Rho,
	LPXLOPER XL_Nu,
	LPXLOPER XL_SabrType,
	LPXLOPER XL_GridSize)
{
	ADD_LOG("Local_SABRDensityFunctor");
	bool PersistentInXL = true;
	return Local_SABRDensityFunctorCommon(XL_Alpha, XL_Beta, XL_Rho, XL_Nu, XL_SabrType, XL_GridSize, PersistentInXL);
}

// Bi Sabr density functor
LPXLOPER Local_BiSABRDensityFunctorCommon(
    LPXLOPER XL_Alpha1,
	LPXLOPER XL_Beta1,
	LPXLOPER XL_Rho1,
	LPXLOPER XL_Nu1,
    LPXLOPER XL_Alpha2,
	LPXLOPER XL_Beta2,
	LPXLOPER XL_Rho2,
	LPXLOPER XL_Nu2,
	LPXLOPER XL_RhoS1S2,
	LPXLOPER XL_RhoS1V2,
	LPXLOPER XL_RhoS2V1,
	LPXLOPER XL_RhoV1V2,
	LPXLOPER XL_SabrType,
	LPXLOPER XL_GridSize,
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

	    double C_Alpha1;
	    XL_readNumCell(XL_Alpha1,C_Alpha1,	" ARM_ERR: numeric expected",C_result);	

	    double C_Beta1;
	    XL_readNumCell(XL_Beta1,C_Beta1,		" ARM_ERR: numeric expected",C_result);	

		double C_Rho1;
	    XL_readNumCell(XL_Rho1,C_Rho1,		" ARM_ERR: numeric expected",C_result);	

		double C_Nu1;
	    XL_readNumCell(XL_Nu1,C_Nu1,			" ARM_ERR: numeric expected",C_result);	

	    double C_Alpha2;
	    XL_readNumCell(XL_Alpha2,C_Alpha2,	" ARM_ERR: numeric expected",C_result);	

	    double C_Beta2;
	    XL_readNumCell(XL_Beta2,C_Beta2,		" ARM_ERR: numeric expected",C_result);	

		double C_Rho2;
	    XL_readNumCell(XL_Rho2,C_Rho2,		" ARM_ERR: numeric expected",C_result);	

		double C_Nu2;
	    XL_readNumCell(XL_Nu2,C_Nu2,			" ARM_ERR: numeric expected",C_result);	

		double C_RhoS1S2;
	    XL_readNumCell(XL_RhoS1S2,C_RhoS1S2,			" ARM_ERR: numeric expected",C_result);	

		double C_RhoS1V2;
	    XL_readNumCell(XL_RhoS1V2,C_RhoS1V2,			" ARM_ERR: numeric expected",C_result);	

		double C_RhoS2V1;
	    XL_readNumCell(XL_RhoS2V1,C_RhoS2V1,			" ARM_ERR: numeric expected",C_result);	

		double C_RhoV1V2;
	    XL_readNumCell(XL_RhoV1V2,C_RhoV1V2,			" ARM_ERR: numeric expected",C_result);	

		long C_SabrType;
		CCString DefaultSabrType = "SABR_IMPLNVOL";
		XL_GETCONVSABRFLAGWD(XL_SabrType,C_SabrType, DefaultSabrType, " ARM_ERR: SabrType flag string expected",C_result);

		double C_GridSize;
		double C_DefaultGridSize = 251;
	    XL_readNumCellWD(XL_GridSize, C_GridSize, C_DefaultGridSize, " ARM_ERR: numeric expected",C_result);	

		exportFunc14Args< double, double, double, double, double, double, double, double, double, double, double, double, long, double > 
			ourFunc(C_Alpha1,C_Beta1, C_Rho1, C_Nu1, C_Alpha2,C_Beta2, C_Rho2, C_Nu2, C_RhoS1S2, C_RhoS1V2, C_RhoS2V1, C_RhoV1V2, C_SabrType, C_GridSize, ARMLOCAL_BiSABRDensityFunctor_Create);

		/// call the general function
		fillXL_Result_withName(ourFunc, C_result, XL_result, PersistentInXL );

	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_NumericalCalibMethodCommon" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;

}

__declspec(dllexport) LPXLOPER WINAPI Local_PXL_BiSABRDensityFunctor(
    LPXLOPER XL_Alpha1,
	LPXLOPER XL_Beta1,
	LPXLOPER XL_Rho1,
	LPXLOPER XL_Nu1,
    LPXLOPER XL_Alpha2,
	LPXLOPER XL_Beta2,
	LPXLOPER XL_Rho2,
	LPXLOPER XL_Nu2,
	LPXLOPER XL_RhoS1S2,
	LPXLOPER XL_RhoS1V2,
	LPXLOPER XL_RhoS2V1,
	LPXLOPER XL_RhoV1V2,
	LPXLOPER XL_SabrType,
	LPXLOPER XL_GridSize)
{
	ADD_LOG("Local_PXL_BiSABRDensityFunctor");
	bool PersistentInXL = false;
	return Local_BiSABRDensityFunctorCommon(XL_Alpha1, XL_Beta1, XL_Rho1, XL_Nu1, XL_Alpha2, XL_Beta2, XL_Rho2, XL_Nu2, XL_RhoS1S2, XL_RhoS1V2, XL_RhoS2V1, XL_RhoV1V2, XL_SabrType, XL_GridSize, PersistentInXL);
}

__declspec(dllexport) LPXLOPER WINAPI Local_BiSABRDensityFunctor(
    LPXLOPER XL_Alpha1,
	LPXLOPER XL_Beta1,
	LPXLOPER XL_Rho1,
	LPXLOPER XL_Nu1,
    LPXLOPER XL_Alpha2,
	LPXLOPER XL_Beta2,
	LPXLOPER XL_Rho2,
	LPXLOPER XL_Nu2,
	LPXLOPER XL_RhoS1S2,
	LPXLOPER XL_RhoS1V2,
	LPXLOPER XL_RhoS2V1,
	LPXLOPER XL_RhoV1V2,
	LPXLOPER XL_SabrType,
	LPXLOPER XL_GridSize)
{
	ADD_LOG("Local_BiSABRDensityFunctor");
	bool PersistentInXL = true;
	return Local_BiSABRDensityFunctorCommon(XL_Alpha1, XL_Beta1, XL_Rho1, XL_Nu1, XL_Alpha2, XL_Beta2, XL_Rho2, XL_Nu2, XL_RhoS1S2, XL_RhoS1V2, XL_RhoS2V1, XL_RhoV1V2, XL_SabrType, XL_GridSize, PersistentInXL);
}

LPXLOPER WINAPI Local_NormalHestonDensityFunctor_Common(
	LPXLOPER XL_Fwd,
	LPXLOPER XL_V0,
	LPXLOPER XL_Kappa,
	LPXLOPER XL_Theta,
	LPXLOPER XL_VVol,
	LPXLOPER XL_Rho,
	LPXLOPER XL_Level,
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

	    double C_Fwd, C_V0, C_Kappa, C_Theta, C_VVol, C_Rho, C_Level, defLevel = 1.;

	    XL_readNumCell(XL_Fwd, C_Fwd," ARM_ERR: forward numeric expected",C_result);	
		XL_readNumCell(XL_V0, C_V0," ARM_ERR: v0 numeric expected",C_result);	
		XL_readNumCell(XL_Kappa, C_Kappa," ARM_ERR: kappa numeric expected",C_result);	
		XL_readNumCell(XL_Theta, C_Theta," ARM_ERR: theta numeric expected",C_result);	
		XL_readNumCell(XL_VVol, C_VVol," ARM_ERR: vvol numeric expected",C_result);	
		XL_readNumCell(XL_Rho, C_Rho," ARM_ERR: rho numeric expected",C_result);	
		XL_readNumCellWD(XL_Level, C_Level, defLevel, " ARM_ERR: level numeric expected",C_result);	

		exportFunc7Args< double, double, double, double, double, double, double > 
			ourFunc(
				C_Fwd,
				C_V0,
				C_Kappa,
				C_Theta,
				C_VVol,
				C_Rho,
				C_Level,
				ARMLOCAL_NormalHestonDensityFunctor_Create);

		/// call the general function
		fillXL_Result_withName(ourFunc, C_result, XL_result, PersistentInXL );

	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_NumericalCalibMethodCommon" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}

__declspec(dllexport) LPXLOPER WINAPI Local_NormalHestonDensityFunctor(
	LPXLOPER XL_Fwd,
	LPXLOPER XL_V0,
	LPXLOPER XL_Kappa,
	LPXLOPER XL_Theta,
	LPXLOPER XL_VVol,
	LPXLOPER XL_Rho,
	LPXLOPER XL_Level)
{
	ADD_LOG("Local_NormalHestonDensityFunctor");
	bool PersistentInXL = true;

	return Local_NormalHestonDensityFunctor_Common(XL_Fwd, XL_V0, XL_Kappa, XL_Theta, XL_VVol, XL_Rho, XL_Level, PersistentInXL);
}

__declspec(dllexport) LPXLOPER WINAPI Local_PXL_NormalHestonDensityFunctor(
	LPXLOPER XL_Fwd,
	LPXLOPER XL_V0,
	LPXLOPER XL_Kappa,
	LPXLOPER XL_Theta,
	LPXLOPER XL_VVol,
	LPXLOPER XL_Rho,
	LPXLOPER XL_Level)
{
	ADD_LOG("Local_PXL_NormalHestonDensityFunctor");
	bool PersistentInXL = false;

	return Local_NormalHestonDensityFunctor_Common(XL_Fwd, XL_V0, XL_Kappa, XL_Theta, XL_VVol, XL_Rho, XL_Level, PersistentInXL);
}


// SPLINE DENSITY FUNCTOR
LPXLOPER Local_SplineDensityFunctorCommon(
    LPXLOPER XL_Money,
	LPXLOPER XL_Vol,
	LPXLOPER XL_VolType,
	LPXLOPER XL_SmileId,
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

		vector<double> defaultVector(0);
		defaultVector.clear();

		CCString C_VolTypeStr;
		XL_readStrCellWD( XL_VolType,C_VolTypeStr,"GAUSS"," ARM_ERR: VolType: string expected",C_result);
		string C_VolType = CCSTringToSTLString(C_VolTypeStr);

		vector<double> C_Money;
		XL_readNumVectorWD(XL_Money,C_Money,defaultVector," ARM_ERR: Moneyness: array of numeric expected",C_result);

		vector<double> C_Vol;
		XL_readNumVectorWD(XL_Vol,C_Vol,defaultVector," ARM_ERR: Vol: array of numeric expected",C_result);

		if (C_Vol.size()!=C_Money.size())
			throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"vol and moneyness must have the same size");

		CCString C_SmileStrId;
		XL_readStrCellWD(XL_SmileId, C_SmileStrId, GETDEFAULTVALUESTR," ARM_ERR: smile: string expected", C_result);
		long C_SmileId ;
		if (C_SmileStrId == GETDEFAULTVALUESTR)
			C_SmileId = ARM_NULL_OBJECT;    
		else
			C_SmileId = LocalGetNumObjectId(C_SmileStrId);
		
		exportFunc4Args< vector<double>, vector<double>, string, long > 
			ourFunc(C_Money,C_Vol, C_VolType, C_SmileId,ARMLOCAL_SplineDensityFunctor_Create);
		
		/// call the general function
		fillXL_Result_withName(ourFunc, C_result, XL_result, PersistentInXL );

	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_NumericalCalibMethodCommon" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;

}

__declspec(dllexport) LPXLOPER WINAPI Local_PXL_SplineDensityFunctor(
    LPXLOPER XL_Money,
	LPXLOPER XL_Vol,
	LPXLOPER XL_VolType,
	LPXLOPER XL_Smile)
{
	ADD_LOG("Local_PXL_SplineDensityFunctor");
	bool PersistentInXL = false;
	return Local_SplineDensityFunctorCommon(XL_Money, XL_Vol, XL_VolType, XL_Smile, PersistentInXL);
}

__declspec(dllexport) LPXLOPER WINAPI Local_SplineDensityFunctor(
    LPXLOPER XL_Money,
	LPXLOPER XL_Vol,
	LPXLOPER XL_VolType,
	LPXLOPER XL_Smile)
{
	ADD_LOG("Local_SplineDensityFunctor");
	bool PersistentInXL = true;
	return Local_SplineDensityFunctorCommon(XL_Money, XL_Vol, XL_VolType, XL_Smile, PersistentInXL);
}

__declspec(dllexport) LPXLOPER WINAPI Local_DensityFunctor_CallOption(
	LPXLOPER XL_DensityFunctorId,
    LPXLOPER XL_Forward,
	LPXLOPER XL_Strike,
	LPXLOPER XL_Maturity)
{
	ADD_LOG("Local_DensityFunctor_CallOption");
	/// this is defined first because it is used in XL macros
	static XLOPER XL_result;
	/// Get the variables from the XLOper variables
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{/// to avoid computation if called by the wizard
	    ARM_NOCALCIFWIZ();
	    
	    /// this is used by macros 
	    /// and therefore this has to be defined
	    static int error;
		static char* reason = "";	
		
		double C_Forward;
	    XL_readNumCell(XL_Forward,C_Forward,		" ARM_ERR: forward: number expected",C_result);	

	    double C_Strike;
	    XL_readNumCell(XL_Strike,C_Strike,			" ARM_ERR: strike: number expected",C_result);	
		
		double C_Maturity;
	    XL_readNumCell(XL_Maturity,C_Maturity,		" ARM_ERR: maturity: number expected",C_result);	
		
		CCString C_DensityFunctorStrId;
		XL_readStrCell( XL_DensityFunctorId, C_DensityFunctorStrId,	" ARM_ERR: Density Functor Id: Object expected", C_result);
		long C_DensityFunctorId = LocalGetNumObjectId(C_DensityFunctorStrId);


		long retCode = ARMLOCAL_DensityFunctor_CallOption(
				C_DensityFunctorId,
				C_Forward,
				C_Strike,
				C_Maturity,
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_VanillaOption_Normal" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}

__declspec(dllexport) LPXLOPER WINAPI Local_DensityFunctor_Quantile(
	LPXLOPER XL_DensityFunctorId,
    LPXLOPER XL_Forward,
	LPXLOPER XL_Proba,
	LPXLOPER XL_Maturity)
{
	ADD_LOG("Local_DensityFunctor_Quantile");
	/// this is defined first because it is used in XL macros
	static XLOPER XL_result;
	/// Get the variables from the XLOper variables
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{/// to avoid computation if called by the wizard
	    ARM_NOCALCIFWIZ();
	    
	    /// this is used by macros 
	    /// and therefore this has to be defined
	    static int error;
		static char* reason = "";	
		
		double C_Forward;
	    XL_readNumCell(XL_Forward,C_Forward,		" ARM_ERR: forward: number expected",C_result);	

	    double C_Proba;
	    XL_readNumCell(XL_Proba,C_Proba,			" ARM_ERR: strike: number expected",C_result);	
		
		double C_Maturity;
	    XL_readNumCell(XL_Maturity,C_Maturity,		" ARM_ERR: maturity: number expected",C_result);	
		
		CCString C_DensityFunctorStrId;
		XL_readStrCell( XL_DensityFunctorId, C_DensityFunctorStrId,	" ARM_ERR: Density Functor Id: Object expected", C_result);
		long C_DensityFunctorId = LocalGetNumObjectId(C_DensityFunctorStrId);


		long retCode = ARMLOCAL_DensityFunctor_Quantile(
				C_DensityFunctorId,
				C_Forward,
				C_Proba,
				C_Maturity,
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_VanillaOption_Normal" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}
////////////////////////////////////////////////////////
/// Create Vanilla Caplet Density (HK)
////////////////////////////////////////////////////////

LPXLOPER WINAPI Local_VanillaSecurityDensityCommon(
    LPXLOPER XL_ResetDate,
    LPXLOPER XL_StartDate,
    LPXLOPER XL_EndDate,
	LPXLOPER XL_DensityFunctorId,
	LPXLOPER XL_Frequency,
	LPXLOPER XL_DayCount,
	LPXLOPER XL_StubRule,
	LPXLOPER XL_Weight,
	LPXLOPER XL_AdjFwdAdd,
	LPXLOPER XL_AdjFwdMult,
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

		double C_ResetDate;
	    XL_readNumCell(XL_ResetDate,C_ResetDate,		" ARM_ERR: reset date: date expected",C_result);	

	    double C_StartDate;
	    XL_readNumCell(XL_StartDate,C_StartDate,		" ARM_ERR: start date: date expected",C_result);	

	    double C_EndDate;
	    XL_readNumCell(XL_EndDate,C_EndDate,			" ARM_ERR: end date: date expected",C_result);	
		
		CCString C_DensityFunctorStrId;
		XL_readStrCell( XL_DensityFunctorId, C_DensityFunctorStrId,	" ARM_ERR: Density Functor Id: Object expected", C_result);
		long C_DensityFunctorId = LocalGetNumObjectId(C_DensityFunctorStrId);

		CCString C_FrequencyStr;
		CCString C_DayCountStr;
		CCString C_StubRuleStr;
	
		XL_readStrCellWD(XL_Frequency, C_FrequencyStr,GETDEFAULTVALUESTR," ARM_ERR: frequency: string expected", C_result);
		XL_readStrCellWD(XL_DayCount,  C_DayCountStr, GETDEFAULTVALUESTR," ARM_ERR: day count: string expected", C_result);
		XL_readStrCellWD(XL_StubRule,  C_StubRuleStr, GETDEFAULTVALUESTR," ARM_ERR: stub rule: string expected", C_result);
		
		double C_Weight, C_DefaultWeight = 1., C_FwdAdd, C_DefaultAdd = 0., C_FwdMult, C_DefaultMult = 1.;
		XL_readNumCellWD(XL_Weight, C_Weight, C_DefaultWeight, "ARM_ERR : weight : numeric expected", C_result);
		XL_readNumCellWD(XL_AdjFwdAdd, C_FwdAdd, C_DefaultAdd, "ARM_ERR : fwd add : numeric expected", C_result);
		XL_readNumCellWD(XL_AdjFwdMult, C_FwdMult, C_DefaultMult, "ARM_ERR : fwd mult : numeric expected", C_result);

		long C_Frequency;
		long C_DayCount;
		long C_StubRule;
	
		if (C_FrequencyStr == GETDEFAULTVALUESTR)
			C_Frequency = GETDEFAULTVALUE;
		else if ( (C_Frequency = ARM_ConvFrequency (C_FrequencyStr, C_result)) == ARM_DEFAULT_ERR )
		{
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}

		if (C_DayCountStr == GETDEFAULTVALUESTR)
			C_DayCount = GETDEFAULTVALUE;
		else 
			C_DayCount = ARM_ConvDayCount (C_DayCountStr);

		if (C_StubRuleStr == GETDEFAULTVALUESTR)
			C_StubRule = GETDEFAULTVALUE;
		else 
			C_StubRule = ARM_ConvStubRule (C_StubRuleStr);

		exportFunc10Args< double, double, double, long, long, long, long, double, double, double> 
			ourFunc(C_ResetDate,C_StartDate,C_EndDate,C_DensityFunctorId,C_Frequency,C_DayCount,C_StubRule,C_Weight,C_FwdAdd, C_FwdMult, ARMLOCAL_VanillaSecurityDensity_Create);	

		/// call the general function
		fillXL_Result_withName(ourFunc, C_result, XL_result, PersistentInXL );
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_VanillaSecurityDensityCommon" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


__declspec(dllexport) LPXLOPER WINAPI Local_PXL_VanillaSecurityDensity(
    LPXLOPER XL_ResetDate,
    LPXLOPER XL_StartDate,
    LPXLOPER XL_EndDate,
	LPXLOPER XL_DensityFunctorId,
	LPXLOPER XL_Frequency,
	LPXLOPER XL_DayCount,
	LPXLOPER XL_StubRule,
	LPXLOPER XL_Weight,
	LPXLOPER XL_AdjFwdAdd,
	LPXLOPER XL_AdjFwdMult)

{
	ADD_LOG("Local_PXL_VanillaSecurityDensity");
	bool PersistentInXL = false;
	return Local_VanillaSecurityDensityCommon(XL_ResetDate,XL_StartDate,XL_EndDate, XL_DensityFunctorId, XL_Frequency,XL_DayCount,XL_StubRule,XL_Weight,XL_AdjFwdAdd, XL_AdjFwdMult, PersistentInXL );

}

__declspec(dllexport) LPXLOPER WINAPI Local_VanillaSecurityDensity(
    LPXLOPER XL_ResetDate,
    LPXLOPER XL_StartDate,
    LPXLOPER XL_EndDate,
	LPXLOPER XL_DensityFunctorId,
	LPXLOPER XL_Frequency,
	LPXLOPER XL_DayCount,
	LPXLOPER XL_StubRule,
	LPXLOPER XL_Weight,
	LPXLOPER XL_AdjFwdAdd,
	LPXLOPER XL_AdjFwdMult)
{
	ADD_LOG("Local_VanillaSecurityDensity");
	bool PersistentInXL = true;
	return Local_VanillaSecurityDensityCommon(XL_ResetDate,XL_StartDate,XL_EndDate,XL_DensityFunctorId,XL_Frequency,XL_DayCount,XL_StubRule,XL_Weight,XL_AdjFwdAdd, XL_AdjFwdMult, PersistentInXL );
}

LPXLOPER WINAPI Local_VanillaSecurityDensitySpreadCommon(
    LPXLOPER XL_ResetDate,
    LPXLOPER XL_StartDate1,
    LPXLOPER XL_EndDate1,
    LPXLOPER XL_StartDate2,
    LPXLOPER XL_EndDate2,
	LPXLOPER XL_DensityFunctorId,
	LPXLOPER XL_Frequency1,
	LPXLOPER XL_DayCount1,
	LPXLOPER XL_Frequency2,
	LPXLOPER XL_DayCount2,
	LPXLOPER XL_StubRule,
	LPXLOPER XL_Weight,
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

		double C_ResetDate;
	    XL_readNumCell(XL_ResetDate,C_ResetDate,		" ARM_ERR: reset date: date expected",C_result);	

	    double C_StartDate1;
	    XL_readNumCell(XL_StartDate1,C_StartDate1,		" ARM_ERR: start date: date expected",C_result);	

	    double C_EndDate1;
	    XL_readNumCell(XL_EndDate1,C_EndDate1,			" ARM_ERR: end date: date expected",C_result);	
		
	    double C_StartDate2;
	    XL_readNumCell(XL_StartDate2,C_StartDate2,		" ARM_ERR: start date: date expected",C_result);	

	    double C_EndDate2;
	    XL_readNumCell(XL_EndDate2,C_EndDate2,			" ARM_ERR: end date: date expected",C_result);	

		CCString C_DensityFunctorStrId;
		XL_readStrCell( XL_DensityFunctorId, C_DensityFunctorStrId,	" ARM_ERR: Density Functor Id: Object expected", C_result);
		long C_DensityFunctorId = LocalGetNumObjectId(C_DensityFunctorStrId);

		CCString C_FrequencyStr1;
		CCString C_DayCountStr1;
		CCString C_FrequencyStr2;
		CCString C_DayCountStr2;
		CCString C_StubRuleStr;
	
		XL_readStrCellWD(XL_Frequency1, C_FrequencyStr1,GETDEFAULTVALUESTR," ARM_ERR: frequency: string expected", C_result);
		XL_readStrCellWD(XL_DayCount1,  C_DayCountStr1, GETDEFAULTVALUESTR," ARM_ERR: day count: string expected", C_result);
		XL_readStrCellWD(XL_Frequency2, C_FrequencyStr2,GETDEFAULTVALUESTR," ARM_ERR: frequency: string expected", C_result);
		XL_readStrCellWD(XL_DayCount2,  C_DayCountStr2, GETDEFAULTVALUESTR," ARM_ERR: day count: string expected", C_result);
		XL_readStrCellWD(XL_StubRule,  C_StubRuleStr, GETDEFAULTVALUESTR," ARM_ERR: stub rule: string expected", C_result);

		double C_Weight, C_DefaultWeight = 1.;
		XL_readNumCellWD(XL_Weight, C_Weight, C_DefaultWeight, "ARM_ERR : weight : numeric expected", C_result);

		long C_Frequency1;
		long C_DayCount1;
		long C_Frequency2;
		long C_DayCount2;
		long C_StubRule;
	
		if (C_FrequencyStr1 == GETDEFAULTVALUESTR)
			C_Frequency1 = GETDEFAULTVALUE;
		else if ( (C_Frequency1 = ARM_ConvFrequency (C_FrequencyStr1, C_result)) == ARM_DEFAULT_ERR )
		{
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}

		if (C_DayCountStr1 == GETDEFAULTVALUESTR)
			C_DayCount1 = GETDEFAULTVALUE;
		else 
			C_DayCount1 = ARM_ConvDayCount (C_DayCountStr1);

		if (C_FrequencyStr2 == GETDEFAULTVALUESTR)
			C_Frequency2 = GETDEFAULTVALUE;
		else if ( (C_Frequency2 = ARM_ConvFrequency (C_FrequencyStr2, C_result)) == ARM_DEFAULT_ERR )
		{
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}

		if (C_DayCountStr2 == GETDEFAULTVALUESTR)
			C_DayCount2 = GETDEFAULTVALUE;
		else 
			C_DayCount2 = ARM_ConvDayCount (C_DayCountStr2);

		if (C_StubRuleStr == GETDEFAULTVALUESTR)
			C_StubRule = GETDEFAULTVALUE;
		else 
			C_StubRule = ARM_ConvStubRule (C_StubRuleStr);

		exportFunc12Args< double, double, double, double, double, long, long, long, long, long, long, double > 
			ourFunc(C_ResetDate,C_StartDate1,C_EndDate1,C_StartDate2, C_EndDate2, C_DensityFunctorId,C_Frequency1,C_DayCount1,C_Frequency2,C_DayCount2, C_StubRule,C_Weight,ARMLOCAL_VanillaSecurityDensitySpread_Create);	

		/// call the general function
		fillXL_Result_withName(ourFunc, C_result, XL_result, PersistentInXL );
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_VanillaSecurityDensitySpreadCommon" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


__declspec(dllexport) LPXLOPER WINAPI Local_PXL_VanillaSecurityDensitySpread(
    LPXLOPER XL_ResetDate,
    LPXLOPER XL_StartDate1,
    LPXLOPER XL_EndDate1,
    LPXLOPER XL_StartDate2,
    LPXLOPER XL_EndDate2,
	LPXLOPER XL_DensityFunctorId,
	LPXLOPER XL_Frequency1,
	LPXLOPER XL_DayCount1,
	LPXLOPER XL_Frequency2,
	LPXLOPER XL_DayCount2,
	LPXLOPER XL_StubRule,
	LPXLOPER XL_Weight)
{
	ADD_LOG("Local_PXL_VanillaSecurityDensitySpread");
	bool PersistentInXL = false;
	return Local_VanillaSecurityDensitySpreadCommon(XL_ResetDate,XL_StartDate1,XL_EndDate1, XL_StartDate2,XL_EndDate2, XL_DensityFunctorId, XL_Frequency1,XL_DayCount1,XL_Frequency2,XL_DayCount2,XL_StubRule,XL_Weight,PersistentInXL );

}

__declspec(dllexport) LPXLOPER WINAPI Local_VanillaSecurityDensitySpread(
    LPXLOPER XL_ResetDate,
    LPXLOPER XL_StartDate1,
    LPXLOPER XL_EndDate1,
    LPXLOPER XL_StartDate2,
    LPXLOPER XL_EndDate2,
	LPXLOPER XL_DensityFunctorId,
	LPXLOPER XL_Frequency1,
	LPXLOPER XL_DayCount1,
	LPXLOPER XL_Frequency2,
	LPXLOPER XL_DayCount2,
	LPXLOPER XL_StubRule,
	LPXLOPER XL_Weight)
{
	ADD_LOG("Local_VanillaSecurityDensitySpread");
	bool PersistentInXL = true;
	return Local_VanillaSecurityDensitySpreadCommon(XL_ResetDate,XL_StartDate1,XL_EndDate1, XL_StartDate2,XL_EndDate2, XL_DensityFunctorId, XL_Frequency1,XL_DayCount1,XL_Frequency2,XL_DayCount2,XL_StubRule,XL_Weight,PersistentInXL );
}

LPXLOPER WINAPI Local_VanillaSecurityDensityFXCommon(
    LPXLOPER XL_ResetDate,
	LPXLOPER XL_DensityFunctorId,
	LPXLOPER XL_DomCurveId,
	LPXLOPER XL_ForCurveId,
	LPXLOPER XL_FXSpot,
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

		double C_ResetDate;
	    XL_readNumCell(XL_ResetDate,C_ResetDate,		" ARM_ERR: reset date: date expected",C_result);
		
		CCString C_DensityFunctorStrId;
		XL_readStrCell( XL_DensityFunctorId, C_DensityFunctorStrId,	" ARM_ERR: Density Functor Id: Object expected", C_result);
		long C_DensityFunctorId = LocalGetNumObjectId(C_DensityFunctorStrId);

		CCString C_DomCurveStrId;
		XL_readStrCell( XL_DomCurveId, C_DomCurveStrId,	" ARM_ERR: Domestic Curve Id: Object expected", C_result);
		long C_DomCurveId = LocalGetNumObjectId(C_DomCurveStrId);

		CCString C_ForCurveStrId;
		XL_readStrCell( XL_ForCurveId, C_ForCurveStrId,	" ARM_ERR: Foreign Curve Id: Object expected", C_result);
		long C_ForCurveId = LocalGetNumObjectId(C_ForCurveStrId);

		double C_FXSpot;
		XL_readNumCell(XL_FXSpot,C_FXSpot," ARM_ERR: FX Spot: numerical expected",C_result);

		exportFunc5Args< double, long, long, long, double > 
			ourFunc(C_ResetDate,C_DensityFunctorId,C_DomCurveId,C_ForCurveId,C_FXSpot,ARMLOCAL_VanillaSecurityDensityFX_Create);

		/// call the general function
		fillXL_Result_withName(ourFunc, C_result, XL_result, PersistentInXL );
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_VanillaSecurityDensityCommon" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


__declspec(dllexport) LPXLOPER WINAPI Local_PXL_VanillaSecurityDensityFX(
    LPXLOPER XL_ResetDate,
	LPXLOPER XL_DensityFunctorId,
	LPXLOPER XL_DomCurveId,
	LPXLOPER XL_ForCurveId,
	LPXLOPER XL_FXSpot)
{
	ADD_LOG("Local_PXL_VanillaSecurityDensityFX");
	bool PersistentInXL = false;
	return Local_VanillaSecurityDensityFXCommon(XL_ResetDate,XL_DensityFunctorId,XL_DomCurveId,XL_ForCurveId,XL_FXSpot,PersistentInXL );

}

__declspec(dllexport) LPXLOPER WINAPI Local_VanillaSecurityDensityFX(
    LPXLOPER XL_ResetDate,
	LPXLOPER XL_DensityFunctorId,
	LPXLOPER XL_DomCurveId,
	LPXLOPER XL_ForCurveId,
	LPXLOPER XL_FXSpot)
{
	ADD_LOG("Local_VanillaSecurityDensityFX");
	bool PersistentInXL = true;
	return Local_VanillaSecurityDensityFXCommon(XL_ResetDate,XL_DensityFunctorId,XL_DomCurveId,XL_ForCurveId,XL_FXSpot,PersistentInXL );
}

/////////////////////////////////////////////////////////////////////////////::
/////////////////////////////////////////////////////////////////////////////::
/////////   For HW2F only
/////////////////////////////////////////////////////////////////////////////::
/////////////////////////////////////////////////////////////////////////////::

/////////////////////////////////////////////////////////////
/// central function that does the creation of the XL function
/////////////////////////////////////////////////////////////


LPXLOPER Local_CalibMethodHW2F_CreateCommon(
	LPXLOPER XL_PortfolioId1,
    LPXLOPER XL_ParamId1,
    LPXLOPER XL_PortfolioId2,
	LPXLOPER XL_ParamId2,
    LPXLOPER XL_PortfolioId3,
	LPXLOPER XL_ParamId3,
	LPXLOPER XL_FlagOptim,
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
		
		CCString C_PortfolioStrId1;
		XL_readStrCellWD(XL_PortfolioId1, C_PortfolioStrId1, GETDEFAULTVALUESTR," ARM_ERR: porfolio1: string expected", C_result);
		long C_PortfolioId1 ;
		if (C_PortfolioStrId1 == GETDEFAULTVALUESTR)
			C_PortfolioId1 = ARM_NULL_OBJECT;    
		else
			C_PortfolioId1 = LocalGetNumObjectId(C_PortfolioStrId1);

		CCString C_PortfolioStrId2;
		XL_readStrCellWD(XL_PortfolioId2, C_PortfolioStrId2, GETDEFAULTVALUESTR," ARM_ERR: porfolio2: string expected", C_result);
		long C_PortfolioId2 ;
		if (C_PortfolioStrId2 == GETDEFAULTVALUESTR)
			C_PortfolioId2 = ARM_NULL_OBJECT;    
		else
			C_PortfolioId2 = LocalGetNumObjectId(C_PortfolioStrId2);

		CCString C_PortfolioStrId3;
		XL_readStrCellWD(XL_PortfolioId3, C_PortfolioStrId3, GETDEFAULTVALUESTR," ARM_ERR: porfolio3: string expected", C_result);
		long C_PortfolioId3 ;
		if (C_PortfolioStrId3 == GETDEFAULTVALUESTR)
			C_PortfolioId3 = ARM_NULL_OBJECT;    
		else
			C_PortfolioId3 = LocalGetNumObjectId(C_PortfolioStrId3);

		CCString C_ParamStrId1;
		XL_readStrCellWD(XL_ParamId1, C_ParamStrId1, GETDEFAULTVALUESTR," ARM_ERR: porfolio1: string expected", C_result);
		long C_ParamId1 ;
		if (C_ParamStrId1 == GETDEFAULTVALUESTR)
			C_ParamId1 = ARM_NULL_OBJECT;    
		else
			C_ParamId1 = LocalGetNumObjectId(C_ParamStrId1);

		CCString C_ParamStrId2;
		XL_readStrCellWD(XL_ParamId2, C_ParamStrId2, GETDEFAULTVALUESTR," ARM_ERR: porfolio2: string expected", C_result);
		long C_ParamId2 ;
		if (C_ParamStrId2 == GETDEFAULTVALUESTR)
			C_ParamId2 = ARM_NULL_OBJECT;    
		else
			C_ParamId2 = LocalGetNumObjectId(C_ParamStrId2);

		CCString C_ParamStrId3;
		XL_readStrCellWD(XL_ParamId3, C_ParamStrId3, GETDEFAULTVALUESTR," ARM_ERR: porfolio3: string expected", C_result);
		long C_ParamId3 ;
		if (C_ParamStrId3 == GETDEFAULTVALUESTR)
			C_ParamId3 = ARM_NULL_OBJECT;    
		else
			C_ParamId3 = LocalGetNumObjectId(C_ParamStrId3);

		double C_WithOptim;
		double withOptimDefault=0.;
		XL_readNumCellWD( XL_FlagOptim, C_WithOptim, withOptimDefault," ARM_ERR: Nb of factors: numeric expected",C_result);
		int C_WithOptimInt = C_WithOptim;

		/// use the concept of Functor to transfer the knowledge of
		/// a function with a context
		exportFunc7Args< long, long, long, long, long, long, int > 
			ourFunc(C_PortfolioId1,C_ParamId1,C_PortfolioId2,C_ParamId2,C_PortfolioId3,C_ParamId3,C_WithOptimInt,ARMLOCAL_CalibMethodHW2F_Create);

		/// call the general function
		fillXL_Result( LOCAL_CALIBMETHOD_CLASS, ourFunc, C_result, XL_result, PersistentInXL );
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_NumericalCalibMethodCommon" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}



///////////////////////////////////
/// Create a Calib Method fo HW2F
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_CalibMethodHW2F_Create(
    LPXLOPER XL_PortfolioId1,
    LPXLOPER XL_ParamId1,
    LPXLOPER XL_PortfolioId2,
	LPXLOPER XL_ParamId2,
    LPXLOPER XL_PortfolioId3,
	LPXLOPER XL_ParamId3,
	LPXLOPER XL_FlagOptim)
{
	ADD_LOG("Local_CalibMethodHW2F_Create");
	bool PersistentInXL = true;
	return Local_CalibMethodHW2F_CreateCommon(
                                    XL_PortfolioId1,
									XL_ParamId1,
									XL_PortfolioId2, 
									XL_ParamId2,
									XL_PortfolioId3,
									XL_ParamId3,
									XL_FlagOptim,
									PersistentInXL );
}

__declspec(dllexport) LPXLOPER WINAPI Local_PXL_CalibMethodHW2F_Create(
	LPXLOPER XL_PortfolioId1,
    LPXLOPER XL_ParamId1,
    LPXLOPER XL_PortfolioId2,
	LPXLOPER XL_ParamId2,
    LPXLOPER XL_PortfolioId3,
	LPXLOPER XL_ParamId3,
	LPXLOPER XL_FlagOptim)
{
	ADD_LOG("Local_PXL_CalibMethodHW2F_Create");
	bool PersistentInXL = false;
	return Local_CalibMethodHW2F_CreateCommon(
                                    XL_PortfolioId1,
									XL_ParamId1,
									XL_PortfolioId2, 
									XL_ParamId2,
									XL_PortfolioId3,
									XL_ParamId3,
									XL_FlagOptim,
									PersistentInXL );
}



/////////////////////////////////////////////////////////////
/// central function that computes basket decomposition of a security
/////////////////////////////////////////////////////////////


LPXLOPER Local_BasketDecomp_CreateCommon(
	LPXLOPER XL_SecuritiesId,
	LPXLOPER XL_ModelsId,
    LPXLOPER XL_Weights,
	LPXLOPER XL_NotionalId,
    LPXLOPER XL_ExerDateStripId,
	LPXLOPER XL_ExerFeesId,
	LPXLOPER XL_Side,
	LPXLOPER XL_MkmoId,
	LPXLOPER XL_Method,
	LPXLOPER XL_Strike,
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
		
		VECTOR<CCString> C_SecuritiesId;
		XL_readStrVector (XL_SecuritiesId,C_SecuritiesId," ARM_ERR: Securities: array of object expected",DOUBLE_TYPE,C_result);
		
		VECTOR<CCString> C_ModelsId;
		XL_readStrVector (XL_ModelsId,C_ModelsId," ARM_ERR: Models: array of object expected",DOUBLE_TYPE,C_result);

		VECTOR < double > C_Weights;
		XL_readNumVector(XL_Weights,C_Weights," ARM_ERR: Weights is not of a good type",C_result);

		size_t size = C_SecuritiesId.size();

		if ( C_Weights.size()!=size || C_ModelsId.size() !=size )
			throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"securities, models and weights must have teh same size!!");
		
		VECTOR<long> C_SecuritiesIdLong (size);
		VECTOR<long> C_ModelsIdLong (size);
		size_t i;
		for(i = 0; i < size; ++i )
		{
			C_SecuritiesIdLong[i] = LocalGetNumObjectId(C_SecuritiesId[i]); 
			C_ModelsIdLong[i] = LocalGetNumObjectId(C_ModelsId[i]); 
		}

		CCString C_ExerDateStripStrId;
		XL_readStrCell( XL_ExerDateStripId, C_ExerDateStripStrId,	" ARM_ERR: Datestrip Id: Object expected",		C_result);
		long C_ExerDateStripId = LocalGetNumObjectId(C_ExerDateStripStrId);
		

		//NotionalValue or NotionalCurve
		double notional;
		CCString notionalStr;
		long notionalId;
		XL_readStrOrNumCell(XL_NotionalId, notionalStr, notional, notionalId,
			   " ARM_ERR: notional: numeric or refValue Id expected",C_result);	

		if (notionalId == XL_TYPE_STRING)
			notionalId = LocalGetNumObjectId(notionalStr);
		else
			notionalId = ARM_NULL_OBJECT;

		//ExerFeesValue or ExerFeesCurve
		double ExerFees;
		CCString ExerFeesStr;
		long ExerFeesId;
		XL_readStrOrNumCell(XL_ExerFeesId, ExerFeesStr, ExerFees, ExerFeesId,
			   " ARM_ERR: ExerFees: numeric or refValue Id expected",C_result);	

		if (ExerFeesId == XL_TYPE_STRING)
			ExerFeesId = LocalGetNumObjectId(ExerFeesStr);
		else
			ExerFeesId = ARM_NULL_OBJECT;

		//payer / receiver
		CCString C_SideStr;
		XL_readStrCellWD( XL_Side,C_SideStr,"R"," ARM_ERR: Side: string expected",C_result);
		C_SideStr.toUpper();
		double C_side;
		if (C_SideStr == "P" )
			C_side = 1.;
		else if (C_SideStr == "R")
			C_side = -1.;
		else
			throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"String \"P\" or \"R\" Expected for Side");

		CCString C_MkmoStrId;
		XL_readStrCell( XL_MkmoId, C_MkmoStrId,	" ARM_ERR: MkMo Id: Object expected",		C_result);
		long C_MkmoId = LocalGetNumObjectId(C_MkmoStrId);    

		CCString C_MethodStr;
		XL_readStrCellWD( XL_Method,C_MethodStr,"BASKET"," ARM_ERR: Method: string expected",C_result);
		string C_Method = CCSTringToSTLString(C_MethodStr);

		CCString C_StrikeStr;
		XL_readStrCellWD( XL_Strike,C_StrikeStr,"EQUIVALENT"," ARM_ERR: Strike: string expected",C_result);
		string C_Strike = CCSTringToSTLString(C_StrikeStr);

		/// use the concept of Functor to transfer the knowledge of
		/// a function with a context
		exportFunc12Args< VECTOR<long>, VECTOR<long>, long, long, VECTOR<double>, double, string, string, double, long, double, long> 
			ourFunc(C_SecuritiesIdLong,C_ModelsIdLong,C_ExerDateStripId, C_MkmoId, C_Weights,C_side,C_Method,C_Strike,notional,notionalId,ExerFees,ExerFeesId,ARMLOCAL_BasketDecomp_Create);

		/// call the general function
		fillXL_Result( LOCAL_BASKET_CLASS, ourFunc, C_result, XL_result, PersistentInXL );
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_BasketDecomp_CreateCommon" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


// Basket decomposition of a security
__declspec(dllexport) LPXLOPER WINAPI Local_BasketDecomp_Create(
	LPXLOPER XL_SecuritiesId,
	LPXLOPER XL_ModelsId,
    LPXLOPER XL_Weights,
	LPXLOPER XL_NotionalId,
    LPXLOPER XL_ExerDateStripId,
	LPXLOPER XL_ExerFeesId,
	LPXLOPER XL_Side,
	LPXLOPER XL_MkmoId,
	LPXLOPER XL_Method,
	LPXLOPER XL_Strike)
{
	ADD_LOG("Local_BasketDecomp_Create");
	bool PersistentInXL = true;
	return Local_BasketDecomp_CreateCommon(
                                    XL_SecuritiesId,
									XL_ModelsId,
									XL_Weights, 
									XL_NotionalId,
									XL_ExerDateStripId,
									XL_ExerFeesId,
									XL_Side,
									XL_MkmoId,
									XL_Method,
									XL_Strike,
									PersistentInXL );
}

__declspec(dllexport) LPXLOPER WINAPI Local_PXL_BasketDecomp_Create(
	LPXLOPER XL_SecuritiesId,
	LPXLOPER XL_ModelsId,
    LPXLOPER XL_Weights,
	LPXLOPER XL_NotionalId,
    LPXLOPER XL_ExerDateStripId,
	LPXLOPER XL_ExerFeesId,
	LPXLOPER XL_Side,
	LPXLOPER XL_MkmoId,
	LPXLOPER XL_Method,
	LPXLOPER XL_Strike)
{
	ADD_LOG("Local_PXL_BasketDecomp_Create");
	bool PersistentInXL = false;
	return Local_BasketDecomp_CreateCommon(
                                    XL_SecuritiesId,
									XL_ModelsId,
									XL_Weights, 
									XL_NotionalId,
									XL_ExerDateStripId,
									XL_ExerFeesId,
									XL_Side,
									XL_MkmoId,
									XL_Method,
									XL_Strike,
									PersistentInXL );
}

__declspec(dllexport) LPXLOPER WINAPI Local_GetData_FromBasket_Common(
	LPXLOPER XL_BasketId,
    LPXLOPER XL_KeyId, 
	bool PersistentInXL )
{
	ADD_LOG("Local_GetData_FromBasket_Common");
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
		
		CCString C_BasketString;
		XL_readStrCell( XL_BasketId, C_BasketString, " ARM_ERR: Basket Id: Object expected",C_result);
		long bsktId = LocalGetNumObjectId(C_BasketString);
		
		CCString C_KeyString;
		XL_readStrCellWD( XL_KeyId, C_KeyString, "PORTFOLIO", " ARM_ERR: Key Id: String expected",C_result);
		string keyString(CCSTringToSTLString(C_KeyString));

		CCString bsktGetClass(BsktGetTypeToClass(keyString,bsktId).c_str());

        exportFunc2Args< long, string > ourFunc(bsktId,keyString, ARMLOCAL_Basket_Get);

		/// call the general function
		fillXL_Result( bsktGetClass, ourFunc, C_result, XL_result, PersistentInXL );

	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_GetData_FromBasket_Common" )

	return (LPXLOPER)&XL_result;
}


__declspec(dllexport) LPXLOPER WINAPI Local_GetData_FromBasket(
	LPXLOPER XL_BasketId,
    LPXLOPER XL_KeyId)
{
	ADD_LOG("Local_GetData_FromBasket");
	return Local_GetData_FromBasket_Common( XL_BasketId, XL_KeyId, true );
}


__declspec(dllexport) LPXLOPER WINAPI Local_PXL_GetData_FromBasket(
	LPXLOPER XL_BasketId,
    LPXLOPER XL_KeyId)
{
	ADD_LOG("Local_PXL_GetData_FromBasket");
	return Local_GetData_FromBasket_Common( XL_BasketId, XL_KeyId, false );
}

/////////////////////////////////////////////////////////////
/// central function that computes basket decomposition of a security
/////////////////////////////////////////////////////////////


LPXLOPER Local_SmileViewer_CreateCommon(
	LPXLOPER XL_SecurityId,
	LPXLOPER XL_ModelId,
	LPXLOPER XL_Moneyness,
	LPXLOPER XL_MoneyType,
	LPXLOPER XL_Strike,
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
		
		CCString C_SecString;
		XL_readStrCell( XL_SecurityId, C_SecString, " ARM_ERR: Security Id: Object expected",C_result);
		long C_SecurityIdLong = LocalGetNumObjectId(C_SecString);

		CCString C_ModString;
		XL_readStrCell( XL_ModelId, C_ModString, " ARM_ERR: Model Id: Object expected",C_result);
		long C_ModelIdLong = LocalGetNumObjectId(C_ModString);
		
		VECTOR < double > C_Moneyness;
		VECTOR < double > C_MoneynessDef(0);
		XL_readNumVectorWD(XL_Moneyness,C_Moneyness,C_MoneynessDef," ARM_ERR: Moneyness is not of a good type",C_result);

		CCString C_MethodStr;
		XL_readStrCellWD( XL_MoneyType,C_MethodStr,"GAUSS"," ARM_ERR: Method: string expected",C_result);
		string C_Method = CCSTringToSTLString(C_MethodStr);
		
		VECTOR < double > C_Strike;
		VECTOR < double > C_StrikeDef(0);
		XL_readNumVectorWD(XL_Strike,C_Strike,C_StrikeDef," ARM_ERR: Strike is not of a good type",C_result);

		/// use the concept of Functor to transfer the knowledge of
		/// a function with a context
		exportFunc5Args< long, long,VECTOR < double >,string, VECTOR <double> > 
			ourFunc(C_SecurityIdLong,C_ModelIdLong,C_Moneyness,C_Method,C_Strike,ARMLOCAL_SmileViewer_Create);

		/// call the general function
		fillXL_Result( LOCAL_SMILEVIEWER_CLASS, ourFunc, C_result, XL_result, PersistentInXL );
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_SmileViewer_CreateCommon" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


// Smile Viewer of a security
__declspec(dllexport) LPXLOPER WINAPI Local_SmileViewer_Create(
	LPXLOPER XL_SecurityId,
	LPXLOPER XL_ModelId,
	LPXLOPER XL_Moneyness,
	LPXLOPER XL_MoneyType,
	LPXLOPER XL_Strike)
{
	ADD_LOG("Local_SmileViewer_Create");
	bool PersistentInXL = true;
	return Local_SmileViewer_CreateCommon(
                                    XL_SecurityId,
									XL_ModelId,
									XL_Moneyness,
									XL_MoneyType,
									XL_Strike,
									PersistentInXL );
}

__declspec(dllexport) LPXLOPER WINAPI Local_PXL_SmileViewer_Create(
	LPXLOPER XL_SecurityId,
	LPXLOPER XL_ModelId,
	LPXLOPER XL_Moneyness,
	LPXLOPER XL_MoneyType,
	LPXLOPER XL_Strike)
{
	ADD_LOG("Local_PXL_SmileViewer_Create");
	bool PersistentInXL = false;
	return Local_SmileViewer_CreateCommon(
                                    XL_SecurityId,
									XL_ModelId,
									XL_Moneyness,
									XL_MoneyType,
									XL_Strike,
									PersistentInXL );
}

LPXLOPER Local_DensityFunctorGen_CreateCommon(
	LPXLOPER XL_SecurityId,
	LPXLOPER XL_ModelId,
	LPXLOPER XL_DecStrike,
	LPXLOPER XL_IsDirect,
    LPXLOPER XL_MinProba,
    LPXLOPER XL_MaxProba,
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
		
		CCString C_SecString;
		XL_readStrCell( XL_SecurityId, C_SecString, " ARM_ERR: Security Id: Object expected",C_result);
		long C_SecurityIdLong = LocalGetNumObjectId(C_SecString);

		CCString C_ModString;
		XL_readStrCell( XL_ModelId, C_ModString, " ARM_ERR: Model Id: Object expected",C_result);
		long C_ModelIdLong = LocalGetNumObjectId(C_ModString);
		
		double C_DecStrike;
		XL_readNumCell(XL_DecStrike,C_DecStrike," ARM_ERR: Dec Strike is not of a good type",C_result);

		CCString C_IsDirectStr;
		XL_readStrCellWD( XL_IsDirect,C_IsDirectStr,"N"," ARM_ERR: Is Direct: string expected",C_result);
		bool C_IsDirectBool;
		C_IsDirectStr.toUpper();
		if (C_IsDirectStr == "Y" )
			C_IsDirectBool = true;
		else if (C_IsDirectStr == "N")
			C_IsDirectBool = false;
		else
			throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"String \"Y\" or \"N\" Expected for Is Direct");

		double C_minProba;
		double C_minProbaDef=1e-5;
		XL_readNumCellWD(XL_MinProba,C_minProba,C_minProbaDef," ARM_ERR: Min Proba is not of a good type",C_result);

		double C_maxProba;
		double C_maxProbaDef=1.-1e-5;
		XL_readNumCellWD(XL_MaxProba,C_maxProba,C_maxProbaDef," ARM_ERR: Max Proba is not of a good type",C_result);

		/// use the concept of Functor to transfer the knowledge of
		/// a function with a context
		exportFunc6Args< long, long, double, bool, double, double > 
			ourFunc(C_SecurityIdLong,C_ModelIdLong,C_DecStrike,C_IsDirectBool,C_minProba,C_maxProba,ARMLOCAL_DensityFunctorGen_Create);

		/// call the general function
		fillXL_Result_withName(ourFunc, C_result, XL_result, PersistentInXL );
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_SmileViewer_CreateCommon" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


// Smile Viewer of a security
__declspec(dllexport) LPXLOPER WINAPI Local_DensityFunctorGen_Create(
	LPXLOPER XL_SecurityId,
	LPXLOPER XL_ModelId,
	LPXLOPER XL_DecStrike,
	LPXLOPER XL_IsDirect,
    LPXLOPER XL_MinProba,
    LPXLOPER XL_MaxProba)
{
	ADD_LOG("Local_DensityFunctorGen_Create");
	bool PersistentInXL = true;
	return Local_DensityFunctorGen_CreateCommon(
                                    XL_SecurityId,
									XL_ModelId,
									XL_DecStrike,
									XL_IsDirect,
									XL_MinProba,
									XL_MaxProba,
									PersistentInXL );
}

__declspec(dllexport) LPXLOPER WINAPI Local_PXL_DensityFunctorGen_Create(
	LPXLOPER XL_SecurityId,
	LPXLOPER XL_ModelId,
	LPXLOPER XL_DecStrike,
	LPXLOPER XL_IsDirect,
    LPXLOPER XL_MinProba,
    LPXLOPER XL_MaxProba)
{
	ADD_LOG("Local_PXL_DensityFunctorGen_Create");
	bool PersistentInXL = false;
	return Local_DensityFunctorGen_CreateCommon(
                                    XL_SecurityId,
									XL_ModelId,
									XL_DecStrike,
									XL_IsDirect,
									XL_MinProba,
									XL_MaxProba,
									PersistentInXL );
}
