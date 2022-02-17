/*
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 * $Log: ARM_xl_gp_nummethod_local.cpp,v $
 * Revision 1.1  2003/13/07 15:08:43  ebenhamou
 * Initial version
 *
 */

/*! \file ARM_xl_gp_nummethod_local.cpp
 *
 *  \brief file for the numerical method part in the generic pricer
 *
 *	\author  E Benhamou
 *	\version 1.0
 *	\date September 2003
 */

#include <ARM\libarm_local\firstToBeIncluded.h>
#include <functional>
#include <libCCxll\CCxll.h>

#include <ARM\libarm_local\ARM_local_gp_nummethod.h>
#include "ARM_xl_gp_nummethod_local.h"
#include "ARM_xl_wrapper_local.h"
#include "ARM_gp_local_interglob.h"
#include "ARM_xl_trycatch_local.h"
#include "ARM_xl_gp_fctorhelper.h"

#include <GP_Infra\gpinfra\nummethod.h>
#include <GP_Base\gpbase\eventviewerfwd.h>

#include "util\tech_macro.h"

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


////////////////////////////////////////////////
/// very rapid definition of ourlonglongFunctor
////////////////////////////////////////////////
struct BiNumFunc : ARMResultLong2LongFunc
{
	long operator()( ARM_result& result, long objId )
	{
		return ARMLOCAL_BINumMethod_Create( result, objId );			
	}
};

////////////////////////////////////////////////
/// very rapid definition of ourlonglongFunctor
////////////////////////////////////////////////
struct FiNumFunc : ARMResultLong2LongFunc
{
	long operator()( ARM_result& result, long objId )
	{
		return ARMLOCAL_FINumMethod_Create( result, objId );			
	}
};

////////////////////////////////////////////////
/// very rapid definition of ourlonglongFunctor
////////////////////////////////////////////////
struct MixteNumFunc : ARMResultLong2LongFunc
{
	long operator()( ARM_result& result, long objId )
	{
		return ARMLOCAL_MixteNumMethod_Create( result, objId );
	}
};


/////////////////////////////////////////////////////////////
/// central function that the creation of a backward induct
///	method
/////////////////////////////////////////////////////////////
template<typename ProtoNumMethod >
	LPXLOPER Local_NumMethodCommon(
		LPXLOPER XL_StepsNb,
		LPXLOPER XL_TruncationPolicy,
		bool PersistentInXL,
		const CCString& LOCAL_OBJECT_CLASS,
		ProtoNumMethod& fctor )
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

		/// call the general function
		fillXL_Result( LOCAL_OBJECT_CLASS, fctor, C_result, XL_result, PersistentInXL );
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_NumMethodCommon" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}
							 

///////////////////////////////////
/// version that takes into account 
/// previous creation of object
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_BINumMethod_Create(
	LPXLOPER XL_StepsNb,
	LPXLOPER XL_TruncationPolicy )
{
	ADD_LOG("Local_BINumMethod_Create");
	bool PersistentInXL = true;
	return Local_NumMethodCommon( XL_StepsNb, XL_TruncationPolicy, PersistentInXL, LOCAL_BINUMMETHOD_CLASS, BiNumFunc() );
}


///////////////////////////////////
/// version for non persistentInXl
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_PXL_BINumMethod_Create(
	LPXLOPER XL_StepsNb,
	LPXLOPER XL_TruncationPolicy )
{
	ADD_LOG("Local_PXL_BINumMethod_Create");
	bool PersistentInXL = false;
	return Local_NumMethodCommon( XL_StepsNb, XL_TruncationPolicy, PersistentInXL, LOCAL_BINUMMETHOD_CLASS, BiNumFunc() );
}


///////////////////////////////////
/// version that takes into account 
/// previous creation of object
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_FINumMethod_Create(
	LPXLOPER XL_StepsNb,
	LPXLOPER XL_TruncationPolicy )
{
	ADD_LOG("Local_FINumMethod_Create");
	bool PersistentInXL = true;
	return Local_NumMethodCommon( XL_StepsNb, XL_TruncationPolicy, PersistentInXL, LOCAL_FINUMMETHOD_CLASS, FiNumFunc() );
}


///////////////////////////////////
/// version for non persistentInXl
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_PXL_FINumMethod_Create(
	LPXLOPER XL_StepsNb,
	LPXLOPER XL_TruncationPolicy )
{
	ADD_LOG("Local_PXL_FINumMethod_Create");
	bool PersistentInXL = false;
	return Local_NumMethodCommon( XL_StepsNb, XL_TruncationPolicy, PersistentInXL, LOCAL_MIXTENUMMETHOD_CLASS, FiNumFunc() );
}


///////////////////////////////////
/// version that takes into account 
/// previous creation of object
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_MixteNumMethod_Create(
	LPXLOPER XL_StepsNb,
	LPXLOPER XL_TruncationPolicy )
{
	ADD_LOG("Local_MixteNumMethod_Create");
	bool PersistentInXL = true;
	return Local_NumMethodCommon( XL_StepsNb, XL_TruncationPolicy, PersistentInXL, LOCAL_MIXTENUMMETHOD_CLASS, MixteNumFunc() );
}


///////////////////////////////////
/// version for non persistentInXl
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_PXL_MixteNumMethod_Create(
	LPXLOPER XL_StepsNb,
	LPXLOPER XL_TruncationPolicy )
{
	ADD_LOG("Local_PXL_MixteNumMethod_Create");
	bool PersistentInXL = false;
	return Local_NumMethodCommon( XL_StepsNb, XL_TruncationPolicy, PersistentInXL, LOCAL_FINUMMETHOD_CLASS, MixteNumFunc() );
}


////////////////////////////////////////////////
/// very rapid definition of ourlonglongFunctor
////////////////////////////////////////////////
class SetObj2ToModelFunc : public ARMResultLong2LongFunc
{
public:
	typedef long (*Function)(long,long,ARM_result&, long);

	SetObj2ToModelFunc(long modelId, long obj2Id, Function f)
	:	C_modelId(modelId), C_Obj2Id(obj2Id), itsFunction(f)  {};

	long operator()( ARM_result& result, long objId )
	{
		return (*itsFunction)(C_modelId, C_Obj2Id, result, objId );			
	}

private:
	long C_modelId;
	long C_Obj2Id;
	Function itsFunction;

};


/////////////////////////////////////////////////////////////
/// central function to set a numerical method to a model
/////////////////////////////////////////////////////////////
LPXLOPER Local_SetSomethingToModelCommon(
	LPXLOPER XL_modelId,
	LPXLOPER XL_obj2Id,
	const string& Obj2Name,
	long (*Function)(long,long,ARM_result&, long ),
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
		XL_readStrCell( XL_modelId, C_modelStrId,	" ARM_ERR: Model Id: Object expected",		C_result);
		long C_modelId = LocalGetNumObjectId(C_modelStrId);

		long C_obj2Id;
		string errorMsg = " ARM_ERR: " + Obj2Name + " Id: Object expected";
		XL_GETOBJID( XL_obj2Id, C_obj2Id, errorMsg.c_str(),	C_result);
		
		/// use the concept of Functor to transfer the knowledge of
		/// a function with a context
		SetObj2ToModelFunc  ourFunc(C_modelId,C_obj2Id,Function);
		CCString stringId = GetLastCurCellEnvValue ();
		CCString curClass = LocalGetStringObjectClassAndError(C_modelStrId);

		/// call the general function
		if( curClass != CCString(" ") )
			fillXL_Result( curClass, ourFunc, C_result, XL_result, PersistentInXL );
		/// bracket are necessary as ARM_ERR_AND_EXIT is a macro
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_SetSomethingToModelCommon" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}
							 


///////////////////////////////////
/// version that takes into account 
/// previous creation of object
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_SetNumMethodtoModel(
	LPXLOPER XL_modelId,
	LPXLOPER XL_numMethodId )
{
	ADD_LOG("Local_SetNumMethodtoModel");
	bool PersistentInXL = true;
	return Local_SetSomethingToModelCommon( XL_modelId, XL_numMethodId, 
		"Numerical Method", ARMLOCAL_SetNumMethod, PersistentInXL );
}


///////////////////////////////////
/// version for non persistentInXl
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_PXL_SetNumMethodtoModel(
	LPXLOPER XL_modelId,
	LPXLOPER XL_numMethodId )
{
	ADD_LOG("Local_PXL_SetNumMethodtoModel");
	bool PersistentInXL = false;
	return Local_SetSomethingToModelCommon( XL_modelId, XL_numMethodId, 
		"Numerical Method", ARMLOCAL_SetNumMethod, PersistentInXL );
}



///////////////////////////////////
/// version that takes into account 
/// previous creation of object
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_SetNumerairetoModel(
	LPXLOPER XL_modelId,
	LPXLOPER XL_numeraireId )
{
	ADD_LOG("Local_SetNumerairetoModel");
	bool PersistentInXL = true;
	return Local_SetSomethingToModelCommon( XL_modelId, XL_numeraireId, 
		"Numeraire", ARMLOCAL_SetNumeraire, PersistentInXL );
}


///////////////////////////////////
/// version for non persistentInXl
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_PXL_SetNumerairetoModel(
	LPXLOPER XL_modelId,
	LPXLOPER XL_numeraireId )
{
	ADD_LOG("Local_PXL_SetNumerairetoModel");
	bool PersistentInXL = false;
	return Local_SetSomethingToModelCommon( XL_modelId, XL_numeraireId, 
		"Numeraire", ARMLOCAL_SetNumeraire, PersistentInXL );
}


/////////////////////////////////////////////////////////////
//// Get the numerical method of a GP model
/////////////////////////////////////////////////////////////
class getNumMethodFunc : public ARMResultLong2LongFunc
{
public:
	getNumMethodFunc(long modelId) : C_modelId(modelId) {};
	
	long operator()( ARM_result& result, long objId ){
		return ARMLOCAL_GetNumMethodFromModel(C_modelId,result,objId);
    }

private:
	long    C_modelId;
};

LPXLOPER Local_GetNumMethodFromModelCommon(
	LPXLOPER XL_modelId,
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

		long modelId;
		XL_GETOBJID(XL_modelId,modelId," ARM_ERR: Model id: object expected",C_result);

        /// Get the type to return
		CCString numMethodClass(ModelNumMethodToClass(modelId).c_str());
		
		/// use the concept of Functor to transfer the knowledge of
		/// a function with a context
		getNumMethodFunc  ourFunc(modelId);

		fillXL_Result( numMethodClass, ourFunc, C_result, XL_result, PersistentInXL );
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_GetNumMethodFromModelCommon" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}
							 
__declspec(dllexport) LPXLOPER WINAPI Local_GetNumMethodFromModel(
	LPXLOPER XL_modelId)
{
	ADD_LOG("Local_GetNumMethodFromModel");
	bool PersistentInXL = true;
	return Local_GetNumMethodFromModelCommon( XL_modelId,PersistentInXL );
}

__declspec(dllexport) LPXLOPER WINAPI Local_PXL_GetNumMethodFromModel(
	LPXLOPER XL_modelId)
{
	ADD_LOG("Local_PXL_GetNumMethodFromModel");
	bool PersistentInXL = false;
	return Local_GetNumMethodFromModelCommon( XL_modelId,PersistentInXL );
}


/////////////////////////////////////////////////////////////
//// Set proba computation flag to a GP tree
/////////////////////////////////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_Tree_SetProbaFlag(
	LPXLOPER XL_treeId,
	LPXLOPER XL_probaFlag)
{
	ADD_LOG("Local_Tree_SetProbaFlag");
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
		long treeId;
		XL_GETOBJID(XL_treeId,treeId," ARM_ERR: Tree id: object expected",C_result);

		double probaFlag;
		XL_readNumCell(XL_probaFlag, probaFlag,	" ARM_ERR: Proba Computation Flag: boolean compatible expected",	C_result);
		bool isSpotProba = probaFlag != 0;
			
		long retCode = ARMLOCAL_Local_Tree_SetProbaFlag(
				treeId,
				isSpotProba,
				C_result );

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
            ARM_ARG_ERR();
        }
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_Tree_SetProbaFlag" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


///----------------------------------------------
///----------------------------------------------
///             Numeraire
/// Inputs :
///     Type
///     Time lags vector
///----------------------------------------------
///----------------------------------------------
////////////////////////////////////////////////
/// very rapid definition of ourlonglongFunctor
////////////////////////////////////////////////
class NumeraireFunc : public ARMResultLong2LongFunc
{
public:
	NumeraireFunc(
        long numeraireType,
        const VECTOR<double>&  numeraireTimes)
    :
    C_NumeraireType(numeraireType),
    C_NumeraireTimes(numeraireTimes)
    {};
	
	long operator()( ARM_result& result, long objId ){
		return ARMLOCAL_Numeraire_Create(
            C_NumeraireType,
            C_NumeraireTimes,
            result, objId);			
	}

private:
    long            C_NumeraireType;
    VECTOR<double>  C_NumeraireTimes;
};


/////////////////////////////////////////////////////////////
/// central function that does the creation of the XL function
/////////////////////////////////////////////////////////////
LPXLOPER Local_NumeraireCommon(
	LPXLOPER XL_NumeraireType,
	LPXLOPER XL_NumeraireTimes,
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

		CCString C_NumeraireTypeStr;
		XL_readStrCell(XL_NumeraireType,C_NumeraireTypeStr," ARM_ERR: Numeraire Type: String expected",C_result);
		long C_NumeraireType;
		if( (C_NumeraireType = ARM_ConvGPNumeraire( C_NumeraireTypeStr, C_result)) == ARM_DEFAULT_ERR )
		{
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}

		VECTOR<double> C_NumeraireTimes;
		VECTOR<double> C_defautNumeraireTimes;
		XL_readNumVectorWD(XL_NumeraireTimes,C_NumeraireTimes,C_defautNumeraireTimes," ARM_ERR: Numeraire Times: array of numeric expected",C_result);


		/// use the concept of Functor to transfer the knowledge of
		/// a function with a context
		NumeraireFunc ourFunc(
			C_NumeraireType,
			C_NumeraireTimes);

		/// call the general function
		fillXL_Result( LOCAL_NUMERAIRE_CLASS, ourFunc, C_result, XL_result, PersistentInXL );
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_NumeraireCommon" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}
							 

///////////////////////////////////
/// version that takes into account 
/// previous creation of object
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_Numeraire_Create(
	LPXLOPER XL_NumeraireType,
	LPXLOPER XL_NumeraireTimes)
{
	ADD_LOG("Local_Numeraire_Create");
	bool PersistentInXL = true;
	return Local_NumeraireCommon(
	    XL_NumeraireType,
	    XL_NumeraireTimes,
        PersistentInXL );
}


///////////////////////////////////
/// version for non persistentInXl
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_PXL_Numeraire_Create(
	LPXLOPER XL_NumeraireType,
	LPXLOPER XL_NumeraireTimes)
{
	ADD_LOG("Local_PXL_Numeraire_Create");
	bool PersistentInXL = false;
	return Local_NumeraireCommon(
	    XL_NumeraireType,
	    XL_NumeraireTimes,
        PersistentInXL );
}



///----------------------------------------------
///----------------------------------------------
///             Tree Method
/// Inputs :
///     Number of steps
///     Truncation level (optional)
///----------------------------------------------
///----------------------------------------------
////////////////////////////////////////////////
/// very rapid definition of ourlonglongFunctor
////////////////////////////////////////////////
class TreeMethodFunc : public ARMResultLong2LongFunc
{
public:
	TreeMethodFunc(
        int     nbSteps,
        double  nbStdDev,
        double  minStdDev,
        int     nbMinSteps,
        bool    isTree1GForced)
	:
    C_NbSteps(nbSteps),
    C_NbStdDev(nbStdDev),
    C_MinStdDev(minStdDev),
    C_NbMinSteps(nbMinSteps),
    C_IsTree1GForced(isTree1GForced)
    {};
	
	long operator()( ARM_result& result, long objId ){
		return ARMLOCAL_TreeMethod_Create(
            C_NbSteps,
            C_NbStdDev,
            C_MinStdDev,
            C_NbMinSteps,
            C_IsTree1GForced,
            result,
            objId);			
	}

private:
	int         C_NbSteps;
	double      C_NbStdDev;
    double      C_MinStdDev;
    int         C_NbMinSteps;
    bool        C_IsTree1GForced;
};


/////////////////////////////////////////////////////////////
/// central function that does the creation of the XL function
/////////////////////////////////////////////////////////////
LPXLOPER Local_TreeMethodCommon(
	LPXLOPER XL_NbSteps,
	LPXLOPER XL_NbStdDev,
	LPXLOPER XL_MinStdDev,
	LPXLOPER XL_NbMinSteps,
    bool IsTree1GForced,
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

		double C_NbSteps;
		XL_readNumCell(XL_NbSteps,C_NbSteps," ARM_ERR: nb of steps: int expected",C_result);

		double C_NbStdDev;
		double defaultNbStdDev = 5.0;
		XL_readNumCellWD(XL_NbStdDev,C_NbStdDev,defaultNbStdDev," ARM_ERR: nb of global StdDev: double expected",C_result);

		double C_MinStdDev;
		double defaultNbMinStdDev = 1.0e-3;
		XL_readNumCellWD(XL_MinStdDev,C_MinStdDev,defaultNbMinStdDev," ARM_ERR: min local StdDev without space step resizing : double expected",C_result);

		double C_NbMinSteps;
		double defaultNbMinSteps = 0;
		XL_readNumCellWD(XL_NbMinSteps,C_NbMinSteps,defaultNbMinSteps," ARM_ERR: min nb of steps before 1st event date : int expected",C_result);

		/// use the concept of Functor to transfer the knowledge of
		/// a function with a context
        IsTree1GForced = IsTree1GForced || !ARM::isTree2G;
		TreeMethodFunc ourFunc(
			(int)(floor(C_NbSteps)),
			C_NbStdDev,
			C_MinStdDev,
			C_NbMinSteps,
            IsTree1GForced);

		/// call the general function
        fillXL_Result( (IsTree1GForced ? LOCAL_TREEMETHOD_CLASS : LOCAL_TREEND_CLASS), ourFunc, C_result, XL_result, PersistentInXL );
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_TreeMethodCommon" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}
							 

///////////////////////////////////
/// version that takes into account 
/// previous creation of object
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_TreeMethod_Create(
    LPXLOPER XL_NbSteps,
	LPXLOPER XL_NbStdDev,
    LPXLOPER XL_MinStdDev,
    LPXLOPER XL_NbMinSteps)
{
	ADD_LOG("Local_TreeMethod_Create");
	bool PersistentInXL = true;
	return Local_TreeMethodCommon(
        XL_NbSteps,
        XL_NbStdDev,
        XL_MinStdDev,
        XL_NbMinSteps,
        false,
        PersistentInXL );
}


///////////////////////////////////
/// version for non persistentInXl
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_PXL_TreeMethod_Create(
    LPXLOPER XL_NbSteps,
	LPXLOPER XL_NbStdDev,
    LPXLOPER XL_MinStdDev,
    LPXLOPER XL_NbMinSteps)
{
	ADD_LOG("Local_PXL_TreeMethod_Create");
	bool PersistentInXL = false;
	return Local_TreeMethodCommon(
        XL_NbSteps,
        XL_NbStdDev,
        XL_MinStdDev,
        XL_NbMinSteps,
        false,
        PersistentInXL );
}



///----------------------------------------------
///----------------------------------------------
///             MC Method
/// Inputs :
///     Number of iteration
///		Fix time steps
///----------------------------------------------
///----------------------------------------------
////////////////////////////////////////////////
/// very rapid definition of ourlonglongFunctor
////////////////////////////////////////////////


/////////////////////////////////////////////////////////////
/// central function that does the creation of the XL function
/////////////////////////////////////////////////////////////
LPXLOPER Local_MCMethod_Common(
	LPXLOPER XL_ItersNb,
	LPXLOPER XL_FixStep,
	LPXLOPER XL_RandGenIds,
	LPXLOPER XL_SamplerType,
	LPXLOPER XL_SamplerDatas,
	LPXLOPER XL_SchedulerType,
	LPXLOPER XL_SchedulerDatas,
	LPXLOPER XL_ExercBoundCalcId,
	LPXLOPER XL_MaxBucketSize,
	LPXLOPER XL_ImpSamplerType,
	LPXLOPER XL_ImpSamplerDatas,
	LPXLOPER XL_PathSchemeType,
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

		double C_ItersNb;
		size_t C_ItersNbInt;
		double defaultItersNb = 2;			/// 2 nb of iteration
		XL_readNumCellWD(XL_ItersNb,C_ItersNb,defaultItersNb, " ARM_ERR: iters nb : double expected",C_result);
		C_ItersNbInt = C_ItersNb;

		double C_FixStep;
		int C_FixStepInt;
		double defaultFixStep = 1;			///  1 days
		XL_readNumCellWD(XL_FixStep,C_FixStep,defaultFixStep, " ARM_ERR: fix steps: double expected", C_result);
		C_FixStepInt = C_FixStep;

		   
		VECTOR<CCString> C_RandGens;
		VECTOR<CCString> C_RandGensDef(0);
		XL_readStrVectorWD (XL_RandGenIds,C_RandGens,C_RandGensDef," ARM_ERR: RandGen: array of object expected",DOUBLE_TYPE,C_result);
		size_t size = C_RandGens.size();
		VECTOR<long> C_RandGensId;
		C_RandGensId.resize(size);
		size_t i;
		for(i = 0; i < size; ++i )    
			C_RandGensId[i] = LocalGetNumObjectId(C_RandGens[i]);

		CCString C_SchedulerTypeStr;
		CCString defaultSchedulerTypeStr = "TimeStepPerYearScheduler";
		XL_readStrCellWD(XL_SchedulerType,C_SchedulerTypeStr,defaultSchedulerTypeStr," ARM_ERR: Scheduler type: String expected",C_result);
		long C_SchedulerType;
		if( (C_SchedulerType = ARM_ConvGPSchedulerType( C_SchedulerTypeStr, C_result)) == ARM_DEFAULT_ERR )
		{
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}

		VECTOR<double> C_SchedulerDatas;
		VECTOR<double> defaultSchedulerDatas(1);
		defaultSchedulerDatas[0] = C_FixStep;
		XL_readNumVectorWD(XL_SchedulerDatas,C_SchedulerDatas,defaultSchedulerDatas," ARM_ERR: Scheduler datas: array of numeric expected",C_result);


		CCString C_SamplerTypeStr;
		CCString defaultSamplerTypeStr = "NormalCentredSampler";
		XL_readStrCellWD(XL_SamplerType,C_SamplerTypeStr,defaultSamplerTypeStr," ARM_ERR: Sampler type: String expected",C_result);
		long C_SamplerType;
		if( (C_SamplerType = ARM_ConvGPSamplerType( C_SamplerTypeStr, C_result)) == ARM_DEFAULT_ERR )
		{
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}

		VECTOR<double> C_SamplerDatas;
		VECTOR<double> defaultSamplerDatas(0);
		XL_readNumVectorWD(XL_SamplerDatas,C_SamplerDatas,defaultSamplerDatas," ARM_ERR: Scheduler datas: array of numeric expected",C_result);

		CCString C_ImpSamplerTypeStr;
		CCString defaultImpSamplerTypeStr = "Dummy";
		XL_readStrCellWD(XL_ImpSamplerType,C_ImpSamplerTypeStr,defaultImpSamplerTypeStr," ARM_ERR: Importance Sampler type: String expected",C_result);
		long C_ImpSamplerType;
		if( (C_ImpSamplerType = ARM_ConvGPImpSamplerType( C_ImpSamplerTypeStr, C_result)) == ARM_DEFAULT_ERR )
		{
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}

		VECTOR<CCString> C_ImpSamplerStrDatas;
		VECTOR<double> C_ImpSamplerDatas;
		VECTOR<CCString> defaultImpSamplerStrDatas(0);
		XL_readStrVectorWD(XL_ImpSamplerDatas,C_ImpSamplerStrDatas,defaultImpSamplerStrDatas," ARM_ERR: Importance Sampler datas: array of id expected",DOUBLE_TYPE,C_result);

		size = C_ImpSamplerStrDatas.size();
		C_ImpSamplerDatas.resize(size);
		for(i = 0; i < size; ++i )    
			C_ImpSamplerDatas[i] = LocalGetNumObjectId(C_ImpSamplerStrDatas[i]);

		CCString C_PathSchemeTypeStr;
		CCString defaultPathSchemeTypeStr = "Incremental";
		XL_readStrCellWD(XL_PathSchemeType,C_PathSchemeTypeStr,defaultPathSchemeTypeStr," ARM_ERR: Path Scheme type: String expected",C_result);
		long C_PathSchemeType;
		if( (C_PathSchemeType = ARM_ConvGPPathSchemeType( C_PathSchemeTypeStr, C_result)) == ARM_DEFAULT_ERR )
		{
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}

		long C_ExercBoundCalcId;
		XL_GETOBJIDWD( XL_ExercBoundCalcId, C_ExercBoundCalcId,	"NULL OBJECT",	" ARM_ERR: Exercise Boundary Calculator: object expected",	C_result);

		double MaxBucketSize;
		size_t C_MaxBucketSize;
		double defaultBucketSize = 100000;			/// 100000 nb of iteration per bucket
		XL_readNumCellWD(XL_MaxBucketSize,MaxBucketSize,defaultBucketSize, " ARM_ERR: Max Bucket Size : double expected",C_result);

		C_MaxBucketSize = MaxBucketSize;

		/// use the concept of Functor to transfer the knowledge of
		/// a function with a context
		exportFunc12Args<size_t,int,VECTOR<long>,long,VECTOR<double>,long,VECTOR<double>,long,VECTOR<double>,long,long,size_t> ourFunc(
			C_ItersNbInt,
			C_FixStepInt,
			C_RandGensId,
			C_SamplerType,
			C_SamplerDatas,
			C_SchedulerType,
			C_SchedulerDatas,
			C_ImpSamplerType,
			C_ImpSamplerDatas,
			C_PathSchemeType,
			C_ExercBoundCalcId,
			C_MaxBucketSize,
			ARMLOCAL_MCMethod_Create);

		/// call the general function
		fillXL_Result( LOCAL_MCMETHOD_CLASS, ourFunc, C_result, XL_result, PersistentInXL );
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_MCMethod_Common" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}
							 

///////////////////////////////////
/// version that takes into account 
/// previous creation of object
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_MCMethod_Create(
	LPXLOPER XL_ItersNb,
	LPXLOPER XL_FixStep,
	LPXLOPER XL_RandGenIds,
	LPXLOPER XL_SamplerType,
	LPXLOPER XL_SamplerDatas,
	LPXLOPER XL_SchedulerType,
	LPXLOPER XL_SchedulerDatas,
	LPXLOPER XL_ExercBoundCalcId, 
	LPXLOPER XL_MaxBucketSize,
	LPXLOPER XL_ImpSamplerType,
	LPXLOPER XL_ImpSamplerDatas,
	LPXLOPER XL_PathSchemeType)
{
	ADD_LOG("Local_MCMethod_Create");
	bool PersistentInXL = true;
	return Local_MCMethod_Common(
        XL_ItersNb,
        XL_FixStep,
		XL_RandGenIds,
		XL_SamplerType,
		XL_SamplerDatas,
		XL_SchedulerType,
		XL_SchedulerDatas,
		XL_ExercBoundCalcId,
		XL_MaxBucketSize,
		XL_ImpSamplerType,
		XL_ImpSamplerDatas,
		XL_PathSchemeType,
		PersistentInXL );
}


///////////////////////////////////
/// version for non persistentInXl
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_PXL_MCMethod_Create(
	LPXLOPER XL_ItersNb,
	LPXLOPER XL_FixStep,
	LPXLOPER XL_RandGenIds,
	LPXLOPER XL_SamplerType,
	LPXLOPER XL_SamplerDatas,
	LPXLOPER XL_SchedulerType,
	LPXLOPER XL_SchedulerDatas,
	LPXLOPER XL_ExercBoundCalcId, 
	LPXLOPER XL_MaxBucketSize,
	LPXLOPER XL_ImpSamplerType,
	LPXLOPER XL_ImpSamplerDatas,
	LPXLOPER XL_PathSchemeType)
{
	ADD_LOG("Local_PXL_MCMethod_Create");
	bool PersistentInXL = false;
	return Local_MCMethod_Common(
        XL_ItersNb,
        XL_FixStep,
		XL_RandGenIds,
		XL_SamplerType,
		XL_SamplerDatas,
		XL_SchedulerType,
		XL_SchedulerDatas,
		XL_ExercBoundCalcId,
		XL_MaxBucketSize,
		XL_ImpSamplerType,
		XL_ImpSamplerDatas,
		XL_PathSchemeType,
		PersistentInXL );
}

/////////////////////////////////////////////////////////////
/// central function that does the creation of the XL function
/////////////////////////////////////////////////////////////
LPXLOPER Local_PDEMethod_Common( 
		LPXLOPER XL_MethodName,
        LPXLOPER XL_SchedulerData,
		LPXLOPER XL_SpaceDiscretization,
		LPXLOPER XL_GridData,
		LPXLOPER XL_NY,
		LPXLOPER XL_NZ,
		LPXLOPER XL_Theta1,
		LPXLOPER XL_Theta2,
		LPXLOPER XL_Theta3,
		LPXLOPER XL_BoundaryLimitName,
		LPXLOPER XL_Lambda,
		LPXLOPER XL_GridType,
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

		/// FIX FIX discounting curve ignored at this stage!
		CCString MethodNameXL;
		XL_readStrCellWD(XL_MethodName,MethodNameXL,""," ARM_ERR: discount curve name: string expected",C_result);
        string MethodName = CCSTringToSTLString(MethodNameXL);

		VECTOR<double> C_SchedulerData;
		long C_SchedulerDataRows, C_SchedulerDataCols;
		XL_readNumVectorAndSize( XL_SchedulerData, C_SchedulerDataRows, C_SchedulerDataCols, C_SchedulerData, " ARM_ERR: scheduler data vector or matrix expected", C_result);
		int func_SchedulerDataRows=C_SchedulerDataRows;
		int func_SchedulerDataCols=C_SchedulerDataCols;

		double C_SpaceDiscretization;
		double defaultSpaceDiscretization = 101;			
		XL_readNumCellWD(XL_SpaceDiscretization,C_SpaceDiscretization,defaultSpaceDiscretization, " ARM_ERR: Space Discretisation (for all the scheme excepts CraigSneyd) Nb : double expected",C_result);
		int func_SpaceDiscretization=C_SpaceDiscretization;

		CCString GridTypeXL;
		XL_readStrCellWD(XL_GridType,GridTypeXL,"StdDev"," ARM_ERR: grid type: string expected",C_result);
        string C_GridType = CCSTringToSTLString(GridTypeXL);

		VECTOR<double> C_GridData;
		long C_GridDataRows, C_GridDataCols;
		XL_readNumVectorAndSize( XL_GridData, C_GridDataRows, C_GridDataCols, C_GridData, " ARM_ERR: grid data vector or matrix expected", C_result);
		int func_GridDataRows=C_GridDataRows;
		int func_GridDataCols=C_GridDataCols;
	
		double C_NY;
		double defaultNY = 3;		
		XL_readNumCellWD(XL_NY,C_NY,defaultNY, " ARM_ERR: Y Grid Nb : double expected",C_result);
		int func_NY=C_NY;


		double C_NZ;
		double defaultNZ = 3;		
		XL_readNumCellWD(XL_NZ,C_NZ,defaultNZ, " ARM_ERR: Z Grid Nb : double expected",C_result);
		int func_NZ=C_NZ;


		double C_Theta1;
		double defaultTheta1 = 0.5;			
		XL_readNumCellWD(XL_Theta1,C_Theta1,defaultTheta1, " ARM_ERR: Theta1 : double expected",C_result);

		double C_Theta2;
		double defaultTheta2 = 0.5;			
		XL_readNumCellWD(XL_Theta2,C_Theta2,defaultTheta2, " ARM_ERR: Theta2 : double expected",C_result);

		double C_Theta3;
		double defaultTheta3 = 0.5;		
		XL_readNumCellWD(XL_Theta3,C_Theta3,defaultTheta3, " ARM_ERR: Theta3 : double expected",C_result);

		CCString BoundaryLimitNameXL;
		XL_readStrCellWD(XL_BoundaryLimitName,BoundaryLimitNameXL,"VonNeumann"," ARM_ERR: Boundary condition name: string expected",C_result);
        string BoundaryConditionName= CCSTringToSTLString(BoundaryLimitNameXL);

		double C_Lambda;
		double defaultLambda = 0.0;		
		XL_readNumCellWD(XL_Lambda,C_Lambda,defaultLambda, " ARM_ERR: Lambda : double expected",C_result);

		/// use the concept of Functor to transfer the knowledge of
		/// a function with a context
		exportFunc16Args<
			string,		// 1 Method Name
			vector<double>,	// 2 Scheduler data
			int,		// 3 SchedulerNbRows
			int,		// 4 SchedulerNbCols
			int,		// 5 Space Iter Nb  
			string,		// 6 Grid Type
			vector<double>,// 7 Grid Data
			int,		// 8 GridNbRows
			int,		// 9 GridNbCols
			int,		// 10 YGridIterNb
			int,		// 11 ZGridIterNb
			double,		// 12 Theta1
			double,		// 13 Theta2
			double,		// 14 Theta3
			string,		// 15 Bound Cond
			double>		// 16 Lambda
			ourFunc(
				MethodName,						// Method Name
				C_SchedulerData,		// Scheduler Data Nb
				func_SchedulerDataRows,		// SchedulerNbRows
				func_SchedulerDataCols,		// SchedulerNbCols	
				func_SpaceDiscretization, // Space Iter Nb 
				C_GridType,						// Grid Type
				C_GridData,						// Grid Data
				func_GridDataRows,		// GridNbRows
				func_GridDataCols,		// GridNbCols	
				func_NY,					// YGridIterNb
				func_NZ,					// ZGridIterNb
				 C_Theta1,						// Theta1
				 C_Theta2,						// Theta2
				 C_Theta3,						// Theta3
				 BoundaryConditionName,			// Bound Cond
				 C_Lambda,						// Lambda
				 ARMLOCAL_PDEMethod_Create );

		/// call the general function
		fillXL_Result_withName( ourFunc, C_result, XL_result, PersistentInXL );
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_MCMethod_Common" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}
							 

///////////////////////////////////
/// version that takes into account 
/// previous creation of object
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_PDEMethod_Create(
		LPXLOPER XL_MethodName,
		LPXLOPER XL_SchedulerData,
		LPXLOPER XL_SpaceDiscretization,
		LPXLOPER XL_GridData,
		LPXLOPER XL_NY,
		LPXLOPER XL_NZ,
		LPXLOPER XL_Theta1,
		LPXLOPER XL_Theta2,
		LPXLOPER XL_Theta3,
		LPXLOPER XL_BoundaryLimitName,
		LPXLOPER XL_Lambda,
		LPXLOPER XL_GridType)
{
	ADD_LOG("Local_PDEMethod_Create");
	bool PersistentInXL = true;
	return Local_PDEMethod_Common(
		XL_MethodName,
		XL_SchedulerData,
		XL_SpaceDiscretization,
		XL_GridData,
		XL_NY,
		XL_NZ,
		XL_Theta1,
		XL_Theta2,
		XL_Theta3,
		XL_BoundaryLimitName,
		XL_Lambda,
		XL_GridType,
		PersistentInXL);
}


///////////////////////////////////
/// version for non persistentInXl
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_PXL_PDEMethod_Create(
		LPXLOPER XL_MethodName,
        LPXLOPER XL_AddTimeStepsNb,
		LPXLOPER XL_SpaceDiscretization,
		LPXLOPER XL_GridData,
		LPXLOPER XL_NY,
		LPXLOPER XL_NZ,
		LPXLOPER XL_Theta1,
		LPXLOPER XL_Theta2,
		LPXLOPER XL_Theta3,
		LPXLOPER XL_BoundaryLimitName,
		LPXLOPER XL_Lambda,
		LPXLOPER XL_GridType)
{
	ADD_LOG("Local_PXL_PDEMethod_Create");
	bool PersistentInXL = false;
	return Local_PDEMethod_Common(
		XL_MethodName,
		XL_AddTimeStepsNb,
		XL_SpaceDiscretization,
		XL_GridData,
		XL_NY,
		XL_NZ,
		XL_Theta1,
		XL_Theta2,
		XL_Theta3,
		XL_BoundaryLimitName,
		XL_Lambda,
		XL_GridType,
		PersistentInXL);
}





/////////////////////////////////////////////////////////////
/// central function that does the creation of the XL function
/////////////////////////////////////////////////////////////
LPXLOPER Local_PdeND_Common( 
		LPXLOPER XL_MethodName,
        LPXLOPER XL_SchedulerType,
        LPXLOPER XL_SchedulerDataId,
		LPXLOPER XL_SpaceDataId,
		LPXLOPER XL_SchemeDataId,
		LPXLOPER XL_BoundCondName,
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

		CCString MethodNameXL;
		XL_readStrCellWD(XL_MethodName,MethodNameXL,""," ARM_ERR: method name: string expected",C_result);
        string C_MethodName = CCSTringToSTLString(MethodNameXL);

		CCString C_SchedulerTypeStr;
		CCString defaultSchedulerTypeStr = "ConstantVarianceMeanRevertingScheduler";
		XL_readStrCellWD(XL_SchedulerType,C_SchedulerTypeStr,defaultSchedulerTypeStr," ARM_ERR: Scheduler type: String expected",C_result);
		long C_SchedulerType;
		if( (C_SchedulerType = ARM_ConvGPSchedulerType( C_SchedulerTypeStr, C_result)) == ARM_DEFAULT_ERR )
		{
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}

		CCString C_SchedulerDataStrId;
		XL_readStrCell( XL_SchedulerDataId, C_SchedulerDataStrId," ARM_ERR: SchedulerData Id: Object expected",C_result);
		long C_SchedulerDataId = LocalGetNumObjectId(C_SchedulerDataStrId);

		CCString C_SpaceDataStrId;
		XL_readStrCell( XL_SpaceDataId, C_SpaceDataStrId,	" ARM_ERR: SpaceData Id: Object expected", C_result);
		long C_SpaceDataId = LocalGetNumObjectId(C_SpaceDataStrId);

		CCString C_SchemeDataStrId;
		XL_readStrCell( XL_SchemeDataId, C_SchemeDataStrId,	" ARM_ERR: SchemeData Id: Object expected",	C_result);
		long C_SchemeDataId = LocalGetNumObjectId(C_SchemeDataStrId);

		CCString BoundCondNameXL;
		XL_readStrCellWD(XL_BoundCondName,BoundCondNameXL,"VonNeumann"," ARM_ERR: Boundary condition name: string expected",C_result);
        string C_BoundCondName= CCSTringToSTLString(BoundCondNameXL);

		/// a function with a context
		//exportFunc6Args< string, long, long, long, long, string > ourFunc(C_MethodName, C_SchedulerType, C_SchedulerDataId, C_SpaceDataId, C_SchemeDataId, C_BoundCondName, ARMLOCAL_PdeND_Create );//pour le moment

		/// call the general function
		//fillXL_Result_withName( ourFunc, C_result, XL_result, PersistentInXL );//pour le moment
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_MCMethod_Common" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}
							 

///////////////////////////////////
/// version that takes into account 
/// previous creation of object
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_PdeND_Create(
		LPXLOPER XL_MethodName,
		LPXLOPER XL_SchedulerType,
        LPXLOPER XL_SchedulerDataId,
		LPXLOPER XL_SpaceDataId,
		LPXLOPER XL_SchemeDataId,
		LPXLOPER XL_BoundCondName)
{
	ADD_LOG("Local_PdeND_Create");
	bool PersistentInXL = true;
	return Local_PdeND_Common(
		XL_MethodName,
		XL_SchedulerType,
        XL_SchedulerDataId,
		XL_SpaceDataId,
		XL_SchemeDataId,
		XL_BoundCondName,
		PersistentInXL);
}


///////////////////////////////////
/// version for non persistentInXl
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_PXL_PdeND_Create(
		LPXLOPER XL_MethodName,
        LPXLOPER XL_SchedulerType,
        LPXLOPER XL_SchedulerDataId,
		LPXLOPER XL_SpaceDataId,
		LPXLOPER XL_SchemeDataId,
		LPXLOPER XL_BoundCondName,
		LPXLOPER XL_Lambda)
{
	ADD_LOG("Local_PXL_PdeND_Create");
	bool PersistentInXL = false;
	return Local_PdeND_Common(
		XL_MethodName,
		XL_SchedulerType,
        XL_SchedulerDataId,
		XL_SpaceDataId,
		XL_SchemeDataId,
		XL_BoundCondName,
		PersistentInXL);
}





///----------------------------------------------
///----------------------------------------------
///             Random Generator
/// Inputs :
///     
///		genType,
///		algoType,
///		baseGenId,
///		seed
///----------------------------------------------
///----------------------------------------------

/////////////////////////////////////////////////////////////
/// central function that does the creation of the XL function
/////////////////////////////////////////////////////////////
LPXLOPER Local_RandomGen_Common(
	LPXLOPER XL_genType,
	LPXLOPER XL_algo,
	LPXLOPER XL_baseGen1Id,
	LPXLOPER XL_seed,
	LPXLOPER XL_dim,
	LPXLOPER XL_factorDim,
	LPXLOPER XL_nbOfPoints,
	LPXLOPER XL_nbStdDevs,
	LPXLOPER XL_baseGen2Id,
	LPXLOPER XL_firstNbTimes,
	LPXLOPER XL_order,
	LPXLOPER XL_firstSimulations,
	LPXLOPER XL_firstNbDims,
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

		CCString C_genTypeCCStr;
		XL_readStrCellWD(XL_genType,C_genTypeCCStr,"UnknownBaseGenAlgorithm", " ARM_ERR: gentype: String expected",C_result);
		string C_genType = CCSTringToSTLString(C_genTypeCCStr);

		CCString C_algoCCStr;
		XL_readStrCellWD(XL_algo,C_algoCCStr,"UnknownTransformAlgo", " ARM_ERR: algo: String expected",C_result);
		string C_algo = CCSTringToSTLString(C_algoCCStr);

		long C_BaseGen1Id;
		XL_GETOBJIDWD( XL_baseGen1Id, C_BaseGen1Id,	"NULL OBJECT", " ARM_ERR: Base gen 1 id: object expected",	C_result);

		long C_BaseGen2Id;
		XL_GETOBJIDWD( XL_baseGen2Id, C_BaseGen2Id,	"NULL OBJECT", " ARM_ERR: Base gen 2 id: object expected",	C_result);

		double C_seed;
		double C_seedDefault=-1;
		XL_readNumCellWD(XL_seed, C_seed, C_seedDefault, " ARM_ERR: seed : numeric expected", C_result);
		int C_seedInt = C_seed;

		double C_dim;
		double C_dimDefault=1;
		XL_readNumCellWD(XL_dim, C_dim, C_dimDefault, " ARM_ERR: dim : numeric expected", C_result);
		int C_dimInt = C_dim;

		double C_factorDim;
		double C_factorDimDefault=1;
		XL_readNumCellWD(XL_factorDim, C_factorDim, C_factorDimDefault, " ARM_ERR: dim : numeric expected", C_result);
		int C_factorDimInt = C_factorDim;

		double C_nbOfPoints;
		double C_nbOfPointsDefault=10;
		XL_readNumCellWD(XL_nbOfPoints, C_nbOfPoints, C_nbOfPointsDefault, " ARM_ERR: nbOfPoints : numeric expected", C_result);
		int C_nbOfPointsInt = C_nbOfPoints;

		double C_nbStdDevs;
		double C_nbStdDevsDefault=4.0;
		XL_readNumCellWD(XL_nbStdDevs, C_nbStdDevs, C_nbStdDevsDefault, " ARM_ERR: nbStdDevs : numeric expected", C_result);

		double C_firstNbTimes;
		int C_firstNbTimesInt;
		double C_firstNbTimesDefault=0;
		XL_readNumCellWD(XL_firstNbTimes, C_firstNbTimes, C_firstNbTimesDefault, " ARM_ERR: firstNbTimes : numeric expected", C_result);
		C_firstNbTimesInt = C_firstNbTimes;

		CCString C_orderStr;
		XL_readStrCellWD(XL_order,C_orderStr,"BucketOrder", " ARM_ERR: order: string expected",C_result);
		string C_order = CCSTringToSTLString(C_orderStr);

		double C_firstSimulations;
		int C_firstSimulationsInt;
		double C_firstSimulationsDef=0;
		XL_readNumCellWD(XL_firstSimulations, C_firstSimulations, C_firstSimulationsDef, " ARM_ERR: firstSimulations : numeric expected", C_result);
		C_firstSimulationsInt = C_firstSimulations;

		double C_firstNbDims;
		int C_firstNbDimsInt;
		double C_firstNbDimsDefault=0;
		XL_readNumCellWD(XL_firstNbDims, C_firstNbDims, C_firstNbDimsDefault, " ARM_ERR: firstNbDims : numeric expected", C_result);
		C_firstNbDimsInt = C_firstNbDims;

		/// use the concept of Functor to transfer the knowledge of
		/// a function with a context
		exportFunc13Args< string, string, long, long, int, int, int, int, double, int, int, string, int>
			ourFunc(
			C_genType,
			C_algo,
			C_BaseGen1Id,
			C_BaseGen2Id,
			(int)C_seedInt,
			(int)C_dimInt,
			(int)C_factorDimInt,
			(int)C_nbOfPointsInt, 
			C_nbStdDevs,
			C_firstNbTimesInt,
			C_firstNbDimsInt,
			C_order,
			C_firstSimulationsInt,
			ARMLOCAL_RandGen_Create);

		/// call the general function
		fillXL_Result( LOCAL_RANDGEN_CLASS, ourFunc, C_result, XL_result, PersistentInXL );
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_RandomGen_Common" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}
							 

///////////////////////////////////
/// version that takes into account 
/// previous creation of object
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_RandGen_Create(
	LPXLOPER XL_genType,
	LPXLOPER XL_algo,
	LPXLOPER XL_baseGen1Id,
	LPXLOPER XL_seed,
	LPXLOPER XL_dim,
	LPXLOPER XL_factorDim,
	LPXLOPER XL_nbOfPoints,
	LPXLOPER XL_nbStdDevs,
	LPXLOPER XL_baseGen2Id,
	LPXLOPER XL_firstNbTimes,
	LPXLOPER XL_order,
	LPXLOPER XL_firstSimulations,
	LPXLOPER XL_firstNbDims)
{
	ADD_LOG("Local_RandGen_Create");
	bool PersistentInXL = true;
	return Local_RandomGen_Common(
		XL_genType,
		XL_algo,
		XL_baseGen1Id,
		XL_seed,
		XL_dim,
		XL_factorDim,
		XL_nbOfPoints,
		XL_nbStdDevs,
		XL_baseGen2Id,
		XL_firstNbTimes,
		XL_order,
		XL_firstSimulations,
		XL_firstNbDims,
		PersistentInXL);
}


///////////////////////////////////
/// version for non persistentInXl
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_PXL_RandGen_Create(
	LPXLOPER XL_genType,
	LPXLOPER XL_algo,
	LPXLOPER XL_baseGen1Id,
	LPXLOPER XL_seed,
	LPXLOPER XL_dim,
	LPXLOPER XL_factorDim,
	LPXLOPER XL_nbOfPoints,
	LPXLOPER XL_nbStdDevs,
	LPXLOPER XL_baseGen2Id,
	LPXLOPER XL_firstNbTimes,
	LPXLOPER XL_order,
	LPXLOPER XL_firstSimulations,
	LPXLOPER XL_firstNbDims)
{
	ADD_LOG("Local_PXL_RandGen_Create");
	bool PersistentInXL = false;
	return Local_RandomGen_Common(
		XL_genType,
		XL_algo,
		XL_baseGen1Id,
		XL_seed,
		XL_dim,
		XL_factorDim,
		XL_nbOfPoints,
		XL_nbStdDevs,
		XL_baseGen2Id,
		XL_firstNbTimes,
		XL_order,
		XL_firstSimulations,
		XL_firstNbDims,
		PersistentInXL );
}

///----------------------------------------------
///----------------------------------------------
///             Simple Random Generator
/// Inputs :
///     
///		GenMode (Traditional, Fast)
///		FirstNbSimulations
///		FirstNbDims
///		SkipNbStdDevs
///----------------------------------------------
///----------------------------------------------

/////////////////////////////////////////////////////////////
/// central function that does the creation of the XL function
/////////////////////////////////////////////////////////////
LPXLOPER Local_SimpleRandomGen_Common(
	LPXLOPER XL_genType1,
	LPXLOPER XL_genType2,
	LPXLOPER XL_algo1,
	LPXLOPER XL_algo2,
	LPXLOPER XL_firstNbTimes,
	LPXLOPER XL_firstNbDims,
	LPXLOPER XL_isAntithetic,
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

		CCString C_genType1CCStr;
		XL_readStrCellWD(XL_genType1,C_genType1CCStr,"NR_Ran2", " ARM_ERR: gentype1: String expected",C_result);
		string C_genType1 = CCSTringToSTLString(C_genType1CCStr);

		CCString C_genType2CCStr;
		XL_readStrCellWD(XL_genType2,C_genType2CCStr,"Sobol", " ARM_ERR: gentype2: String expected",C_result);
		string C_genType2 = CCSTringToSTLString(C_genType2CCStr);

		CCString C_algo1CCStr;
		XL_readStrCellWD(XL_algo1,C_algo1CCStr,"BoxMuller", " ARM_ERR: algo1: String expected",C_result);
		string C_algo1 = CCSTringToSTLString(C_algo1CCStr);

		CCString C_algo2CCStr;
		XL_readStrCellWD(XL_algo2,C_algo2CCStr,"InvNormCum", " ARM_ERR: algo2: String expected",C_result);
		string C_algo2 = CCSTringToSTLString(C_algo2CCStr);

		double C_firstNbTimes;
		int C_firstNbTimesInt;
		double C_firstNbTimesDefault=0;
		XL_readNumCellWD(XL_firstNbTimes, C_firstNbTimes, C_firstNbTimesDefault, " ARM_ERR: firstNbTimes : numeric expected", C_result);
		C_firstNbTimesInt = C_firstNbTimes;

		double C_firstNbDims;
		int C_firstNbDimsInt;
		double C_firstNbDimsDefault=0;
		XL_readNumCellWD(XL_firstNbDims, C_firstNbDims, C_firstNbDimsDefault, " ARM_ERR: firstNbDims : numeric expected", C_result);
		C_firstNbDimsInt = C_firstNbDims;

		CCString C_isAntitheticCCStr;
		XL_readStrCellWD(XL_isAntithetic,C_isAntitheticCCStr,"Y", " ARM_ERR: isAntithetic: String expected",C_result);
		string C_isAntithetic = CCSTringToSTLString(C_isAntitheticCCStr);

		/// use the concept of Functor to transfer the knowledge of
		/// a function with a context
		exportFunc7Args< string, string, string, string, int, int, string>
			ourFunc(
			C_genType1,
			C_genType2,
			C_algo1,
			C_algo2,
			C_firstNbTimesInt,
			C_firstNbDimsInt,
			C_isAntithetic,
			ARMLOCAL_SimpleRandGen_Create);

		/// call the general function
		fillXL_Result( LOCAL_RANDGEN_CLASS, ourFunc, C_result, XL_result, PersistentInXL );
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_RandomGen_Common" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}
							 

///////////////////////////////////
/// version that takes into account 
/// previous creation of object
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_SimpleRandomGen_Create(
	LPXLOPER XL_genType1,
	LPXLOPER XL_genType2,
	LPXLOPER XL_algo1,
	LPXLOPER XL_algo2,
	LPXLOPER XL_firstNbTimes,
	LPXLOPER XL_firstNbDims,
	LPXLOPER XL_isAntithetic)
{
	ADD_LOG("Local_SimpleRandomGen_Create");
	bool PersistentInXL = true;
	return Local_SimpleRandomGen_Common(
		XL_genType1,
		XL_genType2,
		XL_algo1,
		XL_algo2,
		XL_firstNbTimes,
		XL_firstNbDims,
		XL_isAntithetic,
		PersistentInXL);
}


///////////////////////////////////
/// version for non persistentInXl
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_PXL_SimpleRandomGen_Create(
	LPXLOPER XL_genType1,
	LPXLOPER XL_genType2,
	LPXLOPER XL_algo1,
	LPXLOPER XL_algo2,
	LPXLOPER XL_firstNbTimes,
	LPXLOPER XL_firstNbDims,
	LPXLOPER XL_isAntithetic)
{
	ADD_LOG("Local_PXL_SimpleRandomGen_Create");
	bool PersistentInXL = false;
	return Local_SimpleRandomGen_Common(
		XL_genType1,
		XL_genType2,
		XL_algo1,
		XL_algo2,
		XL_firstNbTimes,
		XL_firstNbDims,
		XL_isAntithetic,
		PersistentInXL);
}


///----------------------------------------------
///----------------------------------------------
///             Random Generator Draw
/// Inputs :
///     
///		genId,
///		size,
///----------------------------------------------
///----------------------------------------------

/////////////////////////////////////////////////////////////
/// central function that does the creation of the XL function
/////////////////////////////////////////////////////////////
LPXLOPER Local_RandomGen_DrawVector(
	LPXLOPER XL_RandGenId,
	LPXLOPER XL_size )
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

		long C_RandGenId;
		XL_GETOBJID( XL_RandGenId, C_RandGenId,	" ARM_ERR: random nb gen id: object expected",	C_result);

		double C_size;
		double C_sizeDefault=10.;
		XL_readNumCellWD(XL_size, C_size, C_sizeDefault, " ARM_ERR: size: numeric expected", C_result);

		VECTOR<double> C_DataResult;
		if( ARMLOCAL_RandGen_DrawVector( C_RandGenId, C_size, C_DataResult, C_result ) == ARM_OK )
		{
			/// add these additional lines 
			/// to display blank lines
			const int additionalLinesNb = 100;
			bool fillWithBlank = true;
			FreeCurCellContent ();
			XL_writeNumVectorWithOptions( XL_result, C_DataResult, " ARM_ERR: Could not get result data", C_result, additionalLinesNb, fillWithBlank );
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_RandomGen_DrawVector" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}

///----------------------------------------------
///----------------------------------------------
///             AMC Andersen
/// Inputs :
///     
///		Nothing Yet
///		But it will change for sure
///----------------------------------------------
///----------------------------------------------							 

/////////////////////////////////////////////////////////////
/// function to create an Andersen Method
/////////////////////////////////////////////////////////////
LPXLOPER Local_AMCAndersen_Common( LPXLOPER XL_ItersNb, 
	LPXLOPER XL_sortedMaximisation,
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

		double C_ItersNb;
		XL_readNumCell( XL_ItersNb, C_ItersNb, " ARM_ERR: ItersNb: numeric expected", C_result);

		double C_sortedMaximisationDble;
		double C_sortedMaximisationDefaultDble=0;
		XL_readNumCellWD( XL_sortedMaximisation, C_sortedMaximisationDble, C_sortedMaximisationDefaultDble, " ARM_ERR: sorted maximisation: numeric expected", C_result);
		bool C_sortedMaximisation= C_sortedMaximisationDble!=0;

		/// a function with a context
		exportFunc2Args< double, bool > ourFunc(  C_ItersNb, C_sortedMaximisation, ARMLOCAL_AMCAndersen_Create );

		/// call the general function
		fillXL_Result( LOCAL_AMCANDERSEN_CLASS, ourFunc, C_result, XL_result, PersistentInXL );

	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_LinSurface_Create_Common" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}



///////////////////////////////////
/// Create an AMCAndersen
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_AMCAndersen_Create( LPXLOPER XL_ItersNb, LPXLOPER XL_sortedMaximisation )
{
	ADD_LOG("Local_AMCAndersen_Create");
	bool PersistentInXL = true;
	return Local_AMCAndersen_Common( XL_ItersNb, XL_sortedMaximisation, PersistentInXL );
}

__declspec(dllexport) LPXLOPER WINAPI Local_PXL_AMCAndersen_Create( LPXLOPER XL_ItersNb, LPXLOPER XL_sortedMaximisation )
{
	ADD_LOG("Local_PXL_AMCAndersen_Create");
	bool PersistentInXL = false;
	return Local_AMCAndersen_Common( XL_ItersNb, XL_sortedMaximisation, PersistentInXL );
}


///----------------------------------------------
///----------------------------------------------
///             AMC LongstaffSchwartz
/// Inputs :
///     
///		Nothing Yet
///		But it will change for sure
///----------------------------------------------
///----------------------------------------------							 

/////////////////////////////////////////////////////////////
/// function to create an LongstaffSchwartz Method
/////////////////////////////////////////////////////////////
LPXLOPER Local_AMCLongstaffSchwartz_Common( 
	LPXLOPER XL_ItersNb, 
	LPXLOPER XL_RegMode, 
	LPXLOPER XL_Span,
	LPXLOPER XL_IsAutomatic,
	LPXLOPER XL_Degree,
	bool PersistentInXL )
{	
	/// this is defined first because it is used in XL macros
	static XLOPER XL_result;

	/// Get the variables from the XLOper variables
	ARM_result C_result;;
	
	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
		/// to avoid computation if called by the wizard
		ARM_NOCALCIFWIZ();
		
		/// this is used by macros 
		/// and therefore this has to be defined
		static int error;
		static char* reason = "";

		double C_ItersNb;
		XL_readNumCell( XL_ItersNb, C_ItersNb, " ARM_ERR: ItersNb: numeric expected", C_result);

		CCString C_RegModeCCStr;
		XL_readStrCellWD(XL_RegMode,C_RegModeCCStr,"LS", " ARM_ERR: RegMode: String expected",C_result);
		string C_RegMode = CCSTringToSTLString(C_RegModeCCStr);

		double spanDef=0.8;
		double C_Span;
		XL_readNumCellWD( XL_Span, C_Span, spanDef, " ARM_ERR: Span: numeric expected", C_result);


		CCString C_IsAutomaticStr, C_IsAutomaticDefStr = "N";
		XL_readStrCellWD(XL_IsAutomatic, C_IsAutomaticStr, C_IsAutomaticDefStr, " ARM_ERR: IsAutomatic: String expected",C_result);
		string C_IsAutomatic = CCSTringToSTLString(C_IsAutomaticStr);

		double C_DegreeDbl, C_DegreeDefDbl=3;
		XL_readNumCellWD( XL_Degree, C_DegreeDbl, C_DegreeDefDbl, " ARM_ERR: Degree: numeric expected", C_result);
		int C_Degree = C_DegreeDbl;

		/// a function with a context
		exportFunc5Args< double, string, double, string, int > ourFunc(  C_ItersNb, C_RegMode, C_Span, C_IsAutomatic, C_Degree, ARMLOCAL_AMCLongstaffSchwartz_Create );

		/// call the general function
		fillXL_Result( LOCAL_AMCLS_CLASS, ourFunc, C_result, XL_result, PersistentInXL );

	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_LinSurface_Create_Common" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}



///////////////////////////////////
/// Create an AMCLongstaffSchwartz
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_AMCLongstaffSchwartz_Create( 
LPXLOPER XL_ItersNb, 
LPXLOPER XL_RegMode, 
LPXLOPER XL_Span, 
LPXLOPER XL_IsAutomatic, 
LPXLOPER XL_Degree )
{
	ADD_LOG("Local_AMCLongstaffSchwartz_Create");
	bool PersistentInXL = true;
	return Local_AMCLongstaffSchwartz_Common( XL_ItersNb, XL_RegMode, XL_Span, XL_IsAutomatic, XL_Degree, PersistentInXL );
}

__declspec(dllexport) LPXLOPER WINAPI Local_PXL_AMCLongstaffSchwartz_Create( 
LPXLOPER XL_ItersNb, 
LPXLOPER XL_RegMode, 
LPXLOPER XL_Span, 
LPXLOPER XL_IsAutomatic, 
LPXLOPER XL_Degree )
{
	ADD_LOG("Local_PXL_AMCLongstaffSchwartz_Create");
	bool PersistentInXL = false;
	return Local_AMCLongstaffSchwartz_Common( XL_ItersNb, XL_RegMode, XL_Span, XL_IsAutomatic, XL_Degree, PersistentInXL );
}


///----------------------------------------------
///----------------------------------------------
///             Tree ND (2nd generation)
///----------------------------------------------
///----------------------------------------------

//////////////////////////////////////////////////////////////////////////////
/// central function that does the creation of the XL function with Description
/////////////////////////////////////////////////////////////////////////////
LPXLOPER Local_TreeND_Common(
	LPXLOPER XL_NbDims,
	LPXLOPER XL_SchedulerType,
    LPXLOPER XL_SchedulerDatas,
	LPXLOPER XL_SamplerType,
    LPXLOPER XL_SamplerDatas,
    LPXLOPER XL_TruncatorType,
	LPXLOPER XL_TruncatorDatas,
    LPXLOPER XL_ProbasFlag,
	LPXLOPER XL_ReconnectorType,
	LPXLOPER XL_SmootherType,
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

		double C_NbDimsDble;
		XL_readNumCell(XL_NbDims,C_NbDimsDble," ARM_ERR: nb of dimensions: int expected",C_result);
        long C_NbDims = (long)C_NbDimsDble;


        CCString C_SchedulerTypeStr;
		CCString defaultSchedulerTypeStr = "ConstantVarianceMeanRevertingScheduler";
		XL_readStrCellWD(XL_SchedulerType,C_SchedulerTypeStr,defaultSchedulerTypeStr," ARM_ERR: Scheduler type: String expected",C_result);
		long C_SchedulerType;
		if( (C_SchedulerType = ARM_ConvGPSchedulerType( C_SchedulerTypeStr, C_result)) == ARM_DEFAULT_ERR )
		{
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}

		VECTOR<double> C_SchedulerDatas;
		VECTOR<double> defaultSchedulerDatas(3);
        defaultSchedulerDatas[0]=25;    /// required number of steps
        defaultSchedulerDatas[1]=1;     /// minimum number of steps before 1st event
        defaultSchedulerDatas[2]=0.001; /// minimum of StdDev per year for default schedule construction
		XL_readNumVectorWD(XL_SchedulerDatas,C_SchedulerDatas,defaultSchedulerDatas," ARM_ERR: Scheduler datas: array of numeric expected",C_result);


		CCString C_SamplerTypeStr;
		CCString defaultSamplerTypeStr = "MeanRevertingSampler";
		XL_readStrCellWD(XL_SamplerType,C_SamplerTypeStr,defaultSamplerTypeStr," ARM_ERR: Sampler type: String expected",C_result);
		long C_SamplerType;
		if( (C_SamplerType = ARM_ConvGPSamplerType( C_SamplerTypeStr, C_result)) == ARM_DEFAULT_ERR )
		{
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}

		VECTOR<double> C_SamplerDatas;
		VECTOR<double> defaultSamplerDatas(1);
        defaultSamplerDatas[0]=0.001; /// minimum of StdDev per year for space steps
		XL_readNumVectorWD(XL_SamplerDatas,C_SamplerDatas,defaultSamplerDatas," ARM_ERR: Sampler datas: array of numeric expected",C_result);


		CCString C_TruncatorTypeStr;
		CCString defaultTruncatorTypeStr = "StdDevTruncator";
		XL_readStrCellWD(XL_TruncatorType,C_TruncatorTypeStr,defaultTruncatorTypeStr," ARM_ERR: Truncator type: String expected",C_result);
		long C_TruncatorType;
		if( (C_TruncatorType = ARM_ConvGPTruncatorType( C_TruncatorTypeStr, C_result)) == ARM_DEFAULT_ERR )
		{
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}

		VECTOR<double> C_TruncatorDatas;
		VECTOR<double> defaultTruncatorDatas(1);
		defaultTruncatorDatas[0]=5.0;   /// number of StdDev
		XL_readNumVectorWD(XL_TruncatorDatas,C_TruncatorDatas,defaultTruncatorDatas," ARM_ERR: Truncator datas: array of numeric expected",C_result);

		double C_ProbasFlag;
        double defaultProbasFlag=0;
		XL_readNumCellWD(XL_ProbasFlag,C_ProbasFlag,defaultProbasFlag," ARM_ERR: probas computation flag : double expected",C_result);
        bool probasFlag = (C_ProbasFlag != 0);

        CCString C_ReconnectorTypeStr;
		CCString defaultReconnectorTypeStr = "MeanReconnector";
		XL_readStrCellWD(XL_ReconnectorType,C_ReconnectorTypeStr,defaultReconnectorTypeStr," ARM_ERR: Reconnector type: String expected",C_result);
		long C_ReconnectorType;
		if( (C_ReconnectorType = ARM_ConvGPReconnectorType( C_ReconnectorTypeStr, C_result)) == ARM_DEFAULT_ERR )
		{
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}

        CCString C_SmootherTypeStr;
		CCString defaultSmootherTypeStr = "DoNothingSmoother";
		XL_readStrCellWD(XL_SmootherType,C_SmootherTypeStr,defaultSmootherTypeStr," ARM_ERR: Smoother type: String expected",C_result);
		long C_SmootherType;
		if( (C_SmootherType = ARM_ConvGPSmootherType( C_SmootherTypeStr, C_result)) == ARM_DEFAULT_ERR )
		{
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}

		/// use the concept of Functor to transfer the knowledge of
		/// a function with a context
		exportFunc10Args< long, long, VECTOR<double>, long, VECTOR<double>, long,  VECTOR< double >, bool, long, long > ourFunc(C_NbDims,C_SchedulerType,C_SchedulerDatas,C_SamplerType,C_SamplerDatas,C_TruncatorType,C_TruncatorDatas,probasFlag,C_ReconnectorType,C_SmootherType,ARMLOCAL_TreeND_Create);

		/// call the general function
		fillXL_Result( LOCAL_TREEND_CLASS, ourFunc, C_result, XL_result, PersistentInXL );
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_Tree1D_Common" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}
							 

///////////////////////////////////
/// version that takes into account 
/// previous creation of object
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_TreeND_Create(
	LPXLOPER XL_NbDims,
	LPXLOPER XL_SchedulerType,
    LPXLOPER XL_SchedulerDatas,
	LPXLOPER XL_SamplerType,
	LPXLOPER XL_SamplerDatas,
    LPXLOPER XL_TruncatorType,
	LPXLOPER XL_TruncatorDatas,
    LPXLOPER XL_ProbasFlag,
	LPXLOPER XL_ReconnectorType,
	LPXLOPER XL_SmootherType)
{
	ADD_LOG("Local_TreeND_Create");
	bool PersistentInXL = true;
	return Local_TreeND_Common(
        XL_NbDims,
        XL_SchedulerType,
        XL_SchedulerDatas,
        XL_SamplerType,
        XL_SamplerDatas,
        XL_TruncatorType,
        XL_TruncatorDatas,
        XL_ProbasFlag,
	    XL_ReconnectorType,
	    XL_SmootherType,
        PersistentInXL );
}


///////////////////////////////////
/// version for non persistentInXl
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_PXL_TreeND_Create(
	LPXLOPER XL_NbDims,
	LPXLOPER XL_SchedulerType,
    LPXLOPER XL_SchedulerDatas,
	LPXLOPER XL_SamplerType,
	LPXLOPER XL_SamplerDatas,
    LPXLOPER XL_TruncatorType,
	LPXLOPER XL_TruncatorDatas,
    LPXLOPER XL_ProbasFlag,
	LPXLOPER XL_ReconnectorType,
	LPXLOPER XL_SmootherType)
{
	ADD_LOG("Local_PXL_TreeND_Create");
	bool PersistentInXL = false;
	return Local_TreeND_Common(
        XL_NbDims,
        XL_SchedulerType,
        XL_SchedulerDatas,
        XL_SamplerType,
        XL_SamplerDatas,
        XL_TruncatorType,
        XL_TruncatorDatas,
        XL_ProbasFlag,
	    XL_ReconnectorType,
	    XL_SmootherType,
        PersistentInXL );
}

///----------------------------------------------
///----------------------------------------------
///             CF Num Method
/// Inputs :
///     
///		Closed Form Num Method
///----------------------------------------------
///----------------------------------------------							 

/////////////////////////////////////////////////////////////
/// function to create a Closed Form Num Method
/////////////////////////////////////////////////////////////
LPXLOPER Local_CFMethod_Common(
	LPXLOPER XL_CFMethodName,
	LPXLOPER XL_IntegralParameters,
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

		//methodName
		CCString CFMethodNameXL;
		XL_readStrCellWD(XL_CFMethodName,CFMethodNameXL,"Unknown","ARM_ERR: MethodName: string expected",C_result);
		CFMethodNameXL.toUpper();
        string CFMethodName= CCSTringToSTLString(CFMethodNameXL);


		double integralParams;
		double default_IntegralParams=100.0;

		CCString IntegralParametersStr;
		long     IntegralParametersId;
		XL_readStrOrNumCellWD(XL_IntegralParameters,IntegralParametersStr,integralParams,default_IntegralParams,IntegralParametersId,
			" ARM_ERR: IntegralParameters:  Gen matrix Id expected",C_result);
		IntegralParametersId = IntegralParametersId == XL_TYPE_STRING ? LocalGetNumObjectId(IntegralParametersStr) : ARM_NULL_OBJECT;

		/// a function with a context
		exportFunc2Args<string,long> ourFunc(CFMethodName,
												IntegralParametersId,
												ARMLOCAL_CFMethod_Create);

		/// call the general function
		fillXL_Result( LOCAL_CFMETHOD_CLASS, ourFunc, C_result, XL_result, PersistentInXL );

	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_LinSurface_Create_Common" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}



///////////////////////////////////
/// Create an CFNumMethod
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_CFMethod_Create(LPXLOPER XL_CFMethodName, LPXLOPER XL_IntegralParameters)
{
	ADD_LOG("Local_CFMethod_Create");
	bool PersistentInXL = true;
	return Local_CFMethod_Common( XL_CFMethodName, XL_IntegralParameters, PersistentInXL );
}

__declspec(dllexport) LPXLOPER WINAPI Local_PXL_CFMethod_Create(LPXLOPER XL_CFMethodName, LPXLOPER XL_IntegralParameters)
{
	ADD_LOG("Local_PXL_CFMethod_Create");
	bool PersistentInXL = false;
	return Local_CFMethod_Common( XL_CFMethodName, XL_IntegralParameters, PersistentInXL );
}

///----------------------------------------------
///----------------------------------------------
///             Imp Sampler Optimizer
///----------------------------------------------
///----------------------------------------------

//////////////////////////////////////////////////////////////////////////////
/// central function that does the creation of the XL function with Description
/////////////////////////////////////////////////////////////////////////////
LPXLOPER Local_ImpSampler_Optimize_Common(
	LPXLOPER XL_GenSecId,
	LPXLOPER XL_ModelId,
    LPXLOPER XL_InitGuess,
	LPXLOPER XL_LowerBound,
    LPXLOPER XL_UpperBound,
	LPXLOPER XL_WithMC,
	LPXLOPER XL_NbSteps,
	LPXLOPER XL_Bootstrap,
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
		XL_GETOBJID( XL_GenSecId,	C_genSecId,	" ARM_ERR: Generic Security id: object expected",	C_result);

		CCString C_modelStrId;
		XL_readStrCell( XL_ModelId, C_modelStrId,	" ARM_ERR: Model Id: Object expected",		C_result);
		long C_modelId = LocalGetNumObjectId(C_modelStrId);

		// Init Guess
		double initGuess;
		CCString initGuessStr;
		long     initGuessId;
		XL_readStrOrNumCell(XL_InitGuess, initGuessStr, initGuess, initGuessId,
			   " ARM_ERR: init guess: numeric or curve Id expected",C_result);	
		if(initGuessId == XL_TYPE_STRING)
			initGuessId = LocalGetNumObjectId(initGuessStr);
		else
			initGuessId = ARM_NULL_OBJECT;

		// Lower Bound
		double lowerBound;
		CCString lowerBoundStr;
		long     lowerBoundId;
		XL_readStrOrNumCell(XL_LowerBound, lowerBoundStr, lowerBound, lowerBoundId,
			   " ARM_ERR: lower bound: numeric or curve Id expected",C_result);	
		if(lowerBoundId == XL_TYPE_STRING)
			lowerBoundId = LocalGetNumObjectId(lowerBoundStr);
		else
			lowerBoundId = ARM_NULL_OBJECT;

		// Lower Bound
		double upperBound;
		CCString upperBoundStr;
		long     upperBoundId;
		XL_readStrOrNumCell(XL_UpperBound, upperBoundStr, upperBound, upperBoundId,
			   " ARM_ERR: upper bound: numeric or curve Id expected",C_result);	
		if(upperBoundId == XL_TYPE_STRING)
			upperBoundId = LocalGetNumObjectId(upperBoundStr);
		else
			upperBoundId = ARM_NULL_OBJECT;

		// With MC
		CCString C_withMCStr, C_withMCDef="Y";
		XL_readStrCellWD(XL_WithMC, C_withMCStr, C_withMCDef, " ARM_ERR: with MC: string expected",C_result);
		string C_withMC = CCSTringToSTLString(C_withMCStr);

		// Nb Steps
		double C_NbStepsDbl, C_NbStepsDef=2.0;
		XL_readNumCellWD(XL_NbSteps,C_NbStepsDbl,C_NbStepsDef," ARM_ERR: nb of steps: int expected",C_result);
        long C_NbSteps = (long)C_NbStepsDbl;

		// Bootstrap
		CCString C_bootsrapStr, C_bootsrapDef="N";
		XL_readStrCellWD(XL_Bootstrap, C_bootsrapStr, C_bootsrapDef, " ARM_ERR: bootstrap: string expected",C_result);
		string C_bootstrap = CCSTringToSTLString(C_bootsrapStr);

		/// use the concept of Functor to transfer the knowledge of
		/// a function with a context
		exportFunc11Args< long, long, long, double, long, double, long, double, string, long, string > ourFunc(
			C_genSecId,
			C_modelId,
			initGuessId,
			initGuess,
			lowerBoundId,
			lowerBound,
			upperBoundId,
			upperBound,
			C_withMC,
			C_NbSteps,
			C_bootstrap,
			ARMLOCAL_ImpSampler_Optimize);

		CCString curClass = LocalGetStringObjectClassAndError(C_modelStrId);

		/// call the general function
		if( curClass != CCString(" ") )
			fillXL_Result( curClass, ourFunc, C_result, XL_result, PersistentInXL );
		/// bracket are necessary as ARM_ERR_AND_EXIT is a macro
		else
		{
			ARM_ERR_AND_EXIT( "Could not load the model ", C_result);
		}

		/// call the general function
		
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_ImpSampler_Optimize" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}
							 

///////////////////////////////////
/// version that takes into account 
/// previous creation of object
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_ImpSampler_Optimize(
	LPXLOPER XL_GenSecId,
	LPXLOPER XL_ModelId,
    LPXLOPER XL_InitGuess,
	LPXLOPER XL_LowerBound,
    LPXLOPER XL_UpperBound,
	LPXLOPER XL_WithMC,
	LPXLOPER XL_NbSteps,
	LPXLOPER XL_Bootstrap)
{
	ADD_LOG("Local_ImpSampler_Optimize");
	bool PersistentInXL = true;
	return Local_ImpSampler_Optimize_Common(
        XL_GenSecId,
		XL_ModelId,
		XL_InitGuess,
		XL_LowerBound,
		XL_UpperBound,
		XL_WithMC,
		XL_NbSteps,
		XL_Bootstrap,
        PersistentInXL );
}


///////////////////////////////////
/// version for non persistentInXl
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_PXL_ImpSampler_Optimize(
	LPXLOPER XL_GenSecId,
	LPXLOPER XL_ModelId,
    LPXLOPER XL_InitGuess,
	LPXLOPER XL_LowerBound,
    LPXLOPER XL_UpperBound,
	LPXLOPER XL_WithMC,
	LPXLOPER XL_NbSteps,
	LPXLOPER XL_Bootstrap)
{
	ADD_LOG("Local_PXL_ImpSampler_Optimize");
	bool PersistentInXL = false;
	return Local_ImpSampler_Optimize_Common(
        XL_GenSecId,
		XL_ModelId,
		XL_InitGuess,
		XL_LowerBound,
		XL_UpperBound,
		XL_WithMC,
		XL_NbSteps,
		XL_Bootstrap,
        PersistentInXL );
}
