/*
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 * $Log: ARM_xl_gp_calculators_local.cpp,v $
 * Revision 1.1  2004/02/07 15:08:43  ebenhamou
 * Initial version
 *
 */

#include <ARM\libarm_local\firstToBeIncluded.h>
#include "ARM_xl_gp_calculators_local.h"
#include <libCCxll\CCxll.h>
#include <ARM\libarm_local\ARM_local_gp_calculators.h>
#include <GP_Base\gpbase\datestripcombiner.h>
#include <GP_Infra\gpinfra\argconvdefault.h>
#include "ARM_xl_wrapper_local.h"
#include "ARM_xl_trycatch_local.h"
#include "ARM_local_interface.h"
#include <util\fromto.h>
#include <ARM_xl_gp_common.h>
#include <gpbase\gpmatrix.h>
#include "ARM_xl_gp_fctorhelper.h"
#include <gpbase\gpvector.h>
#include "gpbase/stringmanip.h"
#include <string>
#include "ARM_gp_local_interglob.h"

#include "util\tech_macro.h"

CC_USING_NS(std,string)

using ARM::ARM_DateStripCombiner;
using ARM::ARM_ArgConvReverse_MatFrequency;
using ARM::ARM_ArgConv_MatFrequency;
using ARM::ARM_ArgConv_IndexClass;
using ARM::ARM_ArgConv_InOut;
using ARM::stringGetUpper;


////////////////////////////////////////////////
///  THIS INTERFACE IS MORE FACTORISED THAN THE
///  REST OF ARM .... IN ORDER TO RECODE
///  ONLY THE PART TO DEFINE A FUNCTION CALL
///  WE PROVIDE FOR EACH A FUNCTOR
///  WE ALSO USE THE SAME API FOR PERSISTENT 
///  AND NON PERSISTENT ADDIN
///  BY COMMON DECISION.... WE KEEP THIS METHOD
///  TO THIS FILES AND DO NOT SPREAD ELSEWHERE
////////////////////////////////////////////////

///----------------------------------------------
///----------------------------------------------
///             date strip combiner
/// Inputs :
///     vector of date strip
///     field on which to merge
///----------------------------------------------
///----------------------------------------------
////////////////////////////////////////////////
/// very rapid definition of ourlonglongFunctor
////////////////////////////////////////////////


class dateStripCombinerFunc : public ARMResultLong2LongFunc
{
public:

	dateStripCombinerFunc( const VECTOR<long>& dateStripCombinerIds,
                    const string& mergeFuncName )
	:
	  C_dateStripCombinerIds(dateStripCombinerIds),
	  C_mergeFuncName(mergeFuncName)
    {};
	
	long operator()( ARM_result& result, long objId )
    {
		return ARMLOCAL_DateStripCombiner_Create(C_dateStripCombinerIds,
                                                 C_mergeFuncName,
                                                 result,
                                                 objId);			
	}

private:

	VECTOR<long> C_dateStripCombinerIds;
	string	C_mergeFuncName;
};


/////////////////////////////////////////////////////////////
/// central function that does the creation of the XL function
/////////////////////////////////////////////////////////////
LPXLOPER Local_DateStripCombiner_Common(
	LPXLOPER XL_dateStripCombinerIds,
	LPXLOPER XL_mergeFuncName,
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

    
		VECTOR<CCString> C_dateStripCombinerIdsString;
		XL_readStrVector (XL_dateStripCombinerIds,C_dateStripCombinerIdsString," ARM_ERR: date strips: array of object expected",DOUBLE_TYPE,C_result);
    
		size_t size = C_dateStripCombinerIdsString.size();
		VECTOR<long> C_dateStripCombinerIds;
		C_dateStripCombinerIds.resize(size);
		size_t i;
		for(i = 0; i < size; ++i )    
			C_dateStripCombinerIds[i] = LocalGetNumObjectId(C_dateStripCombinerIdsString[i]); 

		CCString C__mergeFuncName;
		CCString defaultmergeFuncName = "RESETDATE";
		XL_readStrCellWD(XL_mergeFuncName,C__mergeFuncName,defaultmergeFuncName," ARM_ERR: merge Func Name: String expected",C_result);


		/// use the concept of Functor to transfer the knowledge of
		/// a function with a context
		dateStripCombinerFunc ourFunc(C_dateStripCombinerIds,CCSTringToSTLString(C__mergeFuncName));

		/// call the general function
		fillXL_Result( LOCAL_DATESTRIPCOMBINER_CLASS, ourFunc, C_result, XL_result, PersistentInXL );
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
							 

///////////////////////////////////
/// version that takes into account 
/// previous creation of object
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_DateStripCombiner_Create(
	LPXLOPER XL_dateStripCombinerIds,
	LPXLOPER XL_mergeFuncName )
{
	ADD_LOG("Local_DateStripCombiner_Create");
	bool PersistentInXL = true;
	return Local_DateStripCombiner_Common(
		XL_dateStripCombinerIds,
		XL_mergeFuncName,
		PersistentInXL );
}


///////////////////////////////////
/// version for non persistentInXl
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_PXL_DateStripCombiner_Create(
	LPXLOPER XL_dateStripCombinerIds,
	LPXLOPER XL_mergeFuncName )
{
	ADD_LOG("Local_PXL_DateStripCombiner_Create");
	bool PersistentInXL = false;
	return Local_DateStripCombiner_Common(
		XL_dateStripCombinerIds,
		XL_mergeFuncName,
		PersistentInXL );
}


/////////////////////////////////////////
/// Local_DateStripGetData function
/////////////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_DataStripCombiner_GetData(
    LPXLOPER XL_dateStripCombinerId,
    LPXLOPER XL_dataType,
	LPXLOPER XL_dateStripNb )
{
	ADD_LOG("Local_DataStripCombiner_GetData");
	/// this is defined first because it is used in XL macros
	static XLOPER_Holder XL_resultHolder;
	XLOPER& XL_result = XL_resultHolder.GetResult();
	
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

		long C_dateStripCombinerId;
		XL_GETOBJID( XL_dateStripCombinerId,	C_dateStripCombinerId,	" ARM_ERR: date Strip combiner: Object expected",	C_result);

		long C_dataType;
		XL_GETMETHOD(XL_dataType,				C_dataType,				" ARM_ERR: data type: string expected",		C_result, ARM_ConvDateStripDataType );

		double C_dateStripNb;
		XL_readNumCell( XL_dateStripNb,			C_dateStripNb,			" ARM_ERR: date Strip nb: numeric expected",	C_result);
		VECTOR<double> C_DataResult;

		/// call the local function 
		long retCode;
		if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 1)
			retCode = ARMLOCAL_DateStripCombiner_GetData( C_dateStripCombinerId, C_dataType, C_dateStripNb, C_DataResult, C_result );
		else
			retCode = ARM_KO;

		/// feed the LPXLOPER object result 
		if (retCode == ARM_OK)
		{
			/// add these additional lines 
			/// to display blank lines
			const int additionalLinesNb = 100;
			bool fillWithBlank			= true;
			bool filterSpecificValue	= true;
			double specificValue		= ARM_DateStripCombiner::DateStripCombiner_BlankData;
			FreeCurCellContent ();
			XL_writeNumVectorWithOptionsAndFilter( XL_result, C_DataResult, " ARM_ERR: Could not get result data", C_result, additionalLinesNb, fillWithBlank, filterSpecificValue, specificValue );
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_DataStripCombiner_GetData" )
	
	return (LPXLOPER)&XL_result;
}


/////////////////////////////////////////
/// Local_DateStripGetData function
/////////////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_DataStripCombiner_GetMergeData(
    LPXLOPER XL_dateStripCombinerId )
{
	ADD_LOG("Local_DataStripCombiner_GetMergeData");
	/// this is defined first because it is used in XL macros
	static XLOPER_Holder XL_resultHolder;
	XLOPER& XL_result = XL_resultHolder.GetResult();
	
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

		long C_dateStripCombinerId;
		XL_GETOBJID( XL_dateStripCombinerId,	C_dateStripCombinerId,	" ARM_ERR: date Strip combiner: Object expected",	C_result);

		/// call the local function 
		VECTOR<double> C_DataResult;
		long retCode;
		if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 1)
			retCode = ARMLOCAL_DateStripCombiner_GetMergeData( C_dateStripCombinerId, C_DataResult, C_result );
		else
			retCode = ARM_KO;

		/// feed the LPXLOPER object result 
		if (retCode == ARM_OK)
		{
			/// add these additional lines 
			/// to display blank lines
			const int additionalLinesNb = 100;
			bool fillWithBlank			= true;
			bool filterSpecificValue	= true;
			double specificValue		= ARM_DateStripCombiner::DateStripCombiner_BlankData;
			FreeCurCellContent ();
			XL_writeNumVectorWithOptionsAndFilter( XL_result, C_DataResult, " ARM_ERR: Could not get result data", C_result, additionalLinesNb, fillWithBlank, filterSpecificValue, specificValue );
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_DataStripCombiner_GetMergeData" )
	
	return (LPXLOPER)&XL_result;
}

///----------------------------------------------
///----------------------------------------------
///             Gen Calculator
/// Inputs :
///     
///----------------------------------------------
///----------------------------------------------

__declspec(dllexport) LPXLOPER WINAPI Local_GenCalculator_GetPricingData(
	LPXLOPER XL_GenCalculatorId,
    LPXLOPER XL_KeyId )
{
	ADD_LOG("Local_GenCalculator_GetPricingData");
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
		XL_readStrCell( XL_GenCalculatorId, C_PricerString, " ARM_ERR: TARN Calculator Id: Object expected",C_result);
		long C_GenCalcId = LocalGetNumObjectId(C_PricerString);
		
		CCString C_KeyString;
		XL_readStrCell( XL_KeyId, C_KeyString, " ARM_ERR: Key Id: Object expected",C_result);
		
		/// call the function
		ARM_GramFctorArg argResult;
		
		long retCode = ARMLOCAL_GC_GetPricingData(
			C_GenCalcId,
			CCSTringToSTLString( C_KeyString ),
			argResult,
			C_result );
		
		if( retCode == ARM_KO )
		{
			ARM_ERR();
		}
		else
		{
			ARM_GramFunctorToXLOPER(argResult, XL_result, C_result);
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

///----------------------------------------------
///----------------------------------------------
///             Calculator Accessor
/// Inputs :
///     
///----------------------------------------------
///----------------------------------------------

class calculatorGetDataFunc : public ARMResultLong2LongFunc
{
public:
	calculatorGetDataFunc(
        long calcId,
        const string& getType)
    :
    C_calcId(calcId),
    C_getType(getType)
    {};
	
	long operator()( ARM_result& result, long objId ){
		return ARMLOCAL_GC_Get(
            C_calcId,
            C_getType,
            result,
            objId);
    }

private:
	long    C_calcId;
	string  C_getType;
};

LPXLOPER Local_GC_Get_Common(
	LPXLOPER XL_calcId,
	LPXLOPER XL_getType,
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

		long calcId;
		XL_GETOBJID( XL_calcId, calcId,	" ARM_ERR: Calculator: Object expected",C_result);

		CCString getTypeStr;
		XL_readStrCell(XL_getType,getTypeStr," ARM_ERR: Accessor Type: string expected",C_result);
        char* type=getTypeStr.c_str(); // à cause du new !!
		string getType(type);
        delete type;

		CCString calcGetClass(GCGetTypeToClass(getType,calcId).c_str());

		calculatorGetDataFunc ourFunc(calcId,getType);

		/// call the general function
		fillXL_Result( calcGetClass, ourFunc, C_result, XL_result, PersistentInXL );
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_GC_Get_Common" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}

__declspec(dllexport) LPXLOPER WINAPI Local_GenCalculator_GetData(
	LPXLOPER XL_calcId,
	LPXLOPER XL_getType)
{
	ADD_LOG("Local_GenCalculator_GetData");
	bool PersistentInXL = true;
	return Local_GC_Get_Common(
	    XL_calcId,
	    XL_getType,
        PersistentInXL );
}

__declspec(dllexport) LPXLOPER WINAPI Local_PXL_GenCalculator_GetData(
	LPXLOPER XL_calcId,
	LPXLOPER XL_getType)
{
	ADD_LOG("Local_PXL_GenCalculator_GetData");
	bool PersistentInXL = false;
	return Local_GC_Get_Common(
	    XL_calcId,
	    XL_getType,
        PersistentInXL );
}


class calculatorSetDataFunc : public ARMResultLong2LongFunc
{
public:
	calculatorSetDataFunc(
        long calculatorId,
        long dataToSetId,
		const vector< string >& mktDataKeys,
        bool isUpdated)
    :
    C_calculatorId(calculatorId),
    C_dataToSetId(dataToSetId),
	C_mktDataKeys(mktDataKeys),
    C_isUpdated(isUpdated)
    {};
	
	long operator()( ARM_result& result, long objId ){
		return ARMLOCAL_GC_Set(
            C_calculatorId,
            C_dataToSetId,
			C_mktDataKeys,
            C_isUpdated,
            result,
            objId);
	}

private:
	long			C_calculatorId;
	long			C_dataToSetId;
	vector<string>	C_mktDataKeys;
    bool			C_isUpdated;
};


LPXLOPER Local_GenCalculatorSet_Common(
	LPXLOPER XL_calculatorId,
	LPXLOPER XL_dataId,
	LPXLOPER XL_MktDataKeys,
	bool isUpdated,
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

		long calculatorId;
		XL_GETOBJID( XL_calculatorId, calculatorId, " ARM_ERR: Calculator: Object expected", C_result);

		long dataId;
		XL_GETOBJID( XL_dataId, dataId,	" ARM_ERR: Object expected", C_result);

		VECTOR<CCString> mktDataKeys;
		VECTOR<CCString> mktDataKeysDef(0);
		XL_readStrVectorWD(XL_MktDataKeys, mktDataKeys, mktDataKeysDef, " ARM_ERR: Market datas keys: array of string expected", DOUBLE_TYPE, C_result);

		vector<string> mktDataKeysSTL(mktDataKeys.size());
		for (size_t i=0; i<mktDataKeys.size(); ++i)
        {
            mktDataKeys[i].toUpper();
			mktDataKeysSTL[i] = CCSTringToSTLString(mktDataKeys[i]);
        }

		long objId = 0;
		string calculatorShortName;
		calculatorSetDataFunc ourFunc(calculatorId, dataId, mktDataKeysSTL, isUpdated);

        if (isUpdated)
        {
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
        else
        {
            /// General call through the functor with an object creation
			bool getNameFromResult = true;
		    fillXL_Result( CCString(calculatorShortName.c_str()), ourFunc, C_result, XL_result, PersistentInXL, getNameFromResult );
        }
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_GenCalculatorSet_Common" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}
							 
__declspec(dllexport) LPXLOPER WINAPI Local_GenCalculator_SetData(
	LPXLOPER XL_calculatorId,
	LPXLOPER XL_dataId,
	LPXLOPER XL_mktDataKeys)
{
	ADD_LOG("Local_GenCalculator_SetData");
	bool PersistentInXL = true;
    bool isUpdated = false;
	return Local_GenCalculatorSet_Common(
	    XL_calculatorId,
	    XL_dataId,
		XL_mktDataKeys,
        isUpdated,
        PersistentInXL );
}

__declspec(dllexport) LPXLOPER WINAPI Local_PXL_GenCalculator_SetData(
	LPXLOPER XL_calculatorId,
	LPXLOPER XL_dataId,
	LPXLOPER XL_mktDataKeys)
{
	ADD_LOG("Local_PXL_GenCalculator_SetData");
	bool PersistentInXL = false;
    bool isUpdated = false;
	return Local_GenCalculatorSet_Common(
	    XL_calculatorId,
	    XL_dataId,
		XL_mktDataKeys,
        isUpdated,
        PersistentInXL );
}
///----------------------------------------------
///----------------------------------------------
///             CRF Calculator
/// Inputs :
///     
///----------------------------------------------
///----------------------------------------------

__declspec(dllexport) LPXLOPER WINAPI Local_ARM_GetCcyFromGenCalculator(
          LPXLOPER XL_calcId,
          LPXLOPER XL_CcyType)
{
	ADD_LOG("Local_ARM_GetCcyFromGenCalculator");
	/// this is defined first because it is used in XL macros
	static XLOPER XL_result;
	
	/// Get the variables from the XLOper variables
	ARM_result C_result;
	
	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
		/// to avoid computation if called by the wizard
		ARM_NOCALCIFWIZ();
		
		static int error;
		static char* reason = "";

		CCString C_calcId;
		XL_readStrCell(XL_calcId,C_calcId," ARM_ERR: Gen Calculator : Object expected",C_result);
		long calcId = LocalGetNumObjectId (C_calcId);

        CCString C_CcyStr;
		XL_readStrCell( XL_CcyType, C_CcyStr, " ARM_ERR: Ccy Type: string expected",C_result);
		
		long retCode = ARMLOCAL_ARM_GetCcyFromGenCalculator(calcId,C_CcyStr, C_result);

		if(retCode == ARM_OK)
		{
			FreeCurCellErr ();
			XL_result.xltype = xltypeStr;
			XL_result.val.str = XL_StrC2StrPascal (C_result.getString ());
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_ARM_GetCcyFromGenCalculator" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


///----------------------------------------------
///----------------------------------------------
///             CRF Calculator
/// Inputs :
///     
///----------------------------------------------
///----------------------------------------------
class crfCalculatorFunc : public ARMResultLong2LongFunc
{
public:
	crfCalculatorFunc(
        double startDate,
        double endDate,
        double strike,
        long strikeId,
        long payRec,
        long cpnDayCount,
        long cpnFreq,
        double fixEndDate,
        long  fixDayCount,
        long cpnTiming,
        const string& cpnIndexTerm,
        long cpnIndexDayCount,
        const string& cpnResetCal,
        const string& cpnPayCal,
		long stubRule,
        long cpnResetGap,
        double leverage,
        long leverageId,
        double cpnMin,
        long cpnMinId,
        double cpnMax,
        long cpnMaxId,
        double fundSpread,
        long fundSpreadId,
        long  fundFreq,
        long  fundDayCount,
        double nominal,
        long nominalId,
        long exerGap,
        long nbNonCall,
        double exerFee,
        long exerFeeId,
        const vector< string >& flags,
        long mktDataManagerId,
        const vector< string >& keys,
        double fundnominal,
        long fundnominalId)
    :
    C_startDate(startDate),
    C_endDate(endDate),
    C_strike(strike),
    C_strikeId(strikeId),
    C_payRec(payRec),
    C_cpnDayCount(cpnDayCount),
    C_cpnFreq(cpnFreq),
    C_fixEndDate(fixEndDate),
    C_fixDayCount(fixDayCount),
    C_cpnTiming(cpnTiming),
    C_cpnIndexTerm(cpnIndexTerm),
    C_cpnIndexDayCount(cpnIndexDayCount),
    C_cpnResetCal(cpnResetCal),
    C_cpnPayCal(cpnPayCal),
	C_stubRule(stubRule),
    C_cpnResetGap(cpnResetGap),
    C_leverage(leverage),
    C_leverageId(leverageId),
    C_cpnMin(cpnMin),
    C_cpnMinId(cpnMinId),
    C_cpnMax(cpnMax),
    C_cpnMaxId(cpnMaxId),
    C_fundSpread(fundSpread),
    C_fundSpreadId(fundSpreadId),
    C_fundFreq(fundFreq),
    C_fundDayCount(fundDayCount),
    C_nominal(nominal),
    C_nominalId(nominalId),
    C_exerGap(exerGap),
    C_nbNonCall(nbNonCall),
    C_exerFee(exerFee),
    C_exerFeeId(exerFeeId),
    C_flags(flags),
    C_mktDataManagerId(mktDataManagerId),
    C_keys(keys),
    C_fundnominal(fundnominal),
    C_fundnominalId(fundnominalId)
    {};
	
	long operator()( ARM_result& result, long objId ){
		return ARMLOCAL_CRFCalculator_Create(
            C_startDate,
            C_endDate,
            C_strike,
            C_strikeId,
            C_payRec,
            C_fixEndDate,
            C_fixDayCount,
            C_cpnDayCount,
            C_cpnFreq,
            C_cpnTiming,
            C_cpnIndexTerm,
            C_cpnIndexDayCount,
            C_cpnResetCal,
            C_cpnPayCal,
			C_stubRule,
            C_cpnResetGap,
            C_leverage,
            C_leverageId,
            C_cpnMin,
            C_cpnMinId,
            C_cpnMax,
            C_cpnMaxId,
            C_fundSpread,
            C_fundSpreadId,
            C_fundFreq,
            C_fundDayCount,
            C_nominal,
            C_nominalId,
            C_exerGap,
            C_nbNonCall,
            C_exerFee,
            C_exerFeeId,
            C_flags,
            C_mktDataManagerId,
            C_keys,
            C_fundnominal,
            C_fundnominalId,
            result,
            objId);
    }

private:
    double              C_startDate;
    double              C_endDate;
    double              C_strike;
    long                C_strikeId;
    long                C_payRec;
    long                C_cpnDayCount;
    long                C_cpnFreq;
    double              C_fixEndDate;
    long                C_fixDayCount;
    long                C_cpnTiming;
    string              C_cpnIndexTerm;
    long                C_cpnIndexDayCount;
    string              C_cpnResetCal;
    string              C_cpnPayCal;
	long              C_stubRule;
    long                C_cpnResetGap;
    double              C_leverage;
    long                C_leverageId;
    double              C_cpnMin;
    long                C_cpnMinId;
    double              C_cpnMax;
    long                C_cpnMaxId;
    double              C_fundSpread;
    long                C_fundSpreadId;
    long                C_fundFreq;
    long                C_fundDayCount;
    double              C_nominal;
    long                C_nominalId;
    long                C_exerGap;
    long                C_nbNonCall;
    double              C_exerFee;
    long                C_exerFeeId;
    vector< string >    C_flags;
    long                C_mktDataManagerId;
    vector< string >    C_keys;
	string				C_exerCal;
    double              C_fundnominal;
    long                C_fundnominalId;
};


/////////////////////////////////////////////////////////////
/// central function that does the creation of the XL function
/////////////////////////////////////////////////////////////
LPXLOPER Local_CRFCalculator_Common(LPXLOPER XL_startDate,
									LPXLOPER XL_endDate,
									LPXLOPER XL_strike,
									LPXLOPER XL_payRec,
									LPXLOPER XL_cpnDatas,
									LPXLOPER XL_mktDatas,
									LPXLOPER XL_fixEndDate,
									LPXLOPER XL_fixDayCount,
									LPXLOPER XL_cpnResetGap,
									LPXLOPER XL_leverage,
									LPXLOPER XL_cpnMin,
									LPXLOPER XL_cpnMax,
									LPXLOPER XL_fundSpread,
									LPXLOPER XL_fundDatas,
									LPXLOPER XL_nominal,
									LPXLOPER XL_exerGap,
									LPXLOPER XL_nbNonCall,
									LPXLOPER XL_exerFee,
									LPXLOPER XL_flags,
                                    LPXLOPER XL_FundNominal,
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

    
		/// General CRF datas
		/// -----------------
		/// Start & End
		double startDate;
		XL_readNumCell(XL_startDate,startDate," ARM_ERR: start date: date expected",C_result);
		double endDate;
		XL_readNumCell(XL_endDate,endDate," ARM_ERR: end date: date expected",C_result);

		/// Fix rate or RF Strike
		double strike;
		CCString strikeStr;
		long     strikeId;
		XL_readStrOrNumCell(XL_strike, strikeStr, strike, strikeId,
			   " ARM_ERR: strike: numeric or refValue Id expected",C_result);	
		if(strikeId == XL_TYPE_STRING)
			strikeId = LocalGetNumObjectId(strikeStr);
		else
			strikeId = ARM_NULL_OBJECT;

		/// RF Leg Pay/Rec
		CCString payRecStr;
		long payRec;
		XL_readStrCell(XL_payRec,payRecStr," ARM_ERR: payer/receiver: string expected",C_result);
		if((payRec = ARM_ConvRecOrPay (payRecStr, C_result)) == ARM_DEFAULT_ERR)
		{
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}


		/// Required RF coupon datas
		/// ------------------------
		VECTOR<CCString> cpnDatas;
		XL_readStrVector (XL_cpnDatas,cpnDatas," ARM_ERR: Coupon datas: array of string expected",DOUBLE_TYPE,C_result);
		if(cpnDatas.size() < 2)
		{
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}

		/// RF coupon frequency (reset & payment)
		CCString cpnFreqXl=cpnDatas[0];
		long cpnFreq;
		if((cpnFreq = ARM_ConvFrequency (cpnFreqXl, C_result)) == ARM_DEFAULT_ERR)
		{
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}

		/// RF coupon day count
		CCString cpnDayCountXl=cpnDatas[1];
		long cpnDayCount = ARM_ConvDayCount (cpnDayCountXl);


		/// Market datas : curve , MRS and market models for OSW & CF pricings
		/// ------------------------------------------------------------------
		VECTOR<CCString> mktDatas;
		XL_readStrVector(XL_mktDatas,mktDatas," ARM_ERR: Market datas: array of string expected",DOUBLE_TYPE,C_result);
        if(mktDatas.size()<5)
		{
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}

		long mktDataManagerId = LocalGetNumObjectId (mktDatas[0]);
		vector< string > keys(mktDatas.size()-1);
        size_t i;
		for(i=1;i<mktDatas.size();++i)
			keys[i-1]=CCSTringToSTLString(mktDatas[i]);



		//// -------------------------------
		///  ----- OPTIONAL ARGUMENTS ------
		//// -------------------------------



		/// RF coupon datas
		/// ---------------
		/// Advance/Arrears RF coupon reset
		CCString cpnTimingXl("ADV");
		if(cpnDatas.size() >= 3)
			cpnTimingXl=cpnDatas[2];
		long cpnTiming = ARM_ConvPayResetRule (cpnTimingXl);

		/// RF index term
		CCString cpnIndexTermXl(ARM_ArgConvReverse_MatFrequency.GetString(cpnFreq).c_str());
		if(cpnDatas.size() >= 4)
			cpnIndexTermXl=cpnDatas[3];
		string cpnIndexTerm=CCSTringToSTLString(cpnIndexTermXl);

		/// RF index day count
		CCString cpnIndexDayCountXl("A360");
		if(cpnDatas.size() >= 5)
			cpnIndexDayCountXl=cpnDatas[4];
		long cpnIndexDayCount = ARM_ConvDayCount (cpnIndexDayCountXl);

		/// RF reset calendar
		CCString cpnResetCalXl("EUR");
		if(cpnDatas.size() >= 6)
			cpnResetCalXl=cpnDatas[5];
		string cpnResetCal = CCSTringToSTLString(cpnResetCalXl);

		/// RF payment calendar
		CCString cpnPayCalXl("EUR");
		if(cpnDatas.size() >= 7)
			cpnPayCalXl=cpnDatas[6];
		string cpnPayCal = CCSTringToSTLString(cpnPayCalXl);

		/// RF Stub Rule

		long stubRuleId;		
		CCString stubRuleXl("SS");
		if(cpnDatas.size() >= 8)
			stubRuleXl=cpnDatas[7];
		stubRuleId = ARM_ConvStubRule (stubRuleXl);

		/// RF reset gap
		double cpnResetGapD;
		double cpnResetGapDef=2;
		XL_readNumCellWD(XL_cpnResetGap,cpnResetGapD,cpnResetGapDef," ARM_ERR: coupon reset gap: numerical expected",C_result);
		long cpnResetGap = -(long)floor(cpnResetGapD);


		/// Fixed coupon datas
		/// ------------------
		/// Fixed leg end
		double fixEndDate;
		XL_readNumCellWD(XL_fixEndDate,fixEndDate,startDate," ARM_ERR: fixed leg end date: date expected",C_result);

		/// Fixed coupon day count
		CCString fixDayCountXl;
		XL_readStrCellWD(XL_fixDayCount,fixDayCountXl,cpnDayCountXl," ARM_ERR: fix cpn day count: string expected",C_result);
		long fixDayCount = ARM_ConvDayCount (fixDayCountXl);


		/// RF coupon profiles
		/// ------------------
		/// RF coupon index leverage
		double leverage;
		double leverageDef=1.0;
		CCString leverageStr;
		long     leverageId;
		XL_readStrOrNumCellWD(XL_leverage,leverageStr,leverage,leverageDef,leverageId,
			" ARM_ERR: strike: numeric or refValue Id expected",C_result);
		if(leverageId == XL_TYPE_STRING)
			leverageId = LocalGetNumObjectId(leverageStr);
		else
			leverageId = ARM_NULL_OBJECT;

		/// RF min coupon
		double cpnMin;
		double cpnMinDef=0.0;
		CCString cpnMinStr;
		long     cpnMinId;
		XL_readStrOrNumCellWD(XL_cpnMin,cpnMinStr,cpnMin,cpnMinDef,cpnMinId,
			" ARM_ERR: cpn min: numerical or refValue Id expected",C_result);
		if(cpnMinId == XL_TYPE_STRING)
			cpnMinId = LocalGetNumObjectId(cpnMinStr);
		else
			cpnMinId = ARM_NULL_OBJECT;

		/// RF max coupon
		double cpnMax;
		double cpnMaxDef=10000.0;
		CCString cpnMaxStr;
		long     cpnMaxId;
		XL_readStrOrNumCellWD(XL_cpnMax,cpnMaxStr,cpnMax,cpnMaxDef,cpnMaxId,
			" ARM_ERR: cpn max: numerical or refValue Id expected",C_result);
		if(cpnMaxId == XL_TYPE_STRING)
			cpnMaxId = LocalGetNumObjectId(cpnMaxStr);
		else
			cpnMaxId = ARM_NULL_OBJECT;


		/// Funding leg datas
		/// -----------------
		/// Funding coupon spread profile
		double fundSpread;
		double fundSpreadDef=0.0;
		CCString fundSpreadXl;
		long     fundSpreadId;
		XL_readStrOrNumCellWD(XL_fundSpread,fundSpreadXl,fundSpread,fundSpreadDef,fundSpreadId,
			" ARM_ERR: funding spread: numerical or refValue Id expected",C_result);
		if(fundSpreadId == XL_TYPE_STRING)
			fundSpreadId = LocalGetNumObjectId(fundSpreadXl);
		else
			fundSpreadId = ARM_NULL_OBJECT;

		VECTOR<CCString> fundDatas;
		XL_readStrVector (XL_fundDatas,fundDatas," ARM_ERR: Funding datas: array of string expected",DOUBLE_TYPE,C_result);

		/// Funding coupon frequency (reset & payment)
		CCString fundFreqXl(cpnFreqXl);
		if(fundDatas.size() >= 1)
			fundFreqXl=fundDatas[0];
		long fundFreq;
		if((fundFreq = ARM_ConvFrequency (fundFreqXl, C_result)) == ARM_DEFAULT_ERR)
		{
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}

		/// Funding coupon day count
		CCString fundDayCountXl("A360");
		if(fundDatas.size() >= 2)
			fundDayCountXl=fundDatas[1];
		long fundDayCount = ARM_ConvDayCount (fundDayCountXl);

		/// Funding & RF legs notional
		/// --------------------------
		double nominal;
		double nominalDef=100.0;
		CCString nominalStr;
		long     nominalId;
		XL_readStrOrNumCellWD(XL_nominal,nominalStr,nominal,nominalDef,nominalId,
			" ARM_ERR: nominal: numerical or refValue Id expected",C_result);
		if(nominalId == XL_TYPE_STRING)
			nominalId = LocalGetNumObjectId(nominalStr);
		else
			nominalId = ARM_NULL_OBJECT;

        double fundnominal;
		double fundnominalDef = ARM_MISSING_VALUE;
		CCString fundnominalStr;
		long     fundnominalId;
		XL_readStrOrNumCellWD(XL_FundNominal,fundnominalStr,fundnominal,fundnominalDef,fundnominalId,
			" ARM_ERR: nominal: numerical or refValue Id expected",C_result);
		if(fundnominalId == XL_TYPE_STRING)
			fundnominalId = LocalGetNumObjectId(fundnominalStr);
		else
			fundnominalId = ARM_NULL_OBJECT;



		/// Exercise datas
		/// --------------
		/// Notification lag
		double exerGapD;
		XL_readNumCellWD(XL_exerGap,exerGapD,cpnResetGapD," ARM_ERR: notification gap: numerical expected",C_result);
		long exerGap = -(long)floor(exerGapD);
		if (cpnTiming == K_ADVANCE && cpnResetGap < exerGap)
		{
			C_result.setMsg ("ARM_ERR: exercise must be anterior or equal to coupon reset");
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}

		/// Non callable period (w.r.t. RF coupon frequency)
		double nbNonCallD;
		double nbNonCallDDef=0.0;
		XL_readNumCellWD(XL_nbNonCall,nbNonCallD,nbNonCallDDef," ARM_ERR: nb of non call periods: numerical expected",C_result);
		long nbNonCall = (long)floor(nbNonCallD);

		/// Fees profile
		double exerFee;
		double exerFeeDef=0.0;
		CCString exerFeeStr;
		long     exerFeeId;
		XL_readStrOrNumCellWD(XL_exerFee,exerFeeStr,exerFee,exerFeeDef,exerFeeId,
			" ARM_ERR: exercise fees: numerical or refValue Id expected",C_result);
		if(exerFeeId == XL_TYPE_STRING)
			exerFeeId = LocalGetNumObjectId(exerFeeStr);
		else
			exerFeeId = ARM_NULL_OBJECT;


		/// Flags datas
		/// -----------
		VECTOR<CCString> C_flags;
		VECTOR<CCString> flagsDef(8,"");
        flagsDef[0]=CCString("EQUIVALENT");	/// Vol calib
        flagsDef[1]=CCString("UNKNOWN");	/// MRS calib
        flagsDef[2]=CCString("Y");			/// Cap calib
        flagsDef[3]=CCString("Y");			/// floor calib
        flagsDef[4]=CCString("CRF");		/// Product name
		flagsDef[5]=CCString("HWM1F");		/// modelType
		flagsDef[6]=CCString("Y");			/// Skew Calib for QGM Pricing
		flagsDef[7]=CCString("EQUIVALENT"); ///  Which strike have we use to calibrate MRS Parameter


		XL_readStrVectorWD(XL_flags,C_flags,flagsDef," ARM_ERR: flags: array of string expected",DOUBLE_TYPE,C_result);

		vector< string > flags(C_flags.size());
		for(i=0;i<C_flags.size();++i)
        {
            (C_flags[i]).toUpper();
			flags[i]=CCSTringToSTLString(C_flags[i]);
        }

		/// use the concept of Functor to transfer the knowledge of
		/// a function with a context
		crfCalculatorFunc ourFunc(
				startDate,
				endDate,
				strike,
				strikeId,
				payRec,
				cpnDayCount,
				cpnFreq,
				fixEndDate,
				fixDayCount,
				cpnTiming,
				cpnIndexTerm,
				cpnIndexDayCount,
				cpnResetCal,
				cpnPayCal,
				stubRuleId,
				cpnResetGap,
				leverage,
				leverageId,
				cpnMin,
				cpnMinId,
				cpnMax,
				cpnMaxId,
				fundSpread,
				fundSpreadId,
				fundFreq,
				fundDayCount,
				nominal,
				nominalId,
				exerGap,
				nbNonCall,
				exerFee,
				exerFeeId,
                flags,
				mktDataManagerId,
				keys,
                fundnominal,
                fundnominalId);

		/// call the general function
		fillXL_Result( LOCAL_GC_CRF_CLASS, ourFunc, C_result, XL_result, PersistentInXL );
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_CRFCalculator_Common" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}
							 

///////////////////////////////////
/// version that takes into account 
/// previous creation of object
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_CRFCalculator_Create(
    LPXLOPER XL_startDate,
    LPXLOPER XL_endDate,
    LPXLOPER XL_strike,
    LPXLOPER XL_payRec,
    LPXLOPER XL_cpnDatas,
    LPXLOPER XL_mktDatas,
    LPXLOPER XL_fixEndDate,
    LPXLOPER XL_fixDayCount,
    LPXLOPER XL_cpnResetGap,
    LPXLOPER XL_leverage,
    LPXLOPER XL_cpnMin,
    LPXLOPER XL_cpnMax,
    LPXLOPER XL_fundSpread,
    LPXLOPER XL_fundDatas,
    LPXLOPER XL_nominal,
    LPXLOPER XL_exerGap,
    LPXLOPER XL_nbNonCall,
    LPXLOPER XL_exerFee,
    LPXLOPER XL_flags,
    LPXLOPER XL_FundNominal)
{
	ADD_LOG("Local_CRFCalculator_Create");
	bool PersistentInXL = true;

	return Local_CRFCalculator_Common(
        XL_startDate,
        XL_endDate,
        XL_strike,
        XL_payRec,
        XL_cpnDatas,
        XL_mktDatas,
        XL_fixEndDate,
        XL_fixDayCount,
        XL_cpnResetGap,
        XL_leverage,
        XL_cpnMin,
        XL_cpnMax,
        XL_fundSpread,
        XL_fundDatas,
        XL_nominal,
        XL_exerGap,
        XL_nbNonCall,
        XL_exerFee,
        XL_flags,
        XL_FundNominal,
		PersistentInXL );
}



///////////////////////////////////
/// version for non persistentInXl
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_PXL_CRFCalculator_Create(
    LPXLOPER XL_startDate,
    LPXLOPER XL_endDate,
    LPXLOPER XL_strike,
    LPXLOPER XL_payRec,
    LPXLOPER XL_cpnDatas,
    LPXLOPER XL_mktDatas,
    LPXLOPER XL_fixEndDate,
    LPXLOPER XL_fixDayCount,
    LPXLOPER XL_cpnResetGap,
    LPXLOPER XL_leverage,
    LPXLOPER XL_cpnMin,
    LPXLOPER XL_cpnMax,
    LPXLOPER XL_fundSpread,
    LPXLOPER XL_fundDatas,
    LPXLOPER XL_nominal,
    LPXLOPER XL_exerGap,
    LPXLOPER XL_nbNonCall,
    LPXLOPER XL_exerFee,
    LPXLOPER XL_flags,
    LPXLOPER XL_FundNominal)
{
	ADD_LOG("Local_PXL_CRFCalculator_Create");
	bool PersistentInXL = false;
	return Local_CRFCalculator_Common(
        XL_startDate,
        XL_endDate,
        XL_strike,
        XL_payRec,
        XL_cpnDatas,
        XL_mktDatas,
        XL_fixEndDate,
        XL_fixDayCount,
        XL_cpnResetGap,
        XL_leverage,
        XL_cpnMin,
        XL_cpnMax,
        XL_fundSpread,
        XL_fundDatas,
        XL_nominal,
        XL_exerGap,
        XL_nbNonCall,
        XL_exerFee,
        XL_flags,
        XL_FundNominal,
		PersistentInXL );
}

/////////////////////////////////////////////////////////////
/// central function that does the creation of the XL function
/////////////////////////////////////////////////////////////
LPXLOPER Local_Crude_CRFCalculator_Common(LPXLOPER XL_startDate,
									LPXLOPER XL_endDate,
									LPXLOPER XL_strike,
									LPXLOPER XL_payRec,
									LPXLOPER XL_cpnDatas,
									LPXLOPER XL_fixEndDate,
									LPXLOPER XL_fixDayCount,
									LPXLOPER XL_cpnResetGap,
									LPXLOPER XL_leverage,
									LPXLOPER XL_cpnMin,
									LPXLOPER XL_cpnMax,
									LPXLOPER XL_fundSpread,
									LPXLOPER XL_fundDatas,
									LPXLOPER XL_nominal,
									LPXLOPER XL_exerGap,
									LPXLOPER XL_nbNonCall,
									LPXLOPER XL_exerFee,
                                    LPXLOPER XL_FundNominal,
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

    
		/// General CRF datas
		/// -----------------
		/// Start & End
		double startDate;
		XL_readNumCell(XL_startDate,startDate," ARM_ERR: start date: date expected",C_result);
		double endDate;
		XL_readNumCell(XL_endDate,endDate," ARM_ERR: end date: date expected",C_result);

		/// Fix rate or RF Strike
		double strike;
		CCString strikeStr;
		long     strikeId;
		XL_readStrOrNumCell(XL_strike, strikeStr, strike, strikeId,
			   " ARM_ERR: strike: numeric or refValue Id expected",C_result);	
		if(strikeId == XL_TYPE_STRING)
			strikeId = LocalGetNumObjectId(strikeStr);
		else
			strikeId = ARM_NULL_OBJECT;

		/// RF Leg Pay/Rec
		CCString payRecStr;
		long payRec;
		XL_readStrCell(XL_payRec,payRecStr," ARM_ERR: payer/receiver: string expected",C_result);
		if((payRec = ARM_ConvRecOrPay (payRecStr, C_result)) == ARM_DEFAULT_ERR)
		{
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}


		/// Required RF coupon datas
		/// ------------------------
		VECTOR<CCString> cpnDatas;
		XL_readStrVector (XL_cpnDatas,cpnDatas," ARM_ERR: Coupon datas: array of string expected",DOUBLE_TYPE,C_result);
		if(cpnDatas.size() < 2)
		{
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}

		/// RF coupon frequency (reset & payment)
		CCString cpnFreqXl=cpnDatas[0];
		long cpnFreq;
		if((cpnFreq = ARM_ConvFrequency (cpnFreqXl, C_result)) == ARM_DEFAULT_ERR)
		{
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}

		/// RF coupon day count
		CCString cpnDayCountXl=cpnDatas[1];
		long cpnDayCount = ARM_ConvDayCount (cpnDayCountXl);




		//// -------------------------------
		///  ----- OPTIONAL ARGUMENTS ------
		//// -------------------------------



		/// RF coupon datas
		/// ---------------
		/// Advance/Arrears RF coupon reset
		CCString cpnTimingXl("ADV");
		if(cpnDatas.size() >= 3)
			cpnTimingXl=cpnDatas[2];
		long cpnTiming = ARM_ConvPayResetRule (cpnTimingXl);

		/// RF index term
		CCString cpnIndexTermXl(ARM_ArgConvReverse_MatFrequency.GetString(cpnFreq).c_str());
		if(cpnDatas.size() >= 4)
			cpnIndexTermXl=cpnDatas[3];
		string cpnIndexTerm=CCSTringToSTLString(cpnIndexTermXl);

		/// RF index day count
		CCString cpnIndexDayCountXl("A360");
		if(cpnDatas.size() >= 5)
			cpnIndexDayCountXl=cpnDatas[4];
		long cpnIndexDayCount = ARM_ConvDayCount (cpnIndexDayCountXl);

		/// RF reset calendar
		CCString cpnResetCalXl("EUR");
		if(cpnDatas.size() >= 6)
			cpnResetCalXl=cpnDatas[5];
		string cpnResetCal = CCSTringToSTLString(cpnResetCalXl);

		/// RF payment calendar
		CCString cpnPayCalXl("EUR");
		if(cpnDatas.size() >= 7)
			cpnPayCalXl=cpnDatas[6];
		string cpnPayCal = CCSTringToSTLString(cpnPayCalXl);

		/// RF Stub Rule

		long stubRuleId;		
		CCString stubRuleXl("SS");
		if(cpnDatas.size() >= 8)
			stubRuleXl=cpnDatas[7];
		stubRuleId = ARM_ConvStubRule (stubRuleXl);

		/// RF reset gap
		double cpnResetGapD;
		double cpnResetGapDef=2;
		XL_readNumCellWD(XL_cpnResetGap,cpnResetGapD,cpnResetGapDef," ARM_ERR: coupon reset gap: numerical expected",C_result);
		long cpnResetGap = -(long)floor(cpnResetGapD);


		/// Fixed coupon datas
		/// ------------------
		/// Fixed leg end
		double fixEndDate;
		XL_readNumCellWD(XL_fixEndDate,fixEndDate,startDate," ARM_ERR: fixed leg end date: date expected",C_result);

		/// Fixed coupon day count
		CCString fixDayCountXl;
		XL_readStrCellWD(XL_fixDayCount,fixDayCountXl,cpnDayCountXl," ARM_ERR: fix cpn day count: string expected",C_result);
		long fixDayCount = ARM_ConvDayCount (fixDayCountXl);


		/// RF coupon profiles
		/// ------------------
		/// RF coupon index leverage
		double leverage;
		double leverageDef=1.0;
		CCString leverageStr;
		long     leverageId;
		XL_readStrOrNumCellWD(XL_leverage,leverageStr,leverage,leverageDef,leverageId,
			" ARM_ERR: strike: numeric or refValue Id expected",C_result);
		if(leverageId == XL_TYPE_STRING)
			leverageId = LocalGetNumObjectId(leverageStr);
		else
			leverageId = ARM_NULL_OBJECT;

		/// RF min coupon
		double cpnMin;
		double cpnMinDef=0.0;
		CCString cpnMinStr;
		long     cpnMinId;
		XL_readStrOrNumCellWD(XL_cpnMin,cpnMinStr,cpnMin,cpnMinDef,cpnMinId,
			" ARM_ERR: cpn min: numerical or refValue Id expected",C_result);
		if(cpnMinId == XL_TYPE_STRING)
			cpnMinId = LocalGetNumObjectId(cpnMinStr);
		else
			cpnMinId = ARM_NULL_OBJECT;

		/// RF max coupon
		double cpnMax;
		double cpnMaxDef=10000.0;
		CCString cpnMaxStr;
		long     cpnMaxId;
		XL_readStrOrNumCellWD(XL_cpnMax,cpnMaxStr,cpnMax,cpnMaxDef,cpnMaxId,
			" ARM_ERR: cpn max: numerical or refValue Id expected",C_result);
		if(cpnMaxId == XL_TYPE_STRING)
			cpnMaxId = LocalGetNumObjectId(cpnMaxStr);
		else
			cpnMaxId = ARM_NULL_OBJECT;


		/// Funding leg datas
		/// -----------------
		/// Funding coupon spread profile
		double fundSpread;
		double fundSpreadDef=0.0;
		CCString fundSpreadXl;
		long     fundSpreadId;
		XL_readStrOrNumCellWD(XL_fundSpread,fundSpreadXl,fundSpread,fundSpreadDef,fundSpreadId,
			" ARM_ERR: funding spread: numerical or refValue Id expected",C_result);
		if(fundSpreadId == XL_TYPE_STRING)
			fundSpreadId = LocalGetNumObjectId(fundSpreadXl);
		else
			fundSpreadId = ARM_NULL_OBJECT;

		VECTOR<CCString> fundDatas;
		XL_readStrVector (XL_fundDatas,fundDatas," ARM_ERR: Funding datas: array of string expected",DOUBLE_TYPE,C_result);

		/// Funding coupon frequency (reset & payment)
		CCString fundFreqXl(cpnFreqXl);
		if(fundDatas.size() >= 1)
			fundFreqXl=fundDatas[0];
		long fundFreq;
		if((fundFreq = ARM_ConvFrequency (fundFreqXl, C_result)) == ARM_DEFAULT_ERR)
		{
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}

		/// Funding coupon day count
		CCString fundDayCountXl("A360");
		if(fundDatas.size() >= 2)
			fundDayCountXl=fundDatas[1];
		long fundDayCount = ARM_ConvDayCount (fundDayCountXl);

		/// Funding & RF legs notional
		/// --------------------------
		double nominal;
		double nominalDef=100.0;
		CCString nominalStr;
		long     nominalId;
		XL_readStrOrNumCellWD(XL_nominal,nominalStr,nominal,nominalDef,nominalId,
			" ARM_ERR: nominal: numerical or refValue Id expected",C_result);
		if(nominalId == XL_TYPE_STRING)
			nominalId = LocalGetNumObjectId(nominalStr);
		else
			nominalId = ARM_NULL_OBJECT;

        double fundnominal;
		double fundnominalDef = ARM_MISSING_VALUE;
		CCString fundnominalStr;
		long     fundnominalId;
		XL_readStrOrNumCellWD(XL_FundNominal,fundnominalStr,fundnominal,fundnominalDef,fundnominalId,
			" ARM_ERR: nominal: numerical or refValue Id expected",C_result);
		if(fundnominalId == XL_TYPE_STRING)
			fundnominalId = LocalGetNumObjectId(fundnominalStr);
		else
			fundnominalId = ARM_NULL_OBJECT;



		/// Exercise datas
		/// --------------
		/// Notification lag
		double exerGapD;
		XL_readNumCellWD(XL_exerGap,exerGapD,cpnResetGapD," ARM_ERR: notification gap: numerical expected",C_result);
		long exerGap = -(long)floor(exerGapD);
		if (cpnTiming == K_ADVANCE && cpnResetGap < exerGap)
		{
			C_result.setMsg ("ARM_ERR: exercise must be anterior or equal to coupon reset");
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}

		/// Non callable period (w.r.t. RF coupon frequency)
		double nbNonCallD;
		double nbNonCallDDef=0.0;
		XL_readNumCellWD(XL_nbNonCall,nbNonCallD,nbNonCallDDef," ARM_ERR: nb of non call periods: numerical expected",C_result);
		long nbNonCall = (long)floor(nbNonCallD);

		/// Fees profile
		double exerFee;
		double exerFeeDef=0.0;
		CCString exerFeeStr;
		long     exerFeeId;
		XL_readStrOrNumCellWD(XL_exerFee,exerFeeStr,exerFee,exerFeeDef,exerFeeId,
			" ARM_ERR: exercise fees: numerical or refValue Id expected",C_result);
		if(exerFeeId == XL_TYPE_STRING)
			exerFeeId = LocalGetNumObjectId(exerFeeStr);
		else
			exerFeeId = ARM_NULL_OBJECT;

		/// use the concept of Functor to transfer the knowledge of
		/// a function with a context
		exportFunc34Args< double,
        double,
        double,
        long,
        long,
        double,
        long,
        long,
        long,
        long,
        string,
        long ,
        string,
        string,
		long,
        long,
        double,
        long,
        double,
        long,
        double,
        long,
        double,
        long,
        long,
        long,
        double,
        long,
        long,
        long,
        double,
        long,
        double,
        long > ourFunc (startDate, 
			endDate, strike, strikeId, payRec, cpnDayCount, cpnFreq, 
			fixEndDate, fixDayCount, cpnTiming, cpnIndexTerm, cpnIndexDayCount, cpnResetCal, cpnPayCal, 
			stubRuleId, cpnResetGap, leverage, leverageId, cpnMin, cpnMinId, cpnMax, cpnMaxId, fundSpread, 
			fundSpreadId, fundFreq,fundDayCount, nominal, nominalId, exerGap, nbNonCall, exerFee, exerFeeId, 
			fundnominal, fundnominalId, ARMLOCAL_Crude_CRFCalculator_Create);

				

		/// call the general function
		fillXL_Result( LOCAL_GC_CRF_CLASS, ourFunc, C_result, XL_result, PersistentInXL );
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_CRFCalculator_Common" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}
							 

///////////////////////////////////
/// version that takes into account 
/// previous creation of object
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_Crude_CRFCalculator_Create(
    LPXLOPER XL_startDate,
    LPXLOPER XL_endDate,
    LPXLOPER XL_strike,
    LPXLOPER XL_payRec,
    LPXLOPER XL_cpnDatas,
    LPXLOPER XL_fixEndDate,
    LPXLOPER XL_fixDayCount,
    LPXLOPER XL_cpnResetGap,
    LPXLOPER XL_leverage,
    LPXLOPER XL_cpnMin,
    LPXLOPER XL_cpnMax,
    LPXLOPER XL_fundSpread,
    LPXLOPER XL_fundDatas,
    LPXLOPER XL_nominal,
    LPXLOPER XL_exerGap,
    LPXLOPER XL_nbNonCall,
    LPXLOPER XL_exerFee,
    LPXLOPER XL_FundNominal)
{
	ADD_LOG("Local_Crude_CRFCalculator_Create");
	bool PersistentInXL = true;

	return Local_Crude_CRFCalculator_Common(
        XL_startDate,
        XL_endDate,
        XL_strike,
        XL_payRec,
        XL_cpnDatas,
        XL_fixEndDate,
        XL_fixDayCount,
        XL_cpnResetGap,
        XL_leverage,
        XL_cpnMin,
        XL_cpnMax,
        XL_fundSpread,
        XL_fundDatas,
        XL_nominal,
        XL_exerGap,
        XL_nbNonCall,
        XL_exerFee,
        XL_FundNominal,
		PersistentInXL );
}



///----------------------------------------------
///----------------------------------------------
///             CRF Calculator Accessor
/// Inputs :
///     
///----------------------------------------------
///----------------------------------------------
class crfGetDataFunc : public ARMResultLong2LongFunc
{
public:
	crfGetDataFunc(
        long crfId,
        const string& getType)
    :
    C_crfId(crfId),
    C_getType(getType)
    {};
	
	long operator()( ARM_result& result, long objId ){
		return ARMLOCAL_CRF_Get(
            C_crfId,
            C_getType,
            result,
            objId);
    }

private:
	long    C_crfId;
	string  C_getType;
};

class crfSetDataFunc : public ARMResultLong2LongFunc
{
public:
	crfSetDataFunc(
        long crfId,
        long dataToSetId,
        const string& setPortfolioType,
		const vector< string >& mktDataKeys,
        bool isUpdated)
    :
    C_crfId(crfId),
    C_dataToSetId(dataToSetId),
    C_setPortfolioType(setPortfolioType),
	C_mktDataKeys(mktDataKeys),
    C_isUpdated(isUpdated)
    {};
	
	long operator()( ARM_result& result, long objId ){
		return ARMLOCAL_CRF_Set(
            C_crfId,
            C_dataToSetId,
            C_setPortfolioType,
			C_mktDataKeys,
            C_isUpdated,
            result,
            objId);
	}

private:
	long    C_crfId;
	long    C_dataToSetId;
    string  C_setPortfolioType;
	vector< string > C_mktDataKeys;
    bool    C_isUpdated;
};


/////////////////////////////////////////////////////////////
/// central function that does the creation of the XL function
/////////////////////////////////////////////////////////////
LPXLOPER Local_CRFGet_Common(
	LPXLOPER XL_crfId,
	LPXLOPER XL_getType,
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

		long crfId;
		XL_GETOBJID( XL_crfId, crfId,	" ARM_ERR: CRF Calculator: Object expected",C_result);

		CCString getTypeStr;
		XL_readStrCell(XL_getType,getTypeStr," ARM_ERR: Accessor Type: string expected",C_result);
        char* type=getTypeStr.c_str(); // à cause du new !!
		string getType(type);
        delete type;

		CCString crfGetClass(GCGetTypeToClass(getType,crfId).c_str());

		crfGetDataFunc ourFunc(crfId,getType);

		/// call the general function
		fillXL_Result( crfGetClass, ourFunc, C_result, XL_result, PersistentInXL );
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_CRFGet_Common" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}
							 
LPXLOPER Local_CRFSet_Common(
	LPXLOPER XL_crfId,
	LPXLOPER XL_dataId,
	LPXLOPER XL_setPortfolioType,
	LPXLOPER XL_MktDataKeys,
	bool isUpdated,
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

		long crfId;
		XL_GETOBJID( XL_crfId, crfId,	" ARM_ERR: CRF Calculator: Object expected",C_result);

		long dataId;
		XL_GETOBJID( XL_dataId, dataId,	" ARM_ERR: Object expected",C_result);

		CCString setPortfolioTypeStr;
		XL_readStrCellWD(XL_setPortfolioType,setPortfolioTypeStr,"DEFAULT"," ARM_ERR: Portfolio Type: string expected",C_result);
        char * portType = setPortfolioTypeStr.c_str(); // à cause du new !!
		string setPortfolioType(portType);
        delete portType;

		VECTOR<CCString> mktDataKeys;
		VECTOR<CCString> mktDataKeysDef(0);
		XL_readStrVectorWD(XL_MktDataKeys,mktDataKeys,mktDataKeysDef," ARM_ERR: Market datas keys: array of string expected",DOUBLE_TYPE,C_result);

		vector< string > mktDataKeysSTL(mktDataKeys.size());
		for(size_t i=0;i<mktDataKeys.size();++i)
        {
            mktDataKeys[i].toUpper();
			mktDataKeysSTL[i]=CCSTringToSTLString(mktDataKeys[i]);
        }
        crfSetDataFunc ourFunc(crfId,
            dataId,
            setPortfolioType,
            mktDataKeysSTL,
            isUpdated);

        if(isUpdated)
        {
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
        else
        {
            /// General call through the functor with an object creation		    
		    fillXL_Result( LOCAL_GC_CRF_CLASS, ourFunc, C_result, XL_result, PersistentInXL );
        }
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_CRFSet_Common" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}
							 


///////////////////////////////////
/// version that takes into account 
/// previous creation of object
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_CRFCalculator_GetData(
	LPXLOPER XL_crfId,
	LPXLOPER XL_getType)
{
	ADD_LOG("Local_CRFCalculator_GetData");
	bool PersistentInXL = true;
	return Local_CRFGet_Common(
	    XL_crfId,
	    XL_getType,
        PersistentInXL );
}

__declspec(dllexport) LPXLOPER WINAPI Local_CRFCalculator_GetMRS(LPXLOPER XL_crfId)
{
	ADD_LOG("Local_CRFCalculator_GetMRS");
	static XLOPER XL_result;
	
	ARM_result C_result;
	
	ARM_XL_TRY_BLOCK_BEGIN
	{
		ARM_NOCALCIFWIZ();
		
		static int error;
		static char* reason = "";

		long crfId;
		XL_GETOBJID(XL_crfId,crfId," ARM_ERR: CRF calculator id: object expected",C_result);
	
	    long retCode;
	    long objId;
	    CCString prevClass;
	    
	    CCString curClass = LOCAL_REFVAL_CLASS;
	    CCString stringId = GetLastCurCellEnvValue ();
	    
	    if (!stringId)
	    {
			retCode = ARMLOCAL_CRF_GetMRS( crfId, C_result );

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
				retCode = ARMLOCAL_CRF_GetMRS( crfId, C_result, objId );

			    if(retCode == ARM_OK)
			    {			
				    LocalSetCurCellEnvValue (curClass, objId); 
				    stringId = LocalMakeObjectId (objId, curClass);
			    }
		    }
		    else
		    {
			    FreeCurCellContent ();

				retCode = ARMLOCAL_CRF_GetMRS( crfId, C_result );

			    if(retCode == ARM_OK)
			    {
				    objId = C_result.getLong ();
				    LocalSetCurCellEnvValue (curClass, objId); 
				    stringId = LocalMakeObjectId (objId, curClass);
			    }
		    }
	    }

		if (retCode == ARM_OK)
		{
			FreeCurCellErr ();
			XL_result.xltype = xltypeStr;
			XL_result.val.str = XL_StrC2StrPascal (stringId );
			XL_result.xltype |= xlbitDLLFree;
		}
		else
		{
			ARM_ERR();
		}
	}
	ARM_XL_TRY_BLOCK_END

	ARM_XL_CATCH_ARM_EXPT

	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_CRFCalculator_GetMRS" )

	return (LPXLOPER)&XL_result;
}

__declspec(dllexport) LPXLOPER WINAPI Local_CRFCalculator_SetData(
	LPXLOPER XL_crfId,
	LPXLOPER XL_dataId,
    LPXLOPER XL_setPortfolioType,
	LPXLOPER XL_mktDataKeys)
{
	ADD_LOG("Local_CRFCalculator_SetData");
	bool PersistentInXL = true;
    bool isUpdated = false;
	return Local_CRFSet_Common(
	    XL_crfId,
	    XL_dataId,
        XL_setPortfolioType,
		XL_mktDataKeys,
        isUpdated,
        PersistentInXL );
}


///////////////////////////////////
/// version for non persistentInXl
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_PXL_CRFCalculator_GetData(
	LPXLOPER XL_crfId,
	LPXLOPER XL_getType)
{
	ADD_LOG("Local_PXL_CRFCalculator_GetData");
	bool PersistentInXL = false;
	return Local_CRFGet_Common(
	    XL_crfId,
	    XL_getType,
        PersistentInXL );
}

__declspec(dllexport) LPXLOPER WINAPI Local_PXL_CRFCalculator_GetMRS(LPXLOPER XL_crfId)
{
	ADD_LOG("Local_PXL_CRFCalculator_GetMRS");
	static XLOPER XL_result;
	
	ARM_result C_result;
	
	ARM_XL_TRY_BLOCK_BEGIN
	{
		ARM_NOCALCIFWIZ();
		
		static int error;
		static char* reason = "";

		long crfId;
		XL_GETOBJID(XL_crfId,crfId," ARM_ERR: CRF calculator id: object expected",C_result);
	
	    long retCode;
	    long objId;
	    CCString prevClass;
	    
	    CCString curClass = LOCAL_REFVAL_CLASS;
	    CCString stringId = GetLastCurCellEnvValue ();
	    
		retCode = ARMLOCAL_CRF_GetMRS( crfId, C_result );

		if(retCode == ARM_OK)
		{
			objId = C_result.getLong ();

			stringId = LocalMakeObjectId (objId, curClass);
		}

		if (retCode == ARM_OK)
		{
			FreeCurCellErr ();
			XL_result.xltype = xltypeStr;
			XL_result.val.str = XL_StrC2StrPascal (stringId );
			XL_result.xltype |= xlbitDLLFree;
		}
		else
		{
			ARM_ERR();
		}
	}
	ARM_XL_TRY_BLOCK_END

	ARM_XL_CATCH_ARM_EXPT

	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_CRFCalculator_GetMRS" )

	return (LPXLOPER)&XL_result;
}

__declspec(dllexport) LPXLOPER WINAPI Local_PXL_CRFCalculator_SetData(
	LPXLOPER XL_crfId,
	LPXLOPER XL_dataId,
    LPXLOPER XL_setPortfolioType,
	LPXLOPER XL_mktDataKeys)
{
	ADD_LOG("Local_PXL_CRFCalculator_SetData");
	bool PersistentInXL = false;
    bool isUpdated = false;
	return Local_CRFSet_Common(
	    XL_crfId,
	    XL_dataId,
        XL_setPortfolioType,
		XL_mktDataKeys,
        isUpdated,
        PersistentInXL );
}

LPXLOPER Local_Initcalculator_Common(
	LPXLOPER XL_crfId,
	LPXLOPER XL_mktManagerId,
    LPXLOPER XL_toCalSigma,
    LPXLOPER XL_toCalMrs,
	LPXLOPER XL_strikeTypeToCalMRS,
    LPXLOPER XL_toAdjKcap,
    LPXLOPER XL_toAdjKfloor,
    LPXLOPER XL_pricingModel,
    LPXLOPER XL_toCalSkew,
    LPXLOPER XL_Kshift,
	LPXLOPER XL_KFrontier,
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

		long crfId;
		XL_GETOBJID( XL_crfId, crfId,	" ARM_ERR: CRF Calculator: Object expected",C_result);

		long mktManagerId;
		XL_GETOBJID( XL_mktManagerId, mktManagerId,	" ARM_ERR: Object expected",C_result);

        CCString C_toCalSigma;
		CCString defaultCalSigma = "YES";
		XL_readStrCellWD(XL_toCalSigma,C_toCalSigma,defaultCalSigma," ARM_ERR: To calibrate Sigma: String expected",C_result);

        CCString C_toCalMrs;
		CCString defaultCalMrs = "NO";
		XL_readStrCellWD(XL_toCalMrs,C_toCalMrs,defaultCalMrs," ARM_ERR: To calibrate Mean Reversion: String expected",C_result);

		CCString C_strikeTypeToCalMrs;
		CCString defaultstrikeTypeToCal = "STRIKE_EQUIVALENT";
		XL_readStrCellWD(XL_strikeTypeToCalMRS,C_strikeTypeToCalMrs,defaultstrikeTypeToCal," ARM_ERR: Strike Type to calibrate Mean Reversion: String expected",C_result);

        CCString C_AdjKcap;
		CCString defaultAdjKcap = "YES";
		XL_readStrCellWD(XL_toAdjKcap,C_AdjKcap,defaultAdjKcap," ARM_ERR: To adjust Caplet strikes: String expected",C_result);

        CCString C_AdjKfloor;
		CCString defaultAdjKfloor = "NO";
		XL_readStrCellWD(XL_toAdjKfloor,C_AdjKfloor,defaultAdjKfloor," ARM_ERR: To adjust floorlet strikes: String expected",C_result);

        CCString C_pricingModel;
		CCString defaultpricingModel = "HWM1F";
		XL_readStrCellWD(XL_pricingModel,C_pricingModel,defaultpricingModel," ARM_ERR: Model pricing name: String expected",C_result);

        CCString C_CalSkew;
		CCString defaultCalSkew = "YES";
		XL_readStrCellWD(XL_toCalSkew,C_CalSkew,defaultCalSkew," ARM_ERR: To calibrate skew: String expected",C_result);

        double C_Kshift;
		double defaultKshift=0.25;
		XL_readNumCellWD(XL_Kshift,C_Kshift,defaultKshift," ARM_ERR: percentage to shift: numerical expected",C_result);

		double C_Kfrontier;
		double defaultKfrontier=0.0;
		XL_readNumCellWD(XL_KFrontier,C_Kfrontier,defaultKfrontier," ARM_ERR: Frontier iter: numerical expected",C_result);
		long frontier = (long) C_Kfrontier;
        
        exportFunc11Args<long,long,CCString,CCString,CCString,CCString,
                    CCString,CCString,CCString,double, long>  
                    ourFunc(crfId,mktManagerId,C_toCalSigma,C_toCalMrs,C_strikeTypeToCalMrs,C_AdjKcap,
                    C_AdjKfloor,C_pricingModel,C_CalSkew,C_Kshift,frontier,ARMLOCAL_Calculator_Initialize);
		
        /// General call through the functor with an object creation		    
		fillXL_Result( LOCAL_GC_CRF_CLASS, ourFunc, C_result, XL_result, PersistentInXL );
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_Initcalculator_Common" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}

///////////////////////////////////
/// To initialise a crude calculator with 
///  a given Mkt Data Manager ( A PXL version is avalaible below)
//////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_CRFCalculator_Initialize(
	LPXLOPER XL_crfId,
	LPXLOPER XL_mktmanagerId,
    LPXLOPER XL_toCalSigma,
    LPXLOPER XL_toCalMrs,
	LPXLOPER XL_strikeTypeToCalMRS,
    LPXLOPER XL_toAdjKcap,
    LPXLOPER XL_toAdjKfloor,
    LPXLOPER XL_pricingModel,
    LPXLOPER XL_toCalSkew,
    LPXLOPER XL_Kshift,
	LPXLOPER XL_KFrontier
    )
{
	ADD_LOG("Local_CRFCalculator_Initialize");
	bool PersistentInXL = false;
	return Local_Initcalculator_Common(
	    XL_crfId,
	    XL_mktmanagerId,
        XL_toCalSigma,
        XL_toCalMrs,
		XL_strikeTypeToCalMRS,
        XL_toAdjKcap,
        XL_toAdjKfloor,
        XL_pricingModel,
        XL_toCalSkew,
        XL_Kshift,
		XL_KFrontier,
        PersistentInXL );
}

__declspec(dllexport) LPXLOPER WINAPI PXL_Local_CRFCalculator_Initialize(
	LPXLOPER XL_crfId,
	LPXLOPER XL_mktmanagerId,
    LPXLOPER XL_toCalSigma,
    LPXLOPER XL_toCalMrs,
	LPXLOPER XL_strikeTypeToCalMRS,
    LPXLOPER XL_toAdjKcap,
    LPXLOPER XL_toAdjKfloor,
    LPXLOPER XL_pricingModel,
    LPXLOPER XL_toCalSkew,
    LPXLOPER XL_Kshift,
	LPXLOPER XL_KFrontier
    )
{
	ADD_LOG("PXL_Local_CRFCalculator_Initialize");
	bool PersistentInXL = false;
	return Local_Initcalculator_Common(
	    XL_crfId,
	    XL_mktmanagerId,
        XL_toCalSigma,
        XL_toCalMrs,
		XL_strikeTypeToCalMRS,
        XL_toAdjKcap,
        XL_toAdjKfloor,
        XL_pricingModel,
        XL_toCalSkew,
        XL_Kshift,
		XL_KFrontier,
        PersistentInXL );
}

///////////////////////////////////
/// same as SetData but the previous CRF
/// is not cloned simply updated (=> no PXL version)
/// Only GC datas are updattable
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_CRFCalculator_Update(
	LPXLOPER XL_gcId,
	LPXLOPER XL_dataId,
    LPXLOPER XL_setPortfolioType,
	LPXLOPER XL_mktDataKeys)
{
	ADD_LOG("Local_CRFCalculator_Update");
	bool PersistentInXL = true;
    bool isUpdated = true;
	return Local_CRFSet_Common(
	    XL_gcId,
	    XL_dataId,
        XL_setPortfolioType,
		XL_mktDataKeys,
        isUpdated,
        PersistentInXL );
}

/////////////////////////////////////////////////////////////
/// central function that does the creation of the XL function
/////////////////////////////////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_CRFCalculator_SetAutoCalFlags(
	LPXLOPER XL_crfId,
	LPXLOPER XL_flags )
{
	ADD_LOG("Local_CRFCalculator_SetAutoCalFlags");
	
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

		long crfId;
		XL_GETOBJID(XL_crfId,crfId," ARM_ERR: CRF calculator id: object expected",C_result);

		VECTOR<CCString> C_flags;
		VECTOR<CCString> flagsDef(4,"Y");
        flagsDef[3]=CCString("N"); // floor calib
		XL_readStrVectorWD(XL_flags,C_flags,flagsDef," ARM_ERR: flags: array of string expected",DOUBLE_TYPE,C_result);
		vector< string > flags(C_flags.size());
		for(size_t i=0;i<C_flags.size();++i)
        {
            (C_flags[i]).toUpper();
			flags[i]=CCSTringToSTLString(C_flags[i]);
        }
			
		long retCode = ARMLOCAL_CRF_SetAutoCalFlags(
				crfId,
				flags,
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_CRFCalculator_SetAutoCalFlags" )

	
	return (LPXLOPER)&XL_result;
}

/////////////////////////////////////////////////////////////
/// central function that does the creation of the XL function
/////////////////////////////////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_CRFCalculator_SetOneCalFlag(
	LPXLOPER XL_crfId,
	LPXLOPER XL_flag ,
	LPXLOPER XL_type)
{
	ADD_LOG("Local_CRFCalculator_SetOneCalFlag");
	
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

		long crfId;
		XL_GETOBJID(XL_crfId,crfId," ARM_ERR: CRF calculator id: object expected",C_result);

		CCString C_flag;
		XL_readStrCellWD(XL_flag,C_flag,"NO"," ARM_ERR: Faly string Type: string expected",C_result);

		CCString C_type;
		XL_readStrCellWD(XL_type,C_type,"SKEW"," ARM_ERR: Faly string Type: string expected",C_result);

		C_flag.toUpper();
		string falgStr=CCSTringToSTLString(C_flag);
		C_type.toUpper();
		string typeStr=CCSTringToSTLString(C_type);

		exportFunc3Args< long, string, string > ourFunc(crfId,falgStr,typeStr, ARMLOCAL_CRF_SetOneCalFlag);

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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_CRFCalculator_SetOneCalFlag" )

	
	return (LPXLOPER)&XL_result;
}


/////////////////////////////////////////////////////////////
/// central function that does the creation of the XL function
/////////////////////////////////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_CRFCalculator_SetAutoCalFlagsAndClone_Common(
	LPXLOPER XL_crfId,
	LPXLOPER XL_flags,
	bool PersistentInXL)
{
	ADD_LOG("Local_CRFCalculator_SetAutoCalFlagsAndClone_Common");
	
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

		long crfId;
		XL_GETOBJID(XL_crfId,crfId," ARM_ERR: CRF calculator id: object expected",C_result);

		VECTOR<CCString> C_flags;
		VECTOR<CCString> flagsDef(4,"Y");
        flagsDef[3]=CCString("N"); // floor calib
		XL_readStrVectorWD(XL_flags,C_flags,flagsDef," ARM_ERR: flags: array of string expected",DOUBLE_TYPE,C_result);
		vector< string > flags(C_flags.size());
		for(size_t i=0;i<C_flags.size();++i)
        {
            (C_flags[i]).toUpper();
			flags[i]=CCSTringToSTLString(C_flags[i]);
        }
			
		exportFunc2Args< long, vector < string > > ourFunc(crfId,flags,ARMLOCAL_CRF_SetAutoCalFlagsAndClone);

		/// call the general function
		fillXL_Result(LOCAL_GC_CRF_CLASS, ourFunc, C_result, XL_result, PersistentInXL );
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_CRFCalculator_SetAutoCalFlagsAndClone" )

	
	return (LPXLOPER)&XL_result;
}

///////////////////////////////////
/// version that takes into account 
/// previous creation of object
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_CRFCalculator_SetAutoCalFlagsAndClone(
	LPXLOPER XL_crfId,
	LPXLOPER XL_flags )
{
	ADD_LOG("Local_CRFCalculator_SetAutoCalFlagsAndClone");
	bool PersistentInXL = true;
	return Local_CRFCalculator_SetAutoCalFlagsAndClone_Common(
		XL_crfId,
		XL_flags,
        PersistentInXL );
}

__declspec(dllexport) LPXLOPER WINAPI Local_PXL_CRFCalculator_SetAutoCalFlagsAndClone(
	LPXLOPER XL_crfId,
	LPXLOPER XL_flags )
{
	ADD_LOG("Local_PXL_CRFCalculator_SetAutoCalFlagsAndClone");
	bool PersistentInXL = false;
	return Local_CRFCalculator_SetAutoCalFlagsAndClone_Common(
		XL_crfId,
		XL_flags,
        PersistentInXL );
}

/////////////////////////////////////////////////////////////
/// central function that does the creation of the XL function
/////////////////////////////////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_CRFCalculator_SetProductToPrice(
	LPXLOPER XL_crfId,
	LPXLOPER XL_ProdToPrice )
{
	ADD_LOG("Local_CRFCalculator_SetProductToPrice");
	
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

		long crfId;
		CCString C_prodToPrice;
		XL_GETOBJID(XL_crfId,crfId," ARM_ERR: CRF calculator id: object expected",C_result);
		XL_readStrCell(XL_ProdToPrice,C_prodToPrice," ARM_ERR: flags: string expected",C_result);
		string prodToPrice = CCSTringToSTLString(C_prodToPrice);
			
		long retCode = ARMLOCAL_CRF_SetProductToPrice(
				crfId,
				prodToPrice,
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_CRFCalculator_SetProductToPrice" )

	
	return (LPXLOPER)&XL_result;
}


///----------------------------------------------
///----------------------------------------------
///             Bermuda Swaption Calculator   
///----------------------------------------------
///----------------------------------------------
class bermudaswaptionCalculatorFunc : public ARMResultLong2LongFunc
{
public:
	bermudaswaptionCalculatorFunc(
		ARM_Currency ccy,
	    double startDate,
        double endDate,
        double notional,
		long notionalId,
		double strike,
		long strikeId,
		double fees,
		long feesId,
		double spread,
		long spreadId,
		int payReceive,
		int callFreq,
		int callNotice,
		string callCal,
		double firstCallDate,
		double lastCallDate,
		int fixFreq,
		int fixBasis,
		string fixPayCal,
		int fixAdjRule,
		int fixRule,
		int fixPayGap,
		bool isZc,
		int varFreq,
		int varBasis,
		string varResetCal,
		string varPayCal,
		string varIndexTerm,
		int varAdjRule,
		int varRule,
		int varResetGap,
		int varPayGap,
		int stubRule,
		int genSecType,
		vector<int>* controlVariates,
		vector<double>* controlPrices,
		vector<string> mdmKeys,
		long mktDataManagerId,
		int modelType,
		vector<string> modelParams,
		vector<string> calibFlags,
		int numMethodType,
		int amcIter,
		int mcIter,
		int maxBucketSize,
		string genType1,
		string genType2,
		string pathOrder,
		string pathScheme,
		int firstNbTimes,
		int treeSteps,
		vector<int> portfolioMode,
		bool boundaryFlag,
		bool approxMarginFlag,
		bool freezeBetasFlag,
		bool calculateProbaFlag)
    :
	C_currency(ccy),
    C_startDate(startDate),
    C_endDate(endDate),
    C_notional(notional),
	C_notionalId(notionalId),
	C_strike(strike),
    C_strikeId(strikeId),
	C_fees(fees),
    C_feesId(feesId),
	C_spread(spread),
    C_spreadId(spreadId),
    C_payReceive(payReceive),
	C_callFreq(callFreq),
	C_callNotice(callNotice),
	C_callCal(callCal),
	C_firstCallDate(firstCallDate),
	C_lastCallDate(lastCallDate),
	C_fixFreq(fixFreq),
	C_fixBasis(fixBasis),
	C_fixPayCal(fixPayCal),
	C_fixAdjRule(fixAdjRule),
	C_fixRule(fixRule),
	C_fixPayGap(fixPayGap),
	C_isZc(isZc),
	C_varFreq(varFreq),
	C_varBasis(varBasis),
	C_varResetCal(varResetCal),
	C_varPayCal(varPayCal),
	C_varIndexTerm(varIndexTerm),
	C_varAdjRule(varAdjRule),
	C_varRule(varRule),
	C_varResetGap(varResetGap),
	C_varPayGap(varPayGap),
	C_stubRule(stubRule),
	C_genSecType(genSecType),
	C_controlVariates(controlVariates),
	C_controlPrices(controlPrices),
	C_mdmKeys(mdmKeys),
	C_mktDataManagerId(mktDataManagerId),
	C_modelType(modelType),
	C_modelParams(modelParams),
	C_calibFlags(calibFlags),
	C_NumMethodType(numMethodType),
	C_AmcIter(amcIter),
	C_McIter(mcIter),
	C_MaxBucketSize(maxBucketSize),
	C_GenType1(genType1),
	C_GenType2(genType2),
	C_PathOrder(pathOrder),
	C_PathScheme(pathScheme),
	C_FirstNbTimes(firstNbTimes),
	C_TreeSteps(treeSteps),
	C_portfolioMode(portfolioMode),
	C_boundaryFlag(boundaryFlag),
	C_approxMarginFlag(approxMarginFlag),
	C_freezeBetasFlag(freezeBetasFlag),
	C_calculateProbaFlag(calculateProbaFlag)
   
    {};
	
	long operator()( ARM_result& result, long objId ){
		return ARMLOCAL_BERMUDASWAPTIONCalculator_Create(
			C_currency,
	  		C_startDate,
			C_endDate,
			C_notional,
			C_notionalId,
			C_strike,
			C_strikeId,
			C_fees,
			C_feesId,
			C_spread,
			C_spreadId,
			C_payReceive,
			C_callFreq,
			C_callNotice,
			C_callCal,
			C_firstCallDate,
			C_lastCallDate,
			C_fixFreq,
			C_fixBasis,
			C_fixPayCal,
			C_fixAdjRule,
			C_fixRule,
			C_fixPayGap,
			C_isZc,
			C_varFreq,
			C_varBasis,
			C_varResetCal,
			C_varPayCal,
			C_varIndexTerm,
			C_varAdjRule,
			C_varRule,
			C_varResetGap,
			C_varPayGap,
			C_stubRule,
			C_genSecType,
			C_controlVariates,
			C_controlPrices,
			C_mdmKeys,
			C_mktDataManagerId,
			C_modelType,
			C_modelParams,
			C_calibFlags,
			C_NumMethodType,
			C_AmcIter,
			C_McIter,
			C_MaxBucketSize,
			C_GenType1,
			C_GenType2,
			C_PathScheme,
			C_PathOrder,
			C_FirstNbTimes,
			C_TreeSteps,
			C_portfolioMode,
			C_boundaryFlag,
			C_approxMarginFlag,
			C_freezeBetasFlag,
			C_calculateProbaFlag,
            result,
            objId);
    }

private:
		ARM_Currency		C_currency;
		double				C_startDate;
        double				C_endDate;
        double				C_notional;
		long				C_notionalId;
		double				C_strike;
		long				C_strikeId;
		double				C_fees;
		long				C_feesId;
		double				C_spread;
		long				C_spreadId;
		int					C_payReceive;
		int					C_callFreq;
		int					C_callNotice;
		string				C_callCal;
		double				C_firstCallDate;
		double				C_lastCallDate;
		int					C_fixFreq;
		int					C_fixBasis;
		string				C_fixPayCal;
		int					C_fixAdjRule;
		int					C_fixRule;
		int					C_fixPayGap;
		bool				C_isZc;
		int					C_varFreq;
		int					C_varBasis;
		string				C_varResetCal;
		string				C_varPayCal;
		string				C_varIndexTerm;
		int					C_varAdjRule;
		int					C_varRule;
		int					C_varResetGap;
		int					C_varPayGap;
		int					C_stubRule;
		int					C_genSecType;
		vector<int>*		C_controlVariates;
		vector<double>*		C_controlPrices;
		vector<string>		C_mdmKeys;
		long				C_mktDataManagerId;
		vector<string>		C_modelParams;
		vector<string>		C_calibFlags;
		int					C_NumMethodType;
		int					C_AmcIter;
		int					C_McIter;
		int					C_MaxBucketSize;
		string				C_GenType1;
		string				C_GenType2;
		string				C_PathOrder;
		string				C_PathScheme;
		int					C_FirstNbTimes;
		int					C_TreeSteps;
		vector<int>			C_portfolioMode;
		bool				C_boundaryFlag;
		bool				C_approxMarginFlag;
		bool				C_freezeBetasFlag;
		int					C_modelType;
		bool				C_calculateProbaFlag;
};


class bermudaswaptionGetDataFunc : public ARMResultLong2LongFunc
{
public:
	bermudaswaptionGetDataFunc(
        long bsId,
        const string& getType)
    :
    C_bsId(bsId),
    C_getType(getType)
    {};
	
	long operator()( ARM_result& result, long objId ){
		return ARMLOCAL_BERMUDASWAPTION_Get(
            C_bsId,
            C_getType,
            result,
            objId);
    }

private:
	long    C_bsId;
	string  C_getType;
};

class bermudaswaptionSetDataFunc : public ARMResultLong2LongFunc
{
public:
	bermudaswaptionSetDataFunc(
        long bsId,
		long dataToSetId,
        const string& setPortfolioType,
		const vector< string >& mktDataKeys,
        bool isUpdated)
    :
    C_bsId(bsId),
    C_dataToSetId(dataToSetId),
    C_setPortfolioType(setPortfolioType),
	C_mktDataKeys(mktDataKeys),
    C_isUpdated(isUpdated)
    {};
	
	long operator()( ARM_result& result, long objId ){
		return ARMLOCAL_BERMUDASWAPTION_Set(
            C_bsId,
            C_dataToSetId,
            C_setPortfolioType,
			C_mktDataKeys,
            C_isUpdated,
            result,
            objId);
	}

private:
	long    C_bsId;
	long    C_dataToSetId;
    string  C_setPortfolioType;
	vector< string > C_mktDataKeys;
    bool    C_isUpdated;
};



LPXLOPER Local_BermudaSwaptionGet_Common(
	LPXLOPER XL_bsId,
	LPXLOPER XL_getType,
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

		long bsId;
		XL_GETOBJID( XL_bsId, bsId,	" ARM_ERR: Bermuda Swaption Calculator: Object expected",C_result);

		CCString getTypeStr;
		XL_readStrCell(XL_getType,getTypeStr," ARM_ERR: Accessor Type: string expected",C_result);
        char* type=getTypeStr.c_str(); // à cause du new !!
		string getType(type);
        delete type;

		CCString bsGetClass(GCGetTypeToClass(getType,bsId).c_str());

		bermudaswaptionGetDataFunc ourFunc(bsId,getType);

		/// call the general function
		fillXL_Result( bsGetClass, ourFunc, C_result, XL_result, PersistentInXL );
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_BermudaSwaptionGet_Common" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}
							 
LPXLOPER Local_BermudaSwaptionSet_Common(
	LPXLOPER XL_bsId,
	LPXLOPER XL_dataId,
	LPXLOPER XL_setPortfolioType,
	LPXLOPER XL_MktDataKeys,
	bool isUpdated,
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

		long bsId;
		XL_GETOBJID( XL_bsId, bsId,	" ARM_ERR: Bermuda Swaption Calculator: Object expected",C_result);

		long dataId;
		XL_GETOBJID( XL_dataId, dataId,	" ARM_ERR: Object expected",C_result);

		CCString setPortfolioTypeStr;
		XL_readStrCellWD(XL_setPortfolioType,setPortfolioTypeStr,"DEFAULT"," ARM_ERR: Portfolio Type: string expected",C_result);
        char * portType = setPortfolioTypeStr.c_str(); // à cause du new !!
		string setPortfolioType(portType);
        delete portType;

		VECTOR<CCString> mktDataKeys;
		VECTOR<CCString> mktDataKeysDef(0);
		XL_readStrVectorWD(XL_MktDataKeys,mktDataKeys,mktDataKeysDef," ARM_ERR: Market datas keys: array of string expected",DOUBLE_TYPE,C_result);

		vector< string > mktDataKeysSTL(mktDataKeys.size());
		for(size_t i=0;i<mktDataKeys.size();++i)
        {
            mktDataKeys[i].toUpper();
			mktDataKeysSTL[i]=CCSTringToSTLString(mktDataKeys[i]);
        }

        bermudaswaptionSetDataFunc ourFunc(bsId,
										   dataId,
										   setPortfolioType,
           								   mktDataKeysSTL,
										   isUpdated);

		if(isUpdated)
        {
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
        else
        {
            /// General call through the functor with an object creation		    
		    fillXL_Result( LOCAL_GC_BERMUDASWAPTION_CLASS, ourFunc, C_result, XL_result, PersistentInXL );
        }
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_BermudaSwaptionSet_Common" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}



__declspec(dllexport) LPXLOPER WINAPI Local_BermudaSwaptionCalculator_GetData(
	LPXLOPER XL_bsId,
	LPXLOPER XL_getType)
{
	ADD_LOG("Local_BermudaSwaptionCalculator_GetData");
	bool PersistentInXL = true;
	return Local_BermudaSwaptionGet_Common(
	    XL_bsId,
	    XL_getType,
        PersistentInXL );
}

__declspec(dllexport) LPXLOPER WINAPI Local_BermudaSwaptionCalculator_SetData(
	LPXLOPER XL_bsId,
	LPXLOPER XL_dataId,
	LPXLOPER XL_setPortfolioType,
	LPXLOPER XL_mktDataKeys)
{
	ADD_LOG("Local_BermudaSwaptionCalculator_SetData");
	bool PersistentInXL = true;
	bool isUpdated = false;
	return Local_BermudaSwaptionSet_Common(
	    XL_bsId,
	    XL_dataId,
		XL_setPortfolioType,
		XL_mktDataKeys,
		isUpdated,
        PersistentInXL );
}


///////////////////////////////////
/// version for non persistentInXl
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_PXL_BermudaSwaptionCalculator_GetData(
	LPXLOPER XL_bsId,
	LPXLOPER XL_getType)
{
	ADD_LOG("Local_PXL_BermudaSwaptionCalculator_GetData");
	bool PersistentInXL = false;
	return Local_BermudaSwaptionGet_Common(
	    XL_bsId,
	    XL_getType,
        PersistentInXL );
}

__declspec(dllexport) LPXLOPER WINAPI Local_PXL_BermudaSwaptionCalculator_SetData(
	LPXLOPER XL_bsId,
	LPXLOPER XL_dataId,
	LPXLOPER XL_setPortfolioType,
	LPXLOPER XL_mktDataKeys)
{
	ADD_LOG("Local_PXL_BermudaSwaptionCalculator_SetData");
	bool PersistentInXL = false;
	bool isUpdated = false;
	return Local_BermudaSwaptionSet_Common(
	    XL_bsId,
	    XL_dataId,
		XL_setPortfolioType,
		XL_mktDataKeys,
		isUpdated,
        PersistentInXL );
}



LPXLOPER Local_BermudaSwaptionRootMrs_Common(
	LPXLOPER XL_bsId,
	LPXLOPER XL_targetPrice,
	LPXLOPER XL_fTolerance,
	LPXLOPER XL_maxIter,
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

		long bsId;
		XL_GETOBJID( XL_bsId, bsId,	"ARM_ERR: Bermuda Swaption Calculator: Object expected",C_result);

		double C_targetPrice;
		XL_readNumCell(XL_targetPrice, C_targetPrice, "ARM_ERR: targetPrice: double expected", C_result);

		double C_fTolerance;
		XL_readNumCell(XL_fTolerance, C_fTolerance, "ARM_ERR: F or X tolerance: double expected", C_result);

		double C_maxIter;
		XL_readNumCell(XL_maxIter, C_maxIter, "ARM_ERR: maximum iterations: double expected", C_result);

		long retCode = ARMLOCAL_ARM_BermudaRootMrs (bsId, C_targetPrice, C_fTolerance, (int)(C_maxIter), C_result);

		if(retCode == ARM_OK)
		{
			FreeCurCellErr ();
			XL_result.xltype  = xltypeNum;
			XL_result.val.num = C_result.getDouble ();
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_BermudaSwaptionRootMrs_Common" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


__declspec(dllexport) LPXLOPER WINAPI Local_BermudaSwaptionCalculator_RootMrs(
	LPXLOPER XL_bsId,
	LPXLOPER XL_targetPrice,
	LPXLOPER XL_fTolerance,
	LPXLOPER XL_maxIter)
{
	ADD_LOG("Local_BermudaSwaptionCalculator_RootMrs");
	bool PersistentInXL = true;
	return Local_BermudaSwaptionRootMrs_Common(
	    XL_bsId,
	    XL_targetPrice,
		XL_fTolerance,
		XL_maxIter,
        PersistentInXL );
}

__declspec(dllexport) LPXLOPER WINAPI Local_PXL_BermudaSwaptionCalculator_RootMrs(
	LPXLOPER XL_bsId,
	LPXLOPER XL_targetPrice,
	LPXLOPER XL_fTolerance,
	LPXLOPER XL_maxIter)
{
	ADD_LOG("Local_PXL_BermudaSwaptionCalculator_RootMrs");
	bool PersistentInXL = false;
	return Local_BermudaSwaptionRootMrs_Common(
	    XL_bsId,
	    XL_targetPrice,
		XL_fTolerance,
		XL_maxIter,
        PersistentInXL );
}

/////////////////////////////////////////////////////////////
/// central function that does the creation of the XL function
/////////////////////////////////////////////////////////////
LPXLOPER Local_BermudaSwaptionCalculator_Common(LPXLOPER XL_currency,
												LPXLOPER XL_startDate,
												LPXLOPER XL_endDate,
												LPXLOPER XL_notionalCurve,
												LPXLOPER XL_strikeCurve,
												LPXLOPER XL_feesCurve,
												LPXLOPER XL_spreadCurve,
												LPXLOPER XL_payReceive,
												LPXLOPER XL_callDatas,
												LPXLOPER XL_fixDatas,
												LPXLOPER XL_varDatas,
												LPXLOPER XL_genSecType,
												LPXLOPER XL_mktDataManager,
												LPXLOPER XL_calibParams,
												LPXLOPER XL_controlVariates,
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

		//Currency
		CCString ccyStr;
		XL_readStrCell(XL_currency,ccyStr," ARM_ERR: currency: string expected",C_result);
		string sCcy = CCSTringToSTLString(ccyStr);
		ARM_Currency ccy(sCcy.c_str());

		/// Start & End
		double startDate;
		XL_readNumCell(XL_startDate,startDate," ARM_ERR: start date: date expected",C_result);
		double endDate;
		XL_readNumCell(XL_endDate,endDate," ARM_ERR: end date: date expected",C_result);

		//NotionalValue or NotionalCurve
		double notional;
		CCString notionalStr;
		long     notionalId;
		XL_readStrOrNumCell(XL_notionalCurve, notionalStr, notional, notionalId,
			   " ARM_ERR: notional: numeric or refValue Id expected",C_result);	
		if(notionalId == XL_TYPE_STRING)
			notionalId = LocalGetNumObjectId(notionalStr);
		else
			notionalId = ARM_NULL_OBJECT;

		/// StrikeValue or StrikeCurve
		double strike;
		CCString strikeStr;
		long     strikeId;
		XL_readStrOrNumCell(XL_strikeCurve, strikeStr, strike, strikeId,
			   " ARM_ERR: strike: numeric or refValue Id expected",C_result);	
		if(strikeId == XL_TYPE_STRING)
			strikeId = LocalGetNumObjectId(strikeStr);
		else
			strikeId = ARM_NULL_OBJECT;

		/// FeesValue or FeesCurve
		double fees;
		CCString feesStr;
		long     feesId;
		XL_readStrOrNumCell(XL_feesCurve, feesStr, fees, feesId,
			   " ARM_ERR: fees: numeric or refValue Id expected",C_result);	
		if(feesId == XL_TYPE_STRING)
			feesId = LocalGetNumObjectId(feesStr);
		else
			feesId = ARM_NULL_OBJECT;
	
		/// SpreadValue or SpreadCurve
		double spread;
		CCString spreadStr;
		long     spreadId;
		XL_readStrOrNumCell(XL_spreadCurve, spreadStr, spread, spreadId,
			   " ARM_ERR: spread: numeric or refValue Id expected",C_result);	
		if(spreadId == XL_TYPE_STRING)
			spreadId = LocalGetNumObjectId(spreadStr);
		else
			spreadId = ARM_NULL_OBJECT;

		/// Pay/Receive
		CCString payReceiveStr;
		long lPayReceive;
		XL_readStrCell(XL_payReceive,payReceiveStr," ARM_ERR: payer/receiver: string expected",C_result);
		if((lPayReceive = ARM_ConvRecOrPay (payReceiveStr, C_result)) == ARM_DEFAULT_ERR)
		{
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}
		int payReceive = (int)lPayReceive;

		/// Call Datas
		/// ----------
		VECTOR<CCString> callDatas;
		XL_readStrVector (XL_callDatas,callDatas," ARM_ERR: Call datas: array of string expected",DOUBLE_TYPE,C_result);
		if(callDatas.size() < 5)
		{
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}

		/// Call frequency 
		CCString callFreqXl=callDatas[0];
		long lCallFreq;
		if((lCallFreq = ARM_ConvFrequency (callFreqXl, C_result)) == ARM_DEFAULT_ERR)
		{
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}
		int callFreq = (int)lCallFreq;

		//CallNotice
		CCString callNoticeXl = callDatas[1];
		string sCallNotice = CCSTringToSTLString(callNoticeXl);
		int callNotice = atoi(sCallNotice.c_str());
		
		//CallCal
		CCString callCalXl = callDatas[2];
		string callCal = CCSTringToSTLString(callCalXl);
	
		//FirstCallDate
		CCString firstCallDateXl = callDatas[3];
		string sFirstCallDate = CCSTringToSTLString(firstCallDateXl);
		double firstCallDate = atoi(sFirstCallDate.c_str());
	
		//LastCallDate
		CCString lastCallDateXl = callDatas[4];
		string sLastCallDate = CCSTringToSTLString(lastCallDateXl);
		double lastCallDate = atoi(sLastCallDate.c_str());

		/// Underlying Fix Leg Datas
		/// ------------------------
		VECTOR<CCString> fixDatas;
		XL_readStrVector (XL_fixDatas, fixDatas," ARM_ERR: Fix datas: array of string expected",DOUBLE_TYPE,C_result);
		if(fixDatas.size() < 7)
		{
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}

		//Fix Freq
		CCString fixFreqXl = fixDatas[0];
		long lFixFreq;
		if((lFixFreq = ARM_ConvFrequency (fixFreqXl, C_result)) == ARM_DEFAULT_ERR)
		{
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}
		int fixFreq = (int)lFixFreq;

		//Fix Basis
		CCString fixBasisXl = fixDatas[1];
		int fixBasis = (int)ARM_ConvDayCount (fixBasisXl);

		//Fix Pay Cal
		CCString fixPayCalXl = fixDatas[2];
		string 	fixPayCal = CCSTringToSTLString(fixPayCalXl);

		//Fix Adj Rule
		CCString fixAdjRuleXl = fixDatas[3];
		int fixAdjRule = (int)ARM_ConvFwdRule(fixAdjRuleXl);

		//Fix Adj Rule
		CCString fixRuleXl = fixDatas[4];
		int fixRule = (int)ARM_ConvIntRule(fixRuleXl);

		
		//Fix Pay Gap
		CCString fixPayGapXl = fixDatas[5];
		string sFixPayGap = CCSTringToSTLString(fixPayGapXl);
		int fixPayGap = atoi(sFixPayGap.c_str());

		//Is ZC
		CCString isZcXl = fixDatas[6];
		string sisZcXl = CCSTringToSTLString(isZcXl);
		string sisZc(stringGetUpper(sisZcXl));
		bool isZc; 
		if(sisZc == "Y")
		{
			isZc = true;
		}
		else if (sisZc == "N")
		{
			isZc = false;
		}

		/// Underlying Var Leg Datas
		/// ------------------------
		VECTOR<CCString> varDatas;
		XL_readStrVector (XL_varDatas, varDatas," ARM_ERR: Var datas: array of string expected",DOUBLE_TYPE,C_result);
		if(varDatas.size() < 9)
		{
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}

		//Var Freq
		CCString varFreqXl = varDatas[0];
		long lVarFreq;
		if((lVarFreq = ARM_ConvFrequency (varFreqXl, C_result)) == ARM_DEFAULT_ERR)
		{
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}
		int varFreq = (int)lVarFreq;

		//Var Basis
		CCString varBasisXl = varDatas[1];
		int varBasis = (int)ARM_ConvDayCount (varBasisXl);

		//Var Reset Cal
		CCString varResetCalXl = varDatas[2];
		string 	varResetCal = CCSTringToSTLString(varResetCalXl);

		//Var Pay Cal
		CCString varPayCalXl = varDatas[3];
		string 	varPayCal = CCSTringToSTLString(varPayCalXl);

		//Var Index Term
		CCString varIndexTermXl = varDatas[4];
		string 	varIndexTerm = CCSTringToSTLString(varIndexTermXl);

		//Var Adj Rule
		CCString varAdjRuleXl = varDatas[5];
		int varAdjRule = (int)ARM_ConvFwdRule(varAdjRuleXl);

		
		//Var  Rule
		CCString varRuleXl = varDatas[6];
		int varRule = (int)ARM_ConvIntRule(varRuleXl);


		//Var Reset Gap
		CCString varResetGapXl = varDatas[7];
		string sVarResetGap = CCSTringToSTLString(varResetGapXl);
		int varResetGap = atoi(sVarResetGap.c_str());

		//Var Pay Gap
		CCString varPayGapXl = varDatas[8];
		string sVarPayGap = CCSTringToSTLString(varPayGapXl);
		int varPayGap = atoi(sVarPayGap.c_str());

		//Underlying Swap Stub Rule
		CCString stubRuleXl = varDatas[9];
		int stubRule;		
		stubRule = (int)ARM_ConvStubRule(stubRuleXl);

		//Generic Security Type
		double dgenSecType;
		int genSecType;
		XL_readNumCell(XL_genSecType,dgenSecType," ARM_ERR: gen sec type: integer expected",C_result);
		genSecType = (int)dgenSecType;

		//ControlVariates
		VECTOR<CCString> controlVariates;
		VECTOR<CCString> controlVariatesDef(0);
		XL_readStrVectorWD(XL_controlVariates, controlVariates,controlVariatesDef, " ARM_ERR: control variates: array of string expected",DOUBLE_TYPE,C_result);
		int ctrlSize = controlVariates.size();
		vector<int>* ctrlVariateVec = new vector<int>(ctrlSize);
		vector<double>* ctrlVariatePrices = new vector<double>(ctrlSize);
		CCString cCurrCtrl;
		string sCurrCtrl;
		size_t i;
		for (i=0; i<ctrlSize;i++)
		{
			cCurrCtrl = controlVariates[i];
			sCurrCtrl = CCSTringToSTLString(cCurrCtrl);
			(*ctrlVariateVec)[i] = atoi(sCurrCtrl.c_str());
			(*ctrlVariatePrices)[i] = 0;
		}

		//MktDataManager: persistent object
		VECTOR<CCString> mktDataManager;
		XL_readStrVector(XL_mktDataManager,mktDataManager," ARM_ERR: Market datas: array of string expected",DOUBLE_TYPE,C_result);
        if(mktDataManager.size()<4)
		{
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}

		long mktDataManagerId = LocalGetNumObjectId (mktDataManager[0]);
		vector< string > mdmKeys(mktDataManager.size()-1);
		for(i=1;i<mktDataManager.size();++i)
			mdmKeys[i-1]=CCSTringToSTLString(mktDataManager[i]);

		//Model Parameters
		VECTOR<CCString> calibParamsXl;
		XL_readStrVector(XL_calibParams, calibParamsXl," ARM_ERR: Calib Parameters: array expected",DOUBLE_TYPE,C_result);
		if((calibParamsXl.size()<25) || (calibParamsXl.size()>29))

		{
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}
		
		//Model Type
		CCString modelTypeXl = calibParamsXl[0];
		string smodelTypeXl = CCSTringToSTLString(modelTypeXl);
		string smodelType(stringGetUpper(smodelTypeXl));
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
	
		vector< string > calibFlags(2);
		calibFlags[0] = calibMrs;
		calibFlags[1] = calibBeta;	

		//Num Method Type
		CCString methodTypeXl = calibParamsXl[10];
		string smethodTypeXl = CCSTringToSTLString(methodTypeXl);
		string smethodType(stringGetUpper(smethodTypeXl));
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
		string sptflMode(stringGetUpper(sptflModeXl));
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
			throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"Portfolio Mode should be SUMMIT, NEWMODE, MANUAL or FASTER (or BASKET_ATM or BASKET_EQUIV)");
		
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
		string sptflStrike(stringGetUpper(sptflStrikeXl));
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
		string sboundaryFlag(stringGetUpper(sboundaryFlagXl));
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
		string sapproxMarginFlag(stringGetUpper(sapproxMarginFlagXl));
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
		string sfreezeBetasFlag(stringGetUpper(sfreezeBetasFlagXl));
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
		CCString genType1Xl = calibParamsXl[21];
		string 	genType1 = CCSTringToSTLString(genType1Xl);
		genType1 = stringGetUpper(genType1);

		//STM Calib Mode
		CCString mrsptflModeXl = calibParamsXl[22];
		string smrsptflModeXl = CCSTringToSTLString(mrsptflModeXl);
		string smrsptflMode(stringGetUpper(smrsptflModeXl));
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
		string scalculateProbaFlag(stringGetUpper(scalculateProbaFlagXl));
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
		
		// Lag for volatility calibration
		CCString volLagXl = calibParamsXl[24];
		string svolLag = CCSTringToSTLString(volLagXl );
		int volLag = atoi(svolLag .c_str());
		portfolioMode[4] = volLag;


		// Gen Type 2
		CCString genType2XL = "Sobol";
		if (calibParamsXl.size() >= 26)
		{
			genType2XL = calibParamsXl[25];
		}
		string genType2 = CCSTringToSTLString(genType2XL );
		genType2 = stringGetUpper(genType2);


		// Path Order
		CCString pathOrderXL = "BucketOrder";
		if (calibParamsXl.size() >= 27)
		{
			pathOrderXL = calibParamsXl[26];
		}
		string pathOrder = CCSTringToSTLString(pathOrderXL );
		pathOrder = stringGetUpper(pathOrder);

		// Path Scheme
		CCString pathSchemeXL = "Incremental";
		if (calibParamsXl.size() >= 28)
		{
			pathSchemeXL = calibParamsXl[27];
		}
		string pathScheme = CCSTringToSTLString(pathSchemeXL );
		pathScheme = stringGetUpper(pathScheme);

		// firstNbDims
		CCString firstNbTimesXL = "0";
		if (calibParamsXl.size() >= 29)
		{
			firstNbTimesXL = calibParamsXl[28];
		}
		string firstNbTimesStr = CCSTringToSTLString(firstNbTimesXL );
		int firstNbTimes = atoi(firstNbTimesXL .c_str());

		//Functor call
			bermudaswaptionCalculatorFunc ourFunc(
				ccy,
				startDate,
				endDate,
				notional,
				notionalId,
				strike,
				strikeId,
				fees,
				feesId,
				spread,
				spreadId,
				payReceive,
				callFreq,
				callNotice,
				callCal,
				firstCallDate,
				lastCallDate,
				fixFreq,
				fixBasis,
				fixPayCal,
				fixAdjRule,
				fixRule,
				fixPayGap,
				isZc,
				varFreq,
				varBasis,
				varResetCal,
				varPayCal,
				varIndexTerm,
				varAdjRule,
				varRule,
				varResetGap,
				varPayGap,
				stubRule,
				genSecType,
				ctrlVariateVec,
				ctrlVariatePrices,
				mdmKeys,
				mktDataManagerId,
				modelType,
				modelParams,
				calibFlags,
				numMethodType,
				amcIter,
				mcIter,
				maxBucketSize,
				genType1,
				genType2,
				pathOrder,
				pathScheme,
				firstNbTimes,
				treeSteps,
				portfolioMode,
				boundaryFlag,
				approxMarginFlag,
				freezeBetasFlag,
				calculateProbaFlag);

		/// call the general function
		fillXL_Result( LOCAL_GC_BERMUDASWAPTION_CLASS, ourFunc, C_result, XL_result, PersistentInXL );
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_BermudaSwaptionCalculator_Common" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


/////////////////////////////////////////////////////////////
/// central function that does the creation of the XL function
/////////////////////////////////////////////////////////////
LPXLOPER Local_BermudaSwaptionCalc_Common(  LPXLOPER XL_currency,
											LPXLOPER XL_startDate,
											LPXLOPER XL_endDate,
											LPXLOPER XL_notionalCurve,
											LPXLOPER XL_strikeCurve,
											LPXLOPER XL_feesCurve,
											LPXLOPER XL_spreadCurve,
											LPXLOPER XL_payReceive,
											LPXLOPER XL_callDatas,
											LPXLOPER XL_fixDatas,
											LPXLOPER XL_varDatas,
											LPXLOPER XL_genSecType,
											LPXLOPER XL_mktDataManager,
											LPXLOPER XL_calibParams,
											LPXLOPER XL_controlVariates,
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

		//Currency
		CCString ccyStr;
		XL_readStrCell(XL_currency,ccyStr," ARM_ERR: currency: string expected",C_result);
		string sCcy = CCSTringToSTLString(ccyStr);
		ARM_Currency ccy(sCcy.c_str());

		/// Start & End
		double startDate;
		XL_readNumCell(XL_startDate,startDate," ARM_ERR: start date: date expected",C_result);
		double endDate;
		XL_readNumCell(XL_endDate,endDate," ARM_ERR: end date: date expected",C_result);

		//NotionalValue or NotionalCurve
		double notional;
		CCString notionalStr;
		long     notionalId;
		XL_readStrOrNumCell(XL_notionalCurve, notionalStr, notional, notionalId,
			   " ARM_ERR: notional: numeric or refValue Id expected",C_result);	
		if(notionalId == XL_TYPE_STRING)
			notionalId = LocalGetNumObjectId(notionalStr);
		else
			notionalId = ARM_NULL_OBJECT;

		/// StrikeValue or StrikeCurve
		double strike;
		CCString strikeStr;
		long     strikeId;
		XL_readStrOrNumCell(XL_strikeCurve, strikeStr, strike, strikeId,
			   " ARM_ERR: strike: numeric or refValue Id expected",C_result);	
		if(strikeId == XL_TYPE_STRING)
			strikeId = LocalGetNumObjectId(strikeStr);
		else
			strikeId = ARM_NULL_OBJECT;

		/// FeesValue or FeesCurve
		double fees;
		CCString feesStr;
		long     feesId;
		XL_readStrOrNumCell(XL_feesCurve, feesStr, fees, feesId,
			   " ARM_ERR: fees: numeric or refValue Id expected",C_result);	
		if(feesId == XL_TYPE_STRING)
			feesId = LocalGetNumObjectId(feesStr);
		else
			feesId = ARM_NULL_OBJECT;
	
		/// SpreadValue or SpreadCurve
		double spread;
		CCString spreadStr;
		long     spreadId;
		XL_readStrOrNumCell(XL_spreadCurve, spreadStr, spread, spreadId,
			   " ARM_ERR: spread: numeric or refValue Id expected",C_result);	
		if(spreadId == XL_TYPE_STRING)
			spreadId = LocalGetNumObjectId(spreadStr);
		else
			spreadId = ARM_NULL_OBJECT;

		/// Pay/Receive
		CCString payReceiveStr;
		long lPayReceive;
		XL_readStrCell(XL_payReceive,payReceiveStr," ARM_ERR: payer/receiver: string expected",C_result);
		if((lPayReceive = ARM_ConvRecOrPay (payReceiveStr, C_result)) == ARM_DEFAULT_ERR)
		{
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}
		int payReceive = (int)lPayReceive;

		/// Call Datas
		/// ----------
		VECTOR<CCString> callDatas;
		XL_readStrVector (XL_callDatas,callDatas," ARM_ERR: Call datas: array of string expected",DOUBLE_TYPE,C_result);
		if(callDatas.size() < 5)
		{
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}

		/// Call frequency 
		CCString callFreqXl=callDatas[0];
		long lCallFreq;
		if((lCallFreq = ARM_ConvFrequency (callFreqXl, C_result)) == ARM_DEFAULT_ERR)
		{
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}
		int callFreq = (int)lCallFreq;

		//CallNotice
		CCString callNoticeXl = callDatas[1];
		string sCallNotice = CCSTringToSTLString(callNoticeXl);
		int callNotice = atoi(sCallNotice.c_str());
		
		//CallCal
		CCString callCalXl = callDatas[2];
		string callCal = CCSTringToSTLString(callCalXl);
	
		//FirstCallDate
		CCString firstCallDateXl = callDatas[3];
		string sFirstCallDate = CCSTringToSTLString(firstCallDateXl);
		double firstCallDate = atoi(sFirstCallDate.c_str());
	
		//LastCallDate
		CCString lastCallDateXl = callDatas[4];
		string sLastCallDate = CCSTringToSTLString(lastCallDateXl);
		double lastCallDate = atoi(sLastCallDate.c_str());

		/// Underlying Fix Leg Datas
		/// ------------------------
		VECTOR<CCString> fixDatas;
		XL_readStrVector (XL_fixDatas, fixDatas," ARM_ERR: Fix datas: array of string expected",DOUBLE_TYPE,C_result);
	
		int fixFreq;
		int fixBasis;
		string fixPayCal;
		int fixAdjRule;
		int fixRule;
		int fixPayGap;
		bool isZc;
		
		//FixFreq + IsZC + default parameters.
		if(fixDatas.size() == 3) 
		{
			//Fix Freq
			CCString fixFreqXl = fixDatas[0];
			long lFixFreq;
			if((lFixFreq = ARM_ConvFrequency (fixFreqXl, C_result)) == ARM_DEFAULT_ERR)
			{
				ARM_ARG_ERR();
				return (LPXLOPER)&XL_result;
			}
			fixFreq = (int)lFixFreq;

			//Default: Fix Basis
			CCString fixBasisXl = fixDatas[1];
			fixBasis = (int)ARM_ConvDayCount (fixBasisXl);

			//Default: Fix Pay Cal
			fixPayCal = "ZGU";

			//Default: Fix Adj Rule
			fixAdjRule = 2; //Modified Following.
			
			//Default: Fix Rule
			fixRule = K_UNADJUSTED;

			//Default: Fix Pay Gap
			fixPayGap = 0;

			//Is ZC
			CCString isZcXl = fixDatas[2];
			string sisZcXl = CCSTringToSTLString(isZcXl);
			string sisZc(stringGetUpper(sisZcXl));
			if(sisZc == "Y")
			{
				isZc = true;
			}
			else if (sisZc == "N")
			{
				isZc = false;
			}
		}
		//All fix datas are inputed.
		else if (fixDatas.size() == 7)
		{
			//Fix Freq
			CCString fixFreqXl = fixDatas[0];
			long lFixFreq;
			if((lFixFreq = ARM_ConvFrequency (fixFreqXl, C_result)) == ARM_DEFAULT_ERR)
			{
				ARM_ARG_ERR();
				return (LPXLOPER)&XL_result;
			}
			fixFreq = (int)lFixFreq;

			//Fix Basis
			CCString fixBasisXl = fixDatas[1];
			fixBasis = (int)ARM_ConvDayCount (fixBasisXl);

			//Fix Pay Cal
			CCString fixPayCalXl = fixDatas[2];
			fixPayCal = CCSTringToSTLString(fixPayCalXl);

			//Fix Adj Rule
			CCString fixAdjRuleXl = fixDatas[3];
			fixAdjRule = (int)ARM_ConvFwdRule(fixAdjRuleXl);
			
			//Fix Rule
			CCString fixRuleXl = fixDatas[4];
			fixRule = (int)ARM_ConvIntRule(fixRuleXl);

			//Fix Pay Gap
			CCString fixPayGapXl = fixDatas[5];
			string sFixPayGap = CCSTringToSTLString(fixPayGapXl);
			fixPayGap = atoi(sFixPayGap.c_str());

			//Is ZC
			CCString isZcXl = fixDatas[6];
			string sisZcXl = CCSTringToSTLString(isZcXl);
			string sisZc(stringGetUpper(sisZcXl));
			if(sisZc == "Y")
			{
				isZc = true;
			}
			else if (sisZc == "N")
			{
				isZc = false;
			}
		}
		else
		{
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}

		/// Underlying Var Leg Datas
		/// ------------------------
		VECTOR<CCString> varDatas;
		XL_readStrVector (XL_varDatas, varDatas," ARM_ERR: Var datas: array of string expected",DOUBLE_TYPE,C_result);
		int varFreq;
		int varBasis;
		string varResetCal;
		string varPayCal;
		string varIndexTerm;
		int varAdjRule;
		int varRule;
		int varResetGap;
		int varPayGap;
		int stubRule;

		if(varDatas.size() == 1)
		{
			//Var Freq
			CCString varFreqXl = varDatas[0];
			long lVarFreq;
			if((lVarFreq = ARM_ConvFrequency (varFreqXl, C_result)) == ARM_DEFAULT_ERR)
			{
				ARM_ARG_ERR();
				return (LPXLOPER)&XL_result;
			}
			varFreq = (int)lVarFreq;
			
			//Var Basis
			varBasis = KACTUAL_360;

			//Var ResetCal
			varResetCal = "GBP";

			//Var Pay Cal
			varPayCal = "ZGU";

			//Var Index Term
			if (1 == varFreq)
				varIndexTerm = "12M";
			else if (2 == varFreq)
				varIndexTerm = "6M";
			else if (4 == varFreq)
				varIndexTerm = "3M";
			else if (12 == varFreq)
				varIndexTerm = "1M";

			//Var Adj Rule
			varAdjRule = 2; //Modified Following

			//Var Rule
			varRule = K_UNADJUSTED;

			//Var Reset Gap
			varResetGap = -2;

			//Var Pay Gap
			varPayGap = 0;

			//Stub Rule
			stubRule = K_SHORTSTART;
		}
		else if (varDatas.size() == 10)
		{
			//Var Freq
			CCString varFreqXl = varDatas[0];
			long lVarFreq;
			if((lVarFreq = ARM_ConvFrequency (varFreqXl, C_result)) == ARM_DEFAULT_ERR)
			{
				ARM_ARG_ERR();
				return (LPXLOPER)&XL_result;
			}
			varFreq = (int)lVarFreq;
		
			//Var Basis
			CCString varBasisXl = varDatas[1];
			varBasis = (int)ARM_ConvDayCount (varBasisXl);	
		
			//Var Reset Cal
			CCString varResetCalXl = varDatas[2];
			varResetCal = CCSTringToSTLString(varResetCalXl);

			//Var Pay Cal
			CCString varPayCalXl = varDatas[3];
			varPayCal = CCSTringToSTLString(varPayCalXl);

			//Var Index Term
			CCString varIndexTermXl = varDatas[4];
			varIndexTerm = CCSTringToSTLString(varIndexTermXl);

			//Var Adj Rule
			CCString varAdjRuleXl = varDatas[5];
			varAdjRule = (int)ARM_ConvFwdRule(varAdjRuleXl);
		
			//Var  Rule
			CCString varRuleXl = varDatas[6];
			varRule = (int)ARM_ConvIntRule(varRuleXl);
		
			//Var Reset Gap
			CCString varResetGapXl = varDatas[7];
			string sVarResetGap = CCSTringToSTLString(varResetGapXl);
			varResetGap = atoi(sVarResetGap.c_str());

			//Var Pay Gap
			CCString varPayGapXl = varDatas[8];
			string sVarPayGap = CCSTringToSTLString(varPayGapXl);
			varPayGap = atoi(sVarPayGap.c_str());

			//Underlying Swap Stub Rule
			CCString stubRuleXl = varDatas[9];
			stubRule = (int)ARM_ConvStubRule(stubRuleXl);
		}
		else
		{
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}

		//Generic Security Type
		double dgenSecType = 0;
		double dDgenSecType = 0;
		int genSecType;
		 XL_readNumCellWD(XL_genSecType,dgenSecType, dDgenSecType," ARM_ERR: Gen Sec Type: integer expected",C_result);
		if (dgenSecType == 0)
			genSecType = (int)-1;
		else
			genSecType = (int)dgenSecType;

		//ControlVariates
		VECTOR<CCString> controlVariates;
		VECTOR<CCString> controlVariatesDef(0);
		XL_readStrVectorWD(XL_controlVariates, controlVariates,controlVariatesDef, " ARM_ERR: control variates: array of string expected",DOUBLE_TYPE,C_result);
		int ctrlSize = controlVariates.size();
		vector<int>* ctrlVariateVec = new vector<int>(ctrlSize);
		vector<double>* ctrlVariatePrices = new vector<double>(ctrlSize);
		CCString cCurrCtrl;
		string sCurrCtrl;
		size_t i;
		for (i=0; i<ctrlSize;i++)
		{
			cCurrCtrl = controlVariates[i];
			sCurrCtrl = CCSTringToSTLString(cCurrCtrl);
			(*ctrlVariateVec)[i] = atoi(sCurrCtrl.c_str());
			(*ctrlVariatePrices)[i] = 0;
		}

		//MktDataManager: persistent object
		VECTOR<CCString> mktDataManager;
		XL_readStrVector(XL_mktDataManager,mktDataManager," ARM_ERR: Market datas: array of string expected",DOUBLE_TYPE,C_result);
        if(mktDataManager.size()<4)
		{
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}

		long mktDataManagerId = LocalGetNumObjectId (mktDataManager[0]);
		vector< string > mdmKeys(mktDataManager.size()-1);
		for(i=1;i<mktDataManager.size();++i)
			mdmKeys[i-1]=CCSTringToSTLString(mktDataManager[i]);

		//Model Parameters
		VECTOR<CCString> calibParamsXl;
		XL_readStrVector(XL_calibParams, calibParamsXl," ARM_ERR: Calib Parameters: array expected",DOUBLE_TYPE,C_result);
		vector<string> modelParams(7);
		int modelType;						//0
		double initTheta;					//4
		double initSkew;					//5
		double initMrsL;					//6
		double initMrsU;					//7
		string calibMrs;					//8
		string calibBeta;					//9
		vector< string > calibFlags(2);
		int numMethodType;					//10
		int amcIter;						//11
		int mcIter;							//12 
		int maxBucketSize;					//13
		int treeSteps;						//14
		//Portfolio Mode
		vector<int> portfolioMode(5);
		int ptflMode;						//15
		int ptflFreq;						//16
		int ptflStrike;						//17 
		bool boundaryFlag;					//18 
		bool approxMarginFlag;				//19
		bool freezeBetasFlag;				//20 
		string genType1;					//21
		int mrsptflMode;					//22
		bool calculateProbaFlag;			//23
		int volLag;							//24
		string genType2 = "Sobol";			//25
		string pathOrder = "BucketOrder";	//26
		string pathScheme = "Incremental";	//27
		int firstNbDims = 0;				//28


		if(calibParamsXl.size() == 8)
		{
			//Model Type
			CCString modelTypeXl = calibParamsXl[0];
			string smodelTypeXl = CCSTringToSTLString(modelTypeXl);
			string smodelType(stringGetUpper(smodelTypeXl));
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
		
			//Init Vol
			CCString initVolXl = calibParamsXl[1];
			string sinitVol = CCSTringToSTLString(initVolXl);
			modelParams[0] = sinitVol;
			
			//Init Mrs
			CCString initMrsXl = calibParamsXl[2];
			string sinitMrs = CCSTringToSTLString(initMrsXl);
			modelParams[1] = sinitMrs;
		
			//Init Beta
			CCString initBetaXl = calibParamsXl[3];
			string sinitBeta = CCSTringToSTLString(initBetaXl);
			modelParams[2] = sinitBeta;

			//Init Theta
			initTheta = 0.0;
			modelParams[3] = string("DEFAULT");

			//Init Skew
			initSkew = 0.0;
			modelParams[4] = string("DEFAULT");

			//Init Mrs L
			initMrsL = -0.1;
			modelParams[5] = string("DEFAULT");
			
			//Init Mrs U
			initMrsU = 0.1;
			modelParams[6] = string("DEFAULT");
		
			//Calib MRS/Beta ? (Y/N
			CCString calibMrsXl  = calibParamsXl[4];
			calibMrs  = CCSTringToSTLString(calibMrsXl);
			calibBeta = "Y";
			calibFlags[0] = calibMrs;
			calibFlags[1] = calibBeta;	

			//Num Method Type
			CCString methodTypeXl = calibParamsXl[5];
			string smethodTypeXl = CCSTringToSTLString(methodTypeXl);
			string smethodType(stringGetUpper(smethodTypeXl));
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
			CCString amcIterXl = calibParamsXl[6];
			string samcIter = CCSTringToSTLString(amcIterXl);
			amcIter = atoi(samcIter.c_str());

			//MC Iter
			mcIter = amcIter;
	
			//Num Max Bucket Size
			maxBucketSize = amcIter;

			//Num Tree Steps
			CCString treeStepsXl = calibParamsXl[7];
			string streeSteps = CCSTringToSTLString(treeStepsXl);
			treeSteps = atoi(streeSteps.c_str());

			//Calib Mode
			ptflMode = 3; //NEWMODE

			portfolioMode[0] = ptflMode;
	
			//Calib Freq
			ptflFreq = 1;  //Annual
			portfolioMode[1] = ptflFreq;

			//Calib Strike
			ptflStrike = 1; //ATM
			portfolioMode[2] = ptflStrike;

			//Fix Boundary Flag
			boundaryFlag = true;  //Y
		
			//Approx Margin Flag
			approxMarginFlag = true; //Y
			
			//Freeze Betas Flag
			freezeBetasFlag = true;  //Y
		
			//Generator Type
			genType1 = "MRGK5";

			//STM Calib Mode
			mrsptflMode = 1;  //CORREL
			portfolioMode[3] = mrsptflMode;

			//Calculate Proba Flag
			//OB: if we are in mode Tree 
			//and genSecType = -1 ou 2 
			//--> force calculateProbaFlag to True
			if((numMethodType == 1) && ( genSecType == -1 || genSecType == 2 ) )	
				calculateProbaFlag = true;
			else
				calculateProbaFlag = false;  //N
		
			// Lag for volatility calibration
			volLag = 1;
			portfolioMode[4] = volLag;
		}
		else if (calibParamsXl.size() == 25)
		{
			//Model Type
			CCString modelTypeXl = calibParamsXl[0];
			string smodelTypeXl = CCSTringToSTLString(modelTypeXl);
			string smodelType(stringGetUpper(smodelTypeXl));
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
		
			//Init Vol
			CCString initVolXl = calibParamsXl[1];
			string sinitVol = CCSTringToSTLString(initVolXl);
			modelParams[0] = sinitVol;
			
			//Init Mrs
			CCString initMrsXl = calibParamsXl[2];
			string sinitMrs = CCSTringToSTLString(initMrsXl);
			modelParams[1] = sinitMrs;

			//Init Beta
			CCString initBetaXl = calibParamsXl[3];
			string sinitBeta = CCSTringToSTLString(initBetaXl);
			modelParams[2] = sinitBeta;

			//Init Theta
			CCString initThetaXl = calibParamsXl[4];
			string sinitTheta = CCSTringToSTLString(initThetaXl);
			modelParams[3] = sinitTheta;

			//Init Skew
			CCString initSkewXl = calibParamsXl[5];
			string sinitSkew = CCSTringToSTLString(initSkewXl);
			modelParams[4] = sinitSkew;

			//Init Mrs L
			CCString initMrsLXl = calibParamsXl[6];
			string sinitMrsL = CCSTringToSTLString(initMrsLXl);
			modelParams[5] = sinitMrsL;
			
			//Init Mrs U
			CCString initMrsUXl = calibParamsXl[7];
			string sinitMrsU = CCSTringToSTLString(initMrsUXl);
			modelParams[6] = sinitMrsU;
		
			//Calib MRS/Beta ? (Y/N
			CCString calibMrsXl  = calibParamsXl[8];
			CCString calibBetaXl = calibParamsXl[9];
			calibMrs  = CCSTringToSTLString(calibMrsXl);
			calibBeta = CCSTringToSTLString(calibBetaXl);
	
			calibFlags[0] = calibMrs;
			calibFlags[1] = calibBeta;	

			//Num Method Type
			CCString methodTypeXl = calibParamsXl[10];
			string smethodTypeXl = CCSTringToSTLString(methodTypeXl);
			string smethodType(stringGetUpper(smethodTypeXl));
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
			amcIter = atoi(samcIter.c_str());
	
			//Num Mc Iter
			CCString mcIterXl = calibParamsXl[12];
			string smcIter = CCSTringToSTLString(mcIterXl);
			mcIter = atoi(smcIter.c_str());

			//Num Max Bucket Size
			CCString maxBucketSizeXl = calibParamsXl[13];
			string smaxBucket = CCSTringToSTLString(maxBucketSizeXl );
			maxBucketSize = atoi(smaxBucket.c_str());

			//Num Tree Steps
			CCString treeStepsXl = calibParamsXl[14];
			string streeSteps = CCSTringToSTLString(treeStepsXl);
			treeSteps = atoi(streeSteps.c_str());

			//Calib Mode
			CCString ptflModeXl = calibParamsXl[15];
			string sptflModeXl = CCSTringToSTLString(ptflModeXl);
			string sptflMode(stringGetUpper(sptflModeXl));
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
				throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"Portfolio Mode should be SUMMIT, NEWMODE, MANUAL or FASTER (or BASKET_ATM or BASKET_EQUIV)");
		

			portfolioMode[0] = ptflMode;
	
			//Calib Freq
			CCString ptflFreqXl = calibParamsXl[16];
			long lPtflFreq;
			if((lPtflFreq = ARM_ConvFrequency (ptflFreqXl, C_result)) == ARM_DEFAULT_ERR)
			{
				ARM_ARG_ERR();
				return (LPXLOPER)&XL_result;
			}
			ptflFreq = (int)lPtflFreq;
			portfolioMode[1] = ptflFreq;

			//Calib Strike
			CCString ptflStrikeXl = calibParamsXl[17];
			string sptflStrikeXl = CCSTringToSTLString(ptflStrikeXl);
			string sptflStrike(stringGetUpper(sptflStrikeXl));
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
			string sboundaryFlag(stringGetUpper(sboundaryFlagXl));
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
			string sapproxMarginFlag(stringGetUpper(sapproxMarginFlagXl));
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
			string sfreezeBetasFlag(stringGetUpper(sfreezeBetasFlagXl));
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
			CCString genType1Xl = calibParamsXl[21];
			genType1 = CCSTringToSTLString(genType1Xl);
			genType1 = stringGetUpper(genType1);

			//STM Calib Mode
			CCString mrsptflModeXl = calibParamsXl[22];
			string smrsptflModeXl = CCSTringToSTLString(mrsptflModeXl);
			string smrsptflMode(stringGetUpper(smrsptflModeXl));
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
			string scalculateProbaFlag(stringGetUpper(scalculateProbaFlagXl));
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
		
			// Lag for volatility calibration
			CCString volLagXl = calibParamsXl[24];
			string svolLag = CCSTringToSTLString(volLagXl );
			volLag = atoi(svolLag .c_str());
			portfolioMode[4] = volLag;
		}
		else
		{
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}

		//Functor call
		bermudaswaptionCalculatorFunc ourFunc(
				ccy,
				startDate,
				endDate,
				notional,
				notionalId,
				strike,
				strikeId,
				fees,
				feesId,
				spread,
				spreadId,
				payReceive,
				callFreq,
				callNotice,
				callCal,
				firstCallDate,
				lastCallDate,
				fixFreq,
				fixBasis,
				fixPayCal,
				fixAdjRule,
				fixRule,
				fixPayGap,
				isZc,
				varFreq,
				varBasis,
				varResetCal,
				varPayCal,
				varIndexTerm,
				varAdjRule,
				varRule,
				varResetGap,
				varPayGap,
				stubRule,
				genSecType,
				ctrlVariateVec,
				ctrlVariatePrices,
				mdmKeys,
				mktDataManagerId,
				modelType,
				modelParams,
				calibFlags,
				numMethodType,
				amcIter,
				mcIter,
				maxBucketSize,
				genType1,
				genType2,
				pathOrder,
				pathScheme,
				firstNbDims,
				treeSteps,
				portfolioMode,
				boundaryFlag,
				approxMarginFlag,
				freezeBetasFlag,
				calculateProbaFlag);

		/// call the general function
		fillXL_Result( LOCAL_GC_BERMUDASWAPTION_CLASS, ourFunc, C_result, XL_result, PersistentInXL );
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_BermudaSwaptionCalc_Common" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


__declspec(dllexport) LPXLOPER WINAPI Local_BermudaSwaptionCalculator_Create(
			LPXLOPER XL_currency,
			LPXLOPER XL_startDate,
			LPXLOPER XL_endDate,
			LPXLOPER XL_notionalCurve,
			LPXLOPER XL_strikeCurve,
			LPXLOPER XL_feesCurve,
			LPXLOPER XL_spreadCurve,
			LPXLOPER XL_payReceive,
			LPXLOPER XL_callDatas,
			LPXLOPER XL_fixDatas,
			LPXLOPER XL_varDatas,
			LPXLOPER XL_genSecType,
			LPXLOPER XL_mktDataManager,
			LPXLOPER XL_calibParams,
			LPXLOPER XL_controlVariates)
{
	ADD_LOG("Local_BermudaSwaptionCalculator_Create");
	bool PersistentInXL = true;

	return Local_BermudaSwaptionCalculator_Common(
		XL_currency,
	    XL_startDate,
        XL_endDate,
        XL_notionalCurve,
        XL_strikeCurve,
		XL_feesCurve,
		XL_spreadCurve,
        XL_payReceive,
        XL_callDatas,
        XL_fixDatas,
        XL_varDatas,
		XL_genSecType,
		XL_mktDataManager,
		XL_calibParams,
		XL_controlVariates,
		PersistentInXL);
}


///////////////////////////////////
/// version for non persistentInXl
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_PXL_BermudaSwaptionCalculator_Create(
    		LPXLOPER XL_currency,
    		LPXLOPER XL_startDate,
			LPXLOPER XL_endDate,
			LPXLOPER XL_notionalCurve,
			LPXLOPER XL_strikeCurve,
			LPXLOPER XL_feesCurve,
			LPXLOPER XL_spreadCurve,
			LPXLOPER XL_payReceive,
			LPXLOPER XL_callDatas,
			LPXLOPER XL_fixDatas,
			LPXLOPER XL_varDatas,
			LPXLOPER XL_genSecType,
			LPXLOPER XL_mktDataManager,
			LPXLOPER XL_calibParams,
			LPXLOPER XL_controlVariates)
{
	ADD_LOG("Local_PXL_BermudaSwaptionCalculator_Create");
	bool PersistentInXL = false;
	return Local_BermudaSwaptionCalculator_Common(
		XL_currency,
		XL_startDate,
	    XL_endDate,
        XL_notionalCurve,
        XL_strikeCurve,
		XL_feesCurve,
		XL_spreadCurve,
        XL_payReceive,
        XL_callDatas,
        XL_fixDatas,
        XL_varDatas,
		XL_genSecType,
		XL_mktDataManager,
		XL_calibParams,
		XL_controlVariates,
		PersistentInXL);
}

__declspec(dllexport) LPXLOPER WINAPI Local_BermudaSwaptionCalc_Create(
			LPXLOPER XL_currency,
			LPXLOPER XL_startDate,
			LPXLOPER XL_endDate,
			LPXLOPER XL_notionalCurve,
			LPXLOPER XL_strikeCurve,
			LPXLOPER XL_feesCurve,
			LPXLOPER XL_spreadCurve,
			LPXLOPER XL_payReceive,
			LPXLOPER XL_callDatas,
			LPXLOPER XL_fixDatas,
			LPXLOPER XL_varDatas,
			LPXLOPER XL_genSecType,
			LPXLOPER XL_mktDataManager,
			LPXLOPER XL_calibParams,
			LPXLOPER XL_controlVariates)
{
	ADD_LOG("Local_BermudaSwaptionCalc_Create");
	bool PersistentInXL = true;

	return Local_BermudaSwaptionCalc_Common(
		XL_currency,
	    XL_startDate,
        XL_endDate,
        XL_notionalCurve,
        XL_strikeCurve,
		XL_feesCurve,
		XL_spreadCurve,
        XL_payReceive,
        XL_callDatas,
        XL_fixDatas,
        XL_varDatas,
		XL_genSecType,
		XL_mktDataManager,
		XL_calibParams,
		XL_controlVariates,
		PersistentInXL);
}


///////////////////////////////////
/// version for non persistentInXl
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_PXL_BermudaSwaptionCalc_Create(
    		LPXLOPER XL_currency,
    		LPXLOPER XL_startDate,
			LPXLOPER XL_endDate,
			LPXLOPER XL_notionalCurve,
			LPXLOPER XL_strikeCurve,
			LPXLOPER XL_feesCurve,
			LPXLOPER XL_spreadCurve,
			LPXLOPER XL_payReceive,
			LPXLOPER XL_callDatas,
			LPXLOPER XL_fixDatas,
			LPXLOPER XL_varDatas,
			LPXLOPER XL_genSecType,
			LPXLOPER XL_mktDataManager,
			LPXLOPER XL_calibParams,
			LPXLOPER XL_controlVariates)
{
	ADD_LOG("Local_PXL_BermudaSwaptionCalc_Create");
	bool PersistentInXL = false;
	return Local_BermudaSwaptionCalc_Common(
		XL_currency,
		XL_startDate,
	    XL_endDate,
        XL_notionalCurve,
        XL_strikeCurve,
		XL_feesCurve,
		XL_spreadCurve,
        XL_payReceive,
        XL_callDatas,
        XL_fixDatas,
        XL_varDatas,
		XL_genSecType,
		XL_mktDataManager,
		XL_calibParams,
		XL_controlVariates,
		PersistentInXL);
}



__declspec(dllexport) LPXLOPER WINAPI Local_BermudaSwaptionCalculator_GetPricingData(
	LPXLOPER XL_BermudaSwaptionCalculatorId,
    LPXLOPER XL_KeyId )
{
	ADD_LOG("Local_BermudaSwaptionCalculator_GetPricingData");
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
		XL_readStrCell( XL_BermudaSwaptionCalculatorId, C_PricerString, " ARM_ERR: Bermuda Swaption Calculator Id: Object expected",C_result);
		long C_BermudaSwaptionCalcId = LocalGetNumObjectId(C_PricerString);
		
		CCString C_KeyString;
		XL_readStrCell( XL_KeyId, C_KeyString, " ARM_ERR: Key Id: Object expected",C_result);
		
		/// call the function
		ARM_GramFctorArg argResult;
		
		long retCode = ARMLOCAL_GC_GetPricingData(
			C_BermudaSwaptionCalcId,
			CCSTringToSTLString( C_KeyString ),
			argResult,
			C_result );
		
		if( retCode == ARM_KO )
		{
			ARM_ERR();
		}
		else
		{
			ARM_GramFunctorToXLOPER(argResult, XL_result, C_result);
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

///------------------------------------------------------------
///------------------------------------------------------------
///       Local Callable Range Accrual Calculator
/// Inputs :
///------------------------------------------------------------
///------------------------------------------------------------
class cfcraCalculatorFunc : public ARMResultLong2LongFunc
{
public:
	cfcraCalculatorFunc(	ARM_Currency ccy,
							double startDate,
							double endDate,
							int payReceive,
							int callFreq,
							int callNotice,
							string callCal,
							int fundFreq,
							int fundDayCount,
							int cpnPayFreq,
							string cpnResetCal,
							string cpnPayCal,
							int boostedIndexType,
							string boostedVarTerm,
							int boostedResetGap,
							int boostedResetTiming,
							int boostedDayCount,
							int boostedAdjRule,
							int boostedRule,
							int cpnResetFreq,
							int cpnResetTiming,
							int cpnResetGap,
							int refIndexType,
							string  refTerm,
							int refDayCount,
							double refCoeff,
							double notional,
							long notionalId,
							long callFeesId,
							double fundSpread,
							long fundSpreadId,
							double cpnSpread,
							long cpnSpreadId,
							double boostedFix,
							long boostedFixId,
							double bDown,
							long bDownId,
							double bUp,
							long bUpId,
							bool localModel,
							int localModelType,
							vector<double> mrsBeta,
							vector<double> calibSecPFParams,
							int nbSteps,
							int flagToGenerateOSWATM,
							int localResetFreq,
							vector<string> mdmKeys,
							long mktDataManagerId,
							vector<string> productsToPrice,
							bool isStdCalib)
    :
	C_currency(ccy),
    C_startDate(startDate),
    C_endDate(endDate),
	C_payReceive(payReceive),
	C_callFreq(callFreq),
	C_callNotice(callNotice),
	C_callCal(callCal),
	C_fundFreq(fundFreq),
	C_fundDayCount(fundDayCount),
	C_cpnPayFreq(cpnPayFreq),
	C_cpnResetCal(cpnResetCal),
	C_cpnPayCal(cpnPayCal),
	C_boostedIndexType(boostedIndexType),
	C_boostedVarTerm(boostedVarTerm),
	C_boostedResetGap(boostedResetGap),
	C_boostedResetTiming(boostedResetTiming),
	C_boostedDayCount(boostedDayCount),
	C_boostedAdjRule(boostedAdjRule),
	C_boostedRule(boostedRule),
	C_cpnResetFreq(cpnResetFreq),
	C_cpnResetTiming(cpnResetTiming),
	C_cpnResetGap(cpnResetGap),
	C_refIndexType(refIndexType),
	C_refTerm(refTerm),
	C_refDayCount(refDayCount),
	C_refCoeff(refCoeff),
	C_notional(notional),
	C_notionalId(notionalId),
	C_callFeesId(callFeesId),
	C_fundSpread(fundSpread),
	C_fundSpreadId(fundSpreadId),
	C_cpnSpread(cpnSpread),
	C_cpnSpreadId(cpnSpreadId),
	C_boostedFix(boostedFix),
	C_boostedFixId(boostedFixId),
	C_bDown(bDown),
	C_bDownId(bDownId),
	C_bUp(bUp),
	C_bUpId(bUpId),
	C_localModel(localModel),
	C_localModelType(localModelType),
	C_mrsBeta(mrsBeta),
	C_calibSecPFParams(calibSecPFParams),
	C_nbSteps(nbSteps),
	C_flagToGenerateOSWATM(flagToGenerateOSWATM),
	C_localResetFreq(localResetFreq),
	C_mdmKeys(mdmKeys),
	C_mktDataManagerId(mktDataManagerId),
	C_productsToPrice(productsToPrice),
	C_isStdCalib(isStdCalib)
    {
	};
	
	long operator()( ARM_result& result, long objId )
	{
		return ARMLOCAL_CRALocalCalculator_Create(
			C_currency,
			C_startDate,
			C_endDate,
			C_payReceive,
			C_callFreq,
			C_callNotice,
			C_callCal,
			C_fundFreq,
			C_fundDayCount,
			C_cpnPayFreq,
			C_cpnResetCal,
			C_cpnPayCal,
			C_boostedIndexType,
			C_boostedVarTerm,
			C_boostedResetGap,
			C_boostedResetTiming,
			C_boostedDayCount,
			C_boostedAdjRule,
			C_boostedRule,
			C_cpnResetFreq,
			C_cpnResetTiming,
			C_cpnResetGap,
			C_refIndexType,
			C_refTerm,
			C_refDayCount,
			C_refCoeff,
			C_notional,
			C_notionalId,
			C_callFeesId,
			C_fundSpread,
			C_fundSpreadId,
			C_cpnSpread,
			C_cpnSpreadId,
			C_boostedFix,
			C_boostedFixId,
			C_bDown,
			C_bDownId,
			C_bUp,
			C_bUpId,
			C_localModel,
			C_localModelType,
			C_mrsBeta,
			C_calibSecPFParams,
			C_nbSteps,
			C_flagToGenerateOSWATM,
			C_localResetFreq,
			C_mdmKeys,
			C_mktDataManagerId,
			C_productsToPrice,
			C_isStdCalib,
			result,
			objId);
    }

private:
ARM_Currency		C_currency;
double				C_startDate;
double				C_endDate;
int					C_payReceive;
int					C_callFreq;
int					C_callNotice;
string				C_callCal;
int					C_fundFreq;
int					C_fundDayCount;
int					C_cpnPayFreq;
string				C_cpnResetCal;
string				C_cpnPayCal;
int					C_boostedIndexType;
string				C_boostedVarTerm;
int					C_boostedResetGap;
int					C_boostedResetTiming;
int					C_boostedDayCount;
int					C_boostedIndexDayCount;
int					C_boostedAdjRule;
int					C_boostedRule;
int					C_cpnResetFreq;
int					C_cpnResetTiming;
int					C_cpnResetGap;
int					C_refIndexType;
string				C_refTerm;
int					C_refDayCount;
double				C_refCoeff;
double				C_notional;
long				C_notionalId;
long				C_callFeesId;
double				C_fundSpread;
long				C_fundSpreadId;
double				C_cpnSpread;
long				C_cpnSpreadId;
double				C_boostedFix;
long				C_boostedFixId;
double				C_bDown;
long				C_bDownId;
double				C_bUp;
long				C_bUpId;
bool				C_localModel;
int					C_localModelType;
vector<double>		C_mrsBeta;
vector<double>		C_calibSecPFParams;
int					C_nbSteps;
int					C_flagToGenerateOSWATM;
int					C_localResetFreq;
vector<string>		C_mdmKeys;
long				C_mktDataManagerId;
vector<string>		C_productsToPrice;
bool				C_isStdCalib;
};

LPXLOPER Local_CRALocalCalculator_Common(	LPXLOPER XL_generalDatas,
											LPXLOPER XL_callDatas,
											LPXLOPER XL_fundDatas,
											LPXLOPER XL_cpnDatas,
											LPXLOPER XL_notionalCurve,
											LPXLOPER XL_callFeesCurve,
											LPXLOPER XL_fundSpreadCurve,
											LPXLOPER XL_cpnSpreadCurve,
											LPXLOPER XL_boostedFixCurve,
											LPXLOPER XL_barrierDownCurve,
											LPXLOPER XL_barrierUpCurve,
											LPXLOPER XL_localModelParams,
											LPXLOPER XL_mrsBeta,
											LPXLOPER XL_calibSecPFParams,
											LPXLOPER XL_nbSteps,
											LPXLOPER XL_flagToGenerateOSWATM,
											LPXLOPER XL_mktDataManager,
											LPXLOPER XL_productsToPrice,
											LPXLOPER XL_isSfrmStdCalib,
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
	
		//General Datas
		//-------------
		VECTOR<CCString> generalDatas;
		XL_readStrVector (XL_generalDatas, generalDatas, "ARM_ERR: General datas: array of string expected", DOUBLE_TYPE, C_result);
		if (generalDatas.size() < 4)
		{
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}

		//Currency
		CCString ccyStr = generalDatas[0];
		string sCcy = CCSTringToSTLString(ccyStr);
		ARM_Currency ccy(sCcy.c_str());

		//Start & End
		CCString startDateXl = generalDatas[1];
		string sStartDate = CCSTringToSTLString(startDateXl);
		double startDate = atoi(sStartDate.c_str());
		CCString endDateXl = generalDatas[2];
		string sEndDate = CCSTringToSTLString(endDateXl);
		double endDate = atoi(sEndDate.c_str());
	    
		//PayReceive
		CCString payReceiveStr = generalDatas[3];
		long lPayReceive;
		if ((lPayReceive = ARM_ConvRecOrPay (payReceiveStr, C_result)) == ARM_DEFAULT_ERR)
		{
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}
		int payReceive = (int)lPayReceive;

		//CallDatas
		//---------
		VECTOR<CCString> callDatas;
		XL_readStrVector (XL_callDatas, callDatas, "ARM_ERR: Call datas: array of string expected", DOUBLE_TYPE, C_result);
		if (callDatas.size() < 3)
		{
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}

		/// Call frequency 
		CCString callFreqXl = callDatas[0];
		long lCallFreq;
		if ((lCallFreq = ARM_ConvFrequency (callFreqXl, C_result)) == ARM_DEFAULT_ERR)
		{
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}
		int callFreq = (int)lCallFreq;

		//CallNotice
		CCString callNoticeXl = callDatas[1];
		string sCallNotice = CCSTringToSTLString(callNoticeXl);
		int callNotice = atoi(sCallNotice.c_str());
		
		//CallCal
		CCString callCalXl = callDatas[2];
		string callCal = CCSTringToSTLString(callCalXl);

		//FundDatas
		//----------
		VECTOR<CCString> fundDatas;
		XL_readStrVector (XL_fundDatas, fundDatas," ARM_ERR: Fund datas: array of string expected",DOUBLE_TYPE,C_result);
		if (fundDatas.size() < 2)
		{
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}

		//Fund Freq
		CCString fundFreqXl = fundDatas[0];
		long lFundFreq;
		if ((lFundFreq = ARM_ConvFrequency (fundFreqXl, C_result)) == ARM_DEFAULT_ERR)
		{
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}
		int fundFreq = (int)lFundFreq;

		//Fund Day Count
		CCString fundDayCountXl = fundDatas[1];
		int fundDayCount = (int)ARM_ConvDayCount (fundDayCountXl);

		//CpnDatas
		//------------
		VECTOR<CCString> cpnDatas;
		XL_readStrVector (XL_cpnDatas, cpnDatas," ARM_ERR: Cpn datas: array of string expected",DOUBLE_TYPE,C_result);
		if (cpnDatas.size() < 17)
		{
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}

		//Cpn Pay Freq
		CCString cpnPayFreqXl = cpnDatas[0];
		long lCpnPayFreq;
		if ((lCpnPayFreq = ARM_ConvFrequency (cpnPayFreqXl, C_result)) == ARM_DEFAULT_ERR)
		{
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}
		int cpnPayFreq = (int)lCpnPayFreq;

		//Cpn Pay Gap
		CCString cpnResetGapXl = cpnDatas[1];
		string sCpnResetGap = CCSTringToSTLString(cpnResetGapXl);
		int cpnResetGap = atoi(sCpnResetGap.c_str());

		//Cpn Reset Cal
		CCString cpnResetCalXl = cpnDatas[2];
		string 	cpnResetCal = CCSTringToSTLString(cpnResetCalXl);

		//Cpn Pay Cal
		CCString cpnPayCalXl = cpnDatas[3];
		string 	cpnPayCal = CCSTringToSTLString(cpnPayCalXl);
		
		//Boosted Index Type 3
		CCString boostedIndexTypeXl = cpnDatas[4];
		string sboostedIndexTypeXl = CCSTringToSTLString(boostedIndexTypeXl);
		string sboostedIndexType(stringGetUpper(sboostedIndexTypeXl));
		int boostedIndexType;
		if (sboostedIndexType == "FIXED")
		{
			boostedIndexType = K_FIXED;
		}
		else if (sboostedIndexType == "LIBOR")
		{
			boostedIndexType = K_LIBOR;
		}
		else if (sboostedIndexType == "CMS")
		{
			boostedIndexType = K_CMS;
		}
		else
		{
			throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"only FIXED, LIBOR or CMS are Valid ");
		}

		//Boosted Var Term
		CCString boostedVarTermXl = cpnDatas[5];
		string 	boostedVarTerm = CCSTringToSTLString(boostedVarTermXl);

		//Boosted Reset Gap
		CCString boostedResetGapXl = cpnDatas[6];
		string sBoostedResetGap = CCSTringToSTLString(boostedResetGapXl);
		int boostedResetGap = atoi(sBoostedResetGap.c_str());

		//Boosted Reset Timing
		CCString boostedResetTimingXl = cpnDatas[7];
		int boostedResetTiming = (int)ARM_ConvPayResetRule(boostedResetTimingXl);

		//Boosted Day Count
		CCString boostedDayCountXl = cpnDatas[8];
		int boostedDayCount = (int)ARM_ConvDayCount(boostedDayCountXl);

		//Boosted Adj Rule
		CCString boostedAdjRuleXl = cpnDatas[9];
		int boostedAdjRule = (int)ARM_ConvFwdRule(boostedAdjRuleXl);

		//Boosted Rule
		CCString boostedRuleXl = cpnDatas[10];
		int boostedRule = (int)ARM_ConvIntRule(boostedRuleXl);

		//CpnResetFreq
		CCString cpnResetFreqXl = cpnDatas[11];
		long lCpnResetFreq;
		if ((lCpnResetFreq = ARM_ConvFrequency (cpnResetFreqXl, C_result)) == ARM_DEFAULT_ERR)
		{
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}
		int cpnResetFreq = (int)lCpnResetFreq;

		//Cpn Reset Timing
		CCString cpnResetTimingXl = cpnDatas[12];
		int cpnResetTiming = (int)ARM_ConvPayResetRule(cpnResetTimingXl);

		//Ref Index Type : Libor, Cms 
		CCString refIndexTypeXl = cpnDatas[13];
		string srefIndexTypeXl = CCSTringToSTLString(refIndexTypeXl);
		string srefIndexType(stringGetUpper(srefIndexTypeXl));
		int refIndexType;
		if (srefIndexType == "LIBOR")
		{
			refIndexType = K_LIBOR;
		}
		else if (srefIndexType == "CMS")
		{
			refIndexType = K_CMS;
		}
		else
		{
			throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"only LIBOR or CMS are Valid ");
		}

		//Ref Term 
		CCString refTermXl = cpnDatas[14];
		string 	refTerm = CCSTringToSTLString(refTermXl);

		//Ref Day Count  
		CCString refDayCountXl = cpnDatas[15];
		int refDayCount = (int)ARM_ConvDayCount (refDayCountXl);

		//Ref Coeff  
		CCString refCoeffXl = cpnDatas[16];
		string srefCoeff = CCSTringToSTLString(refCoeffXl);
		double refCoeff = atof(srefCoeff.c_str());

		//NotionalValue or NotionalCurve
		double notional;
		CCString notionalStr;
		long notionalId;
		XL_readStrOrNumCell(XL_notionalCurve, notionalStr, notional, notionalId,
			   " ARM_ERR: notional: numeric or refValue Id expected",C_result);	

		if (notionalId == XL_TYPE_STRING)
			notionalId = LocalGetNumObjectId(notionalStr);
		else
			notionalId = ARM_NULL_OBJECT;

		//Call Fees Curve 
		CCString callFeesStr;
		long callFeesId;
		XL_readStrCell(XL_callFeesCurve, callFeesStr," ARM_ERR: Call Fees: String Id expected", C_result);
		
		callFeesId = LocalGetNumObjectId(callFeesStr);

		//FundSpreadValue or FundSpreadCurve
		double fundSpread;
		CCString fundSpreadStr;
		long fundSpreadId;
		XL_readStrOrNumCell(XL_fundSpreadCurve, fundSpreadStr, fundSpread, fundSpreadId,
			   " ARM_ERR: fundSpread: numeric or refValue Id expected",C_result);	
		
		if (fundSpreadId == XL_TYPE_STRING)
			fundSpreadId = LocalGetNumObjectId(fundSpreadStr);
		else
			fundSpreadId = ARM_NULL_OBJECT;

		//CpnSpreadValue or CpnSpreadCurve
		double cpnSpread;
		CCString cpnSpreadStr;
		long cpnSpreadId;
		XL_readStrOrNumCell(XL_cpnSpreadCurve, cpnSpreadStr, cpnSpread, cpnSpreadId,
			   " ARM_ERR: cpnSpread: numeric or refValue Id expected",C_result);	
		
		if (cpnSpreadId == XL_TYPE_STRING)
			cpnSpreadId = LocalGetNumObjectId(cpnSpreadStr);
		else
			cpnSpreadId = ARM_NULL_OBJECT;

		//BoostedFixValue or BoostedFixCurve
		double boostedFix;
		CCString boostedFixStr;
		long boostedFixId;
		XL_readStrOrNumCell(XL_boostedFixCurve, boostedFixStr, boostedFix, boostedFixId,
			   " ARM_ERR: boosted fix rate: numeric or refValue Id expected",C_result);	
		
		if (boostedFixId == XL_TYPE_STRING)
			boostedFixId = LocalGetNumObjectId(boostedFixStr);
		else
			boostedFixId = ARM_NULL_OBJECT;

		//BarrierDownValue or BarrierDownCurve
		double bDown;
		CCString bDownStr;
		long bDownId;
		XL_readStrOrNumCell(XL_barrierDownCurve, bDownStr, bDown, bDownId,
			   " ARM_ERR: BarrierDown: numeric or refValue Id expected",C_result);	
		
		if (bDownId == XL_TYPE_STRING)
			bDownId = LocalGetNumObjectId(bDownStr);
		else
			bDownId = ARM_NULL_OBJECT;

		//BarrierUpValue or BarrierUpCurve
		double bUp;
		CCString bUpStr;
		long bUpId;
		XL_readStrOrNumCell(XL_barrierUpCurve, bUpStr, bUp, bUpId,
			   " ARM_ERR: BarrierUp: numeric or refValue Id expected",C_result);	
		
		if (bUpId == XL_TYPE_STRING)
			bUpId = LocalGetNumObjectId(bUpStr);
		else
			bUpId = ARM_NULL_OBJECT;

		//Local Model Params
		VECTOR<CCString> localModelParams;

		XL_readStrVector(XL_localModelParams,localModelParams," ARM_ERR: Local Model Params: array of string expected",DOUBLE_TYPE,C_result);
        if (localModelParams.size() > 3)
		{
			throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"0 to 5 products to price");
		}

		string localModelStr = CCSTringToSTLString(localModelParams[0]);
		bool localModel;
		if (localModelStr == "Y")
			localModel = true;
		else
			localModel = false;

		string localModelTypeStr = CCSTringToSTLString(localModelParams[1]);
		int localModelType;
		if (localModelTypeStr == "DOWN")
			localModelType = 0; //LocalDown;
		else if (localModelTypeStr == "UP")
			localModelType = 1; //LocalUp;
		else
			localModelType = 2; //LocalDownUp;

		//Local Reset Freq
		int localResetFreq; 
		if (localModelParams.size() == 3)
		{
			CCString localResetFreqXl = localModelParams[2];
			long lLocalResetFreq;
			if ((lLocalResetFreq = ARM_ConvFrequency (localResetFreqXl, C_result)) == ARM_DEFAULT_ERR)
			{
				ARM_ARG_ERR();
				return (LPXLOPER)&XL_result;
			}
			localResetFreq = (int)lLocalResetFreq;
		}
		else
		{
			localResetFreq = cpnResetFreq;
		}

		//Mean reversion and Beta
		VECTOR<double> mrsBeta;
		XL_readNumVector(XL_mrsBeta,mrsBeta," ARM_ERR: mrsBeta: vector of numbers expected",C_result);
		if ((mrsBeta.size() < 4) && (mrsBeta.size() > 6))
		{
			throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"Mrs & Beta: 4 values and 2 optional flags (0/1)");
		}

		//Calib Sec PF Params
		VECTOR<double> calibSecPFParams;
        XL_readNumVector(XL_calibSecPFParams,calibSecPFParams," ARM_ERR: calibSecPFParams: array of numeric expected",C_result);

		//Nb Steps
		double nbSteps;
		XL_readNumCell(XL_nbSteps,nbSteps," ARM_ERR: nbSteps: number expected",C_result);

		//Generate ATM Swaption flag
		CCString flagToGenerateOSWATMStr;
		int flagToGenerateOSWATM;
		XL_readStrCellWD(XL_flagToGenerateOSWATM,flagToGenerateOSWATMStr,"N"," ARM_ERR: flagToGenerateOSWATM: string (Y/N) expected",C_result);

		if (flagToGenerateOSWATMStr == "N")
			flagToGenerateOSWATM = 0;
		else if (flagToGenerateOSWATMStr == "Y")
			flagToGenerateOSWATM = 1;
		else
		{
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}

		//MktDataManager: persistent object
		VECTOR<CCString> mktDataManager;
		XL_readStrVector(XL_mktDataManager,mktDataManager," ARM_ERR: Market datas: array of string expected",DOUBLE_TYPE,C_result);
        
		if (mktDataManager.size()<4)
		{
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}

		long mktDataManagerId = LocalGetNumObjectId (mktDataManager[0]);
		vector<string> mdmKeys(mktDataManager.size()-1);		
		for (int i=1; i<mktDataManager.size(); ++i)
			mdmKeys[i-1] = CCSTringToSTLString(mktDataManager[i]);
	
		//Pricing flags
		VECTOR<CCString> pricingFlags;

		XL_readStrVector(XL_productsToPrice,pricingFlags," ARM_ERR: Pricing flags: array of string expected",DOUBLE_TYPE,C_result);
        if (pricingFlags.size()>7)
		{
			throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"maximum: 7 products for local deal and 5 for classic deal.");
		}

		vector< string > productsToPrice(pricingFlags.size());
		for (i=0; i<pricingFlags.size(); i++)
			productsToPrice[i] = CCSTringToSTLString(pricingFlags[i]);

		//Is SFRM standard calibration ? 
		bool isStdCalib;
		CCString isStdCalibStr;
		XL_readStrCellWD(XL_isSfrmStdCalib, isStdCalibStr,"N"," ARM_ERR: Is Std Calib: string (Y/N) expected",C_result);
		if (isStdCalibStr == "N")
			isStdCalib = false;
		else if (isStdCalibStr == "Y")
			isStdCalib = true;
		else
		{
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}

		//Functor call
		cfcraCalculatorFunc ourFunc(ccy,
									startDate,
									endDate,
									payReceive,
									callFreq,
									callNotice,
									callCal,
									fundFreq,
									fundDayCount,
									cpnPayFreq,
									cpnResetCal,
									cpnPayCal,
									boostedIndexType,
									boostedVarTerm,
									boostedResetGap,
									boostedResetTiming,
									boostedDayCount,
									boostedAdjRule,
									boostedRule,
									cpnResetFreq,
									cpnResetTiming,
									cpnResetGap,
									refIndexType,
									refTerm,
									refDayCount,
									refCoeff,
									notional,
									notionalId,
									callFeesId,
									fundSpread,
									fundSpreadId,
									cpnSpread,
									cpnSpreadId,
									boostedFix,
									boostedFixId,
									bDown,
									bDownId,
									bUp,
									bUpId,
									localModel,
									localModelType,
									mrsBeta,
									calibSecPFParams,
									nbSteps,
									flagToGenerateOSWATM,
									localResetFreq,
									mdmKeys,
									mktDataManagerId,
									productsToPrice,
									isStdCalib);

		/// call the general function
		fillXL_Result( LOCAL_GC_CRA_CLASS, ourFunc, C_result, XL_result, PersistentInXL );
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_CRALocalCalculator_Common" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}

///----------------------------------------------
///----------------------------------------------
///         Callable Range Accrual Calculator
/// Inputs :
///----------------------------------------------
///----------------------------------------------
class craCalculatorFunc : public ARMResultLong2LongFunc
{
public:
	craCalculatorFunc(	ARM_Currency ccy,
						double startDate,
						double endDate,
						int payReceive,
						int callFreq,
						int callNotice,
						string callCal,
						int fundFreq,
						int fundDayCount,
						int cpnPayFreq,
						string cpnResetCal,
						string cpnPayCal,
						int boostedIndexType,
						string boostedVarTerm,
						int boostedResetGap,
						int boostedResetTiming,
						int boostedDayCount,
						int boostedAdjRule,
						int boostedRule,
						int cpnResetFreq,
						int cpnResetTiming,
						int cpnResetGap,
						int refIndexType,
						string  refTerm,
						int refDayCount,
						double refCoeff,
						double notional,
						long notionalId,
						long callFeesId,
						double fundSpread,
						long fundSpreadId,
						double cpnSpread,
						long cpnSpreadId,
						double boostedFix,
						long boostedFixId,
						double bDown,
						long bDownId,
						double bUp,
						long bUpId,
						vector<double> mrsBeta,
						vector<double> calibSecPFParams,
						int nbSteps,
						int flagToGenerateOSWATM,
						vector<string> mdmKeys,
						long mktDataManagerId,
				        vector<string> productsToPrice,
						bool isStdCalib)
    :
	C_currency(ccy),
    C_startDate(startDate),
    C_endDate(endDate),
	C_payReceive(payReceive),
	C_callFreq(callFreq),
	C_callNotice(callNotice),
	C_callCal(callCal),
	C_fundFreq(fundFreq),
	C_fundDayCount(fundDayCount),
	C_cpnPayFreq(cpnPayFreq),
	C_cpnResetCal(cpnResetCal),
	C_cpnPayCal(cpnPayCal),
	C_boostedIndexType(boostedIndexType),
	C_boostedVarTerm(boostedVarTerm),
	C_boostedResetGap(boostedResetGap),
	C_boostedResetTiming(boostedResetTiming),
	C_boostedDayCount(boostedDayCount),
	C_boostedAdjRule(boostedAdjRule),
	C_boostedRule(boostedRule),
	C_cpnResetFreq(cpnResetFreq),
	C_cpnResetTiming(cpnResetTiming),
	C_cpnResetGap(cpnResetGap),
	C_refIndexType(refIndexType),
	C_refTerm(refTerm),
	C_refDayCount(refDayCount),
	C_refCoeff(refCoeff),
	C_notional(notional),
	C_notionalId(notionalId),
	C_callFeesId(callFeesId),
	C_fundSpread(fundSpread),
	C_fundSpreadId(fundSpreadId),
	C_cpnSpread(cpnSpread),
	C_cpnSpreadId(cpnSpreadId),
	C_boostedFix(boostedFix),
	C_boostedFixId(boostedFixId),
	C_bDown(bDown),
	C_bDownId(bDownId),
	C_bUp(bUp),
	C_bUpId(bUpId),
	C_mrsBeta(mrsBeta),
	C_calibSecPFParams(calibSecPFParams),
	C_nbSteps(nbSteps),
	C_flagToGenerateOSWATM(flagToGenerateOSWATM),
	C_mdmKeys(mdmKeys),
	C_mktDataManagerId(mktDataManagerId),
	C_productsToPrice(productsToPrice),
	C_isStdCalib(isStdCalib),
	C_optionPfId(-1)
    {
	};
	
	craCalculatorFunc(	int optionPfId,
						vector<double> mrsBeta,
						vector<double> calibSecPFParams,
						int nbSteps,
						int flagToGenerateOSWATM,
						vector<string> mdmKeys,
						long mktDataManagerId,
				        vector<string> productsToPrice,
						bool localModel,
						int localModelType,
						int  localResetFreq,
						bool isStdCalib)
    :
	C_optionPfId(optionPfId),
	C_mrsBeta(mrsBeta),
	C_calibSecPFParams(calibSecPFParams),
	C_nbSteps(nbSteps),
	C_flagToGenerateOSWATM(flagToGenerateOSWATM),
	C_mdmKeys(mdmKeys),
	C_mktDataManagerId(mktDataManagerId),
	C_productsToPrice(productsToPrice),
	C_localModel(localModel),
	C_localModelType(localModelType),
	C_localResetFreq(localResetFreq),
	C_isStdCalib(isStdCalib)
    {
	};
	
	long operator()( ARM_result& result, long objId )
	{
		if (C_optionPfId != -1)
		{
			return ARMLOCAL_CRACalculator_Create(
				C_optionPfId,
				C_mrsBeta,
				C_calibSecPFParams,
				C_nbSteps,
				C_flagToGenerateOSWATM,
				C_mdmKeys,
				C_mktDataManagerId,
				C_productsToPrice,
				C_localModel,
				C_localModelType,
				C_localResetFreq,
				C_isStdCalib,
				result,
				objId);
		}
		else
		{
			return ARMLOCAL_CRACalculator_Create(
				C_currency,
				C_startDate,
				C_endDate,
				C_payReceive,
				C_callFreq,
				C_callNotice,
				C_callCal,
				C_fundFreq,
				C_fundDayCount,
				C_cpnPayFreq,
				C_cpnResetCal,
				C_cpnPayCal,
				C_boostedIndexType,
				C_boostedVarTerm,
				C_boostedResetGap,
				C_boostedResetTiming,
				C_boostedDayCount,
				C_boostedAdjRule,
				C_boostedRule,
				C_cpnResetFreq,
				C_cpnResetTiming,
				C_cpnResetGap,
				C_refIndexType,
				C_refTerm,
				C_refDayCount,
				C_refCoeff,
				C_notional,
				C_notionalId,
				C_callFeesId,
				C_fundSpread,
				C_fundSpreadId,
				C_cpnSpread,
				C_cpnSpreadId,
				C_boostedFix,
				C_boostedFixId,
				C_bDown,
				C_bDownId,
				C_bUp,
				C_bUpId,
				C_mrsBeta,
				C_calibSecPFParams,
				C_nbSteps,
				C_flagToGenerateOSWATM,
				C_mdmKeys,
				C_mktDataManagerId,
				C_productsToPrice,
				C_isStdCalib,
				result,
				objId);
		}
    }

private:
int					C_optionPfId;
ARM_Currency		C_currency;
double				C_startDate;
double				C_endDate;
int					C_payReceive;
int					C_callFreq;
int					C_callNotice;
string				C_callCal;
int					C_fundFreq;
int					C_fundDayCount;
int					C_cpnPayFreq;
string				C_cpnResetCal;
string				C_cpnPayCal;
int					C_boostedIndexType;
string				C_boostedVarTerm;
int					C_boostedResetGap;
int					C_boostedResetTiming;
int					C_boostedDayCount;
int					C_boostedIndexDayCount;
int					C_boostedAdjRule;
int					C_boostedRule;
int					C_cpnResetFreq;
int					C_cpnResetTiming;
int					C_cpnResetGap;
int					C_refIndexType;
string				C_refTerm;
int					C_refDayCount;
double				C_refCoeff;
double				C_notional;
long				C_notionalId;
long				C_callFeesId;
double				C_fundSpread;
long				C_fundSpreadId;
double				C_cpnSpread;
long				C_cpnSpreadId;
double				C_boostedFix;
long				C_boostedFixId;
double				C_bDown;
long				C_bDownId;
double				C_bUp;
long				C_bUpId;
vector<double>		C_mrsBeta;
vector<double>		C_calibSecPFParams;
int					C_nbSteps;
int					C_flagToGenerateOSWATM;
vector<string>		C_mdmKeys;
long				C_mktDataManagerId;
vector<string>		C_productsToPrice;
bool				C_localModel;
int					C_localModelType;
int					C_localResetFreq;
bool				C_isStdCalib;
};


LPXLOPER Local_CRACalculator_Common(LPXLOPER XL_generalDatas,
									LPXLOPER XL_callDatas,
									LPXLOPER XL_fundDatas,
									LPXLOPER XL_cpnDatas,
									LPXLOPER XL_notionalCurve,
									LPXLOPER XL_callFeesCurve,
									LPXLOPER XL_fundSpreadCurve,
									LPXLOPER XL_cpnSpreadCurve,
									LPXLOPER XL_boostedFixCurve,
									LPXLOPER XL_barrierDownCurve,
									LPXLOPER XL_barrierUpCurve,
									LPXLOPER XL_mrsBeta,
									LPXLOPER XL_calibSecPFParams,
									LPXLOPER XL_nbSteps,
									LPXLOPER XL_flagToGenerateOSWATM,
									LPXLOPER XL_mktDataManager,
									LPXLOPER XL_productsToPrice,
									LPXLOPER XL_isSfrmStdCalib,
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
	
		//General Datas
		//-------------
		VECTOR<CCString> generalDatas;
		XL_readStrVector (XL_generalDatas, generalDatas, "ARM_ERR: General datas: array of string expected", DOUBLE_TYPE, C_result);
		if (generalDatas.size() < 4)
		{
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}

		//Currency
		CCString ccyStr = generalDatas[0];
		string sCcy = CCSTringToSTLString(ccyStr);
		ARM_Currency ccy(sCcy.c_str());

		//Start & End
		CCString startDateXl = generalDatas[1];
		string sStartDate = CCSTringToSTLString(startDateXl);
		double startDate = atoi(sStartDate.c_str());
		CCString endDateXl = generalDatas[2];
		string sEndDate = CCSTringToSTLString(endDateXl);
		double endDate = atoi(sEndDate.c_str());
	    
		//PayReceive
		CCString payReceiveStr = generalDatas[3];
		long lPayReceive;
		if ((lPayReceive = ARM_ConvRecOrPay (payReceiveStr, C_result)) == ARM_DEFAULT_ERR)
		{
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}
		int payReceive = (int)lPayReceive;

		//CallDatas
		//---------
		VECTOR<CCString> callDatas;
		XL_readStrVector (XL_callDatas, callDatas, "ARM_ERR: Call datas: array of string expected", DOUBLE_TYPE, C_result);
		if (callDatas.size() < 3)
		{
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}

		/// Call frequency 
		CCString callFreqXl = callDatas[0];
		long lCallFreq;
		if ((lCallFreq = ARM_ConvFrequency (callFreqXl, C_result)) == ARM_DEFAULT_ERR)
		{
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}
		int callFreq = (int)lCallFreq;

		//CallNotice
		CCString callNoticeXl = callDatas[1];
		string sCallNotice = CCSTringToSTLString(callNoticeXl);
		int callNotice = atoi(sCallNotice.c_str());
		
		//CallCal
		CCString callCalXl = callDatas[2];
		string callCal = CCSTringToSTLString(callCalXl);

		//FundDatas
		//----------
		VECTOR<CCString> fundDatas;
		XL_readStrVector (XL_fundDatas, fundDatas," ARM_ERR: Fund datas: array of string expected",DOUBLE_TYPE,C_result);
		if (fundDatas.size() < 2)
		{
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}

		//Fund Freq
		CCString fundFreqXl = fundDatas[0];
		long lFundFreq;
		if ((lFundFreq = ARM_ConvFrequency (fundFreqXl, C_result)) == ARM_DEFAULT_ERR)
		{
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}
		int fundFreq = (int)lFundFreq;

		//Fund Day Count
		CCString fundDayCountXl = fundDatas[1];
		int fundDayCount = (int)ARM_ConvDayCount (fundDayCountXl);

		//CpnDatas
		//------------
		VECTOR<CCString> cpnDatas;
		XL_readStrVector (XL_cpnDatas, cpnDatas," ARM_ERR: Cpn datas: array of string expected",DOUBLE_TYPE,C_result);
		if (cpnDatas.size() < 17)
		{
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}

		//Cpn Pay Freq
		CCString cpnPayFreqXl = cpnDatas[0];
		long lCpnPayFreq;
		if ((lCpnPayFreq = ARM_ConvFrequency (cpnPayFreqXl, C_result)) == ARM_DEFAULT_ERR)
		{
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}
		int cpnPayFreq = (int)lCpnPayFreq;

		//Cpn Pay Gap
		CCString cpnResetGapXl = cpnDatas[1];
		string sCpnResetGap = CCSTringToSTLString(cpnResetGapXl);
		int cpnResetGap = atoi(sCpnResetGap.c_str());

		//Cpn Reset Cal
		CCString cpnResetCalXl = cpnDatas[2];
		string 	cpnResetCal = CCSTringToSTLString(cpnResetCalXl);

		//Cpn Pay Cal
		CCString cpnPayCalXl = cpnDatas[3];
		string 	cpnPayCal = CCSTringToSTLString(cpnPayCalXl);
		
		//Boosted Index Type 3
		CCString boostedIndexTypeXl = cpnDatas[4];
		string sboostedIndexTypeXl = CCSTringToSTLString(boostedIndexTypeXl);
		string sboostedIndexType(stringGetUpper(sboostedIndexTypeXl));
		int boostedIndexType;
		if (sboostedIndexType == "FIXED")
		{
			boostedIndexType = K_FIXED;
		}
		else if (sboostedIndexType == "LIBOR")
		{
			boostedIndexType = K_LIBOR;
		}
		else if (sboostedIndexType == "CMS")
		{
			boostedIndexType = K_CMS;
		}
		else
		{
			throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"only FIXED, LIBOR or CMS are Valid ");
		}

		//Boosted Var Term
		CCString boostedVarTermXl = cpnDatas[5];
		string 	boostedVarTerm = CCSTringToSTLString(boostedVarTermXl);

		//Boosted Reset Gap
		CCString boostedResetGapXl = cpnDatas[6];
		string sBoostedResetGap = CCSTringToSTLString(boostedResetGapXl);
		int boostedResetGap = atoi(sBoostedResetGap.c_str());

		//Boosted Reset Timing
		CCString boostedResetTimingXl = cpnDatas[7];
		int boostedResetTiming = (int)ARM_ConvPayResetRule(boostedResetTimingXl);

		//Boosted Day Count
		CCString boostedDayCountXl = cpnDatas[8];
		int boostedDayCount = (int)ARM_ConvDayCount(boostedDayCountXl);

		//Boosted Adj Rule
		CCString boostedAdjRuleXl = cpnDatas[9];
		int boostedAdjRule = (int)ARM_ConvFwdRule(boostedAdjRuleXl);

		//Boosted Rule
		CCString boostedRuleXl = cpnDatas[10];
		int boostedRule = (int)ARM_ConvIntRule(boostedRuleXl);

		//CpnResetFreq
		CCString cpnResetFreqXl = cpnDatas[11];
		long lCpnResetFreq;
		if ((lCpnResetFreq = ARM_ConvFrequency (cpnResetFreqXl, C_result)) == ARM_DEFAULT_ERR)
		{
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}
		int cpnResetFreq = (int)lCpnResetFreq;

		//Cpn Reset Timing
		CCString cpnResetTimingXl = cpnDatas[12];
		int cpnResetTiming = (int)ARM_ConvPayResetRule(cpnResetTimingXl);

		//Ref Index Type: Libor, Cms 
		CCString refIndexTypeXl = cpnDatas[13];
		string srefIndexTypeXl = CCSTringToSTLString(refIndexTypeXl);
		string srefIndexType(stringGetUpper(srefIndexTypeXl));
		int refIndexType;
		if (srefIndexType == "LIBOR")
		{
			refIndexType = K_LIBOR;
		}
		else if (srefIndexType == "CMS")
		{
			refIndexType = K_CMS;
		}
		else
		{
			throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"only LIBOR or CMS are Valid ");
		}

		//Ref Term 
		CCString refTermXl = cpnDatas[14];
		string 	refTerm = CCSTringToSTLString(refTermXl);

		//Ref Day Count 
		CCString refDayCountXl = cpnDatas[15];
		int refDayCount = (int)ARM_ConvDayCount (refDayCountXl);

		//Ref Coeff		
		CCString refCoeffXl = cpnDatas[16];
		string srefCoeff = CCSTringToSTLString(refCoeffXl);
		double refCoeff = atof(srefCoeff.c_str());

		//NotionalValue or NotionalCurve
		double notional;
		CCString notionalStr;
		long notionalId;
		XL_readStrOrNumCell(XL_notionalCurve, notionalStr, notional, notionalId,
			   " ARM_ERR: notional: numeric or refValue Id expected",C_result);	

		if (notionalId == XL_TYPE_STRING)
			notionalId = LocalGetNumObjectId(notionalStr);
		else
			notionalId = ARM_NULL_OBJECT;

		//Call Fees Curve //To change
		CCString callFeesStr;
		long callFeesId;
		XL_readStrCell(XL_callFeesCurve, callFeesStr," ARM_ERR: Call Fees: String Id expected", C_result);
		
		callFeesId = LocalGetNumObjectId(callFeesStr);

		//FundSpreadValue or FundSpreadCurve
		double fundSpread;
		CCString fundSpreadStr;
		long fundSpreadId;
		XL_readStrOrNumCell(XL_fundSpreadCurve, fundSpreadStr, fundSpread, fundSpreadId,
			   " ARM_ERR: fundSpread: numeric or refValue Id expected",C_result);	
		
		if (fundSpreadId == XL_TYPE_STRING)
			fundSpreadId = LocalGetNumObjectId(fundSpreadStr);
		else
			fundSpreadId = ARM_NULL_OBJECT;

		//CpnSpreadValue or CpnSpreadCurve
		double cpnSpread;
		CCString cpnSpreadStr;
		long cpnSpreadId;
		XL_readStrOrNumCell(XL_cpnSpreadCurve, cpnSpreadStr, cpnSpread, cpnSpreadId,
			   " ARM_ERR: cpnSpread: numeric or refValue Id expected",C_result);	
		
		if (cpnSpreadId == XL_TYPE_STRING)
			cpnSpreadId = LocalGetNumObjectId(cpnSpreadStr);
		else
			cpnSpreadId = ARM_NULL_OBJECT;

		//BoostedFixValue or BoostedFixCurve
		double boostedFix;
		CCString boostedFixStr;
		long boostedFixId;
		XL_readStrOrNumCell(XL_boostedFixCurve, boostedFixStr, boostedFix, boostedFixId,
			   " ARM_ERR: boosted fix rate: numeric or refValue Id expected",C_result);	
		
		if (boostedFixId == XL_TYPE_STRING)
			boostedFixId = LocalGetNumObjectId(boostedFixStr);
		else
			boostedFixId = ARM_NULL_OBJECT;

		//BarrierDownValue or BarrierDownCurve
		double bDown;
		CCString bDownStr;
		long bDownId;
		XL_readStrOrNumCell(XL_barrierDownCurve, bDownStr, bDown, bDownId,
			   " ARM_ERR: BarrierDown: numeric or refValue Id expected",C_result);	
		
		if (bDownId == XL_TYPE_STRING)
			bDownId = LocalGetNumObjectId(bDownStr);
		else
			bDownId = ARM_NULL_OBJECT;

		//BarrierUpValue or BarrierUpCurve
		double bUp;
		CCString bUpStr;
		long bUpId;
		XL_readStrOrNumCell(XL_barrierUpCurve, bUpStr, bUp, bUpId,
			   " ARM_ERR: BarrierUp: numeric or refValue Id expected",C_result);	
		
		if (bUpId == XL_TYPE_STRING)
			bUpId = LocalGetNumObjectId(bUpStr);
		else
			bUpId = ARM_NULL_OBJECT;

		//Mean reversion and Beta
		VECTOR<double> mrsBeta;
		XL_readNumVector(XL_mrsBeta,mrsBeta," ARM_ERR: mrsBeta: vector of numbers expected",C_result);
		if ((mrsBeta.size() < 4) && (mrsBeta.size() > 6))
		{
			throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"Mrs & Beta: 4 values and 2 optional flags (0/1)");
		}

		//Calib Sec PF Params
		VECTOR<double> calibSecPFParams;
        XL_readNumVector(XL_calibSecPFParams,calibSecPFParams," ARM_ERR: calibSecPFParams: array of numeric expected",C_result);

		//Nb Steps
		double nbSteps;
		XL_readNumCell(XL_nbSteps,nbSteps," ARM_ERR: nbSteps: number expected",C_result);

		//Generate ATM Swaption flag
		CCString flagToGenerateOSWATMStr;
		int flagToGenerateOSWATM;
		XL_readStrCellWD(XL_flagToGenerateOSWATM,flagToGenerateOSWATMStr,"N"," ARM_ERR: flagToGenerateOSWATM: string (Y/N) expected",C_result);

		if (flagToGenerateOSWATMStr == "N")
			flagToGenerateOSWATM = 0;
		else if (flagToGenerateOSWATMStr == "Y")
			flagToGenerateOSWATM = 1;
		else
		{
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}

		//MktDataManager: persistent object
		VECTOR<CCString> mktDataManager;
		XL_readStrVector(XL_mktDataManager,mktDataManager," ARM_ERR: Market datas: array of string expected",DOUBLE_TYPE,C_result);
        
		if (mktDataManager.size()<4)
		{
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}

		long mktDataManagerId = LocalGetNumObjectId (mktDataManager[0]);
		vector<string> mdmKeys(mktDataManager.size()-1);		
		for (int i=1; i<mktDataManager.size(); ++i)
			mdmKeys[i-1] = CCSTringToSTLString(mktDataManager[i]);
	
		//Pricing flags
		VECTOR<CCString> pricingFlags;

		XL_readStrVector(XL_productsToPrice,pricingFlags," ARM_ERR: Pricing flags: array of string expected",DOUBLE_TYPE,C_result);
        if (pricingFlags.size()>7)
		{
			throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"maximum: 7 products to price");
		}

		vector< string > productsToPrice(pricingFlags.size());
		for (i=0; i<pricingFlags.size(); i++)
			productsToPrice[i] = CCSTringToSTLString(pricingFlags[i]);

		//Is SFRM standard calibration ? 
		bool isStdCalib;
		CCString isStdCalibStr;
		XL_readStrCellWD(XL_isSfrmStdCalib, isStdCalibStr,"Y"," ARM_ERR: Is Std Calib: string (Y/N) expected",C_result);
		if (isStdCalibStr == "N")
			isStdCalib = false;
		else if (isStdCalibStr == "Y")
			isStdCalib = true;
		else
		{
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}

		//Functor call
		craCalculatorFunc ourFunc(  ccy,
									startDate,
									endDate,
									payReceive,
									callFreq,
									callNotice,
									callCal,
									fundFreq,
									fundDayCount,
									cpnPayFreq,
									cpnResetCal,
									cpnPayCal,
									boostedIndexType,
									boostedVarTerm,
									boostedResetGap,
									boostedResetTiming,
									boostedDayCount,
									boostedAdjRule,
									boostedRule,
									cpnResetFreq,
									cpnResetTiming,
									cpnResetGap,
									refIndexType,
									refTerm,
									refDayCount,
									refCoeff,
									notional,
									notionalId,
									callFeesId,
									fundSpread,
									fundSpreadId,
									cpnSpread,
									cpnSpreadId,
									boostedFix,
									boostedFixId,
									bDown,
									bDownId,
									bUp,
									bUpId,
									mrsBeta,
									calibSecPFParams,
									nbSteps,
									flagToGenerateOSWATM,
									mdmKeys,
									mktDataManagerId,
									productsToPrice,
									isStdCalib);
	
		/// call the general function
		fillXL_Result( LOCAL_GC_CRA_CLASS, ourFunc, C_result, XL_result, PersistentInXL );	
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_CRACalculator_Common" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


LPXLOPER Local_CRACalculator_Common(LPXLOPER XL_optionPortfolio,
									LPXLOPER XL_mrsBeta,
									LPXLOPER XL_calibSecPFParams,
									LPXLOPER XL_nbSteps,
									LPXLOPER XL_flagToGenerateOSWATM,
									LPXLOPER XL_mktDataManager,
									LPXLOPER XL_productsToPrice,
									LPXLOPER XL_localModelParams,
									LPXLOPER XL_isSfrmStdCalib,
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
	
		CCString C_optionPortfolio;
		XL_readStrCell( XL_optionPortfolio, C_optionPortfolio, " ARM_ERR: Option Portfolio Id: Object expected",C_result);
		long C_optionPortfolioId = LocalGetNumObjectId(C_optionPortfolio);

		//Mean reversion and Beta
		VECTOR<double> mrsBeta;
		VECTOR<double> mrsBetaDef(0);
		XL_readNumVectorWD(XL_mrsBeta,mrsBeta,mrsBetaDef," ARM_ERR: mrsBeta: vector of numbers expected",C_result);
		if ((mrsBeta.size() > 0) && (mrsBeta.size() < 4) && (mrsBeta.size() > 6))
		{
			throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"Mrs & Beta: 4 values and 2 optional flags (0/1)");
		}

		//Calib Sec PF Params
		VECTOR<double> calibSecPFParams;
		VECTOR<double> calibSecPFParamsDef(0);
        XL_readNumVectorWD(XL_calibSecPFParams,calibSecPFParams,calibSecPFParamsDef," ARM_ERR: calibSecPFParams: array of numeric expected",C_result);

		//Nb Steps
		double nbSteps;
		double nbStepsDef = -1.0;
		XL_readNumCellWD(XL_nbSteps,nbSteps,nbStepsDef," ARM_ERR: nbSteps: number expected",C_result);

		//Generate ATM Swaption flag
		CCString flagToGenerateOSWATMStr;
		int flagToGenerateOSWATM;
		XL_readStrCellWD(XL_flagToGenerateOSWATM,flagToGenerateOSWATMStr,"NULL"," ARM_ERR: flagToGenerateOSWATM: string (Y/N) expected",C_result);

		if (flagToGenerateOSWATMStr == "N")
			flagToGenerateOSWATM = 0;
		else if (flagToGenerateOSWATMStr == "Y")
			flagToGenerateOSWATM = 1;
		else
			flagToGenerateOSWATM = -1;

		//MktDataManager: persistent object
		VECTOR<CCString> mktDataManager;
		XL_readStrVector(XL_mktDataManager,mktDataManager," ARM_ERR: Market datas: array of string expected",DOUBLE_TYPE,C_result);
        
		if (mktDataManager.size()<4)
		{
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}

		long mktDataManagerId = LocalGetNumObjectId (mktDataManager[0]);
		vector<string> mdmKeys(mktDataManager.size()-1);		
		for (int i=1; i<mktDataManager.size(); ++i)
			mdmKeys[i-1] = CCSTringToSTLString(mktDataManager[i]);
	
		//Pricing flags
		VECTOR<CCString> pricingFlags;

		XL_readStrVector(XL_productsToPrice,pricingFlags," ARM_ERR: Pricing flags: array of string expected",DOUBLE_TYPE,C_result);
        if (pricingFlags.size()>7)
		{
			throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"maximum: 7 products for classic deal and 5 for local deal");
		}

		vector< string > productsToPrice(pricingFlags.size());
		for (i=0; i<pricingFlags.size(); i++)
			productsToPrice[i] = CCSTringToSTLString(pricingFlags[i]);

		//Local Model Params: 
		VECTOR<CCString> localModelParams;
		VECTOR<CCString> localModelParamsDef(0);
		XL_readStrVectorWD(XL_localModelParams, localModelParams, localModelParamsDef, " ARM_ERR: Local Model Params: array of string expected", DOUBLE_TYPE, C_result);
		
		if (localModelParams.size() > 3)
		{
			throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"0 to 5 products to price");
		}
		bool localModel = false;
		int localModelType;
		int localResetFreq; 
		
		if (localModelParams.size() != 0)
		{
			string localModelStr = CCSTringToSTLString(localModelParams[0]);
			if (localModelStr == "Y")
				localModel = true;
			else
				localModel = false;

			string localModelTypeStr = CCSTringToSTLString(localModelParams[1]);
			if (localModelTypeStr == "DOWN")
				localModelType = 0; //LocalDown;
			else if(localModelTypeStr == "UP")
				localModelType = 1; //LocalUp;
			else
				localModelType = 2; //LocalDownUp;

			//Local ResetFreq
			if (localModelParams.size() == 3)
			{
				CCString localResetFreqXl = localModelParams[2];
				long lLocalResetFreq;
				if ((lLocalResetFreq = ARM_ConvFrequency (localResetFreqXl, C_result)) == ARM_DEFAULT_ERR)
				{
					ARM_ARG_ERR();
					return (LPXLOPER)&XL_result;
				}
				localResetFreq = (int)lLocalResetFreq;
			}
		}
		else
		{
			localResetFreq = 1;
		}

		//Is SFRM standard calibration ? 
		bool isStdCalib;
		CCString isStdCalibStr;
		XL_readStrCellWD(XL_isSfrmStdCalib, isStdCalibStr,"N"," ARM_ERR: Is Std Calib: string (Y/N) expected",C_result);
		if (isStdCalibStr == "N")
			isStdCalib = false;
		else if (isStdCalibStr == "Y")
			isStdCalib = true;
		else
		{
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}

		//Functor call
		craCalculatorFunc  ourFunc( C_optionPortfolioId,
									mrsBeta,
									calibSecPFParams,
									nbSteps,
									flagToGenerateOSWATM,
									mdmKeys,
									mktDataManagerId,
									productsToPrice,
									localModel,
									localModelType,
									localResetFreq,
									isStdCalib);

		/// call the general function
		fillXL_Result( LOCAL_GC_CRA_CLASS, ourFunc, C_result, XL_result, PersistentInXL );
		
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_CRACalculator_Common" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}



///----------------------------------------------
///----------------------------------------------
///     CRA Local Calculator Accessor
/// Inputs :
///     
///----------------------------------------------
///----------------------------------------------
__declspec(dllexport) LPXLOPER WINAPI Local_CRALocalCalculator_Create(
			LPXLOPER XL_generalDatas,
			LPXLOPER XL_callDatas,
			LPXLOPER XL_fundDatas,
			LPXLOPER XL_cpnDatas,
			LPXLOPER XL_notionalCurve,
			LPXLOPER XL_callFeesCurve,
			LPXLOPER XL_fundSpreadCurve,
			LPXLOPER XL_cpnSpreadCurve,
			LPXLOPER XL_boostedFixCurve,
			LPXLOPER XL_barrierDownCurve,
			LPXLOPER XL_barrierUpCurve,
			LPXLOPER XL_localModelParams,
			LPXLOPER XL_mrsBeta,
			LPXLOPER XL_calibSecPFParams,
			LPXLOPER XL_nbSteps,
			LPXLOPER XL_flagToGenerateOSWATM,
			LPXLOPER XL_mktDataManager,
			LPXLOPER XL_productsToPrice,
			LPXLOPER XL_isSfrmStdCalib)
{
	ADD_LOG("Local_CRALocalCalculator_Create");
	bool PersistentInXL = true;

	return Local_CRALocalCalculator_Common(
			XL_generalDatas,
			XL_callDatas,
			XL_fundDatas,
			XL_cpnDatas,
			XL_notionalCurve,
			XL_callFeesCurve,
			XL_fundSpreadCurve,
			XL_cpnSpreadCurve,
			XL_boostedFixCurve,
			XL_barrierDownCurve,
			XL_barrierUpCurve,
			XL_localModelParams,
			XL_mrsBeta,
			XL_calibSecPFParams,
			XL_nbSteps,
			XL_flagToGenerateOSWATM,
			XL_mktDataManager,
			XL_productsToPrice,
			XL_isSfrmStdCalib,
			PersistentInXL);
}

__declspec(dllexport) LPXLOPER WINAPI Local_PXL_CRALocalCalculator_Create(
    		LPXLOPER XL_generalDatas,
			LPXLOPER XL_callDatas,
			LPXLOPER XL_fundDatas,
			LPXLOPER XL_cpnDatas,
			LPXLOPER XL_notionalCurve,
			LPXLOPER XL_callFeesCurve,
			LPXLOPER XL_fundSpreadCurve,
			LPXLOPER XL_cpnSpreadCurve,
			LPXLOPER XL_boostedFixCurve,
			LPXLOPER XL_barrierDownCurve,
			LPXLOPER XL_barrierUpCurve,
			LPXLOPER XL_localModelParams,
			LPXLOPER XL_mrsBeta,
			LPXLOPER XL_calibSecPFParams,
			LPXLOPER XL_nbSteps,
			LPXLOPER XL_flagToGenerateOSWATM,
			LPXLOPER XL_mktDataManager,
			LPXLOPER XL_productsToPrice,
			LPXLOPER XL_isSfrmStdCalib)
{
	ADD_LOG("Local_PXL_CRALocalCalculator_Create");
	bool PersistentInXL = false;
	return Local_CRALocalCalculator_Common(
			XL_generalDatas,
			XL_callDatas,
			XL_fundDatas,
			XL_cpnDatas,
			XL_notionalCurve,
			XL_callFeesCurve,
			XL_fundSpreadCurve,
			XL_cpnSpreadCurve,
			XL_boostedFixCurve,
			XL_barrierDownCurve,
			XL_barrierUpCurve,
			XL_localModelParams,
			XL_mrsBeta,
			XL_calibSecPFParams,
			XL_nbSteps,
			XL_flagToGenerateOSWATM,
			XL_mktDataManager,
			XL_productsToPrice,
			XL_isSfrmStdCalib,
			PersistentInXL);
}


//Local CRA Calculator
class localcraGetDataFunc : public ARMResultLong2LongFunc
{
public:
	localcraGetDataFunc(
        long craId,
        const string& getType)
    :
    C_craId(craId),
    C_getType(getType)
    {};
	
	long operator()( ARM_result& result, long objId ){
		return ARMLOCAL_LocalCRA_Get(
            C_craId,
            C_getType,
            result,
            objId);
    }

private:
	long    C_craId;
	string  C_getType;
};

class localcraSetDataFunc : public ARMResultLong2LongFunc
{
public:
	localcraSetDataFunc(
        long craId,
        long dataToSetId,
        const string& setPortfolioType,
		const vector< string >& mktDataKeys,
        bool isUpdated)
    :
    C_craId(craId),
    C_dataToSetId(dataToSetId),
    C_setPortfolioType(setPortfolioType),
	C_mktDataKeys(mktDataKeys),
    C_isUpdated(isUpdated)
    {};
	
	long operator()( ARM_result& result, long objId ){
		return ARMLOCAL_LocalCRA_Set(
            C_craId,
            C_dataToSetId,
            C_setPortfolioType,
			C_mktDataKeys,
            C_isUpdated,
            result,
            objId);
	}

private:
	long    C_craId;
	long    C_dataToSetId;
    string  C_setPortfolioType;
	vector< string > C_mktDataKeys;
    bool    C_isUpdated;
};


LPXLOPER Local_LocalCRAGet_Common(	LPXLOPER XL_craId,
									LPXLOPER XL_getType,
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
		XL_GETOBJID( XL_craId, craId, " ARM_ERR: Local CRA Calculator: Object expected", C_result);

		CCString getTypeStr;
		XL_readStrCell(XL_getType, getTypeStr, " ARM_ERR: Accessor Type: string expected", C_result);
        char* type = getTypeStr.c_str(); // à cause du new !!
		string getType(type);
        delete type;

		CCString craGetClass(GCGetTypeToClass(getType,craId).c_str());

		localcraGetDataFunc ourFunc(craId,getType);

		/// call the general function
		fillXL_Result( craGetClass, ourFunc, C_result, XL_result, PersistentInXL );
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_LocalCRAGet_Common" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}
							 
LPXLOPER Local_LocalCRASet_Common(
	LPXLOPER XL_craId,
	LPXLOPER XL_dataId,
	LPXLOPER XL_setPortfolioType,
	LPXLOPER XL_MktDataKeys,
	bool isUpdated,
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
		XL_GETOBJID( XL_craId, craId, " ARM_ERR: Local CRA Calculator: Object expected", C_result);

		long dataId;
		XL_GETOBJID( XL_dataId, dataId,	" ARM_ERR: Object expected", C_result);

		CCString setPortfolioTypeStr;
		XL_readStrCellWD(XL_setPortfolioType, setPortfolioTypeStr, "DEFAULT", " ARM_ERR: Portfolio Type: string expected", C_result);
        char * portType = setPortfolioTypeStr.c_str(); // à cause du new !!
		string setPortfolioType(portType);
        delete portType;

		VECTOR<CCString> mktDataKeys;
		VECTOR<CCString> mktDataKeysDef(0);
		XL_readStrVectorWD(XL_MktDataKeys, mktDataKeys, mktDataKeysDef, " ARM_ERR: Market datas keys: array of string expected", DOUBLE_TYPE, C_result);

		vector<string> mktDataKeysSTL(mktDataKeys.size());
		for (size_t i=0; i<mktDataKeys.size(); ++i)
        {
            mktDataKeys[i].toUpper();
			mktDataKeysSTL[i] = CCSTringToSTLString(mktDataKeys[i]);
        }

        localcraSetDataFunc ourFunc( craId,
									 dataId,
									 setPortfolioType,
									 mktDataKeysSTL,
									 isUpdated);

        if (isUpdated)
        {
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
        else
        {
            /// General call through the functor with an object creation		    
		    fillXL_Result( LOCAL_GC_CRA_CLASS, ourFunc, C_result, XL_result, PersistentInXL );
        }
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_LocalCRASet_Common" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


__declspec(dllexport) LPXLOPER WINAPI Local_LocalCRACalculator_GetData(
	LPXLOPER XL_craId,
	LPXLOPER XL_getType)
{
	ADD_LOG("Local_LocalCRACalculator_GetData");
	bool PersistentInXL = true;
	return Local_LocalCRAGet_Common(
	    XL_craId,
	    XL_getType,
        PersistentInXL );
}

__declspec(dllexport) LPXLOPER WINAPI Local_LocalCRACalculator_SetData(
	LPXLOPER XL_craId,
	LPXLOPER XL_dataId,
    LPXLOPER XL_setPortfolioType,
	LPXLOPER XL_mktDataKeys)
{
	ADD_LOG("Local_LocalCRACalculator_SetData");
	bool PersistentInXL = true;
    bool isUpdated = false;
	return Local_LocalCRASet_Common(
	    XL_craId,
	    XL_dataId,
        XL_setPortfolioType,
		XL_mktDataKeys,
        isUpdated,
        PersistentInXL );
}

__declspec(dllexport) LPXLOPER WINAPI Local_PXL_LocalCRACalculator_GetData(
	LPXLOPER XL_craId,
	LPXLOPER XL_getType)
{
	ADD_LOG("Local_PXL_LocalCRACalculator_GetData");
	bool PersistentInXL = false;
	return Local_LocalCRAGet_Common(
	    XL_craId,
	    XL_getType,
        PersistentInXL );
}
							 

__declspec(dllexport) LPXLOPER WINAPI Local_PXL_LocalCRACalculator_SetData(
	LPXLOPER XL_craId,
	LPXLOPER XL_dataId,
    LPXLOPER XL_setPortfolioType,
	LPXLOPER XL_mktDataKeys)
{
	ADD_LOG("Local_PXL_LocalCRACalculator_SetData");
	bool PersistentInXL = false;
    bool isUpdated = false;
	return Local_LocalCRASet_Common(
	    XL_craId,
	    XL_dataId,
        XL_setPortfolioType,
		XL_mktDataKeys,
        isUpdated,
        PersistentInXL );
}




///----------------------------------------------
///----------------------------------------------
///             CRA Calculator Accessor
/// Inputs :
///     
///----------------------------------------------
///----------------------------------------------
class craGetDataFunc : public ARMResultLong2LongFunc
{
public:
	craGetDataFunc(
        long craId,
        const string& getType)
    :
    C_craId(craId),
    C_getType(getType)
    {};
	
	long operator()( ARM_result& result, long objId ){
		return ARMLOCAL_CRA_Get(
            C_craId,
            C_getType,
            result,
            objId);
    }

private:
	long    C_craId;
	string  C_getType;
};

class craSetDataFunc : public ARMResultLong2LongFunc
{
public:
	craSetDataFunc(
        long craId,
        long dataToSetId,
        const string& setPortfolioType,
		const vector< string >& mktDataKeys,
        bool isUpdated)
    :
    C_craId(craId),
    C_dataToSetId(dataToSetId),
    C_setPortfolioType(setPortfolioType),
	C_mktDataKeys(mktDataKeys),
    C_isUpdated(isUpdated)
    {};
	
	long operator()( ARM_result& result, long objId ){
		return ARMLOCAL_CRA_Set(
            C_craId,
            C_dataToSetId,
            C_setPortfolioType,
			C_mktDataKeys,
            C_isUpdated,
            result,
            objId);
	}

private:
	long    C_craId;
	long    C_dataToSetId;
    string  C_setPortfolioType;
	vector< string > C_mktDataKeys;
    bool    C_isUpdated;
};


/////////////////////////////////////////////////////////////
/// central function that does the creation of the XL function
/////////////////////////////////////////////////////////////
LPXLOPER Local_CRAGet_Common(
	LPXLOPER XL_craId,
	LPXLOPER XL_getType,
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
		XL_GETOBJID( XL_craId, craId, " ARM_ERR: CRA Calculator: Object expected", C_result);

		CCString getTypeStr;
		XL_readStrCell(XL_getType, getTypeStr, " ARM_ERR: Accessor Type: string expected", C_result);
        char* type = getTypeStr.c_str(); // à cause du new !!
		string getType(type);
        delete type;

		CCString craGetClass(GCGetTypeToClass(getType,craId).c_str());

		craGetDataFunc ourFunc(craId,getType);

		/// call the general function
		fillXL_Result( craGetClass, ourFunc, C_result, XL_result, PersistentInXL );
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_CRAGet_Common" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}
							 
LPXLOPER Local_CRASet_Common(
	LPXLOPER XL_craId,
	LPXLOPER XL_dataId,
	LPXLOPER XL_setPortfolioType,
	LPXLOPER XL_MktDataKeys,
	bool isUpdated,
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
		XL_GETOBJID( XL_craId, craId, " ARM_ERR: CRA Calculator: Object expected", C_result);

		long dataId;
		XL_GETOBJID( XL_dataId, dataId,	" ARM_ERR: Object expected", C_result);

		CCString setPortfolioTypeStr;
		XL_readStrCellWD(XL_setPortfolioType, setPortfolioTypeStr, "DEFAULT", " ARM_ERR: Portfolio Type: string expected", C_result);
        char * portType = setPortfolioTypeStr.c_str(); // à cause du new !!
		string setPortfolioType(portType);
        delete portType;

		VECTOR<CCString> mktDataKeys;
		VECTOR<CCString> mktDataKeysDef(0);
		XL_readStrVectorWD(XL_MktDataKeys, mktDataKeys, mktDataKeysDef, " ARM_ERR: Market datas keys: array of string expected", DOUBLE_TYPE, C_result);

		vector<string> mktDataKeysSTL(mktDataKeys.size());
		for (size_t i=0; i<mktDataKeys.size(); ++i)
        {
            mktDataKeys[i].toUpper();
			mktDataKeysSTL[i] = CCSTringToSTLString(mktDataKeys[i]);
        }

        craSetDataFunc ourFunc( craId,
								dataId,
								setPortfolioType,
								mktDataKeysSTL,
								isUpdated);

        if (isUpdated)
        {
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
        else
        {
            /// General call through the functor with an object creation		    
		    fillXL_Result( LOCAL_GC_CRA_CLASS, ourFunc, C_result, XL_result, PersistentInXL );
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


__declspec(dllexport) LPXLOPER WINAPI Local_CRACalculator_SetRecalibFlags(
	LPXLOPER XL_craId,
	LPXLOPER XL_mrsFlag,
	LPXLOPER XL_betaFlag)
{
	ADD_LOG("Local_CRACalculator_SetRecalibFlags");
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
		XL_GETOBJID( XL_craId, craId, " ARM_ERR: CRA Calculator: Object expected", C_result);

		double mrsFlag = 0.0;
		XL_readNumCell(XL_mrsFlag, mrsFlag, " ARM_ERR: MRS: integer (0/1) expected", C_result);
		if ((int)mrsFlag != 0 && (int)mrsFlag != 1)
			throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"MRS flag should be equal to 0 or 1");
		
		double betaFlag = 0.0;
		XL_readNumCell(XL_betaFlag, betaFlag, " ARM_ERR: Beta: integer (0/1) expected", C_result);
		if ((int)betaFlag != 0 && (int)betaFlag != 1)
			throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"Beta flag should be equal to 0 or 1");

		long retCode, objId;
		CCString curClass = LOCAL_GC_CRA_CLASS;
		CCString stringId = GetLastCurCellEnvValue ();

		if(!stringId)
		{
			retCode = ARMLOCAL_CRA_SetRecalibFlags(
					craId,
					(int)mrsFlag,
					(int)betaFlag,
					C_result);

			if (retCode == ARM_OK)
			{
				objId = C_result.getLong();

				LocalSetCurCellEnvValue (curClass, objId);

				stringId = LocalMakeObjectId (objId, curClass);
			}
		}
		else
		{
			CCString prevClass = LocalGetStringObjectClass (stringId);

			objId = LocalGetNumObjectId (stringId);
				
			if(curClass == prevClass)
			{
				retCode = ARMLOCAL_CRA_SetRecalibFlags(
						craId,
						(int)mrsFlag,
						(int)betaFlag,
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

				retCode = ARMLOCAL_CRA_SetRecalibFlags(
						craId,
						(int)mrsFlag,
						(int)betaFlag,
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_CRACalculator_SetRecalibFlags" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}

__declspec(dllexport) LPXLOPER WINAPI Local_PXL_CRACalculator_SetRecalibFlags(
	LPXLOPER XL_craId,
	LPXLOPER XL_mrsFlag,
	LPXLOPER XL_betaFlag)
{
	ADD_LOG("Local_PXL_CRACalculator_SetRecalibFlags");
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
		XL_GETOBJID( XL_craId, craId, " ARM_ERR: CRA Calculator: Object expected", C_result);

		double mrsFlag = 0.0;
		XL_readNumCell(XL_mrsFlag, mrsFlag, " ARM_ERR: MRS: integer (0/1) expected", C_result);
		if ((int)mrsFlag != 0 && (int)mrsFlag != 1)
			throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"MRS flag should be equal to 0 or 1");
		
		double betaFlag = 0.0;
		XL_readNumCell(XL_betaFlag, betaFlag, " ARM_ERR: Beta: integer (0/1) expected", C_result);
		if ((int)betaFlag != 0 && (int)betaFlag != 1)
			throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"Beta flag should be equal to 0 or 1");

		long retCode, objId;
		CCString curClass = LOCAL_GC_CRA_CLASS;
		CCString stringId;

		retCode = ARMLOCAL_CRA_SetRecalibFlags(
						craId,
						(int)mrsFlag,
						(int)betaFlag,
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_CRACalculator_SetRecalibFlags" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


//CRA Calculator
__declspec(dllexport) LPXLOPER WINAPI Local_CRACalculator_GetData(
	LPXLOPER XL_craId,
	LPXLOPER XL_getType)
{
	ADD_LOG("Local_CRACalculator_GetData");
	bool PersistentInXL = true;
	return Local_CRAGet_Common(
	    XL_craId,
	    XL_getType,
        PersistentInXL );
}

__declspec(dllexport) LPXLOPER WINAPI Local_CRACalculator_SetData(
	LPXLOPER XL_craId,
	LPXLOPER XL_dataId,
    LPXLOPER XL_setPortfolioType,
	LPXLOPER XL_mktDataKeys)
{
	ADD_LOG("Local_CRACalculator_SetData");
	bool PersistentInXL = true;
    bool isUpdated = false;
	return Local_CRASet_Common(
	    XL_craId,
	    XL_dataId,
        XL_setPortfolioType,
		XL_mktDataKeys,
        isUpdated,
        PersistentInXL );
}

///////////////////////////////////
/// PXL version 
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_PXL_CRACalculator_GetData(
	LPXLOPER XL_craId,
	LPXLOPER XL_getType)
{
	ADD_LOG("Local_PXL_CRACalculator_GetData");
	bool PersistentInXL = false;
	return Local_CRAGet_Common(
	    XL_craId,
	    XL_getType,
        PersistentInXL );
}

__declspec(dllexport) LPXLOPER WINAPI Local_PXL_CRACalculator_SetData(
	LPXLOPER XL_craId,
	LPXLOPER XL_dataId,
    LPXLOPER XL_setPortfolioType,
	LPXLOPER XL_mktDataKeys)
{
	ADD_LOG("Local_PXL_CRACalculator_SetData");
	bool PersistentInXL = false;
    bool isUpdated = false;
	return Local_CRASet_Common(
	    XL_craId,
	    XL_dataId,
        XL_setPortfolioType,
		XL_mktDataKeys,
        isUpdated,
        PersistentInXL );
}





__declspec(dllexport) LPXLOPER WINAPI Local_CRACalculator_Create(
			LPXLOPER XL_generalDatas,
			LPXLOPER XL_callDatas,
			LPXLOPER XL_fundDatas,
			LPXLOPER XL_cpnDatas,
			LPXLOPER XL_notionalCurve,
			LPXLOPER XL_callFeesCurve,
			LPXLOPER XL_fundSpreadCurve,
			LPXLOPER XL_cpnSpreadCurve,
			LPXLOPER XL_boostedFixCurve,
			LPXLOPER XL_barrierDownCurve,
			LPXLOPER XL_barrierUpCurve,
			LPXLOPER XL_mrsBeta,
			LPXLOPER XL_calibSecPFParams,
			LPXLOPER XL_nbSteps,
			LPXLOPER XL_flagToGenerateOSWATM,
			LPXLOPER XL_mktDataManager,
			LPXLOPER XL_productsToPrice,
			LPXLOPER XL_isSfrmStdCalib)
{
	ADD_LOG("Local_CRACalculator_Create");
	bool PersistentInXL = true;

	return Local_CRACalculator_Common(
			XL_generalDatas,
			XL_callDatas,
			XL_fundDatas,
			XL_cpnDatas,
			XL_notionalCurve,
			XL_callFeesCurve,
			XL_fundSpreadCurve,
			XL_cpnSpreadCurve,
			XL_boostedFixCurve,
			XL_barrierDownCurve,
			XL_barrierUpCurve,
			XL_mrsBeta,
			XL_calibSecPFParams,
			XL_nbSteps,
			XL_flagToGenerateOSWATM,
			XL_mktDataManager,
			XL_productsToPrice,
			XL_isSfrmStdCalib,
			PersistentInXL);
}

__declspec(dllexport) LPXLOPER WINAPI Local_PXL_CRACalculator_Create(
    		LPXLOPER XL_generalDatas,
			LPXLOPER XL_callDatas,
			LPXLOPER XL_fundDatas,
			LPXLOPER XL_cpnDatas,
			LPXLOPER XL_notionalCurve,
			LPXLOPER XL_callFeesCurve,
			LPXLOPER XL_fundSpreadCurve,
			LPXLOPER XL_cpnSpreadCurve,
			LPXLOPER XL_boostedFixCurve,
			LPXLOPER XL_barrierDownCurve,
			LPXLOPER XL_barrierUpCurve,
			LPXLOPER XL_mrsBeta,
			LPXLOPER XL_calibSecPFParams,
			LPXLOPER XL_nbSteps,
			LPXLOPER XL_flagToGenerateOSWATM,
			LPXLOPER XL_mktDataManager,
			LPXLOPER XL_productsToPrice,
			LPXLOPER XL_isSfrmStdCalib)
{
	ADD_LOG("Local_PXL_CRACalculator_Create");
	bool PersistentInXL = false;
	return Local_CRACalculator_Common(
			XL_generalDatas,
			XL_callDatas,
			XL_fundDatas,
			XL_cpnDatas,
			XL_notionalCurve,
			XL_callFeesCurve,
			XL_fundSpreadCurve,
			XL_cpnSpreadCurve,
			XL_boostedFixCurve,
			XL_barrierDownCurve,
			XL_barrierUpCurve,
			XL_mrsBeta,
			XL_calibSecPFParams,
			XL_nbSteps,
			XL_flagToGenerateOSWATM,
			XL_mktDataManager,
			XL_productsToPrice,
			XL_isSfrmStdCalib,
			PersistentInXL);
}


__declspec(dllexport) LPXLOPER WINAPI Local_CRACalculator_CreateFromPf(
			LPXLOPER XL_optionPortfolio,
			LPXLOPER XL_mrsBeta,
			LPXLOPER XL_calibSecPFParams,
			LPXLOPER XL_nbSteps,
			LPXLOPER XL_flagToGenerateOSWATM,
			LPXLOPER XL_mktDataManager,
			LPXLOPER XL_productsToPrice,
			LPXLOPER XL_localModelParams,
			LPXLOPER XL_isSfrmStdCalib)
{
	ADD_LOG("Local_CRACalculator_CreateFromPf");
	bool PersistentInXL = true;
	return Local_CRACalculator_Common(
			XL_optionPortfolio,
			XL_mrsBeta,
			XL_calibSecPFParams,
			XL_nbSteps,
			XL_flagToGenerateOSWATM,
			XL_mktDataManager,
			XL_productsToPrice,
			XL_localModelParams,
			XL_isSfrmStdCalib,
			PersistentInXL);
}


__declspec(dllexport) LPXLOPER WINAPI Local_PXL_CRACalculator_CreateFromPf(
			LPXLOPER XL_optionPortfolio,
			LPXLOPER XL_mrsBeta,
			LPXLOPER XL_calibSecPFParams,
			LPXLOPER XL_nbSteps,
			LPXLOPER XL_flagToGenerateOSWATM,
			LPXLOPER XL_mktDataManager,
			LPXLOPER XL_productsToPrice,
			LPXLOPER XL_localModelParams,
			LPXLOPER XL_isSfrmStdCalib)
{
	ADD_LOG("Local_PXL_CRACalculator_CreateFromPf");
	bool PersistentInXL = false;
	return Local_CRACalculator_Common(
			XL_optionPortfolio,
			XL_mrsBeta,
			XL_calibSecPFParams,
			XL_nbSteps,
			XL_flagToGenerateOSWATM,
			XL_mktDataManager,
			XL_productsToPrice,
			XL_localModelParams,
			XL_isSfrmStdCalib,
			PersistentInXL);
}

/****************************************************************
					Global Cap Calculator
*****************************************************************/
class globalCapCalculatorFunc : public ARMResultLong2LongFunc
{
public:
	globalCapCalculatorFunc(ARM_Currency ccy,
						double startDate,
						double endDate,
						int payReceive,
						double notional,
						long notionalId,
						string fundIndexType,
						int fundDayCount,
						int fundFreq,
						int fundResetGap,
						int fundPayGap,
						int fundResetTiming,
						int fundPayTiming,
						string fundResetCal,
						string fundPayCal,
						int fundAdjRule,
						int fundIntRule,
						double fundLev,
						long fundLevId,
						ARM_Vector* globalCapParams,
						long pastFixingsId,
						double capLev,
						long capLevId,
						double capFixed,
						long capFixedId,
						double capStrike,
						long capStrikeId,
						double capSpread,
						long capSpreadId,
						int	nbSteps,
						vector<string> randGenerator,
						int samplerType,
						ARM_Vector* calibParams,
						vector<string> mdmKeys,
						long mktDataManagerId,
				        vector<string> productsToPrice)
    :
	C_currency(ccy),
    C_startDate(startDate),
    C_endDate(endDate),
	C_payReceive(payReceive),
	C_notional(notional),
	C_notionalId(notionalId),
	C_fundIndexType(fundIndexType),
	C_fundDayCount(fundDayCount),
	C_fundFreq(fundFreq),
	C_fundResetGap(fundResetGap),
	C_fundPayGap(fundPayGap),
	C_fundResetTiming(fundResetTiming),
	C_fundPayTiming(fundPayTiming),
	C_fundResetCal(fundResetCal),
	C_fundPayCal(fundPayCal),
	C_fundAdjRule(fundAdjRule),
	C_fundIntRule(fundIntRule),
	C_fundLev(fundLev),
	C_fundLevId(fundLevId),
	C_globalCapParams(globalCapParams),
	C_pastFixingsId(pastFixingsId),
	C_capLev(capLev),
	C_capLevId(capLevId),	
	C_capFixed(capFixed),
	C_capFixedId(capFixedId),
	C_capStrike(capStrike),
	C_capStrikeId(capStrikeId),
	C_capSpread(capSpread),
	C_capSpreadId(capSpreadId),
	C_nbSteps(nbSteps),
	C_randGenerator(randGenerator),
	C_samplerType(samplerType),
	C_calibParams(calibParams),
	C_mdmKeys(mdmKeys),
	C_mktDataManagerId(mktDataManagerId),
	C_productsToPrice(productsToPrice),
	C_globalCapId(-1),
	C_asOfDate(0.0)
    {
	};
	
	globalCapCalculatorFunc(
						long globalCapId,
						double fundLev,
						long fundLevId,
						double capLev,
						long capLevId,
						int	nbSteps,
						vector<string> randGenerator,
						int samplerType,
						ARM_Vector* calibParams,
						vector<string> mdmKeys,
						long mktDataManagerId,
				        vector<string> productsToPrice)
    :
	C_globalCapId(globalCapId),
	C_fundLev(fundLev),
	C_fundLevId(fundLevId),
	C_capLev(capLev),
	C_capLevId(capLevId),
	C_nbSteps(nbSteps),
	C_randGenerator(randGenerator),
	C_samplerType(samplerType),
	C_calibParams(calibParams),
	C_mdmKeys(mdmKeys),
	C_mktDataManagerId(mktDataManagerId),
	C_productsToPrice(productsToPrice),
	C_asOfDate(0.0)
    {
	};
	
	globalCapCalculatorFunc(
						double asOfDate,
						long globalCapId,
						double fundLev,
						long fundLevId,
						double capLev,
						long capLevId)
    :
	C_asOfDate(asOfDate),
	C_globalCapId(globalCapId),
	C_fundLev(fundLev),
	C_fundLevId(fundLevId),
	C_capLev(capLev),
	C_capLevId(capLevId)
    {
	};

	long operator()( ARM_result& result, long objId )
	{
		if (C_asOfDate != 0.0)
		{
			return ARMLOCAL_GlobalCapCalculator_Create_WithoutMktData(
				C_asOfDate,
				C_globalCapId,
				C_fundLev,
				C_fundLevId,
				C_capLev,
				C_capLevId,
				result,
				objId);
		}
		if (C_globalCapId != -1)
		{
			return ARMLOCAL_GlobalCapCalculator_Create(
				C_globalCapId,
				C_fundLev,
				C_fundLevId,
				C_capLev,
				C_capLevId,
				C_nbSteps,
				C_randGenerator,
				C_samplerType,
				C_calibParams,
				C_mdmKeys,
				C_mktDataManagerId,
				C_productsToPrice,
				result,
				objId);
		}
		else
		{
			return ARMLOCAL_GlobalCapCalculator_Create(
				C_currency,
				C_startDate,
				C_endDate,
				C_payReceive,
				C_notional,
				C_notionalId,
				C_fundIndexType,
				C_fundDayCount,
				C_fundFreq,
				C_fundResetGap,
				C_fundPayGap,
				C_fundResetTiming,
				C_fundPayTiming,
				C_fundResetCal,
				C_fundPayCal,
				C_fundAdjRule,
				C_fundIntRule,
				C_fundLev,
				C_fundLevId,
				C_capLev,
				C_capLevId,
				C_globalCapParams,
				C_pastFixingsId,
				C_capFixed,
				C_capFixedId,
				C_capStrike,
				C_capStrikeId,
				C_capSpread,
				C_capSpreadId,
				C_nbSteps,
				C_randGenerator,
				C_samplerType,
				C_calibParams,
				C_mdmKeys,
				C_mktDataManagerId,
				C_productsToPrice,
				result,
				objId);
		}
    }

private:
double				C_asOfDate;
int					C_globalCapId;
ARM_Currency		C_currency;
double				C_startDate;
double				C_endDate;
int					C_payReceive;
string				C_fundIndexType;
int					C_fundDayCount;
int					C_fundFreq;
int					C_fundResetGap;
int					C_fundPayGap;
int					C_fundResetTiming;
int					C_fundPayTiming;
string				C_fundResetCal;
string				C_fundPayCal;
int					C_fundAdjRule;
int					C_fundIntRule;
ARM_Vector*			C_globalCapParams;
long				C_pastFixingsId;

double				C_notional;
long				C_notionalId;
double				C_fundLev;
long				C_fundLevId;
double				C_capLev;
long				C_capLevId;
double				C_capFixed;
long				C_capFixedId;
double				C_capStrike;
long				C_capStrikeId;
double				C_capSpread;
long				C_capSpreadId;

int					C_nbSteps;
vector<string>		C_randGenerator;
int					C_samplerType;

ARM_Vector*			C_calibParams;
vector<string>		C_mdmKeys;
long				C_mktDataManagerId;
vector<string>		C_productsToPrice;
};

LPXLOPER Local_GlobalCapCalculator_Common(LPXLOPER XL_generalData,
									LPXLOPER XL_notionalCurve,
									LPXLOPER XL_fundData,
									LPXLOPER XL_fundLevCurve,
									LPXLOPER XL_globalCapParams,
									LPXLOPER XL_pastFixings,
									LPXLOPER XL_capLevCurve,
									LPXLOPER XL_capFixedCurve,
									LPXLOPER XL_capStrikeCurve,
									LPXLOPER XL_capSpreadCurve,
									LPXLOPER XL_modelParams,
									LPXLOPER XL_calibParams,
									LPXLOPER XL_mktDataManager,
									LPXLOPER XL_productsToPrice,
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
	
		//General Data
		//-------------
		VECTOR<CCString> generalData;
		XL_readStrVector (XL_generalData, generalData, "ARM_ERR: General data: array of string expected", DOUBLE_TYPE, C_result);
		if (generalData.size() != 4)
		{
			throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"General data: 4 parameters expected");
		}

		//Currency
		CCString ccyStr = generalData[0];
		string sCcy = CCSTringToSTLString(ccyStr);
		ARM_Currency ccy(sCcy.c_str());

		//Start & End
		CCString startDateXl = generalData[1];
		string sStartDate = CCSTringToSTLString(startDateXl);
		double startDate = atoi(sStartDate.c_str());
		CCString endDateXl = generalData[2];
		string sEndDate = CCSTringToSTLString(endDateXl);
		double endDate = atoi(sEndDate.c_str());
	    
		//PayReceive
		CCString payReceiveStr = generalData[3];
		long lPayReceive;
		if ((lPayReceive = ARM_ConvRecOrPay (payReceiveStr, C_result)) == ARM_DEFAULT_ERR)
		{
			throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"PayRec: 'P' or 'R' expected");
		}
		int payReceive = (int)lPayReceive;

		//Notional Value or Notional Curve
		double notional;
		CCString notionalStr;
		long notionalId;
		XL_readStrOrNumCell(XL_notionalCurve, notionalStr, notional, notionalId,
			   " ARM_ERR: notional: numeric or refValue Id expected",C_result);	

		if (notionalId == XL_TYPE_STRING)
			notionalId = LocalGetNumObjectId(notionalStr);
		else
			notionalId = ARM_NULL_OBJECT;

		//FundData
		VECTOR<CCString> fundData;
		XL_readStrVector (XL_fundData, fundData, " ARM_ERR: Fund data: array of string expected",DOUBLE_TYPE,C_result);
		int fundSize = fundData.size();
		if ((fundSize < 3) || (fundSize > 11))
		{
			throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"Funding data: 3 to 11 parameters expected");
		}

		//Fund Index Type
		string fundIndexType = CCSTringToSTLString(fundData[0]);

		//Fund Day Count
		CCString fundDayCountXl = fundData[1];
		int fundDayCount = (int)ARM_ConvDayCount (fundDayCountXl);

		//Fund Freq
		CCString fundFreqXl = fundData[2];
		long lFundFreq;
		if ((lFundFreq = ARM_ConvFrequency (fundFreqXl, C_result)) == ARM_DEFAULT_ERR)
		{
			throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"Funding data: invalid frequency");
		}
		int fundFreq = (int)lFundFreq;

		//Fund Reset Gap
		int fundResetGap = -ccy.GetSpotDays();
		if (fundSize > 3)
		{
			CCString fundResetGapXl = fundData[3];
			string sFundResetGap = CCSTringToSTLString(fundResetGapXl);
			fundResetGap = atoi(sFundResetGap.c_str());
		}
		
		//Fund Pay Gap
		int fundPayGap = 0;
		if (fundSize > 4)
		{
			CCString fundPayGapXl = fundData[4];
			string sFundPayGap = CCSTringToSTLString(fundPayGapXl);
			fundPayGap = atoi(sFundPayGap.c_str());
		}

		//Fund Reset Timing
		int fundResetTiming = K_ADVANCE;
		if (fundSize > 5)
		{
			CCString fundResetTimingXl = fundData[5];
			fundResetTiming = (int)ARM_ConvPayResetRule(fundResetTimingXl);
		}

		//Fund Pay Timing
		int fundPayTiming = K_ARREARS;
		if (fundSize > 5)
		{
			CCString fundPayTimingXl = fundData[6];
			fundPayTiming = (int)ARM_ConvPayResetRule(fundPayTimingXl);
		}

		//Fund Reset Calendar
		char* resetCal = NULL;
		resetCal = ccy.GetResetCalName(ccy.GetVanillaIndexType());
		string fundResetCal(resetCal);
		if (fundSize > 6)
		{
			CCString fundResetCalXl = fundData[7];
			fundResetCal = CCSTringToSTLString(fundResetCalXl);
		}

		//Fund Pay Calendar
		char* payCal = NULL;
		payCal = ccy.GetPayCalName(ccy.GetVanillaIndexType());
		string fundPayCal(payCal);
		if (fundSize > 7)
		{
			CCString fundPayCalXl = fundData[8];
			fundPayCal = CCSTringToSTLString(fundPayCalXl);
		}

		//Fund Adjustment Rule
		int fundAdjRule = K_MOD_FOLLOWING;
		if (fundSize > 8)
		{
			CCString fundAdjRuleXl = fundData[9];
			fundAdjRule = (int)ARM_ConvFwdRule(fundAdjRuleXl);
		}

		//Fund Int Rule
		int fundIntRule = K_ADJUSTED;
		if (fundSize > 9)
		{
			CCString fundIntRuleXl = fundData[10];
			fundIntRule = (int)ARM_ConvIntRule(fundIntRuleXl);
		}

		//FundLev Value or FundLev Curve
		double fundLev;
		double fundLevDef = 1.0;
		CCString fundLevStr;
		long fundLevId;
		XL_readStrOrNumCellWD(XL_fundLevCurve, fundLevStr, fundLev, fundLevDef, fundLevId,
			   " ARM_ERR: fundLev: numeric or refValue Id expected",C_result);	

		if (fundLevId == XL_TYPE_STRING)
			fundLevId = LocalGetNumObjectId(fundLevStr);
		else
			fundLevId = ARM_NULL_OBJECT;


		//Global Cap Params
		//-------------
		VECTOR<double> globalCapParamsVect;
		XL_readNumVector (XL_globalCapParams, globalCapParamsVect, "ARM_ERR: Global Cap Params: array of double expected", C_result);
		if (globalCapParamsVect.size() != 3)
		{
			throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"Global cap params: 3 parameters expected");
		}
		ARM_Vector* globalCapParams = CreateARMVectorFromVECTOR(globalCapParamsVect);
		//Test if nbPeriod is bigger then 0
		if ((int)globalCapParams->Elt(0) <= 0)
			throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT," Global Cap Params : nbPeriod must be > 0\n");

		//Past fixings
		//-------------
		CCString C_pastFixings_str;
		long pastFixingsId;
		
		XL_readStrCellWD (XL_pastFixings, C_pastFixings_str, "DEFAULT", "ARM_ERR: Past Fixings: Reference Value expected", C_result);
		pastFixingsId = (long) LocalGetNumObjectId(C_pastFixings_str);

		//CapLev Value or CapLev Curve
		double capLev;
		double capLevDef = 0.0;
		CCString capLevStr;
		long capLevId;
		XL_readStrOrNumCellWD(XL_capLevCurve, capLevStr, capLev, capLevDef, capLevId, "ARM_ERR: capLev: numeric or refValue Id expected",C_result);	

		if (capLevId == XL_TYPE_STRING)
			capLevId = LocalGetNumObjectId(capLevStr);
		else
			capLevId = ARM_NULL_OBJECT;
		
		//CapFixed Value or CapFixed Curve
		double capFixed;
		CCString capFixedStr;
		long capFixedId;
		XL_readStrOrNumCell(XL_capFixedCurve, capFixedStr, capFixed, capFixedId,
			   " ARM_ERR: capFixed: numeric or refValue Id expected",C_result);	

		if (capFixedId == XL_TYPE_STRING)
			capFixedId = LocalGetNumObjectId(capFixedStr);
		else
			capFixedId = ARM_NULL_OBJECT;

		//CapStrike Value or CapStrike Curve
		double capStrike;
		CCString capStrikeStr;
		long capStrikeId;
		XL_readStrOrNumCell(XL_capStrikeCurve, capStrikeStr, capStrike, capStrikeId,
			   " ARM_ERR: capStrike: numeric or refValue Id expected",C_result);	

		if (capStrikeId == XL_TYPE_STRING)
			capStrikeId = LocalGetNumObjectId(capStrikeStr);
		else
			capStrikeId = ARM_NULL_OBJECT;

		//CapSpread Value or CapSpread Curve
		double capSpread;
		CCString capSpreadStr;
		long capSpreadId;
		XL_readStrOrNumCell(XL_capSpreadCurve, capSpreadStr, capSpread, capSpreadId,
			   " ARM_ERR: capSpread: numeric or refValue Id expected",C_result);	

		if (capSpreadId == XL_TYPE_STRING)
			capSpreadId = LocalGetNumObjectId(capSpreadStr);
		else
			capSpreadId = ARM_NULL_OBJECT;

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

		//Calib Params
		//-------------
		VECTOR<double> calibParamsVect;
		XL_readNumVector (XL_calibParams, calibParamsVect, "ARM_ERR: Calib Params: array of double expected", C_result);
		if (calibParamsVect.size() != 5)
		{
			throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"Calib Params: 5 parameters expected");
		}
		ARM_Vector* calibParams = CreateARMVectorFromVECTOR(calibParamsVect);

		//MktDataManager: persistent object
		VECTOR<CCString> mktDataManager;
		XL_readStrVector(XL_mktDataManager,mktDataManager," ARM_ERR: Market data: array of string expected",DOUBLE_TYPE,C_result);
        
		if (mktDataManager.size() != 4)
		{
			throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"Market Data: MDM + 3 model parameters expected");
		}
		
		int i = 0;
		long mktDataManagerId = LocalGetNumObjectId (mktDataManager[0]);
		vector<string> mdmKeys(mktDataManager.size()-1);		
		for (i=1; i<mktDataManager.size(); ++i)
			mdmKeys[i-1] = CCSTringToSTLString(mktDataManager[i]);
	
		//Pricing flags
		VECTOR<CCString> pricingFlags;

		XL_readStrVector(XL_productsToPrice,pricingFlags," ARM_ERR: Pricing flags: array of string expected",DOUBLE_TYPE,C_result);
        if (pricingFlags.size()>9)
		{
			throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"maximum: 9 products to price");
		}

		vector< string > productsToPrice(pricingFlags.size(), string("N"));
		for (i=0; i<pricingFlags.size(); i++)
			productsToPrice[i] = CCSTringToSTLString(pricingFlags[i]);

		//Functor call
		globalCapCalculatorFunc ourFunc(ccy,
										startDate,
										endDate,
										payReceive,
										notional,
										notionalId,
										fundIndexType,
										fundDayCount,
										fundFreq,
										fundResetGap,
										fundPayGap,
										fundResetTiming,
										fundPayTiming,
										fundResetCal,
										fundPayCal,
										fundAdjRule,
										fundIntRule,
										fundLev,
										fundLevId,
										globalCapParams,
										pastFixingsId,
										capLev,
										capLevId,
										capFixed,
										capFixedId,
										capStrike,
										capStrikeId,
										capSpread,
										capSpreadId,
										nbSteps,
										randGenerator,
										samplerType,
										calibParams,
										mdmKeys,
										mktDataManagerId,
										productsToPrice);
	
		/// call the general function
		fillXL_Result( LOCAL_GC_GLOBALCAP_CLASS, ourFunc, C_result, XL_result, PersistentInXL );	
		
		delete resetCal;
		delete payCal;
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_GlobalCapCalculator_Common" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}

LPXLOPER Local_GlobalCapCalculator_Common(
			LPXLOPER XL_globalCap,
			LPXLOPER XL_fundLevCurve,
			LPXLOPER XL_capLevCurve,
			LPXLOPER XL_modelParams,
			LPXLOPER XL_calibParams,
			LPXLOPER XL_mktDataManager,
			LPXLOPER XL_productsToPrice,
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
	
		CCString C_globalCap;

		XL_readStrCell( XL_globalCap, C_globalCap, " ARM_ERR: Global Cap Id: Object expected",C_result);
		long C_globalCapId = LocalGetNumObjectId(C_globalCap);

		//FundLev Value or FundLev Curve
		double fundLev;
		double fundLevDef = 1.0;
		CCString fundLevStr;
		long fundLevId;
		XL_readStrOrNumCellWD(XL_fundLevCurve, fundLevStr, fundLev, fundLevDef, fundLevId,
			   " ARM_ERR: fundLev: numeric or refValue Id expected",C_result);	

		if (fundLevId == XL_TYPE_STRING)
			fundLevId = LocalGetNumObjectId(fundLevStr);
		else
			fundLevId = ARM_NULL_OBJECT;

		//CapLev Value or CapLev Curve
		double capLev;
		double capLevDef = 0.0;
		CCString capLevStr;
		long capLevId;
		XL_readStrOrNumCellWD(XL_capLevCurve, capLevStr, capLev, capLevDef, capLevId,
			   " ARM_ERR: capLev: numeric or refValue Id expected",C_result);	

		if (capLevId == XL_TYPE_STRING)
			capLevId = LocalGetNumObjectId(capLevStr);
		else
			capLevId = ARM_NULL_OBJECT;
		
	
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

		//Calib Params
		//-------------
		VECTOR<double> calibParamsVect;
		XL_readNumVector (XL_calibParams, calibParamsVect, "ARM_ERR: Calib Params: array of double expected", C_result);
		if (calibParamsVect.size() != 5)
		{
			throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"Calib Params: 5 parameters expected");
		}
		ARM_Vector* calibParams = CreateARMVectorFromVECTOR(calibParamsVect);

		//MktDataManager: persistent object
		VECTOR<CCString> mktDataManager;
		XL_readStrVector(XL_mktDataManager,mktDataManager," ARM_ERR: Market data: array of string expected",DOUBLE_TYPE,C_result);
        
		if (mktDataManager.size()<4)
		{
			throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"Market Data: MDM + 3 model parameters expected");
		}

		long mktDataManagerId = LocalGetNumObjectId (mktDataManager[0]);
		vector<string> mdmKeys(mktDataManager.size()-1);
		
		int i;

		for (i=1; i<mktDataManager.size(); ++i)
			mdmKeys[i-1] = CCSTringToSTLString(mktDataManager[i]);
	
		//Pricing flags
		VECTOR<CCString> pricingFlags;

		XL_readStrVector(XL_productsToPrice,pricingFlags," ARM_ERR: Pricing flags: array of string expected",DOUBLE_TYPE,C_result);
        if (pricingFlags.size()>9)
		{
			throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"maximum: 9 products to price");
		}

		vector< string > productsToPrice(pricingFlags.size(), string("N"));
		for (i=0; i<pricingFlags.size(); i++)
			productsToPrice[i] = CCSTringToSTLString(pricingFlags[i]);

		//Functor call
		globalCapCalculatorFunc ourFunc(C_globalCapId,
										fundLev,
										fundLevId,
										capLev,
										capLevId,
										nbSteps,
										randGenerator,
										samplerType,
										calibParams,
										mdmKeys,
										mktDataManagerId,
										productsToPrice);
		
		/// call the general function
		fillXL_Result( LOCAL_GC_GLOBALCAP_CLASS, ourFunc, C_result, XL_result, PersistentInXL );	
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_GlobalCapCalculator_Common" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


LPXLOPER Local_GlobalCapCalculator_Common_WithoutMktData(
			LPXLOPER XL_asOfDate,
			LPXLOPER XL_globalCap,
			LPXLOPER XL_fundLevCurve,
			LPXLOPER XL_capLevCurve,
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
		
		double C_asOfDate;
		CCString C_globalCap;
		double fundLev;
		double fundLevDef = 1.0;
		CCString fundLevStr;
		long fundLevId;
		double capLev;
		double capLevDef = 0.0;
		CCString capLevStr;
		long capLevId;

		CCString C_pastFixings_str;
		
		XL_readNumCell( XL_asOfDate, C_asOfDate, " ARM_ERR: As Of Date : Double expected",C_result);		
		XL_readStrCell( XL_globalCap, C_globalCap, " ARM_ERR: globalCap Id: Object expected",C_result);
		XL_readStrOrNumCellWD(XL_fundLevCurve, fundLevStr, fundLev, fundLevDef, fundLevId, " ARM_ERR: fundLev: numeric or refValue Id expected",C_result);	
		XL_readStrOrNumCellWD(XL_capLevCurve, capLevStr, capLev, capLevDef, capLevId, " ARM_ERR: capLev: numeric or refValue Id expected",C_result);	

		long C_globalCapId = LocalGetNumObjectId(C_globalCap);
		
		if (fundLevId == XL_TYPE_STRING)
			fundLevId = LocalGetNumObjectId(fundLevStr);
		else
			fundLevId = ARM_NULL_OBJECT;

		if (capLevId == XL_TYPE_STRING)
			capLevId = LocalGetNumObjectId(capLevStr);
		else
			capLevId = ARM_NULL_OBJECT;

		//Functor call
		globalCapCalculatorFunc ourFunc(C_asOfDate,
										C_globalCapId,
										fundLev,
										fundLevId,
										capLev,
										capLevId);
		
		/// call the general function
		fillXL_Result( LOCAL_GC_GLOBALCAP_CLASS, ourFunc, C_result, XL_result, PersistentInXL );	
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_GlobalCapCalculator_Common" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}

__declspec(dllexport) LPXLOPER WINAPI Local_GlobalCapCalculator_Create(
			LPXLOPER XL_generalData,
			LPXLOPER XL_notionalCurve,
			LPXLOPER XL_fundData,
			LPXLOPER XL_fundLevCurve,
			LPXLOPER XL_globalCapParams,
			LPXLOPER XL_pastFixings,
			LPXLOPER XL_capLevCurve,
			LPXLOPER XL_capFixedCurve,
			LPXLOPER XL_capStrikeCurve,
			LPXLOPER XL_capSpreadCurve,
			LPXLOPER XL_modelParams,
			LPXLOPER XL_calibParams,
			LPXLOPER XL_mktDataManager,
			LPXLOPER XL_productsToPrice)
{
	ADD_LOG("Local_GlobalCapCalculator_Create");
	bool PersistentInXL = true;

	return Local_GlobalCapCalculator_Common(
			XL_generalData,
			XL_notionalCurve,
			XL_fundData,
			XL_fundLevCurve,
			XL_globalCapParams,
			XL_pastFixings,
			XL_capLevCurve,
			XL_capFixedCurve,
			XL_capStrikeCurve,
			XL_capSpreadCurve,
			XL_modelParams,
			XL_calibParams,
			XL_mktDataManager,
			XL_productsToPrice,
			PersistentInXL);
}

__declspec(dllexport) LPXLOPER WINAPI Local_PXL_GlobalCapCalculator_Create(
			LPXLOPER XL_generalData,
			LPXLOPER XL_notionalCurve,
			LPXLOPER XL_fundData,
			LPXLOPER XL_fundLevCurve,
			LPXLOPER XL_globalCapParams,
			LPXLOPER XL_pastFixings,
			LPXLOPER XL_capLevCurve,
			LPXLOPER XL_capFixedCurve,
			LPXLOPER XL_capStrikeCurve,
			LPXLOPER XL_capSpreadCurve,
			LPXLOPER XL_modelParams,
			LPXLOPER XL_calibParams,
			LPXLOPER XL_mktDataManager,
			LPXLOPER XL_productsToPrice)
{
	ADD_LOG("Local_PXL_GlobalCapCalculator_Create");
	bool PersistentInXL = false;

	return Local_GlobalCapCalculator_Common(
			XL_generalData,
			XL_notionalCurve,
			XL_fundData,
			XL_fundLevCurve,
			XL_globalCapParams,
			XL_pastFixings,
			XL_capLevCurve,
			XL_capFixedCurve,
			XL_capStrikeCurve,
			XL_capSpreadCurve,
			XL_modelParams,
			XL_calibParams,
			XL_mktDataManager,
			XL_productsToPrice,
			PersistentInXL);
}

__declspec(dllexport) LPXLOPER WINAPI Local_GlobalCapCalculator_CreateFromPf(
			LPXLOPER XL_globalCap,
			LPXLOPER XL_fundLevCurve,
			LPXLOPER XL_capLevCurve,
			LPXLOPER XL_modelParams,
			LPXLOPER XL_calibParams,
			LPXLOPER XL_mktDataManager,
			LPXLOPER XL_productsToPrice)
{
	ADD_LOG("Local_GlobalCapCalculator_CreateFromPf");
	bool PersistentInXL = true;

	return Local_GlobalCapCalculator_Common(
			XL_globalCap,
			XL_fundLevCurve,
			XL_capLevCurve,
			XL_modelParams,
			XL_calibParams,
			XL_mktDataManager,
			XL_productsToPrice,
			PersistentInXL);
}

__declspec(dllexport) LPXLOPER WINAPI Local_PXL_GlobalCapCalculator_CreateFromPf(
			LPXLOPER XL_globalCap,
			LPXLOPER XL_fundLevCurve,
			LPXLOPER XL_capLevCurve,
			LPXLOPER XL_modelParams,
			LPXLOPER XL_calibParams,
			LPXLOPER XL_mktDataManager,
			LPXLOPER XL_productsToPrice)
{
	ADD_LOG("Local_PXL_GlobalCapCalculator_CreateFromPf");
	bool PersistentInXL = false;

	return Local_GlobalCapCalculator_Common(
			XL_globalCap,
			XL_fundLevCurve,
			XL_capLevCurve,
			XL_modelParams,
			XL_calibParams,
			XL_mktDataManager,
			XL_productsToPrice,
			PersistentInXL);
}

__declspec(dllexport) LPXLOPER WINAPI Local_GlobalCapCalculator_CreateFromPfWithoutMktData(
			LPXLOPER XL_asOf,
			LPXLOPER XL_globalCap,
			LPXLOPER XL_fundLevCurve,
			LPXLOPER XL_capLevCurve)
{
	ADD_LOG("Local_GlobalCapCalculator_CreateFromPfWithoutMktData");
	bool PersistentInXL = false;

	return Local_GlobalCapCalculator_Common_WithoutMktData(
			XL_asOf,
			XL_globalCap,
			XL_fundLevCurve,
			XL_capLevCurve,
			PersistentInXL);
}


class globalCapSetDataFunc : public ARMResultLong2LongFunc
{
public:
	globalCapSetDataFunc(
        long globalCapId,
        long dataToSetId,
        const string& setPortfolioType,
		const vector< string >& mktDataKeys,
        bool isUpdated)
    :
    C_globalCapId(globalCapId),
    C_dataToSetId(dataToSetId),
    C_setPortfolioType(setPortfolioType),
	C_mktDataKeys(mktDataKeys),
    C_isUpdated(isUpdated)
    {};
	
	long operator()( ARM_result& result, long objId ){
		return ARMLOCAL_GlobalCap_Set(
            C_globalCapId,
            C_dataToSetId,
            C_setPortfolioType,
			C_mktDataKeys,
            C_isUpdated,
            result,
            objId);
	}

private:
	long    C_globalCapId;
	long    C_dataToSetId;
    string  C_setPortfolioType;
	vector< string > C_mktDataKeys;
    bool    C_isUpdated;
};

LPXLOPER Local_GlobalCap_Set_Common(
	LPXLOPER XL_globalCapId,
	LPXLOPER XL_dataId,
	LPXLOPER XL_setPortfolioType,
	LPXLOPER XL_MktDataKeys,
	bool isUpdated,
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

		long globalCapId;
		XL_GETOBJID( XL_globalCapId, globalCapId, " ARM_ERR: Global Cap Calculator: Object expected", C_result);

		long dataId;
		XL_GETOBJID( XL_dataId, dataId,	" ARM_ERR: Object expected", C_result);

		CCString setPortfolioTypeStr;
		XL_readStrCellWD(XL_setPortfolioType, setPortfolioTypeStr, "DEFAULT", " ARM_ERR: Portfolio Type: string expected", C_result);
        char * portType = setPortfolioTypeStr.c_str(); // à cause du new !!
		string setPortfolioType(portType);
        delete portType;

		VECTOR<CCString> mktDataKeys;
		VECTOR<CCString> mktDataKeysDef(0);
		XL_readStrVectorWD(XL_MktDataKeys, mktDataKeys, mktDataKeysDef, " ARM_ERR: Market data keys: array of string expected", DOUBLE_TYPE, C_result);

		vector<string> mktDataKeysSTL(mktDataKeys.size());
		for (size_t i=0; i<mktDataKeys.size(); ++i)
        {
            mktDataKeys[i].toUpper();
			mktDataKeysSTL[i] = CCSTringToSTLString(mktDataKeys[i]);
        }

		globalCapSetDataFunc ourFunc( globalCapId,
									 dataId,
									 setPortfolioType,
									 mktDataKeysSTL,
									 isUpdated);

        if (isUpdated)
        {
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
        else
        {
            /// General call through the functor with an object creation		    
		    fillXL_Result( LOCAL_GC_GLOBALCAP_CLASS, ourFunc, C_result, XL_result, PersistentInXL );
        }
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_GlobalCap_Set_Common" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}

__declspec(dllexport) LPXLOPER WINAPI Local_GlobalCapCalculator_SetData(
	LPXLOPER XL_calculatorId,
	LPXLOPER XL_dataId,
    LPXLOPER XL_setPortfolioType,
	LPXLOPER XL_mktDataKeys)
{
	ADD_LOG("Local_GlobalCapCalculator_SetData");
	bool PersistentInXL = true;
    bool isUpdated = false;
	return Local_GlobalCap_Set_Common(
	    XL_calculatorId,
	    XL_dataId,
        XL_setPortfolioType,
		XL_mktDataKeys,
        isUpdated,
        PersistentInXL );
}
		
__declspec(dllexport) LPXLOPER WINAPI Local_PXL_GlobalCapCalculator_SetData(
	LPXLOPER XL_calculatorId,
	LPXLOPER XL_dataId,
    LPXLOPER XL_setPortfolioType,
	LPXLOPER XL_mktDataKeys)
{
	ADD_LOG("Local_PXL_GlobalCapCalculator_SetData");
	bool PersistentInXL = false;
    bool isUpdated = false;
	return Local_GlobalCap_Set_Common(
	    XL_calculatorId,
	    XL_dataId,
        XL_setPortfolioType,
		XL_mktDataKeys,
        isUpdated,
        PersistentInXL );
}
		
///----------------------------------------------
///----------------------------------------------
///             TARN Calculator
/// Inputs :
///     
///----------------------------------------------
///----------------------------------------------
class tarnCalculatorFunc : public ARMResultLong2LongFunc
{
public:
	tarnCalculatorFunc(
        double startDate,
        double endDate,
        double strike,
        long strikeId,
        long payRec,
        long cpnDayCount,
        long cpnFreq,
        long cpnTiming,
        const string& cpnIndexTerm,
        long cpnIndexDayCount,
        const string& cpnResetCal,
        const string& cpnPayCal,
        long cpnResetGap,
		long intRule,
        double leverage,
        long leverageId,
		double cpnMin,
		long cpnMinId,
		double cpnMax,
		long cpnMaxId,
		double lifeTimeCapTarget,
		bool globalCapFlag,
		double lifeTimeFloorTarget,
        double fundSpread,
        long fundSpreadId,
        long  fundFreq,
        long  fundDayCount,
        double nominal,
        long nominalId,
		double fees,
        long feesId,
		long fundNomId,
		const vector<double>& nbIterations,
		const string& genType1,
		const string& genType2,
		int firstNbTimes,
		int firstNbDims,
		const string& pathScheme,
		const string& pathOrder,
		const string& antithetic,
        const vector< string >& calibFlags,
        const vector< string >& outputFlags,
        long mktDataManagerId,
        const vector< string >& keys,
		bool isCustomResetFlag,
		long resetDatesId,
        double asOfDate)
    :
    C_startDate(startDate),
    C_endDate(endDate),
    C_strike(strike),
    C_strikeId(strikeId),
    C_payRec(payRec),
    C_cpnDayCount(cpnDayCount),
    C_cpnFreq(cpnFreq),
    C_cpnTiming(cpnTiming),
    C_cpnIndexTerm(cpnIndexTerm),
    C_cpnIndexDayCount(cpnIndexDayCount),
    C_cpnResetCal(cpnResetCal),
    C_cpnPayCal(cpnPayCal),
    C_cpnResetGap(cpnResetGap),
	C_intRule(intRule),
    C_leverage(leverage),
	C_cpnMin(cpnMin),
	C_cpnMinId(cpnMinId),
	C_cpnMax(cpnMax),
	C_cpnMaxId(cpnMaxId),
    C_leverageId(leverageId),
	C_lifeTimeCapTarget(lifeTimeCapTarget),
	C_globalCapFlag(globalCapFlag),
	C_lifeTimeFloorTarget(lifeTimeFloorTarget),
    C_fundSpread(fundSpread),
    C_fundSpreadId(fundSpreadId),
    C_fundFreq(fundFreq),
    C_fundDayCount(fundDayCount),
    C_nominal(nominal),
    C_nominalId(nominalId),
	C_fees(fees),
    C_feesId(feesId),
    C_fundNomId(fundNomId),
	C_nbIterations(nbIterations),
	C_genType1(genType1),
	C_genType2(genType2),
	C_firstNbTimes(firstNbTimes),
	C_firstNbDims(firstNbDims),
	C_pathScheme(pathScheme),
	C_pathOrder(pathOrder),
	C_antithetic(antithetic),
    C_calibFlags(calibFlags),
    C_outputFlags(outputFlags),
    C_mktDataManagerId(mktDataManagerId),
    C_keys(keys),
	C_isCustomResetFlag(isCustomResetFlag),
	C_resetDatesId(resetDatesId),
    C_asOfDate(asOfDate)
    {};
	
	long operator()( ARM_result& result, long objId ){
		return ARMLOCAL_TARNCalculator_Create(
            C_startDate,
            C_endDate,
            C_strike,
            C_strikeId,
            C_payRec,
            C_cpnDayCount,
            C_cpnFreq,
            C_cpnTiming,
            C_cpnIndexTerm,
            C_cpnIndexDayCount,
            C_cpnResetCal,
            C_cpnPayCal,
            C_cpnResetGap,
			C_intRule,
            C_leverage,
            C_leverageId,
			C_cpnMin,
			C_cpnMinId,
			C_cpnMax,
			C_cpnMaxId,
			C_lifeTimeCapTarget,
			C_globalCapFlag,
			C_lifeTimeFloorTarget,
            C_fundSpread,
            C_fundSpreadId,
            C_fundFreq,
            C_fundDayCount,
            C_nominal,
            C_nominalId,
			C_fees,
            C_feesId,
            C_fundNomId,
			C_nbIterations,
			C_genType1,
			C_genType2,
			C_firstNbTimes,
			C_firstNbDims,
			C_pathScheme,
			C_pathOrder,
			C_antithetic,
            C_calibFlags,
            C_outputFlags,
            C_mktDataManagerId,
            C_keys,
			C_isCustomResetFlag,
			C_resetDatesId,
            C_asOfDate,
            result,
            objId);
    }

private:

    double              C_startDate;
    double              C_endDate;
    double              C_strike;
    long                C_strikeId;
    long                C_payRec;
    long                C_cpnDayCount;
	long				C_cpnFreq;
    long                C_fixDayCount;
    long                C_cpnTiming;
    string              C_cpnIndexTerm;
    long                C_cpnIndexDayCount;
    string              C_cpnResetCal;
    string              C_cpnPayCal;
    long                C_cpnResetGap;
	long				C_intRule;
    double              C_leverage;
    long                C_leverageId;
	double				C_cpnMin;
	long				C_cpnMinId;
	double				C_cpnMax;
	long				C_cpnMaxId;
	double				C_lifeTimeCapTarget;
	bool				C_globalCapFlag;
	double				C_lifeTimeFloorTarget;
    double              C_fundSpread;
    long                C_fundSpreadId;
    long                C_fundFreq;
    long                C_fundDayCount;
    double              C_nominal;
    long                C_nominalId;
	double              C_fees;
    long                C_feesId;
    long                C_fundNomId;
	vector<double>		C_nbIterations;
	string				C_genType1;
	string				C_genType2;
	int					C_firstNbTimes;
	int					C_firstNbDims;
	string				C_pathScheme;
	string				C_pathOrder;
	string				C_antithetic;
    vector< string >    C_calibFlags;
    vector< string >    C_outputFlags;
    long                C_mktDataManagerId;
    vector< string >    C_keys;
	bool				C_isCustomResetFlag;
	long				C_resetDatesId;
    double              C_asOfDate;
};


/////////////////////////////////////////////////////////////
/// central function that does the creation of the XL function
/////////////////////////////////////////////////////////////
LPXLOPER Local_TARNCalculator_Common(
    LPXLOPER XL_startDate,
    LPXLOPER XL_endDate,
    LPXLOPER XL_strike,
    LPXLOPER XL_payRec,
    LPXLOPER XL_cpnDatas,
    LPXLOPER XL_mktDatas,
    LPXLOPER XL_cpnResetGap,
	LPXLOPER XL_intRule,
    LPXLOPER XL_cpnCurves,
	LPXLOPER XL_lifeTimeCapDatas,
    LPXLOPER XL_lifeTimeFloorTarget,
    LPXLOPER XL_fundSpread,
    LPXLOPER XL_fundDatas,
    LPXLOPER XL_nominal,
	LPXLOPER XL_numDatas,
    LPXLOPER XL_calibFlags,
    LPXLOPER XL_outputFlags,
    LPXLOPER XL_fundNominal,
	LPXLOPER XL_fees,
    LPXLOPER XL_asOfDate,
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

    
		/// General TARN datas
		/// -----------------
		/// Start & End
		double startDate;
		XL_readNumCell(XL_startDate,startDate," ARM_ERR: start date: date expected",C_result);
		double endDate;
		XL_readNumCell(XL_endDate,endDate," ARM_ERR: end date: date expected",C_result);

		/// Fix rate or RF Strike
		double strike;
		CCString strikeStr;
		long     strikeId;
		XL_readStrOrNumCell(XL_strike, strikeStr, strike, strikeId,
			   " ARM_ERR: strike: numeric or refValue Id expected",C_result);	
		if(strikeId == XL_TYPE_STRING)
			strikeId = LocalGetNumObjectId(strikeStr);
		else
			strikeId = ARM_NULL_OBJECT;

		/// RF Leg Pay/Rec
		CCString payRecStr;
		long payRec;
		XL_readStrCell(XL_payRec,payRecStr," ARM_ERR: payer/receiver: string expected",C_result);
		if((payRec = ARM_ConvRecOrPay (payRecStr, C_result)) == ARM_DEFAULT_ERR)
		{
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}


		/// Required RF coupon datas
		/// ------------------------
		VECTOR<CCString> cpnDatas;
		XL_readStrVector (XL_cpnDatas,cpnDatas," ARM_ERR: Coupon datas: array of string expected",DOUBLE_TYPE,C_result);
		if(cpnDatas.size() < 2)
		{
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}

		/// RF coupon frequency (reset & payment)
		CCString cpnFreqXl=cpnDatas[0];
		long cpnFreq;
		if((cpnFreq = ARM_ConvFrequency (cpnFreqXl, C_result)) == ARM_DEFAULT_ERR)
		{
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}

		/// RF coupon day count
		CCString cpnDayCountXl=cpnDatas[1];
		long cpnDayCount = ARM_ConvDayCount (cpnDayCountXl);


		/// Market datas : curve , MRS and market models for OSW & CF pricings
		/// ------------------------------------------------------------------
		VECTOR<CCString> mktDatas;
        VECTOR<CCString> mktDatas_Def;
		
        XL_readStrVectorWD(XL_mktDatas, mktDatas, mktDatas_Def,
                           " ARM_ERR: Market datas: array of string expected", 
                           DOUBLE_TYPE, C_result);
        
        long mktDataManagerId = ARM_NULL_OBJECT;
        vector< string > keys;

        if ( mktDatas.size() >= 1 )
        {
		   mktDataManagerId = LocalGetNumObjectId(mktDatas[0]);
		
           vector< string > theKeys(mktDatas.size()-1);
        
           size_t i;
		
           for (i = 1; i < mktDatas.size(); ++i)
			   theKeys[i-1] = CCSTringToSTLString(mktDatas[i]);

           keys = theKeys;
        }

		//// -------------------------------
		///  ----- OPTIONAL ARGUMENTS ------
		//// -------------------------------

		/// RF coupon datas
		/// ---------------
		/// Advance/Arrears RF coupon reset
		CCString cpnTimingXl("ADV");
		if(cpnDatas.size() >= 3)
			cpnTimingXl=cpnDatas[2];
		long cpnTiming = ARM_ConvPayResetRule (cpnTimingXl);

		/// RF index term
		CCString cpnIndexTermXl(ARM_ArgConvReverse_MatFrequency.GetString(cpnFreq).c_str());
		if(cpnDatas.size() >= 4)
			cpnIndexTermXl=cpnDatas[3];
		string cpnIndexTerm=CCSTringToSTLString(cpnIndexTermXl);

		/// RF index day count
		CCString cpnIndexDayCountXl("A360");
		if(cpnDatas.size() >= 5)
			cpnIndexDayCountXl=cpnDatas[4];
		long cpnIndexDayCount = ARM_ConvDayCount (cpnIndexDayCountXl);

		/// RF reset calendar
		CCString cpnResetCalXl("EUR");
		if(cpnDatas.size() >= 6)
			cpnResetCalXl=cpnDatas[5];
		string cpnResetCal = CCSTringToSTLString(cpnResetCalXl);

		/// RF payment calendar
		CCString cpnPayCalXl("EUR");
		if(cpnDatas.size() >= 7)
			cpnPayCalXl=cpnDatas[6];
		string cpnPayCal = CCSTringToSTLString(cpnPayCalXl);

		/// RF reset gap & Custom Reset dates
		double cpnResetGapD;
		long cpnResetGap;
		CCString resetDatesStr;
		long     resetDatesId;
		bool isCustomResetFlag = false;
		XL_readStrOrNumCell(XL_cpnResetGap, resetDatesStr, cpnResetGapD, resetDatesId,
			 " ARM_ERR: resetDates: numeric or gp_vector Id expected",C_result);	
		//Custom reset dates
		if(resetDatesId == XL_TYPE_STRING)
		{
			resetDatesId = LocalGetNumObjectId(resetDatesStr);
			isCustomResetFlag = true;
			cpnResetGap = 0;
		}
		//Non custom reset dates
		else
		{
			resetDatesId = ARM_NULL_OBJECT;
			cpnResetGap = cpnResetGapD;
		}

		/*double cpnResetGapD;
		double cpnResetGapDef=2;
		XL_readNumCellWD(XL_cpnResetGap,cpnResetGapD,cpnResetGapDef," ARM_ERR: coupon reset gap: numerical expected",C_result);
		long cpnResetGap = cpnResetGapD;*/

		long C_intRule;
		XL_GETINTRULEWD( XL_intRule, C_intRule, "ADJ",		" ARM_ERR: fwdRule : string expected",		C_result );

		/// RF coupon profiles
		/// ------------------
		/// RF coupon index leverage
		VECTOR<CCString> cpnCurves;
		XL_readStrVector(XL_cpnCurves,cpnCurves," ARM_ERR: cpn curves: vector of string or double expected",DOUBLE_TYPE,C_result);

		double leverage;
		long     leverageId;

		if(cpnCurves.size() >= 1)
		{
			if(cpnCurves[0][0] == 'L')
				leverageId = LocalGetNumObjectId(cpnCurves[0]);
			else
			{
				leverageId = ARM_NULL_OBJECT;
				leverage = atof(cpnCurves[0].c_str());
			}
		}
		else
		{
			CCString local_msg ("ARM_ERR: cpn curves [1] : leverage: string or double expected");
			C_result.setMsg (local_msg);
			ARM_ARG_ERR();
		}

		double cpnMin;
		long     cpnMinId;

		if(cpnCurves.size() >= 2)
		{
			if(cpnCurves[0][1] == 'L')
				cpnMinId = LocalGetNumObjectId(cpnCurves[1]);
			else
			{
				cpnMinId = ARM_NULL_OBJECT;
				cpnMin = atof(cpnCurves[1].c_str());
			}
		}
		else
		{
			cpnMinId = ARM_NULL_OBJECT;
			cpnMin = 0.0;
		}

		double cpnMax;
		long     cpnMaxId;

		if(cpnCurves.size() >= 3)
		{
			if(cpnCurves[0][2] == 'L')
				cpnMaxId = LocalGetNumObjectId(cpnCurves[2]);
			else
			{
				cpnMaxId = ARM_NULL_OBJECT;
				cpnMax = atof(cpnCurves[2].c_str());
			}
		}
		else
		{
			cpnMaxId = ARM_NULL_OBJECT;
			cpnMax = 1000.0;
		}

		/// TARN life time cap
		/// TARN life time cap
		VECTOR<CCString> lifeTimeCapDatas;
		XL_readStrVector(XL_lifeTimeCapDatas,lifeTimeCapDatas," ARM_ERR: lifeTime cap data: string expected",DOUBLE_TYPE,C_result);
		double lifeTimeCapTarget;
		bool globalCapFlag = true;

		if(lifeTimeCapDatas.size() >= 1)
		{
			lifeTimeCapTarget = atof(lifeTimeCapDatas[0]);
		}
		else
		{
			CCString local_msg ("ARM_ERR: life time cap datas [1] : LifeTimeCapTarget: double expected");
			C_result.setMsg (local_msg);
			ARM_ARG_ERR();
		}

		if(lifeTimeCapDatas.size() >= 2)
		{
			globalCapFlag = ((lifeTimeCapDatas[1]=="Y")||(lifeTimeCapDatas[1]=="YES"));
		}
		else
		{
			globalCapFlag = true;
		}

		/// TARN life time Floor
		double lifeTimeFloorTarget;
		double lifeTimeFloorTargetDef=0.0;
		XL_readNumCellWD(XL_lifeTimeFloorTarget,lifeTimeFloorTarget,lifeTimeFloorTargetDef, " ARM_ERR: lifeTime Floor: numerical or refValue Id expected",C_result);

		/// Funding leg datas
		/// -----------------
		/// Funding coupon spread profile
		double fundSpread;
		double fundSpreadDef=0.0;
		CCString fundSpreadXl;
		long     fundSpreadId;
		XL_readStrOrNumCellWD(XL_fundSpread,fundSpreadXl,fundSpread,fundSpreadDef,fundSpreadId,
			" ARM_ERR: funding spread: numerical or refValue Id expected",C_result);
		if(fundSpreadId == XL_TYPE_STRING)
			fundSpreadId = LocalGetNumObjectId(fundSpreadXl);
		else
			fundSpreadId = ARM_NULL_OBJECT;

		VECTOR<CCString> fundDatas;
		XL_readStrVector (XL_fundDatas,fundDatas," ARM_ERR: Funding datas: array of string expected",DOUBLE_TYPE,C_result);

		/// Funding coupon frequency (reset & payment)
		CCString fundFreqXl(cpnFreqXl);
		if(fundDatas.size() >= 1)
			fundFreqXl=fundDatas[0];
		long fundFreq;
		if((fundFreq = ARM_ConvFrequency (fundFreqXl, C_result)) == ARM_DEFAULT_ERR)
		{
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}

		/// Funding coupon day count
		CCString fundDayCountXl("A360");
		if(fundDatas.size() >= 2)
			fundDayCountXl=fundDatas[1];
		long fundDayCount = ARM_ConvDayCount (fundDayCountXl);


		/// Funding & RF legs notional
		/// --------------------------
		double nominal, fundNom, fees;
		double nominalDef=100.0, feesDef=1.0;
		CCString nominalStr;
		CCString fundNomStr;
		CCString feesStr;
		long     feesId;
		long     nominalId;
		long     fundNomId;

		XL_readStrOrNumCellWD(XL_nominal,nominalStr,nominal,nominalDef,nominalId,
			" ARM_ERR: nominal: numerical or refValue Id expected",C_result);
		if(nominalId == XL_TYPE_STRING)
			nominalId = LocalGetNumObjectId(nominalStr);
		else
			nominalId = ARM_NULL_OBJECT;

		XL_readStrOrNumCellWD(XL_fees,feesStr,fees,feesDef,feesId,
			" ARM_ERR: fees: numerical or refValue Id expected",C_result);
		if(feesId == XL_TYPE_STRING)
			feesId = LocalGetNumObjectId(feesStr);
		else
			feesId = ARM_NULL_OBJECT;

		XL_readStrOrNumCellWD(XL_fundNominal,fundNomStr,fundNom,nominalDef,fundNomId,
			" ARM_ERR: funding nominal: numerical or refValue Id expected",C_result);
		if(fundNomId == XL_TYPE_STRING)
			fundNomId = LocalGetNumObjectId(fundNomStr);
		else
			fundNomId = nominalId;

		/// Numerical Datas
		VECTOR<CCString> numDatas;
		VECTOR<CCString> numDatasDef(2,"10000.0");
		XL_readStrVectorWD(XL_numDatas,numDatas,numDatasDef," ARM_ERR: num Datas: vector of string expected",DOUBLE_TYPE,C_result);

		if (numDatas.size() > 16)
			throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"Num Data size is required to be less than 16");

		// Nb Iterations for MC
		VECTOR<double> nbIterations;
		if( numDatas.size() >= 1 )
			nbIterations.push_back( atof(numDatas[0]) );
		else
			nbIterations.push_back( 10000.0 );

		// Nb Iterations for Exercise Strike (calibration)
		if( numDatas.size() >= 2 )
			nbIterations.push_back( atof(numDatas[1]) );
		else
			nbIterations.push_back( 10000.0 );

		// Bucket Size
		if( numDatas.size() >= 3 )
			nbIterations.push_back( atof(numDatas[2]) );
		else
			nbIterations.push_back( 10000.0 );

		// Inter step
		if( numDatas.size() >= 4 )
			nbIterations.push_back( atof(numDatas[3]) );
		else
			nbIterations.push_back( 1.0 );

		// PDE_nbSteps
		if( numDatas.size() >= 5 )
			nbIterations.push_back( atof(numDatas[4]) );
		else
			nbIterations.push_back( 1000.0 );

		// PDE_gridSize
		if( numDatas.size() >= 6 )
			nbIterations.push_back( atof(numDatas[5]) );
		else
			nbIterations.push_back( 901.0 );

		// PDE_nbStdDev
		if( numDatas.size() >= 7 )
			nbIterations.push_back( atof(numDatas[6]) );
		else
			nbIterations.push_back( 6.0 );

		// Factor Max
		if( numDatas.size() >= 8 )
			nbIterations.push_back( atof(numDatas[7]) );
		else
			nbIterations.push_back( 2.0 );

		// SABR Grid Size
		if( numDatas.size() >= 9 )
			nbIterations.push_back( atof(numDatas[8]) );
		else
			nbIterations.push_back( 901 );


		/// default generator 1 = MRGK5
		string GenType1 = "MRGK5"; 

		if (numDatas.size() >= 10)
		{
			GenType1 = CCSTringToSTLString(numDatas[9]);
			GenType1 = stringGetUpper(GenType1);
		}

		/// default generator 2 = Sobol
		string GenType2 = "Sobol"; 

		if (numDatas.size() >= 11)
		{
			GenType2 = CCSTringToSTLString(numDatas[10]);
			GenType2 = stringGetUpper(GenType2);
		}

		/// default first nb times = 0
		double FirstNbTimes = 0; 

		if (numDatas.size() >= 12)
		{
			FirstNbTimes = atof(numDatas[11]);
		}

		/// default first nb dims = 0
		double FirstNbDims = 0; 

		if (numDatas.size() >= 13)
		{
			FirstNbDims = atof(numDatas[12]);
		}

		/// default path scheme = Incremental
		string PathScheme = "Incremental"; 

		if (numDatas.size() >= 14)
		{
			PathScheme = CCSTringToSTLString(numDatas[13]);
			PathScheme = stringGetUpper(PathScheme);
		}

		/// default path order = Bucket Order
		string PathOrder = "BucketOrder";

		if (numDatas.size() >= 15)
		{
			PathOrder = CCSTringToSTLString(numDatas[14]);
			PathOrder = stringGetUpper(PathOrder);
		}

		/// default antithetic = "Y"
		string Antithetic = "Y";

		if (numDatas.size() >= 16)
		{
			Antithetic = CCSTringToSTLString(numDatas[15]);
			Antithetic = stringGetUpper(Antithetic);
		}

		/// Calibration Flags datas
		/// -----------
		VECTOR<CCString> C_calibFlags;
		VECTOR<CCString> calibFlagsDef(6,"RF");
        calibFlagsDef[1]=CCString("N"); // Beta/Digital Calib
		calibFlagsDef[2]=CCString("N"); // Swaption/MRS calib
		calibFlagsDef[3]=CCString("N"); // Control Variable
		calibFlagsDef[4]=CCString("N"); // Digital Smoothing
		XL_readStrVectorWD(XL_calibFlags,C_calibFlags,calibFlagsDef," ARM_ERR: calib flags: array of string expected",DOUBLE_TYPE,C_result);

		vector< string > calibFlags(C_calibFlags.size());
        int i;
		for(i=0;i<C_calibFlags.size();++i)
        {
            (C_calibFlags[i]).toUpper();
			calibFlags[i]=CCSTringToSTLString(C_calibFlags[i]);
        }

		/// output Flags datas
		/// -----------
		VECTOR<CCString> C_outputFlags;
		VECTOR<CCString> outputFlagsDef(11,"N");
        outputFlagsDef[0] = CCString("Y"); // TARN Price

		XL_readStrVectorWD(XL_outputFlags,C_outputFlags,outputFlagsDef," ARM_ERR: output flags: array of string expected",DOUBLE_TYPE,C_result);

		vector< string > outputFlags(C_outputFlags.size());
		for(i = 0; i < C_outputFlags.size(); ++i)
        {
            (C_outputFlags[i]).toUpper();
			outputFlags[i]=CCSTringToSTLString(C_outputFlags[i]);
        }

        double C_asOf = -1.0;;
		double C_asOf_default = -1.0;
		
		XL_readNumCellWD(XL_asOfDate, C_asOf,
                         C_asOf_default," ARM_ERR: As Of Date: date expected", C_result);
		
        /// use the concept of Functor to transfer the knowledge of
		/// a function with a context

		tarnCalculatorFunc ourFunc(startDate,
				                   endDate,
				                   strike,
				                   strikeId,
				                   payRec,
				                   cpnDayCount,
				                   cpnFreq,
				                   cpnTiming,
				                   cpnIndexTerm,
				                   cpnIndexDayCount,
				                   cpnResetCal,
				                   cpnPayCal,
				                   cpnResetGap,
				                   C_intRule,
				                   leverage,
				                   leverageId,
								   cpnMin,
								   cpnMinId,
								   cpnMax,
								   cpnMaxId,	
				                   lifeTimeCapTarget,
								   globalCapFlag,
				                   lifeTimeFloorTarget,
				                   fundSpread,
				                   fundSpreadId,
				                   fundFreq,
				                   fundDayCount,
				                   nominal,
				                   nominalId,
				                   fees,
				                   feesId,
				                   fundNomId,
				                   nbIterations,
				GenType1,
				GenType2,
				FirstNbTimes,
				FirstNbDims,
				PathScheme,
				PathOrder,
				Antithetic,
                                   calibFlags,
				                   outputFlags,
				                   mktDataManagerId,
				                   keys,
				                   isCustomResetFlag,
				                   resetDatesId,
                                   C_asOf);

		/// call the general function
		fillXL_Result(LOCAL_GC_TARN_CLASS, ourFunc, C_result, XL_result, PersistentInXL );
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_TARNCalculator_Common" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}
							 

///////////////////////////////////
/// version that takes into account 
/// previous creation of object
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_TARNCalculator_Create(
    LPXLOPER XL_startDate,
    LPXLOPER XL_endDate,
    LPXLOPER XL_strike,
    LPXLOPER XL_payRec,
    LPXLOPER XL_cpnDatas,
    LPXLOPER XL_mktDatas,
    LPXLOPER XL_cpnResetGap,
	LPXLOPER XL_intRule,
    LPXLOPER XL_cpnCurves,
	LPXLOPER XL_lifeTimeCapDatas,
    LPXLOPER XL_lifeTimeFloorTarget,
    LPXLOPER XL_fundSpread,
    LPXLOPER XL_fundDatas,
    LPXLOPER XL_nominal,
	LPXLOPER XL_nbIter,
    LPXLOPER XL_calibFlags,
    LPXLOPER XL_outputFlags,
	LPXLOPER XL_fundNominal,
	LPXLOPER XL_fees,
    LPXLOPER XL_asOfDate)
{
	ADD_LOG("Local_TARNCalculator_Create");
	bool PersistentInXL = true;

	return Local_TARNCalculator_Common(
        XL_startDate,
        XL_endDate,
        XL_strike,
        XL_payRec,
        XL_cpnDatas,
        XL_mktDatas,
        XL_cpnResetGap,
		XL_intRule,
        XL_cpnCurves,
		XL_lifeTimeCapDatas,
		XL_lifeTimeFloorTarget,
        XL_fundSpread,
        XL_fundDatas,
        XL_nominal,
		XL_nbIter,
        XL_calibFlags,
		XL_outputFlags,
		XL_fundNominal,
		XL_fees,
        XL_asOfDate,
		PersistentInXL );
}

LPXLOPER Local_TarnSetOutputFlags_Common(
        LPXLOPER	XL_TarnId,
        LPXLOPER	XL_ProductFlags,
		bool		XL_PersistentInXL)

{
	// this is defined first because it is used in XL macros
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

		//Get the TARN Id
		long tarnId;
		long objId;
		long retCode;
		CCString C_ProductFlags;
		XL_GETOBJID(XL_TarnId,tarnId," ARM_ERR: TARN calculator id: object expected",C_result);
		
		//Get the OUTPUT flags
		VECTOR<CCString> C_outputFlags;
		VECTOR<CCString> outputFlagsDef(10,"N");
        outputFlagsDef[0]=CCString("Y"); // TARN Price
		XL_readStrVectorWD(XL_ProductFlags,C_outputFlags,outputFlagsDef," ARM_ERR: output flags: array of string expected",DOUBLE_TYPE,C_result);

		vector< string > outputFlags(C_outputFlags.size());
		for(int i=0;i<C_outputFlags.size();++i)
        {
            (C_outputFlags[i]).toUpper();
			outputFlags[i]=CCSTringToSTLString(C_outputFlags[i]);
        }
	
		//Is it a first cell calculation ? 
		CCString stringId = GetLastCurCellEnvValue ();
		//First cell calculation: we create a new Tarn object.
		if(!stringId)
		{
			const long CreateBrandNew = ARM_NULL_OBJECT_ID;
			retCode = ARMLOCAL_TARN_SetProductToPrice(tarnId, outputFlags, C_result, CreateBrandNew );
			CCString className = "LTARN"; //C_result.getShortName() ;

			if (retCode == ARM_OK)
			{
				objId = C_result.getLong();
				LocalSetCurCellEnvValue (className, objId); 
				stringId = LocalMakeObjectId (objId, className);
			}
		}
		//Not a first cell calculation: we replace the existing object.
		else
		{
			CCString prevClass = LocalGetStringObjectClass (stringId);
			objId	= LocalGetNumObjectId (stringId);

			retCode = ARMLOCAL_TARN_SetProductToPrice(tarnId, outputFlags, C_result, objId);
			CCString className = "LTARN";

			//Same type as the previous object.
			if( className == prevClass )
			{
				if( retCode == ARM_OK )
				{ 
					LocalSetCurCellEnvValue (className, objId);
					stringId = LocalMakeObjectId (objId, className);
				}
			}
			else 
			{
				FreeCurCellContent ();
				if( retCode == ARM_OK )
				{
					objId = C_result.getLong ();
					LocalSetCurCellEnvValue (className, objId);
					stringId = LocalMakeObjectId (objId, className);
				}
			}
		}

		/// feed the LPXLOPER object result 
		if (retCode == ARM_OK)
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_TarnCalculator_SetProductToPrice" )

	return (LPXLOPER)&XL_result;
}



//Product: TARN / Aim: Set the pricing outputs
__declspec(dllexport) LPXLOPER WINAPI Local_TarnSetOutputFlags(
    LPXLOPER XL_TarnId,
    LPXLOPER XL_ProductFlags)
{
	ADD_LOG("Local_TarnSetOutputFlags");
	bool PersistentInXL = true;
	return Local_TarnSetOutputFlags_Common(
        XL_TarnId,
        XL_ProductFlags,
		PersistentInXL);
}



///////////////////////////////////
/// version for non persistentInXl
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_PXL_TARNCalculator_Create(
    LPXLOPER XL_startDate,
    LPXLOPER XL_endDate,
    LPXLOPER XL_strike,
    LPXLOPER XL_payRec,
    LPXLOPER XL_cpnDatas,
    LPXLOPER XL_mktDatas,
    LPXLOPER XL_cpnResetGap,
	LPXLOPER XL_intRule,
    LPXLOPER XL_cpnCurves,
	LPXLOPER XL_lifeTimeCapTarget,
    LPXLOPER XL_lifeTimeFloorTarget,
    LPXLOPER XL_fundSpread,
    LPXLOPER XL_fundDatas,
    LPXLOPER XL_nominal,
	LPXLOPER XL_nbIter,
    LPXLOPER XL_calibFlags,
    LPXLOPER XL_outputFlags,
	LPXLOPER XL_fundNominal,
	LPXLOPER XL_fees,
    LPXLOPER XL_asOfDate)
{
	ADD_LOG("Local_PXL_TARNCalculator_Create");
	bool PersistentInXL = false;
	return Local_TARNCalculator_Common(
        XL_startDate,
        XL_endDate,
        XL_strike,
        XL_payRec,
        XL_cpnDatas,
        XL_mktDatas,
        XL_cpnResetGap,
		XL_intRule,
        XL_cpnCurves,
		XL_lifeTimeCapTarget,
		XL_lifeTimeFloorTarget,
        XL_fundSpread,
        XL_fundDatas,
        XL_nominal,
		XL_nbIter,
        XL_calibFlags,
        XL_outputFlags,
		XL_fundNominal,
		XL_fees,
        XL_asOfDate,
		PersistentInXL);
}



///----------------------------------------------
///----------------------------------------------
///             TARN SnowBallCalculator
/// Inputs :
///     
///----------------------------------------------
///----------------------------------------------
class tarnSBCalculatorFunc : public ARMResultLong2LongFunc
{
public:
	tarnSBCalculatorFunc(
        double startDate,
        double endDate,
        double strike,
        long strikeId,
        double coupon0,
        long payRec,
        long cpnDayCount,
        long cpnFreq,
        long cpnTiming,
        const string& cpnIndexTerm,
        long cpnIndexDayCount,
        const string& cpnResetCal,
        const string& cpnPayCal,
        long cpnResetGap,
		long intRule,
		string stubRuleStr,
        double leverage,
        long leverageId,
		double cpnMin,
		long cpnMinId,
		double cpnMax,
		long cpnMaxId,
		double levPrev,
        long levPrevId,
		long reverse,
		double lifeTimeCapTarget,
		bool globalCapFlag,
		double lifeTimeFloorTarget,
        double fundSpread,
        long fundSpreadId,
        long  fundFreq,
        long  fundDayCount,
        double nominal,
        long nominalId,
		double fees,
		long feesId,
		double fundNominal,
		long fundNominalId,
		const vector<double>& nbIterations,
		const string& genType1,
		const string& genType2,
		int firstNbTimes,
		int firstNbDims,
		const string& pathScheme,
		const string& pathOrder,
		const string& antithetic,
        const vector< string >& calibFlags,
        const vector< string >& outputFlags,
        long mktDataManagerId,
        const vector< string >& keys)
    :
    C_startDate(startDate),
    C_endDate(endDate),
    C_strike(strike),
    C_strikeId(strikeId),
    C_coupon0(coupon0),
    C_payRec(payRec),
    C_cpnDayCount(cpnDayCount),
    C_cpnFreq(cpnFreq),
    C_cpnTiming(cpnTiming),
    C_cpnIndexTerm(cpnIndexTerm),
    C_cpnIndexDayCount(cpnIndexDayCount),
    C_cpnResetCal(cpnResetCal),
    C_cpnPayCal(cpnPayCal),
    C_cpnResetGap(cpnResetGap),
	C_intRule(intRule),
	C_intStubRule(stubRuleStr),
    C_leverage(leverage),
    C_leverageId(leverageId),
	C_cpnMin(cpnMin),
	C_cpnMinId(cpnMinId),
	C_cpnMax(cpnMax),
	C_cpnMaxId(cpnMaxId),
	C_levPrev(levPrev),
    C_levPrevId(levPrevId),
	C_reverse(reverse),
	C_lifeTimeCapTarget(lifeTimeCapTarget),
	C_globalCapFlag(globalCapFlag),
	C_lifeTimeFloorTarget(lifeTimeFloorTarget),
    C_fundSpread(fundSpread),
    C_fundSpreadId(fundSpreadId),
    C_fundFreq(fundFreq),
    C_fundDayCount(fundDayCount),
    C_nominal(nominal),
    C_nominalId(nominalId),
	C_fees(fees),
	C_feesId(feesId),
	C_fundNominal(fundNominal),
	C_fundNominalId(fundNominalId),
	C_nbIterations(nbIterations),
	C_genType1(genType1),
	C_genType2(genType2),
	C_firstNbTimes(firstNbTimes),
	C_firstNbDims(firstNbDims),
	C_pathScheme(pathScheme),
	C_pathOrder(pathOrder),
	C_antithetic(antithetic),
    C_calibFlags(calibFlags),
    C_outputFlags(outputFlags),
    C_mktDataManagerId(mktDataManagerId),
    C_keys(keys)
    {};
	
	long operator()( ARM_result& result, long objId ){
		return ARMLOCAL_TARNSBCalculator_Create(
            C_startDate,
            C_endDate,
            C_strike,
            C_strikeId,
            C_coupon0,
            C_payRec,
            C_cpnDayCount,
            C_cpnFreq,
            C_cpnTiming,
            C_cpnIndexTerm,
            C_cpnIndexDayCount,
            C_cpnResetCal,
            C_cpnPayCal,
            C_cpnResetGap,
			C_intRule,
			C_intStubRule,
            C_leverage,
            C_leverageId,
			C_cpnMin,
			C_cpnMinId,
			C_cpnMax,
			C_cpnMaxId,
			C_levPrev,
            C_levPrevId,
			C_reverse,
			C_lifeTimeCapTarget,
			C_globalCapFlag,
			C_lifeTimeFloorTarget,
            C_fundSpread,
            C_fundSpreadId,
            C_fundFreq,
            C_fundDayCount,
            C_nominal,
            C_nominalId,
			C_fees,
			C_feesId,
			C_fundNominal,
			C_fundNominalId,
			C_nbIterations,
			C_genType1,
			C_genType2,
			C_firstNbTimes,
			C_firstNbDims,
			C_pathScheme,
			C_pathOrder,
			C_antithetic,
            C_calibFlags,
            C_outputFlags,
            C_mktDataManagerId,
            C_keys,
            result,
            objId);
    }

private:
    double              C_startDate;
    double              C_endDate;
    double              C_strike;
    long                C_strikeId;
	double				C_coupon0;
    long                C_payRec;
    long                C_cpnDayCount;
	long				C_cpnFreq;
    long                C_fixDayCount;
    long                C_cpnTiming;
    string              C_cpnIndexTerm;
    long                C_cpnIndexDayCount;
    string              C_cpnResetCal;
    string              C_cpnPayCal;
    long                C_cpnResetGap;
	long				C_intRule;
	string				C_intStubRule;
    double              C_leverage;
    long                C_leverageId;
	double				C_cpnMin;
	long				C_cpnMinId;
	double				C_cpnMax;
	long				C_cpnMaxId;
	double              C_levPrev;
    long                C_levPrevId;
	long				C_reverse;
	double				C_lifeTimeCapTarget;
	bool				C_globalCapFlag;
	double				C_lifeTimeFloorTarget;
    double              C_fundSpread;
    long                C_fundSpreadId;
    long                C_fundFreq;
    long                C_fundDayCount;
    double              C_nominal;
    long                C_nominalId;
	double				C_fundNominal;
	long				C_fundNominalId;
	double				C_fees;
	long				C_feesId;
	vector<double>		C_nbIterations;
	string				C_genType1;
	string				C_genType2;
	int					C_firstNbTimes;
	int					C_firstNbDims;
	string				C_pathScheme;
	string				C_pathOrder;
	string				C_antithetic;
    vector< string >    C_calibFlags;
    vector< string >    C_outputFlags;
    long                C_mktDataManagerId;
    vector< string >    C_keys;
};


/////////////////////////////////////////////////////////////
/// central function that does the creation of the XL function
/////////////////////////////////////////////////////////////
LPXLOPER Local_TARNSBCalculator_Common(
    LPXLOPER XL_startDate,
    LPXLOPER XL_endDate,
    LPXLOPER XL_strike,
	LPXLOPER XL_coupon0,
    LPXLOPER XL_payRec,
    LPXLOPER XL_cpnDatas,
    LPXLOPER XL_mktDatas,
    LPXLOPER XL_cpnResetGap,
	LPXLOPER XL_intRule,
    LPXLOPER XL_cpnCurves,
    LPXLOPER XL_levPrev,
	LPXLOPER XL_lifeTimeCapDatas,
    LPXLOPER XL_lifeTimeFloorTarget,
    LPXLOPER XL_fundSpread,
    LPXLOPER XL_fundDatas,
    LPXLOPER XL_nominals,
	LPXLOPER XL_numDatas,
    LPXLOPER XL_calibFlags,
    LPXLOPER XL_outputFlags,
    LPXLOPER XL_fundingCcyName,
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

    
		/// General TARN datas
		/// -----------------
		/// Start & End
		double startDate;
		XL_readNumCell(XL_startDate,startDate," ARM_ERR: start date: date expected",C_result);
		double endDate;
		XL_readNumCell(XL_endDate,endDate," ARM_ERR: end date: date expected",C_result);

		/// Fix rate or RF Strike
		double strike;
		CCString strikeStr;
		long     strikeId;
		XL_readStrOrNumCell(XL_strike, strikeStr, strike, strikeId,
			   " ARM_ERR: strike: numeric or refValue Id expected",C_result);	
		if(strikeId == XL_TYPE_STRING)
			strikeId = LocalGetNumObjectId(strikeStr);
		else
			strikeId = ARM_NULL_OBJECT;

		double coupon0;
		XL_readNumCell(XL_coupon0,coupon0," ARM_ERR: coupon 0: numeric expected",C_result);

		/// RF Leg Pay/Rec
		CCString payRecStr;
		long payRec;
		XL_readStrCell(XL_payRec,payRecStr," ARM_ERR: payer/receiver: string expected",C_result);
		if((payRec = ARM_ConvRecOrPay (payRecStr, C_result)) == ARM_DEFAULT_ERR)
		{
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}


		/// Required RF coupon datas
		/// ------------------------
		VECTOR<CCString> cpnDatas;
		XL_readStrVector (XL_cpnDatas,cpnDatas," ARM_ERR: Coupon datas: array of string expected",DOUBLE_TYPE,C_result);
		if(cpnDatas.size() < 2)
		{
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}

		/// RF coupon frequency (reset & payment)
		CCString cpnFreqXl=cpnDatas[0];
		long cpnFreq;
		if((cpnFreq = ARM_ConvFrequency (cpnFreqXl, C_result)) == ARM_DEFAULT_ERR)
		{
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}

		/// RF coupon day count
		CCString cpnDayCountXl=cpnDatas[1];
		long cpnDayCount = ARM_ConvDayCount (cpnDayCountXl);


		/// Market datas : curve , MRS and market models for OSW & CF pricings
		/// ------------------------------------------------------------------
		VECTOR<CCString> mktDatas;
		XL_readStrVector(XL_mktDatas,mktDatas," ARM_ERR: Market datas: array of string expected",DOUBLE_TYPE,C_result);
        if(mktDatas.size()<1)
		{
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}

		long mktDataManagerId = LocalGetNumObjectId (mktDatas[0]);
		vector< string > keys(mktDatas.size()-1);
        size_t i;
		for(i=1;i<mktDatas.size();++i)
			keys[i-1]=CCSTringToSTLString(mktDatas[i]);



		//// -------------------------------
		///  ----- OPTIONAL ARGUMENTS ------
		//// -------------------------------



		/// RF coupon datas
		/// ---------------
		/// Advance/Arrears RF coupon reset
		CCString cpnTimingXl("ADV");
		if(cpnDatas.size() >= 3)
			cpnTimingXl=cpnDatas[2];
		long cpnTiming = ARM_ConvPayResetRule (cpnTimingXl);

		/// RF index term
		CCString cpnIndexTermXl(ARM_ArgConvReverse_MatFrequency.GetString(cpnFreq).c_str());
		if(cpnDatas.size() >= 4)
			cpnIndexTermXl=cpnDatas[3];
		string cpnIndexTerm=CCSTringToSTLString(cpnIndexTermXl);

		/// RF index day count
		CCString cpnIndexDayCountXl("A360");
		if(cpnDatas.size() >= 5)
			cpnIndexDayCountXl=cpnDatas[4];
		long cpnIndexDayCount = ARM_ConvDayCount (cpnIndexDayCountXl);

		/// RF reset calendar
		CCString cpnResetCalXl("EUR");
		if(cpnDatas.size() >= 6)
			cpnResetCalXl=cpnDatas[5];
		string cpnResetCal = CCSTringToSTLString(cpnResetCalXl);

		/// RF payment calendar
		CCString cpnPayCalXl("EUR");
		if(cpnDatas.size() >= 7)
			cpnPayCalXl=cpnDatas[6];
		string cpnPayCal = CCSTringToSTLString(cpnPayCalXl);

		/// Reverse Yes/No
		CCString cpnReverseXl("Y");
		long reverse;
		if(cpnDatas.size() >= 8)
			cpnReverseXl= cpnDatas[7];
		if((reverse = ARM_ConvYesOrNo (cpnReverseXl, C_result)) == ARM_DEFAULT_ERR)
		{
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}

		/// TARN StubRule
		CCString strubRuleXl("SS");
		if(cpnDatas.size() >= 9)
			strubRuleXl=cpnDatas[8];
		string stubRuleStr = CCSTringToSTLString(strubRuleXl);

		/// RF reset gap
		double cpnResetGapD;
		double cpnResetGapDef=2;
		XL_readNumCellWD(XL_cpnResetGap,cpnResetGapD,cpnResetGapDef," ARM_ERR: coupon reset gap: numerical expected",C_result);
		long cpnResetGap = cpnResetGapD;

		long C_intRule;
		XL_GETINTRULEWD( XL_intRule, C_intRule, "ADJ",		" ARM_ERR: fwdRule : string expected",		C_result );


		/// RF coupon profiles
		/// ------------------
		/// RF coupon index leverage

		VECTOR<CCString> cpnCurves;
		XL_readStrVector(XL_cpnCurves,cpnCurves," ARM_ERR: cpn curves: vector of string or double expected",DOUBLE_TYPE,C_result);

		double leverage;
		long     leverageId;

		if(cpnCurves.size() >= 1)
		{
			if(cpnCurves[0][0] == 'L')
				leverageId = LocalGetNumObjectId(cpnCurves[0]);
			else
			{
				leverageId = ARM_NULL_OBJECT;
				leverage = atof(cpnCurves[0].c_str());
			}
		}
		else
		{
			CCString local_msg ("ARM_ERR: cpn curves [1] : leverage: string or double expected");
			C_result.setMsg (local_msg);
			ARM_ARG_ERR();
		}

		double cpnMin;
		long     cpnMinId;

		if(cpnCurves.size() >= 2)
		{
			if(cpnCurves[1][0] == 'L')
				cpnMinId = LocalGetNumObjectId(cpnCurves[1]);
			else
			{
				cpnMinId = ARM_NULL_OBJECT;
				cpnMin = atof(cpnCurves[1].c_str());
			}
		}
		else
		{
			cpnMinId = ARM_NULL_OBJECT;
			cpnMin = 0.0;
		}

		double cpnMax;
		long     cpnMaxId;

		if(cpnCurves.size() >= 3)
		{
			if(cpnCurves[2][0] == 'L')
				cpnMaxId = LocalGetNumObjectId(cpnCurves[2]);
			else
			{
				cpnMaxId = ARM_NULL_OBJECT;
				cpnMax = atof(cpnCurves[2].c_str());
			}
		}
		else
		{
			cpnMaxId = ARM_NULL_OBJECT;
			cpnMax = 1000.0;
		}

		double levPrev;
		double levPrevDef=1.0;
		CCString levPrevStr;
		long     levPrevId;
		XL_readStrOrNumCellWD(XL_levPrev,levPrevStr,levPrev,levPrevDef,levPrevId,
			" ARM_ERR: levPrev: numeric or refValue Id expected",C_result);
		if(levPrevId == XL_TYPE_STRING)
			levPrevId = LocalGetNumObjectId(levPrevStr);
		else
			levPrevId = ARM_NULL_OBJECT;

		/// TARN life time cap
		VECTOR<CCString> lifeTimeCapDatas;
		XL_readStrVector(XL_lifeTimeCapDatas,lifeTimeCapDatas," ARM_ERR: lifeTime cap data: string expected",DOUBLE_TYPE,C_result);
		double lifeTimeCapTarget;
		bool globalCapFlag = true;

		if(lifeTimeCapDatas.size() >= 1)
		{
			lifeTimeCapTarget = atof(lifeTimeCapDatas[0]);
		}
		else
		{
			CCString local_msg ("ARM_ERR: life time cap datas [1] : LifeTimeCapTarget: double expected");
			C_result.setMsg (local_msg);
			ARM_ARG_ERR();
		}

		if(lifeTimeCapDatas.size() >= 2)
		{
			globalCapFlag = ((lifeTimeCapDatas[1]=="Y")||(lifeTimeCapDatas[1]=="YES"));
		}
		else
		{
			globalCapFlag = true;
		}

		/// TARN life time Floor
		double lifeTimeFloorTarget;
		double lifeTimeFloorTargetDef=0.0;
		XL_readNumCellWD(XL_lifeTimeFloorTarget,lifeTimeFloorTarget,lifeTimeFloorTargetDef, " ARM_ERR: lifeTime Floor: numerical or refValue Id expected",C_result);

		/// Funding leg datas
		/// -----------------
		/// Funding coupon spread profile
		double fundSpread;
		double fundSpreadDef=0.0;
		CCString fundSpreadXl;
		long     fundSpreadId;
		XL_readStrOrNumCellWD(XL_fundSpread,fundSpreadXl,fundSpread,fundSpreadDef,fundSpreadId,
			" ARM_ERR: funding spread: numerical or refValue Id expected",C_result);
		if(fundSpreadId == XL_TYPE_STRING)
			fundSpreadId = LocalGetNumObjectId(fundSpreadXl);
		else
			fundSpreadId = ARM_NULL_OBJECT;

		VECTOR<CCString> fundDatas;
		XL_readStrVector (XL_fundDatas,fundDatas," ARM_ERR: Funding datas: array of string expected",DOUBLE_TYPE,C_result);

		/// Funding coupon frequency (reset & payment)
		CCString fundFreqXl(cpnFreqXl);
		if(fundDatas.size() >= 1)
			fundFreqXl=fundDatas[0];
		long fundFreq;
		if((fundFreq = ARM_ConvFrequency (fundFreqXl, C_result)) == ARM_DEFAULT_ERR)
		{
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}

		/// Funding coupon day count
		CCString fundDayCountXl("A360");
		if(fundDatas.size() >= 2)
			fundDayCountXl=fundDatas[1];
		long fundDayCount = ARM_ConvDayCount (fundDayCountXl);


		/// Funding & RF legs notional
		/// --------------------------
		VECTOR<CCString> nominals;
		XL_readStrVector(XL_nominals,nominals, " ARM_ERR: Nominals: array of string or numeric expected",DOUBLE_TYPE,C_result);

		long     nominalId;
		double nominal;

		CCString firstCc("L");
		if(nominals.size() >= 1)
		{
			if(nominals[0].Contain(firstCc))
				nominalId = LocalGetNumObjectId(nominals[0]);
			else
			{
				nominalId = ARM_NULL_OBJECT;
				nominal = atof(nominals[0].c_str());
			}
		}
		else
		{
			nominalId = ARM_NULL_OBJECT;
			nominal = 100.0;
		}

		long     fundNominalId;
		double fundNominal;

		if(nominals.size() >= 2)
		{
			if(nominals[1].Contain(firstCc))
				fundNominalId = LocalGetNumObjectId(nominals[1]);
			else
			{
				fundNominalId = ARM_NULL_OBJECT;
				fundNominal = atof(nominals[1].c_str());
			}
		}
		else
		{
			fundNominalId = nominalId;
			fundNominal = nominal;
		}

		double fees=0.0;
		double feesDef=1.0;
		CCString feesStr;
		long feesId = ARM_NULL_OBJECT;
		
		/// Numerical Datas
		VECTOR<CCString> numDatas;
		VECTOR<CCString> numDatasDef(2,"10000.0");
		XL_readStrVectorWD(XL_numDatas,numDatas,numDatasDef," ARM_ERR: num Datas: vector of string expected",DOUBLE_TYPE,C_result);

		if (numDatas.size() > 16)
			throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"Num Data size is required to be less than 16");

		// Nb Iterations for MC
		VECTOR<double> nbIterations;
		if( numDatas.size() >= 1 )
			nbIterations.push_back( atof(numDatas[0]) );
		else
			nbIterations.push_back( 10000.0 );

		// Nb Iterations for Exercise Strike (calibration)
		if( numDatas.size() >= 2 )
			nbIterations.push_back( atof(numDatas[1]) );
		else
			nbIterations.push_back( 10000.0 );

		// Bucket Size
		if( numDatas.size() >= 3 )
			nbIterations.push_back( atof(numDatas[2]) );
		else
			nbIterations.push_back( 10000.0 );

		// Inter step
		if( numDatas.size() >= 4 )
			nbIterations.push_back( atof(numDatas[3]) );
		else
			nbIterations.push_back( 1.0 );

		// PDE_nbSteps
		if( numDatas.size() >= 5 )
			nbIterations.push_back( atof(numDatas[4]) );
		else
			nbIterations.push_back( 1000.0 );

		// PDE_gridSize
		if( numDatas.size() >= 6 )
			nbIterations.push_back( atof(numDatas[5]) );
		else
			nbIterations.push_back( 901.0 );

		// PDE_nbStdDev
		if( numDatas.size() >= 7 )
			nbIterations.push_back( atof(numDatas[6]) );
		else
			nbIterations.push_back( 6.0 );

		// Factor Max
		if( numDatas.size() >= 8 )
			nbIterations.push_back( atof(numDatas[7]) );
		else
			nbIterations.push_back( 2.0 );

		// SABR Grid Size
		if( numDatas.size() >= 9 )
			nbIterations.push_back( atof(numDatas[8]) );
		else
			nbIterations.push_back( 251.0 );


		/// default generator 1 = MRGK5
		string GenType1 = "MRGK5"; 

		if (numDatas.size() >= 10)
		{
			GenType1 = CCSTringToSTLString(numDatas[9]);
			GenType1 = stringGetUpper(GenType1);
		}

		/// default generator 2 = Sobol
		string GenType2 = "Sobol"; 

		if (numDatas.size() >= 11)
		{
			GenType2 = CCSTringToSTLString(numDatas[10]);
			GenType2 = stringGetUpper(GenType2);
		}

		/// default first nb times = 0
		double FirstNbTimes = 0; 

		if (numDatas.size() >= 12)
		{
			FirstNbTimes = atof(numDatas[11]);
		}

		/// default first nb dims = 0
		double FirstNbDims = 0; 

		if (numDatas.size() >= 13)
		{
			FirstNbDims = atof(numDatas[12]);
		}

		/// default path scheme = Incremental
		string PathScheme = "Incremental"; 

		if (numDatas.size() >= 14)
		{
			PathScheme = CCSTringToSTLString(numDatas[13]);
			PathScheme = stringGetUpper(PathScheme);
		}

		/// default path order = Bucket Order
		string PathOrder = "BucketOrder";

		if (numDatas.size() >= 15)
		{
			PathOrder = CCSTringToSTLString(numDatas[14]);
			PathOrder = stringGetUpper(PathOrder);
		}

		/// default antithetic = "Y"
		string Antithetic = "Y";

		if (numDatas.size() >= 16)
		{
			Antithetic = CCSTringToSTLString(numDatas[15]);
			Antithetic = stringGetUpper(Antithetic);
		}
		

		/// calibration Flags datas
		/// -----------
		VECTOR<CCString> C_calibFlags;
		VECTOR<CCString> calibFlagsDef(5,"RF");
        calibFlagsDef[1]=CCString("N"); // Beta/Digital Calib
		calibFlagsDef[2]=CCString("N"); // Swaption/MRS calib
		calibFlagsDef[3]=CCString("N"); // Control Variable
		calibFlagsDef[4]=CCString("N"); // Digital Smoothing
		XL_readStrVectorWD(XL_calibFlags,C_calibFlags,calibFlagsDef," ARM_ERR: calibration flags: array of string expected",DOUBLE_TYPE,C_result);

		vector< string > calibFlags(C_calibFlags.size());
		for(i=0;i<C_calibFlags.size();++i)
        {
            (C_calibFlags[i]).toUpper();
			calibFlags[i]=CCSTringToSTLString(C_calibFlags[i]);
        }

		/// output Flags datas
		/// -----------
		VECTOR<CCString> C_outputFlags;
		VECTOR<CCString> outputFlagsDef(10,"N");
		outputFlagsDef[0] = CCString("Y"); // TARN Price

		XL_readStrVectorWD(XL_outputFlags,C_outputFlags,outputFlagsDef," ARM_ERR: output flags: array of string expected",DOUBLE_TYPE,C_result);

		vector< string > outputFlags(C_outputFlags.size());
		for(i=0;i<C_outputFlags.size();++i)
        {
            (C_outputFlags[i]).toUpper();
			outputFlags[i]=CCSTringToSTLString(C_outputFlags[i]);
        }

		/// use the concept of Functor to transfer the knowledge of
		/// a function with a context
		tarnSBCalculatorFunc ourFunc(
				startDate,
				endDate,
				strike,
				strikeId,
				coupon0,
				payRec,
				cpnDayCount,
				cpnFreq,
				cpnTiming,
				cpnIndexTerm,
				cpnIndexDayCount,
				cpnResetCal,
				cpnPayCal,
				cpnResetGap,
				C_intRule,
				stubRuleStr,
				leverage,
				leverageId,
				cpnMin,
				cpnMinId,
				cpnMax,
				cpnMaxId,
				levPrev,
				levPrevId,
				reverse,
				lifeTimeCapTarget,
				globalCapFlag,
				lifeTimeFloorTarget,
				fundSpread,
				fundSpreadId,
				fundFreq,
				fundDayCount,
				nominal,
				nominalId,
				fees,
				feesId,
				fundNominal,
				fundNominalId,
				nbIterations,
				GenType1,
				GenType2,
				FirstNbTimes,
				FirstNbDims,
				PathScheme,
				PathOrder,
				Antithetic,
                calibFlags,
				outputFlags,
				mktDataManagerId,
				keys);

		/// call the general function
		fillXL_Result( LOCAL_GC_TARN_CLASS, ourFunc, C_result, XL_result, PersistentInXL );
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_TARNCalculator_Common" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;

}
							 

///////////////////////////////////
/// version that takes into account 
/// previous creation of object
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_TARNSnowBallCalculator_Create(
    LPXLOPER XL_startDate,
    LPXLOPER XL_endDate,
    LPXLOPER XL_strike,
    LPXLOPER XL_coupon0,
    LPXLOPER XL_payRec,
    LPXLOPER XL_cpnDatas,
    LPXLOPER XL_mktDatas,
    LPXLOPER XL_cpnResetGap,
	LPXLOPER XL_intRule,
    LPXLOPER XL_cpnCurves,
	LPXLOPER XL_levPrev,
	LPXLOPER XL_lifeTimeCapDatas,
    LPXLOPER XL_lifeTimeFloorTarget,
    LPXLOPER XL_fundSpread,
    LPXLOPER XL_fundDatas,
    LPXLOPER XL_nominal,
	LPXLOPER XL_nbIter,
    LPXLOPER XL_calibFlags,
    LPXLOPER XL_outputFlags,
    LPXLOPER XL_fundingCcyName)
{
	ADD_LOG("Local_TARNSnowBallCalculator_Create");
	bool PersistentInXL = true;
	return Local_TARNSBCalculator_Common(
        XL_startDate,
        XL_endDate,
        XL_strike,
		XL_coupon0,
        XL_payRec,
        XL_cpnDatas,
        XL_mktDatas,
        XL_cpnResetGap,
		XL_intRule,
        XL_cpnCurves,
		XL_levPrev,
		XL_lifeTimeCapDatas,
		XL_lifeTimeFloorTarget,
        XL_fundSpread,
        XL_fundDatas,
        XL_nominal,
		XL_nbIter,
        XL_calibFlags,
        XL_outputFlags,
        XL_fundingCcyName,
		PersistentInXL );
}


///////////////////////////////////
/// version for non persistentInXl
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_PXL_TARNSnowBallCalculator_Create(
    LPXLOPER XL_startDate,
    LPXLOPER XL_endDate,
    LPXLOPER XL_strike,
    LPXLOPER XL_coupon0,
    LPXLOPER XL_payRec,
    LPXLOPER XL_cpnDatas,
    LPXLOPER XL_mktDatas,
    LPXLOPER XL_cpnResetGap,
	LPXLOPER XL_intRule,
    LPXLOPER XL_cpnCurves,
	LPXLOPER XL_levPrev,
	LPXLOPER XL_lifeTimeCapDatas,
    LPXLOPER XL_lifeTimeFloorTarget,
    LPXLOPER XL_fundSpread,
    LPXLOPER XL_fundDatas,
    LPXLOPER XL_nominal,
	LPXLOPER XL_nbIter,
    LPXLOPER XL_calibFlags,
    LPXLOPER XL_outputFlags,
    LPXLOPER XL_fundingCcyName)
{
	ADD_LOG("Local_PXL_TARNSnowBallCalculator_Create");
	bool PersistentInXL = false;
	return Local_TARNSBCalculator_Common(
        XL_startDate,
        XL_endDate,
        XL_strike,
		XL_coupon0,
        XL_payRec,
        XL_cpnDatas,
        XL_mktDatas,
        XL_cpnResetGap,
		XL_intRule,
        XL_cpnCurves,
		XL_levPrev,
		XL_lifeTimeCapDatas,
		XL_lifeTimeFloorTarget,
        XL_fundSpread,
        XL_fundDatas,
        XL_nominal,
		XL_nbIter,
        XL_calibFlags,
        XL_outputFlags,
        XL_fundingCcyName,
		PersistentInXL );
}


///--------------------------------------------------------
///--------------------------------------------------------
///             TARN (or TARN SnowBall) Calculator Accessor
/// Inputs :
///     
///--------------------------------------------------------
///--------------------------------------------------------
class tarnGetDataFunc : public ARMResultLong2LongFunc
{
public:
	tarnGetDataFunc(
        long tarnId,
        const string& getType)
    :
    C_tarnId(tarnId),
    C_getType(getType)
    {};
	
	long operator()( ARM_result& result, long objId ){
		return ARMLOCAL_TARN_Get(
            C_tarnId,
            C_getType,
            result,
            objId);
    }

private:
	long    C_tarnId;
	string  C_getType;
};

class tarnSetDataFunc : public ARMResultLong2LongFunc
{
public:
	tarnSetDataFunc(
        long tarnId,
        long dataToSetId,
        const string& setPortfolioType,
		const vector<string>& mktDataKeys,
        bool isUpdated)
    :
    C_tarnId(tarnId),
    C_dataToSetId(dataToSetId),
    C_setPortfolioType(setPortfolioType),
	C_mktDataKeys(mktDataKeys),
    C_isUpdated(isUpdated)
    {};
	
	long operator()( ARM_result& result, long objId ){
		return ARMLOCAL_TARN_Set(
            C_tarnId,
            C_dataToSetId,
            C_setPortfolioType,
			C_mktDataKeys,
            C_isUpdated,
            result,
            objId);
	}

private:
	long    C_tarnId;
	long    C_dataToSetId;
    string  C_setPortfolioType;
	vector<string> C_mktDataKeys;
    bool    C_isUpdated;
};


/////////////////////////////////////////////////////////////
/// central function that does the creation of the XL function
/////////////////////////////////////////////////////////////
LPXLOPER Local_TARNGet_Common(
	LPXLOPER XL_tarnId,
	LPXLOPER XL_getType,
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

		long tarnId;
		XL_GETOBJID( XL_tarnId, tarnId,	" ARM_ERR: TARN Calculator: Object expected",C_result);

		CCString getTypeStr;
		XL_readStrCell(XL_getType,getTypeStr," ARM_ERR: Accessor Type: string expected",C_result);
        char* type=getTypeStr.c_str(); // à cause du new !!
		string getType(type);
        delete type;

		CCString tarnGetClass(GCGetTypeToClass(getType,tarnId).c_str());

		tarnGetDataFunc ourFunc(tarnId,getType);

		/// call the general function
		fillXL_Result( tarnGetClass, ourFunc, C_result, XL_result, PersistentInXL );
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_TARNGet_Common" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}
							 
LPXLOPER Local_TARNSet_Common(
	LPXLOPER XL_tarnId,
	LPXLOPER XL_dataId,
	LPXLOPER XL_setPortfolioType,
	LPXLOPER XL_MktDataKeys,
	bool isUpdated,
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

		long tarnId;
		XL_GETOBJID( XL_tarnId, tarnId,	" ARM_ERR: TARN Calculator: Object expected",C_result);

		long dataId;
		XL_GETOBJID( XL_dataId, dataId,	" ARM_ERR: Object expected",C_result);

		CCString setPortfolioTypeStr;
		XL_readStrCellWD(XL_setPortfolioType,setPortfolioTypeStr,"DEFAULT"," ARM_ERR: Portfolio Type: string expected",C_result);
        char * portType = setPortfolioTypeStr.c_str(); // à cause du new !!
		string setPortfolioType(portType);
        delete portType;

		VECTOR<CCString> mktDataKeys;
		VECTOR<CCString> mktDataKeysDef(0);
		XL_readStrVectorWD(XL_MktDataKeys,mktDataKeys,mktDataKeysDef," ARM_ERR: Market datas keys: array of string expected",DOUBLE_TYPE,C_result);

		vector< string > mktDataKeysSTL(mktDataKeys.size());
		for(size_t i=0;i<mktDataKeys.size();++i)
        {
            mktDataKeys[i].toUpper();
			mktDataKeysSTL[i]=CCSTringToSTLString(mktDataKeys[i]);
        }

        if(isUpdated)
        {
            /// Simple updating of the calculator
		    long retCode = ARMLOCAL_TARN_Set(tarnId,dataId,setPortfolioType,
											mktDataKeysSTL,
											isUpdated,
                                            C_result,ARM_NULL_OBJECT_ID);

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
        else
        {
            /// General call through the functor with an object creation
		    tarnSetDataFunc ourFunc(tarnId,dataId,setPortfolioType,mktDataKeysSTL,isUpdated);
		    fillXL_Result( LOCAL_GC_TARN_CLASS, ourFunc, C_result, XL_result, PersistentInXL );
        }
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_CRFSet_Common" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}
							 


///////////////////////////////////
/// version that takes into account 
/// previous creation of object
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_TARNCalculator_GetData(
	LPXLOPER XL_tarnId,
	LPXLOPER XL_getType)
{
	ADD_LOG("Local_TARNCalculator_GetData");
	bool PersistentInXL = true;
	return Local_TARNGet_Common(
	    XL_tarnId,
	    XL_getType,
        PersistentInXL );
}

__declspec(dllexport) LPXLOPER WINAPI Local_TARNCalculator_SetData(
	LPXLOPER XL_tarnId,
	LPXLOPER XL_dataId,
    LPXLOPER XL_setPortfolioType,
	LPXLOPER XL_MktDataKeys)
{
	ADD_LOG("Local_TARNCalculator_SetData");
	bool PersistentInXL = true;
    bool isUpdated = false;
	return Local_TARNSet_Common(
	    XL_tarnId,
	    XL_dataId,
        XL_setPortfolioType,
		XL_MktDataKeys,
        isUpdated,
        PersistentInXL );
}


///////////////////////////////////
/// version for non persistentInXl
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_PXL_TARNCalculator_GetData(
	LPXLOPER XL_tarnId,
	LPXLOPER XL_getType)
{
	ADD_LOG("Local_PXL_TARNCalculator_GetData");
	bool PersistentInXL = false;
	return Local_TARNGet_Common(
	    XL_tarnId,
	    XL_getType,
        PersistentInXL );
}

__declspec(dllexport) LPXLOPER WINAPI Local_PXL_TARNCalculator_SetData(
	LPXLOPER XL_tarnId,
	LPXLOPER XL_dataId,
    LPXLOPER XL_setPortfolioType,
	LPXLOPER XL_MktDataKeys)
{
	ADD_LOG("Local_PXL_TARNCalculator_SetData");
	bool PersistentInXL = false;
    bool isUpdated = false;
	return Local_TARNSet_Common(
	    XL_tarnId,
	    XL_dataId,
        XL_setPortfolioType,
		XL_MktDataKeys,
        isUpdated,
        PersistentInXL );
}

///////////////////////////////////
/// same as SetData but the previous TARN
/// is not cloned simply updated (=> no PXL version)
/// Only GC datas are updattable
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_TARNCalculator_Update(
	LPXLOPER XL_tarnId,
	LPXLOPER XL_dataId,
    LPXLOPER XL_setPortfolioType,
	LPXLOPER XL_mktDataKeys)
{
	ADD_LOG("Local_TARNCalculator_Update");
	bool PersistentInXL = true;
    bool isUpdated = true;
	return Local_TARNSet_Common(
	    XL_tarnId,
	    XL_dataId,
        XL_setPortfolioType,
		XL_mktDataKeys,
        isUpdated,
        PersistentInXL );
}

__declspec(dllexport) LPXLOPER WINAPI Local_TARNCalculator_GetPricingData(
	LPXLOPER XL_TARNCalculatorId,
    LPXLOPER XL_KeyId )
{
	ADD_LOG("Local_TARNCalculator_GetPricingData");
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
		XL_readStrCell( XL_TARNCalculatorId, C_PricerString, " ARM_ERR: TARN Calculator Id: Object expected",C_result);
		long C_TarnCalcId = LocalGetNumObjectId(C_PricerString);
		
		CCString C_KeyString;
		XL_readStrCell( XL_KeyId, C_KeyString, " ARM_ERR: Key Id: Object expected",C_result);
		
		/// call the function
		ARM_GramFctorArg argResult;
		
		long retCode = ARMLOCAL_GC_GetPricingData(
			C_TarnCalcId,
			CCSTringToSTLString( C_KeyString ),
			argResult,
			C_result );
		
		if( retCode == ARM_KO )
		{
			ARM_ERR();
		}
		else
		{
			ARM_GramFunctorToXLOPER(argResult, XL_result, C_result);
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

///----------------------------------------------
///----------------------------------------------
///             Maturity Cap Calculator
/// Inputs :
///     
///----------------------------------------------
///----------------------------------------------
class maturityCapCalculatorFunc : public ARMResultLong2LongFunc
{
public:
	maturityCapCalculatorFunc(
        double startDate,
        double endDate,
		double underlyingEndDate,
        long longShort,
		long capFloor,
        long resetFreq,
		long payFreq,
        const string& indexTerm,
        long dayCount,
		long intRule,
		double spread,
		double initNominal,
		double initTRI,
		double annuity,
		long maturityCapMode,
		double coeff,
		double amortizing,
		long amortizingId,
		long resetGap,
        const string& resetCal,
        const string& payCal,
        long calibrationMode,
		long nbIterations,
        const vector< string >& flags,
        long mktDataManagerId,
        const vector< string >& keys)
    :
    C_startDate(startDate),
    C_endDate(endDate),
	C_underlyingEndDate(underlyingEndDate),
    C_longShort(longShort),
	C_capFloor(capFloor),
	C_resetFreq(resetFreq),
	C_payFreq(payFreq),
	C_indexTerm(indexTerm),
	C_dayCount(dayCount),
	C_intRule(intRule),
	C_spread(spread),
	C_initNominal(initNominal),
	C_initTRI(initTRI),
	C_annuity(annuity),
	C_maturityCapMode(maturityCapMode),
	C_coeff(coeff),
	C_amortizing(amortizing),
	C_amortizingId(amortizingId),
	C_resetGap(resetGap),
	C_resetCal(resetCal),
	C_payCal(payCal),
	C_calibrationMode(calibrationMode),
	C_nbIterations(nbIterations),
	C_flags(flags),
	C_mktDataManagerId(mktDataManagerId),
	C_keys(keys)
    {};
	
	long operator()( ARM_result& result, long objId ){
		return ARMLOCAL_MaturityCapCalculator_Create(
            C_startDate,
			C_endDate,
			C_underlyingEndDate,
			C_longShort,
			C_capFloor,
			C_resetFreq,
			C_payFreq,
			C_indexTerm,
			C_dayCount,
			C_intRule,
			C_spread,
			C_initNominal,
			C_initTRI,
			C_annuity,
			C_maturityCapMode,
			C_coeff,
			C_amortizing,
			C_amortizingId,
			C_resetGap,
			C_resetCal,
			C_payCal,
			C_calibrationMode,
			C_nbIterations,
			C_flags,
			C_mktDataManagerId,
			C_keys,
            result,
            objId);
    }

private:
    double C_startDate;
    double C_endDate;
	double C_underlyingEndDate;
    long C_longShort;
	long C_capFloor;
	long C_resetFreq;
	long C_payFreq;
	const string& C_indexTerm;
	long C_dayCount;
	long C_intRule;
	double C_spread;
	double C_initNominal;
	double C_initTRI;
	double C_annuity;
	long C_maturityCapMode;
	double C_coeff;
	double C_amortizing;
	long C_amortizingId;
	long C_resetGap;
	const string& C_resetCal;
	const string& C_payCal;
	long C_calibrationMode;
	long C_nbIterations;
	const vector< string >& C_flags;
	long C_mktDataManagerId;
	const vector< string >& C_keys;
};


/////////////////////////////////////////////////////////////
/// central function that does the creation of the XL function
/////////////////////////////////////////////////////////////
LPXLOPER Local_MaturityCapCalculator_Common(
    LPXLOPER XL_startDate,
    LPXLOPER XL_endDate,
	LPXLOPER XL_underlyingEndDate,
	LPXLOPER XL_longShort,
	LPXLOPER XL_capFloor,
	LPXLOPER XL_loanDatas,
    LPXLOPER XL_spread,
    LPXLOPER XL_initNominal,
	LPXLOPER XL_initTRI,
    LPXLOPER XL_annuity,
    LPXLOPER XL_maturityCapMode,
    LPXLOPER XL_coeff,
	LPXLOPER XL_amortizing,
    LPXLOPER XL_resetGap,
    LPXLOPER XL_resetCal,
    LPXLOPER XL_payCal,
	LPXLOPER XL_calibrationMode,
	LPXLOPER XL_nbIterations,
    LPXLOPER XL_flags,
	LPXLOPER XL_mktDatas,
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

    
		/// Maturity Cap datas
		/// -----------------

		/// Start, End, Underlying End Date
		double startDate;
		XL_readNumCell(XL_startDate,startDate," ARM_ERR: start date: date expected",C_result);
		double endDate;
		XL_readNumCell(XL_endDate,endDate," ARM_ERR: end date: date expected",C_result);
		double underlyingEndDate;
		XL_readNumCell(XL_underlyingEndDate,underlyingEndDate," ARM_ERR: underlying end date: date expected",C_result);

		/// Long Short
		CCString longShortStr;
		long longShort;
		XL_readStrCell(XL_longShort,longShortStr," ARM_ERR: long/short: string expected",C_result);
		if((longShort = ARM_ConvLongOrShort (longShortStr, C_result)) == ARM_DEFAULT_ERR)
		{
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}

		/// Cap Floor
		CCString capFloorStr;
		long capFloor;
		XL_readStrCell(XL_capFloor,capFloorStr," ARM_ERR: cap/floor: string expected",C_result);
		if((capFloor = ARM_ConvCapOrFloor (capFloorStr, C_result)) == ARM_DEFAULT_ERR)
		{
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}

		// Loan Datas
		VECTOR<CCString> loanDatas;
		XL_readStrVector (XL_loanDatas, loanDatas, " ARM_ERR: Loan datas: array of string expected", DOUBLE_TYPE, C_result);
		if(loanDatas.size() < 5)
		{
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}

		/// reset frequency 
		CCString resetFreqStr = loanDatas[0];
		long resetFreq;
		if((resetFreq = ARM_ConvFrequency (resetFreqStr, C_result)) == ARM_DEFAULT_ERR)
		{
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}

		/// payment frequency 
		CCString payFreqFreqStr = loanDatas[1];
		long payFreq;
		if((payFreq = ARM_ConvFrequency (payFreqFreqStr, C_result)) == ARM_DEFAULT_ERR)
		{
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}

		/// index term
		CCString indexTermXl = loanDatas[2];
		string indexTerm=CCSTringToSTLString(indexTermXl);

		/// day count
		CCString dayCountStr = loanDatas[3];
		long dayCount;
		if((dayCount = ARM_ConvDayCount (dayCountStr)) == ARM_DEFAULT_ERR)
		{
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}

		/// interest rule
		CCString intRuleStr = loanDatas[4];
		long intRule;
		if ((intRule = ARM_ConvIntRule (intRuleStr)) == ARM_DEFAULT_ERR)
		{
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}

		// spread
		double spread;
		XL_readNumCell(XL_spread,spread," ARM_ERR: spread: double expected",C_result);

		// Initial Nominal
		double initNominal;
		XL_readNumCell(XL_initNominal,initNominal," ARM_ERR: init nominal: double expected",C_result);

		// Initial TRI
		double initTRI;
		XL_readNumCell(XL_initTRI,initTRI," ARM_ERR: init TRI: double expected",C_result);

		// Constant Annuity
		double annuity;
		XL_readNumCell(XL_annuity,annuity," ARM_ERR: annuity: double expected",C_result);

		// Maturity Cap Mode
		CCString maturityCapModeStr;
		XL_readStrCell(XL_maturityCapMode,maturityCapModeStr," ARM_ERR: maturity cap mode: string expected",C_result);
		long maturityCapMode;
		if((maturityCapMode = ARM_ConvMaturityCapMode(maturityCapModeStr, C_result)) == ARM_DEFAULT_ERR)
		{
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}

		// Coefficient
		double coeff;
		double coeffDef = 1.0/resetFreq;
		XL_readNumCellWD(XL_coeff,coeff,coeffDef," ARM_ERR: coeff: double expected",C_result);

		// Notional amortizing
		double amortizing;
		double amortizingDef=100.0;
		CCString amortizingStr;
		long     amortizingId;
		XL_readStrOrNumCellWD(XL_amortizing,amortizingStr,amortizing,amortizingDef,amortizingId,
			" ARM_ERR: amortizing: numerical or curve Id expected",C_result);
		if(amortizingId == XL_TYPE_STRING)
			amortizingId = LocalGetNumObjectId(amortizingStr);
		else
			amortizingId = ARM_NULL_OBJECT;

		/// reset gap
		double resetGapD;
		double resetGapDef=2;
		XL_readNumCellWD(XL_resetGap,resetGapD,resetGapDef," ARM_ERR: reset gap: numerical expected",C_result);
		long resetGap = resetGapD;

		/// reset calendar
		CCString resetCalXl;
		CCString resetCalDef("EUR");
		XL_readStrCellWD(XL_resetCal,resetCalXl,resetCalDef," ARM_ERR: reset calendar: string expected",C_result);
		string resetCal = CCSTringToSTLString(resetCalXl);

		/// payment calendar
		CCString payCalXl;
		CCString payCalDef("EUR");
		XL_readStrCellWD(XL_payCal,payCalXl,payCalDef," ARM_ERR: payment calendar: string expected",C_result);
		string payCal = CCSTringToSTLString(payCalXl);

		// calibration mode
		CCString calibrationModeStr;
		XL_readStrCell(XL_calibrationMode,calibrationModeStr," ARM_ERR: calibration mode: string expected",C_result);
		long calibrationMode;
		if((calibrationMode = ARM_ConvMaturityCapCalibrationMode(calibrationModeStr, C_result)) == ARM_DEFAULT_ERR)
		{
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}

		// nb iterations
		double nbIterations;
		XL_readNumCell(XL_nbIterations,nbIterations," ARM_ERR: nb iterations: double expected",C_result);


		/// Market datas : curve , MRS and market models for OSW & CF pricings
		/// ------------------------------------------------------------------
		VECTOR<CCString> mktDatas;
		XL_readStrVector(XL_mktDatas,mktDatas," ARM_ERR: Market datas: array of string expected",DOUBLE_TYPE,C_result);
        if(mktDatas.size()<1)
		{
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}

		long mktDataManagerId = LocalGetNumObjectId (mktDatas[0]);
		vector< string > keys(mktDatas.size()-1);
        size_t i;
		for(i=1;i<mktDatas.size();++i)
			keys[i-1]=CCSTringToSTLString(mktDatas[i]);

		/// Flags datas
		/// -----------
		VECTOR<CCString> C_flags;
		VECTOR<CCString> flagsDef(4,"Y");
        flagsDef[1]=CCString("N"); // Ref Std Cap
		flagsDef[2]=CCString("N"); // Estimated Strikes
		flagsDef[3]=CCString("N"); // Estimated Nominal
		XL_readStrVectorWD(XL_flags,C_flags,flagsDef," ARM_ERR: flags: array of string expected",DOUBLE_TYPE,C_result);

		vector< string > flags(C_flags.size());
		for(i=0;i<C_flags.size();++i)
        {
            (C_flags[i]).toUpper();
			flags[i]=CCSTringToSTLString(C_flags[i]);
        }

		/// use the concept of Functor to transfer the knowledge of
		/// a function with a context
		maturityCapCalculatorFunc ourFunc(
				startDate,
				endDate,
				underlyingEndDate,
				longShort,
				capFloor,
				resetFreq,
				payFreq,
				indexTerm,
				dayCount,
				intRule,
				spread,
				initNominal,
				initTRI,
				annuity,
				maturityCapMode,
				coeff,
				amortizing,
				amortizingId,
				resetGap,
				resetCal,
				payCal,
				calibrationMode,
				nbIterations,
                flags,
				mktDataManagerId,
				keys);

		/// call the general function
		fillXL_Result( LOCAL_GC_MATURITYCAP_CLASS, ourFunc, C_result, XL_result, PersistentInXL );
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_MaturityCapCalculator_Common" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}
							 

///////////////////////////////////
/// version that takes into account 
/// previous creation of object
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_MaturityCapCalculator_Create(
    LPXLOPER XL_startDate,
    LPXLOPER XL_endDate,
	LPXLOPER XL_underlyingEndDate,
	LPXLOPER XL_longShort,
	LPXLOPER XL_capFloor,
	LPXLOPER XL_loanDatas,
    LPXLOPER XL_spread,
    LPXLOPER XL_initNominal,
	LPXLOPER XL_initTRI,
    LPXLOPER XL_annuity,
    LPXLOPER XL_maturityCapMode,
    LPXLOPER XL_coeff,
	LPXLOPER XL_amortizing,
    LPXLOPER XL_resetGap,
    LPXLOPER XL_resetCal,
    LPXLOPER XL_payCal,
	LPXLOPER XL_calibrationMode,
	LPXLOPER XL_nbIterations,
    LPXLOPER XL_flags,
	LPXLOPER XL_mktDatas)
{
	ADD_LOG("Local_MaturityCapCalculator_Create");
	bool PersistentInXL = true;
	return Local_MaturityCapCalculator_Common(
        XL_startDate,
		XL_endDate,
		XL_underlyingEndDate,
		XL_longShort,
		XL_capFloor,
		XL_loanDatas,
		XL_spread,
		XL_initNominal,
		XL_initTRI,
		XL_annuity,
		XL_maturityCapMode,
		XL_coeff,
		XL_amortizing,
		XL_resetGap,
		XL_resetCal,
		XL_payCal,
		XL_calibrationMode,
		XL_nbIterations,
		XL_flags,
		XL_mktDatas,
		PersistentInXL );
}


///////////////////////////////////
/// version for non persistentInXl
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_PXL_MaturityCapCalculator_Create(
    LPXLOPER XL_startDate,
    LPXLOPER XL_endDate,
	LPXLOPER XL_underlyingEndDate,
	LPXLOPER XL_longShort,
	LPXLOPER XL_capFloor,
	LPXLOPER XL_loanDatas,
    LPXLOPER XL_spread,
    LPXLOPER XL_initNominal,
	LPXLOPER XL_initTRI,
    LPXLOPER XL_annuity,
    LPXLOPER XL_maturityCapMode,
    LPXLOPER XL_coeff,
	LPXLOPER XL_amortizing,
    LPXLOPER XL_resetGap,
    LPXLOPER XL_resetCal,
    LPXLOPER XL_payCal,
	LPXLOPER XL_calibrationMode,
	LPXLOPER XL_nbIterations,
    LPXLOPER XL_flags,
	LPXLOPER XL_mktDatas)
{
	ADD_LOG("Local_PXL_MaturityCapCalculator_Create");
	bool PersistentInXL = false;
	return Local_MaturityCapCalculator_Common(
        XL_startDate,
		XL_endDate,
		XL_underlyingEndDate,
		XL_longShort,
		XL_capFloor,
		XL_loanDatas,
		XL_spread,
		XL_initNominal,
		XL_initTRI,
		XL_annuity,
		XL_maturityCapMode,
		XL_coeff,
		XL_amortizing,
		XL_resetGap,
		XL_resetCal,
		XL_payCal,
		XL_calibrationMode,
		XL_nbIterations,
		XL_flags,
		XL_mktDatas,
		PersistentInXL);
}

///----------------------------------------------
///----------------------------------------------
///             Maturity Cap Calculator Accessor
/// Inputs :
///     
///----------------------------------------------
///----------------------------------------------
class maturityCapGetDataFunc : public ARMResultLong2LongFunc
{
public:
	maturityCapGetDataFunc(
        long maturityCapId,
        const string& getType)
    :
    C_maturityCapId(maturityCapId),
    C_getType(getType)
    {};
	
	long operator()( ARM_result& result, long objId ){
		return ARMLOCAL_MaturityCap_Get(
            C_maturityCapId,
            C_getType,
            result,
            objId);
    }

private:
	long    C_maturityCapId;
	string  C_getType;
};

class maturityCapSetDataFunc : public ARMResultLong2LongFunc
{
public:
	maturityCapSetDataFunc(
        long maturityCapId,
        long dataToSetId,
        const string& setPortfolioType,
		const vector<string>& mktDataKeys,
        bool isUpdated)
    :
    C_maturityCapId(maturityCapId),
    C_dataToSetId(dataToSetId),
    C_setPortfolioType(setPortfolioType),
	C_mktDataKeys(mktDataKeys),
    C_isUpdated(isUpdated)
    {};
	
	long operator()( ARM_result& result, long objId ){
		return ARMLOCAL_MaturityCap_Set(
            C_maturityCapId,
            C_dataToSetId,
            C_setPortfolioType,
			C_mktDataKeys,
            C_isUpdated,
            result,
            objId);
	}

private:
	long    C_maturityCapId;
	long    C_dataToSetId;
    string  C_setPortfolioType;
	vector<string> C_mktDataKeys;
    bool    C_isUpdated;
};


/////////////////////////////////////////////////////////////
/// central function that does the creation of the XL function
/////////////////////////////////////////////////////////////
LPXLOPER Local_MaturityCapGet_Common(
	LPXLOPER XL_maturityCapId,
	LPXLOPER XL_getType,
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

		long maturityCapId;
		XL_GETOBJID( XL_maturityCapId, maturityCapId,	" ARM_ERR: Maturity Cap Calculator: Object expected",C_result);

		CCString getTypeStr;
		XL_readStrCell(XL_getType,getTypeStr," ARM_ERR: Accessor Type: string expected",C_result);
        char* type=getTypeStr.c_str(); // à cause du new !!
		string getType(type);
        delete type;

		CCString maturityCapGetClass(GCGetTypeToClass(getType,maturityCapId).c_str());

		maturityCapGetDataFunc ourFunc(maturityCapId,getType);

		/// call the general function
		fillXL_Result( maturityCapGetClass, ourFunc, C_result, XL_result, PersistentInXL );
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_MaturityCapGet_Common" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}
							 
LPXLOPER Local_MaturityCapSet_Common(
	LPXLOPER XL_maturityCapId,
	LPXLOPER XL_dataId,
	LPXLOPER XL_setPortfolioType,
	LPXLOPER XL_MktDataKeys,
	bool isUpdated,
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

		long maturityCapId;
		XL_GETOBJID( XL_maturityCapId, maturityCapId,	" ARM_ERR: Maturity Cap Calculator: Object expected",C_result);

		long dataId;
		XL_GETOBJID( XL_dataId, dataId,	" ARM_ERR: Object expected",C_result);

		CCString setPortfolioTypeStr;
		XL_readStrCellWD(XL_setPortfolioType,setPortfolioTypeStr,"DEFAULT"," ARM_ERR: Portfolio Type: string expected",C_result);
        char * portType = setPortfolioTypeStr.c_str(); // à cause du new !!
		string setPortfolioType(portType);
        delete portType;

		VECTOR<CCString> mktDataKeys;
		VECTOR<CCString> mktDataKeysDef(0);
		XL_readStrVectorWD(XL_MktDataKeys,mktDataKeys,mktDataKeysDef," ARM_ERR: Market datas keys: array of string expected",DOUBLE_TYPE,C_result);

		vector< string > mktDataKeysSTL(mktDataKeys.size());
		for(size_t i=0;i<mktDataKeys.size();++i)
        {
            mktDataKeys[i].toUpper();
			mktDataKeysSTL[i]=CCSTringToSTLString(mktDataKeys[i]);
        }

        if(isUpdated)
        {
            /// Simple updating of the calculator
		    long retCode = ARMLOCAL_MaturityCap_Set(maturityCapId,dataId,setPortfolioType,
											mktDataKeysSTL,
											isUpdated,
                                            C_result,ARM_NULL_OBJECT_ID);

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
        else
        {
            /// General call through the functor with an object creation
		    maturityCapSetDataFunc ourFunc(maturityCapId,dataId,setPortfolioType,mktDataKeysSTL,isUpdated);
		    fillXL_Result( LOCAL_GC_MATURITYCAP_CLASS, ourFunc, C_result, XL_result, PersistentInXL );
        }
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_CRFSet_Common" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}
							 


///////////////////////////////////
/// version that takes into account 
/// previous creation of object
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_MaturityCapCalculator_GetData(
	LPXLOPER XL_maturityCapId,
	LPXLOPER XL_getType)
{
	ADD_LOG("Local_MaturityCapCalculator_GetData");
	bool PersistentInXL = true;
	return Local_MaturityCapGet_Common(
	    XL_maturityCapId,
	    XL_getType,
        PersistentInXL );
}

__declspec(dllexport) LPXLOPER WINAPI Local_MaturityCapCalculator_SetData(
	LPXLOPER XL_maturityCapId,
	LPXLOPER XL_dataId,
    LPXLOPER XL_setPortfolioType,
	LPXLOPER XL_MktDataKeys)
{
	ADD_LOG("Local_MaturityCapCalculator_SetData");
	bool PersistentInXL = true;
    bool isUpdated = false;
	return Local_MaturityCapSet_Common(
	    XL_maturityCapId,
	    XL_dataId,
        XL_setPortfolioType,
		XL_MktDataKeys,
        isUpdated,
        PersistentInXL );
}


///////////////////////////////////
/// version for non persistentInXl
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_PXL_MaturityCapCalculator_GetData(
	LPXLOPER XL_maturityCapId,
	LPXLOPER XL_getType)
{
	ADD_LOG("Local_PXL_MaturityCapCalculator_GetData");
	bool PersistentInXL = false;
	return Local_MaturityCapGet_Common(
	    XL_maturityCapId,
	    XL_getType,
        PersistentInXL );
}

__declspec(dllexport) LPXLOPER WINAPI Local_PXL_MaturityCapCalculator_SetData(
	LPXLOPER XL_maturityCapId,
	LPXLOPER XL_dataId,
    LPXLOPER XL_setPortfolioType,
	LPXLOPER XL_MktDataKeys)
{
	ADD_LOG("Local_PXL_MaturityCapCalculator_SetData");
	bool PersistentInXL = false;
    bool isUpdated = false;
	return Local_MaturityCapSet_Common(
	    XL_maturityCapId,
	    XL_dataId,
        XL_setPortfolioType,
		XL_MktDataKeys,
        isUpdated,
        PersistentInXL );
}

///////////////////////////////////
/// same as SetData but the previous MaturityCap
/// is not cloned simply updated (=> no PXL version)
/// Only GC datas are updattable
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_MaturityCapCalculator_Update(
	LPXLOPER XL_maturityCapId,
	LPXLOPER XL_dataId,
    LPXLOPER XL_setPortfolioType,
	LPXLOPER XL_mktDataKeys)
{
	ADD_LOG("Local_MaturityCapCalculator_Update");
	bool PersistentInXL = true;
    bool isUpdated = true;
	return Local_MaturityCapSet_Common(
	    XL_maturityCapId,
	    XL_dataId,
        XL_setPortfolioType,
		XL_mktDataKeys,
        isUpdated,
        PersistentInXL );
}

__declspec(dllexport) LPXLOPER WINAPI Local_MaturityCapCalculator_GetPricingData(
	LPXLOPER XL_maturityCapId,
    LPXLOPER XL_KeyId )
{
	ADD_LOG("Local_MaturityCapCalculator_GetPricingData");
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
		XL_readStrCell( XL_maturityCapId, C_PricerString, " ARM_ERR: Maturity Cap Calculator Id: Object expected",C_result);
		long C_maturityCapId = LocalGetNumObjectId(C_PricerString);
		
		CCString C_KeyString;
		XL_readStrCell( XL_KeyId, C_KeyString, " ARM_ERR: Key Id: Object expected",C_result);
		
		/// call the function
		ARM_GramFctorArg argResult;
		
		long retCode = ARMLOCAL_GC_GetPricingData(
			C_maturityCapId,
			CCSTringToSTLString( C_KeyString ),
			argResult,
			C_result );
		
		if( retCode == ARM_KO )
		{
			ARM_ERR();
		}
		else
		{
			ARM_GramFunctorToXLOPER(argResult, XL_result, C_result);
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

///////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////Caption Calculator/////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////
///----------------------------------------------
///----------------------------------------------
///             Caption Calculator
/// Inputs :
///     
///----------------------------------------------
///----------------------------------------------
class captionCalculatorFunc : public ARMResultLong2LongFunc
{
public:
	captionCalculatorFunc(
        double startDate,
        double endDate,
		const string& cpnIdxTerm,
		long payRec,
		long CF,
		double coupon,
        long couponId,
		const string& FundIdxTerm,
		long NotifDays,
		long NonCall,
		double exercise,
		long exerStyleId,
		double notional,
        long notionalId,
		long cpnDayCount,
		long cpnResetTiming,
		const string& cpnResetCal,
        const string& cpnPayCal,
		double cpnSpread,
        long cpnSpreadId,
		long fundDayCount,
		const string& fundResetCal,
        const string& fundPayCal,
		double fundSpread,
        long fundSpreadId,
		long factorNb,
		const vector< string >& CalibMod,
		const vector< string >& flags,
		long mktDataManagerId,
        const vector< string >& keys)
    :
	C_startDate(startDate),
    C_endDate(endDate),
	C_cpnIdxTerm(cpnIdxTerm),
	C_payRec(payRec),
	C_CF(CF),
	C_coupon(coupon),
    C_couponId(couponId),
	C_FundIdxTerm(FundIdxTerm),
	C_NotifDays(NotifDays),
	C_NonCall(NonCall),
	C_exercise(exercise),
	C_exerStyleId(exerStyleId),
	C_notional(notional),
    C_notionalId(notionalId),
	C_cpnDayCount(cpnDayCount),
	C_cpnResetTiming(cpnResetTiming),
	C_cpnResetCal(cpnResetCal),
    C_cpnPayCal(cpnPayCal),
	C_cpnSpread(cpnSpread),
    C_cpnSpreadId(cpnSpreadId),
	C_fundDayCount(fundDayCount),
	C_fundResetCal(fundResetCal),
    C_fundPayCal(fundPayCal),
	C_fundSpread(fundSpread),
    C_fundSpreadId(fundSpreadId),
	C_factorNb(factorNb),
	C_CalibMod(CalibMod),
	C_flags(flags),
	C_mktDataManagerId(mktDataManagerId),
    C_keys(keys)
    {};
	long operator()( ARM_result& result, long objId ){
		return ARMLOCAL_CaptionCalculator_Create(
            C_startDate,
			C_endDate,
			C_cpnIdxTerm,
			C_payRec,
			C_CF,
			C_coupon,
			C_couponId,
			C_FundIdxTerm,
			C_NotifDays,
			C_NonCall,
			C_exercise,
			C_exerStyleId,
			C_notional,
			C_notionalId,
			C_cpnDayCount,
			C_cpnResetTiming,
			C_cpnResetCal,
			C_cpnPayCal,
			C_cpnSpread,
			C_cpnSpreadId,
			C_fundDayCount,
			C_fundResetCal,
			C_fundPayCal,
			C_fundSpread,
			C_fundSpreadId,
			C_factorNb,
			C_CalibMod,
			C_flags,
			C_mktDataManagerId,
			C_keys,
            result,
            objId);
    }

private:

	double				C_startDate;
	double				C_endDate;
	string				C_cpnIdxTerm;
	long				C_payRec;
	long				C_CF;
	double				C_coupon;
	long				C_couponId;
	string				C_FundIdxTerm;
	long				C_NotifDays;
	long				C_NonCall;
	double				C_exercise;
	long				C_exerStyleId;
	double				C_notional;
	long				C_notionalId;
	long				C_cpnDayCount;
	long				C_cpnResetTiming;
	string				C_cpnResetCal;
	string				C_cpnPayCal;
	double				C_cpnSpread;
	long				C_cpnSpreadId;
	long				C_fundDayCount;
	string				C_fundResetCal;
	string				C_fundPayCal;
	double				C_fundSpread;
	long				C_fundSpreadId;
	long				C_factorNb;
	vector< string >	C_CalibMod;
	vector< string >	C_flags;
	long				C_mktDataManagerId;
	vector< string >	C_keys;
};


/////////////////////////////////////////////////////////////
/// central function that does the creation of the XL function
/////////////////////////////////////////////////////////////
LPXLOPER Local_CaptionCalculator_Common(
    LPXLOPER XL_startDate,
    LPXLOPER XL_endDate,
	LPXLOPER XL_mktDatas,
	LPXLOPER XL_RequiredCpnDatas,
	LPXLOPER XL_coupon,
	LPXLOPER XL_FundIdxTerm,
	LPXLOPER XL_ExerciseData,
	LPXLOPER XL_exerStyle,
	LPXLOPER XL_Notional,
	LPXLOPER XL_cpnDatas,
	LPXLOPER XL_cpnSpread,
	LPXLOPER XL_fundDatas,
	LPXLOPER XL_fundSpread,
	LPXLOPER XL_factorNb,
	LPXLOPER XL_CalibMod,
	LPXLOPER XL_flags,
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

		///-------------------------------------------------------------------
		///                      Required Caption DATA					   ///
		/// ------------------------------------------------------------------
    
		/// General Caption datas
		/// -----------------
		
		/// Start & End
		double startDate;
		XL_readNumCell(XL_startDate,startDate," ARM_ERR: start date: date expected",C_result);
		double endDate;
		XL_readNumCell(XL_endDate,endDate," ARM_ERR: end date: date expected",C_result);

		
		/// Market datas : zc curve , MRS, beta, volatility curve and market models for OSW & CF pricings
		VECTOR<CCString> mktDatas;
		XL_readStrVector(XL_mktDatas,mktDatas," ARM_ERR: Market datas: array of string expected",DOUBLE_TYPE,C_result);
        if(mktDatas.size()<1)
		{
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}

		long mktDataManagerId = LocalGetNumObjectId (mktDatas[0]);
		vector< string > keys(mktDatas.size()-1);
        size_t i;
		for(i=1;i<mktDatas.size();++i)
			keys[i-1]=CCSTringToSTLString(mktDatas[i]);


		VECTOR<CCString> RequiredCpnDatas;
		XL_readStrVector(XL_RequiredCpnDatas,RequiredCpnDatas," ARM_ERR: Required Cpn Datas : array of string expected",DOUBLE_TYPE,C_result);
        
		if(RequiredCpnDatas.size()!=3)
		{
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}

		/// cpnIdxTerm
		CCString cpnIdxTermStr = RequiredCpnDatas[0];
		string cpnIdxTerm=CCSTringToSTLString(cpnIdxTermStr);

		/// Pay/Rec
		CCString payRecStr = RequiredCpnDatas[1];
		long payRec;
		if((payRec = ARM_ConvRecOrPay (payRecStr, C_result)) == ARM_DEFAULT_ERR)
		{
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}

		//Cap/Floor
		CCString CFStr = RequiredCpnDatas[2];
		long CF;
		if((CF = ARM_ConvCapOrFloor (CFStr, C_result)) == ARM_DEFAULT_ERR)
		{
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}
		
		//coupon
		double coupon;
		CCString couponStr;
		long     couponId;
		XL_readStrOrNumCell(XL_coupon, couponStr, coupon, couponId,
			   " ARM_ERR: coupon: numeric or curve Id expected",C_result);	
		if(couponId == XL_TYPE_STRING)
			couponId = LocalGetNumObjectId(couponStr);
		else
			couponId = ARM_NULL_OBJECT;

		/// fundIdxTerm
		CCString fundIdxTermStr;
		XL_readStrCell(XL_FundIdxTerm,fundIdxTermStr," ARM_ERR: Funding Index Term: string expected",C_result);
		string FundIdxTerm=CCSTringToSTLString(fundIdxTermStr);

		
		///-------------------------------------------------------------------
		///                      Optional Caption DATA					   ///
		/// ------------------------------------------------------------------
    	
		VECTOR<double> ExerciseData;
		VECTOR<double> ExerciseData_default(2,0.0);
		ExerciseData_default[0]=GETDEFAULTVALUE;
		XL_readNumVectorWD(XL_ExerciseData,ExerciseData,ExerciseData_default," ARM_ERR: Required Exercise datas: array of numeric expected", C_result);

		long NotifDays = ExerciseData[0];
		long NonCall = ExerciseData[1];

		/// exerStyle
		double exercise;
		double default_exerStyle=0.0;
		CCString exerStyleStr;
		long     exerStyleId;
		XL_readStrOrNumCellWD(XL_exerStyle,exerStyleStr,exercise,default_exerStyle,exerStyleId,
			" ARM_ERR: exercise Style: numerical or refValue Id expected",C_result);
		if(exerStyleId == XL_TYPE_STRING)
			exerStyleId = LocalGetNumObjectId(exerStyleStr);
		else
			exerStyleId = ARM_NULL_OBJECT;
		
		//Notional		
		double notional;
		double default_notional=100.0;
		CCString notionalStr;
		long     notionalId;
		XL_readStrOrNumCellWD(XL_Notional,notionalStr,notional,default_notional,notionalId,
			" ARM_ERR: notional: numerical or Curve Id expected",C_result);
		if(notionalId == XL_TYPE_STRING)
			notionalId = LocalGetNumObjectId(notionalStr);
		else
			notionalId = ARM_NULL_OBJECT;
		
	
		long default_Longdata = GETDEFAULTVALUE;
		string default_Stringdata = GETDEFAULTVALUESTR;
	/// ------------------------
	/// coupon datas
	/// ------------------------
		VECTOR<CCString> cpnDatas;
		VECTOR<CCString> cpnDatas_default(0);
		XL_readStrVectorWD (XL_cpnDatas,cpnDatas,cpnDatas_default," ARM_ERR: Coupon datas: array of string expected",DOUBLE_TYPE,C_result);
		
		/// coupon daycount
		long cpnDayCount = default_Longdata;
		if(cpnDatas.size() >= 1)
		{
			CCString cpnDayCountXl=cpnDatas[0];
			cpnDayCount = ARM_ConvDayCount (cpnDayCountXl);
		}
			
		// coupon Reset Timing
		long cpnResetTiming = default_Longdata;
		if(cpnDatas.size() >= 2)
		{
			CCString cpnResetTimingXl=cpnDatas[1];
			cpnResetTiming = ARM_ConvPayResetRule (cpnResetTimingXl);
		}

		//coupon reset calendar
		string cpnResetCal = default_Stringdata;
		if(cpnDatas.size() >= 3)
		{
			CCString cpnResetCalXl=cpnDatas[2];
			cpnResetCal = CCSTringToSTLString (cpnResetCalXl);
		}

		//coupon pay calendar
		string cpnPayCal = default_Stringdata;
		if(cpnDatas.size() >= 4)
		{
			CCString cpnPayCalXl=cpnDatas[3];
			cpnPayCal = CCSTringToSTLString (cpnPayCalXl);
		}

		//coupon spread 
		double cpnSpread;
		double default_cpnSpread=0.0;
		CCString cpnSpreadStr;
		long     cpnSpreadId;
		XL_readStrOrNumCellWD(XL_cpnSpread,cpnSpreadStr,cpnSpread,default_cpnSpread,cpnSpreadId,
			" ARM_ERR: cpn spread: numerical or curve Id expected",C_result);

		if(cpnSpreadId == XL_TYPE_STRING)
			cpnSpreadId = LocalGetNumObjectId(cpnSpreadStr);
		else
			cpnSpreadId = ARM_NULL_OBJECT;
	
	/// ------------------------
	/// Funding datas
	/// ------------------------
		VECTOR<CCString> fundDatas;
		VECTOR<CCString> fundDatas_default(0);
		XL_readStrVectorWD (XL_fundDatas,fundDatas,fundDatas_default," ARM_ERR: Funding datas: array of string expected",DOUBLE_TYPE,C_result);
	
		/// funding daycount
		long fundDayCount = default_Longdata;
		if(fundDatas.size() >= 1)
		{
			CCString fundDayCountXl=fundDatas[0];
			fundDayCount = ARM_ConvDayCount (fundDayCountXl);
		}
			
		//Funding reset Calendar
		string fundResetCal = default_Stringdata;
		if(fundDatas.size() >= 2)
		{
			CCString fundResetCalXl=fundDatas[1];
			fundResetCal = CCSTringToSTLString (fundResetCalXl);
		}

		//funding pay calendar
		string fundPayCal = default_Stringdata;
		if(fundDatas.size() >= 3)
		{
			CCString fundPayCalXl=fundDatas[2];
			fundPayCal = CCSTringToSTLString (fundPayCalXl);
		}

		//funding spread 
		double fundSpread;
		double default_fundSpread=0.0;
		CCString fundSpreadStr;
		long     fundSpreadId;
		XL_readStrOrNumCellWD(XL_fundSpread,fundSpreadStr,fundSpread,default_fundSpread,fundSpreadId,
			" ARM_ERR: cpn spread: numerical or curve Id expected",C_result);
		if(fundSpreadId == XL_TYPE_STRING)
			fundSpreadId = LocalGetNumObjectId(fundSpreadStr);
		else
			fundSpreadId = ARM_NULL_OBJECT;

	/// ------------------------
	/// Model datas
	/// ------------------------
	
		//FactorNb
		double factorNbD;
		double default_factorNb = 1.0;
		XL_readNumCellWD(XL_factorNb,factorNbD,default_factorNb," ARM_ERR: fatctorNb : numerical expected",C_result);
		long factorNb = factorNbD;
	


		//CalibMod
		VECTOR<CCString> CalibModXL;
		VECTOR<CCString> CalibMod_default(3);
		CalibMod_default[0]= CCString("DIAG");
		CalibMod_default[1]= CCString("EXSWOPT"); //calibrate mean reversion only on swaption at exercise
		CalibMod_default[2]= CCString("N"); //calibration of beta

		XL_readStrVectorWD (XL_CalibMod,CalibModXL,CalibMod_default," ARM_ERR: Funding datas: array of string expected",DOUBLE_TYPE,C_result);
		vector< string > CalibMod(CalibModXL.size());
		for(i=0;i<CalibModXL.size();++i)
        {            
			CalibMod[i]=CCSTringToSTLString(CalibModXL[i]);
        }
				



	/// ------------------------
	/// Flags 
	/// ------------------------
		VECTOR<CCString> flagsXL;
		VECTOR<CCString> default_flags(5,"N");
        XL_readStrVectorWD(XL_flags,flagsXL,default_flags," ARM_ERR: flags: array of string expected",DOUBLE_TYPE,C_result);

		vector< string > flags(flagsXL.size());
		for(i=0;i<flagsXL.size();++i)
        {
            (flagsXL[i]).toUpper();
			flags[i]=CCSTringToSTLString(flagsXL[i]);
        }

		captionCalculatorFunc ourFunc(
				startDate,
				endDate,
				cpnIdxTerm,
				payRec,
				CF,
				coupon,
				couponId,
				FundIdxTerm,
				NotifDays,
				NonCall,
				exercise,
				exerStyleId,
				notional,
				notionalId,
				cpnDayCount,
				cpnResetTiming,
				cpnResetCal,
				cpnPayCal,
				cpnSpread,
				cpnSpreadId,
				fundDayCount,
				fundResetCal,
				fundPayCal,
				fundSpread,
				fundSpreadId,
				factorNb,
				CalibMod,
				flags,
				mktDataManagerId,
				keys);

		/// call the general function
		fillXL_Result( LOCAL_GC_CAPTION_CLASS, ourFunc, C_result, XL_result, PersistentInXL );
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_CaptionCalculator_Common" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}
							 

///////////////////////////////////
/// version that takes into account 
/// previous creation of object
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_CaptionCalculator_Create(
    LPXLOPER XL_startDate,
    LPXLOPER XL_endDate,
	LPXLOPER XL_mktDatas,
	LPXLOPER XL_RequiredCpnDatas,
	LPXLOPER XL_coupon,
	LPXLOPER XL_FundIdxTerm,
	LPXLOPER XL_ExerciseData,
	LPXLOPER XL_exerStyle,
	LPXLOPER XL_Notional,
	LPXLOPER XL_cpnDatas,
	LPXLOPER XL_cpnSpread,
	LPXLOPER XL_fundDatas,
	LPXLOPER XL_fundSpread,
	LPXLOPER XL_factorNb,
	LPXLOPER XL_CalibMod,
	LPXLOPER XL_flags)
{	bool PersistentInXL = true;
ADD_LOG("Local_CaptionCalculator_Create");
	return Local_CaptionCalculator_Common(
		XL_startDate,
		XL_endDate,
		XL_mktDatas,
		XL_RequiredCpnDatas,
		XL_coupon,
		XL_FundIdxTerm,
		XL_ExerciseData,
		XL_exerStyle,
		XL_Notional,
		XL_cpnDatas,
		XL_cpnSpread,
		XL_fundDatas,
		XL_fundSpread,
		XL_factorNb,
		XL_CalibMod,
		XL_flags,
		PersistentInXL );
}


///////////////////////////////////
/// version for non persistentInXl
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_PXL_CaptionCalculator_Create(
    LPXLOPER XL_startDate,
    LPXLOPER XL_endDate,
	LPXLOPER XL_mktDatas,
	LPXLOPER XL_RequiredCpnDatas,
	LPXLOPER XL_coupon,
	LPXLOPER XL_FundIdxTerm,
	LPXLOPER XL_ExerciseData,
	LPXLOPER XL_exerStyle,
	LPXLOPER XL_Notional,
	LPXLOPER XL_cpnDatas,
	LPXLOPER XL_cpnSpread,
	LPXLOPER XL_fundDatas,
	LPXLOPER XL_fundSpread,
	LPXLOPER XL_factorNb,
	LPXLOPER XL_CalibMod,
	LPXLOPER XL_flags)
{
	ADD_LOG("Local_PXL_CaptionCalculator_Create");
	bool PersistentInXL = false;
	return Local_CaptionCalculator_Common(
        XL_startDate,
		XL_endDate,
		XL_mktDatas,
		XL_RequiredCpnDatas,
		XL_coupon,
		XL_FundIdxTerm,
		XL_ExerciseData,
		XL_exerStyle,
		XL_Notional,
		XL_cpnDatas,
		XL_cpnSpread,
		XL_fundDatas,
		XL_fundSpread,
		XL_factorNb,
		XL_CalibMod,
		XL_flags,
		PersistentInXL );
}

///----------------------------------------------
///----------------------------------------------
///             Caption Calculator Accessor
/// Inputs :
///     
///----------------------------------------------
///----------------------------------------------
class captionGetDataFunc : public ARMResultLong2LongFunc
{
public:
	captionGetDataFunc(
        long captionId,
        const string& getType)
    :
    C_captionId(captionId),
    C_getType(getType)
    {};
	
	long operator()( ARM_result& result, long objId ){
		return ARMLOCAL_Caption_Get(
            C_captionId,
            C_getType,
            result,
            objId);
    }

private:
	long    C_captionId;
	string  C_getType;
};

class captionSetDataFunc : public ARMResultLong2LongFunc
{
public:
	captionSetDataFunc(
        long captionId,
        long dataToSetId,
        const string& setPortfolioType,
		const vector<string>& mktDataKeys,
        bool isUpdated)
    :
    C_captionId(captionId),
    C_dataToSetId(dataToSetId),
    C_setPortfolioType(setPortfolioType),
	C_mktDataKeys(mktDataKeys),
    C_isUpdated(isUpdated)
    {};
	
	long operator()( ARM_result& result, long objId ){
		return ARMLOCAL_Caption_Set(
            C_captionId,
            C_dataToSetId,
            C_setPortfolioType,
			C_mktDataKeys,
            C_isUpdated,
            result,
            objId);
	}

private:
	long    C_captionId;
	long    C_dataToSetId;
    string  C_setPortfolioType;
	vector<string> C_mktDataKeys;
    bool    C_isUpdated;
};


/////////////////////////////////////////////////////////////
/// central function that does the creation of the XL function
/////////////////////////////////////////////////////////////
LPXLOPER Local_CaptionGet_Common(
	LPXLOPER XL_captionId,
	LPXLOPER XL_getType,
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

		long captionId;
		XL_GETOBJID( XL_captionId, captionId,	" ARM_ERR: caption Calculator: Object expected",C_result);

		CCString getTypeStr;
		XL_readStrCell(XL_getType,getTypeStr," ARM_ERR: Accessor Type: string expected",C_result);
        char* type=getTypeStr.c_str(); // à cause du new !!
		string getType(type);
        delete type;

		CCString captionGetClass(GCGetTypeToClass(getType,captionId).c_str());

		captionGetDataFunc ourFunc(captionId,getType);

		/// call the general function
		fillXL_Result( captionGetClass, ourFunc, C_result, XL_result, PersistentInXL );
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_CaptionGet_Common" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}
							 
LPXLOPER Local_CaptionSet_Common(
	LPXLOPER XL_captionId,
	LPXLOPER XL_dataId,
	LPXLOPER XL_setPortfolioType,
	LPXLOPER XL_MktDataKeys,
	bool isUpdated,
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

		long captionId;
		XL_GETOBJID( XL_captionId, captionId,	" ARM_ERR: Caption Calculator: Object expected",C_result);

		long dataId;
		XL_GETOBJID( XL_dataId, dataId,	" ARM_ERR: Object expected",C_result);

		CCString setPortfolioTypeStr;
		XL_readStrCellWD(XL_setPortfolioType,setPortfolioTypeStr,"DEFAULT"," ARM_ERR: Portfolio Type: string expected",C_result);
        char * portType = setPortfolioTypeStr.c_str(); // à cause du new !!
		string setPortfolioType(portType);
        delete portType;

		VECTOR<CCString> mktDataKeys;
		VECTOR<CCString> mktDataKeysDef(0);
		XL_readStrVectorWD(XL_MktDataKeys,mktDataKeys,mktDataKeysDef," ARM_ERR: Market datas keys: array of string expected",DOUBLE_TYPE,C_result);

		vector< string > mktDataKeysSTL(mktDataKeys.size());
		for(size_t i=0;i<mktDataKeys.size();++i)
        {
            mktDataKeys[i].toUpper();
			mktDataKeysSTL[i]=CCSTringToSTLString(mktDataKeys[i]);
        }

        if(isUpdated)
        {
            /// Simple updating of the calculator
		    long retCode = ARMLOCAL_Caption_Set(captionId,dataId,setPortfolioType,
											mktDataKeysSTL,
											isUpdated,
                                            C_result,ARM_NULL_OBJECT_ID);

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
        else
        {
            /// General call through the functor with an object creation
		    captionSetDataFunc ourFunc(captionId,dataId,setPortfolioType,mktDataKeysSTL,isUpdated);
		    fillXL_Result( LOCAL_GC_CAPTION_CLASS, ourFunc, C_result, XL_result, PersistentInXL );
        }
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_CRFSet_Common" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}
							 


///////////////////////////////////
/// version that takes into account 
/// previous creation of object
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_CaptionCalculator_GetData(
	LPXLOPER XL_captionId,
	LPXLOPER XL_getType)
{
	ADD_LOG("Local_CaptionCalculator_GetData");
	bool PersistentInXL = true;
	return Local_CaptionGet_Common(
	    XL_captionId,
	    XL_getType,
        PersistentInXL );
}

__declspec(dllexport) LPXLOPER WINAPI Local_CaptionCalculator_SetData(
	LPXLOPER XL_captionId,
	LPXLOPER XL_dataId,
    LPXLOPER XL_setPortfolioType,
	LPXLOPER XL_MktDataKeys)
{
	ADD_LOG("Local_CaptionCalculator_SetData");
	bool PersistentInXL = true;
    bool isUpdated = false;
	return Local_CaptionSet_Common(
	    XL_captionId,
	    XL_dataId,
        XL_setPortfolioType,
		XL_MktDataKeys,
        isUpdated,
        PersistentInXL );
}


///////////////////////////////////
/// version for non persistentInXl
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_PXL_CaptionCalculator_GetData(
	LPXLOPER XL_captionId,
	LPXLOPER XL_getType)
{
	ADD_LOG("Local_PXL_CaptionCalculator_GetData");
	bool PersistentInXL = false;
	return Local_CaptionGet_Common(
	    XL_captionId,
	    XL_getType,
        PersistentInXL );
}

__declspec(dllexport) LPXLOPER WINAPI Local_PXL_CaptionCalculator_SetData(
	LPXLOPER XL_captionId,
	LPXLOPER XL_dataId,
    LPXLOPER XL_setPortfolioType,
	LPXLOPER XL_MktDataKeys)
{
	ADD_LOG("Local_PXL_CaptionCalculator_SetData");
	bool PersistentInXL = false;
    bool isUpdated = false;
	return Local_CaptionSet_Common(
	    XL_captionId,
	    XL_dataId,
        XL_setPortfolioType,
		XL_MktDataKeys,
        isUpdated,
        PersistentInXL );
}

///////////////////////////////////
/// same as SetData but the previous Caption
/// is not cloned simply updated (=> no PXL version)
/// Only GC datas are updattable
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_CaptionCalculator_Update(
	LPXLOPER XL_captionId,
	LPXLOPER XL_dataId,
    LPXLOPER XL_setPortfolioType,
	LPXLOPER XL_mktDataKeys)
{
	ADD_LOG("Local_CaptionCalculator_Update");
	bool PersistentInXL = true;
    bool isUpdated = true;
	return Local_CaptionSet_Common(
	    XL_captionId,
	    XL_dataId,
        XL_setPortfolioType,
		XL_mktDataKeys,
        isUpdated,
        PersistentInXL );
}

__declspec(dllexport) LPXLOPER WINAPI Local_CaptionCalculator_GetPricingData(
	LPXLOPER XL_CaptionCalculatorId,
    LPXLOPER XL_KeyId )
{
	ADD_LOG("Local_CaptionCalculator_GetPricingData");
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
		XL_readStrCell( XL_CaptionCalculatorId, C_PricerString, " ARM_ERR: Caption Calculator Id: Object expected",C_result);
		long C_CaptionCalcId = LocalGetNumObjectId(C_PricerString);
		
		CCString C_KeyString;
		XL_readStrCell( XL_KeyId, C_KeyString, " ARM_ERR: Key Id: Object expected",C_result);
		
		/// call the function
		ARM_GramFctorArg argResult;
		
		long retCode = ARMLOCAL_GC_GetPricingData(
			C_CaptionCalcId,
			CCSTringToSTLString( C_KeyString ),
			argResult,
			C_result );
		
		if( retCode == ARM_KO )
		{
			ARM_ERR();
		}
		else
		{
			ARM_GramFunctorToXLOPER(argResult, XL_result, C_result);
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


/////////////////////////////////////////////////////////////
/// Creation of a PRDC calculator
/////////////////////////////////////////////////////////////
LPXLOPER Local_PRDCCalculator_Common(
    LPXLOPER XL_prdcId,
	LPXLOPER XL_modelId,
    LPXLOPER XL_otherMktDataIds,
    LPXLOPER XL_schedulerDatas,
    LPXLOPER XL_truncatorDatas,
    LPXLOPER XL_columnsToPrice,
    LPXLOPER XL_markovianDriftSamplerFlag,
    LPXLOPER XL_fxLocalModelFlag,
    LPXLOPER XL_calibType,
    LPXLOPER XL_calibDatas,
    LPXLOPER XL_basisIRCalibFlag,
	LPXLOPER XL_marginFlag,
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

        size_t i;

        /// PRDC Id
		CCString prdcStr;
		XL_readStrCell( XL_prdcId, prdcStr, " ARM_ERR: PRDC Calculator Id: Object expected",C_result);
		long prdcId = LocalGetNumObjectId(prdcStr);

        /// Model Id to price PRDC underlying
		CCString modelStr;
		XL_readStrCell( XL_modelId, modelStr, " ARM_ERR: Model Id: Object expected",C_result);
		long modelId = LocalGetNumObjectId(modelStr);

        /// OtherMktDatas
		vector<CCString> otherMktDataStr;
		XL_readStrVector(XL_otherMktDataIds,otherMktDataStr," ARM_ERR: OtherMktData Ids: array of objects expected",DOUBLE_TYPE,C_result);
		vector<long> otherMktDataIds(otherMktDataStr.size());
        for(i=0;i<otherMktDataIds.size();++i)
		    otherMktDataIds[i] = LocalGetNumObjectId(otherMktDataStr[i]);

        /// SchedulerDatas
		vector<double> schedulerDatas;
		vector<double> defaultSchedulerDatas(0);
		XL_readNumVectorWD(XL_schedulerDatas,schedulerDatas,defaultSchedulerDatas," ARM_ERR: SchedulerDatas: array of numeric expected",C_result);

        /// TruncatorDatas
		vector<double> truncatorDatas;
		vector<double> defaultTruncatorDatas(0);
		XL_readNumVectorWD(XL_truncatorDatas,truncatorDatas,defaultTruncatorDatas," ARM_ERR: TruncatorDatas: array of numeric expected",C_result);

		/// Flags to say what to price
		vector<CCString> columnsToPriceStr;
		vector<CCString> columnsToPriceDef(1,"PRDCOption");
		XL_readStrVectorWD(XL_columnsToPrice,columnsToPriceStr,columnsToPriceDef," ARM_ERR: columnsToPrice: array of string expected",DOUBLE_TYPE,C_result);
		vector< string > columnsToPrice(columnsToPriceStr.size());
		for(i=0;i<columnsToPriceStr.size();++i)
        {
            (columnsToPriceStr[i]).toUpper();
			columnsToPrice[i]=CCSTringToSTLString(columnsToPriceStr[i]);
        }

        vector< bool > prdcFlags(3);
		CCString markovianDriftSamplerStr;
		CCString markovianDriftSamplerDef("Y");
		XL_readStrCellWD( XL_markovianDriftSamplerFlag, markovianDriftSamplerStr,markovianDriftSamplerDef, " ARM_ERR: markovianDriftSamplerFlag: string expected",C_result);
        markovianDriftSamplerStr.toUpper();
        prdcFlags[0] = (CCSTringToSTLString(markovianDriftSamplerStr) == "Y");
            
		CCString fxLocalModelFlagStr;
		CCString fxLocalModelFlagDef("N");
		XL_readStrCellWD( XL_fxLocalModelFlag, fxLocalModelFlagStr,fxLocalModelFlagDef, " ARM_ERR: fxLocalModelFlag: string expected",C_result);
        fxLocalModelFlagStr.toUpper();
        prdcFlags[1] = (CCSTringToSTLString(fxLocalModelFlagStr) == "Y");

		CCString basisIRCalibFlagStr;
		CCString basisIRCalibFlagDef("Y");
		XL_readStrCellWD( XL_basisIRCalibFlag, basisIRCalibFlagStr,basisIRCalibFlagDef, " ARM_ERR: basisIRCalibFlag: string expected",C_result);
        basisIRCalibFlagStr.toUpper();
        prdcFlags[2] = (CCSTringToSTLString(basisIRCalibFlagStr) == "Y");

		CCString marginFlagStr;
		XL_readStrCellWD(XL_marginFlag, marginFlagStr,"FlowByFlow"," ARM_ERR: margin flag: string (Average/FlowByFlow)expected",C_result);
		marginFlagStr.toUpper();
		if(marginFlagStr == "EACHFLOW") marginFlagStr = "FLOWBYFLOW";
		else if(marginFlagStr == "START_END") marginFlagStr = "AVERAGE";
		string marginConvertType = CCSTringToSTLString(marginFlagStr);

		CCString calibTypeStr;
		CCString calibTypeDef("ATM");
		XL_readStrCellWD( XL_calibType, calibTypeStr,calibTypeDef, " ARM_ERR: calibType: string expected",C_result);
        calibTypeStr.toUpper();
        string calibType = CCSTringToSTLString(calibTypeStr);

		vector<double> calibDatas;
		vector<double> defaultCalibDatas(0);
		XL_readNumVectorWD(XL_calibDatas,calibDatas,defaultCalibDatas," ARM_ERR: calibDatas: array of numeric expected",C_result);

        exportFunc10Args< long, long, vector<long>, vector<double>, vector<double>, vector< string >, vector< bool>, string, vector<double>, string >
			ourFunc(prdcId,
			modelId,
			otherMktDataIds,
			schedulerDatas,
			truncatorDatas,
			columnsToPrice,
			prdcFlags,calibType,
			calibDatas,
			marginConvertType,
			ARMLOCAL_PRDCCalculator_Create);

		/// call the general function
		fillXL_Result_withName( ourFunc, C_result, XL_result, PersistentInXL );
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_PRDCCalculator_Common" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}
							 
__declspec(dllexport) LPXLOPER WINAPI Local_PRDCCalculator_Create(
    LPXLOPER XL_prdcId,
	LPXLOPER XL_modelId,
    LPXLOPER XL_otherMktDataIds,
    LPXLOPER XL_schedulerDatas,
    LPXLOPER XL_truncatorDatas,
    LPXLOPER XL_columnsToPrice,
    LPXLOPER XL_markovianDriftSamplerFlag,
    LPXLOPER XL_fxLocalModelFlag,
    LPXLOPER XL_fxCalibType,
    LPXLOPER XL_calibDatas,
    LPXLOPER XL_basisIRCalibFlag,
	LPXLOPER XL_marginFlag)
{
	ADD_LOG("Local_PRDCCalculator_Create");
	bool PersistentInXL = true;
	return Local_PRDCCalculator_Common(
            XL_prdcId,XL_modelId,
			XL_otherMktDataIds,
            XL_schedulerDatas,
			XL_truncatorDatas,
			XL_columnsToPrice,
            XL_markovianDriftSamplerFlag,
			XL_fxLocalModelFlag,
            XL_fxCalibType,
			XL_calibDatas,
			XL_basisIRCalibFlag,
			XL_marginFlag,
            PersistentInXL );
}

__declspec(dllexport) LPXLOPER WINAPI Local_PXL_PRDCCalculator_Create(
    LPXLOPER XL_prdcId,
	LPXLOPER XL_modelId,
    LPXLOPER XL_otherMktDataIds,
    LPXLOPER XL_schedulerDatas,
    LPXLOPER XL_truncatorDatas,
    LPXLOPER XL_columnsToPrice,
    LPXLOPER XL_markovianDriftSamplerFlag,
    LPXLOPER XL_fxLocalModelFlag,
    LPXLOPER XL_fxCalibType,
    LPXLOPER XL_calibDatas,
    LPXLOPER XL_basisIRCalibFlag,
	LPXLOPER XL_marginFlag)
{
	ADD_LOG("Local_PXL_PRDCCalculator_Create");
	bool PersistentInXL = false;
	return Local_PRDCCalculator_Common(
            XL_prdcId,XL_modelId,
			XL_otherMktDataIds,
            XL_schedulerDatas,
			XL_truncatorDatas,
			XL_columnsToPrice,
            XL_markovianDriftSamplerFlag,
			XL_fxLocalModelFlag,
            XL_fxCalibType,
			XL_calibDatas,
			XL_basisIRCalibFlag,
			XL_marginFlag,
            PersistentInXL );
}
/////////////////////////////////////////////////////////////
/// Creation of a PRDC calculator
/////////////////////////////////////////////////////////////
LPXLOPER Local_PRCSCalculator_Common(
			LPXLOPER XL_startDate,
			LPXLOPER XL_fixEndDate,
			LPXLOPER XL_endDate,
			LPXLOPER XL_Ccys,
			LPXLOPER XL_mktDataMger,
			LPXLOPER XL_couponStructure,
			LPXLOPER XL_fundingStructure,
			LPXLOPER XL_exerciseStructure,
			LPXLOPER XL_redemptionDatas,			
			LPXLOPER XL_Notionals,
			LPXLOPER XL_initialFx,
			LPXLOPER XL_Cpns,
			LPXLOPER XL_FundMargin,
			LPXLOPER XL_fees,
			LPXLOPER XL_CalibTypes,
			LPXLOPER XL_CalibDatas,
			LPXLOPER XL_Products,
			LPXLOPER XL_SchedulerDatas,
			LPXLOPER XL_TruncatorDatas,
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
		
		/// Required PRDC DATA						   
		
		/// Start & End
		double startDate;
		XL_readNumCell(XL_startDate,startDate," ARM_ERR: start date: date expected",C_result);
		double fixEndDate;
		XL_readNumCell(XL_fixEndDate,fixEndDate," ARM_ERR: Fix end date: date expected",C_result);
		double endDate;
		XL_readNumCell(XL_endDate,endDate," ARM_ERR: end date: date expected",C_result);

		VECTOR<CCString> Ccys;
		XL_readStrVector(XL_Ccys,Ccys," ARM_ERR: Currencies: array of string expected",DOUBLE_TYPE,C_result);
        if(Ccys.size() > 3)
			throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT," Currencies size should be less or equal than 3");
		string domCcyStr	= CCSTringToSTLString(Ccys[0]);
		string forCcyStr	= CCSTringToSTLString(Ccys[1]);
		string fundCcyStr  = CCSTringToSTLString(Ccys[2]);

		/// Market datas : zc curve , MRS, volatility curve and market models for OSW & CF pricings
		CCString mktDataMger;
		XL_readStrCell(XL_mktDataMger,mktDataMger," ARM_ERR: Mkt Data Manager: mkt Data Manager Id expected",C_result);
		long mktDataManagerId = LocalGetNumObjectId (mktDataMger);

		VECTOR<CCString> RequiredCpnDatas;
		XL_readStrVector(XL_couponStructure,RequiredCpnDatas," ARM_ERR: Required structure Datas : array of string expected",DOUBLE_TYPE,C_result);
        
		if(RequiredCpnDatas.size() > 7) 
			throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"Cpn Datas size != 7, the last data is resetTiming))");
		/// cpnFreq
		string cpnFreqStr = CCSTringToSTLString(RequiredCpnDatas[0]);
		///cpnDaycount
		string cpnDaycountStr = CCSTringToSTLString(RequiredCpnDatas[1]) ;
		/// fx reset lag
		long fxResetGagLg= (long)atof(RequiredCpnDatas[2]);
        ///pay and reset calendar
        string cpnResetCalStr = (RequiredCpnDatas.size() > 4) ? CCSTringToSTLString(RequiredCpnDatas[3]) : string ("");
        string cpnPayCalStr = (RequiredCpnDatas.size() > 5) ? CCSTringToSTLString(RequiredCpnDatas[4]) : string("");

		///stub Rule
		string stubRuleStr = CCSTringToSTLString((RequiredCpnDatas.size() > 5) ? RequiredCpnDatas[5] : CCString("SS"));	
		///reset timing
		string resetTimingStr = CCSTringToSTLString((RequiredCpnDatas.size() > 6) ? RequiredCpnDatas[6] : CCString("ARREARS"));	

		///end cpn structure
		VECTOR<CCString> RequiredFundDatas;
		XL_readStrVector(XL_fundingStructure,RequiredFundDatas," ARM_ERR: Required Funding Datas : array of string expected",DOUBLE_TYPE,C_result);
        
		if(RequiredFundDatas.size()!=2 )
			throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"Funding Datas size != 2))");

		///fundFreq
		string fundFreqStr = CCSTringToSTLString(RequiredFundDatas[0]);

		///fundDaycount
		string fundDaycountStr = CCSTringToSTLString(RequiredFundDatas[1]);

		///Redemption features
		VECTOR<CCString> redemptionDatas;
		XL_readStrVector(XL_redemptionDatas,redemptionDatas," ARM_ERR: Redemption Datas : array of string expected",DOUBLE_TYPE,C_result);
        
		if(redemptionDatas.size()!=3)
			throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"Redemption Datas size != 3))");
		/// redemption type
		string redemTypeStr = CCSTringToSTLString(redemptionDatas[0]);

		/// redemption gap
		long redemGapLg = (long)atof(redemptionDatas[1]);
		
		/// redemption strike
		double redemStrikeDble = (double)atof(redemptionDatas[2]);
		
		///Exercise features
		VECTOR<CCString> RequiredExerDatas;
		XL_readStrVector(XL_exerciseStructure,RequiredExerDatas," ARM_ERR: Required Exercise Datas : array of string expected",DOUBLE_TYPE,C_result);
        
		if(RequiredExerDatas.size()<3)
			throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"Exercise Datas size != 3))");

		/// exer Freq
		string exerFreqStr = CCSTringToSTLString(RequiredExerDatas[0]);	
		
		/// Notification Gap
		long exerNoticeGapLg = (long)atof(RequiredExerDatas[1]);

		/// payer or receiver
		string payRecStr = CCSTringToSTLString( RequiredExerDatas[2]);

		/// number of non call
		long nbNCallLg = (long)atof(RequiredExerDatas[3]);

	/// Market datas : zc curve , MRS, volatility curve and market models for OSW & CF pricings
		CCString initialFxStr;
		XL_readStrCell(XL_initialFx,initialFxStr," ARM_ERR: Initial Fx : curve Id expected",C_result);
		long initialFxId = LocalGetNumObjectId (initialFxStr);

		//cpnNotional		
		VECTOR<CCString> NotionalIds;
		XL_readStrVector(XL_Notionals,NotionalIds," ARM_ERR: nationals Id : array of Ids expected",DOUBLE_TYPE,C_result);
		if(NotionalIds.size() > 2)
			throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT," notional Curve Ids size should be less or equal than 2");
		if(NotionalIds.size() ==1 && domCcyStr != fundCcyStr)
			throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"  ARM_ERR: Notional Ids expected, PRCS needs two notionals");
		 
		long cpnNotionalId = LocalGetNumObjectId(NotionalIds[0]);
		long FundNotionalId = (NotionalIds.size() == 1 || domCcyStr == fundCcyStr) ? cpnNotionalId : LocalGetNumObjectId(NotionalIds[1]);

		//Cpn min andd Cpn max		
		VECTOR<CCString> CpnIds;
		XL_readStrVector(XL_Cpns,CpnIds," ARM_ERR: Cpn Ids : array of Ids expected",DOUBLE_TYPE,C_result);
		if(CpnIds.size() != 4)
			throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT," cpn Ids size should be equal 4");
		long domesticCpnId		= LocalGetNumObjectId(CpnIds[0]);
		long foreignCpnId		= LocalGetNumObjectId(CpnIds[1]);
		long minCpnId			= LocalGetNumObjectId(CpnIds[2]);
		long maxCpnId			= LocalGetNumObjectId(CpnIds[3]);

		VECTOR<CCString> FundCurveIds;
		XL_readStrVector(XL_FundMargin,FundCurveIds," ARM_ERR: Funding Curves Ids : array of Ids expected",DOUBLE_TYPE,C_result);
		if(FundCurveIds.size() > 2)
			throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT," funding curves ( margin and leverage) Ids size should be less than 2");
		long FundMargindId = LocalGetNumObjectId(FundCurveIds[0]);
		long FundLvgeId	= FundCurveIds.size() == 2 ? LocalGetNumObjectId(FundCurveIds[1]) : ARM_NULL_OBJECT;

		//Fees		
		CCString FeesStr;
		XL_readStrCell(XL_fees,FeesStr," ARM_ERR: Fees: Curve Id expected",C_result);
		long FeesId = LocalGetNumObjectId(FeesStr);

		///CalibTypes
		VECTOR<CCString> calibTypesCCStr;
		XL_readStrVector(XL_CalibTypes,calibTypesCCStr," ARM_ERR: Calib Types : array of string expected",DOUBLE_TYPE,C_result);

		vector< string > calibTypesStr(calibTypesCCStr.size());
		for(int i=0;i<calibTypesCCStr.size();++i)
			calibTypesStr[i]=CCSTringToSTLString(calibTypesCCStr[i]);

		///CalibDatas
		VECTOR<double> calibDatas;
		VECTOR<double> calibDatas_default;
		XL_readNumVectorWD(XL_CalibDatas,calibDatas,calibDatas_default," ARM_ERR: Calib Datas : array of numeric expected",C_result);
        
		/// ------------------------
		VECTOR<CCString> ProdXL;
		VECTOR<CCString> Prod_default(1,"PRDCOption");
        XL_readStrVectorWD(XL_Products,ProdXL,Prod_default," ARM_ERR: Product names: array of string expected",DOUBLE_TYPE,C_result);

		vector< string > ProdsStr;
		if ((ProdXL.size() != 1 ) || !(ProdXL[0] == "N"))
		{
			ProdsStr.resize(ProdXL.size());
			for(i=0;i<ProdXL.size();++i)
				ProdsStr[i]=CCSTringToSTLString(ProdXL[i]);
		}

		/// Model Datas 
		VECTOR<double> SchedulerDatasDble;
		VECTOR<double> Scheduler_default;
        XL_readNumVectorWD(XL_SchedulerDatas,SchedulerDatasDble,Scheduler_default," ARM_ERR: Scheduler Datas: array of numeric expected",C_result);
		
		VECTOR<double> TruncatorDatasDble;
		VECTOR<double> Truncator_default;
        XL_readNumVectorWD(XL_TruncatorDatas,TruncatorDatasDble,Truncator_default," ARM_ERR: Truncator Datas: array of numeric expected",C_result);

        exportFunc38Args< double,double,double,string,string,string,string,string,long,string,string,string,string,long,
			long,long,long,long,long,string,string,long,long,long,string,long,string,long,long, 
			string,long, double ,vector < string>, vector< double >,vector< string >,vector< double >, vector< double >,long >
		ourFunc(
				startDate,	
				fixEndDate, 
				endDate,				
				domCcyStr,	
				forCcyStr,		
				fundCcyStr,	
                cpnFreqStr, 
				cpnDaycountStr,
				fxResetGagLg,
				stubRuleStr,
				resetTimingStr,
                cpnResetCalStr,
                cpnPayCalStr,
				cpnNotionalId,
				domesticCpnId,
				foreignCpnId,
				initialFxId,
				minCpnId,
				maxCpnId,
				fundFreqStr,
				fundDaycountStr,
				FundNotionalId,
				FundMargindId,
				FundLvgeId,
				exerFreqStr,
				exerNoticeGapLg,
				payRecStr,
				nbNCallLg,                
				FeesId,
				redemTypeStr,
				redemGapLg,
				redemStrikeDble,
				calibTypesStr,
				calibDatas,
				ProdsStr,
				SchedulerDatasDble,
				TruncatorDatasDble,
				mktDataManagerId,
				ARMLOCAL_PRCSCalculator_Create);

	/// call the general function
		fillXL_Result_withName( ourFunc, C_result, XL_result, PersistentInXL );
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_BasisPRCSCalculator_Common" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;


}

///////////////////////////////////
/// version for non persistentInXl
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_PRCSCalculator_Create(
	LPXLOPER XL_startDate,
    LPXLOPER XL_fixEndDate,
	LPXLOPER XL_endDate,
	LPXLOPER XL_Ccys,
	LPXLOPER XL_mktDatas,
	LPXLOPER XL_couponStructure,
	LPXLOPER XL_fundingStructure,
	LPXLOPER XL_exerciseStructure,
	LPXLOPER XL_redemptionDatas,
	LPXLOPER XL_Notionals,
	LPXLOPER XL_initialFx,
	LPXLOPER XL_Cpns,	
	LPXLOPER XL_FundMargin,
	LPXLOPER XL_fees,
	LPXLOPER XL_CalibTypes,
	LPXLOPER XL_CalibDatas,
	LPXLOPER XL_Products,
	LPXLOPER XL_SchedulerDatas,
	LPXLOPER XL_TruncatorDatas)
{
	ADD_LOG("Local_PRCSCalculator_Create");
	bool PersistentInXL = true;

	return Local_PRCSCalculator_Common(
			XL_startDate,
			XL_fixEndDate,
			XL_endDate,
			XL_Ccys,
			XL_mktDatas,
			XL_couponStructure,
			XL_fundingStructure,
			XL_exerciseStructure,
			XL_redemptionDatas,
			XL_Notionals,
			XL_initialFx,
			XL_Cpns,
			XL_FundMargin,
			XL_fees,
			XL_CalibTypes,
			XL_CalibDatas,
			XL_Products,
			XL_SchedulerDatas,
			XL_TruncatorDatas,
			PersistentInXL );
}

///////////////////////////////////
/// version that takes into account 
/// previous creation of object
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_PXL_PRCSCalculator_Create(
	LPXLOPER XL_startDate,
    LPXLOPER XL_fixEndDate,
	LPXLOPER XL_endDate,
	LPXLOPER XL_Ccys,
	LPXLOPER XL_mktDatas,
	LPXLOPER XL_couponStructure,
	LPXLOPER XL_fundingStructure,
	LPXLOPER XL_exerciseStructure,
	LPXLOPER XL_redemptionDatas,
	LPXLOPER XL_Notionals,
	LPXLOPER XL_initialFx,
	LPXLOPER XL_Cpns,	
	LPXLOPER XL_FundMargin,
	LPXLOPER XL_fees,
	LPXLOPER XL_CalibTypes,
	LPXLOPER XL_CalibDatas,
	LPXLOPER XL_Products,
	LPXLOPER XL_SchedulerDatas,
	LPXLOPER XL_TruncatorDatas)
{
	ADD_LOG("Local_PXL_PRCSCalculator_Create");
	bool PersistentInXL = false;

	return Local_PRCSCalculator_Common(
			XL_startDate,
			XL_fixEndDate,
			XL_endDate,
			XL_Ccys,
			XL_mktDatas,
			XL_couponStructure,
			XL_fundingStructure,
			XL_exerciseStructure,
			XL_redemptionDatas,
			XL_Notionals,
			XL_initialFx,
			XL_Cpns,
			XL_FundMargin,
			XL_fees,
			XL_CalibTypes,
			XL_CalibDatas,
			XL_Products,
			XL_SchedulerDatas,
			XL_TruncatorDatas,
			PersistentInXL );
}
/////////////////////////////////////////////////////////////
/// Creation of a PRDC calculator
/////////////////////////////////////////////////////////////
LPXLOPER Local_PRDKOCalculator_Common(
			LPXLOPER XL_Dates,
			LPXLOPER XL_Ccys,
			LPXLOPER XL_mktDataMger,
			LPXLOPER XL_couponStructure,
			LPXLOPER XL_fundingStructure,
			LPXLOPER XL_exerciseStructure,
			LPXLOPER XL_redemptionDatas,			
			LPXLOPER XL_Notionals,
			LPXLOPER XL_CpnCurves,
			LPXLOPER XL_FundMargin,
			LPXLOPER XL_fees,
			LPXLOPER XL_CalibTypes,
			LPXLOPER XL_CalibDatas,
			LPXLOPER XL_Products,
			LPXLOPER XL_SchedulerDatas,
			LPXLOPER XL_TruncatorDatas,
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
		
		/// Required PRDC DATA						   
		
		/// Start, End, end of fix period & switch dates
		VECTOR<double> vdates;
		XL_readNumVector(XL_Dates,vdates," ARM_ERR: Dates: vector of numbers expected",C_result);
		if(vdates.size() < 2 || vdates.size() > 4)
			throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"vector of dates need at least 2 dates ( start,end) and at most 4 dates(..., fixEnd and switch)");
		double startDate = vdates[0];
		double endDate   = vdates[1];
		double fixEndDate = startDate;
		double switchDate = startDate; 
		if(vdates.size() > 2) fixEndDate= vdates[2];
		if(vdates.size() == 4) switchDate= vdates[3];

		VECTOR<CCString> Ccys;
		XL_readStrVector(XL_Ccys,Ccys," ARM_ERR: Currencies: array of string expected",DOUBLE_TYPE,C_result);
        if(Ccys.size() > 3)
			throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT," Currencies size should be less or equal than 3");
		string domCcyStr	= CCSTringToSTLString(Ccys[0]);
		string forCcyStr	= CCSTringToSTLString(Ccys[1]);
		string fundCcyStr  = CCSTringToSTLString(Ccys[2]);

		/// Market datas : zc curve , MRS, volatility curve and market models for OSW & CF pricings
		CCString mktDataMger;
		XL_readStrCell(XL_mktDataMger,mktDataMger," ARM_ERR: Mkt Data Manager: mkt Data Manager Id expected",C_result);
		long mktDataManagerId = LocalGetNumObjectId (mktDataMger);

		VECTOR<CCString> RequiredCpnDatas;
		XL_readStrVector(XL_couponStructure,RequiredCpnDatas," ARM_ERR: Required structure Datas : array of string expected",DOUBLE_TYPE,C_result);
        
		if(RequiredCpnDatas.size() > 7) 
			throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"Cpn Datas size != 7, the last data is resetTiming))");
		/// cpnFreq
		string cpnFreqStr = CCSTringToSTLString(RequiredCpnDatas[0]);
		///cpnDaycount
		string cpnDaycountStr = CCSTringToSTLString(RequiredCpnDatas[1]) ;
		/// fx reset lag
		long fxResetGagLg= (long)atof(RequiredCpnDatas[2]);
        ///pay and reset calendar
        string cpnResetCalStr = (RequiredCpnDatas.size() > 4) ? CCSTringToSTLString(RequiredCpnDatas[3]) : string ("");
        string cpnPayCalStr = (RequiredCpnDatas.size() > 5) ? CCSTringToSTLString(RequiredCpnDatas[4]) : string("");

		///stub Rule
		string stubRuleStr = CCSTringToSTLString((RequiredCpnDatas.size() > 5) ? RequiredCpnDatas[5] : CCString("SS"));	
		///reset timing
		string resetTimingStr = CCSTringToSTLString((RequiredCpnDatas.size() > 6) ? RequiredCpnDatas[6] : CCString("ARREARS"));	

		///end cpn structure
		VECTOR<CCString> RequiredFundDatas;
		XL_readStrVector(XL_fundingStructure,RequiredFundDatas," ARM_ERR: Required Funding Datas : array of string expected",DOUBLE_TYPE,C_result);
        
		if(RequiredFundDatas.size()!=2 )
			throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"Funding Datas size != 2))");

		///fundFreq
		string fundFreqStr = CCSTringToSTLString(RequiredFundDatas[0]);

		///fundDaycount
		string fundDaycountStr = CCSTringToSTLString(RequiredFundDatas[1]);

		///Redemption features
		VECTOR<CCString> redemptionDatas;
		XL_readStrVector(XL_redemptionDatas,redemptionDatas," ARM_ERR: Redemption Datas : array of string expected",DOUBLE_TYPE,C_result);
        
		if(redemptionDatas.size()!=3)
			throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"Redemption Datas size != 3))");
		/// redemption type
		string redemTypeStr = CCSTringToSTLString(redemptionDatas[0]);

		/// redemption gap
		long redemGapLg = (long)atof(redemptionDatas[1]);
		
		/// redemption strike
		double redemStrikeDble = (double)atof(redemptionDatas[2]);
		
		///Exercise features
		VECTOR<CCString> RequiredExerDatas;
		XL_readStrVector(XL_exerciseStructure,RequiredExerDatas," ARM_ERR: Required Exercise Datas : array of string expected",DOUBLE_TYPE,C_result);
        
		if(RequiredExerDatas.size()<3)
			throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"Exercise Datas size != 3))");

		/// exer Freq
		string exerFreqStr = CCSTringToSTLString(RequiredExerDatas[0]);	
		
		/// Notification Gap
		long exerNoticeGapLg = (long)atof(RequiredExerDatas[1]);

		/// payer or receiver
		string payRecStr = CCSTringToSTLString( RequiredExerDatas[2]);

		/// number of non call
		long nbNCallLg = (long)atof(RequiredExerDatas[3]);

		//cpnNotional		
		VECTOR<CCString> NotionalIds;
		XL_readStrVector(XL_Notionals,NotionalIds," ARM_ERR: nationals Id : array of Ids expected",DOUBLE_TYPE,C_result);
		if(NotionalIds.size() > 2)
			throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT," notional Curve Ids size should be less or equal than 2");
		if(NotionalIds.size() ==1 && domCcyStr != fundCcyStr)
			throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"  ARM_ERR: Notional Ids expected, PRCS needs two notionals");
		 
		long cpnNotionalId = LocalGetNumObjectId(NotionalIds[0]);
		long FundNotionalId = (NotionalIds.size() == 1 || domCcyStr == fundCcyStr) ? cpnNotionalId : LocalGetNumObjectId(NotionalIds[1]);

		//Cpn min andd Cpn max		
		VECTOR<CCString> CpnIds;
		XL_readStrVector(XL_CpnCurves,CpnIds," ARM_ERR: Cpn Ids : array of Ids expected",DOUBLE_TYPE,C_result);
		if(CpnIds.size() != 6)
			throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT," cpn Ids size should be equal 6 ( initial, domestic Cpn, foriegn Cpn, min Cpn , max Cpn and barrier)");
		long initialFxId		= LocalGetNumObjectId(CpnIds[0]);
		long domesticCpnId		= LocalGetNumObjectId(CpnIds[1]);
		long foreignCpnId		= LocalGetNumObjectId(CpnIds[2]);
		long minCpnId			= LocalGetNumObjectId(CpnIds[3]);
		long maxCpnId			= LocalGetNumObjectId(CpnIds[4]);
		long barrierId			= LocalGetNumObjectId(CpnIds[5]);

		VECTOR<CCString> FundCurveIds;
		XL_readStrVector(XL_FundMargin,FundCurveIds," ARM_ERR: Funding Curves Ids : array of Ids expected",DOUBLE_TYPE,C_result);
		if(FundCurveIds.size() > 2)
			throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT," funding curves ( margin and leverage) Ids size should be less than 2");
		long FundMargindId = LocalGetNumObjectId(FundCurveIds[0]);
		long FundLvgeId	= FundCurveIds.size() == 2 ? LocalGetNumObjectId(FundCurveIds[1]) : ARM_NULL_OBJECT;

		//Fees		
		CCString FeesStr;
		XL_readStrCell(XL_fees,FeesStr," ARM_ERR: Fees: Curve Id expected",C_result);
		long FeesId = LocalGetNumObjectId(FeesStr);

		///CalibTypes
		VECTOR<CCString> calibTypesCCStr;
		XL_readStrVector(XL_CalibTypes,calibTypesCCStr," ARM_ERR: Calib Types : array of string expected",DOUBLE_TYPE,C_result);

		vector< string > calibTypesStr(calibTypesCCStr.size());
		for(int i=0;i<calibTypesCCStr.size();++i)
			calibTypesStr[i]=CCSTringToSTLString(calibTypesCCStr[i]);

		///CalibDatas
		VECTOR<double> calibDatas;
		VECTOR<double> calibDatas_default;
		XL_readNumVectorWD(XL_CalibDatas,calibDatas,calibDatas_default," ARM_ERR: Calib Datas : array of numeric expected",C_result);
        
		/// ------------------------
		VECTOR<CCString> ProdXL;
		VECTOR<CCString> Prod_default(1,"PRDCOption");
        XL_readStrVectorWD(XL_Products,ProdXL,Prod_default," ARM_ERR: Product names: array of string expected",DOUBLE_TYPE,C_result);

		vector< string > ProdsStr;
		if ((ProdXL.size() != 1 ) || !(ProdXL[0] == "N"))
		{
			ProdsStr.resize(ProdXL.size());
			for(i=0;i<ProdXL.size();++i)
				ProdsStr[i]=CCSTringToSTLString(ProdXL[i]);
		}

		/// Model Datas 
		VECTOR<double> SchedulerDatasDble;
		VECTOR<double> Scheduler_default;
        XL_readNumVectorWD(XL_SchedulerDatas,SchedulerDatasDble,Scheduler_default," ARM_ERR: Scheduler Datas: array of numeric expected",C_result);
		
		VECTOR<double> TruncatorDatasDble;
		VECTOR<double> Truncator_default;
        XL_readNumVectorWD(XL_TruncatorDatas,TruncatorDatasDble,Truncator_default," ARM_ERR: Truncator Datas: array of numeric expected",C_result);

        exportFunc40Args< double,double,double,double,string,string,string,string,string,long,string,string,string,string,long,
			long,long,long,long,long,long,string,string,long,long,long,string,long,string,long,long, 
			string,long, double ,vector < string>, vector< double >,vector< string >,vector< double >, vector< double >,long >
		ourFunc(
				startDate,	
				fixEndDate,
				switchDate,
				endDate,				
				domCcyStr,	
				forCcyStr,		
				fundCcyStr,	
                cpnFreqStr, 
				cpnDaycountStr,
				fxResetGagLg,
				stubRuleStr,
				resetTimingStr,
                cpnResetCalStr,
                cpnPayCalStr,
				cpnNotionalId,
				domesticCpnId,
				foreignCpnId,
				initialFxId,
				minCpnId,
				maxCpnId,
				barrierId,
				fundFreqStr,
				fundDaycountStr,
				FundNotionalId,
				FundMargindId,
				FundLvgeId,
				exerFreqStr,
				exerNoticeGapLg,
				payRecStr,
				nbNCallLg,                
				FeesId,
				redemTypeStr,
				redemGapLg,
				redemStrikeDble,
				calibTypesStr,
				calibDatas,
				ProdsStr,
				SchedulerDatasDble,
				TruncatorDatasDble,
				mktDataManagerId,
				ARMLOCAL_PRDKOCalculator_Create);

	/// call the general function
		fillXL_Result_withName( ourFunc, C_result, XL_result, PersistentInXL );
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_BasisPRCSCalculator_Common" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;


}

///////////////////////////////////
//// Power Reverse Dual Currencies Knock out
/// version for non persistentInXl
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_PRDKOCalculator_Create(
	LPXLOPER XL_Dates,
	LPXLOPER XL_Ccys,
	LPXLOPER XL_mktDatas,
	LPXLOPER XL_couponStructure,
	LPXLOPER XL_fundingStructure,
	LPXLOPER XL_exerciseStructure,
	LPXLOPER XL_redemptionDatas,
	LPXLOPER XL_Notionals,
	LPXLOPER XL_CpnCurves,	
	LPXLOPER XL_FundMargin,
	LPXLOPER XL_fees,
	LPXLOPER XL_CalibTypes,
	LPXLOPER XL_CalibDatas,
	LPXLOPER XL_Products,
	LPXLOPER XL_SchedulerDatas,
	LPXLOPER XL_TruncatorDatas)
{
	ADD_LOG("Local_PRDKOCalculator_Create");
	bool PersistentInXL = true;

	return Local_PRDKOCalculator_Common(
			XL_Dates,
			XL_Ccys,
			XL_mktDatas,
			XL_couponStructure,
			XL_fundingStructure,
			XL_exerciseStructure,
			XL_redemptionDatas,
			XL_Notionals,
			XL_CpnCurves,
			XL_FundMargin,
			XL_fees,
			XL_CalibTypes,
			XL_CalibDatas,
			XL_Products,
			XL_SchedulerDatas,
			XL_TruncatorDatas,
			PersistentInXL );
}

///////////////////////////////////
/// version that takes into account 
/// previous creation of object
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_PXL_PRDKOCalculator_Create(
	LPXLOPER XL_Dates,
	LPXLOPER XL_Ccys,
	LPXLOPER XL_mktDatas,
	LPXLOPER XL_couponStructure,
	LPXLOPER XL_fundingStructure,
	LPXLOPER XL_exerciseStructure,
	LPXLOPER XL_redemptionDatas,
	LPXLOPER XL_Notionals,
	LPXLOPER XL_CpnCurves,	
	LPXLOPER XL_FundMargin,
	LPXLOPER XL_fees,
	LPXLOPER XL_CalibTypes,
	LPXLOPER XL_CalibDatas,
	LPXLOPER XL_Products,
	LPXLOPER XL_SchedulerDatas,
	LPXLOPER XL_TruncatorDatas)
{
	ADD_LOG("Local_PXL_PRDKOCalculator_Create");
	bool PersistentInXL = false;

	return Local_PRDKOCalculator_Common(
			XL_Dates,
			XL_Ccys,
			XL_mktDatas,
			XL_couponStructure,
			XL_fundingStructure,
			XL_exerciseStructure,
			XL_redemptionDatas,
			XL_Notionals,
			XL_CpnCurves,
			XL_FundMargin,
			XL_fees,
			XL_CalibTypes,
			XL_CalibDatas,
			XL_Products,
			XL_SchedulerDatas,
			XL_TruncatorDatas,
			PersistentInXL );
}



/////////////////////////////////////////////////////////////
/// Get datas from a PRDC calculator
/////////////////////////////////////////////////////////////
LPXLOPER Local_PRDCGet_Common(
	LPXLOPER XL_prdcId,
	LPXLOPER XL_getType,
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

		long prdcId;
		XL_GETOBJID( XL_prdcId,prdcId," ARM_ERR: PRDC Calculator: Object expected",C_result);

		CCString getTypeStr;
		XL_readStrCell(XL_getType,getTypeStr," ARM_ERR: Accessor Type: string expected",C_result);
        string getType(CCSTringToSTLString(getTypeStr));

		CCString prdcGetClass(GCGetTypeToClass(getType,prdcId).c_str());

        exportFunc2Args< long, string > ourFunc(prdcId,getType, ARMLOCAL_PRDC_Get);

		/// call the general function
		fillXL_Result( prdcGetClass, ourFunc, C_result, XL_result, PersistentInXL );
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_PRDCGet_Common" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}

__declspec(dllexport) LPXLOPER WINAPI Local_PRDCCalculator_GetData(
    LPXLOPER XL_prdcId,
    LPXLOPER XL_getType)
{
	ADD_LOG("Local_PRDCCalculator_GetData");
	bool PersistentInXL = true;
	return Local_PRDCGet_Common(XL_prdcId,XL_getType,PersistentInXL );
}

__declspec(dllexport) LPXLOPER WINAPI Local_PXL_PRDCCalculator_GetData(
    LPXLOPER XL_prdcId,
    LPXLOPER XL_getType)
{
	ADD_LOG("Local_PXL_PRDCCalculator_GetData");
	bool PersistentInXL = false;
	return Local_PRDCGet_Common(XL_prdcId,XL_getType,PersistentInXL );
}





///////////////////////////////////////////////////////////////////////////////////////////////
////////////////Callable SnowBall Calculator///////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////
///----------------------------------------------
///----------------------------------------------
///         Callable SnowBall Calculator
/// Inputs :
///     
///----------------------------------------------
///----------------------------------------------
class callableSBCalculatorFunc : public ARMResultLong2LongFunc
{
public:
	callableSBCalculatorFunc(
        double startDate,
        double endDate,
		long cpnFreq,
		long cpnDaycount,
		long cpnIndexTiming,
		const string& cpnIdxTerm,
		long cpnIndexDaycount,
		const string& cpnResetCal,
		const string& cpnPayCal,
		long cpnIntRule,
		long cpnIndexResetLag,
		long payRec,
		long CF,
		double notional,
		long notionalId,
		double fundnotional,
		long fundnotionalId,
		double dConst,
		long constId,
		double lPrevCpn,
		long lPrevCpnId,
		double lNewOpt,
		long lNewOptId,
		double strikeOpt,
		long strikeOptId,
		double minCpn,
		long minCpnId,
		double maxCpn,
		long maxCpnId,
		long fundFreq,
		long fundDaycount,
		double fundCoeff,
		double fundCoeffId,
		double fundMargin,
		long fundMarginId,
		long NotifDays,
		long NonCall,
		const string& exerciseCal,
		bool callSBOrStacky,
		long feesId,
		long NbPathBounding,
		long NbPathPricing,
		long NbMaxBucket,
		long FixSteps,
		const string& USMethod,
		const vector< string >& CalibMod,
		const vector< string >& flags,
		long mktDataManagerId,
        const vector< string >& keys,
		const string& GenType1,
		const string& GenType2,
		int FirstNbTimes,
		int FirstNbDims,
		const string& PathScheme,
		const string& PathOrder,
		const string& Antithetic,
		const string& ModelType,
		const string& TriggerMode,
		int calibSwoptFreq,
		const string& regressors,
		const vector< double >& hkDatas)
    :
	C_startDate(startDate),
    C_endDate(endDate),
	C_cpnFreq(cpnFreq),
	C_cpnDaycount(cpnDaycount),
	C_cpnIndexTiming(cpnIndexTiming),
	C_cpnIdxTerm(cpnIdxTerm),
	C_cpnIndexDaycount(cpnIndexDaycount),
	C_cpnResetCal(cpnResetCal),
	C_cpnPayCal(cpnPayCal),
	C_cpnIntRule(cpnIntRule),
	C_cpnIndexResetLag(cpnIndexResetLag),
	C_payRec(payRec),
	C_CF(CF),
	C_notional(notional),
	C_notionalId(notionalId),
	C_Fundnotional(fundnotional),
	C_FundnotionalId(fundnotionalId),
	C_dConst(dConst),
	C_constId(constId),
	C_lPrevCpn(lPrevCpn),
	C_lPrevCpnId(lPrevCpnId),
	C_lNewOpt(lNewOpt),
	C_lNewOptId(lNewOptId),
	C_strikeOpt(strikeOpt),
	C_strikeOptId(strikeOptId),
	C_minCpn(minCpn),
	C_minCpnId(minCpnId),
	C_maxCpn(maxCpn),
	C_maxCpnId(maxCpnId),
	C_fundFreq(fundFreq),
	C_fundDaycount(fundDaycount),
	C_fundCoeff(fundCoeff),
	C_fundCoeffId(fundCoeffId),
	C_fundMargin(fundMargin),
	C_fundMarginId(fundMarginId),
	C_NotifDays(NotifDays),
	C_NonCall(NonCall),
	C_exerciseCal(exerciseCal),
	C_CallSBOrStacky(callSBOrStacky),
	C_feesId(feesId),
	C_NbPathBounding(NbPathBounding),
	C_NbPathPricing(NbPathPricing),
	C_NbMaxBucket(NbMaxBucket),
	C_FixSteps(FixSteps),
	C_USMethod(USMethod),
	C_CalibMod(CalibMod),
	C_flags(flags),
	C_mktDataManagerId(mktDataManagerId),
    C_keys(keys),
	C_GenType1(GenType1),
	C_GenType2(GenType2),
	C_FirstNbTimes(FirstNbTimes),
	C_FirstNbDims(FirstNbDims),
	C_PathScheme(PathScheme),
	C_PathOrder(PathOrder),
	C_Antithetic(Antithetic),
	C_ModelType(ModelType),
	C_TriggerMode(TriggerMode),
	C_CalibSwoptFreq(calibSwoptFreq),
	C_Regressors(regressors),
	C_hkDatas(hkDatas)
    {};

	long operator()( ARM_result& result, long objId )
	{
		return ARMLOCAL_CallableSBCalculator_Create(
			C_startDate,
			C_endDate,
			C_cpnFreq,
			C_cpnDaycount,
			C_cpnIndexTiming,
			C_cpnIdxTerm,
			C_cpnIndexDaycount,
			C_cpnResetCal,
			C_cpnPayCal,
			C_cpnIntRule,
			C_cpnIndexResetLag,
			C_payRec,
			C_CF,
			C_notional,
			C_notionalId,
			C_Fundnotional,
			C_FundnotionalId,
			C_dConst,
			C_constId,
			C_lPrevCpn,
			C_lPrevCpnId,
			C_lNewOpt,
			C_lNewOptId,
			C_strikeOpt,
			C_strikeOptId,
			C_minCpn,
			C_minCpnId,
			C_maxCpn,
			C_maxCpnId,
			C_fundFreq,
			C_fundDaycount,
			C_fundCoeff,
			C_fundCoeffId,
			C_fundMargin,
			C_fundMarginId,
			C_NotifDays,
			C_NonCall,
			C_exerciseCal,
			C_CallSBOrStacky,
			C_feesId,
			C_NbPathBounding,
			C_NbPathPricing,
			C_NbMaxBucket,
			C_FixSteps,
			C_USMethod,
			C_CalibMod,
			C_flags,
			C_mktDataManagerId,
			C_keys,
			C_GenType1,
			C_GenType2,
			C_FirstNbTimes,
			C_FirstNbDims,
			C_PathScheme,
			C_PathOrder,
			C_Antithetic,
			C_ModelType,
			C_TriggerMode,
			C_CalibSwoptFreq,
			C_Regressors,
			C_hkDatas,
            result,
            objId);
    }

private:

	double				C_startDate;
	double				C_endDate;
	long				C_cpnFreq;
	long				C_cpnDaycount;
	long				C_cpnIndexTiming;
	string				C_cpnIdxTerm;
	long				C_cpnIndexDaycount;
	string				C_cpnResetCal;
	string				C_cpnPayCal;
	long				C_cpnIntRule;
	long				C_cpnIndexResetLag;
	long				C_payRec;
	long				C_CF;
	double				C_notional;
	long				C_notionalId;
	double				C_Fundnotional;
	long				C_FundnotionalId;
	double				C_dConst;
	long				C_constId;
	double				C_lPrevCpn;
	long				C_lPrevCpnId;
	double				C_lNewOpt;
	long				C_lNewOptId;
	double				C_strikeOpt;
	long				C_strikeOptId;
	double				C_minCpn;
	long				C_minCpnId;
	double				C_maxCpn;
	long				C_maxCpnId;
	long				C_fundFreq;
	long				C_fundDaycount;
	double				C_fundCoeff;
	long				C_fundCoeffId;
	double				C_fundMargin;
	long				C_fundMarginId;
	long				C_NotifDays;
	long				C_NonCall;
	bool				C_CallSBOrStacky;
	string				C_exerciseCal;
	long				C_feesId;
	long				C_NbPathBounding;
	long				C_NbPathPricing;
	long				C_NbMaxBucket;
	long				C_FixSteps;
	string				C_USMethod;
	vector< string >	C_CalibMod;

	vector< string >	C_flags;
	long				C_mktDataManagerId;
	vector< string >	C_keys;
	string				C_GenType1;
	string				C_GenType2;
	int					C_FirstNbTimes;
	int					C_FirstNbDims;
	string				C_PathScheme;
	string				C_PathOrder;
	string				C_Antithetic;
	string				C_ModelType;
	string				C_TriggerMode;
	int					C_CalibSwoptFreq;
	string				C_Regressors;
	vector< double >	C_hkDatas;
};



LPXLOPER Local_CallableSBCalculator_Common(
    LPXLOPER XL_startDate,
    LPXLOPER XL_endDate,
	LPXLOPER XL_mktDatas,
	LPXLOPER XL_payOrRec,
	LPXLOPER XL_coupon,
	LPXLOPER XL_Notional,
	LPXLOPER XL_CpnCurves,
	LPXLOPER XL_CapOrFloor,
	LPXLOPER XL_fundDatas,
	LPXLOPER XL_ExerciseData,
	LPXLOPER XL_fees,
	LPXLOPER XL_CalibFlags,
	LPXLOPER XL_ProductsFlags,
	LPXLOPER XL_ModelDatas,
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

		///-------------------------------------------------------------------
		///                     Required Callable SnowBall DATA			   ///
		/// ------------------------------------------------------------------
    
		/// General datas
		/// -----------------
		
		/// Start & End
		double startDate;
		XL_readNumCell(XL_startDate,startDate," ARM_ERR: start date: date expected",C_result);
		double endDate;
		XL_readNumCell(XL_endDate,endDate," ARM_ERR: end date: date expected",C_result);

		/// Market datas : zc curve , MRS, beta, volatility curve and market models for OSW & CF pricings
		VECTOR<CCString> mktDatas;
		XL_readStrVector(XL_mktDatas,mktDatas," ARM_ERR: Market datas: array of string expected",DOUBLE_TYPE,C_result);
        if(mktDatas.size()<1)
		{
			throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"Mkt Data size should be greater than 1");
		}

		long mktDataManagerId = LocalGetNumObjectId (mktDatas[0]);
		vector< string > keys(mktDatas.size()-1);
        size_t i;
		for(i=1;i<mktDatas.size();++i)
			keys[i-1]=CCSTringToSTLString(mktDatas[i]);

		/// Pay/Rec
		CCString payRecStr;
		XL_readStrCell(XL_payOrRec,payRecStr," ARM_ERR: pay or receive: string expected",C_result);
		long payRec;
		if((payRec = ARM_ConvRecOrPay (payRecStr, C_result)) == ARM_DEFAULT_ERR)
		{
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}

		VECTOR<CCString> RequiredCpnDatas;
		XL_readStrVector(XL_coupon,RequiredCpnDatas," ARM_ERR: Required Cpn Datas : array of string expected",DOUBLE_TYPE,C_result);
        
		if(RequiredCpnDatas.size()!=9)
		{
			throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"Cpn Datas size != ))");
		}

		/// cpnFreq
		CCString cpnFreqStr = RequiredCpnDatas[0];
		long cpnFreq;
		if((cpnFreq = ARM_ConvFrequency (cpnFreqStr, C_result)) == ARM_DEFAULT_ERR)
		{
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}

		/// cpnDaycount
		CCString cpnDaycountStr = RequiredCpnDatas[1];
		long cpnDaycount = ARM_ConvDayCount (cpnDaycountStr);

		/// cpnIndexTiming
		CCString cpnIndexTimingStr = RequiredCpnDatas[2];
		long cpnIndexTiming = ARM_ConvPayResetRule (cpnIndexTimingStr);

		/// cpnIndexTerm
		CCString cpnIdxTermStr = RequiredCpnDatas[3];
		string cpnIdxTerm=CCSTringToSTLString(cpnIdxTermStr);

		/// cpnIndexDaycount
		CCString cpnIndexDaycountStr = RequiredCpnDatas[4];
		long cpnIndexDaycount = ARM_ConvDayCount (cpnIndexDaycountStr);

		/// Reset Calendar
		CCString cpnResetCalStr = RequiredCpnDatas[5];
		string cpnResetCal = CCSTringToSTLString(cpnResetCalStr);

		/// pay Calendar
		CCString cpnPayCalStr = RequiredCpnDatas[6];
		string cpnPayCal = CCSTringToSTLString(cpnPayCalStr);

		/// cpnIntRule
		CCString cpnIntRuleStr = RequiredCpnDatas[7];
		long cpnIntRule = ARM_ConvIntRule (cpnIntRuleStr);

		/// cpnIntRule
		double cpnIndexResetLag = atof(RequiredCpnDatas[8]);

		// end cpn structure
		

		// nominals
		VECTOR<CCString> nominals;
		XL_readStrVector(XL_Notional,nominals, " ARM_ERR: Nominals: array of string or numeric expected",DOUBLE_TYPE,C_result);

		long     nominalId;
		double nominal;

		CCString firstCc("L");
		if(nominals.size() >= 1)
		{
			if(nominals[0].Contain(firstCc))
				nominalId = LocalGetNumObjectId(nominals[0]);
			else
			{
				nominalId = ARM_NULL_OBJECT;
				nominal = atof(nominals[0].c_str());
			}
		}
		else
		{
			nominalId = ARM_NULL_OBJECT;
			nominal = 100.0;
		}

		long     fundNominalId;
		double fundNominal;

		if(nominals.size() >= 2)
		{
			if(nominals[1].Contain(firstCc))
				fundNominalId = LocalGetNumObjectId(nominals[1]);
			else
			{
				fundNominalId = ARM_NULL_OBJECT;
				fundNominal = atof(nominals[1].c_str());
			}
		}
		else
		{
			fundNominalId = nominalId;
			fundNominal = nominal;
		}

		long nbRows, nbCols;
		VECTOR<CCString> CpnCurves;
		VECTOR<long> CpnCurvesType;
		XL_readStrVectorSizeAndType(XL_CpnCurves,nbRows,nbCols,CpnCurves,CpnCurvesType," ARM_ERR: Required Cpn Curves : array of string or double expected",DOUBLE_TYPE,C_result);
        
		if(CpnCurves.size()!=6)
		{
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}

		//Const
		double dConst;
		long constId;
		if(CpnCurvesType[0] == XL_TYPE_STRING)
			constId = LocalGetNumObjectId(CpnCurves[0]);
		else
		{
			dConst = atof(CpnCurves[0].c_str());
			constId = ARM_NULL_OBJECT;
		}

		//LPrevCpn
		double lPrevCpn;
		long lPrevCpnId;
		if(CpnCurvesType[1] == XL_TYPE_STRING)
			lPrevCpnId = LocalGetNumObjectId(CpnCurves[1]);
		else
		{
			lPrevCpn = atof(CpnCurves[1].c_str());
			lPrevCpnId = ARM_NULL_OBJECT;
		}

		//LNewOpt
		double lNewOpt;
		long     lNewOptId;
		if(CpnCurvesType[2] == XL_TYPE_STRING)
			lNewOptId = LocalGetNumObjectId(CpnCurves[2]);
		else
		{
			lNewOpt = atof(CpnCurves[2].c_str());
			lNewOptId = ARM_NULL_OBJECT;
		}

		//strikeOpt
		double strikeOpt;
		long     strikeOptId;
		if(CpnCurvesType[3] == XL_TYPE_STRING)
			strikeOptId = LocalGetNumObjectId(CpnCurves[3]);
		else
		{
			strikeOpt = atof(CpnCurves[3].c_str());
			strikeOptId = ARM_NULL_OBJECT;
		}

		//min Cpn
		double minCpn;
		long     minCpnId;
		if(CpnCurvesType[4] == XL_TYPE_STRING)
			minCpnId = LocalGetNumObjectId(CpnCurves[4]);
		else
		{
			minCpn = atof(CpnCurves[4].c_str());
			minCpnId = ARM_NULL_OBJECT;
		}

		//max Cpn
		double maxCpn;
		long     maxCpnId;
		if(CpnCurvesType[5] == XL_TYPE_STRING)
			maxCpnId = LocalGetNumObjectId(CpnCurves[5]);
		else
		{
			maxCpn = atof(CpnCurves[5].c_str());
			maxCpnId = ARM_NULL_OBJECT;
		}

		//Cap/Floor
		CCString CFStr;
		XL_readStrCell(XL_CapOrFloor,CFStr,
			" ARM_ERR: Cap Or Floor: string expected",C_result);
		long CF;
		if((CF = ARM_ConvCapOrFloor (CFStr, C_result)) == ARM_DEFAULT_ERR)
		{
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}

		// FUNDING STRUCTURE
		VECTOR<CCString> RequiredFundDatas;
		XL_readStrVector(XL_fundDatas,RequiredFundDatas," ARM_ERR: Required Funding Datas : array of string expected",DOUBLE_TYPE,C_result);
        
		if(RequiredFundDatas.size()!=4)
		{
			throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"Funding Data size != 4.");
		}

		/// fundIdxTerm
		CCString fundFreqStr = RequiredFundDatas[0];
		long fundFreq;
		if((fundFreq = ARM_ConvFrequency (fundFreqStr, C_result)) == ARM_DEFAULT_ERR)
		{
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}

		/// fundDayCount
		CCString fundDaycountStr = RequiredFundDatas[1];
		long fundDaycount = ARM_ConvDayCount (fundDaycountStr);
		
		// Funding Coeff
		CCString fundCoeffStr = RequiredFundDatas[2];
		double fundCoeff;
		long     fundCoeffId;
		if(fundCoeffStr[0] == 'L')
			fundCoeffId = LocalGetNumObjectId(fundCoeffStr);
		else
		{
			fundCoeffId = ARM_NULL_OBJECT;
			fundCoeff = atof(fundCoeffStr);
		}

		// Funding Margin
		double fundMargin;
		CCString fundMarginStr = RequiredFundDatas[3];
		long     fundMarginId;

		if(fundMarginStr[0] == 'L')
			fundMarginId = LocalGetNumObjectId(fundMarginStr);
		else
		{
			fundMarginId = ARM_NULL_OBJECT;
			fundMargin = atof(fundMarginStr);
		}

    	
		VECTOR<CCString> ExerciseData;
		XL_readStrVector(XL_ExerciseData,ExerciseData," ARM_ERR: Required Exercise datas: array of string expected", DOUBLE_TYPE, C_result);

		long NotifDays = atof(ExerciseData[0]);
		long NonCall = atof(ExerciseData[1]);

		// exercise calendar
		string exerciseCal = CCSTringToSTLString(ExerciseData[2]);

		bool callSBOrStacky = true;
		if (ExerciseData.size() >= 4)
		{
			string callSBOrStackyStr = CCSTringToSTLString(ExerciseData[3]);
			stringGetUpper(callSBOrStackyStr);

			if (callSBOrStackyStr == "SB")
				callSBOrStacky = true;
			if (callSBOrStackyStr == "STACKY")
				callSBOrStacky = false;
		}

		// refvalue des fees
		CCString feesStr;
		XL_readStrCell(XL_fees,feesStr,
			" ARM_ERR: fees : refvalue expected",C_result);

		
	/// ------------------------
	/// Model datas
	/// ------------------------
	
		//Model Datas
		VECTOR<CCString> ModelDatas;
		XL_readStrVector(XL_ModelDatas,ModelDatas," ARM_ERR: Model Datas: array of string expected",DOUBLE_TYPE, C_result);

		long NbPathBounding;
		long NbPathPricing;
		long NbMaxBucket;
		long FixSteps;
		string USMethod;
				
		if (ModelDatas.size() < 5 || ModelDatas.size() > 20)
			throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"Model Data size is required to be between 5 and 20");
		
		NbPathBounding = atof(ModelDatas[0]);
		NbPathPricing = atof(ModelDatas[1]);
		NbMaxBucket = atof(ModelDatas[2]);
		FixSteps = atof(ModelDatas[3]);
		USMethod = CCSTringToSTLString(ModelDatas[4]);
		
		/// default generator 1 = MRGK5
		string GenType1 = "MRGK5"; 

		if (ModelDatas.size() >= 6)
		{
			GenType1 = CCSTringToSTLString(ModelDatas[5]);
			GenType1 = stringGetUpper(GenType1);
		}

		/// default model = SFRM
		string modelType = "SFRM2F";
		string triggerMode = "Coupon";

		string calibSwoptFreqStr;
		int calibSwoptFreq=-1;
		string regressors;

		if (ModelDatas.size() >= 7)
		{
			modelType = CCSTringToSTLString(ModelDatas[6]);
			modelType = stringGetUpper(modelType);
		}

		// TriggerMode
		if (ModelDatas.size() >= 8)
		{
			triggerMode = CCSTringToSTLString(ModelDatas[7]);
			triggerMode = stringGetUpper(triggerMode);
		}

		/// HK datas : PDE space & state steps, standard dev etc.
		vector< double > hkDatas(0);
		
		/// PDE NB TIME STEPS
		if (ModelDatas.size() >= 9)
		{
			hkDatas.push_back( atof(ModelDatas[8]) );
		}

		/// PDE NB STATE STEPS
		if (ModelDatas.size() >= 10)
		{
			hkDatas.push_back( atof(ModelDatas[9]) );
		}

		/// PDE NB STDEV
		if (ModelDatas.size() >= 11)
		{
			hkDatas.push_back( atof(ModelDatas[10]) );
		}

		/// SBGM NB FACTORS
		if (ModelDatas.size() >= 12)
		{
			hkDatas.push_back( atof(ModelDatas[11]) );
		}

		// CalibSwoptFreq
		if (ModelDatas.size() >= 13)
		{
			calibSwoptFreqStr = CCSTringToSTLString(ModelDatas[12]);
			calibSwoptFreqStr = stringGetUpper(calibSwoptFreqStr);
			if (calibSwoptFreqStr != "DEF")
				calibSwoptFreq = ARM_ArgConv_MatFrequency.GetNumber(calibSwoptFreqStr);
			else
				calibSwoptFreq = -1;
		}

		// Regressors
		if (ModelDatas.size() >= 14)
		{
			regressors = CCSTringToSTLString(ModelDatas[13]);
		}

		/// default generator 2 = Sobol
		string GenType2 = "Sobol"; 

		if (ModelDatas.size() >= 15)
		{
			GenType2 = CCSTringToSTLString(ModelDatas[14]);
			GenType2 = stringGetUpper(GenType2);
		}

		/// default first nb times = 0
		double FirstNbTimes; 

		if (ModelDatas.size() >= 16)
		{
			FirstNbTimes = atof(ModelDatas[15]);
		}

		/// default first nb dims = 0
		double FirstNbDims; 

		if (ModelDatas.size() >= 17)
		{
			FirstNbDims = atof(ModelDatas[16]);
		}

		/// default path scheme = Incremental
		string PathScheme = "Incremental"; 

		if (ModelDatas.size() >= 18)
		{
			PathScheme = CCSTringToSTLString(ModelDatas[17]);
			PathScheme = stringGetUpper(PathScheme);
		}

		/// default path order = Bucket Order
		string PathOrder = "BucketOrder";

		if (ModelDatas.size() >= 19)
		{
			PathOrder = CCSTringToSTLString(ModelDatas[18]);
			PathOrder = stringGetUpper(PathOrder);
		}

		/// default antithetic = "Y"
		string Antithetic = "Y";

		if (ModelDatas.size() >= 20)
		{
			Antithetic = CCSTringToSTLString(ModelDatas[19]);
			Antithetic = stringGetUpper(Antithetic);
		}
		
		//Calib Flags
		VECTOR<CCString> CalibFlagsXL;
		VECTOR<CCString> CalibFlags_default(3);
		CalibFlags_default[0]= CCString("CAP"); //calibrate on the CAP Volatility
		CalibFlags_default[1]= CCString("N"); //beta calibration
		CalibFlags_default[2]= CCString("N"); //Fix Boundary

		XL_readStrVectorWD (XL_CalibFlags,CalibFlagsXL,CalibFlags_default," ARM_ERR: Calib Mode: array of string expected",DOUBLE_TYPE,C_result);
		vector< string > CalibFlags(CalibFlagsXL.size());
		for(i=0;i<CalibFlagsXL.size();++i)
        {            
			CalibFlags[i]=CCSTringToSTLString(CalibFlagsXL[i]);
        }
		
		if (CalibFlags.size() == 5)
		{
			// By default no fix beta
			CalibFlags.push_back("N");
		}
		else if( CalibFlags.size() > 6 )
			throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"Calib Flags size != 5");

		/// ------------------------
		/// Flags 
		/// ------------------------
		VECTOR<CCString> ProdFlagsXL;
		VECTOR<CCString> ProdFlags_default(10,"N");
        XL_readStrVectorWD(XL_ProductsFlags,ProdFlagsXL,ProdFlags_default," ARM_ERR: flags: array of string expected",DOUBLE_TYPE,C_result);

		vector< string > ProdFlags(ProdFlagsXL.size());
		for(i=0;i<ProdFlagsXL.size();++i)
        {
            (ProdFlagsXL[i]).toUpper();
			ProdFlags[i]=CCSTringToSTLString(ProdFlagsXL[i]);
        }

		if (ProdFlagsXL.size() > 10)
			throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"Prod Flags size != 10");

		callableSBCalculatorFunc ourFunc(
				startDate,
				endDate,
				cpnFreq,
				cpnDaycount,
				cpnIndexTiming,
				cpnIdxTerm,
				cpnIndexDaycount,
				cpnResetCal,
				cpnPayCal,
				cpnIntRule,
				cpnIndexResetLag,
				payRec,
				CF,
				nominal,
				nominalId,
				fundNominal,
				fundNominalId,
				dConst,
				constId,
				lPrevCpn,
				lPrevCpnId,
				lNewOpt,
				lNewOptId,
				strikeOpt,
				strikeOptId,
				minCpn,
				minCpnId,
				maxCpn,
				maxCpnId,
				fundFreq,
				fundDaycount,
				fundCoeff,
				fundCoeffId,
				fundMargin,
				fundMarginId,
				NotifDays,
				NonCall,
				exerciseCal,
				callSBOrStacky,
				LocalGetNumObjectId(feesStr),
				NbPathBounding,
				NbPathPricing,
				NbMaxBucket,
				FixSteps,
				USMethod,
				CalibFlags,
				ProdFlags,
				mktDataManagerId,
				keys,
				GenType1,
				GenType2,
				FirstNbTimes,
				FirstNbDims,
				PathScheme,
				PathOrder,
				Antithetic,
				modelType,
				triggerMode,
				calibSwoptFreq,
				regressors,
				hkDatas);

		/// call the general function
		fillXL_Result( LOCAL_GC_CALLABLE_SB_CLASS, ourFunc, C_result, XL_result, PersistentInXL );
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_CallableSBCalculator_Common" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


///////////////////////////////////
/// version that takes into account 
/// previous creation of object
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_CallableSBCalculator_Create(
    LPXLOPER XL_startDate,
    LPXLOPER XL_endDate,
	LPXLOPER XL_mktDatas,
	LPXLOPER XL_payOrRec,
	LPXLOPER XL_coupon,
	LPXLOPER XL_Notional,
	LPXLOPER XL_CpnCurves,
	LPXLOPER XL_CapOrFloor,
	LPXLOPER XL_fundDatas,
	LPXLOPER XL_ExerciseData,
	LPXLOPER XL_fees,
	LPXLOPER XL_Calibflags,
	LPXLOPER XL_ProductsFlags,
	LPXLOPER XL_ModelDatas)
{
	ADD_LOG("Local_CallableSBCalculator_Create");
	bool PersistentInXL = true;

	return Local_CallableSBCalculator_Common(
    XL_startDate,
    XL_endDate,
	XL_mktDatas,
	XL_payOrRec,
	XL_coupon,
	XL_Notional,
	XL_CpnCurves,
	XL_CapOrFloor,
	XL_fundDatas,
	XL_ExerciseData,
	XL_fees,
	XL_Calibflags,
	XL_ProductsFlags,
	XL_ModelDatas,
	PersistentInXL );
}


///////////////////////////////////
/// version for non persistentInXl
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_PXL_CallableSBCalculator_Create(
	LPXLOPER XL_startDate,
    LPXLOPER XL_endDate,
	LPXLOPER XL_mktDatas,
	LPXLOPER XL_payOrRec,
	LPXLOPER XL_coupon,
	LPXLOPER XL_Notional,
	LPXLOPER XL_CpnCurves,
	LPXLOPER XL_CapOrFloor,
	LPXLOPER XL_fundDatas,
	LPXLOPER XL_ExerciseData,
	LPXLOPER XL_fees,
	LPXLOPER XL_Calibflags,
	LPXLOPER XL_ProductsFlags,
	LPXLOPER XL_ModelDatas)
{
	ADD_LOG("Local_PXL_CallableSBCalculator_Create");
	bool PersistentInXL = false;
	return Local_CallableSBCalculator_Common(
    XL_startDate,
    XL_endDate,
	XL_mktDatas,
	XL_payOrRec,
	XL_coupon,
	XL_Notional,
	XL_CpnCurves,
	XL_CapOrFloor,
	XL_fundDatas,
	XL_ExerciseData,
	XL_fees,
	XL_Calibflags,
	XL_ProductsFlags,
	XL_ModelDatas,
	PersistentInXL );
}



///--------------------------------------------------------
///--------------------------------------------------------
///             Callable SB Calculator Accessor
/// Inputs :
///     
///--------------------------------------------------------
///--------------------------------------------------------
class callableSBGetDataFunc : public ARMResultLong2LongFunc
{
public:
	callableSBGetDataFunc(
        long callableSBId,
        const string& getType)
    :
    C_callableSBId(callableSBId),
    C_getType(getType)
    {};
	
	long operator()( ARM_result& result, long objId ){
		return ARMLOCAL_CALLABLESB_Get(
            C_callableSBId,
            C_getType,
            result,
            objId);
    }

private:
	long    C_callableSBId;
	string  C_getType;
};

/////////////////////////////////////////////////////////////
/// central function that does the creation of the XL function
/////////////////////////////////////////////////////////////
LPXLOPER Local_CallableSBGet_Common(
	LPXLOPER XL_callableSBId,
	LPXLOPER XL_getType,
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

		long callableSBId;
		XL_GETOBJID( XL_callableSBId, callableSBId,	" ARM_ERR: Callable SnowBall Calculator: Object expected",C_result);

		CCString getTypeStr;
		XL_readStrCell(XL_getType,getTypeStr," ARM_ERR: Accessor Type: string expected",C_result);
        char* type=getTypeStr.c_str(); // à cause du new !!
		string getType(type);
        delete type;

		CCString callableSBGetClass(GCGetTypeToClass(getType,callableSBId).c_str());

		callableSBGetDataFunc ourFunc(callableSBId,getType);

		/// call the general function
		fillXL_Result( callableSBGetClass, ourFunc, C_result, XL_result, PersistentInXL );
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_CallableSBGet_Common" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


///////////////////////////////////
/// version that takes into account 
/// previous creation of object
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_CallableSBCalculator_GetData(
	LPXLOPER XL_tarnId,
	LPXLOPER XL_getType)
{
	ADD_LOG("Local_CallableSBCalculator_GetData");
	bool PersistentInXL = true;
	return Local_CallableSBGet_Common(
	    XL_tarnId,
	    XL_getType,
        PersistentInXL );
}


__declspec(dllexport) LPXLOPER WINAPI Local_PXL_CallableSBCalculator_GetData(
	LPXLOPER XL_tarnId,
	LPXLOPER XL_getType)
{
	ADD_LOG("Local_PXL_CallableSBCalculator_GetData");
	bool PersistentInXL = false;
	return Local_CallableSBGet_Common(
	    XL_tarnId,
	    XL_getType,
        PersistentInXL );
}

LPXLOPER Local_CallableSBSet_Common(
	LPXLOPER XL_csbId,
	LPXLOPER XL_dataId,
	LPXLOPER XL_setPortfolioType,
	LPXLOPER XL_MktDataKeys,
	bool isUpdated,
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

		long csbId;
		XL_GETOBJID( XL_csbId, csbId,	" ARM_ERR: CSB Calculator: Object expected",C_result);

		long dataId;
		XL_GETOBJID( XL_dataId, dataId,	" ARM_ERR: Object expected",C_result);

		CCString setPortfolioTypeStr;
		XL_readStrCellWD(XL_setPortfolioType,setPortfolioTypeStr,"DEFAULT"," ARM_ERR: Portfolio Type: string expected",C_result);
        char * portType = setPortfolioTypeStr.c_str(); // à cause du new !!
		string setPortfolioType(portType);
        delete portType;

		VECTOR<CCString> mktDataKeys;
		VECTOR<CCString> mktDataKeysDef(0);
		XL_readStrVectorWD(XL_MktDataKeys,mktDataKeys,mktDataKeysDef," ARM_ERR: Market datas keys: array of string expected",DOUBLE_TYPE,C_result);

		vector< string > mktDataKeysSTL(mktDataKeys.size());
		for(size_t i=0;i<mktDataKeys.size();++i)
        {
            mktDataKeys[i].toUpper();
			mktDataKeysSTL[i]=CCSTringToSTLString(mktDataKeys[i]);
        }

		exportFunc5Args< long, long, string, vector<string>, bool  >  ourFunc( csbId, dataId, setPortfolioType, mktDataKeysSTL, isUpdated, ARMLOCAL_CALLABLESB_Set );

        if(isUpdated)
        {
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
        else
        {
            /// General call through the functor with an object creation		    
		    fillXL_Result( LOCAL_GC_CRF_CLASS, ourFunc, C_result, XL_result, PersistentInXL );
        }
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_CRFSet_Common" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}

__declspec(dllexport) LPXLOPER WINAPI Local_CallableSBCalculator_SetData(
	LPXLOPER XL_csbId,
	LPXLOPER XL_dataId,
    LPXLOPER XL_setPortfolioType,
	LPXLOPER XL_mktDataKeys)
{
	ADD_LOG("Local_CallableSBCalculator_SetData");
	bool PersistentInXL = true;
    bool isUpdated = false;
	return Local_CallableSBSet_Common(
	    XL_csbId,
	    XL_dataId,
        XL_setPortfolioType,
		XL_mktDataKeys,
        isUpdated,
        PersistentInXL );
}

__declspec(dllexport) LPXLOPER WINAPI Local_PXL_CallableSBCalculator_SetData(
	LPXLOPER XL_csbId,
	LPXLOPER XL_dataId,
    LPXLOPER XL_setPortfolioType,
	LPXLOPER XL_mktDataKeys)
{
	ADD_LOG("Local_PXL_CallableSBCalculator_SetData");
	bool PersistentInXL = false;
    bool isUpdated = false;
	return Local_CallableSBSet_Common(
	    XL_csbId,
	    XL_dataId,
        XL_setPortfolioType,
		XL_mktDataKeys,
        isUpdated,
        PersistentInXL );
}

///////////////////////////////////
/// same as SetData but the previous CSB
/// is not cloned simply updated (=> no PXL version)
/// Only GC datas are updattable
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_CallableSBCalculator_Update(
	LPXLOPER XL_csbId,
	LPXLOPER XL_dataId,
    LPXLOPER XL_setPortfolioType,
	LPXLOPER XL_mktDataKeys)
{
	ADD_LOG("Local_CallableSBCalculator_Update");
	bool PersistentInXL = true;
    bool isUpdated = true;
	return Local_CallableSBSet_Common(
	    XL_csbId,
	    XL_dataId,
        XL_setPortfolioType,
		XL_mktDataKeys,
        isUpdated,
        PersistentInXL );
}


__declspec(dllexport) LPXLOPER WINAPI Local_CallableSBCalculator_GetPricingData (
	LPXLOPER XL_CallableSBCalculatorId,
    LPXLOPER XL_KeyId )
{
	ADD_LOG("Local_CallableSBCalculator_GetPricingData ");
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
		XL_readStrCell( XL_CallableSBCalculatorId, C_PricerString, " ARM_ERR: Caption Calculator Id: Object expected",C_result);
		long C_CallableSBCalculatorId = LocalGetNumObjectId(C_PricerString);
		
		CCString C_KeyString;
		XL_readStrCell( XL_KeyId, C_KeyString, " ARM_ERR: Key Id: Object expected",C_result);
		
		/// call the function
		ARM_GramFctorArg argResult;
		
		long retCode = ARMLOCAL_GC_GetPricingData(
			C_CallableSBCalculatorId,
			CCSTringToSTLString( C_KeyString ),
			argResult,
			C_result );
		
		if( retCode == ARM_KO )
		{
			ARM_ERR();
		}
		else
		{
			ARM_GramFunctorToXLOPER(argResult, XL_result, C_result);
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

///////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////MEPI Calculator///////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////

class CallOnMepiVanillaArgFunc : public ARMResultLong2LongFunc
{
public:
	CallOnMepiVanillaArgFunc ( const string& curveName, const string& EquityName, double startDate, double endDate, 
		long resetFreq, double riskFactor, double strike, double maxBorrow,  double protectionCurveStart, double protectionCurveEnd, 
		double startingPortfolio, double startingCash, double minInvested, double leverageCost, double cashSpread, double fees, double alreadyAsianed, long asianDatesNb) :
		C_CurveName( curveName ),
		C_EquityName( EquityName ),
		C_startDate( startDate ),
		C_endDate( endDate ),
		C_resetFreq( resetFreq ),
		C_strike( strike ),
		C_fees( fees ),
		C_riskFactor( riskFactor ),
		C_leverageCost( leverageCost ),
		C_startingCash( startingCash ),
		C_startingPortfolio( startingPortfolio ),
		C_protectionCurveStart( protectionCurveStart ),
		C_protectionCurveEnd( protectionCurveEnd ),
		C_cashSpread( cashSpread ),
		C_maxBorrow( maxBorrow ),
		C_minInvested( minInvested ),
		C_AlreadyAsianed( alreadyAsianed ),
		C_AsianDatesNb( asianDatesNb )
    {}
	
	long operator()( ARM_result& result, long objId ){
		return ARMLOCAL_CallOnMepiVanillaArgCreate(
			C_CurveName,
			C_EquityName,
			C_startDate,
			C_endDate,
			C_resetFreq,
			C_riskFactor,
			C_strike,
			C_maxBorrow,
			C_protectionCurveStart,
			C_protectionCurveEnd,
			C_startingPortfolio,
			C_startingCash,
			C_minInvested,
			C_leverageCost,
			C_cashSpread,
			C_fees,
			C_AlreadyAsianed,
			C_AsianDatesNb,
            result,
            objId);
    }

private:
	string C_CurveName;
	string C_EquityName;
	double C_startDate;
	double C_endDate;
	long C_resetFreq;
	double C_strike;
	double C_fees;
	double C_riskFactor;
	double C_leverageCost;
	double C_startingCash;
	double C_startingPortfolio;
	double C_protectionCurveStart;
	double C_protectionCurveEnd;
	double C_cashSpread;
	double C_maxBorrow;
	double C_minInvested;
	double C_AlreadyAsianed;
	long C_AsianDatesNb;
};


LPXLOPER Local_CallOnMepiVanillaArg_Common(
	LPXLOPER XL_CurveName,
	LPXLOPER XL_EquityModelName,
	LPXLOPER XL_startDate,
	LPXLOPER XL_endDate,
	LPXLOPER XL_ResetFreq,
	LPXLOPER XL_RiskFactor,
	LPXLOPER XL_Strike,
	LPXLOPER XL_MaxBorrow,
	LPXLOPER XL_ProtectionCurveStart,
	LPXLOPER XL_ProtectionCurveEnd,
	LPXLOPER XL_StartingPf,
	LPXLOPER XL_StartingCash,
	LPXLOPER XL_MinInvested,
	LPXLOPER XL_LvgCost,
	LPXLOPER XL_CashSpread,
	LPXLOPER XL_Fees,
	LPXLOPER XL_AlreadyAsianed,
	LPXLOPER XL_AsianDatesNb,
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

		///-------------------------------------------------------------------
		///                     Required Call On Mepi DATA				   ///
		/// ------------------------------------------------------------------
    
		/// General datas
		/// -----------------
		
		/// Start & End
		double startDate;
		XL_readNumCell(XL_startDate,startDate," ARM_ERR: start date: date expected",C_result);
		double endDate;
		XL_readNumCell(XL_endDate,endDate," ARM_ERR: end date: date expected",C_result);
		/// ResetFreq
		long resetFreq;
		XL_GETFREQUENCYWD(XL_ResetFreq,resetFreq, "A",	" ARM_ERR: reset freq: string expected",	C_result );
		/// Strike
		double strike;
		XL_readNumCell(XL_Strike,strike," ARM_ERR: strike: double expected",C_result);
		/// Fees
		double fees;
		XL_readNumCell(XL_Fees,fees," ARM_ERR: fees: double expected",C_result);
		/// RiskFactor
		double riskFactor;
		XL_readNumCell(XL_RiskFactor,riskFactor," ARM_ERR: riskFactor: double expected",C_result);
		/// LeverageCost
		double leverageCost;
		XL_readNumCell(XL_LvgCost,leverageCost," ARM_ERR: leverageCost: double expected",C_result);
		/// StartingCash
		double startingCash;
		XL_readNumCell(XL_StartingCash,startingCash," ARM_ERR: startingCash: double expected",C_result);
		/// StartingPf
		double startingPortfolio;
		XL_readNumCell(XL_StartingPf,startingPortfolio," ARM_ERR: startingPortfolio: double expected",C_result);
		/// cashSpread
		double cashSpread;
		XL_readNumCell(XL_CashSpread,cashSpread," ARM_ERR: cashSpread: double expected",C_result);
		/// MaxBorrow
		double maxBorrow;
		XL_readNumCell(XL_MaxBorrow,maxBorrow," ARM_ERR: maxBorrow: double expected",C_result);
		/// MinInvested
		double minInvested;
		XL_readNumCell(XL_MinInvested,minInvested," ARM_ERR: minInvested: double expected",C_result);
		/// Protection Curve Start
		double protectionCurveStart;
		XL_readNumCell(XL_ProtectionCurveStart,protectionCurveStart," ARM_ERR: protectionCurveStart: double expected",C_result);
		double protectionCurveEnd;
		XL_readNumCell(XL_ProtectionCurveEnd,protectionCurveEnd," ARM_ERR: protectionCurveEnd: double expected",C_result);
		double alreadyAsianed;
		XL_readNumCell(XL_AlreadyAsianed,alreadyAsianed," ARM_ERR: AlreadyAsianed: double expected",C_result);
		double AsianDatesNbDouble;
		XL_readNumCell(XL_AsianDatesNb,AsianDatesNbDouble," ARM_ERR: AlreadyAsianed: double expected",C_result);
		long AsianDatesNb = (long) AsianDatesNbDouble;


//	CallOnMepiVanillaArgFunc ( const string& curveName, const string& EquityName, double startDate, double endDate, 
//		long resetFreq, double riskFactor, double strike, double maxBorrow,  double protectionCurveStart, double startingPortfolio, 
//		double startingCash, double minInvested, double leverageCost, double cashSpread, double fees)

		CCString GensecType(LOCAL_GENSEC_CLASS);

		CallOnMepiVanillaArgFunc ourFunc( "EUR", "NTHG", startDate, endDate,resetFreq,  riskFactor,  strike, maxBorrow,  protectionCurveStart, protectionCurveEnd, startingPortfolio, startingCash, minInvested, leverageCost, cashSpread, fees, alreadyAsianed, AsianDatesNb);
		fillXL_Result( GensecType, ourFunc, C_result, XL_result, PersistentInXL );
	}

	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_CallOnMepiVanillaArg_Common" )

	return (LPXLOPER)&XL_result;
}


__declspec(dllexport) LPXLOPER WINAPI Local_CallOnMepiVanillaArg_Create(
	LPXLOPER XL_CurveName,
	LPXLOPER XL_EquityModelName,
	LPXLOPER XL_startDate,
	LPXLOPER XL_endDate,
	LPXLOPER XL_ResetFreq,
	LPXLOPER XL_RiskFactor,
	LPXLOPER XL_Strike,
	LPXLOPER XL_MaxBorrow,
	LPXLOPER XL_ProtectionCurveStart,
	LPXLOPER XL_ProtectionCurveEnd,
	LPXLOPER XL_StartingPf,
	LPXLOPER XL_StartingCash,
	LPXLOPER XL_MinInvested,
	LPXLOPER XL_LvgCost,
	LPXLOPER XL_CashSpread,
	LPXLOPER XL_Fees,
	LPXLOPER XL_AlreadyAsianed,
	LPXLOPER XL_AsianDatesNb)
{
	ADD_LOG("Local_CallOnMepiVanillaArg_Create");
	bool PersistentInXL = true;

	return Local_CallOnMepiVanillaArg_Common(
	XL_CurveName,
	XL_EquityModelName,
	XL_startDate,
	XL_endDate,
	XL_ResetFreq,
	XL_RiskFactor,
	XL_Strike,
	XL_MaxBorrow,
	XL_ProtectionCurveStart,
	XL_ProtectionCurveEnd,
	XL_StartingPf,
	XL_StartingCash,
	XL_MinInvested,
	XL_LvgCost,
	XL_CashSpread,
	XL_Fees, 
	XL_AlreadyAsianed,
	XL_AsianDatesNb,
	PersistentInXL );
}



///////////////////////////////////////////////////////////////////////////////////////////////
////////////////Callable SpreadOption///////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////
///----------------------------------------------
///----------------------------------------------
///         Callable SpreadOption
/// Inputs :
///     
///----------------------------------------------
///----------------------------------------------
class CSOCalculatorFunc : public ARMResultLong2LongFunc
{
public:
	CSOCalculatorFunc(
        double startDate,
        double endDate,
		long cpnFreq,
		long CMS1Type,
		long CMS2Type,
		long cpnDaycount,
		long fundFreq,
		long fundDaycount,
		long exerFreq,
		long NotifDays,
		double notional,
		long notionalId,
		double minCpn,
		long minCpnId,
		double maxCpn,
		long maxCpnId,
		double leverage,
		long leverageId,
		double fundSpread,
		long fundSpreadId,
		long FixCpnId,
		long feesId,
		const vector< string >& CalibMod,
		const vector< string >& ProductFlags,
		const vector< double >& ModelDatas,
		long mktDataManagerId,
        const vector< string >& keys)
    :
	C_startDate(startDate),
    C_endDate(endDate),
	C_CMS1Type(CMS1Type),
	C_CMS2Type(CMS2Type),
	C_cpnFreq(cpnFreq),
	C_cpnDaycount(cpnDaycount),
	C_fundFreq(fundFreq),
	C_fundDaycount(fundDaycount),
	C_exerFreq(exerFreq),
	C_NotifDays(NotifDays),
	C_notional(notional),
	C_notionalId(notionalId),
	C_minCpn(minCpn),
	C_minCpnId(minCpnId),
	C_maxCpn(maxCpn),
	C_maxCpnId(maxCpnId),
	C_leverage(leverage),
	C_leverageId(leverageId),
	C_fundMargin(fundSpread),
	C_fundMarginId(fundSpreadId),
	C_fixCpnId(FixCpnId),
	C_feesId(feesId),
	C_CalibMod(CalibMod),
	C_ProductFlags(ProductFlags),
	C_ModelDatas(ModelDatas),
	C_mktDataManagerId(mktDataManagerId),
    C_keys(keys)
    {};

	long operator()( ARM_result& result, long objId )
	{
		return ARMLOCAL_CSOCalculator_Create(
						C_startDate,
						C_endDate,
						C_CMS1Type,
						C_CMS2Type,
						C_cpnFreq,
						C_cpnDaycount,
						C_fundFreq,
						C_fundDaycount,
						C_exerFreq,
						C_NotifDays,
						C_notional,
						C_notionalId,
						C_minCpn,
						C_minCpnId,
						C_maxCpn,
						C_maxCpnId,
						C_leverage,
						C_leverageId,
						C_fundMargin,
						C_fundMarginId,
						C_fixCpnId,
						C_feesId,
						C_CalibMod,
						C_ProductFlags,
						C_ModelDatas,
						C_mktDataManagerId,
						C_keys,
						result,
						objId);
    }

private:

	double				C_startDate;
	double				C_endDate;
	long				C_CMS1Type;
	long				C_CMS2Type;
	long				C_cpnFreq;
	long				C_cpnDaycount;
	long				C_fundFreq;
	long				C_fundDaycount;
	long				C_exerFreq;
	long				C_NotifDays;
	double				C_notional;
	long				C_notionalId;
	double				C_minCpn;
	long				C_minCpnId;
	double				C_maxCpn;
	long				C_maxCpnId;
	double				C_leverage;
	long				C_leverageId;
	double				C_fundMargin;
	long				C_fundMarginId;
	long				C_fixCpnId;
	long				C_feesId;
	vector< string >	C_CalibMod;
	vector< string >	C_ProductFlags;
	vector< double >	C_ModelDatas;
	long				C_mktDataManagerId;
	vector< string >	C_keys;
};




LPXLOPER Local_CSOCalculator_Common(
	LPXLOPER XL_startDate,
    LPXLOPER XL_endDate,
	LPXLOPER XL_mktDatas,
	LPXLOPER XL_couponStructure,
	LPXLOPER XL_fundingStructure,
	LPXLOPER XL_exerciseStructure,
	LPXLOPER XL_Notional,
	LPXLOPER XL_CpnMin,
	LPXLOPER XL_CpnMax,
	LPXLOPER XL_Leverage,
	LPXLOPER XL_FixCpnCurve,
	LPXLOPER XL_FundSpread,
	LPXLOPER XL_Fees,
	LPXLOPER XL_Calibflags,
	LPXLOPER XL_ProductsFlags,
	LPXLOPER XL_ModelDatas,
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

		///-------------------------------------------------------------------
		///                     Required CSO DATA						   ///
		/// ------------------------------------------------------------------
    
		/// General datas
		/// -----------------
		
		/// Start & End
		double startDate;
		XL_readNumCell(XL_startDate,startDate," ARM_ERR: start date: date expected",C_result);
		double endDate;
		XL_readNumCell(XL_endDate,endDate," ARM_ERR: end date: date expected",C_result);

		/// Market datas : zc curve , MRS, beta, volatility curve and market models for OSW & CF pricings
		VECTOR<CCString> mktDatas;
		XL_readStrVector(XL_mktDatas,mktDatas," ARM_ERR: Market datas: array of string expected",DOUBLE_TYPE,C_result);
        if(mktDatas.size()<1)
		{
			throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"Mkt Data size should be greater than 1");
		}

		long mktDataManagerId = LocalGetNumObjectId (mktDatas[0]);
		vector< string > keys(mktDatas.size()-1);
        size_t i;
		for(i=1;i<mktDatas.size();++i)
			keys[i-1]=CCSTringToSTLString(mktDatas[i]);


		VECTOR<CCString> RequiredCpnDatas;
		XL_readStrVector(XL_couponStructure,RequiredCpnDatas," ARM_ERR: Required Cpn Datas : array of string expected",DOUBLE_TYPE,C_result);
        
		if(RequiredCpnDatas.size()!=4)
		{
			throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"Cpn Datas size != ))");
		}
		/// cpnFreq
		CCString cpnFreqStr = RequiredCpnDatas[0];
		long cpnFreq;
		if((cpnFreq = ARM_ConvFrequency (cpnFreqStr, C_result)) == ARM_DEFAULT_ERR)
		{
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}

		/// CMSType1
		CCString CMSType1Str = RequiredCpnDatas[1];
		long CMS1Type = ARM_ConvIrType (CMSType1Str);

		/// CMSType2
		CCString CMSType2Str = RequiredCpnDatas[2];
		long CMS2Type = ARM_ConvIrType (CMSType2Str);

		/// cpnDaycount
		CCString cpnDaycountStr = RequiredCpnDatas[3];
		long cpnDaycount = ARM_ConvDayCount (cpnDaycountStr);

		// end cpn structure

		VECTOR<CCString> RequiredFundDatas;
		XL_readStrVector(XL_fundingStructure,RequiredFundDatas," ARM_ERR: Required Funding Datas : array of string expected",DOUBLE_TYPE,C_result);
        
		if(RequiredFundDatas.size()!=2)
		{
			throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"Funding Datas size != ))");
		}

		/// fundFreq
		CCString fundFreqStr = RequiredFundDatas[0];
		long fundFreq;
		if((fundFreq = ARM_ConvFrequency (fundFreqStr, C_result)) == ARM_DEFAULT_ERR)
		{
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}

		/// fundDaycount
		CCString fundDaycountStr = RequiredFundDatas[1];
		long fundDaycount = ARM_ConvDayCount (fundDaycountStr);


		VECTOR<CCString> RequiredExerDatas;
		XL_readStrVector(XL_exerciseStructure,RequiredExerDatas," ARM_ERR: Required Exercise Datas : array of string expected",DOUBLE_TYPE,C_result);
        
		if(RequiredExerDatas.size()!=2)
		{
			throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"Exercise Datas size != ))");
		}

		/// fundFreq
		CCString exerFreqStr = RequiredExerDatas[0];
		long exerFreq;
		if((exerFreq = ARM_ConvFrequency (exerFreqStr, C_result)) == ARM_DEFAULT_ERR)
		{
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}

		long exerNoticeGap = (long)atof(RequiredExerDatas[1]);

		//Notional		
		double notional;
		double default_notional=100.0;
		CCString notionalStr;
		long     notionalId;
		XL_readStrOrNumCellWD(XL_Notional,notionalStr,notional,default_notional,notionalId,
			" ARM_ERR: notional: numerical or Curve Id expected",C_result);
		if(notionalId == XL_TYPE_STRING)
			notionalId = LocalGetNumObjectId(notionalStr);
		else
			notionalId = ARM_NULL_OBJECT;

		//CpnMin		
		double CpnMin;
		double default_cpnmin=0.0;
		CCString cpnminStr;
		long     cpnminId;
		XL_readStrOrNumCellWD(XL_CpnMin,cpnminStr,CpnMin,default_cpnmin,cpnminId,
			" ARM_ERR: CpnMin: numerical or Curve Id expected",C_result);
		if(cpnminId == XL_TYPE_STRING)
			cpnminId = LocalGetNumObjectId(cpnminStr);
		else
			cpnminId = ARM_NULL_OBJECT;

		//CpnMax		
		double CpnMax;
		double default_cpnmax=0.0;
		CCString cpnmaxStr;
		long     cpnmaxId;
		XL_readStrOrNumCellWD(XL_CpnMax,cpnmaxStr,CpnMax,default_cpnmax,cpnmaxId,
			" ARM_ERR: CpnMax: numerical or Curve Id expected",C_result);
		if(cpnmaxId == XL_TYPE_STRING)
			cpnmaxId = LocalGetNumObjectId(cpnmaxStr);
		else
			cpnmaxId = ARM_NULL_OBJECT;

		//Leverage		
		double Leverage;
		double default_Leverage=1.0;
		CCString LeverageStr;
		long     LeverageId;
		XL_readStrOrNumCellWD(XL_Leverage,LeverageStr,Leverage,default_Leverage,LeverageId,
			" ARM_ERR: Leverage: numerical or Curve Id expected",C_result);
		if(LeverageId == XL_TYPE_STRING)
			LeverageId = LocalGetNumObjectId(LeverageStr);
		else
			LeverageId = ARM_NULL_OBJECT;

		//FixCpnCurve		
		CCString FixCpnCurveStr;
		long     FixCpnCurveId;
		XL_readStrCell(XL_FixCpnCurve,FixCpnCurveStr," ARM_ERR: FixCpnCurve: Curve Id expected",C_result);
		FixCpnCurveId = LocalGetNumObjectId(FixCpnCurveStr);

		//fundingNominal		
		double FundSpread;
		double default_FundSpread=0.0;
		CCString FundSpreadStr;
		long     FundSpreadId;
		XL_readStrOrNumCellWD(XL_FundSpread,FundSpreadStr,FundSpread,default_FundSpread,FundSpreadId,
			" ARM_ERR: FundSpread: numerical or Curve Id expected",C_result);
		if(FundSpreadId == XL_TYPE_STRING)
			FundSpreadId = LocalGetNumObjectId(FundSpreadStr);
		else
			FundSpreadId = ARM_NULL_OBJECT;

		//Fees		
		CCString FeesStr;
		long     FeesId;
		XL_readStrCell(XL_Fees,FeesStr," ARM_ERR: Fees: Curve Id expected",C_result);
		FeesId = LocalGetNumObjectId(FeesStr);

		//Calib Flags
		VECTOR<CCString> CalibFlagsXL;
		VECTOR<CCString> CalibFlags_default(4);
		CalibFlags_default[0]= CCString("BESTFIT"); //calibrate on the CAP Volatility
		CalibFlags_default[1]= CCString("1D"); //beta calibration
		CalibFlags_default[2]= CCString("Y"); //Underlying Adj
		CalibFlags_default[3]= CCString("Y"); //Underlying Adj

		XL_readStrVectorWD (XL_Calibflags,CalibFlagsXL,CalibFlags_default," ARM_ERR: Calib Mode: array of string expected",DOUBLE_TYPE,C_result);
		vector< string > CalibFlags(CalibFlagsXL.size());
		for(i=0;i<CalibFlagsXL.size();++i)
        {            
			CalibFlags[i]=CCSTringToSTLString(CalibFlagsXL[i]);
        }
		
		//
		// already tested in ARM_loca_gp_calculator.cpp
		// the size of CalibFlags will tell what model we use (SFRM or HW)
		// this will has to be modified ...
		// if( CalibFlags.size() != 4 )
		//	throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"Calib Flags size != 4");


		/// ------------------------
		/// Flags 
		/// ------------------------
		VECTOR<CCString> ProdFlagsXL;
		VECTOR<CCString> ProdFlags_default(2,"N");
        XL_readStrVectorWD(XL_ProductsFlags,ProdFlagsXL,ProdFlags_default," ARM_ERR: flags: array of string expected",DOUBLE_TYPE,C_result);

		vector< string > ProdFlags(ProdFlagsXL.size());
		for(i=0;i<ProdFlagsXL.size();++i)
        {
            (ProdFlagsXL[i]).toUpper();
			ProdFlags[i]=CCSTringToSTLString(ProdFlagsXL[i]);
        }

		if (ProdFlagsXL.size() > 5)
			throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"Prod Flags : invalid ");

		/// ------------------------
		/// Model Datas 
		/// ------------------------
		VECTOR<double> ModelDatasXL;
		VECTOR<double> ModelDatas_default(5);
		ModelDatas_default[0] = 5;
		ModelDatas_default[1] = 0.000001;
		ModelDatas_default[2] = 5;
		ModelDatas_default[3] = 0.00001;
		ModelDatas_default[4] = 0.1;

        XL_readNumVectorWD(XL_ModelDatas,ModelDatasXL,ModelDatas_default," ARM_ERR: ModelDatas: array of numeric expected",C_result);

		vector< double > ModelDatas(ModelDatasXL.size());
		for(i=0;i<ModelDatasXL.size();++i)
        {
			ModelDatas[i]=ModelDatasXL[i];
        }

		if (ModelDatasXL.size() > 5)
			throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"ModelDatas size != 5");


		CSOCalculatorFunc ourFunc(
				startDate,
				endDate,
				cpnFreq,
				CMS1Type,
				CMS2Type,
				cpnDaycount,
				fundFreq,
				fundDaycount,
				exerFreq,
				exerNoticeGap,
				notional,
				notionalId,
				CpnMin,
				cpnminId,
				CpnMax,
				cpnmaxId,
				Leverage,
				LeverageId,
				FundSpread,
				FundSpreadId,
				FixCpnCurveId,
				FeesId,
				CalibFlags,
				ProdFlags,
				ModelDatas,
				mktDataManagerId,
				keys);

		/// call the general function
		fillXL_Result( LOCAL_GC_CALLABLE_SO_CLASS, ourFunc, C_result, XL_result, PersistentInXL );
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_CSOCalculator_Common" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}





///////////////////////////////////
/// version that takes into account 
/// previous creation of object
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_CSOCalculator_Create(
	LPXLOPER XL_startDate,
    LPXLOPER XL_endDate,
	LPXLOPER XL_mktDatas,
	LPXLOPER XL_couponStructure,
	LPXLOPER XL_fundingStructure,
	LPXLOPER XL_exerciseStructure,
	LPXLOPER XL_Notional,
	LPXLOPER XL_CpnMin,
	LPXLOPER XL_CpnMax,
	LPXLOPER XL_Leverage,
	LPXLOPER XL_FixCpnCurve,
	LPXLOPER XL_FundSpread,
	LPXLOPER XL_fees,
	LPXLOPER XL_Calibflags,
	LPXLOPER XL_ProductsFlags,
	LPXLOPER XL_ModelDatas)
{
	ADD_LOG("Local_CSOCalculator_Create");
	bool PersistentInXL = true;

	return Local_CSOCalculator_Common(
    XL_startDate,
    XL_endDate,
	XL_mktDatas,
	XL_couponStructure,
	XL_fundingStructure,
	XL_exerciseStructure,
	XL_Notional,
	XL_CpnMin,
	XL_CpnMax,
	XL_Leverage,
	XL_FixCpnCurve,
	XL_FundSpread,
	XL_fees,
	XL_Calibflags,
	XL_ProductsFlags,
	XL_ModelDatas,
	PersistentInXL );
}


///////////////////////////////////
/// version for non persistentInXl
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_PXL_CSOCalculator_Create(
	LPXLOPER XL_startDate,
    LPXLOPER XL_endDate,
	LPXLOPER XL_mktDatas,
	LPXLOPER XL_couponStructure,
	LPXLOPER XL_fundingStructure,
	LPXLOPER XL_exerciseStructure,
	LPXLOPER XL_Notional,
	LPXLOPER XL_CpnMin,
	LPXLOPER XL_CpnMax,
	LPXLOPER XL_Leverage,
	LPXLOPER XL_FixCpnCurve,
	LPXLOPER XL_FundSpread,
	LPXLOPER XL_fees,
	LPXLOPER XL_Calibflags,
	LPXLOPER XL_ProductsFlags,
	LPXLOPER XL_ModelDatas)
{
	ADD_LOG("Local_PXL_CSOCalculator_Create");
	bool PersistentInXL = false;

	return Local_CSOCalculator_Common(
    XL_startDate,
    XL_endDate,
	XL_mktDatas,
	XL_couponStructure,
	XL_fundingStructure,
	XL_exerciseStructure,
	XL_Notional,
	XL_CpnMin,
	XL_CpnMax,
	XL_Leverage,
	XL_FixCpnCurve,
	XL_FundSpread,
	XL_fees,
	XL_Calibflags,
	XL_ProductsFlags,
	XL_ModelDatas,
	PersistentInXL );
}


///-------------------------------------------------------
///-------------------------------------------------------
///    Callable SpreadOption : Extended Version
///		--> leverage long != leverage short	
///		--> strike
///		--> possible non standard reset timing on both legs
///
/// Inputs :
///     
///-------------------------------------------------------
///------------------------------------------------------///
/// Excel add-in common for extended CSO calculator
///
LPXLOPER Local_ExtendedCSOCalculator_Common(
	LPXLOPER XL_startDate,
    LPXLOPER XL_endDate,
	LPXLOPER XL_mktDatas,
	LPXLOPER XL_couponStructure,
	LPXLOPER XL_fundingStructure,
	LPXLOPER XL_exerciseStructure,
	LPXLOPER XL_cpnNotional,
	LPXLOPER XL_CpnMin,
	LPXLOPER XL_CpnMax,
	LPXLOPER XL_LeverageLong,
	LPXLOPER XL_LeverageShort,
	LPXLOPER XL_Strike,
	LPXLOPER XL_FixCpnCurve,
	LPXLOPER XL_FundSpread,
	LPXLOPER XL_Fees,
	LPXLOPER XL_Calibflags,
	LPXLOPER XL_ProductsFlags,
	LPXLOPER XL_ModelDatas,
	LPXLOPER XL_FundLeverage,
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

		///-------------------------------------------------------------------
		///                     Required CSO DATA						   ///
		/// ------------------------------------------------------------------
    
		/// General datas
		/// -----------------
		
		/// Start & End
		double startDate;
		XL_readNumCell(XL_startDate,startDate," ARM_ERR: start date: date expected",C_result);
		double endDate;
		XL_readNumCell(XL_endDate,endDate," ARM_ERR: end date: date expected",C_result);

		/// Market datas : zc curve , MRS, beta, volatility curve and market models for OSW & CF pricings
		VECTOR<CCString> mktDatas;
		XL_readStrVector(XL_mktDatas,mktDatas," ARM_ERR: Market datas: array of string expected",DOUBLE_TYPE,C_result);
        if(mktDatas.size()<1)
		{
			throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"Mkt Data size should be greater than 1");
		}

		long mktDataManagerId = LocalGetNumObjectId (mktDatas[0]);
		vector< string > keys(mktDatas.size()-1);
        size_t i;
		for(i=1;i<mktDatas.size();++i)
			keys[i-1]=CCSTringToSTLString(mktDatas[i]);


		VECTOR<CCString> RequiredCpnDatas;
		XL_readStrVector(XL_couponStructure,RequiredCpnDatas," ARM_ERR: Required Cpn Datas : array of string expected",DOUBLE_TYPE,C_result);
        
		if(RequiredCpnDatas.size() > 7) 
			throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"Cpn Datas size != 8, the last data is stubRule))");
		/// cpnFreq
		CCString cpnFreqStr = RequiredCpnDatas[0];
		long cpnFreq;
		if((cpnFreq = ARM_ConvFrequency (cpnFreqStr, C_result)) == ARM_DEFAULT_ERR)
		{
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}

		///CMSType1
		CCString CMSType1Str = RequiredCpnDatas[1];
		long CMS1Type = ARM_ConvIrType (CMSType1Str);

		///CMSType2
		CCString CMSType2Str = RequiredCpnDatas[2];
		long CMS2Type = ARM_ConvIrType (CMSType2Str);

		///cpnDaycount
		CCString cpnDaycountStr = RequiredCpnDatas[3];
		long cpnDaycount = ARM_ConvDayCount (cpnDaycountStr);

		///cpn reset timing (default = K_ADVANCE)
        long cpnResetTiming = (RequiredCpnDatas.size() == 5) ? ARM_ConvPayResetRule(RequiredCpnDatas[4]) : K_ADVANCE;

        ///pay and reset calendar
        string cpnResetCal = (RequiredCpnDatas.size() == 6) ? CCSTringToSTLString(RequiredCpnDatas[5]) : string ("");
        string cpnPayCal = (RequiredCpnDatas.size() == 7) ? CCSTringToSTLString(RequiredCpnDatas[6]) : string("");

		///end cpn structure
		VECTOR<CCString> RequiredFundDatas;
		XL_readStrVector(XL_fundingStructure,RequiredFundDatas," ARM_ERR: Required Funding Datas : array of string expected",DOUBLE_TYPE,C_result);
        
		if(RequiredFundDatas.size()!=2 && RequiredFundDatas.size()!=3 && RequiredFundDatas.size()!=4)
			throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"Funding Datas size != ))");

		///fundFreq
		CCString fundFreqStr = RequiredFundDatas[0];
		long fundFreq;
		if((fundFreq = ARM_ConvFrequency (fundFreqStr, C_result)) == ARM_DEFAULT_ERR)
		{
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}

		///fundDaycount
		CCString fundDaycountStr = RequiredFundDatas[1];
		long fundDaycount = ARM_ConvDayCount (fundDaycountStr);

		///funding reset timing (default = K_ADVANCE)
		long fundResetTiming = K_ADVANCE;		
		if(RequiredFundDatas.size() >= 3)
		{	
			CCString fundResetTimingStr = RequiredFundDatas[2];
			fundResetTiming = ARM_ConvPayResetRule(fundResetTimingStr);
		}

		///fundingType
		bool switchFlag = false;
		long fundingType = 0;
		if(RequiredFundDatas.size() == 4)
		{
			CCString fundingTypeStr = RequiredFundDatas[3];
			string fundingTypeStrStl = CCSTringToSTLString(fundingTypeStr);
			if (fundingTypeStrStl == "STANDARD")
				switchFlag = false;
			else
			{
				switchFlag = true;
				fundingType = ARM_ConvIrType (fundingTypeStr);
			}
		}

		///Exercise features
		VECTOR<CCString> RequiredExerDatas;
		XL_readStrVector(XL_exerciseStructure,RequiredExerDatas," ARM_ERR: Required Exercise Datas : array of string expected",DOUBLE_TYPE,C_result);
        
		if(RequiredExerDatas.size()<2)
			throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"Exercise Datas size != ))");

		/// exer Freq
		CCString exerFreqStr = RequiredExerDatas[0];
		long exerFreq;
		if((exerFreq = ARM_ConvFrequency (exerFreqStr, C_result)) == ARM_DEFAULT_ERR)
		{
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}
	
		long exerNoticeGap = (long)atof(RequiredExerDatas[1]);

		/// payer or receiver (default = RCV)
		long payRec = K_RCV;
		if(RequiredExerDatas.size()==3)
		{
			CCString payRecStr = RequiredExerDatas[2];
			if((payRec = ARM_ConvRecOrPay (payRecStr, C_result)) == ARM_DEFAULT_ERR)
			{
				ARM_ARG_ERR();
				return (LPXLOPER)&XL_result;
			}
		}

		//cpnNotional		
		double cpnNotional;
		double default_cpnNotional=100.0;
		CCString cpnNotionalStr;
		long     cpnNotionalId;
		XL_readStrOrNumCellWD(XL_cpnNotional,cpnNotionalStr,cpnNotional,default_cpnNotional,cpnNotionalId,
			" ARM_ERR: cpnNotional: numerical or Curve Id expected",C_result);
        cpnNotionalId = cpnNotionalId == XL_TYPE_STRING ? LocalGetNumObjectId(cpnNotionalStr) : ARM_NULL_OBJECT;

		//CpnMin		
		double CpnMin;
		double default_cpnmin=0.0;
		CCString cpnminStr;
		long     cpnminId;
		XL_readStrOrNumCellWD(XL_CpnMin,cpnminStr,CpnMin,default_cpnmin,cpnminId,
			" ARM_ERR: CpnMin: numerical or Curve Id expected",C_result);
        cpnminId = cpnminId == XL_TYPE_STRING ? LocalGetNumObjectId(cpnminStr) : ARM_NULL_OBJECT;

		//CpnMax		
		double CpnMax;
		double default_cpnmax=1.0e+10;
		CCString cpnmaxStr;
		long     cpnmaxId;
		XL_readStrOrNumCellWD(XL_CpnMax,cpnmaxStr,CpnMax,default_cpnmax,cpnmaxId,
			" ARM_ERR: CpnMax: numerical or Curve Id expected",C_result);
        cpnmaxId = cpnmaxId == XL_TYPE_STRING ? LocalGetNumObjectId(cpnmaxStr) : ARM_NULL_OBJECT;

		// LeverageLong		
		double LeverageLong;
		double default_LeverageLong=1.0;
		CCString LeverageLongStr;
		long     LeverageLongId;
		XL_readStrOrNumCellWD(XL_LeverageLong,LeverageLongStr,LeverageLong,default_LeverageLong,LeverageLongId,
			" ARM_ERR: LeverageLong: numerical or Curve Id expected",C_result);
        LeverageLongId = LeverageLongId == XL_TYPE_STRING ? LocalGetNumObjectId(LeverageLongStr):ARM_NULL_OBJECT;

		// LeverageShort		
		double LeverageShort;
		double default_LeverageShort=1.0;
		CCString LeverageShortStr;
		long     LeverageShortId;
		XL_readStrOrNumCellWD(XL_LeverageShort,LeverageShortStr,LeverageShort,default_LeverageShort,LeverageShortId,
			" ARM_ERR: LeverageShort: numerical or Curve Id expected",C_result);
        LeverageShortId = LeverageShortId == XL_TYPE_STRING ? LocalGetNumObjectId(LeverageShortStr): ARM_NULL_OBJECT;

        // Strike		
		double Strike;
		double default_Strike=0.0;
		CCString StrikeStr;
		long     StrikeId;
		XL_readStrOrNumCellWD(XL_Strike,StrikeStr,Strike,default_Strike,StrikeId,
			" ARM_ERR: Strike: numerical or Curve Id expected",C_result);
        StrikeId = StrikeId == XL_TYPE_STRING ? LocalGetNumObjectId(StrikeStr): ARM_NULL_OBJECT;

		//FixCpnCurve		
		CCString FixCpnCurveStr;
		long     FixCpnCurveId;
		XL_readStrCell(XL_FixCpnCurve,FixCpnCurveStr," ARM_ERR: FixCpnCurve: Curve Id expected",C_result);
		FixCpnCurveId = LocalGetNumObjectId(FixCpnCurveStr);

		//fundingNominal		
		double FundSpread;
		double default_FundSpread=0.0;
		CCString FundSpreadStr;
		long     FundSpreadId;
		XL_readStrOrNumCellWD(XL_FundSpread,FundSpreadStr,FundSpread,default_FundSpread,FundSpreadId,
			" ARM_ERR: FundSpread: numerical or Curve Id expected",C_result);
        FundSpreadId = FundSpreadId == XL_TYPE_STRING ? LocalGetNumObjectId(FundSpreadStr): ARM_NULL_OBJECT;

		//FundLeverage
		double FundLeverage;
		double default_FundLeverage=1.0;
		CCString FundLeverageStr;
		long     FundLeverageId;
		XL_readStrOrNumCellWD(XL_FundLeverage,FundLeverageStr,FundLeverage,default_FundLeverage,FundLeverageId,
			" ARM_ERR: FundLeverage: numerical or Curve Id expected",C_result);
        FundLeverageId = FundLeverageId == XL_TYPE_STRING ? LocalGetNumObjectId(FundLeverageStr): ARM_NULL_OBJECT;

		//Fees		
		CCString FeesStr;
		long     FeesId;
		XL_readStrCell(XL_Fees,FeesStr," ARM_ERR: Fees: Curve Id expected",C_result);
		FeesId = LocalGetNumObjectId(FeesStr);

		// Calib Flags
		VECTOR<CCString> CalibFlagsXL;
		VECTOR<CCString> CalibFlags_default(3);
		CalibFlags_default[0]= CCString("DIAG"); // diagonal calibration
		CalibFlags_default[1]= CCString("ATM");  // atm swaptions
		CalibFlags_default[2]= CCString("HWM1F");  // hw1f model
		
		XL_readStrVectorWD (XL_Calibflags,CalibFlagsXL,CalibFlags_default," ARM_ERR: Calib Mode: array of string expected",DOUBLE_TYPE,C_result);
		vector< string > CalibFlags(CalibFlagsXL.size());
		for(i=0;i<CalibFlagsXL.size();++i)
			CalibFlags[i]=CCSTringToSTLString(CalibFlagsXL[i]);
		
		if( CalibFlags.size() > 8 )
			throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"Bad Calib Flags size ");

		/// ------------------------
		/// Flags 
		/// ------------------------
		VECTOR<CCString> ProdFlagsXL;
		VECTOR<CCString> ProdFlags_default(6,"N");
        XL_readStrVectorWD(XL_ProductsFlags,ProdFlagsXL,ProdFlags_default," ARM_ERR: flags: array of string expected",DOUBLE_TYPE,C_result);

		vector< string > ProdFlags(ProdFlagsXL.size());
		for(i=0;i<ProdFlagsXL.size();++i)
        {
            (ProdFlagsXL[i]).toUpper();
			ProdFlags[i]=CCSTringToSTLString(ProdFlagsXL[i]);
        }

		if (ProdFlagsXL.size() > 6)
			throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"Prod Flags : invalid size");


		/// ------------------------
		/// Model Datas 
		/// ------------------------
		VECTOR<double> ModelDatasXL;
		VECTOR<double> ModelDatas_default(1,0);
        XL_readNumVectorWD(XL_ModelDatas,ModelDatasXL,ModelDatas_default," ARM_ERR: ModelDatas: array of numeric expected",C_result);
		vector< double > ModelDatas(ModelDatasXL);

		///if (ModelDatasXL.size() != 1)
		///	throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"ModelDatas size != 1");

        exportFunc38Args< double,double,long,long,long,long,long,long,long,long,long,long,long,double,long,double,long,double,long,double,long,
                          double,long,double,long,double,long,double,long,long,long,bool,long,vector< string >,vector< string >,vector< double >,long,vector< string > >
		ourFunc(
				startDate,
				endDate,
				CMS1Type,
				CMS2Type,
                cpnFreq,
				cpnDaycount,
				cpnResetTiming,
				fundFreq,
				fundDaycount,
				fundResetTiming,
				exerFreq,
				exerNoticeGap,
				payRec,
				cpnNotional,
				cpnNotionalId,
				CpnMin,
				cpnminId,
				CpnMax,
				cpnmaxId,
				LeverageLong,
				LeverageLongId,
				LeverageShort,
				LeverageShortId,
				Strike,
				StrikeId,
				FundSpread,
				FundSpreadId,
				FundLeverage,
				FundLeverageId,
				FixCpnCurveId,
				FeesId,
				switchFlag,
				fundingType,
				CalibFlags,
				ProdFlags,
				ModelDatas,
				mktDataManagerId,
				keys,ARMLOCAL_ExtendedCSOCalculator_Create);

		/// call the general function
		fillXL_Result( LOCAL_GC_CALLABLE_SO_CLASS, ourFunc, C_result, XL_result, PersistentInXL );
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_ExtendedCSOCalculator_Common" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}
///////////////////////////////////
/// version that takes into account 
/// previous creation of object
///////////////////////////////////

///////////////////////////////////
/// version that takes into account 
/// previous creation of object
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_ExtendedCSOCalculator_Create(
	LPXLOPER XL_startDate,
    LPXLOPER XL_endDate,
	LPXLOPER XL_mktDatas,
	LPXLOPER XL_couponStructure,
	LPXLOPER XL_fundingStructure,
	LPXLOPER XL_exerciseStructure,
	LPXLOPER XL_Notional,
	LPXLOPER XL_CpnMin,
	LPXLOPER XL_CpnMax,
	LPXLOPER XL_LeverageLong,
	LPXLOPER XL_LeverageShort,
	LPXLOPER XL_Strike,
	LPXLOPER XL_FixCpnCurve,
	LPXLOPER XL_FundSpread,
	LPXLOPER XL_fees,
	LPXLOPER XL_Calibflags,
	LPXLOPER XL_ProductsFlags,
	LPXLOPER XL_ModelDatas,
	LPXLOPER XL_FundLeverage)
{
	ADD_LOG("Local_ExtendedCSOCalculator_Create");
	bool PersistentInXL = true;

	return Local_ExtendedCSOCalculator_Common(
    XL_startDate,
    XL_endDate,
	XL_mktDatas,
	XL_couponStructure,
	XL_fundingStructure,
	XL_exerciseStructure,
	XL_Notional,
	XL_CpnMin,
	XL_CpnMax,
	XL_LeverageLong,
	XL_LeverageShort,
	XL_Strike,
	XL_FixCpnCurve,
	XL_FundSpread,
	XL_fees,
	XL_Calibflags,
	XL_ProductsFlags,
	XL_ModelDatas,
	XL_FundLeverage,
	PersistentInXL );
}

///////////////////////////////////
/// version for non persistentInXl
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_PXL_ExtendedCSOCalculator_Create(LPXLOPER XL_startDate,
                                                                             LPXLOPER XL_endDate,
	                                                                         LPXLOPER XL_mktDatas,
	                                                                         LPXLOPER XL_couponStructure,
	                                                                         LPXLOPER XL_fundingStructure,
	                                                                         LPXLOPER XL_exerciseStructure,
	                                                                         LPXLOPER XL_cpnNotional,
	                                                                         LPXLOPER XL_CpnMin,
	                                                                         LPXLOPER XL_CpnMax,
	                                                                         LPXLOPER XL_LeverageLong,
	                                                                         LPXLOPER XL_LeverageShort,
	                                                                         LPXLOPER XL_Strike,
	                                                                         LPXLOPER XL_FixCpnCurve,
	                                                                         LPXLOPER XL_FundSpread,
	                                                                         LPXLOPER XL_fees,
	                                                                         LPXLOPER XL_Calibflags,
	                                                                         LPXLOPER XL_ProductsFlags,
	                                                                         LPXLOPER XL_ModelDatas,
	                                                                         LPXLOPER XL_FundLeverage)
{
	ADD_LOG("Local_PXL_ExtendedCSOCalculator_Create");
	bool PersistentInXL = false;

	return Local_ExtendedCSOCalculator_Common(XL_startDate,
                                              XL_endDate,
	                                          XL_mktDatas,
	                                          XL_couponStructure,
	                                          XL_fundingStructure,
	                                          XL_exerciseStructure,
	                                          XL_cpnNotional,
	                                          XL_CpnMin,
	                                          XL_CpnMax,
	                                          XL_LeverageLong,
	                                          XL_LeverageShort,
	                                          XL_Strike,
	                                          XL_FixCpnCurve,
	                                          XL_FundSpread,
	                                          XL_fees,
	                                          XL_Calibflags,
	                                          XL_ProductsFlags,
	                                          XL_ModelDatas,
	                                          XL_FundLeverage,
	                                          PersistentInXL);
}



class basicCSOCalculatorFunc : public ARMResultLong2LongFunc
{
public:

    basicCSOCalculatorFunc(double asOfDate,
                           double startDate,
                           double fixEndDate,
                           double endDate,
                           long CMS1Type,
                           long CMS2Type,
                           long cpnFreq,
                           long cpnDaycount,
                           long cpnResetTiming,
                           long fundFreq,
                           long fundDaycount,
                           long fundResetTiming,
                           long exerFreq,
                           long NotifDays,
                           long payRec,
                           double cpnNotional,
                           long cpnNotionalId,
                           double minCpn,
                           long minCpnId,
                           double maxCpn,
                           long maxCpnId,
                           double leverageLong,
                           long leverageLongId,
                           double leverageShort,
                           long leverageShortId,
                           double strike,
                           long strikeId,
                           double fundNotional,
                           long fundNotionalId,
                           double fundMargin,
                           long fundMarginId,
		                   double fundLeverage,
                           long fundLeverageId,
                           long fixCpnId,
                           long feesId,
                           CCString CpnCcy,
                           CCString FundCcy,
						   long nbNoCall)
    :
	    C_asOfDate(asOfDate),
	    C_startDate(startDate),
        C_fixEndDate(fixEndDate),
        C_endDate(endDate),
	    C_CMS1Type(CMS1Type),
	    C_CMS2Type(CMS2Type),
	    C_cpnFreq(cpnFreq),
	    C_cpnDaycount(cpnDaycount),
	    C_cpnResetTiming(cpnResetTiming),
	    C_fundFreq(fundFreq),
	    C_fundDaycount(fundDaycount),
	    C_fundResetTiming(fundResetTiming),
	    C_exerFreq(exerFreq),
	    C_NotifDays(NotifDays),
	    C_payRec(payRec),
	    C_cpnNotional(cpnNotional),
	    C_cpnNotionalId(cpnNotionalId),
	    C_minCpn(minCpn),
	    C_minCpnId(minCpnId),
	    C_maxCpn(maxCpn),
	    C_maxCpnId(maxCpnId),
	    C_leverageLong(leverageLong),
	    C_leverageLongId(leverageLongId),
	    C_leverageShort(leverageShort),
	    C_leverageShortId(leverageShortId),
	    C_strike(strike),
	    C_strikeId(strikeId),
	    C_fundNotional(fundNotional),
	    C_fundNotionalId(fundNotionalId),
	    C_fundMargin(fundMargin),
	    C_fundMarginId(fundMarginId),
	    C_fundLeverage(fundLeverage),
	    C_fundLeverageId(fundLeverageId),
	    C_fixCpnId(fixCpnId),
	    C_feesId(feesId),
        C_CpnCcy(CpnCcy),
        C_FundCcy(FundCcy),
		C_nbNoCall(nbNoCall)
    {};

	long operator() (ARM_result& result, long objId )
	{
		return ARMLOCAL_BasicCSOCalculator_Create(C_asOfDate,
						                          C_startDate,
                                                  C_fixEndDate,
						                          C_endDate,
						                          C_CMS1Type,
						                          C_CMS2Type,
						                          C_cpnFreq,
						                          C_cpnDaycount,
						                          C_cpnResetTiming,
						                          C_fundFreq,
						                          C_fundDaycount,
						                          C_fundResetTiming,
						                          C_exerFreq,
						                          C_NotifDays,
						                          C_payRec,
						                          C_cpnNotional,
						                          C_cpnNotionalId,
						                          C_minCpn,
						                          C_minCpnId,
						                          C_maxCpn,
						                          C_maxCpnId,
						                          C_leverageLong,
						                          C_leverageLongId,
						                          C_leverageShort,
						                          C_leverageShortId,
						                          C_strike,
						                          C_strikeId,
						                          C_fundNotional,
						                          C_fundNotionalId,
						                          C_fundMargin,
						                          C_fundMarginId,
						                          C_fundLeverage,
						                          C_fundLeverageId,
						                          C_fixCpnId,
						                          C_feesId,
                                                  C_CpnCcy,
                                                  C_FundCcy,
												  C_nbNoCall,
						                          result,
						                          objId);
    }

private:

	double				C_asOfDate;
	double				C_startDate;
    double              C_fixEndDate;
	double				C_endDate;
	long				C_CMS1Type;
	long				C_CMS2Type;
	long				C_cpnFreq;
	long				C_cpnDaycount;
	long				C_cpnResetTiming;
	long				C_fundFreq;
	long				C_fundDaycount;
	long				C_fundResetTiming;
	long				C_exerFreq;
	long				C_NotifDays;
	long				C_payRec;
	double				C_cpnNotional;
	long				C_cpnNotionalId;
	double				C_minCpn;
	long				C_minCpnId;
	double				C_maxCpn;
	long				C_maxCpnId;
	double				C_leverageLong;
	long				C_leverageLongId;
	double				C_leverageShort;
	long				C_leverageShortId;
	double				C_strike;
	long				C_strikeId;
	double				C_fundNotional;
	long				C_fundNotionalId;
	double				C_fundMargin;
	long				C_fundMarginId;
	double				C_fundLeverage;
	long				C_fundLeverageId;
	long				C_fixCpnId;
	long				C_feesId;
	// basis case :
    CCString			C_CpnCcy;
    CCString			C_FundCcy;
	long				C_nbNoCall;
};



LPXLOPER Local_BasicCSOCalculator_Common(LPXLOPER XL_asOfDate,
	                                     LPXLOPER XL_startDate,
                                         LPXLOPER XL_endDate,
	                                     LPXLOPER XL_couponStructure,
	                                     LPXLOPER XL_fundingStructure,
	                                     LPXLOPER XL_exerciseStructure,
	                                     LPXLOPER XL_Notionals,
	                                     LPXLOPER XL_CpnMin,
	                                     LPXLOPER XL_CpnMax,
	                                     LPXLOPER XL_LeverageLong,
	                                     LPXLOPER XL_LeverageShort,
	                                     LPXLOPER XL_Strike,
	                                     LPXLOPER XL_FixCpnCurve,
	                                     LPXLOPER XL_FundSpread,
	                                     LPXLOPER XL_Fees,
	                                     LPXLOPER XL_FundLeverage,
                                         LPXLOPER XL_CpnCcy,
                                         LPXLOPER XL_FundCcy,
                                         LPXLOPER XL_FixEndDate,
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

		///-------------------------------------------------------------------
		///                     Required CSO DATA						   ///
		/// ------------------------------------------------------------------
    
		/// General datas
		/// -----------------
		
		/// AsOf, Start & End
		double asOfDate;
		XL_readNumCell(XL_asOfDate,asOfDate," ARM_ERR: as of date: date expected",C_result);
		double startDate;
		XL_readNumCell(XL_startDate,startDate," ARM_ERR: start date: date expected",C_result);
		double endDate;
		XL_readNumCell(XL_endDate,endDate," ARM_ERR: end date: date expected",C_result);

		VECTOR<CCString> RequiredCpnDatas;
		XL_readStrVector(XL_couponStructure,RequiredCpnDatas," ARM_ERR: Required Cpn Datas : array of string expected",DOUBLE_TYPE,C_result);
        
		if(RequiredCpnDatas.size() > 8) 
			throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"Cpn Datas size != 8, the last data is stubRule))");
		/// cpnFreq
		CCString cpnFreqStr = RequiredCpnDatas[0];
		long cpnFreq;
		if((cpnFreq = ARM_ConvFrequency (cpnFreqStr, C_result)) == ARM_DEFAULT_ERR)
		{
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}

		///CMSType1
		CCString CMSType1Str = RequiredCpnDatas[1];
		long CMS1Type = ARM_ConvIrType (CMSType1Str);

		///CMSType2
		CCString CMSType2Str = RequiredCpnDatas[2];
		long CMS2Type = ARM_ConvIrType (CMSType2Str);

		///cpnDaycount
		CCString cpnDaycountStr = RequiredCpnDatas[3];
		long cpnDaycount = ARM_ConvDayCount (cpnDaycountStr);

		///cpn reset timing (default = K_ADVANCE)
        long cpnResetTiming = (RequiredCpnDatas.size() == 5) ? ARM_ConvPayResetRule(RequiredCpnDatas[4]) : K_ADVANCE;

        ///pay and reset calendar
        string cpnResetCal = (RequiredCpnDatas.size() == 6) ? CCSTringToSTLString(RequiredCpnDatas[5]) : string ("");
        string cpnPayCal = (RequiredCpnDatas.size() == 7) ? CCSTringToSTLString(RequiredCpnDatas[6]) : string("");

		///end cpn structure
		VECTOR<CCString> RequiredFundDatas;
		XL_readStrVector(XL_fundingStructure,RequiredFundDatas," ARM_ERR: Required Funding Datas : array of string expected",DOUBLE_TYPE,C_result);
        
		if(RequiredFundDatas.size()!=2 && RequiredFundDatas.size()!=3)
			throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"Funding Datas size != ))");

		///fundFreq
		CCString fundFreqStr = RequiredFundDatas[0];
		long fundFreq;
		if((fundFreq = ARM_ConvFrequency (fundFreqStr, C_result)) == ARM_DEFAULT_ERR)
		{
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}

		///fundDaycount
		CCString fundDaycountStr = RequiredFundDatas[1];
		long fundDaycount = ARM_ConvDayCount (fundDaycountStr);

		///funding reset timing (default = K_ADVANCE)
		long fundResetTiming = K_ADVANCE;		
		if(RequiredFundDatas.size() == 3)
		{	
			CCString fundResetTimingStr = RequiredFundDatas[2];
			fundResetTiming = ARM_ConvPayResetRule(fundResetTimingStr);
		}

		///Exercise features
		VECTOR<CCString> RequiredExerDatas;
		XL_readStrVector(XL_exerciseStructure,RequiredExerDatas," ARM_ERR: Required Exercise Datas : array of string expected",DOUBLE_TYPE,C_result);
        
		if(RequiredExerDatas.size()<2)
			throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"Exercise Datas size should be > 2");

		/// exer Freq
		CCString exerFreqStr = RequiredExerDatas[0];
		long exerFreq;
		if((exerFreq = ARM_ConvFrequency (exerFreqStr, C_result)) == ARM_DEFAULT_ERR)
		{
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}
	
		long exerNoticeGap = (long)atof(RequiredExerDatas[1]);

		/// payer or receiver (default = RCV)
		long payRec = K_RCV;
		if(RequiredExerDatas.size() > 2)
		{
			CCString payRecStr = RequiredExerDatas[2];
			if((payRec = ARM_ConvRecOrPay (payRecStr, C_result)) == ARM_DEFAULT_ERR)
			{
				ARM_ARG_ERR();
				return (LPXLOPER)&XL_result;
			}
		}
		// number of non callable flows (at beginning of the deal)
		long nbNoCall = 0;
		if(RequiredExerDatas.size() > 3)
		{
			CCString nbNoCallStr = RequiredExerDatas[3];
			nbNoCall = atoi(nbNoCallStr);
		}

        // Manage Currencies
        CCString C_CpnCcy;
        XL_readStrCellWD(XL_CpnCcy, C_CpnCcy, "DEFAULT", " ARM_ERR: CpnCurrency: string expected",C_result);

        CCString C_FundCcy; // = coupon currency by default
        XL_readStrCellWD(XL_FundCcy, C_FundCcy, C_CpnCcy, " ARM_ERR: FundCurrency: string expected",C_result);

		// Notional (2 notionals for basis)
		VECTOR<CCString> NotionalIds;
		VECTOR<CCString> NotionalIds_default(2, "NULL");
		XL_readStrVectorWD(XL_Notionals,NotionalIds,NotionalIds_default," ARM_ERR: nationals Id : array of Ids expected",DOUBLE_TYPE,C_result);
		if(NotionalIds.size() > 2)
			throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT," notional Curve Ids size should be less or equal than 2");
		if(NotionalIds.size() ==1 && C_CpnCcy != C_FundCcy)
			throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"  ARM_ERR: Notional Ids expected, CSO needs two notionals");
		 
		long cpnNotionalId = ARM_NULL_OBJECT;
		double cpnNotional = 100.0;
		if (NotionalIds[0] != CCString("NULL"))
		{
			cpnNotionalId = LocalGetNumObjectId(NotionalIds[0]);
			if (cpnNotionalId == -1)
			{
				cpnNotionalId = ARM_NULL_OBJECT;
				cpnNotional = atoi(NotionalIds[0]);
			}
		}

		long fundNotionalId = ARM_NULL_OBJECT;
		double fundNotional = 100.0;
		if ( (NotionalIds.size() == 2) && (NotionalIds[1] != CCString("NULL")) )
		{
			fundNotionalId = LocalGetNumObjectId(NotionalIds[1]);
			if (fundNotionalId == -1)
			{
				fundNotionalId = ARM_NULL_OBJECT;
				fundNotional = atoi(NotionalIds[1]);
			}
		}

		//CpnMin		
		double CpnMin;
		double default_cpnmin=0.0;
		CCString cpnminStr;
		long     cpnminId;
		XL_readStrOrNumCellWD(XL_CpnMin,cpnminStr,CpnMin,default_cpnmin,cpnminId,
			" ARM_ERR: CpnMin: numerical or Curve Id expected",C_result);
        cpnminId = cpnminId == XL_TYPE_STRING ? LocalGetNumObjectId(cpnminStr) : ARM_NULL_OBJECT;

		//CpnMax		
		double CpnMax;
		double default_cpnmax=1.0e+10;
		CCString cpnmaxStr;
		long     cpnmaxId;
		XL_readStrOrNumCellWD(XL_CpnMax,cpnmaxStr,CpnMax,default_cpnmax,cpnmaxId,
			" ARM_ERR: CpnMax: numerical or Curve Id expected",C_result);
        cpnmaxId = cpnmaxId == XL_TYPE_STRING ? LocalGetNumObjectId(cpnmaxStr) : ARM_NULL_OBJECT;

		// LeverageLong		
		double LeverageLong;
		double default_LeverageLong=1.0;
		CCString LeverageLongStr;
		long     LeverageLongId;
		XL_readStrOrNumCellWD(XL_LeverageLong,LeverageLongStr,LeverageLong,default_LeverageLong,LeverageLongId,
			" ARM_ERR: LeverageLong: numerical or Curve Id expected",C_result);
        LeverageLongId = LeverageLongId == XL_TYPE_STRING ? LocalGetNumObjectId(LeverageLongStr):ARM_NULL_OBJECT;

		// LeverageShort		
		double LeverageShort;
		double default_LeverageShort=1.0;
		CCString LeverageShortStr;
		long     LeverageShortId;
		XL_readStrOrNumCellWD(XL_LeverageShort,LeverageShortStr,LeverageShort,default_LeverageShort,LeverageShortId,
			" ARM_ERR: LeverageShort: numerical or Curve Id expected",C_result);
        LeverageShortId = LeverageShortId == XL_TYPE_STRING ? LocalGetNumObjectId(LeverageShortStr): ARM_NULL_OBJECT;

        // Strike		
		double Strike;
		double default_Strike=0.0;
		CCString StrikeStr;
		long     StrikeId;
		XL_readStrOrNumCellWD(XL_Strike,StrikeStr,Strike,default_Strike,StrikeId,
			" ARM_ERR: Strike: numerical or Curve Id expected",C_result);
        StrikeId = StrikeId == XL_TYPE_STRING ? LocalGetNumObjectId(StrikeStr): ARM_NULL_OBJECT;

		//FixCpnCurve		
		CCString FixCpnCurveStr;
		long     FixCpnCurveId;
		XL_readStrCell(XL_FixCpnCurve,FixCpnCurveStr," ARM_ERR: FixCpnCurve: Curve Id expected",C_result);
		FixCpnCurveId = LocalGetNumObjectId(FixCpnCurveStr);

		//fundingSpread		
		double FundSpread;
		double default_FundSpread=0.0;
		CCString FundSpreadStr;
		long     FundSpreadId;
		XL_readStrOrNumCellWD(XL_FundSpread,FundSpreadStr,FundSpread,default_FundSpread,FundSpreadId,
			" ARM_ERR: FundSpread: numerical or Curve Id expected",C_result);
        FundSpreadId = FundSpreadId == XL_TYPE_STRING ? LocalGetNumObjectId(FundSpreadStr): ARM_NULL_OBJECT;

		//FundLeverage
		double FundLeverage;
		double default_FundLeverage=1.0;
		CCString FundLeverageStr;
		long     FundLeverageId;
		XL_readStrOrNumCellWD(XL_FundLeverage,FundLeverageStr,FundLeverage,default_FundLeverage,FundLeverageId,
			" ARM_ERR: FundLeverage: numerical or Curve Id expected",C_result);
        FundLeverageId = FundLeverageId == XL_TYPE_STRING ? LocalGetNumObjectId(FundLeverageStr): ARM_NULL_OBJECT;

		//Fees		
		CCString FeesStr;
		long     FeesId;
		XL_readStrCell(XL_Fees,FeesStr," ARM_ERR: Fees: Curve Id expected",C_result);
		FeesId = LocalGetNumObjectId(FeesStr);

        double fixEndDate;
        double fixEndDateDef = -1.0;
        XL_readNumCellWD(XL_FixEndDate, fixEndDate, fixEndDateDef,
                         " ARM_ERR: Fix. End Date: date expected", C_result);

		basicCSOCalculatorFunc ourFunc(asOfDate,
									   startDate,
                                       fixEndDate,
									   endDate,
									   CMS1Type,
									   CMS2Type,
									   cpnFreq,
									   cpnDaycount,
									   cpnResetTiming,
									   fundFreq,
									   fundDaycount,
									   fundResetTiming,
									   exerFreq,
									   exerNoticeGap,
									   payRec,
									   cpnNotional,
									   cpnNotionalId,
									   CpnMin,
									   cpnminId,
									   CpnMax,
									   cpnmaxId,
									   LeverageLong,
									   LeverageLongId,
									   LeverageShort,
									   LeverageShortId,
									   Strike,
									   StrikeId,
									   fundNotional,
									   fundNotionalId,
								       FundSpread,
									   FundSpreadId,
									   FundLeverage,
									   FundLeverageId,
									   FixCpnCurveId,
									   FeesId,
                                       C_CpnCcy,
                                       C_FundCcy,
									   nbNoCall);

		/// call the general function
		fillXL_Result(LOCAL_GC_CALLABLE_SO_CLASS, ourFunc, C_result, XL_result, PersistentInXL );
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_BasicCSOCalculator_Common" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}



__declspec(dllexport) LPXLOPER WINAPI Local_BasicCSOCalculator_Create(LPXLOPER XL_asOfDate,
	                                                                  LPXLOPER XL_startDate,
                                                                      LPXLOPER XL_endDate,
	                                                                  LPXLOPER XL_couponStructure,
	                                                                  LPXLOPER XL_fundingStructure,
	                                                                  LPXLOPER XL_exerciseStructure,
	                                                                  LPXLOPER XL_Notional,
	                                                                  LPXLOPER XL_CpnMin,
	                                                                  LPXLOPER XL_CpnMax,
	                                                                  LPXLOPER XL_LeverageLong,
	                                                                  LPXLOPER XL_LeverageShort,
	                                                                  LPXLOPER XL_Strike,
	                                                                  LPXLOPER XL_FixCpnCurve,
	                                                                  LPXLOPER XL_FundSpread,
	                                                                  LPXLOPER XL_fees,
	                                                                  LPXLOPER XL_FundLeverage,
                                                                      LPXLOPER XL_CpnCcy,
                                                                      LPXLOPER XL_FundCcy,
                                                                      LPXLOPER XL_FixEndDate)
{
	ADD_LOG("Local_BasicCSOCalculator_Create");
	bool PersistentInXL = true;

	return Local_BasicCSOCalculator_Common(XL_asOfDate,
                                           XL_startDate,
                                           XL_endDate,
	                                       XL_couponStructure,
	                                       XL_fundingStructure,
	                                       XL_exerciseStructure,
	                                       XL_Notional,
	                                       XL_CpnMin,
	                                       XL_CpnMax,
	                                       XL_LeverageLong,
	                                       XL_LeverageShort,
	                                       XL_Strike,
	                                       XL_FixCpnCurve,
	                                       XL_FundSpread,
	                                       XL_fees,
	                                       XL_FundLeverage,
                                           XL_CpnCcy,
                                           XL_FundCcy,
                                           XL_FixEndDate,
	                                       PersistentInXL);
}

__declspec(dllexport) LPXLOPER WINAPI Local_PXL_BasicCSOCalculator_Create(LPXLOPER XL_asOfDate,
																		  LPXLOPER XL_startDate,
																		  LPXLOPER XL_endDate,
																		  LPXLOPER XL_couponStructure,
																		  LPXLOPER XL_fundingStructure,
																		  LPXLOPER XL_exerciseStructure,
																		  LPXLOPER XL_Notional,
																		  LPXLOPER XL_CpnMin,
																		  LPXLOPER XL_CpnMax,
																		  LPXLOPER XL_LeverageLong,
																		  LPXLOPER XL_LeverageShort,
																		  LPXLOPER XL_Strike,
																		  LPXLOPER XL_FixCpnCurve,
																		  LPXLOPER XL_FundSpread,
																		  LPXLOPER XL_fees,
																		  LPXLOPER XL_FundLeverage,
																		  LPXLOPER XL_CpnCcy,
																		  LPXLOPER XL_FundCcy,
																		  LPXLOPER XL_FixEndDate)
{
	ADD_LOG("Local_PXL_BasicCSOCalculator_Create");
	bool PersistentInXL = false;

	return Local_BasicCSOCalculator_Common(XL_asOfDate,
                                           XL_startDate,
                                           XL_endDate,
	                                       XL_couponStructure,
	                                       XL_fundingStructure,
	                                       XL_exerciseStructure,
	                                       XL_Notional,
	                                       XL_CpnMin,
	                                       XL_CpnMax,
	                                       XL_LeverageLong,
	                                       XL_LeverageShort,
	                                       XL_Strike,
	                                       XL_FixCpnCurve,
	                                       XL_FundSpread,
	                                       XL_fees,
	                                       XL_FundLeverage,
                                           XL_CpnCcy,
                                           XL_FundCcy,
                                           XL_FixEndDate,
	                                       PersistentInXL);
}





///-------------------------------------------------------
///-------------------------------------------------------
///    Callable SpreadOption : Extended Version
///		--> leverage long != leverage short	
///		--> strike
///		--> possible non standard reset timing on both legs
///
/// Inputs :
///     
///-------------------------------------------------------
///------------------------------------------------------///
/// Excel add-in common for extended CSO calculator
///
LPXLOPER Local_BasisCSOCalculator_Common(
	LPXLOPER XL_startDate,
    LPXLOPER XL_fixEndDate,
	LPXLOPER XL_endDate,
	LPXLOPER XL_Ccys,
	LPXLOPER XL_mktDataMger,
	LPXLOPER XL_couponStructure,
	LPXLOPER XL_fundingStructure,
	LPXLOPER XL_exerciseStructure,
	LPXLOPER XL_Notionals,
	LPXLOPER XL_Cpns,
	LPXLOPER XL_Leverages,
	LPXLOPER XL_Strike,
	LPXLOPER XL_FundSpread,
	LPXLOPER XL_fees,
	LPXLOPER XL_Calibflags,
	LPXLOPER XL_ProductsFlags,
	LPXLOPER XL_ModelDatas,
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

		///-------------------------------------------------------------------
		///                     Required CSO DATA						   ///
		/// ------------------------------------------------------------------
    
		/// General datas
		/// -----------------
		
		/// Start & End
		double startDate;
		XL_readNumCell(XL_startDate,startDate," ARM_ERR: start date: date expected",C_result);
		double fixEndDate;
		XL_readNumCell(XL_fixEndDate,fixEndDate," ARM_ERR: Fix end date: date expected",C_result);
		double endDate;
		XL_readNumCell(XL_endDate,endDate," ARM_ERR: end date: date expected",C_result);

		VECTOR<CCString> Ccys;
		XL_readStrVector(XL_Ccys,Ccys," ARM_ERR: Currencies: array of string expected",DOUBLE_TYPE,C_result);
        if(Ccys.size() > 2)
			throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT," Currencies size should be less or equal than 2");
		if(Ccys.size() == 1) Ccys.push_back(Ccys[0]);
		string cpnCcy = CCSTringToSTLString(Ccys[0]);
		string fundCcy = CCSTringToSTLString(Ccys[1]);

		/// Market datas : zc curve , MRS, beta, volatility curve and market models for OSW & CF pricings
		CCString mktDataMger;
		XL_readStrCell(XL_mktDataMger,mktDataMger," ARM_ERR: Mkt Data Manager: mkt Data Manager Id expected",C_result);
		long mktDataManagerId = LocalGetNumObjectId (mktDataMger);

		VECTOR<CCString> RequiredCpnDatas;
		XL_readStrVector(XL_couponStructure,RequiredCpnDatas," ARM_ERR: Required Cpn Datas : array of string expected",DOUBLE_TYPE,C_result);
        
		if(RequiredCpnDatas.size() > 8) 
			throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"Cpn Datas size != 8, the last data is stubRule))");
		/// cpnFreq
		CCString cpnFreqStr = RequiredCpnDatas[0];
		long cpnFreq;
		if((cpnFreq = ARM_ConvFrequency (cpnFreqStr, C_result)) == ARM_DEFAULT_ERR)
		{
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}

		///CMSType1
		CCString CMSType1Str = RequiredCpnDatas[1];
		long CMS1Type = ARM_ConvIrType (CMSType1Str);

		///CMSType2
		CCString CMSType2Str = RequiredCpnDatas[2];
		long CMS2Type = ARM_ConvIrType (CMSType2Str);

		///cpnDaycount
		CCString cpnDaycountStr = RequiredCpnDatas[3];
		long cpnDaycount = ARM_ConvDayCount (cpnDaycountStr);

		///cpn reset timing (default = K_ADVANCE)
        long cpnResetTiming = (RequiredCpnDatas.size() > 5) ? ARM_ConvPayResetRule(RequiredCpnDatas[4]) : K_ADVANCE;

        ///pay and reset calendar
        string cpnResetCal = (RequiredCpnDatas.size() > 6) ? CCSTringToSTLString(RequiredCpnDatas[5]) : string ("");
        string cpnPayCal = (RequiredCpnDatas.size() > 7) ? CCSTringToSTLString(RequiredCpnDatas[6]) : string("");

		///stub Rule
		long stubRule = ARM_ConvStubRule((RequiredCpnDatas.size() > 8) ? RequiredCpnDatas[7] : CCString("SS"));		

		///end cpn structure
		VECTOR<CCString> RequiredFundDatas;
		XL_readStrVector(XL_fundingStructure,RequiredFundDatas," ARM_ERR: Required Funding Datas : array of string expected",DOUBLE_TYPE,C_result);
        
		if(RequiredFundDatas.size()!=2 && RequiredFundDatas.size()!=3)
			throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"Funding Datas size != ))");

		///fundFreq
		CCString fundFreqStr = RequiredFundDatas[0];
		long fundFreq;
		if((fundFreq = ARM_ConvFrequency (fundFreqStr, C_result)) == ARM_DEFAULT_ERR)
		{
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}

		///fundDaycount
		CCString fundDaycountStr = RequiredFundDatas[1];
		long fundDaycount = ARM_ConvDayCount (fundDaycountStr);

		///funding reset timing (default = K_ADVANCE)
		long fundResetTiming = K_ADVANCE;		
		if(RequiredFundDatas.size() == 3)
		{	
			CCString fundResetTimingStr = RequiredFundDatas[2];
			fundResetTiming = ARM_ConvPayResetRule(fundResetTimingStr);
		}

		///Exercise features
		VECTOR<CCString> RequiredExerDatas;
		XL_readStrVector(XL_exerciseStructure,RequiredExerDatas," ARM_ERR: Required Exercise Datas : array of string expected",DOUBLE_TYPE,C_result);
        
		if(RequiredExerDatas.size()<2)
			throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"Exercise Datas size != ))");

		/// exer Freq
		CCString exerFreqStr = RequiredExerDatas[0];
		long exerFreq;
		if((exerFreq = ARM_ConvFrequency (exerFreqStr, C_result)) == ARM_DEFAULT_ERR)
		{
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}
	
		long exerNoticeGap = (long)atof(RequiredExerDatas[1]);

		/// payer or receiver (default = RCV)
		long payRec = K_RCV;
		if(RequiredExerDatas.size()==3)
		{
			CCString payRecStr = RequiredExerDatas[2];
			if((payRec = ARM_ConvRecOrPay (payRecStr, C_result)) == ARM_DEFAULT_ERR)
			{
				ARM_ARG_ERR();
				return (LPXLOPER)&XL_result;
			}
		}
		long nbNCall = 0;
		if(RequiredExerDatas.size()==4)
			nbNCall = (long)atof(RequiredExerDatas[3]);

		//cpnNotional		
		VECTOR<CCString> NotionalIds;
		XL_readStrVector(XL_Notionals,NotionalIds," ARM_ERR: nationals Id : array of Ids expected",DOUBLE_TYPE,C_result);
		if(NotionalIds.size() > 2)
			throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT," notional Curve Ids size should be less or equal than 2");
		if(NotionalIds.size() ==1 && cpnCcy != fundCcy)
			throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"  ARM_ERR: Notional Ids expected, CSO needs two notionals");
		 
		long cpnNotionalId = LocalGetNumObjectId(NotionalIds[0]);
		long FundNotionalId = (NotionalIds.size() == 1 || cpnCcy == fundCcy) ? cpnNotionalId : LocalGetNumObjectId(NotionalIds[1]);


		//Cpn min andd Cpn max		
		VECTOR<CCString> CpnIds;
		XL_readStrVector(XL_Cpns,CpnIds," ARM_ERR: cpn min & cpn max Ids : array of Ids expected",DOUBLE_TYPE,C_result);
		if(CpnIds.size() != 2)
			throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT," cpn min & cpn max Curve Ids size should be equal than 2");
		long cpnminId = LocalGetNumObjectId(CpnIds[0]);
		long cpnmaxId =	LocalGetNumObjectId(CpnIds[1]);

		// Leverage long and short leverage		
		VECTOR<CCString> LeverageStrs;
		XL_readStrVector(XL_Leverages,LeverageStrs," ARM_ERR: Leverage max & Leverage short Ids : array of Ids expected",DOUBLE_TYPE,C_result);
		if(LeverageStrs.size() != 2)
			throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT," Leverage max & Leverage short Curve Ids size should be equal than 2");
		long LeverageLongId = LocalGetNumObjectId(LeverageStrs[0]);
		long LeverageShortId =	LocalGetNumObjectId(LeverageStrs[1]);

        // Strike		
		CCString StrikeStr;
		XL_readStrCell(XL_Strike,StrikeStr," ARM_ERR: Strike:Curve Id expected",C_result);
        long StrikeId = LocalGetNumObjectId(StrikeStr);

		//fundingNominal		
		double FundSpread;
		double default_FundSpread=0.0;
		CCString FundSpreadStr;
		long     FundSpreadId;
		XL_readStrOrNumCellWD(XL_FundSpread,FundSpreadStr,FundSpread,default_FundSpread,FundSpreadId,
			" ARM_ERR: FundSpread: numerical or Curve Id expected",C_result);
        FundSpreadId = FundSpreadId == XL_TYPE_STRING ? LocalGetNumObjectId(FundSpreadStr): ARM_NULL_OBJECT;

		//Fees		
		CCString FeesStr;
		long     FeesId;
		XL_readStrCell(XL_fees,FeesStr," ARM_ERR: Fees: Curve Id expected",C_result);
		FeesId = LocalGetNumObjectId(FeesStr);

		// Calib Flags
		VECTOR<CCString> CalibFlagsXL;
		VECTOR<CCString> CalibFlags_default(5);
		CalibFlags_default[0]= CCString("DIAG");		// diagonal calibration
		CalibFlags_default[1]= CCString("ATM");			// atm swaptions
		CalibFlags_default[2]= CCString("HWM1F");		// hw1f model
		CalibFlags_default[3]= CCString("MOYENESS");	// vns method
		CalibFlags_default[4]= CCString("1.0");			// moyeness level
		
		XL_readStrVectorWD (XL_Calibflags,CalibFlagsXL,CalibFlags_default," ARM_ERR: Calib Mode: array of string expected",DOUBLE_TYPE,C_result);
		vector< string > CalibDatas(CalibFlagsXL.size());
		for(size_t i=0;i<CalibFlagsXL.size();++i)
			CalibDatas[i]=CCSTringToSTLString(CalibFlagsXL[i]);
		
		/*if( CalibDatas.size() > 5 )
			throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"Bad Calib Flags size  ( size =! 5)");*/

		/// ------------------------
		/// Flags 
		/// ------------------------
		VECTOR<CCString> ProdFlagsXL;
		VECTOR<CCString> ProdFlags_default(5,"N");
        XL_readStrVectorWD(XL_ProductsFlags,ProdFlagsXL,ProdFlags_default," ARM_ERR: flags: array of string expected",DOUBLE_TYPE,C_result);

		vector< string > ProdFlags(ProdFlagsXL.size());
		for(i=0;i<ProdFlagsXL.size();++i)
        {
            (ProdFlagsXL[i]).toUpper();
			ProdFlags[i]=CCSTringToSTLString(ProdFlagsXL[i]);
        }

		if (ProdFlagsXL.size() > 5)
			throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"Prod Flags : invalid size");


		/// ------------------------
		/// Model Datas 
		/// ------------------------
		VECTOR<double> ModelDatasXL;
		VECTOR<double> ModelDatas_default(1,0);
        XL_readNumVectorWD(XL_ModelDatas,ModelDatasXL,ModelDatas_default," ARM_ERR: ModelDatas: array of numeric expected",C_result);
		vector< double > ModelDatas(ModelDatasXL);

        exportFunc33Args< double,double,double,string, string,long,long,long,long,long,long,string, string,long,long,long,long,long,long,long,long,long,long,long,long,
                          long,long,long,long,vector< string >,vector< string >,vector< double >,long >
		ourFunc(
				startDate,
				fixEndDate,
				endDate,
				cpnCcy,
				fundCcy,
				CMS1Type,
				CMS2Type,
                cpnFreq,
				cpnDaycount,
				cpnResetTiming,
				stubRule,
                cpnResetCal,
                cpnPayCal,
				cpnNotionalId,
				StrikeId,
				LeverageLongId,
				LeverageShortId,
				cpnminId,
				cpnmaxId,
				fundFreq,
				fundDaycount,
				fundResetTiming,
				FundNotionalId,
				FundSpreadId,
				exerFreq,
				exerNoticeGap,
				payRec,
				nbNCall,                
				FeesId,
				CalibDatas,
				ProdFlags,
				ModelDatas,
				mktDataManagerId,
				ARMLOCAL_BasisCSOCalculator_Create);

		/// call the general function
		fillXL_Result( LOCAL_GC_CALLABLE_SO_CLASS, ourFunc, C_result, XL_result, PersistentInXL );
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_BasisCSOCalculator_Common" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


///////////////////////////////////
/// version for non persistentInXl
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_BasisCSOCalculator_Create(
	LPXLOPER XL_startDate,
    LPXLOPER XL_fixEndDate,
	LPXLOPER XL_endDate,
	LPXLOPER XL_Ccys,
	LPXLOPER XL_mktDatas,
	LPXLOPER XL_couponStructure,
	LPXLOPER XL_fundingStructure,
	LPXLOPER XL_exerciseStructure,
	LPXLOPER XL_Notionals,
	LPXLOPER XL_Cpns,
	LPXLOPER XL_Leverages,
	LPXLOPER XL_Strike,
	LPXLOPER XL_FundSpread,
	LPXLOPER XL_fees,
	LPXLOPER XL_Calibflags,
	LPXLOPER XL_ProductsFlags,
	LPXLOPER XL_ModelDatas)
{
	ADD_LOG("Local_BasisCSOCalculator_Create");
	bool PersistentInXL = true;

	return Local_BasisCSOCalculator_Common(
			XL_startDate,
			XL_fixEndDate,
			XL_endDate,
			XL_Ccys,
			XL_mktDatas,
			XL_couponStructure,
			XL_fundingStructure,
			XL_exerciseStructure,
			XL_Notionals,
			XL_Cpns,
			XL_Leverages,
			XL_Strike,
			XL_FundSpread,
			XL_fees,
			XL_Calibflags,
			XL_ProductsFlags,
			XL_ModelDatas,
			PersistentInXL );
}


///////////////////////////////////
/// version that takes into account 
/// previous creation of object
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_PXL_BasisCSOCalculator_Create(
	LPXLOPER XL_startDate,
    LPXLOPER XL_fixEndDate,
	LPXLOPER XL_endDate,
	LPXLOPER XL_Ccys,
	LPXLOPER XL_mktDatas,
	LPXLOPER XL_couponStructure,
	LPXLOPER XL_fundingStructure,
	LPXLOPER XL_exerciseStructure,
	LPXLOPER XL_Notionals,
	LPXLOPER XL_Cpns,
	LPXLOPER XL_Leverages,
	LPXLOPER XL_Strike,
	LPXLOPER XL_FundSpread,
	LPXLOPER XL_fees,
	LPXLOPER XL_Calibflags,
	LPXLOPER XL_ProductsFlags,
	LPXLOPER XL_ModelDatas)
{
	ADD_LOG("Local_PXL_BasisCSOCalculator_Create");
	bool PersistentInXL = false;

	return Local_BasisCSOCalculator_Common(
			XL_startDate,
			XL_fixEndDate,
			XL_endDate,
			XL_Ccys,
			XL_mktDatas,
			XL_couponStructure,
			XL_fundingStructure,
			XL_exerciseStructure,
			XL_Notionals,
			XL_Cpns,
			XL_Leverages,
			XL_Strike,
			XL_FundSpread,
			XL_fees,
			XL_Calibflags,
			XL_ProductsFlags,
			XL_ModelDatas,
			PersistentInXL );
}



/////////////////////////////////////////////////////////////
/// Get datas from a CSO calculator
/////////////////////////////////////////////////////////////
LPXLOPER Local_CSOGet_Common(
	LPXLOPER XL_csoId,
	LPXLOPER XL_getType,
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

		long csoId;
		XL_GETOBJID( XL_csoId,csoId," ARM_ERR: CSO Calculator: Object expected",C_result);

		CCString getTypeStr;
		XL_readStrCell(XL_getType,getTypeStr," ARM_ERR: Accessor Type: string expected",C_result);
        string getType(CCSTringToSTLString(getTypeStr));

		CCString csoGetClass(GCGetTypeToClass(getType,csoId).c_str());

        exportFunc2Args< long, string > ourFunc(csoId,getType, ARMLOCAL_CSO_Get);

		/// call the general function
		fillXL_Result( csoGetClass, ourFunc, C_result, XL_result, PersistentInXL );
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_CSOGet_Common" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


///////////////////////////////////
/// version that takes into account 
/// previous creation of object
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_CSOCalculator_GetData(
	LPXLOPER XL_csoId,
	LPXLOPER XL_getType)
{
	ADD_LOG("Local_CSOCalculator_GetData");
	bool PersistentInXL = true;
	return Local_CSOGet_Common(
	    XL_csoId,
	    XL_getType,
        PersistentInXL );
}


__declspec(dllexport) LPXLOPER WINAPI Local_PXL_CSOCalculator_GetData(
	LPXLOPER XL_csoId,
	LPXLOPER XL_getType)
{
	ADD_LOG("Local_PXL_CSOCalculator_GetData");
	bool PersistentInXL = false;
	return Local_CSOGet_Common(
	    XL_csoId,
	    XL_getType,
        PersistentInXL );
}

/////////////////////////////////////////////////////////////
/// Set datas to a CSO calculator
/////////////////////////////////////////////////////////////
LPXLOPER Local_GenCalculatorSet_Common(
	LPXLOPER XL_calculatorId,
	LPXLOPER XL_dataId,
	LPXLOPER XL_setPortfolioType,
	LPXLOPER XL_MktDataKeys,
	bool isUpdated,
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

		long calculatorId;
		XL_GETOBJID( XL_calculatorId, calculatorId,	" ARM_ERR: Calculator: Object expected",C_result);

		long dataId;
		XL_GETOBJID( XL_dataId, dataId,	" ARM_ERR: Object expected",C_result);

		CCString setPortfolioTypeStr;
		XL_readStrCellWD(XL_setPortfolioType,setPortfolioTypeStr,"DEFAULT"," ARM_ERR: Portfolio Type: string expected",C_result);
		string setPortfolioType = CCSTringToSTLString(setPortfolioTypeStr);

		VECTOR<CCString> mktDataKeys;
		VECTOR<CCString> mktDataKeysDef(0);
		XL_readStrVectorWD(XL_MktDataKeys,mktDataKeys,mktDataKeysDef," ARM_ERR: Market datas keys: array of string expected",DOUBLE_TYPE,C_result);

		vector< string > mktDataKeysSTL(mktDataKeys.size());
		for(size_t i=0;i<mktDataKeys.size();++i)
        {
            mktDataKeys[i].toUpper();
			mktDataKeysSTL[i]=CCSTringToSTLString(mktDataKeys[i]);
        }

			exportFunc5Args<long, long, string, vector<string>, bool> ourFunc( 
								calculatorId,
								dataId,
								setPortfolioType,
								mktDataKeysSTL,
								isUpdated,
								ARMLOCAL_Calculator_Set);

        if(isUpdated)
        {
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
        else
        {
            /// General call through the functor with an object creation		    
		    fillXL_Result( LOCAL_GC_CRF_CLASS, ourFunc, C_result, XL_result, PersistentInXL );
        }
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_GenCalculatorSet_Common" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}

///////////////////////////////////
/// version that takes into account 
/// previous creation of object
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_GenCalculator_Update(
	LPXLOPER XL_calculatorId,
	LPXLOPER XL_dataId,
    LPXLOPER XL_setPortfolioType,
	LPXLOPER XL_mktDataKeys)
{
	ADD_LOG("Local_GenCalculator_Update");
	bool PersistentInXL = true;
	bool isUpdated = true;
	return Local_GenCalculatorSet_Common(
	    XL_calculatorId,
	    XL_dataId,
		XL_setPortfolioType,
		XL_mktDataKeys,
		isUpdated,
        PersistentInXL );
}

///////////////////////////////////
/// warnings
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_GetWarning_OnObject_Common(
	LPXLOPER XL_ObjId,
	bool persistent )
{
	ADD_LOG("Local_GetWarning_OnObject_Common");
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
		
		CCString C_ObjIdStr;
		XL_readStrCell( XL_ObjId, C_ObjIdStr, " ARM_ERR: Object Id: Object expected",C_result);
		long C_ObjIdId = LocalGetNumObjectId(C_ObjIdStr);

		if( ARMLOCAL_IsWarning_OnObject( C_ObjIdId, C_result ) )
		{
			exportFunc1Arg< long >  ourFunc( C_ObjIdId, ARMLOCAL_GetWarning_OnObject);

			/// call the general function
			fillXL_Result_withName( ourFunc, C_result, XL_result, persistent );
		}
		else
		{
			FreeCurCellErr ();
			XL_result.xltype  = xltypeStr;
			XL_result.val.str = XL_StrC2StrPascal ( "NO WARNING");
			XL_result.xltype |= xlbitDLLFree;
		}
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_GetWarning_OnObject_Common" )

	return (LPXLOPER)&XL_result;
}


__declspec(dllexport) LPXLOPER WINAPI Local_GetWarning_OnObject(
	LPXLOPER XL_ObjId )
{
	ADD_LOG("Local_GetWarning_OnObject");
	return Local_GetWarning_OnObject_Common( XL_ObjId, true );
}


__declspec(dllexport) LPXLOPER WINAPI Local_PXL_GetWarning_OnObject(
	LPXLOPER XL_ObjId )
{
	ADD_LOG("Local_PXL_GetWarning_OnObject");
	return Local_GetWarning_OnObject_Common( XL_ObjId, false );
}



///----------------------------------------------
///----------------------------------------------
///         Callable Range Accrual Spread Calculator
/// Inputs :
///----------------------------------------------
///----------------------------------------------

class craSpreadCalculatorFunc : public ARMResultLong2LongFunc
{
	public:
		craSpreadCalculatorFunc(int optionPfId,
								int refResetFreq,
								vector<string> mdmKeys,
								vector<string> calibParams,
								double payIndexMultValue,
								long payIndexMultId,
								long mktDataManagerId,
								vector<string> productsToPrice)
		:
		C_optionPfId(optionPfId),
		C_refResetFreq(refResetFreq),
		C_mdmKeys(mdmKeys),
		C_calibParams(calibParams),
		C_payIndexMultValue(payIndexMultValue),
		C_payIndexMultId(payIndexMultId),
		C_mktDataManagerId(mktDataManagerId),
		C_productsToPrice(productsToPrice)
		{
		};
		
		long operator()(ARM_result& result, long objId)
		{
			return	ARMLOCAL_CRASpreadCalculator_Create(C_optionPfId,
														C_refResetFreq,
														C_mdmKeys,
														C_calibParams,
														C_payIndexMultValue,
														C_payIndexMultId,
														C_mktDataManagerId,
														C_productsToPrice,
														result,
														objId);
		}


	private:
		int				C_optionPfId;
		int				C_refResetFreq;
		vector<string>	C_mdmKeys;
		vector<string>	C_calibParams;
		double			C_payIndexMultValue;
		long			C_payIndexMultId;
		long			C_mktDataManagerId;
		vector<string>	C_productsToPrice;
};


LPXLOPER Local_CRASpreadCalculator_Common(LPXLOPER XL_generalDatas,
										  LPXLOPER XL_callDatas,
										  LPXLOPER XL_fundDatas,
										  LPXLOPER XL_cpnDatas,
										  LPXLOPER XL_notionalCurve,
										  LPXLOPER XL_callFeesCurve,
										  LPXLOPER XL_fundSpreadCurve,
										  LPXLOPER XL_boostedFixCurve,
										  LPXLOPER XL_payIndexMultCurve,
										  LPXLOPER XL_barrierDownCurve,
										  LPXLOPER XL_barrierUpCurve,
										  LPXLOPER XL_refCoeff1Curve,
										  LPXLOPER XL_refCoeff2Curve,
										  bool IsDbleCorridor,
										  LPXLOPER XL_barrierCond1DownCurve,
										  LPXLOPER XL_barrierCond1UpCurve,
										  bool IsTripleRange,
										  bool IsVms,
										  LPXLOPER XL_ModelDatas,
										  LPXLOPER XL_mktDataManager,
										  LPXLOPER XL_productsToPrice,
										  LPXLOPER XL_localCalibFlags,
										  LPXLOPER XL_miscDatas,
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
	
		//General Datas
		//-------------
		VECTOR<CCString> generalDatas;
		XL_readStrVector (XL_generalDatas, generalDatas, "ARM_ERR: General datas: array of string expected", DOUBLE_TYPE, C_result);
		if (generalDatas.size() < 4)
		{
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}

		//Currency
		CCString ccyStr = generalDatas[0];
		string sCcy = CCSTringToSTLString(ccyStr);
		ARM_Currency ccy(sCcy.c_str());

		//Start & End
		CCString startDateXl = generalDatas[1];
		string sStartDate = CCSTringToSTLString(startDateXl);
		double startDate = atoi(sStartDate.c_str());
		CCString endDateXl = generalDatas[2];
		string sEndDate = CCSTringToSTLString(endDateXl);
		double endDate = atoi(sEndDate.c_str());
	    
		//PayReceive
		CCString payReceiveStr = generalDatas[3];
		long lPayReceive;
		if ((lPayReceive = ARM_ConvRecOrPay (payReceiveStr, C_result)) == ARM_DEFAULT_ERR)
		{
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}
		int payReceive = (int)lPayReceive;

		//CallDatas
		//---------
		VECTOR<CCString> callDatas;
		XL_readStrVector (XL_callDatas, callDatas, "ARM_ERR: Call datas: array of string expected", DOUBLE_TYPE, C_result);
		if (callDatas.size() < 3)
		{
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}

		/// Call frequency 
		CCString callFreqXl = callDatas[0];
		long lCallFreq;
		if ((lCallFreq = ARM_ConvFrequency (callFreqXl, C_result)) == ARM_DEFAULT_ERR)
		{
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}
		int callFreq = (int)lCallFreq;

		//CallNotice
		CCString callNoticeXl = callDatas[1];
		string sCallNotice = CCSTringToSTLString(callNoticeXl);
		int callNotice = atoi(sCallNotice.c_str());
		
		//CallCal
		CCString callCalXl = callDatas[2];
		string callCal = CCSTringToSTLString(callCalXl);

		//FundDatas
		//----------
		VECTOR<CCString> fundDatas;
		XL_readStrVector (XL_fundDatas, fundDatas," ARM_ERR: Fund datas: array of string expected",DOUBLE_TYPE,C_result);
		if (fundDatas.size() < 2)
		{
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}

		//Fund Freq
		CCString fundFreqXl = fundDatas[0];
		long lFundFreq;
		if ((lFundFreq = ARM_ConvFrequency (fundFreqXl, C_result)) == ARM_DEFAULT_ERR)
		{
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}
		int fundFreq = (int)lFundFreq;

		//Fund Day Count
		CCString fundDayCountXl = fundDatas[1];
		int fundDayCount = (int)ARM_ConvDayCount (fundDayCountXl);

		//Fund Leverage (funding = fundLeverage*Libor + fundSpread) :*** NOT USED AT THE MOMENT ***
		double unusedFundLeverage=1.0;
		if(fundDatas.size() > 2)
		{
			CCString fundLeverageXl = fundDatas[2];
			string sfundLeverage = CCSTringToSTLString(fundLeverageXl);
			unusedFundLeverage = atof(sfundLeverage.c_str());
		}

		ARM_Currency fundCcy(ccy);
		if(fundDatas.size() > 3)
		{
			string ccyName = CCSTringToSTLString(fundDatas[3]);
			fundCcy = ARM_Currency(ccyName.c_str());
		}

		//CpnDatas
		//------------
		VECTOR<CCString> cpnDatas;
		XL_readStrVector (XL_cpnDatas, cpnDatas," ARM_ERR: Cpn datas: array of string expected",DOUBLE_TYPE,C_result);
		if (cpnDatas.size() < 9)
		{
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}

		//Cpn Day Count
		CCString cpnDayCountXl = cpnDatas[0];
		int cpnDayCount = (int)ARM_ConvDayCount (cpnDayCountXl);

		//Cpn Pay Freq
		CCString cpnPayFreqXl = cpnDatas[1];
		long lCpnPayFreq;
		if ((lCpnPayFreq = ARM_ConvFrequency (cpnPayFreqXl, C_result)) == ARM_DEFAULT_ERR)
		{
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}
		int cpnPayFreq = (int)lCpnPayFreq;

		//Cpn Pay Gap
		CCString cpnResetGapXl = cpnDatas[2];
		string sCpnResetGap = CCSTringToSTLString(cpnResetGapXl);
		int cpnResetGap = atoi(sCpnResetGap.c_str());

		//Cpn Reset Cal
		CCString cpnResetCalXl = cpnDatas[3];
		string 	cpnResetCal = CCSTringToSTLString(cpnResetCalXl);

		//Cpn Pay Cal
		CCString cpnPayCalXl = cpnDatas[4];
		string 	cpnPayCal = CCSTringToSTLString(cpnPayCalXl);

		//CpnResetFreq
		CCString cpnResetFreqXl = cpnDatas[5];
		long lCpnResetFreq;
		if ((lCpnResetFreq = ARM_ConvFrequency (cpnResetFreqXl, C_result)) == ARM_DEFAULT_ERR)
		{
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}
		int cpnResetFreq = (int)lCpnResetFreq;

		//Cpn Reset Timing
		CCString cpnResetTimingXl = cpnDatas[6];
		int cpnResetTiming = (int)ARM_ConvPayResetRule(cpnResetTimingXl);

		int offsetIdx=7;
		int refCond1Index = K_FIXED;
		if(IsDbleCorridor)
		{
			// 1st condition Ref Index
			refCond1Index = -1000;  // required
			if((refCond1Index = ARM_ConvIrType(cpnDatas[offsetIdx]))==ARM_DEFAULT_ERR)
			{
				ARM_ARG_ERR();
				return (LPXLOPER)&XL_result;
			}
			++offsetIdx;
		}

		//Ref Index 1 (or 1st index of 2nd condition)
		int refIndex1 = ARM_ConvIrType(cpnDatas[offsetIdx]);
		++offsetIdx;

		//Ref Index 2 (or 2nd index of 2nd condition)
		int refIndex2 = ARM_ConvIrType(cpnDatas[offsetIdx]);
		++offsetIdx;

		// pay index (default = FIXED)
		int payIndex = K_FIXED;
		if (cpnDatas.size() > offsetIdx)
		{
			payIndex = ARM_ConvIrType(cpnDatas[offsetIdx]);
		}
		++offsetIdx;

		// pay index reset timing (default = K_ADVANCE)
		int payIndexResetTiming = K_ADVANCE;
		if (cpnDatas.size() > offsetIdx)
		{
			CCString payIndexResetTimingXl = cpnDatas[offsetIdx];
			payIndexResetTiming = (int)ARM_ConvPayResetRule(payIndexResetTimingXl);
		}

		//NotionalValue or NotionalCurve
		double notional;
		CCString notionalStr;
		long notionalId;
		XL_readStrOrNumCell(XL_notionalCurve, notionalStr, notional, notionalId,
			   " ARM_ERR: notional: numeric or refValue Id expected",C_result);	

		if (notionalId == XL_TYPE_STRING)
			notionalId = LocalGetNumObjectId(notionalStr);
		else
			notionalId = ARM_NULL_OBJECT;
		
		long fundNotionaldble = ARM_NULL_OBJECT ;
		if(fundDatas.size() > 3)
		{
			string fundNotionalstr = CCSTringToSTLString(fundDatas[3]);
			fundNotionaldble = atof(fundNotionalstr.c_str());
		}
		else
			fundNotionaldble = notional; ///bad bad 

		//Call Fees Curve //To change
		CCString callFeesStr;
		long callFeesId;
		XL_readStrCell(XL_callFeesCurve, callFeesStr," ARM_ERR: Call Fees: String Id expected", C_result);
		
		callFeesId = LocalGetNumObjectId(callFeesStr);

		//FundSpreadValue or FundSpreadCurve
		double fundSpread;
		CCString fundSpreadStr;
		long fundSpreadId;
		XL_readStrOrNumCell(XL_fundSpreadCurve, fundSpreadStr, fundSpread, fundSpreadId,
			   " ARM_ERR: fundSpread: numeric or refValue Id expected",C_result);	
		
		if (fundSpreadId == XL_TYPE_STRING)
			fundSpreadId = LocalGetNumObjectId(fundSpreadStr);
		else
			fundSpreadId = ARM_NULL_OBJECT;

		//PayIndexMult Value or PayIndexMultCurve
		double payIndexMult;
		CCString payIndexMultStr;
		long payIndexMultId;
		XL_readStrOrNumCell(XL_payIndexMultCurve, payIndexMultStr, payIndexMult, payIndexMultId,
			   " ARM_ERR: PayIndexMult: numeric or refValue Id expected",C_result);	
		
		if (payIndexMultId == XL_TYPE_STRING)
			payIndexMultId = LocalGetNumObjectId(payIndexMultStr);
		else
			payIndexMultId = ARM_NULL_OBJECT;

		double boostedFix;
		long boostedFixId;
		long boostedFix2Id=ARM_NULL_OBJECT;
		long boostedFix3Id=ARM_NULL_OBJECT;
		double bDown;
		long bDownId;
		long bDown2Id=ARM_NULL_OBJECT;
		long bDown3Id=ARM_NULL_OBJECT;
		double bUp;
		long bUpId;
		long bUp2Id=ARM_NULL_OBJECT;
		long bUp3Id=ARM_NULL_OBJECT;

		if (IsTripleRange)
		{
			vector<CCString> C_boostedFixIds;
			XL_readStrVector (XL_boostedFixCurve,C_boostedFixIds," ARM_ERR: boost fix: array of object expected",DOUBLE_TYPE,C_result);
			boostedFixId = LocalGetNumObjectId(C_boostedFixIds[0]); 
			boostedFix2Id = LocalGetNumObjectId(C_boostedFixIds[1]); 
			boostedFix3Id = LocalGetNumObjectId(C_boostedFixIds[2]); 

			vector<CCString> C_bDownIds;
			XL_readStrVector (XL_barrierDownCurve,C_bDownIds," ARM_ERR: down barrier: array of object expected",DOUBLE_TYPE,C_result);
			bDownId = LocalGetNumObjectId(C_bDownIds[0]); 
			bDown2Id = LocalGetNumObjectId(C_bDownIds[1]); 
			bDown3Id = LocalGetNumObjectId(C_bDownIds[2]); 

			vector<CCString> C_bUpIds;
			XL_readStrVector (XL_barrierUpCurve,C_bUpIds," ARM_ERR: up barrier: array of object expected",DOUBLE_TYPE,C_result);
			bUpId = LocalGetNumObjectId(C_bUpIds[0]); 
			bUp2Id = LocalGetNumObjectId(C_bUpIds[1]); 
			bUp3Id = LocalGetNumObjectId(C_bUpIds[2]); 

		}
		else
		{
			//BoostedFixValue or BoostedFixCurve
			CCString boostedFixStr;
			XL_readStrOrNumCell(XL_boostedFixCurve, boostedFixStr, boostedFix, boostedFixId,
				   " ARM_ERR: boosted fix rate: numeric or refValue Id expected",C_result);	
			
			if (boostedFixId == XL_TYPE_STRING)
				boostedFixId = LocalGetNumObjectId(boostedFixStr);
			else
				boostedFixId = ARM_NULL_OBJECT;

			//BarrierDownValue or BarrierDownCurve
			CCString bDownStr;
			XL_readStrOrNumCell(XL_barrierDownCurve, bDownStr, bDown, bDownId,
				   " ARM_ERR: BarrierDown: numeric or refValue Id expected",C_result);	
			
			if (bDownId == XL_TYPE_STRING)
				bDownId = LocalGetNumObjectId(bDownStr);
			else
				bDownId = ARM_NULL_OBJECT;

			//BarrierUpValue or BarrierUpCurve
			CCString bUpStr;
			XL_readStrOrNumCell(XL_barrierUpCurve, bUpStr, bUp, bUpId,
				   " ARM_ERR: BarrierUp: numeric or refValue Id expected",C_result);	
			
			if (bUpId == XL_TYPE_STRING)
				bUpId = LocalGetNumObjectId(bUpStr);
			else
				bUpId = ARM_NULL_OBJECT;
		}

		//RefCoeff1Value or RefCoeff1Curve
		double coeff1, coeff1Def=1.0;
		XL_readNumCellWD(XL_refCoeff1Curve, coeff1, coeff1Def,
			   " ARM_ERR: Ref Coeff 1: numeric expected",C_result);	
		

		double coeff2 = 1.;
		vector< string > TenorVec;
		if (IsVms)
		{
			int k;
			VECTOR<CCString> TenorVecXL;
			XL_readStrVector(XL_refCoeff2Curve,TenorVecXL," ARM_ERR: array of string expected for VMS",DOUBLE_TYPE,C_result);
			TenorVec.resize(TenorVecXL.size());
			for (k=0;k<TenorVecXL.size();++k)
				TenorVec[k]=CCSTringToSTLString(TenorVecXL[k]);
			coeff2 = 0.;
		}
		else
		{
			//RefCoeff2Value or RefCoeff2Curve
			double coeff2Def=1.0;
			XL_readNumCellWD(XL_refCoeff2Curve, coeff2, coeff2Def,
				   " ARM_ERR: Ref Coeff 2: numeric expected",C_result);
		}

		double bCond1Down=0.0;
		long bCond1DownId=ARM_NULL_OBJECT;
		double bCond1Up=0.0;
		long bCond1UpId=ARM_NULL_OBJECT;

		/// Local calib flags
		VECTOR<CCString> LocalCalibFlagsXL;
		VECTOR<CCString> LocalCalibFlagsDefault(5,"N");
		XL_readStrVectorWD(XL_localCalibFlags,LocalCalibFlagsXL,LocalCalibFlagsDefault," ARM_ERR: Local Calib Flag: array of string expected",DOUBLE_TYPE,C_result);
		vector< string > localCalibFlags(5,"N");
		size_t flagOffset=0;
		if(IsDbleCorridor)
		{
			//BarrierDownValue or BarrierDownCurve on the 1st condition
			CCString bCond1DownStr;
			XL_readStrOrNumCell(XL_barrierCond1DownCurve, bCond1DownStr, bCond1Down, bCond1DownId,
				   " ARM_ERR: BarrierDown on 1st Condition : numeric or refValue Id expected",C_result);	
			
			if (bCond1DownId == XL_TYPE_STRING)
				bCond1DownId = LocalGetNumObjectId(bCond1DownStr);
			else
				bCond1DownId = ARM_NULL_OBJECT;

			//BarrierUpValue or BarrierUpCurve on the 1st condition
			CCString bCond1UpStr;
			XL_readStrOrNumCell(XL_barrierCond1UpCurve, bCond1UpStr, bCond1Up, bCond1UpId,
				   " ARM_ERR: BarrierUp on 1st condition : numeric or refValue Id expected",C_result);	
			
			if (bCond1UpId == XL_TYPE_STRING)
				bCond1UpId = LocalGetNumObjectId(bCond1UpStr);
			else
				bCond1UpId = ARM_NULL_OBJECT;

			/// Correl Unsqueezer
			localCalibFlags[0]=CCSTringToSTLString(LocalCalibFlagsXL[1]);

			flagOffset = 1;
		}

		/// Vol Unsqueezer
		localCalibFlags[1]=CCSTringToSTLString(LocalCalibFlagsXL[0]);

		/// PV adjuster
		localCalibFlags[2]=CCSTringToSTLString(LocalCalibFlagsXL[flagOffset+1]);

		/// Bootstrap Optimizer
		localCalibFlags[3]=CCSTringToSTLString(LocalCalibFlagsXL[flagOffset+2]);

		/// Switch to old calib
		localCalibFlags[4]=CCSTringToSTLString(LocalCalibFlagsXL[flagOffset+3]);


		int i;
		// Calib Flags
		VECTOR<CCString> CalibFlagsXL;
		VECTOR<CCString> CalibFlags_default(1);
		CalibFlags_default[0]= CCString("Y");
		
		XL_readStrVectorWD (XL_ModelDatas,CalibFlagsXL,CalibFlags_default," ARM_ERR: Calib Mode: array of string expected",DOUBLE_TYPE,C_result);
		vector< string > CalibFlags(CalibFlagsXL.size());
		for(i=0;i<CalibFlagsXL.size();++i)
			CalibFlags[i]=CCSTringToSTLString(CalibFlagsXL[i]);
		
		if( CalibFlags.size() > 8 )
			throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"Bad Calib Flags size ");


		//MktDataManager: persistent object
		VECTOR<CCString> mktDataManager;
		XL_readStrVector(XL_mktDataManager,mktDataManager," ARM_ERR: Market datas: array of string expected",DOUBLE_TYPE,C_result);
        
		if (mktDataManager.size()<1)
		{
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}

		long mktDataManagerId = LocalGetNumObjectId (mktDataManager[0]);
		vector<string> mdmKeys(mktDataManager.size()-1);		
		for (i=1; i<mktDataManager.size(); ++i)
			mdmKeys[i-1] = CCSTringToSTLString(mktDataManager[i]);
	
		//Pricing flags
		VECTOR<CCString> pricingFlags;

		XL_readStrVector(XL_productsToPrice,pricingFlags," ARM_ERR: Pricing flags: array of string expected",DOUBLE_TYPE,C_result);
		vector< string > productsToPrice(pricingFlags.size());
		for (i=0; i<pricingFlags.size(); i++)
			productsToPrice[i] = CCSTringToSTLString(pricingFlags[i]);

		/// Miscenalleous Datas : ExerProbas + OptimResetDatas
		VECTOR<double> miscDatas;
		VECTOR<double> miscDatasDefault(3,0.0);
		XL_readNumVectorWD(XL_miscDatas,miscDatas,miscDatasDefault," ARM_ERR: MiscDatas: array of double expected",C_result);

		size_t nbMisc = miscDatas.size();
		VECTOR<double> exerProbas(2);
		if(nbMisc>0) exerProbas[0] = miscDatas[0];
		else exerProbas[0]=0.0;
		if(nbMisc>1) exerProbas[1] = miscDatas[1];
		else exerProbas[1]=0.0;
		VECTOR<double> optimResetData(nbMisc-2>0 ? nbMisc-2 : 0);
		for(i=2;i<nbMisc;++i) optimResetData[i-2]=miscDatas[i];

	exportFunc56Args<
			ARM_Currency,				// 1
			ARM_Currency,				// 1'
			double,						// 2
			double,						// 3
			int,						// 4
			int,						// 5
			int,						// 6
			string,						// 7
			int,						// 8
			int,						// 9
			int,						// 10
			int,						// 11
			string,						// 12
			string,						// 13
			int,						// 14
			int,						// 15
			int,						// 16
			int,						// 17
			int,						// 18
			int,						// 19
			int,						// 20
			double,						// 21
			long,						// 22
			double,                      // 22'
			long,						// 23
			double,						// 24
			long,						// 25
			double,						// 26
			long,						// 27
			double,						// 28
			long,						// 29
			double,						// 30
			long,						// 31
			double,						// 32
			long,						// 33
			double,						// 34
			double,						// 35
			int,						// 36
			double,						// 37
			long,						// 39
			double,						// 39
			long,						// 40
			long,						// 41
			long,						// 42
			long,						// 43
			long,						// 44
			long,						// 45
			long,						// 46
			vector<string>,				// 47
			vector<string>,				// 48
			vector<string>,				// 49
			long,						// 50
			vector<string>,				// 51
			vector<string>,				// 52
			vector<double>,				// 53
			vector<double> >			// 54	
			ourFunc(	
			ccy,						// 1
			fundCcy,                    // 1'
			startDate,					// 2
			endDate,					// 3
			payReceive,					// 4
			callFreq,					// 5
			callNotice,					// 6
			callCal,					// 7
			fundFreq,					// 8
			fundDayCount,				// 9
			cpnDayCount,				// 10
			cpnPayFreq,					// 11	
			cpnResetCal,				// 12	
			cpnPayCal,					// 13	
			cpnResetFreq,				// 14	
			cpnResetTiming,				// 15
			cpnResetGap,				// 16
			refIndex1,					// 17
			refIndex2,					// 18
			payIndex,					// 19
			payIndexResetTiming,		// 20
			notional,					// 21
			notionalId,					// 22
			fundNotionaldble,           // 22'           
			callFeesId,					// 23
			fundSpread,					// 24
			fundSpreadId,				// 25
			boostedFix,					// 26
			boostedFixId,				// 27
			payIndexMult,				// 28
			payIndexMultId,				// 29
			bDown,						// 30
			bDownId,					// 31
			bUp,						// 32
			bUpId,						// 33
			coeff1,						// 34
			coeff2,						// 35
			refCond1Index,				// 36
			bCond1Down,					// 37
			bCond1DownId,				// 38
			bCond1Up,					// 39
			bCond1UpId,					// 40
			boostedFix2Id,				// 41
			bDown2Id,					// 42
			bUp2Id,						// 43
			boostedFix3Id,				// 44
			bDown3Id,					// 45
			bUp3Id,						// 46
			TenorVec,					// 47
			CalibFlags,					// 48
			mdmKeys,					// 49
			mktDataManagerId,			// 50
			productsToPrice,			// 51
			localCalibFlags,			// 52
			optimResetData,				// 53
			exerProbas,					// 54
			ARMLOCAL_CRASpreadCalculator_Create);			


		/// call the general function
		fillXL_Result( LOCAL_GC_CCSO_CLASS, ourFunc, C_result, XL_result, PersistentInXL );
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_CRASpreadCalculator_Common" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


LPXLOPER Local_CRASpreadCalculator_Common(	LPXLOPER XL_optionPortfolio,
											LPXLOPER XL_ModelDatas,
											LPXLOPER XL_payIndexMult,
											LPXLOPER XL_mktDataManager,
											LPXLOPER XL_productsToPrice,
											LPXLOPER XL_refResetFreq,
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
	
		CCString C_optionPortfolio;
		XL_readStrCell( XL_optionPortfolio, C_optionPortfolio, " ARM_ERR: Option Portfolio Id: Object expected",C_result);
		long C_optionPortfolioId = LocalGetNumObjectId(C_optionPortfolio);

/*		//ModelType
		VECTOR<CCString> ModelDatas;
		VECTOR<CCString> ModelDatasDef(1);
		ModelDatasDef[0] = "Y";
		XL_readStrVectorWD(XL_ModelDatas,ModelDatas,ModelDatasDef," ARM_ERR: Model Datas: array of string expected",DOUBLE_TYPE, C_result);
		string modelType;
		string calibType;
		string calibStrikeType;
		string vnsPricingMethod = "MONEYNESS"; /// Default value

		if (ModelDatas.size() < 3)
			throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"Invalid Model Data size");

		modelType = CCSTringToSTLString(ModelDatas[0]);
		calibType = CCSTringToSTLString(ModelDatas[1]);
		calibStrikeType = CCSTringToSTLString(ModelDatas[2]);

		if (ModelDatas.size() > 3)
			vnsPricingMethod = CCSTringToSTLString(ModelDatas[3]);*/

		VECTOR<CCString> CalibFlagsXL;
		VECTOR<CCString> CalibFlags_default(1);
		CalibFlags_default[0]= CCString("Y");
		/*CalibFlags_default[1]= CCString("DIAG");
		CalibFlags_default[2]= CCString("ATM");
		CalibFlags_default[2]= CCString("MONEYNESS");*/
		
		XL_readStrVectorWD (XL_ModelDatas,CalibFlagsXL,CalibFlags_default," ARM_ERR: Calib Mode: array of string expected",DOUBLE_TYPE,C_result);
		vector< string > CalibFlags(CalibFlagsXL.size());
		for(int i=0;i<CalibFlagsXL.size();++i)
			CalibFlags[i]=CCSTringToSTLString(CalibFlagsXL[i]);
		
		if( CalibFlags.size() > 8 )
			throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"Bad Calib Flags size ");

		//PayIndexMultValue or PayIndexMultId
		double payIndexMultValue;
		CCString payIndexMultStr;
		long payIndexMultId;
		XL_readStrOrNumCell(XL_payIndexMult, payIndexMultStr, payIndexMultValue, payIndexMultId,
			   " ARM_ERR: PayIndexMult: numeric or refValue Id expected",C_result);	
		
		if (payIndexMultId == XL_TYPE_STRING)
			payIndexMultId = LocalGetNumObjectId(payIndexMultStr);
		else
			payIndexMultId = ARM_NULL_OBJECT;

		//MktDataManager: persistent object
		VECTOR<CCString> mktDataManager;
		XL_readStrVector(XL_mktDataManager,mktDataManager," ARM_ERR: Market datas: array of string expected",DOUBLE_TYPE,C_result);
        
		if (mktDataManager.size()<4)
		{
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}

		long mktDataManagerId = LocalGetNumObjectId (mktDataManager[0]);
		vector<string> mdmKeys(mktDataManager.size()-1);		
		for (i=1; i<mktDataManager.size(); ++i)
			mdmKeys[i-1] = CCSTringToSTLString(mktDataManager[i]);
	
		//Pricing flags
		VECTOR<CCString> pricingFlags;

		XL_readStrVector(XL_productsToPrice,pricingFlags," ARM_ERR: Pricing flags: array of string expected",DOUBLE_TYPE,C_result);
        if (pricingFlags.size()>7)
		{
			throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"maximum: 7 products for classic deal and 5 for local deal");
		}

		vector< string > productsToPrice(pricingFlags.size());
		for (i=0; i<pricingFlags.size(); i++)
			productsToPrice[i] = CCSTringToSTLString(pricingFlags[i]);

		//Reset frequency of reference index
		CCString C_refResetFreq;
		long refResetFreq;
		XL_readStrCellWD(XL_refResetFreq, C_refResetFreq, "-1", "ARM_ERR: Reset Frequency: string expected", C_result);
		if ( C_refResetFreq == "-1" )
		{
		   refResetFreq = K_DEF_FREQ;
		}
		else
		{
			if ((refResetFreq = ARM_ConvFrequency(C_refResetFreq, C_result)) == ARM_DEFAULT_ERR)
			{
				throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"Reset frequency : invalid frequency");
			}
		}

		//Functor call
		craSpreadCalculatorFunc ourFunc(C_optionPortfolioId,
										refResetFreq,
										mdmKeys,
										CalibFlags,
										payIndexMultValue,
										payIndexMultId,
										mktDataManagerId,
										productsToPrice);

		/// call the general function
		fillXL_Result( LOCAL_GC_CCSO_CLASS, ourFunc, C_result, XL_result, PersistentInXL );
		
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_CRASpreadCalculator_Common" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}

__declspec(dllexport) LPXLOPER WINAPI Local_CRASpreadCalculator_CreateWithoutMarketData(LPXLOPER XL_generalDatas,
																					   LPXLOPER XL_callDatas,
																					   LPXLOPER XL_fundDatas,
																					   LPXLOPER XL_cpnDatas,
																					   LPXLOPER XL_notionalCurve,
																					   LPXLOPER XL_callFeesCurve,
																					   LPXLOPER XL_fundSpreadCurve,
																					   LPXLOPER XL_boostedFixCurve,
																					   LPXLOPER XL_payIndexMultCurve,
																					   LPXLOPER XL_barrierDownCurve,
																					   LPXLOPER XL_barrierUpCurve,
																					   LPXLOPER XL_refCoeff1Curve,
																					   LPXLOPER XL_refCoeff2Curve)
{
	ADD_LOG("Local_CRASpreadCalculator_CreateWithoutMarketData");
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
	
		//General Datas
		//-------------
		VECTOR<CCString> generalDatas;
		XL_readStrVector (XL_generalDatas, generalDatas, "ARM_ERR: General datas: array of string expected", DOUBLE_TYPE, C_result);
		if (generalDatas.size() < 4)
		{
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}

		//Currency
		CCString ccyStr = generalDatas[0];
		string sCcy = CCSTringToSTLString(ccyStr);
		ARM_Currency ccy(sCcy.c_str());

		//Start & End
		CCString startDateXl = generalDatas[1];
		string sStartDate = CCSTringToSTLString(startDateXl);
		double startDate = atoi(sStartDate.c_str());
		CCString endDateXl = generalDatas[2];
		string sEndDate = CCSTringToSTLString(endDateXl);
		double endDate = atoi(sEndDate.c_str());
	    
		//PayReceive
		CCString payReceiveStr = generalDatas[3];
		long lPayReceive;
		if ((lPayReceive = ARM_ConvRecOrPay (payReceiveStr, C_result)) == ARM_DEFAULT_ERR)
		{
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}
		int payReceive = (int)lPayReceive;

		//CallDatas
		//---------
		VECTOR<CCString> callDatas;
		XL_readStrVector (XL_callDatas, callDatas, "ARM_ERR: Call datas: array of string expected", DOUBLE_TYPE, C_result);
		if (callDatas.size() < 3)
		{
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}

		/// Call frequency 
		CCString callFreqXl = callDatas[0];
		long lCallFreq;
		if ((lCallFreq = ARM_ConvFrequency (callFreqXl, C_result)) == ARM_DEFAULT_ERR)
		{
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}
		int callFreq = (int)lCallFreq;

		//CallNotice
		CCString callNoticeXl = callDatas[1];
		string sCallNotice = CCSTringToSTLString(callNoticeXl);
		int callNotice = atoi(sCallNotice.c_str());
		
		//CallCal
		CCString callCalXl = callDatas[2];
		string callCal = CCSTringToSTLString(callCalXl);

		//FundDatas
		//----------
		VECTOR<CCString> fundDatas;
		XL_readStrVector (XL_fundDatas, fundDatas," ARM_ERR: Fund datas: array of string expected",DOUBLE_TYPE,C_result);
		if (fundDatas.size() < 2)
		{
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}

		//Fund Freq
		CCString fundFreqXl = fundDatas[0];
		long lFundFreq;
		if ((lFundFreq = ARM_ConvFrequency (fundFreqXl, C_result)) == ARM_DEFAULT_ERR)
		{
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}
		int fundFreq = (int)lFundFreq;

		//Fund Day Count
		CCString fundDayCountXl = fundDatas[1];
		int fundDayCount = (int)ARM_ConvDayCount (fundDayCountXl);

		//CpnDatas
		//------------
		VECTOR<CCString> cpnDatas;
		XL_readStrVector (XL_cpnDatas, cpnDatas," ARM_ERR: Cpn datas: array of string expected",DOUBLE_TYPE,C_result);
		if (cpnDatas.size() < 9)
		{
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}

		//Cpn Day Count
		CCString cpnDayCountXl = cpnDatas[0];
		int cpnDayCount = (int)ARM_ConvDayCount (cpnDayCountXl);

		//Cpn Pay Freq
		CCString cpnPayFreqXl = cpnDatas[1];
		long lCpnPayFreq;
		if ((lCpnPayFreq = ARM_ConvFrequency (cpnPayFreqXl, C_result)) == ARM_DEFAULT_ERR)
		{
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}
		int cpnPayFreq = (int)lCpnPayFreq;

		//Cpn Pay Gap
		CCString cpnResetGapXl = cpnDatas[2];
		string sCpnResetGap = CCSTringToSTLString(cpnResetGapXl);
		int cpnResetGap = atoi(sCpnResetGap.c_str());

		//Cpn Reset Cal
		CCString cpnResetCalXl = cpnDatas[3];
		string 	cpnResetCal = CCSTringToSTLString(cpnResetCalXl);

		//Cpn Pay Cal
		CCString cpnPayCalXl = cpnDatas[4];
		string 	cpnPayCal = CCSTringToSTLString(cpnPayCalXl);

		//CpnResetFreq
		CCString cpnResetFreqXl = cpnDatas[5];
		long lCpnResetFreq;
		if ((lCpnResetFreq = ARM_ConvFrequency (cpnResetFreqXl, C_result)) == ARM_DEFAULT_ERR)
		{
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}
		int cpnResetFreq = (int)lCpnResetFreq;

		//Cpn Reset Timing
		CCString cpnResetTimingXl = cpnDatas[6];
		int cpnResetTiming = (int)ARM_ConvPayResetRule(cpnResetTimingXl);

		int offsetIdx=7;

		//Ref Index 1 (or 1st index of 2nd condition)
		int refIndex1 = ARM_ConvIrType(cpnDatas[offsetIdx]);
		++offsetIdx;

		//Ref Index 2 (or 2nd index of 2nd condition)
		int refIndex2 = ARM_ConvIrType(cpnDatas[offsetIdx]);
		++offsetIdx;

		// pay index (default = FIXED)
		int payIndex = K_FIXED;
		if (cpnDatas.size() > offsetIdx)
		{
			payIndex = ARM_ConvIrType(cpnDatas[offsetIdx]);
		}
		++offsetIdx;

		// pay index reset timing (default = K_ADVANCE)
		int payIndexResetTiming = K_ADVANCE;
		if (cpnDatas.size() > offsetIdx)
		{
			CCString payIndexResetTimingXl = cpnDatas[offsetIdx];
			payIndexResetTiming = (int)ARM_ConvPayResetRule(payIndexResetTimingXl);
		}


		//NotionalValue or NotionalCurve
		double notional;
		CCString notionalStr;
		long notionalId;
		XL_readStrOrNumCell(XL_notionalCurve, notionalStr, notional, notionalId,
			   " ARM_ERR: notional: numeric or refValue Id expected",C_result);	

		if (notionalId == XL_TYPE_STRING)
			notionalId = LocalGetNumObjectId(notionalStr);
		else
			notionalId = ARM_NULL_OBJECT;

		//Call Fees Curve //To change
		CCString callFeesStr;
		long callFeesId;
		XL_readStrCell(XL_callFeesCurve, callFeesStr," ARM_ERR: Call Fees: String Id expected", C_result);
		
		callFeesId = LocalGetNumObjectId(callFeesStr);

		//FundSpreadValue or FundSpreadCurve
		double fundSpread;
		CCString fundSpreadStr;
		long fundSpreadId;
		XL_readStrOrNumCell(XL_fundSpreadCurve, fundSpreadStr, fundSpread, fundSpreadId,
			   " ARM_ERR: fundSpread: numeric or refValue Id expected",C_result);	
		
		if (fundSpreadId == XL_TYPE_STRING)
			fundSpreadId = LocalGetNumObjectId(fundSpreadStr);
		else
			fundSpreadId = ARM_NULL_OBJECT;

		//BoostedFixValue or BoostedFixCurve
		double boostedFix;
		CCString boostedFixStr;
		long boostedFixId;
		XL_readStrOrNumCell(XL_boostedFixCurve, boostedFixStr, boostedFix, boostedFixId,
			   " ARM_ERR: boosted fix rate: numeric or refValue Id expected",C_result);	
		
		if (boostedFixId == XL_TYPE_STRING)
			boostedFixId = LocalGetNumObjectId(boostedFixStr);
		else
			boostedFixId = ARM_NULL_OBJECT;

		//PayIndexMult Value or PayIndexMultCurve
		double payIndexMult;
		CCString payIndexMultStr;
		long payIndexMultId;
		XL_readStrOrNumCell(XL_payIndexMultCurve, payIndexMultStr, payIndexMult, payIndexMultId,
			   " ARM_ERR: PayIndexMult: numeric or refValue Id expected",C_result);	
		
		if (payIndexMultId == XL_TYPE_STRING)
			payIndexMultId = LocalGetNumObjectId(payIndexMultStr);
		else
			payIndexMultId = ARM_NULL_OBJECT;

		//BarrierDownValue or BarrierDownCurve
		double bDown;
		CCString bDownStr;
		long bDownId;
		XL_readStrOrNumCell(XL_barrierDownCurve, bDownStr, bDown, bDownId,
			   " ARM_ERR: BarrierDown: numeric or refValue Id expected",C_result);	
		
		if (bDownId == XL_TYPE_STRING)
			bDownId = LocalGetNumObjectId(bDownStr);
		else
			bDownId = ARM_NULL_OBJECT;

		//BarrierUpValue or BarrierUpCurve
		double bUp;
		CCString bUpStr;
		long bUpId;
		XL_readStrOrNumCell(XL_barrierUpCurve, bUpStr, bUp, bUpId,
			   " ARM_ERR: BarrierUp: numeric or refValue Id expected",C_result);	
		
		if (bUpId == XL_TYPE_STRING)
			bUpId = LocalGetNumObjectId(bUpStr);
		else
			bUpId = ARM_NULL_OBJECT;

		//RefCoeff1Value or RefCoeff1Curve
		double coeff1, coeff1Def=1.0;
		XL_readNumCellWD(XL_refCoeff1Curve, coeff1, coeff1Def,
			   " ARM_ERR: Ref Coeff 1: numeric expected",C_result);	
		

		//RefCoeff2Value or RefCoeff2Curve
		double coeff2, coeff2Def=1.0;
		XL_readNumCellWD(XL_refCoeff2Curve, coeff2, coeff2Def,
			   " ARM_ERR: Ref Coeff 2: numeric expected",C_result);
	

	exportFunc35Args<
			ARM_Currency,				// 1
			double,						// 2
			double,						// 3
			int,						// 4
			int,						// 5
			int,						// 6
			string,						// 7
			int,						// 8
			int,						// 9
			int,						// 10
			int,						// 11
			string,						// 12
			string,						// 13
			int,						// 14
			int,						// 15
			int,						// 16
			int,						// 17
			int,						// 18
			int,						// 19
			int,						// 20
			double,						// 21
			long,						// 22
			long,						// 23
			double,						// 24
			long,						// 25
			double,						// 26
			long,						// 27
			double,						// 28
			long,						// 29
			double,						// 30
			long,						// 31
			double,						// 32
			long,						// 33
			double,						// 34
			double >					// 35

			ourFunc(	
			ccy,						// 1
			startDate,					// 2
			endDate,					// 3
			payReceive,					// 4
			callFreq,					// 5
			callNotice,					// 6
			callCal,					// 7
			fundFreq,					// 8
			fundDayCount,				// 9
			cpnDayCount,				// 10
			cpnPayFreq,					// 11	
			cpnResetCal,				// 12	
			cpnPayCal,					// 13	
			cpnResetFreq,				// 14	
			cpnResetTiming,				// 15
			cpnResetGap,				// 16
			refIndex1,					// 17
			refIndex2,					// 18
			payIndex,					// 19
			payIndexResetTiming,		// 20
			notional,					// 21
			notionalId,					// 22
			callFeesId,					// 23
			fundSpread,					// 24
			fundSpreadId,				// 25
			boostedFix,					// 26
			boostedFixId,				// 27
			payIndexMult,				// 28
			payIndexMultId,				// 29
			bDown,						// 30
			bDownId,					// 31
			bUp,						// 32
			bUpId,						// 33
			coeff1,						// 34
			coeff2,						// 35
			ARMLOCAL_CRASpreadCalculator_Create);			

		bool PersistentInXL = true;

		/// call the general function
		fillXL_Result( LOCAL_GC_CCSO_CLASS, ourFunc, C_result, XL_result, PersistentInXL );
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_CRASpreadCalculator_CreateWithoutMarketData" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}

class basicCraSpreadCalculatorFunc : public ARMResultLong2LongFunc
{
	public:
		basicCraSpreadCalculatorFunc(double asOfDate,
									int securityId,
									int refResetFreq,
									double payIndexMultValue,
									long payIndexMultId)
		:
		C_asOfDate(asOfDate),
		C_securityId(securityId),
		C_refResetFreq(refResetFreq),
		C_payIndexMultValue(payIndexMultValue),
		C_payIndexMultId(payIndexMultId)
		{
		};
		
		long operator()(ARM_result& result, long objId)
		{
			return	ARMLOCAL_CRASpreadCalculator_Create(C_asOfDate,
														C_securityId,
														C_refResetFreq,
														C_payIndexMultValue,
														C_payIndexMultId,
														result,
														objId);
		}

	private:
		double			C_asOfDate;
		int				C_securityId;
		int				C_refResetFreq;
		double			C_payIndexMultValue;
		long			C_payIndexMultId;
};

LPXLOPER Local_CRASpreadCalculator_Common(	LPXLOPER XL_asOfDate,
										    LPXLOPER XL_security,
											LPXLOPER XL_payIndexMult,
											LPXLOPER XL_refResetFreq,
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
	
		double asOfDate;
		XL_readNumCell(XL_asOfDate,asOfDate," ARM_ERR: as of date: date expected",C_result);

		CCString C_security;
		XL_readStrCell( XL_security, C_security, " ARM_ERR: Security Id: Object expected (Option PF or Swaption)",C_result);
		long C_securityId = LocalGetNumObjectId(C_security);

		//PayIndexMultValue or PayIndexMultId
		double payIndexMultValue;
		CCString payIndexMultStr;
		long payIndexMultId;
		XL_readStrOrNumCell(XL_payIndexMult, payIndexMultStr, payIndexMultValue, payIndexMultId,
			   " ARM_ERR: PayIndexMult: numeric or refValue Id expected",C_result);	
		
		if (payIndexMultId == XL_TYPE_STRING)
			payIndexMultId = LocalGetNumObjectId(payIndexMultStr);
		else
			payIndexMultId = ARM_NULL_OBJECT;

		//Reset frequency of reference index
		CCString C_refResetFreq;
		long refResetFreq;
		XL_readStrCellWD(XL_refResetFreq, C_refResetFreq, "-1", "ARM_ERR: Reset Frequency: string expected", C_result);
		if ( C_refResetFreq == "-1" )
		{
		   refResetFreq = K_DEF_FREQ;
		}
		else
		{
			if ((refResetFreq = ARM_ConvFrequency(C_refResetFreq, C_result)) == ARM_DEFAULT_ERR)
			{
				throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"Reset frequency : invalid frequency");
			}
		}

		//Functor call
		basicCraSpreadCalculatorFunc ourFunc(asOfDate,
											C_securityId,
											refResetFreq,
											payIndexMultValue,
											payIndexMultId);

		/// call the general function
		fillXL_Result( LOCAL_GC_CCSO_CLASS, ourFunc, C_result, XL_result, PersistentInXL );
		
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_CRASpreadCalculator_Common" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}

///----------------------------------------------
///----------------------------------------------
///             CRA Spread Calculator Accessor
/// Inputs :
///     
///----------------------------------------------
///----------------------------------------------

/////////////////////////////////////////////////////////////
/// central function that does the creation of the XL function
/////////////////////////////////////////////////////////////
LPXLOPER Local_CRASpreadGet_Common(
	LPXLOPER XL_craId,
	LPXLOPER XL_getType,
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
		XL_GETOBJID( XL_craId, craId, " ARM_ERR: CRA Calculator: Object expected", C_result);

		CCString getTypeStr;
		XL_readStrCell(XL_getType, getTypeStr, " ARM_ERR: Accessor Type: string expected", C_result);
        char* type = getTypeStr.c_str(); // à cause du new !!
		string getType(type);
        delete type;

		CCString craGetClass(GCGetTypeToClass(getType,craId).c_str());

		exportFunc2Args<long, string> ourFunc(craId,getType, ARMLOCAL_CRASpread_Get);

		/// call the general function
		fillXL_Result( craGetClass, ourFunc, C_result, XL_result, PersistentInXL );
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_CRAGet_Common" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}
							 
LPXLOPER Local_CRASpreadSet_Common(
	LPXLOPER XL_craId,
	LPXLOPER XL_dataId,
	LPXLOPER XL_setPortfolioType,
	LPXLOPER XL_MktDataKeys,
	LPXLOPER XL_update,
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
		XL_GETOBJID( XL_craId, craId, " ARM_ERR: CRA Calculator: Object expected", C_result);

		long dataId;
		XL_GETOBJID( XL_dataId, dataId,	" ARM_ERR: Object expected", C_result);

		CCString setPortfolioTypeStr;
		XL_readStrCellWD(XL_setPortfolioType, setPortfolioTypeStr, "DEFAULT", " ARM_ERR: Portfolio Type: string expected", C_result);
        char * portType = setPortfolioTypeStr.c_str(); // à cause du new !!
		string setPortfolioType(portType);
        delete portType;

		VECTOR<CCString> mktDataKeys;
		VECTOR<CCString> mktDataKeysDef(0);
		XL_readStrVectorWD(XL_MktDataKeys, mktDataKeys, mktDataKeysDef, " ARM_ERR: Market datas keys: array of string expected", DOUBLE_TYPE, C_result);

		vector<string> mktDataKeysSTL(mktDataKeys.size());
		for (size_t i=0; i<mktDataKeys.size(); ++i)
        {
            mktDataKeys[i].toUpper();
			mktDataKeysSTL[i] = CCSTringToSTLString(mktDataKeys[i]);
        }

		CCString updateStr;
		CCString updateDefaultStr("N");
		XL_readStrCellWD(XL_update, updateStr, updateDefaultStr, " ARM_ERR: update flag: string expected", C_result);
		updateStr.toUpper();
		bool isUpdated = (updateStr=="Y" || updateStr=="YES");
		exportFunc5Args<long, long, string, vector<string>, bool> ourFunc( 
								craId,
								dataId,
								setPortfolioType,
								mktDataKeysSTL,
								isUpdated,
								ARMLOCAL_CRASpread_Set);

        if (isUpdated)
        {
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
        else
        {
            /// General call through the functor with an object creation		    
		    fillXL_Result( LOCAL_GC_CCSO_CLASS, ourFunc, C_result, XL_result, PersistentInXL );
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
__declspec(dllexport) LPXLOPER WINAPI Local_CRASpreadCalculator_GetData(
	LPXLOPER XL_craId,
	LPXLOPER XL_getType)
{
	ADD_LOG("Local_CRASpreadCalculator_GetData");
	bool PersistentInXL = true;
	return Local_CRASpreadGet_Common(
	    XL_craId,
	    XL_getType,
        PersistentInXL );
}

__declspec(dllexport) LPXLOPER WINAPI Local_CRASpreadCalculator_SetData(
	LPXLOPER XL_craId,
	LPXLOPER XL_dataId,
    LPXLOPER XL_setPortfolioType,
	LPXLOPER XL_mktDataKeys,
	LPXLOPER XL_update)
{
	ADD_LOG("Local_CRASpreadCalculator_SetData");
	bool PersistentInXL = true;
	return Local_CRASpreadSet_Common(
	    XL_craId,
	    XL_dataId,
        XL_setPortfolioType,
		XL_mktDataKeys,
        XL_update,
        PersistentInXL );
}

///////////////////////////////////
/// PXL version 
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_PXL_CRASpreadCalculator_GetData(
	LPXLOPER XL_craId,
	LPXLOPER XL_getType)
{
	ADD_LOG("Local_PXL_CRASpreadCalculator_GetData");
	bool PersistentInXL = false;
	return Local_CRASpreadGet_Common(
	    XL_craId,
	    XL_getType,
        PersistentInXL );
}

__declspec(dllexport) LPXLOPER WINAPI Local_PXL_CRASpreadCalculator_SetData(
	LPXLOPER XL_craId,
	LPXLOPER XL_dataId,
    LPXLOPER XL_setPortfolioType,
	LPXLOPER XL_mktDataKeys,
	LPXLOPER XL_update)
{
	ADD_LOG("Local_PXL_CRASpreadCalculator_SetData");
	bool PersistentInXL = false;
	return Local_CRASpreadSet_Common(
	    XL_craId,
	    XL_dataId,
        XL_setPortfolioType,
		XL_mktDataKeys,
        XL_update,
        PersistentInXL );
}

__declspec(dllexport) LPXLOPER WINAPI Local_CRASpreadCalculator_Create(
	LPXLOPER XL_generalDatas,
	LPXLOPER XL_callDatas,
	LPXLOPER XL_fundDatas,
	LPXLOPER XL_cpnDatas,
	LPXLOPER XL_notionalCurve,
	LPXLOPER XL_callFeesCurve,
	LPXLOPER XL_fundSpreadCurve,
	LPXLOPER XL_boostedFixCurve,
	LPXLOPER XL_payIndexMultCurve,
	LPXLOPER XL_barrierDownCurve,
	LPXLOPER XL_barrierUpCurve,
	LPXLOPER XL_refCoeff1Curve,
	LPXLOPER XL_refCoeff2Curve,
	LPXLOPER XL_ModelDatas,
	LPXLOPER XL_mktDataManager,
	LPXLOPER XL_productsToPrice,
	LPXLOPER XL_localCalibFlags,
	LPXLOPER XL_miscDatas)
{
	ADD_LOG("Local_CRASpreadCalculator_Create");
	bool PersistentInXL = true;
	bool IsDoubleCorridor = false;
	bool IsTripleRange = false;
	bool IsVms =false;

	return Local_CRASpreadCalculator_Common(
			XL_generalDatas,
			XL_callDatas,
			XL_fundDatas,
			XL_cpnDatas,
			XL_notionalCurve,
			XL_callFeesCurve,
			XL_fundSpreadCurve,
			XL_boostedFixCurve,
			XL_payIndexMultCurve,
			XL_barrierDownCurve,
			XL_barrierUpCurve,
			XL_refCoeff1Curve,
			XL_refCoeff2Curve,
			IsDoubleCorridor,
			XL_barrierDownCurve,
			XL_barrierUpCurve,
			IsTripleRange,
			IsVms,
			XL_ModelDatas,
			XL_mktDataManager,
			XL_productsToPrice,
			XL_localCalibFlags,
			XL_miscDatas,
			PersistentInXL);
}

__declspec(dllexport) LPXLOPER WINAPI Local_PXL_CRASpreadCalculator_Create(
	LPXLOPER XL_generalDatas,
	LPXLOPER XL_callDatas,
	LPXLOPER XL_fundDatas,
	LPXLOPER XL_cpnDatas,
	LPXLOPER XL_notionalCurve,
	LPXLOPER XL_callFeesCurve,
	LPXLOPER XL_fundSpreadCurve,
	LPXLOPER XL_boostedFixCurve,
	LPXLOPER XL_payIndexMultCurve,
	LPXLOPER XL_barrierDownCurve,
	LPXLOPER XL_barrierUpCurve,
	LPXLOPER XL_refCoeff1Curve,
	LPXLOPER XL_refCoeff2Curve,
	LPXLOPER XL_ModelDatas,
	LPXLOPER XL_mktDataManager,
	LPXLOPER XL_productsToPrice,
	LPXLOPER XL_localCalibFlags,
	LPXLOPER XL_miscDatas)
{
	ADD_LOG("Local_PXL_CRASpreadCalculator_Create");
	bool PersistentInXL = false;
	bool IsDoubleCorridor = false;
	bool IsTripleRange =false;
	bool IsVms =false;
			
	return Local_CRASpreadCalculator_Common(
			XL_generalDatas,
			XL_callDatas,
			XL_fundDatas,
			XL_cpnDatas,
			XL_notionalCurve,
			XL_callFeesCurve,
			XL_fundSpreadCurve,
			XL_boostedFixCurve,
			XL_payIndexMultCurve,
			XL_barrierDownCurve,
			XL_barrierUpCurve,
			XL_refCoeff1Curve,
			XL_refCoeff2Curve,
			IsDoubleCorridor,
			XL_barrierDownCurve,
			XL_barrierUpCurve,
			IsTripleRange,
			IsVms,
			XL_ModelDatas,
			XL_mktDataManager,
			XL_productsToPrice,
			XL_localCalibFlags,
			XL_miscDatas,
			PersistentInXL);
}

__declspec(dllexport) LPXLOPER WINAPI Local_CRAVMSCalculator_Create(
	LPXLOPER XL_generalDatas,
	LPXLOPER XL_callDatas,
	LPXLOPER XL_fundDatas,
	LPXLOPER XL_cpnDatas,
	LPXLOPER XL_notionalCurve,
	LPXLOPER XL_callFeesCurve,
	LPXLOPER XL_fundSpreadCurve,
	LPXLOPER XL_boostedFixCurve,
	LPXLOPER XL_payIndexMultCurve,
	LPXLOPER XL_barrierDownCurve,
	LPXLOPER XL_barrierUpCurve,
	LPXLOPER XL_refCoeff1Curve,
	LPXLOPER XL_refTenorCurve,
	LPXLOPER XL_ModelDatas,
	LPXLOPER XL_mktDataManager,
	LPXLOPER XL_productsToPrice,
	LPXLOPER XL_localCalibFlags,
	LPXLOPER XL_miscDatas)
{
	ADD_LOG("Local_CRAVMSCalculator_Create");
	bool PersistentInXL = true;
	bool IsDoubleCorridor = false;
	bool IsTripleRange = false;
	bool IsVms = true;
	

	return Local_CRASpreadCalculator_Common(
			XL_generalDatas,
			XL_callDatas,
			XL_fundDatas,
			XL_cpnDatas,
			XL_notionalCurve,
			XL_callFeesCurve,
			XL_fundSpreadCurve,
			XL_boostedFixCurve,
			XL_payIndexMultCurve,
			XL_barrierDownCurve,
			XL_barrierUpCurve,
			XL_refCoeff1Curve,
			XL_refTenorCurve,
			IsDoubleCorridor,
			XL_barrierDownCurve,
			XL_barrierUpCurve,
			IsTripleRange,
			IsVms,
			XL_ModelDatas,
			XL_mktDataManager,
			XL_productsToPrice,
			XL_localCalibFlags,
			XL_miscDatas,
			PersistentInXL);
}

__declspec(dllexport) LPXLOPER WINAPI Local_PXL_CRAVMSCalculator_Create(
	LPXLOPER XL_generalDatas,
	LPXLOPER XL_callDatas,
	LPXLOPER XL_fundDatas,
	LPXLOPER XL_cpnDatas,
	LPXLOPER XL_notionalCurve,
	LPXLOPER XL_callFeesCurve,
	LPXLOPER XL_fundSpreadCurve,
	LPXLOPER XL_boostedFixCurve,
	LPXLOPER XL_payIndexMultCurve,
	LPXLOPER XL_barrierDownCurve,
	LPXLOPER XL_barrierUpCurve,
	LPXLOPER XL_refCoeff1Curve,
	LPXLOPER XL_refTenorCurve,
	LPXLOPER XL_ModelDatas,
	LPXLOPER XL_mktDataManager,
	LPXLOPER XL_productsToPrice,
	LPXLOPER XL_localCalibFlags,
	LPXLOPER XL_miscDatas)
{
	ADD_LOG("Local_PXL_CRAVMSCalculator_Create");
	bool PersistentInXL = false;
	bool IsDoubleCorridor = false;
	bool IsTripleRange =false;
	bool IsVms =true;
			
	return Local_CRASpreadCalculator_Common(
			XL_generalDatas,
			XL_callDatas,
			XL_fundDatas,
			XL_cpnDatas,
			XL_notionalCurve,
			XL_callFeesCurve,
			XL_fundSpreadCurve,
			XL_boostedFixCurve,
			XL_payIndexMultCurve,
			XL_barrierDownCurve,
			XL_barrierUpCurve,
			XL_refCoeff1Curve,
			XL_refTenorCurve,
			IsDoubleCorridor,
			XL_barrierDownCurve,
			XL_barrierUpCurve,
			IsTripleRange,
			IsVms,
			XL_ModelDatas,
			XL_mktDataManager,
			XL_productsToPrice,
			XL_localCalibFlags,
			XL_miscDatas,
			PersistentInXL);
}

__declspec(dllexport) LPXLOPER WINAPI Local_CMRASpreadCalculator_Create(
	LPXLOPER XL_generalDatas,
	LPXLOPER XL_callDatas,
	LPXLOPER XL_fundDatas,
	LPXLOPER XL_cpnDatas,
	LPXLOPER XL_notionalCurve,
	LPXLOPER XL_callFeesCurve,
	LPXLOPER XL_fundSpreadCurve,
	LPXLOPER XL_boostedFixCurve,
	LPXLOPER XL_payIndexMultCurve,
	LPXLOPER XL_barrierDownCurve,
	LPXLOPER XL_barrierUpCurve,
	LPXLOPER XL_refCoeff1Curve,
	LPXLOPER XL_refCoeff2Curve,
	LPXLOPER XL_ModelDatas,
	LPXLOPER XL_mktDataManager,
	LPXLOPER XL_productsToPrice,
	LPXLOPER XL_localCalibFlags,
	LPXLOPER XL_miscDatas)

{
	ADD_LOG("Local_CMRASpreadCalculator_Create");
	bool PersistentInXL = true;
	bool IsDoubleCorridor = false;
	bool IsTripleRange = true;
	bool IsVms =false;

	return Local_CRASpreadCalculator_Common(
			XL_generalDatas,
			XL_callDatas,
			XL_fundDatas,
			XL_cpnDatas,
			XL_notionalCurve,
			XL_callFeesCurve,
			XL_fundSpreadCurve,
			XL_boostedFixCurve,
			XL_payIndexMultCurve,
			XL_barrierDownCurve,
			XL_barrierUpCurve,
			XL_refCoeff1Curve,
			XL_refCoeff2Curve,
			IsDoubleCorridor,
			XL_barrierDownCurve,
			XL_barrierUpCurve,
			IsTripleRange,
			IsVms,
			XL_ModelDatas,
			XL_mktDataManager,
			XL_productsToPrice,
			XL_localCalibFlags,
			XL_miscDatas,
			PersistentInXL);
}

__declspec(dllexport) LPXLOPER WINAPI Local_PXL_CMRASpreadCalculator_Create(
	LPXLOPER XL_generalDatas,
	LPXLOPER XL_callDatas,
	LPXLOPER XL_fundDatas,
	LPXLOPER XL_cpnDatas,
	LPXLOPER XL_notionalCurve,
	LPXLOPER XL_callFeesCurve,
	LPXLOPER XL_fundSpreadCurve,
	LPXLOPER XL_boostedFixCurve,
	LPXLOPER XL_payIndexMultCurve,
	LPXLOPER XL_barrierDownCurve,
	LPXLOPER XL_barrierUpCurve,
	LPXLOPER XL_refCoeff1Curve,
	LPXLOPER XL_refCoeff2Curve,
	LPXLOPER XL_ModelDatas,
	LPXLOPER XL_mktDataManager,
	LPXLOPER XL_productsToPrice,
	LPXLOPER XL_localCalibFlags,
	LPXLOPER XL_miscDatas)
{
	ADD_LOG("Local_PXL_CMRASpreadCalculator_Create");
	bool PersistentInXL = false;
	bool IsDoubleCorridor = false;
	bool IsTripleRange =true;
	bool IsVms =false;
			
	return Local_CRASpreadCalculator_Common(
			XL_generalDatas,
			XL_callDatas,
			XL_fundDatas,
			XL_cpnDatas,
			XL_notionalCurve,
			XL_callFeesCurve,
			XL_fundSpreadCurve,
			XL_boostedFixCurve,
			XL_payIndexMultCurve,
			XL_barrierDownCurve,
			XL_barrierUpCurve,
			XL_refCoeff1Curve,
			XL_refCoeff2Curve,
			IsDoubleCorridor,
			XL_barrierDownCurve,
			XL_barrierUpCurve,
			IsTripleRange,
			IsVms,
			XL_ModelDatas,
			XL_mktDataManager,
			XL_productsToPrice,
			XL_localCalibFlags,
			XL_miscDatas,
			PersistentInXL);
}


__declspec(dllexport) LPXLOPER WINAPI Local_CRASpreadCalculator_CreateFromPf(
			LPXLOPER XL_optionPortfolio,
			LPXLOPER XL_ModelDatas,
			LPXLOPER XL_payIndexMultId,
			LPXLOPER XL_mktDataManager,
			LPXLOPER XL_productsToPrice,
			LPXLOPER XL_refResetFreq)
{
	ADD_LOG("Local_CRASpreadCalculator_CreateFromPf");
	bool PersistentInXL = true;
	return Local_CRASpreadCalculator_Common(XL_optionPortfolio,
											XL_ModelDatas,
											XL_payIndexMultId,
											XL_mktDataManager,
											XL_productsToPrice,
											XL_refResetFreq,
											PersistentInXL);
}

__declspec(dllexport) LPXLOPER WINAPI Local_PXL_CRASpreadCalculator_CreateFromPf(
			LPXLOPER XL_optionPortfolio,
			LPXLOPER XL_ModelDatas,
			LPXLOPER XL_payIndexMultId,
			LPXLOPER XL_mktDataManager,
			LPXLOPER XL_productsToPrice,
			LPXLOPER XL_refResetFreq)
{
	ADD_LOG("Local_PXL_CRASpreadCalculator_CreateFromPf");
	bool PersistentInXL = false;
	return Local_CRASpreadCalculator_Common(XL_optionPortfolio,
											XL_ModelDatas,
											XL_payIndexMultId,
											XL_mktDataManager,
											XL_productsToPrice,
											XL_refResetFreq,
											PersistentInXL);
}

__declspec(dllexport) LPXLOPER WINAPI Local_CRADoubleCalculator_Create(
	LPXLOPER XL_generalDatas,
	LPXLOPER XL_callDatas,
	LPXLOPER XL_fundDatas,
	LPXLOPER XL_cpnDatas,
	LPXLOPER XL_notionalCurve,
	LPXLOPER XL_callFeesCurve,
	LPXLOPER XL_fundSpreadCurve,
	LPXLOPER XL_boostedFixCurve,
	LPXLOPER XL_payIndexMultCurve,
	LPXLOPER XL_rateBarrierDownCurve,
	LPXLOPER XL_rateBarrierUpCurve,
	LPXLOPER XL_spreadBarrierDownCurve,
	LPXLOPER XL_spreadBarrierUpCurve,
	LPXLOPER XL_spreadCoeff1,
	LPXLOPER XL_spreadCoeff2,
	LPXLOPER XL_ModelDatas,
	LPXLOPER XL_mktDataManager,
	LPXLOPER XL_productsToPrice,
	LPXLOPER XL_localCalibFlags,
	LPXLOPER XL_miscDatas)
{
	ADD_LOG("Local_CRADoubleCalculator_Create");
	bool PersistentInXL = true;
	bool IsDoubleCorridor = true;
	bool IsTripleRange = false;
	bool IsVms = false;

	return Local_CRASpreadCalculator_Common(
			XL_generalDatas,
			XL_callDatas,
			XL_fundDatas,
			XL_cpnDatas,
			XL_notionalCurve,
			XL_callFeesCurve,
			XL_fundSpreadCurve,
			XL_boostedFixCurve,
			XL_payIndexMultCurve,
			XL_spreadBarrierDownCurve,
			XL_spreadBarrierUpCurve,
			XL_spreadCoeff1,
			XL_spreadCoeff2,
			IsDoubleCorridor,
			XL_rateBarrierDownCurve,
			XL_rateBarrierUpCurve,
			IsTripleRange,
			IsVms,
			XL_ModelDatas,
			XL_mktDataManager,
			XL_productsToPrice,
			XL_localCalibFlags,
			XL_miscDatas,
			PersistentInXL);
}

__declspec(dllexport) LPXLOPER WINAPI Local_PXL_CRADoubleCalculator_Create(
	LPXLOPER XL_generalDatas,
	LPXLOPER XL_callDatas,
	LPXLOPER XL_fundDatas,
	LPXLOPER XL_cpnDatas,
	LPXLOPER XL_notionalCurve,
	LPXLOPER XL_callFeesCurve,
	LPXLOPER XL_fundSpreadCurve,
	LPXLOPER XL_boostedFixCurve,
	LPXLOPER XL_payIndexMultCurve,
	LPXLOPER XL_rateBarrierDownCurve,
	LPXLOPER XL_rateBarrierUpCurve,
	LPXLOPER XL_spreadBarrierDownCurve,
	LPXLOPER XL_spreadBarrierUpCurve,
	LPXLOPER XL_spreadCoeff1,
	LPXLOPER XL_spreadCoeff2,
	LPXLOPER XL_ModelDatas,
	LPXLOPER XL_mktDataManager,
	LPXLOPER XL_productsToPrice,
	LPXLOPER XL_localCalibFlags,
	LPXLOPER XL_miscDatas)

{
	ADD_LOG("Local_PXL_CRADoubleCalculator_Create");
	bool PersistentInXL = false;
	bool IsDoubleCorridor = true;
	bool IsTripleRange = false;
	bool IsVms = false;

	return Local_CRASpreadCalculator_Common(
			XL_generalDatas,
			XL_callDatas,
			XL_fundDatas,
			XL_cpnDatas,
			XL_notionalCurve,
			XL_callFeesCurve,
			XL_fundSpreadCurve,
			XL_boostedFixCurve,
			XL_payIndexMultCurve,
			XL_spreadBarrierDownCurve,
			XL_spreadBarrierUpCurve,
			XL_spreadCoeff1,
			XL_spreadCoeff2,
			IsDoubleCorridor,
			XL_rateBarrierDownCurve,
			XL_rateBarrierUpCurve,
			IsTripleRange,
			IsVms,
			XL_ModelDatas,
			XL_mktDataManager,
			XL_productsToPrice,
			XL_localCalibFlags,
			XL_miscDatas,
			PersistentInXL);
}


__declspec(dllexport) LPXLOPER WINAPI Local_BasicCRASpreadCalculator_CreateFromPf(
			LPXLOPER XL_asOfDate,
			LPXLOPER XL_optionPortfolio,
			LPXLOPER XL_payIndexMultId,
			LPXLOPER XL_refResetFreq)
{
	ADD_LOG("Local_BasicCRASpreadCalculator_CreateFromPf");
	bool PersistentInXL = true;
	return Local_CRASpreadCalculator_Common(XL_asOfDate,
											XL_optionPortfolio,
											XL_payIndexMultId,
											XL_refResetFreq,
											PersistentInXL);
}

__declspec(dllexport) LPXLOPER WINAPI Local_PXL_BasicCRASpreadCalculator_CreateFromPf(
			LPXLOPER XL_asOfDate,
			LPXLOPER XL_optionPortfolio,
			LPXLOPER XL_payIndexMultId,
			LPXLOPER XL_refResetFreq)
{
	ADD_LOG("Local_PXL_BasicCRASpreadCalculator_CreateFromPf");
	bool PersistentInXL = false;
	return Local_CRASpreadCalculator_Common(XL_asOfDate,
											XL_optionPortfolio,
											XL_payIndexMultId,
											XL_refResetFreq,
											PersistentInXL);
}

__declspec(dllexport) LPXLOPER WINAPI Local_BasicCRASpreadCalculator_CreateFromSwaption(
			LPXLOPER XL_asOfDate,
			LPXLOPER XL_swaption,
			LPXLOPER XL_payIndexMultId,
			LPXLOPER XL_refResetFreq)
{
	ADD_LOG("Local_BasicCRASpreadCalculator_CreateFromSwaption");
	bool PersistentInXL = true;
	return Local_CRASpreadCalculator_Common(XL_asOfDate,
											XL_swaption,
											XL_payIndexMultId,
											XL_refResetFreq,
											PersistentInXL);
}

__declspec(dllexport) LPXLOPER WINAPI Local_PXL_BasicCRASpreadCalculator_CreateFromSwaption(
			LPXLOPER XL_asOfDate,
			LPXLOPER XL_swaption,
			LPXLOPER XL_payIndexMultId,
			LPXLOPER XL_refResetFreq)
{
	ADD_LOG("Local_PXL_BasicCRASpreadCalculator_CreateFromSwaption");
	bool PersistentInXL = false;
	return Local_CRASpreadCalculator_Common(XL_asOfDate,
											XL_swaption,
											XL_payIndexMultId,
											XL_refResetFreq,
											PersistentInXL);
}
/****************************************************************
					Snow Range Calculator
*****************************************************************/
class snowRangeCalculatorFunc : public ARMResultLong2LongFunc
{
public:
	snowRangeCalculatorFunc(ARM_Currency ccy,
							double startDate,
							double endDate,
							int payReceive,
							double notional,
							long notionalId,
							string fundingIndexTerm,
							int fundingDayCount,
							string couponIndexTerm,
							int couponDayCount,
							int resetFreq,
							int payFreq,
							int resetTiming,
							int payTiming,
							int resetGap,
							string resetCal,
							string payCal,
							int adjRule,
							int intRule,
							double spread,
							long spreadId,
							double strike,
							long strikeId,
							double ratchet,
							long ratchetId,
							double cashFlow,
							long cashFlowId,
							double fixedRate,
							long fixedRateId,
							double leverage,
							long leverageId,
							ARM_Vector* snowRangeParams,
							ARM_Vector* calibParams,
							string modelName,
							ARM_Vector* modelParams,
							int	nbSteps,
							string generatorType,
							string inversionMethod,
							bool antithetic,
							int samplerType,
							vector<string> mdmKeys,
							long mktDataManagerId,
							vector<string> productsToPrice)
    :
	C_currency(ccy),
    C_startDate(startDate),
    C_endDate(endDate),
	C_payReceive(payReceive),
	C_notional(notional),
	C_notionalId(notionalId),
	C_fundingIndexTerm(fundingIndexTerm),
	C_fundingDayCount(fundingDayCount),
	C_couponIndexTerm(couponIndexTerm),
	C_couponDayCount(couponDayCount),
	C_resetFreq(resetFreq),
	C_payFreq(payFreq),
	C_resetTiming(resetTiming),
	C_payTiming(payTiming),
	C_resetGap(resetGap),
	C_resetCal(resetCal),
	C_payCal(payCal),
	C_adjRule(adjRule),
	C_intRule(intRule),
	C_spread(spread),
	C_spreadId(spreadId),
	C_strike(strike),
	C_strikeId(strikeId),	
	C_ratchet(ratchet),
	C_ratchetId(ratchetId),
	C_cashFlow(cashFlow),
	C_cashFlowId(cashFlowId),
	C_fixedRate(fixedRate),
	C_fixedRateId(fixedRateId),
	C_leverage(leverage),
	C_leverageId(leverageId),
	C_snowRangeParams(snowRangeParams),
	C_calibParams(calibParams),
	C_modelName(modelName),
	C_modelParams(modelParams),
	C_nbSteps(nbSteps),
	C_generatorType(generatorType),
	C_inversionMethod(inversionMethod),
	C_antithetic(antithetic),
	C_samplerType(samplerType),
	C_mdmKeys(mdmKeys),
	C_mktDataManagerId(mktDataManagerId),
	C_productsToPrice(productsToPrice),
	C_portfolioId(-1)
    {
	};
	
/*	snowRangeCalculatorFunc(
						long portfolioId,
						double fundLev,
						long fundLevId,
						double capLev,
						long capLevId,
						ARM_Vector* snowRangeParams,
						int	nbSteps,
						string generatorType,
						string inversionMethod,
						bool antithetic,
						int samplerType,
						ARM_Vector* calibParams,
						vector<string> mdmKeys,
						long mktDataManagerId,
				        vector<string> productsToPrice)
    :
	C_portfolioId(portfolioId),
	C_fundLev(fundLev),
	C_fundLevId(fundLevId),
	C_capLev(capLev),
	C_capLevId(capLevId),	
	C_snowRangeParams(snowRangeParams),
	C_nbSteps(nbSteps),
	C_generatorType(generatorType),
	C_inversionMethod(inversionMethod),
	C_antithetic(antithetic),
	C_samplerType(samplerType),
	C_calibParams(calibParams),
	C_mdmKeys(mdmKeys),
	C_mktDataManagerId(mktDataManagerId),
	C_productsToPrice(productsToPrice)
    {
	};
*/	
	long operator()( ARM_result& result, long objId )
	{
/*		if (C_portfolioId != -1)
		{
			return ARMLOCAL_SnowRangeCalculator_Create(
				C_portfolioId,
				C_fundLev,
				C_fundLevId,
				C_capLev,
				C_capLevId,	
				C_snowRangeParams,
				C_nbSteps,
				C_generatorType,
				C_inversionMethod,
				C_antithetic,
				C_samplerType,
				C_calibParams,
				C_mdmKeys,
				C_mktDataManagerId,
				C_productsToPrice,
				result,
				objId);
		}
		else
		{
*/			return ARMLOCAL_SnowRangeCalculator_Create( 
				C_currency,
				C_startDate,
				C_endDate,
				C_payReceive,
				C_notional,
				C_notionalId,
				C_fundingIndexTerm,
				C_fundingDayCount,
				C_couponIndexTerm,
				C_couponDayCount,
				C_resetFreq,
				C_payFreq,
				C_resetTiming,
				C_payTiming,
				C_resetGap,
				C_resetCal,
				C_payCal,
				C_adjRule,
				C_intRule,
				C_spread,
				C_spreadId,
				C_strike,
				C_strikeId,	
				C_ratchet,
				C_ratchetId,
				C_cashFlow,
				C_cashFlowId,
				C_fixedRate,
				C_fixedRateId,
				C_leverage,
				C_leverageId,
				C_snowRangeParams,
				C_calibParams,
				C_modelName,
				C_modelParams,
				C_nbSteps,
				C_generatorType,
				C_inversionMethod,
				C_antithetic,
				C_samplerType,
				C_mdmKeys,
				C_mktDataManagerId,
				C_productsToPrice,
				result,
				objId);
//		}
    }

private:
int					C_portfolioId;
ARM_Currency		C_currency;
double				C_startDate;
double				C_endDate;
int					C_payReceive;
string				C_fundingIndexTerm;
int					C_fundingDayCount;
string				C_couponIndexTerm;
int					C_couponDayCount;
int					C_resetFreq;
int					C_payFreq;
int					C_resetTiming;
int					C_payTiming;
int					C_resetGap;
string				C_resetCal;
string				C_payCal;
int					C_adjRule;
int					C_intRule;
double				C_spread;
long				C_spreadId;
double				C_notional;
long				C_notionalId;
double				C_strike;
long				C_strikeId;
double				C_ratchet;
long				C_ratchetId;
double				C_cashFlow;
long				C_cashFlowId;
double				C_fixedRate;
long				C_fixedRateId;
double				C_leverage;
long				C_leverageId;

ARM_Vector*			C_snowRangeParams;
ARM_Vector*			C_calibParams;
string				C_modelName;
ARM_Vector*			C_modelParams;

int					C_nbSteps;
string				C_generatorType;
string				C_inversionMethod;
bool				C_antithetic;
int					C_samplerType;

vector<string>		C_mdmKeys;
long				C_mktDataManagerId;
vector<string>		C_productsToPrice;
};

LPXLOPER Local_SnowRangeCalculator_Common ( LPXLOPER XL_generalData,
											LPXLOPER XL_notionalCurve,
											LPXLOPER XL_fundingData,
											LPXLOPER XL_couponData,
											LPXLOPER XL_spreadCurve,
											LPXLOPER XL_strikeCurve,			
											LPXLOPER XL_ratchetCurve,
											LPXLOPER XL_cashFlowCurve,
											LPXLOPER XL_fixedRateCurve,
											LPXLOPER XL_leverageCurve,
											LPXLOPER XL_snowRangeParams,
											LPXLOPER XL_calibParams,
											LPXLOPER XL_modelParams,
											LPXLOPER XL_mcParams,
											LPXLOPER XL_mktDataManager,
											LPXLOPER XL_productsToPrice,
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
	
		//General Data
		//-------------
		VECTOR<CCString> generalData;
		XL_readStrVector (XL_generalData, generalData, "ARM_ERR: General data: Ccy, Start, End, P/R", DOUBLE_TYPE, C_result);
		if (generalData.size() != 4)
		{
			throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"General data: 4 parameters expected");
		}

		//Currency
		CCString ccyStr = generalData[0];
		string sCcy = CCSTringToSTLString(ccyStr);
		ARM_Currency ccy(sCcy.c_str());

		//Start & End
		CCString startDateXl = generalData[1];
		string sStartDate = CCSTringToSTLString(startDateXl);
		double startDate = atoi(sStartDate.c_str());
		CCString endDateXl = generalData[2];
		string sEndDate = CCSTringToSTLString(endDateXl);
		double endDate = atoi(sEndDate.c_str());
	    
		//PayReceive
		CCString payReceiveStr = generalData[3];
		long lPayReceive;
		if ((lPayReceive = ARM_ConvRecOrPay (payReceiveStr, C_result)) == ARM_DEFAULT_ERR)
		{
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}
		int payReceive = (int)lPayReceive;

		//Notional Value or Notional Curve
		double notional;
		CCString notionalStr;
		long notionalId;
		XL_readStrOrNumCell(XL_notionalCurve, notionalStr, notional, notionalId,
			   " ARM_ERR: notional: numeric or refValue Id expected",C_result);	

		if (notionalId == XL_TYPE_STRING)
			notionalId = LocalGetNumObjectId(notionalStr);
		else
			notionalId = ARM_NULL_OBJECT;

		// FUNDING DATA
		//-------------
		VECTOR<CCString> fundingData;
		XL_readStrVector (XL_fundingData, fundingData," ARM_ERR: Funding data: array of string expected",DOUBLE_TYPE,C_result);
		if (fundingData.size() != 2)
		{
			throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"Funding data: 2 parameters expected");
		}

		//Index Term
		CCString fundingIndexTermStr = fundingData[0];
		string fundingIndexTerm = CCSTringToSTLString(fundingIndexTermStr);

		//Day Count
		CCString fundingDayCountStr = fundingData[1];
		int fundingDayCount = (int)ARM_ConvDayCount(fundingDayCountStr);

		// COUPON DATA
		//-------------
		VECTOR<CCString> couponData;
		XL_readStrVector (XL_couponData, couponData," ARM_ERR: Coupon data: array of string expected",DOUBLE_TYPE,C_result);
		int couponSize = couponData.size();
		if ((couponSize < 6) || (couponSize > 11))
		{
			throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"Coupon data: 6 to 11 parameters expected");
		}

		//Index Term
		CCString couponIndexTermStr = couponData[0];
		string couponIndexTerm = CCSTringToSTLString(couponIndexTermStr);

		//Day Count
		CCString couponDayCountStr = couponData[1];
		int couponDayCount = (int)ARM_ConvDayCount(couponDayCountStr);

		//Reset Freq
		CCString resetFreqStr = couponData[2];
		long lcouponResetFreq;
		if ((lcouponResetFreq = ARM_ConvFrequency (resetFreqStr, C_result)) == ARM_DEFAULT_ERR)
		{
			throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"Coupon data: invalid reset frequency");
		}
		int couponResetFreq = (int)lcouponResetFreq;

		//Pay Freq
		CCString payFreqStr = couponData[3];
		long lcouponPayFreq;
		if ((lcouponPayFreq = ARM_ConvFrequency (payFreqStr, C_result)) == ARM_DEFAULT_ERR)
		{
			throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"Coupon data: invalid pay frequency");
		}
		int couponPayFreq = (int)lcouponPayFreq;

		//coupon Reset Timing
		CCString couponResetTimingStr = couponData[4];
		int couponResetTiming = (int)ARM_ConvPayResetRule(couponResetTimingStr);

		//coupon Pay Timing
		CCString couponPayTimingStr = couponData[5];
		int couponPayTiming = (int)ARM_ConvPayResetRule(couponPayTimingStr);

		//coupon Reset Gap
		int couponResetGap = -ccy.GetSpotDays();
		if (couponSize > 6)
		{
			CCString couponResetGapXL = couponData[6];
			string couponResetGapStr = CCSTringToSTLString(couponResetGapXL);
			couponResetGap = atoi(couponResetGapStr.c_str());
		}

		//coupon Reset Calendar
		char* resetCal = NULL;
		resetCal = ccy.GetResetCalName(ccy.GetVanillaIndexType());
		string couponResetCal(resetCal);
		if (couponSize > 7)
		{
			CCString couponResetCalStr = couponData[7];
			couponResetCal = CCSTringToSTLString(couponResetCalStr);
		}

		//coupon Pay Calendar
		char* payCal = NULL;
		payCal = ccy.GetPayCalName(ccy.GetVanillaIndexType());
		string couponPayCal(payCal);
		if (couponSize > 8)
		{
			CCString couponPayCalStr = couponData[8];
			couponPayCal = CCSTringToSTLString(couponPayCalStr);
		}

		//coupon Adjustment Rule
		int couponAdjRule = K_MOD_FOLLOWING;
		if (couponSize > 9)
		{
			CCString couponAdjRuleStr = couponData[9];
			couponAdjRule = (int)ARM_ConvFwdRule(couponAdjRuleStr);
		}

		//coupon Int Rule
		int couponIntRule = K_ADJUSTED;
		if (couponSize > 10)
		{
			CCString couponIntRuleStr = couponData[10];
			couponIntRule = (int)ARM_ConvIntRule(couponIntRuleStr);
		}

		//spread Value or spread Curve
		double spread;
		CCString spreadStr;
		long spreadId;
		XL_readStrOrNumCell(XL_spreadCurve, spreadStr, spread, spreadId,
			   " ARM_ERR: spread: numeric or refValue Id expected",C_result);	

		if (spreadId == XL_TYPE_STRING)
			spreadId = LocalGetNumObjectId(spreadStr);
		else
			spreadId = ARM_NULL_OBJECT;

		//Strike Value or Strike Curve
		double strike;
		CCString strikeStr;
		long strikeId;
		XL_readStrOrNumCell(XL_strikeCurve, strikeStr, strike, strikeId,
			   " ARM_ERR: strike: numeric or refValue Id expected",C_result);	

		if (strikeId == XL_TYPE_STRING)
			strikeId = LocalGetNumObjectId(strikeStr);
		else
			strikeId = ARM_NULL_OBJECT;

		//Ratchet Value or Ratchet Curve
		double ratchet;
		CCString ratchetStr;
		long ratchetId;
		XL_readStrOrNumCell(XL_ratchetCurve, ratchetStr, ratchet, ratchetId,
			   " ARM_ERR: ratchet: numeric or refValue Id expected",C_result);	

		if (ratchetId == XL_TYPE_STRING)
			ratchetId = LocalGetNumObjectId(ratchetStr);
		else
			ratchetId = ARM_NULL_OBJECT;
		
		//Cash Flow Value or Cash Flow Curve
		double cashFlow;
		CCString cashFlowStr;
		long cashFlowId;
		XL_readStrOrNumCell(XL_cashFlowCurve, cashFlowStr, cashFlow, cashFlowId,
			   " ARM_ERR: cashFlow: numeric or refValue Id expected",C_result);	

		if (cashFlowId == XL_TYPE_STRING)
			cashFlowId = LocalGetNumObjectId(cashFlowStr);
		else
			cashFlowId = ARM_NULL_OBJECT;

		//Fixed Rate Value or Fixed Rate Curve
		double fixedRate;
		CCString fixedRateStr;
		long fixedRateId;
		XL_readStrOrNumCell(XL_fixedRateCurve, fixedRateStr, fixedRate, fixedRateId,
			   " ARM_ERR: fixedRate: numeric or refValue Id expected",C_result);	

		if (fixedRateId == XL_TYPE_STRING)
			fixedRateId = LocalGetNumObjectId(fixedRateStr);
		else
			fixedRateId = ARM_NULL_OBJECT;

		//Leverage Value or Leverage Curve
		double leverage;
		CCString leverageStr;
		long leverageId;
		XL_readStrOrNumCell(XL_leverageCurve, leverageStr, leverage, leverageId,
			   " ARM_ERR: leverage: numeric or refValue Id expected",C_result);	

		if (leverageId == XL_TYPE_STRING)
			leverageId = LocalGetNumObjectId(leverageStr);
		else
			leverageId = ARM_NULL_OBJECT;

		// SNOW RANGE PARAMS
		//------------------
		VECTOR<CCString> snowRangeParamsXL;
		XL_readStrVector (XL_snowRangeParams, snowRangeParamsXL, "ARM_ERR: Snow Range Params: array of string", DOUBLE_TYPE, C_result);
		if (snowRangeParamsXL.size() != 4)
		{
			throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"Snow range params: 4 parameters expected");
		}

		ARM_Vector* snowRangeParams = new ARM_Vector();

		//Barrier Shift
		CCString barrierShiftStr = snowRangeParamsXL[0];
		string sBarrierShiftStr = CCSTringToSTLString(barrierShiftStr);
		snowRangeParams->push_back(atof(sBarrierShiftStr.c_str()));

		//callable
		CCString callableStr = snowRangeParamsXL[1];
		string callable = CCSTringToSTLString(callableStr);
		snowRangeParams->push_back((callable == "YES" || callable == "Y") ? 1.0 : 0.0);

		//Up Front
		CCString upFrontStr = snowRangeParamsXL[2];
		string sUpFront = CCSTringToSTLString(upFrontStr);
		snowRangeParams->push_back(atof(sUpFront.c_str()));

		//Theta
		CCString thetaStr = snowRangeParamsXL[3];
		string sThetaStr = CCSTringToSTLString(thetaStr);
		snowRangeParams->push_back(atof(sThetaStr.c_str()));

		// CALIBRATION PARAMS
		//-------------------
		VECTOR<CCString> calibParamsVect;
		XL_readStrVector (XL_calibParams, calibParamsVect, "ARM_ERR: Calibration Params: array of string expected", DOUBLE_TYPE, C_result);
		if (calibParamsVect.size() != 4)
		{
			throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"Calib params: 4 parameters expected");
		}

		ARM_Vector* calibParams = new ARM_Vector();

		//Time Steps Nb
		CCString timeStepNb = calibParamsVect[0];
		string sTimeStepNbStr = CCSTringToSTLString(timeStepNb);
		calibParams->push_back(atof(sTimeStepNbStr.c_str()));

		//Space Steps Nb
		CCString spaceStepNb = calibParamsVect[1];
		string sSpaceStepNbStr = CCSTringToSTLString(spaceStepNb);
		calibParams->push_back(atof(sSpaceStepNbStr.c_str()));

		//Std Dev Nb
		CCString stdDevNb = calibParamsVect[2];
		string sStdDevNbStr = CCSTringToSTLString(stdDevNb);
		calibParams->push_back(atof(sStdDevNbStr.c_str()));

		//Calibrate Swaption
		CCString calibSwoptStr = calibParamsVect[3];
		string sCalibSwoptStr = CCSTringToSTLString(calibSwoptStr);
		calibParams->push_back(stringGetUpper(sCalibSwoptStr) == "YES");

		//MODEL Params
		//-------------
		VECTOR<CCString> modelParamsVect;
		XL_readStrVector (XL_modelParams, modelParamsVect, "ARM_ERR: Model Params: array of string expected", DOUBLE_TYPE, C_result);
		string modelName = CCSTringToSTLString(modelParamsVect[0]);

		ARM_Vector* modelParams = new ARM_Vector();
		int size = modelParamsVect.size();
		for (int i=1; i<size; i++)
		{
			string param = CCSTringToSTLString(modelParamsVect[i]);
			modelParams->push_back(atof(param.c_str()));
		}

		// MC PARAMS
		//-------------
		VECTOR<CCString> mcParams;
		XL_readStrVector (XL_mcParams, mcParams, "ARM_ERR: Monte Carlo Params: array of string expected", DOUBLE_TYPE, C_result);
		if (mcParams.size() != 5)
		{
			throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"MC : 5 parameters expected");
		}

		//NbSteps
		CCString nbStepsStr = mcParams[0];
		string sNbStepsStr = CCSTringToSTLString(nbStepsStr);
		int nbSteps = atoi(sNbStepsStr.c_str());

		//Generator Type
		CCString generatorTypeStr = mcParams[1];
		string generatorType = CCSTringToSTLString(generatorTypeStr);

		//Inversion Method
		CCString inversionMethodStr = mcParams[2];
		string inversionMethod = CCSTringToSTLString(inversionMethodStr);

		//Antithetic
		CCString antitheticStr = mcParams[3];
		string sAntithetic = CCSTringToSTLString(antitheticStr);
		bool antithetic = (stringGetUpper(sAntithetic) == "YES");

		//Sampler Type
		CCString samplerTypeStr = mcParams[4];
		int samplerType;
		if( (samplerType = ARM_ConvGPSamplerType( samplerTypeStr, C_result)) == ARM_DEFAULT_ERR )
		{
			throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"invalid sampler type");
		}

		// MARKET DATA MANAGER
		//--------------------
		VECTOR<CCString> mktDataManager;
		XL_readStrVector(XL_mktDataManager,mktDataManager," ARM_ERR: Market data: array of string expected",DOUBLE_TYPE,C_result);
        
		if (mktDataManager.size()<5)
		{
			throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"Market data: at least MDM + 4 model parameters");
		}

		long mktDataManagerId = LocalGetNumObjectId (mktDataManager[0]);
		vector<string> mdmKeys(mktDataManager.size()-1);		
		for (i=1; i<mktDataManager.size(); ++i)
			mdmKeys[i-1] = CCSTringToSTLString(mktDataManager[i]);
	
		// PRICING FLAGS
		//--------------
		VECTOR<CCString> pricingFlags;

		XL_readStrVector(XL_productsToPrice,pricingFlags," ARM_ERR: Pricing flags: array of string expected",DOUBLE_TYPE,C_result);
		if (pricingFlags.size()>14)
		{
			throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"maximum: 14 products to price");
		}

		vector< string > productsToPrice(pricingFlags.size(), string("N"));
		for (i=0; i<pricingFlags.size(); i++)
			productsToPrice[i] = CCSTringToSTLString(pricingFlags[i]);

		//Functor call
		snowRangeCalculatorFunc ourFunc(ccy,
										startDate,
										endDate,
										payReceive,
										notional,
										notionalId,
										fundingIndexTerm,
										fundingDayCount,
										couponIndexTerm,
										couponDayCount,
										couponResetFreq,
										couponPayFreq,
										couponResetTiming,
										couponPayTiming,
										couponResetGap,
										couponResetCal,
										couponPayCal,
										couponAdjRule,
										couponIntRule,
										spread,
										spreadId,
										strike,
										strikeId,
										ratchet,
										ratchetId,
										cashFlow,
										cashFlowId,
										fixedRate,
										fixedRateId,
										leverage,
										leverageId,
										snowRangeParams,
										calibParams,
										modelName,
										modelParams,
										nbSteps,
										generatorType,
										inversionMethod,
										antithetic,
										samplerType,
										mdmKeys,
										mktDataManagerId,
										productsToPrice);

		/// call the general function
		fillXL_Result( LOCAL_GC_CALLABLE_SRG_CLASS, ourFunc, C_result, XL_result, PersistentInXL );	
	
		delete resetCal;
		delete payCal;
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_SnowRangeCalculator_Common" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}

__declspec(dllexport) LPXLOPER WINAPI Local_SnowRangeCalculator_Create(
			LPXLOPER XL_generalData,
			LPXLOPER XL_notionalCurve,
			LPXLOPER XL_fundingData,
			LPXLOPER XL_couponData,
			LPXLOPER XL_spreadCurve,
			LPXLOPER XL_strikeCurve,			
			LPXLOPER XL_ratchetCurve,			
			LPXLOPER XL_cashFlowCurve,
			LPXLOPER XL_fixedRateCurve,
			LPXLOPER XL_leverageCurve,
			LPXLOPER XL_snowRangeParams,
			LPXLOPER XL_calibParams,
			LPXLOPER XL_modelParams,
			LPXLOPER XL_mcParams,
			LPXLOPER XL_mktDataManager,
			LPXLOPER XL_productsToPrice)
{
	ADD_LOG("Local_SnowRangeCalculator_Create");
	bool PersistentInXL = true;

	return Local_SnowRangeCalculator_Common(
			XL_generalData,
			XL_notionalCurve,
			XL_fundingData,
			XL_couponData,
			XL_spreadCurve,
			XL_strikeCurve,	
			XL_ratchetCurve,
			XL_cashFlowCurve,
			XL_fixedRateCurve,
			XL_leverageCurve,
			XL_snowRangeParams,
			XL_calibParams,
			XL_modelParams,
			XL_mcParams,
			XL_mktDataManager,
			XL_productsToPrice,
			PersistentInXL);
}

__declspec(dllexport) LPXLOPER WINAPI Local_PXL_SnowRangeCalculator_Create(
			LPXLOPER XL_generalData,
			LPXLOPER XL_notionalCurve,
			LPXLOPER XL_fundingData,
			LPXLOPER XL_couponData,
			LPXLOPER XL_spreadCurve,
			LPXLOPER XL_strikeCurve,			
			LPXLOPER XL_ratchetCurve,			
			LPXLOPER XL_cashFlowCurve,
			LPXLOPER XL_fixedRateCurve,
			LPXLOPER XL_leverageCurve,
			LPXLOPER XL_snowRangeParams,
			LPXLOPER XL_calibParams,
			LPXLOPER XL_modelParams,
			LPXLOPER XL_mcParams,
			LPXLOPER XL_mktDataManager,
			LPXLOPER XL_productsToPrice)
{
	ADD_LOG("Local_PXL_SnowRangeCalculator_Create");
	bool PersistentInXL = false;

	return Local_SnowRangeCalculator_Common(
			XL_generalData,
			XL_notionalCurve,
			XL_fundingData,
			XL_couponData,
			XL_spreadCurve,
			XL_strikeCurve,	
			XL_ratchetCurve,
			XL_cashFlowCurve,
			XL_fixedRateCurve,
			XL_leverageCurve,
			XL_snowRangeParams,
			XL_calibParams,
			XL_modelParams,
			XL_mcParams,
			XL_mktDataManager,
			XL_productsToPrice,
			PersistentInXL);
}

/****************************************************************
	General Swaption Converter to Variable Notional Swaption
*****************************************************************/
class convertToVarNotionalSwaptionFunc : public ARMResultLong2LongFunc
{
public:
	convertToVarNotionalSwaptionFunc(long yieldCurveId,long swaptionId)
    : C_yieldCurveId(yieldCurveId), C_swaptionId(swaptionId) {};
	
	long operator()( ARM_result& result, long objId )
	{
			return ARMLOCAL_ConvertToVarNotionalSwaption(
				C_yieldCurveId,C_swaptionId,result,objId);
    }

private:
	long C_yieldCurveId;
	long C_swaptionId;
};

LPXLOPER Local_ConvertToVarNotionalSwaption_Common( LPXLOPER XL_yieldCurve,
											LPXLOPER XL_swaption,
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
	
	    CCString yieldCurveStr;
	    XL_readStrCell(XL_yieldCurve,yieldCurveStr," ARM_ERR: Yield Curve id: object expected",C_result);
	    long yieldCurveId = LocalGetNumObjectId(yieldCurveStr);

		CCString swaptionStr;
		XL_readStrCell( XL_swaption, swaptionStr, " ARM_ERR: Swaption Id: Object expected",C_result);
		long swaptionId = LocalGetNumObjectId(swaptionStr);

		//Functor call
		convertToVarNotionalSwaptionFunc ourFunc(yieldCurveId,swaptionId);

		/// call the general function
		fillXL_Result( LOCAL_SWAPTION_CLASS, ourFunc, C_result, XL_result, PersistentInXL );	
	
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_ConvertToVarNotionalSwaption_Common" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


__declspec(dllexport) LPXLOPER WINAPI Local_ConvertToVarNotionalSwaption(
			LPXLOPER XL_yieldCurve,
			LPXLOPER XL_swaption)
{
	ADD_LOG("Local_ConvertToVarNotionalSwaption");
	bool PersistentInXL = true;

	return Local_ConvertToVarNotionalSwaption_Common(
			XL_yieldCurve,
			XL_swaption,
			PersistentInXL);
}
__declspec(dllexport) LPXLOPER WINAPI Local_PXL_ConvertToVarNotionalSwaption(
			LPXLOPER XL_yieldCurve,
			LPXLOPER XL_swaption)
{
	ADD_LOG("Local_PXL_ConvertToVarNotionalSwaption");
	bool PersistentInXL = false;

	return Local_ConvertToVarNotionalSwaption_Common(
			XL_yieldCurve,
			XL_swaption,
			PersistentInXL);
}

/////////////////////////////////////////
/// Local_DateStripGetData function
/////////////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_BasisConverter(
	LPXLOPER XL_AsOfDate,
	LPXLOPER XL_DomCcy,
	LPXLOPER XL_ForCcy,
    LPXLOPER XL_DomDateStripId,
	LPXLOPER XL_ForDateStripId,
	LPXLOPER XL_FundDateStripId,
    LPXLOPER XL_DomDayCount,
	LPXLOPER XL_DomFreq,
	LPXLOPER XL_ForDayCount,
	LPXLOPER XL_ForFreq,
	LPXLOPER XL_DomZcId,
	LPXLOPER XL_ForZcId,
	LPXLOPER XL_DomDiscZcId,
	LPXLOPER XL_ForDiscZcId,
	LPXLOPER XL_ForexId,
	LPXLOPER XL_DomNotionalId,
	LPXLOPER XL_ForNotionalId,
	LPXLOPER XL_ForSpreadId)
{
	ADD_LOG("Local_BasisConverter");
	/// to remove memory leak put a result holder!
	static XLOPER_Holder XL_resultHolder;
	XLOPER& XL_result = XL_resultHolder.GetResult();

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

	// AsOfDate
	double C_AsOfDate;
	XL_readNumCell( XL_AsOfDate, C_AsOfDate, " ARM_ERR: AsOfDate: date expected",	C_result);

	// Domestic currency
	CCString C_DomCcy;
	XL_readStrCell( XL_DomCcy, C_DomCcy, " ARM_ERR: Domestic Currency : string expected", C_result);

	// Foreign currency
	CCString C_ForCcy;
	XL_readStrCell( XL_ForCcy, C_ForCcy, " ARM_ERR: Foreign Currency : string expected", C_result);

	/// Domestic DateStrip
	long C_DomDateStripId;
	XL_GETOBJID( XL_DomDateStripId,	C_DomDateStripId,	" ARM_ERR: domestic date strip: Object expected",	C_result);

	/// Foreign DateStrip
	long C_ForDateStripId;
	XL_GETOBJID( XL_ForDateStripId,	C_ForDateStripId,	" ARM_ERR: foreign date strip: Object expected",	C_result);

	/// Funding DateStrip
	long C_FundDateStripId;
	XL_GETOBJID( XL_FundDateStripId, C_FundDateStripId,	" ARM_ERR: funding date strip: Object expected",	C_result);

	/// Domestic DayCount
	CCString C_DomDayCount;
	XL_readStrCell( XL_DomDayCount, C_DomDayCount, " ARM_ERR: Domestic Day Count : string expected", C_result);

	/// Domestic Freq
	CCString C_DomFreq;
	XL_readStrCell( XL_DomFreq, C_DomFreq, " ARM_ERR: Domestic Frequency : string expected", C_result);

	/// Foreign DayCount
	CCString C_ForDayCount;
	XL_readStrCell( XL_ForDayCount, C_ForDayCount, " ARM_ERR: Foreign Day Count : string expected", C_result);

	/// Foreign Frequency
	CCString C_ForFreq;
	XL_readStrCell( XL_ForFreq, C_ForFreq, " ARM_ERR: Foreign Frequency : string expected", C_result);

	/// Domestic Zero Curve
	long C_DomZcId;
	XL_GETOBJID( XL_DomZcId,	C_DomZcId,	" ARM_ERR: domestic curve: Object expected",	C_result);

	/// Foreign Zero Curve
	long C_ForZcId;
	XL_GETOBJID( XL_ForZcId,	C_ForZcId,	" ARM_ERR: foreign curve: Object expected",	C_result);

	/// Domestic Zero Discount Curve
	long C_DomDiscZcId;
	XL_GETOBJID( XL_DomDiscZcId,	C_DomDiscZcId,	" ARM_ERR: domestic discount curve: Object expected",	C_result);

	/// Foreign Zero Discount Curve
	long C_ForDiscZcId;
	XL_GETOBJID( XL_ForDiscZcId,	C_ForDiscZcId,	" ARM_ERR: foreign discount curve: Object expected",	C_result);

	/// Forex
	long C_ForexId;
	XL_GETOBJID( XL_ForexId,	C_ForexId,	" ARM_ERR: forex: Object expected",	C_result);

	/// Domestic Notional
	long C_DomNotionalId;
	XL_GETOBJID( XL_DomNotionalId,	C_DomNotionalId,	" ARM_ERR: domestic notional: Object expected",	C_result);

	/// Foreign Notional
	long C_ForNotionalId;
	XL_GETOBJID( XL_ForNotionalId,	C_ForNotionalId,	" ARM_ERR: foreign notional: Object expected",	C_result);

	/// Foreign Spread
	long C_ForSpreadId;
	XL_GETOBJID( XL_ForSpreadId,	C_ForSpreadId,	" ARM_ERR: foreign spread: Object expected",	C_result);


	VECTOR<double> C_DataResult;
	long retCode;
	retCode = ARMLOCAL_BasisConverter( 
		C_AsOfDate,
		CCSTringToSTLString(C_ForCcy),
		CCSTringToSTLString(C_DomCcy),
		C_DomDateStripId,
		C_ForDateStripId,
		C_FundDateStripId,
		CCSTringToSTLString(C_DomDayCount),
		CCSTringToSTLString(C_DomFreq),
		CCSTringToSTLString(C_ForDayCount),
		CCSTringToSTLString(C_ForFreq),
		C_DomZcId,
		C_ForZcId,
		C_DomDiscZcId,
		C_ForDiscZcId,
		C_ForexId,
		C_DomNotionalId,
		C_ForNotionalId,
		C_ForSpreadId,
		C_DataResult,
		C_result);
	
		
	/// feed the LPXLOPER object result 
	if (retCode == ARM_OK)
	{
		XL_writeNumVector( XL_result, C_DataResult, " ARM_ERR: Could not get result data", C_result);
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_BasisConverter" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}

__declspec(dllexport) LPXLOPER WINAPI Local_CRAQuantoCalculator_Common(LPXLOPER XL_generalDatas,
																	   LPXLOPER XL_callDatas,
																	   LPXLOPER XL_fundDatas,
																	   LPXLOPER XL_cpnDatas,
																	   LPXLOPER XL_notionalCurve,
																	   LPXLOPER XL_callFeesCurve,
																	   LPXLOPER XL_fundSpreadCurve,
																	   LPXLOPER XL_fixCurve,
																	   LPXLOPER XL_barrierDownCurve,
																	   LPXLOPER XL_barrierUpCurve,
																	   LPXLOPER XL_productsToPrice,
																	   bool PersistentInXL)
{
	ADD_LOG("Local_CRAQuantoCalculator_Common");
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
	
	//	General Datas
	//-----------------------------------------------------------------------------------------
		VECTOR<CCString> generalDatas;
		XL_readStrVector (XL_generalDatas, generalDatas, "ARM_ERR: General datas: array of string expected", DOUBLE_TYPE, C_result);
		if (generalDatas.size() < 5)
		{
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}

		CCString ccyDomStr = generalDatas[0];
		string sCcyDom = CCSTringToSTLString(ccyDomStr);
		ARM_Currency ccyDom(sCcyDom.c_str());

		CCString ccyForStr = generalDatas[1];
		string sCcyFor = CCSTringToSTLString(ccyForStr);
		ARM_Currency ccyFor(sCcyFor.c_str());

		CCString startDateXl = generalDatas[2];
		string sStartDate = CCSTringToSTLString(startDateXl);
		double startDate = atoi(sStartDate.c_str());

		CCString endDateXl = generalDatas[3];
		string sEndDate = CCSTringToSTLString(endDateXl);
		double endDate = atoi(sEndDate.c_str());
	    
		CCString payReceiveStr = generalDatas[4];
		long lPayReceive;
		if ((lPayReceive = ARM_ConvRecOrPay (payReceiveStr, C_result)) == ARM_DEFAULT_ERR)
		{
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}
		int payReceive = (int)lPayReceive;

	//	Call Datas
	//-----------------------------------------------------------------------------------------
		VECTOR<CCString> callDatas;
		XL_readStrVector (XL_callDatas, callDatas, "ARM_ERR: Call datas: array of string expected", DOUBLE_TYPE, C_result);
		if (callDatas.size() < 3)
		{
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}

		CCString callFreqXl = callDatas[0];
		long lCallFreq;
		if ((lCallFreq = ARM_ConvFrequency (callFreqXl, C_result)) == ARM_DEFAULT_ERR)
		{
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}
		int callFreq = (int)lCallFreq;

		CCString callNoticeXl = callDatas[1];
		string sCallNotice = CCSTringToSTLString(callNoticeXl);
		int callNotice = atoi(sCallNotice.c_str());

		//CallCal
		CCString callCalXl = callDatas[2];
		string callCal = CCSTringToSTLString(callCalXl);

		
	//	Funding Datas
	//-----------------------------------------------------------------------------------------
		VECTOR<CCString> fundDatas;
		XL_readStrVector (XL_fundDatas, fundDatas," ARM_ERR: Fund datas: array of string expected",DOUBLE_TYPE,C_result);
		if (fundDatas.size() < 2)
		{
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}

		CCString fundFreqXl = fundDatas[0];
		long lFundFreq;
		if ((lFundFreq = ARM_ConvFrequency (fundFreqXl, C_result)) == ARM_DEFAULT_ERR)
		{
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}
		int fundFreq = (int)lFundFreq;

		CCString fundDayCountXl = fundDatas[1];
		int fundDayCount = (int)ARM_ConvDayCount (fundDayCountXl);

	//	Coupon Datas
	//-----------------------------------------------------------------------------------------
		VECTOR<CCString> cpnDatas;
		XL_readStrVector (XL_cpnDatas, cpnDatas," ARM_ERR: Cpn datas: array of string expected",DOUBLE_TYPE,C_result);
		if (cpnDatas.size() < 9)
		{
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}

		/// standard
		CCString cpnDayCountXl = cpnDatas[0];
		int cpnDayCount = (int)ARM_ConvDayCount (cpnDayCountXl);

		CCString cpnPayFreqXl = cpnDatas[1];
		long lCpnPayFreq;
		if ((lCpnPayFreq = ARM_ConvFrequency (cpnPayFreqXl, C_result)) == ARM_DEFAULT_ERR)
		{
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}
		int cpnPayFreq = (int)lCpnPayFreq;

		/// corridor
		CCString cpnResetFreqXl = cpnDatas[2];
		long lCpnResetFreq;
		if ((lCpnResetFreq = ARM_ConvFrequency (cpnResetFreqXl, C_result)) == ARM_DEFAULT_ERR)
		{
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}
		int cpnResetFreq = (int)lCpnResetFreq;

		//Cpn Reset Cal
		CCString cpnResetCalXl = cpnDatas[3];
		string 	cpnResetCal = CCSTringToSTLString(cpnResetCalXl);

		//Cpn Pay Cal
		CCString cpnPayCalXl = cpnDatas[4];
		string 	cpnPayCal = CCSTringToSTLString(cpnPayCalXl);

		CCString cpnResetTimingXl = cpnDatas[5];
		int cpnResetTiming = (int)ARM_ConvPayResetRule(cpnResetTimingXl);

		int refIndex = ARM_ConvIrType(cpnDatas[6]);
		
		/// pay index
		int payIndex = K_FIXED;
		payIndex = ARM_ConvIrType(cpnDatas[7]);
		
		int payIndexResetTiming = K_ADVANCE;
		CCString payIndexResetTimingXl = cpnDatas[8];
		payIndexResetTiming = (int)ARM_ConvPayResetRule(payIndexResetTimingXl);
		
		/// and other
		double notional;
		CCString notionalStr;
		long notionalId;
		XL_readStrOrNumCell(XL_notionalCurve, notionalStr, notional, notionalId,
			   " ARM_ERR: notional: numeric or refValue Id expected",C_result);	
		if (notionalId == XL_TYPE_STRING)
			notionalId = LocalGetNumObjectId(notionalStr);
		else
			notionalId = ARM_NULL_OBJECT;

		CCString callFeesStr;
		long callFeesId;
		XL_readStrCell(XL_callFeesCurve, callFeesStr," ARM_ERR: Call Fees: String Id expected", C_result);
		callFeesId = LocalGetNumObjectId(callFeesStr);

		double fundSpread;
		CCString fundSpreadStr;
		long fundSpreadId;
		XL_readStrOrNumCell(XL_fundSpreadCurve, fundSpreadStr, fundSpread, fundSpreadId,
			   " ARM_ERR: fundSpread: numeric or refValue Id expected",C_result);	
		if (fundSpreadId == XL_TYPE_STRING)
			fundSpreadId = LocalGetNumObjectId(fundSpreadStr);
		else
			fundSpreadId = ARM_NULL_OBJECT;

		double fix;
		CCString fixStr;
		long fixId;
		XL_readStrOrNumCell(XL_fixCurve, fixStr, fix, fixId,
			   " ARM_ERR: boosted fix rate: numeric or refValue Id expected",C_result);	
		if (fixId == XL_TYPE_STRING)
			fixId = LocalGetNumObjectId(fixStr);
		else
			fixId = ARM_NULL_OBJECT;

		double bDown;
		CCString bDownStr;
		long bDownId;
		XL_readStrOrNumCell(XL_barrierDownCurve, bDownStr, bDown, bDownId,
			   " ARM_ERR: BarrierDown: numeric or refValue Id expected",C_result);	
		if (bDownId == XL_TYPE_STRING)
			bDownId = LocalGetNumObjectId(bDownStr);
		else
			bDownId = ARM_NULL_OBJECT;

		double bUp;
		CCString bUpStr;
		long bUpId;
		XL_readStrOrNumCell(XL_barrierUpCurve, bUpStr, bUp, bUpId,
			   " ARM_ERR: BarrierUp: numeric or refValue Id expected",C_result);	
		if (bUpId == XL_TYPE_STRING)
			bUpId = LocalGetNumObjectId(bUpStr);
		else
			bUpId = ARM_NULL_OBJECT;

		//Pricing flags
		VECTOR<CCString> pricingFlags;
		XL_readStrVector(XL_productsToPrice,pricingFlags," ARM_ERR: Pricing flags: array of string expected",DOUBLE_TYPE,C_result);
		vector< string > productsToPrice(pricingFlags.size());
		for (size_t i=0; i<pricingFlags.size(); i++)
			productsToPrice[i] = CCSTringToSTLString(pricingFlags[i]);

		exportFunc31Args<
					ARM_Currency, ARM_Currency, double, double, int, 
					int, int, string, int, int, 
					int, int, string, string, int, 
					int, int, int, int, double, 
					long, long, double, long, double, 
					long, double, long, double, long, 
					vector<string>	>
			ourFunc(
					ccyDom,	ccyFor,	startDate, endDate,	payReceive,	
					callFreq, callNotice, callCal,fundFreq, fundDayCount,
					cpnDayCount, cpnPayFreq, cpnResetCal, cpnPayCal, cpnResetFreq,
					cpnResetTiming, refIndex, payIndex, payIndexResetTiming, notional, 
					notionalId, callFeesId, fundSpread, fundSpreadId, fix, 
					fixId, bDown, bDownId, bUp, bUpId, 
					productsToPrice,
					ARMLOCAL_CRAQuantoCalculator_Create );			

		bool PersistentInXL = true;

		/// call the general function
		fillXL_Result( LOCAL_GC_CRAQUANTO_CLASS, ourFunc, C_result, XL_result, PersistentInXL );
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_CRAQuantoCalculator_Create" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


__declspec(dllexport) LPXLOPER WINAPI Local_CRAQuantoCalculator_Create(
	LPXLOPER XL_generalDatas,
	LPXLOPER XL_callDatas,
	LPXLOPER XL_fundDatas,
	LPXLOPER XL_cpnDatas,
	LPXLOPER XL_notionalCurve,
	LPXLOPER XL_callFeesCurve,
	LPXLOPER XL_fundSpreadCurve,
	LPXLOPER XL_fixCurve,
	LPXLOPER XL_barrierDownCurve,
	LPXLOPER XL_barrierUpCurve,
	LPXLOPER XL_productsToPrice)
{
	ADD_LOG("Local_CRAQuantoCalculator_Create");
	bool PersistentInXL = true;
	
	return Local_CRAQuantoCalculator_Common(
			XL_generalDatas,
			XL_callDatas,
			XL_fundDatas,
			XL_cpnDatas,
			XL_notionalCurve,
			XL_callFeesCurve,
			XL_fundSpreadCurve,
			XL_fixCurve,
			XL_barrierDownCurve,
			XL_barrierUpCurve,
			XL_productsToPrice,
			PersistentInXL);
}

__declspec(dllexport) LPXLOPER WINAPI Local_PXL_CRAQuantoCalculator_Create(
	LPXLOPER XL_generalDatas,
	LPXLOPER XL_callDatas,
	LPXLOPER XL_fundDatas,
	LPXLOPER XL_cpnDatas,
	LPXLOPER XL_notionalCurve,
	LPXLOPER XL_callFeesCurve,
	LPXLOPER XL_fundSpreadCurve,
	LPXLOPER XL_fixCurve,
	LPXLOPER XL_barrierDownCurve,
	LPXLOPER XL_barrierUpCurve,
	LPXLOPER XL_productsToPrice)
{
	ADD_LOG("Local_PXL_CRAQuantoCalculator_Create");
	bool PersistentInXL = false;
			
	return Local_CRAQuantoCalculator_Common(
			XL_generalDatas,
			XL_callDatas,
			XL_fundDatas,
			XL_cpnDatas,
			XL_notionalCurve,
			XL_callFeesCurve,
			XL_fundSpreadCurve,
			XL_fixCurve,
			XL_barrierDownCurve,
			XL_barrierUpCurve,
			XL_productsToPrice,
			PersistentInXL);
}

__declspec(dllexport) LPXLOPER WINAPI Local_FXVanillaCalculator_CreateFromSecurity(
	LPXLOPER XL_SecurityId,
	LPXLOPER XL_BasketType,
	LPXLOPER XL_DigitType,
    LPXLOPER XL_VanillaType)
{
	ADD_LOG("Local_FXVanillaCalculator_CreateFromSecurity");

	/// to remove memory leak put a result holder!
	static XLOPER_Holder XL_resultHolder;
	XLOPER& XL_result = XL_resultHolder.GetResult();

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

	/// Security Id
	long C_FxStripId;
	XL_GETOBJID( XL_SecurityId, C_FxStripId, " ARM_ERR: Security: Object expected",	C_result);

	/// Domestic DayCount
	CCString C_BasketType;
	XL_readStrCell( XL_BasketType, C_BasketType, " ARM_ERR: Basket Type : string expected", C_result);

	/// Domestic Freq
	CCString C_DigitType;
	XL_readStrCell( XL_DigitType, C_DigitType, " ARM_ERR: Digital Type : string expected", C_result);

	/// Foreign DayCount
	CCString C_VanillaType;
	XL_readStrCell( XL_VanillaType, C_VanillaType, " ARM_ERR: Vanilla Type : string expected", C_result);

	long retCode;
	long objId;
	CCString prevClass;
	
	CCString curClass = LOCAL_GC_FXVANILLA_CLASS;
	CCString stringId = GetLastCurCellEnvValue ();
	
	if (!stringId)
	{
		retCode = ARMLOCAL_FXVanillaCalculator_CreateFromSecurity(	C_FxStripId,
																	CCSTringToSTLString(C_BasketType),
																	CCSTringToSTLString(C_DigitType), 
																	CCSTringToSTLString(C_VanillaType),
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
			retCode = ARMLOCAL_FXVanillaCalculator_CreateFromSecurity(	C_FxStripId,
																		CCSTringToSTLString(C_BasketType),
																		CCSTringToSTLString(C_DigitType), 
																		CCSTringToSTLString(C_VanillaType), 
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
			FreeCurCellContent();

			retCode = ARMLOCAL_FXVanillaCalculator_CreateFromSecurity(	C_FxStripId,
																		CCSTringToSTLString(C_BasketType),
																		CCSTringToSTLString(C_DigitType), 
																		CCSTringToSTLString(C_VanillaType),
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
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_FXVanillaCalculator_CreateFromSecurity" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}

__declspec(dllexport) LPXLOPER WINAPI Local_PXL_FXVanillaCalculator_CreateFromSecurity(
	LPXLOPER XL_SecurityId,
	LPXLOPER XL_BasketType,
	LPXLOPER XL_DigitType,
    LPXLOPER XL_VanillaType)
{
	ADD_LOG("Local_FXVanillaCalculator_CreateFromSecurity");

	/// to remove memory leak put a result holder!
	static XLOPER_Holder XL_resultHolder;
	XLOPER& XL_result = XL_resultHolder.GetResult();

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

	/// Security Id
	long C_FxStripId;
	XL_GETOBJID( XL_SecurityId, C_FxStripId, " ARM_ERR: Security: Object expected",	C_result);

	/// Domestic DayCount
	CCString C_BasketType;
	XL_readStrCell( XL_BasketType, C_BasketType, " ARM_ERR: Basket Type : string expected", C_result);

	/// Domestic Freq
	CCString C_DigitType;
	XL_readStrCell( XL_DigitType, C_DigitType, " ARM_ERR: Digital Type : string expected", C_result);

	/// Foreign DayCount
	CCString C_VanillaType;
	XL_readStrCell( XL_VanillaType, C_VanillaType, " ARM_ERR: Vanilla Type : string expected", C_result);

	long retCode;
	long objId;

	CCString curClass = LOCAL_GC_FXVANILLA_CLASS;
	CCString stringId;
	
	retCode = ARMLOCAL_FXVanillaCalculator_CreateFromSecurity(	C_FxStripId,
																CCSTringToSTLString(C_BasketType),
																CCSTringToSTLString(C_DigitType), 
																CCSTringToSTLString(C_VanillaType),
																C_result);

	if (retCode == ARM_OK)
	{
		objId	 = C_result.getLong ();
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
	
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_FXVanillaCalculator_CreateFromSecurity" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


// VOLBOND CALCULATOR

class volbondcalculatorFunctor : public ARMResultLong2LongFunc
{
public:
	volbondcalculatorFunctor(
        long calculatorId,
        long dataToSetId,
		const vector< string >& mktDataKeys,
        bool isUpdated)
    :
    C_calculatorId(calculatorId),
    C_dataToSetId(dataToSetId),
	C_mktDataKeys(mktDataKeys),
    C_isUpdated(isUpdated)
    {};
	
	long operator()( ARM_result& result, long objId ){
		return ARMLOCAL_GC_Set(
            C_calculatorId,
            C_dataToSetId,
			C_mktDataKeys,
            C_isUpdated,
            result,
            objId);
	}

private:
	long			C_calculatorId;
	long			C_dataToSetId;
	vector<string>	C_mktDataKeys;
    bool			C_isUpdated;
};

__declspec(dllexport) LPXLOPER WINAPI Local_VolBondCalculator_Create(
	LPXLOPER XL_NominalStartEndDate,
	LPXLOPER XL_PayFreq,
	LPXLOPER XL_ResetFreq,	
	LPXLOPER XL_DayCount,
	LPXLOPER XL_Tenor,	
	LPXLOPER XL_IntRule,
	LPXLOPER XL_StubRule,
	LPXLOPER XL_ResetGap,	
	LPXLOPER XL_PayCalendar,
	LPXLOPER XL_ResetCalendar,
	LPXLOPER XL_OdeSolvers,
	LPXLOPER XL_RKParameters,	
	LPXLOPER XL_MCParameters,
	LPXLOPER XL_RandomGenerator,
	LPXLOPER XL_PayOffType,	
	LPXLOPER XL_MarketDataManager,
	LPXLOPER XL_MarketDataManagerKeys,
	LPXLOPER XL_ProductsToPrice)
{
	ADD_LOG("Local_VolBondCalculator_Create");

	/// to remove memory leak put a result holder!
	static XLOPER_Holder XL_resultHolder;
	XLOPER& XL_result = XL_resultHolder.GetResult();

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
		
		vector<double> datageneral;
		XL_readNumVector(XL_NominalStartEndDate,datageneral," ARM_ERR: RKParameters: array of objects expected",C_result);
		
		if (datageneral.size() != 3)
			throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"need Nominal, Start Date, End Date");

		double nominal = datageneral[0];
		double startDate = datageneral[1];
		double endDate = datageneral[2];

		CCString payfreq;
		XL_readStrCell(XL_PayFreq,payfreq," ARM_ERR: pay freq: string expected",C_result);
		CCString resetfreq;
		XL_readStrCell(XL_ResetFreq, resetfreq," ARM_ERR: resetfreq: string expected",C_result);

		CCString daycount;
		XL_readStrCell(XL_DayCount,daycount," ARM_ERR: daycount: string expected",C_result);
		CCString tenor;
		XL_readStrCell(XL_Tenor,tenor," ARM_ERR: tenor: string expected",C_result);

		CCString intrule;
		XL_readStrCell(XL_IntRule, intrule," ARM_ERR: intrule: string expected",C_result);
		CCString stubrule;
		XL_readStrCell(XL_StubRule,stubrule," ARM_ERR: stubrule: string expected",C_result);

		double resetgap;
		XL_readNumCell(XL_ResetGap,resetgap," ARM_ERR: resetgap: double expected",C_result);

		CCString paycalendar;
		XL_readStrCell(XL_PayCalendar,paycalendar," ARM_ERR: paycalendar: string expected",C_result);
		CCString resetcalendar;
		XL_readStrCell(XL_ResetCalendar,resetcalendar," ARM_ERR: resetcalendar: string expected",C_result);

		CCString type_odesolver;
		XL_readStrCell(XL_OdeSolvers,type_odesolver," ARM_ERR: type odesolver: string expected",C_result);

		vector<double> RKParameters;
		XL_readNumVector(XL_RKParameters,RKParameters," ARM_ERR: RKParameters: array of objects expected",C_result);
		
		if (RKParameters.size() != 9)
			throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"9 = 2(RK5) + 7(RK4) needed parameters for OdeSolver");

		vector<double> RK4Parameters;
		vector<double> RK5Parameters;

		size_t i = 0;
		for (; i < 2; ++i)
			RK5Parameters.push_back(RKParameters[i]);
		for (; i < RKParameters.size(); ++i)
			RK4Parameters.push_back(RKParameters[i]);


		vector<double> mcparameters;
		XL_readNumVector(XL_MCParameters,mcparameters," ARM_ERR: mcparameters: array of objects expected",C_result);

		if (mcparameters.size() != 3)
			throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"need MCnbsteps, MCBucket, MCNbStepsByYear");
		
		double MCnbsteps = mcparameters[0];
		double MCBucket = mcparameters[1];
		double MCNbStepsByYear = mcparameters[2];

		vector<CCString> CCrandomgenerator;
		vector<string> randomgenerator;
		XL_readStrVector(XL_RandomGenerator,CCrandomgenerator," ARM_ERR: randomgenerator: array of objects expected",DOUBLE_TYPE,C_result);
		for (i = 0; i < CCrandomgenerator.size(); ++i)
			randomgenerator.push_back(string(CCrandomgenerator[i].c_str()));

		CCString payofftype;
		XL_readStrCell( XL_PayOffType, payofftype, " ARM_ERR: PayOffType : string expected", C_result);

		CCString marketdatamanager;
		XL_readStrCell(XL_MarketDataManager,marketdatamanager," ARM_ERR: market data manager: string expected",C_result);
		long marketdatamanager_id = LocalGetNumObjectId (marketdatamanager);

		vector<CCString> CCmarketdatamanagerkeys;
		vector<string> marketdatamanagerkeys;
		XL_readStrVector(XL_MarketDataManagerKeys,CCmarketdatamanagerkeys," ARM_ERR: market data manager keys: array of objects expected",DOUBLE_TYPE,C_result);
		for (i = 0; i < CCmarketdatamanagerkeys.size(); ++i)
			marketdatamanagerkeys.push_back(string(CCmarketdatamanagerkeys[i].c_str()));

		vector<CCString> CCproductstoprice;
		vector<string> productstoprice;
		XL_readStrVector(XL_ProductsToPrice,CCproductstoprice," ARM_ERR: products to price: array of objects expected",DOUBLE_TYPE,C_result);
		for (i = 0; i < CCproductstoprice.size(); ++i)
			productstoprice.push_back(string(CCproductstoprice[i].c_str()));


		long retCode;		

		CCString curClass = LOCAL_GC_VOLBOND_CLASS;
		CCString stringId;
		
		retCode = ARMLOCAL_VolBondCalculator_Create(	nominal, 
														startDate, 
														endDate,
														CCSTringToSTLString(payfreq),
														CCSTringToSTLString(resetfreq), 
														CCSTringToSTLString(daycount),
														CCSTringToSTLString(tenor),
														CCSTringToSTLString(intrule),
														CCSTringToSTLString(stubrule),
														resetgap,
														CCSTringToSTLString(paycalendar),
														CCSTringToSTLString(resetcalendar),
														CCSTringToSTLString(type_odesolver),
														RK4Parameters,
														RK5Parameters,
														MCnbsteps,
														MCBucket,
														MCNbStepsByYear,
														randomgenerator,														
														marketdatamanager_id,
														marketdatamanagerkeys,
														CCSTringToSTLString(payofftype),
														productstoprice,
														C_result);

	

		if (retCode == ARM_OK)
		{
			long objId	 = C_result.getLong ();
			LocalSetCurCellEnvValue (curClass, objId); 
			stringId = LocalMakeObjectId (objId, curClass);
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_VolBondCalculator_Create" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}
