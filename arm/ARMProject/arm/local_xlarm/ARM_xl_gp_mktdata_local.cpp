/*
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 * $Log: ARM_xl_gp_mktdata_local.cpp,v $
 * Revision 1.1  2004/02/07 15:08:43  ebenhamou
 * Initial version
 *
 */

/*! \file ARM_xl_gp_mktdata_local.cpp
 *
 *  \brief file for the mkt data part in the generic pricer
 *
 *	\author  E Benhamou
 *	\version 1.0
 *	\date March 2004
 */

#include <ARM\libarm_local\firstToBeIncluded.h>
#include <libCCxll\CCxll.h>
#include <ARM\libarm_local\ARM_local_gp_mktdata.h>
#include "ARM_xl_gp_mktdata_local.h"
#include "ARM_xl_wrapper_local.h"
#include "ARM_xl_trycatch_local.h"
#include <glob\dates.h>	/// for ARM_DATE

#include "util\tech_macro.h"

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
///             mkt data manager
/// Inputs : as of Date
///----------------------------------------------
///----------------------------------------------

////////////////////////////////////////////////
/// very rapid definition of ourlonglongFunctor
////////////////////////////////////////////////


class MktDataManagerFunc : public ARMResultLong2LongFunc
{
public:

	MktDataManagerFunc( const ARM_Date& asOf)
	:	C_asOf(asOf)
	{};
	
	long operator()( ARM_result& result, long objId ){
		return ARMLOCAL_MktDataManager_Create( C_asOf, result, objId);			
	}
private:
	ARM_Date C_asOf;
};


/////////////////////////////////////////////////////////////
/// central function that does the creation of the XL function
/////////////////////////////////////////////////////////////
LPXLOPER Local_MktDataManager_Common(
	LPXLOPER XL_asOf,
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

		double C_asOfDble;
		XL_readNumCell(	 XL_asOf,		C_asOfDble,		" ARM_ERR: asOfDate: date expected",				C_result);
		ARM_Date C_asOf = ConvertToARMDATE(C_asOfDble);
    
		/// use the concept of Functor to transfer the knowledge of
		/// a function with a context
		MktDataManagerFunc ourFunc(C_asOf);

		/// call the general function
		fillXL_Result( LOCAL_MKTDATAMANAGER_CLASS, ourFunc, C_result, XL_result, PersistentInXL );

	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_MktDataManager_Common" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}
							 

///////////////////////////////////
/// version that takes into account 
/// previous creation of object
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_MktDataManager_Create(
	LPXLOPER XL_asOf )
{
	ADD_LOG("Local_MktDataManager_Create");
	bool PersistentInXL = true;
	return Local_MktDataManager_Common( XL_asOf, PersistentInXL );
}


///////////////////////////////////
/// version for non persistentInXl
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_PXL_MktDataManager_Create(
	LPXLOPER XL_asOf )
{
	ADD_LOG("Local_PXL_MktDataManager_Create");
	bool PersistentInXL = false;
	return Local_MktDataManager_Common( XL_asOf, PersistentInXL );
}




//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
///
/// #    #  #    #   #####
/// ##  ##  #   #      #
/// # ## #  ####       #
/// #    #  #  #       #
/// #    #  #   #      #
/// #    #  #    #     #
///
///
///
/// #####     ##     #####    ##
/// #    #   #  #      #     #  #
/// #    #  #    #     #    #    #
/// #    #  ######     #    ######
/// #    #  #    #     #    #    #
/// #####   #    #     #    #    #
///
///
///   ##     ####    ####   ######   ####    ####    ####   #####
///  #  #   #    #  #    #  #       #       #       #    #  #    #
/// #    #  #       #       #####    ####    ####   #    #  #    #
/// ######  #       #       #            #       #  #    #  #####
/// #    #  #    #  #    #  #       #    #  #    #  #    #  #   #
/// #    #   ####    ####   ######   ####    ####    ####   #    #
///
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////



////////////////////////////////////////////////
/// ZC CURVE
/// functor to get and set zc curve
////////////////////////////////////////////////

class ZCCurveGetSetFunc : public ARMResultLong2LongFunc
{
public:

	ZCCurveGetSetFunc(
		long mktDataManagerId,
		const string& indexName,
		const string& ccy,
		const string& cvName, 
		const ARM_Date& asOf,
		const string& source,
		long zcCurveId = ARM_NULL_OBJECT_ID )
	:
		C_mktDataManagerId(mktDataManagerId),
		C_indexName(indexName),
		C_ccy(ccy),
		C_cvName(cvName),
		C_asOf(asOf),
		C_source(source),
		C_zcCurveId(zcCurveId)
	{};

	long operator()( ARM_result& result, long objId = ARM_NULL_OBJECT_ID )
	{
		if( ARM_NULL_OBJECT_ID == C_zcCurveId )
			return ARMLOCAL_MktDataManager_GetZCCurve( 
				C_mktDataManagerId,
				C_indexName, 
				C_ccy, 
				C_cvName, 
				C_asOf,
				C_source, 
				result, 
				objId);	
		else
			return ARMLOCAL_MktDataManager_RegisterZCCurve(
				C_mktDataManagerId,
				C_indexName,
				C_ccy,
				C_cvName,
				C_asOf,
				C_source,
				C_zcCurveId,
				result );
	}
			
private:
	long C_mktDataManagerId;
	string C_indexName;
	string C_ccy;
	string C_cvName;
	ARM_Date C_asOf;
	string C_source; 
	long C_zcCurveId;
};


/////////////////////////////////////////////////////////////
/// central function that does the creation of the XL function
/////////////////////////////////////////////////////////////
LPXLOPER Local_MktDataManager_ZCCurveGetSet_Common(
	LPXLOPER XL_mktDataManagerId,
	LPXLOPER XL_indexName,
	LPXLOPER XL_ccy,
	LPXLOPER XL_cvName,
	LPXLOPER XL_asOf,
	LPXLOPER XL_source,
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

		long C_mktDataManagerId;
		XL_GETOBJID( XL_mktDataManagerId, C_mktDataManagerId,	" ARM_ERR: Mkt Data Manager id: object expected",	C_result);

		CCString C_indexNameCCStr;
		XL_readStrCell(XL_indexName,C_indexNameCCStr," ARM_ERR: index name: String expected",C_result);
		string C_indexName = CCSTringToSTLString(C_indexNameCCStr);

		CCString C_ccyCCStr;
		XL_readStrCell(XL_ccy,C_ccyCCStr," ARM_ERR: ccy: String expected",C_result);
		string C_ccy = CCSTringToSTLString(C_ccyCCStr);

		CCString C_cvNameCCStr;
		XL_readStrCell(XL_cvName,C_cvNameCCStr," ARM_ERR: curve name: String expected",C_result);
		string C_cvName = CCSTringToSTLString(C_cvNameCCStr);

		double C_asOfDble;
		XL_readNumCell(	 XL_asOf,		C_asOfDble,		" ARM_ERR: asOfDate: date expected",				C_result);
		ARM_Date C_asOf = ConvertToARMDATE(C_asOfDble);

		CCString C_sourceCCStr;
		XL_readStrCell(XL_source,C_sourceCCStr," ARM_ERR: source: String expected",C_result);
		string C_source = CCSTringToSTLString(C_sourceCCStr);

		long C_zcCurveId;
		if( XL_zcCurveId == NULL )
		{
			C_zcCurveId = ARM_NULL_OBJECT_ID;
		}
		else
		{
			XL_GETOBJID( XL_zcCurveId, C_zcCurveId,	" ARM_ERR: zero curve id: object expected",	C_result);
		}

		/// use the concept of Functor to transfer the knowledge of
		/// a function with a context
		ZCCurveGetSetFunc ourFunc(
			C_mktDataManagerId,
			C_indexName,
			C_ccy,
			C_cvName, 
			C_asOf,
			C_source,
			C_zcCurveId );

		if( XL_zcCurveId == NULL )
		{
			/// this is a zero curve
			fillXL_Result( LOCAL_ZERO_CURVE_LIN_CLASS, ourFunc, C_result, XL_result, PersistentInXL );
		}
		else
		{
			/// this is a market data manager
			fillXL_Result( LOCAL_MKTDATAMANAGER_CLASS, ourFunc, C_result, XL_result, PersistentInXL );
		}
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_MktDataManager_ZCCurveGetSet_Common" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


///////////////////////////////////
/// version that takes into account 
/// previous creation of object
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_MktDataManager_ZCCurveGet(
	LPXLOPER XL_mktDataManagerId,
	LPXLOPER XL_indexName,
	LPXLOPER XL_ccy,
	LPXLOPER XL_cvName,
	LPXLOPER XL_asOf,
	LPXLOPER XL_source )
{
	ADD_LOG("Local_MktDataManager_ZCCurveGet");
	bool PersistentInXL = true;
	return Local_MktDataManager_ZCCurveGetSet_Common(
		XL_mktDataManagerId,
		XL_indexName,
		XL_ccy,
		XL_cvName,
		XL_asOf,
		XL_source,
		NULL,
		PersistentInXL );
}


///////////////////////////////////
/// version for non persistentInXl
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_PXL_MktDataManager_ZCCurveGet(
	LPXLOPER XL_mktDataManagerId,
	LPXLOPER XL_indexName,
	LPXLOPER XL_ccy,
	LPXLOPER XL_cvName,
	LPXLOPER XL_asOf,
	LPXLOPER XL_source )
{
	ADD_LOG("Local_PXL_MktDataManager_ZCCurveGet");
	bool PersistentInXL = false;
	return Local_MktDataManager_ZCCurveGetSet_Common(
		XL_mktDataManagerId,
		XL_indexName,
		XL_ccy,
		XL_cvName,
		XL_asOf,
		XL_source,
		NULL,
		PersistentInXL );
}


///////////////////////////////////
/// version that takes into account 
/// previous creation of object
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_MktDataManager_ZCCurveSet(
	LPXLOPER XL_mktDataManagerId,
	LPXLOPER XL_indexName,
	LPXLOPER XL_ccy,
	LPXLOPER XL_cvName,
	LPXLOPER XL_asOf,
	LPXLOPER XL_source,
	LPXLOPER XL_zcCurveId )
{
	ADD_LOG("Local_MktDataManager_ZCCurveSet");
	bool PersistentInXL = true;
	return Local_MktDataManager_ZCCurveGetSet_Common(
		XL_mktDataManagerId,
		XL_indexName,
		XL_ccy,
		XL_cvName,
		XL_asOf,
		XL_source,
		XL_zcCurveId ,
		PersistentInXL );
}


///////////////////////////////////
/// version for non persistentInXl
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_PXL_MktDataManager_ZCCurveSet(
	LPXLOPER XL_mktDataManagerId,
	LPXLOPER XL_indexName,
	LPXLOPER XL_ccy,
	LPXLOPER XL_cvName,
	LPXLOPER XL_asOf,
	LPXLOPER XL_source,
	LPXLOPER XL_zcCurveId )
{
	ADD_LOG("Local_PXL_MktDataManager_ZCCurveSet");
	bool PersistentInXL = false;
	return Local_MktDataManager_ZCCurveGetSet_Common(
		XL_mktDataManagerId,
		XL_indexName,
		XL_ccy,
		XL_cvName,
		XL_asOf,
		XL_source,
		XL_zcCurveId,
		PersistentInXL );
}



/////////////////////////////////////////////////////////////
/// central function that does the creation of the XL function
/////////////////////////////////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_MktDataManager_ZCCurveGetKey(
	LPXLOPER XL_indexName,
	LPXLOPER XL_ccy,
	LPXLOPER XL_cvName,
	LPXLOPER XL_asOf,
	LPXLOPER XL_source )
{
	ADD_LOG("Local_MktDataManager_ZCCurveGetKey");
	
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

		CCString C_indexNameCCStr;
		XL_readStrCell(XL_indexName,C_indexNameCCStr," ARM_ERR: index name: String expected",C_result);
		string C_indexName = CCSTringToSTLString(C_indexNameCCStr);

		CCString C_ccyCCStr;
		XL_readStrCell(XL_ccy,C_ccyCCStr," ARM_ERR: ccy: String expected",C_result);
		string C_ccy = CCSTringToSTLString(C_ccyCCStr);

		CCString C_cvNameCCStr;
		XL_readStrCell(XL_cvName,C_cvNameCCStr," ARM_ERR: curve name: String expected",C_result);
		string C_cvName = CCSTringToSTLString(C_cvNameCCStr);

		double C_asOfDble;
		XL_readNumCell(	 XL_asOf,		C_asOfDble,		" ARM_ERR: asOfDate: date expected",				C_result);
		ARM_Date C_asOf = ConvertToARMDATE(C_asOfDble);

		CCString C_sourceCCStr;
		XL_readStrCell(XL_source,C_sourceCCStr," ARM_ERR: source: String expected",C_result);
		string C_source = CCSTringToSTLString(C_sourceCCStr);

		long retCode = ARMLOCAL_MktDataManager_GetZCCurveKey(
			C_indexName,
			C_ccy,
			C_cvName, 
			C_asOf,
			C_source,
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_MktDataManager_ZCCurveGetKey" )

	return (LPXLOPER)&XL_result;
}



////////////////////////////////////////////////
/// Vol CURVE
/// functor to get and set vol curve
////////////////////////////////////////////////

class VolCurveGetSetFunc : public ARMResultLong2LongFunc
{
public:

	VolCurveGetSetFunc(
		long mktDataManagerId,
		const string& indexName,
		const string& ccy,
		const string& cvName, 
		const ARM_Date& asOf,
		const string& volMktType,
		const string& volType,
		const string& source,
		long volCurveId = ARM_NULL_OBJECT_ID )
	:
		C_mktDataManagerId(mktDataManagerId),
		C_indexName(indexName),
		C_ccy(ccy),
		C_cvName(cvName),
		C_asOf(asOf),
		C_volMktType(volMktType),
		C_volType(volType),
		C_source(source),
		C_volCurveId(volCurveId)
	{};

	long operator()( ARM_result& result, long objId = ARM_NULL_OBJECT_ID  )
	{
		if( ARM_NULL_OBJECT_ID == C_volCurveId )
			return ARMLOCAL_MktDataManager_GetVolCurve( 
				C_mktDataManagerId,
				C_indexName,
				C_ccy,
				C_cvName,
				C_asOf,
				C_volMktType,
				C_volType,
				C_source,
				result,
				objId );
		else
			return ARMLOCAL_MktDataManager_RegisterVolCurve(
				C_mktDataManagerId,
				C_indexName,
				C_ccy,
				C_cvName,
				C_asOf,
				C_volMktType,
				C_volType,
				C_source,
				C_volCurveId,
				result );
	}
private:
	long C_mktDataManagerId;
	string C_indexName;
	string C_ccy;
	string C_cvName;
	ARM_Date C_asOf;
	string C_volMktType;
	string C_volType;
	string C_source;
	long C_volCurveId;
};



/////////////////////////////////////////////////////////////
/// central function that does the creation of the XL function
/////////////////////////////////////////////////////////////
LPXLOPER Local_MktDataManager_VolCurveGetSet_Common(
	LPXLOPER XL_mktDataManagerId,
	LPXLOPER XL_indexName,
	LPXLOPER XL_ccy,
	LPXLOPER XL_cvName,
	LPXLOPER XL_asOf,
	LPXLOPER XL_volMktType,
	LPXLOPER XL_volType,
	LPXLOPER XL_source,
	LPXLOPER XL_volCurveId,
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

		long C_mktDataManagerId;
		XL_GETOBJID( XL_mktDataManagerId, C_mktDataManagerId,	" ARM_ERR: Mkt Data Manager id: object expected",	C_result);

		CCString C_indexNameCCStr;
		XL_readStrCell(XL_indexName,C_indexNameCCStr," ARM_ERR: index name: String expected",C_result);
		string C_indexName = CCSTringToSTLString(C_indexNameCCStr);

		CCString C_ccyCCStr;
		XL_readStrCell(XL_ccy,C_ccyCCStr," ARM_ERR: ccy: String expected",C_result);
		string C_ccy = CCSTringToSTLString(C_ccyCCStr);

		CCString C_cvNameCCStr;
		XL_readStrCell(XL_cvName,C_cvNameCCStr," ARM_ERR: curve name: String expected",C_result);
		string C_cvName = CCSTringToSTLString(C_cvNameCCStr);

		double C_asOfDble;
		XL_readNumCell(	 XL_asOf,		C_asOfDble,		" ARM_ERR: asOfDate: date expected",				C_result);
		ARM_Date C_asOf = ConvertToARMDATE(C_asOfDble);

		CCString C_volMktTypeCCStr;
		XL_readStrCell(XL_volMktType,C_volMktTypeCCStr," ARM_ERR: vol mkt type: String expected",C_result);
		string C_volMktType = CCSTringToSTLString(C_volMktTypeCCStr);

		CCString C_volTypeCCStr;
		XL_readStrCell(XL_volType,C_volTypeCCStr," ARM_ERR: vol type: String expected",C_result);
		string C_volType= CCSTringToSTLString(C_volTypeCCStr);

		CCString C_sourceCCStr;
		XL_readStrCell(XL_source,C_sourceCCStr," ARM_ERR: source: String expected",C_result);
		string C_source = CCSTringToSTLString(C_sourceCCStr);

		long C_volCurveId;
		if( XL_volCurveId == NULL )
		{
			C_volCurveId = ARM_NULL_OBJECT_ID;
		}
		else
		{
			XL_GETOBJID( XL_volCurveId, C_volCurveId,	" ARM_ERR: vol curve id: object expected",	C_result);
		}

		/// use the concept of Functor to transfer the knowledge of
		/// a function with a context
		VolCurveGetSetFunc ourFunc(
			C_mktDataManagerId,
			C_indexName,
			C_ccy,
			C_cvName, 
			C_asOf,
			C_volMktType,
			C_volType,
			C_source,
			C_volCurveId );

		if( XL_volCurveId == NULL )
		{
			/// this is a vol curve
			fillXL_Result( LOCAL_VOL_CURVE_LIN_CLASS, ourFunc, C_result, XL_result, PersistentInXL );
		}
		else
		{
			/// this is a market data manager
			fillXL_Result( LOCAL_MKTDATAMANAGER_CLASS, ourFunc, C_result, XL_result, PersistentInXL );
		}
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_MktDataManager_VolCurveGetSet_Common" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


///////////////////////////////////
/// version that takes into account 
/// previous creation of object
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_MktDataManager_VolCurveGet(
	LPXLOPER XL_mktDataManagerId,
	LPXLOPER XL_indexName,
	LPXLOPER XL_ccy,
	LPXLOPER XL_cvName,
	LPXLOPER XL_asOf,
	LPXLOPER XL_volMktType,
	LPXLOPER XL_volType,
	LPXLOPER XL_source )
{
	ADD_LOG("Local_MktDataManager_VolCurveGet");
	bool PersistentInXL = true;
	return Local_MktDataManager_VolCurveGetSet_Common(
		XL_mktDataManagerId,
		XL_indexName,
		XL_ccy,
		XL_cvName,
		XL_asOf,
		XL_volMktType,
		XL_volType,
		XL_source,
		NULL,
		PersistentInXL );
}


///////////////////////////////////
/// version for non persistentInXl
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_PXL_MktDataManager_VolCurveGet(
	LPXLOPER XL_mktDataManagerId,
	LPXLOPER XL_indexName,
	LPXLOPER XL_ccy,
	LPXLOPER XL_cvName,
	LPXLOPER XL_asOf,
	LPXLOPER XL_volMktType,
	LPXLOPER XL_volType,
	LPXLOPER XL_source )
{
	ADD_LOG("Local_PXL_MktDataManager_VolCurveGet");
	bool PersistentInXL = false;
	return Local_MktDataManager_VolCurveGetSet_Common(
		XL_mktDataManagerId,
		XL_indexName,
		XL_ccy,
		XL_cvName,
		XL_asOf,
		XL_volMktType,
		XL_volType,
		XL_source,
		NULL,
		PersistentInXL );
}


///////////////////////////////////
/// version that takes into account 
/// previous creation of object
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_MktDataManager_VolCurveSet(
	LPXLOPER XL_mktDataManagerId,
	LPXLOPER XL_indexName,
	LPXLOPER XL_ccy,
	LPXLOPER XL_cvName,
	LPXLOPER XL_asOf,
	LPXLOPER XL_volMktType,
	LPXLOPER XL_volType,
	LPXLOPER XL_source,
	LPXLOPER XL_volCurveId
)
{
	ADD_LOG("Local_MktDataManager_VolCurveSet");
	bool PersistentInXL = true;
	return Local_MktDataManager_VolCurveGetSet_Common(
		XL_mktDataManagerId,
		XL_indexName,
		XL_ccy,
		XL_cvName,
		XL_asOf,
		XL_volMktType,
		XL_volType,
		XL_source,
		XL_volCurveId,
		PersistentInXL );
}


///////////////////////////////////
/// version for non persistentInXl
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_PXL_MktDataManager_VolCurveSet(
	LPXLOPER XL_mktDataManagerId,
	LPXLOPER XL_indexName,
	LPXLOPER XL_ccy,
	LPXLOPER XL_cvName,
	LPXLOPER XL_asOf,
	LPXLOPER XL_volMktType,
	LPXLOPER XL_volType,
	LPXLOPER XL_source,
	LPXLOPER XL_volCurveId
)
{
	ADD_LOG("Local_PXL_MktDataManager_VolCurveSet");
	bool PersistentInXL = false;
	return Local_MktDataManager_VolCurveGetSet_Common(
		XL_mktDataManagerId,
		XL_indexName,
		XL_ccy,
		XL_cvName,
		XL_asOf,
		XL_volMktType,
		XL_volType,
		XL_source,
		XL_volCurveId,
		PersistentInXL );
}




/////////////////////////////////////////////////////////////
/// central function that does the creation of the XL function
/////////////////////////////////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_MktDataManager_VolCurveGetKey(
	LPXLOPER XL_indexName,
	LPXLOPER XL_ccy,
	LPXLOPER XL_cvName,
	LPXLOPER XL_asOf,
	LPXLOPER XL_volMktType,
	LPXLOPER XL_volType,
	LPXLOPER XL_source )
{
	ADD_LOG("Local_MktDataManager_VolCurveGetKey");
	
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

		CCString C_indexNameCCStr;
		XL_readStrCell(XL_indexName,C_indexNameCCStr," ARM_ERR: index name: String expected",C_result);
		string C_indexName = CCSTringToSTLString(C_indexNameCCStr);

		CCString C_ccyCCStr;
		XL_readStrCell(XL_ccy,C_ccyCCStr," ARM_ERR: ccy: String expected",C_result);
		string C_ccy = CCSTringToSTLString(C_ccyCCStr);

		CCString C_cvNameCCStr;
		XL_readStrCell(XL_cvName,C_cvNameCCStr," ARM_ERR: curve name: String expected",C_result);
		string C_cvName = CCSTringToSTLString(C_cvNameCCStr);

		double C_asOfDble;
		XL_readNumCell(	 XL_asOf,		C_asOfDble,		" ARM_ERR: asOfDate: date expected",				C_result);
		ARM_Date C_asOf = ConvertToARMDATE(C_asOfDble);

		CCString C_volMktTypeCCStr;
		XL_readStrCell(XL_volMktType,C_volMktTypeCCStr," ARM_ERR: vol mkt type: String expected",C_result);
		string C_volMktType = CCSTringToSTLString(C_volMktTypeCCStr);

		CCString C_volTypeCCStr;
		XL_readStrCell(XL_volType,C_volTypeCCStr," ARM_ERR: vol type: String expected",C_result);
		string C_volType= CCSTringToSTLString(C_volTypeCCStr);

		CCString C_sourceCCStr;
		XL_readStrCell(XL_source,C_sourceCCStr," ARM_ERR: source: String expected",C_result);
		string C_source = CCSTringToSTLString(C_sourceCCStr);

		long retCode = ARMLOCAL_MktDataManager_GetVolCurveKey(
			C_indexName,
			C_ccy,
			C_cvName, 
			C_asOf,
			C_source,
			C_volMktType,
			C_volType,
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_MktDataManager_VolCurveGetKey" )

	return (LPXLOPER)&XL_result;
}



////////////////////////////////////////////////
/// Vol Mkt Model
/// functor to get and set vol mkt model
////////////////////////////////////////////////

class VolMktModelGetSetFunc : public ARMResultLong2LongFunc
{
public:

	VolMktModelGetSetFunc (
		long mktDataManagerId,
		const string& indexName,
		const string& ccy,
		const string& cvName, 
		const ARM_Date& asOf,
		const string& volMktType,
		const string& source,
		long volMktModelId = ARM_NULL_OBJECT_ID )
	:
		C_mktDataManagerId(mktDataManagerId),
		C_indexName(indexName),
		C_ccy(ccy),
		C_cvName(cvName),
		C_asOf(asOf),
		C_volMktType(volMktType),
		C_source(source),
		C_volMktModelId(volMktModelId)
	{};

	long operator()( ARM_result& result, long objId = ARM_NULL_OBJECT_ID )
	{
		if( ARM_NULL_OBJECT_ID == C_volMktModelId )
			return ARMLOCAL_MktDataManager_GetVolMktModel( 
				C_mktDataManagerId,
				C_indexName,
				C_ccy,
				C_cvName,
				C_asOf,
				C_volMktType,
				C_source,
				result,
				objId );
		else
			return ARMLOCAL_MktDataManager_RegisterVolMktModel(
				C_mktDataManagerId,
				C_indexName,
				C_ccy,
				C_cvName,
				C_asOf,
				C_volMktType,
				C_source,
				C_volMktModelId,
				result );
	}
private:
	long C_mktDataManagerId;
	string C_indexName;
	string C_ccy;
	string C_cvName;
	ARM_Date C_asOf;
	string C_volMktType;
	string C_source; 
	long C_volMktModelId;
};




/////////////////////////////////////////////////////////////
/// central function that does the creation of the XL function
/////////////////////////////////////////////////////////////
LPXLOPER Local_MktDataManager_VolMktModelGetSet_Common(
	LPXLOPER XL_mktDataManagerId,
	LPXLOPER XL_indexName,
	LPXLOPER XL_ccy,
	LPXLOPER XL_cvName,
	LPXLOPER XL_asOf,
	LPXLOPER XL_volMktType,
	LPXLOPER XL_source,
	LPXLOPER XL_volMktModelId,
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

		long C_mktDataManagerId;
		XL_GETOBJID( XL_mktDataManagerId, C_mktDataManagerId,	" ARM_ERR: Mkt Data Manager id: object expected",	C_result);

		CCString C_indexNameCCStr;
		XL_readStrCell(XL_indexName,C_indexNameCCStr," ARM_ERR: index name: String expected",C_result);
		string C_indexName = CCSTringToSTLString(C_indexNameCCStr);

		CCString C_ccyCCStr;
		XL_readStrCell(XL_ccy,C_ccyCCStr," ARM_ERR: ccy: String expected",C_result);
		string C_ccy = CCSTringToSTLString(C_ccyCCStr);

		CCString C_cvNameCCStr;
		XL_readStrCell(XL_cvName,C_cvNameCCStr," ARM_ERR: curve name: String expected",C_result);
		string C_cvName = CCSTringToSTLString(C_cvNameCCStr);

		double C_asOfDble;
		XL_readNumCell(	 XL_asOf,		C_asOfDble,		" ARM_ERR: asOfDate: date expected",				C_result);
		ARM_Date C_asOf = ConvertToARMDATE(C_asOfDble);

		CCString C_volMktTypeCCStr;
		XL_readStrCell(XL_volMktType,C_volMktTypeCCStr," ARM_ERR: vol mkt type: String expected",C_result);
		string C_volMktType = CCSTringToSTLString(C_volMktTypeCCStr);

		CCString C_sourceCCStr;
		XL_readStrCell(XL_source,C_sourceCCStr," ARM_ERR: source: String expected",C_result);
		string C_source = CCSTringToSTLString(C_sourceCCStr);

		long C_volMktModelId;
		if( XL_volMktModelId == NULL )
		{
			C_volMktModelId = ARM_NULL_OBJECT_ID;
		}
		else
		{
			XL_GETOBJID( XL_volMktModelId, C_volMktModelId,	" ARM_ERR: vol mkt model id: object expected",	C_result);
		}

		/// use the concept of Functor to transfer the knowledge of
		/// a function with a context
		VolMktModelGetSetFunc ourFunc(
			C_mktDataManagerId,
			C_indexName,
			C_ccy,
			C_cvName, 
			C_asOf,
			C_volMktType,
			C_source,
			C_volMktModelId );

		if( XL_volMktModelId == NULL )
		{
			/// this is a bs model
			fillXL_Result( LOCAL_BSMODEL_CLASS, ourFunc, C_result, XL_result, PersistentInXL );
		}
		else
		{
			/// this is a market data manager
			fillXL_Result( LOCAL_MKTDATAMANAGER_CLASS, ourFunc, C_result, XL_result, PersistentInXL );
		}

	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_MktDataManager_VolMktModelGetSet_Common" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


///////////////////////////////////
/// version that takes into account 
/// previous creation of object
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_MktDataManager_VolMktModelGet(
	LPXLOPER XL_mktDataManagerId,
	LPXLOPER XL_indexName,
	LPXLOPER XL_ccy,
	LPXLOPER XL_cvName,
	LPXLOPER XL_asOf,
	LPXLOPER XL_volMktType,
	LPXLOPER XL_source )
{
	ADD_LOG("Local_MktDataManager_VolMktModelGet");
	bool PersistentInXL = true;
	return Local_MktDataManager_VolMktModelGetSet_Common(
		XL_mktDataManagerId,
		XL_indexName,
		XL_ccy,
		XL_cvName,
		XL_asOf,
		XL_volMktType,
		XL_source,
		NULL,
		PersistentInXL );
}


///////////////////////////////////
/// version for non persistentInXl
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_PXL_MktDataManager_VolMktModelGet(
	LPXLOPER XL_mktDataManagerId,
	LPXLOPER XL_indexName,
	LPXLOPER XL_ccy,
	LPXLOPER XL_cvName,
	LPXLOPER XL_asOf,
	LPXLOPER XL_volMktType,
	LPXLOPER XL_source )
{
	ADD_LOG("Local_PXL_MktDataManager_VolMktModelGet");
	bool PersistentInXL = false;
	return Local_MktDataManager_VolMktModelGetSet_Common(
		XL_mktDataManagerId,
		XL_indexName,
		XL_ccy,
		XL_cvName,
		XL_asOf,
		XL_volMktType,
		XL_source,
		NULL,
		PersistentInXL );
}

///////////////////////////////////
/// version that takes into account 
/// previous creation of object
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_MktDataManager_VolMktModelSet(
	LPXLOPER XL_mktDataManagerId,
	LPXLOPER XL_indexName,
	LPXLOPER XL_ccy,
	LPXLOPER XL_cvName,
	LPXLOPER XL_asOf,
	LPXLOPER XL_volMktType,
	LPXLOPER XL_source,
	LPXLOPER XL_volMktModelId )
{
	ADD_LOG("Local_MktDataManager_VolMktModelSet");
	bool PersistentInXL = true;
	return Local_MktDataManager_VolMktModelGetSet_Common(
		XL_mktDataManagerId,
		XL_indexName,
		XL_ccy,
		XL_cvName,
		XL_asOf,
		XL_volMktType,
		XL_source,
		XL_volMktModelId,
		PersistentInXL );
}



///////////////////////////////////
/// version for non persistentInXl
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_PXL_MktDataManager_VolMktModelSet(
	LPXLOPER XL_mktDataManagerId,
	LPXLOPER XL_indexName,
	LPXLOPER XL_ccy,
	LPXLOPER XL_cvName,
	LPXLOPER XL_asOf,
	LPXLOPER XL_volMktType,
	LPXLOPER XL_source,
	LPXLOPER XL_volMktModelId )
{
	ADD_LOG("Local_PXL_MktDataManager_VolMktModelSet");
	bool PersistentInXL = false;
	return Local_MktDataManager_VolMktModelGetSet_Common(
		XL_mktDataManagerId,
		XL_indexName,
		XL_ccy,
		XL_cvName,
		XL_asOf,
		XL_volMktType,
		XL_source,
		XL_volMktModelId,
		PersistentInXL );
}


/////////////////////////////////////////////////////////////
/// central function that does the creation of the XL function
/////////////////////////////////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_MktDataManager_VolMktModelGetKey(
	LPXLOPER XL_indexName,
	LPXLOPER XL_ccy,
	LPXLOPER XL_cvName,
	LPXLOPER XL_asOf,
	LPXLOPER XL_volMktType,
	LPXLOPER XL_source )
{
	ADD_LOG("Local_MktDataManager_VolMktModelGetKey");
	
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

		CCString C_indexNameCCStr;
		XL_readStrCell(XL_indexName,C_indexNameCCStr," ARM_ERR: index name: String expected",C_result);
		string C_indexName = CCSTringToSTLString(C_indexNameCCStr);

		CCString C_ccyCCStr;
		XL_readStrCell(XL_ccy,C_ccyCCStr," ARM_ERR: ccy: String expected",C_result);
		string C_ccy = CCSTringToSTLString(C_ccyCCStr);

		CCString C_cvNameCCStr;
		XL_readStrCell(XL_cvName,C_cvNameCCStr," ARM_ERR: curve name: String expected",C_result);
		string C_cvName = CCSTringToSTLString(C_cvNameCCStr);

		double C_asOfDble;
		XL_readNumCell(	 XL_asOf,		C_asOfDble,		" ARM_ERR: asOfDate: date expected",				C_result);
		ARM_Date C_asOf = ConvertToARMDATE(C_asOfDble);

		CCString C_volMktTypeCCStr;
		XL_readStrCell(XL_volMktType,C_volMktTypeCCStr," ARM_ERR: vol mkt type: String expected",C_result);
		string C_volMktType = CCSTringToSTLString(C_volMktTypeCCStr);

		CCString C_sourceCCStr;
		XL_readStrCell(XL_source,C_sourceCCStr," ARM_ERR: source: String expected",C_result);
		string C_source = CCSTringToSTLString(C_sourceCCStr);

		long retCode = ARMLOCAL_MktDataManager_GetVolMktModelKey(
			C_indexName,
			C_ccy,
			C_cvName, 
			C_asOf,
			C_volMktType,
			C_source,

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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_MktDataManager_VolMktModelGetKey" )

	return (LPXLOPER)&XL_result;
}



////////////////////////////////////////////////
/// Mkt Data general
/// functor to get and set data
////////////////////////////////////////////////

class DataGetSetFunc : public ARMResultLong2LongFunc
{
public:

	DataGetSetFunc (
		long mktDataManagerId,
		const string& objectKey,
        bool isFill = false,
		long objectId = ARM_NULL_OBJECT_ID )
	:
		C_mktDataManagerId(mktDataManagerId),
		C_objectKey(objectKey),
        C_isFill(isFill),
		C_objectId(objectId)
	{};

	long operator()( ARM_result& result, long objId = ARM_NULL_OBJECT_ID)
	{
		if( ARM_NULL_OBJECT_ID == C_objectId )
			return ARMLOCAL_MktDataManager_GetData( 
				C_mktDataManagerId,
				C_objectKey,
				result,
				objId );
		else
			return ARMLOCAL_MktDataManager_RegisterData(
				C_mktDataManagerId,
				C_objectKey,
				C_objectId,
                C_isFill,
				result,
                objId );
	}
private:
	long C_mktDataManagerId;
	string C_objectKey;
	long C_objectId;
    bool C_isFill;
};



/////////////////////////////////////////////////////////////
/// central function that does the creation of the XL function
/////////////////////////////////////////////////////////////
LPXLOPER Local_MktDataManager_MktDataGetSet_Common(
	LPXLOPER XL_mktDataManagerId,
	LPXLOPER XL_objectKey,
	LPXLOPER XL_objectId,
    bool isFill,
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

		long C_mktDataManagerId;
		XL_GETOBJID( XL_mktDataManagerId, C_mktDataManagerId,	" ARM_ERR: Mkt Data Manager id: object expected",	C_result);

		CCString C_objectKeyCCStr;
		XL_readStrCell(XL_objectKey,C_objectKeyCCStr," ARM_ERR: object key: String expected",C_result);
		string C_objectKey = CCSTringToSTLString(C_objectKeyCCStr);

		long C_objectId;
		if( XL_objectId == NULL )
		{
			C_objectId = ARM_NULL_OBJECT_ID;
		}
		else
		{
			XL_GETOBJID( XL_objectId, C_objectId,	" ARM_ERR: object id: object expected",	C_result);
		}

		/// use the concept of Functor to transfer the knowledge of
		/// a function with a context
		DataGetSetFunc ourFunc(
			C_mktDataManagerId,
			C_objectKey,
            isFill,
			C_objectId );

		if( XL_objectId == NULL )
		{
			/// this is a get!
			/// call the general function
			/// because the object can be of any type, we use the previous class to handle persistance!
			fillXL_Result( LOCAL_ANY_CLASS,	ourFunc, C_result, XL_result, PersistentInXL );
		}
		else
		{
            if(isFill)
            {
                /// Simple updating of the calculator
		        long retCode = ourFunc(C_result);

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
                /// this is a market data manager
                fillXL_Result( LOCAL_MKTDATAMANAGER_CLASS, ourFunc, C_result, XL_result, PersistentInXL );
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


///////////////////////////////////
/// version that takes into account 
/// previous creation of object
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_MktDataManager_MktDataGet(
	LPXLOPER XL_mktDataManagerId,
	LPXLOPER XL_objectKey )
{
	ADD_LOG("Local_MktDataManager_MktDataGet");
	bool PersistentInXL = true;
	return Local_MktDataManager_MktDataGetSet_Common(
		XL_mktDataManagerId,
		XL_objectKey,
		NULL,
        false,
		PersistentInXL );
}


///////////////////////////////////
/// version for non persistentInXl
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_PXL_MktDataManager_MktDataGet(
	LPXLOPER XL_mktDataManagerId,
	LPXLOPER XL_objectKey )
{
	ADD_LOG("Local_PXL_MktDataManager_MktDataGet");
	bool PersistentInXL = false;
	return Local_MktDataManager_MktDataGetSet_Common(
		XL_mktDataManagerId,
		XL_objectKey,
		NULL,
        false,
		PersistentInXL );
}


///////////////////////////////////
/// version that takes into account 
/// previous creation of object
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_MktDataManager_MktDataSet(
	LPXLOPER XL_mktDataManagerId,
	LPXLOPER XL_objectKey,
	LPXLOPER XL_objectId )
{
	ADD_LOG("Local_MktDataManager_MktDataSet");
	bool PersistentInXL = true;
    bool isFill = false;
	return Local_MktDataManager_MktDataGetSet_Common(
		XL_mktDataManagerId,
		XL_objectKey,
		XL_objectId,
        isFill,
		PersistentInXL );
}


///////////////////////////////////
/// version for non persistentInXl
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_PXL_MktDataManager_MktDataSet(
	LPXLOPER XL_mktDataManagerId,
	LPXLOPER XL_objectKey,
	LPXLOPER XL_objectId )
{
	ADD_LOG("Local_PXL_MktDataManager_MktDataSet");
	bool PersistentInXL = false;
    bool isFill = false;
	return Local_MktDataManager_MktDataGetSet_Common(
		XL_mktDataManagerId,
		XL_objectKey,
		XL_objectId,
        isFill,
		PersistentInXL );
}

///////////////////////////////////
/// version that takes into account 
/// previous creation of object and just 
/// fills a new Mkt Data and no cloned object
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_MktDataManager_MktDataFill(
	LPXLOPER XL_mktDataManagerId,
	LPXLOPER XL_objectKey,
	LPXLOPER XL_objectId )
{
	ADD_LOG("Local_MktDataManager_MktDataFill");
	bool PersistentInXL = true;
    bool isFill = true;
	return Local_MktDataManager_MktDataGetSet_Common(
		XL_mktDataManagerId,
		XL_objectKey,
		XL_objectId,
        isFill,
		PersistentInXL );
}

/////////////////////////////////////////////////////////////
/// central function that does the creation of the XL function
/////////////////////////////////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_MktDataManager_SetDetailMode(
	LPXLOPER XL_mktDataManagerId,
	LPXLOPER XL_detailMode )
{
	ADD_LOG("Local_MktDataManager_SetDetailMode");
	
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

		long C_mktDataManagerId;
		XL_GETOBJID( XL_mktDataManagerId, C_mktDataManagerId,	" ARM_ERR: Mkt Data Manager id: object expected",	C_result);

		double C_detailMode;
		XL_readNumCell( XL_detailMode, C_detailMode,	" ARM_ERR: detail mode Flag: boolean compatible expected",	C_result);
		bool C_detailModeFlag = C_detailMode != 0;
			
		long retCode = ARMLOCAL_MktDataManager_SetDetailMode(
				C_mktDataManagerId,
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_MktDataManager_SetDetailMode" )
		
	return (LPXLOPER)&XL_result;
}


/////////////////////////////////////////////////////////////
/// central function that does the creation of the XL function
/////////////////////////////////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_MktDataManager_ResetAllData(
	LPXLOPER XL_mktDataManagerId )
{
	ADD_LOG("Local_MktDataManager_ResetAllData");
	
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

		long C_mktDataManagerId;
		XL_GETOBJID( XL_mktDataManagerId, C_mktDataManagerId,	" ARM_ERR: Mkt Data Manager id: object expected",	C_result);

		long retCode = ARMLOCAL_MktDataManager_ResetAllData(
				C_mktDataManagerId,
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_MktDataManager_ResetAllData" )

	return (LPXLOPER)&XL_result;
}




////////////////////////////////////////////////
/// function to return an empty mkt data manager
/////////////////////////////////////////////////////////////
class ResetMyDataFunc: public ARMResultLong2LongFunc
{
public:
	ResetMyDataFunc( long mktDataManagerId )
	:	C_mktDataManagerId(mktDataManagerId)
	{};

	long operator()( ARM_result& result, long objId = ARM_NULL_OBJECT_ID )
	{
		return ARMLOCAL_MktDataManager_ResetMyData( C_mktDataManagerId, result, objId);
	}

private:
	long C_mktDataManagerId;
};



/////////////////////////////////////////////////////////////
/// central function that does the creation of the XL function
/////////////////////////////////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_MktDataManager_ResetMyData(
	LPXLOPER XL_mktDataManagerId,
	bool PersistentInXL )
{
	ADD_LOG("Local_MktDataManager_ResetMyData");
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

		long C_mktDataManagerId;
		XL_GETOBJID( XL_mktDataManagerId, C_mktDataManagerId,	" ARM_ERR: Mkt Data Manager id: object expected",	C_result);

		/// use the concept of Functor to transfer the knowledge of
		/// a function with a context
		ResetMyDataFunc ourFunc(C_mktDataManagerId);

		/// call the general function
		fillXL_Result( LOCAL_MKTDATAMANAGER_CLASS, ourFunc, C_result, XL_result, PersistentInXL );

	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_MktDataManager_ResetMyData" )

	return (LPXLOPER)&XL_result;
}



class MktDataChangeAsOfFunc : public ARMResultLong2LongFunc
{
public:

	MktDataChangeAsOfFunc(
		long mktDataManagerId,
		const ARM_Date& asOfDate )
	:
		C_mktDataManagerId(mktDataManagerId),
		C_asOfDate(asOfDate)
	{};

	long operator()( ARM_result& result, long objId )
	{
		return ARMLOCAL_MktDataManager_ChangeDate( 
			C_mktDataManagerId,
			C_asOfDate,
			result,
			objId );
	}
private:
	long C_mktDataManagerId;
	ARM_Date C_asOfDate;
};



/////////////////////////////////////////////////////////////
/// central function that does the creation of the XL function
/////////////////////////////////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_MktDataManager_ChangeAsOf_Common(
	LPXLOPER XL_mktDataManagerId,
	LPXLOPER XL_asOfDate,
	bool PersistentInXL )
{
	ADD_LOG("Local_MktDataManager_ChangeAsOf_Common");
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

		long C_mktDataManagerId;
		XL_GETOBJID( XL_mktDataManagerId, C_mktDataManagerId,	" ARM_ERR: Mkt Data Manager id: object expected",	C_result);

		double C_asOfDble;
		XL_readNumCell(	 XL_asOfDate,	C_asOfDble,				" ARM_ERR: asOfDate: date expected",				C_result);
		ARM_Date C_asOf = ConvertToARMDATE(C_asOfDble);

		MktDataChangeAsOfFunc ourFunc(C_mktDataManagerId,C_asOf);

		/// call the general function
		fillXL_Result( LOCAL_MKTDATAMANAGER_CLASS, ourFunc, C_result, XL_result, PersistentInXL );
		

	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_MktDataManager_ResetMyData" )

	return (LPXLOPER)&XL_result;
}


///////////////////////////////////
/// version that takes into account 
/// previous creation of object
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_MktDataManager_ChangeAsOf(
	LPXLOPER XL_mktDataManagerId,
	LPXLOPER XL_asOfDate )
{
	ADD_LOG("Local_MktDataManager_ChangeAsOf");
	bool PersistentInXL = true;
	return Local_MktDataManager_ChangeAsOf_Common(
		XL_mktDataManagerId,
		XL_asOfDate,
		PersistentInXL );
}


///////////////////////////////////
/// version for non persistentInXl
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_PXL_MktDataManager_ChangeAsOf(
	LPXLOPER XL_mktDataManagerId,
	LPXLOPER XL_asOfDate )
{
	ADD_LOG("Local_PXL_MktDataManager_ChangeAsOf");
	bool PersistentInXL = false;
	return Local_MktDataManager_ChangeAsOf_Common(
		XL_mktDataManagerId,
		XL_asOfDate,
		PersistentInXL );
}

