/*! \file ARM_xl_xxxproject_local.cpp
 *
 *  \brief file for the XXX project
 *
 *	\author  R. Guillemot
 *	\version 1.0
 *	\date March 2004
 */

#include <ARM\libarm_local\firstToBeIncluded.h>
#include <libCCxll\CCxll.h>
#include <GP_Infra\gpinfra\gramfunctorarg.h>
#include <ARM\libarm_local\ARM_local_xxxproject.h>
#include <ARM\libarm_local\ARM_local_class.h>
#include "ARM_xl_gp_common.h"
#include "ARM_xl_wrapper_local.h"
#include "ARM_xl_trycatch_local.h"

#include "ARM_local_interglob.h"
#include "ARM_local_interface.h"
#include "ARM_xl_trycatch_local.h"
#include "ARM_xl_gp_fctorhelper.h"

#include "util\tech_macro.h"

using ARM::ARM_GramFctorArg;


LPXLOPER Local_MktDatas_Create_Common(
	LPXLOPER XL_asOfDate,
	LPXLOPER XL_mktDatasKeys,
	LPXLOPER XL_mktDatasId,
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
		XL_readNumCell( XL_asOfDate, asOfDate, " ARM_ERR: As Of Date : Double expected",C_result);

		VECTOR<CCString> mktDatasKeys;
		vector<string> mktDatasKeysSTL;
		XL_readStrVector (XL_mktDatasKeys,mktDatasKeys," ARM_ERR: Market datas keys : array of string expected",DOUBLE_TYPE,C_result);

		VECTOR<CCString> mktDatasIdStr;
		vector<long> mktDatasId; 
		XL_readStrVector (XL_mktDatasId,mktDatasIdStr," ARM_ERR: Market datas keys : array of string expected",DOUBLE_TYPE,C_result);

		mktDatasKeysSTL.resize(mktDatasKeys.size());
		mktDatasId.resize(mktDatasIdStr.size());
		for (int i=0; i < mktDatasIdStr.size(); ++i)
		{	
			mktDatasKeysSTL[i] = CCSTringToSTLString(mktDatasKeys[i]);
			mktDatasId[i] = LocalGetNumObjectId(mktDatasIdStr[i]);
		}

		/// use the concept of Functor to transfer the knowledge of
		/// a function with a context
		exportFunc3Args<double,vector<string>,vector<long> > ourFunc(
			asOfDate,
			mktDatasKeysSTL,
			mktDatasId,
			ARMLOCAL_MktDatas_Create);

		/// call the general function
		fillXL_Result( LOCAL_XXX_MKTDATAS_CLASS, ourFunc, C_result, XL_result, PersistentInXL );

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
__declspec(dllexport) LPXLOPER WINAPI Local_MktDatas_Create(
	LPXLOPER XL_asOfDate,
	LPXLOPER XL_mktDatasKeys,
	LPXLOPER XL_mktDatasId)
{
	ADD_LOG("Local_MktDatas_Create");
	bool PersistentInXL = true;
	return Local_MktDatas_Create_Common( XL_asOfDate, XL_mktDatasKeys, XL_mktDatasId, PersistentInXL );
}


///////////////////////////////////
/// version for non persistentInXl
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_PXL_MktDatas_Create(
	LPXLOPER XL_asOfDate,
	LPXLOPER XL_mktDatasKeys,
	LPXLOPER XL_mktDatasId)
{
	ADD_LOG("Local_PXL_MktDatas_Create");
	bool PersistentInXL = false;
	return Local_MktDatas_Create_Common( XL_asOfDate, XL_mktDatasKeys, XL_mktDatasId, PersistentInXL );
}


// price des XXX object ( security + market data)
__declspec(dllexport) LPXLOPER WINAPI Local_ARM_XXX_Price (	
	LPXLOPER XL_secId,
	LPXLOPER XL_mktDt)
{
	ADD_LOG("Local_ARM_XXX_Price ");
	///	ARM_BEGIN();
	// return
	static XLOPER	XL_result;
	ARM_result		C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
		ARM_NOCALCIFWIZ();

		// C variable
		CCString C_secId;
		CCString C_mktDt;
		long mktDt;
		
		// error
		static int		error;
		static char*	reason = "";

		XL_readStrCell(XL_secId,C_secId," ARM_ERR: security id: object expected",C_result);
		XL_readStrCellWD(XL_mktDt,C_mktDt,"DEFAULT"," ARM_XXX_ERR: model id: object expected",C_result);

		if(C_mktDt == "DEFAULT")	mktDt = ARM_NULL_OBJECT;
		else						mktDt = LocalGetNumObjectId (C_mktDt);

		long retCode = ARMLOCAL_ARM_XXX_Price (LocalGetNumObjectId (C_secId), mktDt, C_result);

		if(retCode == ARM_OK)
		{
			FreeCurCellErr ();
			XL_result.xltype	= xltypeNum;
			XL_result.val.num	= C_result.getDouble ();
		}
		else{
			ARM_ERR();
		}
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_ARM_Price" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


// establish the XXW hedge scenario ( serie of criteras)
__declspec(dllexport) LPXLOPER WINAPI Local_ARM_XXX_Hedge_Create (	LPXLOPER XL_secId,
																	LPXLOPER XL_scenId,
																	LPXLOPER XL_mktDt)
{
	ADD_LOG("Local_ARM_XXX_Hedge_Create ");
	///	ARM_BEGIN();
	// return
	static XLOPER	XL_result;
	ARM_result		C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
		ARM_NOCALCIFWIZ();

		// C variable
		CCString	C_secId;
		CCString	C_scenId;
		CCString	C_mktDt;
		long		secId;
		long		scenId;
		long		mktDt;

		// error
		static int		error;
		static char*	reason = "";

		XL_readStrCellWD(XL_secId,	C_secId,	"DEFAULT"," ARM_XXX_ERR: model id: object expected",C_result);
		XL_readStrCellWD(XL_scenId,	C_scenId,	"DEFAULT"," ARM_XXX_ERR: model id: object expected",C_result);
		XL_readStrCellWD(XL_mktDt,	C_mktDt,	"DEFAULT"," ARM_XXX_ERR: model id: object expected",C_result);

		if(C_secId == "DEFAULT")	secId	= ARM_NULL_OBJECT;
		else						secId	= LocalGetNumObjectId (C_secId);
		if(C_scenId== "DEFAULT")	scenId	= ARM_NULL_OBJECT;
		else						scenId	= LocalGetNumObjectId (C_scenId);
		if(C_mktDt == "DEFAULT")	mktDt	= ARM_NULL_OBJECT;
		else						mktDt	= LocalGetNumObjectId (C_mktDt);

		
		exportFunc3Args<long, long, long > ourFunc(
			secId,
			scenId,
			mktDt,
			ARMLOCAL_Hedge_Create);

		/// call the general function
		fillXL_Result_withName( ourFunc, C_result, XL_result, true );

	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_ARM_Price" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}



__declspec(dllexport) LPXLOPER WINAPI Local_ARM_XXX_Hedge_GetData (	LPXLOPER XL_HedgeId,
																	LPXLOPER XL_Key)
{
	ADD_LOG("Local_ARM_XXX_Hedge_GetData ");
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

				// C variable
		CCString	CC_hedgeId;
		CCString	CC_Key;

		long		hedgeId;
		
		XL_readStrCellWD(XL_HedgeId, CC_hedgeId, "DEFAULT"," ARM_XXX_ERR: hedge id: object expected",C_result);
		if(CC_hedgeId == "DEFAULT")	hedgeId	= ARM_NULL_OBJECT;
		else						hedgeId	= LocalGetNumObjectId (CC_hedgeId);
	
		XL_readStrCell	(XL_Key, CC_Key," ARM_XXX_ERR: string expected",C_result);
		
		/// call the function
		ARM_GramFctorArg argResult;
		
		long retCode = ARMLOCAL_Hedge_GetData(hedgeId,	CCSTringToSTLString( CC_Key ),
															argResult,	C_result );
		
		if( retCode == ARM_KO )
		{
			ARM_ERR();
		}
		else
		{
			ARM_GramFunctorToXLOPER(argResult, XL_result, C_result, true);
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


// establish the XXW hedge scenario ( serie of criteras)
__declspec(dllexport) LPXLOPER WINAPI Local_ARM_XXX_Scenario (	LPXLOPER XL_Shift,
																LPXLOPER XL_Currency,
																LPXLOPER XL_Type_Scenario,
																LPXLOPER XL_SubType_Scenario,
																LPXLOPER XL_Stress_Order,
																LPXLOPER XL_Relative,
																LPXLOPER XL_CumulInv,
																LPXLOPER XL_Perturbative
																)
{
	ADD_LOG("Local_ARM_XXX_Scenario ");
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

		CCString	CC_Currency;
		string		currency;
		XL_readStrCell (XL_Currency, CC_Currency," ARM_ERR: string expected",C_result);
		currency		=	CCSTringToSTLString(CC_Currency);

		CCString	CC_Type_Scenario;
		string		typeScenario;
		XL_readStrCell (XL_Type_Scenario, CC_Type_Scenario," ARM_ERR: string expected",C_result);
		typeScenario	=	CCSTringToSTLString(CC_Type_Scenario);

		CCString	CC_SubType_Scenario;
		string		subTypeScenario;
		XL_readStrCellWD (XL_SubType_Scenario, CC_SubType_Scenario,""," ARM_ERR: string expected",C_result);
		subTypeScenario	=	CCSTringToSTLString(CC_SubType_Scenario);

		CCString	CC_Stress_Order;
		string		stressOrder;
		XL_readStrCellWD (XL_Stress_Order, CC_Stress_Order,""," ARM_ERR: string expected",C_result);
		stressOrder		=	CCSTringToSTLString(CC_Stress_Order);

		CCString	CC_Relative;
		XL_readStrCellWD (XL_Relative, CC_Relative, "N"," ARM_ERR: expected (Y/N)",C_result);
		long relative = (CC_Relative == "Y" || CC_Relative == "YES") ? K_YES : K_NO ;

		CCString		CC_CumulInv;
	    XL_readStrCellWD(XL_CumulInv, CC_CumulInv, "N", " ARM_ERR: expected (Y/N)",C_result);
		long cumulInv = (CC_CumulInv == "Y" || CC_CumulInv == "YES") ? K_YES : K_NO ;

		CCString		CC_Perturbative;
	    XL_readStrCellWD(XL_Perturbative, CC_Perturbative, "N", " ARM_ERR: expected (Y/N)",C_result);
		long perturbative = (CC_Perturbative == "Y" || CC_Perturbative == "YES") ? K_YES : K_NO ;

		double		shift;
		XL_readNumCell (XL_Shift, shift," ARM_ERR: double expected",C_result);



		/// use the concept of Functor to transfer the knowledge of
		/// a function with a context

		exportFunc8Args<double, string, string, string, string, long, long , long > ourFunc(
			shift,
			currency,
			typeScenario,
			subTypeScenario,
			stressOrder,
			relative,
			cumulInv,
			perturbative,
			ARMLOCAL_ARM_XXX_Scenario);

		/// call the general function
		fillXL_Result_withName( ourFunc, C_result, XL_result, true );
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_Arm_Scenario" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


__declspec(dllexport) LPXLOPER WINAPI Local_ARM_XXX_Scenari_Compose(	LPXLOPER XL_scenId1,
																		LPXLOPER XL_scenId2)
{
	ADD_LOG("Local_ARM_XXX_Scenari_Compose");
	static XLOPER	XL_result;
	ARM_result		C_result;

	ARM_XL_TRY_BLOCK_BEGIN
	{
		ARM_NOCALCIFWIZ();

		// C variable
		CCString		C_scenId1;
		CCString		C_scenId2;
		long			scenId1;
		long			scenId2;

		static int		error;
		static char*	reason = "";

		XL_readStrCellWD(XL_scenId1,	C_scenId1,	"DEFAULT"," ARM_XXX_ERR: model id: object expected",C_result);
		XL_readStrCellWD(XL_scenId2,	C_scenId2,	"DEFAULT"," ARM_XXX_ERR: model id: object expected",C_result);

		if(C_scenId1 == "DEFAULT")	scenId1		= ARM_NULL_OBJECT;
		else						scenId1		= LocalGetNumObjectId (C_scenId1);

		if(C_scenId2  == "DEFAULT")	scenId2		= ARM_NULL_OBJECT;
		else						scenId2		= LocalGetNumObjectId (C_scenId2);
		
		exportFunc2Args<long, long > ourFunc(
			scenId1,
			scenId2,
			ARMLOCAL_Scenari_Compose);

		/// call the general function
		fillXL_Result_withName( ourFunc, C_result, XL_result, true );

	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Local_ARM_Price" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}

