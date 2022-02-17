
#pragma warning(disable :4005 4786 4800)

#include <libCCxll\CCxll.h>

#include <ARM\libarm_local\ARM_local_glob.h>
#include <ARM\libarm_local\ARM_local_class.h>

#include <ARM\libicm_local\ICM_local_convertible.h>

#include "ARM_local_interface.h"
#include "ARM_local_interglob.h"

__declspec(dllexport) LPXLOPER WINAPI ARM_Credit_CBOption(LPXLOPER XL_FirstDateIn,
														  LPXLOPER XL_LastDateIn,
														  LPXLOPER XL_OptionTypeIn,
														  LPXLOPER XL_StrikeIn,
														  LPXLOPER XL_StrikeTypeIn,
														  LPXLOPER XL_AccruedOnExIn,
														  LPXLOPER XL_BarrierTypeIn,
														  LPXLOPER XL_BarrierStrikeIn,
														  LPXLOPER XL_BarrierStrikeTypeIn)
{
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	ARM_NOCALCIFWIZ();

	// C variable
	double FirstDate;
	double LastDate;

	CCString C_OptionType;
	double Strike;
	CCString C_StrikeType;
	CCString C_AccruedOnEx;
	CCString C_BarrierType;
	CCString C_BarrierStrikeType;
	double BarrierStrike;
	double BarrierStrike_default=0.0;

	// error
	static int error;
	static char* reason = "";
	
	// Read LPXLOPER
	XL_readNumCell(XL_FirstDateIn,FirstDate,"ARM_ERR: First date: numeric expected",C_result);
	XL_readNumCell(XL_LastDateIn,LastDate,"ARM_ERR: Last date: numeric expected",C_result);

	XL_readStrCell(XL_OptionTypeIn,C_OptionType,"ARM_ERR: Option Type: string expected",C_result);
	XL_readNumCell(XL_StrikeIn,Strike,"ARM_ERR: Strike: numeric expected",C_result);
	XL_readStrCellWD(XL_StrikeTypeIn,C_StrikeType,"BONDPRICE","ARM_ERR: Strike Type: string expected",C_result);
	XL_readStrCellWD(XL_AccruedOnExIn,C_AccruedOnEx,"ACC","ARM_ERR: Accrued On Exercise: string expected",C_result);

	XL_readStrCellWD(XL_BarrierTypeIn,C_BarrierType,"UPANDIN","ARM_ERR: Barrier Type: String expected",C_result);
	XL_readNumCellWD(XL_BarrierStrikeIn,BarrierStrike,BarrierStrike_default,"ARM_ERR: Barrier's strike: numeric expected",C_result);
	XL_readStrCellWD(XL_BarrierStrikeTypeIn,C_BarrierStrikeType,"CONVOVERACCVAL","ARM_ERR: Barrier Strike Type: String expected",C_result);
	
	if(C_OptionType=="CONVERSION")
		C_StrikeType="NBSHARES";

	int qType,qStrikeType,qAccrued,qBarrierType,qBarrierStrikeType;

	if((qType = ARM_ConvCBOptionType (C_OptionType, C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	if((qStrikeType	=ARM_ConvCBOptionStrikeType(C_StrikeType,C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	if((qAccrued	=ARM_ConvCBOptionAccruedOnEx(C_AccruedOnEx,C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	if((qBarrierType=ARM_ConvCBOptionBarrierType(C_BarrierType,C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	if((qBarrierStrikeType=ARM_ConvCBOptionBarrierStrikeType(C_BarrierStrikeType,C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	if((BarrierStrike==0.0)||(qBarrierType==0)||(qBarrierStrikeType==0))
	{
		BarrierStrike=0.0;
		qBarrierType=0;
		qBarrierStrikeType=0;
	};

	// Call liblocal functions
	long retCode;
	long objId;
	CCString prevClass;
	CCString curClass=LOCAL_CBOPTION_CLASS;
	CCString stringId=GetLastCurCellEnvValue();

	if(!stringId)
	{
		retCode=ICMLOCAL_CBOPTION(FirstDate,LastDate,(int)qType,Strike,(int)qStrikeType,(int)qAccrued,
			(int)qBarrierType,BarrierStrike,qBarrierStrikeType,C_result);

		if(retCode==ARM_OK)
		{
			objId=C_result.getLong();
			LocalSetCurCellEnvValue(curClass,objId);
			stringId=LocalMakeObjectId(objId,curClass);
		};
	}
	else
	{
		prevClass=LocalGetStringObjectClass(stringId);
		objId=LocalGetNumObjectId(stringId);

		if(curClass==prevClass)
		{
			retCode=ICMLOCAL_CBOPTION(FirstDate,LastDate,(int)qType,Strike,(int)qStrikeType,(int)qAccrued,
				(int)qBarrierType,BarrierStrike,qBarrierStrikeType,C_result,objId);

			if(retCode==ARM_OK)
			{
				LocalSetCurCellEnvValue(curClass,objId);
				stringId=LocalMakeObjectId(objId,curClass);
			};
		}
		else
		{
			FreeCurCellContent();
			retCode=ICMLOCAL_CBOPTION(FirstDate,LastDate,(int)qType,Strike,(int)qStrikeType,(int)qAccrued,
				(int)qBarrierType,BarrierStrike,qBarrierStrikeType,C_result);
			if(retCode==ARM_OK)
			{
				objId=C_result.getLong();
				LocalSetCurCellEnvValue(curClass,objId);
				stringId=LocalMakeObjectId(objId,curClass);
			};
		};
	};

	// return
	if(retCode==ARM_OK)
	{
		FreeCurCellErr();
		XL_result.xltype=xltypeStr;
		XL_result.val.str=XL_StrC2StrPascal(stringId);
		XL_result.xltype|=xlbitDLLFree;
	}
	else
	{
		ARM_ERR();
	};

	return (LPXLOPER)&XL_result;
};

__declspec(dllexport) LPXLOPER WINAPI ARM_Credit_Convertible(LPXLOPER XL_DefBondIn,
															 LPXLOPER XL_StockIn,
															 LPXLOPER XL_CBOptions,
															 LPXLOPER XL_YieldIn)
{
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	ARM_NOCALCIFWIZ();

	// C variable
	CCString C_DefBondId;
	CCString C_StockId;
	double	 YieldIn;
	double YieldIn_Default=0.;
	VECTOR<CCString> C_CBOptions;
	VECTOR<long> C_CBOptionsIds;

	// error
	static int error;
	static char* reason = "";
	
	// Read LPXLOPER
	XL_readStrCell(XL_DefBondIn,C_DefBondId,"ARM_ERR: Defaultable bond id: string expected",C_result);
	XL_readStrCell(XL_StockIn,C_StockId,"ARM_ERR: Stock id: string expected",C_result);
	XL_readStrVector(XL_CBOptions,C_CBOptions,"ARM_ERR: CB Options: array of string expected",XL_TYPE_STRING,C_result);
	XL_readNumCellWD(XL_YieldIn,YieldIn,YieldIn_Default,"ARM_ERR: Reference Yield: numeric expected",C_result);

	for(int i=0;i<C_CBOptions.size();i++)
	{
		C_CBOptionsIds.push_back(LocalGetNumObjectId(C_CBOptions[i]));
	};

	// Call liblocal functions
	long retCode;
	long objId;
	CCString prevClass;
	CCString curClass=LOCAL_CONVERTIBLE_CLASS;
	CCString stringId=GetLastCurCellEnvValue();

	if(!stringId)
	{
		retCode=ICMLOCAL_CONVERTIBLE(LocalGetNumObjectId(C_DefBondId),LocalGetNumObjectId(C_StockId),YieldIn,C_CBOptionsIds,C_result);

		if(retCode==ARM_OK)
		{
			objId=C_result.getLong();
			LocalSetCurCellEnvValue(curClass,objId);
			stringId=LocalMakeObjectId(objId,curClass);
		};
	}
	else
	{
		prevClass=LocalGetStringObjectClass(stringId);
		objId=LocalGetNumObjectId(stringId);

		if(curClass==prevClass)
		{
			retCode=ICMLOCAL_CONVERTIBLE(LocalGetNumObjectId(C_DefBondId),LocalGetNumObjectId(C_StockId),YieldIn,C_CBOptionsIds,C_result,objId);

			if(retCode==ARM_OK)
			{
				LocalSetCurCellEnvValue(curClass,objId);
				stringId=LocalMakeObjectId(objId,curClass);
			};
		}
		else
		{
			FreeCurCellContent();
			retCode=ICMLOCAL_CONVERTIBLE(LocalGetNumObjectId(C_DefBondId),LocalGetNumObjectId(C_StockId),YieldIn,C_CBOptionsIds,C_result);
			if(retCode==ARM_OK)
			{
				objId=C_result.getLong();
				LocalSetCurCellEnvValue(curClass,objId);
				stringId=LocalMakeObjectId(objId,curClass);
			};
		};
	};

	// return
	if(retCode==ARM_OK)
	{
		FreeCurCellErr();
		XL_result.xltype=xltypeStr;
		XL_result.val.str=XL_StrC2StrPascal(stringId);
		XL_result.xltype|=xlbitDLLFree;
	}
	else
	{
		ARM_ERR();
	};

	return (LPXLOPER)&XL_result;
};

__declspec(dllexport) LPXLOPER WINAPI ARM_Credit_CbCallableAssetSwap(LPXLOPER XL_ConvertibleIn,
																	 LPXLOPER XL_CallFeatures,
																	 LPXLOPER XL_SwapStartDate,
																	 LPXLOPER XL_Index,
																	 LPXLOPER XL_Spread,
																	 LPXLOPER XL_EarlyCallDate,
																	 LPXLOPER XL_RecallSpread,
																	 LPXLOPER XL_HasEuribidClause,
																	 LPXLOPER XL_CallIfITM)
{
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	ARM_NOCALCIFWIZ();

	// C variable
	CCString C_ConvertId;
	VECTOR<CCString> C_CBOptions;
	VECTOR<long> C_CBOptionsIds;
	double C_StartDate;
	CCString C_IndexId;
	double C_Spread;
	double C_EarlyCall;
	double C_EarlyCall_default=0.0;
	double C_RecallSpread;
	double C_RecallSpread_default=0.0;

	CCString C_HasEuribidClause;
	CCString C_HasEuribidClause_default="TRUE";
	CCString C_CallIfITM;
	CCString C_CallIfITM_default="TRUE";
	
	long	hasEuribidClause;
	long	CallIfITM;

	// error
	static int error;
	static char* reason = "";
	
	// Read LPXLOPER
	XL_readStrCell(XL_ConvertibleIn,C_ConvertId,"ARM_ERR: Convertible bond id: string expected",C_result);
	XL_readStrVector(XL_CallFeatures,C_CBOptions,"ARM_ERR: Call features: array of string expected",XL_TYPE_STRING,C_result);
	XL_readNumCell(XL_SwapStartDate,C_StartDate,"ARM_ERR: Swap start date: numeric expected",C_result);
	XL_readStrCell(XL_Index,C_IndexId,"ARM_ERR: Reference index id: string expected",C_result);
	XL_readNumCell(XL_Spread,C_Spread,"ARM_ERR: Asset swap spread: numeric expected",C_result);
	XL_readNumCellWD(XL_EarlyCallDate,C_EarlyCall,C_EarlyCall_default,"ARM_ERR: Early Call Date: numeric expected",C_result);
	XL_readNumCellWD(XL_RecallSpread,C_RecallSpread,C_RecallSpread_default,"ARM_ERR: Recall spread: numeric expected",C_result);
	XL_readStrCellWD(XL_HasEuribidClause,C_HasEuribidClause,C_HasEuribidClause_default,"ARM_ERR: Euribid clause: 'TRUE', 'FALSE' or empty expected",C_result);
	XL_readStrCellWD(XL_CallIfITM,C_CallIfITM,C_CallIfITM_default,"ARM_ERR: Call if in-the-money: 'TRUE', 'FALSE' or empty expected",C_result);

	for(int i=0;i<C_CBOptions.size();i++)
	{
		C_CBOptionsIds.push_back(LocalGetNumObjectId(C_CBOptions[i]));
	};

	C_CBOptions.clear();

	if((hasEuribidClause = ARM_ConvStringToBool(C_HasEuribidClause,C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}
	if((CallIfITM = ARM_ConvStringToBool(C_CallIfITM,C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	// Call liblocal functions
	long retCode;
	long objId;
	CCString prevClass;
	CCString curClass=LOCAL_CBASW_CLASS;
	CCString stringId=GetLastCurCellEnvValue();

	if(!stringId)
	{
		retCode=ICMLOCAL_CBCALLABLEASW(LocalGetNumObjectId(C_ConvertId),C_CBOptionsIds,C_StartDate,LocalGetNumObjectId(C_IndexId),C_Spread,C_EarlyCall,C_RecallSpread,(bool)hasEuribidClause,(bool)CallIfITM,C_result);
		C_CBOptionsIds.clear();

		if(retCode==ARM_OK)
		{
			objId=C_result.getLong();
			LocalSetCurCellEnvValue(curClass,objId);
			stringId=LocalMakeObjectId(objId,curClass);
		};
	}
	else
	{
		prevClass=LocalGetStringObjectClass(stringId);
		objId=LocalGetNumObjectId(stringId);

		if(curClass==prevClass)
		{
			retCode=ICMLOCAL_CBCALLABLEASW(LocalGetNumObjectId(C_ConvertId),C_CBOptionsIds,C_StartDate,LocalGetNumObjectId(C_IndexId),C_Spread,C_EarlyCall,C_RecallSpread,(bool)hasEuribidClause,(bool)CallIfITM,C_result,objId);
			C_CBOptionsIds.clear();

			if(retCode==ARM_OK)
			{
				LocalSetCurCellEnvValue(curClass,objId);
				stringId=LocalMakeObjectId(objId,curClass);
			};
		}
		else
		{
			FreeCurCellContent();
			retCode=ICMLOCAL_CBCALLABLEASW(LocalGetNumObjectId(C_ConvertId),C_CBOptionsIds,C_StartDate,LocalGetNumObjectId(C_IndexId),C_Spread,C_EarlyCall,C_RecallSpread,(bool)hasEuribidClause,(bool)CallIfITM,C_result);
			C_CBOptionsIds.clear();

			if(retCode==ARM_OK)
			{
				objId=C_result.getLong();
				LocalSetCurCellEnvValue(curClass,objId);
				stringId=LocalMakeObjectId(objId,curClass);
			};
		};
	};

	// return
	if(retCode==ARM_OK)
	{
		FreeCurCellErr();
		XL_result.xltype=xltypeStr;
		XL_result.val.str=XL_StrC2StrPascal(stringId);
		XL_result.xltype|=xlbitDLLFree;
	}
	else
	{
		ARM_ERR();
	};

	return (LPXLOPER)&XL_result;
};
