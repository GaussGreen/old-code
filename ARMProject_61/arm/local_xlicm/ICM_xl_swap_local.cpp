
#pragma warning(disable :4005 4786)

#include <libCCxll\CCxll.h>

#include <ARM\libarm_local\ARM_local_glob.h>
#include <ARM\libarm_local\ARM_local_class.h>

#include <ARM\libicm_local\ICM_local_swap.h>
#include <ARM\libicm_local\ICM_local_glob.h>

#include "ARM_local_interface.h"
#include "ARM_local_interglob.h"

#include "ARM_xl_trycatch_local.h"
#include <util\fromto.h>
#include "ICMKernel/glob/icm_enums.h"

#include "ExcelTools.h"

__declspec(dllexport) LPXLOPER WINAPI ARM_Credit_CDS ( LPXLOPER XL_EffectiveDate,
													   LPXLOPER XL_EndDate,
													   LPXLOPER XL_Rate,
													   LPXLOPER XL_FixingFreq,
													   LPXLOPER XL_DayCountFrq,
													   LPXLOPER XL_First_period_refdate, 	
													   LPXLOPER XL_FixedPayerAmount,	
													   LPXLOPER XL_FloatingPayerAmount,
													   LPXLOPER XL_StubRule, 
													   LPXLOPER XL_Currency,
													   LPXLOPER XL_CreditLag,
													   LPXLOPER XL_Adjusted,
													   LPXLOPER XL_IncludeMaturity,
													   LPXLOPER XL_ProtectionStartDate,
													   LPXLOPER XL_ProtectionEndDate,
													   LPXLOPER XL_Name,
													   LPXLOPER XL_Binary,
													   LPXLOPER XL_First_Cpn_EffDate,
													   LPXLOPER XL_StartAdjusted,
													   LPXLOPER XL_OthersParams	
													   )

{
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	// to make the infrastructure very robust we put a try catch!
	try {
	double EffectiveDate;
	double EndDate;
	double Rate;

	double C_ProtectionStartDate;
	double C_ProtectionEndDate;

	double C_ProtectionStartDate_default = -1.;
	double C_ProtectionEndDate_default = -1.;

	CCString C_FixingFreq;
	int FixingFreq;

	int DayCountFrq;

	double First_period_refdate;
	double First_period_refdate_default = -1.;

	double FixedPayerAmount;	
	double FixedPayerAmount_default = 1000000.;	

	double FloatingPayerAmount;	
	double FloatingPayerAmount_default = -1.;	

	CCString Currency;

	CCString C_Currency;
	CCString C_DayCountFrq;

	CCString C_Adjusted;
	CCString C_StartAdjusted;

	// error
	static int error;
	static char* reason = "";

	CCString C_StubRule;
	CCString C_StubRule_default = "SS";
	int l_Stub;

	double C_CreditLag;
	double C_CreditLag_default = 0.;

	CCString C_IncludeMaturity;
	bool b_IncludeMaturity= false;

	CCString C_Name;
	CCString C_Name_default = "undefine";

	double C_Binary;
	double C_Binary_default = -999.;

	double First_Cpn_EffDate;
	double First_Cpn_EffDate_default = -1.;

	int l_Adjusted = 0;
	int l_StartAdjusted = 1; // pour l'ancier cas K_ADJUSTED ds le ICM_CDS

	vector<string> vOthersParams;

	XL_readNumCell(XL_EffectiveDate,EffectiveDate," ARM_ERR: Effective Date date expected",C_result);
	XL_readNumCell(XL_EndDate,EndDate," ARM_ERR: EndDate date: date expected",C_result);
	XL_readNumCell(XL_Rate,Rate," ARM_ERR: Rate: numeric expected",C_result);

	XL_readStrCellWD(XL_FixingFreq,C_FixingFreq,"Q"," ARM_ERR: FixingFreq: string expected",C_result);
	XL_readStrCellWD(XL_DayCountFrq,C_DayCountFrq,"A360"," ARM_ERR: DayCountFrq: string expected",C_result);
	
	XL_readNumCellWD(XL_First_period_refdate,First_period_refdate,First_period_refdate_default," ARM_ERR: First_period_refdate Rate: date expected",C_result);
	XL_readNumCellWD(XL_FixedPayerAmount,FixedPayerAmount,FixedPayerAmount_default," ARM_ERR: FixedPayerAmount: numeric expected",C_result);
	XL_readNumCellWD(XL_FloatingPayerAmount,FloatingPayerAmount,FloatingPayerAmount_default," ARM_ERR: FloatingPayerAmount: numeric expected",C_result);
	XL_readStrCellWD(XL_StubRule,C_StubRule,C_StubRule_default," ARM_ERR: Stub Rule: string expected",C_result);

	XL_readStrCellWD(XL_Currency,C_Currency,"DEFAULT"," ARM_ERR: currency: string expected",C_result);
	XL_readStrCellWD(XL_Adjusted,C_Adjusted,"ADJ"," ARM_ERR: Adjusted: string expected",C_result);
	
	XL_readNumCellWD(XL_CreditLag,C_CreditLag,C_CreditLag_default," ARM_ERR: Credit Lag: string expected",C_result);

	XL_readNumCellWD(XL_ProtectionStartDate,C_ProtectionStartDate,C_ProtectionStartDate_default," ARM_ERR: Protection start date: numeric expected",C_result);
	XL_readNumCellWD(XL_ProtectionEndDate,C_ProtectionEndDate,C_ProtectionEndDate_default," ARM_ERR: Protection End date: numeric expected",C_result);

	XL_readStrCellWD(XL_IncludeMaturity,C_IncludeMaturity,"Y"," ARM_ERR: Include Maturity: string expected",C_result);
	XL_readStrCellWD(XL_StartAdjusted,C_StartAdjusted,"ADJ"," ARM_ERR: StartAdjusted: string expected",C_result);
	
	XL_readStrCellWD(XL_Name,C_Name,C_Name_default," ARM_ERR: Name: string expected",C_result);

	XL_readNumCellWD(XL_Binary,C_Binary,C_Binary_default," ARM_ERR: Binary: numeric expected",C_result);

	XL_readNumCellWD(XL_First_Cpn_EffDate,First_Cpn_EffDate,First_Cpn_EffDate_default," ARM_ERR: First_Cpn_EffDate : date expected",C_result);

	qPAYMENT_PREMIUM_LEG q_accuredOnDef = qACCRUED_SETTLED; //
	string NotionalEchange;
	long l_NotionalEch_Type= K_NX_NONE; // NotionalEch Type default value
	long l_NotionalEchange = -2; 
	vector<string> vOthersParams_Default(2);
	vOthersParams_Default[0] = "ACC";
	vOthersParams_Default[1] = "NXNONE";
	
	ExcelTools::convert(XL_OthersParams,vOthersParams_Default,vOthersParams);

	if (C_IncludeMaturity == "Y") b_IncludeMaturity = true;

	if (FloatingPayerAmount == FloatingPayerAmount_default)
		FloatingPayerAmount = FixedPayerAmount;

	ARM_result currencyres;
	ARMLOCAL_ARM_GetDefaultCurrency (currencyres);

	if(C_Currency == "DEFAULT")
	{
		if(currencyres.getRetCode () != ARM_OK)
		{
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}
		else
		{
			C_Currency = currencyres.getString ();
		}
	}

	if((l_Stub = ARM_ConvStubRule (C_StubRule)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	if((FixingFreq = ARM_ConvFrequency (C_FixingFreq, C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	if((DayCountFrq = ARM_ConvDayCount (C_DayCountFrq)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	if((l_Adjusted = ARM_ConvIntRule (C_Adjusted)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}
	if((l_StartAdjusted = ARM_ConvStartAdjRule ((const string)C_StartAdjusted.GetStr())) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}
	// case others params :
	/*
		Accrued on def
		Notional echange (ref value)
		notional echange type
	*/
	
	int vSize = vOthersParams.size();
	bool ret = false;
	std::string list;
	
	switch (vSize)
	{
		case 3:
		{
			// Notional echange
			NotionalEchange = vOthersParams[2];
			l_NotionalEchange = LocalGetNumObjectId(NotionalEchange.c_str());
			
		}
		case 2 :
		{
			if((l_NotionalEch_Type = ARM_NotionalExchange(vOthersParams[1].c_str())) == ARM_DEFAULT_ERR)
			{
				ARM_ARG_ERR();
				return (LPXLOPER)&XL_result;
			}
			if(l_NotionalEch_Type != K_NX_NONE && l_NotionalEchange == -1)
			{
				C_result.setMsg("ERROR: Notional exchange must be a reference value in case of this Notional exchange type\n ");
				ARM_ARG_ERR();
				return (LPXLOPER)&XL_result;
			}
		}
		case 1 :
		{
			ICM_EnumsCnv::cnv(vOthersParams[0],q_accuredOnDef);	
		}
		break;
	}

	long prevId = ExcelCaller::get().getObjectId();

	long newId =  ICMLOCAL_CDS(Rate,
							   EffectiveDate,
							   EndDate,
							   First_period_refdate,
							   First_Cpn_EffDate,
							   FixingFreq,
							   DayCountFrq,
							   FixedPayerAmount,
							   FloatingPayerAmount,
							   l_Stub,
							   C_Currency,
							   l_Adjusted,
							   (int) C_CreditLag,
							   b_IncludeMaturity,
							   C_ProtectionStartDate,
							   C_ProtectionEndDate,
							   C_Name,
							   C_Binary,
							   l_StartAdjusted,
								q_accuredOnDef,
								l_NotionalEch_Type,
								l_NotionalEchange,		 
							   C_result);
		

		string objectLabel = ExcelCaller::get().setObject(newId,LOCAL_CDS_CLASS);
		ExcelTools::convert(objectLabel, &XL_result);	
	
	} 
	catch (Exception& e)
	{
		ExcelCaller::get().setError(e.GetErrorString());
		ExcelTools::convert("ARM_ERR", &XL_result);
	}
	catch (std::exception& e)
	{
		ExcelCaller::get().setError(e.what() );
		ExcelTools::convert("ARM_ERR", &XL_result);
	}
	catch (...)
	{
		ExcelCaller::get().setError("unknown error");
		ExcelTools::convert("ARM_ERR", &XL_result);
	}

    return (LPXLOPER)&XL_result;
}


__declspec(dllexport) LPXLOPER WINAPI ARM_Credit_CDSGEN(LPXLOPER XL_FeeLeg,
														LPXLOPER XL_DefLeg,
														LPXLOPER XL_TradedNot)
{
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	try {

	ARM_NOCALCIFWIZ();

	// C variable
	
	CCString C_FeeLeg;
	CCString C_DefLeg;
	double TradedNot = 0.;

	long retCode;

	// error
	static int error;
	static char* reason = "";

	XL_readStrCell(XL_FeeLeg,C_FeeLeg," ARM_ERR: FeeLeg: string expected",C_result);
	XL_readStrCell(XL_DefLeg,C_DefLeg," ARM_ERR: DefLeg: string expected",C_result);
	XL_readNumCell(XL_TradedNot,TradedNot," ARM_ERR: TradedNot: numeric expected",C_result);

	long objId;
	CCString prevClass;
	
	CCString curClass = LOCAL_CDS_CLASS;
	CCString stringId = GetLastCurCellEnvValue ();
	
	if(!stringId)
	{	
		retCode = ICMLOCAL_CDSGEN(LocalGetNumObjectId (C_FeeLeg),
								  LocalGetNumObjectId (C_DefLeg),
								  1., //RcvFeeLeg
								  TradedNot,
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
		retCode = ICMLOCAL_CDSGEN(LocalGetNumObjectId (C_FeeLeg),
								  LocalGetNumObjectId (C_DefLeg),
								  1., //RcvFeeLeg
								  TradedNot,
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

		retCode = ICMLOCAL_CDSGEN(LocalGetNumObjectId (C_FeeLeg),
							   LocalGetNumObjectId (C_DefLeg),
							   1., //RcvFeeLeg
							   TradedNot,
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

//	ARM_END();
	}
	catch (Exception&e)
	{
		ExcelCaller::get().setError(e.GetErrorString()); 
		ExcelTools::convert("ARM_ERR",&XL_result) ;
	}
	catch (std::exception&e)
	{
		ExcelCaller::get().setError(e.what()); 
		ExcelTools::convert("ARM_ERR",&XL_result) ;
	}
	catch(...)
	{
		ExcelCaller::get().setError("Unknown Exception");
		ExcelTools::convert("ARM_ERR",&XL_result) ;
	}

	return (LPXLOPER)&XL_result;
}


__declspec(dllexport) LPXLOPER WINAPI ARM_Credit_FTD ( LPXLOPER XL_EffectiveDate,
													   LPXLOPER XL_EndDate,
													   LPXLOPER XL_Rate,
													   LPXLOPER XL_labels,
													   LPXLOPER XL_FixingFreq,
													   LPXLOPER XL_DayCountFrq,
													   LPXLOPER XL_First_period_refdate, 	
													   LPXLOPER XL_IssuerNotional,	
													   LPXLOPER XL_AccruedOnDefault,
													   LPXLOPER XL_Currency,
													   LPXLOPER XL_CreditLag,
													   LPXLOPER XL_stub,				   
													   LPXLOPER XL_First_Cpn_EffDate,
													   LPXLOPER XL_intRule,
													   LPXLOPER XL_startAdj)

{
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{

	ARM_NOCALCIFWIZ();

	// C variable
	

	double EffectiveDate;
	double EndDate;
	double Rate;

	CCString C_FixingFreq;
	int FixingFreq;

	int DayCountFrq;

	double First_period_refdate;
	double First_period_refdate_default = -1.;

	double IssuerNotional;	
	double IssuerNotional_default = 1000000.;	

	CCString Currency;

	// CCString C_AccruedOnDefault;
	// double AccruedOnDefault = 0.;
	
	CCString C_Currency;
	CCString C_DayCountFrq;

	double C_CreditLag;
	double C_CreditLag_default = 0.;

	//double C_PayCreditLag;
	double C_PayCreditLag_default = 0.;

	// error
	static int error;
	static char* reason = "";

	CCString C_AmortizationOnDefault;
	double AmortizationOnDefault = 0.;

	CCString C_InterestOnDefault;
	double InterestOnDefault = 0.;

	CCString C_Stub;
	int l_Stub ;
	CCString C_intRule;
	int l_intRule =1 ;
	CCString C_startAdj;
	int l_startAdj =1 ;

	VECTOR<CCString> C_Labels;

	double First_Cpn_EffDate;
	double First_Cpn_EffDate_default = -1.;

	XL_readNumCell(XL_EffectiveDate,EffectiveDate," ARM_ERR: Effective Date date expected",C_result);
	XL_readNumCell(XL_EndDate,EndDate," ARM_ERR: EndDate date: date expected",C_result);
	XL_readNumCell(XL_Rate,Rate," ARM_ERR: Rate: numeric expected",C_result);

	XL_readStrVector(XL_labels, C_Labels ," ARM_ERR: Issuers Labels : string vector expected",DOUBLE_TYPE,C_result);
	XL_readStrCellWD(XL_FixingFreq,C_FixingFreq,"Q"," ARM_ERR: FixingFreq: string expected",C_result);
	XL_readStrCellWD(XL_DayCountFrq,C_DayCountFrq,"A360"," ARM_ERR: DayCountFrq: string expected",C_result);
	
	XL_readNumCellWD(XL_First_period_refdate,First_period_refdate,First_period_refdate_default," ARM_ERR: First_period_refdate Rate: date expected",C_result);
	XL_readNumCellWD(XL_IssuerNotional,IssuerNotional,IssuerNotional_default," ARM_ERR: FixedPayerAmount: numeric expected",C_result);
	// XL_readStrCellWD(XL_AccruedOnDefault,C_AccruedOnDefault,"ACC"," ARM_ERR: AccruedOnDefault: string expected",C_result);
	qPAYMENT_PREMIUM_LEG  AccruedOnDefault ;
	ExcelTools::econvert(XL_AccruedOnDefault,"ACC",AccruedOnDefault ); 
	XL_readStrCellWD(XL_Currency,C_Currency,"DEFAULT"," ARM_ERR: currency: string expected",C_result);

	XL_readNumCellWD(XL_CreditLag,C_CreditLag,C_CreditLag_default," ARM_ERR: Credit Lag: string expected",C_result);
	XL_readStrCellWD(XL_stub,C_Stub,"LE"," ARM_ERR: Stub: string expected",C_result);
	XL_readNumCellWD(XL_First_Cpn_EffDate,First_Cpn_EffDate,First_Cpn_EffDate_default," ARM_ERR: First_Cpn_EffDate : date expected",C_result);
	XL_readStrCellWD(XL_intRule,C_intRule,"Y"," ARM_ERR: intRule: string expected",C_result);
	XL_readStrCellWD(XL_startAdj,C_startAdj,"Y"," ARM_ERR: StartAdj : string expected",C_result);
	
	ARM_result currencyres;
	ARMLOCAL_ARM_GetDefaultCurrency (currencyres);

	if(C_Currency == "DEFAULT")
	{
		if(currencyres.getRetCode () != ARM_OK)
		{
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}
		else
		{
			C_Currency = currencyres.getString ();
		}
	}


	long retCode;

	if((FixingFreq = ARM_ConvFrequency (C_FixingFreq, C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	if((DayCountFrq = ARM_ConvDayCount (C_DayCountFrq)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	/**
	if((AccruedOnDefault = ARM_ConvAccOnDef (C_AccruedOnDefault, C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}
	**/ 

/*	if((l_CreditLag = ARM_ConvTypePayLag (C_CreditLag, C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}
*/
	if((l_Stub = ARM_ConvStubRule (C_Stub)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}
	if((l_intRule = ARM_ConvIntRule (C_intRule)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}
	if((l_startAdj = ARM_ConvStartAdjRule ((const string)C_startAdj.GetStr())) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}
	long objId;
	CCString prevClass;
	
	CCString curClass = LOCAL_FTD_CLASS;
	CCString stringId = GetLastCurCellEnvValue ();
	
	if(!stringId)
	{	
		retCode = ICMLOCAL_FTD(Rate,
							   EffectiveDate,
							   EndDate,
							   First_period_refdate,
							   First_Cpn_EffDate,
							   C_Labels,
							   FixingFreq,
							   DayCountFrq,
							   IssuerNotional,
							   AccruedOnDefault,
							   C_Currency,
							   (int)C_CreditLag,
							   l_Stub,	
							   l_intRule,
							   l_startAdj,
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
		retCode = ICMLOCAL_FTD(Rate,
							   EffectiveDate,
							   EndDate,
							   First_period_refdate,
							   First_Cpn_EffDate,
							   C_Labels,
							   FixingFreq,
							   DayCountFrq,
							   IssuerNotional,
							   AccruedOnDefault,
							   C_Currency,
							   (int)C_CreditLag,
							   l_Stub,
							   l_intRule,
							   l_startAdj,
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

		retCode = ICMLOCAL_FTD(Rate,
							   EffectiveDate,
							   EndDate,
							   First_period_refdate,
							   First_Cpn_EffDate,
							   C_Labels,
							   FixingFreq,
							   DayCountFrq,
							   IssuerNotional,
							   AccruedOnDefault,
							   C_Currency,
							   (int)C_CreditLag,
							   l_Stub,
							   l_intRule,
							   l_startAdj,
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

//	ARM_END();
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in ARM_Credit_FTD" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}



__declspec(dllexport) LPXLOPER WINAPI ARM_Credit_NTHTD ( LPXLOPER XL_EffectiveDate,
													   LPXLOPER XL_EndDate,
													   LPXLOPER XL_Rate,
													   LPXLOPER XL_FirstNumDefault,
													   LPXLOPER XL_LastNumDefault,
													   LPXLOPER XL_labels,
													   LPXLOPER XL_FixingFreq,
													   LPXLOPER XL_DayCountFrq,
													   LPXLOPER XL_First_period_refdate, 	
													   LPXLOPER XL_IssuerNotional,	
													   LPXLOPER XL_AccruedOnDefault,
													   LPXLOPER XL_Currency,
													   LPXLOPER XL_CreditLag,
													   LPXLOPER XL_stub,
													   LPXLOPER XL_FrequencyDefLeg,
													   LPXLOPER XL_Binary,
													   LPXLOPER XL_PayCal,
													   LPXLOPER XL_IncludeMaturity,
													   LPXLOPER XL_First_Cpn_EffDate,
													   LPXLOPER XL_AdjParam)

{
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();

	// C variable
	

	double EffectiveDate;
	double EndDate;
	double Rate;

	CCString C_FixingFreq;
	int FixingFreq;
	int DayCountFrq;

	CCString C_FrequencyDefLeg;
	CCString C_FrequencyDefLeg_default = "NONE";

	int FrequencyDefLeg = 0;
	int FrequencyDefLeg_default = -1;

	double First_period_refdate;
	double First_period_refdate_default = -1.;

	double IssuerNotional;	
	double IssuerNotional_default = 1000000.;	

	CCString Currency;

	// CCString C_AccruedOnDefault;
	// double AccruedOnDefault = 0.;

	double FirstNumDefault = 0.;
	double LastNumDefault = 0.;
	double Binary = 0.;
	double Binary_default = -999.;

	CCString C_Currency;
	CCString C_DayCountFrq;

	double C_CreditLag;
	double C_CreditLag_default = 0.;

	//double C_PayCreditLag;
	double C_PayCreditLag_default = 0.;

	// error
	static int error;
	static char* reason = "";

	CCString C_AmortizationOnDefault;
	double AmortizationOnDefault = 0.;

	CCString C_InterestOnDefault;
	double InterestOnDefault = 0.;
	double TradedNotional = 0.;

	CCString C_payCal;

	CCString C_Stub;
	int l_Stub;
	CCString C_intRule;
	int l_intRule =1;
	CCString C_startAdj;
	int l_startAdj =1 ;

	CCString C_IncludeMaturity;
	bool b_IncludeMaturity= false;

	VECTOR<CCString> C_Labels;

	double First_Cpn_EffDate;
	double First_Cpn_EffDate_default = -1.;


	VECTOR<CCString> scheduleAdj;
	scheduleAdj.resize(0);
	CCString C_Adjusted = "Y";
	CCString C_StartAdjusted = "Y";
	int l_Adjusted = 1; // pour les cas par default K_ADJUSTED
	int l_StartAdjusted = 1; // pour les cas par default K_ADJUSTED
	

	char** ppreason = NULL;
	int res = XL_getStrVector (XL_AdjParam, ppreason, XL_ERROR_VALUE_MISSING, XL_ERROR_VALUE_INCONSISTENCY, XL_ERROR_VALUE_IS_ERR, scheduleAdj, NULL, DOUBLE_TYPE) ;
	if ( res == xlerrValue) {		
	}

	if ( scheduleAdj.size() > 2){
		C_result.setMsg("ERROR: there is no more than 2 params for EndAdj-StartAdj ");
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}
	if (scheduleAdj.size() > 0){
		C_Adjusted = scheduleAdj[0];
		if ( scheduleAdj.size() > 1){
			C_StartAdjusted = scheduleAdj[1];
		}
	}

	XL_readNumCell(XL_EffectiveDate,EffectiveDate," ARM_ERR: Effective Date date expected",C_result);
	XL_readNumCell(XL_EndDate,EndDate," ARM_ERR: EndDate date: date expected",C_result);
	XL_readNumCell(XL_Rate,Rate," ARM_ERR: Rate: numeric expected",C_result);
	XL_readNumCell(XL_FirstNumDefault,FirstNumDefault," ARM_ERR: NumDefault: numeric expected",C_result);
	XL_readNumCell(XL_LastNumDefault,LastNumDefault," ARM_ERR: NumDefault: numeric expected",C_result);
	XL_readStrVector(XL_labels, C_Labels ," ARM_ERR: Issuers Labels : string vector expected",DOUBLE_TYPE,C_result);
	XL_readStrCellWD(XL_FixingFreq,C_FixingFreq,"Q"," ARM_ERR: FixingFreq: string expected",C_result);
	XL_readStrCellWD(XL_DayCountFrq,C_DayCountFrq,"A360"," ARM_ERR: DayCountFrq: string expected",C_result);
	XL_readNumCellWD(XL_First_period_refdate,First_period_refdate,First_period_refdate_default," ARM_ERR: First_period_refdate Rate: date expected",C_result);
	XL_readNumCellWD(XL_IssuerNotional,IssuerNotional,IssuerNotional_default," ARM_ERR: FixedPayerAmount: numeric expected",C_result);
	TradedNotional = IssuerNotional;
	// XL_readStrCellWD(XL_AccruedOnDefault,C_AccruedOnDefault,"ACC"," ARM_ERR: AccruedOnDefault: string expected",C_result);
	qPAYMENT_PREMIUM_LEG  AccruedOnDefault ;
	ExcelTools::econvert(XL_AccruedOnDefault,"ACC",AccruedOnDefault ); 

	XL_readStrCellWD(XL_Currency,C_Currency,"DEFAULT"," ARM_ERR: currency: string expected",C_result);
	XL_readNumCellWD(XL_CreditLag,C_CreditLag,C_CreditLag_default," ARM_ERR: Credit Lag: string expected",C_result);
	XL_readStrCellWD(XL_stub,C_Stub,"LE"," ARM_ERR: Stub: string expected",C_result);
	XL_readStrCellWD(XL_DayCountFrq,C_DayCountFrq,"A360"," ARM_ERR: DayCountFrq: string expected",C_result);
	XL_readStrCellWD(XL_FrequencyDefLeg,C_FrequencyDefLeg,C_FrequencyDefLeg_default," ARM_ERR: DefLegFrequency: string expected",C_result);
	XL_readNumCellWD(XL_Binary,Binary,Binary_default," ARM_ERR: Binary: numeric expected",C_result);
	XL_readStrCellWD(XL_PayCal,C_payCal,""," ARM_ERR: pay Calendar: string expected",C_result);
	XL_readStrCellWD(XL_IncludeMaturity,C_IncludeMaturity,"Y"," ARM_ERR: Include Maturity: string expected",C_result);
	XL_readNumCellWD(XL_First_Cpn_EffDate,First_Cpn_EffDate,First_Cpn_EffDate_default," ARM_ERR: First_Cpn_EffDate : date expected",C_result);
	XL_readStrCellWD(XL_IncludeMaturity,C_IncludeMaturity,"Y"," ARM_ERR: Include Maturity: string expected",C_result);
	
	ARM_result currencyres;
	ARMLOCAL_ARM_GetDefaultCurrency (currencyres);

	if (C_IncludeMaturity == "Y") b_IncludeMaturity = true;

	if(C_Currency == "DEFAULT")
	{
		if(currencyres.getRetCode () != ARM_OK)
		{
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}
		else
		{
			C_Currency = currencyres.getString ();
		}
	}


	long retCode;

	if((FixingFreq = ARM_ConvFrequency (C_FixingFreq, C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	if (C_FrequencyDefLeg == "NONE") FrequencyDefLeg = FrequencyDefLeg_default;
	else
	if((FrequencyDefLeg = ARM_ConvFrequency (C_FrequencyDefLeg, C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}


	if((DayCountFrq = ARM_ConvDayCount (C_DayCountFrq)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

//	if((AccruedOnDefault = ARM_ConvAccOnDef (C_AccruedOnDefault, C_result)) == ARM_DEFAULT_ERR)
//	{
//		ARM_ARG_ERR();
//		return (LPXLOPER)&XL_result;
//	}

/*	if((l_CreditLag = ARM_ConvTypePayLag (C_CreditLag, C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}*/

	if((l_Stub = ARM_ConvStubRule (C_Stub)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}
	
	if((l_Adjusted = ARM_ConvIntRule (C_Adjusted)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}
	if((l_StartAdjusted = ARM_ConvStartAdjRule ((const string)C_StartAdjusted.GetStr())) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}
	long objId;
	CCString prevClass;
	
	CCString curClass = LOCAL_NTHTD_CLASS;
	CCString stringId = GetLastCurCellEnvValue ();
	
	if(!stringId)
	{	
		retCode = ICMLOCAL_NTHTD(Rate,
							   EffectiveDate,
							   EndDate,
							   First_period_refdate,
							   First_Cpn_EffDate,
							   (int)FirstNumDefault,
							   (int)LastNumDefault,
							   C_Labels,
							   FixingFreq,
							   DayCountFrq,
							   IssuerNotional,
							   AccruedOnDefault,
							   C_Currency,
							   (int)C_CreditLag,
							   l_Stub,
							   FrequencyDefLeg,
							   Binary,
							   C_payCal,
							   1., //RcvFeeLeg
							   IssuerNotional,
							   b_IncludeMaturity,
							   l_Adjusted,
							   l_StartAdjusted,
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
		retCode = ICMLOCAL_NTHTD(Rate,
							   EffectiveDate,
							   EndDate,
							   First_period_refdate,
							   First_Cpn_EffDate,
							   (int)FirstNumDefault,
							   (int)LastNumDefault,
							   C_Labels,
							   FixingFreq,
							   DayCountFrq,
							   IssuerNotional,
							   AccruedOnDefault,
							   C_Currency,
							   (int)C_CreditLag,
							   l_Stub,
							   FrequencyDefLeg,
							   Binary,
							   C_payCal,
							   1., //RcvFeeLeg
							   IssuerNotional,
							   b_IncludeMaturity,
							   l_Adjusted,
							   l_StartAdjusted,
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

		retCode = ICMLOCAL_NTHTD(Rate,
							   EffectiveDate,
							   EndDate,
							   First_period_refdate,
							   First_Cpn_EffDate,
							   (int)FirstNumDefault,
							   (int)LastNumDefault,
							   C_Labels,
							   FixingFreq,
							   DayCountFrq,
							   IssuerNotional,
							   AccruedOnDefault,
							   C_Currency,
							   (int)C_CreditLag,
							   FrequencyDefLeg,
							   Binary,
							   l_Stub,
							   C_payCal,
							   1., //RcvFeeLeg
							   IssuerNotional,
							   b_IncludeMaturity,
							   l_Adjusted,
							   l_StartAdjusted,
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

//	ARM_END();
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in ARM_Credit_NTHTD" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


__declspec(dllexport) LPXLOPER WINAPI ARM_Credit_Mezzanine ( LPXLOPER XL_EffectiveDate,
													   LPXLOPER XL_EndDate,
													   LPXLOPER XL_Rate,
													   LPXLOPER XL_MezzAmount,
													   LPXLOPER XL_SubAmount,
													   LPXLOPER XL_labels,
													   LPXLOPER XL_notionals,
													   LPXLOPER XL_FreqFeeLeg,
													   LPXLOPER XL_FreqDefLeg,
													   LPXLOPER XL_DayCountFrq,
													   LPXLOPER XL_First_period_refdate, 	
													   LPXLOPER XL_AccruedOnDefault,
													   LPXLOPER XL_Currency,
													   LPXLOPER XL_CreditLag,
													   LPXLOPER XL_stub,
													   LPXLOPER XL_Binary,
													   LPXLOPER XL_PayCal,
													   LPXLOPER XL_TypeFeeLeg,
													   LPXLOPER XL_TypeDefLeg,
													   LPXLOPER XL_SheduleAdj
													   )

{
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();

	// C variable
	

	double EffectiveDate;
	double EndDate;
	double Rate;

	CCString C_FreqFeeLeg;
	int FreqFeeLeg;

	CCString C_FreqDefLeg;
	int FreqDefLeg;

	int DayCountFrq;

	double First_period_refdate;
	double First_period_refdate_default = -1.;

	double FixedPayerAmount_default = 1000000.;	
	double FloatingPayerAmount_default = -1.;	

	double MezzAmount =0.;
	double SubAmount =0.;

	CCString Currency;

	// CCString C_AccruedOnDefault;
	// double AccruedOnDefault = 0.;
	
	CCString C_Currency;
	CCString C_DayCountFrq;
	CCString C_payCal;

	long l_VarSpreadsId =0;

	VECTOR<double> C_Notionals;
	VECTOR<CCString> C_Labels;
	
	VECTOR<CCString> scheduleAdj;
	scheduleAdj.resize(0);
	CCString C_Adjusted = "Y";
	CCString C_StartAdjusted = "Y";
	int l_Adjusted = 1; // pour les cas par default K_ADJUSTED
	int l_StartAdjusted = 1; // pour les cas par default K_ADJUSTED
	CCString C_IncludeMaturity = "Y";
	bool b_IncludeMaturity = true;

	// error
	static int error;
	static char* reason = "";

	CCString C_AmortizationOnDefault;
	double AmortizationOnDefault = 0.;

	CCString C_InterestOnDefault;
	double InterestOnDefault = 0.;

	CCString C_Stub;
	int l_Stub;

	long objId;
	CCString prevClass;
	
	double C_CreditLag;
	double C_CreditLag_default = 0.;

	//double C_PayCreditLag;
	double C_PayCreditLag_default = 0.;

	double C_Binary =0.;
	double C_Binary_default = -999.;

	double TradedNotional = 0.;

	CCString curClass = LOCAL_MEZ_CLASS;
	CCString stringId = GetLastCurCellEnvValue ();
	long retCode;

	CCString C_TypeFeeLeg;
	CCString C_TypeDefLeg;

	long l_TypeFeeLeg;
	long l_TypeDefLeg;

	XL_readNumCell(XL_SubAmount,SubAmount," ARM_ERR: SubAmount: numeric expected",C_result);
	TradedNotional = SubAmount;

    XL_readStrVector(XL_labels, C_Labels ," ARM_ERR: Issuers Labels : string vector expected",DOUBLE_TYPE,C_result);
    XL_readNumVector(XL_notionals, C_Notionals ," ARM_ERR: Notionals : numeric vector expected",C_result);
	
	if (C_Labels.size() != C_Notionals.size())
	{
		C_result.setMsg("ERROR: nb of Issuers labels differ from nb of nominals");
		ARM_ERR();
	}

	char** ppreason = NULL;
	int res = XL_getStrVector (XL_SheduleAdj, ppreason, XL_ERROR_VALUE_MISSING, XL_ERROR_VALUE_INCONSISTENCY, XL_ERROR_VALUE_IS_ERR, scheduleAdj, NULL, DOUBLE_TYPE) ;
	if ( res == xlerrValue) {
		
	}

	if ( scheduleAdj.size() > 3){
		C_result.setMsg("ERROR: there is no more than 3 params for InclMat-EndAdj-StartAdj ");
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}
	if (scheduleAdj.size() > 0){
		C_IncludeMaturity = scheduleAdj[0];
		if ( scheduleAdj.size() > 1){
			C_Adjusted = scheduleAdj[1];
			if ( scheduleAdj.size() > 2){
				C_StartAdjusted = scheduleAdj[2];
			}
		}
	}

	double First_Cpn_EffDate;
	double First_Cpn_EffDate_default = -1.;

	

	XL_readNumCell(XL_EffectiveDate,EffectiveDate," ARM_ERR: Effective Date date expected",C_result);
	XL_readNumCell(XL_EndDate,EndDate," ARM_ERR: EndDate date: date expected",C_result);
	XL_readNumCell(XL_Rate,Rate," ARM_ERR: Rate: numeric expected",C_result);
	XL_readNumCell(XL_MezzAmount,MezzAmount," ARM_ERR: MezzAmount: numeric expected",C_result);
	XL_readStrCellWD(XL_FreqFeeLeg,C_FreqFeeLeg,"Q"," ARM_ERR: FreqFeeLeg: string expected",C_result);
	XL_readStrCellWD(XL_FreqDefLeg,C_FreqDefLeg,"Q"," ARM_ERR: FreqFeeLeg: string expected",C_result);
	XL_readStrCellWD(XL_DayCountFrq,C_DayCountFrq,"A360"," ARM_ERR: DayCountFrq: string expected",C_result);
	XL_readNumCellWD(XL_First_period_refdate,First_period_refdate,First_period_refdate_default," ARM_ERR: First_period_refdate Rate: date expected",C_result);
	// XL_readStrCellWD(XL_AccruedOnDefault,C_AccruedOnDefault,"ACC"," ARM_ERR: AccruedOnDefault: string expected",C_result);
	qPAYMENT_PREMIUM_LEG  AccruedOnDefault ;
	ExcelTools::econvert(XL_AccruedOnDefault,"ACC",AccruedOnDefault ); 

	XL_readStrCellWD(XL_Currency,C_Currency,"DEFAULT"," ARM_ERR: currency: string expected",C_result);
	XL_readNumCellWD(XL_CreditLag,C_CreditLag,C_CreditLag_default," ARM_ERR: Credit Lag: string expected",C_result);
	XL_readStrCellWD(XL_stub,C_Stub,"SS"," ARM_ERR: Stub: string expected",C_result);
	XL_readNumCellWD(XL_Binary,C_Binary,C_Binary_default," ARM_ERR: Binary, numeric expected",C_result);
	XL_readStrCellWD(XL_PayCal,C_payCal,""," ARM_ERR: pay Calendar: string expected",C_result);
	XL_readStrCellWD(XL_TypeFeeLeg,C_TypeFeeLeg,"RUNNING"," ARM_ERR: TypeFeeLeg: string expected",C_result);
	XL_readStrCellWD(XL_TypeDefLeg,C_TypeDefLeg,"RECOVERY"," ARM_ERR: TypeFeeLeg: string expected",C_result);
	//XL_readStrCellWD(XL_IncludeMaturity,C_IncludeMaturity,"Y"," ARM_ERR: Include Maturity: string expected",C_result);
	//XL_readStrCellWD(XL_StartAdjusted,C_StartAdjusted,"ADJ"," ARM_ERR: StartAdjusted: string expected",C_result);
	//XL_readStrCellWD(XL_Adjusted,C_Adjusted,"ADJ"," ARM_ERR: Adjusted: string expected",C_result);
	

//	XL_readNumCellWD(XL_First_Cpn_EffDate,First_Cpn_EffDate,First_Cpn_EffDate_default," ARM_ERR: First_Cpn_EffDate : date expected",C_result);

	First_Cpn_EffDate	=	First_Cpn_EffDate_default;

	ARM_result currencyres;
	ARMLOCAL_ARM_GetDefaultCurrency (currencyres);

	if (C_IncludeMaturity == "Y") {
		b_IncludeMaturity = true;
	} else if ( C_IncludeMaturity == "N"){
		b_IncludeMaturity = false;
	} else {
		C_result.setMsg("ERROR: Y or N expected for IncludeMaturity");
		ARM_ERR();
	}

	if(C_Currency == "DEFAULT")
	{
		if(currencyres.getRetCode () != ARM_OK)
		{
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}
		else
		{
			C_Currency = currencyres.getString ();
		}
	}

	if((FreqFeeLeg = ARM_ConvFrequency (C_FreqFeeLeg, C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	if((FreqDefLeg = ARM_ConvFrequency (C_FreqDefLeg, C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	if((DayCountFrq = ARM_ConvDayCount (C_DayCountFrq)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

/**	if((AccruedOnDefault = ARM_ConvAccOnDef (C_AccruedOnDefault, C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}
**/ 
/*	if((l_CreditLag = ARM_ConvTypePayLag (C_CreditLag, C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}*/

	if((l_Stub = ARM_ConvStubRule (C_Stub)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	if((l_TypeFeeLeg = ARM_ConvLegType(C_TypeFeeLeg, C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	if((l_TypeDefLeg = ARM_ConvLegType(C_TypeDefLeg, C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}
	if((l_Adjusted = ARM_ConvIntRule (C_Adjusted)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}
	if((l_StartAdjusted = ARM_ConvStartAdjRule ((const string)C_StartAdjusted.GetStr())) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	if(!stringId)
	{	
		retCode = ICMLOCAL_MEZZANINE(Rate,
							   EffectiveDate,
							   EndDate,
							   First_period_refdate,
							   First_Cpn_EffDate,
							   MezzAmount,
							   SubAmount,
							   C_Labels,
							   C_Notionals,
							   FreqFeeLeg,
							   FreqDefLeg,
							   DayCountFrq,
							   FixedPayerAmount_default,
							   FloatingPayerAmount_default,
							   AccruedOnDefault,
							   C_Currency,
							   (int)C_CreditLag,
							   l_Stub,
							   C_Binary,
							   C_payCal,
							   1., //RcvFeeLeg
							   TradedNotional,
							   l_TypeFeeLeg,
							   l_TypeDefLeg,
							   b_IncludeMaturity,
							   l_Adjusted,
							   l_StartAdjusted,
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
		retCode = ICMLOCAL_MEZZANINE(Rate,
							   EffectiveDate,
							   EndDate,
							   First_period_refdate,
							   First_Cpn_EffDate,
							   MezzAmount,
							   SubAmount,
							   C_Labels,
							   C_Notionals,
							   FreqFeeLeg,
							   FreqDefLeg,
							   DayCountFrq,
							   FixedPayerAmount_default,
							   FloatingPayerAmount_default,
							   AccruedOnDefault,
							   C_Currency,
							   (int)C_CreditLag,
							   l_Stub,
							   C_Binary,
							   C_payCal,
							   1., //RcvFeeLeg
							   TradedNotional,
							   l_TypeFeeLeg,
							   l_TypeDefLeg,
							   b_IncludeMaturity,
							   l_Adjusted,
							   l_StartAdjusted,
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

		retCode = ICMLOCAL_MEZZANINE(Rate,
							   EffectiveDate,
							   EndDate,
							   First_period_refdate,
							   First_Cpn_EffDate,
							   MezzAmount,
							   SubAmount,
							   C_Labels,
							   C_Notionals,
							   FreqFeeLeg,
							   FreqDefLeg,
							   DayCountFrq,
							   FixedPayerAmount_default,
							   FloatingPayerAmount_default,
							   AccruedOnDefault,
							   C_Currency,
							   (int)C_CreditLag,
							   l_Stub,
							   C_Binary,
							   C_payCal,
							   1., //RcvFeeLeg
							   TradedNotional,
							   l_TypeFeeLeg,
							   l_TypeDefLeg,
							   b_IncludeMaturity,
							   l_Adjusted,
							   l_StartAdjusted,
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

//	ARM_END();
	C_Notionals.clear();
	C_Labels.clear();
}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in ARM_Credit_Mezzanine" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}



__declspec(dllexport) LPXLOPER WINAPI ARM_Credit_CMTranche ( LPXLOPER XL_EffectiveDate,
													   LPXLOPER XL_EndDate,
													   LPXLOPER XL_PartRate,
													   LPXLOPER XL_MezzAmount,
													   LPXLOPER XL_SubAmount,
													   LPXLOPER XL_labels,
													   LPXLOPER XL_notionals,
													   LPXLOPER XL_IndexId,				
													   LPXLOPER XL_FreqFeeLeg,
													   LPXLOPER XL_FreqDefLeg,
													   LPXLOPER XL_DayCountFrq,
													   LPXLOPER XL_First_period_refdate, 	
													   LPXLOPER XL_AccruedOnDefault,
													   LPXLOPER XL_Currency,
													   LPXLOPER XL_PayCreditLag,
													   LPXLOPER XL_stub,
													   LPXLOPER XL_Binary,
													   LPXLOPER XL_PayCal,
													   LPXLOPER XL_FwdFixedDate,
													   LPXLOPER XL_IncludeMaturity)

{
	static XLOPER XL_result;
	ARM_result C_result;

	// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();

	// C variable
	

	double EffectiveDate;
	double EndDate;

	CCString C_FreqFeeLeg;
	int FreqFeeLeg;

	CCString C_FreqDefLeg;
	int FreqDefLeg;

	int DayCountFrq;

	double First_period_refdate;
	double First_period_refdate_default = -1.;

	double FwdFixedDate;
	double FwdFixedDate_default = -1.;

	double FixedPayerAmount_default = 1000000.;	
	double FloatingPayerAmount_default = -1.;	

	double MezzAmount =0.;
	double SubAmount =0.;

	CCString Currency;

	CCString C_IndexId;
	double	 C_ParticipationRate =0.;
	double	 C_ParticipationRate_default =0.5;

	// CCString C_AccruedOnDefault;
	// double AccruedOnDefault = 0.;
	
	CCString C_Currency;
	CCString C_DayCountFrq;
	CCString C_payCal;

	VECTOR<double> C_Notionals;
	VECTOR<CCString> C_Labels;

	// error
	static int error;
	static char* reason = "";

	CCString C_AmortizationOnDefault;
	double AmortizationOnDefault = 0.;

	CCString C_InterestOnDefault;
	double InterestOnDefault = 0.;

	CCString C_Stub;
	int l_Stub;

	long objId;
	CCString prevClass;
	
	double C_CreditLag =0.;

//	int l_CreditLag;

	double C_PayCreditLag;
	double C_PayCreditLag_default = 0.;

	double C_Binary =0.;
	double C_Binary_default = -999.;

	CCString C_IncludeMaturity;
	bool b_IncludeMaturity= false;

	double TradedNotional = 0.;

	CCString curClass = LOCAL_MEZ_CLASS;
	CCString stringId = GetLastCurCellEnvValue ();
	long retCode;

	XL_readNumCell(XL_EffectiveDate,EffectiveDate," ARM_ERR: Effective Date date expected",C_result);
	XL_readNumCell(XL_EndDate,EndDate," ARM_ERR: EndDate date: date expected",C_result);
	XL_readNumCell(XL_PartRate,C_ParticipationRate," ARM_ERR: Participation Rate: numeric expected",C_result);
	XL_readNumCell(XL_MezzAmount,MezzAmount," ARM_ERR: MezzAmount: numeric expected",C_result);
	XL_readNumCell(XL_SubAmount,SubAmount," ARM_ERR: SubAmount: numeric expected",C_result);
	TradedNotional = SubAmount;

    XL_readStrVector(XL_labels, C_Labels ," ARM_ERR: Issuers Labels : string vector expected",DOUBLE_TYPE,C_result);
    XL_readNumVector(XL_notionals, C_Notionals ," ARM_ERR: Notionals : numeric vector expected",C_result);

	if (C_Labels.size() != C_Notionals.size())
	{
		C_result.setMsg("ERROR: nb of Issuers labels differ from nb of nominals");
		ARM_ERR();
	}

	XL_readStrCellWD(XL_IndexId, C_IndexId,"NONE","ARM_ERR: Index : Object Index expected",C_result);
	XL_readStrCellWD(XL_FreqFeeLeg,C_FreqFeeLeg,"Q"," ARM_ERR: FreqFeeLeg: string expected",C_result);
	XL_readStrCellWD(XL_FreqDefLeg,C_FreqDefLeg,"Q"," ARM_ERR: FreqFeeLeg: string expected",C_result);
	XL_readStrCellWD(XL_DayCountFrq,C_DayCountFrq,"A360"," ARM_ERR: DayCountFrq: string expected",C_result);
	XL_readNumCellWD(XL_First_period_refdate,First_period_refdate,First_period_refdate_default," ARM_ERR: First_period_refdate Date: date expected",C_result);
	// XL_readStrCellWD(XL_AccruedOnDefault,C_AccruedOnDefault,"ACC"," ARM_ERR: AccruedOnDefault: string expected",C_result);
	qPAYMENT_PREMIUM_LEG  AccruedOnDefault ;
	ExcelTools::econvert(XL_AccruedOnDefault,"ACC",AccruedOnDefault ); 

	XL_readStrCellWD(XL_Currency,C_Currency,"DEFAULT"," ARM_ERR: currency: string expected",C_result);
	XL_readNumCellWD(XL_PayCreditLag,C_PayCreditLag,C_PayCreditLag_default," ARM_ERR: PayCreditLag: date expected",C_result);
	XL_readStrCellWD(XL_stub,C_Stub,"LE"," ARM_ERR: Stub: string expected",C_result);
	XL_readNumCellWD(XL_Binary,C_Binary,C_Binary_default," ARM_ERR: Binary, numeric expected",C_result);
	XL_readStrCellWD(XL_PayCal,C_payCal,"NULL"," ARM_ERR: pay Calendar: string expected",C_result);
	XL_readNumCellWD(XL_FwdFixedDate,FwdFixedDate,FwdFixedDate_default," ARM_ERR: FwdFixedDate Date: date expected",C_result);
	XL_readStrCellWD(XL_IncludeMaturity,C_IncludeMaturity,"Y"," ARM_ERR: Include Maturity: string expected",C_result);

	ARM_result currencyres;
	ARMLOCAL_ARM_GetDefaultCurrency (currencyres);

	if (C_IncludeMaturity == "Y") b_IncludeMaturity = true;

	if(C_Currency == "DEFAULT")
	{
		if(currencyres.getRetCode () != ARM_OK)
		{
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}
		else
		{
			C_Currency = currencyres.getString ();
		}
	}

	if((FreqFeeLeg = ARM_ConvFrequency (C_FreqFeeLeg, C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	if((FreqDefLeg = ARM_ConvFrequency (C_FreqDefLeg, C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	if((DayCountFrq = ARM_ConvDayCount (C_DayCountFrq)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

/** 	if((AccruedOnDefault = ARM_ConvAccOnDef (C_AccruedOnDefault, C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	} **/ 

/*	if((l_CreditLag = ARM_ConvTypePayLag (C_CreditLag, C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}*/

	if((l_Stub = ARM_ConvStubRule (C_Stub)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	if(!stringId)
	{	
		retCode = ICMLOCAL_CMTranche(
							   EffectiveDate,
							   EndDate,
							   First_period_refdate,
							   MezzAmount,
							   SubAmount,
							   C_Labels,
							   C_Notionals,
							   LocalGetNumObjectId(C_IndexId),	
							   C_ParticipationRate,
							   FreqFeeLeg,
							   FreqDefLeg,
							   DayCountFrq,
							   FixedPayerAmount_default,
							   FloatingPayerAmount_default,
							   AccruedOnDefault,
							   C_Currency,
							   (int)C_CreditLag,
							   l_Stub,
							   C_Binary,
							   C_payCal,
							   1., //RcvFeeLeg
							   TradedNotional,
							   FwdFixedDate,
							   b_IncludeMaturity,
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
		retCode = ICMLOCAL_CMTranche(
							   EffectiveDate,
							   EndDate,
							   First_period_refdate,
							   MezzAmount,
							   SubAmount,
							   C_Labels,
							   C_Notionals,
							   LocalGetNumObjectId(C_IndexId),	
							   C_ParticipationRate,
							   FreqFeeLeg,
							   FreqDefLeg,
							   DayCountFrq,
							   FixedPayerAmount_default,
							   FloatingPayerAmount_default,
							   AccruedOnDefault,
							   C_Currency,
							   (int)C_CreditLag,
							   l_Stub,
							   C_Binary,
							   C_payCal,
							   1., //RcvFeeLeg
							   TradedNotional,
							   FwdFixedDate,
							   b_IncludeMaturity,
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

		retCode = ICMLOCAL_CMTranche(
							   EffectiveDate,
							   EndDate,
							   First_period_refdate,
							   MezzAmount,
							   SubAmount,
							   C_Labels,
							   C_Notionals,
							   LocalGetNumObjectId(C_IndexId),	
							   C_ParticipationRate,
							   FreqFeeLeg,
							   FreqDefLeg,
							   DayCountFrq,
							   FixedPayerAmount_default,
							   FloatingPayerAmount_default,
							   AccruedOnDefault,
							   C_Currency,
							   (int)C_CreditLag,
							   l_Stub,
							   C_Binary,
							   C_payCal,
							   1., //RcvFeeLeg
							   TradedNotional,
							   FwdFixedDate,
							   b_IncludeMaturity,
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

//	ARM_END();
	C_Notionals.clear();
	C_Labels.clear();
}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in ARM_Credit_Mezzanine" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


__declspec(dllexport) LPXLOPER WINAPI ARM_Credit_Cdo2 ( LPXLOPER XL_EffectiveDate,
													   LPXLOPER XL_EndDate,
													   LPXLOPER XL_Portfolio,
													   LPXLOPER XL_Rate,
													   LPXLOPER XL_MezzAmount,
													   LPXLOPER XL_SubAmount,
													   LPXLOPER XL_FreqFeeLeg,
													   LPXLOPER XL_DayCountFrq,
													   LPXLOPER XL_First_period_refdate, 	
													   LPXLOPER XL_AccruedOnDefault,
													   LPXLOPER XL_Currency,
													   LPXLOPER XL_CreditLag,
													   LPXLOPER XL_stub,
													   LPXLOPER XL_FreqDefLeg,
													   LPXLOPER XL_Binary,
													   LPXLOPER XL_PayCal,
													   LPXLOPER XL_CrossSubordination,
													   LPXLOPER XL_IncludeMaturity,
													   LPXLOPER XL_First_Cpn_EffDate)

{
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();

	// C variable
	

	double EffectiveDate;
	double EndDate;
	double Rate;

	CCString C_FreqFeeLeg;
	int FreqFeeLeg;

	CCString C_FreqDefLeg;
	int FreqDefLeg;

	int DayCountFrq;

	double First_period_refdate;
	double First_period_refdate_default = -1.;

	double MezzAmount =0.;
	double SubAmount =0.;

	CCString Currency;

	// CCString C_AccruedOnDefault;
	// double AccruedOnDefault = 0.;
	
	CCString C_Currency;
	CCString C_DayCountFrq;
	CCString C_Ptf;
	CCString C_payCal;
	CCString C_CrossSubordination;
	bool b_CrossSubordination = false;

	// error
	static int error;
	static char* reason = "";

	CCString C_AmortizationOnDefault;
	double AmortizationOnDefault = 0.;

	CCString C_InterestOnDefault;
	double InterestOnDefault = 0.;

	long objId;
	CCString prevClass;
	
	double C_CreditLag=0.;
	double C_CreditLag_default = 0.;

	double C_PayCreditLag_default = 0.;

	double C_Binary =0.;
	double C_Binary_default = -999.;

	double TradedNotional = 0.;

	CCString C_Stub;
	int l_Stub;

	CCString C_IncludeMaturity;
	bool b_IncludeMaturity= false;

	CCString curClass = LOCAL_CDO2_CLASS;
	CCString stringId = GetLastCurCellEnvValue ();
	long retCode;

	double First_Cpn_EffDate;
	double First_Cpn_EffDate_default = -1.;

	XL_readNumCell(XL_SubAmount,SubAmount," ARM_ERR: SubAmount: numeric expected",C_result);
	TradedNotional = SubAmount;

	XL_readNumCell(XL_EffectiveDate,EffectiveDate," ARM_ERR: Effective Date date expected",C_result);
	XL_readNumCell(XL_EndDate,EndDate," ARM_ERR: EndDate date: date expected",C_result);
	XL_readNumCell(XL_Rate,Rate," ARM_ERR: Rate: numeric expected",C_result);
	XL_readNumCell(XL_MezzAmount,MezzAmount," ARM_ERR: MezzAmount: numeric expected",C_result);
	XL_readStrCell(XL_Portfolio,C_Ptf," ARM_ERR: Portfolio id expected",C_result);
	XL_readStrCellWD(XL_FreqFeeLeg,C_FreqFeeLeg,"Q"," ARM_ERR: FreqFeeLeg : string expected",C_result);
	XL_readStrCellWD(XL_FreqDefLeg,C_FreqDefLeg,"Q"," ARM_ERR: FreqDefLeg : string expected",C_result);
	XL_readStrCellWD(XL_DayCountFrq,C_DayCountFrq,"A360"," ARM_ERR: DayCountFrq: string expected",C_result);
	XL_readNumCellWD(XL_First_period_refdate,First_period_refdate,First_period_refdate_default," ARM_ERR: First_period_refdate Rate: date expected",C_result);
	// XL_readStrCellWD(XL_AccruedOnDefault,C_AccruedOnDefault,"ACC"," ARM_ERR: AccruedOnDefault: string expected",C_result);
	qPAYMENT_PREMIUM_LEG  AccruedOnDefault ;
	ExcelTools::econvert(XL_AccruedOnDefault,"ACC",AccruedOnDefault ); 

	XL_readStrCellWD(XL_Currency,C_Currency,"DEFAULT"," ARM_ERR: currency: string expected",C_result);
	XL_readNumCellWD(XL_CreditLag,C_CreditLag,C_CreditLag_default," ARM_ERR: Credit Lag: string expected",C_result);
	XL_readStrCellWD(XL_stub,C_Stub,"LE"," ARM_ERR: Stub: string expected",C_result);
	XL_readNumCellWD(XL_Binary,C_Binary,C_Binary_default," ARM_ERR: Binary, numeric expected",C_result);
	XL_readStrCellWD(XL_PayCal,C_payCal,"NULL"," ARM_ERR: pay Calendar: string expected",C_result);
	XL_readStrCellWD(XL_CrossSubordination,C_CrossSubordination,"N"," ARM_ERR: CrossSub: string expected",C_result);
	XL_readStrCellWD(XL_IncludeMaturity,C_IncludeMaturity,"Y"," ARM_ERR: Include Maturity: string expected",C_result);

	XL_readNumCellWD(XL_First_Cpn_EffDate,First_Cpn_EffDate,First_Cpn_EffDate_default," ARM_ERR: First_Cpn_EffDate : date expected",C_result);

	ARM_result currencyres;
	ARMLOCAL_ARM_GetDefaultCurrency (currencyres);

	if (C_CrossSubordination=="Y") b_CrossSubordination=true;
	if (C_IncludeMaturity == "Y") b_IncludeMaturity = true;

	if(C_Currency == "DEFAULT")
	{
		if(currencyres.getRetCode () != ARM_OK)
		{
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}
		else
		{
			C_Currency = currencyres.getString ();
		}
	}

	if((FreqFeeLeg = ARM_ConvFrequency (C_FreqFeeLeg, C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	if((FreqDefLeg = ARM_ConvFrequency (C_FreqDefLeg, C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	if((DayCountFrq = ARM_ConvDayCount (C_DayCountFrq)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

/**	if((AccruedOnDefault = ARM_ConvAccOnDef (C_AccruedOnDefault, C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}
	**/ 

/*	if((l_CreditLag = ARM_ConvTypePayLag (C_CreditLag, C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}*/

	if((l_Stub = ARM_ConvStubRule (C_Stub)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	if(!stringId)
	{	
		retCode = ICMLOCAL_CDO2(Rate,
							   EffectiveDate,
							   EndDate,
							   First_period_refdate,
							   First_Cpn_EffDate,
							   MezzAmount,
							   SubAmount,
							   FreqFeeLeg,
							   DayCountFrq,
							   AccruedOnDefault,
							   C_Currency,
							   LocalGetNumObjectId(C_Ptf),
							   (int)C_CreditLag,
							   l_Stub,
							   FreqDefLeg,
							   C_Binary,
							   C_payCal,
							   1., //RcvFeeLeg
							   SubAmount,
							   b_CrossSubordination,
							   b_IncludeMaturity,
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
		retCode = ICMLOCAL_CDO2(Rate,
							   EffectiveDate,
							   EndDate,
							   First_period_refdate,
							   First_Cpn_EffDate,
							   MezzAmount,
							   SubAmount,
							   FreqFeeLeg,
							   DayCountFrq,
							   AccruedOnDefault,
							   C_Currency,
							   LocalGetNumObjectId(C_Ptf),
							   (int)C_CreditLag,
							   l_Stub,
							   FreqDefLeg,
							   C_Binary,
							   C_payCal,
							   1., //RcvFeeLeg
							   SubAmount,
							   b_CrossSubordination,
							   b_IncludeMaturity,
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

		retCode = ICMLOCAL_CDO2(Rate,
							   EffectiveDate,
							   EndDate,
							   First_period_refdate,
							   First_Cpn_EffDate,
							   MezzAmount,
							   SubAmount,
							   FreqFeeLeg,
							   DayCountFrq,
							   AccruedOnDefault,
							   C_Currency,
							   LocalGetNumObjectId(C_Ptf),
							   (int)C_CreditLag,
							   l_Stub,
							   FreqDefLeg,
							   C_Binary,
							   C_payCal,
							   1., //RcvFeeLeg
							   SubAmount,
							   b_CrossSubordination,
							   b_IncludeMaturity,
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

//	ARM_END();
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in ARM_Credit_Cdo2" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}

__declspec(dllexport) LPXLOPER WINAPI ARM_Credit_CDSIndex (LPXLOPER XL_EffectiveDate,
														   LPXLOPER XL_EndDate,
														   LPXLOPER XL_Rate,
														   LPXLOPER XL_IndexId,
														   LPXLOPER XL_FixingFreq,
														   LPXLOPER XL_DayCountFrq,
														   LPXLOPER XL_First_period_refdate, 	
														   LPXLOPER XL_FixedPayerAmount,	
														   LPXLOPER XL_FloatingPayerAmount,
														   LPXLOPER XL_StubRule, 
														   LPXLOPER XL_Currency,
														   LPXLOPER XL_CreditLag,
														   LPXLOPER XL_Adjusted,
														   LPXLOPER XL_IncludeMaturity,
														   LPXLOPER XL_ProtectionStartDate,
														   LPXLOPER XL_ProtectionEndDate)

{
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();

	// C variable
	

	double EffectiveDate;
	double EndDate;
	double Rate;

	double C_ProtectionStartDate;
	double C_ProtectionEndDate;

	double C_ProtectionStartDate_default = -1.;
	double C_ProtectionEndDate_default = -1.;

	CCString C_IndexId;

	CCString C_FixingFreq;
	int FixingFreq;

	int DayCountFrq;

	double First_period_refdate;
	double First_period_refdate_default = -1.;

	double FixedPayerAmount;	
	double FixedPayerAmount_default = 1000000.;	

	double FloatingPayerAmount;	
	double FloatingPayerAmount_default = -1.;	

	CCString Currency;
	CCString C_Currency;
	CCString C_DayCountFrq;
	CCString C_Adjusted;

	// error
	static int error;
	static char* reason = "";

	CCString C_StubRule;
	CCString C_StubRule_default = "SS";
	int l_Stub;

	double C_CreditLag;
	double C_CreditLag_default = 0.;

	CCString C_IncludeMaturity;
	bool b_IncludeMaturity= false;

	int l_Adjusted = 0;

	XL_readStrCellWD(XL_IndexId, C_IndexId,"NONE","ARM_ERR: Index : Object Index expected",C_result);
	XL_readNumCell(XL_EffectiveDate,EffectiveDate," ARM_ERR: Effective Date date expected",C_result);
	XL_readNumCell(XL_EndDate,EndDate," ARM_ERR: EndDate date: date expected",C_result);
	XL_readNumCell(XL_Rate,Rate," ARM_ERR: Rate: numeric expected",C_result);

	XL_readStrCellWD(XL_FixingFreq,C_FixingFreq,"Q"," ARM_ERR: FixingFreq: string expected",C_result);
	XL_readStrCellWD(XL_DayCountFrq,C_DayCountFrq,"A360"," ARM_ERR: DayCountFrq: string expected",C_result);
	
	XL_readNumCellWD(XL_First_period_refdate,First_period_refdate,First_period_refdate_default," ARM_ERR: First_period_refdate Rate: date expected",C_result);
	XL_readNumCellWD(XL_FixedPayerAmount,FixedPayerAmount,FixedPayerAmount_default," ARM_ERR: FixedPayerAmount: numeric expected",C_result);
	XL_readNumCellWD(XL_FloatingPayerAmount,FloatingPayerAmount,FloatingPayerAmount_default," ARM_ERR: FloatingPayerAmount: numeric expected",C_result);
	XL_readStrCellWD(XL_StubRule,C_StubRule,C_StubRule_default," ARM_ERR: Stub Rule: string expected",C_result);

	XL_readStrCellWD(XL_Currency,C_Currency,"DEFAULT"," ARM_ERR: currency: string expected",C_result);
	XL_readStrCellWD(XL_Adjusted,C_Adjusted,"ADJ"," ARM_ERR: Adjusted: string expected",C_result);
	XL_readNumCellWD(XL_CreditLag,C_CreditLag,C_CreditLag_default," ARM_ERR: Credit Lag: string expected",C_result);

	XL_readNumCellWD(XL_ProtectionStartDate,C_ProtectionStartDate,C_ProtectionStartDate_default," ARM_ERR: Protection start date: numeric expected",C_result);
	XL_readNumCellWD(XL_ProtectionEndDate,C_ProtectionEndDate,C_ProtectionEndDate_default," ARM_ERR: Protection End date: numeric expected",C_result);

	XL_readStrCellWD(XL_IncludeMaturity,C_IncludeMaturity,"Y"," ARM_ERR: Include Maturity: string expected",C_result);

	if (C_IncludeMaturity == "Y") b_IncludeMaturity = true;

	if (FloatingPayerAmount == FloatingPayerAmount_default)
		FloatingPayerAmount = FixedPayerAmount;

	ARM_result currencyres;
	ARMLOCAL_ARM_GetDefaultCurrency (currencyres);

	if(C_Currency == "DEFAULT")
	{
		if(currencyres.getRetCode () != ARM_OK)
		{
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}
		else
		{
			C_Currency = currencyres.getString ();
		}
	}

	long retCode;

	if((l_Stub = ARM_ConvStubRule (C_StubRule)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	if((FixingFreq = ARM_ConvFrequency (C_FixingFreq, C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	if((DayCountFrq = ARM_ConvDayCount (C_DayCountFrq)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	if((l_Adjusted = ARM_ConvIntRule (C_Adjusted)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	long objId;
	CCString prevClass;
	
	CCString curClass = LOCAL_CDSINDEX_CLASS;
	CCString stringId = GetLastCurCellEnvValue ();
	
	if(!stringId)
	{	
		retCode = ICMLOCAL_CDSIndex(Rate,
									LocalGetNumObjectId(C_IndexId),	
									EffectiveDate,
									EndDate,
									First_period_refdate,
									FixingFreq,
									DayCountFrq,
									FixedPayerAmount,
									FloatingPayerAmount,
									l_Stub,
									C_Currency,
									l_Adjusted,
									(int) C_CreditLag,
									b_IncludeMaturity,
									C_ProtectionStartDate,
									C_ProtectionEndDate,
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
			retCode = ICMLOCAL_CDSIndex(Rate,
										LocalGetNumObjectId(C_IndexId),
										EffectiveDate,
										EndDate,
										First_period_refdate,
										FixingFreq,
										DayCountFrq,
										FixedPayerAmount,
										FloatingPayerAmount,
										l_Stub,
										C_Currency,
										l_Adjusted,
										(int) C_CreditLag,
										b_IncludeMaturity,
										C_ProtectionStartDate,
										C_ProtectionEndDate,
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

			retCode = ICMLOCAL_CDSIndex(Rate,
										LocalGetNumObjectId(C_IndexId),
										EffectiveDate,
										EndDate,
										First_period_refdate,
										FixingFreq,
										DayCountFrq,
										FixedPayerAmount,
										FloatingPayerAmount,
										l_Stub,
										C_Currency,
										l_Adjusted,
										(int) C_CreditLag,
										b_IncludeMaturity,
										C_ProtectionStartDate,
										C_ProtectionEndDate,
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

	//	ARM_END();
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END
		
	/// to catch arm exception
	
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in ARM_Credit_CDS" )
	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


__declspec(dllexport) LPXLOPER WINAPI ARM_Credit_Option(LPXLOPER XL_UnderlyingMaturity,
														LPXLOPER XL_ExpiryDate,
														LPXLOPER XL_CCy,
														LPXLOPER XL_cdsAdj,
														LPXLOPER XL_endAdj,
														LPXLOPER XL_Strike,
														LPXLOPER XL_OptionType,
														LPXLOPER XL_KoType,
														LPXLOPER XL_Quantity,
														LPXLOPER XL_UnderlyingType)
{
// return
    static XLOPER XL_result;
	ARM_result C_result;
	try
	{
		ARM_NOCALCIFWIZ();
		long prevId = ExcelCaller::get().getObjectId();
	
		string UnderMaturity = "";
		string ccy="";
		string cdsAdj="";
		qCDS_ADJ qCdsAdj=(qCDS_ADJ)0;
		string EndAdj="";
		bool boolEndAdj=false;
		double expiryDate = -1;
		double Strike = 0.;
		double defaultStrike = 0.;
		string C_OptionType ="";
		string C_OptionType_default = "CALL";
    
		string C_KoType="";
		string C_KoType_default = "KO";

		long l_OptionType = 0.;
		long l_KoType = 0. ;
		
		long retCode = 0;

		double notional = 1;
		const double notional_default = 1; 
		double xlDate =0.;
		ARM_Date MaturityDate ;
		try {
			 ExcelTools::convert(XL_UnderlyingMaturity,xlDate); // try double conversion for date
			 ExcelTools::convert(XL_UnderlyingMaturity,MaturityDate);
			 UnderMaturity = "";
		}catch (...){
			ExcelTools::convert(XL_UnderlyingMaturity,UnderMaturity);
		}

		//ExcelTools::convert( XL_UnderlyingMaturity,UnderMaturity);	
		ExcelTools::convert( XL_ExpiryDate,expiryDate);	
		ExcelTools::convert(XL_CCy,ccy);
		ExcelTools::convert(XL_CCy,ccy);
		ExcelTools::convert(XL_cdsAdj,"STDCDS",cdsAdj);
		ExcelTools::convert(XL_endAdj,"Y",EndAdj);
		ExcelTools::convert(XL_Strike,defaultStrike, Strike);
		ExcelTools::convert(XL_OptionType,"CALL",C_OptionType);
		ExcelTools::convert(XL_KoType,"KO",C_KoType);
		ExcelTools::convert(XL_Quantity,notional_default, notional);

		if(EndAdj=="Y") boolEndAdj = true;
		if((l_OptionType = ARM_ConvCallOrPut(CCString(C_OptionType.c_str()), C_result)) == ARM_DEFAULT_ERR){	ARM_ARG_ERR();	return (LPXLOPER)&XL_result;}

		if((l_KoType = ARM_ConvKO_Or_NoKO(CCString(C_KoType.c_str()), C_result)) == ARM_DEFAULT_ERR){	ARM_ARG_ERR();	return (LPXLOPER)&XL_result;}

		qCdsAdj = ARM_ConvAdjCalCDS(CCString(cdsAdj.c_str()),C_result);

		long newId  = ICMLOCAL_Option(UnderMaturity,MaturityDate,expiryDate, ccy, qCdsAdj,boolEndAdj,
										Strike,l_OptionType,(qDEF_MAT)l_KoType, notional, C_result);
		string objectLabel = ExcelCaller::get().setObject(newId, LOCAL_OPTION_CLASS);
		//ExcelCaller::get().setObject(objectLabel);
		ExcelTools::convert(objectLabel, &XL_result);	
	} 
	catch (Exception& e)
	{
		ExcelCaller::get().setError(e.GetErrorString());
		ExcelTools::convert("ARM_ERR", &XL_result);
	}
	catch (std::exception& e)
	{
		ExcelCaller::get().setError(e.what() );
		ExcelTools::convert("ARM_ERR", &XL_result);
	}
	catch (...)
	{
		ExcelCaller::get().setError("unknown error");
		ExcelTools::convert("ARM_ERR", &XL_result);
	}
    return (LPXLOPER)&XL_result;

}
//	-------------------------------------------------------------------------------------
//	JLA. 
__declspec(dllexport) LPXLOPER WINAPI ARM_Credit_SpreadOption(LPXLOPER XL_UnderlyingInstrument,													   
														LPXLOPER XL_Strike,
														LPXLOPER XL_OptionType,
														LPXLOPER XL_KoType,
														LPXLOPER XL_AccType,
														LPXLOPER XL_ExerciseDates,
														LPXLOPER XL_ExerciseFrequency,
														LPXLOPER XL_ExerciseStyle,
														LPXLOPER XL_UndMatuStyle)
{
	static XLOPER XL_result;
	ARM_result C_result;

	try {

	ARM_NOCALCIFWIZ();

	double Strike;
	CCString C_OptionType;
	// CCString C_KoType;
	CCString C_KoStyle ;
	CCString C_AccStyle;
	CCString C_UnderlyingInst; 
	CCString C_ExerciseStyle ;
	CCString C_UnderlyingMatuStyle ; 
	std::vector<double> C_ExerciseDates;  
	long l_OptionType;
	long koStyle;
	long accStyle;
	long exerciseStyle; 
	double exerciseFrequency; 
	long  matuStyle ;
	long retCode = 0;
	long objId;
	
	//ARM_Vector exerciseDates; 

	// error
	static int error;
	static char* reason = "";

	XL_readStrCell(XL_UnderlyingInstrument,C_UnderlyingInst," ARM_ERR: Security Id: Object expected",C_result);
	XL_readNumCell(XL_Strike,Strike," ARM_ERR: numeric expected",C_result);
	XL_readStrCellWD(XL_OptionType,C_OptionType,"CALL"," ARM_ERR: Option Type : CALL ou PUT expected",C_result);
	// XL_readStrCellWD(XL_KoType,C_KoType,"KO"," ARM_ERR: KO, No_KO_ACC or No_KO_No_ACC expected",C_result);
	XL_readStrCellWD(XL_KoType,C_KoStyle,"KO"," ARM_ERR: KO, NKO expected",C_result);
	XL_readStrCellWD(XL_AccType,C_AccStyle,"NACC"," ARM_ERR: ACC, NACC expected",C_result);
	XL_readStrCellWD(XL_ExerciseStyle,C_ExerciseStyle,"EUROPEAN"," ARM_ERR: ExerciseStyle ?",C_result);
    XL_readNumVector(XL_ExerciseDates, C_ExerciseDates," ARM_ERR: ExerciseDates : numeric vector expected",C_result);
	XL_readNumCell(XL_ExerciseFrequency,exerciseFrequency," ARM_ERR: numeric expected",C_result);
	XL_readStrCellWD(XL_UndMatuStyle,C_UnderlyingMatuStyle ,"CONSTANT"," ARM_ERR: CONSTANT,RESIDUAL expected",C_result);

	if((l_OptionType = ARM_ConvCallOrPut(C_OptionType, C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}
	
	// if((l_KoType = ARM_ConvKO_Or_NoKO(C_KoType, C_result)) == ARM_DEFAULT_ERR)
	// {
	// 	ARM_ARG_ERR();
	// 	return (LPXLOPER)&XL_result;
	// }
	if((exerciseStyle = ARM_ConvExerciseType(C_ExerciseStyle, C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}
	if((matuStyle= ARM_ConvUnderlyingMatuStyle(C_UnderlyingMatuStyle, C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}
	if((koStyle= ARM_ConvKOStyle(C_KoStyle, C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}
	if((accStyle= ARM_ConvAccStyle(C_AccStyle, C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}
	
	
	CCString prevClass;
	
	CCString curClass = LOCAL_OPTION_CLASS;
	CCString stringId = GetLastCurCellEnvValue ();

	bool isCall = (l_OptionType==K_CALL); 
	// bool isKO = ( l_KoType==K_KO ) ;
	// bool isAccelerated = ( l_KoType==K_No_KO_ACC ) ; 

	if (!isCall)
	{
		koStyle = 0 ; // Always KO with puts
		accStyle = 1 ;
	}

	if(!stringId)
	{
		retCode = ICMLOCAL_SpreadOption( LocalGetNumObjectId (C_UnderlyingInst),
							 Strike,
							 isCall,
							 (int)/*(qKoStyle)*/koStyle,
							 (int)/*(qAccelerationStyle)*/accStyle,
							 C_ExerciseDates,
							 exerciseStyle,
							 (int) exerciseFrequency,
							 (int)/*(qUnderlying_Maturity_Style)*/matuStyle,
							 C_result,
							 -1) ;


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
					retCode = ICMLOCAL_SpreadOption( LocalGetNumObjectId (C_UnderlyingInst),
							 Strike,
							 isCall,
							 (int)/*(qKoStyle)*/koStyle,
							 (int)/*(qAccelerationStyle)*/accStyle,
							 C_ExerciseDates,
							 exerciseStyle,
							 (int)exerciseFrequency,
							 (int)/*(qUnderlying_Maturity_Style)*/matuStyle,
							 C_result,
							 objId) ;


/*			retCode = ICMLOCAL_Option(LocalGetNumObjectId (C_UnderlyingPricer),
									  C_Strike,
									  l_OptionType,
									  l_KoType,
									  C_result,
									  objId); */ 
			
			if(retCode == ARM_OK)
			{			
				LocalSetCurCellEnvValue (curClass, objId); 

				stringId = LocalMakeObjectId (objId, curClass);
			}
		}
		else
		{
			FreeCurCellContent ();

					retCode = ICMLOCAL_SpreadOption( LocalGetNumObjectId (C_UnderlyingInst),
							 Strike,
							 isCall,
							 (int)/*(qKoStyle)*/koStyle,
							 (int)/*(qAccelerationStyle)*/accStyle,
							 C_ExerciseDates,
							 exerciseStyle,
							 (int)exerciseFrequency,
							 (int)/*(qUnderlying_Maturity_Style)*/matuStyle,
							 C_result,
							 objId) ;


					/*retCode = ICMLOCAL_Option(LocalGetNumObjectId (C_UnderlyingPricer),
									  C_Strike,
									  l_OptionType,
									  l_KoType,
									  C_result);
									  */ 
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
	catch (Exception&e)
	{
		ExcelCaller::get().setError(e.GetErrorString()); 
		ExcelTools::convert("ARM_ERR",&XL_result) ;
	}
	catch (std::exception&e)
	{
		ExcelCaller::get().setError(e.what()); 
		ExcelTools::convert("ARM_ERR",&XL_result) ;
	}
	catch(...)
	{
		ExcelCaller::get().setError("Unknown Exception");
		ExcelTools::convert("ARM_ERR",&XL_result) ;
	}

	return (LPXLOPER)&XL_result;

}

// -------------------------------------------------------------------------------
// CMCDS
// -------------------------------------------------------------------------------
__declspec(dllexport) LPXLOPER WINAPI ARM_Credit_CMCDS(LPXLOPER XL_EffectiveDate,
													   LPXLOPER XL_EndDate,
													   LPXLOPER XL_PartRate,
													   LPXLOPER XL_IndexId,
													   LPXLOPER XL_FixingFreq,
													   LPXLOPER XL_DayCountFrq,
													   LPXLOPER XL_First_period_refdate, 	
													   LPXLOPER XL_FixedPayerAmount,	
													   LPXLOPER XL_FloatingPayerAmount,
													   LPXLOPER XL_StubRule,
													   LPXLOPER XL_Currency,
													   LPXLOPER XL_CreditLag,
													   LPXLOPER XL_Adjusted,
													   LPXLOPER XL_IncludeMaturity,
													   LPXLOPER XL_ProtectionStartDate,
													   LPXLOPER XL_ProtectionEndDate,
													   LPXLOPER XL_First_Cpn_EffDate)

{
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();

	// C variable
	

	double EffectiveDate;
	double EndDate;
	double PartRate;

	double C_ProtectionStartDate;
	double C_ProtectionEndDate;

	double C_ProtectionStartDate_default = -1.;
	double C_ProtectionEndDate_default = -1.;

	CCString C_FixingFreq;
	int FixingFreq;

	int DayCountFrq;

	double First_period_refdate;
	double First_period_refdate_default = -1.;

	double FixedPayerAmount;	
	double FixedPayerAmount_default = 1000000.;	

	double FloatingPayerAmount;	
	double FloatingPayerAmount_default = -1.;	

	CCString Currency;

	CCString C_Currency;
	CCString C_DayCountFrq;

	CCString C_Adjusted;

	CCString C_IndexId;

	// error
	static int error;
	static char* reason = "";

	CCString C_StubRule;
	CCString C_StubRule_default = "SS";
	int l_Stub;

	double C_CreditLag;
	double C_CreditLag_default = 0.;

	CCString C_IncludeMaturity;
	bool b_IncludeMaturity= false;

	int l_Adjusted = 0;

	double First_Cpn_EffDate;
	double First_Cpn_EffDate_default = -1.;

	XL_readNumCell(XL_EffectiveDate,EffectiveDate," ARM_ERR: Effective Date date expected",C_result);
	XL_readNumCell(XL_EndDate,EndDate," ARM_ERR: EndDate date: date expected",C_result);
	XL_readNumCell(XL_PartRate,PartRate," ARM_ERR: PartRate: numeric expected",C_result);
	XL_readStrCellWD(XL_IndexId, C_IndexId,"NONE","ARM_ERR: Index : Object Index expected",C_result);

	XL_readStrCellWD(XL_FixingFreq,C_FixingFreq,"Q"," ARM_ERR: FixingFreq: string expected",C_result);
	XL_readStrCellWD(XL_DayCountFrq,C_DayCountFrq,"A360"," ARM_ERR: DayCountFrq: string expected",C_result);
	
	XL_readNumCellWD(XL_First_period_refdate,First_period_refdate,First_period_refdate_default," ARM_ERR: First_period_refdate Rate: date expected",C_result);
	XL_readNumCellWD(XL_FixedPayerAmount,FixedPayerAmount,FixedPayerAmount_default," ARM_ERR: FixedPayerAmount: numeric expected",C_result);
	XL_readNumCellWD(XL_FloatingPayerAmount,FloatingPayerAmount,FloatingPayerAmount_default," ARM_ERR: FloatingPayerAmount: numeric expected",C_result);
	XL_readStrCellWD(XL_StubRule,C_StubRule,C_StubRule_default," ARM_ERR: Stub Rule: string expected",C_result);
	
	XL_readStrCellWD(XL_Currency,C_Currency,"DEFAULT"," ARM_ERR: currency: string expected",C_result);
	XL_readStrCellWD(XL_Adjusted,C_Adjusted,"ADJ"," ARM_ERR: Adjusted: string expected",C_result);
	XL_readNumCellWD(XL_CreditLag,C_CreditLag,C_CreditLag_default," ARM_ERR: Credit Lag: string expected",C_result);

	XL_readNumCellWD(XL_ProtectionStartDate,C_ProtectionStartDate,C_ProtectionStartDate_default," ARM_ERR: Protection start date: numeric expected",C_result);
	XL_readNumCellWD(XL_ProtectionEndDate,C_ProtectionEndDate,C_ProtectionEndDate_default," ARM_ERR: Protection End date: numeric expected",C_result);

	XL_readStrCellWD(XL_IncludeMaturity,C_IncludeMaturity,"Y"," ARM_ERR: Include Maturity: string expected",C_result);

	XL_readNumCellWD(XL_First_Cpn_EffDate,First_Cpn_EffDate,First_Cpn_EffDate_default," ARM_ERR: First_Cpn_EffDate : date expected",C_result);
	
	if (C_IncludeMaturity == "Y") b_IncludeMaturity = true;

	if (FloatingPayerAmount == FloatingPayerAmount_default)
		FloatingPayerAmount = FixedPayerAmount;

	ARM_result currencyres;
	ARMLOCAL_ARM_GetDefaultCurrency (currencyres);

	if(C_Currency == "DEFAULT")
	{
		if(currencyres.getRetCode () != ARM_OK)
		{
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}
		else
		{
			C_Currency = currencyres.getString ();
		}
	}

	long retCode;

	if((l_Stub = ARM_ConvStubRule (C_StubRule)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	if((FixingFreq = ARM_ConvFrequency (C_FixingFreq, C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	if((DayCountFrq = ARM_ConvDayCount (C_DayCountFrq)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	if((l_Adjusted = ARM_ConvIntRule (C_Adjusted)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	long objId;
	CCString prevClass;
	
	CCString curClass = LOCAL_CDS_CLASS;
	CCString stringId = GetLastCurCellEnvValue ();
	
	if(!stringId)
	{	
		retCode = ICMLOCAL_CMCDS ( PartRate,
								   EffectiveDate,
								   EndDate,
								   First_Cpn_EffDate,
								   LocalGetNumObjectId(C_IndexId),
								   First_period_refdate,
								   FixingFreq,
								   DayCountFrq,
								   FixedPayerAmount,
								   FloatingPayerAmount,
								   l_Stub,
								   C_Currency,
								   l_Adjusted,
								   (int) C_CreditLag,
								   b_IncludeMaturity,
								   C_ProtectionStartDate,
								   C_ProtectionEndDate,
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
		retCode = ICMLOCAL_CMCDS ( PartRate,
								   EffectiveDate,
								   EndDate,
								   First_Cpn_EffDate,
								   LocalGetNumObjectId(C_IndexId),
								   First_period_refdate,
								   FixingFreq,
								   DayCountFrq,
								   FixedPayerAmount,
								   FloatingPayerAmount,
								   l_Stub,
								   C_Currency,
								   l_Adjusted,
								   (int) C_CreditLag,
								   b_IncludeMaturity,
								   C_ProtectionStartDate,
								   C_ProtectionEndDate,
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

		retCode = ICMLOCAL_CMCDS(  PartRate,
							       EffectiveDate,
  							       EndDate,
								   First_Cpn_EffDate,
								   LocalGetNumObjectId(C_IndexId),
								   First_period_refdate,
								   FixingFreq,
								   DayCountFrq,
								   FixedPayerAmount,
								   FloatingPayerAmount,
								   l_Stub,
								   C_Currency,
								   l_Adjusted,
								   (int) C_CreditLag,
								   b_IncludeMaturity,
								   C_ProtectionStartDate,
								   C_ProtectionEndDate,
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

//	ARM_END();
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in ARM_Credit_CMCDS" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}
/*---- End Of File ----*/

// EOF %M% 

__declspec(dllexport) LPXLOPER WINAPI ARM_Credit_CAPFLOORCMCDS(LPXLOPER XL_EffectiveDate,
															   LPXLOPER XL_EndDate,
															   LPXLOPER XL_CapLevel,
															   LPXLOPER XL_FloorLevel,
															   LPXLOPER XL_PartRate,
															   LPXLOPER XL_IndexId,
															   LPXLOPER XL_FixingFreq,
															   LPXLOPER XL_DayCountFrq,
															   LPXLOPER XL_First_period_refdate, 	
															   LPXLOPER XL_FixedPayerAmount,	
															   LPXLOPER XL_FloatingPayerAmount,
															   LPXLOPER XL_StubRule,
															   LPXLOPER XL_Currency,
															   LPXLOPER XL_CreditLag,
															   LPXLOPER XL_Adjusted,
															   LPXLOPER XL_IncludeMaturity,
															   LPXLOPER XL_ProtectionStartDate,
															   LPXLOPER XL_ProtectionEndDate)
{
	
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();

	// C variable
	

	double EffectiveDate;
	double EndDate;
	double PartRate;
	double CapLevel ;
	double CapLevel_Default = -1.;
	double FloorLevel ;
	double FloorLevel_Default = -1 ;

	double C_ProtectionStartDate;
	double C_ProtectionEndDate;

	double C_ProtectionStartDate_default = -1.;
	double C_ProtectionEndDate_default = -1.;

	CCString C_FixingFreq;
	int FixingFreq;

	int DayCountFrq;

	double First_period_refdate;
	double First_period_refdate_default = -1.;

	double FixedPayerAmount;	
	double FixedPayerAmount_default = 1000000.;	

	double FloatingPayerAmount;	
	double FloatingPayerAmount_default = -1.;	

	CCString Currency;

	CCString C_Currency;
	CCString C_DayCountFrq;

	CCString C_Adjusted;

	CCString C_IndexId;

	// error
	static int error;
	static char* reason = "";

	CCString C_StubRule;
	CCString C_StubRule_default = "SS";
	int l_Stub;

	double C_CreditLag;
	double C_CreditLag_default = 0.;

	CCString C_IncludeMaturity;
	bool b_IncludeMaturity= false;

	int l_Adjusted = 0;

	XL_readNumCell(XL_EffectiveDate,EffectiveDate," ARM_ERR: Effective Date date expected",C_result);
	XL_readNumCell(XL_EndDate,EndDate," ARM_ERR: EndDate date: date expected",C_result);
	XL_readNumCell(XL_PartRate,PartRate," ARM_ERR: PartRate: numeric expected",C_result);
	XL_readNumCellWD(XL_CapLevel,CapLevel,CapLevel_Default," ARM_ERR: Protection End date: numeric expected",C_result);
	XL_readNumCellWD(XL_FloorLevel,FloorLevel,FloorLevel_Default," ARM_ERR: Protection End date: numeric expected",C_result);
	XL_readStrCellWD(XL_IndexId, C_IndexId,"NONE","ARM_ERR: Index : Object Index expected",C_result);

	XL_readStrCellWD(XL_FixingFreq,C_FixingFreq,"Q"," ARM_ERR: FixingFreq: string expected",C_result);
	XL_readStrCellWD(XL_DayCountFrq,C_DayCountFrq,"A360"," ARM_ERR: DayCountFrq: string expected",C_result);
	
	XL_readNumCellWD(XL_First_period_refdate,First_period_refdate,First_period_refdate_default," ARM_ERR: First_period_refdate Rate: date expected",C_result);
	XL_readNumCellWD(XL_FixedPayerAmount,FixedPayerAmount,FixedPayerAmount_default," ARM_ERR: FixedPayerAmount: numeric expected",C_result);
	XL_readNumCellWD(XL_FloatingPayerAmount,FloatingPayerAmount,FloatingPayerAmount_default," ARM_ERR: FloatingPayerAmount: numeric expected",C_result);
	XL_readStrCellWD(XL_StubRule,C_StubRule,C_StubRule_default," ARM_ERR: Stub Rule: string expected",C_result);
	
	XL_readStrCellWD(XL_Currency,C_Currency,"DEFAULT"," ARM_ERR: currency: string expected",C_result);
	XL_readStrCellWD(XL_Adjusted,C_Adjusted,"ADJ"," ARM_ERR: Adjusted: string expected",C_result);
	XL_readNumCellWD(XL_CreditLag,C_CreditLag,C_CreditLag_default," ARM_ERR: Credit Lag: string expected",C_result);

	XL_readNumCellWD(XL_ProtectionStartDate,C_ProtectionStartDate,C_ProtectionStartDate_default," ARM_ERR: Protection start date: numeric expected",C_result);
	XL_readNumCellWD(XL_ProtectionEndDate,C_ProtectionEndDate,C_ProtectionEndDate_default," ARM_ERR: Protection End date: numeric expected",C_result);

	XL_readStrCellWD(XL_IncludeMaturity,C_IncludeMaturity,"Y"," ARM_ERR: Include Maturity: string expected",C_result);

	if (C_IncludeMaturity == "Y") b_IncludeMaturity = true;

	if (FloatingPayerAmount == FloatingPayerAmount_default)
		FloatingPayerAmount = FixedPayerAmount;

	ARM_result currencyres;
	ARMLOCAL_ARM_GetDefaultCurrency (currencyres);

	if(C_Currency == "DEFAULT")
	{
		if(currencyres.getRetCode () != ARM_OK)
		{
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}
		else
		{
			C_Currency = currencyres.getString ();
		}
	}

	long retCode;

	if((l_Stub = ARM_ConvStubRule (C_StubRule)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	if((FixingFreq = ARM_ConvFrequency (C_FixingFreq, C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	if((DayCountFrq = ARM_ConvDayCount (C_DayCountFrq)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	if((l_Adjusted = ARM_ConvIntRule (C_Adjusted)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	long objId;
	CCString prevClass;
	
	CCString curClass = LOCAL_CMCDS_CLASS;
	CCString stringId = GetLastCurCellEnvValue ();
	
	if(!stringId)
	{	
		retCode = ICMLOCAL_CAPFLOORCMCDS (PartRate,
										  EffectiveDate,
										  EndDate,
										  CapLevel,
										  FloorLevel,
										  LocalGetNumObjectId(C_IndexId),
										  First_period_refdate,
										  FixingFreq,
										  DayCountFrq,
										  FixedPayerAmount,
										  FloatingPayerAmount,
										  l_Stub,
										  C_Currency,
										  l_Adjusted,
										  (int) C_CreditLag,
										  b_IncludeMaturity,
										  C_ProtectionStartDate,
										  C_ProtectionEndDate,
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
		retCode = ICMLOCAL_CAPFLOORCMCDS (PartRate,
										  EffectiveDate,
										  EndDate,
										  CapLevel,
										  FloorLevel,
										  LocalGetNumObjectId(C_IndexId),
										  First_period_refdate,
										  FixingFreq,
										  DayCountFrq,
										  FixedPayerAmount,
										  FloatingPayerAmount,
										  l_Stub,
										  C_Currency,
										  l_Adjusted,
										  (int) C_CreditLag,
										  b_IncludeMaturity,
										  C_ProtectionStartDate,
										  C_ProtectionEndDate,
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

		retCode = ICMLOCAL_CAPFLOORCMCDS(PartRate,
										 EffectiveDate,
										 EndDate,
										 CapLevel,
										 FloorLevel,	
										 LocalGetNumObjectId(C_IndexId),
										 First_period_refdate,
										 FixingFreq,
										 DayCountFrq,
										 FixedPayerAmount,
										 FloatingPayerAmount,
										 l_Stub,
										 C_Currency,
										 l_Adjusted,
										 (int) C_CreditLag,
										 b_IncludeMaturity,
										 C_ProtectionStartDate,
										 C_ProtectionEndDate,
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

//	ARM_END();
	}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in ARM_Credit_CAPFLOORCMCDS" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}



__declspec(dllexport) LPXLOPER WINAPI ARM_Credit_CPPI (LPXLOPER XL_StartDate,
													   LPXLOPER XL_EndDate,
													   LPXLOPER XL_SecurityId,
													   LPXLOPER XL_CorrelName,
													   LPXLOPER XL_Min,
													   LPXLOPER XL_Max,
													   LPXLOPER XL_ValueMin,
													   LPXLOPER XL_ValueMax,
													   LPXLOPER XL_Notional,
													   LPXLOPER XL_ProtectedAmount,
													   LPXLOPER XL_AdditionalLeverage,
													   LPXLOPER XL_DesactivateCushion,
													   LPXLOPER XL_ManagementCost,
													   LPXLOPER XL_Currency)
{
	// return
	static XLOPER XL_result;
	ARM_result C_result;

	// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
		ARM_NOCALCIFWIZ();

			// error
			static int error;
			static char* reason = "";
			long objId;
			CCString prevClass;

		// C variable
		double StartDate;
		double EndDate;

		CCString C_SecurityId;

		vector<double> Min;
		Min.clear();
		vector<double> Max;
		Max.clear();
		vector<double> ValueMin;
		ValueMin.clear();
		vector<double> ValueMax;
		ValueMax.clear();

		double Notional;
		
		double ProtectedAmount;
		
		double AdditionalLeverage;
		double AdditionalLeverageDefaultValue = 0.0;
		
		double DesactivateCushion;
		double DesactivateCushionDefaultValue = 0.0;
		
		double ManagementCost;
		double ManagementCostDefaultValue = 0.0;
		
		CCString C_Currency;

		CCString CorrelName;

		CCString curClass = LOCAL_CPPI_CLASS;
		CCString stringId = GetLastCurCellEnvValue ();
		long retCode;

		XL_readNumCell(XL_StartDate,StartDate,"ARM_ERR : StartDate date expected",C_result);
		XL_readNumCell(XL_EndDate,EndDate,"ARM_ERR : EndDate date expected",C_result);
		XL_readStrCellWD(XL_SecurityId,C_SecurityId,"NONE","ARM_ERR: Security : Object Security expected",C_result);
		XL_readNumVector(XL_Min,Min,"ARM_ERR : Min object vector expected",C_result);
		XL_readNumVector(XL_Max,Max,"ARM_ERR : Max object vector expected",C_result);
		XL_readNumVector(XL_ValueMin,ValueMin,"ARM_ERR : Value Min object vector expected",C_result);
		XL_readNumVector(XL_ValueMax,ValueMax,"ARM_ERR : Value Max object vector expected",C_result);
		XL_readNumCell(XL_Notional,Notional,"ARM_ERR : Notional numeric expected",C_result);
		XL_readNumCell(XL_ProtectedAmount,ProtectedAmount,"ARM_ERR : ProtectedAmount numeric expected",C_result);
		XL_readNumCellWD(XL_AdditionalLeverage,AdditionalLeverage,AdditionalLeverageDefaultValue,"ARM_ERR : AdditionalLeverage numeric expected",C_result);
		XL_readNumCellWD(XL_DesactivateCushion,DesactivateCushion,DesactivateCushionDefaultValue,"ARM_ERR : Desactivate Cushion numeric expected",C_result);
		XL_readNumCellWD(XL_ManagementCost,ManagementCost,ManagementCostDefaultValue,"ARM_ERR : Management Cost numeric expected",C_result);
		XL_readStrCellWD(XL_Currency,C_Currency,"DEFAULT"," ARM_ERR: currency: string expected",C_result);
		XL_readStrCellWD(XL_CorrelName,CorrelName,"TRAXX"," ARM_ERR: Correl Name: string expected",C_result);

		ARM_result currencyres;
		ARMLOCAL_ARM_GetDefaultCurrency (currencyres);
		if(C_Currency == "DEFAULT")
		{
			if(currencyres.getRetCode () != ARM_OK)
			{
				ARM_ARG_ERR();
				return (LPXLOPER)&XL_result;
			}
			else
			{
				C_Currency = currencyres.getString ();
			}
		}

	if(!stringId)
	{	
		retCode = ICMLOCAL_CPPI(StartDate,
								EndDate,
								LocalGetNumObjectId(C_SecurityId),
								C_Currency,
								Min,
								Max,
								ValueMin,
								ValueMax,
								Notional,
								ProtectedAmount,
								ManagementCost,
								AdditionalLeverage,
								DesactivateCushion,
								CorrelName,
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
			retCode = ICMLOCAL_CPPI(StartDate,
									EndDate,
									LocalGetNumObjectId(C_SecurityId),
									C_Currency,
									Min,
									Max,
									ValueMin,
									ValueMax,
									Notional,
									ProtectedAmount,
									ManagementCost,
									AdditionalLeverage,
									DesactivateCushion,
									CorrelName,
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
			retCode = ICMLOCAL_CPPI(StartDate,
									EndDate,
									LocalGetNumObjectId(C_SecurityId),
									C_Currency,
									Min,
									Max,
									ValueMin,
									ValueMax,
									Notional,
									ProtectedAmount,
									ManagementCost,
									AdditionalLeverage,
									DesactivateCushion,
									CorrelName,
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

	Min.clear();
	Max.clear();
	ValueMin.clear();
	ValueMax.clear();
}
	/// end of try block
	ARM_XL_TRY_BLOCK_END

	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in ARM_Credit_CPPI" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


__declspec(dllexport) LPXLOPER WINAPI ARM_Credit_NTDGEN(LPXLOPER XL_CdsId,
														LPXLOPER XL_FirstNumDefault,
														LPXLOPER XL_LastNumDefault,
														LPXLOPER XL_CollateralId,
														LPXLOPER XL_Binary,
														LPXLOPER XL_RcvFee)
{
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	try {

	ARM_NOCALCIFWIZ();

	// C variable
	
	CCString C_CdsId;
	CCString C_CollateralId;
	CCString C_RcvFee;
	CCString C_RcvFee_default="R";

	double Binary = 0.;
	double Binary_default = -999.;
	double FirstNumDefault = 0.;
	double LastNumDefault = 0.;
	double d_RcvFee=-1.;

	long retCode;

	// error
	static int error;
	static char* reason = "";

	XL_readStrCell(XL_CdsId,C_CdsId," ARM_ERR: Cds: string expected",C_result);
	XL_readNumCell(XL_FirstNumDefault,FirstNumDefault," ARM_ERR: FirstNumDefault: numeric expected",C_result);
	XL_readNumCell(XL_LastNumDefault,LastNumDefault," ARM_ERR: LastNumDefault: numeric expected",C_result);
	XL_readStrCell(XL_CollateralId,C_CollateralId," ARM_ERR: Collateral: string expected",C_result);
	XL_readNumCellWD(XL_Binary,Binary,Binary_default," ARM_ERR: Binary: numeric expected",C_result);
	XL_readStrCellWD(XL_RcvFee,C_RcvFee,C_RcvFee_default," ARM_ERR: RcvFee: string expected",C_result);
	
	if (C_RcvFee=="R") d_RcvFee = 1.;

	long objId;
	CCString prevClass;
	
	CCString curClass = LOCAL_NTHTD_CLASS;
	CCString stringId = GetLastCurCellEnvValue ();
	
	if(!stringId)
	{	
		retCode = ICMLOCAL_NTDGEN(LocalGetNumObjectId (C_CdsId),
								  (int) FirstNumDefault,
								  (int) LastNumDefault,
								  LocalGetNumObjectId (C_CollateralId),	
								  Binary,
								  d_RcvFee,
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
		retCode = ICMLOCAL_NTDGEN(LocalGetNumObjectId (C_CdsId),
								  (int) FirstNumDefault,
								  (int) LastNumDefault,
								  LocalGetNumObjectId (C_CollateralId),	
								  Binary,
								  d_RcvFee,
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

		retCode = ICMLOCAL_NTDGEN(LocalGetNumObjectId (C_CdsId),
								  (int) FirstNumDefault,
								  (int) LastNumDefault,
								  LocalGetNumObjectId (C_CollateralId),	
								  Binary,
								  d_RcvFee,
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

//	ARM_END();
	}
	catch (Exception&e)
	{
		ExcelCaller::get().setError(e.GetErrorString()); 
		ExcelTools::convert("ARM_ERR",&XL_result) ;
	}
	catch (std::exception&e)
	{
		ExcelCaller::get().setError(e.what()); 
		ExcelTools::convert("ARM_ERR",&XL_result) ;
	}
	catch(...)
	{
		ExcelCaller::get().setError("Unknown Exception");
		ExcelTools::convert("ARM_ERR",&XL_result) ;
	}

	return (LPXLOPER)&XL_result;
}



__declspec(dllexport) LPXLOPER WINAPI ARM_Credit_CDOGEN(LPXLOPER XL_CdsId,
														LPXLOPER XL_SubAmount,
														LPXLOPER XL_CollateralId,
														LPXLOPER XL_Binary,
														LPXLOPER XL_RcvFee)
{
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	try {
	ARM_NOCALCIFWIZ();

	// C variable
	
	CCString C_CdsId;
	CCString C_CollateralId;
	CCString C_RcvFee;
	CCString C_RcvFee_default="R";

	double Binary = 0.;
	double Binary_default = -999.;
	long SubAmount_type;
	double C_SubAmount_double = 0.;
	CCString C_SubAmount_str = 0.;
	double d_RcvFee=-1.;

	long retCode;

	// error
	static int error;
	static char* reason = "";

	XL_readStrCell(XL_CdsId,C_CdsId," ARM_ERR: Cds: string expected",C_result);
	XL_readStrOrNumCell(XL_SubAmount,C_SubAmount_str,C_SubAmount_double,SubAmount_type," ARM_ERR: subamount: string or numeric expected",C_result);
	XL_readStrCell(XL_CollateralId,C_CollateralId," ARM_ERR: Collateral: string expected",C_result);
	XL_readNumCellWD(XL_Binary,Binary,Binary_default," ARM_ERR: Binary: numeric expected",C_result);
	XL_readStrCellWD(XL_RcvFee,C_RcvFee,C_RcvFee_default," ARM_ERR: RcvFee: string expected",C_result);
	
	if ( SubAmount_type == XL_TYPE_STRING )
	{
	   C_SubAmount_double = (double) LocalGetNumObjectId(C_SubAmount_str);
	   SubAmount_type = 1;
	}
	else
	{  SubAmount_type = 0;}
	
	if (C_RcvFee=="R") d_RcvFee = 1.;

	long objId;
	CCString prevClass;
	
	CCString curClass = LOCAL_MEZ_CLASS;
	CCString stringId = GetLastCurCellEnvValue ();
	
	if(!stringId)
	{	
		retCode = ICMLOCAL_CDOGEN(LocalGetNumObjectId (C_CdsId),
								  C_SubAmount_double,
								  SubAmount_type,
								  LocalGetNumObjectId (C_CollateralId),	
								  Binary,
								  d_RcvFee,
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
		retCode = ICMLOCAL_CDOGEN(LocalGetNumObjectId (C_CdsId),
								  C_SubAmount_double,
								  SubAmount_type,
								  LocalGetNumObjectId (C_CollateralId),	
								  Binary,
								  d_RcvFee,
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

		retCode = ICMLOCAL_CDOGEN(LocalGetNumObjectId (C_CdsId),
								  C_SubAmount_double,
								  SubAmount_type,
								  LocalGetNumObjectId (C_CollateralId),	
								  Binary,
								  d_RcvFee,
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

//	ARM_END();
	}
	catch (Exception&e)
	{
		ExcelCaller::get().setError(e.GetErrorString()); 
		ExcelTools::convert("ARM_ERR",&XL_result) ;
	}
	catch (std::exception&e)
	{
		ExcelCaller::get().setError(e.what()); 
		ExcelTools::convert("ARM_ERR",&XL_result) ;
	}
	catch(...)
	{
		ExcelCaller::get().setError("Unknown Exception");
		ExcelTools::convert("ARM_ERR",&XL_result) ;
	}

	return (LPXLOPER)&XL_result;
}

__declspec(dllexport) LPXLOPER WINAPI ARM_Credit_Customized_CDO(
														LPXLOPER XL_Labels,
														LPXLOPER XL_Notionals,
														LPXLOPER XL_Currency,
														LPXLOPER XL_CreditProductDefaultLeg, 
														LPXLOPER XL_CreditProductPremiumLeg,
														LPXLOPER XL_CreditProductParameters
														)

{
	// return
	static XLOPER XL_result;
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
		ARM_NOCALCIFWIZ();

		// C variable

		CCString	C_CreditProductDefaultLegId;
		CCString	C_CreditProductPremiumLegId;
		CCString	C_CreditProductParametersId;
		CCString	C_Currency;


		VECTOR<CCString> C_Labels;
		VECTOR<double> C_Notionals;
		
		// error
		static int error;
		static char* reason = "";

		// LABELS
		XL_readStrVector(XL_Labels, C_Labels," ARM_ERR: Issuers Labels : string vector expected", DOUBLE_TYPE, C_result);

		// Notionals
		XL_readNumVector(XL_Notionals, C_Notionals ," ARM_ERR: Notionals : numeric vector expected",C_result);
	
		//Currency
		XL_readStrCellWD(XL_Currency,C_Currency,"DEFAULT"," ARM_ERR: currency: string expected",C_result);

		// Credit Product Default Leg
		XL_readStrCell(XL_CreditProductDefaultLeg, C_CreditProductDefaultLegId ," ARM_ERR: Credit Product Default Leg Id: object expected",C_result);

		// Credit Product Premium Leg
		XL_readStrCell(XL_CreditProductPremiumLeg, C_CreditProductPremiumLegId ," ARM_ERR: Credit Product Premium Leg Id: object expected",C_result);

		// Credit Product Pricing Parameters
		XL_readStrCell(XL_CreditProductParameters, C_CreditProductParametersId ," ARM_ERR: Credit Product Pricing Parameters Id: object expected",C_result);

		//Currency
		ARM_result currencyres;
		ARMLOCAL_ARM_GetDefaultCurrency (currencyres);

		if(C_Currency == "DEFAULT")
		{
			if(currencyres.getRetCode () != ARM_OK)
			{
				ARM_ARG_ERR();
				return (LPXLOPER)&XL_result;
			}
			else
			{
				C_Currency = currencyres.getString ();
			}
		}

		long retCode;

		long objId;
		CCString prevClass;

		CCString curClass = LOCAL_CUSTOMIZED_CDO_CLASS;
		CCString stringId = GetLastCurCellEnvValue ();

		if(!stringId)
		{	
			retCode = ICMLOCAL_Customized_CDO(
										C_Labels,
										C_Notionals,
										C_Currency,
										LocalGetNumObjectId (C_CreditProductDefaultLegId),
										LocalGetNumObjectId (C_CreditProductPremiumLegId),
										LocalGetNumObjectId (C_CreditProductParametersId),
										C_result);
	
			if (retCode == ARM_OK)
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
				retCode = ICMLOCAL_Customized_CDO(
										C_Labels,
										C_Notionals,
										C_Currency,
										LocalGetNumObjectId (C_CreditProductDefaultLegId),
										LocalGetNumObjectId (C_CreditProductPremiumLegId),
										LocalGetNumObjectId (C_CreditProductParametersId),
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

				retCode = ICMLOCAL_Customized_CDO(
										C_Labels,
										C_Notionals,
										C_Currency,
										LocalGetNumObjectId (C_CreditProductDefaultLegId),
										LocalGetNumObjectId (C_CreditProductPremiumLegId),
										LocalGetNumObjectId (C_CreditProductParametersId),
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG("Unrecognized failure in ARM_Credit_Customized_CDO")

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


/**
__declspec(dllexport) LPXLOPER WINAPI ARM_Credit_Index_Option_Gen(
														LPXLOPER XL_CreditScheduleParameters, 
														LPXLOPER XL_CreditDataParameters
																)
{
	// return
	static XLOPER XL_result;
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
		ARM_NOCALCIFWIZ();

		// C variable

		CCString	C_CreditScheduleParametersId;
		CCString	C_CreditDataParametersId;

		// error
		static int error;
		static char* reason = "";
	
		// Credit Schedule Parameters
		XL_readStrCell(XL_CreditScheduleParameters, C_CreditScheduleParametersId ," ARM_ERR: Credit Index Option Gen Id - Schedule: object expected",C_result);

		// Credit Data Parameters
		XL_readStrCell(XL_CreditDataParameters, C_CreditDataParametersId ," ARM_ERR: Credit Index Option Gen Id - Data: object expected",C_result);

		long retCode;

		long objId;
		CCString prevClass;

		CCString curClass = LOCAL_CREDIT_INDEX_OPTION_GEN;
		CCString stringId = GetLastCurCellEnvValue ();

		if(!stringId)
		{	
			retCode = ICMLOCAL_Option_Index_Gen(
										LocalGetNumObjectId (C_CreditScheduleParametersId),
										LocalGetNumObjectId (C_CreditDataParametersId),
										C_result);
	
			if (retCode == ARM_OK)
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
				retCode = ICMLOCAL_Option_Index_Gen(
										LocalGetNumObjectId (C_CreditScheduleParametersId),
										LocalGetNumObjectId (C_CreditDataParametersId),
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

				retCode = ICMLOCAL_Option_Index_Gen(
										LocalGetNumObjectId (C_CreditScheduleParametersId),
										LocalGetNumObjectId (C_CreditDataParametersId),
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG("Unrecognized failure in ARM_Credit_Index_Option_Gen")

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}

__declspec(dllexport) LPXLOPER WINAPI ARM_Credit_Index_Corridor_Gen(
														LPXLOPER XL_CreditScheduleParameters, 
														LPXLOPER XL_CreditDataParameters,
														LPXLOPER XL_CreditSubScheduleParameters
																)
{
	// return
	static XLOPER XL_result;
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
		ARM_NOCALCIFWIZ();

		// C variable
		CCString	C_CreditScheduleParametersId;
		CCString	C_CreditDataParametersId;
		CCString	C_CreditSubScheduleParametersId;

		// error
		static int error;
		static char* reason = "";
	
		// Credit Schedule Parameters
		XL_readStrCell(XL_CreditScheduleParameters, C_CreditScheduleParametersId ," ARM_ERR: Credit Index Corridor Gen Id - Schedule: object expected",C_result);

		// Credit Data Parameters
		XL_readStrCell(XL_CreditDataParameters, C_CreditDataParametersId ," ARM_ERR: Credit Index Corridor Gen Id - Data: object expected",C_result);

		// Credit Sub Schedule Parameters
		XL_readStrCell(XL_CreditSubScheduleParameters, C_CreditSubScheduleParametersId ," ARM_ERR: Credit Index Corridor Gen Id - Sub Schedule: object expected",C_result);

		long retCode;

		long objId;
		CCString prevClass;

		CCString curClass = LOCAL_CREDIT_INDEX_CORRIDOR_GEN;
		CCString stringId = GetLastCurCellEnvValue ();

		if(!stringId)
		{	
			retCode = ICMLOCAL_Corridor_Index_Gen(
										LocalGetNumObjectId (C_CreditScheduleParametersId),
										LocalGetNumObjectId (C_CreditDataParametersId),
										LocalGetNumObjectId (C_CreditSubScheduleParametersId),
										C_result);
	
			if (retCode == ARM_OK)
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
				retCode = ICMLOCAL_Corridor_Index_Gen(
										LocalGetNumObjectId (C_CreditScheduleParametersId),
										LocalGetNumObjectId (C_CreditDataParametersId),
										LocalGetNumObjectId (C_CreditSubScheduleParametersId),
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

				retCode = ICMLOCAL_Corridor_Index_Gen(
										LocalGetNumObjectId (C_CreditScheduleParametersId),
										LocalGetNumObjectId (C_CreditDataParametersId),
										LocalGetNumObjectId (C_CreditSubScheduleParametersId),
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG("Unrecognized failure in ARM_Credit_Index_Corridor_Gen")

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}
**/ 

_declspec(dllexport) LPXLOPER WINAPI ARM_Credit_Instrument_GetDataMatrix(
														LPXLOPER XL_InstrumentId
														)
{
	// return
	static XLOPER XL_result;
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
		ARM_NOCALCIFWIZ();

	// C variable
	CCString	C_InstrumentId;

	// error
	static int error;
	static char* reason = "";

	/// end of try block
	}
	ARM_XL_TRY_BLOCK_END
	
	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG("Unrecognized failure in ARM_Credit_Instrument_GetDataMatrix")

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}
 

_declspec(dllexport) LPXLOPER WINAPI ARM_Credit_GetLastFixingDate(
														LPXLOPER XL_InstrumentId,
														LPXLOPER XL_AsOfDate
														)
{
	// return
	static XLOPER XL_result;
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
		ARM_NOCALCIFWIZ();

		// C variable
		CCString	C_InstrumentId;
		double AsofDate ;

		// error
		static int error;
		static char* reason = "";

		// 
		XL_readStrCell(XL_InstrumentId,C_InstrumentId," ARM_ERR: Instrument Id: object expected",C_result);
		XL_readNumCell(XL_AsOfDate,AsofDate," ARM_ERR: AsofDate date: date expected",C_result);

		long retCode;
		long instId=LocalGetNumObjectId(C_InstrumentId); 
		retCode = ICMLOCAL_getLastFixingDate(instId,AsofDate,C_result) ;
		if (retCode==ARM_OK) 
		{
			FreeCurCellErr ();
			XL_result.xltype = xltypeNum;
			XL_result.val.num = Local_ARMDATE2XLDATE(C_result.getString ()) ; 
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG("Unrecognized failure in ARM_Credit_GetLastFixingDate")

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}



_declspec(dllexport) LPXLOPER WINAPI ARM_Credit_SetPastFixing(
														LPXLOPER XL_InstrumentId,
														LPXLOPER XL_resetDate,
														LPXLOPER XL_fixingValue
														)
{
	// return
	static XLOPER XL_result;
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
		ARM_NOCALCIFWIZ();

	// C variable
	CCString	C_InstrumentId;
	double		resetDate ;
	double		fixingValue ;

	// error
	static int error;
	static char* reason = "";

	// 
	XL_readStrCell(XL_InstrumentId,C_InstrumentId," ARM_ERR: Instrument Id: object expected",C_result);
	XL_readNumCell(XL_resetDate,resetDate," ARM_ERR: resetDate date: date expected",C_result);
	XL_readNumCell(XL_fixingValue,fixingValue," ARM_ERR: fixingValue : numeric expected",C_result);

	long retCode;
	long instId = LocalGetNumObjectId(C_InstrumentId);

	retCode = ICMLOCAL_setPastFixing(instId, resetDate,fixingValue,C_result);
	if (retCode == ARM_OK)
	{
		FreeCurCellErr ();
		XL_result.xltype = xltypeStr;
		XL_result.val.str = XL_StrC2StrPascal (C_InstrumentId);
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG("Unrecognized failure in ARM_Credit_SetPastFixing")

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


_declspec(dllexport) LPXLOPER WINAPI ARM_Credit_GetEqStrikeUp(
														LPXLOPER XL_CorrelId,
														LPXLOPER XL_IndexName
														)
{
	// return
	static XLOPER XL_result;
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
		ARM_NOCALCIFWIZ();

		// C variable
		CCString	C_CorrelId;
		CCString	C_IndexName;

		// error
		static int error;
		static char* reason = "";

		// 
		XL_readStrCell(XL_CorrelId,C_CorrelId," ARM_ERR: Instrument Id: object expected",C_result);
		XL_readStrCell(XL_IndexName,C_IndexName," ARM_ERR: Instrument Id: object expected",C_result);
		

		long retCode;
		long pricerId=LocalGetNumObjectId(C_CorrelId); 
		std::vector<double> matu,strikes; 
		retCode = ICMLOCAL_GetEqStrike ( pricerId,
							C_IndexName,
							1, /** is Up **/ 
							matu,
							strikes,
							C_result) ;

		if (retCode==ARM_OK) 
		{
			FreeCurCellErr ();
			XL_result.xltype = xltypeNum;
			XL_result.val.num = strikes[0] ; // C_result.getDouble() ; 
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG("Unrecognized failure in ARM_Credit_GetEqStrikeUp")

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}

_declspec(dllexport) LPXLOPER WINAPI ARM_Credit_GetEqStrikeDown(
														LPXLOPER XL_PricerId,
														LPXLOPER XL_IndexName
														)
{
	// return
	static XLOPER XL_result;
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
		ARM_NOCALCIFWIZ();

		// C variable
		CCString	C_PricerId;
		CCString	C_IndexName;

		// error
		static int error;
		static char* reason = "";

		// 
		XL_readStrCell(XL_PricerId,C_PricerId," ARM_ERR: Instrument Id: object expected",C_result);
		XL_readStrCell(XL_IndexName,C_IndexName," ARM_ERR: Instrument Id: object expected",C_result);
		

		long retCode;
		long pricerId=LocalGetNumObjectId(C_PricerId); 
		// retCode = ICMLOCAL_GetEqStrikeDown(pricerId,C_IndexName,C_result) ;
		std::vector<double> matu,strikes; 
		retCode = ICMLOCAL_GetEqStrike ( pricerId,
							C_IndexName,
							0, /** is Up **/ 
							matu,
							strikes,
							C_result) ;
		if (retCode==ARM_OK) 
		{
			FreeCurCellErr ();
			XL_result.xltype = xltypeNum;
			XL_result.val.num = strikes[0] ; // C_result.getDouble() ; 
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG("Unrecognized failure in ARM_Credit_GetEqStrikeDown")

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}

_declspec(dllexport) LPXLOPER WINAPI ARM_Credit_LssGapOption(LPXLOPER XL_MezzId,
															LPXLOPER XL_SpreadsTrigger,
															LPXLOPER XL_Defaulttrigger,
															LPXLOPER XL_Mtmtrigger,
															LPXLOPER XL_mtm_single_cond)
{
	// return
	static XLOPER XL_result;
	ARM_result C_result;

	long retCode;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
		ARM_NOCALCIFWIZ();

		// C variable
		CCString	C_MezzId;
		VECTOR<double>	V_SpreadsTrigger;
		VECTOR<double>	V_Defaulttrigger;
		VECTOR<double>	V_Mtmtrigger;
		double	d_mtm_single_cond ;
		// error
		static int error;
		static char* reason = "";

		// 
		XL_readStrCell(XL_MezzId,C_MezzId," ARM_ERR: Instrument Id: object expected",C_result);

		long nbrows_spreads,nbcolumns_spreads;
		XL_readNumVectorAndSize(XL_SpreadsTrigger,nbrows_spreads,nbcolumns_spreads,V_SpreadsTrigger," ARM_ERR: SpreadsTrigger: object expected",C_result);

		long nbrows_default,nbcolumns_default;
		XL_readNumVectorAndSize(XL_Defaulttrigger,nbrows_default,nbcolumns_default,V_Defaulttrigger," ARM_ERR: Defaulttrigger: object expected",C_result);

		long nbrows_Mtmtrigger,nbcolumns_Mtmtrigger;
		XL_readNumVectorAndSize(XL_Mtmtrigger,nbrows_Mtmtrigger,nbcolumns_Mtmtrigger,V_Mtmtrigger," ARM_ERR: Mtmtrigger: object expected",C_result);

		XL_readNumCell(XL_mtm_single_cond,d_mtm_single_cond," ARM_ERR: mtm_single_cond : numeric",C_result);
		
		CCString prevClass;
		CCString curClass = LOCAL_CREDIT_INDEX_CORRIDOR_GEN;
		CCString stringId = GetLastCurCellEnvValue ();
		long objId;

		if(!stringId)
		{	
			retCode = ICMLOCAL_LssGapOption ( LocalGetNumObjectId(C_MezzId),
							nbrows_spreads,
							nbcolumns_spreads,
							V_SpreadsTrigger,
							nbrows_Mtmtrigger,
							nbcolumns_Mtmtrigger,
							V_Mtmtrigger,
							nbrows_Mtmtrigger,
							nbcolumns_Mtmtrigger,
							V_Mtmtrigger,
							d_mtm_single_cond,
							C_result) ;

	
			if (retCode == ARM_OK)
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
				retCode = ICMLOCAL_LssGapOption ( LocalGetNumObjectId(C_MezzId),
							nbrows_spreads,
							nbcolumns_spreads,
							V_SpreadsTrigger,
							nbrows_Mtmtrigger,
							nbcolumns_Mtmtrigger,
							V_Mtmtrigger,
							nbrows_Mtmtrigger,
							nbcolumns_Mtmtrigger,
							V_Mtmtrigger,
							d_mtm_single_cond,
							C_result,
							objId) ;

			
				if(retCode == ARM_OK)
				{
					LocalSetCurCellEnvValue (curClass, objId); 

					stringId = LocalMakeObjectId (objId, curClass);
				}
			}
			else
			{
				FreeCurCellContent ();

				retCode = ICMLOCAL_LssGapOption ( LocalGetNumObjectId(C_MezzId),
							nbrows_spreads,
							nbcolumns_spreads,
							V_SpreadsTrigger,
							nbrows_Mtmtrigger,
							nbcolumns_Mtmtrigger,
							V_Mtmtrigger,
							nbrows_Mtmtrigger,
							nbcolumns_Mtmtrigger,
							V_Mtmtrigger,
							d_mtm_single_cond,
							C_result) ;
			
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG("Unrecognized failure in ARM_Credit_GetEqStrikeDown")

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}

_declspec(dllexport) LPXLOPER WINAPI ARM_Credit_GetEqStrike(
														LPXLOPER XL_PricerId,
														LPXLOPER XL_IndexName,
														LPXLOPER XL_UpOrLow
														)
{
	// return
	static XLOPER XL_result;
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
		ARM_NOCALCIFWIZ();

		// C variable
		CCString	C_PricerId;
		CCString	C_IndexName;
		CCString	C_UpOrLow ;
		// error
		static int error;
		static char* reason = "";

		// 
		XL_readStrCell(XL_PricerId,C_PricerId," ARM_ERR: Instrument Id: object expected",C_result);
		XL_readStrCell(XL_IndexName,C_IndexName," ARM_ERR: IndexName : string ",C_result);
		XL_readStrCell(XL_UpOrLow,C_UpOrLow," ARM_ERR: UpOrLow : string",C_result);

		int iUpOrLow = 0; 
		if ((iUpOrLow = ARM_ConvTypeUpOrLow(C_UpOrLow,C_result)) == ARM_DEFAULT_ERR ) 
		{
			ARM_ERR() 
		}
		if ( iUpOrLow == 0) iUpOrLow =1;
		else iUpOrLow =0;
		long retCode;
		long pricerId=LocalGetNumObjectId(C_PricerId); 
		// retCode = ICMLOCAL_GetEqStrikeDown(pricerId,C_IndexName,C_result) ;
		std::vector<double> matu,strikes; 
		retCode = ICMLOCAL_GetEqStrike ( pricerId,
							C_IndexName,
							iUpOrLow, /** is Up **/ 
							matu,
							strikes,
							C_result) ;

		if (retCode==ARM_OK) 
		{
			FreeCurCellErr ();
			ICM_QMatrix<double> ret(matu.size(),2); 
			if (!ret.IsEmpty()) 
			{
				ARM_Vector tmpMat(matu); ret.SetCol(0,&tmpMat); 
				ARM_Vector tmpStrikes(strikes); ret.SetCol(1,&tmpStrikes); 
			}
			
			ExcelTools::convert(ret,&XL_result); 
			// XL_result.xltype = xltypeNum;
			// XL_result.val.num = strikes[0] ; // C_result.getDouble() ; 
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
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG("Unrecognized failure in ARM_Credit_GetEqStrikeDown")

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}

_declspec(dllexport) LPXLOPER WINAPI ARM_Credit_Restrikable_CDO(LPXLOPER XL_TriggerStartDate,
																LPXLOPER XL_Expiry,
																LPXLOPER XL_Strike,
																LPXLOPER XL_InitSpread,
																LPXLOPER XL_OptionType,
																LPXLOPER XL_UnderlyingId,
																LPXLOPER XL_Rehauss,
																LPXLOPER XL_TriggerFreq,
																LPXLOPER XL_DiffCDO,
																LPXLOPER XL_CMSpread,
																LPXLOPER XL_CMSpreadMatu)
{
	static XLOPER XL_result;
	ARM_result C_result;

	try
	{
		ARM_NOCALCIFWIZ();
		std::string C_Underlying ; ExcelTools::convert(XL_UnderlyingId,C_Underlying);
		//double UnderlyingMatu ; ExcelTools::convert(XL_UnderlyingMatu,-1.,UnderlyingMatu);
		double ExpiryDate; ExcelTools::convert(XL_Expiry,-1.,ExpiryDate);
		double TriggerStartDate ; ExcelTools::convert(XL_TriggerStartDate,-1.,TriggerStartDate);
		double Strike ; ExcelTools::convert(XL_Strike,-1.,Strike);
		int OptionType ; ExcelTools::convert(XL_OptionType,-999,OptionType);
		double Rehauss ; ExcelTools::convert(XL_Rehauss,-1.,Rehauss);
		double InitSpread ; ExcelTools::convert(XL_InitSpread,-1.,InitSpread);
		std::string  TriggerFreq ; ExcelTools::convert(XL_TriggerFreq,TriggerFreq);
		int DiffCDO ; ExcelTools::convert(XL_CMSpread,-999,DiffCDO);
		int IsCMSpread;ExcelTools::convert(XL_CMSpread,-999,IsCMSpread);
		double CMSpreadMatu;ExcelTools::convert(XL_CMSpreadMatu,0.0,CMSpreadMatu);

		int l_TriggerFreq ;

		if((l_TriggerFreq = ARM_ConvFrequency (CCString(TriggerFreq.c_str()), C_result)) == ARM_DEFAULT_ERR)
		{
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}

		long NewId = ICMLOCAL_RESTRIKABLE_CDO(&ARM_Date(XLDateToJulian(TriggerStartDate)),
											&ARM_Date(XLDateToJulian(ExpiryDate)),
											Strike,
											InitSpread,
											OptionType,
											LocalPersistent::get().getObjectId(C_Underlying),
											Rehauss,
											l_TriggerFreq,
											IsCMSpread,
											CMSpreadMatu,
											DiffCDO);
		std::string ObjName = ExcelCaller::get().setObject(NewId,LOCAL_OPTION_CLASS);
		
		ExcelTools::convert(ObjName,&XL_result);


	}
	catch (Exception&e)
	{
	ExcelCaller::get().setError(e.GetErrorString()); 
	ExcelTools::convert("ARM_ERR",&XL_result) ;
	}
	 catch (std::exception&e)
	 {
	  ExcelCaller::get().setError(e.what()); 
	  ExcelTools::convert("ARM_ERR",&XL_result) ;
	 }
	 catch(...)
	 {
	  ExcelCaller::get().setError("Unknown Exception");
	  ExcelTools::convert("ARM_ERR",&XL_result) ;
	 }
	 return (LPXLOPER)&XL_result;


}


