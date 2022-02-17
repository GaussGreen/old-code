// ARMModule.cpp : Implementation of CLocal_DLLARMApp and DLL registration.

#include "firsttobeincluded.h"
#include "ARM\libarm_frometk\VariantTools.h"
#include "ActiveXModule.h"

#include "ARM_local_class.h"
#include "CCatl.h"


#include "ARM_result.h"
#include "ARM_local_glob.h"

#include "ARM_local_zccurve.h"
#include "ARM_local_mod.h"
#include "ARM_local_swap.h"
#include "ARM_local_class.h"
#include "ARM_local_ccy.h"
#include "ARM_local_volcrv.h"
#include "ICM_local_pwccurve.h" 
#include "ICM_local_glob.h"
#include "ICM_local_swap.h"
#include "ICM_local_mod.h"
#include "ICM_local_pricer.h"
#include "ICM_local_pf.h"
#include "ICM_local_leg.h"
#include "ICM_local_frn.h"

#include "ARM_local_interglob.h"

#include <fromto.h>
#include <currency.h>


STDMETHODIMP ActiveXModule::ARM_Credit_CDS(double	pEffectiveDate,
									   double	pEndDate, 
									   double	pSpread, 
									   BSTR		pFixingFreq, 
									   BSTR		pDayCount, 
									   double	pFirst_period_refdate, 
									   double	pFixedPayerAmount, 
									   double	pFloatingPayerAmount, 
									   BSTR		StubRule,
									   BSTR     pCurrency,
									   BSTR		Adjusted,
									   int		l_CreditLag,
									   BSTR	    pIncludeMaturity,  
									   double	pProtectionStartDate,
									   double	pProtectionEndDate, 
									   BSTR		IssuerName,
									   double   Binary,
									   double	pFstCpnEffDate,
									   BSTR		StartAdj,
									   VARIANT *pRet)
{
	// TODO: Add your implementation code here

	ARM_result C_result;
	try{
	long size =0;

	_bstr_t b_Adjusted(Adjusted);
	CCString C_Adjusted = b_Adjusted;
	int l_Adjusted = ARM_ConvIntRule (C_Adjusted);

	_bstr_t b_StartAdj(StartAdj);
	CCString C_StartAdj = b_StartAdj;
	int l_StartAdj = ARM_ConvStartAdjRule ((const string) C_StartAdj.c_str());

	_bstr_t b_StubRule(StubRule);
	CCString C_StubRule = b_StubRule;
	int l_StubRule = ARM_ConvStubRule(C_StubRule);

	_bstr_t b_IncludeMaturity(pIncludeMaturity);
	CCString C_IncludeMaturity = b_IncludeMaturity;
	bool B_IncludeMaturity = false;
	if (C_IncludeMaturity == "Y") B_IncludeMaturity = true;

	_bstr_t b_FixingFreq(pFixingFreq);
	CCString C_FixingFreq = b_FixingFreq;

	_bstr_t b_IssuerName(IssuerName);
	CCString C_IssuerName = b_IssuerName;

	int l_FixingFreq;
	if((l_FixingFreq = ARM_ConvFrequency (C_FixingFreq, C_result)) == ARM_DEFAULT_ERR)
	{
		ERROR_MSG("Invalid FixingFreq",pRet,ARM_ERROR_FREQ);
		return S_OK;
	}

	_bstr_t b_DayCount(pDayCount);
	CCString C_DayCount = b_DayCount;

	long l_DayCountFrq = ARM_ConvDayCount (C_DayCount);

	_bstr_t b_Currency(pCurrency);
	CCString C_Currency = b_Currency;

	CCString Res;
	CCString curClass = LOCAL_CDS_CLASS;


	long newId = ICMLOCAL_CDS	(pSpread,
					 pEffectiveDate,
					 pEndDate,
					 pFirst_period_refdate,
					 pFstCpnEffDate,
					 l_FixingFreq,
					 l_DayCountFrq,
					 pFixedPayerAmount,
					 pFloatingPayerAmount,
					 l_StubRule,
					 C_Currency,
					 l_Adjusted,
					 l_CreditLag,
					 B_IncludeMaturity,
					 pProtectionStartDate,
					 pProtectionEndDate,
					 C_IssuerName,
					 Binary,
					 l_StartAdj,
					 qACCRUED_SETTLED,
						K_NX_NONE,
					-2,
					 C_result);
	string newLabel= LocalPersistent::get().getStringId(newId, LOCAL_CDS_CLASS);
	VariantTools::convert(newLabel, *pRet); 
	return S_OK ;
	}
	catch(Exception&e)
	{
		return createErrorInfo("",e.GetErrorString()); 
	}
	catch(std::exception&e)
	{
		return createErrorInfo("",e.what()); 
	}
	catch(...) 
	{ 
		return createErrorInfo("","Unknown exception"); 
	}
	return E_FAIL;
}


STDMETHODIMP ActiveXModule::ARM_Credit_FTD(VARIANT *pEffectiveDate,
									   VARIANT *pEndDate, 
									   VARIANT *pSpread, 
									   VARIANT *pLabels, 
									   VARIANT *pFixingFreq, 
									   VARIANT *pDayCountFrq, 
									   VARIANT *pFirst_period_refdate, 
									   VARIANT *pIssuerNotional, 
									   VARIANT *pAccruedOnDefault, 	
									   VARIANT *pCurrency,
									   VARIANT *pPayCreditLag,
									   VARIANT *pStub,	
									   double	pFstCpnEffDate,
									   VARIANT *pintRule,
									   VARIANT *pstartAdj,
									   VARIANT *pRet)
{
	// TODO: Add your implementation code here
	try
	{
		_variant_t v_EffectiveDate (pEffectiveDate);
		v_EffectiveDate.ChangeType( VT_R8 );
		double C_EffectiveDate = (double) v_EffectiveDate;

		_variant_t v_EndDate (pEndDate);
		v_EndDate.ChangeType( VT_R8 );
		double C_EndDate = (double) v_EndDate;

		_variant_t v_Spread (pSpread);
		v_Spread.ChangeType( VT_R8 );
		double C_Spread = (double) v_Spread;

		ARM_result C_result;

		VECTOR<CCString> V_labels;
		long size =0;

		CCString C_FixingFreq;
		if(VARIANT2CCString (*pFixingFreq, C_FixingFreq) != S_OK)
			return S_FALSE;
		int l_FixingFreq;
		if((l_FixingFreq = ARM_ConvFrequency (C_FixingFreq, C_result)) == ARM_DEFAULT_ERR)
		{
			ERROR_MSG("Invalid FixingFreq",pRet,ARM_ERROR_FREQ);
			return S_OK;
		}

		CCString C_intRule;
		if(VARIANT2CCString (*pintRule, C_intRule) != S_OK)
			return S_FALSE;
		int l_intRule;
		if((l_intRule = ARM_ConvIntRule (C_intRule)) == ARM_DEFAULT_ERR)
		{
			ERROR_MSG("Invalid intRule",pRet,ARM_ERROR_FREQ);
			return S_OK;
		}

		CCString C_startAdj;
		if(VARIANT2CCString (*pstartAdj, C_startAdj) != S_OK)
			return S_FALSE;
		int l_startAdj;
		if((l_startAdj = ARM_ConvStartAdjRule ((const string) C_startAdj.c_str())) == ARM_DEFAULT_ERR)
		{
			ERROR_MSG("Invalid startAdj",pRet,ARM_ERROR_FREQ);
			return S_OK;
		}

		CCString C_DayCountFrq;
		if(VARIANT2CCString (*pDayCountFrq, C_DayCountFrq) != S_OK)
			return S_FALSE;

		long l_DayCountFrq = ARM_ConvDayCount (C_DayCountFrq);

		_variant_t v_First_period_refdate(pFirst_period_refdate);
		v_First_period_refdate.ChangeType( VT_R8 );
		double C_First_period_refdate = (double) v_First_period_refdate;

		_variant_t v_IssuerNotional(pIssuerNotional);
		v_IssuerNotional.ChangeType( VT_R8 );
		double C_IssuerNotional = (double) v_IssuerNotional;

		// CCString C_AccruedOnDefault;
		// int l_AccruedOnDefault;
		// if(VARIANT2CCString (*pAccruedOnDefault, C_AccruedOnDefault) != S_OK)
		// 	return S_FALSE;
		// 
		// if((l_AccruedOnDefault = ARM_ConvAccOnDef (C_AccruedOnDefault, C_result)) == ARM_DEFAULT_ERR)
		// {
		// 	ERROR_MSG("Invalid AccruedOnDefault",pRet,ARM_ERROR_ACCONDEF);
		// 	return S_OK;
		// }
		qPAYMENT_PREMIUM_LEG l_AccruedOnDefault; 
		std::string strAccOnDef; 
		VariantTools::convert(*pAccruedOnDefault,strAccOnDef); 
		ICM_EnumsCnv::cnv(strAccOnDef,l_AccruedOnDefault); 

		CCString C_Currency;

		if(VARIANT2CCString (*pCurrency, C_Currency) != S_OK)
			return S_FALSE;

		if(VARIANT2VECTORCCSTRING(*pLabels, V_labels,size) != S_OK)
			return S_FALSE;

		_variant_t v_PayCreditLag(pPayCreditLag);
		v_PayCreditLag.ChangeType( VT_R8 );
		double C_PayCreditLag = (double) v_PayCreditLag;

		CCString C_Stub;
		if(VARIANT2CCString (*pStub, C_Stub) != S_OK)
			return S_FALSE;
		
		int l_Stub = ARM_ConvStubRule(C_Stub);

		CCString Res;
		CCString curClass = LOCAL_FTD_CLASS;

		if(ICMLOCAL_FTD	(C_Spread,
						 C_EffectiveDate,
						 C_EndDate,
						 C_First_period_refdate,
						 pFstCpnEffDate,
						 V_labels,
						 l_FixingFreq,
						 l_DayCountFrq,
						 C_IssuerNotional,
						 l_AccruedOnDefault,
						 C_Currency,
						 C_PayCreditLag,
						 l_Stub,
						 l_intRule,
						 l_startAdj,
						 C_result) == ARM_OK)
		{
			Res = LocalMakeObjectId (C_result.getLong (), curClass);
		}
		else
		{
		ERROR_MSG(C_result.getMsg(),pRet);
		return S_OK;
		}

		_variant_t wrap_Res;
		wrap_Res.Attach(*pRet);
		wrap_Res.SetString(Res);
		*pRet=wrap_Res.Detach();

		pRet->vt=VT_BSTR;

		return S_OK;
	}
	catch(Exception&e)
	{
		return createErrorInfo("",e.GetErrorString()); 
	}
	catch(std::exception&e)
	{
		return createErrorInfo("",e.what()); 
	}
	catch(...) 
	{ 
		return createErrorInfo("","Unknown exception"); 
	}
	return E_FAIL;
}


STDMETHODIMP ActiveXModule::ARM_Credit_NTD(double pEffectiveDate,
									   double pEndDate, 
									   double pSpread, 
									   int	  pFirstNumDefault, 	
									   int	  pLastNumDefault, 	
									   VARIANT *pLabels, 	
									   BSTR		pFixingFreq, 
									   BSTR		pDayCountFrq, 
									   double	pFirst_period_refdate, 
									   double	pIssuerNotional, 
									   BSTR	    pAccruedOnDefault, 	
									   BSTR		pCurrency,
									   double	pPayCreditLag,	
									   BSTR		pStub,
									   BSTR		pFrequencyDefLeg,
									   double	pBinary,
									   BSTR		pPayCal,	
									   BSTR		pRcvFeeLeg,
									   double   TradedNotional,
									   BSTR		InclMaturity,
									   double	pFstCpnEffDate,
									   BSTR		intRule,
									   BSTR		startAdj,
									   VARIANT *pRet)
{
	// TODO: Add your implementation code here
	try
	{
		_bstr_t bRcvFeeLeg(pRcvFeeLeg);
		CCString C_RcvFeeLeg = bRcvFeeLeg;
		double d_RcvFeeLeg = 1.;
		if (C_RcvFeeLeg == "L") d_RcvFeeLeg = 1.; else d_RcvFeeLeg = -1.;

		_bstr_t bInclMaturity(InclMaturity);
		CCString C_InclMaturity = bInclMaturity;
		bool d_InclMaturity = false;
		if (C_InclMaturity == "Y") d_InclMaturity = true;

		_bstr_t b_Adjusted(intRule);
		CCString C_Adjusted = b_Adjusted;
		int l_Adjusted = ARM_ConvIntRule (C_Adjusted);

		_bstr_t b_StartAdj(startAdj);
		CCString C_StartAdj = b_StartAdj;
		int l_StartAdj = ARM_ConvStartAdjRule ((const string) C_StartAdj.c_str());

		ARM_result C_result;

		long size = 0;
		VECTOR<CCString> V_labels;

		_bstr_t bPayCal(pPayCal);
		CCString C_PayCal = bPayCal;

		_bstr_t bFixingFreq(pFixingFreq);
		CCString C_FixingFreq = bFixingFreq;

		_bstr_t bFixingFreq_defleg(pFrequencyDefLeg);
		CCString C_FixingFreq_defleg = bFixingFreq_defleg;

		int l_FixingFreq;
		int l_FixingFreq_defleg;

		if((l_FixingFreq = ARM_ConvFrequency (C_FixingFreq, C_result)) == ARM_DEFAULT_ERR)
		{
			ERROR_MSG("Invalid FixingFreq",pRet,ARM_ERROR_FREQ);
			return S_OK;
		}

		if (C_FixingFreq_defleg == "NONE") l_FixingFreq_defleg = l_FixingFreq;
		else
		if((l_FixingFreq_defleg = ARM_ConvFrequency (C_FixingFreq_defleg, C_result)) == ARM_DEFAULT_ERR)
		{
			ERROR_MSG("Invalid FixingFreq Def Leg",pRet,ARM_ERROR_FREQ);
			return S_OK;
		}

		_bstr_t bDayCountFrq(pDayCountFrq);
		CCString C_DayCountFrq = bDayCountFrq;

		long l_DayCountFrq = ARM_ConvDayCount (C_DayCountFrq);

		double C_TradedNotional = 0.;
		if (!TradedNotional) 
			C_TradedNotional = pIssuerNotional;
		else
			C_TradedNotional = TradedNotional;

		// _bstr_t bAccruedOnDefault(pAccruedOnDefault);
		// CCString C_AccruedOnDefault = bAccruedOnDefault;
		// int l_AccruedOnDefault;
		// 
		// if((l_AccruedOnDefault = ARM_ConvAccOnDef (C_AccruedOnDefault, C_result)) == ARM_DEFAULT_ERR)
		// {
		// 	ERROR_MSG("Invalid AccruedOnDefault",pRet,ARM_ERROR_ACCONDEF);
		// 	return S_OK;
		// }

		qPAYMENT_PREMIUM_LEG l_AccruedOnDefault; 
		std::string strAccOnDef; 
		VariantTools::convert(pAccruedOnDefault,strAccOnDef); 
		ICM_EnumsCnv::cnv(strAccOnDef,l_AccruedOnDefault); 

		_bstr_t bCurrency(pCurrency);
		CCString C_Currency = bCurrency;

		if(VARIANT2VECTORCCSTRING(*pLabels, V_labels,size) != S_OK)
			return S_FALSE;

		_bstr_t bStub(pStub);
		CCString C_Stub = bStub;

		int l_Stub = ARM_ConvStubRule(C_Stub);

		CCString Res;
		CCString curClass = LOCAL_NTHTD_CLASS;


		if(ICMLOCAL_NTHTD(pSpread,
						 pEffectiveDate,
						 pEndDate,
						 pFirst_period_refdate,
						 pFstCpnEffDate,
						 pFirstNumDefault,
						 pLastNumDefault,
						 V_labels,
						 l_FixingFreq,
						 l_DayCountFrq,
						 pIssuerNotional,
						 l_AccruedOnDefault,
						 C_Currency,
						 pPayCreditLag,
						 l_Stub,
						 l_FixingFreq_defleg,
						 pBinary,
						 C_PayCal,
						 d_RcvFeeLeg,
						 C_TradedNotional,
						 d_InclMaturity,
						 l_Adjusted,
						 l_StartAdj,
						 C_result) == ARM_OK)
		{
			Res = LocalMakeObjectId (C_result.getLong (), curClass);
		}
		else
		{
		ERROR_MSG(C_result.getMsg(),pRet);
		return S_OK;
		}

		_variant_t wrap_Res;
		wrap_Res.Attach(*pRet);
		wrap_Res.SetString(Res);
		*pRet=wrap_Res.Detach();

		pRet->vt=VT_BSTR;

		return S_OK;
	}
	catch(Exception&e)
	{
		return createErrorInfo("",e.GetErrorString()); 
	}
	catch(std::exception&e)
	{
		return createErrorInfo("",e.what()); 
	}
	catch(...) 
	{ 
		return createErrorInfo("","Unknown exception"); 
	}
	return E_FAIL;
}



STDMETHODIMP ActiveXModule::ARM_Credit_Mezzanine(double pEffectiveDate,
											 double pEndDate, 
											 double pSpread, 
											 double pMezzAmount, 
											 double pSubAmount,
											 VARIANT *pLabels, 
											 VARIANT *pNotionals, 
											 BSTR	pFreqFeeLeg, 
											 BSTR	pDayCountFrq, 
											 double pFirst_period_refdate, 
											 BSTR	pAccruedOnDefault, 	
											 BSTR	pCurrency,
											 double pPayCreditLag,	
											 BSTR	pStub,	
											 BSTR	pFreqDefLeg, 
											 double	pBinary, 
											 BSTR	pPayCal,	
											 BSTR	pRcvFeeLeg,
											 double TradedNotional,
											 BSTR	InclMaturity,
											 double	pFstCpnEffDate,
											 BSTR intRule,
											 BSTR adjStartDate,
											 VARIANT *pRet)
{
	// TODO: Add your implementation code here
	try
	{
		_bstr_t bRcvFeeLeg(pRcvFeeLeg);
		CCString C_RcvFeeLeg = bRcvFeeLeg;
		double d_RcvFeeLeg = 1.;
		if (C_RcvFeeLeg == "L") d_RcvFeeLeg = 1.; else d_RcvFeeLeg = -1.;

		_bstr_t b_Adjusted(intRule);
		CCString C_Adjusted = b_Adjusted;
		int l_Adjusted = ARM_ConvIntRule (C_Adjusted);

		_bstr_t b_StartAdj(adjStartDate);
		CCString C_StartAdj = b_StartAdj;
		int l_StartAdj = ARM_ConvStartAdjRule ((const string) C_StartAdj.c_str());

		_bstr_t bInclMaturity(InclMaturity);
		CCString C_InclMaturity = bInclMaturity;
		bool d_InclMaturity = true;
		if (C_InclMaturity == "N") d_InclMaturity = false;

		ARM_result C_result;

		VECTOR<CCString> V_labels;
		VECTOR<double> V_amounts;

		_bstr_t bFreqFeeLeg(pFreqFeeLeg);
		CCString C_FreqFeeLeg = bFreqFeeLeg;

		int l_FreqFeeLeg;
		if((l_FreqFeeLeg = ARM_ConvFrequency (C_FreqFeeLeg, C_result)) == ARM_DEFAULT_ERR)
		{
			ERROR_MSG("Invalid FreqFeeLeg",pRet,ARM_ERROR_FREQ);
			return S_OK;
		}

		_bstr_t bPayCal(pPayCal);
		CCString C_PayCal = bPayCal;

		_bstr_t bFreqDefLeg(pFreqDefLeg);
		CCString C_FreqDefLeg = bFreqDefLeg;

		int l_FreqDefLeg;
		if((l_FreqDefLeg = ARM_ConvFrequency (C_FreqDefLeg, C_result)) == ARM_DEFAULT_ERR)
		{
			ERROR_MSG("Invalid FreqDefLeg",pRet,ARM_ERROR_FREQ);
			return S_OK;
		}

		_bstr_t bDayCountFrq(pDayCountFrq);
		CCString C_DayCountFrq = bDayCountFrq;

		long l_DayCountFrq = ARM_ConvDayCount (C_DayCountFrq);

		// _bstr_t bAccruedOnDefault(pAccruedOnDefault);
		// CCString C_AccruedOnDefault = bAccruedOnDefault;

		// int l_AccruedOnDefault;
		// if((l_AccruedOnDefault = ARM_ConvAccOnDef (C_AccruedOnDefault, C_result)) == ARM_DEFAULT_ERR)
		// {
		// 	ERROR_MSG("Invalid AccruedOnDefault",pRet,ARM_ERROR_ACCONDEF);
		// 	return S_OK;
		// }

		qPAYMENT_PREMIUM_LEG l_AccruedOnDefault; 
		std::string strAccOnDef; 
		VariantTools::convert(pAccruedOnDefault,strAccOnDef); 
		ICM_EnumsCnv::cnv(strAccOnDef,l_AccruedOnDefault); 

		_bstr_t bCurrency(pCurrency);
		CCString C_Currency = bCurrency;

		long size =0;

		if(VARIANT2VECTORCCSTRING(*pLabels, V_labels,size) != S_OK)
			return S_FALSE;

		if(VARIANT2VECTORDOUBLE (*pNotionals, V_amounts,size) != S_OK)
			return S_FALSE;

		_bstr_t bStub(pStub);
		CCString C_Stub = bStub;

		int l_Stub = ARM_ConvStubRule(C_Stub);

		CCString Res;
		CCString curClass = LOCAL_MEZ_CLASS;

		double C_TradedNotional = 0.;
		if (!TradedNotional) 
			C_TradedNotional = pSubAmount;
		else
			C_TradedNotional = TradedNotional;

		if(ICMLOCAL_MEZZANINE(pSpread,
								   pEffectiveDate,
								   pEndDate,
								   pFirst_period_refdate,
								   pFstCpnEffDate,
								   pMezzAmount,
								   pSubAmount,
								   V_labels,
								   V_amounts,
								   l_FreqFeeLeg,
								   l_FreqDefLeg,
								   l_DayCountFrq,
								   pSubAmount,
								   pSubAmount,
								   l_AccruedOnDefault,
								   C_Currency,
								   pPayCreditLag,
								   l_Stub,
								   pBinary,
								   C_PayCal,
								   d_RcvFeeLeg,
								   C_TradedNotional,
								   1,//qRunning_Leg,
								   6,//qStandart_Recovery_Leg,
								   d_InclMaturity,
								   l_Adjusted,
								   l_StartAdj,
								   C_result) == ARM_OK)
		{
			Res = LocalMakeObjectId (C_result.getLong (), curClass);
		}
		else
		{
		ERROR_MSG(C_result.getMsg(),pRet);
		return S_OK;
		}
		
		_variant_t wrap_Res;
		wrap_Res.Attach(*pRet);
		wrap_Res.SetString(Res);
		*pRet=wrap_Res.Detach();

		pRet->vt=VT_BSTR;

		return S_OK;
}
		catch(Exception&e)
	 {
	  return createErrorInfo("ARM_Credit_Mezzanine",e.GetErrorString()); 
	 }
	 catch(std::exception&e)
	 {
	  return createErrorInfo("ARM_Credit_Mezzanine",e.what()); 
	 }
	 catch(...) 
	 { 
	  return createErrorInfo("ARM_Credit_Mezzanine","Unknown exception"); 
	 }
	 return E_FAIL;
}


STDMETHODIMP ActiveXModule::ARM_Credit_CDO2(double pEffectiveDate,
										double pEndDate,
										BSTR   pPortfolio,
										double pSpread,
										double pMezzAmount,
										double pSubAmount,
										BSTR pFreqFeeLeg,
										BSTR pFreqDefLeg,
										BSTR pDayCountFrq,
										double pFirst_period_refdate,
										BSTR	pAccruedOnDefault,
										BSTR	pCurrency,
										double	pPayCreditLag,
										BSTR	pStub,
										double	pBinary,
										BSTR	pPayCal,	
										BSTR	pRcvFeeLeg,
										double  TradedNotional,
										BSTR	CrossSubordination,
										BSTR	InclMaturity,
										double	pFstCpnEffDate,
										VARIANT *pRet)
{
	// TODO: Add your implementation code here
	
	try
	{

		_bstr_t bRcvFeeLeg(pRcvFeeLeg);
		CCString C_RcvFeeLeg = bRcvFeeLeg;
		double d_RcvFeeLeg = 1.;
		if (C_RcvFeeLeg == "L") d_RcvFeeLeg = 1.; else d_RcvFeeLeg = -1.;

		_bstr_t bInclMaturity(InclMaturity);
		CCString C_InclMaturity = bInclMaturity;
		bool d_InclMaturity = false;
		if (C_InclMaturity == "Y") d_InclMaturity = true;

		_bstr_t bPortfolio(pPortfolio);
		CCString l_Portfolio = bPortfolio;
		long Portfolio = LocalGetNumObjectId (l_Portfolio);
		
		ARM_result C_result;

		_bstr_t bPayCal(pPayCal);
		CCString C_PayCal = bPayCal;

		_bstr_t bCrossSubordination(CrossSubordination);
		CCString C_CrossSubordination = bCrossSubordination;
		bool b_CrossSub = false;
		if (C_CrossSubordination=="Y") b_CrossSub=true;

		_bstr_t bFreqFeeLeg(pFreqFeeLeg);
		CCString C_FreqFeeLeg = bFreqFeeLeg;

		_bstr_t bFreqDefLeg(pFreqDefLeg);
		CCString C_FreqDefLeg = bFreqDefLeg;

		int l_FreqFeeLeg;
		if((l_FreqFeeLeg = ARM_ConvFrequency (C_FreqFeeLeg, C_result)) == ARM_DEFAULT_ERR)
		{
			ERROR_MSG("Invalid FreqFeeLeg",pRet,ARM_ERROR_FREQ);
			return S_OK;
		}

		int l_FreqDefLeg;
		if((l_FreqDefLeg = ARM_ConvFrequency (C_FreqDefLeg, C_result)) == ARM_DEFAULT_ERR)
		{
			ERROR_MSG("Invalid FreqDefLeg",pRet,ARM_ERROR_FREQ);
			return S_OK;
		}

		_bstr_t bDayCountFrq(pDayCountFrq);
		CCString C_DayCountFrq = bDayCountFrq;

		long l_DayCountFrq = ARM_ConvDayCount (C_DayCountFrq);

		/** _bstr_t bAccruedOnDefault(pAccruedOnDefault);
		CCString C_AccruedOnDefault = bAccruedOnDefault;
		int l_AccruedOnDefault;
		if((l_AccruedOnDefault = ARM_ConvAccOnDef (C_AccruedOnDefault, C_result)) == ARM_DEFAULT_ERR)
		{
			ERROR_MSG("Invalid AccruedOnDefault",pRet,ARM_ERROR_ACCONDEF);
			return S_OK;
		} **/ 

		qPAYMENT_PREMIUM_LEG l_AccruedOnDefault; 
		std::string strAccOnDef; 
		VariantTools::convert(pAccruedOnDefault,strAccOnDef); 
		ICM_EnumsCnv::cnv(strAccOnDef,l_AccruedOnDefault); 

		_bstr_t bCurrency(pCurrency);
		CCString C_Currency = bCurrency;

		long size =0;

		_bstr_t bStub(pStub);
		CCString C_Stub = bStub;
		int l_Stub = ARM_ConvStubRule(C_Stub);

		double C_TradedNotional = 0.;
		if (!TradedNotional) 
			C_TradedNotional = pSubAmount;
		else
			C_TradedNotional = TradedNotional;

		CCString Res;
		CCString curClass = LOCAL_MEZ_CLASS;

		if(ICMLOCAL_CDO2(pSpread,
						 pEffectiveDate,
						 pEndDate,
						 pFirst_period_refdate,
						 pFstCpnEffDate,
						 pMezzAmount,
						 pSubAmount,
						 l_FreqFeeLeg,
						 l_DayCountFrq,
						 l_AccruedOnDefault,
						 C_Currency,
						 Portfolio,
						 pPayCreditLag,
						 l_Stub,
						 l_FreqDefLeg,
						 pBinary,
						 C_PayCal,
						 d_RcvFeeLeg,
						 C_TradedNotional,
						 b_CrossSub,
						 d_InclMaturity,
						 C_result) == ARM_OK)
		{
			Res = LocalMakeObjectId (C_result.getLong (), curClass);
		}
		else
		{
			ERROR_MSG(C_result.getMsg(),pRet);
			return S_OK;
		}

		_variant_t wrap_Res;
		wrap_Res.Attach(*pRet);
		wrap_Res.SetString(Res);
		*pRet=wrap_Res.Detach();

		pRet->vt=VT_BSTR;

		return S_OK;
}
	catch(Exception&e)
	{
		return createErrorInfo("",e.GetErrorString()); 
	}
	catch(std::exception&e)
	{
		return createErrorInfo("",e.what()); 
	}
	catch(...) 
	{ 
		return createErrorInfo("","Unknown exception"); 
	}
	return E_FAIL;

}

STDMETHODIMP ActiveXModule:: ARM_Credit_Portfolio(VARIANT *pSecuritiesID, 
											  BSTR pParameters,
										   	  VARIANT *pRet)
{

try
{
	ARM_result C_result;

	VECTOR<CCString> C_SecuritiesID;
	VECTOR<long> l_SecuritiesID;

	_bstr_t b_Parameters(pParameters);
	CCString C_Parameters = b_Parameters;

	long size_SecuritiesID;

	if(VARIANT2VECTORCCSTRING(*pSecuritiesID, C_SecuritiesID,size_SecuritiesID) != S_OK)
		return S_FALSE;

	long ParameterslId = -1;
	if (pParameters) ParameterslId = LocalGetNumObjectId (C_Parameters);

	long security =0;

	for (int i=0; i<size_SecuritiesID; i++)
	{
		security = LocalGetNumObjectId(C_SecuritiesID[i]);

		if (security == ARM_KO)
		{
		ERROR_MSG("Invalid Security",pRet,ARM_ERROR_SECURITY);
		return S_OK;
		}

		l_SecuritiesID.push_back(security);
	}


	CCString Res;
	
	CCString curClass = LOCAL_PRICER_CLASS;

	if (ICMLOCAL_Portfolio(l_SecuritiesID,
						   ParameterslId,	
						   C_result) == ARM_OK)

	{
		Res = LocalMakeObjectId (C_result.getLong (), curClass);
	}
	else
	{
	ERROR_MSG(C_result.getMsg(),pRet);
	return S_OK;
	}

	_variant_t wrap_Res;
	wrap_Res.Attach(*pRet);
	wrap_Res.SetString(Res);
	*pRet=wrap_Res.Detach();

	pRet->vt=VT_BSTR;

	return S_OK;
}
	catch(Exception&e)
	{
		return createErrorInfo("",e.GetErrorString()); 
	}
	catch(std::exception&e)
	{
		return createErrorInfo("",e.what()); 
	}
	catch(...) 
	{ 
		return createErrorInfo("","Unknown exception"); 
	}
	return E_FAIL;
}



STDMETHODIMP ActiveXModule::ARM_Credit_CMTranche(VARIANT *pEffectiveDate,
											 VARIANT *pEndDate, 
											 VARIANT *pSpread, 
											 VARIANT *pMezzAmount, 
											 VARIANT *pSubAmount,
											 VARIANT *pLabels, 
											 VARIANT *pNotionals, 
											 VARIANT *pIndex,
											 VARIANT *pFreqFeeLeg, 
											 VARIANT *pDayCountFrq, 
											 VARIANT *pFirst_period_refdate, 
											 VARIANT *pAccruedOnDefault, 	
											 VARIANT *pCurrency,
											 VARIANT *pPayCreditLag,	
											 VARIANT *pStub,	
											 VARIANT *pFreqDefLeg, 
											 VARIANT *pBinary, 
											 VARIANT *pPayCal,	
											 BSTR	 pRcvFeeLeg,
											 double  TradedNotional,
											 double  FwdFixedDate,
											 BSTR	 InclMaturity,
											 double	 pFstCpnEffDate,
											 VARIANT *pRet)
{
	// TODO: Add your implementation code here
try
{
	_bstr_t bRcvFeeLeg(pRcvFeeLeg);
	CCString C_RcvFeeLeg = bRcvFeeLeg;
	double d_RcvFeeLeg = 1.;
	if (C_RcvFeeLeg == "L") d_RcvFeeLeg = 1.; else d_RcvFeeLeg = -1.;

	_bstr_t bInclMaturity(InclMaturity);
	CCString C_InclMaturity = bInclMaturity;
	bool d_InclMaturity = false;
	if (C_InclMaturity == "Y") d_InclMaturity = true;

	_variant_t v_EffectiveDate (pEffectiveDate);
	v_EffectiveDate.ChangeType( VT_R8 );
	double C_EffectiveDate = (double) v_EffectiveDate;

	_variant_t v_EndDate (pEndDate);
	v_EndDate.ChangeType( VT_R8 );
	double C_EndDate = (double) v_EndDate;

	_variant_t v_Spread (pSpread);
	v_Spread.ChangeType( VT_R8 );
	double C_Spread = (double) v_Spread;

	_variant_t v_MezzAmount (pMezzAmount);
	v_MezzAmount.ChangeType( VT_R8 );
	double C_MezzAmount = (double) v_MezzAmount;

	_variant_t v_SubAmount (pSubAmount);
	v_SubAmount.ChangeType( VT_R8 );
	double C_SubAmount = (double) v_SubAmount;

	_variant_t v_Binary (pBinary);
	v_Binary.ChangeType( VT_R8 );
	double C_Binary = (double) v_Binary;

	ARM_result C_result;

	VECTOR<CCString> V_labels;
	VECTOR<double> V_amounts;

	CCString C_PayCal;
	if(VARIANT2CCString (*pPayCal, C_PayCal) != S_OK)
		return S_FALSE;

	CCString C_FreqFeeLeg;
	if(VARIANT2CCString (*pFreqFeeLeg, C_FreqFeeLeg) != S_OK)
		return S_FALSE;

	CCString C_IndexId;
	if(VARIANT2CCString (*pIndex, C_IndexId) != S_OK)
		return S_FALSE;

	int l_FreqFeeLeg;

	if((l_FreqFeeLeg = ARM_ConvFrequency (C_FreqFeeLeg, C_result)) == ARM_DEFAULT_ERR)
	{
		ERROR_MSG("Invalid FreqFeeLeg",pRet,ARM_ERROR_FREQ);
		return S_OK;
	}

	CCString C_FreqDefLeg;
	if(VARIANT2CCString (*pFreqDefLeg, C_FreqDefLeg) != S_OK)
		return S_FALSE;

	int l_FreqDefLeg;

	if((l_FreqDefLeg = ARM_ConvFrequency (C_FreqDefLeg, C_result)) == ARM_DEFAULT_ERR)
	{
		ERROR_MSG("Invalid FreqDefLeg",pRet,ARM_ERROR_FREQ);
		return S_OK;
	}

	CCString C_DayCountFrq;
	if(VARIANT2CCString (*pDayCountFrq, C_DayCountFrq) != S_OK)
		return S_FALSE;

	long l_DayCountFrq = ARM_ConvDayCount (C_DayCountFrq);

	_variant_t v_First_period_refdate(pFirst_period_refdate);
	v_First_period_refdate.ChangeType( VT_R8 );
	double C_First_period_refdate = (double) v_First_period_refdate;

	/** 
	CCString C_AccruedOnDefault;
	int l_AccruedOnDefault;

	if(VARIANT2CCString (*pAccruedOnDefault, C_AccruedOnDefault) != S_OK)
		return S_FALSE;

		 if((l_AccruedOnDefault = ARM_ConvAccOnDef (C_AccruedOnDefault, C_result)) == ARM_DEFAULT_ERR)
	{
		ERROR_MSG("Invalid AccruedOnDefault",pRet,ARM_ERROR_ACCONDEF);
		return S_OK;
	}
	**/ 

	qPAYMENT_PREMIUM_LEG l_AccruedOnDefault; 
	std::string strAccOnDef; 
	VariantTools::convert(*pAccruedOnDefault,strAccOnDef); 
	ICM_EnumsCnv::cnv(strAccOnDef,l_AccruedOnDefault); 

	CCString C_Currency;

	if(VARIANT2CCString (*pCurrency, C_Currency) != S_OK)
		return S_FALSE;

	long size =0;

	if(VARIANT2VECTORCCSTRING(*pLabels, V_labels,size) != S_OK)
		return S_FALSE;

	if(VARIANT2VECTORDOUBLE (*pNotionals, V_amounts,size) != S_OK)
		return S_FALSE;

	_variant_t v_PayCreditLag(pPayCreditLag);
	v_PayCreditLag.ChangeType( VT_R8 );
	double C_PayCreditLag = (double) v_PayCreditLag;

	CCString C_Stub;
	if(VARIANT2CCString (*pStub, C_Stub) != S_OK)
		return S_FALSE;
	
	int l_Stub = ARM_ConvStubRule(C_Stub);

	CCString Res;
	CCString curClass = LOCAL_MEZ_CLASS;

	double C_TradedNotional = 0.;
	if (!TradedNotional) 
		C_TradedNotional = C_SubAmount;
	else
		C_TradedNotional = TradedNotional;


	if(ICMLOCAL_CMTranche(C_EffectiveDate,
							   C_EndDate,
							   C_First_period_refdate,
							   C_MezzAmount,
							   C_SubAmount,
							   V_labels,
							   V_amounts,
							   LocalGetNumObjectId(C_IndexId),
							   C_Spread,
							   l_FreqFeeLeg,
							   l_FreqDefLeg,
							   l_DayCountFrq,
							   C_SubAmount,
							   C_SubAmount,
							   l_AccruedOnDefault,
							   C_Currency,
							   C_PayCreditLag,
							   l_Stub,
							   C_Binary,
							   C_PayCal,
							   d_RcvFeeLeg,
							   C_TradedNotional,
							   FwdFixedDate,
							   d_InclMaturity,
							   C_result) == ARM_OK)
	{
		Res = LocalMakeObjectId (C_result.getLong (), curClass);
	}
	else
	{
	ERROR_MSG(C_result.getMsg(),pRet);
	return S_OK;
	}

	_variant_t wrap_Res;
	wrap_Res.Attach(*pRet);
	wrap_Res.SetString(Res);
	*pRet=wrap_Res.Detach();

	pRet->vt=VT_BSTR;

	return S_OK;
}
	catch(Exception&e)
	{
		return createErrorInfo("",e.GetErrorString()); 
	}
	catch(std::exception&e)
	{
		return createErrorInfo("",e.what()); 
	}
	catch(...) 
	{ 
		return createErrorInfo("","Unknown exception"); 
	}
	return E_FAIL;
}

STDMETHODIMP ActiveXModule::ARM_Credit_Index(VARIANT *pLabels,
										 double MaturityDate, 
										 double Spread, 
										 BSTR Method, 
										 BSTR Basis,
										 BSTR ResetFreq,
										 BSTR PayFreq,
										 BSTR ccy,
										 BSTR DefCurve, 
										 BSTR fwdRule,
										 BSTR resetTiming,
										 int resetGap,
										 BSTR payTiming,
										 int payGap,
										 BSTR intRule,
										 BSTR AdjCalType,
										 int cm_resetWeekDay,
										 int cm_resetOccur, 
										 BSTR *pRet)
{
	// TODO: Add your implementation code here
try
{
	ARM_result C_result;
	VECTOR<CCString> V_labels;
	long size =0;
	_bstr_t bCcy(ccy);
	string l_ccy = bCcy;
	ARM_result currencyres;
	ARMLOCAL_ARM_GetDefaultCurrency (currencyres);
	if(l_ccy == "DEFAULT")
	{
		if(currencyres.getRetCode () != ARM_OK)
		{
			return S_FALSE;	
		}
		else
		{
			l_ccy = currencyres.getString ();
		}
	}
	_bstr_t b_Basis(Basis);
	CCString C_Basis = b_Basis;
	long l_Basis = ARM_ConvDayCount (C_Basis);

	_bstr_t b_ResetFreq(ResetFreq);
	CCString C_ResetFreq = b_ResetFreq;
	int l_ResetFreq;
	if((l_ResetFreq = ARM_ConvFrequency (C_ResetFreq, C_result)) == ARM_DEFAULT_ERR)
		return S_FALSE;	

	_bstr_t b_PayFreq(PayFreq);
	CCString C_PayFreq = b_PayFreq;
	int l_PayFreq;
	if((l_PayFreq = ARM_ConvFrequency (C_PayFreq, C_result)) == ARM_DEFAULT_ERR)
		return S_FALSE;

	_bstr_t b_Method(Method);
	CCString C_Method = b_Method;
	long l_Method ;

	_bstr_t b_DefCurve(DefCurve);
	CCString C_DefCurve = b_DefCurve;

	_bstr_t b_fwdRule(fwdRule);
	CCString C_fwdRule = b_fwdRule;
	long l_fwdRule;

	_bstr_t b_resetTiming(resetTiming);
	CCString C_resetTiming = b_resetTiming;
	long l_resetTiming;

	_bstr_t b_payTiming(payTiming);
	CCString C_payTiming = b_payTiming;
	long l_payTiming;

	_bstr_t b_intRule(intRule);
	CCString C_intRule = b_intRule;
	long l_intRule;

	_bstr_t b_AdjCalType(AdjCalType);
	CCString C_AdjCalType = b_AdjCalType;
	qCDS_ADJ l_AdjCalType;

	if((l_Method = ARM_ConvIndexMethod(C_Method, C_result)) == ARM_DEFAULT_ERR)
		return S_FALSE;


	if(VARIANT2VECTORCCSTRING(*pLabels, V_labels,size) != S_OK)
		return S_FALSE;

	if((l_fwdRule = ARM_ConvFwdRule (C_fwdRule, C_result)) == ARM_DEFAULT_ERR)
		return S_FALSE;

	l_resetTiming = ARM_ConvPayResetRule (C_resetTiming);
	l_payTiming = ARM_ConvPayResetRule (C_payTiming);
	l_intRule = ARM_ConvIntRule (C_intRule);

	if((l_AdjCalType = ARM_ConvAdjCalCDS (C_AdjCalType, C_result)) == ARM_DEFAULT_ERR)
		return S_FALSE;

	if((l_Method = ARM_ConvIndexMethod(C_Method, C_result)) == ARM_DEFAULT_ERR)
		return S_FALSE;

	vector<string> strVLabels(size);
	for (int i=0; i<size; i++) 
		strVLabels[i] =string( V_labels[i]);
	
	CCString Res;
	CCString curClass = LOCAL_IRINDEX_CLASS;
	long retCode ;
	ARM_Vector vMatDate(1);
	vMatDate.Elt(0) = MaturityDate;
	ARM_Vector vSpread(1);
	vSpread.Elt(0) = Spread;
	retCode =ICMLOCAL_Index(strVLabels[0],strVLabels, 
						l_Basis, 
						l_ResetFreq,
						l_PayFreq, 
						vMatDate, 
						vSpread ,
						l_ccy, 
						l_Method,
						LocalGetNumObjectId(C_DefCurve),  
						l_fwdRule,
						l_resetTiming,
						resetGap,
						l_payTiming,
						payGap,
						l_intRule,
						(qCDS_ADJ)l_AdjCalType,
						cm_resetWeekDay,
						cm_resetOccur) ; 
	
	Res = LocalMakeObjectId (retCode, curClass);
	

	_bstr_t tmpChaine = Res;

	*pRet = SysAllocString(tmpChaine);

	return S_OK;
}
	catch(Exception&e)
	{
		return createErrorInfo("",e.GetErrorString()); 
	}
	catch(std::exception&e)
	{
		return createErrorInfo("",e.what()); 
	}
	catch(...) 
	{ 
		return createErrorInfo("","Unknown exception"); 
	}
	return E_FAIL;
}

STDMETHODIMP ActiveXModule::ARM_Credit_IndexCompo(BSTR IndexName,
										 VARIANT *pLabels,
										 VARIANT* YearFrac, 
										 VARIANT* Spread, 
										 BSTR Method, 
										 BSTR Basis,
										 BSTR ResetFreq,
										 BSTR PayFreq, 
										 BSTR ccy,
										 BSTR fwdRule,
										 BSTR resetTiming,
										 int resetGap,
										 BSTR payTiming,
										 int payGap,
										 BSTR intRule,
										 BSTR AdjCalType,
										 int cm_resetWeekDay,
										 int cm_resetOccur, 
										 BSTR *pRet)
{
	// TODO: Add your implementation code here
try
{
	ARM_result C_result;
	vector<string> v_labels;
	long size =0;
	string sIndexName((_bstr_t)IndexName);
	_bstr_t bCcy(ccy);
	string l_ccy = bCcy;
	ARM_result currencyres;
	ARMLOCAL_ARM_GetDefaultCurrency (currencyres);
	ARM_Vector vYearTerms;
	ARM_Vector vSpreads;
	VariantTools::convert(*YearFrac,vYearTerms);
	VariantTools::convert(*Spread,vSpreads);
	ARM_Date aDate;
	// YearTerm or Date ?!!
		try {
			for (int j = 0; j < vYearTerms.size(); j++)
			{
				ARM_Date aDate;
				Local_XLDATE2ARMDATE(vYearTerms.Elt(j), aDate);
				vYearTerms.Elt(j) = aDate.GetJulian();
			}
		}catch(...) {
			// nothing
		}
	if(l_ccy == "DEFAULT")
	{
		if(currencyres.getRetCode () != ARM_OK)
		{
			return S_FALSE;	
		}
		else
		{
			l_ccy = currencyres.getString ();
		}
	}
	_bstr_t b_Basis(Basis);
	CCString C_Basis = b_Basis;
	long l_Basis = ARM_ConvDayCount (C_Basis);

	_bstr_t b_ResetFreq(ResetFreq);
	CCString C_ResetFreq = b_ResetFreq;
	int l_ResetFreq;
	if((l_ResetFreq = ARM_ConvFrequency (C_ResetFreq, C_result)) == ARM_DEFAULT_ERR)
		return S_FALSE;	

	_bstr_t b_PayFreq(PayFreq);
	CCString C_PayFreq = b_PayFreq;
	int l_PayFreq;
	if((l_PayFreq = ARM_ConvFrequency (C_PayFreq, C_result)) == ARM_DEFAULT_ERR)
		return S_FALSE;

	_bstr_t b_Method(Method);
	CCString C_Method = b_Method;
	long l_Method ;

	_bstr_t b_fwdRule(fwdRule);
	CCString C_fwdRule = b_fwdRule;
	long l_fwdRule;

	_bstr_t b_resetTiming(resetTiming);
	CCString C_resetTiming = b_resetTiming;
	long l_resetTiming;

	_bstr_t b_payTiming(payTiming);
	CCString C_payTiming = b_payTiming;
	long l_payTiming;

	_bstr_t b_intRule(intRule);
	CCString C_intRule = b_intRule;
	long l_intRule;

	CCString C_DefCurve("");

	_bstr_t b_AdjCalType(AdjCalType);
	CCString C_AdjCalType = b_AdjCalType;
	qCDS_ADJ l_AdjCalType;

	if((l_Method = ARM_ConvIndexMethod(C_Method, C_result)) == ARM_DEFAULT_ERR)
		return S_FALSE;

	VariantTools::convert(*pLabels, v_labels);

	if((l_fwdRule = ARM_ConvFwdRule (C_fwdRule, C_result)) == ARM_DEFAULT_ERR)
		return S_FALSE;

	l_resetTiming = ARM_ConvPayResetRule (C_resetTiming);
	l_payTiming = ARM_ConvPayResetRule (C_payTiming);
	l_intRule = ARM_ConvIntRule (C_intRule);

	if((l_AdjCalType = ARM_ConvAdjCalCDS (C_AdjCalType, C_result)) == ARM_DEFAULT_ERR)
		return S_FALSE;

	if((l_Method = ARM_ConvIndexMethod(C_Method, C_result)) == ARM_DEFAULT_ERR)
		return S_FALSE;
	
	CCString Res;
	CCString curClass = LOCAL_IRINDEX_CLASS;
	long retCode ;
	retCode =ICMLOCAL_Index(sIndexName,
						v_labels, 
						l_Basis, 
						l_ResetFreq,
						l_PayFreq, 
						vYearTerms, 
						vSpreads ,
						l_ccy, 
						l_Method,
						LocalGetNumObjectId(C_DefCurve), 
						l_fwdRule,
						l_resetTiming,
						resetGap,
						l_payTiming,
						payGap,
						l_intRule,
						l_AdjCalType,
						cm_resetWeekDay,
						cm_resetOccur) ; 
	
	Res = LocalMakeObjectId (retCode, curClass);
	
	_bstr_t tmpChaine = Res;

	*pRet = SysAllocString(tmpChaine);

	return S_OK;
	}
	catch(Exception&e)
	{
		return createErrorInfo("",e.GetErrorString()); 
	}
	catch(std::exception&e)
	{
		return createErrorInfo("",e.what()); 
	}
	catch(...) 
	{ 
		return createErrorInfo("","Unknown exception"); 
	}
	return E_FAIL;
}
STDMETHODIMP ActiveXModule::ARM_Credit_CDSIndex(VARIANT *pEffectiveDate,
											VARIANT *pEndDate, 
											VARIANT *pSpread,
											VARIANT *pIndex,
											VARIANT *pFixingFreq, 
											VARIANT *pDayCountFrq, 
											VARIANT *pFirst_period_refdate, 
											VARIANT *pFixedPayerAmount, 
											VARIANT *pFloatingPayerAmount, 
											BSTR	StubRule,
											VARIANT *pCurrency,
											BSTR	Adjusted,
											int		pCreditLag,											
											BSTR	pIncludeMaturity,
											VARIANT *pProtectionStartDate,
											VARIANT *pProtectionEndDate, 
											VARIANT *pRet)
{
	// TODO: Add your implementation code here
try
{
	_variant_t v_EffectiveDate (pEffectiveDate);
	v_EffectiveDate.ChangeType( VT_R8 );
	double C_EffectiveDate = (double) v_EffectiveDate;

	_variant_t v_EndDate (pEndDate);
	v_EndDate.ChangeType( VT_R8 );
	double C_EndDate = (double) v_EndDate;

	_variant_t v_ProtectionStartDate (pProtectionStartDate);
	v_ProtectionStartDate.ChangeType( VT_R8 );
	double C_ProtectionStartDate = (double) v_ProtectionStartDate;

	_variant_t v_ProtectionEndDate (pProtectionEndDate);
	v_ProtectionEndDate.ChangeType( VT_R8 );
	double C_ProtectionEndDate = (double) v_ProtectionEndDate;

	_variant_t v_Spread (pSpread);
	v_Spread.ChangeType( VT_R8 );
	double C_Spread = (double) v_Spread;

	ARM_result C_result;

	VECTOR<CCString> V_labels;
	long size =0;

	_bstr_t b_Adjusted(Adjusted);
	CCString C_Adjusted = b_Adjusted;
	int l_Adjusted = ARM_ConvIntRule (C_Adjusted);

	_bstr_t b_StubRule(StubRule);
	CCString C_StubRule = b_StubRule;
	int l_StubRule = ARM_ConvStubRule(C_StubRule);

	_bstr_t b_IncludeMaturity(pIncludeMaturity);
	CCString C_IncludeMaturity = b_IncludeMaturity;
	bool B_IncludeMaturity = false;
	if (C_IncludeMaturity == "Y") B_IncludeMaturity = true;

	CCString C_FixingFreq;
	if(VARIANT2CCString (*pFixingFreq, C_FixingFreq) != S_OK)
		return S_FALSE;
	
	int l_FixingFreq;

	if((l_FixingFreq = ARM_ConvFrequency (C_FixingFreq, C_result)) == ARM_DEFAULT_ERR)
	{
		ERROR_MSG("Invalid FixingFreq",pRet,ARM_ERROR_FREQ);
		return S_OK;
	}

	CCString C_IndexId;
	if(VARIANT2CCString (*pIndex, C_IndexId) != S_OK)
		return S_FALSE;


	CCString C_DayCountFrq;
	if(VARIANT2CCString (*pDayCountFrq, C_DayCountFrq) != S_OK)
		return S_FALSE;

	long l_DayCountFrq = ARM_ConvDayCount (C_DayCountFrq);

	_variant_t v_First_period_refdate(pFirst_period_refdate);
	v_First_period_refdate.ChangeType( VT_R8 );
	double C_First_period_refdate = (double) v_First_period_refdate;

	_variant_t v_FixedPayerAmount(pFixedPayerAmount);
	v_FixedPayerAmount.ChangeType( VT_R8 );
	double C_FixedPayerAmount = (double) v_FixedPayerAmount;

	_variant_t v_FloatingPayerAmount(pFloatingPayerAmount);
	v_FloatingPayerAmount.ChangeType( VT_R8 );
	double C_FloatingPayerAmount = (double) v_FloatingPayerAmount;

	CCString C_Currency;

	if(VARIANT2CCString (*pCurrency, C_Currency) != S_OK)
		return S_FALSE;

	CCString Res;
	CCString curClass = LOCAL_CDSINDEX_CLASS;


	if(ICMLOCAL_CDSIndex(C_Spread,							
						 LocalGetNumObjectId(C_IndexId),
						 C_EffectiveDate,
						 C_EndDate,
						 C_First_period_refdate,
						 l_FixingFreq,
						 l_DayCountFrq,
						 C_FixedPayerAmount,
						 C_FloatingPayerAmount,
						 l_StubRule,
						 C_Currency,
						 l_Adjusted,
						 pCreditLag,						 
						 B_IncludeMaturity,
						 C_ProtectionStartDate,
						 C_ProtectionEndDate,
						 C_result) == ARM_OK)
	{
		Res = LocalMakeObjectId (C_result.getLong (), curClass);
	}
	else
	{
	ERROR_MSG(C_result.getMsg(),pRet);
	return S_OK;
	}

	_variant_t wrap_Res;
	wrap_Res.Attach(*pRet);
	wrap_Res.SetString(Res);
	*pRet=wrap_Res.Detach();

	pRet->vt=VT_BSTR;

	return S_OK;
}
	catch(Exception&e)
	{
		return createErrorInfo("",e.GetErrorString()); 
	}
	catch(std::exception&e)
	{
		return createErrorInfo("",e.what()); 
	}
	catch(...) 
	{ 
		return createErrorInfo("","Unknown exception"); 
	}
	return E_FAIL;
}

STDMETHODIMP ActiveXModule::ARM_Credit_Option(VARIANT * pUnderlyingMaturity,
												double OptionExpiry,
												BSTR Currency,
												BSTR CdsAdj,
												BSTR EndAdj,
												double pStrike,
											    VARIANT *pOptionType,
												VARIANT *pKoType,
												double Notional,
												VARIANT *pRet )
{

	// TODO: Add your implementation code here

	try
	{
		ARM_result C_result;

		string UnderMaturityTenor="";
		string ccy ="";
		string sCDSAdj ="";
		qCDS_ADJ qCdsAdj=(qCDS_ADJ)0;
		string sEndAdj="";
		bool boolEndAdj=false;
		CCString C_OptionType;
		long optiontype;
		CCString C_KoType;
		long kotype;


		double xlDate =0.;
		ARM_Date UnderlyingMaturityDate ;
		try {
			 VariantTools::convert(*pUnderlyingMaturity,xlDate); // try double conversion for date
			 VariantTools::convert(*pUnderlyingMaturity,UnderlyingMaturityDate);
			 UnderMaturityTenor = "";
		}catch (...){
			VariantTools::convert(pUnderlyingMaturity,UnderMaturityTenor);
		}
		VariantTools::convert(Currency, ccy);
		VariantTools::convert(CdsAdj,sCDSAdj);
		VariantTools::convert(EndAdj,sEndAdj);
		
		qCdsAdj = ARM_ConvAdjCalCDS(CCString(sCDSAdj.c_str()),C_result);
		if(sEndAdj=="Y") boolEndAdj = true;
		
		if(VARIANT2CCString (*pOptionType,C_OptionType) != S_OK)
			return S_FALSE;
		if((optiontype = ARM_ConvCallOrPut (C_OptionType, C_result)) == ARM_DEFAULT_ERR)
			return S_FALSE;
		if(VARIANT2CCString (*pKoType,C_KoType) != S_OK)
			return S_FALSE;

		if((kotype = ARM_ConvKO_Or_NoKO (C_KoType, C_result)) == ARM_DEFAULT_ERR)
			return S_FALSE;
		
		long newId = ICMLOCAL_Option(UnderMaturityTenor,UnderlyingMaturityDate, OptionExpiry, ccy, qCdsAdj,boolEndAdj,
										pStrike,optiontype,(qDEF_MAT)kotype,Notional, C_result);
		string newLabel= LocalPersistent::get().getStringId(newId, LOCAL_OPTION_CLASS);
		VariantTools::convert(newLabel, *pRet); 
		return S_OK ;
	} 
	catch(Exception&e)
	 {
	  return createErrorInfo("",e.GetErrorString()); 
	 }
	 catch(std::exception&e)
	 {
	  return createErrorInfo("",e.what()); 
	 }
	 catch(...) 
	 { 
	  return createErrorInfo("","Unknown exception"); 
	 }
	 return E_FAIL;
}
	

STDMETHODIMP ActiveXModule::ARM_GETINSTRUMENTFROMCALYPSO(BSTR CalypsoId,
													BSTR Type,
													double AsOf,
													BSTR ModelType,
													VARIANT* pRet)
{
	try
{

	ARM_result C_result;

	_bstr_t bCalypsoId(CalypsoId);
	CCString l_CalypsoId = bCalypsoId;

	_bstr_t bType(Type);
	CCString l_Type = bType;

	_bstr_t bModelType(ModelType);
	CCString l_ModelType = bModelType;

	CCString Res;
	CCString curClass = LOCAL_CDS_CLASS;


	if(ARMLOCAL_GETINSTRUMENTFROMCALYPSO (l_CalypsoId,
										 l_Type,
										 l_ModelType,
										 AsOf,
										 C_result) == ARM_OK)
	{
		Res = LocalMakeObjectId (C_result.getLong (), curClass);
	}
	else
	{
	ERROR_MSG(C_result.getMsg(),pRet);
	return S_OK;
	}

	_variant_t wrap_Res;
	wrap_Res.Attach(*pRet);
	wrap_Res.SetString(Res);
	*pRet=wrap_Res.Detach();

	pRet->vt=VT_BSTR;

	return S_OK;
}
	catch(Exception&e)
	{
		return createErrorInfo("",e.GetErrorString()); 
	}
	catch(std::exception&e)
	{
		return createErrorInfo("",e.what()); 
	}
	catch(...) 
	{ 
		return createErrorInfo("","Unknown exception"); 
	}
	return E_FAIL;

}


STDMETHODIMP ActiveXModule::ARM_GETINSTRUMENTFROMSUMMIT(BSTR SummitId,
													BSTR Type,
													double AsOf,
													BSTR Exoticfilter,
													VARIANT* pRet)
{
	// TODO: Add your implementation code here
try
{

	ARM_result C_result;

	_bstr_t bSummitId(SummitId);
	CCString l_SummitId = bSummitId;

	_bstr_t bType(Type);
	CCString l_Type = bType;

	_bstr_t bExoticfilter(Exoticfilter);
	CCString l_Exoticfilter = bExoticfilter;

	CCString Res;
	CCString curClass = LOCAL_CDS_CLASS;


	if(ARMLOCAL_GETINSTRUMENTFROMSUMMIT (l_SummitId,
										 l_Type,
										 AsOf,
										 l_Exoticfilter,
										 C_result) == ARM_OK)
	{
		Res = LocalMakeObjectId (C_result.getLong (), curClass);
	}
	else
	{
	ERROR_MSG(C_result.getMsg(),pRet);
	return S_OK;
	}

	_variant_t wrap_Res;
	wrap_Res.Attach(*pRet);
	wrap_Res.SetString(Res);
	*pRet=wrap_Res.Detach();

	pRet->vt=VT_BSTR;

	return S_OK;
}
	catch(Exception&e)
	{
		return createErrorInfo("",e.GetErrorString()); 
	}
	catch(std::exception&e)
	{
		return createErrorInfo("",e.what()); 
	}
	catch(...) 
	{ 
		return createErrorInfo("","Unknown exception"); 
	}
	return E_FAIL;
}


STDMETHODIMP ActiveXModule::ARM_Credit_EmptyLeg(VARIANT *pRet)
{
	// TODO: Add your implementation code here
try
{

	ARM_result C_result;

	long size =0;

	CCString Res;
	CCString curClass = LOCAL_CDS_CLASS;


	if(ICMLOCAL_GenLeg(47000.,47365.,0.,0.,-1,-1,-1,0,1,4,1,1,-1,"EUR",47365.,
					  true,1,0,-1,0,-999,"UNDEF",0,qACCRUED_SETTLED,C_result) == ARM_OK)
	{
		Res = LocalMakeObjectId (C_result.getLong (), curClass);
	}
	else
	{
	ERROR_MSG(C_result.getMsg(),pRet);
	return S_OK;
	}

	_variant_t wrap_Res;
	wrap_Res.Attach(*pRet);
	wrap_Res.SetString(Res);
	*pRet=wrap_Res.Detach();

	pRet->vt=VT_BSTR;

	return S_OK;
}
	catch(Exception&e)
	{
		return createErrorInfo("",e.GetErrorString()); 
	}
	catch(std::exception&e)
	{
		return createErrorInfo("",e.what()); 
	}
	catch(...) 
	{ 
		return createErrorInfo("","Unknown exception"); 
	}
	return E_FAIL;
}


STDMETHODIMP ActiveXModule::ARM_Credit_IRLEGTOCREDITLEG(BSTR SwapLegId,
													BSTR LegType,
													BSTR creditindexId,
													BSTR PricerId,
													BSTR *pRet)
{
	// TODO: Add your implementation code here
try
{

	ARM_result C_result;

	_bstr_t bSwapLegId(SwapLegId);
	CCString l_SwapLegId = bSwapLegId;

	_bstr_t bcreditindexId(creditindexId);
	CCString l_creditindexId = bcreditindexId;

	_bstr_t bLegType(LegType);
	CCString C_LegType = bLegType;

	_bstr_t bPricerId(PricerId);
	CCString l_PricerId = bPricerId;

	long size =0;

	CCString Res;
	CCString curClass = LOCAL_CDS_CLASS;

	int l_LegType = ARM_ConvLegType(C_LegType, C_result);

	if(ICMLOCAL_IRLEGTOCREDITLEG(LocalGetNumObjectId(l_SwapLegId),
							     l_LegType,
							     LocalGetNumObjectId(l_creditindexId),
							     LocalGetNumObjectId(l_PricerId),C_result) == ARM_OK)
	{
		Res = LocalMakeObjectId (C_result.getLong (), curClass);
	}

	_bstr_t tmpChaine = Res;
	*pRet = SysAllocString(tmpChaine);

	return S_OK;
}
	catch(Exception&e)
	{
		return createErrorInfo("",e.GetErrorString()); 
	}
	catch(std::exception&e)
	{
		return createErrorInfo("",e.what()); 
	}
	catch(...) 
	{ 
		return createErrorInfo("","Unknown exception"); 
	}
	return E_FAIL;
}



STDMETHODIMP ActiveXModule::ARM_Credit_Collateral(VARIANT *pLabels, 
											 VARIANT *pNotionals, 
											 VARIANT *pRet)
{
	// TODO: Add your implementation code here
	try
	{
		std::vector<string> Labels; 
		ARM_Vector notionals; 
		VariantTools::convert(*pLabels,Labels); 
		VariantTools::convert(*pNotionals,notionals);
		long objId = ICMLOCAL_Collateral(Labels,notionals); 
		std::string name = LocalPersistent::get().getStringId(objId,LOCAL_CREDIT_PORTFOLIO_CLASS); 
		VariantTools::convert(name,*pRet); 
		return S_OK;
	}
	catch(Exception&e)
	{
		return createErrorInfo("",e.GetErrorString()); 
	}
	catch(std::exception&e)
	{
		return createErrorInfo("",e.what()); 
	}
	catch(...) 
	{ 
		return createErrorInfo("","Unknown exception"); 
	}
	return E_FAIL;
}
STDMETHODIMP ActiveXModule::ARM_Credit_VariableCollateral(VARIANT pLabels, 
											 VARIANT pNotionalsIds, 
											 VARIANT *pRet)
{
	// TODO: Add your implementation code here
	try
	{
		std::vector<string> Labels; 
		std::vector<string> notionalsStringIds; 
		VariantTools::convert(pLabels,Labels); 
		VariantTools::convert(pNotionalsIds,notionalsStringIds); 
		ICM_Vector<long> notionalIds(notionalsStringIds.size()); 
		for(unsigned int i=0;i<notionalsStringIds.size();i++) 
			notionalIds[i]= LocalPersistent::get().getObjectId(notionalsStringIds[i]); 
		long objId = ICMLOCAL_VariableCollateral(Labels,notionalIds); 
		std::string name = LocalPersistent::get().getStringId(objId,LOCAL_CREDIT_PORTFOLIO_CLASS); 
		VariantTools::convert(name,*pRet); 
		return S_OK;
	}
	catch(Exception&e)
	{
		return createErrorInfo("",e.GetErrorString()); 
	}
	catch(std::exception&e)
	{
		return createErrorInfo("",e.what()); 
	}
	catch(...) 
	{ 
		return createErrorInfo("","Unknown exception"); 
	}
	return E_FAIL;
}


STDMETHODIMP ActiveXModule::ARM_Credit_CDSGEN(BSTR FeeLegId,
										  BSTR DefLegId,
										  double RcvFee,
										  double TradedNot,
										  BSTR *pRet)
{
	// TODO: Add your implementation code here
try
{

	ARM_result C_result;

	_bstr_t bFeeLegId(FeeLegId);
	CCString l_FeeLegId = bFeeLegId;

	_bstr_t bDefLegId(DefLegId);
	CCString l_DefLegId = bDefLegId;

	long size =0;

	CCString Res;
	CCString curClass = LOCAL_CDS_CLASS;

	if(ICMLOCAL_CDSGEN(LocalGetNumObjectId(l_FeeLegId),
					   LocalGetNumObjectId(l_DefLegId),
					   RcvFee,TradedNot,C_result) == ARM_OK)
	{
		Res = LocalMakeObjectId (C_result.getLong (), curClass);
	}

	_bstr_t tmpChaine = Res;
	*pRet = SysAllocString(tmpChaine);

	return S_OK;
}
	catch(Exception&e)
	{
		return createErrorInfo("",e.GetErrorString()); 
	}
	catch(std::exception&e)
	{
		return createErrorInfo("",e.what()); 
	}
	catch(...) 
	{ 
		return createErrorInfo("","Unknown exception"); 
	}
	return E_FAIL;
}

STDMETHODIMP ActiveXModule::ARM_Credit_NTDGEN(BSTR CdsId,
										int firstnumdef,
										int lastnumdef,
										BSTR CollateralId,
										double binary,
										double rcvfee,
										BSTR *pRet)
{
	// TODO: Add your implementation code here
try
{

	ARM_result C_result;

	_bstr_t bCdsId(CdsId);
	CCString l_CdsId = bCdsId;

	_bstr_t bCollateralId(CollateralId);
	CCString l_CollateralId = bCollateralId;

	long size =0;

	CCString Res;
	CCString curClass = LOCAL_CDS_CLASS;

	if(ICMLOCAL_NTDGEN(LocalGetNumObjectId(l_CdsId),
					   firstnumdef,
					   lastnumdef,
					   LocalGetNumObjectId(l_CollateralId),
					   binary,rcvfee,C_result) == ARM_OK)
	{
		Res = LocalMakeObjectId (C_result.getLong (), curClass);
	}

	_bstr_t tmpChaine = Res;
	*pRet = SysAllocString(tmpChaine);

	return S_OK;
}
	catch(Exception&e)
	{
		return createErrorInfo("",e.GetErrorString()); 
	}
	catch(std::exception&e)
	{
		return createErrorInfo("",e.what()); 
	}
	catch(...) 
	{ 
		return createErrorInfo("","Unknown exception"); 
	}
	return E_FAIL;
}

STDMETHODIMP ActiveXModule::ARM_Credit_CDOGEN(BSTR CdsId,
										double subamount,
										BSTR CollateralId,
										double binary,
										double rcvfee,
										BSTR *pRet)
{
	// TODO: Add your implementation code here
try
{

	ARM_result C_result;

	_bstr_t bCdsId(CdsId);
	CCString l_CdsId = bCdsId;

	_bstr_t bCollateralId(CollateralId);
	CCString l_CollateralId = bCollateralId;

	long size =0;

	CCString Res;
	CCString curClass = LOCAL_CDS_CLASS;

	if(ICMLOCAL_CDOGEN(LocalGetNumObjectId(l_CdsId),
					   subamount,0,
					   LocalGetNumObjectId(l_CollateralId),
					   binary,rcvfee,C_result) == ARM_OK)
	{
		Res = LocalMakeObjectId (C_result.getLong (), curClass);
	}

	_bstr_t tmpChaine = Res;
	*pRet = SysAllocString(tmpChaine);

	return S_OK;
}
	catch(Exception&e)
	{
		return createErrorInfo("",e.GetErrorString()); 
	}
	catch(std::exception&e)
	{
		return createErrorInfo("",e.what()); 
	}
	catch(...) 
	{ 
		return createErrorInfo("","Unknown exception"); 
	}
	return E_FAIL;
}


STDMETHODIMP ActiveXModule::ARM_Credit_CDO_SQUARE_GEN(BSTR CdsId,
												  double subamount,
												  BSTR PortfolioId,
												  double binary,
												  double rcvfee,
												  BSTR *pRet)
{
	// TODO: Add your implementation code here
try
{

	ARM_result C_result;

	_bstr_t bCdsId(CdsId);
	CCString l_CdsId = bCdsId;

	_bstr_t bPortfolioId(PortfolioId);
	CCString l_PortfolioId = bPortfolioId;

	long size =0;

	CCString Res;
	CCString curClass = LOCAL_CDS_CLASS;

	if(ICMLOCAL_CDO_SQUARE_GEN(LocalGetNumObjectId(l_CdsId),
								subamount,
								LocalGetNumObjectId(l_PortfolioId),
								binary,rcvfee,C_result) == ARM_OK)
	{
		Res = LocalMakeObjectId (C_result.getLong (), curClass);
	}

	_bstr_t tmpChaine = Res;
	*pRet = SysAllocString(tmpChaine);

	return S_OK;
}
	catch(Exception&e)
	{
		return createErrorInfo("",e.GetErrorString()); 
	}
	catch(std::exception&e)
	{
		return createErrorInfo("",e.what()); 
	}
	catch(...) 
	{ 
		return createErrorInfo("","Unknown exception"); 
	}
	return E_FAIL;
}


STDMETHODIMP ActiveXModule::ARM_Credit_GenLeg(double	StartDate,
										  double	EndDate, 
										  double	FixedRate, 
										  double	FixedNotional, 
										  BSTR		RefValNotional, 
										  BSTR		RefValRate, 
										  BSTR		XChangeNotional, 
										  BSTR		Frequency, 
									      BSTR		Basis,
									      BSTR      payTiming,
									      BSTR		intrule,
									      BSTR		stubrule,
									      BSTR	    ccyid,
									      BSTR		paycalname,
									      double	refdate, 
									      BSTR		includematurity,
									      BSTR	    adjstartdate,
										  BSTR	    legtype,
										  BSTR		indexobj,
										  int		creditlag,
										  double    binary,
										  BSTR      Name,
										  BSTR      Nxchange,
										  BSTR		baccruedOnDef,
									      VARIANT	*pRet)
{
	// TODO: Add your implementation code here
try
{

	ARM_result C_result;

	long size =0;

	_bstr_t b_Name(Name);
	CCString C_Name = b_Name;

	_bstr_t b_RefValNotional(RefValNotional);
	CCString C_RefValNotional = b_RefValNotional;

	_bstr_t b_RefValRate(RefValRate);
	CCString C_RefValRate = b_RefValRate;

	_bstr_t b_XChangeNotional(XChangeNotional);
	CCString C_XChangeNotional = b_XChangeNotional;

	_bstr_t b_Frequency(Frequency);
	CCString C_Frequency = b_Frequency;


	/** _bstr_t b_accruedOnDef(baccruedOnDef);
	CCString C_accruedOnDef = b_accruedOnDef;
	long iaccruedOnDef;
	if ((iaccruedOnDef = ARM_ConvAccOnDef(C_accruedOnDef, C_result)) == ARM_DEFAULT_ERR)
	{
		ERROR_MSG("Invalid accruedOnDef",pRet,ARM_ERROR_FREQ);
		return S_OK;
	}
	**/ 

	qPAYMENT_PREMIUM_LEG iaccruedOnDef; 
	std::string strAccOnDef; 
	VariantTools::convert(baccruedOnDef,strAccOnDef); 
	ICM_EnumsCnv::cnv(strAccOnDef,iaccruedOnDef); 

	int l_Frequency;
	if((l_Frequency = ARM_ConvFrequency (C_Frequency, C_result)) == ARM_DEFAULT_ERR)
	{
		ERROR_MSG("Invalid FixingFreq",pRet,ARM_ERROR_FREQ);
		return S_OK;
	}

	_bstr_t b_Basis(Basis);
	CCString C_Basis = b_Basis;
	long l_Basis = ARM_ConvDayCount(C_Basis);

	_bstr_t b_payTiming(payTiming);
	CCString C_payTiming = b_payTiming;
	int l_payTiming;
	
	if((l_payTiming = ARM_ConvPayResetRule (C_payTiming)) == ARM_DEFAULT_ERR)
	{
		ERROR_MSG("Invalid ARM_ConvPayResetRule",pRet,ARM_ERROR_FREQ);
		return S_OK;
	}

	_bstr_t b_intrule(intrule);
	CCString C_intrule = b_intrule;
	int l_intrule;

	if((l_intrule = ARM_ConvIntRule (C_intrule)) == ARM_DEFAULT_ERR)
	{
		ERROR_MSG("Invalid ARM_ConvIntRule",pRet,ARM_ERROR_FREQ);
		return S_OK;
	}

	_bstr_t b_stubrule(stubrule);
	CCString C_stubrule = b_stubrule;
	int l_stubrule;

	if((l_stubrule = ARM_ConvStubRule (C_stubrule)) == ARM_DEFAULT_ERR)
	{
		ERROR_MSG("Invalid ARM_ConvStubRule",pRet,ARM_ERROR_FREQ);
		return S_OK;
	}

	_bstr_t b_ccyid(ccyid);
	CCString C_ccyid = b_ccyid;

	_bstr_t b_paycalname(paycalname);
	CCString C_paycalname = b_paycalname;

	_bstr_t b_IncludeMaturity(includematurity);
	CCString C_IncludeMaturity = b_IncludeMaturity;
	bool B_IncludeMaturity = false;
	if (C_IncludeMaturity == "Y") B_IncludeMaturity = true;

	_bstr_t b_adjstartdate(adjstartdate);
	CCString C_adjstartdate = b_adjstartdate;
	int l_AdjStartDate = 0;
	if (C_adjstartdate == "Y") l_AdjStartDate = 1;

	_bstr_t b_legtype(legtype);
	CCString C_legtype = b_legtype;
	int l_LegType;
	if((l_LegType = ARM_ConvLegType(C_legtype, C_result)) == ARM_DEFAULT_ERR)
	{
		ERROR_MSG("Invalid ARM_ConvLegType",pRet,ARM_ERROR_FREQ);
		return S_OK;
	}

	_bstr_t b_Nxchange(Nxchange);
	CCString C_Nxchange = b_Nxchange;
	int l_Nxchange = ARM_NotionalExchange((const char*)C_Nxchange);

	_bstr_t b_indexobj(indexobj);
	CCString C_indexobj = b_indexobj;

	CCString Res;
	CCString curClass = LOCAL_CDS_CLASS;


	if(ICMLOCAL_GenLeg (StartDate,
										EndDate,
										FixedRate,
										FixedNotional,
										LocalGetNumObjectId(C_RefValNotional),
										LocalGetNumObjectId(C_RefValRate),
										LocalGetNumObjectId(C_XChangeNotional),
										l_Frequency,
										l_Basis,
										l_payTiming,
										l_intrule,
										l_stubrule,
										LocalGetNumObjectId(C_ccyid),
										C_paycalname,
										refdate,
										B_IncludeMaturity,
										l_AdjStartDate,
										l_LegType,
										LocalGetNumObjectId(C_indexobj),
										creditlag,
										binary,
										C_Name,
										l_Nxchange,
										iaccruedOnDef,
										C_result) == ARM_OK)
	{
		Res = LocalMakeObjectId (C_result.getLong (), curClass);
	}
	else
	{
	ERROR_MSG(C_result.getMsg(),pRet);
	return S_OK;
	}

	_variant_t wrap_Res;
	wrap_Res.Attach(*pRet);
	wrap_Res.SetString(Res);
	*pRet=wrap_Res.Detach();

	pRet->vt=VT_BSTR;

	return S_OK;
}

	catch(Exception&e)
	{
		return createErrorInfo("",e.GetErrorString()); 
	}
	catch(std::exception&e)
	{
		return createErrorInfo("",e.what()); 
	}
	catch(...) 
	{ 
		return createErrorInfo("","Unknown exception"); 
	}
	return E_FAIL;
}


//	------------------------------------------------------------------------------------
STDMETHODIMP 
ActiveXModule::ARM_Credit_GetLastFixingDate(
										BSTR instId_,
										VARIANT xlAsofDate,
										VARIANT*pRet)
{
	try 
	{
		ARM_result result ; 
		long instId=LocalGetNumObjectId(CCString(_bstr_t (instId_))); 
		if( ICMLOCAL_getLastFixingDate(instId,_variant_t(xlAsofDate).dblVal,result)==ARM_KO )
			ERROR_MSG("ERROR",pRet);
		else 
		{
			double dateRes = Local_ARMDATE2XLDATE(result.getString ());
			pRet->dblVal=dateRes; 
			pRet->vt=VT_R8;
			return S_OK ;
		}
	}
	catch(Exception&e)
	{
		return createErrorInfo("",e.GetErrorString()); 
	}
	catch(std::exception&e)
	{
		return createErrorInfo("",e.what()); 
	}
	catch(...) 
	{ 
		return createErrorInfo("","Unknown exception"); 
	}
	return E_FAIL;
}

//	------------------------------------------------------------------------------------
STDMETHODIMP 
ActiveXModule::ARM_Credit_SetPastFixing(
									BSTR instId_,
									VARIANT resetDate,
									VARIANT fixingValue,
									VARIANT*pRet)
{
	try 
	{
		ARM_result result ;
		long instId=LocalGetNumObjectId(CCString(_bstr_t (instId_))); 
		if ( ICMLOCAL_setPastFixing(instId,
				_variant_t(resetDate),
				_variant_t(fixingValue),
				result) == ARM_KO)
			ERROR_MSG("ERROR",pRet);
		return S_OK; 
	}
	catch(Exception&e)
	{
		return createErrorInfo("",e.GetErrorString()); 
	}
	catch(std::exception&e)
	{
		return createErrorInfo("",e.what()); 
	}
	catch(...) 
	{ 
		return createErrorInfo("","Unknown exception"); 
	}
	return E_FAIL;
}


STDMETHODIMP 
ActiveXModule::ARM_Credit_GetBounds(
									BSTR instId_,
									double* low,
									double* up)
{
	try 
	{
		ARM_result result ;
		long instId=LocalGetNumObjectId(CCString(_bstr_t (instId_))); 
		if ( ICMLOCAL_GetBounds(instId,
				*low,
				*up,
				result) != ARM_KO)
		 
		return S_OK; 
	 
	}
	catch(Exception&e)
	{
		return createErrorInfo("",e.GetErrorString()); 
	}
	catch(std::exception&e)
	{
		return createErrorInfo("",e.what()); 
	}
	catch(...) 
	{ 
		return createErrorInfo("","Unknown exception"); 
	}
	return E_FAIL;
}


STDMETHODIMP ActiveXModule::ARM_Credit_Customized_CDO(
									VARIANT* pLabels, 
									VARIANT *pNotionals,
									BSTR pCurrency,
									BSTR pDefaultLeg, 
									BSTR pPremiumLeg, 
									BSTR pParameters, 
									BSTR *pRet)
{

	try 
	{
	ARM_result C_result;

	_bstr_t bPremiumLeg(pPremiumLeg);
	CCString l_PremiumLeg = bPremiumLeg;

	_bstr_t bDefaultLeg(pDefaultLeg);
	CCString l_DefaultLeg = bDefaultLeg;

	_bstr_t bParameters(pParameters);
	CCString l_Parameters = bParameters;

	VECTOR<CCString> V_labels;
	VECTOR<double> V_amounts;

	_bstr_t bCurrency(pCurrency);
	CCString C_Currency = bCurrency;

	long size =0;

	if(VARIANT2VECTORCCSTRING(*pLabels, V_labels,size) != S_OK)
		return S_FALSE;

	if(VARIANT2VECTORDOUBLE (*pNotionals, V_amounts,size) != S_OK)
		return S_FALSE;
	
	CCString Res;
	CCString curClass = LOCAL_CUSTOMIZED_CDO_CLASS;

	if(ICMLOCAL_Customized_CDO(V_labels,
						V_amounts,
						C_Currency,
						LocalGetNumObjectId(l_DefaultLeg),
						LocalGetNumObjectId(l_PremiumLeg),
						LocalGetNumObjectId(l_Parameters),
						C_result) == ARM_OK)
	{
		Res = LocalMakeObjectId (C_result.getLong (), curClass);
	}

	_bstr_t tmpChaine = Res;
	*pRet = SysAllocString(tmpChaine);

	return S_OK;

	
	}
	catch(Exception&e)
	{
		return createErrorInfo("",e.GetErrorString()); 
	}
	catch(std::exception&e)
	{
		return createErrorInfo("",e.what()); 
	}
	catch(...) 
	{ 
		return createErrorInfo("","Unknown exception"); 
	}
	return E_FAIL;


}


STDMETHODIMP ActiveXModule::ARM_Credit_CLN(double	pEffectiveDate,
									   double	pEndDate, 
									   double	pSpread, 
									   BSTR		pIndexId, 
									   double	pFirst_period_refdate, 
									   double	pFstCpnEffDate,
									   double	pNotional, 
									   BSTR		AccOnDef,
									   BSTR		pDayCount, 
									   BSTR		pDecompFreq, 
									   BSTR		StubRule,
									   double	resetgap,
									   BSTR     pCurrency,
									   BSTR     ResetCal,
									   BSTR     PayCal,
									   BSTR     Nxchange,
									   BSTR	    pIncludeMaturity,
									   BSTR		AdjustedStartDate,
									   double   Binary,
									   BSTR		*pRet)
{
	// TODO: Add your implementation code here
try
{

	ARM_result C_result;

	long size =0;

	_bstr_t bIndexId(pIndexId);
	CCString C_IndexId = bIndexId;

	// _bstr_t bAccruedOnDefault(AccOnDef);
	// CCString C_AccruedOnDefault = bAccruedOnDefault;
	// int l_AccruedOnDefault;
	// if((l_AccruedOnDefault = ARM_ConvAccOnDef (C_AccruedOnDefault, C_result)) == ARM_DEFAULT_ERR)
	// {return S_OK; }
	qPAYMENT_PREMIUM_LEG l_AccruedOnDefault; 
	std::string strAccOnDef; 
	VariantTools::convert(AccOnDef,strAccOnDef); 
	ICM_EnumsCnv::cnv(strAccOnDef,l_AccruedOnDefault); 

	_bstr_t b_DayCount(pDayCount);
	CCString C_DayCount = b_DayCount;
	long l_DayCountFrq = ARM_ConvDayCount (C_DayCount);

	_bstr_t b_DecompFreq(pDecompFreq);
	CCString C_DecompFreq = b_DecompFreq;
	long l_DecompFreq = ARM_ConvDecompFrequency (C_DecompFreq);

	_bstr_t b_StubRule(StubRule);
	CCString C_StubRule = b_StubRule;
	int l_StubRule = ARM_ConvStubRule(C_StubRule);

	_bstr_t b_Currency(pCurrency);
	CCString C_Currency = b_Currency;

	_bstr_t b_ResetCal(ResetCal);
	CCString C_ResetCal = b_ResetCal;

	_bstr_t b_PayCal(PayCal);
	CCString C_PayCal = b_PayCal;

	_bstr_t b_Nxchange(Nxchange);
	CCString C_Nxchange = b_Nxchange;
	int l_Nxchange = ARM_NotionalExchange(C_Nxchange);

	_bstr_t b_IncludeMaturity(pIncludeMaturity);
	CCString C_IncludeMaturity = b_IncludeMaturity;
	bool B_IncludeMaturity = false;
	if (C_IncludeMaturity == "Y") B_IncludeMaturity = true;

	_bstr_t b_Adjusted(AdjustedStartDate);
	CCString C_Adjusted = b_Adjusted;
	int l_Adjusted = ARM_ConvIntRule (C_Adjusted);

//	_bstr_t b_IssuerName(IssuerName);
//	CCString C_IssuerName = b_IssuerName;

//	int l_FixingFreq;
//	if((l_FixingFreq = ARM_ConvFrequency (C_FixingFreq, C_result)) == ARM_DEFAULT_ERR)
//	{
//		ERROR_MSG("Invalid FixingFreq",pRet,ARM_ERROR_FREQ);
//		return S_OK;
//	}

	long LegType = ARM_ConvLegType ("SWAPLEG",C_result);
	CCString Res;
	CCString curClass = LOCAL_CDS_CLASS;


	if(ICMLOCAL_CLN(pEffectiveDate,
					pEndDate,
					pFirst_period_refdate,
					pFstCpnEffDate,
					LocalGetNumObjectId(C_IndexId),
					pSpread,
					pNotional,
					l_AccruedOnDefault,
					l_DayCountFrq,
					l_DecompFreq,
					l_StubRule,
					resetgap,
					C_Currency,
					C_ResetCal,
					C_PayCal,
					l_Nxchange,
					B_IncludeMaturity,	
					l_Adjusted,
					LegType,
					Binary,
					C_result) == ARM_OK)
	{
		Res = LocalMakeObjectId (C_result.getLong (), curClass);
	}
	else
	{return S_OK;}

	_bstr_t tmpChaine = Res;
	*pRet = SysAllocString(tmpChaine);

	return S_OK;
}
	catch(Exception&e)
	{
		return createErrorInfo("",e.GetErrorString()); 
	}
	catch(std::exception&e)
	{
		return createErrorInfo("",e.what()); 
	}
	catch(...) 
	{ 
		return createErrorInfo("","Unknown exception"); 
	}
	return E_FAIL;

}

STDMETHODIMP ActiveXModule::ARM_Credit_CorridorLeg_Sche(double Notional,
														BSTR RecieveOrPay,
														BSTR RefValueSpreads,
									                    BSTR floatingIdx,
														double leverageFloatIdx,
														BSTR creditIdx,
														BSTR refvalueKUP,
														BSTR refvalueKDW,
														BSTR ScheduleInfoId,
														VARIANT* accondef,
														BSTR disc_ccy,
														BSTR Name,
														VARIANT		*pRet)
{
	ARM_result C_result;

	try
	{
		long  prevId = -1;
	
		string C_RefValueSpreads = "";
		long l_RefValueSpreads = 0;
		string strRecieveOrPay="";
		long ConvertRorP = 1;
		string C_floatingIdx="";
		long l_floatingIdx =0;		

		string C_refvalueKUP ="";
		long l_refvalueKUP=0;

		string C_refvalueKDW="";
		long l_refvalueKDW=0;

		string strSchedule_info = "";
		long l_schedule_Info = 0;
		string C_disc_ccy = "";
		string C_Name = "";
		string C_creditIdx="";
		int l_creditIdx = -1;


		qPAYMENT_PREMIUM_LEG l_AccruedOnDefault; 
		std::string strAccOnDef; 
		VariantTools::convert(*accondef,strAccOnDef); 
		ICM_EnumsCnv::cnv(strAccOnDef,l_AccruedOnDefault); 

		if (Notional == 0 || Notional < 0) {	ERROR_MSG("Notional must be positive",pRet,ARM_ERROR_FREQ); return S_OK;	}
		VariantTools::convert(RecieveOrPay,strRecieveOrPay);
		if (strRecieveOrPay.empty()) { ERROR_MSG("Recieve or pay expected ",pRet,ARM_ERROR_FREQ); return S_OK;}
		if(strRecieveOrPay == "R" )
			ConvertRorP = 1;
		else if (strRecieveOrPay == "P")
			ConvertRorP = -1;
		else {	ERROR_MSG("recieveOrPay  : P or R ",pRet,ARM_ERROR_FREQ); return S_OK;}

		
		VariantTools::convert(RefValueSpreads,C_RefValueSpreads); 
		if (C_RefValueSpreads.empty()) {
				ERROR_MSG("RefValueSpreads expected",pRet,ARM_ERROR_FREQ); return S_OK;}
		VariantTools::convert(floatingIdx, C_floatingIdx);
		if (C_floatingIdx.empty()) {
			ERROR_MSG("floatingIdx expected",pRet,ARM_ERROR_FREQ); return S_OK; }

		VariantTools::convert(creditIdx,C_creditIdx);
		if (C_creditIdx.empty()) {
			ERROR_MSG("creditIdx expected",pRet,ARM_ERROR_FREQ); return S_OK; }

		VariantTools::convert(refvalueKUP, C_refvalueKUP);
		if (C_refvalueKUP.empty()) {
			ERROR_MSG("refvalueKUP expected",pRet,ARM_ERROR_FREQ); return S_OK;}
		
		VariantTools::convert(refvalueKDW, C_refvalueKDW);
		if (C_refvalueKDW.empty()) {
			ERROR_MSG("refvalueKDW expected",pRet,ARM_ERROR_FREQ); return S_OK; }
		
		VariantTools::convert(ScheduleInfoId, strSchedule_info);
		if (strSchedule_info.empty()) {
			ERROR_MSG("Schedule_info object Id is expected",pRet,ARM_ERROR_FREQ); return S_OK; }

		VariantTools::convert(Name, C_Name);
		if (C_Name.empty()) {
			ERROR_MSG("SName is expected",pRet,ARM_ERROR_FREQ); return S_OK; }

		VariantTools::convert(disc_ccy, C_disc_ccy);

		l_RefValueSpreads = LocalGetNumObjectId(CCString(C_RefValueSpreads.c_str()));
		l_floatingIdx = LocalGetNumObjectId(CCString(C_floatingIdx.c_str()));
		l_creditIdx = LocalGetNumObjectId(CCString(C_creditIdx.c_str()));
		l_refvalueKUP = LocalGetNumObjectId(CCString(C_refvalueKUP.c_str()));
		l_refvalueKDW = LocalGetNumObjectId(CCString(C_refvalueKDW.c_str()));
		l_schedule_Info = LocalGetNumObjectId(CCString(strSchedule_info.c_str()));

		long newId = ICMLOCAL_CORRIDORLEG_SCHE(C_Name, 
								 Notional,
								 ConvertRorP,
								 l_RefValueSpreads,
								 l_floatingIdx,
								 leverageFloatIdx,
								 l_creditIdx,
								 l_refvalueKUP,
								 l_refvalueKDW,
								 l_schedule_Info,
								 l_AccruedOnDefault,
								 C_disc_ccy,
								 C_result);
		string newLabel= LocalPersistent::get().getStringId(newId, LOCAL_ICM_CORRIDOR_LEG_CLASS);
		VariantTools::convert(newLabel, *pRet); 
		return S_OK ;
	} 
	catch(Exception&e)
	 {
	  return createErrorInfo("",e.GetErrorString()); 
	 }
	 catch(std::exception&e)
	 {
	  return createErrorInfo("",e.what()); 
	 }
	 catch(...) 
	 { 
	  return createErrorInfo("","Unknown exception"); 
	 }
	 return E_FAIL;
}


STDMETHODIMP ActiveXModule::ARM_Credit_CPDO(BSTR pRiskyLeg,
											 BSTR pRollLeg, 
											 BSTR pNoRiskyLeg, 
											 double pInitialValo, 
											 double pTarget,
											 double pMaturity, 
											 BSTR pCpnType, 
											 double	pUFFees, 
											 double	pRunningFees, 
											 double pVExpo, 
											 double	pV0Expo, 	
											 double	pAlpha,
											 double pBeta,	
											 double	pDesactivation,	
											 int	pNbAssets,
											 VARIANT *pRet)
{
	// TODO: Add your implementation code here
	try
	{
		_bstr_t bRiskyLeg(pRiskyLeg);
		CCString C_RiskyLeg = bRiskyLeg;

		_bstr_t bRollLeg(pRollLeg);
		CCString C_RollLeg = bRollLeg;

		_bstr_t bNoRiskyLeg(pNoRiskyLeg);
		CCString C_NoRiskyLeg = bNoRiskyLeg;

		_bstr_t bCpnType(pCpnType);
		CCString C_CpnType = bCpnType;
		
		ARM_result C_result;

		//_bstr_t bCurrency(pCurrency);
		//CCString C_Currency = bCurrency;

		CCString Res;
		CCString curClass = LOCAL_CPDO_CLASS;

		if(ICMLOCAL_CPDO(LocalGetNumObjectId(C_RiskyLeg),
							LocalGetNumObjectId(C_RollLeg),
							LocalGetNumObjectId(C_NoRiskyLeg),
							pMaturity,
							pInitialValo,
							pTarget,
							C_CpnType,
							pUFFees,
							pRunningFees,
							pVExpo,
							pV0Expo, 
							pAlpha,
							pBeta,
							pDesactivation, 
							pNbAssets,
    						C_result) == ARM_OK)
		{
			Res = LocalMakeObjectId (C_result.getLong (), curClass);
		}
		else
		{
		ERROR_MSG(C_result.getMsg(),pRet);
		return S_OK;
		}
		
		_variant_t wrap_Res;
		wrap_Res.Attach(*pRet);
		wrap_Res.SetString(Res);
		*pRet=wrap_Res.Detach();

		pRet->vt=VT_BSTR;

		return S_OK;
}
	 catch(Exception&e)
	 {
	  return createErrorInfo("ARM_Credit_CPDO",e.GetErrorString()); 
	 }
	 catch(std::exception&e)
	 {
	  return createErrorInfo("ARM_Credit_CPDO",e.what()); 
	 }
	 catch(...) 
	 { 
	  return createErrorInfo("ARM_Credit_CPDO","Unknown exception"); 
	 }
	 return E_FAIL;
}
STDMETHODIMP ActiveXModule::ARM_Credit_Restrikable_CDO(double TriggerStartDate,
													   double Expiry,
													   double Strike,
													   double InitSpread,
													   int	OptionType,
													   BSTR pUnderlying,
													   double Rehauss,
													   BSTR TriggerFreq,
													   int DiffCDO,
													   int IsCMSpread,
													   double CMSpreadMatu,
													   VARIANT *pRet)
{
	try
	{
		ARM_result C_result;

		_bstr_t b_Underlying(pUnderlying);
		CCString C_Underlying = b_Underlying;

		_bstr_t bFreqTrigger(TriggerFreq);
		CCString C_FreqTrigger = bFreqTrigger;

		int l_FreqTrigger;
		if((l_FreqTrigger = ARM_ConvFrequency (C_FreqTrigger, C_result)) == ARM_DEFAULT_ERR)
		{
			ERROR_MSG("Invalid FreqFeeLeg",pRet,ARM_ERROR_FREQ);
			return S_OK;
		}
		long newId = ICMLOCAL_RESTRIKABLE_CDO(&ARM_Date(XLDateToJulian(TriggerStartDate)),
											&ARM_Date(XLDateToJulian(Expiry)),
											Strike,
											InitSpread,
											OptionType,
											LocalGetNumObjectId (C_Underlying),
											Rehauss,
											l_FreqTrigger,
											IsCMSpread,
											CMSpreadMatu,
											DiffCDO);
		// TODO? créer une nouvelle classe option  pour l'affichage
		std::string objName = LocalPersistent::get().getStringId(newId,LOCAL_OPTION_CLASS) ;

		VariantTools::convert(objName, *pRet); 

		return S_OK;
	}

	catch(Exception&e)
	 {
	  return createErrorInfo("ARM_Credit_Restrikable_CDO",e.GetErrorString()); 
	 }
	 catch(std::exception&e)
	 {
	  return createErrorInfo("ARM_Credit_Restrikable_CDO",e.what()); 
	 }
	 catch(...) 
	 { 
	  return createErrorInfo("ARM_Credit_Restrikable_CDO","Unknown exception"); 
	 }
	 return E_FAIL;

}
