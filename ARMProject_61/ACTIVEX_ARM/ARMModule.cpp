// ARMModule.cpp : Implementation of CLocal_DLLARMApp and DLL registration.
#include "firsttobeincluded.h"
#include "stdafx.h"
#include <ARM/libarm_frometk/VariantTools.h>
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
#include "ARM_local_irindex.h"
#include "ARM_local_irfut.h"
#include "ARM_local_swtion.h"
#include "ARM_local_assetswap.h"
#include "ICM_local_pwccurve.h"
#include "ICM_local_glob.h"
#include "ICM_local_swap.h"
#include "ARM_local_refval.h"
#include "ARM_local_gp_inflation.h"
#include "ICM_local_mod.h"
#include "ICM_local_pricer.h"
#include <ARM\libarm_frometk\ARM_local_etoolkit.h>
#include <ARM\libarm_frometk\ARM_local_etoolkitX.h>
#include <ARM\libarm_frometk\ARM_local_xgiga.h>
#include <ARM_local_bond.h>
#include <ARM\libarm_local\ARM_local_gp_inflation.h>
#include <ARM\libarm_frometk\arm_local_paesexml_calypso.h> 
#include "ARM_local_interglob.h"
#include <ARM\libarm_frometk\PaserManagerUtilities.h>
#include <ARM\libicm_local\icm_local_summit.h>
#include <ARM\libarm_local\ARM_local_xstyle.h>
#include <ARM\libarm_local\ARM_local_pf.h>
#include <ARM\libarm_local\ARM_local_gp_model.h>
#include <ARM\libarm_local\ARM_local_capfl.h>

#include <gpinflation\infcurv.h>
#include <GP_Infra\gpinfra\argconvdefault.h>

using ARM::ARM_InfCurv;


#include <fromto.h>
#include <calend.h>
#include<string>


#ifdef _ACTIVEX_OLD_ 	
	#include "ARMModule.h"
#endif 

ARM_Date TMP_DATE;


/////////////////////////////////////////////////////////////////////////////
//
#ifdef _ACTIVEX_OLD_ 
STDMETHODIMP ARMCommonModule::InterfaceSupportsErrorInfo(REFIID riid)
{
	static const IID* arr[] = 
	{
		&IID_IARMModule,
	};

	for (int i=0;i<sizeof(arr)/sizeof(arr[0]);i++)
	{
		if (InlineIsEqualGUID(*arr[i],riid))
			return S_OK;
	}
	return S_FALSE;
}

#endif _ACTIVEX_OLD_


STDMETHODIMP ActiveXModule::ARMNextBusinessDay(double pDate, BSTR pCalendrier, long pNbDays, double *pDate2)
{
try
{
	// TODO: Add your implementation code here
	ARM_result C_result;

	_bstr_t bCalendrier(pCalendrier);
	CCString l_calendrier = bCalendrier;

	double dateRes;

	if(ARMLOCAL_NextBusinessDay(pDate,l_calendrier,pNbDays,C_result) == ARM_OK)
	{
		dateRes = Local_ARMDATE2XLDATE(C_result.getString ());
	}

	*pDate2 = dateRes;

	return S_OK;
}
catch(...)
{
	return E_FAIL;
}

}

STDMETHODIMP ActiveXModule::ARMAdjustToBusDate(double pDate, BSTR pCalendrier, BSTR pRule, double *pDate2)
{
	// TODO: Add your implementation code here
try
{
	ARM_result C_result;

	_bstr_t bCalendrier(pCalendrier);
	CCString l_calendrier = bCalendrier;

	_bstr_t bRule(pRule);
	CCString l_rule = bRule;

	long ruleId;

	if((ruleId = ARM_ConvRule (l_rule, C_result)) == ARM_DEFAULT_ERR)
		return S_FALSE;

	double dateRes;

	if(ARMLOCAL_ADJUSTTOBUSDATE(pDate,l_calendrier,ruleId,C_result) == ARM_OK)
	{
		dateRes = Local_ARMDATE2XLDATE(C_result.getString ());
	}

	*pDate2 = dateRes;

	return S_OK;
}
catch(...)
{
	return E_FAIL;
}
}


STDMETHODIMP ActiveXModule::ARMFreeObject(BSTR pId, long *pRet)
{
try
{
	_bstr_t bId(pId);
	CCString l_id = bId;

	ARM_result C_result;
	
	long idRes = -1;

	if(ARMLOCAL_FreeObject(LocalGetNumObjectId(l_id),C_result) == ARM_OK)
	{
		idRes = C_result.getLong ();
	}

	*pRet = idRes;

	return S_OK;
}
catch(...)
{
	return E_FAIL;
}
}


STDMETHODIMP ActiveXModule::ARMIsBusinessDay(double pDate, BSTR pCalendrier, long *pRes)
{
	// TODO: Add your implementation code here
try
{

	ARM_result C_result;
	
	_bstr_t bCalendrier(pCalendrier);
	CCString l_calendrier = bCalendrier;

	long Res;

	if(ARMLOCAL_IsBusinessDay(pDate,l_calendrier,C_result) == ARM_OK)
	{
		Res = C_result.getLong ();
	}

	*pRes = Res;

	return S_OK;
}
catch(...)
{
	return E_FAIL;
}
}


STDMETHODIMP ActiveXModule::ARMGetZCFromSummit(BSTR pIndex, BSTR pCurrency, BSTR pCvName, double pDate, BSTR pInterpMethod, BSTR *pRet)
{
	// TODO: Add your implementation code here
try
{
	ARM_result C_result;

	_bstr_t bIndex(pIndex);
	CCString l_index = bIndex;

	_bstr_t bCurrency(pCurrency);
	CCString l_currency = bCurrency;

	_bstr_t bCvName(pCvName);
	CCString l_cvname = bCvName;

	_bstr_t bInterpMethod(pInterpMethod);
	CCString l_InterpMethod = bInterpMethod;
	long InterpId;

	if((InterpId = ARM_ConvInterpMethod (l_InterpMethod, C_result)) == ARM_DEFAULT_ERR)
		return S_FALSE;

	CCString Res;
	CCString curClass = LOCAL_ZERO_CURVE_LIN_CLASS;

	if(ARMLOCAL_GetZCFromSummit(l_index,l_currency,l_cvname,pDate,InterpId,C_result) == ARM_OK)
	{
		Res = LocalMakeObjectId (C_result.getLong (), curClass);
	}

	_bstr_t tmpChaine = Res;

	*pRet = SysAllocString(tmpChaine);

	return S_OK;
}
catch(...)
{
	return E_FAIL;
}
}



STDMETHODIMP ActiveXModule::ARMGetZCFromCalypso(BSTR pIndex, 
												BSTR pCurrency, 
												BSTR pTerm,
												BSTR pPricingEnv, 
												double pDate, 
												BSTR pInterpMethod, 
												BSTR pForceCurveName,
												BSTR pXmlFileName,
												BSTR *pRet)
{
	// TODO: Add your implementation code here
	try
	{
		std::string index ; VariantTools::convert(pIndex,index); 
		std::string ccy ; VariantTools::convert(pCurrency,ccy); 
		std::string term ; VariantTools::convert(pTerm,term); 
		std::string pricingEnv; VariantTools::convert(pPricingEnv,pricingEnv); 
		ARM_Date AsOf; VariantTools::convertXLDate(pDate,AsOf); 
		std::string forceCurveName; VariantTools::convert(pForceCurveName,forceCurveName); 
		std::string xmlFileName; VariantTools::convert(pXmlFileName,xmlFileName); 
		long InterpId; 
		{
			// TODO: 
			// 		ARM_ConvInterpMethod uses const string& 
			//		and throw on failure... 
			//		and returns enumerates
			//
			_bstr_t bInterpMethod(pInterpMethod);
			CCString l_InterpMethod = bInterpMethod;
			ARM_result C_result; 

			if((InterpId = ARM_ConvInterpMethod (l_InterpMethod, C_result)) == ARM_DEFAULT_ERR)
				return S_FALSE;
		}

		CCString Res;
		long objId = ARMLOCAL_GetZCFromCalypso(AsOf,index,ccy,term,pricingEnv,forceCurveName,InterpId,xmlFileName) ;
		Res = LocalMakeObjectId (objId, LOCAL_ZERO_CURVE_LIN_CLASS);
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
	return S_OK;
}


/**

STDMETHODIMP ActiveXModule::ARMGetZCFromCalypso(BSTR pIndex, BSTR pCurrency, BSTR pTerm, BSTR pPricingEnv, double pDate, BSTR pInterpMethod, BSTR forceCurveName, BSTR xmlFile, BSTR *pRet)
{
	// TODO: Add your implementation code here
try
{
	ARM_result C_result;

	_bstr_t bIndex(pIndex);
	CCString l_index = bIndex;

	_bstr_t bCurrency(pCurrency);
	CCString l_currency = bCurrency;

	_bstr_t bCvName(pCvName);
	CCString l_cvname = bCvName;

	_bstr_t bInterpMethod(pInterpMethod);
	CCString l_InterpMethod = bInterpMethod;
	long InterpId;

	if((InterpId = ARM_ConvInterpMethod (l_InterpMethod, C_result)) == ARM_DEFAULT_ERR)
		return S_FALSE;

	CCString Res;
	CCString curClass = LOCAL_ZERO_CURVE_LIN_CLASS;

	if(ARMLOCAL_GetZCFromSummit(l_index,l_currency,l_cvname,pDate,InterpId,C_result) == ARM_OK)
	{
		Res = LocalMakeObjectId (C_result.getLong (), curClass);
	}

	_bstr_t tmpChaine = Res;

	*pRet = SysAllocString(tmpChaine);

	return S_OK;
}
catch(...)
{
	return E_FAIL;
}
}
**/ 


STDMETHODIMP ActiveXModule::ARMCreateZCFromSummit(BSTR pIndex, BSTR pCurrency, BSTR pCvName, double pDate,BSTR pAdj,BSTR pRaw, BSTR *pRet)
{
	// TODO: Add your implementation code here
try
{
	ARM_result C_result;

	_bstr_t bIndex(pIndex);
	CCString l_index = bIndex;

	_bstr_t bCurrency(pCurrency);
	CCString l_currency = bCurrency;

	_bstr_t bCvName(pCvName);
	CCString l_cvname = bCvName;

	_bstr_t bAdj(pAdj);
	CCString l_Adj = bAdj;
	long adjId;
	
	_bstr_t bRaw(pRaw);
	CCString l_Raw = bRaw;

	if((adjId = ARM_ConvYesOrNo (l_Adj, C_result)) == ARM_DEFAULT_ERR)
		return S_FALSE;

	CCString Res;
	CCString curClass = LOCAL_ZERO_CURVE_LIN_CLASS;

	if(ARMLOCAL_CreateZCFromSummit(l_index,l_currency,l_cvname,pDate,adjId,l_Raw,K_DEF_FREQ,C_result) == ARM_OK)
	{
		Res = LocalMakeObjectId (C_result.getLong (), curClass);
	}

	_bstr_t tmpChaine = Res;

	*pRet = SysAllocString(tmpChaine);

	return S_OK;
}
catch(...)
{
	return E_FAIL;
}
}


STDMETHODIMP ActiveXModule::ARMBumpCurve(BSTR pZc, double pEpsilon, long pMethod, BSTR pPlot, BSTR *pRet)
{
	// TODO: Add your implementation code here
try
{
	ARM_result C_result;

	_bstr_t bZc(pZc);
	CCString l_Zc = bZc;

	_bstr_t bPlot(pPlot);
	CCString l_Plot = bPlot;

	ARM_ZeroCurve* inCurve = NULL;
	ARM_InfCurv* inInfCurve = NULL;

	inCurve = (ARM_ZeroCurve *) LOCAL_PERSISTENT_OBJECTS->GetObjectA(LocalGetNumObjectId(l_Zc));
	if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(inCurve, ARM_ZERO_CURVE) == 0)
	{
		inInfCurve = (ARM_InfCurv *) LOCAL_PERSISTENT_OBJECTS->GetObjectA(LocalGetNumObjectId(l_Zc));
	}

	vector<CCString> psMatu;
	vector<double> epsilon;

	if (strcmp((const char*) l_Plot,"") == 0)
	{
		epsilon.push_back(pEpsilon);
	}
	else
	{
		psMatu.push_back(l_Plot);
		epsilon.push_back(pEpsilon);
	}

	CCString Res;
	CCString curClass = LocalGetStringObjectClass(l_Zc);

	if(ARMLOCAL_BumpCurve(LocalGetNumObjectId(l_Zc),psMatu,epsilon,pMethod,C_result) == ARM_OK)
	{
		Res = LocalMakeObjectId (C_result.getLong (), curClass);
	}

	_bstr_t tmpChaine = Res;

	*pRet = SysAllocString(tmpChaine);

	return S_OK;
}
catch(...)
{
	return E_FAIL;
}
}


STDMETHODIMP ActiveXModule::ARMFreeAllObjects(long *pRet)
{
	// TODO: Add your implementation code here

try
{
	ARM_result C_result;
	
	long idRes = -1;

	if(ARMLOCAL_FreeAllObjects(C_result) == ARM_OK)
	{
		idRes = C_result.getLong ();
	}

	*pRet = idRes;

	return S_OK;
}
catch(...)
{
	return E_FAIL;
}
}


STDMETHODIMP ActiveXModule::ARMYcMod(BSTR pZc, BSTR pZcDiscount, BSTR *pRet)
{
	// TODO: Add your implementation code here

try
{
	ARM_result C_result;
	
	_bstr_t bZc(pZc);
	CCString l_zc = bZc;

	_bstr_t bZcDiscount(pZcDiscount);
	CCString l_zcDiscount = bZcDiscount;

	if ( strcmp((const char*)l_zcDiscount,"NULL") == 0)
		l_zcDiscount = l_zc;

	CCString Res;
	CCString curClass = LOCAL_YIELD_CURVE_BASIC_CLASS;

	if(ARMLOCAL_ycmod(LocalGetNumObjectId(l_zc),LocalGetNumObjectId(l_zcDiscount),C_result) == ARM_OK)
	{
		Res = LocalMakeObjectId (C_result.getLong(), curClass);
	}

	_bstr_t tmpChaine = Res;

	*pRet = SysAllocString(tmpChaine);

	return S_OK;
}
catch(...)
{
	return E_FAIL;
}
}


STDMETHODIMP ActiveXModule::ARMForwardYield(BSTR pZc, double pMatu1, double pMatu2, BSTR pMeth, BSTR pAdj, double *pRet)
{
	// TODO: Add your implementation code here
try
{
	ARM_result C_result;

	_bstr_t bZc(pZc);
	CCString l_zc = bZc;

	_bstr_t bMeth(pMeth);
	CCString l_meth = bMeth;
	long methId;

	if((methId = ARM_ConvForwardYieldMethod (l_meth, C_result)) == ARM_DEFAULT_ERR)
		return S_FALSE;

	_bstr_t bAdj(pAdj);
	CCString l_adj = bAdj;
	long adjId;

	if((adjId = ARM_ConvYesOrNo (l_adj, C_result)) == ARM_DEFAULT_ERR)
		return S_FALSE;

	double Res;

	if(ARMLOCAL_ForwardYield(LocalGetNumObjectId(l_zc), pMatu1, pMatu2, methId, adjId, -9999, -1, C_result) == ARM_OK)
	{
		Res = C_result.getDouble();
	}

	*pRet = Res;

	return S_OK;
}
catch(...)
{
	return E_FAIL;
}
}


STDMETHODIMP ActiveXModule::ARMDiscountYield(VARIANT *pZc, VARIANT *pMatu, VARIANT *pMeth, VARIANT *pRet)
{
	// TODO: Add your implementation code here
try
{
	_variant_t valMatu( pMatu );

	valMatu.ChangeType( VT_R8 );

	ARM_result C_result;

	CCString l_zcId;
	
	CCString l_meth;
	long methId;

	if(VARIANT2CCString (*pZc, l_zcId) != S_OK)
		return S_FALSE;

	if(VARIANT2CCString (*pMeth, l_meth) != S_OK)
		return S_FALSE;

	if((methId = ARM_ConvForwardYieldMethod (l_meth, C_result)) == ARM_DEFAULT_ERR)
		return S_FALSE;

	double Res;

	if(ARMLOCAL_DiscountYield(LocalGetNumObjectId(l_zcId), valMatu, methId, C_result) == ARM_OK)
	{
		Res = C_result.getDouble();
	}

	_variant_t wrap_Res;
	wrap_Res.Attach(*pRet);
	wrap_Res.dblVal = Res;
	*pRet=wrap_Res.Detach();

	pRet->vt=VT_R8;

	return S_OK;
}
catch(...)
{
	return E_FAIL;
}
}


STDMETHODIMP ActiveXModule::ARMLiborSwap(VARIANT *pStartDate, VARIANT *pEndDate, VARIANT *pLiborType, VARIANT *pRecOrPay, VARIANT *pFixedRate, VARIANT *pSpread, VARIANT *pCcy, BSTR pDaycount, BSTR pFloatingDaycount, VARIANT *pRet)
{
	// TODO: Add your implementation code here
try
{
	_variant_t valStartDate( pStartDate );
	_variant_t valEndDate( pEndDate );
	_variant_t valFixedRate( pFixedRate );
	_variant_t valSpread( pSpread );

	valStartDate.ChangeType( VT_R8 );
	valEndDate.ChangeType( VT_R8 );
	valFixedRate.ChangeType( VT_R8 );
	valSpread.ChangeType( VT_R8 );

	_bstr_t bDaycount(pDaycount);
	CCString l_Daycount = bDaycount;

	_bstr_t bFloatingDaycount(pFloatingDaycount);
	CCString l_FloatingDaycount = bFloatingDaycount;

	ARM_result C_result;

	CCString l_liborType;
	long liborTypeId;
	
	CCString l_recOrPay;
	long recOrPayId;

	CCString l_ccy;

	if(VARIANT2CCString (*pLiborType, l_liborType) != S_OK)
		return S_FALSE;

	if(VARIANT2CCString (*pRecOrPay, l_recOrPay) != S_OK)
		return S_FALSE;

	if(VARIANT2CCString (*pCcy, l_ccy) != S_OK)
		return S_FALSE;

	if((liborTypeId = ARM_ConvIrIndName (l_liborType, C_result)) == ARM_DEFAULT_ERR)
		return S_FALSE;

	if((recOrPayId = ARM_ConvRecOrPay (l_recOrPay, C_result)) == ARM_DEFAULT_ERR)
		return S_FALSE;

	long dayCountId = ARM_ConvDayCount(l_Daycount);
	long floatingDayCountId = ARM_ConvDayCount(l_FloatingDaycount);

	CCString Res;
	CCString curClass = LOCAL_SWAP_CLASS;

	if(ARMLOCAL_LIBORSWAP (valStartDate, valEndDate, liborTypeId, recOrPayId, valFixedRate, 0, valSpread, true, l_ccy, dayCountId, floatingDayCountId, C_result) == ARM_OK)
	{
		Res = LocalMakeObjectId (C_result.getLong (), curClass);
	}

	_variant_t wrap_Res;
	wrap_Res.Attach(*pRet);
	wrap_Res.SetString(Res);
	*pRet=wrap_Res.Detach();

	pRet->vt=VT_BSTR;

	return S_OK;
}
catch(...)
{
	return E_FAIL;
}
}


STDMETHODIMP ActiveXModule::ARMSwapPriceToRate(VARIANT *pSwap, VARIANT *pDate, VARIANT *pPrice, VARIANT *pModel, VARIANT *pRet)
{
	// TODO: Add your implementation code here
try
{
	_variant_t valDate( pDate );
	_variant_t valPrice( pPrice );

	valDate.ChangeType( VT_R8 );
	valPrice.ChangeType( VT_R8 );

	ARM_result C_result;

	CCString l_swap;
	
	CCString l_model;

	if(VARIANT2CCString (*pSwap, l_swap) != S_OK)
		return S_FALSE;

	if(VARIANT2CCString (*pModel, l_model) != S_OK)
		return S_FALSE;

	double Res;

	if(ARMLOCAL_SWAP_PRICE_TO_RATE (LocalGetNumObjectId (l_swap), valDate, valPrice, LocalGetNumObjectId (l_model), C_result) == ARM_OK)
	{
		Res = C_result.getDouble ();
	}

	_variant_t wrap_Res;
	wrap_Res.Attach(*pRet);
	wrap_Res.dblVal = Res;
	*pRet=wrap_Res.Detach();

	pRet->vt=VT_R8;

	return S_OK;
}
catch(...)
{
	return E_FAIL;
}
}


STDMETHODIMP ActiveXModule::ARMPrice(VARIANT *pSec, VARIANT *pModel, VARIANT *pRet)
{
	// TODO: Add your implementation code here

try
{
	ARM_result C_result;

	CCString l_sec;
	
	CCString l_model;

	if(VARIANT2CCString (*pSec, l_sec) != S_OK)
		return S_FALSE;

	if(VARIANT2CCString (*pModel, l_model) != S_OK)
		return S_FALSE;

	double Res;

	if(ARMLOCAL_ARM_Price (LocalGetNumObjectId (l_sec), LocalGetNumObjectId (l_model), C_result) == ARM_OK)
	{
		Res = C_result.getDouble ();
	}

	_variant_t wrap_Res;
	wrap_Res.Attach(*pRet);
	wrap_Res.dblVal = Res;
	*pRet=wrap_Res.Detach();

	pRet->vt=VT_R8;

	return S_OK;
}
catch(...)
{
	return E_FAIL;
}
}


STDMETHODIMP ActiveXModule::ARMBetweenDates(VARIANT *pDate1, VARIANT *pDate2, VARIANT *pDaycount, VARIANT *pIsYearFrac, VARIANT *pRet)
{
	// TODO: Add your implementation code here
try
{
	_variant_t valDate1( pDate1 );
	_variant_t valDate2( pDate2 );
	_variant_t valIsYearFrac( pIsYearFrac );

	valDate1.ChangeType( VT_R8 );
	valDate2.ChangeType( VT_R8 );
	valIsYearFrac.ChangeType( VT_I4 );

	ARM_result C_result;

	CCString l_daycount;
	long daycountId;

	if(VARIANT2CCString (*pDaycount, l_daycount) != S_OK)
		return S_FALSE;

	daycountId = ARM_ConvDayCount (l_daycount);
	
	double Res;

	if(ARMLOCAL_ARM_BetweenDates (valDate1, valDate2, daycountId, valIsYearFrac, C_result) == ARM_OK)
	{
		Res = C_result.getDouble ();
	}

	_variant_t wrap_Res;
	wrap_Res.Attach(*pRet);
	wrap_Res.dblVal = Res;
	*pRet=wrap_Res.Detach();

	pRet->vt=VT_R8;

	return S_OK;
}
catch(...)
{
	return E_FAIL;
}
}


STDMETHODIMP ActiveXModule::ARMAddPeriod(VARIANT *pDate, VARIANT *pFreq, VARIANT *pCcy, VARIANT *pNbPeriod, VARIANT *pAdjRule, VARIANT *pRet)
{
	// TODO: Add your implementation code here
try
{
	_variant_t valDate( pDate );
	_variant_t valNbPeriods( pNbPeriod );

	valDate.ChangeType( VT_R8 );
	valNbPeriods.ChangeType( VT_I4 );

	ARM_result C_result;

	CCString l_freq;
	long freqId;

	CCString l_ccy;
	CCString l_adjRule;
	long adjRuleId;

	if(VARIANT2CCString (*pFreq, l_freq) != S_OK)
		return S_FALSE;
	
	if(VARIANT2CCString (*pCcy, l_ccy) != S_OK)
		return S_FALSE;

	if(VARIANT2CCString (*pAdjRule, l_adjRule) != S_OK)
		return S_FALSE;

	if((adjRuleId = ARM_ConvRule (l_adjRule, C_result)) == ARM_DEFAULT_ERR)
		return S_FALSE;

	double dateRes;
	long retCode;

	if ( (freqId = ARM_ConvFrequency(l_freq,C_result)) == ARM_DEFAULT_ERR)
	{
		int Nb;
		char matu;

		sscanf(l_freq, "%d%c", &Nb, &matu);

        matu = toupper(matu);

        if ( matu == 'D' ) // Ex : "1D"
        {    
			retCode = ARMLOCAL_ARM_ADDPERIOD(valDate, K_DAILY, l_ccy, (long) (valNbPeriods.intVal * Nb), adjRuleId, 0, C_result);
        }
        else if ( matu == 'W' )  
        {   //  Ex : "1W"    

			retCode = ARMLOCAL_ARM_ADDPERIOD(valDate, K_WEEKLY, l_ccy, (long) (valNbPeriods.intVal * Nb), adjRuleId, 0, C_result);
        }
        else if ( matu == 'M' ) 
        {   //  Ex : "9M"
			retCode = ARMLOCAL_ARM_ADDPERIOD(valDate, K_MONTHLY, l_ccy, (long) (valNbPeriods.intVal * Nb), adjRuleId, 0, C_result);
        }
        else if ( matu == 'Y')  // ->implicitement ce sont des taux de swap
        {   
			retCode = ARMLOCAL_ARM_ADDPERIOD(valDate, K_ANNUAL, l_ccy, (long) (valNbPeriods.intVal * Nb), adjRuleId, 0, C_result);
		}
		else
			return S_FALSE;
	}
	else
	{
		retCode = ARMLOCAL_ARM_ADDPERIOD(valDate, freqId, l_ccy, (long) valNbPeriods.intVal, adjRuleId, 0, C_result);
	}

	dateRes = Local_ARMDATE2XLDATE(C_result.getString ());

	_variant_t wrap_Res;
	wrap_Res.Attach(*pRet);
	wrap_Res.dblVal = dateRes;
	*pRet=wrap_Res.Detach();

	pRet->vt=VT_R8;

	return S_OK;
}
catch(...)
{
	return E_FAIL;
}
}


STDMETHODIMP ActiveXModule::ARMIsoCcy(BSTR pCcy, BSTR pRefObj, BSTR *pRet)
{
	// TODO: Add your implementation code here
try
{
	ARM_result C_result;

	_bstr_t bCcy(pCcy);
	CCString l_ccy = bCcy;

	_bstr_t bRefObj(pRefObj);
	CCString l_refXL = bRefObj;

	CCString curClass = LOCAL_CCY_CLASS;
	CCString stringId = GetEnvVar (l_refXL);
	
	long retCode;
	long objId;

	if(!stringId)
	{
		retCode = ARMLOCAL_ISOCCY (l_ccy, C_result) ;

		if(retCode == ARM_OK)
		{
			objId = C_result.getLong ();

			if (l_refXL.GetLen()>0)
						LocalSetCurCellEnvValue (curClass, objId,l_refXL); 

			stringId = LocalMakeObjectId (objId, curClass);
		}
	}
	else
	{
		CCString prevClass = LocalGetStringObjectClass (stringId);
		
		objId = LocalGetNumObjectId (stringId);
			
		if(curClass == prevClass)
		{
			retCode = ARMLOCAL_ISOCCY (l_ccy, C_result,objId) ;	
 		}
		else
		{
			FreeCurCellContent (stringId);
			retCode = ARMLOCAL_ISOCCY (l_ccy, C_result) ;

			if(retCode == ARM_OK)
			{
				objId = C_result.getLong ();
			
				LocalSetCurCellEnvValue (curClass, objId,l_refXL); 

				stringId = LocalMakeObjectId (objId, curClass);
			}
		}
	}

	_bstr_t tmpChaine = stringId;

	*pRet = SysAllocString(tmpChaine);

	return S_OK;
}
catch(...)
{
	return E_FAIL;
}
}

STDMETHODIMP ActiveXModule::ARMGetSpotDays(VARIANT *pCcy, VARIANT *pRet)
{
	// TODO: Add your implementation code here
try
{
	ARM_result C_result;

	CCString l_ccy;
	
	if(VARIANT2CCString (*pCcy, l_ccy) != S_OK)
		return S_FALSE;

	int Res;

	if(ARMLOCAL_GetSpotDays (l_ccy, C_result) == ARM_OK)
	{
		Res = (int) (C_result.getDouble());
	}

	_variant_t wrap_Res;
	wrap_Res.Attach(*pRet);
	wrap_Res.iVal = Res;
	*pRet=wrap_Res.Detach();

	pRet->vt=VT_R8;

	return S_OK;
}
catch(...)
{
	return E_FAIL;
}
}


STDMETHODIMP ActiveXModule::ARMGetLiborIndexDaycount(VARIANT *pCcy, VARIANT *pRet)
{
	// TODO: Add your implementation code here
try
{
	ARM_result C_result;

	CCString l_ccy;
	
	if(VARIANT2CCString (*pCcy, l_ccy) != S_OK)
		return S_FALSE;

	char* sCcy = l_ccy.GetStr();

	ARM_Currency ccy (sCcy);

	delete sCcy;

	int Res = ccy.GetLiborIndexDayCount();

	_variant_t wrap_Res;
	wrap_Res.Attach(*pRet);
	wrap_Res.intVal = Res;
	*pRet=wrap_Res.Detach();

	pRet->vt=VT_I4;

	return S_OK;
}
catch(...)
{
	return E_FAIL;
}
}

STDMETHODIMP ActiveXModule::ARMGetLiborTerm(VARIANT *pCcy, VARIANT *pRet)
{
	// TODO: Add your implementation code here
try
{
	ARM_result C_result;

	CCString l_ccy;
	
	if(VARIANT2CCString (*pCcy, l_ccy) != S_OK)
		return S_FALSE;

	char* sCcy = l_ccy.GetStr();

	ARM_Currency ccy (sCcy);

	delete sCcy;

	int Res = ccy.GetLiborTerm();

	_variant_t wrap_Res;
	wrap_Res.Attach(*pRet);
	wrap_Res.intVal = Res;
	*pRet=wrap_Res.Detach();

	pRet->vt=VT_I4;

	return S_OK;
}
catch(...)
{
	return E_FAIL;
}
}

STDMETHODIMP ActiveXModule::ARMGetFixedDayCount(VARIANT *pCcy, VARIANT *pRet)
{
	// TODO: Add your implementation code here
try
{
	ARM_result C_result;

	CCString l_ccy;
	
	if(VARIANT2CCString (*pCcy, l_ccy) != S_OK)
		return S_FALSE;

	char* sCcy = l_ccy.GetStr();

	ARM_Currency ccy (sCcy);

	delete sCcy;

	int Res = ccy.GetFixedDayCount();

	_variant_t wrap_Res;
	wrap_Res.Attach(*pRet);
	wrap_Res.intVal = Res;
	*pRet=wrap_Res.Detach();

	pRet->vt=VT_I4;

	return S_OK;
}
catch(...)
{
	return E_FAIL;
}
}

STDMETHODIMP ActiveXModule::ARMGetFixedPayFreq(VARIANT *pCcy, VARIANT *pRet)
{
	// TODO: Add your implementation code here
try
{
	ARM_result C_result;

	CCString l_ccy;
	
	if(VARIANT2CCString (*pCcy, l_ccy) != S_OK)
		return S_FALSE;

	char* sCcy = l_ccy.GetStr();

	ARM_Currency ccy (sCcy);

	delete sCcy;

	int Res = ccy.GetFixedPayFreq();

	_variant_t wrap_Res;
	wrap_Res.Attach(*pRet);
	wrap_Res.intVal = Res;
	*pRet=wrap_Res.Detach();

	pRet->vt=VT_I4;

	return S_OK;
}
catch(...)
{
	return E_FAIL;
}
}


STDMETHODIMP ActiveXModule::ARMComputeVolatility(VARIANT *pVol, VARIANT* pMatu, VARIANT* pStrike, VARIANT* pTenor, VARIANT *pRet)
{
	// TODO: Add your implementation code here
try
{
	ARM_result C_result;

	double Res;
	
	_variant_t valMatu( pMatu );
	valMatu.ChangeType( VT_R8 );

	_variant_t valStrike( pStrike );
	valStrike.ChangeType( VT_R8 );

	_variant_t valTenor( pTenor);
	valTenor.ChangeType( VT_R8 );

	CCString l_vol;

	if(VARIANT2CCString (*pVol, l_vol) != S_OK)
		return S_FALSE;

	if(ARMLOCAL_ComputeVolatility (LocalGetNumObjectId (l_vol), valMatu, valStrike, valTenor, C_result) == ARM_OK)
	{
		Res = C_result.getDouble ();
	}

	_variant_t wrap_Res;
	wrap_Res.Attach(*pRet);
	wrap_Res.dblVal = Res;
	*pRet=wrap_Res.Detach();

	pRet->vt=VT_R8;

	return S_OK;
}
catch(...)
{
	return E_FAIL;
}
}

STDMETHODIMP ActiveXModule::ARMVolCurv(VARIANT *pMatu, VARIANT* pStrikes, VARIANT* pVols, double pAsOf, BSTR pStkType, BSTR pVolType, BSTR pCcy, BSTR pIndexId,BSTR *pRet)
{
	try
	{
		ARM_result C_result;

		CCString Res;
		CCString curClass = LOCAL_VOL_CURVE_LIN_CLASS;

		long matuSize;
		long strikesSize;
		long volsSize;

		VECTOR<double> vMatu;
		VECTOR<double> vStrikes;
		VECTOR<double> vVols;

		_bstr_t bStkType(pStkType);
		CCString l_stkType = bStkType;
		_bstr_t bVolType(pVolType);
		CCString l_volType = bVolType;
		_bstr_t bCcy(pCcy);
		CCString l_ccy = bCcy;
		std::string indexId; VariantTools::convert(pIndexId,indexId); 
		long indexId_= LocalPersistent::get().getObjectId(indexId);


		long stkTypeId;
		long volTypeId;

		if(VARIANT2VECTORDOUBLE (*pMatu,vMatu,matuSize) != S_OK)
			return S_FALSE;

		if(VARIANT2VECTORDOUBLE (*pStrikes,vStrikes,strikesSize) != S_OK)
			return S_FALSE;

		if(VARIANT2VECTORDOUBLE (*pVols,vVols,volsSize) != S_OK)
			return S_FALSE;

		if((volTypeId = ARM_ConvVolType (l_volType, C_result)) == ARM_DEFAULT_ERR)
		{
			return S_FALSE;
		}

		if((stkTypeId = StrikeCode (l_stkType, C_result)) == ARM_DEFAULT_ERR)
		{
			return S_FALSE;
		}

		if(ARMLOCAL_volcurv (vMatu,vStrikes,vVols,pAsOf,stkTypeId,volTypeId,K_LINEAR,l_ccy,"",indexId_,C_result) == ARM_OK)
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
	return S_OK;
}



STDMETHODIMP ActiveXModule::ARMGetVolCubeFromSummit(BSTR pIndex, BSTR pCcy, BSTR pCvName, double pAsOf, BSTR pType, VARIANT* pSmiles, BSTR pTypeCube, BSTR pIndexId,BSTR *pRet)
{
try
{
	ARM_result C_result;

	_bstr_t bIndex(pIndex);
	CCString l_Index = bIndex;

	_bstr_t bCcy(pCcy);
	CCString l_Ccy = bCcy;

	_bstr_t bCvName(pCvName);
	CCString l_CvName = bCvName;

	_bstr_t bType(pType);
	CCString l_Type = bType;

	_bstr_t bTypeCube(pTypeCube);
	CCString l_TypeCube = bTypeCube;

	VECTOR<CCString> vSmiles;
	CCString defaultSmile("");
	long size;

	try{
		if(VARIANT2CCString (*pSmiles,defaultSmile) != S_OK)
			return S_FALSE;
			
	}
	catch (...)
	{
		if(VARIANT2VECTORCCSTRING (*pSmiles,vSmiles,size) != S_OK)
			return S_FALSE;
	}
	
	if(strcmp((const char*)defaultSmile,"") != 0)
	{
		if(VARIANT2VECTORCCSTRING (*pSmiles,vSmiles,size) != S_OK)
			return S_FALSE;
	}

	CCString Res;
	CCString curClass = LOCAL_VOL_CUBE_CLASS;

	std::string indexId; VariantTools::convert(pIndexId,indexId); 
	long indexId_= LocalPersistent::get().getObjectId(indexId);

	if(ARMLOCAL_GetVolCubeFromSummit (l_Index, l_Ccy, l_CvName, pAsOf, l_Type, vSmiles, l_TypeCube, indexId_, 0, C_result) == ARM_OK)
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


STDMETHODIMP ActiveXModule::ARMParallelShift(BSTR pZc, double pBump, BSTR *pRet)
{
	// TODO: Add your implementation code here
try
{
	ARM_result C_result;

	_bstr_t bZc(pZc);
	CCString l_zc = bZc;

	CCString Res;
	CCString curClass = LocalGetStringObjectClass(l_zc);

	if(ARMLOCAL_ParallelShift (LocalGetNumObjectId (l_zc), pBump, C_result) == ARM_OK)
	{
		Res = LocalMakeObjectId (C_result.getLong (), curClass);
	}

	_bstr_t tmpChaine = Res;

	*pRet = SysAllocString(tmpChaine);

	return S_OK;
}
catch(...)
{
	return E_FAIL;
}
}


STDMETHODIMP ActiveXModule::ARMBumpVolatility(BSTR pVol, double pValue, long pNthLine, long pNthCol, BSTR pIsCumul, BSTR *pRet)
{
	// TODO: Add your implementation code here
try
{
	ARM_result C_result;

	_bstr_t bVol(pVol);
	CCString l_vol = bVol;
	
	_bstr_t bIsCumul(pIsCumul);
	CCString l_iscumul = bIsCumul;
	long isCumulId;

	if((isCumulId = ARM_ConvYesOrNo (l_iscumul, C_result)) == ARM_DEFAULT_ERR)
	{
		return S_FALSE;
	}

	CCString Res;
	CCString curClass = LocalGetStringObjectClass(l_vol);

	if(ARMLOCAL_ARM_BumpVolatility (LocalGetNumObjectId (l_vol), pValue, pNthLine, pNthCol, isCumulId, K_YES, C_result) == ARM_OK)
	{
		Res = LocalMakeObjectId (C_result.getLong (), curClass);
	}

	_bstr_t tmpChaine = Res;

	*pRet = SysAllocString(tmpChaine);

	return S_OK;
}
catch(...)
{
	return E_FAIL;
}
}


STDMETHODIMP ActiveXModule::ARMGlobDFBS(BSTR pDomBSId, BSTR pDomCurrId, BSTR pFrgBSId, BSTR pFrgCurrId, BSTR pFxVolCrvId, BSTR pFFxCorrId, BSTR pRatesCorrId, BSTR pFxVolModelId, BSTR *pRet)
{

try
{
	ARM_result C_result;

	_bstr_t bDomBSId(pDomBSId);
	CCString l_DomBSId = bDomBSId;

	_bstr_t bDomCurrId(pDomCurrId);
	CCString l_DomCurrId = bDomCurrId;
	
	_bstr_t bFrgBSId(pFrgBSId);
	CCString l_FrgBSId = bFrgBSId;
	
	_bstr_t bFrgCurrId(pFrgCurrId);
	CCString l_FrgCurrId = bFrgCurrId;
	
	_bstr_t bFxVolCrvId(pFxVolCrvId);
	CCString l_FxVolCrvId= bFxVolCrvId;
	
	_bstr_t bFFxCorrId(pFFxCorrId);
	CCString l_FFxCorrId= bFFxCorrId;

	_bstr_t bRatesCorrId(pRatesCorrId);
	CCString l_RatesCorrId= bRatesCorrId;

	_bstr_t bFxVolModelId(pFxVolModelId);
	CCString l_FxVolModelId= bFxVolModelId;

	long fxVolModelId = ARM_NULL_OBJECT;
	if ( strcmp((const char*)l_FxVolModelId,"DEFAULT") != 0)
		fxVolModelId = LocalGetNumObjectId(l_FxVolModelId);

	CCString Res;
	CCString curClass = LOCAL_DFFXBS_CLASS;

	if(ARMLOCAL_GLOBDFBS(LocalGetNumObjectId(l_DomBSId),
						 LocalGetNumObjectId(l_DomCurrId),
						 LocalGetNumObjectId(l_FrgBSId),
						 LocalGetNumObjectId(l_FrgCurrId),
						 LocalGetNumObjectId(l_FxVolCrvId),
						 LocalGetNumObjectId(l_FFxCorrId),
						 LocalGetNumObjectId(l_RatesCorrId),
						 fxVolModelId,
						 C_result) == ARM_OK)
	{
		Res = LocalMakeObjectId (C_result.getLong (), curClass);
	}

	_bstr_t tmpChaine = Res;

	*pRet = SysAllocString(tmpChaine);

	return S_OK;
}
catch(...)
{
	return E_FAIL;
}

}

STDMETHODIMP ActiveXModule::ARMDFFXBS(BSTR pDVolId,BSTR pFVolId,BSTR pDZcId,BSTR pFZcId,BSTR pDFxCorrId,BSTR pFFxCorrId,BSTR pFxVolId,double pRatesCorr,BSTR *pRet)
{

try
{
	ARM_result C_result;

	_bstr_t bDVolId(pDVolId);
	CCString l_DVolId = bDVolId;

	_bstr_t bFVolId(pFVolId);
	CCString l_FVolId = bFVolId;
	
	_bstr_t bDZcId(pDZcId);
	CCString l_DZcId = bDZcId;
	
	_bstr_t bFZcId(pFZcId);
	CCString l_FZcId = bFZcId;

	_bstr_t bFxVolId(pFxVolId);
	CCString l_FxVolId= bFxVolId;

	_bstr_t bDFxCorrId(pDFxCorrId);
	CCString l_DFxCorrId= bDFxCorrId;
	
	_bstr_t bFFxCorrId(pFFxCorrId);
	CCString l_FFxCorrId= bFFxCorrId;

	CCString Res;
	CCString curClass = LOCAL_DFFXBS_CLASS;

	if(ARMLOCAL_DFFXBS (LocalGetNumObjectId (l_DVolId),
						LocalGetNumObjectId (l_FVolId),
						LocalGetNumObjectId (l_DZcId),
						LocalGetNumObjectId (l_FZcId),
						LocalGetNumObjectId (l_DFxCorrId),
						LocalGetNumObjectId (l_FFxCorrId),
						LocalGetNumObjectId (l_FxVolId),
						pRatesCorr,
						ARM_NULL_OBJECT,
						ARM_NULL_OBJECT,
						0.0,
						ARM_NULL_OBJECT,
						ARM_NULL_OBJECT,
						0.0,
						ARM_NULL_OBJECT,
						ARM_NULL_OBJECT,
						ARM_NULL_OBJECT,
						ARM_NULL_OBJECT,
						ARM_NULL_OBJECT,
						ARM_NULL_OBJECT,
						ARM_NULL_OBJECT,
						ARM_NULL_OBJECT,
						ARM_NULL_OBJECT,
						C_result) == ARM_OK)
	{
		Res = LocalMakeObjectId (C_result.getLong (), curClass);
	}

	_bstr_t tmpChaine = Res;

	*pRet = SysAllocString(tmpChaine);

	return S_OK;
}
catch(...)
{
	return E_FAIL;
}

}


STDMETHODIMP ActiveXModule::ARMTRIBSMODEL(BSTR pModel1,BSTR pModel2,BSTR pDiscModel,BSTR pFX1DiscVol,BSTR pFX2DiscVol,BSTR pIdx1Idx2Corr,BSTR pIdx1DiscIdxCorr,BSTR pIdx2DiscIdxCorr,BSTR pIdx1FxCorr,BSTR pIdx2FxCorr,int pQuantoFlag,BSTR *pRet)
{

try
{
	ARM_result C_result;

	_bstr_t bModel1(pModel1);
	CCString l_Model1 = bModel1;

	_bstr_t bModel2(pModel2);
	CCString l_Model2 = bModel2;
	
	_bstr_t bDiscModel(pDiscModel);
	CCString l_DiscModel = bDiscModel;
	
	_bstr_t bFX1DiscVol(pFX1DiscVol);
	CCString l_FX1DiscVol = bFX1DiscVol;

	_bstr_t bFX2DiscVol(pFX2DiscVol);
	CCString l_FX2DiscVol = bFX2DiscVol;

	_bstr_t bIdx1Idx2Corr(pIdx1Idx2Corr);
	CCString l_Idx1Idx2Corr = bIdx1Idx2Corr;
	
	_bstr_t bIdx1DiscIdxCorr(pIdx1DiscIdxCorr);
	CCString l_Idx1DiscIdxCorr = bIdx1DiscIdxCorr;

	_bstr_t bIdx2DiscIdxCorr(pIdx2DiscIdxCorr);
	CCString l_Idx2DiscIdxCorr = bIdx2DiscIdxCorr;

	_bstr_t bIdx1FxCorr(pIdx1FxCorr);
	CCString l_Idx1FxCorr = bIdx1FxCorr;

	_bstr_t bIdx2FxCorr(pIdx2FxCorr);
	CCString l_Idx2FxCorr = bIdx2FxCorr;

	CCString Res;
	CCString curClass = LOCAL_TRIBSMOD_CLASS;

	if(ARMLOCAL_TriBSModel (LocalGetNumObjectId (l_Model1),
							LocalGetNumObjectId (l_Model2),
							LocalGetNumObjectId (l_DiscModel),
							LocalGetNumObjectId (l_FX1DiscVol),
							LocalGetNumObjectId (l_FX2DiscVol),
							LocalGetNumObjectId (l_Idx1Idx2Corr),
							LocalGetNumObjectId (l_Idx1DiscIdxCorr),
							LocalGetNumObjectId (l_Idx2DiscIdxCorr),
							LocalGetNumObjectId (l_Idx1FxCorr),
							LocalGetNumObjectId (l_Idx2FxCorr),
							pQuantoFlag,
							C_result) == ARM_OK)
	{
		Res = LocalMakeObjectId (C_result.getLong (), curClass);
	}

	_bstr_t tmpChaine = Res;

	*pRet = SysAllocString(tmpChaine);

	return S_OK;
}
catch(...)
{
	return E_FAIL;
}

}


STDMETHODIMP ActiveXModule::ARMTRIBSDUAL(BSTR pModel1,BSTR pModel2,BSTR pDiscModel,BSTR pFX1DiscVol,BSTR pFX2DiscVol,BSTR pIdx1Idx2Corr,BSTR pIdx1DiscIdxCorr,BSTR pIdx2DiscIdxCorr,BSTR pIdx1FxCorr,BSTR pIdx2FxCorr,int pQuantoFlag,double pCorrelForAdj,int pWithslopeflag,BSTR *pRet)
{

try
{
	ARM_result C_result;

	_bstr_t bModel1(pModel1);
	CCString l_Model1 = bModel1;

	_bstr_t bModel2(pModel2);
	CCString l_Model2 = bModel2;
	
	_bstr_t bDiscModel(pDiscModel);
	CCString l_DiscModel = bDiscModel;
	
	_bstr_t bFX1DiscVol(pFX1DiscVol);
	CCString l_FX1DiscVol = bFX1DiscVol;

	_bstr_t bFX2DiscVol(pFX2DiscVol);
	CCString l_FX2DiscVol = bFX2DiscVol;

	_bstr_t bIdx1Idx2Corr(pIdx1Idx2Corr);
	CCString l_Idx1Idx2Corr = bIdx1Idx2Corr;
	
	_bstr_t bIdx1DiscIdxCorr(pIdx1DiscIdxCorr);
	CCString l_Idx1DiscIdxCorr = bIdx1DiscIdxCorr;

	_bstr_t bIdx2DiscIdxCorr(pIdx2DiscIdxCorr);
	CCString l_Idx2DiscIdxCorr = bIdx2DiscIdxCorr;

	_bstr_t bIdx1FxCorr(pIdx1FxCorr);
	CCString l_Idx1FxCorr = bIdx1FxCorr;

	_bstr_t bIdx2FxCorr(pIdx2FxCorr);
	CCString l_Idx2FxCorr = bIdx2FxCorr;

	CCString Res;
	CCString curClass = LOCAL_TRIBSMOD_CLASS;

	if(ARMLOCAL_TriBSDualModel (LocalGetNumObjectId (l_Model1),
								LocalGetNumObjectId (l_Model2),
								LocalGetNumObjectId (l_DiscModel),
								LocalGetNumObjectId (l_FX1DiscVol),
								LocalGetNumObjectId (l_FX2DiscVol),
								LocalGetNumObjectId (l_Idx1Idx2Corr),
								LocalGetNumObjectId (l_Idx1DiscIdxCorr),
								LocalGetNumObjectId (l_Idx2DiscIdxCorr),
								LocalGetNumObjectId (l_Idx1FxCorr),
								LocalGetNumObjectId (l_Idx2FxCorr),
								pQuantoFlag,
								pCorrelForAdj,
								pWithslopeflag,
								C_result) == ARM_OK)
	{
		Res = LocalMakeObjectId (C_result.getLong (), curClass);
	}

	_bstr_t tmpChaine = Res;

	*pRet = SysAllocString(tmpChaine);

	return S_OK;
}
catch(...)
{
	return E_FAIL;
}

}


STDMETHODIMP ActiveXModule::ARMBsSmiledModel(double pDate, double pSpot, BSTR pDividend, BSTR pDiscrate, BSTR pVolATM, BSTR pRo, BSTR pNu, BSTR pIsSABR, BSTR pBeta, BSTR *pRet)
{

try
{
	ARM_result C_result;

	_bstr_t bDiv(pDividend);
	CCString l_dividend = bDiv;
	_bstr_t bDiscrate(pDiscrate);
	CCString l_discrate = bDiscrate;
	_bstr_t bVol(pVolATM);
	CCString l_volATM = bVol;
	_bstr_t bRo(pRo);
	CCString l_ro = bRo;
	_bstr_t bNu(pNu);
	CCString l_nu= bNu;
	_bstr_t bIsSABR(pIsSABR);
	CCString l_issabr= bIsSABR;
	_bstr_t bBeta(pBeta);
	CCString l_beta= bBeta;

	long isSabrId;

	VECTOR<double> maturities;
	VECTOR<double> ro;
	VECTOR<double> nu;
	VECTOR<double> beta;

	if((isSabrId = ARM_ConvSmiledModelFlag (l_issabr, C_result)) == ARM_DEFAULT_ERR)
		return S_FALSE;

	long betaId = ARM_NULL_OBJECT;
	if (strcmp((const char*)l_beta,"NULL") != 0)
		betaId = LocalGetNumObjectId(l_beta);

	CCString Res;
	CCString curClass = LOCAL_BSMODEL_CLASS;

	if(ARMLOCAL_bssmiledmodel(pDate,
							  pSpot,
							  1,
							  LocalGetNumObjectId(l_dividend),
							  1,
							  LocalGetNumObjectId(l_discrate),
							  1,
							  LocalGetNumObjectId(l_volATM),
							  K_YIELD,
							  maturities,
							  1,
							  LocalGetNumObjectId(l_ro),
							  ro,
							  1,
							  LocalGetNumObjectId(l_nu),
							  nu,
							  isSabrId,
							  1,
							  betaId,
							  beta,
							  0.0,
							  1,
							  ARM_NULL_OBJECT,
							  1,
							  ARM_NULL_OBJECT,
							  ARM_NULL_OBJECT,
							  1,
							  C_result) == ARM_OK)
	{
		Res = LocalMakeObjectId (C_result.getLong (), curClass);
	}

	_bstr_t tmpChaine = Res;

	*pRet = SysAllocString(tmpChaine);

	return S_OK;
}
catch(...)
{
	return E_FAIL;
}

}



STDMETHODIMP ActiveXModule::ARMGetVolFromSummit(VARIANT *pIndex, VARIANT* pCcy, VARIANT* pCvName, VARIANT* pAsOf, VARIANT* pType, VARIANT* pMatIndex, VARIANT* pImpOrHist, BSTR pIndexId,VARIANT *pRet)
{
try
{
	ARM_result C_result;

	CCString l_index;
	CCString l_ccy;
	CCString l_cvname;
	CCString l_type;
	CCString l_matIndex;
	CCString l_imporhist;

	_variant_t valDate( pAsOf );
	valDate.ChangeType(VT_R8);
	
	if(VARIANT2CCString (*pIndex, l_index) != S_OK)
		return S_FALSE;

	if(VARIANT2CCString (*pCcy, l_ccy) != S_OK)
		return S_FALSE;

	if(VARIANT2CCString (*pCvName, l_cvname) != S_OK)
		return S_FALSE;

	if(VARIANT2CCString (*pType, l_type) != S_OK)
		return S_FALSE;

	if(VARIANT2CCString (*pMatIndex, l_matIndex) != S_OK)
		return S_FALSE;

	if(VARIANT2CCString (*pImpOrHist, l_imporhist) != S_OK)
		return S_FALSE;

	CCString Res;
	CCString curClass = LOCAL_VOL_CURVE_LIN_CLASS;

	std::string indexId; VariantTools::convert(pIndexId,indexId); 
	long indexId_= LocalPersistent::get().getObjectId(indexId);

	if(ARMLOCAL_GetVolFromSummit (l_index, l_ccy, l_cvname, valDate, l_type, l_matIndex, l_imporhist, indexId_,C_result) == ARM_OK)
	{
		Res = LocalMakeObjectId (C_result.getLong (), curClass);
	}

	_variant_t wrap_Res;
	wrap_Res.Attach(*pRet);
	wrap_Res.SetString(Res);
	*pRet=wrap_Res.Detach();

	pRet->vt=VT_BSTR;

	return S_OK;
}
catch(...)
{
	return E_FAIL;
}
}


STDMETHODIMP ActiveXModule::ARMGetFXVolFromSummit(BSTR pCcy1, BSTR pCcy2, double pDate, BSTR pCvName, BSTR pType, BSTR *pRet)
{
try
{
	ARM_result C_result;

	CCString Res;

	_bstr_t bCcy1(pCcy1);
	CCString l_Ccy1 = bCcy1;

	_bstr_t bCcy2(pCcy2);
	CCString l_Ccy2 = bCcy2;

	_bstr_t bCvName(pCvName);
	CCString l_CvName = bCvName;

	_bstr_t bType(pType);
	CCString l_Type = bType;

	CCString curClass;

	if ((strcmp((const char*)l_Type,"CSMILE")==0) || (strcmp((const char*)l_Type,"CFXSPI")==0))
		curClass = LOCAL_VOL_CUBE_CLASS;
	else
		curClass = LOCAL_VOL_CURVE_LIN_CLASS;

	if(ARMLOCAL_GetFXVolFromSummit (l_Ccy1, l_Ccy2, pDate, l_CvName, "FXVOL", l_Type, C_result) == ARM_OK)
	{
		Res = LocalMakeObjectId (C_result.getLong (), curClass);
	}

	_bstr_t tmpChaine = Res;

	*pRet = SysAllocString(tmpChaine);

	return S_OK;
}
catch(...)
{
	return E_FAIL;
}
}


STDMETHODIMP ActiveXModule::ARMGetFXCorrelFromSummit(BSTR pCcy1, BSTR pIndex, BSTR pCcy2, double pDate, BSTR pCvName, VARIANT* pTenors, BSTR *pRet)
{
try
{
	ARM_result C_result;

	CCString Res;

	_bstr_t bCcy1(pCcy1);
	CCString l_Ccy1 = bCcy1;

	_bstr_t bIndex(pIndex);
	CCString l_Index = bIndex;

	_bstr_t bCcy2(pCcy2);
	CCString l_Ccy2 = bCcy2;

	_bstr_t bCvName(pCvName);
	CCString l_CvName = bCvName;

	long tenorsSize;
	VECTOR<CCString> vTenors;

	if(VARIANT2VECTORCCSTRING (*pTenors,vTenors,tenorsSize) != S_OK)
		return S_FALSE;

	CCString curClass = LOCAL_VOL_CURVE_LIN_CLASS;

	if(ARMLOCAL_GetFXCorrelFromSummit (l_Ccy1, l_Index, l_Ccy2, pDate, l_CvName, vTenors, C_result) == ARM_OK)
	{
		Res = LocalMakeObjectId (C_result.getLong (), curClass);
	}

	_bstr_t tmpChaine = Res;

	*pRet = SysAllocString(tmpChaine);

	return S_OK;
}
catch(...)
{
	return E_FAIL;
}
}


STDMETHODIMP ActiveXModule::ARMGetCorrelFromSummit(BSTR pCcy1, BSTR pIndex1, BSTR pCcy2, BSTR pIndex2, double pDate, BSTR pCvName, BSTR *pRet)
{
try
{
	ARM_result C_result;

	CCString Res;

	_bstr_t bCcy1(pCcy1);
	CCString l_Ccy1 = bCcy1;

	_bstr_t bIndex1(pIndex1);
	CCString l_Index1 = bIndex1;

	_bstr_t bCcy2(pCcy2);
	CCString l_Ccy2 = bCcy2;

	_bstr_t bIndex2(pIndex2);
	CCString l_Index2 = bIndex2;

	_bstr_t bCvName(pCvName);
	CCString l_CvName = bCvName;

	CCString curClass = LOCAL_VOL_CURVE_LIN_CLASS;

	if(ARMLOCAL_GetCorrelFromSummit (l_Ccy1, l_Index1, l_Ccy2, l_Index2, pDate, l_CvName, C_result) == ARM_OK)
	{
		Res = LocalMakeObjectId (C_result.getLong (), curClass);
	}

	_bstr_t tmpChaine = Res;

	*pRet = SysAllocString(tmpChaine);

	return S_OK;
}
catch(...)
{
	return E_FAIL;
}
}


STDMETHODIMP ActiveXModule::ARMVolFlat(double pVol, double pDate, BSTR pCcy, BSTR* pRet)
{
try
{
	ARM_result C_result;

	_bstr_t bCcy(pCcy);
	CCString l_ccy = bCcy;

	CCString Res;
	CCString curClass = LOCAL_VOL_FLAT_CLASS;

	if(ARMLOCAL_volflat (pVol, pDate, l_ccy, C_result) == ARM_OK)
	{
		Res = LocalMakeObjectId (C_result.getLong (), curClass);
	}

	_bstr_t tmpChaine = Res;

	*pRet = SysAllocString(tmpChaine);

	return S_OK;
}
catch(...)
{
	return E_FAIL;
}
}


STDMETHODIMP ActiveXModule::ARMVolCube(BSTR pATMVol,VARIANT *pSmileCurveIds,VARIANT *pTenors,BSTR pVolType,BSTR pRefObj,BSTR *pRet)
{
try
{
	ARM_result C_result;

	CCString Res;

	_bstr_t bATMVol(pATMVol);
	CCString l_ATMVol = bATMVol;

	_bstr_t bVolType(pVolType);
	CCString l_volType = bVolType;
	long volTypeId;

	long tenorsSize;
	VECTOR<double> vTenors;

	long smileSize;
	VECTOR<CCString> vSmiles;
	VECTOR<long> vSmilesId;

	if(VARIANT2VECTORDOUBLE (*pTenors,vTenors,tenorsSize) != S_OK)
		return S_FALSE;

	if(VARIANT2VECTORCCSTRING (*pSmileCurveIds,vSmiles,smileSize) != S_OK)
		return S_FALSE;

	for (int i= 0; i < smileSize; i++)
		vSmilesId.push_back(LocalGetNumObjectId(vSmiles[i]));

	if((volTypeId = ARM_ConvVolType(l_volType, C_result)) == ARM_DEFAULT_ERR)
		return S_FALSE;

	_bstr_t bRefObj(pRefObj);
	CCString l_refXL = bRefObj;

	CCString curClass = LOCAL_VOL_CUBE_CLASS;
	CCString stringId = GetEnvVar (l_refXL);

	long retCode;
	long objId;
	CCString prevClass;

	if (!stringId)
	{

		retCode = ARMLOCAL_VolCube(LocalGetNumObjectId(l_ATMVol),vSmilesId,vTenors,volTypeId, 1, C_result);

		if ( retCode == ARM_OK )
		{
			objId = C_result.getLong ();

			if (l_refXL.GetLen()>0)
						LocalSetCurCellEnvValue (curClass, objId,l_refXL); 


			stringId = LocalMakeObjectId (objId, curClass);
		}
	}
	else
	{
		prevClass = LocalGetStringObjectClass (stringId);
		
		objId = LocalGetNumObjectId (stringId);
			
		if ( curClass == prevClass )
		{
			retCode = ARMLOCAL_VolCube(LocalGetNumObjectId(l_ATMVol),vSmilesId,vTenors,volTypeId, 1,C_result,objId);
		}
		else
		{
			FreeCurCellContent (stringId);
			retCode = ARMLOCAL_VolCube(LocalGetNumObjectId(l_ATMVol),vSmilesId,vTenors,volTypeId,1,C_result);
		
			if ( retCode == ARM_OK )
			{
				objId = C_result.getLong ();

				LocalSetCurCellEnvValue (curClass, objId,l_refXL); 

				stringId = LocalMakeObjectId (objId, curClass);
			}
		}
	}

	_bstr_t tmpChaine = stringId;

	*pRet = SysAllocString(tmpChaine);

	return S_OK;
}
catch(...)
{
	return E_FAIL;
}
}


STDMETHODIMP ActiveXModule::ARMZcFlat(double pZc, double pDate, BSTR pCcy, BSTR* pRet)
{
try
{
	ARM_result C_result;

	_bstr_t bCcy(pCcy);
	CCString l_ccy = bCcy;

	CCString Res;
	CCString curClass = LOCAL_ZERO_CURVE_FLAT_CLASS;

	if(ARMLOCAL_zcflat (pZc, pDate, l_ccy, C_result) == ARM_OK)
	{
		Res = LocalMakeObjectId (C_result.getLong (), curClass);
	}

	_bstr_t tmpChaine = Res;

	*pRet = SysAllocString(tmpChaine);

	return S_OK;
}
catch(...)
{
	return E_FAIL;
}
}


STDMETHODIMP ActiveXModule::ARMBsModel(double pDate, double pSpot, BSTR pDividend, BSTR pDiscrate, BSTR pVol, BSTR pTypeStk, BSTR *pRet)
{

try
{
	ARM_result C_result;

	_bstr_t bDiv(pDividend);
	CCString l_dividend = bDiv;
	_bstr_t bDiscrate(pDiscrate);
	CCString l_discrate = bDiscrate;
	_bstr_t bVol(pVol);
	CCString l_vol = bVol;
	
	_bstr_t bTypeStk(pTypeStk);
	CCString l_typeStk = bTypeStk;
	long typeStkId;

	if((typeStkId = ARM_ConvPriceYield (l_typeStk, C_result)) == ARM_DEFAULT_ERR)
		return S_FALSE;

	CCString Res;
	CCString curClass = LOCAL_BSMODEL_CLASS;

	if(ARMLOCAL_bsmodel(pDate,
						pSpot,
						1,
						LocalGetNumObjectId(l_dividend),
						1,
						LocalGetNumObjectId(l_discrate),
						1,
						LocalGetNumObjectId(l_vol),
						typeStkId,
						C_result) == ARM_OK)
	{
		Res = LocalMakeObjectId (C_result.getLong (), curClass);
	}

	_bstr_t tmpChaine = Res;

	*pRet = SysAllocString(tmpChaine);

	return S_OK;
}
catch(...)
{
	return E_FAIL;
}

}


STDMETHODIMP ActiveXModule::ARMBsSlModel(double pDate, BSTR pZc, BSTR pVolSpreadLock, BSTR pCvCapVol, BSTR pCvIndexVol, BSTR *pRet)
{

try
{
	ARM_result C_result;

	_bstr_t bZc(pZc);
	CCString l_zc = bZc;
	_bstr_t bVolSpreadLock(pVolSpreadLock);
	CCString l_volspreadlock = bVolSpreadLock;
	_bstr_t bCvCapVol(pCvCapVol);
	CCString l_cvcapvol = bCvCapVol;
	_bstr_t bCvIndexVol(pCvIndexVol);
	CCString l_cvindexvol = bCvIndexVol;

	long indexvol_type;
	long indexvol;
	long indexvoldefault=0;

	if (l_cvindexvol=="DEFAULT")
	{
		indexvol_type=0L;
		indexvol=indexvoldefault;
	}
	else
	{
		indexvol_type=1L;
		indexvol=LocalGetNumObjectId(l_cvindexvol);
	}
	
	CCString Res;
	CCString curClass = LOCAL_BSMODEL_CLASS;

	if(ARMLOCAL_bsslmodel(pDate,
						  1,
						  LocalGetNumObjectId(l_zc),
						  1,
						  LocalGetNumObjectId(l_volspreadlock),
						  1,
						  LocalGetNumObjectId(l_cvcapvol),
						  indexvol_type,
						  indexvol,
						  C_result) == ARM_OK)
	{
		Res = LocalMakeObjectId (C_result.getLong (), curClass);
	}

	_bstr_t tmpChaine = Res;

	*pRet = SysAllocString(tmpChaine);

	return S_OK;
}
catch(...)
{
	return E_FAIL;
}

}


STDMETHODIMP ActiveXModule::ARMSwitchToETK()
{
	try {
		switchToETK();
	} catch(...) { return E_FAIL; }

	return S_OK;
}

STDMETHODIMP ActiveXModule::ARMSwitchToWSETK()
{
	try{
		switchToWSETK();
	} catch(...) { return E_FAIL; }

	return S_OK;
}

STDMETHODIMP ActiveXModule::ARMShutdownETK()
{
	try {
		shutdown_etoolkit();
	} catch(...) { return E_FAIL; }

	return S_OK;
}


STDMETHODIMP ActiveXModule::ARMSwitchToFLATFILE()
{
	try {
		switchToFLATFILE();
	} catch(...) { return E_FAIL; }

	return S_OK;
}


STDMETHODIMP ActiveXModule::ARMInfocentreConnect()
{
	try {
		switchToETK();

	connection_etoolkit(SUMMIT_INFOC_CONNEXION_USERNAME,
						SUMMIT_INFOC_CONNEXION_PASSWD,
						SUMMIT_INFOC_CONNEXION_CONTEXT,
						SUMMIT_INFOC_IT_CONFIG_DOMAINSDIR,
						SUMMIT_INFOC_IT_DOMAIN_NAME);
	} catch(...) { return E_FAIL; }
	return S_OK;
}


STDMETHODIMP ActiveXModule::ARMBaseReplicationConnect()
{
	try {
		switchToETK();

	connection_etoolkit(SUMMIT_REPLI_CONNEXION_USERNAME,
						SUMMIT_REPLI_CONNEXION_PASSWD,
						SUMMIT_REPLI_CONNEXION_CONTEXT,
						SUMMIT_REPLI_IT_CONFIG_DOMAINSDIR,
						SUMMIT_REPLI_IT_DOMAIN_NAME);
	} catch(...) { return E_FAIL; }
	return S_OK;
}

STDMETHODIMP ActiveXModule::ARMZCLINT(VARIANT* pMatu, VARIANT* pRate, BSTR pMeth, double pDate, BSTR pCurrency, BSTR pInterpMeth,double pMatuAreDoubles,VARIANT *pSTerms, BSTR *pRet)
{
	// TODO: Add your implementation code here
try
{
	ARM_result C_result;
	ARM_CRV_TERMS sTerms;

	VECTOR<double> vMatu;
	long matuSize;
	VECTOR<double> vRate;
	long rateSize;
	VECTOR<CCString> vTerms;
	long termsSize;
	CCString defaultTerms("");


	if(VARIANT2VECTORDOUBLE (*pMatu,vMatu,matuSize) != S_OK)
		return S_FALSE;

	if(VARIANT2VECTORDOUBLE (*pRate,vRate,rateSize) != S_OK)
		return S_FALSE;



	try{
		if(VARIANT2CCString (*pSTerms,defaultTerms) != S_OK)
		{
			return S_FALSE;
		}
		else
		{
			vTerms.clear();
		}

	}
	catch (...)
	{
		if(VARIANT2VECTORCCSTRING (*pSTerms,vTerms,termsSize) != S_OK)
		{
			return S_FALSE;	
		}
		else
		{
			for (int i = 0; i < vTerms.size(); i++)
			{
				strcpy(sTerms[i], vTerms[i]);
			}
		}
		
	}

	_bstr_t bMeth(pMeth);
	CCString l_meth = bMeth;
	long methId;

	if((methId = ARM_ConvForwardYieldMethod (l_meth, C_result)) == ARM_DEFAULT_ERR)
		return S_FALSE;

	_bstr_t bCurrency(pCurrency);
	CCString l_currency = bCurrency;

	_bstr_t bInterpMeth(pInterpMeth);
	CCString l_interpmeth = bInterpMeth;
	long interpMethId;
	
	if((interpMethId = ARM_ConvInterpMethod (l_interpmeth, C_result)) == ARM_DEFAULT_ERR)
		return S_FALSE;

		CCString Res;
	CCString curClass = LOCAL_ZERO_CURVE_LIN_CLASS;

	int matuAreDoubles = pMatuAreDoubles;

	if(ARMLOCAL_zclint(vMatu,vRate,methId,pDate,l_currency,interpMethId, matuAreDoubles, sTerms, C_result) == ARM_OK)
	{
		Res = LocalMakeObjectId (C_result.getLong (), curClass);
	}

	_bstr_t tmpChaine = Res;

	*pRet = SysAllocString(tmpChaine);

	return S_OK;
}
catch(...)
{
	return E_FAIL;
}
}


STDMETHODIMP ActiveXModule::ARM_zcspreaded(BSTR zcSprId,
										BSTR zcInitId,
										double date,
										BSTR MMFreq,
										BSTR SwapFreq,
										BSTR ccyId,
										BSTR *pRet)
{

try
{
	ARM_result C_result;

	_bstr_t bfreq(MMFreq);
	CCString C_MMFreq = bfreq;

	_bstr_t bSwapFreq(SwapFreq);
	CCString C_SwapFreq = bSwapFreq;

	_bstr_t bccyId(ccyId);
	CCString C_ccyId = bccyId;

	_bstr_t bzcSprId(zcSprId);
	CCString C_zcSprId = bzcSprId;

	_bstr_t bzcInitId(zcInitId);
	CCString C_zcInitId = bzcInitId;

	bool ccyIsObject = false;

	if ((C_ccyId.GetLen() > 3)
		&&
		!(C_ccyId == "DEFAULT"))
		ccyIsObject = true;

	long mmFreqId=0;
	long swapFreqId=0;

	if((mmFreqId = ARM_ConvFrequency (C_MMFreq, C_result)) == ARM_DEFAULT_ERR)
		return S_FALSE;

	if((swapFreqId = ARM_ConvFrequency (C_SwapFreq, C_result)) == ARM_DEFAULT_ERR)
		return S_FALSE;

	CCString Res;
	CCString curClass = LOCAL_ZERO_CURVE_LIN_CLASS;

	if(ARMLOCAL_zcspreaded (LocalGetNumObjectId(C_zcSprId),
							LocalGetNumObjectId(C_zcInitId),
							date,
						    mmFreqId,
							swapFreqId,
							ccyIsObject,
							C_ccyId,
						    C_result) == ARM_OK)
	{
		Res = LocalMakeObjectId (C_result.getLong (), curClass);
	}
	else
	{
		return createErrorInfo("",CCSTringToSTLString(C_result.getMsg())); 
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
	return S_OK;

}


STDMETHODIMP ActiveXModule::ARMCreateZCSwapInt(double pDate,VARIANT *pMatu,VARIANT *pRate,BSTR pMMVsFut,
										   BSTR pSwapVsFut,BSTR pRaw,BSTR pInterp,BSTR pCcy,
										   BSTR pRefObj,BSTR *pRet)
{
try
{
	ARM_result C_result;

	CCString Res;

	long rateSize;

	_bstr_t bCcy(pCcy);
	CCString l_currency = bCcy;
	
	_bstr_t bMMVsFut(pMMVsFut);
	CCString C_MMVsFut = bMMVsFut;
	long MMVsFutId;

	_bstr_t bSwapVsFut(pSwapVsFut);
	CCString C_SwapVsFut = bSwapVsFut;
	long SwapVsFutId;

	_bstr_t braw(pRaw);
	CCString C_raw = braw;
	long RawId;

	_bstr_t binterp(pInterp);
	CCString C_interp = binterp;
	long InterpId;

	_bstr_t bRefObj(pRefObj);
	CCString l_refXL = bRefObj;

	VECTOR<CCString> vMatu;
	VECTOR<double> vRate;
	VECTOR<CCString> vMatuRate;

	long sizeMatu;

	if(VARIANT2VECTORCCSTRING (*pMatu,vMatu,sizeMatu) != S_OK)
		return S_FALSE;

	if(VARIANT2VECTORDOUBLE (*pRate,vRate,rateSize) != S_OK)
		return S_FALSE;

	if (sizeMatu!=rateSize)
		return S_FALSE;

	int j = 0;

	if(l_currency== "DEFAULT")
	{
		ARM_result currencyres;
		ARMLOCAL_ARM_GetDefaultCurrency (currencyres);
		if(currencyres.getRetCode () != ARM_OK)
		{
			return S_FALSE;
		}
		else
		{
			l_currency = currencyres.getString ();
		}
	}

	if((MMVsFutId = ARM_ConvMktType (C_MMVsFut, C_result)) == ARM_DEFAULT_ERR)
			return S_FALSE;
	
	if((SwapVsFutId = ARM_ConvMktType (C_SwapVsFut, C_result)) == ARM_DEFAULT_ERR)
		return S_FALSE;

	RawId = ARM_ConvCvMethod (C_raw);

	if((InterpId = ARM_ConvInterpMethod (C_interp, C_result)) == ARM_DEFAULT_ERR)
		return S_FALSE;

	long retCode;
	long objId;
	CCString prevClass;

	CCString curClass = LOCAL_ZERO_CURVE_LIN_CLASS;
	CCString stringId = GetEnvVar (l_refXL);
	
	if(!stringId)
	{
		retCode = ARMLOCAL_CreateZCSwapInt(pDate, vMatu, vRate, (long)MMVsFutId, 
									   (long)SwapVsFutId, (long)RawId,
									   (long)InterpId, l_currency, K_DEF_FREQ,KNOBASE,
									   C_result);
		if(retCode == ARM_OK)
		{
			objId = C_result.getLong ();

			if (l_refXL.GetLen()>0)
						LocalSetCurCellEnvValue (curClass, objId,l_refXL); 


			stringId = LocalMakeObjectId (objId, curClass);
		}
	}
	else
	{
		prevClass = LocalGetStringObjectClass (stringId);
		
		objId = LocalGetNumObjectId (stringId);
			
		if(curClass == prevClass)
		{
			retCode = ARMLOCAL_CreateZCSwapInt(pDate, vMatu, vRate, (long)MMVsFutId, 
									   (long)SwapVsFutId, (long)RawId,
									   (long)InterpId, l_currency, K_DEF_FREQ,KNOBASE,
									   C_result,objId);

		
		}
		else
		{
			FreeCurCellContent (stringId);

			retCode = ARMLOCAL_CreateZCSwapInt(pDate, vMatu, vRate, (long)MMVsFutId, 
									   (long)SwapVsFutId, (long)RawId,
									   (long)InterpId, l_currency, K_DEF_FREQ,KNOBASE,
									   C_result);
		
			if(retCode == ARM_OK)
			{
				objId = C_result.getLong ();

				LocalSetCurCellEnvValue (curClass, objId,l_refXL); 

				stringId = LocalMakeObjectId (objId, curClass);
			}
		}
	}

	_bstr_t tmpChaine = stringId;

	*pRet = SysAllocString(tmpChaine);

	return S_OK;
}
catch(...)
{
	return E_FAIL;
}
}


STDMETHODIMP ActiveXModule::ARMGetInitialCurveFromSummit( BSTR pIndex, BSTR pCurrency,  
													  BSTR pCvName,  double pDate, BSTR pAdjOrNot,  
													  VARIANT *pRetMat, VARIANT *pRetRate)
{
	// TODO: Add your implementation code here
try
{
	ARM_result C_result;

	_bstr_t bIndex(pIndex);
	CCString l_index = bIndex;

	_bstr_t bCurrency(pCurrency);
	CCString l_currency = bCurrency;

	_bstr_t bCvName(pCvName);
	CCString l_cvname = bCvName;

	_bstr_t bAdjOrNot(pAdjOrNot);
	CCString l_adjornot = bAdjOrNot;
	
	long adjOrNotId;

	if((adjOrNotId = ARM_ConvYesOrNo (l_adjornot, C_result)) == ARM_DEFAULT_ERR)
	{
		return S_FALSE;
	}

	VECTOR<CCString> matu;
	VECTOR<double> rate;

	long retCode;

	retCode = ARMLOCAL_GetInitialCurveFromSummit (l_index,l_currency,l_cvname,pDate,adjOrNotId,&matu,&rate,C_result);

	if (C_result.getLong() == FFRETRIEVER)
	{
		if (retCode == ARM_KO)
			return E_FAIL;

		retCode = LocalExtractCurveFromFileMO (C_result.getString(), matu, rate, adjOrNotId);
		matu.erase(matu.begin() + matu.size() - 1);
		rate.erase(rate.begin() + rate.size() - 1);
	}

	if (retCode == ARM_KO)
		return E_FAIL;

	long Res;
	
	Res = VECTORCCSTRING2VARIANT(matu, pRetMat);
	
	Res = VECTORDOUBLE2VARIANT(rate, pRetRate);


	return S_OK;
}
catch(...)
{
	return E_FAIL;
}
}



STDMETHODIMP ActiveXModule::ARMGetInitialCurveFromCalypso(BSTR pIndex, 
												BSTR pCurrency, 
												BSTR pTerm,
												BSTR pPricingEnv, 
												double pDate, 
												BSTR pForceCurveName,
												BSTR pXmlFileName,
												BSTR pDoAdj, 
												VARIANT *pRetMat,
												VARIANT *pRetRate) 
													  
{
	// TODO: Add your implementation code here
	try
	{
		std::string index ; VariantTools::convert(pIndex,index) ;
		std::string ccy; VariantTools::convert(pCurrency,ccy) ;
		std::string term; VariantTools::convert(pTerm,term) ;
		std::string pricingEnv; VariantTools::convert(pPricingEnv,pricingEnv) ;
		std::string forceCurveName; VariantTools::convert(pForceCurveName,forceCurveName) ;
		std::string xmlFileName; VariantTools::convert(pXmlFileName,xmlFileName) ;
		std::string doAdjStr; VariantTools::convert(pDoAdj,doAdjStr) ;
		ARM_Date AsOf ; VariantTools::convertXLDate(pDate,AsOf); 

		bool doAdj; 
		{
			ARM_result res; 
			doAdj = ARM_ConvYesOrNo (doAdjStr) ; 
		}

		std::vector<std::string> matus; 
		std::vector<double> yields; 
		ARMLOCAL_GetInitialCurveFromCalypso(AsOf,pricingEnv,index,ccy,term,forceCurveName,xmlFileName,doAdj,matus,yields); 

		VariantTools::convert(matus,*pRetMat); 
		VariantTools::convert(yields,*pRetRate); 

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
	return S_OK;
}



STDMETHODIMP ActiveXModule::ARMGetInitialVolFromSummit(BSTR pIndex, BSTR pCurrency,
												   BSTR pCvName,  double pDate, BSTR pType,
												   BSTR pMatIndex, VARIANT *pRetMat, VARIANT *pRetTenor,
												   VARIANT *pRetVol)
{
	// TODO: Add your implementation code here
try
{
	ARM_result C_result;

	_bstr_t bIndex(pIndex);
	CCString l_index = bIndex;

	_bstr_t bCurrency(pCurrency);
	CCString l_currency = bCurrency;

	_bstr_t bCvName(pCvName);
	CCString l_cvname = bCvName;

	_bstr_t bType(pType);
	CCString l_type = bType;

	_bstr_t bMatIndex(pMatIndex);
	CCString l_MatIndex = bMatIndex;

	VECTOR<CCString> matu;
	VECTOR<CCString> tenor;
	VECTOR<double> vol;

	long retCode;

	retCode = ARMLOCAL_GetInitialVolFromSummit (l_index,l_currency,l_cvname,pDate,l_type,l_MatIndex,&matu,&tenor,&vol,C_result);

	if (retCode == ARM_KO)
		return E_FAIL;

	long Res;
	
	Res = VECTORCCSTRING2VARIANT(matu, pRetMat);
	
	Res = VECTORCCSTRING2VARIANT(tenor, pRetTenor);

	Res = VECTORDOUBLE2VARIANT(vol, pRetVol);

	return S_OK;
}
catch(...)
{
	return E_FAIL;
}
}

STDMETHODIMP ActiveXModule::ARM_Credit_CreateBasketCorrelMkDataFromCalypso(BSTR pricingEnv_,
																		   double date_,
																		   BSTR forceCurveName_,
																		   BSTR Ccy_ ,
																		   BSTR xmlFileName_,
																		   BSTR IndexId_,
																		   BSTR *pRet)
{
	try
	{
		std::string pricingEnv ; VariantTools::convert(pricingEnv_,pricingEnv); 
		std::string forceCurveName ; VariantTools::convert(forceCurveName_,forceCurveName); 
		std::string xmlFileName ; VariantTools::convert(xmlFileName_,xmlFileName);
		std::string Ccy ; VariantTools::convert(Ccy_,Ccy);
		std::string IndexId ; VariantTools::convert(IndexId_,IndexId);

		ARM_Date date ; VariantTools::convertXLDate(date_,date); 

		long indexId_= LocalPersistent::get().getObjectId(IndexId);

		long id= ICMLOCAL_CreateBasketCorrelMkDataFromCalypso(pricingEnv,date,forceCurveName,Ccy,xmlFileName,indexId_) ;
		std::string tem = LocalPersistent::get().getStringId(id,LOCAL_VOL_CURVE_LIN_CLASS); 
		VariantTools::convert(tem,*pRet); 
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
	return S_OK;

}

STDMETHODIMP ActiveXModule::ARM_Credit_GetBasketCorrelMkDataFromCalypso(BSTR pricingEnv_,
																  double date_,
																  BSTR forceCurveName_, 
																  BSTR xmlFileName_,
																  VARIANT *pRetMat, 
																  VARIANT *pRetTenor,
																  VARIANT* pRetVol)
{
	// TODO: Add your implementation code here
	try
	{
		std::string pricingEnv ; VariantTools::convert(pricingEnv_,pricingEnv); 
		std::string forceCurveName ; VariantTools::convert(forceCurveName_,forceCurveName); 
		std::string xmlFileName ; VariantTools::convert(xmlFileName_,xmlFileName); 
		ARM_Date date ; VariantTools::convertXLDate(date_,date); 

		std::vector<std::string> matus; 
		std::vector<double> attach; 
		ICM_QMatrix<double> correls; 
		
		ICMLOCAL_GetBasketCorrelMkDataFromCalypso(pricingEnv,date,forceCurveName,xmlFileName,matus,attach,correls) ;
		VariantTools::convert(matus,*pRetMat);
		VariantTools::convert(attach,*pRetTenor);
		VariantTools::convert(correls,*pRetVol);
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
	return S_OK;
}



STDMETHODIMP ActiveXModule::ARMGetInitialFXVolFromSummit(BSTR pCcy1,BSTR pCcy2,double pDate,BSTR pCvName,
													 BSTR pImpOrHist,BSTR pVolType,VARIANT *pRetMat,
													 VARIANT *pRetTenor,VARIANT *pRetVol)
{
	// TODO: Add your implementation code here
try
{
	ARM_result C_result;

	_bstr_t bCcy1(pCcy1);
	CCString l_Ccy1 = bCcy1;

	_bstr_t bCcy2(pCcy2);
	CCString l_Ccy2 = bCcy2;

	_bstr_t bCvName(pCvName);
	CCString l_cvname = bCvName;

	_bstr_t bImpOrHist(pImpOrHist);
	CCString l_ImpOrHist = bImpOrHist;

	_bstr_t bVolType(pVolType);
	CCString l_VolType = bVolType;

	VECTOR<CCString> matu;
	VECTOR<double> tenor;
	VECTOR<double> vol;

	long retCode;

	retCode = ARMLOCAL_GetInitialFXVolFromSummit (l_Ccy1,l_Ccy2,pDate,l_cvname,l_ImpOrHist,l_VolType,&matu,&tenor,&vol,C_result);

	if (retCode == ARM_KO)
		return E_FAIL;

	long Res;
	
	Res = VECTORCCSTRING2VARIANT(matu, pRetMat);
	
	Res = VECTORDOUBLE2VARIANT(tenor, pRetTenor);

	Res = VECTORDOUBLE2VARIANT(vol, pRetVol);

	return S_OK;
}
catch(...)
{
	return E_FAIL;
}
}


STDMETHODIMP ActiveXModule::ARMTHREEMONTHFUT(BSTR pDelivery,
										 long pMarket,
										 BSTR pCcy,
										 BSTR pRefObj,
										 BSTR *pRet)
{
try
{
	ARM_result C_result;

	_bstr_t bDelivery(pDelivery);
	CCString l_delivery = bDelivery;

	_bstr_t bCcy(pCcy);
	CCString l_ccy = bCcy;

	long retCode;
	long objId;
	CCString prevClass;
	
	_bstr_t bRefObj(pRefObj);
	CCString l_refXL = bRefObj;

	CCString curClass = LOCAL_IRFUT_CLASS;
	CCString stringId = GetEnvVar (l_refXL);
	
	if(!stringId)
	{
		retCode = ARMLOCAL_THREE_MONTH_FUT (l_delivery, pMarket, l_ccy, C_result);

		if(retCode == ARM_OK)
		{
			objId = C_result.getLong ();

			if (l_refXL.GetLen()>0)
				LocalSetCurCellEnvValue (curClass, objId,l_refXL); 

			stringId = LocalMakeObjectId (objId, curClass);
		}
	}
	else
	{
		prevClass = LocalGetStringObjectClass (stringId);
		
		objId = LocalGetNumObjectId (stringId);
			
		if(curClass == prevClass)
		{
			retCode = ARMLOCAL_THREE_MONTH_FUT (l_delivery, pMarket, l_ccy, C_result, objId);
		}
		else
		{
			FreeCurCellContent (stringId);
			retCode = ARMLOCAL_THREE_MONTH_FUT (l_delivery, pMarket, l_ccy, C_result);
		
			if(retCode == ARM_OK)
			{
				objId = C_result.getLong ();
			
				LocalSetCurCellEnvValue (curClass, objId,l_refXL); 

				stringId = LocalMakeObjectId (objId, curClass);
			}
		}
	}

	_bstr_t tmpChaine = stringId;

	*pRet = SysAllocString(tmpChaine);

	return S_OK;
}
catch(...)
{
	return E_FAIL;
}
}


STDMETHODIMP ActiveXModule::ARMFutPibor(BSTR pDelivery,BSTR *pRet)
{
try
{
	ARM_result C_result;

	_bstr_t bDelivery(pDelivery);
	CCString l_delivery = bDelivery;

	CCString prevClass;

	CCString curClass = LOCAL_IRFUT_CLASS;
	CCString Res;

	if(ARMLOCAL_THREE_MONTH_FUT (l_delivery, 0, "FRF", C_result) == ARM_OK)
	{
		Res = LocalMakeObjectId (C_result.getLong (), curClass);
	}

	_bstr_t tmpChaine = Res;

	*pRet = SysAllocString(tmpChaine);

	return S_OK;
}
catch(...)
{
	return E_FAIL;
}
}


STDMETHODIMP ActiveXModule::ARMIRFUT(double pDelivery,BSTR pIdUnderlying,BSTR pRefObj,BSTR *pRet)
{
try
{
	ARM_result C_result;

	_bstr_t bIdUnderlying(pIdUnderlying);
	CCString lIdUnderlying = bIdUnderlying;

	_bstr_t bRefObj(pRefObj);
	CCString l_refXL = bRefObj;

	long retCode;
	long objId;
	CCString prevClass;
	
	CCString curClass = LOCAL_IRFUT_CLASS;
	CCString stringId = GetEnvVar (l_refXL);
	
	if(!stringId)
	{
		retCode = ARMLOCAL_IRFUT (pDelivery, LocalGetNumObjectId (lIdUnderlying), C_result);

		if(retCode == ARM_OK)
		{
			objId = C_result.getLong ();

			if (l_refXL.GetLen()>0)
						LocalSetCurCellEnvValue (curClass, objId,l_refXL); 

			stringId = LocalMakeObjectId (objId, curClass);
		}
	}
	else
	{
		prevClass = LocalGetStringObjectClass (stringId);
		
		objId = LocalGetNumObjectId (stringId);
			
		if(curClass == prevClass)
		{
			retCode = ARMLOCAL_IRFUT (pDelivery, LocalGetNumObjectId (lIdUnderlying), C_result, objId);
		}
		else
		{
			FreeCurCellContent (stringId);
			retCode = ARMLOCAL_IRFUT (pDelivery, LocalGetNumObjectId (lIdUnderlying), C_result);
		
			if(retCode == ARM_OK)
			{
				objId = C_result.getLong ();
			
				LocalSetCurCellEnvValue (curClass, objId,l_refXL); 

				stringId = LocalMakeObjectId (objId, curClass);
			}
		}
	}	

	_bstr_t tmpChaine = stringId;

	*pRet = SysAllocString(tmpChaine);

	return S_OK;
}
catch(...)
{
	return E_FAIL;
}
}


STDMETHODIMP ActiveXModule::ARMLibor(BSTR pLiborTypeId,BSTR pCcyId, BSTR pResetFreqId, 
								 BSTR pPayFreqId,BSTR pRefObj,BSTR pBasis,BSTR pIntRule,BSTR *pRet)
{
try
{
	ARM_result C_result;

	_bstr_t bLiborTypeId(pLiborTypeId);
	CCString l_liborType = bLiborTypeId;
	long liborTypeId;

	_bstr_t bCcyId(pCcyId);
	CCString l_ccy = bCcyId;

	_bstr_t bResetFreqId(pResetFreqId);
	CCString l_resetFreq = bResetFreqId;
	long resetFreqId;

	_bstr_t bPayFreqId(pPayFreqId);
	CCString l_payFreq = bPayFreqId;
	long payFreqId;

	_bstr_t bBasis(pBasis);
	CCString l_Basis = bBasis;
	long BasisId;

	_bstr_t bIntRule(pIntRule);
	CCString l_intrule = bIntRule;
	long intruleId;

	if((BasisId = ARM_ConvDayCount (l_Basis)) == ARM_DEFAULT_ERR)
	{
		return S_FALSE;
	}

	if((liborTypeId = ARM_ConvIrIndName (l_liborType, C_result)) == ARM_DEFAULT_ERR)
	{
		return S_FALSE;
	}

	if((resetFreqId = ARM_ConvFrequency (l_resetFreq, C_result)) == ARM_DEFAULT_ERR)
	{
		return S_FALSE;
	}

	if((payFreqId = ARM_ConvFrequency (l_payFreq, C_result)) == ARM_DEFAULT_ERR)
	{
		return S_FALSE;
	}
	
	intruleId = ARM_ConvIntRule (l_intrule);

	long retCode;
	long objId;
	CCString prevClass;
	
	_bstr_t bRefObj(pRefObj);
	CCString l_refXL = bRefObj;

	CCString curClass = LOCAL_IRINDEX_CLASS;
	CCString stringId = GetEnvVar (l_refXL);
	
	if(!stringId)
	{
		retCode = ARMLOCAL_LIBOR (liborTypeId, true, l_ccy, resetFreqId, payFreqId,BasisId,intruleId, C_result);

		if(retCode == ARM_OK)
		{
			objId = C_result.getLong ();

			if (l_refXL.GetLen()>0)
						LocalSetCurCellEnvValue (curClass, objId,l_refXL); 


			stringId = LocalMakeObjectId (objId, curClass);
		}
	}
	else
	{
		prevClass = LocalGetStringObjectClass (stringId);
		
		objId = LocalGetNumObjectId (stringId);
			
		if(curClass == prevClass)
		{
			retCode = ARMLOCAL_LIBOR (liborTypeId, true, l_ccy, resetFreqId, payFreqId,BasisId,intruleId, C_result, objId);

		
 		}
		else
		{
			FreeCurCellContent (stringId);
			retCode = ARMLOCAL_LIBOR (liborTypeId, true, l_ccy, resetFreqId, payFreqId,BasisId,intruleId, C_result);

			if(retCode == ARM_OK)
			{
				objId = C_result.getLong ();
			
				LocalSetCurCellEnvValue (curClass, objId,l_refXL); 

				stringId = LocalMakeObjectId (objId, curClass);
			}
		}
	}
	
	_bstr_t tmpChaine = stringId;

	*pRet = SysAllocString(tmpChaine);

	return S_OK;
}
catch(...)
{
	return E_FAIL;
}
}


STDMETHODIMP ActiveXModule::ARMLiborSwaption(double pStartDate,
										 double pEndDate,
										 BSTR pReceiveOrPay,
										 double pStrike, 
										 double pMaturity,
										 BSTR pLiborType,
										 double pSpread,
										 BSTR pExerciseType,
										 BSTR pResetFreq,
										 BSTR pPayFreq,
										 BSTR pCcyId,
										 BSTR pRefObj,
										 BSTR *pRet)
{
try
{
	ARM_result C_result;

	_bstr_t bLiborType(pLiborType);
	CCString l_liborType = bLiborType;
	long liborTypeId;
	
	_bstr_t bReceiveOrPay(pReceiveOrPay);
	CCString l_receiveOrPay = bReceiveOrPay;
	long receiveOrPayId;

	_bstr_t bExerciseType(pExerciseType);
	CCString l_exerciseType = bExerciseType;
	long exerciseTypeId;

	_bstr_t bPayFreq(pPayFreq);
	CCString l_payFreq = bPayFreq;
	long payFreqId ;

	_bstr_t bResetFreq(pResetFreq);
	CCString l_resetFreq = bResetFreq;
	long resetFreqId ;

	_bstr_t bCcyId(pCcyId);
	CCString l_ccy = bCcyId;

	if((liborTypeId = ARM_ConvIrIndName (l_liborType, C_result)) == ARM_DEFAULT_ERR)
		return S_FALSE;

	if((receiveOrPayId = ARM_ConvRecOrPay (l_receiveOrPay, C_result)) == ARM_DEFAULT_ERR)
		return S_FALSE;

	if((exerciseTypeId = ARM_ConvExerciseType (l_exerciseType, C_result)) == ARM_DEFAULT_ERR)
		return S_FALSE;

	if((payFreqId = ARM_ConvFrequency (l_payFreq, C_result)) == ARM_DEFAULT_ERR)
		return S_FALSE;

	if((resetFreqId = ARM_ConvFrequency (l_resetFreq, C_result)) == ARM_DEFAULT_ERR)
		return S_FALSE;

	long retCode;
	long objId;
	CCString prevClass;

	_bstr_t bRefObj(pRefObj);
	CCString l_refXL = bRefObj;

	CCString curClass = LOCAL_SWAPTION_CLASS;
	CCString stringId = GetEnvVar (l_refXL);

	if(!stringId)
	{
		retCode = ARMLOCAL_LIBORSWAPTION (pStartDate, pEndDate, receiveOrPayId,
								     pStrike, pMaturity, liborTypeId,
									 0, pSpread, exerciseTypeId, resetFreqId,
									 payFreqId, true, l_ccy, C_result);

		if(retCode == ARM_OK)
		{
			objId = C_result.getLong ();

			if (l_refXL.GetLen()>0)
						LocalSetCurCellEnvValue (curClass, objId,l_refXL); 



			stringId = LocalMakeObjectId (objId, curClass);
		}
	}
	else
	{
		prevClass = LocalGetStringObjectClass (stringId);
		
		objId = LocalGetNumObjectId (stringId);
			
		if(curClass == prevClass)
		{
			retCode = ARMLOCAL_LIBORSWAPTION (pStartDate, pEndDate, receiveOrPayId,
										 pStrike, pMaturity, liborTypeId,
										 0, pSpread, exerciseTypeId, resetFreqId,
										 payFreqId, true, l_ccy, C_result, objId);
		
		}
		else
		{
			FreeCurCellContent (stringId);
			retCode = ARMLOCAL_LIBORSWAPTION (pStartDate, pEndDate, receiveOrPayId,
								         pStrike, pMaturity, liborTypeId,
									     0, pSpread, exerciseTypeId, resetFreqId,
									     payFreqId, true, l_ccy, C_result);
		
			if(retCode == ARM_OK)
			{
				objId = C_result.getLong ();
			
				LocalSetCurCellEnvValue (curClass, objId,l_refXL); 

				stringId = LocalMakeObjectId (objId, curClass);
			}
		}
	}

	_bstr_t tmpChaine = stringId;

	*pRet = SysAllocString(tmpChaine);

	return S_OK;
}
catch(...)
{
	return E_FAIL;
}
}



STDMETHODIMP ActiveXModule::ARMFixedLeg(double pStartDate, double pEndDate, BSTR pReceiveOrPay, double pFixRate, BSTR pDayCount, BSTR pFreq, BSTR pDecompFreq, BSTR pPayTiming, BSTR pIntRule, BSTR pStubRule, BSTR pCcyId, BSTR pPayCalName, BSTR pNxChange, double pRefDate, BSTR pAdjStartDate,BSTR *pRet)
{
try
{
	ARM_result C_result;

	_bstr_t bReceiveOrPay(pReceiveOrPay);
	CCString l_receiveOrPay = bReceiveOrPay;
	long receiveOrPayId;

	long fixedRateType = 0L;;

	_bstr_t bFreq(pFreq);
	CCString l_freq = bFreq;
	long freqId;
	
	_bstr_t bCcyId(pCcyId);
	CCString l_ccy = bCcyId;

	_bstr_t bDecompFreq(pDecompFreq);
	CCString l_decompfreq = bDecompFreq;
	long decompFreqId;

	_bstr_t bPayTiming(pPayTiming);
	CCString l_payTiming = bPayTiming;
	long payTimingId;

	_bstr_t bIntRule(pIntRule);
	CCString l_intRule = bIntRule;
	long intRuleId;

	_bstr_t bStubRule(pStubRule);
	CCString l_stubRule = bStubRule;
	long stubRuleId;
	
	_bstr_t bPayCalName(pPayCalName);
	CCString l_payCalName = bPayCalName;

	_bstr_t bNxChange(pNxChange);
	CCString l_nxChange = bNxChange;
    long nxChangeId;

	_bstr_t bDayCount(pDayCount);
	CCString l_daycount = bDayCount;
	long dayCountId ;

	_bstr_t bAdjStartDate(pAdjStartDate);
	CCString l_AdjStartDate = bAdjStartDate;
	long adjStartDateId ;

	if((receiveOrPayId = ARM_ConvRecOrPay (l_receiveOrPay, C_result)) == ARM_DEFAULT_ERR)
		return S_FALSE;

	dayCountId = ARM_ConvDayCount (l_daycount);

	if((freqId = ARM_ConvFrequency (l_freq, C_result)) == ARM_DEFAULT_ERR)
		return S_FALSE;

	decompFreqId = ARM_ConvDecompFrequency (l_decompfreq);

	payTimingId = ARM_ConvPayResetRule (l_payTiming);

	intRuleId = ARM_ConvIntRule (l_intRule);

	stubRuleId = ARM_ConvStubRule (l_stubRule);

	nxChangeId = ARM_NotionalExchange(l_nxChange);
	
	if((adjStartDateId = ARM_ConvYesOrNo (l_AdjStartDate, C_result)) == ARM_DEFAULT_ERR)
		return S_FALSE;

	bool ccyIsObject = false;
	if (( l_ccy.GetLen() > 3 )
		&& 
		( !(l_ccy == "DEFAULT"))
	   )
	   ccyIsObject = true;

	CCString curClass = LOCAL_SWAPLEG_CLASS;
	CCString Res;

	if(ARMLOCAL_FIXEDLEG (pStartDate,
								 pEndDate,
								 receiveOrPayId,
								 fixedRateType,
								 pFixRate,
								 dayCountId,
								 freqId,
								 decompFreqId,
								 payTimingId,
								 intRuleId,
								 stubRuleId,
								 ccyIsObject,
								 l_ccy,
								 l_payCalName,
								 nxChangeId,
								 pRefDate,
								 adjStartDateId,
								 -1,
								 C_result) == ARM_OK)
	{
		Res = LocalMakeObjectId (C_result.getLong (), curClass);
	}

	_bstr_t tmpChaine = Res;

	*pRet = SysAllocString(tmpChaine);

	return S_OK;
}
catch(...)
{
	return E_FAIL;
}
}




STDMETHODIMP ActiveXModule::ARMAswPrice(double pMaturity, double pCpn, BSTR pFreq, BSTR pBase, double pMargin, double pRedemptionPrice, double pAsOf, double pDelivery, BSTR pFixDecompfreq, BSTR pCcy1, BSTR pIndex1, BSTR pFwdCurve1, BSTR pDiscCurve1, BSTR pCcy2, BSTR pIndex2, BSTR pFwdCurve2, BSTR pDiscCurve2, BSTR pAmortizationId, long pSolve, double pMinValue, double pMaxValue, double *pRet)
{
	// TODO: Add your implementation code here

try
{
	ARM_result C_result;

	_bstr_t bFreq(pFreq);
	CCString l_freq = bFreq;
	long freqId;

	_bstr_t bBase(pBase);
	CCString l_base = bBase;
	long bondBaseId;

	_bstr_t bFixDecompFreq(pFixDecompfreq);
	CCString l_fixDecompFreq = bFixDecompFreq;
	long fixdecompfreqId;

	_bstr_t bCcy1(pCcy1);
	CCString l_ccy1 = bCcy1;

	_bstr_t bIndex1(pIndex1);
	CCString l_index1 = bIndex1;

	_bstr_t bFwdCurve1(pFwdCurve1);
	CCString l_fwdCurve1 = bFwdCurve1;

	_bstr_t bDiscCurve1(pDiscCurve1);
	CCString l_discCurve1 = bDiscCurve1;
	long discCurve1Id;

	_bstr_t bCcy2(pCcy2);
	CCString l_ccy2 = bCcy2;

	_bstr_t bIndex2(pIndex2);
	CCString l_index2 = bIndex2;

	_bstr_t bFwdCurve2(pFwdCurve2);
	CCString l_fwdCurve2 = bFwdCurve2;

	_bstr_t bDiscCurve2(pDiscCurve2);
	CCString l_discCurve2 = bDiscCurve2;

	_bstr_t bAmort(pAmortizationId);
	CCString l_amort= bAmort;

	if (strcmp((const char*) l_index2,"DEFAULT") == 0)
		l_index2 = l_index1;

	long liborType1Id;
	long liborType2Id;

	if (!(strcmp((const char*)l_ccy2,"NONE") == 0) && (strcmp((const char*)l_ccy1,(const char*)l_ccy2) != 0) )
	{
		if((liborType2Id = ARM_ConvIrIndName (l_index2, C_result)) == ARM_DEFAULT_ERR)
			return S_FALSE;

		discCurve1Id = LocalGetNumObjectId (l_discCurve1);
	}
	else
	{
		if((liborType2Id = ARM_ConvIrIndName (l_index1, C_result)) == ARM_DEFAULT_ERR)
			return S_FALSE;
	}

	long floatResetFreqId;
	long floatPayFreqId;

	if((freqId = ARM_ConvFrequency (l_freq, C_result)) == ARM_DEFAULT_ERR)
		return S_FALSE;

	if((fixdecompfreqId = ARM_ConvFrequency (l_fixDecompFreq, C_result)) == ARM_DEFAULT_ERR)
		return S_FALSE;

	bondBaseId = ARM_ConvDayCount (l_base);

	if((floatResetFreqId = ARM_ConvIrIndNameToFreq (l_index1, C_result)) == ARM_DEFAULT_ERR)
		return S_FALSE;

	if((floatPayFreqId = ARM_ConvIrIndNameToFreq (l_index2, C_result)) == ARM_DEFAULT_ERR)
		return S_FALSE;

	if((liborType1Id = ARM_ConvIrIndName (l_index1, C_result)) == ARM_DEFAULT_ERR)
		return S_FALSE;

	if(ARMLOCAL_ASWPriceNew(pMaturity,0,pCpn,freqId,bondBaseId,pMargin,pRedemptionPrice,pAsOf,pDelivery,
							fixdecompfreqId,floatResetFreqId,floatPayFreqId,l_ccy1,liborType1Id,
							LocalGetNumObjectId (l_fwdCurve1),discCurve1Id,l_ccy2,liborType2Id,
							LocalGetNumObjectId (l_fwdCurve2),LocalGetNumObjectId (l_discCurve2),
							LocalGetNumObjectId (l_amort),pSolve,pMinValue,pMaxValue,C_result) == ARM_OK)
	{
		*pRet = C_result.getDouble ();
	}

	return S_OK;
}
catch(...)
{
	return E_FAIL;
}
}


STDMETHODIMP ActiveXModule::ARMAswMargin(double pMaturity, double pCpn, BSTR pFreq, BSTR pBase, double pPrice, double pRedemptionPrice, double pAsOf, double pDelivery, BSTR pFixDecompfreq, BSTR pCcy1, BSTR pIndex1, BSTR pFwdCurve1, BSTR pDiscCurve1, BSTR pCcy2, BSTR pIndex2, BSTR pFwdCurve2, BSTR pDiscCurve2, BSTR pAmortizationId, long pSolve, double pMinValue, double pMaxValue, double *pRet)
{
	// TODO: Add your implementation code here

try
{
	ARM_result C_result;

	_bstr_t bFreq(pFreq);
	CCString l_freq = bFreq;
	long freqId;

	_bstr_t bBase(pBase);
	CCString l_base = bBase;
	long bondBaseId;

	_bstr_t bFixDecompFreq(pFixDecompfreq);
	CCString l_fixDecompFreq = bFixDecompFreq;
	long fixdecompfreqId;

	_bstr_t bCcy1(pCcy1);
	CCString l_ccy1 = bCcy1;

	_bstr_t bIndex1(pIndex1);
	CCString l_index1 = bIndex1;

	_bstr_t bFwdCurve1(pFwdCurve1);
	CCString l_fwdCurve1 = bFwdCurve1;

	_bstr_t bDiscCurve1(pDiscCurve1);
	CCString l_discCurve1 = bDiscCurve1;
	long discCurve1Id;

	_bstr_t bCcy2(pCcy2);
	CCString l_ccy2 = bCcy2;

	_bstr_t bIndex2(pIndex2);
	CCString l_index2 = bIndex2;

	_bstr_t bFwdCurve2(pFwdCurve2);
	CCString l_fwdCurve2 = bFwdCurve2;

	_bstr_t bDiscCurve2(pDiscCurve2);
	CCString l_discCurve2 = bDiscCurve2;

	_bstr_t bAmort(pAmortizationId);
	CCString l_amort= bAmort;

	if (strcmp((const char*) l_index2,"DEFAULT") == 0)
		l_index2 = l_index1;

	long liborType1Id;
	long liborType2Id;

	if (!(strcmp((const char*)l_ccy2,"NONE") == 0) && (strcmp((const char*)l_ccy1,(const char*)l_ccy2) != 0) )
	{
		if((liborType2Id = ARM_ConvIrIndName (l_index2, C_result)) == ARM_DEFAULT_ERR)
			return S_FALSE;

		discCurve1Id = LocalGetNumObjectId (l_discCurve1);
	}
	else
	{
		if((liborType2Id = ARM_ConvIrIndName (l_index1, C_result)) == ARM_DEFAULT_ERR)
			return S_FALSE;
	}

	long floatResetFreqId;
	long floatPayFreqId;

	if((freqId = ARM_ConvFrequency (l_freq, C_result)) == ARM_DEFAULT_ERR)
		return S_FALSE;

	if((fixdecompfreqId = ARM_ConvFrequency (l_fixDecompFreq, C_result)) == ARM_DEFAULT_ERR)
		return S_FALSE;

	bondBaseId = ARM_ConvDayCount (l_base);

	if((floatResetFreqId = ARM_ConvIrIndNameToFreq (l_index1, C_result)) == ARM_DEFAULT_ERR)
		return S_FALSE;

	if((floatPayFreqId = ARM_ConvIrIndNameToFreq (l_index2, C_result)) == ARM_DEFAULT_ERR)
		return S_FALSE;

	if((liborType1Id = ARM_ConvIrIndName (l_index1, C_result)) == ARM_DEFAULT_ERR)
		return S_FALSE;

	if(ARMLOCAL_ASWMarginNew(pMaturity,0,pCpn,freqId,bondBaseId,pPrice,pRedemptionPrice,pAsOf,pDelivery,
							 fixdecompfreqId,floatResetFreqId,floatPayFreqId,l_ccy1,liborType1Id,
							 LocalGetNumObjectId (l_fwdCurve1),discCurve1Id,l_ccy2,liborType2Id,
							 LocalGetNumObjectId (l_fwdCurve2),LocalGetNumObjectId (l_discCurve2),
							 LocalGetNumObjectId (l_amort),pSolve,pMinValue,pMaxValue,C_result) == ARM_OK)
	{
		*pRet = C_result.getDouble ();
	}

	return S_OK;
}
catch(...)
{
	return E_FAIL;
}
}

STDMETHODIMP ActiveXModule::ARMFrnMargin(double pAsOf, double pDelivery, double pMaturity, BSTR pCcy1, BSTR pIndex1, BSTR pFwdCurve1, BSTR pDiscCurve1, double pFacialMargin, double pPrice, BSTR pCcy2, BSTR pIndex2, BSTR pFwdCurve2, BSTR pDiscCurve2, double pFixing, double pSpread, double pOutMode, long pSolve, BSTR pAmortizationId, double *pRet)
{
	// TODO: Add your implementation code here

try
{
	ARM_result C_result;

	_bstr_t bCcy1(pCcy1);
	CCString l_ccy1 = bCcy1;

	_bstr_t bIndex1(pIndex1);
	CCString l_index1 = bIndex1;

	_bstr_t bFwdCurve1(pFwdCurve1);
	CCString l_fwdCurve1 = bFwdCurve1;

	_bstr_t bDiscCurve1(pDiscCurve1);
	CCString l_discCurve1 = bDiscCurve1;
	long discCurve1Id;

	_bstr_t bCcy2(pCcy2);
	CCString l_ccy2 = bCcy2;

	_bstr_t bIndex2(pIndex2);
	CCString l_index2 = bIndex2;

	_bstr_t bFwdCurve2(pFwdCurve2);
	CCString l_fwdCurve2 = bFwdCurve2;

	_bstr_t bDiscCurve2(pDiscCurve2);
	CCString l_discCurve2 = bDiscCurve2;

	_bstr_t bAmort(pAmortizationId);
	CCString l_amort= bAmort;

	if (strcmp((const char*) l_index2,"DEFAULT") == 0)
		l_index2 = l_index1;

	long liborType1Id;
	long liborType2Id;
	long frequencyId;
	long frequencyId2;

	if (!(strcmp((const char*)l_ccy2,"NONE") == 0) && (strcmp((const char*)l_ccy1,(const char*)l_ccy2) != 0) )
	{
		if((liborType2Id = ARM_ConvIrIndName (l_index2, C_result)) == ARM_DEFAULT_ERR)
			return S_FALSE;

		discCurve1Id = LocalGetNumObjectId (l_discCurve1);
	}
	else
	{
		if((liborType2Id = ARM_ConvIrIndName (l_index1, C_result)) == ARM_DEFAULT_ERR)
			return S_FALSE;
	}

	if((frequencyId2 = ARM_ConvIrIndNameToFreq (l_index2, C_result)) == ARM_DEFAULT_ERR)
		return S_FALSE;
	
	if((liborType1Id = ARM_ConvIrIndName (l_index1, C_result)) == ARM_DEFAULT_ERR)
		return S_FALSE;

	if((frequencyId = ARM_ConvIrIndNameToFreq (l_index1, C_result)) == ARM_DEFAULT_ERR)
		return S_FALSE;

	if(ARMLOCAL_FRNMarginNew (pAsOf,
							  pDelivery,
							  pMaturity,
							  l_ccy1,
							  liborType1Id,
							  LocalGetNumObjectId (l_fwdCurve1),
							  discCurve1Id,
							  pFacialMargin,
							  pPrice,
							  frequencyId,
							  l_ccy2,
							  liborType2Id,
							  LocalGetNumObjectId (l_fwdCurve2),
							  LocalGetNumObjectId (l_discCurve2),
							  LocalGetNumObjectId (l_amort),
							  frequencyId2,
							  pFixing,
							  pSpread,
							  pSolve,
							  C_result) == ARM_OK)
	{
		*pRet = C_result.getDouble ();
	}

	return S_OK;
}
catch(...)
{
	return E_FAIL;
}
}


STDMETHODIMP ActiveXModule::ARMFrnPrice(double pAsOf, double pDelivery, double pMaturity, BSTR pCcy1, BSTR pIndex1, BSTR pFwdCurve1, BSTR pDiscCurve1, double pFacialMargin, double pValoMargin, BSTR pCcy2, BSTR pIndex2, BSTR pFwdCurve2, BSTR pDiscCurve2, double pFixing, double pSpread, double pOutMode, long pSolve, BSTR pAmortizationId, double *pRet)
{
	// TODO: Add your implementation code here

try
{
	ARM_result C_result;

	_bstr_t bCcy1(pCcy1);
	CCString l_ccy1 = bCcy1;

	_bstr_t bIndex1(pIndex1);
	CCString l_index1 = bIndex1;

	_bstr_t bFwdCurve1(pFwdCurve1);
	CCString l_fwdCurve1 = bFwdCurve1;

	_bstr_t bDiscCurve1(pDiscCurve1);
	CCString l_discCurve1 = bDiscCurve1;
	long discCurve1Id;

	_bstr_t bCcy2(pCcy2);
	CCString l_ccy2 = bCcy2;

	_bstr_t bIndex2(pIndex2);
	CCString l_index2 = bIndex2;

	_bstr_t bFwdCurve2(pFwdCurve2);
	CCString l_fwdCurve2 = bFwdCurve2;

	_bstr_t bDiscCurve2(pDiscCurve2);
	CCString l_discCurve2 = bDiscCurve2;

	_bstr_t bAmort(pAmortizationId);
	CCString l_amort= bAmort;

	if (strcmp((const char*) l_index2,"DEFAULT") == 0)
		l_index2 = l_index1;

	long liborType1Id;
	long liborType2Id;
	long frequencyId;
	long frequencyId2;

	if (!(strcmp((const char*)l_ccy2,"NONE") == 0) && (strcmp((const char*)l_ccy1,(const char*)l_ccy2) != 0) )
	{
		if((liborType2Id = ARM_ConvIrIndName (l_index2, C_result)) == ARM_DEFAULT_ERR)
			return S_FALSE;

		discCurve1Id = LocalGetNumObjectId (l_discCurve1);
	}
	else
	{
		if((liborType2Id = ARM_ConvIrIndName (l_index1, C_result)) == ARM_DEFAULT_ERR)
			return S_FALSE;
	}

	if((frequencyId2 = ARM_ConvIrIndNameToFreq (l_index2, C_result)) == ARM_DEFAULT_ERR)
		return S_FALSE;
	
	if((liborType1Id = ARM_ConvIrIndName (l_index1, C_result)) == ARM_DEFAULT_ERR)
		return S_FALSE;

	if((frequencyId = ARM_ConvIrIndNameToFreq (l_index1, C_result)) == ARM_DEFAULT_ERR)
		return S_FALSE;

	if(ARMLOCAL_FRNPriceNew (pAsOf,
							 pDelivery,
							 pMaturity,
							 l_ccy1,
							 liborType1Id,
							 LocalGetNumObjectId (l_fwdCurve1),
							 discCurve1Id,
							 pFacialMargin,
							 pValoMargin,
							 frequencyId,
							 l_ccy2,
							 liborType2Id,
							 LocalGetNumObjectId (l_fwdCurve2),
							 LocalGetNumObjectId (l_discCurve2),
							 LocalGetNumObjectId (l_amort),
							 frequencyId2,
							 pFixing,
							 pSpread,
							 pSolve,
							 C_result) == ARM_OK)
	{
		*pRet = C_result.getDouble ();
	}

	return S_OK;
}
catch(...)
{
	return E_FAIL;
}
}

STDMETHODIMP ActiveXModule::ARMRefValue(VARIANT* pdates,
									VARIANT* pvalues,
									VARIANT* pvalues2,
									long valueType,
									long conversion,
									BSTR calcMethod,
									BSTR *pRefVal)
{
try
{
	// TODO: Add your implementation code here
	// return
//	static XLOPER XL_result;
	ARM_result C_result;

	long datesSize;
	long valueSize;
	long value2Size;
	long value2Vide;
	long calcModId;

	VECTOR<double> vDates;
	VECTOR<double> vValue;
	VECTOR<double> vValue2;

	CCString curclass = LOCAL_REFVAL_CLASS;
	CCString  refvalRes;

	_bstr_t bcalcMethod(calcMethod);
	CCString l_calcMethod = bcalcMethod;

	if(VARIANT2VECTORDOUBLE (*pdates,vDates,datesSize) != S_OK)
		return S_FALSE;

	if(VARIANT2VECTORDOUBLE (*pvalues,vValue,valueSize) != S_OK)
		return S_FALSE;

	if(VARIANT2Long (*pvalues2,&value2Vide) != S_OK) 
	{
		if(VARIANT2VECTORDOUBLE (*pvalues2,vValue2,value2Size) != S_OK)
			return S_FALSE;
	}

	if((calcModId = ARM_ConvCalculationMethod(l_calcMethod, C_result)) == ARM_DEFAULT_ERR)
		return S_FALSE;

	if(ARMLOCAL_REFVALUE(vDates,vValue,vValue2, valueType, conversion, calcModId,C_result) == ARM_OK)
	{
		refvalRes = LocalMakeObjectId (C_result.getLong (),curclass);
	}

	_bstr_t tmpChaine = refvalRes;

	*pRefVal = SysAllocString(tmpChaine);



	return S_OK;
}
catch(...)
{
	return E_FAIL;
}

}

STDMETHODIMP ActiveXModule::ARMCreateGenCorrelManager(VARIANT* pMktTags,
												  VARIANT* pIntraMktTags,
												  VARIANT* pCorrelCurveIds,
												  BSTR *pRet)
{
try
{
	// TODO: Add your implementation code here
	// return
//	static XLOPER XL_result;
	ARM_result C_result;

	long MktTagsSize;
	long IntraMktTagsSize;
	long CorrelCurveIdsSize;

	VECTOR<CCString> vMktTags;
	VECTOR<CCString> vIntraMktTags;
	VECTOR<CCString> vCorrelCurveIds;

	CCString curclass = LOCAL_CORRELMANAGER_CLASS;
	CCString  corrmgrRes;

	if(VARIANT2VECTORCCSTRING (*pMktTags,vMktTags,MktTagsSize) != S_OK)
		return S_FALSE;

	if(VARIANT2VECTORCCSTRING (*pIntraMktTags,vIntraMktTags,IntraMktTagsSize) != S_OK)
		return S_FALSE;

	if(VARIANT2VECTORCCSTRING (*pCorrelCurveIds,vCorrelCurveIds,CorrelCurveIdsSize) != S_OK)
		return S_FALSE;

	vector<long> CorrelCurveIds(CorrelCurveIdsSize);
	for(int i=0; i < CorrelCurveIdsSize; i++)
	{
		CorrelCurveIds[i] = ( strcmp(vCorrelCurveIds[i], "NULL") == 0 )
							? ARM_NULL_OBJECT : LocalGetNumObjectId ( vCorrelCurveIds[i] );
	}

	if(ARMLOCAL_CreateGenCorrelManager(vMktTags,vIntraMktTags,CorrelCurveIds,C_result) == ARM_OK)
	{
		corrmgrRes = LocalMakeObjectId (C_result.getLong (),curclass);
	}

	_bstr_t tmpChaine = corrmgrRes;

	*pRet = SysAllocString(tmpChaine);

	return S_OK;
}
catch(...)
{
	return E_FAIL;
}

}

STDMETHODIMP ActiveXModule::ARMBSConvAdjust(BSTR pSUMMITFormulaeUsed,
											BSTR pUseSABRCMS,
										BSTR *pRet)
{
try
{
	// TODO: Add your implementation code here
	// return
//	static XLOPER XL_result;
	ARM_result C_result;

	CCString curclass = LOCAL_BSCONVADJUST_CLASS;
	CCString  convAdjRes;

	_bstr_t bSUMMITFormulaeUsed(pSUMMITFormulaeUsed);
	CCString l_SUMMITFormulaeUsed = bSUMMITFormulaeUsed;

	_bstr_t bUseSABRCMS(pUseSABRCMS);
	CCString l_UseSABRCMS = bUseSABRCMS;

	long SUMMITFormulaeUsedId;

	if((SUMMITFormulaeUsedId = ARM_ConvSummitFormulae(l_SUMMITFormulaeUsed, C_result)) == ARM_DEFAULT_ERR)
		return S_FALSE;

	long UseSABRCMSId;
	if((UseSABRCMSId = ARM_ConvYesOrNo (l_UseSABRCMS, C_result)) == ARM_DEFAULT_ERR)
		return S_FALSE;

	if(ARMLOCAL_BSConvAdjust_Create(SUMMITFormulaeUsedId,UseSABRCMSId,C_result) == ARM_OK)
	{
		convAdjRes = LocalMakeObjectId (C_result.getLong (),curclass);
	}

	_bstr_t tmpChaine = convAdjRes;

	*pRet = SysAllocString(tmpChaine);

	return S_OK;
}
catch(...)
{
	return E_FAIL;
}

}


STDMETHODIMP ActiveXModule::ARMBsModelGen(BSTR pYieldCurve,
									  BSTR pVolatility,
									  BSTR pCorrMgr,
									  BSTR pCnvxManager,
									  BSTR pCapletVol,
									  BSTR pSpreadLock,
									  BSTR pDiscCurve,
									  BSTR pCorrel,
									  BSTR pCashVol,
									  BSTR pSpreadVol,
									  BSTR pModelType,
									  BSTR pSpreadVolType,
									  BSTR pSabrMod,
									  BSTR pLnOrNorVol,
									  long pNumSteps,
									  BSTR *pRet)
{
try
{
	// TODO: Add your implementation code here
	// return
//	static XLOPER XL_result;
	ARM_result C_result;

	CCString curclass = LOCAL_BSMODELGEN_CLASS;
	CCString bsmodRes;

	_bstr_t bYieldCurve(pYieldCurve);
	CCString l_YieldCurve = bYieldCurve;

	_bstr_t bVolatility(pVolatility);
	CCString l_Volatility = bVolatility;

	_bstr_t bCorrMgr(pCorrMgr);
	CCString l_CorrMgr = bCorrMgr;

	_bstr_t bCnvxManager(pCnvxManager);
	CCString l_CnvxManager = bCnvxManager;

	_bstr_t bCapletVol(pCapletVol);
	CCString l_CapletVol = bCapletVol;

	_bstr_t bSpreadLock(pSpreadLock);
	CCString l_SpreadLock = bSpreadLock;

	_bstr_t bDiscCurve(pDiscCurve);
	CCString l_DiscCurve = bDiscCurve;

	_bstr_t bCorrel(pCorrel);
	CCString l_Correl = bCorrel;

	_bstr_t bCashVol(pCashVol);
	CCString l_CashVol = bCashVol;

	_bstr_t bSpreadVol(pSpreadVol);
	CCString l_SpreadVol = bSpreadVol;

	_bstr_t bModelType(pModelType);
	CCString l_ModelType = bModelType;

	_bstr_t bSpreadVolType(pSpreadVolType);
	CCString l_SpreadVolType = bSpreadVolType;

	_bstr_t bSabrMod(pSabrMod);
	CCString l_SabrMod = bSabrMod;

	_bstr_t blnOrNorVol(pLnOrNorVol);
	CCString l_lnOrNorVol = blnOrNorVol;

	bool isLnVol=true;
	l_lnOrNorVol.toUpper();
    if(CCSTringToSTLString(l_lnOrNorVol)!="Y")
		isLnVol=false;

	long sabrModId;
	long modelTypeId;
	long spreadVolTypeId;

	if(l_SabrMod == "DEFAULT")
	{
		sabrModId = ARM_NULL_OBJECT;
	}
	else
	{
		sabrModId = LocalGetNumObjectId(l_SabrMod);
	}

	if ((modelTypeId = ARM_ConvModelType (l_ModelType, C_result)) == ARM_DEFAULT_ERR)
		return S_FALSE;

	if ((spreadVolTypeId = ARM_ConvVolType2 (l_SpreadVolType, C_result)) == ARM_DEFAULT_ERR)
		return S_FALSE;

	if(ARMLOCAL_BSMODELGEN(LocalGetNumObjectId(l_YieldCurve),
						   (l_SpreadLock == "NULL")? ARM_NULL_OBJECT : LocalGetNumObjectId (l_SpreadLock),
						   (l_CapletVol == "NULL")? ARM_NULL_OBJECT : LocalGetNumObjectId (l_CapletVol),
						   LocalGetNumObjectId(l_Volatility),
						   (l_CorrMgr == "NULL")? ARM_NULL_OBJECT : LocalGetNumObjectId (l_CorrMgr),
						   (l_CnvxManager == "NULL")? ARM_NULL_OBJECT : LocalGetNumObjectId (l_CnvxManager),
						   (l_DiscCurve == "NULL")? ARM_NULL_OBJECT : LocalGetNumObjectId (l_DiscCurve),
						   (l_Correl == "NULL")? ARM_NULL_OBJECT : LocalGetNumObjectId (l_Correl),
						   (l_CashVol == "NULL")? ARM_NULL_OBJECT : LocalGetNumObjectId (l_CashVol),
						   (l_SpreadVol == "NULL")? ARM_NULL_OBJECT : LocalGetNumObjectId (l_SpreadVol),
						   modelTypeId,
						   spreadVolTypeId,
						   sabrModId,
						   isLnVol,
						   pNumSteps,
						   C_result) == ARM_OK)
	{
		bsmodRes = LocalMakeObjectId (C_result.getLong (),curclass);
	}

	_bstr_t tmpChaine = bsmodRes;

	*pRet = SysAllocString(tmpChaine);

	return S_OK;
}
catch(...)
{
	return E_FAIL;
}

}


STDMETHODIMP ActiveXModule::ARMDisplaySchedule(BSTR pLegId,
										   BSTR pDataType,
										   VARIANT *pRet)
{
try
{
	// TODO: Add your implementation code here
	// return
//	static XLOPER XL_result;
	ARM_result C_result;

	VECTOR<double> dResult;

	_bstr_t bLegId(pLegId);
	CCString l_LegId = bLegId;

	_bstr_t bDataType(pDataType);
	CCString l_DataType = bDataType;

	l_DataType.toUpper();

// FIXMEFRED: mig.vc8 (31/05/2007 11:08:08):use string methods
	if (strstr(l_DataType.c_str(),"DATE") != NULL)
	{
		long typeDatesId;

		if ((typeDatesId = ARM_ConvTypeDates (l_DataType, C_result)) == ARM_DEFAULT_ERR)
			return S_FALSE;

		if(ARMLOCAL_ARM_DisplayScheduleDates(LocalGetNumObjectId(l_LegId),
											 typeDatesId,
											 K_RCV,
											 K_NO,
											 C_result) == ARM_OK)
		{
			VECTOR<CCString> dDate;
			long vecSize;
			long retCode = ExtractVectorDateFromFile("123",dDate,vecSize);

			dDate.resize(vecSize);

			if ( (retCode == ARM_OK) )
			{
				int i;
				for(i=0;i<dDate.size(); ++i)
				dResult.push_back(Local_ARMDATE2XLDATE( dDate[i] ));
			}
		}
	}
	else
	{
		long typeValuesId;

		if ((typeValuesId = ARM_ConvTypeValues (l_DataType, C_result)) == ARM_DEFAULT_ERR)
			return S_FALSE;

		if(ARMLOCAL_ARM_DisplayScheduleValues(LocalGetNumObjectId(l_LegId),
											  typeValuesId,
											  K_RCV,
											  K_NO,
											  C_result) == ARM_OK)
		{
			VECTOR<double> dVal;
			long vecSize;
			long retCode = ExtractVectorDoubleFromFile("123",dVal,vecSize);

			dVal.resize(vecSize);

			if ( (retCode == ARM_OK) )
			{
				int i;
				for(i=0;i<dVal.size(); ++i)
				dResult.push_back(dVal[i]);
			}
		}
	}

	long Res;

	Res = VECTORDOUBLE2VARIANT(dResult, pRet);

	return S_OK;
}
catch(...)
{
	return E_FAIL;
}

}


STDMETHODIMP ActiveXModule::ARMIrIndex(BSTR pDaycount,BSTR pPayFreq,double pMaturity,BSTR pCompMethod,BSTR pFwdRule,BSTR pResetTiming,double pResetGap,BSTR pPayTiming,double pPayGap,BSTR pCcy,BSTR pIndexType,double pDecompFreq,BSTR pIntRule,BSTR pResetFreq, BSTR *pRet)
{
try
{
	// TODO: Add your implementation code here
	// return
	ARM_result C_result;

	_bstr_t bDaycount(pDaycount);
	CCString l_Daycount = bDaycount;

	_bstr_t bPayFreq(pPayFreq);
	CCString l_PayFreq = bPayFreq;

	_bstr_t bCompMethod(pCompMethod);
	CCString l_CompMethod = bCompMethod;

	_bstr_t bFwdRule(pFwdRule);
	CCString l_FwdRule = bFwdRule;

	_bstr_t bResetTiming(pResetTiming);
	CCString l_ResetTiming = bResetTiming;

	_bstr_t bPayTiming(pPayTiming);
	CCString l_PayTiming = bPayTiming;

	_bstr_t bCcy(pCcy);
	CCString l_Ccy = bCcy;

	_bstr_t bIndexType(pIndexType);
	CCString l_IndexType = bIndexType;

	_bstr_t bIntRule(pIntRule);
	CCString l_IntRule = bIntRule;

	_bstr_t bResetFreq(pResetFreq);
	CCString l_ResetFreq = bResetFreq;

	bool ccyIsObject = false;
	if (( l_Ccy.GetLen() > 3 )
		&& 
		( !(l_Ccy == "DEFAULT"))
	   )
	   ccyIsObject = true;

	long indexTypeId = ARM_ConvIrType (l_IndexType);

	long dayCountId;
	if (strcmp(l_Daycount,"-1") == 0)
	{
		dayCountId = -1;
	}
	else
	{
		dayCountId = ARM_ConvDayCount (l_Daycount);
	}

	if( strcmp(l_ResetFreq , "-1") == 0 )
	{
		l_ResetFreq = l_PayFreq;
	}

	long payFreqId;
	if((payFreqId = ARM_ConvFrequency (l_PayFreq, C_result)) == ARM_DEFAULT_ERR)
		return S_FALSE;

	long resetFreqId;
	if((resetFreqId = ARM_ConvFrequency (l_ResetFreq, C_result)) == ARM_DEFAULT_ERR)
		return S_FALSE;

	long compMethodId;
	if((compMethodId = ARM_ConvCompMeth (l_CompMethod, C_result)) == ARM_DEFAULT_ERR)
		return S_FALSE;

	long fwdRuleId;
	if((fwdRuleId = ARM_ConvFwdRule (l_FwdRule, C_result)) == ARM_DEFAULT_ERR)
		return S_FALSE;

	long resetTimingId = ARM_ConvPayResetRule (l_ResetTiming);

	long payTimingId = ARM_ConvPayResetRule (l_PayTiming);

	long intRuleId = ARM_ConvIntRule (l_IntRule);

	CCString Res;
	CCString curClass = LOCAL_IRINDEX_CLASS;

	if(ARMLOCAL_IRINDEX (dayCountId, payFreqId, pMaturity, compMethodId,
						 fwdRuleId, resetTimingId, pResetGap, payTimingId,
						 pPayGap, ccyIsObject, l_Ccy, indexTypeId, (long)pDecompFreq,
						 intRuleId, resetFreqId, C_result) == ARM_OK)
	{
		Res = LocalMakeObjectId (C_result.getLong (), curClass);
	}

	_bstr_t tmpChaine = Res;

	*pRet = SysAllocString(tmpChaine);


	return S_OK;
}
catch(...)
{
	return E_FAIL;
}

}



STDMETHODIMP ActiveXModule::ARMSwapleg(BSTR pIndexId, double pStartDate, double pEndDate, BSTR pRecOrPay, VARIANT pSpread, BSTR pCcy, BSTR pDayCount, double pResetGap, BSTR pResetCal, BSTR pPayCal, double pDecompPricingFlag, BSTR pNxChange, BSTR pStubRule, double pRefDate, BSTR pAdjStartDate, BSTR *pRet)
{
try
{
	// TODO: Add your implementation code here
	// return
	ARM_result C_result;

	_bstr_t bIndexId(pIndexId);
	CCString l_IndexId = bIndexId;

	_bstr_t bRecOrPay(pRecOrPay);
	CCString l_RecOrPay = bRecOrPay;

	_bstr_t bCcy(pCcy);
	CCString l_Ccy = bCcy;

	_bstr_t bDayCount(pDayCount);
	CCString l_DayCount = bDayCount;

	_bstr_t bResetCal(pResetCal);
	CCString l_ResetCal = bResetCal;

	_bstr_t bPayCal(pPayCal);
	CCString l_PayCal = bPayCal;

	_bstr_t bNxChange(pNxChange);
	CCString l_NxChange = bNxChange;

	_bstr_t bStubRule(pStubRule);
	CCString l_StubRule = bStubRule;

	_bstr_t bAdjStartDate(pAdjStartDate);
	CCString l_AdjStartDate = bAdjStartDate;

	bool ccyIsObject = false;
	if (( l_Ccy.GetLen() > 3 )
		&& 
		( !(l_Ccy == "DEFAULT"))
	   )
	   ccyIsObject = true;

	long receiveOrPayId;
	if((receiveOrPayId = ARM_ConvRecOrPay (l_RecOrPay, C_result)) == ARM_DEFAULT_ERR)
		return S_FALSE;

	double dSpread;
	long   spreadType;

	if (VARIANT2Double(pSpread,&dSpread) == S_FALSE)
	{
		_bstr_t bSpread(pSpread);
		CCString l_Spread = bSpread;

		dSpread = (double) LocalGetNumObjectId(l_Spread);
		spreadType = 1L;
	}
	else
	{
	   spreadType = 0L;
	}

	long dayCountId = ARM_ConvDayCount (l_DayCount);

	long nxChange = ARM_NotionalExchange(l_NxChange);

	long stubRuleId = ARM_ConvStubRule (l_StubRule);

	long adjStartDateId;
	if((adjStartDateId = ARM_ConvYesOrNo (l_AdjStartDate, C_result)) == ARM_DEFAULT_ERR)
		return S_FALSE;

	CCString Res;
	CCString curClass = LOCAL_SWAPLEG_CLASS;

	if(ARMLOCAL_SWAPLEG (LocalGetNumObjectId (l_IndexId),
						 pStartDate,
						 pEndDate,
						 receiveOrPayId,
						 spreadType,
						 dSpread,
						 ccyIsObject,
						 l_Ccy,
						 dayCountId,
						 (long) pResetGap,
						 l_ResetCal,
						 l_PayCal,
						 (long) pDecompPricingFlag,
						 nxChange,
						 stubRuleId,
						 pRefDate,
						 adjStartDateId,
						 -1,
						 C_result) == ARM_OK)
	{
		Res = LocalMakeObjectId (C_result.getLong (), curClass);
	}

	_bstr_t tmpChaine = Res;

	*pRet = SysAllocString(tmpChaine);

	return S_OK;
}
catch(...)
{
	return E_FAIL;
}

}


STDMETHODIMP ActiveXModule::ARMConstRefvalue(double pValue, BSTR *pRet)
{
try
{
	// TODO: Add your implementation code here
	// return
	ARM_result C_result;

	CCString Res;
	CCString curClass = LOCAL_REFVAL_CLASS;

	if(ARMLOCAL_CONSTREFVALUE (pValue, C_result) == ARM_OK)
	{
		Res = LocalMakeObjectId (C_result.getLong (), curClass);
	}

	_bstr_t tmpChaine = Res;

	*pRet = SysAllocString(tmpChaine);

	return S_OK;
}
catch(...)
{
	return E_FAIL;
}

}



STDMETHODIMP ActiveXModule::ARMBond(double pIssueDate,double pMaturityDate,double pFirstCpnDate,double pCpnRate,double pRedempPrice,double pPeriodicity,VARIANT pDaycount,double pSettleGap,double pCpnDateFlag,BSTR pCcy,BSTR *pRet)
{
try
{
	// TODO: Add your implementation code here
	// return
	ARM_result C_result;

	double dDaycount;

	if (VARIANT2Double(pDaycount,&dDaycount) == S_FALSE)
	{
		_bstr_t bDaycount(pDaycount);
		CCString l_Daycount = bDaycount;

		dDaycount = (double)ARM_ConvDayCount (l_Daycount);
	}

	_bstr_t bCcy(pCcy);
	CCString l_Ccy = bCcy;

	long l_ccyId;

	if (l_Ccy == "DEFAULT")
		l_ccyId = -1;
	else
		l_ccyId = LocalGetNumObjectId (l_Ccy);

	CCString Res;
	CCString curClass = LOCAL_BOND_CLASS;

	if(ARMLOCAL_bond (pIssueDate, pMaturityDate, pFirstCpnDate, pCpnRate,
							pRedempPrice, (long)pPeriodicity,
							(long)dDaycount, (long)pSettleGap,
							(long)pCpnDateFlag, l_ccyId , C_result) == ARM_OK)
	{
		Res = LocalMakeObjectId (C_result.getLong (), curClass);
	}

	_bstr_t tmpChaine = Res;

	*pRet = SysAllocString(tmpChaine);

	return S_OK;
}
catch(...)
{
	return E_FAIL;
}

}


STDMETHODIMP ActiveXModule::ARMINFCreateOATLeg(double pStartDate,double pEndDate,BSTR pInfIdx,BSTR pRcvOrPay,BSTR pInterpType,double pLeverage,double pSpread,BSTR pResetFreq,BSTR pDaycount,BSTR pResetCal,BSTR pFwdRule,BSTR pIntRule,BSTR pStubRule,double pResetNumGap,double pResetDenomGap,BSTR pPayFreq,double pPayGap,BSTR pPayCal,BSTR pFinalNotionalType,double pFirstReset,double pCoMultiple,BSTR *pRet)
{
try
{
	// TODO: Add your implementation code here
	// return
	ARM_result C_result;

	_bstr_t bInfIdx(pInfIdx);
	CCString l_InfIdx = bInfIdx;

	_bstr_t bRcvOrPay(pRcvOrPay);
	CCString l_RcvOrPay = bRcvOrPay;

	long receiveOrPayId;
	if((receiveOrPayId = ARM_ConvRecOrPay (l_RcvOrPay, C_result)) == ARM_DEFAULT_ERR)
		return S_FALSE;

	_bstr_t bInterpType(pInterpType);
	CCString l_InterpType = bInterpType;

	long interpTypeId;
	if((interpTypeId = ARM_ConvCPIDailyInterpMethod (l_InterpType, C_result)) == ARM_DEFAULT_ERR)
		return S_FALSE;

	_bstr_t bResetFreq(pResetFreq);
	CCString l_ResetFreq = bResetFreq;

	long resetFreqId;
	if((resetFreqId = ARM_ConvFrequency (l_ResetFreq, C_result)) == ARM_DEFAULT_ERR)
		return S_FALSE;
	
	_bstr_t bDaycount(pDaycount);
	CCString l_Daycount = bDaycount;

	long daycountId;
	daycountId = ARM_ConvDayCount (l_Daycount);

	_bstr_t bResetCal(pResetCal);
	CCString l_ResetCal = bResetCal;

	_bstr_t bFwdRule(pFwdRule);
	CCString l_FwdRule = bFwdRule;

	long fwdRuleId;
	if((fwdRuleId = ARM_ConvFwdRule (l_FwdRule, C_result)) == ARM_DEFAULT_ERR)
		return S_FALSE;

	_bstr_t bIntRule(pIntRule);
	CCString l_IntRule = bIntRule;

	long intRuleId;
	intRuleId = ARM_ConvIntRule (l_IntRule);

	_bstr_t bStubRule(pStubRule);
	CCString l_StubRule = bStubRule;

	long stubRuleId;
	stubRuleId = ARM_ConvStubRule (l_StubRule);

	if (pResetDenomGap == -1111)
		pResetDenomGap = pResetNumGap;

	_bstr_t bPayFreq(pPayFreq);
	CCString l_PayFreq = bPayFreq;

	long payFreqId;
	if((payFreqId = ARM_ConvFrequency (l_PayFreq, C_result)) == ARM_DEFAULT_ERR)
		return S_FALSE;

	_bstr_t bPayCal(pPayCal);
	CCString l_PayCal = bPayCal;

	_bstr_t bFinalNotionalType(pFinalNotionalType);
	CCString l_FinalNotionalType = bFinalNotionalType;

	long finalNotionalTypeId;
	if((finalNotionalTypeId = ARM_NotionalType (l_FinalNotionalType, C_result)) == ARM_DEFAULT_ERR)
		return S_FALSE;

	CCString Res;
	CCString curClass = LOCAL_INFSWAPLEG_CLASS;

	if(ARMLOCAL_InfLegAllInputs_Create (pStartDate, pEndDate, l_InfIdx, K_OATTYPE_LEG,
										receiveOrPayId,interpTypeId,pLeverage,pCoMultiple,pSpread,
										resetFreqId,daycountId,l_ResetCal,fwdRuleId,
										intRuleId,stubRuleId,pResetNumGap,pResetDenomGap,
										payFreqId,pPayGap,l_PayCal,K_ADJUSTED,
										finalNotionalTypeId,pFirstReset,C_result) == ARM_OK)
	{
		Res = LocalMakeObjectId (C_result.getLong (), curClass);
	}

	_bstr_t tmpChaine = Res;

	*pRet = SysAllocString(tmpChaine);

	return S_OK;
}
catch(...)
{
	return E_FAIL;
}

}


STDMETHODIMP ActiveXModule::ARMSwap(BSTR pSwapleg1,BSTR pSwapleg2,double pMinPay,BSTR *pRet)
{
try
{
	// TODO: Add your implementation code here
	// return
	ARM_result C_result;

	_bstr_t bSwapleg1(pSwapleg1);
	CCString l_Swapleg1 = bSwapleg1;

	_bstr_t bSwapleg2(pSwapleg2);
	CCString l_Swapleg2 = bSwapleg2;

	VECTOR<double> fixedRates;

	CCString Res;
	CCString curClass = LOCAL_SWAP_CLASS;

	if(ARMLOCAL_SWAP (LocalGetNumObjectId(l_Swapleg1),
					  LocalGetNumObjectId(l_Swapleg2),
					  pMinPay,
					  fixedRates,
					  C_result) == ARM_OK)
	{
		Res = LocalMakeObjectId (C_result.getLong (), curClass);
	}

	_bstr_t tmpChaine = Res;

	*pRet = SysAllocString(tmpChaine);

	return S_OK;
}
catch(...)
{
	return E_FAIL;
}

}


STDMETHODIMP ActiveXModule::ARMPToYield(BSTR pBond,double pSettleDate,double pPrice,VARIANT *pRet)
{
	// TODO: Add your implementation code here

try
{
	ARM_result C_result;

	_bstr_t bBond(pBond);
	CCString l_Bond = bBond;

	double Res;

	if(ARMLOCAL_PTOYIELD (LocalGetNumObjectId (l_Bond), pSettleDate, pPrice, C_result) == ARM_OK)
	{
		Res = C_result.getDouble ();
	}

	_variant_t wrap_Res;
	wrap_Res.Attach(*pRet);
	wrap_Res.dblVal = Res;
	*pRet=wrap_Res.Detach();

	pRet->vt=VT_R8;

	return S_OK;
}
catch(...)
{
	return E_FAIL;
}
}


STDMETHODIMP ActiveXModule::ARMYToPrice(BSTR pBond,double pSettleDate,double pYield,VARIANT *pRet)
{
	// TODO: Add your implementation code here

try
{
	ARM_result C_result;

	_bstr_t bBond(pBond);
	CCString l_Bond = bBond;

	double Res;

	if(ARMLOCAL_PTOYIELD (LocalGetNumObjectId (l_Bond), pSettleDate, pYield, C_result) == ARM_OK)
	{
		Res = C_result.getDouble ();
	}

	_variant_t wrap_Res;
	wrap_Res.Attach(*pRet);
	wrap_Res.dblVal = Res;
	*pRet=wrap_Res.Detach();

	pRet->vt=VT_R8;

	return S_OK;
}
catch(...)
{
	return E_FAIL;
}
}


STDMETHODIMP ActiveXModule::ARMYToDuration(BSTR pBond,double pSettleDate,double pActuRate,double pFlagCpn,VARIANT *pRet)
{
	// TODO: Add your implementation code here

try
{
	ARM_result C_result;

	_bstr_t bBond(pBond);
	CCString l_Bond = bBond;

	double Res;

	if(ARMLOCAL_YTODURATION (LocalGetNumObjectId (l_Bond), pSettleDate, pActuRate, pFlagCpn, C_result) == ARM_OK)
	{
		Res = C_result.getDouble ();
	}

	_variant_t wrap_Res;
	wrap_Res.Attach(*pRet);
	wrap_Res.dblVal = Res;
	*pRet=wrap_Res.Detach();

	pRet->vt=VT_R8;

	return S_OK;
}
catch(...)
{
	return E_FAIL;
}
}


STDMETHODIMP ActiveXModule::ARMLiborleg(double pStartDate,double pEndDate,BSTR pLiborType,BSTR pRecOrPay,VARIANT pSpread,BSTR pResetFreq,BSTR pPayFreq,BSTR pResetTiming,BSTR pPayTiming,BSTR pCcy,BSTR pIntRule,double pResetGap,BSTR pResetCal,BSTR pPayCal,double pDecompPricingFlag,BSTR pNxChange,BSTR pStubRule,double pRefDate,BSTR pAdjStartDate,BSTR pCpnDaycount,BSTR *pRet){
	// TODO: Add your implementation code here

try
{
	ARM_result C_result;

	_bstr_t bLiborType(pLiborType);
	CCString l_LiborType = bLiborType;
	long liborTypeId; 

	if((liborTypeId = ARM_ConvIrIndName (l_LiborType, C_result)) == ARM_DEFAULT_ERR)
		return S_FALSE;

	_bstr_t bRecOrPay(pRecOrPay);
	CCString l_RecOrPay = bRecOrPay;
	long recOrPayId; 

	if((recOrPayId = ARM_ConvRecOrPay (l_RecOrPay, C_result)) == ARM_DEFAULT_ERR)
		return S_FALSE;

	double dSpread;
	long   spreadType;

	if (VARIANT2Double(pSpread,&dSpread) == S_FALSE)
	{
		_bstr_t bSpread(pSpread);
		CCString l_Spread = bSpread;

		dSpread = (double) LocalGetNumObjectId(l_Spread);
		spreadType = 1L;
	}
	else
	{
	   spreadType = 0L;
	}

	_bstr_t bResetFreq(pResetFreq);
	CCString l_ResetFreq = bResetFreq;

	long resetFreqId;
	if((resetFreqId = ARM_ConvFrequency (l_ResetFreq, C_result)) == ARM_DEFAULT_ERR)
		return S_FALSE;

	_bstr_t bPayFreq(pPayFreq);
	CCString l_PayFreq = bPayFreq;

	long payFreqId;
	if((payFreqId = ARM_ConvFrequency (l_PayFreq, C_result)) == ARM_DEFAULT_ERR)
		return S_FALSE;

	_bstr_t bResetTiming(pResetTiming);
	CCString l_ResetTiming = bResetTiming;
	long resetTimingId = ARM_ConvPayResetRule (l_ResetTiming);

	_bstr_t bPayTiming(pPayTiming);
	CCString l_PayTiming = bPayTiming;
	long payTimingId = ARM_ConvPayResetRule (l_PayTiming);

	_bstr_t bCcy(pCcy);
	CCString l_Ccy = bCcy;

	bool ccyIsObject = false;
	if (( l_Ccy.GetLen() > 3 )
		&& 
		( !(l_Ccy == "DEFAULT"))
	   )
	   ccyIsObject = true;

	_bstr_t bIntRule(pIntRule);
	CCString l_IntRule = bIntRule;

	long intRuleId;
	intRuleId = ARM_ConvIntRule (l_IntRule);

	_bstr_t bResetCal(pResetCal);
	CCString l_ResetCal = bResetCal;

	_bstr_t bPayCal(pPayCal);
	CCString l_PayCal = bPayCal;

	_bstr_t bNxChange(pNxChange);
	CCString l_NxChange = bNxChange;
	long nxChange = ARM_NotionalExchange(l_NxChange);

	_bstr_t bStubRule(pStubRule);
	CCString l_StubRule = bStubRule;
	long stubRuleId = ARM_ConvStubRule (l_StubRule);

	_bstr_t bAdjStartDate(pAdjStartDate);
	CCString l_AdjStartDate = bAdjStartDate;
	long adjStartDateId;
	if((adjStartDateId = ARM_ConvYesOrNo (l_AdjStartDate, C_result)) == ARM_DEFAULT_ERR)
		return S_FALSE;

	_bstr_t bCpnDaycount(pCpnDaycount);
	CCString l_CpnDaycount = bCpnDaycount;
	long dayCountId = ARM_ConvDayCount(l_CpnDaycount);

	CCString Res;
	CCString curClass = LOCAL_SWAPLEG_CLASS;

	if(ARMLOCAL_LIBORLEG (pStartDate,
						  pEndDate,
						  liborTypeId,
						  recOrPayId,
						  spreadType,
						  dSpread,
						  resetFreqId,
						  payFreqId,
						  resetTimingId,
						  payTimingId,
						  ccyIsObject,
						  l_Ccy,
						  intRuleId,
						  (long) pResetGap,
						  l_ResetCal,
						  l_PayCal,
						  (long) pDecompPricingFlag,
						  nxChange,
						  stubRuleId,
						  pRefDate,
						  adjStartDateId,
						  dayCountId,
						  C_result) == ARM_OK)
	{
		Res = LocalMakeObjectId (C_result.getLong (), curClass);
	}

	_bstr_t tmpChaine = Res;

	*pRet = SysAllocString(tmpChaine);

	return S_OK;
}
catch(...)
{
	return E_FAIL;
}
}


STDMETHODIMP ActiveXModule::ARMImpliedSpread(BSTR pSwap,BSTR pModel,double pPrice,double pLeg1or2,double *pRet)
{
	// TODO: Add your implementation code here

try
{
	ARM_result C_result;

	_bstr_t bSwap(pSwap);
	CCString l_Swap = bSwap;

	_bstr_t bModel(pModel);
	CCString l_Model = bModel;

	if(ARMLOCAL_IMPLIEDSPREAD (LocalGetNumObjectId (l_Swap), LocalGetNumObjectId (l_Model), pPrice, (long)pLeg1or2, C_result) == ARM_OK)
	{
		*pRet = C_result.getDouble ();
	}

	return S_OK;
}
catch(...)
{
	return E_FAIL;
}
}


STDMETHODIMP ActiveXModule::ARMDiscountPrice(BSTR pZeroCurve,double pMatu,double *pRet)
{
	// TODO: Add your implementation code here

try
{
	ARM_result C_result;

	_bstr_t bZeroCurve(pZeroCurve);
	CCString l_ZeroCurve = bZeroCurve;

	if(ARMLOCAL_DiscountPrice (LocalGetNumObjectId (l_ZeroCurve), pMatu, C_result) == ARM_OK)
	{
		*pRet = C_result.getDouble ();
	}

	return S_OK;
}
catch(...)
{
	return E_FAIL;
}
}


STDMETHODIMP ActiveXModule::ARMINFCreateCurve(double pAsOf,BSTR pIndexName,double pCPIIndexValue,double pCPIIndexDate,VARIANT* pMatu,VARIANT* pRate,BSTR pMonthlyInterpType,BSTR pDailyInterpType,BSTR pDCFMonthly,BSTR pDCFDaily,BSTR pExtrapolType,BSTR pResetManager,BSTR pSeasonManager,BSTR *pRet)
{
	// TODO: Add your implementation code here

try
{
	ARM_result C_result;

	_bstr_t bIndexName(pIndexName);
	CCString l_IndexName = bIndexName;

	VECTOR<CCString> vMatu;
	long size;
	if(VARIANT2VECTORCCSTRING (*pMatu,vMatu,size) != S_OK)
		return S_FALSE;

	VECTOR<double> vRate;
	if(VARIANT2VECTORDOUBLE (*pRate,vRate,size) != S_OK)
		return S_FALSE;

	_bstr_t bMonthlyInterpType(pMonthlyInterpType);
	CCString l_MonthlyInterpType = bMonthlyInterpType;

	long monthlyInterpTypeId;
	if((monthlyInterpTypeId = ARM_ConvCPIMonthlyInterpMethod (l_MonthlyInterpType, C_result)) == ARM_DEFAULT_ERR)
		return S_FALSE;

	_bstr_t bDailyInterpType(pDailyInterpType);
	CCString l_DailyInterpType = bDailyInterpType;

	long dailyInterpTypeId;
	if((dailyInterpTypeId = ARM_ConvCPIDailyInterpMethod (l_DailyInterpType, C_result)) == ARM_DEFAULT_ERR)
		return S_FALSE;

	_bstr_t bExtrapolType(pExtrapolType);
	CCString l_ExtrapolType = bExtrapolType;

	long extrapolTypeId;
	if((extrapolTypeId = ARM_ConvCPIExtrapolMethod (l_ExtrapolType, C_result)) == ARM_DEFAULT_ERR)
		return S_FALSE;

	_bstr_t bDCFMonthly(pDCFMonthly);
	CCString l_DCFMonthly = bDCFMonthly;
	long DCFMonthlyId = ARM_ConvDayCount( l_DCFMonthly);
		
	_bstr_t bDCFDaily(pDCFDaily);
	CCString l_DCFDaily = bDCFDaily;
	long DCFDailyId = ARM_ConvDayCount( l_DCFDaily);

	_bstr_t bResetManager(pResetManager);
	CCString l_ResetManager = bResetManager;
	long resetMgrId;
	if (l_ResetManager == "DEFAULT")
		resetMgrId = ARM_NULL_OBJECT;
	else
		resetMgrId = LocalGetNumObjectId(l_ResetManager);

	_bstr_t bSeasonManager(pSeasonManager);
	CCString l_SeasonManager = bSeasonManager;
	long seasonMgrId;
	if (l_SeasonManager == "DEFAULT")
		seasonMgrId = ARM_NULL_OBJECT;
	else
		seasonMgrId = LocalGetNumObjectId(l_SeasonManager);

	CCString Res;
	CCString curClass = LOCAL_INFCURV_CLASS;
	
	if(ARMLOCAL_InfCurv_Create (pAsOf,
								l_IndexName,
								pCPIIndexValue,
								pCPIIndexDate,
								vMatu,
								vRate,
								monthlyInterpTypeId,
								dailyInterpTypeId,
								DCFMonthlyId,
								DCFDailyId,
								extrapolTypeId,
								resetMgrId,
								seasonMgrId,
								C_result) == ARM_OK)
	{
		Res = LocalMakeObjectId (C_result.getLong (), curClass);
	}

	_bstr_t tmpChaine = Res;

	*pRet = SysAllocString(tmpChaine);

	return S_OK;
}
catch(...)
{
	return E_FAIL;
}
}


STDMETHODIMP ActiveXModule::ARMINFInterpCPI(BSTR pZc,double pCPIDate,BSTR pDCFlag,BSTR pDailyInterpType,BSTR pCPIlag,double pWeight,double *pRet)
{
	// TODO: Add your implementation code here

try
{
	ARM_result C_result;

	_bstr_t bZc(pZc);
	CCString l_Zc = bZc;

	_bstr_t bDCFlag(pDCFlag);
	CCString l_DCFlag = bDCFlag;

	_bstr_t bDailyInterpType(pDailyInterpType);
	CCString l_DailyInterpType = bDailyInterpType;

	long dailyInterpTypeId;
	if((dailyInterpTypeId = ARM_ConvCPIDailyInterpMethod (l_DailyInterpType, C_result)) == ARM_DEFAULT_ERR)
		return S_FALSE;

	_bstr_t bCPIlag(pCPIlag);
	CCString l_CPIlag = bCPIlag;

	if(ARMLOCAL_InfCurv_Interp (LocalGetNumObjectId (l_Zc), pCPIDate, C_result, l_DCFlag, dailyInterpTypeId, l_CPIlag, pWeight, CPI) == ARM_OK)
	{
		*pRet = C_result.getDouble ();
	}

	return S_OK;
}
catch(...)
{
	return E_FAIL;
}
}


STDMETHODIMP ActiveXModule::ARMINFSeasonManager(VARIANT* pMonthList,VARIANT* pValues,BSTR pSeasonAdjMode,BSTR *pRet)
{
	// TODO: Add your implementation code here

try
{
	ARM_result C_result;

	_bstr_t bSeasonAdjMode(pSeasonAdjMode);
	CCString l_SeasonAdjMode = bSeasonAdjMode;

	long SeasonAdjModeId;
	if((SeasonAdjModeId = ARM_ConvSeasonalityMode (l_SeasonAdjMode, C_result)) == ARM_DEFAULT_ERR)
		return S_FALSE;

	VECTOR<CCString> vMonthList;
	long size;
	if(VARIANT2VECTORCCSTRING (*pMonthList,vMonthList,size) != S_OK)
		return S_FALSE;

	VECTOR<double> vValues;
	if(VARIANT2VECTORDOUBLE (*pValues,vValues,size) != S_OK)
		return S_FALSE;

	VECTOR<double>	 C_seasonHorizonList;

	CCString Res;
	CCString curClass = LOCAL_SEASONMANAGER_CLASS;
	
	if(ARMLOCAL_SeasonalityManager_Create (vMonthList,
										   vValues,
										   C_seasonHorizonList,
										   SeasonAdjModeId,
										   C_result) == ARM_OK)
	{
		Res = LocalMakeObjectId (C_result.getLong (), curClass);
	}

	_bstr_t tmpChaine = Res;

	*pRet = SysAllocString(tmpChaine);

	return S_OK;
}
catch(...)
{
	return E_FAIL;
}
}


STDMETHODIMP ActiveXModule::ARMINFResetManager(VARIANT* pDatas,double pNbIndex,BSTR *pRet)
{
	// TODO: Add your implementation code here

try
{
	ARM_result C_result;

	VECTOR<CCString> vDatas;
	long size;
	if(VARIANT2VECTORCCSTRING (*pDatas,vDatas,size) != S_OK)
		return S_FALSE;

	CCString Res;
	CCString curClass = LOCAL_RESETMANAGER_CLASS;
	
	long nbCol = (long)pNbIndex+1;
	if(ARMLOCAL_GetData (vDatas, size / nbCol, nbCol, C_result) == ARM_OK)
	{
		Res = LocalMakeObjectId (C_result.getLong (), curClass);
	}

	_bstr_t tmpChaine = Res;

	*pRet = SysAllocString(tmpChaine);

	return S_OK;
}
catch(...)
{
	return E_FAIL;
}
}



STDMETHODIMP ActiveXModule::ARMINFYcMod(BSTR pYieldCurve,BSTR pInfCurve,BSTR *pRet)
{
	// TODO: Add your implementation code here

try
{
	ARM_result C_result;

	_bstr_t bYieldCurve(pYieldCurve);
	CCString l_YieldCurve = bYieldCurve;

	_bstr_t bInfCurve(pInfCurve);
	CCString l_InfCurve = bInfCurve;

	CCString Res;
	CCString curClass = LOCAL_INFCURVMODEL_CLASS;
	
	if(ARMLOCAL_infYCmod (LocalGetNumObjectId(l_YieldCurve), LocalGetNumObjectId(l_InfCurve), C_result) == ARM_OK)
	{
		Res = LocalMakeObjectId (C_result.getLong (), curClass);
	}

	_bstr_t tmpChaine = Res;

	*pRet = SysAllocString(tmpChaine);

	return S_OK;
}
catch(...)
{
	return E_FAIL;
}
}



STDMETHODIMP ActiveXModule::ARMAccrued(BSTR pSec,double pDate,BSTR pModel,double *pRet)
{
	// TODO: Add your implementation code here
try
{
	_bstr_t bSec(pSec);
	CCString l_Sec = bSec;

	_bstr_t bModel(pModel);
	CCString l_Model = bModel;

	long modelId;
	if (l_Model == "DEFAULT")
		modelId = ARM_NULL_OBJECT;
	else
		modelId = LocalGetNumObjectId(l_Model);

	ARM_result C_result;
	
	if(ARMLOCAL_ARM_Accrued (LocalGetNumObjectId(l_Sec), pDate, modelId, C_result) == ARM_OK)
	{
		*pRet = C_result.getDouble ();
	}

	return S_OK;
}
catch(...)
{
	return E_FAIL;
}
}


STDMETHODIMP ActiveXModule::ARMClonedAndSetNotional(BSTR bLegId,BSTR bAmortId,BSTR *pRet)
{
	// TODO: Add your implementation code here
try
{
	_bstr_t b_bLegId(bLegId);
	CCString l_bLegId = b_bLegId;

	_bstr_t b_bAmortId(bAmortId);
	CCString l_bAmortId = b_bAmortId;

	ARM_result C_result;
	CCString Res;
	CCString curClass = LOCAL_RESETMANAGER_CLASS;
	
	if(ARMLOCAL_ClonedAndSetNotional(LocalGetNumObjectId(l_bLegId),
									LocalGetNumObjectId(l_bAmortId), 
									100., 
									C_result) == ARM_OK)
	{
		Res = LocalMakeObjectId (C_result.getLong (), curClass);
	}

	_bstr_t tmpChaine = Res;

	*pRet = SysAllocString(tmpChaine);

	return S_OK;
}
catch(...)
{
	return E_FAIL;
}
}


STDMETHODIMP ActiveXModule::ARM_INF_GetZcFromSummit(BSTR Index,
											    BSTR Ccy,
												BSTR cvname,
												double date,
												BSTR seasonAdj,
												BSTR seasonAdjMode,
												BSTR *pRet)
{
	// TODO: Add your implementation code here
try
{
	_bstr_t b_Index(Index);
	CCString l_Index = b_Index;

	_bstr_t b_Ccy(Ccy);
	CCString l_Ccy = b_Ccy;

	_bstr_t b_cvname(cvname);
	CCString l_cvname = b_cvname;

	_bstr_t b_seasonAdj(seasonAdj);
	CCString l_seasonAdj = b_seasonAdj;

	_bstr_t b_seasonAdjMode(seasonAdjMode);
	CCString l_seasonAdjMode = b_seasonAdjMode;

	ARM_result C_result;
	CCString Res;
	CCString curClass = LOCAL_INFCURV_CLASS;

	int i_seasonAdjId = ARM_ConvYesOrNo (l_seasonAdj, C_result);
	int i_seasonAdjModeId = ARM_ConvSeasonalityMode (l_seasonAdjMode, C_result);
	
	if(ARMLOCAL_GetInfZcFromSummit_Create(l_Index,
										l_Ccy,
										l_cvname,
										date,
										i_seasonAdjId,
										i_seasonAdjModeId,
										C_result) == ARM_OK)
	{
		Res = LocalMakeObjectId (C_result.getLong (), curClass);
	}

	_bstr_t tmpChaine = Res;

	*pRet = SysAllocString(tmpChaine);

	return S_OK;
}
catch(...)
{
	return E_FAIL;
}
}

/*
STDMETHODIMP ActiveXModule::ARM_INF_CreateGenericLeg(double pStartDate,double pEndDate,BSTR pInfIdx,BSTR pRcvOrPay,BSTR pInterpType,double pLeverage,double pSpread,BSTR pResetFreq,BSTR pDaycount,BSTR pResetCal,BSTR pFwdRule,BSTR pIntRule,BSTR pStubRule,double pResetNumGap,double pResetDenomGap,BSTR pPayFreq,double pPayGap,BSTR pPayCal,BSTR pFinalNotionalType,double pFirstReset,BSTR *pRet)
{
try
{
	// TODO: Add your implementation code here
	// return
	ARM_result C_result;

	_bstr_t bInfIdx(pInfIdx);
	CCString l_InfIdx = bInfIdx;

	_bstr_t bRcvOrPay(pRcvOrPay);
	CCString l_RcvOrPay = bRcvOrPay;

	long receiveOrPayId;
	if((receiveOrPayId = ARM_ConvRecOrPay (l_RcvOrPay, C_result)) == ARM_DEFAULT_ERR)
		return S_FALSE;

	_bstr_t bInterpType(pInterpType);
	CCString l_InterpType = bInterpType;

	long interpTypeId;
	if((interpTypeId = ARM_ConvCPIDailyInterpMethod (l_InterpType, C_result)) == ARM_DEFAULT_ERR)
		return S_FALSE;

	_bstr_t bResetFreq(pResetFreq);
	CCString l_ResetFreq = bResetFreq;

	long resetFreqId;
	if((resetFreqId = ARM_ConvFrequency (l_ResetFreq, C_result)) == ARM_DEFAULT_ERR)
		return S_FALSE;
	
	_bstr_t bDaycount(pDaycount);
	CCString l_Daycount = bDaycount;

	long daycountId;
	daycountId = ARM_ConvDayCount (l_Daycount);

	_bstr_t bResetCal(pResetCal);
	CCString l_ResetCal = bResetCal;

	_bstr_t bFwdRule(pFwdRule);
	CCString l_FwdRule = bFwdRule;

	long fwdRuleId;
	if((fwdRuleId = ARM_ConvFwdRule (l_FwdRule, C_result)) == ARM_DEFAULT_ERR)
		return S_FALSE;

	_bstr_t bIntRule(pIntRule);
	CCString l_IntRule = bIntRule;

	long intRuleId;
	intRuleId = ARM_ConvIntRule (l_IntRule);

	_bstr_t bStubRule(pStubRule);
	CCString l_StubRule = bStubRule;

	long stubRuleId;
	stubRuleId = ARM_ConvStubRule (l_StubRule);

	if (pResetDenomGap == -1111)
		pResetDenomGap = pResetNumGap;

	_bstr_t bPayFreq(pPayFreq);
	CCString l_PayFreq = bPayFreq;

	long payFreqId;
	if((payFreqId = ARM_ConvFrequency (l_PayFreq, C_result)) == ARM_DEFAULT_ERR)
		return S_FALSE;

	_bstr_t bPayCal(pPayCal);
	CCString l_PayCal = bPayCal;

	_bstr_t bFinalNotionalType(pFinalNotionalType);
	CCString l_FinalNotionalType = bFinalNotionalType;

	long finalNotionalTypeId;
	if((finalNotionalTypeId = ARM_NotionalType (l_FinalNotionalType, C_result)) == ARM_DEFAULT_ERR)
		return S_FALSE;

	CCString Res;
	CCString curClass = LOCAL_INFSWAPLEG_CLASS;

	if(ARMLOCAL_InfLegAllInputs_Create (pStartDate, pEndDate, l_InfIdx, K_YEARTOYEAR_LEG,
										receiveOrPayId,interpTypeId,pLeverage,pSpread,
										resetFreqId,daycountId,l_ResetCal,fwdRuleId,
										intRuleId,stubRuleId,pResetNumGap,pResetDenomGap,
										payFreqId,pPayGap,l_PayCal,K_ADJUSTED,
										finalNotionalTypeId,pFirstReset,C_result) == ARM_OK)
	{
		Res = LocalMakeObjectId (C_result.getLong (), curClass);
	}

	_bstr_t tmpChaine = Res;

	*pRet = SysAllocString(tmpChaine);

	return S_OK;
}
catch(...)
{
	return E_FAIL;
}

}
*/


STDMETHODIMP ActiveXModule::ARMRiskyBond(double pIssueDate,double pMaturityDate,double pFirstCpnDate,double pCpnRate,double pRedemptionPrice,long pPeriodicity,VARIANT pDaycount,long pSettleGap,long pCpnDateFlag,BSTR pCcyId,double pRepo,double pSsl,double pRecoveryRate,BSTR *pRet)
{
try
{
	// TODO: Add your implementation code here
	// return
	ARM_result C_result;
	
	double dDaycount;

	if (VARIANT2Double(pDaycount,&dDaycount) == S_FALSE)
	{
		_bstr_t bDaycount(pDaycount);
		CCString l_Daycount = bDaycount;

		dDaycount = (double)ARM_ConvDayCount (l_Daycount);
	}

	_bstr_t bCcy(pCcyId);
	CCString l_Ccy = bCcy;

	long l_ccyId;

	if (l_Ccy == "DEFAULT")
		l_ccyId = -1;
	else
		l_ccyId = LocalGetNumObjectId (l_Ccy);

	CCString Res;
	CCString curClass = LOCAL_BOND_CLASS;

	if(ARMLOCAL_RiskyBond (pIssueDate, pMaturityDate,
						   pFirstCpnDate, pCpnRate,
						   pRedemptionPrice, pPeriodicity,
						   (long)dDaycount, pSettleGap,
						   pCpnDateFlag, l_ccyId,
						   pRepo, pSsl, pRecoveryRate, C_result) == ARM_OK)
	{
		Res = LocalMakeObjectId (C_result.getLong (), curClass);
	}

	_bstr_t tmpChaine = Res;

	*pRet = SysAllocString(tmpChaine);

	return S_OK;
}
catch(...)
{
	return E_FAIL;
}

}


STDMETHODIMP ActiveXModule::ARMRiskyBondWithCF(double pAsOfDate,double pRedemptionPrice,long pPeriodicity,VARIANT pDaycount,VARIANT *pYearTerms,VARIANT *pCashFlows,long pSettleGap,long pCpnDateFlag,BSTR pCcyId,double pRepo,double pSsl,double pRecoveryRate,BSTR *pRet)
{
try
{
	// TODO: Add your implementation code here
	// return
	ARM_result C_result;
	
	double dDaycount;

	if (VARIANT2Double(pDaycount,&dDaycount) == S_FALSE)
	{
		_bstr_t bDaycount(pDaycount);
		CCString l_Daycount = bDaycount;

		dDaycount = (double)ARM_ConvDayCount (l_Daycount);
	}

	_bstr_t bCcy(pCcyId);
	CCString l_Ccy = bCcy;

	long l_ccyId;

	if (l_Ccy == "DEFAULT")
		l_ccyId = -1;
	else
		l_ccyId = LocalGetNumObjectId (l_Ccy);

	VECTOR<double> vYearTerms;
	VECTOR<double> vCashFlows;
	long rateSize;

	if(VARIANT2VECTORDOUBLE (*pYearTerms,vYearTerms,rateSize) != S_OK)
		return S_FALSE;

	if(VARIANT2VECTORDOUBLE (*pCashFlows,vCashFlows,rateSize) != S_OK)
		return S_FALSE;

	CCString Res;
	CCString curClass = LOCAL_BOND_CLASS;

	if(ARMLOCAL_RiskyBondWithCF (pAsOfDate, pRedemptionPrice, pPeriodicity,
								 (long)dDaycount, pSettleGap,
								 pCpnDateFlag, l_ccyId ,
								 vYearTerms, vCashFlows, pRepo,
								 pSsl, pRecoveryRate, C_result) == ARM_OK)
	{
		Res = LocalMakeObjectId (C_result.getLong (), curClass);
	}

	_bstr_t tmpChaine = Res;

	*pRet = SysAllocString(tmpChaine);

	return S_OK;
}
catch(...)
{
	return E_FAIL;
}

}


STDMETHODIMP ActiveXModule::ARMGetFixing(BSTR source,
									 BSTR index,
									 BSTR term,
									 BSTR ccy,
									 double date,
									 double *pRet)
{
	// TODO: Add your implementation code here
try
{
	_bstr_t b_source(source);
	CCString s_source = b_source;

	_bstr_t b_index(index);
	CCString s_index = b_index;

	_bstr_t b_term(term);
	CCString s_term = b_term;

	_bstr_t b_ccy(ccy);
	CCString s_ccy = b_ccy;
	
	ARM_result C_result;

	if(ARMLOCAL_ARM_GetFixing(s_source,
							  s_index,
							  s_term,
							  s_ccy,
							  date,
							  C_result) == ARM_OK)
	{
		*pRet = C_result.getDouble ();
	}

	return S_OK;
}
catch(...)
{
	return E_FAIL;
}
}

STDMETHODIMP ActiveXModule::ARMGetFixingFromCalypso(BSTR source_,
									 BSTR index_,
									 BSTR term_,
									 BSTR ccy_,
									 BSTR curveName_,
									 DATE  date_,
									 VARIANT *pRet)
{
	// TODO: Add your implementation code here
	try
	{
		std::string source; VariantTools::convert(source_,source); 
		std::string index; VariantTools::convert(index_,index); 
		std::string term; VariantTools::convert(term_,term); 
		std::string ccy; VariantTools::convert(ccy_,ccy); 
		std::string curveName; VariantTools::convert(curveName_,curveName); 
		if (curveName=="") curveName="MO"; 
		ARM_Date date; VariantTools::convertXLDate(date_,date); 
		
		double value ; 
		ARM_CalypsoToolkit::GetFixing(index,term,ccy,source,curveName,date,value); 
		VariantTools::convert(value,*pRet); 	
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

STDMETHODIMP ActiveXModule::ARMHyperCube(VARIANT *pVolCurvId,
									 VARIANT* pKeys,
									 BSTR *pRet)
{
	// TODO: Add your implementation code here
try
{
	VECTOR<CCString> vVolCurveId;
	VECTOR<CCString> vKeys;
	VECTOR<long> vIds;
	long size;

	ARM_result C_result;

	if(VARIANT2VECTORCCSTRING (*pVolCurvId,vVolCurveId,size) != S_OK)
		return S_FALSE;

	if(VARIANT2VECTORCCSTRING (*pKeys,vKeys,size) != S_OK)
		return S_FALSE;

	for (int i=0; i < size; i++)
		vIds.push_back(LocalGetNumObjectId(vVolCurveId[i]));

	CCString Res;
	CCString curClass = LOCAL_HYPER_CUBE_CLASS;

	if(ARMLOCAL_HyperCube (vIds, vKeys, C_result) == ARM_OK)
	{
		Res = LocalMakeObjectId (C_result.getLong (), curClass);
	}

	_bstr_t tmpChaine = Res;

	*pRet = SysAllocString(tmpChaine);

	return S_OK;
}
catch(...)
{
	return E_FAIL;
}
}

STDMETHODIMP ActiveXModule::ARMIndexIndexCorrelCube(VARIANT* pVolCurveId,
									 VARIANT* pTenors1List,
									 VARIANT* pTenors2List, 
									 BSTR pInterSurfInterp,
									 BSTR *pRet)
{
	// TODO: Add your implementation code here
try
{
	VECTOR<CCString> vVolCurveId;
	VECTOR<CCString> vTenors1List;
	VECTOR<CCString> vTenors2List;
	VECTOR<long> vIds;

	long size;

	ARM_result C_result;

	if(VARIANT2VECTORCCSTRING (*pVolCurveId,vVolCurveId,size) != S_OK)
		return S_FALSE;

	if(VARIANT2VECTORCCSTRING (*pTenors1List,vTenors1List,size) != S_OK)
		return S_FALSE;

	if(VARIANT2VECTORCCSTRING (*pTenors2List,vTenors2List,size) != S_OK)
		return S_FALSE;

	for (int i=0; i < size; i++)
		vIds.push_back(LocalGetNumObjectId(vVolCurveId[i]));

	_bstr_t bInterSurfInterp(pInterSurfInterp);
	CCString C_InterSurfInterp = bInterSurfInterp;

	CCString Res;
	CCString curClass = LOCAL_HYPER_CUBE_CLASS;

	if(ARMLOCAL_IndexIndexCorrelCube (vIds, vTenors1List, vTenors2List, C_InterSurfInterp, C_result) == ARM_OK)
	{
		Res = LocalMakeObjectId (C_result.getLong (), curClass);
	}

	_bstr_t tmpChaine = Res;

	*pRet = SysAllocString(tmpChaine);

	return S_OK;
}
catch(...)
{
	return E_FAIL;
}
}




STDMETHODIMP ActiveXModule::ARMCmsLeg(double startDate,
					  double endDate,
					  BSTR cmsTypeId,
					  BSTR receiveOrPay,
					  BSTR yieldDecompFreq,
					  BSTR swapLegDayCount,
					  BSTR resetFreq,
					  BSTR payFreq,
					  long resetGap,		 
					  BSTR intRule,
					  BSTR ccyName,
					  BSTR resetTiming,
					  BSTR stubRule,
					  BSTR adjStartDate,
					  BSTR *pRet)
{
try
{
	ARM_result C_result;

	_bstr_t bcmsTypeId(cmsTypeId);
	CCString C_cmsTypeId = bcmsTypeId;
	long l_cmsTypeId = ARM_ConvIrType (C_cmsTypeId);

	_bstr_t breceiveOrPay(receiveOrPay);
	CCString C_receiveOrPay = breceiveOrPay;
	long l_receiveOrPay = ARM_ConvRecOrPay (C_receiveOrPay, C_result);

	_bstr_t byieldDecompFreq(yieldDecompFreq);
	CCString C_yieldDecompFreq = byieldDecompFreq;
	long l_yieldDecompFreq = ARM_ConvDecompFrequency(C_yieldDecompFreq);

	_bstr_t bswapLegDayCount(swapLegDayCount);
	CCString C_swapLegDayCount = bswapLegDayCount;
	long l_swapLegDayCount = ARM_ConvDayCount(C_swapLegDayCount);

	_bstr_t bresetFreq(resetFreq);
	CCString C_resetFreq = bresetFreq;
	long l_resetFreq = ARM_ConvFrequency (C_resetFreq,C_result);

	_bstr_t bpayFreq(payFreq);
	CCString C_payFreq = bpayFreq;
	long l_payFreq = ARM_ConvFrequency (C_payFreq,C_result);

	_bstr_t bintRule(intRule);
	CCString C_intRule = bintRule;
	long l_intRule = ARM_ConvIntRule (C_intRule);

	_bstr_t bccyName(ccyName);
	CCString C_ccyName = bccyName;
	
	_bstr_t bresetTiming(resetTiming);
	CCString C_resetTiming = bresetTiming;
	long l_resetTiming = ARM_ConvPayResetRule(C_resetTiming);

	_bstr_t bstubRule(stubRule);
	CCString C_stubRule = bstubRule;
	long l_stubRule = ARM_ConvStubRule(C_stubRule);

	_bstr_t badjStartDate(adjStartDate);
	CCString C_adjStartDate = badjStartDate;
	long l_adjStartDate = 1;
	if( C_adjStartDate != CCString("YES") ) 	{ l_adjStartDate = 0;}

	CCString Res;
	
	CCString prevClass;

	CCString curClass = LOCAL_SWAPLEG_CLASS;
	
	if(ARMLOCAL_CMSLEG (startDate,
					  endDate,
					  l_cmsTypeId,
					  l_receiveOrPay,
					  0,
					  0,
					  l_yieldDecompFreq,
					  l_swapLegDayCount,
					  l_resetFreq,
					  l_payFreq,
					  resetGap,		 
					  l_intRule,
					  false,
					  C_ccyName,
					  l_resetTiming,
					  l_stubRule,
					  l_adjStartDate,
					  C_result) == ARM_OK)
	{
		Res = LocalMakeObjectId (C_result.getLong (), curClass);
	}
	else
	{
	return S_OK;
	}

	_bstr_t tmpChaine = Res;

	*pRet = SysAllocString(tmpChaine);

	return S_OK;

}
catch(...)
{
	return E_FAIL;
}
}



STDMETHODIMP ActiveXModule::ARM_ReplicConvAdjust_Create(
										BSTR Payoff_ReplicMode,
										double Payoff_StepOrReplicPrecision,
										BSTR Payoff_StopMode,
										double Payoff_StopThreshold,
										BSTR Sensi_ReplicMode,
										double Sensi_StepOrReplicPrecision,
										BSTR Sensi_StopMode,
										double Sensi_StopThreshold,
										BSTR UsedModelId,
										double StrikeMinReplic,
										double StrikeMaxReplic,
										BSTR *pRet)
{
try
{
	ARM_result C_result;

	_bstr_t bPayoff_ReplicMode(Payoff_ReplicMode);
	CCString C_Payoff_ReplicMode = bPayoff_ReplicMode;
	long l_Payoff_ReplicMode = ARM_ConvReplicMode (C_Payoff_ReplicMode, C_result);

	_bstr_t bPayoff_StopMode(Payoff_StopMode);
	CCString C_Payoff_StopMode = bPayoff_StopMode;
	long l_Payoff_StopMode = ARM_ConvStopMode (C_Payoff_StopMode, C_result);

	_bstr_t bSensi_ReplicMode(Sensi_ReplicMode);
	CCString C_Sensi_ReplicMode = bSensi_ReplicMode;
	long l_Sensi_ReplicMode = ARM_ConvReplicMode (C_Sensi_ReplicMode, C_result);

	_bstr_t bSensi_StopMode(Sensi_StopMode);
	CCString C_Sensi_StopMode = bSensi_StopMode;
	long l_Sensi_StopMode = ARM_ConvStopMode (C_Sensi_StopMode, C_result);

	_bstr_t bUsedModelId(UsedModelId);
	CCString C_UsedModelId = bUsedModelId;
	long l_UsedModelId;
	
	if( (l_Payoff_ReplicMode== ARM_DEFAULT_ERR)
		|| (l_Payoff_StopMode == ARM_DEFAULT_ERR)
		|| (l_Sensi_ReplicMode == ARM_DEFAULT_ERR)
		|| (l_Sensi_StopMode == ARM_DEFAULT_ERR)
		)
		return S_FALSE;

	if (C_UsedModelId == "NULL")
		l_UsedModelId = ARM_NULL_OBJECT;
	else
		l_UsedModelId = LocalGetNumObjectId (C_UsedModelId);

	CCString Res;
	
	CCString curClass = LOCAL_REPLICCONVADJUST_CLASS;
	
	if(ARMLOCAL_ReplicConvAdjust_Create (l_Payoff_ReplicMode,
					Payoff_StepOrReplicPrecision,
					l_Payoff_StopMode,
					Payoff_StopThreshold,
					l_Sensi_ReplicMode,
					Sensi_StepOrReplicPrecision,
					l_Sensi_StopMode,
					Sensi_StopThreshold,
					l_UsedModelId,
					StrikeMinReplic,
					StrikeMaxReplic,
					C_result) == ARM_OK)
	{
		Res = LocalMakeObjectId (C_result.getLong (), curClass);
	}
	else
	{
	return S_OK;
	}

	_bstr_t tmpChaine = Res;

	*pRet = SysAllocString(tmpChaine);

	return S_OK;

}
catch(...)
{
	return E_FAIL;
}
}

STDMETHODIMP ActiveXModule::ARM_MapConvAdjust_Create(
										BSTR LiborArrearAdj,
										BSTR NaturalCMSAdj,
										BSTR PaymentLagAdj,
										BSTR *pRet)
{
try
{
	ARM_result C_result;

	_bstr_t bLiborArrearAdj(LiborArrearAdj);
	CCString C_LiborArrearAdj = bLiborArrearAdj;
	long l_LiborArrearAdj;

	if (C_LiborArrearAdj == "NULL")
		l_LiborArrearAdj = ARM_NULL_OBJECT;
	else
		l_LiborArrearAdj = LocalGetNumObjectId (C_LiborArrearAdj);

	_bstr_t bNaturalCMSAdj(NaturalCMSAdj);
	CCString C_NaturalCMSAdj = bNaturalCMSAdj;
	long l_NaturalCMSAdj;

	if (C_NaturalCMSAdj == "NULL")
		l_NaturalCMSAdj = ARM_NULL_OBJECT;
	else
		l_NaturalCMSAdj = LocalGetNumObjectId (C_NaturalCMSAdj);

	_bstr_t bPaymentLagAdj(PaymentLagAdj);
	CCString C_PaymentLagAdj = bPaymentLagAdj;
	long l_PaymentLagAdj;

	if (C_PaymentLagAdj == "NULL")
		l_PaymentLagAdj = ARM_NULL_OBJECT;
	else
		l_PaymentLagAdj = LocalGetNumObjectId (C_PaymentLagAdj);

	CCString Res;
	
	CCString curClass = LOCAL_MAPCONVADJUST_CLASS;
	
	if(ARMLOCAL_MapConvAdjust_Create (l_LiborArrearAdj,
					l_NaturalCMSAdj,
					l_PaymentLagAdj,
					C_result,-1) == ARM_OK)
	{
		Res = LocalMakeObjectId (C_result.getLong (), curClass);
	}
	else
	{
	return S_OK;
	}

	_bstr_t tmpChaine = Res;

	*pRet = SysAllocString(tmpChaine);

	return S_OK;

}
catch(...)
{
	return E_FAIL;
}
}


STDMETHODIMP ActiveXModule::ARMCreateGenCorrelatorManager(VARIANT *pMktTags,
									 VARIANT* pHyperDiagVol,
									 VARIANT *pIndexIndexVol,
									 VARIANT* pCorrelVol,
									 VARIANT* pIndexVol,
									 BSTR *pRet)
{
	// TODO: Add your implementation code here
try
{	
	VECTOR<CCString> C_mktTags;
	VECTOR<CCString> C_HyperDiagVol;
	VECTOR<CCString> C_IndexIndexVol;
	VECTOR<CCString> C_CorrelVol;
	VECTOR<CCString> C_IndexVol;

	VECTOR<CCString> vVolCurveId;
	VECTOR<CCString> vKeys;

	VECTOR<long> C_HyperDiagVolIds;
	VECTOR<long> C_IndexIndexVolIds;
	VECTOR<long> C_CorrelVolIds;
	VECTOR<long> C_IndexVolIds;

	long size;

	ARM_result C_result;

	if(VARIANT2VECTORCCSTRING (*pMktTags,C_mktTags,size) != S_OK)
		return S_FALSE;

	if(VARIANT2VECTORCCSTRING (*pHyperDiagVol,C_HyperDiagVol,size) != S_OK)
		return S_FALSE;

	for (int i=0; i < size; i++)
		C_HyperDiagVolIds.push_back(LocalGetNumObjectId(C_HyperDiagVol[i]));


	if(VARIANT2VECTORCCSTRING (*pIndexIndexVol,C_IndexIndexVol,size) != S_OK)
		return S_FALSE;

	for (i=0; i < size; i++)
		C_IndexIndexVolIds.push_back(LocalGetNumObjectId(C_IndexIndexVol[i]));

	if(VARIANT2VECTORCCSTRING (*pCorrelVol,C_CorrelVol,size) != S_OK)
		return S_FALSE;

	for (i=0; i < size; i++)
		C_CorrelVolIds.push_back(LocalGetNumObjectId(C_CorrelVol[i]));

	if(VARIANT2VECTORCCSTRING (*pIndexVol,C_IndexVol,size) != S_OK)
		return S_FALSE;
	for ( i=0; i < size; i++)
		C_IndexVolIds.push_back(LocalGetNumObjectId(C_IndexVol[i]));

	CCString Res;
	CCString curClass = LOCAL_CORRELMANAGER_CLASS;

	if (0)// ARMLOCAL_CreateGenCorrelatorManager(C_mktTags,
		// 									C_HyperDiagVolIds,
		// 									C_IndexIndexVolIds,
		// 									C_CorrelVolIds,
		// 									C_IndexVolIds,
		// 									C_result) == ARM_OK)
	{
		Res = LocalMakeObjectId (C_result.getLong (), curClass);
	}

	_bstr_t tmpChaine = Res;

	*pRet = SysAllocString(tmpChaine);

	return S_OK;



}
catch(...)
{
	return E_FAIL;
}
}


STDMETHODIMP ActiveXModule::ARMGetResetMgrFromSummit(double pAsOf,
													  BSTR pIndex,
													  BSTR pSource,
													  BSTR pCcy,
													  BSTR pIsInflationIndex,
													  BSTR pTerm,
													  BSTR *pRet)
{
	// TODO: Add your implementation code here
	ARM_result C_result;
try
{

	_bstr_t bIndex(pIndex);
	CCString C_Index = bIndex;

	_bstr_t bSource(pSource);
	CCString C_Source = bSource;

	_bstr_t bCcy(pCcy);
	CCString C_Ccy = bCcy;

	_bstr_t bIsInflationIndex(pIsInflationIndex);
	CCString C_IsInflationIndex = bIsInflationIndex;
	long IsInflationIndexId;
	
	if((IsInflationIndexId = ARM_ConvYesOrNo (C_IsInflationIndex, C_result)) == ARM_DEFAULT_ERR)
		return S_FALSE;

	_bstr_t bTerm(pTerm);
	CCString C_Term = bTerm;

	CCString Res;
	CCString curClass = LOCAL_RESETMANAGER_CLASS;

	if (ARMLOCAL_GetResetMgrFromSummit(pAsOf,
									   C_Index,
									   C_Source,
									   C_Ccy,
									   IsInflationIndexId,
									   C_Term,
									   C_result) == ARM_OK)
	{
		Res = LocalMakeObjectId (C_result.getLong (), curClass);
	}
	else
	{
		return createErrorInfo("",CCSTringToSTLString(C_result.getMsg())); 
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


STDMETHODIMP ActiveXModule::ARMGetReset(BSTR pResetMgr, double pDate, double *pRet)
{
	// TODO: Add your implementation code here

try
{
	ARM_result C_result;

	*pRet = ARM_NULL_OBJECT;

	_bstr_t bResetMgr(pResetMgr);
	CCString C_ResetMgr = bResetMgr;

	if(ARMLOCAL_GetReset(LocalGetNumObjectId(C_ResetMgr), pDate, C_result) == ARM_OK)
	{
		*pRet = C_result.getDouble ();
	}

	return S_OK;
}
catch(...)
{
	return E_FAIL;
}
}


STDMETHODIMP ActiveXModule::ARMSetLastFixing(BSTR pSecurityId, 
												double pRate, 
												double pAsOf, 
												double pBeforeLastFxingDate, 
												double pResetDate, 
												BSTR *pRet)
{
	// TODO: Add your implementation code here

try
{
	ARM_result C_result;

	CCString Res;

	_bstr_t bSecurityId(pSecurityId);
	CCString C_SecurityId = bSecurityId;

	CCString curClass = LocalGetStringObjectClass(C_SecurityId);

	if(ARMLOCAL_ARM_SetLastFixing(LocalGetNumObjectId(C_SecurityId), 
									pRate,
									pBeforeLastFxingDate, 
									pAsOf,
									pResetDate,
									C_result) == ARM_OK)
	{
		Res = LocalMakeObjectId (C_result.getLong (), curClass);
	}

	_bstr_t tmpChaine = Res;

	*pRet = SysAllocString(tmpChaine);

	return S_OK;
}
catch(...)
{
	return E_FAIL;
}
}


STDMETHODIMP ActiveXModule::ARMGenAmortization(BSTR pSwaplegId,
											   BSTR pAmortMethod,
											   BSTR pAmortFreq,
											   double pAmortAmount,
											   BSTR pDaycount,
											   double pLegNotional,
											   double pAmortRate,
											   double pReducedMaturity,
											   BSTR pModelId,
											   double pCleanUp,
											   BSTR *pRet)
{
	try
{
	ARM_result C_result;

	CCString Res;

	_bstr_t bSwaplegId(pSwaplegId);
	CCString C_SwaplegId = bSwaplegId;

	_bstr_t bAmortMethod(pAmortMethod);
	CCString C_AmortMethod = bAmortMethod;

	_bstr_t bAmortFreq(pAmortFreq);
	CCString C_AmortFreq = bAmortFreq;

	_bstr_t bDaycount(pDaycount);
	CCString C_Daycount = bDaycount;

	_bstr_t bModelId(pModelId);
	CCString C_ModelId = bModelId;

	long dayCountId = ARM_ConvDayCount (C_Daycount);

	long amortMethodId;
	if((amortMethodId = ARM_ConvAmortMethod (C_AmortMethod, C_result)) == ARM_DEFAULT_ERR)
		return S_FALSE;

	long amortFreqId;
	if((amortFreqId = ARM_ConvFrequency (C_AmortFreq, C_result)) == ARM_DEFAULT_ERR)
		return S_FALSE;

	long modelId;
	if(C_ModelId == "DEFAULT")
	{
		modelId = ARM_NULL_OBJECT;
	}
	else
	{
		modelId = LocalGetNumObjectId (C_ModelId);
	}

	CCString curClass = LOCAL_REFVAL_CLASS;

	if(ARMLOCAL_ARM_GenAmortization(LocalGetNumObjectId(C_SwaplegId),
									amortMethodId,
									amortFreqId,
									pAmortAmount,
									dayCountId,
									pLegNotional,
									pAmortRate,
									pReducedMaturity,
									modelId,
									pCleanUp,
									C_result) == ARM_OK)
	{
		Res = LocalMakeObjectId (C_result.getLong (), curClass);
	}

	_bstr_t tmpChaine = Res;

	*pRet = SysAllocString(tmpChaine);

	return S_OK;
}
catch(...)
{
	return E_FAIL;
}
}


STDMETHODIMP ActiveXModule::ARMCptRefvalue(BSTR pRefvalId, double pDate, double *pRet)
{
	// TODO: Add your implementation code here

try
{
	ARM_result C_result;

	*pRet = ARM_NULL_OBJECT;

	_bstr_t bRefvalId(pRefvalId);
	CCString C_RefvalId = bRefvalId;

	if(ARMLOCAL_CptRefValue(LocalGetNumObjectId(C_RefvalId), pDate, C_result) == ARM_OK)
	{
		*pRet = C_result.getDouble ();
	}

	return S_OK;
}
catch(...)
{
	return E_FAIL;
}
}

STDMETHODIMP ActiveXModule::ARMCalypsoDevConnect()
{
	try 
	{
		char* user = getenv ("USERNAME");
		std::stringstream sstr ;
		sstr << "-user "<<user<<" -password "<<user<<" -env calypso_dev" ;
		ARM_CalypsoToolkit::init(sstr.str()); 
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
	return S_OK;	
}
STDMETHODIMP ActiveXModule::ARMCalypsoProdConnect()
{
	try 
	{
		char* user_ = getenv ("USERNAME");
		std::string user; 
		if(user_)user=user_;else user="unknown"; 

		std::stringstream sstr ;
		sstr << "-user "<<user<<" -password "<<user<<" -env calypso_prod" ;
		ARM_CalypsoToolkit::init(sstr.str()); 
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
	return S_OK;
}
STDMETHODIMP ActiveXModule::ARMCalypsoRecConnect()
{
	try 
	{
		char* user = getenv ("USERNAME");
		std::stringstream sstr ;
		sstr << "-user "<<user<<" -password "<<user<<" -env calypso_rec" ;
		ARM_CalypsoToolkit::init(sstr.str()); 
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
	return S_OK;
}

STDMETHODIMP ActiveXModule::Local_ARM_ProdConnect()
{
	try 
	{
		if (GetDataRetrieverVersion() == WSETKRETRIEVER)
		{
			connection_wsetoolkit("PROD","OTC","wsotc");
		}
		else
		{
			switchToETK();

			connection_etoolkit(SUMMIT_PROD_CONNEXION_USERNAME,
								SUMMIT_PROD_CONNEXION_PASSWD,
								SUMMIT_PROD_CONNEXION_CONTEXT,
								SUMMIT_PROD_IT_CONFIG_DOMAINSDIR,
								SUMMIT_PROD_IT_DOMAIN_NAME);
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
	return S_OK;
}

STDMETHODIMP ActiveXModule::ARM_SetDefaultCurrency (BSTR isoCCy, VARIANT *pRet)
{
	// return
	ARM_result C_result;
	
	try {
		// C variable
		string currency;
		VariantTools::convert(isoCCy, currency);
		long retCode = ARMLOCAL_ARM_SetDefaultCurrency (CCString(currency.c_str()), C_result);
		string newLabel = "OK";
		if (retCode == -1) newLabel = "KO"; //LocalPersistent::get().getStringId(retCode, LOCAL_CCY_CLASS);
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
	return S_OK;
}


STDMETHODIMP ActiveXModule::ARMLivretALeg(double pStartDate,
											double pEndDate,
											BSTR pRcvOrPay,
											VARIANT pSpread,
											BSTR pResetFreq,
											BSTR pPayFreq,
											BSTR pResetTiming,
											BSTR pPayTiming,
											BSTR pCcy,
											BSTR pIntRule,
											double pResetGap,
											BSTR pResetCal,
											BSTR pPayCal,
											double pDecompPricingFlag,
											BSTR pNxChange,
											BSTR pStubRule,
											double pRefDate,
											BSTR pAdjStartDate,
											BSTR pDayCount,
											BSTR *pRet)
{
	// TODO: Add your implementation code here
try
{
	ARM_result C_result;

	_bstr_t bRcvOrPay(pRcvOrPay);
	CCString l_RcvOrPay = bRcvOrPay;

	_bstr_t bResetFreq(pResetFreq);
	CCString l_ResetFreq = bResetFreq;

	_bstr_t bPayFreq(pPayFreq);
	CCString l_PayFreq = bPayFreq;

	_bstr_t bResetTiming(pResetTiming);
	CCString l_ResetTiming = bResetTiming;

	_bstr_t bPayTiming(pPayTiming);
	CCString l_PayTiming = bPayTiming;

	_bstr_t bCcy(pCcy);
	CCString l_Ccy = bCcy;

	_bstr_t bIntRule(pIntRule);
	CCString l_IntRule = bIntRule;

	_bstr_t bResetCal(pResetCal);
	CCString l_ResetCal = bResetCal;

	_bstr_t bPayCal(pPayCal);
	CCString l_PayCal = bPayCal;

	_bstr_t bNxChange(pNxChange);
	CCString l_NxChange = bNxChange;

	_bstr_t bStubRule(pStubRule);
	CCString l_StubRule = bStubRule;

	_bstr_t bAdjStartDate(pAdjStartDate);
	CCString l_AdjStartDate = bAdjStartDate;

	_bstr_t bDayCount(pDayCount);
	CCString l_DayCount = bDayCount;

	long receiveOrPayId;
	if((receiveOrPayId = ARM_ConvRecOrPay (l_RcvOrPay, C_result)) == ARM_DEFAULT_ERR)
		return S_FALSE;

	long payFreqId; 
	if((payFreqId = ARM_ConvFrequency (l_PayFreq, C_result)) == ARM_DEFAULT_ERR)
		return S_FALSE;

	long resetFreqId; 
	if((resetFreqId = ARM_ConvFrequency (l_ResetFreq, C_result)) == ARM_DEFAULT_ERR)
		return S_FALSE;

	long payTimingId = ARM_ConvPayResetRule (l_PayTiming);
	long resetTimingId = ARM_ConvPayResetRule (l_ResetTiming);

	bool ccyIsObject = false;
	if (( l_Ccy.GetLen() > 3 ) && ( !(l_Ccy == "DEFAULT")) )
	   ccyIsObject = true;

	long intRuleId = ARM_ConvIntRule (l_IntRule);
	long nxChange = ARM_NotionalExchange(l_NxChange);
	long stubRuleId = ARM_ConvStubRule (l_StubRule);

	double dSpread;
	long   spreadType;
	if (VARIANT2Double(pSpread,&dSpread) == S_FALSE)
	{
		_bstr_t bSpread(pSpread);
		CCString l_Spread = bSpread;

		dSpread = (double) LocalGetNumObjectId(l_Spread);
		spreadType = 1L;
	}
	else
	{
	   spreadType = 0L;
	}

	long adjStartDateId; 
	if ((adjStartDateId = ARM_ConvYesOrNo(l_AdjStartDate, C_result)) == ARM_DEFAULT_ERR)
		return S_FALSE;

	long daycountId = ARM_ConvDayCount(l_DayCount);

	CCString Res;
	CCString curClass = LOCAL_SWAPLEG_CLASS;

	if(ARMLOCAL_LIVRETALEG (pStartDate,
							pEndDate,
							K_LIVRET_A,
							receiveOrPayId,
							spreadType,
							dSpread,
							resetFreqId,
							payFreqId,
							resetTimingId,
							payTimingId,
							ccyIsObject,
							l_Ccy,
							intRuleId,
							(long) pResetGap,
							l_ResetCal,
							l_PayCal,
							(long) pDecompPricingFlag,
							nxChange,
							stubRuleId,
							pRefDate,
							adjStartDateId,
							daycountId,
							C_result) == ARM_OK)
	{
		Res = LocalMakeObjectId (C_result.getLong (), curClass);
	}

	_bstr_t tmpChaine = Res;

	*pRet = SysAllocString(tmpChaine);

	return S_OK;
}
catch(...)
{
	return E_FAIL;
}
}


STDMETHODIMP ActiveXModule::ARMLivretACurve(double pAsOf,
											BSTR pInfCurvId,
											BSTR pEuribCurvId,
											double pFlagRouding,
											BSTR pInfResetMgrId,
											BSTR pFixingLivretAId,
											BSTR pFixingEuribId,
											BSTR pMonthForAugust,
											BSTR pMonthForFebruary,
											BSTR *pRet)
{
	// TODO: Add your implementation code here
try
{
	ARM_result C_result;

	_bstr_t bInfCurvId(pInfCurvId);
	CCString l_InfCurvId = bInfCurvId;

	_bstr_t bEuribCurvId(pEuribCurvId);
	CCString l_EuribCurvId = bEuribCurvId;

	_bstr_t bInfResetMgrId(pInfResetMgrId);
	CCString l_InfResetMgrId = bInfResetMgrId;

	_bstr_t bFixingLivretAId(pFixingLivretAId);
	CCString l_FixingLivretAId = bFixingLivretAId;

	_bstr_t bFixingEuribId(pFixingEuribId);
	CCString l_FixingEuribId = bFixingEuribId;

	_bstr_t bMonthForAugust(pMonthForAugust);
	CCString l_MonthForAugust = bMonthForAugust;

	_bstr_t bMonthForFebruary(pMonthForFebruary);
	CCString l_MonthForFebruary = bMonthForFebruary;

	long monthForAugustId = ARM_ConvMonth(l_MonthForAugust);
	long monthForFebruaryId = ARM_ConvMonth(l_MonthForFebruary);


	CCString Res;
	CCString curClass = LOCAL_LIVRET_A_CURVE_CLASS;

	if(ARMLOCAL_LIVRETACURVE (pAsOf,
							  LocalGetNumObjectId (l_InfCurvId),
							  LocalGetNumObjectId (l_EuribCurvId),
							  pFlagRouding,
							  LocalGetNumObjectId(l_InfResetMgrId),
							  LocalGetNumObjectId(l_FixingLivretAId),
							  LocalGetNumObjectId(l_FixingEuribId), 
							  monthForAugustId,
							  monthForFebruaryId,
							  C_result) == ARM_OK)
	{
		Res = LocalMakeObjectId (C_result.getLong (), curClass);
	}

	_bstr_t tmpChaine = Res;

	*pRet = SysAllocString(tmpChaine);

	return S_OK;
}
catch(...)
{
	return E_FAIL;
}
}

STDMETHODIMP ActiveXModule::ARMFutDelivery(BSTR pFut, BSTR pCcy, double *pRet)
{
try
{
	// TODO: Add your implementation code here
	ARM_result C_result;

	_bstr_t bFut(pFut);
	CCString l_Fut = bFut;

	_bstr_t bCcy(pCcy);
	CCString l_Ccy = bCcy;

	double dateRes;

	if(ARMLOCAL_FutDelivery(l_Fut,l_Ccy,C_result) == ARM_OK)
	{
		dateRes = Local_ARMDATE2XLDATE(C_result.getString ());
	}

	*pRet = dateRes;

	return S_OK;
}
catch(...)
{
	return E_FAIL;
}

}

STDMETHODIMP ActiveXModule::ARMInfCurveSetResetMgr(BSTR pInfCurve, BSTR pResetMgr,BSTR *pRet)
{
try
{
	// TODO: Add your implementation code here
	ARM_result C_result;

	_bstr_t bInfCurve(pInfCurve);
	CCString l_InfCurve = bInfCurve;

	_bstr_t bResetMgr(pResetMgr);
	CCString l_ResetMgr = bResetMgr;

	CCString Res;
	CCString curClass = LOCAL_INFCURV_CLASS;

	if(ARMLOCAL_InfCurv_SetResetManager(LocalGetNumObjectId(l_InfCurve),
										LocalGetNumObjectId(l_ResetMgr),
										C_result) == ARM_OK)
	{
		Res = LocalMakeObjectId (C_result.getLong (), curClass);
	}

	_bstr_t tmpChaine = Res;

	*pRet = SysAllocString(tmpChaine);

	return S_OK;
}
catch(...)
{
	return E_FAIL;
}
}

STDMETHODIMP ActiveXModule::ARMSetCalendar(BSTR pFileName,double *pRet)
{
try
{
	// TODO: Add your implementation code here
	ARM_result C_result;

	_bstr_t bFileName(pFileName);
	CCString l_FileName = bFileName;

	if (TMP_DATE.GetCalendar())
		delete TMP_DATE.GetCalendar();

	TMP_DATE.InitCalendar((char*)l_FileName,false);

	*pRet = 1.;

	return S_OK;
}
catch(...)
{
	return E_FAIL;
}
}


STDMETHODIMP ActiveXModule::ARMInitGigaSpaces(BSTR pUrl,double *pRet)
{
try
{
	// TODO: Add your implementation code here
	ARM_result C_result;

	_bstr_t bUrl(pUrl);
	CCString l_Url = bUrl;

	ARM_XGigaToolKit::init((string)l_Url);

	*pRet = 1.;

	return S_OK;
}
catch(...)
{
	return E_FAIL;
}
}

STDMETHODIMP ActiveXModule::ARMBermudanXStyle(VARIANT *pxDates,
											  VARIANT *pexpiryDates,
											  BSTR *pRet)
{
	try
	{
		ARM_result C_result;

		CCString Res;

		long xDatesSize;
		long expiryDatesSize;
		
		VECTOR<double> vxDates;
		VECTOR<double> vexpiryDates;

		CCString defaultExpiryDate("");
		
		CCString curClass = LOCAL_XSTYLE_CLASS;

		if(VARIANT2VECTORDOUBLE (*pxDates,vxDates,xDatesSize) != S_OK)
			return S_FALSE;

		try
		{
			if(VARIANT2CCString (*pexpiryDates,defaultExpiryDate) != S_OK)
				return S_FALSE;
		}
		catch (...)
		{
			if(VARIANT2VECTORDOUBLE (*pexpiryDates,vexpiryDates,expiryDatesSize) != S_OK)
				return S_FALSE;
		}

		if(strcmp((const char*)defaultExpiryDate,"") != 0)
		{
			if(VARIANT2VECTORDOUBLE (*pexpiryDates,vexpiryDates,expiryDatesSize) != S_OK)
				return S_FALSE;
		}


		if(ARMLOCAL_BERMUDANXSTYLE(vxDates,vexpiryDates,C_result) == ARM_OK)
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
	return S_OK;
}

STDMETHODIMP ActiveXModule::ARMSetDiscountPricingMode(BSTR pModelId,
													  int pDiscountPricingMode,
													  BSTR *pRet)
{
	try
	{
		ARM_result C_result;

		CCString Res;

		_bstr_t bModelId(pModelId);
		CCString l_ModelId = bModelId;

		CCString curClass = LocalGetStringObjectClass(l_ModelId);

		if(ARMLOCAL_SetDiscountPricingMode(LocalGetNumObjectId(l_ModelId),
										   pDiscountPricingMode,
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
	return S_OK;
}

STDMETHODIMP ActiveXModule::ARMPF(VARIANT *pinsts,
								  VARIANT *pcoeffs,
								  VARIANT *pmarketPrices,
								  VARIANT *pprecisions,
								  BSTR *pRet)
{
	try
	{
		ARM_result C_result;

		CCString Res;

		long instsSize;
		long coeffsSize;
		long marketPricesSize;
		long precisionsSize;
				
		VECTOR<CCString> vinsts;
		VECTOR<double> vcoeffs;
		VECTOR<double> vmarketPrices;
		VECTOR<double> vprecisions;

	
		CCString curClass = LOCAL_PF_CLASS;

		if(VARIANT2VECTORCCSTRING (*pinsts,vinsts,instsSize) != S_OK)
			return S_FALSE;

		if(VARIANT2VECTORDOUBLE (*pcoeffs,vcoeffs,coeffsSize) != S_OK)
			return S_FALSE;

		if(VARIANT2VECTORDOUBLE (*pmarketPrices,vmarketPrices,marketPricesSize) != S_OK)
			return S_FALSE;
	
		
		if(VARIANT2VECTORDOUBLE (*pprecisions,vprecisions,precisionsSize) != S_OK)
				return S_FALSE;


		VECTOR<long> vinstsId;
		long sz = vinsts.size ();

		for (int i = 0; i < sz; i++)
		{
			vinstsId.push_back (LocalGetNumObjectId (vinsts[i]));
		}

	
		if(ARMLOCAL_PF(vinstsId,vcoeffs,vmarketPrices,vprecisions,C_result) == ARM_OK)
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
	return S_OK;
}


STDMETHODIMP ActiveXModule::ARMBondTEC(double pIssueDate,
									   double pMaturityDate,
									   double pFirstCpnDate,
									   double pCpnRate,
									   double pRedempPrice,
									   long pPeriodicity,
									   VARIANT pDaycount,
									   long pSettleGap,
									   long pCpnDateFlag,
									   BSTR pCcyId,
									   double ptec,
									   BSTR pPFTecId,
									   BSTR pModTecId,
									   BSTR *pRet)


{
try
{
	// TODO: Add your implementation code here
	// return
	ARM_result C_result;

	double dDaycount;

	if (VARIANT2Double(pDaycount,&dDaycount) == S_FALSE)
	{
		_bstr_t bDaycount(pDaycount);
		CCString l_Daycount = bDaycount;

		dDaycount = (double)ARM_ConvDayCount( l_Daycount);
	}

	_bstr_t bCcy(pCcyId);
	CCString l_Ccy = bCcy;

	_bstr_t bPF(pPFTecId);
	CCString l_PF = bPF;

	_bstr_t bMod(pModTecId);
	CCString l_Mod = bMod;


	long l_ccyId;
	long l_pfId;
	long l_modId;
	
	if (l_Ccy == "DEFAULT")
		l_ccyId = -1;
	else
		l_ccyId = LocalGetNumObjectId (l_Ccy);

	if (l_PF == "DEFAULT")
		l_pfId = -1;
	else
		l_pfId = LocalGetNumObjectId (l_PF);

	if (l_Mod == "DEFAULT")
		l_modId = -1;
	else
		l_modId = LocalGetNumObjectId (l_Mod);

	CCString Res;
	CCString curClass = LOCAL_BOND_CLASS;

	if(ARMLOCAL_bondTEC (pIssueDate, pMaturityDate, pFirstCpnDate, pCpnRate,
							pRedempPrice, pPeriodicity,
							(long)dDaycount, pSettleGap,
							pCpnDateFlag, l_ccyId , ptec, l_pfId,l_modId, C_result) == ARM_OK)
	{
		Res = LocalMakeObjectId (C_result.getLong (), curClass);
	}

	_bstr_t tmpChaine = Res;

	*pRet = SysAllocString(tmpChaine);

	return S_OK;
}
catch(...)
{
	return E_FAIL;
}

}


STDMETHODIMP ActiveXModule::ARMPFModFit(BSTR pmodName,
										BSTR ppf,
										double psettlement,
										BSTR pzc,
										VARIANT *pvList,
										VARIANT *pfList,
										long nag_algo,
										long pstep,
										double phorizon,
										BSTR *pRet)

{
try
{
	// TODO: Add your implementation code here
	// return
	ARM_result C_result;

	 double C_settlement;

	VECTOR<double> vvList;
	VECTOR<double> vfList;

	long vListSize;
	long fListSize;

	_bstr_t bmodName(pmodName);
	CCString l_modname = bmodName;
	
	if(VARIANT2VECTORDOUBLE (*pvList,vvList,vListSize) != S_OK)
		return S_FALSE;

	if(VARIANT2VECTORDOUBLE (*pfList,vfList,fListSize) != S_OK)
		return S_FALSE;

	long l_zcId;

	_bstr_t bzc(pzc);
	CCString l_zc = bzc;

	_bstr_t bpf(ppf);
	CCString l_pf = bpf;

	if (l_zc == "DEFAULT")
		l_zcId = ARM_NULL_OBJECT;
	else
		l_zcId = LocalGetNumObjectId (l_zc);

	CCString Res;
	CCString curClass;
    
    if (l_modname == "GYCM")
    {
        curClass = LOCAL_GYCMODEL_CLASS;
    }
    else
    if (l_modname == "HW2F")
    {
        curClass = LOCAL_HW2FMODEL_CLASS;
    }
    else
    if (l_modname == "ZCVS")
    {
        curClass = LOCAL_YIELD_CURVE_BASIC_CLASS;
    }
    else
    if (l_modname == "BKTR")
    {
        curClass = LOCAL_BKIRTREE_CLASS;
    }
    else
    if (l_modname == "ZCSP")
    {
        curClass = LOCAL_YIELD_CURVE_BASIC_CLASS;
    }
    else
    if (l_modname == "GYLS")
    {
        curClass = LOCAL_GYCLSMODEL_CLASS;
    }
    else
    {
        return S_FALSE;
    }

	if(ARMLOCAL_PFMODFIT(l_modname,LocalGetNumObjectId(l_pf), psettlement, l_zcId, 
						vvList, vfList, (long)pstep, phorizon, nag_algo , C_result) == ARM_OK)
	{
		Res = LocalMakeObjectId (C_result.getLong (), curClass);
	}

	_bstr_t tmpChaine = Res;

	*pRet = SysAllocString(tmpChaine);

	return S_OK;
}
catch(...)
{
	return E_FAIL;
}
}


STDMETHODIMP ActiveXModule::ARMTMLeg(BSTR ptmIxType, 
									   double pstartDate, 
									   double pendDate, 
									   BSTR pPorR, 
									   double pspread,
									   BSTR ppayFrequency, 
									   BSTR presetFrequency, 
									   BSTR pinterestRule, 
									   BSTR pfwdRule, 
									   BSTR pstubRule, 
									   BSTR pccy,
									   BSTR *pRet)

{
try
{

	ARM_result C_result;
	CCString Res;

	long tmIxTypeId;
	_bstr_t btmIxType(ptmIxType);
	CCString l_tmIxType = btmIxType;

	tmIxTypeId = ARM_ConvIrType(l_tmIxType);

	long PorRId;
	_bstr_t bPorR(pPorR);
	CCString l_PorR = bPorR;

	if((PorRId = ARM_ConvRecOrPay (l_PorR, C_result)) == ARM_DEFAULT_ERR)
	return S_FALSE;

	long payFrequencyId;
	_bstr_t bpayFrequency(ppayFrequency);
	CCString l_payFrequency = bpayFrequency;

	if((payFrequencyId = ARM_ConvFrequency (l_payFrequency, C_result)) == ARM_DEFAULT_ERR)
		return S_FALSE;

	long resetFrequencyId;
	_bstr_t bresetFrequency(presetFrequency);
	CCString l_resetFrequency = bresetFrequency;

 	if((resetFrequencyId = ARM_ConvFrequency (l_resetFrequency, C_result)) == ARM_DEFAULT_ERR)
		return S_FALSE;

	long interestRuleId;
	_bstr_t binterestRule(pinterestRule);
	CCString l_interestRule = binterestRule;

	interestRuleId = ARM_ConvIntRule (l_interestRule);

	long fwdRuleId;
	_bstr_t bfwdRule(pfwdRule);
	CCString l_fwdRule = bfwdRule;

	if((fwdRuleId = ARM_ConvFwdRule (l_fwdRule, C_result)) == ARM_DEFAULT_ERR)
		return S_FALSE;

	long stubRuleId;
	_bstr_t bstubRule(pstubRule);
	CCString l_stubRule = bstubRule;

	stubRuleId = ARM_ConvStubRule (l_stubRule);

	long ccyId;
	bool ccyIsObject = false;
	_bstr_t bccy(pccy);
	CCString l_ccy = bccy;

	CCString curClass = LOCAL_T4MLEG_CLASS;

	if (( l_ccy.GetLen() > 3 ) && ( !(l_ccy == "DEFAULT")) )
	   ccyIsObject = true;

	if(ARMLOCAL_TMLEG(tmIxTypeId, pstartDate, pendDate, PorRId, 
                                     pspread,
                                     payFrequencyId,
                                     resetFrequencyId,
                                     interestRuleId,
                                     fwdRuleId,
									 stubRuleId,
                                     ccyIsObject,
								     l_ccy,
                                     C_result) == ARM_OK)
	{
		Res = LocalMakeObjectId (C_result.getLong (), curClass);
	}
	
	_bstr_t tmpChaine = Res;

	*pRet = SysAllocString(tmpChaine);

	return S_OK;
}
catch(...)
{
	return E_FAIL;
}
}

STDMETHODIMP ActiveXModule::ARMGetInfoFromFxVolatility(BSTR pCurveId, VARIANT *pRetMatData)
{
	try
	{
		ARM_result C_result;
		int nbColumns, nbRows;
		double* res = NULL;
		VECTOR<double> vMatData;

		_bstr_t bCurveId(pCurveId);
		CCString l_CurveId = bCurveId;

		if (ARMLOCAL_GetInfoFromFxVolatility(LocalGetNumObjectId(l_CurveId), 
                                                        res, 
                                                        nbColumns, nbRows, 
														C_result) == ARM_OK)
		{
			 
			int rowIdx, colIdx;
			int idx = 0;
			
			for (rowIdx = 0; rowIdx < nbRows; ++rowIdx)
			{
				for (colIdx = 0; colIdx < nbColumns; ++colIdx)
				{
					if (colIdx==0)
					{
						vMatData.push_back ( res[idx]);
					}

				
					++idx;
				}
			}
			if (res)
				delete [] res;
		}

		long Res;
		Res = VECTORDOUBLE2VARIANT(vMatData, pRetMatData);
	return S_OK;
}
catch(...)
{
	return E_FAIL;
}
}

STDMETHODIMP ActiveXModule::ARMGetNewFXVolFromSummit(BSTR pCcy1, 
													 BSTR pCcy2, 
													 double pDate, 
													 BSTR pCvName, 
													 BSTR pdomZcId, 
													 BSTR pforZcId,
													 double pfxSpot,
													 VARIANT *pForwards, 
													 BSTR pwhatIsInterpolated, 
													 double pcorrectSplineWithLinear,
													 BSTR pisATM,
													 BSTR *pRet)
{
try
{
	ARM_result C_result;

	CCString Res;
	CCString curClass;

	_bstr_t bCcy1(pCcy1);
	CCString l_Ccy1 = bCcy1;

	_bstr_t bCcy2(pCcy2);
	CCString l_Ccy2 = bCcy2;

	_bstr_t bCvName(pCvName);
	CCString l_cvname = bCvName;

	_bstr_t bdomZcId(pdomZcId);
	CCString l_domZcId = bdomZcId;
	long domZcId;
	domZcId = LocalGetNumObjectId (l_domZcId);

	_bstr_t bforZcId(pforZcId);
	CCString l_forZcId = bforZcId;
	long forZcId;
	forZcId = LocalGetNumObjectId (l_forZcId);

	long forwardSize;
	VECTOR<double> vForward;
	CCString defaultForward("");

	try{
		if(VARIANT2CCString (*pForwards,defaultForward) != S_OK)
			return S_FALSE;
			
	}
	catch (...)
	{
		if(VARIANT2VECTORDOUBLE (*pForwards,vForward,forwardSize) != S_OK)
		return S_FALSE;
	}


		

	_bstr_t bwhatIsInterpolated(pwhatIsInterpolated);
	CCString l_whatIsInterpolated = bwhatIsInterpolated;
	long WhatIsInterpolatedId;
	if (( WhatIsInterpolatedId = ARM_ConvWhatIsInterp(l_whatIsInterpolated, C_result)) == ARM_DEFAULT_ERR )
	{
		return S_FALSE;
	}

	_bstr_t bisATM(pisATM);
	CCString l_isATM = bisATM;
	long isATMId;
    if (( isATMId = ARM_ConvYesOrNo(l_isATM, C_result)) == ARM_DEFAULT_ERR )
	{
	   return S_FALSE;
	}

	if(ARMLOCAL_GetNewFXVolFromSummit(l_Ccy1, l_Ccy2, pDate, l_cvname, domZcId,forZcId,pfxSpot,vForward, WhatIsInterpolatedId,pcorrectSplineWithLinear,isATMId,curClass,C_result) == ARM_OK)
	{
		Res = LocalMakeObjectId (C_result.getLong (), curClass);
	}

	_bstr_t tmpChaine = Res;

	*pRet = SysAllocString(tmpChaine);

	return S_OK;
}
catch(...)
{
	return E_FAIL;
}
}

STDMETHODIMP ActiveXModule::ARMGPModelParamCreate(BSTR pModelParamType, 
													 VARIANT *pParamTimes, 
													 VARIANT *pParamValues, 
													 BSTR pModelParamName, 
													 VARIANT *pLowerBoundaries, 
													 VARIANT *pUpperBoundaries, 
													 BSTR pInterpolMethod,
													 double pAdviseBreakPointTimes,
													 BSTR pCurrency,
													 BSTR *pRet)
{
try
{
	ARM_result C_result;
	ARM_result C_result2;

	CCString Res;
	CCString curClass;

	_bstr_t bModelParamType(pModelParamType);
	CCString l_ModelParamType = bModelParamType;

	long ModelParamTypeId;
	/*if( (ModelParamTypeId = ARM_ConvGPModelParam( l_ModelParamType, C_result)) == ARM_DEFAULT_ERR )
	{
		return S_FALSE;
	}*/

	ModelParamTypeId = ARM::ARM_ArgConv_ModelParam.GetNumber(CCSTringToSTLString(l_ModelParamType));

	

	VECTOR<double> defaultVector(0);
	defaultVector.clear();

	VECTOR<double> ParamTimes;
	VECTOR<double> ParamValues;
	VECTOR<double> LowerBoundaries;
	VECTOR<double> UpperBoundaries;
	long ParamTimesSize;
	long ParamValuesSize;
	long LowerBoundariesSize;
	long UpperBoundariesSize;
	CCString defaultParamTimes("DEFAULT");
	CCString defaultLowerBoundaries("DEFAULT");

	try{
		if(VARIANT2CCString (*pParamTimes,defaultParamTimes) != S_OK)
		{
			return S_FALSE;
		}
		else
		{
			ParamTimes.clear();
		}

	}
	catch (...)
	{
		if(VARIANT2VECTORDOUBLE (*pParamTimes,ParamTimes,ParamTimesSize) != S_OK)
			return S_FALSE;

		if(VARIANT2VECTORDOUBLE (*pParamValues,ParamValues,ParamValuesSize) != S_OK)
			return S_FALSE;
	}

	_bstr_t bModelParamName(pModelParamName);
	CCString l_ModelParamName = bModelParamName;

	try{
		if(VARIANT2CCString (*pLowerBoundaries,defaultLowerBoundaries) != S_OK)
		{
			return S_FALSE;
		}
		else
		{
			LowerBoundaries.clear();
		}

	}
	catch (...)
	{
		if(VARIANT2VECTORDOUBLE (*pLowerBoundaries,LowerBoundaries,LowerBoundariesSize) != S_OK)
			return S_FALSE;
	}

	if(LowerBoundaries.size() > 0)
    {
		if(VARIANT2VECTORDOUBLE (*pUpperBoundaries,UpperBoundaries,UpperBoundariesSize) != S_OK)
			return S_FALSE;
    }
	else 
		UpperBoundaries.clear(); 

	_bstr_t bInterpolMethod(pInterpolMethod);
	CCString l_InterpolMethod = bInterpolMethod;

	bool AdviseBreakPointTimesBool = pAdviseBreakPointTimes != 0;

	_bstr_t bCurrency(pCurrency);
	CCString l_Currency = bCurrency;

	long objId;
			
	if(ARMLOCAL_ModelParam_Create(ModelParamTypeId,ParamTimes,ParamValues,l_ModelParamName,LowerBoundaries, UpperBoundaries, l_InterpolMethod,  AdviseBreakPointTimesBool, l_Currency,C_result,ARM_NULL_OBJECT_ID) == ARM_OK)
	{
		curClass = C_result.getShortName();
		objId = C_result.getLong ();
		Res = LocalMakeObjectId (objId, curClass);	
	}

	_bstr_t tmpChaine = Res;

	*pRet = SysAllocString(tmpChaine);

	return S_OK;
}
catch(...)
{
	return E_FAIL;
}
}


STDMETHODIMP ActiveXModule::ARMFXModelCreate(BSTR pZeroCurveId, 
											 VARIANT *pParamsId,
											 double pSpot,
											 BSTR pForCurveId,
											 long pModelType,
											 BSTR *pRet)
{
try
{
	ARM_result C_result;

	_bstr_t bZeroCurveId(pZeroCurveId);
	CCString l_ZeroCurveId = bZeroCurveId;

	VECTOR<CCString> vParams;
	VECTOR<long> vParamsId;
	long ParamsSize;
		

	if(VARIANT2VECTORCCSTRING (*pParamsId,vParams,ParamsSize) != S_OK)
		return S_FALSE;

	
	for (int i = 0; i < ParamsSize; i++)
	{
		vParamsId.push_back (LocalGetNumObjectId (vParams[i]));
	}

	_bstr_t bForCurveId(pForCurveId);
	CCString l_ForCurveId = bForCurveId;

	CCString Res;
	CCString curClass = LOCAL_INFCURV_CLASS;
	long objId;

	if(ARMLOCAL_FXModel_Create(LocalGetNumObjectId(l_ZeroCurveId),
								vParamsId,
								pSpot,
								LocalGetNumObjectId(l_ForCurveId),
								"ANDREASEN",
								pModelType,
								C_result,
								ARM_NULL_OBJECT_ID) == ARM_OK)
	{
		//Res = LocalMakeObjectId (C_result.getLong (), curClass);
		curClass = C_result.getShortName();
		objId = C_result.getLong ();
		Res = LocalMakeObjectId (objId, curClass);	
	}

	_bstr_t tmpChaine = Res;

	*pRet = SysAllocString(tmpChaine);

	return S_OK;
}
catch(...)
{
	return E_FAIL;
}
}

STDMETHODIMP ActiveXModule::ARMCapFloor(BSTR pSwapLegId, 
											BSTR pIsItCapOrFloor,
											VARIANT pStrike,
											BSTR *pRet)
{
try
{
	ARM_result C_result;
	CCString Res;

	long isItCapOrFloorId;
	double dStrike;
	long   strikeType;

	_bstr_t bSwapLegId(pSwapLegId);
	CCString l_SwapLegId = bSwapLegId;

	_bstr_t bIsItCapOrFloor(pIsItCapOrFloor);
	CCString l_IsItCapOrFloor = bIsItCapOrFloor;

	if((isItCapOrFloorId = ARM_ConvCapOrFloor (l_IsItCapOrFloor, C_result)) == ARM_DEFAULT_ERR)
	{
		return S_FALSE;
	}


	if (VARIANT2Double(pStrike,&dStrike) == S_FALSE)
	{
		_bstr_t bStrike(pStrike);
		CCString l_Strike = bStrike;

		dStrike = (double) LocalGetNumObjectId(l_Strike);
		strikeType = 1L;
	}
	else
	{
	   strikeType = 0L;
	}

	CCString curClass = CorrespondingCapFloorClasswId( LocalGetNumObjectId (l_SwapLegId) );

	if(ARMLOCAL_CAPFLOOR (LocalGetNumObjectId (l_SwapLegId),
									 isItCapOrFloorId,
									 strikeType,
									 dStrike,
									 C_result)== ARM_OK)
	{
		Res = LocalMakeObjectId (C_result.getLong (), curClass);
	}

	_bstr_t tmpChaine = Res;

	*pRet = SysAllocString(tmpChaine);

	return S_OK;
}
catch(...)
{
	return E_FAIL;
}
}


STDMETHODIMP ActiveXModule::ARMSpreadDigital(double pStartDate,
											 double pEndDate,
											 BSTR pCapOrFloor,
											 VARIANT pStrike,
											 VARIANT *pSpread,
											 VARIANT pPayOff,
											 BSTR pLiborType1,
											 BSTR pLiborType2,
											 VARIANT *pWeight,
											 BSTR pDayCount,
											 BSTR pResetFreq,
											 BSTR pPayFreq,
											 BSTR pResetTiming,
											 BSTR pPayTiming,
											 BSTR pCurrency,
											 double pResetGap,
											 BSTR pIntRule,
											 BSTR pStubRule,
											 VARIANT *pFixing1,
											 VARIANT *pFixing2,
											 BSTR *pRet)
{
try
{
	ARM_result C_result;
	CCString Res;

	long capOrFloorId;
	double dStrike;
	double dPayOff;
	long   strike_type;
	long   payoff_type;
	long liborType1Id;
	long liborType2Id;
	double dWeight1;
	double dWeight2;
	double dSpread1;
	double dSpread2;
	long l_slopeFlag;
	long l_slopeFlagDefault=1;
	double dCptStrikeMethod;
	double dCptStrikeMethodDefault = 1.0;
	double dComputedFormula;
	double dComputedFormulaDefault = 1.0;
	long dayCountId;
	long resetFreqId;
	long payFreqId;
	long resetTimingId;
	long payTimingId;
	long l_ccyId;
	long fixing1_type;
	long fixing2_type;
	double dFixing1;
	double dFixing2;
	long intRuleId;
	long stubRuleId;

	CCString defaultFixing("DEFAULT");
	
	long WeightSize;
	long SpreadSize;
	long Fixing1Size;
	long Fixing2Size;

	VECTOR<CCString> vWeight;
	VECTOR<CCString> vWeightDefault;
	VECTOR<double> vCalibInfo;
	VECTOR<double> vSpread;
	VECTOR<double> vFixing1;
	VECTOR<double> vFixing2;
	

	_bstr_t bCapOrFloor(pCapOrFloor);
	CCString l_CapOrFloor = bCapOrFloor;
	if((capOrFloorId = ARM_ConvCapOrFloor (l_CapOrFloor, C_result)) == ARM_DEFAULT_ERR)
	{
		return S_FALSE;
	}

	if (VARIANT2Double(pStrike,&dStrike) == S_FALSE)
	{
		_bstr_t bStrike(pStrike);
		CCString l_Strike = bStrike;

		dStrike = (double) LocalGetNumObjectId(l_Strike);
		strike_type = 1L;
	}
	else
	{
		strike_type = 0L;
	}

	_bstr_t bLiborType1(pLiborType1);
	CCString l_LiborType1 = bLiborType1;
	liborType1Id = ARM_ConvIrType (l_LiborType1);

	_bstr_t bLiborType2(pLiborType2);
	CCString l_LiborType2 = bLiborType2;
	liborType2Id = ARM_ConvIrType (l_LiborType2);

	
	if (VARIANT2Double(pPayOff,&dPayOff) == S_FALSE)
	{
		_bstr_t bPayOff(pPayOff);
		CCString l_PayOff = bPayOff;	
		dPayOff = (double) LocalGetNumObjectId(l_PayOff);
		payoff_type = 1L;
	}
	else
	{
		payoff_type = 0L;
	}

	CCString l_Weight;
	
	try
	{
		if(VARIANT2CCString (*pWeight,l_Weight) != S_OK)
		{
			return S_FALSE;
		}
		else
		{
			if (l_Weight=="DEFAULT")
			{
				vWeight=vWeightDefault;
			}
			else
			{
				return S_FALSE;
			}				
		}
	}
	catch (...)
	{
		if(VARIANT2VECTORCCSTRING (*pWeight,vWeight,WeightSize) != S_OK)
			return S_FALSE;
	}

	

	if(vWeight.size()<2)
		return S_FALSE;

	sscanf ((const char*)vWeight[0], "%lf", &dWeight1);
	sscanf ((const char*)vWeight[1], "%lf", &dWeight2);

	if(vWeight.size() == 2)
	{
		l_slopeFlag = l_slopeFlagDefault;
		dCptStrikeMethod = dCptStrikeMethodDefault;
		dComputedFormula = dComputedFormulaDefault;
	}
	else if(vWeight.size() == 3)
	{
		sscanf ((const char*)vWeight[2], "%ld", &l_slopeFlag);
		dCptStrikeMethod = dCptStrikeMethodDefault;
		dComputedFormula = dComputedFormulaDefault;
	}
	else if(vWeight.size() == 4)
	{
		sscanf ((const char*)vWeight[2], "%ld", &l_slopeFlag);
		sscanf ((const char*)vWeight[3], "%lf", &dCptStrikeMethod);
		dComputedFormula = dComputedFormulaDefault;
	}
	else
	{
		sscanf ((const char*)vWeight[2], "%ld", &l_slopeFlag);
		sscanf ((const char*)vWeight[3], "%lf", &dCptStrikeMethod);
		sscanf ((const char*)vWeight[4], "%lf", &dComputedFormula);

		// get Sabr calib infos for each index

		double tmp;
		for (int j=0;j<vWeight.size()-5;j++)
		{
			sscanf ((const char*)vWeight[j+5], "%lf", &tmp);
			vCalibInfo.push_back (tmp);
		}
	}

	_bstr_t bDayCount(pDayCount);
	CCString l_DayCount = bDayCount;
	dayCountId = ARM_ConvDayCount (l_DayCount);


	_bstr_t bResetFreq(pResetFreq);
	CCString l_ResetFreq = bResetFreq;
	if((resetFreqId = ARM_ConvFrequency (l_ResetFreq, C_result)) == ARM_DEFAULT_ERR)
	{
		return S_FALSE;
	}

	_bstr_t bPayFreq(pPayFreq);
	CCString l_PayFreq = bPayFreq;
	if((payFreqId = ARM_ConvFrequency (l_PayFreq, C_result)) == ARM_DEFAULT_ERR)
	{
		return S_FALSE;
	}

	_bstr_t bResetTiming(pResetTiming);
	CCString l_ResetTiming = bResetTiming;
	resetTimingId = ARM_ConvPayResetRule (l_ResetTiming);

	_bstr_t bPayTiming(pPayTiming);
	CCString l_PayTiming = bPayTiming;
	payTimingId = ARM_ConvPayResetRule (l_PayTiming);

	_bstr_t bCcy(pCurrency);
	CCString l_Ccy = bCcy;

	if (l_Ccy == "DEFAULT")
		l_ccyId = ARM_NULL_OBJECT;
	else
		l_ccyId = LocalGetNumObjectId (l_Ccy);

	if(VARIANT2VECTORDOUBLE (*pSpread,vSpread,SpreadSize) != S_OK)
		return S_FALSE;

	if(vSpread.size()<2)
		return S_FALSE;

	dSpread1 = vSpread[0];
	dSpread2 = vSpread[1];

	CCString l_Fixing1;
	CCString l_Fixing2;

	try
	{
		if(VARIANT2CCString (*pFixing1,l_Fixing1) != S_OK)
		{
			return S_FALSE;
		}
		else
		{
			if (l_Fixing1=="DEFAULT")
			{
				fixing1_type=1;
				dFixing1=ARM_NULL_OBJECT;
				vFixing1.clear();
			}
			else
			{
				fixing1_type=1;
				dFixing1=(double) LocalGetNumObjectId(l_Fixing1);
			}				
		}

	}
	catch (...)
	{
		if (VARIANT2VECTORDOUBLE(*pFixing1,vFixing1,Fixing1Size) != S_OK)
		{
			return S_FALSE;
		}
		else
		{
			fixing1_type=0;
		}
	}

	try
	{
		if(VARIANT2CCString (*pFixing2,l_Fixing2) != S_OK)
		{
			return S_FALSE;
		}
		else
		{
			if (l_Fixing2=="DEFAULT")
			{
				fixing2_type=1;
				dFixing2=ARM_NULL_OBJECT;
				vFixing2.clear();
			}
			else
			{
				fixing2_type=1;
				dFixing2=(double) LocalGetNumObjectId(l_Fixing2);
			}				
		}

	}
	catch (...)
	{
		if (VARIANT2VECTORDOUBLE(*pFixing2,vFixing2,Fixing1Size) != S_OK)
		{
			return S_FALSE;
		}
		else
		{
			fixing2_type=0;
		}
	}

	_bstr_t bIntRule(pIntRule);
	CCString l_IntRule = bIntRule;
	intRuleId = ARM_ConvIntRule(l_IntRule);

	_bstr_t bStubRule(pStubRule);
	CCString l_StubRule = bStubRule;
	stubRuleId = ARM_ConvStubRule(l_StubRule);


	CCString curClass = LOCAL_SPREAD_OPTION_CLASS;

	if(ARMLOCAL_SPREADDIGITAL (pStartDate,
								 pEndDate,
								 capOrFloorId,
								 strike_type,
								 dStrike,
								payoff_type,
								 dPayOff,
								 liborType1Id,
								 liborType2Id,
								 dWeight1,
								 dWeight2,
								 dayCountId,
								 resetFreqId,
								 payFreqId,
								 resetTimingId,
								 payTimingId,
								 l_ccyId,
                                 long(pResetGap),
								 dSpread1,
								 dSpread2,
								 (long) fixing1_type,
								 dFixing1,
								 vFixing1,
								 (long) fixing2_type,
								 dFixing2,
								 vFixing2,
								 intRuleId,
								 stubRuleId,
								 l_slopeFlag,
								 (int) dCptStrikeMethod,
								 (int) dComputedFormula,
								 vCalibInfo,
								 C_result)== ARM_OK)
	{
		Res = LocalMakeObjectId (C_result.getLong (), curClass);
	}

	_bstr_t tmpChaine = Res;
	*pRet = SysAllocString(tmpChaine);
	return S_OK;
}
catch(...)
{
	return E_FAIL;
}
}


STDMETHODIMP ActiveXModule::ARMSpreadDigitalFlt(double pStartDate,
											 double pEndDate,
											 BSTR pCapOrFloor,
											 VARIANT pStrike,
											 VARIANT *pSpread,
											 BSTR pPayOffLiborType,
											 BSTR pLiborType1,
											 BSTR pLiborType2,
											 VARIANT *pWeight,
											 BSTR pDayCount,
											 BSTR pResetFreq,
											 BSTR pPayFreq,
											 BSTR pResetTiming,
											 BSTR pPayTiming,
											 BSTR pCurrency,
											 double pResetGap,
											 BSTR pIntRule,
											 BSTR pStubRule,
											 VARIANT *pFixing1,
											 VARIANT *pFixing2,
											 BSTR *pRet)
{
try
{
	ARM_result C_result;
	CCString Res;

	long capOrFloorId;
	double dStrike;
	double dPayOff;
	long   strike_type;
	long liborType1Id;
	long liborType2Id;
	double dWeight1;
	double dWeight2;
	double dSpread1;
	double dSpread2;
	long l_slopeFlag;
	long l_slopeFlagDefault=1;
	double dCptStrikeMethod;
	double dCptStrikeMethodDefault = 1.0;
	double dComputedFormula;
	double dComputedFormulaDefault = 1.0;
	long dayCountId;
	long resetFreqId;
	long payFreqId;
	long resetTimingId;
	long payTimingId;
	long l_ccyId;
	long fixing1_type;
	long fixing2_type;
	double dFixing1;
	double dFixing2;
	long intRuleId;
	long stubRuleId;

	CCString defaultFixing("DEFAULT");
	
	long WeightSize;
	long SpreadSize;
	long Fixing1Size;
	long Fixing2Size;
	long payoffLiborTypeId;

	VECTOR<double> vWeight;
	VECTOR<double> vWeightDefault(1.0);

	VECTOR<double> vSpread;
	VECTOR<double> vFixing1;
	VECTOR<double> vFixing2;

	CCString l_Weight;
	

	_bstr_t bCapOrFloor(pCapOrFloor);
	CCString l_CapOrFloor = bCapOrFloor;
	if((capOrFloorId = ARM_ConvCapOrFloor (l_CapOrFloor, C_result)) == ARM_DEFAULT_ERR)
	{
		return S_FALSE;
	}

	if (VARIANT2Double(pStrike,&dStrike) == S_FALSE)
	{
		_bstr_t bStrike(pStrike);
		CCString l_Strike = bStrike;

		dStrike = (double) LocalGetNumObjectId(l_Strike);
		strike_type = 1L;
	}
	else
	{
		strike_type = 0L;
	}

	_bstr_t bPayOffLiborType(pPayOffLiborType);
	CCString l_PayOffLiborType = bPayOffLiborType;
	payoffLiborTypeId = ARM_ConvIrType (l_PayOffLiborType);

	_bstr_t bLiborType1(pLiborType1);
	CCString l_LiborType1 = bLiborType1;
	liborType1Id = ARM_ConvIrType (l_LiborType1);

	_bstr_t bLiborType2(pLiborType2);
	CCString l_LiborType2 = bLiborType2;
	liborType2Id = ARM_ConvIrType (l_LiborType2);


	
	/*if (VARIANT2Double(pPayOff,&dPayOff) == S_FALSE)
	{
		_bstr_t bPayOff(pPayOff);
		CCString l_PayOff = bPayOff;	
		dPayOff = (double) LocalGetNumObjectId(l_PayOff);
		payoff_type = 1L;
	}
	else
	{
		payoff_type = 0L;
	}*/

	try
	{
		if(VARIANT2CCString (*pWeight,l_Weight) != S_OK)
		{
			return S_FALSE;
		}
		else
		{
			if (l_Weight=="DEFAULT")
			{
				vWeight=vWeightDefault;
			}
			else
			{
				return S_FALSE;
			}				
		}
	}
	catch (...)
	{
		if (VARIANT2VECTORDOUBLE(*pWeight,vWeight,WeightSize) != S_OK)
		{
			return S_FALSE;
		}
	}
	
	if(vWeight.size()<2)
		return S_FALSE;

	dWeight1=vWeight[0];
	dWeight2=vWeight[1];

	if(vWeight.size() == 2)
	{
		l_slopeFlag = l_slopeFlagDefault;
		dCptStrikeMethod = dCptStrikeMethodDefault;
		dComputedFormula = dComputedFormulaDefault;
	}
	else if(vWeight.size() == 3)
	{
		l_slopeFlag=vWeight[2];
		dCptStrikeMethod = dCptStrikeMethodDefault;
		dComputedFormula = dComputedFormulaDefault;
	}
	else if(vWeight.size() == 4)
	{
		l_slopeFlag=vWeight[2];
		dCptStrikeMethod=vWeight[3];
		dComputedFormula = dComputedFormulaDefault;
	}
	else
	{
		l_slopeFlag=vWeight[2];
		dCptStrikeMethod=vWeight[3];
		dComputedFormula=vWeight[4];
	}

	_bstr_t bDayCount(pDayCount);
	CCString l_DayCount = bDayCount;
	dayCountId = ARM_ConvDayCount (l_DayCount);


	_bstr_t bResetFreq(pResetFreq);
	CCString l_ResetFreq = bResetFreq;
	if((resetFreqId = ARM_ConvFrequency (l_ResetFreq, C_result)) == ARM_DEFAULT_ERR)
	{
		return S_FALSE;
	}

	_bstr_t bPayFreq(pPayFreq);
	CCString l_PayFreq = bPayFreq;
	if((payFreqId = ARM_ConvFrequency (l_PayFreq, C_result)) == ARM_DEFAULT_ERR)
	{
		return S_FALSE;
	}

	_bstr_t bResetTiming(pResetTiming);
	CCString l_ResetTiming = bResetTiming;
	resetTimingId = ARM_ConvPayResetRule (l_ResetTiming);

	_bstr_t bPayTiming(pPayTiming);
	CCString l_PayTiming = bPayTiming;
	payTimingId = ARM_ConvPayResetRule (l_PayTiming);

	_bstr_t bCcy(pCurrency);
	CCString l_Ccy = bCcy;

	if (l_Ccy == "DEFAULT")
		l_ccyId = ARM_NULL_OBJECT;
	else
		l_ccyId = LocalGetNumObjectId (l_Ccy);

	if(VARIANT2VECTORDOUBLE (*pSpread,vSpread,SpreadSize) != S_OK)
		return S_FALSE;

	if(vSpread.size()<2)
		return S_FALSE;

	dSpread1 = vSpread[0];
	dSpread2 = vSpread[1];

	CCString l_Fixing1;
	CCString l_Fixing2;

	try
	{
		if(VARIANT2CCString (*pFixing1,l_Fixing1) != S_OK)
		{
			return S_FALSE;
		}
		else
		{
			if (l_Fixing1=="DEFAULT")
			{
				fixing1_type=1;
				dFixing1=ARM_NULL_OBJECT;
				vFixing1.clear();
			}
			else
			{
				fixing1_type=1;
				dFixing1=(double) LocalGetNumObjectId(l_Fixing1);
			}				
		}

	}
	catch (...)
	{
		if (VARIANT2VECTORDOUBLE(*pFixing1,vFixing1,Fixing1Size) != S_OK)
		{
			return S_FALSE;
		}
		else
		{
			fixing1_type=0;
		}
	}

	try
	{
		if(VARIANT2CCString (*pFixing2,l_Fixing2) != S_OK)
		{
			return S_FALSE;
		}
		else
		{
			if (l_Fixing2=="DEFAULT")
			{
				fixing2_type=1;
				dFixing2=ARM_NULL_OBJECT;
				vFixing2.clear();
			}
			else
			{
				fixing2_type=1;
				dFixing2=(double) LocalGetNumObjectId(l_Fixing2);
			}				
		}

	}
	catch (...)
	{
		if (VARIANT2VECTORDOUBLE(*pFixing2,vFixing2,Fixing1Size) != S_OK)
		{
			return S_FALSE;
		}
		else
		{
			fixing2_type=0;
		}
	}

	_bstr_t bIntRule(pIntRule);
	CCString l_IntRule = bIntRule;
	intRuleId = ARM_ConvIntRule(l_IntRule);

	_bstr_t bStubRule(pStubRule);
	CCString l_StubRule = bStubRule;
	stubRuleId = ARM_ConvStubRule(l_StubRule);


	CCString curClass = LOCAL_SPREAD_OPTION_CLASS;

	if(ARMLOCAL_SPREADDIGITALFLT (pStartDate,
								 pEndDate,
								 capOrFloorId,
								 strike_type,
								 dStrike,
								 payoffLiborTypeId,
								 liborType1Id,
								 liborType2Id,
								 dWeight1,
								 dWeight2,
								 dayCountId,
								 resetFreqId,
								 payFreqId,
								 resetTimingId,
								 payTimingId,
								 l_ccyId,
                                 long(pResetGap),
								 dSpread1,
								 dSpread2,
								 (long) fixing1_type,
								 dFixing1,
								 vFixing1,
								 (long) fixing2_type,
								 dFixing2,
								 vFixing2,
								 intRuleId,
								 stubRuleId,
								 l_slopeFlag,
								 (int) dCptStrikeMethod,
								 (int) dComputedFormula,
								 C_result)== ARM_OK)
	{
		Res = LocalMakeObjectId (C_result.getLong (), curClass);
	}

	_bstr_t tmpChaine = Res;
	*pRet = SysAllocString(tmpChaine);
	return S_OK;
}
catch(...)
{
	return E_FAIL;
}
}







