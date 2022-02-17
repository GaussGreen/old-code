#include "firsttobeincluded.h"
#include "ActiveXModule.h"
#include <libCCatl\CCatl.h>


// #include <ARM\libarm_local\ARM_local_class.h>

#include <ARM\libarm\ARM_result.h>
#include <ARM\libarm_local\ARM_local_glob.h>
#include <ARM\libarm_local\ARM_local_gp_inflation.h>
#include <ARM\libarm_local\ARM_local_swap.h>
#include <ARM\libarm_local\ARM_local_ccy.h>
#include <ARM\libarm_local\ARM_local_mod.h>
#include <ARM\libarm_local\ARM_local_zccurve.h>
#include <ARM\libarm_local\ARM_local_capfl.h>
#include <ARM\libarm_local\ARM_local_refval.h>
#include <ARM\libarm_local\ARM_local_volcrv.h>
#include <ARM\libarm_local\ARM_local_irindex.h>
#include <ARM\libarm_frometk\ARM_local_etoolkit.h>
#include <ARM\libarm_local\ARM_local_pf.h>
#include <ARM\libarm_local\ARM_local_option.h>
#include <ARM\libarm_local\ARM_local_forex.h>

#include "ARM_local_interglob.h"
#include <util\fromto.h>


STDMETHODIMP ActiveXModule::ARMSetEtoolkit(BSTR pUserName,
									   BSTR pPassWord,
									   BSTR pDatabaseContext,
									   BSTR pItConfigDomainDir,
									   BSTR pItDomainName,
									   long *pRet)
{
	// TODO: Add your implementation code here
try
{
	_bstr_t bUserName(pUserName);
	CCString l_username = bUserName;
	_bstr_t bPassWord(pPassWord);
	CCString l_password = bPassWord;
	_bstr_t bDatabaseContext(pDatabaseContext);
	CCString l_databasecontext = bDatabaseContext;
	_bstr_t bItConfigDomainDir(pItConfigDomainDir);
	CCString l_itconfigdomaindir = bItConfigDomainDir;
	_bstr_t bItDomainName(pItDomainName);
	CCString l_itdomainname = bItDomainName;

	set_etoolkit(l_username, l_password, l_databasecontext,l_itconfigdomaindir, l_itdomainname);

	*pRet=1;

	return S_OK;
}
catch(...)
{
	return E_FAIL;
}

}


STDMETHODIMP ActiveXModule::ARMConnectionEtoolkit(VARIANT *pRet)
{
	// TODO: Add your implementation code here
try
{
	connection_etoolkit();

	_variant_t wrap_Res;
	wrap_Res.Attach(*pRet);
	wrap_Res.intVal = 1;
	*pRet=wrap_Res.Detach();

	pRet->vt=VT_I4;

	return S_OK;
}
catch(...)
{
	return E_FAIL;
}

}

STDMETHODIMP ActiveXModule::ARMDeconnectionEtoolkit(VARIANT *pRet)
{
	// TODO: Add your implementation code here
try
{
	deconnection_etoolkit();

	_variant_t wrap_Res;
	wrap_Res.Attach(*pRet);
	wrap_Res.intVal = 1;
	*pRet=wrap_Res.Detach();

	pRet->vt=VT_I4;

	return S_OK;
}
catch(...)
{
	return E_FAIL;
}

}


double compute_le_i(double asOf, double startDate, double endDate, long LiborTypeId, const CCString& ccyStr, long amortMethodId, long resetFreqId, CCString bsmodId)
{

	ARM_result C_result;

	double le_i = 0.0;
	long floatId, fixId, swapId, genAmortId;
	long floatAmortId, fixAmortId;

	if (ARMLOCAL_LIBORLEG(startDate,
						  endDate,
						  LiborTypeId,
						  K_PAY,
						  0,
						  0.0,
						  K_DEF_FREQ,
						  K_DEF_FREQ,
						  K_ADVANCE,
						  K_ARREARS,
						  false,
						  ccyStr,
						  K_ADJUSTED,
						  10000.,
						  "NULL",
						  "NULL",
						  1,
						  K_NX_NONE,
						  K_SHORTSTART,
						  -1.0,
						  K_YES,
						  -1,
						  C_result) == ARM_OK) 
	{
		floatId = C_result.getLong();
	}

	if (ARMLOCAL_ARM_GenAmortization(floatId,
									 amortMethodId,
									 resetFreqId,
									 0.0,
									 -1,
									 -9999.,
									 0.,
									 0.,
									 LocalGetNumObjectId(bsmodId),
									 0.0,
									 C_result) == ARM_OK) 
	{
		genAmortId = C_result.getLong();
	}

	if (ARMLOCAL_ClonedAndSetNotional(floatId,
									  genAmortId,
									  100.,
									  C_result) == ARM_OK)
	{
		floatAmortId = C_result.getLong();
	}

	if (ARMLOCAL_FIXEDLEG(startDate,
						  endDate,
						  K_RCV,
						  0,
						  1.0,
						  K30_360,
						  K_ANNUAL,
						  K_COMP_PROP,
						  K_ARREARS,
						  K_ADJUSTED,
						  K_SHORTSTART,
						  false,
						  ccyStr,
						  "NULL",
						  K_NX_NONE,
						  -1.0,
						  K_YES,
						  -1,
						  C_result) == ARM_OK) 
	{
		fixId = C_result.getLong();
	}

	if (ARMLOCAL_ClonedAndSetNotional(fixId,
									  genAmortId,
									  100.,
									  C_result) == ARM_OK)
	{
		fixAmortId = C_result.getLong();
	}

	vector<double> fixedRates;

	if (ARMLOCAL_SWAP(floatAmortId,
					  fixAmortId,
					  -1.0,
					  fixedRates,
					  C_result) == ARM_OK)
	{
		swapId = C_result.getLong();
	}

	if (ARMLOCAL_SWAP_PRICE_TO_RATE(swapId,
									asOf,
									0.0,
									LocalGetNumObjectId(bsmodId),
									C_result) == ARM_OK)
	{
		le_i = C_result.getDouble();
	}

	long retCode = ARMLOCAL_FreeObject(floatId,C_result);
	retCode = ARMLOCAL_FreeObject(floatAmortId,C_result);
	retCode = ARMLOCAL_FreeObject(fixId,C_result);
	retCode = ARMLOCAL_FreeObject(fixAmortId,C_result);
	retCode = ARMLOCAL_FreeObject(swapId,C_result);
	retCode = ARMLOCAL_FreeObject(genAmortId,C_result);

	return le_i;
}




STDMETHODIMP ActiveXModule::ARMcptBonibor(double pDate,
									  BSTR pDepart,
									  BSTR pMatStruct,
									  BSTR pMatTot,
									  BSTR pAmort,
									  BSTR pFreq,
									  BSTR pSjUSD,
									  BSTR pTiming,
									  double pBarriere,
									  double pSpdPostBar,
									  double pMarge,
									  double pFunding,
									  BSTR pFundingFreq,
									  double pSpd2phase,
									  double pSoulte,
									  BSTR pYcModId,
									  BSTR pBsModId,
									  BSTR pBsModVolUSDId,
									  BSTR pBsModCorrPlusId,
									  BSTR pBsModCorrMoinsId,
									  BSTR pCrossModId,
									  BSTR pProbaMarge,
									  double pInt,
									  VARIANT *pRet)
{
	// TODO: Add your implementation code here
try
{

	ARM_result C_result;

	_bstr_t bdepart(pDepart);
	CCString l_depart = bdepart;

	_bstr_t bmatstruct(pMatStruct);
	CCString l_matstruct = bmatstruct;

	_bstr_t bmattot(pMatTot);
	CCString l_mattot = bmattot;

	_bstr_t bamort(pAmort);
	CCString l_amort = bamort;

	_bstr_t bfreq(pFreq);
	CCString l_freq = bfreq;

	_bstr_t bfundfreq(pFundingFreq);
	CCString l_fundfreq = bfundfreq;

	_bstr_t bresettiming(pTiming);
	CCString l_resettiming = bresettiming;

	_bstr_t bSjUSD(pSjUSD);
	CCString l_SjUSD = bSjUSD;

	_bstr_t bYcModId(pYcModId);
	CCString l_YcModId = bYcModId;

	_bstr_t bBsModId(pBsModId);
	CCString l_BsModId = bBsModId;

	_bstr_t bBsModVolUSDId(pBsModVolUSDId);
	CCString l_BsModVolUSDId = bBsModVolUSDId;

	_bstr_t bBsModCorrPlusId(pBsModCorrPlusId);
	CCString l_BsModCorrPlusId = bBsModCorrPlusId;

	_bstr_t bBsModCorrMoinsId(pBsModCorrMoinsId);
	CCString l_BsModCorrMoinsId = bBsModCorrMoinsId;

	_bstr_t bCrossModId(pCrossModId);
	CCString l_CrossModId = bCrossModId;

	_bstr_t bProbaMarge(pProbaMarge);
	CCString l_ProbaMarge = bProbaMarge;

	double Res = 0.;
	double Res2 = 0.;
	long retCode;
	double the_i;

	// startDate
	double startDate;
	int Nb;
	long freqId;

	GetFrqAndNbFromString((const char*)l_depart,Nb,freqId);

	if(ARMLOCAL_ARM_ADDPERIOD(pDate,K_DAILY,"EUR",2,0,0,C_result) == ARM_OK)
	{
		startDate = Local_ARMDATE2XLDATE(C_result.getString ());
	}
	else
	{
		ERROR_MSG(C_result.getMsg(),pRet);
		return S_OK;
	}

	if(ARMLOCAL_ARM_ADDPERIOD(startDate,freqId,"EUR",Nb,0,0,C_result) == ARM_OK)
	{
		startDate = Local_ARMDATE2XLDATE(C_result.getString ());
	}
	else
	{
		ERROR_MSG(C_result.getMsg(),pRet);
		return S_OK;
	}
	
	// MatStruct
	double endMaturityStruct;

	GetFrqAndNbFromString(l_matstruct,Nb,freqId);

	if(ARMLOCAL_ARM_ADDPERIOD(startDate,freqId,"EUR",Nb,0,0,C_result) == ARM_OK)
	{
		endMaturityStruct = Local_ARMDATE2XLDATE(C_result.getString ());
	}
	else
	{
		ERROR_MSG(C_result.getMsg(),pRet);
		return S_OK;
	}

	// matTot
	double endMaturityTot;

	GetFrqAndNbFromString(l_mattot,Nb,freqId);

	if(ARMLOCAL_ARM_ADDPERIOD(startDate,freqId,"EUR",Nb,0,0,C_result) == ARM_OK)
	{
		endMaturityTot = Local_ARMDATE2XLDATE(C_result.getString ());
	}
	else
	{
		ERROR_MSG(C_result.getMsg(),pRet);
		return S_OK;
	}

	// Amortissement
	double le_i;

	long liborTypeId;
	long sjacLiborTypeId;
	long amortMethodId;
	long amortId;
	long liborlegId;
	long clonedLiborlegid;
	long resetFreqId;
	long resetFundFreqId;
	long fixedlegId;
	long sjUSDTypeId;

	if (strcmp(l_freq,"Q") == 0)
		liborTypeId = K_EURIBOR3M;
	else
		liborTypeId = K_EURIBOR1Y;

	if (strcmp(l_fundfreq,"Q") == 0)
		sjacLiborTypeId = K_EURIBOR3M;
	else
		sjacLiborTypeId = K_EURIBOR1Y;

	long timingId = ARM_ConvPayResetRule (l_resettiming);

	double gap;
	if (timingId == K_ADVANCE)
		gap = -2;
	else
		gap = -15;


	if((resetFreqId = ARM_ConvFrequency(l_freq, C_result)) == ARM_DEFAULT_ERR)
	{
		ERROR_MSG("Invalid Frq",pRet,ARM_ERROR_FREQ);
		return S_OK;
	}

	if((resetFundFreqId = ARM_ConvFrequency(l_fundfreq, C_result)) == ARM_DEFAULT_ERR)
	{
		ERROR_MSG("Invalid Frq",pRet,ARM_ERROR_FREQ);
		return S_OK;
	}

	if((sjUSDTypeId = ARM_ConvIrIndName(l_SjUSD, C_result)) == ARM_DEFAULT_ERR)
	{
		ERROR_MSG("Invalid Libor Type",pRet,ARM_ERROR_FREQ);
		return S_OK;
	}

	//*************************************************************
	//**********    Create Armort Refvalue ID   **********
	//*************************************************************
	long amortIdToDestroy = ARM_NULL_OBJECT;

	if (LocalGetNumObjectId(l_amort) == ARM_KO)
	{
		if((amortMethodId = ARM_ConvAmortMethod (l_amort, C_result)) == ARM_DEFAULT_ERR)
		{
			ERROR_MSG("Invalid AmortMethod",pRet,ARM_ERROR_FREQ);
			return S_OK;
		}

		le_i = compute_le_i(pDate,startDate,endMaturityTot,liborTypeId,"EUR",amortMethodId,resetFreqId,l_BsModId);

		if (pInt == -1.0)
			the_i = le_i;
		else
			the_i = pInt;

		if (ARMLOCAL_FIXEDLEG(startDate,
							  endMaturityTot,
							  K_PAY,
							  0,
							  the_i,
							  K30_360,
							  K_ANNUAL,
							  K_COMP_PROP,
							  K_ARREARS,
							  K_ADJUSTED,
							  K_SHORTSTART,
							  false,
							  "EUR",
							  "NULL",
							  K_NX_NONE,
							  -1.0,
							  K_YES,
						  -1,
							  C_result) == ARM_OK) 
		{
			fixedlegId = C_result.getLong();
		}
		else
		{
			ERROR_MSG(C_result.getMsg(),pRet);
			return S_OK;
		}

		if (ARMLOCAL_ARM_GenAmortization(fixedlegId,
										 amortMethodId,
										 resetFreqId,
										 0.0,
										 -1,
										 -9999.,
										 0.,
										 0.,
										 LocalGetNumObjectId(l_BsModId),
										 0.0,
										 C_result) == ARM_OK) 
		{
			amortId = C_result.getLong();
			amortIdToDestroy = amortId;
		}
		else
		{
			ERROR_MSG(C_result.getMsg(),pRet);
			return S_OK;
		}

		retCode = ARMLOCAL_FreeObject(fixedlegId,C_result);
	}
	else
	{
		amortId = LocalGetNumObjectId(l_amort);
	}
	// FUNDING

	double fund;
	double funding2Price = 0.0;

	if (ARMLOCAL_LIBORLEG(startDate,
						  endMaturityStruct,
						  sjacLiborTypeId,
						  K_PAY,
						  0,
						  (pFunding + pMarge) / 100.0,
						  resetFundFreqId,
						  resetFundFreqId,
						  K_ADVANCE,
						  K_ARREARS,
						  false,
						  "EUR",
						  K_ADJUSTED,
						  -2.0,
						  "NULL",
						  "NULL",
						  1,
						  K_NX_NONE,
						  K_SHORTSTART,
						  -1.0,
						  K_YES,
						  -1,
						  C_result) == ARM_OK) 
	{
		liborlegId = C_result.getLong();
	}
	else
	{
		ERROR_MSG(C_result.getMsg(),pRet);
		return S_OK;
	}

	if (ARMLOCAL_ClonedAndSetNotional(liborlegId,
									  amortId,
									  100.,
									  C_result) == ARM_OK)
	{
		clonedLiborlegid = C_result.getLong();
	}
	else
	{
		ERROR_MSG(C_result.getMsg(),pRet);
		return S_OK;
	}

	if (ARMLOCAL_ARM_Price(clonedLiborlegid,
						   LocalGetNumObjectId(l_YcModId),
						   C_result) == ARM_OK)
	{
		fund = C_result.getDouble();
	}
	else
	{
		ERROR_MSG(C_result.getMsg(),pRet);
		return S_OK;
	}

	retCode = ARMLOCAL_FreeObject(liborlegId,C_result);
	retCode = ARMLOCAL_FreeObject(clonedLiborlegid,C_result);

	if (endMaturityTot > endMaturityStruct)
	{
		if (ARMLOCAL_LIBORLEG(endMaturityStruct,
							  endMaturityTot,
							  sjacLiborTypeId,
							  K_PAY,
							  0,
							  pFunding / 100.0,
							  resetFundFreqId,
							  resetFundFreqId,
							  K_ADVANCE,
							  K_ARREARS,
							  false,
							  "EUR",
							  K_ADJUSTED,
							  -2.0,
							  "NULL",
							  "NULL",
							  1,
							  K_NX_NONE,
							  K_SHORTSTART,
							  -1.0,
							  K_YES,
							  -1,
							  C_result) == ARM_OK) 
		{
			liborlegId = C_result.getLong();
		}
		else
		{
			ERROR_MSG(C_result.getMsg(),pRet);
			return S_OK;
		}

		if (ARMLOCAL_ClonedAndSetNotional(liborlegId,
										  amortId,
										  100.,
										  C_result) == ARM_OK)
		{
			clonedLiborlegid = C_result.getLong();
		}
		else
		{
			ERROR_MSG(C_result.getMsg(),pRet);
			return S_OK;
		}

		if (ARMLOCAL_ARM_Price(clonedLiborlegid,
							   LocalGetNumObjectId(l_YcModId),
							   C_result) == ARM_OK)
		{
			funding2Price = C_result.getDouble();
		}
		else
		{
			ERROR_MSG(C_result.getMsg(),pRet);
			return S_OK;
		}

		retCode = ARMLOCAL_FreeObject(liborlegId,C_result);
		retCode = ARMLOCAL_FreeObject(clonedLiborlegid,C_result);
	}

	// STRUCTURE
	long ccyUSDId;
	long ccyEURId;
	long irIndexId;
	
	long capBId;
	double capB;
	double capB_volUSD;
	double capB_corrplus;
	double capB_corrmoins;
	
	long digitId;
	double digit;
	double digit_volUSD;
	double digit_corrplus;
	double digit_corrmoins;

	long dualcapId;
	long dualclonedId;
	double dual;

	double euribor;

	long fixId;
	long fixAmortId;
	double sensi;

	if (ARMLOCAL_ISOCCY("USD",C_result) == ARM_OK)
	{
		ccyUSDId = C_result.getLong();
	}
	else
	{
		ERROR_MSG(C_result.getMsg(),pRet);
		return S_OK;
	}

	if (ARMLOCAL_ISOCCY("EUR",C_result) == ARM_OK)
	{
		ccyEURId = C_result.getLong();
	}
	else
	{
		ERROR_MSG(C_result.getMsg(),pRet);
		return S_OK;
	}

	if (ARMLOCAL_IRINDEX(KACTUAL_360,
						 resetFreqId,
						 -1.0,
						 K_COMP_PROP,
						 K_MOD_FOLLOWING,
						 timingId,
						 gap,
						 K_ARREARS,
						 10000.,
						 false,
						 "USD",
						 sjUSDTypeId,
						 K_COMP_PROP,
						 K_ADJUSTED,
						 resetFreqId,
						 C_result) == ARM_OK) 
	{
		irIndexId = C_result.getLong();
	}
	else
	{
		ERROR_MSG(C_result.getMsg(),pRet);
		return S_OK;
	}

	if (ARMLOCAL_SWAPLEG(irIndexId,
						 startDate,
						 endMaturityStruct,
						 K_RCV,
						 0,
						 pSpdPostBar / 100.,
						 false,
						 "EUR",
						 KACTUAL_360,
						 gap,
						 "GBP",
						 "GBP",
						 1,
						 K_NX_NONE,
						 K_SHORTSTART,
						 -1.0,
						 K_YES,
						  -1,
						 C_result) == ARM_OK) 
	{
		liborlegId = C_result.getLong();
	}
	else
	{
		ERROR_MSG(C_result.getMsg(),pRet);
		return S_OK;
	}

	if (ARMLOCAL_ClonedAndSetNotional(liborlegId,
									  amortId,
									  100.,
									  C_result) == ARM_OK)
	{
		clonedLiborlegid = C_result.getLong();
	}
	else
	{
		ERROR_MSG(C_result.getMsg(),pRet);
		return S_OK;
	}

	if (ARMLOCAL_CAPFLOOR(clonedLiborlegid,
						  K_CAP,
						  0,
						  pBarriere,
						  C_result) == ARM_OK)
	{
		capBId = C_result.getLong();
	}
	else
	{
		ERROR_MSG(C_result.getMsg(),pRet);
		return S_OK;
	}

	if (ARMLOCAL_ARM_Price(capBId,
						   LocalGetNumObjectId(l_BsModId),
						   C_result) == ARM_OK)
	{
		capB = C_result.getDouble();
	}
	else
	{
		ERROR_MSG(C_result.getMsg(),pRet);
		return S_OK;
	}

	if (ARMLOCAL_ARM_Price(capBId,
						   LocalGetNumObjectId(l_BsModVolUSDId),
						   C_result) == ARM_OK)
	{
		capB_volUSD = C_result.getDouble() - capB;
	}
	else
	{
		ERROR_MSG(C_result.getMsg(),pRet);
		return S_OK;
	}

	if (ARMLOCAL_ARM_Price(capBId,
						   LocalGetNumObjectId(l_BsModCorrPlusId),
						   C_result) == ARM_OK)
	{
		capB_corrplus = C_result.getDouble() - capB;
	}
	else
	{
		ERROR_MSG(C_result.getMsg(),pRet);
		return S_OK;
	}

	if (ARMLOCAL_ARM_Price(capBId,
						   LocalGetNumObjectId(l_BsModCorrMoinsId),
						   C_result) == ARM_OK)
	{
		capB_corrmoins = C_result.getDouble() - capB;
	}
	else
	{
		ERROR_MSG(C_result.getMsg(),pRet);
		return S_OK;
	}

	retCode = ARMLOCAL_FreeObject(capBId,C_result);

	if (ARMLOCAL_DIGITAL(clonedLiborlegid,
						 K_CAP,
						 0,
						 pBarriere,
						 -0.01,
						 0.01,
						 0,
						 1.0,
						 C_result) == ARM_OK)
	{
		digitId = C_result.getLong();
	}
	else
	{
		ERROR_MSG(C_result.getMsg(),pRet);
		return S_OK;
	}

	if (ARMLOCAL_ARM_Price(digitId,
						   LocalGetNumObjectId(l_BsModId),
						   C_result) == ARM_OK)
	{
		digit = C_result.getDouble();
	}
	else
	{
		ERROR_MSG(C_result.getMsg(),pRet);
		return S_OK;
	}

	if (ARMLOCAL_ARM_Price(digitId,
						   LocalGetNumObjectId(l_BsModVolUSDId),
						   C_result) == ARM_OK)
	{
		digit_volUSD = C_result.getDouble() - digit;
	}
	else
	{
		ERROR_MSG(C_result.getMsg(),pRet);
		return S_OK;
	}

	if (ARMLOCAL_ARM_Price(digitId,
						   LocalGetNumObjectId(l_BsModCorrPlusId),
						   C_result) == ARM_OK)
	{
		digit_corrplus = C_result.getDouble() - digit;
	}
	else
	{
		ERROR_MSG(C_result.getMsg(),pRet);
		return S_OK;
	}

	if (ARMLOCAL_ARM_Price(digitId,
						   LocalGetNumObjectId(l_BsModCorrMoinsId),
						   C_result) == ARM_OK)
	{
		digit_corrmoins = C_result.getDouble() - digit;
	}
	else
	{
		ERROR_MSG(C_result.getMsg(),pRet);
		return S_OK;
	}

	retCode = ARMLOCAL_FreeObject(digitId,C_result);
	retCode = ARMLOCAL_FreeObject(clonedLiborlegid,C_result);
	retCode = ARMLOCAL_FreeObject(liborlegId,C_result);
	retCode = ARMLOCAL_FreeObject(irIndexId,C_result);

	if (ARMLOCAL_DUALCAP (startDate,
						  endMaturityStruct,
						  K_CAP,
						  0,
						  pBarriere,
						  liborTypeId,
						  sjUSDTypeId,
						  KACTUAL_360,
						  resetFreqId,
						  timingId,
						  timingId,
						  gap,
						  gap,
						  ccyEURId,
						  ccyUSDId,
						  ccyEURId,
						  "NULL",
						  "NULL",
						  C_result)  == ARM_OK)
	{
		dualcapId = C_result.getLong();
	}
	else
	{
		ERROR_MSG(C_result.getMsg(),pRet);
		return S_OK;
	}

	if (ARMLOCAL_ClonedAndSetNotional(dualcapId,
									  amortId,
									  100.,
									  C_result) == ARM_OK)
	{
		dualclonedId = C_result.getLong();
	}
	else
	{
		ERROR_MSG(C_result.getMsg(),pRet);
		return S_OK;
	}

	if (ARMLOCAL_ARM_Price(dualclonedId,
						   LocalGetNumObjectId(l_CrossModId),
						   C_result) == ARM_OK)
	{
		dual = C_result.getDouble();
	}
	else
	{
		ERROR_MSG(C_result.getMsg(),pRet);
		return S_OK;
	}

	retCode = ARMLOCAL_FreeObject(dualclonedId,C_result);
	retCode = ARMLOCAL_FreeObject(dualcapId,C_result);

	if (ARMLOCAL_LIBORLEG(startDate,
						  endMaturityStruct,
						  liborTypeId,
						  K_RCV,
						  0,
						  0.,
						  resetFreqId,
						  resetFreqId,
						  timingId,
						  K_ARREARS,
						  false,
						  "EUR",
						  K_ADJUSTED,
						  gap,
						  "NULL",
						  "NULL",
						  1,
						  K_NX_NONE,
						  K_SHORTSTART,
						  -1.0,
						  K_YES,
						  -1,
						  C_result) == ARM_OK) 
	{
		liborlegId = C_result.getLong();
	}
	else
	{
		ERROR_MSG(C_result.getMsg(),pRet);
		return S_OK;
	}

	if (ARMLOCAL_ClonedAndSetNotional(liborlegId,
									  amortId,
									  100.,
									  C_result) == ARM_OK)
	{
		clonedLiborlegid = C_result.getLong();
	}
	else
	{
		ERROR_MSG(C_result.getMsg(),pRet);
		return S_OK;
	}

	if (ARMLOCAL_ARM_Price(clonedLiborlegid,
						   LocalGetNumObjectId(l_BsModId),
						   C_result) == ARM_OK)
	{
		euribor = C_result.getDouble();
	}
	else
	{
		ERROR_MSG(C_result.getMsg(),pRet);
		return S_OK;
	}

	retCode = ARMLOCAL_FreeObject(clonedLiborlegid,C_result);
	retCode = ARMLOCAL_FreeObject(liborlegId,C_result);

	if (ARMLOCAL_FIXEDLEG(startDate,
						  endMaturityStruct,
						  K_RCV,
						  0,
						  1.0,
						  KACTUAL_360,
						  resetFreqId,
						  resetFreqId,
						  timingId,
						  K_ADJUSTED,
						  K_SHORTSTART,
						  false,
						  "EUR",
						  "NULL",
						  K_NX_NONE,
						  -1.0,
						  K_YES,
						  -1,
						  C_result) == ARM_OK) 
	{
		fixId = C_result.getLong();
	}

	if (ARMLOCAL_ClonedAndSetNotional(fixId,
									  amortId,
									  100.,
									  C_result) == ARM_OK)
	{
		fixAmortId = C_result.getLong();
	}

	if (ARMLOCAL_ARM_Price(fixAmortId,
						   LocalGetNumObjectId(l_BsModId),
						   C_result) == ARM_OK)
	{
		sensi = C_result.getDouble();
	}
	else
	{
		ERROR_MSG(C_result.getMsg(),pRet);
		return S_OK;
	}

	retCode = ARMLOCAL_FreeObject(fixId,C_result);
	retCode = ARMLOCAL_FreeObject(fixAmortId,C_result);

	double periode2_price = 0.0;

	// 2eme période
	if (endMaturityTot > endMaturityStruct)
	{
		if (ARMLOCAL_LIBORLEG(endMaturityStruct,
							  endMaturityTot,
							  liborTypeId,
							  K_RCV,
							  0,
							  pSpd2phase / 100.,
							  resetFreqId,
							  resetFreqId,
							  K_ADVANCE,
							  K_ARREARS,
							  false,
							  "EUR",
							  K_ADJUSTED,
							  -2.0,
							  "NULL",
							  "NULL",
							  1,
							  K_NX_NONE,
							  K_SHORTSTART,
							  -1.0,
							  K_YES,
							  -1,
							  C_result) == ARM_OK) 
		{
			liborlegId = C_result.getLong();
		}
		else
		{
			ERROR_MSG(C_result.getMsg(),pRet);
			return S_OK;
		}

		if (ARMLOCAL_ClonedAndSetNotional(liborlegId,
										  amortId,
										  100.,
										  C_result) == ARM_OK)
		{
			clonedLiborlegid = C_result.getLong();
		}
		else
		{
			ERROR_MSG(C_result.getMsg(),pRet);
			return S_OK;
		}

		if (ARMLOCAL_ARM_Price(clonedLiborlegid,
							   LocalGetNumObjectId(l_YcModId),
							   C_result) == ARM_OK)
		{
			periode2_price = C_result.getDouble();
		}
		else
		{
			ERROR_MSG(C_result.getMsg(),pRet);
			return S_OK;
		}

		retCode = ARMLOCAL_FreeObject(liborlegId,C_result);
		retCode = ARMLOCAL_FreeObject(clonedLiborlegid,C_result);

	}

	retCode = ARMLOCAL_FreeObject(amortIdToDestroy,C_result);
	retCode = ARMLOCAL_FreeObject(ccyUSDId,C_result);
	retCode = ARMLOCAL_FreeObject(ccyEURId,C_result);

	// calcul spr post
	double delta_NPV;
	double fund_total;
	double probaM;
	double marge_capB;
	double marge_digit;
	double prix;

	delta_NPV = funding2Price + periode2_price;
	fund_total = fund + delta_NPV - pSoulte;

	if (ARMLOCAL_ComputeVolatility(LocalGetNumObjectId(l_ProbaMarge),
								   digit / sensi * 100.0,
								   0.0,
								   0.0,
								   C_result) == ARM_OK)
	{
		probaM = C_result.getDouble();
	}
	else
	{
		ERROR_MSG(C_result.getMsg(),pRet);
		return S_OK;
	}

	marge_capB = capB_volUSD + max(capB_corrplus, capB_corrmoins);
	marge_digit = digit_volUSD + max(digit_corrplus, digit_corrmoins);

	prix = -((fund_total + euribor + capB - dual - marge_capB) / (sensi - (digit - marge_digit - probaM))) * 100.;

	VECTOR<double> vRes;

	vRes.push_back(prix);
	vRes.push_back(the_i);
	vRes.push_back(startDate);
	vRes.push_back(endMaturityStruct);
	vRes.push_back(endMaturityTot);

	Res = VECTORDOUBLE2VARIANT(vRes, pRet);

	_variant_t wrap_Res;
	wrap_Res.Attach(*pRet);
	*pRet=wrap_Res.Detach();

	return S_OK;
}
catch(...)
{
	return E_FAIL;
}
}




STDMETHODIMP ActiveXModule::ARMcptDigital(double pAsOf,
								   double pStartDate,
								   double pEndDate,
								   double pNtl,
								   BSTR pIndex,
								   BSTR pBsmod,
								   BSTR pBsmodDelta,
								   BSTR pBsmodVega,
								   BSTR pFreqP,
								   BSTR pResetTiming,
								   BSTR pPorR,
								   BSTR pCcy,
								   BSTR pCcyIdx,
								   BSTR pDayCount,
								   BSTR pCapOrFloor,
								   BSTR pAmort,
								   BSTR pStrike,
								   BSTR pPayOff,
								   BSTR pSpd,
								   double pResetGap,
								   double pSpreadBelow,
								   double pSpreadAbove,
								   BSTR pFwdRule,
								   BSTR pIntRule,
								   BSTR pStubRule,
								   BSTR pFreqAmort,
								   double pTxAmort,
								   double pAmountAmort,
								   double pRefDate,
								   VARIANT *pRet)

{
	// TODO: Add your implementation code here
try
{
//	double vDate = pDate;

	ARM_result C_result;

	double AsOf = pAsOf;
	double startdate = pStartDate;
	double enddate = pEndDate;
	double Ntl = pNtl;
	double resetgap = pResetGap;
	double spreadbelow = pSpreadBelow;
	double spreadabove = pSpreadAbove;
	double tauxAmort = pTxAmort;
	double amountAmort = pAmountAmort;
	double refdate = pRefDate;

// en attendant le traitement des StepUp

	_bstr_t bspread = pSpd;
	CCString spread = bspread;
		
	_bstr_t bstrike = pStrike;
	CCString strike = bstrike;
		
	_bstr_t bpayoff = pPayOff;
	CCString payoff = bpayoff;
		

// *************************

	long retCode;
	double Res = 0.;

	double StrikeProba;


	_bstr_t bindex(pIndex);
	CCString l_index = bindex;

	_bstr_t bfreqp(pFreqP);
	CCString l_freqp = bfreqp;

	_bstr_t bresettiming(pResetTiming);
	CCString l_resettiming = bresettiming;

	_bstr_t bporr(pPorR);
	CCString l_porr = bporr;

	_bstr_t bccy(pCcy);
	CCString l_ccy = bccy;

	_bstr_t bccyidx(pCcyIdx);
	CCString l_ccyidx = bccyidx;

	_bstr_t bdaycount(pDayCount);
	CCString l_daycount = bdaycount;

	_bstr_t bcaporfloor(pCapOrFloor);
	CCString l_caporfloor = bcaporfloor;

	_bstr_t bamort(pAmort);
	CCString l_amort = bamort;

	_bstr_t bbsmod(pBsmod);
	CCString l_bsmod = bbsmod;

	_bstr_t bbsmodDelta(pBsmodDelta);
	CCString l_bsmodDelta = bbsmodDelta;

	_bstr_t bbsmodVega(pBsmodVega);
	CCString l_bsmodVega = bbsmodVega;

	_bstr_t bfwdrule(pFwdRule);
	CCString l_fwdrule = bfwdrule;

	_bstr_t bintrule(pIntRule);
	CCString l_intrule = bintrule;

	_bstr_t bstubrule(pStubRule);
	CCString l_stubrule = bstubrule;

	_bstr_t bFreqAmort(pFreqAmort);
	CCString l_FreqAmort = bFreqAmort;




	long freqId;

	if((freqId = ARM_ConvIrIndNameToFreq(l_index , C_result)) == ARM_DEFAULT_ERR)
	{
		ERROR_MSG("Invalid Index",pRet);
		return S_OK;
	}

	long indexId;

/*	if((indexId = ARM_ConvIrIndName(l_index , C_result)) == ARM_DEFAULT_ERR)
	{
		ERROR_MSG("Invalid Index",pRet);
		return S_OK;
	}
*/	
	indexId = ARM_ConvIrType(l_index);
	
	long resetFreqId;

	if((resetFreqId = ARM_ConvFrequency(l_freqp, C_result)) == ARM_DEFAULT_ERR)
	{
		ERROR_MSG("Invalid Frq",pRet,ARM_ERROR_FREQ);
		return S_OK;
	}


	long porrId;

	if((porrId = ARM_ConvRecOrPay(l_porr, C_result)) == ARM_DEFAULT_ERR)
	{
		ERROR_MSG("Invalid PorR",pRet);
		return S_OK;
	}

	long caporfloorId;

	if((caporfloorId = ARM_ConvCapOrFloor(l_caporfloor, C_result)) == ARM_DEFAULT_ERR)
	{
		ERROR_MSG("Invalid CaporFloor",pRet);
		return S_OK;
	}

	long amortFreqId;

	if((amortFreqId = ARM_ConvFrequency(l_FreqAmort, C_result)) == ARM_DEFAULT_ERR)
	{
		ERROR_MSG("Invalid Amort Frq",pRet,ARM_ERROR_FREQ);
		return S_OK;
	}

	

	long daycountId= ARM_ConvDayCount(l_daycount);

	long intruleId= ARM_ConvIntRule(l_intrule);

	long resettimingId = ARM_ConvPayResetRule (l_resettiming);

	long stubruleId = ARM_ConvStubRule(l_stubrule);


	// AmortMethodId
	long amortMethodId;

	if((amortMethodId = ARM_ConvAmortMethod (l_amort, C_result)) == ARM_DEFAULT_ERR)
	{
		ERROR_MSG("Invalid AmortMethod",pRet,ARM_ERROR_FREQ);
		return S_OK;
	}

	long irindexId;

	if (ARMLOCAL_IRINDEX(KACTUAL_360,
						 resetFreqId,
						 -1,
						 K_COMP_PROP,
						 K_MOD_FOLLOWING,
						 resettimingId,
						 resetgap,
						 K_ARREARS,
						 10000.,
						 false,
						 l_ccyidx,
						 indexId,
						 K_COMP_PROP,
						 intruleId,
						 resetFreqId,
						 C_result) == ARM_OK) 
	{
		irindexId = C_result.getLong();
	}
	else
	{
		ERROR_MSG(C_result.getMsg(),pRet);
		return S_OK;
	}

	long swaplegId;
	if (ARMLOCAL_SWAPLEG(irindexId,
						 startdate,
						 enddate,
						 porrId,
						 1,
						 LocalGetNumObjectId(spread),
						 false,
						 l_ccy,
						 daycountId,
						 resetgap,
						 "NULL",
						 "NULL",
						 0,
						 K_NX_NONE,
						 stubruleId,
						 refdate,
						 K_NO,
						  -1,
						 C_result) == ARM_OK) 
	{
		swaplegId = C_result.getLong();
	}
	else
	{
		ERROR_MSG(C_result.getMsg(),pRet);
		return S_OK;
	}


	long genAmortId;

	if (ARMLOCAL_ARM_GenAmortization(swaplegId,
									 amortMethodId,
									 amortFreqId,
									 amountAmort,
									 -1,
									 Ntl,
									 tauxAmort,
									 0.,
									 LocalGetNumObjectId(l_bsmod),
									 0.0,
									 C_result) == ARM_OK) 
	{
		genAmortId = C_result.getLong();
	}
	else
	{
		ERROR_MSG(C_result.getMsg(),pRet);
		return S_OK;
	}

	long cloneId;
	if (ARMLOCAL_ClonedAndSetNotional(swaplegId,
									  genAmortId,
									  100.,
									  C_result) == ARM_OK)
	{
		cloneId = C_result.getLong();
	}
	else
	{
		ERROR_MSG(C_result.getMsg(),pRet);
		return S_OK;
	}

	long digitalId;
	if (ARMLOCAL_DIGITAL(cloneId,
					     caporfloorId,
						 1,
						 LocalGetNumObjectId(strike),
						 spreadbelow,
						 spreadabove,
						 1,
						 LocalGetNumObjectId(payoff),
						 C_result) == ARM_OK)
	{
		digitalId = C_result.getLong();
	}
	else
	{
		ERROR_MSG(C_result.getMsg(),pRet);
		return S_OK;
	}

	long digitalIdproba;

	// 1 pour le CAP
	if ( caporfloorId = 1 )
	{
		StrikeProba = 0.01;
	}
	else
	{
		StrikeProba = 1000;
	}
	if (ARMLOCAL_DIGITAL(cloneId,
						 caporfloorId,
						 0,
						 StrikeProba,
						 spreadbelow,
						 spreadabove,
						 1,
						 LocalGetNumObjectId(payoff),
						 C_result) == ARM_OK)
	{
		digitalIdproba = C_result.getLong();
	}
	else
	{
		ERROR_MSG(C_result.getMsg(),pRet);
		 S_OK;
	}


	/**************** PRICE *******************/


	double pxDigit;
	double pxDigitProba;
	double pxDelta;
	double pxVega;
	double proba;
	double Sensi;

	if (ARMLOCAL_ARM_Price(digitalId,
						   LocalGetNumObjectId(l_bsmod),
						   C_result) == ARM_OK)
	{
		pxDigit = C_result.getDouble();
	}
	else
	{
		ERROR_MSG(C_result.getMsg(),pRet);
		return S_OK;
	}


	if (ARMLOCAL_ARM_Price(digitalIdproba,
						   LocalGetNumObjectId(l_bsmod),
						   C_result) == ARM_OK)
	{
		pxDigitProba = C_result.getDouble();
	}
	else
	{
		ERROR_MSG(C_result.getMsg(),pRet);
		return S_OK;
	}


	if (ARMLOCAL_ARM_Price(digitalId,
						   LocalGetNumObjectId(l_bsmodDelta),
						   C_result) == ARM_OK)
	{
		pxDelta = C_result.getDouble();
	}
	else
	{
		ERROR_MSG(C_result.getMsg(),pRet);
		return S_OK;
	}


	if (ARMLOCAL_ARM_Price(digitalId,
						   LocalGetNumObjectId(l_bsmodVega),
						   C_result) == ARM_OK)
	{
		pxVega = C_result.getDouble();
	}
	else
	{
		ERROR_MSG(C_result.getMsg(),pRet);
		return S_OK;
	}

	Sensi = abs( pxVega - pxDigit ) + abs( pxDelta - pxDigit ) ;


	//**********************************************
	//   *******  calcul marge Technique ******
	//**********************************************

	double minProba;
	double maxProba;
	double margetechnique;

	proba = pxDigit / pxDigitProba ;

	if ( proba < (1 - proba ) )
	{
		minProba = proba;
	}
	else
	{
		minProba = 1 - proba;
	}

	if ( minProba - 0.1  < 0 )
	{
		maxProba = 0;
	}
	else
	{
		maxProba = minProba - 0.1;
	}

	if ( pxDigitProba * ( 0.01 + 0.05 * maxProba) < 0 )
	{
		margetechnique =  - pxDigitProba * ( 0.01 + 0.05 * maxProba);
	}
	else
	{
		margetechnique =   pxDigitProba * ( 0.01 + 0.05 * maxProba);
	}

	
	//**********************************************
	//**********************************************
		
		
	retCode = ARMLOCAL_FreeObject(irindexId,C_result);
	retCode = ARMLOCAL_FreeObject(swaplegId,C_result);
	retCode = ARMLOCAL_FreeObject(digitalId,C_result);
	retCode = ARMLOCAL_FreeObject(digitalIdproba,C_result);
	retCode = ARMLOCAL_FreeObject(cloneId,C_result);
	retCode = ARMLOCAL_FreeObject(genAmortId,C_result);

	VECTOR<double> vRes;

	vRes.push_back(pxDigit);
	vRes.push_back(margetechnique);
	vRes.push_back(Sensi);

	Res = VECTORDOUBLE2VARIANT(vRes, pRet);

	

	_variant_t wrap_Res;
	wrap_Res.Attach(*pRet);
	*pRet=wrap_Res.Detach();

//	*pRet = NpvDigit



	return S_OK;
}
catch(...)
{
	return E_FAIL;
}
}



VECTOR<double> computeQTF(double pAsOf,
				  double StartDate,
				  double EndDate,
				  double Ntl,
				  long IndexId,
				  CCString BsmodId,
				  CCString BsmoddeltaIdCcy1,
				  CCString BsmodvegaIdCcy1,
				  CCString BsmoddeltaIdCcy2,
				  CCString BsmodvegaIdCcy2,
				  CCString BsmodFxCorrelId,
	   			  long FreqPId,
	   			  long FreqRId,
				  long ResetTimingId,
				  const CCString& CcyStr,
				  const CCString& CcyIdxStr,
				  long DayCountId,
				  long isRateFixed,
				  double FixedRate,
				  double TauxBoni,
				  CCString Barrier,
				  CCString SpdTf,
				  double ResetGap,
				  double SpreadBelow,
				  double SpreadAbove,
				  long FwdRuleId,
				  long IntRuleId,
				  long StubRuleId,
				  long AmortRefVal,
				  double RefDate)

{

	ARM_result C_result;

	VECTOR<double> resultat;
	resultat.clear();

// *************************
	double prix ;

	double pxFixedLeg;
	double pxFixedLegDelta;

	double pxCap;
	double pxCapDeltaCcy1;
	double pxCapDeltaCcy2;
	double pxCapVegaCcy1;
	double pxCapVegaCcy2;
	double pxCapFxCorrel;

	double pxDigit;
	double pxDigitDeltaCcy1;
	double pxDigitDeltaCcy2;
	double pxDigitVegaCcy1;
	double pxDigitVegaCcy2;
	double pxDigitFxCorrel;

	double pxDigitproba;
	double margetechnique;
	double margedeltaCcy1;
	double margedeltaCcy2;
	double margevega;
	double margeFxCorrel;
	double maxProba;
	double minProba;
	double proba;



	long retCode;
	double Res = 0.;

	//  Index de reference	
	long irindexId;

	if (ARMLOCAL_IRINDEX(KACTUAL_360,
						 FreqPId,
						 -1,
						 K_COMP_PROP,
						 K_MOD_FOLLOWING,
						 ResetTimingId,
						 ResetGap,
						 K_ARREARS,
						 10000.,
						 false,
						 CcyIdxStr,
						 IndexId,
						 K_COMP_PROP,
						 IntRuleId,
						 FreqRId,
						 C_result) == ARM_OK) 
	{
		irindexId = C_result.getLong();
	}
	else
	{
		return resultat;
	}

	

	//  Index de pour la fixedleg  pas de QUANTO pour la pate fixe
	
	long irindexIdfixed;

	if (ARMLOCAL_IRINDEX(KACTUAL_360,
						 FreqPId,
						 -1,
						 K_COMP_PROP,
						 K_MOD_FOLLOWING,
						 ResetTimingId,
						 ResetGap,
						 K_ARREARS,
						 10000.,
						 false,
						 CcyStr,
						 K_FIXED,
						 K_COMP_PROP,
						 IntRuleId,
						 FreqRId,
						 C_result) == ARM_OK) 
	{
		irindexIdfixed = C_result.getLong();
	}
	else
	{
		return resultat;
	}

	// fixed leg

	long fixedlegId;
	if (ARMLOCAL_SWAPLEG(irindexIdfixed,
						 StartDate,
						 EndDate,
						 K_RCV,
						 0,
						 TauxBoni,
						 false,
						 CcyStr,
						 DayCountId,
						 ResetGap,
						 "NULL",
						 "NULL",
						 0,
						 K_NX_NONE,
						 StubRuleId,
						 RefDate,
						 K_NO,
						  -1,
						 C_result) == ARM_OK) 
	{
		fixedlegId = C_result.getLong();
	}
	else
	{
		return resultat;
	}

	//   swapleg de reference

	long swaplegId;
	if (ARMLOCAL_SWAPLEG(irindexId,
						 StartDate,
						 EndDate,
						 K_RCV,
						 0,
						 0,
						 false,
						 CcyStr,
						 DayCountId,
						 ResetGap,
						 "NULL",
						 "NULL",
						 0,
						 K_NX_NONE,
						 StubRuleId,
						 RefDate,
						 K_NO,
						  -1,
						 C_result) == ARM_OK) 
	{
		swaplegId = C_result.getLong();
	}
	else
	{
		return resultat;
	}




	// clone de ref

	long cloneId;
	if (ARMLOCAL_ClonedAndSetNotional(swaplegId,
									  AmortRefVal,
									  100.,
									  C_result) == ARM_OK)
	{
		cloneId = C_result.getLong();
	}
	else
	{
		return resultat;
	}

	long cloneIdfixed;
	if (ARMLOCAL_ClonedAndSetNotional(fixedlegId,
									  AmortRefVal,
									  100.,
									  C_result) == ARM_OK)
	{
		cloneIdfixed = C_result.getLong();
	}
	else
	{
		return resultat;
	}

	long capBId;


	//pas de CAP si on a un taux fixe sur la phase structurée
	if (isRateFixed==K_NO)
	{
		if (ARMLOCAL_CAPFLOOR(cloneId,
							  K_CAP,
							  1,
							  LocalGetNumObjectId(Barrier),
							  C_result) == ARM_OK)
		{
			capBId = C_result.getLong(); 
		}
		else
		{
			return resultat;
		}

		
	}

	// ***** compute payoff ****

	long payoffrebatesId1;
	long payoffrebatesId2;
	long refvalTxBoni;
	long refvalFixedRate;

	if (ARMLOCAL_CONSTREFVALUE(TauxBoni,
						  C_result) == ARM_OK)
	{
		refvalTxBoni = C_result.getLong();
	}
	else
	{
		return resultat;
	}

	if (isRateFixed==K_YES)
	{

		if (ARMLOCAL_CONSTREFVALUE(FixedRate,
					  C_result) == ARM_OK)
		{
			refvalFixedRate = C_result.getLong();
		}
		else
		{
			return resultat;
		}

		if (ARMLOCAL_SumRefValue(refvalFixedRate,
						  refvalTxBoni,
						  -1,
						  C_result) == ARM_OK)
		{
			payoffrebatesId1 = C_result.getLong();
		}
		else
		{
			return resultat;
		}

	}
	else
	{
		if (ARMLOCAL_SumRefValue(LocalGetNumObjectId(Barrier),
						  refvalTxBoni,
						  -1,
						  C_result) == ARM_OK)
		{
			payoffrebatesId1 = C_result.getLong();
		}
		else
		{
			return resultat;
		}
	}

	

	if (ARMLOCAL_SumRefValue(payoffrebatesId1,
						  LocalGetNumObjectId(SpdTf),
						  1,
						  C_result) == ARM_OK)
	{
		payoffrebatesId2 = C_result.getLong();
	}
	else
	{
		return resultat;
	}


	// digital leg de ref
			
	long digitalId;
	if (ARMLOCAL_DIGITAL(cloneId,
					     K_CAP,
						 1,
						 LocalGetNumObjectId(Barrier),
						 SpreadBelow,
						 SpreadAbove,
						 1,
						 payoffrebatesId2,
						C_result) == ARM_OK)
	{
		digitalId = C_result.getLong();
	}
	else
	{
		return resultat;
	}
			
	// digital leg pour la proba
			
	long digitalIdproba;
	if (ARMLOCAL_DIGITAL(cloneId,
					     K_CAP,
						 0,
						 -100,
						 SpreadBelow,
						 SpreadAbove,
						 1,
						 payoffrebatesId2,
						C_result) == ARM_OK)
	{
		digitalIdproba = C_result.getLong();
	}
	else
	{
		return resultat;
	}


			

	// ******** price Cap leg
	//pas de CAP si on a un taux fixe sur la phase structurée
	if (isRateFixed==K_NO)
	{
		if (ARMLOCAL_ARM_Price(capBId ,
							   LocalGetNumObjectId(BsmodId),
							   C_result) == ARM_OK)
		{
			pxCap = C_result.getDouble();
		}
		else
		{
			return resultat;
		}

		#ifdef _DEBUG
			Res = ARMLOCAL_ARM_View(capBId,C_result);
		#endif

		#ifdef _DEBUG
			Res = ARMLOCAL_ARM_View(LocalGetNumObjectId(BsmodId),C_result);
		#endif

		#ifdef _DEBUG
			Res = ARMLOCAL_ARM_View(cloneId,C_result);
		#endif

		#ifdef _DEBUG
			Res = ARMLOCAL_ARM_View(swaplegId,C_result);
		#endif

		#ifdef _DEBUG
			Res = ARMLOCAL_ARM_View(AmortRefVal,C_result);
		#endif

		#ifdef _DEBUG
			Res = ARMLOCAL_ARM_View(irindexId,C_result);
		#endif

		if (ARMLOCAL_ARM_Price(capBId ,
							   LocalGetNumObjectId(BsmoddeltaIdCcy1),
							   C_result) == ARM_OK)
		{
			pxCapDeltaCcy1 = C_result.getDouble();
		}
		else
		{
			return resultat;
		}
		if (ARMLOCAL_ARM_Price(capBId ,
							   LocalGetNumObjectId(BsmoddeltaIdCcy2),
							   C_result) == ARM_OK)
		{
			pxCapDeltaCcy2 = C_result.getDouble();
		}
		else
		{
			return resultat;
		}

		if (ARMLOCAL_ARM_Price(capBId ,
							   LocalGetNumObjectId(BsmodvegaIdCcy1),
							   C_result) == ARM_OK)
		{
			pxCapVegaCcy1 = C_result.getDouble();
		}
		else
		{
			return resultat;
		}

		if (ARMLOCAL_ARM_Price(capBId ,
							   LocalGetNumObjectId(BsmodvegaIdCcy2),
							   C_result) == ARM_OK)
		{
			pxCapVegaCcy2 = C_result.getDouble();
		}
		else
		{
			return resultat;
		}

		if (ARMLOCAL_ARM_Price(capBId ,
							   LocalGetNumObjectId(BsmodFxCorrelId ),
							   C_result) == ARM_OK)
		{
			pxCapFxCorrel = C_result.getDouble();
		}
		else
		{
			return resultat;
		}
	}
	else
	{
		pxCap=0;
		pxCapDeltaCcy1=0;
		pxCapDeltaCcy2=0;
		pxCapVegaCcy1=0;
		pxCapVegaCcy2=0;
		pxCapFxCorrel=0;
	}

	// ******** price digit leg

	if (ARMLOCAL_ARM_Price(digitalId,
						   LocalGetNumObjectId(BsmodId),
						   C_result) == ARM_OK)
	{
		pxDigit = C_result.getDouble();
	}
	else
	{
		return resultat;
	}

	#ifdef _DEBUG
		Res = ARMLOCAL_ARM_View(digitalId,C_result);
	#endif

	#ifdef _DEBUG
		Res = ARMLOCAL_ARM_View(LocalGetNumObjectId(BsmodId),C_result);
	#endif

	#ifdef _DEBUG
		Res = ARMLOCAL_ARM_View(cloneId,C_result);
	#endif

	#ifdef _DEBUG
		Res = ARMLOCAL_ARM_View(payoffrebatesId2,C_result);
	#endif


	if (ARMLOCAL_ARM_Price(digitalId,
						   LocalGetNumObjectId(BsmoddeltaIdCcy1),
						   C_result) == ARM_OK)
	{
		pxDigitDeltaCcy1 = C_result.getDouble();
	}
	else
	{
		return resultat;
	}

	if (ARMLOCAL_ARM_Price(digitalId,
						   LocalGetNumObjectId(BsmoddeltaIdCcy2),
						   C_result) == ARM_OK)
	{
		pxDigitDeltaCcy2 = C_result.getDouble();
	}
	else
	{
		return resultat;
	}

	if (ARMLOCAL_ARM_Price(digitalId,
						   LocalGetNumObjectId(BsmodvegaIdCcy1),
						   C_result) == ARM_OK)
	{
		pxDigitVegaCcy1 = C_result.getDouble();
	}
	else
	{
		return resultat;
	}

	if (ARMLOCAL_ARM_Price(digitalId,
						   LocalGetNumObjectId(BsmodvegaIdCcy2),
						   C_result) == ARM_OK)
	{
		pxDigitVegaCcy2 = C_result.getDouble();
	}
	else
	{
		return resultat;
	}

	if (ARMLOCAL_ARM_Price(digitalId,
						   LocalGetNumObjectId(BsmodFxCorrelId),
						   C_result) == ARM_OK)
	{
		pxDigitFxCorrel = C_result.getDouble();
	}
	else
	{
		return resultat;
	}

	// ******** price digit leg proba

	if (ARMLOCAL_ARM_Price(digitalIdproba,
						   LocalGetNumObjectId(BsmodId),
						   C_result) == ARM_OK)
	{
		pxDigitproba = C_result.getDouble();
	}
	else
	{
		return resultat;
	}

	#ifdef _DEBUG
		Res = ARMLOCAL_ARM_View(digitalIdproba,C_result);
	#endif
	
	// ******** price fixed leg

	if (ARMLOCAL_ARM_Price(cloneIdfixed,
						   LocalGetNumObjectId(BsmodId),
						   C_result) == ARM_OK)
	{
		pxFixedLeg = C_result.getDouble();
	}
	else
	{
		return resultat;
	}

	if (ARMLOCAL_ARM_Price(cloneIdfixed,
						   LocalGetNumObjectId(BsmoddeltaIdCcy1),
						   C_result) == ARM_OK)
	{
		pxFixedLegDelta = C_result.getDouble();
	}
	else
	{
		return resultat;
	}


	// *****  calcul marge ******

	proba = pxDigit / pxDigitproba ;

	if ( proba < (1 - proba ) )
	{
		minProba = proba;
	}
	else
	{
		minProba = 1 - proba;
	}
	if 	(IsCMSIndex((ARM_INDEX_TYPE)IndexId) == 1)
	{


		if ( minProba - 0.1  < 0 )
		{
			maxProba = 0;
		}
		else
		{
			maxProba = minProba - 0.1;
		}

		if ( pxDigitproba * ( 0.01 + 0.05 * maxProba) < 0 )
		{
			margetechnique =  - pxDigitproba * ( 0.01 + 0.05 * maxProba);
		}
		else
		{
			margetechnique =   pxDigitproba * ( 0.01 + 0.05 * maxProba);
		}
	}
	else
	{


		if ( minProba  < 0. )
		{
			maxProba = 0. ;
		}
		else
		{
			maxProba = minProba;
		}

		if ( pxDigitproba * ( 3. /50. * maxProba) < 0 )
		{
			margetechnique =  - pxDigitproba * ( 3. /50. * maxProba);
		}
		else
		{
			margetechnique =   pxDigitproba * ( 3. / 50. * maxProba);
		}

	}
	prix = pxFixedLeg + pxCap + pxDigit ;

	margedeltaCcy1 = ( pxCapDeltaCcy1 - pxCap  +  pxDigitDeltaCcy1 - pxDigit +  pxFixedLegDelta - pxFixedLeg ) ;
	margedeltaCcy2 = ( pxCapDeltaCcy2 - pxCap  +  pxDigitDeltaCcy2 - pxDigit ) ;

	margevega = fabs( pxCapVegaCcy1 - pxCap + pxDigitVegaCcy1 - pxDigit ) + fabs( pxCapVegaCcy2 - pxCap + pxDigitVegaCcy2 - pxDigit );

	margeFxCorrel = ( pxCapFxCorrel - pxCap ) + ( pxDigitFxCorrel - pxDigit );

/**************** PRICE *******************/


	retCode = ARMLOCAL_FreeObject(irindexIdfixed,C_result);
	retCode = ARMLOCAL_FreeObject(fixedlegId,C_result);
	retCode = ARMLOCAL_FreeObject(irindexId,C_result);
	retCode = ARMLOCAL_FreeObject(swaplegId,C_result);
	retCode = ARMLOCAL_FreeObject(capBId,C_result);
	retCode = ARMLOCAL_FreeObject(digitalId,C_result);
	retCode = ARMLOCAL_FreeObject(digitalIdproba,C_result);
	retCode = ARMLOCAL_FreeObject(cloneId,C_result);
	retCode = ARMLOCAL_FreeObject(cloneIdfixed,C_result);
	retCode = ARMLOCAL_FreeObject(refvalTxBoni,C_result);
	retCode = ARMLOCAL_FreeObject(payoffrebatesId1,C_result);
	retCode = ARMLOCAL_FreeObject(payoffrebatesId2,C_result);

	VECTOR<double> vRes;

	vRes.push_back(prix);
	vRes.push_back(margetechnique);
	vRes.push_back(margedeltaCcy1);
	vRes.push_back(margedeltaCcy2);
	vRes.push_back(margevega);
	vRes.push_back(margeFxCorrel);


	return vRes;


}


STDMETHODIMP ActiveXModule::ARMcptSPTQTF(double pAsOf,
										   double pStartDate,
										   double pStartDatePhase2,
										   double pStartDatePhase3, 
										   double pEndDate, 
										   double pNtl, 
										   BSTR pIndexPhase1, 
										   BSTR pIndexPhase2, 
										   BSTR pIndexFund,
										   BSTR pIndexPhase3,  
										   BSTR pFreqPPhase1,
										   BSTR pFreqPPhase2, 
										   BSTR pFreqPFund, 
										   BSTR pFreqPPhase3, 
										   BSTR pFreqR, 
										   BSTR pResetTimingPhase1,
										   BSTR pResetTimingPhase2,
										   BSTR pResetTimingPhase3, 
										   BSTR pCcy, 
										   BSTR pCcyIdx, 
										   BSTR pDayCount,
										   double pFee, 
										   BSTR pIsRateFixedPhase2, 
										   double pFixedRatePhase2,
										   BSTR pBarrier, 
										   BSTR pSpdPhase1, 
										   BSTR  pSpdPhase1fund, 
										   BSTR pSpdPhase2Tf,
										   BSTR pSpdPhase2fund,
										   BSTR pSpdPhase3, 
										   BSTR pSpdPhase3fund, 
										   double pResetGapPhase1, 
										   double pResetGapPhase2, 
										   double pResetGapPhase3,
										   BSTR pAmort,
										   BSTR pBsmod,
										   BSTR pBsmodDeltaCcy1,
										   BSTR pBsmodVegaCcy1,
										   BSTR pBsmodDeltaCcy2,
										   BSTR pBsmodVegaCcy2,
										   BSTR pBsmodFxCorrel,
										   BSTR pFwdRule, 
										   BSTR pIntRule,
										   BSTR pStubRule,
										   BSTR pFreqAmort, 
										   double pTxAmort, 
										   double pAmountAmort, 
										   double pRefDate, 
										   VARIANT *pRet)
{
	// TODO: Add your implementation code here
try
{
	ARM_result C_result;

	double AsOf = pAsOf;
	double startdate = pStartDate;
	double startdatephase2 = pStartDatePhase2;
	double startdatephase3 = pStartDatePhase3;
	double enddate = pEndDate;

	double Ntl = pNtl;
	double resetgapphase1 = pResetGapPhase1;
	double resetgapphase2 = pResetGapPhase2;
	double resetgapphase3 = pResetGapPhase3;

	double refdate = pRefDate;
	double fee = pFee;

	double pxFundingLeg;
	double pxFundingLegDelta;
	double SensiFundingLeg;
	double pxPhase3Leg;
	double pxPhase3LegDelta;
	double pxPhase3fundLeg;
	double pxPhase3fundLegDelta;

	double pxPhase1Leg;
	double pxPhase1LegDelta;
	double pxPhase1fundLeg;
	double pxPhase1fundLegDelta;
	
	double pxQTFLeg;
	double margeQTFTech;
	double margeQTFDeltaCcy1;
	double margeQTFDeltaCcy2;
	double margeQTFVega;
	double margeQTFFxCorrel;
	double margeSensi;
	double marge;

	double Xincr;
	double prixAct;
	double prixAct2;
	double prixInter;

	double fixedRatePhase2 = pFixedRatePhase2;

	bool isphase1;
	bool isphase3;

	isphase1=false;
	isphase3=false;

	VECTOR<double> pvQTFLeg;

// en attendant le traitement des StepUp

	_bstr_t bspdphase1 = pSpdPhase1;
	CCString spdphase1 = bspdphase1;
		
	_bstr_t bspdphase1fund = pSpdPhase1fund;
	CCString spdphase1fund = bspdphase1fund;

	_bstr_t bspreadphase2fund = pSpdPhase2fund;
	CCString spreadphase2fund = bspreadphase2fund;
		
	_bstr_t bspdphase3 = pSpdPhase3;
	CCString spdphase3 = bspdphase3;
		
	_bstr_t bspdphase3fund = pSpdPhase3fund;
	CCString spdphase3fund = bspdphase3fund;
		
	_bstr_t bspdphase2tf = pSpdPhase2Tf;
	CCString spdphase2tf = bspdphase2tf;
		
	_bstr_t bbarrier = pBarrier;
	CCString barrier = bbarrier;
		
	long retCode;
	double Res = 0.;

	_bstr_t bindexphase1(pIndexPhase1);
	CCString l_indexphase1 = bindexphase1;

	_bstr_t bindexphase2(pIndexPhase2);
	CCString l_indexphase2 = bindexphase2;

	_bstr_t bindexFund(pIndexFund);
	CCString l_indexFund = bindexFund;

	_bstr_t bindexphase3(pIndexPhase3);
	CCString l_indexphase3 = bindexphase3;

	_bstr_t bfreqpphase1(pFreqPPhase1);
	CCString l_freqpphase1 = bfreqpphase1;

	_bstr_t bfreqpphase2(pFreqPPhase2);
	CCString l_freqpphase2 = bfreqpphase2;

	_bstr_t bfreqpFund(pFreqPFund);
	CCString l_freqpFund = bfreqpFund;

	_bstr_t bfreqpphase3(pFreqPPhase3);
	CCString l_freqpphase3 = bfreqpphase3;

	_bstr_t bfreqr(pFreqR);
	CCString l_freqr = bfreqr;

	_bstr_t bresettimingphase1(pResetTimingPhase1);
	CCString l_resettimingphase1 = bresettimingphase1;

	_bstr_t bresettimingphase2(pResetTimingPhase2);
	CCString l_resettimingphase2 = bresettimingphase2;

	_bstr_t bresettimingphase3(pResetTimingPhase3);
	CCString l_resettimingphase3 = bresettimingphase3;

	_bstr_t bccy(pCcy);
	CCString l_ccy = bccy;

	_bstr_t bccyidx(pCcyIdx);
	CCString l_ccyidx = bccyidx;


	_bstr_t bdaycount(pDayCount);
	CCString l_daycount = bdaycount;


	_bstr_t bisratefixedphase2(pIsRateFixedPhase2);
	CCString l_isratefixedphase2 = bisratefixedphase2;

	long isratefixedphase2Id;
	if ((isratefixedphase2Id = ARM_ConvYesOrNo (l_isratefixedphase2, C_result)) == ARM_DEFAULT_ERR)
	{
	   ERROR_MSG("Invalid IsRateFixedPhase2",pRet);
	   return S_OK;
	}


	_bstr_t bamort(pAmort);
	CCString l_amort = bamort;

	_bstr_t bbsmod(pBsmod);
	CCString l_bsmod = bbsmod;

	_bstr_t bbsmodDeltaCcy1(pBsmodDeltaCcy1);
	CCString l_bsmodDeltaCcy1 = bbsmodDeltaCcy1;

	_bstr_t bbsmodVegaCcy1(pBsmodVegaCcy1);
	CCString l_bsmodVegaCcy1 = bbsmodVegaCcy1;

	_bstr_t bbsmodDeltaCcy2(pBsmodDeltaCcy2);
	CCString l_bsmodDeltaCcy2 = bbsmodDeltaCcy2;
	
	if ( l_bsmodDeltaCcy2 == "DEFAULT")
	{
		l_bsmodDeltaCcy2 = l_bsmod;
	}		
	
	_bstr_t bbsmodVegaCcy2(pBsmodVegaCcy2);
	CCString l_bsmodVegaCcy2 = bbsmodVegaCcy2;

	if ( l_bsmodVegaCcy2 == "DEFAULT")
	{
		l_bsmodVegaCcy2 = l_bsmod;
	}

	_bstr_t bbsmodFxCorrel(pBsmodFxCorrel);
	CCString l_bsmodFxCorrel = bbsmodFxCorrel;

	if ( l_bsmodFxCorrel == "DEFAULT")
	{
		l_bsmodFxCorrel = l_bsmod;
	}

	_bstr_t bfwdrule(pFwdRule);
	CCString l_fwdrule = bfwdrule;

	_bstr_t bintrule(pIntRule);
	CCString l_intrule = bintrule;

	_bstr_t bstubrule(pStubRule);
	CCString l_stubrule = bstubrule;

	_bstr_t bFreqAmort(pFreqAmort);
	CCString l_FreqAmort = bFreqAmort;

	long freqPayIdPhase1;

	if((freqPayIdPhase1 = ARM_ConvFrequency(l_freqpphase1, C_result)) == ARM_DEFAULT_ERR)
	{
		ERROR_MSG("Invalid Frq",pRet,ARM_ERROR_FREQ);
		return S_OK;
	}

	long freqPayIdPhase2;

	if((freqPayIdPhase2 = ARM_ConvFrequency(l_freqpphase2, C_result)) == ARM_DEFAULT_ERR)
	{
		ERROR_MSG("Invalid Frq",pRet,ARM_ERROR_FREQ);
		return S_OK;
	}

	long freqPayIdFund;

	if((freqPayIdFund = ARM_ConvFrequency(l_freqpFund, C_result)) == ARM_DEFAULT_ERR)
	{
		ERROR_MSG("Invalid Frq",pRet,ARM_ERROR_FREQ);
		return S_OK;
	}

	long freqPayIdPhase3;

	if((freqPayIdPhase3 = ARM_ConvFrequency(l_freqpphase3, C_result)) == ARM_DEFAULT_ERR)
	{
		ERROR_MSG("Invalid SF Frq",pRet,ARM_ERROR_FREQ);
		return S_OK;
	}


	long indexIdPhase1 = ARM_ConvIrType(l_indexphase1);
	long indexIdPhase2 = ARM_ConvIrType(l_indexphase2);
	long indexIdFund = ARM_ConvIrType(l_indexFund);
	long indexIdPhase3 = ARM_ConvIrType(l_indexphase3);
	
	long resetFreqId;

	if((resetFreqId = ARM_ConvFrequency(l_freqr, C_result)) == ARM_DEFAULT_ERR)
	{
		ERROR_MSG("Invalid Frq",pRet,ARM_ERROR_FREQ);
		return S_OK;
	}

	long fwdRuleId;

	if((fwdRuleId = ARM_ConvFwdRule(l_fwdrule, C_result)) == ARM_DEFAULT_ERR)
	{
		ERROR_MSG("Invalid Fwd Rule",pRet,ARM_ERROR_FREQ);
		return S_OK;
	}


	long daycountId= ARM_ConvDayCount(l_daycount);


	long intruleId= ARM_ConvIntRule(l_intrule);

	long resettimingIdPhase1 = ARM_ConvPayResetRule (l_resettimingphase1);
	long resettimingIdPhase2 = ARM_ConvPayResetRule (l_resettimingphase2);
	long resettimingIdPhase3 = ARM_ConvPayResetRule (l_resettimingphase3);

	long stubruleId = ARM_ConvStubRule(l_stubrule);

	//*************************************************************
	//**********    Create Armort Refvalue ID   **********
	//*************************************************************
	long AmortRefvalue;
	long AmortRefvalueToDestroy = ARM_NULL_OBJECT;

	if (LocalGetNumObjectId(l_amort) == ARM_KO)
	{
		double tauxAmort = pTxAmort;
		double amountAmort = pAmountAmort;

		long amortFreqId;

		if((amortFreqId = ARM_ConvFrequency(l_FreqAmort, C_result)) == ARM_DEFAULT_ERR)
		{
			ERROR_MSG("Invalid Amort Frq",pRet,ARM_ERROR_FREQ);
			return S_OK;
		}

		// AmortMethodId
		long amortMethodId;
		if((amortMethodId = ARM_ConvAmortMethod (l_amort, C_result)) == ARM_DEFAULT_ERR)
		{
			ERROR_MSG("Invalid AmortMethod",pRet,ARM_ERROR_FREQ);
			return S_OK;
		}

		// la enddate est la SFDate pour un "une phase" 
		long FixedlegIdamort;
		if (ARMLOCAL_FIXEDLEG(startdate,
							  enddate,
							  K_RCV,
							  0,
							  1.0,
							  daycountId,
							  amortFreqId,
							  K_COMP_PROP,
							  K_ARREARS,
							  intruleId,
							  stubruleId,
							  false,
							  l_ccy,
							  "NULL",
							  K_NX_NONE,
							  -1.0,
							  K_NO,
						  -1,
							  C_result) == ARM_OK) 
							  
		{
			FixedlegIdamort = C_result.getLong();
		}
		else
		{
			ERROR_MSG(C_result.getMsg(),pRet);
			return S_OK;
		}

		if (ARMLOCAL_ARM_GenAmortization(FixedlegIdamort,
										 amortMethodId,
										 amortFreqId,
										 amountAmort,
										 -1,
										 Ntl,
										 tauxAmort,
										 0.,
										 ARM_NULL_OBJECT,
										 0.0,
										 C_result) == ARM_OK) 
		{
			AmortRefvalue = C_result.getLong();
			AmortRefvalueToDestroy = AmortRefvalue;
		}
		else
		{
			ERROR_MSG(C_result.getMsg(),pRet);
			return S_OK;
		}

		#ifdef _DEBUG
			Res = ARMLOCAL_ARM_View(FixedlegIdamort,C_result);
		#endif

		#ifdef _DEBUG
			Res = ARMLOCAL_ARM_View(AmortRefvalue,C_result);
		#endif
	}
	else
	{
		AmortRefvalue = LocalGetNumObjectId(l_amort);
	
		// calcul du 1er nominal
		double dPaymentDate;

		if (ARMLOCAL_DisplayRefValue(AmortRefvalue,true,C_result) == ARM_OK)
		{
			dPaymentDate = C_result.getArray(0);
		}
		else
		{
			ERROR_MSG(C_result.getMsg(),pRet);
			return S_OK;
		}
		
		if (ARMLOCAL_CptRefValue(AmortRefvalue,dPaymentDate,C_result) == ARM_OK)
		{
			Ntl = C_result.getDouble();
		}
		else
		{
			ERROR_MSG(C_result.getMsg(),pRet);
			return S_OK;
		} 
	}
	//************************************************************
	//**********   END Create Amort table    ************
	//************************************************************



	//**********************************************************
	//**********     funding leg & Snd Phase fund    ***********
	//*********      funding en dur en ADV et gap -2  **********
	//**********************************************************
	
	// ici la PayFreq est dite egale a la ResetFreq
	long irindexIdfund;
	if (ARMLOCAL_IRINDEX(KACTUAL_360,
						 freqPayIdFund,
						 -1,
						 K_COMP_PROP,
						 K_MOD_FOLLOWING,
						 K_ADVANCE,
						 -2,
						 K_ARREARS,
						 10000.,
						 false,
						 l_ccy,
						 indexIdFund,
						 K_COMP_PROP,
						 intruleId,
						 freqPayIdFund,
						 C_result) == ARM_OK) 
	{
		irindexIdfund =   C_result.getLong();
	}
	else
	{
		ERROR_MSG(C_result.getMsg(),pRet);
		return S_OK;
	}

	#ifdef _DEBUG
		Res = ARMLOCAL_ARM_View(irindexIdfund,C_result);
	#endif

	long swaplegId;
	if (ARMLOCAL_SWAPLEG(irindexIdfund,
						 startdatephase2,
						 startdatephase3,
						 K_RCV,
						 1,
						 LocalGetNumObjectId(spreadphase2fund),
						 false,
						 l_ccy,
						 daycountId,
						 -2,
						 "NULL",
						 "NULL",
						 0,
						 K_NX_NONE,
						 stubruleId,
						 refdate,
						 K_NO,
						  -1,
						 C_result) == ARM_OK) 
	{
		swaplegId = C_result.getLong();
	}
	else
	{
		ERROR_MSG(C_result.getMsg(),pRet);
		return S_OK;
	}


	
	// ********** Funding Leg ************
	long cloneId;
	if (ARMLOCAL_ClonedAndSetNotional(swaplegId,
									  AmortRefvalue,
									  100.,
									  C_result) == ARM_OK)
	{
		cloneId = C_result.getLong();
	}
	else
	{
		ERROR_MSG(C_result.getMsg(),pRet);
		return S_OK;
	}

	//*********   Price funding Leg   *************************
	
	if (ARMLOCAL_ARM_Price(cloneId,
						   LocalGetNumObjectId(l_bsmod),
						   C_result) == ARM_OK)
	{
		pxFundingLeg = C_result.getDouble();
	}
	else
	{
		ERROR_MSG(C_result.getMsg(),pRet);
		return S_OK;
	}

	#ifdef _DEBUG
		Res = ARMLOCAL_ARM_View(LocalGetNumObjectId(l_bsmod),C_result);
	#endif

	#ifdef _DEBUG
		Res = ARMLOCAL_ARM_View(cloneId,C_result);
	#endif

	#ifdef _DEBUG
		Res = ARMLOCAL_ARM_View(irindexIdfund,C_result);
	#endif

	#ifdef _DEBUG
		Res = ARMLOCAL_ARM_View(swaplegId,C_result);
	#endif

	// shift seulement sur l'euro !?
	if (ARMLOCAL_ARM_Price(cloneId,
						   LocalGetNumObjectId(l_bsmodDeltaCcy1),
						   C_result) == ARM_OK)
	{
		pxFundingLegDelta = C_result.getDouble();
	}
	else
	{
		ERROR_MSG(C_result.getMsg(),pRet);
		return S_OK;
	}

	#ifdef _DEBUG
		Res = ARMLOCAL_ARM_View(LocalGetNumObjectId(l_bsmodDeltaCcy1),C_result);
	#endif

	if ( startdate != startdatephase2 )
	{
		isphase1=true;

	// ********** First Phase Funding Leg ************
		long swaplegIdPhase1fund;
		if (ARMLOCAL_SWAPLEG(irindexIdfund,
							 startdate ,
							 startdatephase2,
							 K_RCV,
							 1,
							 LocalGetNumObjectId(spdphase1fund),
							 false,
							 l_ccy,
							 daycountId,
							 -2,
							 "NULL",
							 "NULL",
							 0,
							 K_NX_NONE,
							 stubruleId,
							 refdate,
							 K_NO,
								-1,
							 C_result) == ARM_OK) 
		{
			swaplegIdPhase1fund = C_result.getLong();
		}
		else
		{
			ERROR_MSG(C_result.getMsg(),pRet);
			return S_OK;
		}

		long cloneIdPhase1fund;
		if (ARMLOCAL_ClonedAndSetNotional(swaplegIdPhase1fund,
										  AmortRefvalue,
										  100.,
										  C_result) == ARM_OK)
		{
			cloneIdPhase1fund = C_result.getLong();
		}
		else
		{
			ERROR_MSG(C_result.getMsg(),pRet);
			return S_OK;
		}

		//*********   Price First Phase Funding Leg   *************************
		if (ARMLOCAL_ARM_Price(cloneIdPhase1fund,
							   LocalGetNumObjectId(l_bsmod),
							   C_result) == ARM_OK)
		{
			pxPhase1fundLeg = C_result.getDouble();
		}
		else
		{
			ERROR_MSG(C_result.getMsg(),pRet);
			return S_OK;
		}

		if (ARMLOCAL_ARM_Price(cloneIdPhase1fund,
							   LocalGetNumObjectId(l_bsmodDeltaCcy1),
							   C_result) == ARM_OK)
		{
			pxPhase1fundLegDelta = C_result.getDouble();
		}
		else
		{
			ERROR_MSG(C_result.getMsg(),pRet);
			return S_OK;
		}

		retCode = ARMLOCAL_FreeObject(swaplegIdPhase1fund,C_result);
		retCode = ARMLOCAL_FreeObject(cloneIdPhase1fund,C_result);

		// ********** ********** ************
		// ******** First Phase Leg **********
		// ********** ********** ************

		long irindexIdPhase1leg;
		
		if (ARMLOCAL_IRINDEX(KACTUAL_360,
							 freqPayIdPhase1 ,
							 -1,
							 K_COMP_PROP,
							 K_MOD_FOLLOWING,
							 resettimingIdPhase1,
							 resetgapphase1,
							 K_ARREARS,
							 10000.,
							 false,
							 l_ccy,
							 indexIdPhase1,
							 K_COMP_PROP,
							 intruleId,
							 freqPayIdPhase1,
							 C_result) == ARM_OK) 
		{
			irindexIdPhase1leg = C_result.getLong();
		}
		else
		{
			ERROR_MSG(C_result.getMsg(),pRet);
			return S_OK;
		}

		long swaplegIdPhase1;
		if (ARMLOCAL_SWAPLEG(irindexIdPhase1leg,
							 startdate ,
							 startdatephase2,
							 K_RCV,
							 1,
							 LocalGetNumObjectId(spdphase1),
							 false,
							 l_ccy,
							 daycountId,
							 resetgapphase1,
							 "NULL",
							 "NULL",
							 0,
							 K_NX_NONE,
							 stubruleId,
							 refdate,
							 K_NO,
						  -1,
							 C_result) == ARM_OK) 
		{
			swaplegIdPhase1 = C_result.getLong();
		}
		else
		{
			ERROR_MSG(C_result.getMsg(),pRet);
			return S_OK;
		}

		long cloneIdPhase1;
		if (ARMLOCAL_ClonedAndSetNotional(swaplegIdPhase1,
										  AmortRefvalue,
										  100.,
										  C_result) == ARM_OK)
		{
			cloneIdPhase1 = C_result.getLong();
		}
		else
		{
			ERROR_MSG(C_result.getMsg(),pRet);
			return S_OK;
		}

		//*********   Price Second Phase Leg   *************************
		if (ARMLOCAL_ARM_Price(cloneIdPhase1,
							   LocalGetNumObjectId(l_bsmod),
							   C_result) == ARM_OK)
		{
			pxPhase1Leg = C_result.getDouble();
		}
		else
		{
			ERROR_MSG(C_result.getMsg(),pRet);
			return S_OK;
		}

		if (ARMLOCAL_ARM_Price(cloneIdPhase1,
							   LocalGetNumObjectId(l_bsmodDeltaCcy1),
							   C_result) == ARM_OK)
		{
			pxPhase1LegDelta = C_result.getDouble();
		}
		else
		{
			ERROR_MSG(C_result.getMsg(),pRet);
			return S_OK;
		}

		retCode = ARMLOCAL_FreeObject(irindexIdPhase1leg,C_result);
		retCode = ARMLOCAL_FreeObject(swaplegIdPhase1,C_result);
		retCode = ARMLOCAL_FreeObject(cloneIdPhase1,C_result);
	}
	else
	{
		pxPhase1fundLeg=0;
		pxPhase1fundLegDelta=0;
		pxPhase1Leg=0;
		pxPhase1LegDelta=0;
	}


	if ( startdatephase3 != enddate )
	{

		isphase3=true;

		// ********** Third Phase Funding Leg ************
		long swaplegIdPhase3fund;
		if (ARMLOCAL_SWAPLEG(irindexIdfund,
							 startdatephase3 ,
							 enddate,
							 K_RCV,
							 1,
							 LocalGetNumObjectId(spdphase3fund),
							 false,
							 l_ccy,
							 daycountId,
							 -2,
							 "NULL",
							 "NULL",
							 0,
							 K_NX_NONE,
							 stubruleId,
							 refdate,
							 K_NO,
						  -1,
							 C_result) == ARM_OK) 
		{
			swaplegIdPhase3fund = C_result.getLong();
		}
		else
		{
			ERROR_MSG(C_result.getMsg(),pRet);
			return S_OK;
		}

		long cloneIdPhase3fund;
		if (ARMLOCAL_ClonedAndSetNotional(swaplegIdPhase3fund,
										  AmortRefvalue,
										  100.,
										  C_result) == ARM_OK)
		{
			cloneIdPhase3fund = C_result.getLong();
		}
		else
		{
			ERROR_MSG(C_result.getMsg(),pRet);
			return S_OK;
		}

		//*********   Price Third Phase Funding Leg   *************************
		if (ARMLOCAL_ARM_Price(cloneIdPhase3fund,
							   LocalGetNumObjectId(l_bsmod),
							   C_result) == ARM_OK)
		{
			pxPhase3fundLeg = C_result.getDouble();
		}
		else
		{
			ERROR_MSG(C_result.getMsg(),pRet);
			return S_OK;
		}

		if (ARMLOCAL_ARM_Price(cloneIdPhase3fund,
							   LocalGetNumObjectId(l_bsmodDeltaCcy1),
							   C_result) == ARM_OK)
		{
			pxPhase3fundLegDelta = C_result.getDouble();
		}
		else
		{
			ERROR_MSG(C_result.getMsg(),pRet);
			return S_OK;
		}

		retCode = ARMLOCAL_FreeObject(irindexIdfund,C_result);
		retCode = ARMLOCAL_FreeObject(swaplegIdPhase3fund,C_result);
		retCode = ARMLOCAL_FreeObject(cloneIdPhase3fund,C_result);

		// ********** ********** ************
		// ******** Third Phase Leg **********
		// ********** ********** ************

		long irindexIdPhase3leg;
		
		if (ARMLOCAL_IRINDEX(KACTUAL_360,
							 freqPayIdPhase3 ,
							 -1,
							 K_COMP_PROP,
							 K_MOD_FOLLOWING,
							 resettimingIdPhase3,
							 resetgapphase3,
							 K_ARREARS,
							 10000.,
							 false,
							 l_ccy,
							 indexIdPhase3,
							 K_COMP_PROP,
							 intruleId,
							 freqPayIdPhase3,
							 C_result) == ARM_OK) 
		{
			irindexIdPhase3leg = C_result.getLong();
		}
		else
		{
			ERROR_MSG(C_result.getMsg(),pRet);
			return S_OK;
		}

		long swaplegIdPhase3;
		if (ARMLOCAL_SWAPLEG(irindexIdPhase3leg,
							 startdatephase3 ,
							 enddate,
							 K_RCV,
							 1,
							 LocalGetNumObjectId(spdphase3),
							 false,
							 l_ccy,
							 daycountId,
							 resetgapphase3,
							 "NULL",
							 "NULL",
							 0,
							 K_NX_NONE,
							 stubruleId,
							 refdate,
							 K_NO,
						  -1,
							 C_result) == ARM_OK) 
		{
			swaplegIdPhase3 = C_result.getLong();
		}
		else
		{
			ERROR_MSG(C_result.getMsg(),pRet);
			return S_OK;
		}

		long cloneIdPhase3;
		if (ARMLOCAL_ClonedAndSetNotional(swaplegIdPhase3,
										  AmortRefvalue,
										  100.,
										  C_result) == ARM_OK)
		{
			cloneIdPhase3 = C_result.getLong();
		}
		else
		{
			ERROR_MSG(C_result.getMsg(),pRet);
			return S_OK;
		}

		//*********   Price Third Phase Leg   *************************
		if (ARMLOCAL_ARM_Price(cloneIdPhase3,
							   LocalGetNumObjectId(l_bsmod),
							   C_result) == ARM_OK)
		{
			pxPhase3Leg = C_result.getDouble();
		}
		else
		{
			ERROR_MSG(C_result.getMsg(),pRet);
			return S_OK;
		}

		if (ARMLOCAL_ARM_Price(cloneIdPhase3,
							   LocalGetNumObjectId(l_bsmodDeltaCcy1),
							   C_result) == ARM_OK)
		{
			pxPhase3LegDelta = C_result.getDouble();
		}
		else
		{
			ERROR_MSG(C_result.getMsg(),pRet);
			return S_OK;
		}

		retCode = ARMLOCAL_FreeObject(irindexIdPhase3leg,C_result);
		retCode = ARMLOCAL_FreeObject(swaplegIdPhase3,C_result);
		retCode = ARMLOCAL_FreeObject(cloneIdPhase3,C_result);
	}
	else
	{
		pxPhase3fundLeg=0;
		pxPhase3fundLegDelta=0;
		pxPhase3Leg=0;
		pxPhase3LegDelta=0;
	}

//******************** Calcul premiere PV ******************
	int iteration =1;
	double prix ;
	double TauxFixe = 10.01;
	double TauxBoniFinal = 10.01;
	double TauxBonifie = 10.01;  // PUREMENT ARBITRAIRE
	fee = Ntl * fee / 100 ;
//******************** Calcul premiere PV ******************

	SensiFundingLeg = pxFundingLegDelta - pxFundingLeg ;

	pvQTFLeg = computeQTF(AsOf ,startdatephase2,startdatephase3,Ntl , indexIdPhase2 ,l_bsmod ,l_bsmodDeltaCcy1 ,l_bsmodVegaCcy1 ,l_bsmodDeltaCcy2 ,l_bsmodVegaCcy2 ,l_bsmodFxCorrel ,freqPayIdPhase2,resetFreqId,resettimingIdPhase2,l_ccy,l_ccyidx ,daycountId, isratefixedphase2Id ,fixedRatePhase2,TauxFixe ,barrier ,spdphase2tf , resetgapphase2 ,-0.01,0.01,fwdRuleId,intruleId,stubruleId,AmortRefvalue ,refdate);
	
	if (pvQTFLeg.size() == 0)	
	{
		ERROR_MSG("Pb in computeQTF",pRet);
		return S_OK;
	}

	pxQTFLeg = pvQTFLeg[0];
	margeQTFTech = pvQTFLeg[1];
	margeQTFDeltaCcy1 = pvQTFLeg[2];
	margeQTFDeltaCcy2 = pvQTFLeg[3];
	margeQTFVega = pvQTFLeg[4];
	margeQTFFxCorrel = pvQTFLeg[5];

	// les deltas se compensent et on sommes en valeur absolue les deltas et vegas

	margeSensi = fabs(fabs( margeQTFDeltaCcy1 - ( pxFundingLegDelta - pxFundingLeg ) + ( pxPhase1LegDelta - pxPhase1Leg ) - ( pxPhase1fundLegDelta - pxPhase1fundLeg ) + ( pxPhase3LegDelta - pxPhase3Leg ) - ( pxPhase3fundLegDelta - pxPhase3fundLeg )) + fabs(margeQTFDeltaCcy2));

	margeQTFVega = fabs(margeQTFVega) ;
		
	margeSensi += margeQTFVega;

	margeQTFFxCorrel = fabs(margeQTFFxCorrel) ;

	marge = max(margeQTFTech,margeSensi) + margeQTFFxCorrel ;
	
	prixAct = pxQTFLeg - pxFundingLeg + pxPhase1Leg - pxPhase1fundLeg + pxPhase3Leg - pxPhase3fundLeg - marge - fee ;
	
	prix = prixAct ;
	prixAct2 =1;
	
//********************                    ******************
//******************** BOUCLE SUR LE PRICE ******************
	while ((abs(prixAct2) > 0.001) && (iteration < 10 ))
	{
		TauxFixe = TauxBonifie;
		
		//*************************
		// **** Premier recalcul ****
		//****************************		
		SensiFundingLeg = pxFundingLegDelta - pxFundingLeg ;

		pvQTFLeg = computeQTF(AsOf ,startdatephase2,startdatephase3,Ntl , indexIdPhase2 ,l_bsmod ,l_bsmodDeltaCcy1 ,l_bsmodVegaCcy1 ,l_bsmodDeltaCcy2 ,l_bsmodVegaCcy2 ,l_bsmodFxCorrel ,freqPayIdPhase2,resetFreqId,resettimingIdPhase2,l_ccy,l_ccyidx ,daycountId,isratefixedphase2Id ,fixedRatePhase2,TauxBonifie ,barrier ,spdphase2tf , resetgapphase2 ,-0.01,0.01,fwdRuleId,intruleId,stubruleId,AmortRefvalue ,refdate);

		if (pvQTFLeg.size() == 0)
		{
			ERROR_MSG("Pb in computeQTF",pRet);
			return S_OK;
		}

		pxQTFLeg = pvQTFLeg[0];
		margeQTFTech = pvQTFLeg[1];
		margeQTFDeltaCcy1 = pvQTFLeg[2];
		margeQTFDeltaCcy2 = pvQTFLeg[3];
		margeQTFVega = pvQTFLeg[4];
		margeQTFFxCorrel = pvQTFLeg[5];

		// les deltas se compensent et on sommes en valeur absolue les deltas et vegas

		margeSensi = fabs(fabs( margeQTFDeltaCcy1 - ( pxFundingLegDelta - pxFundingLeg ) + ( pxPhase1LegDelta - pxPhase1Leg ) - ( pxPhase1fundLegDelta - pxPhase1fundLeg ) + ( pxPhase3LegDelta - pxPhase3Leg ) - ( pxPhase3fundLegDelta - pxPhase3fundLeg )) + fabs(margeQTFDeltaCcy2));
		margeQTFVega = fabs(margeQTFVega) ;
		margeSensi += margeQTFVega;
		margeQTFFxCorrel = fabs(margeQTFFxCorrel) ;
		marge = max(margeQTFTech,margeSensi) + margeQTFFxCorrel ;

		prixAct = pxQTFLeg - pxFundingLeg + pxPhase1Leg - pxPhase1fundLeg + pxPhase3Leg - pxPhase3fundLeg - marge- fee ;

	
		//******************************
		// **** affectation valeurs ****
		//******************************
		prixAct2 = prixAct ;
		Xincr = TauxFixe / 1000000;
		TauxBonifie = TauxFixe + Xincr;		

		//****************************
		// **** Deuxieme recalcul ****
		//****************************
		SensiFundingLeg = pxFundingLegDelta - pxFundingLeg ;

		pvQTFLeg = computeQTF(AsOf ,startdatephase2,startdatephase3,Ntl , indexIdPhase2 ,l_bsmod ,l_bsmodDeltaCcy1 ,l_bsmodVegaCcy1 ,l_bsmodDeltaCcy2 ,l_bsmodVegaCcy2 ,l_bsmodFxCorrel ,freqPayIdPhase2,resetFreqId,resettimingIdPhase2,l_ccy,l_ccyidx ,daycountId,isratefixedphase2Id ,fixedRatePhase2,TauxBonifie ,barrier ,spdphase2tf , resetgapphase2 ,-0.01,0.01,fwdRuleId,intruleId,stubruleId,AmortRefvalue ,refdate);

		if (pvQTFLeg.size() == 0)
		{
			ERROR_MSG("Pb in computeQTF",pRet);
			return S_OK;
		}

		pxQTFLeg = pvQTFLeg[0];
		margeQTFTech = pvQTFLeg[1];
		margeQTFDeltaCcy1 = pvQTFLeg[2];
		margeQTFDeltaCcy2 = pvQTFLeg[3];
		margeQTFVega = pvQTFLeg[4];
		margeQTFFxCorrel = pvQTFLeg[5];

		// les deltas se compensent et on sommes en valeur absolue les deltas et vegas

		margeSensi = fabs(fabs( margeQTFDeltaCcy1 - ( pxFundingLegDelta - pxFundingLeg ) + ( pxPhase1LegDelta - pxPhase1Leg ) - ( pxPhase1fundLegDelta - pxPhase1fundLeg ) + ( pxPhase3LegDelta - pxPhase3Leg ) - ( pxPhase3fundLegDelta - pxPhase3fundLeg )) + fabs(margeQTFDeltaCcy2));
		margeQTFVega = fabs(margeQTFVega) ;
		margeSensi += margeQTFVega;
		margeQTFFxCorrel = fabs(margeQTFFxCorrel) ;
		marge = max(margeQTFTech,margeSensi) + margeQTFFxCorrel ;
		prixAct = pxQTFLeg - pxFundingLeg + pxPhase1Leg - pxPhase1fundLeg + pxPhase3Leg - pxPhase3fundLeg - marge- fee ;
	
		iteration++;

		TauxBoniFinal = TauxBonifie;
		prixInter = ( prixAct - prixAct2 ) / Xincr ;
		TauxBonifie = TauxFixe - prixAct / prixInter ; 
		prixAct2 = prixAct ;		
	}

	//**************** PRICE *******************
	retCode = ARMLOCAL_FreeObject(swaplegId,C_result);
	retCode = ARMLOCAL_FreeObject(cloneId,C_result);
	retCode = ARMLOCAL_FreeObject(AmortRefvalueToDestroy,C_result);

	VECTOR<double> vRes;

	vRes.push_back(TauxBoniFinal);
	vRes.push_back(margeQTFTech);
	vRes.push_back(margeSensi);
	vRes.push_back(margeQTFFxCorrel);

	Res = VECTORDOUBLE2VARIANT(vRes, pRet);

	_variant_t wrap_Res;
	wrap_Res.Attach(*pRet);
	*pRet=wrap_Res.Detach();

	return S_OK;
}
catch(...)
{
	return E_FAIL;
}
}



VECTOR<double> computeBilibor(long irindexIdPF2leg, double startdate, double SFDate, double SprSeekBonifie, CCString pccy, long daycountId, double resetgap, long stubruleId, 
						double prefdate, long pAmortRefvalue, CCString pbsmod, CCString pbsmodDeltaCcy1, CCString pbsmodDeltaCcy2, CCString pbsmodFxCorrel)
;

//Ajout TD le 04/07/2005
STDMETHODIMP ActiveXModule::ARMcomputeBilibor(double pAsOf,
								   double pStartDate,
								   double pDateSecondPhase,
								   double pEndDate,
								   double pNtl,
								   BSTR pIndex,
								   BSTR pIndexFund,
								   BSTR pIndexSF,
								   BSTR pBsmod,
								   BSTR pBsmodFund,
								   BSTR pBsmodDeltaCcy1,
								   BSTR pBsmodDeltaFund,
								   BSTR pBsmodDeltaCcy2,
								   BSTR pBsmodFxCorrel,
								   BSTR pFreqP,
								   BSTR pFreqR,
								   BSTR pFreqPFund,
								   BSTR pFreqRFund, 
								   BSTR pFreqPSF,
								   BSTR pFreqRSF,
								   BSTR pResetTiming,
								   BSTR pResetTimingSF,
								   BSTR pCcy1,
								   BSTR pCcy2,
								   BSTR pDayCount,
								   BSTR pDayCountSF,
								   double pSpdPF,
								   double pSpdSF,
								   double pSpdfund,
								   double pSpdfund2, 
								   double pResetGap,
								   double pResetGapSF,
								   BSTR pAmort,
								   double pRefDate,
								   double pFee,
								   BSTR pFwdRule,
								   BSTR pIntRule,
								   BSTR pStubRule,
								   BSTR pFreqAmort,
								   double pTxAmort,
								   double pAmountAmort,
								   VARIANT *pRet)

{
	// TODO: Add your implementation code here
try
{
	ARM_result C_result;
	 
	double AsOf = pAsOf;
	double startdate = pStartDate;
	double SFDate = pDateSecondPhase;
	double enddate = pEndDate;
	double Ntl = pNtl;
	double resetgap = pResetGap;
	double resetgapSF = pResetGapSF;
	double refdate = pRefDate;
	double fee = pFee;

	double pxFundingLeg;
	double pxFundingLegDelta;

	double pxFundingLeg2;
	double pxFundingLegDelta2;

	double pxSFLeg;
	double pxSFLegDelta;

	double Xincr;
	double prixAct;
	double prixAct2;
	double prixInter;

	double pxPF2Leg;
	double pxPF2LegDeltaCcy1;
	double pxPF2LegDeltaCcy2;
	double pxPF2LegFxCorrel;

	double pxPF1Leg;
	double pxPF1LegDeltaCcy1;
	double pxPF1LegDeltaCcy2;
	double pxPF1LegFxCorrel;

	double pxTotal;
	double pxDeltaCcy1;
	double pxDeltaCcy2;
	double pxFxCorrel; 
	double SensiTotal;

	double pxSFLegLisse;
	double pxSFLegDeltaLisse;

	VECTOR<double> pvBiliborLeg;

	double spread = pSpdPF;
	double spdSF = pSpdSF;
	double spdSFfund = pSpdfund;
	double spdSFfund2 = pSpdfund2;

	double SprSeekFinalLisse; 

	long retCode;
	double Res = 0.;

	_bstr_t bindex(pIndex);
	CCString l_index = bindex;

	_bstr_t bindexFund(pIndexFund);
	CCString l_indexFund = bindexFund;

	_bstr_t bindexSF(pIndexSF);
	CCString l_indexSF = bindexSF;

	_bstr_t bfreqp(pFreqP);
	CCString l_freqp = bfreqp;

	_bstr_t bfreqpFund(pFreqPFund);
	CCString l_freqpFund = bfreqpFund;

	_bstr_t bfreqrFund(pFreqRFund);
	CCString l_freqrFund = bfreqrFund;

	_bstr_t bfreqr(pFreqR);
	CCString l_freqr = bfreqr;
	
	_bstr_t bfreqpSF(pFreqPSF);
	CCString l_freqpSF = bfreqpSF;

	_bstr_t bfreqrSF(pFreqRSF);
	CCString l_freqrSF = bfreqrSF;

	_bstr_t bresettiming(pResetTiming);
	CCString l_resettiming = bresettiming;

	_bstr_t bresettimingSF(pResetTimingSF);
	CCString l_resettimingSF = bresettimingSF;

	_bstr_t bccy(pCcy1);
	CCString l_ccy = bccy;

	_bstr_t bccy2(pCcy2);
	CCString l_ccy2 = bccy2;

	_bstr_t bdaycount(pDayCount);
	CCString l_daycount = bdaycount;

	_bstr_t bdaycountSF(pDayCountSF);
	CCString l_daycountSF = bdaycountSF;

	_bstr_t bamort(pAmort);
	CCString l_amort = bamort;

	_bstr_t bbsmod(pBsmod);
	CCString l_bsmod = bbsmod;

	_bstr_t bbsmodFund(pBsmodFund);
	CCString l_bsmodFund = bbsmodFund;

	if ( l_bsmodFund == "DEFAULT")
	{
		l_bsmodFund = l_bsmod;
	}

	_bstr_t bbsmodDeltaCcy1(pBsmodDeltaCcy1);
	CCString l_bsmodDeltaCcy1 = bbsmodDeltaCcy1;
	if ( l_bsmodDeltaCcy1 == "DEFAULT")
	{
		l_bsmodDeltaCcy1 = l_bsmod;
	}

	_bstr_t bbsmodDeltaCcy2(pBsmodDeltaCcy2);
	CCString l_bsmodDeltaCcy2 = bbsmodDeltaCcy2;
	
	if ( l_bsmodDeltaCcy2 == "DEFAULT")
	{
		l_bsmodDeltaCcy2 = l_bsmod;
	}
	
	_bstr_t bbsmodDeltaFund(pBsmodDeltaFund);
	CCString l_bsmodDeltaFund = bbsmodDeltaFund;
	
	if ( l_bsmodDeltaFund == "DEFAULT")
	{
		l_bsmodDeltaFund = l_bsmod;
	}

	_bstr_t bbsmodFxCorrel(pBsmodFxCorrel);
	CCString l_bsmodFxCorrel = bbsmodFxCorrel;

	if ( l_bsmodFxCorrel == "DEFAULT")
	{
		l_bsmodFxCorrel = l_bsmod;
	}

	_bstr_t bfwdrule(pFwdRule);
	CCString l_fwdrule = bfwdrule;

	_bstr_t bintrule(pIntRule);
	CCString l_intrule = bintrule;

	_bstr_t bstubrule(pStubRule);
	CCString l_stubrule = bstubrule;

	long freqPayId;
	if((freqPayId = ARM_ConvFrequency(l_freqp, C_result)) == ARM_DEFAULT_ERR)
	{
		ERROR_MSG("Invalid Frq",pRet,ARM_ERROR_FREQ);
		return S_OK;
	}

	long resetFreqId;
	if((resetFreqId = ARM_ConvFrequency(l_freqr, C_result)) == ARM_DEFAULT_ERR)
	{
		ERROR_MSG("Invalid Frq",pRet,ARM_ERROR_FREQ);
		return S_OK;
	}

	long freqPayIdSF;
	if((freqPayIdSF = ARM_ConvFrequency(l_freqpSF, C_result)) == ARM_DEFAULT_ERR)
	{
		ERROR_MSG("Invalid Frq",pRet,ARM_ERROR_FREQ);
		return S_OK;
	}

	long resetFreqIdSF;
	if((resetFreqIdSF = ARM_ConvFrequency(l_freqrSF, C_result)) == ARM_DEFAULT_ERR)
	{
		ERROR_MSG("Invalid Frq",pRet,ARM_ERROR_FREQ);
		return S_OK;
	}

	long freqPayIdFund;
	if((freqPayIdFund = ARM_ConvFrequency(l_freqpFund, C_result)) == ARM_DEFAULT_ERR)
	{
		ERROR_MSG("Invalid Frq",pRet,ARM_ERROR_FREQ);
		return S_OK;
	}

	long resetFreqIdFund;
	if((resetFreqIdFund = ARM_ConvFrequency(l_freqrFund, C_result)) == ARM_DEFAULT_ERR)
	{
		ERROR_MSG("Invalid Frq",pRet,ARM_ERROR_FREQ);
		return S_OK;
	}

	long indexId = ARM_ConvIrType(l_index);

	long indexIdFund = ARM_ConvIrType(l_indexFund);

	long indexIdSF = ARM_ConvIrType(l_indexSF);
	
	long fwdRuleId;
	if((fwdRuleId = ARM_ConvFwdRule(l_fwdrule, C_result)) == ARM_DEFAULT_ERR)
	{
		ERROR_MSG("Invalid Fwd Rule",pRet,ARM_ERROR_FREQ);
		return S_OK;
	}

	long daycountId= ARM_ConvDayCount(l_daycount);
	long daycountIdSF= ARM_ConvDayCount(l_daycountSF);
	long intruleId= ARM_ConvIntRule(l_intrule);
	long resettimingId = ARM_ConvPayResetRule (l_resettiming);
	long resettimingIdSF = ARM_ConvPayResetRule (l_resettimingSF);
	long stubruleId = ARM_ConvStubRule(l_stubrule);

	//*************************************************************
	//**********    Create Armort Refvalue ID   **********
	//*************************************************************
	long AmortRefvalue;
	long AmortRefvalueToDestroy = ARM_NULL_OBJECT;

	if (LocalGetNumObjectId(l_amort) == ARM_KO)
	{
		double tauxAmort = pTxAmort;
		double amountAmort = pAmountAmort;

		_bstr_t bFreqAmort(pFreqAmort);
		CCString l_FreqAmort = bFreqAmort;

		long amortFreqId;
		if((amortFreqId = ARM_ConvFrequency(l_FreqAmort, C_result)) == ARM_DEFAULT_ERR)
		{
			ERROR_MSG("Invalid Amort Frq",pRet,ARM_ERROR_FREQ);
			return S_OK;
		}

		// AmortMethodId
		long amortMethodId;
		if((amortMethodId = ARM_ConvAmortMethod (l_amort, C_result)) == ARM_DEFAULT_ERR)
		{
			ERROR_MSG("Invalid AmortMethod",pRet,ARM_ERROR_FREQ);
			return S_OK;
		}

		// la enddate est la SFDate pour un "une phase" 
		long FixedlegIdamort;
		if (ARMLOCAL_FIXEDLEG(startdate,
							  enddate,
							  K_RCV,
							  0,
							  1.0,
							  daycountId,
							  amortFreqId,
							  K_COMP_PROP,
							  K_ARREARS,
							  intruleId,
							  stubruleId,
							  false,
							  l_ccy,
							  "NULL",
							  K_NX_NONE,
							  -1.0,
							  K_NO,
						  -1,
							  C_result) == ARM_OK) 
							  
		{
			FixedlegIdamort = C_result.getLong();
		}
		else
		{
			ERROR_MSG(C_result.getMsg(),pRet);
			return S_OK;
		}

		if (ARMLOCAL_ARM_GenAmortization(FixedlegIdamort,
										 amortMethodId,
										 amortFreqId,
										 amountAmort,
										 -1,
										 Ntl,
										 tauxAmort,
										 0.,
										 ARM_NULL_OBJECT,
										 0.0,
										 C_result) == ARM_OK) 
		{
			AmortRefvalue = C_result.getLong();
			AmortRefvalueToDestroy = AmortRefvalue;
		}
		else
		{
			ERROR_MSG(C_result.getMsg(),pRet);
			return S_OK;
		}
	}
	else
	{
		AmortRefvalue = LocalGetNumObjectId(l_amort);

		// calcul du 1er nominal
		double dPaymentDate;

		if (ARMLOCAL_DisplayRefValue(AmortRefvalue,true,C_result) == ARM_OK)
		{
			dPaymentDate = C_result.getArray(0);
		}
		else
		{
			ERROR_MSG(C_result.getMsg(),pRet);
			return S_OK;
		}
		
		if (ARMLOCAL_CptRefValue(AmortRefvalue,dPaymentDate,C_result) == ARM_OK)
		{
			Ntl = C_result.getDouble();
		}
		else
		{
			ERROR_MSG(C_result.getMsg(),pRet);
			return S_OK;
		} 
	}

	//************************************************************
	//**********   END Create Amort table    ************
	//************************************************************


	//**********************************************************
	//**********     funding leg & Snd Phase fund    ***********
	//*********      funding en dur en ADV et gap -2  **********
	//**********************************************************
	long irindexIdfund;
	if (ARMLOCAL_IRINDEX(KACTUAL_360,
						 freqPayIdFund,
						 -1,
						 K_COMP_PROP,
						 K_MOD_FOLLOWING,
						 K_ADVANCE,
						 -2,
						 K_ARREARS,
						 10000.,
						 false,
						 l_ccy,
						 indexIdFund,
						 K_COMP_PROP,
						 intruleId,
						 resetFreqIdFund,
						 C_result) == ARM_OK) 
	{
		irindexIdfund = C_result.getLong();
	}
	else
	{
		ERROR_MSG(C_result.getMsg(),pRet);
		return S_OK;
	}
	// la enddate est la SFDate pour un "une phase" 
	long swaplegId;
	if (ARMLOCAL_SWAPLEG(irindexIdfund,
						 startdate,
						 SFDate,
						 K_RCV,
						 /*1,
						 LocalGetNumObjectId(spdSFfund),*/
						 0, 
						 spdSFfund, 
						 false,
						 l_ccy,
						 daycountId,
						 -2,
						 "NULL",
						 "NULL",
						 0,
						 K_NX_NONE,
						 stubruleId,
						 refdate,
						 K_NO,
						  -1,
						 C_result) == ARM_OK) 
	{
		swaplegId = C_result.getLong();
	}
	else
	{
		ERROR_MSG(C_result.getMsg(),pRet);
		return S_OK;
	}

	// ********** Funding Leg ************
	long cloneId;
	if (ARMLOCAL_ClonedAndSetNotional(swaplegId,
									  AmortRefvalue,
									  100.,
									  C_result) == ARM_OK)
	{
		cloneId = C_result.getLong();
	}
	else
	{
		ERROR_MSG(C_result.getMsg(),pRet);
		return S_OK;
	}

	//*********   Price funding Leg   *************************
	if (ARMLOCAL_ARM_Price(cloneId,
						   LocalGetNumObjectId(l_bsmodFund),
						   C_result) == ARM_OK)
	{
		pxFundingLeg = C_result.getDouble();
	}
	else
	{
		ERROR_MSG(C_result.getMsg(),pRet);
		return S_OK;
	}

	// shift seulement sur l'euro !?
	if (ARMLOCAL_ARM_Price(cloneId,
						   LocalGetNumObjectId(l_bsmodDeltaFund),
						   C_result) == ARM_OK)
	{
		pxFundingLegDelta = C_result.getDouble();
	}
	else
	{
		ERROR_MSG(C_result.getMsg(),pRet);
		return S_OK;
	}

	long swaplegId2;
	long cloneId2;
	pxFundingLeg2= 0.;
	pxFundingLegDelta2 = 0.;
	if ( SFDate != enddate ) 
	{
		if (ARMLOCAL_SWAPLEG(irindexIdfund,
							 SFDate,
							 enddate,
							 K_RCV,
							 /*1,
							 LocalGetNumObjectId(spdSFfund),*/
							 0, 
							 spdSFfund2, 
							 false,
							 l_ccy,
							 daycountId,
							 -2,
							 "NULL",
							 "NULL",
							 0,
							 K_NX_NONE,
							 stubruleId,
							 refdate,
							 K_NO,
						  -1,
							 C_result) == ARM_OK) 
		{
			swaplegId2 = C_result.getLong();
		}
		else
		{
			ERROR_MSG(C_result.getMsg(),pRet);
			return S_OK;
		}

		if (ARMLOCAL_ClonedAndSetNotional(swaplegId2,
										  AmortRefvalue,
										  100.,
										  C_result) == ARM_OK)
		{
			cloneId2 = C_result.getLong();
		}
		else
		{
			ERROR_MSG(C_result.getMsg(),pRet);
			return S_OK;
		}

		if (ARMLOCAL_ARM_Price(cloneId2,
							   LocalGetNumObjectId(l_bsmodFund),
							   C_result) == ARM_OK)
		{
			pxFundingLeg2 = C_result.getDouble();
		}
		else
		{
			ERROR_MSG(C_result.getMsg(),pRet);
			return S_OK;
		}

		if (ARMLOCAL_ARM_Price(cloneId2,
							   LocalGetNumObjectId(l_bsmodDeltaFund),
							   C_result) == ARM_OK)
		{
			pxFundingLegDelta2 = C_result.getDouble();
		}
		else
		{
			ERROR_MSG(C_result.getMsg(),pRet);
			return S_OK;
		}

		retCode = ARMLOCAL_FreeObject(swaplegId2,C_result);
		retCode = ARMLOCAL_FreeObject(cloneId2,C_result);
	}

	pxFundingLeg += pxFundingLeg2; 
	pxFundingLegDelta += pxFundingLegDelta2;

	retCode = ARMLOCAL_FreeObject(swaplegId,C_result);
	retCode = ARMLOCAL_FreeObject(cloneId,C_result);
		
	// Fin du calcul du prix de la funding Leg


//********************Premiere Phase Patte structure 2 TD****************	
	long irindexIdPF2leg;
	
	if (ARMLOCAL_IRINDEX(KACTUAL_360,
						 freqPayId,
						 -1,
						 K_COMP_PROP,
						 K_MOD_FOLLOWING,
						 resettimingId,
						 resetgap,
						 K_ARREARS,
						 10000.,
						 false,
						 l_ccy2,
						 indexId,
						 K_COMP_PROP,
						 intruleId,
						 resetFreqId,
						 C_result) == ARM_OK) 
	{
		irindexIdPF2leg = C_result.getLong();
	}
	else
	{
		ERROR_MSG(C_result.getMsg(),pRet);
		return S_OK;
	}

	long swaplegIdPF2;

	if (ARMLOCAL_SWAPLEG(irindexIdPF2leg,
						 startdate,
						 SFDate,
						 K_RCV,
//						 1,
//						 LocalGetNumObjectId(spread),
						 0,
						 spread,							
						 false,
						 l_ccy,
						 daycountId,
						 resetgap,
						 "NULL",
						 "NULL",
						 0,
						 K_NX_NONE,
						 stubruleId,
						 refdate,
						 K_NO,
						  -1,
						 C_result) == ARM_OK)
	{
		swaplegIdPF2 = C_result.getLong();
	}
	else
	{
		ERROR_MSG(C_result.getMsg(),pRet);
		return S_OK;
	}

	long cloneIdPF2;
	if (ARMLOCAL_ClonedAndSetNotional(swaplegIdPF2,
									  AmortRefvalue,
									  100.,
									  C_result) == ARM_OK)
	{
		cloneIdPF2 = C_result.getLong();
	}
	else
	{
		ERROR_MSG(C_result.getMsg(),pRet);
		return S_OK;
	}

//*   Price Première Phase structure2  TD*
	if (ARMLOCAL_ARM_Price(cloneIdPF2,
						   LocalGetNumObjectId(l_bsmod),
						   C_result) == ARM_OK)
	{
		pxPF2Leg = C_result.getDouble();
	}
	else
	{
		ERROR_MSG(C_result.getMsg(),pRet);
		return S_OK;
	}

	if (ARMLOCAL_ARM_Price(cloneIdPF2,
						   LocalGetNumObjectId(l_bsmodDeltaCcy1),
						   C_result) == ARM_OK)
	{
		pxPF2LegDeltaCcy1 = C_result.getDouble();
	}
	else
	{
		ERROR_MSG(C_result.getMsg(),pRet);
		return S_OK;
	}

	if (ARMLOCAL_ARM_Price(cloneIdPF2,
						   LocalGetNumObjectId(l_bsmodDeltaCcy2),
						   C_result) == ARM_OK)
	{
		pxPF2LegDeltaCcy2 = C_result.getDouble();
	}
	else
	{
		ERROR_MSG(C_result.getMsg(),pRet);
		return S_OK;
	}

	if (ARMLOCAL_ARM_Price(cloneIdPF2,
						   LocalGetNumObjectId(l_bsmodFxCorrel),
						   C_result) == ARM_OK)
	{
		pxPF2LegFxCorrel = C_result.getDouble();
	}
	else
	{
		ERROR_MSG(C_result.getMsg(),pRet);
		return S_OK;
	}

	retCode = ARMLOCAL_FreeObject(swaplegIdPF2,C_result);
	retCode = ARMLOCAL_FreeObject(cloneIdPF2,C_result);

	//********************Deuxième Phase Patte structure TD****************
	long irindexIdSFleg;

	if ( SFDate != enddate )
	{	
		if (ARMLOCAL_IRINDEX(KACTUAL_360,
							 freqPayIdSF,
							 -1,
							 K_COMP_PROP,
							 K_MOD_FOLLOWING,
							 resettimingIdSF,
							 resetgapSF,
							 K_ARREARS,
							 10000.,
							 false,
							 l_ccy,
							 indexIdSF,
							 K_COMP_PROP,
							 intruleId,
							 resetFreqIdSF,
							 C_result) == ARM_OK) 
		{
			irindexIdSFleg = C_result.getLong();
		}
		else
		{
			ERROR_MSG(C_result.getMsg(),pRet);
			return S_OK;
		}

		long swaplegISF;
		// A remplacer par le spread mis en argument
		if (ARMLOCAL_SWAPLEG(irindexIdSFleg,
							 SFDate,
							 enddate,
							 K_RCV,
//							 1,
//							 LocalGetNumObjectId(spdSF),
							 0,
							 spdSF, 
							 false,
							 l_ccy,
							 daycountId,
							 resetgapSF,
							 "NULL",
							 "NULL",
							 0,
							 K_NX_NONE,
							 stubruleId,
							 refdate,
							 K_NO,
						  -1,
							 C_result) == ARM_OK)
		{
			swaplegISF = C_result.getLong();
		}
		else
		{
			ERROR_MSG(C_result.getMsg(),pRet);
			return S_OK;
		}

		long cloneIdSF;
		if (ARMLOCAL_ClonedAndSetNotional(swaplegISF,
										  AmortRefvalue,
										  100.,
										  C_result) == ARM_OK)
		{
			cloneIdSF = C_result.getLong();
		}
		else
		{
			ERROR_MSG(C_result.getMsg(),pRet);
			return S_OK;
		}

	//* Price Deuxième phase Phase structure TD*
		if (ARMLOCAL_ARM_Price(cloneIdSF,
							   LocalGetNumObjectId(l_bsmodFund),
							   C_result) == ARM_OK)
		{
			pxSFLeg = C_result.getDouble();
		}
		else
		{
			ERROR_MSG(C_result.getMsg(),pRet);
			return S_OK;
		}

		if (ARMLOCAL_ARM_Price(cloneIdSF,
							   LocalGetNumObjectId(l_bsmodDeltaFund),
							   C_result) == ARM_OK)
		{
			pxSFLegDelta = C_result.getDouble();
		}
		else
		{
			ERROR_MSG(C_result.getMsg(),pRet);
			return S_OK;
		}

		retCode = ARMLOCAL_FreeObject(swaplegISF,C_result);
		retCode = ARMLOCAL_FreeObject(cloneIdSF,C_result);	
	}
	else
	{
		pxSFLeg = 0.;
		pxSFLegDelta = 0.;
	}

//******************** Calcul premiere PV ******************
	int iteration =1;
	double prix ;
//	double TauxFixe = 10.01;
	double SprSeek = 10.01;
	double SprSeekFinal= 10.01;
	double SprSeekBonifie = 10.01;  // PUREMENT ARBITRAIRE
	fee = Ntl * fee / 100 ;
//******************** Calcul premiere PV ******************

	pvBiliborLeg = computeBilibor(irindexIdPF2leg, startdate, SFDate, SprSeekBonifie, l_ccy, daycountId, resetgap, stubruleId, 
						refdate, AmortRefvalue, l_bsmod, l_bsmodDeltaCcy1, l_bsmodDeltaCcy2, l_bsmodFxCorrel );

	if (pvBiliborLeg.size() == 0)
	{
		ERROR_MSG("Pb in computeBilibor",pRet);
		return S_OK;
	}

	pxPF1Leg = pvBiliborLeg[0];
	pxPF1LegDeltaCcy1 = pvBiliborLeg[1];
	pxPF1LegDeltaCcy2 = pvBiliborLeg[2];
	pxPF1LegFxCorrel = pvBiliborLeg[3];

	pxTotal = pxPF1Leg + pxPF2Leg + pxSFLeg - pxFundingLeg;
	pxDeltaCcy1 = (pxPF1LegDeltaCcy1 - pxPF1Leg) + (pxPF2LegDeltaCcy2 - pxPF2Leg) + (pxSFLegDelta-pxSFLeg) - (pxFundingLegDelta-pxFundingLeg); 
	if (pxDeltaCcy1 <0.) pxDeltaCcy1 = - pxDeltaCcy1;
	pxDeltaCcy2 = (pxPF1LegDeltaCcy2 - pxPF1Leg) + (pxPF2LegDeltaCcy2 - pxPF2Leg);
	if (pxDeltaCcy2 <0.) pxDeltaCcy2 = - pxDeltaCcy2;
	pxFxCorrel = (pxPF1LegFxCorrel - pxPF1Leg) + (pxPF2LegFxCorrel - pxPF2Leg);
	if (pxFxCorrel <0.) pxFxCorrel = - pxFxCorrel;

	SensiTotal = pxDeltaCcy1 + pxDeltaCcy2 + pxFxCorrel;
	prixAct = pxTotal - SensiTotal - fee; 

	prix = prixAct ;
	prixAct2 =1;

//********************                    ******************
//******************** BOUCLE SUR LE PRICE ******************
	while ((abs(prixAct2) > 0.0001) && (iteration < 10 ))
	{
		SprSeek = SprSeekBonifie;
		
		//*************************
		// **** Premier recalcul ****
		//****************************

		pvBiliborLeg = computeBilibor(irindexIdPF2leg, startdate, SFDate, SprSeekBonifie, l_ccy, daycountId, resetgap, stubruleId, 
						refdate, AmortRefvalue, l_bsmod, l_bsmodDeltaCcy1, l_bsmodDeltaCcy2, l_bsmodFxCorrel );

		if (pvBiliborLeg.size() == 0)
		{
			ERROR_MSG("Pb in computeBilibor",pRet);
			return S_OK;
		}

		pxPF1Leg = pvBiliborLeg[0];
		pxPF1LegDeltaCcy1 = pvBiliborLeg[1];
		pxPF1LegDeltaCcy2 = pvBiliborLeg[2];
		pxPF1LegFxCorrel = pvBiliborLeg[3];

		pxTotal = pxPF1Leg + pxPF2Leg + pxSFLeg - pxFundingLeg;
		pxDeltaCcy1 = (pxPF1LegDeltaCcy1 - pxPF1Leg) + (pxPF2LegDeltaCcy1 - pxPF2Leg) + (pxSFLegDelta-pxSFLeg) - (pxFundingLegDelta-pxFundingLeg); 
		if (pxDeltaCcy1 <0.) pxDeltaCcy1 = - pxDeltaCcy1;
		pxDeltaCcy2 = (pxPF1LegDeltaCcy2 - pxPF1Leg) + (pxPF2LegDeltaCcy2 - pxPF2Leg);
		if (pxDeltaCcy2 <0.) pxDeltaCcy2 = - pxDeltaCcy2;
		pxFxCorrel = (pxPF1LegFxCorrel - pxPF1Leg) + (pxPF2LegFxCorrel - pxPF2Leg);
		if (pxFxCorrel <0.) pxFxCorrel = - pxFxCorrel;

		SensiTotal = pxDeltaCcy1 + pxDeltaCcy2 + pxFxCorrel;
		prixAct = pxTotal - SensiTotal - fee; 
		// les deltas se compensent et on somme en valeur absolue les deltas et vegas

		//******************************
		// **** affectation valeurs ****
		//******************************		
		prixAct2 = prixAct ;
		Xincr = SprSeek / 1000000;
		SprSeekBonifie = SprSeek + Xincr;		

		//****************************
		// **** Deuxieme recalcul ****
		//****************************
		pvBiliborLeg = computeBilibor(irindexIdPF2leg, startdate, SFDate, SprSeekBonifie, l_ccy, daycountId, resetgap, stubruleId, 
						refdate, AmortRefvalue, l_bsmod, l_bsmodDeltaCcy1, l_bsmodDeltaCcy2, l_bsmodFxCorrel );

		if (pvBiliborLeg.size() == 0)
		{
			ERROR_MSG("Pb in computeBilibor",pRet);
			return S_OK;
		}

		pxPF1Leg = pvBiliborLeg[0];
		pxPF1LegDeltaCcy1 = pvBiliborLeg[1];
		pxPF1LegDeltaCcy2 = pvBiliborLeg[2];
		pxPF1LegFxCorrel = pvBiliborLeg[3];

		pxTotal = pxPF1Leg + pxPF2Leg + pxSFLeg - pxFundingLeg;
		pxDeltaCcy1 = (pxPF1LegDeltaCcy1 - pxPF1Leg) + (pxPF2LegDeltaCcy1 - pxPF2Leg) + (pxSFLegDelta-pxSFLeg) - (pxFundingLegDelta-pxFundingLeg); 
		if (pxDeltaCcy1 <0.) pxDeltaCcy1 = - pxDeltaCcy1;
		pxDeltaCcy2 = (pxPF1LegDeltaCcy2 - pxPF1Leg) + (pxPF2LegDeltaCcy2 - pxPF2Leg);
		if (pxDeltaCcy2 <0.) pxDeltaCcy2 = - pxDeltaCcy2;
		pxFxCorrel = (pxPF1LegFxCorrel - pxPF1Leg) + (pxPF2LegFxCorrel - pxPF2Leg);
		if (pxFxCorrel <0.) pxFxCorrel = - pxFxCorrel;

		SensiTotal = pxDeltaCcy1 + pxDeltaCcy2 + pxFxCorrel;
		prixAct = pxTotal - SensiTotal - fee; 

		iteration++;

		SprSeekFinal = SprSeekBonifie;
		prixInter = ( prixAct - prixAct2 ) / Xincr ;
		// if prixinter =  0 sortir erreur
		SprSeekBonifie = SprSeek - prixAct / prixInter ; 
		prixAct2 = prixAct ;		
	}
		
//**************** PRICE *******************
	VECTOR<double> vRes;

	vRes.push_back(SprSeekFinal);

	//**************** xe lissée************
	iteration =1;
	prix ;
//	double TauxFixe = 10.01;
	SprSeek = 10.01;
	SprSeekFinalLisse= 10.01;
	SprSeekBonifie = 10.01;  // PUREMENT ARBITRAIRE
	//fee = Ntl * fee / 100 ;
//******************** Calcul premiere PV ******************

	if ( SFDate != enddate )
	{
		long swaplegISFLisse;
		// A remplacer par le spread mis en argument
		if (ARMLOCAL_SWAPLEG(irindexIdSFleg,
							 SFDate,
							 enddate,
							 K_RCV,
//							 1,
//							 LocalGetNumObjectId(spdSF),
							 0,
							 SprSeekBonifie, 
							 false,
							 l_ccy,
							 daycountId,
							 resetgapSF,
							 "NULL",
							 "NULL",
							 0,
							 K_NX_NONE,
							 stubruleId,
							 refdate,
							 K_NO,
						  -1,
							 C_result) == ARM_OK)
		{
			swaplegISFLisse = C_result.getLong();
		}
		else
		{
			ERROR_MSG(C_result.getMsg(),pRet);
			return S_OK;
		}

		long cloneIdSFLisse;
		if (ARMLOCAL_ClonedAndSetNotional(swaplegISFLisse,
										  AmortRefvalue,
										  100.,
										  C_result) == ARM_OK)
		{
			cloneIdSFLisse = C_result.getLong();
		}
		else
		{
			ERROR_MSG(C_result.getMsg(),pRet);
			return S_OK;
		}

		//* Price Deuxième phase Phase structure TD*
		if (ARMLOCAL_ARM_Price(cloneIdSFLisse,
							   LocalGetNumObjectId(l_bsmodFund),
							   C_result) == ARM_OK)
		{
			pxSFLegLisse = C_result.getDouble();
		}
		else
		{
			ERROR_MSG(C_result.getMsg(),pRet);
			return S_OK;
		}

		if (ARMLOCAL_ARM_Price(cloneIdSFLisse,
							   LocalGetNumObjectId(l_bsmodDeltaFund),
							   C_result) == ARM_OK)
		{
			pxSFLegDeltaLisse = C_result.getDouble();
		}
		else
		{
			ERROR_MSG(C_result.getMsg(),pRet);
			return S_OK;
		}

		retCode = ARMLOCAL_FreeObject(swaplegISFLisse,C_result);
		retCode = ARMLOCAL_FreeObject(cloneIdSFLisse,C_result);		
	}
	else
	{
		pxSFLegLisse = 0.;
		pxSFLegDeltaLisse = 0.;

	}

	pvBiliborLeg = computeBilibor(irindexIdPF2leg, startdate, SFDate, SprSeekBonifie, l_ccy, daycountId, resetgap, stubruleId, 
						refdate, AmortRefvalue, l_bsmod, l_bsmodDeltaCcy1, l_bsmodDeltaCcy2, l_bsmodFxCorrel );

	if (pvBiliborLeg.size() == 0)
	{
		ERROR_MSG("Pb in computeBilibor",pRet);
		return S_OK;
	}

	pxPF1Leg = pvBiliborLeg[0];
	pxPF1LegDeltaCcy1 = pvBiliborLeg[1];
	pxPF1LegDeltaCcy2 = pvBiliborLeg[2];
	pxPF1LegFxCorrel = pvBiliborLeg[3];


	pxTotal = pxPF1Leg + pxPF2Leg + pxSFLegLisse - pxFundingLeg;
	pxDeltaCcy1 = (pxPF1LegDeltaCcy1 - pxPF1Leg) + (pxPF2LegDeltaCcy2 - pxPF2Leg) + (pxSFLegDeltaLisse-pxSFLegLisse) - (pxFundingLegDelta-pxFundingLeg); 
	if (pxDeltaCcy1 <0.) pxDeltaCcy1 = - pxDeltaCcy1;
	pxDeltaCcy2 = (pxPF1LegDeltaCcy2 - pxPF1Leg) + (pxPF2LegDeltaCcy2 - pxPF2Leg);
	if (pxDeltaCcy2 <0.) pxDeltaCcy2 = - pxDeltaCcy2;
	pxFxCorrel = (pxPF1LegFxCorrel - pxPF1Leg) + (pxPF2LegFxCorrel - pxPF2Leg);
	if (pxFxCorrel <0.) pxFxCorrel = - pxFxCorrel;

	SensiTotal = pxDeltaCcy1 + pxDeltaCcy2 + pxFxCorrel;
	prixAct = pxTotal - SensiTotal - fee; 


	prix = prixAct ;
	prixAct2 =1;
	
//********************                    ******************
//******************** BOUCLE SUR LE PRICE ******************
	while ((abs(prixAct2) > 0.0001) && (iteration < 10 ))
	{
		SprSeek = SprSeekBonifie;
		
		//*************************
		// **** Premier recalcul ****
		//****************************
		if ( SFDate != enddate )
		{
			long swaplegISFLisse;
			// A remplacer par le spread mis en argument
			if (ARMLOCAL_SWAPLEG(irindexIdSFleg,
								 SFDate,
								 enddate,
								 K_RCV,
	//							 1,
	//							 LocalGetNumObjectId(spdSF),
								 0,
								 SprSeekBonifie, 
								 false,
								 l_ccy,
								 daycountId,
								 resetgapSF,
								 "NULL",
								 "NULL",
								 0,
								 K_NX_NONE,
								 stubruleId,
								 refdate,
								 K_NO,
						  -1,
								 C_result) == ARM_OK)
			{
				swaplegISFLisse = C_result.getLong();
			}
			else
			{
				ERROR_MSG(C_result.getMsg(),pRet);
				return S_OK;
			}

			long cloneIdSFLisse;
			if (ARMLOCAL_ClonedAndSetNotional(swaplegISFLisse,
											  AmortRefvalue,
											  100.,
											  C_result) == ARM_OK)
			{
				cloneIdSFLisse = C_result.getLong();
			}
			else
			{
				ERROR_MSG(C_result.getMsg(),pRet);
				return S_OK;
			}

			//* Price Deuxième phase Phase structure TD*
			if (ARMLOCAL_ARM_Price(cloneIdSFLisse,
								   LocalGetNumObjectId(l_bsmodFund),
								   C_result) == ARM_OK)
			{
				pxSFLegLisse = C_result.getDouble();
			}
			else
			{
				ERROR_MSG(C_result.getMsg(),pRet);
				return S_OK;
			}

			if (ARMLOCAL_ARM_Price(cloneIdSFLisse,
								   LocalGetNumObjectId(l_bsmodDeltaFund),
								   C_result) == ARM_OK)
			{
				pxSFLegDeltaLisse = C_result.getDouble();
			}
			else
			{
				ERROR_MSG(C_result.getMsg(),pRet);
				return S_OK;
			}

			retCode = ARMLOCAL_FreeObject(swaplegISFLisse,C_result);
			retCode = ARMLOCAL_FreeObject(cloneIdSFLisse,C_result);				
		}
		else
		{
			pxSFLegLisse = 0.;
			pxSFLegDeltaLisse = 0.;
		}

		pvBiliborLeg = computeBilibor(irindexIdPF2leg, startdate, SFDate, SprSeekBonifie, l_ccy, daycountId, resetgap, stubruleId, 
						refdate, AmortRefvalue, l_bsmod, l_bsmodDeltaCcy1, l_bsmodDeltaCcy2, l_bsmodFxCorrel );

		if (pvBiliborLeg.size() == 0)
		{
			ERROR_MSG("Pb in computeBilibor",pRet);
			return S_OK;
		}

		pxPF1Leg = pvBiliborLeg[0];
		pxPF1LegDeltaCcy1 = pvBiliborLeg[1];
		pxPF1LegDeltaCcy2 = pvBiliborLeg[2];
		pxPF1LegFxCorrel = pvBiliborLeg[3];

		pxTotal = pxPF1Leg + pxPF2Leg + pxSFLegLisse - pxFundingLeg;
		pxDeltaCcy1 = (pxPF1LegDeltaCcy1 - pxPF1Leg) + (pxPF2LegDeltaCcy1 - pxPF2Leg) + (pxSFLegDeltaLisse-pxSFLegLisse) - (pxFundingLegDelta-pxFundingLeg); 
		if (pxDeltaCcy1 <0.) pxDeltaCcy1 = - pxDeltaCcy1;
		pxDeltaCcy2 = (pxPF1LegDeltaCcy2 - pxPF1Leg) + (pxPF2LegDeltaCcy2 - pxPF2Leg);
		if (pxDeltaCcy2 <0.) pxDeltaCcy2 = - pxDeltaCcy2;
		pxFxCorrel = (pxPF1LegFxCorrel - pxPF1Leg) + (pxPF2LegFxCorrel - pxPF2Leg);
		if (pxFxCorrel <0.) pxFxCorrel = - pxFxCorrel;

		SensiTotal = pxDeltaCcy1 + pxDeltaCcy2 + pxFxCorrel;
		prixAct = pxTotal - SensiTotal - fee; 
		// les deltas se compensent et on somme en valeur absolue les deltas et vegas

		//******************************
		// **** affectation valeurs ****
		//******************************		
		prixAct2 = prixAct ;
		Xincr = SprSeek / 1000000;
		SprSeekBonifie = SprSeek + Xincr;
		
		//****************************
		// **** Deuxieme recalcul ****
		//****************************

		if ( SFDate != enddate )
		{
			long swaplegISFLisse;
			
			// A remplacer par le spread mis en argument
			
			if (ARMLOCAL_SWAPLEG(irindexIdSFleg,
								 SFDate,
								 enddate,
								 K_RCV,
	//							 1,
	//							 LocalGetNumObjectId(spdSF),
								 0,
								 SprSeekBonifie, 
								 false,
								 l_ccy,
								 daycountId,
								 resetgapSF,
								 "NULL",
								 "NULL",
								 0,
								 K_NX_NONE,
								 stubruleId,
								 refdate,
								 K_NO,
						  -1,
								 C_result) == ARM_OK)
			{
				swaplegISFLisse = C_result.getLong();
			}
			else
			{
				ERROR_MSG(C_result.getMsg(),pRet);
				return S_OK;
			}

			long cloneIdSFLisse;
			if (ARMLOCAL_ClonedAndSetNotional(swaplegISFLisse,
											  AmortRefvalue,
											  100.,
											  C_result) == ARM_OK)
			{
				cloneIdSFLisse = C_result.getLong();
			}
			else
			{
				ERROR_MSG(C_result.getMsg(),pRet);
				return S_OK;
			}

			//* Price Deuxième phase Phase structure TD*
			if (ARMLOCAL_ARM_Price(cloneIdSFLisse,
								   LocalGetNumObjectId(l_bsmodFund),
								   C_result) == ARM_OK)
			{
				pxSFLegLisse = C_result.getDouble();
			}
			else
			{
				ERROR_MSG(C_result.getMsg(),pRet);
				return S_OK;
			}

			if (ARMLOCAL_ARM_Price(cloneIdSFLisse,
								   LocalGetNumObjectId(l_bsmodDeltaFund),
								   C_result) == ARM_OK)
			{
				pxSFLegDeltaLisse = C_result.getDouble();
			}
			else
			{
				ERROR_MSG(C_result.getMsg(),pRet);
				return S_OK;
			}

			retCode = ARMLOCAL_FreeObject(swaplegISFLisse,C_result);
			retCode = ARMLOCAL_FreeObject(cloneIdSFLisse,C_result);		
		
		}
		else
		{
			pxSFLegLisse = 0.;
			pxSFLegDeltaLisse = 0.;

		}

		pvBiliborLeg = computeBilibor(irindexIdPF2leg, startdate, SFDate, SprSeekBonifie, l_ccy, daycountId, resetgap, stubruleId, 
						refdate, AmortRefvalue, l_bsmod, l_bsmodDeltaCcy1, l_bsmodDeltaCcy2, l_bsmodFxCorrel );


		if (pvBiliborLeg.size() == 0)
		{
			ERROR_MSG("Pb in computeBilibor",pRet);
			return S_OK;
		}

		pxPF1Leg = pvBiliborLeg[0];
		pxPF1LegDeltaCcy1 = pvBiliborLeg[1];
		pxPF1LegDeltaCcy2 = pvBiliborLeg[2];
		pxPF1LegFxCorrel = pvBiliborLeg[3];

		pxTotal = pxPF1Leg + pxPF2Leg + pxSFLegLisse - pxFundingLeg;
		pxDeltaCcy1 = (pxPF1LegDeltaCcy1 - pxPF1Leg) + (pxPF2LegDeltaCcy1 - pxPF2Leg) + (pxSFLegDeltaLisse-pxSFLegLisse) - (pxFundingLegDelta-pxFundingLeg); 
		if (pxDeltaCcy1 <0.) pxDeltaCcy1 = - pxDeltaCcy1;
		pxDeltaCcy2 = (pxPF1LegDeltaCcy2 - pxPF1Leg) + (pxPF2LegDeltaCcy2 - pxPF2Leg);
		if (pxDeltaCcy2 <0.) pxDeltaCcy2 = - pxDeltaCcy2;
		pxFxCorrel = (pxPF1LegFxCorrel - pxPF1Leg) + (pxPF2LegFxCorrel - pxPF2Leg);
		if (pxFxCorrel <0.) pxFxCorrel = - pxFxCorrel;

		SensiTotal = pxDeltaCcy1 + pxDeltaCcy2 + pxFxCorrel;
		prixAct = pxTotal - SensiTotal - fee; 
	
		iteration++;
		
		SprSeekFinalLisse = SprSeekBonifie;
		prixInter = ( prixAct - prixAct2 ) / Xincr ;
		// if prixinter =  0 sortir erreur
		SprSeekBonifie = SprSeek - prixAct / prixInter ; 
		prixAct2 = prixAct ;
		
	}
	//***************fin Marge Lissée
	
	if ( SFDate != enddate )
	{		
		retCode = ARMLOCAL_FreeObject(irindexIdSFleg,C_result);
	}

	retCode = ARMLOCAL_FreeObject(irindexIdPF2leg,C_result);
	retCode = ARMLOCAL_FreeObject(AmortRefvalueToDestroy,C_result);
	vRes.push_back(SprSeekFinalLisse);
	Res = VECTORDOUBLE2VARIANT(vRes, pRet);

	_variant_t wrap_Res;
	wrap_Res.Attach(*pRet);
	*pRet=wrap_Res.Detach();

	return S_OK;
}
catch(...)
{
	return E_FAIL;
}
}



VECTOR<double> computeBilibor(long irindexIdPF2leg, double startdate, double SFDate, double SprSeekBonifie, CCString pccy, long daycountId, double resetgap, long stubruleId, 
						double prefdate, long pAmortRefvalue, CCString pbsmod, CCString pbsmodDeltaCcy1, CCString pbsmodDeltaCcy2, CCString pbsmodFxCorrel)
{
	VECTOR<double> resultat;
	resultat.clear();
	
try
{
	double pxPF1Leg;
	double pxPF1LegDeltaCcy1;
	double pxPF1LegDeltaCcy2; 
	double pxPF1LegFxCorrel;
	long swaplegIdPF1;

	ARM_result C_result;

	if (ARMLOCAL_SWAPLEG(irindexIdPF2leg,
						 startdate,
						 SFDate,
						 K_RCV, 
						 0,
						 SprSeekBonifie,
						 false,
						 pccy,
						 daycountId,
						 resetgap,
						 "NULL",
						 "NULL",
						 0,
						 K_NX_NONE,
						 stubruleId,
						 prefdate,
						 K_NO,
						  -1,
						 C_result) == ARM_OK)
	{
		swaplegIdPF1 = C_result.getLong();
	}
	else
	{
		return resultat;
	}


	long cloneIdPF1;
	if (ARMLOCAL_ClonedAndSetNotional(swaplegIdPF1,
									  pAmortRefvalue,
									  100.,
									  C_result) == ARM_OK)
	{
		cloneIdPF1 = C_result.getLong();
	}
	else
	{
		return resultat;
	}

	//*********   Price Première Phase structure 1  TD*************************


	if (ARMLOCAL_ARM_Price(cloneIdPF1,
						   LocalGetNumObjectId(pbsmod),
						   C_result) == ARM_OK)
	{
		pxPF1Leg = C_result.getDouble();
	}
	else
	{
		return resultat;
	}

	if (ARMLOCAL_ARM_Price(cloneIdPF1,
						   LocalGetNumObjectId(pbsmodDeltaCcy1),
						   C_result) == ARM_OK)
	{
		pxPF1LegDeltaCcy1 = C_result.getDouble();
	}
	else
	{
		return resultat;
	}

	if (ARMLOCAL_ARM_Price(cloneIdPF1,
						   LocalGetNumObjectId(pbsmodDeltaCcy2),
						   C_result) == ARM_OK)
	{
		pxPF1LegDeltaCcy2 = C_result.getDouble();
	}
	else
	{
		return resultat;
	}

	if (ARMLOCAL_ARM_Price(cloneIdPF1,
						   LocalGetNumObjectId(pbsmodFxCorrel),
						   C_result) == ARM_OK)
	{
		pxPF1LegFxCorrel = C_result.getDouble();
	}
	else
	{
		return resultat;
	}

	VECTOR<double> vRes;

	vRes.push_back(pxPF1Leg);
	vRes.push_back(pxPF1LegDeltaCcy1);
	vRes.push_back(pxPF1LegDeltaCcy2);
	vRes.push_back(pxPF1LegFxCorrel);

	long retCode = ARMLOCAL_FreeObject(swaplegIdPF1,C_result);
	retCode = ARMLOCAL_FreeObject(cloneIdPF1,C_result);

	return vRes;
}
catch(...)
{
	return resultat;
}
}

VECTOR<double> computeOptilix(long irindexIdPFleg, double startdate, double SFDate, double SprSeekBonifie, CCString pccy, long daycountId, double resetgap, long stubruleId, 
						double prefdate, long pAmortRefvalue, CCString pbsmod, CCString pbsmodDeltaCcy);


//Ajout TD le 16/08/2005
STDMETHODIMP ActiveXModule::ARMcomputeOptilix(double pAsOf,
								   double pStartDate,
								   double pDateSecondPhase,
								   double pEndDate,
								   double pNtl,
								   BSTR pIndex,
								   BSTR pIndexFund,
								   BSTR pIndexSF,
								   BSTR pBsmod,
								   BSTR pBsmodFund,
								   BSTR pBsmodDeltaCcy,
								   BSTR pBsmodDeltaFund,
								   BSTR pFreqP,
								   BSTR pFreqR,
								   BSTR pFreqPFund,
								   BSTR pFreqRFund, 
								   BSTR pFreqPSF,
								   BSTR pFreqRSF,
								   BSTR pResetTiming,
								   BSTR pResetTimingSF,
								   BSTR pCcy,
								   BSTR pDayCount,
								   BSTR pDayCountSF,
								   double pSpdSF,
								   VARIANT pSpdfund,
								   double pResetGap,
								   double pResetGapSF,
								   BSTR pAmort,
								   double pRefDate,
								   double pFee,
								   BSTR pFwdRule,
								   BSTR pIntRule,
								   BSTR pStubRule,
								   BSTR pFreqAmort,
								   double pTxAmort,
								   double pAmountAmort,
								   VARIANT *pRet)

{
	// TODO: Add your implementation code here
try
{
	ARM_result C_result;

	double AsOf = pAsOf;
	double startdate = pStartDate;
	double SFDate = pDateSecondPhase;
	double enddate = pEndDate;
	double Ntl = pNtl;
	double resetgap = pResetGap;
	double resetgapSF = pResetGapSF;
	double refdate = pRefDate;
	double fee = pFee;

	double pxFundingLeg;
	double pxFundingLegDelta;
	double pxSFLeg;
	double pxSFLegDelta;

	double prixInter;
	double Xincr;
	double prixAct;
	double prixAct2;

	double pxPF1Leg;
	double pxPF1LegDeltaCcy1;

	double pxTotal;
	double pxDeltaCcy1;

	double SensiTotal;

	double pxSFLegLisse;
	double pxSFLegDeltaLisse;

	VECTOR<double> pvOptilixLeg;
		
	double spdSF = pSpdSF;

	double spdFund;
	long   spdType;

	if (VARIANT2Double(pSpdfund,&spdFund) == S_FALSE)
	{
		_bstr_t bSpdFund(pSpdfund);
		CCString l_SpdFund = bSpdFund;
		spdFund = (double) LocalGetNumObjectId(l_SpdFund);
		spdType = 1;
	}
	else
	   spdType = 0;

	double SprSeekFinalLisse; 

	long retCode;
	double Res = 0.;

	_bstr_t bindex(pIndex);
	CCString l_index = bindex;

	_bstr_t bindexFund(pIndexFund);
	CCString l_indexFund = bindexFund;

	_bstr_t bindexSF(pIndexSF);
	CCString l_indexSF = bindexSF;

	_bstr_t bfreqp(pFreqP);
	CCString l_freqp = bfreqp;

	_bstr_t bfreqpFund(pFreqPFund);
	CCString l_freqpFund = bfreqpFund;

	_bstr_t bfreqrFund(pFreqRFund);
	CCString l_freqrFund = bfreqrFund;

	_bstr_t bfreqr(pFreqR);
	CCString l_freqr = bfreqr;
	
	_bstr_t bfreqpSF(pFreqPSF);
	CCString l_freqpSF = bfreqpSF;

	_bstr_t bfreqrSF(pFreqRSF);
	CCString l_freqrSF = bfreqrSF;

	_bstr_t bresettiming(pResetTiming);
	CCString l_resettiming = bresettiming;

	_bstr_t bresettimingSF(pResetTimingSF);
	CCString l_resettimingSF = bresettimingSF;

	_bstr_t bccy(pCcy);
	CCString l_ccy = bccy;

	_bstr_t bdaycount(pDayCount);
	CCString l_daycount = bdaycount;

	_bstr_t bdaycountSF(pDayCountSF);
	CCString l_daycountSF = bdaycountSF;

	_bstr_t bamort(pAmort);
	CCString l_amort = bamort;

	_bstr_t bbsmod(pBsmod);
	CCString l_bsmod = bbsmod;

	_bstr_t bbsmodFund(pBsmodFund);
	CCString l_bsmodFund = bbsmodFund;

	if ( l_bsmodFund == "DEFAULT")
	{
		l_bsmodFund = l_bsmod;
	}

	_bstr_t bbsmodDeltaCcy1(pBsmodDeltaCcy);
	CCString l_bsmodDeltaCcy1 = bbsmodDeltaCcy1;
	if ( l_bsmodDeltaCcy1 == "DEFAULT")
	{
		l_bsmodDeltaCcy1 = l_bsmod;
	}

	_bstr_t bbsmodDeltaFund(pBsmodDeltaFund);
	CCString l_bsmodDeltaFund = bbsmodDeltaFund;	
	if ( l_bsmodDeltaFund == "DEFAULT")
	{
		l_bsmodDeltaFund = l_bsmod;
	}

	_bstr_t bfwdrule(pFwdRule);
	CCString l_fwdrule = bfwdrule;

	_bstr_t bintrule(pIntRule);
	CCString l_intrule = bintrule;

	_bstr_t bstubrule(pStubRule);
	CCString l_stubrule = bstubrule;

	long freqPayId;
	if((freqPayId = ARM_ConvFrequency(l_freqp, C_result)) == ARM_DEFAULT_ERR)
	{
		ERROR_MSG("Invalid Frq",pRet,ARM_ERROR_FREQ);
		return S_OK;
	}

	long resetFreqId;
	if((resetFreqId = ARM_ConvFrequency(l_freqr, C_result)) == ARM_DEFAULT_ERR)
	{
		ERROR_MSG("Invalid Frq",pRet,ARM_ERROR_FREQ);
		return S_OK;
	}
	
	long freqPayIdSF;
	if((freqPayIdSF = ARM_ConvFrequency(l_freqpSF, C_result)) == ARM_DEFAULT_ERR)
	{
		ERROR_MSG("Invalid Frq",pRet,ARM_ERROR_FREQ);
		return S_OK;
	}

	long resetFreqIdSF;
	if((resetFreqIdSF = ARM_ConvFrequency(l_freqrSF, C_result)) == ARM_DEFAULT_ERR)
	{
		ERROR_MSG("Invalid Frq",pRet,ARM_ERROR_FREQ);
		return S_OK;
	}

	long freqPayIdFund;
	if((freqPayIdFund = ARM_ConvFrequency(l_freqpFund, C_result)) == ARM_DEFAULT_ERR)
	{
		ERROR_MSG("Invalid Frq",pRet,ARM_ERROR_FREQ);
		return S_OK;
	}

	long resetFreqIdFund;
	if((resetFreqIdFund = ARM_ConvFrequency(l_freqrFund, C_result)) == ARM_DEFAULT_ERR)
	{
		ERROR_MSG("Invalid Frq",pRet,ARM_ERROR_FREQ);
		return S_OK;
	}

	long indexId;
	indexId = ARM_ConvIrType(l_index);

	long indexIdFund;
	indexIdFund = ARM_ConvIrType(l_indexFund);

	long indexIdSF;
	indexIdSF = ARM_ConvIrType(l_indexSF);
	

	long fwdRuleId;
	if((fwdRuleId = ARM_ConvFwdRule(l_fwdrule, C_result)) == ARM_DEFAULT_ERR)
	{
		ERROR_MSG("Invalid Fwd Rule",pRet,ARM_ERROR_FREQ);
		return S_OK;
	}

	long daycountId= ARM_ConvDayCount(l_daycount);
	long daycountIdSF= ARM_ConvDayCount(l_daycountSF);
	long intruleId= ARM_ConvIntRule(l_intrule);
	long resettimingId = ARM_ConvPayResetRule (l_resettiming);
	long resettimingIdSF = ARM_ConvPayResetRule (l_resettimingSF);
	long stubruleId = ARM_ConvStubRule(l_stubrule);

	//*************************************************************
	//**********    Create Armort Refvalue ID   **********
	//*************************************************************

	long amortRefValue;
	long amortRefValueToDestroy = ARM_NULL_OBJECT;

	if (LocalGetNumObjectId(l_amort) == ARM_KO)
	{
		double tauxAmort = pTxAmort;
		double amountAmort = pAmountAmort;
		
		_bstr_t bFreqAmort(pFreqAmort);
		CCString l_FreqAmort = bFreqAmort;

		long amortFreqId;
		if((amortFreqId = ARM_ConvFrequency(l_FreqAmort, C_result)) == ARM_DEFAULT_ERR)
		{
			ERROR_MSG("Invalid Amort Frq",pRet,ARM_ERROR_FREQ);
			return S_OK;
		}

		// AmortMethodId
		long amortMethodId;

		if((amortMethodId = ARM_ConvAmortMethod (l_amort, C_result)) == ARM_DEFAULT_ERR)
		{
			ERROR_MSG("Invalid AmortMethod",pRet,ARM_ERROR_FREQ);
			return S_OK;
		}

		// la enddate est la SFDate pour un "une phase" 
		long FixedlegIdamort;
		if (ARMLOCAL_FIXEDLEG(startdate,
							  enddate,
							  K_RCV,
							  0,
							  1.0,
							  daycountId,
							  amortFreqId,
							  K_COMP_PROP,
							  K_ARREARS,
							  intruleId,
							  stubruleId,
							  false,
							  l_ccy,
							  "NULL",
							  K_NX_NONE,
							  -1.0,
							  K_NO,
						  -1,
							  C_result) == ARM_OK) 
							  
		{
			FixedlegIdamort = C_result.getLong();
		}
		else
		{
			ERROR_MSG(C_result.getMsg(),pRet);
			return S_OK;
		}

		// ********** Funding Leg ************
		if (ARMLOCAL_ARM_GenAmortization(FixedlegIdamort,
										 amortMethodId,
										 amortFreqId,
										 amountAmort,
										 -1,
										 Ntl,
										 tauxAmort,
										 0.,
										 ARM_NULL_OBJECT,
										 0.0,
										 C_result) == ARM_OK) 
		{
			amortRefValue = C_result.getLong();
			amortRefValueToDestroy = amortRefValue;
		}
		else
		{
			ERROR_MSG(C_result.getMsg(),pRet);
			return S_OK;
		}
	}
	else
	{
		amortRefValue = LocalGetNumObjectId(l_amort);

		// calcul du 1er nominal
		double dPaymentDate;

		if (ARMLOCAL_DisplayRefValue(amortRefValue,true,C_result) == ARM_OK)
		{
			dPaymentDate = C_result.getArray(0);
		}
		else
		{
			ERROR_MSG(C_result.getMsg(),pRet);
			return S_OK;
		}
		
		if (ARMLOCAL_CptRefValue(amortRefValue,dPaymentDate,C_result) == ARM_OK)
		{
			Ntl = C_result.getDouble();
		}
		else
		{
			ERROR_MSG(C_result.getMsg(),pRet);
			return S_OK;
		} 
	}
	//************************************************************
	//**********   END Create Amort table    ************
	//************************************************************



	//**********************************************************
	//**********     funding leg & Snd Phase fund    ***********
	//*********      funding en dur en ADV et gap -2  **********
	//**********************************************************
	long irindexIdfund;
	if (ARMLOCAL_IRINDEX(KACTUAL_360,
						 freqPayIdFund,
						 -1,
						 K_COMP_PROP,
						 K_MOD_FOLLOWING,
						 K_ADVANCE,
						 -2,
						 K_ARREARS,
						 10000.,
						 false,
						 l_ccy,
						 indexIdFund,
						 K_COMP_PROP,
						 intruleId,
						 resetFreqIdFund,
						 C_result) == ARM_OK) 
	{
		irindexIdfund = C_result.getLong();
	}
	else
	{
		ERROR_MSG(C_result.getMsg(),pRet);
		return S_OK;
	}
// la enddate est la SFDate pour un "une phase" 
	long swaplegId;
	if (ARMLOCAL_SWAPLEG(irindexIdfund,
						 startdate,
						 enddate,
						 K_RCV,
						 spdType, 
						 spdFund,
						 false,
						 l_ccy,
						 daycountId,
						 -2,
						 "NULL",
						 "NULL",
						 0,
						 K_NX_NONE,
						 stubruleId,
						 refdate,
						 K_NO,
						  -1,
						 C_result) == ARM_OK) 
	{
		swaplegId = C_result.getLong();
	}
	else
	{
		ERROR_MSG(C_result.getMsg(),pRet);
		return S_OK;
	}
	
	// ********** Funding Leg ************
	long cloneId;
	if (ARMLOCAL_ClonedAndSetNotional(swaplegId,
									  amortRefValue,
									  100.,
									  C_result) == ARM_OK)
	{
		cloneId = C_result.getLong();
	}
	else
	{
		ERROR_MSG(C_result.getMsg(),pRet);
		return S_OK;
	}

	//*********   Price funding Leg   *************************
	if (ARMLOCAL_ARM_Price(cloneId,
						   LocalGetNumObjectId(l_bsmodFund),
						   C_result) == ARM_OK)
	{
		pxFundingLeg = C_result.getDouble();
	}
	else
	{
		ERROR_MSG(C_result.getMsg(),pRet);
		return S_OK;
	}

	// shift seulement sur l'euro !?
	if (ARMLOCAL_ARM_Price(cloneId,
						   LocalGetNumObjectId(l_bsmodDeltaFund),
						   C_result) == ARM_OK)
	{
		pxFundingLegDelta = C_result.getDouble();
	}
	else
	{
		ERROR_MSG(C_result.getMsg(),pRet);
		return S_OK;
	}
	// Fin du calcul du prix de la funding Leg

//********************Premiere Phase Patte structure 2 TD****************
	long irindexIdPF2leg;	
	if (ARMLOCAL_IRINDEX(KACTUAL_360,
						 freqPayId,
						 -1,
						 K_COMP_PROP,
						 K_MOD_FOLLOWING,
						 resettimingId,
						 resetgap,
						 K_ARREARS,
						 10000.,
						 false,
						 l_ccy,
						 indexId,
						 K_COMP_PROP,
						 intruleId,
						 resetFreqId,
						 C_result) == ARM_OK) 
	{
		irindexIdPF2leg = C_result.getLong();
	}
	else
	{
		ERROR_MSG(C_result.getMsg(),pRet);
		return S_OK;
	}

	//********************Deuxième Phase Patte structure TD****************
	long irindexIdSFleg;

	if ( SFDate != enddate )
	{		
		if (ARMLOCAL_IRINDEX(KACTUAL_360,
							 freqPayIdSF,
							 -1,
							 K_COMP_PROP,
							 K_MOD_FOLLOWING,
							 resettimingIdSF,
							 resetgapSF,
							 K_ARREARS,
							 10000.,
							 false,
							 l_ccy,
							 indexIdSF,
							 K_COMP_PROP,
							 intruleId,
							 resetFreqIdSF,
							 C_result) == ARM_OK) 
		{
			irindexIdSFleg = C_result.getLong();
		}
		else
		{
			ERROR_MSG(C_result.getMsg(),pRet);
			return S_OK;
		}

		long swaplegISF;

		// A remplacer par le spread mis en argument		
		if (ARMLOCAL_SWAPLEG(irindexIdSFleg,
							 SFDate,
							 enddate,
							 K_RCV,
							 0,
							 spdSF, 
							 false,
							 l_ccy,
							 daycountId,
							 resetgapSF,
							 "NULL",
							 "NULL",
							 0,
							 K_NX_NONE,
							 stubruleId,
							 refdate,
							 K_NO,
						  -1,
							 C_result) == ARM_OK)
		{
			swaplegISF = C_result.getLong();
		}
		else
		{
			ERROR_MSG(C_result.getMsg(),pRet);
			return S_OK;
		}

		long cloneIdSF;
		if (ARMLOCAL_ClonedAndSetNotional(swaplegISF,
										  amortRefValue,
										  100.,
										  C_result) == ARM_OK)
		{
			cloneIdSF = C_result.getLong();
		}
		else
		{
			ERROR_MSG(C_result.getMsg(),pRet);
			return S_OK;
		}

		//* Price Deuxième phase Phase structure TD*
		if (ARMLOCAL_ARM_Price(cloneIdSF,
							   LocalGetNumObjectId(l_bsmodFund),
							   C_result) == ARM_OK)
		{
			pxSFLeg = C_result.getDouble();
		}
		else
		{
			ERROR_MSG(C_result.getMsg(),pRet);
			return S_OK;
		}

		if (ARMLOCAL_ARM_Price(cloneIdSF,
							   LocalGetNumObjectId(l_bsmodDeltaFund),
							   C_result) == ARM_OK)
		{
			pxSFLegDelta = C_result.getDouble();
		}
		else
		{
			ERROR_MSG(C_result.getMsg(),pRet);
			return S_OK;
		}

		retCode = ARMLOCAL_FreeObject(swaplegISF,C_result);
		retCode = ARMLOCAL_FreeObject(cloneIdSF,C_result);	
	}
	else
	{
		pxSFLeg = 0.;
		pxSFLegDelta = 0.;
	}

//******************** Calcul premiere PV ******************
	int iteration =1;
	double prix ;
	double SprSeek = 10.01;
	double SprSeekFinal= 10.01;
	double SprSeekBonifie = 10.01;  // PUREMENT ARBITRAIRE
	fee = Ntl * fee / 100 ;
//******************** Calcul premiere PV ******************
	
	pvOptilixLeg = computeOptilix(irindexIdPF2leg, startdate, SFDate, SprSeekBonifie, l_ccy, daycountId, resetgap, stubruleId, 
						refdate, amortRefValue, l_bsmod, l_bsmodDeltaCcy1);

	if (pvOptilixLeg.size() == 0)
	{
		ERROR_MSG("Pb in computeOptilix",pRet);
		return S_OK;
	}

	pxPF1Leg = pvOptilixLeg[0];
	pxPF1LegDeltaCcy1 = pvOptilixLeg[1];

	pxTotal = pxPF1Leg + pxSFLeg - pxFundingLeg;
	pxDeltaCcy1 = (pxPF1LegDeltaCcy1 - pxPF1Leg)  + (pxSFLegDelta-pxSFLeg) - (pxFundingLegDelta-pxFundingLeg); 
	if (pxDeltaCcy1 <0.) pxDeltaCcy1 = - pxDeltaCcy1;

	SensiTotal = pxDeltaCcy1;
	prixAct = pxTotal - SensiTotal - fee; 

	prix = prixAct ;
	prixAct2 =1;
	
//********************                    ******************
//******************** BOUCLE SUR LE PRICE ******************
	while ((abs(prixAct2) > 0.0001) && (iteration < 10 ))
	{
		SprSeek = SprSeekBonifie;
		
		//*************************
		// **** Premier recalcul ****
		//****************************

		pvOptilixLeg = computeOptilix(irindexIdPF2leg, startdate, SFDate, SprSeekBonifie, l_ccy, daycountId, resetgap, stubruleId, 
						refdate, amortRefValue, l_bsmod, l_bsmodDeltaCcy1 );

		if (pvOptilixLeg.size() == 0)
		{
			ERROR_MSG("Pb in computeOptilix",pRet);
			return S_OK;
		}

		pxPF1Leg = pvOptilixLeg[0];
		pxPF1LegDeltaCcy1 = pvOptilixLeg[1];

		pxTotal = pxPF1Leg + pxSFLeg - pxFundingLeg;
		pxDeltaCcy1 = (pxPF1LegDeltaCcy1 - pxPF1Leg)  + (pxSFLegDelta-pxSFLeg) - (pxFundingLegDelta-pxFundingLeg); 
		if (pxDeltaCcy1 <0.) pxDeltaCcy1 = - pxDeltaCcy1;

		SensiTotal = pxDeltaCcy1;
		prixAct = pxTotal - SensiTotal - fee; 
		// les deltas se compensent et on somme en valeur absolue les deltas et vegas

		//******************************
		// **** affectation valeurs ****
		//******************************		
		prixAct2 = prixAct ;
		Xincr = SprSeek / 1000000;
		SprSeekBonifie = SprSeek + Xincr;		

		//****************************
		// **** Deuxieme recalcul ****
		//****************************
		pvOptilixLeg = computeOptilix(irindexIdPF2leg, startdate, SFDate, SprSeekBonifie, l_ccy, daycountId, resetgap, stubruleId, 
						refdate, amortRefValue, l_bsmod, l_bsmodDeltaCcy1);

		if (pvOptilixLeg.size() == 0)
		{
			ERROR_MSG("Pb in computeOptilix",pRet);
			return S_OK;
		}

		pxPF1Leg = pvOptilixLeg[0];
		pxPF1LegDeltaCcy1 = pvOptilixLeg[1];

		pxTotal = pxPF1Leg + pxSFLeg - pxFundingLeg;
		pxDeltaCcy1 = (pxPF1LegDeltaCcy1 - pxPF1Leg) + (pxSFLegDelta-pxSFLeg) - (pxFundingLegDelta-pxFundingLeg); 
		if (pxDeltaCcy1 <0.) pxDeltaCcy1 = - pxDeltaCcy1;

		SensiTotal = pxDeltaCcy1;
		prixAct = pxTotal - SensiTotal - fee; 	

		iteration++;
		
		SprSeekFinal = SprSeekBonifie;
		prixInter = ( prixAct - prixAct2 ) / Xincr ;
		// if prixinter =  0 sortir erreur
		SprSeekBonifie = SprSeek - prixAct / prixInter ; 
		prixAct2 = prixAct ;		
	}
		
//**************** PRICE *******************

	VECTOR<double> vRes;
	vRes.push_back(SprSeekFinal);

	//**************** xe lissée************
	iteration =1;
	prix ;
	SprSeek = 10.01;
	SprSeekFinalLisse= 10.01;
	SprSeekBonifie = 10.01;  // PUREMENT ARBITRAIRE
//******************** Calcul premiere PV ******************

	if ( SFDate != enddate )
	{
		long swaplegISFLisse;
		
		// A remplacer par le spread mis en argument	
		if (ARMLOCAL_SWAPLEG(irindexIdSFleg,
							 SFDate,
							 enddate,
							 K_RCV,
							 0,
							 SprSeekBonifie, 
							 false,
							 l_ccy,
							 daycountId,
							 resetgapSF,
							 "NULL",
							 "NULL",
							 0,
							 K_NX_NONE,
							 stubruleId,
							 refdate,
							 K_NO,
						  -1,
							 C_result) == ARM_OK)
		{
			swaplegISFLisse = C_result.getLong();
		}
		else
		{
			ERROR_MSG(C_result.getMsg(),pRet);
			return S_OK;
		}

		long cloneIdSFLisse;
		if (ARMLOCAL_ClonedAndSetNotional(swaplegISFLisse,
										  amortRefValue,
										  100.,
										  C_result) == ARM_OK)
		{
			cloneIdSFLisse = C_result.getLong();
		}
		else
		{
			ERROR_MSG(C_result.getMsg(),pRet);
			return S_OK;
		}

	//* Price Deuxième phase Phase structure TD*
		if (ARMLOCAL_ARM_Price(cloneIdSFLisse,
							   LocalGetNumObjectId(l_bsmodFund),
							   C_result) == ARM_OK)
		{
			pxSFLegLisse = C_result.getDouble();
		}
		else
		{
			ERROR_MSG(C_result.getMsg(),pRet);
			return S_OK;
		}

		if (ARMLOCAL_ARM_Price(cloneIdSFLisse,
							   LocalGetNumObjectId(l_bsmodDeltaFund),
							   C_result) == ARM_OK)
		{
			pxSFLegDeltaLisse = C_result.getDouble();
		}
		else
		{
			ERROR_MSG(C_result.getMsg(),pRet);
			return S_OK;
		}

		retCode = ARMLOCAL_FreeObject(swaplegISFLisse,C_result);
		retCode = ARMLOCAL_FreeObject(cloneIdSFLisse,C_result);		
	}
	else
	{
		pxSFLegLisse = 0.;
		pxSFLegDeltaLisse = 0.;
	}

	pvOptilixLeg = computeOptilix(irindexIdPF2leg, startdate, SFDate, SprSeekBonifie, l_ccy, daycountId, resetgap, stubruleId, 
						refdate, amortRefValue, l_bsmod, l_bsmodDeltaCcy1);


	if (pvOptilixLeg.size() == 0)
	{
		ERROR_MSG("Pb in computeOptilix",pRet);
		return S_OK;
	}

	pxPF1Leg = pvOptilixLeg[0];
	pxPF1LegDeltaCcy1 = pvOptilixLeg[1];

	pxTotal = pxPF1Leg + pxSFLegLisse - pxFundingLeg;
	pxDeltaCcy1 = (pxPF1LegDeltaCcy1 - pxPF1Leg) + (pxSFLegDeltaLisse-pxSFLegLisse) - (pxFundingLegDelta-pxFundingLeg); 
	if (pxDeltaCcy1 <0.) pxDeltaCcy1 = - pxDeltaCcy1;

	SensiTotal = pxDeltaCcy1;
	prixAct = pxTotal - SensiTotal - fee; 

	prix = prixAct ;
	prixAct2 =1;

//********************                    ******************
//******************** BOUCLE SUR LE PRICE ******************
	while ((abs(prixAct2) > 0.0001) && (iteration < 10 ))
	{
		SprSeek = SprSeekBonifie;
		
		//*************************
		// **** Premier recalcul ****
		//****************************
		if ( SFDate != enddate)
		{
			long swaplegISFLisse;
			
			// A remplacer par le spread mis en argument
			if (ARMLOCAL_SWAPLEG(irindexIdSFleg,
								 SFDate,
								 enddate,
								 K_RCV,
								 0,
								 SprSeekBonifie, 
								 false,
								 l_ccy,
								 daycountId,
								 resetgapSF,
								 "NULL",
								 "NULL",
								 0,
								 K_NX_NONE,
								 stubruleId,
								 refdate,
								 K_NO,
						  -1,
								 C_result) == ARM_OK)
			{
				swaplegISFLisse = C_result.getLong();
			}
			else
			{
				ERROR_MSG(C_result.getMsg(),pRet);
				return S_OK;
			}

			long cloneIdSFLisse;
			if (ARMLOCAL_ClonedAndSetNotional(swaplegISFLisse,
											  amortRefValue,
											  100.,
											  C_result) == ARM_OK)
			{
				cloneIdSFLisse = C_result.getLong();
			}
			else
			{
				ERROR_MSG(C_result.getMsg(),pRet);
				return S_OK;
			}

		//* Price Deuxième phase Phase structure TD*
			if (ARMLOCAL_ARM_Price(cloneIdSFLisse,
								   LocalGetNumObjectId(l_bsmodFund),
								   C_result) == ARM_OK)
			{
				pxSFLegLisse = C_result.getDouble();
			}
			else
			{
				ERROR_MSG(C_result.getMsg(),pRet);
				return S_OK;
			}

			if (ARMLOCAL_ARM_Price(cloneIdSFLisse,
								   LocalGetNumObjectId(l_bsmodDeltaFund),
								   C_result) == ARM_OK)
			{
				pxSFLegDeltaLisse = C_result.getDouble();
			}
			else
			{
				ERROR_MSG(C_result.getMsg(),pRet);
				return S_OK;
			}

			retCode = ARMLOCAL_FreeObject(swaplegISFLisse,C_result);
			retCode = ARMLOCAL_FreeObject(cloneIdSFLisse,C_result);		
		}
		else
		{
			pxSFLegLisse = 0.;
			pxSFLegDeltaLisse = 0.;

		}

		pvOptilixLeg = computeOptilix(irindexIdPF2leg, startdate, SFDate, SprSeekBonifie, l_ccy, daycountId, resetgap, stubruleId, 
						refdate, amortRefValue, l_bsmod, l_bsmodDeltaCcy1);

		if (pvOptilixLeg.size() == 0)
		{
			ERROR_MSG("Pb in computeOptilix",pRet);
			return S_OK;
		}

		pxPF1Leg = pvOptilixLeg[0];
		pxPF1LegDeltaCcy1 = pvOptilixLeg[1];

		pxTotal = pxPF1Leg +pxSFLegLisse - pxFundingLeg;
		pxDeltaCcy1 = (pxPF1LegDeltaCcy1 - pxPF1Leg) + (pxSFLegDeltaLisse-pxSFLegLisse) - (pxFundingLegDelta-pxFundingLeg); 
		if (pxDeltaCcy1 <0.) pxDeltaCcy1 = - pxDeltaCcy1;

		SensiTotal = pxDeltaCcy1;
		prixAct = pxTotal - SensiTotal - fee; 
		// les deltas se compensent et on somme en valeur absolue les deltas et vegas
		
		//******************************
		// **** affectation valeurs ****
		//******************************
		prixAct2 = prixAct ;
		Xincr = SprSeek / 1000000;
		SprSeekBonifie = SprSeek + Xincr;		

		//****************************
		// **** Deuxieme recalcul ****
		//****************************

		if ( SFDate != enddate )
		{
			long swaplegISFLisse;
			
			// A remplacer par le spread mis en argument
			
			if (ARMLOCAL_SWAPLEG(irindexIdSFleg,
								 SFDate,
								 enddate,
								 K_RCV,
								 0,
								 SprSeekBonifie, 
								 false,
								 l_ccy,
								 daycountId,
								 resetgapSF,
								 "NULL",
								 "NULL",
								 0,
								 K_NX_NONE,
								 stubruleId,
								 refdate,
								 K_NO,
								 -1,
								 C_result) == ARM_OK)
			{
				swaplegISFLisse = C_result.getLong();
			}
			else
			{
				ERROR_MSG(C_result.getMsg(),pRet);
				return S_OK;
			}

			long cloneIdSFLisse;
			if (ARMLOCAL_ClonedAndSetNotional(swaplegISFLisse,
											  amortRefValue,
											  100.,
											  C_result) == ARM_OK)
			{
				cloneIdSFLisse = C_result.getLong();
			}
			else
			{
				ERROR_MSG(C_result.getMsg(),pRet);
				return S_OK;
			}

		//* Price Deuxième phase Phase structure TD*	
			if (ARMLOCAL_ARM_Price(cloneIdSFLisse,
								   LocalGetNumObjectId(l_bsmodFund),
								   C_result) == ARM_OK)
			{
				pxSFLegLisse = C_result.getDouble();
			}
			else
			{
				ERROR_MSG(C_result.getMsg(),pRet);
				return S_OK;
			}

			if (ARMLOCAL_ARM_Price(cloneIdSFLisse,
								   LocalGetNumObjectId(l_bsmodDeltaFund),
								   C_result) == ARM_OK)
			{
				pxSFLegDeltaLisse = C_result.getDouble();
			}
			else
			{
				ERROR_MSG(C_result.getMsg(),pRet);
				return S_OK;
			}

			retCode = ARMLOCAL_FreeObject(swaplegISFLisse,C_result);
			retCode = ARMLOCAL_FreeObject(cloneIdSFLisse,C_result);	
		}
		else
		{
			pxSFLegLisse = 0.;
			pxSFLegDeltaLisse = 0.;
		}

		pvOptilixLeg = computeOptilix(irindexIdPF2leg, startdate, SFDate, SprSeekBonifie, l_ccy, daycountId, resetgap, stubruleId, 
						refdate, amortRefValue, l_bsmod, l_bsmodDeltaCcy1 );


		if (pvOptilixLeg.size() == 0)
		{
			ERROR_MSG("Pb in computeOptilix",pRet);
			return S_OK;
		}

		pxPF1Leg = pvOptilixLeg[0];
		pxPF1LegDeltaCcy1 = pvOptilixLeg[1];

		pxTotal = pxPF1Leg + pxSFLegLisse - pxFundingLeg;
		pxDeltaCcy1 = (pxPF1LegDeltaCcy1 - pxPF1Leg)  + (pxSFLegDeltaLisse-pxSFLegLisse) - (pxFundingLegDelta-pxFundingLeg); 
		if (pxDeltaCcy1 <0.) pxDeltaCcy1 = - pxDeltaCcy1;

		SensiTotal = pxDeltaCcy1;
		prixAct = pxTotal - SensiTotal - fee; 

		iteration++;

		SprSeekFinalLisse = SprSeekBonifie;
		prixInter = ( prixAct - prixAct2 ) / Xincr ;
		// if prixinter =  0 sortir erreur
		SprSeekBonifie = SprSeek - prixAct / prixInter ; 
		prixAct2 = prixAct ;	
	}

	//***************fin Marge Lissée
	if ( SFDate != enddate )
	{		
		retCode = ARMLOCAL_FreeObject(irindexIdSFleg,C_result);
	}

	retCode = ARMLOCAL_FreeObject(amortRefValueToDestroy,C_result);
	retCode = ARMLOCAL_FreeObject(irindexIdPF2leg,C_result);
	vRes.push_back(SprSeekFinalLisse);
	Res = VECTORDOUBLE2VARIANT(vRes, pRet);

	_variant_t wrap_Res;
	wrap_Res.Attach(*pRet);
	*pRet=wrap_Res.Detach();

	return S_OK;
}
catch(...)
{
	return E_FAIL;
}
}


VECTOR<double> computeOptilix(long irindexIdPFleg, double startdate, double SFDate, double SprSeekBonifie, CCString pccy, long daycountId, double resetgap, long stubruleId, 
						double prefdate, long pAmortRefvalue, CCString pbsmod, CCString pbsmodDeltaCcy)
{
	VECTOR<double> resultat;
	resultat.clear();

try
{
	double pxPFLeg;
	double pxPFLegDeltaCcy;
	long swaplegIdPF;

	ARM_result C_result;

	if (ARMLOCAL_SWAPLEG(irindexIdPFleg,
						 startdate,
						 SFDate,
						 K_RCV, 
						 0,
						 SprSeekBonifie,
						 false,
						 pccy,
						 daycountId,
						 resetgap,
						 "NULL",
						 "NULL",
						 0,
						 K_NX_NONE,
						 stubruleId,
						 prefdate,
						 K_NO,
								 -1,
						 C_result) == ARM_OK)
	{
		swaplegIdPF = C_result.getLong();
	}
	else
	{
		return resultat;
	}

	long cloneIdPF;
	if (ARMLOCAL_ClonedAndSetNotional(swaplegIdPF,
									  pAmortRefvalue,
									  100.,
									  C_result) == ARM_OK)
	{
		cloneIdPF = C_result.getLong();
	}
	else
	{
		return resultat;
	}

	//*********   Price Première Phase structure 1  TD*************************
	if (ARMLOCAL_ARM_Price(cloneIdPF,
						   LocalGetNumObjectId(pbsmod),
						   C_result) == ARM_OK)
	{
		pxPFLeg = C_result.getDouble();
	}
	else
	{
		return resultat;
	}

	if (ARMLOCAL_ARM_Price(cloneIdPF,
						   LocalGetNumObjectId(pbsmodDeltaCcy),
						   C_result) == ARM_OK)
	{
		pxPFLegDeltaCcy = C_result.getDouble();
	}
	else
	{
		return resultat;
	}


	VECTOR<double> vRes;

	vRes.push_back(pxPFLeg);
	vRes.push_back(pxPFLegDeltaCcy);

	long retCode = ARMLOCAL_FreeObject(swaplegIdPF,C_result);
	retCode = ARMLOCAL_FreeObject(cloneIdPF,C_result);

	return vRes;
}
catch(...)
{
	return resultat;
}
}

VECTOR<double> computeDeltaByPlot(double pxProduct, long cloneProductId, const VECTOR<CCString>& p_vBumpBsGenMod, bool isCum = true)
{

VECTOR<double> resultat;
resultat.clear();
try 
{
	VECTOR<double> vRes;
	double pxBumpProduct;
	
	ARM_result C_result;

	VECTOR<double> resultat;
	resultat.clear();

	double previouspxProduct = pxProduct;

	long sizeMatu = p_vBumpBsGenMod.size(); 
	for (int i = 0; i < sizeMatu; i++)	
	{
		if (ARMLOCAL_ARM_Price(cloneProductId,
							   LocalGetNumObjectId(p_vBumpBsGenMod[i]),
							   C_result) == ARM_OK)
		{
			pxBumpProduct = C_result.getDouble();
			
			if (isCum)
			{
				vRes.push_back(pxBumpProduct - previouspxProduct);
				previouspxProduct = pxBumpProduct;
			}
			else
			{
				vRes.push_back(pxBumpProduct - pxProduct);
			}
		}
		else
		{
			return resultat;
		}
	}	
	return vRes;
}
catch(...)
{
	return resultat;
}

}



VECTOR<double> computePentifix(const CCString& index, const CCString& indexLong, double startDate, double endDate,
							   const CCString& floorOrCap, const CCString& strike, double payOff, long ccyId, long dayCountId,
							   long resetTimingId, long payTimingId, double resetGap, long resetFreqId, long payFreqId,
							   long intRuleId, long stubRuleId, long amortRefvalue,
								VECTOR<double> vDeltaOther, const VECTOR<CCString>& vBumpBsGenMod, const VECTOR<CCString>& vBumpVolBsGenMod, 
							   const CCString& hyperCubeCorrel, const CCString& VolCurvFromMatriceShift,
							   const CCString& vol, const CCString& volCUB, const CCString& CorrManager,
							   const CCString& ConvexityManager, const CCString& zc, const CCString& smiledMod);


STDMETHODIMP ActiveXModule::ARMcomputePentifix(double pNtl, 
										   double pStartDatePhase1,  
										   BSTR pCcy,
										   BSTR pIndexPhase1, 										   
										   double pSpreadPhase1,
										   BSTR pDayCountPhase1,
										   BSTR pPayFreqPhase1,
										   BSTR pResetFreqPhase1,
										   BSTR pResetTimingPhase1,
										   BSTR pRoll, 
										   BSTR pAdjPhase1, 
										   BSTR pStub,
										   BSTR pIndexPhase2DIG,
										   BSTR pIndexLongPhase2DIG,
										   BSTR pStrikePhase2DIG,
										   BSTR pResetTimingPhase2DIG,
										   BSTR pAdjPhase2DIG, 
										   double pStartDatePhase2,
										   double pSpreadPhase2,										   
										   BSTR pDayCountPhase2, 
										   BSTR pPayFreqPhase2, 
										   BSTR pResetFreqPhase2,
										   BSTR pAdjPhase2,
										   double pStartDatePhase3,
										   double pEndDatePhase3,
										   BSTR pIndexPhase3,
										   double pSpreadPhase3, 
										   BSTR pDayCountPhase3, 
										   BSTR pPayFreqPhase3, 
										   BSTR pResetFreqPhase3,
										   BSTR pResetTimingPhase3,
										   BSTR pAdjPhase3,
										   BSTR pIndexFund,
										   VARIANT pSpreadFund,
										   BSTR pDayCountFund, 
										   BSTR pPayFreqFund, 
										   BSTR pResetFreqFund,
										   BSTR pResetTimingFund,
										   BSTR pAdjFund,  
										   double pEndDateAmort,
										   BSTR pDayCountAmort,
										   BSTR pIntRuleAmort, 
										   double pTxAmort,
										   BSTR pFreqAmort,
										   double pAmountAmort, 
										   BSTR pTypeAmort, 
										   BSTR pFloorOrCap, 
										   double pFee, 
										   BSTR pVolCurvFromMatriceShift,
										   BSTR pVol, 
										   BSTR pVolCub, 
										   BSTR pCorrManager, 
										   BSTR pConvexityManager, 
										   BSTR pZc, 
										   BSTR pSmiledMod, 
										   BSTR pSmiledModBump, 					
										   BSTR pHyperCubeCorrel,
										   VARIANT *pBumpBsGenMod,
										   VARIANT *pBumpVolBsGenMod, 
										   VARIANT *pRet
										   )

{	
	// TODO: Add your implementation code here
try
{									   
	ARM_result C_result;

	VECTOR<double> vDate; 
	VECTOR<double> vSpread; 

	VECTOR<double> deltaPhase1; 
	VECTOR<double> deltaPhase2; 
	VECTOR<double> deltaPhase3; 
	VECTOR<double> deltaPhase2NOTDIG;
	VECTOR<double> deltaOther; 
	VECTOR<double> deltaFund; 
	double res = 0.;
	

	VECTOR<CCString> vBumpBsGenMod; 
	VECTOR<CCString>  vBumpVolBsGenMod; 

	double startDatePhase1 = pStartDatePhase1;  

	double spreadPhase1 = pSpreadPhase1;
	
	double startDatePhase2 = pStartDatePhase2;
	double spreadPhase2 = pSpreadPhase2;

	double startDatePhase3 = pStartDatePhase3;
	double endDatePhase3 = pEndDatePhase3;
	double spreadPhase3 = pSpreadPhase3;

	double spreadFund;
	long   spreadType;
	long sizeMatu;

	if (VARIANT2Double(pSpreadFund,&spreadFund) == S_FALSE)
	{
		_bstr_t bSpreadFund(pSpreadFund);
		CCString l_SpreadFund= bSpreadFund;
		spreadFund = (double) LocalGetNumObjectId(l_SpreadFund);
		spreadType = 1;
	}
	else
	   spreadType = 0;

	if(VARIANT2VECTORCCSTRING (*pBumpBsGenMod,vBumpBsGenMod,sizeMatu) != S_OK)
	{
		ERROR_MSG("Error conversion Model Delta",pRet,ARM_ERROR_FREQ);
		return S_OK;
	}

	if(VARIANT2VECTORCCSTRING (*pBumpVolBsGenMod,vBumpVolBsGenMod,sizeMatu) != S_OK)
	{
		ERROR_MSG("Error conversion Model Vega",pRet,ARM_ERROR_FREQ);
		return S_OK;
	}

	double Ntl = pNtl;
	double fee = pFee; 

	long retCode;
	double Res = 0.;

	_bstr_t ccy(pCcy);
	CCString l_ccy = ccy;

	_bstr_t indexPhase1(pIndexPhase1);
	CCString l_indexPhase1 = indexPhase1;

	_bstr_t dayCountPhase1(pDayCountPhase1);
	CCString l_dayCountPhase1 = dayCountPhase1;

	_bstr_t resetFreqPhase1(pResetFreqPhase1);
	CCString l_resetFreqPhase1 = resetFreqPhase1;

	_bstr_t payFreqPhase1(pPayFreqPhase1);
	CCString l_payFreqPhase1 = payFreqPhase1;

	_bstr_t resetTimingPhase1(pResetTimingPhase1);
	CCString l_resetTimingPhase1 = resetTimingPhase1;

	_bstr_t roll(pRoll);
	CCString l_roll = roll;

	_bstr_t adjPhase1(pAdjPhase1);
	CCString l_adjPhase1 = adjPhase1;
	
	_bstr_t stub(pStub);
	CCString l_stub= stub;

	_bstr_t indexPhase2DIG(pIndexPhase2DIG);
	CCString l_indexPhase2DIG = indexPhase2DIG;

	_bstr_t indexLongPhase2DIG(pIndexLongPhase2DIG);
	CCString l_indexLongPhase2DIG = indexLongPhase2DIG;

	//Strike is double??
	_bstr_t strikePhase2DIG(pStrikePhase2DIG);
	CCString l_strikePhase2DIG= strikePhase2DIG;

	_bstr_t resetTimingPhase2DIG(pResetTimingPhase2DIG);
	CCString l_resetTimingPhase2DIG = resetTimingPhase2DIG;

	_bstr_t adjPhase2DIG(pAdjPhase2DIG);
	CCString l_adjPhase2DIG = adjPhase2DIG;

	_bstr_t dayCountPhase2(pDayCountPhase2);
	CCString l_dayCountPhase2= dayCountPhase2;

	_bstr_t payFreqPhase2(pPayFreqPhase2);
	CCString l_payFreqPhase2= payFreqPhase2;
	
	_bstr_t resetFreqPhase2(pResetFreqPhase2);
	CCString l_resetFreqPhase2 = resetFreqPhase2;

	_bstr_t adjPhase2(pAdjPhase2);
	CCString l_adjPhase2 = adjPhase2;

	_bstr_t indexPhase3(pIndexPhase3);
	CCString l_indexPhase3 = indexPhase3;

	_bstr_t dayCountPhase3(pDayCountPhase3);
	CCString l_dayCountPhase3 = dayCountPhase3;

	_bstr_t payFreqPhase3(pPayFreqPhase3);
	CCString l_payFreqPhase3 = payFreqPhase3;

	_bstr_t resetFreqPhase3(pResetFreqPhase3);
	CCString l_resetFreqPhase3 = resetFreqPhase3;

	_bstr_t resetTimingPhase3(pResetTimingPhase3);
	CCString l_resetTimingPhase3 = resetTimingPhase3;

	_bstr_t adjPhase3(pAdjPhase3);
	CCString l_adjPhase3 = adjPhase3;

	_bstr_t indexFund(pIndexFund);
	CCString l_indexFund = indexFund;

	_bstr_t dayCountFund(pDayCountFund);
	CCString l_dayCountFund= dayCountFund;

	_bstr_t payFreqFund(pPayFreqFund);
	CCString l_payFreqFund= payFreqFund;
	
	_bstr_t resetFreqFund(pResetFreqFund);
	CCString l_resetFreqFund = resetFreqFund;

	_bstr_t resetTimingFund(pResetTimingFund);
	CCString l_resetTimingFund = resetTimingFund;

	_bstr_t adjFund(pAdjFund);
	CCString l_adjFund= adjFund;
	
	_bstr_t typeAmort(pTypeAmort);
	CCString l_typeAmort = typeAmort;

	_bstr_t floorOrCap (pFloorOrCap);
	CCString l_floorOrCap = floorOrCap; 

	_bstr_t volCurvFromMatriceShift(pVolCurvFromMatriceShift);
	CCString l_volCurvFromMatriceShift= volCurvFromMatriceShift;
	
	_bstr_t vol(pVol);
	CCString l_vol = vol;

	_bstr_t volCub(pVolCub);
	CCString l_volCub = volCub;

	_bstr_t corrManager(pCorrManager);
	CCString l_corrManager = corrManager;

	_bstr_t convexityManager(pConvexityManager);
	CCString l_convexityManager = convexityManager;

	_bstr_t zc(pZc);
	CCString l_zc = zc;

	_bstr_t smiledMod (pSmiledMod);
	CCString l_smiledMod = smiledMod;

	_bstr_t smiledModBump(pSmiledModBump); 
	CCString l_smiledModBump = smiledModBump;	

	_bstr_t hyperCubeCorrel (pHyperCubeCorrel);
	CCString l_hyperCubeCorrel = hyperCubeCorrel; 

	long payFreqPhase1Id;	
	if((payFreqPhase1Id = ARM_ConvFrequency(l_payFreqPhase1, C_result)) == ARM_DEFAULT_ERR)
	{
		ERROR_MSG("Invalid Frq",pRet,ARM_ERROR_FREQ);
		return S_OK;
	}

	long resetFreqPhase1Id;
	if((resetFreqPhase1Id = ARM_ConvFrequency(l_resetFreqPhase1, C_result)) == ARM_DEFAULT_ERR)
	{
		ERROR_MSG("Invalid Frq",pRet,ARM_ERROR_FREQ);
		return S_OK;
	}

	long payFreqPhase2Id;
	if((payFreqPhase2Id = ARM_ConvFrequency(l_payFreqPhase2, C_result)) == ARM_DEFAULT_ERR)
	{
		ERROR_MSG("Invalid Frq",pRet,ARM_ERROR_FREQ);
		return S_OK;
	}

	long resetFreqPhase2Id;
	if((resetFreqPhase2Id = ARM_ConvFrequency(l_resetFreqPhase2, C_result)) == ARM_DEFAULT_ERR)
	{
		ERROR_MSG("Invalid Frq",pRet,ARM_ERROR_FREQ);
		return S_OK;
	}

	long payFreqPhase3Id;
	if((payFreqPhase3Id = ARM_ConvFrequency(l_payFreqPhase3, C_result)) == ARM_DEFAULT_ERR)
	{
		ERROR_MSG("Invalid Frq",pRet,ARM_ERROR_FREQ);
		return S_OK;
	}

	long resetFreqPhase3Id;
	if((resetFreqPhase3Id = ARM_ConvFrequency(l_resetFreqPhase3, C_result)) == ARM_DEFAULT_ERR)
	{
		ERROR_MSG("Invalid Frq",pRet,ARM_ERROR_FREQ);
		return S_OK;
	}

	long payFreqFundId;
	if((payFreqFundId = ARM_ConvFrequency(l_payFreqFund, C_result)) == ARM_DEFAULT_ERR)
	{
		ERROR_MSG("Invalid Frq",pRet,ARM_ERROR_FREQ);
		return S_OK;
	}

	long resetFreqFundId;
	if((resetFreqFundId = ARM_ConvFrequency(l_resetFreqFund, C_result)) == ARM_DEFAULT_ERR)
	{
		ERROR_MSG("Invalid Frq",pRet,ARM_ERROR_FREQ);
		return S_OK;
	}

	long indexPhase1Id; 
	indexPhase1Id = ARM_ConvIrType(l_indexPhase1); 

	long indexPhase3Id;
	indexPhase3Id = ARM_ConvIrType(l_indexPhase3);
	
	long indexFundId;
	indexFundId = ARM_ConvIrType(l_indexFund);

	long rollId;
	if((rollId = ARM_ConvFwdRule(l_roll, C_result)) == ARM_DEFAULT_ERR)
	{
		ERROR_MSG("Invalid Int Roll",pRet,ARM_ERROR_FREQ);
		return S_OK;
	}	
	
	long adjPhase1Id = ARM_ConvIntRule(l_adjPhase1);
	long adjPhase2Id = ARM_ConvIntRule(l_adjPhase2);
	long adjPhase2DIGId = ARM_ConvIntRule(l_adjPhase2DIG);
	long adjPhase3Id = ARM_ConvIntRule(l_adjPhase3);
	long adjFundId = ARM_ConvIntRule(l_adjFund);

	long resetTimingPhase1Id = ARM_ConvPayResetRule (l_resetTimingPhase1);
	long resetTimingPhase2DIGId = ARM_ConvPayResetRule (l_resetTimingPhase2DIG);
	long resetTimingPhase3Id = ARM_ConvPayResetRule (l_resetTimingPhase3);
	long resetTimingFundId = ARM_ConvPayResetRule (l_resetTimingFund);

	CCString l_payTimingPhase2DIG = "ARR"; 
	long payTimingPhase2DIGId = ARM_ConvPayResetRule(l_resetTimingPhase2DIG);

	long dayCountPhase1Id= ARM_ConvDayCount(l_dayCountPhase1);
	//long dayCountPhase2DIGId= ARM_ConvDayCount(l_dayCountPhase2DIG);
	long dayCountPhase2Id= ARM_ConvDayCount(l_dayCountPhase2);
	long dayCountPhase3Id= ARM_ConvDayCount(l_dayCountPhase3);
	long dayCountFundId= ARM_ConvDayCount(l_dayCountFund);

	long stubId = ARM_ConvStubRule(l_stub);
	long ccyId; 

	long resetTimingPhase2Id = 1;

	if (ARMLOCAL_ISOCCY(l_ccy,C_result) == ARM_OK)
	{
		ccyId = C_result.getLong();
	}
	else
	{
		ERROR_MSG(C_result.getMsg(),pRet);
		return S_OK;
	}

	//************************************************************************************************
	//*********************************** Création Amort Refvalue ************************************
	//************************************************************************************************
	long amortRefValue;
	long amortRefValueToDestroy = ARM_NULL_OBJECT;

	if (LocalGetNumObjectId(l_typeAmort) == ARM_KO)
	{
		double endDateAmort = pEndDateAmort;
		double txAmort = pTxAmort;
		double amountAmort = pAmountAmort;

		_bstr_t intRuleAmort(pIntRuleAmort);
		CCString l_intRuleAmort =  intRuleAmort;

		_bstr_t freqAmort(pFreqAmort);
		CCString l_freqAmort = freqAmort;

		_bstr_t dayCountAmort(pDayCountAmort);
		CCString l_dayCountAmort = dayCountAmort;

		long freqAmortId; 
		if((freqAmortId = ARM_ConvFrequency(l_freqAmort, C_result)) == ARM_DEFAULT_ERR)
		{
			ERROR_MSG("Invalid Frq",pRet,ARM_ERROR_FREQ);
			return S_OK;
		}

		long intRuleAmortId= ARM_ConvIntRule(l_intRuleAmort);
		long dayCountAmortId = ARM_ConvDayCount(l_dayCountAmort);	

		// AmortMethodId
		long amortMethodId;

		if((amortMethodId = ARM_ConvAmortMethod (l_typeAmort, C_result)) == ARM_DEFAULT_ERR)
		{
			ERROR_MSG("Invalid TypeAmort",pRet,ARM_ERROR_FREQ);
			return S_OK;
		}

		long FixedLegAmortId;
		if (ARMLOCAL_FIXEDLEG(startDatePhase1,
							  endDateAmort,
							  K_RCV,
							  0,
							  1.0,
							  dayCountAmortId,
							  freqAmortId,
							  K_COMP_PROP,
							  K_ARREARS,
							  intRuleAmortId,
							  stubId,
							  false,
							  l_ccy,
							  "NULL",
							  K_NX_NONE,
							  -1.0,
							  K_NO,
								 -1,
							  C_result) == ARM_OK) 
							  
		{
			FixedLegAmortId = C_result.getLong();
		}
		else
		{
			ERROR_MSG(C_result.getMsg(),pRet);
			return S_OK;
		}

		if (ARMLOCAL_ARM_GenAmortization(FixedLegAmortId,
										 amortMethodId,
										 freqAmortId,
										 amountAmort,
										 -1,
										 Ntl,
										 txAmort,
										 0.,
										 ARM_NULL_OBJECT,
										 0.0,
										 C_result) == ARM_OK) 
		{
			amortRefValue = C_result.getLong();
			amortRefValueToDestroy = amortRefValue;
		}
		else
		{
			ERROR_MSG(C_result.getMsg(),pRet);
			return S_OK;
		}

		retCode = ARMLOCAL_FreeObject(FixedLegAmortId,C_result);
	}
	else
	{
		amortRefValue = LocalGetNumObjectId(l_typeAmort);

		// calcul du 1er nominal
		double dPaymentDate;

		if (ARMLOCAL_DisplayRefValue(amortRefValue,true,C_result) == ARM_OK)
		{
			dPaymentDate = C_result.getArray(0);
		}
		else
		{
			ERROR_MSG(C_result.getMsg(),pRet);
			return S_OK;
		}
		
		if (ARMLOCAL_CptRefValue(amortRefValue,dPaymentDate,C_result) == ARM_OK)
		{
			Ntl = C_result.getDouble();
		}
		else
		{
			ERROR_MSG(C_result.getMsg(),pRet);
			return S_OK;
		} 
	}
	//************************************************************************************************
	//*********************************** Fin Création Amort Refvalue ************************************
	//************************************************************************************************


	//************************************************************************************************
	//*********************************** Création Funding Leg Refvalue ******************************
	//************************************************************************************************


	long irIndexFundId;
	if (ARMLOCAL_IRINDEX(KACTUAL_360,
						 payFreqFundId,
						 -1,
						 K_COMP_PROP,
						 K_MOD_FOLLOWING,
						 K_ADVANCE,
						 -2,
						 K_ARREARS,
						 10000.,
						 false,
						 l_ccy,
						 indexFundId,
						 K_COMP_PROP,
						 adjFundId,
						 resetFreqFundId,
						 C_result) == ARM_OK)
	{
		irIndexFundId = C_result.getLong();
	}
	else
	{
		ERROR_MSG(C_result.getMsg(),pRet);
		return S_OK;
	}

	long swapLegId;
	if (ARMLOCAL_SWAPLEG(irIndexFundId,
						 startDatePhase1,
						 endDatePhase3,
						 K_RCV,
						 spreadType, 
						 spreadFund, 
						 false,
						 l_ccy,
						 dayCountFundId,
						 -2,
						 "NULL",
						 "NULL",
						 0,
						 K_NX_NONE,
						 stubId,
						 -1.0,
						 K_NO,
						  -1,
						 C_result) == ARM_OK) 
	{
		swapLegId = C_result.getLong();
	}
	else
	{
		ERROR_MSG(C_result.getMsg(),pRet);
		return S_OK;
	}
	
	
	
	// ********** Funding Leg ************

	long cloneId;
	if (ARMLOCAL_ClonedAndSetNotional(swapLegId,
									  amortRefValue,
									  100.,
									  C_result) == ARM_OK)
	{
		cloneId = C_result.getLong();
	}
	else
	{
		ERROR_MSG(C_result.getMsg(),pRet);
		return S_OK;
	}

	
	//*********   Price funding Leg   *************************
	double pxFund; 
	//double pxFundDelta;

	if (ARMLOCAL_ARM_Price(cloneId,
						   LocalGetNumObjectId(l_smiledMod),
						   C_result) == ARM_OK)
	{
		pxFund = C_result.getDouble();
	}
	else
	{
		ERROR_MSG(C_result.getMsg(),pRet);
		return S_OK;
	}

	#ifdef _DEBUG
		res = ARMLOCAL_ARM_View(irIndexFundId,C_result);
	#endif

	#ifdef _DEBUG
		res = ARMLOCAL_ARM_View(swapLegId,C_result);
	#endif

	#ifdef _DEBUG
		res = ARMLOCAL_ARM_View(cloneId,C_result);
	#endif

	#ifdef _DEBUG
		res = ARMLOCAL_ARM_View(LocalGetNumObjectId(l_smiledMod),C_result);
	#endif


	// shift seulement sur l'euro !?
	//if (ARMLOCAL_ARM_Price(cloneId,
	//					   LocalGetNumObjectId(l_smiledModBump),
	//					   C_result) == ARM_OK)
	//{
	//	pxFundDelta = C_result.getDouble();
	//}
	//else
	//{
	//	ERROR_MSG(C_result.getMsg(),pRet);
	//	return S_OK;
	//}

	deltaFund = computeDeltaByPlot(pxFund, cloneId, vBumpBsGenMod);


	retCode = ARMLOCAL_FreeObject(irIndexFundId, C_result);
	retCode = ARMLOCAL_FreeObject(swapLegId, C_result);
	retCode = ARMLOCAL_FreeObject(cloneId,C_result);
	
	// Fin du calcul du prix de la funding Leg

	//************************************************************************************************
	//*********************************** Fin Création Funding Leg Refvalue **************************
	//************************************************************************************************





	
	//************************************************************************************************
	//*********************************** Création Leg Phase1 Refvalue *******************************
	//************************************************************************************************

	long irIndexPhase1Id;
	double pxPhase1; 
	//double pxPhase1Delta;
	long swapLegPhase1Id;
	long clonePhase1Id;
	double resetGapPhase1; 


	if (resetTimingPhase1Id == 1)
	{
		resetGapPhase1 = -2.; 
	}
	else
	{
		resetGapPhase1 = -15.; 
	}

	if ( startDatePhase1 != startDatePhase2 )
	{
		if (ARMLOCAL_IRINDEX(KACTUAL_360,
							 payFreqPhase1Id,
							 -1,
							 K_COMP_PROP,
							 K_MOD_FOLLOWING,
							 resetTimingPhase1Id,	
							 resetGapPhase1,
							 K_ARREARS,
							 10000.,
							 false,
							 l_ccy,
							 indexPhase1Id,
							 K_COMP_PROP,
							 adjPhase1Id,
							 resetFreqPhase1Id,
							 C_result) == ARM_OK) 
		{
			irIndexPhase1Id = C_result.getLong();
		}
		else
		{
			ERROR_MSG(C_result.getMsg(),pRet);
			return S_OK;
		}

		
		if (ARMLOCAL_SWAPLEG(irIndexPhase1Id,
							 startDatePhase1,
							 startDatePhase2,
							 K_RCV,
							 0, 
							 spreadPhase1, 
							 false,
							 l_ccy,
							 dayCountPhase1Id,
							 resetGapPhase1,
							 "NULL",
							 "NULL",
							 0,
							 K_NX_NONE,
							 stubId,
							 -1.0,
							 K_NO,
						  -1,
							 C_result) == ARM_OK) 
		{
			swapLegPhase1Id = C_result.getLong();
		}
		else
		{
			ERROR_MSG(C_result.getMsg(),pRet);
			return S_OK;
		}
		
		
		
		// ********** Phase1 Leg ************

		
		if (ARMLOCAL_ClonedAndSetNotional(swapLegPhase1Id,
										  amortRefValue,
										  100.,
										  C_result) == ARM_OK)
		{
			clonePhase1Id = C_result.getLong();
		}
		else
		{
			ERROR_MSG(C_result.getMsg(),pRet);
			return S_OK;
		}

		
		//*********   Price Phase1 *************************


		if (ARMLOCAL_ARM_Price(clonePhase1Id,
							   LocalGetNumObjectId(l_smiledMod),
							   C_result) == ARM_OK)
		{
			pxPhase1 = C_result.getDouble();
		}
		else
		{
			ERROR_MSG(C_result.getMsg(),pRet);
			return S_OK;
		}

		//// shift seulement sur l'euro !?
		//if (ARMLOCAL_ARM_Price(clonePhase1Id,
		//					   LocalGetNumObjectId(l_smiledModBump),
		//					   C_result) == ARM_OK)
		//{
		//	pxPhase1Delta = C_result.getDouble();
		//}
		//else
		//{
		//	ERROR_MSG(C_result.getMsg(),pRet);
		//	return S_OK;
		//}

		deltaPhase1 = computeDeltaByPlot(pxPhase1, clonePhase1Id, vBumpBsGenMod);

		retCode = ARMLOCAL_FreeObject(irIndexPhase1Id,C_result);
		retCode = ARMLOCAL_FreeObject(swapLegPhase1Id,C_result);
		retCode = ARMLOCAL_FreeObject(clonePhase1Id,C_result);
	}
	else

	{
		pxPhase1 = 0.;
		//pxPhase1Delta = 0.;

	}
	//************************************************************************************************
	//************************** Fin Création Leg Phase1 Refvalue ************************************
	//************************************************************************************************


	//************************************************************************************************
	//*********************************** Création Leg Phase2 Refvalue *******************************
	//************************************************************************************************

	long irIndexPhase2Id;
	double resetGapPhase2; 

	if (resetTimingPhase2Id == 1)
	{
		resetGapPhase2 = -2.; 
	}
	else
	{
		resetGapPhase2 = -15.; 
	}


	if (ARMLOCAL_IRINDEX(KACTUAL_360,
						 payFreqPhase2Id,
						 -1,
						 K_COMP_PROP,
						 K_MOD_FOLLOWING,
						 resetTimingPhase2Id,	
						 resetGapPhase2,
						 K_ARREARS,
						 10000.,
						 false,
						 l_ccy,
						 K_FIXED,
						 K_COMP_PROP,
						 adjPhase2Id,
						 resetFreqPhase2Id,
						 C_result) == ARM_OK) 
	{
		irIndexPhase2Id = C_result.getLong();
	}
	else
	{
		ERROR_MSG(C_result.getMsg(),pRet);
		return S_OK;
	}

	long swapLegPhase2Id;
	if (ARMLOCAL_SWAPLEG(irIndexPhase2Id,
						 startDatePhase2,
						 startDatePhase3,
						 K_RCV,
						 0, 
						 spreadPhase2, 
						 false,
						 l_ccy,
						 dayCountPhase2Id,
						 resetGapPhase2,
						 "NULL",
						 "NULL",
						 0,
						 K_NX_NONE,
						 stubId,
						 -1.0,
						 K_NO,
						  -1,
						 C_result) == ARM_OK) 
	{
		swapLegPhase2Id = C_result.getLong();
	}
	else
	{
		ERROR_MSG(C_result.getMsg(),pRet);
		return S_OK;
	}
	
	
	
	// ********** Phase2 Leg ************

	long clonePhase2Id;
	if (ARMLOCAL_ClonedAndSetNotional(swapLegPhase2Id,
									  amortRefValue,
									  100.,
									  C_result) == ARM_OK)
	{
		clonePhase2Id = C_result.getLong();
	}
	else
	{
		ERROR_MSG(C_result.getMsg(),pRet);
		return S_OK;
	}

	
	//*********   Price Phase2 *************************
	double pxPhase2;
	//double pxPhase2Delta; 

	if (ARMLOCAL_ARM_Price(clonePhase2Id,
						   LocalGetNumObjectId(l_smiledMod),
						   C_result) == ARM_OK)
	{
		pxPhase2 = C_result.getDouble();
	}
	else
	{
		ERROR_MSG(C_result.getMsg(),pRet);
		return S_OK;
	}


	#ifdef _DEBUG
		res = ARMLOCAL_ARM_View(irIndexPhase2Id,C_result);
	#endif

	#ifdef _DEBUG
		res = ARMLOCAL_ARM_View(swapLegPhase2Id,C_result);
	#endif

	#ifdef _DEBUG
		res = ARMLOCAL_ARM_View(clonePhase2Id,C_result);
	#endif

	#ifdef _DEBUG
		res = ARMLOCAL_ARM_View(LocalGetNumObjectId(l_smiledMod),C_result);
	#endif

	//// shift seulement sur l'euro !?
	//if (ARMLOCAL_ARM_Price(clonePhase2Id,
	//					   LocalGetNumObjectId(l_smiledModBump),
	//					   C_result) == ARM_OK)
	//{
	//	pxPhase2Delta = C_result.getDouble();
	//}
	//else
	//{
	//	ERROR_MSG(C_result.getMsg(),pRet);
	//	return S_OK;
	//}

	deltaPhase2 = computeDeltaByPlot(pxPhase2, clonePhase2Id, vBumpBsGenMod);

	retCode = ARMLOCAL_FreeObject(irIndexPhase2Id,C_result);
	retCode = ARMLOCAL_FreeObject(swapLegPhase2Id,C_result);
	retCode = ARMLOCAL_FreeObject(clonePhase2Id,C_result);

	//************************************************************************************************
	//************************** Fin Création Leg Phase2 Refvalue ************************************
	//************************************************************************************************



	//************************************************************************************************
	//*************************** Création Leg Phase3 Refvalue ***************************************
	//************************************************************************************************

	long irIndexPhase3Id;
	double pxPhase3; 
	//double pxPhase3Delta;
	long swapLegPhase3Id;
	long clonePhase3Id;
	double resetGapPhase3; 

	if ( startDatePhase3 != endDatePhase3 )
	{		
		if (resetTimingPhase3Id == 1)
		{
			resetGapPhase3 = -2.; 
		}
		else
		{
			resetGapPhase3 = -15.; 
		}


		if (ARMLOCAL_IRINDEX(KACTUAL_360,
							 payFreqPhase3Id,
							 -1,
							 K_COMP_PROP,
							 K_MOD_FOLLOWING,
							 resetTimingPhase3Id,	
							 resetGapPhase3,
							 K_ARREARS,
							 10000.,
							 false,
							 l_ccy,
							 indexPhase3Id,
							 K_COMP_PROP,
							 adjPhase3Id,
							 resetFreqPhase3Id,
							 C_result) == ARM_OK) 
		{
			irIndexPhase3Id = C_result.getLong();
		}
		else
		{
			ERROR_MSG(C_result.getMsg(),pRet);
			return S_OK;
		}

		

		if (ARMLOCAL_SWAPLEG(irIndexPhase3Id,
							 startDatePhase3,
							 endDatePhase3,
							 K_RCV,
							 0, 
							 spreadPhase3, 
							 false,
							 l_ccy,
							 dayCountPhase3Id,
							 resetGapPhase3,
							 "NULL",
							 "NULL",
							 0,
							 K_NX_NONE,
							 stubId,
							 -1.0,
							 K_NO,
						  -1,
							 C_result) == ARM_OK) 
		{
			swapLegPhase3Id = C_result.getLong();
		}
		else
		{
			ERROR_MSG(C_result.getMsg(),pRet);
			return S_OK;
		}
		
		
		
		// ********** Phase3 Leg ************

		
		if (ARMLOCAL_ClonedAndSetNotional(swapLegPhase3Id,
										  amortRefValue,
										  100.,
										  C_result) == ARM_OK)
		{
			clonePhase3Id = C_result.getLong();
		}
		else
		{
			ERROR_MSG(C_result.getMsg(),pRet);
			return S_OK;
		}


		//*********   Price Phase3 *************************	
		if (ARMLOCAL_ARM_Price(clonePhase3Id,
							   LocalGetNumObjectId(l_smiledMod),
							   C_result) == ARM_OK)
		{
			pxPhase3 = C_result.getDouble();
		}
		else
		{
			ERROR_MSG(C_result.getMsg(),pRet);
			return S_OK;
		}

		//// shift seulement sur l'euro !?
		//if (ARMLOCAL_ARM_Price(clonePhase3Id,
		//					   LocalGetNumObjectId(l_smiledModBump),
		//					   C_result) == ARM_OK)
		//{
		//	pxPhase3Delta = C_result.getDouble();
		//}
		//else
		//{
		//	ERROR_MSG(C_result.getMsg(),pRet);
		//	return S_OK;
		//}

		deltaPhase3 = computeDeltaByPlot(pxPhase3, clonePhase3Id, vBumpBsGenMod);

		retCode = ARMLOCAL_FreeObject(irIndexPhase3Id,C_result);
		retCode = ARMLOCAL_FreeObject(swapLegPhase3Id,C_result);
		retCode = ARMLOCAL_FreeObject(clonePhase3Id,C_result);
	}
	else
	{
		pxPhase3 = 0.;
		//pxPhase3Delta = 0.; 
	}


	sizeMatu = vBumpBsGenMod.size(); 
	long i; 

	for (i =0; i< sizeMatu; i++)
	{
		deltaOther.push_back(-deltaFund[i]); 
	}

	for (i =0; i< sizeMatu; i++)
	{
		deltaOther[i] = deltaOther[i] + deltaPhase2[i]; 
	}

	if ( startDatePhase1 != startDatePhase2 )
	{
		for (i =0; i< sizeMatu; i++)
		{
			deltaOther[i] = deltaOther[i] + deltaPhase1[i]; 
		}
	}

	if (startDatePhase3 != endDatePhase3 )
	{
		for (i =0; i< sizeMatu; i++)
		{
			deltaOther[i] = deltaOther[i] + deltaPhase3[i]; 
		}
	}

	double pxPhase2DIG; 
	double pxTotal; 
	double pxDelta; 
	double SensiTotal; 


	
	double prixInter;
	double Xincr;
	
	double resetGapPhase2DIG; 

	if (resetTimingPhase2DIGId == 1)
	{
		resetGapPhase2DIG = -2.; 
	}
	else
	{
		resetGapPhase2DIG = -15.; 
	}

//******************** Calcul premiere PV ******************
	int iteration =1;
	double prix;
	double prixAct; 
	double prixAct2;
//	double TauxFixe = 10.01;
	double SprSeek = 10.01;
	double SprSeekFinal= 10.01;
	double SprSeekBonifie = 10.01;  // PUREMENT ARBITRAIRE
	fee = Ntl * fee / 100 ;
//******************** Calcul premiere PV ******************

	
	VECTOR<double> pxPentifix = computePentifix(l_indexPhase2DIG,l_indexLongPhase2DIG, startDatePhase2, startDatePhase3, 
			l_floorOrCap, l_strikePhase2DIG, SprSeekBonifie, ccyId, 
			dayCountPhase2Id, 
			resetTimingPhase2DIGId, payTimingPhase2DIGId,  resetGapPhase2DIG, resetFreqPhase2Id, payFreqPhase2Id, 
			adjPhase2DIGId, stubId,amortRefValue, 
			deltaOther, vBumpBsGenMod, vBumpVolBsGenMod,
			l_hyperCubeCorrel, l_volCurvFromMatriceShift,
			l_vol, l_volCub, l_corrManager,
			l_convexityManager, l_zc, l_smiledMod);


	if (pxPentifix.size() == 0)
	{
		ERROR_MSG("Pb in computePentifix",pRet);
		return S_OK;
	}

	pxPhase2DIG = pxPentifix[0];
	double pxPhase2DIGBumpCorrel = pxPentifix[1];
	double vega = pxPentifix[2]; 
	double deltaTotal = pxPentifix[3];
	double MT = pxPentifix[4];


	pxTotal = pxPhase1 + pxPhase2 + pxPhase3 - pxFund + pxPhase2DIG;

	//pxDelta = fabs(pxPhase2DIGBumpCorrel-pxPhase2DIG) + MT; 

	if ((MT - (fabs(vega) + fabs(deltaTotal) + fabs(pxPhase2DIGBumpCorrel-pxPhase2DIG))) > 0) 
	{
		pxDelta =MT ; 
	}
	else
	{
		pxDelta = fabs(pxPhase2DIGBumpCorrel-pxPhase2DIG) + fabs(vega) + fabs(deltaTotal); 

	}

	if (pxDelta <0.) pxDelta = - pxDelta;


	SensiTotal = pxDelta;
	prixAct = pxTotal - SensiTotal - fee; 


	prix = prixAct ;
	prixAct2 =1;
	
//********************                    ******************
//******************** BOUCLE SUR LE PRICE ******************
	while ((abs(prixAct2) > 0.0001) && (iteration < 10 ))
	{
		SprSeek = SprSeekBonifie;
		
		//*************************
		// **** Premier recalcul ****
		//****************************

	pxPentifix = computePentifix(l_indexPhase2DIG,l_indexLongPhase2DIG, startDatePhase2, startDatePhase3, 
			l_floorOrCap, l_strikePhase2DIG, SprSeekBonifie , ccyId, 
			dayCountPhase2Id, 
			resetTimingPhase2DIGId, payTimingPhase2DIGId,  resetGapPhase2DIG, resetFreqPhase2Id, payFreqPhase2Id, 
			adjPhase2DIGId, stubId, amortRefValue,
			deltaOther, vBumpBsGenMod, vBumpVolBsGenMod,
			l_hyperCubeCorrel, l_volCurvFromMatriceShift,
			l_vol, l_volCub, l_corrManager,
			l_convexityManager, l_zc, l_smiledMod);

		if (pxPentifix.size() == 0)
		{
			ERROR_MSG("Pb in computePentifix",pRet);
			return S_OK;
		}

		pxPhase2DIG = pxPentifix[0];
		pxPhase2DIGBumpCorrel = pxPentifix[1];
		vega = pxPentifix[2]; 
		deltaTotal = pxPentifix[3];
		MT = pxPentifix[4];

		pxTotal = pxPhase1 + pxPhase2 + pxPhase3 - pxFund + pxPhase2DIG;

		if ((MT - (fabs(vega) + fabs(deltaTotal) + fabs(pxPhase2DIGBumpCorrel-pxPhase2DIG))) > 0) 
		{
			pxDelta =MT ; 
		}
		else
		{
			pxDelta = fabs(pxPhase2DIGBumpCorrel-pxPhase2DIG) + fabs(vega) + fabs(deltaTotal); 

		}


		if (pxDelta <0.) pxDelta = - pxDelta;


		SensiTotal = pxDelta;
		prixAct = pxTotal - SensiTotal - fee; 

		// les deltas se compensent et on somme en valeur absolue les deltas et vegas

		//******************************
		// **** affectation valeurs ****
		//******************************

		prixAct2 = prixAct ;
		Xincr = SprSeek / 1000000;
		SprSeekBonifie = SprSeek + Xincr;		

		
		//****************************
		// **** Deuxieme recalcul ****
		//****************************
		pxPentifix = computePentifix(l_indexPhase2DIG,l_indexLongPhase2DIG, startDatePhase2, startDatePhase3, 
									l_floorOrCap, l_strikePhase2DIG, SprSeekBonifie , ccyId, 
									dayCountPhase2Id, 
			resetTimingPhase2DIGId, payTimingPhase2DIGId,  resetGapPhase2DIG, resetFreqPhase2Id, payFreqPhase2Id, 
			adjPhase2DIGId, stubId, amortRefValue, 
			deltaOther, vBumpBsGenMod, vBumpVolBsGenMod, 
			l_hyperCubeCorrel, l_volCurvFromMatriceShift, 
			l_vol, l_volCub, l_corrManager, 
			l_convexityManager, l_zc, l_smiledMod);


		if (pxPentifix.size() == 0)
		{
			ERROR_MSG("Pb in computePentifix",pRet);
			return S_OK;
		}

		pxPhase2DIG = pxPentifix[0];
		pxPhase2DIGBumpCorrel = pxPentifix[1];
		vega = pxPentifix[2]; 
		deltaTotal = pxPentifix[3];
		MT = pxPentifix[4];

		pxTotal = pxPhase1 + pxPhase2 + pxPhase3 - pxFund + pxPhase2DIG;

		if ((MT - (fabs(vega) + fabs(deltaTotal) + fabs(pxPhase2DIGBumpCorrel-pxPhase2DIG))) > 0) 
		{
			pxDelta =MT ; 
		}
		else
		{
			pxDelta = fabs(pxPhase2DIGBumpCorrel-pxPhase2DIG) + fabs(vega) + fabs(deltaTotal); 

		}
		
		if (pxDelta <0.) pxDelta = - pxDelta;

		SensiTotal = pxDelta;
		prixAct = pxTotal - SensiTotal - fee; 
	

		iteration++;
		
		SprSeekFinal = SprSeekBonifie;
		prixInter = ( prixAct - prixAct2 ) / Xincr;
		// if prixinter =  0 sortir erreur
		SprSeekBonifie = SprSeek - prixAct / prixInter ; 
		prixAct2 = prixAct;

		
	}
		
//**************** PRICE *******************


	retCode = ARMLOCAL_FreeObject(amortRefValueToDestroy,C_result);
	retCode = ARMLOCAL_FreeObject(ccyId,C_result);

	VECTOR<double> vRes;

	vRes.push_back(spreadPhase2 + SprSeekFinal);

	vRes.push_back(spreadPhase2 + SprSeekFinal);


	Res = VECTORDOUBLE2VARIANT(vRes, pRet);

	

	_variant_t wrap_Res;
	wrap_Res.Attach(*pRet);
	*pRet=wrap_Res.Detach();

	return S_OK;
}
catch(...)
{
	return E_FAIL;
}
}


VECTOR<double> computePentifix(const CCString& index, const CCString& indexLong, double startDate, double endDate,
							   const CCString& floorOrCap, const CCString& strike, double payOff, long ccyId, long dayCountId,
							   long resetTimingId, long payTimingId, double resetGap, long resetFreqId, long payFreqId,
							   long intRuleId, long stubRuleId, long amortRefvalue,
								VECTOR<double> vDeltaOther, const VECTOR<CCString>& vBumpBsGenMod, const VECTOR<CCString>& vBumpVolBsGenMod, 
							   const CCString& hyperCubeCorrel, const CCString& VolCurvFromMatriceShift,
							   const CCString& vol, const CCString& volCUB, const CCString& CorrManager,
							   const CCString& ConvexityManager, const CCString& zc, const CCString& smiledMod)
{
	VECTOR<double> resultat;
	resultat.clear();

try
{
	double pxDigital;
	double pxBumpDigital; 
	double MT;
	double ProbaExercice; 
	long digitalId;
	long cloneDigitalId; 

	long res;

	double shiftCorrel; 
	
	long bumpHyperCubeCorrelId;
	ARM_result C_result;

	// Modèles
	long bsModGenId; 
	long BumpBsModGenId;

	VECTOR <double> spreadVector;
	VECTOR <double> weightVector;

   	VECTOR<double> C_fixing1_double_vect;
   	VECTOR<double> C_fixing2_double_vect;

	long slopeFlag; 
	int cptStrikeMethod;
	int computedFormula;

	slopeFlag = 1; 
	cptStrikeMethod = 1.; 
	computedFormula = 1.; 

	// Ss-Jacent
	long ssJacentId;
	long cloneSsJacentId;
	double pxSsJacent;

	VECTOR <double> deltaDigital;
	
	double d_index; 
	long indexId; 
	indexId = ARM_ConvIrType(index); 

	if (IsCMSIndex((ARM_INDEX_TYPE)indexId) == 1)
	{
		d_index = (indexId - K_CMS1 +1 ); 
	}	
	else
	{
		d_index = 1./FromLiborTypeToFrequency(indexId); 
	}
	
	double d_indexLong; 
	long indexLongId; 
	indexLongId = ARM_ConvIrType(indexLong); 

	if (IsCMSIndex((ARM_INDEX_TYPE)indexLongId) ==1)
	{
		d_indexLong = indexLongId - K_CMS1 + 1; 
	}
	else
	{
		d_indexLong = 1./FromLiborTypeToFrequency(indexLongId); 
	}

	//Définition du sreadVector et du WeightVector
	spreadVector.push_back(-0.01);
	spreadVector.push_back(0.01); 

	weightVector.push_back(1);
	weightVector.push_back(1);

	long capOrFloorId; 

	if((capOrFloorId = ARM_ConvCapOrFloor(floorOrCap, C_result)) == ARM_DEFAULT_ERR)
	{
		return resultat;
	}
		
	long modelTypeId;
	long spreadVolTypeId;
	CCString modelType; 
	CCString spreadVolType;
	bool isLnVol=true;

	modelType = "2LOG";
	spreadVolType = "C"; 

	if ((modelTypeId = ARM_ConvModelType (modelType, C_result)) == ARM_DEFAULT_ERR)
	{
		return resultat;
	}

	if ((spreadVolTypeId = ARM_ConvVolType2 (spreadVolType, C_result)) == ARM_DEFAULT_ERR)
	{
		return resultat;
	}

	VECTOR<double> C_rate; 
	VECTOR<CCString> C_matu;
	double C_meth = 0;	

	//************************************************************************************************
	//*************************** Création Digitale **************************************************
	//************************************************************************************************

	vector<double> calibInfos;

	//Calcul du prix de la Digitale
	if (ARMLOCAL_SPREADDIGITAL(startDate,
						 endDate,
						 capOrFloorId,
						 1, 
						 LocalGetNumObjectId(strike), 
						 0, 
						 payOff,
						 indexId, 
						 indexLongId,
						 weightVector[0], 
						 weightVector[1],
						 dayCountId,
						 resetFreqId, 
						 payFreqId, 
						 resetTimingId, 
						 payTimingId, 					
						 ccyId,
						 resetGap, 
						 spreadVector[0], 
						 spreadVector[1],
						 1,
						 ARM_NULL_OBJECT,
						 C_fixing1_double_vect,
						 1,
						 ARM_NULL_OBJECT,
						 C_fixing2_double_vect,
						 intRuleId, 
						 stubRuleId, 
						 slopeFlag,
						 cptStrikeMethod, 
						 computedFormula,
						 calibInfos,
						 C_result) == ARM_OK)
	{
		digitalId = C_result.getLong();
	}
	else
	{
		return resultat;
	}

	if (ARMLOCAL_ClonedAndSetNotional(digitalId,
									  amortRefvalue,
									  100.,
									  C_result) == ARM_OK)
	{
		cloneDigitalId = C_result.getLong();
	}
	else
	{
		return resultat;
	}

	if (ARMLOCAL_BSMODELGEN(LocalGetNumObjectId(zc), ARM_NULL_OBJECT, LocalGetNumObjectId(vol), LocalGetNumObjectId(volCUB), 
				LocalGetNumObjectId(CorrManager), LocalGetNumObjectId(ConvexityManager), 
				LocalGetNumObjectId(zc), LocalGetNumObjectId(hyperCubeCorrel), ARM_NULL_OBJECT, ARM_NULL_OBJECT,
                    modelTypeId, spreadVolTypeId, LocalGetNumObjectId(smiledMod), isLnVol, 100, C_result) == ARM_OK)
	{
		bsModGenId = C_result.getLong();
	}
	else
	{
		return resultat;
	}

	if (ARMLOCAL_ARM_Price(cloneDigitalId,
						   bsModGenId,
						   C_result) == ARM_OK)
	{
		pxDigital = C_result.getDouble();
	}
	else
	{
		return resultat;
	}

	#ifdef _DEBUG
		res = ARMLOCAL_ARM_View(digitalId,C_result);
	#endif


	//************************************************************************************************
	//************************************ Calcul Sensi **********************************************
	//************************************************************************************************
	
	deltaDigital = computeDeltaByPlot(pxDigital, cloneDigitalId, vBumpBsGenMod);

	//************************************************************************************************
	//********************************** Fin Calcul Sensi ********************************************
	//************************************************************************************************

	//  Calcul de la marge Correl

		
	VECTOR <double> vegaDigital; 
	VECTOR <double> vegaGlobal; 
	double vega = 0.0;

	long sizeVega = vBumpVolBsGenMod.size(); 
	vegaDigital= computeDeltaByPlot(pxDigital, cloneDigitalId, vBumpVolBsGenMod);

	for (int i = 0; i < sizeVega; i++)	
	{
		vegaGlobal.push_back(vegaDigital[i]);
		vega += fabs(vegaGlobal[i]);
	}

	//************************************************************************************************
	//***********************************  Calcul Delta   ********************************************
	//************************************************************************************************
	
	long sizeMatu = vDeltaOther.size(); 
	double deltaPositif= 0. ; 
	double  deltaNegatif = 0.; 

	for (i = 0; i < sizeMatu; i++)	
	{
		if ((vDeltaOther[i] + deltaDigital[i]) > 0.)
		{
			deltaPositif = deltaPositif + (vDeltaOther[i] + deltaDigital[i]);
		}
		else 
		{
			deltaNegatif = deltaNegatif + abs(vDeltaOther[i] + deltaDigital[i]);
		}
	}

	double deltaGlobal; 
	if (deltaPositif > deltaNegatif)
	{
		deltaGlobal = deltaPositif; 
	}
	else
	{
		deltaGlobal = deltaNegatif; 
	}

	//************************************************************************************************
	//********************************* Fin Calcul Delta *********************************************
	//************************************************************************************************


	//************************************************************************************************
	//********************************* Fin Calcul Vega *********************************************
	//************************************************************************************************


	//  Calcul de la marge Correl
	if (ARMLOCAL_ARM_BumpVolatility(LocalGetNumObjectId(hyperCubeCorrel),1., 0, 0, K_NO, K_YES,C_result) == ARM_OK)
	{
		bumpHyperCubeCorrelId = C_result.getLong();
	}
	else
	{
		return resultat;
	}

	if (ARMLOCAL_BSMODELGEN(LocalGetNumObjectId(zc), ARM_NULL_OBJECT, LocalGetNumObjectId(vol), LocalGetNumObjectId(volCUB), 
				LocalGetNumObjectId(CorrManager), LocalGetNumObjectId(ConvexityManager), 
				LocalGetNumObjectId(zc), bumpHyperCubeCorrelId, ARM_NULL_OBJECT, ARM_NULL_OBJECT,
                modelTypeId, spreadVolTypeId, LocalGetNumObjectId(smiledMod), isLnVol, 100, C_result) == ARM_OK)
	{
		BumpBsModGenId = C_result.getLong();
	}
	else
	{
		return resultat;
	}

	if (ARMLOCAL_ARM_Price(cloneDigitalId,
						   BumpBsModGenId,
						   C_result) == ARM_OK)
	{
		pxBumpDigital = C_result.getDouble();
	}
	else
	{
		return resultat;
	}

	long retCode = ARMLOCAL_FreeObject(bumpHyperCubeCorrelId,C_result);
	retCode = ARMLOCAL_FreeObject(BumpBsModGenId,C_result);

	if ((pxBumpDigital - pxDigital) > 0 )
	{				
		if (ARMLOCAL_ComputeVolatility(LocalGetNumObjectId(VolCurvFromMatriceShift), d_index, d_indexLong, 0, C_result) == ARM_OK)
		{
			shiftCorrel = C_result.getDouble();
		}
		else
		{
			return resultat;
		}

		/*
		
		Exemple d'utilisation de la nouvelle fonction pour récupérer la correl
		CCString i;
		CCString c;
		CCString cur;
		CCString v;
		CCString m;
		CCString ih;
		CCString iid;
		CCString plot;
		CCString indextest;
		double d;
		double yt;

		i="EURIB";
		c="EUR";
		cur="MARGIN";
		d=39216;
		v="IRG";
		m="ATM";
		ih="IRFWDVOL";
		iid="";
		long indexId_= LocalGetNumObjectId(iid);
		indextest = "CMS2";

		long getvolfromsummit;

		if (ARMLOCAL_GetVolFromSummit(i,c,cur,d,v,m,ih,indexId_,C_result) == ARM_OK)
		{
			getvolfromsummit = C_result.getLong();
		}

		if (ARMLOCAL_ConvIndexInYearTerm(indextest,d,c,C_result) == ARM_OK)
		{
			yt = C_result.getDouble();
		}

		if (ARMLOCAL_ComputeVolatility(getvolfromsummit,yt,d_index, 0, C_result) == ARM_OK)
		{
			shiftCorrel = C_result.getDouble();
		}*/
	}		
	else
	{
		if (ARMLOCAL_ComputeVolatility(LocalGetNumObjectId(VolCurvFromMatriceShift),d_indexLong, d_index, 0, C_result) == ARM_OK)
		{
			shiftCorrel = C_result.getDouble();
		}
		else
		{
			return resultat;
		}
	}


	// Bump de l'hyperCube 
	//C_rate.clear();
	//C_rate.push_back (shiftCorrel);
	if (ARMLOCAL_ARM_BumpVolatility(LocalGetNumObjectId(hyperCubeCorrel),shiftCorrel, 0, 0, K_NO, K_YES,C_result) == ARM_OK)
	{
		bumpHyperCubeCorrelId = C_result.getLong();
	}
	else
	{
		return resultat;
	}
	
	if (ARMLOCAL_BSMODELGEN(LocalGetNumObjectId(zc), ARM_NULL_OBJECT, LocalGetNumObjectId(vol), LocalGetNumObjectId(volCUB), 
				LocalGetNumObjectId(CorrManager), LocalGetNumObjectId(ConvexityManager), 
				LocalGetNumObjectId(zc), bumpHyperCubeCorrelId, ARM_NULL_OBJECT, ARM_NULL_OBJECT,
                modelTypeId, spreadVolTypeId, LocalGetNumObjectId(smiledMod), isLnVol, 100, C_result) == ARM_OK)
	{
		BumpBsModGenId = C_result.getLong();
	}
	else
	{
		return resultat;
	}

	if (ARMLOCAL_ARM_Price(cloneDigitalId,
						   BumpBsModGenId,
						   C_result) == ARM_OK)
	{
		pxBumpDigital = C_result.getDouble();
	}
	else
	{
		return resultat;
	}
	//fin calcul marge correl

	double strikeSsJacent; 
	//calcul de la discontinuité
	if (floorOrCap == "C" )
	{
		strikeSsJacent = -1000.;
	}
	else 
	{
		strikeSsJacent = 1000.;
	}

	if (ARMLOCAL_SPREADDIGITAL(startDate,
							   endDate,
							   capOrFloorId,
							   0, 
							   strikeSsJacent, 
							   0, 
							   payOff,
							   indexId, 
							   indexLongId,
							   weightVector[0], 
							   weightVector[1],
							   dayCountId,
							   resetFreqId,
							   payFreqId, 
							   resetTimingId,
							   payTimingId, 
							   ccyId,
							   resetGap, 
							   spreadVector[0], 
							   spreadVector[1],
							   1,
							   ARM_NULL_OBJECT,
							   C_fixing1_double_vect,
							   1,
							   ARM_NULL_OBJECT,
							   C_fixing2_double_vect,
							   intRuleId,
							   stubRuleId,
							   slopeFlag,
							   cptStrikeMethod, 
							   computedFormula,
							   calibInfos,
							   C_result) == ARM_OK)
	{
		ssJacentId = C_result.getLong();
	}
	else
	{
		return resultat;
	}

	if (ARMLOCAL_ClonedAndSetNotional(ssJacentId,
									  amortRefvalue,
									  100.,
									  C_result) == ARM_OK)
	{
		cloneSsJacentId = C_result.getLong();
	}
	else
	{
		return resultat;
	}

	if (ARMLOCAL_ARM_Price(cloneSsJacentId,
						   bsModGenId,
						   C_result) == ARM_OK)
	{
		pxSsJacent = C_result.getDouble();
	}
	else
	{
		return resultat;
	}

	//Fin de la récup du prix du ssJacent

	//Calcul de la marge technique
	double minMT;
	double maxMT;
	
	ProbaExercice = pxDigital/pxSsJacent;

	if (ProbaExercice > (1- ProbaExercice))
	{
		minMT = 1. - ProbaExercice; 
	}
	else
	{
		minMT = ProbaExercice;
	}

	if (minMT - 0.1 > 0)
	{
		maxMT = minMT - 0.1; 
	}
	else
	{
		maxMT = 0.;
	}

	MT = (0.015 + 0.075 *maxMT)*pxSsJacent;

	//************************************************************************************************
	//*************************** Fin création Digitale **********************************************
	//************************************************************************************************


	//*********   Price Première Phase structure 1  TD*************************
	VECTOR<double> vRes;
	
	vRes.push_back(pxDigital);
	vRes.push_back(pxBumpDigital);
	vRes.push_back(vega);
	vRes.push_back(deltaGlobal);
	vRes.push_back(MT);

	retCode = ARMLOCAL_FreeObject(digitalId,C_result);
	retCode = ARMLOCAL_FreeObject(cloneDigitalId,C_result);
	retCode = ARMLOCAL_FreeObject(bsModGenId,C_result);
	retCode = ARMLOCAL_FreeObject(BumpBsModGenId,C_result);
	retCode = ARMLOCAL_FreeObject(bumpHyperCubeCorrelId,C_result);
	retCode = ARMLOCAL_FreeObject(ssJacentId,C_result);
	retCode = ARMLOCAL_FreeObject(cloneSsJacentId,C_result);
	
	return vRes;
}
catch(...)
{
	return resultat;
}
}


VECTOR<double> computePentibor(const CCString& index, const CCString& indexLong, long indexPayId, double startDate, double endDate,
							   const CCString& floorOrCap, const CCString& strike, double payOffSpread, long ccyId, long dayCountId,
							   long resetTimingId, long payTimingId, double resetGap, long resetFreqId, long payFreqId,
							   long intRuleId, long stubRuleId, long amortRefvalue,
							   const CCString& hyperCubeCorrel, const CCString& VolCurvFromMatriceShift,
							   const CCString& vol, const CCString& volCUB, 
							   const CCString& ConvexityManager, const CCString& zc, const CCString& smiledMod, const CCString& indexIndexCorrelCube, 
							const CCString& corrEUR, const CCString& InterCorr);

STDMETHODIMP ActiveXModule::ARMcomputePentibor(double pNtl, 
										   double pStartDatePhase1,  
										   BSTR pCcy,
										   BSTR pIndexPay, 					
										   BSTR pIndexPhase1, 										   
										   double pSpreadPhase1,
										   BSTR pDayCountPhase1,
										   BSTR pPayFreqPhase1,
										   BSTR pResetFreqPhase1,
										   BSTR pResetTimingPhase1,
										   BSTR pRoll, 
										   BSTR pAdjPhase1, 
										   BSTR pStub,
										   BSTR pIndexPhase2DIG,
										   BSTR pIndexLongPhase2DIG,
										   BSTR pStrikePhase2DIG,
										   BSTR pResetTimingPhase2DIG,
										   BSTR pAdjPhase2DIG, 
										   double pStartDatePhase2,
										   double pSpreadPhase2,										   
										   BSTR pDayCountPhase2, 
										   BSTR pPayFreqPhase2, 
										   BSTR pResetFreqPhase2,
										   BSTR pAdjPhase2,
										   double pStartDatePhase3,
										   double pEndDatePhase3,
										   BSTR pIndexPhase3,
										   double pSpreadPhase3, 
										   BSTR pDayCountPhase3, 
										   BSTR pPayFreqPhase3, 
										   BSTR pResetFreqPhase3,
										   BSTR pResetTimingPhase3,
										   BSTR pAdjPhase3,
										   BSTR pIndexFund,
										   VARIANT pSpreadFund,
										   BSTR pDayCountFund, 
										   BSTR pPayFreqFund, 
										   BSTR pResetFreqFund,
										   BSTR pResetTimingFund,
										   BSTR pAdjFund, 
										   double pEndDateAmort,
										   BSTR pDayCountAmort,
										   BSTR pIntRuleAmort, 
										   double pTxAmort,
										   BSTR pFreqAmort,
										   double pAmountAmort, 
										   BSTR pTypeAmort, 
										   double pFee, 
										   BSTR pVolCurvFromMatriceShift,
										   BSTR pVol, 
										   BSTR pVolCub,  
										   BSTR pConvexityManager, 
										   BSTR pZc, 
										   BSTR pSmiledMod, 
										   BSTR pSmiledModBump, 					
										   BSTR pHyperCubeCorrel,
										   BSTR pIndexIndexCorrelCube, 
										   BSTR pCorrEUR, 
										   BSTR pInterCorr, 
										   VARIANT *pRet
										   )

{	
	// TODO: Add your implementation code here
try
{									   
	ARM_result C_result;
	double startDatePhase1 = pStartDatePhase1;  

	double spreadPhase1 = pSpreadPhase1;
	
	double startDatePhase2 = pStartDatePhase2;
	double spreadPhase2 = pSpreadPhase2;

	double startDatePhase3 = pStartDatePhase3;
	double endDatePhase3 = pEndDatePhase3;
	double spreadPhase3 = pSpreadPhase3;

	double spreadFund;
	long   spreadType;

	if (VARIANT2Double(pSpreadFund,&spreadFund) == S_FALSE)
	{
		_bstr_t bSpreadFund(pSpreadFund);
		CCString l_SpreadFund= bSpreadFund;
		spreadFund = (double) LocalGetNumObjectId(l_SpreadFund);
		spreadType = 1;
	}
	else
	{
	   spreadType = 0;
	}

	double Ntl = pNtl;
	double fee = pFee; 

	long retCode;
	double Res = 0.;

	_bstr_t ccy(pCcy);
	CCString l_ccy = ccy;

	_bstr_t indexPay(pIndexPay);
	CCString l_indexPay = indexPay;

	_bstr_t indexPhase1(pIndexPhase1);
	CCString l_indexPhase1 = indexPhase1;

	_bstr_t dayCountPhase1(pDayCountPhase1);
	CCString l_dayCountPhase1 = dayCountPhase1;

	_bstr_t resetFreqPhase1(pResetFreqPhase1);
	CCString l_resetFreqPhase1 = resetFreqPhase1;

	_bstr_t payFreqPhase1(pPayFreqPhase1);
	CCString l_payFreqPhase1 = payFreqPhase1;

	_bstr_t resetTimingPhase1(pResetTimingPhase1);
	CCString l_resetTimingPhase1 = resetTimingPhase1;

	_bstr_t roll(pRoll);
	CCString l_roll = roll;

	_bstr_t adjPhase1(pAdjPhase1);
	CCString l_adjPhase1 = adjPhase1;
	
	_bstr_t stub(pStub);
	CCString l_stub= stub;

	_bstr_t indexPhase2DIG(pIndexPhase2DIG);
	CCString l_indexPhase2DIG = indexPhase2DIG;

	_bstr_t indexLongPhase2DIG(pIndexLongPhase2DIG);
	CCString l_indexLongPhase2DIG = indexLongPhase2DIG;

	//Strike is double??
	_bstr_t strikePhase2DIG(pStrikePhase2DIG);
	CCString l_strikePhase2DIG= strikePhase2DIG;

	_bstr_t resetTimingPhase2DIG(pResetTimingPhase2DIG);
	CCString l_resetTimingPhase2DIG = resetTimingPhase2DIG;

	_bstr_t adjPhase2DIG(pAdjPhase2DIG);
	CCString l_adjPhase2DIG = adjPhase2DIG;

	_bstr_t dayCountPhase2(pDayCountPhase2);
	CCString l_dayCountPhase2= dayCountPhase2;

	_bstr_t payFreqPhase2(pPayFreqPhase2);
	CCString l_payFreqPhase2= payFreqPhase2;
	
	_bstr_t resetFreqPhase2(pResetFreqPhase2);
	CCString l_resetFreqPhase2 = resetFreqPhase2;

	_bstr_t adjPhase2(pAdjPhase2);
	CCString l_adjPhase2 = adjPhase2;

	_bstr_t indexPhase3(pIndexPhase3);
	CCString l_indexPhase3 = indexPhase3;

	_bstr_t dayCountPhase3(pDayCountPhase3);
	CCString l_dayCountPhase3 = dayCountPhase3;

	_bstr_t payFreqPhase3(pPayFreqPhase3);
	CCString l_payFreqPhase3 = payFreqPhase3;

	_bstr_t resetFreqPhase3(pResetFreqPhase3);
	CCString l_resetFreqPhase3 = resetFreqPhase3;

	_bstr_t resetTimingPhase3(pResetTimingPhase3);
	CCString l_resetTimingPhase3 = resetTimingPhase3;

	_bstr_t adjPhase3(pAdjPhase3);
	CCString l_adjPhase3 = adjPhase3;

	_bstr_t indexFund(pIndexFund);
	CCString l_indexFund = indexFund;

	_bstr_t dayCountFund(pDayCountFund);
	CCString l_dayCountFund= dayCountFund;

	_bstr_t payFreqFund(pPayFreqFund);
	CCString l_payFreqFund= payFreqFund;
	
	_bstr_t resetFreqFund(pResetFreqFund);
	CCString l_resetFreqFund = resetFreqFund;

	_bstr_t resetTimingFund(pResetTimingFund);
	CCString l_resetTimingFund = resetTimingFund;

	_bstr_t adjFund(pAdjFund);
	CCString l_adjFund= adjFund;
	
	_bstr_t typeAmort(pTypeAmort);
	CCString l_typeAmort = typeAmort;

	_bstr_t dayCountAmort(pDayCountAmort);
	CCString l_dayCountAmort = dayCountAmort;

	CCString l_floorOrCap = "C"; 

	_bstr_t volCurvFromMatriceShift(pVolCurvFromMatriceShift);
	CCString l_volCurvFromMatriceShift= volCurvFromMatriceShift;

	_bstr_t vol(pVol);
	CCString l_vol = vol;

	_bstr_t volCub(pVolCub);
	CCString l_volCub = volCub;

	_bstr_t convexityManager(pConvexityManager);
	CCString l_convexityManager = convexityManager;

	_bstr_t zc(pZc);
	CCString l_zc = zc;

	_bstr_t smiledMod (pSmiledMod);
	CCString l_smiledMod = smiledMod;

	_bstr_t smiledModBump(pSmiledModBump); 
	CCString l_smiledModBump = smiledModBump;	

	_bstr_t hyperCubeCorrel (pHyperCubeCorrel);
	CCString l_hyperCubeCorrel = hyperCubeCorrel; 

	_bstr_t indexIndexCorrelCube (pIndexIndexCorrelCube);
	CCString l_indexIndexCorrelCube = indexIndexCorrelCube;	

	_bstr_t corrEUR (pCorrEUR);
	CCString l_corrEUR = corrEUR;	

	_bstr_t interCorr (pInterCorr);
	CCString l_interCorr= interCorr;

	long payFreqPhase1Id;	
	if((payFreqPhase1Id = ARM_ConvFrequency(l_payFreqPhase1, C_result)) == ARM_DEFAULT_ERR)
	{
		ERROR_MSG("Invalid Frq",pRet,ARM_ERROR_FREQ);
		return S_OK;
	}

	long resetFreqPhase1Id;
	if((resetFreqPhase1Id = ARM_ConvFrequency(l_resetFreqPhase1, C_result)) == ARM_DEFAULT_ERR)
	{
		ERROR_MSG("Invalid Frq",pRet,ARM_ERROR_FREQ);
		return S_OK;
	}

	long payFreqPhase2Id;
	if((payFreqPhase2Id = ARM_ConvFrequency(l_payFreqPhase2, C_result)) == ARM_DEFAULT_ERR)
	{
		ERROR_MSG("Invalid Frq",pRet,ARM_ERROR_FREQ);
		return S_OK;
	}

	long resetFreqPhase2Id;
	if((resetFreqPhase2Id = ARM_ConvFrequency(l_resetFreqPhase2, C_result)) == ARM_DEFAULT_ERR)
	{
		ERROR_MSG("Invalid Frq",pRet,ARM_ERROR_FREQ);
		return S_OK;
	}

	long payFreqPhase3Id;
	if((payFreqPhase3Id = ARM_ConvFrequency(l_payFreqPhase3, C_result)) == ARM_DEFAULT_ERR)
	{
		ERROR_MSG("Invalid Frq",pRet,ARM_ERROR_FREQ);
		return S_OK;
	}

	long resetFreqPhase3Id;
	if((resetFreqPhase3Id = ARM_ConvFrequency(l_resetFreqPhase3, C_result)) == ARM_DEFAULT_ERR)
	{
		ERROR_MSG("Invalid Frq",pRet,ARM_ERROR_FREQ);
		return S_OK;
	}

	long payFreqFundId;
	if((payFreqFundId = ARM_ConvFrequency(l_payFreqFund, C_result)) == ARM_DEFAULT_ERR)
	{
		ERROR_MSG("Invalid Frq",pRet,ARM_ERROR_FREQ);
		return S_OK;
	}

	long resetFreqFundId;
	if((resetFreqFundId = ARM_ConvFrequency(l_resetFreqFund, C_result)) == ARM_DEFAULT_ERR)
	{
		ERROR_MSG("Invalid Frq",pRet,ARM_ERROR_FREQ);
		return S_OK;
	}

	long indexPhase1Id; 
	indexPhase1Id = ARM_ConvIrType(l_indexPhase1); 

	long indexPhase3Id;
	indexPhase3Id = ARM_ConvIrType(l_indexPhase3);
	
	long indexFundId;
	indexFundId = ARM_ConvIrType(l_indexFund);

	long rollId;
	if((rollId = ARM_ConvFwdRule(l_roll, C_result)) == ARM_DEFAULT_ERR)
	{
		ERROR_MSG("Invalid Int Roll",pRet,ARM_ERROR_FREQ);
		return S_OK;
	}	
	
	long adjPhase1Id = ARM_ConvIntRule(l_adjPhase1);
	long adjPhase2Id = ARM_ConvIntRule(l_adjPhase2);
	long adjPhase2DIGId = ARM_ConvIntRule(l_adjPhase2DIG);
	long adjPhase3Id = ARM_ConvIntRule(l_adjPhase3);
	long adjFundId = ARM_ConvIntRule(l_adjFund);

	long resetTimingPhase1Id = ARM_ConvPayResetRule (l_resetTimingPhase1);
	long resetTimingPhase2DIGId = ARM_ConvPayResetRule (l_resetTimingPhase2DIG);
	long resetTimingPhase3Id = ARM_ConvPayResetRule (l_resetTimingPhase3);
	long resetTimingFundId = ARM_ConvPayResetRule (l_resetTimingFund);

	CCString l_payTimingPhase2DIG = "ARR"; 
	long payTimingPhase2DIGId = ARM_ConvPayResetRule(l_resetTimingPhase2DIG);

	long dayCountPhase1Id= ARM_ConvDayCount(l_dayCountPhase1);
	long dayCountPhase2Id= ARM_ConvDayCount(l_dayCountPhase2);
	long dayCountPhase3Id= ARM_ConvDayCount(l_dayCountPhase3);
	long dayCountFundId= ARM_ConvDayCount(l_dayCountFund);

	long stubId = ARM_ConvStubRule(l_stub);
	long ccyId; 

	long resetTimingPhase2Id = 1;

	if (ARMLOCAL_ISOCCY(l_ccy,C_result) == ARM_OK)
	{
		ccyId = C_result.getLong();
	}
	else
	{
		ERROR_MSG(C_result.getMsg(),pRet);
		return S_OK;
	}

	//************************************************************************************************
	//*********************************** Création Amort Refvalue ************************************
	//************************************************************************************************
	long amortRefValue;
	long amortRefValueToDestroy = ARM_NULL_OBJECT;

	if (LocalGetNumObjectId(l_typeAmort) == ARM_KO)
	{
		double endDateAmort = pEndDateAmort;
		double txAmort = pTxAmort;
		double amountAmort = pAmountAmort;

		_bstr_t intRuleAmort(pIntRuleAmort);
		CCString l_intRuleAmort =  intRuleAmort;

		_bstr_t freqAmort(pFreqAmort);
		CCString l_freqAmort = freqAmort;

		long freqAmortId; 
		if((freqAmortId = ARM_ConvFrequency(l_freqAmort, C_result)) == ARM_DEFAULT_ERR)
		{
			ERROR_MSG("Invalid Frq",pRet,ARM_ERROR_FREQ);
			return S_OK;
		}

		long intRuleAmortId= ARM_ConvIntRule(l_intRuleAmort);
		long dayCountAmortId = ARM_ConvDayCount(l_dayCountAmort);	

		// AmortMethodId
		long amortMethodId;
		if((amortMethodId = ARM_ConvAmortMethod (l_typeAmort, C_result)) == ARM_DEFAULT_ERR)
		{
			ERROR_MSG("Invalid TypeAmort",pRet,ARM_ERROR_FREQ);
			return S_OK;
		}

		long FixedLegAmortId;
		if (ARMLOCAL_FIXEDLEG(startDatePhase1,
							  endDateAmort,
							  K_RCV,
							  0,
							  1.0,
							  dayCountAmortId,
							  freqAmortId,
							  K_COMP_PROP,
							  K_ARREARS,
							  intRuleAmortId,
							  stubId,
							  false,
							  l_ccy,
							  "NULL",
							  K_NX_NONE,
							  -1.0,
							  K_NO,
						  -1,
							  C_result) == ARM_OK) 
							  
		{
			FixedLegAmortId = C_result.getLong();
		}
		else
		{
			ERROR_MSG(C_result.getMsg(),pRet);
			return S_OK;
		}

		if (ARMLOCAL_ARM_GenAmortization(FixedLegAmortId,
										 amortMethodId,
										 freqAmortId,
										 amountAmort,
										 -1,
										 Ntl,
										 txAmort,
										 0.,
										 ARM_NULL_OBJECT,
										 0.0,
										 C_result) == ARM_OK) 
		{
			amortRefValue = C_result.getLong();
			amortRefValueToDestroy = amortRefValue;
		}
		else
		{
			ERROR_MSG(C_result.getMsg(),pRet);
			return S_OK;
		}

		retCode = ARMLOCAL_FreeObject(FixedLegAmortId,C_result);
	}
	else
	{
		amortRefValue = LocalGetNumObjectId(l_typeAmort);

		// calcul du 1er nominal
		double dPaymentDate;

		if (ARMLOCAL_DisplayRefValue(amortRefValue,true,C_result) == ARM_OK)
		{
			dPaymentDate = C_result.getArray(0);
		}
		else
		{
			ERROR_MSG(C_result.getMsg(),pRet);
			return S_OK;
		}
		
		if (ARMLOCAL_CptRefValue(amortRefValue,dPaymentDate,C_result) == ARM_OK)
		{
			Ntl = C_result.getDouble();
		}
		else
		{
			ERROR_MSG(C_result.getMsg(),pRet);
			return S_OK;
		} 
	}
	//************************************************************************************************
	//*********************************** Fin Création Amort Refvalue ************************************
	//************************************************************************************************

	//************************************************************************************************
	//*********************************** Création Funding Leg Refvalue ******************************
	//************************************************************************************************
	long irIndexFundId;
	if (ARMLOCAL_IRINDEX(KACTUAL_360,
						 payFreqFundId,
						 -1,
						 K_COMP_PROP,
						 K_MOD_FOLLOWING,
						 K_ADVANCE,
						 -2,
						 K_ARREARS,
						 10000.,
						 false,
						 l_ccy,
						 indexFundId,
						 K_COMP_PROP,
						 adjFundId,
						 resetFreqFundId,
						 C_result) == ARM_OK)
	{
		irIndexFundId = C_result.getLong();
	}
	else
	{
		ERROR_MSG(C_result.getMsg(),pRet);
		return S_OK;
	}

	long swapLegId;
	if (ARMLOCAL_SWAPLEG(irIndexFundId,
						 startDatePhase1,
						 endDatePhase3,
						 K_RCV,
						 spreadType, 
						 spreadFund, 
						 false,
						 l_ccy,
						 dayCountFundId,
						 -2,
						 "NULL",
						 "NULL",
						 0,
						 K_NX_NONE,
						 stubId,
						 -1.0,
						 K_NO,
						  -1,
						 C_result) == ARM_OK) 
	{
		swapLegId = C_result.getLong();
	}
	else
	{
		ERROR_MSG(C_result.getMsg(),pRet);
		return S_OK;
	}
	
	// ********** Funding Leg ************
	long cloneId;
	if (ARMLOCAL_ClonedAndSetNotional(swapLegId,
									  amortRefValue,
									  100.,
									  C_result) == ARM_OK)
	{
		cloneId = C_result.getLong();
	}
	else
	{
		ERROR_MSG(C_result.getMsg(),pRet);
		return S_OK;
	}

	//*********   Price funding Leg   *************************
	double pxFund; 
	double pxFundDelta;

	if (ARMLOCAL_ARM_Price(cloneId,
						   LocalGetNumObjectId(l_smiledMod),
						   C_result) == ARM_OK)
	{
		pxFund = C_result.getDouble();
	}
	else
	{
		ERROR_MSG(C_result.getMsg(),pRet);
		return S_OK;
	}

	// shift seulement sur l'euro !?
	if (ARMLOCAL_ARM_Price(cloneId,
						   LocalGetNumObjectId(l_smiledModBump),
						   C_result) == ARM_OK)
	{
		pxFundDelta = C_result.getDouble();
	}
	else
	{
		ERROR_MSG(C_result.getMsg(),pRet);
		return S_OK;
	}

	retCode = ARMLOCAL_FreeObject(irIndexFundId, C_result);
	retCode = ARMLOCAL_FreeObject(swapLegId, C_result);
	retCode = ARMLOCAL_FreeObject(cloneId,C_result);
	
	// Fin du calcul du prix de la funding Leg

	//************************************************************************************************
	//*********************************** Fin Création Funding Leg Refvalue **************************
	//************************************************************************************************


	//************************************************************************************************
	//*********************************** Création Leg Phase1 Refvalue *******************************
	//************************************************************************************************
	long irIndexPhase1Id;
	double pxPhase1; 
	double pxPhase1Delta;
	long swapLegPhase1Id;
	long clonePhase1Id;
	double resetGapPhase1; 

	if (resetTimingPhase1Id == 1)
	{
		resetGapPhase1 = -2.; 
	}
	else
	{
		resetGapPhase1 = -15.; 
	}

	if ( startDatePhase1 != startDatePhase2 )
	{
		if (ARMLOCAL_IRINDEX(KACTUAL_360,
							 payFreqPhase1Id,
							 -1,
							 K_COMP_PROP,
							 K_MOD_FOLLOWING,
							 resetTimingPhase1Id,	
							 resetGapPhase1,
							 K_ARREARS,
							 10000.,
							 false,
							 l_ccy,
							 indexPhase1Id,
							 K_COMP_PROP,
							 adjPhase1Id,
							 resetFreqPhase1Id,
							 C_result) == ARM_OK) 
		{
			irIndexPhase1Id = C_result.getLong();
		}
		else
		{
			ERROR_MSG(C_result.getMsg(),pRet);
			return S_OK;
		}

		if (ARMLOCAL_SWAPLEG(irIndexPhase1Id,
							 startDatePhase1,
							 startDatePhase2,
							 K_RCV,
							 0, 
							 spreadPhase1, 
							 false,
							 l_ccy,
							 dayCountPhase1Id,
							 resetGapPhase1,
							 "NULL",
							 "NULL",
							 0,
							 K_NX_NONE,
							 stubId,
							 -1.0,
							 K_NO,
						  -1,
							 C_result) == ARM_OK) 
		{
			swapLegPhase1Id = C_result.getLong();
		}
		else
		{
			ERROR_MSG(C_result.getMsg(),pRet);
			return S_OK;
		}		
		
		// ********** Phase1 Leg ************

		
		if (ARMLOCAL_ClonedAndSetNotional(swapLegPhase1Id,
										  amortRefValue,
										  100.,
										  C_result) == ARM_OK)
		{
			clonePhase1Id = C_result.getLong();
		}
		else
		{
			ERROR_MSG(C_result.getMsg(),pRet);
			return S_OK;
		}

		//*********   Price Phase1 *************************
		if (ARMLOCAL_ARM_Price(clonePhase1Id,
							   LocalGetNumObjectId(l_smiledMod),
							   C_result) == ARM_OK)
		{
			pxPhase1 = C_result.getDouble();
		}
		else
		{
			ERROR_MSG(C_result.getMsg(),pRet);
			return S_OK;
		}

		// shift seulement sur l'euro !?
		if (ARMLOCAL_ARM_Price(clonePhase1Id,
							   LocalGetNumObjectId(l_smiledModBump),
							   C_result) == ARM_OK)
		{
			pxPhase1Delta = C_result.getDouble();
		}
		else
		{
			ERROR_MSG(C_result.getMsg(),pRet);
			return S_OK;
		}

		retCode = ARMLOCAL_FreeObject(irIndexPhase1Id,C_result);
		retCode = ARMLOCAL_FreeObject(swapLegPhase1Id,C_result);
		retCode = ARMLOCAL_FreeObject(clonePhase1Id,C_result);
	}
	else
	{
		pxPhase1 = 0.;
		pxPhase1Delta = 0.;
	}
	//************************************************************************************************
	//************************** Fin Création Leg Phase1 Refvalue ************************************
	//************************************************************************************************


	//************************************************************************************************
	//*********************************** Création Leg Phase2 Refvalue *******************************
	//************************************************************************************************

	long irIndexPhase2Id;
	double resetGapPhase2; 

	if (resetTimingPhase2Id == 1)
	{
		resetGapPhase2 = -2.; 
	}
	else
	{
		resetGapPhase2 = -15.; 
	}

	if (ARMLOCAL_IRINDEX(KACTUAL_360,
						 payFreqPhase2Id,
						 -1,
						 K_COMP_PROP,
						 K_MOD_FOLLOWING,
						 resetTimingPhase2Id,	
						 resetGapPhase2,
						 K_ARREARS,
						 10000.,
						 false,
						 l_ccy,
						 K_FIXED,
						 K_COMP_PROP,
						 adjPhase2Id,
						 resetFreqPhase2Id,
						 C_result) == ARM_OK) 
	{
		irIndexPhase2Id = C_result.getLong();
	}
	else
	{
		ERROR_MSG(C_result.getMsg(),pRet);
		return S_OK;
	}

	long swapLegPhase2Id;
	if (ARMLOCAL_SWAPLEG(irIndexPhase2Id,
						 startDatePhase2,
						 startDatePhase3,
						 K_RCV,
						 0,
						 spreadPhase2,
						 false,
						 l_ccy,
						 dayCountPhase2Id,
						 resetGapPhase2,
						 "NULL",
						 "NULL",
						 0, 
						 K_NX_NONE,
						 stubId, 
						 -1.0,
						 K_NO,
						  -1,
						 C_result) == ARM_OK) 
	{
		swapLegPhase2Id = C_result.getLong();
	}
	else
	{
		ERROR_MSG(C_result.getMsg(),pRet);
		return S_OK;
	}
	
	// ********** Phase2 Leg ************
	long clonePhase2Id;
	if (ARMLOCAL_ClonedAndSetNotional(swapLegPhase2Id,
									  amortRefValue,
									  100.,
									  C_result) == ARM_OK)
	{
		clonePhase2Id = C_result.getLong();
	}
	else
	{
		ERROR_MSG(C_result.getMsg(),pRet);
		return S_OK;
	}

	//*********   Price Phase2 *************************
	double pxPhase2;
	double pxPhase2Delta; 

	if (ARMLOCAL_ARM_Price(clonePhase2Id,
						   LocalGetNumObjectId(l_smiledMod),
						   C_result) == ARM_OK)
	{
		pxPhase2 = C_result.getDouble();
	}
	else
	{
		ERROR_MSG(C_result.getMsg(),pRet);
		return S_OK;
	}

	// shift seulement sur l'euro !?
	if (ARMLOCAL_ARM_Price(clonePhase2Id,
						   LocalGetNumObjectId(l_smiledModBump),
						   C_result) == ARM_OK)
	{
		pxPhase2Delta = C_result.getDouble();
	}
	else
	{
		ERROR_MSG(C_result.getMsg(),pRet);
		return S_OK;
	}

	retCode = ARMLOCAL_FreeObject(irIndexPhase2Id,C_result);
	retCode = ARMLOCAL_FreeObject(swapLegPhase2Id,C_result);
	retCode = ARMLOCAL_FreeObject(clonePhase2Id,C_result);

	//************************************************************************************************
	//************************** Fin Création Leg Phase2 Refvalue ************************************
	//************************************************************************************************



	//************************************************************************************************
	//*************************** Création Leg Phase3 Refvalue ***************************************
	//************************************************************************************************

	long irIndexPhase3Id;
	double pxPhase3; 
	double pxPhase3Delta;
	long swapLegPhase3Id;
	long clonePhase3Id;
	double resetGapPhase3; 

	if ( startDatePhase3 != endDatePhase3 )
	{		
		if (resetTimingPhase3Id == 1)
		{
			resetGapPhase3 = -2.; 
		}
		else
		{
			resetGapPhase3 = -15.; 
		}

		if (ARMLOCAL_IRINDEX(KACTUAL_360,
							 payFreqPhase3Id,
							 -1,
							 K_COMP_PROP,
							 K_MOD_FOLLOWING,
							 resetTimingPhase3Id,	
							 resetGapPhase3,
							 K_ARREARS,
							 10000.,
							 false,
							 l_ccy,
							 indexPhase3Id,
							 K_COMP_PROP,
							 adjPhase3Id,
							 resetFreqPhase3Id,
							 C_result) == ARM_OK) 
		{
			irIndexPhase3Id = C_result.getLong();
		}
		else
		{
			ERROR_MSG(C_result.getMsg(),pRet);
			return S_OK;
		}

		if (ARMLOCAL_SWAPLEG(irIndexPhase3Id,
							 startDatePhase3,
							 endDatePhase3,
							 K_RCV,
							 0, 
							 spreadPhase3, 
							 false,
							 l_ccy,
							 dayCountPhase3Id,
							 resetGapPhase2,
							 "NULL",
							 "NULL",
							 0,
							 K_NX_NONE,
							 stubId,
							 -1.0,
							 K_NO,
						  -1,
							 C_result) == ARM_OK) 
		{
			swapLegPhase3Id = C_result.getLong();
		}
		else
		{
			ERROR_MSG(C_result.getMsg(),pRet);
			return S_OK;
		}

		// ********** Phase3 Leg ************

		if (ARMLOCAL_ClonedAndSetNotional(swapLegPhase3Id,
										  amortRefValue,
										  100.,
										  C_result) == ARM_OK)
		{
			clonePhase3Id = C_result.getLong();
		}
		else
		{
			ERROR_MSG(C_result.getMsg(),pRet);
			return S_OK;
		}

		//*********   Price Phase3 *************************	
		if (ARMLOCAL_ARM_Price(clonePhase3Id,
							   LocalGetNumObjectId(l_smiledMod),
							   C_result) == ARM_OK)
		{
			pxPhase3 = C_result.getDouble();
		}
		else
		{
			ERROR_MSG(C_result.getMsg(),pRet);
			return S_OK;
		}

		// shift seulement sur l'euro !?
		if (ARMLOCAL_ARM_Price(clonePhase3Id,
							   LocalGetNumObjectId(l_smiledModBump),
							   C_result) == ARM_OK)
		{
			pxPhase3Delta = C_result.getDouble();
		}
		else
		{
			ERROR_MSG(C_result.getMsg(),pRet);
			return S_OK;
		}
		retCode = ARMLOCAL_FreeObject(irIndexPhase3Id,C_result);
		retCode = ARMLOCAL_FreeObject(swapLegPhase3Id,C_result);
		retCode = ARMLOCAL_FreeObject(clonePhase3Id,C_result);
	}
	else
	{
		pxPhase3 = 0.;
		pxPhase3Delta = 0.; 
	}

	double pxPhase2DIG; 
	double pxTotal; 
	double pxDelta; 
	double SensiTotal; 
	
	double prixInter;
	double Xincr;
	
	double resetGapPhase2DIG; 

	if (resetTimingPhase2DIGId == 1)
	{
		resetGapPhase2DIG = -2.; 
	}
	else
	{
		resetGapPhase2DIG = -15.; 
	}

	//Construction de l'index PayOff
	long dayCountPay;
	long irIndexPayId; 
	long indexPayId; 
	indexPayId = ARM_ConvIrType(l_indexPay); 

	if (IsCMSIndex((ARM_INDEX_TYPE)indexPayId) == 1)
	{
		dayCountPay = K30_360 ; 
	}	
	else
	{
		dayCountPay = KACTUAL_360; 
	}

	if (ARMLOCAL_IRINDEX(dayCountPay,
						 payFreqPhase2Id,
						 -1,
						 K_COMP_PROP,
						 K_MOD_FOLLOWING,
				 		 resetTimingPhase2DIGId,	
						 resetGapPhase2DIG,
						 K_ARREARS,
						 10000.,
						 false,
						 l_ccy,
						 indexPayId,
						 K_COMP_PROP,
						 adjPhase2DIGId,
						 resetFreqPhase2Id,
						 C_result) == ARM_OK) 
	{
		irIndexPayId = C_result.getLong();
	}
	else
	{
		ERROR_MSG(C_result.getMsg(),pRet);
		return S_OK;
	}

//******************** Calcul premiere PV ******************
	int iteration =1;
	double prix;
	double prixAct; 
	double prixAct2;
	double SprSeek = 10.01;
	double SprSeekFinal= 10.01;
	double SprSeekBonifie = 10.01;  // PUREMENT ARBITRAIRE
	fee = Ntl * fee / 100;
//******************** Calcul premiere PV ******************

		VECTOR<double> pxPentibor = computePentibor(l_indexPhase2DIG,l_indexLongPhase2DIG, irIndexPayId, startDatePhase2, startDatePhase3, 
			l_floorOrCap, l_strikePhase2DIG, SprSeekBonifie, ccyId, 
			dayCountPhase2Id, 
			resetTimingPhase2DIGId, payTimingPhase2DIGId,  resetGapPhase2DIG, resetFreqPhase2Id, payFreqPhase2Id, 
			adjPhase2DIGId, stubId,amortRefValue, 
			l_hyperCubeCorrel, l_volCurvFromMatriceShift,
			l_vol, l_volCub, 
			l_convexityManager, l_zc, l_smiledMod, l_indexIndexCorrelCube, l_corrEUR,l_interCorr);

	if (pxPentibor.size() == 0)
	{
		ERROR_MSG("Pb in computePentibor",pRet);
		return S_OK;
	}

	pxPhase2DIG = pxPentibor[0];
	double pxPhase2DIGBumpCorrel = pxPentibor[1];
	double MT = pxPentibor[2];

	pxTotal = pxPhase1 + pxPhase2 + pxPhase3 - pxFund + pxPhase2DIG;
	pxDelta = fabs(pxPhase2DIGBumpCorrel-pxPhase2DIG) + MT; 
	if (pxDelta <0.) pxDelta = - pxDelta;

	SensiTotal = pxDelta;
	prixAct = pxTotal - SensiTotal - fee; 

	prix = prixAct ;
	prixAct2 =1;
	
//********************                    ******************
//******************** BOUCLE SUR LE PRICE ******************
	while ((abs(prixAct2) > 0.0001) && (iteration < 10 ))
	{
		SprSeek = SprSeekBonifie;

		//*************************
		// **** Premier recalcul ****
		//****************************

		pxPentibor = computePentibor(l_indexPhase2DIG,l_indexLongPhase2DIG, irIndexPayId, startDatePhase2, startDatePhase3, 
			l_floorOrCap, l_strikePhase2DIG, SprSeekBonifie , ccyId, 
			dayCountPhase2Id, 
			resetTimingPhase2DIGId, payTimingPhase2DIGId,  resetGapPhase2DIG, resetFreqPhase2Id, payFreqPhase2Id, 
			adjPhase2DIGId, stubId, amortRefValue,
			l_hyperCubeCorrel, l_volCurvFromMatriceShift,
			l_vol, l_volCub,
			l_convexityManager, l_zc, l_smiledMod, l_indexIndexCorrelCube, l_corrEUR,l_interCorr);

		if (pxPentibor.size() == 0)
		{
			ERROR_MSG("Pb in computePentibor",pRet);
			return S_OK;
		}

		pxPhase2DIG = pxPentibor[0];
		pxPhase2DIGBumpCorrel = pxPentibor[1];
		MT = pxPentibor[2];

		pxTotal = pxPhase1 + pxPhase2 + pxPhase3 - pxFund + pxPhase2DIG;
		pxDelta = fabs(pxPhase2DIGBumpCorrel-pxPhase2DIG) + MT; 
		if (pxDelta <0.) pxDelta = - pxDelta;

		SensiTotal = pxDelta;
		prixAct = pxTotal - SensiTotal - fee; 

		// les deltas se compensent et on somme en valeur absolue les deltas et vegas

		//******************************
		// **** affectation valeurs ****
		//******************************

		prixAct2 = prixAct ;
		Xincr = SprSeek / 1000000;
		SprSeekBonifie = SprSeek + Xincr;		

		//****************************
		// **** Deuxieme recalcul ****
		//****************************
		pxPentibor = computePentibor(l_indexPhase2DIG,l_indexLongPhase2DIG, irIndexPayId, startDatePhase2, startDatePhase3, 
									l_floorOrCap, l_strikePhase2DIG, SprSeekBonifie , ccyId, 
									dayCountPhase2Id, 
			resetTimingPhase2DIGId, payTimingPhase2DIGId,  resetGapPhase2DIG, resetFreqPhase2Id, payFreqPhase2Id, 
			adjPhase2DIGId, stubId, amortRefValue, 
			l_hyperCubeCorrel, l_volCurvFromMatriceShift, 
			l_vol, l_volCub, 
			l_convexityManager, l_zc, l_smiledMod, l_indexIndexCorrelCube, l_corrEUR,l_interCorr);

		if (pxPentibor.size() == 0)
		{
			ERROR_MSG("Pb in computePentibor",pRet);
			return S_OK;
		}

		pxPhase2DIG = pxPentibor[0];
		pxPhase2DIGBumpCorrel = pxPentibor[1];
		MT = pxPentibor[2];

		pxTotal = pxPhase1 + pxPhase2 + pxPhase3 - pxFund + pxPhase2DIG;
		pxDelta = fabs(pxPhase2DIGBumpCorrel-pxPhase2DIG) + MT; 
		if (pxDelta <0.) pxDelta = - pxDelta;

		SensiTotal = pxDelta;
		prixAct = pxTotal - SensiTotal - fee;

		iteration++;

		SprSeekFinal = SprSeekBonifie;
		prixInter = ( prixAct - prixAct2 ) / Xincr;
		// if prixinter =  0 sortir erreur
		SprSeekBonifie = SprSeek - prixAct / prixInter ; 
		prixAct2 = prixAct;		
	}

//**************** PRICE *******************
	retCode = ARMLOCAL_FreeObject(irIndexPayId,C_result);
	retCode = ARMLOCAL_FreeObject(amortRefValueToDestroy,C_result);
	retCode = ARMLOCAL_FreeObject(ccyId,C_result);

	VECTOR<double> vRes;

	vRes.push_back(spreadPhase2 + SprSeekFinal);
	vRes.push_back(spreadPhase2 + SprSeekFinal);

	Res = VECTORDOUBLE2VARIANT(vRes, pRet);

	_variant_t wrap_Res;
	wrap_Res.Attach(*pRet);
	*pRet=wrap_Res.Detach();

	return S_OK;
}
catch(...)
{
	return E_FAIL;
}
}


VECTOR<double> computePentibor(const CCString& index, const CCString& indexLong, long indexPayId, double startDate, double endDate,
							   const CCString& floorOrCap, const CCString& strike, double payOffSpread, long ccyId, long dayCountId,
							   long resetTimingId, long payTimingId, double resetGap, long resetFreqId, long payFreqId,
							   long intRuleId, long stubRuleId, long amortRefvalue,
							   const CCString& hyperCubeCorrel, const CCString& VolCurvFromMatriceShift,
							   const CCString& vol, const CCString& volCUB, 
							   const CCString& ConvexityManager, const CCString& zc, const CCString& smiledMod, const CCString& indexIndexCorrelCube, 
							const CCString& corrEUR, const CCString& InterCorr)
{
	VECTOR<double> resultat;
	resultat.clear();

try
{
	double pxDigital;
	double pxBumpDigital; 

	long digitalId;
	long cloneDigitalId; 
	long corrManagerId; 

	double shiftCorrel; 
	
	long bumpHyperCubeCorrelId;
	long bumpCorrManagerId; 

	ARM_result C_result;

	// Modèles
	long bsModGenId; 
	long BumpBsModGenId;

	VECTOR <double> spreadVector;
	VECTOR <double> weightVector;

   	VECTOR<double> C_fixing1_double_vect;
   	VECTOR<double> C_fixing2_double_vect;
	VECTOR<double> C_fixing3_double_vect;

	long slopeFlag;
	int cptStrikeMethod;
	int computedFormula;

	slopeFlag = 1;
	cptStrikeMethod = 1.;
	computedFormula = 1.;

	// Ss-Jacent
	long ssJacentId;
	long cloneSsJacentId;
	double pxSsJacent;
	
	double d_index; 
	long indexId; 
	indexId = ARM_ConvIrType(index); 

	if (IsCMSIndex((ARM_INDEX_TYPE)indexId) == 1)
	{
		d_index = (indexId - K_CMS1 +1 ); 
	}	
	else
	{
		d_index = 1./FromLiborTypeToFrequency(indexId); 
	}
	
	double d_indexLong; 
	long indexLongId; 
	indexLongId = ARM_ConvIrType(indexLong); 

	if (IsCMSIndex((ARM_INDEX_TYPE)indexLongId) ==1)
	{
		d_indexLong = indexLongId - K_CMS1 + 1; 
	}	
	else
	{
		d_indexLong = 1./FromLiborTypeToFrequency(indexLongId); 
	}

	//Définition du sreadVector et du WeightVector
	spreadVector.push_back(-0.01);
	spreadVector.push_back(0.01); 

	weightVector.push_back(1);
	weightVector.push_back(1);

	long capOrFloorId; 

	if((capOrFloorId = ARM_ConvCapOrFloor(floorOrCap, C_result)) == ARM_DEFAULT_ERR)
	{
		return resultat;
	}

	// Les variables indexIndexCorrelCube, corrEUR, InterCorr sont à créer
	// La fonction ARMLOCAL_CreateGenCorrelatorManager est en attente de modifcation par JPR
 
	VECTOR<CCString> vCcy; 
	VECTOR<long> vHyperCubeCorrelId; 
	VECTOR<long> vIndexIndexCorrelCubeId; 
	VECTOR<long> vCorrEURId; 
	VECTOR<long> vInterCorrId; 
	VECTOR<long> vIrVolId; 
	VECTOR<long> vVolVolId; 
	VECTOR<long> vFXVolHyperCubeId; 

	vCcy.push_back("EUR"); 
	vHyperCubeCorrelId.push_back(LocalGetNumObjectId(hyperCubeCorrel));
	vIndexIndexCorrelCubeId.push_back(LocalGetNumObjectId(indexIndexCorrelCube));
	vInterCorrId.push_back(LocalGetNumObjectId(InterCorr));

	if (ARMLOCAL_CreateGenCorrelatorManager(vCcy, vHyperCubeCorrelId, 
			vIndexIndexCorrelCubeId, vCorrEURId, vInterCorrId, vIrVolId, vVolVolId, vFXVolHyperCubeId, C_result) ==  ARM_OK)
	{
		corrManagerId = C_result.getLong();
	}
	else
	{
		return resultat;
	}
	
	long modelTypeId;
	long spreadVolTypeId;
	CCString modelType; 
	CCString spreadVolType;
	bool isLnVol=true;
	long payMargin_type; 
	double payMarginId;
	long payIdxWeight; 
	double freezeFixing = 0.0;

	payMargin_type = 0; 
	payMarginId = payOffSpread;
	payIdxWeight = 1; 

	modelType = "2LOG";
	spreadVolType = "C"; 

	if ((modelTypeId = ARM_ConvModelType (modelType, C_result)) == ARM_DEFAULT_ERR)
	{
		return resultat;
	}

	if ((spreadVolTypeId = ARM_ConvVolType2 (spreadVolType, C_result)) == ARM_DEFAULT_ERR)
	{
		return resultat;
	}

	VECTOR<double> C_rate; 
	VECTOR<CCString> C_matu;
	double C_meth = 0;

	//************************************************************************************************
	//*************************** Création Corridor **************************************************
	//************************************************************************************************

	//Calcul du prix de la Digitale

	if (ARMLOCAL_SPREADCORRIDOR(startDate,
						 endDate,
						 capOrFloorId,
						 1, 
						 LocalGetNumObjectId(strike), 
						 indexPayId, 
						 indexId, 
						 indexLongId,
						 weightVector[0], 
						 weightVector[1],
						 dayCountId,
						 resetFreqId, 
						 payFreqId, 
						 resetTimingId,
						 payTimingId,				
						 ccyId,
						 resetGap, 
						 spreadVector[0], 
						 spreadVector[1],
						 1,
						 ARM_NULL_OBJECT,
						 C_fixing1_double_vect,
						 1,
						 ARM_NULL_OBJECT,
						 C_fixing2_double_vect,
						 intRuleId, 
						 stubRuleId, 
						 slopeFlag,
						 cptStrikeMethod, 
						 payMargin_type, 
						 payMarginId,  
						 payIdxWeight, 
						 1,
						 ARM_NULL_OBJECT,
						 C_fixing3_double_vect,
						 freezeFixing,
						 computedFormula, 					
						 C_result) == ARM_OK)
	{
		digitalId = C_result.getLong();
	}
	else
	{
		return resultat;
	}

	if (ARMLOCAL_ClonedAndSetNotional(digitalId,
									  amortRefvalue,
									  100.,
									  C_result) == ARM_OK)
	{
		cloneDigitalId = C_result.getLong();
	}
	else
	{
		return resultat;
	}

	long retCode = ARMLOCAL_FreeObject(digitalId,C_result);

	if (ARMLOCAL_BSMODELGEN(LocalGetNumObjectId(zc), ARM_NULL_OBJECT, LocalGetNumObjectId(vol), LocalGetNumObjectId(volCUB), 
				corrManagerId, LocalGetNumObjectId(ConvexityManager), 
				LocalGetNumObjectId(zc), LocalGetNumObjectId(hyperCubeCorrel), ARM_NULL_OBJECT, ARM_NULL_OBJECT,
				modelTypeId, spreadVolTypeId, LocalGetNumObjectId(smiledMod), isLnVol, 100, C_result) == ARM_OK)
	{
		bsModGenId = C_result.getLong();
	}
	else
	{
		return resultat;
	}

	if (ARMLOCAL_ARM_Price(cloneDigitalId,
						   bsModGenId,
						   C_result) == ARM_OK)
	{
		pxDigital = C_result.getDouble();
	}
	else
	{
		return resultat;
	}

	//  Calcul de la marge Correl
	if (ARMLOCAL_ARM_BumpVolatility(LocalGetNumObjectId(hyperCubeCorrel),1., 0, 0, K_NO, K_YES,C_result) == ARM_OK)
	{
		bumpHyperCubeCorrelId = C_result.getLong();
	}
	else
	{
		return resultat;
	}
	
	// Les variables indexIndexCorrelCube, corrEUR, InterCorr sont à créer
	// La fonction ARMLOCAL_CreateGenCorrelatorManager est en attente de modifcation par JPR
	vHyperCubeCorrelId.clear();
	vHyperCubeCorrelId.push_back(bumpHyperCubeCorrelId); 
	
	if (ARMLOCAL_CreateGenCorrelatorManager(vCcy, vHyperCubeCorrelId, 
			vIndexIndexCorrelCubeId, vCorrEURId, vInterCorrId, vIrVolId, vVolVolId, vFXVolHyperCubeId, C_result) ==  ARM_OK)
    {
		bumpCorrManagerId = C_result.getLong();
	}
	else
	{
		return resultat;
	}

	if (ARMLOCAL_BSMODELGEN(LocalGetNumObjectId(zc), ARM_NULL_OBJECT, LocalGetNumObjectId(vol), LocalGetNumObjectId(volCUB), 
				bumpCorrManagerId, LocalGetNumObjectId(ConvexityManager), 
				LocalGetNumObjectId(zc), bumpHyperCubeCorrelId, ARM_NULL_OBJECT, ARM_NULL_OBJECT,
                modelTypeId, spreadVolTypeId, LocalGetNumObjectId(smiledMod), isLnVol, 100, C_result) == ARM_OK)
	{
		BumpBsModGenId = C_result.getLong();
	}
	else
	{
		return resultat;
	}

	if (ARMLOCAL_ARM_Price(cloneDigitalId,
						   BumpBsModGenId,
						   C_result) == ARM_OK)
	{
		pxBumpDigital = C_result.getDouble();
	}
	else
	{
		return resultat;
	}

	retCode = ARMLOCAL_FreeObject(bumpHyperCubeCorrelId,C_result);
	retCode = ARMLOCAL_FreeObject(bumpCorrManagerId,C_result);
	retCode = ARMLOCAL_FreeObject(BumpBsModGenId,C_result);

	if ((pxBumpDigital - pxDigital) > 0 )
	{				
		if (ARMLOCAL_ComputeVolatility(LocalGetNumObjectId(VolCurvFromMatriceShift), d_index, d_indexLong, 0, C_result) == ARM_OK)
		{
			shiftCorrel = C_result.getDouble();
		}
		else
		{
			return resultat;
		}
	}
	else
	{
		if (ARMLOCAL_ComputeVolatility(LocalGetNumObjectId(VolCurvFromMatriceShift),d_indexLong, d_index, 0, C_result) == ARM_OK)
		{
			shiftCorrel = C_result.getDouble();
		}
		else
		{
			return resultat;
		}
	}

	// Bump de l'hyperCube 
	if (ARMLOCAL_ARM_BumpVolatility(LocalGetNumObjectId(hyperCubeCorrel),shiftCorrel, 0, 0, K_NO, K_YES,C_result) == ARM_OK)
	{
		bumpHyperCubeCorrelId = C_result.getLong();
	}
	else
	{
		return resultat;
	}

	// Les variables indexIndexCorrelCube, corrEUR, InterCorr sont à créer
	// La fonction ARMLOCAL_CreateGenCorrelatorManager est en attente de modifcation par JPR
	vHyperCubeCorrelId.clear();
	vHyperCubeCorrelId.push_back(bumpHyperCubeCorrelId); 

	if (ARMLOCAL_CreateGenCorrelatorManager(vCcy, vHyperCubeCorrelId, 
			vIndexIndexCorrelCubeId, vCorrEURId, vInterCorrId, vIrVolId, vVolVolId, vFXVolHyperCubeId, C_result) ==  ARM_OK)
	{
		bumpCorrManagerId = C_result.getLong();
	}
	else
	{
		return resultat;
	}

	if (ARMLOCAL_BSMODELGEN(LocalGetNumObjectId(zc), ARM_NULL_OBJECT, LocalGetNumObjectId(vol), LocalGetNumObjectId(volCUB), 
				bumpCorrManagerId, LocalGetNumObjectId(ConvexityManager),
				LocalGetNumObjectId(zc), bumpHyperCubeCorrelId, ARM_NULL_OBJECT, ARM_NULL_OBJECT,
				modelTypeId, spreadVolTypeId, LocalGetNumObjectId(smiledMod), isLnVol, 100, C_result) == ARM_OK)
	{
		BumpBsModGenId = C_result.getLong();
	}
	else
	{
		return resultat;
	}

	if (ARMLOCAL_ARM_Price(cloneDigitalId,
						   BumpBsModGenId,

						   C_result) == ARM_OK)
	{
		pxBumpDigital = C_result.getDouble();
	}
	else
	{
		return resultat;
	}

	retCode = ARMLOCAL_FreeObject(bumpCorrManagerId,C_result);
	retCode = ARMLOCAL_FreeObject(bumpHyperCubeCorrelId,C_result);
	retCode = ARMLOCAL_FreeObject(BumpBsModGenId,C_result);
	retCode = ARMLOCAL_FreeObject(cloneDigitalId,C_result);
	
	//fin calcul marge correl


	//calcul de la discontinuité
	vector<double> calibInfos;

	if (ARMLOCAL_SPREADDIGITAL(startDate,
							   endDate,
							   capOrFloorId,
							   1, 
							   LocalGetNumObjectId(strike),								
							   0,
							   payOffSpread,
							   indexId, 
							   indexLongId,
							   weightVector[0], 
							   weightVector[1],
							   dayCountId,
							   resetFreqId,
							   payFreqId, 
							   resetTimingId,
							   payTimingId, 
							   ccyId,
							   resetGap, 
							   spreadVector[0], 
							   spreadVector[1],
							   1,
							   ARM_NULL_OBJECT,
							   C_fixing1_double_vect,
							   1,
							   ARM_NULL_OBJECT,
							   C_fixing2_double_vect,
							   intRuleId,
							   stubRuleId,
							   slopeFlag,
							   cptStrikeMethod, 
							   computedFormula,
							   calibInfos,
							   C_result) == ARM_OK)
	{
		ssJacentId = C_result.getLong();
	}
	else
	{
		return resultat;
	}

	if (ARMLOCAL_ClonedAndSetNotional(ssJacentId,
									  amortRefvalue,
									  100.,
									  C_result) == ARM_OK)
	{
		cloneSsJacentId = C_result.getLong();
	}
	else
	{
		return resultat;
	}

	retCode = ARMLOCAL_FreeObject(ssJacentId,C_result);


	if (ARMLOCAL_ARM_Price(cloneSsJacentId,
						   bsModGenId,
						   C_result) == ARM_OK)
	{
		pxSsJacent = C_result.getDouble();
	}
	else
	{
		return resultat;
	}

	retCode = ARMLOCAL_FreeObject(cloneSsJacentId,C_result);
	
	long irIndexIdfixed;
	if (ARMLOCAL_IRINDEX(dayCountId,
						 resetFreqId,
						 -1.0,
						 K_COMP_PROP,
						 K_MOD_FOLLOWING,
						 resetTimingId,
						 resetGap,
						 K_ARREARS,
						 10000.,
						 false,
						 "EUR",
						 K_FIXED,
						 K_COMP_PROP,
						 intRuleId,
						 resetFreqId,
						 C_result) == ARM_OK) 
	{
		irIndexIdfixed = C_result.getLong();
	}
	else
	{
		return resultat;
	}

	long swapLegFixedId;
	if (ARMLOCAL_SWAPLEG(irIndexIdfixed,
						 startDate,
						 endDate,
						 K_RCV,
						 0,
						 payOffSpread,
						 false,
						 "EUR",
						 dayCountId,
						 resetGap,
						 "NULL",
						 "NULL",
						 0, 
						 K_NX_NONE,
						 stubRuleId, 
						 -1.0,
						 K_NO,
						  -1,
						 C_result) == ARM_OK) 
	{
		swapLegFixedId = C_result.getLong();
	}
	else
	{
		return resultat;
	}	

	// ********** Fixed Leg ************
	long cloneFixedId;
	if (ARMLOCAL_ClonedAndSetNotional(swapLegFixedId,
									  amortRefvalue,
									  100.,
									  C_result) == ARM_OK)
	{
		cloneFixedId = C_result.getLong();
	}
	else
	{
		return resultat;
	}

	retCode = ARMLOCAL_FreeObject(swapLegFixedId,C_result);
	
	//*********   Price Fixed *************************
	double pxFixed;

	if (ARMLOCAL_ARM_Price(cloneFixedId,
						   LocalGetNumObjectId(smiledMod),
						   C_result) == ARM_OK)
	{
		pxFixed = C_result.getDouble();
	}
	else
	{
		return resultat;
	}

	retCode = ARMLOCAL_FreeObject(cloneFixedId,C_result);

	double ProbaExerciceFixed; 
	double MTFixed; 
	ProbaExerciceFixed = pxSsJacent/pxFixed; 

	//Calcul de la marge technique
	double minMTFixed;
	double maxMTFixed;
	double maxPriceFixed; 

	maxPriceFixed = pxFixed;

	if (ProbaExerciceFixed > (1- ProbaExerciceFixed))
	{
		minMTFixed = 1. - ProbaExerciceFixed;
	}
	else
	{
		minMTFixed = ProbaExerciceFixed;
	}
	
	if (minMTFixed - 0.1 > 0)
	{
		maxMTFixed = minMTFixed - 0.1; 
	}
	else
	{
		maxMTFixed = 0.;
	}

	MTFixed = (0.015 + 0.075 *maxMTFixed)*maxPriceFixed;

// Calcul MT variable

	long ssJacentCorridorId; 
	long cloneSsJacentCorridorId;
	double pxSsJacentCorridor;

	if (ARMLOCAL_SPREADCORRIDOR(startDate,
						 endDate,
						 capOrFloorId,	
						 1,
						 LocalGetNumObjectId(strike),	
						 indexPayId, 
						 indexId, 
						 indexLongId,
						 weightVector[0], 
						 weightVector[1],
						 dayCountId,
						 resetFreqId, 
						 payFreqId, 
						 resetTimingId, 
						 payTimingId, 					
						 ccyId,
						 resetGap, 
						 spreadVector[0], 
						 spreadVector[1],
						 1,
						 ARM_NULL_OBJECT,
						 C_fixing1_double_vect,
						 1,
						 ARM_NULL_OBJECT,
						 C_fixing2_double_vect,
						 intRuleId, 
						 stubRuleId, 
						 slopeFlag,
						 cptStrikeMethod, 
						 payMargin_type, 
						 0.,  
						 payIdxWeight, 
						 1,
						 ARM_NULL_OBJECT,
						 C_fixing3_double_vect,
						 freezeFixing,
						 computedFormula,					
						 C_result) == ARM_OK)
	{
		ssJacentCorridorId = C_result.getLong();
	}
	else
	{
		return resultat;
	}
	
	if (ARMLOCAL_ClonedAndSetNotional(ssJacentCorridorId,
									  amortRefvalue,
									  100.,
									  C_result) == ARM_OK)
	{
		cloneSsJacentCorridorId = C_result.getLong();
	}
	else
	{
		return resultat;
	}

	retCode = ARMLOCAL_FreeObject(ssJacentCorridorId,C_result);

	if (ARMLOCAL_ARM_Price(cloneSsJacentCorridorId,
						   bsModGenId,
						   C_result) == ARM_OK)
	{
		pxSsJacentCorridor = C_result.getDouble();
	}
	else
	{
		return resultat;
	}

	retCode = ARMLOCAL_FreeObject(corrManagerId,C_result);
	retCode = ARMLOCAL_FreeObject(bsModGenId,C_result);	
	retCode = ARMLOCAL_FreeObject(cloneSsJacentCorridorId,C_result);

	long swapLegVariableId;
	if (ARMLOCAL_SWAPLEG(indexPayId,
						 startDate,
						 endDate,
						 K_RCV,
						 0,
						 0.,
						 false,
						 "EUR",
						 dayCountId, 
						 resetGap,
						 "NULL",
						 "NULL",
						 0,
						 K_NX_NONE,
						 stubRuleId, 
						 -1.0,
						 K_NO,
								 -1,
						 C_result) == ARM_OK) 
	{
		swapLegVariableId = C_result.getLong();
	}
	else
	{
		return resultat;
	}	

	// ********** Variable Leg ************
	long cloneVariableId;
	if (ARMLOCAL_ClonedAndSetNotional(swapLegVariableId,
									  amortRefvalue,
									  100.,
									  C_result) == ARM_OK)
	{
		cloneVariableId = C_result.getLong();
	}
	else
	{
		return resultat;
	}

	retCode = ARMLOCAL_FreeObject(swapLegVariableId,C_result);
	
	//*********   Price Variable *************************
	double pxVariable;

	if (ARMLOCAL_ARM_Price(cloneVariableId,
						   LocalGetNumObjectId(smiledMod),
						   C_result) == ARM_OK)
	{
		pxVariable = C_result.getDouble();
	}
	else
	{
		return resultat;
	}

	retCode = ARMLOCAL_FreeObject(cloneVariableId,C_result);

	double ProbaExerciceVariable; 
	ProbaExerciceVariable = pxSsJacentCorridor/pxVariable;
	
	//Fin de la récup du prix du ssJacent
	
	//Calcul de la marge technique
	double minMTVariable;
	double maxMTVariable;
	double maxPriceVariable; 
	double MTVariable; 

	maxPriceVariable = pxVariable;

	if (ProbaExerciceVariable > (1- ProbaExerciceVariable))
	{
		minMTVariable = 1. - ProbaExerciceVariable;
	}
	else
	{
		minMTVariable = ProbaExerciceVariable;
	}

	if (minMTVariable - 0.1 > 0)
	{
		maxMTVariable = minMTVariable - 0.1; 
	}
	else
	{
		maxMTVariable = 0.;
	}

	MTVariable = (0.015 + 0.075 *maxMTVariable)*maxPriceVariable;

	double MTGlobal; 
	MTGlobal = fabs(MTVariable + MTFixed); 

	//Calcul de la marge Minimale
	double SpreadMin; 
	SpreadMin = 0.8; 

	long swapLegMinId;
	if (ARMLOCAL_SWAPLEG(irIndexIdfixed,
						 startDate,
						 endDate,
						 K_RCV,
						 0,
						 SpreadMin,
						 false,
						 "EUR",
						 dayCountId,
						 resetGap,
						 "NULL",
						 "NULL",
						 0, 
						 K_NX_NONE,
						 stubRuleId, 
						 -1.0,
						 K_NO,
								 -1,
						 C_result) == ARM_OK) 
	{
		swapLegMinId = C_result.getLong();
	}
	else
	{
		return resultat;
	}

	retCode = ARMLOCAL_FreeObject(irIndexIdfixed,C_result);

	// ********** Fixed Leg ************
	long cloneMinId;
	if (ARMLOCAL_ClonedAndSetNotional(swapLegMinId,
									  amortRefvalue,
									  100.,
									  C_result) == ARM_OK)
	{
		cloneMinId = C_result.getLong();
	}
	else
	{
		return resultat;
	}

	retCode = ARMLOCAL_FreeObject(swapLegMinId,C_result);
	
	//*********   Price Fixed *************************	
	double pxMin;

	if (ARMLOCAL_ARM_Price(cloneMinId,
						   LocalGetNumObjectId(smiledMod),
						   C_result) == ARM_OK)
	{
		pxMin = C_result.getDouble();
	}
	else
	{
		return resultat;
	}
	
	retCode = ARMLOCAL_FreeObject(cloneMinId,C_result);	
	
	double ProbaExerciceMin; 
	double minMTMin; 
	double maxMTMin; 
	double maxPriceMin; 
	double MTMin; 

	maxPriceMin = fabs(pxMin);

	ProbaExerciceMin = (ProbaExerciceVariable + ProbaExerciceFixed)/2.; 

	if (ProbaExerciceMin > (1- ProbaExerciceMin))
	{
		minMTMin = 1. - ProbaExerciceMin;
	}
	else
	{
		minMTMin = ProbaExerciceMin;
	}
	
	if (minMTMin - 0.1 > 0)
	{
		maxMTMin = minMTMin - 0.1; 
	}
	else
	{
		maxMTMin = 0.;
	}

	MTMin = (0.015 + 0.075 *maxMTMin)*maxPriceMin;

	MTMin = fabs(MTMin); 

	double MTFinal; 

	if (MTMin - MTGlobal> 0)
	{
		MTFinal = MTMin; 
	}
	else
	{
		MTFinal = MTGlobal;
	}
		
	//************************************************************************************************
	//*************************** Fin création Digitale **********************************************
	//************************************************************************************************



	//*********   Price Première Phase structure 1  TD*************************

	VECTOR<double> vRes;
	
	vRes.push_back(pxDigital);
	vRes.push_back(pxBumpDigital);
	vRes.push_back(MTFinal);

	vHyperCubeCorrelId.clear(); 
	vIndexIndexCorrelCubeId.clear(); 
	vCorrEURId.clear(); 
	vInterCorrId.clear(); 

	retCode = ARMLOCAL_FreeObject(swapLegFixedId,C_result);
	retCode = ARMLOCAL_FreeObject(cloneFixedId,C_result);
 	retCode = ARMLOCAL_FreeObject(swapLegVariableId,C_result);
 	retCode = ARMLOCAL_FreeObject(cloneVariableId,C_result);
 	retCode = ARMLOCAL_FreeObject(swapLegMinId,C_result);
 	retCode = ARMLOCAL_FreeObject(cloneMinId,C_result);
 	retCode = ARMLOCAL_FreeObject(digitalId,C_result);
 	retCode = ARMLOCAL_FreeObject(cloneDigitalId,C_result);
 	retCode = ARMLOCAL_FreeObject(ssJacentCorridorId,C_result);
 	retCode = ARMLOCAL_FreeObject(cloneSsJacentCorridorId,C_result);
 	retCode = ARMLOCAL_FreeObject(corrManagerId,C_result);
 	retCode = ARMLOCAL_FreeObject(bumpCorrManagerId,C_result);
 	retCode = ARMLOCAL_FreeObject(bsModGenId,C_result);
 	retCode = ARMLOCAL_FreeObject(BumpBsModGenId,C_result);
 	retCode = ARMLOCAL_FreeObject(bumpHyperCubeCorrelId,C_result);
 	retCode = ARMLOCAL_FreeObject(ssJacentId,C_result);
 	retCode = ARMLOCAL_FreeObject(cloneSsJacentId,C_result);
 	retCode = ARMLOCAL_FreeObject(irIndexIdfixed,C_result);

	return vRes;
}
catch(...)
{
	return resultat;
}
}


VECTOR<double> computePentilix(const CCString& index, const CCString& indexLong, long indexPayId, double startDate, double endDate,
							   const CCString& floorOrCap, const CCString& strike, double payOffSpread, long ccyId, long dayCountId,
							   long resetTimingId, long payTimingId, double resetGap, long resetFreqId, long payFreqId,
							   long intRuleId, long stubRuleId, long amortRefvalue,
							   const CCString& hyperCubeCorrel, const CCString& VolCurvFromMatriceShift,
							   const CCString& vol, const CCString& volCUB, 
							   const CCString& ConvexityManager, const CCString& zc, const CCString& smiledMod, const CCString& indexIndexCorrelCube, 
							const CCString& corrEUR, const CCString& InterCorr);

STDMETHODIMP ActiveXModule::ARMcomputePentilix(double pNtl, 
										   double pStartDatePhase1,  
										   BSTR pCcy,
										   BSTR pIndexPay, 					
										   BSTR pIndexPhase1, 										   
										   double pSpreadPhase1,
										   BSTR pDayCountPhase1,
										   BSTR pPayFreqPhase1,
										   BSTR pResetFreqPhase1,
										   BSTR pResetTimingPhase1,
										   BSTR pRoll, 
										   BSTR pAdjPhase1, 
										   BSTR pStub,
										   BSTR pIndexPhase2DIG,
										   BSTR pIndexLongPhase2DIG,
										   BSTR pStrikePhase2DIG,
										   BSTR pResetTimingPhase2DIG,
										   BSTR pAdjPhase2DIG, 
										   double pStartDatePhase2,
										   double pSpreadPhase2,										   
										   BSTR pDayCountPhase2, 
										   BSTR pPayFreqPhase2, 
										   BSTR pResetFreqPhase2,
										   BSTR pAdjPhase2,
										   double pStartDatePhase3,
										   double pEndDatePhase3,
										   BSTR pIndexPhase3,
										   double pSpreadPhase3, 
										   BSTR pDayCountPhase3, 
										   BSTR pPayFreqPhase3, 
										   BSTR pResetFreqPhase3,
										   BSTR pResetTimingPhase3,
										   BSTR pAdjPhase3,
										   BSTR pIndexFund,
										   VARIANT pSpreadFund,
										   BSTR pDayCountFund, 
										   BSTR pPayFreqFund, 
										   BSTR pResetFreqFund,
										   BSTR pResetTimingFund,
										   BSTR pAdjFund, 
										   double pEndDateAmort,
										   BSTR pDayCountAmort,
										   BSTR pIntRuleAmort, 
										   double pTxAmort,
										   BSTR pFreqAmort,
										   double pAmountAmort, 
										   BSTR pTypeAmort, 
										   double pFee, 
										   BSTR pVolCurvFromMatriceShift,
										   BSTR pVol, 
										   BSTR pVolCub,  
										   BSTR pConvexityManager, 
										   BSTR pZc, 
										   BSTR pSmiledMod, 
										   BSTR pSmiledModBump, 					
										   BSTR pHyperCubeCorrel,
										   BSTR pIndexIndexCorrelCube, 
										   BSTR pCorrEUR, 
										   BSTR pInterCorr, 
										   VARIANT *pRet
										   )

{	
	// TODO: Add your implementation code here
try
{								   
	ARM_result C_result;

	double startDatePhase1 = pStartDatePhase1;  
	double spreadPhase1 = pSpreadPhase1;

	double startDatePhase2 = pStartDatePhase2;
	double spreadPhase2 = pSpreadPhase2;

	double startDatePhase3 = pStartDatePhase3;
	double endDatePhase3 = pEndDatePhase3;
	double spreadPhase3 = pSpreadPhase3;

	double spreadFund;
	long   spreadType;

	if (VARIANT2Double(pSpreadFund,&spreadFund) == S_FALSE)
	{
		_bstr_t bSpreadFund(pSpreadFund);
		CCString l_SpreadFund= bSpreadFund;
		spreadFund = (double) LocalGetNumObjectId(l_SpreadFund);
		spreadType = 1;
	}
	else
	{
	   spreadType = 0;
	}

	double Ntl = pNtl;
	double fee = pFee; 

	long retCode;
	double Res = 0.;

	_bstr_t ccy(pCcy);
	CCString l_ccy = ccy;

	_bstr_t indexPay(pIndexPay);
	CCString l_indexPay = indexPay;

	_bstr_t indexPhase1(pIndexPhase1);
	CCString l_indexPhase1 = indexPhase1;

	_bstr_t dayCountPhase1(pDayCountPhase1);
	CCString l_dayCountPhase1 = dayCountPhase1;

	_bstr_t resetFreqPhase1(pResetFreqPhase1);
	CCString l_resetFreqPhase1 = resetFreqPhase1;

	_bstr_t payFreqPhase1(pPayFreqPhase1);
	CCString l_payFreqPhase1 = payFreqPhase1;

	_bstr_t resetTimingPhase1(pResetTimingPhase1);
	CCString l_resetTimingPhase1 = resetTimingPhase1;

	_bstr_t roll(pRoll);
	CCString l_roll = roll;

	_bstr_t adjPhase1(pAdjPhase1);
	CCString l_adjPhase1 = adjPhase1;
	
	_bstr_t stub(pStub);
	CCString l_stub= stub;

	_bstr_t indexPhase2DIG(pIndexPhase2DIG);
	CCString l_indexPhase2DIG = indexPhase2DIG;

	_bstr_t indexLongPhase2DIG(pIndexLongPhase2DIG);
	CCString l_indexLongPhase2DIG = indexLongPhase2DIG;

	//Strike is double??
	_bstr_t strikePhase2DIG(pStrikePhase2DIG);
	CCString l_strikePhase2DIG= strikePhase2DIG;

	_bstr_t resetTimingPhase2DIG(pResetTimingPhase2DIG);
	CCString l_resetTimingPhase2DIG = resetTimingPhase2DIG;

	_bstr_t adjPhase2DIG(pAdjPhase2DIG);
	CCString l_adjPhase2DIG = adjPhase2DIG;

	_bstr_t dayCountPhase2(pDayCountPhase2);
	CCString l_dayCountPhase2= dayCountPhase2;

	_bstr_t payFreqPhase2(pPayFreqPhase2);
	CCString l_payFreqPhase2= payFreqPhase2;
	
	_bstr_t resetFreqPhase2(pResetFreqPhase2);
	CCString l_resetFreqPhase2 = resetFreqPhase2;

	_bstr_t adjPhase2(pAdjPhase2);
	CCString l_adjPhase2 = adjPhase2;

	_bstr_t indexPhase3(pIndexPhase3);
	CCString l_indexPhase3 = indexPhase3;

	_bstr_t dayCountPhase3(pDayCountPhase3);
	CCString l_dayCountPhase3 = dayCountPhase3;

	_bstr_t payFreqPhase3(pPayFreqPhase3);
	CCString l_payFreqPhase3 = payFreqPhase3;

	_bstr_t resetFreqPhase3(pResetFreqPhase3);
	CCString l_resetFreqPhase3 = resetFreqPhase3;

	_bstr_t resetTimingPhase3(pResetTimingPhase3);
	CCString l_resetTimingPhase3 = resetTimingPhase3;

	_bstr_t adjPhase3(pAdjPhase3);
	CCString l_adjPhase3 = adjPhase3;

	_bstr_t indexFund(pIndexFund);
	CCString l_indexFund = indexFund;

	_bstr_t dayCountFund(pDayCountFund);
	CCString l_dayCountFund= dayCountFund;

	_bstr_t payFreqFund(pPayFreqFund);
	CCString l_payFreqFund= payFreqFund;
	
	_bstr_t resetFreqFund(pResetFreqFund);
	CCString l_resetFreqFund = resetFreqFund;

	_bstr_t resetTimingFund(pResetTimingFund);
	CCString l_resetTimingFund = resetTimingFund;

	_bstr_t adjFund(pAdjFund);
	CCString l_adjFund= adjFund;

	_bstr_t typeAmort(pTypeAmort);
	CCString l_typeAmort = typeAmort;

	CCString l_floorOrCap = "C"; 

	_bstr_t volCurvFromMatriceShift(pVolCurvFromMatriceShift);
	CCString l_volCurvFromMatriceShift= volCurvFromMatriceShift;
	
	_bstr_t vol(pVol);
	CCString l_vol = vol;

	_bstr_t volCub(pVolCub);
	CCString l_volCub = volCub;

	_bstr_t convexityManager(pConvexityManager);
	CCString l_convexityManager = convexityManager;

	_bstr_t zc(pZc);
	CCString l_zc = zc;

	_bstr_t smiledMod (pSmiledMod);
	CCString l_smiledMod = smiledMod;

	_bstr_t smiledModBump(pSmiledModBump); 
	CCString l_smiledModBump = smiledModBump;	

	_bstr_t hyperCubeCorrel (pHyperCubeCorrel);
	CCString l_hyperCubeCorrel = hyperCubeCorrel; 

	_bstr_t indexIndexCorrelCube (pIndexIndexCorrelCube);
	CCString l_indexIndexCorrelCube = indexIndexCorrelCube;	
	
	_bstr_t corrEUR (pCorrEUR);
	CCString l_corrEUR = corrEUR;	

	_bstr_t interCorr (pInterCorr);
	CCString l_interCorr= interCorr;

	long payFreqPhase1Id;	
	if((payFreqPhase1Id = ARM_ConvFrequency(l_payFreqPhase1, C_result)) == ARM_DEFAULT_ERR)
	{
		ERROR_MSG("Invalid Frq",pRet,ARM_ERROR_FREQ);
		return S_OK;
	}

	long resetFreqPhase1Id;
	if((resetFreqPhase1Id = ARM_ConvFrequency(l_resetFreqPhase1, C_result)) == ARM_DEFAULT_ERR)
	{
		ERROR_MSG("Invalid Frq",pRet,ARM_ERROR_FREQ);
		return S_OK;
	}
	
	long payFreqPhase2Id;
	if((payFreqPhase2Id = ARM_ConvFrequency(l_payFreqPhase2, C_result)) == ARM_DEFAULT_ERR)
	{
		ERROR_MSG("Invalid Frq",pRet,ARM_ERROR_FREQ);
		return S_OK;
	}

	long resetFreqPhase2Id;
	if((resetFreqPhase2Id = ARM_ConvFrequency(l_resetFreqPhase2, C_result)) == ARM_DEFAULT_ERR)
	{
		ERROR_MSG("Invalid Frq",pRet,ARM_ERROR_FREQ);
		return S_OK;
	}

	long payFreqPhase3Id;
	if((payFreqPhase3Id = ARM_ConvFrequency(l_payFreqPhase3, C_result)) == ARM_DEFAULT_ERR)
	{
		ERROR_MSG("Invalid Frq",pRet,ARM_ERROR_FREQ);
		return S_OK;
	}

	long resetFreqPhase3Id;
	if((resetFreqPhase3Id = ARM_ConvFrequency(l_resetFreqPhase3, C_result)) == ARM_DEFAULT_ERR)
	{
		ERROR_MSG("Invalid Frq",pRet,ARM_ERROR_FREQ);
		return S_OK;
	}

	long payFreqFundId;
	if((payFreqFundId = ARM_ConvFrequency(l_payFreqFund, C_result)) == ARM_DEFAULT_ERR)
	{
		ERROR_MSG("Invalid Frq",pRet,ARM_ERROR_FREQ);
		return S_OK;
	}

	long resetFreqFundId;
	if((resetFreqFundId = ARM_ConvFrequency(l_resetFreqFund, C_result)) == ARM_DEFAULT_ERR)
	{
		ERROR_MSG("Invalid Frq",pRet,ARM_ERROR_FREQ);
		return S_OK;
	}

	long indexPhase1Id; 
	indexPhase1Id = ARM_ConvIrType(l_indexPhase1); 

	long indexPhase3Id;
	indexPhase3Id = ARM_ConvIrType(l_indexPhase3);
	
	long indexFundId;
	indexFundId = ARM_ConvIrType(l_indexFund);

	long rollId;
	if((rollId = ARM_ConvFwdRule(l_roll, C_result)) == ARM_DEFAULT_ERR)
	{
		ERROR_MSG("Invalid Int Roll",pRet,ARM_ERROR_FREQ);
		return S_OK;
	}	
	
	long adjPhase1Id = ARM_ConvIntRule(l_adjPhase1);
	long adjPhase2Id = ARM_ConvIntRule(l_adjPhase2);
	long adjPhase2DIGId = ARM_ConvIntRule(l_adjPhase2DIG);
	long adjPhase3Id = ARM_ConvIntRule(l_adjPhase3);
	long adjFundId = ARM_ConvIntRule(l_adjFund);

	long resetTimingPhase1Id = ARM_ConvPayResetRule (l_resetTimingPhase1);
	long resetTimingPhase2DIGId = ARM_ConvPayResetRule (l_resetTimingPhase2DIG);
	long resetTimingPhase3Id = ARM_ConvPayResetRule (l_resetTimingPhase3);
	long resetTimingFundId = ARM_ConvPayResetRule (l_resetTimingFund);

	CCString l_payTimingPhase2DIG = "ARR"; 
	long payTimingPhase2DIGId = ARM_ConvPayResetRule(l_resetTimingPhase2DIG);

	long dayCountPhase1Id= ARM_ConvDayCount(l_dayCountPhase1);
	long dayCountPhase2Id= ARM_ConvDayCount(l_dayCountPhase2);
	long dayCountPhase3Id= ARM_ConvDayCount(l_dayCountPhase3);
	long dayCountFundId= ARM_ConvDayCount(l_dayCountFund);

	long stubId = ARM_ConvStubRule(l_stub);
	long ccyId; 

	long resetTimingPhase2Id = 1;

	if (ARMLOCAL_ISOCCY(l_ccy,C_result) == ARM_OK)
	{
		ccyId = C_result.getLong();
	}
	else
	{
		ERROR_MSG(C_result.getMsg(),pRet);
		return S_OK;
	}

	//************************************************************************************************
	//*********************************** Création Amort Refvalue ************************************
	//************************************************************************************************
	long amortRefValue;
	long amortRefValueToDestroy = ARM_NULL_OBJECT;

	if (LocalGetNumObjectId(l_typeAmort) == ARM_KO)
	{
		double endDateAmort = pEndDateAmort;
		double txAmort = pTxAmort;
		double amountAmort = pAmountAmort;

		_bstr_t intRuleAmort(pIntRuleAmort);
		CCString l_intRuleAmort =  intRuleAmort;

		_bstr_t freqAmort(pFreqAmort);
		CCString l_freqAmort = freqAmort;

		_bstr_t dayCountAmort(pDayCountAmort);
		CCString l_dayCountAmort = dayCountAmort;

		long freqAmortId; 
		if((freqAmortId = ARM_ConvFrequency(l_freqAmort, C_result)) == ARM_DEFAULT_ERR)
		{
			ERROR_MSG("Invalid Frq",pRet,ARM_ERROR_FREQ);
			return S_OK;
		}

		long intRuleAmortId= ARM_ConvIntRule(l_intRuleAmort);
		long dayCountAmortId = ARM_ConvDayCount(l_dayCountAmort);	

		// AmortMethodId
		long amortMethodId;
		if((amortMethodId = ARM_ConvAmortMethod (l_typeAmort, C_result)) == ARM_DEFAULT_ERR)
		{
			ERROR_MSG("Invalid TypeAmort",pRet,ARM_ERROR_FREQ);
			return S_OK;
		}

		long FixedLegAmortId;
		if (ARMLOCAL_FIXEDLEG(startDatePhase1,
							  endDateAmort,
							  K_RCV,
							  0,
							  1.0,
							  dayCountAmortId,
							  freqAmortId,
							  K_COMP_PROP,
							  K_ARREARS,
							  intRuleAmortId,
							  stubId,
							  false,
							  l_ccy,
							  "NULL",
							  K_NX_NONE,
							  -1.0,
							  K_NO,
								 -1,
							  C_result) == ARM_OK) 
							  
		{
			FixedLegAmortId = C_result.getLong();
		}
		else
		{
			ERROR_MSG(C_result.getMsg(),pRet);
			return S_OK;
		}

		if (ARMLOCAL_ARM_GenAmortization(FixedLegAmortId,
										 amortMethodId,
										 freqAmortId,
										 amountAmort,
										 -1,
										 Ntl,
										 txAmort,
										 0.,
										 ARM_NULL_OBJECT,
										 0.0,
										 C_result) == ARM_OK) 
		{
			amortRefValue = C_result.getLong();
			amortRefValueToDestroy = amortRefValue;
		}
		else
		{
			ERROR_MSG(C_result.getMsg(),pRet);
			return S_OK;
		}

		retCode = ARMLOCAL_FreeObject(FixedLegAmortId,C_result);
	}
	else
	{
		amortRefValue = LocalGetNumObjectId(l_typeAmort);

		// calcul du 1er nominal
		double dPaymentDate;

		if (ARMLOCAL_DisplayRefValue(amortRefValue,true,C_result) == ARM_OK)
		{
			dPaymentDate = C_result.getArray(0);
		}
		else
		{
			ERROR_MSG(C_result.getMsg(),pRet);
			return S_OK;
		}
		
		if (ARMLOCAL_CptRefValue(amortRefValue,dPaymentDate,C_result) == ARM_OK)
		{
			Ntl = C_result.getDouble();
		}
		else
		{
			ERROR_MSG(C_result.getMsg(),pRet);
			return S_OK;
		} 
	}
	
	//************************************************************************************************
	//*********************************** Fin Création Amort Refvalue ************************************
	//************************************************************************************************


	//************************************************************************************************
	//*********************************** Création Funding Leg Refvalue ******************************
	//************************************************************************************************

	long irIndexFundId;
	if (ARMLOCAL_IRINDEX(KACTUAL_360,
						 payFreqFundId,
						 -1,
						 K_COMP_PROP,
						 K_MOD_FOLLOWING,
						 K_ADVANCE,
						 -2,
						 K_ARREARS,
						 10000.,
						 false,
						 l_ccy,
						 indexFundId,
						 K_COMP_PROP,
						 adjFundId,
						 resetFreqFundId,
						 C_result) == ARM_OK)
	{
		irIndexFundId = C_result.getLong();
	}
	else
	{
		ERROR_MSG(C_result.getMsg(),pRet);
		return S_OK;
	}

	long swapLegId;
	if (ARMLOCAL_SWAPLEG(irIndexFundId,
						 startDatePhase1,
						 endDatePhase3,
						 K_RCV,
						 spreadType, 
						 spreadFund, 
						 false,
						 l_ccy,
						 dayCountFundId,
						 -2,
						 "NULL",
						 "NULL",
						 0,
						 K_NX_NONE,
						 stubId,
						 -1.0,
						 K_NO,
						  -1,
						 C_result) == ARM_OK) 
	{
		swapLegId = C_result.getLong();
	}
	else
	{
		ERROR_MSG(C_result.getMsg(),pRet);
		return S_OK;
	}
		
	// ********** Funding Leg ************
	long cloneId;
	if (ARMLOCAL_ClonedAndSetNotional(swapLegId,
									  amortRefValue,
									  100.,
									  C_result) == ARM_OK)
	{
		cloneId = C_result.getLong();
	}
	else
	{
		ERROR_MSG(C_result.getMsg(),pRet);
		return S_OK;
	}

	//*********   Price funding Leg   *************************
	double pxFund; 
	double pxFundDelta;

	if (ARMLOCAL_ARM_Price(cloneId,
						   LocalGetNumObjectId(l_smiledMod),
						   C_result) == ARM_OK)
	{
		pxFund = C_result.getDouble();
	}
	else
	{
		ERROR_MSG(C_result.getMsg(),pRet);
		return S_OK;
	}

	// shift seulement sur l'euro !?
	if (ARMLOCAL_ARM_Price(cloneId,
						   LocalGetNumObjectId(l_smiledModBump),
						   C_result) == ARM_OK)
	{
		pxFundDelta = C_result.getDouble();
	}
	else
	{
		ERROR_MSG(C_result.getMsg(),pRet);
		return S_OK;
	}

	retCode = ARMLOCAL_FreeObject(irIndexFundId, C_result);
	retCode = ARMLOCAL_FreeObject(swapLegId, C_result);
	retCode = ARMLOCAL_FreeObject(cloneId,C_result);
	
	// Fin du calcul du prix de la funding Leg

	//************************************************************************************************
	//*********************************** Fin Création Funding Leg Refvalue **************************
	//************************************************************************************************





	
	//************************************************************************************************
	//*********************************** Création Leg Phase1 Refvalue *******************************
	//************************************************************************************************

	long irIndexPhase1Id;
	double pxPhase1; 
	double pxPhase1Delta;
	long swapLegPhase1Id;
	long clonePhase1Id;
	double resetGapPhase1; 


	if (resetTimingPhase1Id == 1)
	{
		resetGapPhase1 = -2.; 
	}
	else
	{
		resetGapPhase1 = -15.; 
	}

	if ( startDatePhase1 != startDatePhase2 )
	{
		if (ARMLOCAL_IRINDEX(KACTUAL_360,
							 payFreqPhase1Id,
							 -1,
							 K_COMP_PROP,
							 K_MOD_FOLLOWING,
							 resetTimingPhase1Id,	
							 resetGapPhase1,
							 K_ARREARS,
							 10000.,
							 false,
							 l_ccy,
							 indexPhase1Id,
							 K_COMP_PROP,
							 adjPhase1Id,
							 resetFreqPhase1Id,
							 C_result) == ARM_OK) 
		{
			irIndexPhase1Id = C_result.getLong();
		}
		else
		{
			ERROR_MSG(C_result.getMsg(),pRet);
			return S_OK;
		}

		if (ARMLOCAL_SWAPLEG(irIndexPhase1Id,
							 startDatePhase1,
							 startDatePhase2,
							 K_RCV,
							 0, 
							 spreadPhase1, 
							 false,
							 l_ccy,
							 dayCountPhase1Id,
							 resetGapPhase1,
							 "NULL",
							 "NULL",
							 0,
							 K_NX_NONE,
							 stubId,
							 -1.0,
							 K_NO,
						  -1,
							 C_result) == ARM_OK) 
		{
			swapLegPhase1Id = C_result.getLong();
		}
		else
		{
			ERROR_MSG(C_result.getMsg(),pRet);
			return S_OK;
		}
		
		// ********** Phase1 Leg ************
		if (ARMLOCAL_ClonedAndSetNotional(swapLegPhase1Id,
										  amortRefValue,
										  100.,
										  C_result) == ARM_OK)
		{
			clonePhase1Id = C_result.getLong();
		}
		else
		{
			ERROR_MSG(C_result.getMsg(),pRet);
			return S_OK;
		}

		//*********   Price Phase1 *************************
		if (ARMLOCAL_ARM_Price(clonePhase1Id,
							   LocalGetNumObjectId(l_smiledMod),
							   C_result) == ARM_OK)
		{
			pxPhase1 = C_result.getDouble();
		}
		else
		{
			ERROR_MSG(C_result.getMsg(),pRet);
			return S_OK;
		}

		// shift seulement sur l'euro !?
		if (ARMLOCAL_ARM_Price(clonePhase1Id,
							   LocalGetNumObjectId(l_smiledModBump),
							   C_result) == ARM_OK)
		{
			pxPhase1Delta = C_result.getDouble();
		}
		else
		{
			ERROR_MSG(C_result.getMsg(),pRet);
			return S_OK;
		}

		retCode = ARMLOCAL_FreeObject(irIndexPhase1Id,C_result);
		retCode = ARMLOCAL_FreeObject(swapLegPhase1Id,C_result);
		retCode = ARMLOCAL_FreeObject(clonePhase1Id,C_result);
	}
	else
	{
		pxPhase1 = 0.;
		pxPhase1Delta = 0.;
	}
	//************************************************************************************************
	//************************** Fin Création Leg Phase1 Refvalue ************************************
	//************************************************************************************************


	//************************************************************************************************
	//*********************************** Création Leg Phase2 Refvalue *******************************
	//************************************************************************************************

	long irIndexPhase2Id;
	double resetGapPhase2; 

	if (resetTimingPhase2Id == 1)
	{
		resetGapPhase2 = -2.; 
	}
	else
	{
		resetGapPhase2 = -15.; 
	}

	if (ARMLOCAL_IRINDEX(KACTUAL_360,
						 payFreqPhase2Id,
						 -1,
						 K_COMP_PROP,
						 K_MOD_FOLLOWING,
						 resetTimingPhase2Id,	
						 resetGapPhase2,
						 K_ARREARS,
						 10000.,
						 false,
						 l_ccy,
						 K_FIXED,
						 K_COMP_PROP,
						 adjPhase2Id,
						 resetFreqPhase2Id,
						 C_result) == ARM_OK) 
	{
		irIndexPhase2Id = C_result.getLong();
	}
	else
	{
		ERROR_MSG(C_result.getMsg(),pRet);
		return S_OK;
	}

	long swapLegPhase2Id;
	if (ARMLOCAL_SWAPLEG(irIndexPhase2Id,
						 startDatePhase2,
						 startDatePhase3,
						 K_RCV,
						 0,
						 spreadPhase2,
						 false,
						 l_ccy,
						 dayCountPhase2Id,
						 resetGapPhase2,
						 "NULL",
						 "NULL",
						 0, 
						 K_NX_NONE,
						 stubId, 
						 -1.0,
						 K_NO,
						  -1,
						 C_result) == ARM_OK) 
	{
		swapLegPhase2Id = C_result.getLong();
	}
	else
	{
		ERROR_MSG(C_result.getMsg(),pRet);
		return S_OK;
	}
	
	// ********** Phase2 Leg ************
	long clonePhase2Id;
	if (ARMLOCAL_ClonedAndSetNotional(swapLegPhase2Id,
									  amortRefValue,
									  100.,
									  C_result) == ARM_OK)
	{
		clonePhase2Id = C_result.getLong();
	}
	else
	{
		ERROR_MSG(C_result.getMsg(),pRet);
		return S_OK;
	}

	//*********   Price Phase2 *************************
	double pxPhase2;
	double pxPhase2Delta; 

	if (ARMLOCAL_ARM_Price(clonePhase2Id,
						   LocalGetNumObjectId(l_smiledMod),
						   C_result) == ARM_OK)
	{
		pxPhase2 = C_result.getDouble();
	}
	else
	{
		ERROR_MSG(C_result.getMsg(),pRet);
		return S_OK;
	}

	// shift seulement sur l'euro !?
	if (ARMLOCAL_ARM_Price(clonePhase2Id,
						   LocalGetNumObjectId(l_smiledModBump),
						   C_result) == ARM_OK)
	{
		pxPhase2Delta = C_result.getDouble();
	}
	else
	{
		ERROR_MSG(C_result.getMsg(),pRet);
		return S_OK;
	}

	retCode = ARMLOCAL_FreeObject(irIndexPhase2Id,C_result);
	retCode = ARMLOCAL_FreeObject(swapLegPhase2Id,C_result);
	retCode = ARMLOCAL_FreeObject(clonePhase2Id,C_result);

	//************************************************************************************************
	//************************** Fin Création Leg Phase2 Refvalue ************************************
	//************************************************************************************************



	//************************************************************************************************
	//*************************** Création Leg Phase3 Refvalue ***************************************
	//************************************************************************************************

	long irIndexPhase3Id;
	double pxPhase3; 
	double pxPhase3Delta;
	long swapLegPhase3Id;
	long clonePhase3Id;
	double resetGapPhase3; 

	if ( startDatePhase3 != endDatePhase3 )
	{		
		if (resetTimingPhase3Id == 1)
		{
			resetGapPhase3 = -2.; 
		}
		else
		{
			resetGapPhase3 = -15.; 
		}

		if (ARMLOCAL_IRINDEX(KACTUAL_360,
							 payFreqPhase3Id,
							 -1,
							 K_COMP_PROP,
							 K_MOD_FOLLOWING,
							 resetTimingPhase3Id,	
							 resetGapPhase3,
							 K_ARREARS,
							 10000.,
							 false,
							 l_ccy,
							 indexPhase3Id,
							 K_COMP_PROP,
							 adjPhase3Id,
							 resetFreqPhase3Id,
							 C_result) == ARM_OK) 
		{
			irIndexPhase3Id = C_result.getLong();
		}
		else
		{
			ERROR_MSG(C_result.getMsg(),pRet);
			return S_OK;
		}

		if (ARMLOCAL_SWAPLEG(irIndexPhase3Id,
							 startDatePhase3,
							 endDatePhase3,
							 K_RCV,
							 0, 
							 spreadPhase3, 
							 false,
							 l_ccy,
							 dayCountPhase3Id,
							 resetGapPhase3,
							 "NULL",
							 "NULL",
							 0,
							 K_NX_NONE,
							 stubId,
							 -1.0,
							 K_NO,
						  -1,
							 C_result) == ARM_OK) 
		{
			swapLegPhase3Id = C_result.getLong();
		}
		else
		{
			ERROR_MSG(C_result.getMsg(),pRet);
			return S_OK;
		}
			
		// ********** Phase3 Leg ************
		if (ARMLOCAL_ClonedAndSetNotional(swapLegPhase3Id,
										  amortRefValue,
										  100.,
										  C_result) == ARM_OK)
		{
			clonePhase3Id = C_result.getLong();
		}
		else
		{
			ERROR_MSG(C_result.getMsg(),pRet);
			return S_OK;
		}

		//*********   Price Phase3 *************************	
		if (ARMLOCAL_ARM_Price(clonePhase3Id,
							   LocalGetNumObjectId(l_smiledMod),
							   C_result) == ARM_OK)
		{
			pxPhase3 = C_result.getDouble();
		}
		else
		{
			ERROR_MSG(C_result.getMsg(),pRet);
			return S_OK;
		}

		// shift seulement sur l'euro !?
		if (ARMLOCAL_ARM_Price(clonePhase3Id,
							   LocalGetNumObjectId(l_smiledModBump),
							   C_result) == ARM_OK)
		{
			pxPhase3Delta = C_result.getDouble();
		}
		else
		{
			ERROR_MSG(C_result.getMsg(),pRet);
			return S_OK;
		}
		retCode = ARMLOCAL_FreeObject(irIndexPhase3Id,C_result);
		retCode = ARMLOCAL_FreeObject(swapLegPhase3Id,C_result);
		retCode = ARMLOCAL_FreeObject(clonePhase3Id,C_result);
	}
	else
	{
		pxPhase3 = 0.;
		pxPhase3Delta = 0.; 
	}

	double pxPhase2DIG; 
	double pxTotal; 
	double pxDelta; 
	double SensiTotal; 	
	double prixInter;
	double Xincr;
	
	double resetGapPhase2DIG; 

	if (resetTimingPhase2DIGId == 1)
	{
		resetGapPhase2DIG = -2.; 
	}
	else
	{
		resetGapPhase2DIG = -15.; 
	}

	//Construction de l'index PayOff
	long indexPayId; 
	indexPayId = ARM_ConvIrType(l_indexPay); 


//******************** Calcul premiere PV ******************
	int iteration =1;
	double prix;
	double prixAct; 
	double prixAct2;
	double SprSeek = 10.01;
	double SprSeekFinal= 10.01;
	double SprSeekBonifie = 10.01;  // PUREMENT ARBITRAIRE
	fee = Ntl * fee / 100;
//******************** Calcul premiere PV ******************

	VECTOR<double> pxPentilix = computePentilix(l_indexPhase2DIG,l_indexLongPhase2DIG, indexPayId, startDatePhase2, startDatePhase3, 
			l_floorOrCap, l_strikePhase2DIG, SprSeekBonifie, ccyId, 
			dayCountPhase2Id, 
			resetTimingPhase2DIGId, payTimingPhase2DIGId,  resetGapPhase2DIG, resetFreqPhase2Id, payFreqPhase2Id, 
			adjPhase2DIGId, stubId,amortRefValue, 
			l_hyperCubeCorrel, l_volCurvFromMatriceShift,
			l_vol, l_volCub, 
			l_convexityManager, l_zc, l_smiledMod, l_indexIndexCorrelCube, l_corrEUR,l_interCorr);


	if (pxPentilix.size() == 0)
	{
		ERROR_MSG("Pb in computePentilix",pRet);
		return S_OK;
	}

	pxPhase2DIG = pxPentilix[0];
	double pxPhase2DIGBumpCorrel = pxPentilix[1];
	double MT = pxPentilix[2];

	pxTotal = pxPhase1 + pxPhase2 + pxPhase3 - pxFund + pxPhase2DIG;
	pxDelta = fabs(pxPhase2DIGBumpCorrel-pxPhase2DIG) + MT; 
	if (pxDelta <0.) pxDelta = - pxDelta;

	SensiTotal = pxDelta;
	prixAct = pxTotal - SensiTotal - fee; 

	prix = prixAct ;
	prixAct2 =1;
	
//********************                    ******************
//******************** BOUCLE SUR LE PRICE ******************
	while ((abs(prixAct2) > 0.0001) && (iteration < 10 ))
	{
		SprSeek = SprSeekBonifie;
		
		//*************************
		// **** Premier recalcul ****
		//****************************

		pxPentilix = computePentilix(l_indexPhase2DIG,l_indexLongPhase2DIG, indexPayId, startDatePhase2, startDatePhase3, 
			l_floorOrCap, l_strikePhase2DIG, SprSeekBonifie , ccyId, 
			dayCountPhase2Id, 
			resetTimingPhase2DIGId, payTimingPhase2DIGId,  resetGapPhase2DIG, resetFreqPhase2Id, payFreqPhase2Id, 
			adjPhase2DIGId, stubId, amortRefValue,
			l_hyperCubeCorrel, l_volCurvFromMatriceShift,
			l_vol, l_volCub,
			l_convexityManager, l_zc, l_smiledMod, l_indexIndexCorrelCube, l_corrEUR,l_interCorr);

		if (pxPentilix.size() == 0)
		{
			ERROR_MSG("Pb in computePentilix",pRet);
			return S_OK;
		}
  
		pxPhase2DIG = pxPentilix[0];
		pxPhase2DIGBumpCorrel = pxPentilix[1];
		MT = pxPentilix[2];

		pxTotal = pxPhase1 + pxPhase2 + pxPhase3 - pxFund + pxPhase2DIG;
		pxDelta = fabs(pxPhase2DIGBumpCorrel-pxPhase2DIG) + MT; 
		if (pxDelta <0.) pxDelta = - pxDelta;

		SensiTotal = pxDelta;
		prixAct = pxTotal - SensiTotal - fee; 

		// les deltas se compensent et on somme en valeur absolue les deltas et vegas

		//******************************
		// **** affectation valeurs ****
		//******************************

		prixAct2 = prixAct ;
		Xincr = SprSeek / 1000000;
		SprSeekBonifie = SprSeek + Xincr;		

		
		//****************************
		// **** Deuxieme recalcul ****
		//****************************
		pxPentilix = computePentilix(l_indexPhase2DIG,l_indexLongPhase2DIG, indexPayId, startDatePhase2, startDatePhase3, 
			l_floorOrCap, l_strikePhase2DIG, SprSeekBonifie , ccyId, 
			dayCountPhase2Id, 
			resetTimingPhase2DIGId, payTimingPhase2DIGId,  resetGapPhase2DIG, resetFreqPhase2Id, payFreqPhase2Id, 
			adjPhase2DIGId, stubId, amortRefValue, 
			l_hyperCubeCorrel, l_volCurvFromMatriceShift, 
			l_vol, l_volCub, 
			l_convexityManager, l_zc, l_smiledMod, l_indexIndexCorrelCube, l_corrEUR,l_interCorr);

		if (pxPentilix.size() == 0)
		{
			ERROR_MSG("Pb in computePentibor",pRet);
			return S_OK;
		}

		pxPhase2DIG = pxPentilix[0];
		pxPhase2DIGBumpCorrel = pxPentilix[1];
		MT = pxPentilix[2];

		pxTotal = pxPhase1 + pxPhase2 + pxPhase3 - pxFund + pxPhase2DIG;
		pxDelta = fabs(pxPhase2DIGBumpCorrel-pxPhase2DIG) + MT; 
		if (pxDelta <0.) pxDelta = - pxDelta;

		SensiTotal = pxDelta;
		prixAct = pxTotal - SensiTotal - fee; 
	
		iteration++;
		
		SprSeekFinal = SprSeekBonifie;
		prixInter = ( prixAct - prixAct2 ) / Xincr;
		// if prixinter =  0 sortir erreur
		SprSeekBonifie = SprSeek - prixAct / prixInter ; 
		prixAct2 = prixAct;	
	}
		
//**************** PRICE *******************

	retCode = ARMLOCAL_FreeObject(amortRefValueToDestroy,C_result);
	retCode = ARMLOCAL_FreeObject(ccyId,C_result);

	VECTOR<double> vRes;

	vRes.push_back(spreadPhase2 + SprSeekFinal);
	vRes.push_back(spreadPhase2 + SprSeekFinal);

	Res = VECTORDOUBLE2VARIANT(vRes, pRet);

	_variant_t wrap_Res;
	wrap_Res.Attach(*pRet);
	*pRet=wrap_Res.Detach();

	return S_OK;
}
catch(...)
{
	return E_FAIL;
}
}

VECTOR<double> computePentilix(const CCString& index, const CCString& indexLong, long indexPayId, double startDate, double endDate,
							   const CCString& floorOrCap, const CCString& strike, double payOffSpread, long ccyId, long dayCountId,
							   long resetTimingId, long payTimingId, double resetGap, long resetFreqId, long payFreqId,
							   long intRuleId, long stubRuleId, long amortRefvalue,
							   const CCString& hyperCubeCorrel, const CCString& VolCurvFromMatriceShift,
							   const CCString& vol, const CCString& volCUB, 
							   const CCString& ConvexityManager, const CCString& zc, const CCString& smiledMod, const CCString& indexIndexCorrelCube, 
							   const CCString& corrEUR, const CCString& InterCorr)
{
	VECTOR<double> resultat;
	resultat.clear();

try
{
	double pxDigital;
	double pxBumpDigital; 
	double pxPentilix;
	double pxBumpPentilix; 
	double pxTotal; 

	long DigitalId;
	long cloneDigitalId; 

	long PentilixId;
	long clonePentilixId; 

	long corrManagerId; 

	double shiftCorrel; 
	
	long bumpHyperCubeCorrelId;
	long bumpCorrManagerId; 

	ARM_result C_result;

	// Modèles
	long bsModGenId; 
	long BumpBsModGenId;

	VECTOR <double> spreadVector;
	VECTOR <double> weightVector;

   	VECTOR<double> C_fixing1_double_vect;
   	VECTOR<double> C_fixing2_double_vect;
	VECTOR<double> C_fixing3_double_vect;

	long slopeFlag;
	int cptStrikeMethod;
	int computedFormula;

	slopeFlag = 1;
	cptStrikeMethod = 1.;
	computedFormula = 1.;

	double d_index; 
	long indexId; 
	indexId = ARM_ConvIrType(index); 

	if (IsCMSIndex((ARM_INDEX_TYPE)indexId) == 1)
	{
		d_index = (indexId - K_CMS1 +1 ); 
	}	
	else
	{
		d_index = 1./FromLiborTypeToFrequency(indexId); 
	}
	
	double d_indexLong; 
	long indexLongId; 
	indexLongId = ARM_ConvIrType(indexLong); 

	if (IsCMSIndex((ARM_INDEX_TYPE)indexLongId) ==1)
	{
		d_indexLong = indexLongId - K_CMS1 + 1; 
	}	
	else
	{
		d_indexLong = 1./FromLiborTypeToFrequency(indexLongId); 
	}

	//Définition du sreadVector et du WeightVector
	spreadVector.push_back(-0.01);
	spreadVector.push_back(0.01); 

	weightVector.push_back(1);
	weightVector.push_back(1);

	long capOrFloorId; 
	if((capOrFloorId = ARM_ConvCapOrFloor(floorOrCap, C_result)) == ARM_DEFAULT_ERR)
	{
		return resultat;
	}
	
	// Les variables indexIndexCorrelCube, corrEUR, InterCorr sont à créer
	// La fonction ARMLOCAL_CreateGenCorrelatorManager est en attente de modifcation par JPR
 
	VECTOR<CCString> vCcy; 
	VECTOR<long> vHyperCubeCorrelId; 
	VECTOR<long> vIndexIndexCorrelCubeId; 
	VECTOR<long> vCorrEURId; 
	VECTOR<long> vInterCorrId; 
	VECTOR<long> vIrVolId; 
	VECTOR<long> vVolVolId; 
	VECTOR<long> vFXVolHyperCubeId; 

	vCcy.push_back("EUR"); 
	vHyperCubeCorrelId.push_back(LocalGetNumObjectId(hyperCubeCorrel));
	vIndexIndexCorrelCubeId.push_back(LocalGetNumObjectId(indexIndexCorrelCube));
	vInterCorrId.push_back(LocalGetNumObjectId(InterCorr));

	if (ARMLOCAL_CreateGenCorrelatorManager(vCcy, vHyperCubeCorrelId, 
		 vIndexIndexCorrelCubeId, vCorrEURId, vInterCorrId, vIrVolId, vVolVolId, vFXVolHyperCubeId, C_result) ==  ARM_OK)    
	{
		corrManagerId = C_result.getLong();
	}
	else
	{
		return resultat;
	}

	long modelTypeId;
	long spreadVolTypeId;
	CCString modelType; 
	CCString spreadVolType;
	bool isLnVol=true;
	long payMargin_type; 
	double payMarginId;
	double payIdxFixedRateId; 
	long payIdxFixedRate_type; 
	long payIdxWeight; 
	double freezeFixing = 0.0;

	payMargin_type = 0; 
	payMarginId = payOffSpread;
	payIdxFixedRateId = 0.; 
	payIdxFixedRate_type = 0;
	payIdxWeight = 1; 

	modelType = "2LOG";
	spreadVolType = "C"; 

	if ((modelTypeId = ARM_ConvModelType (modelType, C_result)) == ARM_DEFAULT_ERR)
	{
		return resultat;
	}

	if ((spreadVolTypeId = ARM_ConvVolType2 (spreadVolType, C_result)) == ARM_DEFAULT_ERR)
	{
		return resultat;
	}

	VECTOR<double> C_rate; 
	VECTOR<CCString> C_matu;
	double C_meth = 0;

	//************************************************************************************************
	//*************************** Création Corridor **************************************************
	//************************************************************************************************

	//Calcul du prix de la Digitale

	if (ARMLOCAL_SPREADDIGITALFLT(startDate,
						 endDate,
						 capOrFloorId,
						 1, 
						 LocalGetNumObjectId(strike), 
						 indexPayId, 
						 indexId, 
						 indexLongId,
						 weightVector[0], 
						 weightVector[1],
						 dayCountId,
						 resetFreqId, 
						 payFreqId, 
						 resetTimingId,
						 payTimingId,				
						 ccyId,
						 resetGap, 
						 spreadVector[0], 
						 spreadVector[1],
						 1,
						 ARM_NULL_OBJECT,
						 C_fixing1_double_vect,
						 1,
						 ARM_NULL_OBJECT,
						 C_fixing2_double_vect,
						 intRuleId, 
						 stubRuleId, 
						 slopeFlag,
						 cptStrikeMethod, 
						 computedFormula, 					
						 C_result) == ARM_OK)
	{
		PentilixId = C_result.getLong();
	}
	else
	{
		return resultat;
	}

	if (ARMLOCAL_ClonedAndSetNotional(PentilixId,
									  amortRefvalue,
									  100.,
									  C_result) == ARM_OK)
	{
		clonePentilixId = C_result.getLong();
	}
	else
	{
		return resultat;
	}

	long retCode = ARMLOCAL_FreeObject(PentilixId,C_result);

	if (ARMLOCAL_BSMODELGEN(LocalGetNumObjectId(zc), ARM_NULL_OBJECT, LocalGetNumObjectId(vol), LocalGetNumObjectId(volCUB), 
				corrManagerId, LocalGetNumObjectId(ConvexityManager), 
				LocalGetNumObjectId(zc), LocalGetNumObjectId(hyperCubeCorrel), ARM_NULL_OBJECT, ARM_NULL_OBJECT,
				modelTypeId, spreadVolTypeId, LocalGetNumObjectId(smiledMod), isLnVol, 100, C_result) == ARM_OK)
	{
		bsModGenId = C_result.getLong();
	}
	else
	{
		return resultat;
	}

	if (ARMLOCAL_ARM_Price(clonePentilixId,
						   bsModGenId,
						   C_result) == ARM_OK)
	{
		pxPentilix = C_result.getDouble();
	}
	else
	{
		return resultat;
	}

	vector<double> calibInfos;
	if (ARMLOCAL_SPREADDIGITAL(startDate,
							   endDate,
							   capOrFloorId,
							   1, 
							   LocalGetNumObjectId(strike),								
							   0,
							   payOffSpread,
							   indexId, 
							   indexLongId,
							   weightVector[0], 
							   weightVector[1],
							   dayCountId,
							   resetFreqId,
							   payFreqId, 
							   resetTimingId,
							   payTimingId, 
							   ccyId,
							   resetGap, 
							   spreadVector[0], 
							   spreadVector[1],
							   1,
							   ARM_NULL_OBJECT,
							   C_fixing1_double_vect,
							   1,
							   ARM_NULL_OBJECT,
							   C_fixing2_double_vect,
							   intRuleId,
							   stubRuleId,
							   slopeFlag,
							   cptStrikeMethod, 
							   computedFormula,
							   calibInfos,
							   C_result) == ARM_OK)
	{
		DigitalId = C_result.getLong();
	}
	else
	{
		return resultat;
	}

	if (ARMLOCAL_ClonedAndSetNotional(DigitalId,
									  amortRefvalue,
									  100.,
									  C_result) == ARM_OK)
	{
		cloneDigitalId=C_result.getLong();
	}
	else
	{
		return resultat;
	}

	retCode = ARMLOCAL_FreeObject(DigitalId,C_result);

	if (ARMLOCAL_ARM_Price(cloneDigitalId,
						   bsModGenId,
						   C_result) == ARM_OK)
	{
		pxDigital = C_result.getDouble();
	}
	else
	{
		return resultat;
	}
	
	pxTotal = pxPentilix + pxDigital;

	//  Calcul de la marge Correl
	if (ARMLOCAL_ARM_BumpVolatility(LocalGetNumObjectId(hyperCubeCorrel),1., 0, 0, K_NO, K_YES,C_result) == ARM_OK)
	{
		bumpHyperCubeCorrelId = C_result.getLong();
	}
	else
	{
		return resultat;
	}

	// Les variables indexIndexCorrelCube, corrEUR, InterCorr sont à créer
	// La fonction ARMLOCAL_CreateGenCorrelatorManager est en attente de modifcation par JPR

	vHyperCubeCorrelId.clear();
	vHyperCubeCorrelId.push_back(bumpHyperCubeCorrelId); 
	
	if (ARMLOCAL_CreateGenCorrelatorManager(vCcy, vHyperCubeCorrelId, 
		 vIndexIndexCorrelCubeId, vCorrEURId, vInterCorrId, vIrVolId, vVolVolId, vFXVolHyperCubeId, C_result) ==  ARM_OK)    
    {
		bumpCorrManagerId = C_result.getLong();
	}
	else
	{
		return resultat;
	}

	if (ARMLOCAL_BSMODELGEN(LocalGetNumObjectId(zc), ARM_NULL_OBJECT, LocalGetNumObjectId(vol), LocalGetNumObjectId(volCUB), 
				bumpCorrManagerId, LocalGetNumObjectId(ConvexityManager), 
				LocalGetNumObjectId(zc), bumpHyperCubeCorrelId, ARM_NULL_OBJECT, ARM_NULL_OBJECT,
				modelTypeId, spreadVolTypeId, LocalGetNumObjectId(smiledMod), isLnVol, 100, C_result) == ARM_OK)
	{
		BumpBsModGenId = C_result.getLong();
	}
	else
	{
		return resultat;
	}

	if (ARMLOCAL_ARM_Price(clonePentilixId,
						   BumpBsModGenId,
						   C_result) == ARM_OK)
	{
		pxBumpPentilix = C_result.getDouble();
	}
	else
	{
		return resultat;
	}

	if (ARMLOCAL_ARM_Price(cloneDigitalId,
						   BumpBsModGenId,
						   C_result) == ARM_OK)
	{
		pxBumpDigital = C_result.getDouble();
	}
	else
	{
		return resultat;
	}

	double pxBumpTotal; 
	pxBumpTotal = pxBumpPentilix + pxBumpDigital; 

	retCode = ARMLOCAL_FreeObject(bumpHyperCubeCorrelId,C_result);
	retCode = ARMLOCAL_FreeObject(bumpCorrManagerId,C_result);
	retCode = ARMLOCAL_FreeObject(BumpBsModGenId,C_result);

	if ((pxBumpTotal - pxTotal) > 0 )
	{				
		if (ARMLOCAL_ComputeVolatility(LocalGetNumObjectId(VolCurvFromMatriceShift), d_index, d_indexLong, 0, C_result) == ARM_OK)
		{
			shiftCorrel = C_result.getDouble();
		}
		else
		{
			return resultat;
		}
	}
	else
	{
		if (ARMLOCAL_ComputeVolatility(LocalGetNumObjectId(VolCurvFromMatriceShift),d_indexLong, d_index, 0, C_result) == ARM_OK)
		{
			shiftCorrel = C_result.getDouble();
		}
		else
		{
			return resultat;
		}
	}

	// Bump de l'hyperCube 
	if (ARMLOCAL_ARM_BumpVolatility(LocalGetNumObjectId(hyperCubeCorrel),shiftCorrel, 0, 0, K_NO, K_YES,C_result) == ARM_OK)
	{
		bumpHyperCubeCorrelId = C_result.getLong();
	}
	else
	{
		return resultat;
	}

	// Les variables indexIndexCorrelCube, corrEUR, InterCorr sont à créer
	// La fonction ARMLOCAL_CreateGenCorrelatorManager est en attente de modifcation par JPR

	vHyperCubeCorrelId.clear();
	vHyperCubeCorrelId.push_back(bumpHyperCubeCorrelId); 

	if (ARMLOCAL_CreateGenCorrelatorManager(vCcy, vHyperCubeCorrelId, 
		 vIndexIndexCorrelCubeId, vCorrEURId, vInterCorrId, vIrVolId, vVolVolId, vFXVolHyperCubeId, C_result) ==  ARM_OK)    
	{
		bumpCorrManagerId = C_result.getLong();
	}
	else
	{
		return resultat;
	}

	if (ARMLOCAL_BSMODELGEN(LocalGetNumObjectId(zc), ARM_NULL_OBJECT, LocalGetNumObjectId(vol), LocalGetNumObjectId(volCUB), 
				bumpCorrManagerId, LocalGetNumObjectId(ConvexityManager),
				LocalGetNumObjectId(zc), bumpHyperCubeCorrelId, ARM_NULL_OBJECT, ARM_NULL_OBJECT,
				modelTypeId, spreadVolTypeId, LocalGetNumObjectId(smiledMod), isLnVol, 100, C_result) == ARM_OK)
	{
		BumpBsModGenId = C_result.getLong();
	}
	else
	{
		return resultat;
	}

	if (ARMLOCAL_ARM_Price(clonePentilixId,
						   BumpBsModGenId,
						   C_result) == ARM_OK)
	{
		pxBumpPentilix = C_result.getDouble();
	}
	else
	{
		return resultat;
	}

	if (ARMLOCAL_ARM_Price(cloneDigitalId,
						   BumpBsModGenId,
						   C_result) == ARM_OK)
	{
		pxBumpDigital = C_result.getDouble();
	}
	else
	{
		return resultat;
	}

	pxBumpTotal = pxBumpPentilix + pxBumpDigital;

	retCode = ARMLOCAL_FreeObject(bumpCorrManagerId,C_result);
	retCode = ARMLOCAL_FreeObject(bumpHyperCubeCorrelId,C_result);
	retCode = ARMLOCAL_FreeObject(BumpBsModGenId,C_result);
	retCode = ARMLOCAL_FreeObject(clonePentilixId,C_result);
	retCode = ARMLOCAL_FreeObject(cloneDigitalId,C_result);

	long irIndexIdfixed;
	if (ARMLOCAL_IRINDEX(dayCountId,
						 resetFreqId,
						 -1.0,
						 K_COMP_PROP,
						 K_MOD_FOLLOWING,
						 resetTimingId,
						 resetGap,
						 K_ARREARS,
						 10000.,
						 false,
						 "EUR",
						 K_FIXED,
						 K_COMP_PROP,
						 intRuleId,
						 resetFreqId,
						 C_result) == ARM_OK) 
	{
		irIndexIdfixed = C_result.getLong();
	}
	else
	{
		return resultat;
	}

	long swapLegFixedId;
	if (ARMLOCAL_SWAPLEG(irIndexIdfixed,
						 startDate,
						 endDate,
						 K_RCV,
						 0,
						 payOffSpread,
						 false,
						 "EUR",
						 dayCountId,
						 resetGap,
						 "NULL",
						 "NULL",
						 0, 
						 K_NX_NONE,
						 stubRuleId, 
						 -1.0,
						 K_NO,
						  -1,
						 C_result) == ARM_OK) 
	{
		swapLegFixedId = C_result.getLong();
	}
	else
	{
		return resultat;
	}
	
	// ********** Fixed Leg ************
	long cloneFixedId;
	if (ARMLOCAL_ClonedAndSetNotional(swapLegFixedId,
									  amortRefvalue,
									  100.,
									  C_result) == ARM_OK)
	{
		cloneFixedId = C_result.getLong();
	}
	else
	{
		return resultat;
	}

	retCode = ARMLOCAL_FreeObject(swapLegFixedId,C_result);
	
	//*********   Price Fixed *************************
	double pxFixed;
	if (ARMLOCAL_ARM_Price(cloneFixedId,
						   LocalGetNumObjectId(smiledMod),
						   C_result) == ARM_OK)
	{
		pxFixed = C_result.getDouble();
	}
	else
	{
		return resultat;
	}

	retCode = ARMLOCAL_FreeObject(cloneFixedId,C_result);

	double ProbaExerciceFixed; 
	double MTFixed; 
	ProbaExerciceFixed = pxDigital/pxFixed; 

	//Calcul de la marge technique
	double minMTFixed;
	double maxMTFixed;
	double maxPriceFixed; 

	maxPriceFixed = pxFixed;

	if (ProbaExerciceFixed > (1- ProbaExerciceFixed))
	{
		minMTFixed = 1. - ProbaExerciceFixed;
	}
	else
	{
		minMTFixed = ProbaExerciceFixed;
	}
	
	if (minMTFixed - 0.1 > 0)
	{
		maxMTFixed = minMTFixed - 0.1; 
	}
	else
	{
		maxMTFixed = 0.;
	}

	MTFixed = (0.015 + 0.075 *maxMTFixed)*maxPriceFixed;

// Calcul MT variable

	retCode = ARMLOCAL_FreeObject(corrManagerId,C_result);
	retCode = ARMLOCAL_FreeObject(bsModGenId,C_result);	

	long irIndexPayId; 
	long dayCountPay;
	
	if (IsCMSIndex((ARM_INDEX_TYPE)indexPayId) == 1)
	{
		dayCountPay = K30_360 ; 
	}	
	else
	{
		dayCountPay = KACTUAL_360; 
	}

	if (ARMLOCAL_IRINDEX(dayCountId,
						 payFreqId,
						 -1,
						 K_COMP_PROP,
						 K_MOD_FOLLOWING,
				 		 resetTimingId,	
						 resetGap,
						 K_ARREARS,
						 10000.,
						 false,
						 "EUR",
						 indexPayId,
						 K_COMP_PROP,
						 intRuleId,
						 resetFreqId, 
						 C_result) == ARM_OK) 
	{
		irIndexPayId = C_result.getLong();
	}
	else
	{
		return resultat;
	}

	long swapLegVariableId;
	if (ARMLOCAL_SWAPLEG(irIndexPayId,
						 startDate,
						 endDate,
						 K_RCV,
						 0,
						 0.,
						 false,
						 "EUR",
						 dayCountId, 
						 resetGap,
						 "NULL",
						 "NULL",
						 0,
						 K_NX_NONE,
						 stubRuleId, 
						 -1.0,
						 K_NO,
								 -1,
						 C_result) == ARM_OK) 
	{
		swapLegVariableId = C_result.getLong();
	}
	else
	{
		return resultat;
	}

	// ********** Variable Leg ************
	long cloneVariableId;
	if (ARMLOCAL_ClonedAndSetNotional(swapLegVariableId,
									  amortRefvalue,
									  100.,
									  C_result) == ARM_OK)
	{
		cloneVariableId = C_result.getLong();
	}
	else
	{
		return resultat;
	}

	retCode = ARMLOCAL_FreeObject(swapLegVariableId,C_result);
	
	//*********   Price Variable *************************
	double pxVariable;

	if (ARMLOCAL_ARM_Price(cloneVariableId,
						   LocalGetNumObjectId(smiledMod),
						   C_result) == ARM_OK)
	{
		pxVariable = C_result.getDouble();
	}
	else
	{
		return resultat;
	}

	retCode = ARMLOCAL_FreeObject(cloneVariableId,C_result);

	double ProbaExerciceVariable; 
	ProbaExerciceVariable = pxPentilix/pxVariable;
	
	//Fin de la récup du prix du ssJacent

	//Calcul de la marge technique
	double minMTVariable;
	double maxMTVariable;
	double maxPriceVariable; 
	double MTVariable; 

	maxPriceVariable = pxVariable;

	if (ProbaExerciceVariable > (1- ProbaExerciceVariable))
	{
		minMTVariable = 1. - ProbaExerciceVariable;
	}
	else
	{
		minMTVariable = ProbaExerciceVariable;
	}

	if (minMTVariable - 0.1 > 0)
	{
		maxMTVariable = minMTVariable - 0.1; 
	}
	else
	{
		maxMTVariable = 0.;
	}

	MTVariable = (0.015 + 0.075 *maxMTVariable)*maxPriceVariable;

	double MTGlobal; 
	MTGlobal = fabs(MTVariable + MTFixed); 

	//Calcul de la marge Minimale
	double SpreadMin; 
	SpreadMin = 0.8; 

	long swapLegMinId;
	if (ARMLOCAL_SWAPLEG(irIndexIdfixed,
						 startDate,
						 endDate,
						 K_RCV,
						 0,
						 SpreadMin,
						 false,
						 "EUR",
						 dayCountId,
						 resetGap,
						 "NULL",
						 "NULL",
						 0, 
						 K_NX_NONE,
						 stubRuleId, 
						 -1.0,
						 K_NO,
								 -1,
						 C_result) == ARM_OK) 
	{
		swapLegMinId = C_result.getLong();
	}
	else
	{
		return resultat;
	}
	
	retCode = ARMLOCAL_FreeObject(irIndexIdfixed,C_result);

	// ********** Fixed Leg ************
	long cloneMinId;
	if (ARMLOCAL_ClonedAndSetNotional(swapLegMinId,
									  amortRefvalue,
									  100.,
									  C_result) == ARM_OK)
	{
		cloneMinId = C_result.getLong();
	}
	else
	{
		return resultat;
	}

	retCode = ARMLOCAL_FreeObject(swapLegMinId,C_result);
	
	//*********   Price Fixed *************************
	
	double pxMin;

	if (ARMLOCAL_ARM_Price(cloneMinId,
						   LocalGetNumObjectId(smiledMod),
						   C_result) == ARM_OK)
	{
		pxMin = C_result.getDouble();
	}
	else
	{
		return resultat;
	}

	retCode = ARMLOCAL_FreeObject(cloneMinId,C_result);	
	
	double ProbaExerciceMin; 
	double minMTMin; 
	double maxMTMin; 
	double maxPriceMin; 
	double MTMin; 

	maxPriceMin = fabs(pxMin);

	ProbaExerciceMin = (ProbaExerciceVariable + ProbaExerciceFixed)/2.; 

	if (ProbaExerciceMin > (1- ProbaExerciceMin))
	{
		minMTMin = 1. - ProbaExerciceMin;
	}
	else
	{
		minMTMin = ProbaExerciceMin;
	}
	
	if (minMTMin - 0.1 > 0)
	{
		maxMTMin = minMTMin - 0.1; 
	}
	else
	{
		maxMTMin = 0.;
	}

	MTMin = (0.015 + 0.075 *maxMTMin)*maxPriceMin;

	MTMin = fabs(MTMin); 

	double MTFinal; 

	if (MTMin - MTGlobal> 0)
	{
		MTFinal = MTMin; 
	}
	else
	{
		MTFinal = MTGlobal;
	}

	//************************************************************************************************
	//*************************** Fin création Digitale **********************************************
	//************************************************************************************************

	//*********   Price Première Phase structure 1  TD*************************

	VECTOR<double> vRes;
	
	vRes.push_back(pxTotal);
	vRes.push_back(pxBumpTotal);
	vRes.push_back(MTFinal);

	vHyperCubeCorrelId.clear(); 
	vIndexIndexCorrelCubeId.clear(); 
	vCorrEURId.clear(); 
	vInterCorrId.clear(); 

	retCode = ARMLOCAL_FreeObject(swapLegFixedId,C_result);
	retCode = ARMLOCAL_FreeObject(cloneFixedId,C_result);
 	retCode = ARMLOCAL_FreeObject(swapLegVariableId,C_result);
 	retCode = ARMLOCAL_FreeObject(cloneVariableId,C_result);
 	retCode = ARMLOCAL_FreeObject(swapLegMinId,C_result);
 	retCode = ARMLOCAL_FreeObject(cloneMinId,C_result);
 	retCode = ARMLOCAL_FreeObject(cloneDigitalId,C_result);
 	retCode = ARMLOCAL_FreeObject(corrManagerId,C_result);
 	retCode = ARMLOCAL_FreeObject(bumpCorrManagerId,C_result);
 	retCode = ARMLOCAL_FreeObject(bsModGenId,C_result);
 	retCode = ARMLOCAL_FreeObject(BumpBsModGenId,C_result);
 	retCode = ARMLOCAL_FreeObject(bumpHyperCubeCorrelId,C_result);
 	retCode = ARMLOCAL_FreeObject(DigitalId,C_result);
 	retCode = ARMLOCAL_FreeObject(cloneDigitalId,C_result);
 	retCode = ARMLOCAL_FreeObject(irIndexIdfixed,C_result);

	return vRes;
}
catch(...)
{
	return resultat;
}
}






VECTOR<double> computeReviPentix(const CCString& index, const CCString& indexLong, long irIndexPhase2Id, double startDate, double endDate,
							   const CCString& floorOrCap, const CCString& strike, double SprSeekBonifie, 
							   double Levier, double TxFixeMax, 
							   long ccyId, long dayCountId, 
							   long resetTimingId, long payTimingId, double resetGap, double resetGap2, long resetFreqId, long payFreqId,
							   long intRuleId, long stubRuleId, long amortRefvalue,
							   VECTOR<double> vDeltaOther, 
							   long iscappedId, double TxCap,
							   const VECTOR<CCString>& vBumpBsGenMod, const VECTOR<CCString>& vBumpVolBsGenMod, 
							   const CCString& hyperCubeCorrel, const CCString& VolCurvFromMatriceShift,
							   const CCString& vol, const CCString& volCUB, const CCString& CorrManager,
							   const CCString& ConvexityManager, const CCString& zc, const CCString& smiledMod,
							   VECTOR<double>* vegaSpreadOption, VECTOR<double>* deltaSpreadOption, VECTOR<double>* vegaSpreadOptionCap, VECTOR<double>* deltaSpreadOptionCap)
{
	VECTOR<double> resultat;
	resultat.clear();

try
{
	double pxDigital;
	double pxSpreadOption; 
	double pxBumpDigital; 
	double pxBumpSpreadOption; 
	double pxBumpSpreadOptionCap; 
	double MT;
	double MargeCorrel;
	double ProbaExercice; 
	long digitalId;
	long cloneDigitalId; 
	long spreadOptionId; 
	long cloneSpreadOptionId;
	long spreadOptionCapId; 
	long cloneSpreadOptionCapId;
	double pxSpreadOptionCap; 

	long res;
	

	double shiftCorrel; 

	long bumpHyperCubeCorrelId;
	ARM_result C_result;

	// Modèles
	long bsModGenId; 
	long BumpBsModGenId;

	VECTOR <double> spreadVector;
	VECTOR <double> weightVector;

   	VECTOR<double> C_fixing1_double_vect;
   	VECTOR<double> C_fixing2_double_vect;

	long slopeFlag; 
	int cptStrikeMethod;
	int computedFormula;

	slopeFlag = 1; 
	cptStrikeMethod = 1.; 
	computedFormula = 1.; 

	// Ss-Jacent
	long ssJacentId;
	long cloneSsJacentId;
	double pxSsJacent;

	VECTOR <double> deltaDigital; 

	///// Calcul de la patte fixe de la structure//////////////////
	long swapLegPhase2Id; 
	long clonePhase2Id; 

	if (ARMLOCAL_SWAPLEG(irIndexPhase2Id,
						 startDate,
						 endDate,
						 K_RCV,
						 0, 
						 SprSeekBonifie, 
						 false,
						 "EUR",
						 dayCountId,
						 resetGap2,
						 "NULL",
						 "NULL",
						 0,
						 K_NX_NONE,
						 stubRuleId,
						 -1.0,
						 K_NO,
								 -1,
						 C_result) == ARM_OK) 
	{
		swapLegPhase2Id = C_result.getLong();
	}

	// ********** Phase2 Leg ************

	if (ARMLOCAL_ClonedAndSetNotional(swapLegPhase2Id,
									  amortRefvalue,
									  100.,
									  C_result) == ARM_OK)
	{
		clonePhase2Id = C_result.getLong();
	}

	//*********   Price Phase2 *************************
	double pxPhase2;

	VECTOR<double> deltaPhase2NOTDIG; 

	if (ARMLOCAL_ARM_Price(clonePhase2Id,
						   LocalGetNumObjectId(smiledMod),
						   C_result) == ARM_OK)
	{
		pxPhase2 = C_result.getDouble();
	}

	#ifdef _DEBUG
		res = ARMLOCAL_ARM_View(irIndexPhase2Id,C_result);
	#endif

	#ifdef _DEBUG
		res = ARMLOCAL_ARM_View(swapLegPhase2Id,C_result);
	#endif


	deltaPhase2NOTDIG = computeDeltaByPlot(pxPhase2, clonePhase2Id, vBumpBsGenMod);

	long retCode = ARMLOCAL_FreeObject(swapLegPhase2Id,C_result);
	retCode = ARMLOCAL_FreeObject(clonePhase2Id,C_result);

	//////////////////////////////////////////////////////////////////////////////////

	double d_index; 
	long indexId; 
	indexId = ARM_ConvIrType(index); 

	if (IsCMSIndex((ARM_INDEX_TYPE)indexId) == 1)
	{
		d_index = (indexId - K_CMS1 +1 ); 
	}	
	else
	{
		d_index = 1./FromLiborTypeToFrequency(indexId); 
	}
	
	double d_indexLong; 
	long indexLongId; 
	indexLongId = ARM_ConvIrType(indexLong); 

	if (IsCMSIndex((ARM_INDEX_TYPE)indexLongId) ==1)
	{
		d_indexLong = indexLongId - K_CMS1 + 1; 
	}
	else
	{
		d_indexLong = 1./FromLiborTypeToFrequency(indexLongId); 
	}

	//Définition du sreadVector et du WeightVector
	spreadVector.push_back(-0.01);
	spreadVector.push_back(0.01); 

	weightVector.push_back(1);
	weightVector.push_back(1);

	long capOrFloorId; 

	if((capOrFloorId = ARM_ConvCapOrFloor(floorOrCap, C_result)) == ARM_DEFAULT_ERR)
	{
		return resultat;
	}

	long modelTypeId;
	long spreadVolTypeId;
	CCString modelType; 
	CCString spreadVolType;
	bool isLnVol=true;

	modelType = "2LOG";
	spreadVolType = "C"; 

	if ((modelTypeId = ARM_ConvModelType (modelType, C_result)) == ARM_DEFAULT_ERR)
	{
		return resultat;
	}

	if ((spreadVolTypeId = ARM_ConvVolType2 (spreadVolType, C_result)) == ARM_DEFAULT_ERR)
	{
		return resultat;
	}

	VECTOR<double> C_rate; 
	VECTOR<CCString> C_matu;
	double C_meth = 0;	

	long txTrans;  
	if (ARMLOCAL_CONSTREFVALUE(TxFixeMax - SprSeekBonifie,
						  C_result) == ARM_OK)
	{
		txTrans = C_result.getLong();
	}

	long payOff; 
	if (ARMLOCAL_SumRefValue ( txTrans,
								  LocalGetNumObjectId(strike),
								  -1.*Levier,
								  C_result)== ARM_OK)
	{
		payOff = C_result.getLong();
	}

	//************************************************************************************************
	//*************************** Création Digitale **************************************************
	//************************************************************************************************

	//Calcul du prix de la Digitale
	vector<double> calibInfos;

	if (ARMLOCAL_SPREADDIGITAL(startDate,
						 endDate,
						 capOrFloorId,
						 1, 
						 LocalGetNumObjectId(strike), 
						 1, 
						 payOff,
						 indexId,
						 indexLongId,
						 weightVector[0], 
						 weightVector[1],
						 dayCountId,
						 resetFreqId, 
						 payFreqId, 
						 resetTimingId, 
						 payTimingId, 					
						 ccyId,
						 resetGap, 
						 spreadVector[0], 
						 spreadVector[1],
						 1,
						 ARM_NULL_OBJECT,
						 C_fixing1_double_vect,
						 1,
						 ARM_NULL_OBJECT,
						 C_fixing2_double_vect,
						 intRuleId, 
						 stubRuleId, 
						 slopeFlag,
						 cptStrikeMethod, 
						 computedFormula,
						 calibInfos,
						 C_result) == ARM_OK)
	{
		digitalId = C_result.getLong();
	}
	else
	{
		return resultat;
	}

	if (ARMLOCAL_ClonedAndSetNotional(digitalId,
									  amortRefvalue,
									  100.,
									  C_result) == ARM_OK)
	{
		cloneDigitalId = C_result.getLong();
	}
	else
	{
		return resultat;
	}
	
	//************************************************************************************************
	//*************************** Fin Création Digitale **********************************************
	//************************************************************************************************

	if (ARMLOCAL_BSMODELGEN(LocalGetNumObjectId(zc), ARM_NULL_OBJECT, LocalGetNumObjectId(vol), LocalGetNumObjectId(volCUB), 
				LocalGetNumObjectId(CorrManager), LocalGetNumObjectId(ConvexityManager), 
				LocalGetNumObjectId(zc), LocalGetNumObjectId(hyperCubeCorrel), ARM_NULL_OBJECT, ARM_NULL_OBJECT,
                modelTypeId, spreadVolTypeId, LocalGetNumObjectId(smiledMod), isLnVol, 100, C_result) == ARM_OK)
	{
		bsModGenId = C_result.getLong();
	}
	else
	{
		return resultat;
	}

	if (ARMLOCAL_ARM_Price(cloneDigitalId,
						   bsModGenId,
						   C_result) == ARM_OK)
	{
		pxDigital = C_result.getDouble();
	}
	else
	{
		return resultat;
	}

	#ifdef _DEBUG
		res = ARMLOCAL_ARM_View(digitalId,C_result);
	#endif

		
	//************************************************************************************************
	//************************************ Calcul Sensi **********************************************
	//************************************************************************************************
	
	deltaDigital = computeDeltaByPlot(pxDigital, cloneDigitalId, vBumpBsGenMod);

	//************************************************************************************************
	//********************************** Fin Calcul Sensi ********************************************
	//************************************************************************************************
	//  Calcul de la marge Correl

	//************************************************************************************************
	//*************************** Création SpreadOption **********************************************
	//************************************************************************************************
	if (ARMLOCAL_SPREADOPTION (startDate,
							   endDate,
							   capOrFloorId,
							   1,
							   LocalGetNumObjectId(strike),
							   indexId,
							   indexLongId,
							   1.,
							   false,
							   weightVector[0],
							   1.,
							   false,
							   weightVector[1],
							   dayCountId,
							   resetFreqId,
							   payFreqId,
							   resetTimingId,
							   payTimingId,
							   ccyId,
							   resetGap,
							   1,
							   ARM_NULL_OBJECT,
							   C_fixing1_double_vect,
							   1,
							   ARM_NULL_OBJECT,
							   C_fixing2_double_vect,
							   intRuleId,
							   stubRuleId,
							   cptStrikeMethod, 
							   computedFormula,
							   calibInfos,
							   C_result) == ARM_OK)
	{
		spreadOptionId = C_result.getLong();
	}
	else
	{
		return resultat;
	}

	if (ARMLOCAL_ClonedAndSetNotional(spreadOptionId,
									  amortRefvalue,
									  100.,
									  C_result) == ARM_OK)
	{
		cloneSpreadOptionId = C_result.getLong();
	}
	else
	{
		return resultat;
	}

	if (ARMLOCAL_ARM_Price(cloneSpreadOptionId,
						   bsModGenId, 
						   C_result) == ARM_OK)
	{
		pxSpreadOption = C_result.getDouble();
	}
	else
	{
		return resultat;
	}

	#ifdef _DEBUG
		res = ARMLOCAL_ARM_View(spreadOptionId,C_result);
	#endif

	//************************************************************************************************
	//*************************** Fin Création SpreadOption ******************************************
	//************************************************************************************************

	

	if ((*vegaSpreadOption).size() == 0)
	{
		*vegaSpreadOption = computeDeltaByPlot(pxSpreadOption, cloneSpreadOptionId, vBumpVolBsGenMod);
		*deltaSpreadOption = computeDeltaByPlot(pxSpreadOption, cloneSpreadOptionId, vBumpBsGenMod);
	}

	

	//************************************************************************************************
	//*************************** Création SpreadOption Cap*******************************************
	//************************************************************************************************

	double cap;

	if (iscappedId==K_YES)
	{
		cap=(TxFixeMax-TxCap)/Levier;

		if (ARMLOCAL_SPREADOPTION (startDate,
								   endDate,
								   capOrFloorId,
								   0,
								   cap,
								   indexId,
								   indexLongId,
								   1.,
								   false,
								   weightVector[0],
								   1.,
								   false,
								   weightVector[1],
								   dayCountId,
								   resetFreqId,
								   payFreqId,
								   resetTimingId,
								   payTimingId,
								   ccyId,
								   resetGap,
								   1,
								   ARM_NULL_OBJECT,
								   C_fixing1_double_vect,
								   1,
								   ARM_NULL_OBJECT,
								   C_fixing2_double_vect,
								   intRuleId,
								   stubRuleId,
								   cptStrikeMethod, 
								   computedFormula,
								   calibInfos,
								   C_result) == ARM_OK)
		{
			spreadOptionCapId = C_result.getLong();
		}
		else
		{
			return resultat;
		}

		if (ARMLOCAL_ClonedAndSetNotional(spreadOptionCapId,
										  amortRefvalue,
										  100.,
										  C_result) == ARM_OK)
		{
			cloneSpreadOptionCapId = C_result.getLong();
		}
		else
		{
			return resultat;
		}

		if (ARMLOCAL_ARM_Price(cloneSpreadOptionCapId,
							   bsModGenId, 
							   C_result) == ARM_OK)
		{
			pxSpreadOptionCap = C_result.getDouble();
		}
		else
		{
			return resultat;
		}

		#ifdef _DEBUG
			res = ARMLOCAL_ARM_View(spreadOptionCapId,C_result);
		#endif

		if ((*vegaSpreadOptionCap).size() == 0)
		{
			*vegaSpreadOptionCap = computeDeltaByPlot(pxSpreadOptionCap, cloneSpreadOptionCapId, vBumpVolBsGenMod);
			*deltaSpreadOptionCap = computeDeltaByPlot(pxSpreadOptionCap, cloneSpreadOptionCapId, vBumpBsGenMod);
		}
	}
	else
	{
		pxSpreadOptionCap=0;
	}

	//************************************************************************************************
	//*************************** Fin Création SpreadOption Cap***************************************
	//************************************************************************************************
	
	VECTOR <double> vegaDigital; 
	VECTOR <double> vegaGlobal; 
	double vega = 0.0;

	long sizeVega = vBumpVolBsGenMod.size(); 
	vegaDigital= computeDeltaByPlot(pxDigital, cloneDigitalId, vBumpVolBsGenMod);


	if (iscappedId==K_YES)
	{
		for (int i = 0; i < sizeVega; i++)	
		{
			vegaGlobal.push_back(vegaDigital[i] + Levier * (*vegaSpreadOption)[i] - Levier * (*vegaSpreadOptionCap)[i]);
			vega += fabs(vegaGlobal[i]);
		}
	}
	else
	{
		for (int i = 0; i < sizeVega; i++)	
		{
			vegaGlobal.push_back(vegaDigital[i] + Levier * (*vegaSpreadOption)[i]);
			vega += fabs(vegaGlobal[i]);
		}
	}

	//************************************************************************************************
	//***********************************  Calcul Delta   ********************************************
	//************************************************************************************************
	
	long sizeMatu = vDeltaOther.size(); 
	double deltaPositif= 0. ; 
	double  deltaNegatif = 0.;

	
	if (iscappedId==K_YES)
	{
		for (int i = 0; i < sizeMatu; i++)	
		{
			if ((vDeltaOther[i] + (Levier * ( (*deltaSpreadOption)[i] - (*deltaSpreadOptionCap)[i] )) + (deltaDigital[i]) + deltaPhase2NOTDIG[i]) > 0.)
			/*if ((vDeltaOther[i] + (Levier * (*deltaSpreadOption)[i]) + (deltaDigital[i]) + deltaPhase2NOTDIG[i] - Levier * (*deltaSpreadOptionCap)[i])) > 0.)*/
			{
				deltaPositif = deltaPositif + (vDeltaOther[i] + Levier * (*deltaSpreadOption)[i] + deltaDigital[i]+ deltaPhase2NOTDIG[i] - Levier * (*deltaSpreadOptionCap)[i]);
			}
			else 
			{
				deltaNegatif = deltaNegatif + abs(vDeltaOther[i] + Levier * (*deltaSpreadOption)[i] + deltaDigital[i]+ deltaPhase2NOTDIG[i] - Levier * (*deltaSpreadOptionCap)[i]);
			}
		}
	}
	else
	{
		for (int i = 0; i < sizeMatu; i++)	
		{
			if ((vDeltaOther[i] + (Levier * (*deltaSpreadOption)[i]) + (deltaDigital[i]) + deltaPhase2NOTDIG[i]) > 0.)
			{
				deltaPositif = deltaPositif + (vDeltaOther[i] + Levier * (*deltaSpreadOption)[i] + deltaDigital[i]+ deltaPhase2NOTDIG[i]);
			}
			else 
			{
				deltaNegatif = deltaNegatif + abs(vDeltaOther[i] + Levier * (*deltaSpreadOption)[i] + deltaDigital[i]+ deltaPhase2NOTDIG[i]);
			}
		}
	}

	double deltaGlobal; 
	if (deltaPositif > deltaNegatif)
	{
		deltaGlobal = deltaPositif; 
	}
	else
	{
		deltaGlobal = deltaNegatif; 
	}

	//************************************************************************************************
	//********************************* Fin Calcul Delta *********************************************
	//************************************************************************************************


	//************************************************************************************************
	//********************************* Fin Calcul Vega *********************************************
	//************************************************************************************************

	if (ARMLOCAL_ARM_BumpVolatility(LocalGetNumObjectId(hyperCubeCorrel),1., 0, 0, K_NO, K_YES,C_result) == ARM_OK)
	{
		bumpHyperCubeCorrelId = C_result.getLong();
	}
	else
	{
		return resultat;
	}

	if (ARMLOCAL_BSMODELGEN(LocalGetNumObjectId(zc), ARM_NULL_OBJECT, LocalGetNumObjectId(vol), LocalGetNumObjectId(volCUB), 
				LocalGetNumObjectId(CorrManager), LocalGetNumObjectId(ConvexityManager), 
				LocalGetNumObjectId(zc), bumpHyperCubeCorrelId, ARM_NULL_OBJECT, ARM_NULL_OBJECT,
                modelTypeId, spreadVolTypeId, LocalGetNumObjectId(smiledMod), isLnVol, 100, C_result) == ARM_OK)
	{
		BumpBsModGenId = C_result.getLong();
	}
	else
	{
		return resultat;
	}

	if (ARMLOCAL_ARM_Price(cloneDigitalId,
						   BumpBsModGenId,
						   C_result) == ARM_OK)
	{
		pxBumpDigital = C_result.getDouble();
	}
	else
	{
		return resultat;
	}

	if (ARMLOCAL_ARM_Price(cloneSpreadOptionId,
						   BumpBsModGenId,
						   C_result) == ARM_OK)
	{
		pxBumpSpreadOption = C_result.getDouble();
	}
	else
	{
		return resultat;
	}

	if (iscappedId==K_YES)
	{
		if (ARMLOCAL_ARM_Price(cloneSpreadOptionCapId,
							   BumpBsModGenId,
							   C_result) == ARM_OK)
		{
			pxBumpSpreadOptionCap = C_result.getDouble();
		}
		else
		{
			return resultat;
		}
	}
	else
	{
		pxBumpSpreadOptionCap=0;
	}

	retCode = ARMLOCAL_FreeObject(bumpHyperCubeCorrelId,C_result);
	retCode = ARMLOCAL_FreeObject(BumpBsModGenId,C_result);

	//Modif BP pour calcul correl*******************************************************************************
	if (((pxBumpDigital - pxDigital) +  Levier * (pxBumpSpreadOption - pxSpreadOption) -  Levier * (pxBumpSpreadOptionCap - pxSpreadOptionCap)) > 0 )
	{
		if (ARMLOCAL_ComputeVolatility(LocalGetNumObjectId(VolCurvFromMatriceShift), d_index, d_indexLong, 0, C_result) == ARM_OK)
		{
			shiftCorrel = C_result.getDouble();
		}
		else
		{
			return resultat;
		}
	}
	else
	{
		if (ARMLOCAL_ComputeVolatility(LocalGetNumObjectId(VolCurvFromMatriceShift),d_indexLong, d_index, 0, C_result) == ARM_OK)
		{
			shiftCorrel = C_result.getDouble();
		}
		else
		{
			return resultat;
		}
	}
	
	// Bump de l'hyperCube 
	if (ARMLOCAL_ARM_BumpVolatility(LocalGetNumObjectId(hyperCubeCorrel),shiftCorrel, 0, 0, K_NO, K_YES,C_result) == ARM_OK)
	{
		bumpHyperCubeCorrelId = C_result.getLong();
	}
	else
	{
		return resultat;
	}
	
	if (ARMLOCAL_BSMODELGEN(LocalGetNumObjectId(zc), ARM_NULL_OBJECT, LocalGetNumObjectId(vol), LocalGetNumObjectId(volCUB), 
				LocalGetNumObjectId(CorrManager), LocalGetNumObjectId(ConvexityManager), 
				LocalGetNumObjectId(zc), bumpHyperCubeCorrelId, ARM_NULL_OBJECT, ARM_NULL_OBJECT,
				modelTypeId, spreadVolTypeId, LocalGetNumObjectId(smiledMod), isLnVol, 100, C_result) == ARM_OK)
	{
		BumpBsModGenId = C_result.getLong();
	}
	else
	{
		return resultat;
	}


	if (ARMLOCAL_ARM_Price(cloneDigitalId,
							   BumpBsModGenId,
							   C_result) == ARM_OK)
	{
		pxBumpDigital = C_result.getDouble();
	}
	else
	{
		return resultat;
	}
	
	if (ARMLOCAL_ARM_Price(cloneSpreadOptionId,
						   BumpBsModGenId,
						   C_result) == ARM_OK)
	{
		pxBumpSpreadOption = C_result.getDouble();
	}
	else
	{
		return resultat;
	}

	if (iscappedId==K_YES)
	{
		if (ARMLOCAL_ARM_Price(cloneSpreadOptionCapId,
						   BumpBsModGenId,
						   C_result) == ARM_OK)
		{
			pxBumpSpreadOptionCap = C_result.getDouble();
		}
		else
		{
			return resultat;
		}
	}
	else
	{
		pxBumpSpreadOptionCap=0;
	}

	// Marge correl correspond à un maxsumsign

	double margepositive;
	double margenegative;
	margepositive=0;
	margenegative=0;

	if ( (pxBumpDigital-pxDigital) >=0 )
	{
		margepositive+=pxBumpDigital-pxDigital;
	}
	else
	{
		margenegative+=fabs(pxBumpDigital-pxDigital);
	}

	if ( Levier*(pxBumpSpreadOption-pxSpreadOption) >=0 )
	{
		margepositive+=Levier*(pxBumpSpreadOption-pxSpreadOption);
	}
	else
	{
		margenegative+=fabs(Levier*(pxBumpSpreadOption-pxSpreadOption));
	}

	if ( -Levier*(pxBumpSpreadOptionCap-pxSpreadOptionCap) >=0 )
	{
		margepositive+=-Levier*(pxBumpSpreadOptionCap-pxSpreadOptionCap);
	}
	else
	{
		margenegative+=fabs(-Levier*(pxBumpSpreadOptionCap-pxSpreadOptionCap));
	}

	MargeCorrel=max(margepositive,margenegative);


	double strikeSsJacent; 
	//calcul de la discontinuité
	if (floorOrCap == "C" )
	{
		strikeSsJacent = -1000.;
	}
	else 
	{
		strikeSsJacent = 1000.;
	}

	if (ARMLOCAL_SPREADDIGITAL(startDate,
							   endDate,
							   capOrFloorId,
							   0, 
							   strikeSsJacent, 
							   1, 
							   payOff,
							   indexId, 
							   indexLongId,
							   weightVector[0], 
							   weightVector[1],
							   dayCountId,
							   resetFreqId,
							   payFreqId, 
							   resetTimingId,
							   payTimingId, 
							   ccyId,
							   resetGap, 
							   spreadVector[0], 
							   spreadVector[1],
							   1,
							   ARM_NULL_OBJECT,
							   C_fixing1_double_vect,
							   1,
							   ARM_NULL_OBJECT,
							   C_fixing2_double_vect,
							   intRuleId,
							   stubRuleId,
							   slopeFlag,
							   cptStrikeMethod, 
							   computedFormula,
							   calibInfos,
							   C_result) == ARM_OK)
	{
		ssJacentId = C_result.getLong();
	}
	else
	{
		return resultat;
	}

	if (ARMLOCAL_ClonedAndSetNotional(ssJacentId,
									  amortRefvalue,
									  100.,
									  C_result) == ARM_OK)
	{
		cloneSsJacentId = C_result.getLong();
	}
	else
	{
		return resultat;
	}

	if (ARMLOCAL_ARM_Price(cloneSsJacentId,
						   bsModGenId,
						   C_result) == ARM_OK)
	{
		pxSsJacent = C_result.getDouble();
	}
	else
	{
		return resultat;
	}

	//Fin de la récup du prix du ssJacent

	//Calcul de la marge technique
	double minMT;
	double maxMT;
	
	ProbaExercice = pxDigital/pxSsJacent;

	if (ProbaExercice > (1- ProbaExercice))
	{
		minMT = 1. - ProbaExercice; 
	}
	else
	{
		minMT = ProbaExercice;
	}

	if (minMT - 0.1 > 0)
	{
		maxMT = minMT - 0.1; 
	}
	else
	{
		maxMT = 0.;
	}

	MT = fabs((0.015 + 0.075 *maxMT)*pxSsJacent);

	//************************************************************************************************
	//*************************** Fin création Digitale **********************************************
	//************************************************************************************************


	//*********   Price Première Phase structure 1  TD*************************
	VECTOR<double> vRes;
	
	
	vRes.push_back(pxDigital+ Levier * pxSpreadOption + pxPhase2 - Levier * pxSpreadOptionCap);
	vRes.push_back(vega);
	vRes.push_back(deltaGlobal);
	vRes.push_back(MargeCorrel);
	vRes.push_back(MT);

	retCode = ARMLOCAL_FreeObject(digitalId,C_result);
	retCode = ARMLOCAL_FreeObject(cloneDigitalId,C_result);
	retCode = ARMLOCAL_FreeObject(spreadOptionId,C_result);
	retCode = ARMLOCAL_FreeObject(cloneSpreadOptionId,C_result);
	retCode = ARMLOCAL_FreeObject(txTrans,C_result);
	retCode = ARMLOCAL_FreeObject(payOff,C_result);
	retCode = ARMLOCAL_FreeObject(bsModGenId,C_result);
	retCode = ARMLOCAL_FreeObject(BumpBsModGenId,C_result);
	retCode = ARMLOCAL_FreeObject(bumpHyperCubeCorrelId,C_result);
	retCode = ARMLOCAL_FreeObject(ssJacentId,C_result);
	retCode = ARMLOCAL_FreeObject(cloneSsJacentId,C_result);
	if (iscappedId==K_YES)
	{
		retCode = ARMLOCAL_FreeObject(spreadOptionCapId,C_result);
		retCode = ARMLOCAL_FreeObject(cloneSpreadOptionCapId,C_result);
	}
	
	return vRes;
}
catch(...)
{
	return resultat;
}
}



STDMETHODIMP ActiveXModule::ARMcomputeReviPentix(double pNtl, 
										   VARIANT pDate,  
										   BSTR pCcy,
										   BSTR pIndexPhase1, 										   
										   VARIANT pSpread,
										   BSTR pDayCountPhase1, 
										   BSTR pPayFreqPhase1,
										   BSTR pResetFreqPhase1,
										   BSTR pResetTimingPhase1,
										   BSTR pRoll, 
										   BSTR pAdjPhase1, 
										   BSTR pStub,
										   BSTR pIndexPhase2DIG,
										   BSTR pIndexLongPhase2DIG,
										   BSTR pStrikePhase2DIG,
										   BSTR pResetTimingPhase2DIG,
										   BSTR pAdjPhase2DIG,										   
										   BSTR pDayCountPhase2, 
										   BSTR pPayFreqPhase2, 
										   BSTR pResetFreqPhase2,
										   BSTR pAdjPhase2,
										   BSTR pIndexPhase3, 
										   BSTR pDayCountPhase3, 
										   BSTR pPayFreqPhase3, 
										   BSTR pResetFreqPhase3,
										   BSTR pResetTimingPhase3,
										   BSTR pAdjPhase3,
										   BSTR pIndexFund,
										   VARIANT pSpreadFund,
										   BSTR pDayCountFund, 
										   BSTR pPayFreqFund, 
										   BSTR pResetFreqFund,
										   BSTR pResetTimingFund,
										   BSTR pAdjFund, 
										   BSTR pDayCountAmort,
										   BSTR pIntRuleAmort, 
										   double pTxAmort,
										   BSTR pFreqAmort,
										   double pAmountAmort, 
										   BSTR pTypeAmort, 
										   BSTR pFloorOrCap, 
										   double pFee, 
										   double pLevier, 
										   double pTxFixeMax, 
										   BSTR pIsCapped,
										   double pTxCap,
										   BSTR pVolCurvFromMatriceShift,
										   BSTR pVol, 
										   BSTR pVolCub, 
										   BSTR pCorrManager, 
										   BSTR pConvexityManager, 
										   BSTR pZc, 
										   BSTR pSmiledMod, 					
										   BSTR pHyperCubeCorrel,
										   VARIANT *pBumpBsGenMod,
										   VARIANT *pBumpVolBsGenMod, 
										   VARIANT *pRet
										   )

{	
	// TODO: Add your implementation code here
try
{
	double Res = 0.;
	ARM_result C_result;

	double Levier = pLevier; 
	double TxFixeMax = pTxFixeMax;
	double TxCap = pTxCap;

	double spreadFund;
	long  spreadType;
	long sizeDate; 
	long sizeSpread; 

	VECTOR<double> deltaPhase1; 
	VECTOR<double> deltaPhase3; 
	VECTOR<double> deltaPhase2NOTDIG;
	VECTOR<double> deltaOther; 
	VECTOR<double> deltaFund; 
	VECTOR<double> vDate; 
	VECTOR<double> vSpread; 

	VECTOR<CCString> vBumpBsGenMod; 
	VECTOR<CCString>  vBumpVolBsGenMod; 
	
	double startDatePhase1;  
	double startDatePhase2; 
	double startDatePhase3;
	double endDatePhase3;
	double spreadPhase1;
	double spreadPhase3;

	if(VARIANT2VECTORDOUBLE (pDate,vDate,sizeDate) != S_OK)
	{
		ERROR_MSG("Error conversion Dates",pRet,ARM_ERROR_FREQ);
		return S_OK;
	}

	if(VARIANT2VECTORDOUBLE (pSpread,vSpread,sizeSpread) != S_OK)
	{
		ERROR_MSG("Error conversion Spread",pRet,ARM_ERROR_FREQ);
		return S_OK;
	}
	
	if (sizeDate !=0 )
	{
		startDatePhase1 = vDate[0];  
		startDatePhase2 = vDate[1];
		startDatePhase3 = vDate[2];
		endDatePhase3 = vDate[3];	
	}

	if (sizeSpread!=0 )
	{
		spreadPhase1 = vSpread[0];
		spreadPhase3 = vSpread[1];
	}

	if (VARIANT2Double(pSpreadFund,&spreadFund) == S_FALSE)
	{
		_bstr_t bSpreadFund(pSpreadFund);
		CCString l_SpreadFund= bSpreadFund;
		spreadFund = (double) LocalGetNumObjectId(l_SpreadFund);
		spreadType = 1;
	}
	else
	   spreadType = 0;

	long sizeMatu;
	
	if(VARIANT2VECTORCCSTRING (*pBumpBsGenMod,vBumpBsGenMod,sizeMatu) != S_OK)
	{
		ERROR_MSG("Error conversion Model Delta",pRet,ARM_ERROR_FREQ);
		return S_OK;
	}

	if(VARIANT2VECTORCCSTRING (*pBumpVolBsGenMod,vBumpVolBsGenMod,sizeMatu) != S_OK)
	{
		ERROR_MSG("Error conversion Model Vega",pRet,ARM_ERROR_FREQ);
		return S_OK;
	}

	double Ntl = pNtl;
	double fee = pFee; 

	long retCode;

	_bstr_t ccy(pCcy);
	CCString l_ccy = ccy;

	_bstr_t indexPhase1(pIndexPhase1);
	CCString l_indexPhase1 = indexPhase1;

	_bstr_t dayCountPhase1(pDayCountPhase1);
	CCString l_dayCountPhase1 = dayCountPhase1;

	_bstr_t resetFreqPhase1(pResetFreqPhase1);
	CCString l_resetFreqPhase1 = resetFreqPhase1;

	_bstr_t payFreqPhase1(pPayFreqPhase1);
	CCString l_payFreqPhase1 = payFreqPhase1;

	_bstr_t resetTimingPhase1(pResetTimingPhase1);
	CCString l_resetTimingPhase1 = resetTimingPhase1;

	_bstr_t roll(pRoll);
	CCString l_roll = roll;

	_bstr_t adjPhase1(pAdjPhase1);
	CCString l_adjPhase1 = adjPhase1;
	
	_bstr_t stub(pStub);
	CCString l_stub= stub;

	_bstr_t indexPhase2DIG(pIndexPhase2DIG);
	CCString l_indexPhase2DIG = indexPhase2DIG;

	_bstr_t indexLongPhase2DIG(pIndexLongPhase2DIG);
	CCString l_indexLongPhase2DIG = indexLongPhase2DIG;

	//Strike is double??
	_bstr_t strikePhase2DIG(pStrikePhase2DIG);
	CCString l_strikePhase2DIG= strikePhase2DIG;

	_bstr_t resetTimingPhase2DIG(pResetTimingPhase2DIG);
	CCString l_resetTimingPhase2DIG = resetTimingPhase2DIG;

	_bstr_t adjPhase2DIG(pAdjPhase2DIG);
	CCString l_adjPhase2DIG = adjPhase2DIG;

	_bstr_t dayCountPhase2(pDayCountPhase2);
	CCString l_dayCountPhase2= dayCountPhase2;

	_bstr_t payFreqPhase2(pPayFreqPhase2);
	CCString l_payFreqPhase2= payFreqPhase2;
	
	_bstr_t resetFreqPhase2(pResetFreqPhase2);
	CCString l_resetFreqPhase2 = resetFreqPhase2;

	_bstr_t adjPhase2(pAdjPhase2);
	CCString l_adjPhase2 = adjPhase2;

	_bstr_t indexPhase3(pIndexPhase3);
	CCString l_indexPhase3 = indexPhase3;

	_bstr_t dayCountPhase3(pDayCountPhase3);
	CCString l_dayCountPhase3 = dayCountPhase3;

	_bstr_t payFreqPhase3(pPayFreqPhase3);
	CCString l_payFreqPhase3 = payFreqPhase3;

	_bstr_t resetFreqPhase3(pResetFreqPhase3);
	CCString l_resetFreqPhase3 = resetFreqPhase3;

	_bstr_t resetTimingPhase3(pResetTimingPhase3);
	CCString l_resetTimingPhase3 = resetTimingPhase3;

	_bstr_t adjPhase3(pAdjPhase3);
	CCString l_adjPhase3 = adjPhase3;

	_bstr_t indexFund(pIndexFund);
	CCString l_indexFund = indexFund;

	_bstr_t dayCountFund(pDayCountFund);
	CCString l_dayCountFund= dayCountFund;

	_bstr_t payFreqFund(pPayFreqFund);
	CCString l_payFreqFund= payFreqFund;
	
	_bstr_t resetFreqFund(pResetFreqFund);
	CCString l_resetFreqFund = resetFreqFund;

	_bstr_t resetTimingFund(pResetTimingFund);
	CCString l_resetTimingFund = resetTimingFund;

	_bstr_t adjFund(pAdjFund);
	CCString l_adjFund= adjFund;
	
	_bstr_t typeAmort(pTypeAmort);
	CCString l_typeAmort = typeAmort;

	_bstr_t floorOrCap (pFloorOrCap);
	CCString l_floorOrCap = floorOrCap; 

	_bstr_t volCurvFromMatriceShift(pVolCurvFromMatriceShift);
	CCString l_volCurvFromMatriceShift= volCurvFromMatriceShift;
	
	_bstr_t vol(pVol);
	CCString l_vol = vol;

	_bstr_t volCub(pVolCub);
	CCString l_volCub = volCub;

	_bstr_t corrManager(pCorrManager);
	CCString l_corrManager = corrManager;

	_bstr_t convexityManager(pConvexityManager);
	CCString l_convexityManager = convexityManager;

	_bstr_t zc(pZc);
	CCString l_zc = zc;

	_bstr_t smiledMod (pSmiledMod);
	CCString l_smiledMod = smiledMod;


	_bstr_t hyperCubeCorrel (pHyperCubeCorrel);
	CCString l_hyperCubeCorrel = hyperCubeCorrel; 

	long payFreqPhase1Id;	
	if((payFreqPhase1Id = ARM_ConvFrequency(l_payFreqPhase1, C_result)) == ARM_DEFAULT_ERR)
	{
		ERROR_MSG("Invalid Pay Frq",pRet,ARM_ERROR_FREQ);
		return S_OK;
	}

	long resetFreqPhase1Id;
	if((resetFreqPhase1Id = ARM_ConvFrequency(l_resetFreqPhase1, C_result)) == ARM_DEFAULT_ERR)
	{
		ERROR_MSG("Invalid Reset Frq",pRet,ARM_ERROR_FREQ);
		return S_OK;
	}

	long payFreqPhase2Id;
	if((payFreqPhase2Id = ARM_ConvFrequency(l_payFreqPhase2, C_result)) == ARM_DEFAULT_ERR)
	{
		ERROR_MSG("Invalid Pay Frq Phase2",pRet,ARM_ERROR_FREQ);
		return S_OK;
	}

	long resetFreqPhase2Id;
	if((resetFreqPhase2Id = ARM_ConvFrequency(l_resetFreqPhase2, C_result)) == ARM_DEFAULT_ERR)
	{
		ERROR_MSG("Invalid reset Frq Phase2",pRet,ARM_ERROR_FREQ);
		return S_OK;
	}

	long payFreqPhase3Id;
	if((payFreqPhase3Id = ARM_ConvFrequency(l_payFreqPhase3, C_result)) == ARM_DEFAULT_ERR)
	{
		ERROR_MSG("Invalid pay Frq Phase3",pRet,ARM_ERROR_FREQ);
		return S_OK;
	}

	long resetFreqPhase3Id;
	if((resetFreqPhase3Id = ARM_ConvFrequency(l_resetFreqPhase3, C_result)) == ARM_DEFAULT_ERR)
	{
		ERROR_MSG("Invalid reset Frq Phase3",pRet,ARM_ERROR_FREQ);
		return S_OK;
	}

	long payFreqFundId;
	if((payFreqFundId = ARM_ConvFrequency(l_payFreqFund, C_result)) == ARM_DEFAULT_ERR)
	{
		ERROR_MSG("Invalid pay Frq Funding",pRet,ARM_ERROR_FREQ);
		return S_OK;
	}

	long resetFreqFundId;
	if((resetFreqFundId = ARM_ConvFrequency(l_resetFreqFund, C_result)) == ARM_DEFAULT_ERR)
	{
		ERROR_MSG("Invalid reset Frq Funding",pRet,ARM_ERROR_FREQ);
		return S_OK;
	}

	long indexPhase1Id; 
	indexPhase1Id = ARM_ConvIrType(l_indexPhase1); 

	long indexPhase3Id;
	indexPhase3Id = ARM_ConvIrType(l_indexPhase3);
	
	long indexFundId;
	indexFundId = ARM_ConvIrType(l_indexFund);

	long rollId;
	if((rollId = ARM_ConvFwdRule(l_roll, C_result)) == ARM_DEFAULT_ERR)
	{
		ERROR_MSG("Invalid Int Roll",pRet,ARM_ERROR_FREQ);
		return S_OK;
	}	
	
	long adjPhase1Id = ARM_ConvIntRule(l_adjPhase1);
	long adjPhase2Id = ARM_ConvIntRule(l_adjPhase2);
	long adjPhase2DIGId = ARM_ConvIntRule(l_adjPhase2DIG);
	long adjPhase3Id = ARM_ConvIntRule(l_adjPhase3);
	long adjFundId = ARM_ConvIntRule(l_adjFund);

	long resetTimingPhase1Id = ARM_ConvPayResetRule (l_resetTimingPhase1);
	long resetTimingPhase2DIGId = ARM_ConvPayResetRule (l_resetTimingPhase2DIG);
	long resetTimingPhase3Id = ARM_ConvPayResetRule (l_resetTimingPhase3);
	long resetTimingFundId = ARM_ConvPayResetRule (l_resetTimingFund);

	CCString l_payTimingPhase2DIG = "ARR"; 
	long payTimingPhase2DIGId = ARM_ConvPayResetRule(l_resetTimingPhase2DIG);

	long dayCountPhase1Id= ARM_ConvDayCount(l_dayCountPhase1);
	long dayCountPhase2Id= ARM_ConvDayCount(l_dayCountPhase2);
	long dayCountPhase3Id= ARM_ConvDayCount(l_dayCountPhase3);
	long dayCountFundId= ARM_ConvDayCount(l_dayCountFund);

	long stubId = ARM_ConvStubRule(l_stub);
	long ccyId; 

	long resetTimingPhase2Id = 1;

	if (ARMLOCAL_ISOCCY(l_ccy,C_result) == ARM_OK)
	{
		ccyId = C_result.getLong();
	}
	else
	{
		ERROR_MSG(C_result.getMsg(),pRet);
		return S_OK;
	}

	_bstr_t biscapped(pIsCapped);
	CCString l_iscapped = biscapped;

	long iscappedId;
	if ((iscappedId = ARM_ConvYesOrNo (l_iscapped, C_result)) == ARM_DEFAULT_ERR)
	{
	   ERROR_MSG("Invalid IsCapped",pRet);
	   return S_OK;
	}

	//************************************************************************************************
	//*********************************** Création Amort Refvalue ************************************
	//************************************************************************************************
	long amortRefValue;
	long amortRefValueToDestroy = ARM_NULL_OBJECT;

	if (LocalGetNumObjectId(l_typeAmort) == ARM_KO)
	{
		double endDateAmort = vDate[4];
		double txAmort = pTxAmort;
		double amountAmort = pAmountAmort;

		_bstr_t intRuleAmort(pIntRuleAmort);
		CCString l_intRuleAmort =  intRuleAmort;

		_bstr_t freqAmort(pFreqAmort);
		CCString l_freqAmort = freqAmort;

		_bstr_t dayCountAmort(pDayCountAmort);
		CCString l_dayCountAmort = dayCountAmort;

		long freqAmortId; 
		if((freqAmortId = ARM_ConvFrequency(l_freqAmort, C_result)) == ARM_DEFAULT_ERR)
		{
			ERROR_MSG("Invalid Frq",pRet,ARM_ERROR_FREQ);
			return S_OK;
		}

		long intRuleAmortId= ARM_ConvIntRule(l_intRuleAmort);
		long dayCountAmortId = ARM_ConvDayCount(l_dayCountAmort);	

		// AmortMethodId
		long amortMethodId;

		if((amortMethodId = ARM_ConvAmortMethod (l_typeAmort, C_result)) == ARM_DEFAULT_ERR)
		{
			ERROR_MSG("Invalid TypeAmort",pRet,ARM_ERROR_FREQ);
			return S_OK;
		}

		long FixedLegAmortId;
		if (ARMLOCAL_FIXEDLEG(startDatePhase1,
							  endDateAmort,
							  K_RCV,
							  0,
							  1.0,
							  dayCountAmortId,
							  freqAmortId,
							  K_COMP_PROP,
							  K_ARREARS,
							  intRuleAmortId,
							  stubId,
							  false,
							  l_ccy,
							  "NULL",
							  K_NX_NONE,
							  -1.0,
							  K_NO,
						  -1,
							  C_result) == ARM_OK) 
							  
		{
			FixedLegAmortId = C_result.getLong();
		}
		else
		{
			ERROR_MSG(C_result.getMsg(),pRet);
			return S_OK;
		}

		if (ARMLOCAL_ARM_GenAmortization(FixedLegAmortId,
										 amortMethodId,
										 freqAmortId,
										 amountAmort,
										 -1,
										 Ntl,
										 txAmort,
										 0.,
										 LocalGetNumObjectId(l_smiledMod),
										 0.0,
										 C_result) == ARM_OK) 
		{
			amortRefValue = C_result.getLong();
			amortRefValueToDestroy = amortRefValue;
		}
		else
		{
			ERROR_MSG(C_result.getMsg(),pRet);
			return S_OK;
		}

		retCode = ARMLOCAL_FreeObject(FixedLegAmortId,C_result);
	}
	else
	{
		amortRefValue = LocalGetNumObjectId(l_typeAmort);

		// calcul du 1er nominal
		double dPaymentDate;

		if (ARMLOCAL_DisplayRefValue(amortRefValue,true,C_result) == ARM_OK)
		{
			dPaymentDate = C_result.getArray(0);
		}
		else
		{
			ERROR_MSG(C_result.getMsg(),pRet);
			return S_OK;
		}
		
		if (ARMLOCAL_CptRefValue(amortRefValue,dPaymentDate,C_result) == ARM_OK)
		{
			Ntl = C_result.getDouble();
		}
		else
		{
			ERROR_MSG(C_result.getMsg(),pRet);
			return S_OK;
		} 
	}
	//************************************************************************************************
	//*********************************** Fin Création Amort Refvalue ************************************
	//************************************************************************************************


	//************************************************************************************************
	//*********************************** Création Funding Leg Refvalue ******************************
	//************************************************************************************************
	long irIndexFundId;
	if (ARMLOCAL_IRINDEX(KACTUAL_360,
						 payFreqFundId,
						 -1,
						 K_COMP_PROP,
						 K_MOD_FOLLOWING,
						 K_ADVANCE,
						 -2,
						 K_ARREARS,
						 10000.,
						 false,
						 l_ccy,
						 indexFundId,
						 K_COMP_PROP,
						 adjFundId,
						 resetFreqFundId,
						 C_result) == ARM_OK)
	{
		irIndexFundId = C_result.getLong();
	}
	else
	{
		ERROR_MSG(C_result.getMsg(),pRet);
		return S_OK;
	}

	long swapLegId;
	if (ARMLOCAL_SWAPLEG(irIndexFundId,
						 startDatePhase1,
						 endDatePhase3,
						 K_RCV,
						 spreadType, 
						 spreadFund, 
						 false,
						 l_ccy,
						 dayCountFundId,
						 -2,
						 "NULL",
						 "NULL",
						 0,
						 K_NX_NONE,
						 stubId,
						 -1.0,
						 K_NO,
						  -1,
						 C_result) == ARM_OK) 
	{
		swapLegId = C_result.getLong();
	}
	else
	{
		ERROR_MSG(C_result.getMsg(),pRet);
		return S_OK;
	}

	// ********** Funding Leg ************

	long cloneId;
	if (ARMLOCAL_ClonedAndSetNotional(swapLegId,
									  amortRefValue,
									  100.,
									  C_result) == ARM_OK)
	{
		cloneId = C_result.getLong();
	}
	else
	{
		ERROR_MSG(C_result.getMsg(),pRet);
		return S_OK;
	}

	//*********   Price funding Leg   *************************
	double pxFund; 

	if (ARMLOCAL_ARM_Price(cloneId,
						   LocalGetNumObjectId(l_smiledMod),
						   C_result) == ARM_OK)
	{
		pxFund = C_result.getDouble();
	}
	else
	{
		ERROR_MSG(C_result.getMsg(),pRet);
		return S_OK;
	}

	deltaFund= computeDeltaByPlot(pxFund, cloneId, vBumpBsGenMod);

	retCode = ARMLOCAL_FreeObject(irIndexFundId, C_result);
	retCode = ARMLOCAL_FreeObject(swapLegId, C_result);
	retCode = ARMLOCAL_FreeObject(cloneId,C_result);

	// Fin du calcul du prix de la funding Leg

	//************************************************************************************************
	//*********************************** Fin Création Funding Leg Refvalue **************************
	//************************************************************************************************


	
	//************************************************************************************************
	//*********************************** Création Leg Phase1 Refvalue *******************************
	//************************************************************************************************

	long irIndexPhase1Id;
	double pxPhase1; 
	double pxPhase1Delta;
	long swapLegPhase1Id;
	long clonePhase1Id;
	double resetGapPhase1; 

	if (resetTimingPhase1Id == 1)
	{
		resetGapPhase1 = -2.; 
	}
	else
	{
		resetGapPhase1 = -15.; 
	}

	if ( startDatePhase1 != startDatePhase2 )
	{
		if (ARMLOCAL_IRINDEX(KACTUAL_360,
							 payFreqPhase1Id,
							 -1,
							 K_COMP_PROP,
							 K_MOD_FOLLOWING,
							 resetTimingPhase1Id,	
							 resetGapPhase1,
							 K_ARREARS,
							 10000.,
							 false,
							 l_ccy,
							 indexPhase1Id,
							 K_COMP_PROP,
							 adjPhase1Id,
							 resetFreqPhase1Id,
							 C_result) == ARM_OK) 
		{
			irIndexPhase1Id = C_result.getLong();
		}
		else
		{
			ERROR_MSG(C_result.getMsg(),pRet);
			return S_OK;
		}
		
		if (ARMLOCAL_SWAPLEG(irIndexPhase1Id,
							 startDatePhase1,
							 startDatePhase2,
							 K_RCV,
							 0, 
							 spreadPhase1, 
							 false,
							 l_ccy,
							 dayCountPhase1Id,
							 resetGapPhase1,
							 "NULL",
							 "NULL",
							 0,
							 K_NX_NONE,
							 stubId,
							 -1.0,
							 K_NO,
						  -1,
							 C_result) == ARM_OK) 
		{
			swapLegPhase1Id = C_result.getLong();
		}
		else
		{
			ERROR_MSG(C_result.getMsg(),pRet);
			return S_OK;
		}

		// ********** Phase1 Leg ************

		if (ARMLOCAL_ClonedAndSetNotional(swapLegPhase1Id,
										  amortRefValue,
										  100.,
										  C_result) == ARM_OK)
		{
			clonePhase1Id = C_result.getLong();
		}
		else
		{
			ERROR_MSG(C_result.getMsg(),pRet);
			return S_OK;
		}
		
		//*********   Price Phase1 *************************

		if (ARMLOCAL_ARM_Price(clonePhase1Id,
							   LocalGetNumObjectId(l_smiledMod),
							   C_result) == ARM_OK)
		{
			pxPhase1 = C_result.getDouble();
		}
		else
		{
			ERROR_MSG(C_result.getMsg(),pRet);
			return S_OK;
		}
		
		deltaPhase1 = computeDeltaByPlot(pxPhase1, clonePhase1Id, vBumpBsGenMod);

		retCode = ARMLOCAL_FreeObject(irIndexPhase1Id,C_result);
		retCode = ARMLOCAL_FreeObject(swapLegPhase1Id,C_result);
		retCode = ARMLOCAL_FreeObject(clonePhase1Id,C_result);
	}
	else
	{
		pxPhase1 = 0.;
		pxPhase1Delta = 0.;
	}

	//************************************************************************************************
	//************************** Fin Création Leg Phase1 Refvalue ************************************
	//************************************************************************************************


	//************************************************************************************************
	//*********************************** Création Leg Phase2 Refvalue *******************************
	//************************************************************************************************

	long irIndexPhase2Id;
	double resetGapPhase2; 

	if (resetTimingPhase2Id == 1)
	{
		resetGapPhase2 = -2.; 
	}
	else
	{
		resetGapPhase2 = -15.; 
	}

	if (ARMLOCAL_IRINDEX(KACTUAL_360,
						 payFreqPhase2Id,
						 -1,
						 K_COMP_PROP,
						 K_MOD_FOLLOWING,
						 resetTimingPhase2Id,	
						 resetGapPhase2,
						 K_ARREARS,
						 10000.,
						 false,
						 l_ccy,
						 K_FIXED,
						 K_COMP_PROP,
						 adjPhase2Id,
						 resetFreqPhase2Id,
						 C_result) == ARM_OK) 
	{
		irIndexPhase2Id = C_result.getLong();
	}
	else
	{
		ERROR_MSG(C_result.getMsg(),pRet);
		return S_OK;
	}

	//************************************************************************************************
	//************************** Fin Création Leg Phase2 Refvalue ************************************
	//************************************************************************************************



	//************************************************************************************************
	//*************************** Création Leg Phase3 Refvalue ***************************************
	//************************************************************************************************

	long irIndexPhase3Id;
	double pxPhase3; 
	double pxPhase3Delta;
	long swapLegPhase3Id;
	long clonePhase3Id;
	double resetGapPhase3; 

	if ( startDatePhase3 != endDatePhase3 )
	{		
		if (resetTimingPhase3Id == 1)
		{
			resetGapPhase3 = -2.; 
		}
		else
		{
			resetGapPhase3 = -15.; 
		}

		if (ARMLOCAL_IRINDEX(KACTUAL_360,
							 payFreqPhase3Id,
							 -1,
							 K_COMP_PROP,
							 K_MOD_FOLLOWING,
							 resetTimingPhase3Id,	
							 resetGapPhase3,
							 K_ARREARS,
							 10000.,
							 false,
							 l_ccy,
							 indexPhase3Id,
							 K_COMP_PROP,
							 adjPhase3Id,
							 resetFreqPhase3Id,
							 C_result) == ARM_OK) 
		{
			irIndexPhase3Id = C_result.getLong();
		}
		else
		{
			ERROR_MSG(C_result.getMsg(),pRet);
			return S_OK;
		}

		if (ARMLOCAL_SWAPLEG(irIndexPhase3Id,
							 startDatePhase3,
							 endDatePhase3,
							 K_RCV,
							 0, 
							 spreadPhase3, 
							 false,
							 l_ccy,
							 dayCountPhase3Id,
							 resetGapPhase3,
							 "NULL",
							 "NULL",
							 0,
							 K_NX_NONE,
							 stubId,
							 -1.0,
							 K_NO,
						  -1,
							 C_result) == ARM_OK) 
		{
			swapLegPhase3Id = C_result.getLong();
		}
		else
		{
			ERROR_MSG(C_result.getMsg(),pRet);
			return S_OK;
		}		

		// ********** Phase3 Leg ************

		if (ARMLOCAL_ClonedAndSetNotional(swapLegPhase3Id,
										  amortRefValue,
										  100.,
										  C_result) == ARM_OK)
		{
			clonePhase3Id = C_result.getLong();
		}
		else
		{
			ERROR_MSG(C_result.getMsg(),pRet);
			return S_OK;
		}

		//*********   Price Phase3 *************************	
		if (ARMLOCAL_ARM_Price(clonePhase3Id,
							   LocalGetNumObjectId(l_smiledMod),
							   C_result) == ARM_OK)
		{
			pxPhase3 = C_result.getDouble();
		}
		else
		{
			ERROR_MSG(C_result.getMsg(),pRet);
			return S_OK;
		}

		deltaPhase3 = computeDeltaByPlot(pxPhase3, clonePhase3Id, vBumpBsGenMod);

		retCode = ARMLOCAL_FreeObject(irIndexPhase3Id,C_result);
		retCode = ARMLOCAL_FreeObject(swapLegPhase3Id,C_result);
		retCode = ARMLOCAL_FreeObject(clonePhase3Id,C_result);
	}
	else
	{
		pxPhase3 = 0.;
		pxPhase3Delta = 0.; 
	}

	sizeMatu = vBumpBsGenMod.size(); 
	long i; 

	for (i =0; i< sizeMatu; i++)
	{
		deltaOther.push_back(-deltaFund[i]); 
	}

	if ( startDatePhase1 != startDatePhase2 )
	{
		for (i =0; i< sizeMatu; i++)
		{
			deltaOther[i] = deltaOther[i] + deltaPhase1[i]; 
		}
	}

	if (startDatePhase3 != endDatePhase3 )
	{
		for (i =0; i< sizeMatu; i++)
		{
			deltaOther[i] = deltaOther[i] + deltaPhase3[i]; 
		}
	}

	double pxPhase2DIG; 
	double pxTotal; 
	double pxDelta; 
	double SensiTotal; 

	double prixInter;
	double Xincr;
	
	double resetGapPhase2DIG; 

	if (resetTimingPhase2DIGId == 1)
	{
		resetGapPhase2DIG = -2.; 
	}
	else
	{
		resetGapPhase2DIG = -15.; 
	}

//******************** Calcul premiere PV ******************
	int iteration =1;
	double prix;
	double prixAct; 
	double prixAct2;
	double SprSeek = 10.01;
	double SprSeekFinal= 10.01;
	double SprSeekBonifie = 10.01;  // PUREMENT ARBITRAIRE
	fee = Ntl * fee / 100 ;
	VECTOR<double> vegaSpdopt;
	VECTOR<double> deltaSpdopt;
	VECTOR<double> vegaSpdoptCap;
	VECTOR<double> deltaSpdoptCap;
//******************** Calcul premiere PV ******************
	
	VECTOR<double> pxReviPentix = computeReviPentix(l_indexPhase2DIG,l_indexLongPhase2DIG, irIndexPhase2Id, startDatePhase2, startDatePhase3, 
			l_floorOrCap, l_strikePhase2DIG, SprSeekBonifie, Levier, TxFixeMax , ccyId, 
			dayCountPhase2Id, 
			resetTimingPhase2DIGId, payTimingPhase2DIGId,  resetGapPhase2DIG, resetGapPhase2, resetFreqPhase2Id, payFreqPhase2Id, 
			adjPhase2DIGId, stubId,amortRefValue, 
				deltaOther, 
			iscappedId, TxCap,
			vBumpBsGenMod, vBumpVolBsGenMod, 
			l_hyperCubeCorrel, l_volCurvFromMatriceShift,
			l_vol, l_volCub, l_corrManager,
			l_convexityManager, l_zc, l_smiledMod,&vegaSpdopt, &deltaSpdopt,&vegaSpdoptCap, &deltaSpdoptCap);

	if (pxReviPentix.size() == 0)
	{
		ERROR_MSG("Pb in computeReviPentix",pRet);
		return S_OK;
	}

	pxPhase2DIG = pxReviPentix[0];
	double vega =pxReviPentix[1]; 
	double deltaTotal = pxReviPentix[2];
	double pxPhase2DIGBumpCorrel = pxReviPentix[3];
	double MT = pxReviPentix[4];

	pxTotal = pxPhase1 + pxPhase3 - pxFund + pxPhase2DIG;
	
	if ((MT - (fabs(vega) + fabs(deltaTotal) + fabs(pxPhase2DIGBumpCorrel))) > 0) 
	{
		pxDelta =MT ; 
	}
	else
	{
		pxDelta = fabs(pxPhase2DIGBumpCorrel) + fabs(vega) + fabs(deltaTotal); 

	}
	//pxDelta = fabs(pxPhase2DIGBumpCorrel) + MT + vega + deltaTotal; 
	if (pxDelta <0.) pxDelta = - pxDelta;

	SensiTotal = pxDelta;
	prixAct = pxTotal - SensiTotal - fee; 

	prix = prixAct ;
	prixAct2 =1;
	
//********************                    ******************
//******************** BOUCLE SUR LE PRICE ******************
	while ((abs(prixAct2) > 0.0001) && (iteration < 10 ))
	{
		SprSeek = SprSeekBonifie;
		
		//*************************
		// **** Premier recalcul ****
		//****************************

		pxReviPentix = computeReviPentix(l_indexPhase2DIG,l_indexLongPhase2DIG, irIndexPhase2Id, startDatePhase2, startDatePhase3, 
			l_floorOrCap, l_strikePhase2DIG, SprSeekBonifie, Levier, TxFixeMax, ccyId, 
			dayCountPhase2Id, 
			resetTimingPhase2DIGId, payTimingPhase2DIGId,  resetGapPhase2DIG, resetGapPhase2, resetFreqPhase2Id, payFreqPhase2Id, 
			adjPhase2DIGId, stubId, amortRefValue,
			deltaOther, 
			iscappedId, TxCap,
			vBumpBsGenMod, vBumpVolBsGenMod, 
			l_hyperCubeCorrel, l_volCurvFromMatriceShift,
			l_vol, l_volCub, l_corrManager,
			l_convexityManager, l_zc, l_smiledMod,&vegaSpdopt, &deltaSpdopt,&vegaSpdoptCap, &deltaSpdoptCap);

		if (pxReviPentix.size() == 0)
		{
			ERROR_MSG("Pb in computeReviPentix",pRet);
			return S_OK;
		}

		pxPhase2DIG = pxReviPentix[0];
		vega =pxReviPentix[1]; 
		deltaTotal = pxReviPentix[2];
		pxPhase2DIGBumpCorrel = pxReviPentix[3];
		MT = pxReviPentix[4];

		pxTotal = pxPhase1 + pxPhase3 - pxFund + pxPhase2DIG;
		if ((MT - (fabs(vega) + fabs(deltaTotal) + fabs(pxPhase2DIGBumpCorrel))) > 0) 
		{
			pxDelta =  MT ; 
		}
		else
		{
			pxDelta = fabs(pxPhase2DIGBumpCorrel) + fabs(vega) + fabs(deltaTotal); 

		}
		//pxDelta = fabs(pxPhase2DIGBumpCorrel) + MT+ vega + deltaTotal;
		if (pxDelta <0.) pxDelta = - pxDelta;

		SensiTotal = pxDelta;
		prixAct = pxTotal - SensiTotal - fee; 

		// les deltas se compensent et on somme en valeur absolue les deltas et vegas

		//******************************
		// **** affectation valeurs ****
		//******************************

		prixAct2 = prixAct ;
		Xincr = SprSeek / 1000000;
		SprSeekBonifie = SprSeek + Xincr;		

		
		//****************************
		// **** Deuxieme recalcul ****
		//****************************
		pxReviPentix = computeReviPentix(l_indexPhase2DIG,l_indexLongPhase2DIG, irIndexPhase2Id, startDatePhase2, startDatePhase3, 
									l_floorOrCap, l_strikePhase2DIG, SprSeekBonifie, Levier, TxFixeMax, ccyId, 
									dayCountPhase2Id, 
									resetTimingPhase2DIGId, payTimingPhase2DIGId,  resetGapPhase2DIG, resetGapPhase2, resetFreqPhase2Id, payFreqPhase2Id, 
									adjPhase2DIGId, stubId, amortRefValue, 
									deltaOther,
									iscappedId, TxCap,
									vBumpBsGenMod, vBumpVolBsGenMod, 
									l_hyperCubeCorrel, l_volCurvFromMatriceShift, 
									l_vol, l_volCub, l_corrManager, 
									l_convexityManager, l_zc, l_smiledMod,&vegaSpdopt, &deltaSpdopt,&vegaSpdoptCap, &deltaSpdoptCap);

		if (pxReviPentix.size() == 0)
		{
			ERROR_MSG("Pb in computeReviPentix",pRet);
			return S_OK;
		}

		pxPhase2DIG = pxReviPentix[0];
		vega =pxReviPentix[1]; 
		deltaTotal = pxReviPentix[2];
		pxPhase2DIGBumpCorrel = pxReviPentix[3];
		MT = pxReviPentix[4];

		pxTotal = pxPhase1 + pxPhase3 - pxFund + pxPhase2DIG;
		//pxDelta = fabs(pxPhase2DIGBumpCorrel) + MT+ vega + deltaTotal;
		if ((MT - (fabs(vega) + fabs(deltaTotal) + fabs(pxPhase2DIGBumpCorrel))) > 0) 
		{
			pxDelta =  MT; 
		}
		else
		{
			pxDelta = fabs(pxPhase2DIGBumpCorrel) + fabs(vega) + fabs(deltaTotal); 

		}

		if (pxDelta <0.) pxDelta = - pxDelta;

		SensiTotal = pxDelta;
		prixAct = pxTotal - SensiTotal - fee; 
	
		iteration++;
		
		SprSeekFinal = SprSeekBonifie;
		prixInter = ( prixAct - prixAct2 ) / Xincr;
		// if prixinter =  0 sortir erreur
		SprSeekBonifie = SprSeek - prixAct / prixInter ; 
		prixAct2 = prixAct;
	}

//**************** PRICE *******************

	retCode = ARMLOCAL_FreeObject(irIndexPhase2Id,C_result);
	retCode = ARMLOCAL_FreeObject(amortRefValueToDestroy,C_result);
	retCode = ARMLOCAL_FreeObject(ccyId,C_result);

	VECTOR<double> vRes;

	vRes.push_back(SprSeekFinal);

	vRes.push_back(SprSeekFinal);

	Res = VECTORDOUBLE2VARIANT(vRes, pRet);

	_variant_t wrap_Res;
	wrap_Res.Attach(*pRet);
	*pRet=wrap_Res.Detach();

	return S_OK;
}
catch(...)
{
	return E_FAIL;
}
}

VECTOR<double> computeLeg2LivretA(double startDate, double endDate,
							    double spread, CCString l_ccy, long dayCountId,
							   long resetTimingId, double resetGap, long resetFreqId, long payFreqId, 
							   long stubRuleId, long amortRefValue,
							   const CCString& l_LAMod, const CCString& l_LAModBumpInflation, const CCString& l_LAModBump, 
							   VECTOR<CCString>& C_resetMgrId)
{
	VECTOR<double> resultat;
	resultat.clear();

try
{
	long retCode; 
	ARM_result C_result;

	long swapLegLAId; 
	long cloneLAId; 
	double pxLA; 

	double asOf;
	if (ARMLOCAL_GetAsOf(LocalGetNumObjectId(l_LAMod),
						 C_result) == ARM_OK)
	{
		asOf = C_result.getDouble();
	}

	if (ARMLOCAL_LIVRETALEG(startDate,
						 endDate,
						 K_LIVRET_A,
						 K_RCV,
						 0, 
						 spread, 
						 resetFreqId,
						 payFreqId, 
						 resetTimingId,
						 K_ARREARS, 	
						 false,
						 l_ccy,
						 K_UNADJUSTED, 
						 0,
						 "NULL",
						 "NULL",
						 0,
						 0,
						 stubRuleId,
						 -1.0,
						 K_NO,
						 dayCountId,
						 C_result) == ARM_OK) 
	{
		swapLegLAId = C_result.getLong();
	}

	// **********  Leg LA  ************	
	long PastResetcloneLAId; 
	
	if (ARMLOCAL_ClonedAndSetNotional(swapLegLAId,
									  amortRefValue,
									  100.,
									  C_result) == ARM_OK)
	{
		cloneLAId = C_result.getLong();
	}

	retCode = ARMLOCAL_DoPastReset(cloneLAId,
								   C_resetMgrId,
								   JulianToXLDate(asOf),
								   C_result);

	if(retCode == ARM_OK)
	{
		PastResetcloneLAId = C_result.getLong ();
	}
	
	//*********   Price LA *************************


	if (ARMLOCAL_ARM_Price(PastResetcloneLAId,
						   LocalGetNumObjectId(l_LAMod),
						   C_result) == ARM_OK)
	{
		pxLA = C_result.getDouble();
	}

	double pxLADeltaInflation; 
	if (ARMLOCAL_ARM_Price(PastResetcloneLAId,
						   LocalGetNumObjectId(l_LAModBumpInflation),
						   C_result) == ARM_OK)
	{
		pxLADeltaInflation = C_result.getDouble();
	}

	double pxLADeltaEuribor; 
	if (ARMLOCAL_ARM_Price(PastResetcloneLAId,
						   LocalGetNumObjectId(l_LAModBump),
						   C_result) == ARM_OK)
	{
		pxLADeltaEuribor = C_result.getDouble();
	}

	retCode = ARMLOCAL_FreeObject(swapLegLAId,C_result);
	retCode = ARMLOCAL_FreeObject(cloneLAId,C_result);
	retCode = ARMLOCAL_FreeObject(PastResetcloneLAId,C_result);

	//*********   Price Première Phase structure 1  TD*************************
	VECTOR<double> vRes;
	
	vRes.push_back(pxLA);
	vRes.push_back(pxLADeltaInflation - pxLA); 
	vRes.push_back(pxLADeltaEuribor - pxLA);
	
	return vRes;
}
catch(...)
{
	return resultat;
}
}



STDMETHODIMP ActiveXModule::ARMcomputeLivretA(double pNtl, 
										   double pStartDateLeg1,  
										   double pEndDateLeg1,
										   BSTR pCcy,
										   BSTR pIndexLeg1, 										   
										   VARIANT pSpreadLeg1,
										   BSTR pDayCountLeg1,
										   BSTR pPayFreqLeg1,
										   BSTR pResetFreqLeg1,
										   BSTR pResetTimingLeg1,
										   BSTR pAdjLeg1,
										   BSTR pRoll, 
										   BSTR pStub,
										   double pEndDateLA,
										   double pSpreadLeg2,
										   BSTR pDayCountLA,
										   BSTR pPayFreqLA,
										   BSTR pResetFreqLA,
										   BSTR pResetTimingLA,									   
										   BSTR pAdjLA, 
										   BSTR pIndexLeg2, 
										   BSTR pDayCountLeg2, 
										   BSTR pPayFreqLeg2, 
										   BSTR pResetFreqLeg2,
										   BSTR pResetTimingLeg2,
										   BSTR pAdjLeg2, 
										   double pEndDateAmort,
										   BSTR pDayCountAmort,
										   BSTR pIntRuleAmort, 
										   double pTxAmort,
										   BSTR pFreqAmort,
										   double pAmountAmort, 
										   BSTR pTypeAmort, 
										   double pFee, 
										   BSTR pSmiledMod, 
										   BSTR pSmiledModBump, 
										   BSTR pLAMod, 
										   BSTR pLAModBump,
										   BSTR pLAModBumpInflation,
										   VARIANT *pResetMgrIds, 
										   VARIANT *pRet
										   )

{	
	// TODO: Add your implementation code here
try
{									   
	ARM_result C_result;
	double startDateLeg1 = pStartDateLeg1; 
	double endDateLeg1 = pEndDateLeg1;
	long spreadType; 
	double spreadLeg1; 
	
	if (VARIANT2Double(pSpreadLeg1,&spreadLeg1) == S_FALSE)
	{
		_bstr_t bpSpreadLeg1(pSpreadLeg1);
		CCString l_SpreadLeg1= bpSpreadLeg1;
		spreadLeg1 = (double) LocalGetNumObjectId(l_SpreadLeg1);
		spreadType = 1;
	}
	else
	{
	  spreadType = 0;

	}
	double endDateLA = pEndDateLA;
	double spreadLeg2= pSpreadLeg2;

	double Ntl = pNtl;
	double fee = pFee; 

	long retCode;
	double Res = 0.;

	_bstr_t ccy(pCcy);
	CCString l_ccy = ccy;

	_bstr_t indexLeg1(pIndexLeg1);
	CCString l_indexLeg1 = indexLeg1;

	_bstr_t dayCountLeg1(pDayCountLeg1);
	CCString l_dayCountLeg1 = dayCountLeg1;

	_bstr_t resetFreqLeg1(pResetFreqLeg1);
	CCString l_resetFreqLeg1 = resetFreqLeg1;

	_bstr_t payFreqLeg1(pPayFreqLeg1);
	CCString l_payFreqLeg1 = payFreqLeg1;

	_bstr_t resetTimingLeg1(pResetTimingLeg1);
	CCString l_resetTimingLeg1 = resetTimingLeg1;

	_bstr_t roll(pRoll);
	CCString l_roll = roll;

	_bstr_t adjLeg1(pAdjLeg1);
	CCString l_adjLeg1 = adjLeg1;
	
	_bstr_t stub(pStub);
	CCString l_stub= stub;

	//Strike is double??
	CCString l_indexLA = "LIVRETA"; 

	_bstr_t payFreqLA(pPayFreqLA);
	CCString l_payFreqLA = payFreqLA;

	_bstr_t resetFreqLA(pResetFreqLA);
	CCString l_resetFreqLA = resetFreqLA;

	_bstr_t resetTimingLA(pResetTimingLA);
	CCString l_resetTimingLA = resetTimingLA;

	_bstr_t adjLA(pAdjLA);
	CCString l_adjLA= adjLA;

	_bstr_t dayCountLA(pDayCountLA);
	CCString l_dayCountLA = dayCountLA;

	_bstr_t resetTimingLeg2(pResetTimingLeg2);
	CCString l_resetTimingLeg2 = resetTimingLeg2;

	_bstr_t dayCountLeg2(pDayCountLeg2);
	CCString l_dayCountLeg2= dayCountLeg2;

	_bstr_t payFreqLeg2(pPayFreqLeg2);
	CCString l_payFreqLeg2= payFreqLeg2;
	
	_bstr_t resetFreqLeg2(pResetFreqLeg2);
	CCString l_resetFreqLeg2 = resetFreqLeg2;

	_bstr_t adjLeg2(pAdjLeg2);
	CCString l_adjLeg2 = adjLeg2;

	_bstr_t indexLeg2(pIndexLeg2);
	CCString l_indexLeg2 = indexLeg2;

	_bstr_t typeAmort(pTypeAmort);
	CCString l_typeAmort = typeAmort;

	_bstr_t smiledMod (pSmiledMod);
	CCString l_smiledMod = smiledMod;

	_bstr_t smiledModBump(pSmiledModBump); 
	CCString l_smiledModBump = smiledModBump;
	
	_bstr_t LAMod(pLAMod);
	CCString l_LAMod = LAMod;
	
	_bstr_t LAModBump(pLAModBump);
	CCString l_LAModBump = LAModBump;

	_bstr_t LAModBumpInflation(pLAModBumpInflation);
	CCString l_LAModBumpInflation = LAModBumpInflation;
	
	long sizeMgr; 
	vector <CCString> vResetMgrIds; 

	if(VARIANT2VECTORCCSTRING (*pResetMgrIds,vResetMgrIds,sizeMgr) != S_OK)
	{
		ERROR_MSG("Error recuperation resetManager",pRet,ARM_ERROR_FREQ);
		return S_OK;
	}

	long payFreqLeg1Id; 	
	if((payFreqLeg1Id = ARM_ConvFrequency(l_payFreqLeg1, C_result)) == ARM_DEFAULT_ERR)
	{
		ERROR_MSG("Invalid Frq",pRet,ARM_ERROR_FREQ);
		return S_OK;
	}

	long resetFreqLeg1Id;
	if((resetFreqLeg1Id = ARM_ConvFrequency(l_resetFreqLeg1, C_result)) == ARM_DEFAULT_ERR)
	{
		ERROR_MSG("Invalid Frq",pRet,ARM_ERROR_FREQ);
		return S_OK;
	}

	long payFreqLeg2Id;
	if((payFreqLeg2Id = ARM_ConvFrequency(l_payFreqLeg2, C_result)) == ARM_DEFAULT_ERR)
	{
		ERROR_MSG("Invalid Frq",pRet,ARM_ERROR_FREQ);
		return S_OK;
	}

	long resetFreqLeg2Id;
	if((resetFreqLeg2Id = ARM_ConvFrequency(l_resetFreqLeg2, C_result)) == ARM_DEFAULT_ERR)
	{
		ERROR_MSG("Invalid Frq",pRet,ARM_ERROR_FREQ);
		return S_OK;
	}

	long payFreqLAId;
	if((payFreqLAId = ARM_ConvFrequency(l_payFreqLA, C_result)) == ARM_DEFAULT_ERR)
	{
		ERROR_MSG("Invalid Frq",pRet,ARM_ERROR_FREQ);
		return S_OK;
	}

	long resetFreqLAId;
	if((resetFreqLAId = ARM_ConvFrequency(l_resetFreqLA, C_result)) == ARM_DEFAULT_ERR)
	{
		ERROR_MSG("Invalid Frq",pRet,ARM_ERROR_FREQ);
		return S_OK;
	}

	long indexLeg1Id; 
	indexLeg1Id = ARM_ConvIrType(l_indexLeg1); 

	long indexLeg2Id; 
	indexLeg2Id = ARM_ConvIrType(l_indexLeg2); 

	long indexLAId; 
	indexLAId = K_LIVRET_A; 

	long rollId;
	if((rollId = ARM_ConvFwdRule(l_roll, C_result)) == ARM_DEFAULT_ERR)
	{
		ERROR_MSG("Invalid Int Roll",pRet,ARM_ERROR_FREQ);
		return S_OK;
	}	
	
	long adjLeg1Id = ARM_ConvIntRule(l_adjLeg1);
	long adjLeg2Id = ARM_ConvIntRule(l_adjLeg2);
	long adjLAId = ARM_ConvIntRule(l_adjLA);

	long resetTimingLeg1Id = ARM_ConvPayResetRule (l_resetTimingLeg1);
	long resetTimingLeg2Id = ARM_ConvPayResetRule (l_resetTimingLeg2);
	long resetTimingLAId = ARM_ConvPayResetRule (l_resetTimingLA);

	long dayCountLeg1Id= ARM_ConvDayCount(l_dayCountLeg1);
	long dayCountLeg2Id= ARM_ConvDayCount(l_dayCountLeg2);
	long dayCountLAId= ARM_ConvDayCount(l_dayCountLA);

	long stubId = ARM_ConvStubRule(l_stub);
	long ccyId; 

	if (ARMLOCAL_ISOCCY(l_ccy,C_result) == ARM_OK)
	{
		ccyId = C_result.getLong();
	}
	else
	{
		ERROR_MSG(C_result.getMsg(),pRet);
		return S_OK;
	}

	//************************************************************************************************
	//*********************************** Création Amort Refvalue ************************************
	//************************************************************************************************
	long amortRefValue;
	long amortRefValueToDestroy = ARM_NULL_OBJECT;

	if (LocalGetNumObjectId(l_typeAmort) == ARM_KO)
	{
		double endDateAmort = pEndDateAmort;
		double txAmort = pTxAmort;
		double amountAmort = pAmountAmort;

		_bstr_t intRuleAmort(pIntRuleAmort);
		CCString l_intRuleAmort =  intRuleAmort;

		_bstr_t freqAmort(pFreqAmort);
		CCString l_freqAmort = freqAmort;

		_bstr_t dayCountAmort(pDayCountAmort);
		CCString l_dayCountAmort = dayCountAmort;

		long freqAmortId; 
		if((freqAmortId = ARM_ConvFrequency(l_freqAmort, C_result)) == ARM_DEFAULT_ERR)
		{
			ERROR_MSG("Invalid Frq",pRet,ARM_ERROR_FREQ);
			return S_OK;
		}

		long intRuleAmortId= ARM_ConvIntRule(l_intRuleAmort);
		long dayCountAmortId = ARM_ConvDayCount(l_dayCountAmort);	

		// AmortMethodId
		long amortMethodId;

		if((amortMethodId = ARM_ConvAmortMethod (l_typeAmort, C_result)) == ARM_DEFAULT_ERR)
		{
			ERROR_MSG("Invalid TypeAmort",pRet,ARM_ERROR_FREQ);
			return S_OK;
		}

		long FixedLegAmortId;
		if (ARMLOCAL_FIXEDLEG(startDateLeg1,
							  endDateAmort,
							  K_RCV,
							  0,
							  1.0,
							  dayCountAmortId,
							  freqAmortId,
							  K_COMP_PROP,
							  K_ARREARS,
							  intRuleAmortId,
							  stubId,
							  false,
							  l_ccy,
							  "NULL",
							  K_NX_NONE,
							  -1.0,
							  K_NO,
								 -1,
							  C_result) == ARM_OK) 
							  
		{
			FixedLegAmortId = C_result.getLong();
		}
		else
		{
			ERROR_MSG(C_result.getMsg(),pRet);
			return S_OK;
		}

		if (ARMLOCAL_ARM_GenAmortization(FixedLegAmortId,
										 amortMethodId,
										 freqAmortId,
										 amountAmort,
										 -1,
										 Ntl,
										 txAmort,
										 0.,
										 LocalGetNumObjectId(l_smiledMod),
										 0.0,
										 C_result) == ARM_OK) 
		{
			amortRefValue = C_result.getLong();
			amortRefValueToDestroy = amortRefValue;
		}
		else
		{
			ERROR_MSG(C_result.getMsg(),pRet);
			return S_OK;
		}

		retCode = ARMLOCAL_FreeObject(FixedLegAmortId,C_result);
	}
	else
	{
		amortRefValue = LocalGetNumObjectId(l_typeAmort);

		// calcul du 1er nominal
		double dPaymentDate;

		if (ARMLOCAL_DisplayRefValue(amortRefValue,true,C_result) == ARM_OK)
		{
			dPaymentDate = C_result.getArray(0);
		}
		else
		{
			ERROR_MSG(C_result.getMsg(),pRet);
			return S_OK;
		}
		
		if (ARMLOCAL_CptRefValue(amortRefValue,dPaymentDate,C_result) == ARM_OK)
		{
			Ntl = C_result.getDouble();
		}
		else
		{
			ERROR_MSG(C_result.getMsg(),pRet);
			return S_OK;
		} 
	}
	//************************************************************************************************
	//*********************************** Fin Création Amort Refvalue ********************************
	//************************************************************************************************
	
	//************************************************************************************************
	//*********************************** Création Leg Phase1 Refvalue *******************************
	//************************************************************************************************
	long irIndexLeg1Id;
	double pxLeg1; 
	double pxLeg1Delta;
	long swapLeg1Id;
	long cloneLeg1Id;
	double resetGapLeg1; 

	if (resetTimingLeg1Id == 1)
	{
		resetGapLeg1 = -2.; 
	}
	else
	{
		resetGapLeg1 = -15.; 
	}

	if (ARMLOCAL_IRINDEX(KACTUAL_360,
						 payFreqLeg1Id,
						 -1,
						 K_COMP_PROP,
						 K_MOD_FOLLOWING,
						 resetTimingLeg1Id,	
						 resetGapLeg1,
						 K_ARREARS,
						 10000.,
						 false,
						 l_ccy,
						 indexLeg1Id,
						 K_COMP_PROP,
						 adjLeg1Id,
						 resetFreqLeg1Id,
						 C_result) == ARM_OK) 
	{
		irIndexLeg1Id = C_result.getLong();
	}
	else
	{
		ERROR_MSG(C_result.getMsg(),pRet);
		return S_OK;
	}
	
	if (ARMLOCAL_SWAPLEG(irIndexLeg1Id,
						 startDateLeg1,
						 endDateLeg1,
						 K_RCV,
						 spreadType, 
						 spreadLeg1, 
						 false,
						 l_ccy,
						 dayCountLeg1Id,
						 resetGapLeg1,
						 "NULL",
						 "NULL",
						 0,
						 K_NX_NONE,
						 stubId,
						 -1.0,
						 K_NO,
								 -1,
						 C_result) == ARM_OK) 
	{
		swapLeg1Id = C_result.getLong();
	}
	else
	{
		ERROR_MSG(C_result.getMsg(),pRet);
		return S_OK;
	}

	// ********** Phase1 Leg ************
	if (ARMLOCAL_ClonedAndSetNotional(swapLeg1Id,
									  amortRefValue,
									  100.,
									  C_result) == ARM_OK)
	{
		cloneLeg1Id = C_result.getLong();
	}
	else
	{
		ERROR_MSG(C_result.getMsg(),pRet);
		return S_OK;
	}

	//*********   Price Phase1 *************************
	if (ARMLOCAL_ARM_Price(cloneLeg1Id,
						   LocalGetNumObjectId(l_smiledMod),
						   C_result) == ARM_OK)
	{
		pxLeg1 = C_result.getDouble();
	}
	else
	{
		ERROR_MSG(C_result.getMsg(),pRet);
		return S_OK;
	}

	// shift seulement sur l'euro !?
	if (ARMLOCAL_ARM_Price(cloneLeg1Id,
						   LocalGetNumObjectId(l_smiledModBump),
						   C_result) == ARM_OK)
	{
		pxLeg1Delta = C_result.getDouble();
	}
	else
	{
		ERROR_MSG(C_result.getMsg(),pRet);
		return S_OK;
	}

	retCode = ARMLOCAL_FreeObject(irIndexLeg1Id,C_result);
	retCode = ARMLOCAL_FreeObject(swapLeg1Id,C_result);
	retCode = ARMLOCAL_FreeObject(cloneLeg1Id,C_result);

	//************************************************************************************************
	//************************** Fin Création Leg Phase1 Refvalue ************************************
	//************************************************************************************************

	//************************************************************************************************
	//*********************************** Création Leg Livret A Refvalue *****************************
	//************************************************************************************************
	long irIndexLAId;
	double pxLA; 
	double resetGapLA; 

	if (resetTimingLAId == 1)
	{
		resetGapLA = -2.; 
	}
	else
	{
		resetGapLA = -15.; 
	}

	if (ARMLOCAL_IRINDEX(KACTUAL_360,
						 payFreqLAId,
						 -1,
						 K_COMP_PROP,
						 K_MOD_FOLLOWING,
						 resetTimingLAId,	
						 resetGapLA,
						 K_ARREARS,
						 10000.,
						 false,
						 l_ccy,
						 indexLAId,
						 K_COMP_PROP,
						 adjLAId,
						 resetFreqLAId,
						 C_result) == ARM_OK) 
	{
		irIndexLAId = C_result.getLong();
	}
	else
	{
		ERROR_MSG(C_result.getMsg(),pRet);
		return S_OK;
	}

	//************************************************************************************************
	//************************** Fin Création Leg Phase1 Refvalue ************************************
	//************************************************************************************************

	//************************************************************************************************
	//*********************************** Création Leg2 Phase2 Refvalue *******************************
	//************************************************************************************************
	long irIndexLeg2Id;
	double resetGapLeg2;
	long swapLeg2Id; 
	long cloneLeg2Id; 

	double pxLeg2Delta; 
	double pxLeg2; 
	if (resetTimingLeg2Id == 1)
	{
		resetGapLeg2 = -2.; 
	}
	else
	{
		resetGapLeg2 = -15.; 
	}

	if ( endDateLA != endDateLeg1 )
	{
		if (ARMLOCAL_IRINDEX(KACTUAL_360,
							 payFreqLeg2Id,
							 -1,
							 K_COMP_PROP,
							 K_MOD_FOLLOWING,
							 resetTimingLeg2Id,	
							 resetGapLeg2,
							 K_ARREARS,
							 10000.,
							 false,
							 l_ccy,
							 indexLeg2Id,
							 K_COMP_PROP,
							 adjLeg2Id,
							 resetFreqLeg2Id,
							 C_result) == ARM_OK) 
		{
			irIndexLeg2Id = C_result.getLong();
		}
		else
		{
			ERROR_MSG(C_result.getMsg(),pRet);
			return S_OK;
		}

		if (ARMLOCAL_SWAPLEG(irIndexLeg2Id,
							endDateLA, 
							endDateLeg1, 						
							 K_RCV,
							 0, 
							 spreadLeg2, 
							 false,
							 l_ccy,
							 dayCountLeg2Id,
							 resetGapLeg2,
							 "NULL",
							 "NULL",
							 0,
							 K_NX_NONE,
							 stubId,
							 -1.0,
							 K_NO,
								 -1,
							 C_result) == ARM_OK) 
		{
			swapLeg2Id = C_result.getLong();
		}
		else
		{
			ERROR_MSG(C_result.getMsg(),pRet);
			return S_OK;
		}

		// ********** Phase2 Leg ************
		if (ARMLOCAL_ClonedAndSetNotional(swapLeg2Id,
										  amortRefValue,
										  100.,
										  C_result) == ARM_OK)
		{
			cloneLeg2Id = C_result.getLong();
		}
		else
		{
			ERROR_MSG(C_result.getMsg(),pRet);
			return S_OK;
		}

		//*********   Price Phase2 *************************
		if (ARMLOCAL_ARM_Price(cloneLeg2Id,
							   LocalGetNumObjectId(l_smiledMod),
							   C_result) == ARM_OK)
		{
			pxLeg2 = C_result.getDouble();
		}
		else
		{
			ERROR_MSG(C_result.getMsg(),pRet);
			return S_OK;
		}

		// shift seulement sur l'euro !?
		if (ARMLOCAL_ARM_Price(cloneLeg2Id,
							   LocalGetNumObjectId(l_smiledModBump),
							   C_result) == ARM_OK)
		{
			pxLeg2Delta = C_result.getDouble();
		}
		else
		{
			ERROR_MSG(C_result.getMsg(),pRet);
			return S_OK;
		}
	}
	else
	{
		pxLeg2 = 0.0; 
		pxLeg2Delta = 0.0;
	}

	retCode = ARMLOCAL_FreeObject(irIndexLeg2Id,C_result);
	retCode = ARMLOCAL_FreeObject(swapLeg2Id,C_result);
	retCode = ARMLOCAL_FreeObject(cloneLeg2Id,C_result);

	//************************************************************************************************
	//************************** Fin Création Leg Phase2 Refvalue ************************************
	//************************************************************************************************

	double pxTotal; 
	double SensiTotal; 

	double prixInter;
	double Xincr;

//******************** Calcul premiere PV ******************
	int iteration =1;
	double prix;
	double prixAct; 
	double prixAct2;
//	double TauxFixe = 10.01;
	double SprSeek = 10.01;
	double SprSeekFinal= 10.01;
	double SprSeekBonifie = 10.01;  // PUREMENT ARBITRAIRE
	fee = Ntl * fee / 100 ;
//******************** Calcul premiere PV ******************
	VECTOR<double> pxLeg2LivretA = computeLeg2LivretA(startDateLeg1, endDateLA, 
			SprSeekBonifie, l_ccy, 
			dayCountLAId, 
			resetTimingLAId, resetGapLA, resetFreqLAId, payFreqLAId, stubId, amortRefValue,
			l_LAMod, l_LAModBumpInflation, l_LAModBump, vResetMgrIds);

	if (pxLeg2LivretA.size() == 0)
	{
		ERROR_MSG("Pb in computeLeg2LivretA",pRet);
		return S_OK;
	}
	
	pxLA = pxLeg2LivretA[0];
	double pxLADeltaInflation = pxLeg2LivretA[1];
	double pxLADeltaEuribor = pxLeg2LivretA[2];	

	pxTotal = 0.0; 
	pxTotal =  pxLeg2 + pxLA - pxLeg1;
	SensiTotal = fabs((pxLeg2Delta-pxLeg2) + pxLADeltaEuribor - (pxLeg1Delta- pxLeg1)) +fabs(pxLADeltaInflation); 

	prixAct = pxTotal - SensiTotal - fee; 

	prix = prixAct ;
	prixAct2 =1;
	
//********************                    ******************
//******************** BOUCLE SUR LE PRICE ******************
	while ((abs(prixAct2) > 0.0001) && (iteration < 10 ))
	{
		SprSeek = SprSeekBonifie;

		//*************************
		// **** Premier recalcul ****
		//****************************

		pxLeg2LivretA = computeLeg2LivretA(startDateLeg1, endDateLA, 
			SprSeekBonifie, l_ccy, 
			dayCountLAId, 
			resetTimingLAId, resetGapLA, resetFreqLAId, payFreqLAId, stubId, amortRefValue,
			l_LAMod, l_LAModBumpInflation, l_LAModBump, vResetMgrIds);

		if (pxLeg2LivretA.size() == 0)
		{
			ERROR_MSG("Pb in computeLeg2LivretA",pRet);
			return S_OK;
		}

		pxLA = pxLeg2LivretA[0];
		pxLADeltaInflation = pxLeg2LivretA[1];
		pxLADeltaEuribor = pxLeg2LivretA[2];

		pxTotal =  pxLeg2 + pxLA - pxLeg1;
		SensiTotal = fabs((pxLeg2Delta-pxLeg2) + (pxLADeltaEuribor) - (pxLeg1Delta- pxLeg1)) +fabs(pxLADeltaInflation); 

		prixAct = pxTotal - SensiTotal - fee; 

		// les deltas se compensent et on somme en valeur absolue les deltas et vegas

		//******************************
		// **** affectation valeurs ****
		//******************************
		prixAct2 = prixAct ;
		Xincr = SprSeek / 1000000;
		SprSeekBonifie = SprSeek + Xincr;		

		//****************************
		// **** Deuxieme recalcul ****
		//****************************		
		pxLeg2LivretA = computeLeg2LivretA(startDateLeg1, endDateLA, 
			SprSeekBonifie, l_ccy, 
			dayCountLAId, 
			resetTimingLAId, resetGapLA, resetFreqLAId, payFreqLAId, stubId, amortRefValue,
			l_LAMod, l_LAModBumpInflation, l_LAModBump, vResetMgrIds);

		if (pxLeg2LivretA.size() == 0)
		{
			ERROR_MSG("Pb in computeLeg2LivretA",pRet);
			return S_OK;
		}

		pxLA = pxLeg2LivretA[0];
		pxLADeltaInflation = pxLeg2LivretA[1];
		pxLADeltaEuribor = pxLeg2LivretA[2];

		pxTotal =  pxLeg2 + pxLA - pxLeg1;

		SensiTotal = fabs((pxLeg2Delta-pxLeg2) + (pxLADeltaEuribor) -  (pxLeg1Delta- pxLeg1)) +fabs(pxLADeltaInflation); 

		prixAct = pxTotal - SensiTotal - fee; 

		iteration++;

		SprSeekFinal = SprSeekBonifie;
		prixInter = ( prixAct - prixAct2 ) / Xincr;
		// if prixinter =  0 sortir erreur
		SprSeekBonifie = SprSeek - prixAct / prixInter ; 
		prixAct2 = prixAct;	
	}

//**************** PRICE *******************
	retCode = ARMLOCAL_FreeObject(amortRefValueToDestroy,C_result);
	retCode = ARMLOCAL_FreeObject(ccyId,C_result);

	retCode = ARMLOCAL_FreeObject(irIndexLAId,C_result);

	VECTOR<double> vRes;

	vRes.push_back(SprSeekFinal);

	Res = VECTORDOUBLE2VARIANT(vRes, pRet);

	_variant_t wrap_Res;
	wrap_Res.Attach(*pRet);
	*pRet=wrap_Res.Detach();

	return S_OK;
}
catch(...)
{
	return E_FAIL;
}
}



VECTOR<double> computeTxFixed(long indexFixedId, double startDate, double endDate,
							    double spread, CCString l_ccy, long dayCountId,
							   long resetTimingId, double resetGap, long resetFreqId, long payFreqId, 
							   long stubRuleId, long amortRefValue,
							   const CCString& l_Mod, const CCString& l_ModBump)
{
	VECTOR<double> resultat;
	resultat.clear();

try
{
	long retCode; 
	ARM_result C_result;

	long swapLegId; 

	if (ARMLOCAL_SWAPLEG(indexFixedId,
						 startDate,
						 endDate,
						 K_RCV,
						 0,
						 spread,
						 false,
						 l_ccy,
						 dayCountId,
						 resetGap,
						 "NULL",
						 "NULL",
						 0,
						 K_NX_NONE,
						 stubRuleId,
						 -1.0,
						 K_NO,
								 -1,
						 C_result) == ARM_OK) 
	{
		swapLegId = C_result.getLong();
	}

	// ********** Phase2 Leg ************
	long cloneLegId;
	if (ARMLOCAL_ClonedAndSetNotional(swapLegId,
									  amortRefValue,
									  100.,
									  C_result) == ARM_OK)
	{
		cloneLegId = C_result.getLong();
	}
	
	//*********   Price Phase2 *************************
	double pxLeg;
	double pxLegDelta; 

	if (ARMLOCAL_ARM_Price(cloneLegId,
						   LocalGetNumObjectId(l_Mod),
						   C_result) == ARM_OK)
	{
		pxLeg = C_result.getDouble();
	}

	// shift seulement sur l'euro !?
	if (ARMLOCAL_ARM_Price(cloneLegId,
						   LocalGetNumObjectId(l_ModBump),
						   C_result) == ARM_OK)
	{
		pxLegDelta = C_result.getDouble();
	}

	retCode = ARMLOCAL_FreeObject(swapLegId,C_result);
	retCode = ARMLOCAL_FreeObject(cloneLegId,C_result);
	//************************************************************************************************
	//*************************** Fin création Digitale **********************************************
	//************************************************************************************************

	//*********   Price Première Phase structure 1  TD*************************
	VECTOR<double> vRes;
	
	vRes.push_back(pxLeg);
	vRes.push_back(pxLegDelta - pxLeg); 

	return vRes;
}
catch(...)
{
	return resultat;
}
}


STDMETHODIMP ActiveXModule::ARMcomputeTxFixed(double pNtl, 
										   double pStartDateLeg1,  
										   double pEndDateLeg1,
										   BSTR pCcy,
										   BSTR pIndexLeg1, 										   
										   VARIANT pSpreadLeg1,
										   BSTR pDayCountLeg1,
										   BSTR pPayFreqLeg1,
										   BSTR pResetFreqLeg1,
										   BSTR pResetTimingLeg1,
										   BSTR pAdjLeg1,
										   BSTR pRoll, 
										   BSTR pStub,
										   double pEndDateFixed,
										   double pSpreadLeg2,
										   BSTR pDayCountFixed,
										   BSTR pPayFreqFixed,
										   BSTR pResetFreqFixed,
										   BSTR pResetTimingFixed,									   
										   BSTR pAdjFixed, 
										   BSTR pIndexLeg2, 
										   BSTR pDayCountLeg2, 
										   BSTR pPayFreqLeg2, 
										   BSTR pResetFreqLeg2,
										   BSTR pResetTimingLeg2,
										   BSTR pAdjLeg2, 
										   double pEndDateAmort,
										   BSTR pDayCountAmort,
										   BSTR pIntRuleAmort, 
										   double pTxAmort,
										   BSTR pFreqAmort,
										   double pAmountAmort, 
										   BSTR pTypeAmort, 
										   double pFee, 
										   BSTR pSmiledMod, 
										   BSTR pSmiledModBump, 
										   VARIANT *pRet
										   )

{	
	// TODO: Add your implementation code here
try
{									   
	ARM_result C_result;
	double startDateLeg1 = pStartDateLeg1; 
	double endDateLeg1 = pEndDateLeg1;

	long spreadType; 
	double spreadLeg1; 
	
	if (VARIANT2Double(pSpreadLeg1,&spreadLeg1) == S_FALSE)
	{
		_bstr_t bpSpreadLeg1(pSpreadLeg1);
		CCString l_SpreadLeg1= bpSpreadLeg1;
		spreadLeg1 = (double) LocalGetNumObjectId(l_SpreadLeg1);
		spreadType = 1;
	}
	else
	{
	  spreadType = 0;
	}

	double endDateFixed = pEndDateFixed;
	double spreadLeg2= pSpreadLeg2;

	double Ntl = pNtl;
	double fee = pFee; 
	VECTOR<double> resultat;

	long retCode;
	double Res = 0.;

	_bstr_t ccy(pCcy);
	CCString l_ccy = ccy;

	_bstr_t indexLeg1(pIndexLeg1);
	CCString l_indexLeg1 = indexLeg1;

	_bstr_t dayCountLeg1(pDayCountLeg1);
	CCString l_dayCountLeg1 = dayCountLeg1;

	_bstr_t resetFreqLeg1(pResetFreqLeg1);
	CCString l_resetFreqLeg1 = resetFreqLeg1;

	_bstr_t payFreqLeg1(pPayFreqLeg1);
	CCString l_payFreqLeg1 = payFreqLeg1;

	_bstr_t resetTimingLeg1(pResetTimingLeg1);
	CCString l_resetTimingLeg1 = resetTimingLeg1;

	_bstr_t roll(pRoll);
	CCString l_roll = roll;

	_bstr_t adjLeg1(pAdjLeg1);
	CCString l_adjLeg1 = adjLeg1;
	
	_bstr_t stub(pStub);
	CCString l_stub= stub;

	//Strike is double??

	_bstr_t payFreqFixed(pPayFreqFixed);
	CCString l_payFreqFixed = payFreqFixed;

	_bstr_t resetFreqFixed(pResetFreqFixed);
	CCString l_resetFreqFixed = resetFreqFixed;

	_bstr_t resetTimingFixed(pResetTimingFixed);
	CCString l_resetTimingFixed = resetTimingFixed;

	_bstr_t adjFixed(pAdjFixed);
	CCString l_adjFixed= adjFixed;

	_bstr_t dayCountFixed(pDayCountFixed);
	CCString l_dayCountFixed = dayCountFixed;

	_bstr_t resetTimingLeg2(pResetTimingLeg2);
	CCString l_resetTimingLeg2 = resetTimingLeg2;

	_bstr_t dayCountLeg2(pDayCountLeg2);
	CCString l_dayCountLeg2= dayCountLeg2;

	_bstr_t payFreqLeg2(pPayFreqLeg2);
	CCString l_payFreqLeg2= payFreqLeg2;
	
	_bstr_t resetFreqLeg2(pResetFreqLeg2);
	CCString l_resetFreqLeg2 = resetFreqLeg2;

	_bstr_t adjLeg2(pAdjLeg2);
	CCString l_adjLeg2 = adjLeg2;

	_bstr_t indexLeg2(pIndexLeg2);
	CCString l_indexLeg2 = indexLeg2;

	_bstr_t typeAmort(pTypeAmort);
	CCString l_typeAmort = typeAmort;

	_bstr_t smiledMod (pSmiledMod);
	CCString l_smiledMod = smiledMod;

	_bstr_t smiledModBump(pSmiledModBump); 
	CCString l_smiledModBump = smiledModBump;
	
	long payFreqLeg1Id; 	
	if((payFreqLeg1Id = ARM_ConvFrequency(l_payFreqLeg1, C_result)) == ARM_DEFAULT_ERR)
	{
		ERROR_MSG("Invalid Frq",pRet,ARM_ERROR_FREQ);
		return S_OK;
	}

	long resetFreqLeg1Id;
	if((resetFreqLeg1Id = ARM_ConvFrequency(l_resetFreqLeg1, C_result)) == ARM_DEFAULT_ERR)
	{
		ERROR_MSG("Invalid Frq",pRet,ARM_ERROR_FREQ);
		return S_OK;
	}

	long payFreqLeg2Id;
	if((payFreqLeg2Id = ARM_ConvFrequency(l_payFreqLeg2, C_result)) == ARM_DEFAULT_ERR)
	{
		ERROR_MSG("Invalid Frq",pRet,ARM_ERROR_FREQ);
		return S_OK;
	}

	long resetFreqLeg2Id;
	if((resetFreqLeg2Id = ARM_ConvFrequency(l_resetFreqLeg2, C_result)) == ARM_DEFAULT_ERR)
	{
		ERROR_MSG("Invalid Frq",pRet,ARM_ERROR_FREQ);
		return S_OK;
	}

	long payFreqFixedId;
	if((payFreqFixedId = ARM_ConvFrequency(l_payFreqFixed, C_result)) == ARM_DEFAULT_ERR)
	{
		ERROR_MSG("Invalid Frq",pRet,ARM_ERROR_FREQ);
		return S_OK;
	}

	long resetFreqFixedId;
	if((resetFreqFixedId = ARM_ConvFrequency(l_resetFreqFixed, C_result)) == ARM_DEFAULT_ERR)
	{
		ERROR_MSG("Invalid Frq",pRet,ARM_ERROR_FREQ);
		return S_OK;
	}

	long indexLeg1Id; 
	indexLeg1Id = ARM_ConvIrType(l_indexLeg1); 

	long indexLeg2Id; 
	indexLeg2Id = ARM_ConvIrType(l_indexLeg2); 

	long rollId;
	if((rollId = ARM_ConvFwdRule(l_roll, C_result)) == ARM_DEFAULT_ERR)
	{
		ERROR_MSG("Invalid Int Roll",pRet,ARM_ERROR_FREQ);
		return S_OK;
	}	
	
	long adjLeg1Id = ARM_ConvIntRule(l_adjLeg1);
	long adjLeg2Id = ARM_ConvIntRule(l_adjLeg2);
	long adjFixedId = ARM_ConvIntRule(l_adjFixed);

	long resetTimingLeg1Id = ARM_ConvPayResetRule (l_resetTimingLeg1);
	long resetTimingLeg2Id = ARM_ConvPayResetRule (l_resetTimingLeg2);
	long resetTimingFixedId = ARM_ConvPayResetRule (l_resetTimingFixed);

	long dayCountLeg1Id= ARM_ConvDayCount(l_dayCountLeg1);
	long dayCountLeg2Id= ARM_ConvDayCount(l_dayCountLeg2);
	long dayCountFixedId= ARM_ConvDayCount(l_dayCountFixed);

	long stubId = ARM_ConvStubRule(l_stub);
	long ccyId; 

	if (ARMLOCAL_ISOCCY(l_ccy,C_result) == ARM_OK)
	{
		ccyId = C_result.getLong();
	}
	else
	{
		ERROR_MSG(C_result.getMsg(),pRet);
		return S_OK;
	}

	//************************************************************************************************
	//*********************************** Création Amort Refvalue ************************************
	//************************************************************************************************
	long amortRefValue;
	long amortRefValueToDestroy = ARM_NULL_OBJECT;

	if (LocalGetNumObjectId(l_typeAmort) == ARM_KO)
	{
		double endDateAmort = pEndDateAmort;
		double txAmort = pTxAmort;
		double amountAmort = pAmountAmort;

		_bstr_t intRuleAmort(pIntRuleAmort);
		CCString l_intRuleAmort =  intRuleAmort;

		_bstr_t freqAmort(pFreqAmort);
		CCString l_freqAmort = freqAmort;

		_bstr_t dayCountAmort(pDayCountAmort);
		CCString l_dayCountAmort = dayCountAmort;

		long freqAmortId; 
		if((freqAmortId = ARM_ConvFrequency(l_freqAmort, C_result)) == ARM_DEFAULT_ERR)
		{
			ERROR_MSG("Invalid Frq",pRet,ARM_ERROR_FREQ);
			return S_OK;
		}

		long intRuleAmortId= ARM_ConvIntRule(l_intRuleAmort);
		long dayCountAmortId = ARM_ConvDayCount(l_dayCountAmort);	

		// AmortMethodId
		long amortMethodId;

		if((amortMethodId = ARM_ConvAmortMethod (l_typeAmort, C_result)) == ARM_DEFAULT_ERR)
		{
			ERROR_MSG("Invalid TypeAmort",pRet,ARM_ERROR_FREQ);
			return S_OK;
		}

		long FixedLegAmortId;
		if (ARMLOCAL_FIXEDLEG(startDateLeg1,
							  endDateAmort,
							  K_RCV,
							  0,
							  1.0,
							  dayCountAmortId,
							  freqAmortId,
							  K_COMP_PROP,
							  K_ARREARS,
							  intRuleAmortId,
							  stubId,
							  false,
							  l_ccy,
							  "NULL",
							  K_NX_NONE,
							  -1.0,
							  K_NO,
								 -1,
							  C_result) == ARM_OK) 
							  
		{
			FixedLegAmortId = C_result.getLong();
		}
		else
		{
			ERROR_MSG(C_result.getMsg(),pRet);
			return S_OK;
		}

		if (ARMLOCAL_ARM_GenAmortization(FixedLegAmortId,
										 amortMethodId,
										 freqAmortId,
										 amountAmort,
										 -1,
										 Ntl,
										 txAmort,
										 0.,
										 LocalGetNumObjectId(l_smiledMod),
										 0.0,
										 C_result) == ARM_OK) 
		{
			amortRefValue = C_result.getLong();
			amortRefValueToDestroy = amortRefValue;
		}
		else
		{
			ERROR_MSG(C_result.getMsg(),pRet);
			return S_OK;
		}

		retCode = ARMLOCAL_FreeObject(FixedLegAmortId,C_result);
	}
	else
	{
		amortRefValue = LocalGetNumObjectId(l_typeAmort);

		// calcul du 1er nominal
		double dPaymentDate;

		if (ARMLOCAL_DisplayRefValue(amortRefValue,true,C_result) == ARM_OK)
		{
			dPaymentDate = C_result.getArray(0);
		}
		else
		{
			ERROR_MSG(C_result.getMsg(),pRet);
			return S_OK;
		}
		
		if (ARMLOCAL_CptRefValue(amortRefValue,dPaymentDate,C_result) == ARM_OK)
		{
			Ntl = C_result.getDouble();
		}
		else
		{
			ERROR_MSG(C_result.getMsg(),pRet);
			return S_OK;
		} 
	}
	//************************************************************************************************
	//*********************************** Fin Création Amort Refvalue ************************************
	//************************************************************************************************
	
	//************************************************************************************************
	//*********************************** Création Leg Phase1 Refvalue *******************************
	//************************************************************************************************

	long irIndexLeg1Id;
	double pxLeg1; 
	double pxLeg1Delta;
	long swapLeg1Id;
	long cloneLeg1Id;
	double resetGapLeg1; 

	if (resetTimingLeg1Id == 1)
	{
		resetGapLeg1 = -2.; 
	}
	else
	{
		resetGapLeg1 = -15.; 
	}

	if (ARMLOCAL_IRINDEX(KACTUAL_360,
						 payFreqLeg1Id,
						 -1,
						 K_COMP_PROP,
						 K_MOD_FOLLOWING,
						 resetTimingLeg1Id,	
						 resetGapLeg1,
						 K_ARREARS,
						 10000.,
						 false,
						 l_ccy,
						 indexLeg1Id,
						 K_COMP_PROP,
						 adjLeg1Id,
						 resetFreqLeg1Id,
						 C_result) == ARM_OK) 
	{
		irIndexLeg1Id = C_result.getLong();
	}
	else
	{
		ERROR_MSG(C_result.getMsg(),pRet);
		return S_OK;
	}
	
	if (ARMLOCAL_SWAPLEG(irIndexLeg1Id,
						 startDateLeg1,
						 endDateLeg1,
						 K_RCV,
						 spreadType, 
						 spreadLeg1, 
						 false,
						 l_ccy,
						 dayCountLeg1Id,
						 resetGapLeg1,
						 "NULL",
						 "NULL",
						 0,
						 K_NX_NONE,
						 stubId,
						 -1.0,
						 K_NO,
								 -1,
						 C_result) == ARM_OK) 
	{
		swapLeg1Id = C_result.getLong();
	}
	else
	{
		ERROR_MSG(C_result.getMsg(),pRet);
		return S_OK;
	}

	// ********** Phase1 Leg ************	
	if (ARMLOCAL_ClonedAndSetNotional(swapLeg1Id,
									  amortRefValue,
									  100.,
									  C_result) == ARM_OK)
	{
		cloneLeg1Id = C_result.getLong();
	}
	else
	{
		ERROR_MSG(C_result.getMsg(),pRet);
		return S_OK;
	}

	//*********   Price Phase1 *************************
	if (ARMLOCAL_ARM_Price(cloneLeg1Id,
						   LocalGetNumObjectId(l_smiledMod),
						   C_result) == ARM_OK)
	{
		pxLeg1 = C_result.getDouble();
	}
	else
	{
		ERROR_MSG(C_result.getMsg(),pRet);
		return S_OK;
	}

	// shift seulement sur l'euro !?
	if (ARMLOCAL_ARM_Price(cloneLeg1Id,
						   LocalGetNumObjectId(l_smiledModBump),
						   C_result) == ARM_OK)
	{
		pxLeg1Delta = C_result.getDouble();
	}
	else
	{
		ERROR_MSG(C_result.getMsg(),pRet);
		return S_OK;
	}

	retCode = ARMLOCAL_FreeObject(irIndexLeg1Id,C_result);
	retCode = ARMLOCAL_FreeObject(swapLeg1Id,C_result);
	retCode = ARMLOCAL_FreeObject(cloneLeg1Id,C_result);

	//************************************************************************************************
	//************************** Fin Création Leg Phase1 Refvalue ************************************
	//************************************************************************************************

	//************************************************************************************************
	//*********************************** Création Leg Fixed Refvalue *******************************
	//************************************************************************************************
	double resetGapFixed; 

	if (resetTimingFixedId == 1)
	{
		resetGapFixed = -2.; 
	}
	else
	{
		resetGapFixed = -15.; 
	}

	long irIndexIdFixed;
	if (ARMLOCAL_IRINDEX(KACTUAL_360,
						 payFreqFixedId,
						 -1.0,
						 K_COMP_PROP,
						 K_MOD_FOLLOWING,
						 resetTimingFixedId,
						 resetGapFixed,
						 K_ARREARS,
						 10000.,
						 false,
						 l_ccy,
						 K_FIXED,
						 K_COMP_PROP,
						 adjFixedId,
						 resetFreqFixedId,
						 C_result) == ARM_OK) 
	{
		irIndexIdFixed = C_result.getLong();
	}
	else
	{
		ERROR_MSG(C_result.getMsg(),pRet);
		return S_OK;
	}
	//************************************************************************************************
	//************************** Fin Création Leg Phase1 Refvalue ************************************
	//************************************************************************************************

	//************************************************************************************************
	//*********************************** Création Leg2 Phase2 Refvalue *******************************
	//************************************************************************************************
	long irIndexLeg2Id;
	double resetGapLeg2;
	long swapLeg2Id; 
	long cloneLeg2Id; 

	double pxLeg2Delta; 
	double pxLeg2; 
	if (resetTimingLeg2Id == 1)
	{
		resetGapLeg2 = -2.; 
	}
	else
	{
		resetGapLeg2 = -15.; 
	}

	if ( endDateFixed != endDateLeg1 )
	{

		if (ARMLOCAL_IRINDEX(KACTUAL_360,
							 payFreqLeg2Id,
							 -1,
							 K_COMP_PROP,
							 K_MOD_FOLLOWING,
							 resetTimingLeg2Id,	
							 resetGapLeg2,
							 K_ARREARS,
							 10000.,
							 false,
							 l_ccy,
							 indexLeg2Id,
							 K_COMP_PROP,
							 adjLeg2Id,
							 resetFreqLeg2Id,
							 C_result) == ARM_OK) 
		{
			irIndexLeg2Id = C_result.getLong();
		}
		else
		{
			ERROR_MSG(C_result.getMsg(),pRet);
			return S_OK;
		}

		if (ARMLOCAL_SWAPLEG(irIndexLeg2Id,
							endDateFixed, 
							endDateLeg1, 						
							 K_RCV,
							 0, 
							 spreadLeg2, 
							 false,
							 l_ccy,
							 dayCountLeg2Id,
							 resetGapLeg2,
							 "NULL",
							 "NULL",
							 0,
							 K_NX_NONE,
							 stubId,
							 -1.0,
							 K_NO,
								 -1,
							 C_result) == ARM_OK) 
		{
			swapLeg2Id = C_result.getLong();
		}
		else
		{
			ERROR_MSG(C_result.getMsg(),pRet);
			return S_OK;
		}

		// ********** Phase2 Leg ************	
		if (ARMLOCAL_ClonedAndSetNotional(swapLeg2Id,
										  amortRefValue,
										  100.,
										  C_result) == ARM_OK)
		{
			cloneLeg2Id = C_result.getLong();
		}
		else
		{
			ERROR_MSG(C_result.getMsg(),pRet);
			return S_OK;
		}

		//*********   Price Phase2 *************************
		if (ARMLOCAL_ARM_Price(cloneLeg2Id,
							   LocalGetNumObjectId(l_smiledMod),
							   C_result) == ARM_OK)
		{
			pxLeg2 = C_result.getDouble();
		}
		else
		{
			ERROR_MSG(C_result.getMsg(),pRet);
			return S_OK;
		}

		// shift seulement sur l'euro !?
		if (ARMLOCAL_ARM_Price(cloneLeg2Id,
							   LocalGetNumObjectId(l_smiledModBump),
							   C_result) == ARM_OK)
		{
			pxLeg2Delta = C_result.getDouble();
		}
		else
		{
			ERROR_MSG(C_result.getMsg(),pRet);
			return S_OK;
		}
	}
	else
	{
		pxLeg2 = 0.0;
		pxLeg2Delta = 0.0;
	}

	retCode = ARMLOCAL_FreeObject(irIndexLeg2Id,C_result);
	retCode = ARMLOCAL_FreeObject(swapLeg2Id,C_result);
	retCode = ARMLOCAL_FreeObject(cloneLeg2Id,C_result);

	//************************************************************************************************
	//************************** Fin Création Leg Phase2 Refvalue ************************************
	//************************************************************************************************
	double pxTotal; 
	double SensiTotal; 

	double prixInter;
	double Xincr;

//******************** Calcul premiere PV ******************
	int iteration =1;
	double prix;
	double prixAct; 
	double prixAct2;
	double SprSeek = 10.01;
	double SprSeekFinal= 10.01;
	double SprSeekBonifie = 10.01;  // PUREMENT ARBITRAIRE
	fee = Ntl * fee / 100 ;
//******************** Calcul premiere PV ******************
	VECTOR<double> pxLegFixed = computeTxFixed(irIndexIdFixed, startDateLeg1, endDateFixed, 
			SprSeekBonifie, l_ccy, 
			dayCountFixedId, 
			resetTimingFixedId, resetGapFixed, resetFreqFixedId, payFreqFixedId, stubId, amortRefValue,
			l_smiledMod, l_smiledModBump);


	if (pxLegFixed.size() == 0)
	{
		ERROR_MSG("Pb in computeTxFixed",pRet);
		return S_OK;
	}
	
	double pxFixed = pxLegFixed[0];
	double pxFixedDelta = pxLegFixed[1];

	pxTotal = 0.0; 
	pxTotal =  pxLeg2 + pxFixed - pxLeg1;
	SensiTotal = fabs((pxLeg2Delta-pxLeg2) + pxFixedDelta - (pxLeg1Delta- pxLeg1)); 

	prixAct = pxTotal - SensiTotal - fee; 

	prix = prixAct ;
	prixAct2 =1;
	
//********************                    ******************
//******************** BOUCLE SUR LE PRICE ******************
	while ((abs(prixAct2) > 0.0001) && (iteration < 10 ))
	{
		SprSeek = SprSeekBonifie;
		
		//*************************
		// **** Premier recalcul ****
		//****************************
		pxLegFixed = computeTxFixed(irIndexIdFixed, startDateLeg1, endDateFixed, 
			SprSeekBonifie, l_ccy, 
			dayCountFixedId, 
			resetTimingFixedId, resetGapFixed, resetFreqFixedId, payFreqFixedId, stubId, amortRefValue,
			l_smiledMod, l_smiledModBump);

		if (pxLegFixed.size() == 0)
		{
			ERROR_MSG("Pb in computeTxFixed",pRet);
			return S_OK;
		}

		pxFixed = pxLegFixed[0];
		pxFixedDelta = pxLegFixed[1];

		pxTotal =  pxLeg2 + pxFixed - pxLeg1;
		SensiTotal = fabs((pxLeg2Delta-pxLeg2) + pxFixedDelta - (pxLeg1Delta- pxLeg1)); 
 
		prixAct = pxTotal - SensiTotal - fee; 

		// les deltas se compensent et on somme en valeur absolue les deltas et vegas

		//******************************
		// **** affectation valeurs ****
		//******************************
		prixAct2 = prixAct ;
		Xincr = SprSeek / 1000000;
		SprSeekBonifie = SprSeek + Xincr;		
		
		//****************************
		// **** Deuxieme recalcul ****
		//****************************
		pxLegFixed = computeTxFixed(irIndexIdFixed, startDateLeg1, endDateFixed, 
			SprSeekBonifie, l_ccy, 
			dayCountFixedId, 
			resetTimingFixedId, resetGapFixed, resetFreqFixedId, payFreqFixedId, stubId, amortRefValue,
			l_smiledMod, l_smiledModBump);

		if (pxLegFixed.size() == 0)
		{
			ERROR_MSG("Pb in computeTxFixed",pRet);
			return S_OK;
		}

		pxFixed = pxLegFixed[0];
		pxFixedDelta = pxLegFixed[1];

		pxTotal =  pxLeg2 + pxFixed - pxLeg1;

		SensiTotal = fabs((pxLeg2Delta-pxLeg2) + pxFixedDelta - (pxLeg1Delta- pxLeg1)); 
		prixAct = pxTotal - SensiTotal - fee; 

		iteration++;

		SprSeekFinal = SprSeekBonifie;
		prixInter = ( prixAct - prixAct2 ) / Xincr;
		// if prixinter =  0 sortir erreur
		SprSeekBonifie = SprSeek - prixAct / prixInter ; 
		prixAct2 = prixAct;			
	}

//**************** PRICE *******************
	retCode = ARMLOCAL_FreeObject(amortRefValueToDestroy,C_result);
	retCode = ARMLOCAL_FreeObject(ccyId,C_result);
	retCode = ARMLOCAL_FreeObject(irIndexIdFixed,C_result);

	VECTOR<double> vRes;

	vRes.push_back(SprSeekFinal);

	Res = VECTORDOUBLE2VARIANT(vRes, pRet);

	_variant_t wrap_Res;
	wrap_Res.Attach(*pRet);
	*pRet=wrap_Res.Detach();

	return S_OK;
}
catch(...)
{
	return E_FAIL;
}
}




VECTOR<double> computeCRAFix(double rate,
							double asofDate,
							double startDate,
							double endDate,
							CCString ccy,
							long levelUpId,
							long upSpecId,
							long levelDownId,
							long downSpecId,
							long indexRefId,
							long daycountId,
							long porrId,
							long payIndexId,
							long payFreqPayIndexId,
							long resetFreqRefIndexId,
							long paidRstTimingId,
							long refRstTimingId,
							long stubId,
							long decompPricingFlagId,
							double discMarginFactor,
							double startcallDate,
							const CCString& XStyle,
							double KStyleId,
							double liborlegId,
							double meanreversion,
							long preinitflagId,
							VECTOR<double> vCalibParams,
							VECTOR<double> vCalibParamsPF,
							VECTOR<double> vMarkovTreeParams,
							long markovtreepathnumber,
							double kerneltogp,
							bool computevega,
							const CCString& bsmod,
							const CCString& bsmodswopt,
							const CCString& bsmodswoptbump,
							const CCString& zc)

{
	VECTOR<double> resultat;
	resultat.clear();

try
{
	ARM_result C_result;
	double pxCorridorLeg;
	double discMargin;
	double pxOption;
	double vega;
	double pxOptionBump;

	long retCode;
	long res;

	VECTOR <long> instrumentIdVector;
	VECTOR <double> coefVector;
	VECTOR <double> marketpriceVector;
	VECTOR<double> precisionVector;
	VECTOR <double> initsigmaVector;
	VECTOR <double> initbetaVector;

	//********************CORRIDOR LEG***********************
	long LDPricingMethodId;

	LDPricingMethodId = ARM_ConvPricCorridorLD("DIG");

	long corridorLegId;

  	if (ARMLOCAL_CORRIDORLEG(startDate,
							 endDate,
							 porrId,
							 payIndexId,
							 payFreqPayIndexId,
							 0,
							 rate,
							 indexRefId,
							 resetFreqRefIndexId,
							 paidRstTimingId,
							 refRstTimingId,
							 stubId,
							 levelDownId,
							 downSpecId,
							 levelUpId,
							 upSpecId,
							 false,
							 ccy,
							 -1,
							 K_LINEAR,
							 LDPricingMethodId,
							 decompPricingFlagId,
							 "NULL", // reset calendar
							 "NULL", // pay calendar
                              0,
							  C_result) == ARM_OK) 
	{
		corridorLegId = C_result.getLong();
	}


	
	if (ARMLOCAL_ARM_Price(corridorLegId,
						   LocalGetNumObjectId(bsmod),
						   C_result) == ARM_OK)
	{
		pxCorridorLeg = C_result.getDouble();
		pxCorridorLeg = pxCorridorLeg*100;
	}
	//******************************************************

	#ifdef _DEBUG
		res = ARMLOCAL_ARM_View(payIndexId,C_result);
	#endif

	#ifdef _DEBUG
		res = ARMLOCAL_ARM_View(indexRefId,C_result);
	#endif

	#ifdef _DEBUG
		res = ARMLOCAL_ARM_View(corridorLegId,C_result);
	#endif

	//********************CORRIDOR LEG FOR CRA**************
	long corridorforCRALegId;

  	if (ARMLOCAL_CORRIDORLEG(startcallDate,
							 endDate,
							 porrId,
							 payIndexId,
							 payFreqPayIndexId,
							 0,
							 rate,
							 indexRefId,
							 resetFreqRefIndexId,
							 paidRstTimingId,
							 refRstTimingId,
							 stubId,
							 levelDownId,
							 downSpecId,
							 levelUpId,
							 upSpecId,
							 false,
							 ccy,
							 -1,
							 K_LINEAR,
							 LDPricingMethodId,
							 decompPricingFlagId,
							 "NULL", // reset calendar
							 "NULL", // pay calendar
                              0,
							  C_result) == ARM_OK) 
	{
		corridorforCRALegId = C_result.getLong();
	}
	//******************************************************

  	//*****************FIXED LEG******************************
	long fixId;
	double pxfixleg;

	if (ARMLOCAL_FIXEDLEG(startDate,
						  endDate,
						  K_RCV,
						  0,
						  rate,
						  daycountId,
						  payFreqPayIndexId,
						  0,
						  K_ARREARS,
						  K_ADJUSTED,
						  K_SHORTSTART,
						  false,
						  ccy,
						  "NULL",
						  K_NX_NONE,
						  -1.0,
						  K_YES,
						  -1,
						  C_result) == ARM_OK) 
	{
		fixId = C_result.getLong();
	}
	
	
	if (ARMLOCAL_ARM_Price(fixId,
						   LocalGetNumObjectId(bsmod),
						   C_result) == ARM_OK)
	{
		pxfixleg = C_result.getDouble();
		pxfixleg = pxfixleg * 100;
	}
	//*******************************************************

	#ifdef _DEBUG
		res = ARMLOCAL_ARM_View(fixId,C_result);
	#endif

	#ifdef _DEBUG
		res = ARMLOCAL_ARM_View(LocalGetNumObjectId(bsmod),C_result);
	#endif

	//*****************DISC MARGIN***************************
	double proba;
		
	proba = fabs(pxCorridorLeg / pxfixleg);
	
	discMargin =   pxfixleg * ( 0.01 + 0.05 * max(0,min(proba,1-proba)-0.1));
	discMargin = discMargin/discMarginFactor;
	//*******************************************************

	//*******************SWAP********************************
	long swapId;

	instrumentIdVector.push_back(liborlegId);
	instrumentIdVector.push_back(corridorforCRALegId);


	coefVector.push_back(-1);
	coefVector.push_back(-1); 

	marketpriceVector.push_back(1);
	marketpriceVector.push_back(1);

	if (ARMLOCAL_PF(instrumentIdVector,
						  coefVector,
						  marketpriceVector,
						  precisionVector,
						  C_result) == ARM_OK) 
	{
		swapId = C_result.getLong();
	}
	//*******************************************************


	//*******************OPTION PORTFOLIO********************
	long optionportId;
	
	if (ARMLOCAL_ARM_OPTIONPORTFOLIO(swapId,
						   LocalGetNumObjectId(XStyle),
						  KStyleId,
						  K_CALL,
						  C_result) == ARM_OK) 
	{
		optionportId = C_result.getLong();
	}
	//*******************************************************

	//******************CALIBRATOR SFRM**********************
	long calibratorId;


	if (ARMLOCAL_CALIBRATORSFRM(optionportId,
						LocalGetNumObjectId(bsmod),
						LocalGetNumObjectId(bsmodswopt),
						meanreversion,
						vCalibParams,
						preinitflagId,
						initsigmaVector,
						initbetaVector,
						ARM_NULL_OBJECT,
						ARM_NULL_OBJECT,
						ARM_NULL_OBJECT,
						K_DIAG,
						ARM_NULL_OBJECT,
						vCalibParamsPF,
						"NO",
						 C_result) == ARM_OK) 
	{
		calibratorId = C_result.getLong();
	}
	//*******************************************************

		//******************SFRM CALIBRATE***********************
	long sfrmcalibrateId;

	if (ARMLOCAL_SFRMCALIBRATE (calibratorId,
                             ARM_NULL_OBJECT,
                             ARM_NULL_OBJECT,
                             "YES",
                             "YES",
                             kerneltogp,
							  C_result) == ARM_OK) 
	{
		sfrmcalibrateId = C_result.getLong();
	}
	//********************************************************

		//*****************MARKOVTREE*****************************
	long markovtreeId;

	if (ARMLOCAL_FRMMARKOVTREE (asofDate,
                             endDate,
                             LocalGetNumObjectId(zc),
							 markovtreepathnumber,
							 vMarkovTreeParams,
							sfrmcalibrateId,
							  C_result) == ARM_OK) 
	{
		markovtreeId = C_result.getLong();
	}
	//********************************************************
	
		//******************CALCUL PRIX OPTION********************
	if (ARMLOCAL_ARM_Price(optionportId,
						   markovtreeId,
						   C_result) == ARM_OK)
	{
		pxOption= C_result.getDouble();
		pxOption = pxOption * 100;
	}
	//********************************************************

	#ifdef _DEBUG
		res = ARMLOCAL_ARM_View(LocalGetNumObjectId(zc),C_result);
	#endif

	#ifdef _DEBUG
		res = ARMLOCAL_ARM_View(LocalGetNumObjectId(bsmod),C_result);
	#endif

	#ifdef _DEBUG
		res = ARMLOCAL_ARM_View(LocalGetNumObjectId(bsmodswopt),C_result);
	#endif

	#ifdef _DEBUG
		res = ARMLOCAL_ARM_View(optionportId,C_result);
	#endif

	#ifdef _DEBUG
		res = ARMLOCAL_ARM_View(calibratorId,C_result);
	#endif

	#ifdef _DEBUG
		res = ARMLOCAL_ARM_View(markovtreeId,C_result);
	#endif

		long calibratorbumpId;
		long sfrmcalibratebumpId;
		long markovtreebumpId;
	if (computevega==true)
	{

		//******************CALIBRATOR SFRM BUMPED**********************
		
		


		if (ARMLOCAL_CALIBRATORSFRM(optionportId,
							LocalGetNumObjectId(bsmod),
							LocalGetNumObjectId(bsmodswoptbump),
							meanreversion,
							vCalibParams,
							preinitflagId,
							initsigmaVector,
							initbetaVector,
							ARM_NULL_OBJECT,
							ARM_NULL_OBJECT,
							ARM_NULL_OBJECT,
							K_DIAG,
							ARM_NULL_OBJECT,
							vCalibParamsPF,
							"NO",
							 C_result) == ARM_OK) 
		{
			calibratorbumpId = C_result.getLong();
		}
		
		//*******************************************************



		//******************SFRM CALIBRATE BUMPED*****************
		

		if (ARMLOCAL_SFRMCALIBRATE (calibratorbumpId,
								 ARM_NULL_OBJECT,
								 ARM_NULL_OBJECT,
								 "YES",
								 "YES",
								 kerneltogp,
								  C_result) == ARM_OK) 
		{
			sfrmcalibratebumpId = C_result.getLong();
		}
		//********************************************************




	  //*****************MARKOVTREE BUMPED*****************************
		

		if (ARMLOCAL_FRMMARKOVTREE (asofDate,
								 endDate,
								 LocalGetNumObjectId(zc),
								 markovtreepathnumber,
								 vMarkovTreeParams,
								sfrmcalibratebumpId,
								  C_result) == ARM_OK) 
		{
			markovtreebumpId = C_result.getLong();
		}
		//********************************************************



		//******************CALCUL VEGA********************
		if (ARMLOCAL_ARM_Price(optionportId,
							   markovtreebumpId,
							   C_result) == ARM_OK)
		{
			pxOptionBump= C_result.getDouble();
			pxOptionBump = pxOptionBump * 100;
		}

		vega = fabs(pxOptionBump-pxOption);
		//********************************************************
	}

	VECTOR<double> vRes;
	
	vRes.push_back(pxCorridorLeg);
	vRes.push_back(discMargin);
	vRes.push_back(pxOption);
	if (computevega==true)
	{
		vRes.push_back(vega);
	}

	retCode = ARMLOCAL_FreeObject(corridorLegId,C_result);
	retCode = ARMLOCAL_FreeObject(corridorforCRALegId,C_result);
	retCode = ARMLOCAL_FreeObject(fixId,C_result);
	retCode = ARMLOCAL_FreeObject(swapId,C_result);
	retCode = ARMLOCAL_FreeObject(optionportId,C_result);
	retCode = ARMLOCAL_FreeObject(calibratorId,C_result);
	if (computevega==true)
	{
		retCode = ARMLOCAL_FreeObject(calibratorbumpId,C_result);
		retCode = ARMLOCAL_FreeObject(sfrmcalibratebumpId,C_result);
		retCode = ARMLOCAL_FreeObject(markovtreebumpId,C_result);
	}
	retCode = ARMLOCAL_FreeObject(sfrmcalibrateId,C_result);
	retCode = ARMLOCAL_FreeObject(markovtreeId,C_result);
	
	return vRes;
}
catch(...)
{
	return resultat;
}
}


VECTOR<double> computeCRAFloat(double rate,
							double asofDate,
							double startDate,
							double endDate,
							CCString ccy,
							long levelUpId,
							long upSpecId,
							long levelDownId,
							long downSpecId,
							long indexRefId,
							long daycountId,
							long porrId,
							long payIndexId,
							long payIndexFixId,
							long payFreqPayIndexId,
							long resetFreqRefIndexId,
							long paidRstTimingId,
							long refRstTimingId,
							long stubId,
							long decompPricingFlagId,
							double discMarginFactor,
							double startcallDate,
							const CCString& XStyle,
							double KStyleId,
							double liborlegId,
							double meanreversion,
							long preinitflagId,
							VECTOR<double> vCalibParams,
							VECTOR<double> vCalibParamsPF,
							VECTOR<double> vMarkovTreeParams,
							long markovtreepathnumber,
							double kerneltogp,
							double pxFloatingLeg,
							double pxRaPaysFloatIndexLeg,
							bool computevega,
							const CCString& bsmod,
							const CCString& bsmodswopt,
							const CCString& bsmodswoptbump,
							const CCString& zc)

{
	VECTOR<double> resultat;
	resultat.clear();

try
{
	ARM_result C_result;
	double pxCorridorLeg;
	double discMargin;
	double pxOption;
	double vega;
	double pxOptionBump;

	long retCode;
	long res;

	VECTOR<double> vRes;

	VECTOR <long> instrumentIdVector;
	VECTOR <double> coefVector;
	VECTOR <double> marketpriceVector;
	VECTOR<double> precisionVector;
	VECTOR <double> initsigmaVector;
	VECTOR <double> initbetaVector;

	//*****************FIXED LEG******************************
	long fixId;
	double pxfixleg;

	if (ARMLOCAL_FIXEDLEG(startDate,
						  endDate,
						  K_RCV,
						  0,
						  rate,
						  daycountId,
						  payFreqPayIndexId,
						  0,
						  K_ARREARS,
						  K_ADJUSTED,
						  K_SHORTSTART,
						  false,
						  ccy,
						  "NULL",
						  K_NX_NONE,
						  -1.0,
						  K_YES,
						  -1,
						  C_result) == ARM_OK) 
	{
		fixId = C_result.getLong();
	}
	
	
	if (ARMLOCAL_ARM_Price(fixId,
						   LocalGetNumObjectId(bsmod),
						   C_result) == ARM_OK)
	{
		pxfixleg = C_result.getDouble();
	}
	//*******************************************************

	#ifdef _DEBUG
		res = ARMLOCAL_ARM_View(fixId,C_result);
	#endif

	#ifdef _DEBUG
		res = ARMLOCAL_ARM_View(LocalGetNumObjectId(bsmod),C_result);
	#endif

	//********************RA PAY SPREAD ON INDEX****************
	long LDPricingMethodId;

	LDPricingMethodId = ARM_ConvPricCorridorLD("DIG");

	long rapayspreadonindexId;

	if (ARMLOCAL_CORRIDORLEG(startDate,
							 endDate,
							 porrId,
							 payIndexFixId,
							 payFreqPayIndexId,
							 0,
							 rate,
							 indexRefId,
							 resetFreqRefIndexId,
							 paidRstTimingId,
							 refRstTimingId,
							 stubId,
							 levelDownId,
							 downSpecId,
							 levelUpId,
							 upSpecId,
							 false,
							 ccy,
							 -1,
							 K_LINEAR,
							 LDPricingMethodId,
							 decompPricingFlagId,
							 "NULL", // reset calendar
							 "NULL", // pay calendar
                              0,
							  C_result) == ARM_OK) 
	{
		rapayspreadonindexId = C_result.getLong();
	}

	double pxrapayspreadonindexLeg;
	
	if (ARMLOCAL_ARM_Price(rapayspreadonindexId,
						   LocalGetNumObjectId(bsmod),
						   C_result) == ARM_OK)
	{
		pxrapayspreadonindexLeg = C_result.getDouble();
	}
	//******************************************************

	#ifdef _DEBUG
		res = ARMLOCAL_ARM_View(payIndexFixId,C_result);
	#endif

	#ifdef _DEBUG
		res = ARMLOCAL_ARM_View(indexRefId,C_result);
	#endif

	#ifdef _DEBUG
		res = ARMLOCAL_ARM_View(rapayspreadonindexId,C_result);
	#endif

	
	discMargin= pxFloatingLeg * (0.01+0.05*max(min(pxRaPaysFloatIndexLeg/pxFloatingLeg,1-pxRaPaysFloatIndexLeg/pxFloatingLeg)-0.1,0));
	discMargin= discMargin + pxfixleg * (0.01+0.05*max(min(pxrapayspreadonindexLeg/pxfixleg,1-pxrapayspreadonindexLeg/pxfixleg)-0.1,0));
	discMargin=100 * discMargin;
	discMargin= discMargin / discMarginFactor;


	//********************CORRIDOR LEG***********************
	long corridorLegId;

  	if (ARMLOCAL_CORRIDORLEG(startDate,
							 endDate,
							 porrId,
							 payIndexId,
							 payFreqPayIndexId,
							 0,
							 rate,
							 indexRefId,
							 resetFreqRefIndexId,
							 paidRstTimingId,
							 refRstTimingId,
							 stubId,
							 levelDownId,
							 downSpecId,
							 levelUpId,
							 upSpecId,
							 false,
							 ccy,
							 -1,
							 K_LINEAR,
							 LDPricingMethodId,
							 decompPricingFlagId,
							 "NULL", // reset calendar
							 "NULL", // pay calendar
                              0,
							  C_result) == ARM_OK) 
	{
		corridorLegId = C_result.getLong();
	}


	
	if (ARMLOCAL_ARM_Price(corridorLegId,
						   LocalGetNumObjectId(bsmod),
						   C_result) == ARM_OK)
	{
		pxCorridorLeg = C_result.getDouble();
		pxCorridorLeg = pxCorridorLeg*100;
	}
	//******************************************************

	//********************CORRIDOR LEG FOR CRA**************
	long corridorforCRALegId;

  	if (ARMLOCAL_CORRIDORLEG(startcallDate,
							 endDate,
							 porrId,
							 payIndexId,
							 payFreqPayIndexId,
							 0,
							 rate,
							 indexRefId,
							 resetFreqRefIndexId,
							 paidRstTimingId,
							 refRstTimingId,
							 stubId,
							 levelDownId,
							 downSpecId,
							 levelUpId,
							 upSpecId,
							 false,
							 ccy,
							 -1,
							 K_LINEAR,
							 LDPricingMethodId,
							 decompPricingFlagId,
							 "NULL", // reset calendar
							 "NULL", // pay calendar
                              0,
							  C_result) == ARM_OK) 
	{
		corridorforCRALegId = C_result.getLong();
	}
	//******************************************************


	//*******************SWAP********************************
	long swapId;

	instrumentIdVector.push_back(liborlegId);
	instrumentIdVector.push_back(corridorforCRALegId);


	coefVector.push_back(-1);
	coefVector.push_back(-1); 

	marketpriceVector.push_back(1);
	marketpriceVector.push_back(1);

	if (ARMLOCAL_PF(instrumentIdVector,
						  coefVector,
						  marketpriceVector,
						  precisionVector,
						  C_result) == ARM_OK) 
	{
		swapId = C_result.getLong();
	}
	//*******************************************************


	//*******************OPTION PORTFOLIO********************
	long optionportId;
	
	if (ARMLOCAL_ARM_OPTIONPORTFOLIO(swapId,
						   LocalGetNumObjectId(XStyle),
						  KStyleId,
						  K_CALL,
						  C_result) == ARM_OK) 
	{
		optionportId = C_result.getLong();
	}
	//*******************************************************

	//******************CALIBRATOR SFRM**********************
	long calibratorId;


	if (ARMLOCAL_CALIBRATORSFRM(optionportId,
						LocalGetNumObjectId(bsmod),
						LocalGetNumObjectId(bsmodswopt),
						meanreversion,
						vCalibParams,
						preinitflagId,
						initsigmaVector,
						initbetaVector,
						ARM_NULL_OBJECT,
						ARM_NULL_OBJECT,
						ARM_NULL_OBJECT,
						K_DIAG,
						ARM_NULL_OBJECT,
						vCalibParamsPF,
						"YES",
						 C_result) == ARM_OK) 
	{
		calibratorId = C_result.getLong();
	}
	//*******************************************************

		//******************SFRM CALIBRATE***********************
	long sfrmcalibrateId;

	if (ARMLOCAL_SFRMCALIBRATE (calibratorId,
                             ARM_NULL_OBJECT,
                             ARM_NULL_OBJECT,
                             "YES",
                             "YES",
                             kerneltogp,
							  C_result) == ARM_OK) 
	{
		sfrmcalibrateId = C_result.getLong();
	}
	//********************************************************

		//*****************MARKOVTREE*****************************
	long markovtreeId;

	if (ARMLOCAL_FRMMARKOVTREE (asofDate,
                             endDate,
                             LocalGetNumObjectId(zc),
							 markovtreepathnumber,
							 vMarkovTreeParams,
							sfrmcalibrateId,
							  C_result) == ARM_OK) 
	{
		markovtreeId = C_result.getLong();
	}
	//********************************************************
	
	//******************CALCUL PRIX OPTION********************
	if (ARMLOCAL_ARM_Price(optionportId,
						   markovtreeId,
						   C_result) == ARM_OK)
	{
		pxOption= C_result.getDouble();
		pxOption = pxOption * 100;
	}
	//********************************************************

	#ifdef _DEBUG
		res = ARMLOCAL_ARM_View(LocalGetNumObjectId(zc),C_result);
	#endif

	#ifdef _DEBUG
		res = ARMLOCAL_ARM_View(LocalGetNumObjectId(bsmod),C_result);
	#endif

	#ifdef _DEBUG
		res = ARMLOCAL_ARM_View(LocalGetNumObjectId(bsmodswopt),C_result);
	#endif

	#ifdef _DEBUG
		res = ARMLOCAL_ARM_View(LocalGetNumObjectId(swapId),C_result);
	#endif

	#ifdef _DEBUG
		res = ARMLOCAL_ARM_View(LocalGetNumObjectId(XStyle),C_result);
	#endif

	#ifdef _DEBUG
		res = ARMLOCAL_ARM_View(optionportId,C_result);
	#endif

	#ifdef _DEBUG
		res = ARMLOCAL_ARM_View(calibratorId,C_result);
	#endif

	#ifdef _DEBUG
		res = ARMLOCAL_ARM_View(markovtreeId,C_result);
	#endif

	long calibratorbumpId;
	long sfrmcalibratebumpId;
	long markovtreebumpId;
	if (computevega==true)
	{
		//******************CALIBRATOR SFRM BUMPED**********************
		
		if (ARMLOCAL_CALIBRATORSFRM(optionportId,
							LocalGetNumObjectId(bsmod),
							LocalGetNumObjectId(bsmodswoptbump),
							meanreversion,
							vCalibParams,
							preinitflagId,
							initsigmaVector,
							initbetaVector,
							ARM_NULL_OBJECT,
							ARM_NULL_OBJECT,
							ARM_NULL_OBJECT,
							K_DIAG,
							ARM_NULL_OBJECT,
							vCalibParamsPF,
							"YES",
							 C_result) == ARM_OK) 
		{
			calibratorbumpId = C_result.getLong();
		}
		
		//*******************************************************



		//******************SFRM CALIBRATE BUMPED*****************
		

		if (ARMLOCAL_SFRMCALIBRATE (calibratorbumpId,
								 ARM_NULL_OBJECT,
								 ARM_NULL_OBJECT,
								 "YES",
								 "YES",
								 kerneltogp,
								  C_result) == ARM_OK) 
		{
			sfrmcalibratebumpId = C_result.getLong();
		}
		//********************************************************


	  //*****************MARKOVTREE BUMPED*****************************
		

		if (ARMLOCAL_FRMMARKOVTREE (asofDate,
								 endDate,
								 LocalGetNumObjectId(zc),
								 markovtreepathnumber,
								 vMarkovTreeParams,
								sfrmcalibratebumpId,
								  C_result) == ARM_OK) 
		{
			markovtreebumpId = C_result.getLong();
		}
		//********************************************************



		//******************CALCUL VEGA********************
		if (ARMLOCAL_ARM_Price(optionportId,
							   markovtreebumpId,
							   C_result) == ARM_OK)
		{
			pxOptionBump= C_result.getDouble();
			pxOptionBump = pxOptionBump * 100;
		}

		vega = fabs(pxOptionBump-pxOption);
		//********************************************************

		#ifdef _DEBUG
			res = ARMLOCAL_ARM_View(LocalGetNumObjectId(zc),C_result);
		#endif

		#ifdef _DEBUG
			res = ARMLOCAL_ARM_View(LocalGetNumObjectId(bsmod),C_result);
		#endif

		#ifdef _DEBUG
			res = ARMLOCAL_ARM_View(LocalGetNumObjectId(bsmodswopt),C_result);
		#endif

		#ifdef _DEBUG
			res = ARMLOCAL_ARM_View(optionportId,C_result);
		#endif

		#ifdef _DEBUG
			res = ARMLOCAL_ARM_View(calibratorId,C_result);
		#endif

		#ifdef _DEBUG
			res = ARMLOCAL_ARM_View(markovtreeId,C_result);
		#endif
	}



	vRes.push_back(pxCorridorLeg);
	vRes.push_back(discMargin);
	vRes.push_back(pxOption);
	if (computevega==true)
	{
		vRes.push_back(vega);
	}
	
	retCode = ARMLOCAL_FreeObject(fixId,C_result);
	retCode = ARMLOCAL_FreeObject(rapayspreadonindexId,C_result);
	retCode = ARMLOCAL_FreeObject(corridorLegId,C_result);
	retCode = ARMLOCAL_FreeObject(corridorforCRALegId,C_result);
	retCode = ARMLOCAL_FreeObject(swapId,C_result);
	retCode = ARMLOCAL_FreeObject(optionportId,C_result);
	retCode = ARMLOCAL_FreeObject(calibratorId,C_result);

	if (computevega==true)
	{
		retCode = ARMLOCAL_FreeObject(calibratorbumpId,C_result);
		retCode = ARMLOCAL_FreeObject(sfrmcalibratebumpId,C_result);
		retCode = ARMLOCAL_FreeObject(markovtreebumpId,C_result);
	}
	retCode = ARMLOCAL_FreeObject(sfrmcalibrateId,C_result);
	retCode = ARMLOCAL_FreeObject(markovtreeId,C_result);
	return vRes;
}
catch(...)
{
	return resultat;
}
}


VECTOR<double> computeCRAFullFloat(double AsOf, 
							double fee,
							double startdate,
							double enddate,
							CCString l_ccy,
							long levelup,
							long upSpecId,
							long leveldown,
							long downSpecId,
							long indexRef,
							long refIndexId,
							long daycountId,
							long payFreqPayIndexId,
							long resetFreqRefIndexId,
							long paidRstTimingId,
							long refRstTimingId,
							long stubId,
							long porrId,
							double startcalldate,
							CCString l_XStyle,
							double KStyleId,
							double decompricingflag,
							double discMarginFactor,
							long preinitflagId,
							double MeanReversion,
							VECTOR<double> vCalibParams,
							VECTOR<double> vCalibParamsPF,
							double kerneltogp,
							VECTOR<double> vMarkovTreeParams,
							double markovtreepathnumber,
					        const CCString& bsmod,
							const CCString& bsmodswopt,
							const CCString& bsmodswoptbump,
							const CCString& zc,
							double liborlegId,
							double fund)

{	
	// TODO: Add your implementation code here
	
	VECTOR<double> resultat;
	resultat.clear();
try
{	
	ARM_result C_result;

	VECTOR<double> vRes;

 											
	//****************PAY INDEX*******************
	long payIndexFixId;
	long payIndexId;

	if (ARMLOCAL_FixedIndex(daycountId,l_ccy,C_result) == ARM_OK) 
	{
		payIndexFixId = C_result.getLong ();
	}

	if (ARMLOCAL_LIBOR(indexRef,
						false,
						l_ccy,
						payFreqPayIndexId,
						payFreqPayIndexId,
						daycountId,
						K_ADJUSTED,
						C_result) == ARM_OK) 
	{
		payIndexId = C_result.getLong ();
	}
  	//*********************************************

	//********************FLOATING LEG***********************
	long floatinglegId;
	if (ARMLOCAL_LIBORLEG(startdate,
						  enddate,
						  indexRef,
						  K_RCV,
						  0L,
						  0,
						  payFreqPayIndexId,
						  payFreqPayIndexId,
						  K_ADVANCE,
						  K_ARREARS,
						  false,
						  l_ccy,
						  K_ADJUSTED,
						  -2,
						  "NULL",
						  "NULL",
						  (long) decompricingflag,
						  K_NX_NONE,
						  K_SHORTSTART,
						  -1.0,
						  K_YES,
						  -1,
						  C_result) == ARM_OK) 
	{
		floatinglegId = C_result.getLong();
	}
	else
	{
		/*ERROR_MSG(C_result.getMsg(),pRet);*/
// FIXMEFRED: mig.vc8 (27/03/2007 17:08:00): problem S_OK != Vector<double>
		throw;
//		return S_OK;
	}
	
	double pxfloatingLeg;
	if (ARMLOCAL_ARM_Price(floatinglegId,
						   LocalGetNumObjectId(bsmod),
						   C_result) == ARM_OK)
	{
		pxfloatingLeg = C_result.getDouble();
	}
	else
	{
		/*ERROR_MSG(C_result.getMsg(),pRet);*/
// FIXMEFRED: mig.vc8 (27/03/2007 17:08:00): problem S_OK != Vector<double>
		throw;
//		return S_OK;
	}                 
	//*******************************************************

	//********************RA PAY INDEX FLOAT*****************
	long LDPricingMethodId;

	LDPricingMethodId = ARM_ConvPricCorridorLD("DIG");

	long rapayfloatindexlegId;

  	if (ARMLOCAL_CORRIDORLEG(startdate,
							 enddate,
							 porrId,
							 payIndexId,
							 payFreqPayIndexId,
							 0,
							 0,
							 refIndexId,
							 resetFreqRefIndexId,
							 paidRstTimingId,
							 refRstTimingId,
							 stubId,
							 leveldown,
							 downSpecId,
							 levelup,
							 upSpecId,
							 false,
							 l_ccy,
							 -1,
							 K_LINEAR,
							 LDPricingMethodId,
							 (long) decompricingflag,
							 "NULL", // reset calendar
							 "NULL", // pay calendar
                              0,
							  C_result) == ARM_OK) 
	{
		rapayfloatindexlegId = C_result.getLong();
	}
	
	double pxrapayfloatLeg;
	if (ARMLOCAL_ARM_Price(rapayfloatindexlegId,
						   LocalGetNumObjectId(bsmod),
						   C_result) == ARM_OK)
	{
		pxrapayfloatLeg = C_result.getDouble();
	}
	else
	{
		/*ERROR_MSG(C_result.getMsg(),pRet);*/
// FIXMEFRED: mig.vc8 (27/03/2007 17:08:00): problem S_OK != Vector<double>
		throw;
//		return S_OK;
	}                 
	//*******************************************************

	long retCode;

	double Xincr;
	int iteration =1;
	double SprSeek = 10.01;
	double SprSeekFinal= 10.01;
	double SprSeekBonifie = 10.01;  // PUREMENT ARBITRAIRE
	double CCNR;
	double CCNR2;
	double CCNRinter;

	CCNR=1000;

	double pxCorridorLeg;
	double discMargin;
	double pxOption;
	double vega;
	VECTOR<double> pxCRA;

	double pxCorridorLeg2;
	double discMargin2;
	double pxOption2;
	VECTOR<double> pxCRA2;


	double maxNegativeCCNR;
	double minPositiveCCNR;
	double maxNegativeCCNRRate;
	double minPositiveCCNRRate;
	double optimizationMethod;


	CCNR=1000;
	maxNegativeCCNR=-1000;
	minPositiveCCNR=1000;
	maxNegativeCCNRRate=10.01;
	minPositiveCCNRRate=0.01;
	optimizationMethod=0;


	while ((fabs(CCNR) > 2) && (iteration < 10 ))
	{
			optimizationMethod=0;
			SprSeek = SprSeekBonifie;
			pxCRA = computeCRAFloat(SprSeekBonifie,AsOf,startdate,enddate, l_ccy,levelup,upSpecId,
					leveldown,downSpecId,refIndexId,daycountId,porrId,payIndexId,payIndexFixId,
					payFreqPayIndexId,resetFreqRefIndexId,paidRstTimingId,refRstTimingId,
					stubId, (long) decompricingflag, discMarginFactor,startcalldate, l_XStyle,
					KStyleId,liborlegId,MeanReversion,preinitflagId,vCalibParams,
					vCalibParamsPF,vMarkovTreeParams,(long) markovtreepathnumber,
					kerneltogp,pxfloatingLeg,pxrapayfloatLeg,true, bsmod,bsmodswopt,bsmodswoptbump,zc);

			if (pxCRA.size() == 0)
			{
				//ERROR_MSG("Pb in computeCRAFix",pRet);
// FIXMEFRED: mig.vc8 (27/03/2007 17:08:00): problem S_OK != Vector<double>
		throw;
//		return S_OK;
			}

			pxCorridorLeg = pxCRA[0];
			discMargin = pxCRA[1];
			pxOption = pxCRA[2];
			vega=pxCRA[3];

			if (vega==0 && l_ccy=="USD" && discMarginFactor==4)
			{
				discMargin=2*discMargin;
			}	

			CCNR=fund+pxCorridorLeg+pxOption-discMargin-vega-fee;


			if (CCNR>0)
			{
				if (CCNR<minPositiveCCNR)
				{
					minPositiveCCNR=CCNR;
					minPositiveCCNRRate=SprSeekBonifie;
				}
			}
			else
			{
				if (CCNR>maxNegativeCCNR)
				{
					maxNegativeCCNR=CCNR;
					maxNegativeCCNRRate=SprSeekBonifie;
				}
			}

			Xincr = SprSeek / 1000000;
			SprSeekBonifie = SprSeek + Xincr;

			pxCRA2 = computeCRAFloat(SprSeekBonifie,AsOf,startdate,enddate, l_ccy,levelup,upSpecId,
					leveldown,downSpecId,refIndexId,daycountId,porrId,payIndexId,payIndexFixId,
					payFreqPayIndexId,resetFreqRefIndexId,paidRstTimingId,refRstTimingId,
					stubId, (long) decompricingflag, discMarginFactor,startcalldate, l_XStyle,
					KStyleId,liborlegId,MeanReversion,preinitflagId,vCalibParams,
					vCalibParamsPF,vMarkovTreeParams,(long) markovtreepathnumber,
					kerneltogp,pxfloatingLeg,pxrapayfloatLeg,false, bsmod,bsmodswopt,bsmodswoptbump,zc);

			if (pxCRA2.size() == 0)
			{
				//ERROR_MSG("Pb in computeCRAFix",pRet);
// FIXMEFRED: mig.vc8 (27/03/2007 17:08:00): problem S_OK != Vector<double>
		throw;
//		return S_OK;
			}

			pxCorridorLeg2 = pxCRA2[0];
			discMargin2 = pxCRA2[1];
			pxOption2 = pxCRA2[2];
			//vega=pxCRA[3];

			if (vega==0 && l_ccy=="USD" && discMarginFactor==4)
			{
				discMargin2=2*discMargin2;
			}

			CCNR2=fund+pxCorridorLeg2+pxOption2-discMargin2-vega-fee;

				if (CCNR2>0)
			{
				if (CCNR2<minPositiveCCNR)
				{
					minPositiveCCNR=CCNR2;
					minPositiveCCNRRate=SprSeekBonifie;
				}
			}
			else
			{
				if (CCNR2>maxNegativeCCNR)
				{
					maxNegativeCCNR=CCNR2;
					maxNegativeCCNRRate=SprSeekBonifie;
				}
			}
			

			CCNRinter=(CCNR2-CCNR) / Xincr;

			SprSeekBonifie = SprSeek -  (CCNR / CCNRinter) ; 

			iteration++;
	}

	//In this case, no convergence - dichotomy algorithm
	if (fabs(CCNR) > 2)
	{
		optimizationMethod=1;
		iteration=0;

		while ((fabs(CCNR) > 2) && (iteration < 10 ))
		{
			SprSeekBonifie=(minPositiveCCNRRate+maxNegativeCCNRRate)/2;
			pxCRA = computeCRAFloat(SprSeekBonifie,AsOf,startdate,enddate, l_ccy,levelup,upSpecId,
					leveldown,downSpecId,refIndexId,daycountId,porrId,payIndexId,payIndexFixId,
					payFreqPayIndexId,resetFreqRefIndexId,paidRstTimingId,refRstTimingId,
					stubId, (long) decompricingflag, discMarginFactor,startcalldate, l_XStyle,
					KStyleId,liborlegId,MeanReversion,preinitflagId,vCalibParams,
					vCalibParamsPF,vMarkovTreeParams,(long) markovtreepathnumber,
					kerneltogp,pxfloatingLeg,pxrapayfloatLeg,true, bsmod,bsmodswopt,bsmodswoptbump,zc);

			if (pxCRA.size() == 0)
			{
				//ERROR_MSG("Pb in computeCRAFloat",pRet);
// FIXMEFRED: mig.vc8 (27/03/2007 17:08:00): problem S_OK != Vector<double>
		throw;
//		return S_OK;
			}

			pxCorridorLeg = pxCRA[0];
			discMargin = pxCRA[1];
			pxOption = pxCRA[2];
			vega=pxCRA[3];

			if (vega==0 && l_ccy=="USD" && discMarginFactor==4)
			{
				discMargin=2*discMargin;
			}	

			CCNR=fund+pxCorridorLeg+pxOption-discMargin-vega-fee;
			
			if (CCNR>0)
			{
				minPositiveCCNRRate=SprSeekBonifie;
			}
			else
			{
				maxNegativeCCNRRate=SprSeekBonifie;
			}

			iteration++;
		}

	}


	vRes.push_back(SprSeek);
	vRes.push_back(fund);
	vRes.push_back(pxCorridorLeg);
	vRes.push_back(pxOption);
	vRes.push_back(discMargin);
	vRes.push_back(vega);
	vRes.push_back(optimizationMethod);

	retCode = ARMLOCAL_FreeObject(floatinglegId,C_result);
	retCode = ARMLOCAL_FreeObject(payIndexId,C_result);
	retCode = ARMLOCAL_FreeObject(payIndexFixId,C_result);
	retCode = ARMLOCAL_FreeObject(rapayfloatindexlegId,C_result);

	return vRes;


	}
catch(...)
{
	return resultat;
}
}



VECTOR<double> computeCRAFullFix(double AsOf, 
							double fee,
							double startdate,
							double enddate,
							CCString l_ccy,
							long levelup,
							long upSpecId,
							long leveldown,
							long downSpecId,
							long indexRef,
							long refIndexId,
							long daycountId,
							long payFreqPayIndexId,
							long resetFreqRefIndexId,
							long paidRstTimingId,
							long refRstTimingId,
							long stubId,
							long porrId,
							double startcalldate,
							CCString l_XStyle,
							double KStyleId,
							double decompricingflag,
							double discMarginFactor,
							long preinitflagId,
							double MeanReversion,
							VECTOR<double> vCalibParams,
							VECTOR<double> vCalibParamsPF,
							double kerneltogp,
							VECTOR<double> vMarkovTreeParams,
							double markovtreepathnumber,
					        const CCString& bsmod,
							const CCString& bsmodswopt,
							const CCString& bsmodswoptbump,
							const CCString& zc,
							double liborlegId,
							double fund)

{	
	// TODO: Add your implementation code here
	VECTOR<double> resultat;
	resultat.clear();
try
{		
	ARM_result C_result;

	double Res = 0.;
	long retCode;

	VECTOR<double> vRes;

	

	//****************PAY INDEX*******************
	long payIndexId;
	if (ARMLOCAL_FixedIndex(daycountId,l_ccy,C_result) == ARM_OK) 
	{
		payIndexId = C_result.getLong ();
	}
  	//*********************************************

  	double Xincr;
	int iteration =1;
	double SprSeek = 10.01;
	double SprSeekFinal= 10.01;
	double SprSeekBonifie = 10.01;  // PUREMENT ARBITRAIRE
	double CCNR;
	double CCNR2;
	double CCNRinter;

	double pxCorridorLeg;
	double discMargin;
	double pxOption;
	double vega;

	double pxCorridorLeg2;
	double discMargin2;
	double pxOption2;

	VECTOR<double> pxCRA;
	VECTOR<double> pxCRA2;

	double maxNegativeCCNR;
	double minPositiveCCNR;
	double maxNegativeCCNRRate;
	double minPositiveCCNRRate;
	double optimizationMethod;

	CCNR=1000;
	maxNegativeCCNR=-1000;
	minPositiveCCNR=1000;
	maxNegativeCCNRRate=10.01;
	minPositiveCCNRRate=0.01;
	optimizationMethod=0;


	while ((fabs(CCNR) > 2) && (iteration < 10 ))
	{
			optimizationMethod=0;
			SprSeek = SprSeekBonifie;
			pxCRA = computeCRAFix(SprSeekBonifie,AsOf,startdate,enddate, l_ccy,levelup,upSpecId,
					leveldown,downSpecId,refIndexId,daycountId,porrId,payIndexId,
					payFreqPayIndexId,resetFreqRefIndexId,paidRstTimingId,refRstTimingId,
					stubId, (long) decompricingflag, discMarginFactor,startcalldate, l_XStyle,
					KStyleId,liborlegId,MeanReversion,preinitflagId,vCalibParams,
					vCalibParamsPF,vMarkovTreeParams,(long) markovtreepathnumber,
					kerneltogp,true,bsmod,bsmodswopt,bsmodswoptbump,zc);

			if (pxCRA.size() == 0)
			{
				//ERROR_MSG("Pb in computeCRAFix",pRet);
// FIXMEFRED: mig.vc8 (27/03/2007 17:08:00): problem S_OK != Vector<double>
		throw;
//		return S_OK;
			}

			pxCorridorLeg = pxCRA[0];
			discMargin = pxCRA[1];
			pxOption = pxCRA[2];
			vega=pxCRA[3];

			if (vega==0 && l_ccy=="USD" && discMarginFactor==4)
			{
				discMargin=2*discMargin;
			}		
			
			CCNR=fund+pxCorridorLeg+pxOption-discMargin-vega-fee;

			if (CCNR>0)
			{
				if (CCNR<minPositiveCCNR)
				{
					minPositiveCCNR=CCNR;
					minPositiveCCNRRate=SprSeekBonifie;
				}
			}
			else
			{
				if (CCNR>maxNegativeCCNR)
				{
					maxNegativeCCNR=CCNR;
					maxNegativeCCNRRate=SprSeekBonifie;
				}
			}
			
		
			Xincr = SprSeek / 1000000;
			SprSeekBonifie = SprSeek + Xincr;

			pxCRA2 = computeCRAFix(SprSeekBonifie,AsOf,startdate,enddate, l_ccy,levelup,upSpecId,
					leveldown,downSpecId,refIndexId,daycountId,porrId,payIndexId,
					payFreqPayIndexId,resetFreqRefIndexId,paidRstTimingId,refRstTimingId,
					stubId, (long) decompricingflag, discMarginFactor,startcalldate, l_XStyle,
					KStyleId,liborlegId,MeanReversion,preinitflagId,vCalibParams,
					vCalibParamsPF,vMarkovTreeParams,(long) markovtreepathnumber,
					kerneltogp,false,bsmod,bsmodswopt,bsmodswoptbump,zc);

			if (pxCRA2.size() == 0)
			{
				//ERROR_MSG("Pb in computeCRAFix",pRet);
// FIXMEFRED: mig.vc8 (27/03/2007 17:08:00): problem S_OK != Vector<double>
		throw;
//		return S_OK;
			}

			pxCorridorLeg2 = pxCRA2[0];
			discMargin2 = pxCRA2[1];
			pxOption2 = pxCRA2[2];
			//vega=pxCRA2[3];


			if (vega==0 && l_ccy=="USD" && discMarginFactor==4)
			{
				discMargin2=2*discMargin2;
			}

			CCNR2=fund+pxCorridorLeg2+pxOption2-discMargin2-vega-fee;

			if (CCNR2>0)
			{
				if (CCNR2<minPositiveCCNR)
				{
					minPositiveCCNR=CCNR2;
					minPositiveCCNRRate=SprSeekBonifie;
				}
			}
			else
			{
				if (CCNR2>maxNegativeCCNR)
				{
					maxNegativeCCNR=CCNR2;
					maxNegativeCCNRRate=SprSeekBonifie;
				}
			}
			

			CCNRinter=(CCNR2-CCNR) / Xincr;

			SprSeekBonifie = SprSeek -  (CCNR / CCNRinter) ; 

			iteration++;
	}


	//In this case, no convergence - dichotomy algorithm
	if (fabs(CCNR) > 2)
	{
		optimizationMethod=1;
		iteration=0;

		while ((fabs(CCNR) > 2) && (iteration < 10 ))
		{
			SprSeekBonifie=(minPositiveCCNRRate+maxNegativeCCNRRate)/2;
			pxCRA = computeCRAFix(SprSeekBonifie,AsOf,startdate,enddate, l_ccy,levelup,upSpecId,
						leveldown,downSpecId,refIndexId,daycountId,porrId,payIndexId,
						payFreqPayIndexId,resetFreqRefIndexId,paidRstTimingId,refRstTimingId,
						stubId, (long) decompricingflag, discMarginFactor,startcalldate, l_XStyle,
						KStyleId,liborlegId,MeanReversion,preinitflagId,vCalibParams,
						vCalibParamsPF,vMarkovTreeParams,(long) markovtreepathnumber,
						kerneltogp,true,bsmod,bsmodswopt,bsmodswoptbump,zc);

			if (pxCRA.size() == 0)
			{
				/*ERROR_MSG("Pb in computeCRAFix",pRet);*/
// FIXMEFRED: mig.vc8 (27/03/2007 17:08:00): problem S_OK != Vector<double>
		throw;
//		return S_OK;
			}

			pxCorridorLeg = pxCRA[0];
			discMargin = pxCRA[1];
			pxOption = pxCRA[2];
			vega=pxCRA[3];

			if (vega==0 && l_ccy=="USD" && discMarginFactor==4)
			{
				discMargin=2*discMargin;
			}		
			
			CCNR=fund+pxCorridorLeg+pxOption-discMargin-vega-fee;

			if (CCNR>0)
			{
				minPositiveCCNRRate=SprSeekBonifie;
			}
			else
			{
				maxNegativeCCNRRate=SprSeekBonifie;
			}

			iteration++;
		}

	}

	retCode = ARMLOCAL_FreeObject(payIndexId,C_result);

	vRes.push_back(SprSeek);
	vRes.push_back(fund);
	vRes.push_back(pxCorridorLeg);
	vRes.push_back(pxOption);
	vRes.push_back(discMargin);
	vRes.push_back(vega);
	vRes.push_back(optimizationMethod);


	return vRes;
		
	}
catch(...)
{
	return resultat;
}
}

		  
STDMETHODIMP ActiveXModule::ARMcomputeCRA(BSTR pFixorFloat,
										  double pFee,
										  double pAsOf, 
											double pStartDate,
											double pEndDate,
											BSTR pCcy,
											double pLevelUp,
											BSTR pUpSpec,
											double pLevelDown,
											BSTR pDownSpec,
											BSTR pRefIndex,
											BSTR pDayCount,
											BSTR pPayFreqPayIndex,
											BSTR pResetFreqRefIndex,
											BSTR pPaidRstTiming,
											BSTR pRefRstTiming,
											BSTR pStubRule,
											BSTR pPOrR,
											double pStartCallDate,
											BSTR pXStyle,
											BSTR pFundingIndex,
											BSTR pResetFreqFunding,
											BSTR pPayFreqFunding,
											VARIANT pSpreadFunding,
											BSTR pPOrRFunding,
											double pDecompPricingFlag,
											double pdiscMarginFactor,
											BSTR pPreInitFlag,
											double pMeanReversion,
											VARIANT pCalibParams,
											VARIANT pCalibParamsPF,
											double pKernelToGP,
											VARIANT pMarkovTreeParams,
											double pMarkovTreePathNumber,
					                        BSTR pBsmodId,
								            BSTR pBsmodSwoptId,
											BSTR pBsmodSwoptBumpId,
											BSTR pzcId,
											VARIANT *pRet)

{	
	// TODO: Add your implementation code here
try
{		

	double Res = 0.;
	VECTOR<double> SprSeekFinal;

	long retCode;

	double dlevelup = pLevelUp;
    double dleveldown = pLevelDown;
	double startdate = pStartDate;
	double enddate = pEndDate;
	double decompricingflag = pDecompPricingFlag;

	double AsOf = pAsOf;
	double startcalldate = pStartCallDate;
	double fee=pFee;
	double discMarginFactor = pdiscMarginFactor;
  	double MeanReversion = pMeanReversion;
	double kerneltogp = pKernelToGP;
	double markovtreepathnumber=pMarkovTreePathNumber;

	_bstr_t bfixorfloat(pFixorFloat);
	CCString l_fixorfloat = bfixorfloat;
	ARM_result C_result;

	//***************CURRENCY***********************
   	_bstr_t bccy(pCcy);
	CCString l_ccy = bccy;
	//***********FIN CURRENCY********************

	//***********LEVEL UP AND LEVEL DOWN*********
	long levelup;  
	if (ARMLOCAL_CONSTREFVALUE(dlevelup,C_result) == ARM_OK)
	{
		levelup = C_result.getLong();
	}

	long leveldown;  
	if (ARMLOCAL_CONSTREFVALUE(dleveldown,C_result) == ARM_OK)
	{
		leveldown = C_result.getLong();
	}
	//***********FIN LEVEL UP AND LEVEL DOWN*****

	//************LEVEL SPEC**********************
	_bstr_t bupspec(pUpSpec);
	CCString l_upspec = bupspec;

  	long upSpecId;
	if((upSpecId = ARM_ConvStdBoost (l_upspec, C_result)) == ARM_DEFAULT_ERR)
	{
		ERROR_MSG("Invalid Level Spec",pRet,ARM_ERROR_FREQ);
		return S_OK;
	}

	_bstr_t bdownspec(pDownSpec);
	CCString l_downspec = bdownspec;

	long downSpecId;
	if((downSpecId = ARM_ConvStdBoost (l_downspec, C_result)) == ARM_DEFAULT_ERR)
	{
		ERROR_MSG("Invalid Level Spec",pRet,ARM_ERROR_FREQ);
		return S_OK;
	}
	//********************************************

	//************REF INDEX***********************
	_bstr_t brefIndex(pRefIndex);
	CCString l_refindex = brefIndex;

	long indexRef; 
	indexRef = ARM_ConvIrType(l_refindex); 
	//********************************************

	//****************REF INDEX OBJECT************
	long refIndexId;
	if(ARMLOCAL_LIBOR (indexRef,
					   false,
					   l_ccy,
					   -1,
					   -1,
                       -1,
					   K_ADJUSTED,
					   C_result) == ARM_OK)
	{
		refIndexId = C_result.getLong ();
	}
  	//*********************************************

	//***********DAY COUNT************************
	_bstr_t bdaycount(pDayCount);
	CCString l_daycount = bdaycount;

	long daycountId= ARM_ConvDayCount(l_daycount);
	//********************************************

	//************PAY FREQ PAY INDEX**************
	_bstr_t bpayfreqpayindex(pPayFreqPayIndex);
	CCString l_payfreqpayindex = bpayfreqpayindex;

	long payFreqPayIndexId;	
	if((payFreqPayIndexId = ARM_ConvFrequency(l_payfreqpayindex, C_result)) == ARM_DEFAULT_ERR)
	{
		ERROR_MSG("Invalid Frq",pRet,ARM_ERROR_FREQ);
		return S_OK;
	}
	//********************************************

	//*************RESET FREQ REF INDEX***********
	_bstr_t bresetfreqrefindex(pResetFreqRefIndex);
	CCString l_resetfreqrefindex = bresetfreqrefindex;

	long resetFreqRefIndexId;	
	if((resetFreqRefIndexId = ARM_ConvFrequency(l_resetfreqrefindex, C_result)) == ARM_DEFAULT_ERR)
	{
		ERROR_MSG("Invalid Frq",pRet,ARM_ERROR_FREQ);
		return S_OK;
	}
	//********************************************

	//*************PAID RESET TIMING**************
	_bstr_t bpaidrsttiming(pPaidRstTiming);
	CCString l_paidrsttiming = bpaidrsttiming;
	long paidRstTimingId = ARM_ConvPayResetRule (l_paidrsttiming);
	//**********************************************

	//*************REF RESET TIMING**************
	_bstr_t brefrsttiming(pRefRstTiming);
	CCString l_refrsttiming = brefrsttiming;
	long refRstTimingId = ARM_ConvPayResetRule (l_refrsttiming);
	//**********************************************

	//*************STUB RULE*********************
	_bstr_t bstubrule(pStubRule);
	CCString l_stubrule = bstubrule;
	long stubId = ARM_ConvStubRule(l_stubrule);
	//******************************************

	//*************P OR R************************
  	_bstr_t bporr(pPOrR);
	CCString l_porr = bporr;

	long porrId;

	if((porrId = ARM_ConvRecOrPay(l_porr, C_result)) == ARM_DEFAULT_ERR)
	{
		ERROR_MSG("Invalid PorR",pRet);
		return S_OK;
	}
	//**********************************************

	//**************XSTYLE**************************
	_bstr_t bXStyle(pXStyle);
	CCString l_XStyle= bXStyle;
	//**********************************************

	//************FUNDING INDEX*********************
	_bstr_t bfundingIndex(pFundingIndex);
	CCString l_fundingindex = bfundingIndex;

	long indexFundingId; 
	indexFundingId = ARM_ConvIrType(l_fundingindex); 
	//**********************************************

	//*************RESET FREQ FUNDING INDEX*********
	_bstr_t bresetfreqfundingindex(pResetFreqFunding);
	CCString l_resetfreqfundingindex = bresetfreqfundingindex;

	long resetFreqFundingIndexId;	
	if((resetFreqFundingIndexId = ARM_ConvFrequency(l_resetfreqfundingindex, C_result)) == ARM_DEFAULT_ERR)
	{
		ERROR_MSG("Invalid Frq",pRet,ARM_ERROR_FREQ);
		return S_OK;
	}
	//**********************************************

	//************PAY FREQ FUNDING INDEX************
	_bstr_t bpayfreqfundingindex(pPayFreqFunding);
	CCString l_payfreqfundingindex = bpayfreqfundingindex;

	long payFreqFundingIndexId;	
	if((payFreqFundingIndexId = ARM_ConvFrequency(l_payfreqfundingindex, C_result)) == ARM_DEFAULT_ERR)
	{
		ERROR_MSG("Invalid Frq",pRet,ARM_ERROR_FREQ);
		return S_OK;
	}
	//**********************************************

	//*************PRE INIT FLAG********************
	_bstr_t bpreinitflag(pPreInitFlag);
	CCString l_preinitflag = bpreinitflag;

	long preinitflagId;
	if ((preinitflagId = ARM_ConvYesOrNo (l_preinitflag, C_result)) == ARM_DEFAULT_ERR)
	{
	   ERROR_MSG("Invalid PreInitFlag",pRet);
	   return S_OK;
	}
	//************************************************

	//**************CALIB PARAMS**********************
  	VECTOR<double> vCalibParams;

	long sizecalibparams;

  	if(VARIANT2VECTORDOUBLE (pCalibParams,vCalibParams,sizecalibparams) != S_OK)
	{
		ERROR_MSG("Error conversion CalibParams",pRet,ARM_ERROR_FREQ);
		return S_OK;
	}
	//***************************************************

  	//**************CALIB PARAMS PF**********************
  	VECTOR<double> vCalibParamsPF;

	long sizecalibparamspf;

  	if(VARIANT2VECTORDOUBLE (pCalibParamsPF,vCalibParamsPF,sizecalibparamspf) != S_OK)
	{
		ERROR_MSG("Error conversion CalibParamsPF",pRet,ARM_ERROR_FREQ);
		return S_OK;
	}
	//***************************************************
	
	//**************MARKOV TREE PARAMS*********************
  	VECTOR<double> vMarkovTreeParams;

	long sizemarkovtreeparams;

  	if(VARIANT2VECTORDOUBLE (pMarkovTreeParams,vMarkovTreeParams,sizemarkovtreeparams) != S_OK)
	{
		ERROR_MSG("Error conversion MarkovTreeParams",pRet,ARM_ERROR_FREQ);
		return S_OK;
	}
	//***************************************************

	//********************MODELES**********************
	_bstr_t bbsmodId(pBsmodId);
	CCString l_bsmodId = bbsmodId;

	_bstr_t bbsmodswoptId(pBsmodSwoptId);
	CCString l_bsmodswoptId = bbsmodswoptId;

  	_bstr_t bbsmodswoptbumpId(pBsmodSwoptBumpId);
	CCString l_bsmodswoptbumpId = bbsmodswoptbumpId;

  	_bstr_t bzcId(pzcId);
	CCString l_zcId = bzcId;						
	//*************************************************

	//**************SPREAD FUNDING*********************
  	double SpreadFunding;
	long   SpreadFundingType;

	if (VARIANT2Double(pSpreadFunding,&SpreadFunding) == S_FALSE)
	{
		_bstr_t bSpreadFunding(pSpreadFunding);
		CCString l_SpreadFunding = bSpreadFunding;

		SpreadFunding = (double) LocalGetNumObjectId(l_SpreadFunding);
		SpreadFundingType = 1L;
	}
	else
	{
	   SpreadFundingType = 0L;
	}
	//**********************************************

	//*************P OR R FUNDING********************
  	_bstr_t bporrfunding(pPOrRFunding);
	CCString l_porrfunding = bporrfunding;

	long porrfundingId;

	if((porrfundingId = ARM_ConvRecOrPay(l_porrfunding, C_result)) == ARM_DEFAULT_ERR)
	{
		ERROR_MSG("Invalid PorR",pRet);
		return S_OK;
	}
	//**********************************************

	//**************** FUNDING *******************
	long liborlegId;
	if (ARMLOCAL_LIBORLEG(startdate,
						  enddate,
						  indexFundingId,
						  porrfundingId,
						  SpreadFundingType,
						  SpreadFunding,
						  resetFreqFundingIndexId,
						  payFreqFundingIndexId,
						  K_ADVANCE,
						  K_ARREARS,
						  false,
						  l_ccy,
						  K_ADJUSTED,
						  -2,
						  "NULL",
						  "NULL",
						  (long) decompricingflag,
						  K_NX_NONE,
						  K_SHORTSTART,
						  -1.0,
						  K_YES,
						  -1,
						  C_result) == ARM_OK) 
	{
		liborlegId = C_result.getLong();
	}
	else
	{
		ERROR_MSG(C_result.getMsg(),pRet);
		return S_OK;
	}

	
	double fund;
	if (ARMLOCAL_ARM_Price(liborlegId,
						   LocalGetNumObjectId(l_bsmodId),
						   C_result) == ARM_OK)
	{
		fund = C_result.getDouble();
		fund=fund*100;
	}
	else
	{
		ERROR_MSG(C_result.getMsg(),pRet);
		return S_OK;
	}                     
	//******************************************************
	#ifdef _DEBUG
		retCode = ARMLOCAL_ARM_View(liborlegId,C_result);
	#endif

	#ifdef _DEBUG
		retCode = ARMLOCAL_ARM_View(LocalGetNumObjectId(l_bsmodId),C_result);
	#endif

	#ifdef _DEBUG
		retCode = ARMLOCAL_ARM_View(LocalGetNumObjectId(l_bsmodswoptbumpId),C_result);
	#endif

	#ifdef _DEBUG
		retCode = ARMLOCAL_ARM_View(LocalGetNumObjectId(l_zcId),C_result);
	#endif

	#ifdef _DEBUG
		retCode = ARMLOCAL_ARM_View(LocalGetNumObjectId(l_XStyle),C_result);
	#endif
	//********************K STYLE***************************
	long KStyleId;  
	if (ARMLOCAL_CONSTREFVALUE(0, C_result) == ARM_OK)
	{
		KStyleId = C_result.getLong();
	}
  	//******************************************************

	if (l_fixorfloat=="FIX")
	{
			SprSeekFinal=computeCRAFullFix( AsOf, fee,startdate,	enddate,l_ccy,
				levelup,upSpecId,leveldown,downSpecId,indexRef,refIndexId,daycountId,
				payFreqPayIndexId,resetFreqRefIndexId,paidRstTimingId, refRstTimingId,
				 stubId, porrId, startcalldate,	 l_XStyle, KStyleId,
				  decompricingflag,discMarginFactor,preinitflagId,MeanReversion,
				  vCalibParams,vCalibParamsPF,kerneltogp,vMarkovTreeParams,markovtreepathnumber,
				  l_bsmodId,l_bsmodswoptId,l_bsmodswoptbumpId,l_zcId,liborlegId,fund);
	}
	if (l_fixorfloat=="FLOAT")
	{
			SprSeekFinal=computeCRAFullFloat( AsOf, fee,startdate,	enddate,l_ccy,
				levelup,upSpecId,leveldown,downSpecId,indexRef,refIndexId,daycountId,
				payFreqPayIndexId,resetFreqRefIndexId,paidRstTimingId, refRstTimingId,
				 stubId, porrId, startcalldate,	 l_XStyle,KStyleId, 
				  decompricingflag,discMarginFactor,preinitflagId,MeanReversion,
				  vCalibParams,vCalibParamsPF,kerneltogp,vMarkovTreeParams,markovtreepathnumber,
				  l_bsmodId,l_bsmodswoptId,l_bsmodswoptbumpId,l_zcId,liborlegId,fund);
	}

	retCode = ARMLOCAL_FreeObject(levelup,C_result);
	retCode = ARMLOCAL_FreeObject(leveldown,C_result);
	retCode = ARMLOCAL_FreeObject(liborlegId,C_result);
	retCode = ARMLOCAL_FreeObject(KStyleId,C_result);
	retCode = ARMLOCAL_FreeObject(refIndexId,C_result);
	
	Res = VECTORDOUBLE2VARIANT(SprSeekFinal, pRet);

	_variant_t wrap_Res;
	wrap_Res.Attach(*pRet);
	*pRet=wrap_Res.Detach();

	return S_OK;
	}
catch(...)
{
	return E_FAIL;
}
}

double computeHelvetixIIIDiscMargin(long underlyingId, //forexId
									double startDate,
									double endDate,
									double spot,
									CCString& domccy,
									CCString& forccy,
									long optionTypeId,
									long strikeRefValue,
									double epsilon,
									long payOffRefValue,
									long amortRefValue,
									double leverage,
									long freqId,
									long dayCountId,
									long resetTimingId,
									double resetGap,
									long intRuleId,
									long rollRuleId,
									VECTOR<double> vProba,
									VECTOR<double> vMT,
									const CCString& l_smiledModCHF,
									const CCString& l_mixtureModEURCHF)

{
	double Res=0;
	try
	{
		long retCode;
		ARM_result C_result;
		long i;
		//long j;
		//long k;
		

		i=0;
		//j=1;
		//k=1;
		//l=1;
		//m=1;
		//n=1;
		//o=1;


		double endDate_opt_i_for;
		double endDate_opt_i_dom;
		double startDate_opt_i_for;
		double startDate_opt_i_dom;
		double term_opt_i;
		double expiryDate_opt_i_dom;

	
		startDate_opt_i_for=startDate;
		startDate_opt_i_dom=startDate;

		double dPaymentDate;

		VECTOR<double> amortDate;
		VECTOR<double> amortValue;
		VECTOR<double> strikeDate;
		VECTOR<double> strikeValue;
		VECTOR<double> payOffDate;
		VECTOR<double> payOffValue;

		double strike_opt=0;
		double payoff_opt=0;
		double amort_opt=0;

		long nbrowsAmort;
		long nbrowsStrike;
		long nbrowsPayOff;
		long nbcolumns=2;

		long optionKPlusId;
		long optionKMinusId;

		double pxOptionKPlus;
		double pxOptionKMinus;
		double stripDigPrice;
		long fixedLegId;
		double proba;
		double discMargin=0;
		double fixedLegPrice;
		double interpol;

		if (ARMLOCAL_DisplayRefValue(amortRefValue,true,C_result) == ARM_OK)
		{
			nbrowsAmort = C_result.getLong();
			for (int x=0;x<nbrowsAmort;x++)
			{
				amortDate.push_back(C_result.getArray(x));
				amortValue.push_back(C_result.getArray(nbrowsAmort + x));
			}
		}

		if (ARMLOCAL_DisplayRefValue(strikeRefValue,true,C_result) == ARM_OK)
		{
			nbrowsStrike = C_result.getLong();
			for (int y=0;y<nbrowsStrike;y++)
			{
				strikeDate.push_back(C_result.getArray(y));
				strikeValue.push_back(C_result.getArray(nbrowsStrike + y));
			}
		}

  		if (ARMLOCAL_DisplayRefValue(payOffRefValue,true,C_result) == ARM_OK)
		{
			nbrowsPayOff = C_result.getLong();
			for (int z=0;z<nbrowsPayOff;z++)
			{
				payOffDate.push_back(C_result.getArray(z));
				payOffValue.push_back(C_result.getArray(nbrowsPayOff + z));
			}
		}
	
  
		while (endDate_opt_i_for < endDate)
		{
			
			strike_opt=0;
			payoff_opt=0;
			amort_opt=0;

			if (ARMLOCAL_ARM_ADDPERIOD(startDate, freqId, forccy, i+1, rollRuleId, 0, C_result)==ARM_OK)
			{
					endDate_opt_i_for=Local_ARMDATE2XLDATE(C_result.getString ());
					
			}

			if (ARMLOCAL_ARM_BetweenDates(startDate_opt_i_for, endDate_opt_i_for, dayCountId, 1, C_result)==ARM_OK)
			{
					term_opt_i=C_result.getDouble();
			}


			if (ARMLOCAL_ARM_ADDPERIOD(startDate, freqId, domccy, i+1, rollRuleId, 0, C_result)==ARM_OK)
			{
					endDate_opt_i_dom=Local_ARMDATE2XLDATE(C_result.getString ());
			}

			if (resetTimingId==K_ARREARS)
			{
				if (ARMLOCAL_ARM_ADDPERIOD(endDate_opt_i_dom, K_DAILY, domccy, resetGap, rollRuleId, 0, C_result)==ARM_OK)
				{
					expiryDate_opt_i_dom=Local_ARMDATE2XLDATE(C_result.getString ());
				}
				
			}
			else
			{
				if (ARMLOCAL_ARM_ADDPERIOD(startDate_opt_i_dom, K_DAILY, domccy, resetGap, rollRuleId, 0, C_result)==ARM_OK)
				{
					expiryDate_opt_i_dom=Local_ARMDATE2XLDATE(C_result.getString ());
				}
				
			}

			bool wentthere = false;
			
			for (int j=0;j<nbrowsStrike;j++)
			{	
				if (wentthere==false)
				{
					if ((strikeDate[j] == expiryDate_opt_i_dom) || (abs(expiryDate_opt_i_dom-strikeDate[j]) <= 2))
					{
						strike_opt=strikeValue[j];
						wentthere=true;
					}
				}
				
			}

			wentthere = false;

			for (int k=0;k<nbrowsPayOff;k++)
			{	
				if (wentthere==false)
				{
					if ((payOffDate[k] == expiryDate_opt_i_dom) || (abs(expiryDate_opt_i_dom-payOffDate[k]) <= 2))
					{
						payoff_opt=payOffValue[k];
						wentthere=true;
					}
				}
				
			}

			wentthere = false;

			for (int l=0;l<nbrowsAmort;l++)
			{	
				if (wentthere==false)
				{
					if ((amortDate[l] == endDate_opt_i_for) || (abs(endDate_opt_i_for-amortDate[l]) <= 2))
					{
						amort_opt=amortValue[l];
						wentthere=true;
					}
				}
				
			}

			//création de l'option K+
			if (ARMLOCAL_OPTION (underlyingId, expiryDate_opt_i_dom, strike_opt+epsilon,optionTypeId,K_EUROPEAN,0,-1, endDate_opt_i_for, C_result)==ARM_OK)
			{
				optionKPlusId=C_result.getLong();
			}

			if (ARMLOCAL_ARM_Price(optionKPlusId,
							   LocalGetNumObjectId(l_mixtureModEURCHF),
							   C_result) == ARM_OK)
			{
				pxOptionKPlus = C_result.getDouble();
			}

			retCode = ARMLOCAL_FreeObject(optionKPlusId,C_result);

			//Création de l'option K-
			if (ARMLOCAL_OPTION (underlyingId, expiryDate_opt_i_dom, strike_opt-epsilon,optionTypeId,K_EUROPEAN,0,-1, endDate_opt_i_for, C_result)==ARM_OK)
			{
				optionKMinusId=C_result.getLong();
			}

			if (ARMLOCAL_ARM_Price(optionKMinusId,
							   LocalGetNumObjectId(l_mixtureModEURCHF),
							   C_result) == ARM_OK)
			{
				pxOptionKMinus = C_result.getDouble();
			}

			retCode = ARMLOCAL_FreeObject(optionKMinusId,C_result);

			stripDigPrice = leverage * amort_opt * term_opt_i * payoff_opt * (pxOptionKPlus - pxOptionKMinus) / (2 * epsilon);

			if (ARMLOCAL_FIXEDLEG(startDate_opt_i_for,
							  endDate_opt_i_for,
							  K_RCV,
							  0,
							  payoff_opt,
							  dayCountId,
							  freqId,
							  K_COMP_PROP,
							  K_ARREARS,
							  intRuleId,
							  K_SHORTSTART,
							  false,
							  domccy,
							  "NULL",
							  K_NX_NONE,
							  -1.0,
							  K_YES,
							  -1,
							  C_result) == ARM_OK) 
			{
				fixedLegId=C_result.getLong();
			}



  			if (ARMLOCAL_ARM_Price(fixedLegId,
								   LocalGetNumObjectId(l_smiledModCHF),
								   C_result) == ARM_OK)
			{
				fixedLegPrice = C_result.getDouble();
			}

			#ifdef _DEBUG
				retCode = ARMLOCAL_ARM_View(fixedLegId,C_result);
			#endif
		
			
			retCode = ARMLOCAL_FreeObject(fixedLegId,C_result);
	

			fixedLegPrice = leverage * amort_opt * fixedLegPrice;

			proba = stripDigPrice / fixedLegPrice;

			if (ARMLOCAL_INTERPOL(vProba, vMT, proba, K_LINEAR, C_result) == ARM_OK)
			{
				interpol=C_result.getDouble();
			}
			
			discMargin = discMargin + interpol*fixedLegPrice;
			
			
			startDate_opt_i_for = endDate_opt_i_for;
			startDate_opt_i_dom = endDate_opt_i_dom;

			i++;
		}

		Res = discMargin/spot;

		return Res;
	}
	catch(...)
	{
		return Res;
	}
}


VECTOR<double> computeHelvetix(double tauxBonifie,
								double asOfDate,
								double startDate,
								double startDatePhaseStruct,
								double endDatePhaseStruct,
								double endDate,
								double helvetixType,
								double spot,
								long amortRefValue,
								CCString& l_ccy, 
								long stubRuleId,
								bool isSmoothedPhasePre,
								long irIndexPhasePreId,
								long dayCountPhasePreId,
								double resetGapPhasePre,
								bool isSmoothedPhasePost,
								long irIndexPhasePostId,
								long dayCountPhasePostId,
								double resetGapPhasePost,
								long dayCountPhaseStructId,
								double resetGapPhaseStruct,
								long irIndexPhaseStructId,				
								long iscappedId, 
								double tauxCap,
								const CCString& l_strikePhaseStruct,
								long forexId,
								long optionTypePhaseStructId,
								CCString& l_foreignccy, 
								CCString& l_domesticccy, 
								long payFreqOptionPhaseStructId,
								long dayCountOptionPhaseStructId,
								long fwdRuleId,
								long adjOptionPhaseStructId,
								double resetGapOptionPhaseStruct,
								double payGapOptionPhaseStruct,
  								long resetTimingOptionPhaseStructId,
								long payTimingOptionPhaseStructId,
								double levier,
								const CCString& l_smiledModEUR,
								const CCString& l_deltaSmiledModEUR,
  								const CCString& l_mixtureModEURCHF,
								const CCString& l_deltaEURMixtureMod,
								const CCString& l_bSEURMixtureMod,
								const CCString& l_volFxEURMixtureMod,
								const CCString& l_deltaCHFMixtureMod,
								const CCString& l_bSCHFMixtureMod)
{
	VECTOR<double> vRes;
	vRes.clear();
try
{
	long retCode;
	ARM_result C_result;

	double pxFinalPhasePre=0;
	double pxFinalPhaseStruct=0;
	double pxFinalPhasePost=0;
	double pxFinalCap=0;

	double pxBumpPhasePre=0;
	double pxBumpPhaseStruct=0;
	double pxBumpPhasePost=0;
	double pxBumpDeltaEURCap=0;
	double pxBumpDeltaCHFCap=0;
	double pxBumpBsEURCap=0;
	double pxBumpBsCHFCap=0;
	double pxBumpVolFxEURCap=0;

	double deltaFinalPhasePre=0;
	double deltaFinalPhasePost=0;
	double deltaFinalPhaseStruct=0;
	double deltaEURFinalCap=0;
	double deltaCHFFinalCap=0;
	double bsEURFinalCap=0;
	double bsCHFFinalCap=0;
	double volFxFinalCap=0;

	double levierCap;
	double strikeCapHIII;

	VECTOR<double> strikeDate;
	VECTOR<double> strikeValue;
	VECTOR<double> capStrikeValue;
	VECTOR<double> valueVide;

	long capStrikeRefValue;
	
	long nbrowsStrike;

	//Phase Pre
	if (isSmoothedPhasePre==true)
	{
		long swapLegPhasePreId;
		if (ARMLOCAL_SWAPLEG(irIndexPhasePreId,
							startDate,
							startDatePhaseStruct,
							K_RCV,
							0, 
							tauxBonifie, 
							false,
							l_ccy,
							dayCountPhasePreId,
							resetGapPhasePre,
							"NULL",
							"NULL",
							0,
							K_NX_NONE,
							stubRuleId,
							-1.0,
							K_NO,
							-1,
							C_result) == ARM_OK) 
		{
			swapLegPhasePreId = C_result.getLong();
		}
	

		long clonePhasePreId;
		if (ARMLOCAL_ClonedAndSetNotional(swapLegPhasePreId,
									  amortRefValue,
									  100.,
									  C_result) == ARM_OK)
		{
			clonePhasePreId = C_result.getLong();
		}
		

		if (ARMLOCAL_ARM_Price(clonePhasePreId,
							   LocalGetNumObjectId(l_smiledModEUR),
							   C_result) == ARM_OK)
		{
			pxFinalPhasePre = C_result.getDouble();
		}

		//Delta Phase Pre

		if (ARMLOCAL_ARM_Price(clonePhasePreId,
							   LocalGetNumObjectId(l_deltaSmiledModEUR),
							   C_result) == ARM_OK)
		{
			pxBumpPhasePre = C_result.getDouble();
		}
	

		deltaFinalPhasePre=(pxBumpPhasePre-pxFinalPhasePre);
		

		retCode = ARMLOCAL_FreeObject(swapLegPhasePreId,C_result);
		retCode = ARMLOCAL_FreeObject(clonePhasePreId,C_result);
	}

	//PhaseStruct
	long swapLegPhaseStructId;
	if (ARMLOCAL_SWAPLEG(irIndexPhaseStructId,
						startDatePhaseStruct,
						endDatePhaseStruct,
						K_RCV,
						0, 
						tauxBonifie, 
						false,
						l_ccy,
						dayCountPhaseStructId,
						resetGapPhaseStruct,
						"NULL",
						"NULL",
						0,
						K_NX_NONE,
						stubRuleId,
						-1.0,
						K_NO,
						-1,
						C_result) == ARM_OK) 
	{
		swapLegPhaseStructId = C_result.getLong();
	}


	long clonePhaseStructId;
	if (ARMLOCAL_ClonedAndSetNotional(swapLegPhaseStructId,
								  amortRefValue,
								  100.,
								  C_result) == ARM_OK)
	{
		clonePhaseStructId = C_result.getLong();
	}
	

	if (ARMLOCAL_ARM_Price(clonePhaseStructId,
						   LocalGetNumObjectId(l_smiledModEUR),
						   C_result) == ARM_OK)
	{
		pxFinalPhaseStruct = C_result.getDouble();
	}

	//Delta Phase Struct

	if (ARMLOCAL_ARM_Price(clonePhaseStructId,
						   LocalGetNumObjectId(l_deltaSmiledModEUR),
						   C_result) == ARM_OK)
	{
		pxBumpPhaseStruct = C_result.getDouble();
	}
	

	deltaFinalPhaseStruct=(pxBumpPhaseStruct-pxFinalPhaseStruct);


	retCode = ARMLOCAL_FreeObject(swapLegPhaseStructId,C_result);
	retCode = ARMLOCAL_FreeObject(clonePhaseStructId,C_result);

	//Phase Post
	if (isSmoothedPhasePost==true)
	{
		long swapLegPhasePostId;
		if (ARMLOCAL_SWAPLEG(irIndexPhasePostId,
							endDatePhaseStruct,
							endDate,
							K_RCV,
							0, 
							tauxBonifie, 
							false,
							l_ccy,
							dayCountPhasePostId,
							resetGapPhasePost,
							"NULL",
							"NULL",
							0,
							K_NX_NONE,
							stubRuleId,
							-1.0,
							K_NO,
							-1,
							C_result) == ARM_OK) 
		{
			swapLegPhasePostId = C_result.getLong();
		}
		

		long clonePhasePostId;
		if (ARMLOCAL_ClonedAndSetNotional(swapLegPhasePostId,
									  amortRefValue,
									  100.,
									  C_result) == ARM_OK)
		{
			clonePhasePostId = C_result.getLong();
		}
	

		if (ARMLOCAL_ARM_Price(clonePhasePostId,
							   LocalGetNumObjectId(l_smiledModEUR),
							   C_result) == ARM_OK)
		{
			pxFinalPhasePost = C_result.getDouble();
		}

		//Delta Phase Post

		if (ARMLOCAL_ARM_Price(clonePhasePostId,
							   LocalGetNumObjectId(l_deltaSmiledModEUR),
							   C_result) == ARM_OK)
		{
			pxBumpPhasePost = C_result.getDouble();

		}
		

		deltaFinalPhasePost=(pxBumpPhasePost-pxFinalPhasePost);
		

		retCode = ARMLOCAL_FreeObject(swapLegPhasePostId,C_result);
		retCode = ARMLOCAL_FreeObject(clonePhasePostId,C_result);
	}
	
	//Cap éventuel
	if (iscappedId==K_YES)
	{
		levierCap =  levier * ( 1 + ((tauxCap - tauxBonifie) / (100*levier)));
		long fxCapOptionId;

		if (helvetixType==3)
		{
			//Un seul strike dépendant du Taux Fixe Bonifié et indépendant des strikes du put
			strikeCapHIII = spot / (1 + ((tauxCap-tauxBonifie) / (100*levier)));

			if (ARMLOCAL_CONSTREFVALUE(strikeCapHIII,C_result)==ARM_OK)
			{
				capStrikeRefValue =  C_result.getLong();
			}
			
		}

		if (helvetixType==2)
		{
			//Un échéancier de strike fonction du TFB et des strike du put
			//on commence par décomposer le strike du put

			long strikeRefValue = LocalGetNumObjectId(l_strikePhaseStruct);

			if (ARMLOCAL_DisplayRefValue(strikeRefValue,true,C_result) == ARM_OK)
			{
				nbrowsStrike = C_result.getLong();
				for (int y=0;y<nbrowsStrike;y++)
				{
					strikeDate.push_back(C_result.getArray(y));
					strikeValue.push_back(C_result.getArray(nbrowsStrike + y));
				}
			}

			for (int j=0;j<nbrowsStrike;j++)
			{	
				capStrikeValue.push_back(strikeValue[j] / (1 + ((tauxCap-tauxBonifie) / (100*levier))));
			}

			if(ARMLOCAL_REFVALUE(strikeDate,capStrikeValue,valueVide, 1, 1, K_STEPUP_RIGHT,C_result) == ARM_OK)
			{
				capStrikeRefValue = C_result.getLong ();
			}
		}

		if (ARMLOCAL_FxOptionStrip(XLDateToJulian(asOfDate),
							 forexId,
							 capStrikeRefValue,
							 optionTypePhaseStructId,
							 startDatePhaseStruct,
							 endDatePhaseStruct,
							 amortRefValue,
							 l_foreignccy,
							 payFreqOptionPhaseStructId,
							 dayCountOptionPhaseStructId,
							 l_domesticccy,
							 fwdRuleId,
							 adjOptionPhaseStructId,
							 stubRuleId,
							 resetGapOptionPhaseStruct,
							 payFreqOptionPhaseStructId, 
							 payGapOptionPhaseStruct,
							 l_foreignccy,
							 resetTimingOptionPhaseStructId,
							 payTimingOptionPhaseStructId,
							 "RCV",							 
							 ARM_NULL_OBJECT,
							 false,
							 1,
							 0.1,
							 ARM_NULL_OBJECT,
							 ARM_NULL_OBJECT,
							 levierCap,
							 C_result) == ARM_OK)
		{
			fxCapOptionId=C_result.getLong();
		}

	
		if (ARMLOCAL_ARM_Price(fxCapOptionId,
							   LocalGetNumObjectId(l_mixtureModEURCHF),
							   C_result) == ARM_OK)
		{
			pxFinalCap = C_result.getDouble();
		}

		if (ARMLOCAL_ARM_Price(fxCapOptionId,
							   LocalGetNumObjectId(l_deltaEURMixtureMod),
							   C_result) == ARM_OK)
		{
			pxBumpDeltaEURCap = C_result.getDouble();
		}

		deltaEURFinalCap = (pxBumpDeltaEURCap-pxFinalCap);

		if (ARMLOCAL_ARM_Price(fxCapOptionId,
							   LocalGetNumObjectId(l_deltaCHFMixtureMod),
							   C_result) == ARM_OK)
		{
			pxBumpDeltaCHFCap = C_result.getDouble();
		}

		deltaCHFFinalCap = (pxBumpDeltaCHFCap-pxFinalCap);

		if (ARMLOCAL_ARM_Price(fxCapOptionId,
							   LocalGetNumObjectId(l_bSEURMixtureMod),
							   C_result) == ARM_OK)
		{
			pxBumpBsEURCap = C_result.getDouble();
		}

		bsEURFinalCap=(pxBumpBsEURCap-pxFinalCap);

		if (ARMLOCAL_ARM_Price(fxCapOptionId,
							   LocalGetNumObjectId(l_bSCHFMixtureMod),
							   C_result) == ARM_OK)
		{
			pxBumpBsCHFCap = C_result.getDouble();
		}

		bsCHFFinalCap=(pxBumpBsCHFCap-pxFinalCap);

		if (ARMLOCAL_ARM_Price(fxCapOptionId,
							   LocalGetNumObjectId(l_volFxEURMixtureMod),
							   C_result) == ARM_OK)
		{
			pxBumpVolFxEURCap = C_result.getDouble();
		}

		volFxFinalCap=(pxBumpVolFxEURCap-pxFinalCap);

	
		retCode = ARMLOCAL_FreeObject(capStrikeRefValue,C_result);
		retCode = ARMLOCAL_FreeObject(fxCapOptionId,C_result);
	}

		
	vRes.push_back(pxFinalPhasePre);
	vRes.push_back(pxFinalPhaseStruct);
	vRes.push_back(pxFinalCap);
	vRes.push_back(pxFinalPhasePost);
	vRes.push_back(deltaFinalPhasePre);
	vRes.push_back(deltaFinalPhaseStruct);
	vRes.push_back(deltaFinalPhasePost);
	vRes.push_back(deltaEURFinalCap);
	vRes.push_back(deltaCHFFinalCap);
	vRes.push_back(bsEURFinalCap);
	vRes.push_back(bsCHFFinalCap);
	vRes.push_back(volFxFinalCap);

	return vRes;
}
catch(...)
{
	return vRes;
}
}


STDMETHODIMP ActiveXModule::ARMcomputeHelvetix(double pNtl, 
										   double pAsOfDate, 
										   double pStartDate, 
										   double pEndDate,
										   double pHelvetixType,
										   double pFee,
										   BSTR pCcy,
										   BSTR pIndexFund,
										   VARIANT pSpreadFund,
										   BSTR pDayCountFund, 
										   BSTR pPayFreqFund, 
										   BSTR pResetFreqFund,
										   BSTR pResetTimingFund,
										   BSTR pAdjFund, 
										   BSTR pDayCountAmort,
										   BSTR pIntRuleAmort, 
										   double pTxAmort,
										   BSTR pFreqAmort,
										   double pAmountAmort, 
										   BSTR pTypeAmort, 
										   BSTR pIndexPhasePre, 										   
										   double pSpreadPhasePre,
										   BSTR pDayCountPhasePre, 
										   BSTR pPayFreqPhasePre,
										   BSTR pResetFreqPhasePre,
										   BSTR pResetTimingPhasePre,
										   BSTR pAdjPhasePre, 
										   double pStartDatePhaseStruct,
										   double pEndDatePhaseStruct,
										   BSTR pResetTimingPhaseStruct,
										   BSTR pAdjPhaseStruct,										   
										   BSTR pDayCountPhaseStruct, 
										   BSTR pPayFreqPhaseStruct, 
										   BSTR pResetFreqPhaseStruct,
										   double pSpotPhaseStruct,
										   BSTR pOptionTypePhaseStruct,
										   BSTR pForeignCcyPhaseStruct,
										   BSTR pDomesticCcyPhaseStruct,
										   double pLevierPhaseStruct,
										   BSTR pStrikePhaseStruct, //objet représentant les barrières
										   BSTR pPayOffDigitalPhaseStruct, //objet représentant les pay off de l'éventuelle digitale
										   BSTR pResetTimingOptionPhaseStruct,
										   BSTR pAdjOptionPhaseStruct,										   
										   BSTR pDayCountOptionPhaseStruct, 
										   BSTR pPayFreqOptionPhaseStruct, 
										   BSTR pResetFreqOptionPhaseStruct,
										   BSTR pPayTimingOptionPhaseStruct,
										   double pPayGapOptionPhaseStruct,
										   double pResetGapOptionPhaseStruct,								   
										   BSTR pIsCapped,
   										   double pTxCap, 
										   BSTR pIndexPhasePost,
										   double pSpreadPhasePost,
										   BSTR pDayCountPhasePost, 
										   BSTR pPayFreqPhasePost, 
										   BSTR pResetFreqPhasePost,
										   BSTR pResetTimingPhasePost,
										   BSTR pAdjPhasePost,
										   BSTR pRoll, 
										   BSTR pFwdRule,
										   BSTR pIntRule,
										   BSTR pStubRule,
										   VARIANT *pProba,
										   VARIANT *pMT,
										   BSTR pSmiledModEUR,
										   BSTR pSmiledModCHF,
										   BSTR pDeltaSmiledModEUR,
										   BSTR pMixtureModEURCHF,										   
										   BSTR pDeltaEURMixtureMod,
										   BSTR pBSEURMixtureMod,
										   BSTR pVolFxEURMixtureMod,
										   BSTR pDeltaCHFMixtureMod,
										   BSTR pBSCHFMixtureMod,
										   VARIANT *pRet)
{	
try
{									   
	//**************Init Global*******************************
	ARM_result C_result;

	double asOfDate = pAsOfDate; 
	double startDate = pStartDate;  
	double endDate = pEndDate; 
	double startDatePhaseStruct = pStartDatePhaseStruct;
	double endDatePhaseStruct = pEndDatePhaseStruct;
	double helvetixType = pHelvetixType;

	double putDeltaEUR=0;
	double putDeltaCHF=0;
	double capDeltaEUR=0;
	double capDeltaCHF=0;
	double putBSEUR=0;
	double putBSCHF=0;
	double capBSEUR=0;
	double capBSCHF=0;
	double putVolFxEUR=0;
	double capVolFxEUR=0;
	double fundingDeltaEUR=0;
	double phasePreDeltaEUR=0;
	double phasePostDeltaEUR=0;
	double phaseStructDeltaEUR=0;

	double pxPut=0;
	double pxPutBumpDeltaEUR=0;
	double pxPutBumpDeltaCHF=0;
	double pxPutBumpBSEUR=0;
	double pxPutBumpBSCHF=0;
	double pxPutBumpFxVolEUR=0;

	double pxDigital=0;
	double pxDigitalBumpDeltaEUR=0;
	double pxDigitalBumpDeltaCHF=0;
	double pxDigitalBumpBSEUR=0;
	double pxDigitalBumpBSCHF=0;
	double pxDigitalBumpFxVolEUR=0;

	double digitalDeltaEUR=0;
	double digitalDeltaCHF=0;
	double digitalBSEUR=0;
	double digitalBSCHF=0;
	double digitalVolFxEUR=0;

	double pxPhasePost=0; 
	double pxBumpPhasePost=0; 
	double pxPhasePre=0; 
	double pxBumpPhasePre=0;

	double pxcertainleg=0;
	double proba=0;
	double discMargin=0;

	double Ntl = pNtl;
	double fee = pFee; 

	long retCode;
	double Res = 0.;

	//Pay Currency = EUR
	_bstr_t ccy(pCcy);
	CCString l_ccy = ccy;
	long ccyId; 
	if (ARMLOCAL_ISOCCY(l_ccy,C_result) == ARM_OK)
	{
		ccyId = C_result.getLong();
	}
	else
	{
		ERROR_MSG(C_result.getMsg(),pRet);
		return S_OK;
	}

	_bstr_t roll(pRoll);
	CCString l_roll = roll;

	long rollId;
	if((rollId = ARM_ConvFwdRule(l_roll, C_result)) == ARM_DEFAULT_ERR)
	{
		ERROR_MSG("Invalid Int Roll",pRet,ARM_ERROR_FREQ);
		return S_OK;
	}

	_bstr_t stubRule(pStubRule);
	CCString l_stub= stubRule;
	long stubRuleId = ARM_ConvStubRule(l_stub);
	
	_bstr_t fwdRule(pFwdRule);
	CCString l_fwd= fwdRule;
	long fwdRuleId;
	if((fwdRuleId = ARM_ConvFwdRule(l_fwd, C_result)) == ARM_DEFAULT_ERR)
	{
		ERROR_MSG("Invalid Fwd Rule",pRet,ARM_ERROR_FREQ);
		return S_OK;
	}

	_bstr_t intRule(pIntRule);
	CCString l_int= intRule;
	long intRuleId= ARM_ConvIntRule(l_int);
	//***************End init Global*************************

  	//***************Init Market Datas***********************
	_bstr_t smiledModEUR (pSmiledModEUR);
	CCString l_smiledModEUR = smiledModEUR;

	_bstr_t smiledModCHF (pSmiledModCHF);
	CCString l_smiledModCHF = smiledModCHF;

	_bstr_t deltaSmiledModEUR (pDeltaSmiledModEUR);
	CCString l_deltaSmiledModEUR = deltaSmiledModEUR;

  	_bstr_t mixtureModEURCHF (pMixtureModEURCHF);
	CCString l_mixtureModEURCHF = mixtureModEURCHF;

	_bstr_t deltaEURMixtureMod (pDeltaEURMixtureMod);
	CCString l_deltaEURMixtureMod = deltaEURMixtureMod;
	
	_bstr_t bSEURMixtureMod (pBSEURMixtureMod);
	CCString l_bSEURMixtureMod = bSEURMixtureMod;
	
	_bstr_t volFxEURMixtureMod (pVolFxEURMixtureMod);
	CCString l_volFxEURMixtureMod = volFxEURMixtureMod;
	
	_bstr_t deltaCHFMixtureMod (pDeltaCHFMixtureMod);
	CCString l_deltaCHFMixtureMod = deltaCHFMixtureMod;
	
	_bstr_t bSCHFMixtureMod (pBSCHFMixtureMod);
	CCString l_bSCHFMixtureMod = bSCHFMixtureMod;								   
	//***************End Init Market Datas***********************

	//***************Init Amort**************************
	_bstr_t typeAmort(pTypeAmort);
	CCString l_typeAmort = typeAmort;

  	long amortRefValue;
	long amortRefValueToDestroy = ARM_NULL_OBJECT;

	if (LocalGetNumObjectId(l_typeAmort) == ARM_KO)
	{
		//double endDateAmort = pEndDateAmort;
		double txAmort = pTxAmort;
		double amountAmort = pAmountAmort;

		_bstr_t intRuleAmort(pIntRuleAmort);
		CCString l_intRuleAmort =  intRuleAmort;

		_bstr_t freqAmort(pFreqAmort);
		CCString l_freqAmort = freqAmort;

		_bstr_t dayCountAmort(pDayCountAmort);
		CCString l_dayCountAmort = dayCountAmort;

		long freqAmortId; 
		if((freqAmortId = ARM_ConvFrequency(l_freqAmort, C_result)) == ARM_DEFAULT_ERR)
		{
			ERROR_MSG("Invalid Frq",pRet,ARM_ERROR_FREQ);
			return S_OK;
		}

		long intRuleAmortId= ARM_ConvIntRule(l_intRuleAmort);
		long dayCountAmortId = ARM_ConvDayCount(l_dayCountAmort);	

		// AmortMethodId
		long amortMethodId;

		if((amortMethodId = ARM_ConvAmortMethod (l_typeAmort, C_result)) == ARM_DEFAULT_ERR)
		{
			ERROR_MSG("Invalid TypeAmort",pRet,ARM_ERROR_FREQ);
			return S_OK;
		}

		long FixedLegAmortId;
		if (ARMLOCAL_FIXEDLEG(startDate,
							  endDate,
							  K_RCV,
							  0,
							  1.0,
							  dayCountAmortId,
							  freqAmortId,
							  K_COMP_PROP,
							  K_ARREARS,
							  intRuleAmortId,
							  stubRuleId,
							  false,
							  l_ccy,
							  "NULL",
							  K_NX_NONE,
							  -1.0,
							  K_NO,
								 -1,
							  C_result) == ARM_OK) 
							  
		{
			FixedLegAmortId = C_result.getLong();
		}
		else
		{
			ERROR_MSG(C_result.getMsg(),pRet);
			return S_OK;
		}

		if (ARMLOCAL_ARM_GenAmortization(FixedLegAmortId,
										 amortMethodId,
										 freqAmortId,
										 amountAmort,
										 -1,
										 Ntl,
										 txAmort,
										 0.,
										 ARM_NULL_OBJECT,
										 0.0,
										 C_result) == ARM_OK) 
		{
			amortRefValue = C_result.getLong();
			//amortRefValueToDestroy = amortRefValue;
		}
		else
		{
			ERROR_MSG(C_result.getMsg(),pRet);
			return S_OK;
		}

		retCode = ARMLOCAL_FreeObject(FixedLegAmortId,C_result);
	}
	else
	{
		amortRefValue = LocalGetNumObjectId(l_typeAmort);

		// calcul du 1er nominal
		double dPaymentDate;

		if (ARMLOCAL_DisplayRefValue(amortRefValue,true,C_result) == ARM_OK)
		{
			dPaymentDate = C_result.getArray(0);
		}
		else
		{
			ERROR_MSG(C_result.getMsg(),pRet);
			return S_OK;
		}
		
		if (ARMLOCAL_CptRefValue(amortRefValue,dPaymentDate,C_result) == ARM_OK)
		{
			Ntl = C_result.getDouble();
		}
		else
		{
			ERROR_MSG(C_result.getMsg(),pRet);
			return S_OK;
		} 
	}
	//***************End Init Amort**********************

	//***************Init Funding************************
	double spreadFund;
	long   spreadType;
	long sizeMatu;

	if (VARIANT2Double(pSpreadFund,&spreadFund) == S_FALSE)
	{
		_bstr_t bSpreadFund(pSpreadFund);
		CCString l_SpreadFund= bSpreadFund;
		spreadFund = (double) LocalGetNumObjectId(l_SpreadFund);
		spreadType = 1;
	}
	else
	{
	   spreadType = 0;
	}

	_bstr_t indexFund(pIndexFund);
	CCString l_indexFund = indexFund;
	 long indexFundId;
	indexFundId = ARM_ConvIrType(l_indexFund);

	_bstr_t dayCountFund(pDayCountFund);
	CCString l_dayCountFund= dayCountFund;
	long dayCountFundId= ARM_ConvDayCount(l_dayCountFund);

	_bstr_t payFreqFund(pPayFreqFund);
	CCString l_payFreqFund= payFreqFund;
	long payFreqFundId;
	if((payFreqFundId = ARM_ConvFrequency(l_payFreqFund, C_result)) == ARM_DEFAULT_ERR)
	{
		ERROR_MSG("Invalid pay Frq Fund",pRet,ARM_ERROR_FREQ);
		return S_OK;
	}
	
	_bstr_t resetFreqFund(pResetFreqFund);
	CCString l_resetFreqFund = resetFreqFund;
	long resetFreqFundId;
	if((resetFreqFundId = ARM_ConvFrequency(l_resetFreqFund, C_result)) == ARM_DEFAULT_ERR)
	{
		ERROR_MSG("Invalid reset Frq Fund",pRet,ARM_ERROR_FREQ);
		return S_OK;
	}

	_bstr_t resetTimingFund(pResetTimingFund);
	CCString l_resetTimingFund = resetTimingFund;
	long resetTimingFundId = ARM_ConvPayResetRule (l_resetTimingFund);

	_bstr_t adjFund(pAdjFund);
	CCString l_adjFund= adjFund;
	long adjFundId = ARM_ConvIntRule(l_adjFund);

	//Création Funding Leg Refvalue 
	long irIndexFundId;
	if (ARMLOCAL_IRINDEX(KACTUAL_360,
						 payFreqFundId,
						 -1,
						 K_COMP_PROP,
						 K_MOD_FOLLOWING,
						 K_ADVANCE,
						 -2,
						 K_ARREARS,
						 10000.,
						 false,
						 l_ccy,
						 indexFundId,
						 K_COMP_PROP,
						 adjFundId,
						 resetFreqFundId,
						 C_result) == ARM_OK)
	{
		irIndexFundId = C_result.getLong();
	}
	else
	{
		ERROR_MSG(C_result.getMsg(),pRet);
		return S_OK;
	}

	long swapLegId;
	if (ARMLOCAL_SWAPLEG(irIndexFundId,
						 startDate,
						 endDate,
						 K_RCV,
						 spreadType, 
						 spreadFund, 
						 false,
						 l_ccy,
						 dayCountFundId,
						 -2,
						 "NULL",
						 "NULL",
						 0,
						 K_NX_NONE,
						 stubRuleId,
						 -1.0,
						 K_NO,
						  -1,
						 C_result) == ARM_OK) 
	{
		swapLegId = C_result.getLong();
	}
	else
	{
		ERROR_MSG(C_result.getMsg(),pRet);
		return S_OK;
	}
	
	long cloneId;
	if (ARMLOCAL_ClonedAndSetNotional(swapLegId,
									  amortRefValue,
									  100.,
									  C_result) == ARM_OK)
	{
		cloneId = C_result.getLong();
	}
	else
	{
		ERROR_MSG(C_result.getMsg(),pRet);
		return S_OK;
	}

	double pxFund; 
	double pxBumpFund;
	if (ARMLOCAL_ARM_Price(cloneId,
						   LocalGetNumObjectId(l_smiledModEUR),
						   C_result) == ARM_OK)
	{
		pxFund = C_result.getDouble();
	}
	else
	{
		ERROR_MSG(C_result.getMsg(),pRet);
		return S_OK;
	}


	#ifdef _DEBUG
		Res = ARMLOCAL_ARM_View(irIndexFundId,C_result);
	#endif

	#ifdef _DEBUG
		Res = ARMLOCAL_ARM_View(swapLegId,C_result);
	#endif

	#ifdef _DEBUG
		Res = ARMLOCAL_ARM_View(cloneId,C_result);
	#endif

	#ifdef _DEBUG
		Res = ARMLOCAL_ARM_View(LocalGetNumObjectId(l_smiledModEUR),C_result);
	#endif


	//Delta Funding

	if (ARMLOCAL_ARM_Price(cloneId,
						   LocalGetNumObjectId(l_deltaSmiledModEUR),
						   C_result) == ARM_OK)
	{
		pxBumpFund = C_result.getDouble();
	}
	else
	{
		ERROR_MSG(C_result.getMsg(),pRet);
		return S_OK;
	}
	
	fundingDeltaEUR=(pxBumpFund-pxFund);
	
	retCode = ARMLOCAL_FreeObject(irIndexFundId, C_result);
	retCode = ARMLOCAL_FreeObject(swapLegId, C_result);
	retCode = ARMLOCAL_FreeObject(cloneId,C_result);
	//***************End Init Funding********************

	//**************Init Phase Pre************************
	double spreadPhasePre;
	bool b_isSmoothedPhasePre=false;
	
	//Création Leg PhasePre Refvalue 
	long irIndexPhasePreId = ARM_NULL_OBJECT;
	long dayCountPhasePreId;

	

	double resetGapPhasePre; 

	if ( startDate != startDatePhaseStruct )
	{

		//Smoothed fixed index ?
		_bstr_t indexPhasePre(pIndexPhasePre);
		CCString l_indexPhasePre = indexPhasePre;

		if (l_indexPhasePre=="SMOOTHEDFIXED")
		{
			l_indexPhasePre="FIXED";
			spreadPhasePre = 0;
			b_isSmoothedPhasePre = true;
		}
		else
		{
			spreadPhasePre = pSpreadPhasePre;
			b_isSmoothedPhasePre = false;
		}
		
		long indexPhasePreId; 
		indexPhasePreId = ARM_ConvIrType(l_indexPhasePre); 

			
		_bstr_t dayCountPhasePre(pDayCountPhasePre);
		CCString l_dayCountPhasePre = dayCountPhasePre;
		dayCountPhasePreId= ARM_ConvDayCount(l_dayCountPhasePre);

		_bstr_t resetFreqPhasePre(pResetFreqPhasePre);
		CCString l_resetFreqPhasePre = resetFreqPhasePre;
		long resetFreqPhasePreId;
		if((resetFreqPhasePreId = ARM_ConvFrequency(l_resetFreqPhasePre, C_result)) == ARM_DEFAULT_ERR)
		{
			ERROR_MSG("Invalid reset Frq Pre",pRet,ARM_ERROR_FREQ);
			return S_OK;
		}

		_bstr_t payFreqPhasePre(pPayFreqPhasePre);
		CCString l_payFreqPhasePre = payFreqPhasePre;
  		long payFreqPhasePreId;	
		if((payFreqPhasePreId = ARM_ConvFrequency(l_payFreqPhasePre, C_result)) == ARM_DEFAULT_ERR)
		{
			ERROR_MSG("Invalid pay Frq Pre",pRet,ARM_ERROR_FREQ);
			return S_OK;
		}

		_bstr_t resetTimingPhasePre(pResetTimingPhasePre);
		CCString l_resetTimingPhasePre = resetTimingPhasePre;
		long resetTimingPhasePreId = ARM_ConvPayResetRule (l_resetTimingPhasePre);

		_bstr_t adjPhasePre(pAdjPhasePre);
		CCString l_adjPhasePre = adjPhasePre;    	
  		long adjPhasePreId = ARM_ConvIntRule(l_adjPhasePre);

		

		long swapLegPhasePreId;
		long clonePhasePreId;


		if (resetTimingPhasePreId == 1)
		{
			resetGapPhasePre = -2.; 
		}
		else
		{
			resetGapPhasePre = -15.; 
		}


		if (ARMLOCAL_IRINDEX(KACTUAL_360,
							 payFreqPhasePreId,
							 -1,
							 K_COMP_PROP,
							 K_MOD_FOLLOWING,
							 resetTimingPhasePreId,	
							 resetGapPhasePre,
							 K_ARREARS,
							 10000.,
							 false,
							 l_ccy,
							 indexPhasePreId,
							 K_COMP_PROP,
							 adjPhasePreId,
							 resetFreqPhasePreId,
							 C_result) == ARM_OK) 
		{
			irIndexPhasePreId = C_result.getLong();
		}
		else
		{
			ERROR_MSG(C_result.getMsg(),pRet);
			return S_OK;
		}

		if (b_isSmoothedPhasePre==false)
		{
		
			if (ARMLOCAL_SWAPLEG(irIndexPhasePreId,
								 startDate,
								 startDatePhaseStruct,
								 K_RCV,
								 0, 
								 spreadPhasePre, 
								 false,
								 l_ccy,
								 dayCountPhasePreId,
								 resetGapPhasePre,
								 "NULL",
								 "NULL",
								 0,
								 K_NX_NONE,
								 stubRuleId,
								 -1.0,
								 K_NO,
							  -1,
								 C_result) == ARM_OK) 
			{
				swapLegPhasePreId = C_result.getLong();
			}
			else
			{
				ERROR_MSG(C_result.getMsg(),pRet);
				return S_OK;
			}	

			if (ARMLOCAL_ClonedAndSetNotional(swapLegPhasePreId,
										  amortRefValue,
										  100.,
										  C_result) == ARM_OK)
			{
				clonePhasePreId = C_result.getLong();
			}
			else
			{
				ERROR_MSG(C_result.getMsg(),pRet);
				return S_OK;
			}

			if (ARMLOCAL_ARM_Price(clonePhasePreId,
								   LocalGetNumObjectId(l_smiledModEUR),
								   C_result) == ARM_OK)
			{
				pxPhasePre = C_result.getDouble();
			}
			else
			{
				ERROR_MSG(C_result.getMsg(),pRet);
				return S_OK;
			}

			//Delta Phase Pre

			if (ARMLOCAL_ARM_Price(clonePhasePreId,
								   LocalGetNumObjectId(l_deltaSmiledModEUR),
								   C_result) == ARM_OK)
			{
				pxBumpPhasePre= C_result.getDouble();
			}
			else
			{
				ERROR_MSG(C_result.getMsg(),pRet);
				return S_OK;
			}

			phasePreDeltaEUR = (pxBumpPhasePre-pxPhasePre);

			retCode = ARMLOCAL_FreeObject(swapLegPhasePreId,C_result);
			retCode = ARMLOCAL_FreeObject(clonePhasePreId,C_result);
		}
	}
	//*************End Phase Pre*****************************

	///************Init Phase struct***************************
	
  	_bstr_t dayCountPhaseStruct(pDayCountPhaseStruct);
	CCString l_dayCountPhaseStruct= dayCountPhaseStruct;
	long dayCountPhaseStructId= ARM_ConvDayCount(l_dayCountPhaseStruct);

  	_bstr_t payFreqPhaseStruct(pPayFreqPhaseStruct);
	CCString l_payFreqPhaseStruct= payFreqPhaseStruct;
	long payFreqPhaseStructId;
	if((payFreqPhaseStructId = ARM_ConvFrequency(l_payFreqPhaseStruct, C_result)) == ARM_DEFAULT_ERR)
	{
		ERROR_MSG("Invalid pay Frq Struct",pRet,ARM_ERROR_FREQ);
		return S_OK;
	}

	_bstr_t resetFreqPhaseStruct(pResetFreqPhaseStruct);
	CCString l_resetFreqPhaseStruct = resetFreqPhaseStruct;
	long resetFreqPhaseStructId;
	if((resetFreqPhaseStructId = ARM_ConvFrequency(l_resetFreqPhaseStruct, C_result)) == ARM_DEFAULT_ERR)
	{
		ERROR_MSG("Invalid reset Frq Struct",pRet,ARM_ERROR_FREQ);
		return S_OK;
	}

	_bstr_t resetTimingPhaseStruct(pResetTimingPhaseStruct);
	CCString l_resetTimingPhaseStruct = resetTimingPhaseStruct;
	long resetTimingPhaseStructId = ARM_ConvPayResetRule (l_resetTimingPhaseStruct);

	_bstr_t adjPhaseStruct(pAdjPhaseStruct);
	CCString l_adjPhaseStruct = adjPhaseStruct;
	long adjPhaseStructId = ARM_ConvIntRule(l_adjPhaseStruct);

	double resetGapPhaseStruct; 

	if (resetTimingPhaseStructId == 1)
	{
		resetGapPhaseStruct = -2.; 
	}
	else
	{
		resetGapPhaseStruct = -15.; 
	}

  	_bstr_t biscapped(pIsCapped);
	CCString l_iscapped = biscapped;

	long iscappedId;
	if ((iscappedId = ARM_ConvYesOrNo (l_iscapped, C_result)) == ARM_DEFAULT_ERR)
	{
	   ERROR_MSG("Invalid IsCapped",pRet);
	   return S_OK;
	}

	double TxCap = pTxCap;

	//Option
	//Foreign = EUR
	_bstr_t foreignccy(pForeignCcyPhaseStruct);
	CCString l_foreignccy = foreignccy;
	long foreignccyId; 
	if (ARMLOCAL_ISOCCY(l_foreignccy,C_result) == ARM_OK)
	{
		foreignccyId = C_result.getLong();
	}
	else
	{
		ERROR_MSG(C_result.getMsg(),pRet);
		return S_OK;
	}

  	//Domectic=CHF
	_bstr_t domesticccy(pDomesticCcyPhaseStruct);
	CCString l_domesticccy = domesticccy;
	long domesticccyId; 
	if (ARMLOCAL_ISOCCY(l_domesticccy,C_result) == ARM_OK)
	{
		domesticccyId = C_result.getLong();
	}
	else
	{
		ERROR_MSG(C_result.getMsg(),pRet);
		return S_OK;
	}

	double resetGapOptionPhaseStruct = pResetGapOptionPhaseStruct; 
	double payGapOptionPhaseStruct = pPayGapOptionPhaseStruct; 

	_bstr_t dayCountOptionPhaseStruct(pDayCountOptionPhaseStruct);
	CCString l_dayCountOptionPhaseStruct= dayCountOptionPhaseStruct;
	long dayCountOptionPhaseStructId= ARM_ConvDayCount(l_dayCountOptionPhaseStruct);

    _bstr_t payFreqOptionPhaseStruct(pPayFreqOptionPhaseStruct);
	CCString l_payFreqOptionPhaseStruct= payFreqOptionPhaseStruct;
	long payFreqOptionPhaseStructId;
	if((payFreqOptionPhaseStructId = ARM_ConvFrequency(l_payFreqOptionPhaseStruct, C_result)) == ARM_DEFAULT_ERR)
	{
		ERROR_MSG("Invalid pay Frq Struct Option",pRet,ARM_ERROR_FREQ);
		return S_OK;
	}

	_bstr_t resetFreqOptionPhaseStruct(pResetFreqOptionPhaseStruct);
	CCString l_resetFreqOptionPhaseStruct = resetFreqOptionPhaseStruct;
	long resetFreqOptionPhaseStructId;
	if((resetFreqOptionPhaseStructId = ARM_ConvFrequency(l_resetFreqOptionPhaseStruct, C_result)) == ARM_DEFAULT_ERR)
	{
		ERROR_MSG("Invalid reset Frq Struct Option",pRet,ARM_ERROR_FREQ);
		return S_OK;
	}
	
	_bstr_t resetTimingOptionPhaseStruct(pResetTimingOptionPhaseStruct);
	CCString l_resetTimingOptionPhaseStruct = resetTimingOptionPhaseStruct;
	long resetTimingOptionPhaseStructId = ARM_ConvPayResetRule (l_resetTimingOptionPhaseStruct);

	_bstr_t payTimingOptionPhaseStruct(pPayTimingOptionPhaseStruct);
	CCString l_payTimingOptionPhaseStruct = payTimingOptionPhaseStruct;
	long payTimingOptionPhaseStructId = ARM_ConvPayResetRule (l_payTimingOptionPhaseStruct);

  	_bstr_t adjOptionPhaseStruct(pAdjOptionPhaseStruct);
	CCString l_adjOptionPhaseStruct = adjOptionPhaseStruct;
	long adjOptionPhaseStructId = ARM_ConvIntRule(l_adjOptionPhaseStruct);

  	_bstr_t optionTypePhaseStruct(pOptionTypePhaseStruct);
	CCString l_optionTypePhaseStruct = optionTypePhaseStruct;
	long optionTypePhaseStructId;
	if((optionTypePhaseStructId = ARM_ConvCallOrPut (l_optionTypePhaseStruct, C_result)) == ARM_DEFAULT_ERR)
	{
		ERROR_MSG(C_result.getMsg(),pRet);
		return S_OK;
	}

  	_bstr_t strikePhaseStruct(pStrikePhaseStruct);
	CCString l_strikePhaseStruct= strikePhaseStruct;

	_bstr_t payOffDigitalPhaseStruct(pPayOffDigitalPhaseStruct);
	CCString l_payOffDigitalPhaseStruct= payOffDigitalPhaseStruct;
	
	double levier = pLevierPhaseStruct;
	double spot = pSpotPhaseStruct;

	
	////Création Leg PhaseStruct Refvalue 
	long irIndexPhaseStructId;
	
	if (ARMLOCAL_IRINDEX(KACTUAL_360,
						 payFreqPhaseStructId,
						 -1,
						 K_COMP_PROP,
						 K_MOD_FOLLOWING,
						 resetTimingPhaseStructId,	
						 resetGapPhaseStruct,
						 K_ARREARS,
						 10000.,
						 false,
						 l_ccy,
						 K_FIXED,
						 K_COMP_PROP,
						 adjPhaseStructId,
						 resetFreqPhaseStructId,
						 C_result) == ARM_OK) 
	{
		irIndexPhaseStructId = C_result.getLong();
	}
	else
	{
		ERROR_MSG(C_result.getMsg(),pRet);
		return S_OK;
	}

	//Construction de l'option
	long forexId;
	if (ARMLOCAL_FOREX(foreignccyId,
						domesticccyId,
						spot,
						C_result) == ARM_OK)
	{
		forexId = C_result.getLong();
	}
	else
	{
		ERROR_MSG(C_result.getMsg(),pRet);
		return S_OK;
	}

	long fxoptionstripId;
	if (ARMLOCAL_FxOptionStrip(XLDateToJulian(asOfDate),
							 forexId,
							 LocalGetNumObjectId(l_strikePhaseStruct),
							 optionTypePhaseStructId,
							 startDatePhaseStruct,
							 endDatePhaseStruct,
							 amortRefValue,
							 l_foreignccy,
							 payFreqOptionPhaseStructId,
							 dayCountOptionPhaseStructId,
							 l_domesticccy,
							 fwdRuleId,
							 adjOptionPhaseStructId,
							 stubRuleId,
							 resetGapOptionPhaseStruct,
							 payFreqOptionPhaseStructId, 
							 payGapOptionPhaseStruct,
							 l_foreignccy,
							 resetTimingOptionPhaseStructId,
							 payTimingOptionPhaseStructId,
							 "RCV",							 
							 ARM_NULL_OBJECT,
							 false,
							 1,
							 0.1,
							 ARM_NULL_OBJECT,
							 ARM_NULL_OBJECT,
							 levier,
							 C_result) == ARM_OK)
	//long ARMLOCAL_FxOptionStrip( double asOfDate,
	//						 long C_underlyingId,
	//						 long C_strikesCurveId,
	//						 long optionType,
	//						 double C_startDate,
	//						 double C_endDate,
	//						 long C_notionalId,
	//						 CCString C_paymentCcy,
	//						 long C_resetFreq,
	//						 long C_dayCount,
	//						 CCString C_resetCalendar,
	//						 long C_fwdRule,
	//						 long C_intRule,
	//						 long C_stubRule,
	//						 long C_resetGap,
	//						 long C_payFreq, 
	//						 long C_payGap, 
	//						 CCString C_payCalendar, 
	//						 long C_resetTiming, 
	//						 long C_payTiming,
	//						 CCString C_PorS,
	//						 long C_fxFixingsId,
	//						 bool isDigital,
	//						 int C_callSpreadFlag,
	//						 double C_epsilon,
	//						 long C_payoffCurveId,
	//						 long C_leverageId,
	//						 double C_leverageValue,
	//						 ARM_result& result,
	//						 long objId)
	{
		fxoptionstripId=C_result.getLong();
	}
	else
	{
		ERROR_MSG(C_result.getMsg(),pRet);
		return S_OK;
	}

	
	
	

	
	
	if (ARMLOCAL_ARM_Price(fxoptionstripId,
						   LocalGetNumObjectId(l_mixtureModEURCHF),
						   C_result) == ARM_OK)
	{
		pxPut = C_result.getDouble();
	}
	else
	{
		ERROR_MSG(C_result.getMsg(),pRet);
		return S_OK;
	}

	#ifdef _DEBUG
		Res = ARMLOCAL_ARM_View(fxoptionstripId,C_result);
	#endif

	#ifdef _DEBUG
		Res = ARMLOCAL_ARM_View(LocalGetNumObjectId(l_mixtureModEURCHF),C_result);
	#endif

	if (ARMLOCAL_ARM_Price(fxoptionstripId,
						   LocalGetNumObjectId(l_deltaEURMixtureMod),
						   C_result) == ARM_OK)
	{
		pxPutBumpDeltaEUR = C_result.getDouble();
	}
	else
	{
		ERROR_MSG(C_result.getMsg(),pRet);
		return S_OK;
	}

	if (ARMLOCAL_ARM_Price(fxoptionstripId,
						   LocalGetNumObjectId(l_deltaCHFMixtureMod),
						   C_result) == ARM_OK)
	{
		pxPutBumpDeltaCHF = C_result.getDouble();
	}
	else
	{
		ERROR_MSG(C_result.getMsg(),pRet);
		return S_OK;
	}


	putDeltaEUR=(pxPutBumpDeltaEUR-pxPut);
	putDeltaCHF=(pxPutBumpDeltaCHF-pxPut);

	if (ARMLOCAL_ARM_Price(fxoptionstripId,
						   LocalGetNumObjectId(l_bSEURMixtureMod),
						   C_result) == ARM_OK)
	{
		pxPutBumpBSEUR = C_result.getDouble();
	}
	else
	{
		ERROR_MSG(C_result.getMsg(),pRet);
		return S_OK;
	}

	if (ARMLOCAL_ARM_Price(fxoptionstripId,
						   LocalGetNumObjectId(l_bSCHFMixtureMod),
						   C_result) == ARM_OK)
	{
		pxPutBumpBSCHF = C_result.getDouble();
	}
	else
	{
		ERROR_MSG(C_result.getMsg(),pRet);
		return S_OK;
	}


	putBSEUR=pxPutBumpBSEUR-pxPut;
	putBSCHF=pxPutBumpBSCHF-pxPut;

	if (ARMLOCAL_ARM_Price(fxoptionstripId,
						   LocalGetNumObjectId(l_volFxEURMixtureMod),
						   C_result) == ARM_OK)
	{
		pxPutBumpFxVolEUR = C_result.getDouble();
	}
	else
	{
		ERROR_MSG(C_result.getMsg(),pRet);
		return S_OK;
	}

	putVolFxEUR = pxPutBumpFxVolEUR- pxPut;

	
	//suppression des objets qui ne sont plus utilisés
	retCode = ARMLOCAL_FreeObject(fxoptionstripId,C_result);


	long fxoptionstripdigitalId;
	if (helvetixType==3)
	{	
		if (ARMLOCAL_FxOptionStrip(XLDateToJulian(asOfDate),
								 forexId,
								 LocalGetNumObjectId(l_strikePhaseStruct),
								 optionTypePhaseStructId,
								 startDatePhaseStruct,
								 endDatePhaseStruct,
								 amortRefValue,
								 l_foreignccy,
								 payFreqOptionPhaseStructId,
								 dayCountOptionPhaseStructId,
								 l_domesticccy,
								 fwdRuleId,
								 adjOptionPhaseStructId,
								 stubRuleId,
								 resetGapOptionPhaseStruct,
								 payFreqOptionPhaseStructId, 
								 payGapOptionPhaseStruct,
								 l_foreignccy,
								 resetTimingOptionPhaseStructId,
								 payTimingOptionPhaseStructId,
								 "RCV",							 
								 ARM_NULL_OBJECT,
								 true,
								 1,
								 0.0001,
								 LocalGetNumObjectId(l_payOffDigitalPhaseStruct),
								 ARM_NULL_OBJECT,
								 levier,
								 C_result) == ARM_OK)
		{
			fxoptionstripdigitalId=C_result.getLong();
		}
		else
		{
			ERROR_MSG(C_result.getMsg(),pRet);
			return S_OK;
		}

		if (ARMLOCAL_ARM_Price(fxoptionstripdigitalId,
							   LocalGetNumObjectId(l_mixtureModEURCHF),
							   C_result) == ARM_OK)
		{
			pxDigital = C_result.getDouble();
		}
		else
		{
			ERROR_MSG(C_result.getMsg(),pRet);
			return S_OK;
		}

		#ifdef _DEBUG
			Res = ARMLOCAL_ARM_View(fxoptionstripdigitalId,C_result);
		#endif

	
		if (ARMLOCAL_ARM_Price(fxoptionstripdigitalId,
							   LocalGetNumObjectId(l_deltaEURMixtureMod),
							   C_result) == ARM_OK)
		{
			pxDigitalBumpDeltaEUR = C_result.getDouble();
		}
		else
		{
			ERROR_MSG(C_result.getMsg(),pRet);
			return S_OK;
		}

		if (ARMLOCAL_ARM_Price(fxoptionstripdigitalId,
							   LocalGetNumObjectId(l_deltaCHFMixtureMod),
							   C_result) == ARM_OK)
		{
			pxDigitalBumpDeltaCHF = C_result.getDouble();
		}
		else
		{
			ERROR_MSG(C_result.getMsg(),pRet);
			return S_OK;
		}


		digitalDeltaEUR=(pxDigitalBumpDeltaEUR-pxDigital);
		digitalDeltaCHF=(pxDigitalBumpDeltaCHF-pxDigital);

		if (ARMLOCAL_ARM_Price(fxoptionstripdigitalId,
							   LocalGetNumObjectId(l_bSEURMixtureMod),
							   C_result) == ARM_OK)
		{
			pxDigitalBumpBSEUR = C_result.getDouble();
		}
		else
		{
			ERROR_MSG(C_result.getMsg(),pRet);
			return S_OK;
		}

		if (ARMLOCAL_ARM_Price(fxoptionstripdigitalId,
							   LocalGetNumObjectId(l_bSCHFMixtureMod),
							   C_result) == ARM_OK)
		{
			pxDigitalBumpBSCHF = C_result.getDouble();
		}
		else
		{
			ERROR_MSG(C_result.getMsg(),pRet);
			return S_OK;
		}


		digitalBSEUR=pxDigitalBumpBSEUR-pxDigital;
		digitalBSCHF=pxDigitalBumpBSCHF-pxDigital;

		if (ARMLOCAL_ARM_Price(fxoptionstripdigitalId,
							   LocalGetNumObjectId(l_volFxEURMixtureMod),
							   C_result) == ARM_OK)
		{
			pxDigitalBumpFxVolEUR = C_result.getDouble();
		}
		else
		{
			ERROR_MSG(C_result.getMsg(),pRet);
			return S_OK;
		}

		digitalVolFxEUR = pxDigitalBumpFxVolEUR- pxDigital;

		

		//Calcul Disc = > certain leg (fixed leg)
		//On re créer le fixoptionstrip en décomposant pour chaque échéance, en un put long et un put court
		//afin de recréer une digitale
		//Le calcul de la marge est fait dans une fonction séparée

		VECTOR<double> vecProba;
		VECTOR<double> vecMT;
		long ProbaSize;
		long MTSize;

		if(VARIANT2VECTORDOUBLE (*pProba,vecProba,ProbaSize) != S_OK)
			return S_FALSE;

		if(VARIANT2VECTORDOUBLE (*pMT,vecMT,MTSize) != S_OK)
			return S_FALSE;

		discMargin = computeHelvetixIIIDiscMargin(forexId,
												startDatePhaseStruct,
												endDatePhaseStruct,
												spot,
												l_domesticccy,
												l_foreignccy,
												optionTypePhaseStructId,
												LocalGetNumObjectId(l_strikePhaseStruct),
												0.0001,
												LocalGetNumObjectId(l_payOffDigitalPhaseStruct),
												amortRefValue,
												levier,
												payFreqOptionPhaseStructId,
												dayCountOptionPhaseStructId,
												resetTimingOptionPhaseStructId,
												resetGapOptionPhaseStruct,
												adjOptionPhaseStructId,
												fwdRuleId,
												vecProba,
												vecMT,
												l_smiledModCHF,
												l_mixtureModEURCHF);

	
		retCode = ARMLOCAL_FreeObject(fxoptionstripdigitalId,C_result);
	
	}


	//*************End Init Phase Struct************************

	//*************Init Phase Post***************************

	bool b_isSmoothedPhasePost=false;
	long irIndexPhasePostId;
	long swapLegPhasePostId;
	long clonePhasePostId;
	long dayCountPhasePostId;
	double spreadPhasePost;
	double resetGapPhasePost; 

	if ( endDatePhaseStruct != endDate )
	{		

		_bstr_t indexPhasePost(pIndexPhasePost);
		CCString l_indexPhasePost = indexPhasePost;

		//Smoothed fixed index ?
		if (l_indexPhasePost=="SMOOTHEDFIXED")
		{
			l_indexPhasePost="FIXED";
			spreadPhasePost = 0;
			b_isSmoothedPhasePost = true;
		}
		else
		{
			spreadPhasePost = pSpreadPhasePost;
			b_isSmoothedPhasePost = false;
		}
		
		
  		long indexPhasePostId;
		indexPhasePostId = ARM_ConvIrType(l_indexPhasePost);
		
		_bstr_t dayCountPhasePost(pDayCountPhasePost);
		CCString l_dayCountPhasePost = dayCountPhasePost;
		dayCountPhasePostId= ARM_ConvDayCount(l_dayCountPhasePost);

		_bstr_t payFreqPhasePost(pPayFreqPhasePost);
		CCString l_payFreqPhasePost = payFreqPhasePost;
		long payFreqPhasePostId;
		if((payFreqPhasePostId = ARM_ConvFrequency(l_payFreqPhasePost, C_result)) == ARM_DEFAULT_ERR)
		{
			ERROR_MSG("Invalid pay Frq Post",pRet,ARM_ERROR_FREQ);
			return S_OK;
		}

		_bstr_t resetFreqPhasePost(pResetFreqPhasePost);
		CCString l_resetFreqPhasePost = resetFreqPhasePost;
		long resetFreqPhasePostId;
		if((resetFreqPhasePostId = ARM_ConvFrequency(l_resetFreqPhasePost, C_result)) == ARM_DEFAULT_ERR)
		{
			ERROR_MSG("Invalid reset Frq Post",pRet,ARM_ERROR_FREQ);
			return S_OK;
		}

		_bstr_t resetTimingPhasePost(pResetTimingPhasePost);
		CCString l_resetTimingPhasePost = resetTimingPhasePost;
		long resetTimingPhasePostId = ARM_ConvPayResetRule (l_resetTimingPhasePost);

		_bstr_t adjPhasePost(pAdjPhasePost);
		CCString l_adjPhasePost = adjPhasePost;
		long adjPhasePostId = ARM_ConvIntRule(l_adjPhasePost);

  		


		if (resetTimingPhasePostId == 1)
		{
			resetGapPhasePost = -2.; 
		}
		else
		{
			resetGapPhasePost = -15.; 
		}


		if (ARMLOCAL_IRINDEX(KACTUAL_360,
							 payFreqPhasePostId,
							 -1,
							 K_COMP_PROP,
							 K_MOD_FOLLOWING,
							 resetTimingPhasePostId,	
							 resetGapPhasePost,
							 K_ARREARS,
							 10000.,
							 false,
							 l_ccy,
							 indexPhasePostId,
							 K_COMP_PROP,
							 adjPhasePostId,
							 resetFreqPhasePostId,
							 C_result) == ARM_OK) 
		{
			irIndexPhasePostId = C_result.getLong();
		}
		else
		{
			ERROR_MSG(C_result.getMsg(),pRet);
			return S_OK;
		}
		
		if (b_isSmoothedPhasePost==false)
		  {
			if (ARMLOCAL_SWAPLEG(irIndexPhasePostId,
										 endDatePhaseStruct,
										 endDate,
										 K_RCV,
										 0, 
										 spreadPhasePost, 
										 false,
										 l_ccy,
										 dayCountPhasePostId,
										 resetGapPhasePost,
										 "NULL",
										 "NULL",
										 0,
										 K_NX_NONE,
										 stubRuleId,
										 -1.0,
										 K_NO,
									  -1,
										 C_result) == ARM_OK) 
			{
				swapLegPhasePostId = C_result.getLong();
			}
			else
			{
				ERROR_MSG(C_result.getMsg(),pRet);
				return S_OK;
			}
					
					
			if (ARMLOCAL_ClonedAndSetNotional(swapLegPhasePostId,
											  amortRefValue,
											  100.,
											  C_result) == ARM_OK)
			{
				clonePhasePostId = C_result.getLong();
			}
			else
			{
				ERROR_MSG(C_result.getMsg(),pRet);
				return S_OK;
			}


			if (ARMLOCAL_ARM_Price(clonePhasePostId,
								   LocalGetNumObjectId(l_smiledModEUR),
								   C_result) == ARM_OK)
			{
				pxPhasePost = C_result.getDouble();
			}
			else
			{
				ERROR_MSG(C_result.getMsg(),pRet);
				return S_OK;
			}

			//Delta Phase Post
			if (ARMLOCAL_ARM_Price(clonePhasePostId,
								   LocalGetNumObjectId(l_deltaSmiledModEUR),
								   C_result) == ARM_OK)
			{
				pxBumpPhasePost = C_result.getDouble();
			}
			else
			{
				ERROR_MSG(C_result.getMsg(),pRet);
				return S_OK;
			}
		
			phasePostDeltaEUR=(pxBumpPhasePost-pxPhasePost);

			retCode = ARMLOCAL_FreeObject(swapLegPhasePostId,C_result);
			retCode = ARMLOCAL_FreeObject(clonePhasePostId,C_result);
		 }
	}
	

	//*************End Init Phase Post ************************
	

	//******************** Calcul premiere PV ******************
	int iteration =1;
	double prix=0;
	double pxPhaseStruct=0;
	double pxCap=0;
	double pxTotal=0;
	double prixAct=0; 
	double prixAct2=0;
	double SensiTotal=0;
	double Xincr=0;
	double prixInter=0;
	double SprSeek = 10.01;
	double SprSeekFinal= 10.01;
	double SprSeekBonifie = 10.01;  // PUREMENT ARBITRAIRE
	double DeltaEUR=0;
	double DeltaCHF=0;
	double BSEUR=0;
	double BSCHF=0;
	double VolFxEUR=0;
	
	fee = Ntl * fee / 100 ;
	//******************** Calcul premiere PV ******************
	
	VECTOR<double> pxHelvetix = computeHelvetix(SprSeekBonifie,
													asOfDate,
													startDate,
													startDatePhaseStruct,
													endDatePhaseStruct,
													endDate,
													helvetixType,
													spot,
													amortRefValue,
													l_ccy,
													stubRuleId,
													b_isSmoothedPhasePre,
													irIndexPhasePreId,
													dayCountPhasePreId,
													resetGapPhasePre,
													b_isSmoothedPhasePost,
													irIndexPhasePostId,
													dayCountPhasePostId,
													resetGapPhasePost,
													dayCountPhaseStructId,
													resetGapPhaseStruct,
													irIndexPhaseStructId,
													iscappedId,
													TxCap,
													l_strikePhaseStruct,
													forexId,
													optionTypePhaseStructId,
													l_foreignccy,
													l_domesticccy,
													payFreqOptionPhaseStructId,
													dayCountOptionPhaseStructId,
													fwdRuleId,
													adjOptionPhaseStructId,
													resetGapOptionPhaseStruct,
													payGapOptionPhaseStruct,
													resetTimingOptionPhaseStructId,
													payTimingOptionPhaseStructId,
													levier,
													l_smiledModEUR,
													l_deltaSmiledModEUR,
													l_mixtureModEURCHF,
													l_deltaEURMixtureMod,
													l_bSEURMixtureMod,
													l_volFxEURMixtureMod,
													l_deltaCHFMixtureMod,
													l_bSCHFMixtureMod);

	if (pxHelvetix.size() == 0)					
	{
		ERROR_MSG("Pb in computepxHelvetix",pRet);
		return S_OK;
	}

	if (b_isSmoothedPhasePre==true)
	{
		pxPhasePre = pxHelvetix[0];
		phasePreDeltaEUR = pxHelvetix[4];
	}

	pxPhaseStruct = pxHelvetix[1];
	phaseStructDeltaEUR = pxHelvetix[5];

	if (iscappedId==K_YES)
	{
		pxCap =  pxHelvetix[2];
		capDeltaEUR=pxHelvetix[7];
		capDeltaCHF=pxHelvetix[8];
		capBSEUR=pxHelvetix[9];
		capBSCHF=pxHelvetix[10];
		capVolFxEUR=pxHelvetix[11];
	}


	if (b_isSmoothedPhasePost==true)
	{
		pxPhasePost = pxHelvetix[3];
		phasePostDeltaEUR= pxHelvetix[6];
	}

	pxTotal = pxPhasePre + pxPhaseStruct + pxPut + pxDigital + pxPhasePost - pxCap - pxFund;

	DeltaEUR = phasePreDeltaEUR + phasePostDeltaEUR + phaseStructDeltaEUR + putDeltaEUR + digitalDeltaEUR - capDeltaEUR - fundingDeltaEUR;
	DeltaCHF = putDeltaCHF + digitalDeltaCHF - capDeltaCHF;
	BSEUR = putBSEUR + digitalBSEUR - capBSEUR;
	BSCHF = putBSCHF + digitalBSCHF - capBSCHF;
	VolFxEUR = putVolFxEUR+ digitalVolFxEUR - capVolFxEUR;

	SensiTotal = abs(DeltaEUR) + abs(DeltaCHF) + abs(BSEUR) + abs(BSCHF) + abs(VolFxEUR) + abs(discMargin);

	
	 //En attendant
	SensiTotal = 0;

	prixAct = pxTotal - SensiTotal - fee; 

	prixAct2 =1;
	
	//******************** BOUCLE SUR LE PRICE ******************
	while ((abs(prixAct2) > 0.0001) && (iteration < 10 ))
	{
		SprSeek = SprSeekBonifie;
		
		//*************************
		// **** Premier recalcul ****
		//****************************

		pxHelvetix = computeHelvetix(SprSeekBonifie,
										asOfDate,
										startDate,
										startDatePhaseStruct,
										endDatePhaseStruct,
										endDate,
										helvetixType,
										spot,
										amortRefValue,
										l_ccy,
										stubRuleId,
										b_isSmoothedPhasePre,
										irIndexPhasePreId,
										dayCountPhasePreId,
										resetGapPhasePre,
										b_isSmoothedPhasePost,
										irIndexPhasePostId,
										dayCountPhasePostId,
										resetGapPhasePost,
										dayCountPhaseStructId,
										resetGapPhaseStruct,
										irIndexPhaseStructId,
										iscappedId,
										TxCap,
										l_strikePhaseStruct,
										forexId,
										optionTypePhaseStructId,
										l_foreignccy,
										l_domesticccy,
										payFreqOptionPhaseStructId,
										dayCountOptionPhaseStructId,
										fwdRuleId,
										adjOptionPhaseStructId,
										resetGapOptionPhaseStruct,
										payGapOptionPhaseStruct,
										resetTimingOptionPhaseStructId,
										payTimingOptionPhaseStructId,
										levier,
										l_smiledModEUR,
										l_deltaSmiledModEUR,
										l_mixtureModEURCHF,
										l_deltaEURMixtureMod,
										l_bSEURMixtureMod,
										l_volFxEURMixtureMod,
										l_deltaCHFMixtureMod,
										l_bSCHFMixtureMod);

		if (pxHelvetix.size() == 0)					
		{
			ERROR_MSG("Pb in computepxHelvetix",pRet);
			return S_OK;
		}

		if (b_isSmoothedPhasePre==true)
		{
			pxPhasePre = pxHelvetix[0];
			phasePreDeltaEUR = pxHelvetix[4];
		}

		pxPhaseStruct = pxHelvetix[1];
		phaseStructDeltaEUR = pxHelvetix[5];

		if (iscappedId==K_YES)
		{
			pxCap =  pxHelvetix[2];
			capDeltaEUR=pxHelvetix[7];
			capDeltaCHF=pxHelvetix[8];
			capBSEUR=pxHelvetix[9];
			capBSCHF=pxHelvetix[10];
			capVolFxEUR=pxHelvetix[11];
		}


		if (b_isSmoothedPhasePost==true)
		{
			pxPhasePost = pxHelvetix[3];
			phasePostDeltaEUR= pxHelvetix[6];
		}

		pxTotal = pxPhasePre + pxPhaseStruct + pxPut + pxDigital + pxPhasePost - pxCap - pxFund;

		DeltaEUR = phasePreDeltaEUR + phasePostDeltaEUR + phaseStructDeltaEUR + putDeltaEUR + digitalDeltaEUR - capDeltaEUR - fundingDeltaEUR;
		DeltaCHF = putDeltaCHF + digitalDeltaCHF - capDeltaCHF;
		BSEUR = putBSEUR + digitalBSEUR - capBSEUR;
		BSCHF = putBSCHF + digitalBSCHF - capBSCHF;
		VolFxEUR = putVolFxEUR+ digitalVolFxEUR - capVolFxEUR;

		SensiTotal = abs(DeltaEUR) + abs(DeltaCHF) + abs(BSEUR) + abs(BSCHF) + abs(VolFxEUR)  + abs(discMargin);

		prixAct = pxTotal - SensiTotal - fee; 
			

		//******************************
		// **** affectation valeurs ****
		//******************************

		prixAct2 = prixAct ;
		Xincr = SprSeek / 1000000;
		SprSeekBonifie = SprSeek + Xincr;		

		
		//****************************
		// **** Deuxieme recalcul ****
		//****************************
		pxHelvetix = computeHelvetix(SprSeekBonifie,
										asOfDate,
										startDate,
										startDatePhaseStruct,
										endDatePhaseStruct,
										endDate,
										helvetixType,
										spot,
										amortRefValue,
										l_ccy,
										stubRuleId,
										b_isSmoothedPhasePre,
										irIndexPhasePreId,
										dayCountPhasePreId,
										resetGapPhasePre,
										b_isSmoothedPhasePost,
										irIndexPhasePostId,
										dayCountPhasePostId,
										resetGapPhasePost,
										dayCountPhaseStructId,
										resetGapPhaseStruct,
										irIndexPhaseStructId,
										iscappedId,
										TxCap,
										l_strikePhaseStruct,
										forexId,
										optionTypePhaseStructId,
										l_foreignccy,
										l_domesticccy,
										payFreqOptionPhaseStructId,
										dayCountOptionPhaseStructId,
										fwdRuleId,
										adjOptionPhaseStructId,
										resetGapOptionPhaseStruct,
										payGapOptionPhaseStruct,
										resetTimingOptionPhaseStructId,
										payTimingOptionPhaseStructId,
										levier,
										l_smiledModEUR,
										l_deltaSmiledModEUR,
										l_mixtureModEURCHF,
										l_deltaEURMixtureMod,
										l_bSEURMixtureMod,
										l_volFxEURMixtureMod,
										l_deltaCHFMixtureMod,
										l_bSCHFMixtureMod);


			if (pxHelvetix.size() == 0)					
		{
			ERROR_MSG("Pb in computepxHelvetix",pRet);
			return S_OK;
		}

		if (b_isSmoothedPhasePre==true)
		{
			pxPhasePre = pxHelvetix[0];
			phasePreDeltaEUR = pxHelvetix[4];
		}

		pxPhaseStruct = pxHelvetix[1];
		phaseStructDeltaEUR = pxHelvetix[5];

		if (iscappedId==K_YES)
		{
			pxCap =  pxHelvetix[2];
			capDeltaEUR=pxHelvetix[7];
			capDeltaCHF=pxHelvetix[8];
			capBSEUR=pxHelvetix[9];
			capBSCHF=pxHelvetix[10];
			capVolFxEUR=pxHelvetix[11];
		}


		if (b_isSmoothedPhasePost==true)
		{
			pxPhasePost = pxHelvetix[3];
			phasePostDeltaEUR= pxHelvetix[6];
		}

		pxTotal = pxPhasePre + pxPhaseStruct + pxPut + pxDigital + pxPhasePost - pxCap - pxFund;

		DeltaEUR = phasePreDeltaEUR + phasePostDeltaEUR + phaseStructDeltaEUR + putDeltaEUR + digitalDeltaEUR - capDeltaEUR - fundingDeltaEUR;
		DeltaCHF = putDeltaCHF + digitalDeltaCHF - capDeltaCHF;
		BSEUR = putBSEUR + digitalBSEUR - capBSEUR;
		BSCHF = putBSCHF + digitalBSCHF - capBSCHF;
		VolFxEUR = putVolFxEUR+ digitalVolFxEUR - capVolFxEUR;

		SensiTotal = abs(DeltaEUR) + abs(DeltaCHF) + abs(BSEUR) + abs(BSCHF) + abs(VolFxEUR) + abs(discMargin);
	
		prixAct = pxTotal - SensiTotal - fee; 	

		iteration++;
		
		SprSeekFinal = SprSeekBonifie;
		prixInter = ( prixAct - prixAct2 ) / Xincr;
		SprSeekBonifie = SprSeek - prixAct / prixInter ; 
		prixAct2 = prixAct;
	}
		
//**************** PRICE *******************

	//suppression de tous les objets construits...
	retCode = ARMLOCAL_FreeObject(ccyId,C_result);
	retCode = ARMLOCAL_FreeObject(irIndexPhasePostId,C_result);
	retCode = ARMLOCAL_FreeObject(irIndexPhasePreId,C_result);
	retCode = ARMLOCAL_FreeObject(irIndexPhaseStructId,C_result);
	retCode = ARMLOCAL_FreeObject(forexId,C_result);


	//Test de la valeur de retour - no convergence
	if ((abs(prixAct2) > 0.1) && (iteration > 9 ))
	{
		SprSeekFinal=-999;
	}

	VECTOR<double> vRes;

	vRes.push_back(SprSeekFinal);


	Res = VECTORDOUBLE2VARIANT(vRes, pRet);

	_variant_t wrap_Res;
	wrap_Res.Attach(*pRet);
	*pRet=wrap_Res.Detach();

	return S_OK;
}
catch(...)
{
	return E_FAIL;
}
}




