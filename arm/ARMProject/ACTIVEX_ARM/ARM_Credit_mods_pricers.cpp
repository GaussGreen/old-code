// ARMModule.cpp : Implementation of CLocal_DLLARMApp and DLL registration.
#include "firsttobeincluded.h"
#include "ActiveXModule.h"

#include <ARM/libarm_frometk/VariantTools.h>
#include "CCatl.h"
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
#include "ICM_local_summit.h"
#include "ARM_local_interglob.h"



STDMETHODIMP ActiveXModule::ARM_Credit_DefProbModelNew(BSTR pDefCurve,
												   BSTR pIRcurve,
												   BSTR pVolcurve,
												   VARIANT *pRet)
{
	// TODO: Add your implementation code here
try
{
	ARM_result C_result;

	_bstr_t b_IRcurve(pIRcurve);
	CCString C_IRcurve = b_IRcurve;

	long ircurve = LocalGetNumObjectId( C_IRcurve);

	if (ircurve == ARM_KO)
	{
	ERROR_MSG("Invalid IRcurve",pRet,ARM_ERROR_IRCURVE);
	return S_OK;
	}

	_bstr_t b_DefCurve(pDefCurve);
	CCString C_Defcurve = b_DefCurve;

	long defcurve = LocalGetNumObjectId(C_Defcurve);

	_bstr_t b_Volcurve(pVolcurve);
	CCString C_Volcurve = b_Volcurve;

	long Volcurve = LocalGetNumObjectId(C_Volcurve);

	CCString Res;
	CCString curClass = LOCAL_DEFPROBMODEL_CLASS;

	if(ICMLOCAL_DefProbModel(defcurve, 
							 ircurve,
							 Volcurve,
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

STDMETHODIMP ActiveXModule::ARM_Credit_ModelMultiCurves(VARIANT *pIRcurve,
													VARIANT *pDefCurves, 
													VARIANT *pRecoveryRates, 
													BSTR	CorrelId,
													BSTR	pVolcurve,
													BSTR	pCpnInfcurve,
													BSTR	pCpnIRcurve,
													VARIANT *pRet)
{
	// TODO: Add your implementation code here
try
{
	ARM_result C_result;

	CCString C_IRcurve;

	_bstr_t b_CorrelId(CorrelId);
	CCString C_CorrelId = b_CorrelId;

	if(VARIANT2CCString (*pIRcurve, C_IRcurve) != S_OK)
		return S_FALSE;

	long ircurve = LocalGetNumObjectId( C_IRcurve);

	if (ircurve == ARM_KO)
	{
	ERROR_MSG("Invalid IRcurve",pRet,ARM_ERROR_IRCURVE);
	return S_OK;
	}

	_bstr_t b_Volcurve(pVolcurve);
	CCString C_Volcurve = b_Volcurve;

	long Volcurve = LocalGetNumObjectId(C_Volcurve);

	_bstr_t b_CpnInfcurve(pCpnInfcurve);
	CCString C_CpnInfcurve = b_CpnInfcurve;

	long CpnInfcurve = LocalGetNumObjectId(C_CpnInfcurve);

	_bstr_t b_CpnIRcurve(pCpnIRcurve);
	CCString C_CpnIRcurve = b_CpnIRcurve;

	long CpnIRcurve = LocalGetNumObjectId(C_CpnIRcurve);

	VECTOR<CCString> C_DefCurves;
	VECTOR<long> l_DefCurves;

	long size_DefCurves;

	if(VARIANT2VECTORCCSTRING(*pDefCurves, C_DefCurves,size_DefCurves) != S_OK)
		return S_FALSE;

	long defcurve =0;

	for (int i=0; i<size_DefCurves; i++)
	{
		defcurve = LocalGetNumObjectId(C_DefCurves[i]);

		if (defcurve == ARM_KO)
		{
		ERROR_MSG("Invalid DefCurve",pRet,ARM_ERROR_DFCURVE);
		return S_OK;
		}

		l_DefCurves.push_back(defcurve);
	}

	VECTOR<double> C_RecoveryRates;
	long size_RecoveryRates;

	if(VARIANT2VECTORDOUBLE (*pRecoveryRates, C_RecoveryRates,size_RecoveryRates) != S_OK)
		return S_FALSE;

	CCString Res;
	CCString curClass = LOCAL_MULTICURVESMODEL_CLASS;

	if(ICMLOCAL_ModelMultiCurves (size_DefCurves, 
								l_DefCurves,
								ircurve,
								C_RecoveryRates,
								LocalGetNumObjectId(C_CorrelId),
								Volcurve,
								true,
								CpnInfcurve,
								CpnIRcurve,
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
catch(...)
{
	return E_FAIL;
}
}


STDMETHODIMP ActiveXModule::ARM_Credit_ModelMultiCvMktDataMng(VARIANT *pIRcurve,
													VARIANT *pDefCurves, 
													VARIANT *pRecoveryRates, 
													BSTR	CorrelId,
													BSTR	MktDataMngId,
													BSTR	pVolcurve,
													BSTR	cloneOrNot,
													VARIANT *pRet)
{
	ARM_result C_result;
	try {
		// DefCurves
		vector<string> vLabelStr;
		VariantTools::convert(*pDefCurves, vLabelStr);
		int size_DefCurves = vLabelStr.size();
		vector<long> vLabelId(size_DefCurves);
		for (int i=0; i<size_DefCurves; i++)
		{
			long defcurve = LocalGetNumObjectId(CCString(vLabelStr[i].c_str()));
			if (defcurve == ARM_KO)
			{
				ICMTHROW(ERR_INVALID_ARGUMENT, "ARM_Credit_ModelMultiCvMktDataMng : Bad defCurve");
			}
			vLabelId[i] = defcurve;
		}
		// irCurves
		string irCurveStr("");
		VariantTools::convert(pIRcurve, irCurveStr);
		long ircurve = LocalGetNumObjectId( CCString(irCurveStr.c_str()));
		// recov
		vector<double> recovV;
		VariantTools::convert(*pRecoveryRates, recovV);
		string StrcorrelId;
		VariantTools::convert(CorrelId, StrcorrelId);
		long LongcorrelId = LocalGetNumObjectId( CCString(StrcorrelId.c_str()));
		// MktDataMng
		string MktDataMngStr;
		VariantTools::convert(MktDataMngId, MktDataMngStr);
		long MktDataMngId = LocalGetNumObjectId( CCString(MktDataMngStr.c_str()));
		// VolCurve
		string VolCurveStr;
		VariantTools::convert(pVolcurve, VolCurveStr);
		long VoCurveId = LocalGetNumObjectId( CCString(VolCurveStr.c_str()));
		// cloneOrNot
		string cloneOrNotStr;
		VariantTools::convert(cloneOrNot, cloneOrNotStr);
		bool cloneOrNotBool = false;
		cloneOrNotBool = ARM_ConvYesOrNo(cloneOrNotStr);
		

		long  prevId = -1;
		long newId = ICMLOCAL_ModelMultiCurves (size_DefCurves, 
								vLabelId,
								ircurve,
								recovV,
								LongcorrelId,
								MktDataMngId,
								VoCurveId,
								cloneOrNotBool,
								prevId);
		string newLabel= LocalPersistent::get().getStringId(newId, LOCAL_MULTICURVESMODEL_CLASS);
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

STDMETHODIMP ActiveXModule::ARM_Credit_Pricer(BSTR pSecurity, 
										  BSTR pModel,
										  BSTR pPricerType,
										  int		l_nbpaths,
										  BSTR     pParameters,
										  double   AsOfDate,
										  VARIANT *pRet)
{

try
{

	ARM_result C_result;

	CCString Res;

	_bstr_t b_Security(pSecurity);
	CCString C_Security = b_Security;

	_bstr_t b_Model(pModel);
	CCString C_Model = b_Model;

	_bstr_t b_PricerType(pPricerType);
	CCString C_PricerType = b_PricerType;

	_bstr_t b_Parameters(pParameters);
	CCString C_Parameters = b_Parameters;

 
	CCString curClass = LOCAL_PRICER_CLASS;

	int l_PricerType;

	if((l_PricerType = ARM_ConvCreditPricerType(C_PricerType, C_result)) == ARM_DEFAULT_ERR)
	{
		ERROR_MSG("Invalid Pricer Type",pRet,ARM_ERROR_FREQ);
		return S_OK;
	}

	

	long newId = ICMLOCAL_Pricer(
						AsOfDate==-1?(ARM_Date*)0:&ARM_Date(XLDateToJulian(AsOfDate)),
						LocalGetNumObjectId (C_Security),
						LocalGetNumObjectId (C_Model),
						l_PricerType,
						l_nbpaths,
						LocalGetNumObjectId (C_Parameters)
						) ;
 
	std::string objName = LocalPersistent::get().getStringId(newId,LOCAL_PRICER_CLASS) ;
	VariantTools::convert(objName,*pRet); 
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

/**

STDMETHODIMP ActiveXModule::ARM_Credit_GetModelFromSummit(BSTR IRcurve,
													  BSTR IDSummit, 
													  BSTR type,	
													  VARIANT *pRet)
{
	// TODO: Add your implementation code here
try
{
	ARM_result C_result;

	_bstr_t b_IRcurve(IRcurve);
	CCString C_IRcurve = b_IRcurve;
	int l_IRcurve = LocalGetNumObjectId( C_IRcurve);

	if (l_IRcurve == ARM_KO)
	{
	ERROR_MSG("Invalid IRcurve",pRet,ARM_ERROR_IRCURVE);
	return S_OK;
	}

	_bstr_t b_IDSummit(IDSummit);
	CCString C_IDSummit = b_IDSummit;

	_bstr_t b_type(type);
	CCString C_type = b_type;

	CCString CurveId;
	CCString CorrCurveId;

	CCString Res;
	CCString curClass = LOCAL_MULTICURVESMODEL_CLASS;

	if(ICMLOCAL_GetModelFromSummit  (l_IRcurve, 
								C_IDSummit,
								C_type,
								CurveId,
								CorrCurveId,
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
**/ 

STDMETHODIMP ActiveXModule::ARM_Credit_SetVolatility(VARIANT *pPricer , 
												 BSTR pVolcurve,
												 VARIANT *pRet)
{
	// TODO: Add your implementation code here

try
{
	ARM_result C_result;

	CCString l_Pricer;
	
	if(VARIANT2CCString (*pPricer, l_Pricer) != S_OK)
		return S_FALSE;

	double Res =0.;

	long Pricerid = LocalGetNumObjectId (l_Pricer);

	if (Pricerid == ARM_KO)
	{
	ERROR_MSG("Invalid Pricer",pRet,ARM_ERROR_PRICER);
	return S_OK;
	}

	_bstr_t b_Volcurve(pVolcurve);
	CCString C_Volcurve = b_Volcurve;

	long Volcurve = LocalGetNumObjectId(C_Volcurve);

	if(ICMLOCAL_SetVolatility(Pricerid,
							  Volcurve,
							  C_result) == ARM_OK)
	{
		Res = C_result.getDouble ();
	}
	else
	{
	ERROR_MSG(C_result.getMsg(),pRet);
	return S_OK;
	}

	_variant_t wrap_Res;
	wrap_Res.Attach(*pRet);
	wrap_Res.dblVal = Res;
	*pRet=wrap_Res.Detach();

	pRet->vt=VT_R8;

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


STDMETHODIMP ActiveXModule::ARM_Credit_GenerateImpliedCurve (BSTR pricerId, 
														 BSTR *pRet)
{
	// TODO: Add your implementation code here

try
{
	ARM_result C_result;

	_bstr_t bpricerId(pricerId);
	CCString l_pricerId = bpricerId;

	CCString retour;

	double Res =0.;

	CCString curClass = LOCAL_ZERO_CURVE_CDS_CLASS;

	long idpricer = LocalGetNumObjectId (l_pricerId);

	if(ICMLOCAL_GenerateImpliedCurve (idpricer,C_result) != ARM_OK)
	{retour = LocalMakeObjectId (C_result.getLong (), curClass);}

	_bstr_t tmpChaine = retour;
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