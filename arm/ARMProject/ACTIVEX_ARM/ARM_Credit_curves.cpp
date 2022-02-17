// ARMModule.cpp : Implementation of CLocal_DLLARMApp and DLL registration.
#include "firsttobeincluded.h"
#include "ActiveXModule.h"


#include "ARM_local_class.h"
#include <ARM/libarm_frometk/VariantTools.h>
#include "CCatl.h"
#include "CCxll.h"


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

#include "ARM_local_interglob.h"
#include "ARMKernel\ccy\currency.h"

#include <CCString.h>

#include <currency.h>
#include <fromto.h>
#include <fstream>

void ERROR_MSG(CCString Err_mess,VARIANT* pRet,long noerr )
{
	_variant_t wrap_Res;
	wrap_Res.Attach(*pRet);
	char NOERR[4];
	sprintf(NOERR,":%d:",noerr);

	CCString TXT = (CCString)"ERROR"+ (CCString) NOERR + Err_mess;

	wrap_Res.SetString(TXT);
	*pRet=wrap_Res.Detach();
	pRet->vt=VT_BSTR;
}

// -----------------------------------------------------------------------------------------------
// Fonctionalités Credit
// -----------------------------------------------------------------------------------------------

STDMETHODIMP ActiveXModule::ARM_Credit_GetInitialCurveFromSummit(VARIANT *pIndex, VARIANT *pCurrency, VARIANT *pCvName, VARIANT *pDate, VARIANT* pTypeValue, VARIANT *pRet)
{
	// TODO: Add your implementation code here
try
{
	_variant_t valDate( pDate );

	ARM_result C_result;

	CCString l_index;
	CCString l_currency;
	CCString l_cvname;
	CCString l_typevalue;
	
	if(VARIANT2CCString (*pIndex, l_index) != S_OK)
		return S_FALSE;

	if(VARIANT2CCString (*pCurrency, l_currency) != S_OK)
		return S_FALSE;

	if(VARIANT2CCString (*pCvName, l_cvname) != S_OK)
		return S_FALSE;

	if(VARIANT2CCString (*pTypeValue, l_typevalue) != S_OK)
		return S_FALSE;

	valDate.ChangeType( VT_R8 );

	double vDate = (double)valDate;

	VECTOR<CCString> matu;
	VECTOR<double> rate;

	long retCode;

	retCode = ARMLOCAL_GetInitialCurveFromSummit (l_index,l_currency,l_cvname,vDate,0,&matu,&rate,C_result);


	long Res;

	if (l_typevalue == "MATU")
		Res = VECTORCCSTRING2VARIANT(matu, pRet);
	else
		Res = VECTORDOUBLE2VARIANT(rate, pRet);

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

STDMETHODIMP ActiveXModule::ARM_Credit_GetDPFromCalypso(double pDate,
										  BSTR pricingEnv, 
										  BSTR issuer,
										  BSTR seniority,
										  BSTR ccy,
										  BSTR forceCurveName,
										  BSTR xmlFileName,
										  BSTR irCurveId,
										  BSTR label,
										  VARIANT *ret)
{
	try
	{		
		std::string sPricingEnv ; VariantTools::convert(pricingEnv,sPricingEnv); 
		std::string sIssuer ; VariantTools::convert(issuer,sIssuer); 
		std::string sSeniority; VariantTools::convert(seniority,sSeniority); 
		std::string sCcy; VariantTools::convert(ccy,sCcy); 
		std::string sForceCurveName; VariantTools::convert(forceCurveName,sForceCurveName); 
		std::string sXmlFileName; VariantTools::convert(xmlFileName,sXmlFileName); 
		std::string sLabel; VariantTools::convert(label,sLabel); 
		ARM_Date AsOf ; VariantTools::convertXLDate(pDate,AsOf); 
		std::string sIrCurveId; VariantTools::convert(irCurveId,sIrCurveId); 
		long ircurve = LocalPersistent::get().getObjectId(sIrCurveId);
		long id= ICMLOCAL_GetDPFromCalypso(AsOf,sIssuer,sSeniority,sCcy,sPricingEnv,sForceCurveName,ircurve,sLabel,sXmlFileName);	
		std::string tem = LocalPersistent::get().getStringId(id,LOCAL_ZERO_CURVE_CDS_CLASS); 
		VariantTools::convert(tem,*ret); 
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


STDMETHODIMP ActiveXModule::ARM_Credit_GetDPFromSummit(double AsOf_,
												BSTR Issuer, 
												BSTR CurveName, 
												BSTR PWCcurveId, 
												BSTR label, 
												VARIANT *pRet)
{
	// TODO: Add your implementation code here
try
{
	ARM_result C_result;

	_bstr_t bIssuer(Issuer);
	CCString C_Issuer = bIssuer;

	_bstr_t bCurveName(CurveName);
	CCString C_CurveName = bCurveName;

	_bstr_t bPWCcurveId(PWCcurveId);
	CCString C_PWCcurveId = bPWCcurveId;

	_bstr_t blabel(label);
	CCString C_label = blabel;

	CCString Res;
	CCString curClass = LOCAL_ZERO_CURVE_CDS_CLASS;

	long ircurve = LocalGetNumObjectId(C_PWCcurveId);

	ARM_Date AsOf; 
	VariantTools::convertXLDate(AsOf_,AsOf); 
	if(ICMLOCAL_GetDPFromSummit(AsOf,
							   C_Issuer,
							   C_CurveName,
							   ircurve,	
							   C_label,
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


STDMETHODIMP ActiveXModule::ARM_Credit_DPMktDataFromSummit(double AsOfDate,
												BSTR Issuer, 
												BSTR CurveName, 
												BSTR Parameter, 
												VARIANT *pRet)
{
	// TODO: Add your implementation code here
try
{
	ARM_result C_result;

	_bstr_t bIssuer(Issuer);
	CCString C_Issuer = bIssuer;

	_bstr_t bCurveName(CurveName);
	CCString C_CurveName = bCurveName;

	_bstr_t bParameter(Parameter);
	CCString C_Parameter = bParameter;

	CCString C_currency;
	CCString C_IndexName;
	CCString C_StandardCDS;
	long CDSAdj;
	
	CCString Res;
	CCString curClass = LOCAL_ZERO_CURVE_CDS_CLASS;

	VECTOR<CCString> Matu;
	VECTOR<double> Spread;
	VECTOR<double> Recovery;
	VECTOR<CCString> CCY;
	VECTOR<CCString> INDEXNAME;
	VECTOR<CCString> STANDARDCDS;

	if(ICMLOCAL_DPMktDataFromSummit(AsOfDate,
								C_Issuer,
								C_CurveName,
								Matu,
								Spread,
								Recovery,
								C_currency,
								C_IndexName,
								CDSAdj, 
							    C_result) == ARM_OK)
	{
		Res = LocalMakeObjectId (C_result.getLong (), curClass);
	}
	else
	{
	ERROR_MSG(C_result.getMsg(),pRet);
	return S_OK;
	}
	
	std::string tmp; 
	ICM_EnumsCnv::toString((qCDS_ADJ) CDSAdj,tmp) ;	// will throw if !ok
	C_StandardCDS=tmp.c_str(); 


	for (int i=0;i<Matu.size();i++)
	{
		CCY.push_back(C_currency);
		INDEXNAME.push_back(C_IndexName);
		STANDARDCDS.push_back(C_StandardCDS); 
	}

	if (C_Parameter == "MATU")
		Res = VECTORCCSTRING2VARIANT(Matu, pRet);
	else if (C_Parameter == "SPREAD")
		Res = VECTORDOUBLE2VARIANT(Spread, pRet);
	else if (C_Parameter == "CCY")
		Res = VECTORCCSTRING2VARIANT(CCY, pRet);
	else if (C_Parameter == "INDEX")
		Res = VECTORCCSTRING2VARIANT(INDEXNAME, pRet);
	else if (C_Parameter == "ADJCDS") 
		Res = VECTORCCSTRING2VARIANT(STANDARDCDS, pRet);
	else
		Res = VECTORDOUBLE2VARIANT(Recovery, pRet);

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

STDMETHODIMP ActiveXModule::ARM_Credit_DPMktDataFromCalypso(double _AsOfDate,
												BSTR _pricingEnv, 
												BSTR _issuer, 
												BSTR _seniority, 
												BSTR _ccy, 
												BSTR _forceCurveName, 
												BSTR _xmlFileName, 
												VARIANT *_parameters,
												VARIANT *pRet)
{
	// TODO: Add your implementation code here
	try
	{
		ARM_Date AsOfDate; VariantTools::convertXLDate(_AsOfDate,AsOfDate); 
		std::string pricingEnv; VariantTools::convert(_pricingEnv,pricingEnv); 
		std::string issuer; VariantTools::convert(_issuer,issuer); 
		std::string seniority; VariantTools::convert(_seniority,seniority); 
		std::string ccy; VariantTools::convert(_ccy,ccy); 
		std::string forceCurveName; VariantTools::convert(_forceCurveName,forceCurveName); 
		std::string xmlFileName; VariantTools::convert(_xmlFileName,xmlFileName); 
		vector<string> vParameters; VariantTools::convert(*_parameters,vParameters); 

		std::vector<std::string> Matus;
		std::vector<double> Spreads;
		double Recovery;
		std::string outputCcy;
		std::string indexName;
		qCDS_ADJ adj ;

		ICMLOCAL_DPMktDataFromCalypso(AsOfDate,
									pricingEnv,
									issuer,
									seniority,
									ccy,
									forceCurveName,
									xmlFileName,
										//output: 
									Matus,
									Spreads,
									Recovery,
									outputCcy, 
									indexName,
									adj) ;

		std::string tmp; 

		vector<vector<string> > vectorResult;
		vector<string> stringSpread(Spreads.size());
		for ( int i= 0; i< vParameters.size(); i++) {
			ICM_EnumsCnv::toString(adj,tmp) ;	// will throw if !ok
			if  (vParameters[i] =="MATU") {		
					vectorResult.push_back(Matus);
					continue;
			}
			if  (vParameters[i] =="SPREAD"){
					char TMP[20];
					
					for ( int i=0; i<Spreads.size() ; i++) {
						sprintf(TMP,"%8.4lf",Spreads[i]);
						stringSpread[i] = string(TMP);
					}
					vectorResult.push_back(stringSpread);
					continue;
			}	
			if  (vParameters[i] =="CCY") {		
					std::vector<std::string> ccy(Matus.size(),outputCcy); 
					vectorResult.push_back(ccy);
					continue;
			}
			if  (vParameters[i] =="INDEX") {
					std::vector<std::string> index(Matus.size(),indexName); 
					vectorResult.push_back(index);
					continue;
			}
			if  (vParameters[i] =="ADJCDS") {
					std::vector<std::string> adj(Matus.size(),ICM_EnumsCnv::toString(adj)); 
					vectorResult.push_back(adj);
					continue;
			}
			if  (vParameters[i] =="RECOVERY") {
					char TMP[20];
					stringSpread.clear();
					stringSpread.resize(Matus.size());
					for ( int i=0; i<Spreads.size() ; i++) {
						sprintf(TMP,"%8.4lf",Recovery);
						stringSpread[i] = string(TMP);
					}
					vectorResult.push_back(stringSpread);
					continue;
			}else {
				ICMTHROW(ERR_INVALID_ARGUMENT,"ARM_Credit_DPMktDataFromCalypso: parameters are MATU|SPREAD|CCY|INDEX|ADJCDS|MATU|RECOVERY"); 
			}
		}
		VariantTools::convert(vectorResult,*pRet);
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



STDMETHODIMP ActiveXModule::ARM_Credit_ConstantDefaultCurve(DATE AsOf_,
														VARIANT pTenors, 
														VARIANT pRates, 
														double  Recovery, 
														BSTR	IRCurveId,
														BSTR	bCurrency, 
														BSTR	bLabel, 
														BSTR	bAdjCalType, 
														BSTR	bIsSummit, 
														BSTR	calibrationData_,
														int lag,
														BSTR	calibrationAlgo_,
														BSTR	params_,
														BSTR	intRule_,
														BSTR	adjStartRule_,
														VARIANT *pRet)
{
	// TODO: Add your implementation code here
try
{
	ARM_result C_result;

	
	// VECTOR<double> V_Rates;

	// CCString l_IRcurveId;
	// CCString l_currency;
	// CCString l_Label;

	// int lag = 2;
	long size =0;

	std::string C_Currency ; VariantTools::convert(bCurrency,C_Currency); 
	// _bstr_t b_Currency(bCurrency);
	// CCString C_Currency = b_Currency;
	int lag = ARM_Currency(C_Currency.c_str()).GetCreditStartDateLag();

	std::string C_Label  ; VariantTools::convert(bLabel,C_Label ); 
	// _bstr_t b_Label(bLabel);
	// CCString C_Label = b_Label;

	std::string C_IRCurveId  ; VariantTools::convert(IRCurveId,C_IRCurveId ); 
	// _bstr_t b_IRCurveId(IRCurveId);
	// CCString C_IRCurveId = b_IRCurveId;

	std::string C_IsSummit   ; VariantTools::convert(bIsSummit,C_IsSummit  ); 
	// _bstr_t b_IsSummit(bIsSummit);
	// CCString C_IsSummit = b_IsSummit;

	// _bstr_t b_AdjCalType(bAdjCalType);
	// CCString C_AdjCalType = b_AdjCalType;
	// qCDS_ADJ AdjCalType= qCredit_Adjust20;
	
	ARM_Date AsOf ; VariantTools::convertXLDate(AsOf_,AsOf); 
	std::string calibrationData; VariantTools::convert(calibrationData_,"STD",calibrationData); 
	qDEFCURVE_CALIB_ALGO calibrationAlgo ; VariantTools::convert(calibrationAlgo_,"DICHO",calibrationAlgo); 
	std::string sParamId; VariantTools::convert(params_,"",sParamId); 

	qCDS_ADJ AdjCalType ; VariantTools::convert(bAdjCalType,"STDCDS",AdjCalType); 
	// if((AdjCalType = ARM_ConvAdjCalCDS (C_AdjCalType, C_result)) == ARM_DEFAULT_ERR)
	// {
	// 	ERROR_MSG("Invalid ADJ",pRet,ARM_ERROR_FREQ);
	// 	return S_OK;
	// }

	bool issummit = false;
	if (C_IsSummit == "Y") issummit = true; 

	
	
	std::vector<std::string> V_Tenors; 
	VariantTools::convert(pTenors,V_Tenors); 

	std::vector<double> V_Rates ;
	VariantTools::convert(pRates,V_Rates); 

	std::string sIntRule ; VariantTools::convert(intRule_,"ADJ",sIntRule); 
	std::string sAdjStartRule; VariantTools::convert(adjStartRule_,"ADJ",sAdjStartRule); 
	long intRule = ARM_ConvIntRule(sIntRule.c_str()); 
	long adjStartRule= ARM_ConvStartAdjRule(sAdjStartRule); 
	
	// if(VARIANT2VECTORDOUBLE (*pRates, V_Rates,size) != S_OK)
	// 	return S_FALSE;

	// CCString Res;
	// CCString curClass = LOCAL_ZERO_CURVE_CDS_CLASS;

	// long ircurve = LocalGetNumObjectId(C_IRCurveId);

	// if (ircurve == ARM_KO)
	// {
	// 	ERROR_MSG("Invalid IrCurve",pRet,ARM_ERROR_IRCURVE);
	// 	return S_OK;
	// }


	long objId = ICMLOCAL_ConstantDefaultCurve (AsOf, 
									  V_Tenors, 
									  V_Rates,
									  Recovery,
									  LocalPersistent::get().getObjectId(C_IRCurveId),
									  intRule,
									  adjStartRule,
									  C_Currency, 
									  C_Label, 
									  AdjCalType,
									  issummit,
									  calibrationAlgo,
									  -1, //volcurve
										calibrationData,
										lag,
										LocalPersistent::get().getObjectId(sParamId),
										-1
										) ; 

	std::string objName = LocalPersistent::get().getStringId(objId,LOCAL_ZERO_CURVE_CDS_CLASS); 
	VariantTools::convert(objName,*pRet); 
	return S_OK;
}
	catch(Exception&e)
	{
		return createErrorInfo("ARM_Credit_ConstantDefaultCurve",e.GetErrorString()); 
	}
	catch(std::exception&e)
	{
		return createErrorInfo("ARM_Credit_ConstantDefaultCurve",e.what()); 
	}
	catch(...) 
	{ 
		return createErrorInfo("ARM_Credit_ConstantDefaultCurve","Unknown exception"); 
	}
	return E_FAIL ;
}


STDMETHODIMP ActiveXModule::ARM_Credit_ZeroCouponDefaultCurveFromSummit(double AsOfDate,
																	BSTR   bIssuer, 
																	BSTR   bCurrency, 
																	BSTR   bCvName,	
															        BSTR   IRCurveId,
																	BSTR   bLabel,
																	VARIANT *pRet)
{
	// TODO: Add your implementation code here
try
{
	ARM_result C_result;

	long size =0;

	_bstr_t b_Currency(bCurrency);
	CCString C_Currency = b_Currency;

	_bstr_t b_Issuer(bIssuer);
	CCString C_Issuer = b_Issuer;

	_bstr_t b_CvName(bCvName);
	CCString C_CvName = b_CvName;

	_bstr_t b_IRCurveId(IRCurveId);
	CCString C_IRCurveId = b_IRCurveId;

	_bstr_t b_Label(bLabel);
	CCString C_Label = b_Label;

	CCString Res;
	CCString curClass = LOCAL_ZERO_CURVE_CDS_CLASS;

	long ircurve = LocalGetNumObjectId(C_IRCurveId);

	if (ircurve == ARM_KO)
	{
		ERROR_MSG("Invalid IrCurve",pRet,ARM_ERROR_IRCURVE);
		return S_OK;
	}


	if(ICMLOCAL_GetZC_DP_FromSummit (C_Issuer,
									 C_Currency,
									 C_CvName,
									 AsOfDate,
									 ircurve,
									 C_Label,
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



STDMETHODIMP ActiveXModule::ARM_Credit_ZeroCouponDefaultCurveFromCalypso(double AsOfDate_,
																	  BSTR pricingEnv_, 
																	  BSTR issuer_,
																	  BSTR seniority_,
																	  BSTR ccy_,
																	  BSTR forceCurveName_,
																	  BSTR xmlFileName_,
																	  BSTR irCurveId_,
																	  BSTR label_,
																	  VARIANT *pRet)
{
	// TODO: Add your implementation code here
	try
	{
		std::string pricingEnv ; VariantTools::convert(pricingEnv_,pricingEnv); 
		std::string issuer ; VariantTools::convert(issuer_,issuer); 
		std::string seniority; VariantTools::convert(seniority_,seniority); 
		std::string ccy; VariantTools::convert(ccy_,ccy); 
		std::string forceCurveName; VariantTools::convert(forceCurveName_,forceCurveName); 
		std::string xmlFileName; VariantTools::convert(xmlFileName_,xmlFileName); 
		std::string label; VariantTools::convert(label_,label); 
		ARM_Date AsOf ; VariantTools::convertXLDate(AsOfDate_,AsOf); 
		std::string strCurveId; VariantTools::convert(irCurveId_,strCurveId); 

		long ircurve = LocalPersistent::get().getObjectId(strCurveId); 
		long id = ICMLOCAL_GetZC_DP_FromCalypso(AsOf,pricingEnv,issuer,seniority,ccy,forceCurveName,xmlFileName,ircurve,label); 
		std::string item = LocalPersistent::get().getStringId(id,LOCAL_ZERO_CURVE_CDS_CLASS); 
		VariantTools::convert(item,*pRet); 
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




STDMETHODIMP ActiveXModule::ARM_Credit_InputDefaultCurve(VARIANT  AsOfDate,
														VARIANT  pDates, 
														VARIANT  pRates, 
														double  Recovery, 
														BSTR	IRCurveId,
														BSTR	bCurrency, 
														BSTR	bLabel, 
														// BSTR	bInterpolType, 
														VARIANT *pRet)
{
	// TODO: Add your implementation code here
	try
	{
		ARM_Date AsOf; VariantTools::convert(AsOfDate,AsOf); 
		std::vector<ARM_Date> Dates; VariantTools::convert(pDates,Dates); 
		ARM_Vector Rates; VariantTools::convert(pRates,Rates); 
		std::string irCurveId; VariantTools::convert(IRCurveId,irCurveId); 
		long ircurveid = LocalPersistent::get().getObjectId(irCurveId); 
		std::string ccy ; VariantTools::convert(bCurrency,ccy); 
		std::string label ; VariantTools::convert(bLabel,label); 

		long objId = ICMLOCAL_InputDefCurve_Dates(AsOf,
				Dates,
				Rates,
				Recovery,
				ircurveid,
				ccy,
				label); 

		std::string objName = LocalPersistent::get().getStringId(objId,LOCAL_ZERO_CURVE_CDS_CLASS); 
		VariantTools::convert(objName,*pRet); 
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

		// ARM_result C_result;

		// VECTOR<double> V_Dates;
		// VECTOR<double> V_Rates;

		// CCString l_IRcurveId;
		// CCString l_currency;
		// CCString l_Label;

		// long size =0;

		// _bstr_t b_Currency(bCurrency);
		// CCString C_Currency = b_Currency;

		// _bstr_t b_Label(bLabel);
		// CCString C_Label = b_Label;

		// _bstr_t b_IRCurveId(IRCurveId);
		// CCString C_IRCurveId = b_IRCurveId;

		// std::string interpolStr ; VariantTools::convert(bInterpolType,interpolStr);
		
		// qINTERPOL_TYPE InterpolType ; 
		// ICM_EnumsCnv::cnv(interpolStr,InterpolType ); 

		// if(VARIANT2VECTORDOUBLE(*pDates, V_Dates,size) != S_OK)
		// 	return S_FALSE;

		// if(VARIANT2VECTORDOUBLE (*pRates, V_Rates,size) != S_OK)
		// 	return S_FALSE;
/** 
		CCString Res;
		CCString curClass = LOCAL_ZERO_CURVE_CDS_CLASS;

		long ircurve = LocalGetNumObjectId(C_IRCurveId);

		if (ircurve == ARM_KO)
		{
			ERROR_MSG("Invalid IrCurve",pRet,ARM_ERROR_IRCURVE);
			return S_OK;
		}


		if(ICMLOCAL_InputDefCurve_Dates (AsOfDate, 
										  V_Dates, 
										  V_Rates,
										  Recovery,
										  ircurve,
										  C_Currency, 
										  C_Label, 
										  // InterpolType,
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
	**/ 
}


STDMETHODIMP ActiveXModule::ARM_Credit_FixingCurve(VARIANT *pDates, 
														VARIANT *pValues, 
														double  AsOfDate,
														BSTR	B_IndexName,
														BSTR	B_IndexID, 
														VARIANT *pRet)
{
	// TODO: Add your implementation code here
	try
	{
		ARM_result C_result;

		vector<ARM_Date> V_Dates;
		vector<double> V_Values;
		string Name =""; 
		VariantTools::convert(B_IndexName,"", Name);
		string NameID ="";
		VariantTools::convert(B_IndexID,"", NameID);
		ARM_Date ADate; ADate.Today();
		if (AsOfDate != -1) 
			VariantTools::convertXLDate(AsOfDate, ADate);

		long size =0;

		VariantTools::convert(*pDates,V_Dates);
		VariantTools::convert(*pValues, V_Values);
	
		long newId = ICMLOCAL_Fixing_Curve (V_Dates, V_Values,
										  ADate,
										  Name,
										  (long)LocalGetNumObjectId(CCString(NameID.c_str())),
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

STDMETHODIMP ActiveXModule::ARM_Credit_SetVolCurve(BSTR Model , 
												 BSTR pVolcurve,
												 VARIANT *pRet)
{
	// TODO: Add your implementation code here

try
{
	ARM_result C_result;

	string l_Model;
	string l_VolId;
	VariantTools::convert(Model,l_Model);
	VariantTools::convert(pVolcurve,l_VolId);

	double Res =0.;

	long ModelId = LocalGetNumObjectId (CCString(l_Model.c_str()));

	if (ModelId == ARM_KO)
	{
		ERROR_MSG("Invalid ModelId",pRet,ARM_ERROR_PRICER);
		return S_OK;
	}
	long Volcurve = LocalGetNumObjectId(CCString(l_VolId.c_str()));

	if(ICMLOCAL_SetVolCurve(ModelId,
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

STDMETHODIMP ActiveXModule::ARM_Credit_DefCurveIntensityPWC( double AsOfDate,
														VARIANT *pMatuRates, 
														VARIANT *pInputs, 
														double Type,
														double  Recovery, 
														BSTR	IRCurveId,
														BSTR	bCurrency, 
														BSTR	bLabel, 
														BSTR	VolCurveId, 
														BSTR	calibrationAlgo,
														int lag,
														BSTR	params_,
														VARIANT *pRet)
{
	try
	{
		ARM_result C_result;

		ARM_Date	C_AsOfDate ;
		VECTOR<double> C_matu;
		VECTOR<double> C_Inputs;
	
		std::string C_ircurve;
		std::string C_Currency;
		std::string C_label ;
		std::string C_VolCurveId;
		std::string C_CalibrationAlgo ;
		std::string sParamId;
		
		VariantTools::convertXLDate(AsOfDate, C_AsOfDate);
		long size = 0;


		if(VARIANT2VECTORDOUBLE (*pMatuRates, C_matu,size) != S_OK)
			return S_FALSE;

			
		if(VARIANT2VECTORDOUBLE (*pInputs, C_Inputs,size) != S_OK)
			return S_FALSE;

		VariantTools::convert(IRCurveId,C_ircurve);
		VariantTools::convert(bCurrency,"DEFAULT",C_Currency);
		VariantTools::convert(bLabel,"NONE",C_label);
		VariantTools::convert(VolCurveId,C_VolCurveId);
		VariantTools::convert(calibrationAlgo,"STD",C_CalibrationAlgo);
		VariantTools::convert(params_,"",sParamId);

		if ( C_Currency=="DEFAULT")
		{
			ARM_result currencyres;
			ARMLOCAL_ARM_GetDefaultCurrency (currencyres);
			C_Currency = currencyres.getString();

		}
		int lag = 2;
		lag = ARM_Currency(C_Currency.c_str()).GetCreditStartDateLag();
		
		long objId = ICMLOCAL_CstDefCurve_Dates(C_AsOfDate,
											C_matu,
											C_Inputs,
											Type,
											Recovery,
											(long)LocalGetNumObjectId(CCString(C_ircurve.c_str())),
											CCString(C_Currency.c_str()),
											CCString(C_label.c_str()),
											(long)LocalGetNumObjectId(CCString(C_VolCurveId.c_str())),
											C_CalibrationAlgo,
											lag,
											LocalPersistent::get().getObjectId(sParamId),
											-1
											);

		std::string objName = LocalPersistent::get().getStringId(objId,LOCAL_ZERO_CURVE_CDS_CLASS); 
		VariantTools::convert(objName,*pRet); 

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

STDMETHODIMP ActiveXModule::ARM_Credit_DefCurvePWC_ABS(double AsOfDate,
														VARIANT *Tenors, 
														VARIANT *Rates, 
														VARIANT	*UpFront,
														VARIANT	*RefValIds,
														double  Recovery, 
														BSTR	IRCurveId,
														BSTR	bCurrency, 
														BSTR	bLabel, 
														BSTR	adjCalType,
														BSTR	IsSummitCurve,
														BSTR	VolCurveId, 
														BSTR	Accrued, 
														BSTR	calibrationAlgo,
														BSTR	calibrationData,
														int lag,
														VARIANT *pRet)
{

	try
	{
		ARM_result C_result;

		ARM_Date		C_AsOfDate ;
		VECTOR<std::string>	C_matu;
		VECTOR<double>	C_Rates;
		VECTOR<double>	C_UpFronts;

		VECTOR<std::string>	C_RefIds;
		VECTOR<long>	l_RefIds;

	
		std::string C_ircurve;
		std::string C_Currency;
		std::string C_label ;
		std::string C_VolCurveId;

		std::string C_AdjCal ;

		qDEFCURVE_CALIB_ALGO C_CalibrationAlgo ;
		std::string C_CalibrationData;
		std::string C_IsSummit;

		bool C_Cal_Summit=false ;
		qCDS_ADJ l_AdjCal = qCredit_Adjust20; // STDCDS
		
		VariantTools::convertXLDate(AsOfDate, C_AsOfDate);
		long size = 0;


		VariantTools::convert(*Tenors,C_matu);

			
		if(VARIANT2VECTORDOUBLE (*Rates, C_Rates,size) != S_OK)
		{
			ERROR_MSG("Invalid Rates",pRet,ARM_ERROR_FREQ);
			return S_FALSE;
		}
		if ( C_Rates.size() != C_matu.size())
		{
			ERROR_MSG(" check your maturities & spreads array",pRet,ARM_ERROR_FREQ);
			return S_FALSE;

		}



		VariantTools::convert(IRCurveId,C_ircurve);
		VariantTools::convert(bCurrency,"DEFAULT",C_Currency);
		VariantTools::convert(bLabel,"NONE",C_label);
		VariantTools::convert(VolCurveId,C_VolCurveId);
		VariantTools::convert(calibrationAlgo,"DICHO",C_CalibrationAlgo);
		VariantTools::convert(calibrationData,"STD",C_CalibrationData);
		VariantTools::convert(IsSummitCurve,"Y",C_IsSummit);
		VariantTools::convert(adjCalType,"STDCDS",C_AdjCal);

		if (C_IsSummit=="Y") C_Cal_Summit = true ;

		VariantTools::convert(*RefValIds,C_RefIds);
		VariantTools::convert(*UpFront,C_UpFronts);

		l_RefIds.resize(C_RefIds.size());

		for (int i=0;i<C_RefIds.size();i++)
			{l_RefIds[i]=(long)LocalGetNumObjectId(C_RefIds[i].c_str());}


		if ( C_Currency=="DEFAULT")
		{
			ARM_result currencyres;
			ARMLOCAL_ARM_GetDefaultCurrency (currencyres);
			C_Currency = currencyres.getString();

		}

		if((l_AdjCal = ARM_ConvAdjCalCDS (CCString(C_AdjCal.c_str()), C_result)) == ARM_DEFAULT_ERR)
		{
			ERROR_MSG("Invalid ADJ",pRet,ARM_ERROR_FREQ);
			return S_FALSE;
		}

		int lag = 2;
		lag = ARM_Currency(C_Currency.c_str()).GetCreditStartDateLag();

		
		long objId = ICMLOCAL_ABS_PWCDefaultCurve(C_AsOfDate, 
												C_matu, 
												C_Rates, 
												Recovery,
												(long)LocalGetNumObjectId(C_ircurve.c_str()),
												C_Currency,
												C_label,
												l_AdjCal,
												C_Cal_Summit,
												// C_Cal_Brent_Solver,
												C_CalibrationAlgo,
												(long)LocalGetNumObjectId(C_VolCurveId.c_str()),
												C_CalibrationData,
												(int)lag,
												C_UpFronts,
												l_RefIds,
												-1);

		std::string objName = LocalPersistent::get().getStringId(objId,LOCAL_ZERO_CURVE_CDS_CLASS); 
		VariantTools::convert(objName,*pRet); 

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