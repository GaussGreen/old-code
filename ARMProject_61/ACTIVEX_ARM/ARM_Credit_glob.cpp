// ARMModule.cpp : Implementation of CLocal_DLLARMApp and DLL registration.


#include "firsttobeincluded.h"
#include "ICMKernel\glob\icm_mktdatamng.h"
#include "ActiveXModule.h"

#include "ARM_local_class.h"

#include <ARM/libarm_frometk/VariantTools.h>
#include <ARM/libarm_frometk/arm_local_xgiga.h>
#include <ARM\libarm_frometk\arm_local_paesexml_calypso.h> 
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
#include "ICM_local_cf.h"
#include "ICM_local_leg.h"

#include "ARM_local_interglob.h"
#include "ARM_local_assetswap.h"

#include <fromto.h>

#define EDITOR							"notepad "
#define XML_EDITOR							"iexplore "
#define VIEW_FILE_PREFIX				"VF"
#define VIEW_FILE_CLIENT_LOCATION		"c:\\temp\\"



STDMETHODIMP ActiveXModule::ARM_FxConvert(VARIANT *pccy1, 
									  VARIANT *pccy2, 
									  VARIANT *pAsOfDate,
									  VARIANT *pCvName,
									  VARIANT *pRet)
{
	// TODO: Add your implementation code here
try
{
	_variant_t valDate( pAsOfDate );

	ARM_result C_result;

	CCString l_ccy1;
	CCString l_ccy2;
	CCString l_cvName;
	
	if(VARIANT2CCString (*pccy1, l_ccy1) != S_OK)
		return S_FALSE;

	if(VARIANT2CCString (*pccy2, l_ccy2) != S_OK)
		return S_FALSE;

	if(VARIANT2CCString (*pCvName, l_cvName) != S_OK)
		return S_FALSE;

	valDate.ChangeType( VT_R8 );

	double vDate = (double)valDate;

	double Res = 0.;

	if(ARMLOCAL_FxConvert(l_ccy1,l_ccy2,vDate,1.,l_cvName,C_result) == ARM_OK)
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
catch(...)
{
	return E_FAIL;
}
}


STDMETHODIMP ActiveXModule::ARM_FxConvertFromCalypso(BSTR ccy1_, 
												 BSTR ccy2_,
												 BSTR cvName_,
												 DATE asof_,
												 VARIANT*pRet)
{
	// TODO: Add your implementation code here
	try
	{
		std::string ccy1 ; VariantTools::convert(ccy1_,ccy1); 
		std::string ccy2 ; VariantTools::convert(ccy2_,ccy2); 
		std::string cvName ; VariantTools::convert(cvName_,cvName); 
		if (cvName=="") cvName="MO"; 
		ARM_Date asof; VariantTools::convertXLDate(asof_,asof); 
		double ret; 
		ARM_CalypsoToolkit::GetFXRate(ccy1,ccy2,cvName,asof,ret); 
		VariantTools::convert(ret,*pRet) ;
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


STDMETHODIMP ActiveXModule::ARM_DiscountPrice(VARIANT *pCurve, VARIANT *pMatu, VARIANT *pRet)
{
	// TODO: Add your implementation code here

try
{
	ARM_result C_result;

	_variant_t Matu(pMatu);

	CCString l_Curve;
	
	Matu.ChangeType( VT_R8 );
	double vDate = (double)Matu;

	if(VARIANT2CCString (*pCurve, l_Curve) != S_OK)
		return S_FALSE;

	double Res =0.;

	long curve = LocalGetNumObjectId (l_Curve);

	if (curve == ARM_KO)
	{
	ERROR_MSG("Invalid IrCurveId",pRet,ARM_ERROR_IRCURVE);
	return S_OK;
	}


	if(ARMLOCAL_DiscountPrice (curve, vDate, C_result) == ARM_OK)
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
catch(...)
{
	return E_FAIL;
}
}


// -----------------------------------------------------------------------------------------------
// Fonctionalités Credit
// -----------------------------------------------------------------------------------------------



STDMETHODIMP ActiveXModule::ARM_Credit_Delivery(VARIANT *pAsOfDate, 
											VARIANT *pTenorContract, 
											VARIANT *pRet)
{
	// TODO: Add your implementation code here
try
{
	_variant_t valDate(pAsOfDate);

	ARM_result C_result;

	CCString l_TenorContract;

	if(VARIANT2CCString (*pTenorContract, l_TenorContract) != S_OK)
		return S_FALSE;

	valDate.ChangeType( VT_R8 );

	double vDate = (double)valDate;

	double Res = 0.;

	if(ICMLOCAL_Delivery(vDate,l_TenorContract,C_result) == ARM_OK)
	{
		Res = Local_ARMDATE2XLDATE(C_result.getString());
	}
	else
	{
	ERROR_MSG(C_result.getMsg(),pRet);
	return S_OK;
	}

	_variant_t wrap_Res;
	wrap_Res.Attach(*pRet);
	wrap_Res.dblVal =Res;
	*pRet=wrap_Res.Detach();

	pRet->vt=VT_R8;

	return S_OK;
}
catch(...)
{
	return E_FAIL;
}
}


STDMETHODIMP ActiveXModule::ARM_Credit_CptInterpolDefCurve(BSTR pCurve,VARIANT pTenor,VARIANT *pRet)
{
	// TODO: Add your implementation code here

try
{
	std::string defCurveId; VariantTools::convert(pCurve,defCurveId); 
	ARM_Date date; 
	std::string tenor ;
	if (VariantTools::isXLDate(pTenor)) VariantTools::convert(pTenor,date); 
	else VariantTools::convert(pTenor,tenor); 

	double res; 
	if (!tenor.empty()) 
		res=ICMLOCAL_CptImplicitSpreadInterpol(
		LocalPersistent::get().getObjectId(defCurveId),tenor); 
	else 
		res=ICMLOCAL_CptImplicitSpreadInterpol(
		LocalPersistent::get().getObjectId(defCurveId),date); 
	VariantTools::convert(res,*pRet); 
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
 
STDMETHODIMP ActiveXModule::ARM_Credit_createFlatCurve(BSTR pCurve,VARIANT pTenor,VARIANT *pRet)
{
	// TODO: Add your implementation code here

	try
	{
		std::string defCurveId; VariantTools::convert(pCurve,defCurveId); 
		ARM_Date date; 
		std::string tenor ;
		if (VariantTools::isXLDate(pTenor)) VariantTools::convert(pTenor,date); 
		else VariantTools::convert(pTenor,tenor); 

		long objId; 
		if (!tenor.empty()) 
			objId=ICMLOCAL_createFlatCurve(
			LocalPersistent::get().getObjectId(defCurveId),tenor,-1); 
		else 
			objId=ICMLOCAL_createFlatCurve(
			LocalPersistent::get().getObjectId(defCurveId),date,-1); 
		
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

STDMETHODIMP ActiveXModule::ARM_Credit_createDefCurveFromBase(BSTR pCurveCDS,BSTR pCurveIndex,VARIANT vBase,VARIANT *pRet)
{
	// TODO: Add your implementation code here

	try
	{
		std::string defCurveCDS; VariantTools::convert(pCurveCDS,defCurveCDS); 
		std::string defCurveIndex; VariantTools::convert(pCurveIndex,defCurveIndex); 
		ARM_Date date; 
		ARM_Vector base ;
		VariantTools::convert(vBase,base); 
		
		long objId; 
		objId=ICMLOCAL_createDefCurveFromBase(LocalPersistent::get().getObjectId(defCurveCDS),
											  LocalPersistent::get().getObjectId(defCurveIndex),base,-1);
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

STDMETHODIMP ActiveXModule::ARM_Credit_DefaultProba(VARIANT *pCurve, 
												VARIANT *pMatu, 
												VARIANT *pRet)
{
	// TODO: Add your implementation code here

try
{
	ARM_result C_result;

	_variant_t Matu(pMatu);

	CCString l_Curve;
	
	Matu.ChangeType( VT_R8 );

	double vDate = (double)Matu;

	if(VARIANT2CCString (*pCurve, l_Curve) != S_OK)
		return S_FALSE;

	double Res =0.;

	long defcurve = LocalGetNumObjectId (l_Curve);

	if (defcurve == ARM_KO)
	{
	ERROR_MSG("Invalid defCurve",pRet,ARM_ERROR_DFCURVE);
	return S_OK;
	}

	if(ICMLOCAL_DefaultProba(defcurve, vDate, C_result) == ARM_OK)
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
catch(...)
{
	return E_FAIL;
}
}


STDMETHODIMP ActiveXModule::ARM_Credit_GetBeta (VARIANT *pPricer, 
											VARIANT *pLabel,
											VARIANT *pRet)
{
	// TODO: Add your implementation code here

try
{
	ARM_result C_result;

	CCString l_Pricer;
	CCString l_Label;
	
	if(VARIANT2CCString (*pPricer, l_Pricer) != S_OK)
		return S_FALSE;

	if(VARIANT2CCString (*pLabel, l_Label) != S_OK)
		return S_FALSE;

	double Res =0.;

	long pricerid = LocalGetNumObjectId (l_Pricer);

	if (pricerid == ARM_KO)
	{
	ERROR_MSG("Invalid Pricer",pRet,ARM_ERROR_PRICER);
	return S_OK;
	}


	if(ICMLOCAL_GetBeta(pricerid, l_Label, C_result) == ARM_OK)
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

STDMETHODIMP ActiveXModule::ARM_Credit_Price (VARIANT *pPricer , 
										  VARIANT *pAsOfDate,
										  VARIANT *pRet)
{
	// TODO: Add your implementation code here

try
{
	ARM_result C_result;

	_variant_t v_AsOfDate(pAsOfDate);
	v_AsOfDate.ChangeType( VT_R8 );
	double C_AsOfDate = ((double) v_AsOfDate);

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

	if(ICMLOCAL_Price(Pricerid,
					C_AsOfDate,
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
catch(...)
{
	return E_FAIL;
}
}


STDMETHODIMP ActiveXModule::ARM_Credit_Spread (VARIANT *pPricer , 
										  VARIANT *pMTM,
										  VARIANT *pRet)
{
	// TODO: Add your implementation code here

try
{
	ARM_result C_result;

	_variant_t v_MTM(pMTM);
	v_MTM.ChangeType( VT_R8 );
	double C_MTM = ((double) v_MTM);

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

	if(ICMLOCAL_Spread(Pricerid,
					C_MTM,
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
catch(...)
{
	return E_FAIL;
}
}



STDMETHODIMP ActiveXModule::ARM_Credit_RiskyDuration (BSTR pDefCurve, 
													  VARIANT vDate, 
													  VARIANT *pRet)
{
	// TODO: Add your implementation code here

try
{
	std::string defCurve ; VariantTools::convert(pDefCurve,defCurve) ; 
	long defCurveId = LocalPersistent::get().getObjectId(defCurve); 
	double res=0; 
	if( VariantTools::isXLDate(vDate) )
	{
		ARM_Date date ; VariantTools::convert(vDate,date); 
		res = ICMLOCAL_RiskyDuration(defCurveId,date); 
	}
	else
	{
		std::string tenor; VariantTools::convert(vDate,tenor); 
		res = ICMLOCAL_RiskyDuration(defCurveId,tenor); 
	}
	VariantTools::convert(res,*pRet); 
	return S_OK ;
}
/**	ARM_result C_result;

	CCString l_DefCurve;
	
	_bstr_t bTenor(Tenor);
	CCString C_Tenor = bTenor;

	if(VARIANT2CCString (*pDefCurve, l_DefCurve) != S_OK)
		return S_FALSE;

	double Res =0.;

	long defcurveid = LocalGetNumObjectId (l_DefCurve);

	if (defcurveid == ARM_KO)
	{
	ERROR_MSG("Invalid DefCurve",pRet,ARM_ERROR_DFCURVE);
	return S_OK;
	}


	if(ICMLOCAL_RiskyDuration(defcurveid,Date,C_Tenor,C_result) == ARM_OK)
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
}**/ 
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


STDMETHODIMP ActiveXModule::ARM_Credit_GetDefProbTranche (BSTR pPricer, 
													  double yearterm,
													  VARIANT *pRet)
{
	// TODO: Add your implementation code here

try
{
	ARM_result C_result;

	_bstr_t bPricer(pPricer);
	CCString C_Pricer = bPricer;


	double Res =0.;

	long Pricerid = LocalGetNumObjectId (C_Pricer);

	if (Pricerid == ARM_KO)
	{
	ERROR_MSG("Invalid Pricer",pRet,ARM_ERROR_DFCURVE);
	return S_OK;
	}

	if(ICMLOCAL_GetDefProbTranche(Pricerid,yearterm,C_result) == ARM_OK)
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


STDMETHODIMP ActiveXModule::ARM_Credit_GetCleanSpreadTranche (VARIANT *pPricer, 
											VARIANT *pPlot,
											VARIANT *pRet)
{
	// TODO: Add your implementation code here

try
{
	ARM_result C_result;

	CCString l_Pricer;
	CCString l_Plot;
	
	if(VARIANT2CCString (*pPricer,l_Pricer) != S_OK)
		return S_FALSE;

	if(VARIANT2CCString (*pPlot,l_Plot) != S_OK)
		return S_FALSE;

	double Res =0.;

	long Pricerid = LocalGetNumObjectId (l_Pricer);

	if (Pricerid == ARM_KO)
	{
	ERROR_MSG("Invalid Pricer",pRet,ARM_ERROR_DFCURVE);
	return S_OK;
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


STDMETHODIMP ActiveXModule::ARM_Credit_CDONPV (VARIANT *pPricer , 
										  VARIANT *pCPTTYPE,
										  VARIANT *pRet)
{
	// TODO: Add your implementation code here

try
{
	if (!pPricer) 
		ICMTHROW(ERR_INVALID_ARGUMENT,"Pricer not provided"); 
	ARM_result C_result;

	// CCString l_Pricer;
	// if(VARIANT2CCString (*pPricer, l_Pricer) != S_OK)
	// 	return S_FALSE;
	std::string pricerId; 
	VariantTools::convert(*pPricer,pricerId); 
	
	qCMPMETH measure ;
	VariantTools::convert(pCPTTYPE,"",measure); // will throw if measure not ok. 
 
	
	// double Res =0.;

	long Pricerid = LocalPersistent::get().getObjectId(pricerId);

// 	if (Pricerid == ARM_KO)
// 	{
// 	ERROR_MSG("Invalid Pricer",pRet,ARM_ERROR_PRICER);
// 	return S_OK;
// 	}
	
	double Res = ICMLOCAL_NPV(Pricerid,measure); 
	VariantTools::convert(Res,*pRet); 
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


STDMETHODIMP ActiveXModule::ARM_Credit_CorrMatrix(VARIANT *pLabels,
											  VARIANT *pCoefs,
											  double AsOf_,
											  BSTR	 Name,
											  VARIANT *pRet)
{
	// TODO: Add your implementation code here
try
{
	ARM_result C_result;

	std::vector<std::string> V_Labels;
	VECTOR<double> V_Coefs;
	long size ;
	

	std::string C_Name; VariantTools::convert(Name,C_Name); 
	ARM_Date AsOf; VariantTools::convertXLDate(AsOf_,AsOf); 

	// if(VARIANT2VECTORCCSTRING(*pLabels,V_Labels,il) != S_OK)
	// 	return S_FALSE;
	VariantTools::convert(*pLabels,V_Labels); 

	size = V_Labels.size();

	if (!size)
	{	
		ERROR_MSG("Empty labels",pRet);
		return S_OK;
	}

	long il; 
	if(VARIANT2VECTORDOUBLE (*pCoefs,V_Coefs,il) != S_OK)
		return S_FALSE;

	if (size != (V_Coefs.size()/size))
	{	
		ERROR_MSG("Size of labels <> sqrt(Size Matrix)",pRet);
		return S_OK;
	}

	CCString Res;
	CCString curClass = LOCAL_PRICER_CORRMATRIX_CLASS;

	if( ICMLOCAL_CORRMATRIX (AsOf,
							 C_Name,
							 V_Labels,
							 V_Coefs, 
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


STDMETHODIMP ActiveXModule::ARM_Credit_ExtractCorrMatrix(VARIANT *pCorrMatrixId,
													 VARIANT *pLabels,
													 VARIANT *pRet)
{
	// TODO: Add your implementation code here
try
{
	ARM_result C_result;

	vector<string> V_Labels;
	VECTOR<double> V_Coefs;
	long size =0;
	long il = 0;

	CCString C_CorrId;

	if(VARIANT2CCString (*pCorrMatrixId, C_CorrId) != S_OK)
		return S_FALSE;

	long CorrMatrixId = LocalGetNumObjectId(C_CorrId);

	if (CorrMatrixId == ARM_KO)
	{
	ERROR_MSG("Invalid object CorrMatrix",pRet,ARM_ERROR_CORRMATRIX);
	return S_OK;
	}

	
	// if(VARIANT2VECTORCCSTRING(*pLabels,V_Labels,il) != S_OK)
	// 	return S_FALSE;
	VariantTools::convert(*pLabels,V_Labels); 

	size = V_Labels.size();

	if( ICMLOCAL_EXTRACTCORRMATRIX(CorrMatrixId,
						   		   V_Labels,
								   V_Coefs, 
								   C_result) == ARM_OK)
	{
	}
	else
	{
	ERROR_MSG(C_result.getMsg(),pRet);
	return S_OK;
	}

	//VECTORDOUBLE2MATRIXVARIANT(V_Coefs,size,pRet);
	VECTORDOUBLE2VARIANT(V_Coefs,pRet);

	_variant_t wrap_Res;
	wrap_Res.Attach(*pRet);
	*pRet=wrap_Res.Detach();

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


STDMETHODIMP ActiveXModule::ARM_Credit_GetLabel (VARIANT *pCurveId , 
										  	 VARIANT *pRet)
{
	// TODO: Add your implementation code here

try
{
	ARM_result C_result;

	CCString l_Curve;
	
	if(VARIANT2CCString (*pCurveId,l_Curve) != S_OK)
		return S_FALSE;

	CCString Res ;

	long Curveid = LocalGetNumObjectId (l_Curve);

	if (Curveid == ARM_KO)
	{
	ERROR_MSG("Invalid Curve",pRet,ARM_ERROR_DFCURVE);
	return S_OK;
	}

	if(ICMLOCAL_GetLabel(Curveid,
						 C_result) == ARM_OK)
	{
		Res = C_result.getMsg();
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


STDMETHODIMP ActiveXModule::ARM_Credit_SetLabel (VARIANT *pCurveId , 
											 VARIANT *pLabel, 
										  	 VARIANT *pRet)
{
	// TODO: Add your implementation code here

try
{
	ARM_result C_result;

	CCString C_Label;
	
	if(VARIANT2CCString (*pLabel,C_Label) != S_OK)
		return S_FALSE;

	CCString l_Curve;
	
	if(VARIANT2CCString (*pCurveId,l_Curve) != S_OK)
		return S_FALSE;

	CCString Res ;

	long Curveid = LocalGetNumObjectId (l_Curve);

	if (Curveid == ARM_KO)
	{
	ERROR_MSG("Invalid Curve",pRet,ARM_ERROR_DFCURVE);
	return S_OK;
	}

	if(ICMLOCAL_SetLabel(Curveid,
						 C_Label,
						 C_result) == ARM_OK)
	{
		Res = C_result.getMsg();
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


STDMETHODIMP ActiveXModule::ARM_Credit_Sensitivity(VARIANT *pPricer , 
											   VARIANT *pType,
											   VARIANT *pPlot,
											   VARIANT *pLabel,
											   VARIANT *pEpsilon,
											   double epsilonGamma,
											   VARIANT *pRet)
{
	// TODO: Add your implementation code here

	try
	{
		ARM_result C_result;
		_variant_t v_Epsilon(pEpsilon);
		v_Epsilon.ChangeType( VT_R8 );
		double C_Epsilon = (double) v_Epsilon;

		CCString l_Pricer;

		if(VARIANT2CCString (*pPricer, l_Pricer) != S_OK)
			return S_FALSE;

		double Res = 0.;

		long Pricerid = LocalGetNumObjectId (l_Pricer);

		if (Pricerid == ARM_KO)
		{
		ERROR_MSG("Invalid Pricer",pRet,ARM_ERROR_PRICER);
		return S_OK;
		}

		// CCString C_Type;
		// long crvtype;
		
		// if(VARIANT2CCString (*pType,C_Type) != S_OK)
		// 	return S_FALSE;
		std::string sensiTypeStr; VariantTools::convert(pType,sensiTypeStr); 
		qSENSITIVITY_TYPE crvtype ; 
		ICM_EnumsCnv::cnv(sensiTypeStr,crvtype ); 
		// if((crvtype = ICM_ConvCurvType (C_Type, C_result)) == ARM_DEFAULT_ERR)
		// 	return S_FALSE;

		CCString C_plot;

		if(VARIANT2CCString (*pPlot,C_plot) != S_OK)
			return S_FALSE;

		CCString C_Label;

		if(VARIANT2CCString (*pLabel,C_Label) != S_OK)
			return S_FALSE;

		long CurveId = LocalGetNumObjectId (C_Label);
		ARM_Object* curve = NULL;
		curve = (ARM_Object* )LOCAL_PERSISTENT_OBJECTS->GetObject(CurveId);
		std::string s_label("");
		if (curve)
		{	
			deducelabelforobject(*curve, s_label);
			
		} else {
			/*ERROR_MSG("Invalid CurveLabel",pRet,ARM_ERROR_PRICER);
			return S_OK;*/
			s_label = CCSTringToSTLString(C_Label);
		}

		if(ICMLOCAL_Sensitivity(Pricerid,
						(qSENSITIVITY_TYPE)crvtype,
						CCSTringToSTLString(C_plot),
						s_label,
						C_Epsilon,
						epsilonGamma,
						C_result) == ARM_OK)
		{
			Res =  C_result.getDouble();
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


STDMETHODIMP ActiveXModule::ARM_Credit_GenSchedule(VARIANT *pAccStartDate,
											   VARIANT *pAccEndDate,
											   VARIANT *pFixingFreq,
											   VARIANT *pDayCountFrq,
											   VARIANT *prefdate,
											   VARIANT *pCurrency,
											   VARIANT *ptypeDates,
											   VARIANT *pModFol,
											   VARIANT *pGapCredit,
											   VARIANT *pRet)
{
	// TODO: Add your implementation code here
try
{
	ARM_result C_result;
	VECTOR<double> DateOut;

	_variant_t v_AccStartDate (pAccStartDate);
	v_AccStartDate.ChangeType( VT_R8 );
	double C_AccStartDate = (double) v_AccStartDate;

	_variant_t v_AccEndDate (pAccEndDate);
	v_AccEndDate.ChangeType( VT_R8 );
	double C_AccEndDate = (double) v_AccEndDate;

	_variant_t v_GapCredit (pGapCredit);
	v_GapCredit.ChangeType( VT_R8 );
	double C_GapCredit = (double) v_GapCredit;

	CCString C_FixingFreq;
	if(VARIANT2CCString (*pFixingFreq, C_FixingFreq) != S_OK)
		return S_FALSE;

	CCString C_ModFol;
	if(VARIANT2CCString (*pModFol, C_ModFol) != S_OK)
		return S_FALSE;

	int l_FixingFreq;

	if((l_FixingFreq = ARM_ConvFrequency (C_FixingFreq, C_result)) == ARM_DEFAULT_ERR)
	{
		ERROR_MSG("Invalid FixingFreq",pRet,ARM_ERROR_FREQ);
		return S_OK;
	}

	CCString C_DayCountFrq;
	if(VARIANT2CCString (*pDayCountFrq, C_DayCountFrq) != S_OK)
		return S_FALSE;

	long l_DayCountFrq = ARM_ConvDayCount (C_DayCountFrq);

	long ModFollId;
	if((ModFollId = ARM_ConvFwdRule (C_ModFol, C_result)) == ARM_DEFAULT_ERR)
	{
		ERROR_MSG("ModFoll",pRet,ARM_ERROR_DAYCOUNT);
		return S_OK;
	}

	_variant_t v_refdate (prefdate);
	v_refdate.ChangeType( VT_R8 );
	double C_refdate = (double) v_refdate;

	if (C_refdate == -999.) C_refdate= -1.;

	CCString C_Currency;
	if(VARIANT2CCString (*pCurrency, C_Currency) != S_OK)
		return S_FALSE;

	CCString C_typeDates;
	if(VARIANT2CCString (*ptypeDates, C_typeDates) != S_OK)
		return S_FALSE;

	long l_typeDates;

	if((l_typeDates = ARM_ConvTypeGenDates (C_typeDates, C_result)) == ARM_DEFAULT_ERR)
	{
		ERROR_MSG("TypeDate",pRet,ARM_ERROR_DAYCOUNT);
		return S_OK;
	}


	if(ICMLOCAL_GenSchedule(C_AccStartDate,
							 C_AccEndDate,
							 C_refdate,
							 l_FixingFreq,
							 l_DayCountFrq,
							 C_Currency,
							 l_typeDates,
							 ModFollId,
							 (int) C_GapCredit,
							 C_result) == ARM_OK)
	{
		VECTOR<CCString> dDate;
		long vecSize;
		ExtractVectorDateFromFile("123",dDate,vecSize);
		DateOut.clear();

		int nbrows = vecSize;
		int nbcolumns = 1;

		for(int i = 0; i < nbrows; i++)
			DateOut.push_back(Local_ARMDATE2XLDATE(dDate[i]));
	}
	else
	{
	ERROR_MSG(C_result.getMsg(),pRet);
	return S_OK;
	}

	VECTORDOUBLE2VARIANT(DateOut,pRet);

	_variant_t wrap_Res;
	wrap_Res.Attach(*pRet);
	*pRet=wrap_Res.Detach();

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


STDMETHODIMP ActiveXModule::ARM_Credit_CashFlows(VARIANT *pCoefs, 
											  VARIANT *pRet)
{
	// TODO: Add your implementation code here
try
{
	ARM_result C_result;

	VECTOR<CCString> V_Coefs;
	long size = 0;

	if(VARIANT2VECTORCCSTRING(*pCoefs,V_Coefs,size) != S_OK)
		return S_FALSE;

	CCString Res;
	CCString curClass = LOCAL_PRICER_CORRMATRIX_CLASS;

	if( ICMLOCAL_CashFlows(V_Coefs,
						   size/9,
						   9,	
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


STDMETHODIMP ActiveXModule::ARM_Credit_Parameters(VARIANT *pCoefs, 
											  long	   nbcols,
											  VARIANT *pRet)
{
	// TODO: Add your implementation code here
try
{
	ARM_result C_result;

	VECTOR<CCString> V_Coefs;
	long size = 0;

	if(VARIANT2VECTORCCSTRING(*pCoefs,V_Coefs,size) != S_OK)
		return S_FALSE;

	CCString Res;
	CCString curClass = LOCAL_PRICER_CORRMATRIX_CLASS;

	if( ICMLOCAL_Parameters(V_Coefs,
						    size/nbcols,
						    nbcols,	
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


long LocalDisplayArmFile (const CCString& fileName)
{
	CCString	vEditor(EDITOR);
	int	vLen = fileName.GetLen();
	if(vLen > 3)
	{
		char*	vExtension = fileName.c_str() + vLen - 3;
		CCString	cExtension(vExtension);
		if(cExtension == "xml")
		{
			vEditor = XML_EDITOR;	
		}
	}

	/// the start command helps to start on a new thread the editor!
	CCString command = CCString("start ") + vEditor + " " + fileName;

	//_flushall();

	system((const char*) command);

	return(ARM_OK);
}


long LocalGetArmViewFile (const CCString& sockId)
{

    CCString clientViewFileName = CCString(VIEW_FILE_CLIENT_LOCATION)+CCString(VIEW_FILE_PREFIX)+sockId;

	LocalDisplayArmFile (clientViewFileName);

	/// Watch out! because we start in a new thread
	/// we cannot delete the file anymore!!!!
	/// otherwise we will delete while it is processed

	return ARM_OK;
}

long Local_ViewFile (const CCString& C_instId)
{
	ARM_result C_result;
	
	long retCode = ARMLOCAL_ARM_View (LocalGetNumObjectId (C_instId), C_result);

	if(retCode == ARM_OK)
	{
		char username[100];
		DWORD nbChar = sizeof(username);
		GetUserName(username,&nbChar);

		retCode = LocalGetArmViewFile ((CCString)"123" + (CCString)username);
	}

	return retCode;
}

STDMETHODIMP ActiveXModule::ARM_View (VARIANT *pObjet,
								  VARIANT *pRet)
{
	try 
	{
		if (!pObjet) 
			ICMTHROW(ERR_INVALID_ARGUMENT,"ActiveXModule::ARM_View: null argument"); 

		CCString C_Label;		
		if(VARIANT2CCString (*pObjet,C_Label) != S_OK) return S_FALSE;

		long retCode = Local_ViewFile (C_Label);

	}
	catch(...)
	{
		return E_FAIL ;
	}
	return S_OK;
}


STDMETHODIMP ActiveXModule::ARM_Credit_Version(VARIANT *pRet)
{
	// TODO: Add your implementation code here

try
{
	VariantTools::convert(ICMLOCAL_Version(),*pRet); 
	return S_OK; 	
} 
catch(...)
{
	return E_FAIL;
}
}

/** 
STDMETHODIMP ActiveXModule::ARM_Credit_CorrelationSmile(VARIANT *pPricerId,
													double pSmileType,
													double pMktPrice,
													double pSeed,
													double pUpfrontPay,
													int datatype,
													VARIANT *pRet)	
{
	// TODO: Add your implementation code here
try
{
	CCString C_pricerId;

	ARM_result C_result;
	double Res = 0;

	if(VARIANT2CCString (*pPricerId, C_pricerId) != S_OK)
		return S_FALSE;

	if (ICMLOCAL_CorrelationSmile (LocalGetNumObjectId (C_pricerId), 
									pSmileType,
									pMktPrice,
									pSeed,
									pUpfrontPay,
									datatype,
									C_result) == ARM_OK)


	{
		Res = C_result.getDouble();
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
**/ 
STDMETHODIMP ActiveXModule::ARM_Credit_CptBaseCorrelation(double AsOf,
															BSTR parameterId,
															BSTR CalMethod,
															BSTR IndexId,
															VARIANT *pStrikeLow,
															VARIANT *pStrikeHigh,
															VARIANT *pVMktBid,
															VARIANT *pVMktAsk,
															VARIANT *pVUpfBid,
															VARIANT *pVUpfAsk,
															VARIANT *pVInitialCorrel,
															VARIANT *pVDeltaLevrage,
															BSTR ModelId,
															double integrationStep,
															double lagStartDate,
															double creditLag,
															VARIANT *pVectorPrevIndexId,
															VARIANT *pMatrixPrevBC,
															double step,
															BSTR CalMeth,
															VARIANT *pRet)

{
	
	CCString C_pricerId;
	
	ARM_result C_result;
	VECTOR<double> correlOut;
try {
	//
	int l_CalMethod =0;
	long size = 0;

	_bstr_t bvar(parameterId);
	CCString C_parameterId = bvar;

	_bstr_t bvar1(CalMethod);
	CCString C_CalMethod = bvar1;

	_bstr_t bvar2(IndexId);
	CCString C_IndexId = bvar2;

	_bstr_t bvar3(ModelId);
	CCString C_ModelId = bvar3;

	_bstr_t bvar5(CalMeth);
	CCString C_CalMeth = bvar5;

	VECTOR<double> C_VStrikeLow; 
	VECTOR<double> C_VStrikeHigh; 	
	VECTOR<double> C_VMktBid; 
	VECTOR<double> C_VMktAsk; 
	VECTOR<double> C_VUpfBid; 	
	VECTOR<double> C_VUpfAsk; 
	// from optional parameters 
	VECTOR<double> C_VInitialCorrel;
	VECTOR<double> C_VDeltaLevrage; 
	vector<double> DefaultValueDbl;
	DefaultValueDbl.clear();DefaultValueDbl.empty();


	VECTOR<CCString> C_PrevIndexId;

	VECTOR<double> C_PrevMatrixBC;

	VECTOR<double> basesCorrelation;
	

	int l_CalMeth =0;

	if(VARIANT2VECTORDOUBLE(*pStrikeLow,C_VStrikeLow,size) != S_OK)
		return S_FALSE;
	if(VARIANT2VECTORDOUBLE(*pStrikeHigh,C_VStrikeHigh,size) != S_OK)
		return S_FALSE;
	if(VARIANT2VECTORDOUBLE(*pVMktBid,C_VMktBid,size) != S_OK)
		return S_FALSE;
	if(VARIANT2VECTORDOUBLE(*pVMktAsk,C_VMktAsk,size) != S_OK)
		return S_FALSE;

	if (pVInitialCorrel->pvarVal->vt == 0 ) {
		C_VInitialCorrel.clear();
		C_VInitialCorrel.empty();
	} else {
		if(VARIANT2VECTORDOUBLE (*pVInitialCorrel, C_VInitialCorrel,size) != S_OK)
			return S_FALSE;
	}
	if (pVDeltaLevrage->pvarVal->vt == 0 ) {
		C_VDeltaLevrage.clear();
		C_VDeltaLevrage.empty();
	} else {
		if(VARIANT2VECTORDOUBLE (*pVDeltaLevrage, C_VDeltaLevrage,size) != S_OK)
			return S_FALSE;
	}


	if(VARIANT2VECTORDOUBLE(*pVUpfBid,C_VUpfBid,size) != S_OK)
		return S_FALSE;
	if(VARIANT2VECTORDOUBLE(*pVUpfAsk,C_VUpfAsk,size) != S_OK)
		return S_FALSE;
	if (pMatrixPrevBC->pvarVal->vt == 0 ) {
		C_PrevMatrixBC.clear();
		C_PrevMatrixBC.empty();
	} else {
		if(VARIANT2VECTORDOUBLE (*pMatrixPrevBC, C_PrevMatrixBC,size) != S_OK)
			return S_FALSE;
	}
	
	if (pVectorPrevIndexId->pvarVal->vt == 0 ) {
		C_PrevIndexId.clear();
		C_PrevIndexId.empty();
	} else {
		if(VARIANT2VECTORCCSTRING (*pVectorPrevIndexId, C_PrevIndexId,size) != S_OK)
			return S_FALSE;
	}
	VECTOR<long> V_PrevIndexId;V_PrevIndexId.resize(C_PrevIndexId.size());
	for (int kk=0; kk<C_PrevIndexId.size();kk++)
	{	
		V_PrevIndexId[kk]=LocalGetNumObjectId(C_PrevIndexId[kk]);
	}
	long l_params=LocalGetNumObjectId(C_parameterId);
	
	if((l_CalMethod = ARM_ConvCalibType (C_CalMethod, C_result)) == ARM_DEFAULT_ERR)
	{
		return E_FAIL;
	}

	if((l_CalMeth = ARM_ConvCalibMeth (C_CalMeth, C_result)) == ARM_DEFAULT_ERR)
	{
		return E_FAIL;
	}
	
	if (ICMLOCAL_CPT_BASE_CORRELATION(AsOf,
								//C_Name,
								l_CalMethod,
								LocalGetNumObjectId(C_IndexId),
								C_VStrikeLow,
								C_VStrikeHigh,
								C_VMktBid,
								C_VMktAsk,
								C_VUpfBid,
								C_VUpfAsk,
								C_VInitialCorrel,
								C_VDeltaLevrage,
								basesCorrelation,
								LocalGetNumObjectId(C_ModelId),
								(int)integrationStep,
								(int)lagStartDate,
								(int)creditLag,
								V_PrevIndexId,
								C_PrevMatrixBC,
								step,
								(qOPTIMIZE_TYPE)l_CalMeth,
								l_params,
								C_result) == ARM_OK)


	{
		VECTOR<double> correl;
		long vecSize;
		ExtractVectorDoubleFromFile("123",correl,vecSize);
		correlOut.clear();

		int nbrows = vecSize;
		int nbcolumns = 1;

		for(int i = 0; i < nbrows; i++)
			correlOut.push_back(correl[i]);
	}
	else
	{
		ERROR_MSG(C_result.getMsg(),pRet);
		return S_OK;
	}

	VECTORDOUBLE2VARIANT(correlOut,pRet);
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

 
/**
STDMETHODIMP ActiveXModule::ARM_Credit_GetBaseCorrelation(VARIANT *pPricerId,
													  double pEquityAmount,
													  VARIANT *pMktSpreads,
													  double pSmileType,
													  VARIANT *pSeeds,
													  VARIANT *pUpfronts,
													  VARIANT *pRet)
{
	// TODO: Add your implementation code here
try
{
	CCString C_pricerId;
	
	ARM_result C_result;
	double Res = 0;

	if(VARIANT2CCString (*pPricerId, C_pricerId) != S_OK)
		return S_FALSE;

	VECTOR<double> V_MktSpreads;
	long size = 0;

	if(VARIANT2VECTORDOUBLE(*pMktSpreads,V_MktSpreads,size) != S_OK)
		return S_FALSE;

	VECTOR<double> V_Seeds;
	VECTOR<double> V_Upfronts;

	if(VARIANT2VECTORDOUBLE(*pSeeds,V_Seeds,size) != S_OK)
		return S_FALSE;

	if(VARIANT2VECTORDOUBLE(*pUpfronts,V_Upfronts,size) != S_OK)
		return S_FALSE;

	if (ICMLOCAL_GetBaseCorrelation(LocalGetNumObjectId (C_pricerId), 
									pEquityAmount,
									V_MktSpreads,
									pSmileType,
									V_Seeds,
									V_Upfronts,
									C_result) == ARM_OK)


	{
		Res = C_result.getDouble();
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
	catch(...)
{
	return E_FAIL;
}
}
**/ 

/**
STDMETHODIMP ActiveXModule::ARM_Credit_GetImpliedCorrelation(VARIANT *pPricerId,
													  double pRange,
													  VARIANT *pMktSpreads,
													  double pSmileType,
													  VARIANT *pSeeds,
													  VARIANT *pUpfronts,
													  VARIANT *pRet)
{
	// TODO: Add your implementation code here
try
{
	CCString C_pricerId;
	
	ARM_result C_result;
	double Res = 0;

	if(VARIANT2CCString (*pPricerId, C_pricerId) != S_OK)
		return S_FALSE;

	VECTOR<double> V_MktSpreads;
	long size = 0;

	if(VARIANT2VECTORDOUBLE(*pMktSpreads,V_MktSpreads,size) != S_OK)
		return S_FALSE;

	VECTOR<double> V_Seeds;
	VECTOR<double> V_Upfronts;

	if(VARIANT2VECTORDOUBLE(*pSeeds,V_Seeds,size) != S_OK)
		return S_FALSE;

	if(VARIANT2VECTORDOUBLE(*pUpfronts,V_Upfronts,size) != S_OK)
		return S_FALSE;

	if (ICMLOCAL_GetImpliedCorrelation(LocalGetNumObjectId (C_pricerId), 
									pRange,
									V_MktSpreads,
									pSmileType,
									V_Seeds,									
									V_Upfronts,
									C_result) == ARM_OK)


	{
		Res = C_result.getDouble();
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
	catch(...)
{
	return E_FAIL;
}
}
**/ 
STDMETHODIMP ActiveXModule::ARM_Credit_GetDuration (VARIANT *pPricer,VARIANT *pRet)
{
	// TODO: Add your implementation code here
try
{
	ARM_result C_result;
	CCString l_Pricer;
	
	if(VARIANT2CCString (*pPricer,l_Pricer) != S_OK)
		return S_FALSE;

	double Res =0.;
	long Pricerid = LocalGetNumObjectId (l_Pricer);

	if (Pricerid == ARM_KO)
	{	ERROR_MSG("Invalid Pricer",pRet,ARM_ERROR_DFCURVE);
	return S_OK;
	}


	if(ICMLOCAL_GetDuration(Pricerid,C_result) == ARM_OK)
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

STDMETHODIMP ActiveXModule::ARM_Credit_SetCorrelationMatrix (BSTR pModelMultiCurves, 
														 BSTR pCorrMatrix, 
														 BSTR *pRet)
{
	// TODO: Add your implementation code here

try
{
	ARM_result C_result;

	_bstr_t bModelMultiCurves(pModelMultiCurves);
	CCString l_MMC = bModelMultiCurves;

	_bstr_t bCorrMatrix(pCorrMatrix);
	CCString l_CORR = bCorrMatrix;

	CCString retour;

	double Res =0.;

	long MMCid = LocalGetNumObjectId (l_MMC);

	if (MMCid == ARM_KO)
	{
	retour = "Invalid Model Multi Curves";
	_bstr_t tmpChaine = retour;
	*pRet = SysAllocString(tmpChaine);
	return S_OK;
	}

	long Corrid = LocalGetNumObjectId (l_CORR);

	if (Corrid == ARM_KO)
	{
	retour = "Invalid Correlation Matrix";
	_bstr_t tmpChaine = retour;
	*pRet = SysAllocString(tmpChaine);
	return S_OK;
	}

	if(ICMLOCAL_SetCorrelationMatrix(MMCid,Corrid,C_result) != ARM_OK)
		return E_FAIL;

	retour = "SetCorrelationMatrix OK";
	_bstr_t tmpChaine = retour;
	*pRet = SysAllocString(tmpChaine);

	return S_OK;
}
catch(...)
{
	return E_FAIL;
}
}


STDMETHODIMP ActiveXModule:: ARM_Credit_CloneCorrMatrixBary(BSTR pCorrMatrix, 
														double Beta,
														int UpOrDown,
														BSTR *pRet)
{

try
{

	ARM_result C_result;
	CCString Res = "ERROR";

	_bstr_t bCorrMatrix(pCorrMatrix);
	CCString C_CORR = bCorrMatrix;
	
	CCString curClass = LOCAL_PRICER_CLASS;

	if (ICMLOCAL_CloneCorrMatrixBary(LocalGetNumObjectId (C_CORR),
									Beta,
									UpOrDown,
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


STDMETHODIMP ActiveXModule::ARM_Credit_FwdSpreadPricer(VARIANT *pPricer , 
												   double  Mty1,
												   double  Mty2,
												   VARIANT *pRet)
{
	// TODO: Add your implementation code here

try
{
	ARM_result C_result;
	ARM_Date matu1,matu2; 
	VariantTools::convertXLDate(Mty1, matu1);
	VariantTools::convertXLDate(Mty2, matu2);


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

	if(ICMLOCAL_FwdSpreadPricer(Pricerid,
								matu1,
								matu2,
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
catch(...)
{
	return E_FAIL;
}
}


STDMETHODIMP ActiveXModule::ARM_Credit_ImpliedVol (VARIANT *pPricer , 
											   double pMktPrice,
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

	if(ICMLOCAL_ImpliedVol(Pricerid,
						   pMktPrice,
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
catch(...)
{
	return E_FAIL;
}
}



STDMETHODIMP ActiveXModule::ARM_Credit_VirtualCdsSpread (VARIANT *pPricer,
													 VARIANT *pMaturity,
													 VARIANT *pRet)
{
	// TODO: Add your implementation code here

	try
		{
		
		ARM_result C_result;

		CCString l_Pricer;

		_variant_t v_Maturity (pMaturity);
		v_Maturity.ChangeType( VT_R8 );
		double C_Maturity = (double) v_Maturity;
		
		if(VARIANT2CCString (*pPricer, l_Pricer) != S_OK)
			return S_FALSE;
		
		double Res =0.;
		
		long Pricerid = LocalGetNumObjectId (l_Pricer);
		
		
		if (Pricerid == ARM_KO)
		{
			ERROR_MSG("Invalid Pricer",pRet,ARM_ERROR_PRICER);
			return S_OK;
		}

		if(ICMLOCAL_VirtualCdsSpread(Pricerid,
									 C_Maturity,
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
	catch(...)
	{
		return E_FAIL;
	}
}


STDMETHODIMP ActiveXModule::ARM_Credit_BSGreeks(VARIANT *pPricer,
											VARIANT *pGreekType,
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
		
		CCString C_GreekType;
		long greektype;
		
		if(VARIANT2CCString (*pGreekType,C_GreekType) != S_OK)
			return S_FALSE;

		if((greektype = ICM_ConvGreekType (C_GreekType, C_result)) == ARM_DEFAULT_ERR)
			return S_FALSE;

		
		if(ICMLOCAL_BSGreeks(Pricerid,
							 greektype,
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
	catch(...)
	{
		return E_FAIL;
	}
}



STDMETHODIMP ActiveXModule::ARM_Credit_SetCorrelation (BSTR pModelMultiCurvesId , 
												   BSTR pCorrelationId,
												   VARIANT *pRet)
{
	// TODO: Add your implementation code here

try
{
	ARM_result C_result;

	_bstr_t b_ModelMultiCurvesId(pModelMultiCurvesId);
	CCString C_ModelMultiCurvesId = b_ModelMultiCurvesId;

	_bstr_t b_CorrelationId(pCorrelationId);
	CCString C_CorrelationId = b_CorrelationId;

	double Res =0.;

	long ModelMultiCurvesId = LocalGetNumObjectId (C_ModelMultiCurvesId);

	if (ModelMultiCurvesId == ARM_KO)
	{
	ERROR_MSG("Invalid Model",pRet,ARM_ERROR_PRICER);
	return S_OK;
	}

	long CorrelationId = LocalGetNumObjectId (C_CorrelationId);

	if (CorrelationId == ARM_KO)
	{
	ERROR_MSG("Invalid Correlation",pRet,ARM_ERROR_PRICER);
	return S_OK;
	}

	if(ICMLOCAL_SetCorrelation(ModelMultiCurvesId,
							   CorrelationId,
							   C_result) == ARM_OK)
	{
		Res = 1.;
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
catch(...)
{
	return E_FAIL;
}
}


STDMETHODIMP ActiveXModule::ARM_Credit_CorrelationStrike(VARIANT* pLabels, 
													 VARIANT* pVolCurves, 
													 VARIANT* pProportions,
													 VARIANT* pSmileStrikeLow,
													 VARIANT* pSmileStrikeHigh,
													 VARIANT* pVIndex,
													 double AsOf_,
													 BSTR Name,
													 BSTR *pRet)
{
try
{
	ARM_result C_result;

	CCString Res;
	CCString curClass = LOCAL_PRICER_CORRMATRIX_CLASS;

	long size = 0;
	int i=0;

	std::string C_Name ;
	VariantTools::convert(Name,C_Name); 
	

	vector<string> vLabels;
	VECTOR<CCString> vVolCurves;
	VECTOR<long> LVolCurves;
	VECTOR<double> vProportions;
	VECTOR<double> vSmileStrikeLow;
	VECTOR<double> vSmileStrikeHigh;

	VECTOR<CCString> vIndex;
	VECTOR<long> LIndex;LIndex.clear();

	ARM_Date AsOf; 
	VariantTools::convertXLDate(AsOf_,AsOf); 
	// if(VARIANT2VECTORCCSTRING(*pLabels,vLabels,size) != S_OK)
	// 	return S_FALSE;
	VariantTools::convert(*pLabels,vLabels); 

	if(VARIANT2VECTORCCSTRING(*pVolCurves,vVolCurves,size) != S_OK)
		return S_FALSE;

	LVolCurves.resize(size);

	for (i=0; i<size; i++)
		LVolCurves[i] = LocalGetNumObjectId(vVolCurves[i]);

	size = 0;

	if(VARIANT2VECTORCCSTRING(*pVIndex,vIndex,size) != S_OK)
		return S_FALSE;

	LIndex.resize(size);

	for (i=0; i<size; i++)
		LIndex[i] = LocalGetNumObjectId(vIndex[i]);

	if(VARIANT2VECTORDOUBLE(*pProportions,vProportions,size) != S_OK)
		return S_FALSE;

	if(VARIANT2VECTORDOUBLE(*pSmileStrikeLow,vSmileStrikeLow,size) != S_OK)
		return S_FALSE;

	if(VARIANT2VECTORDOUBLE(*pSmileStrikeHigh,vSmileStrikeHigh,size) != S_OK)
		return S_FALSE;

	ICM_QMatrix<double> fullSmileLow, fullSmileUp ; 

	if(ICMLOCAL_CORRELATION_STRIKE(AsOf,
								   C_Name,	
								   vLabels,
								   LVolCurves,
								   vProportions,
								   vSmileStrikeLow,
								   vSmileStrikeHigh,
								   fullSmileLow,
								   fullSmileUp,
								   LIndex,
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

STDMETHODIMP ActiveXModule::ARM_Credit_CorrelationSmileStrike(VARIANT* pLabels, 
													 VARIANT* pVolCurves, 
													 VARIANT* pProportions,
													 double asOfDate,
													 VARIANT* pSmileStrikeLow,
													 VARIANT* pSmileStrikeHigh,
													 VARIANT* pVIndex,
													 BSTR Name,
													 VARIANT * pFullStrikeLow,
													 VARIANT * pFullStrikeUp,
													 BSTR *pRet)
{
try
	{
	ARM_result C_result;

	long size = 0;
	int i=0;
	_bstr_t b_Name(Name);
	CCString C_Name = b_Name;

	VECTOR<string> vLabels;
	VECTOR<string> vVolCurves;
	VECTOR<long> LVolCurves;
	ARM_Vector vProportions;
	ARM_Vector vSmileStrikeLow;
	ARM_Vector vSmileStrikeHigh;
	ICM_QMatrix<double> fullStrikeLow,fullStrikeUp; 

	VECTOR<string> vIndex;
	VECTOR<long> LIndex;LIndex.clear();

	VariantTools::convert(*pLabels,vLabels);
	VariantTools::convert(*pVolCurves,vVolCurves);
	VariantTools::convert(*pProportions,vProportions);
	//Index
	if (pVIndex->vt != VT_ERROR && pVIndex->pvarVal->vt != VT_EMPTY )
		VariantTools::convert(*pVIndex,vIndex);
	// or SmileStrikeLow/High
	if (pSmileStrikeLow->vt != VT_ERROR && pSmileStrikeLow->pvarVal->vt != VT_EMPTY && pSmileStrikeLow->pvarVal->vt != VT_NULL)
		VariantTools::convert(*pSmileStrikeLow,vSmileStrikeLow);
	if (pSmileStrikeHigh->vt != VT_ERROR && pSmileStrikeHigh->pvarVal->vt != VT_EMPTY && pSmileStrikeHigh->pvarVal->vt != VT_NULL)	
		VariantTools::convert(*pSmileStrikeHigh,vSmileStrikeHigh);
	// or FullSmileStrikeLow/High
	if (pFullStrikeLow->vt != VT_ERROR && pFullStrikeLow->pvarVal->vt != VT_EMPTY && pFullStrikeLow->pvarVal->vt != VT_NULL)
		VariantTools::convert(*pFullStrikeLow, fullStrikeLow);
	if (pFullStrikeUp->vt != VT_ERROR && pFullStrikeUp->pvarVal->vt != VT_EMPTY && pFullStrikeUp->pvarVal->vt != VT_NULL)
		VariantTools::convert(*pFullStrikeUp, fullStrikeUp);
	
	size = vVolCurves.size();
	LVolCurves.resize(size);
	for (i=0; i<size; i++)
		LVolCurves[i] = LocalPersistent::get().getObjectId(vVolCurves[i]);
	
	if ( vIndex.size() == size) {
		LIndex.resize(size);
		for (i=0; i<size; i++)
			LIndex[i] = LocalPersistent::get().getObjectId(vIndex[i]);
	}

	long newId = ICMLOCAL_CORRELATION_SMILE_STRIKE(asOfDate,
										   (string)C_Name.c_str(),	
										   vLabels,
										   LVolCurves,
										   ARM_Vector(vProportions),
										   ARM_Vector(vSmileStrikeLow),
										   ARM_Vector(vSmileStrikeHigh),
										   fullStrikeLow,
										   fullStrikeUp,
										   LIndex);
	string newLabel = LocalPersistent::get().getStringId(newId, LOCAL_PRICER_CORRMATRIX_CLASS);
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
	 return E_FAIL ;
}



STDMETHODIMP ActiveXModule::ARM_Credit_Beta_Correlation(VARIANT *pLabels,
													VARIANT *pCoefs, 
													double AsOf_,
													BSTR Name,
													BSTR index1_,
													BSTR index2_,
													VARIANT *pRet)
{
	// TODO: Add your implementation code here
try
{
	ARM_result C_result;

	std::vector<std::string> V_Labels;
	VECTOR<double> V_Coefs;
	long size ;
	long il = 0;

	// _bstr_t b_Name(Name);
	// CCString C_Name = b_Name;
	std::string C_Name; 
	VariantTools::convert(Name,C_Name); 

	// if(VARIANT2VECTORCCSTRING(*pLabels,V_Labels,il) != S_OK)
	// 	return S_FALSE;
	VariantTools::convert(*pLabels,V_Labels); 

	size = V_Labels.size();

	if (!size)
	{	
		ERROR_MSG("Empty labels",pRet);
		return S_OK;
	}

	/*if(VARIANT2VECTORDOUBLE (*pCoefs,V_Coefs,il) != S_OK)
		return S_FALSE;*/
	VariantTools::convert(*pCoefs,V_Coefs);

	if (size != V_Coefs.size())
	{	
		ERROR_MSG("Size of labels <> sqrt(Size Matrix)",pRet);
		return S_OK;
	}

	std::string index1; VariantTools::convert(index1_,index1); 
	std::string index2; VariantTools::convert(index2_,index2); 
	long idIndex1 = LocalPersistent::get().getObjectId(index1); 
	long idIndex2 = LocalPersistent::get().getObjectId(index2); 

	CCString Res;
	CCString curClass = LOCAL_PRICER_CORRMATRIX_CLASS;

	ARM_Date AsOf ;
	VariantTools::convertXLDate(AsOf_,AsOf); 
	if( ICMLOCAL_BETA_CORRELATION (AsOf,
								 C_Name,	
								 V_Labels,
								 V_Coefs, 
								 idIndex1,
								 idIndex2,
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


STDMETHODIMP ActiveXModule::ARM_Credit_Flat_Correlation(DATE AsOf_, 
														  BSTR structName_,
														  double correlValue_,
														  BSTR idIndex1_,
														  BSTR idIndex2_,
														  VARIANT *pRet)
{
	// TODO: Add your implementation code here
try
	{
		std::string structName; VariantTools::convert(structName_,structName); 
		ARM_Date AsOf ; VariantTools::convertXLDate(AsOf_,AsOf); 
		std::string index1; VariantTools::convert(idIndex1_,index1); 
		std::string index2; VariantTools::convert(idIndex2_,index2); 
		long idIndex1 = LocalPersistent::get().getObjectId(index1); 
		long idIndex2 = LocalPersistent::get().getObjectId(index2); 
		
		long oldId=-1; 
		long newId =  ICMLOCAL_FLAT_CORRELATION(AsOf,structName,correlValue_,idIndex1,idIndex1,oldId); 

		std::string newName = LocalPersistent::get().getStringId(newId,LOCAL_PRICER_CORRMATRIX_CLASS); 
		VariantTools::convert(newName,*pRet); 
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

STDMETHODIMP ActiveXModule::ARM_Credit_GetEqStrikeDown(BSTR correlId,
												   BSTR indexname, 
													double *pRet)
{
	// TODO: Add your implementation code here

try
{
	ARM_result C_result;

	_bstr_t b_correlId(correlId);
	CCString C_correlId = b_correlId;

	_bstr_t b_indexname(indexname);
	CCString C_indexname = b_indexname;

	double Res =0.;

	long Pricerid = LocalGetNumObjectId (C_correlId);

	if(ICMLOCAL_GetEqStrikeDown(Pricerid,
								C_indexname,
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

STDMETHODIMP ActiveXModule::ARM_Credit_GetEqStrikeUp(BSTR correlId,
												 BSTR indexname, 
												 double *pRet)
{
	// TODO: Add your implementation code here

try
{
	ARM_result C_result;

	_bstr_t b_correlId(correlId);
	CCString C_correlId = b_correlId;

	_bstr_t b_indexname(indexname);
	CCString C_indexname = b_indexname;

	double Res =0.;

	long Pricerid = LocalGetNumObjectId (C_correlId);

	if(ICMLOCAL_GetEqStrikeUp(Pricerid,
							C_indexname,
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

STDMETHODIMP ActiveXModule::ARM_Credit_GetCorrelStrikeDown(BSTR correlId,
													   double yfmaturity,
													   double *pRet)
{
	// TODO: Add your implementation code here

try
{
	ARM_result C_result;

	_bstr_t b_correlId(correlId);
	CCString C_correlId = b_correlId;

	double Res =0.;

	long Pricerid = LocalGetNumObjectId (C_correlId);

	if(ICMLOCAL_GetCorrelStrikeDown(Pricerid,
									yfmaturity,
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


STDMETHODIMP ActiveXModule::ARM_Credit_GetCorrelStrikeUp(BSTR correlId,
													 double yfmaturity,
													 double *pRet)
{
	// TODO: Add your implementation code here

try
{
	ARM_result C_result;

	_bstr_t b_correlId(correlId);
	CCString C_correlId = b_correlId;

	double Res =0.;

	long Pricerid = LocalGetNumObjectId (C_correlId);

	if(ICMLOCAL_GetCorrelStrikeUp(Pricerid,
								  yfmaturity,
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


STDMETHODIMP ActiveXModule::ARM_Credit_GetCorrelation(BSTR ModelId,
												  BSTR *pRet)
{
try
{
	ARM_result C_result;

	_bstr_t b_ModelId(ModelId);
	CCString C_ModelId = b_ModelId;

	CCString Res;
	CCString curClass = LOCAL_PRICER_CORRMATRIX_CLASS;

	if(ICMLOCAL_GetCorrelation(LocalGetNumObjectId(C_ModelId),
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


STDMETHODIMP ActiveXModule::ARM_NextCpnDate(double AsOfDate,
										double maturity,
										BSTR frequency,
										BSTR rule,
										BSTR currency,
										BSTR intrule,
 									    double *pRet)
{
try
{
	ARM_result C_result;

	_bstr_t b_frequency(frequency);
	CCString C_frequency = b_frequency;

	_bstr_t b_rule(rule);
	CCString C_rule = b_rule;

	_bstr_t b_currency(currency);
	CCString C_currency = b_currency;

	_bstr_t b_intrule(intrule);
	CCString C_intrule = b_intrule;

	long frequencyId;
	long ruleId;
	long intruleId;

	if(C_currency == "DEFAULT")
	{
		ARM_result ccyres;
		ARMLOCAL_ARM_GetDefaultCurrency (ccyres);
		if(ccyres.getRetCode () != ARM_OK)
		{
			return E_FAIL;
		}
		else
		{
			C_currency = ccyres.getString ();
		}
	}

	if((ruleId = ARM_ConvRule (C_rule, C_result)) == ARM_DEFAULT_ERR)
	{
		return E_FAIL;
	}
	
	if((frequencyId = ARM_ConvFrequency (C_frequency, C_result)) == ARM_DEFAULT_ERR)
	{
		return E_FAIL;
	}

	intruleId = ARM_ConvIntRule (C_intrule);

	double Res;

	if(ARMLOCAL_ARM_NextCpnDate (AsOfDate,
											 maturity,
											 frequencyId,
											 ruleId,
											 C_currency,
											 C_result,
											 intruleId) == ARM_OK)
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


STDMETHODIMP ActiveXModule::ARM_Credit_SetProportionsInfos(BSTR correlId,
													   BSTR IndexName,
													   double proportion,
													   double forcedstrikelow,
													   double forcedstrikehigh)
{
	// TODO: Add your implementation code here

try
{
	ARM_result C_result;

	_bstr_t b_correlId(correlId);
	CCString C_correlId = b_correlId;

	_bstr_t b_IndexName(IndexName);
	CCString C_IndexName = b_IndexName;

	double Res =0.;

	long CorrelId = LocalGetNumObjectId (C_correlId);

	if(ICMLOCAL_SetProportionsInfos(CorrelId,
									C_IndexName,
									proportion,
									forcedstrikelow,
									forcedstrikehigh,
									C_result) == ARM_OK)
	{
		double pRet = C_result.getDouble ();
	}

	return S_OK;
}
catch(...)
{
	return E_FAIL;
}
}


STDMETHODIMP ActiveXModule::ARM_Credit_CptImplCvForCDO2(BSTR pricerId,
													BSTR Name,
													BSTR Tenor,
													double *pRet)
{
	// TODO: Add your implementation code here
try
{
	ARM_result C_result;

	long il = 0;

	_bstr_t b_Tenor(Tenor);
	CCString C_Tenor = b_Tenor;

	_bstr_t b_Name(Name);
	CCString C_Name = b_Name;

	_bstr_t b_pricerId(pricerId);
	CCString C_pricerId = b_pricerId;
	long l_pricerId = LocalGetNumObjectId (C_pricerId);

	double result = 0.;
	CCString Res;
	CCString curClass = LOCAL_PRICER_CORRMATRIX_CLASS;

	if( ICMLOCAL_ComputeImplicitCurveForCDO2 (l_pricerId,
								 C_Name,	
								 C_Tenor, 
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


STDMETHODIMP ActiveXModule::ARM_Credit_AddPeriod(double pAsOf,
										BSTR Maturity,
										BSTR pCcy,
										BSTR AdjRule,
										BSTR AdjCDS,
 									    VARIANT *pRet)
{
try
{
	ARM_result C_result;

	_bstr_t b_Maturity(Maturity);
	CCString C_Maturity = b_Maturity;

	_bstr_t b_Currency(pCcy);
	CCString C_Currency = b_Currency;

	_bstr_t b_AdjRule(AdjRule);
	CCString C_AdjRule = b_AdjRule;
	
	_bstr_t b_AdjCDS(AdjCDS);
	CCString C_AdjCDS = b_AdjCDS;


	if(C_Currency == "DEFAULT")
	{
		ARM_result ccyres;
		ARMLOCAL_ARM_GetDefaultCurrency (ccyres);
		if(ccyres.getRetCode () != ARM_OK)
		{
			return E_FAIL;
		}
		else
		{
			C_Currency = ccyres.getString ();
		}
	}


	bool isAdj = false;
	if ((C_AdjRule == "Y") || (C_AdjRule == "Adj"))  isAdj = true; 


	qCDS_ADJ AdjCDS = qCredit_Adjust20; // STDCDS
	AdjCDS = ARM_ConvAdjCalCDS(C_AdjCDS,C_result);

	if ((AdjCDS = ARM_ConvAdjCalCDS(C_AdjCDS,C_result)) == ARM_DEFAULT_ERR ) 
	{
		ERROR_MSG("Invalid Adj CDS : choose STDCDS or NONE ",pRet,ARM_ERROR_FREQ);
		return S_OK;
	}


	double Res;

	if(ICMLOCAL_Credit_AddPeriod (pAsOf,
											 C_Maturity,
											 C_Currency,
											 isAdj,
											 AdjCDS,
											 C_result) == ARM_OK)
	{
		Res = C_result.getDouble();
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
catch(...)
{
	return E_FAIL;
}
}


STDMETHODIMP ActiveXModule::ARM_Credit_SetCoupons (BSTR CdsorCdoId, 
											   BSTR CouponsId, 
											   BSTR TypesId, 
											   BSTR PartId, 
											   BSTR *pRet)
{
	// TODO: Add your implementation code here

try
{
	ARM_result C_result;

	_bstr_t bCdsorCdoId(CdsorCdoId);
	CCString l_CdsorCdoId = bCdsorCdoId;

	_bstr_t bCouponsId(CouponsId);
	CCString l_CouponsId = bCouponsId;

	_bstr_t bTypesId(TypesId);
	CCString l_TypesId = bTypesId;

	_bstr_t bPartId(PartId);
	CCString l_PartId = bPartId;

	CCString retour;

	double Res =0.;

	long m_CdsorCdoId = LocalGetNumObjectId (l_CdsorCdoId);
	long m_CouponsId = LocalGetNumObjectId (l_CouponsId);
	long m_TypesId = LocalGetNumObjectId (l_TypesId);
	long m_PartId = LocalGetNumObjectId (l_PartId);

	if(ICMLOCAL_SetCoupons(m_CdsorCdoId,m_CouponsId,m_TypesId,m_PartId,C_result) != ARM_OK)
		return E_FAIL;

	retour = "SetCorrelationMatrix OK";
	_bstr_t tmpChaine = retour;
	*pRet = SysAllocString(tmpChaine);

	return S_OK;
}
catch(...)
{
	return E_FAIL;
}
}

STDMETHODIMP ActiveXModule::ARM_Credit_SetRiskyProfile (BSTR CdsorCdoId, 
											   BSTR CouponsId, 
											   BSTR TypesId, 
											   BSTR *pRet)
{
	// TODO: Add your implementation code here

try
{
	ARM_result C_result;

	_bstr_t bCdsorCdoId(CdsorCdoId);
	CCString l_CdsorCdoId = bCdsorCdoId;

	_bstr_t bCouponsId(CouponsId);
	CCString l_CouponsId = bCouponsId;

	_bstr_t bTypesId(TypesId);
	CCString l_TypesId = bTypesId;

	CCString retour;

	double Res =0.;

	long m_CdsorCdoId = LocalGetNumObjectId (l_CdsorCdoId);
	long m_CouponsId = LocalGetNumObjectId (l_CouponsId);
	long m_TypesId = LocalGetNumObjectId (l_TypesId);

	if(ICMLOCAL_SetRiskyProfile(m_CdsorCdoId,m_TypesId,m_CouponsId,C_result) != ARM_OK)
		return E_FAIL;

	retour = "SetCorrelationMatrix OK";
	_bstr_t tmpChaine = retour;
	*pRet = SysAllocString(tmpChaine);

	return S_OK;
}
catch(...)
{
	return E_FAIL;
}
}


STDMETHODIMP ActiveXModule::ARM_Credit_DataFromLabel (BSTR	 pPricer,
												  BSTR	 pLabel,
												  double *pRet)
{
	// TODO: Add your implementation code here
try
{
	ARM_result C_result;
	
	_bstr_t bPricer(pPricer);
	CCString l_Pricer = bPricer;

	_bstr_t bLabel(pLabel);
	CCString l_Label = bLabel;

	double Res =0.;
	long Pricerid = LocalGetNumObjectId (l_Pricer);

	if(ICMLOCAL_GetDataFromLabel(Pricerid,l_Label,C_result) == ARM_OK)
	{Res = C_result.getDouble();}

	(*pRet) = Res;

	return S_OK;
}
catch(...)
{
	return E_FAIL;
}
}



STDMETHODIMP ActiveXModule::ARM_Credit_GetEqStrike(BSTR correlId,
												BSTR indexname, 
												BSTR UpOrLow, 
												VARIANT *pRetMatu,
												VARIANT *pRetStrikes)
{
	// TODO: Add your implementation code here

try
{
	ARM_result C_result;

	_bstr_t b_correlId(correlId);
	CCString C_correlId = b_correlId;

	_bstr_t b_indexname(indexname);
	CCString C_indexname = b_indexname;

	_bstr_t b_UpOrLow(UpOrLow);
	CCString C_UpOrLow = b_UpOrLow;

	int iUpOrLow = 0; 
	iUpOrLow = ARM_ConvTypeUpOrLow(C_UpOrLow,C_result);

	if ((iUpOrLow = ARM_ConvTypeUpOrLow(C_UpOrLow,C_result)) == ARM_DEFAULT_ERR ) 
	{
		ERROR_MSG("Invalid Type : choose UP or LOW or DOWN ",pRetStrikes, ARM_ERROR_FREQ);
		return S_OK;
	}

	double Res =0.;

	long Pricerid = LocalGetNumObjectId (C_correlId);

	VECTOR<double> matu;
	VECTOR<double> strikes;

	if (ICMLOCAL_GetEqStrike (Pricerid,
								C_indexname,
								iUpOrLow,
								matu,
								strikes,
								C_result)== ARM_OK)

	{
	Res = VECTORDOUBLE2VARIANT(matu, pRetMatu);
	
	Res = VECTORDOUBLE2VARIANT(strikes, pRetStrikes);
	}


	return S_OK;
}
catch(...)
{
	return E_FAIL;
}
}


STDMETHODIMP ActiveXModule::ARM_Credit_DefaultIntensity(BSTR pricerId,
													double Maturity,
													double *pRet)
{
	// TODO: Add your implementation code here
try
{
	ARM_result C_result;

	_bstr_t b_pricerId(pricerId);
	CCString C_pricerId = b_pricerId;
	long l_pricerId = LocalGetNumObjectId (C_pricerId);


	if( ICMLOCAL_DefaultIntensity (l_pricerId,
								 Maturity,	
								 0, 
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


STDMETHODIMP ActiveXModule::ARM_Credit_SetMatuLabel (BSTR pCurveId , 
											 VARIANT *pMatuLabels, 
										  	 double *pRet)
{
	// TODO: Add your implementation code here

try
{
	ARM_result C_result;

	_bstr_t b_CurveId(pCurveId);
	CCString C_CurveId = b_CurveId;

	long l_CurveId = LocalGetNumObjectId (C_CurveId);

	VECTOR<CCString> V_MatuLabels;
	long il=0;

	if(VARIANT2VECTORCCSTRING(*pMatuLabels,V_MatuLabels,il) != S_OK)
		return S_FALSE;

	if( ICMLOCAL_SetMatuLabel (l_CurveId,
								 V_MatuLabels,	
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

STDMETHODIMP ActiveXModule::ARM_Credit_SetPricerForRatesComputation(BSTR legId,
															  BSTR pricerId,
															  BSTR *pRet)
{
	// TODO: Add your implementation code here
try
{
	ARM_result C_result;

	_bstr_t b_legId(legId);
	CCString C_legId = b_legId;

	_bstr_t b_pricerId(pricerId);
	CCString C_pricerId = b_pricerId;

	CCString retour;

	long l_legId = LocalGetNumObjectId (C_legId);
	long l_pricerId = LocalGetNumObjectId (C_pricerId);


	if( ICMLOCAL_SetPricerForRatesComputation (l_legId,
											l_pricerId,	
											C_result) == ARM_OK)
	{
		retour = "SetPricerOK";
		_bstr_t tmpChaine = retour;
		*pRet = SysAllocString(tmpChaine);
	}

	return S_OK;
}
catch(...)
{
	return E_FAIL;
}
}

STDMETHODIMP ActiveXModule::ARM_Credit_SetRecovCoef (BSTR pSecId , 
											 double RecovCoef, 
										  	 double *pRet)
{
	// TODO: Add your implementation code here

try
{
	ARM_result C_result;

	_bstr_t b_SecId(pSecId);
	CCString C_SecId = b_SecId;

	long l_SecId = LocalGetNumObjectId (C_SecId);

	if( ICMLOCAL_SetRecovCoef (l_SecId,
								 RecovCoef,	
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

STDMETHODIMP ActiveXModule::ARM_Credit_SetFees(BSTR securityId,
										   BSTR RefvalueId)
{
	// TODO: Add your implementation code here
try
{
	ARM_result C_result;

	_bstr_t b_securityId(securityId);
	CCString C_securityId = b_securityId;
	long l_securityId = LocalGetNumObjectId (C_securityId);

	_bstr_t b_RefvalueId(RefvalueId);
	CCString C_RefvalueId = b_RefvalueId;
	long l_RefvalueId = LocalGetNumObjectId (C_RefvalueId);

	if( ICMLOCAL_SetFees (l_securityId,
						  l_RefvalueId,
						  C_result) == ARM_OK)
	{}

	return S_OK;
}
catch(...)
{
	return E_FAIL;
}
}

STDMETHODIMP ActiveXModule::ARM_Credit_SetInterpolationType(VARIANT *pVolCurve,
														VARIANT *pInterpolType,
														VARIANT *pRet)
{
	// TODO: Add your implementation code here

	try
		{
		
		ARM_result C_result;

		CCString l_VolCurve;
		
		if(VARIANT2CCString (*pVolCurve, l_VolCurve) != S_OK)
			return S_FALSE;
		
		double Res =0.;
		
		long VolCurveid = LocalGetNumObjectId (l_VolCurve);
		
		
		if (VolCurveid == ARM_KO)
		{
			ERROR_MSG("Invalid VolCurve",pRet,ARM_ERROR_IRCURVE);
			return S_OK;
		}
		
		CCString C_InterpolType;
		long InterpolType;
		
		if(VARIANT2CCString (*pInterpolType,C_InterpolType) != S_OK)
			return S_FALSE;

		if((InterpolType = ARM_ConvInterpMethod (C_InterpolType, C_result)) == ARM_DEFAULT_ERR)
			return S_FALSE;

		
		if(ICMLOCAL_SetInterpolationType(VolCurveid,
										 InterpolType,
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
	catch(...)
	{
		return E_FAIL;
	}
}


STDMETHODIMP ActiveXModule::ARM_Credit_RiskyPV01 (BSTR	 DefCurve_,
													  VARIANT Date1,
													  VARIANT Date2,
													  VARIANT *pRet)
{
	// TODO: Add your implementation code here
	try
	{
		std::string defCurve ; VariantTools::convert(DefCurve_,defCurve ); 
		double dateValue ; VariantTools::convert(Date1,-1.,dateValue); 
		ARM_Date date1 ; if (dateValue!=-1) VariantTools::convert(Date1,date1); 
		long defcurveId = LocalPersistent::get().getObjectId(defCurve); 
		double res ;
		if (VariantTools::isXLDate(Date2)) 
		{
			ARM_Date date2; VariantTools::convert(Date2,date2); 
			if (dateValue==-1) res= ICMLOCAL_RiskyPV01(defcurveId,NULL,date2) ;
			else res = ICMLOCAL_RiskyPV01(defcurveId,&date1,date2) ;
		}
		else 
		{
			std::string tenor ; VariantTools::convert(Date2,tenor); 
			if (dateValue==-1) res= ICMLOCAL_RiskyPV01(defcurveId,NULL,tenor) ;
			else res = ICMLOCAL_RiskyPV01(defcurveId,&date1,tenor) ;
		}
		VariantTools::convert(res,*pRet); 
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

STDMETHODIMP ActiveXModule::ARM_Credit_RiskyPV01AsSensitivity (BSTR	 DefCurve_,
													  BSTR Tenor,
													  VARIANT *pRet)
{
	// TODO: Add your implementation code here
	try
	{
		std::string defCurve ; VariantTools::convert(DefCurve_,defCurve ); 
		std::string tenor ; VariantTools::convert(Tenor,tenor); 
		long defcurveId = LocalPersistent::get().getObjectId(defCurve); 
		double res = ICMLOCAL_RiskyPV01AsSensitivity(defcurveId,tenor) ;
		VariantTools::convert(res,*pRet); 
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


STDMETHODIMP ActiveXModule::ARM_Credit_Sectorial_Correlation(DATE AsOf_,
															 BSTR structName_,
															 BSTR correlation_Type,
															VARIANT * vLabels,
															VARIANT * vector_Membership,
															double intra_Sector_Correlation,
															double inter_Sector_Correlation,
															VARIANT * vBetas,
															VARIANT * vLambdas,
															VARIANT * vBetas_Down,
															VARIANT * vLambdas_Down,
															VARIANT *pRet){

try{
	ARM_result C_result;
	long size = 0;

	qTWO_FACTORS_CORRELATION_TYPE eCorrelType ; 
	VariantTools::convert(correlation_Type,"no default",eCorrelType); 
	// _bstr_t correlType(correlation_Type);
	// CCString C_CorrelType = correlType;
	// long eCorrelType =0;
	// if((eCorrelType = ARM_Conv_Sectorial_Correlation (C_CorrelType, C_result)) == ARM_DEFAULT_ERR)
	// {
	//	ERROR_MSG("Invalid Correl Type", pRet, ARM_ERROR_CORRMATRIX);
	//	return S_OK;
	// }
	

	ARM_Date AsOf; VariantTools::convertXLDate(AsOf_,AsOf); 
	std::string structName ; VariantTools::convert(structName_,structName); 
	vector<string> labels;
	VariantTools::convert(*vLabels,labels); 
	// if(VARIANT2VECTORCCSTRING(*vLabels, labels,size) != S_OK)
	// 	return S_FALSE;

	VECTOR<int> iMemberShip;
	VECTOR<double> dMemberShip;
	if(VARIANT2VECTORDOUBLE (*vector_Membership, dMemberShip,size) != S_OK)
		return S_FALSE;
	size = dMemberShip.size();
	iMemberShip.resize(size);
	for (int i=0; i<size; i++)
		iMemberShip[i]	=	(int)	dMemberShip[i];

	VECTOR<double> betas;
	if (vBetas->pvarVal->vt == 0){
		betas.clear();
		betas.empty();
	} else {
		if(VARIANT2VECTORDOUBLE (*vBetas, betas,size) != S_OK)
			return S_FALSE;
	}

	VECTOR<double> lambdas;
	if (vLambdas->pvarVal->vt == 0 ) {
		lambdas.clear();
		lambdas.empty();
	} else {
		if(VARIANT2VECTORDOUBLE (*vLambdas, lambdas,size) != S_OK)
			return S_FALSE;
	}

	VECTOR<double> betas_down;
	if (vBetas_Down->pvarVal->vt == 0){
		betas_down.clear();
		betas_down.empty();
	} else {
		if(VARIANT2VECTORDOUBLE (*vBetas_Down, betas_down,size) != S_OK)
			return S_FALSE;
	}

	VECTOR<double> lambdas_down;
	if (vLambdas_Down->pvarVal->vt == 0 ) {
		lambdas_down.clear();
		lambdas_down.empty();
	} else {
		if(VARIANT2VECTORDOUBLE (*vLambdas_Down, lambdas_down,size) != S_OK)
			return S_FALSE;
	}

	
	long retCode =0;
	long objId =0;
	CCString Res;
	CCString prevClass;
	CCString curClass = LOCAL_PRICER_CORRMATRIX_CLASS;


	if(ICMLOCAL_SECTORIAL_CORRELATION(	AsOf,
		structName,
		eCorrelType,
									labels,
									iMemberShip,
									intra_Sector_Correlation,
									inter_Sector_Correlation,
									betas,
									lambdas,
									betas_down,
									lambdas_down,
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


STDMETHODIMP ActiveXModule::ARM_Credit_DefProbInverse (BSTR pCurveId , 
											 double dDefProba, 
										  	 VARIANT *pRet)
{

try
{
	ARM_result C_result;

	_bstr_t b_CurveId(pCurveId);
	CCString C_CurveId = b_CurveId;

	long l_CurveId = LocalGetNumObjectId (C_CurveId);

	double Res=0.;

	if( ICMLOCAL_DefProbInverse (l_CurveId,
								 dDefProba,	
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
catch(...)
{
	return E_FAIL;
}
}

STDMETHODIMP ActiveXModule::ARM_Credit_QMatrix(VARIANT *pQMatrix, 
									  VARIANT *pRet)
{
	// TODO: Add your implementation code here
	ARM_result C_result;

	try {

		ICM_QMatrix<double> lQMatrix;
		long  prevId = -1;
		VariantTools::convert(*pQMatrix, lQMatrix);
		long newId = ICMLOCAL_QMatrix(lQMatrix, prevId);
		string newLabel = LocalPersistent::get().getStringId(newId, LOCAL_CREDIT_UTIL);
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

STDMETHODIMP ActiveXModule::ARM_Credit_MarketDataMng(VARIANT *pstrVect, 
									  VARIANT *pRet)
{
	// TODO: Add your implementation code here
	ARM_result C_result;

	try {

		vector<string> vLabelStr;
		long  prevId = -1;
		VariantTools::convert(*pstrVect, vLabelStr);
		vector<long> vLabelObj;
		vLabelObj.resize(vLabelStr.size());
		for (int i = 0; i < vLabelObj.size(); i++)  
			vLabelObj[i] = LocalGetNumObjectId ((CCString)vLabelStr[i].c_str());
		long newId = ICMLOCAL_MarketDataMng(vLabelObj,prevId);
		string newLabel= LocalPersistent::get().getStringId(newId, LOCAL_MARKETDATAMANAGER_CLASS);
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


STDMETHODIMP ActiveXModule::ARM_Credit_FwdSpread (BSTR defcurveId, 
					  double Maturity1,
					  double Maturity2,
					  double FwdStartDate,
					  double FwdEndDate,
					  BSTR VolId, 
					  double *pRet)
{

// TODO: Add your implementation code here
	ARM_result C_result;

	try {
		long  prevId = -1;
		string strDefcurveId;  VariantTools::convert(defcurveId, strDefcurveId);
		string strVolId;  VariantTools::convert(VolId, strVolId);
		long newId = ICMLOCAL_FwdSpread(LocalGetNumObjectId(CCString(strDefcurveId.c_str())),
									  Maturity1,
									  Maturity2,
									  FwdStartDate,
									  FwdEndDate,
									  LocalGetNumObjectId(CCString(strVolId.c_str())),
									  C_result);
		*pRet = C_result.getDouble ();
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

STDMETHODIMP ActiveXModule::ARM_Credit_FwdSpreadAsIndex (BSTR defcurveId, 
					  double Maturity1,
					  double Maturity2,
					  double *pRet)
{

// TODO: Add your implementation code here
	ARM_result C_result;

	try {
		long  prevId = -1;
		string strDefcurveId;  VariantTools::convert(defcurveId, strDefcurveId);
		long newId = ICMLOCAL_FwdSpreadAsIndex(LocalGetNumObjectId(CCString(strDefcurveId.c_str())),
									  Maturity1,
									  Maturity2,
									  C_result);
		*pRet = C_result.getDouble ();
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

STDMETHODIMP ActiveXModule::ARM_Credit_Schedule_Info (double 	EffectiveDate,
																double 	MaturityDate,			
																BSTR 	payFrequency,
																BSTR 	ResetFreq ,
																BSTR 	DayCount,
																BSTR 	Stubrule,
																BSTR 	intRule,
																BSTR 	payCalName,
																BSTR 	PayTiming,
																BSTR 	ResetTiming,
																BSTR 	fwdRule,
																BSTR 	IncludeMaturity,
																BSTR 	adj,
																BSTR 	intStartAdj,
																BSTR 	AccDayCount,
																double 	ReferenceDate,
																double 	FirstCpnEffDate,
																BSTR 	AdjCal,
																VARIANT		*pRet)
{
	ARM_result C_result;

	try
	{
		 long  prevId = -1;
		 string strPayFreq = "";
		 int	ipayFrequency = 0;
		 string strResetFreq = "";
		 int iResetFreq =0;
		 string Basis = "";
		 int	iDayCount =0;
		 std::string strStubrule = "";
		 int	iStubrule =0;
		 string strIntRule = "";
		 int	iIntrule =0;
		 std::string strpayCalName = "";

		 string strPayTiming ="";
		 int iPayTiming =0;
		 string strResetTiming ="";
		 int iResetTiming =0 ;
		 string strfwdRule = "";
		 int ifwdRule=0;
		 string strIncludeMaturity = "N";
		 bool	IncludeMaturity = false;
		 string strAdj = "";
		 int iAdj =0;
		 string strStartAdj ="";
		 int	iIntStartAdj =0;
		 string strAccDayCount = "";
		 int iAccDayCount =0;
		 string strReferenceDate="";
		 string strAdjCal="";
		 int iAdjCal = 0;

		VariantTools::convert(payFrequency, strPayFreq);
		if((ipayFrequency = ARM_ConvFrequency(strPayFreq.c_str(), C_result)) == ARM_DEFAULT_ERR)	{ERROR_MSG("Invalid PayFreq",pRet,ARM_ERROR_FREQ); return S_OK;}			
		VariantTools::convert( ResetFreq , strResetFreq);
		if((iResetFreq = ARM_ConvFrequency (strResetFreq.c_str(), C_result)) == ARM_DEFAULT_ERR)	{ERROR_MSG("Invalid ResetFreq",pRet,ARM_ERROR_FREQ); return S_OK;}

		VariantTools::convert( DayCount,Basis);
		if((iDayCount = ARM_ConvDayCount (Basis.c_str())) == ARM_DEFAULT_ERR)	{ERROR_MSG("Invalid Basis",pRet,ARM_ERROR_FREQ); return S_OK;}

		VariantTools::convert( Stubrule,strStubrule);
		if((iStubrule = ARM_ConvStubRule (strStubrule.c_str())) == ARM_DEFAULT_ERR)	{ERROR_MSG("Invalid stubrule",pRet,ARM_ERROR_FREQ); return S_OK;}
	
		VariantTools::convert( intRule,strIntRule);
			iIntrule = ARM_ConvIntRule (strIntRule.c_str());

		VariantTools::convert( payCalName,strpayCalName);

		
		VariantTools::convert( PayTiming,strPayTiming);
			iPayTiming = ARM_ConvPayResetRule (strPayTiming.c_str());
			
		VariantTools::convert( ResetTiming,strResetTiming);
			iResetTiming = ARM_ConvPayResetRule (strResetTiming.c_str());

		VariantTools::convert( fwdRule,strfwdRule);
			if((ifwdRule =  ARM_ConvFwdRule (strfwdRule.c_str(), C_result)) == ARM_DEFAULT_ERR)	{ERROR_MSG("Invalid fwd rule index",pRet,ARM_ERROR_FREQ); return S_OK;}

		VariantTools::convert( IncludeMaturity,strIncludeMaturity);
		if (strIncludeMaturity == "Y" ) IncludeMaturity = true;

		VariantTools::convert( adj,strAdj);
			if((iAdj = ARM_ConvIntRule (strAdj.c_str())) == ARM_DEFAULT_ERR) {ERROR_MSG("Invalid adj",pRet,ARM_ERROR_FREQ); return S_OK;}
		
		VariantTools::convert( intStartAdj,strStartAdj);
			if( iIntStartAdj = ARM_ConvStartAdjRule (strStartAdj.c_str()) == ARM_DEFAULT_ERR)	{ERROR_MSG("Invalid start adj",pRet,ARM_ERROR_FREQ); return S_OK;}

		VariantTools::convert( AccDayCount,strAccDayCount);
			if((iAccDayCount = ARM_ConvDayCount (strAccDayCount.c_str())) == ARM_DEFAULT_ERR)  {ERROR_MSG("Invalid AccDayCount",pRet,ARM_ERROR_FREQ); return S_OK;}
		
		VariantTools::convert( AdjCal,strAdjCal);
			if((iAdjCal = ARM_ConvAdjCalCDS (CCString(strAdjCal.c_str()), C_result)) == ARM_DEFAULT_ERR){ERROR_MSG("Invalid AdjCal",pRet,ARM_ERROR_FREQ); return S_OK;}


		long newId = ICMLOCAL_SCHEDULE_INFO(EffectiveDate,
											MaturityDate,			
											ipayFrequency,
											iResetFreq ,
											iDayCount,
											iStubrule,
											iIntrule,
											strpayCalName,
											iPayTiming,
											iResetTiming,
											ifwdRule,
											IncludeMaturity,
											iAdj,
											iIntStartAdj,
											iAccDayCount,
											ReferenceDate,
											FirstCpnEffDate,
											iAdjCal,
											prevId);
	
		string newLabel= LocalPersistent::get().getStringId(newId, LOCAL_SCHEDULE_INFO_CLASS);
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


STDMETHODIMP ActiveXModule::ARM_Credit_PriceVector (VARIANT *pPricer , 
										  BSTR pCPTTYPE,
										  VARIANT *pRetVectorValos)
{
	// TODO: Add your implementation code here

try
{
	CCString l_Pricer;	
	if(VARIANT2CCString (*pPricer, l_Pricer) != S_OK)
		return S_FALSE;
	long Pricerid = LocalGetNumObjectId (l_Pricer);

	string strCptType="";
	VariantTools::convert(pCPTTYPE, strCptType);

	ARM_Vector output; 
	

	ICMLOCAL_PriceVector(Pricerid, strCptType, output);

	VariantTools::convert(output, *pRetVectorValos); 
	
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


STDMETHODIMP ActiveXModule::ARM_Credit_GenPrice (VARIANT *pPricer , 
										  BSTR pCPTTYPE,
										  VARIANT *pRet)
{
	// TODO: Add your implementation code here

try
{
	CCString l_Pricer;	
	if(VARIANT2CCString (*pPricer, l_Pricer) != S_OK)
		return S_FALSE;
	long Pricerid = LocalGetNumObjectId (l_Pricer);

	string strCptType="";
	VariantTools::convert(pCPTTYPE, strCptType);

	double Valo;

	ICMLOCAL_GenPrice(Pricerid, strCptType, Valo);

	VariantTools::convert(Valo, *pRet); 
	
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

STDMETHODIMP ActiveXModule::ARM_Credit_CptInterpolDefCurveOLD(BSTR pCurve,BSTR pTENOR,double pSlope,double pDate,double pInterpDate, VARIANT *pRet)
{
	// TODO: Add your implementation code here

try
{
	ARM_result C_result;

	_bstr_t b_Tenor(pTENOR);
	CCString C_Tenor = b_Tenor;

	_bstr_t b_Curve(pCurve);
	CCString C_Curve = b_Curve;

	double Res =0.;
	VECTOR<double> Output;


	if(ICMLOCAL_CptImplicitSpreadInterpolOLD (LocalGetNumObjectId (C_Curve), C_Tenor,pInterpDate, pSlope,pDate, Output,C_result) == ARM_OK)
	{
		Res = Output[0];
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
catch(...)
{
	return E_FAIL;
}
}

STDMETHODIMP ActiveXModule::ARM_Credit_GetExpectedLoss (BSTR 	pricerId,
														double 	YearTerm ,
														double		*pRet)
{
	ARM_result C_result;
	try
	{	 
		double expectedLoss = -1.;
		string strPricerId;  VariantTools::convert(pricerId, strPricerId);
		double dYearTerm = YearTerm;  
		
		if(ICMLOCAL_GetExpectedLoss (LocalGetNumObjectId (strPricerId.c_str()),
											 dYearTerm,
											C_result) == ARM_OK)
		{
			*pRet = C_result.getDouble ();
		}
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


STDMETHODIMP ActiveXModule::ARM_Credit_FunctionRegister (long address)
{
	ARM_result C_result;
	try
	{	 
	
		ICMLOCAL_Register(address,C_result);
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

STDMETHODIMP ActiveXModule::ARM_Credit_Math_BivNormale (double x, double y, double rho, double *pRet)
{
	ARM_result C_result;
	try
	{
		if(ICMLOCAL_Math_Bivariate_normale (x, y, rho, C_result) == ARM_OK)
		{
			*pRet = C_result.getDouble ();
		}
		else
			ICMTHROW(ERR_INVALID_ARGUMENT,"ARM_Credit_Math_BivNormale failed");
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

STDMETHODIMP ActiveXModule::ARM_Credit_Math_RandUniform(double seed, double *pRet)
{
	ARM_result C_result;
	try
	{
		if(ICMLOCAL_Math_random_uniform (seed, C_result) == ARM_OK)
		{
			*pRet = C_result.getDouble ();
		}
		else
			ICMTHROW(ERR_INVALID_ARGUMENT,"ARM_Credit_Math_RandUniform failed");
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
STDMETHODIMP ActiveXModule::ARM_Credit_Math_Interpol(VARIANT *X ,
													 VARIANT *Y,
													 double value,
													 double type,
													 double smooth,
													 VARIANT *Weights,
													 double modeSpline,
													 double withC1condition,
													 double leftSlope,
													 double rightSlope,
													 double *pRet)
{
	ARM_result C_result;

	try
	{
		VECTOR<double> V_X ;
		VECTOR<double> V_Y ;
		VECTOR<double> V_weights;
		VariantTools::convert(*X,V_X);
		VariantTools::convert(*Y,V_Y);
		VariantTools::convert(V_weights,*Weights);

		if(ICMLOCAL_Math_Interpol(V_X,V_Y,value,(int)type,smooth,V_weights,(int)modeSpline,(int)withC1condition,leftSlope,rightSlope,C_result) == ARM_OK)
		{
			*pRet = C_result.getDouble ();
		}
		else
			ICMTHROW(ERR_INVALID_ARGUMENT,"ARM_Credit_Math_Interpol failed");
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
STDMETHODIMP ActiveXModule::ARM_Credit_Random_Generator (BSTR RandomType, 
														 BSTR ParamId, 
														 VARIANT *pRet)
{
	//ARM_result C_result;
	try
	{
		long prevId = -1 ; 
		string strParamId; 
		string ParameterIdDefault="";
		long ParamLong = -1;
		VariantTools::convert(ParamId, "",strParamId);
		if(ParameterIdDefault != strParamId)
			ParamLong = LocalGetNumObjectId (CCString(strParamId.c_str()));
		
		string strRandomTp;  
		VariantTools::convert(RandomType, strRandomTp);
		qRAN_GEN rType;
		bool res=false;
		string list;
		ICM_EnumsCnv::cnv(strRandomTp,rType,res, list);
		if(!res){
			ICMTHROW(ERR_INVALID_ARGUMENT,"ARM_Credit_Random_Generator: type of random not defined list is" << list) ;
		}
		
		long newId = ICMLOCAL_RandomGenerator(rType,
											ParamLong, prevId);
	
		string newLabel= LocalPersistent::get().getStringId(newId, LOCAL_CREDIT_RANDOM_GENERATOR);
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
STDMETHODIMP ActiveXModule::ARM_Credit_GenerateOneRandom (BSTR RandomId,  
														 double *pRet)
{
	//ARM_result C_result;
	try
	{
		string strRandomId;  
		VariantTools::convert(RandomId, strRandomId);
		long lRandomId = LocalGetNumObjectId (strRandomId.c_str());
		
		double RandomNb;

		ICMLOCAL_GenerateOneRandom(lRandomId, RandomNb);

		*pRet = RandomNb;
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
STDMETHODIMP ActiveXModule::ARM_Credit_GenerateRandoms (BSTR RandomId,  
														int DimVector,
														 VARIANT *pRet)
{
	//ARM_result C_result;
	try
	{
		string strRandomId;  
		VariantTools::convert(RandomId, strRandomId);
		long lRandomId = LocalGetNumObjectId (strRandomId.c_str());
		
		ARM_Vector RandomVector(DimVector,0.); 

		ICMLOCAL_GenerateRandoms(lRandomId, RandomVector);

		VariantTools::convert( RandomVector, *pRet); 
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



STDMETHODIMP ActiveXModule::ARM_Credit_ResetRandom(BSTR RandomId)
{
	//ARM_result C_result;
	try
	{
		string strRandomId;  
		VariantTools::convert(RandomId, strRandomId);
		long lRandomId = LocalGetNumObjectId (strRandomId.c_str());
		
		ICMLOCAL_ResetRandom(lRandomId);

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

STDMETHODIMP ActiveXModule::ARM_Credit_PropertyList(VARIANT attrNames,
													   VARIANT attrValues,
													   VARIANT attrTypes,
													   VARIANT *pRet)
{
	//ARM_result C_result;
	try
	{
		ICM_PropertyList pl;
		VariantTools::convert(attrNames,attrValues,attrTypes,pl); 
		long objId = LocalPersistent::get().adopt(pl.Clone(),-1); 
		std::string objName = LocalPersistent::get().getStringId(objId,"LPLST"); 
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
