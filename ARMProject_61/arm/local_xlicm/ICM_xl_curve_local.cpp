#pragma warning(disable :4005 4786)

#include <libCCxll\CCxll.h>

#include "ARM_local_interglob.h"
#include "ARM_local_interface.h"

#include <ARM\libarm_local\ARM_local_glob.h>
#include <ARM\libarm_local\ARM_local_class.h>

#include <ARM\libicm_local\ICM_local_pwccurve.h>
#include <ARM\libicm_local\ICM_local_glob.h>
#include <ARM\libicm_local\ICM_local_summit.h>

#include "ExcelTools.h"
#include "fromto.h"
#include "ARMKernel\ccy\currency.h"

__declspec(dllexport) LPXLOPER WINAPI ARM_Credit_CptInterpolDefCurve  (LPXLOPER XL_inCv,
																	LPXLOPER XL_plot)
{
	static XLOPER XL_result;
	
	try {
		ARM_NOCALCIFWIZ();

		std::string defCurveId ; ExcelTools::convert(XL_inCv,defCurveId); 
		ARM_Date date; 
		std::string plot ; 
		if (ExcelTools::isXLDate(XL_plot)) ExcelTools::convert(XL_plot,date); 
		else ExcelTools::convert(XL_plot,plot); 
		double res; 
		if (!plot.empty()) 
		{
			res = ICMLOCAL_CptImplicitSpreadInterpol(
				LocalPersistent::get().getObjectId(defCurveId),plot); 
		}
		else 
		{
			res = ICMLOCAL_CptImplicitSpreadInterpol(
				LocalPersistent::get().getObjectId(defCurveId),date); 
		}
		ExcelTools::convert(res,&XL_result); 
		return (LPXLOPER)&XL_result;
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
 

__declspec(dllexport) LPXLOPER WINAPI ARM_Credit_createFlatCurve(LPXLOPER XL_inCv,LPXLOPER XL_plot)
{
	static XLOPER XL_result;
	
	try {
		ARM_NOCALCIFWIZ();

		std::string defCurveId ; ExcelTools::convert(XL_inCv,defCurveId); 
		ARM_Date date; 
		std::string plot ; 
		if (ExcelTools::isXLDate(XL_plot)) ExcelTools::convert(XL_plot,date); 
		else ExcelTools::convert(XL_plot,plot); 
		long objId ; 
		if (!plot.empty()) 
		{
			objId = ICMLOCAL_createFlatCurve(
				LocalPersistent::get().getObjectId(defCurveId),plot,ExcelCaller::get().getObjectId()); 
		}
		else 
		{
			objId = ICMLOCAL_createFlatCurve(
				LocalPersistent::get().getObjectId(defCurveId),date,ExcelCaller::get().getObjectId()); 
		}
		std::string objName = ExcelCaller::get().setObject(objId,LOCAL_ZERO_CURVE_CDS_CLASS); 
		ExcelTools::convert(objName,&XL_result); 
		return (LPXLOPER)&XL_result;
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
__declspec(dllexport) LPXLOPER WINAPI  ARM_Credit_createDefCurveFromBase(LPXLOPER XL_inCvCDS,LPXLOPER XL_inCvIndex,LPXLOPER XL_VBase)
{
	static XLOPER XL_result;
	
	try {
		ARM_NOCALCIFWIZ();

		std::string defCurveIdCDS ; ExcelTools::convert(XL_inCvCDS,defCurveIdCDS); 
		std::string defCurveIdIndex ; ExcelTools::convert(XL_inCvIndex,defCurveIdIndex); 
		ARM_Vector vBase;
		ExcelTools::convert(XL_VBase,vBase); 
		long objId  =0; 
		objId = ICMLOCAL_createDefCurveFromBase(LocalPersistent::get().getObjectId(defCurveIdCDS),
				LocalPersistent::get().getObjectId(defCurveIdIndex),vBase,ExcelCaller::get().getObjectId()); 

		std::string objName = ExcelCaller::get().setObject(objId,LOCAL_ZERO_CURVE_CDS_CLASS); 
		ExcelTools::convert(objName,&XL_result); 
		return (LPXLOPER)&XL_result;
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

__declspec(dllexport) LPXLOPER WINAPI ARM_Credit_DefCurveFromSummit (LPXLOPER XL_AsOfDate,
																	LPXLOPER XL_Issuer,
																	LPXLOPER XL_CurveName,
																	LPXLOPER XL_IRcurve,
																	LPXLOPER XL_Label)
{

	static XLOPER XL_result;
	ARM_result C_result;

	try {
    ARM_NOCALCIFWIZ();

	// C variable
	// double C_AsOfDate;
	CCString C_Issuer;
	CCString C_CurveName;
	CCString C_IRcurve;
	CCString C_Ccy;
	CCString C_label;

	long retCode=0;

	// error
	static int error;
	static char* reason = "";

	// XL_readNumCell(XL_AsOfDate,C_AsOfDate," ARM_ERR: AsOfDate of the curve, numeric expected",C_result);
	ARM_Date AsOf ; ExcelTools::convert(XL_AsOfDate,AsOf); 
	XL_readStrCell(XL_Issuer,C_Issuer," ARM_ERR: Issuer: string expected",C_result);
	XL_readStrCell(XL_CurveName,C_CurveName," ARM_ERR: curve name: string expected",C_result);
	XL_readStrCellWD(XL_IRcurve,C_IRcurve,"NONE"," ARM_ERR: : string expected",C_result);
	XL_readStrCellWD(XL_Label,C_label,"NONE"," ARM_ERR: Label : string expected",C_result);

	long objId;
	CCString prevClass;
	
	CCString curClass = LOCAL_ZERO_CURVE_CDS_CLASS;
	CCString stringId = GetLastCurCellEnvValue ();
	
	if (C_label=="NONE")
		C_label = C_Issuer;

	if(!stringId)
	{
		retCode = ICMLOCAL_GetDPFromSummit (AsOf,
											C_Issuer,
											C_CurveName,
											LocalGetNumObjectId(C_IRcurve),	
											C_label,
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
			retCode = ICMLOCAL_GetDPFromSummit (AsOf,
												C_Issuer,
												C_CurveName,
												LocalGetNumObjectId(C_IRcurve),	
												C_label,
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
			retCode = ICMLOCAL_GetDPFromSummit (AsOf,
												C_Issuer,
												C_CurveName,
												LocalGetNumObjectId(C_IRcurve),	
												C_label,
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

__declspec(dllexport) LPXLOPER WINAPI ARM_Credit_DPMktDataFromCalypso(LPXLOPER XL_AsOfDate,
																	  LPXLOPER XL_PricingEnv,
																	  LPXLOPER XL_Issuer,
																	  LPXLOPER XL_Seniority,
																	  LPXLOPER XL_Ccy,
																	  LPXLOPER XL_Parameter,
																	  LPXLOPER XL_forceCurveName,
																	  LPXLOPER XL_xmlfile)
																	  
{
	static XLOPER XL_result;
	LPXLOPER pxArray;
	
	try
	{
		ARM_NOCALCIFWIZ();

		std::string issuer ; ExcelTools::convert(XL_Issuer,"",issuer);
		std::string ccy ; ExcelTools::convert(XL_Ccy,"",ccy);
		std::string seniority; ExcelTools::convert(XL_Seniority,"",seniority);
		std::string pricingEnv; ExcelTools::convert(XL_PricingEnv,"",pricingEnv);
		std::string forceCurveName ; ExcelTools::convert(XL_forceCurveName,"",forceCurveName) ; 
		std::string xmlFileName ; ExcelTools::convert(XL_xmlfile,"",xmlFileName); 
		vector<string> vParameters; ExcelTools::convert(XL_Parameter,vector<string>(),vParameters);
		ARM_Date AsOfDate; ExcelTools::convert(XL_AsOfDate,AsOfDate); 

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

		int nbcolumns = vParameters.size() ;
		int nbrows = Matus.size() ;

		//FreeCurcellErr();
		FreeCurCellErr ();
		XL_result.xltype= xltypeMulti ;
		XL_result.val.array.columns = nbcolumns;
		XL_result.val.array.rows = nbrows;
		XL_result.val.array.lparray = pxArray = (LPXLOPER)GlobalAlloc(GMEM_ZEROINIT, nbrows * nbcolumns * sizeof(XLOPER));
		
		char * tmp_str ;
		CCString tmpCCstring ;
		

		for (int i =0 ; i<nbcolumns ; i++)
		{
			ICM_EnumsCnv::toString(adj,tmp) ;	// will throw if !ok
			if  (vParameters[i] =="MATU") {
				for (int j = 0 ; j < nbrows ; j++)
				{	
					tmp_str  =(char *) Matus[j].c_str();
					tmpCCstring = CCString(tmp_str);

					pxArray[XL_Coordonnate2Rank (j, i, nbcolumns)].xltype = xltypeStr;
					pxArray[XL_Coordonnate2Rank (j, i, nbcolumns)].val.str= XL_StrC2StrPascal(tmpCCstring);
					pxArray[XL_Coordonnate2Rank (j, i, nbcolumns)].xltype |= xlbitDLLFree;
				}
			}

			if  (vParameters[i] =="SPREAD") {
				for (int j = 0 ; j < nbrows ; j++)
				{	
					pxArray[XL_Coordonnate2Rank (j, i, nbcolumns)].xltype = xltypeNum;
					pxArray[XL_Coordonnate2Rank (j, i, nbcolumns)].val.num = Spreads[j];
				}
			}


			if  (vParameters[i] =="CCY") {
				tmp_str  =(char *) outputCcy.c_str();
				tmpCCstring = CCString(tmp_str);
				for (int j = 0 ; j < nbrows ; j++)
				{
					
					pxArray[XL_Coordonnate2Rank (j, i, nbcolumns)].xltype = xltypeStr;
					pxArray[XL_Coordonnate2Rank (j, i, nbcolumns)].val.str =XL_StrC2StrPascal(tmpCCstring);
				}
			}

			if  (vParameters[i] =="INDEX") {
				tmp_str  =(char *) indexName.c_str();
				tmpCCstring = CCString(tmp_str);
				for (int j = 0 ; j < nbrows ; j++)
				{	
					pxArray[XL_Coordonnate2Rank (j, i, nbcolumns)].xltype = xltypeStr;
					pxArray[XL_Coordonnate2Rank (j, i, nbcolumns)].val.str =XL_StrC2StrPascal(tmpCCstring);
				}
			}

			if  (vParameters[i] =="ADJCDS") {
				tmp_str  =(char *) ICM_EnumsCnv::toString(adj).c_str();
				tmpCCstring = CCString(tmp_str);
				for (int j = 0 ; j < nbrows ; j++)
				{	
					pxArray[XL_Coordonnate2Rank (j, i, nbcolumns)].xltype = xltypeStr;
					pxArray[XL_Coordonnate2Rank (j, i, nbcolumns)].val.str =XL_StrC2StrPascal(tmpCCstring);
				}
			}

			if  (vParameters[i] =="RECOVERY") {
				for (int j = 0 ; j < nbrows ; j++)
				{	
					pxArray[XL_Coordonnate2Rank (j, i, nbcolumns)].xltype = xltypeNum;;
					pxArray[XL_Coordonnate2Rank (j, i, nbcolumns)].val.num = Recovery;
				}
			}	
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


__declspec(dllexport) LPXLOPER WINAPI ARM_Credit_DefCurveFromCalypso(LPXLOPER XL_AsOfDate,
																	LPXLOPER XL_Issuer,
																	LPXLOPER XL_Ccy,
																	LPXLOPER XL_Seniority,
																	LPXLOPER XL_PricingEnv,
																	LPXLOPER XL_IRcurve,
																	LPXLOPER XL_Label,
																	LPXLOPER XL_forceCurveName,
																	LPXLOPER XL_xmlFileName)
{
	static XLOPER XL_result;
	try 
	{
		//	
		if (ExcelCaller::get().isCalledByWizard()) 
		{
			ExcelTools::convert("",&XL_result);
			return &XL_result; 
		}
		// 
		std::string issuer ; ExcelTools::convert(XL_Issuer,"",issuer) ; 
		std::string ccy; ExcelTools::convert(XL_Ccy,"",ccy) ; 
		std::string seniority; ExcelTools::convert(XL_Seniority,"",seniority) ; 
		std::string pricingEnv; ExcelTools::convert(XL_PricingEnv,"MO",pricingEnv) ; 
		std::string ircurve; ExcelTools::convert(XL_IRcurve,"",ircurve) ; 
		std::string label; ExcelTools::convert(XL_Label,"",label) ; 
		std::string forceCurveName ; ExcelTools::convert(XL_forceCurveName,"",forceCurveName) ; 
		std::string xmlFileName ; ExcelTools::convert(XL_xmlFileName,"",xmlFileName); 
		ARM_Date AsOf; ExcelTools::convert(XL_AsOfDate,AsOf); 
		// 
		long prevId = ExcelCaller::get().getObjectId();  
		long retId = ICMLOCAL_GetDPFromCalypso(AsOf,issuer,seniority,ccy,pricingEnv,forceCurveName,
												LocalGetNumObjectId(ircurve.c_str()),	
												label,xmlFileName,prevId ) ;
		std::string objectLabel =ExcelCaller::get().setObject(retId, LOCAL_ZERO_CURVE_CDS_CLASS) ;
		// ExcelCaller::get().setObject(objectLabel) ;
		ExcelTools::convert(objectLabel,&XL_result); 
		return &XL_result;		
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
	return &XL_result;
}

__declspec(dllexport) LPXLOPER WINAPI ARM_Credit_DPMktDataFromSummit (LPXLOPER XL_AsOfDate,
																	LPXLOPER XL_Issuer,
																	LPXLOPER XL_CurveName)
{

	static XLOPER XL_result;
	ARM_result C_result;
	LPXLOPER pxArray;

	try {
    ARM_NOCALCIFWIZ();

	// C variable
	double C_AsOfDate;
	CCString C_Issuer;
	CCString C_CurveName;
	CCString C_Ccy;
	CCString C_IndexName;
	long CDSAdj;

	VECTOR<CCString> Matu;
	VECTOR<double> Spread;
	VECTOR<double> Recovery;

	long retCode=0;

	// error
	static int error;
	static char* reason = "";

	XL_readNumCell(XL_AsOfDate,C_AsOfDate," ARM_ERR: AsOfDate of the curve, numeric expected",C_result);
	XL_readStrCell(XL_Issuer,C_Issuer," ARM_ERR: Issuer: string expected",C_result);
	XL_readStrCell(XL_CurveName,C_CurveName," ARM_ERR: curve name: string expected",C_result);


	CCString prevClass;
	
	CCString curClass = LOCAL_ZERO_CURVE_CDS_CLASS;
	CCString stringId = GetLastCurCellEnvValue ();
	
	retCode = ICMLOCAL_DPMktDataFromSummit (C_AsOfDate,
											C_Issuer,
											C_CurveName,
											Matu,
											Spread,
											Recovery,
											C_Ccy,
											C_IndexName,
											CDSAdj,
											C_result);


	if(retCode == ARM_OK)
	{
		int nbrows = Matu.size ();
		int nbcolumns = 5;
		
		FreeCurCellErr ();
		XL_result.xltype = xltypeMulti;
		XL_result.val.array.columns = nbcolumns;
		XL_result.val.array.rows = nbrows; 
		XL_result.val.array.lparray = pxArray = (LPXLOPER)GlobalAlloc (GMEM_ZEROINIT, nbrows * nbcolumns * sizeof (XLOPER));

		for(int i = 0; i < nbrows; i++)
		{
			pxArray[XL_Coordonnate2Rank (i, 0, nbcolumns)].xltype = xltypeStr;
			pxArray[XL_Coordonnate2Rank (i, 0, nbcolumns)].val.str = XL_StrC2StrPascal (Matu[i]);
			pxArray[XL_Coordonnate2Rank (i, 0, nbcolumns)].xltype |= xlbitDLLFree;
			pxArray[XL_Coordonnate2Rank (i, 1, nbcolumns)].xltype = xltypeNum;
			pxArray[XL_Coordonnate2Rank (i, 1, nbcolumns)].val.num = Spread[i]; 
			pxArray[XL_Coordonnate2Rank (i, 2, nbcolumns)].xltype = xltypeNum;
			pxArray[XL_Coordonnate2Rank (i, 2, nbcolumns)].val.num = Recovery[i]; 
			pxArray[XL_Coordonnate2Rank (i, 3, nbcolumns)].xltype = xltypeStr;
			pxArray[XL_Coordonnate2Rank (i, 3, nbcolumns)].val.str = XL_StrC2StrPascal (C_Ccy);
			pxArray[XL_Coordonnate2Rank (i, 4, nbcolumns)].xltype = xltypeStr;
			pxArray[XL_Coordonnate2Rank (i, 4, nbcolumns)].val.str = XL_StrC2StrPascal (C_IndexName);
		}
		
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



__declspec(dllexport) LPXLOPER WINAPI ARM_Credit_GetMktDataCorrFromSummit(LPXLOPER XL_AsOfDate,
																		 LPXLOPER XL_Issuer1,
																		 LPXLOPER XL_Issuer2,
																		 LPXLOPER XL_CurveName)
{
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	try {
	ARM_NOCALCIFWIZ();

	// C variable
	double C_AsOfDate;
	CCString C_IssuerId1;
	CCString C_IssuerId2;
	CCString C_CurveName;

	// error
	static int error;
	static char* reason = "";

	XL_readNumCell(XL_AsOfDate,C_AsOfDate," ARM_ERR: AsOfDate of the curve, numeric expected",C_result);
	XL_readStrCell(XL_Issuer1,C_IssuerId1," ARM_ERR: Issuer name 1: string expected (Defprob issuer)",C_result);
	XL_readStrCell(XL_Issuer2,C_IssuerId2," ARM_ERR: Issuer name 2: string expected (Defprob issuer)",C_result);
	XL_readStrCell(XL_CurveName,C_CurveName," ARM_ERR: curve name: string expected",C_result);

	long retCode = ICMLOCAL_CorrFromSummit(C_AsOfDate,
										  C_IssuerId1,
										  C_IssuerId2,
										  C_CurveName,
										  C_result);

	if(retCode == ARM_OK)
	{
		FreeCurCellErr ();
		XL_result.xltype = xltypeNum;
		XL_result.val.num = C_result.getDouble ();
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


__declspec(dllexport) LPXLOPER WINAPI ARM_Credit_DefCurvePWC(LPXLOPER XL_date,
															 LPXLOPER XL_Maturities,
															 LPXLOPER XL_Rates,
															 LPXLOPER XL_ircurve,
															 LPXLOPER XL_Recovery,
															 LPXLOPER XL_ccy,
															 LPXLOPER XL_label,
															 LPXLOPER XL_AdjCal,
															 LPXLOPER XL_Cal_Summit,
															 LPXLOPER XL_VolCurve,
															 LPXLOPER XL_Accrued,
															 LPXLOPER XL_CalibrationAlgo,
															 LPXLOPER XL_CalibrationData,
															 LPXLOPER XL_Lag,
															 LPXLOPER XL_Parameters,
															 LPXLOPER XL_intRule,
															 LPXLOPER XL_adjStartRule)
{
	static XLOPER XL_result;
	// ARM_result C_result;

	try {
		ARM_NOCALCIFWIZ();

		// VECTOR<double> C_Rates; 

		// CCString C_ircurve;
		// long ircurveId = 0;

		// CCString C_ccy;
		// CCString C_label;
		// CCString C_AdjCal;
		// CCString C_Parameters;

		// double lag =0;
		// double defaultLag = 2;

		// qCDS_ADJ l_AdjCal = qCredit_Adjust20; // STDCDS
		// double C_Recovery = 0.;

		// CCString C_Summit;
		// bool C_Cal_Summit = false;


		// CCString C_VolCurve;
		// long VolCurveId = 0;
		// CCString C_VolCurveDefaultValue ="NONE";

		// error
		// static int error;
		// static char* reason = "";

		ARM_Date AsOf ; ExcelTools::convert(XL_date,AsOf); 
		std::string C_ircurve; ExcelTools::convert(XL_ircurve,C_ircurve); 
		// XL_readStrCell(XL_ircurve,C_ircurve," ARM_ERR: Curve Id expected ",C_result);
		std::vector<std::string> C_matu ;
		ExcelTools::convert(XL_Maturities,C_matu ); 
		std::vector<double> C_Rates; ExcelTools::convert(XL_Rates,C_Rates); 
		// XL_readNumVector(XL_Rates,C_Rates," ARM_ERR: rates : array expected",C_result);
		double C_Recovery; ExcelTools::convert(XL_Recovery,C_Recovery); // XL_readNumCell(XL_Recovery,C_Recovery," ARM_ERR: Recovery: numeric expected",C_result);
		// XL_readStrCellWD(XL_ccy,C_ccy,"DEFAULT"," ARM_ERR: currency: string expected",C_result);
		std::string C_ccy; ExcelTools::convert(XL_ccy,C_ccy);   
		// XL_readStrCellWD(XL_label,C_label,"NONE"," ARM_ERR: Curve Label: string expected",C_result);
		std::string C_label; ExcelTools::convert(XL_label,"NONE",C_label); 
		// XL_readStrCellWD(XL_AdjCal,C_AdjCal,"STDCDS"," ARM_ERR: Adjust Cal: string expected",C_result);
		qCDS_ADJ C_AdjCal; ExcelTools::econvert(XL_AdjCal,"STDCDS",C_AdjCal); 
		std::string C_Summit; ExcelTools::convert(XL_Cal_Summit,"Y",C_Summit); 
		// XL_readStrCellWD(XL_Cal_Summit,C_Summit,"Y"," ARM_ERR: Calibration with summit: string expected",C_result);
		std::string C_VolCurve; ExcelTools::convert(XL_VolCurve,"NONE",C_VolCurve); 
		// XL_readStrCellWD(XL_VolCurve,C_VolCurve,C_VolCurveDefaultValue," ARM_ERR: Vol Curve Id expected ",C_result);
		// XL_readStrCellWD(XL_Parameters,C_Parameters,"NONE"," ARM_ERR: Parameter Id expected ",C_result);
		std::string C_Parameters; ExcelTools::convert(XL_Parameters,"NONE",C_Parameters); 
		
		std::string sIntRule; ExcelTools::convert(XL_intRule,"ADJ",sIntRule); 
		long intRule = ARM_ConvIntRule(sIntRule.c_str()) ; 
		std::string sAdjStartRule; ExcelTools::convert(XL_adjStartRule,"ADJ",sAdjStartRule); 
		long adjStartRule= ARM_ConvStartAdjRule(sAdjStartRule) ; 

		qDEFCURVE_CALIB_ALGO calibrationAlgo ; ExcelTools::econvert(XL_CalibrationAlgo,"DICHO",calibrationAlgo); 
		std::string calibrationData ; ExcelTools::convert(XL_CalibrationData,"STD",calibrationData) ;

		bool C_Cal_Summit ; 
		if (C_Summit == "Y")  C_Cal_Summit = true;
		else C_Cal_Summit =false; 

		// JLA. if(C_ccy == "DEFAULT")
		// JLA. {
		// JLA. 	ARM_result currencyres;
		// JLA. 	ARMLOCAL_ARM_GetDefaultCurrency (currencyres);
		// JLA. 	if(currencyres.getRetCode () != ARM_OK)
		// JLA. 	{
		// JLA. 		ARM_ARG_ERR();
		// JLA. 		return (LPXLOPER)&XL_result;
		// JLA. 	}
		// JLA. 	else
		// JLA. 	{
		// JLA. 		C_ccy = currencyres.getString ();
		// JLA. 	}
		// JLA. }
		
		int defaultLag= ARM_Currency(C_ccy.c_str()).GetCreditStartDateLag();

		int lag ; ExcelTools::convert(XL_Lag,defaultLag,lag);
		// XL_readNumCellWD(XL_Lag,lag,defaultLag, "ARM_ERR: lag expected ",C_result);
		
		// if((l_AdjCal = ARM_ConvAdjCalCDS (C_AdjCal, C_result)) == ARM_DEFAULT_ERR)
		// {
		// 		ARM_ARG_ERR();
		// 	return (LPXLOPER)&XL_result;
		// }

		// JLA double tmp_rate = -1;
		// JLA CCString tmp_matu ("NULL");

		// JLA JLA if(C_matu.size () != C_Rates.size())
		// JLA {
		// JLA 	C_result.setMsg ("ARM_ERR: check your maturities & spreads array");
		// JLA 	ARM_ARG_ERR();
		// JLA 	return (LPXLOPER)&XL_result;
		// JLA }

		// JLA double tmp = 0.;
		// JLA int j = 0;
 

		//	long retCode;
		// JLA CCString prevClass;
		
		// JLA CCString curClass = LOCAL_ZERO_CURVE_CDS_CLASS;
		// JLA CCString stringId = GetLastCurCellEnvValue ();


		long prevId = ExcelCaller::get().getObjectId(); 
		long objId = ICMLOCAL_ConstantDefaultCurve(AsOf, 
													C_matu, 
													C_Rates, 
													C_Recovery,
													LocalPersistent::get().getObjectId(C_ircurve),
													intRule,
													adjStartRule,
													// (long)LocalGetNumObjectId(C_ircurve),
													C_ccy,
													C_label,
													C_AdjCal,
													C_Cal_Summit,
													calibrationAlgo,
													LocalPersistent::get().getObjectId(C_VolCurve),
													// (long)LocalGetNumObjectId(C_VolCurve),
													calibrationData,
													lag,
													LocalPersistent::get().getObjectId(C_Parameters),
													// (long)LocalGetNumObjectId(C_Parameters),
													prevId);


		std::string objName = ExcelCaller::get().setObject(objId,LOCAL_ZERO_CURVE_CDS_CLASS); 
		ExcelTools::convert(objName,&XL_result); 
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


__declspec(dllexport) LPXLOPER WINAPI ARM_Credit_DefCurvePWLIN(LPXLOPER XL_date,
															 LPXLOPER XL_Maturities,
															 LPXLOPER XL_Rates,
															 LPXLOPER XL_ircurve,
															 LPXLOPER XL_Recovery,
															 LPXLOPER XL_ccy,
															 LPXLOPER XL_label,
															 LPXLOPER XL_AdjCal,
															 LPXLOPER XL_Cal_Summit)
{
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;
	try {

	ARM_NOCALCIFWIZ();

	// C variable
	// double C_date =0.;
	// VECTOR<CCString> C_matu;
	VECTOR<double> C_Rates; 

	CCString C_ircurve;
	long ircurveId = 0;

	CCString C_ccy;
	CCString C_label;
	CCString C_AdjCal;

	qCDS_ADJ l_AdjCal = qCredit_Adjust20;
	double C_Recovery = 0.;

	CCString C_Summit;
	bool C_Cal_Summit = false;

	// error
	static int error;
	static char* reason = "";

	// XL_readNumCell(XL_date,C_date," ARM_ERR: as of date: date expected",C_result);
	ARM_Date AsOf ; ExcelTools::convert(XL_date,AsOf); 
	XL_readStrCell(XL_ircurve,C_ircurve," ARM_ERR: Curve Id expected ",C_result);
	// XL_readStrVector(XL_Maturities,C_matu," ARM_ERR: maturities : array expected",DOUBLE_TYPE,C_result);
	std::vector<std::string> C_matu ;
	ExcelTools::convert(XL_Maturities,C_matu) ; 
	XL_readNumVector(XL_Rates,C_Rates," ARM_ERR: rates : array expected",C_result);
	XL_readNumCell(XL_Recovery,C_Recovery," ARM_ERR: Recovery: numeric expected",C_result);
	XL_readStrCellWD(XL_ccy,C_ccy,"DEFAULT"," ARM_ERR: currency: string expected",C_result);
	XL_readStrCellWD(XL_label,C_label,"NONE"," ARM_ERR: Curve Label: string expected",C_result);
	XL_readStrCellWD(XL_AdjCal,C_AdjCal,"STDCDS"," ARM_ERR: Curve Label: string expected",C_result);
	XL_readStrCellWD(XL_Cal_Summit,C_Summit,"Y"," ARM_ERR: Calibration with summit: string expected",C_result);

	if (C_Summit == "Y")  C_Cal_Summit = true;

	if(C_ccy == "DEFAULT")
	{
		ARM_result currencyres;
		ARMLOCAL_ARM_GetDefaultCurrency (currencyres);
		if(currencyres.getRetCode () != ARM_OK)
		{
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}
		else
		{
			C_ccy = currencyres.getString ();
		}
	}

	if((l_AdjCal = ARM_ConvAdjCalCDS (C_AdjCal, C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	double tmp_rate = -1;
	CCString tmp_matu ("NULL");

	if(C_matu.size () != C_Rates.size())
	{
		C_result.setMsg ("ARM_ERR: check your maturities & spreads array");
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	double tmp = 0.;
	int j = 0;

	/** 
	CCString last = C_matu[C_matu.size () - 1];
	if((!last) || (last[0] == ' '))
	{
		C_result.setMsg ("ARM_ERR: check your maturities and rates array");
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}
	*/ 
	long retCode;
	long objId;
	CCString prevClass;
	
	CCString curClass = LOCAL_ZERO_CURVE_CDS_CLASS;
	CCString stringId = GetLastCurCellEnvValue ();

	if(!stringId)
	{
		retCode = ICMLOCAL_LinearDefaultCurve(AsOf, 
												C_matu, 
												C_Rates, 
												C_Recovery,
												(long)LocalGetNumObjectId(C_ircurve),
												C_ccy,
												C_label,
												l_AdjCal,
												C_Cal_Summit,
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
			retCode = ICMLOCAL_LinearDefaultCurve(AsOf, 
												C_matu, 
												C_Rates, 
												C_Recovery,
												(long)LocalGetNumObjectId(C_ircurve),
												C_ccy,
												C_label,
												l_AdjCal,
												C_Cal_Summit,
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

			retCode = ICMLOCAL_LinearDefaultCurve(AsOf, 
												C_matu, 
												C_Rates, 
												C_Recovery,
												(long)LocalGetNumObjectId(C_ircurve),
												C_ccy,
												C_label,
												l_AdjCal,
												C_Cal_Summit,
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

	// clean memory
	C_matu.clear();
	C_Rates.clear();

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


__declspec(dllexport) LPXLOPER WINAPI ARM_Credit_FwdSpread (LPXLOPER XL_defcurveId, 
															LPXLOPER XL_Maturity1,
															LPXLOPER XL_Maturity2,
															LPXLOPER XL_FwdStartDate,
															LPXLOPER XL_FwdEndDate,
															LPXLOPER XL_VolId)
{

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	try {
	ARM_NOCALCIFWIZ();

	// C variable

	double C_maturity1 = 0.0;
	double C_maturity2 = 0.0;

	double C_fwdstart = 0.0;
	double C_fwdstart_default = -1.0;
	double C_fwdend = 0.0;
	double C_fwdend_default = 0.0;

	CCString C_defcurveId;
	CCString C_Adj;
	CCString C_VolId;

	CCString C_VolId_default = "NONE";
	
	// error
	static int error;
	static char* reason = "";

	XL_readStrCell(XL_defcurveId,C_defcurveId," ARM_ERR: Default curve id: object expected",C_result);
	XL_readNumCell(XL_Maturity1, C_maturity1, " ARM_ERR: maturity1 : numeric expected",C_result);
	XL_readNumCell(XL_Maturity2, C_maturity2, " ARM_ERR: maturity2 : numeric expected",C_result);
	XL_readNumCellWD(XL_FwdStartDate, C_fwdstart,C_fwdstart_default, " ARM_ERR: FwdStartDate : numeric expected",C_result);
	XL_readNumCellWD(XL_FwdEndDate, C_fwdend,C_fwdend_default, " ARM_ERR: FwdEndDate : numeric expected",C_result);
	XL_readStrCellWD(XL_VolId,C_VolId,C_VolId_default," ARM_ERR: volatility curve id: object expected",C_result);

	long retCode = ICMLOCAL_FwdSpread(LocalGetNumObjectId(C_defcurveId),
									  C_maturity1,
									  C_maturity2,
									  C_fwdstart,
									  C_fwdend,
									  LocalGetNumObjectId(C_VolId),
									  C_result);

	if(retCode == ARM_OK)
	{
		FreeCurCellErr ();
		XL_result.xltype = xltypeNum;
		XL_result.val.num = C_result.getDouble ();
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


__declspec(dllexport) LPXLOPER WINAPI ARM_Credit_DefCurveFromZCwithSummit (LPXLOPER XL_AsOfDate,
																			LPXLOPER XL_Issuer,
																			LPXLOPER XL_CurveName,
																			LPXLOPER XL_IRcurve,
																			LPXLOPER XL_Ccy,
																			LPXLOPER XL_Label)
{

	static XLOPER XL_result;
	ARM_result C_result;

	try {
    ARM_NOCALCIFWIZ();

	// C variable
	double C_AsOfDate;
	CCString C_Issuer;
	CCString C_CurveName;
	CCString C_IRcurve;
	CCString C_Ccy;
	CCString C_label;

	long retCode=0;

	// error
	static int error;
	static char* reason = "";

	XL_readNumCell(XL_AsOfDate,C_AsOfDate," ARM_ERR: AsOfDate of the curve, numeric expected",C_result);
	XL_readStrCell(XL_Issuer,C_Issuer," ARM_ERR: Issuer: string expected",C_result);
	XL_readStrCell(XL_CurveName,C_CurveName," ARM_ERR: curve name: string expected",C_result);
	XL_readStrCell(XL_IRcurve,C_IRcurve," ARM_ERR: : string expected",C_result);
	XL_readStrCellWD(XL_Ccy,C_Ccy,"DEFAULT"," ARM_ERR: currency: string expected",C_result);
	XL_readStrCellWD(XL_Label,C_label,""," ARM_ERR: Label: string expected",C_result);

	if(C_Ccy == "DEFAULT")
	{
		ARM_result currencyres;
		ARMLOCAL_ARM_GetDefaultCurrency (currencyres);
		if(currencyres.getRetCode () != ARM_OK)
		{
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}
		else
		{
			C_Ccy = currencyres.getString ();
		}
	}

	long objId;
	CCString prevClass;
	
	CCString curClass = LOCAL_ZERO_CURVE_CDS_CLASS;
	CCString stringId = GetLastCurCellEnvValue ();
	
	if(!stringId)
	{
		retCode = ICMLOCAL_GetZC_DP_FromSummit (C_Issuer,
												C_Ccy,	
												C_CurveName,
												C_AsOfDate,	
												LocalGetNumObjectId(C_IRcurve),	
												C_label,
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
			retCode = ICMLOCAL_GetZC_DP_FromSummit (C_Issuer,
													C_Ccy,	
													C_CurveName,
													C_AsOfDate,	
													LocalGetNumObjectId(C_IRcurve),	
													C_label,
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
			retCode = ICMLOCAL_GetZC_DP_FromSummit (C_Issuer,
													C_Ccy,	
													C_CurveName,
													C_AsOfDate,	
													LocalGetNumObjectId(C_IRcurve),	
													C_label,
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


__declspec(dllexport) LPXLOPER WINAPI ARM_Credit_FwdSpreadAsIndex (LPXLOPER XL_DefCurve,
																   LPXLOPER XL_Maturity1,
																   LPXLOPER XL_Maturity2)
{
	static XLOPER XL_result;
	ARM_result C_result;
	try {
		double C_maturity1 = 0.0;
		double C_maturity2 = 0.0;
		string DefCurveId ="";
	
		ExcelTools::convert(XL_DefCurve, DefCurveId); 
		ExcelTools::convert(XL_Maturity1, C_maturity1);
		ExcelTools::convert(XL_Maturity2, C_maturity2);
	
		long retCode = ICMLOCAL_FwdSpreadAsIndex(LocalGetNumObjectId(CCString(DefCurveId.c_str())),
											C_maturity1,
											C_maturity2,
											C_result);

		if(retCode == ARM_OK)
		{
			FreeCurCellErr ();
			XL_result.xltype = xltypeNum;
			XL_result.val.num = C_result.getDouble ();
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

/** 
__declspec(dllexport) LPXLOPER WINAPI ARM_Credit_DefCurveFromIntensityOrSurvivalProbas_YF(LPXLOPER XL_date,
																						  LPXLOPER XL_matuRates,
																						  LPXLOPER XL_Inputs,
																						  LPXLOPER XL_Type,																					   
																						  LPXLOPER XL_ircurve,
																						  LPXLOPER XL_Recovery,
																						  LPXLOPER XL_ccy,
																						  LPXLOPER XL_label,
																						  LPXLOPER XL_VolCurve,
																						  LPXLOPER XL_CalibrationMethod,
																						  LPXLOPER XL_Lag,
																						  LPXLOPER XL_intRule,
																						  LPXLOPER XL_adjStartRule
																						  )
{
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	// ARM_result C_result;

	try {
		ARM_NOCALCIFWIZ();

		// C variable
		// double C_date =0.;
		// VECTOR<double> C_matu;
		// VECTOR<double> C_Inputs; 
		
		// CCString C_ircurve;
		// long ircurveId = 0;

		// CCString C_ccy;
		// CCString C_label;
		// CCString C_AdjCal;

		// double C_Recovery = 0.;

		// CCString C_VolCurve;
		// long VolCurveId = 0;
		// CCString C_VolCurveDefaultValue ="NONE";

		// double C_Type;
		// double C_Type_default = -1. ;

		// error
		// static int error;
		// static char* reason = "";

		// double lag = 0;
		// double defaultLag = 2;

		// XL_readNumCell(XL_date,C_date," ARM_ERR: as of date: date expected",C_result);
		ARM_Date AsOf ; ExcelTools::convert(XL_date,AsOf); 
		std::string C_ircurve; ExcelTools::convert(XL_ircurve,C_ircurve); // XL_readStrCell(XL_ircurve,C_ircurve," ARM_ERR: Curve Id expected ",C_result);
		std::vector<double> C_matu ; ExcelTools::convert(XL_matuRates,C_matu ); // XL_readNumVector(XL_matuRates, C_matu ," ARM_ERR: Year Fractions : numeric vector expected",C_result);
		std::vector<double> C_Inputs; ExcelTools::convert(XL_Inputs,C_Inputs); // XL_readNumVector(XL_Inputs,C_Inputs," ARM_ERR: Inputs : array expected",C_result);
		double C_Type; ExcelTools::convert(XL_Type,-1.,C_Type); // XL_readNumCellWD(XL_Type,C_Type,C_Type_default," ARM_ERR: 1 or -1 expected",C_result);
		double C_Recovery; ExcelTools::convert(XL_Recovery,C_Recovery); // XL_readNumCell(XL_Recovery,C_Recovery," ARM_ERR: Recovery: numeric expected",C_result);
		std::string C_ccy ; ExcelTools::convert(XL_ccy,C_ccy);//  XL_readStrCellWD(XL_ccy,C_ccy,"DEFAULT"," ARM_ERR: currency: string expected",C_result);
		std::string C_label ; ExcelTools::convert(XL_label,"NONE",C_label) ;// XL_readStrCellWD(XL_label,C_label,"NONE"," ARM_ERR: Curve Label: string expected",C_result);
		std::string C_VolCurve; ExcelTools::convert(XL_VolCurve,"NONE",C_VolCurve); // XL_readStrCellWD(XL_VolCurve,C_VolCurve,C_VolCurveDefaultValue," ARM_ERR: Vol Curve Id expected ",C_result);
		
		std::string calibrationMethod ; ExcelTools::convert(XL_CalibrationMethod,"STD",calibrationMethod); 
		//JLA if (C_matu.size() != C_Inputs.size())
		//JLA {
		//JLA 	C_result.setMsg("ERROR: nb of Inputs differ from nb of maturates");
		//JLA 	ARM_ARG_ERR();
		//JLA 	return (LPXLOPER)&XL_result;
		//JLA }

		// JLA if(C_ccy == "DEFAULT")
		// JLA {
		// JLA 	ARM_result currencyres;
		// JLA 	ARMLOCAL_ARM_GetDefaultCurrency (currencyres);
		// JLA 	if(currencyres.getRetCode () != ARM_OK)
		// JLA 	{
		// JLA 		ARM_ARG_ERR();
		// JLA 		return (LPXLOPER)&XL_result;
		// JLA 	}
		// JLA 	else
		// JLA 	{
		// JLA 		C_ccy = currencyres.getString ();
		// JLA 	}
		// JLA }
		int defaultLag = ARM_Currency(C_ccy.c_str()).GetCreditStartDateLag();
		int lag ; ExcelTools::convert(XL_Lag,defaultLag,lag); // XL_readNumCellWD(XL_Lag,lag,defaultLag) ; // ," ARM_ERR: Lag expected",C_result);
		std::string sIntRule ; ExcelTools::convert(XL_intRule,sIntRule); 
		std::string sAdjStartRule; ExcelTools::convert(XL_adjStartRule,sAdjStartRule); 
		// CCString prevClass;
		
		// CCString curClass = LOCAL_ZERO_CURVE_CDS_CLASS;
		// CCString stringId = GetLastCurCellEnvValue ();

		long prevId = ExcelCaller::get().getObjectId(); 
		long objId = ICMLOCAL_CstDefCurve(AsOf,
			    						   C_matu, 
										   C_Inputs, 
										   C_Type,
										   C_Recovery,
										   LocalPersistent::get().getObjectId(C_ircurve),
										   ARM_ConvIntRule(sIntRule.c_str()),
										   ARM_ConvStartAdjRule(sAdjStartRule),
										   // (long)LocalGetNumObjectId(C_ircurve),
										   C_ccy,
										   C_label,
										   // (long)LocalGetNumObjectId(C_VolCurve),
											LocalPersistent::get().getObjectId(C_VolCurve),
										   calibrationMethod, 
											lag,
											prevId);

		std::string objName = ExcelCaller::get().setObject(objId,LOCAL_ZERO_CURVE_CDS_CLASS); 
		ExcelTools::convert(objName,&XL_result); 

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
**/ 

__declspec(dllexport) LPXLOPER WINAPI ARM_Credit_DefCurveFromIntensityOrSurvivalProbas_Dates(LPXLOPER XL_date,
																							 LPXLOPER XL_matuRates,
																							 LPXLOPER XL_Inputs,
																							 LPXLOPER XL_Type,																					   
																							 LPXLOPER XL_ircurve,
																							 LPXLOPER XL_Recovery,
																							 LPXLOPER XL_ccy,
																							 LPXLOPER XL_label,
																							 LPXLOPER XL_VolCurve,
																							 LPXLOPER XL_CalibrationMethod,
																							 LPXLOPER XL_Lag,
																							 LPXLOPER XL_Parameters)
{
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	try {
	ARM_NOCALCIFWIZ();

	// C variable
	// double C_date =0.;
	VECTOR<double> C_matu;
	VECTOR<double> C_Inputs; 
	
	CCString C_ircurve;
	CCString C_paramid;
	long ircurveId = 0;

	CCString C_ccy;
	CCString C_label;
	CCString C_AdjCal;

	double C_Recovery = 0.;

	CCString C_VolCurve;
	long VolCurveId = 0;
	CCString C_VolCurveDefaultValue ="NONE";

	double C_Type;
	double C_Type_default = -1. ;

	// error
	static int error;
	static char* reason = "";

	double lag = 0;
	double defaultLag = 2;

	// XL_readNumCell(XL_date,C_date," ARM_ERR: as of date: date expected",C_result);
	ARM_Date AsOf; ExcelTools::convert(XL_date,AsOf); 
	XL_readStrCell(XL_ircurve,C_ircurve," ARM_ERR: Curve Id expected ",C_result);
	XL_readStrCellWD(XL_Parameters,C_paramid,"NONE"," ARM_ERR: Param Id expected ",C_result);
	XL_readNumVector(XL_matuRates, C_matu ," ARM_ERR: Dates : Dates expected",C_result);
	XL_readNumVector(XL_Inputs,C_Inputs," ARM_ERR: Inputs : array expected",C_result);
	XL_readNumCellWD(XL_Type,C_Type,C_Type_default," ARM_ERR: 1 or -1 expected",C_result);
	XL_readNumCell(XL_Recovery,C_Recovery," ARM_ERR: Recovery: numeric expected",C_result);
	XL_readStrCellWD(XL_ccy,C_ccy,"DEFAULT"," ARM_ERR: currency: string expected",C_result);
	XL_readStrCellWD(XL_label,C_label,"NONE"," ARM_ERR: Curve Label: string expected",C_result);
	XL_readStrCellWD(XL_VolCurve,C_VolCurve,C_VolCurveDefaultValue," ARM_ERR: Vol Curve Id expected ",C_result);
	

	std::string calibrationMethod ; ExcelTools::convert(XL_CalibrationMethod,"STD",calibrationMethod); 
	if (C_matu.size() != C_Inputs.size())
	{
		C_result.setMsg("ERROR: nb of Inputs differ from nb of maturates");
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	if(C_ccy == "DEFAULT")
	{
		ARM_result currencyres;
		ARMLOCAL_ARM_GetDefaultCurrency (currencyres);
		if(currencyres.getRetCode () != ARM_OK)
		{
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}
		else
		{
			C_ccy = currencyres.getString ();
		}
	}
	defaultLag = ARM_Currency(C_ccy.c_str()).GetCreditStartDateLag();
	XL_readNumCellWD(XL_Lag,lag,defaultLag," ARM_ERR: Lag expected",C_result);
	CCString prevClass;
	
	CCString curClass = LOCAL_ZERO_CURVE_CDS_CLASS;
	CCString stringId = GetLastCurCellEnvValue ();

	long prevId = ExcelCaller::get().getObjectId(); 
	long objId = ICMLOCAL_CstDefCurve_Dates(AsOf,
											 C_matu, 
											 C_Inputs, 
											 C_Type,
											 C_Recovery,
											 (long)LocalGetNumObjectId(C_ircurve),
											 C_ccy,
											 C_label,
											 (long)LocalGetNumObjectId(C_VolCurve),
											 calibrationMethod, (int) lag,
											 (long)LocalGetNumObjectId(C_paramid),
											 prevId);
	std::string objName = ExcelCaller::get().setObject(objId,LOCAL_ZERO_CURVE_CDS_CLASS); 
	ExcelTools::convert(objName,&XL_result); 


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




__declspec(dllexport) LPXLOPER WINAPI ARM_Credit_InputDefaultCurve(LPXLOPER XL_date,
																	LPXLOPER XL_matuRates,
																	LPXLOPER XL_Inputs,																				   
																	LPXLOPER XL_ircurve,
																	LPXLOPER XL_Recovery,
																	LPXLOPER XL_ccy,
																	LPXLOPER XL_label
																	// LPXLOPER XL_InterpolTp
																	)
{

	static XLOPER XL_result;
	try
	{
		ARM_NOCALCIFWIZ();

		ARM_Date AsOf ; ExcelTools::convert(XL_date,AsOf); 
		std::vector<ARM_Date> matuRates ;  ExcelTools::convert(XL_matuRates,matuRates); 
		ARM_Vector inputs ;  ExcelTools::convert(XL_Inputs ,inputs ); 
		std::string ircurve ;  ExcelTools::convert(XL_ircurve,ircurve); 
		long ircurveid = LocalPersistent::get().getObjectId(ircurve); 
		double recovery ; ExcelTools::convert(XL_Recovery,recovery); 
		std::string ccy ; ExcelTools::convert(XL_ccy,ccy);
		std::string label; ExcelTools::convert(XL_label,"NONE",label);
		long prevId = ExcelCaller::get().getObjectId(); 
		long objId = ICMLOCAL_InputDefCurve_Dates(AsOf,matuRates,inputs,recovery,ircurveid,ccy,label,prevId); 
		std::string objName = ExcelCaller::get().setObject(objId,LOCAL_ZERO_CURVE_CDS_CLASS); 
		ExcelTools::convert(objName,&XL_result); 
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
	return &XL_result;
/** 
	// return
	static XLOPER XL_result;
	try 
	{

		ARM_result C_result;

		ARM_NOCALCIFWIZ();

		// C variable
		double C_date =0.;
		VECTOR<double> C_matu;
		VECTOR<double> C_Inputs; 
		
		CCString C_ircurve;
		long ircurveId = 0;

		CCString C_ccy;
		CCString C_label;

		// CCString C_InterpolTp;
		

		double C_Recovery = 0.;

		// error
		static int error;
		static char* reason = "";

		XL_readNumCell(XL_date,C_date," ARM_ERR: as of date: date expected",C_result);
		XL_readStrCell(XL_ircurve,C_ircurve," ARM_ERR: Curve Id expected ",C_result);
		XL_readNumVector(XL_matuRates, C_matu ," ARM_ERR: Dates : Dates expected",C_result);
		XL_readNumVector(XL_Inputs,C_Inputs," ARM_ERR: Inputs : array expected",C_result);
	//	XL_readNumCellWD(XL_Type,C_Type,C_Type_default," ARM_ERR: 1 or -1 expected",C_result);
		XL_readNumCell(XL_Recovery,C_Recovery," ARM_ERR: Recovery: numeric expected",C_result);
		XL_readStrCellWD(XL_ccy,C_ccy,"DEFAULT"," ARM_ERR: currency: string expected",C_result);
		XL_readStrCellWD(XL_label,C_label,"NONE"," ARM_ERR: Curve Label: string expected",C_result);
		// XL_readStrCellWD(XL_InterpolTp,C_InterpolTp,"LINEAR"," ARM_ERR:  string expected",C_result);
		
		qINTERPOL_TYPE l_InterpolTp ;
		ExcelTools::econvert(XL_InterpolTp,"LINEAR",l_InterpolTp ); 
		

		//	XL_readStrCellWD(XL_VolCurve,C_VolCurve,C_VolCurveDefaultValue," ARM_ERR: Vol Curve Id expected ",C_result);

		if (C_matu.size() != C_Inputs.size())
		{
			C_result.setMsg("ERROR: nb of Inputs differ from nb of maturities");
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}

		if(C_ccy == "DEFAULT")
		{
			ARM_result currencyres;
			ARMLOCAL_ARM_GetDefaultCurrency (currencyres);
			if(currencyres.getRetCode () != ARM_OK)
			{
				ARM_ARG_ERR();
				return (LPXLOPER)&XL_result;
			}
			else
			{
				C_ccy = currencyres.getString ();
			}
		}

		// l_InterpolTp = ARM_ConvInterpolType (C_InterpolTp); 
		// if((l_InterpolTp = ARM_ConvInterpolType (C_InterpolTp, C_result)) == ARM_DEFAULT_ERR)
		// {
		// 	ARM_ARG_ERR();
		// 	return (LPXLOPER)&XL_result;
		// }

		long retCode;
		long objId;
		CCString prevClass;
		
		CCString curClass = LOCAL_ZERO_CURVE_CDS_CLASS;
		CCString stringId = GetLastCurCellEnvValue ();

		if(!stringId)
		{
			retCode = ICMLOCAL_InputDefCurve_Dates(C_date,
												 C_matu, 
												 C_Inputs, 
												 C_Recovery,
												 (long)LocalGetNumObjectId(C_ircurve),
												 C_ccy,
												 C_label,
												 l_InterpolTp,
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
			retCode = ICMLOCAL_InputDefCurve_Dates(C_date,
												 C_matu, 
												 C_Inputs, 
												 C_Recovery,
												 (long)LocalGetNumObjectId(C_ircurve),
												 C_ccy,
												 C_label,
												 l_InterpolTp,
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

			retCode = ICMLOCAL_InputDefCurve_Dates(C_date,
												 C_matu, 
												 C_Inputs, 
												 C_Recovery,
												 (long)LocalGetNumObjectId(C_ircurve),
												 C_ccy,
												 C_label,
												 l_InterpolTp,
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

		// clean memory
		C_matu.clear();
		C_Inputs.clear();

		return (LPXLOPER)&XL_result;
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
	return &XL_result;
	**/ 
}




/** 
__declspec(dllexport) LPXLOPER WINAPI ARM_Credit_DefCurveIndex(LPXLOPER XL_date,
															   LPXLOPER XL_Maturities,
															   LPXLOPER XL_Rates,
															   LPXLOPER XL_RefRates,
															   LPXLOPER XL_ircurve,
															   LPXLOPER XL_Recovery,
															   LPXLOPER XL_ccy,
															   LPXLOPER XL_label,
															   LPXLOPER XL_AdjCal,
															   LPXLOPER XL_Cal_Summit,
															   LPXLOPER XL_VolCurve,
															   LPXLOPER XL_Accrued,
															   LPXLOPER XL_CalibrationAlgo,
															   LPXLOPER XL_CalibrationData,
															   LPXLOPER XL_Lag)
{
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
try {

	ARM_result C_result;

	ARM_NOCALCIFWIZ();

	VECTOR<double> C_Rates; 
	VECTOR<double> C_RefRates; 

	CCString C_ircurve;
	long ircurveId = 0;

	CCString C_ccy;
	CCString C_label;
	CCString C_AdjCal;

	qCDS_ADJ l_AdjCal = qCredit_Adjust20;
	double C_Recovery = 0.;

	CCString C_Summit;
	bool C_Cal_Summit = false;

	// CCString C_Brent_Solver;
	// bool C_Cal_Brent_Solver = false;

	CCString C_VolCurve;
	long VolCurveId = 0;
	CCString C_VolCurveDefaultValue ="NONE";

	// error
	static int error;
	static char* reason = "";

	double lag = 0;
	double defaultLag = 2;

	ARM_Date AsOf; ExcelTools::convert(XL_date,AsOf); 
	XL_readStrCell(XL_ircurve,C_ircurve," ARM_ERR: Curve Id expected ",C_result);
	std::vector<std::string> C_matu; 
	ExcelTools::convert(XL_Maturities,C_matu); 
	XL_readNumVector(XL_Rates,C_Rates," ARM_ERR: rates : array expected",C_result);
	XL_readNumVector(XL_RefRates,C_RefRates," ARM_ERR: rates : array expected",C_result);
	XL_readNumCell(XL_Recovery,C_Recovery," ARM_ERR: Recovery: numeric expected",C_result);
	XL_readStrCellWD(XL_ccy,C_ccy,"DEFAULT"," ARM_ERR: currency: string expected",C_result);
	XL_readStrCellWD(XL_label,C_label,"NONE"," ARM_ERR: Curve Label: string expected",C_result);
	XL_readStrCellWD(XL_AdjCal,C_AdjCal,"STDCDS"," ARM_ERR: Adjust Cal: string expected",C_result);
	XL_readStrCellWD(XL_Cal_Summit,C_Summit,"Y"," ARM_ERR: Calibration with summit: string expected",C_result);
	// XL_readStrCellWD(XL_Brent_Solver, C_Brent_Solver, "N", " ARM_ERR: Calibration with Brent Solver: string expected",C_result);
	XL_readStrCellWD(XL_VolCurve,C_VolCurve,C_VolCurveDefaultValue," ARM_ERR: Vol Curve Id expected ",C_result);
	std::string calibrationData; ExcelTools::convert(XL_CalibrationData,"STD",calibrationData); 
	qDEFCURVE_CALIB_ALGO calibrationAlgo ;ExcelTools::econvert(XL_CalibrationAlgo,"DICHO",calibrationAlgo); 
	if (C_Summit == "Y")  C_Cal_Summit = true;
	// if (C_Brent_Solver == "Y")  C_Cal_Brent_Solver = true;

	if(C_ccy == "DEFAULT")
	{
		ARM_result currencyres;
		ARMLOCAL_ARM_GetDefaultCurrency (currencyres);
		if(currencyres.getRetCode () != ARM_OK)
		{
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}
		else
		{
			C_ccy = currencyres.getString ();
		}
	}

	defaultLag = ARM_Currency(C_ccy.c_str()).GetCreditStartDateLag();
	XL_readNumCellWD(XL_Lag,lag,defaultLag," ARM_ERR: Lag expected",C_result);
	if((l_AdjCal = ARM_ConvAdjCalCDS (C_AdjCal, C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	double tmp_rate = -1;
	CCString tmp_matu ("NULL");

	if(C_matu.size () != C_Rates.size())
	{
		C_result.setMsg ("ARM_ERR: check your maturities & spreads array");
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	double tmp = 0.;
	int j = 0;

	
	
	CCString prevClass;
	
	CCString curClass = LOCAL_ZERO_CURVE_CDS_CLASS;
	CCString stringId = GetLastCurCellEnvValue ();

	long prevId = ExcelCaller::get().getObjectId(); 
	long objId = ICMLOCAL_DefaultCurveIndexPWC(AsOf, 
												C_matu, 
												C_Rates, 
												C_RefRates,
												C_Recovery,
												(long)LocalGetNumObjectId(C_ircurve),
												C_ccy,
												C_label,
												l_AdjCal,
												C_Cal_Summit,
												calibrationAlgo,
												(long)LocalGetNumObjectId(C_VolCurve),
												calibrationData, (int) lag,
												prevId);

	std::string objName = ExcelCaller::get().setObject(objId,LOCAL_ZERO_CURVE_CDS_CLASS); 
	ExcelTools::convert(objName,&XL_result); 


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
	return &XL_result;
}

**/ 

__declspec(dllexport) LPXLOPER WINAPI ARM_Credit_ImpliedLossTree(LPXLOPER XL_asofdate,
																 LPXLOPER XL_defcurve,
															     LPXLOPER XL_VolCurve,
																 LPXLOPER XL_correltype,
																 LPXLOPER XL_adjusted,
																 LPXLOPER XL_step)
{
	// return
	static XLOPER XL_result;
	ARM_result C_result;

	try {

	ARM_NOCALCIFWIZ();

	// C variable
	double C_date =0.;
	CCString C_defcurve;
	long defcurveId = 0;

	CCString C_VolCurve;
	long VolCurveId = 0;

	CCString C_Correlationtype;
	long l_Correlationtype = 0;

	// error
	static int error;
	static char* reason = "";

	CCString c_Adjust;
	CCString c_Adjust_default="YES";

	double c_step =0.;
	double c_step_default =30.;


	XL_readNumCell(XL_asofdate,C_date," ARM_ERR: as of date: date expected",C_result);
	XL_readStrCell(XL_defcurve,C_defcurve," ARM_ERR: DefCurve Id expected ",C_result);
	XL_readStrCell(XL_VolCurve,C_VolCurve," ARM_ERR: Vol Curve Id expected ",C_result);
	XL_readStrCellWD(XL_correltype,C_Correlationtype,"BASE"," ARM_ERR: correltype expected ",C_result);
	XL_readStrCellWD(XL_adjusted,c_Adjust,c_Adjust_default," ARM_ERR: correltype expected ",C_result);
	XL_readNumCellWD(XL_step,c_step,c_step_default," ARM_ERR: step: numeric expected",C_result);


	bool adjusted = false;

	if (c_Adjust == c_Adjust_default) 
		adjusted = true;

	long retCode;
	long objId;
	CCString prevClass;
	
	CCString curClass = LOCAL_ZERO_CURVE_CDS_CLASS;
	CCString stringId = GetLastCurCellEnvValue ();

	if((l_Correlationtype = ARM_ConvCalibType (C_Correlationtype, C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	if(!stringId)
	{
		retCode = ICMLOCAL_ImpliedLossTree(C_date, 
										   (long)LocalGetNumObjectId(C_defcurve),
										   (long)LocalGetNumObjectId(C_VolCurve),
										   l_Correlationtype,
										   adjusted,
										   (int)c_step,
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
			retCode = ICMLOCAL_ImpliedLossTree(C_date, 
										   (long)LocalGetNumObjectId(C_defcurve),
										   (long)LocalGetNumObjectId(C_VolCurve),
										   l_Correlationtype,
										   adjusted,
										   c_step,
										   C_result
										   ,objId);

			if(retCode == ARM_OK)
			{			
				LocalSetCurCellEnvValue (curClass, objId); 

				stringId = LocalMakeObjectId (objId, curClass);
			}
		}
		else
		{
			FreeCurCellContent ();

			retCode = ICMLOCAL_ImpliedLossTree(C_date, 
										   (long)LocalGetNumObjectId(C_defcurve),
										   (long)LocalGetNumObjectId(C_VolCurve),
										   l_Correlationtype,
										   adjusted,
										   c_step,
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


__declspec(dllexport) LPXLOPER WINAPI ARM_Credit_ModifyLossTree(LPXLOPER XL_TreeId,
																LPXLOPER XL_IndexNo,
																LPXLOPER XL_VectordatesYF,
																LPXLOPER XL_VectorlossesYF,
																LPXLOPER XL_transproba,
																LPXLOPER XL_stateloss)
{
//    ARM_BEGIN();

    // return
    static XLOPER XL_result;
    ARM_result C_result;

	try {

	ARM_NOCALCIFWIZ();

    // C variable
	int i = 0,j =0;
	double nbcf =0.;
	long nbrows=0,nbcolumns=0;

	VECTOR<double> C_yterms;
	VECTOR<double> C_nolosses;
	VECTOR<double> C_transproba;
	VECTOR<double> C_stateloss;
	VECTOR<double> vectdefault;

    // error
    static int error;
    static char* reason = "";

	CCString C_TreeId;
	double C_IndexNo;

	XL_readStrCell(XL_TreeId,C_TreeId," ARM_ERR: tree : object expected",C_result);
	XL_readNumCell(XL_IndexNo,C_IndexNo," ARM_ERR: index : numeric expected",C_result);
	XL_readNumVector(XL_VectordatesYF,C_yterms," ARM_ERR: year terms : vector expected",C_result);
	XL_readNumVector(XL_VectorlossesYF,C_nolosses," ARM_ERR: losses : vector expected",C_result);

    XL_readNumVectorWD(XL_transproba,C_transproba,vectdefault," ARM_ERR: trans proba : object matrix expected",C_result);
	XL_readNumVectorWD(XL_stateloss,C_stateloss,vectdefault," ARM_ERR: state losses : object matrix expected",C_result);

	double size_transproba = C_transproba.size();
	double size_stateloss = C_stateloss.size();

	if ((size_transproba == 0)&&(size_stateloss == 0))
	{
		C_result.setMsg ("ARM_ERR: check your transition probability Matrix");
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	long retCode;
    CCString prevClass;
    
    CCString curClass = LOCAL_CREDIT_CASHFLOWS_CLASS;
    CCString stringId = "OK";

	retCode = ICMLOCAL_ModifyLossTree (LocalGetNumObjectId(C_TreeId),
									   (long)C_IndexNo,
										C_yterms,
										C_nolosses,
										C_transproba,
										C_stateloss);

    if ( retCode == ARM_OK )
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


__declspec(dllexport) LPXLOPER WINAPI ARM_Credit_QuickDefaultCurve (LPXLOPER XL_Spread,
																	LPXLOPER XL_Recovery,
																	LPXLOPER XL_Label)
{

	static XLOPER XL_result;
	ARM_result C_result;

	try {
    ARM_NOCALCIFWIZ();

	// C variable
	double C_Spread;
	double C_Recovery;
	CCString C_label;

	long retCode=0;

	// error
	static int error;
	static char* reason = "";

	XL_readNumCell(XL_Spread,C_Spread," ARM_ERR: spread, numeric expected",C_result);
	XL_readNumCell(XL_Recovery,C_Recovery," ARM_ERR: recovery: numeric expected",C_result);
	XL_readStrCell(XL_Label,C_label," ARM_ERR: Label : string expected",C_result);

	long objId;
	CCString prevClass;
	
	CCString curClass = LOCAL_ZERO_CURVE_CDS_CLASS;
	CCString stringId = GetLastCurCellEnvValue ();
	
	if(!stringId)
	{
		retCode = ICMLOCAL_QuickDefaultCurve (C_Spread,
											C_Recovery,
											C_label,
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
			retCode = ICMLOCAL_QuickDefaultCurve (C_Spread,
											C_Recovery,
											C_label,
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
			retCode = ICMLOCAL_QuickDefaultCurve (C_Spread,
											C_Recovery,
											C_label,
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

__declspec(dllexport) LPXLOPER WINAPI ARM_Credit_FixingCurve(LPXLOPER XL_VDates,
															 LPXLOPER XL_Vvalue,
															 LPXLOPER XL_AsOfDate,																				   
															 LPXLOPER XL_IndexName,
															 LPXLOPER XL_IndexID)
{
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	try 
	{
		long prevId = ExcelCaller::get().getObjectId();
		ARM_result C_result;
		ARM_NOCALCIFWIZ();
		// C variable
		double C_date =0.;
		vector<ARM_Date> v_Dates;
		vector<double> v_Values; 
		ARM_Date AsOf ;
		ARM_Date defaultDate;
		defaultDate.Today();
		string  indexName = "";
		string  indexNameID = "";
		long indexId = 0;

		// error
		static int error;
		static char* reason = "";

		ExcelTools::convert(XL_VDates, v_Dates);
		ExcelTools::convert(XL_Vvalue, v_Values);
		ExcelTools::convert(XL_AsOfDate,defaultDate, AsOf); 
		ExcelTools::convert(XL_IndexName, "", indexName);
		if (indexName.empty()){
			ExcelTools::convert(XL_IndexID, "", indexNameID);
			if (indexNameID.empty()){
				C_result.setMsg("ERROR: IndexName or IndexId must be set");
				ARM_ARG_ERR();
				return (LPXLOPER)&XL_result;
			}
		}
		
		long retCode = ICMLOCAL_Fixing_Curve(v_Dates,v_Values, 
										AsOf, indexName, 
										(long)LocalGetNumObjectId(CCString(indexNameID.c_str())),
										C_result, prevId);
		if (retCode == ARM_KO){
			ARM_ERR();
		} else {
		
			string objectLabel = ExcelCaller::get().setObject(retCode,LOCAL_FIXING_CURVE_CLASS);
			//ExcelCaller::get().setObject(objectLabel);
			ExcelTools::convert(objectLabel, &XL_result);
		}
		return (LPXLOPER)&XL_result;
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

__declspec(dllexport) LPXLOPER WINAPI ARM_Credit_CreateBasketCorrelMkDataFromCalypso(LPXLOPER XL_PricingEnv,
															 LPXLOPER XL_AsOfDate,	
															 LPXLOPER XL_ForceCurveName,
															 LPXLOPER XL_Currency,
															 LPXLOPER XL_xmlFile,
															 LPXLOPER XL_indexId)
{
	static XLOPER XL_result;
	try
	{
		long objId ;
		std::string pricingEnv ; ExcelTools::convert(XL_PricingEnv,"",pricingEnv);
		std::string forceCurveName ; ExcelTools::convert(XL_ForceCurveName,"",forceCurveName);
		std::string Currency ; ExcelTools::convert(XL_Currency,"",Currency);
		std::string xmlFile ; ExcelTools::convert(XL_xmlFile,"",xmlFile);
		std::string indexId ; ExcelTools::convert(XL_indexId,"",indexId);
		ARM_Date asof ; ExcelTools::convert(XL_AsOfDate,asof);
		
		CCString prevClass;
		CCString stringId = GetLastCurCellEnvValue();
		CCString curClass = LOCAL_VOL_CURVE_LIN_CLASS ;
		long indexId_= LocalPersistent::get().getObjectId(indexId);

		if(!stringId)
		{
			objId = ICMLOCAL_CreateBasketCorrelMkDataFromCalypso(pricingEnv,asof,forceCurveName,Currency,xmlFile,indexId_);
			LocalSetCurCellEnvValue (curClass, objId); 
			stringId = LocalMakeObjectId (objId, curClass);
			ExcelTools::convert(CCSTringToSTLString(stringId),&XL_result); 
			return &XL_result ;
		}

		// cell is not empty
		prevClass = LocalGetStringObjectClass (stringId);
		objId = LocalGetNumObjectId (stringId);
		if ( curClass == prevClass)
		{
			objId = ICMLOCAL_CreateBasketCorrelMkDataFromCalypso(pricingEnv,asof,forceCurveName,Currency,xmlFile,indexId_);
			LocalSetCurCellEnvValue (curClass, objId); 
			stringId = LocalMakeObjectId (objId, curClass);
			ExcelTools::convert(CCSTringToSTLString(stringId),&XL_result); 
			return &XL_result ;
		}

		FreeCurCellContent ();
		objId = ICMLOCAL_CreateBasketCorrelMkDataFromCalypso(pricingEnv,asof,forceCurveName,Currency,xmlFile,indexId_);
		LocalSetCurCellEnvValue (curClass, objId); 
		stringId = LocalMakeObjectId (objId, curClass);
		ExcelTools::convert(CCSTringToSTLString(stringId),&XL_result); 
		return &XL_result ;
		
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
	return &XL_result;

}
__declspec(dllexport) LPXLOPER WINAPI ARM_Credit_GetBasketCorrelMkDataFromCalypso(LPXLOPER XL_PricingEnv,
															 LPXLOPER XL_AsOfDate,	
															 LPXLOPER XL_ForceCurveName,
															 LPXLOPER XL_xmlFile)
{
	static XLOPER XL_result;
	LPXLOPER pxArray;
	try
	{
		ARM_NOCALCIFWIZ();
		
		std::string pricingEnv ; ExcelTools::convert(XL_PricingEnv,"",pricingEnv);
		std::string forceCurveName ; ExcelTools::convert(XL_ForceCurveName,"",forceCurveName);
		std::string xmlFile ; ExcelTools::convert(XL_xmlFile,"",xmlFile);
		ARM_Date asof ; ExcelTools::convert(XL_AsOfDate,asof);

		std::vector<std::string> matus;
		std::vector<double> attach;
		ICM_QMatrix<double> correls;

		ICMLOCAL_GetBasketCorrelMkDataFromCalypso(pricingEnv,asof,forceCurveName,xmlFile,matus,attach,correls) ;
		

		int nbcolumns = attach.size() ;
		int nbrows = matus.size() ;

		//FreeCurcellErr();
		FreeCurCellErr ();
		XL_result.xltype= xltypeMulti ;
		XL_result.val.array.columns = nbcolumns+1;
		XL_result.val.array.rows = nbrows+1;
		XL_result.val.array.lparray = pxArray = (LPXLOPER)GlobalAlloc(GMEM_ZEROINIT, (nbrows+1) * (nbcolumns+1) * sizeof(XLOPER));
		
		char * tmp_str ;
		CCString tmpCCstring ;


		pxArray[XL_Coordonnate2Rank (0, 0, nbcolumns+1)].xltype = xltypeStr;
		pxArray[XL_Coordonnate2Rank (0, 0, nbcolumns+1)].val.str = XL_StrC2StrPascal("");
		pxArray[XL_Coordonnate2Rank (0, 0, nbcolumns+1)].xltype |= xlbitDLLFree;
		
		
		for(int j=1;j<nbcolumns+1;j++)
		{
			pxArray[XL_Coordonnate2Rank (0, j, nbcolumns+1)].xltype = xltypeNum;
			pxArray[XL_Coordonnate2Rank (0, j, nbcolumns+1)].val.num = attach[j-1]/100;

		}

		for(int i = 1; i < nbrows+1; i++)
		{
			tmp_str  =(char *) matus[i-1].c_str();
			tmpCCstring = CCString(tmp_str);
			pxArray[XL_Coordonnate2Rank (i, 0, nbcolumns+1)].xltype = xltypeStr;
			pxArray[XL_Coordonnate2Rank (i, 0, nbcolumns+1)].val.str= XL_StrC2StrPascal(tmpCCstring);
			pxArray[XL_Coordonnate2Rank (i, 0, nbcolumns+1)].xltype |= xlbitDLLFree;
		}


		for( i = 1; i < nbrows+1; i++)
		{
		
			for( j=1;j<nbcolumns+1;j++)
			{
				pxArray[XL_Coordonnate2Rank (i, j, nbcolumns+1)].xltype = xltypeNum;
				pxArray[XL_Coordonnate2Rank (i, j, nbcolumns+1)].val.num = correls(i-1,j-1);

			}
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


__declspec(dllexport) LPXLOPER WINAPI ARM_Credit_CptInterpolDefCurveOLD  (LPXLOPER XL_inCv,
																	LPXLOPER XL_plot,
																	LPXLOPER XL_Date,
																	LPXLOPER XL_slope,
																	LPXLOPER XL_ExactDate)
{
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	LPXLOPER pxArray;

	try {
	ARM_result C_result;
	ARM_NOCALCIFWIZ();

	// C variable
	CCString C_defcurveId;
	CCString C_plot;
	CCString C_plot_default = "NONE";

	double C_slope;
	double C_slope_default = 0;

	double C_date;
	double C_date_default = 0;

	double C_Exactdate;
	double C_Exactdate_default = 0;

	VECTOR<double> output;
	output.clear();

	// error
	static int error;
	static char* reason = "";

	XL_readStrCell(XL_inCv,C_defcurveId," ARM_ERR: curve id: object expected",C_result);
	XL_readStrCellWD(XL_plot,C_plot,C_plot_default," ARM_ERR: Plot: string expected ex: 2Y",C_result);
	XL_readNumCellWD(XL_Date,C_date,C_date_default," ARM_ERR: slope: numeric expected",C_result);
	XL_readNumCellWD(XL_slope,C_slope,C_slope_default," ARM_ERR: slope: numeric expected",C_result);
	XL_readNumCellWD(XL_ExactDate,C_Exactdate,C_Exactdate_default," ARM_ERR: slope: numeric expected",C_result);

	long retCode = ICMLOCAL_CptImplicitSpreadInterpolOLD ((long)LocalGetNumObjectId(C_defcurveId), 
									C_plot,
									C_date,
									C_slope,
									C_Exactdate,	
									output,
									C_result);

	if(retCode == ARM_OK)
	{
		int nbrows = 1;
		int nbcolumns = 2;
		
		FreeCurCellErr ();
		XL_result.xltype = xltypeMulti;
		XL_result.val.array.columns = nbcolumns;
		XL_result.val.array.rows = nbrows; 
		XL_result.val.array.lparray = pxArray = (LPXLOPER)GlobalAlloc (GMEM_ZEROINIT, nbrows * nbcolumns * sizeof (XLOPER));

		pxArray[XL_Coordonnate2Rank (0, 0, nbcolumns)].xltype = xltypeNum;
		pxArray[XL_Coordonnate2Rank (0, 0, nbcolumns)].val.num = output[0];
		pxArray[XL_Coordonnate2Rank (0, 1, nbcolumns)].xltype = xltypeNum;
		pxArray[XL_Coordonnate2Rank (0, 1, nbcolumns)].val.num = output[1]; 
		
	}
	else
	{
		ARM_ERR();
	}

	output.clear();
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


__declspec(dllexport) LPXLOPER WINAPI ARM_Credit_DefCurvePWC_ABS(LPXLOPER XL_date,
															 LPXLOPER XL_Maturities,
															 LPXLOPER XL_Rates,
															 LPXLOPER XL_ircurve,
															 LPXLOPER XL_Recovery,
															 LPXLOPER XL_ccy,
															 LPXLOPER XL_label,
															 LPXLOPER XL_AdjCal,
															 LPXLOPER XL_Cal_Summit,
															 LPXLOPER XL_VolCurve,
															 LPXLOPER XL_Accrued,
															 LPXLOPER XL_CalibrationAlgo,
															 LPXLOPER XL_CalibrationData,
															 LPXLOPER XL_Lag,
															 LPXLOPER XL_Upfront,
															 LPXLOPER XL_RefIds)
{
	static XLOPER XL_result;
	ARM_result C_result;

	try {
	ARM_NOCALCIFWIZ();

	VECTOR<double> C_Rates; 
	VECTOR<double> C_UpFront; 
	VECTOR<double> C_UpFront_default; 

	VECTOR<std::string> C_RefId; 
	VECTOR<std::string> C_RefId_default; 
	VECTOR<long> l_RefId; 

	std::string  C_ircurve;
	long ircurveId = 0;

	std::string  C_ccy;
	std::string C_label;
	std::string C_AdjCal;

	double lag =0;
	double defaultLag = 2;

	qCDS_ADJ l_AdjCal = qCredit_Adjust20; // STDCDS
	double C_Recovery = 0.;

	std::string C_Summit;
	bool C_Cal_Summit = false;

	std::string C_VolCurve;
	long VolCurveId = 0;
	std::string C_VolCurveDefaultValue ="NONE";

	// error
	static int error;
	static char* reason = "";

	ARM_Date AsOf ; ExcelTools::convert(XL_date,AsOf); 
	ExcelTools::convert(XL_ircurve,"",C_ircurve);

	std::vector<std::string> C_matu ;
	ExcelTools::convert(XL_Maturities,C_matu); 
	ExcelTools::convert(XL_Rates,C_Rates);
	ExcelTools::convert(XL_Recovery,0.0,C_Recovery);
	ExcelTools::convert(XL_ccy,"DEFAULT",C_ccy);
	ExcelTools::convert(XL_label,"NONE",C_label);
	ExcelTools::convert(XL_AdjCal,"STDCDS",C_AdjCal);
	ExcelTools::convert(XL_Cal_Summit,"Y",C_Summit);
	ExcelTools::convert(XL_VolCurve,"",C_VolCurve);

	ExcelTools::convert(XL_Upfront,C_UpFront);
	ExcelTools::convert(XL_RefIds,C_RefId);

	qDEFCURVE_CALIB_ALGO calibrationAlgo ; ExcelTools::econvert(XL_CalibrationAlgo,"DICHO",calibrationAlgo); 
	std::string calibrationData ; ExcelTools::convert(XL_CalibrationData,"STD",calibrationData) ;

	l_RefId.resize(C_RefId.size());

	for (int i=0;i<C_RefId.size();i++)
	{l_RefId[i]=(long)LocalGetNumObjectId(C_RefId[i].c_str());}

	if (C_Summit == "Y")  C_Cal_Summit = true;

	if(C_ccy == "DEFAULT")
	{
		ARM_result currencyres;
		ARMLOCAL_ARM_GetDefaultCurrency (currencyres);
		if(currencyres.getRetCode () != ARM_OK)
		{
			ARM_ARG_ERR();
			return (LPXLOPER)&XL_result;
		}
		else
		{
			C_ccy = currencyres.getString ();
		}
	}
	defaultLag= ARM_Currency((const char *)C_ccy.c_str()).GetCreditStartDateLag();

	XL_readNumCellWD(XL_Lag,lag,defaultLag, "ARM_ERR: lag expected ",C_result);
	if((l_AdjCal = ARM_ConvAdjCalCDS (CCString(C_AdjCal.c_str()), C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	double tmp_rate = -1;
	CCString tmp_matu ("NULL");

	if(C_matu.size () != C_Rates.size())
	{
		C_result.setMsg ("ARM_ERR: check your maturities & spreads array");
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	double tmp = 0.;
	int j = 0;
 

//	long retCode;
	CCString prevClass;
	
	CCString curClass = LOCAL_ZERO_CURVE_CDS_CLASS;
	CCString stringId = GetLastCurCellEnvValue ();


	long prevId = ExcelCaller::get().getObjectId(); 
	long objId = ICMLOCAL_ABS_PWCDefaultCurve(AsOf, 
												C_matu, 
												C_Rates, 
												C_Recovery,
												(long)LocalGetNumObjectId(C_ircurve.c_str()),
												C_ccy,
												C_label,
												l_AdjCal,
												C_Cal_Summit,
												// C_Cal_Brent_Solver,
												calibrationAlgo,
												(long)LocalGetNumObjectId(C_VolCurve.c_str()),
												calibrationData,
												(int)lag,
												C_UpFront,
												l_RefId,
												prevId);

	// clean memory
	// C_matu.clear();
	// C_Rates.clear();
	std::string objName = ExcelCaller::get().setObject(objId,LOCAL_ZERO_CURVE_CDS_CLASS); 
	ExcelTools::convert(objName,&XL_result); 
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
