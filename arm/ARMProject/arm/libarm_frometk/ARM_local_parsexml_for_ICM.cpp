
#include <ARM\libarm_local\firstToBeIncluded.h>

#include <ICMKernel\crv\icm_constant_piecewise.h>
//#include <ICMKernel\mod\modelmulticurves.h>
#include <ICMKernel\glob\icm_corrmatrix.h>


#include <atlbase.h>
#include <ARM\libarm_frometk\XMLTools.h>
//#include <ARM\libarm_frometk\arm_local_parsexml_nt_curves.h>
#include <ARM\libarm_frometk\arm_local_parsexml_util.h>
#include "VariantTools.h"

//#include <ICMKernel\inst\icm_cds.h>
//#include <ICMKernel\inst\icm_mez.h>
//#include <ICMKernel\inst\icm_nthtd.h>

//#include <ARM_local_parsexml.h>

//#include <ARM\libarm_local\arm_local_glob.h>
#include <ARM\libarm_frometk\arm_local_parsexml.h>
//#include <ARM\libarm_frometk\arm_local_parsexml_common.h>
//#include <ARM\libarm_frometk\arm_local_parsexml_util.h>
#include <ARM\libarm_frometk\arm_local_etoolkit_for_icm.h>
#include <ARM\libarm_frometk\arm_local_parsexml_for_icm.h>

//#include <ARM\libarm_frometk\arm_local_parsexml_nt_security.h>
#include <ARM\libarm_frometk\PaserManagerUtilities.h>

//#include <ARM\libarm_local\ARM_local_swap.h>

//#include <ARM\local_xlarm\ARM_local_interglob.h>
/// headers to remove definition of va_start and va_end as this is redefined later on!
/// handle this with care after sorting out why this is so!!
//#include <ARM\libarm_local\undef_va_vars.h>
#include <ARM\libarm_frometk\arm_local_paesexml_calypso.h> 

/**
#include <glob\linalg.h>
#include <glob\dates.h>
#include <crv\zeroint.h>
#include <ccy\currency.h>
**/

//#include <util\fromto.h>
//#include <wtypes.h>

//#include <ARM\libarm_frometk\VariantTools.h>

//using namespace MSXML2;

//#include <ARM\libarm_frometk\arm_local_etoolkit.h>


#include <libCCatl\CCatl.h>


#ifndef ARM_LOCAL_ETOOLKIT_DTD
#define ARM_LOCAL_ETOOLKIT_DTD
#define DTD_GETDEFPROBCURVE  "/Response/MKTDATA/MktSpecs/MKTSPEC/MktPoints/MKTPOINT"
#endif
#define PATHFORCED  "/CurveProbabilityList/CurveProbability"

#define ICM_RELEASEXML(argICM_RELEASEXML) \
{ if (argICM_RELEASEXML) argICM_RELEASEXML->Release(); argICM_RELEASEXML=0; }

// --------------------------------------------------------------------------------
// Recuperation de la RecoveryCurve
// --------------------------------------------------------------------------------

ARM_Vector* ARMLOCAL_RecoveryCurve (const char* chaineXML_Recovery)
{

	VARIANT_BOOL bOK;
	HRESULT hr;
	MSXML2::IXMLDOMDocument *XMLDoc = NULL;

	ARM_Vector* pRecovery = NULL;
	char** pPLOT = NULL;

	try
	{
		hr = CoInitialize(NULL); 
//		SUCCEEDED(hr) ? 0 : throw hr;
		hr = CoCreateInstance(__uuidof(MSXML2::DOMDocument30), NULL, CLSCTX_INPROC_SERVER, __uuidof(MSXML2::IXMLDOMDocument), (void**)&XMLDoc);
		SUCCEEDED(hr) ? 0 : throw hr;

		//JLA _bstr_t* tmpChaine = new _bstr_t(chaineXML_Recovery);
		_bstr_t* tmpChaine = new _bstr_t ;
		VariantTools::convert(std::string(chaineXML_Recovery),*tmpChaine); 

		XMLDoc->loadXML(*tmpChaine, &bOK);

		delete tmpChaine;
	}
	catch(...)
	{
		ICM_RELEASEXML(XMLDoc) ;
		hr = S_FALSE;
		throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
         "Pb in creating XML document for RecoveryCurve");
	}

	MSXML2::IXMLDOMNodeList * resultList = NULL;
	MSXML2::IXMLDOMNode * listItem = NULL, * theNode = NULL;
	long nbNodes;

	try
	{
		_bstr_t* Info = new _bstr_t("Response/CURVEHEAD/CurvePts/CURVE");

		if (XMLDoc->selectNodes(*Info, &resultList) == S_OK)
		{
			resultList->get_length(&nbNodes);

			if (nbNodes == 0)
			{
				hr = S_FALSE;
				throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
				 "Invalid XML string for Def Prob Curve");
			}

			pRecovery = new ARM_Vector(nbNodes,0.);

			for (long indexNode=0 ; indexNode<nbNodes ; indexNode++)
			{
				hr=resultList->get_item(indexNode, &listItem);
				if (hr==S_OK && listItem!=NULL)
				{
					_bstr_t* bRate =  new _bstr_t("Rate");
					listItem->selectSingleNode(*bRate,&theNode);
					if (theNode!=NULL)
					{
						BSTR resultat = NULL;
						theNode->get_text(&resultat);	//printf("resultat : %S\n", resultat);

						_bstr_t ff(resultat,false);
						char * ff1=(char *)ff;

						double rate = atof(ff1);	//printf("rate : %.13lf\n", rate);

						pRecovery->InitElt(indexNode,rate);

						theNode->Release();
						theNode=NULL;
						if (resultat) SysFreeString(resultat);
					}
					delete bRate;

					listItem->Release();
					listItem=NULL;
				}
			}

		}
		delete Info;


		if (listItem)
		{
			listItem->Release();
			listItem = NULL;
		}

		if (theNode)
		{
			theNode->Release();
			theNode = NULL;
		}
		if (resultList)
		{
			resultList->Release();
			resultList = NULL;
		}
		ICM_RELEASEXML(XMLDoc); 

		hr = S_OK;

		return pRecovery;
	}

	catch(Exception& )
	{
		if (pRecovery)
			delete[] pRecovery;
		
		if (listItem) listItem->Release();
		if (resultList) resultList->Release();
		if (theNode) theNode->Release();
		ICM_RELEASEXML(XMLDoc); 

		hr = S_FALSE;
		throw; 
		
	}
	catch(...)
	{
		if (pRecovery)
			delete[] pRecovery;
		
		if (listItem) listItem->Release();
		if (resultList) resultList->Release();
		if (theNode) theNode->Release();
		ICM_RELEASEXML(XMLDoc); 

		hr = S_FALSE;
		throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
         "Error in XML parsing for ZC");
	}
}
// definition structure contenant spread + matu
struct DPPoint{
	double spread;
	std::string matu;
	} ;


ICM_DefaultCurve* ARMLOCAL_ParseDefProbCurveCalypso (const std::string& chaineXML_defcurve,
											const ARM_Date&CurveDate, 
											ARM_ZeroCurve *zcpy,
											const std::string& label,
											const std::string& Seniority,
											const std::string& Currency,
											const std::string& ForceCurveName,
											map<string,ARM_ZeroCurve*>* MercureMap)
{


	vector <string> pPLOT ;
	vector <double> pSPREAD ;
	ARM_ZeroCurve* IRCurve = NULL;
	// ARM_Currency* ccy = NULL;
	std::string ccy; 
	long Frequency = 0;
	double Recovery = CREDIT_DEFAULT_VALUE;
	std::string  IndexName = "NONE";
	std::string	 curveName ;

	std::map <ARM_Date,DPPoint>	 DPMap;
	std::map<ARM_Date,DPPoint>::iterator DPMapIterator;
	struct DPPoint quote ;
	ARM_Date keydate = ARM_Date();
	qCDS_ADJ AIMMADJ = qCredit_Adjust20;
	
	ARM_Vector* mktData = NULL;
	ARM_Vector* cdsData = NULL;

	MSXML2::IXMLDOMDocumentPtr XMLDoc;
	XMLDoc = XMLTools::LoadXML(chaineXML_defcurve);
	#ifdef _DEBUG
		_variant_t tmp("c:\\temp\\Calypso_defprobcurve.xml");
		XMLDoc->save(tmp);
	#endif


	std::stringstream sstr;
	if ( ! ForceCurveName.empty())
	{
		sstr << "/CurveProbabilityList/CurveProbability[CurveName='"<< ForceCurveName <<"']";
	}
	else
	{
		sstr << "/CurveProbabilityList/CurveProbability[Issuer/ShortName='"<< label <<"']";
	}

	//sstr << "/CurveProbabilityList/CurveProbability";
	
	MSXML2::IXMLDOMNodePtr theNode;
	theNode = XMLTools::selectSingleNode(XMLDoc,sstr.str());
	// the node contains the curveprobability node

	MSXML2::IXMLDOMNodeListPtr Underlyings ;
	std::stringstream sstrUnderlyings;
	sstrUnderlyings << sstr.str() << "/Underlyings/Underlying" ;
	Underlyings = XMLTools::selectNodes(XMLDoc,sstrUnderlyings.str());

	long NbNodes(0);
	Underlyings ->get_length(&NbNodes);
	// getting frequency from firt item
	std::string sFreq;
	XMLTools::convert(XMLTools::selectSingleNode(theNode,"Underlyings/Underlying[0]/Frequency"),sFreq);
	Frequency = FromSummitFreqToARMFreq(sFreq.c_str());

	pPLOT.resize(NbNodes);
	pSPREAD.resize(NbNodes);

	for ( int i=0;i< NbNodes;i++)
	{
		string theplot;
		double thespread;
		
		XMLTools::convert(XMLTools::selectSingleNode(XMLTools::get_item(Underlyings,i),"Maturity"),theplot);
		
		XMLTools::convert(XMLTools::selectSingleNode(XMLTools::get_item(Underlyings,i),"Quote/Value"),thespread); 
		quote.matu = theplot;
		quote.spread = (thespread/10000.);
		keydate = ARM_Date() ;
		keydate.AddPeriod(theplot);
		DPMap[keydate] = quote;
	}
	
	XMLTools::convert(XMLTools::selectSingleNode(theNode,"CurveInfo/RecoveryRate"),Recovery);
	XMLTools::convert(XMLTools::selectSingleNode(theNode,"CurveInfo/Currency"),ccy);
	XMLTools::convert(XMLTools::selectSingleNode(theNode,"CurveInfo/RisklessCurve"),curveName);
	std::string adjRule;
	XMLTools::convert(XMLTools::selectSingleNode(theNode,"CurveInfo/MaturitiesAdjustment"),adjRule);
	
	if (adjRule == "NONE") 
	{
		XMLTools::convert(XMLTools::selectSingleNode(theNode,"Underlyings/Underlying[0]/RollMaturity"),adjRule);
		if (adjRule=="QTR") AIMMADJ = qCredit_Adjust20 ;
		else if (adjRule=="SA") AIMMADJ = qCredit_Adjust20SA ;
		else AIMMADJ = qCredit_Default ;
	}
	else
	{
		ICM_EnumsCnv::cnv(adjRule,AIMMADJ); 
	}
		// Fill plots and pspread from the map
		int j =0;
		for(DPMapIterator = DPMap.begin(); DPMapIterator != DPMap.end();DPMapIterator++)
		{
			pPLOT[j] = (*DPMapIterator).second.matu;
			pSPREAD[j] =(*DPMapIterator).second.spread;
			j++;
		}

	try
	{
		vector<string> psMatu ;
		for (int i=0; i<NbNodes; i++)
		{
			psMatu.push_back(string((const char*) pPLOT[i].c_str() ));
		}
	
		mktData = new ARM_Vector(pSPREAD);

		if (zcpy)
			IRCurve = (ARM_ZeroCurve*) zcpy->Clone();
		else
		{

			if (MercureMap)
			{
			IRCurve = (ARM_ZeroCurve*) (*MercureMap)[ccy]->Clone();
			}
			else
			{
				std::string xmlResponse ;
				
				ARM_CalypsoToolkit::GetCurveZero("","","",curveName,SUMMIT_DEFAULT_CURVE,CurveDate,"",xmlResponse); 
				IRCurve = ARMLOCAL_ParseXMLForCalypsoZC(xmlResponse, CurveDate,ccy,0);

			}
		}


		ICM_DefaultCurve* pwcdefault = (ICM_DefaultCurve*) new ICM_Constant_Piecewise(CurveDate,
																					psMatu,
																					mktData, 
																					Recovery,
																					IRCurve,
																					K_ADJUSTED,	// intRule
																					K_ADJUSTED,	// adjStartDate
																					AIMMADJ,
																					ccy,
																					label,
																					true,//Summit Curve
																					//2 NULL,//Default Curve
																					NULL,//Volatility Curve
																					Frequency,
																					qDEFCURVE_DICHO,
																					"STD", 
																					ARM_Currency(ccy.c_str()).GetCreditStartDateLag()
																					,ICM_Parameters());
		if (IRCurve)
			delete IRCurve;
		IRCurve=NULL;


		pPLOT.clear();
		pSPREAD.clear();

		if (mktData)
			delete mktData;
		mktData = NULL;

		return (pwcdefault);


	}

	catch(...)
	{
		

		if (IRCurve)
			delete IRCurve;
		IRCurve=NULL;

		pPLOT.clear();
		pSPREAD.clear(); 
		if (mktData)
			delete mktData;
		mktData = NULL;

		ICM_DefaultCurve* pwcdefault=NULL ;
		throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
         "Pb in creating the DefProbCurve");

	}
	
}

// --------------------------------------------------------------------------------
// Recuperation de la DefProbCurve
// --------------------------------------------------------------------------------

ICM_DefaultCurve* ARMLOCAL_ParseDefProbCurve (const char* chaineXML_defcurve,
											ARM_Date CurveDate, 
											ARM_ZeroCurve *zcpy,
											const char* label,
											map<string,ARM_ZeroCurve*>* MercureMap)
{

	VARIANT_BOOL bOK;
	HRESULT hr;
	MSXML2::IXMLDOMDocument *XMLDoc = NULL;

	char** pPLOT = NULL;
	double* pSPREAD = NULL;
	ARM_ZeroCurve* IRCurve = NULL;
	// ARM_Currency* ccy = NULL;
	std::string ccy; 
	long Frequency = 0;
	double Recovery = CREDIT_DEFAULT_VALUE;
	CCString IndexName = "NONE";
	qCDS_ADJ AIMMADJ = qCredit_Adjust20;

	try
	{
		hr = CoInitialize(NULL); 
//		SUCCEEDED(hr) ? 0 : throw hr;
		hr = CoCreateInstance(__uuidof(MSXML2::DOMDocument30), NULL, CLSCTX_INPROC_SERVER, __uuidof(MSXML2::IXMLDOMDocument), (void**)&XMLDoc);
		SUCCEEDED(hr) ? 0 : throw hr;

		_bstr_t  tmpChaine (chaineXML_defcurve);

		XMLDoc->loadXML(tmpChaine, &bOK);

		#ifdef _DEBUG
		CCString chtmp = (CCString)"c:\\temp\\" + (CCString)label + (CCString)".xml";
		VARIANT v_chtmp;
		CCString2VARIANT (chtmp, &v_chtmp);
		XMLDoc->save(v_chtmp);
		#endif

		// delete tmpChaine;
	}
	catch(...)
	{
		ICM_RELEASEXML(XMLDoc); 
		hr = S_FALSE;
		throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
         "Pb in creating XML document for DefProbCurve");
	}

	MSXML2::IXMLDOMNodeList * resultList = NULL;
	MSXML2::IXMLDOMNode * listItem = NULL, * theNode = NULL;
	long nbNodes;

	try
	{
		//---------------------------------------------------------------
		//Récupération des spreads et des Tenors
		//---------------------------------------------------------------
		_bstr_t  Info(DTD_GETDEFPROBCURVE);
		if (XMLDoc->selectNodes(Info, &resultList) == S_OK)
		{
			resultList->get_length(&nbNodes);

			if (nbNodes == 0)
			{
				hr = S_FALSE;
				throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
				 "Invalid XML string for Def Prob Curve");
			}

			pPLOT = new char*[nbNodes];
			pSPREAD = new double[nbNodes];

			for (long indexNode=0 ; indexNode<nbNodes ; indexNode++)
			{
				hr=resultList->get_item(indexNode, &listItem);
				if (hr==S_OK && listItem!=NULL)
				{
					_bstr_t* bDate = new _bstr_t("Date");
					listItem->selectSingleNode(*bDate, &theNode);
					if (theNode!=NULL)
					{
						BSTR resultat = NULL;
						theNode->get_text(&resultat);

						_bstr_t ff(resultat,false);
						char * ff1=(char *)ff;

						pPLOT[indexNode] = new char[5];
						strcpy(pPLOT[indexNode],ff1); 

						theNode->Release();
						theNode=NULL;
						if (resultat) SysFreeString(resultat);
					}
					delete bDate;

					_bstr_t* bRate =  new _bstr_t("MidRate");
					listItem->selectSingleNode(*bRate,&theNode);
					if (theNode!=NULL)
					{
						BSTR resultat = NULL;
						theNode->get_text(&resultat);	//printf("resultat : %S\n", resultat);

						_bstr_t ff(resultat,false);
						char * ff1=(char *)ff;

						double rate =  atof(ff1);	//printf("rate : %.13lf\n", rate);

						pSPREAD[indexNode] = rate;

						theNode->Release();
						theNode=NULL;
						if (resultat) SysFreeString(resultat);
					}
					delete bRate;

					listItem->Release();
					listItem=NULL;
				}
			}

		}
		// delete Info;

		//---------------------------------------------------------------
		//Récupération Recovery/Currency/Frequency
		//---------------------------------------------------------------

 		if (XMLDoc->selectNodes((_bstr_t)("//COMMDATA"), &resultList) == S_OK)
		{

			hr=resultList->get_item(0, &listItem);
			Recovery = GetDoubleFromXMLNode (listItem,"Value");
		}

		if (listItem)
		{
			listItem->Release();
			listItem = NULL;
		}
		if (resultList)
		{
			resultList->Release();
			resultList = NULL;
		}

 		if (XMLDoc->selectNodes((_bstr_t)("//MKTDATA"), &resultList) == S_OK)
		{
			hr=resultList->get_item(0, &listItem);
			CCString Meth = GetCCStringFromXMLNode (listItem,"Method");
			if (Meth == "RAW") AIMMADJ = qCredit_Default;
			else if (Meth == "cCDS") AIMMADJ = qCredit_Adjust20;
		}

		if (listItem)
		{
			listItem->Release();
			listItem = NULL;
		}
		if (resultList)
		{
			resultList->Release();
			resultList = NULL;
		}


 		if (XMLDoc->selectNodes((_bstr_t)("//MKTSPEC"), &resultList) == S_OK)
		{
			hr=resultList->get_item(0, &listItem);
			// CCString Meth = GetCCStringFromXMLNode (listItem,"GenAssetId");
			std::string Meth = GetStringFromXMLNode(listItem,"GenAssetId"); 

			if (AIMMADJ != qCredit_Default)
			{
				ICM_EnumsCnv::cnv(Meth,AIMMADJ); 
				/**
			if (Meth == "CDSDRULE") AIMMADJ = qCredit_Adjust20;
			else if (Meth == "CDSDIND") AIMMADJ = qCredit_CDSDIND;
			else if (Meth == "CDSDTRX") AIMMADJ = qCredit_CDSDTRX;
			else if (Meth == "CDSDIDZ") AIMMADJ = qCredit_CDSINDZ;
			*/ 
			}
		}

		if (listItem)
		{
			listItem->Release();
			listItem = NULL;
		}
		if (resultList)
		{
			resultList->Release();
			resultList = NULL;
		}

 		if (XMLDoc->selectNodes((_bstr_t)("//GASSET"), &resultList) == S_OK)
		{

			hr=resultList->get_item(0, &listItem);
	
			//Recuperation de frequence
			listItem->selectSingleNode(_bstr_t("SCHED_Pay_Freq"), &theNode);
			if (theNode!=NULL)
			{
				BSTR resultat = NULL;
				theNode->get_text(&resultat);

				_bstr_t ff(resultat,false);
				char * ff1=(char *)ff;
				Frequency = FromSummitFreqToARMFreq((const char*) ff1);

				theNode->Release();
				theNode=NULL;
				if (resultat) SysFreeString(resultat);
			}
			if (theNode)
			{
			theNode->Release();
			theNode = NULL;
			}

			//Recuperation de la ccy
			listItem->selectSingleNode(_bstr_t("Ccy"), &theNode);
			if (theNode!=NULL)
			{
				BSTR resultat = NULL;
				theNode->get_text(&resultat);

				_bstr_t ff(resultat,false);
				ccy = ff ; 
				// char * ff1=(char *)ff;
				// ccy = new ARM_Currency((const char*)ff1);

				theNode->Release();
				theNode=NULL;
				if (resultat) SysFreeString(resultat);
			}
			if (theNode)
			{
			theNode->Release();
			theNode = NULL;
			}

			//Recuperation de l'index
			listItem->selectSingleNode(_bstr_t("INTEREST_UnderIndex"), &theNode);
			if (theNode!=NULL)
			{
				BSTR resultat = NULL;
				theNode->get_text(&resultat);

				_bstr_t ff(resultat,false);
				char * ff1=(char *)ff;
				IndexName = (CCString)ff1;

				theNode->Release();
				theNode=NULL;
				if (resultat) SysFreeString(resultat);
			}
			if (theNode)
			{
			theNode->Release();
			theNode = NULL;
			}

		}

		if (listItem)
		{
			listItem->Release();
			listItem = NULL;
		}

		if (theNode)
		{
			theNode->Release();
			theNode = NULL;
		}
		if (resultList)
		{
			resultList->Release();
			resultList = NULL;
		}
		ICM_RELEASEXML(XMLDoc); 
		hr = S_OK;
	}

	catch(...)
	{
		for (int j=0; j<nbNodes; j++)
		{	
			delete[] pPLOT[j];
		}

		if (pPLOT)
			delete[] pPLOT;
		if (pSPREAD)
			delete[] pSPREAD;
		
		if (listItem) listItem->Release();
		if (resultList) resultList->Release();
		if (theNode) theNode->Release();
		ICM_RELEASEXML(XMLDoc); 
		hr = S_FALSE;
		throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
         "Error in XML parsing for ZC");
	}

	ARM_Vector* mktData = NULL;
	ARM_Vector* cdsData = NULL;

	try
	{
		vector<string> psMatu ;
		for (int i=0; i<nbNodes; i++)
		{
			psMatu.push_back(string((const char*) pPLOT[i]));
		}

		//Recuperation du niveau de recovery
		CCString XML_output;
		CCString XML_messagelist;

		mktData = new ARM_Vector(nbNodes ,pSPREAD);

		if (zcpy)
			IRCurve = (ARM_ZeroCurve*) zcpy->Clone();
		else
		{
			CCString xmlResponse;

			if (MercureMap)
			{
			IRCurve = (ARM_ZeroCurve*) (*MercureMap)[ccy]->Clone();
			}
			else
			{
			xmlResponse = etoolkit_getXMLMYAndZCFromSummit(IndexName,ccy.c_str(),SUMMIT_DEFAULT_CURVE,CurveDate);
			int interp = ARMLOCAL_ParseXMLForInterpolator(xmlResponse);
			xmlResponse = etoolkit_getXMLZCFromSummit(IndexName,ccy.c_str(),SUMMIT_DEFAULT_CURVE,CurveDate);
			IRCurve = ARMLOCAL_ParseXMLForZC(xmlResponse, CurveDate,ccy.c_str(),interp);
			}
		}


		ICM_DefaultCurve* pwcdefault = (ICM_DefaultCurve*) new ICM_Constant_Piecewise(CurveDate,
																					psMatu,
																					mktData, 
																					Recovery,
																					IRCurve,
																					K_ADJUSTED,	// intRule
																					K_ADJUSTED,	// adjStartDate
																					AIMMADJ,
																					ccy,
																					label,
																					true,//Summit Curve
																					//2 NULL,//Default Curve
																					NULL,//Volatility Curve
																					Frequency,qDEFCURVE_DICHO,
																					"STD", ARM_Currency(ccy.c_str()).GetCreditStartDateLag(),
																					ICM_Parameters());

//		if (ccy)
//			delete ccy;
//		ccy = NULL;

		if (IRCurve)
			delete IRCurve;
		IRCurve=NULL;

		for (int j=0; j<nbNodes; j++)
		{
			delete[] pPLOT[j];
		}

		delete[] pPLOT;
		delete[] pSPREAD;

		if (mktData)
			delete mktData;
		mktData = NULL;

		return (pwcdefault);
	}
	catch(Exception& )
	{
		if (pPLOT)
		for (int j=0; j<nbNodes; j++)
			if (pPLOT[j])
			delete[] pPLOT[j];

		delete[] pPLOT;
		delete[] pSPREAD;

		if (IRCurve)
			delete IRCurve;
		IRCurve=NULL;

		if (mktData)
			delete mktData;
		mktData = NULL;

//		if (ccy)
//			delete ccy;
//		ccy = NULL;

		hr = S_FALSE;
		throw; 
	}
	catch(...)
	{
		if (pPLOT)
		for (int j=0; j<nbNodes; j++)
			if (pPLOT[j])
			delete[] pPLOT[j];

		if (IRCurve)
			delete IRCurve;
		IRCurve=NULL;

		delete[] pPLOT;
		delete[] pSPREAD;

		if (mktData)
			delete mktData;
		mktData = NULL;

//		if (ccy)
//			delete ccy;
//		ccy = NULL;

		hr = S_FALSE;
		throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
         "Error in creating DefProbCurve");
	}
}


// --------------------------------------------------------------------------------
// Recuperation des MktData de la DefProbCurve
// --------------------------------------------------------------------------------

long  ARMLOCAL_ParseDefProbCurveMktData (const char* chaineXML,
										VECTOR<CCString>& matu,
										VECTOR<double>& spread,
										VECTOR<double>& recovery,
										CCString& currency,
										CCString& indexname,
										long& AIMMADJ)
{

	VARIANT_BOOL bOK;
	HRESULT hr;
	MSXML2::IXMLDOMDocument *XMLDoc = NULL;

	char** pPLOT = NULL;
	double* pSPREAD = NULL;
	ARM_ZeroCurve* IRCurve = NULL;
	ARM_Currency* ccy = NULL;
	long Frequency = 0;
	double Recovery = CREDIT_DEFAULT_VALUE;
	CCString IndexName = "NONE";
	//qCDS_ADJ AIMMADJ = qCredit_Adjust20;

	try
	{
		hr = CoInitialize(NULL); 
//		SUCCEEDED(hr) ? 0 : throw hr;
		hr = CoCreateInstance(__uuidof(MSXML2::DOMDocument30), NULL, CLSCTX_INPROC_SERVER, __uuidof(MSXML2::IXMLDOMDocument), (void**)&XMLDoc);
		SUCCEEDED(hr) ? 0 : throw hr;

		//JLA _bstr_t* tmpChaine = new _bstr_t(chaineXML);
		_bstr_t* tmpChaine = new _bstr_t ;
		VariantTools::convert(std::string(chaineXML),*tmpChaine); 

		XMLDoc->loadXML(*tmpChaine, &bOK);

		#ifdef _DEBUG
		CCString chtmp = (CCString)"c:\\temp\\" + (CCString)"MKT_DATA" + (CCString)".xml";
		VARIANT v_chtmp;
		CCString2VARIANT (chtmp, &v_chtmp);
		XMLDoc->save(v_chtmp);
		#endif

		delete tmpChaine;
	}
	catch(...)
	{
		ICM_RELEASEXML(XMLDoc); 
		hr = S_FALSE;
		throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
         "Pb in creating XML document for DefProbCurve");
	}

	MSXML2::IXMLDOMNodeList * resultList = NULL;
	MSXML2::IXMLDOMNode * listItem = NULL, * theNode = NULL;
	long nbNodes;

	try
	{
		//---------------------------------------------------------------
		//Récupération des spreads et des Tenors
		//---------------------------------------------------------------
		_bstr_t* Info = new _bstr_t(DTD_GETDEFPROBCURVE);
		if (XMLDoc->selectNodes(*Info, &resultList) == S_OK)
		{
			resultList->get_length(&nbNodes);

			if (nbNodes == 0)
			{
				hr = S_FALSE;
				throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
				 "Invalid XML string for Def Prob Curve");
			}

			pPLOT = new char*[nbNodes];
			pSPREAD = new double[nbNodes];

			for (long indexNode=0 ; indexNode<nbNodes ; indexNode++)
			{
				hr=resultList->get_item(indexNode, &listItem);
				if (hr==S_OK && listItem!=NULL)
				{
					_bstr_t* bDate = new _bstr_t("Date");
					listItem->selectSingleNode(*bDate, &theNode);
					if (theNode!=NULL)
					{
						BSTR resultat = NULL;
						theNode->get_text(&resultat);

						_bstr_t ff(resultat,false);
						char * ff1=(char *)ff;

						pPLOT[indexNode] = new char[5];
						strcpy(pPLOT[indexNode],ff1); 

						theNode->Release();
						theNode=NULL;
						if (resultat) SysFreeString(resultat);
					}
					delete bDate;

					_bstr_t* bRate =  new _bstr_t("MidRate");
					listItem->selectSingleNode(*bRate,&theNode);
					if (theNode!=NULL)
					{
						BSTR resultat = NULL;
						theNode->get_text(&resultat);	//printf("resultat : %S\n", resultat);

						_bstr_t ff(resultat,false);
						char * ff1=(char *)ff;

						double rate =  atof(ff1);	//printf("rate : %.13lf\n", rate);

						pSPREAD[indexNode] = rate;

						theNode->Release();
						theNode=NULL;
						if (resultat) SysFreeString(resultat);
					}
					delete bRate;

					listItem->Release();
					listItem=NULL;
				}
			}

		}
		delete Info;

		//---------------------------------------------------------------
		//Récupération Recovery/Currency/Frequency
		//---------------------------------------------------------------

 		if (XMLDoc->selectNodes((_bstr_t)("//COMMDATA"), &resultList) == S_OK)
		{

			hr=resultList->get_item(0, &listItem);
			Recovery = GetDoubleFromXMLNode (listItem,"Value");
		}

		if (listItem)
		{
			listItem->Release();
			listItem = NULL;
		}
		if (resultList)
		{
			resultList->Release();
			resultList = NULL;
		}

 		if (XMLDoc->selectNodes((_bstr_t)("//MKTDATA"), &resultList) == S_OK)
		{
			hr=resultList->get_item(0, &listItem);
			CCString Meth = GetCCStringFromXMLNode (listItem,"Method");
			AIMMADJ = FromSummitAdjCDSToARM(Meth);
			/*if (Meth == "RAW") {AIMMADJ = qCredit_Default; AdjCDS = "NONE";}
			else if (Meth == "cCDS") AIMMADJ = qCredit_Adjust20;*/
		}

		if (listItem)
		{
			listItem->Release();
			listItem = NULL;
		}
		if (resultList)
		{
			resultList->Release();
			resultList = NULL;
		}


 		if (XMLDoc->selectNodes((_bstr_t)("//MKTSPEC"), &resultList) == S_OK)
		{
			hr=resultList->get_item(0, &listItem);
			CCString Meth = GetCCStringFromXMLNode (listItem,"GenAssetId");

			if (AIMMADJ != qCredit_Default)
			{AIMMADJ = FromSummitAdjCDSToARM(Meth);
			/*	if (Meth == "CDSDRULE") {AIMMADJ = qCredit_Adjust20;AdjCDS = "STDCDS";}
			else if (Meth == "CDSDIND") {AIMMADJ = qCredit_CDSDIND;AdjCDS = "STDIDX";}
			else if (Meth == "CDSDTRX") {AIMMADJ = qCredit_CDSDTRX;AdjCDS = "STDTRX";}
			else if (Meth == "CDSDIDZ") {AIMMADJ = qCredit_CDSINDZ;AdjCDS = "STDIDZ";}*/
			}
		}

		if (listItem)
		{
			listItem->Release();
			listItem = NULL;
		}
		if (resultList)
		{
			resultList->Release();
			resultList = NULL;
		}



 		if (XMLDoc->selectNodes((_bstr_t)("//GASSET"), &resultList) == S_OK)
		{

			hr=resultList->get_item(0, &listItem);
	
			//Recuperation de frequence
			listItem->selectSingleNode(_bstr_t("SCHED_Pay_Freq"), &theNode);
			if (theNode!=NULL)
			{
				BSTR resultat = NULL;
				theNode->get_text(&resultat);

				_bstr_t ff(resultat,false);
				char * ff1=(char *)ff;
				Frequency = FromSummitFreqToARMFreq((const char*) ff1);

				theNode->Release();
				theNode=NULL;
				if (resultat) SysFreeString(resultat);
			}
			if (theNode)
			{
			theNode->Release();
			theNode = NULL;
			}

			//Recuperation de la ccy
			listItem->selectSingleNode(_bstr_t("Ccy"), &theNode);
			if (theNode!=NULL)
			{
				BSTR resultat = NULL;
				theNode->get_text(&resultat);

				_bstr_t ff(resultat,false);
				char * ff1=(char *)ff;
				ccy = new ARM_Currency((const char*)ff1);

				theNode->Release();
				theNode=NULL;
				if (resultat) SysFreeString(resultat);
			}
			if (theNode)
			{
			theNode->Release();
			theNode = NULL;
			}

			//Recuperation de l'index
			listItem->selectSingleNode(_bstr_t("INTEREST_UnderIndex"), &theNode);
			if (theNode!=NULL)
			{
				BSTR resultat = NULL;
				theNode->get_text(&resultat);

				_bstr_t ff(resultat,false);
				char * ff1=(char *)ff;
				IndexName = (CCString)ff1;

				theNode->Release();
				theNode=NULL;
				if (resultat) SysFreeString(resultat);
			}
			if (theNode)
			{
			theNode->Release();
			theNode = NULL;
			}

		}

		if (listItem)
		{
			listItem->Release();
			listItem = NULL;
		}

		if (theNode)
		{
			theNode->Release();
			theNode = NULL;
		}
		if (resultList)
		{
			resultList->Release();
			resultList = NULL;
		}
		ICM_RELEASEXML(XMLDoc); 

		hr = S_OK;
	}

	catch(...)
	{
		for (int j=0; j<nbNodes; j++)
		{	
			delete[] pPLOT[j];
		}

		if (pPLOT)
			delete[] pPLOT;
		if (pSPREAD)
			delete[] pSPREAD;
		
		if (listItem) listItem->Release();
		if (resultList) resultList->Release();
		if (theNode) theNode->Release();
		ICM_RELEASEXML(XMLDoc); 

		hr = S_FALSE;
		throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
         "Error in XML parsing for ZC");
	}

	ARM_Vector* mktData = NULL;

	try
	{

		for (int i=0; i<nbNodes; i++)
		{
			matu.push_back(pPLOT[i]);
			spread.push_back(pSPREAD[i]*10000.);
			recovery.push_back(Recovery);
		}

		currency = (CCString)ccy->GetCcyName();
		indexname = (CCString)IndexName;

		if (ccy)
			delete ccy;
		ccy = NULL;

		for (int j=0; j<nbNodes; j++)
			delete[] pPLOT[j];

		delete[] pPLOT;
		delete[] pSPREAD;

		if (mktData)
			delete mktData;
		mktData = NULL;

		return ARM_OK;

	}
	catch(Exception& )
	{
		{
		if (pPLOT)
		for (int j=0; j<nbNodes; j++)
			if (pPLOT[j])
				delete[] pPLOT[j];
		}

		delete[] pPLOT;
		delete[] pSPREAD;

		if (mktData)
			delete mktData;
		mktData = NULL;

		if (ccy)
			delete ccy;
		ccy = NULL;

		hr = S_FALSE;
		throw; 
	}
	catch(...)
	{
		{
		if (pPLOT)
		for (int j=0; j<nbNodes; j++)
			if (pPLOT[j])
				delete[] pPLOT[j];
		}

		delete[] pPLOT;
		delete[] pSPREAD;

		if (mktData)
			delete mktData;
		mktData = NULL;

		if (ccy)
			delete ccy;
		ccy = NULL;

		hr = S_FALSE;
		throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
         "Error in creating DefProbCurve");
	}
}



// --------------------------------------------------------------------------------
// Recuperation des MktData de la DefProbCurve
// --------------------------------------------------------------------------------

void ARMLOCAL_ParseDefProbCurveMktDataCalypso (const std::string& chaineXML,
											   std::vector<std::string>& matus,
											   std::vector<double>& spreads,
											   double & recovery,
											   std::string& currency,
											   std::string& indexname,
											   qCDS_ADJ & adj)
{
	MSXML2::IXMLDOMDocumentPtr xmlDoc = XMLTools::CreateDOMDocument30(); 
	VARIANT_BOOL bOK;
	_bstr_t tmp; VariantTools::convert(chaineXML,tmp); 
	xmlDoc->loadXML(tmp,&bOK); 
	
	std::map <ARM_Date,DPPoint>	 DPMap;
	std::map<ARM_Date,DPPoint>::iterator DPMapIterator;
	struct DPPoint quote ;
	ARM_Date keydate = ARM_Date();

	MSXML2::IXMLDOMNodeListPtr resultList; 
	resultList = XMLTools::selectNodes(xmlDoc,"/CurveProbabilityList/CurveProbability/Underlyings/Underlying") ;
	long nbNodes(0);
	resultList->get_length(&nbNodes); 
	if (nbNodes==0) 
		ICMTHROW(ERR_INVALID_ARGUMENT,"ARMLOCAL_ParseDefProbCurveMktDataCalypso: no underlyings found..."); 
	matus.resize(nbNodes); 
	spreads.resize(nbNodes); 
	
	string theplot = "";
	double thespread = 0.;
	for(long i=0;i<nbNodes;i++) 
	{
		MSXML2::IXMLDOMNodePtr item = XMLTools::get_item(resultList,i); 
		MSXML2::IXMLDOMNodePtr item2 = XMLTools::selectSingleNode(item,"Maturity"); 
		XMLTools::convert(item2,theplot) ;
		quote.matu = theplot;
		item2 = XMLTools::selectSingleNode(item,"Quote/Value"); 
		XMLTools::convert(item2, thespread);
		quote.spread = thespread;
	
		keydate = ARM_Date() ;
		keydate.AddPeriod(theplot);
		DPMap[keydate] = quote;
	}
	int j =0;
	for(DPMapIterator = DPMap.begin(); DPMapIterator != DPMap.end();DPMapIterator++)
	{
		matus[j] = (*DPMapIterator).second.matu;
		spreads[j] =(*DPMapIterator).second.spread;
		j++;
	}

	MSXML2::IXMLDOMNodePtr item= XMLTools::selectSingleNode(xmlDoc,"/CurveProbabilityList/CurveProbability/CurveInfo/RecoveryRate"); 
	XMLTools::convert(item,recovery) ;
	item= XMLTools::selectSingleNode(xmlDoc,"/CurveProbabilityList/CurveProbability/CurveInfo/Currency"); 
	XMLTools::convert(item,currency) ;
	indexname = "EURIB" ;
	adj=qCredit_Adjust20; 
	item= XMLTools::selectSingleNode(xmlDoc,"/CurveProbabilityList/CurveProbability/CurveInfo/Index"); 
	std::string isIndex; 
	XMLTools::convert(item,isIndex); 
	/*if (isIndex=="No")
	{
		item= XMLTools::selectSingleNode(xmlDoc,"/CurveProbabilityList/CurveProbability/Underlyings/Underlying[0]"); 
		std::string adjrule; 
		XMLTools::convert(item,adjrule); 
		if (adjrule=="QTR") 	adj=qCredit_Adjust20;
		else if (adjrule=="SA") 	adj=qCredit_Adjust20SA;
		else 	adj=qCredit_Default;
	}
	else
	{
		item= XMLTools::selectSingleNode(xmlDoc,"/CurveProbabilityList/CurveProbability/CurveInfo/MaturitiesAdjustment"); 
		std::string adjrule; 
		XMLTools::convert(item,adjrule); 
		ICM_EnumsCnv::cnv(adjrule,adj); 
	}*/

	item= XMLTools::selectSingleNode(xmlDoc,"/CurveProbabilityList/CurveProbability/CurveInfo/MaturitiesAdjustment"); 
	std::string adjrule; 
	XMLTools::convert(item,adjrule); 

	if (adjrule == "NONE") 
	{
		item= XMLTools::selectSingleNode(xmlDoc,"/CurveProbabilityList/CurveProbability/Underlyings/Underlying[0]/RollMaturity"); 
		std::string adjrule; 
		XMLTools::convert(item,adjrule); 
		if (adjrule=="QTR") 	adj=qCredit_Adjust20;
		else if (adjrule=="SA") 	adj=qCredit_Adjust20SA;
		else 	adj=qCredit_Default;
	}
	else
	{
		ICM_EnumsCnv::cnv(adjrule,adj); 
	}


	
/*** 
	MSXML2::IXMLDOMNodeList * resultList = NULL;
	MSXML2::IXMLDOMNode * listItem = NULL, * theNode = NULL;
	long nbNodes;

	try
	{
		//---------------------------------------------------------------
		//Récupération des spreads et des Tenors
		//---------------------------------------------------------------
		_bstr_t* Info = new _bstr_t(DTD_GETDEFPROBCURVE);
		if (XMLDoc->selectNodes(*Info, &resultList) == S_OK)
		{
			resultList->get_length(&nbNodes);

			if (nbNodes == 0)
			{
				hr = S_FALSE;
				throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
				 "Invalid XML string for Def Prob Curve");
			}

			pPLOT = new char*[nbNodes];
			pSPREAD = new double[nbNodes];

			for (long indexNode=0 ; indexNode<nbNodes ; indexNode++)
			{
				hr=resultList->get_item(indexNode, &listItem);
				if (hr==S_OK && listItem!=NULL)
				{
					_bstr_t* bDate = new _bstr_t("Date");
					listItem->selectSingleNode(*bDate, &theNode);
					if (theNode!=NULL)
					{
						BSTR resultat = NULL;
						theNode->get_text(&resultat);

						_bstr_t ff(resultat,false);
						char * ff1=(char *)ff;

						pPLOT[indexNode] = new char[5];
						strcpy(pPLOT[indexNode],ff1); 

						theNode->Release();
						theNode=NULL;
						if (resultat) SysFreeString(resultat);
					}
					delete bDate;

					_bstr_t* bRate =  new _bstr_t("MidRate");
					listItem->selectSingleNode(*bRate,&theNode);
					if (theNode!=NULL)
					{
						BSTR resultat = NULL;
						theNode->get_text(&resultat);	//printf("resultat : %S\n", resultat);

						_bstr_t ff(resultat,false);
						char * ff1=(char *)ff;

						double rate =  atof(ff1);	//printf("rate : %.13lf\n", rate);

						pSPREAD[indexNode] = rate;

						theNode->Release();
						theNode=NULL;
						if (resultat) SysFreeString(resultat);
					}
					delete bRate;

					listItem->Release();
					listItem=NULL;
				}
			}

		}
		delete Info;

		//---------------------------------------------------------------
		//Récupération Recovery/Currency/Frequency
		//---------------------------------------------------------------

 		if (XMLDoc->selectNodes((_bstr_t)("//COMMDATA"), &resultList) == S_OK)
		{

			hr=resultList->get_item(0, &listItem);
			Recovery = GetDoubleFromXMLNode (listItem,"Value");
		}

		if (listItem)
		{
			listItem->Release();
			listItem = NULL;
		}
		if (resultList)
		{
			resultList->Release();
			resultList = NULL;
		}

 		if (XMLDoc->selectNodes((_bstr_t)("//MKTDATA"), &resultList) == S_OK)
		{
			hr=resultList->get_item(0, &listItem);
			CCString Meth = GetCCStringFromXMLNode (listItem,"Method");
			AIMMADJ = FromSummitAdjCDSToARM(Meth);
		}

		if (listItem)
		{
			listItem->Release();
			listItem = NULL;
		}
		if (resultList)
		{
			resultList->Release();
			resultList = NULL;
		}


 		if (XMLDoc->selectNodes((_bstr_t)("//MKTSPEC"), &resultList) == S_OK)
		{
			hr=resultList->get_item(0, &listItem);
			CCString Meth = GetCCStringFromXMLNode (listItem,"GenAssetId");

			if (AIMMADJ != qCredit_Default)
			{AIMMADJ = FromSummitAdjCDSToARM(Meth);
			}
		}

		if (listItem)
		{
			listItem->Release();
			listItem = NULL;
		}
		if (resultList)
		{
			resultList->Release();
			resultList = NULL;
		}



 		if (XMLDoc->selectNodes((_bstr_t)("//GASSET"), &resultList) == S_OK)
		{

			hr=resultList->get_item(0, &listItem);
	
			//Recuperation de frequence
			listItem->selectSingleNode(_bstr_t("SCHED_Pay_Freq"), &theNode);
			if (theNode!=NULL)
			{
				BSTR resultat = NULL;
				theNode->get_text(&resultat);

				_bstr_t ff(resultat,false);
				char * ff1=(char *)ff;
				Frequency = FromSummitFreqToARMFreq((const char*) ff1);

				theNode->Release();
				theNode=NULL;
				if (resultat) SysFreeString(resultat);
			}
			if (theNode)
			{
			theNode->Release();
			theNode = NULL;
			}

			//Recuperation de la ccy
			listItem->selectSingleNode(_bstr_t("Ccy"), &theNode);
			if (theNode!=NULL)
			{
				BSTR resultat = NULL;
				theNode->get_text(&resultat);

				_bstr_t ff(resultat,false);
				char * ff1=(char *)ff;
				ccy = new ARM_Currency((const char*)ff1);

				theNode->Release();
				theNode=NULL;
				if (resultat) SysFreeString(resultat);
			}
			if (theNode)
			{
			theNode->Release();
			theNode = NULL;
			}

			//Recuperation de l'index
			listItem->selectSingleNode(_bstr_t("INTEREST_UnderIndex"), &theNode);
			if (theNode!=NULL)
			{
				BSTR resultat = NULL;
				theNode->get_text(&resultat);

				_bstr_t ff(resultat,false);
				char * ff1=(char *)ff;
				IndexName = (CCString)ff1;

				theNode->Release();
				theNode=NULL;
				if (resultat) SysFreeString(resultat);
			}
			if (theNode)
			{
			theNode->Release();
			theNode = NULL;
			}

		}

		if (listItem)
		{
			listItem->Release();
			listItem = NULL;
		}

		if (theNode)
		{
			theNode->Release();
			theNode = NULL;
		}
		if (resultList)
		{
			resultList->Release();
			resultList = NULL;
		}
		ICM_RELEASEXML(XMLDoc); 

		hr = S_OK;
	}

	catch(...)
	{
		for (int j=0; j<nbNodes; j++)
		{	
			delete[] pPLOT[j];
		}

		if (pPLOT)
			delete[] pPLOT;
		if (pSPREAD)
			delete[] pSPREAD;
		
		if (listItem) listItem->Release();
		if (resultList) resultList->Release();
		if (theNode) theNode->Release();
		ICM_RELEASEXML(XMLDoc); 

		hr = S_FALSE;
		throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
         "Error in XML parsing for ZC");
	}

	ARM_Vector* mktData = NULL;

	try
	{

		for (int i=0; i<nbNodes; i++)
		{
			matu.push_back(pPLOT[i]);
			spread.push_back(pSPREAD[i]*10000.);
			recovery.push_back(Recovery);
		}

		currency = (CCString)ccy->GetCcyName();
		indexname = (CCString)IndexName;

		if (ccy)
			delete ccy;
		ccy = NULL;

		for (int j=0; j<nbNodes; j++)
			delete[] pPLOT[j];

		delete[] pPLOT;
		delete[] pSPREAD;

		if (mktData)
			delete mktData;
		mktData = NULL;

		return ARM_OK;

	}
	catch(Exception& x)
	{
		{
		if (pPLOT)
		for (int j=0; j<nbNodes; j++)
			if (pPLOT[j])
				delete[] pPLOT[j];
		}

		delete[] pPLOT;
		delete[] pSPREAD;

		if (mktData)
			delete mktData;
		mktData = NULL;

		if (ccy)
			delete ccy;
		ccy = NULL;

		hr = S_FALSE;
		char Tmp[1000];
		x.GetErrorMessage(Tmp);
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,Tmp);
	}
	catch(...)
	{
		{
		if (pPLOT)
		for (int j=0; j<nbNodes; j++)
			if (pPLOT[j])
				delete[] pPLOT[j];
		}

		delete[] pPLOT;
		delete[] pSPREAD;

		if (mktData)
			delete mktData;
		mktData = NULL;

		if (ccy)
			delete ccy;
		ccy = NULL;

		hr = S_FALSE;
		throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
         "Error in creating DefProbCurve");
	}
	**/ 
}

// ---------------------------------------------------------------------------------------------------
// Generic Tranche Parsing for Correlation 
// ---------------------------------------------------------------------------------------------------

double ICMLOCAL_ParseCorrCredit(const char* chaineXML, const CCString& bookName, const CCString& structureId)
{
	VARIANT_BOOL bOK;
	HRESULT hr;

	MSXML2::IXMLDOMDocument *XMLDoc = NULL;
	MSXML2::IXMLDOMNodeList * resultList = NULL;
	MSXML2::IXMLDOMNode * listItem = NULL, * theNode = NULL;
	MSXML2::IXMLDOMNode * floatingNode = NULL;
	long nbNodes;

	try
	{
	hr = CoInitialize(NULL); 
//	SUCCEEDED(hr) ? 0 : throw hr;
	hr = CoCreateInstance(__uuidof(MSXML2::DOMDocument30), NULL, CLSCTX_INPROC_SERVER, __uuidof(MSXML2::IXMLDOMDocument), (void**)&XMLDoc);
	SUCCEEDED(hr) ? 0 : throw hr;

	_bstr_t tmpChaine ( chaineXML );

	XMLDoc->loadXML(tmpChaine, &bOK);

	}
	catch(...)
	{
	ICM_RELEASEXML(XMLDoc); 
	hr = S_FALSE;
	throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"Pb in creating XML document for getting Tranches informations");
	}

	try
	{
	// ---------------------------------------------------------------------------------
	// Recuperation de la correlation Issuer1/Issuer2
	// ---------------------------------------------------------------------------------

	if (XMLDoc->selectNodes((_bstr_t)("Response/Entity/COMMSET"), &resultList) == S_OK)
	{
		resultList->get_length(&nbNodes);
		if (nbNodes != 1)
		{
			hr = S_FALSE;
			throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT, "Invalid XML string for getting startadte");
		}

		hr=resultList->get_item(0, &listItem);

		//JLA: resultList should be deleted here. 
		if (resultList)
		{
			resultList->Release();
			resultList = NULL;
		}

		double corr = XML_doubleNodeTreating(listItem,"Value");
		if (listItem) {
			listItem->Release(); listItem=0; }
		return (corr);
	}

	if (resultList)
	{
		resultList->Release();
		resultList = NULL;
	}
	if (listItem) { listItem->Release(); listItem=0; }
	} 

	catch(...)
	{		
		if (theNode) theNode->Release();
		if (listItem) listItem->Release();
		if (resultList) resultList->Release();
		ICM_RELEASEXML(XMLDoc); 

		hr = S_FALSE;
		throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
         "Error in XML parsing for getting Tranche");
	}

	return -999.;
}



// ---------------------------------------------------------------------------------------------------
// Generic Tranche Parsing for Model 
// ---------------------------------------------------------------------------------------------------

/** 

ARM_Model* ICMLOCAL_ParseModelCredit(ARM_ZeroCurve* DiscountCurve, 
									 const char* chaineXML, 
									 const CCString& SummitId, 
									 const CCString& Type,
									 const CCString& CurveId,
									 const CCString& CorrCurveId)
{
	VARIANT_BOOL bOK;
	HRESULT hr;

	MSXML2::IXMLDOMDocument *XMLDoc = NULL;
	MSXML2::IXMLDOMNodeList * resultList = NULL;
	MSXML2::IXMLDOMNode * listItem = NULL, * theNode = NULL;
	MSXML2::IXMLDOMNode * floatingNode = NULL;
	long nbNodes;

	char StructId[60];
	int i = 0,j = 0,k = 0;
	int nbissuers = 0;
	// char** Issuers = NULL; 
	std::vector<std::string> Issuers; 
	ARM_Vector Recovery ;
	ARM_Currency* ccy = NULL;
	double* Betas = NULL;
	// double** correlations = NULL;
	int avoidnbr = -1;
	ICM_Correlation* CorrId = NULL;

	ICM_ModelMultiCurves* modelmc = NULL;
	ICM_DefaultCurve** DefaultCurves = NULL;
	char Book[60];
	strcpy(Book,(const char*)SummitId);

	CCString xml_DefCurve;
	CCString messageList;

	try
	{
		hr = CoInitialize(NULL); 
//		SUCCEEDED(hr) ? 0 : throw hr;
		hr = CoCreateInstance(__uuidof(MSXML2::DOMDocument30), NULL, CLSCTX_INPROC_SERVER, __uuidof(MSXML2::IXMLDOMDocument), (void**)&XMLDoc);
		SUCCEEDED(hr) ? 0 : throw hr;

		// _bstr_t tmpChaine = chaineXML;
		wchar_t * xmlWCharText = NULL;
		xmlWCharText = constchar2wchar(chaineXML);

		// XMLDoc->loadXML(tmpChaine, &bOK);
		XMLDoc->loadXML((_bstr_t)xmlWCharText, &bOK);


	}
	catch(...)
	{
		ICM_RELEASEXML(XMLDoc); 
		hr = S_FALSE;
		throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"Pb in creating XML document for getting Tranches informations");
	}

	try
	{

		// ---------------------------------------------------------------------------------
		// Recuperation des infos communes à tous les deals
		// ---------------------------------------------------------------------------------

		if (XMLDoc->selectNodes((_bstr_t)("Response/EXOTIC/Env/ENV"), &resultList) == S_OK)
		{
			resultList->get_length(&nbNodes);

			if (nbNodes != 1)
			{
				hr = S_FALSE;
				throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT, "Invalid XML string for getting startadte");
			}

			hr=resultList->get_item(0, &listItem);
			
			//Recuperation de la Trade Date
			listItem->selectSingleNode(_bstr_t("StructureId"), &theNode);
			if (theNode!=NULL)
			{
				BSTR resultat = NULL;
				theNode->get_text(&resultat);

				_bstr_t ff(resultat,false);
				char * ff1=(char *)ff;
				strcpy(StructId,ff1); 

				theNode->Release();
				theNode=NULL;
				if (resultat) SysFreeString(resultat);
			}
			listItem->Release();
			listItem=NULL;
		}

		if (resultList)
		{
			resultList->Release();
			resultList = NULL;
		}

		if (XMLDoc->selectNodes((_bstr_t)("Response/EXOTIC/Assets/ASSET"), &resultList) == S_OK)
		{
			resultList->get_length(&nbNodes);

			if (nbNodes != 1)
			{
				hr = S_FALSE;
				throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT, "Invalid XML string");
			}

			hr=resultList->get_item(0, &listItem);

			listItem->selectSingleNode(_bstr_t("Ccy"), &theNode);
			if (theNode!=NULL)
			{
				BSTR resultat = NULL;
				theNode->get_text(&resultat);

				_bstr_t ff(resultat,false);
				char * ff1=(char *)ff;
				ccy = new ARM_Currency(ff1);

				theNode->Release();
				theNode=NULL;
				if (resultat) SysFreeString(resultat);
			}

			listItem->Release();
			listItem=NULL;
		}

		// ---------------------------------------------------------------------------------
		// Recuperation des infos sur les Issuers
		// ---------------------------------------------------------------------------------

		if (XMLDoc->selectNodes((_bstr_t)("Response/EXOTIC/Assets/ASSET/ProdData/CREDSWAP/CreditEntList/CREDIT"), &resultList) == S_OK)
		{
			resultList->get_length(&nbNodes);

			if (nbNodes == 0)
			{	
				hr = S_FALSE;
				throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT, "Invalid XML string for getting ");
			}

			//nbissuers = nbNodes - 1;
			nbissuers = nbNodes;
			// Issuers = new char*[nbissuers];
			Issuers.resize(nbissuers); 

			k = 0;
			for (long indexNode=0 ; indexNode<nbNodes ; indexNode++)
			{
				hr=resultList->get_item(indexNode, &listItem);

				listItem->selectSingleNode(_bstr_t("DefpIssuer"), &theNode);
				if (theNode!=NULL)
				{
				BSTR resultat = NULL;
				theNode->get_text(&resultat);

				_bstr_t ff(resultat,false);
				char * ff1=(char *)ff;

				if ((strcmp(ff1,StructId)) && (strcmp(ff1,Book)))
				{
					Issuers[k] = ff1; 
					// new char[_size_zclabel_];
					// memcpy(Issuers[k],ff1,sizeof(char)*_size_zclabel_);
					k++;
				}
				else avoidnbr = indexNode;
		
				theNode->Release();
				theNode=NULL;
				if (resultat) SysFreeString(resultat);
				}
			}

			listItem->Release();
			listItem=NULL;
		}

		if (resultList)
		{
			resultList->Release();
			resultList = NULL;
		}

		if (XMLDoc->selectNodes((_bstr_t)("Response/EXOTIC/Assets/ASSET/ProdData/CREDSWAP/CreditEntList/CREDIT"), &resultList) == S_OK)
		{
			resultList->get_length(&nbNodes);

			if (nbNodes == 0)
			{	
				hr = S_FALSE;
				throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT, "Invalid XML string for getting ");
			}

			Recovery.Resize(nbissuers);
			DefaultCurves = new ICM_DefaultCurve*[nbissuers];

			k=0;
			for (long indexNode=0 ; indexNode<nbNodes ; indexNode++)
			{
				hr=resultList->get_item(indexNode, &listItem);

				ARM_Date irdate = DiscountCurve->GetAsOfDate();

				if (avoidnbr != indexNode)
				{		
				long out = etoolkit_GetDefProbCurve(Issuers[k].c_str(),irdate,CurveId, xml_DefCurve,messageList);

				if (out == ARM_KO)  
					throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"etoolkit_GetDefProbCurve fail");

				DefaultCurves[k] = ARMLOCAL_ParseDefProbCurve(xml_DefCurve,
															  irdate,
															  DiscountCurve,
															   Issuers[k].c_str());

				Recovery[k] = XML_doubleNodeTreating(listItem,"RecRate");
				k++;
				}

			}

			listItem->Release();
			listItem=NULL;
		}


		if (resultList)
		{
			resultList->Release();
			resultList = NULL;
		}


		long retour = 0.;
		CCString xmlResponse;
		CCString messageList;
		double corr = 0;

		ICM_QMatrix<double> correlations (nbissuers,nbissuers); 
		// correlations = new double*[nbissuers];

		// for (i=0; i<nbissuers; i++)
		// {
		// 	correlations[i] = new double[nbissuers];
		// 	memset(correlations[i],'\0',sizeof(double)*nbissuers); 
		// }

		for (i=0; i<nbissuers; i++)
		{
			if (!(CorrCurveId == "NONE"))
			{
				for (j=0; j<i; j++)
				{
				retour =  etoolkit_GetCorrFromSummit (DiscountCurve->GetAsOfDate(),
															  Issuers[i].c_str(),
															  Issuers[j].c_str(),
															  CorrCurveId,
															  xmlResponse,
															  messageList);

				correlations(i,j)  = ICMLOCAL_ParseCorrCredit(xmlResponse,SummitId,Type);

				if (retour == ARM_KO)
				{
					throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
					"Error in correlations");
				}

				correlations(j,i)= correlations(i,j) ;
				}
			}

			correlations(i,i)  = 1.0;
		}

		//------------------------------------------------------------------------------
		//Création du model
		//------------------------------------------------------------------------------

		if (Type == "MMC")
		{
			CorrId = new ICM_CorrMatrix(DiscountCurve->GetAsOfDate(),"ETK",Issuers,correlations);

			modelmc = new ICM_ModelMultiCurves(nbissuers,
											   DefaultCurves,
											   DiscountCurve,
											   Recovery,
											   CorrId);

			if (CorrId) 
				delete CorrId;
			CorrId = NULL;
	
		}
		else
			throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
			     "Unknwon Model");


		for (i = 0; i<nbissuers; i++)
		{
			// if (Issuers[i]) delete[] Issuers[i];
			if (DefaultCurves[i]) delete DefaultCurves[i]; 
			// if (correlations[i]) delete correlations[i];
		}
		
		//if (correlations) delete[] correlations;
		if (DefaultCurves) delete[] DefaultCurves;
		//if (Recovery) delete[] Recovery;		
		//if (Issuers) delete[] Issuers;		

		if (CorrId) delete CorrId;
		CorrId = NULL;

		if (ccy) delete ccy;

		if (theNode)
		{
			theNode->Release();
			theNode = NULL;
		}

		if (resultList)
		{
			resultList->Release();
			resultList = NULL;
		}

		if (listItem)
		{
			listItem->Release();
			listItem = NULL;
		}

		ICM_RELEASEXML(XMLDoc); 

		return (modelmc);
	}

	catch(...)
	{		
		if (theNode) theNode->Release();
		if (listItem) listItem->Release();
		if (resultList) resultList->Release();
		ICM_RELEASEXML(XMLDoc); 

		hr = S_FALSE;
		throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
         "Error in XML parsing for getting Tranche");
	}

	return NULL;
}
**/

// --------------------------------------------------------------------------------
// Recuperation du secid associé à l'issuer --> Recovery
// --------------------------------------------------------------------------------
long  ARMLOCAL_GetSecidForRecovery (const char* chaineXML,CCString& SecId)
{

	VARIANT_BOOL bOK;
	HRESULT hr;
	MSXML2::IXMLDOMDocument *XMLDoc = NULL;

	try
	{
		hr = CoInitialize(NULL); 
//		SUCCEEDED(hr) ? 0 : throw hr;
		hr = CoCreateInstance(__uuidof(MSXML2::DOMDocument30), NULL, CLSCTX_INPROC_SERVER, __uuidof(MSXML2::IXMLDOMDocument), (void**)&XMLDoc);
		SUCCEEDED(hr) ? 0 : throw hr;

		//JLA _bstr_t* tmpChaine = new _bstr_t(chaineXML);
		_bstr_t* tmpChaine = new _bstr_t ;
		VariantTools::convert(std::string(chaineXML),*tmpChaine); 


		XMLDoc->loadXML(*tmpChaine, &bOK);

		delete tmpChaine;
	}
	catch(...)
	{
		ICM_RELEASEXML(XMLDoc); 
		hr = S_FALSE;
		throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
         "Pb in creating XML document for ParseDefProbCurveMktData");
	}

	MSXML2::IXMLDOMNodeList * resultList = NULL;
	MSXML2::IXMLDOMNode * listItem = NULL, * theNode = NULL;
	long nbNodes;

	try
	{
		_bstr_t* Info = new _bstr_t(DTD_GETDEFPROBCURVE);

		if (XMLDoc->selectNodes(*Info, &resultList) == S_OK)
		{
			resultList->get_length(&nbNodes);

			if (nbNodes == 0)
			{
				hr = S_FALSE;
				throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
				 "Invalid XML string for Def Prob Curve");
			}

			hr=resultList->get_item(0, &listItem);

			if (hr==S_OK && listItem!=NULL)
			{

			_bstr_t* bSec =  new _bstr_t("Sec");
			listItem->selectSingleNode(*bSec,&theNode);
			if (theNode!=NULL)
			{
				BSTR resultat = NULL;
				theNode->get_text(&resultat);	//printf("resultat : %S\n", resultat);

				_bstr_t ff(resultat,false);
				char * ff1=(char *)ff;
				
				SecId = (CCString) ff1;

				theNode->Release();
				theNode=NULL;
				if (resultat) SysFreeString(resultat);
			}
			delete bSec;

			}

			listItem->Release();
			listItem=NULL;
		}

		if (listItem)
		{
			listItem->Release();
			listItem = NULL;
		}

		if (theNode)
		{
			theNode->Release();
			theNode = NULL;
		}
		if (resultList)
		{
			resultList->Release();
			resultList = NULL;
		}

		return ARM_OK;
		
	}
	catch(Exception& )
	{
	
		if (listItem) listItem->Release();
		if (resultList) resultList->Release();
		if (theNode) theNode->Release();
		ICM_RELEASEXML(XMLDoc); 

		hr = S_FALSE;
		throw; 
	}
	catch(...)
	{
		if (listItem) listItem->Release();
		if (resultList) resultList->Release();
		if (theNode) theNode->Release();
		ICM_RELEASEXML(XMLDoc) ;

		hr = S_FALSE;
		throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
         "Error in XML parsing for ZC");
	}

	return ARM_KO;
}


// ---------------------------------------------------------------------------------------------------
// Generic Tranche Parsing for Model 
// ---------------------------------------------------------------------------------------------------

ARM_Model* ICMLOCAL_ParseModelForSpreadOptions(ARM_ZeroCurve* DiscountCurve, 
												const char* chaineXML)
{
	VARIANT_BOOL bOK;
	HRESULT hr;

	MSXML2::IXMLDOMDocument *XMLDoc = NULL;
	MSXML2::IXMLDOMNodeList * resultList = NULL;
	MSXML2::IXMLDOMNode * listItem = NULL, * theNode = NULL;
	MSXML2::IXMLDOMNode * floatingNode = NULL;
	long nbNodes;

	int i = 0,j = 0,k = 0;
	ARM_Currency* ccy = new ARM_Currency("EUR");
	ICM_DefaultCurveModel* model = NULL;
	ICM_DefaultCurve* DefaultCurves = NULL;
	ARM_VolCurve* VolCurve = NULL;

	long retour = 0.;
	CCString xmlResponse_Credit;
	CCString xmlResponse;
	CCString messageList;

	ARM_Date irdate;
	char Currency[5];
	char DefProbIssuer[255];
	CCString xml_DefCurve;

	try
	{
		hr = CoInitialize(NULL); 
//		SUCCEEDED(hr) ? 0 : throw hr;
		hr = CoCreateInstance(__uuidof(MSXML2::DOMDocument30), NULL, CLSCTX_INPROC_SERVER, __uuidof(MSXML2::IXMLDOMDocument), (void**)&XMLDoc);
		SUCCEEDED(hr) ? 0 : throw hr;

		// _bstr_t tmpChaine = chaineXML;
		wchar_t * xmlWCharText = NULL;
		xmlWCharText = constchar2wchar(chaineXML);

		XMLDoc->loadXML((_bstr_t)xmlWCharText, &bOK);
	}
	catch(...)
	{
		ICM_RELEASEXML(XMLDoc); 
		hr = S_FALSE;
		throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"Pb in creating XML document for getting Tranches informations");
	}

	try
	{
		irdate = DiscountCurve->GetAsOfDate();

		// ---------------------------------------------------------------------------------
		// Recuperation de la discount curve
		// ---------------------------------------------------------------------------------

		if (XMLDoc->selectNodes((_bstr_t)("//Option/OPTION"), &resultList) == S_OK)
		{
			resultList->get_length(&nbNodes);

			if (nbNodes != 1)
			{
				hr = S_FALSE;
				throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT, "Invalid XML string for getting startadte");
			}

			hr=resultList->get_item(0, &listItem);
			
			//Recuperation de la Trade Date
			listItem->selectSingleNode(_bstr_t("StkCcy"), &theNode);
			if (theNode!=NULL)
			{
				BSTR resultat = NULL;
				theNode->get_text(&resultat);

				_bstr_t ff(resultat,false);
				char * ff1=(char *)ff;
				strcpy(Currency,ff1); 

				theNode->Release();
				theNode=NULL;
				if (resultat) SysFreeString(resultat);
			}
			listItem->Release();
			listItem=NULL;
		}

		if (resultList)
		{
			resultList->Release();
			resultList = NULL;
		}

		// ---------------------------------------------------------------------------------
		// Recuperation des infos sur les Issuers
		// ---------------------------------------------------------------------------------

		if (XMLDoc->selectNodes((_bstr_t)("//CREDSWAP//CreditEntList//CREDIT"), &resultList) == S_OK)
		{
			resultList->get_length(&nbNodes);

			if (nbNodes == 0)
			{	
				hr = S_FALSE;
				throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT, "Invalid XML string for getting ");
			}

			k = 0;
			for (long indexNode=0 ; indexNode<nbNodes ; indexNode++)
			{
				hr=resultList->get_item(indexNode, &listItem);

				listItem->selectSingleNode(_bstr_t("DefpIssuer"), &theNode);
				if (theNode!=NULL)
				{
				BSTR resultat = NULL;
				theNode->get_text(&resultat);

				_bstr_t ff(resultat,false);
				char * ff1=(char *)ff;
				strcpy(DefProbIssuer,ff1);
		
				theNode->Release();
				theNode=NULL;
				if (resultat) SysFreeString(resultat);
				}
			}

			listItem->Release();
			listItem=NULL;
		}

		if (resultList)
		{
			resultList->Release();
			resultList = NULL;
		}

		long out = etoolkit_GetDefProbCurve(DefProbIssuer,irdate, SUMMIT_DEFAULT_CURVE, xmlResponse_Credit,messageList);

		DefaultCurves = ARMLOCAL_ParseDefProbCurve(xmlResponse_Credit,irdate,DiscountCurve,DefProbIssuer);


		VolCurve = etoolkit_GetVolATMFromSummit(DefProbIssuer,
												Currency,
												SUMMIT_DEFAULT_CURVE,
												irdate,
												"ISSUER",
												"ISSRVOL");


		//------------------------------------------------------------------------------
		//Création du model
		//------------------------------------------------------------------------------

		model = new ICM_DefaultCurveModel(DefaultCurves,
										  DiscountCurve,
										  VolCurve);


		if (VolCurve) delete VolCurve;
		if (DefaultCurves) delete DefaultCurves;

		if (theNode)
		{
			theNode->Release();
			theNode = NULL;
		}

		if (resultList)
		{
			resultList->Release();
			resultList = NULL;
		}

		if (listItem)
		{
			listItem->Release();
			listItem = NULL;
		}

		ICM_RELEASEXML(XMLDoc); 

		return (model);
	}

	catch(...)
	{		
		if (theNode) theNode->Release();
		if (listItem) listItem->Release();
		if (resultList) resultList->Release();
		ICM_RELEASEXML(XMLDoc); 

		hr = S_FALSE;
		throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
         "Error in XML parsing for getting Tranche");
	}

	return NULL;
}


//	----------------------------------------------------------------------------------------
void 
ICMLOCAL_ParseBasketCorrelMkDataFromCalypso(const std::string& xml,
										  std::vector<std::string>& matus,
										  std::vector<double>&attach,
										  ICM_QMatrix<double>& correls)
{
	MSXML2::IXMLDOMDocumentPtr xmlDoc; 
	xmlDoc = XMLTools::LoadXML(xml); 
	#ifdef _DEBUG
	_variant_t tmp("c:\\temp\\Calypso_basketCorrel.xml"); 
	xmlDoc->save(tmp); 
	#endif

	MSXML2::IXMLDOMNodeListPtr xmlNodeList ;
	xmlNodeList = XMLTools::selectNodes(xmlDoc,"/ListeBaseCorrel/BaseCorrel/Points/Point"); 
	long nbNodes(0); 
	xmlNodeList->get_length(&nbNodes); 
	std::set<std::string> dateSet; 
	std::set<double> attachSet; 
	for(long i=0;i<nbNodes;i++) 
	{
		std::string tmpDate; 
		double tmpAttach; 
		XMLTools::convert(XMLTools::selectSingleNode(XMLTools::get_item(xmlNodeList,i),"Date"),tmpDate); 
		XMLTools::convert(XMLTools::selectSingleNode(XMLTools::get_item(xmlNodeList,i),"Detachment"),tmpAttach); 
		dateSet.insert(tmpDate) ; 
		attachSet.insert(tmpAttach*10000.) ;
	}
	matus.resize(dateSet.size()); 
	std::copy(dateSet.begin(),dateSet.end(),matus.begin()); 
	attach.resize(attachSet.size()); 
	std::copy(attachSet.begin(),attachSet.end(),attach.begin()); 
	correls.Resize(matus.size(),attach.size()); 
	
	std::set<std::string>::const_iterator itMatus ;
	std::set<double>::const_iterator itAttach;
	long j; 
	for(i=0,itMatus =dateSet.begin(); itMatus!=dateSet.end(); ++itMatus,i++)
		for(j=0,itAttach=attachSet.begin(); itAttach!=attachSet.end(); ++itAttach,j++)
		{
			std::stringstream sstr; 
			sstr<<"/ListeBaseCorrel/BaseCorrel/Points/Point[Date='"<<*itMatus<<"' and Detachment="<<*itAttach/10000.<<"]/Value" ;
			double tmp; 
			XMLTools::convert(XMLTools::selectSingleNode(xmlDoc,sstr.str()),tmp) ; 
			correls(i,j)=tmp*100.; 
		}
	if ( attach.size() > 0 )
	{
		if (attach[0] == 1) attach[0]=0 ;
	}
}
