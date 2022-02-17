
#include <ARM\libarm_local\firstToBeIncluded.h>


#include <atlbase.h>
#include "XMLTools.h"
//#include <ARM\libarm_local\arm_local_glob.h>
#include <ARM\libarm_frometk\arm_local_parsexml_nt_curves.h>
#include <ARM\libarm_frometk\arm_local_parsexml_util.h>
#include "VariantTools.h"
/** 
#include <ARM\libarm_frometk\arm_local_parsexml_common.h>
#include <inst\irindex.h>
#include <glob\linalg.h>
#include <glob\dates.h>
#include <crv\zeroint.h>
#include <ccy\currency.h>
**/ 

#ifndef XML_DEFINE
#define XML_DEFINE
#import <msxml3.dll> raw_interfaces_only
using namespace MSXML2;
#endif


// -----------------------------------------------
// construction de la zcpy avec stripper summit
// -----------------------------------------------

ARM_ZeroCurve* ARMLOCAL_XML_ZCPY_with_summit_stripper(const char* chaineXML,const char* chaineCCY)
{
	VARIANT_BOOL bOK;
	HRESULT hr;
	MSXML2::IXMLDOMDocument* XMLDoc = NULL;

	ARM_Date AsOfDate;
	char MCType[255];
	ARM_Currency* ccy=NULL;
	int InterpolType = 0;

	double* pdRate = NULL;
	double* pdMatu = NULL;

	CCString Input = (CCString)chaineCCY;
	// CCString strtmp;
	wchar_t * xmlWCharText = NULL;
	
	try
	{
		hr = CoInitialize(NULL); 
//		SUCCEEDED(hr) ? 0 : throw hr;
		hr = CoCreateInstance(__uuidof(MSXML2::DOMDocument30), NULL, CLSCTX_INPROC_SERVER, __uuidof(MSXML2::IXMLDOMDocument), (void**)&XMLDoc);
		SUCCEEDED(hr) ? 0 : throw hr;
		
		xmlWCharText = constchar2wchar(chaineXML);

		Local_XML_TRACE("c:\\temp\\log_arm_credit.txt","Lecture XML pour ZCPY Curve");
		XMLDoc->loadXML((_bstr_t)xmlWCharText, &bOK);

	}
	catch(...)
	{
		if (XMLDoc) XMLDoc->Release();
		hr = S_FALSE;

		CCString msg ((CCString)"Pb in creating XML document for ARMLOCAL_XML_ZCPY_with_summit_stripper \n");

		Local_XML_TRACE("c:\\temp\\log_arm_credit.txt",msg);

		if (xmlWCharText)
			delete xmlWCharText;
		xmlWCharText = NULL;

		throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
						 (const char*) msg);
	}

	MSXML2::IXMLDOMNodeList * resultList = NULL;
	MSXML2::IXMLDOMNode * listItem = NULL, * theNode = NULL;
	long nbNodes;

	try
	{
		// ---------------------------------------------------------------------------------
		// Recuperation des infos de la courbe de taux
		// ---------------------------------------------------------------------------------

		Local_XML_TRACE("c:\\temp\\log_arm_credit.txt","Recuperation des infos de la courbe de taux");
		std::string indexName = GetSummitCurveIndexFromCurrency(chaineCCY); 
		std::stringstream sstr; 
		if (Input=="")
			sstr<< "//CURVEHEAD/CVCHAR[MCType = \"ZERO\"]";
		else
			sstr<<"//CURVEHEAD/CVCHAR[MCType = \"ZERO\" and Ccy = \""<<(const char*)Input<<"\" and dmIndex = \"" <<indexName <<"\"]";

 		if (XMLDoc->selectNodes((_bstr_t)sstr.str().c_str(), &resultList) == S_OK)
		{
			hr=resultList->get_item(0, &listItem);
			
			//Recuperation de la Date AsOf
			AsOfDate = GetDateFromXMLNode(listItem,"AsOfDate");

			//Recuperation du type de la courbe
			listItem->selectSingleNode((_bstr_t)"MCType", &theNode);
			if (theNode!=NULL)
			{
				BSTR resultat = NULL;
				theNode->get_text(&resultat);

				_bstr_t ff(resultat,false);
				char * ff1=(char *)ff;
				strcpy(MCType,ff1);

				theNode->Release();
				theNode=NULL;
				if (resultat) SysFreeString(resultat);
			}

			//Recuperation de la devise
			listItem->selectSingleNode((_bstr_t)"Ccy", &theNode);
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

			//Recuperation de la méthode d'interpolation
			listItem->selectSingleNode((_bstr_t)"InterpMethod", &theNode);
			if (theNode!=NULL)
			{
				BSTR resultat = NULL;
				theNode->get_text(&resultat);

				_bstr_t ff(resultat,false);
				char * ff1=(char *)ff;
				InterpolType = Summit_InterpMethod((const char*) ff1);

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
		// Recuperation des Dates & des rates de la courbe de taux
		// ---------------------------------------------------------------------------------

		Local_XML_TRACE("c:\\temp\\log_arm_credit.txt","Recuperation des Dates & des rates de la courbe de taux");

		std::stringstream sstr2 ;
		
		if (Input=="")
			sstr2<<"//CURVEHEAD/ENTITYLIST/CURVE[../../CVCHAR/MCType = \"ZERO\"]";
		else
			sstr2<<"//CURVEHEAD/ENTITYLIST/CURVE[../../CVCHAR/MCType = \"ZERO\" and ../../CVCHAR/Ccy = \"" 
				<<(const char*)Input << "\" and ../../CVCHAR/dmIndex = \"" << indexName  << "\"]";

		if (XMLDoc->selectNodes((_bstr_t)sstr2.str().c_str(), &resultList) == S_OK)
		{
			resultList->get_length(&nbNodes);

			if (nbNodes == 0)
			{
				hr = S_FALSE;

				CCString msg ((CCString)"Invalid XML string for ZC \n" + chaineXML);
		
				throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
								 (const char*) msg);
			}

			pdRate = new double[nbNodes];
			pdMatu = new double[nbNodes];

			for (long indexNode=0 ; indexNode<nbNodes ; indexNode++)
			{
				hr=resultList->get_item(indexNode, &listItem);
				if (hr==S_OK && listItem!=NULL)
				{
					listItem->selectSingleNode((_bstr_t)"Date", &theNode);
					if (theNode!=NULL)
					{
						BSTR resultat = NULL;
						theNode->get_text(&resultat);

						int y;
						int m;
						int d;

						_bstr_t ff(resultat,false);
						char * ff1=(char *)ff;

						sscanf(ff1, "%04d%02d%02d", &y, &m, &d);

						ARM_Date endDate(d, m, y);
						double maturite = (endDate - AsOfDate) / 365.;

						pdMatu[indexNode] = maturite;
						theNode->Release();
						theNode=NULL;
						if (resultat) SysFreeString(resultat);
					}

					listItem->selectSingleNode((_bstr_t)"Rate",&theNode);
					if (theNode!=NULL)
					{
						BSTR resultat = NULL;
						theNode->get_text(&resultat);	//printf("resultat : %S\n", resultat);

						_bstr_t ff(resultat,false);
						char * ff1=(char *)ff;

						double rate = 100. * atof(ff1);	//printf("rate : %.13lf\n", rate);
						pdRate[indexNode] = rate;
						theNode->Release();
						theNode=NULL;
						if (resultat) SysFreeString(resultat);
					}

					listItem->Release();
					listItem=NULL;
				}
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

		if (XMLDoc) XMLDoc->Release();
		XMLDoc = NULL;

	}

	catch(Exception& )
	{
		if (pdRate)
			delete[] pdRate;
		if (pdMatu)
			delete[] pdMatu;
		
		if (listItem) listItem->Release();
		if (resultList) resultList->Release();
		if (theNode) theNode->Release();
		if (XMLDoc) XMLDoc->Release();

		XMLDoc = NULL;

		hr = S_FALSE;

		CCString msg ((CCString)"Error in XML parsing for ARMLOCAL_XML_ZCPY_with_summit_stripper \n");
		Local_XML_TRACE("c:\\temp\\log_arm_credit.txt",msg);

		if (xmlWCharText)
			delete xmlWCharText;
		xmlWCharText = NULL;

		throw; 
	}
	catch(...)
	{
		if (pdRate)
			delete[] pdRate;
		if (pdMatu)
			delete[] pdMatu;
		
		if (listItem) listItem->Release();
		if (resultList) resultList->Release();
		if (theNode) theNode->Release();
		if (XMLDoc) XMLDoc->Release();

		XMLDoc = NULL;

		hr = S_FALSE;

		CCString msg ((CCString)"Error in XML parsing for ARMLOCAL_XML_ZCPY_with_summit_stripper \n");

		Local_XML_TRACE("c:\\temp\\log_arm_credit.txt",msg);

		if (xmlWCharText)
			delete xmlWCharText;
		xmlWCharText = NULL;

		throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
						 (const char*) msg);
	}

	ARM_Vector* mktData = NULL;
	ARM_Vector* mktMatu = NULL;

	ARM_Currency* aCcy = NULL;

	ARM_ZeroLInterpol* newZcLin = NULL;

	try
	{
		mktData = new ARM_Vector(nbNodes,pdRate);
		mktMatu = new ARM_Vector(nbNodes,pdMatu);

		Local_XML_TRACE("c:\\temp\\log_arm_credit.txt","Création de la courbe de taux");

		newZcLin = new ARM_ZeroLInterpol(AsOfDate, mktMatu, mktData, 0, 0, InterpolType);
		newZcLin->SetCurrencyUnit(ccy);

		Local_XML_TRACE("c:\\temp\\log_arm_credit.txt","Fin création de la courbe de taux");

		if (ccy)
			delete ccy;
		ccy = NULL;

		delete[] pdRate;
		delete[] pdMatu;

		if (mktData)
			delete mktData;
		mktData = NULL;

		if (mktMatu)
			delete mktMatu;
		mktMatu = NULL;

		if (xmlWCharText)
			delete xmlWCharText;
		xmlWCharText = NULL;

		return (ARM_ZeroCurve*) newZcLin;
	}
	catch(Exception& )
	{
		delete[] pdRate;
		delete[] pdMatu;

		if (mktData)
			delete mktData;
		mktData = NULL;

		if (mktMatu)
			delete mktMatu;
		mktMatu = NULL;

		if (aCcy)
			delete aCcy;
		aCcy = NULL;

		if (xmlWCharText)
			delete xmlWCharText;
		xmlWCharText = NULL;

		hr = S_FALSE;

		Local_XML_TRACE("c:\\temp\\log_arm_credit.txt","Unknown Error in  ZCPY construction");
		throw; 
	}
	catch(...)
	{
		delete[] pdRate;
		delete[] pdMatu;

		if (mktData)
			delete mktData;
		mktData = NULL;

		if (mktMatu)
			delete mktMatu;
		mktMatu = NULL;

		if (aCcy)
			delete aCcy;
		aCcy = NULL;

		if (xmlWCharText)
			delete xmlWCharText;
		xmlWCharText = NULL;

		hr = S_FALSE;

		Local_XML_TRACE("c:\\temp\\log_arm_credit.txt","Unknown Error in  ZCPY construction");

		throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
         "Error in  ZCPY construction");
	}
}


// -----------------------------------------------
// construction de la defprob avec stripper summit
// -----------------------------------------------

ICM_DefaultCurve* ARMLOCAL_XML_DEFPROB_with_summit_stripper(const char* chaineXML,
														    const char* IssuerName, 
														    ARM_ZeroCurve* ircurve_input)
{
	VARIANT_BOOL bOK;
	HRESULT hr;
	MSXML2::IXMLDOMDocument2 *XMLDoc = NULL;

	ARM_Date AsOfDate;
	char MCType[255];
	ARM_Currency* ccy=NULL;
	int InterpolType = 0;
	char IssuerLabel[255];

	ARM_ZeroCurve* ircurve = NULL;
	double* pYearFractions = NULL;
	ARM_Vector* pIntensity = NULL;
	ARM_Vector* pRecovery = NULL;

	wchar_t * xmlWCharText = NULL;

	CCString strtmp;
	CCString Input = (CCString)IssuerName;
	CCString Titre;

	double Fixed_Recovery_Rate = CREDIT_DEFAULT_VALUE;
	double Fixed_Recovery_Rate_Index = 0.4;
	double Fixed_Junior_Recovery = 0.2;

	try
	{
		hr = CoInitialize(NULL); 
//		SUCCEEDED(hr) ? 0 : throw hr;
		hr = CoCreateInstance(__uuidof(MSXML2::DOMDocument30), NULL, CLSCTX_INPROC_SERVER, __uuidof(MSXML2::IXMLDOMDocument), (void**)&XMLDoc);
		SUCCEEDED(hr) ? 0 : throw hr;

		xmlWCharText = constchar2wchar(chaineXML);

		Local_XML_TRACE("c:\\temp\\log_arm_credit.txt","Lecture XML pour DefProb Curve");
		XMLDoc->loadXML((_bstr_t)xmlWCharText, &bOK);
	}
	catch(...)
	{
		if (XMLDoc) XMLDoc->Release();
		hr = S_FALSE;

		CCString msg ((CCString)"Pb in creating XML document for DefProb Curve : ARMLOCAL_XML_DEFPROB_with_summit_stripper \n");

		Local_XML_TRACE("c:\\temp\\log_arm_credit.txt",msg);

		if (xmlWCharText)
			delete xmlWCharText;
		xmlWCharText = NULL;

		throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
						 (const char*) msg);
	}

	MSXML2::IXMLDOMNodeList * resultList = NULL;
	MSXML2::IXMLDOMNode * listItem = NULL, * theNode = NULL;
	long nbNodes;

	try
	{
		
		// ---------------------------------------------------------------------------------
		// Création de la courbe de défault
		// ---------------------------------------------------------------------------------

		if (Input == "")
		 strtmp = "//CURVEHEAD/CVCHAR[MCType = \"DEFPROB\"]";
		else
		{	
		 strtmp = "//CURVEHEAD/CVCHAR[MCType = \"DEFPROB\" and Issuer = \"";	
		 strtmp = strtmp + Input;
		 strtmp	= strtmp + "\"]";
		}

 		if (XMLDoc->selectNodes((_bstr_t)(const char*)strtmp, &resultList) == S_OK)
		{

			hr=resultList->get_item(0, &listItem);
			
			//Recuperation de la Date AsOf
			AsOfDate = GetDateFromXMLNode(listItem,"AsOfDate");

			//Recuperation du type de la courbe
			listItem->selectSingleNode((_bstr_t)"MCType", &theNode);
			if (theNode!=NULL)
			{
				BSTR resultat = NULL;
				theNode->get_text(&resultat);

				_bstr_t ff(resultat,false);
				char * ff1=(char *)ff;
				strcpy(MCType,ff1);

				theNode->Release();
				theNode=NULL;
				if (resultat) SysFreeString(resultat);
			}

			//Recuperation de la méthode d'interpolation
			listItem->selectSingleNode((_bstr_t)"InterpMethod", &theNode);
			if (theNode!=NULL)
			{
				BSTR resultat = NULL;
				theNode->get_text(&resultat);

				_bstr_t ff(resultat,false);
				char * ff1=(char *)ff;
				InterpolType = Summit_InterpMethod((const char*) ff1);

				theNode->Release();
				theNode=NULL;
				if (resultat) SysFreeString(resultat);
			}

			//Recuperation du label de l'issuer
			listItem->selectSingleNode((_bstr_t)"Issuer", &theNode);
			if (theNode!=NULL)
			{
				BSTR resultat = NULL;
				theNode->get_text(&resultat);

				_bstr_t ff(resultat,false);
				char * ff1=(char *)ff;
				strcpy(IssuerLabel,ff1);

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
		// Création de la courbe de taux
		// ---------------------------------------------------------------------------------

		Local_XML_TRACE("c:\\temp\\log_arm_credit.txt","Création de la courbe de taux dans la defprob");

		if (ircurve_input == NULL)
			ircurve = ARMLOCAL_XML_ZCPY_with_summit_stripper(chaineXML,"");
		else
			ircurve = (ARM_ZeroCurve*) ircurve_input->Clone();

		// ---------------------------------------------------------------------------------
		// Recuperation des Dates & des rates de la courbe de default
		// ---------------------------------------------------------------------------------

		Local_XML_TRACE("c:\\temp\\log_arm_credit.txt","Recuperation des Dates & des rates de la courbe de default");

		if (Input == "")
	 	 strtmp = "//CURVEHEAD/ENTITYLIST/CURVE[../../CVCHAR/MCType = \"DEFPROB\"]";
		else
		{	
		 strtmp = strtmp = "//CURVEHEAD/ENTITYLIST/CURVE[../../CVCHAR/MCType = \"DEFPROB\" and ../../CVCHAR/Issuer = \"";	
		 strtmp = strtmp + Input;
		 strtmp	= strtmp + "\"]";
		}

		if (XMLDoc->selectNodes((_bstr_t)(const char*)strtmp, &resultList) == S_OK)
		{
			resultList->get_length(&nbNodes);

			if (nbNodes == 0)
			{
				hr = S_FALSE;

				CCString msg ((CCString)"Invalid XML string for ZC \n" + chaineXML);
		
				throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
								 (const char*) msg);
			}

			pYearFractions = new double[nbNodes];
			pIntensity = new ARM_Vector(nbNodes,0.);
			pRecovery = new ARM_Vector(nbNodes,0.);

			for (long indexNode=0 ; indexNode<nbNodes ; indexNode++)
			{
				hr=resultList->get_item(indexNode, &listItem);
				if (hr==S_OK && listItem!=NULL)
				{

					listItem->selectSingleNode((_bstr_t)"Date", &theNode);
					if (theNode!=NULL)
					{
						BSTR resultat = NULL;
						theNode->get_text(&resultat);

						int y;
						int m;
						int d;

						_bstr_t ff(resultat,false);
						char * ff1=(char *)ff;

						sscanf(ff1, "%04d%02d%02d", &y, &m, &d);

						ARM_Date endDate(d, m, y);

						double YearFraction = (endDate.GetJulian() - AsOfDate.GetJulian())/K_YEAR_LEN;
						pYearFractions[indexNode] = YearFraction;

						theNode->Release();
						theNode=NULL;
						if (resultat) SysFreeString(resultat);
					}

					listItem->selectSingleNode((_bstr_t)"Rate",&theNode);
					if (theNode!=NULL)
					{
						BSTR resultat = NULL;
						theNode->get_text(&resultat);	//printf("resultat : %S\n", resultat);

						_bstr_t ff(resultat,false);
						char * ff1=(char *)ff;

						double rate = atof(ff1);	//printf("rate : %.13lf\n", rate);
						pIntensity->InitElt(indexNode,rate);
						theNode->Release();
						theNode=NULL;
						if (resultat) SysFreeString(resultat);
					}

					listItem->Release();
					listItem=NULL;
				}
			}
		}

		if (resultList)
		{
			resultList->Release();
			resultList = NULL;
		}

		// -------------------------------------------------------------
		// Recuperation du titre
		// -------------------------------------------------------------

		strtmp = (CCString)"//CREDSWAP/CreditEntList/CREDIT[DefpIssuer = \"" +
			(CCString) IssuerName + (CCString) "\"]";

		if (XMLDoc->selectNodes((_bstr_t)(const char*)strtmp, &resultList) == S_OK)
		{
			resultList->get_length(&nbNodes);


			if (nbNodes != 0)
			{
			hr=resultList->get_item(0, &listItem);

			listItem->selectSingleNode((_bstr_t)"SecID",&theNode);
			if (theNode!=NULL)
			{
				BSTR resultat = NULL;
				theNode->get_text(&resultat);	//printf("resultat : %S\n", resultat);

				_bstr_t ff(resultat,false);

				Titre = (CCString) ((char *)ff);

				theNode->Release();
				theNode=NULL;
				if (resultat) SysFreeString(resultat);
			}

			listItem->Release();
			listItem=NULL;
			}
			else
				Titre = "NONE";

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

		// -------------------------------------------------------------
		// Recuperation des recovery rates
		// -------------------------------------------------------------

		strtmp = (CCString)"//CURVEHEAD/ENTITYLIST/CURVE[../../CVCHAR/MCType = \"RECRATES\" and ../../COM1 = \"RECRATES/BOND/" +
			Titre + (CCString) "\"]";

		if (XMLDoc->selectNodes((_bstr_t)(const char*)strtmp, &resultList) == S_OK)
		{
			resultList->get_length(&nbNodes);

			if (nbNodes != 0)
			{
			hr=resultList->get_item(0, &listItem);

			Fixed_Recovery_Rate = XML_doubleNodeTreating(listItem,"Rate");

 			if (listItem) listItem->Release();
			listItem=NULL;
			}
			
			if (Titre == "NONE") Fixed_Recovery_Rate = Fixed_Recovery_Rate_Index;
			else
			{
				CCString DefPValid = IssuerName;
				DefPValid.toUpper();
				if ((DefPValid[DefPValid.GetLen()-2] == '_') && (DefPValid[DefPValid.GetLen()-1] == 'J') && (Fixed_Recovery_Rate != Fixed_Junior_Recovery))
				{
					Fixed_Recovery_Rate = Fixed_Junior_Recovery;
					//CCString msg ((CCString)"Incompatible Recovery with RefSecurity :" + Titre + " and DefProbIssuer :" + DefPValid);
					//throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,(const char*) msg);
				}
			}

			for (int i=0; i<pRecovery->GetSize(); i++)
				pRecovery->InitElt(i,Fixed_Recovery_Rate);

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

		// -------------------------------------------------------------
		// Recuperation des recovery rates (Ftd Case)
		// -------------------------------------------------------------

		if (Fixed_Recovery_Rate == CREDIT_DEFAULT_VALUE)
		{

		strtmp = "//CURVEHEAD/ENTITYLIST/CURVE[../../CVCHAR/MCType = \"RECRATES\"]";

		if (XMLDoc->selectNodes((_bstr_t)(const char*)strtmp, &resultList) == S_OK)
		{
			resultList->get_length(&nbNodes);

			if (nbNodes == 0)
			{
				hr = S_FALSE;

				CCString msg ((CCString)"Invalid XML for recovery rates \n" + chaineXML);
		
				throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
								 (const char*) msg);
			}

			hr=resultList->get_item(0, &listItem);

			Fixed_Recovery_Rate = XML_doubleNodeTreating(listItem,"Rate");

			for (int i=0; i<pRecovery->GetSize(); i++)
				pRecovery->InitElt(i,Fixed_Recovery_Rate);

			listItem->Release();
			listItem=NULL;
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

		if (XMLDoc)	{ XMLDoc->Release();}
		XMLDoc = NULL;

		if (xmlWCharText)
			delete xmlWCharText;
		xmlWCharText = NULL;

	}
	catch(Exception& )
	{
		if (pYearFractions) 
			delete pYearFractions;
		pYearFractions = NULL;

		if (pIntensity) 
			delete pIntensity;
		pIntensity = NULL;

		if (pRecovery) 
			delete pRecovery;
		pRecovery = NULL;

		if (ircurve)
			delete ircurve;
		ircurve = NULL;

		if (listItem) listItem->Release();
		if (resultList) resultList->Release();
		if (theNode) theNode->Release();
		if (XMLDoc) XMLDoc->Release();

		hr = S_FALSE;

		CCString msg ((CCString)"Unknown Error for parsing DefProb Curve XML : ARMLOCAL_XML_DEFPROB_with_summit_stripper \n");

		Local_XML_TRACE("c:\\temp\\log_arm_credit.txt",msg);

		if (xmlWCharText)
			delete xmlWCharText;
		xmlWCharText = NULL;

		throw; 
	}
	catch(...)
	{
		if (pYearFractions) 
			delete pYearFractions;
		pYearFractions = NULL;

		if (pIntensity) 
			delete pIntensity;
		pIntensity = NULL;

		if (pRecovery) 
			delete pRecovery;
		pRecovery = NULL;

		if (ircurve)
			delete ircurve;
		ircurve = NULL;
		
		if (listItem) listItem->Release();
		if (resultList) resultList->Release();
		if (theNode) theNode->Release();
		if (XMLDoc) XMLDoc->Release();

		hr = S_FALSE;

		CCString msg ((CCString)"Unknown Error for parsing DefProb Curve XML : ARMLOCAL_XML_DEFPROB_with_summit_stripper \n");

		Local_XML_TRACE("c:\\temp\\log_arm_credit.txt",msg);

		if (xmlWCharText)
			delete xmlWCharText;
		xmlWCharText = NULL;

		throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
						 (const char*) msg);
	}

	ARM_Vector* mktData = NULL;
	ARM_Vector* mktMatu = NULL;

	// ARM_Currency* aCcy = NULL;

	try
	{

		Local_XML_TRACE("c:\\temp\\log_arm_credit.txt","Création de la courbe de taux");

		ICM_DefaultCurve* newDefProb = new ICM_Constant_Piecewise(AsOfDate,
																  pYearFractions,
																  pIntensity,
																  1,
																  pRecovery->Elt(0),
																  ircurve,
																  K_ADJUSTED, // intRule
																  K_ADJUSTED, // adjStartDate
																  qCredit_Adjust20,
																  ircurve->GetCurrencyUnit()->GetCcyName(),
																  IssuerLabel,
																  NULL,qDEFCURVE_DICHO,"STD", 
																  ARM_Currency(ircurve->GetCurrencyUnit()->GetCcyName()).GetCreditStartDateLag()
																  );


		Local_XML_TRACE("c:\\temp\\log_arm_credit.txt","Fin création de la courbe de taux");

		if (ccy)
			delete ccy;
		ccy = NULL;

		if (pYearFractions) 
			delete pYearFractions;
		pYearFractions = NULL;

		if (pIntensity) 
			delete pIntensity;
		pIntensity = NULL;

		if (pRecovery) 
			delete pRecovery;
		pRecovery = NULL;

		if (ircurve)
			delete ircurve;
		ircurve = NULL;

		return newDefProb;
	}
	catch(Exception& )
	{
		if (pYearFractions) 
			delete pYearFractions;
		pYearFractions = NULL;

		if (pIntensity) 
			delete pIntensity;
		pIntensity = NULL;

		if (pRecovery) 
			delete pRecovery;
		pRecovery = NULL;

		if (ccy)
			delete ccy;
		ccy = NULL;

		if (ircurve)
			delete ircurve;
		ircurve = NULL;

		if (xmlWCharText)
			delete xmlWCharText;
		xmlWCharText = NULL;

		hr = S_FALSE;

		Local_XML_TRACE("c:\\temp\\log_arm_credit.txt","Unknown Error in  DefProbCurve object construction");

		throw; 
	}
	catch(...)
	{
		if (pYearFractions) 
			delete pYearFractions;
		pYearFractions = NULL;

		if (pIntensity) 
			delete pIntensity;
		pIntensity = NULL;

		if (pRecovery) 
			delete pRecovery;
		pRecovery = NULL;

		if (ccy)
			delete ccy;
		ccy = NULL;

		if (ircurve)
			delete ircurve;
		ircurve = NULL;

		hr = S_FALSE;

		Local_XML_TRACE("c:\\temp\\log_arm_credit.txt","Unknown Error in  DefProbCurve object construction");

		if (xmlWCharText)
			delete xmlWCharText;
		xmlWCharText = NULL;

		throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
         "Error in  ZCPY construction");
	}
}



// -------------------------------------------------------------
// construction de la defprob avec stripper summit via E-toolkit
// -------------------------------------------------------------

ICM_DefaultCurve* ARMLOCAL_XML_DEFPROB_with_summit_stripper_Etk (const char* chaineXML,
																 ARM_ZeroCurve* ircurve_input,
																 CCString OverLabel)
{
	VARIANT_BOOL bOK;
	HRESULT hr;
	MSXML2::IXMLDOMDocument2 *XMLDoc = NULL;

	ARM_Date AsOfDate;
	char Curen[4];
	ARM_Currency* ccy=NULL;
	char IssuerLabel[255];
	double recovery = 0.;

	ARM_ZeroCurve* ircurve = NULL;
	double* pYearFractions = NULL;
	ARM_Vector* pIntensity = NULL;
	ARM_Vector* pRecovery = NULL;

	wchar_t * xmlWCharText = NULL;

	CCString strtmp;

	try
	{
		hr = CoInitialize(NULL); 
//		SUCCEEDED(hr) ? 0 : throw hr;
		hr = CoCreateInstance(__uuidof(MSXML2::DOMDocument30), NULL, CLSCTX_INPROC_SERVER, __uuidof(MSXML2::IXMLDOMDocument), (void**)&XMLDoc);
		SUCCEEDED(hr) ? 0 : throw hr;

		xmlWCharText = constchar2wchar(chaineXML);

		Local_XML_TRACE("c:\\temp\\log_arm_credit.txt","Lecture XML pour DefProb Curve");
		XMLDoc->loadXML((_bstr_t)xmlWCharText, &bOK);

/*
		_variant_t fichier=(_bstr_t)"c:\\test\\curve_summit_strip.xml";
		XMLDoc->save(fichier);
*/

	}
	catch(...)
	{
		if (XMLDoc) XMLDoc->Release();
		hr = S_FALSE;

		CCString msg ((CCString)"Pb in creating XML document for DefProb Curve : ARMLOCAL_XML_DEFPROB_with_summit_stripper \n");

		Local_XML_TRACE("c:\\temp\\log_arm_credit.txt",msg);

		if (xmlWCharText)
			delete xmlWCharText;
		xmlWCharText = NULL;

		throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
						 (const char*) msg);
	}

	MSXML2::IXMLDOMNodeList * resultList = NULL;
	MSXML2::IXMLDOMNode * listItem = NULL, * theNode = NULL;
	long nbNodes;

	try
	{
		
		// ---------------------------------------------------------------------------------
		// Création de la courbe de défault
		// ---------------------------------------------------------------------------------

		 strtmp = "//CURVEHEAD";

 		if (XMLDoc->selectNodes((_bstr_t)(const char*)strtmp, &resultList) == S_OK)
		{

			hr=resultList->get_item(0, &listItem);
			
			//Recuperation de la Date AsOf
			AsOfDate = GetDateFromXMLNode(listItem,"AsOfDate");

			//Recuperation du label de l'issuer
			listItem->selectSingleNode((_bstr_t)"Issuer", &theNode);
			if (theNode!=NULL)
			{
				BSTR resultat = NULL;
				theNode->get_text(&resultat);

				_bstr_t ff(resultat,false);
				char * ff1=(char *)ff;
				strcpy(IssuerLabel,ff1);

				theNode->Release();
				theNode=NULL;
				if (resultat) SysFreeString(resultat);
			}

			//Recuperation du label de l'issuer
			listItem->selectSingleNode((_bstr_t)"Ccy", &theNode);
			if (theNode!=NULL)
			{
				BSTR resultat = NULL;
				theNode->get_text(&resultat);

				_bstr_t ff(resultat,false);
				char * ff1=(char *)ff;
				strcpy(Curen,ff1);

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
		// Création de la courbe de taux
		// ---------------------------------------------------------------------------------

		ircurve = (ARM_ZeroCurve*) ircurve_input->Clone();

		// ---------------------------------------------------------------------------------
		// Recuperation des Dates & des rates de la courbe de default
		// ---------------------------------------------------------------------------------

		Local_XML_TRACE("c:\\temp\\log_arm_credit.txt","Recuperation des Dates & des rates de la courbe de default");

	 	strtmp = "//CURVEHEAD/CurvePts/CURVE";

		if (XMLDoc->selectNodes((_bstr_t)(const char*)strtmp, &resultList) == S_OK)
		{
			resultList->get_length(&nbNodes);

			if (nbNodes == 0)
			{
				hr = S_FALSE;

				CCString msg ((CCString)"Invalid XML string for ZC \n" + chaineXML);
		
				throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
								 (const char*) msg);
			}

			pYearFractions = new double[nbNodes];
			pIntensity = new ARM_Vector(nbNodes,0.);
			pRecovery = new ARM_Vector(nbNodes,0.);

			for (long indexNode=0 ; indexNode<nbNodes ; indexNode++)
			{
				hr=resultList->get_item(indexNode, &listItem);
				if (hr==S_OK && listItem!=NULL)
				{

					listItem->selectSingleNode((_bstr_t)"Date", &theNode);
					if (theNode!=NULL)
					{
						BSTR resultat = NULL;
						theNode->get_text(&resultat);

						int y;
						int m;
						int d;

						_bstr_t ff(resultat,false);
						char * ff1=(char *)ff;

						sscanf(ff1, "%04d%02d%02d", &y, &m, &d);

						ARM_Date endDate(d, m, y);
						//double YearFraction = CountYears(KACTUAL_360,AsOfDate,endDate);
						double YearFraction = (endDate.GetJulian() - AsOfDate.GetJulian())/K_YEAR_LEN;

						pYearFractions[indexNode]=YearFraction;
						theNode->Release();
						theNode=NULL;
						if (resultat) SysFreeString(resultat);
					}

					listItem->selectSingleNode((_bstr_t)"Rate",&theNode);
					if (theNode!=NULL)
					{
						BSTR resultat = NULL;
						theNode->get_text(&resultat);	//printf("resultat : %S\n", resultat);

						_bstr_t ff(resultat,false);
						char * ff1=(char *)ff;

						double rate = atof(ff1);	//printf("rate : %.13lf\n", rate);
						pIntensity->InitElt(indexNode,rate);
						theNode->Release();
						theNode=NULL;
						if (resultat) SysFreeString(resultat);
					}

					listItem->Release();
					listItem=NULL;
				}
			}
		}

		if (resultList)
		{
			resultList->Release();
			resultList = NULL;
		}

		// -------------------------------------------------------------
		// Recuperation des recovery rates
		// -------------------------------------------------------------

		CCString XML_output;
		CCString XML_messagelist;
		long retour = 0;
		CCString SecurityId;

		retour = etoolkit_GetDefProbCurve(IssuerLabel,AsOfDate,"MO",XML_output,XML_messagelist);
		retour = ARMLOCAL_GetSecidForRecovery(XML_output,SecurityId);
		
		retour =etoolkit_GetRecoveryCurve(SecurityId,AsOfDate,"MO",XML_output,XML_messagelist);
		ARM_Vector* vout = ARMLOCAL_RecoveryCurve((const char*)XML_output);
		recovery = vout->Elt(0);
		if (vout) 
			delete vout;
		vout = NULL;
		

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

		if (XMLDoc)	{ XMLDoc->Release();}
		XMLDoc = NULL;

		if (xmlWCharText)
			delete xmlWCharText;
		xmlWCharText = NULL;

	}
	catch(Exception& )
	{
		if (pYearFractions) 
			delete[] pYearFractions;
		pYearFractions = NULL;

		if (pIntensity) 
			delete pIntensity;
		pIntensity = NULL;

		if (pRecovery) 
			delete pRecovery;
		pRecovery = NULL;

		if (ircurve)
			delete ircurve;
		ircurve = NULL;

		if (listItem) listItem->Release();
		if (resultList) resultList->Release();
		if (theNode) theNode->Release();
		if (XMLDoc) XMLDoc->Release();

		hr = S_FALSE;

		CCString msg ((CCString)"Unknown Error for parsing DefProb Curve XML : ARMLOCAL_XML_DEFPROB_with_summit_stripper \n");

		Local_XML_TRACE("c:\\temp\\log_arm_credit.txt",msg);

		if (xmlWCharText)
			delete xmlWCharText;
		xmlWCharText = NULL;

		throw; 
	}
	catch(...)
	{
		if (pYearFractions) 
			delete[] pYearFractions;
		pYearFractions = NULL;

		if (pIntensity) 
			delete pIntensity;
		pIntensity = NULL;

		if (pRecovery) 
			delete pRecovery;
		pRecovery = NULL;

		if (ircurve)
			delete ircurve;
		ircurve = NULL;
		
		if (listItem) listItem->Release();
		if (resultList) resultList->Release();
		if (theNode) theNode->Release();
		if (XMLDoc) XMLDoc->Release();

		hr = S_FALSE;

		CCString msg ((CCString)"Unknown Error for parsing DefProb Curve XML : ARMLOCAL_XML_DEFPROB_with_summit_stripper \n");

		Local_XML_TRACE("c:\\temp\\log_arm_credit.txt",msg);

		if (xmlWCharText)
			delete xmlWCharText;
		xmlWCharText = NULL;

		throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
						 (const char*) msg);
	}

	ARM_Vector* mktData = NULL;
	ARM_Vector* mktMatu = NULL;

	// ARM_Currency* aCcy = NULL;

	try
	{

		Local_XML_TRACE("c:\\temp\\log_arm_credit.txt","Création de la courbe de taux");

		if (!(OverLabel=="")) strcpy(IssuerLabel,OverLabel.GetStr());

		ICM_DefaultCurve* newDefProb = new ICM_Constant_Piecewise(AsOfDate,
																  pYearFractions,
																  pIntensity,
																  1,
																  recovery,
																  ircurve,
																  K_ADJUSTED,	// intRule
																  K_ADJUSTED,	// adjStartDate
																  qCredit_Adjust20,
																  ircurve->GetCurrencyUnit()->GetCcyName(),
															      IssuerLabel,
																  NULL,qDEFCURVE_DICHO,"STD",
																  ARM_Currency(ircurve->GetCurrencyUnit()->GetCcyName()).GetCreditStartDateLag()
																  );


		Local_XML_TRACE("c:\\temp\\log_arm_credit.txt","Fin création de la courbe de taux");

		if (ccy)
			delete ccy;
		ccy = NULL;

		if (pYearFractions) 
			delete[] pYearFractions;
		pYearFractions = NULL;

		if (pIntensity) 
			delete pIntensity;
		pIntensity = NULL;

		if (pRecovery) 
			delete pRecovery;
		pRecovery = NULL;

		if (ircurve)
			delete ircurve;
		ircurve = NULL;

		return newDefProb;
	}
	catch(Exception& )
	{
		if (pYearFractions) 
			delete[] pYearFractions;
		pYearFractions = NULL;

		if (pIntensity) 
			delete pIntensity;
		pIntensity = NULL;

		if (pRecovery) 
			delete pRecovery;
		pRecovery = NULL;

		if (ccy)
			delete ccy;
		ccy = NULL;

		if (ircurve)
			delete ircurve;
		ircurve = NULL;

		if (xmlWCharText)
			delete xmlWCharText;
		xmlWCharText = NULL;

		hr = S_FALSE;

		Local_XML_TRACE("c:\\temp\\log_arm_credit.txt","Unknown Error in  DefProbCurve object construction");

		throw ; 
	}
	catch(...)
	{
		if (pYearFractions) 
			delete[] pYearFractions;
		pYearFractions = NULL;

		if (pIntensity) 
			delete pIntensity;
		pIntensity = NULL;

		if (pRecovery) 
			delete pRecovery;
		pRecovery = NULL;

		if (ccy)
			delete ccy;
		ccy = NULL;

		if (ircurve)
			delete ircurve;
		ircurve = NULL;

		hr = S_FALSE;

		Local_XML_TRACE("c:\\temp\\log_arm_credit.txt","Unknown Error in  DefProbCurve object construction");

		if (xmlWCharText)
			delete xmlWCharText;
		xmlWCharText = NULL;

		throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
         "Error in  ZCPY construction");
	}
}


// -------------------------------------------------------------
// construction de la defprob avec stripper calypso  
// -------------------------------------------------------------

ICM_DefaultCurve* ARMLOCAL_XML_DEFPROB_with_calypso_stripper(const std::string& chaineXML,
															 const ARM_ZeroCurve& ircurve,
															 const std::string& newLabel)
{

	MSXML2::IXMLDOMDocumentPtr xmlDoc = XMLTools::CreateDOMDocument30(); 
	VARIANT_BOOL bOK;
	_bstr_t tmp; VariantTools::convert(chaineXML,tmp); 
	xmlDoc->loadXML(tmp,&bOK); 
	
	
	MSXML2::IXMLDOMNodeListPtr resultList; 
	resultList = XMLTools::selectNodes(xmlDoc,"/CurveProbabilityList/CurveProbability/Points/Point") ;
	long nbNodes;
	resultList->get_length(&nbNodes); 
	if (nbNodes<=1) 
		ICMTHROW(ERR_INVALID_ARGUMENT,"ARMLOCAL_XML_DEFPROB_with_Calypso_stripper: more than 1 point required..."); 
	
	std::vector<double> yearFractions(nbNodes-1); 
	ARM_Vector intensities(nbNodes-1); 
	
	ARM_Date AsOfDate ;
	MSXML2::IXMLDOMNodePtr item= XMLTools::selectSingleNode(item,"/CurveProbabilityList/CurveProbability/Date"); 
	XMLTools::convert(item,AsOfDate,"YYYYMMDD hh:mm") ;

	// skip the first point (asof curveDate) 
	for(long i=1;i<nbNodes;i++) 
	{
		MSXML2::IXMLDOMNodePtr item = XMLTools::get_item(resultList,i); 
		MSXML2::IXMLDOMNodePtr item2 = XMLTools::selectSingleNode(item,"Maturity"); 
		ARM_Date date ; XMLTools::convert(item2,date,"YYYYMMDD"); 
		yearFractions[i-1] = CountYears(KACTUAL_365,AsOfDate,date); 
		item2 = XMLTools::selectSingleNode(item,"Mid"); 
		XMLTools::convert(item2,intensities.Elt(i-1)) ; 
	}
	double recovery; 
	item= XMLTools::selectSingleNode(item,"/CurveProbabilityList/CurveProbability/CurveInfo/RecoveryRate"); 
	XMLTools::convert(item,recovery) ;
	std::string issuerLabel; 
	item = XMLTools::selectSingleNode(item,"/CurveProbabilityList/CurveProbability/Issuer/ShortName"); 
	XMLTools::convert(item,issuerLabel) ;
	if (newLabel!="") issuerLabel=newLabel; 
	
	
	qCDS_ADJ adj ;
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

	
	ICM_DefaultCurve* newDefProb = new ICM_Constant_Piecewise(AsOfDate,
														  (double*)&(*yearFractions.begin()),
														  (ARM_Vector*)(&intensities),
														  -1, // survival proba
														  recovery,
														  (ARM_ZeroCurve*)(&ircurve),	// will be cloned. 
														  K_ADJUSTED,	// intRule
														  K_ADJUSTED,	// adjStartDate
														  adj,
														  ircurve.GetCurrencyUnit()->GetCcyName(),
														  issuerLabel,
														  NULL,qDEFCURVE_DICHO,"STD",
														  ARM_Currency(ircurve.GetCurrencyUnit()->GetCcyName()).GetCreditStartDateLag()
														  );
	return newDefProb ;

}
 