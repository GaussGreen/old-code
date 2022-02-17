
#include <ARM\libarm_local\firstToBeIncluded.h>

#include <ARM\libarm_local\arm_local_glob.h>
#include <ARM\libarm_frometk\arm_local_parsexml_nt_security.h>
#include <ARM\libarm_frometk\arm_local_parsexml_util.h>
#include <ARM\libarm_frometk\arm_local_etoolkit.h>
#include <ARM\libarm_frometk\ARM_local_eToolkit_for_ICM.h>
#include <libCCatl\CCatl.h>

#include <glob\linalg.h>
#include <glob\dates.h>
#include <crv\zeroint.h>
#include <ccy\currency.h>
#include <atlbase.h>
#include <memory>

#include <ARM\libarm_frometk\PaserManagerUtilities.h>
#include "VariantTools.h"

#define ICM_RELEASEXML(argICM_RELEASEXML) \
{ if (argICM_RELEASEXML) argICM_RELEASEXML->Release(); argICM_RELEASEXML=0; }

#include <algorithm> 




// ------------------------------------
// Recuperation du Cust Id & du Deal Id
// ------------------------------------

void ARMLOCAL_XML_STDINFOCREDIT(const char* chaineXML, 
								std::string & bookName, 
								std::string & structureId, 
								std::string & custId, 
								std::string & dealId,
								std::string& summitProdType)
{
	VARIANT_BOOL bOK;
	HRESULT hr;

	
	MSXML2::IXMLDOMDocument *XMLDoc = NULL;
	MSXML2::IXMLDOMNode* aNode = NULL;
	wchar_t * xmlWCharText = NULL;

	bookName = "";
	structureId = "";
	custId = "";
	dealId ="";
	summitProdType=""; 
	try
	{
		hr = CoInitialize(NULL); 
		hr = CoCreateInstance(__uuidof(MSXML2::DOMDocument30), NULL, CLSCTX_INPROC_SERVER, __uuidof(MSXML2::IXMLDOMDocument), (void**)&XMLDoc);
		SUCCEEDED(hr) ? 0 : throw hr;

		xmlWCharText = constchar2wchar(chaineXML);

		Local_XML_TRACE("c:\\temp\\log_arm_credit.txt","Lecture XML pour Infos Credit");
		XMLDoc->loadXML((_bstr_t)xmlWCharText, &bOK);

	}
	catch (std::exception&e)
	{
		ICM_RELEASEXML(XMLDoc); 
		// if (XMLDoc) XMLDoc->Release(); 
		if (xmlWCharText)	delete[] xmlWCharText;
		xmlWCharText = NULL;
		ICMTHROW(ERR_INVALID_ARGUMENT,"ARMLOCAL_XML_STDINFOCREDIT: can't create XMLDocument: "<<e.what()); 
	}
	catch(Exception& )
	{
		ICM_RELEASEXML(XMLDoc); 
		// if (XMLDoc) XMLDoc->Release();
		if (xmlWCharText)	delete[] xmlWCharText;
		xmlWCharText = NULL;

		throw; 
	}
	catch(...)
	{
		if (xmlWCharText)	delete[] xmlWCharText;
		xmlWCharText = NULL;
		// if (XMLDoc) XMLDoc->Release();
		ICM_RELEASEXML(XMLDoc); 
		ICMTHROW(ERR_INVALID_ARGUMENT,"ARMLOCAL_XML_STDINFOCREDIT: can't create XMLDocument: unknown exception"); 
	}

	
	try
	{
		// JLA. Retrieve 	
		std::string chemin = "Response/EXOTIC/Env/ENV";
		XMLDoc->selectSingleNode((_bstr_t)chemin.c_str(), &aNode);
		if (aNode)
		{
			custId= GetStringFromXMLNode(aNode, "Cust");
			dealId = GetStringFromXMLNode(aNode, "DealId") ;
			bookName = GetStringFromXMLNode(aNode, "Book") ;
			structureId = GetStringFromXMLNode(aNode, "StructureId") ;
			summitProdType = GetStringFromXMLNode(aNode,"ProductType") ;
			ICM_RELEASEXML(aNode); 
		}
		ICM_RELEASEXML(XMLDoc); 
		delete [] xmlWCharText ; xmlWCharText =0; 
		return ; 
	}
	catch(std::exception& e)
	{
		ICM_RELEASEXML(aNode); 
		ICM_RELEASEXML(XMLDoc); 
		// if (aNode) aNode->Release();
		// if (XMLDoc) XMLDoc->Release();
		if (xmlWCharText)delete[] xmlWCharText;
		ICMTHROW(ERR_INVALID_ARGUMENT,"ARMLOCAL_XML_STDINFOCREDIT: Parsing Error: "<<e.what()); 
	}
	catch(Exception& )
	{
		ICM_RELEASEXML(aNode); 
		ICM_RELEASEXML(XMLDoc); 
		// if (aNode) aNode->Release();
		// if (XMLDoc) XMLDoc->Release();
		if (xmlWCharText)	delete[] xmlWCharText;
		throw; 
	}
	catch(...)
	{		
		ICM_RELEASEXML(aNode); 
		ICM_RELEASEXML(XMLDoc); 
		//if (aNode) aNode->Release();
		//if (XMLDoc) XMLDoc->Release();
		if (xmlWCharText)	delete[] xmlWCharText;
		ICMTHROW(ERR_INVALID_ARGUMENT,"ARMLOCAL_XML_STDINFOCREDIT: Parsing Error: unknown exception"); 
	}
}


// ------------------------------------
// Recuperation des events
// ------------------------------------

ARM_ReferenceValue* ARMLOCAL_XML_EVENT(const char* chaineXML, CCString EventType)
{
	VARIANT_BOOL bOK;
	HRESULT hr;

	MSXML2::IXMLDOMDocument *XMLDoc = NULL;
	MSXML2::IXMLDOMNodeList * resultList = NULL;
	MSXML2::IXMLDOMNode * listItem = NULL, * theNode = NULL;
	long nbNodes = 0;
	int k =0;	

	ARM_Vector* Dates = NULL;
	ARM_Vector* Values = NULL;
	ARM_ReferenceValue* refval = NULL;

	wchar_t * xmlWCharText = NULL;

	CCString EventTypeMod; 			

	try
	{
		hr = CoInitialize(NULL); 
		hr = CoCreateInstance(__uuidof(MSXML2::DOMDocument30), NULL, CLSCTX_INPROC_SERVER, __uuidof(MSXML2::IXMLDOMDocument), (void**)&XMLDoc);
		SUCCEEDED(hr) ? 0 : throw hr;

		xmlWCharText = constchar2wchar(chaineXML);

		Local_XML_TRACE("c:\\temp\\log_arm_credit.txt","Lecture XML pour Events Credit");
		XMLDoc->loadXML((_bstr_t)xmlWCharText, &bOK);

	}
	catch(Exception& )
	{
		if (XMLDoc) XMLDoc->Release();
		hr = S_FALSE;

		if (xmlWCharText)
			delete[] xmlWCharText;
		xmlWCharText = NULL;

		throw; 
	}
	catch(...)
	{
		if (xmlWCharText)
			delete[] xmlWCharText;
		xmlWCharText = NULL;

		if (XMLDoc) XMLDoc->Release();
		hr = S_FALSE;
		throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"Pb in creating XML document for getting Tranches informations");
	}

	try
	{
		EventTypeMod = (CCString) "//ASSET/Events/EVENT[Type = \""+ EventType + (CCString) "\"]";

		// Récupération des events
		if (XMLDoc->selectNodes((_bstr_t)(EventTypeMod), &resultList) == S_OK)
		{
			resultList->get_length(&nbNodes);

			if (nbNodes == 0)
			{
				if (theNode) theNode->Release();theNode = NULL;
				if (resultList)	resultList->Release();	resultList = NULL;
				if (listItem) listItem->Release(); listItem = NULL;
				if (XMLDoc) XMLDoc->Release(); XMLDoc = NULL;
				if (xmlWCharText) delete[] xmlWCharText;	xmlWCharText = NULL;

				return NULL;
			}		

			ARM_Vector* Dates = new ARM_Vector(nbNodes,0.);
			ARM_Vector* Values = new ARM_Vector(nbNodes,0.);

			int k =0;

			for (long indexNode=0 ; indexNode<nbNodes ; indexNode++)
			{
				hr=resultList->get_item(indexNode, &listItem);

				//Recuperation de la Date
				Dates->InitElt(indexNode,GetDateFromXMLNode(listItem,"ADate").GetJulian());

				double Amount= XML_doubleNodeTreating(listItem,"Amount");
				Values->InitElt(indexNode,Amount);

				if (listItem)
				{
				listItem->Release();
				listItem = NULL;
				}

			}

			if (resultList)
			{
				resultList->Release();
				resultList = NULL;
			}


			refval = new ARM_ReferenceValue(Dates,Values,0);
			refval->SetCalcMethod(K_STEPUP_RIGHT);
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

		if (listItem)
		{
		listItem->Release();
		listItem = NULL;
		}

		if (XMLDoc)	
		{
		XMLDoc->Release();
		XMLDoc = NULL;
		}

		if (xmlWCharText)
			delete[] xmlWCharText;
		xmlWCharText = NULL;


		return (refval);
	}
	catch(Exception& )
	{
		if (theNode) theNode->Release();
		if (listItem) listItem->Release();
		if (resultList) resultList->Release();
		if (XMLDoc) XMLDoc->Release();

		if (xmlWCharText)
			delete[] xmlWCharText;
		xmlWCharText = NULL;

		hr = S_FALSE;
		throw; 
	}
	catch(...)
	{		
		if (theNode) theNode->Release();
		if (listItem) listItem->Release();
		if (resultList) resultList->Release();
		if (XMLDoc) XMLDoc->Release();

		if (xmlWCharText)
			delete[] xmlWCharText;
		xmlWCharText = NULL;

		hr = S_FALSE;
		throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
         "Error in XML parsing for getting Tranche");
	}
}


// ------------------------------------
// Recuperation des events for fixing
// ------------------------------------

void ARMLOCAL_XML_EVENT_FIXING(const char* chaineXML, CCString EventType,vector<double>& FixDates,vector<double>& Fixings)
{
	VARIANT_BOOL bOK;
	HRESULT hr;

	MSXML2::IXMLDOMDocument *XMLDoc = NULL;
	MSXML2::IXMLDOMNodeList * resultList = NULL;
	MSXML2::IXMLDOMNode * listItem = NULL, * theNode = NULL;
	long nbNodes = 0;
	int k =0;	

	FixDates.resize(0);
	Fixings.resize(0);

	wchar_t * xmlWCharText = NULL;

	CCString EventTypeMod; 			

	try
	{
		hr = CoInitialize(NULL); 
//		SUCCEEDED(hr) ? 0 : throw hr;
		hr = CoCreateInstance(__uuidof(MSXML2::DOMDocument30), NULL, CLSCTX_INPROC_SERVER, __uuidof(MSXML2::IXMLDOMDocument), (void**)&XMLDoc);
		SUCCEEDED(hr) ? 0 : throw hr;

		xmlWCharText = constchar2wchar(chaineXML);

		Local_XML_TRACE("c:\\temp\\log_arm_credit.txt","Lecture XML pour Events Credit");
		XMLDoc->loadXML((_bstr_t)xmlWCharText, &bOK);

	}
	catch(Exception& )
	{
		if (XMLDoc) XMLDoc->Release();
		hr = S_FALSE;

		if (xmlWCharText)
			delete[] xmlWCharText;
		xmlWCharText = NULL;

		throw; 
	}
	catch(...)
	{
		if (xmlWCharText)
			delete[] xmlWCharText;
		xmlWCharText = NULL;

		if (XMLDoc) XMLDoc->Release();
		hr = S_FALSE;
		throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"Pb in creating XML document for getting Tranches informations");
	}

	try
	{
		EventTypeMod = (CCString) "//ASSET/Events/EVENT[Type = \""+ EventType + (CCString) "\"]";

		// Récupération des events
		if (XMLDoc->selectNodes((_bstr_t)(EventTypeMod), &resultList) == S_OK)
		{
			resultList->get_length(&nbNodes);

			if (nbNodes == 0)
			{
				if (theNode) theNode->Release();theNode = NULL;
				if (resultList)	resultList->Release();	resultList = NULL;
				if (listItem) listItem->Release(); listItem = NULL;
				if (XMLDoc) XMLDoc->Release(); XMLDoc = NULL;
				if (xmlWCharText) delete[] xmlWCharText;	xmlWCharText = NULL;

				return;
			}		

			FixDates.resize(nbNodes);
			Fixings.resize(nbNodes);

			int k =0;

			for (long indexNode=0 ; indexNode<nbNodes ; indexNode++)
			{
				hr=resultList->get_item(indexNode, &listItem);

				//Recuperation de la Date
				FixDates[indexNode] = GetDateFromXMLNode(listItem,"ADate").GetJulian();

				double Amount= XML_doubleNodeTreating(listItem,"Amount");
				Fixings[indexNode]= Amount;

				if (listItem)
				{
				listItem->Release();
				listItem = NULL;
				}

			}

			if (resultList)
			{
				resultList->Release();
				resultList = NULL;
			}

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

		if (listItem)
		{
		listItem->Release();
		listItem = NULL;
		}

		if (XMLDoc)	
		{
		XMLDoc->Release();
		XMLDoc = NULL;
		}

		if (xmlWCharText)
			delete[] xmlWCharText;
		xmlWCharText = NULL;


		return;
	}
	catch(Exception& )
	{
		if (theNode) theNode->Release();
		if (listItem) listItem->Release();
		if (resultList) resultList->Release();
		if (XMLDoc) XMLDoc->Release();

		if (xmlWCharText)
			delete[] xmlWCharText;
		xmlWCharText = NULL;

		hr = S_FALSE;
		throw; 
	}
	catch(...)
	{		
		if (theNode) theNode->Release();
		if (listItem) listItem->Release();
		if (resultList) resultList->Release();
		if (XMLDoc) XMLDoc->Release();

		if (xmlWCharText)
			delete[] xmlWCharText;
		xmlWCharText = NULL;

		hr = S_FALSE;
		throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
         "Error in XML parsing for getting Tranche");
	}
}



// ------------------------------------
// Generic Tranche Parsing for Ntd 
// ------------------------------------

ICM_Nthtd* ARMLOCAL_XML_NTD(const char* chaineXML, CCString& bookId, CCString& custId, CCString& dealId,
								ARM_ReferenceValue*fees)
{
	VARIANT_BOOL bOK;
	HRESULT hr;

	MSXML2::IXMLDOMDocument *XMLDoc = NULL;
	MSXML2::IXMLDOMNodeList * resultList = NULL;
	MSXML2::IXMLDOMNodeList * resultList0 = NULL;
	MSXML2::IXMLDOMNode * listItem = NULL, * theNode = NULL;
	long nbNodes = 0;

	ARM_Date StartDate;
	ARM_Date EndDate;
	ARM_Date fstCpEffectDate ;
	bool isFstCpEffect(false); 
	double spread = 0.;

	char ReferenceDate[11];
	int payFreq = K_DEF_FREQ;
	int daycount = KACTUAL_360;
	int ACC_daycount = KACTUAL_360;
	int DelivDays=0; 

	double FixedPayerAmount=0.;
	double FloatingPayerAmount=0.;

	// ARM_Currency* ccy = NULL;
	std::string ccy ;

	double inferior = 0.;
	double superior = 0.;
	double Binary = 0.;

	int nbissuers = 0;
	// char** Issuers = NULL; 
	std::vector<std::string> Issuers ; 
	double* Notionals = NULL;
	double* Rate = NULL;
	// char* PayCal = NULL;
	std::string PayCal ;

	int interestRule, adjStartDateRule; 
	double SensDeLaTransaction = 1.;

	//ARM_Currency* discountCcy = NULL;
	//ARM_ReferenceValue* notional = NULL;

	int k =0;	

	wchar_t * xmlWCharText = NULL;

	ARM_Date startStubDate;
	ARM_Date endStubDate;
	int refDay = 0;
	int stub = K_SHORTSTART;

	std::string tradeId ; 
	bool excludeMaturity = EXCLUDE_MATURITY ;

	// not used for FTD ARM_Date Protect_StartDate;
	// not used for FTD char TProtect_StartDate[11];
	ARM_Date Protect_EndDate;
	char TProtect_EndDate[11];

	try
	{
		hr = CoInitialize(NULL); 
//		SUCCEEDED(hr) ? 0 : throw hr;
		hr = CoCreateInstance(__uuidof(MSXML2::DOMDocument30), NULL, CLSCTX_INPROC_SERVER, __uuidof(MSXML2::IXMLDOMDocument), (void**)&XMLDoc);
		SUCCEEDED(hr) ? 0 : throw hr;

		xmlWCharText = constchar2wchar(chaineXML);

		Local_XML_TRACE("c:\\temp\\log_arm_credit.txt","Lecture XML pour Matrice de Correlation");
		XMLDoc->loadXML((_bstr_t)xmlWCharText, &bOK);

		#ifdef _DEBUG
		CCString chtmp = "c:\\temp\\AAA_NTD_OUTPUT.xml";
		VARIANT v_chtmp;
		CCString2VARIANT (chtmp, &v_chtmp);
		XMLDoc->save(v_chtmp);
		#endif

	}
	catch(Exception& )
	{
		if (XMLDoc) XMLDoc->Release();
		hr = S_FALSE;
		throw; 


		if (xmlWCharText)
			delete[] xmlWCharText;
		xmlWCharText = NULL;
	}
	catch(...)
	{
		if (XMLDoc) XMLDoc->Release();
		hr = S_FALSE;

		if (xmlWCharText)
			delete[] xmlWCharText;
		xmlWCharText = NULL;

		throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"Pb in creating XML document for getting Tranches informations");
	}

	try
	{
		{
			// JLA. Retrieve 
			MSXML2::IXMLDOMNode* aNode = NULL;
			std::string chemin = "Response/EXOTIC/Env/ENV";
			XMLDoc->selectSingleNode((_bstr_t)chemin.c_str(), &aNode);
			if (aNode)
			{
				custId= GetStringFromXMLNode(aNode, "Cust").c_str();
				dealId = GetStringFromXMLNode(aNode, "DealId").c_str();
				bookId = GetStringFromXMLNode(aNode, "Book").c_str();
				tradeId = GetStringFromXMLNode(aNode, "TradeId") ;
				aNode->Release();
			}
		}


		// ---------------------------------------------------------------------------------
		// Recuperation des infos communes à tous les deals
		// ---------------------------------------------------------------------------------

		Local_XML_TRACE("c:\\temp\\log_arm_credit.txt","Recuperation des infos communes à tous les deals ");

		if (fees) {
			// JLA. Retrieve 
			MSXML2::IXMLDOMNode* aNode = NULL;
			std::string chemin = "//ASSET";
			XMLDoc->selectSingleNode((_bstr_t)chemin.c_str(), &aNode);
			ARM_ReferenceValue* item= 0 ;
			if (aNode) 
			{ 
				item =GetPrimes(aNode) ; 
				if (item) *fees=*item; 
			}
			if (item) delete item; item=0; 
			aNode->Release(); 
		}


 		if (XMLDoc->selectNodes((_bstr_t)("//ASSET"), &resultList) == S_OK)
		{

			hr=resultList->get_item(0, &listItem);
			
			//Recuperation de la Trade Date
			StartDate = GetStartDate(listItem);

			//Recuperation de la Maturity Date
			EndDate = GetEndDate(listItem);

			//Recuperation de frequence
			listItem->selectSingleNode(_bstr_t("SCHED_Pay_Freq"), &theNode);
			if (theNode!=NULL)
			{
				BSTR resultat = NULL;
				theNode->get_text(&resultat);

				_bstr_t ff(resultat,false);
				char * ff1=(char *)ff;
				payFreq = FromSummitFreqToARMFreq((const char*) ff1);

				theNode->Release();
				theNode=NULL;
				if (resultat) SysFreeString(resultat);
			}

			//Recuperation du calendrier
			listItem->selectSingleNode(_bstr_t("SCHED_Pay_Cal"), &theNode);
			if (theNode!=NULL)
			{
				BSTR resultat = NULL;
				theNode->get_text(&resultat);

				_bstr_t ff(resultat,false);
				PayCal =  ff;
				// if (PayCal)
				//	delete[] PayCal;
				//PayCal = new char[5];
				//st// rcpy(PayCal,ff1);


				theNode->Release();
				theNode=NULL;
				if (resultat) SysFreeString(resultat);
			}

			//Recuperation du basis daycount
			listItem->selectSingleNode(_bstr_t("INTEREST_Basis"), &theNode);
			if (theNode!=NULL)
			{
				BSTR resultat = NULL;
				theNode->get_text(&resultat);

				_bstr_t ff(resultat,false);
				char * ff1=(char *)ff;
				daycount = FromSummitDaycountToARMDaycount((const char*) ff1);

				theNode->Release();
				theNode=NULL;
				if (resultat) SysFreeString(resultat);
			}

			//Recuperation du Accrual basis daycount
			listItem->selectSingleNode(_bstr_t("INT_ACC_Basis"), &theNode);
			if (theNode!=NULL)
			{
				BSTR resultat = NULL;
				theNode->get_text(&resultat);

				_bstr_t ff(resultat,false);
				char * ff1=(char *)ff;
				ACC_daycount = FromSummitDaycountToARMDaycount((const char*) ff1);

				theNode->Release();
				theNode=NULL;
				if (resultat) SysFreeString(resultat);
			}

			//Recuperation de la Reference Date
			listItem->selectSingleNode(_bstr_t("STUB_Date1"), &theNode);
			
			if (theNode!=NULL)
			{
				BSTR resultat = NULL;
				theNode->get_text(&resultat);
				_bstr_t ff(resultat,false);
				char * ff1=(char *)ff;

				if (strcmp((const char*)ff1,"") == 0)
					strcpy(ReferenceDate,"NULL");
				else
				{
					ARM_Date tmpDate (ff1,"YYYYMMDD");
					tmpDate.JulianToStrDate(ReferenceDate); 
				}

				theNode->Release();
				theNode=NULL;
				if (resultat) SysFreeString(resultat);
			}
			else
				strcpy(ReferenceDate,"NULL");

			// reference Date
			listItem->selectSingleNode(_bstr_t((const char *)"STUB_Date2"), &theNode);
			if (theNode!=NULL)
			{
				BSTR resultat = NULL;
				theNode->get_text(&resultat);

				_bstr_t ff(resultat,false);
				char * ff1=(char *)ff;

				if (strcmp((const char*)ff1,"") != 0)
				{
				ARM_Date tmpDate (ff1,"YYYYMMDD");
				endStubDate = tmpDate;
				}

				theNode->Release();
				theNode=NULL;
				if (resultat) SysFreeString(resultat);
			}
			// Fst Cpn Effective Date
			
			listItem->selectSingleNode(_bstr_t((const char *)"FstCpEffect"), &theNode);
			if (theNode!=NULL)
			{
				BSTR resultat = NULL;
				theNode->get_text(&resultat);

				_bstr_t ff(resultat,false);
				char * ff1=(char *)ff;

				if (strcmp((const char*)ff1,"") != 0)
				{
					isFstCpEffect=true ; 
					ARM_Date tmpDate (ff1,"YYYYMMDD");
					fstCpEffectDate= tmpDate;
				}

				theNode->Release();
				theNode=NULL;
				if (resultat) SysFreeString(resultat);
			}

			// interestRule
			interestRule = GetIntRule(listItem) ;
			if (interestRule==K_UNADJUSTED) adjStartDateRule=K_UNADJUSTED; 
			else adjStartDateRule=K_ADJUSTED ; 

			// reference Day1
			refDay  = GetIntFromXMLNode (listItem,"SCHED_Pay_AnnDay"); 

			//Recuperation du stub
			listItem->selectSingleNode(_bstr_t("STUB_StubType"), &theNode);
			if (theNode!=NULL)
			{
				BSTR resultat = NULL;
				theNode->get_text(&resultat);

				_bstr_t ff(resultat,false);
				char * ff1=(char *)ff;
				stub = SummitStub2ARMStub(ff1);

				theNode->Release();
				theNode=NULL;
				if (resultat) SysFreeString(resultat);
			}

			listItem->selectSingleNode(_bstr_t("Ccy"), &theNode);
			if (theNode!=NULL)
			{
				BSTR resultat = NULL;
				theNode->get_text(&resultat);

				_bstr_t ff(resultat,false);
				ccy = ff ; 
				// char * ff1=(char *)ff;
				// ccy = new ARM_Currency(ff1);

				theNode->Release();
				theNode=NULL;
				if (resultat) SysFreeString(resultat);
			}

			if (startStubDate<=StartDate || startStubDate>=EndDate) startStubDate=ARM_Date(); 
			if (endStubDate<=StartDate || endStubDate>=EndDate) endStubDate=ARM_Date(); 
			DeduceRefDateAndStub(StartDate,EndDate,startStubDate,endStubDate,payFreq,refDay,ccy,(char*)ReferenceDate,stub);

			FixedPayerAmount= XML_doubleNodeTreating(listItem,"Notional");
			FloatingPayerAmount=FixedPayerAmount;

			//Recuperation du PorS
			listItem->selectSingleNode(_bstr_t("PorS"), &theNode);
			if (theNode!=NULL)
			{
				BSTR resultat = NULL;
				theNode->get_text(&resultat);

				_bstr_t ff(resultat,false);
				char * ff1=(char *)ff;
				
				if (!strcmp(ff1,"P"))
					SensDeLaTransaction = K_RCV;
				else
					SensDeLaTransaction = K_PAY;

				theNode->Release();
				theNode=NULL;
				if (resultat) SysFreeString(resultat);
			}

			//	Recuperation de INTEREST_Method
			{
				std::string method ;
				if ( GetStringFromXMLNode(listItem, "INTEREST_Method",method) ) 
				{
					if (method =="") excludeMaturity=EXCLUDE_MATURITY ; 
					else if (method =="STAND") excludeMaturity=EXCLUDE_MATURITY ; 
					else if (method =="INCL") excludeMaturity=INCLUDE_MATURITY; 
					// else if (method =="EFF") excludeMaturity=INCLUDE_MATURITY ; 
					else 
					{
						ICMTHROW(ERR_INVALID_ARGUMENT,tradeId<<": Can't understand INTEREST_Method. "<< method<<" for " << tradeId ); 
					}
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

		// ---------------------------------------------------------------------------------
		// Recuperation du spread fix
		// ---------------------------------------------------------------------------------

	
		if (XMLDoc->selectNodes((_bstr_t)("//ASSET/ProdData/CREDSWAP"), &resultList) == S_OK)
		{
			resultList->get_length(&nbNodes);

			if (nbNodes == 0)
			{
				hr = S_FALSE;
				ICMTHROW(ERR_INVALID_ARGUMENT,tradeId<<": Invalid XML string for getting //ASSET/ProdData/CREDSWAP");
			}

			hr=resultList->get_item(0, &listItem);
			spread = XML_doubleNodeTreating(listItem,"Formula");
			spread /= 100.;

			// listItem->Release();
			// listItem=NULL;

			listItem->selectSingleNode(_bstr_t("DelivDays"), &theNode);
			if (theNode!=NULL)
			{
				BSTR resultat = NULL;
				theNode->get_text(&resultat);

				_bstr_t ff(resultat,false);
				char * ff1=(char *)ff;
				
				sscanf(ff1, "%2dD", &DelivDays);

				theNode->Release();
				theNode=NULL;
				if (resultat) SysFreeString(resultat);
			}

			//Recuperation de la Protect End Date
			listItem->selectSingleNode(_bstr_t("MatDate"), &theNode);
			
			if (theNode!=NULL)
			{
				BSTR resultat = NULL;
				theNode->get_text(&resultat);
				_bstr_t ff(resultat,false);
				char * ff1=(char *)ff;

				if (strcmp((const char*)ff1,"") == 0)
					strcpy(TProtect_EndDate,"NULL");
				else
				{
					ARM_Date tmpDate (ff1,"YYYYMMDD");
					tmpDate.JulianToStrDate(TProtect_EndDate); 
				}

				theNode->Release();
				theNode=NULL;
				if (resultat) SysFreeString(resultat);
			}
			else
				strcpy(TProtect_EndDate,"NULL");	

			// LEAK: hr=resultList->get_item(0, &listItem);
			listItem->selectSingleNode(_bstr_t("SettleMethod"), &theNode);
			if (theNode!=NULL)
			{
				BSTR resultat = NULL;
				theNode->get_text(&resultat);

				_bstr_t ff(resultat,false);
				char * ff1=(char *)ff;
				
				if (strcmp(ff1,"BIN") != NULL) Binary = CREDIT_DEFAULT_VALUE;

				theNode->Release();
				theNode=NULL;
				if (resultat) SysFreeString(resultat);
			}
			listItem->Release();
			listItem=NULL;

			if (Binary != CREDIT_DEFAULT_VALUE)
			{
				hr=resultList->get_item(0, &listItem);
				Binary = 1. - XML_doubleNodeTreating(listItem,"BinaryAmt");
				listItem->Release();
				listItem=NULL;
			}
		}

		if (resultList)
		{
			resultList->Release();
			resultList = NULL;
		}

		// ---------------------------------------------------------------------------------
		// Recuperation des infos sur les Issuers
		// ---------------------------------------------------------------------------------

		Local_XML_TRACE("c:\\temp\\log_arm_credit.txt","Recuperation des infos sur les Issuers");

			if (XMLDoc->selectNodes((_bstr_t)("//ASSET/ProdData/CREDSWAP/CreditEntList/CREDIT[DefpIssuer != \"FTD1\"]"), &resultList0) == S_OK)
			{
				resultList0->get_length(&nbNodes);

				if (nbNodes == 0)
				{	
					hr = S_FALSE;
					ICMTHROW(ERR_INVALID_ARGUMENT,tradeId<<": Invalid XML string for getting //ASSET/ProdData/CREDSWAP/CreditEntList/CREDIT[DefpIssuer != \"FTD1\"]");
				}

				// Issuers = new char*[nbNodes];
				Issuers.resize(nbNodes); 
				nbissuers = nbNodes; 

				for (long indexNode=0 ; indexNode<nbNodes ; indexNode++)
				{
					hr=resultList0->get_item(indexNode, &listItem);

					listItem->selectSingleNode(_bstr_t("DefpIssuer"), &theNode);
					if (theNode!=NULL)
					{
					BSTR resultat = NULL;
					theNode->get_text(&resultat);

					_bstr_t ff(resultat,false);
					char * ff1=(char *)ff;

					// Issuers[indexNode] = new char[1000];
					// strcpy(Issuers[indexNode],ff1);
					Issuers[indexNode]=ff1 ;
		
					theNode->Release();
					theNode=NULL;
					if (resultat) SysFreeString(resultat);
					}

					listItem->Release();
					listItem = NULL;

				}

				if (resultList0)
				{
					resultList0->Release();
					resultList0 = NULL;
				}


			if (XMLDoc->selectNodes((_bstr_t)("//ASSET/ProdData/CREDSWAP/CreditEntList/CREDIT"), &resultList) == S_OK)
			{
				resultList->get_length(&nbNodes);

				if (nbNodes == 0)
				{	
					hr = S_FALSE;
					ICMTHROW(ERR_INVALID_ARGUMENT,tradeId<<": Invalid XML string for getting //ASSET/ProdData/CREDSWAP/CreditEntList/CREDIT");
				}

				Notionals = new double[nbNodes];

				for (long indexNode=0 ; indexNode<nbNodes ; indexNode++)
				{
					hr=resultList->get_item(indexNode, &listItem);

					Notionals[indexNode] = XML_doubleNodeTreating(listItem,"Quantity");

					listItem->Release();
					listItem = NULL;

				}

				// Useless ??: 
				if (theNode)
				{
				theNode->Release();
				theNode=NULL;
				}
	
			}

			if (resultList)
			{
			resultList->Release();
			resultList = NULL;
			}

			if (XMLDoc->selectNodes((_bstr_t)("//ASSET/ProdData/CREDSWAP"), &resultList) == S_OK)
			{
				resultList->get_length(&nbNodes);

				if (nbNodes == 0)
				{
				hr = S_FALSE;
				ICMTHROW(ERR_INVALID_ARGUMENT,tradeId<<": Invalid XML string for getting //ASSET/ProdData/CREDSWAP");
				}

				hr=resultList->get_item(0, &listItem);

				listItem->selectSingleNode(_bstr_t("Tranched"), &theNode);
				if (theNode!=NULL)
				{
					BSTR resultat = NULL;
					theNode->get_text(&resultat);

					_bstr_t ff(resultat,false);
					char * ff1=(char *)ff;

					theNode->Release();
					theNode=NULL;
					if (resultat) SysFreeString(resultat);
				}

				inferior = XML_doubleNodeTreating(listItem,"StartDefaultNo");
				superior = XML_doubleNodeTreating(listItem,"EndDefaultNo");

				if (theNode)
				{
				theNode->Release();
				theNode=NULL;
				}

				listItem->Release();
				listItem=NULL;
			}

			if (resultList)
			{
			resultList->Release();
			resultList = NULL;
			}
		}

		if (resultList)
		{
			resultList->Release();
			resultList = NULL;
		}

		//------------------------------------------------------------------------------
		//Création de l'instrument financier
		//------------------------------------------------------------------------------

		Local_XML_TRACE("c:\\temp\\log_arm_credit.txt","Création du NTD");

		if (EndDate <= StartDate)
		{
			Local_XML_TRACE("c:\\temp\\log_arm_credit.txt","Error: EndDate <= StartDate");
			ICMTHROW(ERR_INVALID_ARGUMENT,tradeId<<": Error: EndDate <= StartDate");
		}

		if (inferior > superior)
		{
			Local_XML_TRACE("c:\\temp\\log_arm_credit.txt","Error: inferior > superior");
			ICMTHROW(ERR_INVALID_ARGUMENT,tradeId<<": Error: EndDate <= StartDate");
		}

		if (nbissuers == 0)
		{
			Local_XML_TRACE("c:\\temp\\log_arm_credit.txt","Error: inferior > superior");
			ICMTHROW(ERR_INVALID_ARGUMENT,tradeId<<": Error: nbissuers = 0");
		}

		bool toto (isFstCpEffect) ; 
		if (isFstCpEffect && fstCpEffectDate>=StartDate) toto=false;
		ICM_Nthtd* ntd = new ICM_Nthtd(StartDate,
										EndDate,
										strcmp(ReferenceDate,"NULL")==0 ? 0 : &ARM_Date(ReferenceDate),
										toto?&fstCpEffectDate:0,
										spread,
										interestRule,		// intRule
										adjStartDateRule,	// adjStartDate
										inferior,
										superior,
										// nbissuers,
										Issuers,
										// ReferenceDate,
										payFreq,
										daycount,
										SensDeLaTransaction*FixedPayerAmount, 
										qACCRUED_SETTLED,
										ccy,
										SensDeLaTransaction*FloatingPayerAmount,
										stub,
										// DEFAULT_CREDIT_LAG,
										DelivDays,
										DEFAULT_FRQ_DEFLEG,
										Binary,
										PayCal,
										excludeMaturity);

		ntd->SetPorS(SensDeLaTransaction);

		//ARM_ReferenceValue* XNL = ARMLOCAL_XML_EVENT(chaineXML,"XNL");
		//if (XNL) ntd->GetFeeLeg()->SetAmount(XNL);
		//if (XNL) delete (XNL);

		//if (ccy) delete ccy;

		Local_XML_TRACE("c:\\temp\\log_arm_credit.txt","Fin Création du NTD");

		// for (int i = 0; i<nbissuers; i++)
		// {
		// 	if (Issuers[i]) delete[] Issuers[i];
		// 	Issuers[i] = NULL;
		// }

		// if (Issuers) delete[] Issuers;
		if (Notionals) delete[] Notionals;		
		if (Rate) delete[] Rate;		
		//if (PayCal) delete[] PayCal;

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

		if (XMLDoc)	{ XMLDoc->Release();}
		XMLDoc = NULL;

		if (xmlWCharText)
			delete[] xmlWCharText;
		xmlWCharText = NULL;

		return (ntd);
	}
	catch(Exception& )
	{
		if (theNode) theNode->Release();
		if (listItem) listItem->Release();
		if (resultList) resultList->Release();
		if (XMLDoc) XMLDoc->Release();


		if (xmlWCharText)
			delete[] xmlWCharText;
		xmlWCharText = NULL;

		hr = S_FALSE;
		throw; 
		// throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,Tmp);
	}
	catch(...)
	{		
		if (theNode) theNode->Release();
		if (listItem) listItem->Release();
		if (resultList) resultList->Release();
		if (XMLDoc) XMLDoc->Release();

		if (xmlWCharText)
			delete[] xmlWCharText;
		xmlWCharText = NULL;

		hr = S_FALSE;
		ICMTHROW(ERR_INVALID_ARGUMENT,tradeId<<": Error in XML parsing for getting Tranche");
	}

	return NULL;
}



/**		JLA useless

  
// ------------------------------------
// Generic Tranche Parsing for Cdo 
// ------------------------------------


ICM_Mez* ARMLOCAL_XML_CDO(const char* chaineXML, CCString& bookId, CCString& custId, CCString& dealId,
						  ARM_ReferenceValue*fees)
{
	VARIANT_BOOL bOK;
	HRESULT hr;

	MSXML2::IXMLDOMDocument *XMLDoc = NULL;
	MSXML2::IXMLDOMNodeList * resultList = NULL;
	MSXML2::IXMLDOMNodeList * resultList0 = NULL;
	MSXML2::IXMLDOMNode * listItem = NULL, * theNode = NULL;
	long nbNodes = 0;

	//ARM_Currency* INDEX_ccy = NULL;
	string index_ccy;
	string index_term;
	string index_name;

	bool INDEX_spreadonly = false;

	CCString strtmp;
	
	ARM_Date StartDate;
	ARM_Date EndDate;
	ARM_Date fstCpEffectDate ;
	bool isFstCpEffect(false); 
	bool isStubStart(false);
	bool isStubEnd(false); 
	double spread = 0.;

	char ReferenceDate[11]; strcpy(ReferenceDate,"NULL"); 
	int payFreq = K_DEF_FREQ;
	int daycount = KACTUAL_360;
	int ACC_daycount = KACTUAL_360;
	int DelivDays=0 ;
	int interestRule,adjStartDateRule ;

	double FixedPayerAmount=0.;
	double FloatingPayerAmount=0.;
	// ARM_Currency* ccy = NULL;
	std::string ccy ; 

	double MezzAmount = 0.;
	double SubAmount = 0.;

	int nbissuers = 0;
	char** Issuers = NULL; 
	double* Notionals = NULL;
	double* Rate = NULL;
	double CoefAdjust = 0.;
	std::string PayCal ;

	double SensDeLaTransaction = 1.;
	double Binary = 0.;
	

	//ARM_Currency* discountCcy = NULL;
	// ARM_ReferenceValue* notional = NULL;
	ICM_ProportionsInfo* info = NULL;
	char SmilProp[500];
	memset(SmilProp,'\0',sizeof(SmilProp));

	int k =0;	
	int IntegrationMethod = CREDIT_DEFAULT_VALUE;
	
	wchar_t * xmlWCharText = NULL;

	ARM_Date startStubDate;
	ARM_Date endStubDate;
	int refDay = 0;
	int stub = K_SHORTSTART;
	int stub0(stub); 

	std::string tradeId ; 
	bool excludeMaturity = EXCLUDE_MATURITY ;

	ARM_Date Protect_StartDate;
	char TProtect_StartDate[11];
	ARM_Date Protect_EndDate;
	char TProtect_EndDate[11];

	qCredit_Leg_Type FeeLegType = qRunning_Leg;
	qCredit_Leg_Type DefaultLegType = qStandart_Recovery_Leg;

	ICM_Leg* FeeLeg = NULL;
	ICM_Leg* DefLeg = NULL;

	try
	{
		hr = CoInitialize(NULL); 
		hr = CoCreateInstance(__uuidof(MSXML2::DOMDocument30), NULL, CLSCTX_INPROC_SERVER, __uuidof(MSXML2::IXMLDOMDocument), (void**)&XMLDoc);
		SUCCEEDED(hr) ? 0 : throw hr;

		xmlWCharText = constchar2wchar(chaineXML);

		Local_XML_TRACE("c:\\temp\\log_arm_credit.txt","Lecture XML pour Matrice de Correlation");
		XMLDoc->loadXML((_bstr_t)xmlWCharText, &bOK);

		
		#ifdef _DEBUG
		CCString chtmp = "c:\\temp\\AAA_CDO_OUTPUT.xml";
		VARIANT v_chtmp;
		CCString2VARIANT (chtmp, &v_chtmp);
		XMLDoc->save(v_chtmp);
		#endif
	}
	catch(Exception& x)
	{
		if (XMLDoc) XMLDoc->Release();
		hr = S_FALSE;
		char Tmp[1000];
		x.GetErrorMessage(Tmp);
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,Tmp);


		if (xmlWCharText)
			delete[] xmlWCharText;
		xmlWCharText = NULL;
	}
	catch(...)
	{
		if (XMLDoc) XMLDoc->Release();
		hr = S_FALSE;

		if (xmlWCharText)
			delete[] xmlWCharText;
		xmlWCharText = NULL;

		throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"Pb in creating XML document for getting Tranches informations");
	}

	try
	{
		{
			// JLA. Retrieve 
			MSXML2::IXMLDOMNode* aNode = NULL;
			std::string chemin = "Response/EXOTIC/Env/ENV";
			XMLDoc->selectSingleNode((_bstr_t)chemin.c_str(), &aNode);
			if (aNode)
			{
			custId= GetStringFromXMLNode(aNode, "Cust").c_str();
			dealId = GetStringFromXMLNode(aNode, "DealId").c_str();
			bookId = GetStringFromXMLNode(aNode, "Book").c_str();
			tradeId = GetStringFromXMLNode(aNode, "TradeId").c_str();
			aNode->Release();
			}
		}

		if (fees) {
			// JLA. Retrieve 
			MSXML2::IXMLDOMNode* aNode = NULL;
			std::string chemin = "//ASSET";
			XMLDoc->selectSingleNode((_bstr_t)chemin.c_str(), &aNode);
			ARM_ReferenceValue* item= 0 ;
			if (aNode) 
			{ 
				item =GetPrimes(aNode) ; 
				if (item) *fees=*item; 
			}
			if (item) delete item; item=0; 
			aNode->Release(); 
		}

		// ------------------------------------------------------------------------
		// Pricing avec correl par strike ?
		// ------------------------------------------------------------------------
		strtmp = "//ASSET/ProdData/CREDSWAP";

		if (XMLDoc->selectNodes((_bstr_t)(const char*)strtmp, &resultList) == S_OK)
		{
			resultList->get_length(&nbNodes);

			if (nbNodes == 0)
			{
				hr = S_FALSE;

				CCString msg ((CCString)"Invalid XML string for ZC \n" + chaineXML);
		
				ICMTHROW(ERR_INVALID_ARGUMENT,tradeId<<":"<<msg);
			}

			hr=resultList->get_item(0, &listItem);

			//Recuperation de la Protect Start Date
			listItem->selectSingleNode(_bstr_t("BindingDate"), &theNode);
			
			if (theNode!=NULL)
			{
				BSTR resultat = NULL;
				theNode->get_text(&resultat);
				_bstr_t ff(resultat,false);
				char * ff1=(char *)ff;

				if (strcmp((const char*)ff1,"") == 0)
					ICMTHROW(ERR_INVALID_ARGUMENT,tradeId<<": ProtectionStartDate not found... "); 
				
				ARM_Date tmpDate (ff1,"YYYYMMDD");
				tmpDate.JulianToStrDate(TProtect_StartDate); 
				Protect_StartDate = tmpDate;
				

				theNode->Release();
				theNode=NULL;
				if (resultat) SysFreeString(resultat);
			}

			//Recuperation de la Protect End Date
			listItem->selectSingleNode(_bstr_t("MatDate"), &theNode);
			
			if (theNode!=NULL)
			{
				BSTR resultat = NULL;
				theNode->get_text(&resultat);
				_bstr_t ff(resultat,false);
				char * ff1=(char *)ff;

				if (strcmp((const char*)ff1,"") == 0)
					ICMTHROW(ERR_INVALID_ARGUMENT,tradeId<<": ProtectionEndDate not found..."); 
				ARM_Date tmpDate (ff1,"YYYYMMDD");
				tmpDate.JulianToStrDate(TProtect_EndDate); 
				Protect_EndDate = tmpDate;
				

				theNode->Release();
				theNode=NULL;
				if (resultat) SysFreeString(resultat);
			}
			
			//Recuperation du type de la courbe
			listItem->selectSingleNode((_bstr_t)"cSmileProp", &theNode);
			if (theNode!=NULL)
			{
				BSTR resultat = NULL;
				theNode->get_text(&resultat);

				_bstr_t ff(resultat,false);
				char * ff1=(char *)ff;
				strcpy(SmilProp,ff1);

				theNode->Release();
				theNode=NULL;
				if (resultat) SysFreeString(resultat);
			}

			
			listItem->selectSingleNode(_bstr_t("DelivDays"), &theNode);
			if (theNode!=NULL)
			{
				BSTR resultat = NULL;
				theNode->get_text(&resultat);

				_bstr_t ff(resultat,false);
				char * ff1=(char *)ff;
				
				sscanf(ff1, "%2dD", &DelivDays);

				theNode->Release();
				theNode=NULL;
				if (resultat) SysFreeString(resultat);
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

			if (strcmp(SmilProp,"") == NULL)
			{
				if (XMLDoc) XMLDoc->Release();
				XMLDoc = NULL;

				if (xmlWCharText)
					delete[] xmlWCharText;
				xmlWCharText = NULL;

				return NULL;
			}

			CvtStrPropInfo(SmilProp,info);

		}

		// ---------------------------------------------------------------------------------
		// Recuperation des infos communes à tous les deals
		// ---------------------------------------------------------------------------------

		Local_XML_TRACE("c:\\temp\\log_arm_credit.txt","Recuperation des infos communes à tous les deals");

 		if (XMLDoc->selectNodes((_bstr_t)("//OPSPEC"), &resultList) == S_OK)
		{

			hr=resultList->get_item(0, &listItem);
			
			if (listItem)
			{
			//Recuperation de la Trade Date
			IntegrationMethod = (int) XML_doubleNodeTreating(listItem,"NumIter");

			if (IntegrationMethod == 0) IntegrationMethod = INTEGRATION_STEP;

			listItem->Release();
			listItem=NULL;
			}
		}

		if (resultList)
		{
			resultList->Release();
			resultList = NULL;
		}

		
	
 		if (XMLDoc->selectNodes((_bstr_t)("//ASSET"), &resultList) == S_OK)
		{

			hr=resultList->get_item(0, &listItem);
			
			//Recuperation de la Trade Date
			StartDate = GetStartDate(listItem);

			//Recuperation de la Maturity Date
			EndDate = GetEndDate(listItem);

			//Recuperation de frequence
			listItem->selectSingleNode(_bstr_t("SCHED_Pay_Freq"), &theNode);
			if (theNode!=NULL)
			{
				BSTR resultat = NULL;
				theNode->get_text(&resultat);

				_bstr_t ff(resultat,false);
				char * ff1=(char *)ff;
				payFreq = FromSummitFreqToARMFreq((const char*) ff1);

				theNode->Release();
				theNode=NULL;
				if (resultat) SysFreeString(resultat);
			}

			//Recuperation de frequence
			listItem->selectSingleNode(_bstr_t("SCHED_Pay_Cal"), &theNode);
			if (theNode!=NULL)
			{
				BSTR resultat = NULL;
				theNode->get_text(&resultat);

				_bstr_t ff(resultat,false);
				PayCal = ff;
				// if (PayCal)
				// 	delete[] PayCal;
				// PayCal = new char[5];
				// strcpy(PayCal,ff1);

				theNode->Release();
				theNode=NULL;
				if (resultat) SysFreeString(resultat);
			}

			// interestRule
			interestRule = GetIntRule(listItem) ;
			if (interestRule==K_UNADJUSTED) adjStartDateRule=K_UNADJUSTED; 
			else adjStartDateRule=K_ADJUSTED ; 

			//Recuperation du basis daycount
			listItem->selectSingleNode(_bstr_t("INTEREST_Basis"), &theNode);
			if (theNode!=NULL)
			{
				BSTR resultat = NULL;
				theNode->get_text(&resultat);

				_bstr_t ff(resultat,false);
				char * ff1=(char *)ff;
				daycount = FromSummitDaycountToARMDaycount((const char*) ff1);

				theNode->Release();
				theNode=NULL;
				if (resultat) SysFreeString(resultat);
			}

			//Recuperation du Accrual basis daycount
			listItem->selectSingleNode(_bstr_t("INT_ACC_Basis"), &theNode);
			if (theNode!=NULL)
			{
				BSTR resultat = NULL;
				theNode->get_text(&resultat);

				_bstr_t ff(resultat,false);
				char * ff1=(char *)ff;
				ACC_daycount = FromSummitDaycountToARMDaycount((const char*) ff1);

				theNode->Release();
				theNode=NULL;
				if (resultat) SysFreeString(resultat);
			}

			// reference Date
			isStubStart=false; 
			listItem->selectSingleNode(_bstr_t((const char *)"STUB_Date1"), &theNode);
			if (theNode!=NULL)
			{
				BSTR resultat = NULL;
				theNode->get_text(&resultat);

				_bstr_t ff(resultat,false);
				char * ff1=(char *)ff;

				if (strcmp((const char*)ff1,"") != 0)
				{
					ARM_Date tmpDate (ff1,"YYYYMMDD");
					startStubDate = tmpDate;
					isStubStart=true; 
				}

				theNode->Release();
				theNode=NULL;
				if (resultat) SysFreeString(resultat);
			}

			// reference Date
			isStubEnd=false ;
			listItem->selectSingleNode(_bstr_t((const char *)"STUB_Date2"), &theNode);
			if (theNode!=NULL)
			{
				BSTR resultat = NULL;
				theNode->get_text(&resultat);

				_bstr_t ff(resultat,false);
				char * ff1=(char *)ff;

				if (strcmp((const char*)ff1,"") != 0)
				{
					ARM_Date tmpDate (ff1,"YYYYMMDD");
					endStubDate = tmpDate;
					isStubEnd=true; 
				}

				theNode->Release();
				theNode=NULL;
				if (resultat) SysFreeString(resultat);
			}


			// Fst Cpn Effective Date
			isFstCpEffect=false; 
			listItem->selectSingleNode(_bstr_t((const char *)"FstCpEffect"), &theNode);
			if (theNode!=NULL)
			{
				BSTR resultat = NULL;
				theNode->get_text(&resultat);

				_bstr_t ff(resultat,false);
				char * ff1=(char *)ff;

				if (strcmp((const char*)ff1,"") != 0)
				{
					isFstCpEffect=true ; 
					ARM_Date tmpDate (ff1,"YYYYMMDD");
					fstCpEffectDate= tmpDate;
				}

				theNode->Release();
				theNode=NULL;
				if (resultat) SysFreeString(resultat);
			}

			// reference Day1
			listItem->selectSingleNode(_bstr_t((const char *)"SCHED_Pay_AnnDay"), &theNode);
			if (theNode!=NULL)
			{
				BSTR resultat = NULL;
				theNode->get_text(&resultat);

				_bstr_t ff(resultat,false);
				char * ff1=(char *)ff;

				refDay = atoi((const char*) ff1);

				theNode->Release();
				theNode=NULL;
				if (resultat) SysFreeString(resultat);
			}

			//Recuperation du stub
			listItem->selectSingleNode(_bstr_t("STUB_StubType"), &theNode);
			if (theNode!=NULL)
			{
				BSTR resultat = NULL;
				theNode->get_text(&resultat);

				_bstr_t ff(resultat,false);
				char * ff1=(char *)ff;
				stub = SummitStub2ARMStub(ff1);
				stub0=stub ;

				theNode->Release();
				theNode=NULL;
				if (resultat) SysFreeString(resultat);
			}

			listItem->selectSingleNode(_bstr_t("Ccy"), &theNode);
			if (theNode!=NULL)
			{
				BSTR resultat = NULL;
				theNode->get_text(&resultat);

				_bstr_t ff(resultat,false);
				ccy = ff; 
				

				theNode->Release();
				theNode=NULL;
				if (resultat) SysFreeString(resultat);
			}


			if (startStubDate<=StartDate || startStubDate>=EndDate) startStubDate=ARM_Date(); 
			if (endStubDate<=StartDate || endStubDate>=EndDate) endStubDate=ARM_Date(); 
			if ( payFreq != K_DAILY) 
				DeduceRefDateAndStub(StartDate,EndDate,startStubDate,endStubDate,payFreq,refDay,ccy,(char*)ReferenceDate,stub);
			else 
			{
				startStubDate=ARM_Date(); 
				endStubDate=ARM_Date(); 
				isFstCpEffect=false; 
			}
			
			FixedPayerAmount= XML_doubleNodeTreating(listItem,"Notional");
			FloatingPayerAmount=FixedPayerAmount;

			CoefAdjust = FixedPayerAmount;

			//Recuperation du PorS
			listItem->selectSingleNode(_bstr_t("PorS"), &theNode);
			if (theNode!=NULL)
			{
				BSTR resultat = NULL;
				theNode->get_text(&resultat);

				_bstr_t ff(resultat,false);
				char * ff1=(char *)ff;
				
				if (!strcmp(ff1,"P"))
					SensDeLaTransaction = 1.;
				else
					SensDeLaTransaction = -1.;

				theNode->Release();
				theNode=NULL;
				if (resultat) SysFreeString(resultat);
			}

			//	Recuperation de INTEREST_Method
			{
				std::string method ; 
				if ( GetStringFromXMLNode(listItem, "INTEREST_Method",method) ) 
				{
					if (method =="") excludeMaturity=EXCLUDE_MATURITY ; 
					else if (method =="STAND") excludeMaturity=EXCLUDE_MATURITY ; 
					else if (method =="INCL") excludeMaturity=INCLUDE_MATURITY; 
					// else if (method =="EFF") excludeMaturity=INCLUDE_MATURITY ; 
					else 
					{
						ICMTHROW(ERR_INVALID_ARGUMENT,tradeId<<": Can't understand INTEREST_Method. "<< method<<" for " << tradeId ); 
					}
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



		
		
		// ---------------------------------------------------------------------------------
		// Recuperation du spread fix
		// ---------------------------------------------------------------------------------
		
		

		if (XMLDoc->selectNodes((_bstr_t)("//ASSET/ProdData/CREDSWAP"), &resultList) == S_OK)
		{
			resultList->get_length(&nbNodes);

			if (nbNodes == 0)
			{
				hr = S_FALSE;
				ICMTHROW(ERR_INVALID_ARGUMENT,tradeId<<": Invalid XML string for getting //ASSET/ProdData/CREDSWAP");
			}

			hr=resultList->get_item(0, &listItem);
			//spread = XML_doubleNodeTreating(listItem,"Formula");
			ParseSpreadFollowTxt(listItem, spread, index_ccy, index_term, index_name,INDEX_spreadonly);
			spread /= 100.;


			listItem->Release();
			listItem=NULL;

			hr=resultList->get_item(0, &listItem);
			listItem->selectSingleNode(_bstr_t("SettleMethod"), &theNode);
			if (theNode!=NULL)
			{
				BSTR resultat = NULL;
				theNode->get_text(&resultat);

				_bstr_t ff(resultat,false);
				char * ff1=(char *)ff;
				
				if (strcmp(ff1,"BIN") != NULL) Binary = CREDIT_DEFAULT_VALUE;

				theNode->Release();
				theNode=NULL;
				if (resultat) SysFreeString(resultat);
			}
			listItem->Release();
			listItem=NULL;

			if (Binary != CREDIT_DEFAULT_VALUE)
			{
				hr=resultList->get_item(0, &listItem);
				Binary = 1. - XML_doubleNodeTreating(listItem,"BinaryAmt");
				listItem->Release();
				listItem=NULL;
			}

		}

		if (resultList)
		{
			resultList->Release();
			resultList = NULL;
		}

		// ---------------------------------------------------------------------------------
		// Recuperation des infos sur les Issuers
		// ---------------------------------------------------------------------------------

		
		Local_XML_TRACE("c:\\temp\\log_arm_credit.txt","Recuperation des infos sur les Issuers");

		if (XMLDoc->selectNodes((_bstr_t)("//ASSET/ProdData/CREDSWAP/CreditEntList/CREDIT[DefpIssuer != \"FTD1\"]"), &resultList0) == S_OK)
			{
				resultList0->get_length(&nbNodes);

				if (nbNodes == 0)
				{	
					hr = S_FALSE;
					ICMTHROW(ERR_INVALID_ARGUMENT,tradeId <<": Invalid XML string for getting ");
				}

				Issuers = new char*[nbNodes];
				nbissuers = nbNodes; 

				for (long indexNode=0 ; indexNode<nbNodes ; indexNode++)
				{
					hr=resultList0->get_item(indexNode, &listItem);

					listItem->selectSingleNode(_bstr_t("DefpIssuer"), &theNode);
					if (theNode!=NULL)
					{
					BSTR resultat = NULL;
					theNode->get_text(&resultat);

					_bstr_t ff(resultat,false);
					char * ff1=(char *)ff;

					Issuers[indexNode] = new char[1000];
					strcpy(Issuers[indexNode],ff1);
		
					theNode->Release();
					theNode=NULL;
					if (resultat) SysFreeString(resultat);
					}

					listItem->Release();
					listItem = NULL;

				}

				if (resultList0)
				{
					resultList0->Release();
					resultList0 = NULL;
				}


			if (XMLDoc->selectNodes((_bstr_t)("//ASSET/ProdData/CREDSWAP/CreditEntList/CREDIT"), &resultList) == S_OK)
			{
				resultList->get_length(&nbNodes);

				if (nbNodes == 0)
				{	
					hr = S_FALSE;
					ICMTHROW(ERR_INVALID_ARGUMENT,tradeId<<": Invalid XML string for getting ");
				}

				Notionals = new double[nbNodes];

				for (long indexNode=0 ; indexNode<nbNodes ; indexNode++)
				{
					hr=resultList->get_item(indexNode, &listItem);

					Notionals[indexNode] = XML_doubleNodeTreating(listItem,"Quantity");

					listItem->Release();
					listItem = NULL;

				}

				if (theNode)
				{
				theNode->Release();
				theNode=NULL;
				}
	
			}

			if (resultList)
			{
			resultList->Release();
			resultList = NULL;
			}

			if (XMLDoc->selectNodes((_bstr_t)("//ASSET/ProdData/CREDSWAP"), &resultList) == S_OK)
			{
				resultList->get_length(&nbNodes);

				if (nbNodes == 0)
				{
				hr = S_FALSE;
				ICMTHROW(ERR_INVALID_ARGUMENT,tradeId<<": Invalid XML string for getting ");
				}

				hr=resultList->get_item(0, &listItem);

				listItem->selectSingleNode(_bstr_t("Tranched"), &theNode);
				if (theNode!=NULL)
				{
				BSTR resultat = NULL;
				theNode->get_text(&resultat);

				_bstr_t ff(resultat,false);
				char * ff1=(char *)ff;

				theNode->Release();
				theNode=NULL;
				if (resultat) SysFreeString(resultat);
				}

				MezzAmount = XML_doubleNodeTreating(listItem,"ProtectStartAmt");
				SubAmount = XML_doubleNodeTreating(listItem,"ProtectEndAmt");

				if (theNode)
				{
				theNode->Release();
				theNode=NULL;
				}

				listItem->Release();
				listItem=NULL;
			}

			if (resultList)
			{
			resultList->Release();
			resultList = NULL;
			}
		}

		if (resultList)
		{
			resultList->Release();
			resultList = NULL;
		}

		
		//------------------------------------------------------------------------------
		//Création de l'instrument financier
		//------------------------------------------------------------------------------

		Local_XML_TRACE("c:\\temp\\log_arm_credit.txt","Création du NTD");

		if (EndDate <= StartDate)
		{
			Local_XML_TRACE("c:\\temp\\log_arm_credit.txt","Error: EndDate <= StartDate");
			ICMTHROW(ERR_INVALID_ARGUMENT,tradeId<<": Error: EndDate <= StartDate");
		}

		if (MezzAmount > SubAmount)
		{
			Local_XML_TRACE("c:\\temp\\log_arm_credit.txt","Error: MezzAmount > SubAmount");
			ICMTHROW(ERR_INVALID_ARGUMENT,tradeId<<": Error: EndDate <= StartDate");
		}

		if (nbissuers == 0)
		{
			Local_XML_TRACE("c:\\temp\\log_arm_credit.txt","Error: MezzAmount > SubAmount");
			ICMTHROW(ERR_INVALID_ARGUMENT,tradeId<<": Error: nbissuers = 0");
		}

		double* Notionals_Modified = new double[nbissuers];
		double MezzAmount_Modified = 0.;
		double SubAmount_Modified = 0.;

		ConversionMnts(nbissuers,Notionals,MezzAmount,SubAmount,
					Notionals_Modified, SubAmount_Modified, MezzAmount_Modified);

		
		CoefAdjust = FixedPayerAmount/MezzAmount_Modified;
		

		int deflegPayFreq (payFreq) ; 
		char deflegRefDate[11]; 
		strcpy(deflegRefDate,ReferenceDate); 
		if (payFreq==K_DAILY) 
		{
			deflegPayFreq =K_QUARTERLY ;
			strcpy(deflegRefDate,"NULL"); 
		} 
		//	SUMMIT: for unknown reason, there is an additional summit call
		//		with ProtectionStart=AsOf that can be > End Date. 
		//		We do not control this call and don't want to throw.
		//
		if (Protect_StartDate>Protect_EndDate) Protect_StartDate=Protect_EndDate; 
		//on check si la default leg est NULL
		if (Protect_StartDate==Protect_EndDate) 
			//	either both= EndDate : RISKY_LEG
			//	or both < EndDate : RISKY_LEG + FROZEN NOT
		{ 
			//	This is RISKY_LEG type 
			//	The schedule is the same as the feeleg. 
			DefaultLegType = qNone_Leg; 
			bool toto (isFstCpEffect) ; 
			if (isFstCpEffect && fstCpEffectDate>=StartDate) toto=false;
			int stub=stub0; 
			ARM_Date startStubDate2(startStubDate); 
			ARM_Date endStubDate2(endStubDate); 
			if (startStubDate2<=StartDate || startStubDate2>=Protect_EndDate) startStubDate2=ARM_Date(); 
			if (endStubDate2<=StartDate || endStubDate2>=Protect_EndDate) endStubDate2=ARM_Date(); 
			DeduceRefDateAndStub(StartDate,Protect_EndDate,
					startStubDate2,
					endStubDate2,
					deflegPayFreq,refDay,ccy,(char*)deflegRefDate,stub);

			DefLeg = new ICM_Leg(StartDate // Protect_StartDate, 
									  Protect_EndDate // Matu, 
									  strcmp(deflegRefDate,"NULL")==0 ? 0 : &ARM_Date(deflegRefDate),
									  toto?&fstCpEffectDate:0,
									  spread,
									  SensDeLaTransaction*MezzAmount_Modified,
									  NULL,
									  NULL,
									  NULL,
									  deflegPayFreq , 
									  daycount,
									  K_ARREARS,
									  K_ADJUSTED,	
									  stub ,
									  ccy,
									  PayCal,
									  DefaultLegType,
									  excludeMaturity,
									  adjStartDateRule,
										NULL, // ICM_Credit_Index* irIndex // = NULL,
									CREDIT_DEFAULT_VALUE, // const double& Binary //= CREDIT_DEFAULT_VALUE,
									ISSUER_UNDEFINE, // const string& name //= ISSUER_UNDEFINE,
									K_NX_NONE, // const int& NXchange //= K_NX_NONE) : ARM_SwapLeg()
									qACCRUED_SETTLED 
									  );
		}
		else
		{
			bool toto (isFstCpEffect) ; 
			if (isFstCpEffect && fstCpEffectDate>=StartDate) toto=false;
			int stub=stub0; 
			ARM_Date startStubDate2(startStubDate); 
			ARM_Date endStubDate2(endStubDate); 
			if (startStubDate2<=Protect_StartDate || startStubDate2>=Protect_EndDate) startStubDate2=ARM_Date(); 
			if (endStubDate2<=Protect_StartDate || endStubDate2>=Protect_EndDate) endStubDate2=ARM_Date(); 
			DeduceRefDateAndStub(StartDate	// Protect_StartDate,
					Protect_EndDate,
					startStubDate2,
					endStubDate2,
					deflegPayFreq,refDay,ccy,(char*)deflegRefDate,stub);
			DefLeg = new ICM_Leg(StartDate,	// Protect_StartDate 
										  Protect_EndDate, 
										  strcmp(deflegRefDate,"NULL")==0 ? 0 : &ARM_Date(deflegRefDate),
										  toto?&fstCpEffectDate:0,
										  spread,
										  SensDeLaTransaction*MezzAmount_Modified,
										  NULL,
										  NULL,
										  NULL,
										  deflegPayFreq , 
										  daycount,
										  K_ARREARS,
										  K_MATUNADJUSTED, // K_ADJUSTED, 
										  stub ,
										  ccy,
										  PayCal,
										  DefaultLegType,
										  excludeMaturity,
										  adjStartDateRule,
										NULL, // ICM_Credit_Index* irIndex = NULL,
									CREDIT_DEFAULT_VALUE, // const double& Binary = CREDIT_DEFAULT_VALUE,
									ISSUER_UNDEFINE, // const string& name = ISSUER_UNDEFINE,
									K_NX_NONE, // const int& NXchange = K_NX_NONE) : ARM_SwapLeg()
									qACCRUED_SETTLED 
										  );
		}

		bool toto (isFstCpEffect) ; 
		if (isFstCpEffect && fstCpEffectDate>=StartDate) toto=false;
		FeeLeg = new ICM_Leg(StartDate, 
		  				 EndDate, 
						 strcmp(ReferenceDate,"NULL")==0 ? 0 : &ARM_Date(ReferenceDate),
						 toto?&fstCpEffectDate:0,
						 spread,
						 SensDeLaTransaction*MezzAmount_Modified,
						 NULL,
						 NULL,
						 NULL,
						 payFreq , 
						 daycount,
						 K_ARREARS,
						 interestRule ,
						 stub ,
						 ccy,
						 PayCal,
						 FeeLegType,
						 excludeMaturity,
						 adjStartDateRule,
						NULL, // ICM_Credit_Index* irIndex = NULL,
						CREDIT_DEFAULT_VALUE, // const double& Binary = CREDIT_DEFAULT_VALUE,
						ISSUER_UNDEFINE, // const string& name = ISSUER_UNDEFINE,
						K_NX_NONE ,// const int& NXchange = K_NX_NONE) : ARM_SwapLeg()
						qACCRUED_SETTLED 
						 );

	   	//	This is the Floating Coupon
		if (INDEX_spreadonly==false)
		{
			ARM_ReferenceValue Ref((double)qCPN_FIXED_AND_FWD);
			FeeLeg->SetFwdCalcTypes(&Ref);
			FeeLeg->SetCreditLegType(qSwapLeg);
			int indexType = FromStringToIndexType((char*)index_term.c_str(),(char*)index_name.c_str());
			ARM_IRIndex index((ARM_INDEX_TYPE)indexType, payFreq,payFreq, &ARM_Currency((char*)index_ccy.c_str()));
			FeeLeg->SetIRIndex(&index);

			ARM_ReferenceValue* Fixings = ARMLOCAL_XML_EVENT(chaineXML,"FRC");
			if (Fixings)
			{	for (int ila=0;ila<Fixings->GetDiscreteDates()->GetSize();ila++)
				{	ARM_Date date = (ARM_Date)Fixings->GetDiscreteDates()->Elt(ila);
					double value = 100.*Fixings->GetDiscreteValues()->Elt(ila);
					// FeeLeg->GetCreditInfos()->setPastFixing(date,value) ;}
					FeeLeg->setPastFixing(date,value) ;}

				if (Fixings) delete Fixings;Fixings=NULL;
			}
		}
		// }

		ICM_Collateral collateral(nbissuers,Issuers,Notionals_Modified);
		DefLeg->SetCreditLag(DelivDays); 
		ICM_Cds	cds(FeeLeg,DefLeg,DelivDays);	// to be fixed some day... 
		if (FeeLeg) {delete FeeLeg;FeeLeg=NULL;}
		if (DefLeg) {delete DefLeg;DefLeg=NULL;}

		ICM_Mez* cdo= new ICM_Mez(&cds,SubAmount_Modified,&collateral,Binary);

		if (IntegrationMethod != CREDIT_DEFAULT_VALUE)
		{
				ICM_SummitContainer SummitContainer;
				SummitContainer.itsIntegrationMethod = IntegrationMethod;
				cdo->SetSummitContainer(&SummitContainer);
		}
		
		if (info)
		{
			cdo->SetProportions(info);
			delete info;
		}

		cdo->SetTradedCoef(CoefAdjust);
		cdo->SetPorS(SensDeLaTransaction);
 		if (Protect_EndDate<EndDate) 
		{
			ARM_Vector dates(2), types(2), specdates(2) ; 
			dates.Elt(0) = StartDate.GetJulian() ; 
			dates.Elt(1) = Protect_EndDate.GetJulian() ; 
			types.Elt(0) = qStdRisk; 
			types.Elt(1) = qRiskToDate; 
			specdates.Elt(0) = ARM_Date(01,01,2000).GetJulian(); 
			specdates.Elt(1) = Protect_EndDate.GetJulian() ; 
			ARM_ReferenceValue refcpn((ARM_Vector*)dates.Clone(),(ARM_Vector*)specdates.Clone()); refcpn.SetCalcMethod(K_STEPUP_LEFT); 
			ARM_ReferenceValue refstyle((ARM_Vector*)dates.Clone(),(ARM_Vector*)types.Clone()); refstyle.SetCalcMethod(K_STEPUP_LEFT); 
			cdo->GetFeeLeg()->SetRefRiskyDate(&refcpn,&refstyle);
		}
  
		// if (ccy) delete ccy;

		Local_XML_TRACE("c:\\temp\\log_arm_credit.txt","Fin Création du CDO ");

		for (int i = 0; i<nbissuers; i++)
		{
			if (Issuers[i]) delete[] Issuers[i];
			Issuers[i] = NULL;
		}

		if (Issuers) delete[] Issuers;
		if (Notionals) delete[] Notionals;		
		if (Notionals_Modified) delete[] Notionals_Modified;
		if (Rate) delete[] Rate;		


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

		if (XMLDoc)	{ XMLDoc->Release();}
		XMLDoc = NULL;

		if (xmlWCharText)
			delete[] xmlWCharText;
		xmlWCharText = NULL;

		return (cdo);
	}
	catch(Exception& x)
	{
		if (theNode) theNode->Release();
		if (listItem) listItem->Release();
		if (resultList) resultList->Release();
		if (XMLDoc) XMLDoc->Release();


		if (xmlWCharText)
			delete[] xmlWCharText;
		xmlWCharText = NULL;

		hr = S_FALSE;
		char Tmp[1000];
		x.GetErrorMessage(Tmp);
		ICMTHROW(ERR_INVALID_ARGUMENT,tradeId<<": "<<Tmp);
	}
	catch(...)
	{		
		if (theNode) theNode->Release();
		if (listItem) listItem->Release();
		if (resultList) resultList->Release();
		if (XMLDoc) XMLDoc->Release();

		if (xmlWCharText)
			delete[] xmlWCharText;
		xmlWCharText = NULL;

		hr = S_FALSE;
		ICMTHROW(ERR_INVALID_ARGUMENT,tradeId<<": Error in XML parsing for getting Tranche");
	}

	return NULL;
}
**/ 


// ------------------------------------------------------------------------------
// Amount Conversion to avoid rounding error on LossUnit
// ------------------------------------------------------------------------------

void ConversionMnts(int nbnames, double* not_init,double sub_init, double end_init,
					double*& not_out, double& sub_out, double& tranche_out)
{
	int i=0;
	double size_ptf_init = 0.;
	pgcd t;
	
// FIXMEFRED: mig.vc8 (30/05/2007 17:47:54):std::_cpp_max is not right
	for (i=0;i<nbnames;i++) size_ptf_init += _cpp_max(0.,not_init[i]);

	size_ptf_init = t.round(size_ptf_init);

	double startprotect = sub_init/size_ptf_init;
	double endprotect = end_init/size_ptf_init;
	not_out[0] = 10000000.;

	for (i=1;i<nbnames;i++)	not_out[i] = t.round(not_out[0]*not_init[i]/not_init[0]);

	double size_ptf_out = 0.;
// FIXMEFRED: mig.vc8 (30/05/2007 17:47:54):std::_cpp_max is not right
	for (i=0;i<nbnames;i++) size_ptf_out += _cpp_max(0.,not_out[i]);

	size_ptf_out = t.round(size_ptf_out);
		
	sub_out = t.round(size_ptf_out*startprotect);
	tranche_out = t.round(size_ptf_out*(endprotect-startprotect));
}

// ------------------------------------------------------------------------------
// Création d'un instrument option pour un pricer définissant le sous-jacent
// ------------------------------------------------------------------------------

ICM_Option* ARMLOCAL_XML_CDS_OPTION(const char* chaineXML, 
									   CCString& bookName, 
									   CCString& custId, 
									   CCString& dealId,
									   ICM_Pricer* UnderlyingPricer)
{
	VARIANT_BOOL bOK;
	HRESULT hr;

	MSXML2::IXMLDOMDocument *XMLDoc = NULL;
	MSXML2::IXMLDOMNodeList * resultList = NULL, * resultList2 = NULL;
	MSXML2::IXMLDOMNode * listItem = NULL, * theNode = NULL, *item = NULL;
	MSXML2::IXMLDOMNode * floatingNode = NULL;
	long nbNodes;

	ARM_Date startDate;
	ARM_Date endDate;
	int rcvPay;
	double strike;
	char index[20];
	char term[20];
	char ccy[20];
    ARM_Date maturity;
	//int liborType;
    double swapYearTerm = -1000000.0;
	int resetFreq = K_DEF_FREQ;
	int payFreq = K_DEF_FREQ;

	double dPorS = 1.0;
	int PutOrCall = 0;

	//ARM_Currency* discountCcy = NULL;
	//ARM_ReferenceValue* notional = NULL;

	ICM_Option* Option = NULL;
	std::string tradeId("trade?"); 
	try
	{
		hr = CoInitialize(NULL); 
//		SUCCEEDED(hr) ? 0 : throw hr;
		hr = CoCreateInstance(__uuidof(MSXML2::DOMDocument30), NULL, CLSCTX_INPROC_SERVER, __uuidof(MSXML2::IXMLDOMDocument), (void**)&XMLDoc);
		SUCCEEDED(hr) ? 0 : throw hr;

		//JLA _bstr_t tmpChaine = chaineXML;
		_bstr_t tmpChaine; 
		VariantTools::convert(std::string(chaineXML),tmpChaine); 



		XMLDoc->loadXML(tmpChaine, &bOK);
				
		#ifdef _DEBUG
		CCString chtmp = "c:\\temp\\AAA_CDSOPT_OUTPUT.xml";
		VARIANT v_chtmp;
		CCString2VARIANT (chtmp, &v_chtmp);
		XMLDoc->save(v_chtmp);
		#endif

	}
	catch(...)
	{
		if (XMLDoc) XMLDoc->Release();
		hr = S_FALSE;

		CCString msg((CCString)"Pb in creating XML document for getting Option \n" + chaineXML);

		throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
						 (const char *) msg);
	}

	try
	{
		// Récupération du customer id
		if (XMLDoc->selectNodes((_bstr_t)("Response/SWAPTION/Env/ENV"), &resultList) == S_OK)
		{
			resultList->get_length(&nbNodes);

			if (nbNodes != 1)
			{
				hr = S_FALSE;

				CCString msg((CCString)"Invalid XML string for getting Swaption \n" + chaineXML);

				ICMTHROW(ERR_INVALID_ARGUMENT,tradeId<<":"<<(const char *) msg);
			}

			hr=resultList->get_item(0, &listItem);

			custId = GetCustumerId(listItem);

			dealId = GetDealId(listItem);
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

		// Recuperation de la partie optionnelle
		if (XMLDoc->selectNodes((_bstr_t)("Response/SWAPTION/Option/OPTION"), &resultList) == S_OK)
		{
			resultList->get_length(&nbNodes);

			if (nbNodes != 1)
			{
				hr = S_FALSE;

				CCString msg((CCString)"Invalid XML string for getting Option \n" + chaineXML);

				throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
								 (char*) msg);
			}

			hr=resultList->get_item(0, &listItem);

			maturity = GetDateFromXMLNode(listItem,"ExpDate");

			strike = XML_doubleNodeTreating(listItem,"Strike") * 100.;

			dPorS = GetPorS(listItem);

			// ATTENTION c'est ici l'inverse des cap/floor
			listItem->selectSingleNode(_bstr_t("PorC"), &theNode);
			if (theNode!=NULL)
			{
				BSTR resultat = NULL;
				theNode->get_text(&resultat);

				_bstr_t ff(resultat,false);
				char * ff1=(char *)ff;

				if (strcmp(ff1,"C") == 0)
				{ rcvPay = K_PAY; PutOrCall = K_CALL; } 
				else
				{ rcvPay = K_RCV; PutOrCall = K_PUT; } 

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

		// Recuperation du book
		if (XMLDoc->selectNodes((_bstr_t)("Response/SWAPTION/Env/ENV"), &resultList) == S_OK)
		{
			resultList->get_length(&nbNodes);

			if (nbNodes != 1)
			{
				hr = S_FALSE;

				CCString msg((CCString)"Invalid XML string for getting OPtion \n" + chaineXML);

				throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
								 (const char *) msg);
			}

			hr=resultList->get_item(0, &listItem);

			bookName = GetBook(listItem);
			listItem->Release();
			listItem=NULL;
		}

		if (resultList)
		{
			resultList->Release();
			resultList = NULL;
		}

		// Recuperation de la partie swap
		if (XMLDoc->selectNodes((_bstr_t)("Response/SWAPTION/Assets/ASSET"), &resultList) == S_OK)
		{
			resultList->get_length(&nbNodes);

			if (nbNodes != 2)
			{
				hr = S_FALSE;

				CCString msg((CCString)"Invalid XML string for getting Option \n" + chaineXML);

				throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
								 (char*) msg);
			}

			for (long indexNode=0 ; indexNode<nbNodes ; indexNode++)
			{
				hr=resultList->get_item(indexNode, &listItem);

				listItem->selectSingleNode(_bstr_t("INTEREST_FixFloat"), &theNode);
				if (theNode!=NULL)
				{
					BSTR resultat = NULL;
					theNode->get_text(&resultat);

					_bstr_t ff(resultat,false);
					char * ff1=(char *)ff;

					if (strcmp(ff1,"FLO") == 0)
					{
						floatingNode = listItem;
						indexNode = nbNodes;
					}
					else
					{
						listItem->Release();
						listItem=NULL;
					}

					theNode->Release();
					theNode=NULL;
					if (resultat) SysFreeString(resultat);
				}
			}

			startDate = GetStartDate(listItem);

			endDate = GetEndDate(listItem);

			listItem->selectSingleNode(_bstr_t("INTEREST_dmIndex"), &theNode);
			if (theNode!=NULL)
			{
				BSTR resultat = NULL;
				theNode->get_text(&resultat);

				_bstr_t ff(resultat,false);
				char * ff1=(char *)ff;

				strcpy(index,ff1);

				theNode->Release();
				theNode=NULL;
				if (resultat) SysFreeString(resultat);
			}

			listItem->selectSingleNode(_bstr_t("INTEREST_Term"), &theNode);
			if (theNode!=NULL)
			{
				BSTR resultat = NULL;
				theNode->get_text(&resultat);

				_bstr_t ff(resultat,false);
				char * ff1=(char *)ff;

				strcpy(term,ff1);

				theNode->Release();
				theNode=NULL;
				if (resultat) SysFreeString(resultat);
			}

			//liborType = FromStringToIndexType(term,index);

			listItem->selectSingleNode(_bstr_t("INTEREST_Ccy"), &theNode);
			if (theNode!=NULL)
			{
				BSTR resultat = NULL;
				theNode->get_text(&resultat);

				_bstr_t ff(resultat,false);
				char * ff1=(char *)ff;

				strcpy(ccy,ff1);

				theNode->Release();
				theNode=NULL;
				if (resultat) SysFreeString(resultat);
			}

			//discountCcy = new ARM_Currency(ccy);

			double dNotional = XML_doubleNodeTreating(listItem,"Notional");

			// notional = new ARM_ReferenceValue(dNotional);

			listItem->Release();
			listItem=NULL;

			//unused
			//startDate 
			//endDate 
			//K_EUROPEAN 
			//rcvPay

			ICM_Cds* Cdsunderlying = ARMLOCAL_XML_CDS_SIMPLIFIED(chaineXML,bookName,custId,dealId);

			double FixedRate = Cdsunderlying->GetFeeLeg()->GetFixedRate();
			strike = 10000 * FixedRate/100.;

			Option = new ICM_Option(Cdsunderlying->GetEndDateNA(),
									maturity,
									strike,
									PutOrCall,
									(qDEF_MAT)1);

			//Option->SetExpiry(maturity);

			if (Cdsunderlying)
				delete Cdsunderlying;
			Cdsunderlying = NULL;

			//newSwaption->SetPorS(dPorS);

			//if (discountCcy)
			//	delete discountCcy;
			//discountCcy = NULL;

//			if (notional)
//				delete notional;
//			notional = NULL;
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

		if (listItem)
		{
			listItem->Release();
			listItem = NULL;
		}

		if (XMLDoc)
		{
			XMLDoc->Release();
			XMLDoc = NULL;
		}

		return Option;
	}

	catch(...)
	{		
		if (theNode) theNode->Release();
		if (listItem) listItem->Release();
		if (resultList) resultList->Release();
		if (XMLDoc) XMLDoc->Release();

		hr = S_FALSE;

		CCString msg((CCString)"Error in XML parsing for getting Swaption \n" + chaineXML);

		throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
						 (const char *) msg);
	}

}


// -----------------------------------------------
// Generic Tranche Parsing for Credit Default Swap
// -----------------------------------------------

ICM_Cds* ARMLOCAL_XML_CDS_SIMPLIFIED(const char* chaineXML,CCString& bookId, CCString& custId, CCString& dealId)
{
	VARIANT_BOOL bOK;
	HRESULT hr;

	MSXML2::IXMLDOMDocument *XMLDoc = NULL;
	MSXML2::IXMLDOMNodeList * resultList = NULL;
	MSXML2::IXMLDOMNodeList * resultList0 = NULL;
	MSXML2::IXMLDOMNode * listItem = NULL, * theNode = NULL;
	long nbNodes = 0;

	ARM_Date StartDate;
	ARM_Date EndDate;
	double spread = 0.;

	char ReferenceDate[11];
	int payFreq = K_DEF_FREQ;
	int daycount = KACTUAL_360;
	int ACC_daycount = KACTUAL_360;
	int DelivDays =0 ;

	double FixedPayerAmount=0.;
	double FloatingPayerAmount=0.;

	// ARM_Currency* ccy = NULL;
	std::string ccy ;
	int stub = 0;

	double MezzAmount = 0.;
	double SubAmount = 0.;

	int nbissuers = 0;
	char** Issuers = NULL; 
	double* Notionals = NULL;
	double* Rate = NULL;
	double CoefAdjust = 0.;
	// char* PayCal = NULL;
	std::string PayCal ;

	double SensDeLaTransaction = 1.;
	double Binary = 0.;

	//ARM_Currency* discountCcy = NULL;
	//ARM_ReferenceValue* notional = NULL;

	int k =0;	

	wchar_t * xmlWCharText = NULL;
	std::string tradeId; // for report

	try
	{
		hr = CoInitialize(NULL); 
//		SUCCEEDED(hr) ? 0 : throw hr;
		hr = CoCreateInstance(__uuidof(MSXML2::DOMDocument30), NULL, CLSCTX_INPROC_SERVER, __uuidof(MSXML2::IXMLDOMDocument), (void**)&XMLDoc);
		SUCCEEDED(hr) ? 0 : throw hr;

		xmlWCharText = constchar2wchar(chaineXML);

		Local_XML_TRACE("c:\\temp\\log_arm_credit.txt","Lecture XML pour Matrice de Correlation");
		XMLDoc->loadXML((_bstr_t)xmlWCharText, &bOK);

	}
	catch(Exception& )
	{
		if (XMLDoc) XMLDoc->Release();
		hr = S_FALSE;
		throw; 


		if (xmlWCharText)
			delete[] xmlWCharText;
		xmlWCharText = NULL;
	}
	catch(...)
	{
		if (XMLDoc) XMLDoc->Release();
		hr = S_FALSE;

		if (xmlWCharText)
			delete[] xmlWCharText;
		xmlWCharText = NULL;

		throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"Pb in creating XML document for getting Tranches informations");
	}

	try
	{
		{
			// JLA. Retrieve 
			MSXML2::IXMLDOMNode* aNode = NULL;
			std::string chemin = "//ENV";
			XMLDoc->selectSingleNode((_bstr_t)chemin.c_str(), &aNode);
			custId= GetStringFromXMLNode(aNode, "Cust").c_str();
			dealId = GetStringFromXMLNode(aNode, "DealId").c_str();
			bookId = GetStringFromXMLNode(aNode, "Book").c_str();
			tradeId = GetStringFromXMLNode(aNode, "TradeId") ;
			aNode->Release();
		}

		// ---------------------------------------------------------------------------------
		// Recuperation des infos communes à tous les deals
		// ---------------------------------------------------------------------------------

		Local_XML_TRACE("c:\\temp\\log_arm_credit.txt","Recuperation des infos communes à tous les deals");

 		if (XMLDoc->selectNodes((_bstr_t)("//CREDSWAP"), &resultList) == S_OK)
		{

			hr=resultList->get_item(0, &listItem);
			
			//Recuperation de la Trade Date
			StartDate = GetDateFromXMLNode(listItem,"BindingDate");

			//Recuperation de la Maturity Date
			EndDate = GetEndDate(listItem);

			//Recuperation de frequence
			payFreq = FromSummitFreqToARMFreq("Q");

			//Recuperation de frequence
			//if (PayCal)
			//	delete[] PayCal;
			//PayCal = new char[5];
			//strcpy(PayCal,"EUR");
			PayCal="EUR"; 

			//Recuperation du basis daycount
			daycount = FromSummitDaycountToARMDaycount("A360");

			//Recuperation du Accrual basis daycount
			ACC_daycount = FromSummitDaycountToARMDaycount("A360");

			//Recuperation de la Reference Date
			strcpy(ReferenceDate,"NULL");

			//Recuperation du stub
			stub = SummitStub2ARMStub("SS");

			FixedPayerAmount= XML_doubleNodeTreating(listItem,"Notional");
			FloatingPayerAmount=FixedPayerAmount;

			CoefAdjust = FixedPayerAmount;

			//Recuperation du PorS
			SensDeLaTransaction = 1.;

			// ccy = new ARM_Currency("EUR");
			ccy="EUR"; 

			spread = XML_doubleNodeTreating(listItem,"Formula");
			spread /= 100.;

			Binary = CREDIT_DEFAULT_VALUE;

			
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

		Local_XML_TRACE("c:\\temp\\log_arm_credit.txt","Recuperation des infos sur les Issuers");

			if (XMLDoc->selectNodes((_bstr_t)("//CREDSWAP/CreditEntList/CREDIT[DefpIssuer != \"FTD1\"]"), &resultList0) == S_OK)
			{
				resultList0->get_length(&nbNodes);

				if (nbNodes == 0)
				{	
					hr = S_FALSE;
					throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT, "Invalid XML string for getting ");
				}

				Issuers = new char*[nbNodes];
				nbissuers = nbNodes; 

				for (long indexNode=0 ; indexNode<nbNodes ; indexNode++)
				{
					hr=resultList0->get_item(indexNode, &listItem);

					listItem->selectSingleNode(_bstr_t("DefpIssuer"), &theNode);
					if (theNode!=NULL)
					{
					BSTR resultat = NULL;
					theNode->get_text(&resultat);

					_bstr_t ff(resultat,false);
					char * ff1=(char *)ff;

					Issuers[indexNode] = new char[1000];
					strcpy(Issuers[indexNode],ff1);
		
					theNode->Release();
					theNode=NULL;
					if (resultat) SysFreeString(resultat);
					}

					listItem->Release();
					listItem = NULL;

				}

				if (resultList0)
				{
					resultList0->Release();
					resultList0 = NULL;
				}


			if (XMLDoc->selectNodes((_bstr_t)("//CREDSWAP/CreditEntList/CREDIT"), &resultList) == S_OK)
			{
				resultList->get_length(&nbNodes);

				if (nbNodes == 0)
				{	
					hr = S_FALSE;
					throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT, "Invalid XML string for getting ");
				}

				Notionals = new double[nbNodes];

				for (long indexNode=0 ; indexNode<nbNodes ; indexNode++)
				{
					hr=resultList->get_item(indexNode, &listItem);

					Notionals[indexNode] = XML_doubleNodeTreating(listItem,"Quantity");

					listItem->Release();
					listItem = NULL;

				}

				if (theNode)
				{
				theNode->Release();
				theNode=NULL;
				}
	
			}

			if (resultList)
			{
			resultList->Release();
			resultList = NULL;
			}

			if (XMLDoc->selectNodes((_bstr_t)("//CREDSWAP"), &resultList) == S_OK)
			{
				resultList->get_length(&nbNodes);

				if (nbNodes == 0)
				{
				hr = S_FALSE;
				throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT, "Invalid XML string for getting ");
				}

				hr=resultList->get_item(0, &listItem);

				listItem->selectSingleNode(_bstr_t("Tranched"), &theNode);
				if (theNode!=NULL)
				{
					BSTR resultat = NULL;
					theNode->get_text(&resultat);

					_bstr_t ff(resultat,false);
					char * ff1=(char *)ff;

					theNode->Release();
					theNode=NULL;
					if (resultat) SysFreeString(resultat);
				}

				MezzAmount = XML_doubleNodeTreating(listItem,"ProtectStartAmt");
				SubAmount = XML_doubleNodeTreating(listItem,"ProtectEndAmt");

				if (theNode)
				{
				theNode->Release();
				theNode=NULL;
				}

				listItem->Release();
				listItem=NULL;
			}

			if (resultList)
			{
			resultList->Release();
			resultList = NULL;
			}
		}

		if (resultList)
		{
			resultList->Release();
			resultList = NULL;
		}

		//------------------------------------------------------------------------------
		//Création de l'instrument financier
		//------------------------------------------------------------------------------

		Local_XML_TRACE("c:\\temp\\log_arm_credit.txt","Création du CDS");

		if (EndDate <= StartDate)
		{
			Local_XML_TRACE("c:\\temp\\log_arm_credit.txt","Error: EndDate <= StartDate");
			throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT, "Error: EndDate <= StartDate");
		}

		if (MezzAmount > SubAmount)
		{
			Local_XML_TRACE("c:\\temp\\log_arm_credit.txt","Error: MezzAmount > SubAmount");
			throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT, "Error: EndDate <= StartDate");
		}

		if (nbissuers == 0)
		{
			Local_XML_TRACE("c:\\temp\\log_arm_credit.txt","Error: MezzAmount > SubAmount");
			throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT, "Error: nbissuers = 0");
		}

		double* Notionals_Modified = new double[nbissuers];
		double MezzAmount_Modified = 1000000.;
		
		CoefAdjust = FixedPayerAmount/MezzAmount_Modified;

		ICM_Cds* cds = new ICM_Cds(StartDate,
								   	EndDate,
									strcmp(ReferenceDate,"NULL")==0 ? 0 : &ARM_Date(ReferenceDate),
									0,
									StartDate,
								   	EndDate,
									spread,
									MezzAmount_Modified,MezzAmount_Modified,// 0,0,
									payFreq,
									daycount,
									qACCRUED_SETTLED,
									ccy,
									stub,
									DEFAULT_CREDIT_LAG, 
									DEFAULT_FRQ_DEFLEG,
									K_ADJUSTED, 
									INCLUDE_MATURITY, 
									K_ADJUSTED, 
									PayCal,
									qRunning_Leg,
									qStandart_Recovery_Leg,
									ISSUER_UNDEFINE, // const string& name = ISSUER_UNDEFINE ,
									CREDIT_DEFAULT_VALUE // const double& Binary = CREDIT_DEFAULT_VALUE );
									);

		cds->SetTradedCoef(CoefAdjust);
		cds->SetPorS(SensDeLaTransaction);

		//if (ccy) delete ccy;

		Local_XML_TRACE("c:\\temp\\log_arm_credit.txt","Fin Création du NTD");

		for (int i = 0; i<nbissuers; i++)
		{
			if (Issuers[i]) delete[] Issuers[i];
			Issuers[i] = NULL;
		}

		if (Issuers) delete[] Issuers;
		if (Notionals) delete[] Notionals;		
		if (Notionals_Modified) delete[] Notionals_Modified;
		if (Rate) delete[] Rate;		
		//if (PayCal) delete[] PayCal;

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

		if (XMLDoc)	{ XMLDoc->Release();}
		XMLDoc = NULL;

		if (xmlWCharText)
			delete[] xmlWCharText;
		xmlWCharText = NULL;

		return (cds);
	}
	catch(Exception& )
	{
		if (theNode) theNode->Release();
		if (listItem) listItem->Release();
		if (resultList) resultList->Release();
		if (XMLDoc) XMLDoc->Release();


		if (xmlWCharText)
			delete[] xmlWCharText;
		xmlWCharText = NULL;

		hr = S_FALSE;
		throw; 
	}
	catch(...)
	{		
		if (theNode) theNode->Release();
		if (listItem) listItem->Release();
		if (resultList) resultList->Release();
		if (XMLDoc) XMLDoc->Release();

		if (xmlWCharText)
			delete[] xmlWCharText;
		xmlWCharText = NULL;

		hr = S_FALSE;
		throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
         "Error in XML parsing for getting Tranche");
	}

	return NULL;

}


// ------------------------------------
// Generic Tranche Parsing for Cds 
// ------------------------------------

ICM_Cds* ARMLOCAL_XML_CDS(const char* chaineXML, CCString& bookId, CCString& custId, CCString& dealId,
							ARM_ReferenceValue* fees)
{
	VARIANT_BOOL bOK;
	HRESULT hr;

	MSXML2::IXMLDOMDocument *XMLDoc = NULL;
	MSXML2::IXMLDOMNodeList * resultList = NULL;
	MSXML2::IXMLDOMNodeList * resultList0 = NULL;
	MSXML2::IXMLDOMNode * listItem = NULL, * theNode = NULL;
	long nbNodes = 0;

	ARM_Date StartDate;
	ARM_Date EndDate;
	ARM_Date fstCpEffectDate ;
	bool isFstCpEffect(false); 
	double spread = 0.;

	char ReferenceDate[11];
	int payFreq = K_DEF_FREQ;
	int daycount = KACTUAL_360;
	int ACC_daycount = KACTUAL_360;
	int interestRule, adjStartDateRule ;

	double FixedPayerAmount=0.;
	double FloatingPayerAmount=0.;

	std::string ccy; 
	// ARM_Currency* ccy = NULL;

	double Binary = 0.;

	double* Rate = NULL;
	// char* PayCal = NULL;
	std::string PayCal; 

	double SensDeLaTransaction = 1.;

	// ARM_Currency* discountCcy = NULL;
	// ARM_ReferenceValue* notional = NULL;
	std::auto_ptr<ARM_ReferenceValue> notional ;
	ARM_Date startStubDate;
	ARM_Date endStubDate;
	int refDay = 0;
	int stub = K_SHORTSTART;
	int DelivDays = 0;	
	int k =0;	

	wchar_t * xmlWCharText = NULL;
	std::string tradeId; 
	bool excludeMaturity=EXCLUDE_MATURITY; 

	try
	{
		hr = CoInitialize(NULL); 
//		SUCCEEDED(hr) ? 0 : throw hr;
		hr = CoCreateInstance(__uuidof(MSXML2::DOMDocument30), NULL, CLSCTX_INPROC_SERVER, __uuidof(MSXML2::IXMLDOMDocument), (void**)&XMLDoc);
		SUCCEEDED(hr) ? 0 : throw hr;

		xmlWCharText = constchar2wchar(chaineXML);

		Local_XML_TRACE("c:\\temp\\log_arm_credit.txt","Lecture XML pour Matrice de Correlation");
		XMLDoc->loadXML((_bstr_t)xmlWCharText, &bOK);

		#ifdef _DEBUG
		CCString chtmp = "c:\\temp\\ARMLOCAL_XML_CDS.xml";
		VARIANT v_chtmp;
		CCString2VARIANT (chtmp, &v_chtmp);
		XMLDoc->save(v_chtmp);
		#endif // _DEBUG

	}
	catch(Exception& )
	{
		if (XMLDoc) XMLDoc->Release();
		hr = S_FALSE;
		throw; 


		if (xmlWCharText)
			delete[] xmlWCharText;
		xmlWCharText = NULL;
	}
	catch(...)
	{
		if (XMLDoc) XMLDoc->Release();
		hr = S_FALSE;

		if (xmlWCharText)
			delete[] xmlWCharText;
		xmlWCharText = NULL;

		throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"Pb in creating XML document for getting Tranches informations");
	}

	try
	{
		{
			// JLA. Retrieve 
			MSXML2::IXMLDOMNode* aNode = NULL;
			std::string chemin = "Response/EXOTIC/Env/ENV";
			XMLDoc->selectSingleNode((_bstr_t)chemin.c_str(), &aNode);
			if (aNode)
			{
				custId= GetStringFromXMLNode(aNode, "Cust").c_str();
				dealId = GetStringFromXMLNode(aNode, "DealId").c_str();
				bookId = GetStringFromXMLNode(aNode, "Book").c_str();
				tradeId = GetStringFromXMLNode(aNode, "TradeId") ;
				aNode->Release();
			}
		}


		// ---------------------------------------------------------------------------------
		// Recuperation des infos communes à tous les deals
		// ---------------------------------------------------------------------------------

		Local_XML_TRACE("c:\\temp\\log_arm_credit.txt","Recuperation des infos communes à tous les deals");

		if (fees) {
			// JLA. Retrieve 
			MSXML2::IXMLDOMNode* aNode = NULL;
			std::string chemin = "//ASSET";
			XMLDoc->selectSingleNode((_bstr_t)chemin.c_str(), &aNode);
			ARM_ReferenceValue* item= 0 ;
			if (aNode) 
			{ 
				item =GetPrimes(aNode) ; 
				if (item) *fees=*item; 
			}
			if (item) delete item; item=0; 
			aNode->Release(); 
		}


 		if (XMLDoc->selectNodes((_bstr_t)("//ASSET"), &resultList) == S_OK)
		{

			hr=resultList->get_item(0, &listItem);

			
			//Recuperation de la Trade Date
			StartDate = GetStartDate(listItem);

			//Recuperation de la Maturity Date
			EndDate = GetEndDate(listItem);

			// notional reference value. 
			notional=auto_ptr<ARM_ReferenceValue>(GetNotional(listItem))  ;

			//Recuperation de frequence
			listItem->selectSingleNode(_bstr_t("SCHED_Pay_Freq"), &theNode);
			if (theNode!=NULL)
			{
				BSTR resultat = NULL;
				theNode->get_text(&resultat);

				_bstr_t ff(resultat,false);
				char * ff1=(char *)ff;
				payFreq = FromSummitFreqToARMFreq((const char*) ff1);

				theNode->Release();
				theNode=NULL;
				if (resultat) SysFreeString(resultat);
			}

			//Recuperation de frequence
			listItem->selectSingleNode(_bstr_t("SCHED_Pay_Cal"), &theNode);
			if (theNode!=NULL)
			{
				BSTR resultat = NULL;
				theNode->get_text(&resultat);

				_bstr_t ff(resultat,false);
				PayCal=  ff;
				//if (PayCal)
				//	delete[] PayCal;
				//PayCal = new char[5];
				// strcpy(PayCal,ff1);


				theNode->Release();
				theNode=NULL;
				if (resultat) SysFreeString(resultat);
			}

			//Recuperation du basis daycount
			listItem->selectSingleNode(_bstr_t("INTEREST_Basis"), &theNode);
			if (theNode!=NULL)
			{
				BSTR resultat = NULL;
				theNode->get_text(&resultat);

				_bstr_t ff(resultat,false);
				char * ff1=(char *)ff;
				daycount = FromSummitDaycountToARMDaycount((const char*) ff1);

				theNode->Release();
				theNode=NULL;
				if (resultat) SysFreeString(resultat);
			}

			// interestRule
			interestRule = GetIntRule(listItem) ;
			if (interestRule==K_UNADJUSTED) adjStartDateRule=K_UNADJUSTED; 
			else adjStartDateRule=K_ADJUSTED ; 

			//Recuperation du Accrual basis daycount
			listItem->selectSingleNode(_bstr_t("INT_ACC_Basis"), &theNode);
			if (theNode!=NULL)
			{
				BSTR resultat = NULL;
				theNode->get_text(&resultat);

				_bstr_t ff(resultat,false);
				char * ff1=(char *)ff;
				ACC_daycount = FromSummitDaycountToARMDaycount((const char*) ff1);

				theNode->Release();
				theNode=NULL;
				if (resultat) SysFreeString(resultat);
			}

			// reference Date
			listItem->selectSingleNode(_bstr_t((const char *)"STUB_Date1"), &theNode);
			if (theNode!=NULL)
			{
				BSTR resultat = NULL;
				theNode->get_text(&resultat);

				_bstr_t ff(resultat,false);
				char * ff1=(char *)ff;

				if (strcmp((const char*)ff1,"") == 0)
					strcpy(ReferenceDate,"NULL");
				else
				{
					ARM_Date tmpDate (ff1,"YYYYMMDD");
					startStubDate = tmpDate;
				}

				theNode->Release();
				theNode=NULL;
				if (resultat) SysFreeString(resultat);
			}

			// reference Date
			listItem->selectSingleNode(_bstr_t((const char *)"STUB_Date2"), &theNode);
			if (theNode!=NULL)
			{
				BSTR resultat = NULL;
				theNode->get_text(&resultat);

				_bstr_t ff(resultat,false);
				char * ff1=(char *)ff;

				if (strcmp((const char*)ff1,"") != 0)
				{
				ARM_Date tmpDate (ff1,"YYYYMMDD");
				endStubDate = tmpDate;
				}

				theNode->Release();
				theNode=NULL;
				if (resultat) SysFreeString(resultat);
			}

			// Fst Cpn Effective Date			
			listItem->selectSingleNode(_bstr_t((const char *)"FstCpEffect"), &theNode);
			if (theNode!=NULL)
			{
				BSTR resultat = NULL;
				theNode->get_text(&resultat);

				_bstr_t ff(resultat,false);
				char * ff1=(char *)ff;

				if (strcmp((const char*)ff1,"") != 0)
				{
					isFstCpEffect=true ; 
					ARM_Date tmpDate (ff1,"YYYYMMDD");
					fstCpEffectDate= tmpDate;
				}

				theNode->Release();
				theNode=NULL;
				if (resultat) SysFreeString(resultat);
			}


			// reference Day1
			listItem->selectSingleNode(_bstr_t((const char *)"SCHED_Pay_AnnDay"), &theNode);
			if (theNode!=NULL)
			{
				BSTR resultat = NULL;
				theNode->get_text(&resultat);

				_bstr_t ff(resultat,false);
				char * ff1=(char *)ff;

				refDay = atoi((const char*) ff1);

				theNode->Release();
				theNode=NULL;
				if (resultat) SysFreeString(resultat);
			}


			//Recuperation du stub
			listItem->selectSingleNode(_bstr_t("STUB_StubType"), &theNode);
			if (theNode!=NULL)
			{
				BSTR resultat = NULL;
				theNode->get_text(&resultat);

				_bstr_t ff(resultat,false);
				char * ff1=(char *)ff;
				stub = SummitStub2ARMStub(ff1);

				theNode->Release();
				theNode=NULL;
				if (resultat) SysFreeString(resultat);
			}

			listItem->selectSingleNode(_bstr_t("Ccy"), &theNode);
			if (theNode!=NULL)
			{
				BSTR resultat = NULL;
				theNode->get_text(&resultat);

				_bstr_t ff(resultat,false);
				ccy=ff; 
				// char * ff1=(char *)ff;
				// ccy = new ARM_Currency(ff1);

				theNode->Release();
				theNode=NULL;
				if (resultat) SysFreeString(resultat);
			}

			if (startStubDate<=StartDate || startStubDate>=EndDate) startStubDate=ARM_Date(); 
			if (endStubDate<=StartDate || endStubDate>=EndDate) endStubDate=ARM_Date(); 
			DeduceRefDateAndStub(StartDate,EndDate,startStubDate,endStubDate,payFreq,refDay,ccy,(char*)ReferenceDate,stub);

			FixedPayerAmount= XML_doubleNodeTreating(listItem,"Notional");
			FloatingPayerAmount=FixedPayerAmount;

			//Recuperation du PorS
			listItem->selectSingleNode(_bstr_t("PorS"), &theNode);
			if (theNode!=NULL)
			{
				BSTR resultat = NULL;
				theNode->get_text(&resultat);

				_bstr_t ff(resultat,false);
				char * ff1=(char *)ff;
				
				if (!strcmp(ff1,"P"))
					SensDeLaTransaction = K_RCV;
				else
					SensDeLaTransaction = K_PAY;

				theNode->Release();
				theNode=NULL;
				if (resultat) SysFreeString(resultat);
			}

			{
				// Recupération de l'INTEREST_Method (Last Day Interest) 
				std::string method ;
				if ( GetStringFromXMLNode(listItem, "INTEREST_Method",method) )
				{
					if (method =="") excludeMaturity=EXCLUDE_MATURITY ; 
					else if (method =="STAND") excludeMaturity=EXCLUDE_MATURITY ; 
					else if (method =="INCL") excludeMaturity=INCLUDE_MATURITY; 
					// else if (method =="EFF") excludeMaturity=INCLUDE_MATURITY ; 
					else 
					{
						ICMTHROW(ERR_INVALID_ARGUMENT,"Can't understand INTEREST_Method. "<< method<<" for " << tradeId ); 
					}
				}
			}

			listItem->Release();
			listItem=NULL;

			if (startStubDate<=StartDate || startStubDate>=EndDate) startStubDate=ARM_Date(); 
			if (endStubDate<=StartDate || endStubDate>=EndDate) endStubDate=ARM_Date(); 
			DeduceRefDateAndStub(StartDate,EndDate,startStubDate,endStubDate,payFreq,refDay,ccy,(char*)ReferenceDate,stub);



		}

		if (resultList)
		{
			resultList->Release();
			resultList = NULL;
		}

		// ---------------------------------------------------------------------------------
		// Recuperation du spread fix
		// ---------------------------------------------------------------------------------

		if (XMLDoc->selectNodes((_bstr_t)("//ASSET/ProdData/CREDSWAP"), &resultList) == S_OK)
		{
			resultList->get_length(&nbNodes);

			if (nbNodes == 0)
			{
				hr = S_FALSE;
				throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT, "Invalid XML string for getting ");
			}

			hr=resultList->get_item(0, &listItem);
			spread = XML_doubleNodeTreating(listItem,"Formula");
			spread /= 100.;

			listItem->Release();
			listItem=NULL;

			hr=resultList->get_item(0, &listItem);
			listItem->selectSingleNode(_bstr_t("SettleMethod"), &theNode);
			if (theNode!=NULL)
			{
				BSTR resultat = NULL;
				theNode->get_text(&resultat);

				_bstr_t ff(resultat,false);
				char * ff1=(char *)ff;
				
				if (strcmp(ff1,"BIN") != NULL) Binary = CREDIT_DEFAULT_VALUE;

				theNode->Release();
				theNode=NULL;
				if (resultat) SysFreeString(resultat);
			}
			listItem->Release();
			listItem=NULL;

			if (Binary != CREDIT_DEFAULT_VALUE)
			{
				hr=resultList->get_item(0, &listItem);
				Binary = 1. - XML_doubleNodeTreating(listItem,"BinaryAmt") ;
				listItem->Release();
				listItem=NULL;
			}

			hr=resultList->get_item(0, &listItem);
			listItem->selectSingleNode(_bstr_t("DelivDays"), &theNode);
			if (theNode!=NULL)
			{
				BSTR resultat = NULL;
				theNode->get_text(&resultat);

				_bstr_t ff(resultat,false);
				char * ff1=(char *)ff;
				
				sscanf(ff1, "%2dD", &DelivDays);

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

		// ------------- Recuperation du Defprob issuer (if any) -----------------------
		// 
		XMLDoc->selectSingleNode((_bstr_t)("//ASSET/ProdData/CREDSWAP/CreditEntList/CREDIT/DefpIssuer"), &theNode)  ;
		std::string defprobIssuer ;
		if (theNode!=NULL)
		{
			BSTR resultat = NULL;
			theNode->get_text(&resultat);

			_bstr_t ff(resultat,false);
			defprobIssuer = (char*)ff; 

			theNode->Release();
			theNode=NULL;
			if (resultat) SysFreeString(resultat);
		}

		//------------------------------------------------------------------------------
		//Création de l'instrument financier
		//------------------------------------------------------------------------------

		Local_XML_TRACE("c:\\temp\\log_arm_credit.txt","Création du NTD");

		if (EndDate <= StartDate)
		{
			Local_XML_TRACE("c:\\temp\\log_arm_credit.txt","Error: EndDate <= StartDate");
			throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT, "Error: EndDate <= StartDate");
		}


// 		if (excludeMaturity!=INCLUDE_MATURITY)
// 			ICMTHROW(ERR_INVALID_ARGUMENT,"Wrong handling of INCLUDE MATURITY for " << tradeId) ;

		(*notional) *= SensDeLaTransaction; 
		bool toto (isFstCpEffect) ; 
		if (isFstCpEffect && fstCpEffectDate>=StartDate) toto=false;
		ICM_Cds* cds = new ICM_Cds(StartDate,
									 EndDate,
									 strcmp(ReferenceDate,"NULL")==0 ? 0 : &ARM_Date(ReferenceDate),
									 toto?&fstCpEffectDate:0,
									 StartDate,
									 EndDate,
									 spread,
									 // 0,0,
									 *notional.get(),
									 *notional.get(),
									 payFreq,
									 daycount,
									 qACCRUED_SETTLED,
									 ccy,
									 stub,
									 DelivDays,
									 DEFAULT_FRQ_DEFLEG,
									 interestRule,
									 excludeMaturity,
									 adjStartDateRule,
									 PayCal,
									qRunning_Leg,
									qStandart_Recovery_Leg,
									 ISSUER_UNDEFINE,
									 Binary);

		if (defprobIssuer!="") cds->SetSingleName(defprobIssuer); 
		//ARM_ReferenceValue* XNL = ARMLOCAL_XML_EVENT(chaineXML,"XNL");
		//if (XNL) ntd->GetFeeLeg()->SetAmount(XNL);
		//if (XNL) delete (XNL);

		//if (ccy) delete ccy;

		Local_XML_TRACE("c:\\temp\\log_arm_credit.txt","Fin Création du NTD");

		if (Rate) delete[] Rate;		
		//if (PayCal) delete[] PayCal;

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

		if (XMLDoc)	{ XMLDoc->Release();}
		XMLDoc = NULL;

		if (xmlWCharText)
			delete[] xmlWCharText;
		xmlWCharText = NULL;

		return (cds);
	}
	catch(Exception& )
	{
		if (theNode) theNode->Release();
		if (listItem) listItem->Release();
		if (resultList) resultList->Release();
		if (XMLDoc) XMLDoc->Release();


		if (xmlWCharText)
			delete[] xmlWCharText;
		xmlWCharText = NULL;

		hr = S_FALSE;
		throw; 
	}
	catch(...)
	{		
		if (theNode) theNode->Release();
		if (listItem) listItem->Release();
		if (resultList) resultList->Release();
		if (XMLDoc) XMLDoc->Release();

		if (xmlWCharText)
			delete[] xmlWCharText;
		xmlWCharText = NULL;

		hr = S_FALSE;
		throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
         "Error in XML parsing for getting Tranche");
	}

	return NULL;
}


// ------------------------------------
// Generic Credit Index
// ------------------------------------

/** 
ICM_Credit_Index* ARMLOCAL_XML_CREDIT_INDEX(const char* chaineXML, 
											CCString& bookId, 
											CCString& custId, 
											CCString& dealId)
{
	VARIANT_BOOL bOK;
	HRESULT hr;

	MSXML2::IXMLDOMDocument *XMLDoc = NULL;
	MSXML2::IXMLDOMNodeList * resultList = NULL;
	MSXML2::IXMLDOMNodeList * resultList0 = NULL;
	MSXML2::IXMLDOMNode * listItem = NULL, * theNode = NULL;
	long nbNodes = 0;

	int l_daycount;
	int l_SCHED_Reset_Freq;
	int l_SCHED_Pay_Freq;
	double CST_yearterm = 5.;
	int CST_NbCurves =0;
	std::string ccy("");
	double CST_Spread = 0.;
	int l_SCHED_Reset_Rule;
	int l_SCHED_Pay_IntRule;
	int l_SCHED_Reset_Time;
	int l_SCHED_Reset_Gap;
	int l_SCHED_Pay_Time;
	int l_SCHED_Pay_Gap;

	string str_SCHED_Reset_Freq;

	int k =0;	

	wchar_t * xmlWCharText = NULL;

	try
	{
		hr = CoInitialize(NULL); 
		SUCCEEDED(hr) ? 0 : throw hr;
		hr = CoCreateInstance(__uuidof(MSXML2::DOMDocument30), NULL, CLSCTX_INPROC_SERVER, __uuidof(MSXML2::IXMLDOMDocument), (void**)&XMLDoc);
		SUCCEEDED(hr) ? 0 : throw hr;

		xmlWCharText = constchar2wchar(chaineXML);

		Local_XML_TRACE("c:\\temp\\log_arm_credit.txt","Lecture XML pour Matrice de Correlation");
		XMLDoc->loadXML((_bstr_t)xmlWCharText, &bOK);

	}
	catch(Exception& )
	{
		if (XMLDoc) XMLDoc->Release();
		hr = S_FALSE;
		throw; 


		if (xmlWCharText)
			delete[] xmlWCharText;
		xmlWCharText = NULL;
	}
	catch(...)
	{
		if (XMLDoc) XMLDoc->Release();
		hr = S_FALSE;

		if (xmlWCharText)
			delete[] xmlWCharText;
		xmlWCharText = NULL;

		throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"Pb in creating XML document for getting Tranches informations");
	}

	try
	{
		{
			// JLA. Retrieve 
			MSXML2::IXMLDOMNode* aNode = NULL;
			std::string chemin = "Response/EXOTIC/Env/ENV";
			XMLDoc->selectSingleNode((_bstr_t)chemin.c_str(), &aNode);
			if (aNode)
			{
				custId= GetStringFromXMLNode(aNode, "Cust").c_str();
				dealId = GetStringFromXMLNode(aNode, "DealId").c_str();
				bookId = GetStringFromXMLNode(aNode, "Book").c_str();
				aNode->Release();
			}
		}


		// ---------------------------------------------------------------------------------
		// Recuperation des infos communes à tous les deals
		// ---------------------------------------------------------------------------------

		Local_XML_TRACE("c:\\temp\\log_arm_credit.txt","Recuperation des infos communes à tous les deals");

 		if (XMLDoc->selectNodes((_bstr_t)("//ASSET"), &resultList) == S_OK)
		{

			hr=resultList->get_item(0, &listItem);

			//Recuperation du basis daycount
			listItem->selectSingleNode(_bstr_t("INTEREST_Basis"), &theNode);
			if (theNode!=NULL)
			{
				BSTR resultat = NULL;
				theNode->get_text(&resultat);

				_bstr_t ff(resultat,false);
				char * ff1=(char *)ff;
				l_daycount = FromSummitDaycountToARMDaycount((const char*) ff1);

				theNode->Release();
				theNode=NULL;
				if (resultat) SysFreeString(resultat);
			}

			//Recuperation de frequence
			listItem->selectSingleNode(_bstr_t("SCHED_Reset_Freq"), &theNode);
			if (theNode!=NULL)
			{
				BSTR resultat = NULL;
				theNode->get_text(&resultat);

				_bstr_t ff(resultat,false);
				char * ff1=(char *)ff;
				l_SCHED_Reset_Freq = FromSummitFreqToARMFreq((const char*) ff1);

				theNode->Release();
				theNode=NULL;
				if (resultat) SysFreeString(resultat);
			}

			//Recuperation de frequence
			listItem->selectSingleNode(_bstr_t("SCHED_Pay_Freq"), &theNode);
			if (theNode!=NULL)
			{
				BSTR resultat = NULL;
				theNode->get_text(&resultat);

				_bstr_t ff(resultat,false);
				char * ff1=(char *)ff;
				l_SCHED_Pay_Freq = FromSummitFreqToARMFreq((const char*) ff1);

				theNode->Release();
				theNode=NULL;
				if (resultat) SysFreeString(resultat);
			}

			listItem->selectSingleNode(_bstr_t("Ccy"), &theNode);
			if (theNode!=NULL)
			{
				BSTR resultat = NULL;
				theNode->get_text(&resultat);

				_bstr_t ff(resultat,false);
				char * ff1=(char *)ff;
				//ccy = new ARM_Currency(ff1);
				ccy = string(ff1);
				theNode->Release();
				theNode=NULL;
				if (resultat) SysFreeString(resultat);
			}

			//Recuperation de frequence
			listItem->selectSingleNode(_bstr_t("SCHED_Reset_Rule"), &theNode);
			if (theNode!=NULL)
			{
				BSTR resultat = NULL;
				theNode->get_text(&resultat);

				_bstr_t ff(resultat,false);
				char * ff1=(char *)ff;
				l_SCHED_Reset_Rule = ARM_ConvFwdRule((const char*) ff1);

				theNode->Release();
				theNode=NULL;
				if (resultat) SysFreeString(resultat);
			}


			//Recuperation de frequence
			listItem->selectSingleNode(_bstr_t("SCHED_Pay_IntRule"), &theNode);
			if (theNode!=NULL)
			{
				BSTR resultat = NULL;
				theNode->get_text(&resultat);

				_bstr_t ff(resultat,false);
				char * ff1=(char *)ff;
				l_SCHED_Pay_IntRule = ARM_ConvIntRule((const char*) ff1);

				theNode->Release();
				theNode=NULL;
				if (resultat) SysFreeString(resultat);
			}

			//Recuperation de frequence
			listItem->selectSingleNode(_bstr_t("SCHED_Reset_Time"), &theNode);
			if (theNode!=NULL)
			{
				BSTR resultat = NULL;
				theNode->get_text(&resultat);

				_bstr_t ff(resultat,false);
				char * ff1=(char *)ff;
				l_SCHED_Reset_Time = ARM_ConvPayResetRule((const char*) ff1);

				theNode->Release();
				theNode=NULL;
				if (resultat) SysFreeString(resultat);
			}

			//Recuperation de frequence
			listItem->selectSingleNode(_bstr_t("SCHED_Reset_Gap"), &theNode);
			if (theNode!=NULL)
			{
				BSTR resultat = NULL;
				theNode->get_text(&resultat);

				_bstr_t ff(resultat,false);
				char * ff1=(char *)ff;
				l_SCHED_Reset_Gap = FromSummitGapToARMGap((const char*) ff1);

				theNode->Release();
				theNode=NULL;
				if (resultat) SysFreeString(resultat);
			}

			//Recuperation de frequence
			listItem->selectSingleNode(_bstr_t("SCHED_Pay_Time"), &theNode);
			if (theNode!=NULL)
			{
				BSTR resultat = NULL;
				theNode->get_text(&resultat);

				_bstr_t ff(resultat,false);
				char * ff1=(char *)ff;
				l_SCHED_Pay_Time = ARM_ConvPayResetRule((const char*) ff1);

				theNode->Release();
				theNode=NULL;
				if (resultat) SysFreeString(resultat);
			}

			//Recuperation de frequence
			listItem->selectSingleNode(_bstr_t("SCHED_Pay_Gap"), &theNode);
			if (theNode!=NULL)
			{
				BSTR resultat = NULL;
				theNode->get_text(&resultat);

				_bstr_t ff(resultat,false);
				char * ff1=(char *)ff;
				l_SCHED_Pay_Gap = FromSummitGapToARMGap((const char*) ff1);

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


       ICM_Credit_Index* index = new ICM_Credit_Index(l_daycount,
												  l_SCHED_Reset_Freq, 
												  l_SCHED_Pay_Freq, 
												  CST_yearterm,
											      (ARM_Date*)0,
												  ccy,
												  CST_NbCurves,
												  NULL,
												  qAVERAGE,
												  CST_Spread,
												  // NULL, maturity
												  NULL, //ForcedCurve
												  NULL, //Volatility
												  l_SCHED_Reset_Rule,
												  l_SCHED_Pay_IntRule,	
												  l_SCHED_Reset_Time,	
												  l_SCHED_Reset_Gap,	
											      l_SCHED_Pay_Time,
												  l_SCHED_Pay_Gap,
											      qCredit_Adjust20,5,2);

		Local_XML_TRACE("c:\\temp\\log_arm_credit.txt","Fin Création du NTD");

		//if (ccy) delete ccy;		

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

		if (XMLDoc)	{ XMLDoc->Release();}
		XMLDoc = NULL;

		if (xmlWCharText)
			delete[] xmlWCharText;
		xmlWCharText = NULL;

		return (index);
	}
	catch(Exception& )
	{
		if (theNode) theNode->Release();
		if (listItem) listItem->Release();
		if (resultList) resultList->Release();
		if (XMLDoc) XMLDoc->Release();


		if (xmlWCharText)
			delete[] xmlWCharText;
		xmlWCharText = NULL;

		hr = S_FALSE;
		throw; 
	}
	catch(...)
	{		
		if (theNode) theNode->Release();
		if (listItem) listItem->Release();
		if (resultList) resultList->Release();
		if (XMLDoc) XMLDoc->Release();

		if (xmlWCharText)
			delete[] xmlWCharText;
		xmlWCharText = NULL;

		hr = S_FALSE;
		throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
         "Error in XML parsing for getting Tranche");
	}

	return NULL;
}
**/ 


/** JLA useless 

  // ------------------------------------
// Generic Tranche Parsing for Cdo 2
// ------------------------------------

ICM_Cdo2* ARMLOCAL_XML_CDO2_EMPTY(const char* chaineXML, 
								CCString& bookId, 
								CCString& custId, 
								CCString& dealId,
								ICM_Portfolio* Pf)
{
	VARIANT_BOOL bOK;
	HRESULT hr;

	MSXML2::IXMLDOMDocument *XMLDoc = NULL;
	MSXML2::IXMLDOMNodeList * resultList = NULL;
	MSXML2::IXMLDOMNodeList * resultList0 = NULL;
	MSXML2::IXMLDOMNode * listItem = NULL, * theNode = NULL;
	long nbNodes = 0;

	CCString strtmp;

	ARM_Date StartDate;
	ARM_Date EndDate;
	ARM_Date fstCpEffectDate ;
	bool isFstCpEffect(false); 
	double spread = 0.;

	char ReferenceDate[11];
	int payFreq = K_DEF_FREQ;
	int daycount = KACTUAL_360;
	int ACC_daycount = KACTUAL_360;
	int DelivDays=0 ;

	double FixedPayerAmount=0.;
	double FloatingPayerAmount=0.;
	// ARM_Currency* ccy = NULL;
	std::string ccy; 

	double MezzAmount = 0.;
	double SubAmount = 0.;

	int nbissuers = 0;
	char** Issuers = NULL; 
	double* Notionals = NULL;
	double* Rate = NULL;
	double CoefAdjust = 0.;
	// char* PayCal = NULL;
	std::string PayCal ;


	int interestRule,adjStartDateRule; 
	double SensDeLaTransaction = 1.;
	double Binary = 0.;

	//ARM_Currency* discountCcy = NULL;
	ARM_ReferenceValue* notional = NULL;
	ICM_ProportionsInfo* info = NULL;
	char SmilProp[500];
	memset(SmilProp,'\0',sizeof(SmilProp));

	int k =0;	
	int IntegrationMethod = CREDIT_DEFAULT_VALUE;
	
	wchar_t * xmlWCharText = NULL;

	ARM_Date startStubDate;
	ARM_Date endStubDate;
	int refDay = 0;
	int stub = K_SHORTSTART;

	std::string tradeId; 
	bool excludeMaturity=EXCLUDE_MATURITY; 
	try
	{
		hr = CoInitialize(NULL); 
//		SUCCEEDED(hr) ? 0 : throw hr;
		hr = CoCreateInstance(__uuidof(MSXML2::DOMDocument30), NULL, CLSCTX_INPROC_SERVER, __uuidof(MSXML2::IXMLDOMDocument), (void**)&XMLDoc);
		SUCCEEDED(hr) ? 0 : throw hr;

		xmlWCharText = constchar2wchar(chaineXML);

		Local_XML_TRACE("c:\\temp\\log_arm_credit.txt","Lecture XML pour Matrice de Correlation");
		XMLDoc->loadXML((_bstr_t)xmlWCharText, &bOK);

		
		#ifdef _DEBUG
		CCString chtmp = "c:\\temp\\AAA_CDO_OUTPUT.xml";
		VARIANT v_chtmp;
		CCString2VARIANT (chtmp, &v_chtmp);
		XMLDoc->save(v_chtmp);
		#endif


	}
	catch(Exception& x)
	{
		if (XMLDoc) XMLDoc->Release();
		hr = S_FALSE;
		char Tmp[1000];
		x.GetErrorMessage(Tmp);
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,Tmp);


		if (xmlWCharText)
			delete[] xmlWCharText;
		xmlWCharText = NULL;
	}
	catch(...)
	{
		if (XMLDoc) XMLDoc->Release();
		hr = S_FALSE;

		if (xmlWCharText)
			delete[] xmlWCharText;
		xmlWCharText = NULL;

		throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"Pb in creating XML document for getting Tranches informations");
	}

	try
	{
		{
			// JLA. Retrieve 
			MSXML2::IXMLDOMNode* aNode = NULL;
			std::string chemin = "Response/EXOTIC/Env/ENV";
			XMLDoc->selectSingleNode((_bstr_t)chemin.c_str(), &aNode);
			if (aNode)
			{
			custId= GetStringFromXMLNode(aNode, "Cust").c_str();
			dealId = GetStringFromXMLNode(aNode, "DealId").c_str();
			bookId = GetStringFromXMLNode(aNode, "Book").c_str();
			tradeId= GetStringFromXMLNode(aNode, "TradeId");
			aNode->Release();
			}
		}

		// ------------------------------------------------------------------------
		// Pricing avec correl par strike ?
		// ------------------------------------------------------------------------
		strtmp = "//ASSET/ProdData/CREDSWAP";

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

			hr=resultList->get_item(0, &listItem);
			
			//Recuperation du type de la courbe
			listItem->selectSingleNode((_bstr_t)"cSmileProp", &theNode);
			if (theNode!=NULL)
			{
				BSTR resultat = NULL;
				theNode->get_text(&resultat);

				_bstr_t ff(resultat,false);
				char * ff1=(char *)ff;
				strcpy(SmilProp,ff1);

				theNode->Release();
				theNode=NULL;
				if (resultat) SysFreeString(resultat);
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

			if (strcmp(SmilProp,"") == NULL)
			{
				if (XMLDoc) XMLDoc->Release();
				XMLDoc = NULL;

				if (xmlWCharText)
					delete[] xmlWCharText;
				xmlWCharText = NULL;

				return NULL;
			}

			CvtStrPropInfo(SmilProp,info);

		}

		// ---------------------------------------------------------------------------------
		// Recuperation des infos communes à tous les deals
		// ---------------------------------------------------------------------------------

		Local_XML_TRACE("c:\\temp\\log_arm_credit.txt","Recuperation des infos communes à tous les deals");

 		if (XMLDoc->selectNodes((_bstr_t)("//OPSPEC"), &resultList) == S_OK)
		{

			hr=resultList->get_item(0, &listItem);
			
			if (listItem)
			{
			//Recuperation de la Trade Date
			IntegrationMethod = (int) XML_doubleNodeTreating(listItem,"NumIter");

			if (IntegrationMethod == 0) IntegrationMethod = INTEGRATION_STEP;

			listItem->Release();
			listItem=NULL;
			}
		}

		if (resultList)
		{
			resultList->Release();
			resultList = NULL;
		}

 		if (XMLDoc->selectNodes((_bstr_t)("//ASSET"), &resultList) == S_OK)
		{

			hr=resultList->get_item(0, &listItem);
			
			//Recuperation de la Trade Date
			StartDate = GetStartDate(listItem);

			//Recuperation de la Maturity Date
			EndDate = GetEndDate(listItem);

			//Recuperation de frequence
			listItem->selectSingleNode(_bstr_t("SCHED_Pay_Freq"), &theNode);
			if (theNode!=NULL)
			{
				BSTR resultat = NULL;
				theNode->get_text(&resultat);

				_bstr_t ff(resultat,false);
				char * ff1=(char *)ff;
				payFreq = FromSummitFreqToARMFreq((const char*) ff1);

				theNode->Release();
				theNode=NULL;
				if (resultat) SysFreeString(resultat);
			}

			//Recuperation de frequence
			listItem->selectSingleNode(_bstr_t("SCHED_Pay_Cal"), &theNode);
			if (theNode!=NULL)
			{
				BSTR resultat = NULL;
				theNode->get_text(&resultat);

				_bstr_t ff(resultat,false);
				PayCal =  ff;
				// if (PayCal)
					// delete[] PayCal;
				// PayCal = new char[5];
				// strcpy(PayCal,ff1);

				theNode->Release();
				theNode=NULL;
				if (resultat) SysFreeString(resultat);
			}

			//Recuperation du basis daycount
			listItem->selectSingleNode(_bstr_t("INTEREST_Basis"), &theNode);
			if (theNode!=NULL)
			{
				BSTR resultat = NULL;
				theNode->get_text(&resultat);

				_bstr_t ff(resultat,false);
				char * ff1=(char *)ff;
				daycount = FromSummitDaycountToARMDaycount((const char*) ff1);

				theNode->Release();
				theNode=NULL;
				if (resultat) SysFreeString(resultat);
			}

			//Recuperation du Accrual basis daycount
			listItem->selectSingleNode(_bstr_t("INT_ACC_Basis"), &theNode);
			if (theNode!=NULL)
			{
				BSTR resultat = NULL;
				theNode->get_text(&resultat);

				_bstr_t ff(resultat,false);
				char * ff1=(char *)ff;
				ACC_daycount = FromSummitDaycountToARMDaycount((const char*) ff1);

				theNode->Release();
				theNode=NULL;
				if (resultat) SysFreeString(resultat);
			}

			// reference Date
			listItem->selectSingleNode(_bstr_t((const char *)"STUB_Date1"), &theNode);
			if (theNode!=NULL)
			{
				BSTR resultat = NULL;
				theNode->get_text(&resultat);

				_bstr_t ff(resultat,false);
				char * ff1=(char *)ff;

				if (strcmp((const char*)ff1,"") == 0)
					strcpy(ReferenceDate,"NULL");
				else
				{
					ARM_Date tmpDate (ff1,"YYYYMMDD");
					startStubDate = tmpDate;
				}

				theNode->Release();
				theNode=NULL;
				if (resultat) SysFreeString(resultat);
			}

			// reference Date
			listItem->selectSingleNode(_bstr_t((const char *)"STUB_Date2"), &theNode);
			if (theNode!=NULL)
			{
				BSTR resultat = NULL;
				theNode->get_text(&resultat);

				_bstr_t ff(resultat,false);
				char * ff1=(char *)ff;

				if (strcmp((const char*)ff1,"") != 0)
				{
				ARM_Date tmpDate (ff1,"YYYYMMDD");
				endStubDate = tmpDate;
				}

				theNode->Release();
				theNode=NULL;
				if (resultat) SysFreeString(resultat);
			}
			// Fst Cpn Effective Date
			
			listItem->selectSingleNode(_bstr_t((const char *)"FstCpEffect"), &theNode);
			if (theNode!=NULL)
			{
				BSTR resultat = NULL;
				theNode->get_text(&resultat);

				_bstr_t ff(resultat,false);
				char * ff1=(char *)ff;

				if (strcmp((const char*)ff1,"") != 0)
				{
					isFstCpEffect=true ; 
					ARM_Date tmpDate (ff1,"YYYYMMDD");
					fstCpEffectDate= tmpDate;
				}

				theNode->Release();
				theNode=NULL;
				if (resultat) SysFreeString(resultat);
			}
			// reference Day1
			listItem->selectSingleNode(_bstr_t((const char *)"SCHED_Pay_AnnDay"), &theNode);
			if (theNode!=NULL)
			{
				BSTR resultat = NULL;
				theNode->get_text(&resultat);

				_bstr_t ff(resultat,false);
				char * ff1=(char *)ff;

				refDay = atoi((const char*) ff1);

				theNode->Release();
				theNode=NULL;
				if (resultat) SysFreeString(resultat);
			}

			//Recuperation du stub
			listItem->selectSingleNode(_bstr_t("STUB_StubType"), &theNode);
			if (theNode!=NULL)
			{
				BSTR resultat = NULL;
				theNode->get_text(&resultat);

				_bstr_t ff(resultat,false);
				char * ff1=(char *)ff;
				stub = SummitStub2ARMStub(ff1);

				theNode->Release();
				theNode=NULL;
				if (resultat) SysFreeString(resultat);
			}

			listItem->selectSingleNode(_bstr_t("Ccy"), &theNode);
			if (theNode!=NULL)
			{
				BSTR resultat = NULL;
				theNode->get_text(&resultat);

				_bstr_t ff(resultat,false);
				ccy=ff; 

				//char * ff1=(char *)ff;
				//ccy = new ARM_Currency(ff1);

				theNode->Release();
				theNode=NULL;
				if (resultat) SysFreeString(resultat);
			}

			// interestRule
			interestRule = GetIntRule(listItem) ;
			if (interestRule==K_UNADJUSTED) adjStartDateRule=K_UNADJUSTED; 
			else adjStartDateRule=K_ADJUSTED ; 

			if (startStubDate<=StartDate || startStubDate>=EndDate) startStubDate=ARM_Date(); 
			if (endStubDate<=StartDate || endStubDate>=EndDate) endStubDate=ARM_Date(); 
			DeduceRefDateAndStub(StartDate,EndDate,startStubDate,endStubDate,payFreq,refDay,ccy,(char*)ReferenceDate,stub);

			FixedPayerAmount= XML_doubleNodeTreating(listItem,"Notional");
			FloatingPayerAmount=FixedPayerAmount;

			CoefAdjust = FixedPayerAmount;

			//Recuperation du PorS
			listItem->selectSingleNode(_bstr_t("PorS"), &theNode);
			if (theNode!=NULL)
			{
				BSTR resultat = NULL;
				theNode->get_text(&resultat);

				_bstr_t ff(resultat,false);
				char * ff1=(char *)ff;
				
				if (!strcmp(ff1,"P"))
					SensDeLaTransaction = 1.;
				else
					SensDeLaTransaction = -1.;

				theNode->Release();
				theNode=NULL;
				if (resultat) SysFreeString(resultat);
			}
			//	Recuperation de INTEREST_Method
			{
				std::string method ;
				if ( GetStringFromXMLNode(listItem, "INTEREST_Method",method) ) 
				{
					if (method =="") excludeMaturity=EXCLUDE_MATURITY ; 
					else if (method =="STAND") excludeMaturity=EXCLUDE_MATURITY ; 
					else if (method =="INCL") excludeMaturity=INCLUDE_MATURITY; 
					// else if (method =="EFF") excludeMaturity=INCLUDE_MATURITY ; 
					else 
					{
						ICMTHROW(ERR_INVALID_ARGUMENT,"Can't understand INTEREST_Method. "<< method<<" for " << tradeId ); 
					}
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

		// ---------------------------------------------------------------------------------
		// Recuperation du spread fix
		// ---------------------------------------------------------------------------------

		if (XMLDoc->selectNodes((_bstr_t)("//ASSET/ProdData/CREDSWAP"), &resultList) == S_OK)
		{
			resultList->get_length(&nbNodes);

			if (nbNodes == 0)
			{
				hr = S_FALSE;
				throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT, "Invalid XML string for getting ");
			}

			hr=resultList->get_item(0, &listItem);
			spread = XML_doubleNodeTreating(listItem,"Formula");
			spread /= 100.;

			listItem->Release();
			listItem=NULL;

			hr=resultList->get_item(0, &listItem);
			listItem->selectSingleNode(_bstr_t("SettleMethod"), &theNode);
			if (theNode!=NULL)
			{
				BSTR resultat = NULL;
				theNode->get_text(&resultat);

				_bstr_t ff(resultat,false);
				char * ff1=(char *)ff;
				
				if (strcmp(ff1,"BIN") != NULL) Binary = CREDIT_DEFAULT_VALUE;

				theNode->Release();
				theNode=NULL;
				if (resultat) SysFreeString(resultat);
			}
			listItem->Release();
			listItem=NULL;

			if (Binary != CREDIT_DEFAULT_VALUE)
			{
				hr=resultList->get_item(0, &listItem);
				Binary = 1. - XML_doubleNodeTreating(listItem,"BinaryAmt");
				listItem->Release();
				listItem=NULL;
			}

			hr=resultList->get_item(0, &listItem);
			listItem->selectSingleNode(_bstr_t("DelivDays"), &theNode);
			if (theNode!=NULL)
			{
				BSTR resultat = NULL;
				theNode->get_text(&resultat);

				_bstr_t ff(resultat,false);
				char * ff1=(char *)ff;
				
				sscanf(ff1, "%2dD", &DelivDays);

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

		Local_XML_TRACE("c:\\temp\\log_arm_credit.txt","Recuperation des infos sur les Issuers");

			if (XMLDoc->selectNodes((_bstr_t)("//ASSET/ProdData/CREDSWAP/CreditEntList/CREDIT[DefpIssuer != \"FTD1\"]"), &resultList0) == S_OK)
			{
				resultList0->get_length(&nbNodes);

				if (nbNodes == 0)
				{	
					hr = S_FALSE;
					throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT, "Invalid XML string for getting ");
				}

				Issuers = new char*[nbNodes];
				nbissuers = nbNodes; 

				for (long indexNode=0 ; indexNode<nbNodes ; indexNode++)
				{
					hr=resultList0->get_item(indexNode, &listItem);

					listItem->selectSingleNode(_bstr_t("DefpIssuer"), &theNode);
					if (theNode!=NULL)
					{
					BSTR resultat = NULL;
					theNode->get_text(&resultat);

					_bstr_t ff(resultat,false);
					char * ff1=(char *)ff;

					Issuers[indexNode] = new char[1000];
					strcpy(Issuers[indexNode],ff1);
		
					theNode->Release();
					theNode=NULL;
					if (resultat) SysFreeString(resultat);
					}

					listItem->Release();
					listItem = NULL;

				}

				if (resultList0)
				{
					resultList0->Release();
					resultList0 = NULL;
				}


			if (XMLDoc->selectNodes((_bstr_t)("//ASSET/ProdData/CREDSWAP/CreditEntList/CREDIT"), &resultList) == S_OK)
			{
				resultList->get_length(&nbNodes);

				if (nbNodes == 0)
				{	
					hr = S_FALSE;
					throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT, "Invalid XML string for getting ");
				}

				Notionals = new double[nbNodes];

				for (long indexNode=0 ; indexNode<nbNodes ; indexNode++)
				{
					hr=resultList->get_item(indexNode, &listItem);

					Notionals[indexNode] = XML_doubleNodeTreating(listItem,"Quantity");

					listItem->Release();
					listItem = NULL;

				}

				if (theNode)
				{
				theNode->Release();
				theNode=NULL;
				}
	
			}

			if (resultList)
			{
			resultList->Release();
			resultList = NULL;
			}

			if (XMLDoc->selectNodes((_bstr_t)("//ASSET/ProdData/CREDSWAP"), &resultList) == S_OK)
			{
				resultList->get_length(&nbNodes);

				if (nbNodes == 0)
				{
				hr = S_FALSE;
				throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT, "Invalid XML string for getting ");
				}

				hr=resultList->get_item(0, &listItem);

				listItem->selectSingleNode(_bstr_t("Tranched"), &theNode);
				if (theNode!=NULL)
				{
				BSTR resultat = NULL;
				theNode->get_text(&resultat);

				_bstr_t ff(resultat,false);
				char * ff1=(char *)ff;

				theNode->Release();
				theNode=NULL;
				if (resultat) SysFreeString(resultat);
				}

				MezzAmount = XML_doubleNodeTreating(listItem,"ProtectStartAmt");
				SubAmount = XML_doubleNodeTreating(listItem,"ProtectEndAmt");

				if (theNode)
				{
				theNode->Release();
				theNode=NULL;
				}

				listItem->Release();
				listItem=NULL;
			}

			if (resultList)
			{
			resultList->Release();
			resultList = NULL;
			}
		}

		if (resultList)
		{
			resultList->Release();
			resultList = NULL;
		}

		//------------------------------------------------------------------------------
		//Création de l'instrument financier
		//------------------------------------------------------------------------------

		Local_XML_TRACE("c:\\temp\\log_arm_credit.txt","Création du NTD");

		if (EndDate <= StartDate)
		{
			Local_XML_TRACE("c:\\temp\\log_arm_credit.txt","Error: EndDate <= StartDate");
			throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT, "Error: EndDate <= StartDate");
		}

		if (MezzAmount > SubAmount)
		{
			Local_XML_TRACE("c:\\temp\\log_arm_credit.txt","Error: MezzAmount > SubAmount");
			throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT, "Error: EndDate <= StartDate");
		}

		if (nbissuers == 0)
		{
			Local_XML_TRACE("c:\\temp\\log_arm_credit.txt","Error: MezzAmount > SubAmount");
			throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT, "Error: nbissuers = 0");
		}

		double* Notionals_Modified = new double[nbissuers];
		double MezzAmount_Modified = 0.;
		double SubAmount_Modified = 0.;

		ConversionMnts(nbissuers,Notionals,MezzAmount,SubAmount,
					Notionals_Modified, SubAmount_Modified, MezzAmount_Modified);

		
		CoefAdjust = FixedPayerAmount/MezzAmount_Modified;

		bool toto (isFstCpEffect) ; 
		if (isFstCpEffect && fstCpEffectDate>=StartDate) toto=false;
		ICM_Cdo2* cdo2 = new ICM_Cdo2(StartDate,
										EndDate,
										strcmp(ReferenceDate,"NULL")==0 ? 0 : &ARM_Date(ReferenceDate),
										toto?&fstCpEffectDate:0,
										spread,
										interestRule,	// intRule
										adjStartDateRule,	// adjStartDate
										SubAmount_Modified,
										SensDeLaTransaction*MezzAmount_Modified,
										//ReferenceDate,
										payFreq,
										daycount,
										qACCRUED_SETTLED,
										ccy,
										stub,
										Pf,
										DelivDays,
										DEFAULT_FRQ_DEFLEG,
										Binary,
										PayCal,
										excludeMaturity,
										INCLUDE_MATURITY // const bool& IncludeMaturity = INCLUDE_MATURITY 
										);

		if (IntegrationMethod != CREDIT_DEFAULT_VALUE)
		{
				ICM_SummitContainer SummitContainer;
				SummitContainer.itsIntegrationMethod = IntegrationMethod;
				cdo2->SetSummitContainer(&SummitContainer);
		}
		
		if (info)
		{
			cdo2->SetProportions(info);
			delete info;
		}

		cdo2->SetTradedCoef(CoefAdjust);
		cdo2->SetPorS(SensDeLaTransaction);

		//if (ccy) delete ccy;

		// cdo2->SetIssuersInfos(nbissuers,Issuers,Notionals_Modified); 
		// JLA Added: we update the tranche filles Traded Coeff. ... 
		for (int i=0;i<nbissuers;i++) 
		{
			ICM_Mez * mez = dynamic_cast<ICM_Mez*>( cdo2->GetPortfolio()->GetSecurity(i) ); 
			mez->SetTradedCoef( 
				fabs ( Notionals_Modified[i] ) / mez->GetMezzAmount() 
				) ;
		}

		
		Local_XML_TRACE("c:\\temp\\log_arm_credit.txt","Fin Création du CDO2");

		for (  i = 0; i<nbissuers; i++)
		{
			if (Issuers[i]) delete[] Issuers[i];
			Issuers[i] = NULL;
		}

		if (Issuers) delete[] Issuers;
		if (Notionals) delete[] Notionals;		
		if (Notionals_Modified) delete[] Notionals_Modified;
		if (Rate) delete[] Rate;		
		//if (PayCal) delete[] PayCal;

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

		if (XMLDoc)	{ XMLDoc->Release();}
		XMLDoc = NULL;

		if (xmlWCharText)
			delete[] xmlWCharText;
		xmlWCharText = NULL;

		return (cdo2);
	}
	catch(Exception& x)
	{
		if (theNode) theNode->Release();
		if (listItem) listItem->Release();
		if (resultList) resultList->Release();
		if (XMLDoc) XMLDoc->Release();


		if (xmlWCharText)
			delete[] xmlWCharText;
		xmlWCharText = NULL;

		hr = S_FALSE;
		char Tmp[1000];
		x.GetErrorMessage(Tmp);
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,Tmp);
	}
	catch(...)
	{		
		if (theNode) theNode->Release();
		if (listItem) listItem->Release();
		if (resultList) resultList->Release();
		if (XMLDoc) XMLDoc->Release();

		if (xmlWCharText)
			delete[] xmlWCharText;
		xmlWCharText = NULL;

		hr = S_FALSE;
		throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
         "Error in XML parsing for getting Tranche");
	}

	return NULL;
}

**/ 

// ------------------------------------
// Extract Collateral
// ------------------------------------

void ARMLOCAL_XML_COLLATERAL(const char* chaineXML, 
								  vector<string>& IssuersNames)
{
	VARIANT_BOOL bOK;
	HRESULT hr;

	MSXML2::IXMLDOMDocument *XMLDoc = NULL;
	MSXML2::IXMLDOMNodeList * resultList = NULL;
	MSXML2::IXMLDOMNodeList * resultList0 = NULL;
	MSXML2::IXMLDOMNode * listItem = NULL, * theNode = NULL;
	long nbNodes = 0;

	CCString strtmp;

	int nbissuers = 0;
	char** Issuers = NULL; 
	
	wchar_t * xmlWCharText = NULL;

	try
	{
		hr = CoInitialize(NULL); 
//		SUCCEEDED(hr) ? 0 : throw hr;
		hr = CoCreateInstance(__uuidof(MSXML2::DOMDocument30), NULL, CLSCTX_INPROC_SERVER, __uuidof(MSXML2::IXMLDOMDocument), (void**)&XMLDoc);
		SUCCEEDED(hr) ? 0 : throw hr;

		xmlWCharText = constchar2wchar(chaineXML);

		Local_XML_TRACE("c:\\temp\\log_arm_credit.txt","Lecture XML pour Matrice de Correlation");
		XMLDoc->loadXML((_bstr_t)xmlWCharText, &bOK);

		
		#ifdef _DEBUG
		CCString chtmp = "c:\\temp\\AAA_CDO_OUTPUT.xml";
		VARIANT v_chtmp;
		CCString2VARIANT (chtmp, &v_chtmp);
		XMLDoc->save(v_chtmp);
		#endif


	}
	catch(Exception& )
	{
		if (XMLDoc) XMLDoc->Release();
		hr = S_FALSE;
		throw; 


		if (xmlWCharText)
			delete[] xmlWCharText;
		xmlWCharText = NULL;
	}
	catch(...)
	{
		if (XMLDoc) XMLDoc->Release();
		hr = S_FALSE;

		if (xmlWCharText)
			delete[] xmlWCharText;
		xmlWCharText = NULL;

		throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"Pb in creating XML document for getting Tranches informations");
	}

	try
	{

		// ---------------------------------------------------------------------------------
		// Recuperation des infos sur les Issuers
		// ---------------------------------------------------------------------------------

		Local_XML_TRACE("c:\\temp\\log_arm_credit.txt","Recuperation des infos sur les Issuers");

			if (XMLDoc->selectNodes((_bstr_t)("//ASSET/ProdData/CREDSWAP/CreditEntList/CREDIT[DefpIssuer != \"FTD1\"]"), &resultList0) == S_OK)
			{
				resultList0->get_length(&nbNodes);

				if (nbNodes == 0)
				{	
					hr = S_FALSE;
					throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT, "Invalid XML string for getting ");
				}

				Issuers = new char*[nbNodes];
				nbissuers = nbNodes; 

				for (long indexNode=0 ; indexNode<nbNodes ; indexNode++)
				{
					hr=resultList0->get_item(indexNode, &listItem);

					listItem->selectSingleNode(_bstr_t("DefpIssuer"), &theNode);
					if (theNode!=NULL)
					{
					BSTR resultat = NULL;
					theNode->get_text(&resultat);

					_bstr_t ff(resultat,false);
					char * ff1=(char *)ff;

					Issuers[indexNode] = new char[1000];
					strcpy(Issuers[indexNode],ff1);
		
					theNode->Release();
					theNode=NULL;
					if (resultat) SysFreeString(resultat);
					}

					listItem->Release();
					listItem = NULL;

				}

				if (resultList0)
				{
					resultList0->Release();
					resultList0 = NULL;
				}

				if (theNode)
				{
				theNode->Release();
				theNode=NULL;
				}
	
			}

		if (resultList)
		{
		resultList->Release();
		resultList = NULL;
		}

		//------------------------------------------------------------------------------
		//Création de l'instrument financier
		//------------------------------------------------------------------------------

		for (int i = 0; i<nbissuers; i++)
		{
			IssuersNames.push_back((string)Issuers[i]);
			if (Issuers[i]) delete[] Issuers[i];
			Issuers[i] = NULL;
		}

		if (Issuers) delete[] Issuers;

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

		if (XMLDoc)	{ XMLDoc->Release();}
		XMLDoc = NULL;

		if (xmlWCharText)
			delete[] xmlWCharText;
		xmlWCharText = NULL;

	}
	catch(Exception& )
	{
		if (theNode) theNode->Release();
		if (listItem) listItem->Release();
		if (resultList) resultList->Release();
		if (XMLDoc) XMLDoc->Release();


		if (xmlWCharText)
			delete[] xmlWCharText;
		xmlWCharText = NULL;

		hr = S_FALSE;
		throw; 
	}
	catch(...)
	{		
		if (theNode) theNode->Release();
		if (listItem) listItem->Release();
		if (resultList) resultList->Release();
		if (XMLDoc) XMLDoc->Release();

		if (xmlWCharText)
			delete[] xmlWCharText;
		xmlWCharText = NULL;

		hr = S_FALSE;
		throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
         "Error in XML parsing for getting Tranche");
	}

}


// ------------------------------------
// Generic Tranche Parsing for Cdo 2
// ------------------------------------

/** 
ICM_Cdo2* ARMLOCAL_XML_CDO2(const char* chaineXML, 
								CCString& bookId, 
								CCString& custId, 
								CCString& dealId,
								ARM_ReferenceValue*fees)
{

	vector<string> UnderlyingCdos;

	ARMLOCAL_XML_COLLATERAL(chaineXML,UnderlyingCdos);

	ARM_Security** TAB = new ARM_Security*[UnderlyingCdos.size()];
	CCString bookId2;
	CCString custId2;
	CCString dealId2;
	int i=0,j=0;

	ICM_Mez* cdo = ARMLOCAL_XML_CDO(chaineXML,bookId2,custId2,dealId2,0);

	vector<string> Cdos = etoolkit_getTradeListByBook(dealId2);
	vector<string> TradeIdCdos;
	vector<string> summitTypes;

	SearchTradeIdFromList(Cdos,TradeIdCdos,summitTypes);

	for (i=0;i<UnderlyingCdos.size();i++)
	{
		CCString Result; 
		if (summitTypes[i]=="CDO") Result =etoolkit_getXMLObjectFromSummit(TradeIdCdos[i].c_str(),"CDO");
		else if (summitTypes[i]=="STCDO") Result =etoolkit_getXMLObjectFromSummit(TradeIdCdos[i].c_str(),"CDO");
		else if (summitTypes[i]=="NTD") Result =etoolkit_getXMLObjectFromSummit(TradeIdCdos[i].c_str(),"NTD");
		else if (summitTypes[i]=="FTD") Result =etoolkit_getXMLObjectFromSummit(TradeIdCdos[i].c_str(),"NTD");
		else ICMTHROW(ERR_INVALID_ARGUMENT,"ARMLOCAL_XML_CDO2: Can't load underlying "<<TradeIdCdos[i]<<" of type "<<summitTypes[i]); 

		if (summitTypes[i]=="CDO") TAB[i] = ARMLOCAL_XML_CDO(Result,bookId2,custId2,dealId2,0);	
		else if (summitTypes[i]=="STCDO") TAB[i] = ARMLOCAL_XML_CDO(Result,bookId2,custId2,dealId2,0);	
		else if (summitTypes[i]=="NTD") TAB[i] = ARMLOCAL_XML_NTD(Result,bookId2,custId2,dealId2,0);	
		else if (summitTypes[i]=="FTD") TAB[i] = ARMLOCAL_XML_NTD(Result,bookId2,custId2,dealId2,0);	
	}

	ICM_Portfolio* pf= new ICM_Portfolio(TAB,UnderlyingCdos.size());
	ICM_Parameters Parameters;

	
	//
	//	This will faild for FTD .
	int sizeprop = ((ICM_Ftd*)TAB[0])->GetProportions()->GetProportions().size();
	for (i=0;i<sizeprop;i++)
	{
		ARM_Vector* v = new ARM_Vector(UnderlyingCdos.size(),0.);
		string INDEX_name = ((ICM_Ftd*)TAB[0])->GetProportions()->GetIndexDefProb()[i];
		for (j=0;j<UnderlyingCdos.size();j++)
		{
			v->Elt(j) = ((ICM_Ftd*)TAB[0])->GetProportions()->GetProportions()[i];
		}
		Parameters.Push(v,(char*)(const char*)INDEX_name.c_str());
	}
	
	ICM_QMatrix<string>* TradeIdCdos2 = new ICM_QMatrix<string>(TradeIdCdos.size(),1);
	for (i=0;i<TradeIdCdos.size();i++)
	{
		(*TradeIdCdos2)(i,0)=TradeIdCdos[i];
	}
	Parameters.Push(TradeIdCdos2,"NAME");

	pf->SetParams(&Parameters);

	ICM_Cdo2* cdo2 = ARMLOCAL_XML_CDO2_EMPTY(chaineXML,bookId,custId,dealId,pf);
	delete[] TAB;

	if (cdo) delete cdo;

	return (cdo2);
}

**/ 
void SearchTradeIdFromList(const vector<string>& TradeList,
						   vector<string>& TradeIdList,
						   vector<string>& TradeType)
{
	MSXML2::IXMLDOMDocument* XMLDoc = NULL;
	MSXML2::IXMLDOMNodeList * resultList = NULL;
	MSXML2::IXMLDOMNode * listItem = NULL, * theNode = NULL;
	TradeIdList.resize(TradeList.size());
	TradeType.resize(TradeList.size());

	HRESULT hr;
	VARIANT_BOOL bOK;

	for (int i = 0; i < TradeList.size (); i++)
	{
		try
		{
			hr = CoInitialize(NULL); 
			hr = CoCreateInstance(__uuidof(MSXML2::DOMDocument30), NULL, CLSCTX_INPROC_SERVER, __uuidof(MSXML2::IXMLDOMDocument), (void**)&XMLDoc);
			SUCCEEDED(hr) ? 0 : throw hr;

			// JLA _bstr_t tmpChaine = (TradeList[i]).c_str();
			_bstr_t tmpChaine; 
			VariantTools::convert(TradeList[i],tmpChaine); 


			XMLDoc->loadXML(tmpChaine, &bOK);
		}
		catch(...)
		{
			if (XMLDoc) XMLDoc->Release();
			hr = S_FALSE;
			ICMTHROW(ERR_INVALID_ARGUMENT,"Pb in creating XML document for getting "<<TradeList[i]); 
		}
		std::string tradeId ; 
		TradeIdList[i]= GetStringFromXMLNode(XMLDoc,"//TradeId"); 
		TradeType[i]= GetStringFromXMLNode(XMLDoc,"//ProductType"); 
		XMLDoc->Release();
		XMLDoc = NULL;
	}
}
