
#include <ARM\libarm_local\firstToBeIncluded.h>

#include <ARM\libarm_local\arm_local_glob.h>
#include <ARM\libarm_frometk\arm_local_parsexml_nt_pricer.h>
#include <ARM\libarm_frometk\arm_local_parsexml_util.h>
 
#ifndef XML_DEFINE
#define XML_DEFINE
#import <msxml3.dll> raw_interfaces_only
using namespace MSXML2;
#endif

/** 

ICM_ModelMultiCurves* ARMLOCAL_XML_ModelMultiCurves(const char* chaineXML)
{
	VARIANT_BOOL bOK;
	VARIANT_BOOL bBinary;
	HRESULT hr;
	MSXML2::IXMLDOMDocument *XMLDoc = NULL;

	ARM_Date AsOfDate;
	ARM_Currency* ccy=NULL;
	int InterpolType = 0;

	CCString Currency;
	CCString strtmp;

	double			   binary = 0.;
	int				   NbDefCurves = 0;
	ICM_DefaultCurve**   DefaultCurves = NULL;;
	ARM_Vector RecoveryRates ; 

	ICM_Smile_Correlation* Correlation = NULL;
	ICM_CorrMatrix*   CorrelationMatrix = NULL;

	ARM_ZeroCurve* DiscountCurve = NULL;
	ARM_ZeroCurve* CpnZeroCurve = NULL;
	// double*			   Betas = NULL;	
	int				 il = 0;
	double			BetaFlat = 0.5;

	string index_ccy;
	string index_term;
	string index_name;
	double unusedspread=0.;
	bool index_spreadonly = false;

	wchar_t * xmlWCharText = NULL;

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
	catch(Exception& )
	{
		if (XMLDoc) XMLDoc->Release();
		hr = S_FALSE;

		throw; 
	}
	catch(...)
	{
		if (XMLDoc) XMLDoc->Release();
		hr = S_FALSE;

		CCString msg ((CCString)"Pb in creating XML document for ModelMultiCurves :ARMLOCAL_XML_ModelMultiCurves \n");

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
		// Création de la courbe de taux
		// ---------------------------------------------------------------------------------

		Local_XML_TRACE("c:\\temp\\log_arm_credit.txt","Création de la courbe de taux dans le MMC");

 		if (XMLDoc->selectNodes((_bstr_t)("//ASSET"), &resultList) == S_OK)
		{

			hr=resultList->get_item(0, &listItem);

			listItem->selectSingleNode(_bstr_t("Ccy"), &theNode);
			if (theNode!=NULL)
			{
				BSTR resultat = NULL;
				theNode->get_text(&resultat);

				_bstr_t ff(resultat,false);
				Currency = (CCString) (char *)ff;

				theNode->Release();
				theNode=NULL;
				if (resultat) SysFreeString(resultat);
			}

			ParseSpreadFollowTxt(listItem, unusedspread, index_ccy, index_term, index_name,index_spreadonly);

			listItem->Release();
			listItem=NULL;
		}

		if (resultList)
		{
			resultList->Release();
			resultList = NULL;
		}

		DiscountCurve = ARMLOCAL_XML_ZCPY_with_summit_stripper(chaineXML,Currency);

		if (index_spreadonly==false)
		{ CCString IndexCcy = (CCString)index_ccy.c_str();
		  CpnZeroCurve = ARMLOCAL_XML_ZCPY_with_summit_stripper(chaineXML,IndexCcy);}

		// ---------------------------------------------------------------------------------
		// Création de la matrice de correlation
		// ---------------------------------------------------------------------------------

		Correlation = ARMLOCAL_XML_CORRELATION_SMILE(chaineXML);
		CorrelationMatrix = ARMLOCAL_XML_CORRMATRIX(chaineXML);

		// ---------------------------------------------------------------------------------
		// Recuperation des DefProbCurves
		// ---------------------------------------------------------------------------------

		Local_XML_TRACE("c:\\temp\\log_arm_credit.txt","Recuperation des DefProbCurves");

		strtmp = "//CURVEHEAD/CVCHAR[MCType = \"DEFPROB\"]";	

		if (XMLDoc->selectNodes((_bstr_t)(const char*)(strtmp), &resultList) == S_OK)
		{
			resultList->get_length(&nbNodes);

			if (nbNodes == 0)
			{
				hr = S_FALSE;

				CCString msg ((CCString)"Invalid XML string for ZC \n" + chaineXML);
		
				throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
								 (const char*) msg);
			}

			NbDefCurves = nbNodes;
			DefaultCurves = new ICM_DefaultCurve*[nbNodes];
			RecoveryRates.Resize(nbNodes);
			
			for (il=0; il<nbNodes; il++)
			{	
				DefaultCurves[il] =  NULL;
				RecoveryRates[il] = 0.;
			}		

			for (long indexNode=0 ; indexNode<nbNodes ; indexNode++)
			{
				hr=resultList->get_item(indexNode, &listItem);
				if (hr==S_OK && listItem!=NULL)
				{
					_bstr_t* bIssuers = new _bstr_t("Issuer");
					listItem->selectSingleNode(*bIssuers, &theNode);
					if (theNode!=NULL)
					{
						BSTR resultat = NULL;
						theNode->get_text(&resultat);

						_bstr_t ff(resultat,false);
						char * ff1=(char *)ff;
						
						DefaultCurves[indexNode] =  ARMLOCAL_XML_DEFPROB_with_summit_stripper(chaineXML,(const char*) ff1,DiscountCurve);
						RecoveryRates[indexNode] = DefaultCurves[indexNode]->GetRecovery();

						theNode->Release();
						theNode=NULL;
						if (resultat) SysFreeString(resultat);
					}
					delete bIssuers;

					listItem->Release();
					listItem=NULL;
				}
			}
		}


		if (resultList)
		{
		resultList->Release();
		resultList=NULL;
		}


		// ---------------------------------------------------------------------------------
		// Test binary Amount
		// ---------------------------------------------------------------------------------

		strtmp = "//ASSET/ProdData/CREDSWAP";

		if (XMLDoc->selectNodes((_bstr_t)strtmp, &resultList) == S_OK)
		{
			hr=resultList->get_item(0, &listItem);
			listItem->selectSingleNode(_bstr_t("BinaryAmt"), &theNode);
			if (theNode!=NULL)
			{
			BSTR resultat1 = NULL;
			theNode->get_text(&resultat1);

			_bstr_t ff(resultat1,false);
			char * ff1=(char *)ff;
			
			if (!strcmp(ff1,"")) bBinary = VARIANT_FALSE;

			theNode->Release();
			theNode=NULL;
			if (resultat1) SysFreeString(resultat1);
			}		
		}

		if (resultList)
		{
			resultList->Release();
			resultList = NULL;
		}

		// ---------------------------------------------------------------------------------
		// Récupération du binary Amount
		// ---------------------------------------------------------------------------------

		if (bBinary == VARIANT_TRUE )
		{
		strtmp = "//ASSET/ProdData/CREDSWAP";

		if (XMLDoc->selectNodes((_bstr_t)strtmp, &resultList) == S_OK)
		{
			resultList->get_length(&nbNodes);
		
			if (nbNodes == 0)
			{
				hr = S_FALSE;

				CCString msg ((CCString)"Invalid XML for binary Amount\n" + chaineXML);
		
				throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
								 (const char*) msg);
			}

			hr=resultList->get_item(0, &listItem);

			binary = XML_doubleNodeTreating(listItem,"BinaryAmt");

			for (int i=0; i<NbDefCurves; i++)
				RecoveryRates[i] = binary;

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

		// ---------------------------------------------------------------------------------
		// Fin Récupération du binary Amount
		// ---------------------------------------------------------------------------------

		if (XMLDoc) {XMLDoc->Release();}
		XMLDoc = NULL;


		if (xmlWCharText)
			delete xmlWCharText;
		xmlWCharText = NULL;
	}
	catch(Exception& )
	{
		if (DiscountCurve)
			delete DiscountCurve;
		DiscountCurve = NULL;
		
		if (listItem) listItem->Release();
		if (resultList) resultList->Release();
		if (theNode) theNode->Release();
		if (XMLDoc) XMLDoc->Release();

		hr = S_FALSE;

		if (xmlWCharText)
			delete xmlWCharText;
		xmlWCharText = NULL;

		throw; 
	}
	catch(...)
	{
		if (DiscountCurve)
			delete DiscountCurve;
		DiscountCurve = NULL;
		
		if (listItem) listItem->Release();
		if (resultList) resultList->Release();
		if (theNode) theNode->Release();
		if (XMLDoc) XMLDoc->Release();

		if (xmlWCharText)
			delete xmlWCharText;
		xmlWCharText = NULL;

		hr = S_FALSE;

		CCString msg ((CCString)"Unknown Error for parsing DefProb Curve XML : ARMLOCAL_XML_ModelMultiCurves\n");

		Local_XML_TRACE("c:\\temp\\log_arm_credit.txt",msg);

		throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
						 (const char*) msg);
	}

	ARM_Vector* mktData = NULL;
	ARM_Vector* mktMatu = NULL;

	ARM_Currency* aCcy = NULL;

	try
	{

		Local_XML_TRACE("c:\\temp\\log_arm_credit.txt","Création du model multicurves");

		if (CorrelationMatrix == NULL)
		{
			// char** labels = new char*[NbDefCurves];
			std::vector<std::string> labels(NbDefCurves) ; 
			// Betas = new double[NbDefCurves];
			for (il=0; il<NbDefCurves; il++) 
			{ 
				labels[il] = DefaultCurves[il]->GetLabel();
				//Betas[il] = BetaFlat;
			}
			ICM_QMatrix<double> corr(NbDefCurves,NbDefCurves); 
			corr.Fill(0);
			CorrelationMatrix =  new ICM_CorrMatrix(DiscountCurve->GetAsOfDate(),
													"SERVARM",labels,corr);
			// delete [] labels;
		}

		ICM_ModelMultiCurves* modelmc = new ICM_ModelMultiCurves(NbDefCurves,
											   DefaultCurves,
											   DiscountCurve,
											   RecoveryRates,
											   CorrelationMatrix);
		if (CpnZeroCurve) {modelmc->SetCpnIRCurve(CpnZeroCurve);}
		
		if (Correlation) 
		{ 
			modelmc->SetCorrelation(Correlation);
			delete Correlation;
			Correlation = NULL;
		}
		
		for (int i=0;i<NbDefCurves;i++)
		{
			delete DefaultCurves[i];
			DefaultCurves[i] = NULL;	
		}

		if (DefaultCurves)
			delete[] DefaultCurves;
		DefaultCurves = NULL;

		// if (RecoveryRates)
		//	delete[] RecoveryRates;
		//RecoveryRates = NULL;

		// if (Betas)
		// 	delete[] Betas;
		// Betas = NULL;

		if (CorrelationMatrix)
			delete CorrelationMatrix;
		CorrelationMatrix = NULL;

		if (DiscountCurve)
			delete DiscountCurve;
		DiscountCurve= NULL;

		Local_XML_TRACE("c:\\temp\\log_arm_credit.txt","Fin Création du model multicurves");

		if (ccy)
			delete ccy;
		ccy = NULL;

		if (xmlWCharText)
			delete xmlWCharText;
		xmlWCharText = NULL;

		return (modelmc);
	}
	catch(Exception& )
	{
		if (ccy)
			delete ccy;
		ccy = NULL;

		if (DefaultCurves)
			delete[] DefaultCurves;
		DefaultCurves = NULL;

		//if (RecoveryRates)
		//	delete[] RecoveryRates;
		//RecoveryRates = NULL;

//		if (Betas)
//			delete[] Betas;
//		Betas = NULL;

		if (CorrelationMatrix)
			delete CorrelationMatrix;
		CorrelationMatrix = NULL;

		if (DiscountCurve)
			delete DiscountCurve;
		DiscountCurve= NULL;

		hr = S_FALSE;

		if (xmlWCharText)
			delete xmlWCharText;
		xmlWCharText = NULL;

		throw; 
	}
	catch(...)
	{

		if (ccy)
			delete ccy;
		ccy = NULL;

		if (DefaultCurves)
			delete[] DefaultCurves;
		DefaultCurves = NULL;

		// if (RecoveryRates)
		// 	delete[] RecoveryRates;
		// RecoveryRates = NULL;

// 		if (Betas)
// 			delete[] Betas;
// 		Betas = NULL;

		if (CorrelationMatrix)
			delete CorrelationMatrix;
		CorrelationMatrix = NULL;

		if (DiscountCurve)
			delete DiscountCurve;
		DiscountCurve= NULL;

		if (xmlWCharText)
			delete xmlWCharText;
		xmlWCharText = NULL;

		hr = S_FALSE;

		Local_XML_TRACE("c:\\temp\\log_arm_credit.txt","Unknown Error in  model multicurves object construction");

		throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
         "Error in  ZCPY construction");
	}
}
**/ 

/** 
ICM_DefaultCurveModel * ARMLOCAL_XML_DefPobModelForSpreadOption(const char* chaineXML)
{
	VARIANT_BOOL bOK;
	VARIANT_BOOL bBinary;
	HRESULT hr;
	MSXML2::IXMLDOMDocument *XMLDoc = NULL;

	ARM_Date AsOfDate;
	ARM_Currency* ccy=NULL;
	int InterpolType = 0;

	CCString Currency;
	CCString strtmp;

	ARM_ZeroCurve* DiscountCurve = NULL;
	ICM_DefaultCurve*   DefaultCurves = NULL;;
	ARM_VolCurve* VolCurve = NULL;

	wchar_t * xmlWCharText = NULL;

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
	catch(Exception& )
	{
		if (XMLDoc) XMLDoc->Release();
		hr = S_FALSE;

		throw; 
	}
	catch(...)
	{
		if (XMLDoc) XMLDoc->Release();
		hr = S_FALSE;

		CCString msg ((CCString)"Pb in creating XML document for ModelMultiCurves :ARMLOCAL_XML_ModelMultiCurves \n");

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
		// Création de la courbe de taux
		// ---------------------------------------------------------------------------------

		Local_XML_TRACE("c:\\temp\\log_arm_credit.txt","Création de la courbe de taux dans le MMC");

 		if (XMLDoc->selectNodes((_bstr_t)("//ASSET"), &resultList) == S_OK)
		{

			hr=resultList->get_item(0, &listItem);

			listItem->selectSingleNode(_bstr_t("Ccy"), &theNode);
			if (theNode!=NULL)
			{
				BSTR resultat = NULL;
				theNode->get_text(&resultat);

				_bstr_t ff(resultat,false);
				Currency = (CCString) (char *)ff;

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

		DiscountCurve = ARMLOCAL_XML_ZCPY_with_summit_stripper(chaineXML,Currency);

		// ---------------------------------------------------------------------------------
		// Recuperation des DefProbCurves
		// ---------------------------------------------------------------------------------

		Local_XML_TRACE("c:\\temp\\log_arm_credit.txt","Recuperation des DefProbCurves");

		strtmp = "//CREDSWAP//CreditEntList//CREDIT";	

		if (XMLDoc->selectNodes((_bstr_t)(const char*)(strtmp), &resultList) == S_OK)
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
			if (hr==S_OK && listItem!=NULL)
			{
				_bstr_t* bIssuers = new _bstr_t("DefpIssuer");
				listItem->selectSingleNode(*bIssuers, &theNode);
				if (theNode!=NULL)
				{
					BSTR resultat = NULL;
					theNode->get_text(&resultat);

					_bstr_t ff(resultat,false);
					char * ff1=(char *)ff;
						
					DefaultCurves =  ARMLOCAL_XML_DEFPROB_with_summit_stripper(chaineXML,(const char*) ff1);

					theNode->Release();
					theNode=NULL;
					if (resultat) SysFreeString(resultat);
				}
				delete bIssuers;

				listItem->Release();
				listItem=NULL;
			}
		}

		if (resultList)
		{
		resultList->Release();
		resultList=NULL;
		}


		// ---------------------------------------------------------------------------------
		// Test binary Amount
		// ---------------------------------------------------------------------------------

		strtmp = "//ASSET/ProdData/CREDSWAP";

		if (XMLDoc->selectNodes((_bstr_t)strtmp, &resultList) == S_OK)
		{
			hr=resultList->get_item(0, &listItem);
			listItem->selectSingleNode(_bstr_t("BinaryAmt"), &theNode);
			if (theNode!=NULL)
			{
			BSTR resultat1 = NULL;
			theNode->get_text(&resultat1);

			_bstr_t ff(resultat1,false);
			char * ff1=(char *)ff;
			
			if (!strcmp(ff1,"")) bBinary = VARIANT_FALSE;

			theNode->Release();
			theNode=NULL;
			if (resultat1) SysFreeString(resultat1);
			}		
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

		// ---------------------------------------------------------------------------------
		// Fin Récupération du binary Amount
		// ---------------------------------------------------------------------------------

		if (XMLDoc) {XMLDoc->Release();}
		XMLDoc = NULL;


		if (xmlWCharText)
			delete xmlWCharText;
		xmlWCharText = NULL;
	}
	catch(Exception& )
	{
		if (DiscountCurve)
			delete DiscountCurve;
		DiscountCurve = NULL;
		
		if (listItem) listItem->Release();
		if (resultList) resultList->Release();
		if (theNode) theNode->Release();
		if (XMLDoc) XMLDoc->Release();

		hr = S_FALSE;

		if (xmlWCharText)
			delete xmlWCharText;
		xmlWCharText = NULL;

		throw; 
	}
	catch(...)
	{
		if (DiscountCurve)
			delete DiscountCurve;
		DiscountCurve = NULL;
		
		if (listItem) listItem->Release();
		if (resultList) resultList->Release();
		if (theNode) theNode->Release();
		if (XMLDoc) XMLDoc->Release();

		if (xmlWCharText)
			delete xmlWCharText;
		xmlWCharText = NULL;

		hr = S_FALSE;

		CCString msg ((CCString)"Unknown Error for parsing DefProb Curve XML : ARMLOCAL_XML_ModelMultiCurves\n");

		Local_XML_TRACE("c:\\temp\\log_arm_credit.txt",msg);

		throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
						 (const char*) msg);
	}

	ARM_Vector* mktData = NULL;
	ARM_Vector* mktMatu = NULL;

	ARM_Currency* aCcy = NULL;

	try
	{

		Local_XML_TRACE("c:\\temp\\log_arm_credit.txt","Création du model multicurves");


		ICM_DefaultCurveModel* modelmc = new ICM_DefaultCurveModel(DefaultCurves,
																   DiscountCurve,
																   VolCurve);

	
		if (DefaultCurves) delete DefaultCurves;

 		Local_XML_TRACE("c:\\temp\\log_arm_credit.txt","Fin Création du model multicurves");

		if (ccy)
			delete ccy;
		ccy = NULL;

		if (xmlWCharText)
			delete xmlWCharText;
		xmlWCharText = NULL;

		return (modelmc);
	}
	catch(Exception& )
	{
		if (ccy)
			delete ccy;
		ccy = NULL;

		if (DefaultCurves)
			delete DefaultCurves;
		DefaultCurves = NULL;

		if (DiscountCurve)
			delete DiscountCurve;
		DiscountCurve= NULL;

		hr = S_FALSE;

		if (xmlWCharText)
			delete xmlWCharText;
		xmlWCharText = NULL;

		throw; 
	}
	catch(...)
	{

		if (ccy)
			delete ccy;
		ccy = NULL;

		if (DefaultCurves)
			delete DefaultCurves;
		DefaultCurves = NULL;

		if (DiscountCurve)
			delete DiscountCurve;
		DiscountCurve= NULL;

		if (xmlWCharText)
			delete xmlWCharText;
		xmlWCharText = NULL;

		hr = S_FALSE;

		Local_XML_TRACE("c:\\temp\\log_arm_credit.txt","Unknown Error in  model multicurves object construction");

		throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
         "Error in  ZCPY construction");
	}
}

  **/ 