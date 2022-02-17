
#include <ARM\libarm_local\firstToBeIncluded.h>

#include <ARM\libarm_local\arm_local_glob.h>
#include <ARM\libarm_frometk\arm_local_parsexml_nt_glob.h>
#include <ARM\libarm_frometk\arm_local_parsexml_util.h>
#include <ARM\libarm_frometk\PaserManagerUtilities.h>
#include <ARM\libarm_frometk\arm_local_parsexml_nt_curves.h>
#include <ARM\libarm_frometk\arm_local_parsexml_nt_security.h>
#include <ICMKernel/inst/icm_credit_index.h>

#include <glob\linalg.h>
#include <glob\dates.h>
#include <crv\zeroint.h>
#include <ccy\currency.h>
#include <atlbase.h>

#ifndef XML_DEFINE
#define XML_DEFINE
#import <msxml3.dll> raw_interfaces_only
using namespace MSXML2;
#endif

// ------------------------------------
// Generic Smile Correlation Matrix Using XML
// ------------------------------------

extern ICM_Smile_Correlation* ARMLOCAL_XML_CORRELATION_SMILE(const char* chaineXML)
{
	VARIANT_BOOL bOK;
	HRESULT hr;

	MSXML2::IXMLDOMDocument *XMLDoc = NULL;
	MSXML2::IXMLDOMNodeList * resultList = NULL;
	MSXML2::IXMLDOMNode * listItem = NULL, * theNode = NULL;
	long nbNodes = 0;

	CCString strtmp;
	char** Vector_Index = NULL;
	char** Vector_Currency = NULL;
	char** Vector_Term = NULL;

	char* Index = NULL;
	char* Currency = NULL;
	char* Term = NULL;

	// char** Vector_Index_Ccy_Term = NULL;
	std::vector<std::string> Vector_Index_Ccy_Term ; 
	int size_Vector_Index_Ccy_Term = 0;

	int size_max = 0;
	int il2=0,il3,il4=0;

	ARM_Date AsOfDate;

	ARM_Vector* vYearTerms = NULL;
	ARM_Vector* pdRate = NULL;
	ARM_Matrix* mVolatility = NULL;

	ARM_Vector** vStrikes = NULL;
	vector<ARM_VolCurve*> Vector_VolCurve;
	vector<const ARM_VolCurve*> Vector_VolCurve2;
	wchar_t * xmlWCharText = NULL;

	char SmilProp[500];
	memset(SmilProp,'\0',sizeof(char)*500);

	char** Labels_Prop = NULL;

	ARM_Vector* vProp = NULL;
	ARM_Vector* vProp_sort = NULL;

	ARM_Vector* vCorr_Low = NULL;
	ARM_Vector* vCorr_Low_sort = NULL;

	ARM_Vector* vCorr_Up = NULL;
	ARM_Vector* vCorr_Up_sort = NULL;

	ICM_Credit_Index** vIndex;
	vector<const ICM_Credit_Index*> vIndex_sort;

	ICM_Smile_Correlation* Corr = NULL;

	try
	{
		hr = CoInitialize(NULL); 
//		SUCCEEDED(hr) ? 0 : throw hr;
		hr = CoCreateInstance(__uuidof(MSXML2::DOMDocument30), NULL, CLSCTX_INPROC_SERVER, __uuidof(MSXML2::IXMLDOMDocument), (void**)&XMLDoc);
		SUCCEEDED(hr) ? 0 : throw hr;

		xmlWCharText = constchar2wchar(chaineXML);

		Local_XML_TRACE("c:\\temp\\log_arm_credit.txt","Lecture du XML du smile de correlation");
		XMLDoc->loadXML((_bstr_t)xmlWCharText, &bOK);
	}
	catch(...)
	{

		if (xmlWCharText)
			delete xmlWCharText;
		xmlWCharText = NULL;

		if (XMLDoc) XMLDoc->Release();
		hr = S_FALSE;
		throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"Pb in creating XML document for smile correlation");
	}
	try
	{
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
					delete xmlWCharText;
				xmlWCharText = NULL;

				return NULL;
			}

			CvtStrProp(SmilProp,Labels_Prop, vProp,vCorr_Low,vCorr_Up,vIndex);

		}

		// ------------------------------------------------------------------------
		// 1ere passe pour déterminer les triplets indices/devises/maturités 
		// ------------------------------------------------------------------------
		strtmp = "//CURVEHEAD/CVCHAR[MCType = \"SMILE\"]";

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

			Vector_Index = new char*[nbNodes];
			Vector_Currency = new char*[nbNodes];
			Vector_Term = new char*[nbNodes];

			size_max = 0;  

			for (long indexNode=0 ; indexNode<nbNodes ; indexNode++)
			{
				hr=resultList->get_item(indexNode, &listItem);
				if (hr==S_OK && listItem!=NULL)
				{
					AsOfDate = GetDateFromXMLNode(listItem,"AsOfDate");

					listItem->selectSingleNode((_bstr_t)"Ccy",&theNode);
					if (theNode!=NULL)
					{
						BSTR resultat = NULL;
						theNode->get_text(&resultat);	//printf("resultat : %S\n", resultat);

						_bstr_t ff(resultat,false);
						Currency = new char[20];
						strcpy(Currency,(char *)ff);

						theNode->Release();
						theNode=NULL;
						if (resultat) SysFreeString(resultat);
					}

					listItem->selectSingleNode((_bstr_t)"dmIndex",&theNode);
					if (theNode!=NULL)
					{
						BSTR resultat = NULL;
						theNode->get_text(&resultat);	//printf("resultat : %S\n", resultat);

						_bstr_t ff(resultat,false);
						Index = new char[20];
						strcpy(Index,(char *)ff);

						theNode->Release();
						theNode=NULL;
						if (resultat) SysFreeString(resultat);
					}

					listItem->selectSingleNode((_bstr_t)"Term",&theNode);
					if (theNode!=NULL)
					{
						BSTR resultat = NULL;
						theNode->get_text(&resultat);	//printf("resultat : %S\n", resultat);

						_bstr_t ff(resultat,false);
						Term = new char[20];
						strcpy(Term,(char *)ff);

						theNode->Release();
						theNode=NULL;
						if (resultat) SysFreeString(resultat);
					}

					listItem->Release();
					listItem=NULL;
				}

				bool flag = false;
					
				for (int il=0; il<size_max; il++)
				{
					if ((strcmp(Currency,Vector_Currency[il])==NULL) &&
						(strcmp(Index,Vector_Index[il])==NULL) && 
						(strcmp(Term,Vector_Term[il])==NULL))
					{	  
						  flag = true;
						  break;
					}
				}

				if (!flag) 
				{
					Vector_Currency[size_max] = Currency;
					Vector_Index[size_max] = Index;
					Vector_Term[size_max] = Term;
					size_max++;
				}
				else
				{
					delete[] Currency;
					delete[] Index;
					delete[] Term;
				}
			
				// End Boucle For ----------------------------------------------------
			}
		}

	if (listItem)
	{
		listItem->Release();
		listItem = NULL;
	}

	if (theNode)
	{
		theNode = NULL;
		theNode->Release();
	}
	if (resultList)
	{
		resultList->Release();
		resultList = NULL;
	}


	// ------------------------------------------------------------------------
	// 2eme passe pour déterminer les courbes de vol
	// ------------------------------------------------------------------------

	Vector_Index_Ccy_Term.resize(size_max); 
	// if (Vector_Index_Ccy_Term== NULL) 
	// 	Vector_Index_Ccy_Term = new char*[size_max];
	
	size_Vector_Index_Ccy_Term = size_max;

	vStrikes = new ARM_Vector*[size_Vector_Index_Ccy_Term];

	Vector_VolCurve.resize(size_Vector_Index_Ccy_Term);

	for (il2=0; il2<size_max; il2++)
	{
		Vector_VolCurve[il2]=NULL;
		// Vector_Index_Ccy_Term[il2] = new char[255];
		// memset(Vector_Index_Ccy_Term[il2],'\0',sizeof(char)*255);
		// sprintf(Vector_Index_Ccy_Term[il2],"%s_%s_%s",Vector_Index[il2],Vector_Currency[il2],Vector_Term[il2]);
		std::stringstream sstr; sstr<<Vector_Index[il2]<<"_"<<Vector_Currency[il2]<<"_"<<Vector_Term[il2];
		sstr>> Vector_Index_Ccy_Term[il2]; 
	}		

	for (il2=0; il2<size_max; il2++)
	{

		strtmp = (CCString)"//CURVEHEAD[CVCHAR/MCType = \"SMILE\" and CVCHAR/Ccy = \"" + 
			(CCString) Vector_Currency[il2] + 
			"\" and CVCHAR/dmIndex = \"" + 
			(CCString) Vector_Index[il2] + 
			(CCString) "\" and CVCHAR/Term = \"" + 
			(CCString) Vector_Term[il2] + (CCString) "\"]";
		
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

			vStrikes[il2] = new ARM_Vector(nbNodes,0.);

			for (long indexNode=0 ; indexNode<nbNodes ; indexNode++)
			{
				hr=resultList->get_item(indexNode, &listItem);
				if (hr==S_OK && listItem!=NULL)
				{
					double strike = XML_doubleNodeTreating(listItem,"COM2");
					strike /= 100.;
					vStrikes[il2]->InitElt(indexNode,strike);
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

	}
		(*vStrikes[il2]) = vStrikes[il2]->Sort();
	}
	}
	catch(Exception& )
	{
		if (theNode) theNode->Release();
		if (listItem) listItem->Release();
		if (resultList) resultList->Release();
		if (XMLDoc) XMLDoc->Release();

		hr = S_FALSE;

		if (xmlWCharText)
			delete xmlWCharText;
		xmlWCharText = NULL;

		throw; 
	}
	catch(...)
	{		
		if (theNode) theNode->Release();
		if (listItem) listItem->Release();
		if (resultList) resultList->Release();
		if (XMLDoc) XMLDoc->Release();

		if (xmlWCharText)
			delete xmlWCharText;
		xmlWCharText = NULL;

		hr = S_FALSE;
		throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
         "Error in XML parsing for getting Tranche");
	}


	// ---------------------------------------------------------------------------------
	// Recuperation des Dates & des rates pour chaque strike & chaque index
	// ---------------------------------------------------------------------------------

try{

for (il2=0; il2<size_Vector_Index_Ccy_Term; il2++)
	for (il4=0; il4<vStrikes[il2]->GetSize(); il4++)
	{	

	char TabT[50];
	memset(TabT,'\0',sizeof(char)*50);
	sprintf(TabT,"%f",100.*vStrikes[il2]->Elt(il4));

		strtmp = (CCString)"//CURVEHEAD/ENTITYLIST/CURVE[../../CVCHAR/MCType = \"SMILE\" and ../../CVCHAR/Ccy = \"" + 
				(CCString) Vector_Currency[il2] + 
				"\" and ../../CVCHAR/dmIndex = \"" + 
				(CCString) Vector_Index[il2] + 
				(CCString) "\" and ../../CVCHAR/Term = \"" + 
				(CCString) Vector_Term[il2] + (CCString) "\" and ../../COM2 = " + 
				(CCString) TabT + (CCString) "]";

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

			pdRate = new ARM_Vector(nbNodes,0.);
			vYearTerms = new ARM_Vector(nbNodes,0.);

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
						/*if (strcmp(Vector_Currency[il2],"USD")==NULL)
							endDate.NextBusinessDay(1,"GBP");
						else
							endDate.AdjustToBusDate(Vector_Currency[il2],K_FOLLOWING);*/

						double maturite = (endDate - AsOfDate) / 365.;

						vYearTerms->InitElt(indexNode,maturite);
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
						pdRate->InitElt(indexNode,rate);
						theNode->Release();
						theNode=NULL;
						if (resultat) SysFreeString(resultat);
					}

					listItem->Release();
					listItem=NULL;
				}
			}

		const ARM_VolCurve* Active_VolCurve = NULL;
		int indice_actif = 0;


		for (il3=0; il3<size_Vector_Index_Ccy_Term; il3++)
		{
			char TTab[255];
			memset(TTab,'\0',sizeof(char)*255);
			sprintf(TTab,"%s_%s_%s",Vector_Index[il2],Vector_Currency[il2],Vector_Term[il2]);

			if (strcmp(Vector_Index_Ccy_Term[il3].c_str(),TTab)==NULL)
			{
				Active_VolCurve = Vector_VolCurve[il3];
				indice_actif = il3;
				break;
			}
		}

		if (Active_VolCurve)
			mVolatility = Active_VolCurve->GetVolatilities();
		else
			mVolatility = new ARM_Matrix(vYearTerms->GetSize(),vStrikes[indice_actif]->GetSize());

		for (il3=0; il3<pdRate->GetSize(); il3++)
			mVolatility->Elt(il3,il4) = 100.*pdRate->Elt(il3);

		if (Active_VolCurve == NULL)
		{
			ARM_Currency* Cy = new ARM_Currency(Vector_Currency[il2]);
			Vector_VolCurve[indice_actif] = new ARM_VolLInterpol( AsOfDate,(ARM_Vector*)vYearTerms->Clone(),
												vStrikes[indice_actif], mVolatility, 1, K_SMILE_VOL,Cy);	

			Vector_VolCurve[indice_actif]->SetIndexName((char*)(const char*)Vector_Index[il2]);
			if (Cy) delete Cy;
		}
		Vector_VolCurve2.resize(Vector_VolCurve.size());
		for (int j= 0; j<Vector_VolCurve.size(); j++) { Vector_VolCurve2[j] = const_cast<ARM_VolCurve*>(Vector_VolCurve[j]);}
		if (pdRate) 
			delete pdRate;
		pdRate = NULL;

		if (vYearTerms) 
			delete vYearTerms;
		vYearTerms = NULL;

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

		}
		}

		//On ordonne les proportions
		vProp_sort = new ARM_Vector(size_Vector_Index_Ccy_Term,0.);
		vCorr_Low_sort = new ARM_Vector(size_Vector_Index_Ccy_Term,0.);
		vCorr_Up_sort = new ARM_Vector(size_Vector_Index_Ccy_Term,0.);
		
		if (vIndex)
			vIndex_sort.resize(vProp_sort->GetSize());

		for (il2=0; il2<size_Vector_Index_Ccy_Term; il2++)
		for (int y=0; y<vProp->GetSize(); y++)
		{
			if (strcmp(Vector_Index_Ccy_Term[il2].c_str(),Labels_Prop[y]) == NULL)
			{
				vProp_sort->InitElt(il2,vProp->Elt(y));
				vCorr_Low_sort->InitElt(il2,vCorr_Low->Elt(y));
				vCorr_Up_sort->InitElt(il2,vCorr_Up->Elt(y));
				if (vIndex) vIndex_sort[il2] = vIndex[y];
			}
		}

// FIXMEFRED: mig.vc8 (30/05/2007 17:49:37):cast					
		if (vIndex)
			Corr = new ICM_Smile_Correlation(AsOfDate,
											 "SERVARM",	
											 // size_Vector_Index_Ccy_Term, 
											&(*Vector_VolCurve2.begin()), 
											 Vector_Index_Ccy_Term,
											 vProp_sort,
											&(*vIndex_sort.begin()));
		else
			Corr = new ICM_Smile_Correlation(AsOfDate,
											 "SERVARM",	
											 // size_Vector_Index_Ccy_Term, 
											&(*Vector_VolCurve2.begin()), 
											 Vector_Index_Ccy_Term,
											 vProp_sort,
											 vCorr_Low_sort,
											 vCorr_Up_sort);

		if (vIndex)
		{
			for (il2=0; il2<vProp_sort->GetSize(); il2++)
			{
			if (vIndex[il2]) delete vIndex[il2];
			vIndex[il2] = NULL;
			}
			delete[] vIndex;
			//delete[] vIndex_sort;
		}

		if (Labels_Prop)
		{
		for (int y2=0; y2<vProp->GetSize(); y2++)
			delete[] Labels_Prop[y2];
		}

		if (Labels_Prop) delete[] Labels_Prop;
		if (vProp) delete vProp;
		if (vProp_sort) delete vProp_sort;
		if (vCorr_Low_sort) delete vCorr_Low_sort;
		if (vCorr_Up_sort) delete vCorr_Up_sort;
		if (vCorr_Low) delete vCorr_Low;
		if (vCorr_Up) delete vCorr_Up;

		for (il2=0; il2<size_Vector_Index_Ccy_Term; il2++)
		{
			for ( int i= 0; i< Vector_VolCurve.size() ; i++) delete Vector_VolCurve[i];
			// delete[] Vector_Index_Ccy_Term[il2];
			delete[] Vector_Index[il2];
			delete[] Vector_Currency[il2];
			delete[] Vector_Term[il2];
		}

		for ( int i= 0; i< Vector_VolCurve.size() ; i++) delete Vector_VolCurve[i];
		// if (Vector_Index_Ccy_Term) delete[] Vector_Index_Ccy_Term;
		if (Vector_Index) delete[] Vector_Index;
		if (Vector_Currency) delete[] Vector_Currency;
		if (Vector_Term) delete[] Vector_Term;
		if (vStrikes) delete[] vStrikes;


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

		if (xmlWCharText)
			delete xmlWCharText;
		xmlWCharText = NULL;

		return (Corr);
	}
	catch(Exception& )
	{
		if (theNode) theNode->Release();
		if (listItem) listItem->Release();
		if (resultList) resultList->Release();
		if (XMLDoc) XMLDoc->Release();

		hr = S_FALSE;

		if (xmlWCharText)
			delete xmlWCharText;
		xmlWCharText = NULL;

		throw; 
	}
	catch(...)
	{		
		if (theNode) theNode->Release();
		if (listItem) listItem->Release();
		if (resultList) resultList->Release();
		if (XMLDoc) XMLDoc->Release();

		if (xmlWCharText)
			delete xmlWCharText;
		xmlWCharText = NULL;

		hr = S_FALSE;
		throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
         "Error in XML parsing for getting Tranche");
	}


	return NULL;
}





// ------------------------------------
// Generic Correlation Matrix Using XML
// ------------------------------------

extern ICM_CorrMatrix* ARMLOCAL_XML_CORRMATRIX(const char* chaineXML)
{
	VARIANT_BOOL bOK;
	HRESULT hr;

	MSXML2::IXMLDOMDocument *XMLDoc = NULL;
	MSXML2::IXMLDOMNodeList * resultList = NULL;
	MSXML2::IXMLDOMNode * listItem = NULL, * theNode = NULL;
	long nbNodes = 0;

	int nbissuers = 0;
	char** Issuers = NULL; 
	// double** correlation = NULL;
	ICM_QMatrix<double> correlation; 
	long indexNode = 0;
	ICM_CorrMatrix* Corr = NULL;

	std::vector<std::string> IssuersTMP; 
	std::vector<std::string> IssuersTMP1; 
	std::vector<std::string> IssuersTMP2; 
	// char** IssuersTMP = NULL;
	// char** IssuersTMP1 = NULL;
	// char** IssuersTMP2 = NULL;
	int* StatusCorr = NULL;

	ARM_Date AsOfDate; //to modify

	wchar_t * xmlWCharText = NULL;

	int k =0;	

	try
	{
		hr = CoInitialize(NULL); 
//		SUCCEEDED(hr) ? 0 : throw hr;
		hr = CoCreateInstance(__uuidof(MSXML2::DOMDocument30), NULL, CLSCTX_INPROC_SERVER, __uuidof(MSXML2::IXMLDOMDocument), (void**)&XMLDoc);
		SUCCEEDED(hr) ? 0 : throw hr;

		xmlWCharText = constchar2wchar(chaineXML);

		Local_XML_TRACE("c:\\temp\\log_arm_credit.txt","Lecture du XML de la matrice de correlation");
		XMLDoc->loadXML((_bstr_t)xmlWCharText, &bOK);

	}
	catch(...)
	{

		if (xmlWCharText)
			delete xmlWCharText;
		xmlWCharText = NULL;

		if (XMLDoc) XMLDoc->Release();
		hr = S_FALSE;
		throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"Pb in creating XML document for Correlation Matrix");
	}

	try
	{

			// ---------------------------------------------------------------------------------
			// Recuperation des infos sur les Correlation
			// ---------------------------------------------------------------------------------
			Local_XML_TRACE("c:\\temp\\log_arm_credit.txt","Initialisation des matrices pour la correlation");

			if (XMLDoc->selectNodes((_bstr_t)("//ASSET/ProdData/CREDSWAP/CreditEntList/CREDIT"), &resultList) == S_OK)
			{
				resultList->get_length(&nbNodes);

				if (nbNodes == 0)
				{	
					hr = S_FALSE;
					throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT, "Invalid XML string for getting nbissuers");
				}

				nbissuers = nbNodes;

			}

			if (resultList)
			{
			resultList->Release();
			resultList = NULL;
			}

			if (XMLDoc->selectNodes((_bstr_t)("//CVLLIST/CURVEHEAD"), &resultList) == S_OK)
			{
				resultList->get_length(&nbNodes);

				if (nbNodes == 0)
				{	
					hr = S_FALSE;
					throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT, "Invalid XML string for getting ");
				}

				//Création et initialisation de la matrice de correlation
				// IssuersTMP = new char*[nbissuers];
				IssuersTMP.resize(nbissuers); 
				IssuersTMP1.resize(nbissuers); 
				IssuersTMP2.resize(nbissuers); 
				// IssuersTMP1 = new char*[nbNodes];
				// IssuersTMP2 = new char*[nbNodes];
				int il1=0;

				//memset(IssuersTMP,'\0',nbissuers*sizeof(char*));
				// memset(IssuersTMP1,'\0',nbNodes*sizeof(char*));
				// memset(IssuersTMP2,'\0',nbNodes*sizeof(char*));

				correlation.Resize(nbissuers,nbissuers) ; 
				// correlation = new double*[nbissuers];
				// memset(correlation,'\0',sizeof(double*)*nbissuers);
				StatusCorr = new int[nbNodes];
				memset(StatusCorr,'\0',sizeof(int*)*nbNodes);

				for (il1=0; il1<nbissuers; il1++)
				{
					// correlation[il1] = new double[nbissuers];
					// memset(correlation[il1],'\0',nbissuers*sizeof(double));
					correlation(il1,il1)=1.;
				}

				int pos = 0;
				bool alreadyin = false;

				//----------------------------------------------------------
				// On crée la liste d'issuers initiaux
				//----------------------------------------------------------
				Local_XML_TRACE("c:\\temp\\log_arm_credit.txt","Chargement des issuers 1/2 pour la correlation");

				for (indexNode=0 ; indexNode<nbNodes ; indexNode++)
				{
					hr=resultList->get_item(indexNode, &listItem);

					listItem->selectSingleNode(_bstr_t("COM1"), &theNode);
					if (theNode!=NULL)
					{
					BSTR resultat2 = NULL;
					theNode->get_text(&resultat2);

					_bstr_t ff(resultat2,false);
					char * ff1=(char *)ff;

					char TTMP[1000];
					strcpy(TTMP,ff1);
					TTMP[8]='\0';	

					if (!strcmp(TTMP,"ISSRCORR")) 
					{
						alreadyin = false;
						char* issuer = ExtractCorr(ff1);
						StatusCorr[indexNode] = 1;

						IssuersTMP1[indexNode] = issuer; 
						// new char[1000];
						// strcpy(IssuersTMP1[indexNode],issuer);

						for (int il = 0; il<pos; il++) //On check si l'issuer n'est pas déjà présent dans la liste
						{
							if (!IssuersTMP.empty())
								if (!IssuersTMP[il].empty())
									if (!strcmp(issuer,IssuersTMP[il].c_str()))
									{
									alreadyin = true;	
									break;
									}
						}

						if (alreadyin)
							delete[] issuer;
						else
						{
							IssuersTMP[pos] = issuer;
							pos++;
						}
					}
					theNode->Release();
					theNode=NULL;
					if (resultat2) SysFreeString(resultat2);
					}

					if (listItem)
					{
						listItem->Release();
						listItem = NULL;
					}

				}

				//----------------------------------------------------------
				// On compléte la liste d'issuers si besoin avec COM2
				//----------------------------------------------------------
				Local_XML_TRACE("c:\\temp\\log_arm_credit.txt","Chargement des issuers 2/2 pour la correlation");

				for (indexNode=0 ; indexNode<nbNodes ; indexNode++)
				{
					hr=resultList->get_item(indexNode, &listItem);

					listItem->selectSingleNode(_bstr_t("COM2"), &theNode);
					if (theNode!=NULL)
					{
					BSTR resultat3 = NULL;
					theNode->get_text(&resultat3);

					_bstr_t ff(resultat3,false);
					char * ff1=(char *)ff;

					if (StatusCorr[indexNode]==1)
					{
						alreadyin = false;

						IssuersTMP2[indexNode] =ff1 ;

						for (int il = 0; il<pos; il++) //On check si l'issuer n'est pas déjà présent dans la liste
						{
							if (!IssuersTMP.empty())
								if (!IssuersTMP[il].empty())
									if (!strcmp(ff1,IssuersTMP[il].c_str()))
									{
									alreadyin = true;	
									break;
									}
						}

						if (alreadyin==false)
						{
							IssuersTMP[pos] =ff1; 
							pos++;
						}

					}
					theNode->Release();
					theNode=NULL;
					if (resultat3) SysFreeString(resultat3);
					}

					if (listItem)
					{
						listItem->Release();
						listItem = NULL;
					}


				}

				//------------------------------------------------------------------------------
				//Création de la matrice de correlation
				//------------------------------------------------------------------------------
				Local_XML_TRACE("c:\\temp\\log_arm_credit.txt","Création de l'objet CorrMatrix");

				if (nbissuers <= 0)
				{
					Local_XML_TRACE("c:\\temp\\log_arm_credit.txt","Error: nbissuers <= 0");
					throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT, "Error: nbissuers <= 0");
				}

				if (IssuersTMP.empty())
				{
					Local_XML_TRACE("c:\\temp\\log_arm_credit.txt","Error: Unable to find issuers");
					throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT, "Error: Unable to find issuers");
				}

				if (IssuersTMP[0].empty() == NULL)
				{
					Local_XML_TRACE("c:\\temp\\log_arm_credit.txt","Error: Unable to find issuers");
				}
				else
					Corr = new 	ICM_CorrMatrix(AsOfDate, //to modify
											 "SERVARM",
											 IssuersTMP,correlation
											 );

				//----------------------------------------------------------
				// Allimentation de la matrice de correlation
				//----------------------------------------------------------

				if (Corr)
				{

				for (indexNode=0 ; indexNode<nbNodes ; indexNode++)
				{
					hr=resultList->get_item(indexNode, &listItem);

					if (StatusCorr[indexNode]==1)
					{
					Corr->SetCorrelation(IssuersTMP1[indexNode] ,IssuersTMP2[indexNode] ,XML_doubleNodeTreating(listItem,"ENTITYLIST/CURVE/Rate"));
					}

					listItem->Release();
					listItem=NULL;

				}

				}

				Local_XML_TRACE("c:\\temp\\log_arm_credit.txt","Fin Création de l'objet CorrMatrix");

			}


		//Suppression des variables temporaires
		for (k=0; k <nbissuers; k++)
		{
			// if (IssuersTMP[k])
			// 	delete[] IssuersTMP[k];

			// if (correlation[k])
			// 	delete[] correlation[k];
		}

		//Suppression des variables temporaires
		for (k=0; k <nbNodes; k++)
		{
			// if (IssuersTMP1[k])
			// 	delete[] IssuersTMP1[k];

			// if (IssuersTMP2[k])
			// 	delete[] IssuersTMP2[k];
		}

		//if (IssuersTMP) delete[] IssuersTMP;
		// if (IssuersTMP1) delete[] IssuersTMP1;
		// if (IssuersTMP2) delete[] IssuersTMP2;
		// if (correlation) delete[] correlation;
		if (StatusCorr) delete[] StatusCorr; 

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

		if (XMLDoc) {XMLDoc->Release();}
		XMLDoc = NULL;

		if (xmlWCharText)
			delete xmlWCharText;
		xmlWCharText = NULL;

		return (Corr);
	}
	catch(Exception& )
	{
		if (theNode) theNode->Release();
		if (listItem) listItem->Release();
		if (resultList) resultList->Release();
		if (XMLDoc) XMLDoc->Release();

		hr = S_FALSE;

		if (xmlWCharText)
			delete xmlWCharText;
		xmlWCharText = NULL;

		throw; 
	}
	catch(...)
	{		
		if (theNode) theNode->Release();
		if (listItem) listItem->Release();
		if (resultList) resultList->Release();
		if (XMLDoc) XMLDoc->Release();

		if (xmlWCharText)
			delete xmlWCharText;
		xmlWCharText = NULL;

		hr = S_FALSE;
		throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
         "Error in XML parsing for getting Tranche");
	}

	return NULL;
}


extern ICM_Matrix<ARM_Vector>* ARMLOCAL_XML_CASHFLOWS(const char* chaineXML, 
										const CCString& Path, 
										CCString& bookName, 
										CCString& structureId, 
										CCString& custId, 
										CCString& dealId)
{
	VARIANT_BOOL bOK;
	HRESULT hr;

	MSXML2::IXMLDOMDocument *XMLDoc = NULL;
	MSXML2::IXMLDOMNodeList * resultList = NULL;
	MSXML2::IXMLDOMNode * listItem = NULL, * theNode = NULL;
	long nbNodes = 0;
	int k =0;	
	ICM_Matrix<ARM_Vector>* Matrix = NULL;

	wchar_t * xmlWCharText = NULL;

	vector<string> CHAMPS = SUMMIT_NAMES::list_names();;
	int size = SUMMIT_NAMES::nbParameters();

	bookName = "";
	structureId = "";
	custId = "";
	dealId ="";

	try
	{
		hr = CoInitialize(NULL); 
//		SUCCEEDED(hr) ? 0 : throw hr;
		hr = CoCreateInstance(__uuidof(MSXML2::DOMDocument30), NULL, CLSCTX_INPROC_SERVER, __uuidof(MSXML2::IXMLDOMDocument), (void**)&XMLDoc);
		SUCCEEDED(hr) ? 0 : throw hr;

		xmlWCharText = constchar2wchar(chaineXML);

		Local_XML_TRACE("c:\\temp\\log_arm_credit.txt","Lecture XML pour Infos Credit");
		XMLDoc->loadXML((_bstr_t)xmlWCharText, &bOK);

	}
	catch(Exception& )
	{
		if (XMLDoc) XMLDoc->Release();
		hr = S_FALSE;

		if (xmlWCharText)
			delete xmlWCharText;
		xmlWCharText = NULL;

		throw; 
	}
	catch(...)
	{
		if (xmlWCharText)
			delete xmlWCharText;
		xmlWCharText = NULL;

		if (XMLDoc) XMLDoc->Release();
		hr = S_FALSE;
		throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"Pb in creating XML document for getting Tranches informations");
	}

	try
	{
	long nbNodes;

	CCString tmpChaine;

	Matrix = new ICM_Matrix<ARM_Vector>();

	tmpChaine = Path;

	if ( XMLDoc->selectNodes(_bstr_t((const char *)(tmpChaine)), &resultList) == S_OK )
	{
		resultList->get_length(&nbNodes);

		if ( nbNodes == 0 )
		{
			hr = S_FALSE;

			CCString msg((CCString)"Invalid XML string for getting flows \n");

			throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
							 (char*) msg);
		}

		for ( int indexNode = 0; indexNode<nbNodes; indexNode++)
		{
			hr = resultList->get_item(indexNode, &theNode);

			if( hr == S_OK && theNode != NULL)
			{
				for (k=0;k<size;k++)
				{
					CCString TMP = (CCString)CHAMPS[k].c_str();
					double value = 0.;

					if ((Matrix->GetColVect((char*)CHAMPS[k].c_str())==NULL) && (ExistXMLNode((void*)theNode,CHAMPS[k])))
					{
						ARM_Vector* vector = new ARM_Vector(nbNodes,0.);
						Matrix->Push(vector,(char*)CHAMPS[k].c_str());

						value = GetDoubleFromXMLNode (theNode,(char*)CHAMPS[k].c_str());

						if (TMP.Contain("Date") && value)
							value = GetDateFromXMLNode (theNode,(char*)CHAMPS[k].c_str()).GetJulian();

						Matrix->GetColVect((char*)CHAMPS[k].c_str())->Elt(indexNode) = value;
					}
					else if (Matrix->GetColVect((char*)CHAMPS[k].c_str())!=NULL)
					{
						value = GetDoubleFromXMLNode (theNode,(char*)CHAMPS[k].c_str());

						if (TMP.Contain("Date") && value)
							value = GetDateFromXMLNode (theNode,(char*)CHAMPS[k].c_str()).GetJulian();
						
						Matrix->GetColVect((char*)CHAMPS[k].c_str())->Elt(indexNode) = value;
					}
				}

				theNode->Release();
				theNode = NULL;
			}
			else
			{
				throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
							 "Probleme in get_item \n");
			}
		}

/*		if (indexNode >= nbNodes)
		{
			throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
							 "First Start Date after AsOf not Found \n");
		}
*/
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
		
        if (XMLDoc)
		{
			XMLDoc->Release();
			XMLDoc = NULL;
		}
	}

		delete xmlWCharText; xmlWCharText=NULL; 
		return Matrix;
	}
	catch(Exception& )
	{
		if (theNode) theNode->Release();
		if (listItem) listItem->Release();
		if (resultList) resultList->Release();
		if (XMLDoc) XMLDoc->Release();

		if (xmlWCharText)
			delete xmlWCharText;
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
			delete xmlWCharText;
		xmlWCharText = NULL;

		hr = S_FALSE;
		throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
         "Error in XML parsing for getting Tranche");
	}

	return NULL;
}



/** 
extern ICM_Matrix<ARM_Vector>* ARMLOCAL_XML_PRICING_CMCDS(const char* chaineXML, 
											  const CCString& Path, 
											  CCString& bookName, 
											  CCString& structureId, 
											  CCString& custId, 
											  CCString& dealId)
{
	VARIANT_BOOL bOK;
	HRESULT hr;

	MSXML2::IXMLDOMDocument *XMLDoc = NULL;
	MSXML2::IXMLDOMNodeList * resultList = NULL;
	MSXML2::IXMLDOMNode * listItem = NULL, * theNode = NULL;
	long nbNodes = 0;
	int k =0;	
	ICM_Matrix<ARM_Vector>* Matrix = NULL;

	wchar_t * xmlWCharText = NULL;

	vector<string> CHAMPS = SUMMIT_NAMES::list_names();;
	int size = SUMMIT_NAMES::nbParameters();

	bookName = "";
	structureId = "";
	custId = "";
	dealId ="";

	string FwdDefprob;
	string DiscDefprob;
	string Currency;

	ICM_Matrix<ARM_Vector>* CFMatrix = NULL;
	ARM_ZeroCurve* ZcCurve = NULL;
	ICM_DefaultCurve* pFwdDefprob = NULL; 
	ICM_DefaultCurve* pDiscDefprob = NULL; 
	ARM_VolCurve* VolFwdDefprob = NULL;
	ICM_Credit_Index* index=NULL;

	double Tenor=0.;
	double Floor=0.;
	double Cap=0.;
	string temp;

	ARM_Date BeginCMDate;
	ARM_Date EndCMDate;
	double FixedStartRate = 0.;
	int OnDefAccrued = 0;
	int daycount = 0;
	ARM_Vector* VBasis = NULL;

	try
	{
		hr = CoInitialize(NULL); 
//		SUCCEEDED(hr) ? 0 : throw hr;
		hr = CoCreateInstance(__uuidof(MSXML2::DOMDocument30), NULL, CLSCTX_INPROC_SERVER, __uuidof(MSXML2::IXMLDOMDocument), (void**)&XMLDoc);
		SUCCEEDED(hr) ? 0 : throw hr;

		xmlWCharText = constchar2wchar(chaineXML);

		Local_XML_TRACE("c:\\temp\\log_arm_credit.txt","Lecture XML pour Infos Credit");
		XMLDoc->loadXML((_bstr_t)xmlWCharText, &bOK);

	}
	catch(Exception& )
	{
		if (XMLDoc) XMLDoc->Release();
		hr = S_FALSE;

		if (xmlWCharText)
			delete xmlWCharText;
		xmlWCharText = NULL;

		throw; 
	}
	catch(...)
	{
		if (xmlWCharText)
			delete xmlWCharText;
		xmlWCharText = NULL;

		if (XMLDoc) XMLDoc->Release();
		hr = S_FALSE;
		throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"Pb in creating XML document for getting Tranches informations");
	}

	try
	{
	long nbNodes;

	CFMatrix = ARMLOCAL_XML_CASHFLOWS(chaineXML,"//cBLBLST_I",bookName,custId,dealId,structureId);

	if ( XMLDoc->selectNodes(_bstr_t("//ProdData/cBLOB_I"), &resultList) == S_OK )
	{
		resultList->get_length(&nbNodes);

		if ( nbNodes == 0 )
		{
			hr = S_FALSE;

			if (resultList)
			{
			   resultList->Release();
			
			   resultList = NULL;
			}
		}

		hr = resultList->get_item(0, &theNode);

		if ( hr == S_OK && theNode != NULL )
		{
		FwdDefprob = GetStringFromXMLNode(theNode, "String1");
		DiscDefprob = GetStringFromXMLNode(theNode, "String2");
		}

		theNode->Release();
		theNode = NULL;

		resultList->Release();
		resultList = NULL;
	}

	if ( XMLDoc->selectNodes(_bstr_t("//ASSET"), &resultList) == S_OK )
	{
		resultList->get_length(&nbNodes);

		if ( nbNodes == 0 )
		{
			hr = S_FALSE;

			if (resultList)
			{
			   resultList->Release();
			
			   resultList = NULL;
			}
		}

		hr = resultList->get_item(0, &theNode);

		if ( hr == S_OK && theNode != NULL )
		{
		Currency = GetStringFromXMLNode(theNode, "INTEREST_Ccy");
		}

		//Recuperation du basis daycount
		temp = GetStringFromXMLNode(theNode, "INTEREST_Basis");
		daycount = FromSummitDaycountToARMDaycount(temp.c_str());
		int size = CFMatrix->GetCol(0)->GetSize();
		VBasis = new ARM_Vector(size,(double)daycount);
		CFMatrix->Push(VBasis,"Basis");
		theNode->Release();
		theNode = NULL;

		resultList->Release();
		resultList = NULL;
	}

	

	if ( XMLDoc->selectNodes(_bstr_t("//cBLOB_I"), &resultList) == S_OK )
	{
		resultList->get_length(&nbNodes);

		if ( nbNodes == 0 )
		{
			hr = S_FALSE;

			if (resultList)
			{
			   resultList->Release();
			
			   resultList = NULL;
			}
		}

		hr = resultList->get_item(0, &theNode);

		if ( hr == S_OK && theNode != NULL )
		{
		temp = GetStringFromXMLNode(theNode, "Date2");
		if (temp!="") BeginCMDate = GetDateFromXMLNode (theNode,"Date2");
		temp = GetStringFromXMLNode(theNode, "Date1");
		if (temp!="") EndCMDate = GetDateFromXMLNode (theNode,"Date1");
		FixedStartRate = GetDoubleFromXMLNode (theNode, "Rate4");
		OnDefAccrued = (int) GetDoubleFromXMLNode(theNode, "Num3");
		}

		theNode->Release();
		theNode = NULL;

		resultList->Release();
		resultList = NULL;
	}

	ZcCurve = ARMLOCAL_XML_ZCPY_with_summit_stripper(chaineXML,Currency.c_str());
	pFwdDefprob =  ARMLOCAL_XML_DEFPROB_with_summit_stripper(chaineXML,
															 FwdDefprob.c_str(),
															 ZcCurve);

	pDiscDefprob =  ARMLOCAL_XML_DEFPROB_with_summit_stripper(chaineXML,
															 DiscDefprob.c_str(),
															 ZcCurve);


	VolFwdDefprob = ARMLOCAL_XML_VOLCURVE_DEFPROB(chaineXML,FwdDefprob);
	index = ARMLOCAL_XML_CREDIT_INDEX(chaineXML,bookName,custId,dealId);
	

	#ifdef _DEBUG
	FILE *stream2;
	stream2 = fopen("c:\\temp\\Vol.txt", "w+");
	VolFwdDefprob->View("",stream2);
	fclose(stream2);
	#endif _DEBUG

	vector<double> FixDates;
	vector<double> Fixings;
	ARMLOCAL_XML_EVENT_FIXING(chaineXML,"FRC",FixDates,Fixings);

	if (FixDates.size()>0)
	{
		ARM_Vector* FixingDate = CFMatrix->GetColVect("FixingDate");
		ARM_Vector* VFixings = new ARM_Vector(FixingDate->size(),0.);
		for (int il1=0; il1<FixingDate->GetSize();il1++)
			for (int il2=0; il2<FixDates.size();il2++)
			{
				if CHECK_EQUAL(FixingDate->Elt(il1),FixDates[il2])
					VFixings->Elt(il1) = Fixings[il2];

			}

		CFMatrix->Push(VFixings,"FixingRate");
	}

		PricingCMCDSWithCFMatrix(CFMatrix,
							 ZcCurve,
							 pFwdDefprob,
							 pDiscDefprob,
							 VolFwdDefprob,
							 index,
							 BeginCMDate,
							 EndCMDate,
							 FixedStartRate,
							 OnDefAccrued);
							 


	if (VolFwdDefprob)
		delete VolFwdDefprob;
	VolFwdDefprob = NULL;

	if (pFwdDefprob)
		delete pFwdDefprob;
	pFwdDefprob = NULL;

	if (pDiscDefprob)
		delete pDiscDefprob;
	pDiscDefprob = NULL;

	if (Matrix)
		delete Matrix;
	Matrix = NULL;

	if (ZcCurve)
		delete ZcCurve;
	ZcCurve = NULL;

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
		
    if (XMLDoc)
	{
		XMLDoc->Release();
		XMLDoc = NULL;
	}
	if (xmlWCharText) delete xmlWCharText ; xmlWCharText=0 ;

	return (CFMatrix);

	}
	catch(Exception& )
	{
		if (theNode) theNode->Release();
		if (listItem) listItem->Release();
		if (resultList) resultList->Release();
		if (XMLDoc) XMLDoc->Release();

		if (xmlWCharText)
			delete xmlWCharText;
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
			delete xmlWCharText;
		xmlWCharText = NULL;

		hr = S_FALSE;
		throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
         "Error in XML parsing for getting Tranche");
	}

	return NULL;
}

  **/ 
// ------------------------------------
// Generic Smile Correlation Matrix Using XML
// ------------------------------------
extern ARM_Vector* ARMLOCAL_XML_VOLCURVE_DEFPROB_BY_TENOR(const char* chaineXML,
											  const string& DefProbName,
											  const string& TENOR,
											  string RatesOrDate /*R ou D*/,
											  ARM_Date AsOfDate)
{
	VARIANT_BOOL bOK;
	HRESULT hr;

	MSXML2::IXMLDOMDocument *XMLDoc = NULL;
	MSXML2::IXMLDOMNodeList * resultList = NULL;
	MSXML2::IXMLDOMNode * listItem = NULL, * theNode = NULL;
	long nbNodes = 0;
	CCString strtmp;

	ARM_Vector* Rates = NULL;

	wchar_t * xmlWCharText = NULL;

	try
	{
		hr = CoInitialize(NULL); 
//		SUCCEEDED(hr) ? 0 : throw hr;
		hr = CoCreateInstance(__uuidof(MSXML2::DOMDocument30), NULL, CLSCTX_INPROC_SERVER, __uuidof(MSXML2::IXMLDOMDocument), (void**)&XMLDoc);
		SUCCEEDED(hr) ? 0 : throw hr;

		xmlWCharText = constchar2wchar(chaineXML);

		Local_XML_TRACE("c:\\temp\\log_arm_credit.txt","Lecture du XML du smile de correlation");
		XMLDoc->loadXML((_bstr_t)xmlWCharText, &bOK);
	}
	catch(...)
	{

		if (xmlWCharText)
			delete xmlWCharText;
		xmlWCharText = NULL;

		if (XMLDoc) XMLDoc->Release();
		hr = S_FALSE;
		throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"Pb in creating XML document for smile correlation");
	}
	try
	{
	strtmp = (CCString)"//CURVEHEAD/ENTITYLIST/CURVE[../../CVCHAR/MCType = \"ISSUER\" and ../../CVCHAR/Issuer = \"" + 
				(CCString)DefProbName.c_str() + 
				"\" and ../../CVCHAR/Term = \"" + 
				(CCString)TENOR.c_str() + 
				(CCString) "\"]";

	if (XMLDoc->selectNodes((_bstr_t)(const char*)strtmp, &resultList) == S_OK)
	{
			resultList->get_length(&nbNodes);

			if (nbNodes == 0)
			{
				hr = S_FALSE;

				CCString msg ((CCString)"Invalid XML string for VOLCURVE \n" + chaineXML);
		
				throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
								 (const char*) msg);
			}

			Rates = new ARM_Vector(nbNodes,0.);

			for (long indexNode=0 ; indexNode<nbNodes ; indexNode++)
			{
				hr=resultList->get_item(indexNode, &listItem);
				if (hr==S_OK && listItem!=NULL)
				{
					if (RatesOrDate == "R")
						Rates->InitElt(indexNode,100.*GetDoubleFromXMLNode(listItem, "Rate"));
					else
						Rates->InitElt(indexNode,(GetDateFromXMLNode(listItem, "Date")-AsOfDate)/365.);
				}

				if (listItem)
				{
				listItem->Release();
				listItem = NULL;
				}

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
		
    if (XMLDoc)
	{
		XMLDoc->Release();
		XMLDoc = NULL;
	}

	if (xmlWCharText) delete xmlWCharText; xmlWCharText=0; 
	return (Rates);

	}
	catch(Exception& )
	{
		if (theNode) theNode->Release();
		if (listItem) listItem->Release();
		if (resultList) resultList->Release();
		if (XMLDoc) XMLDoc->Release();

		hr = S_FALSE;

		if (xmlWCharText)
			delete xmlWCharText;
		xmlWCharText = NULL;

		throw; 
	}
	catch(...)
	{		
		if (theNode) theNode->Release();
		if (listItem) listItem->Release();
		if (resultList) resultList->Release();
		if (XMLDoc) XMLDoc->Release();

		if (xmlWCharText)
			delete xmlWCharText;
		xmlWCharText = NULL;

		hr = S_FALSE;
		throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
         "Error in XML parsing for getting Tranche");
	}


	return NULL;
}


// ------------------------------------
// Generic Smile Correlation Matrix Using XML
// ------------------------------------
extern void ARMLOCAL_XML_VOLCURVE_DEFPROB_MATURITIES(const char* chaineXML,
												const string& DefProbName,
											    vector<string>& Maturities,
												ARM_Vector*& YF_sort,
												ARM_Date& DAsOf)
{
	VARIANT_BOOL bOK;
	HRESULT hr;

	MSXML2::IXMLDOMDocument *XMLDoc = NULL;
	MSXML2::IXMLDOMNodeList * resultList = NULL;
	MSXML2::IXMLDOMNode * listItem = NULL, * theNode = NULL;
	long nbNodes = 0;
	CCString strtmp;

	ARM_Date AsOfDate;

	wchar_t * xmlWCharText = NULL;
	vector<string> Maturities1;
	ARM_Vector* YF = NULL;

	try
	{
		hr = CoInitialize(NULL); 
//		SUCCEEDED(hr) ? 0 : throw hr;
		hr = CoCreateInstance(__uuidof(MSXML2::DOMDocument30), NULL, CLSCTX_INPROC_SERVER, __uuidof(MSXML2::IXMLDOMDocument), (void**)&XMLDoc);
		SUCCEEDED(hr) ? 0 : throw hr;

		xmlWCharText = constchar2wchar(chaineXML);

		Local_XML_TRACE("c:\\temp\\log_arm_credit.txt","Lecture du XML du smile de correlation");
		XMLDoc->loadXML((_bstr_t)xmlWCharText, &bOK);
	}
	catch(...)
	{

		if (xmlWCharText)
			delete xmlWCharText;
		xmlWCharText = NULL;

		if (XMLDoc) XMLDoc->Release();
		hr = S_FALSE;
		throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"Pb in creating XML document for smile correlation");
	}
	try
	{
	strtmp = (CCString)"//CURVEHEAD/CVCHAR[MCType = \"ISSUER\" and Issuer = \"" + 
				(CCString)DefProbName.c_str() + 
				"\"]";

	if (XMLDoc->selectNodes((_bstr_t)(const char*)strtmp, &resultList) == S_OK)
	{
			resultList->get_length(&nbNodes);

			if (nbNodes == 0)
			{
				hr = S_FALSE;

				CCString msg ((CCString)"Invalid XML string for VOLCURVE \n" + chaineXML);
		
				throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
								 (const char*) msg);
			}

			Maturities.resize(nbNodes);
			Maturities1.resize(nbNodes);
			YF = new ARM_Vector(nbNodes,0.);
			ARM_Date Tmp;

			for (long indexNode=0 ; indexNode<nbNodes ; indexNode++)
			{
				hr=resultList->get_item(indexNode, &listItem);
				if (hr==S_OK && listItem!=NULL)
				{
					Maturities1[indexNode] = GetStringFromXMLNode(listItem, "Term");
					DAsOf = GetDateFromXMLNode(listItem, "AsOfDate");
					Tmp = AddPeriod(DAsOf, Maturities1[indexNode],ARM_DEFAULT_COUNTRY,false,qCredit_Default);
					YF->InitElt(indexNode,(Tmp-DAsOf)/365.);
				}

				if (listItem)
				{
				listItem->Release();
				listItem = NULL;
				}

			}
	}


	ARM_Vector* TMPV = (ARM_Vector*) YF->Clone();
	YF_sort = TMPV->Sort_Compact();
	delete TMPV;

	for (int i=0; i<nbNodes;i++)
		for (int j=0; j<nbNodes;j++)
			if (YF_sort->Elt(i) == YF->Elt(j)) 
				Maturities[i] = Maturities1[j];


	Maturities1.empty();

	if (YF) delete YF;
	YF = NULL;

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
		
    if (XMLDoc)
	{
		XMLDoc->Release();
		XMLDoc = NULL;
	}

	if (xmlWCharText) delete xmlWCharText; xmlWCharText=0; 
	}
	catch(Exception& )
	{
		if (theNode) theNode->Release();
		if (listItem) listItem->Release();
		if (resultList) resultList->Release();
		if (XMLDoc) XMLDoc->Release();

		hr = S_FALSE;

		if (xmlWCharText)
			delete xmlWCharText;
		xmlWCharText = NULL;

		throw; 
	}
	catch(...)
	{		
		if (theNode) theNode->Release();
		if (listItem) listItem->Release();
		if (resultList) resultList->Release();
		if (XMLDoc) XMLDoc->Release();

		if (xmlWCharText)
			delete xmlWCharText;
		xmlWCharText = NULL;

		hr = S_FALSE;
		throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
         "Error in XML parsing for getting Tranche");
	}
}


// ------------------------------------
// Generic Smile Correlation Matrix Using XML
// ------------------------------------
extern ARM_VolCurve* ARMLOCAL_XML_VOLCURVE_DEFPROB(const char* chaineXML,const string& DefProbName)
{

	ARM_VolCurve* volcurve = NULL;
	vector<string> Maturities;
	ARM_Vector*	YearFrac = NULL;
	ARM_Matrix* Matrix = NULL;

	ARM_Date AsOf;

	ARMLOCAL_XML_VOLCURVE_DEFPROB_MATURITIES(chaineXML,DefProbName,Maturities,YearFrac,AsOf);

	ARM_Vector* Dates = ARMLOCAL_XML_VOLCURVE_DEFPROB_BY_TENOR(chaineXML,
															  DefProbName,
															  Maturities[0],
															  "D",
															  AsOf);

	Matrix = new ARM_Matrix(Dates->GetSize(),YearFrac->GetSize());

	for (int i=0; i<Maturities.size();i++)
	{
		ARM_Vector* vector = ARMLOCAL_XML_VOLCURVE_DEFPROB_BY_TENOR(chaineXML,
															  DefProbName,
															  Maturities[i],
															  "R",
															  AsOf);
		for (int j=0;j<Dates->GetSize();j++)
			Matrix->Elt(j,i) = vector->Elt(j);

		if (vector) delete vector;

	}

	volcurve = new ARM_VolLInterpol(AsOf,Dates,YearFrac, Matrix,1, K_SMILE_VOL,ARM_DEFAULT_CURRENCY);	
	volcurve->SetIndexName((char*)DefProbName.c_str());

	return (volcurve);
}
