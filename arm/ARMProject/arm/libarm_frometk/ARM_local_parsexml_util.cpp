//#include <wtypes.h>

/// headers to remove definition of va_start and va_end as this is redefined later on!
/// handle this with care after sorting out why this is so!!
#include "ARM\libarm_local\undef_va_vars.h"
#include "ARM\libarm_local\firstToBeIncluded.h"

#include "ARM\libarm_local\ARM_local_glob.h"
#include "inst\swapleg.h" // sinon conflit avec #import de arm_local_parsexml_util.h
#include "ARM\libarm_frometk\arm_local_parsexml_util.h"


#include "util\fromto.h"
#include "util\refvalue.h"
#include "util\exercise.h"
#include "ccy\currency.h"
#include "inst\cmsleg.h"
#include "ARM_ccy.h"

#include "libCCtools++\CCString.h"

using ARM::std::vector<double>;


char * Summit_DefaultDate="99991231";

ARM_Date GetDateFromXMLNode (MSXML2::IXMLDOMNode* node, const char* nodeName)
{
	MSXML2::IXMLDOMNode* tmpNode = NULL;
	BSTR resultat = NULL;
	ARM_Date resDate(ARM_DEFAULT_DATE);

	try
	{
		node->selectSingleNode(_bstr_t(nodeName), &tmpNode);
		if (tmpNode != NULL)
		{
			tmpNode->get_text(&resultat);
			_bstr_t ff(resultat,false);
			char * ff1=(char *)ff;

			if ((strcmp(ff1,"") != 0) && (strcmp(ff1,"SUN") != 0))
				resDate = ARM_Date(ff1,"YYYYMMDD");

			if (resultat) SysFreeString(resultat);
			tmpNode->Release();
			tmpNode = NULL;
		}
		else
		{
			CCString msg((CCString)"Pb in getting Date from node " + nodeName);
	
			throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
					 (char*) msg);
		}
	}
	catch(...)
	{
		CCString msg((CCString)"Pb in getting Date from node " + nodeName);

		throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
				 (char*) msg);
	}

	return resDate;
}

/**void populateDate(ARM_Date&res,const MSXML2::IXMLDOMNode&parentNode,const std::string&nodeName)
{
	MSXML2::IXMLDOMNode* tmpNode(0); 
	parentNode.selectSingleNode(_bstr_t(nodeName),&tmpNode); 
	if(!tmpNode) ICMTHROW(ERR_INVALID_ARGUMENT,"populateDate: can't find "<<nodeName); 

}**/


double GetDoubleFromXMLNode (MSXML2::IXMLDOMNode* node, const char* nodeName, bool throwExcept)
{
	MSXML2::IXMLDOMNode* tmpNode = NULL;
	BSTR resultat = NULL;
	double resDouble = 0.0;

	try
	{
		node->selectSingleNode(_bstr_t(nodeName), &tmpNode);
		if (tmpNode != NULL)
		{
			tmpNode->get_text(&resultat);
			_bstr_t ff(resultat,false);
			char * ff1=(char *)ff;

			if ( (strcmp(ff1,"") == 0) && (throwExcept) )
				resDouble = ARM_MISSING_VALUE;
			else
				resDouble = atof((const char*) ff1);

			if (resultat) SysFreeString(resultat);
			tmpNode->Release();
			tmpNode = NULL;
		}
		else
		{
			if (throwExcept)
				return ARM_MISSING_VALUE;
		}
	}
	catch(...)
	{
	}

	return resDouble;
}

int GetIntFromXMLNode (MSXML2::IXMLDOMNode* node, const char* nodeName)
{
	MSXML2::IXMLDOMNode* tmpNode = NULL;
	BSTR resultat = NULL;
	int resInteger = 0;

	try
	{
		node->selectSingleNode(_bstr_t(nodeName), &tmpNode);
		if (tmpNode != NULL)
		{
			tmpNode->get_text(&resultat);
			_bstr_t ff(resultat,false);
			char * ff1=(char *)ff;

			resInteger = atoi((const char*) ff1);

			if (resultat) SysFreeString(resultat);
			tmpNode->Release();
			tmpNode = NULL;
		}
	}
	catch(...)
	{
	}

	return resInteger;
}


CCString GetCCStringFromXMLNode (MSXML2::IXMLDOMNode* node, const char* nodeName)
{
	MSXML2::IXMLDOMNode* tmpNode = NULL;
	BSTR resultat = NULL;
	CCString value;

	try
	{
		node->selectSingleNode(_bstr_t(nodeName), &tmpNode);
		if (tmpNode != NULL)
		{
			tmpNode->get_text(&resultat);
			_bstr_t ff(resultat,false);
			char * ff1=(char *)ff;

			value = ff1;

			if (resultat) SysFreeString(resultat);
			tmpNode->Release();
			tmpNode = NULL;
		}
		else
		{
			value = "";
		}
	}
	catch(...)
	{
	}
	return value;
}


CCString GetCustumerId(MSXML2::IXMLDOMNode* node)
{
	return GetCCStringFromXMLNode(node,"Cust");
}

CCString GetDealId(MSXML2::IXMLDOMNode* node)
{
	return GetCCStringFromXMLNode(node,"DealId");
}

CCString GetBook(MSXML2::IXMLDOMNode* node)
{
	return GetCCStringFromXMLNode(node,"Book");
}


//***************************************************************************************************
// Appel au Niveau OPTION
//***************************************************************************************************

int GetPorC(MSXML2::IXMLDOMNode* node)
{
	CCString res = GetCCStringFromXMLNode(node,"PorC");

	if (strcmp((const char*)res,"P")==0)
			// swaption Payeuse du taux fixe= CALL		
		return K_PUT;
	else
		return K_CALL;
}

ARM_Currency* GetOptionCcy(MSXML2::IXMLDOMNode* node)
{
	MSXML2::IXMLDOMNode* theNode=NULL;
	BSTR resultat = NULL;
	ARM_Currency* ccy=NULL;

	node->selectSingleNode(_bstr_t("PREMDATA_Ccy"), &theNode);
	if (theNode!=NULL)
	{
		theNode->get_text(&resultat);
		_bstr_t ff(resultat,false);
		char * ff1=(char *)ff;
		ccy = new ARM_Currency(ff1);

		theNode->Release();
		theNode=NULL;
		if (resultat) SysFreeString(resultat);
	}	
	return ccy;
}



int GetBLOBTiming(MSXML2::IXMLDOMNode* node, int n)
{
	MSXML2::IXMLDOMNode* theNode=NULL;
	BSTR resultat = NULL;
	int value=K_ARREARS;

	char nodeName[30];
	sprintf(nodeName, "ProdData/cBLOB_I/Timing%i", n);
	
	node->selectSingleNode(_bstr_t(nodeName), &theNode);
	if (theNode!=NULL)
	{
		theNode->get_text(&resultat);
		_bstr_t ff(resultat,false);
		char * ff1=(char *)ff;
		value = ARM_ConvPayResetRule(CCString(ff1));
		
		theNode->Release();
		theNode=NULL;
		if (resultat) SysFreeString(resultat);
	}
	else
	{
		CCString msg((CCString)"Pb in getting Timing from node " + nodeName);

		throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT, (char*) msg);
	}

	return value;
}

ARM_ReferenceValue*  GetCRASpread(MSXML2::IXMLDOMNode* node)
{
	MSXML2::IXMLDOMNode* theNode=NULL;
	MSXML2::IXMLDOMNode* theNode2=NULL;
	MSXML2::IXMLDOMNodeList* resultList = NULL;
	ARM_ReferenceValue* newRefVal = NULL;
	BSTR resultat = NULL;
	HRESULT hr;
	long nb=0;
	
	VECTOR<double> RefDates;
	VECTOR<double> Spreads;

	node->selectNodes(_bstr_t("ProdData/cBLOB_I/Schedule1/cBLBLST_I"), &resultList);
	if (resultList!=NULL)
	{
		long nb;
		resultList->get_length(&nb);
		if (nb==0)
		{
			hr = S_FALSE;

			CCString msg((CCString)"Invalid XML string for getting CRA Spread \n");

			throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
					 (char*) msg);
		}
		for (long i=0 ; i<nb ; i++)
		{	
			hr=resultList->get_item(i, &theNode);
			if (hr==S_OK && theNode!=NULL)
			{
				// voir si le refvalue se cale sur les StartDate ou sur les EndDate...
				RefDates.push_back(GetDateFromXMLNode(theNode,"StartDate").GetJulian());
				
				theNode->selectSingleNode(_bstr_t("Amount1"), &theNode2);
				if (theNode2!=NULL)
				{
					theNode2->get_text(&resultat);
					_bstr_t ff(resultat,false);
					char * ff1=(char *)ff;
					Spreads.push_back(atof(ff1)/100.0);

					if (resultat) SysFreeString(resultat);
					theNode2->Release();
					theNode2=NULL;
				}
			}
			if (theNode) theNode->Release();
			theNode = NULL;
		}
		if (resultList) resultList->Release();
		resultList=NULL;
	}

	if (nb==1)
		newRefVal = new ARM_ReferenceValue(Spreads[0]);
	else
	{
		ARM_Vector* VRefDates= CreateARMVectorFromVECTOR(RefDates);
		ARM_Vector* VSpreads= CreateARMVectorFromVECTOR(Spreads);
		newRefVal = new ARM_ReferenceValue(VRefDates, VSpreads,K_YIELD,1);
		newRefVal->SetCalcMethod(K_STEPUP_LEFT);
	}
	return newRefVal;
}

ARM_ReferenceValue*  GetCRABarrierDown(MSXML2::IXMLDOMNode* node)
{	
	MSXML2::IXMLDOMNode* theNode=NULL;
	MSXML2::IXMLDOMNode* theNode2=NULL;
	MSXML2::IXMLDOMNodeList* resultList = NULL;
	ARM_ReferenceValue* newRefVal = NULL;
	BSTR resultat = NULL;
	HRESULT hr;
	long nb;

	VECTOR<double> RefDates;
	VECTOR<double> Barrieres;
	node->selectNodes(_bstr_t("ProdData/cBLOB_I/Schedule1/cBLBLST_I"), &resultList);
	if (resultList!=NULL)
	{
		resultList->get_length(&nb);
		if (nb==0)
		{
			hr = S_FALSE;

			CCString msg((CCString)"Invalid XML string for getting CRA Barrier Down \n");

			throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
					 (char*) msg);
		}
		for (long i=0 ; i<nb ; i++)
		{	
			hr=resultList->get_item(i, &theNode);
			if (hr==S_OK && theNode!=NULL)
			{
				// voir si le refvalue se cale sur les StartDate ou sur les EndDate...
				RefDates.push_back(GetDateFromXMLNode(theNode,"StartDate").GetJulian());
				
				theNode->selectSingleNode(_bstr_t("Rate1"), &theNode2);
				if (theNode2!=NULL)
				{
					theNode2->get_text(&resultat);
					_bstr_t ff(resultat,false);
					char * ff1=(char *)ff;
					Barrieres.push_back(atof(ff1)*100);

					if (resultat) SysFreeString(resultat);
					theNode2->Release();
					theNode2=NULL;
				}
			}
			if (theNode) theNode->Release();
			theNode = NULL;
		}
	}

	if (resultList) resultList->Release();
	resultList=NULL;
	
	if (nb==1)
		newRefVal = new ARM_ReferenceValue(Barrieres[0]);
	else
	{
		ARM_Vector* VRefDates= CreateARMVectorFromVECTOR(RefDates);
		ARM_Vector* VBarrieres= CreateARMVectorFromVECTOR(Barrieres);
		newRefVal = new ARM_ReferenceValue(VRefDates, VBarrieres,K_YIELD,1);
		newRefVal->SetCalcMethod(K_STEPUP_LEFT);
	}

	return newRefVal;
}

ARM_ReferenceValue*  GetCRABarrierUp(MSXML2::IXMLDOMNode* node)
{
	MSXML2::IXMLDOMNode* theNode=NULL;
	MSXML2::IXMLDOMNode* theNode2=NULL;
	MSXML2::IXMLDOMNodeList* resultList = NULL;
	ARM_ReferenceValue* newRefVal = NULL;
	BSTR resultat = NULL;
	HRESULT hr;
	long nb;

	VECTOR<double> RefDates;
	VECTOR<double> Barrieres;
	node->selectNodes(_bstr_t("ProdData/cBLOB_I/Schedule1/cBLBLST_I"), &resultList);
	if (resultList!=NULL)
	{
		resultList->get_length(&nb);
		if (nb==0)
		{
			hr = S_FALSE;

			CCString msg((CCString)"Invalid XML string for getting CRA Barrier Up \n");

			throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
					 (char*) msg);
		}
		for (long i=0 ; i<nb ; i++)
		{	
			hr=resultList->get_item(i, &theNode);
			if (hr==S_OK && theNode!=NULL)
			{
				// voir si le refvalue se cale sur les StartDate ou sur les EndDate...
				RefDates.push_back(GetDateFromXMLNode(theNode,"StartDate").GetJulian());

				theNode->selectSingleNode(_bstr_t("Rate2"), &theNode2);
				if (theNode2!=NULL)
				{
					theNode2->get_text(&resultat);
					_bstr_t ff(resultat,false);
					char * ff1=(char *)ff;
					Barrieres.push_back(atof(ff1)*100);

					if (resultat) SysFreeString(resultat);
					theNode2->Release();
					theNode2=NULL;
				}
			}
			if (theNode) theNode->Release();
			theNode = NULL;
		}
	}
	if (resultList) resultList->Release();
	resultList=NULL;

	if (nb==1)
		newRefVal = new ARM_ReferenceValue(Barrieres[0]);
	else
	{
		ARM_Vector* VRefDates= CreateARMVectorFromVECTOR(RefDates);
		ARM_Vector* VBarrieres= CreateARMVectorFromVECTOR(Barrieres);
		newRefVal = new ARM_ReferenceValue(VRefDates, VBarrieres,K_YIELD,1);
		newRefVal->SetCalcMethod(K_STEPUP_LEFT);
	}

	return newRefVal;
}

void ConvSummitIndexToARMIndexTerm(char* summitIndex, ARM_Currency* ccy, char* armIndexTerm)
{
	char IndexTerm[5];

	// index du type [6mEUReurib]
	// on recupere la maturité IndexTerm
// FIXMEFRED: mig.vc8 (30/05/2007 17:52:16):const missing
	const char* tmpStr = strstr((const char*)summitIndex, (const char*)ccy->GetCcyName());

	if (tmpStr == NULL)
	{
		CCString msg((CCString)"Currency not found in summitIndex " + summitIndex);

		throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
						 (char*) msg);
	}

	int pDest = (int)(tmpStr - summitIndex + 1);

	strncpy(IndexTerm,summitIndex,pDest);

	// decalage a cause du crochet
	for (int i = 0; i < pDest - 2; i++)
		IndexTerm[i] = IndexTerm[i+1];

	IndexTerm[pDest - 2] = '\0';

	strcpy(armIndexTerm,strupr(IndexTerm));
}

void ConvSummitIndexToARMIndex(char* summitIndex, ARM_Currency* ccy, char* armIndex, int isCMT = 0)
{
	int nbMatu;
	char frq;
	char buffer, buffer2;
	char sccy[20];
	char indexName[20];

	summitIndex = strupr(summitIndex);

	// index du type [6mEUReurib]
	sscanf((const char*)summitIndex,"%c%d%c%3s%[0-9,A-Z]%c",&buffer,&nbMatu,&frq,sccy,indexName,&buffer2);

	if (strcmp(sccy, ccy->GetCcyName()) != 0)
	{
		// QUANTO case !!!
		ccy->SetCcyName(sccy);
		ccy->SetCurrencyProperties();
		/*ARM_Currency* ccytmp = new ARM_Currency(sccy);
		ccy = ccytmp;
		delete ccytmp;
		ccytmp = NULL;*/
	}

	if ( (strcmp(indexName,"EUR3M") == 0)
		|| (strcmp(indexName,"EUR12") == 0) )
	{
		strcpy(armIndex,indexName);
	}
	else
	{
		if (strcmp(indexName,"EURIB") == 0)
			strcat(indexName, "OR");

		if (strcmp(indexName,"STIBO") == 0)
			strcpy(indexName, "LIBOR");

		// Cas d'un EURIBOR10Y a convertir en CMS10...
		char matu[3];
		sprintf(matu,"%d",nbMatu);

		if (frq == 'Y')
		{
			if (isCMT)
				strcpy(indexName,"CMT");
			else
				strcpy(indexName,"CMS");
			strcat(indexName, matu);
			strcpy(armIndex,indexName);
		}
		else
		{
			strcat(indexName, matu);
			strcat(indexName, "M");
			strcpy(armIndex,indexName);
		}
	}
}


int GetBLOBIdx(MSXML2::IXMLDOMNode* node, ARM_Currency* ccy,int n)
{	
	MSXML2::IXMLDOMNode* theNode=NULL;
	BSTR resultat = NULL;
	char summitIndex[15];
	char armIndex[15];

	char nodeName[30];
	sprintf(nodeName, "ProdData/cBLOB_I/Index%i", n);

	node->selectSingleNode(_bstr_t(nodeName), &theNode);
	if (theNode!=NULL)
	{
		BSTR resultat = NULL;
		theNode->get_text(&resultat);

		_bstr_t ff(resultat,false);
		char * ff1=(char *)ff;
		strcpy(summitIndex,ff1);
		
		theNode->Release();
		theNode=NULL;
	}
	if (resultat) SysFreeString(resultat);
		
	if ( (strcmp(summitIndex,"")==0) || (summitIndex[0] != '[') )
		return K_FIXED;
	else
	{
		char tmp[20];
		ConvSummitIndexToARMIndex(summitIndex,ccy,tmp);
		strcpy(armIndex,tmp);

		return ARM_ConvIrType(armIndex);
	}
}



ARM_Date GetOptionExpiry(MSXML2::IXMLDOMNode* node)
{
	return GetDateFromXMLNode(node,"ExpDate");
}

ARM_Date GetOptionStartExpiry(MSXML2::IXMLDOMNode* node)
{
	return GetDateFromXMLNode(node,"FstExpDate");
}

ARM_Vector* GetBermudDates(MSXML2::IXMLDOMNode* node)
{
	MSXML2::IXMLDOMNode* theNode=NULL;
	MSXML2::IXMLDOMNode* theNode1=NULL;
	MSXML2::IXMLDOMNode* theNode2=NULL;
	MSXML2::IXMLDOMNodeList* resultList = NULL;
	BSTR resultat = NULL;
	HRESULT hr;
	ARM_Vector* res = NULL;

	VECTOR<double> RefDates;
	node->selectNodes(_bstr_t("OpEvents/OPEVENT"), &resultList);
	if (resultList!=NULL)
	{
		long nb;
		resultList->get_length(&nb);
		if (nb==0)
		{
			hr = S_FALSE;

			CCString msg((CCString)"Invalid XML string for getting BermudDates\n");

			throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
					 (char*) msg);
		}
		for (long i=0 ; i<nb ; i++)
		{	
			hr=resultList->get_item(i, &theNode);
			if (hr==S_OK && theNode!=NULL)
			{
				// voir si le refvalue se cale sur les StartDate ou sur les EndDate...
				theNode->selectSingleNode(_bstr_t("ADate"), &theNode1);
				theNode->selectSingleNode(_bstr_t("Type"), &theNode2);
				if (theNode2!=NULL)
				{
					theNode2->get_text(&resultat);
					_bstr_t ff(resultat,false);
					char * ff1=(char *)ff;
					if ((theNode1!=NULL) && (strcmp(ff1,"CL")==0))
					{
						theNode1->get_text(&resultat);
						_bstr_t ff(resultat,false);
						char * ff1=(char *)ff;

						ARM_Date tmpDate (ff1,"YYYYMMDD");

						RefDates.push_back(tmpDate.GetJulian());
					}
					if (resultat) SysFreeString(resultat);
					theNode2->Release();
					theNode2=NULL;
				}
				if (theNode1) theNode1->Release();
				theNode1=NULL;

				if (theNode) theNode->Release();
				theNode = NULL;
			}
		}
		if (resultList) resultList->Release();
		resultList=NULL;
	}

	res = CreateARMVectorFromVECTOR(RefDates);
	return res;
}


int GetExerciceStyleType(MSXML2::IXMLDOMNode* node)
{
	MSXML2::IXMLDOMNode* theNode=NULL;
	BSTR resultat = NULL;
	int value= K_EUROPEAN;
	
	node->selectSingleNode(_bstr_t("Style"), &theNode);
	if (theNode!=NULL)
	{
		theNode->get_text(&resultat);

		_bstr_t ff(resultat,false);
		char * ff1=(char *)ff;
		
		if (strcmp(ff1,"EURO") ==0)
			value = K_EUROPEAN;
		else if (strcmp(ff1,"AMER") ==0)
			value = K_AMERICAN;
		else if (strcmp(ff1,"BERM") ==0)
			value = K_BERMUDAN;
		
		theNode->Release();
		theNode=NULL;
		if (resultat) SysFreeString(resultat);
	}
	return value;
}


ARM_Vector* GetExerciseDates(MSXML2::IXMLDOMNode* node)
{
	MSXML2::IXMLDOMNode* theNode=NULL;
	MSXML2::IXMLDOMNode* theNode1=NULL;
	MSXML2::IXMLDOMNode* theNode2=NULL;
	MSXML2::IXMLDOMNodeList* resultList = NULL;
	BSTR resultat = NULL;
	HRESULT hr;

	VECTOR<double> RefDates;
	node->selectNodes(_bstr_t("OpEvents/OPEVENT"), &resultList);
	if (resultList!=NULL)
	{
		long nb;
		resultList->get_length(&nb);
		if (nb==0)
		{
			hr = S_FALSE;

			CCString msg((CCString)"Invalid XML string for getting BermudDates\n");

			throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
					 (char*) msg);
		}
		for (long i=0 ; i<nb ; i++)
		{	
			hr=resultList->get_item(i, &theNode);
			if (hr==S_OK && theNode!=NULL)
			{
				// voir si le refvalue se cale sur les StartDate ou sur les EndDate...
				theNode->selectSingleNode(_bstr_t("Date"), &theNode1);
				theNode->selectSingleNode(_bstr_t("Type"), &theNode2);
				if (theNode2!=NULL)
				{
					theNode2->get_text(&resultat);
					_bstr_t ff(resultat,false);
					char * ff1=(char *)ff;
					if ((theNode1!=NULL) && (strcmp(ff1,"CL")==0))
					{
						theNode1->get_text(&resultat);
						_bstr_t ff(resultat,false);
						char * ff1=(char *)ff;

						ARM_Date tmpDate (ff1,"YYYYMMDD");

						RefDates.push_back(tmpDate.GetJulian());

						theNode1->Release();
						theNode1=NULL;
					}
					if (resultat) SysFreeString(resultat);
					theNode2->Release();
					theNode2=NULL;
				}
			}
			if (theNode) theNode->Release();
			theNode = NULL;
		}
		if (resultList) resultList->Release();
		resultList=NULL;
	}
	return CreateARMVectorFromVECTOR(RefDates);
}


ARM_ExerciseStyle * GetExerciceStyle(MSXML2::IXMLDOMNode* node)
{
	ARM_ExerciseStyle * style = NULL;
	int styletype = GetExerciceStyleType(node);
		
	if (styletype == K_EUROPEAN)
	{
		// Constructor for European Style
		ARM_Date xDate = GetOptionExpiry(node);
		style = new ARM_ExerciseStyle(xDate);  
	}
	else if (styletype == K_AMERICAN)
	{
		// Constructor for American Style
		ARM_Date xStartDate = GetOptionStartExpiry(node);
		ARM_Date xEndDate = GetOptionExpiry(node);
		style = new ARM_ExerciseStyle(xStartDate, xEndDate);  
	}
	else if (styletype == K_BERMUDAN)
	{
		// Constructor for Bermudan Style
		ARM_Vector* NoticeDates = GetBermudDates(node);
		ARM_Vector* ExpiryDates = GetExerciseDates(node);
		// Store notice and corresponding expiry dates
		style = new ARM_ExerciseStyle(NoticeDates, ExpiryDates);
		delete NoticeDates;
		delete ExpiryDates;
	}
	
	return style;

}

ARM_Vector* GetFees(MSXML2::IXMLDOMNode* node)
{
	MSXML2::IXMLDOMNode* theNode=NULL;
	MSXML2::IXMLDOMNode* theNode1=NULL;
	MSXML2::IXMLDOMNode* theNode2=NULL;
	MSXML2::IXMLDOMNodeList* resultList = NULL;
	BSTR resultat = NULL;
	HRESULT hr;
	ARM_Vector* res = NULL;

	VECTOR<double> RefFees;

	node->selectNodes(_bstr_t("OpEvents/OPEVENT"), &resultList);
	if (resultList!=NULL)
	{
		long nb;
		resultList->get_length(&nb);
		if (nb==0)
		{
			hr = S_FALSE;

			CCString msg((CCString)"Invalid XML string for getting Fees\n");

			throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
					 (char*) msg);
		}
		for (long i=0 ; i<nb ; i++)
		{	
			hr=resultList->get_item(i, &theNode);
			if (hr==S_OK && theNode!=NULL)
			{
				// voir si le refvalue se cale sur les StartDate ou sur les EndDate...
				theNode->selectSingleNode(_bstr_t("Amount"), &theNode1);
				theNode->selectSingleNode(_bstr_t("Type"), &theNode2);
				if (theNode2!=NULL)
				{
					theNode2->get_text(&resultat);
					_bstr_t ff(resultat,false);
					char * ff1=(char *)ff;
					if ((theNode1!=NULL) && (strcmp(ff1,"CL")==0))
					{
						theNode1->get_text(&resultat);
						_bstr_t ff(resultat,false);
						char * ff1=(char *)ff;
						double fee= atof((const char*) ff1);
						RefFees.push_back(fee);
					}
					if (resultat) SysFreeString(resultat);
					theNode2->Release();
					theNode2=NULL;
				}
				if (theNode1) theNode1->Release();
				theNode1=NULL;
			}
			if (theNode) theNode->Release();
			theNode = NULL;
		}
		if (resultList) resultList->Release();
		resultList=NULL;
	}

	res = CreateARMVectorFromVECTOR(RefFees);
	return res;
}

double GetStrike(MSXML2::IXMLDOMNode* node)
{	
	MSXML2::IXMLDOMNode* theNode=NULL;
	BSTR resultat = NULL;
	int value=0.0;

	node->selectSingleNode(_bstr_t("Strike"), &theNode);
	if (theNode!=NULL)
	{
		theNode->get_text(&resultat);
		_bstr_t ff(resultat,false);
		char * ff1=(char *)ff;
		
		value = atof((const char*) ff1)*100.0;
		
		theNode->Release();
		theNode=NULL;
		if (resultat) SysFreeString(resultat);
	}
	
	return value;
}


ARM_ReferenceValue*  GetStrikes(MSXML2::IXMLDOMNode* node)
{
	ARM_ReferenceValue* value=NULL;
	int styletype=GetExerciceStyleType(node);
	
	if ( styletype == K_EUROPEAN )
	{
		double strike = GetStrike(node);
		value = new ARM_ReferenceValue(strike);
	}
	else if ((styletype ==K_AMERICAN) || (styletype ==K_BERMUDAN))
	{
		ARM_Vector* Dates = GetBermudDates(node);
		ARM_Vector* Fees  = GetFees(node);

		if ( Dates->GetSize() == 1 )
		{
		   value = new ARM_ReferenceValue(Fees->Elt(0));
		}
		else
		{
		   value = new ARM_ReferenceValue(Dates,Fees,K_YIELD,1);
			
		   // interpol par STEP
	       value->SetCalcMethod(K_STEPUP_LEFT);
		}
	}
	
	return value;
	
}

ARM_ReferenceValue*  GetExerFees(MSXML2::IXMLDOMNode* node)
{
	ARM_ReferenceValue* value=NULL;

	ARM_Vector* Dates=GetBermudDates(node);
	ARM_Vector* Fees=GetFees(node);
	value = new ARM_ReferenceValue(Dates,Fees);
	return value;
	
}

int GetCRFIndexDayCount(MSXML2::IXMLDOMNode* node)
{
	MSXML2::IXMLDOMNode* theNode=NULL;
	BSTR resultat = NULL;
	int value=KNOBASE;
	
	node->selectSingleNode(_bstr_t("ProdData/cCRFARM_I/Basis"), &theNode);
	if (theNode!=NULL)
	{
		theNode->get_text(&resultat);

		_bstr_t ff(resultat,false);
		char * ff1=(char *)ff;
		
		if (strcmp((const char*) ff1,"") != 0)
			value = FromSummitDaycountToARMDaycount((const char*) ff1);
		
		theNode->Release();
		theNode=NULL;
		if (resultat) SysFreeString(resultat);
	}

	return value;
}

int GetCRFDayCount(MSXML2::IXMLDOMNode* node)
{
	MSXML2::IXMLDOMNode* theNode=NULL;
	BSTR resultat = NULL;
	int value=KNOBASE;
	
	node->selectSingleNode(_bstr_t("ProdData/cCRFARM_I/Basis2"), &theNode);
	if (theNode!=NULL)
	{
		theNode->get_text(&resultat);

		_bstr_t ff(resultat,false);
		char * ff1=(char *)ff;
		
		if (strcmp((const char*) ff1,"") != 0)
			value = FromSummitDaycountToARMDaycount((const char*) ff1);
		
		theNode->Release();
		theNode=NULL;
		if (resultat) SysFreeString(resultat);
	}

	return value;
}


int GetCRFTiming(MSXML2::IXMLDOMNode* node)
{
	MSXML2::IXMLDOMNode* theNode=NULL;
	BSTR resultat = NULL;
	int value=K_ADVANCE;
	
	node->selectSingleNode(_bstr_t("ProdData/cCRFARM_I/Timing"), &theNode);
	if (theNode!=NULL)
	{
		BSTR resultat = NULL;
		theNode->get_text(&resultat);

		_bstr_t ff(resultat,false);
		char * ff1=(char *)ff;
		value = ARM_ConvPayResetRule(CCString((char *)ff));

		theNode->Release();
		theNode=NULL;
		if (resultat) SysFreeString(resultat);
	}
	return value;

}

int GetCRFResetGap(MSXML2::IXMLDOMNode* node)
{
	MSXML2::IXMLDOMNode* theNode=NULL;
	BSTR resultat = NULL;
	int value=0;
	
	node->selectSingleNode(_bstr_t("ProdData/cCRFARM_I/ResetGap"), &theNode);
	if (theNode!=NULL)
	{
		BSTR resultat = NULL;
		theNode->get_text(&resultat);

		_bstr_t ff(resultat,false);
		char * ff1=(char *)ff;

		value = FromSummitGapToARMGap((const char*) ff1);

		theNode->Release();
		theNode=NULL;
		if (resultat) SysFreeString(resultat);
	}
	return value;

}


void GetCRFResetCal(MSXML2::IXMLDOMNode* node, char* resetCal)
{
	MSXML2::IXMLDOMNode* theNode=NULL;
	BSTR resultat = NULL;
	
	node->selectSingleNode(_bstr_t("ProdData/cCRFARM_I/ResetCal"), &theNode);
	if (theNode!=NULL)
	{
		BSTR resultat = NULL;
		theNode->get_text(&resultat);

		_bstr_t ff(resultat,false);
		char * ff1=(char *)ff;

		strcpy(resetCal,ff1);

		theNode->Release();
		theNode=NULL;
		if (resultat) SysFreeString(resultat);
	}
}

void GetCRFIndexTerm(MSXML2::IXMLDOMNode* node, ARM_Currency* ccy, char* indexTerm)
{	
	MSXML2::IXMLDOMNode* theNode=NULL;
	BSTR resultat = NULL;
	char summitIndex[15];

	node->selectSingleNode(_bstr_t("ProdData/cCRFARM_I/Index"), &theNode);
	if (theNode!=NULL)
	{
		BSTR resultat = NULL;
		theNode->get_text(&resultat);

		_bstr_t ff(resultat,false);
		char * ff1=(char *)ff;
		strcpy(summitIndex,ff1);
		
		theNode->Release();
		theNode=NULL;
	}
	if (resultat) SysFreeString(resultat);
		
	if (strcmp(summitIndex,"")==0)
		strcpy(indexTerm,"FIXED");
	else
	{
		char tmp[20];
		ConvSummitIndexToARMIndexTerm(summitIndex,ccy,tmp);
		strcpy(indexTerm,tmp);
	}
}

ARM_Date GetCRFFixEndDate(MSXML2::IXMLDOMNode* node)
{
	MSXML2::IXMLDOMNode* theNode=NULL;
	
	node->selectSingleNode(_bstr_t("ProdData/cCRFARM_I"), &theNode);
	if (theNode!=NULL)
	{
		ARM_Date tmpDate = GetDateFromXMLNode(theNode,"FromDate");

		if (theNode) theNode->Release();
		theNode = NULL;

		return tmpDate;
	}
	else
	{
		CCString msg((CCString)"Pb in getting FixEndDate from node ProdData/cCRFARM_I/FromDate");
		throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT, (char*) msg);
	}
}

ARM_ReferenceValue*  GetCRFOneLeverage(MSXML2::IXMLDOMNode* node)
{
	MSXML2::IXMLDOMNode* theNode=NULL;
	ARM_ReferenceValue* newRefVal = NULL;
	BSTR resultat = NULL;
	
	double Spread;
	node->selectSingleNode(_bstr_t("ProdData/cCRFARM_I/Coef"), &theNode);
	if (theNode!=NULL)
	{
		theNode->get_text(&resultat);
		_bstr_t ff(resultat,false);
		char * ff1=(char *)ff;
		Spread=atof(ff1);

		if (resultat) SysFreeString(resultat);
		theNode->Release();
		theNode=NULL;
	}

	newRefVal = new ARM_ReferenceValue(Spread);
	
	return newRefVal;
}

ARM_ReferenceValue*  GetCRFLeverage(MSXML2::IXMLDOMNode* node)
{
	MSXML2::IXMLDOMNode* theNode=NULL;
	MSXML2::IXMLDOMNode* theNode2=NULL;
	MSXML2::IXMLDOMNodeList* resultList = NULL;
	ARM_ReferenceValue* newRefVal = NULL;
	BSTR resultat = NULL;
	HRESULT hr;
	long nb=0;
	
	VECTOR<double> RefDates;
	VECTOR<double> Spreads;
	node->selectNodes(_bstr_t("ProdData/cCRFARM_I/Schedule/cCRFLST_I"), &resultList);
	if (resultList!=NULL)
	{
		long nb;
		resultList->get_length(&nb);
		if (nb==0)
		{
			newRefVal = GetCRFOneLeverage(node);
			if (resultList) resultList->Release();
			resultList = NULL;
			return newRefVal;
		}
		for (long i = 0 ; i < nb ; i++)
		{	
			hr = resultList->get_item(i, &theNode);
			if (hr == S_OK && theNode != NULL)
			{
				// voir si le refvalue se cale sur les StartDate ou sur les EndDate...
				RefDates.push_back(GetDateFromXMLNode(theNode,"StartDate").GetJulian());

				theNode->selectSingleNode(_bstr_t("Coef"), &theNode2);
				if (theNode2 != NULL)
				{
					theNode2->get_text(&resultat);
					_bstr_t ff(resultat,false);
					char * ff1=(char *)ff;
					Spreads.push_back(atof(ff1));

					if (resultat) SysFreeString(resultat);
					theNode2->Release();
					theNode2 = NULL;
				}
			}
			if (theNode) theNode->Release();
			theNode = NULL;
		}
	}
	if (resultList) resultList->Release();
	resultList = NULL;

	if (nb == 1)
		newRefVal = new ARM_ReferenceValue(Spreads[0]);
	else
	{
		ARM_Vector* VRefDates= CreateARMVectorFromVECTOR(RefDates);
		ARM_Vector* VSpreads= CreateARMVectorFromVECTOR(Spreads);
		newRefVal = new ARM_ReferenceValue(VRefDates, VSpreads);
	}
	return newRefVal;
}

ARM_ReferenceValue*  GetCRFOneCpnMin(MSXML2::IXMLDOMNode* node)
{
	MSXML2::IXMLDOMNode* theNode=NULL;
	ARM_ReferenceValue* newRefVal = NULL;
	BSTR resultat = NULL;
	
	double Spread;
	node->selectSingleNode(_bstr_t("ProdData/cCRFARM_I/FloorValue"), &theNode);
	if (theNode!=NULL)
	{
		theNode->get_text(&resultat);
		_bstr_t ff(resultat,false);
		char * ff1=(char *)ff;
		Spread=atof(ff1);

		if (resultat) SysFreeString(resultat);
		theNode->Release();
		theNode=NULL;
	}

	newRefVal = new ARM_ReferenceValue(Spread);
	
	return newRefVal;
}

ARM_ReferenceValue*  GetCRFCpnMin(MSXML2::IXMLDOMNode* node)
{
	MSXML2::IXMLDOMNode* theNode=NULL;
	MSXML2::IXMLDOMNode* theNode2=NULL;
	MSXML2::IXMLDOMNodeList* resultList = NULL;
	ARM_ReferenceValue* newRefVal = NULL;
	BSTR resultat = NULL;
	HRESULT hr;
	long nb=0;
	
	VECTOR<double> RefDates;
	VECTOR<double> Spreads;
	node->selectNodes(_bstr_t("ProdData/cCRFARM_I/Schedule/cCRFLST_I"), &resultList);
	if (resultList!=NULL)
	{
		long nb;
		resultList->get_length(&nb);
		if (nb==0)
		{
			newRefVal = GetCRFOneCpnMin(node);
			if (resultList) resultList->Release();
			resultList=NULL;

			return newRefVal;
		}
		for (long i=0 ; i<nb ; i++)
		{	
			hr=resultList->get_item(i, &theNode);
			if (hr==S_OK && theNode!=NULL)
			{
				// voir si le refvalue se cale sur les StartDate ou sur les EndDate...
				RefDates.push_back(GetDateFromXMLNode(theNode,"StartDate").GetJulian());

				theNode->selectSingleNode(_bstr_t("FloorValue"), &theNode2);
				if (theNode2!=NULL)
				{
					theNode2->get_text(&resultat);
					_bstr_t ff(resultat,false);
					char * ff1=(char *)ff;
					Spreads.push_back(atof(ff1));

					if (resultat) SysFreeString(resultat);
					theNode2->Release();
					theNode2=NULL;
				}
			}
			if (theNode) theNode->Release();
			theNode = NULL;
		}
	}

	if (resultList) resultList->Release();
	resultList=NULL;
	if (nb==1)
		newRefVal = new ARM_ReferenceValue(Spreads[0]);
	else
	{
		ARM_Vector* VRefDates= CreateARMVectorFromVECTOR(RefDates);
		ARM_Vector* VSpreads= CreateARMVectorFromVECTOR(Spreads);
		newRefVal = new ARM_ReferenceValue(VRefDates, VSpreads);
	}
	return newRefVal;
}

ARM_ReferenceValue*  GetCRFOneCpnMax(MSXML2::IXMLDOMNode* node)
{
	MSXML2::IXMLDOMNode* theNode=NULL;
	ARM_ReferenceValue* newRefVal = NULL;
	BSTR resultat = NULL;
	
	double Spread;
	node->selectSingleNode(_bstr_t("ProdData/cCRFARM_I/CapValue"), &theNode);
	if (theNode!=NULL)
	{
		theNode->get_text(&resultat);
		_bstr_t ff(resultat,false);
		char * ff1=(char *)ff;
		Spread=atof(ff1);

		if (resultat) SysFreeString(resultat);
		theNode->Release();
		theNode=NULL;
	}

	newRefVal = new ARM_ReferenceValue(Spread);
	
	return newRefVal;
}

ARM_ReferenceValue*  GetCRFCpnMax(MSXML2::IXMLDOMNode* node)
{
	MSXML2::IXMLDOMNode* theNode=NULL;
	MSXML2::IXMLDOMNode* theNode2=NULL;
	MSXML2::IXMLDOMNodeList* resultList = NULL;
	ARM_ReferenceValue* newRefVal = NULL;
	BSTR resultat = NULL;
	HRESULT hr;
	long nb=0;
	
	VECTOR<double> RefDates;
	VECTOR<double> Spreads;
	node->selectNodes(_bstr_t("ProdData/cCRFARM_I/Schedule/cCRFLST_I"), &resultList);
	if (resultList!=NULL)
	{
		long nb;
		resultList->get_length(&nb);
		if (nb==0)
		{
			newRefVal = GetCRFOneCpnMax(node);				
			if (resultList) resultList->Release();
			resultList=NULL;

			return newRefVal;
		}
		for (long i=0 ; i<nb ; i++)
		{	
			hr=resultList->get_item(i, &theNode);
			if (hr==S_OK && theNode!=NULL)
			{
				// voir si le refvalue se cale sur les StartDate ou sur les EndDate...
				RefDates.push_back(GetDateFromXMLNode(theNode,"StartDate").GetJulian());

				theNode->selectSingleNode(_bstr_t("CapValue"), &theNode2);
				if (theNode2!=NULL)
				{
					theNode2->get_text(&resultat);
					_bstr_t ff(resultat,false);
					char * ff1=(char *)ff;
					Spreads.push_back(atof(ff1));

					if (resultat) SysFreeString(resultat);
					theNode2->Release();
					theNode2=NULL;
				}
			}
			if (theNode) theNode->Release();
			theNode = NULL;
		}
	}
	if (resultList) resultList->Release();
	resultList=NULL;
	
	if (nb==1)
		newRefVal = new ARM_ReferenceValue(Spreads[0]);
	else
	{
		ARM_Vector* VRefDates = CreateARMVectorFromVECTOR(RefDates);
		ARM_Vector* VSpreads  = CreateARMVectorFromVECTOR(Spreads);
		newRefVal = new ARM_ReferenceValue(VRefDates, VSpreads);
	}
	return newRefVal;
}


//***************************************************************************************************
// Appel au Niveau ASSET
//***************************************************************************************************

int GetStubRule(MSXML2::IXMLDOMNode* node)
{	
	MSXML2::IXMLDOMNode* theNode=NULL;
	BSTR resultat = NULL;
	int value;

	node->selectSingleNode(_bstr_t("STUB_StubType"), &theNode);
	if (theNode!=NULL)
	{
		theNode->get_text(&resultat);

		_bstr_t ff(resultat,false);
		char * ff1=(char *)ff;
			
		value = ARM_ConvStubRule(ff1);


		theNode->Release();
		theNode=NULL;
		if (resultat) SysFreeString(resultat);
	}
	else 
		value=K_SHORTSTART;
	return value;
}

int GetAssetDayCount(MSXML2::IXMLDOMNode* node)
{
	MSXML2::IXMLDOMNode* theNode=NULL;
	BSTR resultat = NULL;
	int value=KNOBASE;
	
	node->selectSingleNode(_bstr_t("INTEREST_Basis"), &theNode);
	if (theNode!=NULL)
	{
		theNode->get_text(&resultat);

		_bstr_t ff(resultat,false);
		char * ff1=(char *)ff;
		
		if (strcmp((const char*) ff1,"") != 0)
			value = FromSummitDaycountToARMDaycount((const char*) ff1);
		
		theNode->Release();
		theNode=NULL;
		if (resultat) SysFreeString(resultat);
	}
	return value;
}

int GetInterestDayCount(MSXML2::IXMLDOMNode* node)
{
	MSXML2::IXMLDOMNode* theNode=NULL;
	BSTR resultat = NULL;
	int value=KNOBASE;
	
	node->selectSingleNode(_bstr_t("INT_ACC_Basis"), &theNode);
	if (theNode!=NULL)
	{
		theNode->get_text(&resultat);

		_bstr_t ff(resultat,false);
		char * ff1=(char *)ff;
		
		if (strcmp((const char*) ff1,"") != 0)
			value = FromSummitDaycountToARMDaycount((const char*) ff1);
		
		theNode->Release();
		theNode=NULL;
		if (resultat) SysFreeString(resultat);
	}
	return value;
}

int GetFwdRule(MSXML2::IXMLDOMNode* node)
{
	MSXML2::IXMLDOMNode* theNode=NULL;
	BSTR resultat = NULL;
	int value=KNOBASE;
	
	node->selectSingleNode(_bstr_t("SCHED_Reset_Rule"), &theNode);
	if (theNode!=NULL)
	{
		theNode->get_text(&resultat);

		_bstr_t ff(resultat,false);
		char * ff1=(char *)ff;
		
		if (strcmp((const char*) ff1,"") != 0)
			value = ARM_ConvFwdRule((const char*) ff1);
		
		theNode->Release();
		theNode=NULL;
		if (resultat) SysFreeString(resultat);
	}
	return value;
}

int GetScheduleRule(MSXML2::IXMLDOMNode* node)
{
	MSXML2::IXMLDOMNode* theNode=NULL;
	BSTR resultat = NULL;
	int value=KNOBASE;
	
	node->selectSingleNode(_bstr_t("SCHED_Reset_Rule"), &theNode);
	if (theNode!=NULL)
	{
		theNode->get_text(&resultat);

		_bstr_t ff(resultat,false);
		char * ff1=(char *)ff;
		
		if (strcmp((const char*) ff1,"") != 0)
			value = ARM_ConvFwdRule((const char*) ff1);
		
		theNode->Release();
		theNode=NULL;
		if (resultat) SysFreeString(resultat);
	}
	return value;
}

int GetPayIndexFwdRule(MSXML2::IXMLDOMNode* node)
{
	MSXML2::IXMLDOMNode* theNode	= NULL;
	BSTR resultat			= NULL;
	int value				= KNOBASE;
	
	node->selectSingleNode(_bstr_t("SCHED_Pay_Rule"), &theNode);
	if (theNode != NULL)
	{
		theNode->get_text(&resultat);

		_bstr_t ff(resultat,false);
		char * ff1 = (char *)ff;
		
		if (strcmp((const char*) ff1 , "") != 0)
			value = ARM_ConvFwdRule((const char*) ff1);
		
		theNode->Release();
		theNode = NULL;
		if (resultat) SysFreeString(resultat);
	}
	return value;
}

ARM_ReferenceValue* GetNotionalConst(MSXML2::IXMLDOMNode* node)
{
	MSXML2::IXMLDOMNode* theNode=NULL;
	BSTR resultat = NULL;
	double notional;

	notional = GetDoubleFromXMLNode(node,"Notional");

	return new ARM_ReferenceValue(notional);
}

ARM_Date GetStartDate(MSXML2::IXMLDOMNode* node)
{
	return GetDateFromXMLNode(node,"EffDate");
}

ARM_Date GetEndDate(MSXML2::IXMLDOMNode* node)
{
	return GetDateFromXMLNode(node,"MatDate");
}

ARM_Date GetFstEffDate(MSXML2::IXMLDOMNode* node,ARM_Date defaultDate)
{
	ARM_Date fstCpnEffDate (defaultDate);
	MSXML2::IXMLDOMNode* theNode=NULL;

	node->selectSingleNode(_bstr_t("FstCpEffect"), &theNode);
	if (theNode!=NULL)
	{
		BSTR resultat = NULL;
		theNode->get_text(&resultat);

		_bstr_t ff(resultat,false);
		char * ff1=(char *)ff;

		if (strcmp((const char*)ff1,"") != 0)
		{
			ARM_Date tmpDate (ff1,"YYYYMMDD");

			fstCpnEffDate = tmpDate;
		}

		theNode->Release();
		theNode=NULL;
		if (resultat) SysFreeString(resultat);
	}

	return fstCpnEffDate;
}

ARM_ReferenceValue* GetNotional(MSXML2::IXMLDOMNode* node, bool isKernel)
{
	MSXML2::IXMLDOMNode* theNode=NULL;
	MSXML2::IXMLDOMNode* listItem=NULL;
	MSXML2::IXMLDOMNodeList* resultList = NULL;
	ARM_ReferenceValue* newRefVal = NULL;
	BSTR resultat = NULL;
	HRESULT hr;
	long nbEvents;
	
	VECTOR<double> notionalDate;
	VECTOR<double> notionalVal;
	
	if (! isKernel)
	{
		ARM_Date startDate=GetStartDate(node);
		notionalDate.push_back(startDate.GetJulian());
	}

	ARM_Date endDate=GetEndDate(node);
	node->selectSingleNode(_bstr_t("Notional"), &theNode);
	if (theNode!=NULL)
	{
		theNode->get_text(&resultat);

		_bstr_t ff(resultat,false);
		char * ff1=(char *)ff;
		double notional = atof((const char*) ff1);
		notionalVal.push_back(notional);
		theNode->Release();
		theNode=NULL;
		if (resultat) SysFreeString(resultat);
	}	
	node->selectNodes(_bstr_t("Events/EVENT"), &resultList);
	if (resultList!=NULL)
	{
		resultList->get_length(&nbEvents);
		CCString type;

		for (int i = 0; i < nbEvents; i++)
		{
			hr=resultList->get_item(i, &listItem);
			listItem->selectSingleNode(_bstr_t("Type"), &theNode);
			if (theNode!=NULL)
			{
				BSTR resultat = NULL;
				theNode->get_text(&resultat);

				_bstr_t ff(resultat,false);
				char * ff1=(char *)ff;

				type = (const char*) ff1;

				theNode->Release();
				theNode=NULL;
				if (resultat) SysFreeString(resultat);
			}
//***********************************************************
// WARNING
// Parse XML sur la Prod Type= XNL ie Notional Exchange si Amort.
// Parse XML sur la Recette Type= NTL ie Notional si Amort
// Parse XML sur la Recette Type= XNL ie Notional si Amort + Notional Exchange
// verif contenu du Amount si nominal ou adjutement de nominal...
//***********************************************************

			if ((strcmp(type,"NTL") == 0) || (strcmp(type,"XNL") == 0))
			{
				notionalDate.push_back(GetDateFromXMLNode(listItem,"ADate").GetJulian());

				listItem->selectSingleNode(_bstr_t("Amount"), &theNode);
				if (theNode!=NULL)
				{
					BSTR resultat = NULL;
					theNode->get_text(&resultat);

					_bstr_t ff(resultat,false);
					char * ff1=(char *)ff;

					double strike = atof(ff1);
					notionalVal.push_back(strike);

					theNode->Release();
					theNode=NULL;
					if (resultat) SysFreeString(resultat);
				}
			}
			if (listItem)
				listItem->Release();
		}
		if (resultList)
			resultList->Release();
	}

	if (isKernel)
		notionalDate.push_back(endDate.GetJulian());
	
	if (notionalVal.size()==1)
		newRefVal = new ARM_ReferenceValue(notionalVal[0]);
	else
	{
		ARM_Vector* VRefDates= CreateARMVectorFromVECTOR(notionalDate);
		ARM_Vector* VNotional= CreateARMVectorFromVECTOR(notionalVal);
		newRefVal = new ARM_ReferenceValue((ARM_Vector*)VRefDates->Clone(), (ARM_Vector*)VNotional->Clone());
		newRefVal->SetCalcMethod(K_STEPUP_RIGHT);
		newRefVal->SetExtrapolMeth(1);
		delete VRefDates;
		VRefDates=NULL;
		delete VNotional;
		VNotional=NULL;

	}
	return newRefVal;
}



ARM_ReferenceValue* GetStepUpFixCoupon(MSXML2::IXMLDOMNode* node)
{
	MSXML2::IXMLDOMNode* theNode=NULL;
	MSXML2::IXMLDOMNode* listItem=NULL;
	MSXML2::IXMLDOMNodeList* resultList = NULL;
	ARM_ReferenceValue* newRefVal = NULL;
	BSTR resultat = NULL;
	HRESULT hr;
	long nbEvents;
	
	VECTOR<double> cpnDate;
	VECTOR<double> cpnVal;
	
	char tmpPayCal[4];
	GetPayCalendar(node,tmpPayCal);

	cpnVal.push_back(GetDoubleFromXMLNode(node,"INTEREST_Rate") * 100.);
	cpnDate.push_back(GetStartDate(node).GapBusinessDay(-2, tmpPayCal).GetJulian());

	node->selectNodes(_bstr_t("Events/EVENT"), &resultList);
	if (resultList!=NULL)
	{
		resultList->get_length(&nbEvents);
		CCString type;

		for (int i = 0; i < nbEvents; i++)
		{
			hr=resultList->get_item(i, &listItem);
			listItem->selectSingleNode(_bstr_t("Type"), &theNode);
			if (theNode!=NULL)
			{
				BSTR resultat = NULL;
				theNode->get_text(&resultat);

				_bstr_t ff(resultat,false);
				char * ff1=(char *)ff;

				type = (const char*) ff1;

				theNode->Release();
				theNode=NULL;
				if (resultat) SysFreeString(resultat);
			}

			if (strcmp(type,"FRC") == 0)
			{
				ARM_Date tmpDate = GetDateFromXMLNode(listItem,"Date");

				tmpDate.GapBusinessDay(-2, tmpPayCal); //pour un index fixe, c'est toujours -2
				cpnDate.push_back(tmpDate.GetJulian());

				cpnVal.push_back(GetDoubleFromXMLNode(listItem,"Amount") * 100.);				
			}
			if (listItem)
				listItem->Release();
		}
		if (resultList)
			resultList->Release();
	}

	if (cpnVal.size()==1)
		newRefVal = new ARM_ReferenceValue(cpnVal[0]);
	else
	{
		ARM_Vector* VCpnDates= CreateARMVectorFromVECTOR(cpnDate);
		ARM_Vector* VCpnVal= CreateARMVectorFromVECTOR(cpnVal);
		newRefVal = new ARM_ReferenceValue((ARM_Vector*)VCpnDates->Clone(), (ARM_Vector*)VCpnVal->Clone());
		newRefVal->SetCalcMethod(K_STEPUP_LEFT);
		newRefVal->SetExtrapolMeth(1);
		delete VCpnDates;
		VCpnDates=NULL;
		delete VCpnVal;
		VCpnVal=NULL;

	}
	return newRefVal;
}


ARM_ReferenceValue* GetNotionalCust(MSXML2::IXMLDOMNode* node)
{
	MSXML2::IXMLDOMNode* theNode=NULL;
	MSXML2::IXMLDOMNode* listItem=NULL;
	MSXML2::IXMLDOMNodeList* resultList = NULL;
	ARM_ReferenceValue* newRefVal = NULL;
	BSTR resultat = NULL;
	HRESULT hr;
	long nbEvents;
	
	VECTOR<double> notionalDate;
	VECTOR<double> notionalVal;

	double tmpNotional = 0;
	int isFirstNotional = 1;
	int PorS = GetPorS(node);

	ARM_Date endDate=GetEndDate(node);

	node->selectSingleNode(_bstr_t("Notional"), &theNode);
	if (theNode!=NULL)
	{
		theNode->get_text(&resultat);

		_bstr_t ff(resultat,false);
		char * ff1=(char *)ff;
		double notional = atof((const char*) ff1);
		notionalVal.push_back(notional);
		tmpNotional = notional;
		theNode->Release();
		theNode=NULL;
		if (resultat) SysFreeString(resultat);
	}

	node->selectNodes(_bstr_t("Events/EVENT"), &resultList);
	if (resultList!=NULL)
	{
		resultList->get_length(&nbEvents);
		CCString type;

		for (int i = 0; i < nbEvents; i++)
		{
			hr=resultList->get_item(i, &listItem);
			listItem->selectSingleNode(_bstr_t("Type"), &theNode);
			if (theNode!=NULL)
			{
				BSTR resultat = NULL;
				theNode->get_text(&resultat);

				_bstr_t ff(resultat,false);
				char * ff1=(char *)ff;

				type = (const char*) ff1;

				theNode->Release();
				theNode=NULL;
				if (resultat) SysFreeString(resultat);
			}
//***********************************************************
// WARNING
// Parse XML sur la Prod Type= XNL ie Notional Exchange si Amort.
// Parse XML sur la Recette Type= NTL ie Notional si Amort
// Parse XML sur la Recette Type= XNL ie Notional si Amort + Notional Exchange
// verif contenu du Amount si nominal ou adjutement de nominal...
//***********************************************************

			if ((strcmp(type,"NTL") == 0) || (strcmp(type,"XNL") == 0))
			{
				if (isFirstNotional == 0)
				{
					notionalDate.push_back(GetDateFromXMLNode(listItem,"ADate").GetJulian());

					listItem->selectSingleNode(_bstr_t("Amount"), &theNode);
					if (theNode!=NULL)
					{
						BSTR resultat = NULL;
						theNode->get_text(&resultat);

						_bstr_t ff(resultat,false);
						char * ff1=(char *)ff;

						double notional = atof(ff1);
						notionalVal.push_back(tmpNotional - PorS * notional);
						tmpNotional -= PorS * notional;

						theNode->Release();
						theNode=NULL;
						if (resultat) SysFreeString(resultat);
					}
				}
				isFirstNotional = 0;
			}
			if (listItem)
				listItem->Release();
		}
		if (resultList)
			resultList->Release();
	}
	notionalDate.push_back(endDate.GetJulian());

	if (notionalVal.size() == 2)
		newRefVal = new ARM_ReferenceValue(notionalVal[0]);
	else
	{
		ARM_Vector* VRefDates = CreateARMVectorFromVECTOR(notionalDate,notionalDate.size() - 1);
		ARM_Vector* VNotional = CreateARMVectorFromVECTOR(notionalVal,notionalVal.size() - 1);
		newRefVal = new ARM_ReferenceValue((ARM_Vector*)VRefDates->Clone(), (ARM_Vector*)VNotional->Clone());
		newRefVal->SetCalcMethod(K_STEPUP_RIGHT);
		newRefVal->SetExtrapolMeth(1);
		delete VRefDates;
		VRefDates=NULL;
		delete VNotional;
		VNotional=NULL;

	}
	return newRefVal;
}


ARM_ReferenceValue* GetNotionalWithPayDates(MSXML2::IXMLDOMNode* node, bool isCust)
{
	MSXML2::IXMLDOMNode        *theNode    = NULL;
	MSXML2::IXMLDOMNode        *theNode2    = NULL;
	MSXML2::IXMLDOMNode        *listItem   = NULL;
	MSXML2::IXMLDOMNodeList    *resultList = NULL;
	BSTR               resultat    = NULL;
	HRESULT            hr;

	VECTOR<double>     notionalDate;
	VECTOR<double>     notionalVal;
	double			   firstExchNot;
	double			   varNotional;
	double			   currDate;
	long               nbEvents;
	
	ARM_ReferenceValue *newRefVal  = NULL;
	ARM_Currency	   *ccy        = NULL;
	ARM_Vector         *payDates   = NULL;
	ARM_Vector         *startDates = NULL;
	
	// generation du schedule de payment et des endDates correspondantes
	ccy					= GetCcy(node);
	ARM_Date startDate  = GetStartDate(node);
	ARM_Date endDate    = GetEndDate(node);
	int payFreq		    = GetPayFreq(node);
	int payTiming		= GetPayTiming(node);
	int indexType		= FromFrequencyToXiborType(payFreq, ccy->GetCcyName()); // a changer en ARM_INDEX_TYPE, ex : EURIBOR6M
	int intRule			= GetPayIntRule(node);
	int idx;

	ARM_SwapLeg sched(startDate, endDate, (ARM_INDEX_TYPE)indexType, K_RCV, 0.0, K_DEF_FREQ, payFreq, K_ADVANCE, payTiming, ccy, intRule);
	delete ccy;

	payDates   = sched.GetPaymentDates();
	startDates = sched.GetFlowStartDates();

    if ( !payDates || !startDates )
	{
		CCString msg((CCString)"Error generating payment schedule");
		throw Exception (__LINE__, __FILE__, ERR_OBJECT_NULL, (char*) msg);
	}		

	//Nominal de depart
	node->selectSingleNode(_bstr_t("Notional"), &theNode);
	if (theNode!=NULL)
	{
		theNode->get_text(&resultat);
		_bstr_t ff(resultat,false);
		char* ff1 = (char *)ff;
		double notional = atof((const char*) ff1);
		notionalVal.push_back(notional);
		theNode->Release();
		theNode=NULL;
		if (resultat) SysFreeString(resultat);
	}

	// la 1ere date de paiement correspond au nominal de depart
	notionalDate.push_back(payDates->Elt(0));

	//Nominaux suivants si non constants
	node->selectNodes(_bstr_t("Events/EVENT"), &resultList);
	if (resultList!=NULL)
	{
		resultList->get_length(&nbEvents);
		CCString type;

		for (int i = 0; i < nbEvents; i++)
		{
			hr=resultList->get_item(i, &listItem);
			listItem->selectSingleNode(_bstr_t("Type"), &theNode);
			if (theNode!=NULL)
			{
				BSTR resultat = NULL;
				theNode->get_text(&resultat);
				_bstr_t ff(resultat,false);
				char * ff1=(char *)ff;
				type = (const char*) ff1;
				theNode->Release();
				theNode = NULL;
				if (resultat) SysFreeString(resultat);
			}
//***********************************************************
// WARNING
// Parse XML sur la Prod Type= XNL ie Notional Exchange si Amort.
// Parse XML sur la Recette Type= NTL ie Notional si Amort
// Parse XML sur la Recette Type= XNL ie Notional si Amort + Notional Exchange
// verif contenu du Amount si nominal ou adjutement de nominal...
//***********************************************************
			if ((strcmp(type,"NTL") == 0) || (strcmp(type,"XNL") == 0))
			{
				// date de paiement correspondante a la enddate ds EVENT
				//Le XML contient les évent NTL avec (dateDépart, Ntl) correspondant.
				currDate = GetDateFromXMLNode(listItem,"ADate").GetJulian();

				//We test the notional exchange (not taken in account)
				if (isCust)
				{
					if (currDate == startDate.GetJulian())
					{
						listItem->selectSingleNode(_bstr_t("Amount"), &theNode);
						if (theNode != NULL)
						{
							BSTR resultat = NULL;
							theNode->get_text(&resultat);
							_bstr_t ff(resultat,false);
							char * ff1=(char *)ff;
							firstExchNot = atof(ff1);
							theNode->Release();
							theNode=NULL;
							if (resultat) SysFreeString(resultat);
						}	
					}
					else if ((currDate != startDate.GetJulian()) && (currDate != endDate.GetJulian()))
					{
						idx = startDates->find(currDate);
						if (idx == -1)
						{
							CCString msg((CCString)"Error: payment sched and event sched don't match !");
							throw Exception (__LINE__, __FILE__, ERR_CONTAINER_INDEX, (char*) msg);
						}
						notionalDate.push_back(payDates->Elt(idx));

						listItem->selectSingleNode(_bstr_t("Amount"), &theNode);
						if (theNode != NULL)
						{
							BSTR resultat = NULL;
							theNode->get_text(&resultat);
							_bstr_t ff(resultat,false);
							char * ff1=(char *)ff;
							if (firstExchNot > 0)
							{
								varNotional = firstExchNot + atof(ff1);
								firstExchNot = varNotional;
							}
							else if (firstExchNot < 0)
							{
								varNotional = -firstExchNot - atof(ff1);
								firstExchNot = -varNotional;
							}
							notionalVal.push_back(varNotional);
							theNode->Release();
							theNode=NULL;
							if (resultat) SysFreeString(resultat);
						}	
					}
				}
				else if (!isCust)
				{
					idx = startDates->find(currDate);
					if (idx == -1)
					{
						CCString msg((CCString)"Error: payment sched and event sched don't match !");
						throw Exception (__LINE__, __FILE__, ERR_CONTAINER_INDEX, (char*) msg);
					}
					notionalDate.push_back(payDates->Elt(idx));
				
					listItem->selectSingleNode(_bstr_t("Amount"), &theNode);
					if (theNode != NULL)
					{
						BSTR resultat = NULL;
						theNode->get_text(&resultat);
						_bstr_t ff(resultat,false);
						char * ff1=(char *)ff;
						double currNotional = atof(ff1);
						notionalVal.push_back(currNotional);
						theNode->Release();
						theNode=NULL;
						if (resultat) SysFreeString(resultat);
					}
				}

			}
			if (listItem)
				listItem->Release();
		}
		if (resultList)
			resultList->Release();
	}

	if (notionalVal.size() == 1)
		newRefVal = new ARM_ReferenceValue(notionalVal[0]);
	else
	{
		ARM_Vector* VRefDates= CreateARMVectorFromVECTOR(notionalDate);
		ARM_Vector* VNotional= CreateARMVectorFromVECTOR(notionalVal);
		newRefVal = new ARM_ReferenceValue((ARM_Vector*)VRefDates->Clone(), (ARM_Vector*)VNotional->Clone());
		newRefVal->SetCalcMethod(K_STEPUP_LEFT);
		newRefVal->SetExtrapolMeth(1);
		delete VRefDates;
		VRefDates=NULL;
		delete VNotional;
		VNotional=NULL;

	}
	return newRefVal;
}

ARM_ReferenceValue* GetFundNotionalWithPayDates(MSXML2::IXMLDOMNode* node)
{
	MSXML2::IXMLDOMNode        *theNode    = NULL;
	MSXML2::IXMLDOMNode        *theNode2   = NULL;
	MSXML2::IXMLDOMNode        *listItem   = NULL;
	MSXML2::IXMLDOMNodeList    *resultList = NULL;
	BSTR               resultat    = NULL;
	HRESULT            hr;

	VECTOR<double>     notionalDate;
	VECTOR<double>     notionalVal;
	long               nbEvents;
	ARM_ReferenceValue *newRefVal  = NULL;
	double firstNot	   = 0;
	double currNot     = 0;
	//ARM_Date endDate   = GetEndDate(node);

	// nominal de depart
	node->selectSingleNode(_bstr_t("ProdData/cBLOB_I/Amount4"), &theNode);
	if (theNode!=NULL)
	{
		theNode->get_text(&resultat);

		_bstr_t ff(resultat,false);
		char* ff1 = (char *)ff;
		firstNot = atof((const char*) ff1);
		
		if (resultat) SysFreeString(resultat);
	}

	node->selectNodes(_bstr_t("ProdData/cBLOB_I/Schedule2/cBLBLST_I"), &resultList);
	//Nominal Amorti
	if (resultList != NULL)
	{
		resultList->get_length(&nbEvents);
		for (int i = 0; i < nbEvents; i++)
		{
			//End Date:
			hr = resultList->get_item(i, &listItem);
			listItem->selectSingleNode(_bstr_t("EndDate"), &theNode2);
			if (theNode2 != NULL)
			{
				BSTR resultat = NULL;
				theNode2->get_text(&resultat);

				_bstr_t ff(resultat,false);
				char * ff1=(char *)ff;
				ARM_Date endDate (ff1,"YYYYMMDD");
				notionalDate.push_back(endDate.GetJulian());
				
				theNode2->Release();
				theNode2=NULL;
				if (resultat) SysFreeString(resultat);
			}
			//Corresponding Notional:
			listItem->selectSingleNode(_bstr_t("Amount2"), &theNode2);
			if (theNode2!=NULL)
			{
				BSTR resultat = NULL;
				theNode2->get_text(&resultat);

				_bstr_t ff(resultat,false);
				char * ff1=(char *)ff;
				currNot = atof(ff1);
				notionalVal.push_back(currNot);
				
				theNode2->Release();
				theNode2 = NULL;
				if (resultat) SysFreeString(resultat);
			}

			if (listItem)
				listItem->Release();
		}
		if (resultList)
			resultList->Release();
	}
	//Il y a peut-être un seul nominal
	if ((theNode != NULL) && (nbEvents == 0))
	{
		//notionalDate.push_back(endDate.GetJulian());
		notionalVal.push_back(firstNot);
		theNode->Release();
		theNode=NULL;
	}

	if (notionalVal.size() == 1)
		newRefVal = new ARM_ReferenceValue(notionalVal[0]);
	else if (notionalVal.size() > 1)
	{
		ARM_Vector* VRefDates= CreateARMVectorFromVECTOR(notionalDate);
		ARM_Vector* VNotional= CreateARMVectorFromVECTOR(notionalVal);
		newRefVal = new ARM_ReferenceValue((ARM_Vector*)VRefDates->Clone(), (ARM_Vector*)VNotional->Clone());
		newRefVal->SetCalcMethod(K_STEPUP_LEFT);
		newRefVal->SetExtrapolMeth(1);
		delete VRefDates;
		VRefDates=NULL;
		delete VNotional;
		VNotional=NULL;
	}
	return newRefVal;
}


ARM_Currency* GetCcy(MSXML2::IXMLDOMNode* node)
{
	MSXML2::IXMLDOMNode* theNode=NULL;
	BSTR resultat = NULL;
	ARM_Currency* ccy=NULL;

	node->selectSingleNode(_bstr_t("Ccy"), &theNode);
	if (theNode!=NULL)
	{
		theNode->get_text(&resultat);
		ccy = new ARM_Currency((char*)(_bstr_t)resultat);

		theNode->Release();
		theNode=NULL;
		if (resultat) SysFreeString(resultat);
	}	
	return ccy;
}

// CMT underlying bond frequency
int GetYieldFreq(MSXML2::IXMLDOMNode* node)
{
	MSXML2::IXMLDOMNode* theNode=NULL;
	BSTR resultat = NULL;
	int value= K_DEF_FREQ;

	node->selectSingleNode(_bstr_t("INTEREST_YldFreq"), &theNode);
	if (theNode!=NULL)
	{
		theNode->get_text(&resultat);
		_bstr_t ff(resultat,false);
		char * ff1=(char *)ff;

		value = FromSummitFreqToARMFreq((const char*) ff1);

		if (resultat) SysFreeString(resultat);
		theNode->Release();
		theNode=NULL;
	}
	return value;
}

int GetResetFreq(MSXML2::IXMLDOMNode* node)
{
	MSXML2::IXMLDOMNode* theNode=NULL;
	BSTR resultat = NULL;
	int value= K_DEF_FREQ;

	node->selectSingleNode(_bstr_t("SCHED_Reset_Freq"), &theNode);
	if (theNode!=NULL)
	{
		theNode->get_text(&resultat);
		_bstr_t ff(resultat,false);
		char * ff1=(char *)ff;

		value = FromSummitFreqToARMFreq((const char*) ff1);

		if (resultat) SysFreeString(resultat);
		theNode->Release();
		theNode=NULL;
	}
	return value;
}

int GetPayFreq(MSXML2::IXMLDOMNode* node)
{
	MSXML2::IXMLDOMNode* theNode=NULL;
	BSTR resultat = NULL;
	int value= K_DEF_FREQ;

	node->selectSingleNode(_bstr_t("SCHED_Pay_Freq"), &theNode);
	if (theNode!=NULL)
	{
		theNode->get_text(&resultat);
		_bstr_t ff(resultat,false);
		char * ff1=(char *)ff;

		value = FromSummitFreqToARMFreq((const char*) ff1);

		if (resultat) SysFreeString(resultat);
		theNode->Release();
		theNode=NULL;
	}
	return value;
}

int GetCompoundingFreq(IXMLDOMNode* node)
{
	IXMLDOMNode* theNode=NULL;
	BSTR resultat = NULL;
	int value= K_DEF_FREQ;

	node->selectSingleNode(_bstr_t("COMP_CompFreq"), &theNode);
	if (theNode!=NULL)
	{
		theNode->get_text(&resultat);
		_bstr_t ff(resultat,false);
		char * ff1=(char *)ff;

		if (strcmp((const char*) ff1,"") != 0)
			value = FromSummitFreqToARMFreq((const char*) ff1);

		if (resultat) SysFreeString(resultat);
		theNode->Release();
		theNode=NULL;
	}
	return value;
}

int GetResetGap(MSXML2::IXMLDOMNode* node)
{
	MSXML2::IXMLDOMNode* theNode=NULL;
	BSTR resultat = NULL;
	int value=0;
	
	node->selectSingleNode(_bstr_t("SCHED_Reset_Gap"), &theNode);
	if (theNode!=NULL)
	{
		BSTR resultat = NULL;
		theNode->get_text(&resultat);

		_bstr_t ff(resultat,false);
		char * ff1=(char *)ff;

		value = FromSummitGapToARMGap((const char*) ff1);

		theNode->Release();
		theNode=NULL;
		if (resultat) SysFreeString(resultat);
	}
	return value;

}

int GetPayResetGap(MSXML2::IXMLDOMNode* node)			// for corridor leg
{
	MSXML2::IXMLDOMNode* theNode=NULL;
	BSTR resultat = NULL;
	int value=0;
	
	node->selectSingleNode(_bstr_t("INT_RES_NTL_Gap"), &theNode);
	if (theNode!=NULL)
	{
		BSTR resultat = NULL;
		theNode->get_text(&resultat);

		_bstr_t ff(resultat,false);
		char * ff1=(char *)ff;

		value = FromSummitGapToARMGap((const char*) ff1);

		theNode->Release();
		theNode=NULL;
		if (resultat) SysFreeString(resultat);
	}
	return value;

}

int GetPayGap(MSXML2::IXMLDOMNode* node)
{
	MSXML2::IXMLDOMNode* theNode = NULL;
	BSTR resultat		 = NULL;
	int value			 = 0;
	
	node->selectSingleNode(_bstr_t("SCHED_Pay_Gap"), &theNode);
	if (theNode != NULL)
	{
		BSTR resultat = NULL;
		theNode->get_text(&resultat);

		_bstr_t ff(resultat,false);
		char * ff1 = (char *)ff;
		
		value = FromSummitGapToARMGap((const char*) ff1);

		theNode->Release();
		theNode = NULL;
		if (resultat) SysFreeString(resultat);
	}
	return value;

}

int GetResetTiming(MSXML2::IXMLDOMNode* node)
{
	MSXML2::IXMLDOMNode* theNode=NULL;
	BSTR resultat = NULL;
	int value=K_ADVANCE;
	
	node->selectSingleNode(_bstr_t("SCHED_Reset_Time"), &theNode);
	if (theNode!=NULL)
	{
		BSTR resultat = NULL;
		theNode->get_text(&resultat);

		_bstr_t ff(resultat,false);
		char * ff1=(char *)ff;
		value = ARM_ConvPayResetRule(CCString((char *)ff));

		theNode->Release();
		theNode=NULL;
		if (resultat) SysFreeString(resultat);
	}
	return value;

}

int GetIfResrtikableType(MSXML2::IXMLDOMNode* node, int& payIndexType, int& refIndexType)
{
	MSXML2::IXMLDOMNode* theNode = NULL;
	BSTR result = NULL;
	int value = 1;
	refIndexType = K_CMS;

	node->selectSingleNode(_bstr_t("ProductName"), &theNode);
	if( theNode!=NULL)
	{
		theNode->get_text(&result);
		_bstr_t ff(result, false);
		char* ff1=(char* )ff;
		if(strcmp(ff1,"RESTRIKABLECMSLIBOR") == 0)
			payIndexType = K_LIBOR;
		else if(strcmp(ff1,"RESTRIKABLECMSCMS") == 0)
			payIndexType = K_CMS;
		else if(strcmp(ff1,"RESTRIKABLECMSFIXED") == 0)
			payIndexType = K_FIXED;
		else
			value = 0;

		theNode->Release();
		theNode = NULL;
		if(result)
			SysFreeString(result);
	}
	return value;
}

double GetRestrikableAlpha(MSXML2::IXMLDOMNode* node)
{
	char buffer[15];
	double alpha = 0.5;
	GetNthArgInFormula(node, 7, buffer);
	if(strcmp(buffer,"")!=0)
		alpha = atof(buffer);

	return alpha;
}

double GetRestrikableBeta(MSXML2::IXMLDOMNode* node)
{
	char buffer[15];
	double beta = 0.5;
	GetNthArgInFormula(node, 8, buffer);
	if(strcmp(buffer,"")!=0)
		beta = atof(buffer);

	return beta;
}

double GetRestrikableRange(MSXML2::IXMLDOMNode* node)
{
	char buffer[15];
	double range = 0.;
	GetNthArgInFormula(node, 3, buffer);
	if(strcmp(buffer,"")!=0)
		range = atof(buffer);

	return range;
}

int GetRestrikablePayIndex(MSXML2::IXMLDOMNode* node, ARM_Currency* ccy)
{
	char buffer[50];
	char payIndex[50];
	GetNthArgInFormula(node,2, buffer);
	if(strcmp(buffer, "")==0)
		return K_FIXED;
	ConvSummitIndexToARMIndex(buffer, ccy, payIndex);
	return ARM_ConvIrType(payIndex);
}

int GetRestrikableRefIndex(MSXML2::IXMLDOMNode* node, ARM_Currency* ccy)
{
	char buffer[50];
	char refIndex[50];
	GetNthArgInFormula(node,1, buffer);
	if(strcmp(buffer, "")==0)
		return K_FIXED;
	ConvSummitIndexToARMIndex(buffer, ccy, refIndex);
	return ARM_ConvIrType(refIndex);
}

			



int GetRngProductType(MSXML2::IXMLDOMNode* node)
{
	MSXML2::IXMLDOMNode* theNode=NULL;
	BSTR resultat = NULL;
	int value = K_FIXED;
	
	node->selectSingleNode(_bstr_t("ProductName"), &theNode);
	if (theNode!=NULL)
	{
		BSTR resultat = NULL;
		theNode->get_text(&resultat);

		_bstr_t ff(resultat,false);
		char * ff1=(char *)ff;
		if (strcmp(ff1, "RNGLIBORLIBOR") == 0)
			value = K_LIBOR;

		theNode->Release();
		theNode=NULL;
		if (resultat) SysFreeString(resultat);
	}
	return value;
}

int GetPayTiming(MSXML2::IXMLDOMNode* node)
{
	MSXML2::IXMLDOMNode* theNode=NULL;
	BSTR resultat = NULL;
	int value=K_ARREARS;
	
	node->selectSingleNode(_bstr_t("SCHED_Pay_Time"), &theNode);
	if (theNode!=NULL)
	{
		BSTR resultat = NULL;
		theNode->get_text(&resultat);

		_bstr_t ff(resultat,false);
		char * ff1=(char *)ff;
		value = ARM_ConvPayResetRule(CCString((char *)ff));

		theNode->Release();
		theNode=NULL;
		if (resultat) SysFreeString(resultat);
	}
	return value;

}

int GetPorS(MSXML2::IXMLDOMNode* node)
{
	MSXML2::IXMLDOMNode* theNode=NULL;
	BSTR resultat = NULL;
	int value= K_RCV;

	node->selectSingleNode(_bstr_t("PorS"), &theNode);
	if (theNode!=NULL)
	{
		theNode->get_text(&resultat);
		_bstr_t ff(resultat,false);
		char * ff1=(char *)ff;

		if (strcmp((const char*) ff1,"S") == 0)
			value = K_PAY;

		if (resultat) SysFreeString(resultat);
		theNode->Release();
		theNode=NULL;
	}
	return value;
}

int GetLDPricingMeth(MSXML2::IXMLDOMNode* node)
{	
	MSXML2::IXMLDOMNode* theNode=NULL;
	BSTR resultat = NULL;
	int value=K_DIGITALE;

	node->selectSingleNode(_bstr_t("OPForm"), &theNode);
	if (theNode!=NULL)
	{
		theNode->get_text(&resultat);
		_bstr_t ff(resultat,false);
		char * ff1=(char *)ff;

		value = ARM_ConvPricCorridorLD(ff1);

		theNode->Release();
		theNode=NULL;
		if (resultat) SysFreeString(resultat);
	}
	
	return value;
}

int GetIndexType(MSXML2::IXMLDOMNode* node)
{
	MSXML2::IXMLDOMNode* theNode=NULL;
	BSTR resultat = NULL;
	char summitIndex[15];
	char summitTerm[15];
	char armIndex[15];

	node->selectSingleNode(_bstr_t("INTEREST_dmIndex"), &theNode);
	if (theNode!=NULL)
	{
		BSTR resultat = NULL;
		theNode->get_text(&resultat);

		_bstr_t ff(resultat,false);
		char * ff1=(char *)ff;
		strcpy(summitIndex,ff1);
		
		theNode->Release();
		theNode=NULL;
	}
	if (resultat) SysFreeString(resultat);
	
	if (strcmp((const char*)summitIndex, "FORM") == 0)
	{
		//Index is defined in Formula !
		ARM_Currency* ccy = GetCcy(node);

		return GetSpreadIndexType1(node, ccy);
	}

	node->selectSingleNode(_bstr_t("INTEREST_Term"), &theNode);
	if (theNode!=NULL)
	{
		BSTR resultat = NULL;
		theNode->get_text(&resultat);

		_bstr_t ff(resultat,false);
		char * ff1=(char *)ff;
		strcpy(summitTerm,ff1);
		
		theNode->Release();
		theNode=NULL;
	}
	if (resultat) SysFreeString(resultat);
	
	if (strcmp(summitIndex,"")==0)
		return K_FIXED;
	else
	{
		if ((strcmp(summitIndex,"EUR3M") == 0) || (strcmp(summitIndex,"EUR12") == 0) || (strcmp(summitIndex,"EUR1M") == 0))
		{
			return ARM_ConvIrType(summitIndex);
		}
		if (strcmp(summitIndex,"EURIB") == 0)
		{
			strcpy(summitIndex,strcat(summitIndex, "OR"));
		}
		if ((strcmp(summitIndex,"STIBO") == 0) || (strcmp(summitIndex,"AUBB") == 0))
		{
			strcpy(summitIndex,"LIBOR");
		}

		strcpy(armIndex,strcat(summitIndex,summitTerm));
		return ARM_ConvIrType(armIndex);
	}
	
}


double GetSpread(MSXML2::IXMLDOMNode* node)
{	
	MSXML2::IXMLDOMNode* theNode=NULL;
	BSTR resultat = NULL;
	double value=0.0;

	node->selectSingleNode(_bstr_t("INTEREST_Spread"), &theNode);
	if (theNode!=NULL)
	{
		theNode->get_text(&resultat);
		_bstr_t ff(resultat,false);
		char * ff1=(char *)ff;
		
		value = atof((const char*) ff1);
		
		theNode->Release();
		theNode=NULL;
		if (resultat) SysFreeString(resultat);
	}
	return value;
}


ARM_ReferenceValue* GetSpreadVariable(MSXML2::IXMLDOMNode* node)
{
	MSXML2::IXMLDOMNode* theNode=NULL;
	MSXML2::IXMLDOMNode* theNode2=NULL;
	MSXML2::IXMLDOMNodeList* resultList = NULL;
	ARM_ReferenceValue* newRefVal = NULL;
	BSTR resultat = NULL;
	HRESULT hr;
	long nb=0;

	int trouveSpread = 0;
	
	VECTOR<double> RefDates;
	VECTOR<double> Spreads;

	node->selectNodes(_bstr_t("Events/EVENT"), &resultList);
	if (resultList!=NULL)
	{
		long nb;
		resultList->get_length(&nb);

		if (nb==0)
		{
			if (resultList) resultList->Release();
			resultList=NULL;

			return NULL;
		}

		for (long i=0 ; i<nb ; i++)
		{	
			hr=resultList->get_item(i, &theNode);
			if (hr==S_OK && theNode!=NULL)
			{
				theNode->selectSingleNode(_bstr_t("Type"), &theNode2);
				if (theNode2!=NULL)
				{
					theNode2->get_text(&resultat);
					_bstr_t ff(resultat,false);
					char * ff1=(char *)ff;

					if (strcmp(ff1,"SPR")==0)
					{
						if (trouveSpread == 0)
						{
							Spreads.push_back(GetSpread(node)*100.0);
							RefDates.push_back(GetStartDate(node).GetJulian());
							trouveSpread = 1;
						}

						Spreads.push_back(GetDoubleFromXMLNode(theNode,"Amount")*100.0);
						RefDates.push_back(GetDateFromXMLNode(theNode,"Date").GetJulian());
					}
					theNode2->Release();
					theNode2 = NULL;
				}
				theNode->Release();
				theNode = NULL;
			}
		}
	}

	if (resultList) resultList->Release();
	resultList=NULL;

	if (trouveSpread==0)
		return NULL;
	else
	{
		ARM_Vector* VRefDates= CreateARMVectorFromVECTOR(RefDates);
		ARM_Vector* VSpreads= CreateARMVectorFromVECTOR(Spreads);
		newRefVal = new ARM_ReferenceValue(VRefDates, VSpreads,K_YIELD,1);
		newRefVal->SetCalcMethod(K_STEPUP_LEFT);
	}
	return newRefVal;
}


void GetPayCalendar(MSXML2::IXMLDOMNode* node, char* payCal)
{
	MSXML2::IXMLDOMNode* theNode=NULL;
	BSTR resultat = NULL;
	
	node->selectSingleNode(_bstr_t("SCHED_Pay_Cal"), &theNode);
	if (theNode!=NULL)
	{
		BSTR resultat = NULL;
		theNode->get_text(&resultat);

		_bstr_t ff(resultat,false);
		char * ff1=(char *)ff;

		strcpy(payCal,ff1);

		theNode->Release();
		theNode=NULL;
		if (resultat) SysFreeString(resultat);
	}
}

void GetResetCalendar(MSXML2::IXMLDOMNode* node, char* resetCal)
{
	MSXML2::IXMLDOMNode* theNode=NULL;
	BSTR resultat = NULL;
	
	node->selectSingleNode(_bstr_t("SCHED_Reset_Cal"), &theNode);
	if (theNode!=NULL)
	{
		BSTR resultat = NULL;
		theNode->get_text(&resultat);

		_bstr_t ff(resultat,false);
		char * ff1=(char *)ff;

		strcpy(resetCal,ff1);

		theNode->Release();
		theNode=NULL;
		if (resultat) SysFreeString(resultat);
	}
}

int GetIntRule(MSXML2::IXMLDOMNode* node)
{	
	MSXML2::IXMLDOMNode* theNode=NULL;
	BSTR resultat = NULL;
	int value=K_ADJUSTED;
	// c'est la meme valeur qui se trouve dans SCHED_Pay_IntRule
	node->selectSingleNode(_bstr_t("SCHED_Reset_IntRule"), &theNode);
	if (theNode!=NULL)
	{
		theNode->get_text(&resultat);
		_bstr_t ff(resultat,false);
		char * ff1=(char *)ff;
		
		value = ARM_ConvIntRule(ff1);

		theNode->Release();
		theNode=NULL;
		if (resultat) SysFreeString(resultat);
	}
	
	return value;
}

int GetAmortFreq(MSXML2::IXMLDOMNode* node)
{
	MSXML2::IXMLDOMNode* theNode=NULL;
	BSTR resultat = NULL;
	int value= K_DEF_FREQ;

	node->selectSingleNode(_bstr_t("AMORT_Freq"), &theNode);
	if (theNode!=NULL)
	{
		theNode->get_text(&resultat);
		_bstr_t ff(resultat,false);
		char * ff1=(char *)ff;

		value = FromSummitFreqToARMFreq((const char*) ff1);

		if (resultat) SysFreeString(resultat);
		theNode->Release();
		theNode=NULL;
	}
	return value;
}


int GetPayIntRule(MSXML2::IXMLDOMNode* node)
{	
	MSXML2::IXMLDOMNode* theNode=NULL;
	BSTR resultat = NULL;
	int value=K_ADJUSTED;
	
	node->selectSingleNode(_bstr_t("SCHED_Pay_IntRule"), &theNode);
	if (theNode!=NULL)
	{
		theNode->get_text(&resultat);
		_bstr_t ff(resultat,false);
		char * ff1=(char *)ff;
		
		value = ARM_ConvIntRule(ff1);

		theNode->Release();
		theNode=NULL;
		if (resultat) SysFreeString(resultat);
	}
	
	return value;
}

int GetDecompFreq(MSXML2::IXMLDOMNode* node)
{	
	MSXML2::IXMLDOMNode* theNode=NULL;
	BSTR resultat = NULL;
	int value=K_COMP_PROP;
	
	node->selectSingleNode(_bstr_t("FreqDecompFreq"), &theNode);
	if (theNode!=NULL)
	{
		theNode->get_text(&resultat);
		_bstr_t ff(resultat,false);
		char * ff1=(char *)ff;
		if (strcmp(ff1,"")!=0)
			value = ARM_ConvDecompFrequency(ff1);

		theNode->Release();
		theNode=NULL;
		if (resultat) SysFreeString(resultat);
	}
	
	return value;
}

int GetNotionalExchangeFlag(MSXML2::IXMLDOMNode* node)
{	
	MSXML2::IXMLDOMNode* theNode=NULL;
	BSTR resultat = NULL;
	int value=K_NX_NONE;
	
	// c'est la meme valeur qui se trouve dans SCHED_Pay_IntRule
	node->selectSingleNode(_bstr_t("CROSS_CCY_NotExch"), &theNode);
	if (theNode!=NULL)
	{
		theNode->get_text(&resultat);
		_bstr_t ff(resultat,false);
		char * ff1=(char *)ff;
		if (strcmp(ff1,"")!=0)
			value = ARM_NotionalExchange(ff1);

		theNode->Release();
		theNode=NULL;
		if (resultat) SysFreeString(resultat);
	}

	return value;
}

double GetMeanRev(MSXML2::IXMLDOMNode* node)
{	
	MSXML2::IXMLDOMNode* theNode=NULL;
	BSTR resultat = NULL;
	double value=10000.0;

	node->selectSingleNode(_bstr_t("INTEREST_FundSprd"), &theNode);
	if (theNode!=NULL)
	{
		theNode->get_text(&resultat);
		_bstr_t ff(resultat,false);
		char * ff1=(char *)ff;
		
		if (strcmp((const char*)ff1,"") != 0)
			value = atof((const char*) ff1)*100.0;
		
		theNode->Release();
		theNode=NULL;
		if (resultat) SysFreeString(resultat);
	}
	
	return value;
}

ARM_ReferenceValue* GetCRFStrike(MSXML2::IXMLDOMNode* node)
{
	MSXML2::IXMLDOMNode* theNode=NULL;
	MSXML2::IXMLDOMNode* theNode1=NULL;
	MSXML2::IXMLDOMNode* theNode2=NULL;
	MSXML2::IXMLDOMNode* theNode3=NULL;
	MSXML2::IXMLDOMNodeList* resultList = NULL;
	BSTR resultat = NULL;
	HRESULT hr;
	ARM_ReferenceValue* strike=NULL;

	VECTOR<double> strikeDate;
	VECTOR<double> strikeVal;

	node->selectSingleNode(_bstr_t("INTEREST_Rate"), &theNode);
	if (theNode!=NULL)
	{
		BSTR resultat = NULL;
		theNode->get_text(&resultat);

		_bstr_t ff(resultat,false);
		char * ff1=(char *)ff;

		double strike = atof((const char*) ff1);
		strikeVal.push_back(strike);
		strikeDate.push_back(GetStartDate(node).GetJulian());

		theNode->Release();
		theNode=NULL;
		if (resultat) SysFreeString(resultat);
	}
	
	// Recuperation des evenements							{
	if( node->selectNodes(_bstr_t("Events/EVENT"), &resultList)== S_OK)
	{
		if (resultList!=NULL)
		{
			long nb;
			resultList->get_length(&nb);
			for (long i=0 ; i<nb ; i++)
			{	
				hr=resultList->get_item(i, &theNode);
				if (hr==S_OK && theNode!=NULL)
				{
					// voir si le refvalue se cale sur les StartDate ou sur les EndDate...
					theNode->selectSingleNode(_bstr_t("Amount"), &theNode1);
					theNode->selectSingleNode(_bstr_t("Date"), &theNode3);
					theNode->selectSingleNode(_bstr_t("Type"), &theNode2);
					if (theNode2!=NULL)
					{
						theNode2->get_text(&resultat);
						_bstr_t ff(resultat,false);
						char * ff1=(char *)ff;
						if ((theNode3!=NULL) && (theNode1!=NULL) && (strcmp(ff1,"FRC")==0))
						{
							theNode1->get_text(&resultat);
							_bstr_t ff(resultat,false);
							char * ff1=(char *)ff;
							double strike= atof((const char*) ff1);
							strikeVal.push_back(strike);

							theNode3->get_text(&resultat);
							_bstr_t ff2(resultat,false);
							char * ff21=(char *)ff2;
							ARM_Date tmpDate (ff21,"YYYYMMDD");
							strikeDate.push_back(tmpDate.GetJulian());
						}

						if (resultat) SysFreeString(resultat);
						theNode2->Release();
						theNode2=NULL;
					}
					if (theNode1) theNode1->Release();
					theNode1=NULL;
					if (theNode3) theNode3->Release();
					theNode3=NULL;
				}
				if (theNode)
					theNode->Release();
			}
			if (resultList) resultList->Release();
			resultList=NULL;
		}
	}
	ARM_Vector* vstrikeVal= CreateARMVectorFromVECTOR(strikeVal);
	ARM_Vector* vstrikeDate= CreateARMVectorFromVECTOR(strikeDate);
	
	if (vstrikeDate->GetSize() == 1)
	{
		strike = new ARM_ReferenceValue(vstrikeVal->Elt(0));

		if (vstrikeVal)
			delete vstrikeVal;
		vstrikeVal = NULL;

		if (vstrikeDate)
			delete vstrikeDate;
		vstrikeDate = NULL;
	}
	else
	{
		strike = new ARM_ReferenceValue(vstrikeDate,vstrikeVal);
		strike->SetCalcMethod(K_STEPUP_LEFT);
	}

	return strike;
}

int GetCapFloor(MSXML2::IXMLDOMNode* node)
{
	MSXML2::IXMLDOMNode* theNode=NULL;
	BSTR resultat = NULL;
	int value=K_CALL;

	node->selectSingleNode(_bstr_t("SubType"), &theNode);
	if (theNode!=NULL)
	{
		theNode->get_text(&resultat);
		_bstr_t ff(resultat,false);
		char * ff1=(char *)ff;
		if (strcmp(ff1,"FLOOR")==0)
			value=K_PUT;
		else
			value=K_CALL;
		
		theNode->Release();
		theNode=NULL;
		if (resultat) SysFreeString(resultat);
	}
	return value;
}

ARM_ReferenceValue* GetCapStrikes(MSXML2::IXMLDOMNode* node)
{
	MSXML2::IXMLDOMNode* theNode=NULL;
	MSXML2::IXMLDOMNode* theNode1=NULL;
	MSXML2::IXMLDOMNode* theNode2=NULL;
	MSXML2::IXMLDOMNode* theNode3=NULL;
	MSXML2::IXMLDOMNodeList* resultList = NULL;
	BSTR resultat = NULL;
	HRESULT hr;
	ARM_ReferenceValue* strike=NULL;

	VECTOR<double> strikeDate;
	VECTOR<double> strikeVal;

	strikeVal.push_back(GetCapStrike(node));
	strikeDate.push_back(GetStartDate(node).GetJulian());
	
	// Recuperation des evenements							{
	if( node->selectNodes(_bstr_t("Events/EVENT"), &resultList)== S_OK)
	{
		if (resultList!=NULL)
		{
			long nb;
			resultList->get_length(&nb);
			for (long i=0 ; i<nb ; i++)
			{	
				hr=resultList->get_item(i, &theNode);
				if (hr==S_OK && theNode!=NULL)
				{
					// voir si le refvalue se cale sur les StartDate ou sur les EndDate...
					theNode->selectSingleNode(_bstr_t("Amount"), &theNode1);
					theNode->selectSingleNode(_bstr_t("Date"), &theNode3);
					theNode->selectSingleNode(_bstr_t("Type"), &theNode2);
					if (theNode2!=NULL)
					{
						theNode2->get_text(&resultat);
						_bstr_t ff(resultat,false);
						char * ff1=(char *)ff;
						if ((theNode3!=NULL) && (theNode1!=NULL) && (strcmp(ff1,"STR")==0))
						{
							theNode3->get_text(&resultat);
							_bstr_t ff2(resultat,false);
							char * ff21=(char *)ff2;
							ARM_Date tmpDate (ff21,"YYYYMMDD");

							// pour eviter d'éventuelles erreurs de saisie dans Summit,
							// chaque date d'évènement récupérée doit être postérieure à la date précédente !
							if (tmpDate.GetJulian() > strikeDate[strikeDate.size()-1])
							{
								strikeDate.push_back(tmpDate.GetJulian());

								theNode1->get_text(&resultat);
								_bstr_t ff(resultat,false);
								char * ff1=(char *)ff;
								double strike= atof((const char*) ff1)*100.0;
								strikeVal.push_back(strike);
							}
						}

						if (resultat) SysFreeString(resultat);
						theNode2->Release();
						theNode2=NULL;	
					}
					if (theNode1) theNode1->Release();
					theNode1=NULL;
					if (theNode3) theNode3->Release();
					theNode3=NULL;
				}
				if (theNode)
					theNode->Release();
			}
			if (resultList) resultList->Release();
			resultList=NULL;
		}
	}

	ARM_Vector* vstrikeVal= CreateARMVectorFromVECTOR(strikeVal);
	ARM_Vector* vstrikeDate= CreateARMVectorFromVECTOR(strikeDate);

	if (vstrikeDate->GetSize() == 1)
	{
		strike = new ARM_ReferenceValue(vstrikeVal->Elt(0));

		if (vstrikeVal)
			delete vstrikeVal;
		vstrikeVal = NULL;
		
		if (vstrikeDate)
			delete vstrikeDate;
		vstrikeDate = NULL;
	}
	else
	{
		strike = new ARM_ReferenceValue(vstrikeDate,vstrikeVal);

		// interpol par STEP
		strike->SetCalcMethod(K_STEPUP_LEFT);
	}

	return strike;

}

double GetCapStrike(MSXML2::IXMLDOMNode* node)
{
	MSXML2::IXMLDOMNode* theNode=NULL;
	MSXML2::IXMLDOMNodeList* resultList = NULL;
	BSTR resultat = NULL;
	double strike;

	node->selectSingleNode(_bstr_t("INTEREST_Rate"), &theNode);
	if (theNode!=NULL)
	{
		BSTR resultat = NULL;
		theNode->get_text(&resultat);

		_bstr_t ff(resultat,false);
		char * ff1=(char *)ff;

		strike = atof((const char*) ff1)*100.0;

		theNode->Release();
		theNode=NULL;
		if (resultat) SysFreeString(resultat);
	}
	return strike;
}

void GetResetRollDate(MSXML2::IXMLDOMNode* node, char* resetRollDate)
{	
	MSXML2::IXMLDOMNode* theNode=NULL;
	BSTR resultat = NULL;

	node->selectSingleNode(_bstr_t("SCHED_Reset_AnnDay"), &theNode);
	if (theNode!=NULL)
	{
		BSTR resultat = NULL;
		theNode->get_text(&resultat);

		_bstr_t ff(resultat,false);
		char * ff1=(char *)ff;
		strcpy(resetRollDate,ff1);
		
		theNode->Release();
		theNode=NULL;
	}
	if (resultat) SysFreeString(resultat);	
}

void GetPayRollDate(MSXML2::IXMLDOMNode* node, char* payRollDate)
{	
	MSXML2::IXMLDOMNode* theNode=NULL;
	BSTR resultat = NULL;

	node->selectSingleNode(_bstr_t("SCHED_Pay_AnnDay"), &theNode);
	if (theNode!=NULL)
	{
		BSTR resultat = NULL;
		theNode->get_text(&resultat);

		_bstr_t ff(resultat,false);
		char * ff1=(char *)ff;
		strcpy(payRollDate,ff1);
		
		theNode->Release();
		theNode=NULL;
	}
	if (resultat) SysFreeString(resultat);
}

int GetDecompFlag(MSXML2::IXMLDOMNode* node)
{
	MSXML2::IXMLDOMNode* theNode=NULL;
	BSTR resultat = NULL;
	int value= 1;

	node->selectSingleNode(_bstr_t("INTEREST_FwdDecomp"), &theNode);
	if (theNode!=NULL)
	{
		theNode->get_text(&resultat);
		_bstr_t ff(resultat,false);
		char * ff1=(char *)ff;

		if (strcmp((const char*) ff1,"N") == 0)
			value = 0;
		if (resultat) SysFreeString(resultat);
		theNode->Release();
		theNode=NULL;
	}
	return value;
}


void GetStubDate1(MSXML2::IXMLDOMNode* node, char* stubdate)
{
	MSXML2::IXMLDOMNode* theNode=NULL;
	BSTR resultat = NULL;

	node->selectSingleNode(_bstr_t("STUB_Date1"), &theNode);
	if (theNode!=NULL)
	{
		BSTR resultat = NULL;
		theNode->get_text(&resultat);

		_bstr_t ff(resultat,false);
		char * ff1=(char *)ff;

		if (strcmp((const char*)ff1,"") == 0)
			strcpy(stubdate,"NULL");
		else
		{
			ARM_Date tmpDate (ff1,"YYYYMMDD");
			tmpDate.JulianToStrDate(stubdate);
		}

		theNode->Release();
		theNode=NULL;
		if (resultat) SysFreeString(resultat);
	}
	else
		strcpy(stubdate,"NULL");
}


void GetStubDate2(MSXML2::IXMLDOMNode* node, char* stubdate)
{
	MSXML2::IXMLDOMNode* theNode=NULL;
	BSTR resultat = NULL;

	node->selectSingleNode(_bstr_t("STUB_Date2"), &theNode);
	if (theNode != NULL)
	{
		BSTR resultat = NULL;
		theNode->get_text(&resultat);

		_bstr_t ff(resultat,false);
		char * ff1=(char *)ff;

		if (strcmp((const char*)ff1,"") == 0)
			strcpy(stubdate,"NULL");
		else
		{
			ARM_Date tmpDate (ff1,"YYYYMMDD");
			tmpDate.JulianToStrDate(stubdate);
		}

		theNode->Release();
		theNode=NULL;
		if (resultat) SysFreeString(resultat);
	}
	else
		strcpy(stubdate,"NULL");
}


void GetCustom(MSXML2::IXMLDOMNode* node, char* cust)
{
	MSXML2::IXMLDOMNode* theNode=NULL;
	BSTR resultat = NULL;

	node->selectSingleNode(_bstr_t("Custom"), &theNode);
	if (theNode!=NULL)
	{
		BSTR resultat = NULL;
		theNode->get_text(&resultat);

		_bstr_t ff(resultat,false);
		char * ff1=(char *)ff;

		if (strcmp((const char*)ff1,"") == 0)
			strcpy(cust,"STND");
		else
			strcpy(cust,ff1);

		theNode->Release();
		theNode=NULL;
		if (resultat) SysFreeString(resultat);
	}
	else
		strcpy(cust,"STND");
}


//***************************************************************************************************
// Appel au Niveau CONFIG_VECTOR
//***************************************************************************************************


bool GetConfigVectorFlag(MSXML2::IXMLDOMNode* node, int pos)
{
	BSTR resultat = NULL;
	int size;
	bool flag = false;

	try
	{
		if (node != NULL)
		{
			node->get_text(&resultat);
			_bstr_t ff(resultat,false);
			char * ff1=(char *)ff;

			size = strlen((const char*)ff1);

			if ((0 <= pos) && (pos < size))
			{
				flag = ('1' == ff1[pos]);
			}

			if (resultat) SysFreeString(resultat);
		}
	}
	catch(...)
	{
	}

	return flag;
}


int ParseCapCashFlows(MSXML2::IXMLDOMNode* node, 
						ARM_Vector * flowStartDates, 
						ARM_Vector * flowEndDates,
						ARM_Vector * paymentDates, 
						ARM_Vector * resetDates, 
						ARM_Vector * intDays, 
						ARM_Vector * strike, 
						ARM_Vector * fwd,
						ARM_ReferenceValue * spread, 
						ARM_ReferenceValue * Notional)
{
	MSXML2::IXMLDOMNode* theNode=NULL;
	MSXML2::IXMLDOMNode* theNode1=NULL;
	MSXML2::IXMLDOMNode* theNode2=NULL;
	MSXML2::IXMLDOMNode* theNode3=NULL;
	MSXML2::IXMLDOMNodeList* resultList = NULL;
	HRESULT hr;
	BSTR resultat = NULL;
	
	VECTOR<double> vstartDate;
	VECTOR<double> vendDate;
	VECTOR<double> vresetDate;
	VECTOR<double> vpayDate;
	VECTOR<double> vintDays;
	VECTOR<double> vnominal;
	VECTOR<double> vstrike;
	VECTOR<double> vspread;
	VECTOR<double> vfwd;

	// Recuperation des Flows							
	if( node->selectNodes(_bstr_t("Flow"), &resultList)== S_OK)
	{
		if (resultList!=NULL)
		{
			long nb;
			resultList->get_length(&nb);
			for (long i=0 ; i<nb ; i++)
			{	
				hr=resultList->get_item(i, &theNode);
				if (hr==S_OK && theNode!=NULL)
				{
					theNode->selectSingleNode(_bstr_t("Type"), &theNode2);
					if (theNode2!=NULL)
					{
						theNode2->get_text(&resultat);
						_bstr_t ff(resultat,false);
						char * ff1=(char *)ff;
						if ((strcmp(ff1,"IPR")==0) || (strcmp(ff1,"INT")==0) || (strcmp(ff1,"")==0))
						{
							vstartDate.push_back(GetDateFromXMLNode(theNode,"StartDate").GetJulian());

							vendDate.push_back(GetDateFromXMLNode(theNode,"EndDate").GetJulian());

							double reset= GetDateFromXMLNode(theNode,"FixingDate").GetJulian();
							if ( reset == ARM_Date(Summit_DefaultDate,"YYYYMMDD").GetJulian() )
								vresetDate.push_back(vresetDate.back());
							else
								vresetDate.push_back(GetDateFromXMLNode(theNode,"FixingDate").GetJulian());

							theNode->selectSingleNode(_bstr_t("IntDays"), &theNode3);
							if (theNode3!=NULL)
							{
								theNode3->get_text(&resultat);
								_bstr_t ff(resultat,false);
								char * ff1=(char *)ff;
								double tmp= atof((const char*) ff1);
								vintDays.push_back(tmp);
								theNode3->Release();
								theNode3=NULL;
							}
							theNode->selectSingleNode(_bstr_t("Spread"), &theNode3);
							if (theNode3!=NULL)
							{
								theNode3->get_text(&resultat);
								_bstr_t ff(resultat,false);
								char * ff1=(char *)ff;
								double tmp= atof((const char*) ff1);
								vspread.push_back(tmp/100.0);
								theNode3->Release();
								theNode3=NULL;
							}

							vpayDate.push_back(GetDateFromXMLNode(theNode,"PayDate").GetJulian());

							theNode->selectSingleNode(_bstr_t("Strike"), &theNode3);
							if (theNode3!=NULL)
							{
								theNode3->get_text(&resultat);
								_bstr_t ff(resultat,false);
								char * ff1=(char *)ff;
								double tmp= atof((const char*) ff1);
								vstrike.push_back(tmp);
								theNode3->Release();
								theNode3=NULL;
							}
							theNode->selectSingleNode(_bstr_t("Notional"), &theNode3);
							if (theNode3!=NULL)
							{
								theNode3->get_text(&resultat);
								_bstr_t ff(resultat,false);
								char * ff1=(char *)ff;
								double tmp= atof((const char*) ff1);
								vnominal.push_back(tmp);
								theNode3->Release();
								theNode3=NULL;
							}
							theNode->selectSingleNode(_bstr_t("Forward"), &theNode3);
							if (theNode3!=NULL)
							{
								theNode3->get_text(&resultat);
								_bstr_t ff(resultat,false);
								char * ff1=(char *)ff;
								double tmp= atof((const char*) ff1);
								vfwd.push_back(tmp);
								theNode3->Release();
								theNode3=NULL;
							}

						}

						if (resultat) SysFreeString(resultat);
						theNode2->Release();
						theNode2=NULL;	
					}
					
					theNode->Release();
					theNode=NULL;
					
				}
			}
			if (resultList) resultList->Release();
			resultList=NULL;
		}
	}
	*flowStartDates = *CreateARMVectorFromVECTOR(vstartDate);
	*flowEndDates= *CreateARMVectorFromVECTOR(vendDate);
	*resetDates= *CreateARMVectorFromVECTOR(vresetDate);
	*paymentDates= *CreateARMVectorFromVECTOR(vpayDate);
	*intDays= *CreateARMVectorFromVECTOR(vintDays);
	*strike= *CreateARMVectorFromVECTOR(vstrike);
	*fwd= *CreateARMVectorFromVECTOR(vfwd);
	
	// Traitement des dates de paiement pour resetFreq > paimentFreq : Summit met 31/12/9999 
	
	ARM_Date defaultDate(31,12,9999);
	double dDefaultDate = defaultDate.GetJulian();

	for (int i = 0; i < paymentDates->GetSize(); i++)
	{
		if (paymentDates->Elt(i) == dDefaultDate)
		{
			int j = i;
			while (j < paymentDates->GetSize())
			{
				if (paymentDates->Elt(j) != dDefaultDate)
				{
					paymentDates->Elt(i) = paymentDates->Elt(j);
					j = paymentDates->GetSize();
				}
				j++;
			}
		}
	}
	
	*spread= ARM_ReferenceValue((ARM_Vector*)paymentDates->Clone(), (ARM_Vector*)CreateARMVectorFromVECTOR(vspread)->Clone());
	spread->SetCalcMethod(K_STEPUP_LEFT);
	*Notional = ARM_ReferenceValue((ARM_Vector*)paymentDates->Clone(), (ARM_Vector*)CreateARMVectorFromVECTOR(vnominal)->Clone());
	Notional->SetCalcMethod(K_STEPUP_RIGHT);
	Notional->SetExtrapolMeth(1);
	
	return 1;

}

ARM_ReferenceValue* GetPrimes(MSXML2::IXMLDOMNode* node)
{
	MSXML2::IXMLDOMNode* theNode=NULL;
	MSXML2::IXMLDOMNode* listItem=NULL;
	MSXML2::IXMLDOMNodeList* resultList = NULL;
	ARM_ReferenceValue* newRefVal = NULL;
	BSTR resultat = NULL;
	HRESULT hr;
	long nbEvents;
	
	VECTOR<double> PrimeDate;
	VECTOR<double> PrimeVal;
	VECTOR<double> PrimeCcy;
	
	node->selectNodes(_bstr_t("Events/EVENT"), &resultList);
	if (resultList!=NULL)
	{
		resultList->get_length(&nbEvents);
		CCString type;

		for (int i = 0; i < nbEvents; i++)
		{
			hr=resultList->get_item(i, &listItem);
			listItem->selectSingleNode(_bstr_t("Type"), &theNode);
			if (theNode!=NULL)
			{
				BSTR resultat = NULL;
				theNode->get_text(&resultat);

				_bstr_t ff(resultat,false);
				char * ff1=(char *)ff;

				type = (const char*) ff1;

				theNode->Release();
				theNode=NULL;
				if (resultat) SysFreeString(resultat);
			}

			if ((strcmp(type,"PRM") == 0) || (strcmp(type,"F/M") == 0)
				|| (strcmp(type,"F/E") == 0) || (strcmp(type,"ADJ") == 0)
				|| (strcmp(type,"FEE") == 0) || (strcmp(type,"F/O") == 0)
				|| (strcmp(type,"PLA") == 0)
				)
			{
				PrimeDate.push_back(GetDateFromXMLNode(listItem,"Date").GetJulian());

				listItem->selectSingleNode(_bstr_t("Amount"), &theNode);
				if (theNode!=NULL)
				{
					BSTR resultat = NULL;
					theNode->get_text(&resultat);
					_bstr_t ff(resultat,false);
					char * ff1=(char *)ff;

					double amount = atof((const char*) ff1);
					PrimeVal.push_back(amount);

					if (resultat) SysFreeString(resultat);
					theNode->Release();
					theNode=NULL;
				}

				listItem->selectSingleNode(_bstr_t("Ccy"), &theNode);
				if (theNode!=NULL)
				{
					BSTR resultat = NULL;
					theNode->get_text(&resultat);
					_bstr_t ff(resultat,false);
					char * ff1=(char *)ff;

					PrimeCcy.push_back((double)ARM_GetCountrySymbol(ff1));

					if (resultat) SysFreeString(resultat);
					theNode->Release();
					theNode=NULL;
				}
			}
			if (listItem) listItem->Release();
			listItem=NULL;
		}
		if (resultList) resultList->Release();
		resultList=NULL;
		
	}
	
	if (PrimeDate.size() > 0)
	{
		ARM_Vector* vPrimeDate = CreateARMVectorFromVECTOR(PrimeDate);
		ARM_Vector* vPrimeVal = CreateARMVectorFromVECTOR(PrimeVal);
		ARM_Vector* vPrimeCcy = CreateARMVectorFromVECTOR(PrimeCcy);

		newRefVal= new ARM_ReferenceValue((ARM_Vector*)vPrimeDate->Clone(),(ARM_Vector*)vPrimeVal->Clone(),(ARM_Vector*)vPrimeCcy->Clone());
		newRefVal->SetCalcMethod(K_STEPUP_RIGHT);
		newRefVal->SetExtrapolMeth(1);
		if (vPrimeDate)
			delete vPrimeDate;
		vPrimeDate = NULL;

		if (vPrimeVal)
			delete vPrimeVal;
		vPrimeVal = NULL;

		if (vPrimeCcy)
			delete vPrimeCcy;
		vPrimeCcy = NULL;
	}
	
	return newRefVal;

}

ARM_Vector * GetFixing(ARM_Vector * fwd, ARM_Vector * fixingDates, ARM_Date date)
{
	VECTOR<double> Fixing;
	
	int nb = fwd->GetSize();
	for (long i=0 ; i<nb ; i++)
	{
		if (fixingDates->Elt(i) <= date.GetJulian())
			Fixing.push_back(fwd->Elt(i));
	}

	return CreateARMVectorFromVECTOR(Fixing);
}




char* filterBlank( const char* input, char c = ' ' )
{
 char* result = new char[strlen(input)+1];
 char* tmp = result;
 while( *input )
 {
  if( *input != c )
   *tmp++ = *input;
  ++input;
 }
 *tmp ='\0';
 return result;
}


void GetFormulaInMemorySO(MSXML2::IXMLDOMNode* node, char* formula)
{
	//SPREADOPTIONLOGFLTDIGITAL ( [2yEUReurib], [10yEUReurib], [10yEUReurib] )
	
	MSXML2::IXMLDOMNode* theNode=NULL;
	BSTR resultat = NULL;

	node->selectSingleNode(_bstr_t("Formula/FORMULA/Formula"), &theNode);
	if (theNode!=NULL)
	{
		BSTR resultat = NULL;
		theNode->get_text(&resultat);

		_bstr_t ff(resultat,false);
		char * ff1=(char *)ff;

		// on supprime les espaces dans la formule
		char* tmp = NULL;
		char* tmptmp = NULL;
		tmp = filterBlank(ff1);
		tmptmp = filterBlank(tmp,'"');
		strcpy(ff1,tmptmp);
		if (tmp)
			delete [] tmp;
		if (tmptmp)
			delete [] tmptmp;

		int c ='(';
		char* tmp2= strchr(ff1, c);
		char* tmp3= strchr(tmp2+1, c);
		if (tmp3!=NULL)
			strncpy(formula,tmp2+1,(int) (tmp3-tmp2-1));
	
		formula[(int) (tmp3-tmp2-1)] = '\0';

		theNode->Release();
		theNode=NULL;
	}
	if (resultat) SysFreeString(resultat);
	strcpy(formula, strupr(formula));
}


void GetFormula(MSXML2::IXMLDOMNode* node, char* formula)
{
	//SPREADOPTIONLOGFLTDIGITAL ( [2yEUReurib], [10yEUReurib], [10yEUReurib] )
	
	MSXML2::IXMLDOMNode* theNode=NULL;
	BSTR resultat = NULL;

	node->selectSingleNode(_bstr_t("Formula/FORMULA/Formula"), &theNode);
	if (theNode!=NULL)
	{
		BSTR resultat = NULL;
		theNode->get_text(&resultat);

		_bstr_t ff(resultat,false);
		char * ff1=(char *)ff;

		// on supprime les espaces dans la formule
		char* tmp = NULL;
		tmp = filterBlank(ff1);
		strcpy(ff1,tmp);
		if (tmp)
			delete [] tmp;

		int c ='(';
		char* tmp2= strchr(ff1, c);
		if (tmp2!=NULL)
			strncpy(formula,ff1,(int) (tmp2-ff1));
	
		formula[(int) (tmp2-ff1)] = '\0';

		theNode->Release();
		theNode=NULL;
	}
	if (resultat) SysFreeString(resultat);
	strcpy(formula, strupr(formula));
}


void GetNthArgInFormula(MSXML2::IXMLDOMNode* node, int n, char* Arg, int mode)
{
	MSXML2::IXMLDOMNode* theNode=NULL;
	BSTR resultat = NULL;

	if (mode == 0)
		node->selectSingleNode(_bstr_t("Formula/FORMULA/Formula"), &theNode);
	else
		node->selectSingleNode(_bstr_t("Formula"), &theNode);
	if (theNode)
	{
		BSTR resultat = NULL;
		theNode->get_text(&resultat);

		_bstr_t ff(resultat,false);
		char * ff1=(char *)ff;

		// on supprime les espaces dans la formule
		char* tmp = NULL;
		tmp = filterBlank(ff1);
		strcpy(ff1,tmp);
		if (tmp)
			delete [] tmp;

		// les arguments commencent apres '('
		// pour les Memory, il peut y en avoir plusieurs
		ff1 =  strchr(ff1, '(');
		char* index = strchr(ff1+1, '(');

		if (index != NULL)
			ff1 = index;

		// le separateur d'argument
		int c =',';
		for (int i= 1; (i<n) && (ff1); i++)
		{
			if (ff1)
				ff1= strchr(ff1+1, c);
		}
		if (ff1)
		{
			char * end= strchr(ff1+2, c);
			if (end)
			{
				strncpy(Arg,ff1+1,(int) (end-ff1-1));
				Arg[(int) (end-ff1-1)] = '\0';
			}
			else // c'est le dernier argument
			{
				char * end= strchr(ff1+1, ')');
				if (end)
				{
					strncpy(Arg,ff1+1,(int) (end-ff1-1));
					Arg[(int) (end-ff1-1)] = '\0';
				}
			}
		}
		else
		{
			strcpy(Arg,"");
		}

		theNode->Release();
		theNode=NULL;
	}
	if (resultat) SysFreeString(resultat);
}

ARM_ReferenceValue* GetRefValueFromFormula(MSXML2::IXMLDOMNode* node, int nthArg)
{
	MSXML2::IXMLDOMNode* theNode=NULL;
	MSXML2::IXMLDOMNodeList* listItem=NULL;
	BSTR resultat = NULL;
	long nbFormulas = 0;
	ARM_Vector* dates = NULL;
	ARM_Vector* values = NULL;
	ARM_ReferenceValue* refValue = NULL;
	bool isValid = false;

	node->selectNodes(_bstr_t("Formula/FORMULA"), &listItem);
	if (listItem!=NULL)
	{
		listItem->get_length(&nbFormulas);
		dates = new ARM_Vector(nbFormulas);
		values = new ARM_Vector(nbFormulas);

		for (int i = 0; i < nbFormulas; i++)
		{
			listItem->get_item(i, &theNode);
			dates->Elt(i) = (GetDateFromXMLNode(theNode,"Date")).GetJulian();

			char buffer[15];
			GetNthArgInFormula(theNode,nthArg,buffer,1);
			if (strcmp(buffer, "") != 0)
				isValid = true;

			values->Elt(i) = atof(buffer);

			if (theNode)
				theNode->Release();
		}
		listItem->Release();

		if (!isValid) //no value found, return NULL
		{
			delete dates;
			delete values;
		}
		else
		{
			if (nbFormulas == 1)
			{
				refValue = new ARM_ReferenceValue(values->Elt(0));
			}
			else
			{
				refValue = new ARM_ReferenceValue(dates, values);
				refValue->SetCalcMethod(K_STEPUP_LEFT);
			}
		}
	}
	return refValue;
}

/// In formula "FORM ( Arg1[, Arg2...] )", returns Arg1 that represents Index
int GetSpreadIndexType1(MSXML2::IXMLDOMNode* node, ARM_Currency* ccy, int isCMT)
{
	char summitIndex[50];
	char armIndex[50];
	
	GetNthArgInFormula(node,1,summitIndex);

	if (strcmp(summitIndex,"")==0)
		return K_FIXED;
	else
	{
		ConvSummitIndexToARMIndex(summitIndex,ccy,armIndex, isCMT);

		return ARM_ConvIrType(armIndex);
	}
}

/// In formula "FORM ( Arg1, Arg2... )", returns Arg2
int GetSpreadIndexType2(MSXML2::IXMLDOMNode* node, ARM_Currency* ccy)
{
	char summitIndex[50];
	char armIndex[50];
	
	GetNthArgInFormula(node,2,summitIndex);

	if (strcmp(summitIndex,"")==0)
		return K_FIXED;
	else
	{
		ConvSummitIndexToARMIndex(summitIndex,ccy,armIndex);

		return ARM_ConvIrType(armIndex);
	}
}

/// In formula "FORM1 ( FORM2 (Arg1) )"
/// returns index type represented by FORM2
double GetSubFormulaIndexType(MSXML2::IXMLDOMNode* node, ARM_Currency* ccy)
{
	BSTR resultat = NULL;
	node->get_text(&resultat);

	_bstr_t ff(resultat,false);
	char* ff1 = filterBlank((char*)ff);
	char* subform = strchr(ff1+1, '(');
	char* index = strchr(subform+1, '(');
	char* Argument = "";

	if (index)
	{
		char* end = strchr(index+1, ')');
		if (end)
		{
			strncpy(Argument, index+1, (int)(end-index-1));
			Argument[(int)(end-index-1)] = '\0';
		}
	}
	else
	{
		strcpy(Argument,"");
	}

	char armIndex[15];
	if (strcmp(Argument,"")==0)
		return K_FIXED;
	else
	{
		ConvSummitIndexToARMIndex(Argument, ccy, armIndex);
		return ARM_ConvIrType(armIndex);
	}
}

///	Returns: 
///    K_LIVRET_A if livretA
///    -1 if not quanto
///    index type otherwise
///	Action : find quanto index in formula
///    formula can be like : 
///    - QUANTO ( Index ( [TermCcyIbor] ) )
///    - QUANTOG ( [TermCcyIbor] )
///    - QUANTOG2 ( [TermCcyIbor] )
int GetFormulaIndexType(MSXML2::IXMLDOMNode* node, ARM_Currency* ccy, int mode)
{
	MSXML2::IXMLDOMNode* theNode=NULL;
	BSTR resultat = NULL;

	try
	{
		if (mode == 0)
			node->selectSingleNode(_bstr_t("Formula/FORMULA/Formula"), &theNode);
		else
			node->selectSingleNode(_bstr_t("Formula"), &theNode);

		if (theNode)
		{
			BSTR resultat = NULL;
			theNode->get_text(&resultat);

			_bstr_t ff(resultat,false);
			char* ff1 = filterBlank((char*)ff);
			char* index = strchr(ff1, '(');
			char summitIndex[15];
			char armIndex[15];
			if (index)
			{
				if (strncmp(ff1, "LIVRET_A", 8) == 0)
				{
					if (resultat) SysFreeString(resultat);
					theNode->Release();
					theNode = NULL;

					return K_LIVRET_A;
				}
				else if ( (strncmp(ff1, "QUANTOG2", 8) == 0) || (strncmp(ff1, "QUANTOG", 7) == 0) )
				{
					GetNthArgInFormula(node,1,summitIndex);
					ConvSummitIndexToARMIndex(summitIndex, ccy, armIndex);

					if (resultat) SysFreeString(resultat);
					theNode->Release();
					theNode = NULL;

					return ARM_ConvIrType(armIndex);
				}
				else
				{
					// Formula can be : FORM(QUANTOG())
					char* subFormula = strchr(index+1, '(');
					if ( (subFormula) && ((strncmp(index+1, "QUANTOG2", 8) == 0) || (strncmp(index+1, "QUANTOG", 7) == 0)) )
					{
						int res = GetSubFormulaIndexType(theNode, ccy);

						if (resultat) SysFreeString(resultat);
						theNode->Release();
						theNode = NULL;

						return res;
					}
				}
				if (strncmp(ff1, "QUANTO", 6)==0)
				{
					int res = GetSubFormulaIndexType(theNode, ccy);

					if (resultat) SysFreeString(resultat);
					theNode->Release();
					theNode = NULL;

					return res;
				}
			}
			if (resultat) SysFreeString(resultat);
			theNode->Release();
			theNode = NULL;
		}

		return -1;
	}
	catch (...)
	{
		CCString msg((CCString)"Error in GetFormulaIndexType() ");

		throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT, (char*) msg);
	}
}
		
double GetSpreadPayOff(MSXML2::IXMLDOMNode* node)
{
	char buffer[15];
	
	GetNthArgInFormula(node,3,buffer);
	return atof(buffer);
}

double GetSpreadSlopeFlag(MSXML2::IXMLDOMNode* node)
{
	char buffer[15];
	
	GetNthArgInFormula(node,6,buffer);
	if(strcmp(buffer,"")==0)
		return 1.0;
	else
		return atof(buffer);
}

ARM_ReferenceValue* GetVariableSpreadPayOff(MSXML2::IXMLDOMNode* node)
{
	MSXML2::IXMLDOMNode* theNode=NULL;
	MSXML2::IXMLDOMNodeList* listItem=NULL;
	BSTR resultat = NULL;
	long nbFormulas;
	ARM_ReferenceValue* newRefVal = NULL;
	ARM_Vector* vDates = NULL;
	ARM_Vector* vVal = NULL;

	node->selectNodes(_bstr_t("Formula/FORMULA"), &listItem);
	if (listItem!=NULL)
	{
		listItem->get_length(&nbFormulas);
		vDates = new ARM_Vector(nbFormulas);
		vVal = new ARM_Vector(nbFormulas);

		for (int i = 0; i < nbFormulas; i++)
		{
			listItem->get_item(i, &theNode);
			vDates->Elt(i) = (GetDateFromXMLNode(theNode,"Date")).GetJulian();
			char buffer[15];

			GetNthArgInFormula(theNode,3,buffer,1);
			vVal->Elt(i) = atof(buffer);

			if (theNode)
				theNode->Release();
		}
		listItem->Release();

		if (nbFormulas == 1)
		{
			newRefVal = new ARM_ReferenceValue(vVal->Elt(0));
			if (vVal)
				delete vVal;
			if (vDates)
				delete vDates;
		}
		else
		{
			newRefVal = new ARM_ReferenceValue(vDates, vVal, K_YIELD, 1);
			newRefVal->SetCalcMethod(K_STEPUP_LEFT);
		}
	}
	return newRefVal;
}

ARM_ReferenceValue* GetLowBarrier(MSXML2::IXMLDOMNode* node, int resetGap, char* resetCal)
{
	MSXML2::IXMLDOMNode* theNode			= NULL;
	MSXML2::IXMLDOMNodeList* listItem		= NULL;
	BSTR resultat					= NULL;
	ARM_ReferenceValue* newRefVal	= NULL;
	ARM_Vector* vDates				= NULL;
	ARM_Vector* vVal				= NULL;
	long nbFormulas;
	double dCurrDate;
	ARM_Date currDate;

	node->selectNodes(_bstr_t("Formula/FORMULA"), &listItem);
	if (listItem != NULL)
	{
		listItem->get_length(&nbFormulas);
		vDates	= new ARM_Vector(nbFormulas);
		vVal	= new ARM_Vector(nbFormulas);

		for (int i = 0; i < nbFormulas; i++)
		{
			listItem->get_item(i, &theNode);
			dCurrDate	= (GetDateFromXMLNode(theNode,"Date")).GetJulian();
			currDate	= ARM_Date(dCurrDate);
			//we adjust the date on the first ref fixing date.
			currDate.GapBusinessDay(resetGap, resetCal);
			vDates->Elt(i) = currDate.GetJulian();
			
			char buffer[15];
			GetNthArgInFormula(theNode,2,buffer,1);
			vVal->Elt(i) = atof(buffer);

			if (theNode)
				theNode->Release();
		}
		listItem->Release();

		if (nbFormulas == 1)
		{
			newRefVal = new ARM_ReferenceValue(vVal->Elt(0));
		}
		else
		{
			newRefVal = new ARM_ReferenceValue(vDates, vVal, K_YIELD, 1);
			newRefVal->SetCalcMethod(K_STEPUP_LEFT);
		}
	}
	return newRefVal;
}

ARM_ReferenceValue* GetUpBarrier(MSXML2::IXMLDOMNode* node, int resetGap, char* resetCal)
{
	MSXML2::IXMLDOMNode* theNode			= NULL;
	MSXML2::IXMLDOMNodeList* listItem		= NULL;
	BSTR resultat					= NULL;
	ARM_ReferenceValue* newRefVal	= NULL;
	ARM_Vector* vDates				= NULL;
	ARM_Vector* vVal				= NULL;
	
	long nbFormulas;
	double dCurrDate;
	ARM_Date currDate;

	node->selectNodes(_bstr_t("Formula/FORMULA"), &listItem);

	if (listItem != NULL)
	{
		listItem->get_length(&nbFormulas);
		vDates	= new ARM_Vector(nbFormulas);
		vVal	= new ARM_Vector(nbFormulas);

		for (int i = 0; i < nbFormulas; i++)
		{
			listItem->get_item(i, &theNode);
			dCurrDate = (GetDateFromXMLNode(theNode,"Date")).GetJulian();
			currDate = ARM_Date(dCurrDate);
			//we adjust the date on the first ref fixing date.
			currDate.GapBusinessDay(resetGap, resetCal);
			vDates->Elt(i) = currDate.GetJulian();

			char buffer[15];
			GetNthArgInFormula(theNode,3,buffer,1);
			vVal->Elt(i) = atof(buffer);

			if (theNode)
				theNode->Release();
		}
		listItem->Release();

		if (nbFormulas == 1)
		{
			newRefVal = new ARM_ReferenceValue(vVal->Elt(0));
		}
		else
		{
			newRefVal = new ARM_ReferenceValue(vDates, vVal, K_YIELD, 1);
			newRefVal->SetCalcMethod(K_STEPUP_LEFT);
		}
	}
	return newRefVal;
}

ARM_ReferenceValue* GetBoostedRate(MSXML2::IXMLDOMNode* node)
{
	MSXML2::IXMLDOMNode* theNode			= NULL;
	MSXML2::IXMLDOMNodeList* listItem		= NULL;
	BSTR resultat					= NULL;
	long nbFormulas					= 0;
	ARM_ReferenceValue* newRefVal	= NULL;
	ARM_Vector* vDates				= NULL;
	ARM_Vector* vVal				= NULL;
	int check						= 0.0;
	double oneBoostedRate			= 0.0;

	node->selectNodes(_bstr_t("Formula/FORMULA"), &listItem);
	if (listItem!=NULL)
	{
		listItem->get_length(&nbFormulas);
		vDates = new ARM_Vector(nbFormulas);
		vVal = new ARM_Vector(nbFormulas);

		for (int i = 0; i < nbFormulas; i++)
		{
			listItem->get_item(i, &theNode);
			vDates->Elt(i) = (GetDateFromXMLNode(theNode,"Date")).GetJulian(); 
			char buffer[15];
			GetNthArgInFormula(theNode,5,buffer,1);
			check = atof(buffer);

			//If the boosted rate does not appear: it is not step-up and then unique.   
			if (check == 0)
			{
				MSXML2::IXMLDOMNode* ayadiNode = NULL;
				node->selectSingleNode(_bstr_t("INTEREST_Rate"), &ayadiNode);
				if (ayadiNode != NULL)
				{
					BSTR resultat = NULL;
					ayadiNode->get_text(&resultat);
					_bstr_t ff(resultat,false);
					char * ff1 = (char *)ff;
					oneBoostedRate = (atof(ff1))*100;
					if (resultat) SysFreeString(resultat);
					ayadiNode->Release();
					ayadiNode = NULL;	
				}
				newRefVal = new ARM_ReferenceValue(oneBoostedRate);
				return newRefVal;
			}

			vVal->Elt(i) = atof(buffer);
		
			if (theNode)
				theNode->Release();
		}
		listItem->Release();
		
		//Only one boosted rate in the formula. 
		if ((nbFormulas == 1) || (check == 0.0)) 
		{
			newRefVal = new ARM_ReferenceValue(vVal->Elt(0));
		}
		//Step-Up boosted rates in the fomulasssssssssssss.
		else
		{
			newRefVal = new ARM_ReferenceValue(vDates, vVal, K_YIELD, 1);
			newRefVal->SetCalcMethod(K_STEPUP_LEFT);
		}
	}
	return newRefVal;
}


//The FRC events concerns the Reference Index and the NRC events concerns the Payment Index.
ARM_ReferenceValue* GetCorridorPastFixings(MSXML2::IXMLDOMNode* node, int indexType, const ARM_Date& asOfDate)
{
	MSXML2::IXMLDOMNode* theNode			= NULL;
	MSXML2::IXMLDOMNode* theNode2			= NULL;
	MSXML2::IXMLDOMNode* listItem			= NULL;
	MSXML2::IXMLDOMNodeList* resultList		= NULL;
	BSTR resultat					= NULL;
	long nbEvents					= 0;
	ARM_ReferenceValue* newRefVal	= NULL;
	ARM_Vector* vDates				= NULL;
	ARM_Vector* vVal				= NULL;
	VECTOR<double> vecDates;
	VECTOR<double> vecVal;
	double currFixing;
	ARM_Date currFixingDate;

	node->selectNodes(_bstr_t("Events/EVENT"), &resultList);

	if (resultList != NULL)
	{
		resultList->get_length(&nbEvents);

		for (int i = 0; i < nbEvents; i++)
		{
			HRESULT hr;
			hr = resultList->get_item(i, &listItem);
			listItem->selectSingleNode(_bstr_t("Type"), &theNode);
			
			if (theNode != NULL)
			{
				BSTR resultat = NULL;
				theNode->get_text(&resultat);
				_bstr_t ff(resultat, false);
				char * eventType =(char *)ff;

				//Reference Index: FRC events.
				if ((strcmp(eventType,"FRC") == 0) && (indexType == 0)) 
				{
					currFixing		= 100.0 * GetDoubleFromXMLNode(listItem, "Amount");
					currFixingDate	= GetDateFromXMLNode(listItem, "ADate");
					if (currFixingDate	<= asOfDate)
					{
						vecDates.push_back(currFixingDate.GetJulian());
						vecVal.push_back(currFixing);
					}
				}
				//Payment Index: NRC events.
				else if ((strcmp(eventType,"NRC") == 0) && (indexType == 1))
				{
					currFixing		= 100.0 * GetDoubleFromXMLNode(listItem, "Amount");
					currFixingDate	= GetDateFromXMLNode(listItem, "ADate");
					if (currFixingDate	<= asOfDate)
					{
						vecDates.push_back(currFixingDate.GetJulian());
						vecVal.push_back(currFixing);
					}
				}
			}

			if (theNode)
				theNode->Release();

			if (theNode2)
				theNode2->Release();

			if (listItem)
				listItem->Release();
		}

		if (resultList)
			resultList->Release();
	}

	//No past fixings available.
	if (vecDates.size() == 0)
	{
		return newRefVal;
	}
	//Past fixings available.
	else
	{

		ARM_Vector* VRefDates = CreateARMVectorFromVECTOR(vecDates);
		ARM_Vector* VNotional = CreateARMVectorFromVECTOR(vecVal);
		newRefVal = new ARM_ReferenceValue((ARM_Vector*)VRefDates->Clone(), (ARM_Vector*)VNotional->Clone());
		newRefVal->SetCalcMethod(K_STEPUP_LEFT);
		newRefVal->SetExtrapolMeth(1);
		
		//Free memory.
		delete VRefDates;
		VRefDates = NULL;
		delete VNotional;
		VNotional = NULL;

		return newRefVal;
	}	
}

//The SPR events concern the step-up spread on CorridorLeg or even on SwapLeg. 
ARM_ReferenceValue* GetCorridorStepUpSpread(MSXML2::IXMLDOMNode* node)
{
	MSXML2::IXMLDOMNode* theNode			= NULL;
	MSXML2::IXMLDOMNode* theNode2			= NULL;
	MSXML2::IXMLDOMNode* listItem			= NULL;
	MSXML2::IXMLDOMNodeList* resultList		= NULL;
	BSTR resultat					= NULL;
	long nbEvents					= 0;
	ARM_ReferenceValue* newRefVal	= NULL;
	ARM_Vector* vDates				= NULL;
	ARM_Vector* vVal				= NULL;
	VECTOR<double> vecDates;
	VECTOR<double> vecVal;
	double currSpread;
	ARM_Date currSpreadDate;

	node->selectNodes(_bstr_t("Events/EVENT"), &resultList);

	if (resultList != NULL)
	{
		resultList->get_length(&nbEvents);

		for (int i = 0; i < nbEvents; i++)
		{
			HRESULT hr;
			hr = resultList->get_item(i, &listItem);
			listItem->selectSingleNode(_bstr_t("Type"), &theNode);
			
			if (theNode != NULL)
			{
				BSTR resultat = NULL;
				theNode->get_text(&resultat);
				_bstr_t ff(resultat, false);
				char * eventType =(char *)ff;

				if (strcmp(eventType,"SPR") == 0)  
				{
					currSpread		= 100.0 * GetDoubleFromXMLNode(listItem, "Amount");
					currSpreadDate	= GetDateFromXMLNode(listItem, "ADate") + 1;
					vecDates.push_back(currSpreadDate.GetJulian());
					vecVal.push_back(currSpread);
				}
			}

			if (theNode)
				theNode->Release();

			if (theNode2)
				theNode2->Release();

			if (listItem)
				listItem->Release();
		}

		if (resultList)
			resultList->Release();
	}

	//No past fixings available.
	if (vecDates.size() == 0)
	{
		return newRefVal;
	}
	//Past fixings available.
	else
	{
		ARM_Vector* VRefDates = CreateARMVectorFromVECTOR(vecDates);
		ARM_Vector* VNotional = CreateARMVectorFromVECTOR(vecVal);
		newRefVal = new ARM_ReferenceValue((ARM_Vector*)VRefDates->Clone(), (ARM_Vector*)VNotional->Clone());
		newRefVal->SetCalcMethod(K_STEPUP_LEFT);
		newRefVal->SetExtrapolMeth(1);
		
		//Free memory.
		delete VRefDates;
		VRefDates = NULL;
		delete VNotional;
		VNotional = NULL;

		return newRefVal;
	}	
}


//The SPR events concern the step-up spread on CorridorLeg or even on SwapLeg. 
ARM_ReferenceValue* GetSwaplegStepUpSpread(MSXML2::IXMLDOMNode* node, bool adjustStartDate)
{
	MSXML2::IXMLDOMNode* theNode			= NULL;
	MSXML2::IXMLDOMNode* theNode2			= NULL;
	MSXML2::IXMLDOMNode* listItem			= NULL;
	MSXML2::IXMLDOMNodeList* resultList		= NULL;
	BSTR resultat					= NULL;
	long nbEvents					= 0;
	ARM_ReferenceValue* newRefVal	= NULL;
	ARM_Vector* vDates				= NULL;
	ARM_Vector* vVal				= NULL;
	VECTOR<double> vecDates;
	VECTOR<double> vecVal;
	double currSpread;
	ARM_Date currSpreadDate;

	double firstSpd = GetSpread(node);
	if (firstSpd)
	{
		ARM_Date startDate = GetStartDate(node);
		if ((adjustStartDate) && (GetIntRule(node) == K_ADJUSTED))
		{
			char calendar[4];
			GetResetCalendar(node, calendar);
			startDate.GoodBusinessDay(K_FOLLOWING, calendar); //first start date is always adjusted with FOLLOWING !
		}
		vecDates.push_back(startDate.GetJulian());
		vecVal.push_back(firstSpd * 100.0);
	}
	
	node->selectNodes(_bstr_t("Events/EVENT"), &resultList);

	if (resultList != NULL)
	{
		resultList->get_length(&nbEvents);

		for (int i = 0; i < nbEvents; i++)
		{
			HRESULT hr;
			hr = resultList->get_item(i, &listItem);
			listItem->selectSingleNode(_bstr_t("Type"), &theNode);
			
			if (theNode != NULL)
			{
				BSTR resultat = NULL;
				theNode->get_text(&resultat);
				_bstr_t ff(resultat, false);
				char * eventType =(char *)ff;

				if (strcmp(eventType,"SPR") == 0)  
				{
					currSpread		= 100.0 * GetDoubleFromXMLNode(listItem, "Amount");
					currSpreadDate	= GetDateFromXMLNode(listItem, "Date");
					vecDates.push_back(currSpreadDate.GetJulian());
					vecVal.push_back(currSpread);
				}
			}

			if (theNode)
				theNode->Release();

			if (theNode2)
				theNode2->Release();

			if (listItem)
				listItem->Release();
		}

		if (resultList)
			resultList->Release();
	}

	//No past fixings available.
	if (vecDates.size() == 0)
	{
		return newRefVal;
	}
	//Past fixings available.
	else
	{
		ARM_Vector* VRefDates = CreateARMVectorFromVECTOR(vecDates);
		ARM_Vector* VNotional = CreateARMVectorFromVECTOR(vecVal);
		newRefVal = new ARM_ReferenceValue((ARM_Vector*)VRefDates->Clone(), (ARM_Vector*)VNotional->Clone());
		newRefVal->SetCalcMethod(K_STEPUP_LEFT);
		newRefVal->SetExtrapolMeth(1);
		
		//Free memory.
		delete VRefDates;
		VRefDates = NULL;
		delete VNotional;
		VNotional = NULL;

		return newRefVal;
	}	
}



int GetCorridorIndexType(MSXML2::IXMLDOMNode* node, ARM_Currency* ccy, int refOrPaid)
{
	char summitIndex[15];
	char armIndex[15];
	
	GetNthArgInFormula(node,refOrPaid,summitIndex);

	if (strcmp(summitIndex,"")==0)
		return K_FIXED;
	else
	{
		ConvSummitIndexToARMIndex(summitIndex,ccy,armIndex);

		return ARM_ConvIrType(armIndex);
	}
}

int GetSpreadPayoffLibor(MSXML2::IXMLDOMNode* node, ARM_Currency* ccy)
{
	char summitIndex[20];
	char armIndex[20];
	
	GetNthArgInFormula(node,3,summitIndex);
	if (strcmp(summitIndex,"")==0)
		return K_FIXED;
	else
	{
		char tmp[20];
		ConvSummitIndexToARMIndex(summitIndex,ccy,tmp);
		strcpy(armIndex,tmp);

		return ARM_ConvIrType(armIndex);
	}
}

double GetWeight1(MSXML2::IXMLDOMNode* node)
{
	char arg[15];
	
	GetNthArgInFormula(node, 4, arg);
	return atof((const char*) arg);
	
}

double GetWeight2(MSXML2::IXMLDOMNode* node)
{
	char arg[15];
	
	GetNthArgInFormula(node, 5, arg);
	return atof((const char*) arg);
	
}

ARM_SwapLeg* GetPayLeg(MSXML2::IXMLDOMNode* node, int indexType, ARM_Currency* discountCcy, int basisDaycount, int payFreq, const ARM_Date& endDateNA, ARM_Date asOf)
{
	ARM_SwapLeg* payLeg = NULL;

	ARM_Vector* startDates = NULL;
	ARM_Vector* endDates = NULL;
	ARM_Vector* resetDates = NULL;
	ARM_Vector* paymentDates = NULL;
	ARM_Vector* spread = NULL;
	ARM_Vector* amount = NULL;
	VECTOR<double> rate;
	ARM_Vector* vRate = NULL;
	ARM_Vector* interestDays = NULL;
	ARM_ReferenceValue* notional = NULL;
	ARM_ReferenceValue* variableSpread = NULL;

	IXMLDOMNode* schedLineNode=NULL;
	IXMLDOMNodeList* resultList=NULL;
	HRESULT hr;

	node->selectNodes(_bstr_t("ProdData/cBLOB_I/Schedule1/cBLBLST_I"), &resultList);
	if (resultList!=NULL)
	{
		long nb;
		resultList->get_length(&nb);
		if (nb > 0)
		{
			startDates = new ARM_Vector(nb);
			endDates = new ARM_Vector(nb);
			resetDates = new ARM_Vector(nb);
			paymentDates = new ARM_Vector(nb);
			spread = new ARM_Vector(nb);
			amount = new ARM_Vector(nb);
			//rate = new ARM_Vector(nb);

			for (long i=0 ; i<nb ; i++)
			{	
				hr=resultList->get_item(i, &schedLineNode);
				if (hr==S_OK && schedLineNode!=NULL)
				{
					startDates->Elt(i) = GetDateFromXMLNode(schedLineNode,"StartDate").GetJulian();
					endDates->Elt(i) = GetDateFromXMLNode(schedLineNode,"EndDate").GetJulian();
					resetDates->Elt(i) = GetDateFromXMLNode(schedLineNode,"FixingDate").GetJulian();
					paymentDates->Elt(i) = GetDateFromXMLNode(schedLineNode,"PayDate").GetJulian();

					amount->Elt(i) = GetDoubleFromXMLNode(schedLineNode,"Amount1");
					spread->Elt(i) = GetDoubleFromXMLNode(schedLineNode,"Rate2");
					double dRate;
					
					if (IsFixedIndex((ARM_INDEX_TYPE)indexType))
					{
						dRate = GetDoubleFromXMLNode(schedLineNode,"Rate1");
						rate.push_back(dRate*100.);
					}
					else
					{
						dRate = GetDoubleFromXMLNode(schedLineNode,"Rate1",true);
						if (dRate != ARM_MISSING_VALUE && resetDates->Elt(i) < asOf.GetJulian())
							rate.push_back(dRate*100.);
					}
				}
				schedLineNode->Release();
				schedLineNode = NULL;
			}
		}
		resultList->Release();
		resultList = NULL;
	}

	if (rate.size() > 0)
		vRate = CreateARMVectorFromVECTOR(rate);

	int payBasis = basisDaycount;
	int isBondBasis = GetBLOBNum(node,3);
	if (isBondBasis == 1)
		payBasis = K30_360;

	interestDays = new ARM_Vector(startDates->GetSize());
	for (int i = 0; i < startDates->GetSize(); i++)
		interestDays->Elt(i) = DaysBetweenDates(payBasis, startDates->Elt(i), endDates->Elt(i));

	notional = new ARM_ReferenceValue((ARM_Vector*)paymentDates->Clone(),amount);
	variableSpread = new ARM_ReferenceValue((ARM_Vector*)resetDates->Clone(),spread);

	int resetTiming = GetBLOBTiming(node,2);
	char sFreq[20];
	ltoa(GetBLOBNum(node,2),sFreq,10);
	int resetFreq = FromSummitFreqToARMFreq(sFreq);

	if (IsCMSIndex((ARM_INDEX_TYPE)indexType))
	{
		payLeg = (ARM_SwapLeg *) new ARM_CMSLeg(startDates,
												endDates,
												resetDates,
												paymentDates,
												interestDays,
												(ARM_INDEX_TYPE)indexType,
												notional,
												K_RCV,
												variableSpread,
												K_COMP_PROP,
												payBasis,
												K_ADJUSTED,
												resetFreq,
												payFreq,
												10000,
												discountCcy,
												resetTiming,
												K_SHORTSTART,
												"NULL",
												"NULL");

		payLeg->SetEndDateNA(endDateNA);

		payLeg->SetAmount(notional);

		payLeg->SetIndexType((ARM_INDEX_TYPE)indexType);

		payLeg->SetFixRates(vRate);
	}
	else if (IsLiborIndex((ARM_INDEX_TYPE)indexType))
	{
		ARM_IRIndex newIndex((ARM_INDEX_TYPE)indexType,resetFreq,payFreq,discountCcy);

		payLeg = new ARM_SwapLeg(startDates,
								 endDates,
								 paymentDates,
								 resetDates,
								 interestDays,
								 NULL,
								 notional,
								 &newIndex,
								 K_RCV,
								 0.0,
								 K_FLOAT_RATE,
								 discountCcy,
								 K_NX_NONE,
								 0);

		payLeg->SetEndDateNA(endDateNA);

		payLeg->SetAmount(notional);

		payLeg->SetDayCount(payBasis);

		payLeg->SetVariableSpread(variableSpread);

//		payLeg->CptRealCashFlowDates();

		payLeg->SetFixRates(vRate);
	}
	else if (IsFixedIndex((ARM_INDEX_TYPE)indexType))
	{
		ARM_IRIndex newIndex(discountCcy->GetCcyName(),-1);

		ARM_Vector * fixRates = new ARM_Vector(resetDates->GetSize());
		for (int i = 0; i < fixRates->GetSize(); i++)
			fixRates->Elt(i) = vRate->Elt(i) + spread->Elt(i);

		if (variableSpread)
			delete variableSpread;
		variableSpread = new ARM_ReferenceValue((ARM_Vector*)resetDates->Clone(),fixRates);

		payLeg = new ARM_SwapLeg(startDates,
								 endDates,
								 paymentDates,
								 resetDates,
								 interestDays,
								 NULL,
								 notional,
								 &newIndex,
								 K_RCV,
								 0.0,
								 K_FIXED_RATE,
								 discountCcy);

		payLeg->SetEndDateNA(endDateNA);

		payLeg->SetAmount(notional);

		payLeg->SetDayCount(payBasis);

		payLeg->SetVariableSpread(variableSpread);
	}


	if (variableSpread)
		delete variableSpread;
	variableSpread = NULL;

	if (notional)
		delete notional;
	notional = NULL;

	if (startDates)
		delete startDates;
	startDates = NULL;

	if (endDates)
		delete endDates;
	endDates = NULL;

	if (paymentDates)
		delete paymentDates;
	paymentDates = NULL;

	if (resetDates)
		delete resetDates;
	resetDates = NULL;

	if (interestDays)
		delete interestDays;
	interestDays = NULL;

	if (vRate)
		delete vRate;
	vRate = NULL;

	return payLeg;
}


// Parse the payleg on corridor spread option
void GetPayIdxPastFixings(IXMLDOMNode* node, const ARM_Date& asOf, int indexType, ARM_ReferenceValue* payIdxPastFixings)
{
	ARM_Vector* resetDates = NULL;
	VECTOR<double> rate;
	VECTOR<double> date;
	ARM_Vector* vRate = NULL;
	ARM_Vector* vDate = NULL;

	IXMLDOMNode* schedLineNode=NULL;
	IXMLDOMNodeList* resultList=NULL;
	HRESULT hr;

	node->selectNodes(_bstr_t("ProdData/cBLOB_I/Schedule1/cBLBLST_I"), &resultList);
	if (resultList!=NULL)
	{
		long nb;
		resultList->get_length(&nb);
		if (nb > 0)
		{
			resetDates = new ARM_Vector(nb);

			for (long i=0 ; i<nb ; i++)
			{	
				hr=resultList->get_item(i, &schedLineNode);
				if (hr==S_OK && schedLineNode!=NULL)
				{
					resetDates->Elt(i) = GetDateFromXMLNode(schedLineNode,"FixingDate").GetJulian();

					double dRate;
					
					if (IsFixedIndex((ARM_INDEX_TYPE)indexType))
					{
						dRate = GetDoubleFromXMLNode(schedLineNode,"Rate1");
						{
							rate.push_back(dRate*100.);
							date.push_back(resetDates->Elt(i));
						}
					}
					else
					{
						dRate = GetDoubleFromXMLNode(schedLineNode,"Rate1",true);
						if (dRate != ARM_MISSING_VALUE && resetDates->Elt(i) < asOf.GetJulian())
						{
							rate.push_back(dRate*100.);
							date.push_back(resetDates->Elt(i));
						}
					}
				}
				schedLineNode->Release();
				schedLineNode = NULL;
			}
		}
		resultList->Release();
		resultList = NULL;
	}

	if (rate.size() > 0)
	{
		vRate = CreateARMVectorFromVECTOR(rate);
		vDate = CreateARMVectorFromVECTOR(date);

        *payIdxPastFixings = ARM_ReferenceValue(vDate, vRate); 
	}

	if (resetDates)
	   delete resetDates;
	resetDates = NULL;
}



ARM_IRIndex* GetPayIndex(MSXML2::IXMLDOMNode* node, int indexType, ARM_Currency* indexCcy, int basisDaycount, int payFreq)
{
	ARM_IRIndex* indexP = NULL;

	int payBasis = basisDaycount;
	int isBondBasis = GetBLOBNum(node,3);
	if (isBondBasis == 1)
		payBasis = K30_360;

	int resetTiming = GetBLOBTiming(node,2);

	char sFreq[20];
	ltoa(GetBLOBNum(node,2),sFreq,10);
	int resetFreq = FromSummitFreqToARMFreq(sFreq);

	if (IsCMSIndex((ARM_INDEX_TYPE)indexType))
	{
		indexP = new ARM_IRIndex((ARM_INDEX_TYPE)indexType, indexCcy->GetVanillaIndexType(), resetFreq, K_ANNUAL, indexCcy);
	}
	else if (IsLiborIndex((ARM_INDEX_TYPE)indexType))
	{
		indexP = new ARM_IRIndex((ARM_INDEX_TYPE)indexType, resetFreq, payFreq, indexCcy, basisDaycount);
	}
	else if (IsFixedIndex((ARM_INDEX_TYPE)indexType))
	{
		indexP = new ARM_IRIndex(indexCcy->GetCcyName(),payBasis);

		indexP->SetResetFrequency(resetFreq);
		indexP->SetPayFrequency(payFreq);
		indexP->SetResetTiming(resetTiming);
	}

	return indexP;

}

ARM_Vector* GetSpreadFirstFixing(MSXML2::IXMLDOMNode* node, const ARM_Date& date)
{
	MSXML2::IXMLDOMNode* theNode=NULL;
	MSXML2::IXMLDOMNode* theNode2=NULL;
	MSXML2::IXMLDOMNodeList* resultList = NULL;
	BSTR resultat = NULL;
	HRESULT hr;
	ARM_Date FixingDate;

	VECTOR<double> VVal;
	
	// Recuperation des evenements							{
	if( node->selectNodes(_bstr_t("Events/EVENT"), &resultList)== S_OK)
	{
		if (resultList!=NULL)
		{
			long nb;
			resultList->get_length(&nb);
			for (long i=0 ; i<nb ; i++)
			{	
				hr=resultList->get_item(i, &theNode);
				if (hr==S_OK && theNode!=NULL)
				{
					theNode->selectSingleNode(_bstr_t("Type"), &theNode2);
					if (theNode2!=NULL)
					{
						theNode2->get_text(&resultat);
						_bstr_t ff(resultat,false);
						char * ff1=(char *)ff;
						if (strcmp(ff1,"FRC")==0)
						{
							double fix= 100.0*GetDoubleFromXMLNode(theNode,"Amount");
							FixingDate = GetDateFromXMLNode(theNode, "ADate");
							if (FixingDate <= date)
									VVal.push_back(fix);
						}
						if (resultat) SysFreeString(resultat);
						theNode2->Release();
						theNode2=NULL;	
					}
					theNode->Release();
					theNode=NULL;
				}
			}
			if (resultList) resultList->Release();
			resultList=NULL;
		}
	}
	ARM_Vector* vectorVal= CreateARMVectorFromVECTOR(VVal);
	
	return vectorVal;
}


ARM_ReferenceValue* GetRefValSpreadFirstFixing(MSXML2::IXMLDOMNode* node, const ARM_Date& date)
{
	MSXML2::IXMLDOMNode* theNode=NULL;
	MSXML2::IXMLDOMNode* theNode2=NULL;
	MSXML2::IXMLDOMNodeList* resultList = NULL;
	BSTR resultat = NULL;
	HRESULT hr;
	ARM_Date FixingDate;

	VECTOR<double> VVal;
	VECTOR<double> VDate;
	
	// Recuperation des evenements							{
	if( node->selectNodes(_bstr_t("Events/EVENT"), &resultList)== S_OK)
	{
		if (resultList!=NULL)
		{
			long nb;
			resultList->get_length(&nb);
			for (long i=0 ; i<nb ; i++)
			{	
				hr=resultList->get_item(i, &theNode);
				if (hr==S_OK && theNode!=NULL)
				{
					theNode->selectSingleNode(_bstr_t("Type"), &theNode2);
					if (theNode2!=NULL)
					{
						theNode2->get_text(&resultat);
						_bstr_t ff(resultat,false);
						char * ff1=(char *)ff;
						if (strcmp(ff1,"FRC")==0)
						{
							double fix= 100.0*GetDoubleFromXMLNode(theNode,"Amount");
							FixingDate = GetDateFromXMLNode(theNode, "ADate");
							if (FixingDate <= date)
							{
								VVal.push_back(fix);
								VDate.push_back(FixingDate.GetJulian());
							}
						}
						if (resultat) SysFreeString(resultat);
						theNode2->Release();
						theNode2=NULL;	
					}
					theNode->Release();
					theNode=NULL;
				}
			}
			if (resultList) resultList->Release();
			resultList=NULL;
		}
	}

	ARM_Vector* vectorVal= CreateARMVectorFromVECTOR(VVal);
	ARM_Vector* vectorDate= CreateARMVectorFromVECTOR(VDate);

	ARM_ReferenceValue* newRefVal = NULL;

	if (vectorVal)
		newRefVal = new ARM_ReferenceValue(vectorDate,vectorVal);
	
	return newRefVal;
}



ARM_Vector* GetSpreadSecondFixing(MSXML2::IXMLDOMNode* node, const ARM_Date& date)
{
	MSXML2::IXMLDOMNode* theNode=NULL;
	MSXML2::IXMLDOMNode* theNode2=NULL;
	MSXML2::IXMLDOMNodeList* resultList = NULL;
	BSTR resultat = NULL;
	HRESULT hr;
	ARM_Date FixingDate;

	VECTOR<double> VVal;
	
	// Recuperation des evenements							{
	if( node->selectNodes(_bstr_t("Events/EVENT"), &resultList)== S_OK)
	{
		if (resultList!=NULL)
		{
			long nb;
			resultList->get_length(&nb);
			for (long i=0 ; i<nb ; i++)
			{	
				hr=resultList->get_item(i, &theNode);
				if (hr==S_OK && theNode!=NULL)
				{
					
					theNode->selectSingleNode(_bstr_t("Type"), &theNode2);
					
					if (theNode2!=NULL)
					{
						theNode2->get_text(&resultat);
						_bstr_t ff(resultat,false);
						char * ff1=(char *)ff;
						if (strcmp(ff1,"NRC")==0)
						{
							double fix= 100.0*GetDoubleFromXMLNode(theNode,"Amount");
							FixingDate = GetDateFromXMLNode(theNode, "ADate");
							if (FixingDate <= date)
								VVal.push_back(fix);
						}

						if (resultat) SysFreeString(resultat);
						theNode2->Release();
						theNode2=NULL;						
					}
					theNode->Release();
					theNode=NULL;
				}
			}
			if (resultList) resultList->Release();
			resultList=NULL;
		}
	}
	ARM_Vector* vectorVal= CreateARMVectorFromVECTOR(VVal);
	
	return vectorVal;
}


ARM_ReferenceValue* GetRefValSpreadSecondFixing(MSXML2::IXMLDOMNode* node, const ARM_Date& date)
{
	MSXML2::IXMLDOMNode* theNode=NULL;
	MSXML2::IXMLDOMNode* theNode2=NULL;
	MSXML2::IXMLDOMNodeList* resultList = NULL;
	BSTR resultat = NULL;
	HRESULT hr;
	ARM_Date FixingDate;

	VECTOR<double> VVal;
	VECTOR<double> VDate;
	
	// Recuperation des evenements							{
	if( node->selectNodes(_bstr_t("Events/EVENT"), &resultList)== S_OK)
	{
		if (resultList!=NULL)
		{
			long nb;
			resultList->get_length(&nb);
			for (long i=0 ; i<nb ; i++)
			{	
				hr=resultList->get_item(i, &theNode);
				if (hr==S_OK && theNode!=NULL)
				{
					theNode->selectSingleNode(_bstr_t("Type"), &theNode2);
					if (theNode2!=NULL)
					{
						theNode2->get_text(&resultat);
						_bstr_t ff(resultat,false);
						char * ff1=(char *)ff;
						if (strcmp(ff1,"NRC")==0)
						{
							double fix= 100.0*GetDoubleFromXMLNode(theNode,"Amount");
							FixingDate = GetDateFromXMLNode(theNode, "ADate");
							if (FixingDate <= date)
							{
								VVal.push_back(fix);
								VDate.push_back(FixingDate.GetJulian());
							}
						}
						if (resultat) SysFreeString(resultat);
						theNode2->Release();
						theNode2=NULL;	
					}
					theNode->Release();
					theNode=NULL;
				}
			}
			if (resultList) resultList->Release();
			resultList=NULL;
		}
	}

	ARM_Vector* vectorVal= CreateARMVectorFromVECTOR(VVal);
	ARM_Vector* vectorDate= CreateARMVectorFromVECTOR(VDate);

	ARM_ReferenceValue* newRefVal = NULL;

	if (vectorVal)
		newRefVal = new ARM_ReferenceValue(vectorDate,vectorVal);
	
	return newRefVal;
}



void GetEcheancier(MSXML2::IXMLDOMNode* node, int resetFreq, int payFreq, ARM_Vector** startDates, ARM_Vector** endDates, ARM_Vector** resetDates, ARM_Vector** paymentDates, ARM_Date asOf)
{
	MSXML2::IXMLDOMNode* theNode=NULL;
	MSXML2::IXMLDOMNode* theNode2=NULL;
	MSXML2::IXMLDOMNodeList* resultList = NULL;
	BSTR resultat = NULL;
	HRESULT hr;

	VECTOR<double> VStartDates;
	VECTOR<double> VEndDates;
	VECTOR<double> VResetDates;

	VECTOR<double> VEndDatesForPaymentDates;
	VECTOR<double> VPaymentDates;

	char ProductName[30];
	GetFormula(node,ProductName);

	// Recuperation des evenements							{
	if( node->selectNodes(_bstr_t("Events/EVENT"), &resultList)== S_OK)
	{
		if (resultList!=NULL)
		{
			long nb;
			resultList->get_length(&nb);

			int firstROF = 1;
			int previousFRC = 0;

			for (long i=0 ; i<nb ; i++)
			{	
				hr=resultList->get_item(i, &theNode);
				if (hr==S_OK && theNode!=NULL)
				{
					theNode->selectSingleNode(_bstr_t("Type"), &theNode2);
					
					if (theNode2!=NULL)
					{
						theNode2->get_text(&resultat);
						_bstr_t ff(resultat,false);
						char * ff1=(char *)ff;

						if (strcmp(ff1,"ROF")==0)
						{
							ARM_Date tmpDate = GetDateFromXMLNode(theNode, "Date");
							ARM_Date tmpDate1 = GetDateFromXMLNode(theNode, "ADate");
							ARM_Date tmpDate2 = GetDateFromXMLNode(theNode, "BDate");

							if (tmpDate1 != ARM_DEFAULT_DATE) 
							{
								if (tmpDate2 != ARM_DEFAULT_DATE)
								{
									if (strcmp(ProductName, "GLOBALCAP") || (tmpDate2.GetJulian() >= asOf.GetJulian()))
									{
										VResetDates.push_back(tmpDate2.GetJulian());
										VStartDates.push_back(tmpDate.GetJulian());
									}
								}
								else
								{
									if (strcmp(ProductName, "GLOBALCAP") || (tmpDate1.GetJulian() >= asOf.GetJulian()))
									{
										VResetDates.push_back(tmpDate1.GetJulian());
										VStartDates.push_back(tmpDate.GetJulian());
									}
								}
							}

							if (firstROF == 0)
							{
								if (strcmp(ProductName, "GLOBALCAP") || (tmpDate.GetJulian() >= asOf.GetJulian()))
								{
									VEndDatesForPaymentDates.push_back(tmpDate.GetJulian());
									VPaymentDates.push_back(tmpDate.GetJulian());
								}
							}

							previousFRC = 0;
							firstROF = 0;
						}
						else if ( (strcmp(ff1,"RES")==0) || (strcmp(ff1,"FRC")==0) )
						{
							ARM_Date tmpDate = GetDateFromXMLNode(theNode, "ADate");
							if (tmpDate != ARM_DEFAULT_DATE)
							{
								if (strcmp(ProductName, "GLOBALCAP") || (tmpDate.GetJulian() >= asOf.GetJulian()))
								{
									VResetDates.push_back(tmpDate.GetJulian());
									VStartDates.push_back(GetDateFromXMLNode(theNode, "Date").GetJulian());
								}
							}
							previousFRC = 1;
						}
						else if (strcmp(ff1,"OFF")==0)
						{
							ARM_Date tmpDate = GetDateFromXMLNode(theNode, "ADate");
							if (tmpDate != ARM_DEFAULT_DATE)
							{
								if (strcmp(ProductName, "GLOBALCAP") || (tmpDate.GetJulian() >= asOf.GetJulian()))
								{
									VEndDatesForPaymentDates.push_back(GetDateFromXMLNode(theNode, "Date").GetJulian());
									VPaymentDates.push_back(GetDateFromXMLNode(theNode, "ADate").GetJulian());
								}
							}
							else
							{
								// WARNING : ADate peut ne pas etre définie
								// dans ce cas, il faut récupérer le champ Date
								double tmpPayDate = GetDateFromXMLNode(theNode, "Date").GetJulian();
								if (tmpPayDate != GetStartDate(node).GetJulian())
								{
									if (strcmp(ProductName, "GLOBALCAP") || (tmpPayDate >= asOf.GetJulian()))
									{
										VPaymentDates.push_back(tmpPayDate);
										VEndDatesForPaymentDates.push_back(tmpPayDate);
									}
									
									// pour les cas reset moins fréquents que pay
									if ( (resetFreq < payFreq) && (payFreq != K_ZEROCOUPON) && (previousFRC == 0) )
									{
										if (strcmp(ProductName, "GLOBALCAP") || (VResetDates[VResetDates.size()-1] >= asOf.GetJulian()))
										{
											VResetDates.push_back(VResetDates[VResetDates.size()-1]);
											VStartDates.push_back(tmpPayDate);
										}
									}
								}
							}
							previousFRC = 0;
							firstROF = 0;
						}

						if (resultat) SysFreeString(resultat);
						theNode2->Release();
						theNode2=NULL;
					}
					theNode->Release();
					theNode=NULL;
				}
			}
			if (resultList) resultList->Release();
			resultList=NULL;
		}
	}

	if ((resetFreq < payFreq) && (payFreq != K_ZEROCOUPON) )
	{
		*startDates = CreateARMVectorFromVECTOR(VStartDates,VStartDates.size()-1);
		*resetDates = CreateARMVectorFromVECTOR(VResetDates,VResetDates.size()-1);
	}
	else
	{
		*startDates = CreateARMVectorFromVECTOR(VStartDates);
		*resetDates = CreateARMVectorFromVECTOR(VResetDates);
	}

	*endDates = new ARM_Vector((*startDates)->GetSize());
	for (int i = 0; i < (*endDates)->GetSize()-1; i++)
		(*endDates)->Elt(i) = VStartDates[i+1];

	(*endDates)->Elt((*endDates)->GetSize()-1) = VEndDatesForPaymentDates[VEndDatesForPaymentDates.size()-1];

	*paymentDates = new ARM_Vector((*startDates)->GetSize(), 0.0);

	if (payFreq == K_ZEROCOUPON)
	{
		// Il ne suffit pas de stocker la seule date de paiement !
		// Il faut egalement récupérer les autres dates
		// qui permettront de gérer le cas d'un éventuel amortissement de notionnel
		for (i = 0; i < (*paymentDates)->GetSize()-1; i++)
			(*paymentDates)->Elt(i) = (*endDates)->Elt(i);

		(*paymentDates)->Elt((*paymentDates)->GetSize()-1) = VPaymentDates[VPaymentDates.size()-1];
	}
	else
	{
		// WARNING : case payFreq < resetFreq
		// VPaymentDates size may be smaller than the other vectors size
		// => complete paymentDates !
		if (payFreq < resetFreq)
		{
			// suppose that payment gap = 0
			int paySize = VPaymentDates.size();
			int payTiming = GetPayTiming(node);
			int j = 0;
			if (payTiming == K_ADVANCE)
			{
				for (int i = 0; (i < VStartDates.size()) && (j < paySize); i++)
				{
					if (VStartDates[i] <= VPaymentDates[j])
						(*paymentDates)->Elt(i) = VPaymentDates[j];
					else
						j++;
				}
			}
			else
			{
				for (int i = 0; (i < VEndDates.size()) && (j < paySize); i++)
				{
					if (VEndDates[i] <= VPaymentDates[j])
						(*paymentDates)->Elt(i) = VPaymentDates[j];
					else
						j++;
				}
			}
		}
		else
		{
			for (i = 0; i < (*paymentDates)->GetSize(); i++)
			{
				(*paymentDates)->Elt(i) = VPaymentDates[i];
			}
		}
	}
}


ARM_Date GetBLOBDate(MSXML2::IXMLDOMNode* node, int n)
{	
	char nodeName[30];
	sprintf(nodeName, "ProdData/cBLOB_I/Date%i", n);
	
	return GetDateFromXMLNode(node, nodeName).GetJulian();
}


double GetBLOBAmount(MSXML2::IXMLDOMNode* node, int n)
{
	char nodeName[30];
	sprintf(nodeName, "ProdData/cBLOB_I/Amount%i", n);
	
	return GetDoubleFromXMLNode(node,nodeName);
}


int GetBLOBNum(MSXML2::IXMLDOMNode* node, int n)
{
	char nodeName[30];
	sprintf(nodeName, "ProdData/cBLOB_I/Num%i", n);
	
	return GetIntFromXMLNode(node,nodeName);
}


double GetBLOBRate(MSXML2::IXMLDOMNode* node, int n)
{
	char nodeName[30];
	sprintf(nodeName, "ProdData/cBLOB_I/Rate%i", n);
	
	return GetDoubleFromXMLNode(node,nodeName);	
}



char* GetBLOBIdxTerm(MSXML2::IXMLDOMNode* node, ARM_Currency* ccy, int n)
{
	MSXML2::IXMLDOMNode* theNode=NULL;
	BSTR resultat = NULL;
	char summitIndex[15];

	char nodeName[30];
	sprintf(nodeName, "ProdData/cBLOB_I/Index%i", n);
	
	node->selectSingleNode(_bstr_t(nodeName), &theNode);
	if (theNode!=NULL)
	{
		BSTR resultat = NULL;
		theNode->get_text(&resultat);

		_bstr_t ff(resultat,false);
		char * ff1=(char *)ff;
		strcpy(summitIndex,ff1);
		
		theNode->Release();
		theNode=NULL;
	}
	if (resultat) SysFreeString(resultat);
		
	char* tmp = new char[20];
	ConvSummitIndexToARMIndexTerm(summitIndex,ccy,tmp);
	return tmp;

}

int GetBLOBBasis(MSXML2::IXMLDOMNode* node, int n)
{
	MSXML2::IXMLDOMNode* theNode = NULL;
	BSTR resultat        = NULL;
	int value            = KNOBASE;
	
	char nodeName[30];
	sprintf(nodeName, "ProdData/cBLOB_I/Basis%i", n);
	
	node->selectSingleNode(_bstr_t(nodeName), &theNode);
	if (theNode!=NULL)
	{
		theNode->get_text(&resultat);

		_bstr_t ff(resultat,false);
		char * ff1=(char *)ff;
		
		if (strcmp((const char*) ff1,"") != 0)
			value = FromSummitDaycountToARMDaycount((const char*) ff1);
		
		theNode->Release();
		theNode=NULL;
		if (resultat) SysFreeString(resultat);
	}
	else
	{
		CCString msg((CCString)"Pb in getting Date from node " + nodeName);

		throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT, (char*) msg);
	}

	return value;
}

char* GetBLOBResetCal(MSXML2::IXMLDOMNode* node, int n)
{
	MSXML2::IXMLDOMNode* theNode = NULL;
	BSTR resultat        = NULL;
	char* resetCal       = NULL;
	
	char nodeName[30];
	sprintf(nodeName, "ProdData/cBLOB_I/ResetCal%i", n);
	
	node->selectSingleNode(_bstr_t(nodeName), &theNode);
	if (theNode!=NULL)
	{
		BSTR resultat = NULL;
		theNode->get_text(&resultat);

		_bstr_t ff(resultat,false);
		char * ff1=(char *)ff;

		resetCal = new char[4];
		strcpy(resetCal,ff1);

		theNode->Release();
		theNode=NULL;
		if (resultat) SysFreeString(resultat);
	}
	else
	{
		CCString msg((CCString)"Pb in getting Date from node " + nodeName);

		throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT, (char*) msg);
	}

	return resetCal;
}
					  
char* GetBLOBString(MSXML2::IXMLDOMNode* node, int n)
{
	MSXML2::IXMLDOMNode* theNode = NULL;
	BSTR resultat        = NULL;
	char* resetCal       = NULL;
	
	char nodeName[30];
	sprintf(nodeName, "ProdData/cBLOB_I/String%i", n);
	
	node->selectSingleNode(_bstr_t(nodeName), &theNode);
	if (theNode!=NULL)
	{
		BSTR resultat = NULL;
		theNode->get_text(&resultat);

		_bstr_t ff(resultat,false);
		char * ff1=(char *)ff;

		resetCal = new char[4];
		strcpy(resetCal,ff1);

		theNode->Release();
		theNode=NULL;
		if (resultat) SysFreeString(resultat);
	}
	else
	{
		CCString msg((CCString)"Pb in getting Date from node " + nodeName);

		throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT, (char*) msg);
	}

	return resetCal;
}
					  
int GetBLOBResetGap(MSXML2::IXMLDOMNode* node, int n)
{
	MSXML2::IXMLDOMNode* theNode = NULL;
	BSTR resultat        = NULL;
	int value            = 0;
	
	char nodeName[30];
	sprintf(nodeName, "ProdData/cBLOB_I/ResetGap%i", n);
	
	node->selectSingleNode(_bstr_t(nodeName), &theNode);
	if (theNode!=NULL)
	{
		BSTR resultat = NULL;
		theNode->get_text(&resultat);

		_bstr_t ff(resultat,false);
		char * ff1=(char *)ff;

		value = FromSummitGapToARMGap((const char*) ff1);

		theNode->Release();
		theNode=NULL;
		if (resultat) SysFreeString(resultat);
	}
	else
	{
		CCString msg((CCString)"Pb in getting Date from node " + nodeName);

		throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT, (char*) msg);
	}

	return value;
}

// ARM_ReferenceValue is given with startdate
ARM_ReferenceValue* GetBLOBSchedAmount(MSXML2::IXMLDOMNode* node, int schedIdx, int amountIdx)
{
	MSXML2::IXMLDOMNode*     schedLineNode = NULL;
	MSXML2::IXMLDOMNode*     amountNode    = NULL;
	MSXML2::IXMLDOMNodeList* resultList    = NULL;

	ARM_ReferenceValue* refVal = NULL;

	BSTR resultat = NULL;
	HRESULT hr;

	long nb=0;
	
	VECTOR<double> dates;
	VECTOR<double> amounts;

	char nodeName[40];
	sprintf(nodeName, "ProdData/cBLOB_I/Schedule%i/cBLBLST_I", schedIdx);
	
	node->selectNodes(_bstr_t(nodeName), &resultList);
	if (resultList!=NULL)
	{
		long nb;
		resultList->get_length(&nb);
		if (0 < nb)
		{
			sprintf(nodeName, "Amount%i", amountIdx);

			for (long i=0 ; i<nb ; i++)
			{	
				hr=resultList->get_item(i, &schedLineNode);
				if (hr==S_OK && schedLineNode!=NULL)
				{
					dates.push_back(GetDateFromXMLNode(schedLineNode,"StartDate").GetJulian());

					schedLineNode->selectSingleNode(_bstr_t(nodeName), &amountNode);
					if (amountNode!=NULL)
					{
						amountNode->get_text(&resultat);
						_bstr_t ff(resultat,false);
						char * ff1=(char *)ff;

						amounts.push_back(atof(ff1));

						if (resultat) SysFreeString(resultat);
						amountNode->Release();
						amountNode = NULL;
					}

					schedLineNode->Release();
					schedLineNode = NULL;
				}
			}

			if (nb == 1)
				refVal = new ARM_ReferenceValue(amounts[0]);
			else
			{
				ARM_Vector* VDates   = CreateARMVectorFromVECTOR(dates);
				ARM_Vector* VAmounts = CreateARMVectorFromVECTOR(amounts);

				refVal = new ARM_ReferenceValue(VDates, VAmounts);
			}
		}
		resultList->Release();
		resultList = NULL;
	}

	return refVal;
}

// ARM_ReferenceValue
ARM_ReferenceValue* GetBLOBSchedRate(MSXML2::IXMLDOMNode* node, const char* dateName,int schedIdx, int rateIdx)
{
	MSXML2::IXMLDOMNode*     schedLineNode = NULL;
	MSXML2::IXMLDOMNode*     amountNode    = NULL;
	MSXML2::IXMLDOMNodeList* resultList    = NULL;

	ARM_ReferenceValue* refVal = NULL;

	BSTR resultat = NULL;
	HRESULT hr;

	long nb=0;
	
	VECTOR<double> dates;
	VECTOR<double> amounts;

	char nodeName[40];
	sprintf(nodeName, "ProdData/cBLOB_I/Schedule%i/cBLBLST_I", schedIdx);
	
	node->selectNodes(_bstr_t(nodeName), &resultList);
	if (resultList!=NULL)
	{
		long nb;
		resultList->get_length(&nb);
		if (0 < nb)
		{
			sprintf(nodeName, "Rate%i", rateIdx);

			for (long i=0 ; i<nb ; i++)
			{	
				hr=resultList->get_item(i, &schedLineNode);
				if (hr==S_OK && schedLineNode!=NULL)
				{
					dates.push_back(GetDateFromXMLNode(schedLineNode,dateName).GetJulian());

					schedLineNode->selectSingleNode(_bstr_t(nodeName), &amountNode);
					if (amountNode!=NULL)
					{
						amountNode->get_text(&resultat);
						_bstr_t ff(resultat,false);
						char * ff1=(char *)ff;

						amounts.push_back(atof(ff1));

						if (resultat) SysFreeString(resultat);
						amountNode->Release();
						amountNode = NULL;
					}

					schedLineNode->Release();
					schedLineNode = NULL;
				}
			}

			if (nb == 1)
				refVal = new ARM_ReferenceValue(amounts[0]);
			else
			{
				ARM_Vector* VDates   = CreateARMVectorFromVECTOR(dates);
				ARM_Vector* VAmounts = CreateARMVectorFromVECTOR(amounts);

				refVal = new ARM_ReferenceValue(VDates, VAmounts);
			}
		}
		resultList->Release();
		resultList = NULL;
	}

	return refVal;
}

// ARM_ReferenceValue
ARM_ReferenceValue* GetBLOBSchedRateMulti(MSXML2::IXMLDOMNode* node, const char* dateName,int schedIdx, int rateIdx, double factorMulti)
{
	MSXML2::IXMLDOMNode*     schedLineNode = NULL;
	MSXML2::IXMLDOMNode*     amountNode    = NULL;
	MSXML2::IXMLDOMNodeList* resultList    = NULL;

	ARM_ReferenceValue* refVal = NULL;

	BSTR resultat = NULL;
	HRESULT hr;

	long nb=0;
	
	VECTOR<double> dates;
	VECTOR<double> amounts;

	char nodeName[40];
	sprintf(nodeName, "ProdData/cBLOB_I/Schedule%i/cBLBLST_I", schedIdx);
	
	node->selectNodes(_bstr_t(nodeName), &resultList);
	if (resultList!=NULL)
	{
		long nb;
		resultList->get_length(&nb);
		if (0 < nb)
		{
			sprintf(nodeName, "Rate%i", rateIdx);

			for (long i=0 ; i<nb ; i++)
			{	
				hr=resultList->get_item(i, &schedLineNode);
				if (hr==S_OK && schedLineNode!=NULL)
				{
					dates.push_back(GetDateFromXMLNode(schedLineNode,dateName).GetJulian());

					schedLineNode->selectSingleNode(_bstr_t(nodeName), &amountNode);
					if (amountNode!=NULL)
					{
						amountNode->get_text(&resultat);
						_bstr_t ff(resultat,false);
						char * ff1=(char *)ff;

						amounts.push_back(atof(ff1)*factorMulti);

						if (resultat) SysFreeString(resultat);
						amountNode->Release();
						amountNode = NULL;
					}

					schedLineNode->Release();
					schedLineNode = NULL;
				}
			}

			if (nb == 1)
				refVal = new ARM_ReferenceValue(amounts[0]);
			else
			{
				ARM_Vector* VDates   = CreateARMVectorFromVECTOR(dates);
				ARM_Vector* VAmounts = CreateARMVectorFromVECTOR(amounts);

				refVal = new ARM_ReferenceValue(VDates, VAmounts);
			}
		}
		resultList->Release();
		resultList = NULL;
	}

	return refVal;
}

// ARM_ReferenceValue
ARM_ReferenceValue* GetBLOBSchedRateMultiWithLastDate(MSXML2::IXMLDOMNode* node, const char* dateName,const ARM_Date& date,int schedIdx, int rateIdx, double factorMulti)
{
	MSXML2::IXMLDOMNode*     schedLineNode = NULL;
	MSXML2::IXMLDOMNode*     amountNode    = NULL;
	MSXML2::IXMLDOMNodeList* resultList    = NULL;

	ARM_ReferenceValue* refVal = NULL;

	BSTR resultat = NULL;
	HRESULT hr;

	ARM_Date SelectedDate;

	long nb=0;
	
	VECTOR<double> dates;
	VECTOR<double> amounts;
	long nbElement=0;

	char nodeName[40];
	sprintf(nodeName, "ProdData/cBLOB_I/Schedule%i/cBLBLST_I", schedIdx);
	
	node->selectNodes(_bstr_t(nodeName), &resultList);
	if (resultList!=NULL)
	{
		long nb;
		resultList->get_length(&nb);
		if (0 < nb)
		{
			sprintf(nodeName, "Rate%i", rateIdx);

			for (long i=0 ; i<nb ; i++)
			{	
				hr=resultList->get_item(i, &schedLineNode);
				if (hr==S_OK && schedLineNode!=NULL)
				{
					SelectedDate = GetDateFromXMLNode(schedLineNode,dateName);
					if (SelectedDate<=date)
					{
						schedLineNode->selectSingleNode(_bstr_t(nodeName), &amountNode);
						if (amountNode!=NULL)
						{
							nbElement=nbElement+1;
							amountNode->get_text(&resultat);
							_bstr_t ff(resultat,false);
							char * ff1=(char *)ff;

							amounts.push_back(atof(ff1)*factorMulti);
							dates.push_back(SelectedDate.GetJulian());

							if (resultat) SysFreeString(resultat);
							amountNode->Release();
							amountNode = NULL;
						}	
					}
					schedLineNode->Release();
					schedLineNode = NULL;
				}
			}

			if (nbElement >0)
			{
				if (nbElement == 1)
					refVal = new ARM_ReferenceValue(amounts[0]);
				else
				{
					ARM_Vector* VDates   = CreateARMVectorFromVECTOR(dates);
					ARM_Vector* VAmounts = CreateARMVectorFromVECTOR(amounts);

					refVal = new ARM_ReferenceValue(VDates, VAmounts);
				}
			}
		}
		resultList->Release();
		resultList = NULL;
	}

	return refVal;
}

// Specific Data 
ARM_Date GetStartDateSpecific(MSXML2::IXMLDOMNode* node)
{	
	return GetDateFromXMLNode(node,"SPECIFIC_DATA/STARTDATE");
}

double GetStrikeSpecific(MSXML2::IXMLDOMNode* node)
{
	return GetDoubleFromXMLNode(node, "SPECIFIC_DATA/STRIKE");
}

double GetNotionalSpecific(MSXML2::IXMLDOMNode* node)
{
	return GetDoubleFromXMLNode(node, "SPECIFIC_DATA/NOTIONAL");
}

void GetSwapCustomSchedule(MSXML2::IXMLDOMNode* node, 
						   ARM_Vector** startDates, 
						   ARM_Vector** endDates, 
						   ARM_Vector** resetDates, 
						   ARM_Vector** paymentDates)
{
	MSXML2::IXMLDOMNode* theNode=NULL;
	MSXML2::IXMLDOMNode* theNode2=NULL;
	MSXML2::IXMLDOMNodeList* resultList = NULL;
	BSTR resultat = NULL;
	HRESULT hr;

	VECTOR<double> VStartDates;
	VECTOR<double> VEndDates;
	VECTOR<double> VResetDates;
	VECTOR<double> VPaymentDates;
	VECTOR<double> VFixings1;
	VECTOR<double> VFixings2;

	// Recuperation des evenements							{
	if( node->selectNodes(_bstr_t("Events/EVENT"), &resultList)== S_OK)
	{
		if (resultList!=NULL)
		{
			long nb;
			resultList->get_length(&nb);

			int firstROF = 1;
			for (long i=0 ; i<nb ; i++)
			{	
				hr=resultList->get_item(i, &theNode);
				if (hr==S_OK && theNode!=NULL)
				{
					theNode->selectSingleNode(_bstr_t("Type"), &theNode2);
					
					if (theNode2!=NULL)
					{
						theNode2->get_text(&resultat);
						_bstr_t ff(resultat,false);
						char * ff1=(char *)ff;

						if (strcmp(ff1,"ROF")==0)
						{
							ARM_Date tmpDate = GetDateFromXMLNode(theNode, "ADate");

							if (tmpDate != ARM_DEFAULT_DATE)
							{
								VStartDates.push_back(GetDateFromXMLNode(theNode, "Date").GetJulian());
								VResetDates.push_back(tmpDate.GetJulian());
							}
							if (firstROF == 0)
							{
								//EndDate[i]=StartDate[i+1]
								VEndDates.push_back(GetDateFromXMLNode(theNode, "Date").GetJulian());
								//Payment=End
								VPaymentDates.push_back(GetDateFromXMLNode(theNode, "Date").GetJulian());
							}
							firstROF = 0;
						}
						else if ( (strcmp(ff1,"RES")==0) || (strcmp(ff1,"FRC")==0) )
						{
							ARM_Date tmpDate = GetDateFromXMLNode(theNode, "ADate");
							if (tmpDate != ARM_DEFAULT_DATE)
							{
								VResetDates.push_back(tmpDate.GetJulian());
								VStartDates.push_back(GetDateFromXMLNode(theNode, "Date").GetJulian());
							}
						}
						else if (strcmp(ff1,"OFF")==0)
						{
							//Last end date
							VEndDates.push_back(GetDateFromXMLNode(theNode, "Date").GetJulian());
							//Payment=End
							VPaymentDates.push_back(GetDateFromXMLNode(theNode, "Date").GetJulian());
						}

						if (resultat) SysFreeString(resultat);
						theNode2->Release();
						theNode2=NULL;
					}
					theNode->Release();
					theNode=NULL;
				}
			}
			if (resultList) resultList->Release();
			resultList=NULL;
		}
	}

	// WARNING : case payFreq < resetFreq
	// VPaymentDates size may be smaller than the other vectors size
	// => complete paymentDates !
	if (VPaymentDates < VResetDates)
	{
		// suppose that payment gap = 0
		VECTOR<double> VPaymentDatesNew(VResetDates.size());
		int paySize = VPaymentDates.size();
		int payTiming = GetPayTiming(node);
		int j = 0;
		if (payTiming == K_ADVANCE)
		{
			for (int i = 0; (i < VStartDates.size()) && (j < paySize); i++)
			{
				if (VStartDates[i] <= VPaymentDates[j])
					VPaymentDatesNew[i] = VPaymentDates[j];
				else
					j++;
			}
		}
		else
		{
			for (int i = 0; (i < VEndDates.size()) && (j < paySize); i++)
			{
				if (VEndDates[i] <= VPaymentDates[j])
					VPaymentDatesNew[i] = VPaymentDates[j];
				else
					j++;
			}
		}
		*paymentDates = CreateARMVectorFromVECTOR(VPaymentDatesNew);
	}
	else
		*paymentDates = CreateARMVectorFromVECTOR(VPaymentDates);

	*startDates = CreateARMVectorFromVECTOR(VStartDates);
	*resetDates = CreateARMVectorFromVECTOR(VResetDates);
	*endDates = CreateARMVectorFromVECTOR(VEndDates);
}


/// Cette fonction permet dans le cas d'un schedule STND d'une spread option
/// de prendre les eventuels fixings qui ont pu etre forcés dans Summit
/// Le traitement se fait sur chacune des 2 jambes du spreadLeg
void GetSpreadCustomFixings(MSXML2::IXMLDOMNode* node, 
							ARM_Vector** startDates1, 
							ARM_Vector** resetDates1, 
							ARM_Vector** fixings1, 
							ARM_Vector** startDates2,
							ARM_Vector** resetDates2, 
							ARM_Vector** fixings2)
{
	MSXML2::IXMLDOMNode* theNode=NULL;
	MSXML2::IXMLDOMNode* theNode2=NULL;
	MSXML2::IXMLDOMNodeList* resultList = NULL;
	BSTR resultat = NULL;
	HRESULT hr;

	VECTOR<double> VStartDates1;
	VECTOR<double> VStartDates2;
	VECTOR<double> VResetDates1;
	VECTOR<double> VResetDates2;
	VECTOR<double> VFixings1;
	VECTOR<double> VFixings2;

	// Recuperation des evenements							{
	if( node->selectNodes(_bstr_t("Events/EVENT"), &resultList)== S_OK)
	{
		if (resultList!=NULL)
		{
			long nb;
			resultList->get_length(&nb);

			int firstROF = 1;
			for (long i=0 ; i<nb ; i++)
			{	
				hr=resultList->get_item(i, &theNode);
				if (hr==S_OK && theNode!=NULL)
				{
					theNode->selectSingleNode(_bstr_t("Type"), &theNode2);
					
					if (theNode2!=NULL)
					{
						theNode2->get_text(&resultat);
						_bstr_t ff(resultat,false);
						char * ff1=(char *)ff;

						if ( strcmp(ff1,"FRC")==0 ) // Jambe 1
						{
							ARM_Date tmpDate = GetDateFromXMLNode(theNode, "ADate");
							if (tmpDate != ARM_DEFAULT_DATE)
							{
								VStartDates1.push_back(GetDateFromXMLNode(theNode, "Date").GetJulian());
								VResetDates1.push_back(tmpDate.GetJulian());
								VFixings1.push_back(100.0*GetDoubleFromXMLNode(theNode,"Amount"));
							}
						}
						if ( strcmp(ff1,"NRC")==0 ) // Jambe 2
						{
							ARM_Date tmpDate = GetDateFromXMLNode(theNode, "ADate");
							if (tmpDate != ARM_DEFAULT_DATE)
							{
								VStartDates2.push_back(GetDateFromXMLNode(theNode, "Date").GetJulian());
								VResetDates2.push_back(tmpDate.GetJulian());
								VFixings2.push_back(100.0*GetDoubleFromXMLNode(theNode,"Amount"));
							}
						}

						if (resultat) SysFreeString(resultat);
						theNode2->Release();
						theNode2=NULL;
					}
					theNode->Release();
					theNode=NULL;
				}
			}
			if (resultList) resultList->Release();
			resultList=NULL;
		}
	}

	*startDates1 = CreateARMVectorFromVECTOR(VStartDates1);
	*resetDates1 = CreateARMVectorFromVECTOR(VResetDates1);
	*fixings1 = CreateARMVectorFromVECTOR(VFixings1);
	*startDates2 = CreateARMVectorFromVECTOR(VStartDates2);
	*resetDates2 = CreateARMVectorFromVECTOR(VResetDates2);
	*fixings2 = CreateARMVectorFromVECTOR(VFixings2);
}


// Convertion

ARM_Curve* FromRefValueToARM_Curve(const ARM_ReferenceValue* refVal, const ARM_Date& date, ARM_Interpolator<double,double>* interpolator)
{
	int size = 0;

	ARM_Vector*    dates    = NULL;
	ARM_Vector*    values   = NULL;
	std::vector<double>& gpDates  = NULL;
	std::vector<double>& gpValues = NULL;
	ARM_Curve*     curve    = NULL;

	dates = refVal->GetDiscreteDates();
	if (dates != NULL)
		size = dates->GetSize();

	if (size == 0)
	{
		gpDates = new std::vector<double>(1, 0.0);
	}
	else
	{
		gpDates  = new std::vector<double>(dates->GetSize(),dates->GetElt());
		*gpDates -= date.GetJulian();
	}

	values   = refVal->GetDiscreteValues();
	gpValues = new std::vector<double>(values->GetSize(),values->GetElt());

	// steptUpLeft car les strikes sont donnés avec les start dates
	curve = new ARM_Curve(*gpDates,*gpValues, interpolator);

	if (gpDates)
	   delete gpDates;
	gpDates = NULL;

	if (gpValues)
	   delete gpValues;
	gpValues = NULL;

	return curve;
}


//TARN Calculator: Get customized reset dates.
std::vector<double>& GetTarnCustomDates(MSXML2::IXMLDOMNode* node)
{
	MSXML2::IXMLDOMNode        *theNode    = NULL;
	MSXML2::IXMLDOMNode        *theNode2   = NULL;
	MSXML2::IXMLDOMNode        *listItem   = NULL;
	MSXML2::IXMLDOMNodeList    *resultList = NULL;
	BSTR               resultat    = NULL;
	HRESULT            hr;

	std::vector<double>&     genDates  = NULL;
	std::vector<double>&     gen       = NULL;

	long               nbEvents;
	ARM_Date		   currGenDate;

	// nominaux suivants si non constants
	node->selectNodes(_bstr_t("Events/EVENT"), &resultList);
	if (resultList!=NULL)
	{
		resultList->get_length(&nbEvents);
		CCString type;
		CCString sun;
		int nbDates = 0;
		genDates = new std::vector<double>(nbEvents,0.0);
		for (int i = 0; i < nbEvents; i++)
		{
			hr=resultList->get_item(i, &listItem);
		
			listItem->selectSingleNode(_bstr_t("Type"), &theNode);
			if (theNode!=NULL)
			{
				resultat = NULL;
				theNode->get_text(&resultat);
				_bstr_t ff(resultat,false);
				char* ff1=(char *)ff;
				type = (const char*) ff1;
				theNode->Release();
				theNode=NULL;
				if (resultat) SysFreeString(resultat);
			}

			if ((strcmp(type,"ROF") == 0))
			{
				listItem->selectSingleNode(_bstr_t("ADate"), &theNode2);
			
				if (theNode2!=NULL)
				{
					resultat = NULL;
					theNode2->get_text(&resultat);
					_bstr_t ff(resultat,false);
					char* ff1=(char *)ff;
					sun = (const char*)ff1;
					if((strcmp(sun, "SUN") !=0))
					{	
						nbDates++;
						currGenDate = ARM_Date(ff1,"YYYYMMDD");
						genDates->Elt(i) = currGenDate.GetJulian();
						theNode2->Release();
						theNode2=NULL;
						if (resultat) SysFreeString(resultat);
					}
				}
			}
			if (listItem)
				listItem->Release();
		}

		//Fullfill the result
		if (nbDates != 0)
		{
			gen = new std::vector<double>(nbDates, 0.0);
			int loop = 0;
			double currVal=0;
			for (i=0; i<nbEvents; i++)
			{
				currVal = genDates->Elt(i);
				if (currVal != 0.0)
				{
					gen->Elt(loop) = currVal;
					loop++;
				}
			}
		}

		if (resultList)
			resultList->Release();
	
	delete genDates;
	genDates = NULL;
	}

	return gen;
}

//TARN Calculator: Get past fixings.
std::vector<double>& GetPastFixings(MSXML2::IXMLDOMNode* node, ARM_Date asOfDate)
{
	MSXML2::IXMLDOMNode        *theNode    = NULL;
	MSXML2::IXMLDOMNode        *theNode2   = NULL;
	MSXML2::IXMLDOMNode        *theNode3   = NULL;
	MSXML2::IXMLDOMNode        *listItem   = NULL;
	MSXML2::IXMLDOMNodeList    *resultList = NULL;
	BSTR               resultat    = NULL;
	HRESULT            hr;

	std::vector<double>&     pastFixings = NULL;
	std::vector<double>&	   fixings	   = NULL;
	long               nbEvents;
	double			   currFixing;
	ARM_Date		   currFixDate;

	node->selectNodes(_bstr_t("Events/EVENT"), &resultList);
	if (resultList!=NULL)
	{
		resultList->get_length(&nbEvents);
		CCString type;
		int nbFixing = 0;
		pastFixings = new std::vector<double>(nbEvents,0.0);
		for (int i = 0; i < nbEvents; i++)
		{
			hr=resultList->get_item(i, &listItem);
			listItem->selectSingleNode(_bstr_t("Type"), &theNode);
			if (theNode!=NULL)
			{
				resultat = NULL;
				theNode->get_text(&resultat);
				_bstr_t ff(resultat,false);
				char* ff1=(char *)ff;
				type = (const char*) ff1;
				theNode->Release();
				theNode=NULL;
				if (resultat) SysFreeString(resultat);
			}

			if ((strcmp(type,"FRC") == 0))
			{
				listItem->selectSingleNode(_bstr_t("ADate"), &theNode3);
				if (theNode3 != NULL)
				{
					resultat = NULL;
					theNode3->get_text(&resultat);
					_bstr_t ff(resultat,false);
					char* ff1=(char *)ff;
					currFixDate = ARM_Date(ff1,"YYYYMMDD");
					
					//Fixing Date <= AsOfDate
					if (currFixDate.GetJulian() <= asOfDate.GetJulian())
					{
						nbFixing++;
						listItem->selectSingleNode(_bstr_t("Amount"), &theNode2);
						if (theNode2!=NULL)
						{
							resultat = NULL;
							theNode2->get_text(&resultat);
							_bstr_t ff(resultat,false);
							char* ff1=(char *)ff;
							currFixing = atof(ff1);
							pastFixings->Elt(i) = currFixing;
							theNode2->Release();
							theNode2=NULL;
							if (resultat) SysFreeString(resultat);
						}
					}

					theNode3->Release();
					theNode3=NULL;
					if (resultat) SysFreeString(resultat);
				}
			}
			if (listItem)
				listItem->Release();
		}

		//Fullfill the result
		if (nbFixing != 0)
		{
			fixings = new std::vector<double>(nbFixing, 0.0);
			int loop = 0;
			double currVal=0;
			for (i=0; i<nbEvents; i++)
			{
				currVal = pastFixings->Elt(i);
				if (currVal != 0.0)
				{
					fixings->Elt(loop) = currVal;
					loop++;
				}
			}
		}
		if (resultList)
			resultList->Release();
	
	delete pastFixings;
	pastFixings = NULL;
	}

	return fixings;
}

//TARN Calculator: Get past fixings.
std::vector<double>& GetPastStart(MSXML2::IXMLDOMNode* node, ARM_Date asOfDate)
{
	MSXML2::IXMLDOMNode        *theNode    = NULL;
	MSXML2::IXMLDOMNode        *theNode2   = NULL;
	MSXML2::IXMLDOMNode        *theNode3   = NULL;
	MSXML2::IXMLDOMNode        *listItem   = NULL;
	MSXML2::IXMLDOMNodeList    *resultList = NULL;
	BSTR               resultat    = NULL;
	HRESULT            hr;

	std::vector<double>&     pastStart = NULL;
	std::vector<double>&	   start     = NULL;
	long               nbEvents;
	ARM_Date		   currStartDate;
	ARM_Date		   currFixDate;

	node->selectNodes(_bstr_t("Events/EVENT"), &resultList);
	if (resultList!=NULL)
	{
		resultList->get_length(&nbEvents);
		CCString type;
		int nbStart = 0;
		pastStart = new std::vector<double>(nbEvents,0.0);
		for (int i = 0; i < nbEvents; i++)
		{
			hr=resultList->get_item(i, &listItem);
			listItem->selectSingleNode(_bstr_t("Type"), &theNode);
			if (theNode != NULL)
			{
				resultat = NULL;
				theNode->get_text(&resultat);
				_bstr_t ff(resultat,false);
				char* ff1=(char *)ff;
				type = (const char*) ff1;
				theNode->Release();
				theNode = NULL;
				if (resultat) SysFreeString(resultat);
			}

			if ((strcmp(type,"FRC") == 0))
			{
				listItem->selectSingleNode(_bstr_t("ADate"), &theNode2);
				if (theNode2 != NULL)
				{
					resultat = NULL;
					theNode2->get_text(&resultat);
					_bstr_t ff(resultat,false);
					char* ff1=(char *)ff;
					currFixDate = ARM_Date(ff1,"YYYYMMDD");
					
					theNode2->Release();
					theNode2 = NULL;
					if (resultat) SysFreeString(resultat);
					         
					//Fixing Date <= AsOfDate
					if (currFixDate.GetJulian() <= asOfDate.GetJulian())
					{
						nbStart++;
						//We select the start date
						listItem->selectSingleNode(_bstr_t("Date"), &theNode3);
						if (theNode3 != NULL)
						{
							resultat = NULL;
							theNode3->get_text(&resultat);
							_bstr_t ff(resultat,false);
							char* ff1=(char *)ff;
							currStartDate = ARM_Date(ff1,"YYYYMMDD");
							pastStart->Elt(i) = currStartDate.GetJulian();
							
							theNode3->Release();
							theNode3 = NULL;
							if (resultat) SysFreeString(resultat);
						}
					}	
				}
			}
			if (listItem)
				listItem->Release();
		}

		//Fullfill the result
		if (nbStart != 0)
		{
			start = new std::vector<double>(nbStart, 0.0);
			int loop = 0;
			double currVal=0;
			for (i=0; i<nbEvents; i++)
			{
				currVal = pastStart->Elt(i);
				if (currVal != 0.0)
				{
					start->Elt(loop) = currVal;
					loop++;
				}
			}
		}
		if (resultList)
			resultList->Release();
	
	delete pastStart;
	pastStart = NULL;
	}

	return start;
}


string GetTarnNodeName(MSXML2::IXMLDOMNode* node, string nodeName)
{

	MSXML2::IXMLDOMNode* theNode = NULL;
	string sname;
	node->selectSingleNode(_bstr_t(nodeName.c_str()), &theNode);
	if (theNode!=NULL)
	{
		BSTR resultat = NULL;
		theNode->get_text(&resultat);
		_bstr_t ff(resultat,false);
		sname = string((char *)ff);
		theNode->Release();
		theNode=NULL;
		if (resultat) SysFreeString(resultat);
	}

	return sname;
}

void GetFormulaForCredit(MSXML2::IXMLDOMNode* node, char* formula)
{
	//SPREADOPTIONLOGFLTDIGITAL ( [2yEUReurib], [10yEUReurib], [10yEUReurib] )
	
	MSXML2::IXMLDOMNode* theNode=NULL;
	BSTR resultat = NULL;

	node->selectSingleNode(_bstr_t("//CREDSWAP/Formula"), &theNode);
	if (theNode!=NULL)
	{
		BSTR resultat = NULL;
		theNode->get_text(&resultat);

		_bstr_t ff(resultat,false);
		char * ff1=(char *)ff;

		// on supprime les espaces dans la formule
		char* tmp = NULL;
		tmp = filterBlank(ff1);
		strcpy(ff1,tmp);
		if (tmp)
			delete [] tmp;

		//int c ='(';
		//char* tmp2= strchr(ff1, c);
		//if (tmp2!=NULL)
		//	strncpy(formula,ff1,(int) (tmp2-ff1));
		//	formula[(int) (tmp2-ff1)] = '\0';
		strcpy(formula,ff1);

		theNode->Release();
		theNode=NULL;
	}
	if (resultat) SysFreeString(resultat);
	strcpy(formula, strupr(formula));
}

void ParseSpreadFollowTxt(MSXML2::IXMLDOMNode* node, 
						  double& spread,
						  string& ccy,
						  string& Term,
						  string& IndexName,
						  bool& spreadonly)
{	
	MSXML2::IXMLDOMNode* theNode=NULL;
	BSTR resultat = NULL;
	//ccy=NULL;
	spread = 0.;

	char Formula[255];
	GetFormulaForCredit(node,Formula);

	string s_formula = (string)Formula;
	
	int val = 0;

	if ((val=s_formula.find("["))<0) 
	{		
		spreadonly=true;
		spread = atof((char*)s_formula.c_str());
		return;
	}

	int begin = s_formula.find("[");
	int end = s_formula.find("]");
	string formula_spread_begin(s_formula,0,begin-1);
	string formula_spread_end(s_formula,end+1,s_formula.size());
	spread += atof((char*)formula_spread_begin.c_str());
	spread += atof((char*)formula_spread_end.c_str());

	string real_formula(s_formula,begin+1,end-begin-1);

	//ConvSummitIndexToARMIndexTerm((char*)real_formula.c_str(), 
	//							  ccy, armIndexTerm);
	begin = real_formula.find("EUR");
	ccy="EUR";

	if (begin<0) { begin = real_formula.find("USD");ccy="USD";}
	if (begin<0) { begin = real_formula.find("JPY");ccy="JPY";}
	if (begin<0) { begin = real_formula.find("CHF");ccy="CHF";}
	if (begin<0) { begin = real_formula.find("GBP");ccy="GBP";}
	if (begin<0) { begin = real_formula.find("AUD");ccy="AUD";}
	if (begin<0) { begin = real_formula.find("NOK");ccy="NOK";}
	if (begin<0) { begin = real_formula.find("SEK");ccy="SEK";}
	if (begin<0) { begin = real_formula.find("DKK");ccy="DKK";}
	if (begin<0) { begin = real_formula.find("CZK");ccy="CZK";}

	Term = string(real_formula,0,begin);
	IndexName = string(real_formula,begin+3,real_formula.size()-1);
}

////////////////////////////////////////////////////////////////////////
// CAPTION CALCULATOR DATA
////////////////////////////////////////////////////////////////////////
int GetCaptionIndex(MSXML2::IXMLDOMNode* node, ARM_Currency* ccy)
{
	MSXML2::IXMLDOMNode* theNode=NULL;
	char summitIndex[15];
	char armIndex[15];

	node->selectSingleNode(_bstr_t("ProdData/cCAPT_I/Index"), &theNode);
	if (theNode!=NULL)
	{
		BSTR resultat = NULL;
		theNode->get_text(&resultat);

		_bstr_t ff(resultat, false);
		char* ff1 = (char *)ff;
		strcpy(summitIndex, ff1);
		
		if (resultat) 
			SysFreeString(resultat);
		theNode->Release();
		theNode = NULL;
	}
		
	if ( (strcmp(summitIndex,"") == 0) || (summitIndex[0] != '[') )
		return K_FIXED;
	else
	{
		char tmp[20];
		ConvSummitIndexToARMIndex(summitIndex, ccy, tmp);
		strcpy(armIndex,tmp);

		return ARM_ConvIrType(armIndex);
	}
}

int GetCaptionNotifDays(MSXML2::IXMLDOMNode* node)
{
	MSXML2::IXMLDOMNode* theNode = NULL;
	int value = 0;

	node->selectSingleNode(_bstr_t("ProdData/cCAPT_I/NotifDays"), &theNode);
	if (theNode!=NULL)
	{
		BSTR resultat = NULL;
		theNode->get_text(&resultat);

		_bstr_t ff(resultat, false);
		char* ff1 = (char*)ff;
		value = atoi((const char*) ff1);
		
		if (resultat) 
			SysFreeString(resultat);
		theNode->Release();
		theNode = NULL;
	}
	return value;
}

int GetCaptionUnderPorS(MSXML2::IXMLDOMNode* node)
{
	MSXML2::IXMLDOMNode* theNode = NULL;
	int value = K_PAY;

	node->selectSingleNode(_bstr_t("ProdData/cCAPT_I/UnderPorS"), &theNode);
	if (theNode!=NULL)
	{
		BSTR resultat = NULL;
		theNode->get_text(&resultat);
		_bstr_t ff(resultat, false);
		char* ff1 = (char *)ff;

		if (strcmp((const char*) ff1, "S") == 0)
			value = K_RCV;

		if (resultat) 
			SysFreeString(resultat);
		theNode->Release();
		theNode = NULL;
	}
	return value;
}

ARM_Date GetCaptionLockout(MSXML2::IXMLDOMNode* node)
{
	// corresponds to first expiry
	MSXML2::IXMLDOMNode* theNode = NULL;
	ARM_Date resDate(ARM_DEFAULT_DATE);

	node->selectSingleNode(_bstr_t("ProdData/cCAPT_I/Lockout"), &theNode);
	if (theNode!=NULL)
	{
		BSTR resultat = NULL;
		theNode->get_text(&resultat);
		_bstr_t ff(resultat, false);
		char* ff1 = (char *)ff;

		if ((strcmp(ff1,"") != 0) && (strcmp(ff1,"SUN") != 0))
			resDate = ARM_Date(ff1,"YYYYMMDD");

		if (resultat) 
			SysFreeString(resultat);
		theNode->Release();
		theNode = NULL;
	}
	return resDate;
}

ARM_Date GetCaptionFinalExpiry(MSXML2::IXMLDOMNode* node)
{
	MSXML2::IXMLDOMNode* theNode = NULL;
	ARM_Date resDate(ARM_DEFAULT_DATE);

	node->selectSingleNode(_bstr_t("ProdData/cCAPT_I/FinalExpiry"), &theNode);
	if (theNode!=NULL)
	{
		BSTR resultat = NULL;
		theNode->get_text(&resultat);
		_bstr_t ff(resultat, false);
		char* ff1 = (char *)ff;

		if ((strcmp(ff1,"") != 0) && (strcmp(ff1,"SUN") != 0))
			resDate = ARM_Date(ff1,"YYYYMMDD");

		if (resultat) 
			SysFreeString(resultat);
		theNode->Release();
		theNode = NULL;
	}
	return resDate;
}

ARM_Vector* GetCaptionExerciseDates(MSXML2::IXMLDOMNode* node, 
									ARM_Date AsOf, 
									ARM_INDEX_TYPE indexType,
									int resetTiming,
									int payTiming,
									ARM_Currency* ccy,
									int notifDays,
									char* resetCal,
									char* payCal)
{
	// build a leg using:
	// - start date,
	// - last expry,
	// - notification days (as reset gap)
	// and return the RESET DATES generated by ARM_SwapLeg

	ARM_Date startDate = GetStartDate(node);
	// end date = last expiry + 1 period :
	ARM_Date endDate = GetCaptionFinalExpiry(node);
	int freq = K_DEF_FREQ;

	MSXML2::IXMLDOMNode* theNode = NULL;
	node->selectSingleNode(_bstr_t("ProdData/cCAPT_I/ExerciseFreq"), &theNode);
	if (theNode!=NULL)
	{
		BSTR resultat = NULL;
		theNode->get_text(&resultat);
		_bstr_t ff(resultat, false);
		char* ff1 = (char *)ff;
		freq = FromSummitFreqToARMFreq(ff1);
		endDate = endDate.AddPeriod(freq, ccy->GetCcyName());
	}
	
	ARM_SwapLeg* temp = new ARM_SwapLeg(startDate, endDate, indexType,
										K_RCV, 0.0, freq, freq, 
										resetTiming, payTiming,
										ccy, K_ADJUSTED, notifDays, 
										resetCal, payCal);

	// NB : take ever date, do not filter forward dates
	ARM_Vector* ExerciseDates = new ARM_Vector();
	for (int i=0; i<temp->GetResetDates()->GetSize(); i++)
	{
		ExerciseDates->push_back(temp->GetResetDates()->Elt(i));
	}

	if (temp)
		delete temp;

	return ExerciseDates;
}

int GetCaptionCpnDayCount(MSXML2::IXMLDOMNode* node)
{
	int value = KNOBASE;

	MSXML2::IXMLDOMNode* theNode = NULL;
	node->selectSingleNode(_bstr_t("ProdData/cCAPT_I/BasisBSpread"), &theNode);
	if (theNode!=NULL)
	{
		BSTR resultat = NULL;
		theNode->get_text(&resultat);
		_bstr_t ff(resultat, false);
		char* ff1 = (char *)ff;

		if (strcmp((const char*) ff1,"") != 0)
			value = FromSummitDaycountToARMDaycount((const char*) ff1);

		if (resultat) 
			SysFreeString(resultat);
		theNode->Release();
		theNode = NULL;
	}

	return value;
}

ARM_ReferenceValue* GetCaptionCpnSpread(MSXML2::IXMLDOMNode* node)
{
	MSXML2::IXMLDOMNode* theNode = NULL;
	node->selectSingleNode(_bstr_t("ProdData/cCAPT_I/FloatBSpread"), &theNode);
	ARM_ReferenceValue* spreads = NULL;
	double spread = 0.0;

	if (theNode!=NULL)
	{
		BSTR resultat = NULL;
		theNode->get_text(&resultat);
		_bstr_t ff(resultat, false);
		char* ff1 = (char *)ff;

		if ((strcmp(ff1, "") != 0) && (strcmp(ff1, "SUN") != 0))
			spread = atof(ff1)/100.;

		spreads = new ARM_ReferenceValue(spread);
		
		if (resultat) 
			SysFreeString(resultat);
		theNode->Release();
		theNode = NULL;
	}

	return spreads;
}



ARM_ReferenceValue* GetCaptionFundSpread(MSXML2::IXMLDOMNode* node)
{
	MSXML2::IXMLDOMNode* theNode = NULL;
	
    node->selectSingleNode(_bstr_t("ProdData/cCAPT_I/FloatASpread"), &theNode);
	
    ARM_ReferenceValue* spreads = NULL;

	if ( theNode != NULL )
	{
		BSTR resultat = NULL;
		
        theNode->get_text(&resultat);
		_bstr_t ff(resultat, false);
		
        char* ff1 = (char *)ff;

		if (strcmp(ff1, ""))
		{
			double spread = atof(ff1)/100.;

            if ( spread == 0.0 )
            {
               return(NULL);
            }

			spreads = new ARM_ReferenceValue(spread);
		}

		if (resultat) 
			SysFreeString(resultat);

		theNode->Release();
		theNode = NULL;
	}

	return spreads;
}



int GetCaptionFundDayCount(MSXML2::IXMLDOMNode* node)
{
	int value = KNOBASE;

	MSXML2::IXMLDOMNode* theNode = NULL;
	node->selectSingleNode(_bstr_t("ProdData/cCAPT_I/BasisASpread"), &theNode);
	if (theNode!=NULL)
	{
		BSTR resultat = NULL;
		theNode->get_text(&resultat);
		_bstr_t ff(resultat, false);
		char* ff1 = (char *)ff;

		if (strcmp((const char*) ff1,"") != 0)
			value = FromSummitDaycountToARMDaycount((const char*) ff1);

		if (resultat) 
			SysFreeString(resultat);
		theNode->Release();
		theNode = NULL;
	}

	return value;
}

int GetGlobalCapIndexType(MSXML2::IXMLDOMNode* node, ARM_Currency* ccy)
{
	MSXML2::IXMLDOMNode* theNode=NULL;
	char summitIndex[15];
	char armIndex[15];

	node->selectSingleNode(_bstr_t("ProdData/cCRFARM_I/Index"), &theNode);
	if (theNode!=NULL)
	{
		BSTR resultat = NULL;
		theNode->get_text(&resultat);

		_bstr_t ff(resultat, false);
		char* ff1 = (char *)ff;
		strcpy(summitIndex, ff1);
		
		if (resultat) 
			SysFreeString(resultat);
		theNode->Release();
		theNode = NULL;
	}
		
	if ( (strcmp(summitIndex,"") == 0) || (summitIndex[0] != '[') )
		return K_FIXED;
	else
	{
		char tmp[20];
		ConvSummitIndexToARMIndex(summitIndex, ccy, tmp);
		strcpy(armIndex,tmp);

		return ARM_ConvIrType(armIndex);
	}
}

void GetGlobalCapInfo(	MSXML2::IXMLDOMNode* node,
						ARM_ReferenceValue*& spreads,
						ARM_ReferenceValue*& fixedRates,
						ARM_ReferenceValue*& barriers,
						ARM_SwapLeg* leg,
						double& numerator,
						double& denominator,
						double& nbIter,
						ARM_Date asOf)
{
	MSXML2::IXMLDOMNode			*theNode	= NULL;
	MSXML2::IXMLDOMNode			*theNode2	= NULL;
	MSXML2::IXMLDOMNode			*theNode3	= NULL;
	MSXML2::IXMLDOMNode			*theNode4	= NULL;
	MSXML2::IXMLDOMNodeList		*resultList = NULL;
	HRESULT				hr = S_OK;
	BSTR				resultat	= NULL;
	long				nb;
	VECTOR<double>		vFixedRates;
	VECTOR<double>		vSpreads;
	VECTOR<double>		vBarriers;
	ARM_Vector			*myFixedRatesVector = NULL;
	ARM_Vector			*mySpreadsVector	= NULL;
	ARM_Vector			*myBarriersVector	= NULL;

	node->selectNodes(_bstr_t("ProdData/cCRFARM_I/Schedule/cCRFLST_I"), &resultList);
	
	if (resultList!=NULL)
	{
		resultList->get_length(&nb);
		CCString type;
		
		if (nb==0)
		{
			hr = S_FALSE;

			CCString msg((CCString)"Invalid XML string for getting Global Cap Info \n");

			throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
					 (char*) msg);
		}
		
		for (int i = 0;i < nb; i++)
		{
			hr=resultList->get_item(i, &theNode);
			
			if (hr==S_OK && theNode!=NULL)
			{
				theNode->selectSingleNode(_bstr_t("FixedCoupon"),&theNode2);
				theNode->selectSingleNode(_bstr_t("CapValue"),&theNode3);
				theNode->selectSingleNode(_bstr_t("FloorValue"),&theNode4);
				
				if (theNode2!=NULL)
				{
					resultat = NULL;
					theNode2->get_text(&resultat);
					_bstr_t ff(resultat,false);
					char * ff1=(char *)ff;
					vFixedRates.push_back(atof(ff1));

					if (resultat) SysFreeString(resultat);
					theNode2->Release();
					theNode2=NULL;
				}
				if (theNode3!=NULL)
				{
					resultat = NULL;
					theNode3->get_text(&resultat);
					_bstr_t ff(resultat,false);
					char * ff1=(char *)ff;
					vSpreads.push_back(atof(ff1));

					if (resultat) SysFreeString(resultat);
					theNode3->Release();
					theNode3=NULL;
				}
				if (theNode4!=NULL)
				{
					resultat = NULL;
					theNode4->get_text(&resultat);
					_bstr_t ff(resultat,false);
					char * ff1=(char *)ff;
					vBarriers.push_back(atof(ff1));

					if (resultat) SysFreeString(resultat);
					theNode4->Release();
					theNode4=NULL;
				}
			}
			if (theNode) theNode->Release();
			theNode = NULL;
		}

		if (resultList) resultList->Release();
		resultList=NULL;
		
		if (nb==1)
		{
			fixedRates = new ARM_ReferenceValue(fixedRates[0]);
			spreads = new ARM_ReferenceValue(spreads[0]);
			barriers = new ARM_ReferenceValue(barriers[0]);
		}
		else
		{
			ARM_Vector* myDatesVector= leg->GetResetDates();
			ARM_Vector* myFixedRatesVector= CreateARMVectorFromVECTOR(vFixedRates);
			ARM_Vector* mySpreadsVector= CreateARMVectorFromVECTOR(vSpreads);
			ARM_Vector* myBarriersVector= CreateARMVectorFromVECTOR(vBarriers);
			
			fixedRates = new ARM_ReferenceValue((ARM_Vector*) myDatesVector->Clone(), myFixedRatesVector);
			fixedRates->SetCalcMethod(K_STEPUP_LEFT);
			spreads = new ARM_ReferenceValue((ARM_Vector*) myDatesVector->Clone(), mySpreadsVector);
			spreads->SetCalcMethod(K_STEPUP_LEFT);
			barriers = new ARM_ReferenceValue((ARM_Vector*) myDatesVector->Clone(), myBarriersVector);
			barriers->SetCalcMethod(K_STEPUP_LEFT);
		}
		
		node->selectSingleNode(_bstr_t("ProdData/cCRFARM_I/FloorValue"), &theNode);
		if (theNode!=NULL)
		{
			BSTR resultat = NULL;
			theNode->get_text(&resultat);

			_bstr_t ff(resultat, false);
			char* ff1 = (char *)ff;
			denominator = atof(ff1);
			
			if (resultat) 
			SysFreeString(resultat);
			theNode->Release();
			theNode = NULL;
		}
		
		node->selectSingleNode(_bstr_t("ProdData/cCRFARM_I/CapValue"), &theNode);
		if (theNode!=NULL)
		{
			BSTR resultat = NULL;
			theNode->get_text(&resultat);

			_bstr_t ff(resultat, false);
			char* ff1 = (char *)ff;
			numerator = atof(ff1);
			
			if (resultat) 
			SysFreeString(resultat);
			theNode->Release();
			theNode = NULL;
		}
	}
}


ARM_ReferenceValue* GetFxFixings(MSXML2::IXMLDOMNode* node)
{
	MSXML2::IXMLDOMNode* theNode = NULL;
	MSXML2::IXMLDOMNode* listItem = NULL;
	MSXML2::IXMLDOMNodeList* resultList = NULL;
	ARM_ReferenceValue* newRefVal = NULL;
	BSTR resultat = NULL;
	HRESULT hr;
	long nbEvents;
	
	VECTOR<double> fxFixingDate;
	VECTOR<double> fxFixingVal;

	if (node->selectNodes(_bstr_t("Events/EVENT"), &resultList) == S_OK)
	{
		resultList->get_length(&nbEvents);

		CCString type;

		for(int i = 0; i < nbEvents; i++)
		{
			hr = resultList->get_item(i, &listItem);
			
			listItem->selectSingleNode(_bstr_t((const char *)"Type"), &theNode);
			
			if(theNode != NULL)
			{
				BSTR resultat = NULL;
				theNode->get_text(&resultat);

				_bstr_t ff(resultat,false);
				char * ff1=(char *)ff;

				type = (const char*) ff1;

				theNode->Release();
				theNode = NULL;
				if(resultat) 
					SysFreeString(resultat);
			}

			if ( strcmp(type,"FX") == 0 )
			{
				listItem->selectSingleNode(_bstr_t((const char*)"Amount"), &theNode);
				if (theNode != NULL)
				{
					BSTR resultat = NULL;
					theNode->get_text(&resultat);

					_bstr_t ff(resultat,false);
					char * ff1=(char *)ff;

					double dfxfixing = atof(ff1);
					fxFixingVal.push_back(dfxfixing);

					theNode->Release();
					theNode = NULL;
					if (resultat) 
						SysFreeString(resultat);
				}

				fxFixingDate.push_back(GetDateFromXMLNode(listItem,"ADate").GetJulian());
			
			}
		}

		if (fxFixingDate.size() > 0)
		{
			ARM_Vector*	vFxFixingDates = CreateARMVectorFromVECTOR(fxFixingDate);
			ARM_Vector* vFxFixingVal = CreateARMVectorFromVECTOR(fxFixingVal);
			newRefVal = new ARM_ReferenceValue(vFxFixingDates, vFxFixingVal);
			newRefVal->SetCalcMethod(K_PERFECT_DISCRETE_REF);
			newRefVal->SetExtrapolMeth(1);
		}
	}

	return	newRefVal;
}


string calypso2ARMCalendar(string& cal){
	string ARM_cal;
	int len = cal.size();

	for(int i =0;i<len;i+=4){
	 ARM_cal +=cal.substr(i,3);
	}

	return ARM_cal;
}

int GetCalypsoPorS(MSXML2::IXMLDOMNodePtr& node)
{
	string proS; XMLTools::convert(XMLTools::selectSingleNode(node,"/BuyOrSell"),proS);
	return (strcmp(proS.c_str(),"Sell")==0)?K_PAY:K_RCV;
}


void parseFormula(string formula,double& strike,double& leverage,string& index) {
		
		strike =0;
		leverage =1;
		index="";

		if(formula.find_first_not_of(' ') != -1){
			int minusIndex = formula.find((int) '-');
		
			string strikeStr;
			if (minusIndex > 0) 
				strikeStr = formula.substr(0, minusIndex);
			else 
				strikeStr = formula;
		
		
			strike = atof(strikeStr.c_str());
			if(strike == 0) 
				throw Exception(__LINE__, __FILE__, ERR_INVALID_INPUT,"Invalid strike in  formula: "+formula);
			
			if (minusIndex <= 0)
				return;
		
			string rightSide = formula.substr(minusIndex + 1);
			int timesIndex = rightSide.find((int) '*');
			if (timesIndex < 0) {
				index = rightSide;
			}else{
				string leftFromTimes = rightSide.substr(0, timesIndex - 1);
				string rightFromTimes = rightSide.substr(timesIndex + 1);
			
				leverage = atof(leftFromTimes.c_str());
				index = rightFromTimes;
				if(leverage==0){
					leverage = atof(rightFromTimes.c_str());
					index = leftFromTimes;
					if(leverage ==0)
						throw Exception(__LINE__, __FILE__, ERR_INVALID_INPUT,"Invalid leverage in formula: "+formula);
				} 
			}
		}
}

int GetCompFreq(MSXML2::IXMLDOMNode* node)
{
	MSXML2::IXMLDOMNode* theNode=NULL;
	BSTR resultat = NULL;
	int value = K_ANNUAL;
	
	node->selectSingleNode(_bstr_t("COMP_CompFreq"), &theNode);
	if (theNode!=NULL)
	{
		theNode->get_text(&resultat);
		_bstr_t ff(resultat,false);
		char * ff1=(char *)ff;
		if (strcmp(ff1,"")!=0)
			value = FromSummitFreqToARMFreq((const char*) ff1);

		theNode->Release();
		theNode=NULL;
		if (resultat) SysFreeString(resultat);
	}
	
	return value;
}

int GetCompMode(MSXML2::IXMLDOMNode* node)
{
	MSXML2::IXMLDOMNode* theNode=NULL;
	BSTR resultat = NULL;
	int value = K_COMP_NONE;
	
	node->selectSingleNode(_bstr_t("COMP_Mode"), &theNode);
	if (theNode!=NULL)
	{
		theNode->get_text(&resultat);
		_bstr_t ff(resultat,false);
		char * ff1=(char *)ff;

		if (strcmp(ff1,"INC")==0)
			value = K_SPREAD_INC;

		theNode->Release();
		theNode=NULL;
		if (resultat) SysFreeString(resultat);
	}
	
	return value;
}


ARM_ReferenceValue* GetPayFixedRate(MSXML2::IXMLDOMNode* node)
{
	ARM_ReferenceValue* fixedRates = NULL;

	VECTOR<double> rate;
	ARM_Vector* vRate = NULL;
	ARM_Vector* spread = NULL;
	
	ARM_Vector* resetDates = NULL;

	MSXML2::IXMLDOMNode* schedLineNode=NULL;
	MSXML2::IXMLDOMNodeList* resultList=NULL;
	HRESULT hr;

	node->selectNodes(_bstr_t("ProdData/cBLOB_I/Schedule1/cBLBLST_I"), &resultList);
	if (resultList!=NULL)
	{
		long nb;
		resultList->get_length(&nb);
		if (nb > 0)
		{
			resetDates = new ARM_Vector(nb);
			spread = new ARM_Vector(nb);
			
			for (long i=0 ; i<nb ; i++)
			{	
				hr=resultList->get_item(i, &schedLineNode);
				if (hr==S_OK && schedLineNode!=NULL)
				{
					spread->Elt(i) = GetDoubleFromXMLNode(schedLineNode,"Rate2");
					resetDates->Elt(i) = GetDateFromXMLNode(schedLineNode,"FixingDate").GetJulian();

					double dRate = GetDoubleFromXMLNode(schedLineNode,"Rate1");
					rate.push_back(dRate*100.);
									}
				schedLineNode->Release();
				schedLineNode = NULL;
			}
		}
		resultList->Release();
		resultList = NULL;
	}

	if (rate.size() > 0)
		vRate = CreateARMVectorFromVECTOR(rate);
	
	ARM_Vector* fixRates = new ARM_Vector(resetDates->GetSize());
		for (int i = 0; i < fixRates->GetSize(); i++)
			fixRates->Elt(i) = vRate->Elt(i) + spread->Elt(i);

	fixedRates = new ARM_ReferenceValue((ARM_Vector*)resetDates,fixRates);

	if (spread)
		delete spread;
	spread = NULL;

	return(fixedRates);
}

void GetPayIndexResetCal(MSXML2::IXMLDOMNode* node, char* resetCal)
{
	MSXML2::IXMLDOMNode* theNode=NULL;
	BSTR resultat = NULL;
	
	node->selectSingleNode(_bstr_t("ProdData/cBLOB_I/ResetCal1"), &theNode);
	if (theNode!=NULL)
	{
		BSTR resultat = NULL;
		theNode->get_text(&resultat);

		_bstr_t ff(resultat,false);
		char * ff1=(char *)ff;

		strcpy(resetCal,ff1);

		theNode->Release();
		theNode=NULL;
		if (resultat) SysFreeString(resultat);
	}
}


ARM_ReferenceValue*  createRefValue(vector<double> dates, vector<double> values, int calcMeth, int extrapolMeth, bool checkIfConstant, bool useDateEvenIfConstant) {
	ARM_ReferenceValue* refValue = NULL;
	ARM_Vector* armDates = NULL;
	ARM_Vector* armValues = NULL;

	if(values.size() > 0) {
		double firstValue = values[0];
		bool isConstant = checkIfConstant;
		
		std::vector<double>::iterator it = values.begin();
		while (it != values.end()) {
			if(*it != firstValue) {
				isConstant = false;
				break;
			}
			++it;
		}
		if(isConstant) {
			if(useDateEvenIfConstant) {
				vector<double> newDates;
				vector<double> newValues;
				newDates.push_back(dates[0]);
				newValues.push_back(firstValue);
				armDates = new ARM_Vector(newDates);
				armValues = new ARM_Vector(newValues);
				refValue = new ARM_ReferenceValue((ARM_Vector*)armDates->Clone(), (ARM_Vector*)armValues->Clone());
				refValue->SetCalcMethod(calcMeth);
				refValue->SetExtrapolMeth(extrapolMeth);
			} else {
				refValue = new ARM_ReferenceValue(firstValue);
			}
		} else {
			armDates = new ARM_Vector(dates);
			armValues = new ARM_Vector(values);
			refValue = new ARM_ReferenceValue((ARM_Vector*)armDates->Clone(), (ARM_Vector*)armValues->Clone());
			refValue->SetCalcMethod(calcMeth);
			refValue->SetExtrapolMeth(extrapolMeth);
		}
	}
	if(armDates) {
		delete armDates;
		armDates = NULL;
	}
	if(armValues) {
		delete armValues;
		armValues = NULL;
	}
	return refValue;
}
