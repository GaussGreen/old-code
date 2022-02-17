
#include <ARM\libarm_local\firstToBeIncluded.h>
#include <ARM\libarm_local\arm_local_glob.h>
#include <ARM\libarm_local\arm_local_persistent.h>

#include <ARM\libarm_frometk\arm_local_etoolkit.h>

#include <ARM\libarm_local\undef_va_vars.h>

#include "PaserManagerUtilities.h"
#include "ParserManager.h"
#include "VariantTools.h"

#import <msxml3.dll> raw_interfaces_only
using namespace MSXML2;
using namespace etoolkit;




long GetTradeIdFromXML(const CCString& xmlResponse, CCString& tradeId, CCString& typeDeal)
{

	VARIANT_BOOL bOK;
	HRESULT hr;
	MSXML2::IXMLDOMDocument *XMLDoc = NULL;
	MSXML2::IXMLDOMNode *theNode = NULL;

	BSTR resultat = NULL;

	try
	{
		hr = CoInitialize(NULL); 
//		SUCCEEDED(hr) ? 0 : throw hr;
		hr = CoCreateInstance(__uuidof(MSXML2::DOMDocument30), NULL, CLSCTX_INPROC_SERVER, __uuidof(MSXML2::IXMLDOMDocument), (void**)&XMLDoc);
		SUCCEEDED(hr) ? 0 : throw hr;

		//JLA _bstr_t tmpChaine = xmlResponse;
		_bstr_t tmpChaine; 
		VariantTools::convert(CCSTringToSTLString(xmlResponse),tmpChaine); 


		XMLDoc->loadXML(tmpChaine, &bOK);
	}
	catch(...)
	{
		if (XMLDoc) XMLDoc->Release();
		hr = S_FALSE;

		CCString msg((CCString)"Pb in creating XML document for getting List of deals (" + xmlResponse + ")\n");

		throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
						 (const char *) msg);
	}

	XMLDoc->selectSingleNode(_bstr_t("*/*/TradeId"), &theNode);
	if (theNode!=NULL)
	{
		BSTR resultat = NULL;
		theNode->get_text(&resultat);

		_bstr_t ff(resultat,false);
		char * ff1=(char *)ff;

		tradeId = ff1;

		theNode->Release();
		theNode=NULL;
		if (resultat) SysFreeString(resultat);
	}
	
	ConvertManager::Init();
	typeDeal = ConvertManager::GetETKKey(XMLDoc).c_str();

	return ARM_OK;
}


VECTOR<CCString> ARMLOCAL_GetDealByFilter(const CCString& filter, VECTOR<CCString>& listTypes)
{
	CCString xmlResponse;
	CCString msglist;

	long retCode;

	retCode = etoolkit_getFilter(filter, xmlResponse, msglist);

	VARIANT_BOOL bOK;
	HRESULT hr;
	MSXML2::IXMLDOMDocument *XMLDoc = NULL;

	BSTR resultat = NULL;

	try
	{
		hr = CoInitialize(NULL); 
//		SUCCEEDED(hr) ? 0 : throw hr;
		hr = CoCreateInstance(__uuidof(MSXML2::DOMDocument30), NULL, CLSCTX_INPROC_SERVER, __uuidof(MSXML2::IXMLDOMDocument), (void**)&XMLDoc);
		SUCCEEDED(hr) ? 0 : throw hr;

		//JLA _bstr_t tmpChaine = xmlResponse;
		_bstr_t tmpChaine; 
		VariantTools::convert(CCSTringToSTLString(xmlResponse),tmpChaine); 


		XMLDoc->loadXML(tmpChaine, &bOK);
	}
	catch(...)
	{
		if (XMLDoc) XMLDoc->Release();
		hr = S_FALSE;

		CCString msg((CCString)"Pb in creating XML document for getting List of deals (" + xmlResponse + ")\n");

		throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
						 (const char *) msg);
	}

	MSXML2::IXMLDOMNodeList * resultList = NULL;
	MSXML2::IXMLDOMNode * listItem = NULL, * theNode = NULL;

	CCString tradeListHandle;
	int nbTrades;

	XMLDoc->selectSingleNode((_bstr_t)("/Response/Iterator"), &theNode);
	if (theNode!=NULL)
	{
		theNode->get_text(&resultat);
		tradeListHandle = _com_util::ConvertBSTRToString(resultat);
		if (resultat) SysFreeString(resultat);
		theNode->Release();
		theNode=NULL;
	}

	XMLDoc->selectSingleNode((_bstr_t)("/Response/Length"), &theNode);
	if (theNode!=NULL)
	{
		theNode->get_text(&resultat);
		_bstr_t tmp (resultat,false);
		sscanf((char*)tmp, "%d", &nbTrades);
		if (resultat) SysFreeString(resultat);
		theNode->Release();
		theNode=NULL;
	}

	if (XMLDoc)
	{
		XMLDoc->Release();
		XMLDoc = NULL;
	}

	VECTOR <CCString> listDeals;
	listDeals.clear();

	if (nbTrades == 0)
	{
		return listDeals;
	}

	retCode = etoolkit_iterateTradeList(tradeListHandle,
										"BEGINNING",
										xmlResponse,
										msglist);

	for (int i = 0; i < nbTrades; i++)
	{
		retCode = etoolkit_iterateTradeList(tradeListHandle,
											"NEXT",
											xmlResponse,
											msglist);

		CCString tradeId;
		CCString typeDeal;

		retCode = GetTradeIdFromXML(xmlResponse,tradeId,typeDeal);

		if (retCode == ARM_OK)
		{
			listDeals.push_back(tradeId);
			listTypes.push_back(typeDeal);
		}
		else
		{
			i = nbTrades;
			listDeals.clear();
			listTypes.clear();
		}
	}

	retCode = etoolkit_releaseIterateTradeList(tradeListHandle,
											   xmlResponse,
											   msglist);

	return listDeals;
}
