#include <ARM\libarm_local\firstToBeIncluded.h>

#include <ARM\libarm_frometk\ARM_local_wsetk.h>
#include <ARM\local_xlarm\ARM_local_interglob.h>
#include "VariantTools.h"

#define CHECK_HR(hr) { if (FAILED(hr)) { throw -1; } }

#import <msxml3.dll> raw_interfaces_only


#import <mssoap30.dll>  \
   exclude("IStream", "IErrorInfo", "ISequentialStream", "_LARGE_INTEGER", \
   "_ULARGE_INTEGER", "tagSTATSTG", "_FILETIME")

// using namespace MSSOAPLib30;
using namespace std;



ARM_EtkSoapClient::ARM_EtkSoapClient()
{
      itsURL = _T(armlocal_init->wsetkprod.c_str());
	itsBase = _T("OTC");

	strcpy(itsUserName,"wsotc");

    HRESULT hr = itsSoapClient.CoCreateInstance( __uuidof(MSSOAPLib30::SoapClient30));

	hr = itsSoapClient->MSSoapInit(_bstr_t(itsURL), L"", L"", L"");

}


ARM_EtkSoapClient::ARM_EtkSoapClient(const CCString&  pProd, const CCString&  pBase, const CCString& pUserName)
{
if (strcmp((const char*)pProd, "PROD") == 0)
	{
		itsURL = _T(armlocal_init->wsetkprod.c_str());
	}
	else
	{
		itsURL = _T(armlocal_init->wsetkrec.c_str());
	}

	itsBase = _T((const char*)pBase);
	strcpy(itsUserName,(const char*)pUserName);

	HRESULT hr = itsSoapClient.CoCreateInstance(__uuidof(MSSOAPLib30::SoapClient30));

	hr = itsSoapClient->MSSoapInit(_bstr_t(itsURL), L"", L"", L"");

}


void ARM_EtkSoapClient::Execute(CCString command, CCString xmlRequest, CCString & xmlResponse, CCString & messageList)
{

	
	USES_CONVERSION;
	
	WCHAR *pwcMethodName = L"OneShotExecute";
	DISPID dispidFn = 0;
	HRESULT hr  = itsSoapClient->GetIDsOfNames(IID_NULL, &pwcMethodName, 1, 
	LOCALE_SYSTEM_DEFAULT, &dispidFn);
	CHECK_HR(hr);

	_bstr_t bAppName("ARM");
	_bstr_t bBase (itsBase);
	_bstr_t bLogin(itsUserName);
	_bstr_t bCommand ((const char*)command);
	_bstr_t bRequest ((const char*)xmlRequest);

	unsigned int uArgErr;
	VARIANT varg[5];
	varg[0].vt = VT_BSTR;
	varg[0].bstrVal = bRequest;
	varg[1].vt = VT_BSTR;
	varg[1].bstrVal = bCommand;
	varg[2].vt = VT_BSTR;
	varg[2].bstrVal = bLogin;
	varg[3].vt = VT_BSTR;
	varg[3].bstrVal = bBase;
	varg[4].vt = VT_BSTR;
	varg[4].bstrVal = bAppName;

	DISPPARAMS params;
	params.cArgs = 5;
	params.rgvarg = varg;
	params.cNamedArgs    = 0;
	params.rgdispidNamedArgs = NULL;

	_variant_t result;

	uArgErr = -1;

	EXCEPINFO excepInfo;
	memset(&excepInfo, 0, sizeof(excepInfo));
	hr  = itsSoapClient->Invoke(dispidFn, IID_NULL, LOCALE_SYSTEM_DEFAULT, 
	DISPATCH_METHOD, &params, &result, &excepInfo, &uArgErr);
	CHECK_HR(hr);

	VARIANT_BOOL bOK;
	MSXML2::IXMLDOMDocument *XMLDoc = NULL;
	MSXML2::IXMLDOMNode *node= NULL;
	
	hr = CoCreateInstance(__uuidof(MSXML2::DOMDocument30), NULL, 
			              CLSCTX_INPROC_SERVER, __uuidof(MSXML2::IXMLDOMDocument), (void**)&XMLDoc);

	SUCCEEDED(hr) ? 0 : throw hr;

	XMLDoc->loadXML(result.bstrVal, &bOK);

	if ( XMLDoc->selectSingleNode(_bstr_t("/data/XMLResponse"), &node) == S_OK )
	{
		BSTR resultat = NULL;
		node->get_text(&resultat);

//		_bstr_t ff(resultat,false);
//		char * ff1=(char *)ff;
		//JLA xmlResponse = W2A(resultat);
		std::string toto; 
		VariantTools::convert(resultat,toto); 
		xmlResponse=toto.c_str() ;

		node->Release();
		node=NULL;
		if (resultat) SysFreeString(resultat);
	}

	if (XMLDoc)
		XMLDoc->Release();
	XMLDoc = NULL;

}




VECTOR<string> ARM_EtkSoapClient::GetTradeList(CCString xmlRequest)
{
	USES_CONVERSION;

	VECTOR<string> res;
	


	WCHAR *pwcMethodName = L"OneShotGetTradeList";
	DISPID dispidFn = 0;
	HRESULT hr  = itsSoapClient->GetIDsOfNames(IID_NULL, &pwcMethodName, 1, 
	LOCALE_SYSTEM_DEFAULT, &dispidFn);
	CHECK_HR(hr);

	_bstr_t bAppName("ARM");
	_bstr_t bBase (itsBase);
	_bstr_t bLogin(itsUserName);
	_bstr_t bRequest ((const char*)xmlRequest);
	long bMaxTrade = 0;

	unsigned int uArgErr;
	VARIANT varg[5];
	varg[0].vt = VT_I4;
	varg[0].intVal = bMaxTrade;
	varg[1].vt = VT_BSTR;
	varg[1].bstrVal = bRequest;
	varg[2].vt = VT_BSTR;
	varg[2].bstrVal = bLogin;
	varg[3].vt = VT_BSTR;
	varg[3].bstrVal = bBase;
	varg[4].vt = VT_BSTR;
	varg[4].bstrVal = bAppName;

	DISPPARAMS params;
	params.cArgs = 5;
	params.rgvarg = varg;
	params.cNamedArgs    = 0;
	params.rgdispidNamedArgs = NULL;

	_variant_t result;

	uArgErr = -1;

	EXCEPINFO excepInfo;
	memset(&excepInfo, 0, sizeof(excepInfo));
	hr  = itsSoapClient->Invoke(dispidFn, IID_NULL, LOCALE_SYSTEM_DEFAULT, 
	DISPATCH_METHOD, &params, &result, &excepInfo, &uArgErr);
	CHECK_HR(hr);

	VARIANT_BOOL bOK;
	MSXML2::IXMLDOMDocument *XMLDoc = NULL;
	MSXML2::IXMLDOMNodeList *nodeList = NULL;
	MSXML2::IXMLDOMNode *node= NULL, *node2 = NULL;
		
	hr = CoInitialize(NULL); 
//	SUCCEEDED(hr) ? 0 : throw hr;
	hr = CoCreateInstance(__uuidof(MSXML2::DOMDocument30), NULL, 
			              CLSCTX_INPROC_SERVER, __uuidof(MSXML2::IXMLDOMDocument), (void**)&XMLDoc);
	SUCCEEDED(hr) ? 0 : throw hr;

	XMLDoc->loadXML(result.bstrVal, &bOK);

	long nbNodes;
	CCString xmlResponse;

	if ( XMLDoc->selectSingleNode(_bstr_t("/data/XMLResponse"), &node) == S_OK )
	{
		BSTR resultat = NULL;
		node->get_text(&resultat);

//		xmlResponse = W2A(resultat);

//		_bstr_t ff(resultat,false);

		node->Release();
		node=NULL;

		if (XMLDoc)
			XMLDoc->Release();
		XMLDoc = NULL;

		hr = CoInitialize(NULL); 
//		SUCCEEDED(hr) ? 0 : throw hr;
		hr = CoCreateInstance(__uuidof(MSXML2::DOMDocument30), NULL, 
							  CLSCTX_INPROC_SERVER, __uuidof(MSXML2::IXMLDOMDocument), (void**)&XMLDoc);
		SUCCEEDED(hr) ? 0 : throw hr;

		XMLDoc->loadXML(resultat, &bOK);

		if (XMLDoc->selectNodes(_bstr_t("tradeList/Response"), &nodeList) == S_OK)
		{
			nodeList->get_length(&nbNodes);

			for (int i = 0; i < nbNodes; i++)
			{
				hr=nodeList->get_item(i, &node2);

				BSTR resultat2 = NULL;
				hr=node2->get_xml(&resultat2);

//				_bstr_t ff(resultat2,false);
//				char * ff1=(char *)ff;
				{
					std::string toto; 
					VariantTools::convert(resultat2,toto); 
					res.push_back(toto);
					// JLA res.push_back(W2A(resultat2));
				}


				if (node2)
					node2->Release();
				node2=NULL;
				if (resultat2) SysFreeString(resultat2);
			}
			if (nodeList)
				nodeList->Release();
			nodeList = NULL;
		}
		if (node)
			node->Release();
		node=NULL;

		if (resultat) SysFreeString(resultat);
	}

	if (XMLDoc)
		XMLDoc->Release();
	XMLDoc = NULL;
 
	return res;
}




void parser(char* pChaine,std::string &pResult)
{
	int separator='//';
	
	char* p = strchr (pChaine, separator);
	char *temp=pChaine;
	char strNoeud[10]="";

	char strCopy[500]="";

	pResult.append("<Name>");
	pResult.append("<Val>");
	pResult.append(pChaine);
	pResult.append("</Val>");

	int d = 0,e=1;
	while(p)
	{
		d = strlen (temp) - strlen (p); 
		if(d > 0)
		{
			sprintf(strNoeud,"<P%i>",e);
			pResult.append(strNoeud);
			strncpy (strCopy, temp, d);
			strCopy[d]='\0';
			pResult.append(CCtrim_right(strCopy));

			sprintf(strNoeud,"</P%i>",e);
			pResult.append(strNoeud);

			e++;
		}

		p++;
		temp = p;
		p = strchr (temp, separator);
	}

	if((temp) && (strlen (temp) > 0))
	{
		sprintf(strNoeud,"<P%i>",e);
		pResult.append(strNoeud);
		pResult.append(CCtrim_right(temp));

		sprintf(strNoeud,"</P%i>",e);
		pResult.append(strNoeud);
	}
	pResult.append("</Name>");

}



void parser2(char* pChaine,std::string &pResult)
{
	pResult.append("<Name2><Val>");
	pResult.append(pChaine);
	pResult.append("</Val></Name2>");
}



HRESULT ParseXml(BSTR* pXml,std::string &result)
{
	HRESULT hr;
	MSXML2::IXMLDOMDocument * pXMLDoc=NULL;
	MSXML2::IXMLDOMNode * pXDN=NULL;
	VARIANT_BOOL bOK;
	MSXML2::IXMLDOMNodeList *resultList=NULL;
	
	MSXML2::IXMLDOMNode *listItem=NULL,*theNode=NULL;
	MSXML2::IXMLDOMNamedNodeMap *attributeMap=NULL;
	BSTR valeur=NULL;
	VARIANT varValue;

	try
	{
		hr = CoInitialize(NULL); 
		hr = CoCreateInstance(__uuidof(MSXML2::DOMDocument30), NULL, 
                              CLSCTX_INPROC_SERVER, __uuidof(MSXML2::IXMLDOMDocument), 
                              (void**) &pXMLDoc);
		SUCCEEDED(hr) ? 0 : throw hr;

		if (pXMLDoc->loadXML(*pXml,&bOK)==S_OK)
		{
			if (pXMLDoc->selectNodes((_bstr_t)("/Response/SummitQuery/xml/rs:data/z:row"),&resultList)==S_OK)
			{
				long nb_nodes;
				resultList->get_length(&nb_nodes);

				for (long index_node=0;index_node<nb_nodes;index_node++)
				{
					hr=resultList->get_item(index_node,&listItem);

					hr = listItem->get_attributes(&attributeMap);
					  if(SUCCEEDED(hr) )
					  {
						 result.append("<Commdata>");
						 hr = attributeMap->getNamedItem((_bstr_t)("c1"), &theNode);
						 if(SUCCEEDED(hr) && theNode)
						 {
							theNode->get_nodeValue(&varValue);

							parser((_bstr_t)varValue,result);

							theNode->Release();
							theNode = NULL;
						 }
						 hr = attributeMap->getNamedItem((_bstr_t)("c2"), &theNode);
						 if(SUCCEEDED(hr) && theNode)
						 {
							theNode->get_nodeValue(&varValue);

							parser2((_bstr_t)varValue,result);

							theNode->Release();
							theNode = NULL;
						 }
						 attributeMap->Release();
						 attributeMap = NULL;
						 result.append("</Commdata>");
					  }
				}
			}
		}
		if (attributeMap) attributeMap->Release();
		if (listItem) listItem->Release();
		if (theNode) theNode->Release();
		if (resultList) resultList->Release();
		if (pXDN) pXDN->Release();
		if (pXMLDoc) pXMLDoc->Release();
		hr=S_OK;
	}
	catch(...)
	{
		if (attributeMap) attributeMap->Release();
		if (listItem) listItem->Release();
		if (resultList) resultList->Release();
		if (pXDN) pXDN->Release();
		if (theNode) theNode->Release();
		if (pXMLDoc) pXMLDoc->Release();
		hr=S_FALSE;
	}

	return hr;
}




void ARM_EtkSoapClient::GetCommsetName(CCString name, CCString name2, CCString asOf, CCString type, CCString cvname, CCString & xmlResponse, CCString & messageList)
{

	USES_CONVERSION;
	
	WCHAR *pwcMethodName = L"OneShotExecute";
	DISPID dispidFn = 0;
	HRESULT hr  = itsSoapClient->GetIDsOfNames(IID_NULL, &pwcMethodName, 1, 
	LOCALE_SYSTEM_DEFAULT, &dispidFn);
	CHECK_HR(hr);

	_bstr_t bAppName("ARM");
	_bstr_t bBase (itsBase);
	_bstr_t bLogin(itsUserName);
	_bstr_t bCommand ("s_base:DBQuery");

	CCString request = (CCString)"<Request><SummitSQL>select COMMSET.Name, COMMSET.Name2 from COMMSET where COMMSET.Name like '";
	request = request + name;
	request = request + (CCString)"' and COMMSET.AsOfDate ='";
	request = request + asOf;
	request = request + (CCString)"' and COMMSET.Type = '";
	request = request + type;
	request = request + (CCString)"' and COMMSET.Id='";
	request = request + cvname;
	
	if (strcmp(name2,"") != 0)
	{
		request = request + (CCString)"' and COMMSET.Name2 like '";
		request = request + name2;
	}

	request = request + (CCString)"'</SummitSQL></Request>";

	_bstr_t bRequest ((const char*)request);

	unsigned int uArgErr;
	VARIANT varg[5];
	varg[0].vt = VT_BSTR;
	varg[0].bstrVal = bRequest;
	varg[1].vt = VT_BSTR;
	varg[1].bstrVal = bCommand;
	varg[2].vt = VT_BSTR;
	varg[2].bstrVal = bLogin;
	varg[3].vt = VT_BSTR;
	varg[3].bstrVal = bBase;
	varg[4].vt = VT_BSTR;
	varg[4].bstrVal = bAppName;

	DISPPARAMS params;
	params.cArgs = 5;
	params.rgvarg = varg;
	params.cNamedArgs    = 0;
	params.rgdispidNamedArgs = NULL;

	_variant_t result;

	uArgErr = -1;

	EXCEPINFO excepInfo;
	memset(&excepInfo, 0, sizeof(excepInfo));
	hr  = itsSoapClient->Invoke(dispidFn, IID_NULL, LOCALE_SYSTEM_DEFAULT, 
	DISPATCH_METHOD, &params, &result, &excepInfo, &uArgErr);
	CHECK_HR(hr);

	VARIANT_BOOL bOK;
	MSXML2::IXMLDOMDocument *XMLDoc = NULL;
	MSXML2::IXMLDOMNode *node= NULL;
	
	hr = CoCreateInstance(__uuidof(MSXML2::DOMDocument30), NULL, 
			              CLSCTX_INPROC_SERVER, __uuidof(MSXML2::IXMLDOMDocument), (void**)&XMLDoc);

	SUCCEEDED(hr) ? 0 : throw hr;

	XMLDoc->loadXML(result.bstrVal, &bOK);

	if ( XMLDoc->selectSingleNode(_bstr_t("/data/XMLResponse"), &node) == S_OK )
	{
		BSTR resultat = NULL;
		node->get_text(&resultat);

		std::string strRes="<Response>";
		ParseXml(&resultat,strRes);
		strRes.append("</Response>");
		xmlResponse = strRes.c_str();

		node->Release();
		node=NULL;
		if (resultat) SysFreeString(resultat);
	}

	if (XMLDoc)
		XMLDoc->Release();
	XMLDoc = NULL;

}



void ARM_EtkSoapClient::GetRefRate(CCString source, CCString ccy, CCString index, CCString tenor, CCString & xmlResponse, CCString & messageList)
{


    USES_CONVERSION;
	
	WCHAR *pwcMethodName = L"OneShotExecute";
	DISPID dispidFn = 0;
	HRESULT hr  = itsSoapClient->GetIDsOfNames(IID_NULL, &pwcMethodName, 1, 
	LOCALE_SYSTEM_DEFAULT, &dispidFn);
	CHECK_HR(hr);

	_bstr_t bAppName("ARM");
	_bstr_t bBase (itsBase);
	_bstr_t bLogin(itsUserName);
	_bstr_t bCommand ("s_base:DBQuery");

	CCString request = (CCString)"<Request><SummitSQL>select AsOfDate, Rate, Term from dmREFRATE where Source = '";
	request = request + source;
	request = request + (CCString)"' and Ccy = '";
	request = request + ccy;
	request = request + (CCString)"' and dmIndex = '";
	request = request + index;
	
	if (strcmp(tenor,"") != 0)
	{
		request = request + (CCString)"' and Term like '%";
		request = request + tenor;
	}

	request = request + (CCString)"'</SummitSQL></Request>";

	_bstr_t bRequest ((const char*)request);

	unsigned int uArgErr;
	VARIANT varg[5];
	varg[0].vt = VT_BSTR;
	varg[0].bstrVal = bRequest;
	varg[1].vt = VT_BSTR;
	varg[1].bstrVal = bCommand;
	varg[2].vt = VT_BSTR;
	varg[2].bstrVal = bLogin;
	varg[3].vt = VT_BSTR;
	varg[3].bstrVal = bBase;
	varg[4].vt = VT_BSTR;
	varg[4].bstrVal = bAppName;

	DISPPARAMS params;
	params.cArgs = 5;
	params.rgvarg = varg;
	params.cNamedArgs    = 0;
	params.rgdispidNamedArgs = NULL;

	_variant_t result;

	uArgErr = -1;

	EXCEPINFO excepInfo;
	memset(&excepInfo, 0, sizeof(excepInfo));
	hr  = itsSoapClient->Invoke(dispidFn, IID_NULL, LOCALE_SYSTEM_DEFAULT, 
	DISPATCH_METHOD, &params, &result, &excepInfo, &uArgErr);
	CHECK_HR(hr);

	VARIANT_BOOL bOK;
	MSXML2::IXMLDOMDocument *XMLDoc = NULL;
	MSXML2::IXMLDOMNode *node= NULL;
	
	hr = CoCreateInstance(__uuidof(MSXML2::DOMDocument30), NULL, 
			              CLSCTX_INPROC_SERVER, __uuidof(MSXML2::IXMLDOMDocument), (void**)&XMLDoc);

	SUCCEEDED(hr) ? 0 : throw hr;

	XMLDoc->loadXML(result.bstrVal, &bOK);

	if ( XMLDoc->selectSingleNode(_bstr_t("/data/XMLResponse"), &node) == S_OK )
	{
		BSTR resultat = NULL;
		node->get_text(&resultat);

//		_bstr_t ff(resultat,false);
//		char * ff1=(char *)ff;
		// JLA xmlResponse = W2A(resultat);
		std::string toto; 
		VariantTools::convert(resultat,toto); 
		xmlResponse=toto.c_str() ;


		node->Release();
		node=NULL;
		if (resultat) SysFreeString(resultat);
	}

	if (XMLDoc)
		XMLDoc->Release();
	XMLDoc = NULL;

}

