#include <ARM\libarm_local\firstToBeIncluded.h>

#include <atlbase.h>
#include <glob\dates.h>
//#include <CCDate.h>
#include "XMLTools.h"
#include "VariantTools.h"
using namespace MSXML2 ;

//	------------------------------------------------------------------------------------------------
void 
XMLTools::convert(MSXML2::IXMLDOMNode*node,std::string&ret)
{
	if (!node) ICMTHROW(ERR_INVALID_ARGUMENT,"XMLTools::convert: no node"); 
	
	CComBSTR tmp; 
	node->get_text(&tmp);
	ret=_bstr_t(tmp); 
}
//	------------------------------------------------------------------------------------------------
void 
XMLTools::convert(MSXML2::IXMLDOMNode*node,double&ret)
{
	std::string tmp; 
	convert(node,tmp); 
	ret=atof(tmp.c_str()); 
} 
//	------------------------------------------------------------------------------------------------
void 
XMLTools::convert(MSXML2::IXMLDOMNode*node,int& ret)
{
	std::string tmp; 
	convert(node,tmp); 
	ret=atoi(tmp.c_str()); 
}
//	------------------------------------------------------------------------------------------------
void 
XMLTools::convert(MSXML2::IXMLDOMNode*node,ARM_Date&ret,const std::string& format)
{
	std::string tmp;
	convert(node,tmp); 
	ret=ARM_Date(tmp,format);
}


//	----------------------------------------------------------
//
//		smart pointer node manipulations
//
//	------------------------------------------------------------------------------------------------
// static 
MSXML2::IXMLDOMDocumentPtr 
XMLTools::CreateDOMDocument30() 
{
	HRESULT hr ; 
	hr = CoInitialize(NULL); 
	//if (hr!=S_OK) 
	//	ICMTHROW(ERR_INVALID_ARGUMENT,"XMLTools::CreateMSXML2::DOMDocument30: failed"); 
	MSXML2::IXMLDOMDocumentPtr ptr; 
	ptr.CreateInstance(__uuidof(MSXML2::DOMDocument30)); 
	return ptr; 
}
//	------------------------------------------------------------------------------------------------
// static 
MSXML2::IXMLDOMNodeListPtr 
XMLTools::selectNodes(MSXML2::IXMLDOMDocumentPtr&xmlDoc,const std::string& path_) 
{
	HRESULT hr ;
	_bstr_t path ; 
	VariantTools::convert(path_,path); 
	MSXML2::IXMLDOMNodeList * item =0 ; 
	hr = xmlDoc->selectNodes(path,&item); 
	if (hr!=S_OK) 
		ICMTHROW(ERR_INVALID_ARGUMENT,"XMLTools::selectNodes: failed: path="<<path_); 
	MSXML2::IXMLDOMNodeListPtr ret ;
	ret.Attach(item) ; // AddRef is not called
	return ret; 
}

//	------------------------------------------------------------------------------------------------
// static 
MSXML2::IXMLDOMNodeListPtr 
XMLTools::selectNodes(MSXML2::IXMLDOMNodePtr&xmlDoc,const std::string& path_) 
{
	HRESULT hr ;
	_bstr_t path ; 
	VariantTools::convert(path_,path); 
	MSXML2::IXMLDOMNodeList * item =0 ; 
	hr = ((MSXML2::IXMLDOMNode *)xmlDoc)->selectNodes(path,&item); 
	if (hr!=S_OK) 
		ICMTHROW(ERR_INVALID_ARGUMENT,"XMLTools::selectNodes: failed: path="<<path_); 
	MSXML2::IXMLDOMNodeListPtr ret ;
	ret.Attach(item) ; // AddRef is not called
	return ret; 
}
//	------------------------------------------------------------------------------------------------
// static 
MSXML2::IXMLDOMNodePtr	  
XMLTools::get_item(MSXML2::IXMLDOMNodeListPtr&xmlList,long indexNode ) 
{

	HRESULT hr ;
	MSXML2::IXMLDOMNode * item =0 ; 
	hr = xmlList->get_item(indexNode,&item); 
	if (hr!=S_OK) 
		ICMTHROW(ERR_INVALID_ARGUMENT,"XMLTools::get_item: failed: node#"<<indexNode); 
	MSXML2::IXMLDOMNodePtr ret ;
	ret.Attach(item) ; // AddRef is not called
	return ret; 
}
//	------------------------------------------------------------------------------------------------
// static 
MSXML2::IXMLDOMNodePtr	  
XMLTools::selectSingleNode(MSXML2::IXMLDOMDocumentPtr&xmlDoc,const std::string& path_) 
{
	HRESULT hr ;
	MSXML2::IXMLDOMNode * item =0 ; 
	_bstr_t path ; 
	VariantTools::convert(path_,path); 
	hr = xmlDoc->selectSingleNode(path,&item); 
	if (hr!=S_OK) 
		ICMTHROW(ERR_INVALID_ARGUMENT,"XMLTools::selectSingleNode: failed: path="<<path_); 
	MSXML2::IXMLDOMNodePtr ret ;
	ret.Attach(item) ; // AddRef is not called
	return ret; 
}
//	------------------------------------------------------------------------------------------------
// static 
MSXML2::IXMLDOMNodePtr	  
XMLTools::selectSingleNode(MSXML2::IXMLDOMNodePtr&xmlDoc,const std::string& path_) 
{
	HRESULT hr ;
	MSXML2::IXMLDOMNode * item =0 ; 
	_bstr_t path ; 
	VariantTools::convert(path_,path); 
	hr = xmlDoc->selectSingleNode(path,&item); 
	if (hr!=S_OK) 
		ICMTHROW(ERR_INVALID_ARGUMENT,"XMLTools::selectSingleNode: failed: path="<<path_); 
	MSXML2::IXMLDOMNodePtr ret ;
	ret.Attach(item) ; // AddRef is not called
	return ret; 
}
//	------------------------------------------------------------------------------------------------
// static 
MSXML2::IXMLDOMDocumentPtr 
XMLTools::LoadXML(const std::string& xmlContent)
{
	VARIANT_BOOL bOK;
	MSXML2::IXMLDOMDocumentPtr xmlDoc=XMLTools::CreateDOMDocument30(); 
								   
	_bstr_t b;
	VariantTools::convert(xmlContent,b); 
	if (xmlDoc->loadXML(b,&bOK)!=S_OK)
		ICMTHROW(ERR_INVALID_ARGUMENT,"Can't load XML from string"); 
	return xmlDoc; 
}