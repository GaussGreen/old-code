
#ifndef _XMLTOOLS_H_
#define _XMLTOOLS_H_


#include <string>
#include <ICMKernel/util/icm_qmatrix.h>
#import <msxml3.dll> raw_interfaces_only


//	
class XMLTools
{
public:	
	//	----------------------------------------------------------
	//
	//		unitary conversion functions 
	//		all might throw 
	//
	static void convert(MSXML2::IXMLDOMNode*node,std::string&ret); 
	static void convert(MSXML2::IXMLDOMNode*node,double&ret); 
	static void convert(MSXML2::IXMLDOMNode*node,int& ret); 
	static void convert(MSXML2::IXMLDOMNode*node,ARM_Date&,const std::string& format); 
public:
	//	----------------------------------------------------------
	//
	//		smart pointer node manipulations
	//
	static MSXML2::IXMLDOMDocumentPtr CreateDOMDocument30() ;
	static MSXML2::IXMLDOMDocumentPtr LoadXML(const std::string& xmlContent) ;
	static MSXML2::IXMLDOMNodeListPtr selectNodes(MSXML2::IXMLDOMDocumentPtr&xmlDoc,const std::string& path) ;
	static MSXML2::IXMLDOMNodeListPtr selectNodes(MSXML2::IXMLDOMNodePtr&xmlDoc,const std::string& path_);
    static MSXML2::IXMLDOMNodePtr	  get_item(MSXML2::IXMLDOMNodeListPtr&xmlList,long indexNode );
	static MSXML2::IXMLDOMNodePtr	  selectSingleNode(MSXML2::IXMLDOMDocumentPtr&xmlDoc,const std::string& path) ;
	static MSXML2::IXMLDOMNodePtr	  selectSingleNode(MSXML2::IXMLDOMNodePtr&xmlNode,const std::string& path) ;
}; 

#endif	// _XMLTOOLS_H_
