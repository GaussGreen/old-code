/*
 *
 * Copyright (c) CDC IXIS CM October 2004 Paris
 *
 * $Log: $
 *
 */

/*! \file ParserManagerUtilities.h
 *
 *  \brief
 *	\author  M. Campet
 *	\version 1.0
 *	\date October 2004
 */

#ifndef _PARSERMANAGER_UTIL_H
#define _PARSERMANAGER_UTIL_H


#include <string>
#include <ICMKernel/util/icm_qmatrix.h>
using namespace std;


#import <msxml3.dll> raw_interfaces_only
using namespace MSXML2;


class ARM_Date;

ARM_Date GetDateFromXMLNode (MSXML2::IXMLDOMNode* node, string nodeName);
double GetDoubleFromXMLNode (MSXML2::IXMLDOMNode* node, string nodeName);
int GetIntFromXMLNode (MSXML2::IXMLDOMNode* node, string nodeName);
string GetNodeName (MSXML2::IXMLDOMNode* node);

std::string GetStringFromXMLNode (MSXML2::IXMLDOMNode* node, const std::string & nodeName) ;
bool GetStringFromXMLNode(MSXML2::IXMLDOMNode* node, const std::string & nodeName,std::string& result) ; 



#endif

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/
