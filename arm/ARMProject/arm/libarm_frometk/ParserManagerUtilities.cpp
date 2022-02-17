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

#include <atlbase.h>
#include <glob\dates.h>
#include <CCDate.h>
#include "PaserManagerUtilities.h"
#include "ICMKernel\util\icm_macro.h"

ARM_Date GetDateFromXMLNode (MSXML2::IXMLDOMNode* node, string nodeName)
{
	MSXML2::IXMLDOMNode* tmpNode = NULL;
	BSTR resultat = NULL;
	ARM_Date resDate;

	try
	{
		node->selectSingleNode((_bstr_t)nodeName.c_str(), &tmpNode);
		if (tmpNode != NULL)
		{
			tmpNode->get_text(&resultat);
			_bstr_t ff(resultat,false);
			char * ff1=(char *)ff;

			resDate = ARM_Date(ff1,"YYYYMMDD");

			if (resultat) SysFreeString(resultat);
			tmpNode->Release();
			tmpNode = NULL;
		}
		else
		{
			string msg="Pb in getting Date from node " + nodeName;
	
			throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
					 (char*) msg.c_str());
		}
	}
	catch(...)
	{
		string msg="Pb in getting Date from node " + nodeName;

		throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
				 (char*) msg.c_str());
	}

	return resDate;
}



double GetDoubleFromXMLNode (MSXML2::IXMLDOMNode* node, string nodeName)
{
	MSXML2::IXMLDOMNode* tmpNode = NULL;
	BSTR resultat = NULL;
	double resDouble = 0.0;

	try
	{
		node->selectSingleNode((_bstr_t)nodeName.c_str(), &tmpNode);
		if (tmpNode != NULL)
		{
			tmpNode->get_text(&resultat);
			_bstr_t ff(resultat,false);
			char * ff1=(char *)ff;

			resDouble = atof((const char*) ff1);

			if (resultat) SysFreeString(resultat);
			tmpNode->Release();
			tmpNode = NULL;
		}
	}
	catch(...)
	{
		string msg="Pb in getting Double from node " + nodeName;

		throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
				 (char*) msg.c_str());
	}

	return resDouble;
}

int GetIntFromXMLNode (MSXML2::IXMLDOMNode* node, string nodeName)
{
	MSXML2::IXMLDOMNode* tmpNode = NULL;
	BSTR resultat = NULL;
	int resInteger = 0;

	try
	{
		node->selectSingleNode((_bstr_t)nodeName.c_str(), &tmpNode);
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
		string msg="Pb in getting Int from node " + nodeName;

		throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
				 (char*) msg.c_str());
	}

	return resInteger;
}


string GetNodeName (MSXML2::IXMLDOMNode* node)
{
	BSTR resultat = NULL;
	string resString = "";

	try
	{
		node->get_nodeName(&resultat);
		_bstr_t ff(resultat,false);

		resString = string((char *)ff);
	}
	catch(...)
	{
		string msg="Pb in getting NodeName";

		throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
				 (char*) msg.c_str());
	}

	return resString;
}


// JLA. Same version with const args
//		+ catch standard exception
//		+ log message clickable
//
std::string GetStringFromXMLNode (MSXML2::IXMLDOMNode* node, const std::string & nodeName)
{
	MSXML2::IXMLDOMNode* tmpNode = NULL;
	BSTR resultat = NULL;
	string resString = "";

	try
	{
		node->selectSingleNode((_bstr_t)nodeName.c_str(), &tmpNode);
		if (tmpNode != NULL)
		{
			tmpNode->get_text(&resultat);
			_bstr_t ff(resultat,false);

			resString = string((char *)ff);

			if (resultat) SysFreeString(resultat);
			tmpNode->Release();
			tmpNode = NULL;
		}
	}
	catch (std::exception&e) 
	{
		ICMTHROW(ERR_INVALID_ARGUMENT,"GetStringFromXMLNode "<<nodeName<<" : " << e.what()) 
	}
	catch(...)
	{
		ICMTHROW(ERR_INVALID_ARGUMENT,"GetStringFromXMLNode "<<nodeName<<" : unknown exception") 
	}

	return resString;
}


// The same function will not throw on failure. 
bool GetStringFromXMLNode(MSXML2::IXMLDOMNode* node, const std::string & nodeName,std::string& result)
{
	try 
	{
		result=GetStringFromXMLNode(node,nodeName); 
	}
	catch(...)
	{
		return false; 
	}
	return true; 
}

