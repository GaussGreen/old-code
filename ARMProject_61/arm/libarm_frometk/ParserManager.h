/*
 *
 * Copyright (c) CDC IXIS CM October 2004 Paris
 *
 * $Log: $
 *
 */

/*! \file ParserManager.h
 *
 *  \brief
 *	\author  M. Campet
 *	\version 1.0
 *	\date October 2004
 */

#ifndef _PARSERMANAGER_H
#define _PARSERMANAGER_H


#include <string>
#include <map>
using namespace std;

#import <msxml3.dll> raw_interfaces_only
using namespace MSXML2;

class ARM_Security;
class ARM_Date;

namespace etoolkit
{

#define SUMMITSTRING "ARMNT_CALL/ASSETINFO/"
#define ETKSTRING "Response/"

//! pure virtual class that defines a Summit Instrument loader interface
class AbstractParser 
{
	public:
		virtual ~AbstractParser();

		//! retrieves a ARM Security from summit via etoolkit 
		virtual ARM_Object* BuildInstrument(MSXML2::IXMLDOMDocument *XMLDoc,
											  const ARM_Date& date,
											  const string filter,
											  string& bookName,
											  string& structureId,
											  string& custId,
											  string& dealId,
											  vector<string>& listAssetId,
											  int isEtkOrNot = 1) = 0;
		inline void SetSummitType(string type) { itsSummitType=type;};
		inline string GetSummitType(void) { return itsSummitType;};

		int GetResetFreqForCorridorOptim(void)		{ return itsResetFreqForCorridorOptim;}
		void SetResetFreqForCorridorOptim(int data)	{ itsResetFreqForCorridorOptim = data;}

	private:
		string itsSummitType;
		int itsResetFreqForCorridorOptim;
};

class ParserManager
{
	public:
		static void Init(int isEtkOrNot = 1);
		static void Release(void);
		static AbstractParser* GetParser(string type);
		static ARM_Object* BuildInstrument(const string type,
											 MSXML2::IXMLDOMDocument *XMLDoc,
											 const ARM_Date& date,
											 const string filter,
											 string& bookName,
											 string& structureId,
											 string& custId,
											 string& dealId,
											 vector<string>& listAssetId,
											 int aResetFreqForCorridorOptim = 365);

		static bool MustBeParsed(string& type, vector<string> listFilterPlus, vector<string> listFilterMinus);
		static vector<string> GetListFilter(const string filter, vector<string>& listFilterMinus);

		// TMP static string GetKey(MSXML2::IXMLDOMDocument*XMLDoc);

		static int IsEtkOrSummit(){ return isEtkOrSummit; };
	private:
		//! \brief registers Instrument loaders
		static bool RegisterParser(string& type, AbstractParser* parser);

		static map<string, AbstractParser*> itsParsers;
		static bool initDone;

		static int isEtkOrSummit;
};


class ConvertManager
{
	public:
		static void Init(void);
		static void Release(void);
		static string BuildSummitKey(MSXML2::IXMLDOMNode *theNode, int isEtkOrNot = 1);
		static string ConvertETKKey(string& SummitType);
		static string GetETKKey(MSXML2::IXMLDOMNode *theNode, int isEtkOrNot = 1);
		static string GetAssetKey(MSXML2::IXMLDOMNode *XMLNode);
		
	private:
		//! \brief registers Summit Type corresponding to ARM Type
		static bool RegisterKey(string& SummitType, string& Type);

		static map<string, string> itsConvertTable;
		static bool initDone;
};

} // namespace

#endif

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/
