/*
 *
 * Copyright (c) CDC IXIS CM October 2004 Paris
 *
 * $Log: $
 *
 */

/*! \file AssetParserManager.h
 *
 *  \brief
 *	\author  M. Campet
 *	\version 1.0
 *	\date October 2004
 */

#ifndef _ASSETPARSERMANAGER_H
#define _ASSETPARSERMANAGER_H


#include <string>
#include <map>
using namespace std;

#import <msxml3.dll> raw_interfaces_only
using namespace MSXML2;

//class ARM_Security;
//class ARM_Object;
class ARM_Date;
class ARM_ReferenceValue;
//class MSXML2::IXMLDOMDocument;





namespace etoolkit
{

//! pure virtual class that defines a Market Data loader interface
class AbstractAssetParser 
{
	public:
		virtual ~AbstractAssetParser();

		//! retrieves a ARM Security from summit via etoolkit 
		virtual ARM_Object* BuildAsset(MSXML2::IXMLDOMNode* XMLDoc,
									   const ARM_Date& date,
									   string& assetId,
									   int isEtkOrNot = 1) = 0;

		inline void SetSummitType(string type) { itsSummitType = type;};
		inline void SetTradeId(string tradeId) { itsTradeId = tradeId; };
		inline void SetTradeType(string tradeType) { itsTradeType = tradeType; };

		inline string GetSummitType(void) { return itsSummitType;};		
		inline string GetTradeId(void) { return itsTradeId; };
		inline string GetTradeType(void) { return itsTradeType; };

	
		//Used for TARN
		virtual ARM_Object* BuildTARNAsset( MSXML2::IXMLDOMNode* XMLDoc,
											const ARM_Date& date,
											string& assetId,
											int isEtkOrNot,
											bool oswCalibFlag){ return NULL; };

		inline double GetNbIterations(){return itsNbIterations;};
		inline void   SetNbIterations(double nbIter){itsNbIterations = nbIter;};

		int GetResetFreqForCorridorOptim(){return itsResetFreqForCorridorOptim;}
		void SetResetFreqForCorridorOptim(int data){itsResetFreqForCorridorOptim = data;}

	private:

		static string itsSummitType;		
		string itsTradeId;
		string itsTradeType;
		
		//Used for TARN
		double itsNbIterations;

		// use for Mercure optim
		int itsResetFreqForCorridorOptim;
};



class AssetParserManager
{
	public:

		static void Init(int isEtkOrNot = 1);
		static void Release(void);
		static AbstractAssetParser* GetAssetParser(string type);
		static ARM_Object* BuildAsset(const string type,
									 MSXML2::IXMLDOMDocument* XMLDoc,
									 string& assetId,
									 const ARM_Date& date);
		
		static void CorrectNotional(ARM_ReferenceValue* initNotional,
								    ARM_Security* object);

		static int IsEtkOrSummit(){ return isEtkOrSummit; };

	private:

		//! \brief registers Market Data loaders
		static bool RegisterParser(string& type, AbstractAssetParser* parser);

		static map<string, AbstractAssetParser*> itsAssetParsers;
		static bool initDone;
		static int isEtkOrSummit;
};

} // namespace



#endif
/*---------------------------------------------------------------------------*/
/*---- End of file ----*/
