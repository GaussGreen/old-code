/*
 *
 * Copyright (c) CDC IXIS CM October 2004 Paris
 *
 * $Log: $
 *
 */

/*! \file ParserManager.cpp
 *
 *  \brief
 *	\author  M. Campet
 *	\version 1.0
 *	\date October 2004
 */

#include "firsttoinc.h"
#include <glob\expt.h>
#include <inst\security.h>
#include <inst\spreadoption.h>
#include <inst\swap.h>
#include <inst\portfolio.h>
#include "ParserManager.h"
#include "AssetParserManager.h"
#include "PaserManagerUtilities.h"
#include "gpcalculators\callablesnowballcalculator.h"
#include "gpcalculators\tarncalculator.h"

// GP_Infra
#include "gpinfra\gensecurity.h"
#include "gpinfra\cstmanager.h"
#include "gpinfra\pricingmodel.h"
#include "gpinfra\pricingadviser.h"

// GP_Base
#include "gpbase\globalportfolio.h"

// GP_Calib
#include "gpcalib\calibmethod.h"


#if defined( va_start)
	#undef va_start
#endif

#if defined( va_end)
	#undef va_end
#endif

#include "ARM_local_parsexml.h"
#include "arm_local_parsexml_util.h"

using namespace etoolkit;
using namespace std;

#define SUMMITSTRING2 "ARMNT_CALL/OPSPEC"

using ARM::ARM_TARNCalculator;
using ARM::ARM_GlobalPortfolio;

bool ParserManager::initDone = false;
int ParserManager::isEtkOrSummit = 1;
bool ConvertManager::initDone = false;
map<string, AbstractParser*> ParserManager::itsParsers;
map<string, string> ConvertManager::itsConvertTable;



void ParserManager::Release(void) 
{
	if ( !initDone ) 
	{
		throw Exception (__LINE__, __FILE__, ERR_INVALID_OPERATION,
			 "ParserManager::Release, ParserManager::Init not done");
	}

	map<string, AbstractParser*>::iterator iter = itsParsers.begin();
	for (; iter != itsParsers.end(); iter++) 
	{
		delete iter->second;
	}
}



AbstractParser* ParserManager::GetParser(string type)
{
	if ( !initDone ) 
	{
		throw Exception (__LINE__, __FILE__, ERR_INVALID_OPERATION,
			 "ParserManager::GetParser, ParserManager::Init not done");
	}

	map<string, AbstractParser*>::iterator iter = itsParsers.find(type);	
	if ( iter != itsParsers.end() )
	{
		return iter->second;
	}
	else 
	{	
		CCString msg("ParserManager::GetParser, no loader for Instrument type ");
		msg += type.c_str();
		throw Exception (__LINE__, __FILE__, ERR_INVALID_OPERATION, msg);		
	}	
}



bool ParserManager::RegisterParser(string& key, AbstractParser* loader)
{
	pair<string, AbstractParser*> p(key, loader);
	pair<map<string, AbstractParser*>::iterator, bool> r = itsParsers.insert(p);

	return r.second;
}



bool ParserManager::MustBeParsed(string& type, vector<string> listFilterPlus, vector<string> listFilterMinus)
{
	string filter;

	for (int i=0; i < listFilterMinus.size(); i++)
	{
		filter = listFilterMinus[i];

		string::size_type pos = type.find(filter,0);
		if (pos != string::npos)
			return 0;
	}

	for (i=0; i < listFilterPlus.size(); i++)
	{
		filter = listFilterPlus[i];

		if ( filter == "ALL" )
			return 1;
		else
		{
			string::size_type pos = type.find(filter,0);
			if (pos != string::npos)
				return 1;
		}
	}

	return 0;
}



vector<string> ParserManager::GetListFilter(const string filter,vector<string>& listFilterMinus)
{
	vector<string> listFilterPlus;

	long sepPos=0;
	long sepPosPlus=0;
	long sepPosMinus=0;
	long sepPosPlusOld=0;
	long sepPosMinusOld=0;
	long sepPosOld=0;

	while ( (sepPosPlus != -1) || (sepPosMinus != -1) )
	{
		if (sepPosPlusOld == sepPosOld)
			sepPosPlus = filter.find_first_of("+",sepPosOld);

		if (sepPosMinusOld == sepPosOld)
			sepPosMinus = filter.find_first_of("-",sepPosOld);

		if (sepPosPlus == -1)
			sepPos = sepPosMinus;
		else if (sepPosMinus == -1)
			sepPos = sepPosPlus;
		else
			sepPos = MIN(sepPosPlus,sepPosMinus);

		string sFilter = string(filter,sepPosOld,sepPos-sepPosOld);

		if ( (sepPosOld == 0) || (sepPosPlusOld == sepPosOld) )
		{
			listFilterPlus.push_back(sFilter);
			sepPosPlusOld = sepPosPlus+1;
			if (sepPosOld == 0)
				sepPosMinusOld = sepPosMinus+1;
		}
		else
		{
			listFilterMinus.push_back(sFilter);
			sepPosMinusOld = sepPosMinus+1;
		}

		sepPosOld = sepPos+1;
	}

	return listFilterPlus;
}



ARM_Object* ParserManager::BuildInstrument(string type,
										   MSXML2::IXMLDOMDocument *XMLDoc,
										   const ARM_Date& date,
										   const string filter,
										   string& bookName,
										   string& structureId,
										   string& custId,
										   string& dealId,
										   VECTOR<string>& listAssetId,
										   int aResetFreqForCorridorOptim)
{
	ARM_Object* sec = NULL;

	string key = ConvertManager::GetETKKey(XMLDoc, IsEtkOrSummit());

	AbstractParser* p= GetParser(key);

    if ( p ) 
	{
		p->SetResetFreqForCorridorOptim(aResetFreqForCorridorOptim);
		sec = p->BuildInstrument(XMLDoc, date, filter, bookName, 
                                 structureId, custId, dealId, listAssetId, IsEtkOrSummit());
	}
	else 
	{
		CCString msg("ParserManager::BuildInstrument, no loader for Instrument type ");
		msg += type.c_str();
		throw Exception (__LINE__, __FILE__, ERR_INVALID_OPERATION, msg);		
	}

	return sec;
}



AbstractParser::~AbstractParser() 
{
}



//! \brief loader for Spread Option
class SpreadOptionLoader : public AbstractParser
{
public:
	SpreadOptionLoader()
	{
		SetSummitType("CAPTR");
	}

	~SpreadOptionLoader()
	{
	}

	ARM_Object* BuildInstrument(MSXML2::IXMLDOMDocument *XMLDoc,
								  const ARM_Date& date,
								  const string filter,
								  string& bookName,
								  string& structureId,
								  string& custId,
								  string& dealId,
								  vector<string>& listAssetId,
								  int isEtkOrNot = 1)
	{
		HRESULT hr;
		MSXML2::IXMLDOMNode* aNode = NULL;
		MSXML2::IXMLDOMNodeList* resultList = NULL;
		long nbNodes=0;
		ARM_Object* sec = NULL;

		try
		{
			// Informations contenues dans l'enveloppe du deal
			string chemin = ETKSTRING + GetSummitType() +"/Env/ENV";
			XMLDoc->selectSingleNode((_bstr_t)chemin.c_str(), &aNode);
			custId= GetStringFromXMLNode(aNode, "Cust");
			dealId = GetStringFromXMLNode(aNode, "DealId");
			bookName = GetStringFromXMLNode(aNode, "Book");
			if (aNode)
				aNode->Release();
			aNode = NULL;
			
			chemin = ETKSTRING + GetSummitType() +"/Assets/ASSET";

            if (XMLDoc->selectNodes((_bstr_t)chemin.c_str(), &resultList) == S_OK)
			{
				AssetParserManager::Init();

				hr=resultList->get_item(0, &aNode);
				if (hr==S_OK && aNode!=NULL)
				{
					string assetType = ConvertManager::GetAssetKey(aNode);

					AbstractAssetParser* p= AssetParserManager::GetAssetParser(assetType);
					if ( p ) 
					{
						string assetId;
						p->SetResetFreqForCorridorOptim(GetResetFreqForCorridorOptim());
						sec = p->BuildAsset(aNode, date, assetId);
						listAssetId.push_back(assetId);

						if (!sec)
						{
							throw Exception (__LINE__, __FILE__, ERR_INVALID_OPERATION,
							 "ParserManager::BuildInstrument,no loader for Instrument type %s", GetSummitType());
						}
					}
			
                    if (aNode)
						aNode->Release();
					aNode = NULL;
				}
			}

			if (resultList)
				resultList->Release();
			resultList = NULL;
		}

		catch (...) 
		{
			throw;
		} 

		return (ARM_Object*) sec;
	}
};



//! \brief loader for Exotic
class ExoticLoader : public AbstractParser
{
public:

	ExoticLoader()
	{
		SetSummitType("EXOTIC");
	}

	~ExoticLoader()
	{
	}

	ARM_Object* BuildInstrument(MSXML2::IXMLDOMDocument* XMLDoc,
								const ARM_Date& date,
								const string filter,
								string& bookName,
								string& structureId,
								string& custId,
								string& dealId,
								vector<string>& listAssetId,
								int isEtkOrNot = 1)
	{
		HRESULT hr;
		
        MSXML2::IXMLDOMNode* aNode          = NULL;
		MSXML2::IXMLDOMNodeList* resultList = NULL;
		
        long nbNodes = 0;

		ARM_StdPortfolio* exotic    = NULL;

		vector<ARM_Security*> assets;
		
        vector<double> weights;
		vector<double> prices;
		vector<double> vegas;
		
		// Just to get a calculator object.
		ARM_Object* calculator = NULL;
		bool isCalculator = false;

		

		try
		{
			vector<string> listFilterMinus;
			vector<string> listFilterPlus = ParserManager::GetListFilter(filter,listFilterMinus);

			string chemin;
			string status;

			if (isEtkOrNot)
			{
				// Informations contenues dans l'enveloppe du deal
				chemin = ETKSTRING + GetSummitType() +"/Env/ENV";
				XMLDoc->selectSingleNode((_bstr_t)chemin.c_str(), &aNode);

				custId= GetStringFromXMLNode(aNode, "Cust"); 
				dealId = GetStringFromXMLNode(aNode, "DealId");
				bookName = GetStringFromXMLNode(aNode, "Book");

				if (aNode)
					aNode->Release();
				aNode = NULL;

				chemin = ETKSTRING + GetSummitType() +"/Back/BACK";
				XMLDoc->selectSingleNode((_bstr_t)chemin.c_str(), &aNode);

				status = GetStringFromXMLNode(aNode, "TermAssignStatus");
				if (aNode)
					aNode->Release();
				aNode = NULL;
			}
			else
			{
				custId = "ZZZ"; 
				dealId = "ZZZ"; 
				bookName = "ZZZ"; 
				status="";
			}
           

			if ( strcmp(status.c_str(), "TERM") == 0 )
			{
			   exotic = new ARM_StdPortfolio(assets, weights, prices, (const vector<double>)NULL);
                
               // Free Cloned Assets

               int i;
               int sz = assets.size();

               for (i = 0; i < sz; i++)
               {
                   delete assets[i];

                   assets[i] = NULL;
               }

			   return((ARM_Object *) exotic);
			}
			
			//Chemin must depend on isEtkOrNot (serviceARM / Mercure)
			if (isEtkOrNot)
				chemin = ETKSTRING + GetSummitType() +"/Assets/ASSET";
			else
				chemin = SUMMITSTRING + string((const char *) "ASSET");
			
			
            if (XMLDoc->selectNodes((_bstr_t)chemin.c_str(), &resultList) == S_OK)
			{
				AssetParserManager::Init();

				resultList->get_length(&nbNodes);
				for (long indexAsset = 0 ; indexAsset < nbNodes ; indexAsset++)
				{
					hr=resultList->get_item(indexAsset, &aNode);
					if (hr == S_OK && aNode != NULL)
					{	
						string assetType = ConvertManager::GetAssetKey(aNode);
						if (ParserManager::MustBeParsed(assetType,listFilterPlus,listFilterMinus))
						{
							AbstractAssetParser* p= AssetParserManager::GetAssetParser(assetType);

							//We put a "rustine" for an exotic containing just one TARN Calculator and other non calculator products.
							if (strcmp(assetType.c_str(),"IRG.RFTARN.FORM") == 0)
							{
								isCalculator  = true;
								double nbIter = 100000;
								p->SetNbIterations(nbIter);
								string assetId;

								// get beta calib flags from CONFIG_VECTOR if SUMMIT
								bool configVectFlag = false;
								if (isEtkOrNot != 1)
								{
									long nbNodes = 0;
									MSXML2::IXMLDOMNodeList* resultList = NULL;
									MSXML2::IXMLDOMNode* item = NULL;
									if (XMLDoc->selectNodes(_bstr_t((const char *)"ARMNT_CALL/CONFIG_VECTOR"), &resultList) == S_OK)
									{
										resultList->get_length(&nbNodes);
										if (nbNodes == 0)
										{		
											CCString msg((CCString)"Invalid XML string for getting Config Vector");
											throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT, (char*) msg);
										}

										HRESULT hr=resultList->get_item(0, &item);
										configVectFlag = GetConfigVectorFlag(item, 0);

										if (item) 
											item->Release();
									}

									if (resultList)
										resultList->Release();
									resultList = NULL;
								}
								calculator = p->BuildTARNAsset(aNode, date, assetId, isEtkOrNot, configVectFlag);
							}
							else if (strcmp(assetType.c_str(),"IRG.ALMCAPTION.FORM") == 0)
							{
								//////////////////////////////////////////////////////////////////////////////
								// WARNING : Caption calculator doesn't inherit from ARM_Security !!
								// We can not add a caption leg in an ARM_Portfolio
								// We add a test:
								// - if we specify an "ALMCAPTION" filter, we will return only the caption leg
								// - otherwise, return a portfolio WITHOUT the caption leg
								//////////////////////////////////////////////////////////////////////////////
								for (int i = 0; i < listFilterPlus.size(); i++)
								{
									if (strcmp(listFilterPlus[i].c_str(), "ALMCAPTION") == 0)
									{
										isCalculator  = true;
										string assetId;
										calculator = p->BuildAsset(aNode, date, assetId);
									}
								}
							}
							else
							{
								if ( p ) 
								{
									ARM_Object* sec = NULL;
									string assetId;
			
									p->SetResetFreqForCorridorOptim(GetResetFreqForCorridorOptim());
                                    sec = p->BuildAsset(aNode, date, assetId);

									if (sec)
									{
										assets.push_back(dynamic_cast<ARM_Security*>(sec));
										weights.push_back(1.0);
										prices.push_back(1.0);
										vegas.push_back(0.0);
										listAssetId.push_back(assetId);
									}
								}
							}
						}

						if (aNode)
						   aNode->Release();
						aNode = NULL;
					}
				}
			}

			if (resultList)
			   resultList->Release();
			resultList = NULL;

			if (!isCalculator)
            {	
               exotic = new ARM_StdPortfolio(assets, weights, prices, vegas);

               // Free Cloned Assets

               int i;
               int sz = assets.size();

               for (i = 0; i < sz; i++)
               {
                   delete assets[i];

                   assets[i] = NULL;
               }
            }
		}

		catch (Exception& x)
		{
			throw x;
		}

		catch (...) 
		{
			throw;
		}

		if (isCalculator)
		   return (ARM_Object *) calculator;
		else
		   return (ARM_Object *) exotic;
	}
	
};


//! \brief loader for Global Exotic (security & calculator)
class GlobalExoticLoader : public AbstractParser
{
public:

	GlobalExoticLoader()
	{
		SetSummitType("GLOBAL_EXOTIC");
	}

	~GlobalExoticLoader()
	{
	}

	ARM_Object* BuildInstrument(MSXML2::IXMLDOMDocument* XMLDoc,
								const ARM_Date& date,
								const string filter,
								string& bookName,
								string& structureId,
								string& custId,
								string& dealId,
								vector<string>& listAssetId,
								int isEtkOrNot = 1)
	{
		HRESULT hr;
		
        MSXML2::IXMLDOMNode* aNode          = NULL;
		MSXML2::IXMLDOMNodeList* resultList = NULL;
		
        long nbNodes = 0;

		ARM_GlobalPortfolio* exotic    = NULL;

		vector<ARM_Object*> assets;
		
        std::vector<double> weights;
		std::vector<double> prices;

		try
		{
			vector<string> listFilterMinus;
			vector<string> listFilterPlus = ParserManager::GetListFilter(filter,listFilterMinus);

			// Informations contenues dans l'enveloppe du deal
			string chemin = ETKSTRING + GetSummitType() +"/Env/ENV";
			XMLDoc->selectSingleNode((_bstr_t)chemin.c_str(), &aNode);
			structureId = GetStringFromXMLNode(aNode, "StructureId");
			custId   = GetStringFromXMLNode(aNode, "Cust");
			dealId   = GetStringFromXMLNode(aNode, "DealId");
			bookName = GetStringFromXMLNode(aNode, "Book");
		
            if (aNode)
			   aNode->Release();
			aNode = NULL;

			// Informations contenues dans le Back du deal
			chemin = ETKSTRING + GetSummitType() +"/Back/BACK";
			XMLDoc->selectSingleNode((_bstr_t)chemin.c_str(), &aNode);
			string status = GetStringFromXMLNode(aNode, "TermAssignStatus");

            if (aNode)
				aNode->Release();
			aNode = NULL;

			if ( strcmp(status.c_str(), "TERM") == 0 )
			{
			   exotic = new ARM_GlobalPortfolio(assets, weights, prices);
                
               // Free Cloned Assets

               int i;
               int sz = assets.size();

               for (i = 0; i < sz; i++)
               {
                   delete assets[i];

                   assets[i] = NULL;
               }

			   return((ARM_Object*) exotic);
			}
			
			chemin = ETKSTRING + GetSummitType() +"/Assets/ASSET";
			
            if (XMLDoc->selectNodes((_bstr_t)chemin.c_str(), &resultList) == S_OK)
			{
				AssetParserManager::Init();

				resultList->get_length(&nbNodes);
				for (long indexAsset = 0 ; indexAsset < nbNodes ; indexAsset++)
				{
					hr=resultList->get_item(indexAsset, &aNode);
					if (hr == S_OK && aNode != NULL)
					{	
						string assetType = ConvertManager::GetAssetKey(aNode);
						if (ParserManager::MustBeParsed(assetType,listFilterPlus,listFilterMinus))
						{
							AbstractAssetParser* p= AssetParserManager::GetAssetParser(assetType);

							//Special case for TARN calculator
							if (strcmp(assetType.c_str(),"IRG.RFTARN.FORM") == 0)
							{
								ARM_Object* sec = NULL;
								double nbIter = 100000;
								p->SetNbIterations(nbIter);
								string assetId;

								// get beta calib flags from CONFIG_VECTOR if SUMMIT
								bool configVectFlag = false;
								if (isEtkOrNot != 1)
								{
									MSXML2::IXMLDOMNodeList* resultList = NULL;
									MSXML2::IXMLDOMNode* item = NULL;
									ARM_Object* sec = NULL;
									long nbNodes = 0;

									if (XMLDoc->selectNodes(_bstr_t((const char *)"ARMNT_CALL/CONFIG_VECTOR"), &resultList) == S_OK)
									{
										resultList->get_length(&nbNodes);
										if (nbNodes == 0)
										{		
											CCString msg((CCString)"Invalid XML string for getting Config Vector");
											throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT, (char*) msg);
										}

										HRESULT hr=resultList->get_item(0, &item);
										configVectFlag = GetConfigVectorFlag(item, 0);

										if (item) 
											item->Release();
									}

									if (resultList)
										resultList->Release();
									resultList = NULL;
								}

								sec = p->BuildTARNAsset(aNode, date, assetId, isEtkOrNot, configVectFlag);

								if (sec)
								{
									assets.push_back(dynamic_cast<ARM_Object*>(sec));
									weights.push_back(1.0);
									prices.push_back(1.0);
									listAssetId.push_back(assetId);
								}
							}
							else
							{
								if ( p ) 
								{
									ARM_Object* sec = NULL;
									string assetId;
			
									p->SetResetFreqForCorridorOptim(GetResetFreqForCorridorOptim());
                                    sec = p->BuildAsset(aNode, date, assetId);

									if (sec)
									{
										assets.push_back(dynamic_cast<ARM_Object*>(sec));
										weights.push_back(1.0);
										prices.push_back(1.0);
										listAssetId.push_back(assetId);
									}
								}
							}
						}

						if (aNode)
						   aNode->Release();
						aNode = NULL;
					}
				}
			}

			if (resultList)
			   resultList->Release();
			resultList = NULL;

	
           exotic = new ARM_GlobalPortfolio(assets, weights, prices);

           // Free Cloned Assets

           int i;
           int sz = assets.size();

           for (i = 0; i < sz; i++)
           {
               delete assets[i];

               assets[i] = NULL;
           }

		}

		catch (Exception& x)
		{
			throw x;
		}

		catch (...) 
		{
			throw;
		}
		
		return (ARM_Object*) exotic;
	}	
};

class MaturityCapLoader : public AbstractParser
{
public:
	MaturityCapLoader()
	{
		SetSummitType("CAPTR");
	}

	~MaturityCapLoader()
	{
	}

	ARM_Object* BuildInstrument(MSXML2::IXMLDOMDocument *XMLDoc,
								  const ARM_Date& date,
								  const string filter,
								  string& bookName,
								  string& structureId,
								  string& custId,
								  string& dealId,
								  vector<string>& listAssetId,
								  int isEtkOrNot = 1)
	{
		HRESULT hr;
		MSXML2::IXMLDOMNode* aNode = NULL;
		MSXML2::IXMLDOMNodeList* resultList = NULL;
		long nbNodes=0;
		ARM_Object* sec = NULL;

		try
		{
			// Informations contenues dans l'enveloppe du deal
			string chemin;
			
			if (isEtkOrNot)
			{
				chemin = ETKSTRING + GetSummitType() +"/Env/ENV";
				
				XMLDoc->selectSingleNode((_bstr_t)chemin.c_str(), &aNode);

				custId= GetStringFromXMLNode(aNode, "Cust"); 
				dealId = GetStringFromXMLNode(aNode, "DealId");
				bookName = GetStringFromXMLNode(aNode, "Book");

				if (aNode)
					aNode->Release();
				aNode = NULL;
			}
			else
			{
				custId = "ZZZ"; 
				dealId = "ZZZ"; 
				bookName = "ZZZ"; 
			}

			string tradeId;
			string tradeType;

			if (isEtkOrNot)
			{
				chemin = ETKSTRING + GetSummitType();
				XMLDoc->selectSingleNode((_bstr_t)chemin.c_str(), &aNode);

				tradeId = GetStringFromXMLNode(aNode, "TradeId"); 
				tradeType = GetStringFromXMLNode(aNode, "TradeType"); 
			}
			else
			{
				chemin = SUMMITSTRING + string((const char *) "ASSET");
				XMLDoc->selectSingleNode((_bstr_t)chemin.c_str(), &aNode);

				tradeId = GetStringFromXMLNode(aNode, "TradeId"); 
				tradeType = GetStringFromXMLNode(aNode, "Type");
			}

			if (aNode)
				aNode->Release();
			aNode = NULL;

			if (isEtkOrNot)
				chemin = ETKSTRING + GetSummitType() +"/Assets/ASSET";
			else
				chemin = SUMMITSTRING + string((const char *) "ASSET");
		
			if ( XMLDoc->selectNodes((_bstr_t)chemin.c_str(), &resultList) == S_OK )
			{
				AssetParserManager::Init(isEtkOrNot);

				hr = resultList->get_item(0, &aNode);
				
				if (( hr == S_OK ) && ( aNode != NULL))
				{
					string assetType = ConvertManager::GetAssetKey(aNode);

					AbstractAssetParser* p= AssetParserManager::GetAssetParser(assetType);
					
					if ( p ) 
					{
						p->SetTradeId(tradeId);
						p->SetTradeType(tradeType);
				
						string assetId;

						sec = p->BuildAsset(aNode, date,assetId,isEtkOrNot);
				
						listAssetId.push_back(assetId);

						if (!sec)
						{
							CCString msg("ParserManager::BuildInstrument,no loader for Instrument type ");
							msg += GetSummitType().c_str();
							throw Exception (__LINE__, __FILE__, ERR_INVALID_OPERATION, msg);
						}
					}
					if (aNode)
						aNode->Release();
					aNode = NULL;
				}
			}

			if (resultList)
			   resultList->Release();
			resultList = NULL;
		}

		catch (...) 
		{
			throw;
		} 

		return((ARM_Object *) sec);
	}
	
};



class MemorySOLoader : public AbstractParser
{
public:

    MemorySOLoader()
	{
		SetSummitType("CAPTR");
	}

   ~MemorySOLoader()
	{
	}

	ARM_Object* BuildInstrument(MSXML2::IXMLDOMDocument *XMLDoc,
								const ARM_Date& date,
								const string filter,
								string& bookName,
								string& structureId,
								string& custId,
								string& dealId,
								vector<string>& listAssetId,
								 int isEtkOrNot = 1)
	{
		HRESULT hr;
		MSXML2::IXMLDOMNode* aNode = NULL;
		MSXML2::IXMLDOMNodeList* resultList = NULL;
		
        long nbNodes = 0;

		ARM_Object* sec = NULL;

		try
		{
			// Informations contenues dans l'enveloppe du deal
			string chemin;
			
			if (isEtkOrNot)
			{
				chemin = ETKSTRING + GetSummitType() +"/Env/ENV";
				
				XMLDoc->selectSingleNode((_bstr_t)chemin.c_str(), &aNode);

				custId= GetStringFromXMLNode(aNode, "Cust"); 
				dealId = GetStringFromXMLNode(aNode, "DealId");
				bookName = GetStringFromXMLNode(aNode, "Book");

				if (aNode)
					aNode->Release();
				aNode = NULL;
			}
			else
			{
				custId = "ZZZ"; 
				dealId = "ZZZ"; 
				bookName = "ZZZ"; 
			}

			string tradeId;
			string tradeType;

			if (isEtkOrNot)
			{
				chemin = ETKSTRING + GetSummitType();
				XMLDoc->selectSingleNode((_bstr_t)chemin.c_str(), &aNode);

				tradeId = GetStringFromXMLNode(aNode, "TradeId"); 
				tradeType = GetStringFromXMLNode(aNode, "TradeType"); 
			}
			else
			{
				chemin = SUMMITSTRING + string((const char *) "ASSET");
				XMLDoc->selectSingleNode((_bstr_t)chemin.c_str(), &aNode);

				tradeId = GetStringFromXMLNode(aNode, "TradeId"); 
				tradeType = GetStringFromXMLNode(aNode, "Type");
			}

			if (aNode)
				aNode->Release();
			aNode = NULL;

			if (isEtkOrNot)
				chemin = ETKSTRING + GetSummitType() +"/Assets/ASSET";
			else
				chemin = SUMMITSTRING + string((const char *) "ASSET");
		
			if ( XMLDoc->selectNodes((_bstr_t)chemin.c_str(), &resultList) == S_OK )
			{
				AssetParserManager::Init(isEtkOrNot);

				hr = resultList->get_item(0, &aNode);
				
				if (( hr == S_OK ) && ( aNode != NULL))
				{
					string assetType = ConvertManager::GetAssetKey(aNode);

					AbstractAssetParser* p= AssetParserManager::GetAssetParser(assetType);
					
					if ( p ) 
					{
						p->SetTradeId(tradeId);
						p->SetTradeType(tradeType);
				
						string assetId;

						sec = p->BuildAsset(aNode, date,assetId,isEtkOrNot);
				
						listAssetId.push_back(assetId);

						if (!sec)
						{
							CCString msg("ParserManager::BuildInstrument,no loader for Instrument type ");
							msg += GetSummitType().c_str();
							throw Exception (__LINE__, __FILE__, ERR_INVALID_OPERATION, msg);
						}
					}
					
                    if (aNode)
					   aNode->Release();
					aNode = NULL;
				}
			}

			if (resultList)
			   resultList->Release();
			resultList = NULL;
		}

		catch (...) 
		{
			throw;
		}

		return((ARM_Object *) sec);
	}
	
};




class MemoryCapLoader : public AbstractParser
{
public:
	MemoryCapLoader()
	{
		SetSummitType("CAPTR");
	}

	~MemoryCapLoader()
	{
	}

	ARM_Object* BuildInstrument(MSXML2::IXMLDOMDocument *XMLDoc,
								const ARM_Date& date,
								const string filter,
								string& bookName,
								string& structureId,
								string& custId,
								string& dealId,
								vector<string>& listAssetId,
								int isEtkOrNot = 1)
	{
		HRESULT hr;
		MSXML2::IXMLDOMNode* aNode          = NULL;
		MSXML2::IXMLDOMNodeList* resultList = NULL;

		long nbNodes = 0;
		ARM_Object* sec = NULL;

		try
		{
			// Informations contenues dans l'enveloppe du deal
			string chemin;
			
			if (isEtkOrNot)
			{
				chemin = ETKSTRING + GetSummitType() +"/Env/ENV";
				
				XMLDoc->selectSingleNode((_bstr_t)chemin.c_str(), &aNode);

				custId= GetStringFromXMLNode(aNode, "Cust"); 
				dealId = GetStringFromXMLNode(aNode, "DealId");
				bookName = GetStringFromXMLNode(aNode, "Book");

				if (aNode)
					aNode->Release();
				aNode = NULL;
			}
			else
			{
				custId = "ZZZ"; 
				dealId = "ZZZ"; 
				bookName = "ZZZ"; 
			}

			string tradeId;
			string tradeType;

			if (isEtkOrNot)
			{
				chemin = ETKSTRING + GetSummitType();
				XMLDoc->selectSingleNode((_bstr_t)chemin.c_str(), &aNode);

				tradeId = GetStringFromXMLNode(aNode, "TradeId"); 
				tradeType = GetStringFromXMLNode(aNode, "TradeType"); 
			}
			else
			{
				chemin = SUMMITSTRING + string((const char *) "ASSET");
				XMLDoc->selectSingleNode((_bstr_t)chemin.c_str(), &aNode);

				tradeId = GetStringFromXMLNode(aNode, "TradeId"); 
				tradeType = GetStringFromXMLNode(aNode, "Type");
			}

			if (aNode)
				aNode->Release();
			aNode = NULL;

			if (isEtkOrNot)
				chemin = ETKSTRING + GetSummitType() +"/Assets/ASSET";
			else
				chemin = SUMMITSTRING + string((const char *) "ASSET");
		
			if ( XMLDoc->selectNodes((_bstr_t)chemin.c_str(), &resultList) == S_OK )
			{
				AssetParserManager::Init(isEtkOrNot);

				hr = resultList->get_item(0, &aNode);
				
				if (( hr == S_OK ) && ( aNode != NULL))
				{
					string assetType = ConvertManager::GetAssetKey(aNode);

					AbstractAssetParser* p= AssetParserManager::GetAssetParser(assetType);
					
					if ( p ) 
					{
						p->SetTradeId(tradeId);
						p->SetTradeType(tradeType);
				
						string assetId;

						sec = p->BuildAsset(aNode, date,assetId,isEtkOrNot);
				
						listAssetId.push_back(assetId);

						if (!sec)
						{
							CCString msg("ParserManager::BuildInstrument,no loader for Instrument type ");
							msg += GetSummitType().c_str();
							throw Exception (__LINE__, __FILE__, ERR_INVALID_OPERATION, msg);
						}
					}
					if (aNode)
						aNode->Release();
					aNode = NULL;
				}
			}
			if (resultList)
				resultList->Release();
			resultList = NULL;
		}

		catch (...) 
		{
			throw;
		}

		return((ARM_Object *) sec);
	}
	
};


//class to read the Xml string. CSB booked like a swaption.
/*class CallableSnowBallLoader : public AbstractParser
{
public:
	CallableSnowBallLoader()
	{
		SetSummitType("SWOPT"); //Vérify with Summit team
	}

	~CallableSnowBallLoader()
	{
	}

	ARM_Object* BuildInstrument(MSXML2::IXMLDOMDocument *XMLDoc,
								  const ARM_Date& date,
								  const string filter,
								  string& bookName,
								  string& structureId,
								  string& custId,
								  string& dealId,
								  vector<string>& listAssetId,
								  int isEtkOrNot = 1)
	{
		HRESULT hr;
		MSXML2::IXMLDOMNode* aNode = NULL;
		MSXML2::IXMLDOMNode* aNode2 = NULL;
		MSXML2::IXMLDOMNodeList* assetList = NULL;
		MSXML2::IXMLDOMNodeList* optionList = NULL;
		long nbAssets=0;
		ARM_Object* sec = NULL;
		BSTR resultat;
		ARM_CallableSnowBallCalculator* newCsb = NULL;

		try
		{
			// Informations contenues dans l'enveloppe du deal
			string chemin;
			
			if (isEtkOrNot)
			{
				chemin = ETKSTRING + GetSummitType() +"/Env/ENV";
				
				XMLDoc->selectSingleNode((_bstr_t)chemin.c_str(), &aNode);

				custId= GetStringFromXMLNode(aNode, "Cust"); 
				dealId = GetStringFromXMLNode(aNode, "DealId");
				bookName = GetStringFromXMLNode(aNode, "Book");

				if (aNode)
					aNode->Release();
				aNode = NULL;
			}
			else
			{
				custId = "ZZZ"; 
				dealId = "ZZZ"; 
				bookName = "ZZZ"; 
			}

			//Récupération de la partie ASSET'S
			if (isEtkOrNot)
				chemin = ETKSTRING + GetSummitType() +"/Assets/ASSET";
			else
				chemin = SUMMITSTRING + string((const char *) "ASSET");
		
			if ( XMLDoc->selectNodes((_bstr_t)chemin.c_str(), &assetList) == S_OK )
			{
				assetList->get_length(&nbAssets);
				if (nbAssets != 2)
				{
					CCString msg((CCString)"Option XML string is not valid \n" );
					throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT, (const char *) msg);
				}
				AssetParserManager::Init(isEtkOrNot);

				//Récupération des 2 ASSETS
				for (long indexAsset=0 ; indexAsset<nbAssets ; indexAsset++)
				{
					hr = assetList->get_item(indexAsset, &aNode);
					if (( hr == S_OK ) && ( aNode != NULL))
					{
						aNode->selectSingleNode(_bstr_t((const char *)"INTEREST_FixFloat"), &aNode2);
						if(aNode2 != NULL)
						{
							aNode2->get_text(&resultat);
							aNode2->Release();
							aNode2 = NULL;
							//Jambe variable de la swaption
							if ( strcmp((const char*)(_bstr_t)resultat,"FLO") == 0 )	
							{
								//Appeler le reader de la jambe variable swaption
							}
							//Jambe fixe de la swaption
							else //Jambe fixe
							{
								//Appeler le reader de la jambe fixe swaption
							}
							if (resultat) SysFreeString(resultat);

						}
					}
				}
			}
			if (assetList) 
				assetList->Release();
			assetList = NULL;


			//Récupération de la partie OPTION
			if (isEtkOrNot)
				chemin = ETKSTRING + GetSummitType() +"/Option/OPTION";
			else
				chemin = SUMMITSTRING + string((const char *) "OPTION");			
			if ( XMLDoc->selectNodes((_bstr_t)chemin.c_str(), &optionList) == S_OK )	
			{
				optionList->get_length(&nbAssets);
				if ( nbAssets != 1 )
				{
					CCString msg((CCString)"Option XML string is not valid \n" );
					throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,(char*) msg);
				}
				hr = optionList->get_item(0, &aNode);
				
				//Appeler le reader de la partie optionnelle.

			}

			//Appel du constructeur du callable snowball.
				
		}

		catch (...) 
		{
			throw;
		} 

		return((ARM_Object *) sec);
	}
	
};*/



// anciennement dérivé de MaturityCap.
//On séparer car le nombre d'itérations est différent.
class TarnLoader : public AbstractParser
{
public:
	TarnLoader()
	{
		SetSummitType("CAPTR");
	}

	~TarnLoader()
	{
	}

	ARM_Object* BuildInstrument(MSXML2::IXMLDOMDocument *XMLDoc,
								  const ARM_Date& date,
								  const string filter,
								  string& bookName,
								  string& structureId,
								  string& custId,
								  string& dealId,
								  vector<string>& listAssetId,
								  int isEtkOrNot = 1)
	{
		HRESULT hr;
		MSXML2::IXMLDOMNode* aNode = NULL;
		MSXML2::IXMLDOMNodeList* resultList = NULL;
		long nbNodes=0;
		ARM_Object* sec = NULL;

		try
		{
			//Informations contenues dans l'enveloppe du deal
			string chemin;
			
			if (isEtkOrNot)
			{
				chemin = ETKSTRING + GetSummitType() +"/Env/ENV";
				
				XMLDoc->selectSingleNode((_bstr_t)chemin.c_str(), &aNode);

				custId= GetStringFromXMLNode(aNode, "Cust"); 
				dealId = GetStringFromXMLNode(aNode, "DealId");
				bookName = GetStringFromXMLNode(aNode, "Book");

				if (aNode)
					aNode->Release();
				aNode = NULL;
			}
			else
			{
				custId		= "ZZZ"; 
				dealId		= "ZZZ"; 
				bookName	= "ZZZ"; 
			}

			string tradeId;
			string tradeType;
			double nbIter;

			if (isEtkOrNot)
			{
				chemin = ETKSTRING + GetSummitType();
				XMLDoc->selectSingleNode((_bstr_t)chemin.c_str(), &aNode);

				tradeId		= GetStringFromXMLNode(aNode, "TradeId"); 
				tradeType	= GetStringFromXMLNode(aNode, "TradeType"); 

				nbIter = 100000;
			}
			else
			{
				//NbIterations on node OPSPEC
				chemin = SUMMITSTRING2;
				XMLDoc->selectSingleNode((_bstr_t)chemin.c_str(), &aNode);
				nbIter = GetDoubleFromXMLNode(aNode, "NumIter");
				if (nbIter <= 0) //Field not existing
					nbIter = 100000; //Default value
			
				chemin = SUMMITSTRING + string((const char *) "ASSET");
				XMLDoc->selectSingleNode((_bstr_t)chemin.c_str(), &aNode);

				tradeId = GetStringFromXMLNode(aNode, "TradeId"); 
				tradeType = GetStringFromXMLNode(aNode, "Type");
			}

			if (aNode)
				aNode->Release();
			aNode = NULL;

			if (isEtkOrNot)
				chemin = ETKSTRING + GetSummitType() +"/Assets/ASSET";
			else
				chemin = SUMMITSTRING + string((const char *) "ASSET");
		
			if ( XMLDoc->selectNodes((_bstr_t)chemin.c_str(), &resultList) == S_OK )
			{
				AssetParserManager::Init(isEtkOrNot);

				hr = resultList->get_item(0, &aNode);
				
				if (( hr == S_OK ) && ( aNode != NULL))
				{
					string assetType = ConvertManager::GetAssetKey(aNode);

					AbstractAssetParser* p= AssetParserManager::GetAssetParser(assetType);
					
					if ( p ) 
					{
						p->SetTradeId(tradeId);
						p->SetTradeType(tradeType);
						p->SetNbIterations(nbIter);

						string assetId;

						bool configVectFlag = false;
						if (isEtkOrNot != 1)
						{
							long nbNodes = 0;
							MSXML2::IXMLDOMNodeList* resultList = NULL;
							MSXML2::IXMLDOMNode* item = NULL;
							if (XMLDoc->selectNodes(_bstr_t((const char *)"ARMNT_CALL/CONFIG_VECTOR"), &resultList) == S_OK)
							{
								resultList->get_length(&nbNodes);
								if (nbNodes == 0)
								{		
									CCString msg((CCString)"Invalid XML string for getting Config Vector");
									throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT, (char*) msg);
								}

								HRESULT hr=resultList->get_item(0, &item);
								configVectFlag = GetConfigVectorFlag(item, 0);

								if (item) 
									item->Release();
							}

							if (resultList)
								resultList->Release();
							resultList = NULL;
						}
						sec = p->BuildTARNAsset(aNode, date, assetId, isEtkOrNot, configVectFlag);
				
						listAssetId.push_back(assetId);

						if (!sec)
						{
							CCString msg("ParserManager::BuildInstrument,no loader for Instrument type ");
							msg += GetSummitType().c_str();
							throw Exception (__LINE__, __FILE__, ERR_INVALID_OPERATION, msg);
						}
					}
					if (aNode)
						aNode->Release();
					aNode = NULL;
				}
			}
			if (resultList)
				resultList->Release();
			resultList = NULL;
		}

		catch (...) 
		{
			throw;
		} 

		return((ARM_Object *) sec);
	}
};



class CaptionLoader : public AbstractParser
{

public:

    CaptionLoader()
	{
		SetSummitType("CAPTR");
	}

   ~CaptionLoader()
	{
	}

	ARM_Object* BuildInstrument(MSXML2::IXMLDOMDocument *XMLDoc,
								  const ARM_Date& date,
								  const string filter,
								  string& bookName,
								  string& structureId,
								  string& custId,
								  string& dealId,
								  vector<string>& listAssetId,
								  int isEtkOrNot = 1)
	{
		HRESULT hr;
		MSXML2::IXMLDOMNode* aNode = NULL;
		MSXML2::IXMLDOMNodeList* resultList = NULL;
		long nbNodes=0;
		ARM_Object* sec = NULL;

		try
		{
			//Informations contenues dans l'enveloppe du deal
			string chemin;
			
			if (isEtkOrNot)
			{
				chemin = ETKSTRING + GetSummitType() +"/Env/ENV";
				
				XMLDoc->selectSingleNode((_bstr_t)chemin.c_str(), &aNode);

				custId= GetStringFromXMLNode(aNode, "Cust"); 
				dealId = GetStringFromXMLNode(aNode, "DealId");
				bookName = GetStringFromXMLNode(aNode, "Book");

				if (aNode)
					aNode->Release();
				aNode = NULL;
			}
			else
			{
				custId		= "ZZZ"; 
				dealId		= "ZZZ"; 
				bookName	= "ZZZ"; 
			}

			string tradeId;
			string tradeType;
			double nbIter;

			if (isEtkOrNot)
			{
				chemin = ETKSTRING + GetSummitType();
				XMLDoc->selectSingleNode((_bstr_t)chemin.c_str(), &aNode);

				tradeId		= GetStringFromXMLNode(aNode, "TradeId"); 
				tradeType	= GetStringFromXMLNode(aNode, "TradeType"); 

				nbIter = 100000;
			}
			else
			{
				//NbIterations on node OPSPEC
				chemin = SUMMITSTRING2;
				XMLDoc->selectSingleNode((_bstr_t)chemin.c_str(), &aNode);
				nbIter = GetDoubleFromXMLNode(aNode, "NumIter");
				if (nbIter <= 0) //Field not existing
					nbIter = 100000; //Default value
			
				chemin = SUMMITSTRING + string((const char *) "ASSET");
				XMLDoc->selectSingleNode((_bstr_t)chemin.c_str(), &aNode);

				tradeId = GetStringFromXMLNode(aNode, "TradeId"); 
				tradeType = GetStringFromXMLNode(aNode, "Type");
			}

			if (aNode)
				aNode->Release();
			aNode = NULL;

			if (isEtkOrNot)
				chemin = ETKSTRING + GetSummitType() +"/Assets/ASSET";
			else
				chemin = SUMMITSTRING + string((const char *) "ASSET");
		
			if ( XMLDoc->selectNodes((_bstr_t)chemin.c_str(), &resultList) == S_OK )
			{
				AssetParserManager::Init(isEtkOrNot);

				hr = resultList->get_item(0, &aNode);
				
				if (( hr == S_OK ) && ( aNode != NULL))
				{
					string assetType = ConvertManager::GetAssetKey(aNode);

					AbstractAssetParser* p= AssetParserManager::GetAssetParser(assetType);
					
					if ( p ) 
					{
						p->SetTradeId(tradeId);
						p->SetTradeType(tradeType);
						p->SetNbIterations(nbIter);

						string assetId;

						sec = p->BuildAsset(aNode, date,assetId,isEtkOrNot);
				
						listAssetId.push_back(assetId);

						if (!sec)
						{
							CCString msg("ParserManager::BuildInstrument,no loader for Instrument type ");
							msg += GetSummitType().c_str();
							throw Exception (__LINE__, __FILE__, ERR_INVALID_OPERATION, msg);
						}
					}
					if (aNode)
						aNode->Release();
					aNode = NULL;
				}
			}
			if (resultList)
				resultList->Release();
			resultList = NULL;
		}

		catch (...) 
		{
			throw;
		} 

		return((ARM_Object *) sec);
	}
};



// Cap Loader
class CapLoader : public AbstractParser
{

public:

    CapLoader()
	{
		SetSummitType("CAPTR");
	}

   ~CapLoader()
	{
	}

	ARM_Object* BuildInstrument(MSXML2::IXMLDOMDocument *XMLDoc,
								const ARM_Date& date,
								const string filter,
								string& bookName,
								string& structureId,
								string& custId,
								string& dealId,
								vector<string>& listAssetId,
								int isEtkOrNot = 1)
	{
		HRESULT hr;
		MSXML2::IXMLDOMNode* aNode = NULL;
		MSXML2::IXMLDOMNodeList* resultList = NULL;
		long nbNodes = 0;
		ARM_Object* sec = NULL;

		try
		{
			//Informations contenues dans l'enveloppe du deal
			string chemin;

			if (isEtkOrNot)
			{
				chemin = ETKSTRING + GetSummitType() +"/Env/ENV";
				
				XMLDoc->selectSingleNode((_bstr_t)chemin.c_str(), &aNode);

				custId= GetStringFromXMLNode(aNode, "Cust"); 
				dealId = GetStringFromXMLNode(aNode, "DealId");
				bookName = GetStringFromXMLNode(aNode, "Book");

				if (aNode)
					aNode->Release();
				aNode = NULL;
			}
			else
			{
				custId		= "ZZZ"; 
				dealId		= "ZZZ"; 
				bookName	= "ZZZ"; 
			}

			string tradeId;
			string tradeType;

			if (isEtkOrNot)
			{
				chemin = ETKSTRING + GetSummitType();
				XMLDoc->selectSingleNode((_bstr_t)chemin.c_str(), &aNode);

				tradeId		= GetStringFromXMLNode(aNode, "TradeId"); 
				tradeType	= GetStringFromXMLNode(aNode, "TradeType");
			}
			else
			{			
				chemin = SUMMITSTRING + string((const char *) "ASSET");
				XMLDoc->selectSingleNode((_bstr_t)chemin.c_str(), &aNode);

				tradeId = GetStringFromXMLNode(aNode, "TradeId"); 
				tradeType = GetStringFromXMLNode(aNode, "Type");
			}

			if (aNode)
				aNode->Release();
			aNode = NULL;

			if (isEtkOrNot)
				chemin = ETKSTRING + GetSummitType() +"/Assets/ASSET";
			else
				chemin = SUMMITSTRING + string((const char *) "ASSET");
		
			if ( XMLDoc->selectNodes((_bstr_t)chemin.c_str(), &resultList) == S_OK )
			{
				AssetParserManager::Init(isEtkOrNot);

				hr = resultList->get_item(0, &aNode);
				
				if (( hr == S_OK ) && ( aNode != NULL))
				{
					string assetType = ConvertManager::GetAssetKey(aNode);

					AbstractAssetParser* p= AssetParserManager::GetAssetParser(assetType);
					
					if ( p ) 
					{
						p->SetTradeId(tradeId);
						p->SetTradeType(tradeType);
						string assetId;

						sec = p->BuildAsset(aNode, date,assetId,isEtkOrNot);
				
						listAssetId.push_back(assetId);

						if (!sec)
						{
							CCString msg("ParserManager::BuildInstrument,no loader for Instrument type ");
							msg += GetSummitType().c_str();
							throw Exception (__LINE__, __FILE__, ERR_INVALID_OPERATION, msg);
						}
					}
					if (aNode)
						aNode->Release();
					aNode = NULL;
				}
			}
			if (resultList)
				resultList->Release();
			resultList = NULL;
		}

		catch (...) 
		{
			throw;
		} 

		return((ARM_Object *) sec);
	}
};



//! \brief loader for Swap
class SwapLoader : public AbstractParser
{

public:

    SwapLoader()
	{
		SetSummitType("SWAP");
	}

   ~SwapLoader()
	{
	}

	ARM_Object* BuildInstrument(MSXML2::IXMLDOMDocument *XMLDoc,
								const ARM_Date& date,
								const string filter,
								string& bookName,
								string& structureId,
								string& custId,
								string& dealId,
								vector<string>& listAssetId,
								int isEtkOrNot = 1)
	{
		HRESULT hr;
		MSXML2::IXMLDOMNode* aNode = NULL;
		MSXML2::IXMLDOMNodeList* resultList = NULL;
		long nbNodes = 0;
		long nbAssets = 0;
		ARM_Swap* newSwap = NULL;
		vector<ARM_Security*> assets;

		try
		{
			// Informations contenues dans l'enveloppe du deal
			string chemin = ETKSTRING + GetSummitType() +"/Env/ENV";
			XMLDoc->selectSingleNode((_bstr_t)chemin.c_str(), &aNode);
			custId= GetStringFromXMLNode(aNode, "Cust");
			dealId = GetStringFromXMLNode(aNode, "DealId");
			bookName = GetStringFromXMLNode(aNode, "Book");
			if (aNode)
				aNode->Release();
			aNode = NULL;
			
			chemin = ETKSTRING + GetSummitType() +"/Assets/ASSET";
			if (XMLDoc->selectNodes((_bstr_t)chemin.c_str(), &resultList) == S_OK)
			{
				AssetParserManager::Init();

				resultList->get_length(&nbAssets);
				if (nbAssets != 2)
				{
					CCString msg((CCString)"Option XML string is not valid \n");

					throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
									 (const char *) msg);
				}

				for (long indexAsset=0 ; indexAsset<nbAssets ; indexAsset++)
				{
					hr=resultList->get_item(indexAsset, &aNode);
					if (hr==S_OK && aNode!=NULL)
					{	
						string assetType = ConvertManager::GetAssetKey(aNode);
						AbstractAssetParser* p= AssetParserManager::GetAssetParser(assetType);
						if ( p ) 
						{
							ARM_Object* sec = NULL;

							string assetId;
							sec = p->BuildAsset(aNode, date, assetId);
							if (sec)
							{
								assets.push_back(dynamic_cast<ARM_SwapLeg*>(sec));
								listAssetId.push_back(assetId);
							}
							else
							{
								throw Exception (__LINE__, __FILE__, ERR_INVALID_OPERATION,
								 "ParserManager::BuildInstrument,no loader for Instrument type %s", GetSummitType());
							}
						}
					}
					if (aNode)
						aNode->Release();
					aNode = NULL;
				}

				//CorridorLeg contained in the first asset.
				if (ARM_CorridorLeg* corridorleg = dynamic_cast<ARM_CorridorLeg*>(assets[0]))
				{
					newSwap = new ARM_Swap(dynamic_cast<ARM_SwapLeg*>(assets[0]), dynamic_cast<ARM_SwapLeg*>(assets[1]));
				}
				//CorridorLeg contained in the second asset.
				else if (ARM_CorridorLeg* corridorleg = dynamic_cast<ARM_CorridorLeg*>(assets[1]))
				{
					newSwap = new ARM_Swap(dynamic_cast<ARM_SwapLeg*>(assets[1]), dynamic_cast<ARM_SwapLeg*>(assets[0]));
				}
				//No corridor leg: just standard legs.
				else
				{
					newSwap = new ARM_Swap(dynamic_cast<ARM_SwapLeg*>(assets[0]), dynamic_cast<ARM_SwapLeg*>(assets[1]));
				}
			}

			if (resultList)
				resultList->Release();
			resultList = NULL;
		}

		catch (Exception& x) 
		{
	        x.DebugPrint();

			throw x;
		} 
		catch (...) 
		{
			throw;
		} 

		return (ARM_Object*) newSwap;
	}
};


//! \brief loader for FXOptionStrip
class FXOptionStripLoader : public AbstractParser
{

public:

    FXOptionStripLoader()
	{
		SetSummitType("CAPTR");
	}

	~FXOptionStripLoader()
	{
	}

	ARM_Object* BuildInstrument(MSXML2::IXMLDOMDocument *XMLDoc,
								const ARM_Date& date,
								const string filter,
								string& bookName,
								string& structureId,
								string& custId,
								string& dealId,
								vector<string>& listAssetId,
								int isEtkOrNot = 1)
	{
		HRESULT hr;
		MSXML2::IXMLDOMNode* aNode = NULL;
		MSXML2::IXMLDOMNodeList* resultList = NULL;
		long nbNodes = 0;
		long nbAssets = 0;
		ARM_Object* sec = NULL;

		try
		{
			// Informations contenues dans l'enveloppe du deal : BOOK, STRUCTURE, CUSTOMER, DEAL
			string chemin;

			if (isEtkOrNot)
			{
				chemin = ETKSTRING + GetSummitType() +"/Env/ENV";
				
				XMLDoc->selectSingleNode((_bstr_t)chemin.c_str(), &aNode);

				bookName	= GetStringFromXMLNode(aNode, "Book");
				structureId = GetStringFromXMLNode(aNode, "StructureId");
				custId		= GetStringFromXMLNode(aNode, "Cust"); 
				dealId		= GetStringFromXMLNode(aNode, "DealId");

				if (aNode)
					aNode->Release();
				aNode = NULL;
			}
			else
			{
				bookName	= "ZZZ"; 
				structureId = "ZZZ";
				custId		= "ZZZ"; 
				dealId		= "ZZZ"; 
			}

			string tradeId;
			string tradeType;

			if (isEtkOrNot)
			{
				chemin = ETKSTRING + GetSummitType();
				XMLDoc->selectSingleNode((_bstr_t)chemin.c_str(), &aNode);

				tradeId		= GetStringFromXMLNode(aNode, "TradeId"); 
				tradeType	= GetStringFromXMLNode(aNode, "TradeType"); 
			}
			else
			{			
				chemin = SUMMITSTRING + string((const char *) "ASSET");
				XMLDoc->selectSingleNode((_bstr_t)chemin.c_str(), &aNode);

				tradeId = GetStringFromXMLNode(aNode, "TradeId"); 
				tradeType = GetStringFromXMLNode(aNode, "Type");
			}

			if (aNode)
				aNode->Release();
			aNode = NULL;

			if (isEtkOrNot)
				chemin = ETKSTRING + GetSummitType() +"/Assets/ASSET";
			else
				chemin = SUMMITSTRING + string((const char *) "ASSET");

			if ( XMLDoc->selectNodes((_bstr_t)chemin.c_str(), &resultList) == S_OK )
			{
				AssetParserManager::Init(isEtkOrNot);

				hr = resultList->get_item(0, &aNode);
				
				if (( hr == S_OK ) && ( aNode != NULL))
				{
					string assetType = ConvertManager::GetAssetKey(aNode);

					AbstractAssetParser* p= AssetParserManager::GetAssetParser(assetType);
					
					if ( p ) 
					{
						p->SetTradeId(tradeId);
						p->SetTradeType(tradeType);

						string assetId;

						sec = p->BuildAsset(aNode, date,assetId,isEtkOrNot);
				
						listAssetId.push_back(assetId);

						if (!sec)
						{
							CCString msg("ParserManager::BuildInstrument,no loader for Instrument type ");
							msg += GetSummitType().c_str();
							throw Exception (__LINE__, __FILE__, ERR_INVALID_OPERATION, msg);
						}
					}
					if (aNode)
						aNode->Release();
					aNode = NULL;
				}
			}

			if (resultList)
				resultList->Release();
			resultList = NULL;
		}

		catch (Exception& x) 
		{
	        x.DebugPrint();

			throw x;
		} 
		catch (...) 
		{
			throw;
		} 

		return (ARM_Object*) sec;
	}
};


//! \brief loader for FXOptionStrip
class FXStripLoader : public AbstractParser
{

public:

    FXStripLoader()
	{
		SetSummitType("CAPTR");
	}

	~FXStripLoader()
	{
	}

	ARM_Object* BuildInstrument(IXMLDOMDocument *XMLDoc,
								const ARM_Date& date,
								const string filter,
								string& bookName,
								string& structureId,
								string& custId,
								string& dealId,
								vector<string>& listAssetId,
								int isEtkOrNot = 1)
	{
		HRESULT hr;
		IXMLDOMNode* aNode = NULL;
		IXMLDOMNodeList* resultList = NULL;
		long nbNodes = 0;
		long nbAssets = 0;
		ARM_Object* sec = NULL;

		try
		{
			// Informations contenues dans l'enveloppe du deal : BOOK, STRUCTURE, CUSTOMER, DEAL
			string chemin;

			if (isEtkOrNot)
			{
				chemin = ETKSTRING + GetSummitType() +"/Env/ENV";
				
				XMLDoc->selectSingleNode((_bstr_t)chemin.c_str(), &aNode);

				bookName	= GetStringFromXMLNode(aNode, "Book");
				structureId = GetStringFromXMLNode(aNode, "StructureId");
				custId		= GetStringFromXMLNode(aNode, "Cust"); 
				dealId		= GetStringFromXMLNode(aNode, "DealId");

				if (aNode)
					aNode->Release();
				aNode = NULL;
			}
			else
			{
				bookName	= "ZZZ"; 
				structureId = "ZZZ";
				custId		= "ZZZ"; 
				dealId		= "ZZZ"; 
			}

			string tradeId;
			string tradeType;

			if (isEtkOrNot)
			{
				chemin = ETKSTRING + GetSummitType();
				XMLDoc->selectSingleNode((_bstr_t)chemin.c_str(), &aNode);

				tradeId		= GetStringFromXMLNode(aNode, "TradeId"); 
				tradeType	= GetStringFromXMLNode(aNode, "TradeType"); 
			}
			else
			{			
				chemin = SUMMITSTRING + string((const char *) "ASSET");
				XMLDoc->selectSingleNode((_bstr_t)chemin.c_str(), &aNode);

				tradeId = GetStringFromXMLNode(aNode, "TradeId"); 
				tradeType = GetStringFromXMLNode(aNode, "Type");
			}

			if (aNode)
				aNode->Release();
			aNode = NULL;

			if (isEtkOrNot)
				chemin = ETKSTRING + GetSummitType() +"/Assets/ASSET";
			else
				chemin = SUMMITSTRING + string((const char *) "ASSET");

			if ( XMLDoc->selectNodes((_bstr_t)chemin.c_str(), &resultList) == S_OK )
			{
				AssetParserManager::Init(isEtkOrNot);

				hr = resultList->get_item(0, &aNode);
				
				if (( hr == S_OK ) && ( aNode != NULL))
				{
					string assetType = ConvertManager::GetAssetKey(aNode);

					AbstractAssetParser* p= AssetParserManager::GetAssetParser(assetType);
					
					if ( p ) 
					{
						p->SetTradeId(tradeId);
						p->SetTradeType(tradeType);

						string assetId;

						sec = p->BuildAsset(aNode, date,assetId,isEtkOrNot);
				
						listAssetId.push_back(assetId);

						if (!sec)
						{
							CCString msg("ParserManager::BuildInstrument,no loader for Instrument type ");
							msg += GetSummitType().c_str();
							throw Exception (__LINE__, __FILE__, ERR_INVALID_OPERATION, msg);
						}
					}
					if (aNode)
						aNode->Release();
					aNode = NULL;
				}
			}

			if (resultList)
				resultList->Release();
			resultList = NULL;
		}

		catch (Exception& x) 
		{
	        x.DebugPrint();

			throw x;
		} 
		catch (...) 
		{
			throw;
		} 

		return (ARM_Object*) sec;
	}
};


//! \brief loader for FXOptionStrip with calculator
class FXStripCalculatorLoader : public AbstractParser
{

public:

    FXStripCalculatorLoader()
	{
		SetSummitType("CAPTR");
	}

	~FXStripCalculatorLoader()
	{
	}

	ARM_Object* BuildInstrument(IXMLDOMDocument *XMLDoc,
								const ARM_Date& date,
								const string filter,
								string& bookName,
								string& structureId,
								string& custId,
								string& dealId,
								vector<string>& listAssetId,
								int isEtkOrNot = 1)
	{
		HRESULT hr;
		IXMLDOMNode* aNode = NULL;
		IXMLDOMNodeList* resultList = NULL;
		long nbNodes = 0;
		long nbAssets = 0;
		ARM_Object* sec = NULL;

		try
		{
			// Informations contenues dans l'enveloppe du deal : BOOK, STRUCTURE, CUSTOMER, DEAL
			string chemin;

			if (isEtkOrNot)
			{
				chemin = ETKSTRING + GetSummitType() +"/Env/ENV";
				
				XMLDoc->selectSingleNode((_bstr_t)chemin.c_str(), &aNode);

				bookName	= GetStringFromXMLNode(aNode, "Book");
				structureId = GetStringFromXMLNode(aNode, "StructureId");
				custId		= GetStringFromXMLNode(aNode, "Cust"); 
				dealId		= GetStringFromXMLNode(aNode, "DealId");

				if (aNode)
					aNode->Release();
				aNode = NULL;
			}
			else
			{
				bookName	= "ZZZ"; 
				structureId = "ZZZ";
				custId		= "ZZZ"; 
				dealId		= "ZZZ"; 
			}

			string tradeId;
			string tradeType;

			if (isEtkOrNot)
			{
				chemin = ETKSTRING + GetSummitType();
				XMLDoc->selectSingleNode((_bstr_t)chemin.c_str(), &aNode);

				tradeId		= GetStringFromXMLNode(aNode, "TradeId"); 
				tradeType	= GetStringFromXMLNode(aNode, "TradeType"); 
			}
			else
			{			
				chemin = SUMMITSTRING + string((const char *) "ASSET");
				XMLDoc->selectSingleNode((_bstr_t)chemin.c_str(), &aNode);

				tradeId = GetStringFromXMLNode(aNode, "TradeId"); 
				tradeType = GetStringFromXMLNode(aNode, "Type");
			}

			if (aNode)
				aNode->Release();
			aNode = NULL;

			if (isEtkOrNot)
				chemin = ETKSTRING + GetSummitType() +"/Assets/ASSET";
			else
				chemin = SUMMITSTRING + string((const char *) "ASSET");

			if ( XMLDoc->selectNodes((_bstr_t)chemin.c_str(), &resultList) == S_OK )
			{
				AssetParserManager::Init(isEtkOrNot);

				hr = resultList->get_item(0, &aNode);
				
				if (( hr == S_OK ) && ( aNode != NULL))
				{
					string assetType = ConvertManager::GetAssetKey(aNode);
					
					assetType = "C." + assetType;

					AbstractAssetParser* p= AssetParserManager::GetAssetParser(assetType);
					
					if ( p ) 
					{
						p->SetTradeId(tradeId);
						p->SetTradeType(tradeType);

						string assetId;

						sec = p->BuildAsset(aNode, date,assetId,isEtkOrNot);
				
						listAssetId.push_back(assetId);

						if (!sec)
						{
							CCString msg("ParserManager::BuildInstrument,no loader for Instrument type ");
							msg += GetSummitType().c_str();
							throw Exception (__LINE__, __FILE__, ERR_INVALID_OPERATION, msg);
						}
					}
					if (aNode)
						aNode->Release();
					aNode = NULL;
				}
			}

			if (resultList)
				resultList->Release();
			resultList = NULL;
		}

		catch (Exception& x) 
		{
	        x.DebugPrint();

			throw x;
		} 
		catch (...) 
		{
			throw;
		} 

		return (ARM_Object*) sec;
	}
};


//! \brief loader for RNGDouble
class CorridorDblConditionLoader : public AbstractParser
{

public:

    CorridorDblConditionLoader()
	{
		//BPU : A checker !!!!
		SetSummitType("CAPTR");
	}

	~CorridorDblConditionLoader()
	{
	}

	ARM_Object* BuildInstrument(IXMLDOMDocument *XMLDoc,
								const ARM_Date& date,
								const string filter,
								string& bookName,
								string& structureId,
								string& custId,
								string& dealId,
								vector<string>& listAssetId,
								int isEtkOrNot = 1)
	{
		HRESULT hr;
		IXMLDOMNode* aNode = NULL;
		IXMLDOMNodeList* resultList = NULL;
		long nbNodes = 0;
		long nbAssets = 0;
		ARM_Object* sec = NULL;

		try
		{
			// Informations contenues dans l'enveloppe du deal : BOOK, STRUCTURE, CUSTOMER, DEAL

			//BPU : A FAIRE
			string chemin;

			if (isEtkOrNot)
			{
				chemin = ETKSTRING + GetSummitType() +"/Env/ENV";
				
				XMLDoc->selectSingleNode((_bstr_t)chemin.c_str(), &aNode);

				bookName	= GetStringFromXMLNode(aNode, "Book");
				structureId = GetStringFromXMLNode(aNode, "StructureId");
				custId		= GetStringFromXMLNode(aNode, "Cust"); 
				dealId		= GetStringFromXMLNode(aNode, "DealId");

				if (aNode)
					aNode->Release();
				aNode = NULL;
			}
			else
			{
				bookName	= "ZZZ"; 
				structureId = "ZZZ";
				custId		= "ZZZ"; 
				dealId		= "ZZZ"; 
			}

			string tradeId;
			string tradeType;

			if (isEtkOrNot)
			{
				chemin = ETKSTRING + GetSummitType();
				XMLDoc->selectSingleNode((_bstr_t)chemin.c_str(), &aNode);

				tradeId		= GetStringFromXMLNode(aNode, "TradeId"); 
				tradeType	= GetStringFromXMLNode(aNode, "TradeType"); 
			}
			else
			{			
				chemin = SUMMITSTRING + string((const char *) "ASSET");
				XMLDoc->selectSingleNode((_bstr_t)chemin.c_str(), &aNode);

				tradeId = GetStringFromXMLNode(aNode, "TradeId"); 
				tradeType = GetStringFromXMLNode(aNode, "Type");
			}

			if (aNode)
				aNode->Release();
			aNode = NULL;

			if (isEtkOrNot)
				chemin = ETKSTRING + GetSummitType() +"/Assets/ASSET";
			else
				chemin = SUMMITSTRING + string((const char *) "ASSET");

			if ( XMLDoc->selectNodes((_bstr_t)chemin.c_str(), &resultList) == S_OK )
			{
				AssetParserManager::Init(isEtkOrNot);

				hr = resultList->get_item(0, &aNode);
				
				if (( hr == S_OK ) && ( aNode != NULL))
				{
					string assetType = ConvertManager::GetAssetKey(aNode);

					AbstractAssetParser* p= AssetParserManager::GetAssetParser(assetType);
					
					if ( p ) 
					{
						p->SetTradeId(tradeId);
						p->SetTradeType(tradeType);

						string assetId;

						sec = p->BuildAsset(aNode, date,assetId,isEtkOrNot);
				
						listAssetId.push_back(assetId);

						if (!sec)
						{
							CCString msg("ParserManager::BuildInstrument,no loader for Instrument type ");
							msg += GetSummitType().c_str();
							throw Exception (__LINE__, __FILE__, ERR_INVALID_OPERATION, msg);
						}
					}
					if (aNode)
						aNode->Release();
					aNode = NULL;
				}
			}

			if (resultList)
				resultList->Release();
			resultList = NULL;
		}

		catch (Exception& x) 
		{
	        x.DebugPrint();

			throw x;
		} 
		catch (...) 
		{
			throw;
		} 

		return (ARM_Object*) sec;
	}
};

void ParserManager::Init(int isEtkOrNot) 
{
	if ( !initDone ) 
	{
		isEtkOrSummit = isEtkOrNot;
		RegisterParser(string("IRG"), new CapLoader());
		RegisterParser(string("SPDOPT"), new SpreadOptionLoader());
		RegisterParser(string("EXOTIC"), new ExoticLoader());
		RegisterParser(string("GLOBALEXOTIC"), new GlobalExoticLoader());
		RegisterParser(string("MATURITYCAP"), new MaturityCapLoader());
		RegisterParser(string("RFTARN"), new TarnLoader());
		//RegisterParser(string("CSB"), new CallableSnowBallLoader());
		RegisterParser(string("SWAP"), new SwapLoader());
		RegisterParser(string("MEMORYSO"), new MemorySOLoader());
		RegisterParser(string("MEMORYCAP"), new MemoryCapLoader());
		RegisterParser(string("ALMCAPTION"), new CaptionLoader());
		RegisterParser(string("FXOPTSTRIP"), new FXOptionStripLoader());
		RegisterParser(string("FXSTRIP"), new FXStripLoader());
		RegisterParser(string("CFXSTRIP"), new FXStripCalculatorLoader());
		RegisterParser(string("RNG_DOUBLE"), new CorridorDblConditionLoader());
		initDone = true;
		ConvertManager::Init();
	}
	isEtkOrSummit = isEtkOrNot;
}


void ConvertManager::Init(void)
{
	if ( !initDone ) 
	{
		RegisterKey(string("SWOPT.RNGLIBOR_OPTION"), string("CRA"));
		RegisterKey(string("SWOPT.OPT_PRCS"), string("PRCS"));
		RegisterKey(string("SWOPT.OPT_PRCS_DUALE"), string("PRCS"));
		RegisterKey(string("SWOPT.OPT_PRCS_DUALE_MANDAT"), string("PRCS"));
		RegisterKey(string("SWOPT.OPT_PRCS_DUALE_SPOT"), string("PRCS"));
		RegisterKey(string("SWOPT.OPT_PRDC"), string("PRCS"));
		RegisterKey(string("SWOPT.OPT_PRDC_DUALE"), string("PRCS"));
		RegisterKey(string("SWOPT.OPT_PRDC_DUALE_MANDAT"), string("PRCS"));
		RegisterKey(string("SWOPT.OPT_PRDC_DUALE_SPOT"), string("PRCS"));
		RegisterKey(string("SWOPT"), string("SWOPT"));
		RegisterKey(string("FXOPT"), string("FXOPT"));
		RegisterKey(string("FXOPT.FX_OPT_CDC"), string("FXOPT"));
		RegisterKey(string("EXOTIC"), string("EXOTIC"));
		RegisterKey(string("IRG.SPREADOPTIONLOG"), string("SPDOPT"));
		RegisterKey(string("IRG.SPREADOPTIONLOGDIGITAL"), string("SPDOPT"));
		RegisterKey(string("IRG.SPREADOPTIONLOGFLTDIGITAL"), string("SPDOPT"));
		RegisterKey(string("IRG.SPREADOPTIONLOGCORRIDOR"), string("SPDOPT"));
		RegisterKey(string("MATURITYCAP"), string("MATURITYCAP"));
		RegisterKey(string("IRG.MATURITYCAP"), string("MATURITYCAP"));
		RegisterKey(string("IRG.MATURITY_CAP"), string("MATURITYCAP"));
		RegisterKey(string("IRG.MEMORYSO"), string("MEMORYSO"));
		RegisterKey(string("IRG.MEMORYCAP"), string("MEMORYCAP"));
		// 
		RegisterKey(string("EXOTIC.CREDIT_SWAP.FTD"), string("FTD"));
		RegisterKey(string("EXOTIC.CREDIT_SWAP.STCDO"), string("TRANCHE"));
		RegisterKey(string("EXOTIC.CREDIT_SWAP.STCDOCFWD"), string("TRANCHE"));
		RegisterKey(string("EXOTIC.CREDIT_SWAP.STCDOFWD"), string("TRANCHE"));
		RegisterKey(string("EXOTIC.CREDIT_SWAP.RISKY_LEG"), string("TRANCHE"));
		RegisterKey(string("EXOTIC.CREDIT_SWAP.INDEXCDO"), string("INDEXCDO"));
		RegisterKey(string("EXOTIC.CREDIT_SWAP.CDOSQUARE"), string("CDO2"));
		// 
		RegisterKey(string("EXOTIC.CREDIT_SWAP.CDSCEURO"), string("CDS"));
		RegisterKey(string("EXOTIC.CREDIT_SWAP.CDSCUS"), string("CDS"));
		RegisterKey(string("EXOTIC.CREDIT_SWAP.CDSCUSHY"), string("CDS"));
		RegisterKey(string("EXOTIC.CREDIT_SWAP.CDSCASIA"), string("CDS"));
		RegisterKey(string("EXOTIC.CREDIT_SWAP.CDSCJAPAN"), string("CDS"));
		RegisterKey(string("EXOTIC.CREDIT_SWAP.CDSCAUS"), string("CDS"));
		RegisterKey(string("EXOTIC.CREDIT_SWAP.CDSCSING"), string("CDS"));
		RegisterKey(string("EXOTIC.CREDIT_SWAP.CDSCNEWZ"), string("CDS"));
		RegisterKey(string("EXOTIC.CREDIT_SWAP.CDSCASIA"), string("CDS"));
		RegisterKey(string("EXOTIC.CREDIT_SWAP.CDSCISUB"), string("CDS"));
		RegisterKey(string("EXOTIC.CREDIT_SWAP.CDSCMONO"), string("CDS"));
		RegisterKey(string("EXOTIC.CREDIT_SWAP.CDSSASIA"), string("CDS"));
		RegisterKey(string("EXOTIC.CREDIT_SWAP.CDSCLATAM"), string("CDS"));
		RegisterKey(string("EXOTIC.CREDIT_SWAP.CDSSCEEUR"), string("CDS"));
		RegisterKey(string("EXOTIC.CREDIT_SWAP.CDSSJAPAN"), string("CDS"));
		RegisterKey(string("EXOTIC.CREDIT_SWAP.CDSSLATAME"), string("CDS"));
		RegisterKey(string("EXOTIC.CREDIT_SWAP.CDSSME"), string("CDS"));
		RegisterKey(string("EXOTIC.CREDIT_SWAP.CDSSWESEUR"), string("CDS"));
		RegisterKey(string("EXOTIC.CREDIT_SWAP.CDSSAUS"), string("CDS"));
		RegisterKey(string("EXOTIC.CREDIT_SWAP.CDSSNZ"), string("CDS"));
		RegisterKey(string("EXOTIC.CREDIT_SWAP.CDSSSING"), string("CDS"));
		RegisterKey(string("EXOTIC.CREDIT_SWAP.CDS_ITRAXX"), string("CDS"));
		RegisterKey(string("EXOTIC.CREDIT_SWAP.CDS_CDX"), string("CDS"));
		RegisterKey(string("EXOTIC.CREDIT_SWAP.EDS"), string("CDS"));
		RegisterKey(string("EXOTIC.CREDIT_SWAP.CDS_RECOV"), string("CDS"));
		RegisterKey(string("EXOTIC.CREDIT_SWAP.CDSSTRUCT"), string("CDS"));
		// 
		RegisterKey(string("SWAP"), string("SWAP"));
		RegisterKey(string("SWAP.RCMS"), string("SWAP"));
		RegisterKey(string("SWAP.RCMSD"), string("SWAP"));
		RegisterKey(string("SWAP.RCMT"), string("SWAP"));
		RegisterKey(string("SWAP.RCMTD"), string("SWAP"));
		RegisterKey(string("SWAP.TEC_SPLK"), string("SWAP"));
		RegisterKey(string("SWAP.RNGLIBORFIX"), string("SWAP"));
		RegisterKey(string("SWAP.RNGLIBORLIBOR"), string("SWAP"));
		RegisterKey(string("SWAP.LIVRET_A"), string("SWAP"));
		RegisterKey(string("SWAP.LIVRET_A_ARRONDI"), string("SWAP"));
		RegisterKey(string("IRG.RNGLIBORFIX"), string("IRG"));
		RegisterKey(string("IRG.RNGLIBORLIBOR"), string("IRG"));	
		RegisterKey(string("IRG.TEC_SPLK_CF"), string("IRG"));
		RegisterKey(string("IRG.RCMS"), string("IRG"));
		RegisterKey(string("IRG.RCMSD"), string("IRG"));
		RegisterKey(string("IRG.RCMT"), string("IRG"));
		RegisterKey(string("IRG.RCMTD"), string("IRG"));
		RegisterKey(string("IRG"), string("IRG"));
		RegisterKey(string("IRG.LIVRET_A"), string("IRG"));
		RegisterKey(string("IRG.LIVRET_A_ARRONDI"), string("IRG"));
		RegisterKey(string("IRG.CUO_SMILE"), string("IRG"));
		RegisterKey(string("IRG.CUO_SPLK"), string("IRG"));
		RegisterKey(string("IRG.CUI_SMILE"), string("IRG"));
		RegisterKey(string("IRG.CUI_SPLK"), string("IRG"));
		RegisterKey(string("IRG.PDO_SMILE"), string("IRG"));
		RegisterKey(string("IRG.PDO_SPLK"), string("IRG"));
		RegisterKey(string("IRG.PDI_SMILE"), string("IRG"));
		RegisterKey(string("IRG.PDI_SPLK"), string("IRG"));
		RegisterKey(string("IRG.DIGITAL_SMILE"), string("IRG"));
		RegisterKey(string("IRG.DIGITAL_SPLK"), string("IRG"));
		RegisterKey(string("IRG.GLOBALCAP"), string("IRG"));

		RegisterKey(string("RFTARN"), string("RFTARN"));
		RegisterKey(string("IRG.RFTARN"), string("RFTARN"));

		RegisterKey(string("IRG.ALMCAPTION"), string("ALMCAPTION"));
		RegisterKey(string("IRG.FX_CALL_STRIP"), string("FXOPTSTRIP"));
		RegisterKey(string("IRG.FX_CALL_STRIP_DIGITAL"), string("FXOPTSTRIP"));

		RegisterKey(string("IRG.FX_QUANTO_STRIP"), string("FXSTRIP"));
		RegisterKey(string("IRG.FX_QUANTO_STRIP_DIGITAL"), string("FXSTRIP"));
		RegisterKey(string("IRG.FX_SPREAD_STRIP"), string("FXSTRIP"));
		RegisterKey(string("IRG.FX_SPREAD_STRIP_DIGITAL"), string("FXSTRIP"));

		//RegisterKey(string("IRG.FX_QUANTO_STRIP"), string("CFXSTRIP"));
		//RegisterKey(string("IRG.FX_QUANTO_STRIP_DIGITAL"), string("CFXSTRIP"));
		//RegisterKey(string("IRG.FX_SPREAD_STRIP"), string("CFXSTRIP"));
		//RegisterKey(string("IRG.FX_SPREAD_STRIP_DIGITAL"), string("CFXSTRIP"));
		
		RegisterKey(string("RNG_DOUBLE"), string("RNG_DOUBLE"));
		RegisterKey(string("SWAP.RNG_DOUBLE.FORM"), string("RNG_DOUBLE"));
		RegisterKey(string("EXOTIC.RNG_DOUBLE"), string("RNG_DOUBLE"));

		initDone = true;
	}
}


void ConvertManager::Release(void)
{
	if ( !initDone ) 
	{
		throw Exception (__LINE__, __FILE__, ERR_INVALID_OPERATION,
			 "ConvertManager::Release, ConvertManager::Init not done");
	}

}

string ConvertManager::BuildSummitKey(MSXML2::IXMLDOMNode *theNode, int isEtkOrNot)
{
	string key;
	MSXML2::IXMLDOMNode* aNode = NULL;
	BSTR resultat = NULL;
	
	key= GetStringFromXMLNode(theNode,"*/*/TradeType");

	if ( strcmp(key.c_str(),"EXOTIC") != 0 )
	{
		string formulae;

		if (isEtkOrNot)
		   formulae = GetStringFromXMLNode(theNode,"*/*/*/*/ProductName");
		else
		   formulae = GetStringFromXMLNode(theNode,"*/PRODUCTNAME");

		if ( formulae != "" )
		   key = key + "." + formulae;
	}
	else 
	{
		// JLA: additional Credit Typing Credit typin,g		
		if (GetStringFromXMLNode(theNode,"*/*/*/*/ProductName")=="CREDIT_SWAP") 
			key += std::string(".CREDIT_SWAP.") +  GetStringFromXMLNode(theNode,"*/*/*/*/ProductType") ;
		
	}

	return key;
}


string ConvertManager::GetAssetKey(MSXML2::IXMLDOMNode* XMLNode)
{
	string key;
	string sType;
	string sProductName;
	string sInterestDmIndex;

	BSTR resultat = NULL;

	sType = GetStringFromXMLNode(XMLNode,"Type");
	sProductName = GetStringFromXMLNode(XMLNode,"ProductName");
	sInterestDmIndex = GetStringFromXMLNode(XMLNode,"INTEREST_dmIndex");

	key = sType + "." + sProductName + "." + sInterestDmIndex;

	return key;
}


string ConvertManager::ConvertETKKey(string& SummitType)
{
	if ( !initDone ) 
	{
		throw Exception (__LINE__, __FILE__, ERR_INVALID_OPERATION,
			 "ConvertManager::ConvertKey, ConvertManager::Init not done");
	}

	map<string, string>::iterator iter = itsConvertTable.find(SummitType);	
	if ( iter != itsConvertTable.end() )
	{
		return iter->second;
	}
	else 
	{
		throw Exception (__LINE__, __FILE__, ERR_INVALID_OPERATION,
			 "ConvertManager::ConvertKey, no ETK type for Summit type %s",SummitType);		
	}	
}

		
string ConvertManager::GetETKKey(MSXML2::IXMLDOMNode* theNode, int isEtkOrNot)
{
	return ConvertETKKey(BuildSummitKey(theNode,isEtkOrNot));
}


bool ConvertManager::RegisterKey(string& SummitType, string& Type)
{
	pair<string, string> p(SummitType, Type);
	pair<map<string, string>::iterator, bool> r = itsConvertTable.insert(p);

	return r.second;
}




/*---------------------------------------------------------------------------*/
/*---- End of file ----*/
