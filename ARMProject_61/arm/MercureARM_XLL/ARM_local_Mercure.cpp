#pragma warning(disable : 4541)

#include "firstToBeIncluded.h"
#include <ARM\libarm\ARM_result.h>
#include <ARM\libarm_local\ARM_local_class.h>
#include <ARM\local_xlarm\ARM_local_interglob.h>

#include "ARM_local_persistent.h"
#include "ARM_local_glob.h"
#include "ARM_local_class.h"

#include "ARM_MercureResult.h"
#include "ARM_MercureMarketDataManager.h"
#include "ARM_local_Mercure.h"

#include "Mercure.h"
#include "MarketDataManager.h"
#include "MarketDataDictionary.h"
#include "MetaMarketData.h"
#include "MetaZeroCurve.h"
#include "MetaVolatilityCurve.h"
#include "MetaVolatilityCube.h"
#include "MetaScalarData.h"
#include "MetaFXVolatilityCurve.h"
#include "Utility.h"

#include "zeroint.h"
#include "currency.h"
#include "volint.h"
#include "volcube.h"
#include "armscalardata.h"
#include "portfolio.h"

#ifdef __STDC__
#undef __STDC__
#endif
#include "BasicParser.h"
//#define __STDC__

// Ne pas mettre au-dessus
#include "MercureServerUtil.h"
#include "ARM_local_etoolkit.h"


using namespace mercure;


long ARMMercure_Hedge(vector<CCString>& C_TradeId,
					  CCString C_TradeType,
					  double   C_AsOfDate,
					  vector<CCString>& C_FallBack,
					  CCString C_MarketDataManagerId,
					  CCString C_HedgeRatioFilesDirectory,
					  CCString C_ConfigFile,
					  CCString C_MarketDataFile,
					  bool traceSecurityAndModel,
					  CCString C_HedgeRatioFile,
					  ARM_result& result,
					  long objId)
{
	long	MercureResultId;

    ARM_MercureResult*	vMercureResult = NULL;
	ARM_MercureResult*	vNewMercureResult = NULL;
 
	if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0)
	{
		result.setMsg ("ARM_ERR: Pb accessing objects");
		return ARM_KO;
	}

	char	sAsOfDate[11];

	CCString msg ("");

    try
    {
		Local_XLDATE2ARMDATE(C_AsOfDate, sAsOfDate);
		string		vAsOfDate(sAsOfDate);
		ARM_Date	vAsOfDateARM(sAsOfDate);

		MercureKernelService	vMercure;
		string	vHedgeRatioFilesDirectory(C_HedgeRatioFilesDirectory);
		LoadAndSetHedgeRatioFiles(vHedgeRatioFilesDirectory, vMercure, string(C_TradeType), string(C_HedgeRatioFile));

		vector<string>	vFallBack;

		vector<CCString>::iterator	it = C_FallBack.begin();
		while( it != C_FallBack.end() )
		{
			vFallBack.push_back(string(*it));
			it++;
		}

		string	vConfigFile(C_ConfigFile);
		string	configFile = LoadConfigFileAsXMLString(vConfigFile);
		vMercure.SetConfigFile(configFile);

		string	vMarketDataFile(C_MarketDataFile);
		string	marketDataFile = LoadInitialMarketDataFileAsXMLString(vMarketDataFile);
		vMercure.SetInitialMarketDataFile(marketDataFile);
		vMercure.GetMarketDataManager().LoadInitialMarketData(vAsOfDateARM);

		if(C_MarketDataManagerId.GetLen() != 0)
		{
			long	vMktDataMgrId = LocalGetNumObjectId(C_MarketDataManagerId);
			ARM_MercureMarketDataManager*	vMktDataMgr = (ARM_MercureMarketDataManager*) 
														LOCAL_PERSISTENT_OBJECTS->GetARMObject(vMktDataMgrId);

			if(vMktDataMgr)
			{
				map<string, MarketDataDictionary*>&	vDictionaries = vMktDataMgr->GetMarketDataManager()->GetMktDictionaries();
				map<string, MarketDataDictionary*>::iterator	vDicoIter = vDictionaries.find(vAsOfDate);

				if( vDicoIter != vDictionaries.end() )
				{
					MarketDataDictionary*	vMktDataDico = vDicoIter->second;
					MarketDataDictionary*	vMercureMDD = vMercure.GetMarketDataManager().GetMktDictionaries().find(vAsOfDate)->second;
					map<string, MetaMarketData*>::iterator	iter = vMktDataDico->GetDictionary().begin();
					map<string, MetaMarketData*>::iterator	end = vMktDataDico->GetDictionary().end();

					for(; iter != end; iter++)
					{
						string	vKey(iter->first);
						vMercureMDD->InsertMarketData(vKey, iter->second);
					}
				}
			}
		}

		string	calculationNumber("MercureARM_XLL");

		vector<string>*	vXMLOutput = new vector<string>;
		vector<string>*	vTXTOutput = new vector<string>;

		if(C_TradeType == "ARM_Security")	// ARM Object
		{
			string	sXMLOutput;
			string	sTXTOutput;
			long vARMObjectId;
			ARM_Object* vSecurity = NULL;
			string vTradeId = string("Excel");
			if (C_TradeId.size()==1)
			{
				vARMObjectId = LocalGetNumObjectId(C_TradeId[0]);
				vSecurity = (ARM_Object*) LOCAL_PERSISTENT_OBJECTS->GetARMObject(vARMObjectId)->Clone();
			}
			else if (C_TradeId.size()==2)
			{
				vTradeId = CCSTringToSTLString(C_TradeId[0]);
				vARMObjectId = LocalGetNumObjectId(C_TradeId[1]);
				vSecurity = (ARM_Object*) LOCAL_PERSISTENT_OBJECTS->GetARMObject(vARMObjectId)->Clone();
			}

			vMercure.Compute(vAsOfDateARM, vTradeId, vSecurity, vFallBack, calculationNumber, sXMLOutput, true, &sTXTOutput);
			vXMLOutput->push_back(sXMLOutput);
			vTXTOutput->push_back(sTXTOutput);
		}
		else if(C_TradeType == "Portfolio")
		{
			vector<string> tradeId;
			ARM_Portfolio* vPortfolio = new ARM_Portfolio(0);
			for(int i=0; i<C_TradeId.size(); i++)
			{
				tradeId.push_back(CCSTringToSTLString(C_TradeId[i++]));

				long vARMObjectId = LocalGetNumObjectId(C_TradeId[i]);
				ARM_Security* vSecurity = (ARM_Security*) LOCAL_PERSISTENT_OBJECTS->GetARMObject(vARMObjectId)->Clone();
				vPortfolio->AddInstrument(vSecurity, 0, 0, 0.0001);
			}

			vMercure.Compute(vAsOfDateARM, tradeId, vPortfolio, vFallBack, calculationNumber, vXMLOutput, true, vTXTOutput);
			delete vPortfolio;
		}
		else
		{
			vector<Deal>	deals;
			vector<CCString>::iterator iter_TradeId = C_TradeId.begin();
			vector<CCString>::iterator iter_TradeId_end = C_TradeId.end();
			for( ; iter_TradeId != iter_TradeId_end; iter_TradeId++)
			{
				deals.push_back(Deal(string(*iter_TradeId), string(C_TradeType)));
			}

			vMercure.Compute(vAsOfDateARM, deals, vFallBack, calculationNumber, vXMLOutput, true, vTXTOutput, traceSecurityAndModel);
		}
		// Fin Calcul


		if(objId == -1)
		{
			// Création d'un nouveau MercureResult
			vNewMercureResult = new ARM_MercureResult(vXMLOutput, vTXTOutput);

			if(vNewMercureResult == NULL)
			{
				result.setMsg ("ARM_ERR: Mercure Result is null");
				return ARM_KO;
			}

			CREATE_GLOBAL_OBJECT();
			
			MercureResultId = LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)vNewMercureResult);

			if (MercureResultId == RET_KO)
			{
				delete	vNewMercureResult;
				vNewMercureResult = NULL;

				result.setMsg ("ARM_ERR: Pb inserting object");				
				return ARM_KO;
			}

			result.setLong(MercureResultId);
			result.setString( vXMLOutput->begin()->c_str() );

			return ARM_OK;
		}
		else if(objId == -2)	// Hedge array : pas besoin de MercureResult
		{
			vector<string>::iterator iter_XML = vXMLOutput->begin();
			vector<string>::iterator iter_XML_end = vXMLOutput->end();
			for( ;iter_XML != iter_XML_end ; iter_XML++)
				result.setStringInVect( iter_XML->c_str() );

			delete	vXMLOutput;
			delete	vTXTOutput;

			return ARM_OK;
		}
		else
		{
			vMercureResult = (ARM_MercureResult*) LOCAL_PERSISTENT_OBJECTS->GetARMObject(objId);

			if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(vMercureResult, ARM_MERCURE_RESULT) == 1)
			{
				if (vMercureResult)
				{
					delete	vMercureResult;
					vMercureResult = NULL;
				}

				// Création d'un nouveau MercureResult
				vNewMercureResult = new ARM_MercureResult(vXMLOutput, vTXTOutput);

				if(vNewMercureResult == NULL)
				{
					result.setMsg ("ARM_ERR: Mercure Result is null");
					return ARM_KO;
				}

				LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_MercureResult*)vNewMercureResult, objId);

				result.setString( vXMLOutput->begin()->c_str() );
				return ARM_OK;
			}

			else
			{
				result.setMsg ("ARM_ERR: previous object is not of a good type");
				return	ARM_KO;
			}
		}
    }

    catch(Exception& x)
    {
        x.DebugPrint();
 
		delete	vNewMercureResult;
		vNewMercureResult = NULL;

		ARM_RESULT();
    } 

	catch(ParsingException& vException)
	{
		result.setMsg (vException.GetMessage().c_str());
		return	ARM_KO;
	}

	catch(...)
	{
		result.setMsg("Unknown exception caught in ARMMercure_Hedge");
		return	ARM_KO;
	}
}


long ARMCreateMarketDataManager(vector<CCString>& C_MarketDataIds,
								double   C_AsOfDate,
								vector<CCString>& C_FallBack,
								ARM_result& result,
								long objId)
{
	long	vMDMId;

    ARM_MercureMarketDataManager*	vMMDM = NULL;
	ARM_MercureMarketDataManager*	vNewMMDM = NULL;
 
	if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0)
	{
		result.setMsg ("ARM_ERR: Pb accessing objects");
		return ARM_KO;
	}

	CCString	msg ("");
	char	sAsOfDate[11];
	string	vKey;

	MetaMarketData*	vMetaMarketData = NULL;

    try
    {
		Local_XLDATE2ARMDATE(C_AsOfDate, sAsOfDate);
		string		vAsOfDate(sAsOfDate);
		ARM_Date	vAsOfDateARM(sAsOfDate);

		vector<string>	vFallBack;
		vector<CCString>::iterator	iter = C_FallBack.begin();
		while( iter != C_FallBack.end() )
		{
			vFallBack.push_back(string(*iter));
			iter++;
		}

		MarketDataManager*	vMktDataMgr = new MarketDataManager();
		vMktDataMgr->SetDefaultFallback(vFallBack);

		MarketDataDictionary*	vMktDataDico = new MarketDataDictionary(vAsOfDateARM, vMktDataMgr);

		vector<CCString>::iterator	vIdIter = C_MarketDataIds.begin();
		while( vIdIter != C_MarketDataIds.end() )
		{
			long	vMarketDataId = LocalGetNumObjectId(*vIdIter);
			ARM_Object*	vMarketData = LOCAL_PERSISTENT_OBJECTS->GetARMObject(vMarketDataId)->Clone();

			vMetaMarketData = EncapsulateMarketData(vMarketData);

			vMktDataDico->InsertMarketData(vMetaMarketData->GetKey(), vMetaMarketData); 

			vIdIter++;
		}

		(vMktDataMgr->GetMktDictionaries())[vAsOfDate] = vMktDataDico;
		vMktDataMgr->SetInitDone();


		if(objId == -1)
		{
			vNewMMDM = new ARM_MercureMarketDataManager(vMktDataMgr, vAsOfDate);

			if(vNewMMDM == NULL)
			{
				result.setMsg ("ARM_ERR: Mercure MarketDataManager is null");
				return ARM_KO;
			}

			CREATE_GLOBAL_OBJECT();
			
			vMDMId = LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)vNewMMDM);

			if(vMDMId == RET_KO)
			{
				delete	vNewMMDM;
				vNewMMDM = NULL;

				result.setMsg ("ARM_ERR: Pb inserting object");				
				return ARM_KO;
			}

			// Garder une copie originale pour pouvoir comparer si les MktData ont changé
			map<string, ARM_Object*>&	vNewMarketDataList = vNewMMDM->GetMarketDataList();
			map<string, MetaMarketData*>::iterator	vMktDataIter4 = vMktDataDico->GetDictionary().begin();
			for(; vMktDataIter4 != vMktDataDico->GetDictionary().end(); vMktDataIter4++)
			{
				vNewMarketDataList[vMktDataIter4->first] = &(vMktDataIter4->second->GetMarketData());
			}


			result.setLong(vMDMId);

			return ARM_OK;
		}
		else
		{
			vMMDM = (ARM_MercureMarketDataManager*) LOCAL_PERSISTENT_OBJECTS->GetARMObject(objId);

			if(LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(vMMDM, ARM_MERCURE_MARKETDATAMANAGER) == 1)
			{
				vNewMMDM = new ARM_MercureMarketDataManager(vMktDataMgr, vAsOfDate);

				if(vNewMMDM == NULL)
				{
					result.setMsg ("ARM_ERR: Mercure MarketDataManager is null");
					return ARM_KO;
				}

				if (vMMDM)
				{
					map<string, ARM_Object*>&	vPreviousMarketDataList = vMMDM->GetMarketDataList();
					map<string, ARM_Object*>&	vNewMarketDataList = vNewMMDM->GetMarketDataList();

					MarketDataManager*	vPreviousMDM = vMMDM->GetMarketDataManager();
					map<string, MarketDataDictionary*>::iterator	vPreviousDicoIter = vPreviousMDM->GetMktDictionaries().find(vAsOfDate);
					if(vPreviousDicoIter != vPreviousMDM->GetMktDictionaries().end())
					{
						MarketDataDictionary*	vPreviousDico = vPreviousDicoIter->second;
						map<string, MetaMarketData*>::iterator	vMktDataIter = vPreviousDico->GetDictionary().begin();
						map<string, MetaMarketData*>::iterator	vMktDataEnd = vPreviousDico->GetDictionary().end();
						for(; vMktDataIter != vMktDataEnd; vMktDataIter++)
						{
							MetaMarketData*	vMMD = vMktDataIter->second;
							ARM_Object*		vARMMD = &(vMMD->GetMarketData());

							map<string, MetaMarketData*>::iterator	vMktDataIter2 = vMktDataDico->GetDictionary().begin();
							for(; vMktDataIter2 != vMktDataDico->GetDictionary().end(); vMktDataIter2++)
							{
								// Ne pas détruire les objets communs aux 2 MarketDataManager
								if( vARMMD == &(vMktDataIter2->second->GetMarketData()) )
								{
									vMktDataIter->second = NULL;
									break;
								}
							}

							// Ne pas détruire les MarketData qui ont changé.
							// Ils ont déjà été détruits lors du rafraîchissement
							if( (vMktDataIter->second != NULL) &&
								( vNewMarketDataList[vMktDataIter->first] != &(vMktDataIter->second->GetMarketData()) )
							  )
							{
								vMktDataIter->second = NULL;
							}
						}
					}

					// Garder une copie originale pour pouvoir comparer si les MktData ont changé
					map<string, MetaMarketData*>::iterator	vMktDataIter3 = vMktDataDico->GetDictionary().begin();
					for(; vMktDataIter3 != vMktDataDico->GetDictionary().end(); vMktDataIter3++)
					{
						vNewMarketDataList[vMktDataIter3->first] = &(vMktDataIter3->second->GetMarketData());
					}

					delete	vMMDM;
					vMMDM = NULL;
				}

				LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_MercureMarketDataManager*)vNewMMDM, objId);

				return ARM_OK;
			}

			else
			{
				result.setMsg ("ARM_ERR: previous object is not of a good type");
				return	ARM_KO;
			}
		}
    }

    catch(Exception& x)
    {
        x.DebugPrint();
 
		delete	vNewMMDM;
		vNewMMDM = NULL;

		ARM_RESULT();
    } 

	catch(ParsingException& vException)
	{
		result.setMsg (vException.GetMessage().c_str());
		return	ARM_KO;
	}

	catch(...)
	{
		result.setMsg("Unknow exception caught in CreateMarketDataManager");
		return	ARM_KO;
	}
}


MetaMarketData*	EncapsulateMarketData(ARM_Object* aMarketData)
{
	MetaMarketData*	vMetaMarketData = NULL;

	ARM_AbstractMarketClass*	vMarketData = (ARM_AbstractMarketClass*)aMarketData;
	ARM_Date	vAsOfDate = vMarketData->GetAsOfDate();
	
	string	vType = vMarketData->GetStrType();		// "ZERO", "BASIS USD", "VOL SWOPT HISTO", ...
	string	vCurrency = vMarketData->GetStrCurrency();		// "EUR"
	string	vIndex = " " + vMarketData->GetExternalIndex();	// "LIBOR"
	string	vCurveId = vMarketData->GetExternalCrvId() == "NONE" ? "" : vMarketData->GetExternalCrvId();	
						// "MO", "MOSMILE", ...

	if(vType == "BASIS USD" || vMarketData->GetExternalIndex() == "")
		vIndex = "";

	char str[100];
	sprintf(str, "%s %s %s %s", vType.c_str(), vCurrency.c_str(), vIndex.c_str(), vCurveId.c_str());
	string	vKey = string(str);


	switch(vMarketData->GetName())
	{
		case ARM_ZERO_LIN_INTERPOL:
		{
			ARM_ZeroLInterpol*	vCurve = (ARM_ZeroLInterpol*)vMarketData;

			string	vSmoothingMethod = vCurve->GetSmoothingMethod().empty() ? "" : vCurve->GetSmoothingMethod();

			sprintf(str, "%s %s %s %s", vType.c_str(), vCurrency.c_str(), vIndex.c_str(), vSmoothingMethod.c_str(), vCurveId.c_str());
			vKey = string(str);

			vMetaMarketData = new MetaZeroCurve(vAsOfDate, vKey, vCurve, vSmoothingMethod);
		}

		break;


		case ARM_VOL_LIN_INTERPOL:
		{
			ARM_VolLInterpol*	vCurve = (ARM_VolLInterpol*)vMarketData;

			if(vCurve->GetStrType() == "VOL LIBOR CMS")
			{
				sprintf(str, "%s %s CORR %s", vType.c_str(), vCurrency.c_str(), vCurveId.c_str());
				vKey = string(str);
			}
			
			if(vCurve->GetExternalIndex() == "ROLIB" || vCurve->GetExternalIndex() == "NULIB")	// SABR
			{
				sprintf(str, "%s %s %s %s", vType.c_str(), vCurrency.c_str(), vIndex.c_str(), vCurve->GetExternalIndex().c_str(), vCurveId.c_str());
				vKey = string(str);
			}

			vMetaMarketData = new MetaVolatilityCurve(vAsOfDate, vKey, vCurve);
		}

		break;


		case ARM_VOL_CUBE:
		{
			ARM_VolCube*	vCurve = (ARM_VolCube*)vMarketData;

			vMetaMarketData = new MetaVolatilityCube(vAsOfDate, vKey, vCurve);
		}
			
		break;


		case ARM_SCALAR_DATA:
		{
			ARM_ScalarData*	vCurve = (ARM_ScalarData*)vMarketData;

			vMetaMarketData = new MetaScalarData(vAsOfDate, vKey, vCurve);
		}

		break;


		case ARM_VOL_CURVE:
		{
			ARM_VolCurve*	vCurve = (ARM_VolCurve*)vMarketData;

			vMetaMarketData = new MetaFXVolatilityCurve(vCurve->GetAsOfDate(), vKey, vCurve);
		}

		break;
	}

	return	vMetaMarketData;
}


long	ARMCreateScalarData(double	 C_Value,
							CCString C_Type,
							CCString C_Currency,
							CCString C_Index,
							double   C_AsOfDate,
							CCString C_CurveId,
							ARM_result& result,
							long objId)
{
	long	ARMScalarDataId;

    ARM_ScalarData*	vARMScalarData = NULL;
	ARM_ScalarData*	vNewARMScalarData = NULL;
 
	if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0)
	{
		result.setMsg ("ARM_ERR: Pb accessing objects");
		return ARM_KO;
	}

	char	sAsOfDate[11];

	CCString msg ("");

    try
    {
		if( (C_Type == "FX") || (C_Type == "CUTOFF 2F") || (C_Type == "CUTOFF 3F") ||
			(C_Type == "MEANREV 2F") || (C_Type == "MEANREV 3F") || (C_Type == "QMOD0") || (C_Type == "QMOD1") ||
			(C_Type == "FXIR CORR 2F") || (C_Type == "FXIR CORR 3F") || 
			(C_Type == "IRIR CORR 2F") || (C_Type == "IRIR CORR 3F") )
		{
			Local_XLDATE2ARMDATE(C_AsOfDate, sAsOfDate);
			string		vAsOfDate(sAsOfDate);
			ARM_Date	vAsOfDateARM(sAsOfDate);

			string	vType(C_Type);
			string	vCurrency(C_Currency);
			string	vIndex(C_Index);
			string	vCurveId(C_CurveId);

			if(objId == -1)
			{
				// Création d'un nouveau ARM_ScalarData
				vNewARMScalarData = new ARM_ScalarData(C_Value);
				vNewARMScalarData->SetStrType(vType);
				vNewARMScalarData->SetStrCurrency(vCurrency);
				vNewARMScalarData->SetExternalIndex(vIndex);
				vNewARMScalarData->SetExternalCrvId(vCurveId);
				vNewARMScalarData->SetAsOfDate(vAsOfDateARM);

				if(vNewARMScalarData == NULL)
				{
					result.setMsg ("ARM_ERR: ARMScalarData is null");
					return ARM_KO;
				}

				CREATE_GLOBAL_OBJECT();
				
				ARMScalarDataId = LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)vNewARMScalarData);

				if (ARMScalarDataId == RET_KO)
				{
					delete	vNewARMScalarData;
					vNewARMScalarData = NULL;

					result.setMsg ("ARM_ERR: Pb inserting object");				
					return ARM_KO;
				}

				result.setLong(ARMScalarDataId);

				return ARM_OK;
			}
			else
			{
				vARMScalarData = (ARM_ScalarData*) LOCAL_PERSISTENT_OBJECTS->GetARMObject(objId);

				if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(vARMScalarData, ARM_SCALAR_DATA) == 1)
				{
					if (vARMScalarData)
					{
						delete	vARMScalarData;
						vARMScalarData = NULL;
					}

					// Création d'un nouveau ARM_ScalarData
					vNewARMScalarData = new ARM_ScalarData(C_Value);
					vNewARMScalarData->SetStrType(vType);
					vNewARMScalarData->SetStrCurrency(vCurrency);
					vNewARMScalarData->SetExternalIndex(vIndex);
					vNewARMScalarData->SetExternalCrvId(vCurveId);
					vNewARMScalarData->SetAsOfDate(vAsOfDateARM);

					if(vNewARMScalarData == NULL)
					{
						result.setMsg ("ARM_ERR: ARMScalarData is null");
						return ARM_KO;
					}

					LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_ScalarData*)vNewARMScalarData, objId);

					return ARM_OK;
				}

				else
				{
					result.setMsg ("ARM_ERR: previous object is not of a good type");
					return	ARM_KO;
				}
			}
		}
		else
		{
			result.setMsg ("ARM_ERR: Type should be one of the following strings : \
							FX, CUTOFF 2F, CUTOFF 3F, MEANREV 2F, MEANREV 3F, QMOD0, QMOD1, \
							FXIR CORR 2F, FXIR CORR 3F, IRIR CORR 2F, IRIR CORR 3F");
			return ARM_KO;
		}
    }

    catch(Exception& x)
    {
        x.DebugPrint();
 
		delete	vNewARMScalarData;
		vNewARMScalarData = NULL;

		ARM_RESULT();
    } 

	catch(ParsingException& vException)
	{
		result.setMsg (vException.GetMessage().c_str());
		return	ARM_KO;
	}

	catch(...)
	{
		result.setMsg("Unknow exception caught in ARMCreateScalarData");
		return	ARM_KO;
	}
}


long ARMGetPostProcessedData(long C_HedgeId, 
							 CCString C_DataName, 
							 CCString C_PlotName,
							 ARM_result& result)
{
	CCString msg ("");
	try
	{
		ARM_MercureResult*	vMercureResult = NULL;
		vMercureResult = (ARM_MercureResult*) LOCAL_PERSISTENT_OBJECTS->GetARMObject(C_HedgeId);

		if ( (strcmp(C_DataName.c_str(), "PV") != 0) && (C_PlotName == "") )
		{
			result.setMsg ("ARM_ERR: a PlotName must be specified in 1st order sensitivity");
			return ARM_KO;
		}

		if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(vMercureResult, ARM_MERCURE_RESULT) == 1)
		{
			double dataValue = vMercureResult->GetPostProcessedData(string(C_DataName), string(C_PlotName));
			result.setDouble(dataValue);
			return ARM_OK;
		}
		else
		{
			result.setMsg ("ARM_ERR: Mercure object is not of good type");
			return ARM_KO;
		}
	}
    catch(Exception& x)
    {
        x.DebugPrint();
 		ARM_RESULT();
    } 
}


long ARMViewModelParams(CCString MetaModelName, 
						ARM_result& result,
						long objId)
{
	CCString msg ("");
	long MercureHelpId;
	ARM_MercureHelp* mercureHelp = NULL;

	try
	{
		if ( (IsMetaModel(string(MetaModelName))) || MetaModelName == "ALL" )
		{
			if(objId == -1)
			{
				mercureHelp = new ARM_MercureHelp(string(MetaModelName));

				if(mercureHelp == NULL)
				{
					result.setMsg ("ARM_ERR: Mercure Help object is null");
					return ARM_KO;
				}

				CREATE_GLOBAL_OBJECT();
				
				MercureHelpId = LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)mercureHelp);

				if (MercureHelpId == RET_KO)
				{
					delete	mercureHelp;
					mercureHelp = NULL;

					result.setMsg ("ARM_ERR: Pb inserting object");				
					return ARM_KO;
				}

				result.setLong(MercureHelpId);

				return ARM_OK;
			}
			else
			{
				mercureHelp = (ARM_MercureHelp*) LOCAL_PERSISTENT_OBJECTS->GetARMObject(objId);

				if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(mercureHelp, ARM_MERCURE_HELP) == 1)
				{
					if (mercureHelp)
					{
						delete	mercureHelp;
						mercureHelp = NULL;
					}

					// Création d'un nouveau MercureHelp
					mercureHelp = new ARM_MercureHelp(string(MetaModelName));

					if(mercureHelp == NULL)
					{
						result.setMsg ("ARM_ERR: Mercure Help object is null");
						return ARM_KO;
					}

					LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_MercureHelp*)mercureHelp, objId);

					return ARM_OK;
				}
				else
				{
					result.setMsg ("ARM_ERR: previous object is not of a good type");
					return	ARM_KO;
				}
			}
		}
		else
		{
			char msg[200];
			sprintf(msg, "ARM_ERR: %s is not a MetaModel in Mercure. \
						  Expected : 2IRFXModel, BS, BSGen, BSSmiled, CRACalculator, \
						  Mixture, Multi3F, FXOption3F, QModel", MetaModelName);

			result.setMsg(msg);
			return ARM_KO;
		}
	}
    catch(Exception& x)
    {
        x.DebugPrint();
 		ARM_RESULT();
    } 
	/// catch the rest
	catch (...)
	{
		result.setMsg ("ARM_ERR: unrecognized failure");
		return ARM_KO;
	}

	return ARM_OK;
}
