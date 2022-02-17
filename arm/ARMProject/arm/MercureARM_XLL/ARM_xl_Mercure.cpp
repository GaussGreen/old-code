#pragma warning(disable : 4786)
#pragma warning(disable : 4005)

#ifdef __STDC__
#undef __STDC__
#endif

#include "ARM_xl_trycatch_local.h"
#include <util\fromto.h>

#include "ARM_local_interglob.h"
#include "ARM_local_interface.h"

#include <ARM\libarm_local\ARM_local_class.h>

#include "ARM_local_Mercure.h"

#include "ARM_xl_Mercure.h"

#include <string>
#include "BasicParser.h"
#include "ARM_local_etoolkit.h"
#include "ARM_local_etoolkitX.h"


using namespace mercure;



__declspec(dllexport) LPXLOPER WINAPI	Mercure_Hedge (	LPXLOPER XL_TradeId,
														LPXLOPER XL_TradeType,
														LPXLOPER XL_AsOfDate,
														LPXLOPER XL_FallBack,
														LPXLOPER XL_MarketDataType,
														LPXLOPER XL_MarketDataManagerId,
														LPXLOPER XL_HedgeRatioFilesDirectory,
														LPXLOPER XL_ConfigFile,
														LPXLOPER XL_MarketDataFile,
														LPXLOPER XL_TraceSecurityAndModel,
														LPXLOPER XL_HedgeRatioFile )
{
	static XLOPER	XL_result;
	ARM_result		C_result;

	ARM_XL_TRY_BLOCK_BEGIN
	{
		ARM_NOCALCIFWIZ();

		CCString C_TradeId;
		vector<CCString>	vC_TradeId;
		CCString	C_TradeType;
		double		C_AsOfDate;
		vector<CCString>	C_FallBack;
		CCString	C_MarketDataType;
		CCString	C_MarketDataManagerId;
		CCString	C_HedgeRatioFilesDirectory;
		CCString	C_ConfigFile;
		CCString	C_MarketDataFile;
		CCString	C_TraceSecurityAndModel;
		bool		traceSecurityAndModel;
		CCString	C_HedgeRatioFile;

		static int		error;
		static char*	reason = "";

		XL_readStrCell(XL_TradeId, C_TradeId," ARM_ERR: TradeId: string expected", C_result);
		XL_readStrCellWD(XL_TradeType, C_TradeType, "ARM_Security"," ARM_ERR: TradeType: string expected", C_result);
		XL_readNumCell(XL_AsOfDate, C_AsOfDate," ARM_ERR: AsOfDate: string expected", C_result);
		XL_readStrVector(XL_FallBack, C_FallBack," ARM_ERR: FallBack: string expected", XL_TYPE_STRING, C_result);
		XL_readStrCellWD(XL_MarketDataType, C_MarketDataType, "", " ARM_ERR: MarketDataType: string expected", C_result);
		XL_readStrCellWD(XL_MarketDataManagerId, C_MarketDataManagerId, "", " ARM_ERR: MarketDataManager: string expected", C_result);
		XL_readStrCellWD(XL_HedgeRatioFilesDirectory, C_HedgeRatioFilesDirectory, "C:\\Mercure\\Parameters_Weekly", " ARM_ERR: HedgeRatioFilesDirectory: string expected", C_result);
		XL_readStrCellWD(XL_ConfigFile, C_ConfigFile, "C:\\Mercure\\SERVER\\MercureServer\\Debug\\MercureConfigFile.xml", " ARM_ERR: ConfigFile: string expected", C_result);
		XL_readStrCellWD(XL_MarketDataFile, C_MarketDataFile, "C:\\Mercure\\SERVER\\MercureServer\\Debug\\MercureInitialMarketDataFile.xml", " ARM_ERR: MarketDataFile: string expected", C_result);
		XL_readStrCellWD(XL_TraceSecurityAndModel, C_TraceSecurityAndModel, "NO", " ARM_ERR: Trace security and model: string expected", C_result);
		XL_readStrCellWD(XL_HedgeRatioFile, C_HedgeRatioFile, "", " ARM_ERR: HedgeRatioFile: string expected", C_result);

		C_TraceSecurityAndModel.toUpper();
		traceSecurityAndModel = (C_TraceSecurityAndModel == "YES") ? true : false;

		switchToETK();
		
/*		connection_etoolkit(SUMMIT_REC342_CONNEXION_USERNAME,
							SUMMIT_REC342_CONNEXION_PASSWD,
							SUMMIT_REC342_CONNEXION_CONTEXT,
							SUMMIT_REC342_IT_CONFIG_DOMAINSDIR,
							SUMMIT_REC342_IT_DOMAIN_NAME);
*/
		if(C_MarketDataType.GetLen() != 0)
		{
			C_HedgeRatioFilesDirectory += CCString("\\") + C_MarketDataType;
		}

		vC_TradeId.push_back(C_TradeId);

		long	retCode;
		long	objId;

		CCString	prevClass;	
		CCString	curClass = LOCAL_MERCURE_CLASS;
		CCString	stringId = GetLastCurCellEnvValue();
		
		if(!stringId)
		{
			retCode = ARMMercure_Hedge(	vC_TradeId,
										C_TradeType,
										C_AsOfDate,
										C_FallBack,
										C_MarketDataManagerId,
										C_HedgeRatioFilesDirectory,
										C_ConfigFile,
										C_MarketDataFile,
										traceSecurityAndModel,
										C_HedgeRatioFile,
										C_result );

			if(retCode == ARM_OK)
			{
				objId = C_result.getLong();
				LocalSetCurCellEnvValue(curClass, objId); 
				stringId = LocalMakeObjectId (objId, curClass);
			}
		}
		else
		{
			prevClass = LocalGetStringObjectClass(stringId);
			
			objId = LocalGetNumObjectId(stringId);
				
			if(curClass == prevClass)
			{
				retCode = ARMMercure_Hedge(	vC_TradeId,
											C_TradeType,
											C_AsOfDate,
											C_FallBack,
											C_MarketDataManagerId,
											C_HedgeRatioFilesDirectory,
											C_ConfigFile,
											C_MarketDataFile,
											traceSecurityAndModel,
											C_HedgeRatioFile,
											C_result,
											objId );

				if(retCode == ARM_OK)
				{			
					LocalSetCurCellEnvValue(curClass, objId); 
					stringId = LocalMakeObjectId (objId, curClass);
				}
			}
			else
			{
				FreeCurCellContent();
				retCode = ARMMercure_Hedge(	vC_TradeId,
											C_TradeType,
											C_AsOfDate,
											C_FallBack,
											C_MarketDataManagerId,
											C_HedgeRatioFilesDirectory,
											C_ConfigFile,
											C_MarketDataFile,
											traceSecurityAndModel,
											C_HedgeRatioFile,
											C_result );
			
				if(retCode == ARM_OK)
				{
					objId = C_result.getLong();
					LocalSetCurCellEnvValue(curClass, objId); 
					stringId = LocalMakeObjectId(objId, curClass);
				}
			}
		}

		if(retCode == ARM_OK)
		{
			FreeCurCellErr();
			XL_result.xltype = xltypeStr;
			XL_result.val.str = XL_StrC2StrPascal(stringId);
			XL_result.xltype |= xlbitDLLFree;
		}
		else
		{
			ARM_ERR();
		}
	}
	ARM_XL_TRY_BLOCK_END

	ARM_XL_CATCH_ARM_EXPT

	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Mercure_Hedge" )

	return (LPXLOPER)&XL_result;
}

__declspec(dllexport) LPXLOPER WINAPI	Mercure_Hedge_Array(LPXLOPER XL_TradeId,
															LPXLOPER XL_TradeType,
															LPXLOPER XL_AsOfDate,
															LPXLOPER XL_FallBack,
															LPXLOPER XL_MarketDataType,
															LPXLOPER XL_MarketDataManagerId,
															LPXLOPER XL_HedgeRatioFilesDirectory,
															LPXLOPER XL_ConfigFile,
															LPXLOPER XL_MarketDataFile,
															LPXLOPER XL_TraceSecurityAndModel,
															LPXLOPER XL_HedgeRatioFile)
{
	static XLOPER	XL_result;
	LPXLOPER		pxArray;
	ARM_result		C_result;

	ARM_XL_TRY_BLOCK_BEGIN
	{
		ARM_NOCALCIFWIZ();

		vector<CCString>	C_TradeId;
		CCString	C_TradeType;
		double		C_AsOfDate;
		vector<CCString>	C_FallBack;
		CCString	C_MarketDataType;
		CCString	C_MarketDataManagerId;
		CCString	C_HedgeRatioFilesDirectory;
		CCString	C_ConfigFile;
		CCString	C_MarketDataFile;
		CCString	C_TraceSecurityAndModel;
		bool		traceSecurityAndModel;
		CCString	C_HedgeRatioFile;

		static int		error;
		static char*	reason = "";

		XL_readStrVector(XL_TradeId, C_TradeId," ARM_ERR: TradeId: string expected", XL_TYPE_STRING, C_result);
		XL_readStrCellWD(XL_TradeType, C_TradeType, "ARM_Security", " ARM_ERR: TradeType: string expected", C_result);
		XL_readNumCell(XL_AsOfDate, C_AsOfDate," ARM_ERR: AsOfDate: string expected", C_result);
		XL_readStrVector(XL_FallBack, C_FallBack," ARM_ERR: FallBack: string expected", XL_TYPE_STRING, C_result);
		XL_readStrCellWD(XL_MarketDataType, C_MarketDataType, "", " ARM_ERR: MarketDataType: string expected", C_result);
		XL_readStrCellWD(XL_MarketDataManagerId, C_MarketDataManagerId, "", " ARM_ERR: MarketDataManager: string expected", C_result);
		XL_readStrCellWD(XL_HedgeRatioFilesDirectory, C_HedgeRatioFilesDirectory, "C:\\Mercure\\Parameters_Weekly", " ARM_ERR: HedgeRatioFilesDirectory: string expected", C_result);
		XL_readStrCellWD(XL_ConfigFile, C_ConfigFile, "C:\\Mercure\\SERVER\\MercureServer\\Debug\\MercureConfigFile.xml", " ARM_ERR: ConfigFile: string expected", C_result);
		XL_readStrCellWD(XL_MarketDataFile, C_MarketDataFile, "C:\\Mercure\\SERVER\\MercureServer\\Debug\\MercureInitialMarketDataFile.xml", " ARM_ERR: MarketDataFile: string expected", C_result);
		XL_readStrCellWD(XL_TraceSecurityAndModel, C_TraceSecurityAndModel, "NO", " ARM_ERR: Trace security and model: string expected", C_result);
		XL_readStrCellWD(XL_HedgeRatioFile, C_HedgeRatioFile, "", " ARM_ERR: HedgeRatioFile: string expected", C_result);

		C_TraceSecurityAndModel.toUpper();
		traceSecurityAndModel = (C_TraceSecurityAndModel == "YES") ? true : false;

		switchToETK();
/*
		connection_etoolkit(SUMMIT_REC342_CONNEXION_USERNAME,
							SUMMIT_REC342_CONNEXION_PASSWD,
							SUMMIT_REC342_CONNEXION_CONTEXT,
							SUMMIT_REC342_IT_CONFIG_DOMAINSDIR,
							SUMMIT_REC342_IT_DOMAIN_NAME);
*/
		if(C_MarketDataType.GetLen() != 0)
		{
			C_HedgeRatioFilesDirectory += CCString("\\") + C_MarketDataType;
		}

		long	retCode = ARMMercure_Hedge(	C_TradeId,
											C_TradeType,
											C_AsOfDate,
											C_FallBack,
											C_MarketDataManagerId,
											C_HedgeRatioFilesDirectory,
											C_ConfigFile,
											C_MarketDataFile,
											traceSecurityAndModel,
											C_HedgeRatioFile,
											C_result,
											-2);	// Pour détruire l'ARMMercureResult

		string tradeType = CCSTringToSTLString(C_TradeType);

		bool isGlobal	= (tradeType.find("DTL")  != string::npos ? false : true);	 

		if(retCode == ARM_OK)
		{
			FreeCurCellErr ();

			vector<CCString>	vName;
			vector<double>		vValue;
			vector<CCString>	vValueCcy;
			vector<CCString>    vTradeId;

			if(ExtractHedgeResut(C_result, vName, vValue, vValueCcy, vTradeId, isGlobal) == false)
			{
				ARM_ARG_ERR();
				return (LPXLOPER)&XL_result;
			}

			int	vNbRow = vName.size();
			int	vNbCol = 4;
			
			XL_result.xltype = xltypeMulti;
			XL_result.val.array.columns = vNbCol;
			XL_result.val.array.rows = vNbRow; 
			XL_result.val.array.lparray = pxArray = (LPXLOPER)GlobalAlloc(	GMEM_ZEROINIT, vNbRow * vNbCol * sizeof (XLOPER) );

			for(int i = 0; i < vNbRow; i++)
			{
				pxArray[XL_Coordonnate2Rank (i, 0, vNbCol)].xltype = xltypeStr;
				pxArray[XL_Coordonnate2Rank (i, 0, vNbCol)].val.str = XL_StrC2StrPascal(vTradeId[i]);
				pxArray[XL_Coordonnate2Rank (i, 0, vNbCol)].xltype |= xlbitDLLFree;

				pxArray[XL_Coordonnate2Rank (i, 1, vNbCol)].xltype = xltypeStr;
				pxArray[XL_Coordonnate2Rank (i, 1, vNbCol)].val.str = XL_StrC2StrPascal(vName[i]);
				pxArray[XL_Coordonnate2Rank (i, 1, vNbCol)].xltype |= xlbitDLLFree;

				pxArray[XL_Coordonnate2Rank (i, 2, vNbCol)].xltype = xltypeNum;
				pxArray[XL_Coordonnate2Rank (i, 2, vNbCol)].val.num = vValue[i]; 
				pxArray[XL_Coordonnate2Rank (i, 2, vNbCol)].xltype |= xlbitDLLFree;

				pxArray[XL_Coordonnate2Rank (i, 3, vNbCol)].xltype = xltypeStr;
				pxArray[XL_Coordonnate2Rank (i, 3, vNbCol)].val.str = XL_StrC2StrPascal(vValueCcy[i]);
				pxArray[XL_Coordonnate2Rank (i, 3, vNbCol)].xltype |= xlbitDLLFree;
			}
		}
		else
		{
			ARM_ERR();
		}
	}
	ARM_XL_TRY_BLOCK_END

	ARM_XL_CATCH_ARM_EXPT

	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Mercure_Hedge" )

	return (LPXLOPER)&XL_result;
}


__declspec(dllexport) LPXLOPER WINAPI	Mercure_ARM_Hedge (	LPXLOPER XL_SecurityId,
															LPXLOPER XL_AsOfDate,
															LPXLOPER XL_FallBack,
															LPXLOPER XL_MarketDataType,
															LPXLOPER XL_MarketDataManagerId,
															LPXLOPER XL_HedgeRatioFilesDirectory,
															LPXLOPER XL_ConfigFile,
															LPXLOPER XL_MarketDataFile,
															LPXLOPER XL_TraceSecurityAndModel,
															LPXLOPER XL_HedgeRatioFile)
{
	static XLOPER	XL_result;
	ARM_result		C_result;

	ARM_XL_TRY_BLOCK_BEGIN
	{
		ARM_NOCALCIFWIZ();

		CCString	C_SecurityId;
		vector<CCString> vC_SecurityId;
		double		C_AsOfDate;
		vector<CCString>	C_FallBack;
		CCString	C_MarketDataType;
		CCString	C_MarketDataManagerId;
		CCString	C_HedgeRatioFilesDirectory;
		CCString	C_ConfigFile;
		CCString	C_MarketDataFile;
		CCString	C_TraceSecurityAndModel;
		bool		traceSecurityAndModel;
		CCString	C_HedgeRatioFile;

		static int		error;
		static char*	reason = "";

		XL_readStrCell(XL_SecurityId, C_SecurityId," ARM_ERR: SecurityId: string expected", C_result);
		XL_readNumCell(XL_AsOfDate, C_AsOfDate," ARM_ERR: AsOfDate: string expected", C_result);
		XL_readStrVector(XL_FallBack, C_FallBack," ARM_ERR: FallBack: string expected", XL_TYPE_STRING, C_result);
		XL_readStrCellWD(XL_MarketDataType, C_MarketDataType, "", " ARM_ERR: MarketDataType: string expected", C_result);
		XL_readStrCellWD(XL_MarketDataManagerId, C_MarketDataManagerId, "", " ARM_ERR: MarketDataManager: string expected", C_result);
		XL_readStrCellWD(XL_HedgeRatioFilesDirectory, C_HedgeRatioFilesDirectory, "C:\\Mercure\\Parameters_Weekly", " ARM_ERR: HedgeRatioFilesDirectory: string expected", C_result);
		XL_readStrCellWD(XL_ConfigFile, C_ConfigFile, "C:\\Mercure\\SERVER\\MercureServer\\Debug\\MercureConfigFile.xml", " ARM_ERR: ConfigFile: string expected", C_result);
		XL_readStrCellWD(XL_MarketDataFile, C_MarketDataFile, "C:\\Mercure\\SERVER\\MercureServer\\Debug\\MercureInitialMarketDataFile.xml", " ARM_ERR: MarketDataFile: string expected", C_result);
		XL_readStrCellWD(XL_TraceSecurityAndModel, C_TraceSecurityAndModel, "NO", " ARM_ERR: Trace security and model: string expected", C_result);
		XL_readStrCellWD(XL_HedgeRatioFile, C_HedgeRatioFile, "", " ARM_ERR: HedgeRatioFile: string expected", C_result);

		C_TraceSecurityAndModel.toUpper();
		traceSecurityAndModel = (C_TraceSecurityAndModel == "YES") ? true : false;

		switchToETK();
/*
		connection_etoolkit(SUMMIT_REC342_CONNEXION_USERNAME,
							SUMMIT_REC342_CONNEXION_PASSWD,
							SUMMIT_REC342_CONNEXION_CONTEXT,
							SUMMIT_REC342_IT_CONFIG_DOMAINSDIR,
							SUMMIT_REC342_IT_DOMAIN_NAME);
*/
		if(C_MarketDataType.GetLen() != 0)
		{
			C_HedgeRatioFilesDirectory += CCString("\\") + C_MarketDataType;
		}

		vC_SecurityId.push_back(C_SecurityId);

		long	retCode;
		long	objId;

		CCString	prevClass;	
		CCString	curClass = LOCAL_MERCURE_CLASS;
		CCString	stringId = GetLastCurCellEnvValue();
		
		if(!stringId)
		{
			retCode = ARMMercure_Hedge(	vC_SecurityId,
										"",
										C_AsOfDate,
										C_FallBack,
										C_MarketDataManagerId,
										C_HedgeRatioFilesDirectory,
										C_ConfigFile,
										C_MarketDataFile,
										traceSecurityAndModel,
										C_HedgeRatioFile,
										C_result );

			if(retCode == ARM_OK)
			{
				objId = C_result.getLong();
				LocalSetCurCellEnvValue(curClass, objId); 
				stringId = LocalMakeObjectId (objId, curClass);
			}
		}
		else
		{
			prevClass = LocalGetStringObjectClass(stringId);
			
			objId = LocalGetNumObjectId(stringId);
				
			if(curClass == prevClass)
			{
				retCode = ARMMercure_Hedge(	vC_SecurityId,
											"",
											C_AsOfDate,
											C_FallBack,
											C_MarketDataManagerId,
											C_HedgeRatioFilesDirectory,
											C_ConfigFile,
											C_MarketDataFile,
											traceSecurityAndModel,
											C_HedgeRatioFile,
											C_result,
											objId );

				if(retCode == ARM_OK)
				{			
					LocalSetCurCellEnvValue(curClass, objId); 
					stringId = LocalMakeObjectId (objId, curClass);
				}
			}
			else
			{
				FreeCurCellContent();
				retCode = ARMMercure_Hedge(	vC_SecurityId,
											"",
											C_AsOfDate,
											C_FallBack,
											C_MarketDataManagerId,
											C_HedgeRatioFilesDirectory,
											C_ConfigFile,
											C_MarketDataFile,
											traceSecurityAndModel,
											C_HedgeRatioFile,
											C_result );
			
				if(retCode == ARM_OK)
				{
					objId = C_result.getLong();
					LocalSetCurCellEnvValue(curClass, objId); 
					stringId = LocalMakeObjectId(objId, curClass);
				}
			}
		}

		if(retCode == ARM_OK)
		{
			FreeCurCellErr();
			XL_result.xltype = xltypeStr;
			XL_result.val.str = XL_StrC2StrPascal(stringId);
			XL_result.xltype |= xlbitDLLFree;
		}
		else
		{
			ARM_ERR();
		}
	}
	ARM_XL_TRY_BLOCK_END

	ARM_XL_CATCH_ARM_EXPT

	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Mercure_ARM_Hedge" )

	return (LPXLOPER)&XL_result;
}


bool ExtractHedgeResut(ARM_result& C_result, vector<CCString>& aName, 
					   vector<double>& aValue, vector<CCString>& aValueCcy, vector<CCString>& aTradeId,
					   bool isGlobal)
{
	if (!isGlobal)
	{
		return ExtractPortfolioHedgeResult(C_result, aName, aValue, aValueCcy, aTradeId);
	}

	vector<CCString> vXMLOutput = C_result.getStringVector();

	string codeError;
	string valueError;
	string tradeId;
	XMLFile	vXMLFile;

	vector<CCString>::iterator iter_vXMLOutput = vXMLOutput.begin();
	for( ; iter_vXMLOutput != vXMLOutput.end(); iter_vXMLOutput++)
	{
		vXMLFile.LoadXMLString(iter_vXMLOutput->c_str());

		IXMLDOMNode* vDealNode = vXMLFile.getNode("/Deal");
		tradeId = vXMLFile.getString(vDealNode, "TradeId");
		IXMLDOMNode* vHedgeRatioClassNode = vXMLFile.getNode("/Deal/Hedge_Ratio_Class");

		if(vHedgeRatioClassNode == NULL)
		{
			codeError = vXMLFile.getString(vDealNode, "ErrorType");
			valueError = vXMLFile.getString(vDealNode, "ErrorMessage");

			aName.push_back(codeError.c_str());
			aValue.push_back(0.0);
			aValueCcy.push_back(valueError.c_str());
			aTradeId.push_back(tradeId.c_str());
		}
		else if ((valueError = vXMLFile.getString(vHedgeRatioClassNode, "Error")) != "")
		{
			aName.push_back(CCString("ARM_ERR"));
			aValue.push_back(0.0);
			aValueCcy.push_back(valueError.c_str());
			aTradeId.push_back(tradeId.c_str());
		}
		else
		{
			long	vNbGPP = 0;
			IXMLDOMNodeList*	vGPPList = vXMLFile.getNodeList(vHedgeRatioClassNode, "GlobalPostProcessed");
			if(vGPPList == NULL)
			{
				aName.push_back(CCString("ARM_ERR"));
				aValue.push_back(0.0);
				aValueCcy.push_back(CCString("Error getting GlobalPostProcessed node list from XML result file"));
				aTradeId.push_back(tradeId.c_str());
			
			}
			else
			{
				vGPPList->get_length(&vNbGPP);
				for(int indexGPP = 0; indexGPP < vNbGPP; indexGPP++)
				{
					IXMLDOMNode*	vGPPNode = NULL;
					vGPPList->get_item(indexGPP, &vGPPNode);
					if(vGPPNode == NULL)
					{
						aName.push_back(CCString("ARM_ERR"));
						aValue.push_back(indexGPP);
						aValueCcy.push_back(CCString("Error getting GlobalPostProcessed node from XML result file"));
						aTradeId.push_back(tradeId.c_str());
					}
					else
					{
						string	codeMaturite = vXMLFile.getString(vGPPNode, "Name");
						double	sensi = vXMLFile.getDouble(vGPPNode, "Value");
						string	valueCcy = vXMLFile.getString(vGPPNode, "ValueCcy");

						aName.push_back(codeMaturite.c_str());
						aValue.push_back(sensi);
						aValueCcy.push_back(valueCcy.c_str());
						aTradeId.push_back(tradeId.c_str());
					}
					vGPPNode->Release();
				}
				vGPPList->Release();
			}	
			vHedgeRatioClassNode->Release();
		}
		vDealNode->Release();
	}
	return true;
}


bool ExtractPortfolioHedgeResult(ARM_result& C_result, vector<CCString>& aName, 
								 vector<double>& aValue, vector<CCString>& aValueCcy,
								 vector<CCString>& aTradeId)
{
	vector<CCString> vXMLOutput = C_result.getStringVector();

	string codeError;
	string valueError;
	string tradeId;
	XMLFile	vXMLFile;

	string rootNodeName("/Exotic_Portfolio");
	string mainNodeName("/Exotic_Portfolio/Deals");

	vector<CCString>::iterator iter_vXMLOutput = vXMLOutput.begin();
	for( ; iter_vXMLOutput != vXMLOutput.end(); iter_vXMLOutput++)
	{
		vXMLFile.LoadXMLString(iter_vXMLOutput->c_str());

		IXMLDOMNode* vPortfolioNode = vXMLFile.getNode(rootNodeName);

		long vNbDN = 0;
		IXMLDOMNode* vDealsNode = vXMLFile.getNode(mainNodeName);
		IXMLDOMNodeList* vDNList = vXMLFile.getNodeList(vDealsNode, "Deal");

		vDNList->get_length(&vNbDN);
		
		for (long indexDeal=0;indexDeal<vNbDN;indexDeal++)
		{
			IXMLDOMNode*	vDealNode = NULL;
			vDNList->get_item(indexDeal, &vDealNode);

			tradeId = vXMLFile.getString(vDealNode, "TradeId");

			IXMLDOMNodeList* vHRCList = vXMLFile.getNodeList(vDealNode, "Hedge_Ratio_Class");

			IXMLDOMNode* vHedgeRatioClassNode = NULL;
			vHRCList->get_item(0, &vHedgeRatioClassNode);
			
			if(vHedgeRatioClassNode == NULL)
			{
				codeError = vXMLFile.getString(vDealNode, "ErrorType");
				valueError = vXMLFile.getString(vDealNode, "ErrorMessage");

				aName.push_back(codeError.c_str());
				aValue.push_back(0.0);
				aValueCcy.push_back(valueError.c_str());
				aTradeId.push_back(tradeId.c_str());
			}
			else if ((valueError = vXMLFile.getString(vHedgeRatioClassNode, "Error")) != "")
			{
				aName.push_back(CCString("ARM_ERR"));
				aValue.push_back(0.0);
				aValueCcy.push_back(valueError.c_str());
				aTradeId.push_back(tradeId.c_str());
			}
			else
			{
				long	vNbGPP = 0;
				IXMLDOMNodeList*	vGPPList = vXMLFile.getNodeList(vHedgeRatioClassNode, "GlobalPostProcessed");
				if(vGPPList == NULL)
				{
					aName.push_back(CCString("ARM_ERR"));
					aValue.push_back(0.0);
					aValueCcy.push_back(CCString("Error getting GlobalPostProcessed node list from XML result file"));
					aTradeId.push_back(tradeId.c_str());
				
				}
				else
				{
					vGPPList->get_length(&vNbGPP);
					for(int indexGPP = 0; indexGPP < vNbGPP; indexGPP++)
					{
						IXMLDOMNode*	vGPPNode = NULL;
						vGPPList->get_item(indexGPP, &vGPPNode);
						if(vGPPNode == NULL)
						{
							aName.push_back(CCString("ARM_ERR"));
							aValue.push_back(indexGPP);
							aValueCcy.push_back(CCString("Error getting GlobalPostProcessed node from XML result file"));
							aTradeId.push_back(tradeId.c_str());
						}
						else
						{
							string	codeMaturite = vXMLFile.getString(vGPPNode, "Name");
							double	sensi = vXMLFile.getDouble(vGPPNode, "Value");
							string	valueCcy = vXMLFile.getString(vGPPNode, "ValueCcy");

							aName.push_back(codeMaturite.c_str());
							aValue.push_back(sensi);
							aValueCcy.push_back(valueCcy.c_str());
							aTradeId.push_back(tradeId.c_str());
						}
						vGPPNode->Release();
					}
					vGPPList->Release();
				}	
				vHedgeRatioClassNode->Release();
			}
			vDealNode->Release();
		}

	}
	return true;
}


__declspec(dllexport) LPXLOPER WINAPI	CreateMarketDataManager(LPXLOPER XL_MarketDataIds,
																LPXLOPER XL_AsOfDate,
																LPXLOPER XL_FallBack,
																LPXLOPER XL_SwitchToETK)
{
	static XLOPER	XL_result;
	ARM_result		C_result;

	ARM_XL_TRY_BLOCK_BEGIN
	{
		ARM_NOCALCIFWIZ();

		vector<CCString>	C_MarketDataIds;
		double				C_AsOfDate;
		vector<CCString>	C_FallBack;
		CCString			C_SwitchToETK;

		static int		error;
		static char*	reason = "";

		XL_readStrVector(XL_MarketDataIds, C_MarketDataIds,
						" ARM_ERR: MarketDataIds: string vector expected", XL_TYPE_STRING, C_result);
		XL_readNumCell(XL_AsOfDate, C_AsOfDate," ARM_ERR: AsOfDate: string expected", C_result);
		XL_readStrVector(XL_FallBack, C_FallBack," ARM_ERR: FallBack: string expected", XL_TYPE_STRING, C_result);
		XL_readStrCellWD(XL_SwitchToETK, C_SwitchToETK, "NO", " ARM_ERR: SwitchToETK: string expected", C_result);

		C_SwitchToETK.toUpper();
		if(C_SwitchToETK == "YES")
		{
			switchToETK();
		}
		else if(C_SwitchToETK == "RECETTE")
		{
			switchToETK();
		
			connection_etoolkit(SUMMIT_REC342_CONNEXION_USERNAME,
								SUMMIT_REC342_CONNEXION_PASSWD,
								SUMMIT_REC342_CONNEXION_CONTEXT,
								SUMMIT_REC342_IT_CONFIG_DOMAINSDIR,
								SUMMIT_REC342_IT_DOMAIN_NAME);
		}

		long	retCode;
		long	objId;

		CCString	prevClass;	
		CCString	curClass = LOCAL_MARKETDATAMANAGER_CLASS;
		CCString	stringId = GetLastCurCellEnvValue();
		
		if(!stringId)
		{
			retCode = ARMCreateMarketDataManager(C_MarketDataIds,
												 C_AsOfDate,
												 C_FallBack,
												 C_result);

			if(retCode == ARM_OK)
			{
				objId = C_result.getLong();
				LocalSetCurCellEnvValue(curClass, objId); 
				stringId = LocalMakeObjectId (objId, curClass);
			}
		}
		else
		{
			prevClass = LocalGetStringObjectClass(stringId);
			objId = LocalGetNumObjectId(stringId);
				
			if(curClass == prevClass)
			{
				retCode = ARMCreateMarketDataManager(C_MarketDataIds,
													 C_AsOfDate,
													 C_FallBack,
													 C_result,
													 objId);

				if(retCode == ARM_OK)
				{			
					LocalSetCurCellEnvValue(curClass, objId); 
					stringId = LocalMakeObjectId (objId, curClass);
				}
			}
			else
			{
				FreeCurCellContent();
				retCode = ARMCreateMarketDataManager(C_MarketDataIds,
													 C_AsOfDate,
													 C_FallBack,
													 C_result);
			
				if(retCode == ARM_OK)
				{
					objId = C_result.getLong();
					LocalSetCurCellEnvValue(curClass, objId); 
					stringId = LocalMakeObjectId(objId, curClass);
				}
			}
		}

		if(retCode == ARM_OK)
		{
			FreeCurCellErr();
			XL_result.xltype = xltypeStr;
			XL_result.val.str = XL_StrC2StrPascal(stringId);
			XL_result.xltype |= xlbitDLLFree;
		}
		else
		{
			ARM_ERR();
		}
	}

	ARM_XL_TRY_BLOCK_END

	ARM_XL_CATCH_ARM_EXPT

	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in CreateMarketDataManager" )

	return (LPXLOPER)&XL_result;
}


__declspec(dllexport) LPXLOPER WINAPI	CreateARMScalarData(LPXLOPER XL_Value,
															LPXLOPER XL_Type,
															LPXLOPER XL_Currency,
															LPXLOPER XL_Index,
															LPXLOPER XL_AsOfDate,
															LPXLOPER XL_CurveId,
															LPXLOPER XL_SwitchToETK)
{
	static XLOPER	XL_result;
	ARM_result		C_result;

	ARM_XL_TRY_BLOCK_BEGIN
	{
		ARM_NOCALCIFWIZ();

		double		C_Value;
		CCString	C_Type;
		CCString	C_Currency;
		CCString	C_Index;
		double		C_AsOfDate;
		CCString	C_CurveId;
		CCString	C_SwitchToETK;

		static int		error;
		static char*	reason = "";

		XL_readNumCell(XL_Value, C_Value," ARM_ERR: Value: string expected", C_result);
		XL_readStrCell(XL_Type, C_Type," ARM_ERR: Type: string expected", C_result);
		XL_readStrCell(XL_Currency, C_Currency," ARM_ERR: Type: string expected", C_result);
		XL_readStrCell(XL_Index, C_Index," ARM_ERR: Type: string expected", C_result);	
		XL_readNumCell(XL_AsOfDate, C_AsOfDate," ARM_ERR: AsOfDate: string expected", C_result);
		XL_readStrCellWD(XL_CurveId, C_CurveId, "MO", " ARM_ERR: CurveId: string expected", C_result);
		XL_readStrCellWD(XL_SwitchToETK, C_SwitchToETK, "NO", " ARM_ERR: SwitchToETK: string expected", C_result);

		C_SwitchToETK.toUpper();
		if(C_SwitchToETK == "YES")
		{
			switchToETK();
		}
		else if(C_SwitchToETK == "RECETTE")
		{
			switchToETK();
		
			connection_etoolkit(SUMMIT_REC342_CONNEXION_USERNAME,
								SUMMIT_REC342_CONNEXION_PASSWD,
								SUMMIT_REC342_CONNEXION_CONTEXT,
								SUMMIT_REC342_IT_CONFIG_DOMAINSDIR,
								SUMMIT_REC342_IT_DOMAIN_NAME);
		}

		long	retCode;
		long	objId;

		CCString	prevClass;	
		CCString	curClass = LOCAL_ARMSCALARDATA_CLASS;
		CCString	stringId = GetLastCurCellEnvValue();
		
		if(!stringId)
		{
			retCode = ARMCreateScalarData(C_Value,
										  C_Type,
										  C_Currency,
										  C_Index,
										  C_AsOfDate,
										  C_CurveId,
										  C_result);

			if(retCode == ARM_OK)
			{
				objId = C_result.getLong();
				LocalSetCurCellEnvValue(curClass, objId); 
				stringId = LocalMakeObjectId (objId, curClass);
			}
		}
		else
		{
			prevClass = LocalGetStringObjectClass(stringId);
			
			objId = LocalGetNumObjectId(stringId);
				
			if(curClass == prevClass)
			{
				retCode = ARMCreateScalarData(C_Value,
											  C_Type,
											  C_Currency,
											  C_Index,
											  C_AsOfDate,
											  C_CurveId,
											  C_result,
											  objId);

				if(retCode == ARM_OK)
				{			
					LocalSetCurCellEnvValue(curClass, objId); 
					stringId = LocalMakeObjectId (objId, curClass);
				}
			}
			else
			{
				FreeCurCellContent();
				retCode = ARMCreateScalarData(C_Value,
											  C_Type,
											  C_Currency,
											  C_Index,
											  C_AsOfDate,
											  C_CurveId,
											  C_result);			
				if(retCode == ARM_OK)
				{
					objId = C_result.getLong();
					LocalSetCurCellEnvValue(curClass, objId); 
					stringId = LocalMakeObjectId(objId, curClass);
				}
			}
		}

		if(retCode == ARM_OK)
		{
			FreeCurCellErr();
			XL_result.xltype = xltypeStr;
			XL_result.val.str = XL_StrC2StrPascal(stringId);
			XL_result.xltype |= xlbitDLLFree;
		}
		else
		{
			ARM_ERR();
		}
	}

	ARM_XL_TRY_BLOCK_END

	ARM_XL_CATCH_ARM_EXPT

	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in CreateARMScalarData" )

	return (LPXLOPER)&XL_result;
}


// Affichage d'une valeur calculée dans le postprocessing
// Sortie : nombre ou erreur
__declspec(dllexport) LPXLOPER WINAPI	Mercure_GetPostProcessedData(LPXLOPER XL_HedgeObject,
																	 LPXLOPER XL_DataName,
																	 LPXLOPER XL_PlotName)
{
	static XLOPER	XL_result;
	ARM_result		C_result;

	ARM_XL_TRY_BLOCK_BEGIN
	{
		CCString C_HedgeObject;
		CCString C_DataName, C_PlotName;

		static int		error;
		static char*	reason = "";

		XL_readStrCell(XL_HedgeObject, C_HedgeObject," ARM_ERR: Hedge object expected", C_result);
		XL_readStrCellWD(XL_DataName, C_DataName, "PV", " ARM_ERR: DataName : string expected", C_result);
		XL_readStrCellWD(XL_PlotName, C_PlotName, "", " ARM_ERR: DataName : string expected", C_result);
		C_DataName.toUpper();
		C_PlotName.toUpper();

		long hedgeId = LocalGetNumObjectId (C_HedgeObject);

		long retCode = ARMGetPostProcessedData(hedgeId, C_DataName, C_PlotName, C_result);

		if(retCode == ARM_OK)
		{
			FreeCurCellErr ();
			XL_result.xltype = xltypeNum;
			XL_result.val.num = C_result.getDouble();
		}
		else
		{
			ARM_ERR();
		}
	}

	ARM_XL_TRY_BLOCK_END

	ARM_XL_CATCH_ARM_EXPT

	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Mercure_GetPostProcessedData" )

	return (LPXLOPER)&XL_result;
}


// Affichage des paramètres d'un MetaModel donné
// renvoie un objet
__declspec(dllexport) LPXLOPER WINAPI	Mercure_ViewModelParams(LPXLOPER XL_MetaModelName)
{
	static XLOPER	XL_result;
	ARM_result		C_result;

	ARM_XL_TRY_BLOCK_BEGIN
	{
		ARM_NOCALCIFWIZ();

		static int		error;
		static char*	reason = "";

		CCString C_MetaModelName;

		XL_readStrCellWD(XL_MetaModelName, C_MetaModelName, "ALL", " ARM_ERR: Model name : string expected", C_result);

		long	retCode;
		long	objId;

		CCString	prevClass;	
		CCString	curClass = LOCAL_MERCURE_HELP_CLASS;
		CCString	stringId = GetLastCurCellEnvValue();

		if(!stringId)
		{
			retCode = ARMViewModelParams(C_MetaModelName, C_result);

			if (retCode == ARM_OK)
				objId = C_result.getLong();
		}
		else
		{
			prevClass = LocalGetStringObjectClass(stringId);			
			objId = LocalGetNumObjectId(stringId);
				
			if(curClass == prevClass)
			{
				retCode = ARMViewModelParams(C_MetaModelName, C_result, objId);
			}
			else
			{
				FreeCurCellContent();
				retCode = ARMViewModelParams(C_MetaModelName, C_result );
			
				if(retCode == ARM_OK)
					objId = C_result.getLong();
			}
		}

		if(retCode == ARM_OK)
		{
			LocalSetCurCellEnvValue(curClass, objId); 
			stringId = LocalMakeObjectId (objId, curClass);

			FreeCurCellErr();
			XL_result.xltype = xltypeStr;
			XL_result.val.str = XL_StrC2StrPascal(stringId);
			XL_result.xltype |= xlbitDLLFree;
		}
		else
		{
			ARM_ERR();
		}
	}

	ARM_XL_TRY_BLOCK_END
	ARM_XL_CATCH_ARM_EXPT
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in Mercure_ViewModelParams" )

	return (LPXLOPER)&XL_result;
}