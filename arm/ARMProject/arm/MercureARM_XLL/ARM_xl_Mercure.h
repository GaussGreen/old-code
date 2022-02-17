#ifndef ARM_XL_MERCURE_H
#define ARM_XL_MERCURE_H

#include <libCCxll\CCxll.h>
#include <vector>

class ARM_result;
class CCString;



// Hedge � partir d'un tradeId
// Sortie : objet ARM_MercureResult
__declspec(dllexport) LPXLOPER WINAPI	Mercure_Hedge (	LPXLOPER XL_TradeId,
														LPXLOPER XL_TradeType,
														LPXLOPER XL_AsOfDate,
														LPXLOPER XL_FallBack,
														LPXLOPER XL_MarketDataType,
														LPXLOPER XL_HedgeRatioFilesDirectory,
														LPXLOPER XL_ConfigFile,
														LPXLOPER XL_MarketDataFile,
														LPXLOPER XL_TraceSecurityAndModel,
														LPXLOPER XL_HedgeRatioFile);


// Hedge � partir d'un tradeId
// Sortie : matrice Excel 
__declspec(dllexport) LPXLOPER WINAPI	Mercure_Hedge_Array(LPXLOPER XL_TradeId,
															LPXLOPER XL_TradeType,
															LPXLOPER XL_AsOfDate,
															LPXLOPER XL_FallBack,
															LPXLOPER XL_HedgeRatioFilesDirectory,
															LPXLOPER XL_ConfigFile,
															LPXLOPER XL_MarketDataFile,
															LPXLOPER XL_TraceSecurityAndModel,
															LPXLOPER XL_HedgeRatioFile);


// Extraction de la cha�nde XML de r�sultat sous forme de vecteurs

bool ExtractHedgeResut(ARM_result& C_result, std::vector<CCString>& aName,std::vector<double>& aValue,
					   std::vector<CCString>& aValueCcy, std::vector<CCString>& aTradeId,
					   bool isGlobal = false);

bool ExtractPortfolioHedgeResult(ARM_result& C_result, vector<CCString>& aName, 
								 vector<double>& aValue, vector<CCString>& aValueCcy,
								 vector<CCString>& aTradeId);


// Hedge � partir d'un ARM Object d�j� cr�� sur la spread sheet
// Sortie : objet ARM_MercureResult
__declspec(dllexport) LPXLOPER WINAPI	Mercure_ARM_Hedge (	LPXLOPER XL_SecurityId,
															LPXLOPER XL_AsOfDate,
															LPXLOPER XL_FallBack,
															LPXLOPER XL_MarketDataType,
															LPXLOPER XL_HedgeRatioFilesDirectory,
															LPXLOPER XL_ConfigFile,
															LPXLOPER XL_MarketDataFile,
															LPXLOPER XL_TraceSecurityAndModel,
															LPXLOPER XL_HedgeRatioFile);


// Cr�ation d'un MarketDataManager � partir d'Excel
// Sortie : MarketDataManager
__declspec(dllexport) LPXLOPER WINAPI	CreateMarketDataManager(LPXLOPER XL_MarketDataIds,
																LPXLOPER XL_AsOfDate,
																LPXLOPER XL_FallBack,
																LPXLOPER XL_SwitchToETK);

// Cr�ation d'un ARMScalarData � partir d'Excel
// Sortie : ARMScalarData
__declspec(dllexport) LPXLOPER WINAPI	CreateARMScalarData(LPXLOPER XL_Value,
															LPXLOPER XL_Type,
															LPXLOPER XL_Currency,
															LPXLOPER XL_Index,
															LPXLOPER XL_AsOfDate,
															LPXLOPER XL_CurveId,
															LPXLOPER XL_SwitchToETK);

// Affichage d'une valeur calcul�e dans le postprocessing
// Sortie : valeur (ou erreur)
__declspec(dllexport) LPXLOPER WINAPI	Mercure_GetPostProcessedData(LPXLOPER XL_HedgeObject,
																	 LPXLOPER XL_DataName,
																	 LPXLOPER XL_PlotName);

// Affichage des param�tres d'un MetaModel donn�
// renvoie un objet
__declspec(dllexport) LPXLOPER WINAPI	Mercure_ViewModelParams(LPXLOPER XL_MetaModelName);

#endif