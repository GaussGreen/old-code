#ifndef ARM_LOCAL_MERCURE_H
#define ARM_LOCAL_MERCURE_H


#include "MetaMarketData.h"


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
					  long objId = -1);

long ARMCreateMarketDataManager(vector<CCString>& C_MarketDataIds,
								double   C_AsOfDate,
								vector<CCString>& C_FallBack,
								ARM_result& result,
								long objId = -1);

mercure::MetaMarketData* EncapsulateMarketData(ARM_Object* vMarketData);

long ARMCreateScalarData(double	 C_Value,
						 CCString C_Type,
						 CCString C_Currency,
						 CCString C_Index,
						 double   C_AsOfDate,
						 CCString C_CurveId,
						 ARM_result& result,
						 long objId = -1);

long ARMGetPostProcessedData(long C_HedgeId, 
							 CCString C_DataName, 
							 CCString C_PlotName,
							 ARM_result& result);

long ARMViewModelParams(CCString C_MetaModelName, 
						ARM_result& result,
						long objId = -1);

#endif /* ARM_LOCAL_MERCURE_H */