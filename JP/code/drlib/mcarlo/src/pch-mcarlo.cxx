#include "../credit/src/pch-credit.cxx"
/*
#include "../mcarlo/edginc/AssetData.hpp"
#include "../mcarlo/edginc/CCMCopulaModel.hpp"
#include "../mcarlo/edginc/CondensedVector.hpp"
#include "../mcarlo/edginc/CopulaModel.hpp"
#include "../mcarlo/edginc/Dependence.hpp"
#include "../mcarlo/edginc/DependenceGauss.hpp"
#include "../mcarlo/edginc/DependenceGaussTerm.hpp"
#include "../mcarlo/edginc/DependentCopulaModel.hpp"
#include "../mcarlo/edginc/FastTDateMap.hpp"
#include "../mcarlo/edginc/HitSample.hpp"
#include "../mcarlo/edginc/HyperTrigUtils.hpp"
#include "../mcarlo/edginc/IndependentCopulaModel.hpp"
#include "../mcarlo/edginc/IndexedPerfList.hpp"
#include "../mcarlo/edginc/IQMCAssetRNG.hpp"
#include "../mcarlo/edginc/IQMCDiffusibleAsset.hpp"
#include "../mcarlo/edginc/IQMCDiffusibleAssetBase.hpp"
#include "../mcarlo/edginc/IQMCHelperBoundedDiffusion.hpp"
#include "../mcarlo/edginc/IQMCHelperDateTimeCache.hpp"
#include "../mcarlo/edginc/IQMCHelperMaxDiffusionDate.hpp"
#include "../mcarlo/edginc/IQMCHelperTimeLogic.hpp"
#include "../mcarlo/edginc/IQMCPureJumps.hpp"
#include "../mcarlo/edginc/IQMCRNGManager.hpp"
#include "../mcarlo/edginc/SVQmcImplemented.hpp"
#include "../mcarlo/edginc/IQMCStateVariableBase.hpp"
#include "../mcarlo/edginc/IRGridPointCache.hpp"
#include "../mcarlo/edginc/LocalVolGrid.hpp"
#include "../mcarlo/edginc/LRGenerator.hpp"
#include "../mcarlo/edginc/MarketDataFetcherCDS.hpp"
#include "../mcarlo/edginc/MarketDataFetcherCIS.hpp"
#include "../mcarlo/edginc/MarketDataFetcherDemo.hpp"
#include "../mcarlo/edginc/MarketDataFetcherE3F.hpp"
#include "../mcarlo/edginc/MarketDataFetcherSRM.hpp"
#include "../mcarlo/edginc/MCCache.hpp"
#include "../mcarlo/edginc/MCPathBase.hpp"
#include "../mcarlo/edginc/MCPathConfig.hpp"
#include "../mcarlo/edginc/MCPathConfigCCM.hpp"
#include "../mcarlo/edginc/MCPathConfigParametric.hpp"
#include "../mcarlo/edginc/MCPathConfigSRM.hpp"
#include "../mcarlo/edginc/MCPathConfigSRMGen.hpp"
#include "../mcarlo/edginc/MCPathConfigSRMGenSV.hpp"
#include "../mcarlo/edginc/MCPrices.hpp"
#include "../mcarlo/edginc/MCPricing.hpp"
#include "../mcarlo/edginc/MCProduct.hpp"
#include "../mcarlo/edginc/MCProductClient.hpp"
#include "../mcarlo/edginc/MCProductEngineClient.hpp"
#include "../mcarlo/edginc/MCProductIQuicks.hpp"
#include "../mcarlo/edginc/MCRandom.hpp"
#include "../mcarlo/edginc/MCRandomLite.hpp"
#include "../mcarlo/edginc/MCReader.hpp"
#include "../mcarlo/edginc/MCWriter.hpp"
#include "../mcarlo/edginc/MonteCarlo.hpp"
#include "../mcarlo/edginc/MonteCarloLib.hpp"
#include "../mcarlo/edginc/PastPathGenerator.hpp"
#include "../mcarlo/edginc/QMCAssetRNG.hpp"
#include "../mcarlo/edginc/QMCBasisSpreadDiffuse.hpp"
#include "../mcarlo/edginc/QMCCreditCIDDiffuse.hpp"
#include "../mcarlo/edginc/QMCCreditCIDJumps.hpp"
#include "../mcarlo/edginc/QMCCreditComposite.hpp"
#include "../mcarlo/edginc/QMCCreditDiffuse.hpp"
#include "../mcarlo/edginc/QMCEnergyDiffuse.hpp"
#include "../mcarlo/edginc/QMCEquityDiffuse.hpp"
#include "../mcarlo/edginc/QMCFXBaseDiffuse.hpp"
#include "../mcarlo/edginc/QMCGenDiffusibleAsset.hpp"
#include "../mcarlo/edginc/QMCHelperBoundedDiffusion.hpp"
#include "../mcarlo/edginc/QMCHelperCachingTimeLogic.hpp"
#include "../mcarlo/edginc/QMCHelperDateTimeCache.hpp"
#include "../mcarlo/edginc/QMCHelperFastDateTimeCache.hpp"
#include "../mcarlo/edginc/QMCHelperVectorDiffusionDate.hpp"
#include "../mcarlo/edginc/QMCPureJumps.hpp"
#include "../mcarlo/edginc/QMCRatesDiffuse.hpp"
#include "../mcarlo/edginc/QMCRNGManager.hpp"
//!!!#include "../mcarlo/edginc/QuantoCDSAlgorithm.hpp"
#include "../mcarlo/edginc/RefLevel.hpp"
#include "../mcarlo/edginc/Reprice.hpp"
#include "../mcarlo/edginc/RFLCopulaModel.hpp"
#include "../mcarlo/edginc/RFLCopulaParameters.hpp"
#include "../mcarlo/edginc/SampleList.hpp"
#include "../mcarlo/edginc/SimpathI.hpp"
#include "../mcarlo/edginc/SparseMatrix.hpp"
#include "../mcarlo/edginc/SRMBasisSpreadHJMDiffuse.hpp"
#include "../mcarlo/edginc/SRMBasisSpreadHJMUtil.hpp"
#include "../mcarlo/edginc/SRMBlack.hpp"
#include "../mcarlo/edginc/SRMConstants.hpp"
#include "../mcarlo/edginc/SRMCorrelation.hpp"
#include "../mcarlo/edginc/SRMCreditCIRDiffuse.hpp"
#include "../mcarlo/edginc/SRMCreditCIRUtil.hpp"
#include "../mcarlo/edginc/SRMCreditDiffuse.hpp"
#include "../mcarlo/edginc/SRMCreditHJMDiffuse.hpp"
#include "../mcarlo/edginc/SRMCreditHJMUtil.hpp"
#include "../mcarlo/edginc/SRMCreditLiborDiffuse.hpp"
#include "../mcarlo/edginc/SRMCreditLiborUtil.hpp"
#include "../mcarlo/edginc/SRMCreditUtil.hpp"
#include "../mcarlo/edginc/SRMEnergyDiffuse.hpp"
#include "../mcarlo/edginc/SRMEnergyUtil.hpp"
#include "../mcarlo/edginc/SRMEquityDiffuse.hpp"
#include "../mcarlo/edginc/SRMEquityDiffuseEuler.hpp"
#include "../mcarlo/edginc/SRMEquityDiffuseMappingMethod.hpp"
#include "../mcarlo/edginc/SRMEquityDiffuseSecondOrder.hpp"
#include "../mcarlo/edginc/SRMEquityDiffuseStrictSecondOrder.hpp"
#include "../mcarlo/edginc/SRMEquityUtil.hpp"
#include "../mcarlo/edginc/SRMFXDiffuse.hpp"
#include "../mcarlo/edginc/SRMFXUtil.hpp"
#include "../mcarlo/edginc/SRMFXVolSimple.hpp"
#include "../mcarlo/edginc/SRMGenDiffusibleAsset.hpp"
#include "../mcarlo/edginc/SRMICE.hpp"
#include "../mcarlo/edginc/SRMRatesFactor.hpp"
#include "../mcarlo/edginc/SRMRatesHJMDiffuse.hpp"
#include "../mcarlo/edginc/SRMRatesHJMUtil.hpp"
#include "../mcarlo/edginc/SRMRatesLiborDiffuse.hpp"
#include "../mcarlo/edginc/SRMRatesLiborUtil.hpp"
#include "../mcarlo/edginc/SRMRatesUtil.hpp"
#include "../mcarlo/edginc/SRMSwap.hpp"
#include "../mcarlo/edginc/SRMSwaption.hpp"
#include "../mcarlo/edginc/SRMSwaptionPricer.hpp"
#include "../mcarlo/edginc/SRMUtil.hpp"
#include "../mcarlo/edginc/SubReprice.hpp"
#include "../mcarlo/edginc/SVGenAggregatedSurvivalDiscFactor.hpp"
#include "../mcarlo/edginc/SVGenBarrier.hpp"
#include "../mcarlo/edginc/SVGenDateOfDefault.hpp"
#include "../mcarlo/edginc/SVGenDemo.hpp"
#include "../mcarlo/edginc/SVGenDiscFactor.hpp"
#include "../mcarlo/edginc/SVGenEnergyFuturePrice.hpp"
#include "../mcarlo/edginc/SVGenExpectedBasisForward.hpp"
#include "../mcarlo/edginc/SVGenExpectedDiscFactor.hpp"
#include "../mcarlo/edginc/SVGenExpectedEnergyFuture.hpp"
#include "../mcarlo/edginc/SVGenExpectedSpot.hpp"
#include "../mcarlo/edginc/SVGenExpectedSurvivalDiscFactor.hpp"
#include "../mcarlo/edginc/SVGenIRFloatRate.hpp"
#include "../mcarlo/edginc/SVGenIRSwap.hpp"
#include "../mcarlo/edginc/SVGenSpot.hpp"
#include "../mcarlo/edginc/SVGenSurvivalDiscFactor.hpp"
#include "../mcarlo/edginc/TemplateIdx.hpp"
#include "../mcarlo/edginc/SVPath.hpp"
#include "../mcarlo/edginc/ISVBase.hpp"
#include "../mcarlo/edginc/SVAllInterfaces.hpp"
*/
