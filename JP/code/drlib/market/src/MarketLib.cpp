//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : MarketLib.cpp
//
//   Description : Force linker to include all relevant files/symbols
//
//   Author      : Mark A Robson
//
//   Date        : 2 March 2001
//
//

#include "edginc/config.hpp"
#include "edginc/MarketLib.hpp"
#include "edginc/Actual360.hpp"
#include "edginc/Actual365F.hpp"
#include "edginc/FixedSettlement.hpp"
#include "edginc/FlatFXVol.hpp"
#include "edginc/LinearStrikeVolRequest.hpp"
#include "edginc/PseudoSimpleEquity.hpp"
#include "edginc/ProtEquity.hpp"
#include "edginc/StruckEquity.hpp"
#include "edginc/UnitXCBWithVol.hpp"
#include "edginc/PercXCBWithVol.hpp"
#include "edginc/XCB.hpp"
#include "edginc/Future.hpp"
#include "edginc/VolSpline.hpp"
#include "edginc/MixedSettlement.hpp"
#include "edginc/PhysicalAndCashSettlement.hpp"
#include "edginc/BadDayNone.hpp"
#include "edginc/BadDayFollowing.hpp"
#include "edginc/BadDayModified.hpp"
#include "edginc/BadDayPrevious.hpp"
#include "edginc/BootstrappedYieldCurve.hpp"
#include "edginc/CashSwapCurve.hpp"
#include "edginc/FourPlusI.hpp"
#include "edginc/ZeroCurve3.hpp"
#include "edginc/B30360.hpp"
#include "edginc/B30E360.hpp"
#include "edginc/Actual365.hpp"
#include "edginc/ActualActual.hpp"
#include "edginc/Business252.hpp"
#include "edginc/Actual365FJ.hpp"
#include "edginc/B30E360I.hpp"
#include "edginc/B30EP360.hpp"
#include "edginc/CashSettleDate.hpp"
#include "edginc/CashSettlePeriod.hpp"
#include "edginc/PhysicalSettlement.hpp"
#include "edginc/AtMaturity.hpp"
#include "edginc/ESWDividend.hpp"
#include "edginc/ESWCashFlow.hpp"
#include "edginc/Fund.hpp"
#include "edginc/IlliquidStock.hpp"
#include "edginc/MultiFactors.hpp"
#include "edginc/IRVol.hpp"
#include "edginc/CEVJ.hpp"
#include "edginc/CEVJProcessed.hpp"
#include "edginc/LocVolRequest.hpp"
#include "edginc/LiquiditySpreadCurve.hpp"
#include "edginc/ClosedFormCDSPSandFA.hpp"
#include "edginc/ResetSchedule.hpp"
#include "edginc/AssetDDE.hpp"
#include "edginc/StochasticYieldCurve.hpp"
#include "edginc/PreciousMetal.hpp"
#include "edginc/CommodityIndex.hpp"
#include "edginc/DeltaStrikeVolSurface.hpp"
#include "edginc/CDSParSpreadsLegalBasis.hpp"
#include "edginc/AdjustedCDSParSpreads.hpp"
#include "edginc/BespokeCreditIndexMap.hpp"
#include "edginc/CreditIndexPreferred.hpp"
#include "edginc/Inflation.hpp"
#include "edginc/SimpleCashFlowStream.hpp"
#include "edginc/ParCDS.hpp"
#include "edginc/AdhocCashFlow.hpp"
#include "edginc/FixedCashFlow.hpp"
#include "edginc/CreditCashFlow.hpp"
#include "edginc/BaseMetal.hpp"
#include "edginc/DeltaToStrike.hpp"
#include "edginc/EnergyUnderlyer.hpp"
#include "edginc/EnergyCurve.hpp"
#include "edginc/EnergyFuturesCurve.hpp"
#include "edginc/EnergyInstVolExplicit.hpp"
#include "edginc/EnergyInstVolExplicitRegular.hpp"
#include "edginc/EnergyInstVolExplicitTier2.hpp"
#include "edginc/EnergyImpliedVolSurface.hpp"
#include "edginc/EnergyInstVolCalibrated.hpp"
#include "edginc/EnergyInstVolRegular.hpp"
#include "edginc/EnergyInstVolSeasonal.hpp"
#include "edginc/EnergyVolCurve.hpp"
#include "edginc/CreditIndexMap.hpp"
#include "edginc/CDSVolATMMatrix.hpp"
#include "edginc/CDSVolCubeBSImpliedSmile.hpp"
#include "edginc/CDSVolCubeMultiQSmile.hpp"
#include "edginc/VolVSCurve.hpp"
#include "edginc/SRMEQVol.hpp"
#include "edginc/SRMFXVol.hpp"
#include "edginc/IDecretionCurve.hpp"
#include "edginc/PrepayCurve.hpp"
#include "edginc/MarketObservable.hpp"
#include "edginc/IndexSpec.hpp"
#include "edginc/IndexSpecEQ.hpp"
#include "edginc/IndexSpecIR.hpp"
#include "edginc/IndexSpecFX.hpp"
#include "edginc/PiecewiseLinearIncrementalMappingFunction.hpp"
#include "edginc/PiecewiseFlatIncrementalMappingFunction.hpp"
#include "edginc/ZCBrzFI.hpp"
#include "edginc/ToolkitDebug.hpp"
#include "edginc/MQQuasiIRVol.hpp"
#include "edginc/SCIDparameters.hpp"
#include "edginc/SCIDCalibparameters.hpp"
#include "edginc/MQQuasiIRVolCMS.hpp"
#include "edginc/MQCMSPricer.hpp"
#include "edginc/VolProcessedMQ.hpp"
#include "edginc/VolProcessedMQCMS.hpp"
#include "edginc/CmRflParameters.hpp"
#include "edginc/CmCcmParameters.hpp"
#include "edginc/CmOnlyParameters.hpp"
#include "edginc/CompositeCreditEngineParameters.hpp"
#include "edginc/CmCcmRflParameters.hpp"
#include "edginc/CmCcmBaseCorrelationParameters.hpp"
#include "edginc/RflOnlyParameters.hpp"
#include "edginc/BaseCorrelationOnlyParameters.hpp"
#include "edginc/CmBaseCorrelationParameters.hpp"
#include "edginc/CorrelationCommon.hpp"
#include "edginc/FactorCorrelation.hpp"
#include "edginc/ProxyVol.hpp"
#include "edginc/IREngineTree.hpp"
#include "edginc/IREngineMC.hpp"
#include "edginc/IRSmile2Q.hpp"
#include "edginc/IRSmileQuasi2Q.hpp"
#include "edginc/IRSmileMQ.hpp"
#include "edginc/IRModelVNFM.hpp"
#include "edginc/ExplicitDates.hpp"
#include "edginc/RegimeFactor.hpp"
#include "edginc/RateRegime.hpp"
#include "edginc/HazardRegime.hpp"
#include "edginc/VolRegime.hpp"
#include "edginc/MultiRegimeFactor.hpp"
#include "edginc/AbsCdoParameters.hpp"

DRLIB_BEGIN_NAMESPACE
// to avoid header file
extern bool VolNormalLogLinkIn();
extern bool VolNormalLogLegacyLinkIn();
extern bool VolPreferredLinkIn();
extern bool VolQuadExpLinkIn();
extern bool VolQuadLinkIn();
extern bool VolSVLinkIn();
extern bool VolSVJLinkIn();
extern bool VolSVJJLinkIn();
extern bool VolSplineLinkIn();
extern bool VolExpSmileSkewLinkIn();
extern bool VolGammaOULinkIn();
extern bool VolIGOULinkIn();
extern bool VolCGMYHestonLinkIn();
extern bool VolAJDSuperLinkIn();
extern bool VolNormalLogModLinkIn();
extern bool VolMertonLinkIn();
extern bool VolHyperTrigLoad();
extern bool VolLogLinearPlusLinkIn();
extern bool VolMertonLVLinkIn();
extern bool UntweakableYCLinkIn();
extern bool InstrumentAssetLoad();
extern bool VolSVCJLinkIn();
extern bool VolSVCJLVLinkIn();
extern bool VolStochGarfLinkIn();
extern bool PDFLogNormalLinkIn();
extern bool UntweakableCDSParSpreadLinkIn();
extern bool FlatCDSSpotVolLinkIn();
extern bool CDOContingentLegLinkIn();
extern bool CDOFullContingentLegLinkIn();
extern bool CIDParametersLinkIn();
extern bool VarSwapBasisLinkIn();
extern bool CDSIndexParSpreadsLinkIn();
extern bool CreditIndexLinkIn();
extern bool FutureVanillaPriceVolSurfaceLinkIn();
extern bool FutureAssetLoad();
extern bool ICreditIndexMapLinkIn();
extern bool SRMEQVolLoad();
extern bool SRMFXVolLoad();
extern bool IRVolPairLoad();
extern bool CorrSwapSamplingAdjLoad();
extern bool CorrSwapBasisAdjLoad();
extern bool IndexSpecIRLoad();
extern bool EnergyFuturesCurveLoad();
extern bool BootstrappedBasisIndexCurveLoad();
extern bool UntweakableBasisIndexCurveLinkIn();
extern bool BrazilCDILoad();
extern bool ESWAverageLoad();
extern bool DeltaToStrikeLoad();
extern bool SRMFXVolSpotLoad();
extern bool ExplicitDatesLoad();
extern bool ISDAConventionLoad();
extern bool VSCurveDeltaPointwiseLinkIn();
extern bool VSCurveDeltaParallelLinkIn();
extern bool VSCurveCrossGammaLinkIn();
extern bool VSCurveVVegaLinkIn();
extern bool VSCurveMRRTweakLinkIn();
extern bool MQQuasiIRVolLoad();
extern bool MQQuasiIRVolCMSLoad();
extern bool CreditIndexSpreadRhoPointwiseLinkIn();
extern bool CreditIndexSpreadRhoParallelLinkIn();
extern bool IModelConfigMapperLoad();
extern bool MQCMSPricerLoad();
extern bool ModelConfigMapperLoad();
extern bool SCIDparametersLoad();
extern bool SCIDCalibparametersLoad();
extern bool VolProcessedMQLoad();
extern bool FloatCashFlowLoad();
extern bool VolProcessedMQCMSLoad();
extern bool PortfolioNameLinkIn();
extern bool CDOPortfolioLinkIn();
extern bool CreditTrancheLossConfigLinkIn();
extern bool DateTimeLiteLinkIn();
extern bool ICreditFeeLegLoad();
extern bool CDSVolATMMatrixLoad();
extern bool CDSVolCubeBSImpliedSmileLoad();
extern bool CDSVolCubeMultiQSmileLoad();
extern bool ProxyVolLoad();
extern bool VolTypeSensitiveStrikesLoad();
extern bool IREngineTreeLoad();
extern bool IREngineMCLoad();
extern bool IRSmile2QLoad();
extern bool IRSmileQuasi2QLoad();
extern bool IRSmileMQLoad();
extern bool IRModelVNFMLoad();
extern bool RegimeFactorLinkIn();
extern bool RateRegimeLinkIn();
extern bool HazardRegimeLinkIn();
extern bool VolRegimeLinkIn();
extern bool MultiRegimeFactorLinkIn();
extern bool MultiAssetLinkIn();
extern bool NToDefaultLossConfigLinkIn();
extern bool DateAdjustmentLoad();
extern bool ScheduleGeneratorLoad();

void CMarketLib::linkInClasses(){
    /* The sole purpose of this method is to ensure that the linker
       includes all symbols out of the market directory. Many symbols
       are automatically linked because they are used by other classes
       which are already included.

       An example of symbols that could be dropped would be an entire class
       representing volatility. For example FlatVol might be dropped since
       the only references might be throught the abstract parent class. */


    /* there is no order dependency - this could be automated */

    // Please don't add any new ...::TYPE entries here ---
    // use FFFLoad() where FFF.cpp is the source filename.

    bool success =
        AssetDDE::TYPE &&
        Actual360::TYPE &&
        Actual365F::TYPE &&
        Equity::TYPE &&
        FixedSettlement::TYPE &&
        RollingSettlement::TYPE &&
        YieldCurve::TYPE &&
        FlatVol::TYPE &&
        FlatFXVol::TYPE &&
        LinearStrikeVolRequest::TYPE &&
        EquityBase::TYPE &&
        SimpleEquity::TYPE &&
        PseudoSimpleEquity::TYPE &&
        ProtEquity::TYPE &&
        StruckEquity::TYPE &&
        UnitXCBWithVol::TYPE &&
        PercXCBWithVol::TYPE &&
        XCB::TYPE &&
        Future::TYPE &&
        IVolatilityBS::TYPE &&
        IVolatilityDVF::TYPE &&
        CVolBase::TYPE &&
        CVolProcessed::TYPE &&
        CVolProcessedBS::TYPE &&
        CVolRequest::TYPE &&
        VolSurface::TYPE &&
        MixedSettlement::TYPE &&
        PhysicalAndCashSettlement::TYPE &&
        BadDayNone::TYPE &&
        BadDayFollowing::TYPE   &&
        BadDayModified::TYPE   &&
        BadDayPrevious::TYPE &&
        BootstrappedYieldCurve::TYPE  &&
        CashSwapCurve::TYPE  &&
        IZeroCurveFactory::TYPE &&
        FourPlusI::TYPE	 &&
        ZeroCurve3::TYPE	 &&
        ZeroCurve::TYPE  &&
        FourPlusIZeroCurve::TYPE  &&
        ZC3ZeroCurve::TYPE  &&
        Correlation::TYPE &&
        CorrelationCommon::TYPE &&
        FactorCorrelation::TYPE &&
        DividendList::TYPE &&
        B30360::TYPE  &&
        B30E360::TYPE &&
        Actual365::TYPE &&
        ActualActual::TYPE &&
        FXAsset::TYPE &&
        BorrowCurve::TYPE &&
        Business252::TYPE &&
        Actual365FJ::TYPE  &&
        B30E360I::TYPE &&
        B30EP360::TYPE &&
        CashSettleDate::TYPE &&
        CashSettlePeriod::TYPE &&
        PhysicalSettlement::TYPE &&
        AtMaturity::TYPE &&
        CValueDateCollector::TYPE &&
        PayStream::load() &&
        FloatRate::TYPE &&
        CFXRateCollector::TYPE &&
        IDecretionCurve::TYPE &&
        PrepayCurve::TYPE &&
        DecretionCurve::TYPE && IABCDSDecretion::TYPE &&
        NullDecretionCurve::TYPE &&
        // probably don't need these - should get pulled in by Equity Swap
        ESWAverageLoad() &&
        ESWEquity::TYPE &&
        ESWLibor::TYPE &&
        ESWDividend::TYPE &&
        ESWCashFlow::TYPE &&
        Fund::TYPE &&
        INextStrike::TYPE &&
        IlliquidStock::TYPE &&
        IRVol::TYPE &&
        IRVolSpot::TYPE &&
        IMultiFactors::TYPE &&
        ResetSchedule::TYPE &&
        // param vols
        VolNormalLogLinkIn() &&
        VolNormalLogModLinkIn() &&
        VolNormalLogLegacyLinkIn() &&
        VolPreferredLinkIn() &&
        VolQuadExpLinkIn() &&
        VolQuadLinkIn() &&
        VolSVLinkIn() &&
        VolSVJLinkIn() &&
        VolSVJJLinkIn() &&
        VolSplineLinkIn() &&
        VolExpSmileSkewLinkIn() &&
        LocVolRequest::TYPE &&
        CEVJ::TYPE &&
        CEVJProcessed::TYPE &&
        CVolProcessedDVF::TYPE &&
        ICanBeRisky::TYPE &&
        LiquiditySpreadCurve::TYPE &&
        ClosedFormCDSPSandFA::TYPE &&
        VolGammaOULinkIn() &&
        VolIGOULinkIn() &&
        VolCGMYHestonLinkIn() &&
        VolAJDSuperLinkIn() &&
        VolMertonLinkIn() &&
        VolHyperTrigLoad() &&
        VolLogLinearPlusLinkIn() &&
        VolMertonLVLinkIn() &&
        VolSVCJLinkIn() &&
        VolSVCJLVLinkIn() &&
        VolStochGarfLinkIn() &&
        InstrumentAssetLoad() && 
        IStochasticYieldCurve::TYPE &&
        CurrencyBasis::TYPE &&
        PreciousMetal::TYPE &&
        CommodityIndex::TYPE &&
        UntweakableYCLinkIn() &&
        UntweakableCDSParSpreadLinkIn() &&
        UntweakableBasisIndexCurveLinkIn() &&
        FlatCDSSpotVolLinkIn() &&
        PDFLogNormalLinkIn() &&
        SRMEQVolLoad() &&
        SRMFXVolLoad() &&
        CDOContingentLegLinkIn() &&
        CDOFullContingentLegLinkIn() &&
        CIDParametersLinkIn() &&
		EnergyFuturesCurveLoad() &&
        DeltaStrikeVolSurface::TYPE &&
        Inflation::TYPE &&
        AdjustedCDSParSpreads::TYPE &&
        SimpleCashFlowStream::TYPE &&
        ParCDS::TYPE &&
        AbstractCashFlow::TYPE &&
        AdhocCashFlow::TYPE &&
        FixedCashFlow::TYPE &&
        CreditCashFlow::TYPE &&
        CDSParSpreadsLegalBasis::TYPE &&
        BespokeCreditIndexMap::TYPE &&
        ICreditIndexMap::TYPE &&
        CreditIndexPreferred::TYPE &&
        Duration::TYPE &&
        BaseMetal::TYPE &&
        VarSwapBasisLinkIn() &&
        CDSIndexParSpreadsLinkIn() &&
        CreditIndexLinkIn() &&
        DeltaToStrikeLoad() &&
        FutureAssetLoad() &&
        FutureVanillaPriceVolSurfaceLinkIn() &&
        EnergyUnderlyer::TYPE &&
        EnergyCurve::TYPE &&
        EnergyFuturesCurve::TYPE &&
        EnergyImpliedVolSurface::TYPE &&
        EnergyInstVolBase::TYPE &&
		EnergyInstVolExplicit::TYPE &&
        EnergyInstVolExplicitRegular::TYPE &&
        EnergyInstVolExplicitTier2::TYPE &&
        EnergyInstVolCalibrated::TYPE &&
        EnergyInstVolRegular::TYPE &&
        EnergyInstVolSeasonal::TYPE &&
        EnergyVolCurve::TYPE &&
        ICreditIndexMapLinkIn() &&
        VolVSCurve::TYPE &&
        IRVolPairLoad() &&
        CDSVolATMMatrixLoad() &&
        CDSVolCubeBSImpliedSmileLoad() &&
        CDSVolCubeMultiQSmileLoad() &&
		CorrSwapSamplingAdjLoad() &&
		CorrSwapBasisAdjLoad() &&
        VolProcessedMQLoad() &&
        VolProcessedMQCMSLoad() &&
        MQQuasiIRVolLoad() &&
        MQQuasiIRVolCMSLoad() &&
        SRMFXVolSpotLoad() &&
		IMarketObservable::TYPE &&
        IndexSpec::TYPE &&
        IndexSpecEQ::TYPE &&
        IndexSpecIR::TYPE &&
        IndexSpecFX::TYPE &&
        BootstrappedBasisIndexCurveLoad() &&
        BrazilCDILoad() &&
        PiecewiseLinearMappingFunction::TYPE &&
        PiecewiseFlatMappingFunction::TYPE &&
        PiecewiseMappingFunction::TYPE &&
        IMappingFunction::TYPE &&
        PiecewiseLinearIncrementalMappingFunction::TYPE &&
        PiecewiseFlatIncrementalMappingFunction::TYPE &&
        PiecewiseIncrementalMappingFunction::TYPE &&
        ZCBrzFI::TYPE &&
        ExplicitDatesLoad() &&
        ISDAConventionLoad() &&
		VSCurveDeltaPointwiseLinkIn() &&
        VSCurveDeltaParallelLinkIn() &&
		VSCurveCrossGammaLinkIn() &&
        VSCurveVVegaLinkIn() &&
        VSCurveMRRTweakLinkIn() &&
        CreditIndexSpreadRhoParallelLinkIn() &&
        CreditIndexSpreadRhoPointwiseLinkIn() &&
        IModelConfigMapperLoad() &&
        ModelConfigMapperLoad() &&
        FloatCashFlowLoad() &&
        MQCMSPricerLoad() &&
		SCIDparametersLoad() &&
		SCIDCalibparametersLoad() &&
        CmRflParameters::TYPE &&
        CmCcmParameters::TYPE &&
        CmOnlyParameters::TYPE &&
        CompositeCreditEngineParameters::TYPE &&
        CmCcmRflParameters::TYPE &&
        CmCcmBaseCorrelationParameters::TYPE &&
        RflOnlyParameters::TYPE &&
        BaseCorrelationOnlyParameters::TYPE &&
        CmBaseCorrelationParameters::TYPE &&
        PortfolioNameLinkIn() &&
		CDOPortfolioLinkIn() &&
		CreditTrancheLossConfigLinkIn() &&
		DateTimeLiteLinkIn() &&
		ICreditFeeLegLoad() &&
        ProxyVolLoad() &&
        VolTypeSensitiveStrikesLoad() && 
		IREngineTreeLoad() &&
		IREngineMCLoad() &&
		IRSmile2QLoad() &&
		IRSmileQuasi2QLoad() &&
		IRSmileMQLoad() &&
		IRModelVNFMLoad() &&
		RegimeFactorLinkIn() &&
		RateRegimeLinkIn() &&
		HazardRegimeLinkIn() &&
		VolRegimeLinkIn() &&
		MultiRegimeFactorLinkIn() &&
		MultiAssetLinkIn() &&
        NToDefaultLossConfigLinkIn() &&
		DateAdjustmentLoad() &&
		ScheduleGeneratorLoad() &&
		AbsCdoParameters::TYPE &&
        true;

    if (!success){
        throw ModelException("CMarketLib::registerClasses",
                             "Registration error");
    }
}

DRLIB_END_NAMESPACE
