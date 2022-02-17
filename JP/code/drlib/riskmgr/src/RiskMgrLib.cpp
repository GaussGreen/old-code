
#include "edginc/config.hpp"
#include "edginc/RiskMgrLib.hpp"
#include "edginc/FXCrossGamma.hpp"
#include "edginc/CrossGamma.hpp"
#include "edginc/VegaMatrix.hpp"
#include "edginc/RiskMgrInterface.hpp"
#include "edginc/SensitiveStrikes.hpp"
#include "edginc/ScenarioInterface.hpp"
#include "edginc/SpotLevel.hpp"
#include "edginc/SpotPrice.hpp"
#include "edginc/ThetaFwdSpot.hpp"
#include "edginc/ThetaFwdRate.hpp"
#include "edginc/ThetaNIE.hpp"
#include "edginc/DeltaNextDay.hpp"
#include "edginc/BasketDelta.hpp"
#include "edginc/RollingTheta.hpp"
#include "edginc/ImpliedScalarShift.hpp"
#include "edginc/DeltaProxy.hpp"
#include "edginc/VegaProxyParallel.hpp"
#include "edginc/DeltaProxyNextDay.hpp"
#include "edginc/SpotPriceProxy.hpp"
#include "edginc/ImpliedVol.hpp"
#include "edginc/ImpliedYTM.hpp"
#include "edginc/ImpliedYTP.hpp"
#include "edginc/CreditSpreadRhoParallel.hpp"
#include "edginc/CreditSpreadRhoPointwise.hpp"
#include "edginc/CreditSpreadLevel.hpp"
#include "edginc/NakedBondRhoParallel.hpp"
#include "edginc/NakedBondRhoPointwise.hpp"
#include "edginc/OutputRequestOCBInt.hpp"
#include "edginc/ScaleOutputs.hpp"
#include "edginc/HandlePaymentEvents.hpp"
#include "edginc/LiquiditySpreadRhoParallel.hpp"
#include "edginc/DeltaToCredit.hpp"
#include "edginc/DeltaToCreditNextDay.hpp"
#include "edginc/AssetVegaParallel.hpp"
#include "edginc/AssetVegaPointwise.hpp"
#include "edginc/E2CModel.hpp"
#include "edginc/AdjCreditSpreadRhoParallel.hpp"
#include "edginc/AdjCreditSpreadRhoPointwise.hpp"
#include "edginc/OptionStrikeRhoParallel.hpp"
#include "edginc/OptionStrikeRhoPointwise.hpp"
#include "edginc/IntrinsicMTM.hpp"
#include "edginc/EquityVegaEquivalentParallel.hpp"
#include "edginc/CompositeModel.hpp"
#include "edginc/VegaParallel2Sided.hpp"
#include "edginc/DeltaDDE.hpp"
#include "edginc/DDeltaDVol.hpp"
#include "edginc/VegaWeightedPhi.hpp"
#include "edginc/DeltaTPlusNHack.hpp"
#include "edginc/CreditDefaultSensBase.hpp"
#include "edginc/ClosedForm.hpp"
#include "edginc/VHTSkewFactorsPointwise.hpp"
#include "edginc/VHTSkewCelerityFactorsPointwise.hpp"
#include "edginc/VHTConvexityFactorsPointwise.hpp"
#include "edginc/VHTConvexityCelerityFactorsPointwise.hpp"
#include "edginc/IBootstrapper.hpp"
#include "edginc/PrepayParallel.hpp"
#include "edginc/ParSpreadUpfrontParallel.hpp"
#include "edginc/ToolkitDebug.hpp"

DRLIB_BEGIN_NAMESPACE

extern bool CalibratorLinkIn();
extern bool ParamSensLinkIn();
extern bool BasketDeltaNextDayLinkIn();
extern bool ShortTermSqueezeParallelLinkIn();
extern bool LongTermSqueezeParallelLinkIn();
extern bool CorrParallelShiftLinkIn();
extern bool CorrelationShiftLinkIn();
extern bool ParSpreadRhoParallelLinkIn();
extern bool ParSpreadRhoParallelTwoSidedLinkIn();
extern bool ParSpreadRhoPointwiseLinkIn();
extern bool ParSpreadUpfrontPointwiseLinkIn();
extern bool PrepayParallelLinkIn();
extern bool ParSpreadUpfrontParallelLinkIn();
extern bool LegalBasisAdditiveParallelLinkIn();
extern bool LegalBasisAdditivePointwiseLinkIn();
extern bool LegalBasisMultiplierParallelLinkIn();
extern bool LegalBasisMultiplierPointwiseLinkIn();
extern bool LegalBasisMultiplierRecoveryLinkIn();
extern bool LegalBasisAdditiveRecoveryLinkIn();
extern bool CCMDBetaDSpreadLinkIn();
extern bool CCMDBasketBetaDSpreadLinkIn();
extern bool CompositeVegaParallelLinkIn();
extern bool CompositeVegaPointwiseLinkIn();
extern bool SpotVegaParallelLinkIn();
extern bool SpotVegaPointwiseLinkIn();
extern bool CurrencyBasisRhoParallelLinkIn();
extern bool CurrencyBasisRhoPointwiseLinkIn();
extern bool CRSpotVegaParallelLinkIn();
extern bool CurrencyBasisSpreadLevelLinkIn();
extern bool CRMeanReversionParallelLinkIn();
extern bool IRMeanReversionParallelLinkIn();
extern bool SqueezeParallelTweakLinkIn();
extern bool CompressionRatioTweakLinkIn();
extern bool StrikeMappingTweakLinkIn();
extern bool StrikeMappingOverrideTweakLinkIn();
extern bool RecoveryTweakLinkIn();
extern bool RecoveryShiftLinkIn();
extern bool RateShiftLinkIn();
extern bool YCWeightedAdditiveShiftLinkIn();
extern bool YCWeightedMultiplicativeShiftLinkIn();
extern bool ParSpreadWeightedAdditiveShiftLinkIn();
extern bool ParSpreadWeightedMultiplicativeShiftLinkIn();
extern bool BetaSkewMatrixTweakLinkIn();
extern bool IRVegaMatrixLinkIn();
extern bool MultiRiskMgrInterfaceLinkIn();
extern bool RecoverySetLinkIn();
extern bool BCStrikeMappingTweakBaseLinkIn();
extern bool BCStrikeMappingTweakLinkIn();
extern bool BCStrikeMappingSetLinkIn();
extern bool BCStrikeMappingOverrideTweakLinkIn();
extern bool RecoveryTweakWithShiftLinkIn();
extern bool CreditDeltaPointwiseWithShiftLinkIn();
extern bool EnergyDeltaPointwiseWithShiftLinkIn();
extern bool LiquiditySpreadRhoPointwiseLinkIn();
extern bool CreditIndexBasisRhoParallelLinkIn();
extern bool CreditIndexBasisRhoPointwiseLinkIn();
extern bool ExposureReporterLinkIn();
extern bool FlexibleSensitivityLinkIn();
extern bool PerEntryFlexiblePerturbationLinkIn();
extern bool ScalarFlexiblePerturbationLinkIn();
extern bool PerturbationTestLinkIn();
extern bool RhoParallelTwoSidedLinkIn();
extern bool RhoPointwiseTwoSidedLinkIn();
extern bool IRDeltaPointwiseLinkIn();
extern bool CorrSwapBasisAdjTweakLinkIn();
extern bool CorrSwapBasisAdjParallelLinkIn();
extern bool CorrSwapSamplingAdjTweakLinkIn();
extern bool HolidayLoad();
extern bool VegaMatrixLoad();
extern bool QuasiContractualBetaSkewParallelLinkIn();
extern bool QuasiContractualBetaSkewMatrixLinkIn();
extern bool CRVegaParallelLinkIn();
extern bool CRVegaPointwiseLinkIn();
extern bool ModelFilterLinkIn();
extern bool VegaAtmParallelConstConvxLinkIn();
extern bool VegaAtmPointwiseConstConvxLinkIn();
extern bool MarketDataNameSubstitutionLinkIn();
extern bool CreditMultiBetaLevelLinkIn();
extern bool ImpliedScalarShiftMultiLinkIn();
extern bool CreditNameNotionalLevelLinkIn();
extern bool WeightedInstrumentCollectionLinkIn();
extern bool MultiModelLinkIn();
extern bool WeightedInstrumentShiftLinkIn();
extern bool CDOParallelStrikeShiftLinkIn();
extern bool CreditFeeLegCouponShiftLinkIn();
extern bool Smile2QElementwiseLinkIn();
extern bool CompositeForwardRatePricerLoad();
extern bool ClosedFormForwardRatePricerLoad();
//extern bool CreditIndexSpreadRhoPointwiseLinkIn();
extern bool LegalBasisAdditiveRelativeParallelLinkIn();
extern bool CorridorVarSwapSensLinkIn();
extern bool TierLoad();
extern bool ParSpreadRhoPointwiseLoad();
extern bool SPIGapRiskSensLoad();
       bool MOQTestMarketObjectLoad();
       bool MOQTestModelLoad();
       bool PhiLinkIn();
       bool PhiParallelLinkIn();

extern bool MarketObjectQualifiersLinkIn();

void CRiskMgrLib::linkInClasses(){
    /* The sole purpose of this method is to ensure that the linker
       includes all symbols out of the riskmr directory. Many symbols
       are automatically linked because they are used by other classes
       which are already included.

       An example of symbols that could be dropped would be an entire class
       representing a product which was referenced by no other classes */

    /* Since the classes in the riskmgr are naturally referenced by many other
       classes, the list of symbols to include is very short */

    /* there is no order dependency - this could be automated */

    // Please don't add any new ...::TYPE entries here ---
    // use FFFLoad() where FFF.cpp is the source filename.

    bool success =
        RiskMgrInterface::TYPE    &&
        MultiRiskMgrInterfaceLinkIn() &&
        CrossGamma::TYPE          &&
        FXCrossGamma::TYPE        &&
        VegaMatrixLoad()          &&
        ScenarioInterface::TYPE   &&
        SpotLevel::TYPE           &&
        VegaMatrix::TYPE          &&
        ThetaFwdSpot::TYPE        &&
        ThetaFwdRate::TYPE        &&
        ThetaNIE::TYPE            &&
        DeltaNextDay::TYPE        &&
        SpotPrice::TYPE           &&
        BasketDelta::TYPE         &&
        RollingTheta::TYPE        &&
        DeltaProxy::TYPE          &&
        VegaProxyParallel::TYPE   &&
        DeltaProxyNextDay::TYPE   &&
        SpotPriceProxy::TYPE      &&
        ImpliedScalarShift::TYPE  &&
        ImpliedVol::TYPE          &&
        ImpliedYTM::TYPE                  &&
        ImpliedYTP::TYPE                  &&
        CreditSpreadRhoParallel::TYPE     &&
        CreditSpreadRhoPointwise::TYPE    &&
        ParSpreadRhoParallelLinkIn()      &&
        ParSpreadRhoPointwiseLinkIn()     &&
        ParSpreadUpfrontPointwiseLinkIn()    &&
        ParSpreadUpfrontParallelLinkIn()    &&
        PrepayParallelLinkIn()               &&
        ParSpreadRhoParallelTwoSidedLinkIn() &&
        LegalBasisAdditiveParallelLinkIn()   &&
        LegalBasisAdditivePointwiseLinkIn()  &&
        LegalBasisMultiplierParallelLinkIn() &&
        LegalBasisMultiplierPointwiseLinkIn()&&
        LegalBasisMultiplierRecoveryLinkIn() &&
        LegalBasisAdditiveRecoveryLinkIn()   &&
        NakedBondRhoParallel::TYPE        &&
        NakedBondRhoPointwise::TYPE       &&
        CreditSpreadLevel::TYPE           &&
        OutputRequestOCBInt::TYPE         &&
        IHandlePaymentEvents::TYPE        &&
        CalibratorLinkIn()                &&
        ParamSensLinkIn()                 &&
        IScaleOutputs::TYPE               &&
        LiquiditySpreadRhoParallel::TYPE  &&
        LiquiditySpreadRhoPointwiseLinkIn() &&
        DeltaToCredit::TYPE               &&
        DeltaToCreditNextDay::TYPE        &&
        AssetVegaParallel::TYPE           &&
        AssetVegaPointwise::TYPE          &&
        LiquiditySpreadRhoParallel::TYPE  &&
        AdjCreditSpreadRhoParallel::TYPE  &&
        AdjCreditSpreadRhoPointwise::TYPE &&
        IE2CModel::TYPE                   &&
        OptionStrikeRhoParallel::TYPE     &&
        OptionStrikeRhoPointwise::TYPE    &&
        IIntrinsicMTM::TYPE               &&
        EquityVegaEquivalentParallel::TYPE &&
        CompositeModel::TYPE              &&
        VegaParallel2Sided::TYPE          &&
        DeltaDDE::TYPE                    &&
        DDeltaDVol::TYPE                  &&
        DeltaTPlusN::TYPE                 &&
        DeltaTPlusNHack::TYPE             &&
        BasketDeltaNextDayLinkIn()        &&
        ShortTermSqueezeParallelLinkIn()  &&
        LongTermSqueezeParallelLinkIn()   &&
        CorrParallelShiftLinkIn()         &&
        VegaWeightedPhi::TYPE             &&
        CreditDefaultSensBase::TYPE       &&
        ClosedForm::TYPE                  &&
        CorrelationShiftLinkIn()          &&
        CCMDBetaDSpreadLinkIn()           &&
        CCMDBasketBetaDSpreadLinkIn()     &&
        CompositeVegaParallelLinkIn()     &&
        CompositeVegaPointwiseLinkIn()    &&
        SpotVegaParallelLinkIn()          &&
        SpotVegaPointwiseLinkIn()         &&
        CurrencyBasisRhoParallelLinkIn()  &&
        CurrencyBasisRhoPointwiseLinkIn() &&
        CRSpotVegaParallelLinkIn()        &&
        CurrencyBasisSpreadLevelLinkIn()  &&
        CRMeanReversionParallelLinkIn()   &&
        IRMeanReversionParallelLinkIn()   &&
        SqueezeParallelTweakLinkIn()      &&
        CompressionRatioTweakLinkIn()     &&
        StrikeMappingOverrideTweakLinkIn()&&
        StrikeMappingTweakLinkIn()        &&
        RateShiftLinkIn()                 &&
        YCWeightedAdditiveShiftLinkIn()   &&
        YCWeightedMultiplicativeShiftLinkIn()         &&
        ParSpreadWeightedAdditiveShiftLinkIn()        &&
        ParSpreadWeightedMultiplicativeShiftLinkIn()  &&
        BetaSkewMatrixTweakLinkIn()                   &&
        RecoveryTweakLinkIn()                         &&
        RecoveryShiftLinkIn()                         &&
        IRVegaMatrixLinkIn()                          &&
        VHTSkewFactorsPointwise::TYPE                 &&
        VHTSkewCelerityFactorsPointwise::TYPE         &&
        VHTConvexityFactorsPointwise::TYPE            &&
        VHTConvexityCelerityFactorsPointwise::TYPE    &&
        RecoverySetLinkIn()                           &&
        BCStrikeMappingOverrideTweakLinkIn()          &&
        BCStrikeMappingSetLinkIn()                    &&
        BCStrikeMappingTweakLinkIn()                  &&
        RecoveryTweakWithShiftLinkIn()                &&
        CreditDeltaPointwiseWithShiftLinkIn()         &&
        EnergyDeltaPointwiseWithShiftLinkIn()         &&
        CreditIndexBasisRhoParallelLinkIn()           &&
        ExposureReporterLinkIn()                      &&
        CreditIndexBasisRhoPointwiseLinkIn()          &&
        PerEntryFlexiblePerturbationLinkIn()          &&
        ScalarFlexiblePerturbationLinkIn()            &&
        FlexibleSensitivityLinkIn()                   &&
        RhoParallelTwoSidedLinkIn()                   &&
        RhoPointwiseTwoSidedLinkIn()                  &&
        PerturbationTestLinkIn()                      &&
        IRDeltaPointwiseLinkIn()                      &&
        CorrSwapBasisAdjTweakLinkIn()                 &&
        CorrSwapBasisAdjParallelLinkIn()              &&
        CorrSwapSamplingAdjTweakLinkIn()              &&
        IBootstrapper::TYPE                           &&
        HolidayLoad()                                 &&
        QuasiContractualBetaSkewParallelLinkIn()      &&
        QuasiContractualBetaSkewMatrixLinkIn()        &&
        CRVegaParallelLinkIn()                        &&  
        CRVegaPointwiseLinkIn()                       &&  
        ModelFilterLinkIn()                           &&
//        CreditIndexSpreadRhoPointwiseLinkIn()         &&
        VegaAtmParallelConstConvxLinkIn()             &&
        VegaAtmPointwiseConstConvxLinkIn()            &&
        MarketDataNameSubstitutionLinkIn()            &&
        Smile2QElementwiseLinkIn()                    &&
        CreditMultiBetaLevelLinkIn()                  &&
        ImpliedScalarShiftMultiLinkIn()               &&
        CreditNameNotionalLevelLinkIn()               &&
        MultiModelLinkIn()                            &&
        WeightedInstrumentShiftLinkIn()               &&
        CDOParallelStrikeShiftLinkIn()                &&
        CreditFeeLegCouponShiftLinkIn()               &&
        WeightedInstrumentCollectionLinkIn()          &&
        LegalBasisAdditiveRelativeParallelLinkIn()    &&
        ClosedFormForwardRatePricerLoad()             &&
        TierLoad()                                    &&
        CompositeForwardRatePricerLoad()              &&
        CorridorVarSwapSensLinkIn()                   &&
        SPIGapRiskSensLoad()                          &&
        MOQTestMarketObjectLoad()                     &&
        MOQTestModelLoad()                            &&
        PhiLinkIn()                                   &&
        PhiParallelLinkIn()                           &&
        MarketObjectQualifiersLinkIn()                &&
        true;

    if (!success){
        throw ModelException("CRiskMgrLib::registerClasses",
                             "Registration error");
    }

}
DRLIB_END_NAMESPACE
