//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : ProductsLib.cpp
//
//   Description : Force linker to include all relevant files/symbols
//
//   Author      : Mark A Robson
//
//   Date        : 2 March 2001
//
//
//
// This file is not itself compiled and linked into QLib.  In order to support
// selective builds*, the system generates a modified version
// ProductsLib-filtered.cpp, in which references to "FooLoad()" are removed
// unless Foo.cpp is selected for inclusion.  It's ProductsLib-filtered.cpp
// which is actually built into the library.
//
// (NB don't edit ProductsLib-filtered.cpp, since any changes you make to it
// will be overwritten --- edit this file.)
//
// [Technical note: the filtering is performed by
// ../../../makerules/scripts/filterProductsLib.pl, invoked from
// ../../../makerules/gnu/selective-srcs.mkh and from
// ../../QLibSolution.vsmproj:QLibPartial.checkDeps().]
//
// ---
//   *see selectivebuild.example for more info

#include "edginc/config.hpp"
#include "edginc/ProductsLib.hpp"
#include "edginc/ModelException.hpp"
#include "edginc/ToolkitDebug.hpp"


DRLIB_BEGIN_NAMESPACE
extern bool RainbowSPILoad();
extern bool NFBinaryLoad();
extern bool GenericNFactorTesterLoad();
extern bool TaxWrapperLoad();
extern bool CallableEquityKOSwapLoad();
extern bool CorpActLoad();
extern bool AccumulatorLoad();
extern bool AllWeatherTARNLoad();
extern bool CreditMetricsLondonFloorLoad();
extern bool EnergySwapCouponLoad();
extern bool TracerLoad();
extern bool BaskAvLoad();
extern bool CDblBarrierLoad();
extern bool VanillaSmoothStrikeLoad();
extern bool Generic1FactorTesterLoad();
extern bool CredDefSwapLoad();
extern bool CashFlowStreamLoad();
extern bool CPVDivLoad();
extern bool InsuranceAnnuityGMBLoad();
extern bool InsuranceALoad();
extern bool InsuranceAGMWB1Load();
extern bool InsuranceAGMWB2Load();
extern bool ImpliedLossAddInLoad();
extern bool FloatingBondLoad();
extern bool BondFloatNtlLoad();
extern bool EGKnockInLoad();
extern bool PickRainbowLoad();
extern bool LadderAverageLoad();
extern bool SyntheticPortfolioInsuranceLoad();
extern bool OptOnConvBondLoad();
extern bool GMinInvBenefitLoad();
extern bool ClosedFormEnergyLoad();
extern bool BarrierLoad();
extern bool CalendarDropLoad();
extern bool EnergyCapFloorLoad();
extern bool DropRainbowLoad();
extern bool QPVanillaLoad();
extern bool InstallmentLoad();
extern bool ImpliedIntegrationLoad();
extern bool CForwardContractLoad();
extern bool ClosedFormFALoad();
extern bool ECOLoad();
extern bool TrailFeeLoad();
extern bool RealizedCorrelationLoad();
extern bool AUD3YBondFutureLoad();
extern bool VWAPLoad();
extern bool CliquetPerformanceLoad();
extern bool CDOLoad();
extern bool GeneralisedCDOLoad();
extern bool CEqGainKONoteLoad();
extern bool ConvBondLoad();
extern bool FRLoad();
extern bool FRBarrierVariableLoad();
extern bool EscalatorLoad();
extern bool VanillaTSOLoad();
extern bool MomentumLoad();
extern bool Abs2DeltaBasedStrikeVolConverterLoad();
extern bool BaskRebalancedLoad();
extern bool BestOrWorstLoad();
extern bool RiskyBondSeriesLoad();
extern bool TrancheIndexLeastSquareFitLoad();
extern bool PhysicalSettlementOverrideNameLoad();
extern bool SplineTrancheQuoteInterpolatorLoad();
extern bool DividendAdjustedFLLoad();
extern bool CreditEventOverrideLoad();
extern bool CBCashFlowLoad();
extern bool RollingAveAmerLoad();
extern bool RangeNoteLoad();
extern bool EquityStabilitySwapLoad();
extern bool VanillaGridMultiLoad();
extern bool AssetValueLoad();
extern bool CreditEventOverrideNameLoad();
extern bool VanillaGridLoad();
extern bool SwaptionLoad();
extern bool EnergyFixingLoad();
extern bool CForwardContractRiskyLoad();
extern bool VolCheckLoad();
extern bool MCDemoProductLoad();
extern bool EquityDistressSwapLoad();
extern bool StrikeResetLoad();
extern bool CEquitySwapLoad();
extern bool NumericalIntegrationLNLoad();
extern bool CorrCovSwapLoad();
extern bool VanillaCDSLoad();
extern bool OptionOnIRFLoad();
extern bool IndexBasisCalcAddinLoad();
extern bool TriggerECOLoad();
extern bool SyntheticConvertLoad();
extern bool CallableEDSLoad();
extern bool IRCapLoad();
extern bool CExtendableNoteLoad();
extern bool ForwardOptionSeriesLoad();
extern bool CCMBaseCorrelationLoad();
extern bool ClosedFormIRLNLoad();
extern bool ClosedFormCDSBasketLoad();
extern bool LiborStreamLoad();
extern bool GammaSwapLoad();
extern bool SuperRainbowLoad();
extern bool AsianASWStrikeLoad();
extern bool BasketOfLookBacksLoad();
extern bool AUD90DayBondFutureLoad();
extern bool CFDGridPassLoad();
extern bool TimingRuleTesterLoad();
extern bool EntropyExpLossInterpolatorLoad();
extern bool BondFutureLoad();
extern bool CashSettlementOverrideNameLoad();
extern bool DetailedCreditEventOverrideNameLoad();
extern bool EnhancedStragglerLoad();
extern bool VIXHistLoad();
extern bool CollarGammaLoad();
extern bool AverageLockerLoad();
extern bool DividendAdjustedLoad();
extern bool LadderLoad();
extern bool CredDefSwaptionLoad();
extern bool OptimoLoad();
extern bool BondParamsLoad();
extern bool SyntheticAssetLoad();
extern bool EnergySwapScheduleLoad();
extern bool TrendOptionLoad();
extern bool PartCashSettlementOverrideNameLoad();
extern bool COPPERLoad();
extern bool VanillaCreditFeeLegLoad();
extern bool EnergyVanillaLoad();
extern bool CreditIndexSwapLoad();
extern bool AbstractCashFlowStreamLoad();
extern bool RainbowKOLoad();
extern bool VanillaGridCEVJLoad();
extern bool VolatilitySwapLoad();
extern bool AmerSpreadLoad();
extern bool RainbowLoad();
extern bool CVanillaLoad();
extern bool GenericCoreDumpLoad();
extern bool ExplicitBondLoad();
extern bool SingleCashFlowLoad();
extern bool BoostedNFBLoad();
extern bool FlatExpLossPriorLoad();
extern bool CFuturesLoad();
extern bool CallOnCliqLoad();
extern bool AssetCVBLoad();
extern bool CorporateBondLoad();
extern bool KnockInFwdLoad();
extern bool CalendarRainbowLoad();
extern bool EnergyFutureSwapScheduleLoad();
extern bool FastQuoteLoad();
extern bool CreditKOXCcySwapLoad();
extern bool EnergySwapLoad();
extern bool BondCashFlowsLoad();
extern bool CorridorVarSwapLoad();
extern bool RainbowDKOLoad();
extern bool CreditMetricsModelLoad();
extern bool PickAverageLoad();
extern bool IRFutureLoad();
extern bool AverageLoad();
extern bool ClosedFormCDSPSLoad();
extern bool CreditDefaultSwapLoad();
extern bool ExtendibleLoad();
extern bool VanillaGridInstrumentCollectionLoad();
extern bool VanillaCreditContingentLegLoad();
extern bool CDOFineGridLoad();
extern bool VarSwapHedgingSupportLoad();
extern bool RotatorLoad();
extern bool UniformDensityPriorLoad();
extern bool RainbowYenmanLoad();
extern bool ShrinkLoad();
extern bool RainbowRangeKOLoad();
extern bool FRLoad();
extern bool CStepDownBondLoad();
extern bool EGKBondLoad();
extern bool TargetRedemptionNoteLoad();
extern bool CCMLondonFloorLoad();
extern bool EnergyFutureLoad();
extern bool VolVarSwapLoad();
extern bool VarianceIndexForwardLoad();
extern bool MertonLVCalibLoad();
extern bool VanillaMomentLoad();
extern bool CClosedFormLNLoad();
extern bool VDaxLoad();
extern bool VAssetLoad();
extern bool VIXFutureLoad();
extern bool VForwardLoad();
extern bool AmericanLadderLoad();
extern bool RangeAccrueLoad();
extern bool OperaLoad();
extern bool PutCliquetLoad();
extern bool InsuranceAnnuityLoad();
extern bool CallableDepositLoad();
extern bool SimpathIInstrumentLoad();
extern bool BandedBaskAvLoad();
extern bool ASRLoad();
extern bool CCMLoad();
extern bool ForwardBasisIndexSeriesLoad();
extern bool CreditMetricsBaseCorrelationLoad();
extern bool BondLoad();
extern bool TrancheOptionLoad();
extern bool CreditMetricsRFLLoad();
extern bool CCMRFLLoad();
extern bool CreditMetricsABSCDOLoad();
extern bool CDSOptionLoad();
extern bool VarianceOptionLoad();
extern bool Q3MQQuasiPricerLoad();
extern bool CMCDSFloatingLegLoad();
extern bool MultiQQuasiSmileLoad();
extern bool LHPADensityPriorLoad();
extern bool RangeCounterLoad();
extern bool TrancheCashSettlementOverrideLoad();
extern bool TranchePhysicalSettlementOverrideLoad();
extern bool TranchePartCashSettlementOverrideLoad();
extern bool DeliveryDetailsLoad();
extern bool NoticeOfPhysicalSettlementLoad();
extern bool LookbackLoad();
extern bool AbstractionExampleLoad();
extern bool MCStatisticsProductsLoad();
extern bool PyramidLegacyFFXLoad();
extern bool SCIDLoad();
extern bool TranchePricerLoad();
extern bool InstalmentWarrantLoad();
extern bool MCTestingCIDMatricesLoad();
extern bool TranchePricerCalibratorLoad();
extern bool TranchePricerForCalibratorLoad();
extern bool BespokeCDOModelTestLoad();
extern bool BespokeCDOModelBCTestLoad();
extern bool SCIDCreditTARNLoad();
extern bool SCIDRiskyZeroLoad();
extern bool BasicCreditContingentLegLoad();
extern bool CDSLoad();
extern bool ClosedFormForwardRatePricerLoad();
extern bool CDOIndexOptionLoad();
extern bool TrancheOptionLoad();
extern bool CreditLegConventionLoad();
extern bool AccNewLoad();
extern bool ICDOQuotesGeneratorLoad();
extern bool CDOQuotesMaturityFilterLoad();
extern bool CDOQuotesDiscreteWeightsLoad();
extern bool CDOQuotesGeneratorLoad();
extern bool CDOQuotesBCGeneratorLoad();
extern bool HybNVanillaLoad();
extern bool DailyAccumulatorLoad();
extern bool BasketAverageLoad();

void CProductsLib::linkInClasses(){
    /* The sole purpose of this method is to ensure that the linker
       includes all symbols out of the products directory. Many symbols
       are automatically linked because they are used by other classes
       which are already included.

       An example of symbols that could be dropped would be an entire class
       representing a product. For example Vanilla might be dropped since
       the only references might be throught the abstract parent class. */


    /* there is no order dependency - this could be automated
       The load functions below belong to auxiliary classes
       and hence are included regardless of which products are selected to
       build.  The load functions for all other products which are
       specified in the build specification are automatically added
       by the make scripts.
    */

    // Please don't add any new ...::TYPE entries here ---
    // use FFFLoad() where FFF.cpp is the source filename.

    bool success =
        RainbowSPILoad() &&
        NFBinaryLoad() &&
        GenericNFactorTesterLoad() &&
        TaxWrapperLoad() &&
        CallableEquityKOSwapLoad() &&
        CorpActLoad() &&
        AccumulatorLoad() &&
        AllWeatherTARNLoad() &&
        CreditMetricsLondonFloorLoad() &&
        EnergySwapCouponLoad() &&
        TracerLoad() &&
        BaskAvLoad() &&
        CDblBarrierLoad() &&
        VanillaSmoothStrikeLoad() &&
        Generic1FactorTesterLoad() &&
        CredDefSwapLoad() &&
        CashFlowStreamLoad() &&
        CPVDivLoad() &&
        InsuranceAnnuityGMBLoad() &&
        InsuranceALoad() &&
        InsuranceAGMWB1Load() &&
        InsuranceAGMWB2Load() &&
        ImpliedLossAddInLoad() &&
        FloatingBondLoad() &&
        BondFloatNtlLoad() &&
        EGKnockInLoad() &&
        PickRainbowLoad() &&
        LadderAverageLoad() &&
        SyntheticPortfolioInsuranceLoad() &&
        OptOnConvBondLoad() &&
        GMinInvBenefitLoad() &&
        ClosedFormEnergyLoad() &&
        BarrierLoad() &&
        CalendarDropLoad() &&
        EnergyCapFloorLoad() &&
        DropRainbowLoad() &&
        QPVanillaLoad() &&
        InstallmentLoad() &&
        ImpliedIntegrationLoad() &&
        CForwardContractLoad() &&
        ClosedFormFALoad() &&
        ECOLoad() &&
        TrailFeeLoad() &&
        RealizedCorrelationLoad() &&
        AUD3YBondFutureLoad() &&
        VWAPLoad() &&
        CliquetPerformanceLoad() &&
        CDOLoad() &&
        GeneralisedCDOLoad() &&
        CEqGainKONoteLoad() &&
        ConvBondLoad() &&
        FRLoad() &&
        FRBarrierVariableLoad() &&
        EscalatorLoad() &&
        VanillaTSOLoad() &&
        MomentumLoad() &&
        Abs2DeltaBasedStrikeVolConverterLoad() &&
        BaskRebalancedLoad() &&
        BestOrWorstLoad() &&
        RiskyBondSeriesLoad() &&
        TrancheIndexLeastSquareFitLoad() &&
        PhysicalSettlementOverrideNameLoad() &&
        SplineTrancheQuoteInterpolatorLoad() &&
        DividendAdjustedFLLoad() &&
        CreditEventOverrideLoad() &&
        CBCashFlowLoad() &&
        RollingAveAmerLoad() &&
        RangeNoteLoad() &&
        EquityStabilitySwapLoad() &&
        VanillaGridMultiLoad() &&
        AssetValueLoad() &&
        CreditEventOverrideNameLoad() &&
        VanillaGridLoad() &&
        SwaptionLoad() &&
        EnergyFixingLoad() &&
        CForwardContractRiskyLoad() &&
        VolCheckLoad() &&
        MCDemoProductLoad() &&
        EquityDistressSwapLoad() &&
        StrikeResetLoad() &&
        CEquitySwapLoad() &&
        NumericalIntegrationLNLoad() &&
        CorrCovSwapLoad() &&
        VanillaCDSLoad() &&
        OptionOnIRFLoad() &&
        IndexBasisCalcAddinLoad() &&
        TriggerECOLoad() &&
        SyntheticConvertLoad() &&
        CallableEDSLoad() &&
        IRCapLoad() &&
        CExtendableNoteLoad() &&
        ForwardOptionSeriesLoad() &&
        CCMBaseCorrelationLoad() &&
        ClosedFormIRLNLoad() &&
        ClosedFormCDSBasketLoad() &&
        LiborStreamLoad() &&
        GammaSwapLoad() &&
        SuperRainbowLoad() &&
        AsianASWStrikeLoad() &&
        BasketOfLookBacksLoad() &&
        AUD90DayBondFutureLoad() &&
        CFDGridPassLoad() &&
        TimingRuleTesterLoad() &&
        EntropyExpLossInterpolatorLoad() &&
        BondFutureLoad() &&
        CashSettlementOverrideNameLoad() &&
        DetailedCreditEventOverrideNameLoad() &&
        EnhancedStragglerLoad() &&
        VIXHistLoad() &&
        CollarGammaLoad() &&
        AverageLockerLoad() &&
        DividendAdjustedLoad() &&
        LadderLoad() &&
        CredDefSwaptionLoad() &&
        OptimoLoad() &&
        BondParamsLoad() &&
        SyntheticAssetLoad() &&
        EnergySwapScheduleLoad() &&
        TrendOptionLoad() &&
        PartCashSettlementOverrideNameLoad() &&
        COPPERLoad() &&
        VanillaCreditFeeLegLoad() &&
        EnergyVanillaLoad() &&
        CreditIndexSwapLoad() &&
        AbstractCashFlowStreamLoad() &&
        RainbowKOLoad() &&
        VanillaGridCEVJLoad() &&
        VolatilitySwapLoad() &&
        AmerSpreadLoad() &&
        RainbowLoad() &&
        CVanillaLoad() &&
        GenericCoreDumpLoad() &&
        ExplicitBondLoad() &&
        SingleCashFlowLoad() &&
        BoostedNFBLoad() &&
        FlatExpLossPriorLoad() &&
        CFuturesLoad() &&
        CallOnCliqLoad() &&
        AssetCVBLoad() &&
        CorporateBondLoad() &&
        KnockInFwdLoad() &&
        CalendarRainbowLoad() &&
        EnergyFutureSwapScheduleLoad() &&
        FastQuoteLoad() &&
        CreditKOXCcySwapLoad() &&
        EnergySwapLoad() &&
        BondCashFlowsLoad() &&
        CorridorVarSwapLoad() &&
        RainbowDKOLoad() &&
        CreditMetricsModelLoad() &&
        PickAverageLoad() &&
        IRFutureLoad() &&
        AverageLoad() &&
        ClosedFormCDSPSLoad() &&
        CreditDefaultSwapLoad() &&
        ExtendibleLoad() &&
        VanillaGridInstrumentCollectionLoad() &&
        VanillaCreditContingentLegLoad() &&
        CDOFineGridLoad() &&
        VarSwapHedgingSupportLoad() &&
        RotatorLoad() &&
        UniformDensityPriorLoad() &&
        RainbowYenmanLoad() &&
        ShrinkLoad() &&
        RainbowRangeKOLoad() &&
        FRLoad() &&
        CStepDownBondLoad() &&
        EGKBondLoad() &&
        TargetRedemptionNoteLoad() &&
        CCMLondonFloorLoad() &&
        EnergyFutureLoad() &&
        VolVarSwapLoad() &&
        VarianceIndexForwardLoad() &&
        MertonLVCalibLoad() &&
        VanillaMomentLoad() &&
        CClosedFormLNLoad() &&
        VDaxLoad() &&
        VAssetLoad() &&
        VIXFutureLoad() &&
        VForwardLoad() &&
        AmericanLadderLoad() &&
        RangeAccrueLoad() &&
        OperaLoad() &&
        PutCliquetLoad() &&
        InsuranceAnnuityLoad() &&
        CallableDepositLoad() &&
        SimpathIInstrumentLoad() &&
        BandedBaskAvLoad() &&
        ASRLoad() &&
        CCMLoad() &&
        CreditMetricsBaseCorrelationLoad() &&
        BondLoad() &&
        TrancheOptionLoad() &&
        CreditMetricsRFLLoad() &&
        CCMRFLLoad() &&
        ForwardBasisIndexSeriesLoad() &&
        CreditMetricsABSCDOLoad() &&
        CDSOptionLoad() &&
        VarianceOptionLoad() &&
      Q3MQQuasiPricerLoad() &&
       MultiQQuasiSmileLoad() &&
      CMCDSFloatingLegLoad() &&
        LHPADensityPriorLoad() &&
        RangeCounterLoad() &&
        TrancheCashSettlementOverrideLoad() &&
        TranchePhysicalSettlementOverrideLoad() &&
        TranchePartCashSettlementOverrideLoad() &&
        DeliveryDetailsLoad() &&
        NoticeOfPhysicalSettlementLoad() &&
        LookbackLoad() &&
        AbstractionExampleLoad() &&
        MCStatisticsProductsLoad() &&
        PyramidLegacyFFXLoad() &&
        SCIDLoad() &&
        TranchePricerLoad() &&
        InstalmentWarrantLoad() &&
        MCTestingCIDMatricesLoad() &&
		TranchePricerCalibratorLoad() &&
        TranchePricerForCalibratorLoad() &&
        BespokeCDOModelTestLoad() &&
        BespokeCDOModelBCTestLoad() &&
		SCIDCreditTARNLoad() &&
		SCIDRiskyZeroLoad() &&
        BasicCreditContingentLegLoad() &&
        CDSLoad() &&
        ClosedFormForwardRatePricerLoad() &&
        CDOIndexOptionLoad() &&
        TrancheOptionLoad() &&
        CreditLegConventionLoad() &&
        AccNewLoad() &&
        ICDOQuotesGeneratorLoad() &
        CDOQuotesGeneratorLoad() &
        CDOQuotesBCGeneratorLoad() &
        CDOQuotesMaturityFilterLoad() &&
        CDOQuotesDiscreteWeightsLoad() &&
        HybNVanillaLoad() &&
        DailyAccumulatorLoad() &&
        BasketAverageLoad() &&
        true;


    if (!success){
        throw ModelException("CProductsLib::registerClasses",
                             "Registration error");
    }
}

DRLIB_END_NAMESPACE
