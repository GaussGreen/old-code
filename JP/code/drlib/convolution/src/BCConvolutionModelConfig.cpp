#include "edginc/config.hpp"
#include "edginc/BCConvolutionModelConfig.hpp"
#include "edginc/BCPortfolioNameConvolutionLossGen.hpp"
#include "edginc/BCTrancheConvolutionLossGen.hpp"
#include "edginc/SingleCreditAsset.hpp"
#include "edginc/FlatCDO2LossConfig.hpp"
#include "edginc/NToDefaultLossConfig.hpp"
#include "edginc/CreditAsset.hpp"
#include "edginc/CounterPartyCredit.hpp"
#include "edginc/FtDEffectiveCurveGen.hpp"

DRLIB_BEGIN_NAMESPACE

DEFINE_TEMPLATE_TYPE(BCConvolutionModelConfigArray);

/** TYPE (for reflection) */        
CClassConstSP const BCConvolutionModelConfig::TYPE =
CClass::registerClassLoadMethod(
    "BCConvolutionModelConfig",
    typeid(BCConvolutionModelConfig),
    BCConvolutionModelConfig::load);

/** Virtual destructor */
BCConvolutionModelConfig::~BCConvolutionModelConfig()
{
}

/** Constructor */
BCConvolutionModelConfig::BCConvolutionModelConfig():
    CObject(TYPE),
    compressionRatio(DEFAULT_THETA),
    strikeMappingOverride(0),
    lowerBCBetasOverride(0),
    upperBCBetasOverride(0),
    bcBetasTimeline(0),
    spreadRatioOverride(0),
    spreadRatiosRegionalOverride(0),
    spreadRatiosRegions(0),
    historicalBetaAdjustment(0),
    namesOverride(new StringArray(0)),
    indexWeightsOverride(new IndexWeightsArray(0)),
    authoriseNegativeEL(false),
    useExpectedLossRatio(false){}

/** Explicit constructor */
BCConvolutionModelConfig::BCConvolutionModelConfig(
    IConditionalDefaultsModelSP condDefaultsModel,
    IConvolutorSP convolutor,
    bool returnBinaryDistribution,
    double compressionRatio,
    CDoubleSP strikeMappingOverride,
    DoubleArraySP lowerBCBetasOverride,
    DoubleArraySP upperBCBetasOverride,
    DateTimeArraySP bcBetasTimeline,
    CDoubleSP spreadRatioOverride,
    DoubleArraySP spreadRatiosRegionalOverride,
    StringArraySP spreadRatiosRegions,
    double historicalBetaAdjustment,
    StringArraySP namesOverride,
    IndexWeightsArraySP indexWeightsOverride,
    bool authoriseNegativeEL,
    bool useExpectedLossRatio):
        CObject(TYPE),
        convolutor(convolutor),
        returnBinaryDistribution(returnBinaryDistribution),
        condDefaultsModel(condDefaultsModel),
        compressionRatio(compressionRatio),
        strikeMappingOverride(strikeMappingOverride),
        lowerBCBetasOverride(lowerBCBetasOverride),
        upperBCBetasOverride(upperBCBetasOverride),
        bcBetasTimeline(bcBetasTimeline),
        spreadRatioOverride(spreadRatioOverride),
        spreadRatiosRegionalOverride(spreadRatiosRegionalOverride),
        spreadRatiosRegions(spreadRatiosRegions),
        historicalBetaAdjustment(historicalBetaAdjustment),
        namesOverride(namesOverride),
        indexWeightsOverride(indexWeightsOverride),
        authoriseNegativeEL(authoriseNegativeEL),
        useExpectedLossRatio(useExpectedLossRatio){}

IObject* BCConvolutionModelConfig::defaultConstructor()
{
    return new BCConvolutionModelConfig();
}

void BCConvolutionModelConfig::load(CClassSP& clazz)
{
    clazz->setPublic();
    REGISTER(BCConvolutionModelConfig, clazz);
    IMPLEMENTS(ICreditLossModelConfig);
    IMPLEMENTS(IModelConfigMapper);
    SUPERCLASS(CObject);
    EMPTY_SHELL_METHOD(defaultConstructor);
    FIELD(convolutor, "Convolution engine");
    FIELD(returnBinaryDistribution,
        "Flag to control if we want to approximate the returned"
        "loss distribution by a binary distribution");
    FIELD(condDefaultsModel, "Conditional defaults model - used for e.g. CDO2");

    FIELD(compressionRatio,
                 "Compression ratio (also called theta) "
                 "used to rescale the historical betas dispersion");
    FIELD_MAKE_OPTIONAL(compressionRatio);
    FIELD(strikeMappingOverride,
          "Optional global strike mapping parameter (also called 'q')"
          " to override all strike mapping parameters in index skew matrix");
    FIELD_MAKE_OPTIONAL(strikeMappingOverride);
    FIELD(lowerBCBetasOverride, "Optional override for lower strike base correlation betas");
    FIELD_MAKE_OPTIONAL(lowerBCBetasOverride);
    FIELD(upperBCBetasOverride, "Optional override for upper strike base correlation betas");
    FIELD_MAKE_OPTIONAL(upperBCBetasOverride);
    FIELD(bcBetasTimeline,
          "Timeline associated with lowerBCBetasOverride and upperBCBetasOverride."
          "This input is used for checking purpose only (bcBetasTimeline must be"
          "exactly the same as the effective curve timeline computed internally)");
    FIELD_MAKE_OPTIONAL(bcBetasTimeline);
    FIELD(spreadRatioOverride,
          "Optional global parameter to override the spread ratio used "
          "for that trade, typically use 1 for index trades");
    FIELD_MAKE_OPTIONAL(spreadRatioOverride);
    FIELD(spreadRatiosRegionalOverride,
        "Optional override: allow to specify a 'local' spread ratio override per "
        "region (corresponding regions are defined by 'spreadRatiosRegions' field). "
        "'spreadRatioOverride' and 'spreadRatiosRegionalOverride' are exclusive");
    FIELD_MAKE_OPTIONAL(spreadRatiosRegionalOverride);
    FIELD(spreadRatiosRegions,
        "Define regions corresponding to 'spreadRatiosRegionalOverride'");
    FIELD_MAKE_OPTIONAL(spreadRatiosRegions);
    FIELD(historicalBetaAdjustment, "Historical Beta adjustment");
    FIELD_MAKE_OPTIONAL(historicalBetaAdjustment);
    FIELD(namesOverride,
          "Array of names (in the portfolio) for which we want to "
          "override the associated indexWeights");
    FIELD_MAKE_OPTIONAL(namesOverride);
    FIELD(indexWeightsOverride,
          "Array of overriden indexWeights");
    FIELD_MAKE_OPTIONAL(indexWeightsOverride);
    FIELD(authoriseNegativeEL,
        "TRUE=authorise negative expected losses, "
        "FALSE=do not authorise negative expected losses (floors to zero) "
        "[default is FALSE]");
    FIELD_MAKE_OPTIONAL(authoriseNegativeEL);
    FIELD(useExpectedLossRatio,
        "TRUE=use expected loss ratio instead of spread ratio in the spread mapping, "
        "FALSE=use par spread mapping "
        "[default is FALSE]");
    FIELD_MAKE_OPTIONAL(useExpectedLossRatio);
    FIELD(discountCurve,
        "Instrument discount curve - only used to compute duration "
        "weighted average spread when using 'par spreads' methodology");
    FIELD_MAKE_OPTIONAL(discountCurve);
}

/** Default value for compression ratio (theta) */
double const BCConvolutionModelConfig::DEFAULT_THETA = 1.0;

/** Dummy implementation of IModelConfigMapper (just returns "this") */
ICreditLossModelConfigConstSP BCConvolutionModelConfig::innerModel(
    ICreditLossConfigConstSP lossConfig) const
{
    return ICreditLossModelConfigConstSP(this);
}

ICreditLossGenSP BCConvolutionModelConfig::lossGenerator(
    ICreditLossConfigConstSP lossConfig,
    IModelConfigMapperConstSP mapper) const
{
    // Expects lossConfig to be a CreditTrancheLossConfig.
    // We could support more than one implementation of
    // ICreditLossConfig (using double dispatch mechanism) but not sure
    // there is too much to share / factorise between radically different
    // loss configs (eg: PortfolioName vs CreditTrancheLossConfig).
    // We already have lots of genericity in the recursive definition
    // of loss configs.

    if (CreditTrancheLossConfig::TYPE->isInstance(lossConfig) ||
		FlatCDO2LossConfig::TYPE->isInstance(lossConfig) )
    {
        const CreditTrancheLossConfig* tranche =
            dynamic_cast<const CreditTrancheLossConfig*>( lossConfig.get());

		const FlatCDO2LossConfig* flatCDO2 =
            dynamic_cast<const FlatCDO2LossConfig*>( lossConfig.get());

		const ITranchesCombinationPayoff* trancheComb = 
			dynamic_cast<const ITranchesCombinationPayoff*>( lossConfig.get() );

		const ICreditLossConfig * portfolio;
		DateTime today;
		double outstandingNotional;

		if (tranche)
		{
			today = tranche->getToday();
			portfolio = tranche->getPortfolio().get();
			outstandingNotional = tranche->trancheOutstandingNotional(false);

		} else {
			today = flatCDO2->getToday();
			portfolio = flatCDO2->getPortfolio().get();
			outstandingNotional = flatCDO2->trancheOutstandingNotional(false);

		};

        return ICreditLossGenSP(new BCTrancheConvolutionLossGen(
            trancheComb,
            outstandingNotional,
            portfolio,
            today,
            mapper.get(),
            convolutor.get(),
            returnBinaryDistribution,
            condDefaultsModel.get(),
            compressionRatio,
            strikeMappingOverride,
            lowerBCBetasOverride,
            upperBCBetasOverride,
            bcBetasTimeline,
            spreadRatioOverride,
            spreadRatiosRegionalOverride,
            spreadRatiosRegions,
            historicalBetaAdjustment,
            namesOverride,
            indexWeightsOverride,
            authoriseNegativeEL,
            useExpectedLossRatio,
            discountCurve));
    }
    else if (PortfolioName::TYPE->isInstance(lossConfig))
    {
        return ICreditLossGenSP(new BCPortfolioNameConvolutionLossGen(
            DYNAMIC_CONST_CAST(PortfolioName, lossConfig.get()),
            condDefaultsModel.get()));
    }
    else
    {
        throw ModelException(
            "BCConvolutionModelConfig does not handle " +
            lossConfig->getClass()->getName());
    }
}


/** Creates an IEffectiveCurveGen for an 'NtD' */
ICreditLossGenSP BCConvolutionModelConfig::effCurveGenerator(
    ICreditLossConfigConstSP  lossConfig,
    CounterPartyCreditConstSP cpty,
    const bool                recoverNotional,
    IModelConfigMapperConstSP modelConfigMapper) const
{
    static const string method("BCConvolutionModelConfig::effCurveGenerator");
                         

    if (NToDefaultLossConfig::TYPE->isInstance(lossConfig)) {
        throw ModelException (method, "This model config can NOT be used "
                              "to price NtDs.");
    }
    else {
        throw ModelException(method, "BCConvolutionModelConfig does not handle " +
            lossConfig->getClass()->getName());
    }
}


// Hack to pass discount curve - note the input discount curve is cloned
void BCConvolutionModelConfig::setDiscountCurve(
    YieldCurveConstSP discount)
{
    if (discountCurve.get() == 0)
    {
        discountCurve.reset(
            DYNAMIC_CAST(YieldCurve, discount->clone()));
    }
    else
    {
        throw ModelException(
            "BCConvolutionModelConfig::setDiscountCurve",
            "Discount curve already populated: can not override.");
    }
}

/* external symbol to allow class to be forced to be linked in */
bool BCConvolutionModelConfigLoad(){
    return (BCConvolutionModelConfig::TYPE != 0);
}

DRLIB_END_NAMESPACE
