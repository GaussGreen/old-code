#include "edginc/config.hpp"
#include "edginc/BCPortfolioNameConvolutionLossGen.hpp"
#include "edginc/SingleCreditAsset.hpp"
#include "edginc/DiscreteDistribution.hpp"
#include "edginc/BaseCorrelationOnlyParameters.hpp"
#include "edginc/CmOnlyParameters.hpp"
#include "edginc/Format.hpp"
#include "edginc/LinearInterpolator.hpp"

DRLIB_BEGIN_NAMESPACE

#define TINY 3e-15

BCPortfolioNameConvolutionLossGen::~BCPortfolioNameConvolutionLossGen() {}

BCPortfolioNameConvolutionLossGen::BCPortfolioNameConvolutionLossGen(
    const PortfolioName* ptfName,
    const IConditionalDefaultsModel* condDefaultsModel):
        ptfName(ptfName),
        condDefaultsModel(condDefaultsModel)
{
    // Retrieves model parameters
    CreditEngineParametersConstSP params(
        ptfName->getEngineParams(BaseCorrelationOnlyParameters::TYPE));
    BaseCorrelationOnlyParametersConstSP bcParam =
        BaseCorrelationOnlyParametersConstSP::dynamicCast(params); 
    indexWeights = bcParam->getIndexWeights();
    
    // Retrieves underlying "single name" CDS
    SingleCreditAssetConstSP myAsset =
        SingleCreditAssetConstSP::dynamicCast(ptfName->getCreditAsset());
    const ICDSParSpreads* cds = myAsset->cdsParSpreads.get();
    
    // Retrieves clean spreads
    defaultRates = cds->defaultRates();
}

/**
 * Mechanism used to pass the "single name - index" mapping to
 * the tranche level
 * */
IndexWeightsConstSP BCPortfolioNameConvolutionLossGen::getIndexWeights() const
{
    return indexWeights;
}
    
/**
 * Mechanism used to pass the "single name" notional to
 * the tranche level
 * */
double BCPortfolioNameConvolutionLossGen::getNotional() const
{
    return ptfName->getNameNotional();
}

/**
 * Mechanism used to pass the "single name" historical beta to
 * the tranche level
 * */
double BCPortfolioNameConvolutionLossGen::getHistoricalBeta() const
{
    return ptfName->getBeta();
}

/**
 * Mechanism used to pass the "single name" (unconditional)
 * expected loss to the tranche level
 * */
double BCPortfolioNameConvolutionLossGen::getExpectedLoss(const DateTime& time) const
{
    double defaultProba = 1.0 - 
        defaultRates->calcDefaultPV(
            ptfName->getProtectionStartDate(), time);
                
    return defaultProba * ptfName->getNameNotional() * (1.0 - ptfName->getNameRecovery());
}

/** Mechanism used to pass the "default status" to the tranche level */
bool BCPortfolioNameConvolutionLossGen::hasDefaulted() const
{
    return ptfName->hasDefaulted();
}

/** Mechanism used to pass the loss config "single name" name (!) to the tranche level */
string BCPortfolioNameConvolutionLossGen::getName() const
{
    return ptfName->getName();
}

/**
 * Mechanism used to pass the "single name" implied par spreads and 
 * durations to the tranche level.
 * - Needed for "wide spread" with duration weighted spread ratio -
 * */
void BCPortfolioNameConvolutionLossGen::impliedParSpreadsAndDurations(
    const YieldCurveConstSP discount,
    const DateTimeArray& dates,
    DoubleArray& impliedSpreads,  /* (Output) */
    DoubleArray& durations) const /* (Output) */
{
    ptfName->impliedParSpreadsAndDurations(
        discount,
        dates,
        impliedSpreads,
        durations);
}

 /** [Implements IBCCondLossDistributionsGen] */
ICondLossDistributionsGenKeySP BCPortfolioNameConvolutionLossGen::initialise(
    const DateTime& time,
    double adjustedStrike,
    const IBCSkewCalculator& bcSkewCalculator) const
{
    static const string method("BCPortfolioNameConvolutionLossGen::initialise");
    
    if (condDefaultsModel == 0)
    {
        throw ModelException(method,
            "Unable to initialise thresholds "
            "(no conditional defaults model has been provided).");
    }
    
    // Retrieves loss level
    double lossLevel = ptfName->getNameNotional() * (1.0 - ptfName->getNameRecovery());
    
    const DateTime& protectionStartDate = ptfName->getProtectionStartDate();
    const DateTime& protectionEndDate = ptfName->getProtectionEndDate(time);
    double defaultProba = 1.0 - 
        defaultRates->calcDefaultPV(protectionStartDate, protectionEndDate);
    
    // Account for protectionEndDate (maturity cutoff)
    double adjustedBeta = bcSkewCalculator.computeBCadjustedBeta(
        ptfName->getBeta(), adjustedStrike, time.min(protectionEndDate));
    
    // Creates "Credit Metrics" parameters with the adjusted beta
    CreditEngineParametersSP betaCM(new CmOnlyParameters(
        "betaCM_" + Format::toString(adjustedStrike),
        adjustedBeta,   // beta
        CDoubleSP(0))); // no decretion beta

    return condDefaultsModel->initialise(
        defaultProba,
        defaultProba * lossLevel, // expected loss
        ptfName->getNameNotional(),
        betaCM,
        protectionStartDate,
        protectionEndDate);
}

#undef TINY

DRLIB_END_NAMESPACE
