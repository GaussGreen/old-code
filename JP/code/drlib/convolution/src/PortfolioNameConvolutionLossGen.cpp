#include "edginc/config.hpp"
#include "edginc/PortfolioNameConvolutionLossGen.hpp"
#include "edginc/SingleCreditAsset.hpp"
#include "edginc/DiscreteDistribution.hpp"

DRLIB_BEGIN_NAMESPACE

PortfolioNameConvolutionLossGen::~PortfolioNameConvolutionLossGen() {}

PortfolioNameConvolutionLossGen::PortfolioNameConvolutionLossGen(
    const PortfolioName* ptfName,
    const IConditionalDefaultsModel* condDefaultsModel):
        ptfName(ptfName),
        condDefaultsModel(condDefaultsModel){}

/** Implements "ILossDistributionsGen" */
IDistribution1DArraySP PortfolioNameConvolutionLossGen::createLossDistributions(
    const DateTimeArray& timeline) const
{
    IDistribution1DArraySP lossDistributions(
        new IDistribution1DArray(timeline.size()));
    
    double lossLevel = ptfName->getNameNotional() * (1.0 - ptfName->getNameRecovery());

    // Retrieves underlying "single name" CDS
    const ICDSParSpreads* cds = SingleCreditAssetConstSP::dynamicCast(
        ptfName->getCreditAsset())->cdsParSpreads.get();
    
    // Retrieves clean spreads
    DefaultRatesSP defaultRates = cds->defaultRates();

    double defaultProba;
    for (int i = 0; i < timeline.size(); ++i) {
        defaultProba = 1.0 - 
            defaultRates->calcDefaultPV(
                ptfName->getProtectionStartDate(), timeline[i]);
        
        (*lossDistributions)[i].reset(new DiscreteDistribution(lossLevel, defaultProba));
    }

    return lossDistributions;
}
   
/** [Implements ICondLossDistributionsGen] */
ICondLossDistributionsGenKeyArrayConstSP PortfolioNameConvolutionLossGen::initialise(
    const DateTimeArray& timeline) const
{
    if (condDefaultsModel == 0)
    {
        throw ModelException(
            "PortfolioNameConvolutionLossGen::initialise",
            "Unable to compute conditional loss distribution "
            "(no conditional defaults model has been provided).");
    }
    
    // Retrieves model parameters
    CreditEngineParametersConstSP modelParameters =
        ptfName->getEngineParams(condDefaultsModel->engineParamsType());

    // Retrieves loss level
    double lossLevel = ptfName->getNameNotional() * (1.0 - ptfName->getNameRecovery());
    
    // Retrieves underlying "single name" CDS
    const ICDSParSpreads* cds = SingleCreditAssetConstSP::dynamicCast(
        ptfName->getCreditAsset())->cdsParSpreads.get();
    
    // Retrieves clean spreads
    DefaultRatesSP defaultRates = cds->defaultRates();

    int nbDates = timeline.size();

    // Initialises ICondLossDistributionsGenKey
    ICondLossDistributionsGenKeyArraySP keys(
        new ICondLossDistributionsGenKeyArray(nbDates));
    
    const DateTime& protectionStartDate = ptfName->getProtectionStartDate();
    
    double defaultProba;
    for (int t = 0; t < nbDates; ++t) {
        const DateTime& protectionEndDate = ptfName->getProtectionEndDate(timeline[t]);
        defaultProba = 0.0;
        if (protectionStartDate < protectionEndDate)
        {
            defaultProba = 1.0 - 
                defaultRates->calcDefaultPV(protectionStartDate, protectionEndDate);
        }

        (*keys)[t] = condDefaultsModel->initialise(
            defaultProba,
            defaultProba * lossLevel, // expected loss
            ptfName->getNameNotional(),
            modelParameters,
            protectionStartDate,
            protectionEndDate);
    }

    return keys;
}

DRLIB_END_NAMESPACE
