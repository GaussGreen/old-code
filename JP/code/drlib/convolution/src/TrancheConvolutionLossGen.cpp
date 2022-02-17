#include "edginc/config.hpp"
#include "edginc/TrancheConvolutionLossGen.hpp"
#include "edginc/IModelConfigMapper.hpp"
#include "edginc/CDOPortfolio.hpp"
#include "edginc/DiscreteDistribution.hpp"
#include "edginc/DoubleArrayMarketFactor.hpp"
#include "edginc/TrancheCondELIntegrand.hpp"

DRLIB_BEGIN_NAMESPACE

TrancheConvolutionLossGen::~TrancheConvolutionLossGen(){}

/** TrancheConvolutionLossGen constructor */
TrancheConvolutionLossGen::TrancheConvolutionLossGen(
    const ICreditLossConfig* lossConfig,
    double lossLevel,
    const ICreditLossConfig* innerLossConfig,
    const IModelConfigMapper* mapper,
    const IConvolutor* convolutor,
    bool returnBinaryDistribution,
    const IConditionalDefaultsModel* condDefaultsModel,
    const DateTime& valueDate):
        lossConfig(lossConfig),
        lossLevel(lossLevel),
        valueDate(valueDate),
        innerLossConfig(innerLossConfig),
        mapper(mapper),
        convolutor(convolutor),
        returnBinaryDistribution(returnBinaryDistribution),
        condDefaultsModel(condDefaultsModel)
{
}
/** Implements "ILossDistributionsGen" */
IDistribution1DArraySP TrancheConvolutionLossGen::createLossDistributions(
    const DateTimeArray& timeline) const
{
    static const string method = "TrancheConvolutionLossGen::createLossDistributions";
    
    int i,t;
    
    // Use a reduced timeline if 1st point = value date
    // (in which case expected loss = 0)
    int timeLineOffset = (timeline[0] == valueDate ? 1 : 0);
    DateTimeArraySP reducedTimeLine(new DateTimeArray(
        timeline.begin() + timeLineOffset, timeline.end()));

    if (CDOPortfolio::TYPE->isInstance(innerLossConfig))
    {
        int nbNames = innerLossConfig->numInnerLossConfigs();
        int nbDates = reducedTimeLine->size();
        
        // Vector of size nbDates
        // For each date, contains an IKeyArraySP corresponding to names
        // in the portfolio
        vector<ICondLossDistributionsGenKeyArraySP> keysByDate(nbDates);
        for (t = 0; t < nbDates; ++t) {
			keysByDate[t].reset(new ICondLossDistributionsGenKeyArray(nbNames));
		}

        // Temporary variables
        ICreditLossGenConstSP lossGen;
        ICondLossDistributionsGenConstSP condLossGen;
        ICreditLossConfigConstSP name;
        ICreditLossModelConfigConstSP innerModel;
        ICondLossDistributionsGenKeyArrayConstSP keys;
        for (i = 0; i < nbNames; ++i) {
            name = innerLossConfig->getInnerLossConfig(i);
            
            // Retrieves relevant ICreditLossModelConfig
            innerModel = mapper->innerModel(name);
            
            // Retrieves loss generator
            lossGen = innerModel->lossGenerator(
                name, IModelConfigMapperConstSP(mapper));
            
            // Casts lossGen: here wants a ICondLossDistributionsGen
            condLossGen = DYNAMIC_POINTER_CAST<const ICondLossDistributionsGen>(lossGen);
            
            if (condLossGen.get() == 0)
            {
                throw ModelException(
                    method,
                    "Unable to cast asset " +
                    Format::toString(i) +
                    " into a ICondLossDistributionsGen.");
            }
            
            // Initialises loss gen (NB: potentially a time consuming operation)
            keys = condLossGen->initialise(*reducedTimeLine);
                
            // Populates keysByDate
            for (t = 0; t < nbDates; ++t) {
                (*keysByDate[t])[i] = (*keys)[t];
            }
        }

        if (returnBinaryDistribution)
        {
            int nbDates = reducedTimeLine->size();
            DoubleArray expectedLosses(nbDates);
            for (t = 0; t < nbDates; ++t) {
                TrancheCondELIntegrand condELFunction(
                        keysByDate[t],
                        condDefaultsModel->marketFactorDimension(),
                        timeline[t],
                        lossConfig,
                        convolutor);
                
                expectedLosses[t] = condDefaultsModel->integrateCondFunction(
                    &condELFunction, keysByDate[t], timeline[t]);
            }

            // Builds binary distributions for each point in the timeline
            IDistribution1DArraySP binaryDistr(new IDistribution1DArray(timeline.size()));

            if (timeline[0] == valueDate)
            {
                (*binaryDistr)[0].reset(
                    new DiscreteDistribution(lossLevel, 0));
                for (t = 1; t < timeline.size(); ++t) {
                    (*binaryDistr)[t].reset(
                        new DiscreteDistribution(lossLevel, expectedLosses[t-1] / lossLevel));
                }
            }
            else
            {
                for (t = 0; t < timeline.size(); ++t) {
                    (*binaryDistr)[t].reset(
                        new DiscreteDistribution(lossLevel, expectedLosses[t] / lossLevel));
                }
            }
            
            return binaryDistr;
        }
        else
        {
            //TODO
            throw ModelException(
                method,
                "Not implemented (returning non binary distributions)");
        }
    }
    else
    {
        // In that case inner loss config is not a portfolio.
        // Probably just wants to apply the lossConfig payoff to
        // the inner loss config distribution...
        // TODO
        throw ModelException(
            method,
            "Not implemented (inner loss config is not a CDOPortfolio)");
    }
}

/** [Implements ICondLossDistributionsGen] */
ICondLossDistributionsGenKeyArrayConstSP TrancheConvolutionLossGen::initialise(
    const DateTimeArray& timeline) const
{
    if (condDefaultsModel == 0)
    {
        throw ModelException(
            "TrancheConvolutionLossGen::initialise",
            "Unable to compute conditional loss distribution "
            "(no conditional defaults model has been provided).");
    }
    
    // Retrieves model parameters
    CreditEngineParametersConstSP modelParameters =
        innerLossConfig->getEngineParams(condDefaultsModel->engineParamsType());
    
    // Retrieves unconditional loss distributions for each point of the timeline
    // (usually very expensive !)
    IDistribution1DArraySP lossDistributions = createLossDistributions(timeline);

    int nbDates = timeline.size();

    // Initialises ICondLossDistributionsGenKey
    ICondLossDistributionsGenKeyArraySP keys(
        new ICondLossDistributionsGenKeyArray(nbDates));
    
    double expectedLoss;
    for (int t = 0; t < nbDates; ++t) {
        expectedLoss = (*lossDistributions)[t]->expectation();
        (*keys)[t] = condDefaultsModel->initialise(
            expectedLoss / lossLevel, // default proba
            expectedLoss,
            lossLevel, // notional
            modelParameters,
            valueDate,
            timeline[t]);
    }

    return keys;
}

DRLIB_END_NAMESPACE
