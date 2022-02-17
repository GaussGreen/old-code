#include "edginc/config.hpp"
#include "edginc/BCTrancheConvolutionLossGen.hpp"
#include "edginc/CDOPortfolio.hpp"
#include "edginc/DiscreteDistribution.hpp"
#include "edginc/IBCCondLossDistributionsGen.hpp"
#include "edginc/RecursiveConvolutor.hpp"
#include "edginc/BetaConvolutor.hpp"
#include "edginc/DoubleArrayMarketFactor.hpp"
#include "edginc/TrancheCondELIntegrand.hpp"

DRLIB_BEGIN_NAMESPACE

#define TINY 3e-15

/** IBCSkewCalculator for BCTrancheConvolutionLossGen */
class BCSkewCalculator: public virtual IBCSkewCalculator
{
public:
    BCSkewCalculator(
        SkewSurface::SkewType skewType,
        const BCTrancheConvolutionLossGen::IndexDataMap& indexDataMap,
        double histAvgBeta,
        double compressionRatio,
        double historicalBetaAdjustment):
            skewType(skewType),
            indexDataMap(indexDataMap),
            histAvgBeta(histAvgBeta),
            compressionRatio(compressionRatio),
            historicalBetaAdjustment(historicalBetaAdjustment) {}

    /** [Implements IBCSkewCalculator] */
    virtual double computeBCadjustedBeta(
        double nameHistBeta, double adjStrike, const DateTime& time) const;

private:
    SkewSurface::SkewType skewType;
    const BCTrancheConvolutionLossGen::IndexDataMap& indexDataMap;
    double histAvgBeta;
    double compressionRatio;
    double historicalBetaAdjustment;
};

BCTrancheConvolutionLossGen::~BCTrancheConvolutionLossGen(){}

/** BCTrancheConvolutionLossGen constructor */
BCTrancheConvolutionLossGen::BCTrancheConvolutionLossGen(
    const ITranchesCombinationPayoff* tranchesCombinationPayoff,
     double lossLevel,
    const ICreditLossConfig* innerLossConfig,
    const DateTime& valueDate,
    const IModelConfigMapper* mapper,
    const IConvolutor* convolutor,
    bool returnBinaryDistribution,
    const IConditionalDefaultsModel* condDefaultsModel,
    double compressionRatio,
    CDoubleConstSP strikeMappingOverride,
    DoubleArrayConstSP lowerBCBetasOverride,
    DoubleArrayConstSP upperBCBetasOverride,
    DateTimeArrayConstSP bcBetasTimeline,
    CDoubleConstSP spreadRatioOverride,
    DoubleArrayConstSP spreadRatiosRegionalOverride,
    StringArrayConstSP spreadRatiosRegions,
    double historicalBetaAdjustment,
    StringArrayConstSP namesOverride,
    IndexWeightsArrayConstSP indexWeightsOverride,
    bool authoriseNegativeEL,
    bool useExpectedLossRatio,
    YieldCurveConstSP discountCurve):
        tranchesCombinationPayoff(tranchesCombinationPayoff),
        lossLevel(lossLevel),
        innerLossConfig(innerLossConfig),
        mapper(mapper),
        convolutor(convolutor),
        returnBinaryDistribution(returnBinaryDistribution),
        condDefaultsModel(condDefaultsModel),
        valueDate(valueDate),
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
        useExpectedLossRatio(useExpectedLossRatio),
        discountCurve(discountCurve){}

/** Intermediate structure for index related data (for a given timepoint) */
struct BCTrancheConvolutionLossGen::IndexData {
    // Weight of this index in the bespoke portfolio
    double weight; // NB - not timepoint dependent
    
    // Expected loss of this index
    double expectedLoss;
    
    // Expected loss of the single names of the bespoke portfolio 
    // assigned to this index
    double expectedLossBespoke;

    // Weighted sum of the notionals of the single names
    // assigned to this index
    double sumNotionalBespoke;
    
    // Skew associated to this index
    IndexSkewConstSP skew; // NB - not timepoint dependent
    
    // Index strike
    double indexStrike;
    
    // Spread ratio
    double spreadRatio;
    
    // Notional and duration weighted sum for names of the bespoke portfolio
    // assigned to this index
    double sumNotionalDurationBespoke;

    // Strike mapping (also called 'q')
    double strikeMapping; // NB - not timepoint dependent
    
    // Constructor
    IndexData(IndexSkewConstSP skew, double strikeMapping);
    
    // Default constructor (used by map<string, IndexData>[key])
    IndexData();
};

/** Constructor for IndexData */
BCTrancheConvolutionLossGen::IndexData::IndexData(
    IndexSkewConstSP skew, double strikeMapping) : 
    weight(0.0), 
    expectedLoss(0.0),
    expectedLossBespoke(0.0),
    sumNotionalBespoke(0.0), 
    skew(skew),
    indexStrike(0.0),
    spreadRatio(0.0),
    sumNotionalDurationBespoke(0.0),
    strikeMapping(strikeMapping) {}

/** Default constructor for IndexData (used by map<string, IndexData>[key]) */
BCTrancheConvolutionLossGen::IndexData::IndexData() {} 

/** Implements "ILossDistributionsGen" */
IDistribution1DArraySP BCTrancheConvolutionLossGen::createLossDistributions(
    const DateTimeArray& timeline) const
{
    static const string method = "BCTrancheConvolutionLossGen::createLossDistributions";
    
    int i,t;
    
    if (CDOPortfolio::TYPE->isInstance(innerLossConfig))
    {
        int nbNames = innerLossConfig->numInnerLossConfigs();
        
        IBCCondLossDistributionsGenArraySP lossGens(
            new IBCCondLossDistributionsGenArray(nbNames));
        
        for (i = 0; i < nbNames; ++i) {
            ICreditLossConfigConstSP name = innerLossConfig->getInnerLossConfig(i);
            
            // Retrieves relevant ICreditLossModelConfig
            ICreditLossModelConfigConstSP innerModel = mapper->innerModel(name);
            
            // Retrieves loss generator
            ICreditLossGenSP lossGen =
                innerModel->lossGenerator(name, IModelConfigMapperConstSP(mapper));
            
            // Casts lossGen: here wants a IBCCondLossDistributionsGen
            (*lossGens)[i] =
                DYNAMIC_POINTER_CAST<IBCCondLossDistributionsGen>(lossGen);

            if ((*lossGens)[i].get() == 0)
            {
                throw ModelException(
                    method,
                    "Unable to cast asset " +
                    Format::toString(i) +
                    " into a IBCCondLossDistributionsGen.");
            }
        }

        if (returnBinaryDistribution)
        {
            // Pass a reduced timeline to avoid calculation for 1st point (if relevant)
            int timeLineOffset = (timeline[0] == valueDate ? 1 : 0);
            DateTimeArraySP reducedTimeLine(new DateTimeArray(
                timeline.begin() + timeLineOffset, timeline.end()));
            int nbDates = reducedTimeLine->size();

            // Choose skew type (fast / normal)
            SkewSurface::SkewType skewType;
            if (RecursiveConvolutor::TYPE->isInstance(convolutor))
            {
                skewType = SkewSurface::NORMAL;
            }
            else if (BetaConvolutor::TYPE->isInstance(convolutor))
            {
                skewType = SkewSurface::FAST;
            }
            else
            {
                throw ModelException(
                    "BCTrancheConvolutionLossGen::computeBCSkews",
                    "No skew surface corresponds to convolutor of type" +
                    convolutor->getClass()->getName() + ".");
            }
            
            IndexDataMap indexDataMap;
            double pastPortLoss =
                dynamic_cast<const CDOPortfolio*>(innerLossConfig)->portfolioLongRecoveredNotional();

            // Computes sumOutNotional and histAvgBeta
            double histAvgBeta, sumOutPosNotional, sumOutNegNotional, sumOutNetNotional;
            computeSumNotionalAndBeta(
                lossGens, 
                sumOutPosNotional, 
                sumOutNegNotional, 
                sumOutNetNotional, 
                histAvgBeta);
        
            // initialise expected losses
            DoubleArray expectedLosses(nbDates, 0.0);
            
            for (t = 0; t < nbDates; ++t) {

                // initialise "payoff" for timepoint t
                DoubleArray baseStrikes;
                DoubleArray baseStrikesWeights;
                double expectedLossOffset;
                tranchesCombinationPayoff->linearDecomposition(
                    (*reducedTimeLine)[t],
                    baseStrikes,
                    baseStrikesWeights,
                    expectedLossOffset);

                // add expected loss offset
                expectedLosses[t] += expectedLossOffset;

                initIndexDataMap(
                    valueDate,
                    lossGens,
                    discountCurve,
                    (*reducedTimeLine)[t],
                    pastPortLoss,
                    indexDataMap);
                    
                BCSkewCalculator bcSkewCalculator(
                    skewType, indexDataMap, histAvgBeta, compressionRatio, historicalBetaAdjustment);
                    
                // loop over baseStrikes
                for (int k = 0; k < baseStrikes.size(); ++k) {

                    // Adjustment of strike and conversion to percentage
                    double adjStrike = (baseStrikes[k] - pastPortLoss ) / (sumOutPosNotional - sumOutNegNotional);

                    ICondLossDistributionsGenKeyArraySP keys(
                        new ICondLossDistributionsGenKeyArray(nbNames));
                    for (i = 0; i < nbNames; ++i) {
                        // Initialises loss gen (NB: potentially a time consuming operation)
                        (*keys)[i] = (*lossGens)[i]->initialise(
                            (*reducedTimeLine)[t],
                            adjStrike,
                            bcSkewCalculator);
                    }
                    
                    CreditTrancheLossConfigSP tranche(
                        new CreditTrancheLossConfig(0.0, baseStrikes[k]));
        
                    TrancheCondELIntegrand condELFunction(
                            keys,
                            condDefaultsModel->marketFactorDimension(),
                            (*reducedTimeLine)[t],
                            tranche.get(),
                            convolutor);
                        
                    expectedLosses[t] += baseStrikesWeights[k] *
                        condDefaultsModel->integrateCondFunction(&condELFunction, keys, (*reducedTimeLine)[t]);
                }

                if (!authoriseNegativeEL) {
                    expectedLosses[t] = Maths::max(expectedLosses[t], 0.0);
                }
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
        // Probably just wants to apply tranche payoffs to
        // the inner loss config distribution...
        // TODO
        throw ModelException(
            method,
            "Not implemented (inner loss config is not a CDOPortfolio)");
    }
}

/** [Implements ICondLossDistributionsGen] */
ICondLossDistributionsGenKeyArrayConstSP BCTrancheConvolutionLossGen::initialise(
    const DateTimeArray& timeline) const
{
    if (condDefaultsModel == 0)
    {
        throw ModelException(
            "BCTrancheConvolutionLossGen::initialise",
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

/** Initialises indexDataMap and computes sumOutNotional */
//TODO: refactor this method - overly long and complex
void BCTrancheConvolutionLossGen::initIndexDataMap(
    const DateTime& today,
    IBCCondLossDistributionsGenArrayConstSP lossGens,
    YieldCurveConstSP discount,
    const DateTime& timepoint,
    double pastPortLoss,
    IndexDataMap& indexDataMap) const
{
    indexDataMap.clear();
    
    // Name of this method
    static const string method(
        "BCTrancheConvolutionLossGen::initIndexDataMap");
    try {
        IndexWeightsMap indexWeightsMap;
        
        // Builds a (portfolio name, index weights) map with overridden names
        IndexWeightsMap indexWeightsOverrideMap;
        if (namesOverride.get()) {
            for(int i = 0; i < namesOverride->size(); i++) {
                indexWeightsOverrideMap[(*namesOverride)[i]] =
                    (*indexWeightsOverride)[i];
            }    
        }
        int numNames = lossGens->size();
        // Builds the (portfolio name, index weights) map
        int i;
        for (i=0; i < numNames; i++) {
            if (!(*lossGens)[i]->hasDefaulted()){ // inner "name" has not defaulted
                const string& name = (*lossGens)[i]->getName();
                IndexWeightsMap::const_iterator iter = 
                    indexWeightsOverrideMap.find(name); 
                if (iter != indexWeightsOverrideMap.end()) {
                    // use the override
                    indexWeightsMap[name] = iter->second;
                } else {
                    // use the market value
                    indexWeightsMap[name] = (*lossGens)[i]->getIndexWeights();
                }
            }
        }
        
        // Some pre process for spread ratios computation
        // Note that this is NOT used when spreadRatioOverride != NULL
        // We still do it here to keep the code readable
        DoubleArray impliedSpreads(1);
        DoubleArray durations(1);
        DateTimeArray timepointAsArray(1, timepoint);
        
        // Initialise "name to spread ratio override" local map
        map<string, double> nameToSpreadRatioOverride;
        if (spreadRatiosRegionalOverride.get() != 0) {
            for (i = 0; i < spreadRatiosRegionalOverride->size(); ++i) {
                nameToSpreadRatioOverride[(*spreadRatiosRegions)[i]] = (*spreadRatiosRegionalOverride)[i];
            }
        }
        
        // Builds the (index name, IndexData) map
        for(i = 0; i< numNames; i++) {
            double nameNotional    = (*lossGens)[i]->getNotional();
            
            // flag to check that duration calculation is needed (true = needed)
            bool needsDurations = false;
            if (!(*lossGens)[i]->hasDefaulted()){ // inner "name" has not defaulted
                const string& name = (*lossGens)[i]->getName();
                IndexWeightsConstSP indexWeights = indexWeightsMap[name];
                IndexSkewWrapperArraySP indexes(indexWeights->getIndexNames());
                CDoubleArraySP weights = indexWeights->getWeights();
                
                if (spreadRatioOverride.get() == NULL) {
                    // calculate implied spreads and durations for this name
                    
                    // checks that duration calculation is really needed
                    for(int j = 0; j < indexes->size(); j++) {
                        needsDurations =
                            needsDurations || ((*indexes)[j]->useOffSpreads());
                    }
                    // if we use expected loss ratio then we need neither the durations
                    // nor the par spreads
                    needsDurations = needsDurations && !useExpectedLossRatio;
                
                    if (needsDurations) {
                        // NB: supports forward starting names
                        (*lossGens)[i]->impliedParSpreadsAndDurations(
                            discount,
                            timepointAsArray,
                            impliedSpreads,
                            durations);
                    }
                }
                
                for(int j = 0; j < indexes->size(); j++) {
                    const string& indexName = (*indexes)[j].getName();
                    if (indexDataMap.find(indexName) == indexDataMap.end()) {
                        // indexName not found in the map : initialise
                        double strikeMapping =
                            strikeMappingOverride.get() == NULL?
                            (*indexes)[j]->getStrikeMapping():
                            strikeMappingOverride->doubleValue();
                        
                        indexDataMap[indexName].reset(
                            new IndexData((*indexes)[j].getSP(), strikeMapping));
                    }
                    IndexDataSP indexData = indexDataMap[indexName];
                    double weight = (*weights)[j];
                    
                    if (spreadRatioOverride.get() == NULL && // no global override
                        nameToSpreadRatioOverride.find(indexName) == nameToSpreadRatioOverride.end() && // no local override
                        needsDurations)
                    {
                        // update spread ratio
                        double coefSR = weight * nameNotional * durations[0];
                        indexData->sumNotionalDurationBespoke += coefSR;
                        indexData->spreadRatio += coefSR * impliedSpreads[0];
                    }
                    
                    // Update indexData
                    indexData->sumNotionalBespoke += weight * nameNotional;
                        indexData->expectedLossBespoke +=
                            weight * (*lossGens)[i]->getExpectedLoss(timepoint);
                }            
            }
        }


        // we need to compute the sum of absolute value of each weight assigned 
        // to each index. The final weight will be the absolute value of the ones computed above
        // divided by the sum of abs of these weight.
        // This defines the BC weigts to use with shorts :
        // 1) if a portfolio is 100 long US, 10 long Europe, 10 short Europe, BC is 100% US
        // 2) if a portfolio is 110 long US, 10 short Europe, BC is 100/110 US + 10/110 Europe 
        double sumSumAbsNotionalBespoke = 0.;
        IndexDataMap::iterator mapIter = indexDataMap.begin();
        for (; mapIter != indexDataMap.end(); mapIter++){
            IndexDataSP indexData = (*mapIter).second;
            sumSumAbsNotionalBespoke += fabs(indexData->sumNotionalBespoke);
        }


        // Normalises and populate index expected loss
        mapIter = indexDataMap.begin();
        for (; mapIter != indexDataMap.end(); mapIter++){
            IndexDataSP indexData = (*mapIter).second;
            
            // Normalises expectedLossesBespoke
             indexData->weight = fabs(indexData->sumNotionalBespoke) / sumSumAbsNotionalBespoke;
            if (fabs(indexData->sumNotionalBespoke) < TINY) {
                indexData->expectedLossBespoke = 0.0;
            } else {
                indexData->expectedLossBespoke /= indexData->sumNotionalBespoke;
            }

            // Populates index expected loss
            // NB : this is not time consuming because clean spreads are cached
            DefaultRatesSP indexDefaultRates(indexData->skew->defaultRates());
            double indexRecovery = indexData->skew->getRecovery();
            
            indexData->expectedLoss =
                (1.0 - indexRecovery) *
                (1.0 - indexDefaultRates->calcDefaultPV(today,timepoint));
            
            // Update spread ratios
            map<string, double>::iterator nameToSpreadRatioOverrideIter = nameToSpreadRatioOverride.find((*mapIter).first);
            if (spreadRatioOverride.get() == NULL && // no global override
                nameToSpreadRatioOverrideIter == nameToSpreadRatioOverride.end()) // no local override
            {
                if (indexData->skew->useOffSpreads())
                {
                    if (useExpectedLossRatio)
                    {
                        // compute ratio of expected losses
                        if (indexData->expectedLoss != 0.0)
                        {
                            indexData->spreadRatio = indexData->expectedLossBespoke / indexData->expectedLoss;
                        }
                        else
                        {
                            indexData->spreadRatio = 1.0;
                        }
                    }
                    else
                    {
                        // Computes index implied spreads
                        indexData->skew->impliedSpreadsAndDurations(discount,
                                                                    timepointAsArray,
                                                                    impliedSpreads,
                                                                    durations);
                        
                        // Normalises spread ratios
                        double normalisationFactor =
                            indexData->sumNotionalDurationBespoke * impliedSpreads[0];
                        
                        // hack: normalisationFactor(t) can be =0 when we
                        // use wide spreads with forward starting tranches and
                        // t < forward start date.
                        // In that case the survival probability of a 
                        // name conditional of a market factor is always = 1
                        // and hence does not depend on adjusted beta, a fortiori
                        // does not depend on spread ratio (so we can "arbitrarily" choose
                        // spread ratio=1). Even if spread ratio does not affect pricing,
                        // invalid values (=1/0) can lead to failure, that is why we set it to 1.
                        // A non hacky solution would be to use expected loss instead of spreads
                        // and durations and compute an "expected loss ratio" instead of
                        // a "spread ratio".
                        if (normalisationFactor != 0.0) {
                            indexData->spreadRatio /= normalisationFactor;
                        } else {
                            indexData->spreadRatio = 1.0;
                        }
                    }
                } else {
                    // indexData->skew does not use "off spreads" -> defaults
                    // spread ratio to 1.0
                    indexData->spreadRatio = 1.0;
                }
            } else {
                // Uses relevant spread ratio override
                double srOverride;
                if (spreadRatioOverride.get() != NULL)
                {
                    srOverride = spreadRatioOverride->doubleValue();
                }
                else
                {
                    srOverride = nameToSpreadRatioOverrideIter->second;
                }
                indexData->spreadRatio = srOverride;
            }
        }
    } catch (exception& e){
        throw ModelException(e, method);
    }
}

/** Auxiliary function to compute sumOutNotional and histAvgBeta */
void BCTrancheConvolutionLossGen::computeSumNotionalAndBeta(
    IBCCondLossDistributionsGenArrayConstSP lossGens,
    double& sumOutPosNotional,
    double& sumOutNegNotional,
    double& sumOutNetNotional,
    double& histAvgBeta) const
{
    static const string routine(
        "BCTrancheConvolutionLossGen::computeSumNotionalAndBeta");
    try {
        sumOutPosNotional = 0.0;
        sumOutNegNotional = 0.0;
        sumOutNetNotional = 0.0;
        histAvgBeta = 0.0;        
        double sumBetaWeights = 0.0;
        int numNames = lossGens->size();
        for(int i = 0; i < numNames; i++) {
            double nameNotional = (*lossGens)[i]->getNotional();
            double betaWeight = fabs(nameNotional);
            if (!(*lossGens)[i]->hasDefaulted()){ // inner "name" has not defaulted
                histAvgBeta += betaWeight * (*lossGens)[i]->getHistoricalBeta();
                sumBetaWeights += betaWeight;
                sumOutPosNotional += Maths::max( nameNotional, 0.0);
                sumOutNegNotional += Maths::max(-nameNotional, 0.0);
                sumOutNetNotional += nameNotional;
            }
        }
        
        // Checks sumOutPosNotional != 0.0
        if (fabs(sumOutPosNotional) < TINY) {
            throw ModelException(routine, "Sum of positive outstanding"
                                 " name notionals is 0.0");
        }

        // Checks sumBetaWeights != 0.0
        if (fabs(sumBetaWeights) < TINY) {
            throw ModelException(routine, "Sum of abs(outstanding "
                                 "name notionals) is 0.0");
        }
        histAvgBeta /= sumBetaWeights;
    } 
    catch (exception& e){
        throw ModelException(e, routine);
    }
}

double BCSkewCalculator::computeBCadjustedBeta(
    double nameHistBeta, double adjStrike, const DateTime& time) const
{
    // Name of this method
    static const string method =
        "BCSkewCalculator::computeBCSkew";

    try {
        // Local variables
        double q, expLoss, expLossBespoke, indexStrike;
        BCTrancheConvolutionLossGen::IndexDataSP indexData;
        
        double betaSkew = 0.0;

        // Caps and floors
        adjStrike = Maths::min(1.0, Maths::max(adjStrike, 0.0));
        
        // Computes the index strikes
        BCTrancheConvolutionLossGen::IndexDataMap::const_iterator mapIter = indexDataMap.begin();
        for (; mapIter != indexDataMap.end(); mapIter++){
            indexData = (*mapIter).second;
            q = indexData->strikeMapping;
            expLoss = indexData->expectedLoss;
            expLossBespoke = indexData->expectedLossBespoke;
            indexStrike = adjStrike;
            if (fabs(expLossBespoke) >= TINY) {
                indexStrike *= (q + (1.0 - q) * expLoss / expLossBespoke);
            }
            // Caps and floors
            indexData->indexStrike =
                Maths::min(1.0, Maths::max(indexStrike, 0.0));
        }
        
        // Computes Beta Implied Skew
        mapIter = indexDataMap.begin();
        for (; mapIter != indexDataMap.end(); mapIter++){
            indexData = (*mapIter).second;
            
            betaSkew += 
                indexData->weight
                * indexData->skew->getSkew(
                    indexData->indexStrike, time, indexData->spreadRatio, skewType);
        }
        
        // Computes the squeeze tweak applied to the Beta skew
        double betaTwkSqueeze = 0.0;
        mapIter = indexDataMap.begin();
        for (; mapIter != indexDataMap.end(); mapIter++){
            indexData = (*mapIter).second;
            betaTwkSqueeze += indexData->weight * 
                indexData->skew->getSqueeze()->map(
                    histAvgBeta
                    - indexData->skew->getHistoricalBeta()
                    - historicalBetaAdjustment);
        }
        
        // Applies the squeeze
        double betaSkewSqueezed = PortfolioName::betaTweak(betaSkew, betaTwkSqueeze);
        
        // Dispersion adjust beta
        double dispersionAdjBeta =
            histAvgBeta + compressionRatio * (nameHistBeta - histAvgBeta);
        if (dispersionAdjBeta > 1.0 - TINY || dispersionAdjBeta < -1.0 + TINY) {
            throw ModelException(method,
                "Dispersion adjustment leads to"
                " beta outside [-1,1]");
        }
        
        // Adjust for historical average beta
        double betaTwk = (betaSkewSqueezed - histAvgBeta) / (1.0 - histAvgBeta);
        
        // Apply tweak
        return PortfolioName::betaTweak(dispersionAdjBeta, betaTwk);
        
    } catch (exception& e){
        throw ModelException(e, method);
    }
}

#undef TINY

DRLIB_END_NAMESPACE
