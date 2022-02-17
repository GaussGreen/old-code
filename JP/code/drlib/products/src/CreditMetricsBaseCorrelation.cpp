//----------------------------------------------------------------------------
//
//   Group       : Quantitative Research & Development
//
//   Filename    : CreditMetricsBaseCorrelation.cpp
//
//   Description : Credit Metrics model with Base Correlation
//
//   Author      : Antoine Gregoire
//
//   Date        : May 2005
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/AtomicArray.hpp"
#include "edginc/CreditAsset.hpp"
#include "edginc/CreditMetricsBaseCorrelation.hpp"
#include "edginc/SingleCreditAsset.hpp"
#include "edginc/IInstrumentCollection.hpp"
#include "edginc/CDSPricer.hpp"
#include "edginc/BaseCorrelationLossCalculator.hpp"
#include "edginc/LinearInterpolator.hpp"
#include "edginc/CreditMetricsLossCalculatorBase.hpp"
#include "edginc/PortfolioName.hpp"
#include "edginc/IndexSkew.hpp"
#include "edginc/Results.hpp"
#include "edginc/Control.hpp"
#include "edginc/BaseCorrelationOnlyParameters.hpp"
#include "edginc/CreditTrancheLossConfig.hpp"
#include "edginc/CounterPartyCredit.hpp"

DRLIB_BEGIN_NAMESPACE

#define TINY 3e-15

/** Useful function to check range of a field */
void inline checkRange(string fieldName, 
                       double fieldValue, 
                       double minValue, 
                       double maxValue) 
{
    if (fieldValue < minValue || fieldValue > maxValue) {
        throw ModelException(fieldName +
                             " (=" +
                             Format::toString(fieldValue) +
                             ") outside [" +
                             Format::toString(minValue) +
                             "," +
                             Format::toString(maxValue) +
                             "]");
    }
}

/** Intermediate structure for index related datas. */
struct CreditMetricsBaseCorrelation::IndexData {
    // Weight of this index in the bespoke portfolio
    double weight;
    
    // Expected loss of this index across the timeline
    DoubleArraySP expectedLosses;
    
    // Expected loss of the single names of the bespoke portfolio 
    // assigned to this index
    DoubleArraySP expectedLossesBespoke;

    // Weighted sum of the notionals of the single names
    // assigned to this index
    double sumNotionalBespoke;
    
    // Skew associated to this index
    IndexSkewConstSP skew;
    
    // Index strikes
    DoubleArraySP indexStrikes;
    
    // Spread ratios
    DoubleArraySP spreadRatios;
    
    // Notional and duration weighted sum for names of the bespoke portfolio
    // assigned to this index
    DoubleArraySP sumNotionalDurationBespoke;

    // Strike mapping (also called 'q')
    double strikeMapping;
    
    // Constructor
    IndexData(int size, IndexSkewConstSP skew, double strikeMapping);

    // Constructor
    IndexData(int size);
    
    // Default constructor (used by map<string, IndexData>[key])
    IndexData();
};

/** Constructor for IndexData */
CreditMetricsBaseCorrelation::IndexData::IndexData(
        int size, IndexSkewConstSP skew, double strikeMapping) : 
    weight(0.0), 
    expectedLosses(new DoubleArray(size, 0.0)),
    expectedLossesBespoke(new DoubleArray(size, 0.0)),
    sumNotionalBespoke(0.0), 
    skew(skew),
    indexStrikes(new DoubleArray(size, 0.0)),
    spreadRatios(new DoubleArray(size, 0.0)),
    sumNotionalDurationBespoke(new DoubleArray(size)),
    strikeMapping(strikeMapping) 
{}

/** Constructor for IndexData, with no skew or strike mapping */
CreditMetricsBaseCorrelation::IndexData::IndexData(int size) : 
    weight(0.0), 
    expectedLosses(new DoubleArray(size, 0.0)),
    expectedLossesBespoke(new DoubleArray(size, 0.0)),
    sumNotionalBespoke(0.0), 
    indexStrikes(new DoubleArray(size, 0.0)),
    spreadRatios(new DoubleArray(size, 0.0)),
    sumNotionalDurationBespoke(new DoubleArray(size)),
    strikeMapping(0.0) 
{}

/** Default constructor for IndexData (used by map<string, IndexData>[key]) */
CreditMetricsBaseCorrelation::IndexData::IndexData() 
{} 

/** Auxiliary function to compute sumOutNotional and histAvgBeta */
void computeSumNotionalAndBeta2(CreditTrancheLossConfigConstSP tranche,
                                double&                        sumOutPosNotional,
								double&					       sumOutNegNotional,
								double&				           sumOutNetNotional,
                                double&                        histAvgBeta)
{
    static const string routine("computeSumNotionalAndBeta2");
    try {
        sumOutPosNotional = 0.0;
		sumOutNegNotional = 0.0;
		sumOutNetNotional = 0.0;
        histAvgBeta = 0.0;        
        double sumBetaWeights = 0.0;
        int numNames = tranche->numInnerLossConfigs();
        for(int i = 0; i < numNames; i++) {
            double nameNotional = tranche->nameNotional(i);
            double betaWeight = fabs(nameNotional);
            if (!tranche->nameDefaulted(i)){ // name has not defaulted
                histAvgBeta += betaWeight * tranche->nameBeta(i);
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


/** Adjusts the dispersion of the single name betas around the average beta */
DoubleArraySP CreditMetricsBaseCorrelation::dispersionAdjust2(
    CreditTrancheLossConfigConstSP tranche,
    double                         avgBeta) const
{
    static const string method(
        "CreditMetricsBaseCorrelation::dispersionAdjust2");

    try {
        // reduction of historical betas dispersion
        int nbNames = tranche->numInnerLossConfigs();
        DoubleArraySP adjustedBetas(new DoubleArray(nbNames));
        for (int i = 0; i < nbNames; i++) {
            // adjust single name beta
            double beta = tranche->nameBeta(i);
            double newBeta = avgBeta + compressionRatio * (beta - avgBeta);
            if (newBeta > 1.0 - TINY || newBeta < -1.0 + TINY) {
                throw ModelException(method, "Dispersion adjustment leads to"
                                     " beta outside [-1,1] for name" +
                                     tranche->getInnerLossConfig(i)->getName());
            }
            (*adjustedBetas)[i] = newBeta; // save
        }
        return adjustedBetas;
    } 
    catch (exception& e){
        throw ModelException(e, method);
    }
}

/** Creates betas for each time point, for each name, using Base Correlation
    methodology. Betas are accessed as [time point][name] */
DoubleArrayArraySP CreditMetricsBaseCorrelation::calculateBCBetas(
    const DateTimeArray&                timeline,        // includes today
    CreditTrancheLossConfigConstSP      tranche,         // view onto Inst
    OutputRequest*                      betaRequest,     // may be null
    Results*                            results,         // for output result
    double                              strike,          // eg lowStrike
    DoubleArrayConstSP                  bcBetasOverride, // optional override
    SkewSurface::SkewType               skewType,        // fast or full
    double                              pastPortLoss,    // portfolio losses
    double                              sumOutPosNotional,
	double								sumOutNegNotional,
    double                              histAvgBeta,
    double                              avgBeta,         // average beta
    DoubleArrayConstSP                  adjBetas,        // dispersionAdjust
    const DoubleArray&                  timelineAsDouble,// dates->doubles
    const IndexDataMap&                 indexDataMap) const  // from initMaps
{
    try {
        int nbDates = timeline.size(); // for ease
        // Calculate the bespoke portfolio skew information for given strike
        DoubleArraySP bcBetas(new DoubleArray(nbDates));
        ccmBespokeSkewOneStrike(timeline,
                                pastPortLoss,
                                strike,
                                skewType,
                                indexDataMap,
                                sumOutPosNotional,
								sumOutNegNotional,
                                histAvgBeta,
                                *bcBetas);
        int t;
        for (t = 1; t < nbDates; ++t) {
            // "idx" is such that:
            // a) idx >= 0 and "timeline[t] = bcBetasTimeline[idx]
            // or
            // b) "idx" = -1 (means "no override" is used)
            int idx = -1;
            if (bcBetasOverride.get()) {
                try {
                    idx = timeline[t].find(*bcBetasTimeline);
                } 
                catch (exception&) {
                    // "idx" not found - this is weak. To do: fix
                    idx = -1;
                }
            }
            if (idx != -1 && bcBetasOverride.get()) {
                (*bcBetas)[t] = (*bcBetasOverride)[idx]; // use override
            }
        }        
        
        // if an output request has been supplied then store these betas 
        if (betaRequest){
            results->storeRequestResult(betaRequest, bcBetas);
        }
        LinearInterpolator interpolator;
        Interpolator::InterpolantConstSP betasInterpolator( 
            interpolator.computeInterp(timelineAsDouble, *bcBetas));
        
        // need to compute adjusted Beta for every name at every timepoint
        int numNames = tranche->numInnerLossConfigs();
        // betas accessed as [timepoint][name index]
        DoubleArrayArraySP betas(new DoubleArrayArray(nbDates));
        /* Loop through points of the timeline */
        for (t = 1; t < nbDates; ++t) {
            double betaTwk = ((*bcBetas)[t] - avgBeta) / (1.0 - avgBeta);
            // reserve some space
            (*betas)[t].resize(numNames);
            // then for each name
            for (int i = 0; i < numNames; ++i) {
                double betaTwkCutoff = betaTwk;
                const DateTime& date = timeline[t]; // for ease
                /* maturity cutoff adjustment */
                const DateTime& matCutOff = 
                    tranche->nameProtectionEndDate(i, date);
                if (matCutOff < date){
                    betaTwkCutoff =
                        (betasInterpolator->value(matCutOff.getDate()) 
                         - avgBeta) / (1.0 - avgBeta);
                }
                double adjustedBeta = (*adjBetas)[i]; // for ease
                /* Tweak the betas at lower strike */
                (*betas)[t][i] = PortfolioName::betaTweak(adjustedBeta, betaTwkCutoff);
            }
        }
        return betas;
    } 
    catch (exception& e){
        throw ModelException(e, "CreditMetricsBaseCorrelation::"
                             "calculateBCBetas");
    }
}

/** Generate the timeline on which the effective curve will be calculated.
    This overrides the one CreditMetricsModel in order to do some 
    additional validation */
DateTimeArraySP CreditMetricsBaseCorrelation::generateTimeline(
    const DateTime& today,
    const DateTime& lastObservationDate) const
{
    static const string method("CreditMetricsBaseCorrelation::"
                               "generateTimeline");

    /* 1. Check that lgdFloor = 0.0, lgdCap = 1.0, and
       lgdNotional = 1.0 for each name */
    /* This check is now performed in PortfolioName::validatePop2Object
       for all models, no need to duplicate */
//     int numNames = product->numNames();
//     for (int i = 0; i < numNames; i++) {
//         double lgdNotional, lgdFloor, lgdCap;
//         product->nameLGDRanges(i, lgdNotional, lgdFloor, lgdCap);
//         if (lgdNotional != 1.0 || lgdCap != 1.0 || lgdFloor != 0.0) {
//             SingleCreditAssetConstSP asset(product->nameAsset(i));
//             throw ModelException(method, "Name " + asset->getName() +
//                                  ": Base Correlation only supports "
//                                  "lgdNotional = 1, "
//                                  "lgdCap = 1 and lgdFloor = 0");
//         }
//     }

    // 2. generate timeline
    DateTimeArraySP timeline(
        CreditMetricsModel::generateTimeline(today, lastObservationDate));

    // 3. Checks consistency between bcBetasTimeline and timeline
    // Don't see how this works with theta for example
    if (bcBetasTimeline.get()) {
        // make sure they are increasing
        string fieldDesc("bcBetasTimeline in " + getClass()->getName());
        DateTime::ensureStrictlyIncreasing(*bcBetasTimeline, fieldDesc, false);

        // and that they are a subset of the timeline
        if (!DateTime::isSubset(*timeline, *bcBetasTimeline)) {
            throw ModelException(method, fieldDesc + " must be a subset of "
                                 "the effective curve timeline");
        }
    }
    return timeline;
}


/** Overridden to use "Base Correlation" methodology. Note that a
    IFixedTrancheLossCalculator is returned via the conditionalLossCalc 
    parameter. */
IFixedTrancheLossCalculator* 
CreditMetricsBaseCorrelation::createFixedLossCalculator(
    const DateTimeArray&                timeline,           /* (I) */
    CreditTrancheLossConfigConstSP      tranche,            /* (I) */
    CounterPartyCreditConstSP           cpty,               /* (I) */
    YieldCurveConstSP                   discount,           /* (I) */
    Control*                            control,            /* (I) */
    Results*                            results,            /* (I) */
    IFixedTrancheLossCalculatorConstSP& conditionalLossCalc) const // (O)
{
    // to do: refactor this method: too long and complex

    // create CreditMetrics style loss calculator
    CreditMetricsLossCalculatorBaseSP cmCalculator(
        createLossCalculatorBase(timeline, tranche, cpty));

    // Choose skew type (fast / normal)
    SkewSurface::SkewType skewType = useFastConvolution(tranche) ?
        SkewSurface::FAST : SkewSurface::NORMAL;
    // compute probabilities - to do reuse probs in CM loss calculator constctor
    DoubleArrayArray namesSurvivalProb;
    tranche->computeNameSurvivalProb(timeline, namesSurvivalProb);

    IndexWeightsMap indexWeightsMap;
    IndexDataMap    indexDataMap;
    double pastPortLoss = tranche->portfolioLoss();
    double sumOutPosNotional, sumOutNegNotional, sumOutNetNotional, histAvgBeta; // initialised by initMaps
    initMaps2(tranche, discount, timeline, namesSurvivalProb, pastPortLoss,
              indexWeightsMap, indexDataMap, sumOutPosNotional, 
              sumOutNegNotional, sumOutNetNotional, histAvgBeta);

    // calculate average beta
    double avgBeta = calculateAverageBeta2(tranche);
    // beta compression/dispersion
    DoubleArraySP adjBetas(dispersionAdjust2(tranche, avgBeta));

    int numDates = timeline.size(); // for ease
    /* A bit hacky - build a timeline made up of doubles for interpolation */
    DoubleArray timelineAsDouble(numDates);
    for (int t = 0; t < numDates; ++t) {
        timelineAsDouble[t] = timeline[t].getDate();
    }
    // look up the strikes
    double lowStrike, highStrike, baseStrike;
    tranche->getTrancheStrikes(lowStrike, highStrike);
	baseStrike =  - tranche->portfolioShortNotional();
    // then can calculate lower and upper betas
    OutputRequest* betaRequest = 
        control->requestsOutput(OutputRequest::TRANCHE_LOWER_BC_BETAS);
    DoubleArrayArraySP lowerBetas(
        calculateBCBetas(timeline, tranche, betaRequest, results,
                         lowStrike, lowerBCBetasOverride, // this line varies
                         skewType, pastPortLoss, sumOutPosNotional, 
                         sumOutNegNotional, histAvgBeta, avgBeta, 
                         adjBetas, timelineAsDouble, indexDataMap));

    betaRequest = 
        control->requestsOutput(OutputRequest::TRANCHE_UPPER_BC_BETAS);
    DoubleArrayArraySP upperBetas(
        calculateBCBetas(timeline, tranche, betaRequest, results,
                         highStrike, upperBCBetasOverride, // this line varies
                         skewType, pastPortLoss, sumOutPosNotional, 
                         sumOutNegNotional, histAvgBeta, avgBeta, 
                         adjBetas, timelineAsDouble, indexDataMap));
    // then for conditionalLossCalculator
    if (!!cpty){
        double portfolioLongNotional = tranche->portfolioLongNotional();
        double cccStrike = portfolioLongNotional * EQUITY_STRIKE_FOR_CCC;
        // get holds of the CCC betas
        DoubleArrayArraySP cccBetas(
            calculateBCBetas(timeline, tranche, NULL, results,
                             cccStrike, DoubleArraySP(), // this line varies
                             skewType, pastPortLoss, sumOutPosNotional, 
                             sumOutNegNotional, histAvgBeta, avgBeta, 
                             adjBetas, timelineAsDouble, indexDataMap));
        // now build a CM or CCM style loss calculator
        CreditMetricsLossCalculatorBaseSP 
            lossCalculator(createLossCalculatorBase(timeline, tranche, cpty));

        // then plug it into into a fixed strike one which allows the betas
        // to be overridden
        conditionalLossCalc.reset(
            CreditMetricsLossCalculatorBase::createFixedTrancheLossCalculator(
                lowStrike, highStrike, cccBetas, lossCalculator));
    }
    // now can finally build our loss calculator
    return new BaseCorrelationLossCalculator(cmCalculator, baseStrike,
                                             lowStrike, highStrike,
                                             lowerBetas, upperBetas,
                                             authoriseNegativeEL,
                                             useExpectedLossRatio);
}


/** Overridden to use "Base Correlation" methodology. Note that a
    IFixedTrancheLossCalculator is returned via the conditionalLossCalc 
    parameter. */
IFixedTrancheLossCalculator* 
CreditMetricsBaseCorrelation::createFixedRecoveredNotionalCalculator(
    const DateTimeArray&                timeline,           /* (I) */
    CreditTrancheLossConfigConstSP      tranche,            /* (I) */
    CounterPartyCreditConstSP           cpty,               /* (I) */
    YieldCurveConstSP                   discount,           /* (I) */
    Control*                            control,            /* (I) */
    Results*                            results,            /* (I) */
    IFixedTrancheLossCalculatorConstSP& conditionalLossCalc) const // (O)
{
    // to do: refactor this method: too long and complex

    // create CreditMetrics style loss calculator
    CreditMetricsLossCalculatorBaseSP cmCalculator(
        createRecoveredNotionalCalculatorBase(timeline, tranche, cpty));

    // Choose skew type (fast / normal)
    SkewSurface::SkewType skewType = useFastConvolution(tranche)?
        SkewSurface::FAST : SkewSurface::NORMAL;

    // compute probabilities - to do reuse probs in CM loss calculator constctor
    DoubleArrayArray namesSurvivalProb;
    tranche->computeNameSurvivalProb(timeline, namesSurvivalProb);

    IndexWeightsMap indexWeightsMap;
    IndexDataMap    indexDataMap;
    double pastPortLoss = tranche->portfolioLoss();
    double sumOutPosNotional, sumOutNegNotional, sumOutNetNotional, histAvgBeta; // initialised by initMaps
    initMaps2(tranche, discount, timeline, namesSurvivalProb, pastPortLoss,
              indexWeightsMap, indexDataMap, sumOutPosNotional, 
              sumOutNegNotional, sumOutNetNotional, histAvgBeta);

    // calculate average beta
    double avgBeta = calculateAverageBeta2(tranche);
    // beta compression/dispersion
    DoubleArraySP adjBetas(dispersionAdjust2(tranche, avgBeta));

    int numDates = timeline.size(); // for ease
    /* A bit hacky - build a timeline made up of doubles for interpolation */
    DoubleArray timelineAsDouble(numDates);
    for (int t = 0; t < numDates; ++t) {
        timelineAsDouble[t] = timeline[t].getDate();
    }
    // look up the strikes
    double lowStrike, highStrike;
    tranche->getTrancheStrikes(lowStrike, highStrike);
    // and the total notional
    double prtfNtnl = tranche->portfolioLongNotional();
    //and invert them
    double newLowStrike = prtfNtnl-highStrike;
    double newHighStrike = prtfNtnl-lowStrike;

	double baseStrike =  - tranche->portfolioShortNotional();

    // then can calculate lower and upper betas
    OutputRequest* betaRequest = 
        control->requestsOutput(OutputRequest::TRANCHE_LOWER_BC_BETAS);
    DoubleArrayArraySP lowerBetas(
        calculateBCBetas(timeline, tranche, betaRequest, results,
                         newLowStrike, lowerBCBetasOverride, // this line varies
                         skewType, pastPortLoss, sumOutPosNotional, 
                         sumOutNegNotional, histAvgBeta, avgBeta, 
                         adjBetas, timelineAsDouble, indexDataMap));

    betaRequest = 
        control->requestsOutput(OutputRequest::TRANCHE_UPPER_BC_BETAS);
    DoubleArrayArraySP upperBetas(
        calculateBCBetas(timeline, tranche, betaRequest, results,
                         newHighStrike, upperBCBetasOverride, // this line varies
                         skewType, pastPortLoss, sumOutPosNotional, 
                         sumOutNegNotional, histAvgBeta, avgBeta, 
                         adjBetas, timelineAsDouble, indexDataMap));
    // then for conditionalLossCalculator
    if (!cpty){
        double portfolioLongNotional = tranche->portfolioLongNotional();
        double cccStrike = portfolioLongNotional * EQUITY_STRIKE_FOR_CCC;
        // get holds of the CCC betas
        DoubleArrayArraySP cccBetas(
            calculateBCBetas(timeline, tranche, NULL, results,
                             cccStrike, DoubleArraySP(), // this line varies
                             skewType, pastPortLoss, sumOutPosNotional, 
                             sumOutNegNotional, histAvgBeta, avgBeta, 
                             adjBetas, timelineAsDouble, indexDataMap));
        // now build a CM or CCM style loss calculator
        CreditMetricsLossCalculatorBaseSP 
            lossCalculator(createRecoveredNotionalCalculatorBase(
                timeline, tranche, cpty));

        // then plug it into into a fixed strike one which allows the betas
        // to be overridden
        conditionalLossCalc.reset(
            CreditMetricsLossCalculatorBase::createFixedTrancheLossCalculator(
                newLowStrike, newHighStrike, cccBetas, lossCalculator));
    }
    // now can finally build our loss calculator
    return new BaseCorrelationLossCalculator(cmCalculator, baseStrike,
                                             newLowStrike, newHighStrike,
                                             lowerBetas, upperBetas,
                                             authoriseNegativeEL,
                                             useExpectedLossRatio);
}


/** Returns the indexWeightsMap for all names in the tranche, using
    the override if present */
void CreditMetricsBaseCorrelation::getIndexWeightsMap(
    CreditTrancheLossConfigConstSP tranche,               // (I)
    IndexWeightsMap&               indexWeightsMap) const // (O)
{
    // Builds a (portfolio name, index weights) map with overridden names
    IndexWeightsMap indexWeightsOverrideMap;
    if (namesOverride.get()) {
        for(int i=0; i < namesOverride->size(); ++i) {
            indexWeightsOverrideMap[(*namesOverride)[i]] = (*indexWeightsOverride)[i];
        }    
    }

    int numNames = tranche->numInnerLossConfigs();
    // Builds the (portfolio name, index weights) map
    for (int i=0; i < numNames; ++i) {
        if (!tranche->nameDefaulted(i)){ // name has not defaulted
            ICreditLossConfigConstSP lossCfg(tranche->getInnerLossConfig(i));
            const string& name = lossCfg->getName();
            IndexWeightsMap::const_iterator iter = 
                indexWeightsOverrideMap.find(name); 
            if (iter != indexWeightsOverrideMap.end()) { // override present, use it
                indexWeightsMap[name] = iter->second;
            } 
            else { // use the market value
                CreditEngineParametersConstSP params(lossCfg->getEngineParams(
                    BaseCorrelationOnlyParameters::TYPE));
                BaseCorrelationOnlyParametersConstSP bcParam =
                    BaseCorrelationOnlyParametersConstSP::dynamicCast(params); 

                indexWeightsMap[name] = bcParam->getIndexWeights();
            }
        }
    }
    return;
}



/** Initialises indexWeightsMap, indexDataMap and computes sumOutNotional.
    to do: refactor this method - overly long and complex */
void CreditMetricsBaseCorrelation::initMaps2(
    CreditTrancheLossConfigConstSP tranche,
    YieldCurveConstSP              discount,
    const DateTimeArray&           timeLine,
    const DoubleArrayArray&        namesSurvivalProb,
    double                         pastPortLoss,
    IndexWeightsMap&               indexWeightsMap,
    IndexDataMap&                  indexDataMap,
    double&                        sumOutPosNotional,
	double&                        sumOutNegNotional,
	double&                        sumOutNetNotional,
    double&                        histAvgBeta) const
{
    // Name of this method
    static const string method("CreditMetricsBaseCorrelation::initMaps2");
    try {
        // Computes sumOutNotional and histAvgBeta
        computeSumNotionalAndBeta2(tranche, 
                                   sumOutPosNotional, 
                                   sumOutNegNotional, 
                                   sumOutNetNotional, 
                                   histAvgBeta);

        // Obtain the indexWeightsMap for all names in the tranche
        getIndexWeightsMap(tranche,          // (I)
                           indexWeightsMap); // (O)
        
        int numNames = tranche->numInnerLossConfigs();
        // Builds the (index name, IndexData) map
        for(int i = 0; i< numNames; i++) {
			double nameNotional    = tranche->nameNotional(i);
            double nameRecovery    = tranche->nameRecovery(i);
            
            // flag to check that duration calculation is needed (true = needed)
            if (!tranche->nameDefaulted(i)){ // name has not defaulted
                SingleCreditAssetConstSP myAsset = tranche->nameAsset(i);
                const string& name = myAsset->getName();
                IndexWeightsConstSP indexWeights = indexWeightsMap[name];
                IndexSkewWrapperArraySP indexes(indexWeights->getIndexNames());
                CDoubleArraySP weights = indexWeights->getWeights();
                               
                for(int j = 0; j < indexes->size(); j++) {
                    const string& indexName = (*indexes)[j].getName();
                    if (indexDataMap.find(indexName) == indexDataMap.end()) {
                        // indexName not found in the map : initialise
                        double strikeMapping =
                            strikeMappingOverride == NULL?
                            (*indexes)[j]->getStrikeMapping():
                            (*strikeMappingOverride);
                        
                        indexDataMap[indexName].reset(
                            new IndexData(timeLine.size(),
                                          (*indexes)[j].getSP(),
                                          strikeMapping));
                    }
                    IndexDataSP indexData = indexDataMap[indexName];
                    double weight = (*weights)[j];
                                       
                    // Update indexData
                    double coef = weight * nameNotional * (1.0 - nameRecovery);
                    indexData->sumNotionalBespoke += weight * nameNotional;
                    for( int t = 0; t < timeLine.size(); t++) {
                        (*indexData->expectedLossesBespoke)[t] +=
                            coef * (1.0 - namesSurvivalProb[t][i]);
                    }
                }            
            }
        }

        // checks if duration calculation is really needed
        CBoolArraySP needsDurations(new CBoolArray(numNames, false));
        for(int i=0; i < numNames; ++i) {
            if (!tranche->nameDefaulted(i)){ // name has not defaulted
                // no local override, so need to compute the durations if any of the
                // indices this name maps to uses OffSpreads
                ICreditLossConfigConstSP nameLossCfg = tranche->getInnerLossConfig(i);
                const string& name = nameLossCfg->getName();
                IndexWeightsConstSP indexWeights = indexWeightsMap[name];
                IndexSkewWrapperArraySP indexes(indexWeights->getIndexNames());
                for(int j=0; j < indexes->size(); j++) {
                    (*needsDurations)[i] = (*needsDurations)[i] || ((*indexes)[j]->useOffSpreads());
                }
            }
        }

        if (!useExpectedLossRatio && (spreadRatioOverride == NULL)) {
            populateDurWeightAvgSpreadData(tranche,
                                           needsDurations,
                                           discount,
                                           timeLine, 
                                           indexWeightsMap,
                                           spreadRatiosRegions, // indices to ignore
                                           indexDataMap);
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


        const DateTime& today = tranche->getToday();
        // Some pre process for spread ratios computation
        // Note that this is NOT used when spreadRatioOverride != NULL
        // We still do it here to keep the code readable
        // Use a reduced timeline to avoid computing implied par
        // spread / duration for first point - which is set to "today"
        DateTimeArray reducedTimeLine(timeLine.begin()+1, timeLine.end());
        DoubleArray impliedSpreads(reducedTimeLine.size());
        DoubleArray durations(reducedTimeLine.size());

        // Normalises and populate index expected loss
        mapIter = indexDataMap.begin();
        for (; mapIter != indexDataMap.end(); mapIter++){
            IndexDataSP indexData = (*mapIter).second;
            
            // Normalises expectedLossesBespoke
             indexData->weight = fabs(indexData->sumNotionalBespoke) / sumSumAbsNotionalBespoke;
            if (fabs(indexData->sumNotionalBespoke) < TINY) {
                for(int t = 0; t < timeLine.size(); t++) {
                    (*indexData->expectedLossesBespoke)[t] = 0.0;
                }            
            } else {
                for(int t = 0; t < timeLine.size(); t++) {
                    (*indexData->expectedLossesBespoke)[t] /= 
                        indexData->sumNotionalBespoke;
                }
            }

            // Populates index expected loss
            // NB : this is not time consuming because clean spreads are cached
            DefaultRatesSP indexDefaultRates(indexData->skew->defaultRates());
            double indexRecovery = indexData->skew->getRecovery();
            
            for (int t = 0; t < timeLine.size(); t++) {
                (*indexData->expectedLosses)[t] =
                    (1.0 - indexRecovery) *
                    (1.0 - indexDefaultRates->calcDefaultPV(today,timeLine[t]));
            }
            
            // Initialise "name to spread ratio override" local map
            map<string, double> nameToSpreadRatioOverride;
            if (spreadRatiosRegionalOverride.get() != 0) {
                for (int i=0; i < spreadRatiosRegionalOverride->size(); ++i) {
                    nameToSpreadRatioOverride[(*spreadRatiosRegions)[i]] = 
                        (*spreadRatiosRegionalOverride)[i];
                }
            }
            // Update spread ratios
            map<string, double>::iterator nameToSpreadRatioOverrideIter = 
                nameToSpreadRatioOverride.find((*mapIter).first);
            if (spreadRatioOverride == NULL && // no global override
                nameToSpreadRatioOverrideIter == nameToSpreadRatioOverride.end()) // no local override
            {
                if (indexData->skew->useOffSpreads()) {
                    if (useExpectedLossRatio) {
                        // compute ratio of expected losses
                        for (int t = 1; t<timeLine.size(); t++) {       
                            if ((*indexData->expectedLosses)[t] != 0.0) {
                                (*indexData->spreadRatios)[t] = 
                                    (*indexData->expectedLossesBespoke)[t] / 
                                    (*indexData->expectedLosses)[t];
                            }
                            else {
                                (*indexData->spreadRatios)[t] = 1.0;
                            }
                        }
                    }
                    else {
                        // Computes index implied spreads
                        indexData->skew->impliedSpreadsAndDurations(discount,
                                                                    reducedTimeLine,
                                                                    impliedSpreads,
                                                                    durations);
                        
                        // Normalises spread ratios
                        (*indexData->spreadRatios)[0] = 1.0;
                        double normalisationFactor;
                        for (int t = 1; t<timeLine.size(); t++) {
                            normalisationFactor =
                                (*indexData->sumNotionalDurationBespoke)[t] * 
                                impliedSpreads[t-1];
                            
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
                                (*indexData->spreadRatios)[t] /= normalisationFactor;
                            } else {
                                (*indexData->spreadRatios)[t] = 1.0;
                            }
                        }
                    }
                } else {
                    // indexData->skew does not use "off spreads" -> defaults
                    // spread ratio to 1.0
                    for (int t = 0; t < timeLine.size(); t++) {
                        (*indexData->spreadRatios)[t] = 1.0;
                    }
                }
            } else {
                // Uses relevant spread ratio override
                double srOverride;
                if (spreadRatioOverride != NULL) {
                    srOverride = (*spreadRatioOverride);
                }
                else {
                    srOverride = nameToSpreadRatioOverrideIter->second;
                }
                for (int t=0; t < timeLine.size(); t++) {
                    (*indexData->spreadRatios)[t] = srOverride;
                }
            }
        }
    } catch (exception& e){
        throw ModelException(e, method);
    }
}

map<string, double> CreditMetricsBaseCorrelation::getDurationWeightedAverage(
    CreditTrancheLossConfigConstSP tranche,
    const DateTime&                maturityDate,
    const YieldCurveConstSP        discount) const 
{
    DateTimeArray timeLine(2);
    timeLine[0] = DateTime();   // dummy date - will be ignored
    timeLine[1] = maturityDate; // this is the date we really are interested in

    IndexWeightsMap indexWeightsMap;
    getIndexWeightsMap(tranche, indexWeightsMap /* (O) */);
    
    IndexDataMap indexDataMap;
    populateDurWeightAvgSpreadData(tranche,
                                   CBoolArraySP(), // include all names
                                   discount,
                                   timeLine,
                                   indexWeightsMap,
                                   StringArrayConstSP(), // no indices should be ignored
                                   indexDataMap); // (O)

    map<string, double> dwas;
    IndexDataMap::const_iterator mapIter = indexDataMap.begin();
    for (; mapIter != indexDataMap.end(); mapIter++) {
        // We are only interested in the indexData for the second tenor
        // (ie, [1]) since that is what corresponds to the maturityDate
        IndexDataSP data = (*mapIter).second;
        string name = (*mapIter).first;
        dwas[name] = (*data->spreadRatios)[1] / 
            (*data->sumNotionalDurationBespoke)[1];
    }
    return dwas;
}


// CAUTION: the spread ratios etc will NOT be computed for the first
// point of the timeline (typically today)
// - includeName: if a NULL SP is passed in, ALL (non-defaulted) names 
//   will be included
// - indicesToIgnore: if a NULL SP is passed in, NO indices will be ignored
void CreditMetricsBaseCorrelation::populateDurWeightAvgSpreadData(
    CreditTrancheLossConfigConstSP tranche,
    CBoolArraySP                   includeName,
    const YieldCurveConstSP        discount,
    const DateTimeArray&           timeLine,
    IndexWeightsMap&               indexWeightsMap,
    StringArrayConstSP             indicesToIgnore,
    IndexDataMap&                  indexDataMap) const // (O)
{
    // Some pre process for spread ratios computation
    // Note that this is NOT used when spreadRatioOverride != NULL
    // We still do it here to keep the code readable
    // Use a reduced timeline to avoid computing implied par
    // spread / duration for first point - which is set to "today"
    DateTimeArray reducedTimeLine(timeLine.begin()+1, timeLine.end());
    DoubleArray impliedSpreads(reducedTimeLine.size());
    DoubleArray durations(reducedTimeLine.size());

    int numNames = tranche->numInnerLossConfigs();
    const DateTime& today = tranche->getToday();

    for(int i=0; i < numNames; ++i) {
        if (((!!includeName) && !(*includeName)[i]) // do not include name
            || tranche->nameDefaulted(i)) // name defaulted
        {
            // This name is to be excluded from the computation
            continue;
        }
        SingleCreditAssetConstSP myAsset = tranche->nameAsset(i);
        const string& name = myAsset->getName();
        IndexWeightsConstSP indexWeights = indexWeightsMap[name];
        IndexSkewWrapperArraySP indexes(indexWeights->getIndexNames());
        CDoubleArraySP weights = indexWeights->getWeights();
            
        // support for forward starting names
        DateTimeSP protectionStartDate(new DateTime(
            tranche->nameProtectionStartDate(i)));
            
        // retrieve ICDSParSpreads object directly
        CDSPricer::impliedParSpreadsAndDurations(
            DateTimeSP(new DateTime(today)),
            protectionStartDate,
            myAsset->getParSpreadCurve(),
            discount,
            reducedTimeLine,
            impliedSpreads,
            durations);
        
        double nameNotional = tranche->nameNotional(i);
        int wholeTimelineSize = reducedTimeLine.size() + 1;
        for(int j=0; j < indexes->size(); ++j) {
            const string& indexName = (*indexes)[j].getName();
            bool ignoreIndex = false;
            int numIndices = (!indicesToIgnore ? 0 : indicesToIgnore->size());
            for (int k=0; k < numIndices; ++k) {
                if (indexName == (*indicesToIgnore)[k]) {
                    ignoreIndex = true;
                    break;
                }
            }

            if (!ignoreIndex) { // don't ignore this index
                if (indexDataMap.find(indexName) == indexDataMap.end()) {
                    indexDataMap[indexName].reset(new IndexData(timeLine.size()));
                }
                IndexDataSP indexData = indexDataMap[indexName];
                double weight = (*weights)[j];

                // Do not ignore this index. Compute spread ratio
                for(int t=1; t < wholeTimelineSize; ++t) {
                    double coefSR = weight * nameNotional * durations[t-1];
                    (*indexData->sumNotionalDurationBespoke)[t] += coefSR;
                    (*indexData->spreadRatios)[t] += coefSR * impliedSpreads[t-1];
                }
            }
        }
    }
}


/** Calculate the bespoke portfolio skew information for one given strike */
// NB : we assume lgdFloor = 0.0, lgdCap = 1.0, lgdNotional = 1.0
// This is checked in CDO when pricing with BaseCorrelation
void CreditMetricsBaseCorrelation::ccmBespokeSkewOneStrike(
    const DateTimeArray&  timeLineEff,
    double                pastPortLoss,       
    double                strike,
    SkewSurface::SkewType skewType,
    const IndexDataMap&   indexDataMap,
    double                sumOutPosNotional,
	double			      sumOutNegNotional,
    double                histAvgBeta,
    DoubleArray&          result) const // (O) result[0] not calculated
{
    // Name of this method
    static const string method = "CreditMetricsBaseCorrelation::"
        "ccmBespokeSkewOneStrike";

    try {
        // Local variables
        int t;
        double q, expLoss, expLossBespoke, indexStrike;
        IndexDataSP indexData;
        
		// Adjustment of strike and conversion to percentage
		double adjStrike = (strike - pastPortLoss ) / (sumOutPosNotional - sumOutNegNotional);

        // Caps and floors
        adjStrike = Maths::min(1.0, Maths::max(adjStrike, 0.0));
        
        // Computes the index strikes
        IndexDataMap::const_iterator mapIter = indexDataMap.begin();
        for (; mapIter != indexDataMap.end(); mapIter++){
            indexData = (*mapIter).second;
            q = indexData->strikeMapping;
            for(t=0; t<timeLineEff.size(); t++) {
                expLoss = (*indexData->expectedLosses)[t];
                expLossBespoke = (*indexData->expectedLossesBespoke)[t];
                indexStrike = adjStrike;
                if (fabs(expLossBespoke) >= TINY) {
                    indexStrike *= (q + (1.0 - q) * expLoss / expLossBespoke);
                }
                // Caps and floors
                (*indexData->indexStrikes)[t] =
                    Maths::min(1.0, Maths::max(indexStrike, 0.0));
            }
        }
        
        // Computes Beta Implied Skew
        mapIter = indexDataMap.begin();
        for (; mapIter != indexDataMap.end(); mapIter++){
            indexData = (*mapIter).second;
            
            if (betaBasis.get()) {
                // betaBasis is defined
                for(t=1; t<timeLineEff.size(); t++) {
                    result[t] += 
                        indexData->weight * 
                        (indexData->skew->getSkew(
                            (*indexData->indexStrikes)[t], timeLineEff[t], (*indexData->spreadRatios)[t], skewType)
                         + betaBasis->getSkew(
                             (*indexData->indexStrikes)[t], timeLineEff[t], skewType));
                }
            } else {
                // betaBasis is not defined : don't use it
                for(t=1; t<timeLineEff.size(); t++) {
                    result[t] += 
                        indexData->weight
                        * indexData->skew->getSkew(
                            (*indexData->indexStrikes)[t], timeLineEff[t], (*indexData->spreadRatios)[t], skewType);
                }
            }
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
        for(t=1; t<timeLineEff.size(); t++) {
            result[t] = PortfolioName::betaTweak(result[t], betaTwkSqueeze);
        }

// The following code can be used to create a debug file
// in the same format than CCM2 library
#if 0
        const string filename("H:\\BC_Strike_"+Format::toString(strike));
        FILE* fp = fopen(filename.c_str(),"w");
        if (fp != NULL) {
            
            /* DEAL DATA FIRST */
            fprintf (fp, "####### BESPOKE SKEW INFORMATION ######\n");
            fprintf (fp, "# historical average beta\n");
            fprintf (fp, "%.12f\n", histAvgBeta);
        
            fprintf (fp, "# %% strike used\n");
            fprintf (fp, "%.12f\n", adjStrike);
        
        
            fprintf (fp, "# number of reference index\n");
            fprintf (fp, "%ld\n", indexDataMap.size());
        
            fprintf (fp, "# number of points in the timeline\n");
            fprintf (fp, "%ld\n", timeLineEff.size());
        
            fprintf (fp, "# reference idx SPN and assigned weight\n");
            mapIter = indexDataMap.begin();
            for (; mapIter != indexDataMap.end(); mapIter++){
                indexData = (*mapIter).second;
                fprintf (fp, "%s\t%.12f\n", indexData->skew->getName().c_str(), indexData->weight);
            }
        
            /* expected loss of each index */
            fprintf (fp, "# EL%% of index\n");
            mapIter = indexDataMap.begin();
            for (; mapIter != indexDataMap.end(); mapIter++){
                indexData = (*mapIter).second;
                fprintf (fp, "\t%s", indexData->skew->getName().c_str());
            }
            fprintf (fp, "\n");       
            for (t=0; t<timeLineEff.size(); ++t)
            {
                fprintf (fp, "%s", timeLineEff[t].toString().c_str());        
                mapIter = indexDataMap.begin();
                for (; mapIter != indexDataMap.end(); mapIter++){
                    indexData = (*mapIter).second;
                    fprintf (fp, "\t%.12f", (*indexData->expectedLosses)[t]);
                }
                fprintf (fp, "\n");        
            }
        
            /* expected loss of the part of the bespoke related to the index */
            fprintf (fp, "# EL%% of the bespoke part related to each index\n\t");
            mapIter = indexDataMap.begin();
            for (; mapIter != indexDataMap.end(); mapIter++){
                indexData = (*mapIter).second;
                fprintf (fp, "\t%s", indexData->skew->getName().c_str());
            }
            fprintf (fp, "\n");       
            for (t=0; t<timeLineEff.size(); ++t)
            {
                fprintf (fp, "%s", timeLineEff[t].toString().c_str());        
                mapIter = indexDataMap.begin();
                for (; mapIter != indexDataMap.end(); mapIter++){
                    indexData = (*mapIter).second;
                    fprintf (fp, "\t%.12f", (*indexData->expectedLossesBespoke)[t]);
                }
                fprintf (fp, "\n");        
            }
        
            /* strike looked up for each index */
            fprintf (fp, "# strike looked up for each index\n\t");
            mapIter = indexDataMap.begin();
            for (; mapIter != indexDataMap.end(); mapIter++){
                indexData = (*mapIter).second;
                fprintf (fp, "\t%s", indexData->skew->getName().c_str());
            }
            fprintf (fp, "\n");       
            for (t=0; t<timeLineEff.size(); ++t)
            {
                fprintf (fp, "%s", timeLineEff[t].toString().c_str());        
                mapIter = indexDataMap.begin();
                for (; mapIter != indexDataMap.end(); mapIter++){
                    indexData = (*mapIter).second;
                    fprintf (fp, "\t%.12f", (*indexData->indexStrikes)[t]);
                }
                fprintf (fp, "\n");        
            }
        
            /* beta tweak looked up for each index */
            fprintf (fp, "# beta of strike for each index\n\t");
            mapIter = indexDataMap.begin();
            for (; mapIter != indexDataMap.end(); mapIter++){
                indexData = (*mapIter).second;
                fprintf (fp, "\t%s", indexData->skew->getName().c_str());
            }
            double value;
            fprintf (fp, "\n");       
            for (t=0; t<timeLineEff.size(); ++t)
            {
                fprintf (fp, "%s", timeLineEff[t].toString().c_str());        
                mapIter = indexDataMap.begin();
                for (; mapIter != indexDataMap.end(); mapIter++){
                    indexData = (*mapIter).second;
                    value = indexData->skew->getSkew(
                        (*indexData->indexStrikes)[t], timeLineEff[t], (*indexData->spreadRatios)[t], skewType);
                    fprintf (fp, "\t%.12f", value);
                }
                fprintf (fp, "\n");        
            }
        
            /* beta tweak looked up for the basis */
            fprintf (fp, "# beta of strike for the basis in relation to each index\n\t");
            mapIter = indexDataMap.begin();
            for (; mapIter != indexDataMap.end(); mapIter++){
                indexData = (*mapIter).second;
                fprintf (fp, "\t%s", indexData->skew->getName().c_str());
            }
            fprintf (fp, "\n");       
            for (t=0; t<timeLineEff.size(); ++t)
            {
                fprintf (fp, "%s", timeLineEff[t].toString().c_str());        
                mapIter = indexDataMap.begin();
                for (; mapIter != indexDataMap.end(); mapIter++){
                    indexData = (*mapIter).second;
                    value = 0.0;
                    if (betaBasis.get()) {
                        value = betaBasis->getSkew(
                            (*indexData->indexStrikes)[t], timeLineEff[t], skewType);
                    }
                    fprintf (fp, "\t%.12f", value);
                }
                fprintf (fp, "\n");        
            }
        
            fprintf (fp, "# beta tweak for squeeze\n");
            fprintf (fp, "%.12f\n", betaTwkSqueeze);
        
        
            fprintf (fp, "# Date, bespoke EL%% and beta\n");
            for (t=0; t<timeLineEff.size(); ++t)
            {
                fprintf (fp, "%s\t%.12f\t\n",
                         timeLineEff[t].toString().c_str(),
                         result[t]);
            }
            fclose(fp);
        }
#endif
    } catch (exception& e){
        throw ModelException(e, method);
    }
}


/** Calculate the average historical beta */
double CreditMetricsBaseCorrelation::calculateAverageBeta2(
    CreditTrancheLossConfigConstSP tranche) const
{
    double sumOutPosNotional, sumOutNegNotional, sumOutNetNotional, histAvgBeta;
    computeSumNotionalAndBeta2(tranche, sumOutPosNotional, sumOutNegNotional, 
                               sumOutNetNotional, histAvgBeta);
    return histAvgBeta;
}

/** Returns all the points on the skew surface (of given name) to which this
    product is sensitive. This implementation returns the points used by the
    Base Correlation methodology */
BetaSkewGridPointArrayConstSP 
CreditMetricsBaseCorrelation::getSensitiveBetaSkewPoints(
    OutputNameConstSP              outputName,
    const DateTime&                lastObservationDate, 
    CreditTrancheLossConfigConstSP tranche, 
    YieldCurveConstSP              discount,
    bool                           tweakAll, 
    const int                      numberOfNeighbours) const
{
    static const string method("CreditMetricsBaseCorrelation::"
                               "getSensitiveBetaSkewPoints");
    try {
        int i;
        
        // Special case for beta basis
        if (betaBasis.get() != 0 && outputName->equals(BETA_BASIS_NAME)) {
            return betaBasis->getAllPoints();
        }
        // get index information - need to precompute some items
        IndexWeightsMap indexWeightsMap;
        IndexDataMap indexDataMap;
        double pastPortLoss = tranche->portfolioLoss();
        
        DateTimeArraySP timeline = generateTimeline(tranche->getToday(),
                                                    lastObservationDate);
        DoubleArrayArray namesSurvivalProb;
        tranche->computeNameSurvivalProb(*timeline, namesSurvivalProb);
        double sumOutPosNotional, sumOutNegNotional, sumOutNetNotional, histAvgBeta;
        initMaps2(tranche, discount, *timeline,
                  namesSurvivalProb, pastPortLoss,
                  indexWeightsMap, indexDataMap,
                  sumOutPosNotional, sumOutNegNotional, sumOutNetNotional, histAvgBeta);

        // retrieve skew surface (including "wide spreads" surfaces)
        // corresponding to outputName
        StringArraySP names;
        string surfaceName = outputName->toString();
        bool found = false;
        IndexDataSP indexData;
        for (IndexDataMap::iterator iter = indexDataMap.begin(); 
             !found && iter != indexDataMap.end(); 
             iter++)
        {
			names = (*iter).second->skew->getSkewSurfaceNames();
            for(i=0; i < names->size();i++) {
                if ((*names)[i] == surfaceName) {
                    indexData = (*iter).second;
                    found = true;
                }
            }
		}
            
        if (found) {
            if (tweakAll) {
                // ignore first point of the timeline as it
                // is never used for base correlations computations
                DoubleArraySP reducedSpreadRatios(new DoubleArray(
                    indexData->spreadRatios->begin()+1, indexData->spreadRatios->end()));
                
                return indexData->skew->getAllPoints(
                    surfaceName, reducedSpreadRatios);
            } 
            else {
                double lowStrike, highStrike; // get hold of the strikes
                tranche->getTrancheStrikes(lowStrike, highStrike);
                BetaSkewGridPointSet pointsToTweak;
                // loop over both strikes
                for (i = 0 ; i < 2; i++) {
                    double strike = i == 0? lowStrike: highStrike;
                    BetaSkewGridPointSet pointsOneStrike(
                        getSensitiveBetaSkewPointsOneStrike(*timeline,
                                                            pastPortLoss,
                                                            strike,
                                                            indexData,
                                                            sumOutPosNotional,
                                                            surfaceName, 
                                                            numberOfNeighbours));
                    
                    pointsToTweak.insert(
                        pointsOneStrike.begin(), pointsOneStrike.end()); 
                }
                return BetaSkewGridPoint::toArray(pointsToTweak);
            }
        } 
        else {
            return BetaSkewGridPointArrayConstSP();
        }

    } 
    catch (exception& e){
        throw ModelException(e, method);
    }
}


/** */
BetaSkewGridPointSet 
CreditMetricsBaseCorrelation::getSensitiveBetaSkewPointsOneStrike(
    const DateTimeArray& timeLine,
    double               pastPortLoss,       
    double               strike,
    IndexDataSP          indexData,
    double               sumOutPosNotional,
    string               surfaceName, 
    const int            numberOfNeighbours) const
{
    static const string method("CreditMetricsBaseCorrelation::"
                               "getSensitiveBetaSkewPointsOneStrike");
    try {
        // Adjustment of strike and conversion to percentage
        double adjStrike =
            Maths::max(strike - pastPortLoss, 0.0) / sumOutPosNotional;
        
        // Caps and floors
        adjStrike = Maths::min(1.0, Maths::max(adjStrike, 0.0));
        
        // Computes the index strikes
        BetaSkewGridPointSet pointsSet;
        double q = indexData->strikeMapping;
        for(int t = 1; t < timeLine.size(); t++) {
            double expLoss = (*indexData->expectedLosses)[t];
            double expLossBespoke = (*indexData->expectedLossesBespoke)[t];
            double indexStrike = adjStrike;
            if (fabs(expLossBespoke) >= TINY) {
                indexStrike *= (q + (1.0 - q) * expLoss / expLossBespoke);
            }
            // Caps and floors
            double currentStrike = 
                Maths::min(1.0, Maths::max(indexStrike, 0.0));
            BetaSkewGridPointSet localPoints(
                indexData->skew->getNeighbours(currentStrike, 
                                               timeLine[t], 
                                               surfaceName, 
                                               (*indexData->spreadRatios)[t], 
                                               numberOfNeighbours)); 
            pointsSet.insert(localPoints.begin(), localPoints.end());
        }
        return pointsSet;
    } 
    catch (exception& e){
        throw ModelException(e, method);
    }
}        

/** Compression ratio tweak support */
string CreditMetricsBaseCorrelation::sensName(CompressionRatioTwk* shift) const {
    throw ModelException(
        "CreditMetricsBaseCorrelation::sensName(CompressionRatioTwk* shift)",
        "Not implemented");            
}

/** Compression ratio tweak support */  
bool CreditMetricsBaseCorrelation::sensShift(CompressionRatioTwk* shift) {
    double value = shift->getShiftSize();
    compressionRatio += value;
    
    //Cap and floor tweaked value
    if (compressionRatio < 0.0) {
        compressionRatio = 0.0;
    }
    if (compressionRatio > 1.0) {
        compressionRatio = 1.0;
    }
    return false;
}

/** Strike mapping override tweak support */
string CreditMetricsBaseCorrelation::sensName(StrikeMappingOverrideTwk* shift) const{
    throw ModelException(
        "CreditMetricsBaseCorrelation::sensName(StrikeMappingOverrideTwk* shift)",
        "Not implemented");            
}

/** Strike mapping override tweak support */
bool CreditMetricsBaseCorrelation::sensShift(StrikeMappingOverrideTwk* shift) {
    if (strikeMappingOverride != NULL) {
        double value = shift->getShiftSize();
        *strikeMappingOverride += value;
        
        //Cap and floor tweaked value
        if (*strikeMappingOverride < 0.0) {
            *strikeMappingOverride = 0.0;
        }
        if (*strikeMappingOverride > 1.0) {
            *strikeMappingOverride = 1.0;
        }
    }
    return false;
}


/** Base Correlation Strike mapping override tweak support */
string CreditMetricsBaseCorrelation::sensName(BCStrikeMappingOverrideTwk* shift) const{
    throw ModelException(
        "CreditMetricsBaseCorrelation::sensName(BCStrikeMappingOverrideTwk* shift)",
        "Not implemented");            
}

/** Base Correlation Strike mapping override tweak support */
bool CreditMetricsBaseCorrelation::sensShift(BCStrikeMappingOverrideTwk* shift) {
    static const string method = 
        "CreditMetricsBaseCorrelation::sensShift(BCStrikeMappingOverrideTwk* shift)";

    double newStrikeMappingOverride = shift->applyShift(*strikeMappingOverride);

    // Checks if the new strikeMappingOverride is valid
    if ((newStrikeMappingOverride < 0.0) || (newStrikeMappingOverride > 1.0)) {
        throw ModelException(method,
                             "Shifted strikeMappingOverride (="
                             + Format::toString(newStrikeMappingOverride) + 
                             ") outside of [0,1].");   
    }
    *strikeMappingOverride = newStrikeMappingOverride;
    return false; /* none of our components has a strikeMapping
                   * type sensitivity */
}


/** Compression ratio and betaBasis setting support */
bool CreditMetricsBaseCorrelation::sensShift(QuasiContractualBaseCorrelation* shift) {
    static const string method = 
        "CreditMetricsBaseCorrelation::sensShift(QuasiContractualBaseCorrelation* shift)";

    // Set the compression ratio to the shift value
    double newCompressionRatio = shift->getCompressionRatio();

    // Checks if the new compressionRatio is valid
    try {
        checkRange("compressionRatio", newCompressionRatio, 0.0, 1.0);
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }
    compressionRatio = newCompressionRatio;
    return true;
}


/** Invoked by the containing model after the instrument data is fetched
    ie after CInstrument::GetMarket is invoked. This invokes parent method
    and retrieves indexWeightsOverride data */
void CreditMetricsBaseCorrelation::postFetchMarketData(
    IModel*            model,
    MarketDataConstSP market)
{
    CreditMetricsModel::postFetchMarketData(model, market);
    if (indexWeightsOverride.get() != 0) {
        for(int i=0; i<indexWeightsOverride->size(); i++) {
            (*indexWeightsOverride)[i]->getMarket(model, market.get());
        }
    }
}
/** Called immediately after object constructed */
void CreditMetricsBaseCorrelation::validatePop2Object() {
    CreditMetricsModel::validatePop2Object();
    static const string method = "CreditMetricsBaseCorrelation::validatePop2Object";
    try{
        int i;
        
        // Set name for beta basis (if defined)
        if (betaBasis.get() != 0) {
            betaBasis->setName(BETA_BASIS_NAME);
        }
      
        // Checks range for strikeMappingOverride     
        if (strikeMappingOverride != 0) {
            checkRange("strikeMappingOverride", *strikeMappingOverride, 0.0, 1.0);
        }
        
        // spreadRatioOverride and spreadRatiosRegionalOverride are mutually exclusive
        if (spreadRatioOverride != 0 && spreadRatiosRegionalOverride.get() != 0)
        {
            throw ModelException(method,
                "spreadRatioOverride and spreadRatiosRegionalOverride are mutually exclusive");
        }

        // Checks size of spreadRatiosRegionalOverride and spreadRatiosRegions
        if (spreadRatiosRegionalOverride.get() != 0)
        {
            if (spreadRatiosRegions.get() == 0)
            {
                throw ModelException(method,
                    "No spreadRatiosRegions specified for spreadRatiosRegionalOverride");
            }
            if (spreadRatiosRegionalOverride->size() != spreadRatiosRegions->size())
            {
                throw ModelException(method,
                    "spreadRatiosRegions and spreadRatiosRegionalOverride don't have the same number of elements");
            }
        }
        
        if (lowerBCBetasOverride.get() != 0) {
            // Checks consistency between lowerBCBetasOverride and bcBetasTimeline
            if (bcBetasTimeline.get() == 0) {
                throw ModelException(method,"No timeline associated to lowerBCBetasOverride");
            }
            if (lowerBCBetasOverride->size() != bcBetasTimeline->size()){
                throw ModelException(method,
                                     "Inconsistency: bcBetasTimeline and lowerBCBetasOverride have not the same size");
            }
            
            // Checks range
            for (i = 0; i < lowerBCBetasOverride->size(); i++) {
                checkRange(
                    "lowerBCBetasOverride[" + Format::toString(i) + "]",
                    (*lowerBCBetasOverride)[i], -1.0, 1.0);
            }
        }

        if (upperBCBetasOverride.get() != 0) {
            // Checks consistency between upperBCBetasOverride and bcBetasTimeline
            if (bcBetasTimeline.get() == 0) {
                throw ModelException(method,"No timeline associated to upperBCBetasOverride");
            }
            if (upperBCBetasOverride->size() != bcBetasTimeline->size()){
                throw ModelException(method,
                                     "Inconsistency: bcBetasTimeline and upperBCBetasOverride have not the same size");
            }
            
            // Checks range
            for (i = 0; i < upperBCBetasOverride->size(); i++) {
                checkRange(
                    "upperBCBetasOverride[" + Format::toString(i) + "]",
                    (*upperBCBetasOverride)[i], -1.0, 1.0);
            }
        }
        
        // Checks range for compressionRatio
        checkRange("compressionRatio", compressionRatio, 0.0, 1.0);

        // Checks namesOverride and indexWeightsOverride have same length
        if (namesOverride->size() != indexWeightsOverride->size()) {
            throw ModelException(method,
                                 "namesOverride and indexWeightsOverride have different sizes");            
        }
    } catch (exception& e){
        throw ModelException(e, "CreditMetricsBaseCorrelation::validatePop2Object");
    }
}

/** Destructor */
CreditMetricsBaseCorrelation::~CreditMetricsBaseCorrelation() {
    delete strikeMappingOverride;
    delete spreadRatioOverride;
}

/** Public constructor */
// Initialises compressionRatio and strikeMappingOverride pointer
CreditMetricsBaseCorrelation::CreditMetricsBaseCorrelation() : 
    CreditMetricsModel(TYPE),
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
    useExpectedLossRatio(false)
{
}


/** Private constructor */
// Initialises compressionRatio and strikeMappingOverride pointer
CreditMetricsBaseCorrelation::CreditMetricsBaseCorrelation(const CClassConstSP& clazz) : 
    CreditMetricsModel(clazz),
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
    useExpectedLossRatio(false) {}

/** Default constructor */
IObject* CreditMetricsBaseCorrelation::defaultConstructor(){
    return new CreditMetricsBaseCorrelation(CreditMetricsBaseCorrelation::TYPE);
}

/** Invoked when Class is 'loaded' */
void CreditMetricsBaseCorrelation::load(CClassSP& clazz){
    clazz->setPublic(); // make visible to EAS/spreadsheet
    REGISTER(CreditMetricsBaseCorrelation, clazz);
    SUPERCLASS(CreditMetricsModel);
    IMPLEMENTS(TweakableWith<CompressionRatioTwk>);
    IMPLEMENTS(TweakableWith<StrikeMappingOverrideTwk>);
    IMPLEMENTS(TweakableWith<BCStrikeMappingOverrideTwk>);
    IMPLEMENTS(QuasiContractualBaseCorrelation::IShift);
    EMPTY_SHELL_METHOD(defaultConstructor);
    
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
    
    FIELD(betaBasis,
          "Optional parameter used to offset the skew surface");
    FIELD_MAKE_OPTIONAL(betaBasis);
    
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
}

/** Default value for compression ratio (theta) */
double const CreditMetricsBaseCorrelation::DEFAULT_THETA = 1.0;

/** Default strike used in the calculation of CCC */
double const CreditMetricsBaseCorrelation::EQUITY_STRIKE_FOR_CCC = 0.03;

/** Name of the (optional) beta basis */
string const CreditMetricsBaseCorrelation::BETA_BASIS_NAME = "Beta Basis";

CClassConstSP const CreditMetricsBaseCorrelation::TYPE = CClass::registerClassLoadMethod(
    "CreditMetricsBaseCorrelation", typeid(CreditMetricsBaseCorrelation), CreditMetricsBaseCorrelation::load);

bool CreditMetricsBaseCorrelationLoad() {
    return (CreditMetricsBaseCorrelation::TYPE != 0);
}
DRLIB_END_NAMESPACE
