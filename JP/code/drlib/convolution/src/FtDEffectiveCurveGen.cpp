//----------------------------------------------------------------------------
//
// Group       : Credit Hybrids QR
//
// Description : A generator of FtDEffectiveCurveSPs for the fee and contingent
//               legs, using the information in the instrument and model
//
// Date        : October 2006
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/mathlib.hpp"
#include "edginc/AtomicArray.hpp"
#include "edginc/NToDefaultLossConfig.hpp"
#include "edginc/ICondLossDistributionsGen.hpp"
#include "edginc/ConditionalSurvProbIntegrand.hpp"
#include "edginc/ConditionalFtDLossIntegrand.hpp"
#include "edginc/EffectiveCurve.hpp"
#include "edginc/FtDEffectiveCurve.hpp"
#include "edginc/FtDEffectiveCurveGen.hpp"


DRLIB_BEGIN_NAMESPACE


FtDEffectiveCurveGen::~FtDEffectiveCurveGen()
{}

/** Constructor */
FtDEffectiveCurveGen::FtDEffectiveCurveGen(
        NToDefaultLossConfigConstSP ntdLossCfg,
        IModelConfigMapperConstSP   mapper,
        IConditionalDefaultsModelSP condDefaultsModel) :
    modelConfigMapper(mapper),
    condDefaultsModel(condDefaultsModel)
{
    static const string method("FtDEffectiveCurveGen::FtDEffectiveCurveGen");

    // Verify (even if it has been done already) that this NtD is really an FtD
    if (ntdLossCfg->getDefaultNumber() == 1) {
        ftdLossCfg = ntdLossCfg;
    }
    else {
        throw ModelException(method, "The loss config object is not an FtD "
                             "(i.e., it is not an NToDefaultLossConfig "
                             "with N=1)");
    }
}


void FtDEffectiveCurveGen::getFeeAndCtgEffectiveCurve(
    IDiscountCurveRiskySP&   feeEffCurve, // (O)
    IDiscountCurveRiskySP&   ctgEffCurve, // (O) 
    YieldCurveConstSP        discount,
    const IConvolutionModel* convolutionModel,
    const DateTime&          lastObservationDate)
{
    static const string method("FtDEffectiveCurveGen::getFeeAndCtgEffectiveCurve");
    
    const DateTime& valueDate = ftdLossCfg->getToday();

    DateTimeArraySP timeline = ftdLossCfg->getTimeLine();
    if (!timeline || timeline->size() == 0) {
        // Ask the model to generate the timeline
        timeline = convolutionModel->generateTimeline(valueDate, 
                                                      lastObservationDate);
    }

    // Use a reduced timeline if 1st point = value date
    // (in which case expected loss = 0)
    const int timeLineOffset = ((*timeline)[0] == valueDate ? 1 : 0);
    DateTimeArraySP reducedTimeLine(new DateTimeArray(
        timeline->begin() + timeLineOffset, timeline->end()));


    const int nbReducedDates = reducedTimeLine->size();
    const int nbNames = ftdLossCfg->numInnerLossConfigs();
    
    // Vector of size nbReducedDates
    // For each date, it contains an IKeyArraySP corresponding to names
    // in the portfolio
    vector<ICondLossDistributionsGenKeyArraySP> keysByDate(nbReducedDates);
    for (int t=0; t < nbReducedDates; ++t) {
        keysByDate[t].reset(new ICondLossDistributionsGenKeyArray(nbNames));
    }
    
    // Loss given default for each of the names - used to compute ctgEffCurve
    DoubleArraySP namesLoss(new DoubleArray(nbNames));

    // Temporary objects, to save constructing and destructing them in the
    // for loop below
    ICreditLossGenConstSP lossGen;
    ICondLossDistributionsGenConstSP condLossGen;
    ICreditLossConfigConstSP name; // we know it really is an ISingleDefaultCreditLossConfig
    ICreditLossModelConfigConstSP innerModel;
    ICondLossDistributionsGenKeyArrayConstSP keys;
    for (int i=0; i < nbNames; ++i) {
        name = ftdLossCfg->getInnerLossConfig(i);
            
        // Retrieves relevant ICreditLossModelConfig
        innerModel = modelConfigMapper->innerModel(name);
        
        // Retrieves loss generator
        lossGen = innerModel->lossGenerator(
            name, IModelConfigMapperConstSP(modelConfigMapper));
            
        // Casts lossGen: here want an ICondLossDistributionsGen
        condLossGen = DYNAMIC_POINTER_CAST<const ICondLossDistributionsGen>(lossGen);
        
        if (!condLossGen) {
            throw ModelException(method, "Unable to cast asset " +
                                 Format::toString(i) +
                                 " ('" + name->getName() + "') "
                                 "into a ICondLossDistributionsGen.");
        }
            
        // Initialises keys (NB: potentially a time consuming operation)
        keys = condLossGen->initialise(*reducedTimeLine);
        
        // Populates keysByDate
        for (int t=0; t < nbReducedDates; ++t) {
            (*keysByDate[t])[i] = (*keys)[t];
        }

        // Populate namesLoss
        (*namesLoss)[i] = name->maxPossibleLoss();
    }

    // Now we have all the per-name information, so go and 
    // compute the effective curves

    // Fee leg first
    const int nbDates = nbReducedDates + timeLineOffset; // ie, timeline->size();
    DoubleArray survivalProbabilities(nbDates);
    survivalProbabilities[0] = 1.0; // ie, no default today
    for (int t=timeLineOffset; t < nbDates; ++t) {
        const int reducedT = t-timeLineOffset;
        ConditionalSurvProbIntegrand condSurvProbFunction(
            keysByDate[reducedT],
            condDefaultsModel->marketFactorDimension(),
            (*timeline)[t],
            ftdLossCfg.get());
        
        survivalProbabilities[t] = condDefaultsModel->integrateCondFunction(
            &condSurvProbFunction, keysByDate[reducedT], (*timeline)[t]);
    }

    string interpolationStyle = convolutionModel->getLossInterpolation();

    feeEffCurve.reset(new EffectiveCurve(valueDate, 
                                         discount,
                                         *timeline,
                                         survivalProbabilities,
                                         interpolationStyle));

    // And now ctgEffCurve:
    // It is not possible to pre-compute the data required to create a typical
    // EffectiveCurve in the case of the FtDEffectiveCurveGen - what we do is 
    // to pass in all the data required to compute it on the fly.
    // The issue is that the recoveryDelay (which is only know at the point
    // of pricing the contingent leg) is required in the integration accross
    // market factors... so we need to postpone computing everything until the
    // last minute, when everything is known
    ctgEffCurve.reset(new FtDEffectiveCurve(valueDate, 
                                            discount,
                                            timeline,
                                            keysByDate,
                                            condDefaultsModel,
                                            namesLoss,
                                            ftdLossCfg->notional()));
}

DRLIB_END_NAMESPACE
