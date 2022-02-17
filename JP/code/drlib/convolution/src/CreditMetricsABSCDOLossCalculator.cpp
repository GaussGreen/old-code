//----------------------------------------------------------------------------
//
//   Group       : QR
//
//   Author      : Jay Wang
//
//   Date        : Fed 2006
//
//----------------------------------------------------------------------------
#include "edginc/config.hpp"
#include "edginc/AbsCdoParameters.hpp"
#include "edginc/CCMLossUnit.hpp"
#include "edginc/CCMTrancheUtils.hpp"
#include "edginc/ConvolutionProduct.hpp"
#include "edginc/CreditMetricsABSCDOLossCalculator.hpp"
#include "edginc/Control.hpp"
#include "edginc/Results.hpp"
#include "edginc/SingleCreditAsset.hpp"
#include "edginc/CounterPartyCredit.hpp"
#include "edginc/ABSCDOConvolution.hpp"
#include "edginc/OutputRequest.hpp"
#include "edginc/CreditTrancheLossConfig.hpp"
#include "edginc/CounterPartyCredit.hpp"

DRLIB_BEGIN_NAMESPACE

/** Constructor - takes in full timeline to allow for optimisations. If
    computeCondCurve if false then lossConditional in loss() below will be
    set to zero. */
CreditMetricsABSCDOLossCalculator::CreditMetricsABSCDOLossCalculator(
        const DateTimeArray&           timeline,               /* (I) */
        CreditTrancheLossConfigConstSP tranche,                /* (I) */
        double                         lossUnit,               /* (I) */
        CounterPartyCreditConstSP      counterParty,           /* (I) */
        bool                           useSaddlePoint,         /* (I) */
        int                            numSaddlePoints,        /* (I) */
        int                            lossMarketFactors,      /* (I) */
        int                            decretionMarketFactors, /* (I) */
        bool                           useFastConvolution,     /* (I) */
        const DoubleArray&             lossBetaOverride,       /* (I) */
        const DoubleArray&             decBetaOverride,        /* (I) */
        const DoubleArray&             lossDecBetaOverride,    /* (I) */
        bool                           calcExpectedNotional) : /* (I) */
    CreditMetricsLossCalculatorBase(timeline, tranche, counterParty), 
    timeline(timeline), 
    lossUnit(lossUnit),
    useSaddlePoint(useSaddlePoint), 
    numSaddlePoints(numSaddlePoints),
    lossMarketFactors(lossMarketFactors), 
    decretionMarketFactors(decretionMarketFactors),
    useFastConvolution(useFastConvolution), 
    calcExpectedNotional(calcExpectedNotional)
{
    try {
        // set up name parameters
        pastLoss = tranche->portfolioLoss();

        populateBasketInfoCM(tranche, lossBetaOverride, decBetaOverride, basketInfo);
        cpty = createCptyInfoCM(counterParty);

        /*  renormalize notionals and allocate density */
        int n = basketInfo.size();
        ntl.resize(n);
        for (int i = 0; i < n; ++i) {
            basketInfo[i]->ntl /= lossUnit;
            basketInfo[i]->ntlLoss /= lossUnit;
            basketInfo[i]->ntlRecovery /= lossUnit;
            ntl[i] = basketInfo[i]->ntl;
        }
        
        CCMLossUnit::calcLossDensityBounds(ntl, maxl, maxs);
        
        density.resize(maxs+maxl+1);
        densityCond.resize(maxs+maxl+1);
        
        initializeOutputResults(timeline, tranche);

        createMarketScenarios(lossMarketFactors, decretionMarketFactors);

		AbsCdoParametersConstSP modelParam = 
			AbsCdoParametersConstSP::dynamicCast(tranche->getEngineParams(AbsCdoParameters::TYPE));
        double lossDecretionBeta = lossDecBetaOverride.empty() ? 
            modelParam->getLossDecBeta() : lossDecBetaOverride[0];

        lossMarketWeight = lossDecretionBeta;
        decretionMarketWeight = sqrt(1 - lossMarketWeight * lossMarketWeight);
        createConditionalCurves(timeline);
    } 
    catch (exception& e){
        throw ModelException(e, "CreditMetricsABSCDOLossCalculator::ctor");
    }   
}


void CreditMetricsABSCDOLossCalculator::createMarketScenarios(
    const int lossMarketFactors, const int decretionMarketFactors)
{
    try {
        // create matrix of MarketScenarios
        decIntMethod.calcPoints(decretionMarketFactors);
        lossIntMethod.calcPoints(lossMarketFactors);
        
        scenarios.reserve(decretionMarketFactors);
        for (int i = 0; i < decretionMarketFactors; ++i) {
            ABSCDOConvolution::MarketScenarioArray tmp;
            tmp.reserve(lossMarketFactors);
            for (int j = 0; j < lossMarketFactors; ++j) {
                ABSCDOConvolution::MarketScenarioSP scenario (
                    new ABSCDOConvolution::MarketScenario(lossIntMethod.points[j],
                                                          lossIntMethod.weights[j],
                                                          decIntMethod.points[i],
                                                          decIntMethod.weights[i]));
                tmp.push_back(scenario);
            }
            scenarios.push_back(tmp);
        }
    } catch (exception& e) {
        throw ModelException(e, "CreditMetricsABSCDOLossCalculator::createMarketScenarios");
    }
}
void CreditMetricsABSCDOLossCalculator::createConditionalCurves(const DateTimeArray& timeline)
{
    try {
        // create conditional default and dec curve for each name under each MarketScenario
        // memory-intensive operation
        for (int i = 0; i < decretionMarketFactors; ++i) {
            ABSCDOConvolution::MarketScenarioArray& currentDecScenario = scenarios[i];
            for (int j = 0; j < lossMarketFactors; ++j) {
                ABSCDOConvolution::MarketScenarioSP& currentLossScenario = currentDecScenario[j];
                currentLossScenario->createConditionalCurves(basketInfo, timeline, 
                                                             lossMarketWeight, decretionMarketWeight);
            }
        }
    } catch (exception& e) {
        throw ModelException(e, "CreditMetricsABSCDOLossCalculator::createConditionalCurves");
    }
}

void CreditMetricsABSCDOLossCalculator::initializeOutputResults(
    const DateTimeArray&           timeline,
    CreditTrancheLossConfigConstSP tranche)
{
    initialNotional = tranche->trancheOutstandingNotional(true); // always recover notional
    portNotional = tranche->portfolioNotional();
    int numDates = timeline.size();

    expectedBottomLoss.resize(numDates);
    expectedBottomLoss[0].date = timeline[0];
    expectedBottomLoss[0].amount = 0;

    expectedTopLoss.resize(numDates);
    expectedTopLoss[0].date = timeline[0];
    expectedTopLoss[0].amount = 0;

    expectedNotional.resize(numDates);
    expectedNotional[0].date = timeline[0];
    expectedNotional[0].amount = initialNotional;    
}

void CreditMetricsABSCDOLossCalculator::populateBasketInfoCM(
    CreditTrancheLossConfigConstSP     tranche,      /* (I) */
    const DoubleArray&                 lossBetaOverride,
    const DoubleArray&                 decBetaOverride,
    ABSCDOConvolution::NameParamArray& basketInfo)   // (O)  
{
    int numNames = tranche->numInnerLossConfigs();
    basketInfo.clear();
    basketInfo.reserve(numNames);
    bool lossOverride = !lossBetaOverride.empty();
    bool decOverride = !decBetaOverride.empty();
    // the whole approach here is pretty weak since all but one of the 
    // parameters is time independent yet we copy everything over each time.
    // Should review all of this - not sure how necessary it really is
    for (int i = 0; i < numNames; ++i) {
        ABSCDOConvolution::NameParamSP nm(new ABSCDOConvolution::NameParam());
        // NB Non credit metrics params are set to zero in constructor
        nm->beta          = lossOverride ? lossBetaOverride[0] : tranche->nameBeta(i);
        nm->decBeta       = decOverride ? decBetaOverride[0] : tranche->nameDecBeta(i);
        SingleCreditAssetConstSP asset(tranche->nameAsset(i));
        nm->nameId        = asset->getName();
        nm->ntl           = tranche->nameDefaulted(i) ? 0.0 : tranche->nameNotional(i);
        nm->R             = tranche->nameRecovery(i);
        nm->decretion     = asset->getParSpreadCurve()->getPrepayCurve();
        nm->defaultrate   = asset->getParSpreadCurve()->defaultRates();
        nm->index         = i;
        nm->ntlLoss       = nm->ntl * (1. - nm->R);
        nm->ntlRecovery   = nm->ntl * nm->R;

        //add to the basket
        basketInfo.push_back(nm);
    }
}


void CreditMetricsABSCDOLossCalculator::setupCM(
    const DoubleArray&        survivalProb, /* name survival proba */
    double                    counterPartyProb) const /* name survival proba */
{
    for (size_t i = 0; i < basketInfo.size(); i++){
        basketInfo[i]->survival = survivalProb[i];
    }
    if (cpty.get()){
        cpty->survival = counterPartyProb;
    }
}


ABSCDOConvolution::NameParamSP CreditMetricsABSCDOLossCalculator::createCptyInfoCM(
    CounterPartyCreditConstSP counterParty) /* (I) */
{
    ABSCDOConvolution::NameParamSP cpty;
    if (!!counterParty) {
        cpty.reset(new ABSCDOConvolution::NameParam());
        cpty->beta          = counterParty->getBeta();
        cpty->nameId        = counterParty->getName();
        cpty->R             = counterParty->getRecovery();
    }
    return cpty;
}


/** Calculate the expected loss for specified timepoint and strikes. 
    no caching supported yet */
void CreditMetricsABSCDOLossCalculator::loss(
        int     timePoint,        // (I) do the calculation for this timepoint
        double  k1,               /* (I) lower strike      */
        double  k2,               /* (I) upper strike      */
        double& loss,             /* (O) tranche loss amt  */
        double& lossCond) const  /* (O) tranche loss amt cond on cpty
                                     surviving */
{
    try{
        // set up convolution parameters, to be removed
        setupCM(survivalProb[timePoint], false);
     
        if (useFastConvolution) {
            double lossTop;
            ABSCDOConvolution::calcABSCDOTrancheLossFast(
                                  timeline[timePoint], k1, k2, lossUnit, portNotional,
                                  scenarios, basketInfo, loss, lossTop);
            expectedTopLoss[timePoint].date = timeline[timePoint];
            expectedTopLoss[timePoint].amount = lossTop;
            expectedBottomLoss[timePoint].date = timeline[timePoint];
            expectedBottomLoss[timePoint].amount = loss;
        } else if (useSaddlePoint) {
            double lossTop;
            ABSCDOConvolution::calcABSCDOTrancheLossSaddle(
                                  timeline[timePoint], k1, k2, lossUnit, portNotional,
                                  lossMarketFactors, lossIntMethod.points, lossIntMethod.weights,
                                  scenarios, basketInfo, loss, lossTop, numSaddlePoints);
            expectedTopLoss[timePoint].date = timeline[timePoint];
            expectedTopLoss[timePoint].amount = lossTop;
            expectedBottomLoss[timePoint].date = timeline[timePoint];
            expectedBottomLoss[timePoint].amount = loss;
        } else {
            ABSCDOConvolution::calcABSCDOLossDistribution(
                                  timeline[timePoint], maxl + 1L, scenarios, basketInfo);
            ABSCDOConvolution::calcABSCDOLossDensities(density, LossDistribution::TOPLOSS,
                                                       scenarios, basketInfo.size(), maxl + 1L);
            expectedTopLoss[timePoint].date = timeline[timePoint];
            expectedTopLoss[timePoint].amount = 
                                CCMTrancheUtils::calcExpectedTrancheLoss(
                                                                portNotional-k2,
                                                                portNotional-k1,
                                                                pastLoss,
                                                                lossUnit,
                                                                density,
                                                                maxs,
                                                                maxl);
            ABSCDOConvolution::calcABSCDOLossDensities(density, LossDistribution::BOTTOMLOSS,
                                                       scenarios, basketInfo.size(), 
                                                       maxl + 1L);
            expectedBottomLoss[timePoint].date = timeline[timePoint];
            loss = expectedBottomLoss[timePoint].amount =
                                CCMTrancheUtils::calcExpectedTrancheLoss(
                                                            k1,
                                                            k2,
                                                            pastLoss,
                                                            lossUnit,
                                                            density,
                                                            maxs,
                                                            maxl);
        }
            
        lossCond = loss; // do not support county party default calc

        if (calcExpectedNotional) {
            expectedNotional[timePoint].date = timeline[timePoint];
            expectedNotional[timePoint].amount = 
                initialNotional - expectedBottomLoss[timePoint].amount
                    - expectedTopLoss[timePoint].amount;

            // patch
            if (expectedNotional[timePoint].amount <= 0) {
                expectedTopLoss[timePoint].amount = initialNotional - expectedBottomLoss[timePoint].amount - 1;
                expectedNotional[timePoint].amount = 
                    initialNotional - expectedBottomLoss[timePoint].amount
                        - expectedTopLoss[timePoint].amount;
            }
        }
    } catch (exception& e){
        throw ModelException(e, "CreditMetricsABSCDOLossCalculator::loss");
    }
}

void CreditMetricsABSCDOLossCalculator::getTopLoss(
        int     timePoint,        // (I) do the calculation for this timepoint
        double& loss,             /* (O) tranche loss amt  */
        double& lossCond) const  /* (O) tranche loss amt cond on cpty
                                     surviving */
{
    loss = lossCond = expectedTopLoss[timePoint].amount;
}

/** store calculator results */
void CreditMetricsABSCDOLossCalculator::storeResults(
    Results* results, Control* control) const
{
    if (!control->isPricing())
        return;
    int numDates = timeline.size();

    OutputRequest* request = control->requestsOutput(
            OutputRequest::TRANCHE_ABS_BOTTOM_LOSS_CURVE);
    if (request){
        CashFlowArraySP eloss(new CashFlowArray(numDates));
        for (int i = 0; i < numDates; ++i) {
            (*eloss)[i].date = expectedBottomLoss[i].date;
            (*eloss)[i].amount = expectedBottomLoss[i].amount;
        }
        results->storeRequestResult(request, eloss);
    }

    request = control->requestsOutput(
            OutputRequest::TRANCHE_ABS_TOP_LOSS_CURVE);
    if (request){
        CashFlowArraySP edec(new CashFlowArray(numDates));
        for (int i = 0; i < numDates; ++i) {
            (*edec)[i].date = expectedTopLoss[i].date;
            (*edec)[i].amount = expectedTopLoss[i].amount;
        }
        results->storeRequestResult(request, edec);
    }

    request = control->requestsOutput(
            OutputRequest::TRANCHE_ABS_NOTIONAL_CURVE);
    if (request){
        CashFlowArraySP enot(new CashFlowArray(numDates));
        for (int i = 0; i < numDates; ++i) {
            (*enot)[i].date = expectedNotional[i].date;
            (*enot)[i].amount = expectedNotional[i].amount;
        }
        results->storeRequestResult(request, enot);
    }
}

/** Returns a key used to optimise repeated calculations of losses
    at the same timepoint. Not supported yet */
ITrancheLossCalculator::IKey* CreditMetricsABSCDOLossCalculator::lossKey(
        int timePoint) const // (I) do the calculation for this timepoint
{
    throw ModelException("Not implemented!");
}

/** Same as above but allows the betas used to be overridden. The
    length of the array must be the same as the number of names or
    be empty (=> no beta overrides). This is used by Base Correlation 
    not supported yet */
ITrancheLossCalculator::IKey* CreditMetricsABSCDOLossCalculator::lossKey(
    int                timePoint,   /* (I) do the calculation for 
                                       this timepoint */
    const DoubleArray& betaOverride) const // (I) using these betas
{
    throw ModelException("Not implemented!");
}

/////////////////////////////////////////////////////////////////////////////////

/** Calculate the expected loss for specified timepoint and strikes. 
    no caching supported yet */
void ProxyCalculator::loss(
    int     timePoint,        // (I) do the calculation for this timepoint
    double  k1,               /* (I) lower strike      */
    double  k2,               /* (I) upper strike      */
    double& loss,             /* (O) tranche loss amt  */
    double& lossCond) const  /* (O) tranche loss amt cond on cpty
                                 surviving */
{
    parent->getTopLoss(timePoint, loss, lossCond);
}

/** Returns a key used to optimise repeated calculations of losses
    at the same timepoint. Not supported yet */
ITrancheLossCalculator::IKey* ProxyCalculator::lossKey(
        int timePoint) const // (I) do the calculation for this timepoint
{
    throw ModelException("Not implemented!");
}

/** Same as above but allows the betas used to be overridden. The
    length of the array must be the same as the number of names or
    be empty (=> no beta overrides). This is used by Base Correlation 
    not supported yet */
ITrancheLossCalculator::IKey* ProxyCalculator::lossKey(
    int                timePoint,   /* (I) do the calculation for 
                                       this timepoint */
    const DoubleArray& betaOverride) const // (I) using these betas
{
    throw ModelException("Not implemented!");
}

/////////////////////////////////////////////////////////////////////////////////

ABSCDOBaseCorrelationLossCalculator::ABSCDOBaseCorrelationLossCalculator(
        const DateTimeArray&           timeline,            /* (I) */
        CreditTrancheLossConfigConstSP tranche,             /* (I) */
        double                         lossUnit,            /* (I) */
        CounterPartyCreditConstSP      counterParty,        /* (I) */
        const DateTime&                maturity,           /* (I) */
        bool                           useSaddlePoint,      /* (I) */
        int                            numSaddlePoints,     /* (I) */
        int                            lossMarketFactors,   /* (I) */
        int                            decretionMarketFactors, /* (I) */
        bool                           fastConvolution,      /* (I) */
		bool							useLossSkew,		/* (I) */
		bool							useDecSkew,			/* (I) */
		bool							useLossDecSkew,		/* (I) */
        bool                           authoriseNegativeEL) :
    timeline(timeline), authoriseNegativeEL(authoriseNegativeEL)
{
    try {
        // get strikes
        double portNotional = tranche->portfolioNotional();
        tranche->getTrancheStrikes(lowStrike, highStrike);

        // get base correlations
		AbsCdoParametersConstSP modelParam = 
			AbsCdoParametersConstSP::dynamicCast(tranche->getEngineParams(AbsCdoParameters::TYPE));

        double lowLossBC, highLossBC;
        double lowDecBC, highDecBC;
        double lowLossDecBC, highLossDecBC;
        DoubleArray lowLossBCOverride, lowDecBCOverride, lowLossDecBCOverride;
        DoubleArray highLossBCOverride, highDecBCOverride, highLossDecBCOverride;

        SkewSurface::SkewType st = fastConvolution? SkewSurface::FAST : SkewSurface::NORMAL;

        lowBetas.reset(new DoubleArray(3));
        highBetas.reset(new DoubleArray(3));
        if (useLossSkew) {
			SkewSurfaceConstSP lossBCSkew = modelParam->getLossSkew();
            lowLossBC = lossBCSkew->getSkew(lowStrike/portNotional, maturity, st);
            if (lowLossBC > 1 || lowLossBC < -1)
                throw ModelException("interpolated loss BC at low strike out of bound");
            highLossBC = lossBCSkew->getSkew(highStrike/portNotional, maturity, st);
            if (highLossBC > 1 || highLossBC < -1)
                throw ModelException("interpolated loss BC at high strike out of bound");
            lowLossBCOverride.resize(1);
            lowLossBCOverride[0] = lowLossBC;
            highLossBCOverride.resize(1);
            highLossBCOverride[0] = highLossBC;
            (*lowBetas)[0] = lowLossBC;
            (*highBetas)[0] = highLossBC;
        }

        if (useDecSkew) {
			SkewSurfaceConstSP decBCSkew = modelParam->getDecSkew();
            lowDecBC = decBCSkew->getSkew(lowStrike/portNotional, maturity, st);
            if (lowDecBC > 1 || lowDecBC < -1)
                throw ModelException("interpolated decretion BC at low strike out of bound");
            highDecBC = decBCSkew->getSkew(highStrike/portNotional, maturity, st);
            if (highDecBC > 1 || highDecBC < -1)
                throw ModelException("interpolated decretion BC at high strike out of bound");
            lowDecBCOverride.resize(1);
            lowDecBCOverride[0] = lowDecBC;
            highDecBCOverride.resize(1);
            highDecBCOverride[0] = highDecBC;
            (*lowBetas)[1] = lowDecBC;
            (*highBetas)[1] = highDecBC;
        }

        if (useLossDecSkew) {
			SkewSurfaceConstSP lossDecBCSkew = modelParam->getLossDecSkew();
            lowLossDecBC = lossDecBCSkew->getSkew(lowStrike/portNotional, maturity, st);
            if (lowLossDecBC > 1 || lowLossDecBC < -1)
                throw ModelException("interpolated loss-decretion BC at low strike out of bound");
            highLossDecBC = lossDecBCSkew->getSkew(highStrike/portNotional, maturity, st);
            if (highLossDecBC > 1 || highLossDecBC < -1)
                throw ModelException("interpolated loss-decretion BC at high strike out of bound");
            lowLossDecBCOverride.resize(1);
            lowLossDecBCOverride[0] = lowLossDecBC;
            highLossDecBCOverride.resize(1);
            highLossDecBCOverride[0] = highLossDecBC;
            (*lowBetas)[2] = lowLossDecBC;
            (*highBetas)[2] = highLossDecBC;
        }

        
        // create calculators
        if (lowStrike < 1e-6) {
            calculatorLow = CreditMetricsABSCDOLossCalculatorSP(0);
        } 
        else {
            calculatorLow.reset(new CreditMetricsABSCDOLossCalculator(
                                                        timeline,
                                                        tranche,
                                                        lossUnit,
                                                        counterParty,
                                                        useSaddlePoint,
                                                        numSaddlePoints,
                                                        lossMarketFactors,
                                                        decretionMarketFactors,
                                                        fastConvolution,
                                                        lowLossBCOverride,
                                                        lowDecBCOverride,
                                                        lowLossDecBCOverride,
                                                        false));
        }
        calculatorHigh.reset(new CreditMetricsABSCDOLossCalculator(
                                                        timeline,
                                                        tranche,
                                                        lossUnit,
                                                        counterParty,
                                                        useSaddlePoint,
                                                        numSaddlePoints,
                                                        lossMarketFactors,
                                                        decretionMarketFactors,
                                                        fastConvolution,
                                                        highLossBCOverride,
                                                        highDecBCOverride,
                                                        highLossDecBCOverride,
                                                        false));

        // initialize output results
        initialNotional = tranche->trancheOutstandingNotional(true); // always recover notional
        int numDates = timeline.size();

        expectedBottomLoss.resize(numDates);
        expectedBottomLoss[0].date = timeline[0];
        expectedBottomLoss[0].amount = 0; // XXX should be past loss

        expectedTopLoss.resize(numDates);
        expectedTopLoss[0].date = timeline[0];
        expectedTopLoss[0].amount = 0; // XXX

        expectedNotional.resize(numDates);
        expectedNotional[0].date = timeline[0];
        expectedNotional[0].amount = initialNotional;    
    } 
    catch (exception& e) {
        throw ModelException(e, "ABSCDOBaseCorrelationLossCalculator::"
                             "ABSCDOBaseCorrelationLossCalculator");
    }
}


/** Calculate the expected loss for specified timepoint and strikes. 
    no caching supported yet */
void ABSCDOBaseCorrelationLossCalculator::loss(
    int     timePoint,        // (I) do the calculation for this timepoint
    double  k1,               /* (I) lower strike      */
    double  k2,               /* (I) upper strike      */
    double& loss,             /* (O) tranche loss amt  */
    double& lossCond) const  /* (O) tranche loss amt cond on cpty
                                 surviving */
{
    try {
        double lowLoss, lowLossCond;
        double highLoss, highLossCond;
        if (calculatorLow.get() == 0) {
            lowLoss = lowLossCond = 0;
        } else 
            calculatorLow->loss(timePoint, 0, lowStrike, lowLoss, lowLossCond);
        calculatorHigh->loss(timePoint, 0, highStrike, highLoss, highLossCond);
        loss = highLoss - lowLoss;
        if (!authoriseNegativeEL && loss < 0)
            loss = 0.0;
        lossCond = loss;

        if (loss < 0)
            throw ModelException("negative tranche expected loss, "
                    "most likely this means the base correlation curve is too steep");

        double topLoss;
        if (calculatorLow.get() == 0) {
            lowLoss = lowLossCond = 0;
        } else
            calculatorLow->getTopLoss(timePoint, lowLoss, lowLossCond);
        calculatorHigh->getTopLoss(timePoint, highLoss, highLossCond);
        topLoss = highLoss - lowLoss;
        if (!authoriseNegativeEL && topLoss < 0)
            topLoss = 0.0;

        if (topLoss < 0)
            throw ModelException("negative tranche top expected loss, "
                    "most likely this means the base correlation curve is too steep");

        expectedTopLoss[timePoint].date = timeline[timePoint];
        expectedTopLoss[timePoint].amount = topLoss;
        expectedBottomLoss[timePoint].date = timeline[timePoint];
        expectedBottomLoss[timePoint].amount = loss;

        expectedNotional[timePoint].date = timeline[timePoint];
        expectedNotional[timePoint].amount = 
            initialNotional - expectedBottomLoss[timePoint].amount
                - expectedTopLoss[timePoint].amount;

        // patch
        if (expectedNotional[timePoint].amount <= 0) {
            expectedTopLoss[timePoint].amount = initialNotional - expectedBottomLoss[timePoint].amount - 1;
            expectedNotional[timePoint].amount = 
                initialNotional - expectedBottomLoss[timePoint].amount
                    - expectedTopLoss[timePoint].amount;
        }        
    } catch (exception& e) {
        throw ModelException(e, "ABSCDOBaseCorrelationLossCalculator::loss");
    }
}

/** calculate tranche top loss */
void ABSCDOBaseCorrelationLossCalculator::getTopLoss(
    int     timePoint,
    double& loss,
    double& lossCond) const
{
    loss = lossCond = expectedTopLoss[timePoint].amount;
}


void ABSCDOBaseCorrelationLossCalculator::storeResults(Results* results, 
                                                       Control* control) const
{
    if (!control->isPricing())
        return;
    int numDates = timeline.size();

    OutputRequest* request = control->requestsOutput(
            OutputRequest::TRANCHE_ABS_BOTTOM_LOSS_CURVE);
    if (request){
        CashFlowArraySP eloss(new CashFlowArray(numDates));
        for (int i = 0; i < numDates; ++i) {
            (*eloss)[i].date = expectedBottomLoss[i].date;
            (*eloss)[i].amount = expectedBottomLoss[i].amount;
        }
        results->storeRequestResult(request, eloss);
    }

    request = control->requestsOutput(
            OutputRequest::TRANCHE_ABS_TOP_LOSS_CURVE);
    if (request){
        CashFlowArraySP edec(new CashFlowArray(numDates));
        for (int i = 0; i < numDates; ++i) {
            (*edec)[i].date = expectedTopLoss[i].date;
            (*edec)[i].amount = expectedTopLoss[i].amount;
        }
        results->storeRequestResult(request, edec);
    }

    request = control->requestsOutput(
            OutputRequest::TRANCHE_ABS_NOTIONAL_CURVE);
    if (request){
        CashFlowArraySP enot(new CashFlowArray(numDates));
        for (int i = 0; i < numDates; ++i) {
            (*enot)[i].date = expectedNotional[i].date;
            (*enot)[i].amount = expectedNotional[i].amount;
        }
        results->storeRequestResult(request, enot);
    }

    request = control->requestsOutput(
            OutputRequest::TRANCHE_LOWER_BC_BETAS);
    if (request) {
        results->storeRequestResult(request, lowBetas);
    }

    request = control->requestsOutput(
            OutputRequest::TRANCHE_UPPER_BC_BETAS);
    if (request) {
        results->storeRequestResult(request, highBetas);
    }
}

/** Returns a key used to optimise repeated calculations of losses
    at the same timepoint. Not supported yet */
ITrancheLossCalculator::IKey* ABSCDOBaseCorrelationLossCalculator::lossKey(
        int timePoint) const // (I) do the calculation for this timepoint
{
    throw ModelException("Not implemented!");
}

/** Same as above but allows the betas used to be overridden. The
    length of the array must be the same as the number of names or
    be empty (=> no beta overrides). This is used by Base Correlation 
    not supported yet */
ITrancheLossCalculator::IKey* ABSCDOBaseCorrelationLossCalculator::lossKey(
    int                timePoint,   /* (I) do the calculation for 
                                       this timepoint */
    const DoubleArray& betaOverride) const // (I) using these betas
{
    throw ModelException("Not implemented!");
}

DRLIB_END_NAMESPACE
