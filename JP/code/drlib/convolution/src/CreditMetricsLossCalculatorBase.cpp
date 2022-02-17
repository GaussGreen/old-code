//----------------------------------------------------------------------------
//
//   Group       : QR
//
//   Author      : Mark A Robson
//
//   Date        : 23 Dec 2005
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/CreditMetricsLossCalculatorBase.hpp"
#include "edginc/FixedTrancheLossCalculator.hpp"
#include "edginc/ConvolutionProduct.hpp"
#include "edginc/CCMLossCalculator.hpp"
#include "edginc/TrancheLossCalculatorLegacy.hpp"
#include "edginc/SingleCreditAsset.hpp"
#include "edginc/CounterPartyCredit.hpp"
#include "edginc/RflOnlyParameters.hpp"
#include "edginc/CreditTrancheLossConfig.hpp"
#include "edginc/CounterPartyCredit.hpp"

DRLIB_BEGIN_NAMESPACE
// avoid having another source file
ITrancheLossCalculator::ITrancheLossCalculator(){}
ITrancheLossCalculator::~ITrancheLossCalculator(){}
ITrancheLossCalculator::IKey::~IKey(){}
ITrancheLossCalculator::IKey::IKey(){}
IFixedTrancheLossCalculator::~IFixedTrancheLossCalculator(){}


/** Retrieves RFL stlye engine params from Credit asset. Fails if they are
    not there or are of wrong type */
RflOnlyParametersConstSP CreditMetricsLossCalculatorBase::getRFLEngineParameters(
    SingleCreditAssetConstSP asset)
{
    CreditEngineParametersConstSP params(
        asset->getEngineParams(RflOnlyParameters::TYPE));
    return RflOnlyParametersConstSP::dynamicCast(params); 
}


/** Tiny little utility class to morph a ITrancheLossCalculator into a
    IFixedTrancheLossCalculator */
class ITrancheLossCalculator::FixedStrikes: 
    public virtual IFixedTrancheLossCalculator
{
public:
    FixedStrikes(double                        lowerStrike,  
                 double                        upperStrike,  
                 ITrancheLossCalculatorConstSP lossCalculator):
        lowerStrike(lowerStrike), upperStrike(upperStrike), 
        lossCalculator(lossCalculator){}

    /** Calculate the expected loss for specified time point. */
    virtual void loss(
        int     timePoint,        // (I) do the calculation for this timepoint
        double& loss,             /* (O) tranche loss amt  */
        double& lossCond) const{  /* (O) tranche loss amt cond on cpty
                                         surviving */
        lossCalculator->loss(timePoint, lowerStrike, upperStrike, 
                             loss, lossCond);

    }

    /** store calculator results, empty by default */
    void storeResults(Results* result, Control* control) const {
        lossCalculator->storeResults(result, control); }

private:
    double                        lowerStrike;
    double                        upperStrike;  
    ITrancheLossCalculatorConstSP lossCalculator;
};

/** Utility method to build a IFixedTrancheLossCalculator using a
    ITrancheLossCalculator. Implementation in
    CreditMetricsLossCalculatorBase.cpp */
IFixedTrancheLossCalculator* 
ITrancheLossCalculator::createFixedTrancheLossCalculator(
    double                        lowerStrike,  
    double                        upperStrike,  
    ITrancheLossCalculatorConstSP lossCalculator)
{
    return new FixedStrikes(lowerStrike, upperStrike, lossCalculator);
}

CreditMetricsLossCalculatorBase::~CreditMetricsLossCalculatorBase()
{}

/** Constructor - takes in full timeline to allow for optimisations. If
    computeCondCurve if false then lossConditional in loss() below will be
    set to zero. */
CreditMetricsLossCalculatorBase::CreditMetricsLossCalculatorBase(
        const DateTimeArray&           timeline, /* (I) */
        CreditTrancheLossConfigConstSP tranche,  /* (I) */
        CounterPartyCreditConstSP      cpty):    /* (I) */
    computeCondCurve(!!cpty)
{
    try {
        tranche->computeNameSurvivalProb(timeline, survivalProb);
        // 'true' is there temporarily to match existing tests - to do
        if (true /*computeCondCurve*/) {
            computeCounterPartySurvivalProb(cpty,
                                            timeline, 
                                            counterPartyProb);
        }
    } 
    catch (exception& e) {
        throw ModelException(e, "CreditMetricsLossCalculatorBase::"
                             "CreditMetricsLossCalculatorBase(tranche)");
    }
}


void CreditMetricsLossCalculatorBase::computeCounterPartySurvivalProb(
    CounterPartyCreditConstSP counterParty, // (I)
    const DateTimeArray& timeline,          // (I)
    DoubleArray& counterPartyProb)          // (O)
{
    static const string method("CreditMetricsLossCalculatorBase::"
                               "computeCounterPartySurvivalProb");
    try {
        if (!!counterParty) {
            counterPartyProb.resize(timeline.size()); // reserve space
            // get hold of relevant market data
            SingleCreditAssetConstSP asset(counterParty->getAsset());
            DefaultRatesSP defaultRates(asset->
                                        getParSpreadCurve()->defaultRates());
            const DateTime& today = defaultRates->getValueDate();
            for (int t=0; t < timeline.size(); ++t) {
                counterPartyProb[t] = defaultRates->calcDefaultPV(today,
                                                                  timeline[t]);
            }
        }
    } 
    catch (exception& e) {
        throw ModelException(e, "ConvolutionProduct::"
                             "computeCounterPartySurvivalProb");
    }

}


/** store calculator results, empty by default */
void CreditMetricsLossCalculatorBase::storeResults(
    Results* result,
    Control* control) const
{}

/** Calculate the expected loss for specified time point and strikes. 
    Note if repeated calculations at the same time point are required then
    IKey::calc should be used instead. This implementation just does
    "lossKey()->loss(...)" */
void CreditMetricsLossCalculatorBase::loss(
    int     timePoint,      /* (I) do the calculation for this timepoint */
    double  k1,             /* (I) lower strike      */
    double  k2,             /* (I) upper strike      */
    double& loss,           /* (O) tranche loss amt  */
    double& lossCond) const /* (O) tranche loss amt cond on cpty surviving */
{
    IKeySP key(lossKey(timePoint));
    key->loss(k1, k2, loss, lossCond);
}

/** Returns a key used to optimise repeated calculations of losses
    at the same time point. */
ITrancheLossCalculator::IKey* CreditMetricsLossCalculatorBase::lossKey(
    int timePoint) const // (I) do the calculation for this timepoint
{
    return lossKey(timePoint, DoubleArray());
}

CreditMetricsLossCalculatorBase::Key::Key(
    ITrancheLossCalculatorLegacySP trancheCalc): trancheCalc(trancheCalc){}
CreditMetricsLossCalculatorBase::Key::~Key(){}

/** Calculates the appropriate rate/factor between the two dates */
void CreditMetricsLossCalculatorBase::Key::loss(
    double  k1,               /* (I) lower strike      */
    double  k2,               /* (I) upper strike      */
    double& loss,             /* (O) tranche loss amt  */
    double& lossCond) const   /* (O) tranche loss amt cond on cpty surviving */
{
    trancheCalc->loss(k1, k2, loss, lossCond);
}

/** Allows a CreditMetricsLossCalculatorBase to be used as
    IFixedTrancheLossCalculator with overrides for the betas */
class CreditMetricsLossCalculatorBase::BetaOverride:
    public virtual IFixedTrancheLossCalculator
{
    double                                 lowerStrike;
    double                                 upperStrike;  
    DoubleArrayArraySP                     betas; //override for [time][name]
    CreditMetricsLossCalculatorBaseConstSP lossCalculator;

public:
    BetaOverride(double                      lowerStrike,  
                 double                      upperStrike,
                 DoubleArrayArraySP          betas, //override for [time][name]
                 CreditMetricsLossCalculatorBaseConstSP lossCalculator):
        lowerStrike(lowerStrike), upperStrike(upperStrike), betas(betas),
        lossCalculator(lossCalculator)
    {}

    /** Calculate the expected loss for specified time point. */
    virtual void loss(
        int     timePoint,       // (I) do the calculation for this timepoint
        double& loss,            // (O) tranche loss amt
        double& lossCond) const  // (O) tranche loss amt cond on cpty surviving
    {
        // create IKey with beta override
        IKeySP key(lossCalculator->lossKey(timePoint, (*betas)[timePoint]));
        key->loss(lowerStrike, upperStrike, loss, lossCond);
    }

    /** store calculator results, empty by default */
    void storeResults(Results* result, Control* control) const {
        lossCalculator->storeResults(result, control); 
    }
};
    
/** Utility method to build a IFixedTrancheLossCalculator using a
    CreditMetricsLossCalculatorBase where you can override the betas
    at each time point/for each name. */
IFixedTrancheLossCalculator* 
CreditMetricsLossCalculatorBase::createFixedTrancheLossCalculator(
    double                                 lowerStrike,  
    double                                 upperStrike,
    DoubleArrayArraySP                     betas, //override for [time][name]
    CreditMetricsLossCalculatorBaseConstSP lossCalculator)
{
    return new BetaOverride(lowerStrike, upperStrike, betas, lossCalculator);
}

DRLIB_END_NAMESPACE
