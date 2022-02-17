//----------------------------------------------------------------------------
//
//   Group       : QR
//
//   Author      : Sebastien Gay
//
//   Date        : March 2006
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/CreditAsset.hpp"
#include "edginc/CreditMetricsRFLFastLossCalculator.hpp"
#include "edginc/CCMTrancheFastLossCalculatorLegacy.hpp"
#include "edginc/CreditTrancheLossConfig.hpp"
#include "edginc/CounterPartyCredit.hpp"

DRLIB_BEGIN_NAMESPACE

CreditMetricsRFLFastLossCalculator::~CreditMetricsRFLFastLossCalculator(){}

/** Constructor - takes in full timeline to allow for optimisations. If
    computeCondCurve if false then lossConditional in loss() below will be
    set to zero. */
CreditMetricsRFLFastLossCalculator::CreditMetricsRFLFastLossCalculator(
        const DateTimeArray&           timeline,        /* (I) */
        CreditTrancheLossConfigConstSP tranche,         /* (I) */
        const bool                     recoverNotional, /* (I) */
        CounterPartyCreditConstSP      cpty):           /* (I) */
    CreditMetricsLossCalculatorBase(timeline, tranche, cpty),
    calculatorTemplate(new CCMTrancheFastLossCalculatorLegacy(true,
                                                              true,
                                                              recoverNotional,
                                                              tranche,
                                                              cpty))
{}


/** Same as original lossKey method but allows the betas used to be
    overridden. The length of the array must be the same as the number
    of names or be empty (=> no beta overrides). This is used by Base
    Correlation */
ITrancheLossCalculator::IKey* CreditMetricsRFLFastLossCalculator::lossKey(
    int                timePoint,   /* (I) do the calculation for 
                                       this timepoint */
    const DoubleArray& betaOverride) const // (I) using these betas
{
    // To do: scrap CCMTrancheCalculator. Do everything here instead.
    try{
        // create new calculate from our template - this is inefficient since
        // the Credit Metrics convolution algorithm leaves it unchanged
        CCMTrancheFastLossCalculatorLegacySP trancheCalc(
            new CCMTrancheFastLossCalculatorLegacy(calculatorTemplate.get()));
        trancheCalc->setupCM(survivalProb[timePoint],
                             computeCondCurve? counterPartyProb[timePoint]: 0.0,
                             betaOverride);
        return new Key(trancheCalc);
    } catch (exception& e){
        throw ModelException(e, "CreditMetricsRFLFastLossCalculator::lossKey");
    }
}

DRLIB_END_NAMESPACE
