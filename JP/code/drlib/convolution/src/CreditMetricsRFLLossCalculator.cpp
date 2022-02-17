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
#include "edginc/CreditMetricsRFLLossCalculator.hpp"
#include "edginc/CCMTrancheCalculatorLegacy.hpp"
#include "edginc/CreditTrancheLossConfig.hpp"
#include "edginc/CounterPartyCredit.hpp"

DRLIB_BEGIN_NAMESPACE

CreditMetricsRFLLossCalculator::~CreditMetricsRFLLossCalculator(){}

/** Constructor - takes in full timeline to allow for optimisations. If
    computeCondCurve if false then lossConditional in loss() below will be
    set to zero. */
CreditMetricsRFLLossCalculator::CreditMetricsRFLLossCalculator(
        const DateTimeArray&           timeline,        /* (I) */
        CreditTrancheLossConfigConstSP tranche,         /* (I) */
        double                         lossUnit,        /* (I) */
        const bool                     recoverNotional, /* (I) */
        CounterPartyCreditConstSP      cpty) :           /* (I) */
    CreditMetricsLossCalculatorBase(timeline, tranche, cpty),
    calculatorTemplate(new CCMTrancheCalculatorLegacy(
        true, true, recoverNotional, tranche, lossUnit, cpty))
{}


/** Same as original lossKey method but allows the betas used to be
    overridden. The length of the array must be the same as the number
    of names or be empty (=> no beta overrides). This is used by Base
    Correlation */
ITrancheLossCalculator::IKey* CreditMetricsRFLLossCalculator::lossKey(
    int                timePoint,   /* (I) do the calculation for 
                                       this timepoint */
    const DoubleArray& betaOverride) const // (I) using these betas
{
    try{
        // create new calculate from our template - this is inefficient since
        // the Credit Metrics convolution algorithm leaves it unchanged
        refCountPtr<CCMTrancheCalculatorLegacy> trancheCalc(
            new CCMTrancheCalculatorLegacy(calculatorTemplate.get()));
        trancheCalc->setupCM(survivalProb[timePoint],
                             computeCondCurve? counterPartyProb[timePoint]: 0.0,
                             betaOverride);
        
        /* the following line is the key part of the algorithm */
        // to do: probably remove all of these parameters
        trancheCalc->convolution(0 /* control - to be removed */,
                                 0 /* model to be removed */, timePoint,
                                 1 /* presumably means something? */);
        return new Key(trancheCalc);
    } catch (exception& e){
        throw ModelException(e, "CreditMetricsLossCalculator::lossKey");
    }
}

DRLIB_END_NAMESPACE
