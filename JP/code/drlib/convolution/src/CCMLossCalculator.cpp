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
#include "edginc/CreditAsset.hpp"
#include "edginc/CCMLossCalculator.hpp"
#include "edginc/CCMTrancheCalculatorLegacy.hpp"
#include "edginc/CreditTrancheLossConfig.hpp"
#include "edginc/CounterPartyCredit.hpp"


DRLIB_BEGIN_NAMESPACE

CCMLossCalculator::~CCMLossCalculator(){}

/** Constructor - takes in full timeline to allow for optimisations. If
    computeCondCurve if false then lossConditional in loss() below will be
    set to zero. */
CCMLossCalculator::CCMLossCalculator(
        const DateTimeArray&           timeline,        /* (I) */
        CreditTrancheLossConfigConstSP tranche,         /* (I) */
        double                         lossUnit,        /* (I) */
        const bool                     recoverNotional, /* (I) */    
        CounterPartyCreditConstSP      cpty):           /* (I) */
    CCMLossCalculatorBase(timeline, tranche, cpty),
    calculatorTemplate(new CCMTrancheCalculatorLegacy(
        false, false, recoverNotional, tranche, lossUnit, cpty))
{}

/** Same as original lossKey method but allows the betas used to be
    overridden. The length of the array must be the same as the number
    of names or be empty (=> no beta overrides). This is used by Base
    Correlation */
ITrancheLossCalculator::IKey* CCMLossCalculator::lossKey(
    int                timePoint,   /* (I) do the calculation for 
                                       this timepoint */
    const DoubleArray& betaOverride) const // (I) using these betas
{
    // To do: scrap CCMTrancheCalculator. Do everything here instead.
    try{
        // create new calculate from our template - this could be more efficient
        refCountPtr<CCMTrancheCalculatorLegacy> trancheCalc(
            new CCMTrancheCalculatorLegacy(calculatorTemplate.get()));
        trancheCalc->setupCCM(survivalProb[timePoint],
                              floorSurvivalProb[timePoint],
                              computeCondCurve? 
                              counterPartyProb[timePoint]: 0.0,
                              betaOverride);
        /* the following line is the key part of the algorithm */
        // to do: probably remove all of these parameters
        trancheCalc->convolution(0 /* control - to be removed */,
                                 0 /* model to be removed */, timePoint,
                                 1 /* presumably means something? */);
        return new Key(trancheCalc);
    } catch (exception& e){
        throw ModelException(e, "CCMTrancheCalculator::lossKey");
    }
}

DRLIB_END_NAMESPACE
