//----------------------------------------------------------------------------
//
//   Group       : QR
//
//   Author      : Sebastien Gay
//
//   Date        : March 2006
//
//----------------------------------------------------------------------------

#ifndef QR_CCMRFLLOSSCALCULATOR_HPP
#define QR_CCMRFLLOSSCALCULATOR_HPP

#include "edginc/CCMLossCalculator.hpp"
#include "edginc/CCMConvolution.hpp"
#include "edginc/FORWARD_DECLARE.hpp"

DRLIB_BEGIN_NAMESPACE

FORWARD_DECLARE(CCMTrancheCalculatorLegacy);

/** Implementation of ITrancheLossCalculator using 'Composite Copula Model' */
class CONVOLUTION_DLL CCMRFLLossCalculator: public CCMLossCalculatorBase {
public:
    virtual ~CCMRFLLossCalculator();

    /** Constructor - takes in full timeline to allow for optimisations. If
        computeCondCurve if false then lossConditional in loss() below will be
        set to zero. */
    CCMRFLLossCalculator(
        const DateTimeArray&           timeline,        /* (I) */
        CreditTrancheLossConfigConstSP tranche,         /* (I) */
        double                         lossUnit,        /* (I) */
        const bool                     recoverNotional, /* (I) */    
        CounterPartyCreditConstSP      cpty);           /* (I) */

    /** Same as original lossKey method but allows the betas used to be
        overridden. The length of the array must be the same as the number
        of names or be empty (=> no beta overrides). This is used by Base
        Correlation */
    virtual IKey* lossKey(
        int                timePoint,   /* (I) do the calculation for 
                                           this timepoint */
        const DoubleArray& betaOverride) const; // (I) using these betas

private:
    CCMTrancheCalculatorLegacySP calculatorTemplate;
};


DRLIB_END_NAMESPACE
#endif
