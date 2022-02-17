//----------------------------------------------------------------------------
//
//   Group       : QR
//
//   Author      : Mark A Robson
//
//   Date        : 23 Dec 2005
//
//----------------------------------------------------------------------------

#ifndef QR_CCMFASTLOSSCALCULATOR_HPP
#define QR_CCMFASTLOSSCALCULATOR_HPP

#include "edginc/CCMLossCalculatorBase.hpp"
#include "edginc/DateTime.hpp"

DRLIB_BEGIN_NAMESPACE
class CCMTrancheFastLossCalculatorLegacy;
typedef refCountPtr<
    CCMTrancheFastLossCalculatorLegacy> CCMTrancheFastLossCalculatorLegacySP;

/** Implementation of ITrancheLossCalculator using 'fast' Composite Copula
    Model which is accurate only as the number of names tends to infinity */
class CONVOLUTION_DLL CCMFastLossCalculator: public CCMLossCalculatorBase {
public:
    virtual ~CCMFastLossCalculator();

    /** Constructor - takes in full timeline to allow for optimisations. If
        computeCondCurve if false then lossConditional in loss() below will be
        set to zero. */
    CCMFastLossCalculator(
        const DateTimeArray&           timeline,        /* (I) */
        CreditTrancheLossConfigConstSP tranche,         /* (I) */
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
    CCMTrancheFastLossCalculatorLegacySP calculatorTemplate;
};

DRLIB_END_NAMESPACE
#endif
