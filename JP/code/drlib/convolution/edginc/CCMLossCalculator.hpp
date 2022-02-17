//----------------------------------------------------------------------------
//
//   Group       : QR
//
//   Author      : Mark A Robson
//
//   Date        : 23 Dec 2005
//
//----------------------------------------------------------------------------

#ifndef QR_CCMLOSSCALCULATOR_HPP
#define QR_CCMLOSSCALCULATOR_HPP

#include "edginc/CCMLossCalculatorBase.hpp"
#include "edginc/CCMConvolution.hpp"
#include "edginc/FORWARD_DECLARE.hpp"

DRLIB_BEGIN_NAMESPACE

class ConvolutionProduct;
FORWARD_DECLARE(CCMTrancheCalculatorLegacy);

/** Implementation of ITrancheLossCalculator using 'Composite Copula Model' */
class CONVOLUTION_DLL CCMLossCalculator: public CCMLossCalculatorBase {
public:
    virtual ~CCMLossCalculator();

    /** Constructor - takes in full timeline to allow for optimisations. If
        computeCondCurve if false then lossConditional in loss() below will be
        set to zero. */
    CCMLossCalculator(
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

protected:
    CCMTrancheCalculatorLegacySP calculatorTemplate;
};


DRLIB_END_NAMESPACE
#endif
