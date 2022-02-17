//----------------------------------------------------------------------------
//
//   Group       : QR
//
//   Author      : Sebastien Gay
//
//   Date        : March 2006
//
//----------------------------------------------------------------------------

#ifndef QR_CREDITMETRICSRFLFASTLOSSCALCULATOR_HPP
#define QR_CREDITMETRICSRFLFASTLOSSCALCULATOR_HPP

#include "edginc/CreditMetricsFastLossCalculator.hpp"
#include "edginc/CreditMetricsRFLLossCalculator.hpp"
#include "edginc/CCMRFLFastLossCalculator.hpp"

DRLIB_BEGIN_NAMESPACE

/** Implemenetation of ITrancheLossCalculator using 'fast' 'CreditMetrics' 
    which is accurate only as the number of names tends to infinity */
class CONVOLUTION_DLL CreditMetricsRFLFastLossCalculator: public CreditMetricsLossCalculatorBase {
public:
    virtual ~CreditMetricsRFLFastLossCalculator();

    /** Constructor - takes in full timeline to allow for optimisations. If
        computeCondCurve if false then lossConditional in loss() below will be
        set to zero. */
    CreditMetricsRFLFastLossCalculator(
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
