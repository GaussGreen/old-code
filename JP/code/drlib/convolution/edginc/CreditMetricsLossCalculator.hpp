//----------------------------------------------------------------------------
//
//   Group       : QR
//
//   Author      : Mark A Robson
//
//   Date        : 23 Dec 2005
//
//
//----------------------------------------------------------------------------

#ifndef QR_CREDITMETRICSLOSSCALCULATOR_HPP
#define QR_CREDITMETRICSLOSSCALCULATOR_HPP

#include "edginc/CreditMetricsLossCalculatorBase.hpp"
#include "edginc/FORWARD_DECLARE.hpp"

DRLIB_BEGIN_NAMESPACE
FORWARD_DECLARE(CCMTrancheCalculatorLegacy);

/** Implemenetation of ITrancheLossCalculator using 'CreditMetrics' */
class CONVOLUTION_DLL CreditMetricsLossCalculator: public CreditMetricsLossCalculatorBase{
public:
    virtual ~CreditMetricsLossCalculator();

    /** Constructor - takes in full timeline to allow for optimisations. If
        computeCondCurve if false then lossConditional in loss() below will be
        set to zero. */
    CreditMetricsLossCalculator(
        const DateTimeArray&           timeline,        /* (I) */
        CreditTrancheLossConfigConstSP tranche,         /* (I) */
        double                         lossUnit,        /* (I) */
        const bool                     recoverNotional, /* (I) */
        CounterPartyCreditConstSP      cpty);           /* (I) */

    /** Same as above but allows the betas used to be overridden. The
        length of the array must be the same as the number of names or
        be empty (=> no beta overrides). This is used by Base Correlation */
    virtual IKey* lossKey(
        int                timePoint,   /* (I) do the calculation for 
                                           this timepoint */
        const DoubleArray& betaOverride) const; // (I) using these betas

private:
    CCMTrancheCalculatorLegacySP calculatorTemplate; // holds invariant data 
};

DRLIB_END_NAMESPACE
#endif
