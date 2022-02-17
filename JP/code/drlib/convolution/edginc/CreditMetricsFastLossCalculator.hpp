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

#ifndef QR_CREDITMETRICSFASTLOSSCALCULATOR_HPP
#define QR_CREDITMETRICSFASTLOSSCALCULATOR_HPP
#include "edginc/CreditMetricsLossCalculatorBase.hpp"
#include "edginc/CCMFastLossCalculator.hpp"
#include "edginc/FORWARD_DECLARE.hpp"

DRLIB_BEGIN_NAMESPACE

/* Forward declare the classes refered to from this file if we are not
 * interested in their storage properties (ie, they are used through (smart)
 * pointers and therefore their include files are not required here - they
 * will be required in the .cpp file though)*/
FORWARD_DECLARE(CreditTrancheLossConfig);
FORWARD_DECLARE(CounterPartyCredit);
class ConvolutionProduct;

/** Implemenetation of ITrancheLossCalculator using 'fast' 'CreditMetrics' 
    which is accurate only as the number of names tends to infinity */
class CONVOLUTION_DLL CreditMetricsFastLossCalculator: public CreditMetricsLossCalculatorBase {
public:
    virtual ~CreditMetricsFastLossCalculator();

    /** Constructor - takes in full timeline to allow for optimisations. If
        computeCondCurve if false then lossConditional in loss() below will be
        set to zero. */
    CreditMetricsFastLossCalculator(
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
