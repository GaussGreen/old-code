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

#ifndef QR_BASECORRELATIONLOSSCALCULATOR_HPP
#define QR_BASECORRELATIONLOSSCALCULATOR_HPP
#include "edginc/FixedTrancheLossCalculator.hpp"
#include "edginc/AtomicArray.hpp"
#include "edginc/FORWARD_DECLARE.hpp"
#include "edginc/Control_forward.hpp"
#include "edginc/Results_forward.hpp"

DRLIB_BEGIN_NAMESPACE
class ConvolutionProduct;
FORWARD_DECLARE(CreditMetricsLossCalculatorBase);

/** Implementation of ITrancheLossCalculator using 'Base Correlation' adjustment
    on top of an existing ITrancheLossCalculator */
class CONVOLUTION_DLL BaseCorrelationLossCalculator: public IFixedTrancheLossCalculator{
public:
    virtual ~BaseCorrelationLossCalculator();

    /** Constructor - takes in the loss calculator which does the unadjusted  
        loss together with the lower and upper beta name and time dependent
        overrides. The lowerBetas[0] and upperBetas[0] values aren't actually
        used - the issue is that the timeline has today as the first date which
        makes the first value irrelevant. */
    BaseCorrelationLossCalculator(
        CreditMetricsLossCalculatorBaseSP originalLossCalculator, /* (I) */
		double							  baseStrike,
        double                            lowerStrike,
        double                            upperStrike,
        DoubleArrayArraySP                lowerBetas, // [timepoint][name index]
        DoubleArrayArraySP                upperBetas, // [timepoint][name index]
        bool                              authoriseNegativeEL,
        bool                              useExpectedLossRatio);

    /** Calculate the expected loss for specified timepoint and strikes using
        'base correlation' methodology */
    virtual void loss(
        int     timePoint,        // (I) do the calculation for this timepoint
        double& loss,             /* (O) tranche loss amt  */
        double& lossCond) const;  /* (O) tranche loss amt cond on cpty
                                     surviving */

    /** store calculator results */
    void storeResults(Results* result, Control* control) const {}

private:
    /// fields ////
    CreditMetricsLossCalculatorBaseSP originalLossCalculator;
	double						      baseStrike;
    double                            lowerStrike;
    double                            upperStrike;
    DoubleArrayArraySP                lowerBetas; // [timepoint][name index]
    DoubleArrayArraySP                upperBetas; // [timepoint][name index]
    bool                              authoriseNegativeEL;
    bool                              useExpectedLossRatio;
};

typedef smartConstPtr<
    BaseCorrelationLossCalculator> BaseCorrelationLossCalculatorConstSP;

DRLIB_END_NAMESPACE
#endif
