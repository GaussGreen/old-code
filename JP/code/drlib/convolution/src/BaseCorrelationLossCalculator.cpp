//----------------------------------------------------------------------------
//
//   Group       : QR
//
//   Author      : Mark A Robson (from CDO.cpp/ccm2 lib)
//
//   Date        : 3 Jan 2006
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/BaseCorrelationLossCalculator.hpp"
#include "edginc/CreditMetricsLossCalculator.hpp"
#include "edginc/Maths.hpp"

DRLIB_BEGIN_NAMESPACE
BaseCorrelationLossCalculator::~BaseCorrelationLossCalculator(){}

/** Constructor - takes in the loss calculator which does the unadjusted  
    loss together with the lower and upper beta name and time dependent
    overrides */
BaseCorrelationLossCalculator::BaseCorrelationLossCalculator(
    CreditMetricsLossCalculatorBaseSP originalLossCalculator, /* (I) */
	double						      baseStrike,
    double                            lowerStrike,
    double                            upperStrike,
    DoubleArrayArraySP                lowerBetas, // [timepoint][name index]
    DoubleArrayArraySP                upperBetas, // [timepoint][name index]
    bool                              authoriseNegativeEL,
    bool                              useExpectedLossRatio):
    originalLossCalculator(originalLossCalculator),
	baseStrike(baseStrike),
    lowerStrike(lowerStrike), upperStrike(upperStrike),
    lowerBetas(lowerBetas), upperBetas(upperBetas),
    authoriseNegativeEL(authoriseNegativeEL),
    useExpectedLossRatio(useExpectedLossRatio){}

/** Calculate the expected loss for specified timepoint and strikes using
    'base correlation' methodology */
void BaseCorrelationLossCalculator::loss(
    int     timePoint,        // (I) do the calculation for this timepoint
    double& loss,                 /* (O) tranche loss amt  */
    double& lossCond) const{      /* (O) tranche loss amt cond on cpty
                                     surviving */
    const char routine[] = "BaseCorrelationLossCalculator::loss";
    try {

        double loss1 = 0.0;
        double loss1Cond = 0.0;
        double loss2 = 0.0;
        double loss2Cond= 0.0;

        // Avoid doing useless computation (with equity tranches for example...)
        if (baseStrike != lowerStrike)
        {
            // Create loKey using originalLossCalculator for this timepoint
            // but using the beta overrides
            ITrancheLossCalculator::IKeySP loKey(
                originalLossCalculator->lossKey(timePoint, 
                                                (*lowerBetas)[timePoint]));
    
            // Then compute losses
            loKey->loss(baseStrike, lowerStrike, loss1, loss1Cond);
        }

        // Avoid doing useless computation, just in case
        if (baseStrike != upperStrike)
        {
            // Create loKey using originalLossCalculator for this timepoint
            // but using the beta overrides
            ITrancheLossCalculator::IKeySP hiKey(
                originalLossCalculator->lossKey(timePoint, 
                                                (*upperBetas)[timePoint]));
    
            // Then compute losses
            hiKey->loss(baseStrike, upperStrike, loss2, loss2Cond);
        }
        
        // 3. combine results
        if (authoriseNegativeEL) {
            loss = loss2 - loss1;
            lossCond = loss2Cond - loss1Cond;
        } else {
            loss = Maths::max(loss2 - loss1, 0.0);
            lossCond = Maths::max(loss2Cond - loss1Cond, 0.0);
        }
    } catch (exception& e) {
        throw ModelException(e, routine);
    }
};

DRLIB_END_NAMESPACE
