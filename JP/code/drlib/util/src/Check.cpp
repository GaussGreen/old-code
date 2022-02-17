//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : Check.cpp
//
//   Description : Various useful checking functions etc
//
//   Date        : Oct 2002
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/Check.hpp"
#include "edginc/Maths.hpp"

DRLIB_BEGIN_NAMESPACE


/** Ensures that weights are positive, that weights.size() == numWeights
    and that the weights sum to 1.0 to within reasonable
    tolerance. The difference between 1.0 and the sum of the weights
    is then distributed back among the weights. */
void Check::percWeights(DoubleArray&       weights,
                        int                numToMatch,
                        string             desc) {
    static const string method("Check::percWeights");
    if (numToMatch != weights.size()){
        throw ModelException(method,
                             "Different number of " + desc + " ("+
                             Format::toString(numToMatch)+")"
                             " to weights ("+
                             Format::toString(weights.size())+")");
    }
    double weightSum = 0.0;
    for (int i = 0; i < numToMatch; i++){
        Maths::checkNonNegative(weights[i], "Weight number "+ Format::toString(i+1));
        weightSum += weights[i];
    }
    // check weights add up to 100% (or almost anyway)
    if (!Maths::areEqualWithinTol(weightSum, 1.0,  0.000000001)){
        throw ModelException(method, "Total % weights ("+
                             Format::toString("%.12f", weightSum)+") must "
                             "equal 1.0");
    }
    if (!Maths::equals(weightSum, 1.0)){
        // spread remainder among assets
        for (int i =0; i < numToMatch; i++){
            weights[i] += (1.0 - weightSum)/numToMatch;
        }
    }
}

DRLIB_END_NAMESPACE
