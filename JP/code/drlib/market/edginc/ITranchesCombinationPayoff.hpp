//----------------------------------------------------------------------------
//
//   Group       : CH Quantitative Research
//
//   Date        : 16-Oct-2006
//
//----------------------------------------------------------------------------

#ifndef QLIB_ITRANCHESCOMBINATIONPAYOFF_HPP
#define QLIB_ITRANCHESCOMBINATIONPAYOFF_HPP

#include "edginc/DateTime.hpp"

DRLIB_BEGIN_NAMESPACE

/**
 * Interface for credit payoffs capable to decompose themselves into
 * a (possibly time dependent) linear combination of base tranche payoffs so that:
 * 
 * ExpectedLoss(T) = ExpectedLossOffset(T) + sum(i, baseStrikesWeights_i(T) * ExpectedLoss(baseStrikes_i(T))
 * 
 * where ExpectedLoss(baseStrikes_i(T)) is the expected loss of the base tranche (0, baseStrikes_i(T))
 *
 * */
class MARKET_DLL ITranchesCombinationPayoff
{
public:

    /**
     * Decomposition method.
     * Expects "baseStrikes" to be sorted in increasing order.
     * Each element of index i in "baseStrikesWeights" should correspond to an
     * element of index i in "baseStrikes". 
     * */
    virtual void linearDecomposition(
        const DateTime& time,
        DoubleArray& baseStrikes,        /* output */
        DoubleArray& baseStrikesWeights, /* output */
        double& expectedLossOffset       /* output */) const = 0;
};

DRLIB_END_NAMESPACE

#endif /*QLIB_ITRANCHESCOMBINATIONPAYOFF_HPP*/
