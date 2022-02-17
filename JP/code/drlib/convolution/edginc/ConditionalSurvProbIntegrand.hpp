//----------------------------------------------------------------------------
//
// Group       : Credit Hybrids QR
//
// Description : Computes the conditional survival probability at a given 
//               timepoint. Essentially multiplies the "inner names" 
//               conditional survival probabilities.
//
// Date        : Oct 2006
//
//----------------------------------------------------------------------------

#ifndef QLIB_CONDITIONALSURVPROBINTEGRAND_HPP
#define QLIB_CONDITIONALSURVPROBINTEGRAND_HPP

#include "edginc/ICondLossDistributionsGen.hpp"
#include "edginc/IConvolutor.hpp"
#include "edginc/IConditionalDefaultsModel.hpp"

DRLIB_BEGIN_NAMESPACE

/** Computes the  conditional survival probability at a given timepoint.
    Essentially calls the "convolution algorithm" on "inner names" 
    conditional survival probabilities. */
class ConditionalSurvProbIntegrand : public virtual FunctionNDDouble {

public:
    // destructor
    virtual ~ConditionalSurvProbIntegrand();
    
    // constructor
    ConditionalSurvProbIntegrand(
        ICondLossDistributionsGenKeyArrayConstSP keys,
        int mfDim,
        const DateTime& time,
        const ICreditLossConfig* lossConfig);
    
    // Function
    virtual double operator()(
        const CDoubleArray&  marketFactor) const;

private:
    // Point in the timeline
    DateTime time;

    // Array of keys (one per name)
    ICondLossDistributionsGenKeyArrayConstSP keys;

    const ICreditLossConfig* lossConfig;
};

DRLIB_END_NAMESPACE

#endif
