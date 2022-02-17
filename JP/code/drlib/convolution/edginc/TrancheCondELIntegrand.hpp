//----------------------------------------------------------------------------
//
//   Group       : CH Quantitative Research
//
//   Date        : 09-Oct-2006
//
//----------------------------------------------------------------------------

#ifndef QLIB_TRANCHECONDELINTEGRAND_HPP
#define QLIB_TRANCHECONDELINTEGRAND_HPP

#include "edginc/ICondLossDistributionsGen.hpp"
#include "edginc/IConvolutor.hpp"
#include "edginc/IConditionalDefaultsModel.hpp"

DRLIB_BEGIN_NAMESPACE

/**
 * Computes the tranche conditional expected loss at a given timepoint.
 * Essentially calls the convolution algorithm on "inner names" conditional loss distributions.
 * */
class TrancheCondELIntegrand :
    public virtual FunctionNDDouble
{
public:
    
    // destructor
    virtual ~TrancheCondELIntegrand();
    
    // constructor
    TrancheCondELIntegrand(
        ICondLossDistributionsGenKeyArrayConstSP lossKeys,
        int mfDim,
        DateTime time,
        const ICreditLossConfig* lossConfig,
        const IConvolutor* convolutor);
        
    // function
    virtual double operator()(
        const CDoubleArray&  marketFactor) const;

private:
    // Point in the timeline
    DateTime time;

    // Array of loss keys (one per name)
    ICondLossDistributionsGenKeyArrayConstSP lossKeys;

    const ICreditLossConfig* lossConfig;
    const IConvolutor* convolutor;
    
    // Local variable - avoid allocating a new array for each value of the market factor
    IDistribution1DArraySP distributionsToConvolute;
};

DRLIB_END_NAMESPACE

#endif /*QLIB_TRANCHECONDELINTEGRAND_HPP*/
