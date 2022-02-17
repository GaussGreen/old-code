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

#include "edginc/config.hpp"
#include "edginc/ConditionalSurvProbIntegrand.hpp"
#include "edginc/DoubleArrayMarketFactor.hpp"


DRLIB_BEGIN_NAMESPACE

// destructor
ConditionalSurvProbIntegrand::~ConditionalSurvProbIntegrand() 
{}

// constructor
ConditionalSurvProbIntegrand::ConditionalSurvProbIntegrand(
        ICondLossDistributionsGenKeyArrayConstSP keys,
        int mfDim,
        const DateTime& time,
        const ICreditLossConfig* lossConfig) :
    FunctionNDDouble(mfDim,
                     *InfiniteRange::createInfiniteRangeArray(mfDim)),
    time(time),
    keys(keys),
    lossConfig(lossConfig)
{}

// function
double ConditionalSurvProbIntegrand::operator()(
    const CDoubleArray&  marketFactor) const
{
    IMarketFactorValueSP mfValue(new DoubleArrayMarketFactor(marketFactor));

    double survivalProb = 1.0;
    int nbNames = !keys ? 0 : keys->size();
    for (int i=0; i < nbNames; ++i) {
        survivalProb *= (*keys)[i]->conditionalSurvProb(mfValue);
	}
    return survivalProb;
}

DRLIB_END_NAMESPACE
