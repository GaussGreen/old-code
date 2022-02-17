#include "edginc/config.hpp"
#include "edginc/TrancheCondELIntegrand.hpp"
#include "edginc/DoubleArrayMarketFactor.hpp"
#include "edginc/ICreditLossConfig.hpp"

DRLIB_BEGIN_NAMESPACE

// destructor
TrancheCondELIntegrand::~TrancheCondELIntegrand() {}

// constructor
TrancheCondELIntegrand::TrancheCondELIntegrand(
    ICondLossDistributionsGenKeyArrayConstSP lossKeys,
    int mfDim,
    DateTime time,
    const ICreditLossConfig* lossConfig,
    const IConvolutor* convolutor) :
        FunctionNDDouble(
            mfDim,
            *InfiniteRange::createInfiniteRangeArray(mfDim)),
        time(time),
        lossKeys(lossKeys),
        lossConfig(lossConfig),
        convolutor(convolutor)
{
    distributionsToConvolute.reset(new IDistribution1DArray(lossKeys->size()));
}

// function
double TrancheCondELIntegrand::operator()(
    const CDoubleArray&  marketFactor) const
{
    IMarketFactorValueSP mfValue(new DoubleArrayMarketFactor(marketFactor));

    for (int i = 0; i < distributionsToConvolute->size(); ++i) {
        // Not very nice, ideally distributionsToConvolute should be
        // an array<IDistribution1DConstSP, IDistribution1D>.
        // Note that distributionsToConvolute is then passed as a
        // IDistribution1DArrayConstSP in convolutor->convoluteAndIntegrate
        (*distributionsToConvolute)[i] = IDistribution1DSP::constCast(
            (*lossKeys)[i]->conditionalLossDistribution(mfValue));
	}

	DateTime valueDate;
    // Convolutes & apply lossConfig payoff
    double expectedLoss = convolutor->convoluteAndIntegrate(
        distributionsToConvolute,
        ICreditLossConfigConstSP(lossConfig),
        time);
        
    // Multiply by market factor density
    return expectedLoss;
}

DRLIB_END_NAMESPACE
