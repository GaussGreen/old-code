//----------------------------------------------------------------------------
//
//   Group       : QR&D
//
//   Filename    : SVQmcTrancheLoss.cpp
//
//   Description : Implementing derived state variable for PortfolioLoss object 
//
//
//----------------------------------------------------------------------------


#include "edginc/config.hpp"
#include "edginc/SVQmcTrancheLoss.hpp"
#include "edginc/DiscreteDistribution.hpp"
#include <numeric>

DRLIB_BEGIN_NAMESPACE

SVQmcTrancheLoss::SVQmcTrancheLoss(
								   const DoubleArraySP _notionals,
								   const DateTimeArraySP _portfolioLossDates,
								   const DoubleArraySP _lowerPct,
								   const DoubleArraySP _upperPct)
								   :
notionals(_notionals),
portfolioLossDates(_portfolioLossDates),
lowerLoss(_lowerPct->size()),
rangeOfLoss(_upperPct->size()),
totalNotional(accumulate(_notionals->begin(), _notionals->end(), 0.0)),
portfolioLoss(_portfolioLossDates->size(), 0.0),
trancheLossConfig(CreditTrancheLossConfigArraySP(new CreditTrancheLossConfigArray(_lowerPct->size())))
{
	for(int i=0; i<_lowerPct->size(); ++i)
	{
		lowerLoss[i] = totalNotional*(*_lowerPct)[i];
		rangeOfLoss[i]= totalNotional*(*_upperPct)[i] - lowerLoss[i];
		(*trancheLossConfig)[i] = CreditTrancheLossConfigSP(new CreditTrancheLossConfig((*_lowerPct)[i],(*_upperPct)[i]));
	}

	// NOTE: input checking should already be done by the SV Gen
}


void SVQmcTrancheLoss::prepare(bool /*mm*/ ) // parameter not used here, as no mm exists for this SV
{ 
	trancheLoss.reset(new CDoubleMatrix(portfolioLossDates->size(), lowerLoss.size())); 
}

SVQmcTrancheLossFullMC::SVQmcTrancheLossFullMC(
	const vector<SVDateOfDefaultSP>& _ddSVs,
	const DoubleArraySP _notionals,
	const DateTimeArraySP _portfolioLossDates,
	const DoubleArraySP _lowerPct,
	const DoubleArraySP _upperPct) 
	:
SVQmcTrancheLoss(_notionals, _portfolioLossDates, _lowerPct, _upperPct), 
ddSVs(_ddSVs)
{}



CDoubleMatrixConstSP SVQmcTrancheLossFullMC::getTrancheLosses()
{
	const DateTime& lastDate = portfolioLossDates->back();

	fill(portfolioLoss.begin(), portfolioLoss.end(), 0.0);

	// Go thru each name and check for default
	for (size_t i = 0; i < ddSVs.size(); ++i)
	{
		const DateTime& dDate = ddSVs[i]->getDateOfDefault();

		if (dDate > lastDate) // default happened after the interesting interval
			continue;

		double R = ddSVs[i]->getRecoveryRateAtDefault();
		double loss = (1-R)*(*notionals)[i];

		// If default before the last observation date, this default had 
		// some contribution to the portfolio losses. We start counting from the end
		// as we can stop as soon as we reach the dDate
		for(int t=portfolioLossDates->size() - 1; t>=0 && dDate >= (*portfolioLossDates)[t]; --t)
			portfolioLoss[t] += loss;
	}

	// now fill in the tranches:
	for(int t=0; t<portfolioLossDates->size(); ++t)
		for(int i=0; i<lowerLoss.size(); ++i)
			if ( portfolioLoss[t] > lowerLoss[i]) 
				(*trancheLoss)[t][i] = min(portfolioLoss[t] - lowerLoss[i], rangeOfLoss[i]);

	return trancheLoss;
}


SVQmcTrancheLossFastMC::SVQmcTrancheLossFastMC(
	const vector<SVSurvivalDiscFactorSP>& _sdfSVs,
	const DoubleArraySP _notionals,
	const DateTimeArraySP _portfolioLossDates,
	const DoubleArraySP _lowerPct,
	const DoubleArraySP _upperPct,
	const IConvolutorSP _convolutor)
	:
SVQmcTrancheLoss(_notionals, _portfolioLossDates, _lowerPct, _upperPct), 
sdfSVs(_sdfSVs),
convolutor(_convolutor)
{}


void SVQmcTrancheLossFastMC::prepare(bool /*mm*/ ) // parameter not used here, as no mm exists for this SV
{ 
	SVQmcTrancheLoss::prepare(false);

	// only if necessary -- if not, this function can be removed later
}

CDoubleMatrixConstSP SVQmcTrancheLossFastMC::getTrancheLosses()
{
	IDistribution1DArraySP distributions = IDistribution1DArraySP (new IDistribution1DArray(portfolioLoss.size()));
	for (int tr=0; tr<trancheLossConfig->size(); ++tr)
	{
		for (int t=0; t<portfolioLossDates->size(); ++t)
		{
			for (size_t i = 0; i < sdfSVs.size(); ++i)
			{
				double lossGivenDefault = (*notionals)[i] * ( 1.0 - sdfSVs[i]->getRecoveryRate(t));
				(*distributions)[i] = DiscreteDistributionSP (
								new DiscreteDistribution(lossGivenDefault, 1.0 - sdfSVs[i]->getSDF(t)));
			}
			for (int i=0; i<trancheLossConfig->size(); ++i)
				(*trancheLoss)[t][i] = convolutor->convoluteAndIntegrate(distributions,(*trancheLossConfig)[tr],DateTime());
		}
	}
	return trancheLoss;
}





DRLIB_END_NAMESPACE

