//----------------------------------------------------------------------------
//
//   Group       : QR&D
//
//   Filename    : SVQmcPortfolioLoss.cpp
//
//   Description : Implementing derived state variable for PortfolioLoss object 
//
//
//----------------------------------------------------------------------------


#include "edginc/config.hpp"
#include "edginc/SVQmcPortfolioLoss.hpp"

DRLIB_BEGIN_NAMESPACE

void SVQmcPortfolioLossFullMC::calculateLossEvents() const
{
    for(size_t i=0; i<ddSVs.size(); ++i)
    {
        PortfolioNameSVSubResult sub(
            ddSVs[i]->getDateOfDefault(), 
            ddSVs[i]->getRecoveryRateAtDefault());
        (*mySubResults)[i] = sub;
    }

    SVQmcPortfolioLoss::calculateLossEvents();
}

void SVQmcPortfolioLossFullMC::setSubResults(const PortfolioNameSVSubResultVectorConstSP&  rhsSubResults) const
{
    throw ModelException("SVQmcPortfolioLossFullMC::setSubResults","This function cannot be called in QMC mode.");
}


void SVQmcPortfolioLossFastMC::calculateLossEvents() const
{
    // do this calculation for every point in the timeline

    creditLossConfigSVResults.clear();

    double lastloss = 0, lastaffected = 0;
	for (size_t idx = 0; idx < getTimeline()->size(); ++idx)
	{
		const DateTimeLite& date = getTimeline()->at(idx);

		if ( (date < protectionStartDate) || (maturityCutOff < date) )   //asset has not defaulted
			continue;

        double loss = 0.0, affected = 0.0;
	    for (size_t i = 0; i < sdfSVs.size(); ++i)
	    {
		    double recovery = sdfSVs[i]->getRecoveryRate(idx);

		    if (recoveryOverride.get())
			    recovery = recoveryOverride->doubleValue();

            double sdf = sdfSVs[i]->getSDF(idx);

            loss += getNotional()  * (1.0 - sdf) * ( 1.0 - recovery );
            affected += getNotional() * (1.0 - sdf);
        }
        // storing incremental effect on this date
		creditLossConfigSVResults.push_back(CreditLossConfigSVResult(
																date,
																loss-lastloss,
																affected-lastaffected,
    															idx));
        lastloss = loss;
        lastaffected = affected;
	}

}

void SVQmcPortfolioLossFastMC::setSubResults(const PortfolioNameSVSubResultVectorConstSP&  rhsSubResults) const
{
    throw ModelException("SVQmcPortfolioLossFastMC::setSubResults","This function cannot be called in QMC mode.");
}

//SVQmcPortfolioLoss::SVQmcPortfolioLoss(
//									   const DoubleArraySP _notionals,
//									   const DateTimeArraySP _portfolioLossDates)
//									   :
//notionals(_notionals),
//portfolioLossDates(_portfolioLossDates)
//{
//	// NOTE: input checking should already be done by the SV Gen
//}
//
//
//void SVQmcPortfolioLoss::prepare(bool /*mm*/ ) // parameter not used here, as no mm exists for this SV
//{ 
//	portfolioLoss = DoubleArraySP(new DoubleArray(portfolioLossDates->size(), 0.0)); 
//}
//
//SVQmcPortfolioLossFullMC::SVQmcPortfolioLossFullMC(
//	const vector<SVDateOfDefaultSP>& _ddSVs,
//	const DoubleArraySP _notionals,
//	const DateTimeArraySP _portfolioLossDates) 
//	:
//SVQmcPortfolioLoss(_notionals, _portfolioLossDates), 
//ddSVs(_ddSVs) 
//{}
//
//

//DoubleArrayConstSP SVQmcPortfolioLossFullMC::getPortfolioLoss()
//{
//	const DateTime& lastDate = portfolioLossDates->back();
//
//	fill(portfolioLoss->begin(), portfolioLoss->end(), 0.0);
//
//	// Go thru each name and check for default
//	for (size_t i = 0; i < ddSVs.size(); ++i)
//	{
//		const DateTime& dDate = ddSVs[i]->getDateOfDefault();
//
//		if (dDate > lastDate) // default happened after the interesting interval
//			continue;
//
//		double R = ddSVs[i]->getRecoveryRateAtDefault();
//		double loss = (1-R)*(*notionals)[i];
//
//		// If default before the last observation date, this default had 
//		// some contribution to the portfolio losses. We start counting from the end
//		// as we can stop as soon as we reach the dDate
//		for(int t=portfolioLossDates->size() - 1; t>=0 && dDate >= (*portfolioLossDates)[t]; --t)
//			(*portfolioLoss)[t] += loss;
//	}
//
//	return portfolioLoss;
//}
//
//
//SVQmcPortfolioLossFastMC::SVQmcPortfolioLossFastMC(
//	const vector<SVSurvivalDiscFactorSP>& _sdfSVs,
//	const DoubleArraySP _notionals,
//	const DateTimeArraySP _portfolioLossDates)
//	:
//SVQmcPortfolioLoss(_notionals, _portfolioLossDates), 
//sdfSVs(_sdfSVs) 
//{}
//
//
//void SVQmcPortfolioLossFastMC::prepare(bool /*mm*/ ) // parameter not used here, as no mm exists for this SV
//{ 
//	SVQmcPortfolioLoss::prepare(false);
//
//	// only if necessary -- if not, this function can be removed later
//}
//
//DoubleArrayConstSP SVQmcPortfolioLossFastMC::getPortfolioLoss()
//{
//	fill(portfolioLoss->begin(), portfolioLoss->end(), 0.0); // resetting all values to zero
//
//	for (size_t i = 0; i < sdfSVs.size(); ++i)
//	{
//		for (size_t idx = 0; idx < portfolioLoss->size(); ++idx)
//		{
//			double lossGivenDefault = (*notionals)[i] * ( 1.0 - sdfSVs[i]->getRecoveryRate(idx));
//			(*portfolioLoss)[idx] += lossGivenDefault*(1.0 - sdfSVs[i]->getSDF(idx));
//		}
//	}
//	// do the calculation of the portfolio loss out of the bunch of single name conditional
//	// survival probabilities
//	// i.e.   sdfSVs[i]->getSDF(idx); will get the conditional SP for the time point idx
//	//   (corresponding to the real date = portfolioLossDates[idx] )
//
//
//	// other given things are:
//	// notionals
//	// portfolioLossDates
//	// recovery rates - not in this code yet, but can be easily enabled (use as recovery[idx])
//
//	return portfolioLoss;
//}





DRLIB_END_NAMESPACE

