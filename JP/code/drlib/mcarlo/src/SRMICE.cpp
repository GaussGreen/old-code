//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : SRMICE.cpp
//
//   Description : Code responsible for make ICE adjustments to diffusion
//
//   Date        : 27 May 2004
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/SRMICE.hpp"

USING_DRLIB_NAMESPACE


ICE::ICE(SRMRatesUtilSP _srmRatesUtil, QMCRatesDiffuseSP _srmRatesDiffuse): 
    srmRatesDiffuse(_srmRatesDiffuse)
{
    static const string method("SRMICE::SRMICE");

    srmRatesHJMUtil = SRMRatesHJMUtilSP(dynamic_cast<SRMRatesHJMUtil*>(_srmRatesUtil.get()));
	// the model must be HJM
    if (srmRatesHJMUtil.get() == 0)
	{
        throw ModelException(method, "Failed to cast IR model to HJM");
	}

    const DateTimeArray& extendedTimeLine = 
        srmRatesHJMUtil->getExtendedTimeLine();
    const DateTimeArray& swapStartDates = srmRatesHJMUtil->getSwapStartDates();
    const DateTimeArray& swapMatDates = srmRatesHJMUtil->getSwapMatDates();
    const DateTime& lastSimDate = srmRatesHJMUtil->getSimDates().back();
    const DateTimeArray& swaptionExpiryDates = 
        srmRatesHJMUtil->getSwaptionExpiryDates();
    // while swaption expiry date is within simulation
    for (int i = 0; 
            i < swaptionExpiryDates.size() &&
                lastSimDate.isGreaterOrEqual(swaptionExpiryDates[i]);
            i++){
        const DateTime& prevSwaptionExpiry = (i == 0) ?
            srmRatesHJMUtil->getBaseDate(): swaptionExpiryDates[i-1];
        SRMSwaptionSP swaption(
            new SRMSwaption(
                *srmRatesHJMUtil, 
                swaptionExpiryDates[i],
                swapStartDates[i],
                swapMatDates[i],
                prevSwaptionExpiry,
                extendedTimeLine));
        swaptions.push_back(swaption);
        prices.push_back(0.0);
        SRMSwaptionPricerSP pricer(swaption->createPricer(*srmRatesHJMUtil));
        swaptionPricers.push_back(pricer);
    }
}

/** Ask for state variables */
void ICE::collectStateVars(
    IStateVariableCollectorSP svCollector) const{
    for (unsigned int i = 0; i < swaptionPricers.size(); i++){
        swaptionPricers[i]->collectStateVars(svCollector);
    }
}

/** Tell IMCPrices that pathGen has been updated ie get hold of the
    state variables from the generator */
void ICE::pathGenUpdated(IStateVariableGen::IStateGen* newPathGen)
{
    for (unsigned int i = 0; i < swaptionPricers.size(); i++){
        swaptionPricers[i]->pathGenUpdated(newPathGen);
    }
}

/** recalibrate the supplied SRMRates using the results obtained from 
    pricing the swaptions */
void ICE::recalibrate(int pathIdx)
{
    const double MAX_VOL_ERR =  0.15;
    // get hold of swaption vols
    const DoubleArray& swapVols = srmRatesHJMUtil->getSwaptionVols();
    // get hold of the vols we're recalibrating
    vector<double>& spotVol = srmRatesDiffuse->getOrigSVol();
    for (unsigned int i = 0; i < prices.size(); i++){
        double meanPrice = prices[i]/(pathIdx+1);
        double impVol = swaptions[i]->impliedVol(meanPrice);
        if (impVol < 0.0 || swapVols[i] == 0.0){
            // implied vol method failed, or swapVol is zero. So skip
        } else {
            double bias = Maths::square(impVol) -
                Maths::square(swapVols[i]);
            /* New target to be reached by the calibration */
            double target = sqrt(Maths::square(swapVols[i]) - bias);
            double simErr = fabs(target - swapVols[i])/swapVols[i];
            if (simErr > MAX_VOL_ERR) {
                // too far away - don't recalibrate
            } else {
                swaptions[i]->calibrate(*srmRatesHJMUtil, target, spotVol);
            }
        }
    }
    // then update the path generator
    srmRatesDiffuse->recalibrate(srmRatesHJMUtil);
}

void ICE::updatePriceWithPayoff(size_t pricerId)
{
    prices[pricerId] += swaptionPricers[pricerId]->payoff();
}        

void ICE::updateAllPrices()
{
    for(size_t i=0; i < swaptionPricers.size(); ++i)
        updatePriceWithPayoff(i);
}
