//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : SRMICE.hpp
//
//   Description : Code responsible for make ICE adjustments to diffusion
//
//   Date        : 27 May 2004
//
//
//----------------------------------------------------------------------------
#ifndef SRMICE_HPP
#define SRMICE_HPP

#include "edginc/SRMSwaption.hpp"
#include "edginc/QMCRatesDiffuse.hpp"

DRLIB_BEGIN_NAMESPACE

// little class for storing ICE data per currency
class MCARLO_DLL ICE{
public:
    ICE(SRMRatesUtilSP srmRatesUtil, QMCRatesDiffuseSP srmRatesDiffuse);

    /** Ask for state variables */
    void collectStateVars(IStateVariableCollectorSP svCollector) const;

    /** Tell IMCPrices that pathGen has been updated ie get hold of the
        state variables from the generator */
    void pathGenUpdated(IStateVariableGen::IStateGen* newPathGen);

    /** recalibrate the supplied SRMRates using the results obtained from 
        pricing the swaptions */
    void recalibrate(int pathIdx);
    void updateAllPrices(); // calls updatePriceWithPayoff on each pricer

private:
    void updatePriceWithPayoff(size_t pricerId); // adds payoff to the price

    vector<SRMSwaptionPricerSP>   swaptionPricers; // prices  swaptions
    vector<double>                prices; // sum so far of swaption prices
    vector<SRMSwaptionSP>         swaptions; // the swaptions we're pricing
    SRMRatesHJMUtilSP             srmRatesHJMUtil;
    QMCRatesDiffuseSP             srmRatesDiffuse;
};

DRLIB_END_NAMESPACE
#endif // SRMICE_HPP

