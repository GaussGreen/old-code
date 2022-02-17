//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : SRMSwaptionPricer.cpp
//
//   Description : Helper for SRM - used for calibration against a swaption
//
//   Author      : Mark A Robson
//
//   Date        : 11 June 2004
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/SRMSwaptionPricer.hpp"
#include "edginc/Maths.hpp"
//#include "edginc/SwapTool.hpp"
//#include "edginc/SRMConstants.hpp"

DRLIB_BEGIN_NAMESPACE

/** Ask for state variables */
void SRMSwaptionPricer::collectStateVars(
    IStateVariableCollectorSP svCollector) const
{
    svCollector->append(swapGen.get()); // ask for a swap one
    svCollector->append(dfGen.get()); // and a DiscFactor one
}

/** Get hold of the state variables from the generator */
void SRMSwaptionPricer::pathGenUpdated(
    IStateVariableGen::IStateGen* newPathGen)
{
    swapSV = swapGen->getIRSwapSV(newPathGen);
    dfSV = dfGen->getSVDiscFactor(newPathGen);
}

double SRMSwaptionPricer::payoff()
{
    double parYield, annuity;
    // calculate rate for a par swap and what coupons are worth for a rate
    // of 1.0
    swapSV->parYield(parYield, annuity);
    double myPayoff = Maths::max(0.0, parYield - strike);
    // then multiple reduced rate by annuity
    myPayoff *= annuity;
    myPayoff *= dfSV->firstDF(); // and then discount to today
    return myPayoff;
}

//// simple constructor
SRMSwaptionPricer::SRMSwaptionPricer(
    double                    strike,
    SVGenIRSwapSP                swapGen, // swap state variable generator
    SVGenDiscFactorSP            dfGen):  // df state variable generator
        strike(strike), 
        swapGen(swapGen), 
        dfGen(dfGen){}


        
DRLIB_END_NAMESPACE
