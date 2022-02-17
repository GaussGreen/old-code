//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : SRMSwaptionPricer.hpp
//
//   Description : Helper for SRM - used for calibration against a swaption
//
//   Author      : Mark A Robson
//
//   Date        : 11 June 2004
//
//
//----------------------------------------------------------------------------

#ifndef EDR_SRMSWAPTIONPRICER_HPP
#define EDR_SRMSWAPTIONPRICER_HPP
//#include "edginc/SRMRatesUtil.hpp"
#include "edginc/SVGenIRSwap.hpp"
#include "edginc/SVGenDiscFactor.hpp"


DRLIB_BEGIN_NAMESPACE

/** For pricing [call?] swaptions for use in ICE */
class MCARLO_DLL SRMSwaptionPricer :  public virtual VirtualDestructorBase
{
public:
    /** Ask for state variables */
    void collectStateVars(IStateVariableCollectorSP svCollector) const;
    /** Get hold of the state variables from the generator */
    void pathGenUpdated(IStateVariableGen::IStateGen* newPathGen);

    //// returns price of current path
    double payoff();
    //// simple constructor
    SRMSwaptionPricer(
        double                    strike,
        SVGenIRSwapSP                swapGen,
        SVGenDiscFactorSP            dfGen);

private:
    double                    strike;
    SVGenIRSwap::IStateVarSP     swapSV;  // swap state variable
    SVDiscFactorSP              dfSV;    // df state variable
    SVGenIRSwapSP                swapGen;
    SVGenDiscFactorSP            dfGen;
};
typedef smartPtr<SRMSwaptionPricer> SRMSwaptionPricerSP;


DRLIB_END_NAMESPACE
#endif // EDR_SRMSWAPTIONPRICER_HPP
