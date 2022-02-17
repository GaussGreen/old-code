//----------------------------------------------------------------------------
//
//   Group       : CH Quantitative Research
//
//   Date        : 03-Aug-2006
//
//----------------------------------------------------------------------------

#ifndef QLIB_ICONDLOSSDISTRIBUTIONSGEN_HPP
#define QLIB_ICONDLOSSDISTRIBUTIONSGEN_HPP

#include "edginc/ICreditLossGen.hpp"
#include "edginc/IDistribution1D.hpp"
#include "edginc/DateTime.hpp"
#include "edginc/IMarketFactorValue.hpp"

DRLIB_BEGIN_NAMESPACE

/**
 * Object responsible for actually computing a conditional loss
 * distribution for a given market factor value.
 * IKey are created (once) by ICondLossDistributionsGen::initialise and are
 * then used many times - once for each market factor value (typically 100 times).
 * 
 * Without this mechanism, we would have to use a cache to retrieve parameters
 * calibrated in ICondLossDistributionsGen::initialise() each time we call
 * conditionalLossDistribution, which is usually less efficient.
 * 
 * An example of implementation can be found in CreditMetricsDefaultModel.
 * */
class CONVOLUTION_DLL ICondLossDistributionsGenKey: public virtual IObject
{
public:
    static CClassConstSP const TYPE;

    /**
     * Computes a survival probability conditional on
     * a "market factor" value
     * */
    virtual IDistribution1DConstSP conditionalLossDistribution(
        IMarketFactorValueConstSP marketFactorValue) const = 0;

    /** Computes a survival probability conditional on
        a "market factor" value */
    virtual double conditionalSurvProb(
      IMarketFactorValueConstSP marketFactorValue) const { return 0.0; }
};

DECLARE(ICondLossDistributionsGenKey);

/**
 * Specialisation of ICreditLossGen responsible for creating loss distributions
 * conditional on some "market factor" value.
 * */
class CONVOLUTION_DLL ICondLossDistributionsGen: public virtual ICreditLossGen {
public:

    virtual ~ICondLossDistributionsGen();
    
    /**
     * Give a chance to do some market factor independent
     * initialisation (eg: computing clean spreads and default 
     * probabilities, calibrating thresholds...) and improve
     * performance.
     * Creates an IKey that will be used to compute loss distributions conditional
     * on some market factor value.
     * */
    virtual ICondLossDistributionsGenKeyArrayConstSP initialise(
        const DateTimeArray& timeline) const = 0;

protected:
    ICondLossDistributionsGen();

private:    
    ICondLossDistributionsGen(const ICondLossDistributionsGen& rhs); // don't use
    ICondLossDistributionsGen& operator=(const ICondLossDistributionsGen& rhs); // don't use
};

DECLARE_REF_COUNT(ICondLossDistributionsGen);

DRLIB_END_NAMESPACE

#endif /*QLIB_ICONDLOSSDISTRIBUTIONSGEN_HPP*/
