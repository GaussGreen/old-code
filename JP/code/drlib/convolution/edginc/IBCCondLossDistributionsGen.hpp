//----------------------------------------------------------------------------
//
//   Group       : CH Quantitative Research
//
//   Date        : 18-Aug-2006
//
//----------------------------------------------------------------------------

#ifndef QLIB_IBCCONDLOSSDISTRIBUTIONSGEN_HPP
#define QLIB_IBCCONDLOSSDISTRIBUTIONSGEN_HPP

#include "edginc/ICondLossDistributionsGen.hpp"
#include "edginc/IDistribution1D.hpp"
#include "edginc/DateTime.hpp"
#include "edginc/IMarketFactorValue.hpp"
#include "edginc/IndexWeights.hpp"
#include "edginc/IBCSkewCalculator.hpp"

DRLIB_BEGIN_NAMESPACE

/**
 * Specialisation of ICreditLossGen responsible for creating loss distributions
 * conditional on some "market factor" value and some strikes.
 * Used for "Base Correlation" methodology.
 * */
class CONVOLUTION_DLL IBCCondLossDistributionsGen: public virtual ICreditLossGen {
public:

    virtual ~IBCCondLossDistributionsGen();

    /**
     * Mechanism used to pass the "single name - index" mapping to
     * the tranche level
     * */
     virtual IndexWeightsConstSP getIndexWeights() const = 0;

    /**
     * Mechanism used to pass the "single name" notional to
     * the tranche level
     * */
     virtual double getNotional() const = 0;

    /**
     * Mechanism used to pass the "single name" historical beta to
     * the tranche level
     * */
     virtual double getHistoricalBeta() const = 0;
     
    /**
     * Mechanism used to pass the "single name" (unconditional)
     * expected loss to the tranche level
     * */
     virtual double getExpectedLoss(const DateTime& time) const = 0;
     
    /** Mechanism used to pass the "single name" default status to the tranche level */
     virtual bool hasDefaulted() const = 0;

    /** Mechanism used to pass the "single name" name (!) to the tranche level */
     virtual string getName() const = 0;

    /**
     * Mechanism used to pass the "single name" implied par spreads and 
     * durations to the tranche level.
     * - Needed for "wide spread" with duration weighted spread ratio -
     * */
    virtual void impliedParSpreadsAndDurations(
        const YieldCurveConstSP discount,
        const DateTimeArray& dates,
        DoubleArray& impliedSpreads,       /* (Output) */
        DoubleArray& durations) const = 0; /* (Output) */

    /**
     * Give a chance to do some market factor independent
     * initialisation (eg: computing clean spreads and default 
     * probabilities, calibrating thresholds...) and improve
     * performance.
     * Creates an IKey that will be used to compute loss distributions conditional
     * on some market factor value.
     * 
     * In the case of base correlation, we need to pass all the
     * relevant parameters to adjust the historical betas.
     * 
     * */
    virtual ICondLossDistributionsGenKeySP initialise(
        const DateTime& time,
        double adjustedStrike,
        const IBCSkewCalculator& bcSkewCalculator) const = 0;

protected:
    IBCCondLossDistributionsGen();

private:    
    IBCCondLossDistributionsGen(const IBCCondLossDistributionsGen& rhs); // don't use
    IBCCondLossDistributionsGen& operator=(const IBCCondLossDistributionsGen& rhs); // don't use
};

DECLARE_REF_COUNT(IBCCondLossDistributionsGen);

DRLIB_END_NAMESPACE

#endif /*QLIB_IBCCONDLOSSDISTRIBUTIONSGEN_HPP*/
