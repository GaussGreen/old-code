//----------------------------------------------------------------------------
//
//   Group       : CH Quantitative Research
//
//   Date        : 04-Sep-2006
//
//----------------------------------------------------------------------------

#ifndef QLIB_BCPORTFOLIONAMECONVOLUTIONLOSSGEN_HPP
#define QLIB_BCPORTFOLIONAMECONVOLUTIONLOSSGEN_HPP

#include "edginc/IModelConfigMapper.hpp"
#include "edginc/IBCCondLossDistributionsGen.hpp"
#include "edginc/IConditionalDefaultsModel.hpp"
#include "edginc/PortfolioName.hpp"

DRLIB_BEGIN_NAMESPACE

/**
 * Loss generator implementation for PortfolioName using
 * Base Correlation methodology.
 * Needs all relevant "loss generators" implementations since
 * the "outer model" may use any of them.
 * */
class CONVOLUTION_DLL BCPortfolioNameConvolutionLossGen:
    public virtual IBCCondLossDistributionsGen
{
public:

    virtual ~BCPortfolioNameConvolutionLossGen();

    BCPortfolioNameConvolutionLossGen(
        const PortfolioName* ptfName,
        const IConditionalDefaultsModel* condDefaultsModel);

    /**
     * Mechanism used to pass the "single name - index" mapping to
     * the tranche level
     * */
    virtual IndexWeightsConstSP getIndexWeights() const;

    /**
     * Mechanism used to pass the "single name" notional to
     * the tranche level
     * */
    virtual double getNotional() const;

    /**
     * Mechanism used to pass the "single name" historical beta to
     * the tranche level
     * */
    virtual double getHistoricalBeta() const;
     
     /**
     * Mechanism used to pass the "single name" (unconditional)
     * expected loss to the tranche level
     * */
    virtual double getExpectedLoss(const DateTime& time) const;
     
     /** Mechanism used to pass the "default status" to the tranche level */
    virtual bool hasDefaulted() const;

    /** Mechanism used to pass the loss config "single name" name (!) to the tranche level */
    virtual string getName() const;

    /**
     * Mechanism used to pass the "single name" implied par spreads and 
     * durations to the tranche level.
     * - Needed for "wide spread" with duration weighted spread ratio -
     * */
    virtual void impliedParSpreadsAndDurations(
        const YieldCurveConstSP discount,
        const DateTimeArray& dates,
        DoubleArray& impliedSpreads, /* (Output) */
        DoubleArray& durations) const;     /* (Output) */

    /** [Implements IBCCondLossDistributionsGen] */
    ICondLossDistributionsGenKeySP initialise(
        const DateTime& time,
        double adjustedStrike,
        const IBCSkewCalculator& bcSkewCalculator) const;
    
private:
    BCPortfolioNameConvolutionLossGen(const BCPortfolioNameConvolutionLossGen& rhs); // copy constructor - not defined
    BCPortfolioNameConvolutionLossGen& operator=(const BCPortfolioNameConvolutionLossGen& rhs); // operator= method - not defined

    // fields
    const PortfolioName* ptfName;
    const IConditionalDefaultsModel* condDefaultsModel;
    
    // Index Weights (BC Parameters)
    IndexWeightsConstSP indexWeights;
    
    // DefaultRates stored locally for performance
    DefaultRatesSP defaultRates;
};

DRLIB_END_NAMESPACE

#endif /*BCPORTFOLIONAMECONVOLUTIONLOSSGEN*/
