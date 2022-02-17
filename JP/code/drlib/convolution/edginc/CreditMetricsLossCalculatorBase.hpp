//----------------------------------------------------------------------------
//
//   Group       : QR
//
//   Author      : Mark A Robson
//
//   Date        : 23 Dec 2005
//
//----------------------------------------------------------------------------

#ifndef QR_CREDITMETRICSLOSSCALCULATORBASE_HPP
#define QR_CREDITMETRICSLOSSCALCULATORBASE_HPP

#include "edginc/TrancheLossCalculator.hpp"
#include "edginc/DateTime.hpp"
#include "edginc/FORWARD_DECLARE.hpp"
#include "edginc/Control_forward.hpp"
#include "edginc/Results_forward.hpp"

DRLIB_BEGIN_NAMESPACE

class ConvolutionProduct;
class ITrancheLossCalculatorLegacy;
FORWARD_DECLARE(RflOnlyParameters);
FORWARD_DECLARE(SingleCreditAsset);
FORWARD_DECLARE(CreditMetricsLossCalculatorBase);
FORWARD_DECLARE(CreditTrancheLossConfig);
FORWARD_DECLARE(CounterPartyCredit);
typedef refCountPtr<
    ITrancheLossCalculatorLegacy> ITrancheLossCalculatorLegacySP;

/** Base class to ease implementations of ITrancheLossCalculator */
class CONVOLUTION_DLL CreditMetricsLossCalculatorBase:
    public virtual ITrancheLossCalculator
{
public:
    virtual ~CreditMetricsLossCalculatorBase();

    /** Constructor - takes in full timeline to allow for optimisations. If
        computeCondCurve if false then lossConditional in loss() below will be
        set to zero. */
    CreditMetricsLossCalculatorBase(
        const DateTimeArray&           timeline, /* (I) */
        CreditTrancheLossConfigConstSP tranche,  /* (I) */
        CounterPartyCreditConstSP      cpty);    /* (I) */

    /** Calculate the expected loss for specified timepoint and strikes. 
        Note if repeated calculations at the same timepoint are required then
        IKey::calc should be used instead. This implementation just does
        "lossKey()->loss(...)" */
    virtual void loss(
        int     timePoint,        // (I) do the calculation for this timepoint
        double  k1,               /* (I) lower strike      */
        double  k2,               /* (I) upper strike      */
        double& loss,             /* (O) tranche loss amt  */
        double& lossCond) const;  /* (O) tranche loss amt cond on cpty
                                     surviving */

    /** store calculator results */
    virtual void storeResults(Results* result, Control* control) const;

    /** Returns a key used to optimise repeated calculations of losses
        at the same timepoint. Implementation here just routes through
        lossKey below */
    virtual IKey* lossKey(
        int timePoint) const; // (I) do the calculation for this timepoint

    /** Same as above but allows the betas used to be overridden. The
        length of the array must be the same as the number of names or
        be empty (=> no beta overrides). This is used by Base Correlation */
    virtual IKey* lossKey(
        int                timePoint,   /* (I) do the calculation for 
                                           this timepoint */
        const DoubleArray& betaOverride) const = 0; // (I) using these betas

    /** Utility method to build a IFixedTrancheLossCalculator using a
        CreditMetricsLossCalculatorBase where you can override the betas
        at each time point/for each name. */
    static IFixedTrancheLossCalculator* createFixedTrancheLossCalculator(
        double                                 lowerStrike,  
        double                                 upperStrike,
        DoubleArrayArraySP                     betas,//override for [time][name]
        CreditMetricsLossCalculatorBaseConstSP lossCalculator);

    /** Retrieves RFL stlye engine params from Credit asset. Fails if they are
        not there or are of wrong type */
    static RflOnlyParametersConstSP getRFLEngineParameters(
        SingleCreditAssetConstSP asset);

protected:
    // implementation of IKey when using a CCMTrancheCalculatorLegacy
    class CONVOLUTION_DLL Key: public virtual ITrancheLossCalculator::IKey{
        ITrancheLossCalculatorLegacySP trancheCalc;
    public:
        ~Key();
        Key(ITrancheLossCalculatorLegacySP trancheCalc);

        /** Calculate the expected loss for specified strikes. */
        virtual void loss(
            double  k1,                   /* (I) lower strike      */
            double  k2,                   /* (I) upper strike      */
            double& loss,                 /* (O) tranche loss amt  */
            double& lossCond) const;      /* (O) tranche loss amt cond on cpty
                                             surviving */
    };
    /// fields //////////////
    bool             computeCondCurve; // avoid unnecessary calculations
    DoubleArrayArray survivalProb;
    DoubleArray      counterPartyProb;

private:
    class BetaOverride;

    static void computeCounterPartySurvivalProb(
        CounterPartyCreditConstSP counterParty, // (I)
        const DateTimeArray& timeline,          // (I)
        DoubleArray& counterPartyProb);         // (O)

};

typedef smartPtr<
    CreditMetricsLossCalculatorBase> CreditMetricsLossCalculatorBaseSP;

DRLIB_END_NAMESPACE
#endif
