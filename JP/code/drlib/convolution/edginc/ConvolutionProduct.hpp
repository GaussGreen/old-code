//----------------------------------------------------------------------------
//
//   Group       : QR
//
//   Description : View onto instrument as required by ConvolutionEngine
//
//   Date        : 18th Nov 2005
//
//   Author      : Mark Robson
//
//
//----------------------------------------------------------------------------

#ifndef QR_CONVOLUTIONPRODUCT_HPP
#define QR_CONVOLUTIONPRODUCT_HPP

#include "edginc/IGeneralisedConvolutionProduct.hpp"
#include "edginc/Control_forward.hpp"
#include "edginc/Results_forward.hpp"
#include "edginc/DateTime.hpp"
#include "edginc/FORWARD_DECLARE_WRAPPER.hpp"

DRLIB_BEGIN_NAMESPACE

class OutputRequest;
class IConvolutionModel;
FORWARD_DECLARE(ConvolutionEngine);
FORWARD_DECLARE(FlatFwdZeroCurve); // defined below in this file
FORWARD_DECLARE(CounterPartyCredit);
FORWARD_DECLARE(SingleCreditAsset);
FORWARD_DECLARE(YieldCurve);
FORWARD_DECLARE(IFixedTrancheLossCalculator);
FORWARD_DECLARE(EffectiveCurve);
FORWARD_DECLARE(ICreditLossConfig);
FORWARD_DECLARE_REF_COUNT(IEffectiveCurveLossGen);
FORWARD_DECLARE(IForwardRatePricer);

/** Provides view onto Instrument required by model together with payoff()
    type functionality. Unclear whether we should return some sort of
    object per name which gives per name stuff */
class CONVOLUTION_DLL ConvolutionProduct : 
    public virtual IGeneralisedConvolutionProduct {

public:
    //// Possible values for creditChargeViewType in
    //// computeCreditCharge method below
    static const string CC_CONSERVATIVE;
    static const string CC_AGGRESSIVE;

    virtual ~ConvolutionProduct();

    /** Returns the last 'scheduled' pay date (ie ignoring credit
        events).  Introduced for the output request IND_CDS_PAR_SPREAD
        - not clear if this is right date to use */
    virtual DateTime lastPayDate() const = 0;

    /** Returns number of names contained within instrument */
    virtual int numNames() const = 0;

    /** Returns the representation of the underlying name corresponding to the
        specified index which must lie in range [0, numNames()-1] */
    virtual SingleCreditAssetConstSP nameAsset(int index) const = 0;

    /** Returns the notional for the name corresponding to the
        specified index which must lie in range [0, numNames()-1] */
    virtual double nameNotional(int index) const = 0;

    /** Returns the loss given default for the name corresponding to the
        specified index which must lie in range [0, numNames()-1] */
    virtual double nameLossGivenDefault(int index) const = 0;

    /** Returns the max of valueDate and protection start date
        for specified name (if any). The specified index must lie in
        range [0, numNames()-1]*/
    virtual const DateTime& nameProtectionStartDate(
        int index) const = 0;

    /** Returns the min of last date parameter and protection end date
        for specified name (if any). The specified index must lie in
        range [0, numNames()-1]. This if for name 'cutoff' */
    virtual const DateTime& nameProtectionEndDate(
        int index, const DateTime& lastDate) const = 0;

    /** Returns the recovery for the name corresponding to the
        specified index which must lie in range [0, numNames()-1]. Note that
        the recovery can  be overridden at the trade level hence this method,
        which reflects any trade level overrides, should be used rather than
        the recovery off the par spread curve */
    virtual double nameRecovery(int index) const = 0;

    /** Returns true if the name has defaulted [subject to
        nameProtectionEndDate()] where the name corresponds to the
        specified index which must lie in range [0, numNames()-1] */
    virtual bool nameDefaulted(int index) const = 0;

    /** Returns the 'beta' for the name corresponding to the specified
        index which must lie in range [0, numNames()-1] */
    virtual double nameBeta(int index) const = 0;

    /** Returns the beta between loss market and decretion market */
    virtual double lossDecretionBeta() const = 0;

    /** Returns the ranges for Loss Given Default for the name
        corresponding to the specified index which must lie in range
        [0, numNames()-1]. In the particular, the LGD of a name is given by
        nameNotional * MIN(lgdCap, MAX(lgdFloor, lgdNotional - recovery)) */
    virtual void nameLGDRanges(int     index,       /* (I) */
                               double& lgdNotional, /* (O) */
                               double& lgdFloor,    /* (O) */
                               double& lgdCap)      /* (O) */ const = 0;

    /** Returns a CounterPartyCredit object containing data required for
        pricing CCC. Note may return null if no counter party information is
        present. Going forward we might find that a CounterPartyCredit is too
        strong a type (ie may want something that promises to do less) */
    virtual CounterPartyCreditConstSP getCounterParty() const = 0;

    /** Returns the total notional of the portfolio */
    virtual double portfolioNotional() const = 0;

	/** Returns the total notional of the portfolio */
    virtual double portfolioLongNotional() const = 0;

	/** Returns the total notional of the portfolio */
    virtual double portfolioShortNotional() const = 0;

    /** Returns the sum of the past losses of the portfolio which eat
     * into the total notional returned by portfolioNotional() */
    virtual double portfolioLoss() const = 0;

    /** Returns whether the instument requires notional to be recovered
        from the top of the portfolio
        NB may be later overridden by the model if it is deemed to be unnecessary */
    virtual bool recoverNotional() const = 0;

    /** Returns the sum of the past recovered notional of the portfolio which eat
     * into the total notional returned by portfolioNotional() */
    virtual double portfolioRecoveredNotional() const = 0;

    /** Returns the sum of the past losses of the portfolio which eat
     * into the total notional returned by portfolioNotional() */
    virtual double portfolioLongLoss() const = 0;

    /** Returns the sum of the past losses of the portfolio which eat
     * into the total notional returned by portfolioNotional() */
    virtual double portfolioShortLoss() const = 0;

    /** Returns the sum of the past recovered notional of the portfolio which eat
     * into the total notional returned by portfolioNotional() */
    virtual double portfolioLongRecoveredNotional() const = 0;

	/** Returns the sum of the past recovered notional of the portfolio which eat
     * into the total notional returned by portfolioNotional() */
    virtual double portfolioShortRecoveredNotional() const = 0;

    // jlhp make private
    /** Returns the strikes of the tranche. These are relative to the values
        returned by portfolioRanges() */
    virtual void trancheStrikes(double& lowStrike,
                                double& highStrike) const = 0;

    /** Returns highStrike-lowStrike with highStrike and lowStrike obtained from
        trancheStrikes */
    virtual double trancheNotional() const;

    /** Returns the outstanding notional as seen from today */
    virtual double trancheOutstandingNotional(const bool recoveredNotional) const;

    /** Returns the outstanding risky notional of the tranche */
	virtual double trancheRiskyNotional(const bool recoverNotional) const;

    /** Product specific 'output' class. Useful because one will be created
        for both normal pricing and when calculating CCC */
    class CONVOLUTION_DLL Output: public virtual VirtualDestructorBase {
    public:
        virtual ~Output();
        Output();

        /** Returns the fair value (aka mark to market) given this
            Output object */
        virtual double price() const = 0;
    };

    DECLARE(Output); // sets up OutputSP

    /** Computes payoff for instrument given effectiveCurve.
        The forwardFactor is a relic of Kapital backward compatibility 
        and is a double by which all amounts should be scaled by and
        relects pv'ing things to IR spot date rather than today.
        Note when doing CCC it can't be read off the supplied
        discount curve. */
    virtual OutputSP payoff(
        EffectiveCurveSP     ctgLegEffectiveCurve,
        EffectiveCurveSP     feeLegEffectiveCurve,
        double               forwardFactor,
        IForwardRatePricerSP model) const = 0;  // see above

    /** Compute the credit charge using supplied results and credit charge
        view type. To do: sort out rep of creditChargeViewType - currently it
        can be one of the two strings defined above: CC_CONSERVATIVE or CC_
        AGGRESSIVE */
    virtual double computeCreditCharge(
        //FlatFwdZeroCurveSP   risklessDiscount,
        //FlatFwdZeroCurveSP   cptyCleanSpreadCurve,
        const EffectiveCurveSP cptyCurve,
        const string&          creditChargeViewType,
        OutputSP               output,
        OutputSP               cccOutput) const = 0;

    /** Populate Results object with product specific requests etc. Note
        cccOutput may be null if CCC not being computed. The engine
        will store fair value, CPTY_CREDIT_CHARGE, FV_MINUS_CCC */
    virtual void storeExtraResults(Results*  results,
                                   Control*  control,
                                   OutputSP  output,
                                   OutputSP  cccOutput,  // may be 0
                                   IForwardRatePricerSP model,
                                   const IConvolutionModel* convolutionModel) const = 0;
    
    /** Kapital backward compatibility mode - to be removed eventually.
        Returns whether we pv to cfCutOffDate() */
    bool pvToSpot() const;
    
    /** Kapital backward compatibility mode - to be removed eventually.
        Ignore cashflows as well as protection etc before this date. */
    DateTime cfCutOffDate() const;

    /** Returns today aka value date. Note not yield curve spot date */
    DateTime getToday() const;

    /** Returns the [riskless] curve for discounting (could remove this if
        we had better methods on ICDSParSpreads - only used for calculating
        durations) */
    YieldCurveConstSP getDiscount() const;

    /** Utility method: calculates counter party survival probabilities along
        timeline. Resizes counterPartyProb as required. */
    void computeCounterPartySurvivalProb(
        const DateTimeArray& timeline,          // (I)
        DoubleArray&         counterPartyProb) const; // (O) [time index]

    /** Invoked by model to compute the price. Note not virtual (although
        could be) as don't want (or rather anticipate) products to override */
    void price(const IConvolutionModel* convolutionModel,
               Control*                 control, 
               Results*                 results,
               IForwardRatePricerSP     model) const;

    /** Computes the effective curve along the specified timeline
        storing the result in effectiveCurve. If computeCondCurve is
        true then an effective curve conditional on the counterparty
        not defaulting should be returned too (in this case
        ConvolutionProduct::getCounterParty() will not return null).
        The timeline supplied here is the same as that generated by 
        generateTimeline, which will have been called first. 
        The DoubleArraySP will initially be null. */
    void calculateEffectiveCurve(
    IFixedTrancheLossCalculatorConstSP  lossCalculator,              /* (I) */
    IFixedTrancheLossCalculatorConstSP  recoveredNotionalCalculator, /* (I, may be null) */
    const DateTimeArray&                timeline,                    /* (I) */
    DoubleArraySP&                      ctgEffectiveCurve,           /* (O) */
    DoubleArraySP&                      feeEffectiveCurve,           /* (O) */
    bool                                computeCondCurve,            /* (I) */
    DoubleArraySP&                      ctgEffectiveCurveCond,       /* (O) */ 
    DoubleArraySP&                      feeEffectiveCurveCond)       /* (O) */ const;

protected:
    /** Constructor for derived types - just takes references */
    ConvolutionProduct(ConvolutionEngineConstSP convolutionEngine,
                       YieldCurveConstSP        discount);

private:
    ConvolutionProduct(const ConvolutionProduct& rhs);
    ConvolutionProduct& operator=(const ConvolutionProduct& rhs);
    
    void storeExpectedLossCurve(
        const EffectiveCurveSP effectiveCurve,
        OutputRequest*         request,
        Results*               results) const;

    void storeIndicativeSpreads(OutputRequest* indRequest,
                                OutputRequest* parRequest,
                                Results*       results) const;

    //FlatFwdZeroCurveSP createCptyCleanSpreadCurve() const;
    double calculateForwardFactor(YieldCurveConstSP discountCurve) const;

    void storePrimaryResults(const EffectiveCurveSP ctgEffectiveCurve,
                             const EffectiveCurveSP feeEffectiveCurve,
                             Results*             results,
                             Control*             control,
                             OutputSP             output,
                             OutputRequest*       cccRequest,
                             bool                 doCCC,      
                             double               ccc) const;

    void storeCalculatorResults(
        Results* results, Control* control, 
        IFixedTrancheLossCalculatorConstSP calculator) const;

    /// fields ////
    ConvolutionEngineConstSP convolutionEngine;
    YieldCurveConstSP        discount;
};

DRLIB_END_NAMESPACE

#endif
