//----------------------------------------------------------------------------
//
//   Group       : CH Quantitative Research
//
//   Description : An ICreditLossConfig used for CDO tranche products.
//
//   Date        : July 2006
//
//----------------------------------------------------------------------------


#ifndef QLIB_CREDITTRANCHELOSSCONFIG_HPP
#define QLIB_CREDITTRANCHELOSSCONFIG_HPP

#include "edginc/Atomic.hpp"
#include "edginc/Theta.hpp"
#include "edginc/CashFlow.hpp"
#include "edginc/MarketWrapper.hpp"
#include "edginc/ICreditLossConfig.hpp"
#include "edginc/RationalisedCreditEngineParameters.hpp"
#include "edginc/IEffectiveCurveLossModelConfig.hpp"
#include "edginc/FORWARD_DECLARE_WRAPPER.hpp"
#include "edginc/CDOParallelStrikeShift.hpp"
#include "edginc/CreditLossConfigMC.hpp"
#include "edginc/ITranchesCombinationPayoff.hpp"

DRLIB_BEGIN_NAMESPACE

/* Forward declare the classes refered to from this file if we are not
 * interested in their storage properties (ie, they are used through (smart)
 * pointers and therefore their include files are not required here - they
 * will be required in the .cpp file though)*/
FORWARD_DECLARE(IProtectionProvider);
FORWARD_DECLARE(GeneralisedCDO); // jlhp just here temporarily
FORWARD_DECLARE(CDOPortfolio); // jlhp just here temporarily
FORWARD_DECLARE(SingleCreditAsset);


class MARKET_DLL CreditTrancheLossConfig :
    public CObject,
    virtual public ICreditLossConfig,
    virtual public IEffectiveCurveLossModelConfig::IIntoLossGen,
    virtual public IRestorableWithRespectTo<CDOParallelStrike>,
    virtual public Theta::Shift,
	virtual public ICreditLossConfigSVGenMC,  // so that it supports MC
	virtual public ICreditLossConfigIndexedSVGenMC,
    virtual public ITranchesCombinationPayoff
{
    friend class GeneralisedCDO; // jlhp just temporarily

public:
    /** TYPE (for reflection) */
    static CClassConstSP const TYPE;

    virtual ~CreditTrancheLossConfig();

    ICreditLossConfigConstSP getPortfolio() const {return portfolio;}; // jlhp remove, only exists temporarily

    // Public constructor, with all fields passed in
    CreditTrancheLossConfig(const string& name,
                            const DateTime& valueDate,
                            const double lowStrike,
                            const double highStrike,
                            ICreditLossConfigSP portfolio);

    /** Constructor for tranches without portfolio - i.e. used as pure
        "low strike" - "high strike" payoff */
    CreditTrancheLossConfig(double lowStrike,
                            double highStrike);

    virtual IObject* clone() const;

    virtual void validatePop2Object();

    //++++++++++++++++++++++++++++++++++++++++
    //ICreditLossConfig methods
    //
    /** All instances have a 'name' which can be viewed as the name of
        the portfolio or the name of the single curve. */
    virtual string getName() const;

    /** Returns the historic losses on the contingent leg, due to the
        settlement of this name's default. In the presence of shorts
        the 'losses' can be negative.
        If there are no losses, it may return an empty (0 size)
        CtgLegLossPerDefaultArraySP or a "null" SP.
        Note: the cashflow dates are the CALCULATION DATES (as opposed to
        settlement dates) when the corresponding losses are determined. */
    virtual CtgLegLossPerDefaultArraySP historicContingentLegLosses(
         CIntConstSP triggerDelay,
         CIntConstSP defaultToCalculationDelay,
         const DateTime& lastTriggerDate,
         IBadDayAdjusterConstSP bda,
         const IProtectionProvider* const protect) const;

    /** Returns the notional reductions on the fee leg (due to losses and/or
        recovered notional) and the corresponding fee rebates if any. In the
        presence of shorts the 'losses' can be negative
        The accrualPeriods are potentially required to determine the actual
        date of the notional reductions. */
    virtual FeeLegReductionPerDefaultArraySP historicFeeLegReductions(
        CIntConstSP triggerDelay,
        CIntConstSP defaultToCalculationDelay,
        const double temporaryLossAmount,
        const DateTime& lastTriggerDate,
        AccrualPeriodArrayConstSP accrualPeriods,
        IBadDayAdjusterConstSP bda,
        const bool recoverNotional) const;

    /** Compute the rebate payments caused by historic fee leg reductions */
    virtual CashFlowArraySP historicRebatePayments(
        const IRebateCalculator* const   rebateCalc,
        FeeLegReductionPerDefaultArraySP trancheReductions,
        IForwardRatePricerSP             model,
        const bool                       recoverNotional) const;

    /** populate from market cache */
    virtual void getMarket(const IModel* model, const MarketData* market);

    /** Returns the engine parameters for this loss cfg (NB Model needs to
        have selected appropriate set during getMarket()). The parameter
        can specify what type of engine parameters are required. An
        exception is thrown if the relevant type is not found */
    virtual CreditEngineParametersConstSP getEngineParams(
        CClassConstSP engineParamsType) const;

    /** Returns the number of "inner loss configs" contained in this LossConfig.
        In this case cascade down the call to the portfolio. */
    virtual int numInnerLossConfigs() const;

    /** Returns the inner loss config */
    virtual ICreditLossConfigConstSP getInnerLossConfig() const;

    /** Returns the inner loss config number "index".
        In this case cascade down the call to the portfolio. */
    virtual ICreditLossConfigConstSP getInnerLossConfig(const int index) const;

    /** Returns the maximum (realistic) loss that may be produced by this
        loss config. Equivalent to lossGivenDefault for single names.
        "Realistic" in the sense that this method does NOT assume, eg,
        RR=0, since otherwise this method would be identical to notional() */
    virtual double maxPossibleLoss() const;

    /** Returns the loss config notional, ie, the amount that could be
        (or could have been) lost in the worst possible scenario, even if
        this scenario is no longer possible (e.g., if a name has defaulted
        with a RR > 0 it will still be accounted as if RR=0 if the losses
        are greater in this case).
        Note this notional bears no relation to the trade position */
    virtual double notional() const;

    /** Computes the loss and recovered notional currently produced by this
        loss config. Note both can be 0 no names have defaulted */
    virtual void currentLossAndRecovered(
        double& loss,               // (O)
        double& recoveredNotional,  // (O)
        const bool recoverNotional) const;

    /** Computes how much of this loss config's original notional has been
        prepaid (amortized) before "prepay date" */
    virtual double getPrepaidNotional(const DateTime& prepayDate) const;
    //
    // ICreditLossConfig methods
    //----------------------------------------


    //++++++++++++++++++++++++++++++++++++++++
    // Methods required for convolution purposes (note non-virtual)
    //
    /** Stores the strikes of the tranche in the method arguments (references) */
    void getTrancheStrikes(double& k1, double& k2) const;

    /** Stores the strikes of the tranche in the method arguments (references) */
    void getTrancheXSubAdjustedStrikes(double& k1, double& k2) const;

    /** Utility method: calculates survival probabilities along timeline.
        Resizes survivalProb as required. */
    void computeNameSurvivalProb(
        const DateTimeArray& timeline,            // (I)
        DoubleArrayArray&    survivalProb) const; // (O) [time index][name index]


    // JLHP These methods below should NOT be here, they just do not fit
    // into the new infrastrucutre...

    /** Returns the notional for the name corresponding to the
        specified index which must lie in range [0, numInnerLossConfigs()-1] */
    double nameNotional(int index) const;

    /** Returns the 'beta' for the name corresponding to the specified
        index which must lie in range [0, numNames()-1] */
    double nameBeta(int index) const;

    /** Returns the 'decretion beta' for the name corresponding to the specified
        index which must lie in range [0, numNames()-1] */
    double nameDecBeta(int index) const;

    /** Returns the recovery for the name corresponding to the
        specified index which must lie in range [0, numNames()-1]. Note that
        the recovery can  be overridden at the trade level hence this method,
        which reflects any trade level overrides, should be used rather than
        the recovery off the par spread curve */
    double nameRecovery(int index) const;

    /** Returns true if the name has defaulted where the name corresponds to the
        specified index which must lie in range [0, numNames()-1] */
    bool nameDefaulted(int index) const;

    /** Returns the ranges for Loss Given Default for the name
        corresponding to the specified index which must lie in range
        [0, numNames()-1]. In the particular, the LGD of a name is given by
        nameNotional * MIN(lgdCap, MAX(lgdFloor, lgdNotional - recovery)) */
    void nameLGDRanges(int     index,       /* (I) */
                       double& lgdNotional, /* (O) */
                       double& lgdFloor,    /* (O) */
                       double& lgdCap)      /* (O) */ const;

    /** Returns the max of valueDate and protection start date
        for specified name (if any). The specified index must lie in
        range [0, numNames()-1]*/
    const DateTime& nameProtectionStartDate(int index) const;

    /** Returns the min of last date parameter and protection end date
        for specified name (if any). The specified index must lie in
        range [0, numNames()-1]. This if for name 'cutoff' */
    const DateTime& nameProtectionEndDate(int index,
                                          const DateTime& lastDate) const;

    /** Returns the 'loss decretion beta' */
    double lossDecretionBeta() const;

    /** Returns the representation of the underlying name corresponding to the
        specified index which must lie in range [0, numNames()-1] */
    SingleCreditAssetConstSP nameAsset(int index) const;

    /** Returns the total notional of the portfolio */
    double portfolioNotional() const;

    /** Returns the total long notional of the portfolio */
    double portfolioLongNotional() const;

    /** Returns the total short notional of the portfolio */
    double portfolioShortNotional() const;

    /** Returns the TOTAL historical losses in the underlying portfolio,
        ie, before applying the tranche */
    double portfolioLoss() const;

    /** Returns the sum of the past long recovered notional of the
     * portfolio that ate into the total long notional returned by
     * portfolioLongNotional() */
    double portfolioLongRecoveredNotional() const;

    // jlhp Can (should) this be private?
    /** Returns highStrike-lowStrike with highStrike and lowStrike obtained from
        trancheStrikes */
    double trancheNotional() const;

    /** the outstanding notional as seen from today */
    double trancheOutstandingNotional(const bool recoverNotional) const;

    DateTime getToday() const;

    const double trancheLoss() const;

    const double trancheRecoveredNotional() const;
    //
    // Methods required for convolution purposes (note non-virtual)
    //----------------------------------------


    /** EffectiveLossCurveModelConfig::IIntoLossGen method:
        Create an IEffectiveLossCurveGen as specified by the supplied
        EffectiveLossCurveModelConfig */
    virtual IEffectiveCurveLossGenSP lossGenerator(
        IEffectiveCurveLossModelConfigConstSP effCurveLossModelConfig) const;


    //// For CDOParalleStrikeShift
    TweakOutcome sensShift(const PropertyTweak<CDOParallelStrike>& shift);
    string sensName(const CDOParallelStrike* shift) const;
    void sensRestore(const PropertyTweak<CDOParallelStrike>& shift);

    /** Theta::IShift method */
    virtual bool sensShift(Theta* shift);

//# Methods of the interface, ICreditTrancheLossConfigSVGen follow
    /** Creates the corresponding state variable  */
	virtual ICreditLossConfigSVGenConstSP createSVGen(
		const DateTimeLiteVectorConstSP& timeline,
		const CIntConstSP& triggerDelay,
		const CIntConstSP& defaultToCalculationDelay,
		double temporaryLossAmount,
		const DateTime& lastTriggerDate,
		const AccrualPeriodArrayConstSP& accrualPeriods,
		const IBadDayAdjusterConstSP& bda,
		const IProtectionProviderConstSP& protect,
		const IRebateCalculatorConstSP& rebateCalc,
        const bool recoverNotional
		) const;

//# Methods of the interface, ICreditTrancheIndexedLossConfigSVGen follow
    /** Creates the corresponding indexed state variable  */
	virtual ICreditLossConfigIndexedSVGenConstSP createIndexedSVGen(
		const DateTimeArrayConstSP& timeline,
		const CIntConstSP& triggerDelay,
		const CIntConstSP& defaultToCalculationDelay,
		double temporaryLossAmount,
		const DateTime& lastTriggerDate,
		const AccrualPeriodArrayConstSP& accrualPeriods,
		const IBadDayAdjusterConstSP& bda,
		const IProtectionProviderConstSP& protect,
		const IRebateCalculatorConstSP& rebateCalc,
        const bool recoverNotional
		) const;

    /** [Implements ITranchesCombinationPayoff] */
    virtual void linearDecomposition(
        const DateTime& time,
        DoubleArray& baseStrikes,        /* output */
        DoubleArray& baseStrikesWeights, /* output */
        double& expectedLossOffset       /* output */) const;

	virtual bool isCrossSub() const {return isXSub;};

	bool hasCutOffFlag() const {return hasCutOff;};

	const DateTime& getCutOffDate() const {return cutOffDate;};

	const DateTime& getProtectionStartDate() const;

private:

/** Forward declaring inner classes representing state variables generators representing loss event times
	*/
	FORWARD_DECLARE(SVGen);
	FORWARD_DECLARE(IndexedSVGen);

	/** Invoked once at start up when this class is 'loaded' */
    static void load(CClassSP& clazz);
    static IObject* defaultConstructor();
    CreditTrancheLossConfig();

    CreditTrancheLossConfig(const CreditTrancheLossConfig&);
    CreditTrancheLossConfig& operator=(const CreditTrancheLossConfig&);

    /** Sorts and aggregates the incoming "portfolioLosses" array, and
        returns a array of CtgLegLossPerDefault after applying the tranche 
        as defined by the k1 and k2 (low and high) strikes.
        CAUTION: the entries in the returned array are non-cumulative */
    static CtgLegLossPerDefaultArraySP getTrancheReductions(
        CtgLegLossPerDefaultArraySP portfolioLosses,
        const double                k1,
        const double                k2);

    /** Sorts and aggregates the incoming "portfolioLosses" array, and
        returns a array of FeeLegReductionPerDefault after applying the tranche 
        as defined by the k1 and k2 (low and high) strikes.
        CAUTION: the entries in the returned array are non-cumulative */
    static FeeLegReductionPerDefaultArraySP getTrancheReductions(
        FeeLegReductionPerDefaultArraySP portfolioLosses,
        const double                     k1,
        const double                     k2,
        const bool                       isLoss);

    // Fields
    string   name;           // Name
    DateTime valueDate;      // Today
    double   lowStrike;      // Tranche low strike (absolute value, not a percentage)
    double   highStrike;     // Tranche high strike (absolute value, not a percentage)
    ICreditLossConfigSP portfolio; // Portfolio of names in the tranche
    RationalisedCreditEngineParametersWrapper engineParams; // Parameters for this loss config
    CDOPortfolioConstSP portf; // unregistered JLHP - temporarily, keep a (cast) of the CDOPortfolio

	mutable bool isXSub;
	mutable bool isOuter;
	double sumBufferForXSub;

	bool hasCutOff;
	DateTime cutOffDate;

	DateTimeSP protectionStartDate;
};

DECLARE (CreditTrancheLossConfig);
typedef MarketWrapper<CreditTrancheLossConfig> CreditTrancheLossConfigWrapper;


DRLIB_END_NAMESPACE

#endif
