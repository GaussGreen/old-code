//----------------------------------------------------------------------------
//
//   Group       : CH Quantitative Research
//
//   Description : A portfolio of credit assets for use in a CDO
//
//   Author      : Antoine Gregoire
//
//----------------------------------------------------------------------------

#ifndef QLIB_CDOPORTFOLIO_HPP
#define QLIB_CDOPORTFOLIO_HPP

#include "edginc/TweakableWith.hpp"
#include "edginc/Atomic.hpp"
#include "edginc/CreditAsset.hpp"
#include "edginc/Theta.hpp"
#include "edginc/DateTime.hpp"
#include "edginc/PortfolioName.hpp"
#include "edginc/CreditMultiBetaLevel.hpp"
#include "edginc/CreditLossConfigMC.hpp"
#include "edginc/RationalisedCreditEngineParameters.hpp"
#include "edginc/FORWARD_DECLARE_WRAPPER.hpp"

DRLIB_BEGIN_NAMESPACE

/* Forward declare the classes refered to from this file if we are not
 * interested in their storage properties (ie, they are used through (smart)
 * pointers and therefore their include files are not required here - they
 * will be required in the .cpp file though)*/
FORWARD_DECLARE(SingleCreditAsset);
FORWARD_DECLARE(IProtectionProvider);
FORWARD_DECLARE_REF_COUNT(CtgLegLossPerDefault);
FORWARD_DECLARE_REF_COUNT(FeeLegReductionPerDefault);

class MARKET_DLL CDOPortfolio:
    public CObject,
    virtual public ICreditLossConfig,
    virtual public CreditMultiBetaLevel::IShift,
    virtual public Theta::Shift,
    virtual public ICreditLossConfigSVGenMC,  // so that it supports MC
    virtual public ICreditLossConfigIndexedSVGenMC  // to support indexed SVGen creation
{
public:
    /** TYPE (for reflection) */
    static CClassConstSP const TYPE;
    /** Destructor */
    virtual ~CDOPortfolio();

    //++++++++++++++++++++++++++++++++++++++++
    //ICreditLossConfig methods
    //
    /** Name of the single curve. */
    virtual string getName() const;

    /** Returns valuedate */
    virtual DateTime getToday() const;

    virtual void setEngineParameters(CreditEngineParametersSP engineParams);

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
        CIntConstSP               triggerDelay,
        CIntConstSP               defaultToCalculationDelay,
        const double              temporaryLossAmount,
        const DateTime&           lastTriggerDate,
        AccrualPeriodArrayConstSP accrualPeriods,
        IBadDayAdjusterConstSP    bda,
        const bool                recoverNotional) const;

    /** Compute the rebate payments caused by historic fee leg reductions */
    virtual CashFlowArraySP historicRebatePayments(
        const IRebateCalculator* const   rebateCalc,
        FeeLegReductionPerDefaultArraySP reductions,
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

    /** Returns the engine parameters for this loss cfg (NB Model needs to
        have selected appropriate set during getMarket()). The parameter
        can specify what type of engine parameters are required. An
        exception is thrown if the relevant type is not found */
    CreditEngineParametersConstSP getEngineParams() const { return engineParams.getSP(); }

    /** Returns the number of "inner loss configs" contained in this LossConfig.
        In this case, return the number of elements in the "names" array. */
    virtual int numInnerLossConfigs() const;

    /** Returns the inner loss config number "index" */
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
        are greater in this case). */
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

    virtual void validatePop2Object();

    virtual IObject* clone() const;

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
     * portfolio which eat into the total long notional returned by
     * portfolioLongNotional() */
    double portfolioLongRecoveredNotional() const;

    void getCompositionDetails(
        // Outputs
        double& totalLongNotional,
        double& totalShortNotional,
        double& totalPastLongLoss,
        double& totalPastShortLoss,
        double& totalPastRecoveredNotional,
        double& totalPastLongRecoveredNotional,
        double& totalPastShortRecoveredNotional) const;

    /** Returns number of names in that portfolio */
    int nbName() const;

    /** Returns the PortfolioName for the specified index which must lie in
        in the range [0, nbName()-1] */
    PortfolioNameConstSP getName(int index) const;

    /** Returns the representation of the underlying name corresponding to the
        specified index which must lie in range [0, numNames()-1] */
    virtual SingleCreditAssetConstSP nameAsset(int index) const;

    /** CreditMultiBetaLevel::IShift implementation */
    virtual string sensName(CreditMultiBetaLevel* shift) const;

    /** CreditMultiBetaLevel::IShift implementation */
    virtual bool sensShift(CreditMultiBetaLevel* shift);

    /** CAUTION: Do not allow the calculation dates of defaults which
        should and should not be paid get interleaved, ie, if a default
        should not be paid do not allow later defaults which should
        be paid to have earlier calculation dates, and viceversa.
        This is in order to avoid the wrong default hitting the tranche
        losses and producing default payments when there shouldn't be,
        or not producing them when there should be.
        Note: this will typically only happen if at least one of the names
        has a credit event override. */
    static void validateCtgLegLosses(CtgLegLossPerDefaultArraySP losses,
                                     const IProtectionProvider* const protect);


    //# Methods of the interface, ICreditLossConfigSVGenMC follow

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

    //# Methods of the interface, ICreditLossConfigIndexedSVGenMC follow

	/** Creates the corresponding state variable  */
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

    // ------
    // FIELDS
    // ------
    /* correlation between loss and decretion */
	// this should be moved into engine parameters 
	// it is still here for backward compatibility
    double lossDecretionBeta;

    /** Theta::IShift method */
    virtual bool sensShift(Theta* shift);

private:
	//nested class declaration for the state variables corresponding to this lossConfig
	FORWARD_DECLARE(SVGen);
	FORWARD_DECLARE(IndexedSVGen);

	/** Constructor */
    CDOPortfolio();

    /** Default constructor (for reflection) */
    static IObject* defaultConstructor();

    /** Invoked once at start up when this class is 'loaded' */
    static void load(CClassSP& clazz);

    // ------
    // FIELDS
    // ------
    /* portfolio names */
    ICreditLossConfigArray names;

    /* Name of the CDOPortfolio (used, e.g., to retrieve per name per
       portfolio per model parameters) */
    string portfolioName;
    DateTime valueDate;
    RationalisedCreditEngineParametersWrapper engineParams; // Parameters for this loss config
    PortfolioNameArray portfolioNames; // JLHP - temporarily, keep a cast of the PortfolioNames

	CDoubleSP externalNotional;
	double notionalFactor;

public:

	CDOPortfolio(PortfolioNameArray lclNames, CClassConstSP clazz = TYPE);

};

DECLARE(CDOPortfolio);

DRLIB_END_NAMESPACE

#endif

