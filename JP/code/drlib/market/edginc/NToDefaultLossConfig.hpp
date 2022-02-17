//----------------------------------------------------------------------------
//
//   Group       : Credit Hybrids QR
//
//   Description : An ICreditLossConfig used for N-to-default products.
//                 It implements the ISingleDefaultCreditLossConfig interface 
//                 since an NtD loses all its notional an a single default date.
//
//   Date        : October 2006
//
//----------------------------------------------------------------------------

#ifndef QLIB_NTODEFAULTLOSSCONFIG_HPP
#define QLIB_NTODEFAULTLOSSCONFIG_HPP

#include "edginc/Atomic.hpp"
#include "edginc/Theta.hpp"
#include "edginc/CashFlow.hpp"
#include "edginc/MarketWrapper.hpp"
#include "edginc/ICreditLossConfig.hpp"
#include "edginc/CreditEngineParameters.hpp"
#include "edginc/ISingleDefaultCreditLossConfig.hpp"
#include "edginc/IEffectiveCurveLossModelConfig.hpp"
#include "edginc/FORWARD_DECLARE.hpp"


DRLIB_BEGIN_NAMESPACE

/* Forward declare the classes refered to from this file if we are not
 * interested in their storage properties (ie, they are used through (smart)
 * pointers and therefore their include files are not required here - they
 * will be required in the .cpp file though)*/
FORWARD_DECLARE(IProtectionProvider);


class MARKET_DLL NToDefaultLossConfig : 
    public CObject,
    virtual public IEffectiveCurveLossModelConfig::IIntoEffCurveGen,
    virtual public ISingleDefaultCreditLossConfig,
    virtual public Theta::Shift
{

public:
    /** TYPE (for reflection) */
    static CClassConstSP const TYPE;

    virtual ~NToDefaultLossConfig();

    virtual void validatePop2Object();

    /** Gathers the CtgLegLosses corresponding to the Nth default in
        chronological default order */
    static CtgLegLossPerDefaultArraySP getLossesOfNthDefault(
        CtgLegLossPerDefaultArraySP portfolioLosses, 
        int n);

    //++++++++++++++++++++++++++++++++++++++++
    //ICreditLossConfig methods
    //
    /** All loss config instances have a 'name'. */
    virtual string getName() const;

    /** Returns valuedate */
    DateTime getToday() const;

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
        const bool recoverNotional) const;       // (O)

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

    /** Returns the number of "inner loss configs" contained in this LossConfig.
        In this case return the number of names in the portfolio. */
    virtual int numInnerLossConfigs() const;

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
    // ISingleDefaultCreditLossConfig methods
    //
    /** Returns true if this loss config has defaulted */
    virtual bool hasDefaulted() const;

    /** Returns the default date, or an empty DateTime if the loss config
        has not defaulted */
    virtual DateTime getDefaultDate() const;

    /** Returns the timeline used in this loss config, ie, the timepoints 
        defining the regions at which it is safe to assume that the (flat 
        forwards) default rate is constant - If this concept does not make 
        sense for this loss config, a null SP should be returned.
        In this case, a union of the timelines of the underlying names in
        this NtD */
    virtual DateTimeArraySP getTimeLine() const;
    //
    // ISingleDefaultCreditLossConfig methods
    //----------------------------------------

    /** EffectiveCurveLossModelConfig::IIntoLossGen method:
        Create an IEffectiveLossCurveGen as specified by the supplied
        EffectiveLossCurveModelConfig */
    ICreditLossGenSP effCurveGenerator(
         IEffectiveCurveLossModelConfigConstSP effCurveLossModelConfig,
         CounterPartyCreditConstSP             cpty,
         const bool                            recoverNotional,
         IModelConfigMapperConstSP             mapper) const;

    /** Returns the default number, N, of this NtD */
    const double getDefaultNumber() const;

    /** Theta::IShift method */
    virtual bool sensShift(Theta* shift);


private:
	/** Invoked once at start up when this class is 'loaded' */
    static void load(CClassSP& clazz);
    static IObject* defaultConstructor();
    NToDefaultLossConfig();
    
    NToDefaultLossConfig(const NToDefaultLossConfig&); // Not implemented, dont use
    NToDefaultLossConfig& operator=(const NToDefaultLossConfig&); // Not implemented, dont use

    // Methods
    /** Returns the Nth defaultedName (ie, the one triggering this NtD) if
        available, or null otherwise */
    ISingleDefaultCreditLossConfigSP getNthDefaultedName() const;

    // Fields
    string   name;          // Name
    DateTime valueDate;     // Today
    int      defaultNumber; // Number (index) of default covered for protection
//     bool     payFeesUntilCalcDate;
    ISingleDefaultCreditLossConfigArray basket;  // portfolio names
    CreditEngineParametersWrapper engineParams; // Parameters for this loss config
};

DECLARE(NToDefaultLossConfig);
typedef MarketWrapper<NToDefaultLossConfig> NToDefaultLossConfigWrapper;

DRLIB_END_NAMESPACE

#endif
