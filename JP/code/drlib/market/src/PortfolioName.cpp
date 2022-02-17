//----------------------------------------------------------------------------
//
//   Group       : CH Quantitative Research
//
//   Description : Each of the underlying "names" in, say, a CDO
//
//   Author      : Antoine Gregoire
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/Maths.hpp"
#include "edginc/Format.hpp"
#include "edginc/Model.hpp"
#include "edginc/MarketData.hpp"
#include "edginc/IBadDayAdjuster.hpp"
#include "edginc/ITrancheCreditEventOverride.hpp"
#include "edginc/IRebateCalculator.hpp"
#include "edginc/FeeLegReductionPerDefault.hpp"
#include "edginc/CtgLegLossPerDefault.hpp"
#include "edginc/IProtectionProvider.hpp"
#include "edginc/PortfolioName.hpp"
#include "edginc/CDSPricer.hpp"
#include "edginc/SingleCreditAsset.hpp"
#include "edginc/CmOnlyParameters.hpp"

DRLIB_BEGIN_NAMESPACE

/* Flag used to get an error message when beta is tweaked outside [-1,1] */
#define BETA_TWEAK_DEBUG 0

static double const BETA_LOW_CAP = -0.9999; /** Minimum value for betas */
static double const BETA_UP_CAP  =  1.0;    /** Maximum value for betas */

/** Computes betaTweak = beta + shift * (1 - beta) */
double PortfolioName::betaTweak(double beta, double shift, bool capAndFloor) {
    double result = beta + shift * (1.0 - beta);

    if (capAndFloor) {
        // Caps extreme values
        if (result < BETA_LOW_CAP) {
            return BETA_LOW_CAP;
        }
        if (result > BETA_UP_CAP) {
            return BETA_UP_CAP;
        }
    }
    else {
#if BETA_TWEAK_DEBUG
        // Checks that tweaked beta is inside [-1,1]
        checkRange("betaTweak", result, -1.0, 1.0);
#endif
    }
    return result;
}


/** Gets beta utility method. Throws an exception if this name's model
    parameters are not of type CmOnlyParameters */
double PortfolioName::getBeta() const {
    CreditEngineParametersConstSP engineParams =
        asset->getEngineParams(CmOnlyParameters::TYPE);

    CmOnlyParametersConstSP cmParams(
        dynamic_cast<const CmOnlyParameters*>(engineParams.get()));

    return cmParams->getBeta();
}


/** basic validation */
void PortfolioName::validatePop2Object() {
    static const string method("PortfolioName::validatePop2Object");
    const string& name = asset.getName();

    /* For simplicity in the implementation of credit events handling,
     * only allow the "default" lgdNotional, lgdFloor and lgdCap of
     * 1, 0, 1
    if (lgdNotional < 0.0 || lgdNotional > 1.0){
        throw ModelException(method, "Portfolio: lgdNotional for "
                             + name + " outside [0,1]");
    }
    if (lgdFloor < 0.0 || lgdFloor > 1.0){
        throw ModelException(method, "Portfolio: lgdFloor for "
                             + name + " outside [0,1]");
    }
    if (lgdCap < 0.0 || lgdCap > 1.0) {
        throw ModelException(method, "Portfolio: lgdCap for "
                             + name + " outside [0,1]");
    }
    if (lgdFloor > lgdCap){
        throw ModelException(method, "lgdFloor is greater than lgdCap "
                             "for name "+ name);
    }
    */

    if (!Maths::equals(lgdNotional, 1.0)) {
        throw ModelException(method, "Portfolio: lgdNotional for "
                             + name + " must be 1.0");
    }
    if (!Maths::isZero(lgdFloor)) {
        throw ModelException(method, "Portfolio: lgdFloor for "
                             + name + " must be 0.0");
    }
    if (!Maths::equals(lgdCap, 1.0)) {
        throw ModelException(method, "Portfolio: lgdCap for "
                             + name + " must be 1.0.");
    }

    // check recovery
    if (defaultParamOverride) {
        if (nameRecovery < 0.0 || nameRecovery > 1.0) {
            throw ModelException(method,
                                 "Portfolio: Recovery outside [0,1] for "
                                 + name + " (recovery="
                                 + Format::toString(nameRecovery) + ")");
        }
    }
}

/** populate from market cache */
void PortfolioName::getMarket(const IModel* model, const MarketData* market){
    static const string method("PortfolioName::getMarket");
    try {
        market->GetReferenceDate(valueDate);

        // Fetch the data for the asset, including the engine parameters
        asset.getData(model, market);
        // Now, potentially, replace the old-style engine parameters with the
        // new-style ones, which may include the beta
        asset->convertToNewParamStyle(getName(), beta, decretionBeta);

        if (!nameDefDate.empty()) {
            if (nameDefDate > valueDate)
                throw ModelException("PortfolioName::getMarket",
                                     "Portfolio: DEFAULT DATE ("
                                     + nameDefDate.toString() + ") > TODAY ("
                                     + valueDate.toString() + ") for " +
                                     asset.getName());
        }
        // If there is a credit event override, let it get its market data
        if (!!creditEventOverride) {
            creditEventOverride->getMarket(model, market);
        }
    }
    catch (exception& e) {
        throw ModelException(e, method, "Problem fetching the data for name '" +
                             getName() + "'.");
    }
}

/** CreditNameNotionalLevel::Shift implementation */
bool PortfolioName::sensShift(CreditNameNotionalLevel* shift) {
    nameNotional = shift->getShiftSize();
    return true;
}

string PortfolioName::sensName(CreditNameNotionalLevel* shift) const
{
    return asset->getName();
}


/** Destructor */
PortfolioName::~PortfolioName() {}

/** Returns the notional */
double PortfolioName::getNameNotional() const {
    return nameNotional;
}

/** Indicates whether the name has defaulted before the name-specific
 * "protectionStartDate" - returns false if:
 * - the name has not defaulted, or
 * - there is no name-specific protectionStartDate, or
 * - the name has defaulted on or after the protectionStartDate */
bool PortfolioName::defaultedBeforeProtectionStarts() const {
    const DateTime& defaultDate = getNameDefaultDate();
    if (!defaultDate.empty() && !!protectionStartDate) {
        if (defaultDate < *protectionStartDate) {
            return true;
        }
    }
    return false;
}


/**
 * Returns the recovery:
 * - field "nameRecovery" if defaultParamOverride=TRUE
 * - market recovery otherwise
 *  */
double PortfolioName::getNameRecovery() const {
    // Special case for forward starting names
    // Want to use recovery = 100% if default date < protection start date
    if (defaultedBeforeProtectionStarts()) {
        return 1.0;
    }

    double recRate;
    if (defaultParamOverride) {
        recRate = nameRecovery;
    }
    else {
        recRate = asset->recoveryRate();
    }

    if (!!creditEventOverride) {
        // The override may includes a rec rate
        recRate = creditEventOverride->getOverallRecoveryRate(recRate);
    }
    return recRate;
}

/**
 * Returns the default date:
 * - field "nameDefDate" if specified
 * - market default date (if any) otherwise
 *  */
DateTime PortfolioName::getNameDefaultDate() const {
    if (nameDefDate.empty()) {
        return asset->defaultDate();
    } else {
        return nameDefDate;
    }
}

/** Returns the underlying credit asset */
CreditAssetConstSP PortfolioName::getCreditAsset() const {
    return asset.getSP();
}

//++++++++++++++++++++++++++++++++++++++++
// ISingleDefaultCreditLossConfig methods
//
/** Returns true if this name has defaulted up to and including nameMatCutoff */
bool PortfolioName::hasDefaulted() const{
    const DateTime& defaultDate = getNameDefaultDate();
    if (defaultDate.empty()){
        return false;
    }
    if (nameMatCutoff.empty()){
        return true;
    }
    return (defaultDate <= nameMatCutoff);
}


/** Returns the default date, or an empty DateTime if the loss config
    has not defaulted */
DateTime PortfolioName::getDefaultDate() const {
    return getNameDefaultDate();
}


/** Returns the timeline used in this loss config, ie, the timepoints
    defining the regions at which it is safe to assume that the (flat
    forwards) default rate is constant. */
DateTimeArraySP PortfolioName::getTimeLine() const {
    DateTimeArray dates = asset->getParSpreadCurve()->zeroDates();
    return DateTimeArraySP(new DateTimeArray(dates));
}
//
// ISingleDefaultCreditLossConfig methods
//----------------------------------------

/** Returns protection start date */
const DateTime& PortfolioName::getProtectionStartDate() const
{
    if (!protectionStartDate) {
        return valueDate;
    }
    else {
        return protectionStartDate->max(valueDate);
    }
}

/** Returns the min of last date parameter and protection end date
    for specified name (if any). This if for name 'cutoff' */
const DateTime& PortfolioName::getProtectionEndDate(
    const DateTime& lastDate) const
{
    if (nameMatCutoff.empty()){
        return lastDate;
    }
    return nameMatCutoff.min(lastDate);
}

void
PortfolioName::AdjustNotionalAndDates(
	double newNotional,
	DateTime newStart,
	DateTime newCutOff)
{
	nameNotional  = newNotional;
	nameMatCutoff = newCutOff;
	protectionStartDate = DateTimeSP(new DateTime(newStart));
};


/** Returns the name maturity cut off for the name **/
const DateTime& PortfolioName::getNameMaturityCutOff() const
{
	return nameMatCutoff;
}

//// computes MIN(k2, MAX(k1, x)). Must have k2 >= k1
static double lgdCalc(double x, double k1, double k2) {
    if (x<=k1){
        return k1;
    }
    if (x>=k2){
        return k2;
    }
    return x;
}

/** Returns
    nameNotional * MIN(lgdCap, MAX(lgdFloor, lgdNotional - recovery)) */
double PortfolioName::lossGivenDefault() const{
    double notional = nameNotional;
    if (!!creditEventOverride) {
        notional *= creditEventOverride->getOverallDefaultedNotionalFraction();
    }

    return notional * lgdCalc(lgdNotional-getNameRecovery(), lgdFloor, lgdCap);
}


/** Indicates whether this name has a credit event override */
bool PortfolioName::hasEventOverride() const {
    return !!creditEventOverride;
}

/** Indicates whether this name has a default param override */
bool PortfolioName::defaultParamOverrideExists() const
{
	return defaultParamOverride;
}


//++++++++ ICreditLossConfig methods
/** Name of the underlying asset */
string PortfolioName::getName() const {
    return asset->getName();
}

DateTime PortfolioName::getToday() const {
    return valueDate;
}

void PortfolioName::impliedParSpreadsAndDurations(
    const YieldCurveConstSP discount,
    const DateTimeArray& dates,
    DoubleArray& impliedSpreads,  /* (Output) */
    DoubleArray& durations) const /* (Output) */
{
    CDSPricer::impliedParSpreadsAndDurations(
        DateTimeSP(new DateTime(valueDate)),
        DateTimeSP(new DateTime(getProtectionStartDate())),
        SingleCreditAssetConstSP::dynamicCast(getCreditAsset())->getParSpreadCurve(),
        discount,
        dates,
        impliedSpreads,
        durations);
}

/** Returns an array with the cashflows corresponding to all "real"
    losses due to the settlement of this name's default.
    If there are no losses, it may return an empty (0 size)
    CtgLegLossPerDefaultArraySP or a "null" SP.
    Note: the cashflow dates are the CALCULATION DATES (as opposed to
    settlement dates) when the corresponding losses are determined.
    CAUTION: Assumes that the defaultDate has been verified to fall within
    the protection period. */
CtgLegLossPerDefaultArraySP PortfolioName::historicContingentLegLosses(
    CIntConstSP triggerDelay,
    CIntConstSP defaultToCalculationDelay,
    const DateTime& lastTriggerDate,
    IBadDayAdjusterConstSP bda,
    const IProtectionProvider* const protect) const
{
    static const string method("PortfolioName::historicContingentLegLosses");

    // Should keep track of past deliveries to compute the effective
    // (remaining) lgdNotional, lgdFloor and lgdCap - and adjust each
    // cashflow amount doing something like:
    //   (nameNotional *
    //             lgdCalc(lgdNotional-getNameRecovery(), lgdFloor, lgdCap));
    // Not yet done - Instead verify that these parameters take the "default"
    // values (1, 0, 1) in validatePop2Object

    try {
        CtgLegLossPerDefaultArraySP losses;

        if (protect) {
            throw ModelException(method, "PortfolioName can not be used to "
                                 "produce a CDS-style instrument");
        }

        if (!hasDefaulted() || defaultedBeforeProtectionStarts()) {
            // This is a forward starting tranche and requires special
            // treatment: Assume RR = 100% and therefore there is no loss
            return CtgLegLossPerDefaultArraySP();
        }

        const DateTime& defaultDate = getNameDefaultDate();
        if (!creditEventOverride) {
            // There is no override, so use the values in the instrument

            // Note: Even if triggerDelay is 0, still invoke rollAndAdjustDates
            // in order to take the lastTriggerDate into account and BD-adjust
            // the default date
            const DateTime& triggerDate =
                ITrancheCreditEventOverride::rollAndAdjustDate(
                    defaultDate,
                    triggerDelay,
                    valueDate,
                    valueDate,
                    lastTriggerDate,
                    bda);

            if (triggerDate.empty()) {
                // The default has not been triggered, or it was triggered
                // after the lastTriggerDate -> no losses. Note losses is
                // already "null"
                ;
            }
            else {
                DateTime calculationDate =
                    ITrancheCreditEventOverride::rollAndAdjustDate(
                        defaultDate,
                        defaultToCalculationDelay,
                        valueDate,
                        triggerDate,
                        DateTime(), //lastTriggerDate not applicable for calculationDate
                        bda);

                losses.reset(new CtgLegLossPerDefaultArray(
                    1, CtgLegLossPerDefaultSP(new CtgLegLossPerDefault(
                        defaultDate,
                        calculationDate,
                        lossGivenDefault()))));
            }
        }
        else {
            // There is an override, so use it and ignore triggerDelay and
            // defaultToCalculationDelay
            losses = creditEventOverride->historicContingentLegLosses(
                nameNotional,
                getNameRecovery(),
                lastTriggerDate,
                defaultDate,
                bda);
        }
        return losses;
    }
    catch (exception& e) {
        throw ModelException(e, method, "Error for name " + getName());
    }
}


/** Computes the notional reductions on the fee leg (due to losses and/or
    recovered notional) and the corresponding fee rebates, if any - based
    on the creditEventOverride in this name.
    If there are no losses, it may return an empty (0 size)
    CtgLegLossPerDefaultArraySP or a "null" SP.
    The accrualPeriods are an argument because they are required to
    determine the actual date of the notional reductions (in the fee leg,
    all notional reductions happen either on the determination date or on
    an accrual start date.
    NOTE: The recoverNotional flag is ignored */
FeeLegReductionPerDefaultArraySP PortfolioName::historicFeeLegReductions(
    CIntConstSP triggerDelay,
    CIntConstSP defaultToCalculationDelay,
    const double temporaryLossAmount,
    const DateTime& lastTriggerDate,
    AccrualPeriodArrayConstSP accrualPeriods,
    IBadDayAdjusterConstSP bda,
    const bool recoverNotional) const
{
    static const string method("PortfolioName::historicFeeLegReductions");

    // Should keep track of past deliveries to compute the effective
    // (remaining) lgdNotional, lgdFloor and lgdCap - and adjust each
    // cashflow amount doing something like:
    //   (nameNotional *
    //             lgdCalc(lgdNotional-getNameRecovery(), lgdFloor, lgdCap));
    // not yet done
    try {

        FeeLegReductionPerDefaultArraySP feeReductions;

        const DateTime& defaultDate = getNameDefaultDate();

        if (!hasDefaulted()) {
            // No losses or recovered amounts
        }
        else if (defaultedBeforeProtectionStarts()) {
            // This is a forward starting instrument and requires special
            // treatment: Assume RR = 100% and therefore there will be no
            // losses, all the notional will be recovered
            feeReductions.reset(new FeeLegReductionPerDefaultArray(
                1, FeeLegReductionPerDefaultSP(new FeeLegReductionPerDefault(
                    defaultDate,       // determination date
                    defaultDate,       // effective date = loss date -> no rebate
                    defaultDate,       // calculation date
                    getNameNotional(), // defaulted notional
                    getNameNotional(), // total notional
                    1.0))));           // (assumed) recovery rate
        }
        else if (!creditEventOverride) {
            // Name defaulted after protection starts and there is no override,
            // so use the values in the instrument.
            // Note: Even if triggerDelay is 0, still invoke rollAndAdjustDate
            // in order to take the lastTriggerDate into account and BD-adjust
            // the default date
            const DateTime& triggerDate =
                ITrancheCreditEventOverride::rollAndAdjustDate(
                    defaultDate,
                    triggerDelay,
                    valueDate,
                    valueDate,
                    lastTriggerDate,
                    bda);

            if (triggerDate.empty()) {
                // The default has not been triggered, or it was triggered after
                // the lastTriggerDate -> no losses or recovered amounts.
                // Leave feeReductions as a NULL SP
            }
            else {
                DateTime calculationDate;
                if (!defaultToCalculationDelay) {
                    calculationDate = triggerDate;
                }
                else {
                    calculationDate =
                        ITrancheCreditEventOverride::rollAndAdjustDate(
                            defaultDate,
                            defaultToCalculationDelay,
                            valueDate,
                            triggerDate,
                            bda);
                }
                feeReductions.reset(new FeeLegReductionPerDefaultArray(
                    1, FeeLegReductionPerDefaultSP(new FeeLegReductionPerDefault(
                        triggerDate,           // determination date
                        calculationDate,       // effective date
                        calculationDate,
                        getNameNotional(),     // defaulted notional
                        getNameNotional(),     // total notional
                        getNameRecovery())))); // recovery rate
            }
        }
        else {
            // Name defaulted after protection starts and there is an override
            feeReductions = creditEventOverride->historicFeeLegReductions(
                nameNotional,
                getNameRecovery(),
                lastTriggerDate,
                defaultDate,
                bda);
        }
        return feeReductions;
    }
    catch (exception& e) {
        throw ModelException(e, method, "Error for name " + asset.getName());
    }
}


/** Compute the rebate payments caused by historic fee leg reductions */
CashFlowArraySP PortfolioName::historicRebatePayments(
    const IRebateCalculator* const   rebateCalc,
    FeeLegReductionPerDefaultArraySP reductions,
    IForwardRatePricerSP             model,
    const bool                       recoverNotional) const
{
    if (rebateCalc) {
        throw ModelException("JLHP PortfolioName::historicRebatePayments",
                             "Internal error: cannot compute rebates here yet");
    }
    return CashFlowArraySP();
}

/** Returns the engine parameters for this name (NB Model needs to
    have selected appropriate set during getMarket()). The parameter
    can specify what type of engine parameters are required. An
    exception is thrown if the relevant type is not found */
CreditEngineParametersConstSP PortfolioName::getEngineParams(
    CClassConstSP engineParamsType) const
{
    return asset->getEngineParams(engineParamsType);
}

/** Returns the number of "inner loss configs" contained in this LossConfig.
    In this case returns '1'. */
int PortfolioName::numInnerLossConfigs() const {
    return 1;
}

/** Returns the inner loss config number "index".
    In this case throws an exception. */
ICreditLossConfigConstSP PortfolioName::getInnerLossConfig(
    const int index) const
{
    throw ModelException("PortfolioName::getInnerLossConfig",
                         "Internal error: no inner loss configs in "
                         "a PortfolioName");
}

/** Returns lossGivenDefault() */
double PortfolioName::maxPossibleLoss() const {
    return lossGivenDefault();
}

/** Returns nameNotional */
double PortfolioName::notional() const {
    return nameNotional;
}

/** Computes the loss and recovered notional currently produced by this
    name. Note both can be 0 if this name has not defaulted */
void PortfolioName::currentLossAndRecovered(
    double& loss,               // (O)
    double& recoveredNotional,  // (O)
    const bool recoverNotional) const
{
    // Do not bother checking the recoverNotional flag, since it costs
    // close to nothing to compute the recoveredNotional
    if (hasDefaulted()) {
        if (defaultedBeforeProtectionStarts()) {
            // This is a forward starting tranche and requires special
            // treatment: assume RR = 100% and therefore there is no
            // loss - all the notional is, potentially, recovered.
            loss = 0.0;
            recoveredNotional += nameNotional;
        }
        else {
            loss = lossGivenDefault();
            recoveredNotional = nameNotional - loss;
        }
    }
    else {
        loss = recoveredNotional = 0.0;
    }
    return;
}

/** Computes how much of this loss config's original notional has been
    prepaid (amortized) before "prepay date" */
double PortfolioName::getPrepaidNotional(const DateTime& prepayDate) const {
    double prepaidNotional(0.0);
    if (hasDefaulted()) {
        // JLH. Should we not return what the prepaid notional is, regardless
        // of the default status? we did not use to and this behaviour has
        // been left unchanged
    }
    else {
        IDecretionCurveConstSP decretion =
            asset->getParSpreadCurve()->getPrepayCurve();

        if (!!decretion) {
            double currentFactor = decretion->getFactor(valueDate);
            prepaidNotional = nameNotional * (1 - currentFactor);
        }
    }
    return prepaidNotional;
}
//-------- ICreditLossConfig methods


//-------- ICreditLossConfigSVGenMC methods
ICreditLossConfigSVGenConstSP
PortfolioName::createSVGen(
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
	) const
{
    /** Note recoverNotional is not yet supported in
        PortfolioName::SVGen */
	return ICreditLossConfigSVGenConstSP(
		new PortfolioName::SVGen(
			PortfolioNameConstSP(this),
			timeline,
			triggerDelay,
			defaultToCalculationDelay,
			temporaryLossAmount,
			lastTriggerDate,
			accrualPeriods,
			bda,
			protect,
			rebateCalc));
}

//-------- ICreditLossConfigSVGenMC methods

//-------- ICreditLossConfigIndexedSVGenMC methods
ICreditLossConfigIndexedSVGenConstSP
PortfolioName::createIndexedSVGen(
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
	) const
{
    /** Note recoverNotional is not yet supported in
        PortfolioName::IndexedSVGen */
	return ICreditLossConfigIndexedSVGenConstSP(
		new PortfolioName::IndexedSVGen(
			timeline,
			PortfolioNameConstSP(this),
			triggerDelay,
			defaultToCalculationDelay,
			temporaryLossAmount,
			lastTriggerDate,
			accrualPeriods,
			bda,
			protect,
			rebateCalc));
}

//-------- ICreditLossConfigIndexedSVGenMC methods


/** Theta::IShift method */
bool PortfolioName::sensShift(Theta* shift) {
    // alter immediate data
    valueDate = shift->rollDate(valueDate);
    return true; // then shift components
}


/** Returns the ranges for Loss Given Default for this name.
    In the particular, the LGD of a name is given by
    nameNotional * MIN(lgdCap, MAX(lgdFloor, lgdNotional - recovery)) */
void PortfolioName::nameLGDRanges(double& lgdNotional, /* (O) */
                                  double& lgdFloor,    /* (O) */
                                  double& lgdCap)      /* (O) */ const{
    lgdNotional = this->lgdNotional;
    lgdCap = this->lgdCap;
    lgdFloor = this->lgdFloor;
}


/** Constructor for PortfolioName */
PortfolioName::PortfolioName():
    CObject(TYPE),
    lgdNotional(1.0),
    lgdFloor(0.0),
    lgdCap(1.0),
    nameRecovery(0.0)
{}

/** Default constructor (for reflection) */
IObject* PortfolioName::defaultConstructor(){
    return new PortfolioName();
}

/** Invoked once at start up when this class is 'loaded' */
void PortfolioName::load(CClassSP& clazz){
    clazz->setPublic();
    REGISTER(PortfolioName, clazz);
    SUPERCLASS(CObject);
    IMPLEMENTS(CreditNameNotionalLevel::Shift);
    IMPLEMENTS(ISingleDefaultCreditLossConfig);
    IMPLEMENTS(Theta::Shift);

    EMPTY_SHELL_METHOD(defaultConstructor);
    FIELD(nameNotional, "Notional of the name");
    FIELD(lgdNotional,
        "LGD = Loss Given Default "
        "= NameNotional * max(lgdFloor,min(lgdCap, lgdNotional - recovery)). "
        "Default values are : lgdNotional = 1, lgdFloor = 0, lgdCap = 1 ");
    FIELD_MAKE_OPTIONAL(lgdNotional);
    FIELD(lgdFloor,
        "LGD = Loss Given Default "
        "= NameNotional * max(lgdFloor,min(lgdCap, lgdNotional - recovery)). "
        "Default values are : lgdNotional = 1, lgdFloor = 0, lgdCap = 1 ");
    FIELD_MAKE_OPTIONAL(lgdFloor);
    FIELD(lgdCap,
        "LGD = Loss Given Default "
        "= NameNotional * max(lgdFloor,min(lgdCap, lgdNotional - recovery)). "
        "Default values are : lgdNotional = 1, lgdFloor = 0, lgdCap = 1 ");
    FIELD_MAKE_OPTIONAL(lgdCap);
    FIELD(asset, "credit asset");
    FIELD(beta, "beta");
    FIELD_MAKE_OPTIONAL(beta); // Not present here if using new style model parameters
    FIELD(decretionBeta, "notional decretion correlation"
                 "Default value is 0 (no correlation)");
    FIELD_MAKE_OPTIONAL(decretionBeta);
    FIELD(nameMatCutoff, "Maturity cut off date of the name");
    FIELD_MAKE_OPTIONAL(nameMatCutoff);
    FIELD(defaultParamOverride,
                 "TRUE: get recovery and default date from portfolio "
                 "FALSE: get recovery and default date from market");
    FIELD(nameRecovery, "Recovery in % override");
    FIELD_MAKE_OPTIONAL(nameRecovery);
    FIELD(nameDefDate, "Default date of the name");
    FIELD_MAKE_OPTIONAL(nameDefDate);
    FIELD(protectionStartDate,
        "Protection start date "
        "(used to book tranches on forward starting CDS)");
    FIELD_MAKE_OPTIONAL(protectionStartDate);
    FIELD(creditEventOverride,
          "Information about the credit event settlement for this name");
    FIELD_MAKE_OPTIONAL(creditEventOverride);
    FIELD       (valueDate, "Value date");
    FIELD_MAKE_OPTIONAL(valueDate);
}


CClassConstSP const PortfolioName::TYPE = CClass::registerClassLoadMethod(
    "PortfolioName", typeid(PortfolioName), load);

DEFINE_TEMPLATE_TYPE(PortfolioNameArray);

void PortfolioName::NameView::load(CClassSP& clazz)
{
    clazz->setPublic();
    REGISTER(NameView, clazz);
    SUPERCLASS(CObject);

	FIELD(name, "Name");

	FIELD(notional, "Notional");

	FIELD(protectionStart, "Protection start");

	FIELD(protectionEnd, "Risk cut-off");

	FIELD(recovery, "Recovery rate");

	EMPTY_SHELL_METHOD(defaultNameView);
}

CClassConstSP const PortfolioName::NameView::TYPE = CClass::registerClassLoadMethod(
	"PortfolioName::NameView", typeid(PortfolioName::NameView),  PortfolioName::NameView::load);


typedef  PortfolioName::NameViewArray PortfolioNameNameViewArray;

DEFINE_TEMPLATE_TYPE_WITH_NAME("PortfolioName::NameViewArray",
							   PortfolioNameNameViewArray);

// ################################################################################################################
/** Implementation for classes, PortfolioName::SVGen and SV */
PortfolioName::SVGen::SVGen(
	const PortfolioNameConstSP&  portfolioName,
	const DateTimeLiteVectorConstSP& timeline,
	const CIntConstSP& triggerDelay,
	const CIntConstSP& defaultToCalculationDelay,
	double temporaryLossAmount,
	const DateTime& lastTriggerDate,
	const AccrualPeriodArrayConstSP& accrualPeriods,
	const IBadDayAdjusterConstSP& bda,
	const IProtectionProviderConstSP& protect,
	const IRebateCalculatorConstSP& rebateCalc)
	:timeline(timeline),
	name(portfolioName->getName()),
	portfolioName(portfolioName),
	triggerDelay(triggerDelay),
	defaultToCalculationDelay(defaultToCalculationDelay),
	temporaryLossAmount(temporaryLossAmount),
	lastTriggerDate(lastTriggerDate),
	accrualPeriods(accrualPeriods),
	bda(bda),
	protect(protect),
	rebateCalc(rebateCalc)
{}

const string&
PortfolioName::SVGen::getName() const
{
	return name;
}

ICreditLossConfigSVSP
PortfolioName::SVGen::createNewSV(IStateVariableGen::IStateGen* stateGen) const
{
	try
	{
		//this method gets called twice - think of how to remove it later
		CtgLegLossPerDefaultArraySP histContLosses;

/*		These calls have been suppressed as these methods,  historicContingentLegLosses(...) and historicFeeLegReductions(...)
		are not working as of now.
		Assuming that later, the handling would include whether the default of the asset would be between start and end date!

		CtgLegLossPerDefaultArraySP histContLosses = portfolioName->historicContingentLegLosses(triggerDelay,
																								defaultToCalculationDelay,
																								lastTriggerDate,
																								bda.get(),
																								protect.get());

*/
		FeeLegReductionPerDefaultArraySP feeLosses;
		FeeLegReductionPerDefaultArraySP feeRecovered;

/*		CashFlowArraySP outputRebates;
		portfolioName->historicFeeLegReductions(feeLosses, feeRecovered, portfolioName->getNameNotional(), triggerDelay, defaultToCalculationDelay,
			temporaryLossAmount, lastTriggerDate, accrualPeriods, bda.get(), rebateCalc.get(), outputRebates);
*/
		CDoubleSP recoveryOverride;
		if (portfolioName->defaultParamOverrideExists())
			recoveryOverride = CDoubleSP(CDouble::create(portfolioName->getNameRecovery()));

		return ICreditLossConfigSVSP(
			new SV(
				this->getName(),
				portfolioName->getNameNotional(),
				portfolioName->getCreditAsset(),
				portfolioName->getBeta(),
				recoveryOverride,
				portfolioName->getProtectionStartDate(),
				portfolioName->getProtectionEndDate(timeline->back()),
				histContLosses,
				feeLosses,
				feeRecovered,
				timeline));
	}
	catch(exception& e)
	{
		throw ModelException(e,"PortfolioName::SVGen::createNewSV");
	}
}


ICreditLossConfigSVSP
PortfolioName::SVGen::getLossSV(
	IStateVariableSP oldStateVar,
	IStateVariableGen::IStateGen* stateGen) const
{
	try
	{
		IStateVariableSP stateVariable = stateGen->create(this);
		SV* portfolioNameSV = dynamic_cast<SV*>(stateVariable.get());
		return ICreditLossConfigSVSP(portfolioNameSV);
	}
	catch(exception& e)
	{
		throw ModelException(e,"PortfolioName::SVGen::getLossSV");
	}
}

IStateVariableSP
PortfolioName::SVGen::create(
	IStateVariableSP oldStateVar,
	IStateVariableGen::IStateGen* stateGen) const
{
	try
	{
		return getLossSV(oldStateVar, stateGen);
	}
	catch(exception& e)
	{
		throw ModelException(e,"PortfolioName::SVGen::create");
	}
}

void
PortfolioName::SVGen::collectStateVars(IStateVariableCollectorSP svdb) const
{
	svdb->appendElementary(this);
}


void
PortfolioName::SVGen::attachSVGen(class IElemStateVariableGenVisitor* sv) const
{
	//sv->processSVGen(this);
}


PortfolioName::SVGen::SV::SV(
	const string& name,
	double notional,
	const CreditAssetConstSP& asset,
	double beta,
	const CDoubleConstSP& recoveryOverride,
	const DateTime& protectionStartDate,
	const DateTime& maturityCutOff,
	const CtgLegLossPerDefaultArrayConstSP& histContLosses,
	const FeeLegReductionPerDefaultArrayConstSP& feeLosses,
	const FeeLegReductionPerDefaultArrayConstSP& feeRecovered,
	const DateTimeLiteVectorConstSP& timeline)
	:
	timeline(timeline),
	name(name),
	notional(notional),
	asset(asset),
	beta(beta),
	histContLosses(histContLosses),
	feeLosses(feeLosses),
	feeRecovered(feeRecovered),
	recoveryOverride(recoveryOverride),
	protectionStartDate(protectionStartDate),
	maturityCutOff(maturityCutOff)
{}


const CreditLossConfigSVResultList&
PortfolioName::SVGen::SV::getResults() const
{
	this->calculateLossEvents();
	return creditLossConfigSVResults;
}

const string&
PortfolioName::SVGen::SV::getName() const
{
	return name;
}

CreditLossConfigSVResultList&
PortfolioName::SVGen::SV::results()
{
	this->calculateLossEvents();
	return creditLossConfigSVResults;
}

const PortfolioNameSVSubResultVectorConstSP&
PortfolioName::SVGen::SV::getSubResults()
{
	return subResults;
}

void
PortfolioName::SVGen::SV::setSubResults(const PortfolioNameSVSubResultVectorConstSP&  rhsSubResults) const
{
	subResults = rhsSubResults;
}

void
PortfolioName::SVGen::SV::calculateLossEvents() const
{
/*
	CAUTION: Past defaults not handled. Implement something on the lines of Indexed SV
*/
	creditLossConfigSVResults.clear();

	for (unsigned int i=0; i< subResults->size(); ++i)
	{
		const DateTimeLite& date = (*subResults)[i].defaultTime;

		if ( (date < protectionStartDate) || (maturityCutOff < date) )   //asset has not defaulted
			continue;

		double recovery = (*subResults)[i].recovery;

		if (recoveryOverride.get())
			recovery = recoveryOverride->doubleValue();

		creditLossConfigSVResults.push_back(CreditLossConfigSVResult(
																date,
																notional*(1-recovery),
																notional,
																0));
	}
}

bool
PortfolioName::SVGen::SV::doingPast() const
{
	return false;
}

double
PortfolioName::SVGen::SV::getNotional() const
{
	return notional;
}

const CreditAssetConstSP&
PortfolioName::SVGen::SV::getAsset() const
{
	return asset;
}

double
PortfolioName::SVGen::SV::getBeta() const
{
	return beta;
}

const DateTimeLiteVectorConstSP&
PortfolioName::SVGen::SV::getTimeline() const
{
	return timeline;
}

// ################################################################################################################
/** Implementation for classes, PortfolioName::IndexedSVGen*/

PortfolioName::IndexedSVGen::IndexedSVGen(
	const DateTimeArrayConstSP& timeline,
	const PortfolioNameConstSP&  portfolioName,
	const CIntConstSP& triggerDelay,
	const CIntConstSP& defaultToCalculationDelay,
	double temporaryLossAmount,
	const DateTime& lastTriggerDate,
	const AccrualPeriodArrayConstSP& accrualPeriods,
	const IBadDayAdjusterConstSP& bda,
	const IProtectionProviderConstSP& protect,
	const IRebateCalculatorConstSP& rebateCalc)
	:timeline(timeline),
	name(portfolioName->getName()),
	portfolioName(portfolioName),
	triggerDelay(triggerDelay),
	defaultToCalculationDelay(defaultToCalculationDelay),
	temporaryLossAmount(temporaryLossAmount),
	lastTriggerDate(lastTriggerDate),
	accrualPeriods(accrualPeriods),
	bda(bda),
	protect(protect),
	rebateCalc(rebateCalc)
{}

const string&
PortfolioName::IndexedSVGen::getName() const
{
	return name;
}

ICreditLossConfigIndexedSVSP
PortfolioName::IndexedSVGen::createNewSV(IStateVariableGen::IStateGen* stateGen) const
{
	try
	{
		DateTimeConstSP defaultSettDate;

		CtgLegLossPerDefaultArraySP histContLosses =
			portfolioName->historicContingentLegLosses(
				triggerDelay,
				defaultToCalculationDelay,
				lastTriggerDate,
				bda,
				0);

		if (histContLosses.get()) //check with Jose, if we are expecting more than one default per name
			if (histContLosses->size() > 0)
				defaultSettDate = DateTimeConstSP(new DateTime((*histContLosses)[0]->defaultDate));

		FeeLegReductionPerDefaultArraySP feeLosses;
		FeeLegReductionPerDefaultArraySP feeRecovered;

/*		CashFlowArraySP outputRebates;
		portfolioName->historicFeeLegReductions(feeLosses, feeRecovered, portfolioName->getNameNotional(), triggerDelay, defaultToCalculationDelay,
			temporaryLossAmount, lastTriggerDate, accrualPeriods, bda.get(), rebateCalc.get(), outputRebates);
*/

/*		Past default handling - this will go once, the function calls, historicContingentLegLosses(...) and historicFeeLegReductions(...)
		start working
*/
		double lossGivenDefault =  portfolioName->lossGivenDefault();

		SingleCreditAssetConstSP singleCreditAsset = SingleCreditAssetConstSP::dynamicCast(portfolioName->getCreditAsset());
		ICDSParSpreadsConstSP cdsParSpreads = singleCreditAsset->getParSpreadCurve();

		int protectionStartIndex;
		try
		{
			const DateTime& protectionStartDate = portfolioName->getProtectionStartDate();
			protectionStartIndex = protectionStartDate.find(*timeline);
		}
		catch(exception& e)
		{
			throw ModelException(
				e,
				"PortfolioName::IndexedSVGen::createNewSV",
				"ProtectionStartDate not found on the product timeline");
		}

		int maturityCutOffIndex = -1;
		try
		{
			const DateTime& maturityCutOff = portfolioName->getNameMaturityCutOff();
			if (!maturityCutOff.empty())
				maturityCutOffIndex = maturityCutOff.find(*timeline);
		}
		catch(exception& e)
		{
			throw ModelException(
				e,
				"PortfolioName::IndexedSVGen::createNewSV",
				"MaturityCutOff not found on the product timeline");
		}

//calculate the recovery override
		CDoubleSP recoveryOverride;
		if (portfolioName->defaultParamOverrideExists())
			recoveryOverride = CDoubleSP(CDouble::create(portfolioName->getNameRecovery()));

//return the state variable
		return ICreditLossConfigIndexedSVSP(
			new SV(
				timeline,
				this->getName(),
				portfolioName->getNameNotional(),
				portfolioName->getCreditAsset(),
				portfolioName->getBeta(),
				portfolioName->hasDefaulted(),
				lossGivenDefault,
				defaultSettDate,
				recoveryOverride,
				protectionStartIndex,
				maturityCutOffIndex,
				feeLosses,
				feeRecovered));
	}
	catch(exception& e)
	{
		throw ModelException(e,"PortfolioName::IndexedSVGen::createNewSV");
	}
}


ICreditLossConfigIndexedSVSP
PortfolioName::IndexedSVGen::getLossSV(
	IStateVariableSP oldStateVar,
	IStateVariableGen::IStateGen* stateGen) const
{
	try
	{
		IStateVariableSP stateVariable = stateGen->create(this);
		SV* portfolioNameIndexedSV = dynamic_cast<SV*>(stateVariable.get());
		return ICreditLossConfigIndexedSVSP(portfolioNameIndexedSV);
	}
	catch(exception& e)
	{
		throw ModelException(e,"PortfolioName::IndexedSVGen::getLossSV");
	}
}

IStateVariableSP
PortfolioName::IndexedSVGen::create(
	IStateVariableSP oldStateVar,
	IStateVariableGen::IStateGen* stateGen) const
{
	try
	{
		return getLossSV(oldStateVar, stateGen);
	}
	catch(exception& e)
	{
		throw ModelException(e,"PortfolioName::IndexedSVGen::create");
	}
}

void
PortfolioName::IndexedSVGen::collectStateVars(IStateVariableCollectorSP svdb) const
{
	svdb->appendElementary(this);
}


void
PortfolioName::IndexedSVGen::attachSVGen(class IElemStateVariableGenVisitor* sv) const
{
	//sv->processSVGen(this);
}

// ################################################################################################################
/** Implementation for classes, PortfolioName::IndexedSVGen::SV*/

PortfolioName::IndexedSVGen::SV::SV(
	const DateTimeArrayConstSP& timeline,
	const string& name,
	double notional,
	const CreditAssetConstSP& asset,
	double beta,
	bool hasDefaulted,
	double lossGivenDefault,
	const DateTimeConstSP& defaultSettDate,
	const CDoubleConstSP& recoveryOverride,
	int protectionStartDateIndex,
	int maturityCutOffIndex,
	const FeeLegReductionPerDefaultArrayConstSP& feeLosses,
	const FeeLegReductionPerDefaultArrayConstSP& feeRecovered)
	:resultPathInitialized(false),
	timeline(timeline),
	name(name),
	notional(notional),
	asset(asset),
	beta(beta),
	hasDefaulted(hasDefaulted),
	defaultSettDate(defaultSettDate),
	lossGivenDefault(lossGivenDefault),
	feeLosses(feeLosses),
	feeRecovered(feeRecovered),
	recoveryOverride(recoveryOverride),
	protectionStartDateIndex(protectionStartDateIndex),
	maturityCutOffIndex(maturityCutOffIndex)
{
}

const DateTimeConstSP&
PortfolioName::IndexedSVGen::SV::getDefaultSettDate() const
{
	return defaultSettDate;
}


void
PortfolioName::IndexedSVGen::SV::initializePath(
	const CIntConstSP& pastDefaultIndex,
	int productTimelineIndex,
	int pathSize)
{
	creditLossConfigIndexedSVResultVector.resize(pathSize);
	// JCP was pathSize -1
	creditLossConfigIndexedSVResultPath.initialize(&creditLossConfigIndexedSVResultVector[0], 0, pathSize-1);

//first set losses and notionals to them to zero
		for (unsigned int i=0; i< creditLossConfigIndexedSVResultVector.size(); ++i)
		{
			creditLossConfigIndexedSVResultVector[i].lossChange = 0;
			creditLossConfigIndexedSVResultVector[i].notionalChange = 0;
			creditLossConfigIndexedSVResultVector[i].index = -1;
		}


//Logic for handling past defaults is written here. Going forward, once it is ready,
//we will call methods on PortfolioName
	if (hasDefaulted)
	{
		if (pastDefaultIndex.get())
		{
			creditLossConfigIndexedSVResultVector[pastDefaultIndex->intValue()].lossChange = lossGivenDefault;
			creditLossConfigIndexedSVResultVector[pastDefaultIndex->intValue()].notionalChange = notional;
			creditLossConfigIndexedSVResultVector[pastDefaultIndex->intValue()].index = productTimelineIndex;
		}
	}
	resultPathInitialized = true;
}


const string&
PortfolioName::IndexedSVGen::SV::getName() const
{
	return name;
}

const CreditLossConfigIndexedSVResultPath&
PortfolioName::IndexedSVGen::SV::getResults() const
{
	this->calculateLossEvents();
	return creditLossConfigIndexedSVResultPath;
}

CreditLossConfigIndexedSVResultPath&
PortfolioName::IndexedSVGen::SV::results()
{
	this->calculateLossEvents();
	return creditLossConfigIndexedSVResultPath;
}

const PortfolioNameIndexedSVSubResultArrayConstSP&
PortfolioName::IndexedSVGen::SV::getSubResults()
{
	return subResults;
}

void
PortfolioName::IndexedSVGen::SV::setSubResults(
	const PortfolioNameIndexedSVSubResultArrayConstSP&  rhsSubResults) const
{
	subResults = rhsSubResults;
}

void
PortfolioName::IndexedSVGen::SV::calculateLossEvents() const
{
	if (hasDefaulted) //as this case has been handled by the constructor
		return;

	if (!subResults.get())
		throw ModelException("PortfolioName::IndexedSVGen::SV::calculateLossEvents", "sub results not set");

//first set losses and notionals to them to zero
	for (unsigned int i=0; i< defaultIndices.size(); ++i)
	{
		creditLossConfigIndexedSVResultVector[defaultIndices[i]].lossChange = 0;
		creditLossConfigIndexedSVResultVector[defaultIndices[i]].notionalChange = 0;
		creditLossConfigIndexedSVResultVector[defaultIndices[i]].index = -1;
	}

	defaultIndices.clear();

//Logic for stochastic defaults
	for (unsigned int i=0; i< subResults->size(); ++i)
	{
		int defIndex = (*subResults)[i].defaultTimeIndex;
		if (defIndex == -1) continue; //we return -1 if the asset does not default - maybe change to DEFINE later

		int productTimeIndex = (*subResults)[i].productTimelineIndex;

		if (productTimeIndex < protectionStartDateIndex)
			continue;

		if (maturityCutOffIndex != -1)
		{
			if (productTimeIndex > maturityCutOffIndex)
				continue;
		}

		double lgd = (1 - (*subResults)[i].recovery)*notional;

		if (recoveryOverride.get())
			lgd  = lossGivenDefault;


		defaultIndices.push_back(defIndex);

		creditLossConfigIndexedSVResultVector[defIndex].lossChange = lgd;
		creditLossConfigIndexedSVResultVector[defIndex].notionalChange = notional;
		creditLossConfigIndexedSVResultVector[defIndex].index = productTimeIndex;
	}
}

int
PortfolioName::IndexedSVGen::SV::getResultPathSize() const
{
	if (resultPathInitialized)
		return creditLossConfigIndexedSVResultVector.size();
	else
		throw ModelException(
			"PortfolioName::IndexedSVGen::SV::getResultPathSize",
			"ResultPath not initialized");
}

double
PortfolioName::IndexedSVGen::SV::getNotional() const
{
	return notional;
}

const CreditAssetConstSP&
PortfolioName::IndexedSVGen::SV::getAsset() const
{
	return asset;
}

double
PortfolioName::IndexedSVGen::SV::getBeta() const
{
	return beta;
}


bool
PortfolioName::IndexedSVGen::SV::doingPast() const
{
	return false;
}

bool PortfolioNameLinkIn() {
    return PortfolioName::TYPE != NULL;
}

static const DateTime LARGEDATE = DateTime(10000000,0);

PortfolioName::NameViewSP PortfolioName::getNameView() const
{
	return NameViewSP(new NameView(
		getName(), getNameNotional(), getNameRecovery(), getProtectionStartDate(), getProtectionEndDate(LARGEDATE)));
};

// ################################################################################################################
DRLIB_END_NAMESPACE
