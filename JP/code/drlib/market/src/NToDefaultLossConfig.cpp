//----------------------------------------------------------------------------
//
//   Group       : Credit Hybrids QR
//
//   Description : An ICreditLossConfig used for N-to-default products.
//
//   Date        : October 2006
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/Maths.hpp"
#include "edginc/Atomic.hpp"
#include "edginc/Format.hpp"
#include "edginc/MarketData.hpp"
#include "edginc/IProtectionProvider.hpp"
#include "edginc/CtgLegLossPerDefault.hpp"
#include "edginc/FeeLegReductionPerDefault.hpp"
#include "edginc/IRebateCalculator.hpp"
#include "edginc/CmOnlyParameters.hpp"
#include "edginc/CreditAsset.hpp"
#include "edginc/CounterPartyCredit.hpp"
#include "edginc/IModelConfigMapper.hpp"
#include "edginc/NToDefaultLossConfig.hpp"


DRLIB_BEGIN_NAMESPACE


NToDefaultLossConfig::~NToDefaultLossConfig()
{}

NToDefaultLossConfig::NToDefaultLossConfig() : 
    CObject(TYPE), 
    defaultNumber(1)
{}


void NToDefaultLossConfig::validatePop2Object() {
    static const string method("NToDefaultLossConfig::validatePop2Object");

    if (defaultNumber < 1) {
        throw ModelException(method, "The number of the default 'N' in an "
                             "NtD must be greater than 0");
    }

    int portfolioSize = basket.size();
    if (portfolioSize < 1) {
        throw ModelException(method, "There are no names in the basket");
    }
    if (portfolioSize == 1) {
        // We could just reject this case, but it is usefull if only for testing
        // purposes (FtD on 1 name = CDS).
        // However, need to set verify that this name's beta (if present) is 
        // set to 0
        bool error = false;
        try {
            CreditEngineParametersConstSP engineParams =
                basket[0]->getEngineParams(CmOnlyParameters::TYPE);

            CmOnlyParametersConstSP cmParams(
                dynamic_cast<const CmOnlyParameters*>(engineParams.get()));

            if (!Maths::isZero(cmParams->getBeta())) {
                error = true;
            }
        }
        catch(exception e) {
            ;   // Ignore any issues here
        }
        if (error) {
            throw ModelException(method, "While baskets containing a single "
                                 "name are allowed in NtDs, the name's beta "
                                 "must be set to 0");
        }
    }

    double namesNotional = basket[0]->notional();
    for (int i=0; i < portfolioSize; ++i) {
        if (namesNotional != basket[i]->notional()) {
            throw ModelException(method, "All names in an NtD basket must "
                                 "have the same notional, but the notional of "
                                 "the first name is " +
                                 Format::toString(namesNotional) + " and the "
                                 "notional of name '" + Format::toString(i+1) + 
                                 "' is " + 
                                 Format::toString(basket[i]->notional()) + ".");
        }
    }
}


/** Returns the default number, N, of this NtD */
const double NToDefaultLossConfig::getDefaultNumber() const{
    return defaultNumber;
}


/** Returns the Nth defaultedName (ie, the one triggering this NtD) if
    available, or null otherwise */
ISingleDefaultCreditLossConfigSP NToDefaultLossConfig::getNthDefaultedName() const {
    static const string method("NToDefaultLossConfig::getNthDefaultedName");

    const DateTime& defaultDate = getDefaultDate();

    if (!defaultDate.empty()) {
        // Loop through all names until the one with the correct default date
        // is found
        for (int i=0; i < basket.size(); ++i) {
            if (defaultDate == basket[i]->getDefaultDate()) {
                return basket[i];
            }
        }
        throw ModelException(method, "Internal error: successfully identified "
                             "the default date (" + defaultDate.toString() +
                             ") but failed to find the name it corresponds to.");
    }
    // This NtD has not defaulted, so return null
    return ISingleDefaultCreditLossConfigSP();
}


//++++++++ ICreditLossConfig methods
//
/** All loss config instances have a 'name' */
string NToDefaultLossConfig::getName() const {
    return name;
}


DateTime NToDefaultLossConfig::getToday() const {
    return valueDate;
}


/** Returns an array with the cashflows corresponding to all "real"
    losses due to the settlement of this name's default.
    If there are no losses, it may return an empty (0 size)
    CtgLegLossPerDefaultArraySP or a "null" SP.
    Note: the cashflow dates are the CALCULATION DATES (as opposed to
    settlement dates) when the corresponding losses are determined. */
CtgLegLossPerDefaultArraySP NToDefaultLossConfig::historicContingentLegLosses(
    CIntConstSP triggerDelay,
    CIntConstSP defaultToCalculationDelay,
    const DateTime& lastTriggerDate,
    IBadDayAdjusterConstSP bda,
    const IProtectionProvider* const protect) const
{
    ISingleDefaultCreditLossConfigSP name = getNthDefaultedName();

    CtgLegLossPerDefaultArraySP losses;
    if (!name) {
        // This NtD has not defaulted, so there are no contingent losses.
        // Note the losses array is already set to null
    }
    else {
        // Could use "protect" to verify that this default is covered for 
        // protection. If it is not, there are several options:
        // 1. Increase defaultNumber and loop until we find one that is 
        //    covered or there are no more defaults...
        // 2. Find all defaults actually covered for protection and pick the Nth 
        //    (could be desirable for forward starting NTDs) - However, like for 
        //    tranches, this should be implemented using the PortfolioNames' 
        //    protectionStartDate field.
        // 3. If the Nth default is not covered for protection, treat it as if 
        //    the Nth default produced no losses.
        // 4. Return the losses not bothering about this - it is down to the 
        //    product using this NtD to verify whether to pay for this default 
        //    or not. The "protect" parameter is only for error verification
        //    purposes (see CDOPortfolio::validateCtgLegLosses)
        // Do (4)
        losses = name->historicContingentLegLosses(
            triggerDelay,
            defaultToCalculationDelay,
            lastTriggerDate,
            bda,
            NULL); // Pass null IProtectionProvider - verification is done here
    }

    return losses;
}


/** Compute the rebate payments caused by historic fee leg reductions */
CashFlowArraySP NToDefaultLossConfig::historicRebatePayments(
    const IRebateCalculator* const   rebateCalc,
    FeeLegReductionPerDefaultArraySP reductions,
    IForwardRatePricerSP             model,
    const bool                       recoverNotional) const
{
    static const string method("NToDefaultLossConfig::historicRebatePayments");

    // JLH Need to add support for CDS-style rebates
    throw ModelException(method, "Method not yet implemented!");

    return CashFlowArraySP();
}


/** Computes the notional reductions on the fee leg (due to losses and/or
    recovered notional) and the corresponding fee rebates, if any - based
    on the creditEventOverride in this name.
    If there are no losses, it may return an empty (0 size)
    CtgLegLossPerDefaultArraySP or a "null" SP.
    The accrualPeriods are an argument because they are required to
    determine the actual date of the notional reductions (in the fee leg,
    all notional reductions happen either on the determination date or on
    an accrual start date. */
FeeLegReductionPerDefaultArraySP NToDefaultLossConfig::historicFeeLegReductions(
    CIntConstSP triggerDelay,
    CIntConstSP defaultToCalculationDelay,
    const double temporaryLossAmount,
    const DateTime& lastTriggerDate,
    AccrualPeriodArrayConstSP accrualPeriods,
    IBadDayAdjusterConstSP bda,
    const bool recoverNotional) const
{
    FeeLegReductionPerDefaultArraySP reductions;
    ISingleDefaultCreditLossConfigSP name = getNthDefaultedName();
    if (!name) {
        // This NtD has not defaulted, so there are no losses or recovered
        // notional
    }
    else {
        // JLH Need to add support for CDS-style rebates
        reductions = name->historicFeeLegReductions(triggerDelay,
                                                    defaultToCalculationDelay,
                                                    temporaryLossAmount,
                                                    lastTriggerDate,
                                                    accrualPeriods,
                                                    bda,
                                                    recoverNotional);
    }
    return reductions;
}


/** populate from market cache */
void NToDefaultLossConfig::getMarket(const IModel* model,
                                     const MarketData* market)
{
    market->GetReferenceDate(valueDate);
    int numNames = basket.size();
    for (int i=0; i < numNames; ++i) {
        basket[i]->getMarket(model, market);
    }
}


/** Returns the engine parameters for this loss cfg (NB Model needs to
    have selected appropriate set during getMarket()). The parameter
    can specify what type of engine parameters are required. An
    exception is thrown if the relevant type is not found */
CreditEngineParametersConstSP NToDefaultLossConfig::getEngineParams(
    CClassConstSP engineParamsType) const
{
    static const string method("NToDefaultLossConfig::getEngineParams");
    try {
        if (!engineParamsType){
            throw ModelException(method, "The type of engine parameters "
                                 "requested is Null");
        }
        if (!engineParams){
            throw ModelException(method, "Null engine parameters found but "
                                 "require type: " + engineParamsType->getName());
        }

        return engineParams->getEngineParams(engineParamsType);
    }
    catch (exception& e) {
        throw ModelException(e, method, "Failed when trying to retrieve credit "
                             "engine parameters of type: " +
                             engineParamsType->getName());
    }
}


/** Returns the number of "inner loss configs" contained in this LossConfig.
    In this case cascade down the call to the basket. */
int NToDefaultLossConfig::numInnerLossConfigs() const {
    return basket.size();
}

/** Returns the inner loss config number "index". */
ICreditLossConfigConstSP NToDefaultLossConfig::getInnerLossConfig(
    const int index) const
{
    return basket[index];
}

/** Returns the maximum (realistic) loss that may be produced by this 
    loss config. Equivalent to lossGivenDefault for single names.
    "Realistic" in the sense that this method does NOT assume, eg,
    RR=0, since otherwise this method would be identical to notional() */
double NToDefaultLossConfig::maxPossibleLoss() const {
    // Need to loop through all names in the basket and pick the 
    // maximum maxLossGivenDefault
    throw ModelException("NToDefaultLossConfig::maxPossibleLoss",
                         "Internal error: this method is not yet implemented");
}


/** Returns the loss config notional, ie, the amount that could be 
    (or could have been) lost in the worst possible scenario, even if 
    this scenario is no longer possible (e.g., if a name has defaulted
    with a RR > 0 it will still be accounted as if RR=0 if the losses
    are greater in this case). 
    Note this notional bears no relation to the trade position */
double NToDefaultLossConfig::notional() const {
    // We have enforced all names to have the same notional, so the notional
    // of the NtD is now the notional of any of the underlying names
    return basket[0]->notional();
}


/** Computes the loss and recovered notional currently produced by this 
    loss config. Note both can be 0 no names have defaulted */
void NToDefaultLossConfig::currentLossAndRecovered(
    double& loss,               // (O)
    double& recoveredNotional,  // (O)
    const bool recoverNotional) const
{
    static const string method("NToDefaultLossConfig::currentLossAndRecovered");

    throw ModelException(method, "Internal error: this method is not "
                         "yet implemented");
    return;
}

/** Computes how much of this loss config's original notional has been
    prepaid (amortized) before "prepay date" */
double NToDefaultLossConfig::getPrepaidNotional(const DateTime& prepayDate) const {
    // This method is not yet fully implemented!!
    return 0.0;
}
//
//-------- ICreditLossConfig methods


//++++++++++++++++++++++++++++++++++++++++
// ISingleDefaultCreditLossConfig methods
//
/** Returns true if this loss config has defaulted */
bool NToDefaultLossConfig::hasDefaulted() const {
    int numDefaultedNames = 0;
    for (int i=0; 
         (i < basket.size()) && (numDefaultedNames < defaultNumber); 
         ++i) 
    {
        if (basket[i]->hasDefaulted()) {
            ++numDefaultedNames;
        }
    }
    return (numDefaultedNames >= defaultNumber);
}


/** Returns the default date, or an empty DateTime if the loss config
    has not defaulted */
DateTime NToDefaultLossConfig::getDefaultDate() const {
    static const string method("NToDefaultLossConfig::getDefaultDate");
    // Get and store all default dates in the defaultDates array
    DateTimeArray defaultDates(0);
    for (int i=0; i < basket.size(); ++i) {
        if (basket[i]->hasDefaulted()) {
            defaultDates.push_back(basket[i]->getDefaultDate());
        }
    }
    // If there are less than 'N' defaults, then this NtD has not 
    // defaulted yet, return an empty DateTime
    const int numOfDefaults = defaultDates.size();
    if (numOfDefaults < defaultNumber) {
        return DateTime(); 
    }

    // Sort the array and pick the Nth date (using DateTime::operator <)
    std::sort(defaultDates.begin(), defaultDates.end()); 

    // The Nth default is in element defaultNumber-1. If two or more 
    // defaults occur at the same (date)time as the Nth default, throw
    // an exception
    const DateTime& dateOfNthDefault = defaultDates[defaultNumber-1];
    if (((defaultNumber > 1) && (dateOfNthDefault == defaultDates[defaultNumber-2])) ||
         (defaultNumber <= numOfDefaults-1) && (dateOfNthDefault == defaultDates[defaultNumber]))
    {
        throw ModelException(method,
                             "Two or more defaults occured at the same time on " +
                             dateOfNthDefault.toString() +
                             ", and one of them is the default number " + 
                             Format::toString(defaultNumber) + ". The losses "
                             "associated to default number " + 
                             Format::toString(defaultNumber) +
                             " are therefore ambiguous");
    }
    return dateOfNthDefault;
}


/** Returns the timeline used in this loss config, ie, the timepoints 
    defining the regions at which it is safe to assume that the (flat 
    forwards) default rate is constant - If this concept does not make 
    sense for this loss config, a null SP should be returned.
    In this case, a union of the timelines of the underlying names in
    this NtD */
DateTimeArraySP NToDefaultLossConfig::getTimeLine() const {
    DateTimeArray mergedDates(1, valueDate);

    int numNames = basket.size();
    for (int i=0; i < numNames; ++i) {
        DateTimeArraySP nameDates = basket[i]->getTimeLine();
        if (!nameDates) {
            // This name can not produce a timeline - so forget about
            // trying to generate this timeline, since there is now way
            // it is going to be accurate anyway. Return a null SP
            return DateTimeArraySP();
        }
        mergedDates = DateTime::merge(mergedDates, *nameDates);
    }
    return DateTimeArraySP(new DateTimeArray(mergedDates));
}
//
// ISingleDefaultCreditLossConfig methods
//----------------------------------------


/** EffectiveCurveLossModelConfig::IIntoLossGen method:
    Create an IEffectiveLossCurveGen as specified by the supplied
    EffectiveLossCurveModelConfig */
ICreditLossGenSP NToDefaultLossConfig::effCurveGenerator(
    IEffectiveCurveLossModelConfigConstSP effCurveLossModelConfig,
    CounterPartyCreditConstSP             cpty,
    const bool                            recoverNotional,
    IModelConfigMapperConstSP             mapper) const
{
    return effCurveLossModelConfig->createEffCurveGenerator(
        NToDefaultLossConfigConstSP::attachToRef(this),
        cpty,
        recoverNotional);
}


/** Theta::IShift method */
bool NToDefaultLossConfig::sensShift(Theta* shift) {
    // alter immediate data
    valueDate = shift->rollDate(valueDate);
    return true; // then shift components
}


void NToDefaultLossConfig::load(CClassSP& clazz){
    clazz->setPublic(); // make visible to EAS/spreadsheet
    REGISTER(NToDefaultLossConfig, clazz);
    SUPERCLASS(CObject);
    IMPLEMENTS(ISingleDefaultCreditLossConfig);
    IMPLEMENTS(Theta::Shift);
    EMPTY_SHELL_METHOD(defaultConstructor);

    FIELD(name, "Name of the N-to-default loss config");
    FIELD(valueDate, "Value date");
    FIELD_MAKE_OPTIONAL(valueDate);
    FIELD(defaultNumber, "Number (index) of default covered for protection");
    FIELD_MAKE_OPTIONAL(defaultNumber);
//     FIELD(payFeesUntilCalcDate, "Whether fees continue to be paid after the "
//           "'defaultNumber' default, until the corresponding calculation date. "
//           "If true a rebate will be paid back when the default is settled "
//           "(on the first settlement date in case of multiple settlements). "
//           "Default: false");
//     FIELD_MAKE_OPTIONAL(payFeesUntilCalcDate);
    FIELD(basket, "Underlying loss config to apply the N-to-default to");
    FIELD(engineParams, "Engine parameters");
    FIELD_MAKE_OPTIONAL(engineParams);
}

IObject* NToDefaultLossConfig::defaultConstructor() {
    return new NToDefaultLossConfig();
}


CClassConstSP const NToDefaultLossConfig::TYPE =
    CClass::registerClassLoadMethod("NToDefaultLossConfig",
                                    typeid(NToDefaultLossConfig),
                                    load);

DEFINE_TEMPLATE_TYPE(NToDefaultLossConfigWrapper);

// Array has to have its own type
DEFINE_TEMPLATE_TYPE(NToDefaultLossConfigArray);

bool NToDefaultLossConfigLinkIn() {
    return NToDefaultLossConfig::TYPE != NULL;
}

DRLIB_END_NAMESPACE
