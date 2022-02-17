//----------------------------------------------------------------------------
//
//   Group       : Credit Hybrids QR
//
//   Filename    : DetailedCreditEventOverrideName.cpp
//
//   Description : Simple class used by credit instruments (e.g., CIS) in 
//                 order to override a name's default-related parameters.
//
//   Author      : Jose Hilera
//
//   Date        : February 2006
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/Format.hpp"
#include "edginc/DetailedCreditEventOverrideName.hpp"
#include "edginc/ICreditFeeLeg.hpp"
#include "edginc/IBadDayAdjuster.hpp"
#include "edginc/FeeLegReductionPerDefault.hpp"

DRLIB_BEGIN_NAMESPACE


DetailedCreditEventOverrideName::DetailedCreditEventOverrideName() : 
    CObject(TYPE)
{}

DetailedCreditEventOverrideName::DetailedCreditEventOverrideName(
    CClassConstSP clazz) : CObject(clazz)
{}

/** Constructor with all fields passed in */
DetailedCreditEventOverrideName::DetailedCreditEventOverrideName(
        CClassConstSP clazz,
        ICDSParSpreadsWrapper name,
        DateTime eventDeterminationDate,
        CIntSP triggerDelay,
        CDoubleSP recovery) :
    CObject(clazz),
    name(name),
    eventDeterminationDate(eventDeterminationDate),
    triggerDelay(triggerDelay),
    recovery(recovery)
{}

DetailedCreditEventOverrideName::~DetailedCreditEventOverrideName()
{}


/** Called immediately after object constructed */
void DetailedCreditEventOverrideName::validatePop2Object() {
    static const string method = 
        "DetailedCreditEventOverrideName::validatePop2Object";
    
    if (!!recovery) {
        double dRecovery = recovery->doubleValue();
        if (dRecovery < 0.0 || dRecovery > 1.0) {
            throw ModelException(method,
                                 "Recovery of name " + name.getName() +
                                 " out of [0,1] (recovery = " +
                                 Format::toString(dRecovery) + ")");
        }
    }

    // The eventDeterminationDate or the triggerDelay need to
    // be provided. EventDeterminationDate will be preferred if provided
    if (eventDeterminationDate.empty()) {
        if (!triggerDelay) {
            throw ModelException(method,
                                 "Need to enter eventDeterminationDate or "
                                 "triggerDelay for name " +
                                 name.getName());
        }
        else {
            // Will use the triggerDelay, so verify it is right
            if (triggerDelay->intValue() < 0) {
                throw ModelException(method,
                                     "The triggerDelay for name " +
                                     name.getName() +
                                     " is negative, and negative delays are "
                                     "not accepted.");
            }
        }
    }
}


/** Name of the curve associated to this override */
string DetailedCreditEventOverrideName::getName() const {
    return name.getName();
}


/** Returns the event determination date with the information in this 
 * override - it can be an estimate or the actual eventDeterminationDate,
 * if set. 
 * If the eventDeterminationDate is set but it is after lastTriggerDate, or
 * if it cannot be estimated because it would be after lastTriggerDate, an
 * empty date is returned indicating that this CDS will not be triggered. */
const DateTime DetailedCreditEventOverrideName::getEventDeterminationDate(
    const DateTime& creditEventDate,
    const DateTime& lastTriggerDate,
    IBadDayAdjusterConstSP bda) const
{
    static const string method(
        "DetailedCreditEventOverrideName::getEventDeterminationDate");

    DateTime determinationDate;
    if (eventDeterminationDate.empty()) {
        determinationDate = rollAndAdjustDate(creditEventDate,
                                              triggerDelay,
                                              valueDate,
                                              valueDate,
                                              lastTriggerDate,
                                              bda);
    }
    else {
#ifndef QLIB_REMOVE_EDD_VALIDATION_CDS
        if (creditEventDate > eventDeterminationDate) {
            throw ModelException(method,
                                 name.getName() + 
                                 "'s event determination date (" +
                                 eventDeterminationDate.toString() +
                                 ") can not be before the default date ("+
                                 creditEventDate.toString() + ").");
        }
#endif
        if (eventDeterminationDate > lastTriggerDate) {
            // determinationDate is already an empty DateTime
        }
        else {
            // Return eventDeterminationDate
            determinationDate = eventDeterminationDate;
        }
    }
    return determinationDate;
}


FeeLegReductionPerDefaultArraySP 
DetailedCreditEventOverrideName::historicFeeLegReductions(
    const double notional,
    const double recoveryRate,
    const DateTime& lastTriggerDate,
    const DateTime& creditEventDate,
    IBadDayAdjusterConstSP bda) const
{
    static const string method(
        "DetailedCreditEventOverrideName::getFeeLegCashFlows");

    DateTime accrualPayDate;
    DateTime determinationDate = getEventDeterminationDate(creditEventDate, 
                                                           lastTriggerDate,
                                                           bda);
    
    if (determinationDate.empty()) {
        // This means that the default cannot be triggered. Therefore do not
        //  bother computing the accrual pay date, since it will not be used.
        //    
        // accrualPayDate is already an empty date
    }
    else {
        accrualPayDate = getAccrualPayDate(creditEventDate,
                                           lastTriggerDate,
                                           valueDate,
                                           bda);
        if (accrualPayDate.empty()) {
            // The default will/has been triggered on determinationDate but
            // we failed to determine when the accrual will be paid -> throw
            // an exception
            throw ModelException (method,
                                  "Internal error: determinationDate is " +
                                  determinationDate.toString() +
                                  " but failed to determine the accrualPayDate "
                                  "for name " + name.getName());
        }
    }

    // Obtain the percentage of the notional that has been triggered 
    // for protection
    const double defaultedNotionalFraction = getOverallDefaultedNotionalFraction();
    const double defaultedNotional = notional * defaultedNotionalFraction;

    FeeLegReductionPerDefaultArraySP feeReductions(
        new FeeLegReductionPerDefaultArray(
            1, FeeLegReductionPerDefaultSP(new FeeLegReductionPerDefault(
                determinationDate,
                accrualPayDate, // effective date
                accrualPayDate, // calculation date
                defaultedNotional,
                notional,
                getOverallRecoveryRate(recoveryRate)))));

    return feeReductions;
}


/* Returns the overriden recovery rate */
double DetailedCreditEventOverrideName::getOverridenRecoveryRate(
    double recoveryRate) const
{
    return (!recovery) ? 
        recoveryRate : // No override, use method argument
        recovery->doubleValue();
}


/** GetMarket implementation*/
void DetailedCreditEventOverrideName::getMarket(const IModel* model, 
                                                const MarketData* market)
{
    market->GetReferenceDate(valueDate);
}

void DetailedCreditEventOverrideName::load (CClassSP& clazz) {
    clazz->setPublic(); // make visible to EAS/spreadsheet
    REGISTER(DetailedCreditEventOverrideName, clazz);
    SUPERCLASS(CObject);
    IMPLEMENTS(ICreditEventOverrideName);
    IMPLEMENTS(IGetMarket);

    FIELD(name, "Name of the CDS Par Spreads this override applies to");
    FIELD(eventDeterminationDate, "Date when fees and accrual stop");
    FIELD_MAKE_OPTIONAL(eventDeterminationDate);
    FIELD(triggerDelay, "Delay, in days, between default and "
          "eventDeterminationDate - used to estimate eventDeterminationDate "
          "before it is known.");
    FIELD_MAKE_OPTIONAL(triggerDelay);
    FIELD(recovery, "Recovery override");
    FIELD_MAKE_OPTIONAL(recovery);
    FIELD(valueDate, "ValueDate");
    FIELD_MAKE_OPTIONAL(valueDate);
}


CClassConstSP const DetailedCreditEventOverrideName::TYPE = 
    CClass::registerClassLoadMethod("DetailedCreditEventOverrideName", 
                                    typeid(DetailedCreditEventOverrideName), 
                                    load);

// Array has to have its own type
DEFINE_TEMPLATE_TYPE(DetailedCreditEventOverrideNameArray);


/** Included in ProductsLib-modified::linkInClasses() via the productSrcsMap
 * script to force the linker to include this file */
bool DetailedCreditEventOverrideNameLoad() {
    return (DetailedCreditEventOverrideName::TYPE != 0);
}

DRLIB_END_NAMESPACE
