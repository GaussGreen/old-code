//----------------------------------------------------------------------------
//
//   Group       : Credit Hybrids QR
//
//   Filename    : TrancheCreditEventOverride.cpp
//
//   Description : Simple class used by CDO in 
//                 order to override a name's default-related parameters.
//
//   Author      : Jose Hilera
//
//   Date        : April 2006
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/Format.hpp"
#include "edginc/TrancheCreditEventOverride.hpp"
#include "edginc/IRebateCalculator.hpp"
#include "edginc/IBadDayAdjuster.hpp"

DRLIB_BEGIN_NAMESPACE


TrancheCreditEventOverride::TrancheCreditEventOverride() : 
    CObject(TYPE)
{}

TrancheCreditEventOverride::TrancheCreditEventOverride(
    CClassConstSP clazz) : CObject(clazz)
{}

/** Constructor with all fields passed in */
TrancheCreditEventOverride::TrancheCreditEventOverride(
    CClassConstSP clazz,
    CIntSP triggerDelay,
    const DateTime& eventDeterminationDate,
    CDoubleSP recovery,
    const DateTime& valueDate) :
        CObject(clazz),
        triggerDelay(triggerDelay),
        eventDeterminationDate(eventDeterminationDate),
        recovery(recovery),
        valueDate(valueDate)
{}

TrancheCreditEventOverride::~TrancheCreditEventOverride()
{}


/** Called immediately after object constructed */
void TrancheCreditEventOverride::validatePop2Object() {
    static const string method = 
        "TrancheCreditEventOverride::validatePop2Object";
    
    if (!!recovery) {
        double dRecovery = recovery->doubleValue();
        if (dRecovery < 0.0 || dRecovery > 1.0) {
            throw ModelException(method,
                                 "Recovery out of [0,1] (recovery = " +
                                 Format::toString(dRecovery) + ")");
        }
    }
    
    // Either the eventDeterminationDate or the triggerDelay need to
    // be provided. EventDeterminationDate will be preferred if provided
    if (eventDeterminationDate.empty()) {
        if (!triggerDelay) {
            throw ModelException(method,
                                 "Need to enter eventDeterminationDate or "
                                 "triggerDelay");
        }
        else {
            // Will use the triggerDelay, so verify it is right
            if (triggerDelay->intValue() < 0) {
                throw ModelException(method,
                                     "The triggerDelay is negative, and "
                                     "negative delays are not accepted.");
            }
        }
    }
}


/** GetMarket implementation*/
void TrancheCreditEventOverride::getMarket(const IModel* model, 
                                           const MarketData* market)
{
    market->GetReferenceDate(valueDate);
}


/** Returns the event determination date with the information in this 
 * override - it can be an estimate or the actual eventDeterminationDate,
 * if set. 
 * If the eventDeterminationDate is set but it is after lastTriggerDate, or
 * if it cannot be estimated because it would be after lastTriggerDate, an
 * empty date is returned indicating that this CDS will not be triggered. */
const DateTime TrancheCreditEventOverride::getEventDeterminationDate(
    const DateTime& creditEventDate,
    const DateTime& lastTriggerDate,
    IBadDayAdjusterConstSP bda) const
{
    static const string method = 
        "TrancheCreditEventOverride::getEventDeterminationDate";

    DateTime determinationDate;
    if (eventDeterminationDate.empty()) {
        determinationDate = ITrancheCreditEventOverride::rollAndAdjustDate(
            creditEventDate,
            triggerDelay,
            valueDate,
            valueDate,
            lastTriggerDate,
            bda);
    }
    else {
#ifndef QLIB_REMOVE_EDD_VALIDATION_CDO
        if (creditEventDate > eventDeterminationDate) {
            throw ModelException(method,
                                 "Determination date (" +
                                 eventDeterminationDate.toString() +
                                 ") can not be before the default date ("+
                                 creditEventDate.toString() + ").");
        }
#endif
        if (!lastTriggerDate.empty() &&
            (eventDeterminationDate > lastTriggerDate))
        {
            // determinationDate is already an empty DateTime
        }
        else {
            // Return eventDeterminationDate
            determinationDate = eventDeterminationDate;
        }
    }
    return determinationDate;
}


void TrancheCreditEventOverride::load (CClassSP& clazz) {
    clazz->setPublic(); // make visible to EAS/spreadsheet
    REGISTER(TrancheCreditEventOverride, clazz);
    SUPERCLASS(CObject);
    IMPLEMENTS(ITrancheCreditEventOverride);
    IMPLEMENTS(IGetMarket);

    FIELD(eventDeterminationDate, "Date when fees and accrual stop");
    FIELD_MAKE_OPTIONAL(eventDeterminationDate);

    FIELD(triggerDelay, "Delay between default and eventDeterminationDate - "
          "used to estimate eventDeterminationDate before it is known.");
    FIELD_MAKE_OPTIONAL(triggerDelay);

    FIELD(recovery, "Recovery override");
    FIELD_MAKE_OPTIONAL(recovery);

    FIELD(valueDate, "Value date");
    FIELD_MAKE_OPTIONAL(valueDate);
}


CClassConstSP const TrancheCreditEventOverride::TYPE = 
    CClass::registerClassLoadMethod("TrancheCreditEventOverride", 
                                    typeid(TrancheCreditEventOverride), 
                                    load);

// Array has to have its own type
DEFINE_TEMPLATE_TYPE(TrancheCreditEventOverrideArray);


/** Included in ProductsLib-modified::linkInClasses() via the productSrcsMap
 * script to force the linker to include this file */
bool TrancheCreditEventOverrideLoad() {
    return (TrancheCreditEventOverride::TYPE != 0);
}

DRLIB_END_NAMESPACE
