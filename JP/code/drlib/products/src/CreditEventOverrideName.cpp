//----------------------------------------------------------------------------
//
//   Group       : Quatitative Research
//
//   Filename    : CreditEventOverrideName.cpp
//
//   Description : Simple class used by credit instruments (e.g., CIS) in order 
//                 to override (at instrument level) a name's default-related 
//                 parameters.
//
//   Author      : Jose Hilera
//
//   Date        : December 2005
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/Format.hpp"
#include "edginc/CreditEventOverrideName.hpp"
#include "edginc/ICreditFeeLeg.hpp"
#include "edginc/ICreditContingentLeg.hpp"
#include "edginc/IBadDayAdjuster.hpp"
#include "edginc/IForwardRatePricer.hpp"
#include "edginc/FeeLegReductionPerDefault.hpp"
#include "edginc/CtgLegLossPerDefault.hpp"


DRLIB_BEGIN_NAMESPACE


CreditEventOverrideName::CreditEventOverrideName() : CObject(TYPE)
{}


CreditEventOverrideName::~CreditEventOverrideName()
{}


/** Called immediately after object constructed */
void CreditEventOverrideName::validatePop2Object() {
    static const string method = "CreditEventOverrideName::validatePop2Object";
    
    if (!!recovery) {
        double dRecovery = recovery->doubleValue();
        if (dRecovery < 0.0 || dRecovery > 1.0) {
            throw ModelException(method,
                                 "Recovery out of [0,1] (recovery = " +
                                 Format::toString(dRecovery) + ")");
        }
    }
    if (eventDeterminationDate > defaultPayDate) {
        throw ModelException(method,
                             "DefaultPayDate (" + defaultPayDate.toString() +
                             ") happens before eventDeterminationDate (" +
                             eventDeterminationDate.toString() + ").");
    }
}


/** Name of this market object */
string CreditEventOverrideName::getName() const {
    return name.getName();
}


FeeLegReductionPerDefaultArraySP 
CreditEventOverrideName::historicFeeLegReductions(
    const double notional,
    const double recoveryRate,
    const DateTime& lastTriggerDate,
    const DateTime& creditEventDate,
    IBadDayAdjusterConstSP bda) const
{
    static const string method(
        "CreditEventOverrideName::historicFeeLegReductions");

    if (creditEventDate > eventDeterminationDate) {
        throw ModelException(method,
                             "The event determination date (" +
                             eventDeterminationDate.toString() +
                             ") can not be before the default date ("+
                             creditEventDate.toString() + ").");
    }

    FeeLegReductionPerDefaultArraySP feeReductions;
    if (eventDeterminationDate > lastTriggerDate) {
        // The default was triggered after the lastTriggerDate, so no fee 
        // reductions should be paid.
        // feeReductions is already a NULL sp
    }
    else {
        // Check if we are overriding the recovery rate here
        double recRate = getOverallRecoveryRate(recoveryRate);

        // Since we know the actual determination date, use it (and ignore the
        // lastTriggerDate and creditEventDate)
        feeReductions = FeeLegReductionPerDefaultArraySP(
            new FeeLegReductionPerDefaultArray(
                1, FeeLegReductionPerDefaultSP(new FeeLegReductionPerDefault(
                    eventDeterminationDate,
                    defaultPayDate,
                    defaultPayDate,
                    notional,
                    notional,
                    recRate))));
    }
    return feeReductions;
}


/** Returns the contingent leg loss cashflows corresponding to a defaulted 
    name, assuming that the default settles in one date. */
CtgLegLossPerDefaultArraySP
CreditEventOverrideName::historicContingentLegLosses(
    const double nameNotional,
    const double recoveryRate,
    const DateTime& lastTriggerDate,
    const DateTime& creditEventDate,
    IBadDayAdjusterConstSP bda) const
{
    static const string method( 
        "CreditEventOverrideName::historicContingentLegLosses");

    if (creditEventDate > eventDeterminationDate) {
        throw ModelException(method,
                             "The event determination date (" +
                             eventDeterminationDate.toString() +
                             ") can not be before the default date ("+
                             creditEventDate.toString() + ").");
    }

    CtgLegLossPerDefaultArraySP ctgReductions;
    if (eventDeterminationDate > lastTriggerDate) {
        // The default was triggered after the lastTriggerDate, so no
        // contingent payment will be made.
        // ctgReductions is already a NULL sp.
    }
    else {
        // Check if we are overriding the recovery rate here
        double recRate = getOverallRecoveryRate(recoveryRate);

        // Since we know the actual determination date, use it (and ignore the
        // lastTriggerDate and creditEventDate)
        ctgReductions = CtgLegLossPerDefaultArraySP(new CtgLegLossPerDefaultArray(
            1, CtgLegLossPerDefaultSP(new CtgLegLossPerDefault(
                creditEventDate,
                defaultPayDate,
                nameNotional * (1.0 - recRate)))));
    }
    return ctgReductions;
}


/* Returns the overall recovery rate accross all settlements, such that
 * total loss = overall defaulted notional * (1-overall recovery rate) */
double CreditEventOverrideName::getOverallRecoveryRate(
    double recoveryRate) const
{
    return (!recovery) ? 
        recoveryRate : // No override, use method argument
        recovery->doubleValue();
}

/** Returns the percentage of the notional that has been triggered 
 * for protection */
double CreditEventOverrideName::getOverallDefaultedNotionalFraction() const {
    return 1.0; // 100%
}


IObject* CreditEventOverrideName::defaultCreditEventOverrideName() {
    return new CreditEventOverrideName();
}


/** GetMarket implementation*/
void CreditEventOverrideName::getMarket(const IModel* model, 
                                        const MarketData* market)
{
    // nothing to do here
}


void CreditEventOverrideName::load (CClassSP& clazz) {
    clazz->setPublic(); // make visible to EAS/spreadsheet
    REGISTER(CreditEventOverrideName, clazz);
    SUPERCLASS(CObject);
    IMPLEMENTS(ICreditEventOverrideName);
    EMPTY_SHELL_METHOD(defaultCreditEventOverrideName);

    FIELD(name,                   "Name of the CDS Par Spreads this "
                                  "override applies to");
    FIELD(eventDeterminationDate, "Date when fees and accrual stop");
    FIELD(defaultPayDate,         "Date when fees and (1-R) are paid");
    FIELD(recovery,               "Recovery override");
    FIELD_MAKE_OPTIONAL(recovery);
}


CClassConstSP const CreditEventOverrideName::TYPE = 
    CClass::registerClassLoadMethod("CreditEventOverrideName", 
                                    typeid(CreditEventOverrideName), 
                                    load);

// Array has to have its own type
DEFINE_TEMPLATE_TYPE(CreditEventOverrideNameArray);


/** Included in ProductsLib-modified::linkInClasses() via the productSrcsMap
 * script to force the linker to include this file */
bool CreditEventOverrideNameLoad() {
    return (CreditEventOverrideName::TYPE != 0);
}

DRLIB_END_NAMESPACE
