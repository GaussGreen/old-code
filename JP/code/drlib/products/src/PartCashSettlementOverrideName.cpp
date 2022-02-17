//----------------------------------------------------------------------------
//
//   Group       : Quatitative Research
//
//   Filename    : PartCashSettlementOverrideName.cpp
//
//   Description : Class used by credit instruments (e.g., CIS) in order 
//                 to override a name's default-related parameters, when
//                 the default is settled physically but there is also a 
//                 partial cash settlement.
//
//   Author      : Jose Hilera
//
//   Date        : February 2006
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/IBadDayAdjuster.hpp"
#include "edginc/CashSettlementOverrideName.hpp"
#include "edginc/PhysicalSettlementOverrideName.hpp"
#include "edginc/PartCashSettlementOverrideName.hpp"
#include "edginc/DeliveryDetails.hpp"
#include "edginc/ICreditContingentLeg.hpp"
#include "edginc/CtgLegLossPerDefault.hpp"
#include "edginc/FeeLegReductionPerDefault.hpp"

DRLIB_BEGIN_NAMESPACE


PartCashSettlementOverrideName::~PartCashSettlementOverrideName()
{}


PartCashSettlementOverrideName::PartCashSettlementOverrideName() : 
    DetailedCreditEventOverrideName(TYPE)
{}


/** Called immediately after object constructed */
void PartCashSettlementOverrideName::validatePop2Object() {
    static const string method(
        "PartCashSettlementOverrideName::validatePop2Object");

    // Validate fields in the parent class
    DetailedCreditEventOverrideName::validatePop2Object();

    // Create and validate the cash settlement component
    try {
        cashSettlement = CashSettlementOverrideNameSP(
            new CashSettlementOverrideName(name,
                                           eventDeterminationDate,
                                           triggerDelay,
                                           recovery,
                                           notionalFractionCash,
                                           defaultToSettlementDelayCash,
                                           calculationDate,
                                           calculationToSettlementDelay));
        cashSettlement->validatePop2Object();
    }
    catch (exception& e) {
        throw ModelException(e,
                             method,
                             "Error producing the cash settlement component of "
                             "the part cash settlement name override.");
    }


    // Create and validate the physical settlement component
    try {
        physicalSettlement = PhysicalSettlementOverrideNameSP(
            new PhysicalSettlementOverrideName(name,
                                               eventDeterminationDate,
                                               triggerDelay,
                                               recovery,
                                               NoPS,
                                               defaultToSettlementDelayPhys,
                                               deliveries,
                                               1.0 - notionalFractionCash));
        physicalSettlement->validatePop2Object();
    }
    catch (exception& e) {
        throw ModelException(e,
                             method,
                             "Error producing the physical settlement component "
                             "of the part cash settlement name override.");
    }
}


/** Returns the portfolio reductions for the fee leg (corresponding to losses
    and recovered amounts). */
CtgLegLossPerDefaultArraySP 
PartCashSettlementOverrideName::historicContingentLegLosses(
    const double nameNotional,
    const double recoveryRate,
    const DateTime& lastTriggerDate,
    const DateTime& creditEventDate,
    IBadDayAdjusterConstSP bda) const
{
    static const string method( 
        "PartCashSettlementOverrideName::historicContingentLegLosses");

    try {
        // Get the cash settlement cashflows
        CtgLegLossPerDefaultArraySP cashReductions;
        if (!Maths::isZero(notionalFractionCash)) {
            cashReductions = 
                cashSettlement->historicContingentLegLosses(
                    nameNotional,
                    recoveryRate,
                    lastTriggerDate,
                    creditEventDate,
                    bda);
        }

        // Get the physical settlement cashflows
        CtgLegLossPerDefaultArraySP physicalReductions;
        if (!Maths::equals(notionalFractionCash, 1.0)) {
            physicalReductions = 
                physicalSettlement->historicContingentLegLosses(
                    nameNotional,
                    recoveryRate,
                    lastTriggerDate,
                    creditEventDate,
                    bda);
        }

        // Return the result of merging both arrays
        if (!cashReductions) {
            return physicalReductions;
        }
        if (!physicalReductions) {
            return cashReductions;
        }
        // we have both -> merge them. Arbitrarily, add them to cashReductions
        int numPhysicalReductions = physicalReductions->size();
        cashReductions->reserve(cashReductions->size() + numPhysicalReductions);
        for (int i=0; i < numPhysicalReductions; ++i) {
            cashReductions->push_back((*physicalReductions)[i]);
        }
        return cashReductions;
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }
}


/** Returns the accrual pay date. For part cash settlements this is the 
 * earlies of the accrual dates for the physical and cash components */
const DateTime PartCashSettlementOverrideName::getAccrualPayDate(
    const DateTime& creditEventDate,
    const DateTime& lastTriggerDate,
    const DateTime& valueDate,
    IBadDayAdjusterConstSP bda) const
{
    if (notionalFractionCash == 0.0) {
        return physicalSettlement->getAccrualPayDate(
            creditEventDate, lastTriggerDate, valueDate, bda);
    }
    else if (notionalFractionCash == 1.0) {
        return cashSettlement->getAccrualPayDate(
            creditEventDate, lastTriggerDate, valueDate, bda);
    }
    else {
        const DateTime& physAccrualPayDate = physicalSettlement->getAccrualPayDate(
            creditEventDate, lastTriggerDate, valueDate, bda);

        const DateTime& cashAccrualPayDate = cashSettlement->getAccrualPayDate(
            creditEventDate, lastTriggerDate, valueDate, bda);

        return physAccrualPayDate.min(cashAccrualPayDate);
    }
}



/* Returns the overall recovery rate accross all settlements, such that
 * total loss = overall defaulted notional * (1-overall recovery rate) */
double PartCashSettlementOverrideName::getOverallRecoveryRate(
    double recoveryRate) const
{
    static const string method(
        "PartCashSettlementOverrideName::getOverallRecoveryRate");

    try { 
        const double cashDefaultedNotional = 
            cashSettlement->getOverallDefaultedNotionalFraction();
        const double cashRecovery = 
            cashSettlement->getOverallRecoveryRate(recoveryRate);

        const double physicalDefaultedNotional = 
            physicalSettlement->getOverallDefaultedNotionalFraction();
        const double physicalRecovery = 
            physicalSettlement->getOverallRecoveryRate(recoveryRate);

        const double totalLoss = (cashDefaultedNotional * (1-cashRecovery)) +
            (physicalDefaultedNotional * (1-physicalRecovery));

        const double overallRecovery = 1 - 
            (totalLoss / (cashDefaultedNotional + physicalDefaultedNotional));

        return overallRecovery;
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }
}


/** Returns the percentage of the notional that has been triggered 
 * for protection */
double PartCashSettlementOverrideName::getOverallDefaultedNotionalFraction() const
{
    static const string method(
        "PartCashSettlementOverrideName::getOverallDefaultedNotionalFraction");
    try { 
        const double cashDefaultedNotional = 
            cashSettlement->getOverallDefaultedNotionalFraction();
        const double physicalDefaultedNotional = 
            physicalSettlement->getOverallDefaultedNotionalFraction();

        const double totalDefaultedNotional = 
            cashDefaultedNotional + physicalDefaultedNotional;

        return totalDefaultedNotional;
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }
}


IObject* PartCashSettlementOverrideName::defaultPartCashSettlementOverrideName() {
    return new PartCashSettlementOverrideName();
}


void PartCashSettlementOverrideName::load (CClassSP& clazz) {
    clazz->setPublic(); // make visible to EAS/spreadsheet
    REGISTER(PartCashSettlementOverrideName, clazz);
    SUPERCLASS(DetailedCreditEventOverrideName);
    EMPTY_SHELL_METHOD(defaultPartCashSettlementOverrideName);

    FIELD(notionalFractionCash,         "Percentage of the notional being "
                                               "cash-settled (eg, 0.5 for 50%)");
    FIELD(defaultToSettlementDelayCash,        "Delay between credit event and "
                                               "settlement");
    FIELD(calculationDate,              "Date when the recovery rate is "
                                               "determined");
    FIELD(calculationToSettlementDelay, "Delay between calculation date "
                                               "and settlement");
    FIELD(NoPS,                                "Specification of the amounts "
                                               "that will be physically settled");
    FIELD(defaultToSettlementDelayPhys,        "Delay between credit event and "
                                               "settlement");
    FIELD(deliveries,                          "Details for the scheduled "
                                               "physical deliveries");
    FIELD(physicalSettlement,                  "Transient: Physical settlement");
    FIELD(cashSettlement,                      "Transient: Cash settlement");

    // notionalFractionCash is optional in CashSettlementOverrideName but not
    // on a PartCashSettlementOverrideName
    // FIELD_MAKE_OPTIONAL(notionalFractionCash);
    FIELD_MAKE_OPTIONAL(defaultToSettlementDelayCash);
    FIELD_MAKE_OPTIONAL(calculationDate);
    FIELD_MAKE_OPTIONAL(NoPS);
    FIELD_MAKE_OPTIONAL(defaultToSettlementDelayPhys);
    FIELD_MAKE_OPTIONAL(deliveries);
    FIELD_MAKE_TRANSIENT(physicalSettlement);
    FIELD_MAKE_TRANSIENT(cashSettlement);
}


CClassConstSP const PartCashSettlementOverrideName::TYPE = 
    CClass::registerClassLoadMethod("PartCashSettlementOverrideName", 
                                    typeid(PartCashSettlementOverrideName), 
                                    load);

// Array has to have its own type
DEFINE_TEMPLATE_TYPE(PartCashSettlementOverrideNameArray);


/** Included in ProductsLib-modified::linkInClasses() via the productSrcsMap
 * script to force the linker to include this file */
bool PartCashSettlementOverrideNameLoad() {
    return (PartCashSettlementOverrideName::TYPE != 0);
}

DRLIB_END_NAMESPACE
