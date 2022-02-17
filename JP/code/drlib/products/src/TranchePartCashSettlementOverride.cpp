//----------------------------------------------------------------------------
//
//   Group       : Quatitative Research
//
//   Filename    : TranchePartCashSettlementOverride.cpp
//
//   Description : Class used by CDO in order to override a name's default-
//                 related parameters, when the default is settled physically
//                 with a part cash settlement
//
//   Author      : Jose Hilera
//
//   Date        : June 2006
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/Format.hpp"
#include "edginc/DeliveryDetails.hpp"
#include "edginc/FeeLegReductionPerDefault.hpp"
#include "edginc/TrancheCashSettlementOverride.hpp"
#include "edginc/TranchePhysicalSettlementOverride.hpp"
#include "edginc/TranchePartCashSettlementOverride.hpp"

DRLIB_BEGIN_NAMESPACE


TranchePartCashSettlementOverride::~TranchePartCashSettlementOverride()
{}


TranchePartCashSettlementOverride::TranchePartCashSettlementOverride() : 
    TrancheCreditEventOverride(TYPE)
{}


/** Called immediately after object constructed */
void TranchePartCashSettlementOverride::validatePop2Object() {
    static const string method = 
        "TranchePartCashSettlementOverride::validatePop2Object";

    // Validate fields in the parent class
    TrancheCreditEventOverride::validatePop2Object();
    
    // Create and validate the cash settlement component
    try {
        cashSettlement.reset(new TrancheCashSettlementOverride(
            triggerDelay,
            eventDeterminationDate,
            !cashRecovery ? recovery : cashRecovery,
            valueDate,
            cashDefaultedNotionalFraction,
            cashDefaultToCalculationDelay,
            cashCalculationDate,
            cashNotionalFraction));
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
        physicalSettlement.reset(new TranchePhysicalSettlementOverride(
            triggerDelay,
            eventDeterminationDate,
            recovery,
            valueDate,
            NoPS,
            deliveries,
            settleDeliveriesSeparately,
            physicalDefaultToCalculationDelay,
            physicalCalculationDate,
            cutoffDate,
            moreDeliveriesPending,
            1.0 - cashNotionalFraction));
        physicalSettlement->validatePop2Object();
    }
    catch (exception& e) {
        throw ModelException(e,
                             method,
                             "Error producing the physical settlement component "
                             "of the part cash settlement name override.");
    }
    return;
}


/** Returns the contingent leg cashflows corresponding to a defaulted name.
 * This object contains the credit event parameters */
CtgLegLossPerDefaultArraySP 
TranchePartCashSettlementOverride::historicContingentLegLosses(
    const double nameNotional,
    const double recoveryRate,
    const DateTime& lastTriggerDate,
    const DateTime& creditEventDate,
    IBadDayAdjusterConstSP bda) const
{
    static const string method(
        "TranchePartCashSettlementOverride::historicContingentLegLosses");

    try { 
        // Initialise arrays of cash and physical settlement contingent losses
        CtgLegLossPerDefaultArraySP cashContingentLosses = 
            cashSettlement->historicContingentLegLosses(
                 nameNotional,
                 recoveryRate,
                 lastTriggerDate,
                 creditEventDate,
                 bda);
        CtgLegLossPerDefaultArraySP physicalContingentLosses = 
            physicalSettlement->historicContingentLegLosses(
                 nameNotional,
                 recoveryRate,
                 lastTriggerDate,
                 creditEventDate,
                 bda);

        // Merge the arrays and return
        CtgLegLossPerDefaultArraySP totalLosses;
        if (!!cashContingentLosses) {
            // copy the cash Losses
            totalLosses = cashContingentLosses;
            // and add the physical losses
            if (!!physicalContingentLosses) {
                unsigned int numPhysicalLosses = physicalContingentLosses->size();
                totalLosses->reserve(totalLosses->size() + numPhysicalLosses);
                for (unsigned int i=0; i < numPhysicalLosses; ++i) {
                    totalLosses->push_back((*physicalContingentLosses)[i]);
                }
            }
        }
        else {
            // there are no cash Losses, so just copy the physical losses
            totalLosses = physicalContingentLosses;
        }
        
        return totalLosses;
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }
}



/** Returns the portfolio reductions for the fee leg (corresponding to losses
 * and recovered amounts). */
FeeLegReductionPerDefaultArraySP 
TranchePartCashSettlementOverride::historicFeeLegReductions(
    const double nameNotional,
    const double recoveryRate,
    const DateTime& lastTriggerDate,
    const DateTime& creditEventDate,
    IBadDayAdjusterConstSP bda) const
{
    static const string method(
        "TranchePartCashSettlementOverride::historicFeeLegReductions");

    try { 
        // Get the cash and physical settlement contingent fee reductions
    
        FeeLegReductionPerDefaultArraySP cashReductions =
            cashSettlement->historicFeeLegReductions(nameNotional,
                                                     recoveryRate,
                                                     lastTriggerDate,
                                                     creditEventDate,
                                                     bda);

        FeeLegReductionPerDefaultArraySP physicalReductions = 
            physicalSettlement->historicFeeLegReductions(nameNotional,
                                                         recoveryRate,
                                                         lastTriggerDate,
                                                         creditEventDate,
                                                         bda);

        if (!cashReductions) {
            return physicalReductions;
        }
        if (!physicalReductions) {
            return cashReductions;
        }
        // we have both -> merge them. Arbitrarily, add them to cashReductions
        unsigned int numPhysicalReductions = physicalReductions->size();
        cashReductions->reserve(cashReductions->size() + numPhysicalReductions);
        for (unsigned int i=0; i < numPhysicalReductions; ++i) {
            cashReductions->push_back((*physicalReductions)[i]);
        }
        return cashReductions;
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }
}


/* Returns the overall recovery rate accross all settlements, such that
 * total loss = overall defaulted notional * (1-overall recovery rate) */
double TranchePartCashSettlementOverride::getOverallRecoveryRate(
    double recoveryRate) const 
{
    static const string method(
        "TranchePartCashSettlementOverride::getOverallRecoveryRate");

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
double TranchePartCashSettlementOverride::getOverallDefaultedNotionalFraction() const {
    static const string method = 
        "TranchePartCashSettlementOverride::getOverallDefaultedNotionalFraction";

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


IObject* TranchePartCashSettlementOverride::defaultTranchePartCashSettlementOverride() {
    return new TranchePartCashSettlementOverride();
}


void TranchePartCashSettlementOverride::load (CClassSP& clazz) {
    clazz->setPublic(); // make visible to EAS/spreadsheet
    REGISTER(TranchePartCashSettlementOverride, clazz);
    SUPERCLASS(TrancheCreditEventOverride);
    EMPTY_SHELL_METHOD(defaultTranchePartCashSettlementOverride);

    FIELD(NoPS, "Notice of Physical Settlement: specification of the amounts "
          "being settled");
    FIELD_MAKE_OPTIONAL(NoPS);
    FIELD(deliveries, "Array with the details of all deliveries");
    FIELD_MAKE_OPTIONAL(deliveries);
    FIELD(settleDeliveriesSeparately, "If true, each delivery will "
          "have its own calculation and settlement date; otherwise all "
          "deliveries will share a common calculation and settlement "
          "date");
    FIELD_MAKE_OPTIONAL(settleDeliveriesSeparately);
    FIELD(physicalDefaultToCalculationDelay, "An estimate for the delay between "
          "a credit event and the calculation date for all obligations settling "
          "physically (if field settleDeliveriesSeparately is false) or for the "
          "remaining expected obligations (if true)");
    FIELD_MAKE_OPTIONAL(physicalDefaultToCalculationDelay);
    FIELD(physicalCalculationDate, "Date when the recovery rate for the "
                 "obligations settling physically is determined "
                 "if NOT settling deliveries separately");
    FIELD_MAKE_OPTIONAL(physicalCalculationDate);
    FIELD(cutoffDate, "The last date when deliveries are accepted. "
                 "If the 'Recover notional' flag is set, any outstanding "
                 "notional not delivered on or before this date will "
                 "be deemed recovered.");
    FIELD_MAKE_OPTIONAL(cutoffDate);
    FIELD(moreDeliveriesPending, "Indicates whether more deliveries are "
                 "to be expected before the cut-off date.");
    FIELD_MAKE_OPTIONAL(moreDeliveriesPending);

    FIELD(cashDefaultedNotionalFraction, 
          "Percentage of the notional triggered for protection in cash. "
          "Must be smaller than the cashNotionalFraction.");
    FIELD_MAKE_OPTIONAL(cashDefaultedNotionalFraction);
    FIELD(cashDefaultToCalculationDelay, 
          "Delay between credit event and cash settlement");
    FIELD_MAKE_OPTIONAL(cashDefaultToCalculationDelay);
    FIELD(cashCalculationDate, 
                 "Date when the recovery rate is determined for the "
                 "obligations being cash-settled.");
    FIELD_MAKE_OPTIONAL(cashCalculationDate);
    FIELD(cashNotionalFraction,
                 "Percentage of the notional considered for cash settlement");
    FIELD(cashRecovery, "Recovery rate for the cash settlement. If not "
          "present the general 'recovery' field will be used for both cash "
          "and physical settlement");
    FIELD_MAKE_OPTIONAL(cashRecovery);    

    FIELD(physicalSettlement, "Transient: Physical settlement");
    FIELD_MAKE_TRANSIENT(physicalSettlement);
    FIELD(cashSettlement, "Transient: Cash settlement");
    FIELD_MAKE_TRANSIENT(cashSettlement);
}


CClassConstSP const TranchePartCashSettlementOverride::TYPE = 
    CClass::registerClassLoadMethod("TranchePartCashSettlementOverride", 
                                    typeid(TranchePartCashSettlementOverride), 
                                    load);

/** Included in ProductsLib to force the linker to include this file */
bool TranchePartCashSettlementOverrideLoad() {
    return (TranchePartCashSettlementOverride::TYPE != 0);
}

DRLIB_END_NAMESPACE
