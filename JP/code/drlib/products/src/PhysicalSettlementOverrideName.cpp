//----------------------------------------------------------------------------
//
//   Group       : CH Quantitative Research
//
//   Filename    : PhysicalSettlementOverrideName.cpp
//
//   Description : Class used by credit instruments (e.g., CIS) in order 
//                 to override a name's default-related parameters, when
//                 the default is settled physically.
//
//   Author      : Jose Hilera
//
//   Date        : February 2006
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/Format.hpp"
#include "edginc/PhysicalSettlementOverrideName.hpp"
#include "edginc/IBadDayAdjuster.hpp"
#include "edginc/DeliveryDetails.hpp"
#include "edginc/ICreditContingentLeg.hpp"
#include "edginc/CtgLegLossPerDefault.hpp"
#include "edginc/FeeLegReductionPerDefault.hpp"

DRLIB_BEGIN_NAMESPACE


PhysicalSettlementOverrideName::PhysicalSettlementOverrideName() : 
    DetailedCreditEventOverrideName(TYPE),
    totalExpectedNotional(1.0)
{}

/** Public constructor, with all fields passed in */
PhysicalSettlementOverrideName::PhysicalSettlementOverrideName(
        ICDSParSpreadsWrapper name,
        DateTime eventDeterminationDate,
        CIntSP triggerDelay,
        CDoubleSP recovery,
        NoticeOfPhysicalSettlementSP NoPS,
        CIntSP defaultToSettlementDelay,
        DeliveryDetailsArraySP deliveries,
        double totalExpectedNotional) :
    DetailedCreditEventOverrideName(TYPE,
                                    name, 
                                    eventDeterminationDate,
                                    triggerDelay,
                                    recovery),
    NoPS(NoPS),
    defaultToSettlementDelay(defaultToSettlementDelay),
    deliveries(deliveries),
    totalExpectedNotional(totalExpectedNotional)
{}


PhysicalSettlementOverrideName::~PhysicalSettlementOverrideName()
{}


/** Called immediately after object constructed */
void PhysicalSettlementOverrideName::validatePop2Object() {
    static const string method = 
        "PhysicalSettlementOverrideName::validatePop2Object";

    // Validate fields in the parent class
    DetailedCreditEventOverrideName::validatePop2Object();

    // Validate fields in this class
    double nopsNotional;
    double deliveryNotional = getNotionalScheduledForDelivery();
    if (!NoPS) {
        nopsNotional = 0.0;
        if (deliveryNotional != 0.0) {
            throw ModelException(method,
                                 "For defaulted name " + name.getName() + ", "
                                 "the total notional scheduled for delivery is " +
                                 Format::toString(100.0 * deliveryNotional) +
                                 "%, but there is no NoPS specified.");
        }
    }
    else {
        nopsNotional = NoPS->getNoPSNotionalFraction();
    }

    if (nopsNotional > totalExpectedNotional) {
        // This is to take into account part cash settlements, where
        // only up to x% of the notional is expected to settle physically
        // (the rest being settled in cash)
        throw ModelException(method,
                             "For defaulted name " + name.getName() +
                             ", the NoPS specifies that " +
                             Format::toString(100.0 * nopsNotional) +
                             "% will be triggered, but no more than " +
                             Format::toString(100.0 * totalExpectedNotional) +
                             "% is expected to settle physically.");
    }

    if (nopsNotional < deliveryNotional) {
        throw ModelException(method,
                             "For defaulted name " + name.getName() +
                             ", the total notional scheduled for delivery is " +
                             Format::toString(100.0 * deliveryNotional) +
                             "%, whereas the NoPS specifies that only " +
                             Format::toString(100.0 * nopsNotional) +
                             "% will be delivered.");
    }
    else if (nopsNotional > deliveryNotional) {
        if (!defaultToSettlementDelay) {
            throw ModelException(method,
                                 "For defaulted name " + name.getName() +
                                 ", since not all deliveries are specified "
                                 "(expect to receive " +
                                 Format::toString(100.0 * nopsNotional) +
                                 "% of the notional, but only " +
                                 Format::toString(100.0 * deliveryNotional) +
                                 "% is scheduled), parameter "
                                 "defaultToSettlementDelay must be specified.");
        }
        else {
            // Will use the defaultToSettlementDelay, so verify it is right
            if (defaultToSettlementDelay->intValue() < 0) {
                throw ModelException(method,
                                     "The defaultToSettlementDelay for name " +
                                     name.getName() +
                                     " is negative, and negative delays are "
                                     "not accepted.");
            }
        }
    }

    // Cross-validate fields in the parent class and in this class
    if (!!NoPS && eventDeterminationDate.empty()) {
        throw ModelException(method,
                             "Name " + name.getName() +
                             " has not been triggered (eventDeterminationDate "
                             "is not set) so there can be no NoPS for it. "
                             "Either mark the name as triggered or remove "
                             "its NoPS.");
    }

    if (!!deliveries) {
        if (eventDeterminationDate.empty()) {
            throw ModelException(method,
                                 "Name " + name.getName() +
                                 " has not been triggered (eventDeterminationDate "
                                 "is not set) so there can be no default "
                                 "deliveries scheduled for it. "
                                 "Either mark the name as triggered or remove "
                                 "its deliveries.");
        }
        else {
            int numOfDeliveries = deliveries->size();
            // Check that all deliveries are in chronological order
            for (int i=1; i < numOfDeliveries; ++i) {
                if ((*deliveries)[i-1]->deliveryDate > (*deliveries)[i]->deliveryDate) {
                    throw ModelException(method,
                                         name.getName() +
                                         "'s delivery dates are not in "
                                         "chronological order (delivery[" +
                                         Format::toString(i-1) + "] = " +
                                         (*deliveries)[i-1]->deliveryDate.toString() +
                                         ", delivery[" + 
                                         Format::toString(i) + "] = " +
                                         (*deliveries)[i]->deliveryDate.toString() +
                                         ")");
                }
            }

            // Check if the first delivery happens after the trigger date
            if ((numOfDeliveries > 0) && 
                (*deliveries)[0]->deliveryDate < eventDeterminationDate) 
            {
                throw ModelException(method,
                                     "Name " + name.getName() +
                                     " has been triggered on " +
                                     eventDeterminationDate.toString() +
                                     " so cannot have deliveries scheduled "
                                     "before that (the first delivery is on " +
                                     (*deliveries)[0]->deliveryDate.toString() +
                                     ")");
            }
        }
    }
}


/** Returns the overall notional fraction scheduled for delivery */
double PhysicalSettlementOverrideName::getNotionalScheduledForDelivery() const {
    double deliveryNotional = 0.0;
    if (!!deliveries) {
        for (int i=0; i < deliveries->size(); ++i) {
            deliveryNotional += (*deliveries)[i]->getNotionalFractionInDelivery();
        }
    }
    return deliveryNotional;
}


/* Returns the fraction of the notional which is expected to be settled
 * but still has no settlement (delivery) date */
double PhysicalSettlementOverrideName::getUnscheduledNotional() const {
    double unscheduledNotional;
    if (!NoPS) {
        // If there is no NoPS, assume all the notional will be delivered
        unscheduledNotional = totalExpectedNotional;
    }
    else {
        unscheduledNotional = NoPS->getNoPSNotionalFraction();
        // Substract the notional scheduled for delivery, since it has 
        // already been taken into consideration
        unscheduledNotional -= getNotionalScheduledForDelivery();
    }
    return unscheduledNotional;
}


/** Returns the portfolio reductions for the fee leg (corresponding to losses
    and recovered amounts). */
CtgLegLossPerDefaultArraySP 
PhysicalSettlementOverrideName::historicContingentLegLosses(
    const double nameNotional,
    const double recoveryRate,
    const DateTime& lastTriggerDate,
    const DateTime& creditEventDate,
    IBadDayAdjusterConstSP bda) const
{
    static const string method(
        "PhysicalSettlementOverrideName::getContingentLegCashFlows");

    const DateTime& determinationDate = getEventDeterminationDate(
        creditEventDate, lastTriggerDate, bda);

    if (determinationDate.empty() || (determinationDate > lastTriggerDate)) {
        // The default cannot be triggered, or it was triggered after 
        // the lastTriggerDate -> no contingent payment will be made
        return CtgLegLossPerDefaultArraySP();
    }

    // Check if we are overriding the recovery rate here
    double recRate = getOverridenRecoveryRate(recoveryRate);

    // Initialise an empty array of cash flows - and reserve enough space for
    // deliveries->size() + 1 cashflows (the maximum possible number of entries)
    CashFlowArraySP contingentCashFlows(new CashFlowArray(0));
    int numOfDeliveries = (!deliveries) ? 0 : deliveries->size();
    contingentCashFlows->reserve(numOfDeliveries + 1);

    // Need to estimate when the notional not scheduled for delivery will
    // arrive (if there is any).
    // We need to sort the cashflows in the array so only create the 
    // cashflow here - it will be inserted into the array later
    double unscheduledNotional = getUnscheduledNotional();
    bool insertUnscheduledCashFlow = false;
    CashFlowSP unscheduledCashFlow;  // will be added if insertUnscheduledCashFlow
    if (unscheduledNotional > 0) {
        const DateTime& settlementDate = 
            estimateSettlementDate(creditEventDate, lastTriggerDate, valueDate, bda);
        
        if (settlementDate <= valueDate) {
            // This cashflow should be ignored
        }
        else {
            // Create the new cashflow
            unscheduledCashFlow.reset(new CashFlow (
                settlementDate, 
                nameNotional * unscheduledNotional * (1.0 - recRate)));

            insertUnscheduledCashFlow = true;
        }
    }

    // Insert the delivery and unscheduled cashflows in the array
    if (!!deliveries) {
        for (int i=0; i < numOfDeliveries; ++i) {
            const DateTime& deliveryDate = (*deliveries)[i]->deliveryDate;

            if (deliveryDate <= valueDate) {
                // This cashflow should be ignored
            }
            else {
                // Get the remaining data for this delivery
                const double deliveryNotional = nameNotional * 
                    (*deliveries)[i]->getNotionalFractionInDelivery();

                const double deliveryRecovery = !((*deliveries)[i]->recovery) ? 
                    recRate :
                    (*deliveries)[i]->recovery->doubleValue();

                // If we need to insert the unscheduled cashflow and it happens 
                // before this delivery, insert it first
                if (insertUnscheduledCashFlow && 
                    (unscheduledCashFlow->date < deliveryDate))
                {
                    contingentCashFlows->push_back(*unscheduledCashFlow);
                    insertUnscheduledCashFlow = false; // already inserted
                }

                // Create the new cashflow and insert it into contingentCashFlows
                CashFlow newCashFlow(deliveryDate,
                                     deliveryNotional * (1.0 - deliveryRecovery));

                contingentCashFlows->push_back(newCashFlow);
            }
        }
    }

    // If the unscheduled cashflow has not been inserted yet, do it now
    if (insertUnscheduledCashFlow) {
        contingentCashFlows->push_back(*unscheduledCashFlow);
        // insertUnscheduledCashFlow = false; // Not required here
    }

    CtgLegLossPerDefaultArraySP ctgReductions =
        CtgLegLossPerDefault::produceLossesPerDefault(creditEventDate,
                                                      contingentCashFlows);
    
    return ctgReductions;
}


/** Returns the overall loss resulting from the settlement of this name */
double PhysicalSettlementOverrideName::getOverallLoss(
    const double nameNotional,
    const double recoveryRate) const 
{
    const double unscheduledNotional = getUnscheduledNotional();

    const double totalExpectedNotionalFraction = (!NoPS) ?
        totalExpectedNotional : // no NoPS so assume all the notional will be delivered
        NoPS->getNoPSNotionalFraction();

    // Check if we are overriding the recovery rate here
    double recRate = getOverridenRecoveryRate(recoveryRate);
                                          
    // Initialise totalLoss with the unscheduled notional
    double totalLoss = unscheduledNotional * (1.0 - recRate);

    const int numDeliveries = (!deliveries) ? 0 : deliveries->size();
    double deliveredProportionSoFar = 0.0;
    for (int i=0; i < numDeliveries; ++i) {
        const double deliveryNotionalFraction = Maths::min(
            (*deliveries)[i]->getNotionalFractionInDelivery(),
            totalExpectedNotionalFraction - deliveredProportionSoFar);

        deliveredProportionSoFar += deliveryNotionalFraction;

        if (!Maths::isZero(deliveryNotionalFraction)) {
            const double deliveryRecoveryRate =
                (*deliveries)[i]->getDeliveryRecoveryRate(recRate);
                
            totalLoss += deliveryNotionalFraction * 
                nameNotional * (1.0 - deliveryRecoveryRate);
        }
    }
    return totalLoss;
}


/* Returns the overall recovery rate accross all settlements, such that
 * total loss = overall defaulted notional * (1-overall recovery rate) */
double PhysicalSettlementOverrideName::getOverallRecoveryRate(
    double recoveryRate) const
{
    static const string method(
        "PhysicalSettlementOverrideName::getOverallRecoveryRate");
    double overallNotional = getOverallDefaultedNotionalFraction();
    if (!Maths::isZero(overallNotional)) {
        double overallLoss = getOverallLoss(1.0, // notional: dealing with fractions here
                                            recoveryRate);

        return (1.0 - overallLoss / overallNotional);
    }
    // else the name will not cause any losses, so it's all recovered
    return 1.0;
}


/** Returns the percentage of the notional that has been triggered 
 * for protection */
double PhysicalSettlementOverrideName::getOverallDefaultedNotionalFraction() const
{
    const double notionalFraction = (!NoPS) ? 
        totalExpectedNotional : NoPS->getNoPSNotionalFraction();

    return notionalFraction;
}


/** Returns the accrual pay date. For physical settlements this is the date of
 * the first delivery */
const DateTime PhysicalSettlementOverrideName::getAccrualPayDate(
    const DateTime& creditEventDate,
    const DateTime& lastTriggerDate,
    const DateTime& valueDate,
    IBadDayAdjusterConstSP bda) const
{
    DateTime accPayDate;
    
    const int numOfDeliveries = (!deliveries) ? 0 : deliveries->size();
    if (!NoPS || (!deliveries) || (numOfDeliveries == 0)) {
        // accrual will be paid on settlementDate
        accPayDate = estimateSettlementDate(
            creditEventDate, lastTriggerDate, valueDate, bda);
    }
    else {
        // Get the date of the first scheduled delivery - since the deliveries
        // are in chronological order, can ignore the other deliveries
        accPayDate = (*deliveries)[0]->deliveryDate;

        // If there is still notional with no delivery schedule, consider the 
        // date when it will be delivered
        double unscheduledNotional = getUnscheduledNotional();
        if (unscheduledNotional > 0) {
            const DateTime& settlementDate = estimateSettlementDate(
                creditEventDate, lastTriggerDate, valueDate, bda);
            accPayDate = accPayDate.min(settlementDate);
        }
    }
    return accPayDate;
}


/** Returns the estimated settlement date based on the defaultToSettlementDelay */
const DateTime PhysicalSettlementOverrideName::estimateSettlementDate(
    const DateTime& creditEventDate,
    const DateTime& lastTriggerDate,
    const DateTime& valueDate,
    IBadDayAdjusterConstSP bda) const
{
    static const string method = 
        "PhysicalSettlementOverrideName::estimateSettlementDate";

    DateTime settlementDate;
    const DateTime& determinationDate = getEventDeterminationDate(
        creditEventDate, lastTriggerDate, bda);

    if (determinationDate.empty()) {
        // An empty determinationDate means the default has not been trigered.
        // Therefore return an empty settlement date too.
        // Note settlementDate is an empty date already.
    }
    else {

#ifdef QLIB_REMOVE_EDD_VALIDATION_CDS
        // Need to verify that the determinationDate
        // falls after the creditEventDate
        if (creditEventDate > determinationDate) {
            throw ModelException(method,
                                 name.getName() + 
                                 "'s event determination date (" +
                                 determinationDate.toString() +
                                 ") can not be before the default date ("+
                                 creditEventDate.toString() + 
                                 ") because not all deliveries have "
                                 "been specified.");
        }
#endif
        if (!defaultToSettlementDelay) {
            throw ModelException(method,
                                 "For defaulted name " + name.getName() +
                                 ", since not all deliveries are specified, "
                                 "parameter defaultToSettlementDelay must "
                                 "be present.");
        }

        settlementDate = bda->badDayAdjust(
            creditEventDate.rollDate(defaultToSettlementDelay->intValue()));
        
        // If this date is before valueDate or determinationDate (whichever is 
        // greater) it means that the estimate for defaultToSettlementDelay was
        // short. Will use that date
        const DateTime& adjValueDate = bda->badDayAdjust(valueDate.rollDate(1)); // JLH. make this "1" configurable
        const DateTime& maxDate = determinationDate.max(adjValueDate);
        if (settlementDate < maxDate) {
            settlementDate = maxDate;
        }
    }
    return settlementDate;
}


IObject* PhysicalSettlementOverrideName::defaultPhysicalSettlementOverrideName() {
    return new PhysicalSettlementOverrideName();
}


void PhysicalSettlementOverrideName::load (CClassSP& clazz) {
    clazz->setPublic(); // make visible to EAS/spreadsheet
    REGISTER(PhysicalSettlementOverrideName, clazz);
    SUPERCLASS(DetailedCreditEventOverrideName);
    EMPTY_SHELL_METHOD(defaultPhysicalSettlementOverrideName);

    FIELD(NoPS,                        "Specification of the amounts that will be "
                                       "physically settled");
    FIELD(defaultToSettlementDelay,    "Delay between credit event and settlement, "
                                       "in calendar days");
    FIELD(deliveries,                  "Details for the scheduled physical "
                                       "deliveries - in chronological order");
    FIELD(totalExpectedNotional,"Transient field: Fraction of the notional "
                                       "expected to settle physically");

    FIELD_MAKE_OPTIONAL(NoPS);
    FIELD_MAKE_OPTIONAL(defaultToSettlementDelay);
    FIELD_MAKE_OPTIONAL(deliveries);

    FIELD_MAKE_TRANSIENT(totalExpectedNotional);
}


CClassConstSP const PhysicalSettlementOverrideName::TYPE = 
    CClass::registerClassLoadMethod("PhysicalSettlementOverrideName", 
                                    typeid(PhysicalSettlementOverrideName), 
                                    load);

// Array has to have its own type
DEFINE_TEMPLATE_TYPE(PhysicalSettlementOverrideNameArray);


/** Included in ProductsLib-modified::linkInClasses() via the productSrcsMap
 * script to force the linker to include this file */
bool PhysicalSettlementOverrideNameLoad() {
    return (PhysicalSettlementOverrideName::TYPE != 0);
}

DRLIB_END_NAMESPACE
