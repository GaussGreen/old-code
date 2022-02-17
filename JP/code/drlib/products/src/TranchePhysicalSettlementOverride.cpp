//----------------------------------------------------------------------------
//
//   Group       : Quatitative Research
//
//   Filename    : TranchePhysicalSettlementOverride.cpp
//
//   Description : Class used by CDO in order to override a name's default-
//                 related parameters, when the default is settled physically.
//
//   Author      : Jose Hilera
//
//   Date        : May 2006
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/Format.hpp"
#include "edginc/CreditFeeLeg.hpp"
#include "edginc/IBadDayAdjuster.hpp"
#include "edginc/IRebateCalculator.hpp"
#include "edginc/DeliveryDetails.hpp"
#include "edginc/FeeLegReductionPerDefault.hpp"
#include "edginc/CtgLegLossPerDefault.hpp"
#include "edginc/TranchePhysicalSettlementOverride.hpp"

DRLIB_BEGIN_NAMESPACE

const bool DEFAULT_SETTLE_DELIVERIES_SEPARATELY(true);
const bool DEFAULT_MORE_DELIVERIES_PENDING(true);

TranchePhysicalSettlementOverride::TranchePhysicalSettlementOverride() : 
    TrancheCreditEventOverride(TYPE), 
    settleDeliveriesSeparately(DEFAULT_SETTLE_DELIVERIES_SEPARATELY),
    moreDeliveriesPending(DEFAULT_MORE_DELIVERIES_PENDING),
    notionalFractionSettlingPhysically(1.0)
{}


/** Public constructor, with all fields passed in */
TranchePhysicalSettlementOverride::TranchePhysicalSettlementOverride(
    CIntSP triggerDelay,
    DateTime eventDeterminationDate,
    CDoubleSP recovery,
    DateTime valueDate,
    NoticeOfPhysicalSettlementSP NoPS,
    DeliveryDetailsArraySP deliveries,
    CBoolSP settleDeliveriesSeparatelySP,
    CIntSP defaultToCalculationDelay,
    const DateTime& calculationDate,
    const DateTime& cutoffDate,
    CBoolSP moreDeliveriesPendingSP,
    double notionalFractionSettlingPhysically) :
        TrancheCreditEventOverride(TYPE,
                                   triggerDelay, 
                                   eventDeterminationDate,
                                   recovery,
                                   valueDate),
        NoPS(NoPS),
        deliveries(deliveries),
        defaultToCalculationDelay(defaultToCalculationDelay),
        calculationDate(calculationDate),
        cutoffDate(cutoffDate),
        notionalFractionSettlingPhysically(notionalFractionSettlingPhysically)
{
    // Initialise values using the smart pointers. If null, set the default
    settleDeliveriesSeparately = !settleDeliveriesSeparatelySP ? 
        DEFAULT_SETTLE_DELIVERIES_SEPARATELY :
        settleDeliveriesSeparatelySP->boolValue();

    moreDeliveriesPending = !moreDeliveriesPendingSP ?
        DEFAULT_MORE_DELIVERIES_PENDING :
        moreDeliveriesPendingSP->boolValue();
}


TranchePhysicalSettlementOverride::~TranchePhysicalSettlementOverride()
{}


/** Called immediately after object constructed */
void TranchePhysicalSettlementOverride::validatePop2Object() {
    static const string method = 
        "TranchePhysicalSettlementOverride::validatePop2Object";

    // Validate fields in the parent class
    TrancheCreditEventOverride::validatePop2Object();
    
    if (settleDeliveriesSeparately) {
        if (!calculationDate.empty()) {
            throw ModelException(method,
                                 "Since deliveries are settling separately "
                                 "the overall calculation date cannot be set. "
                                 "Either specify that calculation date for the "
                                 "appropriate delivery, or mark the name as "
                                 "settling in a single date.");
        }
    }
    else if (!cutoffDate.empty() &&        // If both dates are present, fail
             !calculationDate.empty() &&   // if cutoff date > calculation date
             cutoffDate > calculationDate) 
    {
            throw ModelException(method,
                                 "Deliveries are expected to arrive until "
                                 "cut-off date (" +
                                 cutoffDate.toString() +
                                 ") so the calculation date must fall after "
                                 "that (currently it is " + 
                                 calculationDate.toString() + ").");
    }

    if (!defaultToCalculationDelay) {
        // We cannot tell at this point all cases where 
        // defaultToCalculationDelay will be required (maybe all details
        // are not yet known) - Will have to perform other checks where
        // required also
        if (!deliveries) {
            throw ModelException(method,
                                 "If no deliveries are specified the "
                                 "defaultToCalculationDelay must be provided");
        }

        if (!settleDeliveriesSeparately && calculationDate.empty()) {
            throw ModelException(method,
                                 "Since all deliveries are settling together "
                                 "and the calculation date is not yet known, the "
                                 "defaultToCalculationDelay must be provided");
        }
    }
    else {
        // May use the defaultToCalculationDelay, so verify it is right
        if (defaultToCalculationDelay->intValue() < 0) {
            throw ModelException(method,
                                 "defaultToCalculationDelay is " + 
                                 Format::toString(defaultToCalculationDelay->intValue()) +
                                 ", and negative delays are not accepted.");
        }
    }


    // Cross validate NoPS and Deliveries
    double deliveryNotional = getNotionalScheduledForDelivery();
    if (!NoPS) {
        if(!Maths::isZero(deliveryNotional)) {
            throw ModelException(method,
                                 "The total notional scheduled for delivery is " +
                                 Format::toString(100.0 * deliveryNotional) +
                                 "%, but there is no NoPS specified.");
        }
    }
    else {
        double nopsNotional = NoPS->getNoPSNotionalFraction();
        if (nopsNotional > notionalFractionSettlingPhysically) {
            // This is to take into account part cash settlements, where
            // only up to x% of the notional is expected to settle physically
            // (the rest being settled in cash)
            throw ModelException(method,
                                 "The NoPS specifies that " +
                                 Format::toString(100.0 * nopsNotional) +
                                 "% will be triggered, but no more than " +
                                 Format::toString(100.0 * notionalFractionSettlingPhysically) +
                                 "% is expected to settle physically.");
        }
    }

    double totalExpectedNotionalFraction = getTotalExpectedNotional();
    if (totalExpectedNotionalFraction < deliveryNotional) {
        // This is fine - the notional delivered on top of 
        // totalExpectedNotionalFraction will be considered a "gif" to the 
        // protection seller, and will potentially be treated as recovered 
        // notional
        ;
    }

    if (cutoffDate.empty()) {
        double acceptedDeliveryNotional = Maths::min(totalExpectedNotionalFraction, 
                                                     deliveryNotional);
        double unscheduledNotionalFraction = getUnscheduledNotionalFraction();

        if (acceptedDeliveryNotional + unscheduledNotionalFraction <
            notionalFractionSettlingPhysically)
        {
            throw ModelException(method,
                                 Format::toString(notionalFractionSettlingPhysically -
                                                  acceptedDeliveryNotional -
                                                  unscheduledNotionalFraction) +
                                 "% of the notional can not be delivered. "
                                 "The cut-off date needs to be provided for "
                                 "rebate computation purposes (and, potentially, "
                                 "to determine the point at which notional will "
                                 "be considered recovered).");
        }
    }
    else { // cutoffDate present
        if (eventDeterminationDate.empty()) {
            throw ModelException(method,
                                 "The cut-off date can not be present if the "
                                 "default has not been triggered yet ("
                                 "eventDeterminationDate is empty)");            
        }
        if (cutoffDate < eventDeterminationDate) { // eventDeterminationDate not empty
            throw ModelException(method,
                                 "The cut-off date (" + cutoffDate.toString() +
                                 ") must fall after the eventDeterminationDate (" +
                                 eventDeterminationDate.toString() + ").");
        }
    }        

    // Cross-validate fields in the parent class and in this class
    if (!!NoPS && eventDeterminationDate.empty()) {
        throw ModelException(method,
                             "The name has not been triggered "
                             "(eventDeterminationDate is not set) so there can "
                             "be no NoPS for it. Either mark the name as "
                             "triggered or remove its NoPS.");
    }

    if (!!deliveries) {
        if (eventDeterminationDate.empty()) {
            throw ModelException(method,
                                 "The name has not been triggered "
                                 "(eventDeterminationDate is not set) so there "
                                 "can be no default deliveries scheduled for "
                                 "it. Either mark the name as triggered or "
                                 "remove its deliveries.");
        }
        else {
            int numOfDeliveries = deliveries->size();
            // Check that all deliveries are in chronological order
            for (int i=1; i < numOfDeliveries; ++i) {
                if ((*deliveries)[i-1]->getDeliveryDate() > (*deliveries)[i]->getDeliveryDate()) {
                    throw ModelException(method,
                                         "Delivery dates are not in "
                                         "chronological order (delivery[" +
                                         Format::toString(i-1) + "] = " +
                                         (*deliveries)[i-1]->getDeliveryDate().toString() +
                                         ", delivery[" + 
                                         Format::toString(i) + "] = " +
                                         (*deliveries)[i]->getDeliveryDate().toString() +
                                         ")");
                }
            }

            // Check if the first delivery happens after the trigger date
            if ((numOfDeliveries > 0) && 
                (*deliveries)[0]->getDeliveryDate() < eventDeterminationDate) 
            {
                throw ModelException(method,
                                     "The name has been triggered on " +
                                     eventDeterminationDate.toString() +
                                     " so cannot have deliveries scheduled "
                                     "before that (the first delivery is on " +
                                     (*deliveries)[0]->getDeliveryDate().toString() +
                                     ")");
            }

            // Check if the last delivery happens before the cutoffDate
            if (!cutoffDate.empty() && 
                (numOfDeliveries > 0) && 
                (deliveries->back()->getDeliveryDate() > cutoffDate)) 
            {
                throw ModelException(method,
                                     "All deliveries must be scheduled before "
                                     "the cut-off date (" +
                                     cutoffDate.toString() +
                                     ") but there is a delivery on " +
                                     (*deliveries)[numOfDeliveries-1]->getDeliveryDate().toString());
            }
        }
    }

    return;
}


/** Returns the contingent leg cashflows corresponding to a defaulted name.
 * This object contains the credit event parameters */
CtgLegLossPerDefaultArraySP 
TranchePhysicalSettlementOverride::historicContingentLegLosses(
    const double nameNotional,
    const double recoveryRate,
    const DateTime& lastTriggerDate,
    const DateTime& creditEventDate,
    IBadDayAdjusterConstSP bda) const
{
    static const string method = 
        "TranchePhysicalSettlementOverride::historicContingentLegLosses";

    const DateTime& determinationDate = getEventDeterminationDate(
        creditEventDate, lastTriggerDate, bda);

    if (determinationDate.empty()) {
        // The default has not been triggered, or it was triggered after 
        // the lastTriggerDate -> no losses
        return CtgLegLossPerDefaultArraySP();
    }
    
    // Initialise an empty array of cash flows
    CashFlowArraySP contingentLosses;

    if (settleDeliveriesSeparately) {
        double totalExpectedNotionalFraction = getTotalExpectedNotional();

        // Check if we are overriding the recovery rate here
        double recRate = (!recovery) ?
            recoveryRate : // No override, use argument passed in to this method
            recovery->doubleValue();

        // Insert the deliveries' losses in the array
        const int numDeliveries = (!deliveries) ? 0 : deliveries->size();
        double deliveredProportionSoFar = 0.0;
        for (int i=0; i < numDeliveries; ++i) {
            double deliveryNotionalFraction = Maths::min(
                (*deliveries)[i]->getNotionalFractionInDelivery(),
                totalExpectedNotionalFraction - deliveredProportionSoFar);

            deliveredProportionSoFar += deliveryNotionalFraction;

            if (!Maths::isZero(deliveryNotionalFraction)) {
                const double deliveryNotional = 
                    nameNotional * deliveryNotionalFraction;

                const double deliveryRecoveryRate =
                    (*deliveries)[i]->getDeliveryRecoveryRate(recRate);

                const DateTime& deliveryCalcDate = 
                    (*deliveries)[i]->deliveryCalcDate(creditEventDate,
                                                       determinationDate,
                                                       valueDate,
                                                       bda);
                
                // Create the new cashflow and insert it into contingentLosses
                CashFlowArraySP deliveryCashFlow(new CashFlowArray(
                    1, CashFlow(deliveryCalcDate,
                                deliveryNotional * (1.0 - deliveryRecoveryRate))));

                contingentLosses = CashFlow::merge(contingentLosses, 
                                                   deliveryCashFlow);
            }
        }

        // Need to estimate when the notional not scheduled for delivery will
        // arrive (if there is any left). Note getUnscheduledNotionalFraction
        // takes into account that the "delivered notional" cannot be greater 
        // than the "nops Notional"
        double unscheduledNotionalFraction = getUnscheduledNotionalFraction();
        if (unscheduledNotionalFraction > 0) {
            const DateTime& unscheduledCalcDate = 
                commonCalculationDate(creditEventDate, lastTriggerDate, bda);
        
            const double amount = nameNotional * 
                unscheduledNotionalFraction * (1.0 - recRate);

            // Set and insert the unscheduled cash flow
            CashFlowArraySP unscheduledCashFlow(new CashFlowArray(
                1, CashFlow(unscheduledCalcDate, amount)));
            contingentLosses = CashFlow::merge(contingentLosses,
                                               unscheduledCashFlow);
        }
    }
    else {
        // Since all deliveries settle at once, need to determine:
        // - the (one) calculation date
        // - the losses due to each delivery (each can have different notional
        //   and recovery rate)
        // - the losses due to the unscheduled notional
        const DateTime& calcDate = 
            commonCalculationDate(creditEventDate, lastTriggerDate, bda);
       const double totalLoss = getOverallLoss(nameNotional, recoveryRate);
       
       // Produce the actual loss cashflow
       contingentLosses.reset(new CashFlowArray(1, CashFlow(calcDate, 
                                                            totalLoss)));
    }

    CtgLegLossPerDefaultArraySP ctgLegLoss = 
        CtgLegLossPerDefault::produceLossesPerDefault(creditEventDate,
                                                      contingentLosses);

    return ctgLegLoss;
}


/** Estimates the calculation date for the outstanding notional with no
 * delivery information */
DateTime TranchePhysicalSettlementOverride::commonCalculationDate(
    const DateTime& creditEventDate,
    const DateTime& lastTriggerDate,
    IBadDayAdjusterConstSP bda) const
{
    static const string method(
        "TranchePhysicalSettlementOverride::commonCalculationDate");
    
    DateTime defaultCalculationDate;
    if (calculationDate.empty()) {
        const DateTime& determinationDate = getEventDeterminationDate(
            creditEventDate, lastTriggerDate, bda);

        if (determinationDate.empty()) {
            // An empty determinationDate means the default has not been trigered.
            // Therefore return an empty calculation date date too.
            // Note defaultCalculationDate is an empty date already.
        }
        else {
#ifdef QLIB_REMOVE_EDD_VALIDATION_CDO
            // CalculationDate is empty, so need to check that the determinationDate
            // falls after the creditEventDate
            if (creditEventDate > determinationDate) {
                throw ModelException(method,
                                     "Event determination date (" +
                                     determinationDate.toString() +
                                     ") can not be before the default date ("+
                                     creditEventDate.toString() + 
                                     ") because the calculation date has not "
                                     "been specified.");
            }
#endif

            defaultCalculationDate = 
                ITrancheCreditEventOverride::rollAndAdjustDate(
                    creditEventDate,
                    defaultToCalculationDelay,
                    valueDate,
                    determinationDate,
                    bda);
        }
    }
    else {
        // We know what the calculation date is - return it
        defaultCalculationDate = calculationDate;
    }

    return defaultCalculationDate;
}



/** Returns the notional specified in the NoPS, or 
 * notionalFractionSettlingPhysically if no NoPS have been provided */
double TranchePhysicalSettlementOverride::getTotalExpectedNotional() const {
    return (!NoPS) ? 
        notionalFractionSettlingPhysically : NoPS->getNoPSNotionalFraction();
}


/** Returns the overall notional fraction scheduled for delivery. Note this
 * may be greater than the notional indicated in the NoPS  */
double TranchePhysicalSettlementOverride::getNotionalScheduledForDelivery() const {
    double deliveryNotional = 0.0;
    if (!!deliveries) {
        for (int i=0; i < deliveries->size(); ++i) {
            deliveryNotional += 
                (*deliveries)[i]->getNotionalFractionInDelivery();
        }
    }
    return deliveryNotional;
}


/* Returns the fraction of the notional which is expected to settle
 * but still has no settlement (delivery) date. */
double TranchePhysicalSettlementOverride::getUnscheduledNotionalFraction() const 
{
    double unscheduledNotional;
    if (!NoPS) {
        // If there is no NoPS, assume all the notional will be delivered
        unscheduledNotional = notionalFractionSettlingPhysically;
    }
    else if (!moreDeliveriesPending ||      // Nothing else pending
             (!cutoffDate.empty() && (valueDate > cutoffDate)) || // Gone past cutoff date
             (!calculationDate.empty() && (valueDate > calculationDate))) // past calcDate
    {
        unscheduledNotional = 0.0;
    }
    else {
        double totalExpectedNotionalFraction = getTotalExpectedNotional();
        double deliveryNotional = getNotionalScheduledForDelivery();
        // Substract the notional scheduled for delivery, since it has 
        // already been taken into consideration
        unscheduledNotional = 
            Maths::max(totalExpectedNotionalFraction - deliveryNotional, 
                       0.0); 
    }
    return unscheduledNotional;
}



/** Returns the portfolio reductions for the fee leg (corresponding to losses
 * and recovered amounts). This object contains the credit event parameters. */
FeeLegReductionPerDefaultArraySP 
TranchePhysicalSettlementOverride::historicFeeLegReductions(
    const double nameNotional,
    const double recoveryRate,
    const DateTime& lastTriggerDate,
    const DateTime& creditEventDate,
    IBadDayAdjusterConstSP bda) const
{
    // Initialize the output losses and recovered notional arrays
    FeeLegReductionPerDefaultArraySP feeReductions;

    /* Obtain the settlement parameters using the information in
     * this override */
    const DateTime& determinationDate = getEventDeterminationDate(
        creditEventDate, lastTriggerDate, bda);

    if (determinationDate.empty()) {
        // The default can not be triggered, or it was triggered after 
        // the lastTriggerDate -> no losses or recovered amounts
        return feeReductions;
    }     

    if (settleDeliveriesSeparately) {
        // Check if we are overriding the recovery rate here
        double recRate = (!recovery) ?
            recoveryRate : // No override, use argument passed in to this method
            recovery->doubleValue();

        feeReductions.reset(new FeeLegReductionPerDefaultArray(0));

        double totalExpectedNotionalFraction = getTotalExpectedNotional();

        // Insert each delivery's losses
        const int numDeliveries = (!deliveries) ? 0 : deliveries->size();
        double deliveredProportionSoFar = 0.0;
        for (int i=0; i < numDeliveries; ++i) {
            double deliveryNotionalFraction = Maths::min(
                (*deliveries)[i]->getNotionalFractionInDelivery(),
                totalExpectedNotionalFraction - deliveredProportionSoFar);

            deliveredProportionSoFar += deliveryNotionalFraction;

            if (!Maths::isZero(deliveryNotionalFraction)) {
                const double deliveryNotional = 
                    nameNotional * deliveryNotionalFraction;

                const double deliveryRecoveryRate =
                    (*deliveries)[i]->getDeliveryRecoveryRate(recRate);

                const DateTime& deliveryCalcDate = 
                    (*deliveries)[i]->deliveryCalcDate(creditEventDate,
                                                       determinationDate,
                                                       valueDate,
                                                       bda);

                // Add to the total feeReductions this delivery's reductions
                feeReductions->push_back(FeeLegReductionPerDefaultSP(
                    new FeeLegReductionPerDefault(determinationDate,
                                                  deliveryCalcDate, // effective date
                                                  deliveryCalcDate, // calc date
                                                  deliveryNotional, // defaulted notional
                                                  deliveryNotional, // total notional
                                                  deliveryRecoveryRate)));
            }
        }

        // Need to estimate when the notional not scheduled for delivery will
        // arrive (if there is any left). Note getUnscheduledNotionalFraction
        // takes into account that the "delivered notional" cannot be greater 
        // than the "nops Notional"
        const double unscheduledNotionalFraction = 
            getUnscheduledNotionalFraction();
        if (unscheduledNotionalFraction > 0) {
            const DateTime& unscheduledCalcDate = 
                commonCalculationDate(creditEventDate, lastTriggerDate, bda);
            
            const double defaultedNotional = nameNotional * unscheduledNotionalFraction;

            // Add to the total feeReductions the unscheduled reductions
            feeReductions->push_back(FeeLegReductionPerDefaultSP(
                new FeeLegReductionPerDefault(determinationDate,
                                              unscheduledCalcDate, // effective date
                                              unscheduledCalcDate, // calc date
                                              defaultedNotional, // defaulted notional
                                              defaultedNotional, // total notional
                                              recRate)));
        }

        // If there is notional that has not been delivered and that is not
        // expected to be delivered (typically because the cut-off date has 
        // been reached and deliveries were still pending, but maybe because
        // the moreDeliveriesPending flag has been set to false), recover that 
        // notional on cut-off date, as if it had been delivered and RR=1
        const double undeliverableNotionalFraction = 
            notionalFractionSettlingPhysically -
            (deliveredProportionSoFar + unscheduledNotionalFraction);

        if (!Maths::isZero(undeliverableNotionalFraction)) {
            const double undeliverableNotional = 
                nameNotional * undeliverableNotionalFraction;

            // Add to the total feeReductions the undeliverable reductions
            feeReductions->push_back(FeeLegReductionPerDefaultSP(
                new FeeLegReductionPerDefault(determinationDate,
                                              cutoffDate, // effective date
                                              cutoffDate, // calc date
                                              undeliverableNotional, // defaulted notional
                                              undeliverableNotional, // total notional
                                              1.0))); // recovery rate
        }
    }
    else {
        // All deliveries settle at once
        const DateTime& calcDate = 
            commonCalculationDate(creditEventDate, lastTriggerDate, bda);

        const double totalLoss = getOverallLoss(nameNotional, recoveryRate);

        // Arrays of losses/recovered notional 
        feeReductions.reset(new FeeLegReductionPerDefaultArray(
            1, FeeLegReductionPerDefaultSP(
                new FeeLegReductionPerDefault(
                    determinationDate,
                    calcDate,
                    calcDate,
                    totalLoss,  // defaulted notional (see below)
                    nameNotional * notionalFractionSettlingPhysically,
                    0.0)))); // RR = 0 -> defaulted notional = loss

    }

    return feeReductions;
}


/** Returns the overall loss resulting from the settlement of this name */
double TranchePhysicalSettlementOverride::getOverallLoss(
    const double nameNotional,
    const double recoveryRate) const 
{
    const double unscheduledNotionalFraction = 
        getUnscheduledNotionalFraction();

    double totalExpectedNotionalFraction = getTotalExpectedNotional();

    // Check if we are overriding the recovery rate here
    double recRate = (!recovery) ?
        recoveryRate : // No override, use argument passed in to this method
        recovery->doubleValue();
                                          
    // Initialise totalLoss with the unscheduled notional
    double totalLoss = nameNotional * unscheduledNotionalFraction * (1.0 - recRate);

    const int numDeliveries = (!deliveries) ? 0 : deliveries->size();
    double deliveredProportionSoFar = 0.0;
    
    for (int i=0; i < numDeliveries; ++i) {
        double deliveryNotionalFraction = Maths::min(
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
double TranchePhysicalSettlementOverride::getOverallRecoveryRate(
    double recoveryRate) const 
{
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
double TranchePhysicalSettlementOverride::getOverallDefaultedNotionalFraction() const {
    double notional = getNotionalScheduledForDelivery();
    double totalExpectedNotionalFraction = getTotalExpectedNotional();

    if (notional > totalExpectedNotionalFraction) {
        notional = totalExpectedNotionalFraction;
    }
    else {
        double unscheduledNotional = getUnscheduledNotionalFraction();
        // Note getUnscheduledNotionalFraction already takes into account that
        // the "delivered notional" cannot be greater than the "nops Notional"
        notional += unscheduledNotional;
    }
    return notional;
}


IObject* TranchePhysicalSettlementOverride::defaultTranchePhysicalSettlementOverride() {
    return new TranchePhysicalSettlementOverride();
}


void TranchePhysicalSettlementOverride::load (CClassSP& clazz) {
    clazz->setPublic(); // make visible to EAS/spreadsheet
    REGISTER(TranchePhysicalSettlementOverride, clazz);
    SUPERCLASS(TrancheCreditEventOverride);
    EMPTY_SHELL_METHOD(defaultTranchePhysicalSettlementOverride);

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

    FIELD(defaultToCalculationDelay, "An estimate for the delay between a "
          "credit event and the calculation date for all obligations (if "
          "field settleDeliveriesSeparately is false) or for the remaining "
          "expected obligations (if true)");
    FIELD_MAKE_OPTIONAL(defaultToCalculationDelay);

    FIELD(calculationDate, "Date when the recovery rate is determined "
                 "if NOT settling deliveries separately");
    FIELD_MAKE_OPTIONAL(calculationDate);

    FIELD(cutoffDate, "The last date when deliveries are accepted. "
                 "If the 'Recover notional' flag is set, any outstanding "
                 "notional not delivered on or before this date will "
                 "be deemed recovered.");
    FIELD_MAKE_OPTIONAL(cutoffDate);
    
    FIELD(moreDeliveriesPending, "Indicates whether more deliveries are "
                 "to be expected before the cut-off date.");
    FIELD_MAKE_OPTIONAL(moreDeliveriesPending);

    FIELD(notionalFractionSettlingPhysically, "Fraction of the notional "
                 "settling physically");
    FIELD_MAKE_TRANSIENT(notionalFractionSettlingPhysically);
}


CClassConstSP const TranchePhysicalSettlementOverride::TYPE = 
    CClass::registerClassLoadMethod("TranchePhysicalSettlementOverride", 
                                    typeid(TranchePhysicalSettlementOverride), 
                                    load);

// Array has to have its own type
DEFINE_TEMPLATE_TYPE(TranchePhysicalSettlementOverrideArray);


/** Included in ProductsLib to force the linker to include this file */
bool TranchePhysicalSettlementOverrideLoad() {
    return (TranchePhysicalSettlementOverride::TYPE != 0);
}

DRLIB_END_NAMESPACE
