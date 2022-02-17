//----------------------------------------------------------------------------
//
//   Group       : Quatitative Research
//
//   Filename    : TrancheCashSettlementOverride.cpp
//
//   Description : Class used by CDO in order 
//                 to override a name's default-related parameters, when
//                 the default is settled in cash.
//
//   Author      : Jose Hilera
//
//   Date        : April 2006
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/Format.hpp"
#include "edginc/CreditFeeLeg.hpp"
#include "edginc/IBadDayAdjuster.hpp"
#include "edginc/IRebateCalculator.hpp"
#include "edginc/FeeLegReductionPerDefault.hpp"
#include "edginc/CtgLegLossPerDefault.hpp"
#include "edginc/TrancheCashSettlementOverride.hpp"


DRLIB_BEGIN_NAMESPACE


TrancheCashSettlementOverride::TrancheCashSettlementOverride() : 
    TrancheCreditEventOverride(TYPE),
    defaultedNotionalFraction(1.0), 
    notionalFraction(1.0)
{}


/** Public constructor, with all fields passed in */
TrancheCashSettlementOverride::TrancheCashSettlementOverride(
    CIntSP triggerDelay,
    DateTime eventDeterminationDate,
    CDoubleSP recovery,
    DateTime valueDate,
    CDoubleSP defaultedNotionalFractionSP,
    CIntSP defaultToCalculationDelay,
    DateTime calculationDate,
    double notionalFraction) :
        TrancheCreditEventOverride(TYPE,
                                   triggerDelay, 
                                   eventDeterminationDate,
                                   recovery,
                                   valueDate),
        defaultToCalculationDelay(defaultToCalculationDelay),
        calculationDate(calculationDate),
        notionalFraction(notionalFraction)
{
    defaultedNotionalFraction = !defaultedNotionalFractionSP ?
        notionalFraction : // default: all notional settling in cash
        defaultedNotionalFractionSP->doubleValue();
}


TrancheCashSettlementOverride::~TrancheCashSettlementOverride()
{}


/** Called immediately after object constructed */
void TrancheCashSettlementOverride::validatePop2Object() {
    static const string method = 
        "TrancheCashSettlementOverride::validatePop2Object";

    // Validate fields in the parent class
    TrancheCreditEventOverride::validatePop2Object();

    if (notionalFraction < 0.0 || notionalFraction > 1.0) {
        throw ModelException(method,
                             "notionalFraction is " +
                             Format::toString(notionalFraction) +
                             ". Value must be in [0,1].");
    }

    // Validate fields in this class
    if ((defaultedNotionalFraction < 0.0) || 
        (defaultedNotionalFraction > notionalFraction))
    {
        throw ModelException(method,
                             "defaultedNotionalFraction is " +
                             Format::toString(defaultedNotionalFraction) +
                             ". Value must be in [0, " +
                             Format::toString(notionalFraction) + "].");
    }
    
    // Either the calculationDate or the defaultToCalculationDelay need 
    // to be provided. calculationDate will be preferred if provided
    if (calculationDate.empty()) {
        if (!defaultToCalculationDelay) {
            throw ModelException(method,
                                 "Need to enter calculationDate or "
                                 "defaultToCalculationDelay");
        }
        else {
            // Will use the defaultToCalculationDelay, so verify it is right
            if (defaultToCalculationDelay->intValue() < 0) {
                throw ModelException(method,
                                     "DefaultToCalculationDelay is " + 
                                     Format::toString(defaultToCalculationDelay->intValue()) +
                                     ", and negative delays are not accepted.");
            }
        }
    }
    else { // Cross-validate fields in the parent class and in this class
        if (eventDeterminationDate.empty()) {
            throw ModelException(method,
                                 "Name has not been triggered "
                                 "(eventDeterminationDate is not set) so its "
                                 "calculationDate cannot be known. Either mark "
                                 "the name as triggered or remove its "
                                 "calculationDate.");
        }
        if (calculationDate < eventDeterminationDate) {
            throw ModelException(method,
                                 "Name has been triggered on " +
                                 eventDeterminationDate.toString() +
                                 " so its calculationDate cannot be before "
                                 "that (it is " + 
                                 calculationDate.toString() + ").");
        }
    }
}


/** Returns the contingent leg cashflows corresponding to a defaulted name,
 * assuming that the default settles in one date. */
CtgLegLossPerDefaultArraySP TrancheCashSettlementOverride::historicContingentLegLosses(
    const double nameNotional,
    const double recoveryRate,
    const DateTime& lastTriggerDate,
    const DateTime& creditEventDate,
    IBadDayAdjusterConstSP bda) const
{
    static const string method(
        "TrancheCashSettlementOverride::historicContingentLegLosses");
    
    // Estimate the calculation date
    DateTime calcDate = getCalculationDate(creditEventDate,
                                           lastTriggerDate,
                                           bda);

    CtgLegLossPerDefaultArraySP ctgLegLoss;
    if (calcDate.empty()) {
        // The default has not been triggered, or it was triggered after 
        // the lastTriggerDate -> no losses
        ctgLegLoss.reset(new CtgLegLossPerDefaultArray(0));
    }
    else {
        // Obtain the percentage of the notional that has been triggered 
        // for protection
        const double nameLoss = getLossAmount(nameNotional, recoveryRate);

        CashFlowArraySP loss(new CashFlowArray(1, CashFlow(calcDate, nameLoss)));

        ctgLegLoss = 
            CtgLegLossPerDefault::produceLossesPerDefault(creditEventDate, loss);
    }
    return ctgLegLoss;
}


/* Returns the loss amount associated to this name's default.
 * JLH - should use the lbd params, as per PortfolioName::lossGivenDefault()
 * Since it does not, we have had to add a stricter check to 
 * PortfolioName::validatePop2Object() to verify that the lbd params do
 * not enter the equation - if required, this method would have to change */
double TrancheCashSettlementOverride::getLossAmount(
    const double nameNotional,
    const double recoveryRate) const 
{
    // Obtain the percentage of the notional that has been triggered 
    // for protection
    const double defaultedNotionalFraction = getOverallDefaultedNotionalFraction();
    const double defaultedNotional = nameNotional * defaultedNotionalFraction;
    const double recRate = getOverallRecoveryRate(recoveryRate);
    
    return defaultedNotional * (1.0 - recRate);
}


/** Returns the portfolio reductions for the fee leg (corresponding to losses
 * and recovered amounts). */
FeeLegReductionPerDefaultArraySP 
TrancheCashSettlementOverride::historicFeeLegReductions(
    const double nameNotional,
    const double recoveryRate,
    const DateTime& lastTriggerDate,
    const DateTime& creditEventDate,
    IBadDayAdjusterConstSP bda) const
{
    FeeLegReductionPerDefaultArraySP feeReductions;

    /* First, obtain the settlement parameters using the information in
     * this override */
    const DateTime& determinationDate = getEventDeterminationDate(
        creditEventDate, lastTriggerDate, bda);

    if (determinationDate.empty()) {
        // The default can not be triggered, or it was triggered after 
        // the lastTriggerDate -> no losses or recovered amounts
    }
    else {
        // Get the calculation date
        DateTime calcDate = 
            getCalculationDate(creditEventDate, lastTriggerDate, bda);

        // Obtain the percentage of the notional that has been triggered 
        // for protection
        const double defaultedNotionalFraction = getOverallDefaultedNotionalFraction();
        const double recRate = getOverallRecoveryRate(recoveryRate);

        feeReductions.reset(new FeeLegReductionPerDefaultArray(
            1, FeeLegReductionPerDefaultSP(new FeeLegReductionPerDefault(
                determinationDate,  
                calcDate,     // effective date
                calcDate,     // calculation date
                nameNotional * defaultedNotionalFraction,
                nameNotional * notionalFraction,
                recRate))));
    }

    return feeReductions;
}


/** Returns the date of the calculation date based on the information in 
 * this override - it can be an estimate or the actual date if known */
DateTime TrancheCashSettlementOverride::getCalculationDate(
    const DateTime& creditEventDate,
    const DateTime& lastTriggerDate,
    IBadDayAdjusterConstSP bda) const
{
    static const string method(
        "TrancheCashSettlementOverride::getCalculationDate");
    
    DateTime calcDate;
    if (calculationDate.empty()) {
        const DateTime& determinationDate = getEventDeterminationDate(
            creditEventDate, lastTriggerDate, bda);

        if (determinationDate.empty()) {
            // An empty determinationDate means the default has not been trigered.
            // Therefore return an empty calculation date date too.
            // Note calcDate is an empty date already.
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

            calcDate = 
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
        calcDate = calculationDate;
    }

    return calcDate;
}




/* Returns the overall recovery rate accross all settlements, such that
 * total loss = overall defaulted notional * (1-overall recovery rate) */
double TrancheCashSettlementOverride::getOverallRecoveryRate(
    double recoveryRate) const 
{
    return (!recovery ? recoveryRate : recovery->doubleValue());
}


/** Returns the percentage of the notional that has been triggered 
 * for protection */
double TrancheCashSettlementOverride::getOverallDefaultedNotionalFraction() const {
    return defaultedNotionalFraction;
}


IObject* TrancheCashSettlementOverride::defaultTrancheCashSettlementOverride() {
    return new TrancheCashSettlementOverride();
}


void TrancheCashSettlementOverride::load (CClassSP& clazz) {
    clazz->setPublic(); // make visible to EAS/spreadsheet
    REGISTER(TrancheCashSettlementOverride, clazz);
    SUPERCLASS(TrancheCreditEventOverride);
    EMPTY_SHELL_METHOD(defaultTrancheCashSettlementOverride);

    FIELD(defaultedNotionalFraction, 
                 "Percentage of the notional triggered for protection");
    FIELD_MAKE_OPTIONAL(defaultedNotionalFraction);
    FIELD(defaultToCalculationDelay, 
          "Delay between credit event and settlement");
    FIELD_MAKE_OPTIONAL(defaultToCalculationDelay);
    FIELD(calculationDate, 
                 "Date when the recovery rate is determined");
    FIELD_MAKE_OPTIONAL(calculationDate);

    FIELD(notionalFraction,
                 "Percentage of the notional considered for cash settlement");
    FIELD_MAKE_TRANSIENT(notionalFraction);
}


CClassConstSP const TrancheCashSettlementOverride::TYPE = 
    CClass::registerClassLoadMethod("TrancheCashSettlementOverride", 
                                    typeid(TrancheCashSettlementOverride), 
                                    load);

// Array has to have its own type
DEFINE_TEMPLATE_TYPE(TrancheCashSettlementOverrideArray);

/** Included in ProductsLib to force the linker to include this file */
bool TrancheCashSettlementOverrideLoad() {
    return (TrancheCashSettlementOverride::TYPE != 0);
}

DRLIB_END_NAMESPACE
