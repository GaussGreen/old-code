//----------------------------------------------------------------------------
//
//   Group       : Quatitative Research
//
//   Filename    : CashSettlementOverrideName.cpp
//
//   Description : Class used by credit instruments (e.g., CIS) in order 
//                 to override a name's default-related parameters, when
//                 the default is settled in cash.
//
//   Author      : Jose Hilera
//
//   Date        : February 2006
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/Format.hpp"
#include "edginc/CashSettlementOverrideName.hpp"
#include "edginc/IBadDayAdjuster.hpp"
#include "edginc/ICreditContingentLeg.hpp"
#include "edginc/CtgLegLossPerDefault.hpp"

DRLIB_BEGIN_NAMESPACE


CashSettlementOverrideName::CashSettlementOverrideName() : 
    DetailedCreditEventOverrideName(TYPE), notionalFraction(1.0)
{}


/** Public constructor, with all fields passed in */
CashSettlementOverrideName::CashSettlementOverrideName(
        ICDSParSpreadsWrapper name,
        DateTime eventDeterminationDate,
        CIntSP triggerDelay,
        CDoubleSP recovery,
        double notionalFraction,
        CIntSP defaultToSettlementDelay,
        DateTime calculationDate,
        int calculationToSettlementDelay) :
    DetailedCreditEventOverrideName(TYPE,
                                    name, 
                                    eventDeterminationDate,
                                    triggerDelay,
                                    recovery),
    notionalFraction(notionalFraction),
    defaultToSettlementDelay(defaultToSettlementDelay),
    calculationDate(calculationDate),
    calculationToSettlementDelay(calculationToSettlementDelay)
{}


CashSettlementOverrideName::~CashSettlementOverrideName()
{}


/** Called immediately after object constructed */
void CashSettlementOverrideName::validatePop2Object() {
    static const string method = 
        "CashSettlementOverrideName::validatePop2Object";

    // Validate fields in the parent class
    DetailedCreditEventOverrideName::validatePop2Object();

    // Validate fields in this class
    if (notionalFraction < 0.0 || notionalFraction > 1.0) {
        throw ModelException(method,
                             "notionalFraction value for name " +
                             name.getName() + " is " +
                             Format::toString(notionalFraction) +
                             ". Value must be in [0, 1].");
    }
    
    // Either the calculationDate or the defaultToSettlementDelay need 
    // to be provided. calculationDate will be preferred if provided
    if (calculationDate.empty()) {
        if (!defaultToSettlementDelay) {
            throw ModelException(method,
                                 "Need to enter calculationDate or "
                                 "defaultToSettlementDelay for name " +
                                 name.getName());
        }
        else {
            // Will use the defaultToSettlementDelay, so verify it is right
            if (defaultToSettlementDelay->intValue() < 0) {
                throw ModelException(method,
                                     name.getName() + 
                                     "'s defaultToSettlementDelay is " + 
                                     Format::toString(defaultToSettlementDelay->intValue()) +
                                     ", and negative delays are not accepted.");
            }
        }
    }

    // Cross-validate fields in the parent class and in this class
    if (!calculationDate.empty()) {
        if (eventDeterminationDate.empty()) {
            throw ModelException(method,
                                 "Name " + name.getName() +
                                 " has not been triggered (eventDeterminationDate "
                                 "is not set) so its calculationDate cannot be "
                                 "known. Either mark the name as triggered or "
                                 "remove its calculationDate.");
        }
        if (calculationDate < eventDeterminationDate) {
            throw ModelException(method,
                                 "Name " + name.getName() +
                                 " has been triggered on " +
                                 eventDeterminationDate.toString() +
                                 " so its calculationDate cannot be before "
                                 "that (it is " + 
                                 calculationDate.toString() + ").");
        }
    }
}


/** Returns the contingent leg loss cashflows corresponding to a defaulted 
    name, assuming that the default settles in one date. */
CtgLegLossPerDefaultArraySP
CashSettlementOverrideName::historicContingentLegLosses(
    const double nameNotional,
    const double recoveryRate,
    const DateTime& lastTriggerDate,
    const DateTime& creditEventDate,
    IBadDayAdjusterConstSP bda) const
{
    static const string method( 
        "CashSettlementOverrideName::historicContingentLegLosses");

    const DateTime& determinationDate = getEventDeterminationDate(
        creditEventDate, lastTriggerDate, bda);

    if (determinationDate.empty() || (determinationDate > lastTriggerDate)) {
        // The default has not been triggered, or it was triggered after 
        // the lastTriggerDate -> no contingent payment will be made
        return CtgLegLossPerDefaultArraySP();
    }

    // Estimate for the delivery date
    DateTime defaultPayDate = getSettlementDate(creditEventDate,
                                                lastTriggerDate,
                                                valueDate,
                                                bda);

    if (defaultPayDate <= valueDate) {
        // Since this contingent payment has been paid already, return a null SP
        return CtgLegLossPerDefaultArraySP();
    }

    // Check if we are overriding the recovery rate here
    double recRate = getOverallRecoveryRate(recoveryRate);

    // Obtain the percentage of the notional that has been triggered 
    // for protection
    const double defaultedNotionalFraction = getOverallDefaultedNotionalFraction();
    const double defaultedNotional = nameNotional * defaultedNotionalFraction;

    CtgLegLossPerDefaultArraySP contingentReduction(
        new CtgLegLossPerDefaultArray(
            1, CtgLegLossPerDefaultSP(new CtgLegLossPerDefault(
                creditEventDate,
                defaultPayDate,
                defaultedNotional * (1.0 - recRate)))));

    return contingentReduction;
}


/** Returns the accrual pay date. For cash settlements this is the same
 * as the settlement date. */
const DateTime CashSettlementOverrideName::getAccrualPayDate(
    const DateTime& creditEventDate,
    const DateTime& lastTriggerDate,
    const DateTime& valueDate,
    IBadDayAdjusterConstSP bda) const
{
    return getSettlementDate(creditEventDate, lastTriggerDate, valueDate, bda);
}


/** Returns the date of the cash settlement based on the information in 
 * this override - it can be an estimate or the actual date if known */
const DateTime CashSettlementOverrideName::getSettlementDate(
    const DateTime& creditEventDate,
    const DateTime& lastTriggerDate,
    const DateTime& valueDate,
    IBadDayAdjusterConstSP bda) const
{
    static const string method = 
        "CashSettlementOverrideName::getSettlementDate";

    DateTime settlementDate;
    if (calculationDate.empty()) {
        const DateTime& determinationDate = getEventDeterminationDate(
            creditEventDate, lastTriggerDate, bda);

        if (determinationDate.empty()) {
            // An empty determinationDate means the default can not be trigered.
            // Therefore return an empty settlement date too.
            // Note settlementDate is an empty date already.
        }
        else {
#ifdef QLIB_REMOVE_EDD_VALIDATION_CDS
            // CalculationDate is empty, so need to check that the determinationDate
            // falls after the creditEventDate
            if (creditEventDate > determinationDate) {
                throw ModelException(method,
                                     name.getName() + 
                                     "'s event determination date (" +
                                     determinationDate.toString() +
                                     ") can not be before the default date ("+
                                     creditEventDate.toString() + 
                                     ") because the calculation date has not "
                                     "been specified.");
            }
#endif
            settlementDate = bda->badDayAdjust(
                creditEventDate.rollDate(defaultToSettlementDelay->intValue()));
       
            // If this date is before valueDate or determinationDate (whichever is 
            // greater) it means that the estimate for defaultToSettlementDelay was
            // short. Will use that date
            const DateTime& adjValueDate = 
                bda->badDayAdjust(valueDate.rollDate(1)); // JLH. make this "1" configurable
            const DateTime& maxDate = determinationDate.max(adjValueDate);
            if (settlementDate < maxDate) {
                settlementDate = maxDate;
            }
        }
    }
    else {
        // We know what the settlement date is
        // Note: if calculationToSettlementDelay is 0 and calculationDate falls 
        // on a holiday, so will settlementDate. Though bizarre, this is 
        // probably what is intended? Note there is no possible value for
        // calculationToSettlementDelay that will make settlementDate to fall
        // on the date immediately after the holiday.
        settlementDate = bda->addBusinessDays(calculationDate, 
                                              calculationToSettlementDelay);
    }
    return settlementDate;
}


/* Returns the overall recovery rate accross all settlements, such that
 * total loss = overall defaulted notional * (1-overall recovery rate) */
double CashSettlementOverrideName::getOverallRecoveryRate(
    double recoveryRate) const
{
    return getOverridenRecoveryRate(recoveryRate);
}


/** Returns the percentage of the notional that has been triggered 
 * for protection */
double CashSettlementOverrideName::getOverallDefaultedNotionalFraction() const {
    return notionalFraction;
}


IObject* CashSettlementOverrideName::defaultCashSettlementOverrideName() {
    return new CashSettlementOverrideName();
}


void CashSettlementOverrideName::load (CClassSP& clazz) {
    clazz->setPublic(); // make visible to EAS/spreadsheet
    REGISTER(CashSettlementOverrideName, clazz);
    SUPERCLASS(DetailedCreditEventOverrideName);
    EMPTY_SHELL_METHOD(defaultCashSettlementOverrideName);

    FIELD(notionalFraction,"Percentage of the notional being "
          "cash-settled (eg, 0.5 for 50%)");
    FIELD_MAKE_OPTIONAL(notionalFraction);
    FIELD(defaultToSettlementDelay, "Delay between credit event and "
          "settlement, in calendar days");
    FIELD_MAKE_OPTIONAL(defaultToSettlementDelay);
    FIELD(calculationDate, "Date when the recovery rate is determined");
    FIELD_MAKE_OPTIONAL(calculationDate);
    FIELD(calculationToSettlementDelay, "Delay between calculation date "
          "and settlement, in business days");
}


CClassConstSP const CashSettlementOverrideName::TYPE = 
    CClass::registerClassLoadMethod("CashSettlementOverrideName", 
                                    typeid(CashSettlementOverrideName), 
                                    load);

// Array has to have its own type
DEFINE_TEMPLATE_TYPE(CashSettlementOverrideNameArray);


/** Included in ProductsLib-modified::linkInClasses() via the productSrcsMap
 * script to force the linker to include this file */
bool CashSettlementOverrideNameLoad() {
    return (CashSettlementOverrideName::TYPE != 0);
}

DRLIB_END_NAMESPACE
