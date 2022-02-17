//----------------------------------------------------------------------------
//
//   Group       : CH Quantitative Research
//
//   Filename    : DeliveryDetails.cpp
//
//   Description : Class used to contain the information regarding a delivery
//                 for physically settled defaults
//
//   Author      : Jose Hilera
//
//   Date        : May 2006
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/Format.hpp"
#include "edginc/IBadDayAdjuster.hpp"
#include "edginc/DeliveryDetails.hpp"
#include "edginc/ITrancheCreditEventOverride.hpp"


DRLIB_BEGIN_NAMESPACE

DeliveryDetails::DeliveryDetails() : CObject(TYPE)
{}

DeliveryDetails::~DeliveryDetails()
{}


/** Called immediately after object constructed */
void DeliveryDetails::validatePop2Object() {
    static const string method = "DeliveryDetails::validatePop2Object";

    double numOfDeliveryAmounts = deliveryAmountFractions->size();
    if (numOfDeliveryAmounts == 0) {
        throw ModelException(method,
                             "Delivery on date " + deliveryDate.toString() +
                             "has no deliveryAmountFractions. Do not specify "
                             "the delivery if its info is not yet known.");
    }

    for (int i=0; i < numOfDeliveryAmounts; ++i) {
        if ((*deliveryAmountFractions)[i] < 0) {
            throw ModelException(method,
                                 "Delivery on date " + deliveryDate.toString() +
                                 ": deliveryAmountFractions[" + 
                                 Format::toString(i) + "] = " + 
                                 Format::toString((*deliveryAmountFractions)[i]) +
                                 ". Negative amounts are not allowed.");
        }
    }

    if (!!recovery) {
        double dRecovery = recovery->doubleValue();
        if (dRecovery < 0.0 || dRecovery > 1.0) {
            throw ModelException(method,
                                 "Overall recovery for delivery on date " +
                                 deliveryDate.toString() +
                                 " is out of [0,1] (recovery = " +
                                 Format::toString(dRecovery) + ")");
        }
    }

    if (deliveryDate.empty()) {
        throw ModelException(method,
                             "The delivery date is empty. If there is no "
                             "planned delivery do not enter its details, but "
                             "if it is the delivery date must be entered");
    }

    if (!calculationDate.empty() && calculationDate < deliveryDate) {
        throw ModelException(method,
                             "The calculation date (" +
                             calculationDate.toString() +
                             ") must fall after the delivery date (" +
                             deliveryDate.toString() + ")");
    }

    if (!!defaultToCalculationDelay && defaultToCalculationDelay->intValue() < 0) {
        throw ModelException(method,
                             "defaultToCalculationDelay is " + 
                             Format::toString(defaultToCalculationDelay->intValue()) +
                             ", and negative delays are not accepted.");
    }
}


/** Returns the % notional expected to arrive during this 
 * name's settlement period */
double DeliveryDetails::getNotionalFractionInDelivery() const
{
    double notionalFraction = 0.0;
    for (int i=0; i < deliveryAmountFractions->size(); ++i) {
        notionalFraction += (*deliveryAmountFractions)[i];
    }
    return notionalFraction;
}


/** Returns the delivery date of this delivery */
DateTime DeliveryDetails::getDeliveryDate() const {
    return deliveryDate;
}

/** Returns the recovery rate for this delivery - the parameter
 * is the default recovery rate, in case it is not overriden here */
double DeliveryDetails::getDeliveryRecoveryRate(const double recoveryRate) const {
    if (!recovery) {
        return recoveryRate;
    }
    else {
        return recovery->doubleValue();
    }
}

/** Returns the (potentially estimate) calculation date for this
 * recovery, as used for CDO tranches */
DateTime DeliveryDetails::deliveryCalcDate(const DateTime& creditEventDate,
                                           const DateTime& determinationDate,
                                           const DateTime& valueDate,
                                           IBadDayAdjusterConstSP bda) const 
{
    static const string method = "DeliveryDetails::deliveryCalcDate";
    
    DateTime calcDate;
    if (determinationDate.empty()) {
        ; // Not triggered, so no calculation date (calcDate is already empty)
    }
    else if (calculationDate.empty()) {
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

        calcDate = ITrancheCreditEventOverride::rollAndAdjustDate(
            creditEventDate,
            defaultToCalculationDelay,
            valueDate,
            determinationDate,
            bda);

        // Check for invalid configurations
        if (calcDate < deliveryDate) {
            throw ModelException(method,
                                 "The delivery on date " + 
                                 deliveryDate.toString() +
                                 " has a calculation date on "
                                 + calcDate.toString() + 
                                 ". This is not allowed.");
        }
    }
    else {
        calcDate = calculationDate;
    }
    return calcDate;
}



/** Invoked when Class is 'loaded' */
void DeliveryDetails::load(CClassSP& clazz) {
    clazz->setPublic();  // make visible to EAS/spreadsheet
    REGISTER(DeliveryDetails, clazz);
    SUPERCLASS(CObject);
    EMPTY_SHELL_METHOD(defaultDeliveryDetails);

    FIELD(deliveryDate,       "Date for this delivery");
    FIELD(deliveryAmountFractions,   "Notice of Physical Delivery - Percentage "
                                     "of notional amounts that will be delivered");
    FIELD(recovery,                  "Overall RR for this delivery");
    FIELD(defaultToCalculationDelay, "Delay between the credit event date and "
                                     "the calculation date for this delivery "
                                     "(if settling deliveries separately");
    FIELD(calculationDate,    "Date when the recovery rate is determined "
                                     "for the obligations in this delivery (if "
                                     "settling deliveries separately");
    FIELD_MAKE_OPTIONAL(recovery);
    FIELD_MAKE_OPTIONAL(defaultToCalculationDelay);
    FIELD_MAKE_OPTIONAL(calculationDate);
}


IObject* DeliveryDetails::defaultDeliveryDetails() {
    return new DeliveryDetails();
}

CClassConstSP const DeliveryDetails::TYPE = 
    CClass::registerClassLoadMethod("DeliveryDetails", 
                                    typeid(DeliveryDetails),
                                    load);

// Array has to have its own type
DEFINE_TEMPLATE_TYPE(DeliveryDetailsArray);

/** Included in ProductsLib to force the linker to include this file */
bool DeliveryDetailsLoad() {
    return (DeliveryDetails::TYPE != 0);
}

DRLIB_END_NAMESPACE
