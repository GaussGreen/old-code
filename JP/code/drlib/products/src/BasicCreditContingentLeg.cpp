//----------------------------------------------------------------------------
//
//   Filename    : BasicCreditContingentLeg.cpp
//
//   Description : Single period contingent leg
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/BasicCreditContingentLeg.hpp"
#include "edginc/IDecretionCurve.hpp"
#include "edginc/ICreditEventOverrideName.hpp"
#include "edginc/IBadDayAdjuster.hpp"
#include "edginc/CtgLegLossPerDefault.hpp"

DRLIB_BEGIN_NAMESPACE

//---------------------------------
// BasicCreditContingentLeg methods
//---------------------------------

BasicCreditContingentLeg::BasicCreditContingentLeg(
    const  DateTime& protectionStartDate,
    const  DateTime& protectionEndDate)
: CObject(TYPE),
  protectionStartDate(protectionStartDate),
  protectionEndDate(protectionEndDate)
{}

BasicCreditContingentLeg::~BasicCreditContingentLeg()
{}

//----------------
// CObject methods
//----------------

void BasicCreditContingentLeg::validatePop2Object()
{
    static const string method = "BasicCreditContingentLeg::validatePop2Object";

    //ensure end is later than start
    if (protectionStartDate.isGreater(protectionEndDate))
    {
        throw ModelException(method,"protection start must preceed protection end");
    }

    //ensure swap recovery was passed if useSwapRecovery set to true
    if (useSwapRecovery)
    {
        if (!swapRecovery)
        {
            throw ModelException(method,
                "swapRecovery must be provided if useSwapRecovery is set to true");
        }
    }
}

//-----------------------------
// ICreditContingentLeg methods
//-----------------------------

/** Price this contingent leg. */
double BasicCreditContingentLeg::price(
    double                       initialNotional,     // highStrike - lowStrike
    double                       outstandingNotional, /* initialTrancheSize -
                                                            pastTrancheLoss */
    const DateTime&              today,
    const DateTime&              valDateCF, // to be scrapped
    const IDiscountCurveRiskySP  effectiveCurve,
    const CashFlowArray&         pastTrancheLosses,
    const BoolArray&             payPastTrancheLosses,
    bool                         computeDebugPrices, /* true: populate arrays 
                                                        below */
    DoubleArray&                 debugUnitPrice, // price for each leg unit
    DoubleArray&                 debugUnitHistPrice, /* price for each leg unit
                                                        * due to historical default */
    IBadDayAdjusterConstSP        bda) const
{
    static const string method = "BasicCreditContingentLeg::price";
    //for now
    throw ModelException(method,"not yet implemented");
}

/** Returns the earliest observation start date */
DateTime BasicCreditContingentLeg::firstObservationStartDate() const
{
    return protectionStartDate;
}

/** Returns the last pay date */
DateTime BasicCreditContingentLeg::lastPayDate(IBadDayAdjusterConstSP bda) const
{
    return protectionEndDate;
}

/** Returns the last observation date */
DateTime BasicCreditContingentLeg::lastObservationEndDate() const
{
    return protectionEndDate;
}

/** When to stop tweaking */
DateTime BasicCreditContingentLeg::lastYCSensDate(const DateTime& currentLastDate,
                                                  IBadDayAdjusterConstSP bda) const
{
    return protectionEndDate;
}

CashFlowArraySP BasicCreditContingentLeg::generateKnownCashFlows(
    const DateTime&        today,
    double                 initialTrancheSize,
    const CashFlowArray&   pastTrancheLosses,
    const BoolArray&       payPastTrancheLosses,
    IBadDayAdjusterConstSP bda) const
{
    static const string method = "BasicCreditContingentLeg::generateKnownCashFlows";
    //for now
    throw ModelException(method,"not yet implemented");
}

/** CAUTION: HACK!
    * This method is for the benefit of the fee leg only: it returns
    * payment details required in the fee leg (because they should 
    * have been added to the CDO instrument in the first place, but
    * now it is too late to change them). Although it would be possible
    * to add these fields to the fee leg also, it has been estimated 
    * that this would cause confussion among the library users 
    * (specially the start/end arrays) and instead they will be
    * obtained from the contingent leg and passed to the fee leg as
    * required.
    * All parameters are really outputs of the method - the contents of
    * the smart pointers passed in will be discarded */
void BasicCreditContingentLeg::getPaymentInformation(BoolArraySP&     payAsYouGoArray,
                                                     IntArraySP&      numDelayDaysArray,
                                                     DateTimeArraySP& startDate,
                                                     DateTimeArraySP& endDate,
                                                     DateTimeArraySP& paymentDate) const
{
    static const string method = "BasicCreditContingentLeg::getPaymentInformation";

    payAsYouGoArray.reset(new BoolArray(1));
    numDelayDaysArray.reset(new IntArray(1));
    startDate.reset(new DateTimeArray(1));
    endDate.reset(new DateTimeArray(1));
    paymentDate.reset(new DateTimeArray(1)); //empty datetime

    (*payAsYouGoArray)[0] = true;
    (*numDelayDaysArray)[0] = 0;
    (*startDate)[0] = protectionStartDate;
    (*endDate)[0] = protectionEndDate;
}

/* Checks if the input date is covered for protection, i.e., falls in
    * one of the observation periods of this leg */
bool BasicCreditContingentLeg::isDateCoveredForProtection(const DateTime& date) const
{
    return date.within(protectionStartDate, protectionEndDate);
}

/**Compute PV of contingent leg at valuation date, with instrument default settlement and
    payment date behaviour.*/
double BasicCreditContingentLeg::getContingentLegPV(const DateTime&              valuationDate, 
                                                    const IDiscountCurveRisky&   crv,
                                                    IBadDayAdjusterConstSP bda) const
{
    return getContingentLegPV(valuationDate, valuationDate, crv, bda);
}

/**Compute PV of contingent leg at valuationDate, but for unconditional payment at 
    paymentDate, conditional on no new defaults before valuationDate.*/
double BasicCreditContingentLeg::getContingentLegPV(const DateTime&              valuationDate, 
                                                    const DateTime&              paymentDate, 
                                                    const IDiscountCurveRisky&   crv,
                                                    IBadDayAdjusterConstSP bda) const // bda is currently unused
{
    static const string method = "BasicCreditContingentLeg::getContingentLegPV";
    
    double sp = 1.0; // survival probability from valuationDate to paymentDate
    double df = 1.0; // IR discount from valuationDate to paymentDate
    double survivalPV = 0.0;
    double defaultPV  = 0.0;
    double pv         = 0.0;

   if(valuationDate < paymentDate)
   {
       sp = crv.survivalProb(valuationDate, paymentDate);
       if(sp == 0.0)
       {
           throw ModelException (method, "Unexpected zero survival probability.");
       }
       df = crv.pv(valuationDate, paymentDate);
       df /= sp;
   }
   else if (valuationDate > paymentDate)
   {
       throw ModelException(method, "valuationDate must always be on or before paymentDate.");
   }
   else if( valuationDate >= protectionEndDate)
   {   // passed end of protection date so value is zero
       return 0.0;
   }

   // Value of protection leg if credit survives until paymentDate
   survivalPV = crv.protectionPV(paymentDate,
                                 paymentDate.max(protectionStartDate),
                                 paymentDate.max(protectionEndDate),
                                 IDiscountCurveRisky::RECOVER_1, // recover 100% - scaled later
                                 0.0);
                                 
   //// Protection PV from valuationDate to paymentDate
   //defaultPV = crv.protectionPV(valuationDate,
   //                             valuationDate.max(protectionStartDate),
   //                             paymentDate.min(protectionEndDate),
   //                             IDiscountCurveRisky::RECOVER_1, // recover 100% - scaled later
   //                             paymentDate);

   pv = notional * (1- getRecovery(crv)) * (sp * survivalPV + defaultPV) * df;
   return pv;
}

/** return the recovery rate to be used */
double BasicCreditContingentLeg::getRecovery(const IDiscountCurveRisky& crv) const
{
    double recovery;

    if (!useSwapRecovery)
    {
        recovery = crv.getRecovery();
    }
    else
    {
        recovery = swapRecovery->doubleValue();
    }

    return recovery;
}

/** Returns the amount that would be recovered upon default */
double BasicCreditContingentLeg::recoveredValue(const DateTime& valueDate,
                                                const IDiscountCurveRisky& crv,
                                                const IDecretionCurveConstSP prepay,
                                                const IDiscountCurveRisky::RecoveryType recType) const
{
    double factor = prepay->getFactor(valueDate);
    double balance = prepay->pv(valueDate);  

    double recRate;
    switch(recType)
    {
        case IDiscountCurveRisky::RECOVER_R:
            recRate = getRecovery(crv);
            break;
        case IDiscountCurveRisky::RECOVER_1:
            recRate = 1.0;
            break;
        case IDiscountCurveRisky::RECOVER_1_MINUS_R:
            recRate = 1.0-getRecovery(crv);
            break;
        case IDiscountCurveRisky::RECOVER_0:
            recRate = 0.0;
            break;
    }


    double defaultValue = recRate * notional * factor * balance;
    return defaultValue;
}

/** Returns the amount that would be recovered upon default,
    using the supplied recovery rate */
double BasicCreditContingentLeg::recoveredValue(const DateTime& valueDate,
                                                const IDiscountCurveRisky& crv,
                                                const IDecretionCurveConstSP prepay,
                                                const double recoveryToUse,
                                                const IDiscountCurveRisky::RecoveryType recType) const
{
    double factor = prepay->getFactor(valueDate);
    double balance = prepay->pv(valueDate);  

    double recRate;
    switch(recType)
    {
        case IDiscountCurveRisky::RECOVER_R:
            recRate = recoveryToUse;
            break;
        case IDiscountCurveRisky::RECOVER_1:
            recRate = 1.0;
            break;
        case IDiscountCurveRisky::RECOVER_1_MINUS_R:
            recRate = 1.0-recoveryToUse;
            break;
        case IDiscountCurveRisky::RECOVER_0:
            recRate = 0.0;
            break;
    }

    double defaultValue = recRate * notional * factor * balance;
    return defaultValue;
}

/** Prices the leg sufferring a default, under the (optional) credit event */
double BasicCreditContingentLeg::getContingentLegDefaultedPV(
    const DateTime&              valuationDate, 
    const IDiscountCurveRisky&   crv,
    const IDiscountCurveConstSP  discount,
    const IDecretionCurveConstSP prepay,
    const DateTime&              defaultDate,
    IBadDayAdjusterConstSP       badDayAdjuster,
    const bool                   allowIncludingTodaysPayments,
    ICreditEventOverrideNameSP   creditEventOverride,
    CIntSP                       triggerDelay,
    CIntSP                       defaultToSettlementDelay,
    DateTime                     lastTriggerDate) const
{
    static const string method("BasicCreditContingentLeg::getContingentLegDefaultedPV");

    CashFlowArraySP contingentCashFlows(new CashFlowArray(0));
    if ((defaultDate < protectionStartDate) ||
        (defaultDate > protectionEndDate)) 
    {
        // Default happened outside the protection period, so 
        // contingentCashFlows is empty 
    }
    else {
        if (!creditEventOverride) {
            // Flag to indicate whether payments on valueDate should be included
            // in the valuation or not. In general they should not but for
            // compatibility with the old CredDefSwap implementation, if the
            // triggerDelay and  defaultToSettlementDelay are NOT passed in,
            // they will be included
            bool includeTodaysPayments = false;
            DateTime defaultPaymentDate;

            if (!defaultToSettlementDelay) {
                includeTodaysPayments = allowIncludingTodaysPayments;
                
                // No adjustment is done here.
                defaultPaymentDate = defaultDate;
            }
            else {
                // Note this means both triggerDelay and defaultToSettlementDelay
                // are not NULL.
                defaultPaymentDate = 
                    ITrancheCreditEventOverride::rollAndAdjustDate(
                        defaultDate,
                        defaultToSettlementDelay,
                        valuationDate,
                        valuationDate,
                        lastTriggerDate,
                        badDayAdjuster);
            }

            // Generate the one cashflow "manually"
            if ((defaultDate <= lastTriggerDate) &&
                ((defaultPaymentDate > valuationDate) ||
                 ((defaultPaymentDate == valuationDate) &&
                  includeTodaysPayments)))
            {
                // Notional < 0 means short protection
                // Notional > 0 means long protection
                contingentCashFlows->push_back(
                    CashFlow(defaultPaymentDate,
                             notional * (1.0 - getRecovery(crv))));
            }
            else {
                // Since this contingent payment would have been paid before
                // today or the default was triggered too late, return an 
                // empty cash flow array
            }
        }
        else {
            // Use the override
            CtgLegLossPerDefaultArraySP reductions = 
                creditEventOverride->historicContingentLegLosses(
                    notional,
                    getRecovery(crv),
                    lastTriggerDate,
                    defaultDate,
                    badDayAdjuster);

            CashFlowArraySP allCashFlows = 
                CtgLegLossPerDefault::getReductions(reductions);

            // Insert all future cashflows
            int numCashFlows = allCashFlows->size();
            contingentCashFlows->reserve(numCashFlows);
            for (int i=0; i < numCashFlows; ++i) {
                CashFlow cf((*allCashFlows)[i]); // for ease
                if (cf.date > valuationDate) {
                    contingentCashFlows->push_back(cf);
                }
            }
        }
    }

    double pv;
    double contingentLegValue = 0.0;
    auto_ptr<IYieldCurve::IKey> key(discount->logOfDiscFactorKey());
    // pv cashFlows
    for (int i=0; i < contingentCashFlows->size(); ++i) {
        CashFlow cf = (*contingentCashFlows)[i]; // For ease
        pv = exp(key->calc(valuationDate, cf.date));
        contingentLegValue += pv * cf.amount;
    }

    const double factor = prepay->getFactor(valuationDate);
    const double balance = prepay->pv(defaultDate);  

    contingentLegValue *= factor * balance;

    return contingentLegValue;
}

/** Returns the leg notional */
double BasicCreditContingentLeg::getContingentLegNotional() const
{
    return notional;
}

/** Sets the leg notional */
void BasicCreditContingentLeg::setContingentLegNotional(double newNotional)
{
    notional = newNotional;
}

//--------------------------------------
// ICreditContingentLegGenerator methods
//--------------------------------------

ICreditContingentLegSP BasicCreditContingentLeg::generateCreditContingentLeg(
    const DateTime& startDate,
    const DateTime& endDate) const
{
    //just build a new BasicCreditContingentLeg
    BasicCreditContingentLegSP newLeg = BasicCreditContingentLegSP(
        new BasicCreditContingentLeg(startDate, endDate));

    return newLeg;
}

//-------------------------------------------
// BasicCreditContingentLeg methods (private)
//-------------------------------------------

void BasicCreditContingentLeg::load(CClassSP& clazz)
{
    clazz->setPublic();
    REGISTER(BasicCreditContingentLeg, clazz);
    SUPERCLASS(CObject);
    IMPLEMENTS(ICreditContingentLeg);
    IMPLEMENTS(ICreditContingentLegGenerator);
    EMPTY_SHELL_METHOD(defaultConstructor);

    FIELD(notional,            "notional (>0 implies long protection)");
    FIELD(protectionStartDate, "protection start date");
    FIELD(protectionEndDate,   "protection end date");
    FIELD(useSwapRecovery,     "override underlying recovery rate");
    FIELD       (swapRecovery,        "recovery rate to use if overridden");

    FIELD_MAKE_OPTIONAL(swapRecovery);  //required if useSwapRecovery=TRUE
}

IObject* BasicCreditContingentLeg::defaultConstructor()
{
    return new BasicCreditContingentLeg();
}

BasicCreditContingentLeg::BasicCreditContingentLeg()
: CObject(TYPE)
{}

CClassConstSP const BasicCreditContingentLeg::TYPE = 
CClass::registerClassLoadMethod("BasicCreditContingentLeg", typeid(BasicCreditContingentLeg), load);

bool BasicCreditContingentLegLoad(){
    return BasicCreditContingentLeg::TYPE != NULL;
}

DRLIB_END_NAMESPACE
