//----------------------------------------------------------------------------
//
//   Group       : Derivatives Research
//
//   Filename    : VanillaCreditContingentLeg.cpp
//
//   Description : Implementation of a contingent leg on a Credit Vanilla Instrument :
//
//   Author      : Doris Morris
//
//   Date        : February 2006
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/VanillaCreditContingentLeg.hpp"
#include "edginc/SensControl.hpp"
#include "edginc/IDecretionCurve.hpp"
#include "edginc/ICreditEventOverrideName.hpp"
#include "edginc/IBadDayAdjuster.hpp"

DRLIB_BEGIN_NAMESPACE

//-----------------------------------
// VanillaCreditContingentLeg methods
//-----------------------------------

VanillaCreditContingentLeg::VanillaCreditContingentLeg(const  DateTime& valueDate,
                                                       const  DateTime& startDate,
                                                       const  DateTime& endDate,
                                                       double baseNotional,
                                                       double recoveryRate,
                                                       double delay,
                                                       const  YieldCurveWrapper discountCrv,
                                                       const  BadDayConventionSP pmtBdc):CObject(TYPE),
                                                       baseNotional(baseNotional),
                                                       valueDate(valueDate),
                                                       effectiveDate(startDate),
                                                       maturityDate(endDate),
                                                       recoveryRate(recoveryRate),
                                                       delay(delay),
                                                       pmtBdc(pmtBdc),
                                                       discount(discountCrv)
{
    /* A VanillaCreditContingentLeg is always 1-R */
    recoveryTypeString = "RECOVER_1_MINUS_R";
    recovery = IDiscountCurveRisky::RECOVER_1_MINUS_R;
}

VanillaCreditContingentLeg::~VanillaCreditContingentLeg(){}

/* Set the base notional*/
void VanillaCreditContingentLeg::setNotional(double newNotional)
{
   const static string method = "VanillaCreditContingentLeg::setNotional";

   if(baseNotional == 0.0)
   {
       throw ModelException(method, "You cannot reset the notional of a VanillaCreditContingentLeg which has zero notional!");
   }

    baseNotional = newNotional;
}

//----------------
// CObject methods
//----------------

/** Called immediately after object constructed */
void VanillaCreditContingentLeg::validatePop2Object() 
{
    static const string method("VanillaCreditContingentLeg::validatePop2Object");

    if(effectiveDate > maturityDate)
    {
        throw ModelException(method, "protection start date (" +
                             effectiveDate.toString() + 
                             ") should be before protectionEndDate (" +
                             maturityDate.toString() + ").");
    }
}

//-----------------------------
// ICreditContingentLeg methods
//-----------------------------

double VanillaCreditContingentLeg::price(
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
    IBadDayAdjusterConstSP       bda) const
{
    const static string method("VanillaCreditContingentLeg::price");
    //for now
    throw ModelException(method,"not yet implemented");
}

/** Returns the earliest observation start date */
DateTime VanillaCreditContingentLeg::firstObservationStartDate() const
{
    const static string method("VanillaCreditContingentLeg::firstObservationStartDate");
    //for now
    return effectiveDate;
}

/** Returns the last pay date */
DateTime VanillaCreditContingentLeg::lastPayDate(IBadDayAdjusterConstSP bda) const
{
    const static string method("VanillaCreditContingentLeg::lastPayDate");
    //for now
    throw ModelException(method,"not yet implemented");
}

/** Returns the last observation date */
DateTime VanillaCreditContingentLeg::lastObservationEndDate() const
{
    const static string method("VanillaCreditContingentLeg::lastObservationDate");
    
    return maturityDate;
}

/** When to stop tweaking */
DateTime VanillaCreditContingentLeg::lastYCSensDate(
    const DateTime& currentLastDate,
    IBadDayAdjusterConstSP bda) const
{
    const static string method("VanillaCreditContingentLeg::lastYCSensDate");
    //for now
    throw ModelException(method,"not yet implemented");
}

CashFlowArraySP VanillaCreditContingentLeg::generateKnownCashFlows(
    const DateTime&        today,
    double                 initialTrancheSize,
    const CashFlowArray&   pastTrancheLosses,
    const BoolArray&       payPastTrancheLosses,
    IBadDayAdjusterConstSP bda) const
{
    const static string method("VanillaCreditContingentLeg::generateKnownCashFlows");
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
void VanillaCreditContingentLeg::getPaymentInformation(BoolArraySP&     payAsYouGoArray,
                                    IntArraySP&      numDelayDaysArray,
                                    DateTimeArraySP& startDate,
                                    DateTimeArraySP& endDate,
                                    DateTimeArraySP& paymentDate) const
{
    const static string method("VanillaCreditContingentLeg::getPaymentInformation");
    //for now
    throw ModelException(method,"not yet implemented");
}

/* Checks if the input date is covered for protection, i.e., falls in
    * one of the observation periods of this leg */
bool VanillaCreditContingentLeg::isDateCoveredForProtection(const DateTime& date) const
{
    const static string method("VanillaCreditContingentLeg::isDateCoveredForProtection");
    //for now
    throw ModelException(method,"not yet implemented");
}


/**
 * Compute PV of fee leg at valuation date.
 */
double VanillaCreditContingentLeg::getContingentLegPV(const DateTime&            valuationDate, 
                                                      const IDiscountCurveRisky& crv,
                                                      IBadDayAdjusterConstSP     bda) const
{
    return this->getContingentLegPV(valuationDate, valuationDate, crv, bda);
}

/**
 * Compute PV of fee leg at valuationDate, but for unconditional payment at 
 * paymentDate, conditional on no new defaults before valuationDate.
 */
double VanillaCreditContingentLeg::getContingentLegPV(const DateTime&            valuationDate, 
                                                      const DateTime&            paymentDate, 
                                                      const IDiscountCurveRisky& crv,
                                                      IBadDayAdjusterConstSP     bda) const //bda is currently unused
{
    const static string method("VanillaCreditContingentLeg::getContingentLegPV");
    
    double sp = 1.0; /* survival probability from valuationDate to paymentDate */
    double df = 1.0; /* IR discount from valuationDate to paymentDate */
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
   else if( valuationDate >= maturityDate) /* i.e. valuationDate >= protectionEndDate*/
   {   /* passed end of protection date so value is zero */
       return 0.0;
   }

   /* Value of protection leg if credit survives until paymentDate */
   survivalPV = crv.protectionPV(paymentDate,
                                 paymentDate.max(effectiveDate),
                                 paymentDate.max(maturityDate),
                                 recovery,//(!useSwapRecovery ? IDiscountCurveRisky::RECOVER_1_MINUS_R : IDiscountCurveRisky::RECOVER_1) ,
                                 delay);
                                 
   /* Protection PV from valuationDate to paymentDate*/
   //Commented out the below as this was also commented out in BasicCreditContingentLeg
   if (valuationDate.max(effectiveDate) <=paymentDate.min(maturityDate))
   {
   defaultPV = crv.protectionPV(valuationDate,
                                valuationDate.max(effectiveDate),
                                paymentDate.min(maturityDate),
                                recovery,//(!useSwapRecovery ? IDiscountCurveRisky::RECOVER_1_MINUS_R : IDiscountCurveRisky::RECOVER_1) ,
                                paymentDate);
   }

   /*  
    ** CredDefSwap has useSwapRecovery; if true, use 1-R, else 1 as the recovery is 
    ** incorporated in the crv when computing the protectionPV.
    ** @TODO: The underlying risky curve can have a fixed constant recovery rate or a
    ** recovery rate term structure.  Hence the survivalPV should be adjusted in 
    ** the contingentLeg by either its own recovery rate or it's already adjusted 
    ** by the risky curve when asking for protectionPV.
    */
   pv = baseNotional *  1 // //(useSwapRecovery ? IDiscountCurveRisky::RECOVER_1_MINUS_R : IDiscountCurveRisky::RECOVER_1) 
        * (sp * survivalPV + defaultPV) * df;
   return pv;
}

/** Return the recovery rate used
    crv should be the underlying if the leg does not support an override */
double VanillaCreditContingentLeg::getRecovery(const IDiscountCurveRisky& crv) const
{
    return recoveryRate;
}

/** Returns the amount that would be recovered upon default */
double VanillaCreditContingentLeg::recoveredValue(const DateTime& valueDate,
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

    double defaultValue = recRate * baseNotional * factor * balance;
    return defaultValue;
}

/** Returns the amount that would be recovered upon default,
    using the supplied recovery rate */
double VanillaCreditContingentLeg::recoveredValue(const DateTime& valueDate,
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
    double defaultValue = recRate * baseNotional * factor * balance;
    return defaultValue;
}

/** Prices the leg sufferring a default, under the (optional) credit event */
double VanillaCreditContingentLeg::getContingentLegDefaultedPV(
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
    const static string method(
        "VanillaCreditContingentLeg::getContingentLegDefaultedPV");
    //for now
    throw ModelException(method,"not yet implemented");
}

/** Returns the leg notional */
double VanillaCreditContingentLeg::getContingentLegNotional() const
{
    return baseNotional;
}

/** Sets the leg notional */
void VanillaCreditContingentLeg::setContingentLegNotional(double newNotional)
{
    baseNotional = newNotional;
}

//--------------------------------------
// ICreditContingentLegGenerator methods
//--------------------------------------

ICreditContingentLegSP VanillaCreditContingentLeg::generateCreditContingentLeg(const DateTime& startDate,
                                                                               const DateTime& endDate) const
{

    const static string method = "VanillaCreditContingentLeg::generateCreditContingentLeg";
    VanillaCreditContingentLegSP newContingentLeg(dynamic_cast<VanillaCreditContingentLeg*>(this->clone()));
	
    if (!newContingentLeg) {
        throw ModelException(method, "Generating contingent leg failed.");
    }

    /* reset the protection start and end date */
    newContingentLeg->effectiveDate = startDate;
    newContingentLeg->maturityDate  = endDate;

    return newContingentLeg;
}

//---------------------
// Theta::Shift methods
//---------------------

/** Implementation of the Theta Shift interface */
bool VanillaCreditContingentLeg::sensShift(Theta* shift) 
{
    static const string method = "VanillaCreditContingentLeg::sensShift";

    return true;
}

//-------------------
// IGetMarket methods
//-------------------

/** Populate from market cache */
void VanillaCreditContingentLeg::getMarket(const IModel* model, const MarketData* market) 
{
    const static string method = "VanillaCreditContingentLeg::getMarket";

    market->GetReferenceDate(valueDate);

    if(!discount.isEmpty())
    {
        discount.getData(model, market);
    }

    if(!pmtHol.isEmpty())
    {
        pmtHol.getData(model, market);
    }
    else
    {
        // default it to weekends only
        pmtHol.setObject(MarketObjectSP(Holiday::weekendsOnly()));
    }
}

//---------------------------------------------
// VanillaCreditContingentLeg methods (private)
//---------------------------------------------

VanillaCreditContingentLeg::VanillaCreditContingentLeg() :CObject(TYPE)
{
    /* A VanillaCreditContingentLeg is always 1-R */
    recoveryTypeString = "RECOVER_1_MINUS_R";
    recovery = IDiscountCurveRisky::RECOVER_1_MINUS_R;
}

VanillaCreditContingentLeg::VanillaCreditContingentLeg(const DateTime& startDate,
                                                       const DateTime& endDate):CObject(TYPE),effectiveDate(startDate),
                                                       maturityDate(endDate)
{
    /* A VanillaCreditContingentLeg is always 1-R */
    recoveryTypeString = "RECOVER_1_MINUS_R";
    recovery = IDiscountCurveRisky::RECOVER_1_MINUS_R;
}

void VanillaCreditContingentLeg::load(CClassSP& clazz){
    clazz->setPublic();
    REGISTER(VanillaCreditContingentLeg, clazz);
    SUPERCLASS(CObject);
    IMPLEMENTS(Theta::Shift);
    IMPLEMENTS(ICreditContingentLeg);
    IMPLEMENTS(ICreditContingentLegGenerator);
    IMPLEMENTS(IGetMarket);
    EMPTY_SHELL_METHOD(defaultConstructor);

    FIELD(baseNotional,               "initial notional ");
    FIELD(effectiveDate,              "protection start date ");
    FIELD(valueDate,                  "valuation date ");
    FIELD(maturityDate,               "maturity date ");
    FIELD(recoveryRate,               "recovery rate");
    FIELD(delay,                      "payment delay");
    FIELD       (pmtBdc,                     "payment bad day convention");
    FIELD(pmtHol,                     "payment holiday convention");
    FIELD(discount,                   "Discount curve");

    FIELD_MAKE_OPTIONAL(baseNotional);
    FIELD_MAKE_OPTIONAL(valueDate);
    FIELD_MAKE_OPTIONAL(delay);
    FIELD_MAKE_OPTIONAL(discount);

    FIELD(recoveryTypeString,         "recovery type");
    FIELD_MAKE_TRANSIENT(recoveryTypeString);
}

IObject* VanillaCreditContingentLeg::defaultConstructor()
{
    return new VanillaCreditContingentLeg();
}

CClassConstSP const VanillaCreditContingentLeg::TYPE = 
CClass::registerClassLoadMethod("VanillaCreditContingentLeg", typeid(VanillaCreditContingentLeg), load);

/* to ensure class is linked in */
bool VanillaCreditContingentLegLoad(){
    return VanillaCreditContingentLeg::TYPE != NULL;
}


//
//
///****************** ICreditVanillaInstrument ******************/
//
///* Will return the default value at maturity date */
//CashFlowArraySP VanillaCreditContingentLeg::getInstrumentCashFlows() const
//{
//    double defaultValue = baseNotional * (1 - recoveryRate);
//
//    CashFlowArraySP cflArray(new CashFlowArray(1));
//    (*cflArray)[0] = CashFlow(maturityDate, defaultValue);    
//
//    return cflArray;
//}
//
///**Accrued interest (cash amount, not percentage) for settlement on settlementDate.*/
//double VanillaCreditContingentLeg::getAccruedInterest(const DateTime& settlementDate) const
//{
//    return 0.0;
//}
//
///** Calculates PV given a valuation date. 
// *  Settlement is whatever the instrument determines it to be (e.g. T+1 CDS; T+3 bond, etc.)
// *  relative to the valuation date. Value should be conditional on no new defaults 
// *  before valuationDate.
// */
//double VanillaCreditContingentLeg::getPV(const DateTime& valuationDate, 
//                                         const IDiscountCurveRisky& crv) const
//{
//    return this->getPV(valuationDate, valuationDate, crv);
//}
//
///*Value for unconditional settlement on paymentDate conditional on survival to valuationDate.*/
//
///** Calculates the PV on a given valuationDate for unconditional settlement
// *  on paymentDate conditional on survival to valuationDate.
// */
//double VanillaCreditContingentLeg::getPV(const DateTime&                 valuationDate, 
//                                         const DateTime&                 paymentDate, 
//                                         const IDiscountCurveRisky&      crv) const
//{
//    const static string method("VanillaCreditContingentLeg::getPV");
//    
//    double sp = 1.0; /* survival probability from valuationDate to paymentDate */
//    double df = 1.0; /* IR discount from valuationDate to paymentDate */
//    double survivalPV = 0.0;
//    double defaultPV  = 0.0;
//    double pv         = 0.0;
//
//   if(valuationDate < paymentDate)
//   {
//       sp = crv.survivalProb(valuationDate, paymentDate);
//       if(sp == 0.0)
//       {
//           throw ModelException (method, "Unexpected zero survival probability.");
//       }
//       df = crv.pv(valuationDate, paymentDate);
//       df /= sp;
//   }
//   else if (valuationDate > paymentDate)
//   {
//       throw ModelException(method, "valuationDate must always be on or before paymentDate.");
//   }
//   else if( valuationDate >= maturityDate) /* i.e. valuationDate >= protectionEndDate*/
//   {   /* passed end of protection date so value is zero */
//       return 0.0;
//   }
//
//   /* Value of protection leg if credit survives until paymentDate */
//   survivalPV = crv.protectionPV(paymentDate,
//                                 paymentDate.max(effectiveDate),
//                                 paymentDate.max(maturityDate),
//                                 recovery,//(!useSwapRecovery ? IDiscountCurveRisky::RECOVER_1_MINUS_R : IDiscountCurveRisky::RECOVER_1) ,
//                                 delay);
//                                 
//   /* Protection PV from valuationDate to paymentDate*/
//   defaultPV = crv.protectionPV(valuationDate,
//                                valuationDate.max(effectiveDate),
//                                paymentDate.min(maturityDate),
//                                recovery,//(!useSwapRecovery ? IDiscountCurveRisky::RECOVER_1_MINUS_R : IDiscountCurveRisky::RECOVER_1) ,
//                                paymentDate);
//
//   /*  
//    ** CredDefSwap has useSwapRecovery; if true, use 1-R, else 1 as the recovery is 
//    ** incorporated in the crv when computing the protectionPV.
//    ** @TODO: The underlying risky curve can have a fixed constant recovery rate or a
//    ** recovery rate term structure.  Hence the survivalPV should be adjusted in 
//    ** the contingentLeg by either its own recovery rate or it's already adjusted 
//    ** by the risky curve when asking for protectionPV.
//    */
//   pv = baseNotional *  1 // //(useSwapRecovery ? IDiscountCurveRisky::RECOVER_1_MINUS_R : IDiscountCurveRisky::RECOVER_1) 
//        * (sp * survivalPV + defaultPV) * df;
//   return pv;
//}
//
///** Calculates the PV of an instrument on defaultDate, given that default has just occured.
// *  Note that the curve is needed in case there is a delay between default and recovery
// *  payment. This method assumes a complete default, so should be interpreted with care
// *  when the instrument/curve has multiple names.
// */
//double VanillaCreditContingentLeg::getPVGivenDefault(const DateTime&                 defaultDate,
//                                                     const IDiscountCurveRisky&      crv) const
//{
//    const static string method = "VanillaCreditContingentLeg::getPVGivenDefault";
//    DateTime adjDefaultDate;
//
//    /* default date is in the future */
//    if(defaultDate > valueDate)
//    {
//        throw ModelException(method, "Default date is in the future " 
//                             "(Default date is " + defaultDate.toString() +
//                             ", Value date is " + valueDate.toString() + ")");
//    }
//
//    /* This assumes recovery rate is from the curve */
//    /* TODO: incorpate getting the recovery rate from the risky curve */
//    //return baseNotional*(1-crv.getRecovery(defaultDate));
//
//    return baseNotional * (1 - recoveryRate);
//}
//
///* Is not perpetual */
//bool VanillaCreditContingentLeg::hasFiniteMaturity() const
//{
//    return true;
//}
//
///* Returns the last accrual date for the fee payments. */
//DateTime VanillaCreditContingentLeg::getMaturity() const
//{
//    return maturityDate;
//}
//
///**
// *  Returns the earliest possible "interesting" date: this will be the earliest of the
// *  start of the first accrual period, the start of the contingent leg, the first
// *  cash-flow, etc.
//*/
//DateTime VanillaCreditContingentLeg::getStartDate() const
//{
//    return effectiveDate;
//}
//
///** 
// *  Return the notional of the instrument. If the instrument is amortizing, this
// *  will be the base notional. The idea is that getPV(...)/getNotional() should be
// *  the price of the instrument, in some meaningful sense.
// */
//double VanillaCreditContingentLeg::getNotional() const
//{
//    return baseNotional;
//}
//
//
///************* End ICreditVanillaInstrument ****************************/
//
//
//
//
//
//
//
///********* IInstrument **************************************/
//void VanillaCreditContingentLeg::GetMarket(const IModel* model, const CMarketDataSP market)
//{
//    this->getMarket(model, market.get());
//}
//
//
///** Called once before the initial pricing */
//void VanillaCreditContingentLeg::Validate()
//{
//    static const string method = "VanillaCreditContingentLeg::validate";
//
//    /* A VanillaCreditContingentLeg will always be 1-R */
//    /* Should this be set in the default constructor?  */
//    recoveryTypeString == "RECOVER_1_MINUS_R";
//    recovery = IDiscountCurveRisky::RECOVER_1_MINUS_R;
//
//}
//
//CSensControl* VanillaCreditContingentLeg::AlterControl(const IModel*          modelParams,
//                                                       const CSensControl*    sensControl) const
//{      
//    return 0;
//}
//
//
///** Returns the value date (aka today) the instrument is currently
//    pricing for */
//DateTime VanillaCreditContingentLeg::getValueDate() const
//{
//    return valueDate;
//}
//
//bool VanillaCreditContingentLeg::priceDeadInstrument(CControl* control,
//                                                     CResults* results) const
//{
//    return false;
//}
//
//
//string VanillaCreditContingentLeg::discountYieldCurveName() const
//{
//    return discount.getName();
//}
//

DRLIB_END_NAMESPACE
