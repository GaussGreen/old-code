//----------------------------------------------------------------------------
//
//   Group       : Convertibles DR
//
//   Filename    : CredDefSwap.cpp
//
//   Description : Credit default swap
//
//   Author      : Tycho von Rosenvinge
//
//   Date        : August 8, 2001
//
//----------------------------------------------------------------------------
#include "edginc/config.hpp"
#include "math.h"
#include "edginc/Maths.hpp"
#include "edginc/Class.hpp"
#include "edginc/Format.hpp"
#include "edginc/ExpiryResult.hpp"
#include "edginc/CredDefSwap.hpp"
#include "edginc/Nrfns.hpp"
#include "edginc/ZeroCurve.hpp"
#include "edginc/Actual365F.hpp"
#include "edginc/CompoundBasis.hpp"
#include "edginc/OutputRequestUtil.hpp"
#include "edginc/RateConversion.hpp"
#include "edginc/Delta.hpp"
#include "edginc/DayCountConventionFactory.hpp"
#include "edginc/VegaMatrix.hpp"
#include "edginc/RootTimeVega.hpp"
#include "edginc/PowerVega.hpp"
#include "edginc/VegaSkewParallel.hpp"
#include "edginc/VegaSkewPointwise.hpp"
#include "edginc/BadDayConventionFactory.hpp"
#include "edginc/Addin.hpp"
#include "edginc/DoubleMatrix.hpp"
#include "edginc/BenchmarkDate.hpp"
#include "edginc/imsl.h"
#include "edginc/CDSHelper.hpp"
#include "edginc/CDSPricer.hpp"
#include "edginc/SwapTool.hpp"
#include "edginc/StubFactory.hpp"
#include "edginc/Actual365F.hpp"
#include "edginc/CDSCreditSupport.hpp"
#include "edginc/SensMgr.hpp"
#include "edginc/IRestorableWithRespectTo.hpp"
#include "edginc/PropertyTweakHypothesis.hpp"
#include "edginc/RiskProperty.hpp"
#include "edginc/Spot.hpp"
#include "edginc/IDiscountCurveRisky.hpp"
#include "edginc/ICreditEventOverrideName.hpp"
#include "edginc/DetailedCreditEventOverrideName.hpp"
#include "edginc/FeeLegReductionPerDefault.hpp"
#include "edginc/XMLWriter.hpp"
#include "edginc/CtgLegLossPerDefault.hpp"

DRLIB_BEGIN_NAMESPACE

// Allows for margin of error in computing max tweak sizes
const double CredDefSwap::CDS_ADJUSTMENT = .75;

static CFieldConstSP cdsParSpreadsField;
static CFieldConstSP assetField;

//--------------------
// CredDefSwap methods
//--------------------

CredDefSwap::~CredDefSwap(){}

//--------------
// ICDS methods
//--------------

/** Return the fee leg */
ICreditFeeLegSP CredDefSwap::getFeeLeg() const
{
    const static string method = "CredDefSwap::getFeeLeg";
    //we allow the CredDefSwap to behave as a fee leg
    //due to the way the fee inputs are defined on the instrument
    return ICreditFeeLegSP(dynamic_cast<ICreditFeeLeg*>(const_cast<CredDefSwap*>(this)));
}

/** Return the contingent leg */
ICreditContingentLegSP CredDefSwap::getContingentLeg() const
{
    const static string method = "CredDefSwap::getContingentLeg";
    //we allow the CredDefSwap to behave as a ctg leg
    //due to the way the fee inputs are defined on the instrument
    return ICreditContingentLegSP(dynamic_cast<ICreditContingentLeg*>(const_cast<CredDefSwap*>(this)));
}

//----------------------
// ICreditFeeLeg methods
//----------------------

/** Return the value of this leg */
double CredDefSwap::price(
    const DateTime&             today,
    const DateTime&             valDateCF,
    const IDiscountCurveRiskySP effectiveCurve,
    const YieldCurveWrapper&    discount,
    const double&               lowStrike,
    const double&               highStrike,
    const double&               outstandingNotional,
    const CashFlowArray&        pastTrancheLosses,
    double&                     riskyDurationTotal,  // (O) notional weighted risky duration
    double&                     riskyNotionalsMean,  // (O) mean value of notional per period
    double&                     risklessCFPV,        // (O) fair value of riskless payments
    bool                        computeExtra,        // true: populate arrays
    DoubleArray&                debugUnitPrice,      // price for each leg unit
    DoubleArray&                debugUnitHistPrice,  // price for each leg unit due to historical default
    CashFlowArraySP             rebatePayments,
    BoolArrayConstSP            payAsYouGoArray,
    IntArrayConstSP             numDelayDaysArray,
    DateTimeArrayConstSP        startDates,
    DateTimeArrayConstSP        endDates,
    DateTimeArrayConstSP        paymentDates,
    IBadDayAdjusterConstSP      bda,
    IForwardRatePricerSP model) const
{
    static const string method = "CredDefSwap::price";
    //for now
    throw ModelException(method,"not yet implemented");
}

/** Return all cash flow dates */
DateTimeArraySP CredDefSwap::getCashFlowDates() const
{
    int numFees = feePayments->size();
    DateTimeArraySP cfDates = DateTimeArraySP(
        new DateTimeArray(numFees));

    for (int i=0; i<numFees; i++)
    {
        (*cfDates)[i] = ((*feePayments)[i]).date;
    }

    return cfDates;
}

/** Return risky cash flow dates */
DateTimeArraySP CredDefSwap::getRiskyCashFlowDates() const
{
    //all fees are risky in this structure
    return getCashFlowDates();
}

/** Return riskfree cash flow dates */
DateTimeArraySP CredDefSwap::getRisklessCashFlowDates() const
{
    //no riskless fees in this structure
    return DateTimeArraySP();
}

/** Returns the accrual periods (start and end dates) in this fee leg. */
AccrualPeriodArrayConstSP CredDefSwap::getAccrualPeriods() const
{
    static const string method = "CredDefSwap::getAccrualPeriods";
    //requires backing out from the fee amount
    //so for now....
    throw ModelException(method,"not yet implemented");
}

/** Returns risky accrual periods (start and end dates) in this fee leg. */
AccrualPeriodArrayConstSP CredDefSwap::getRiskyAccrualPeriods() const
{
    //all fees are risky in this structure
    return getAccrualPeriods();
}

/** Returns risky accrual periods (start and end dates) in this fee leg. */
AccrualPeriodArrayConstSP CredDefSwap::getRisklessAccrualPeriods() const
{
    //no riskless fees in this structure
    return AccrualPeriodArrayConstSP();
}

/** Returns risky coupon notional types */
CouponNotionalTypesArraySP CredDefSwap::getRiskyCouponNotionalTypes() const
{
    static const string method = "CredDefSwap::getRiskyCouponNotionalTypes";

    //for now
    throw ModelException(method, "Not yet implemented");
}

/** Returns risky observation dates */
DateTimeArraySP CredDefSwap::getRiskyObservationDates() const
{
    static const string method = "CredDefSwap::getRiskyObservationDates";

    //for now
    throw ModelException(method, "Not yet implemented");
}

/** Returns all cash flows */
AbstractCashFlowArrayConstSP CredDefSwap::getCashFlows(IForwardRatePricerSP model) const
{
    static const string method = "CredDefSwap::getCashFlows";
    //for now
    throw ModelException(method,"not yet implemented");
}

/** Returns risky cash flows only */
CashFlowArraySP CredDefSwap::getRiskyCashFlows(IForwardRatePricerSP model) const
{
    static const string method = "CredDefSwap::getRiskyCashFlows";
    //for now
    throw ModelException(method,"not yet implemented");
}

/** Returns risk free cash flows only */
CashFlowArraySP CredDefSwap::getRisklessCashFlows(IForwardRatePricerSP model) const
{
    static const string method = "CredDefSwap::getRisklessCashFlows";
    //for now
    throw ModelException(method,"not yet implemented");
}

/** Returns risky notional dates */
DateTimeArraySP CredDefSwap::getRiskyNotionalDates(IForwardRatePricerSP model) const
{
    static const string method = "CredDefSwap::getRiskyNotionalDates";
    //for now
    throw ModelException(method,"not yet implemented");
}

/** Return known cash flows corresponding to a CDO tranche */
CashFlowArraySP CredDefSwap::generateKnownCashFlows(
        const DateTime       today,
        const double         initialTrancheSize,
        const DateTimeArray  pastTrancheLossDates,
        const DoubleArray    pastTrancheLosses,
        const double         pastTrancheLoss,
        IForwardRatePricerSP model)
{
    static const string method = "CredDefSwap::generateKnownCashFlows";
    //for now
    throw ModelException(method,"not yet implemented");
}

/** Return known cash flows corresponding to a non-defaulted CDS */
CashFlowArraySP CredDefSwap::generateKnownCashFlows(IForwardRatePricerSP model) const
{
    static const string method = "CredDefSwap::generateKnownCashFlows";
    //for now
    throw ModelException(method,"not yet implemented");
}

/** Estimates the known cash flows (so they are not really "known")
    * corresponding to a CDO tranche - takes into account estimated losses
    * in the future */
CashFlowArraySP CredDefSwap::estimateKnownCashFlows(
        const double         initialTrancheSize,
        const DateTimeArray  pastTrancheLossDates,
        const DoubleArray    pastTrancheLosses,
        IForwardRatePricerSP model)
{
    static const string method = "CredDefSwap::estimateKnownCashFlows";
    //for now
    throw ModelException(method,"not yet implemented");
}

/** Returns the known cash flows corresponding to a defaulted CDS
    * (taking accrued payments into consideration), that happen after a specific
    * date (typically used to exclude cashflows in the past, already paid).
    * If excludePaymentsBeforeDate is empty, all cashflows will be returned.
    *
    * CAUTION:
    * It does not aggregate cashflows happening on the same dates into one.
    * If this is required (e.g., to output these cashflows for MiddleOffice use)
    * you should call CashFlow::agregate(...) on the resulting CashFlowArray -
    * It is not done here for performance reasons */
CashFlowArraySP CredDefSwap::generateKnownCashFlowsGivenDefault(
    const DateTime&            valueDate,
    const DateTime&            defaultDate,
    const DateTime&            excludePaymentsBeforeDate,
    const DateTime&            protectionStartDate,
    const DateTime&            protectionEndDate,
    const bool                 allowIncludingTodaysPayments,
    IForwardRatePricerSP       model,
    ICreditEventOverrideNameSP creditEventOverride,
    IBadDayAdjusterConstSP     badDayAdjuster,
    CIntSP                     triggerDelay,
    CIntSP                     defaultToSettlementDelay,
    const DateTime&            lastTriggerDate) const
{
    static const string method = "CredDefSwap::generateKnownCashFlowsDefaulted";
    //for now
    throw ModelException(method, "not yet implemented");
}

/** Returns the pay date which is last in terms of time */
DateTime CredDefSwap::getLastPayDate() const
{
    CashFlow last = feePayments->back();
    return last.date;
}

/** Returns the observation date which is last in terms of time */
DateTime CredDefSwap::getLastObservationDate() const
{
    static const string method = "CredDefSwap::getLastObservationDate";
    //for now
    throw ModelException(method,"not yet implemented");
}

/** When to stop tweaking for Yield Curve type tweaks */
DateTime CredDefSwap::lastYCSensDate(const DateTime& currentLastDate) const
{
    static const string method = "CredDefSwap::lastYCSensDate";
    //for now
    throw ModelException(method,"not yet implemented");
}

/** Feedback method for getting information about a fee cashflow */
void CredDefSwap::getActiveFee(const DateTime&      withRespectTo,       // (I) get the fee whose accrual period contains this date
                               const DateTime&      earliestAccrualDate, // (I) for fee legs that dont specify accrue start dates, and we are interested in the first fee
                               IForwardRatePricerSP model,               // (I) for calculating the amount
                               DateTime&            accrueStartDate,     // (O) when the fee starts accruing
                               DateTime&            accrueEndDate,       // (O) when the fee finishes accruing
                               DateTime&            paymentDate,         // (O) when the fee is paid
                               double&              amount) const        // (O) the cashflow amount
{
    int numFees = feePayments->size();
    int idx;

    for (idx=0; idx<numFees; idx++)
    {
        //accrual period treated as payment to payment
        accrueStartDate = idx==0 ? earliestAccrualDate : (*feePayments)[idx-1].date;
        accrueEndDate   = (*feePayments)[idx].date;

        if (withRespectTo.within(accrueStartDate,accrueEndDate))
        {
            break;
        }
    }

    if (idx == numFees)
    {
        DateTime empty;
        accrueStartDate = empty;
        accrueEndDate   = empty;
        paymentDate     = empty;
        amount = 0.0;
    }
    else
    {
        paymentDate     = accrueEndDate;
        amount          = (*feePayments)[idx].amount;
    }
}

/** Price this fee leg assuming it corresponds to a defaulted CDS */
double CredDefSwap::pv(const DateTime&         valueDate,
                    const DateTime&         defaultDeterminationDate,
                    const DateTime&         accrualPaymentDate,
                    const YieldCurveConstSP discount,
                    const bool              computeAccrual,
                      IForwardRatePricerSP    model) const
{
    static const string method = "CredDefSwap::pv";
    //for now
    throw ModelException(method,"not yet implemented");
}

/**
* Compute PV of fee leg at valuationDate, but for unconditional payment at
* paymentDate, conditional on no new defaults before valuationDate.
*/
double CredDefSwap::getFeeLegPV(const DateTime&              valuationDate,
                                const DateTime&              earliestRiskyDate,
                                const DateTime&              paymentDate,
                                const DateTime&              latestRiskyDate,
                                const IDiscountCurve&        discount,
                                const IDiscountCurveRisky&   crv,
                                const IDecretionCurveConstSP prepay,
                                const bool                   includeAccrued,
                                const DayCountConventionSP   dcc,
                                IForwardRatePricerSP         model) const
{
    static const string method = "CredDefSwap::getFeeLegPV";
    //for now
    throw ModelException(method,"not yet implemented");
}

/** Returns the accrued interest */
double CredDefSwap::getFeeLegAI(const DateTime&              valuationDate, 
                                const DateTime&              paymentDate, 
                                const DateTime&              earliestAccrualStart,
                                const DateTime&              latestAccrualEnd,
                                const DayCountConventionSP   dcc, //allows an override to be specified
                                const IDiscountCurveRisky&   crv,
                                const IDiscountCurveConstSP  discount,
                                const IDecretionCurveConstSP prepay,
                                IForwardRatePricerSP         model) const
{
    const static string method = "VanillaCreditFeeLeg::getFeeLegAI";
    throw ModelException(method,"not yet implemented");
}

/** Prices the leg sufferring a default, under the (optional) credit event */
double CredDefSwap::getFeeLegDefaultedPV(const DateTime&              valuationDate,
                                         const DateTime&              defaultDate,
                                         const DateTime&              protectionStartDate,
                                         const DateTime&              protectionEndDate,
                                         const bool                   allowIncludingTodaysPayments,
                                         const IDiscountCurveConstSP  discount,
                                         const IDecretionCurveConstSP prepay,
                                         IForwardRatePricerSP         model,
                                         IBadDayAdjusterConstSP       badDayAdjuster,
                                         ICreditEventOverrideNameSP   creditEventOverride,
                                         CIntSP                       triggerDelay,
                                         CIntSP                       defaultToSettlementDelay,
                                         DateTime                     lastTriggerDate) const
{
    const static string method = "VanillaCreditFeeLeg::getFeeLegDefaultedPV";

    //should resolve to priceFeeLegGivenDefault
    throw ModelException(method,"not yet implemented");
}

/** Returns the leg notional */
double CredDefSwap::getFeeLegNotional() const
{
    return notional;
}

/** Returns the leg notional */
void CredDefSwap::setFeeLegNotional(double newNotional)
{
    notional = newNotional;
}

ICreditFeeLegSVGenSP
CredDefSwap::createSVGen(
	const DateTime& valueDate,
	const DateTimeArray& modifiedFeeLegObservationDates,					 
	const DateTimeArray& productTimeline, //product timeline
	const IntArray& dateToDiscFactorIndex, //map from a date to an index on the DiscountFactor SV
	double lossConfigNotional
	) const
{
	throw ModelException("VanillaCreditFeeLeg::createSVGen","Method not implemented");
}


//-------------------------------
// IFixedRateCreditFeeLeg methods
//-------------------------------
/* Note this is ugly, and will lead to nasty cumulative rounding errors if you change the
    often, but is all I can do to keep backwards compatibility. This is also !VERY BAD! if you
   defaulted the feeRate to 1.0. */

void CredDefSwap::setRate(double newFee) {
    const static string method = "CredDefSwap::setRate";
    if (newFee==feeRate) {
        return;
    } else if (newFee==0.0 && feeRate!=0.0) {
        throw ModelException(method,
            "The newFee must be non-zero unless the old one is zero, too!");
    }
    double feeRatio = newFee/feeRate;
    //assign new rate
    feeRate = newFee;
    int i=0;
    for (i=0; i<feePayments->size(); i++) {
        //scale input cashflows by the ratio of rates
        (*feePayments)[i].amount *= feeRatio;
    }
    return;
}

double CredDefSwap::getRate() const {
    return feeRate;
}

//-----------------------------
// ICreditContingentLeg methods
//-----------------------------

/** Price this contingent leg. */
double CredDefSwap::price(
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
    static const string method = "CredDefSwap::price";
    //for now
    throw ModelException(method,"not yet implemented");
}

/** Returns the earliest observation start date */
DateTime CredDefSwap::firstObservationStartDate() const
{
    return swapEffectiveDate;
}

/** Returns the last pay date */
DateTime CredDefSwap::lastPayDate(IBadDayAdjusterConstSP bda) const
{
    CashFlow last = feePayments->back();
    return last.date;
}

/** Returns the last observation date */
DateTime CredDefSwap::lastObservationEndDate() const
{
    return protectionEndDate;
}

/** Returns the leg notional */
double CredDefSwap::getContingentLegNotional() const
{
    return notional;
}

/** Sets the leg notional */
// actually applies to both legs
void CredDefSwap::setContingentLegNotional(double newNotional)
{
    notional = newNotional;
}

/** When to stop tweaking */
DateTime CredDefSwap::lastYCSensDate(const DateTime& currentLastDate,
                                     IBadDayAdjusterConstSP bda) const
{
    static const string method = "CredDefSwap::lastYCSensDate";
    //for now
    throw ModelException(method,"not yet implemented");
}

CashFlowArraySP CredDefSwap::generateKnownCashFlows(
    const DateTime&        today,
    double                 initialTrancheSize,
    const CashFlowArray&   pastTrancheLosses,
    const BoolArray&       payPastTrancheLosses,
    IBadDayAdjusterConstSP bda) const
{
    static const string method = "CredDefSwap::generateKnownCashFlows";
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
void CredDefSwap::getPaymentInformation(BoolArraySP&     payAsYouGoArray,
                                        IntArraySP&      numDelayDaysArray,
                                        DateTimeArraySP& startDate,
                                        DateTimeArraySP& endDate,
                                        DateTimeArraySP& paymentDate) const
{
    static const string method = "CredDefSwap::getPaymentInformation";
    //for now
    throw ModelException(method,"not yet implemented");
}
/* Checks if the input date is covered for protection, i.e., falls in
    * one of the observation periods of this leg */
bool CredDefSwap::isDateCoveredForProtection(const DateTime& date) const
{
    static const string method = "CredDefSwap::isDateCoveredForProtection";
    //for now
    throw ModelException(method,"not yet implemented");
}

///**Compute PV of contingent leg at valuation date, with instrument default settlement and
//    payment date behaviour.*/
//double CredDefSwap::getContingentLegPV(const DateTime&                 valuationDate,
//                                       const IDiscountCurveRisky&      crv) const
//{
//    static const string method = "CredDefSwap::getContingentLegPV";
//    //for now
//    throw ModelException(method,"not yet implemented");
//}
//
///**Compute PV of contingent leg at valuationDate, but for unconditional payment at
//    paymentDate, conditional on no new defaults before valuationDate.*/
//double CredDefSwap::getContingentLegPV(const DateTime&                 valuationDate,
//                                       const DateTime&                 paymentDate,
//                                       const IDiscountCurveRisky&      crv) const
//{
//    static const string method = "CredDefSwap::getContingentLegPV";
//    //for now
//    throw ModelException(method,"not yet implemented");
//}


YieldCurveWrapper CredDefSwap::getYieldCurveWrapper() const {
    return this->discount;
}
ICDSParSpreadsWrapper CredDefSwap::getParSpreadsWrapper() const {
    return this->cdsParSpreads;
}

double CredDefSwap::getAccruedInterest(const DateTime& settlementDate,
                                       IForwardRatePricerSP model) const 
{
    int idx=0;
    return getAccruedInterest(settlementDate, model, idx);
}

double CredDefSwap::getAccruedInterest(const DateTime& settlementDate,
                                       IForwardRatePricerSP model,
                                       int& feeIdx) const 
{
    return priceFeeLegGivenDefault(settlementDate, 
                                   settlementDate, 
                                   feeIdx, 
                                   model, 
                                   true); // pricing accrued here
}

/* Value for unconditional settlement on paymentDate conditional on survival 
   to valuationDate. */
double CredDefSwap::getContingentLegPV(
    const DateTime&            valuationDate,
    const DateTime&            paymentDate,
    const IDiscountCurveRisky& crv,
    IBadDayAdjusterConstSP     bda) const { //bda is currently unused

    const static string method("CredDefSwap::getContingentLegPV");

    /*sp is the survival probability from valuationDate to paymentDate.*/
    double sp = 1.0;
    /* pure IR discount from valuationDate to paymentDate */
    double df = 1.0;
    if (valuationDate<paymentDate) {
        sp = crv.survivalProb(valuationDate, paymentDate);
        if (sp==0.0) {
            throw ModelException(method, "Unexpected zero survival probability.");
        }
        df = crv.pv(valuationDate, paymentDate);
        df /= sp;
    } else if (valuationDate>paymentDate) {
        throw ModelException(method, "valuationDate must always be on or before paymentDate.");
    } else if (valuationDate>=getProtectionEndDate()) {
        // passed end of protection date - value is zero
        return 0.0;
    }

    /*Value of protection leg if credit survives until paymentDate */
    double survivalPV = crv.protectionPV(paymentDate,
                        paymentDate.max(swapEffectiveDate),
                        paymentDate.max(getProtectionEndDate()),
                        (!useSwapRecovery ? IDiscountCurveRisky::RECOVER_1_MINUS_R : IDiscountCurveRisky::RECOVER_1),
                        0);
    /*protection PV from valuationDate to paymentDate: Value of (1-R) payment at fixed date paymentDate */
    //This seems a bit bogus & leads to inconsistencies when
    //pricing a forward starting trade, and seems only to impact
    //cds options, and therefore this has been commented out to
    //achieve consitency for the options.
    //double defaultPV = crv.protectionPV(
    //                    valuationDate,
    //                    valuationDate.max(swapEffectiveDate),
    //                    paymentDate.min(getProtectionEndDate()),
    //                    (!useSwapRecovery ? IDiscountCurveRisky::RECOVER_1_MINUS_R : IDiscountCurveRisky::RECOVER_1),
    //                    paymentDate);
    //instead we just remove the effect
    double defaultPV = 0.0;

    return notional
        *(useSwapRecovery ? 1.0 - swapRecovery : 1.0)
        *(sp*survivalPV + defaultPV/df);
}


double CredDefSwap::getContingentLegPV(
    const DateTime&            valuationDate,
    const IDiscountCurveRisky& crv,
    IBadDayAdjusterConstSP     bda) const 
{
    return getContingentLegPV(valuationDate,valuationDate, crv, bda);
}

/** Return the recovery rate used
    crv should be the underlying if the leg does not support an override */
double CredDefSwap::getRecovery(const IDiscountCurveRisky& crv) const
{
    double recovery;

    if (!useSwapRecovery)
    {
        recovery = crv.getRecovery();
    }
    else
    {
        recovery = swapRecovery;
    }

    return recovery;
}

/** Returns the amount that would be recovered upon default */
double CredDefSwap::recoveredValue(const DateTime& valueDate,
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
double CredDefSwap::recoveredValue(
    const DateTime& valueDate,
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
double CredDefSwap::getContingentLegDefaultedPV(
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
    const static string method("CredDefSwap::getContingentLegDefaultedPV");
    //should resolve to priceContingentLegGivenDefault
    //for now
    throw ModelException(method,"not yet implemented");
}

/*Value for unconditional settlement on paymentDate conditional on survival to valuationDate.
  For now, I have ignored the coupon recovery possible in case of default before settlement
  and just assumed that the payment contingent on default between valuationDate and paymentDate
  is just zero.*/
double CredDefSwap::getFeeLegPV(
        const DateTime&              valuationDate,
        const DateTime&              paymentDate, 
        const DateTime&              earliestRiskyDate,
        const DateTime&              latestRiskyDate,
        const IDiscountCurve&        discount,
        const IDiscountCurveRisky&   crv,
        const IDecretionCurveConstSP prepay,
        const bool                   includeAccrued,
        const DayCountConventionSP   dcc,
        bool                         defaultValueOnly,
        IForwardRatePricerSP         model) const {

    const static string method("CredDefSwap::getFeeLegPV");

    /* If default value only and no accrued recovery, then no value at all*/
    if (defaultValueOnly && !payAccruedFee) return 0.0;

    /*sp is the survival probability from valuationDate to paymentDate.*/
    double sp = 1.0;
    /* pure IR discount from valuationDate to paymentDate */
    double df = 1.0;
    if (valuationDate<paymentDate) {
        DateTime toRiskyDate = (paymentDate > latestRiskyDate) ? latestRiskyDate : paymentDate;
        sp = crv.survivalProb(valuationDate, toRiskyDate);
        if (sp==0.0) {
            throw ModelException(method, "Unexpected zero survival probability.");
        }
        df = discount.pv(valuationDate, paymentDate);
    } else if (valuationDate>paymentDate) {
        throw ModelException(method, "valuationDate must always be on or before paymentDate.");
    } else if (valuationDate>=getProtectionEndDate()) {
        // passed end of protection date - value is zero
        return 0.0;
    }

    /*Value for payment at paymentDate contingent on survival until then.*/
    double survivalPV = crv.annuityPV(
        *feePayments,
        paymentDate,
        (payAccruedFee ? IDiscountCurveRisky::RECOVER_1 : IDiscountCurveRisky::RECOVER_0),
        0,getStartDate()) -
        (defaultValueOnly ?
            crv.annuityPV(
            *feePayments,
            paymentDate,
            IDiscountCurveRisky::RECOVER_0,
            0,getStartDate()) : 0.0);

    // only get survivalPV if survives until paymentDate, so * by survival prob.
    return sp*survivalPV;
}
double CredDefSwap::getFeeLegPV(
        const DateTime& valuationDate,
        const DateTime& earliestRiskyDate,
        const DateTime& latestRiskyDate,
        const IDiscountCurve& discount,
        const IDiscountCurveRisky& crv,
        const IDecretionCurveConstSP prepay,
        const bool includeAccrued,
        const DayCountConventionSP dcc,
        IForwardRatePricerSP model) const {
            /*This is just fee PV for natural settlement for the instrument. */
    return getFeeLegPV(valuationDate,valuationDate,
        earliestRiskyDate, latestRiskyDate,
        discount, crv, prepay, 
        includeAccrued, dcc, model);
}

double CredDefSwap::getPVGivenDefault(
    const DateTime& valuationDate, 
    const DateTime& defaultDate,
    const IDiscountCurveRisky& crv,
    IForwardRatePricerSP model) const 
{
    /* There is no payment delay for CredDefSwap, so this is easy. 
       Recovery of AI is always 100%.
       Positive notional represents long protection; the fee will be negative, 
       so the accrued interest will be negative, too. */

    if ((defaultDate < swapEffectiveDate) ||
        (defaultDate > getProtectionEndDate()))
    {
        /* JLH This is wrong if there are riskless fees in the fee leg!
         * Should somehow call the getPriceIfDefaulted method, all the logic
         * is there */

        //outside protection period
        return 0.0;
    }
    else {
        /* JLH Again, this is wrong if there are riskless fees in the fee leg.
         * The fee leg value is the accrued interest plus any potential 
         * riskless fees... the logic is in getPriceIfDefaulted. */

        double ctgVal = notional * (1 - crv.getRecovery(defaultDate));
        double feeVal = payAccruedFee ? getAccruedInterest(defaultDate, model) : 0.0;

        double df = crv.risklessPV(valuationDate, defaultDate);
        return (df * ctgVal) + feeVal;
    }
}

DayCountConventionSP CredDefSwap::getAccrualDcc() const
{
    return swpAccrualDCC;
}

CashFlowArraySP CredDefSwap::getInstrumentCashFlows(
        IForwardRatePricerSP            model) const
{
    return knownCashFlows();
}

/*This is pretty crude, and not really satisfactory - should do for initial testing until
  we have a better CDS object.
  TODO: this is causing crashes that I can't diagnose fince my f***ing debugger will not work! I have not
  fixed this problem! [Charles Morcom February 18, 2006]*/
ICDSSP CredDefSwap::generateCDS(
        const DateTime& startDate,
        const DateTime& endDate,
        double newFeeRate) const {
    const static string method = "CredDefSwap::regenerateWithDifferentDates";

    if (startDate>=endDate) throw ModelException(method,"startDate must be before endDate.");

    // deep copy everything
    CredDefSwapSP newCDS(dynamic_cast<CredDefSwap*>(this->clone()));
    if (!newCDS) {
        throw ModelException(method, "Clone failed for this.");
    }

    CashFlowArray& newCFL = *(newCDS->feePayments);
    // extend the CDS so that the first cash-flow is zero and the date is
    // the accrual start date - makes extending coupons much easier, and
    // shouldn't affect value calculations
    if (newCFL.front().date>newCDS->swapEffectiveDate) {
        newCFL.insert(newCFL.begin(), CashFlow(newCDS->swapEffectiveDate, 0.0));
    }
    /* Try to guess the coupon frequency based on the number of days
       between middle payments (to avoid start/end stubs */
    int p0;
    int p1;
    if (newCFL.size()==2) {
        p0 = 0;
        p1 = 1;
    } else {
        p0 = 1;
        p1 = 2;
    }
    // for now, just use same number of days to determine new coupon dates
    const DateTime::Interval dayDiff = newCFL[p1].date.subtract(newCFL[p0].date);
    const DateTime::Interval negDayDiff = newCFL[p0].date.subtract(newCFL[p1].date);

    /*=========================================================================
     * EXTEND COUPONS AND DATES IF endDate>lastdt or startDate<firstdt
     *=======================================================================*/
    // extend maturity end out as far as endDate
    /* Extend coupons by making them the same number of days long. This
       is not quite correct, but with no frequency or coupon information
       it's difficult to see what else to do */
    while (newCFL.end()->date<endDate) {
        const CashFlow& lastCF = newCFL.back();
        const CashFlow& lastButOneCF = newCFL[newCFL.size()-2];
        if (lastCF.date<=lastButOneCF.date) {
            throw ModelException(method,"Cash-Flow dates are not correctly ordered at back.");
        }

        // next coupon date is lesser of endDate or next full coupon
        DateTime nextCouponDate = lastCF.date.add(dayDiff);
        if (nextCouponDate>endDate) {
            nextCouponDate = endDate;
        }
        double nextCouponAmount = lastCF.amount *
            swpAccrualDCC->years(lastCF.date,nextCouponDate)/swpAccrualDCC->years(lastButOneCF.date,lastCF.date);
        CashFlow newCashFlow(nextCouponDate, nextCouponAmount);
        newCFL.push_back(newCashFlow);
    }
    // extend start end back to startDate
    while (newCFL.begin()->date>startDate) {
        const CashFlow& nextCF = newCFL.front();
        const CashFlow& nextButOneCF = newCFL[1];
        if (nextCF.date>=nextButOneCF.date) {
            throw ModelException(method,"Cash-Flow dates are not correctly ordered at front.");
        }

        DateTime newCouponDate = nextCF.date.add(negDayDiff);
        if (newCouponDate<startDate) {
            newCouponDate = startDate;
        }
        double newCouponAmount = nextButOneCF.amount *
            swpAccrualDCC->years(newCouponDate,nextCF.date)/swpAccrualDCC->years(nextCF.date,nextButOneCF.date);
        newCFL.front().amount = newCouponAmount;
        CashFlow newCashFlow(newCouponDate, 0.0);
        // this is not very efficient
        newCFL.insert(newCFL.begin(),newCashFlow);
    }

    /*=========================================================================
     * TRUNCATE IF endDate<lastdt or startDate>firstdt
     *=======================================================================*/
    while (endDate < newCFL.end()->date) {
        if (endDate<=newCFL[newCFL.size()-2].date) {
            // remove whole coupon
            newCFL.pop_back();
        } else {
            // remove part of coupon
            DateTime& lastDate = newCFL[newCFL.size()-2].date;
            double newCouponAmount = newCFL.end()->amount *
            swpAccrualDCC->years(lastDate, endDate)/swpAccrualDCC->years(lastDate,newCFL.end()->date);
            newCFL.end()->date = endDate;
            newCFL.end()->amount = newCouponAmount;
        }
    }
    while (startDate>newCFL.begin()->date) {
        if (startDate>=newCFL[1].date) {
            // remove whole coupon
            newCFL.erase(newCFL.begin());
            newCFL.front().amount = 0.0;
        } else {
            // remove part of coupon
            DateTime& nextDate = newCFL[1].date;
            double newCouponAmount = newCFL[1].amount *
            swpAccrualDCC->years(startDate, nextDate)/swpAccrualDCC->years(newCFL.front().date, nextDate);
            newCFL.front().date = startDate;
            newCFL[1].amount = newCouponAmount;
        }
    }
    newCDS->swapEffectiveDate = startDate;
    newCDS->protectionEndDate = endDate;
    newCDS->createE2Csensitivities = false;
    newCDS->setRate(newFeeRate);

    return newCDS;
}
/**Returns earlier of first accrual date, or protection start date*/
DateTime CredDefSwap::getStartDate() const {
    return swapEffectiveDate;
}
bool CredDefSwap::hasFiniteMaturity() const {
    return true;
}
/* Return last fee payment date, or protection end date if there is none */
DateTime CredDefSwap::getMaturity() const {
    if (!feePayments || feePayments->size()<=0) {
        return protectionEndDate;
    } else {
        return feePayments->back().date;
    }
}
double CredDefSwap::getNotional() const {
    return notional;
}
void CredDefSwap::setNotional(double newNotional) {
    const static string method = "CredDefSwap::setNotional";
    if (notional==newNotional) {
        return;
    } else if (notional==0.0) {
        throw ModelException(method,"You cannot reset the notional of a CredDefSwap which has zero notional!");
    }
    double ntlRatio = newNotional/notional;
    notional = newNotional;
    int i=0;

    for (i=0; i<feePayments->size(); i++) {
        (*feePayments)[i].amount *= ntlRatio;
    }
    return;

}


/** Obtains the fee leg cashflows using the credit event override if
    available.
    Depending on the value of the pricingAccrue flag, two things can be
    configured:
    - Riskless fees are only included if pricingAccrue is false (ie,
      when pricing a defaulted fee leg).
    - Cashflows on "pretendedValueDate" may get included in the valuation
      for backwards compatibility under certain conditions, but only if
      pricingAccrue is false.
    The "pretendedDefaultDate" and "pretendedValueDate" are
    required to support the (THETA_)ACCRUED_INTEREST output requests */
double CredDefSwap::priceFeeLegGivenDefault(
    const DateTime& pretendedDefaultDate,
    const DateTime& pretendedValueDate,
    int& feeIdx,
    IForwardRatePricerSP model,
    const bool pricingAccrue) const
{
    const static string method("CredDefSwap::priceFeeLegGivenDefault");

    // Flag to indicate whether payments on valueDate should be included in
    // the valuation or not. In general they should not but for compatibility
    // with the old CredDefSwap implementation, if the triggerDelay and
    // defaultToSettlementDelay are NOT passed in, they will be included
    bool includeTodaysPayments = false;

    const DateTime& defaultDate = pretendedDefaultDate;
    int idx = 0;
    double feeLegValue = 0.0;

    if (defaultDate < swapEffectiveDate) { // aka protectionStartDate
        // If the default happenend before the protection start date,
        // ignore all risky fees
    }
    else {
        DateTime determinationDate;
        DateTime accrualPayDate;
        if (defaultDate > getProtectionEndDate()) {
            // All fees need to be accrued riskless since the default cannot be
            // triggered in this contract. No need to look into the override here,
            // since whatever the eventDeterminationDate etc are, those parameters
            // will not be used
            // This is flagged by having an empty determinationDate
    
            // determinationDate and accrualPayDate are already empty
        }
        else {
            if (!creditEventOverride) {
                // There is no override, so use the values in the instrument
                if (!triggerDelay) {
                    includeTodaysPayments = !pricingAccrue;
    
                    // No adjustment is done here.
                    determinationDate = defaultDate;
                    accrualPayDate = defaultDate;
                }
                else {
                    // Note this means both triggerDelay and defaultToSettlementDelay
                    // are not NULL.
                    // Obtain the determination and payment dates using this
                    // instruments delay parameters + adjustments
                    determinationDate =
                        ITrancheCreditEventOverride::rollAndAdjustDate(
                            defaultDate,
                            triggerDelay,
                            pretendedValueDate,
                            pretendedValueDate,
                            lastTriggerDate,
                            IBadDayAdjusterConstSP::attachToRef(this)); // badDayAdjuster
    
                    accrualPayDate =
                        ITrancheCreditEventOverride::rollAndAdjustDate(
                            defaultDate,
                            defaultToSettlementDelay,
                            pretendedValueDate,
                            determinationDate,
                            lastTriggerDate,
                            IBadDayAdjusterConstSP::attachToRef(this)); // badDayAdjuster
                }
            }
            else {
                FeeLegReductionPerDefaultArrayConstSP reductions =
                    creditEventOverride->historicFeeLegReductions(
                        1.0, // notional - not used for single names
                        0.0, // recovery rate override - not used for single names
                        lastTriggerDate,
                        defaultDate,
                        IBadDayAdjusterConstSP::attachToRef(this));
                
                // Sanity-check: verify that all losses have the same determination
                // date and effective date
                int numReductions = reductions->size();
                if (!reductions || (numReductions == 0)) {
                    // Leave determinationDate and accrualPayDate as empty dates
                }
                else {
                    determinationDate = (*reductions)[0]->determinationDate;
                    accrualPayDate = (*reductions)[0]->calculationDate;
                    for (int i=1; i < numReductions; ++i) {
                        if (determinationDate != (*reductions)[i]->determinationDate) {
                            throw ModelException(
                                method, "A single name is not expected to "
                                "produce losses with different determination "
                                "dates (" +
                                determinationDate.toString() + " vs " +
                                (*reductions)[i]->determinationDate.toString() +
                                ").");
                        }
                        // The accrualPayDate is the first of all calculationDates
                        if ((*reductions)[i]->calculationDate < accrualPayDate) {
                            accrualPayDate = (*reductions)[i]->calculationDate;
                        }
                    }
                }
            }
        }
    
        // Skip past fees
        int numFees = feePayments->size();
        for (/*idx*/; idx < numFees; ++idx) {
            if ((((*feePayments)[idx].date) > pretendedValueDate) ||
                (includeTodaysPayments && 
                 ((*feePayments)[idx].date == pretendedValueDate)))
            {
                // This fee is paid after "valueDate" or on valueDate and we should
                // include such fees. So do not skip this fee, get out of the loop.
                break;
            }
        }
    
        IDecretionCurveConstSP psPrepay = cdsParSpreads->getPrepayCurve();
        double factor = psPrepay->getFactor(pretendedValueDate);
        double balance;
        // Accrue all payments up to determinationDate - Note: assumes feePayments
        // are in chronological order
        for (/*idx*/; idx < numFees; ++idx) {
            const DateTime& loDate = (idx == 0)?
                                      swapEffectiveDate :
                                      (*feePayments)[idx-1].date;
    
            if (!(determinationDate.empty()) &&
                (loDate >= determinationDate))  // JLH change to >
            {
                // The default has happened and has been triggered, and have gone
                // past the last accruable fee - we are done computing the accrual
                break;
            }
    
            double pv;
            DateTime hiDate;
            if (swpAccrualDCC->years(loDate, (*feePayments)[idx].date)) {
                // cannot divide by 0
    
                bool lastFee = (idx == (numFees - 1));
                const DateTime& accrueEndDate = lastFee ? 
                    getProtectionEndDate() : (*feePayments)[idx].date;
    
                if (determinationDate.empty() ||
                    accrueEndDate <= determinationDate)
                {
                    // If the fee is to be paid before the determinationDate,
                    // accrue it all and discount it to the fee payment date.
                    // Note this is equivalent to this payment being riskless
                    // since we are stating that this name will not default
                    // between today and fee payment date.
                    // Caution: need to accrue up to accrueEndDate, but the 
                    // payment happens on (*feePayments)[idx].date
                    hiDate = accrueEndDate;
                    pv = discount->pv(pretendedValueDate, (*feePayments)[idx].date);
                    balance = psPrepay->pv((*feePayments)[idx].date); // JLHP should this be hiDate?
                }
                else {
                    if (!payAccruedFee ||
                        (accrualPayDate < pretendedValueDate) ||
                        ((accrualPayDate == pretendedValueDate) &&
                         !(includeTodaysPayments || pricingAccrue)))
                    {
                        // We are done if:
                        // - we are not paying accrued fees
                        // - the accrued fee has been paid in the past
                        // - the accrued fee is to be paid today AND we are not
                        //   + including today's payments, nor
                        //   + computing the accrual
                        break;
                    }
    
                    // Otherwise accrue up to determinationDate. The payment
                    // will happen on accrualPayDate
                    hiDate = determinationDate;
                    pv = discount->pv(pretendedValueDate, accrualPayDate);
                    balance = psPrepay->pv(accrualPayDate); // JLH should this be hiDate?
                }
    
                feeLegValue += pv * (*feePayments)[idx].amount *
                    swpAccrualDCC->years(loDate, hiDate) /
                    swpAccrualDCC->years(loDate, (*feePayments)[idx].date) *
                    factor * balance;
            }
            // else, nothing to accrue here
        }
    }
    
    // Add riskless cashflows if required, including payments on valueDate if 
    // required
    if (!pricingAccrue && !!payAsYouGo) {
        feeLegValue += discount->pv(*payAsYouGo, 
                                    pretendedValueDate,
                                    !includeTodaysPayments);
    }

    // let the calling function know the relevant fee payment index of
    // the last risky fee we dealt with
    feeIdx = (idx > 0) ? idx - 1 : 0;
    return feeLegValue;
}



/** Prices the contingent leg of a defaulted CDS, using the credit event
 * override if available */
double CredDefSwap::priceContingentLegGivenDefault() const {
    const DateTime& curveDefaultDate = cdsParSpreads->getDefaultDate();

    CashFlowArraySP contingentCashFlows(new CashFlowArray(0));
    if((curveDefaultDate < swapEffectiveDate) || // aka protectionStartDate
       (curveDefaultDate > getProtectionEndDate()))
    {
        // Default happened outside the protection period, so
        // contingentCashFlows is empty
    }
    else {
        if (!creditEventOverride) {
            // There is no override, so use the values in the instrument
            // Flag to indicate whether payments on valueDate should be included
            // in the valuation or not. In general they should not but for
            // compatibility with the old CredDefSwap implementation, if the
            // triggerDelay and  defaultToSettlementDelay are NOT passed in,
            // they will be included
            bool includeTodaysPayments = false;
            DateTime defaultPaymentDate;

            if (!defaultToSettlementDelay) {
                includeTodaysPayments = true;

                // No adjustment is done here.
                defaultPaymentDate = curveDefaultDate;
            }
            else {
                // Note this means both triggerDelay and defaultToSettlementDelay
                // are not NULL.
                defaultPaymentDate =
                    ITrancheCreditEventOverride::rollAndAdjustDate(
                        curveDefaultDate,
                        defaultToSettlementDelay,
                        valueDate,
                        valueDate,
                        lastTriggerDate,
                        IBadDayAdjusterConstSP::attachToRef(this)); // badDayAdjuster
            }

            // Generate the one cashflow "manually"
            if ((curveDefaultDate <= lastTriggerDate) &&
                ((defaultPaymentDate > valueDate) ||
                 (includeTodaysPayments && (defaultPaymentDate == valueDate))))
            {
                // Notional < 0 means short protection
                // Notional > 0 means long protection
                contingentCashFlows->push_back(
                    CashFlow(defaultPaymentDate,
                             notional * (1.0 - getRecovery())));
            }
            else {
                // Since this contingent payment would have been paid before
                // today or the default was triggered too late, return an
                // empty cash flow array
            }
        }
        else {
            // Use the override
            CtgLegLossPerDefaultArraySP contingentReductions =
                creditEventOverride->historicContingentLegLosses(
                    notional,
                    getRecovery(),
                    lastTriggerDate,
                    curveDefaultDate,
                    IBadDayAdjusterConstSP::attachToRef(this)); // BadDayAdjuster

            CashFlowArraySP allCashFlows = 
                CtgLegLossPerDefault::getReductions(contingentReductions);

            // Insert all future cashflows
            int numCashFlows = allCashFlows->size();
            contingentCashFlows->reserve(numCashFlows);
            for (int i=0; i < numCashFlows; ++i) {
                CashFlow cf((*allCashFlows)[i]); // for ease
                if (cf.date > valueDate) {
                    contingentCashFlows->push_back(cf);
                }
            }
        }
    }

    double contingentLegValue = 0.0;
    auto_ptr<IYieldCurve::IKey> key(discount->logOfDiscFactorKey());
    // pv feeCashFlows
    for (int i=0; i < contingentCashFlows->size(); ++i) {
        CashFlow cf = (*contingentCashFlows)[i]; // For ease
        const double pv = exp(key->calc(valueDate, cf.date));
        contingentLegValue += pv * cf.amount;
    }

    IDecretionCurveConstSP psPrepay = cdsParSpreads->getPrepayCurve();
    const double factor = psPrepay->getFactor(valueDate);
    const double balance = psPrepay->pv(curveDefaultDate);

    contingentLegValue *= factor * balance;

    return contingentLegValue;
}


/** Returns the CDS price when the underlying name has defaulted */
double CredDefSwap::getPriceIfDefaulted(Control* control,
                                        IForwardRatePricerSP model) const {
    static const string method = "CredDefSwap::getPriceIfDefaulted";

    const DateTime& defaultDate = cdsParSpreads->getDefaultDate();

    int dummy;
    double feeLegValue = 
        priceFeeLegGivenDefault(defaultDate, valueDate, dummy, model, false);

    double contingentLegValue = priceContingentLegGivenDefault();

    const double price = feeLegValue + contingentLegValue;
    return price;
}


bool CredDefSwap::recurse(const CFieldConstSP& field,
                          const CClassConstSP& targetClass) const
{
    // this gets called as part of the tweaking and allows us to specify
    // whether the fields within this class should be tweaked or not. The
    // target class indicates what is being shifted
    if (field == assetField ) {
        if (targetClass == ITweakableWithRespectTo<VolParallel>::TYPE   ||
            targetClass == IRestorableWithRespectTo<VolParallel>::TYPE  ||
            targetClass == ITweakableWithRespectTo<VolPointwise>::TYPE  ||
            targetClass == IRestorableWithRespectTo<VolPointwise>::TYPE ||
            targetClass == VegaMatrix::IShift::TYPE                     ||
            targetClass == VegaMatrix::IRestorableShift::TYPE           ||
            targetClass == RootTimeVega::IShift::TYPE                   ||
            targetClass == RootTimeVega::IRestorableShift::TYPE         ||
            targetClass == VegaSkewParallel::IShift::TYPE               ||
            targetClass == VegaSkewParallel::IRestorableShift::TYPE     ||
            targetClass == VegaSkewPointwise::IShift::TYPE              ||
            targetClass == VegaSkewPointwise::IRestorableShift::TYPE    ||
            targetClass == CreditSpreadRhoParallel::Shift::TYPE         ||
            targetClass == CreditSpreadRhoPointwise::IShift::TYPE       ||
            targetClass == ITweakableWithRespectTo<Spot>::TYPE          ||
            targetClass == IRestorableWithRespectTo<Spot>::TYPE         )
        {
            return false;
        }
    }

    return true;
}

void CredDefSwap::getMarket(const IModel* model, const MarketData* market)
{
    if (ClosedFormCDSPSandFA::TYPE->isInstance(model) ||
        ClosedFormFA::TYPE->isInstance(model))  {
        createE2Csensitivities = true;
    } else {
        createE2Csensitivities = false;
    }
    if (!swpAccrualDCC) {
        // default to par curve DCC
        swpAccrualDCC = DayCountConventionSP(cdsParSpreads->dayCountConv());
    }

#ifdef CDS_BACKWARD_COMPATIBILITY
    cdsParSpreads->setHolidays(HolidaySP(Holiday::weekendsOnly()));
    if (!BDC.get()) {
        BDC = BadDayConventionSP(BadDayConventionFactory::make("None"));
    }
    cdsParSpreads->setBadDayConvention(BDC);
#endif

    if (settlementHols.isEmpty()) {
        // default to no holidays (so business days = calendar days)
        settlementHols.setObject(MarketObjectSP(Holiday::noHolidays()));
    }
    else {
        settlementHols.getData(model, market);
    }
}


/** Validates that the dcc field is not empty (it is flagged as optional) and
    that the bdc is not supplied (it is to be removed) */
void CredDefSwap::validateDCCSuppliedAndNoBDC() const{
    static const string method = "CredDefSwap::validateDCCSuppliedAndNoBDC";
    if (dcc.empty()) {
        throw ModelException(method,
                             "Instrument's Day Count Convention must be "
                             "supplied");
    }
    if (!bdc.empty()){
        throw ModelException(method,
                             "Bad Day Convention for par spreads should not be "
                             "specified in the insrument");
    }
}

void CredDefSwap::validatePop2Object(){
    static const string method("CredDefSwap::validatePop2Object");
    try{
        if (!bdc.empty()){
            BDC.reset(BadDayConventionFactory::make(bdc));
        }
        if (!dcc.empty()) {
            swpAccrualDCC.reset(DayCountConventionFactory::make(dcc));
        }
        if (feePayments->empty()){
            throw ModelException(method,
                                 "There must be at least one fee payment");
        }

        if (lastTriggerDate.empty()) {
            // lastTriggerDate is optional and may not be present. If so,
            // set it to the protectionEndDate
            lastTriggerDate = getProtectionEndDate();
        }
        else if (lastTriggerDate < getProtectionEndDate()) {
            throw ModelException(method,
                                 "lastTriggerDate (" +
                                 lastTriggerDate.toString() +
                                 ") cannot be before protectionEndDate (" +
                                 getProtectionEndDate().toString() + ").");
        }

        if (!!creditEventOverride &&
            (creditEventOverride->getName() != cdsParSpreads.getName()))
        {
            throw ModelException(method,
                                 "A creditEventOverride is present, but its "
                                 "name (" + creditEventOverride->getName() +
                                 ") is not the same as the underlying curve's "
                                 "name (" + cdsParSpreads.getName() + ").");
        }

        if ((!!triggerDelay && !defaultToSettlementDelay) ||
            (!triggerDelay && !!defaultToSettlementDelay))
        {
            throw ModelException(method,
                                 "triggerDelay and defaultToSettlementDelay are "
                                 "optional but both of them have to be provided "
                                 "if one is provided - alternativeley, remove "
                                 "both of them.");
        }

        if (!!triggerDelay) { // both delays are provided
            if (defaultToSettlementDelay->intValue() < triggerDelay->intValue()) {
                throw ModelException(method,
                                     "The defaultToSettlementDelay (" +
                                     Format::toString(defaultToSettlementDelay->intValue()) +
                                     ") cannot be smaller than the triggerDelay(" +
                                     Format::toString(triggerDelay->intValue()) +
                                     ").");
            }

            if (triggerDelay->intValue() < 0) {
                throw ModelException(method,
                                     "The triggerDelay is negative (" +
                                     Format::toString(triggerDelay->intValue()) +
                                     "), and negative delays are not accepted.");
            }
            if(defaultToSettlementDelay->intValue() < 0) {
                throw ModelException(method,
                                     "The defaultToSettlementDelay is negative (" +
                                     Format::toString(defaultToSettlementDelay->intValue()) +
                                     "), and negative delays are not accepted.");
            }
        }

        adjFeePayments = CashFlowArraySP(feePayments.clone());
    }
    catch (exception& e){
        throw ModelException(e, method);
    }
}


void CredDefSwap::Validate() {
    static const string method = "CredDefSwap::Validate";
    int i;

    try {
        // validate the asset specific stuff at the generic level
        validate();

        for (i = 1; i < feePayments->getLength(); i++) {
            DateTime current  = (*feePayments)[i].date;
            DateTime previous = (*feePayments)[i-1].date;
            if (previous.isGreaterOrEqual(current)) {
                throw ModelException(
                    method, "cash flow date (" + current.toString() + ") must be " +
                    "greater than the previous cashflow date " + previous.toString());
            }
        }
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }
}

double CredDefSwap::getRecovery() const {
    double recovery;
    if (!useSwapRecovery) {
        recovery = cdsParSpreads->getRecovery();
    } else {
        recovery = swapRecovery;
    }
    return recovery;
}

/** Returns the protection end date. Should be used in preference to direct
    reference to protectionEndDate field. Note that we avoid populating
    protectionEndDate so that serialised object shows whether protectionEndDate
    was supplied or not */
DateTime CredDefSwap::getProtectionEndDate() const{
    DateTime pe = (protectionEndDate.empty()?
            // not passed on construction, so default to last fee payment date
            feePayments->back().date:
            // return the input date
            protectionEndDate);

    // ensure that this is greater than protection start
    if (swapEffectiveDate.isGreater(pe))
    {
        throw ModelException(
            "CredDefSwap::getProtectionEndDate",
            "Protection end date is not later than protection start");
    }
    return pe;
}


// -- Pricing Methods -- //

void CredDefSwap::priceParSpreads(CResults* results, Control* control, IForwardRatePricerSP model) const{
    static const string method = "CredDefSwap::priceParSpreads";

    // Moving this out here (from "psDefRates = ..." below) seems to avoid
    // an obscure solaris.opt exception-catching crash.

    DefaultRatesSP psDefRates;
    IDecretionCurveConstSP psPrepay;

    try {
        if (cdsParSpreads->defaulted()) {
            double priceDefaulted = getPriceIfDefaulted(control, model);
            results->storePrice(priceDefaulted, discount->getCcy());
        } else {
            double value = 0.0;

            DateTime maturity = (*feePayments)[feePayments->size()-1].date;
            // set the rolling effective date - for bootstrapping
            DateTime effDate = cdsParSpreads->spotDate(valueDate);

            psPrepay = cdsParSpreads->getPrepayCurve();
            psDefRates = cdsParSpreads->defaultRates();

            DateTime protectionStartDate = swapEffectiveDate.max(valueDate);

            CDSPricer cdsPricer = CDSPricer(
                feePayments,
                psDefRates,
                (-1.0) * notional,
                getRecovery(),
                payAccruedFee,
                swpAccrualDCC,
                valueDate,
                protectionStartDate,
                getProtectionEndDate(),
                swapEffectiveDate,
                discount.getSP(),
                psPrepay);

            value = cdsPricer.price();
            if (!!payAsYouGo) { // add riskless cashflows
                value += discount->pv(*payAsYouGo, valueDate);
            }

            results->storePrice(value, discount->getCcy());

            if (control && control->isPricing()) {
                CashFlowArraySP cleanSpreadCurve = psDefRates->getCleanSpreadCurve();

                string ccyName = discount->getName();
                addRequests(control, results, cleanSpreadCurve,
                            cdsParSpreads->getCurrentSpreadOrUntweakable(effDate,
                                                                         maturity),
                            model);

                // IMPLIED_CDS_SPREAD and CDS_RISKY_DURATION - can't do it in addRequests at the moment
                OutputRequest* requestImpliedSpread =
                    control->requestsOutput(OutputRequest::IMPLIED_CDS_SPREAD);

                OutputRequest* requestDuration =
                    control->requestsOutput(OutputRequest::CDS_RISKY_DURATION);

                if (requestImpliedSpread || requestDuration) {
                    try {
                        CDSPricer impliedSpreadPricer = CDSPricer(
                            1.0, // coupon rate
                            CashFlow::dates(*feePayments),
                            psDefRates,
                            -1.0, // unit notional
                            getRecovery(),
                            payAccruedFee,
                            swpAccrualDCC,
                            valueDate,
                            protectionStartDate,
                            getProtectionEndDate(),
                            swapEffectiveDate,
                            discount.getSP(),
                            psPrepay);

                        double parSpread = 0.0;
                        double riskyDuration = 0.0;

                        // do only one call to impliedParSpreadAndDuration(...)
                        impliedSpreadPricer.impliedParSpreadAndDuration(parSpread,
                                                                        riskyDuration);

                        if (requestImpliedSpread) {
                            results->storeRequestResult(requestImpliedSpread,
                                                        parSpread);
                        }
                        if (requestDuration) {
                            results->storeRequestResult(requestDuration,
                                                        riskyDuration);
                        }
                    } catch (exception) {
                        if (requestImpliedSpread) {
                            results->storeNotApplicable(requestImpliedSpread);
                        }
                        if (requestDuration) {
                            results->storeNotApplicable(requestDuration);
                        }
                    }
                }
            }
        }
    }
    catch (exception& e) {
        throw ModelException(&e, method);
    }
}


void CredDefSwap::addRequests(Control*             control,
                              Results*             results,
                              CashFlowArraySP      cleanSpreadCurve,
                              IObjectSP            currentSpread,
                              IForwardRatePricerSP model) const
{
    static const string method = "CredDefSwap::addRequests";

    try {
        // CLEAN_DEFAULT_SPREAD_CURVE
        OutputRequest* request =
            control->requestsOutput(OutputRequest::CLEAN_DEFAULT_SPREAD_CURVE);

        if (request && cleanSpreadCurve.get()) {
            string name = cdsParSpreads->getName();
            OutputNameConstSP defaultOutput(new OutputName(name));
            IObjectSP    cflows(new CashFlowList(cleanSpreadCurve.get()));

            results->storeRequestResult(request, cflows, defaultOutput);
        }

        // PAYMENT_DATES
        request = control->requestsOutput(OutputRequest::PAYMENT_DATES);
        if (request) {
            DateTimeArraySP dates(paymentDates());
            OutputRequestUtil::recordPaymentDates(control,results,dates.get());
        }

        // KNOWN_CASHFLOWS, actually input cashflows
        request = control->requestsOutput(OutputRequest::KNOWN_CASHFLOWS);
        if (request) {
            CashFlowArraySP cfl(knownCashFlows());
            OutputRequestUtil::recordKnownCashflows(control,
                                                    results,
                                                    discount->getCcy(),
                                                    cfl.get());
        }

        // DEFAULT_PROBABILITY
        request = control->requestsOutput(OutputRequest::DEFAULT_PROBABILITY);
        if (request) {
            if ( !!asset) {
                DateTime maturity = (*feePayments)[feePayments->size()-1].date;
                CAsset*    ncAsset   = const_cast<CAsset*>(asset.get());
                FirmAsset* firmAsset = dynamic_cast<FirmAsset*>(ncAsset);
                double defaultProb = 1. - firmAsset->CalcNoDefaultProb(valueDate,
                                                                       maturity);
                results->storeRequestResult(request, defaultProb);
            }
            else {
                results->storeNotApplicable(request);
            }
        }

        // IMPLIED_DEFAULT_PROBABILITY at maturity
        request = control->requestsOutput(OutputRequest::IMPLIED_DEFAULT_PROBABILITY);
        try {
            if (request)  {
                if (!!cdsParSpreads) {
                   DateTime maturity = (*feePayments)[feePayments->size()-1].date;
                   // set the rolling effective date - for bootstrapping
                   DateTime effDate = cdsParSpreads->spotDate(valueDate);
                    DefaultRatesSP psDefRates = cdsParSpreads->defaultRates();

                    double defaultProb = 1.0 - psDefRates->calcDefaultPV(effDate,
                                                                         maturity);
                    results->storeRequestResult(request, defaultProb);
                } else {
                   results->storeNotApplicable(request);
                }
            }
        } catch (exception&) {
            results->storeNotApplicable(request);
        }

        // THETA_ACCRUED_INTEREST
        request = control->requestsOutput(OutputRequest::THETA_ACCRUED_INTEREST);
        if (request) {
            double thetaAI = calcThetaAccruedInterest(model);
            results->storeRequestResult(request, thetaAI);
        }

        // ACCRUED_INTEREST
        request = control->requestsOutput(OutputRequest::ACCRUED_INTEREST);
        if (request) {
            double accInt = getAccruedInterest(valueDate, model);
            results->storeRequestResult(request, accInt);
        }

        // RECOVERY_VALUE
        request = control->requestsOutput(OutputRequest::RECOVERY_VALUE);
        if (request) {
            IDecretionCurveConstSP psPrepay = cdsParSpreads->getPrepayCurve();
            double factor = psPrepay->getFactor(valueDate);
            double balance = psPrepay->pv(valueDate);
            double defaultValue = getRecovery() * notional * factor * balance;
            results->storeRequestResult(request, defaultValue);
        }

        // CURRENT_SPREAD
        request = control->requestsOutput(OutputRequest::CURRENT_SPREAD);
        if (request) {
            results->storeRequestResult(request, currentSpread);
        }

        // IND_CDS_PAR_SPREAD
        // This one is the same as CURRENT_SPREAD
        // but goes in its own packet and is qualified by curve name
        request = control->requestsOutput(OutputRequest::IND_CDS_PAR_SPREAD);
        if (request) {
            results->storeRequestResult(
                request,
                currentSpread,
                OutputNameConstSP(new OutputName(cdsParSpreads->getName())));
        }

        // PAR_SPREAD_CURVE
        request = control->requestsOutput(OutputRequest::PAR_SPREAD_CURVE);
        if (request) {
            if (!cdsParSpreads) {
                results->storeNotApplicable(request);
            }
            else {
                results->storeRequestResult(
                    request,
                    ExpiryResult::createExpiryResultArray(
                        *cdsParSpreads->getParSpreadsExpiries().get(),
                        *cdsParSpreads->getParSpreads().get()),
                    OutputNameSP(new OutputName(cdsParSpreads->getName())));
            }
        }
        // FEE_LEG_RISKLESS_FV - hopefully self-explanatory
        request = control->requestsOutput(OutputRequest::FEE_LEG_RISKLESS_FV);
        if (request) {
            IDecretionCurveConstSP psPrepay = cdsParSpreads->getPrepayCurve();
            double factor = psPrepay->getFactor(valueDate);
            double balance = psPrepay->pv(valueDate);
            double risklessPV = discount->pv(*feePayments, valueDate) * factor * balance;
            results->storeRequestResult(request, risklessPV);
        }
     }
    catch (exception& e) {
        throw ModelException(&e, method);
    }
}


double CredDefSwap::calcThetaAccruedInterest(IForwardRatePricerSP model)const {
    static const string method = "CredDefSwap::calcThetaAccruedInterest";

    try {
        HolidaySP hols(Holiday::weekendsOnly());
        BadDayConventionSP badDay(BadDayConventionFactory::make("Following"));
        int feeIdx = 0;
        double accruedInterestNow = getAccruedInterest(valueDate, model, feeIdx);

        const DateTime& feeDate = (*feePayments)[feeIdx].date;
        double feeAmount = (*feePayments)[feeIdx].amount;

        const DateTime& tomorrow = badDay->adjust(valueDate.rollDate(1), 
                                                  hols.get());

        double accruedInterestTomorrow = getAccruedInterest(tomorrow, model, feeIdx);

        // Account for the case where today and tomorrow are in different fee periods
        double thetaAI = 0.;
        if ((valueDate < feeDate && feeDate <= tomorrow) &&
            !(swapEffectiveDate >= feeDate) &&
            !(valueDate <= swapEffectiveDate && swapEffectiveDate <= feeDate))
        {
            thetaAI = (accruedInterestTomorrow + feeAmount) - accruedInterestNow;
        }
        else {
            thetaAI = accruedInterestTomorrow - accruedInterestNow;
        }

        return thetaAI;

    }
    catch (exception& e) {
        throw ModelException(&e, method);
    }
}


void CredDefSwap::priceFirmAsset(CResults* results, Control* control, IForwardRatePricerSP model) const{
    static const string method = "CredDefSwap::priceFirmAsset";
    try {
        if (cdsParSpreads->defaulted()) {
            double priceDefaulted = getPriceIfDefaulted(control, model);
            results->storePrice(priceDefaulted, discount->getCcy());
        } else {
            double value = 0.0;

            // set the rolling effective date - do it here so Theta also causes effDate to change
            DateTime effDate = cdsParSpreads->spotDate(valueDate);

            CAsset*    ncAsset   = const_cast<CAsset*>(asset.get());
            FirmAsset* firmAsset = dynamic_cast<FirmAsset*>(ncAsset);

            if (!firmAsset) {
                throw ModelException("CredDefSwap::priceFirmAsset",
                        "Underlying credit asset must be a FirmAsset");
            }

            const DateTime& maturity = (*feePayments)[feePayments->size()-1].date;

            // interpolate vol for firm asset
            firmAsset->calculateProcessedVol(maturity);

            // calculate the clean liquidity spread curve
            DefaultRatesSP cleanLiquiditySpreads(
                new CDSHelper::CParSpreadDefaultRates(valueDate));

            try {
                if (!firmAsset->getLiquiditySpreadCurve().getName().empty())  {
                    ICDSParSpreadsSP liquiditySpreads =
                        firmAsset->getLiquiditySpreadCurve().getSP();
#ifdef CDS_BACKWARD_COMPATIBILITY
                    liquiditySpreads->setBadDayConvention(BDC);
                    liquiditySpreads->setHolidays(HolidaySP(Holiday::weekendsOnly()));
#endif
                    cleanLiquiditySpreads =
                        liquiditySpreads->defaultRates();
                }
            } catch (exception &e) {
                throw ModelException(&e,
                                     "CredDefSwap::basePrice: Failed to "
                                     "calculate clean liquidity spread curve");
            }

            CDSHelper::CFirmAssetDefaultRates faDefRates(firmAsset,
                                                         cleanLiquiditySpreads);

            CDSHelper::CredDefSwapFromDefRatesObject(valueDate,
                                                     effDate,
                                                     feePayments,
                                                     (-1.)*notional,
                                                     getRecovery(),
                                                     payAccruedFee,
                                                     swapEffectiveDate,
                                                     maturity,
                                                     discount.getSP(),
                                                     &faDefRates,
                                                     swpAccrualDCC.get(),
                                                     true, /* isE2C */
                                                     &value);

            results->storePrice(value, discount->getCcy());

            if (control && control->isPricing()) {
                string ccyName = discount->getName();
                addRequests(control, results, CashFlowArraySP(   ),
                            IObjectSP(new NotApplicable()),
                            model);

                e2cBasePriceCalculated = true;
                e2cBasePrice           = value;

                // IMPLIED_CDS_SPREAD - can't do it in addRequests at the moment
                OutputRequest* request =
                    control->requestsOutput(OutputRequest::IMPLIED_CDS_SPREAD);
                if (request) {
                    try {
                        int frequency = 4;
                        double parSpread = CDSHelper::CDSParSpread(
                            valueDate,
                            swapEffectiveDate,
                            frequency,
                            (-1.)*notional,
                            getRecovery(),
                            payAccruedFee,
                            swapEffectiveDate,
                            maturity,
                            discount.getSP(),
                            true, /* isE2C */
                            &faDefRates,
                            swpAccrualDCC.get());

                        results->storeRequestResult(request, parSpread);
                    } catch (exception) {
                        results->storeNotApplicable(request);
                    }
                }
            }
        }
    }
    catch (exception& e) {
        throw ModelException(&e, method);
    }
}


// -- The model implementations and product definitions -- //

/** private class. */
class CredDefSwapClosedFormCDSPS: public ClosedFormCDSPS::IProduct{
private:
    const CredDefSwap* cf; // a reference

public:
    CredDefSwapClosedFormCDSPS(const CredDefSwap* cf): cf(cf)
    {}

    void price(ClosedFormCDSPS* model,
               Control*         control,
               CResults*        results) const
    {
        //get the forward rate pricer from the model
        IHasForwardRatePricer* ihfrp = dynamic_cast<IHasForwardRatePricer*>(model);
        if (!ihfrp) {
            throw ModelException("CredDefSwapClosedFormCDSPS::price",
                "Model must implement IHasForwardRatePricer");
        }
        IForwardRatePricerSP ifrp = ihfrp->getForwardRatePricer();

        cf->priceParSpreads(results, control, ifrp);
    }
};

/** Implementation of ClosedFormCDSPS::IntoProduct interface */
ClosedFormCDSPS::IProduct* CredDefSwap::createProduct(
    ClosedFormCDSPS* model) const{
    return new CredDefSwapClosedFormCDSPS(this);
}


/** private class */
class CredDefSwapClosedFormFA: public ClosedFormFA::IProduct{
private:
    const CredDefSwap* cf; // a reference

public:
    CredDefSwapClosedFormFA(const CredDefSwap* cf): cf(cf)
    {}

    void price(ClosedFormFA* model,
               Control*         control,
               CResults*        results) const
    {
        //get the forward rate pricer from the model
        IHasForwardRatePricer* ihfrp = dynamic_cast<IHasForwardRatePricer*>(model);
        if (!ihfrp) {
            throw ModelException("CredDefSwapClosedFormCDSPS::price",
                "Model must implement IHasForwardRatePricer");
        }
        IForwardRatePricerSP ifrp = ihfrp->getForwardRatePricer();

        cf->priceFirmAsset(results, control, ifrp);
    }
};

/** Implementation of ClosedFormCDSPS::IntoProduct interface */
ClosedFormFA::IProduct* CredDefSwap::createProduct(
    ClosedFormFA* model) const
{
    return new CredDefSwapClosedFormFA(this);
}


/** private class. */
class CredDefSwapClosedFormCDSPSandFA: public ClosedFormCDSPSandFA::IProduct{
private:
    const CredDefSwap* cf; // a reference

public:
    CredDefSwapClosedFormCDSPSandFA(const CredDefSwap* cf): cf(cf)
    {}

    void price(ClosedFormCDSPSandFA* model,
               Control*              control,
               CResults*             results) const
    {
        // retrieve the model to be used in the calculation
        // of fee leg forward rates
        IHasForwardRatePricer* ihfrp = dynamic_cast<IHasForwardRatePricer*>(model);
        if (!ihfrp)
        {
            throw ModelException(
                "CredDefSwapClosedFormCDSPSandFA::price",
                "Model must implement IHasForwardRateModel");
        }
        IForwardRatePricerSP frModel = ihfrp->getForwardRatePricer();

        if (model->isPricePar()) {
            cf->priceParSpreads(results, control, frModel);
        } else {
            cf->priceFirmAsset(results, control, frModel);
        }
    }
};

/** Implementation of ClosedFormCDSPS::IntoProduct interface */
ClosedFormCDSPSandFA::IProduct* CredDefSwap::createProduct(
    ClosedFormCDSPSandFA* model) const
{
    return new CredDefSwapClosedFormCDSPSandFA(this);
}


/** Implementation of CreditSupport::Interface interface */
CreditSupportSP CredDefSwap::createCreditSupport(CMarketDataSP market){
    return CreditSupportSP(new CDSCreditSupport(this, market));
}

// -- Stuff every product needs to do -- //

/** what's today ? */
DateTime CredDefSwap::getValueDate() const {
    return valueDate;
}

/** when to stop tweaking */
DateTime CredDefSwap::endDate(const Sensitivity* sensControl) const {
    const DateTime& sensLastDate =
        cdsParSpreads->stopTweaking(feePayments->back().date);
    if( ! sensControl )
        return sensLastDate;
    // this is a mess. Need to route endDate through model not through
    // instrument
    IModel* model = sensControl->getModel();
    ClosedFormCDSPS* cdsModel = dynamic_cast<ClosedFormCDSPS*>(model);
    return cdsModel? cdsModel->endDate(sensControl, this,
                                       sensLastDate): sensLastDate;
}

/** when do payments occur ? */
DateTimeArraySP CredDefSwap::paymentDates() const {
    DateTimeArraySP paydates(new DateTimeArray(feePayments->size()));

    for (int i = 0; i < feePayments->size(); i++) {
        (*paydates)[i] = (*feePayments)[i].date;
    }
    return paydates;
}

/** Returns all known cash flows */
CashFlowArraySP CredDefSwap::knownCashFlows() const
{
    static const string method = "CredDefSwap::knownCashFlows";

    try {
        IDecretionCurveConstSP psPrepay = cdsParSpreads->getPrepayCurve();
        double factor = psPrepay->getFactor(valueDate);

        for (int i = 0; i < feePayments->size(); i++) {
            DateTime date = (*feePayments)[i].date;
            double balance = psPrepay->pv(date.rollDate(-1)); 
            double amount = (*feePayments)[i].amount;
            //calculate payment date
            (*adjFeePayments)[i].date = psPrepay->getSettlement()->settles(date);
            (*adjFeePayments)[i].amount =  amount * factor * balance;
        }

        return adjFeePayments;

        //return feePayments;
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }
}

/** Returns name identifying vol for vega parallel */
string CredDefSwap::sensName(DeltaToCredit* shift) const
{
    static const string method = "CredDefSwap::sensName(DeltaToCredit)";

    if (createE2Csensitivities && !!asset)
    {
        if (!FirmAsset::TYPE->isInstance(asset.get())) {
            throw  ModelException(method, "asset is not of type FirmAsset");
        }

        OutputNameArrayConstSP names(
            RiskProperty<Spot>().subjectNames(asset.getSP()));

        if (names->size() > 1) {
            throw ModelException("CorporateBond::sensName(DeltaToCredit)",
                    "E2C currently not supported for baskets!");
        }

        return names->empty() ? "" : (*names)[0]->toString();
    }
    else {
        return "";
    }
}

/** Shifts the object using given shift */
bool CredDefSwap::sensShift(DeltaToCredit* shift)
{
    static const string method = "CredDefSwap::sensShift(DeltaToCredit)";

    try {
        if (createE2Csensitivities) {
            if (!FirmAsset::TYPE->isInstance(asset.get())) {
                throw  ModelException(method, "asset is not of type FirmAsset");
            }

            TweakOutcome outcome = PropertyTweakHypothesis<Spot>(
                shift->getShiftSize(), shift->getMarketDataName()).
                    applyTo_TweakOutcome(asset.getSP());

            shift->setInitialValue(outcome.oldValue());
        }
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }

    return false; // none of our components has a DeltaToCredit type sensitivity
}

double CredDefSwap::calcParSpreadEquivalent(const double& lsSpreadRho) const
{
    try {
        DateTime maturity = (*feePayments)[feePayments->size()-1].date;
        CAsset*    ncAsset   = const_cast<CAsset*>(asset.get());
        FirmAsset* firmAsset = dynamic_cast<FirmAsset*>(ncAsset);

        HolidaySP hols(Holiday::weekendsOnly());
        double currentParSpread       = cdsParSpreads->getCurrentSpread(valueDate,
                                                                        maturity,
                                                                        BDC.get(),
                                                                        hols.get());
        double currentLiquiditySpread = firmAsset->getLiquiditySpread(valueDate,
                                                                      maturity,
                                                                      BDC.get(),
                                                                      hols.get());
        double adjustmentFactor;
        if ( currentLiquiditySpread > 0.0 ) {
            adjustmentFactor = currentParSpread/currentLiquiditySpread; // this is equivalent to 1/(1-rho)
        } else {
            adjustmentFactor = 0.0;
        }
        return lsSpreadRho * adjustmentFactor;
    } catch (exception& e){
        throw ModelException(e, "CredDefSwap::calcParSpreadEquivalent");
    }
}

/** converts a CDS par spread curve into a clean spread curve. This method lives here, because it uses the CDS pricing
    heavily. The CDS code implements its own default rate classes, which I chose not to move into the market directory,
    as the CDS pricing code and the default rate classes interfere heavily with each other - instead, a new class
    CleanSpreadCurve is born (ASe) */
CleanSpreadCurveSP CredDefSwap::getCleanSpreadCurve(CDSParSpreadsSP parSpreadCurve,
                                                    YieldCurveSP    discountCurve,
                                                    const DateTime& valueDate,
                                                    const bool      returnFwdRates)
{
    DefaultRatesSP psDefRates =
        parSpreadCurve->defaultRates();

    CDoubleArraySP rates;

    // convert forward rates to spot rates if requested
    if (!returnFwdRates) {
        rates = CashFlow::amounts(
            *psDefRates->convertToSpot(
                psDefRates->getValueDate(),
                Actual365F()));
    } else {
        rates = CDoubleArraySP(
            new CDoubleArray(psDefRates->getRates()));
    }

    // convert the default dates into an expiry array
    DateTimeArraySP defaultDates = psDefRates->getDefaultDates();
    ExpiryArraySP expiries(new ExpiryArray(defaultDates->size()));
    for (int i = 0; i < defaultDates->size(); ++i) {
        (*expiries)[i] = ExpirySP(new BenchmarkDate((*defaultDates)[i]));
    }

    CleanSpreadCurveSP cleanSpreadCurve(
        new CleanSpreadCurve(
            valueDate,
            expiries.get(),
            rates.get()));

    return cleanSpreadCurve;
}


//++++++++++++++++++++++++++++++++++++++++
//  IBadDayAdjuster methods
//
/** Returns "date" bad day adjusted using the bad day convention
 * and holidays in this object */
DateTime CredDefSwap::badDayAdjust(const DateTime& date) const {
    return BDC->adjust(date, settlementHols.get());
}

/** Add a number of business days to a date */
DateTime CredDefSwap::addBusinessDays(const DateTime& from, int busDays) const {
    return settlementHols->addBusinessDays(from, busDays);
}
//
//  IBadDayAdjuster methods
//------------------------------------------


// for reflection
CredDefSwap::CredDefSwap():
    Generic1FactorCredit(TYPE), feeRate(1.0), useSwapRecovery(false),
    swapRecovery(0), payAccruedFee(true), payAsYouGo(0),
    createE2Csensitivities(false), e2cBasePriceCalculated(false),
    e2cBasePrice(0.0)
{}

CredDefSwap::CredDefSwap(CClassConstSP clazz):
    Generic1FactorCredit(clazz), feeRate(1.0), useSwapRecovery(false),
    swapRecovery(0), payAccruedFee(true), payAsYouGo(0),
    createE2Csensitivities(false), e2cBasePriceCalculated(false),
    e2cBasePrice(0.0)
{}


CredDefSwap::CredDefSwap(const DateTime&        valueDate,
                         const bool&            oneContract,
                         const double&          notional,
                         const string&          ccyTreatment,
                         InstrumentSettlementSP instSettle,
                         InstrumentSettlementSP premiumSettle,
                         CAssetWrapper          asset,
                         ICDSParSpreadsWrapper  cdsParSpreads,
                         YieldCurveWrapper      discount,
                         const DateTime&        swapEffectiveDate,
                         const DateTime&        protectionEndDate,
                         CashFlowArraySP        feePayments,
                         const bool&            useSwapRecovery,
                         const double&          swapRecovery,
                         const bool&            payAccruedFee,
                         const string&          dcc,
                         const string&          bdc,
                         HolidayWrapper         settlementHols):
    Generic1FactorCredit(TYPE),
    swapEffectiveDate(swapEffectiveDate),
    protectionEndDate(protectionEndDate),
    feeRate(1.0),
    feePayments(feePayments),
    useSwapRecovery(useSwapRecovery),
    swapRecovery(swapRecovery),
    payAccruedFee(payAccruedFee),
    dcc(dcc),
    bdc(bdc),
    payAsYouGo(0),
    settlementHols(settlementHols),
    createE2Csensitivities(false),
    e2cBasePriceCalculated(false),
    e2cBasePrice(0.0)
{
    this->valueDate     = valueDate;
    this->oneContract   = oneContract;
    this->notional      = notional;
    this->ccyTreatment  = ccyTreatment;
    this->instSettle    = instSettle;
    this->premiumSettle = premiumSettle;
    this->asset         = asset;
    this->cdsParSpreads = cdsParSpreads;
    this->discount      = discount;
}

class CredDefSwapHelper {
public:
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(CredDefSwap, clazz);
        SUPERCLASS(Generic1FactorCredit);
        IMPLEMENTS(ClosedFormCDSPS::IIntoProduct);
        IMPLEMENTS(ClosedFormFA::IIntoProduct);
        IMPLEMENTS(ClosedFormCDSPSandFA::IIntoProduct);
        IMPLEMENTS(LastSensDate);
        IMPLEMENTS(DeltaToCredit::IShift);
        IMPLEMENTS(ObjectIteration::IOverride);
        IMPLEMENTS(IGetMarket);
        IMPLEMENTS(ICDS);
        IMPLEMENTS(ICDSConvention);
        IMPLEMENTS(IBadDayAdjuster);
        EMPTY_SHELL_METHOD(defaultCredDefSwap);

        FIELD(swapEffectiveDate,      "protection and accrual start date");
        FIELD(protectionEndDate,      "protection end date");
        FIELD(feePayments,            "cash flows");
        FIELD(useSwapRecovery,        "use passed swap recovery value");
        FIELD(swapRecovery,           "overidden by par recovery rate if "
                                      "useSwapRecovery = false");
        FIELD(payAccruedFee,          "make accrued payment on default?");
        FIELD(dcc,                    "day count conv for accrual periods");
        FIELD(bdc,                    "bad day convention"); // ideally to be removed
        FIELD(swpAccrualDCC,          "swap accrual day count convention");
        FIELD(BDC,                    "bad day convention");
        FIELD(createE2Csensitivities, "whether to create debt/equity sensitivities");
        FIELD(e2cBasePriceCalculated, "internal field");
        FIELD(e2cBasePrice,           "internal field");
        FIELD(adjFeePayments,         "adjusted fee payments by prepay and current factors");
        FIELD(feeRate,                "Fee rate of the CDS. If missing, defaults to 1.0");
        FIELD(creditEventOverride,    "Override for credit event parameters");
        FIELD(triggerDelay,           "Delay, in days, between default and "
                                      "eventDeterminationDate - used to "
                                      "estimate eventDeterminationDate "
                                      "before it is known.");
        FIELD(defaultToSettlementDelay,"Delay between credit event and "
                                      "settlement.");
        FIELD(lastTriggerDate,        "Last date when a default occurred during "
                                      "the protection period can be triggered.");
        FIELD(payAsYouGo,             "Ad-hoc pay-as-you-go riskless payments");
        FIELD(settlementHols,         "Holidays for the tranche. Used to "
                                      "adjust the delays regarding defaults' "
                                      "settlements. Default: weekends only");

        FIELD_MAKE_OPTIONAL(payAsYouGo);
        FIELD_MAKE_OPTIONAL(protectionEndDate);
        FIELD_MAKE_OPTIONAL(swapRecovery);
        FIELD_MAKE_OPTIONAL(dcc);
        FIELD_MAKE_OPTIONAL(bdc); // really belongs on CDSParSpreads
        FIELD_MAKE_OPTIONAL(feeRate);
        FIELD_MAKE_OPTIONAL(creditEventOverride);
        FIELD_MAKE_OPTIONAL(triggerDelay);
        FIELD_MAKE_OPTIONAL(defaultToSettlementDelay);
        FIELD_MAKE_OPTIONAL(lastTriggerDate);
        FIELD_MAKE_OPTIONAL(settlementHols);
        FIELD_MAKE_TRANSIENT(swpAccrualDCC);
        FIELD_MAKE_TRANSIENT(BDC);
        FIELD_MAKE_TRANSIENT(createE2Csensitivities);
        FIELD_MAKE_TRANSIENT(e2cBasePriceCalculated);
        FIELD_MAKE_TRANSIENT(e2cBasePrice);
        FIELD_MAKE_TRANSIENT(adjFeePayments);

        // look up field for use on recurse
        cdsParSpreadsField = clazz->getSuperClass()->getDeclaredField("cdsParSpreads");
        assetField         = clazz->getSuperClass()->getDeclaredField("asset");
    }

    static IObject* defaultCredDefSwap(){
        return new CredDefSwap();
    }
};

CClassConstSP const CredDefSwap::TYPE = CClass::registerClassLoadMethod(
    "CredDefSwap", typeid(CredDefSwap), CredDefSwapHelper::load);



// REDUCED INTERFACE CLASS
// Class wraps around CredDefSwap and generates cashflows from parameters

// for reflection
CredDefSwap::QuickPricer::QuickPricer(): Generic1FactorCredit(TYPE),
                                         swapRecovery(0.0),
                                         dcc("Act/360")
{}

// CredDefSwap::QuickPricer::QuickPricer(const DateTime&         valueDate,
//                                       const bool&             oneContract,
//                                       const double&           notional,
//                                       const string&           ccyTreatment,
//                                       InstrumentSettlementSP  instSettle,
//                                       InstrumentSettlementSP  premiumSettle,
//                                       CAssetWrapper           asset,
//                                       ICDSParSpreadsWrapper    cdsParSpreads,
//                                       YieldCurveWrapper       discount,
//                                       const DateTime&         swapEffectiveDate,
//                                       const bool&             useSwapRecovery,
//                                       const double&           swapRecovery,
//                                       const bool&             payAccruedFee,
//                                       const string&           dcc,
//                                       const string&           bdc,
//                                       const double&           dealSpread,
//                                       const int&              frequency,
//                                       const DateTime&         maturityDate) :
//     Generic1FactorCredit(TYPE),
//     swapEffectiveDate(swapEffectiveDate),
//     useSwapRecovery(useSwapRecovery),
//     swapRecovery(swapRecovery),
//     payAccruedFee(payAccruedFee),
//     dcc(dcc),
//     bdc(bdc)
// {
//     this->valueDate     = valueDate;
//     this->oneContract   = oneContract;
//     this->notional      = notional;
//     this->ccyTreatment  = ccyTreatment;
//     this->instSettle    = instSettle;
//     this->premiumSettle = premiumSettle;
//     this->asset         = asset;
//     this->cdsParSpreads = cdsParSpreads;
//     this->discount      = discount;
// }

void CredDefSwap::QuickPricer::validatePop2Object() {
    static const string method = "CredDefSwap::QuickPricer::validatePop2Object";
    try {
        // Check frequency equates to an integer number of months
        int noMonths = 0;
        if (12 % frequency != 0) {
            throw ModelException("method",
                                 "Payment frequency (" +
                                 Format::toString(frequency) +
                                 ") should equate to an integer number of months.");
        } 
        else {
            noMonths = 12 / frequency;
        }

        // Need the DCC and BDC for cashflow generation
        DayCountConventionSP DCC(DayCountConventionFactory::make(dcc));
        BadDayConventionSP BDC(BadDayConventionFactory::make(bdc));

        // Generate cashflows based on unit notional. Don't have access to holidays so use
        // all days but weekends for now.
        HolidaySP hols(Holiday::weekendsOnly());
        StubSP stub(StubFactory::make("Simple"));
        CashFlowArray cashFlows = SwapTool::cashflows(
            swapEffectiveDate,  // cashflows start from Effective date
            maturityDate,
            stub.get(),         // stub
            true,               // true: stub at end, false: at beginning
            BDC.get(),          // Accrual BDC
            BDC.get(),          // Redemption BDC (n/a)
            hols.get(),         // holidays
            false,              // Keep start date?
            false,              // Add redemption?
            dealSpread,         // rate
            noMonths,           // interval = count periods
            "M",                // period
            DCC.get());

        // Multiply by notional and change sign
        for (int i=0; i<cashFlows.size(); i++) {
            cashFlows[i].amount = -cashFlows[i].amount * notional;
        }

        CashFlowArraySP feePayments(new CashFlowArray(cashFlows));

        // Create fat CDS
        cds = CredDefSwapSP(new CredDefSwap(
            valueDate,
            oneContract,
            notional,
            ccyTreatment,
            instSettle,
            premiumSettle,
            asset,
            cdsParSpreads,
            discount,
            swapEffectiveDate,
            maturityDate,  //TODO : may not be correct
            feePayments,
            useSwapRecovery,
            swapRecovery,
            payAccruedFee,
            dcc,
            bdc,
            settlementHols));

        cds->validatePop2Object();
    }
    catch (exception &e) {
        throw ModelException(e, method);
    }
}

// Wrapper functions around cds
// Generic1FactorCredit interface
void CredDefSwap::QuickPricer::Validate()
{
    cds->Validate();
}

void CredDefSwap::QuickPricer::GetMarket(const IModel* model,
                                         const CMarketDataSP market)
{
    (const_cast<IModel*>(model))->getInstrumentAndModelMarket(market.get(),
                                                              cds.get());
}

DateTime CredDefSwap::QuickPricer::getValueDate() const
{
    return cds->getValueDate();
}

// IIntoProduct interfaces
ClosedFormCDSPS::IProduct* CredDefSwap::QuickPricer::createProduct(
    ClosedFormCDSPS* model) const
{
    return cds->createProduct(model);
}

ClosedFormFA::IProduct* CredDefSwap::QuickPricer::createProduct(
    ClosedFormFA* model) const
{
    return cds->createProduct(model);
}

ClosedFormCDSPSandFA::IProduct* CredDefSwap::QuickPricer::createProduct(
    ClosedFormCDSPSandFA* model) const
{
    return cds->createProduct(model);
}

// CreditSupport interface
CreditSupportSP CredDefSwap::QuickPricer::createCreditSupport(CMarketDataSP market)
{
    return CreditSupportSP(new CDSCreditSupport(cds.get(), market));
}

// LastSensDate interface
DateTime CredDefSwap::QuickPricer::endDate(const Sensitivity* sensControl) const
{
    return cds->endDate(sensControl);
}

// IGetMarket interface
void CredDefSwap::QuickPricer::getMarket(const IModel* model,
                                         const MarketData* market)
{
    cds->getMarket(model, market);
}

CashFlowArraySP CredDefSwap::QuickPricer::knownCashFlows() const
{
    return cds->knownCashFlows();
}

DateTimeArraySP CredDefSwap::QuickPricer::paymentDates() const
{
    return cds->paymentDates();
}

// note needs to have unique name in our namespace
class CDSQuickPricerHelper{
public:
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(CredDefSwap::QuickPricer, clazz);
        SUPERCLASS(Generic1FactorCredit);
        IMPLEMENTS(ClosedFormCDSPS::IIntoProduct);
        IMPLEMENTS(ClosedFormFA::IIntoProduct);
        IMPLEMENTS(ClosedFormCDSPSandFA::IIntoProduct);
        IMPLEMENTS(LastSensDate);
        IMPLEMENTS(IGetMarket);
        EMPTY_SHELL_METHOD(defaultInterface);
        FIELD(swapEffectiveDate,     "effective date");
        FIELD(useSwapRecovery,       "use passed swap recovery value")
        FIELD(swapRecovery,          "overidden by par recovery rate if useSwapRecovery = false");
        FIELD(payAccruedFee,         "make accrued payment on default?");
        FIELD(dcc,                   "day count conv for accrual periods");
        FIELD(bdc,                   "bad day convention");
        FIELD(dealSpread,            "fee payment spread");
        FIELD(frequency,             "payment frequency");
        FIELD(maturityDate,          "swap maturity date");
        FIELD(cds,                          "credit default swap");
        FIELD(settlementHols,        "Holidays for the tranche. Used to "
                                            "adjust the delays regarding defaults' "
                                            "settlements. Default: weekends only");
        FIELD_MAKE_TRANSIENT_BUT_TWEAKABLE(cds);
        FIELD_MAKE_OPTIONAL(swapRecovery);
        FIELD_MAKE_OPTIONAL(dcc);
        FIELD_MAKE_OPTIONAL(settlementHols);

        Addin::registerConstructor("QP_CDS",
                                   Addin::CONV_BOND,
                                   "Creates a CDS instrument",
                                   CredDefSwap::QuickPricer::TYPE);
    }

    static IObject* defaultInterface(){
        return new CredDefSwap::QuickPricer();
    }
};

CClassConstSP const CredDefSwap::QuickPricer::TYPE =
CClass::registerClassLoadMethod(
    "CredDefSwap::QuickPricer",
    typeid(CredDefSwap::QuickPricer),
    CDSQuickPricerHelper::load);


bool CredDefSwapLoad() {
    return (CredDefSwap::TYPE != 0);
}

DRLIB_END_NAMESPACE
