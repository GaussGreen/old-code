//----------------------------------------------------------------------------
//
//   Group       : Derivatives Research
//
//   Filename    : VanillaCreditFeeLeg.cpp
//
//   Description : Implementation of a fee leg a Credit Vanilla Instrument :
//
//   Author      : Doris Morris
//
//   Date        : February 2006
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/VanillaCreditFeeLeg.hpp"
#include "edginc/SensControl.hpp"
#include "edginc/CDSHelper.hpp" /* To generate fee payments */
#include "edginc/MaturityPeriod.hpp"
#include "edginc/CashFlow.hpp"
#include "edginc/DateTime.hpp"
#include "edginc/StubFactory.hpp"
#include "edginc/SwapTool.hpp"
#include "edginc/AdhocCashFlow.hpp"
#include "edginc/IForwardRatePricer.hpp"
#include "edginc/IDecretionCurve.hpp"
#include "edginc/IBadDayAdjuster.hpp"
#include "edginc/ICreditEventOverrideName.hpp"

DRLIB_BEGIN_NAMESPACE

//----------------------------
// VanillaCreditFeeLeg methods
//----------------------------

VanillaCreditFeeLeg::VanillaCreditFeeLeg(
    double baseNotional,                   // initial notional
    const  DateTime& valueDate,            // when to value the fee leg as of
    const  DateTime& effectiveDate,        // protection start date
    const  DateTime& maturityDate,         // protection end date
    double feeRate,                        // fixed fee spread
    const  MaturityPeriodSP paymentFreq,   // payment frequency
    const  DateTime& fcd,                  // first full coupon
    const  DateTime& lcd,                  // last full coupon
    const  DayCountConventionSP pmtDcc,    // payment day count convention
    const  DayCountConventionSP accrualDcc,// accrual day count convention
    const  BadDayConventionSP pmtBdc,      // payment bad day convention
    const  BadDayConventionSP accrualBdc,  // accrual bad day convention
    const  BadDayConventionSP valuationBdc,// valuation bad day convention
    const  StubSP stubPaymentType,         // stub pmt type
    const  HolidayWrapper pmtHol,          // holiday calendar for payments
    const  HolidayWrapper accrualHol,      // accrual holidays
    const  YieldCurveWrapper discountCrv,  // discount curve
    bool   payAccruedFee,
    double delay) 
    : 
    CObject(TYPE),
    baseNotional(baseNotional),
    valueDate(valueDate),
    effectiveDate(effectiveDate),
    maturityDate(maturityDate),
    feeRate(feeRate),
    paymentFreq(paymentFreq),
    fcd(fcd),
    lcd(lcd),
    pmtDcc(pmtDcc),
    accrualDcc(accrualDcc),
    pmtBdc(pmtBdc),
    accrualBdc(accrualBdc),
    valuationBdc(valuationBdc),
    pmtHol(pmtHol),
    accrualHol(accrualHol),
    discount(discountCrv),
    payAccruedFee(payAccruedFee),
    delay(delay),
    stubPaymentType(stubPaymentType)
{
    /* generate the fee cash flow payments */
    feePayments = generateFeePayments();

    /* TODO: Check to see if these optional parameters are populated */
    /* fcd, lcd */

    // ck to see values for baseNotional, delay and payAccruedFee
}

VanillaCreditFeeLeg::~VanillaCreditFeeLeg(){}

/**
 * Accrued interest (cash amount, not percentage) for settlement on settlementDate.
 * Assumes fee payments are in chronological order.
 */
double VanillaCreditFeeLeg::getAccruedInterest(const DateTime& settlementDate) const
{
    const static string method("VanillaCreditFeeLeg::getAccruedInterest");

    DateTime      adjSettlementDate;
    DateTimeArray feeDates = CashFlow::dates(*feePayments);
    int           numDates;
    int           idx;
    double        ai = 0.0;
    double        pv = 0.0;


    // Adjust settlement date if delay
    if(delay > 0)
    {
        // Roll the settlementDate for the delay
        adjSettlementDate = settlementDate.rollDate((int)delay);

        // Adjust for bdc
        adjSettlementDate = accrualBdc->adjust(adjSettlementDate, accrualHol.get());
    }


    //
    // Settlement date is before protection start date or after protection end date.
    // Hence ignore all fees.
    //
    if(settlementDate < effectiveDate || settlementDate > maturityDate )
    {
        return 0.0;
    }


    // Find the coupon dates that straddle the settlement date for the accrued interest
    numDates = feeDates.size();
    if(feeDates.size() < 1)
    {
        throw ModelException(method, "No fee dates.");
    }

    for (idx = 0; idx < (*feePayments).size(); idx ++)
    {
        const DateTime& loDate = (idx == 0)? effectiveDate : (*feePayments)[idx-1].date;

        // Determine if the hiDate ends on the next fee payment date or on the protection end date
        const DateTime& hiDate = (*feePayments)[idx].date;

        // Check the start and end protection dates
        // if the settlement date falls on the a fee payment date, no accrual
        if(idx > 0)
        {
            if(settlementDate == loDate || settlementDate == hiDate)
            {
                return 0.0;
            }
        }


        // does the settlement date lie btw 2 fee payments
        if(loDate < settlementDate && settlementDate < hiDate)
        {
            // get df from value date to settlement date + delay
            if(delay > 0)
            {
                pv = discount->pv(valueDate, adjSettlementDate);
            }
            else
            {
                pv = discount->pv(valueDate, settlementDate);
            }

            // TOADD - type of dccType!!!!!! Account for SHORT_FIRST, etc

            // compute the accrual amount between the previous cpn
            // from the previous coupon date to the effective (i.e. settlement) date
            double accruedFraction = accrualDcc->accruedFactor(settlementDate,
                                                          loDate,
                                                          hiDate,
                                                          paymentFreq->annualFrequency(),
                                                          DayCountConvention::REGULAR,
                                                          false,  // eomAdjSec
                                                          false); // goneEx

            ai = pv * feeRate * accruedFraction * baseNotional;
        }
    }
    // return ai as of valuation date
    return ai;
}

/* Set the base notional*/
void VanillaCreditFeeLeg::setNotional(double newNotional)
{
    const static string method = "VanillaCreditFeeLeg::setNotional";
    double ntlRatio;
    int    i;

    if(baseNotional == 0.0)
    {
        throw ModelException(method, "You cannot reset the notional of a VanillaCreditFeeLeg which has zero notional!");
    }

    ntlRatio = newNotional/baseNotional;

    baseNotional = newNotional;

    // re-adjust fee payments
    for (i=0; i<feePayments->size(); i++)
    {
        (*feePayments)[i].amount *= ntlRatio;
    }
}

//----------------
// CObject methods
//----------------

/** Called immediately after object constructed */
void VanillaCreditFeeLeg::validatePop2Object()
{
    static const string method("VanillaCreditFeeLeg::validatePop2Object");

    /* check to see if feePayments have been generated */
    if (feePayments->empty())
    {
        throw ModelException(method, "There must be at least one fee payment generated!");
    }

    if(effectiveDate > maturityDate)
    {
        throw ModelException(method, "protection start date (" +
                             effectiveDate.toString() +
                             ") should be before protectionEndDate (" +
                             maturityDate.toString() + ").");
    }


   /* set the fcd and lcd */
}

/** Called once before the initial pricing */
void VanillaCreditFeeLeg::Validate()
{
    static const string method = "VanillaCreditFeeLeg::Validate";

    /* checks fee cashflow dates are in increasing order */
    CashFlow::ensureDatesIncreasing(*feePayments, method, true);
}

//----------------------
// ICreditFeeLeg methods
//----------------------

double VanillaCreditFeeLeg::price(
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
        IForwardRatePricerSP        model) const
{
    static const string method = "VanillaCreditFeeLeg::price";

    //for now
    throw ModelException(method, "Not yet implemented");
}

/** Return all cash flow dates */
DateTimeArraySP VanillaCreditFeeLeg::getCashFlowDates() const
{
    static const string method = "VanillaCreditFeeLeg::getCashFlowDates";

    int numFlows = feePayments->size();
    DateTimeArraySP cfDates = DateTimeArraySP(
        new DateTimeArray(numFlows));

    for (int i=0; i<numFlows; i++)
    {
        (*cfDates)[i] = (*feePayments)[i].date;
    }
    return cfDates;
}

/** Return risky cash flow dates */
DateTimeArraySP VanillaCreditFeeLeg::getRiskyCashFlowDates() const
{
    //all fees are risky in this structure
    return getCashFlowDates();
}

/** Return riskfree cash flow dates */
DateTimeArraySP VanillaCreditFeeLeg::getRisklessCashFlowDates() const
{
    //no fees are riskless in this structure
    return DateTimeArraySP(new DateTimeArray(0));
}

/** Returns the accrual periods (start and end dates) in this fee leg. */
AccrualPeriodArrayConstSP VanillaCreditFeeLeg::getAccrualPeriods() const
{
    static const string method = "VanillaCreditFeeLeg::getAccrualPeriods";

    //fee leg provides for contiguous fee periods, with no payment delay
    // (delay field is currently ignored)
    //so can back out from the payment dates
    DateTimeArraySP payDates = getCashFlowDates();
    int numFees = payDates->size();

    AccrualPeriodArraySP accruals = AccrualPeriodArraySP(
        new AccrualPeriodArray(numFees));

    for (int i=0; i<numFees; i++)
    {
        AccrualPeriodSP acc = AccrualPeriodSP(
            new AccrualPeriod(
                i==0 ? effectiveDate : (*payDates)[i-1], // accrue start
                (*payDates)[i]));                        // accrue end

        (*accruals)[i] = acc;
    }

    return accruals;

}

/** Returns risky accrual periods (start and end dates) in this fee leg. */
AccrualPeriodArrayConstSP VanillaCreditFeeLeg::getRiskyAccrualPeriods() const
{
    //all fees are risky in this structure
    return getAccrualPeriods();
}

/** Returns risky accrual periods (start and end dates) in this fee leg. */
AccrualPeriodArrayConstSP VanillaCreditFeeLeg::getRisklessAccrualPeriods() const
{
    //no fees are riskless in this structure
    return AccrualPeriodArraySP(new AccrualPeriodArray(0));
}


/** Returns all cash flows */
AbstractCashFlowArrayConstSP VanillaCreditFeeLeg::getCashFlows(IForwardRatePricerSP model) const
{
    static const string method = "VanillaCreditFeeLeg::getCashFlows";

    //convert feePayments into array of ad-hoc cashflows
    int numCfls = feePayments->size();
    AbstractCashFlowArraySP absCfls = AbstractCashFlowArraySP(
        new AbstractCashFlowArray(numCfls));

    for (int i=0; i<numCfls; i++)
    {
        //for convenience
        CashFlow cfl = (*feePayments)[i];

        //make new adhoc cashflow
        AdhocCashFlowSP cf = AdhocCashFlowSP(
            new AdhocCashFlow(cfl.date, cfl.amount));

        //and assign
        (*absCfls)[i] = cf;
    }

    return absCfls;
}

/** Returns risky cash flows only */
CashFlowArraySP VanillaCreditFeeLeg::getRiskyCashFlows(IForwardRatePricerSP model) const
{
    static const string method = "VanillaCreditFeeLeg::getRiskyCashFlows";

    //this leg type only has risky cashflows
    AbstractCashFlowArrayConstSP cfls = getCashFlows(model);
    //convert to simple cashflow type
    return AbstractCashFlow::asCashFlows(cfls, model);
}

/** Returns risk free cash flows only */
CashFlowArraySP VanillaCreditFeeLeg::getRisklessCashFlows(IForwardRatePricerSP model) const
{
    static const string method = "VanillaCreditFeeLeg::getRisklessCashFlows";

    //this leg type only has risky cashflows
    return CashFlowArraySP(new CashFlowArray(0));
}

/** Returns risky notional dates */
DateTimeArraySP VanillaCreditFeeLeg::getRiskyNotionalDates(IForwardRatePricerSP model) const
{
    static const string method = "VanillaCreditFeeLeg::getRiskyNotionalDates";

    //for now
    throw ModelException(method, "Not yet implemented");
}

/** Returns risky coupon notional types */
CouponNotionalTypesArraySP VanillaCreditFeeLeg::getRiskyCouponNotionalTypes() const
{
    static const string method = "VanillaCreditFeeLeg::getRiskyCouponNotionalTypes";

    //for now
    throw ModelException(method, "Not yet implemented");
}

/** Returns risky observation dates */
DateTimeArraySP VanillaCreditFeeLeg::getRiskyObservationDates() const
{
    static const string method = "VanillaCreditFeeLeg::getRiskyObservationDates";

    //for now
    throw ModelException(method, "Not yet implemented");
}



/** Return known cash flows corresponding to a CDO tranche */
CashFlowArraySP VanillaCreditFeeLeg::generateKnownCashFlows(
        const DateTime       today,
        const double         initialTrancheSize,
        const DateTimeArray  pastTrancheLossDates,
        const DoubleArray    pastTrancheLosses,
        const double         pastTrancheLoss,
        IForwardRatePricerSP model)
{
    static const string method = "VanillaCreditFeeLeg::generateKnownCashFlows";

    //for now
    throw ModelException(method, "Not yet implemented");
}

/** Return known cash flows corresponding to a non-defaulted CDS */
CashFlowArraySP VanillaCreditFeeLeg::generateKnownCashFlows(
        IForwardRatePricerSP model) const
{
    static const string method = "VanillaCreditFeeLeg::generateKnownCashFlows";

    //for now
    throw ModelException(method, "Not yet implemented");
}

/** Estimates the known cash flows (so they are not really "known")
    * corresponding to a CDO tranche - takes into account estimated losses
    * in the future */
CashFlowArraySP VanillaCreditFeeLeg::estimateKnownCashFlows(
        const double         initialTrancheSize,
        const DateTimeArray  pastTrancheLossDates,
        const DoubleArray    pastTrancheLosses,
        IForwardRatePricerSP model)
{
    static const string method = "VanillaCreditFeeLeg::estimateKnownCashFlows";

    //for now
    throw ModelException(method, "Not yet implemented");
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
CashFlowArraySP VanillaCreditFeeLeg::generateKnownCashFlowsGivenDefault(
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
    static const string method(
        "VanillaCreditFeeLeg::generateKnownCashFlowsGivenDefault");

    //for now
    throw ModelException(method, "Not yet implemented");
}

/** Returns the pay date which is last in terms of time */
DateTime VanillaCreditFeeLeg::getLastPayDate() const
{
    //simply the last pay date
    DateTimeArraySP payDates = getCashFlowDates();

    return payDates->back();
}

/** Returns the observation date which is last in terms of time */
DateTime VanillaCreditFeeLeg::getLastObservationDate() const
{
    static const string method = "VanillaCreditFeeLeg::getLastObservationDate";

    //for now
    throw ModelException(method, "Not yet implemented");
}

/** When to stop tweaking for Yield Curve type tweaks */
DateTime VanillaCreditFeeLeg::lastYCSensDate(const DateTime& currentLastDate) const
{
    return getLastPayDate();
}

/** Feedback method for getting information about a fee cashflow */
void VanillaCreditFeeLeg::getActiveFee(
    const DateTime&      withRespectTo,       // (I) get the fee whose accrual period contains this date
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
double VanillaCreditFeeLeg::pv(const DateTime&         valueDate,
                    const DateTime&         defaultDeterminationDate,
                    const DateTime&         accrualPaymentDate,
                    const YieldCurveConstSP discount,
                    const bool              computeAccrual,
                    IForwardRatePricerSP    model) const
{
    static const string method = "VanillaCreditFeeLeg::pv";

    //for now
    throw ModelException(method, "Not yet implemented");
}

/** Compute PV of fee leg at valuation date. */
double VanillaCreditFeeLeg::getFeeLegPV(const DateTime&              valuationDate, 
                                        const DateTime&              earliestRiskyDate,
                                        const DateTime&              latestRiskyDate,
                                        const IDiscountCurve&        discount,
                                        const IDiscountCurveRisky&   crv,
                                        const IDecretionCurveConstSP prepay,
                                        const bool                   includeAccrued,
                                        const DayCountConventionSP   dcc,
                                        IForwardRatePricerSP         model) const
{
    const static string method = "VanillaCreditFeeLeg::getFeeLegPV";

    return this->getFeeLegPV(valuationDate, valuationDate,
                             earliestRiskyDate, latestRiskyDate,
                             discount, crv, prepay,
                             includeAccrued, dcc, model);
}

/**
 * Compute PV of fee leg at valuationDate, but for unconditional payment at
 * paymentDate, conditional on no new defaults before valuationDate.
 */
double VanillaCreditFeeLeg::getFeeLegPV(const DateTime&              valuationDate, 
                                        const DateTime&              paymentDate, 
                                        const DateTime&              earliestRiskyDate,
                                        const DateTime&              latestRiskyDate,
                                        const IDiscountCurve&        discount,
                                        const IDiscountCurveRisky&   crv,
                                        const IDecretionCurveConstSP prepay,
                                        const bool                   includeAccrued,
                                        const DayCountConventionSP   dcc,
                                        IForwardRatePricerSP         model) const
{
    const static string method = "VanillaCreditFeeLeg::getFeeLegPV";

    return this->getFeeLegPV(valuationDate, paymentDate,
                             earliestRiskyDate, latestRiskyDate,
                             discount, crv, prepay,
                             includeAccrued, dcc, false, model);
}

double VanillaCreditFeeLeg::getFeeLegPV(const DateTime&              valuationDate, 
                                        const DateTime&              paymentDate, 
                                        const DateTime&              earliestRiskyDate,
                                        const DateTime&              latestRiskyDate,
                                        const IDiscountCurve&        discount,
                                        const IDiscountCurveRisky&   crv,
                                        const IDecretionCurveConstSP prepay,
                                        const bool                   includeAccrued,
                                        const DayCountConventionSP   dcc,
                                        bool                         defaultValueOnly,
                                        IForwardRatePricerSP         model) const
{
    // NB PREPAY CURVE IS CURRENTLY IGNORED

    double                            sp             = 1.0; // survival probability
    double                            df             = 1.0; // ir discount factor
    double                            survivalPV     = 1.0;
    double                            defaultValuePV = 0.0; // pv if you only want the defaultValue

    const static string               method("VanillaCreditFeeLeg::getFeeLegPV");

    // If default value only and no accrued recovery, then no value at all
    if (defaultValueOnly && !payAccruedFee)
    {
        return 0.0;
    }

    if (valuationDate < paymentDate)
    {
        DateTime toRiskyDate = (paymentDate > latestRiskyDate) ?
                               latestRiskyDate :
                               paymentDate;

        sp = crv.survivalProb(valuationDate, toRiskyDate);
        if (sp==0.0) 
        {
            throw ModelException(method, "Unexpected zero survival probability.");
        }
        //df = discount.pv(valuationDate, paymentDate);
        df = crv.pv(valuationDate, paymentDate);
        df /= sp;
    }
    else if (valuationDate > paymentDate)
    {
        throw ModelException(method, "valuationDate must always be on or before paymentDate.");
    }
    else if (valuationDate >= maturityDate)
    {
        // the last payment date of a fee leg is the protection end date
        // passed end of protection date - value is zero
        return 0.0;
    }

    // pv of defaultValue
    if(defaultValueOnly)
    {
        defaultValuePV = crv.annuityPV(*(feePayments),
                                         paymentDate,
                                         IDiscountCurveRisky::RECOVER_0,
                                         delay);

    }

    survivalPV = crv.annuityPV(*(feePayments),
                                 paymentDate,
                                 payAccruedFee ? IDiscountCurveRisky::RECOVER_1 : IDiscountCurveRisky::RECOVER_0,
                                 delay) - defaultValuePV;


    // only get survivalPV if survives until paymentDate, so * by survival prob.
    // if valuationDate < payment date, return pv as of valuationDate
    return df*sp*survivalPV;
}

/** Returns the accrued interest */
double VanillaCreditFeeLeg::getFeeLegAI(
    const DateTime&              valuationDate, 
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
double VanillaCreditFeeLeg::getFeeLegDefaultedPV(
    const DateTime&              valuationDate,
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
    throw ModelException(method, "not yet implemented");
}

/* Computes the interest accrued up to valuationDate */
double VanillaCreditFeeLeg::getAccruedInterest(const DateTime& valuationDate,
                                               IForwardRatePricerSP model) const
{
    const static string method = "VanillaCreditFeeLeg::getAccruedInterest";
    throw ModelException(method, "not yet implemented");
}

/** Retrieve market data for this fee leg */
void VanillaCreditFeeLeg::getMarket(const IModel* model, const MarketData* market)
{
    const static string method = "VanillaCreditFeeLeg::getMarket";

    market->GetReferenceDate(valueDate);

    /* assumes discount, pmtHol, accrualHol is empty */
    /* Get the interest rate discount curve */
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


    if(!accrualHol.isEmpty())
    {
        accrualHol.getData(model, market);
    }
    else
    {
        // default it to weekends only
        accrualHol.setObject(MarketObjectSP(Holiday::weekendsOnly()));
    }
}

//-------------------------------
// ICreditFeeLegGenerator methods
//-------------------------------

ICreditFeeLegSP VanillaCreditFeeLeg::generateCreditFeeLeg(const DateTime& startDate,
                                                          const DateTime& endDate,
                                                          double feeRate) const
{
    const static string method = "VanillaCreditFeeLeg::generateCreditFeeLeg";
    VanillaCreditFeeLegSP newFeeLeg(dynamic_cast<VanillaCreditFeeLeg*>(this->clone()));

    if (!newFeeLeg)
    {
        throw ModelException(method, "Generating fee leg failed.");
    }

    /* reset effective date, maturity date and fee rate */
    newFeeLeg->effectiveDate = startDate;
    newFeeLeg->maturityDate  = endDate;
    newFeeLeg->feeRate       = feeRate;

    /* regenerate the fee payments */
    newFeeLeg->feePayments   = newFeeLeg->generateFeePayments();

    /* reset fcd */
    newFeeLeg->fcd           = newFeeLeg->getFirstCouponDate();

   /* TODO: reset lcd */

    return newFeeLeg;

}

//-------------------------------
// IFixedRateCreditFeeLeg methods
//-------------------------------

/**Get the fee rate of the fee leg.*/
double VanillaCreditFeeLeg::getRate() const
{
    return feeRate;
}

/**
 * Change the fee rate of the fee leg.
 */
void VanillaCreditFeeLeg::setRate(double newRate)
{
    const static string method = "VanillaCreditFeeLeg::setRate";

    if( newRate == 0.0 && feeRate != 0.0)
    {
        throw ModelException(method, "The new fee must be non-zero unless the old one is zero too!");
    }

    feeRate = newRate;
}

//---------------------
// Theta::Shift methods
//---------------------

/** Implementation of the Theta Shift interface */
bool VanillaCreditFeeLeg::sensShift(Theta* shift)
{
    static const string method = "VanillaCreditFeeLeg::sensShift";

    return true;
}

//--------------------------------------
// VanillaCreditFeeLeg methods (private)
//--------------------------------------

VanillaCreditFeeLeg::VanillaCreditFeeLeg() :CObject(TYPE),
                                            baseNotional(1.0),
                                            feeRate(1.0),
                                            payAccruedFee(true)
{

    /* default stub payment type to SIMPLE */
    StubSP stub(StubFactory::make("Simple"));
    stubPaymentType = stub;

    /* generate the fee cash flow payments */
    //disactivated because the object is empty at this point
    //and the call will fail - MC 21Jun06
    //feePayments = generateFeePayments();

}

VanillaCreditFeeLeg::VanillaCreditFeeLeg(const DateTime& startDate,
                                         const DateTime& endDate,
                                         double feeRate):CObject(TYPE),effectiveDate(startDate),
                                         maturityDate(endDate),
                                         feeRate(feeRate)
{
    /* default stub payment type to SIMPLE */
    StubSP stub(StubFactory::make("Simple"));
    stubPaymentType = stub;

    /* generate the fee cash flow payments */
    feePayments = generateFeePayments();
}

CashFlowArraySP VanillaCreditFeeLeg::generateFeePayments()
{
    static const string method = "VanillaCreditFeeLeg::generateFeePayments";
    //CashFlowArraySP feeCFLs;
    int feePaymentSize = 0;

    // TO ADD
    /* Create a DateListGenerator to take both a front stub and/or back stubs */
    /* Can then use StubPlacement and pass this in as an arg to this
     * DateListGenerator as the field frontStub = "sf", etc.  This would in turn
     * be deciphered by the StubPlacement class where the stub is.
     */

    CashFlowArraySP cfl (new CashFlowArray(SwapTool::cashflows(effectiveDate,
                                           maturityDate,
                                           stubPaymentType.get(),
                                           false, // stub at end
                                           accrualBdc.get(),
                                           pmtBdc.get(),
                                           pmtHol.get(),
                                           true, // used to be false - keep start date,
                                           false, // add redemption/final
                                           feeRate,
                                           paymentFreq->toMonths(),
                                           "M",
                                           pmtDcc.get())));


    // Multiply by notional and change sign
    for (int i=0; i<cfl->size(); i++)
    {
        (*cfl)[i].amount *= baseNotional;
    }

    feePayments = cfl;


    if (!feePayments)
    {
        throw ModelException(method, "Could not generate fee payments.");
    }
    /* there must be at least ONE fee */
    feePaymentSize = feePayments->getLength();
    if(feePaymentSize < 1)
    {
        throw ModelException(method, "Less than one fee payment generated.");
    }

    return feePayments;
}

/* Determine the first full coupon date */
DateTime VanillaCreditFeeLeg::getFirstCouponDate()const
{
    int  numDates;
    DateTimeArray feeDates;

    static const string method = "VanillaCreditFeeLeg::getFirstCouponDate";

    /* if fcd has been set, return it otherwise find the first cpn date */
    if(!fcd.empty())
    {
        return fcd;
    }

    feeDates = CashFlow::dates(*feePayments);

    /* Must have at least one fee cashflow! */
    numDates = feeDates.size();
    if(numDates < 1)
    {
        throw ModelException(method, "Need at least one fee date.");
    }

    /* return the first coupon payment date */
    return feeDates[0];
}

/** Returns the leg notional */
double VanillaCreditFeeLeg::getFeeLegNotional() const
{
    return baseNotional;
}

/** Returns the leg notional */
void VanillaCreditFeeLeg::setFeeLegNotional(double newNotional)
{
    baseNotional = newNotional;
}

void VanillaCreditFeeLeg::load(CClassSP& clazz){
    clazz->setPublic();
    REGISTER(VanillaCreditFeeLeg, clazz);
    SUPERCLASS(CObject);
    IMPLEMENTS(Theta::Shift);
    IMPLEMENTS(ICreditFeeLeg);
    IMPLEMENTS(ICreditFeeLegGenerator);
    IMPLEMENTS(IFixedRateCreditFeeLeg);
    EMPTY_SHELL_METHOD(defaultConstructor);

    FIELD(baseNotional,               "initial notional ");
    FIELD(effectiveDate,              "protection start date ");
    FIELD(valueDate,                  "valuation date ");
    FIELD(maturityDate,               "maturity date ");
    FIELD(feeRate,                    "fee");
    FIELD       (paymentFreq,                "payment frequency");
    FIELD(fcd,                        "first coupon date");
    FIELD(lcd,                        "last coupon date");
    FIELD       (pmtDcc,                     "payment day count convention");
    FIELD       (accrualDcc,                 "accrual day count convention");
    FIELD       (pmtBdc,                     "payment bad day convention");
    FIELD       (accrualBdc,                 "accrual bad day convention");
    FIELD       (valuationBdc,               "valuation bad day convention");
    FIELD(pmtHol,                     "payment holiday convention");
    FIELD(accrualHol,                 "accrual holiday convention");
    FIELD(discount,                   "Discount curve");
    FIELD(payAccruedFee,              "pay accrued fee on default");
    FIELD(delay,                      "payment offset");

    FIELD_MAKE_OPTIONAL(baseNotional);
    FIELD_MAKE_OPTIONAL(valueDate);
    FIELD_MAKE_OPTIONAL(fcd);
    FIELD_MAKE_OPTIONAL(lcd);
    FIELD_MAKE_OPTIONAL(accrualDcc);
    FIELD_MAKE_OPTIONAL(accrualBdc);
    FIELD_MAKE_OPTIONAL(valuationBdc);
    FIELD_MAKE_OPTIONAL(accrualHol);
    FIELD_MAKE_OPTIONAL(discount);
    FIELD_MAKE_OPTIONAL(payAccruedFee);
    FIELD_MAKE_OPTIONAL(delay);
}

IObject* VanillaCreditFeeLeg::defaultConstructor()
{
     return new VanillaCreditFeeLeg();
}


CClassConstSP const VanillaCreditFeeLeg::TYPE =
CClass::registerClassLoadMethod(
    "VanillaCreditFeeLeg", typeid(VanillaCreditFeeLeg), load);

/* to ensure class is linked in */
bool VanillaCreditFeeLegLoad(){
    return VanillaCreditFeeLeg::TYPE != NULL;
}


ICreditFeeLegSVGenSP
VanillaCreditFeeLeg::createSVGen(
	const DateTime& valueDate,
	const DateTimeArray& modifiedFeeLegObservationDates,
	const DateTimeArray& productTimeline, //product timeline
	const IntArray& dateToDiscFactorIndex, //map from a date to an index on the DiscountFactor SV
	double lossConfigNotional
	) const
{
	throw ModelException("VanillaCreditFeeLeg::createSVGen","Method not implemented");
}







/****************** ICreditVanillaInstrument ******************/
/*
CashFlowArraySP VanillaCreditFeeLeg::getCashFlows() const
{
    return feePayments;
}
*/

/** Calculates PV given a valuation date.
 *  Settlement is whatever the instrument determines it to be (e.g. T+1 CDS; T+3 bond, etc.)
 *  relative to the valuation date. Value should be conditional on no new defaults
 *  before valuationDate.
 */
/*
double VanillaCreditFeeLeg::getPV(const DateTime& valuationDate,
                                  const IDiscountCurveRisky& crv) const
{
    return getFeeLegPV(valuationDate, crv);
}
*/

/** Calculates the PV on a given valuationDate for unconditional settlement
 *  and payment on paymentDate.
 */
/*
double VanillaCreditFeeLeg::getPV(const DateTime&                 valuationDate,
                                  const DateTime&                 paymentDate,
                                  const IDiscountCurveRisky&      crv) const
{
    return getFeeLegPV(valuationDate, paymentDate, crv);
}
*/

/* Is not perpetual */
/*
bool VanillaCreditFeeLeg::hasFiniteMaturity() const
{
    return true;
}
*/

/* Returns the last accrual date for the fee payments. */
/*
DateTime VanillaCreditFeeLeg::getMaturity() const
{
    return maturityDate;
}
*/

/**
 *  Returns the earliest possible "interesting" date: this will be the earliest of the
 *  start of the first accrual period, the start of the contingent leg, the first
 *  cash-flow, etc.
*/
/*
DateTime VanillaCreditFeeLeg::getStartDate() const
{
    return effectiveDate;
}
*/

/**
 *  Return the notional of the instrument. If the instrument is amortizing, this
 *  will be the base notional. The idea is that getPV(...)/getNotional() should be
 *  the price of the instrument, in some meaningful sense.
 */
/*
double VanillaCreditFeeLeg::getNotional() const
{
    return baseNotional;
}
*/


/*
void VanillaCreditFeeLeg::GetMarket(const IModel* model, const CMarketDataSP market)
{
    this->getMarket(model, market.get());
}
*/

/*
CSensControl* VanillaCreditFeeLeg::AlterControl(const IModel*          modelParams,
                                                const CSensControl*    sensControl) const
{
    return 0;
}
*/

/** Returns the value date (aka today) the instrument is currently
    pricing for */
/*
DateTime VanillaCreditFeeLeg::getValueDate() const
{
    return valueDate;
}
*/

/*
bool VanillaCreditFeeLeg::priceDeadInstrument(CControl* control,
                                              CResults* results) const
{
    return false;
}
*/
/*
string VanillaCreditFeeLeg::discountYieldCurveName() const
{
	return discount.getName();
}
*/













DRLIB_END_NAMESPACE
